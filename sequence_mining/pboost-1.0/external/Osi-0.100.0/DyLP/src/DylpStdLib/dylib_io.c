/*
  This file is part of the support library for the Dylp LP distribution.

        Copyright (C) 2005 -- 2007 Lou Hafer

        School of Computing Science
        Simon Fraser University
        Burnaby, B.C., V5A 1S6, Canada
        lou@cs.sfu.ca

  This code is licensed under the terms of the Common Public License (CPL).
*/

/*
  This file contains a basic i/o package, designed to isolate the client
  program from ANSI C file i/o primitives and augment their functionality.
  The package uses an integer stream id, and keeps its own records on open
  files. Multiple open/close requests to the same file are tracked. Output to
  a terminal and/or log file (and, in particular, coordinated output to one
  or both) is handled with reasonable sophistication. There is a moderately
  capable lexer, a number of formatted output routines, and various random
  utilities.

  This package requires the set of error reporting/logging utilities in
  errs.c.
*/

#include "dylib_std.h"

static char sccsid[] UNUSED = "@(#)io.c	3.14	11/11/04" ;
static char svnid[] UNUSED = "$Id: dylib_io.c 245 2008-07-08 13:40:02Z lou $" ;

#include <stdio.h>
#include <stdarg.h>
#include <string.h>
#include "dylib_io.h"
#include "dylib_errs.h"


/*
  These declarations are used only by dyio_outfmt_, the version of dyio_outfmt
  that's designed to be called from Fortran. Unless you're the proud
  possessor of some ancient LP software, you likely don't want to know about
  xmp.h.
*/

#ifdef _DYLIB_FORTRAN
#  include "dylib_fortran.h"
#endif

#ifndef MALLOC
#define MALLOC(zz_sze_zz) malloc(zz_sze_zz)
#define CALLOC(zz_cnt_zz,zz_sze_zz) calloc(zz_cnt_zz,zz_sze_zz)
#define FREE(zz_fptr_zz) free(zz_fptr_zz)
#endif

#define MAXPATH FILENAME_MAX	/* Maximum length of full path name */

/*
  The structure filblks keeps track of open files. In addition to the FILE
  structure, a filblks entry stores the channel modes and the full path
  (directory path plus file name) of the file. The i/o id used by the client
  is the index in this array. The size of this array is set in dyio_ioinit
  to be 2 less than the number of file descriptors available to the process.
  The -2 accounts for two files used in the error reporting/logging package
  (the error message text file, and the error message log file). They are
  handled without reference to this package to allow these routines to use
  the error reporting/logging package without creating a potential infinite
  recursion.

  Field		Content
  -----		-------
  stream	ANSI C stdio FILE structure
  modes		Flags specifying the channel mode.
  refcnt	The number of opens that resolved to this file
  dname		Directory path for this file.
  fname		File name.

  Note that filblks is indexed from 1 to maxfiles, so that we can use NULL (0)
  to indicate no stream id.
*/

typedef struct { FILE *stream ;
                 flags modes ;
		 int refcnt ;
                 char *dname ;
                 char *fname ; } filblk_struct ;

static filblk_struct *filblks ;
static ioid maxfiles ;

#define INVALID_IOID(zz_id) ((zz_id < (ioid) 0) || (maxfiles < zz_id))
#define INVALID_STREAMID(zz_id) ((zz_id < (ioid) 1) || (maxfiles < zz_id))

/*
  Flags for the modes field.

  Flag		Description
  ----		-----------
  io_active	This entry is active.
  io_line	Line oriented i/o, with '\n' treated as a delimiter character.
  io_free	Free i/o, with '\n' treated as white space.
  io_read	Stream is opened for input.
  io_write	Stream is opened for output.
*/

#define io_active	1<<0
#define io_line		1<<1
#define io_free		1<<2
#define io_read		1<<3
#define io_write	1<<4



/*
  The chartab provides the character translation and class codes on input.
  nxtchar reads characters from a stream and uses chartab to provide the
  correct character translation and class codes to the various scan routines.
  The description of the chartab_struct structure used for chartab entries is
  as follows:

  Field		Description
  -----		-----------
  chrtrn	Character supplied when this ascii code is read.
  class1	Primary character class.
  class2	Secondary character class, used during construction of
		identifiers.
  
  The values for class1 and class2 are drawn from the following set, defined
  below as the enum chrclass:

  Value		Description
  -----		-----------
  CCINV		(INVtermin) Characters in this class terminate any symbol
		being formed, but are otherwise ignored.
  CCDIG		(DIGit) The obvious set of digits.
  CCIDB		(ID Begin) Characters which can start an identifier.
  CCDEL		(DELimiter) Single character delimiters.
  CCIG		(IGnore) Characters which are alway ignored.
  CCIDC		(ID Continuation) Characters which are valid in an identifier
		but cannot start one.
  CCEOF		Indicates end-of-file was seen.
  CCERR		Indicates i/o error was seen.
*/

typedef enum {CCINV,CCDIG,CCIDB,CCDEL,CCIG,CCIDC,CCEOF,CCERR} chrclass ;

typedef struct { chrclass class1 ;
		 chrclass class2 ;
		 char chrtrn ; } chartab_struct ;

/*
  A few special case chartab_struct's used by nxtchar to handle special cases.
  chartab_eof and chartab_err are returned on EOF and i/o error respectively.
  chartab_nl is returned when '\n' is seen on a stream operating in line mode.
*/

static
chartab_struct chartab_eof = { CCEOF,CCEOF,'\377' },
		      chartab_err = { CCERR,CCERR,'\377' },
		      chartab_nl  = { CCDEL,CCDEL,'\n' } ;

static
chartab_struct chartab[128] = { { CCIG,CCIG,'\0' },	 	/* NULL, ^@ */
				{ CCIG,CCIG,'\0' },		/* SOH,  ^A */
				{ CCIG,CCIG,'\0' },		/* STX,  ^B */
				{ CCIG,CCIG,'\0' },		/* ETX,  ^C */
				{ CCIG,CCIG,'\0' },		/* EOT,  ^D */
				{ CCIG,CCIG,'\0' },		/* ENQ,  ^E */
				{ CCIG,CCIG,'\0' },		/* ACK,  ^F */
				{ CCIG,CCIG,'\0' },		/* BEL,  ^G */
				{ CCIG,CCIG,'\0' },		/* BS,	 ^H */
				{ CCINV,CCINV,'\t' },		/* HT,	 ^I */
				{ CCINV,CCINV,'\n' },		/* LF,   ^J */
				{ CCINV,CCINV,'\013' },		/* VT,   ^K */
				{ CCINV,CCINV,'\f' },		/* FF,   ^L */
				{ CCINV,CCINV,'\r' },		/* CR,   ^M */
				{ CCIG,CCIG,'\0' },		/* SO,   ^N */
				{ CCIG,CCIG,'\0' },		/* SI,   ^O */
				{ CCIG,CCIG,'\0' },		/* DLE,  ^P */
				{ CCIG,CCIG,'\0' },		/* DC1,  ^Q */
				{ CCIG,CCIG,'\0' },		/* DC2,  ^R */
				{ CCIG,CCIG,'\0' },		/* DC3,  ^S */
				{ CCIG,CCIG,'\0' },		/* DC4,  ^T */
				{ CCIG,CCIG,'\0' },		/* NAK,  ^U */
				{ CCIG,CCIG,'\0' },		/* SYN,  ^V */
				{ CCIG,CCIG,'\0' },		/* ETB,  ^W */
				{ CCIG,CCIG,'\0' },		/* CAN,  ^X */
				{ CCIG,CCIG,'\0' },		/* EOM,  ^Y */
				{ CCIG,CCIG,'\0' },		/* SUB,  ^Z */
				{ CCIG,CCIG,'\0' },		/* ESC,  ^[ */
				{ CCIG,CCIG,'\0' },		/* FS,	 ^\ */
				{ CCIG,CCIG,'\0' },		/* GS,   ^] */
				{ CCIG,CCIG,'\0' },		/* RS,   ^^ */
				{ CCIG,CCIG,'\0' },		/* US,   ^_ */
				{ CCINV,CCINV,' ' },		/* SP */
				{ CCDEL,CCDEL,'!' },
				{ CCDEL,CCDEL,'"' },
				{ CCDIG,CCDIG,'#' },
				{ CCDEL,CCDEL,'$' },
				{ CCDEL,CCIDC,'%' },
				{ CCDEL,CCDEL,'&' },
				{ CCDIG,CCDIG,'\047' },		/* ' */
				{ CCDEL,CCDEL,'(' },
				{ CCDEL,CCDEL,')' },
				{ CCDEL,CCDEL,'*' },
				{ CCDIG,CCDIG,'+' },
				{ CCDEL,CCDEL,',' },
				{ CCDIG,CCDIG,'-' },
				{ CCDIG,CCIDC,'.' },
				{ CCDEL,CCDEL,'/' },
				{ CCDIG,CCIDC,'0' },
				{ CCDIG,CCIDC,'1' },
				{ CCDIG,CCIDC,'2' },
				{ CCDIG,CCIDC,'3' },
				{ CCDIG,CCIDC,'4' },
				{ CCDIG,CCIDC,'5' },
				{ CCDIG,CCIDC,'6' },
				{ CCDIG,CCIDC,'7' },
				{ CCDIG,CCIDC,'8' },
				{ CCDIG,CCIDC,'9' },
				{ CCDEL,CCDEL,':' },
				{ CCDEL,CCDEL,';' },
				{ CCDEL,CCDEL,'<' },
				{ CCDEL,CCDEL,'=' },
				{ CCDEL,CCDEL,'>' },
				{ CCDEL,CCDEL,'?' },
				{ CCDEL,CCDEL,'@' },
				{ CCIDB,CCIDC,'A' },
				{ CCIDB,CCIDC,'B' },
				{ CCIDB,CCIDC,'C' },
				{ CCIDB,CCIDC,'D' },
				{ CCIDB,CCIDC,'E' },
				{ CCIDB,CCIDC,'F' },
				{ CCIDB,CCIDC,'G' },
				{ CCIDB,CCIDC,'H' },
				{ CCIDB,CCIDC,'I' },
				{ CCIDB,CCIDC,'J' },
				{ CCIDB,CCIDC,'K' },
				{ CCIDB,CCIDC,'L' },
				{ CCIDB,CCIDC,'M' },
				{ CCIDB,CCIDC,'N' },
				{ CCIDB,CCIDC,'O' },
				{ CCIDB,CCIDC,'P' },
				{ CCIDB,CCIDC,'Q' },
				{ CCIDB,CCIDC,'R' },
				{ CCIDB,CCIDC,'S' },
				{ CCIDB,CCIDC,'T' },
				{ CCIDB,CCIDC,'U' },
				{ CCIDB,CCIDC,'V' },
				{ CCIDB,CCIDC,'W' },
				{ CCIDB,CCIDC,'X' },
				{ CCIDB,CCIDC,'Y' },
				{ CCIDB,CCIDC,'Z' },
				{ CCDEL,CCDEL,'[' },
				{ CCDEL,CCDEL,'\\' },
				{ CCDEL,CCDEL,']' },
				{ CCDEL,CCDEL,'^' },
				{ CCDEL,CCIDC,'_' },
				{ CCDEL,CCDEL,'`' },
				{ CCIDB,CCIDC,'a' },
				{ CCIDB,CCIDC,'b' },
				{ CCIDB,CCIDC,'c' },
				{ CCIDB,CCIDC,'d' },
				{ CCIDB,CCIDC,'e' },
				{ CCIDB,CCIDC,'f' },
				{ CCIDB,CCIDC,'g' },
				{ CCIDB,CCIDC,'h' },
				{ CCIDB,CCIDC,'i' },
				{ CCIDB,CCIDC,'j' },
				{ CCIDB,CCIDC,'k' },
				{ CCIDB,CCIDC,'l' },
				{ CCIDB,CCIDC,'m' },
				{ CCIDB,CCIDC,'n' },
				{ CCIDB,CCIDC,'o' },
				{ CCIDB,CCIDC,'p' },
				{ CCIDB,CCIDC,'q' },
				{ CCIDB,CCIDC,'r' },
				{ CCIDB,CCIDC,'s' },
				{ CCIDB,CCIDC,'t' },
				{ CCIDB,CCIDC,'u' },
				{ CCIDB,CCIDC,'v' },
				{ CCIDB,CCIDC,'w' },
				{ CCIDB,CCIDC,'x' },
				{ CCIDB,CCIDC,'y' },
				{ CCIDB,CCIDC,'z' },
				{ CCDEL,CCDEL,'{' },
				{ CCDEL,CCDEL,'|' },
				{ CCDEL,CCDEL,'}' },
				{ CCIG,CCIG,'\176' },		/* ~ */
				{ CCIG,CCIG,'\177' } } ;	/* DEL */

/*
  We need a few special lexemes too. lex_eof and lex_err are used to indicate
  EOF and lexical scan error, respectively. lex_nil is the null lexeme.
*/

static lex_struct lex_eof = {DY_LCEOF,"end-of-file"},
		  lex_err = {DY_LCERR,"lexical scan error"},
		  lex_nil = {DY_LCNIL,NULL} ;
 


bool dyio_ioinit (void)

/*
  This routine initialises the i/o data structures. It should be called only
  once, at the start of the client's program. It does the following:

  * Initialise the filblks array. We first get the size of the file descriptor
    table and acquire space for filblks. The first three entries are set up
    for stdin, stdout, and stderr, the streams opened by default under Unix.

  * A possible fourth entry is made if the user has opened an error log
    file.  We need this entry here so that a later call to dyio_openfile can
    find it (the particular case of interest would be a desire to log normal
    and error messages to the same file). To maintain information hiding,
    there is a routine, errlogq, which returns the required information.

  Parameters:
    errlogpath:	the path name of the error log file; used only if errlogq
		returns a non-NULL FILE handle

  Returns: TRUE if the initialisation succeeds, FALSE otherwise
*/

{ int len ;
  filblk_struct *filblk ;
  char *fname,*errlogpath,*tmp ;
  const char *rtnnme = "dyio_ioinit" ;

  extern FILE *errlogq(char **elogpath) ;

/*
  See the comments on the filblks array at the beginning of the file, for the
  explanation of the magic number 2. It's possible that the o/s may have
  imposed a more stringent limit on the client process.
*/
  maxfiles = FOPEN_MAX-2 ;
  filblks = (filblk_struct *) CALLOC(maxfiles+1,sizeof(filblk_struct)) ;
  filblk = &filblks[1] ;
  filblk->stream = stdin ;
  filblk->dname = NULL ;
  filblk->fname = "stdin" ;
  setflg(filblk->modes,io_active|io_free|io_read) ;
  filblk->refcnt = 1 ;

  filblk = &filblks[2] ;
  filblk->stream = stdout ;
  filblk->dname = NULL ;
  filblk->fname = "stdout" ;
  setflg(filblk->modes,io_active|io_free|io_write) ;
  filblk->refcnt = 1 ;

  filblk = &filblks[3] ;
  filblk->stream = stderr ;
  filblk->dname = NULL ;
  filblk->fname = "stderr" ;
  setflg(filblk->modes,io_active|io_free|io_write) ;
  filblk->refcnt = 1 ;
/*
  Check for an error log file, and make a filblks entry for it, if necessary.
*/
  filblk = &filblks[4] ;
  filblk->stream = errlogq(&errlogpath) ;
  if (filblk->stream != NULL)
  { if (errlogpath == NULL)
    { errmsg(14,rtnnme) ;
      errmsg(1,rtnnme,__LINE__) ;
      return(FALSE) ; }
    fname = strrchr(errlogpath,'/') ;
    if (fname == NULL)
    { filblk->dname = NULL ;
      fname = errlogpath ; }
    else
    { len = fname-errlogpath ;
      tmp = (char *) MALLOC(len+1) ;
      (void) strncpy(tmp,errlogpath,len) ;
      tmp[len] = '\0' ;
      filblk->dname = tmp ;
      fname++ ; }
  tmp = (char *) MALLOC(strlen(fname)+1) ;
  (void) strcpy(tmp,fname) ;
  filblk->fname = tmp ;
  setflg(filblk->modes,io_active|io_write) ;
  filblk->refcnt = 1 ; }

  return(TRUE) ; }


void dyio_ioterm (void)

/*
  This routine cleans up file management data structures allocated by the
  i/o package.

  Parameters: none

  Returns: undefined.
*/

{ int id ;

/*
  Data structures allocated by dyio_ioinit to track open files.
*/
  for (id = 4 ; id <= maxfiles ; id++)
  { if (filblks[id].dname != NULL) FREE(filblks[id].dname) ;
    if (filblks[id].fname != NULL) FREE(filblks[id].fname) ; }
  FREE(filblks) ;

  return ; }




static bool rwmodecmp (filblk_struct *filblk, const char *mode)

/*
  This routine compares the r/w mode flags in the filblks entry for a stream
  with the r/w mode information in mode.

  Parameters:
    filblk:	filblks entry to be checked.
    mode:	Mode parameters as supplied for the standard i/o routine fopen.

  Returns: TRUE if the r/w mode information matches, FALSE otherwise.
*/

{ const char *rtnnme = "rwmodecmp" ;
  flags modes ;
  bool rw ;

  if (mode == NULL)
  { errmsg(2,rtnnme,"r/w mode") ;
    return (FALSE) ; }
  if (filblk == NULL)
  { errmsg(2,rtnnme,"filblk") ;
    return (FALSE) ; }
  modes = filblk->modes ;

  rw = (mode[1] == '+')?TRUE:FALSE ;
  switch (mode[0])
  { case 'r':
    { if (flgoff(modes,io_read) == TRUE) return (FALSE) ;
      if (rw == TRUE && flgoff(modes,io_write) == TRUE) return (FALSE) ;
      return (TRUE) ; }
    case 'a':
    case 'w':
    { if (flgoff(modes,io_write) == TRUE) return (FALSE) ;
      if (rw == TRUE && flgoff(modes,io_read) == TRUE) return (FALSE) ;
      return (TRUE) ; }
    default:
    { errmsg(4,rtnnme,"r/w mode",mode) ;
      return (FALSE) ; } } }



static bool setrwmode (filblk_struct *filblk, char *mode)

/*
  This routine sets the r/w mode flags in the filblks entry for a stream.

  Parameters:
    filblk:	Filblks entry to be set.
    mode:	Mode parameters as supplied to the standard i/o routine fopen.

  Returns: TRUE if the mode is set successfully, FALSE otherwise.
*/

{ const char *rtnnme = "setrwmode" ;
  bool rw ;

  if (mode == NULL)
  { errmsg(2,rtnnme,"r/w mode") ;
    return (FALSE) ; }
  if (filblk == NULL)
  { errmsg(2,rtnnme,"filblk") ;
    return (FALSE) ; }

  rw = (mode[1] == '+')?TRUE:FALSE ;
  switch (mode[0])
  { case 'r':
    { if (rw == TRUE)
	setflg(filblk->modes,io_read|io_write) ;
      else
	setflg(filblk->modes,io_read) ;
      return (TRUE) ; }
    case 'a':
    case 'w':
    { if (rw == TRUE)
	setflg(filblk->modes,io_read|io_write) ;
      else
	setflg(filblk->modes,io_write) ;
      return (TRUE) ; }
    default:
    { errmsg(4,rtnnme,"r/w mode",mode) ;
      return (FALSE) ; } } }



bool dyio_setmode (ioid id, char mode)

/*
  This routine is provided as an information-hiding function so that the
  client can switch a stream between line and free mode.

  Parameters:
    id:		i/o id.
    mode:	Desired mode, 'f' for free, 'l' for line.

  Returns: TRUE if the mode was successfully set, FALSE otherwise.
*/

{ filblk_struct *filblk ;
  const char *rtnnme = "dyio_setmode" ;
/*
  Check to make sure the stream ID is OK.
*/
  if (INVALID_STREAMID(id))
  { errmsg(5,rtnnme,"stream id",id) ;
    return (FALSE) ; }
  filblk = &filblks[id] ;
  if (flgon(filblk->modes,io_active) == FALSE)
  { errmsg(15,rtnnme,id) ;
    return (FALSE) ; }
  if (flgon(filblk->modes,io_read) == FALSE)
  { errmsg(16,rtnnme,dyio_idtopath(id)) ;
    return (FALSE) ; }
/*
  Now set the mode.
*/
  switch (mode)
  { case 'l':
    case 'L':
    { clrflg(filblk->modes,io_free) ;
      setflg(filblk->modes,io_line) ;
      break ; }
    case 'f':
    case 'F':
    { clrflg(filblk->modes,io_line) ;
      setflg(filblk->modes,io_free) ;
      break ; }
    default:
    { errmsg(3,rtnnme,"scanning mode",mode) ;
      return (FALSE) ; } }
  
  return (TRUE) ; }



const char *dyio_idtopath (ioid id)

/*
  This routine returns the file name associated with the specified stream.

  Parameter:
    id:	i/o id.

  Returns: Path name of the file opened on this stream, or an error string
	   if the stream is not currently active.
*/

{ static char fullpath[MAXPATH] ;
  static const char *badid = "!invalid id!" ;
  filblk_struct *filblk ;
  const char *rtnnme = "dyio_idtopath" ;
/*
  Check to make sure the stream ID is OK.
*/
  if (INVALID_STREAMID(id))
  { errmsg(5,rtnnme,"stream id",id) ;
    return (badid) ; }
  filblk = &filblks[id] ;
  if (flgon(filblk->modes,io_active) == FALSE)
  { errmsg(15,rtnnme,id) ;
    return (badid) ; }
/*
  Now construct the full pathname from the directory and file names.
*/
  fullpath[0] = '\0' ;
  if (filblk->dname != NULL)
  { (void) strcat(fullpath,filblk->dname) ;
    (void) strcat(fullpath,"/") ; }
  (void) strcat(fullpath,filblk->fname) ;

  return (fullpath) ; }



ioid dyio_pathtoid (const char *path, const char *mode)

/*
  This routine searches the filblks array and returns the stream ID associated
  with the full path name given by path. If mode is non-null, the i/o mode of
  the file must also match.

  Parameters:
    path:	Full path name of the file.
    mode:	i/o mode, as passed to the standard i/o routine fopen

  Returns: I/O ID associated with the path name (and mode), or IOID_INV if
	   nothing matched.
*/

{ const char *fname ;
  int ndx,dlen ;
  filblk_struct *filblk ;
  const char *rtnnme = "dyio_pathtoid" ;

  if (path == NULL)
  { errmsg(2,rtnnme,"path") ;
    return (-1) ; }

  fname = strrchr(path,'/') ;
  if (fname == NULL)
  { dlen = 0 ;
    fname = path ; }
  else
  { dlen = fname-path ;
    fname++ ; }

  for (ndx = 1 ; ndx <= maxfiles ; ndx++)
  { filblk = &filblks[ndx] ;
    if (flgoff(filblk->modes,io_active) == TRUE) continue ;
    if (strcmp(filblk->fname,fname) != 0) continue ;
    if (filblk->dname == NULL)
    { if (dlen != 0) continue ; }
    else
    { if (strncmp(filblk->dname,path,dlen) != 0) continue ; }
    if (mode != NULL && rwmodecmp(filblk,mode) == FALSE) continue ;
    return (ndx) ; }
  
  return (IOID_INV) ; }




bool dyio_ttyq (ioid id)

/*
  This routine answers the question of whether or not the stream specified by
  id is attached to a tty. It is simply an interface to isatty.

  isatty and fileno are not ANSI standard functions, but seem pretty common
  among Unix implementations. Your mileage may vary.

  Parameters:
    id:	i/o id

  Returns: TRUE if the stream is attached to a tty, FALSE otherwise; in the
	   event of a error, FALSE is returned.
*/

{ filblk_struct *filblk ;
  const char *rtnnme = "dyio_ttyq" ;

  extern int isatty(int fildes) ;
  /*extern int fileno(FILE *stream) ;*/

/*
  Check to make sure the stream ID is OK.
*/
  if (INVALID_STREAMID(id))
  { errmsg(5,rtnnme,"stream id",id) ;
    return (FALSE) ; }
  filblk = &filblks[id] ;
  if (flgon(filblk->modes,io_active) == FALSE)
  { errmsg(15,rtnnme,id) ;
    return (FALSE) ; }
/*
  Now inquire if the stream is attached to a terminal.
*/
  if (isatty(fileno(filblk->stream)) == 1)
    return (TRUE) ;
  else
    return (FALSE) ; }



ioid dyio_openfile (const char *path, const char *mode)

/*
  This routine sets up the file specified in path for i/o. It acts as the
  interface between the client and the C stdio library. The major differences
  from the standard fopen function are
  
    1): dyio_openfile looks to see if the specified file is already opened with
	compatible r/w mode (see below for 'compatible') and if so it returns
	the ID for the preexisting stream, with no further actions.
    2): Six additional modes "R","R+","W","W+","A","A+" are provided to
	override feature 1, i.e., to force a rewind or truncate operation, as
	appropriate.
    3): The letter `q' (query) is interpreted as `open the file if it exists,
	but it's not an error if it doesn't exist.'
	
  'Compatible' mode means everything has to match except the case of the
  letter.

  Parameters:
    path:	Address of the full path name (directory path plus file name) 
                of the file to be opened.
    mode:	Address of a string giving the file i/o mode, as specified 
                for the standard fopen function, with the extra modes as
		described above.

  Returns: stream ID, or IOID_INV if the file cannot be opened.
*/

{ FILE *handle ;
  filblk_struct *filblk ;
  const char *fname ;
  char *tmp ;
  char mode_var[3] ;
  int len ;
  ioid id ;
  bool mustexist,totop ;
  const char *rtnnme = "dyio_openfile" ;
/*
  Make sure the parameters are ok.
*/
  if (path == NULL)
  { errmsg(2,rtnnme,"path") ;
    return (IOID_INV) ; }
  if (mode == NULL)
  { errmsg(2,rtnnme,"r/w mode") ;
    return (IOID_INV) ; }
/*
  Get a copy of mode so we can play with the string.
*/
  mode_var[0] = mode[0] ;
  mode_var[1] = mode[1] ;
  mode_var[2] = '\0' ;
/*
  If mode_var starts with an upper case letter, change to lower case. For
  modes beginning with 'W', 'R', set the rewind flag. Convert query modes
  to read, but note that the file need not exist.
*/
  mustexist = FALSE ;
  totop = FALSE ;
  switch (mode_var[0])
  { case 'R':
    { mode_var[0] = 'r' ;
      totop = TRUE ;
      mustexist = TRUE ;
      break ; }
    case 'Q':
    { mode_var[0] = 'r' ;
      totop = TRUE ;
      break ; }
    case 'W':
    { mode_var[0] = 'w' ;
      totop = TRUE ;
      break ; }
    case 'A':
    { mode_var[0] = 'a' ;
      break ; }
    case 'r':
    { mustexist = TRUE ;
      break ; }
    case 'q':
    { mode_var[0] = 'r' ;
      break ; }
    case 'w':
    case 'a':
    { break ; }
    default:
    { errmsg(4,rtnnme,"r/w mode",mode) ;
      return (IOID_INV) ; } }

/*
  First check if the file is already open with compatible (lower case
  equivalent) mode. Note that dyio_pathtoid() checks modes. If the file is
  already open, return with the stream id. If additionally the rewind
  flag totop is true, take appropriate action. The rewind flag and/or
  uppercase mode mean nothing if the file is not already open.
*/
  id = dyio_pathtoid(path,mode_var) ; 
  if (id != IOID_INV)
  { if (totop == TRUE) 
    { switch (mode_var[0])
      { case 'r':
	{ rewind(filblks[id].stream) ;
	  break ; }
	case 'w':
	{ fclose(filblks[id].stream) ;
	  filblks[id].stream = fopen(path,mode_var) ; 
	  if (filblks[id].stream == NULL) 
	  { errmsg(10,rtnnme,dyio_idtopath(id),mode_var) ;
	    perror(rtnnme) ;
	    return (IOID_INV) ; }
	  break ; }
	case 'a':
	{ errmsg(1,rtnnme,__LINE__) ;
	  return (IOID_INV) ; } } }
    filblks[id].refcnt++ ;
    return (id) ; }
/*
  The file is not currently open, so find an empty entry in filblks.
*/
  for (id = 1 ;
       id <= maxfiles && flgon(filblks[id].modes,io_active) == TRUE ;
       id++) ;
  if (id > maxfiles)
  { errmsg(13,rtnnme) ;
    return (IOID_INV) ; }
  filblk = &filblks[id] ;
/*
  Now attempt to open the file.
*/
  handle = fopen(path,mode_var) ;
  if (handle == NULL)
  { if (mode_var[0] == 'r' && mustexist == FALSE)
    { return (IOID_NOSTRM) ; }
    else
    { errmsg(10,rtnnme,path,mode_var) ;
      perror(rtnnme) ;
      return (IOID_INV) ; } }
/*
  File is successfully opened, so fill in the information in filblks.
*/
  filblk->stream = handle ;
  setrwmode(filblk,mode_var) ;
  fname = strrchr(path,'/') ;
  if (fname == NULL)
  { filblk->dname = NULL ;
    fname = path ; }
  else
  { fname++ ;
    len = fname-path ;
    tmp = (char *) MALLOC(len+1) ;
    (void) strncpy(tmp,path,len) ;
    tmp[len] = '\0' ;
    filblk->dname = tmp ; }
  tmp = (char *) MALLOC(strlen(fname)+1) ;
  (void) strcpy(tmp,fname) ;
  filblk->fname = tmp ;
  setflg(filblk->modes,io_active|io_free) ;
  filblk->refcnt = 1 ;

  return (id) ; }



bool dyio_isactive (ioid id)

/*
  This routine checks to see if id corresponds to an active stream.

  Parameter:
    id: stream ID obtained from dyio_openfile
  
  Returns: TRUE if the stream is active, FALSE otherwise
*/

{ const char *rtnnme = "dyio_isactive" ;
/*
  We'll still complain if the id isn't valid.
*/
  if (INVALID_IOID(id))
  { errmsg(5,rtnnme,"stream id",id) ;
    return (FALSE) ; }

  if (id == IOID_NOSTRM) return (FALSE) ;

  return (flgon(filblks[id].modes,io_active)) ; }



bool dyio_closefile (ioid id)

/*
  This routine closes a stream and cleans up the filblks entry, iff the
  reference count for the stream has dropped to 0.

  Parameter:
    id:	stream ID obtained from dyio_openfile

  Returns: TRUE if the stream is closed without error, FALSE otherwise.
*/

{ filblk_struct *filblk ;
  bool retval ;
  const char *rtnnme = "dyio_closefile" ;
/*
  Make sure the id is valid.
*/
  if (INVALID_STREAMID(id))
  { errmsg(5,rtnnme,"stream id",id) ;
    return (FALSE) ; }
  filblk = &filblks[id] ;
  if (flgon(filblk->modes,io_active) == FALSE)
  { errmsg(15,rtnnme,id) ;
    return (FALSE) ; }
/*
  Decrement the reference count. If it's non-zero, someone else is still using
  the stream, so just return.
*/
  if (--filblk->refcnt > 0) return (TRUE) ;
/*
  Reference count is 0, so close the stream.
*/
  if (fclose(filblk->stream) == EOF)
  { errmsg(11,rtnnme,dyio_idtopath(id)) ;
    perror(rtnnme) ;
    retval = FALSE ; }
  else
    retval = TRUE ;
/*
  Mark the filblk entry as inactive and free the directory path and file
  name space. We can close the three standard streams (stdin, stdout, stderr)
  but the name strings are static (vid. dyio_ioinit) and should not be freed.
*/
  clrflg(filblk->modes,io_active) ;
  if (id > 3)
  { if (filblk->dname != NULL)
    { FREE(filblk->dname) ;
      filblk->dname = NULL ; }
    FREE(filblk->fname) ;
    filblk->fname = NULL ; }

  return (retval) ; }



bool dyio_chgerrlog (const char *newerrpath, bool echo)

/*
  This routine allows for manipulation of the error logging channel via the
  standard interface provided by io.c. The supplied path is checked against
  open streams, and if one is found it is passed on to the error logging
  package. If there is no open stream for the new log file, it's opened.

  If we had no record of the old error log file, we'll let the error logging
  package do the close. But if we have a record, we'll handle the close here.

  It's possible we've been called solely to change the echo, in which case
  we'll do it up front and bail out.

  Parameters:
    newerrpath:	path name of the file
    echo:	TRUE if errors should be echoed to stderr, FALSE otherwise

  Returns: TRUE if things go well, FALSE otherwise.
*/

{ ioid olderrid,newerrid ;
  bool close ;
  char *olderrpath ;
  FILE *newerrchn ;

  const char *rtnnme = "dyio_chgerrlog" ;

  extern FILE *errlogq(char **elogpath) ;
  extern bool reseterrlogchn(const char *elogpath, FILE *elogchn,
			     bool echo, bool close) ;

/*
  Just here to change the echo?
*/
  if (newerrpath == NULL)
  { (void) reseterrlogchn(NULL,NULL,echo,FALSE) ;
    return (TRUE) ; }
/*
  The first thing we need to do is check for the existence of the present
  error logging file. If the i/o package has no record, we want the error
  logging package to close out the old file. If we have a record, we'll do
  it here.
*/
  (void) errlogq(&olderrpath) ;
  if (olderrpath == NULL)
  { close = FALSE ;
    olderrid = IOID_INV ; }
  else
  { olderrid = dyio_pathtoid(olderrpath,NULL) ;
    if (olderrid == IOID_INV)
    { close = TRUE ; }
    else
    { close = FALSE ; } }
/*
  Now see if the new file is already open. If not, open it.
*/
  newerrid = dyio_pathtoid(newerrpath,"w") ;
  if (newerrid == IOID_INV)
  { newerrid = dyio_openfile(newerrpath,"w") ;
    if (newerrid == IOID_INV)
    { errmsg(10,rtnnme,newerrpath,"w") ;
      return (FALSE) ; } }
  newerrchn = filblks[newerrid].stream ;
/*
  Do the swap. Once we've finished, we may need to deactivate the filblks
  entry for olderrid.
*/
  if (reseterrlogchn(newerrpath,newerrchn,echo,close) == FALSE)
  { errmsg(18,rtnnme,olderrpath,newerrpath) ;
    return (FALSE) ; }
/*
  Did we have a record for this file? If so, call dyio_closefile to tidy up.
*/
  if (olderrid != IOID_INV) (void) dyio_closefile(olderrid) ;
  
  return (TRUE) ; }



long dyio_mark (ioid id)

/*
  This routine marks the current place in the input stream. It uses ftell to
  obtain the current offset. As best I can figure, ftell should always return a
  positive integer, at least on Unix, so the error return here is -1.

  Parameter:
    id: Stream ID from dyio_openfile.

  Returns: A "mark", or -1 in the event of any bogosity.
*/

{ filblk_struct *filblk ;
  long here ;
  const char *rtnnme = "dyio_mark" ;
/*
  Check to make sure the stream ID is OK.
*/
  if (INVALID_STREAMID(id))
  { errmsg(5,rtnnme,"stream id",id) ;
    return (-1) ; }
  filblk = &filblks[id] ;
  if (flgon(filblk->modes,io_active) == FALSE)
  { errmsg(15,rtnnme,id) ;
    return (-1) ; }
/*
  Now try the call to ftell to get the mark.
*/
  here = ftell(filblk->stream) ;
  if (here < 0)
  { errmsg(23,rtnnme,dyio_idtopath(id)) ;
    perror(rtnnme) ; }

  return (here) ; }



bool dyio_backup (ioid id, long there)

/*
  This routine moves the position of the next i/o operation on this stream to
  there.

  Parameters:
    id:		Stream ID from dyio_openfile.
    there:	Mark supplied by dyio_mark.

  Returns: TRUE if the position is successfully set, FALSE otherwise.
*/

{ filblk_struct *filblk ;
  const char *rtnnme = "dyio_backup" ;
/*
  Check to make sure the stream ID is OK.
*/
  if (INVALID_STREAMID(id))
  { errmsg(5,rtnnme,"stream id",id) ;
    return (FALSE) ; }
  filblk = &filblks[id] ;
  if (flgon(filblk->modes,io_active) == FALSE)
  { errmsg(15,rtnnme,id) ;
    return (FALSE) ; }
/*
  Call fseek to position the file pointer to there.
*/
  if (fseek(filblk->stream,there,SEEK_SET) < 0)
  { errmsg(24,rtnnme,dyio_idtopath(id),there) ;
    perror(rtnnme) ;
    return (FALSE) ; }
  
  return (TRUE) ; }



static chartab_struct *nxtchar (FILE *stream, flags modes)

/*
  This routine reads characters from stream and removes characters of class
  CCIG (IGnore).

  Parameters:
    stream:	stream handle
    modes:	modes field from filblk entry for this channel
  
  Returns: chartab_struct for the character read. Note that EOF, i/o error,
	   and '\n' receive special handling.
*/

{ const char *rtnnme = "nxtchar" ;
  int nxt,id ;
/*
  Read characters until not class CCIG or until EOF.
*/
  for (nxt = getc(stream) ;
       nxt != EOF && chartab[nxt].class1 == CCIG ;
       nxt = getc(stream)) ;
/*
  Check for read errors.
*/
  if (nxt == EOF)
  { if (ferror(stream) != 0)
    { for (id = 1 ; id <= maxfiles && filblks[id].stream != stream ; id++) ;
      if (id <= maxfiles)
      { errmsg(12,rtnnme,dyio_idtopath(id)) ; }
      else
      { errmsg(12,rtnnme,"unknown") ;
	errmsg(1,rtnnme,__LINE__) ; }
      perror(rtnnme) ;
      return (&chartab_err) ; }
    else
    { return (&chartab_eof) ; } }
  
  if (nxt == '\n' && flgon(modes,io_line) == TRUE)
    return (&chartab_nl) ;
  else
    return (&chartab[nxt]) ; }




bool dyio_scan (ioid id, const char pattern[], bool rwnd, bool wrap)

/*
  This routine scans the input stream specified by id for the string specified
  in pattern. It leaves the i/o pointer positioned so that the next read from
  the channel will fetch the character following the string matched. rwnd
  specifies whether the stream should be reset to the beginning of the file
  before the scan starts. wrap is valid only if the file is not rewound, and
  specifies whether the search can be continued at EOF by resetting to the 
  beginning of file and resuming the scan. wrap is a bit expensive unless you
  are pretty sure that the desired string will occur between where you're at
  and the end of the file, since there is no cheap way of stopping at the point
  where the scan originally started without digging into the guts of the stdio
  package. Hence the scan continues til the second EOF.

  The algorithm used is the linear substring matching algorithm described in
  Section 9.3 of Aho, Hopcroft, & Ullman, "The Design and Analysis of Computer
  Algorithms". The routine limits the length of pattern arbitrarily to 30 chars
  and will truncate the pattern and complain if the pattern is longer.

  Parameters:
    id:		i/o id.
    pattern:	Character string to be matched.
    rwnd:	TRUE: start the scan from the beginning of the file
		FALSE: start the scan from the present position
    wrap:	TRUE: allow the restart
		FALSE: quit after first EOF
  
  Returns: TRUE if the string is found, FALSE otherwise.
*/

#define MAXPATLEN 30

{ static struct state_struct { char nxtchr ;
			       struct state_struct *fail ; } states[MAXPATLEN] ;
  struct state_struct *state ;
  int patlen,patndx,statendx ;
  const char *patptr ;
  bool retry ; 
  chartab_struct *chr ;
  filblk_struct *filblk ;
  FILE *stream ;
  flags mode ;
  const char *rtnnme = "dyio_scan" ;
/*
  Check to make sure the stream ID is OK and lex points to a string.
*/
  if (INVALID_STREAMID(id))
  { errmsg(5,rtnnme,"stream id",id) ;
    return (FALSE) ; }
  filblk = &filblks[id] ;
  if (flgon(filblk->modes,io_active) == FALSE)
  { errmsg(15,rtnnme,id) ;
    return (FALSE) ; }
  if (flgon(filblk->modes,io_read) == FALSE)
  { errmsg(16,rtnnme,dyio_idtopath(id)) ;
    return (FALSE) ; }
  if (pattern == NULL)
  { errmsg(2,rtnnme,"pattern") ;
    return (FALSE) ; }
  patlen = strlen(pattern) ;
  if (patlen > MAXPATLEN)
  { errmsg(25,rtnnme,pattern,MAXPATLEN) ;
    patlen = MAXPATLEN ; }
/*
  Now set up the states array. For each state, the entry specifies the next
  character and the failure function value. We initialize the entries for
  states[0] and states[1], then enter the main loop to calculate the rest.
*/
  state = states ; patptr = pattern ;
  state->nxtchr = *patptr ;
  state->fail = states ;
  state++ ; patptr++ ;
  state->nxtchr = *patptr ;
  state->fail = states ;
  for (statendx = 2 ; statendx < patlen ; statendx++)
  { patndx = state->fail-states ;
    state++ ;
    while (*patptr != pattern[patndx] && patndx > 0)
      patndx = states[patndx].fail-states ;
    if (*patptr != pattern[patndx] && patndx == 0)
      state->fail = states ;
    else
      state->fail = &states[patndx+1] ;
    patptr++ ;
    state->nxtchr = *patptr ; }
/*
  Set up the control variables and position the file pointer.
*/
  stream = filblk->stream ;
  mode = filblk->modes ;
  if (rwnd == FALSE)
    retry = wrap ;
  else
  { rewind(stream) ;
    retry = FALSE ; }
/*
  Now begin the actual matching. The first action is to read in a character.
  i/o error causes immediate termination; nxtchar will already have printed an
  error message. EOF is handled by rewinding the file, if allowed, otherwise
  it too causes a return. If we rewind, state is forced to 0 and we abort this
  iteration of the loop.
*/
  for (state = states ; state != &states[patlen] ; )
  { chr = nxtchar(stream,mode) ;
    if (chr->class1 == CCERR) return (FALSE) ;
    if (chr->class1 == CCEOF)
    { if (retry == FALSE)
	return (FALSE) ;
      else
      { rewind(stream) ;
	retry = FALSE ;
	state = states ;
	continue ; } }
/*
  OK, we've got a character and we're actually ready to check it against the
  pattern. If it matches, we advance to the next state.
*/
    if (chr->chrtrn == state->nxtchr)
    { state++ ;
      continue ; }
/*
  It didn't match, so we apply the failure function until we get to a state
  where chr == state->nxtchr or state == 0. If we reach a state where nxtchr
  matches chr, we bump the state by one and iterate, otherwise we just iterate
  leaving state at 0.
*/
    for (state = state->fail ;
	 !(chr->chrtrn == state->nxtchr || state == states) ;
	 state = state->fail) ;
    if (chr->chrtrn == state->nxtchr) state++ ; }
/*
  To get here, we've matched. (Otherwise we would have returned due to EOF up
  above. Celebrate by returning TRUE.
*/
  return (TRUE) ; }

#undef MAXPATLEN



lex_struct *dyio_scanlex (ioid id)

/*
  This routine is a general lexical scanner. It recognizes three classes of 
  lexemes: numbers (DY_LCNUM), identifiers (DY_LCID), and delimiters (DY_LCDEL). As
  with character classes, these are augmented by two classes for EOF (DY_LCEOF)
  and i/o error (DY_LCERR). Lexemes are limited to 80 characters in length. Longer
  lexemes are truncated and the user is warned, but no error is indicated in
  the returned value.

  The basic form of the scanner follows the outline presented in Section 3.3
  of Gries, "Compiler Construction for Digital Computers".

  Parameters:
    id:	stream id from dyio_openfile.

  Returns: A lex_struct reflecting the result of the scan. This will be a
	   normal lex_struct or a lex_struct reflecting an exception (EOF or
	   error)
*/

#define MAXLEXLEN 80

{ chartab_struct *chr ;
  char *lexptr ;
  bool ovf ;
  filblk_struct *filblk ;
  FILE *stream ;
  flags mode ;
  static char stringspace[MAXLEXLEN+1] ;
  static lex_struct lex = {DY_LCNIL,stringspace} ;
  const char *rtnnme = "dyio_scanlex" ;
/*
  Check to make sure that the stream id is OK.
*/
  if (INVALID_STREAMID(id))
  { errmsg(5,rtnnme,"stream id",id) ;
    return (&lex_err) ; }
  filblk = &filblks[id] ;
  if (flgon(filblk->modes,io_active) == FALSE)
  { errmsg(15,rtnnme,id) ;
    return (&lex_err) ; }
  if (flgon(filblk->modes,io_read) == FALSE)
  { errmsg(16,rtnnme,dyio_idtopath(id)) ;
    return (&lex_err) ; }
  stream = filblk->stream ;
  mode = filblk->modes ;
/*
  Step over any characters of class CCINV (i.e., whitespace). If we get EOF
  or i/o error here, return appropriately.
*/
  for (chr = nxtchar(stream,mode) ;
       chr->class1 == CCINV ;
       chr = nxtchar(stream,mode)) ;
  if (chr->class1 == CCEOF) return (&lex_eof) ;
  if (chr->class1 == CCERR) return (&lex_err) ;
/*
  Form the lexeme. Note that the algorithm for forming a number does not
  attempt to be picky - it will handle a well-formed number, but it will also
  pass malformed numbers of various sorts.
*/
  lexptr = lex.string ;
  *lexptr++ = chr->chrtrn ;
  ovf = FALSE ;
  switch (chr->class1)
  { case CCDIG:				/* Form a number. */
    { for (chr = nxtchar(stream,mode) ; (1) ; chr = nxtchar(stream,mode))
      { if ((chr->class1 != CCDIG &&
	     chr->chrtrn != 'E' && chr->chrtrn != 'e') ||
	    ((chr->chrtrn == '+' || chr->chrtrn == '-') &&
	     *(lexptr-1) != 'E' && *(lexptr-1) != 'e'))
	  break ;
	if (ovf == FALSE)
	{ if (lexptr < &stringspace[MAXLEXLEN])
	  { *lexptr++ = chr->chrtrn ; }
	  else
	  { ovf = TRUE ;
	    *lexptr = '\0' ;
	    errmsg(26,rtnnme,lex.string,MAXLEXLEN) ; } } }
      if (!(chr->class1 == CCEOF || chr->class1 == CCERR))
	ungetc(chr->chrtrn,stream) ;
      if ((lexptr-1) == lex.string &&
	  (*(lexptr-1) == '+' || *(lexptr-1) == '-'))
	lex.class = DY_LCDEL ;
      else
	lex.class = DY_LCNUM ;
      break ; }
    case CCIDB:				/* Form an identifier. */
    { for (chr = nxtchar(stream,mode) ;
	   chr->class2 == CCIDC ;
	   chr = nxtchar(stream,mode))
	if (ovf == FALSE)
	{ if (lexptr < &stringspace[MAXLEXLEN])
	  { *lexptr++ = chr->chrtrn ; }
	  else
	  { ovf = TRUE ;
	    *lexptr = '\0' ;
	    errmsg(26,rtnnme,lex.string,MAXLEXLEN) ; } }
      if (!(chr->class1 == CCEOF || chr->class1 == CCERR))
	ungetc(chr->chrtrn,stream) ;
      lex.class = DY_LCID ;
      break ; }
    case CCDEL:				/* Form single-character delimiter. */
    { lex.class = DY_LCDEL ;
      break ; }
    default:				/* Should never execute this. */
    { errmsg(1,rtnnme,__LINE__) ;
      return (&lex_err) ; } }
/*
  Clean up and return. We check to see if the scan was terminated by an i/o
  error, and if so return lex_err. Otherwise, we make sure that the string ends
  with a NULL and return. Termination due to EOF isn't necessarily an error;
  an EOF return will be generated on the next call.
*/
  if (chr->class1 == CCERR) return (&lex_err) ;
  *lexptr = '\0' ;
  return (&lex) ; }

#undef MAXLEXLEN



lex_struct *dyio_scanstr (ioid id,
			  lexclass stype, int fslen, char qschr, char qechr)

/*
  This routine is a string scanner. It scans two types of strings, fixed length
  (DY_LCFS) and quoted (DY_LCQS). Before starting the string scan, the input is read
  until a non-CCINV character is encountered (i.e, we skip over whitespace
  characters).
  
  Fixed length strings have the length specified by fslen. The algorithm is
  to copy the next fslen characters into a text string and return.

  Quoted strings are bounded by a start character and an end character, passed
  in qschr and qechr, respectively.

  A check is made to see if the specified start character is present unless the
  start character is '\0', in which case no check is made, or the start
  character is ' ', in which case any CCINV (whitespace) character satisfies
  the check. It is an error if a specific start character is not matched in the
  input.
  
  Any single occurrence of the end character ends the string, unless the end
  character is '\0' or ' '. For '\0', any CCINV or CCDEL character ends the
  string, while for ' ' only CCINV characters will end the string. Double
  occurrences of printing end characters are reduced to a single occurrence in
  the string returned; surrounding quote characters are stripped. For
  non-printing end characters, the double-occurrence escape is not available.
  (This avoids problems with forming strings ended by \n or \f from terminal
  input, where the read may block because there are no more characters.)
 
  NOTE: It is possible to parse a string of length 0 (and hence get a DY_LCNIL
	return type) when a quoted string is parsed with start character '\0'.
	Example: a request for a string quoted with '\0','\n' (to parse off
	entire lines) will return DY_LCNIL if the line contains only a '\n'.

	Further, a string quoted with '\0','\0' is satisfied by any amount of
	whitespace, even the whitespace skipped over at the start of the
	parse.

	Finally, to avoid leaking the buffer allocated for the last non-null
	string, use a call of dyio_scanstr(ioid,DY_LCQS,0,'\0','\0') and free
	lex.string iff lex.class == DY_LCQS.

  Parameters:
    id:		i/o ID from dyio_openfile.
    stype:	String type, allowable values are DY_LCFS and DY_LCQS.
    fslen:	String length (fixed-length strings only).
    qschr:	Start character (quoted strings only).
    qechr:	End character (quoted strings only).
  
  Returns: A lex_struct reflecting the result of the scan. This will be a
	   normal lex_struct or a lex_struct reflecting an exception (EOF or
	   error)
*/

#define MAXQSLEN 250

{ chartab_struct *chr ;
  char *lexptr ;
  bool invseen ;
  int ndx,expcnt ;
  filblk_struct *filblk ;
  FILE *stream ;
  flags mode ;
  static lex_struct lex = {DY_LCNIL,NULL} ;
  const char *rtnnme = "dyio_scanstr" ;
/*
  Release the previous lex space, if necessary.
*/
  if (lex.string != NULL)
  { FREE(lex.string) ;
    lex.string = NULL ; }
/*
  Check to make sure that the stream id is OK.
*/
  if (INVALID_IOID(id))
  { errmsg(5,rtnnme,"stream id",id) ;
    return (&lex_err) ; }
  filblk = &filblks[id] ;
  if (flgon(filblk->modes,io_active) == FALSE)
  { errmsg(15,rtnnme,id) ;
    return (&lex_err) ; }
  if (flgon(filblk->modes,io_read) == FALSE)
  { errmsg(16,rtnnme,dyio_idtopath(id)) ;
    return (&lex_err) ; }
  stream = filblk->stream ;
  mode = filblk->modes ;
/*
  Step over any characters of class CCINV (i.e., whitespace). If we get EOF
  or i/o error here, return appropriately.
*/
  invseen = FALSE ;
  for (chr = nxtchar(stream,mode) ;
       chr->class1 == CCINV ;
       chr = nxtchar(stream,mode))
    invseen = TRUE ;
  if (chr->class1 == CCERR) return (&lex_err) ;
  if (chr->class1 == CCEOF)
  { if (invseen == TRUE && qschr == '\0' && (qechr == '\0' || qechr == ' '))
      return (&lex_nil) ;
    else
      return (&lex_eof) ; }
/*
  Now begin to construct the string.
*/
  switch (stype)
  {
/*
  Fixed length strings. Acquire some space of the right size, then copy the
  required number of characters.
*/
    case DY_LCFS:
    { if (fslen <= 0)
      { errmsg(5,rtnnme,"fslen",fslen) ;
	ungetc(chr->chrtrn,stream) ;
	return (&lex_err) ; }
      lexptr = (char *) MALLOC(fslen+1) ;
      lex.string = lexptr ;
      *lexptr++ = chr->chrtrn ;
      for (ndx = fslen-1 ; ndx != 0 ; ndx--)
      { chr = nxtchar(stream,mode) ;
	if (chr->class1 == CCEOF || chr->class1 == CCERR)
	{ *lexptr = '\0' ;
	  errmsg(27,rtnnme,fslen,lex.string) ;
	  FREE(lex.string) ;
	  lex.string = NULL ;
	  return (&lex_err) ; }
	*lexptr++ = chr->chrtrn ; }
      lex.class = DY_LCFS ;
      break ; }
/*
  Quoted strings. We'll guess that MAXQSLEN characters will normally hold us,
  and expand if necessary. The first switch looks to see if we've satisfied the
  start character conditions. Then we allocate space and proceed to scan off
  the string. We allocate MAXQSLEN+1 to allow space for a closing null.
*/
    case DY_LCQS:
    { switch (qschr)
      { case '\0':
	{ if (invseen == TRUE && (qechr == '\0' || qechr == ' '))
	  { ungetc(chr->chrtrn,stream) ;
	    return (&lex_nil) ; }
	  break ; }
	case ' ':
	{ if (invseen != TRUE)
	  { ungetc(chr->chrtrn,stream) ;
	    return (&lex_err) ; }
	  break ; }
	default:
	{ if (qschr != chr->chrtrn)
	  { ungetc(chr->chrtrn,stream) ;
	    return (&lex_err) ; }
	  chr = nxtchar(stream,mode) ;
	  break ; } }
      lexptr = (char *) MALLOC(MAXQSLEN+1) ;
      lex.string = lexptr ;
      ndx = MAXQSLEN ;
      expcnt = 0 ;
      while (1)
      { if (chr->class1 == CCEOF || chr->class1 == CCERR)
	{ *lexptr = '\0' ;
	  errmsg(28,rtnnme,qschr,qechr,lex.string) ;
	  FREE(lex.string) ;
	  lex.string = NULL ;
	  return (&lex_err) ; }
	if ((qechr == ' ' && chr->class1 == CCINV) ||
	    (qechr == '\0' && (chr->class1 == CCINV || chr->class1 == CCDEL)))
	  break ;
	if (chr->chrtrn == qechr)
	{ if (qechr >= ' ')
	  { chr = nxtchar(stream,mode) ;
	    if (chr->chrtrn != qechr)
	    { if (chr->class1 != CCEOF && chr->class1 != CCERR)
		ungetc(chr->chrtrn,stream) ;
	      break ; } }
	  else
	    break ; }
	*lexptr++ = chr->chrtrn ;
	ndx-- ;
	if (ndx <= 0)
	{ expcnt++ ;
	  lexptr = (char *) MALLOC((expcnt+1)*MAXQSLEN+1) ;
	  (void) strncpy(lexptr,lex.string,expcnt*MAXQSLEN) ;
	  FREE(lex.string) ;
	  lex.string = lexptr ;
	  lexptr += expcnt*MAXQSLEN ;
	  ndx = MAXQSLEN ; }
	chr = nxtchar(stream,mode) ; }
      if (lex.string == lexptr)
      { FREE(lex.string) ;
	lex.string = NULL ;
	return (&lex_nil) ; }
      else
	lex.class = DY_LCQS ;
      break ; }
    default:
    { errmsg(5,rtnnme,"string type",stype) ;
      ungetc(chr->chrtrn,stream) ;
      return (&lex_err) ; } }
/*
  That's it. Terminate the string and  return. Note that this return is used
  only when a non-null string has been scanned.
*/
  *lexptr = '\0' ;
  return (&lex) ; } 

#undef MAXQSLEN



/*
  The following output routines provide some measure of isolation for output,
  but their major function in life is to localize echoing control.
*/

void dyio_flushio (ioid id, bool echo)

/*
  This routine flushes buffered output. It is strictly an interface routine,
  to hide the translation from integer stream id to the stdio stream handle
  and localize echoing control.

  If id == 0 and echo == FALSE, nothing happens.
  If id == 0 and echo == TRUE, only stdout is affected.
                  (etc.)

  Parameters:
    id:		stream id from dyio_openfile
    echo:	specifies if we're echoing to stdout, and hence should flush
		it

  Returns: undefined
*/

{ FILE *strmid ;
  filblk_struct *filblk ;
  const char *rtnnme = "dyio_flushio" ;

/*
  Check for sane id.
*/
  if (INVALID_IOID(id))
  { errmsg(5,rtnnme,"i/o id",id) ;
    return ; }
/*
  Flush the stream specified by id, if it's active and writeable.
*/
  if (id != IOID_NOSTRM)
  { filblk = &filblks[id] ;
    if (flgon(filblk->modes,io_active) == FALSE)
    { errmsg(15,rtnnme,id) ; }
    else
    if (flgon(filblk->modes,io_write) == FALSE)
    { errmsg(17,rtnnme,dyio_idtopath(id)) ; }
    else
    { strmid = filblk->stream ;
      if (fflush(strmid) != 0) perror(rtnnme) ; } }
/*
  Flush stdout if echo is TRUE. stdout is assumed active and writeable.
*/
  if (echo == TRUE)
    if (fflush(stdout) != 0) perror(rtnnme) ;

  return ; } 



void dyio_outfmt (ioid id, bool echo, const char *pattern, ... )

/*
  This routine provides the client with an interface to the vfprintf routine,
  while maintaining a logging facility.

  The form of the call should be
    dyio_outfmt(id,echo,pattern,parm1, ... ,parmn)

  If id == 0 and echo == FALSE, nothing happens.
  If id == 0 and echo == TRUE, only stdout is affected.
                    (etc.)

  Parameters:
    id:		stream id from dyio_openfile
    echo:	TRUE to echo to user's tty, FALSE otherwise
    pattern:	output pattern, passed to printf
    parmi:	the parameters for the pattern
  
  Returns: undefined
*/

{ va_list parms ;

  filblk_struct *filblk ;
  const char *rtnnme = "dyio_outfmt" ;
/*
  Sanity checks on i/o id and pattern.
*/
  if (INVALID_IOID(id))
  { errmsg(5,rtnnme,"i/o id",id) ;
    return ; }
  if (pattern == NULL)
  { errmsg(2,rtnnme,"pattern") ;
    return ; }
/*
  Print to stream id, assuming it's active and writeable.
*/
  va_start(parms,pattern) ;

  if (id != IOID_NOSTRM)
  { filblk = &filblks[id] ;
    if (flgon(filblk->modes,io_active) == FALSE)
    { errmsg(15,rtnnme,id) ; }
    else
    if (flgon(filblk->modes,io_write) == FALSE)
    { errmsg(17,rtnnme,dyio_idtopath(id)) ; }
    else
    { (void) vfprintf(filblk->stream,pattern,parms) ; } }
/*
  Print to stdout if echo is TRUE and the stream id is not stdout.
  stdout is assumed active and writeable.
*/
  if (echo == TRUE && id != dyio_pathtoid("stdout",NULL)) 
    (void) vfprintf(stdout,pattern,parms) ;

  va_end(parms) ;

  return ; }



void dyio_outchr (ioid id, bool echo, char chr)

/*
  This routine interfaces to the fputc function. Note that it traps chr == '\0'
  as an error.

  If id == 0 and echo == FALSE, nothing happens.
  If id == 0 and echo == TRUE, only stdout is affected.

  Parameters:
    id:		stream id from dyio_openfile
    echo:	TRUE to echo to user's tty, FALSE otherwise
    chr:	character to be output

  Returns: undefined
*/

{ FILE *strmid ;
  filblk_struct *filblk ;
  const char *rtnnme = "dyio_outchr" ;

/*
  Sanity checks on the parameters.
*/
  if (INVALID_IOID(id))
  { errmsg(5,rtnnme,"i/o id",id) ;
    return ; }
  if (chr == '\0')
  { errmsg(2,rtnnme,"chr") ;
    return ; }
/*
  Write the character to the stream, if it's active and writeable.
*/
  if (id != IOID_NOSTRM)
  { filblk = &filblks[id] ;
    if (flgon(filblk->modes,io_active) == FALSE)
    { errmsg(15,rtnnme,id) ; }
    else
    if (flgon(filblk->modes,io_write) == FALSE)
    { errmsg(17,rtnnme,dyio_idtopath(id)) ; }
    else
    { strmid = filblk->stream ;
      putc(chr,filblks[id].stream) ; } }
/*
  Write the character to stdout, if echo is TRUE and the stream id is not
  stdout. stdout is assumed to be active and writeable.
*/
    if (echo == TRUE && id != dyio_pathtoid("stdout",NULL)) putc(chr,stdout) ;

  return ; }



int dyio_outfxd (char *buffer, int fldsze, char lcr, const char *pattern, ... )

/*
  The form of the call should be
    dyio_outfxd (buffer,fldsze,lcr,pattern,parm1, ... parmn)

  The input to this routine is a varargs list comprising the name of a buffer,
  a field size, an 'lcr' (left/centre/right) character, and finally a pattern
  and parameter list to be used by vsprintf. 'lcr' is an indication of whether 
  the output is to be left justified, centered, or right justified.  Strings 
  longer than the specified fieldsize are truncated. The maximum fieldsize is 
  limited (arbitrarily) to FLDMAX characters, since we have to keep a hidden 
  buffer (FLDMAX is 240 currently). Centering is the only real work here. 
  Spaces is an array of FLDMAX spaces. A negative fldsze for left justification 
  is interpreted to mean "don't pad, just truncate if necessary".

  Parameters:
    buffer:	buffer area used to construct the field
    fldsze:	desired size of the output field
    lcr:	left-justified ('l'), centered ('c'), or right-justified ('r')
    pattern:	output pattern, passed to printf
    parms:	the parameters for the pattern
  
  Returns: the length of the string placed in the buffer; 0 is returned in the
	   event of an error
*/

#define FLDMAX 240

{ int count ;
  va_list parms ;

  static char ourbuf[FLDMAX+1] ;
  static char spaces[FLDMAX] ;
  int tlen,hlen ;
  bool pad ;
  static bool frstflg = TRUE ;
  const char *rtnnme = "dyio_outfxd" ;

/*
  Load the spaces array just once.
*/
  if (frstflg == TRUE)
  for (count = 0 ; count < FLDMAX ; count++)
  { *(spaces+count) = ' ' ;
    frstflg = FALSE ; }
/*
  Do some sanity checks.
*/
  if (buffer == NULL)
  { errmsg(2,rtnnme,"buffer") ;
    return (0) ; }
  pad = (fldsze < 0)? FALSE : TRUE ;
  fldsze = abs(fldsze) ;
  if (fldsze < 1 || fldsze > FLDMAX)
  { errmsg(5,rtnnme,"fldsze",fldsze) ;
    return (0) ; }
  if (lcr != 'l' && lcr != 'c' && lcr != 'r')
  { errmsg(3,rtnnme,"left/center/right",lcr) ;
    return (0) ; }
  if (pattern == NULL)
  { errmsg(2,rtnnme,"pattern") ;
    return (0) ; }
/*
  vsprintf builds the user's string in ourbuf. Then use fldsze and lcr to
  construct the string we return to the user.
*/
  va_start(parms,pattern) ;
  (void) vsprintf (ourbuf,pattern,parms) ;
  va_end(parms) ;

  tlen = strlen(ourbuf) ;
  switch (lcr)
  { case 'l':
    { if (fldsze > tlen)
      { memcpy((void *) buffer,(void *) ourbuf,tlen) ;
	if (pad == TRUE)
	  memcpy((void *) (buffer+tlen),(void *) spaces,fldsze-tlen) ;
	else
	  fldsze = tlen ; }
      else
      { memcpy((void *) buffer,(void *) ourbuf,fldsze) ; }
      break ; }
    case 'c':
    { if (fldsze > tlen)
      { hlen = fldsze-tlen ;
	memcpy((void *) buffer,(void *) spaces,hlen>>1) ;
	memcpy((void *) (buffer+(hlen>>1)),(void *) ourbuf,tlen) ;
	memcpy((void *) (buffer+(hlen>>1)+tlen),
	      (void *) spaces,hlen-(hlen>>1)) ; }
      else
      { memcpy((void *) buffer,(void *) ourbuf,fldsze) ; }
      break ; }
    case 'r':
    { if (fldsze > tlen)
      { memcpy((void *) buffer,(void *) spaces,fldsze-tlen) ;
	memcpy((void *) (buffer+fldsze-tlen),(void *) ourbuf,tlen) ; }
      else
      { memcpy((void *) buffer,(void *) ourbuf,fldsze) ; }
      break ; } }
  buffer[fldsze] = '\0' ;
  return (fldsze) ; }

#undef FLDMAX



#ifdef _DYLIB_FORTRAN

void dyio_outfmt_ (integer *ftnid, logical *ftnecho, char *pattern, ... )

/*
  This routine provides a Fortran client with an interface to the vfprintf
  routine of ANSI C and the logging facilities of this i/o library. It deals
  with translating the argument types supplied by Fortran-to-C interface
  conventions into the arguments required by vprintf. The method is to
  construct a new varargs block, which is handed over to vprintf. To do
  this, we have to make the fragile assumption that a varargs block is
  constructed in a straightforward manner --- as data items written into
  a contiguous block of storage which we can allocate. There are more
  extensive comments with errmsg_ and warn_ in errs.c.

  The call over in Fortran will look like
    dyio_outfmt(id,echo,pattern,ftnargtype1,arg1, ... ,
		ftnargtypen,argn,ftnargEND)

  The routine deals with a limited set of argument types: the Fortran types
  integer, double_precision, and character, and the special categories of
  variable and constraint names (the last two historical artifacts from
  the initial impetus for this routine -- the Fortran ylp library -- now
  replaced by dylp).  A special type code is used to indicate the end of
  the list of <type,arg> pairs.

  If id == 0 and echo == FALSE, nothing happens.
  If id == 0 and echo == TRUE, only stdout is affected.
                               (etc.)

  Parameters:
    ftnid:	the i/o id
    ftnecho:	true to echo output to stdout
    pattern:	output pattern, passed to printf
    argtype, arg: the parameters for the pattern
  
  Returns: undefined
*/

{ va_list fargs,varargp ;
  ioid id ;
  bool echo ;
  filblk_struct *filblk ;
  int type;

  double varargs[64] ;			/* double avoids alignment problems */
  int intarg ;
  double dblarg ;
  char *chararg ;

  const char *rtnnme = "dyio_outfmt_" ;

/*
  Convert the stream id and echo values.
*/
  id = (ioid) *ftnid ;
  echo = (*ftnecho == TRUEL)?TRUE:FALSE ;
/*
  Sanity checks on stream id and pattern.
*/
  if (INVALID_IOID(id))
  { errmsg(5,rtnnme,"i/o id",id) ;
    return ; }
  if (id != IOID_NOSTRM)
  { filblk = &filblks[id] ;
    if (flgon(filblk->modes,io_active) == FALSE)
    { errmsg(15,rtnnme,id) ;
      return ; }
    if (flgon(filblk->modes,io_write) == FALSE)
    { errmsg(17,rtnnme,dyio_idtopath(id)) ;
      return ; } }
  if (pattern == NULL)
  { errmsg(2,rtnnme,"pattern") ;
    return ; }
/*
  Now start up a loop to process the remainder of the arguments. For each
  <type,arg> pair, we pull off the type and use it to condition a switch
  with a case for each type of argument we're prepared to deal with.
*/
  varargp = (va_list) &varargs[0] ;
  va_start(fargs,pattern) ;
  for (type = (int) *va_arg(fargs,integer *) ;
       type != ftnargEND ;
       type = (int) *va_arg(fargs,integer *))
    switch (type)
    { case ftnargINTEGER:
      { intarg = (int) *va_arg(fargs,integer *) ;
        *((int *) varargp) = intarg ;
	if (intarg != va_arg(varargp,int))
	{ errmsg(1,rtnnme,__LINE__) ;
	  return ; }
	break ; }
      case ftnargDOUBLE_PRECISION:
      { dblarg = (double) *va_arg(fargs,double_precision *) ;
	memcpy(varargp,&dblarg,sizeof(double)) ;
        if (dblarg != va_arg(varargp,double))
	{ errmsg(1,rtnnme,__LINE__) ;
	  return ; }
	break ; }
      case ftnargCHARACTER:
      { chararg = va_arg(fargs,char *) ;
        *((char **) varargp) = chararg ;
	if (chararg != va_arg(varargp,char *))
	{ errmsg(1,rtnnme,__LINE__) ;
	  return ; }
	break ; }
      default:
      { errmsg(7,rtnnme,__LINE__,"Fortran argument type code",type) ;
	return ; } }
  va_end(fargs) ;
/*
  Finally, do the printing. Print to stdout if echo is TRUE. stdout is assumed
  active and writeable.
*/
  if (id != IOID_NOSTRM)
  { (void) vfprintf(filblk->stream,pattern,((va_list) &varargs[0])) ; }
  if (echo == TRUE)
  { (void) vfprintf(stdout,pattern,((va_list) &varargs[0])) ; }

  return ; }

#endif /* _DYLIB_FORTRAN */

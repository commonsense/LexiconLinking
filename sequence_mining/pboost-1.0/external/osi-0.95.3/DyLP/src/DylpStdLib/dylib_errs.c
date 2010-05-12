/*
  This file is part of the support library  for the OsiDylp LP distribution.

        Copyright (C) 2005 Lou Hafer

        School of Computing Science
        Simon Fraser University
        Burnaby, B.C., V5A 1S6, Canada
        lou@cs.sfu.ca

  This program is free software; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the Free
  Software Foundation; either version 2 of the License, or (at your option)
  any later version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
  more details.

  You should have received a copy of the GNU General Public License along
  with this program; if not, write to the Free Software Foundation, Inc.,
  51 Franklin St., Fifth Floor, Boston, MA  02110-1301  USA
*/

/*
  This is an error reporting/logging package. It contains routines to print
  formatted error and warning messages. The client passes a numeric
  identifier, along with any necessary parameters. The routines will look up
  the identifier in a file to find the associated error message and pass the
  message, with parameters, to printf.

  If the conditional compilation variable DYLP_NDEBUG is defined, the warn routine
  is converted to a complete noop.

  The error message file may be passed as a parameter to errinit, or specified
  using the environment variable ERRMSGTXT. Failing anything else, it will
  default to "errmsg.txt".

  In addition to reporting error messages via stderr, the package has
  provision for optionally logging error messages into a separate file,
  specified by elogchn.

  The conditional compilation variable _DYLIB_FORTRAN, when defined, will
  trigger the inclusion of code which provides an errmsg subroutine and
  interface usable from Fortran code.
*/

#include "dylib_std.h"

static char sccsid[] UNUSED = "@(#)errs.c	3.13	11/11/04" ;
static char svnid[] UNUSED = "$Id: dylib_errs.c 71 2006-06-09 04:21:15Z andreasw $" ;

#include <stdio.h>
#include <stdarg.h>

#ifndef MALLOC
#define MALLOC(zz_sze_zz) malloc(zz_sze_zz)
#endif

/*
  emsgchn is the channel used to access the error message text file; elogchn is
  the channel used to access an error log file. errecho controls whether error
  messages go to stderr.
*/

static FILE *emsgchn,*elogchn ;
static char *elogname ;
static bool errecho ;

/*
  The text of an error message (the unexpanded pattern) can be no longer than
  MAXERRTXT.
*/

#define MAXERRTXT 240

#define ERRPFXLEN  13
static char errtxt[MAXERRTXT+ERRPFXLEN] = "\n%s (error): " ;

#ifndef DYLP_NDEBUG
#define WARNPFXLEN 15
static char warntxt[MAXERRTXT+WARNPFXLEN] = "\n%s (warning): " ;
#endif


/*
  Only errmsg_, warn_, and a small bit of initialisation in errinit need the
  declarations of fortran.h. The file xmp.h is peculiar to the bonsai code
  with Fortran ylp and la05 libraries.
*/

#ifdef _DYLIB_FORTRAN
#  include "fortran.h"
#endif




void errinit (const char *emsgpath, char *elogpath, bool echo)

/*
  This routine initializes the error reporting facilities. It attempts to
  open the error message file and an error logging file.

  Parameters:
    emsgpath:	path name for the error message file (defaults to errmsg.txt)
    elogpath:	path name for logging error messages (the default is to not log
		error messages)
    echo:	TRUE if error messages should be echoed to stderr, FALSE
		otherwise.

  Returns: undefined
*/

{ const char *rtnnme = "errinit" ;

/*
  See if we can open a channel to the error message source text file. If we
  can't, we'll only be able to report routine names and error numbers.
*/
  if (emsgpath == NULL)
  { emsgpath = getenv("ERRMSGTXT") ;
    if (emsgpath == NULL) emsgpath = "errmsg.txt" ; }
  emsgchn = fopen(emsgpath,"r") ;
  if (emsgchn == NULL)
  { fprintf(stderr,"\n%s: couldn't open error message text file \"%s\".\n",
	    rtnnme,emsgpath) ;
    perror(rtnnme) ;
    fprintf(stderr,"\n%s: only numeric error codes will be reported.\n",
	    rtnnme) ; }
/*
  If a log file name has been supplied, try to open it. Keep the name locally
  for the benefit of errlogq, should it be called.
*/
  if (elogpath == NULL)
  { elogchn = NULL ;
    elogname = NULL ; }
  else
  { elogname = (char *) MALLOC(strlen(elogpath)+1) ;
    strcpy(elogname,elogpath) ;
    elogchn = fopen(elogname,"w") ;
    if (elogchn == NULL)
    { fprintf(stderr,"\n%s: couldn't open error logging file \"%s\".\n",
	      rtnnme,elogname) ;
      perror(rtnnme) ; } }
  
  errecho = echo ;

#ifdef _DYLIB_FORTRAN
/*
  Initialise a common block with the type codes used by Fortran when it calls
  errmsg_.
*/
  argcod_.integer_code = ftnargINTEGER ;
  argcod_.double_precision_code = ftnargDOUBLE_PRECISION ;
  argcod_.character_code = ftnargCHARACTER ;
  argcod_.varname_code = ftnargVARNAME ;
  argcod_.conname_code = ftnargCONNAME ;
  argcod_.end_code = ftnargEND ;
#endif

  return ; }


void errterm (void)

/*
  This routine cleans up data structures owned by the error package.

  Parameters: none

  Returns: undefined.
*/

{ if (elogname != NULL) FREE(elogname) ;

  return ; }



/*
  These next two routines are exported only to io.c in the default setup.
  They are used to allow the user to control the error log files through
  io.c
*/

FILE *errlogq (char **elogpath)

/*
  This routine is here to keep the compartmentalisation of errs.c. Its sole
  purpose in life is to return the stdio FILE handle and path name for the
  error log file.

  Parameters:
    elogpath	(o) Will be set to point to the name of the error log file,
		    or NULL if no error log file is specified.

  Returns: the file handle, which will be NULL if no error log file exists.
*/

{ *elogpath = elogname ;
  return (elogchn) ; }


bool reseterrlogchn (const char *newpath, FILE *newchn, bool echo, bool close)

/*
  Another routine solely dedicated to compartmentalisation. If newpath is
  non-NULL, the routine will reset the log file, opening it if newchn is NULL,
  otherwise using newchn. The errecho flag is set to echo.

  Parameters:
    newpath:	new path for the log file; ignored if NULL
    newchn:	new handle for the log file channel; ignored if NULL
    echo:	TRUE if errors should be echoed to stderr, FALSE otherwise.
    close:	TRUE of the old file should be closed, FALSE otherwise.
  
  Returns: TRUE if the new file is opened successfully, FALSE otherwise
	   (note that failure to close the old file still gets a TRUE return)
*/

{ bool success ;
  const char *rtnnme = "reseterrlogchn" ;

/*
  If newchn is NULL, try to open the specified path. Retain the previous
  settings if we can't open the new file. If newchn isn't NULL, use it
  without question. There's the possibility the user has called us only to
  change the echo, in which case newpath will be NULL.
*/
  success = TRUE ;
  if (newpath != NULL)
  { if (newchn == NULL)
    { newchn = fopen(newpath,"w") ;
      if (newchn == NULL)
      { fprintf(stderr,"\n%s: couldn't open error logging file \"%s\".\n",
		rtnnme,newpath) ;
	perror(rtnnme) ;
	fprintf(stderr,"\n%s: retaining previous file \"%s\".\n",
		rtnnme,elogname) ;
	success = FALSE ; } }

    if (success == TRUE)
    { if (close == TRUE)
      { if (fclose(elogchn) == EOF)
	{ fprintf(stderr,
		  "\n%s: couldn't close previous error logging file \"%s\".\n",
		  rtnnme,elogname) ;
	  perror(rtnnme) ; } }
      elogchn = newchn ;
      if (elogname != NULL) FREE(elogname) ;
      elogname = (char *) MALLOC(strlen(newpath)+1) ;
      strcpy(elogname,newpath) ; } }

  errecho = echo ;

  return (success) ; }



static char *finderrmsg (int errid, char *buffer)

/*
  This is a utility routine to search the error message file for the error
  message specified by errid.

  The format of an error message is:

    @<errid>@<message text>@

  A line which has '@' as the first character is taken as the start of an
  error message. <errid> can be of the form <number> or <number>:<number> and
  must be on one line. A newline immediately following the '@' which ends the
  error id line is stripped.  Otherwise the message is free format, with "\"
  used to escape an '@' in the message text.  The message text can be any
  number of lines.

  Comments can be inserted between error messages.

  The form <number>:<number> is intended as a convenience for leaving blocks
  of error message ids open for future use. The accompanying message text is
  used for all error numbers in the range.

  Parameter:
    errid:	error message number
    buffer:	(o) buffer to hold the message text; assumed to be of length
		MAXERRTXT or larger

  Returns: the error message, or NULL
*/

{ int argcnt,chrcnt,id,id2,chr ;
  bool sync ;
  char *txtptr ;
  const char *rtnnme = "finderrmsg" ;

/*
  Rewind to get in sync. Then scan for the '@' that introduces the message
  number. When it's found, scan the error id. If it's the right message,
  break out of the loop. If it's the wrong message, just scan till the '@'
  that closes the text, then iterate the search loop. EOF causes a break from
  the loop and a return. Anything else unexpected in the scan causes a sync
`  error and a return.`,

*/
  rewind(emsgchn) ;
  id = 0 ;
  sync = FALSE ;
  while (fgets(buffer,MAXERRTXT,emsgchn) != NULL)
  { if (buffer[0] != '@') continue ;
    argcnt = sscanf(buffer,"@%d%n:%d%n",&id,&chrcnt,&id2,&chrcnt) ;
    if ((argcnt == 1 || argcnt == 2) && buffer[chrcnt] == '@')
    { if ((argcnt == 1 && id == errid) ||
	  (argcnt == 2 && id <= errid && errid <= id2))
      { sync = TRUE ;
	break ; } }
    else
    { fprintf(stderr,"\n%s: bad error message id format; line is:\n%s\n",
	      rtnnme,buffer) ;
      fprintf(stderr,"\tskipping to start of next message.\n") ; } }

  if (sync == FALSE)
  { if (feof(emsgchn))
    { fprintf(stderr,"\n%s: couldn't find error message %d.\n",rtnnme,errid) ; }
    else if (ferror(emsgchn))
    { fprintf(stderr,"\n%s: i/o error.\n",rtnnme) ;
      perror(rtnnme) ; }
    else
    { fprintf(stderr,"\n%s: internal confusion at line %d.\n",
	      rtnnme,__LINE__) ; }
    return (NULL) ; }
/*
  We've found the message. Shift any remaining message text on this line up
  into the front of the buffer.  If the end of the line is '@\n', and the @
  isn't preceded by '\', write a null and return.  Otherwise, start a loop to
  pick up any further lines of text. There's still a chance for EOF errors.
*/
  txtptr = &buffer[chrcnt+1] ;
  if (*txtptr != '\n')
  { chrcnt = 0 ;
    while (*txtptr != '\0')
    { chr = *txtptr++ ;
      if (chr == '\\') chr = *txtptr++ ;
      buffer[chrcnt++] = chr ; }
    buffer[chrcnt] = '\0' ;
    txtptr = &buffer[strlen(buffer)-1] ;
    if (*txtptr == '\n' && *(txtptr-1) == '@' && *(txtptr-2) != '\\')
    { *(txtptr-1) = '\0' ;
      return (buffer) ; }
    else
      txtptr++ ; }
  else
  { txtptr = buffer ; }
  for (chr = getc(emsgchn) ; chr != EOF ; chr = getc(emsgchn))
  { if (chr == '@') break ;
    if (chr == '\\')
    { chr = getc(emsgchn) ;
      if (chr == EOF)
      { fprintf(stderr,
		"\n%s: sync error - EOF following \"\\\" in message %d.\n",
		rtnnme,errid) ;
	return (NULL) ; } }
    *txtptr++ = chr ; }
  if (chr == EOF)
  { fprintf(stderr,"\n%s: sync error - EOF collecting text for message %d.\n",
	    rtnnme,errid) ;
    return (NULL) ; }
  *txtptr = '\0' ;

  return (buffer) ; }



void errmsg (int errid, ... )

/*
  Actual Call:
    errmsg(errid,ident,arg1, ... ,argn)

  errmsg is the error reporting function used by clients of this package.  It
  prints a message with the format "ident: error message".  errid identifies
  a printf-style error message template which is looked up in the error
  message text file.  A variable number of parameters may follow ident; they
  are assumed to be compatible with the way vfprintf will interpret the error
  message template. The message is printed to stderr and to the error log
  file, if it is open. In the event that emsgchn indicates the error message
  file wasn't found, only the error number and ident are printed.

  Parameters:
    errid:	number of the generic error message
    ident:	string used to prefix the error message (usually a routine
		name or similar identifier)
    arg1,...:	the arguments to be spliced into the error message template

  Returns: undefined
*/

{ char *ident ;
  va_list varargs ; 

/*
  Flush stdout and elogchn, so that any buffered normal output appears before
  the error message.
*/
  fflush(stdout) ;
  if (elogchn != NULL) fflush(elogchn) ;
/*
  If there's no error message file, or we can't find the message, just print
  the error number. Otherwise call vfprintf to handle expanding the message
  text.
*/
  va_start(varargs,errid) ; 
  if (emsgchn == NULL || finderrmsg(errid,&errtxt[ERRPFXLEN]) == NULL)
  { ident = va_arg(varargs,char *) ;
    if (errecho == TRUE)
      fprintf(stderr,"\n%s: error %d.\n",ident,errid) ;
    if (elogchn != NULL)
      fprintf(elogchn,"\n%s: error %d.\n",ident,errid) ; }
  else
  { if (errecho == TRUE)
    { vfprintf(stderr,errtxt,varargs) ;
      putc('\n',stderr) ; }
    if (elogchn != NULL) 
    { vfprintf(elogchn,errtxt,varargs) ;
      putc('\n',elogchn) ; } }
  va_end(varargs) ;
/*
  Flush the logged error message, so that the user will definitely see it.
*/
  if (elogchn != NULL) fflush(elogchn) ;

  return ; }



#ifdef DYLP_NDEBUG

void warn (int errid, ... )

{ return ; }

#else

void warn (int errid, ... )

/*
  Actual Call:
    warn(errid,ident,arg1, ... ,argn)

  Warn is functionally identical to errmsg, but prints the message
  "ident (warning): warning message".

  Parameters:
    errid:	number of the generic error message
    ident:	string used to prefix the error message (usually a routine
		name or similar identifier)
    arg1,...:	the arguments to be spliced into the error message template

  Returns: undefined
*/

{ char *ident ;
  va_list varargs ;

/*
  Flush stdout and elogchn, so that any buffered normal output appears before
  the error message.
*/
  fflush(stdout) ;
  if (elogchn != NULL) fflush(elogchn) ;
/*
  If there's no error message file, or we can't find the message, just print
  the error number. Otherwise call vfprintf to handle expanding the message
  text.
*/
  va_start(varargs,errid) ; 
  if (emsgchn == NULL || finderrmsg(errid,&warntxt[WARNPFXLEN]) == NULL)
  { ident = va_arg(varargs,char *) ;
    if (errecho == TRUE)
      fprintf(stderr,"\n%s: error %d.\n",ident,errid) ;
    if (elogchn != NULL)
      fprintf(elogchn,"\n%s: error %d.\n",ident,errid) ; }
  else
  { if (errecho == TRUE)
    { vfprintf(stderr,warntxt,varargs) ;
      putc('\n',stderr) ; }
    if (elogchn != NULL) 
    { vfprintf(elogchn,warntxt,varargs) ;
      putc('\n',elogchn) ; } }
  va_end(varargs) ;
/*
  Flush the logged error message, so that the user will definitely see it.
*/
  if (elogchn != NULL) fflush(elogchn) ;

  return ; }

#endif /* DYLP_NDEBUG */

#ifdef _DYLIB_FORTRAN



/*
  The two routines which follow, errmsg_ and warn_, are intended to be
  called from Fortran code. (The is also a companion routine, io.c:outfmt_.)
  This interface was created because early versions of the bonsai MILP
  code, written in C, used Fortran versions of the ylp and la05 subroutine
  libraries.

  The problem is that the things that Fortran will hand over cannot be passed
  directly to vfprintf - it passes pointers, vfprintf wants values for some
  things, pointers for others. This means we have to interpret the values
  passed from Fortran, writing our own varargs block to hand to vfprintf.
  To do this, we have to make the fragile assumption that a varargs block
  is constructed in a straightforward manner --- as data items written
  into a contiguous block of storage which we can allocate. Older versions
  of the varargs macros actually did this. With more recent versions,
  the assumption has become increasingly tenuous. As of GCC 2.96, gcc
  declares the varargs macros as builtins, with no hint of implementation.

  This hack has worked for some 10 years now, surviving ports to various
  architectures and the transition to ANSI C (varargs.h -> stdarg.h), but I'm
  starting to feel I'm pushing my luck. Fortunately, la05 has joined ylp in
  the dustbin, replaced by C basis maintenance code from glpk. At present
  (2005), this code is historical. Eventually it will disappear.

  The routines deal with a limited set of argument types: the Fortran types
  integer, double_precision, and character, and the special categories of
  variable and constraint names (these last two artifacts from the ylp LP
  code, the predecessor to dylp). A special type code is used to indicate
  the end of the list of arguments.

  Note that we simply ignore the length parameters which accompany the
  character strings (they follow ftnargEND as hidden parameters), under
  the assumption that an explict '\0' was placed at the end of each string
  when it was defined in the Fortran code.

  Older implementations of va_arg were written in such a way that you could
  actually write to va_arg(varargp,<type>). This changed with GCC 2.96,
  necessitating the fancy two-step of writing to varargp, then advancing it.
  We need to use the value to make gcc shut up, so turn it to advantage
  as a check this hack is still working. Speed isn't an issue here.
*/

void errmsg_ (integer *errid, char *ident, ... )

/*
  This routine provides the functionality of errmsg to Fortran code. Its
  primary purpose in life is to take the Fortran arguments and put them into
  a varargs parameter block that is acceptable to vfprintf. As mentioned above,
  the method used is a little dubious.

  The call, over in Fortran, looks something like
    errmsg(errno,rtnnme,ftnargtype1,arg1, ... ,ftnargtypen,argn,ftnargEND)

  By the time it gets here, of course, it's all pointers, and there are
  additional arguments tacked on at the end for character string lengths.
  We ignore the string length arguments, and expect that a null terminator
  has been added explicitly over in the Fortran code.

  Parameters: as explained above

  Returns: undefined
*/

{ va_list fargs,varargp ;
  int type ;

  double varargs[64] ;			/* double avoids alignment problems */
  int intarg ;
  double dblarg ;
  char *chararg ;

  char *rtnnme = "errmsg_" ;

/*
  Flush stdout and elogchn, so that any buffered normal output appears before
  the error message.
*/
  fflush(stdout) ;
  if (elogchn != NULL) fflush(elogchn) ;
/*
  Initialise various pointers and indices.
*/
  varargp = (va_list) &varargs[0] ;
  va_start(fargs,ident) ;
/*
  Put the ident string into the varargs block.
*/
  *((char **) varargp) = ident ;
  if (ident != va_arg(varargp,char *))
  { errmsg(1,rtnnme,__LINE__) ;
    return ; }
/*
  Now start up a loop to process the remainder of the arguments. For each
  <type,arg> pair, we pull off the type and use it to condition a switch
  with a case for each type of argument we're prepared to deal with.
*/
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
      { if (errecho == TRUE)
	  fprintf(stderr,
		  "\n%s: unrecognised Fortran argument type code %d.\n",
		  rtnnme,type) ;
	if (elogchn != NULL)
	  fprintf(elogchn,
		  "\n%s: unrecognised Fortran argument type code %d.\n",
		  rtnnme,type) ;
	return ; } }
  va_end(fargs) ;
/*
  That's it for the preprocessing. The rest of this code is just a copy of
  the standard errmsg used by the C part of the code.

  If there's no error message file, or we can't find the message, just print
  the error number. Otherwise call vfprintf to handle expanding the message
  text.
*/
  if (emsgchn == NULL || finderrmsg((int) *errid,&errtxt[ERRPFXLEN]) == NULL)
  { if (errecho == TRUE)
      fprintf(stderr,"\n%s: error %d.\n",ident,(int) *errid) ;
    if (elogchn != NULL)
      fprintf(elogchn,"\n%s: error %d.\n",ident,(int) *errid) ; }
  else
  { if (errecho == TRUE)
    { vfprintf(stderr,errtxt,((va_list) &varargs[0])) ;
      putc('\n',stderr) ; }
    if (elogchn != NULL) 
    { vfprintf(elogchn,errtxt,((va_list) &varargs[0])) ;
      putc('\n',elogchn) ; } }
/*
  Flush the logged error message, so that the user will definitely see it.
*/
  if (elogchn != NULL) fflush(elogchn) ;

  return ; }



#ifdef DYLP_NDEBUG

void warn_ (integer *errid, char *ident, ... )

{ return ; }

#else

void warn_ (integer *errid, char *ident, ... )

/*
  This routine provides the functionality of warn to Fortran code. Its
  primary purpose in life is to take the Fortran arguments and put them into
  a varargs parameter block that is acceptable to vfprintf. As mentioned
  above, the method used is a little dubious.

  The call, over in Fortran, looks something like
    warn(errno,rtnnme,ftnargtype1,arg1, ... ,ftnargtypen,argn,ftnargEND)

  By the time it gets here, of course, it's all pointers, and there are
  additional arguments tacked on at the end for character string lengths.
  We ignore the string length arguments, and expect that a null terminator
  has been added explicitly over in the Fortran code.

  Parameters: as explained above

  Returns: undefined
*/

{ va_list fargs,varargp ;
  int type ;

  double varargs[64] ;			/* double avoids alignment problems */
  int intarg ;
  double dblarg ;
  char *chararg ;

  char *rtnnme = "warn_" ;

/*
  Flush stdout and elogchn, so that any buffered normal output appears before
  the error message.
*/
  fflush(stdout) ;
  if (elogchn != NULL) fflush(elogchn) ;
/*
  Initialise various pointers and indices.
*/
  varargp = (va_list) &varargs[0] ;
  va_start(fargs,ident) ;
/*
  Put the ident string into the varargs block.
*/
  *((char **) varargp) = ident ;
  if (ident != va_arg(varargp,char *))
  { errmsg(1,rtnnme,__LINE__) ;
    return ; }
/*
  Now start up a loop to process the remainder of the arguments. For each
  <type,arg> pair, we pull off the type and use it to condition a switch
  with a case for each type of argument we're prepared to deal with.
*/
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
      { if (errecho == TRUE)
	  fprintf(stderr,
		  "\n%s: unrecognised Fortran argument type code %d.\n",
		  rtnnme,type) ;
	if (elogchn != NULL)
	  fprintf(elogchn,
		  "\n%s: unrecognised Fortran argument type code %d.\n",
		  rtnnme,type) ;
	return ; } }
  va_end(fargs) ;
/*
  That's it for the preprocessing. The rest of this code is just a copy of
  the standard warn used by the C part of the code.

  If there's no error message file, or we can't find the message, just print
  the error number. Otherwise call vfprintf to handle expanding the message
  text.
*/
  if (emsgchn == NULL || finderrmsg((int) *errid,&warntxt[WARNPFXLEN]) == NULL)
  { if (errecho == TRUE)
      fprintf(stderr,"\n%s: error %d.\n",ident,(int) *errid) ;
    if (elogchn != NULL)
      fprintf(elogchn,"\n%s: error %d.\n",ident,(int) *errid) ; }
  else
  { if (errecho == TRUE)
    { vfprintf(stderr,warntxt,((va_list) &varargs[0])) ;
    putc('\n',stderr) ; }
    if (elogchn != NULL) 
    { vfprintf(elogchn,warntxt,((va_list) &varargs[0])) ;
      putc('\n',elogchn) ; } }
/*
  Flush the logged error message, so that the user will definitely see it.
*/
  if (elogchn != NULL) fflush(elogchn) ;

  return ; }

#endif /* DYLP_NDEBUG */

#endif /* _DYLIB_FORTRAN */

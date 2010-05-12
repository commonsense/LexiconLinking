/*
  This file is a part of the OsiDylp LP distribution.

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
  This file contains a medium strength command interpreter. It's a generic
  piece of code, also used in the bonsaiG MIP code. For bare dylp, or dylp
  embedded in a COIN OSI layer, various bits should be excised. The symbol
  BONSAIG should be defined only if dylp is being built as part of bonsaiG.
*/



#include "dylib_strrtns.h"
#include "dy_cmdint.h"
#include "dylib_keytab.h"

static char sccsid[] UNUSED = "@(#)cmdint.c	3.12	11/11/04" ;
static char svnid[] UNUSED = "$Id: dy_cmdint.c 94 2006-06-29 23:06:51Z lou $" ;

/*
  MAXCMDFILES is pretty much arbitrary. The i/o package enforces the system
  limit. All we need here is something reasonable for sizing the cmdchns
  array.
*/

#define MAXCMDFILES 16

/*
  cmdchns keeps track of the i/o id and echo status for each command file
  that's open. level is the current nesting level.
*/

static int level ;
static bool prompt ;
static struct { ioid chn ;
		bool cecho ;
		bool gecho ;
		bool prompt ; } cmdchns[MAXCMDFILES-1] ;

/*
  User command table definitions. It is important that the commands be in 
  alphabetic order in the usercmds structure, as the keytab routines use
  binary search. The various _NDX values can be changed at will.

  For the dylp OSI layer, only lpcontrol, and lpprint are available.
*/

#ifdef BONSAIG
#  define HELP_NDX 1
#  define LPPRINT_NDX 2
#  define OBJECTIVE_NDX 3
#  define MIPPRINT_NDX 5
#  define PHIC_NDX 6
#  define OBJDELTA_NDX 7
#  define INCUMBENT_NDX 8
#  define RECOVER_NDX 9
#  define PRIORITY_NDX 10
#  define TOURCLASS_NDX 11
#  define LPCTL_NDX 12
#  define MIPCTL_NDX 13
#  define MIPLIM_NDX 14
#else
#  define HELP_NDX 1
#  define LPPRINT_NDX 2
#  define OBJECTIVE_NDX 255
#  define MIPPRINT_NDX 255
#  define PHIC_NDX 255
#  define OBJDELTA_NDX 255
#  define INCUMBENT_NDX 255
#  define RECOVER_NDX 255
#  define PRIORITY_NDX 255
#  define TOURCLASS_NDX 255
#  define LPCTL_NDX 12
#  define MIPCTL_NDX 255
#  define MIPLIM_NDX 255
#endif

#define UNIMP_NDX 255

static keytab_entry usercmds[] = { { "help", 1, HELP_NDX },
				   { "incumbent", 1, INCUMBENT_NDX },
				   { "lpcontrol", 3, LPCTL_NDX },
				   { "lpprint", 3, LPPRINT_NDX },
				   { "mipcontrol", 4, MIPCTL_NDX },
				   { "miplimit", 4, MIPLIM_NDX },
				   { "mipprint", 4, MIPPRINT_NDX },
				   { "objdelta", 4, OBJDELTA_NDX },
				   { "objective", 4, OBJECTIVE_NDX },
				   { "phic", 2, PHIC_NDX },
				   { "priority", 2, PRIORITY_NDX },
				   { "recovery", 1, RECOVER_NDX },
				   { "tourclass", 1, TOURCLASS_NDX }
				 } ;

#define NUMUSERCMDS (sizeof usercmds/sizeof (keytab_entry))




static cmd_retval docmd (lex_struct *txt)

/*
  This routine is just some prep and cleanup code, plus a big case statement.
  It figures out what command is requested, echoes it, dispatches to the
  appropriate routine, and on return cleans off any remaining junk on the
  command line.

  Parameter:
    txt:	the command name
  
  Returns: cmd_retval code supplied by the command execution routine, or
	   cmdHALTERROR if an error occurs here in docmd.
*/

{ lex_struct *lex ;
  const char *keywd ;
  int cmd ;
  cmd_retval retval ;
  const char *rtnnme = "docmd" ;

/* dy_setup.c */

  extern cmd_retval dy_printopt(const char *keywd),
		    dy_ctlopt(const char *keywd) ;

#ifdef BONSAIG

/* mip_setup.c */

  extern cmd_retval mip_objopt(const char *keywd),
		    mip_objdeltaopt(const char *keywd),
		    mip_incumbentopt(const char *keywd),
		    mip_phicopt(const char *keywd),
		    mip_printopt(const char *keywd),
		    mip_ctlopt(const char *keywd),
		    mip_recoveryopt(const char *keywd) ;

/* premature.c */

  extern cmd_retval mip_limit (const char *keywd) ;

/* priority_utils.c */

  extern cmd_retval pri_priorityopt(const char *keywd) ;

/* tourclass_utils.c */

  extern cmd_retval tc_tourclassopt(const char *keywd) ;

#endif


  if (txt == NULL)
  { errmsg(2,rtnnme,"txt") ;
    return (cmdHALTERROR) ; }
  if (txt->class != LCID)
  { errmsg(5,rtnnme,(int) txt->class) ;
    return (cmdHALTERROR) ; }
/*
  Look up the command in the command table. Return code of -1 means we can't
  find the string, -2 that it's ambiguous.
*/
  dyio_outfmt(dy_logchn,dy_cmdecho,txt->string) ;
  dyio_flushio(dy_logchn,dy_cmdecho) ;
  cmd = ambig(txt->string,usercmds,NUMUSERCMDS) ;
  if (cmd < 0) 
  { if (cmd < -1)
      errmsg(233,rtnnme,txt->string) ;
    else
      errmsg(234,rtnnme,txt->string) ;
    return (cmdHALTERROR) ; }
/*
  Call the appropriate execution routine. Go to free input mode while parsing
  the command, then revert to line-oriented mode to remove any trailing junk.
*/
  keywd = STRALLOC(txt->string) ;
  (void) dyio_setmode(dy_cmdchn,'f') ;
  retval = cmdHALTERROR ;
  (void) dyio_setmode(dy_cmdchn,'f') ;
  switch (cmd)
  { case LPPRINT_NDX:
    { retval = dy_printopt(keywd) ;
      break ; }
    case LPCTL_NDX:
    { retval = dy_ctlopt(keywd) ;
      break ; }
#   ifdef BONSAIG
    case MIPPRINT_NDX:
    { retval = mip_printopt(keywd) ;
      break ; }
    case MIPCTL_NDX:
    { retval = mip_ctlopt(keywd) ;
      break ; }
    case MIPLIM_NDX:
    { retval = mip_limit(keywd) ;
      break ; }
    case OBJECTIVE_NDX:
    { retval = mip_objopt(keywd) ;
      break ; }
    case OBJDELTA_NDX:
    { retval = mip_objdeltaopt(keywd) ;
      break ; }
    case PHIC_NDX:
    { retval = mip_phicopt(keywd) ;
      break ; }
    case INCUMBENT_NDX:
    { retval = mip_incumbentopt(keywd) ;
      break ; }
    case RECOVER_NDX:
    { retval = mip_recoveryopt(keywd) ;
      break ; }
    case PRIORITY_NDX:
    { retval = pri_priorityopt(keywd) ;
      break ; }
    case TOURCLASS_NDX:
    { retval = tc_tourclassopt(keywd) ;
      break ; }
#   endif
    case HELP_NDX:
    { dyio_outfmt(dy_logchn,dy_gtxecho,
		  "\n\tHeh heh heh.\t\tSurely you jest?\n") ;
      retval = cmdOK ;
      break ; }
    case UNIMP_NDX:
    { dyio_outfmt(dy_logchn,dy_gtxecho,
	      "\nThe command \"%s\" is unimplemented in this configuration.\n",
	      keywd) ;
      retval = cmdOK ;
      break ; } }
  STRFREE(keywd) ;
  (void) dyio_setmode(dy_cmdchn,'l') ;
/*
  Check what happened, and clean up if necessary. Cleanup consists of
  scanning off whatever is left on the command line and echoing it.
*/
  if (retval != cmdHALTERROR)
  { lex =  dyio_scanstr(dy_cmdchn,LCQS,0,'\0','\n') ;
    if (!(lex->class == LCNIL || lex->class == LCEOF || lex->class == LCERR))
      dyio_outfmt(dy_logchn,dy_cmdecho," %s",lex->string) ;
    if (lex->class == LCERR) retval = cmdHALTERROR ; }
  
  return (retval) ; }



static cmd_retval indcmd (bool silent)

/*
  This routine takes care of opening an indirect command file and doing the
  necessary bookkeeping to change command channels. The form of the command
  to open an indirect command file is '@ "filename"'.

  Parameters:
    silent: TRUE if echoing of commands & responses to ttyout should be
	    suppressed, FALSE otherwise.

  Returns: TRUE if all goes well, FALSE otherwise
*/

{ lex_struct *file ;
  const char *rtnnme = "indcmd" ;

/*
  Get the file name.
*/
  file = dyio_scanstr(dy_cmdchn,LCQS,0,'"','"') ;
  if (file->class != LCQS)
  { errmsg(236,rtnnme,"file name","parameter","@") ;
    return (cmdHALTERROR) ; }
/*
  Stash the old command channel, then try to open the new file.
*/
  cmdchns[level].chn = dy_cmdchn ;
  cmdchns[level].cecho = dy_cmdecho ;
  cmdchns[level].gecho = dy_gtxecho ;
  cmdchns[level].prompt = prompt ;
  dy_cmdchn = dyio_openfile(file->string,"r") ;
  if (dy_cmdchn < 0)
  { dy_cmdchn = cmdchns[level].chn ;
    return (cmdHALTERROR) ; }
  (void) dyio_setmode(dy_cmdchn,'l') ;
/*
  Acknowledge that the file is successfully opened.
*/
  dyio_outfmt(dy_logchn,dy_cmdecho," \"%s\"\n",file->string) ;
  dyio_outfmt(dy_logchn,dy_gtxecho,
	      "\tcommand source file now %s\n",file->string) ;
/*
  Set level and echo control variables. dy_gtxecho stays the same, so no action
  for it. We have to check dy_cmdecho and prompt, in case the input has switched
  to ttyin.
*/
  level++ ;
  if (dyio_ttyq(dy_cmdchn) == TRUE || silent == FALSE)
    prompt = TRUE ;
  else
    prompt = FALSE ;
  if (dyio_ttyq(dy_cmdchn) == FALSE && silent == FALSE)
    dy_cmdecho = TRUE ;
  else
    dy_cmdecho = FALSE ;

  return (cmdOK) ; }



static cmd_retval dobuiltin (lex_struct *txt, bool silent)

/*
  This routine handles the built-in functions of the command interpreter.
  Not much, really -- disposal of empty lines and comments, opening
  indirect command files.

  Parameters:
    txt:    the first lexeme of the command line (a delimiter).
    silent: TRUE if echoing of commands & responses to ttyout should be
	    suppressed, FALSE otherwise.
  
  Returns: the appropriate cmd_retval code.
*/

{ lex_struct *lex ;
  cmd_retval retval ;
  const char *rtnnme = "dobuiltin" ;

  if (txt == NULL)
  { errmsg(2,rtnnme,"txt") ;
    return (cmdHALTERROR) ; }
  if (txt->class != LCDEL)
  { errmsg(5,rtnnme,(int) txt->class) ;
    return (cmdHALTERROR) ; }
/*
  Deal with empty lines, comments, and garbage.
*/
  if (*txt->string == '\n') return cmdOK ;

  if (*txt->string == '!')
  { dyio_outchr(dy_logchn,dy_cmdecho,'!') ;
    lex = dyio_scanstr(dy_cmdchn,LCQS,0,'\0','\n') ;
    if (lex->class != LCNIL)
    { dyio_outfmt(dy_logchn,dy_gtxecho," %s",lex->string) ; }
    return (cmdOK) ; }

  if (*txt->string != '@')
  { errmsg(230,rtnnme,txt->string) ;
    return (cmdHALTERROR) ; }
/*
  The various 'built-in' commands. For now, just opening an indirect command
  file.
*/
  retval = indcmd(silent) ;

  return (retval) ; }



cmd_retval process_cmds (bool silent)

/*
  This routine is a driver routine for processing a command file. dy_cmdchn and
  dy_logchn must be valid before process_cmds is called.

  Parameters:
    silent: TRUE if echoing of commands & responses to ttyout should be
	    suppressed, FALSE otherwise.

  Returns: cmdOK 	  if eof is reached and the nesting level is 0
	   cmdHALTNOERROR if the return code from a command execution
		 	  specified terminate w/o error
	   cmdHALTERROR	  otherwise
*/

{ lex_struct *txt ;
  cmd_retval retval ;
  const char *rtnnme = "process_cmds" ; 

/*
  Initialise to level 0 and begin the command interpretation loop. We need to
  generate a prompt if the command channel is the terminal, or if the user
  wants commands to be echoed. We want to echo commands if the command channel
  is not a terminal and the user wants commands to be echoed. We echo generated
  text if the user wants it.
*/
  level = 0 ;
  txt = NULL ;
  retval = cmdOK ;
  if (dyio_ttyq(dy_cmdchn) == TRUE || silent == FALSE)
    prompt = TRUE ;
  else
    prompt = FALSE ;
  if (dyio_ttyq(dy_cmdchn) == FALSE && silent == FALSE)
    dy_cmdecho = TRUE ;
  else
    dy_cmdecho = FALSE ;
  dy_gtxecho = !silent ;
/*
  Open the command interpretation loop. First action is to prompt for a
  command. Remember that a return value of cmdOK from a command execution
  is interpreted to mean carry on with command interpretation.
*/
  while (retval == cmdOK)
  { dyio_outfmt(dy_logchn,prompt,"\n(%d)%% ",level) ;
    retval = cmdHALTERROR ;
/*
  Read in the first lexeme of the response and act accordingly. With any luck,
  we'll have an id, which we pass along to docmd to process. A delimiter could
  be an empty line, comment, or builtin, and they're handled by dobuiltin. EOF
  means we've finished a command file. Attempt a close and continue at the next
  level up. If we're at level 0, return. All other possibilities are errors of
  one sort or another.
*/
    txt = dyio_scanlex(dy_cmdchn) ;
    switch (txt->class)
    { case LCID:
      { retval = docmd(txt) ;
	break ; }
      case LCDEL:
      { retval = dobuiltin(txt,silent) ;
	break ; }
      case LCEOF:
      { if (level == 0)
	{ retval = cmdHALTNOERROR ; }
	else
	{ if (dyio_closefile(dy_cmdchn) == FALSE)
	    warn(232,rtnnme,dyio_idtopath(dy_cmdchn)) ;
	  dy_cmdchn = cmdchns[--level].chn ;
	  dyio_outfmt(dy_logchn,dy_gtxecho,
		      "\n\treturning to command source file %s\n",
		      dyio_idtopath(dy_cmdchn)) ;
	  dy_cmdecho = cmdchns[level].cecho ;
	  dy_gtxecho = cmdchns[level].gecho ;
	  retval = cmdOK ; }
	break ; }
      case LCERR:
      { break ; }
      default: /* LCNIL, LCNUM, LCFS, LCQS */
      { errmsg(230,rtnnme,(txt->string == NULL)?"<<nil>>":txt->string) ;
	break ; } }
    dyio_flushio(dy_logchn,dy_cmdecho) ; }
/*
  Command interpretation has been stopped. Fix up the return code and return.
*/
  if (retval == cmdHALTERROR) errmsg(235,rtnnme) ;
  if (retval == cmdHALTNOERROR && txt->class == LCEOF) retval = cmdOK ;

  return (retval) ; }

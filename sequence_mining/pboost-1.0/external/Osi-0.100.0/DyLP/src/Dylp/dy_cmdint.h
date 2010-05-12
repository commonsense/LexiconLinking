/*
  This file is a part of the Dylp LP distribution.

        Copyright (C) 2005 -- 2007 Lou Hafer

        School of Computing Science
        Simon Fraser University
        Burnaby, B.C., V5A 1S6, Canada
        lou@cs.sfu.ca

  This code is licensed under the terms of the Common Public License (CPL).
*/

#ifndef _DY_CMDINT_H
#define _DY_CMDINT_H

/*
  @(#)dy_cmdint.h	3.3	06/22/04
  svn/cvs: $Id: dy_cmdint.h 148 2007-06-09 03:15:30Z lou $

  Declarations specific to dylp's command interpreter.
*/

#include "dylib_std.h"
#include "dylib_io.h"
#include "dylib_errs.h"

/*
  Globals for log stream and echo control, command input stream and echo
  control. These must be declared in a main program somewhere.

  dy_cmdchn		i/o id for command input
  dy_cmdecho		controls echoing of command input to stdout

  dy_logchn		i/o id for log file
  dy_gtxecho		controls echoing of generated text to stdout

  dylp.h also contains extern declarations for dy_logchn and dy_gtxecho. Turns out
  that the files related to the command interpreter don't need the main dylp
  structures, so it's useful to duplicate the extern decl's in both .h files.
*/

extern ioid dy_logchn,dy_cmdchn ;
extern bool dy_gtxecho,dy_cmdecho ;


/*
  cmdint.c
*/

/*
  Return codes for command execution routines called from the command
  interpreter:

    cmdOK	execution of the command was adequately successful, further
		command interpretation should continue.
    cmdHALTNOERROR execution of the command was adequately successful, but break
		out of the command interpretation loop.
    cmdHALTERROR an error occurred during execution of the command, break
		out of the command interpretation loop.

  As return codes for process_cmds, the interpretation is slightly different:
    cmdOK	command interpretation was ended by an eof on the top level
		command channel (this is the normal case when command execution
		completes without error).
    cmdHALTNOERROR some command returned a cmdHALTNOERROR return code.
    cmdHALTERROR either a command returned a cmdHALTERROR return code, or a
		fatal error occurred in process_cmds.
*/

typedef enum { cmdOK, cmdHALTERROR, cmdHALTNOERROR } cmd_retval ;

cmd_retval process_cmds(bool silent) ;

#endif	/* _DY_CMDINT_H */

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
  This file contains generic routines to process simple options. In the main,
  these are called by the option processing routines that deal with the
  options file.
*/


#include "dylib_bnfrdr.h"
#include "dylib_strrtns.h"
#include "dylib_keytab.h"
#include "dy_cmdint.h"
#include "dylp.h"

#include <limits.h>
#include <float.h>

static char sccsid[] UNUSED = "@(#)options.c	3.5	09/25/04" ;
static char svnid[] UNUSED = "$Id: dy_options.c 94 2006-06-29 23:06:51Z lou $" ;

/*
  Routines to acquire dylp's defaults and limits for options and tolerances.
*/

extern void dy_exposeOptDefaults(lpopts_struct **opts_lb,
				 lpopts_struct **opts_dflt,
				 lpopts_struct **opts_ub),
	    dy_exposeTolDefaults(lptols_struct **tols_dflt) ;

/*
  We need only the options and tolerances structure from the main milp code,
  so a solitary declaration seems preferable to dragging in all of milp.h
  At best, though, this is a hack, and I need to deal with it properly, to
  retain dylp as a self-contained package.
*/

extern lpopts_struct *main_lpopts ;
extern lptols_struct *main_lptols ;



static bool string_opt (char **str)

/*
  Generic (and fairly trivial) routine to parse a string.

  Parameters:
    str: address where the string is to be placed.

  Returns: TRUE, barring some sort of parsing failure.
*/

{ parse_any result ;
  const char *rtnnme = "string_opt" ;

/*
  BNF to parse a string, returning (char *).
*/
  static tdef(zid,bnfttID,NULL,NULL) ;
  static tref(zgetstring_zid,zid,bnfstore|bnfatsgn,0) ;

  static comphd(zgetstring_alt) = { compcnt(1), mkcref(zgetstring_zid) } ;
  static gdef(zgetstring_int,sizeof(char *),NULL,zgetstring_alt) ;
  static gref(zgetstring,zgetstring_int,NULL,NULL,NULLP) ;

/*
  Initialise the reader, parse the string, and check for errors.
*/
  rdrinit() ;
  if (parse(dy_cmdchn,&zgetstring,&result) == FALSE)
  { rdrclear() ;
    errmsg(240,rtnnme,"zgetstring") ;
    return (FALSE) ; }
  dyio_outfmt(dy_logchn,dy_cmdecho," %s",*(char **) result.g) ;
  dyio_flushio(dy_logchn,dy_cmdecho) ;

  *str = *(char **) result.g ;

/*
  Clean up and return.
*/
  FREE(result.g) ;
  rdrclear() ;

  return (TRUE) ; }



static bool integer_opt (int *iloc)

/*
  Generic (and fairly trivial) routine to parse an integer.

  Parameters:
    iloc: address where the integer is to be placed.

  Returns: TRUE, barring some sort of parsing failure.
*/

{ parse_any result ;
  const char *rtnnme = "integer_opt" ;

/*
  BNF to parse an integer, returning (int *).
*/
  static tdef(zdnum,bnfttN,10,NULL) ;
  static tref(zgetnum_zdnum,zdnum,bnfstore,0) ;
  static comphd(zgetnum_alt) = { compcnt(1), mkcref(zgetnum_zdnum) } ;
  static gdef(zgetnum_int,sizeof(int),NULL,zgetnum_alt) ;
  static gref(zgetnum,zgetnum_int,NULL,NULL,NULLP) ;

/*
  Initialise the reader, parse the integer, and check for errors.
*/
  rdrinit() ;
  if (parse(dy_cmdchn,&zgetnum,&result) == FALSE)
  { rdrclear() ;
    errmsg(240,rtnnme,"zgetnum") ;
    return (FALSE) ; }
  dyio_outfmt(dy_logchn,dy_cmdecho," %d",*(int *) result.g) ;
  dyio_flushio(dy_logchn,dy_cmdecho) ;

  *iloc = *(int *) result.g ;

/*
  Clean up and return.
*/
  FREE(result.g) ;
  rdrclear() ;

  return (TRUE) ; }



static bool double_opt (double *rloc)

/*
  Generic (and fairly trivial) routine to parse a real as a double. A double
  will provide about 15 -- 17 significant decimal digits.

  Parameters:
    rloc: address where the real is to be placed.

  Returns: TRUE, barring some sort of parsing failure.
*/

{ parse_any result ;
  const char *rtnnme = "double_opt" ;

/*
  BNF to parse a real, returning (double *).
*/
  static tdef(zdnum,bnfttN,10,NULL) ;
  static tref(zgetnum_zdnum,zdnum,bnfdbl|bnfstore,0) ;
  static comphd(zgetnum_alt) = { compcnt(1), mkcref(zgetnum_zdnum) } ;
  static gdef(zgetnum_int,sizeof(double),NULL,zgetnum_alt) ;
  static gref(zgetnum,zgetnum_int,NULL,NULL,NULLP) ;

/*
  Initialise the reader, parse the real, and check for errors.
*/
  rdrinit() ;
  if (parse(dy_cmdchn,&zgetnum,&result) == FALSE)
  { rdrclear() ;
    errmsg(240,rtnnme,"zgetnum") ;
    return (FALSE) ; }
  dyio_outfmt(dy_logchn,dy_cmdecho," %g",*(double *) result.g) ;
  dyio_flushio(dy_logchn,dy_cmdecho) ;

  *rloc = *(double *) result.g ;

/*
  Clean up and return.
*/
  FREE(result.g) ;
  rdrclear() ;

  return (TRUE) ; }



static bool real_opt (float *rloc)

/*
  Generic (and fairly trivial) routine to parse a real as a float. A float
  will provide about 6 -- 9 significant decimal digits.

  Parameters:
    rloc: address where the real is to be placed.

  Returns: TRUE, barring some sort of parsing failure.
*/

{ parse_any result ;
  const char *rtnnme = "real_opt" ;

/*
  BNF to parse a real, returning (float *).
*/
  static tdef(zdnum,bnfttN,10,NULL) ;
  static tref(zgetnum_zdnum,zdnum,bnfflt|bnfstore,0) ;
  static comphd(zgetnum_alt) = { compcnt(1), mkcref(zgetnum_zdnum) } ;
  static gdef(zgetnum_int,sizeof(float),NULL,zgetnum_alt) ;
  static gref(zgetnum,zgetnum_int,NULL,NULL,NULLP) ;

/*
  Initialise the reader, parse the real, and check for errors.
*/
  rdrinit() ;
  if (parse(dy_cmdchn,&zgetnum,&result) == FALSE)
  { rdrclear() ;
    errmsg(240,rtnnme,"zgetnum") ;
    return (FALSE) ; }
  dyio_outfmt(dy_logchn,dy_cmdecho," %g",*(float *) result.g) ;
  dyio_flushio(dy_logchn,dy_cmdecho) ;

  *rloc = *(float *) result.g ;

/*
  Clean up and return.
*/
  FREE(result.g) ;
  rdrclear() ;

  return (TRUE) ; }



static bool bool_opt (bool *bloc)

/*
  Generic (and fairly trivial) routine to parse a boolean. But ... we have to
  take some care here because of the way bool is handled. If you look at the
  typedef for bool in loustd.h, you'll see that it can change in size --- this
  is necessary for C++ compatibility in the COIN OSI layer implementation. But
  a bnfIdef_struct (an immediate) holds its value as an int, and when
  doimmediate tries to load a field, it casts to an int. To make a long
  story short, the val field in a boolopt_struct must be an int.

  Parameters:
    bloc: address where the boolean is to be placed.

  Returns: TRUE, barring some sort of parsing failure.
*/

{ struct boolopt_struct { char *str ;
			  int val ; } *boolopt ;
  parse_any result ;
  const char *rtnnme = "bool_opt" ;

/*
  BNF to parse a boolean, returning the string.
*/
  static idef(ziTRUE,TRUE) ;
  static idef(ziFALSE,FALSE) ;
  static tdef(zTRUE,bnfttID,NULL,"TRUE") ;
  static tdef(zFALSE,bnfttID,NULL,"FALSE") ;

  static iref(zparsebool_ziTRUE,ziTRUE,mkoff(struct boolopt_struct,val)) ;
  static iref(zparsebool_ziFALSE,ziFALSE,mkoff(struct boolopt_struct,val)) ;
  static tref(zparsebool_zTRUE,zTRUE,bnfstore|bnfatsgn|bnfmin,
	      mkoff(struct boolopt_struct,str)) ;
  static tref(zparsebool_zFALSE,zFALSE,bnfstore|bnfatsgn|bnfmin,
	      mkoff(struct boolopt_struct,str)) ;
  static comphd(zparsebool_alt1) = { compcnt(2),
      mkcref(zparsebool_zFALSE), mkcref(zparsebool_ziFALSE) } ;
  static comphd(zparsebool_alt2) = { compcnt(2),
      mkcref(zparsebool_zTRUE), mkcref(zparsebool_ziTRUE) } ;
  static althd(zparsebool_alts) = { altcnt(2),
      mkaref(zparsebool_alt1), mkaref(zparsebool_alt2) } ;
  static npdef(zparsebool,zparsebool_alts) ;
  static npref(zparsebool_ref,zparsebool,NULL,NULLP) ;

  static comphd(zgetbool_alt) = { compcnt(1), mkcref(zparsebool_ref) } ;
  static gdef(zgetbool_int,sizeof(struct boolopt_struct),NULL,zgetbool_alt) ;
  static gref(zgetbool,zgetbool_int,NULL,NULL,NULLP) ;

/*
  Initialise the reader, parse the string, and check for errors.
*/
  rdrinit() ;
  if (parse(dy_cmdchn,&zgetbool,&result) == FALSE)
  { rdrclear() ;
    errmsg(240,rtnnme,"zgetbool") ;
    return (FALSE) ; }
  boolopt = (struct boolopt_struct *) result.g ;
  rdrclear() ;
  
  dyio_outfmt(dy_logchn,dy_cmdecho," %s",boolopt->str) ;
  dyio_flushio(dy_logchn,dy_cmdecho) ;

  *bloc = boolopt->val ;

/*
  Clean up and return.
*/
  STRFREE(boolopt->str) ;
  FREE(boolopt) ;

  return (TRUE) ; }



/*
  Option processing routines. These are exported only to cmdint.c:docmds
  where they are called. In turn, they will often make use of one of the set
  of generic option processing routines above.
*/

cmd_retval dy_printopt (const char *keywd)

/*
  This routine handles the manifold parameters that control just how much dylp
  tells you while it's working. This routine doesn't return an error, even
  when it can't parse the command, on the general theory that print levels
  aren't really critical. It will complain, however.

  The bnf for the print option command is:
    <printopt> ::= lpprint <what> <level>
    <what> ::= basis | conmgmt | crash | degen | dual | force | major |
	       phase1 | phase2 | pivoting | pivreject | pricing |
	       scaling | setup | varmgmt
    <level> ::= <integer>

  Parameters:
    keywd:	The command keyword

  Returns: cmdOK
*/

{ char *what,cmdstr[50] ;
  int code,dflt,lb,ub,*opt ;
  lpopts_struct *opts_lb,*opts_dflt,*opts_ub ;
  const char *rtnnme = "dy_printopt" ;

/*
  A lookup table with the various <what> keywords recognised by the print
  command.
*/

  enum prntcodes { poINV = 0, poBASIS, poCONMGMT, poCRASH,
		   poDEGEN, poDUAL, poFORCE, poMAJOR,
		   poPHASE1, poPHASE2, poPIVOTING, poPIVREJ,
		   poPRICING, poSCALING, poSETUP, poVARMGMT } prntcode ;

  static keytab_entry prntkwds[] = { { "basis", 1, (int) poBASIS },
				     { "conmgmt", 2, (int) poCONMGMT },
				     { "crash", 2, (int) poCRASH },
				     { "degen", 2, (int) poDEGEN },
				     { "dual", 2, (int) poDUAL },
				     { "force", 2, (int) poFORCE },
				     { "major", 1, (int) poMAJOR },
				     { "phase1", 6, (int) poPHASE1 },
				     { "phase2", 6, (int) poPHASE2 },
				     { "pivoting", 4, (int) poPIVOTING },
				     { "pivreject", 4, (int) poPIVREJ },
				     { "pricing", 2, (int) poPRICING },
				     { "scaling", 2, (int) poSCALING },
				     { "setup", 2, (int) poSETUP },
				     { "varmgmt", 1, (int) poVARMGMT }
				   } ;

  static int numprntcodes = (sizeof prntkwds/sizeof (keytab_entry)) ;

/*
  Acquire the option and tolerance limits and defaults.
*/
  dy_exposeOptDefaults(&opts_lb,&opts_dflt,&opts_ub) ;
/*
  Now to work. Parse off the <what> keyword and see if we can look it up.
*/
  prntcode = poINV ;
  if (string_opt(&what) == TRUE)
  { code = ambig(what,prntkwds,numprntcodes) ;
    if (code < 0) 
    { if (code < -1)
	errmsg(233,rtnnme,what) ;
      else
	errmsg(234,rtnnme,what) ; }
    else
      prntcode = (enum prntcodes) code ; }
/*
  Set the various variables for each command.
*/
  dyio_outfxd(cmdstr,-((int) (sizeof(cmdstr)-1)),'l',"%s %s",keywd,what) ;
  switch (prntcode)
  { case poCONMGMT:
    { opt = &main_lpopts->print.conmgmt ;
      dflt = opts_dflt->print.conmgmt ;
      lb = opts_lb->print.conmgmt ;
      ub = opts_ub->print.conmgmt ;
      break ; }
    case poCRASH:
    { opt = &main_lpopts->print.crash ;
      dflt = opts_dflt->print.crash ;
      lb = opts_lb->print.crash ;
      ub = opts_ub->print.crash ;
      break ; }
    case poDEGEN:
    { opt = &main_lpopts->print.degen ;
      dflt = opts_dflt->print.degen ;
      lb = opts_lb->print.degen ;
      ub = opts_ub->print.degen ;
      break ; }
    case poBASIS:
    { opt = &main_lpopts->print.basis ;
      dflt = opts_dflt->print.basis ;
      lb = opts_lb->print.basis ;
      ub = opts_ub->print.basis ;
      break ; }
    case poMAJOR:
    { opt = &main_lpopts->print.major ;
      dflt = opts_dflt->print.major ;
      lb = opts_lb->print.major ;
      ub = opts_ub->print.major ;
      break ; }
    case poPHASE1:
    { opt = &main_lpopts->print.phase1 ;
      dflt = opts_dflt->print.phase1 ;
      lb = opts_lb->print.phase1 ;
      ub = opts_ub->print.phase1 ;
      break ; }
    case poPHASE2:
    { opt = &main_lpopts->print.phase2 ;
      dflt = opts_dflt->print.phase2 ;
      lb = opts_lb->print.phase2 ;
      ub = opts_ub->print.phase2 ;
      break ; }
    case poDUAL:
    { opt = &main_lpopts->print.dual ;
      dflt = opts_dflt->print.dual ;
      lb = opts_lb->print.dual ;
      ub = opts_ub->print.dual ;
      break ; }
    case poFORCE:
    { opt = &main_lpopts->print.force ;
      dflt = opts_dflt->print.force ;
      lb = opts_lb->print.force ;
      ub = opts_ub->print.force ;
      break ; }
    case poPIVOTING:
    { opt = &main_lpopts->print.pivoting ;
      dflt = opts_dflt->print.pivoting ;
      lb = opts_lb->print.pivoting ;
      ub = opts_ub->print.pivoting ;
      break ; }
    case poPIVREJ:
    { opt = &main_lpopts->print.pivreject ;
      dflt = opts_dflt->print.pivreject ;
      lb = opts_lb->print.pivreject ;
      ub = opts_ub->print.pivreject ;
      break ; }
    case poPRICING:
    { opt = &main_lpopts->print.pricing ;
      dflt = opts_dflt->print.pricing ;
      lb = opts_lb->print.pricing ;
      ub = opts_ub->print.pricing ;
      break ; }
    case poSCALING:
    { opt = &main_lpopts->print.scaling ;
      dflt = opts_dflt->print.scaling ;
      lb = opts_lb->print.scaling ;
      ub = opts_ub->print.scaling ;
      break ; }
    case poSETUP:
    { opt = &main_lpopts->print.setup ;
      dflt = opts_dflt->print.setup ;
      lb = opts_lb->print.setup ;
      ub = opts_ub->print.setup ;
      break ; }
    case poVARMGMT:
    { opt = &main_lpopts->print.varmgmt ;
      dflt = opts_dflt->print.varmgmt ;
      lb = opts_lb->print.varmgmt ;
      ub = opts_ub->print.varmgmt ;
      break ; }
    default:
    { errmsg(236,rtnnme,"<what>","keyword",keywd) ;
      return (cmdOK) ; } }
/*
  Last but not least, the actual work. A negative value is taken as a request
  from the user to be told the default value.
*/
  if (integer_opt(opt) == TRUE)
  { if (*opt >= 0)
    { if (*opt > ub)
      { warn(241,rtnnme,lb,cmdstr,ub,*opt,ub) ;
	*opt  = ub ; } }
    else
    { warn(243,rtnnme,cmdstr,dflt) ; } }
  else
  { errmsg(236,rtnnme,"<level>","parameter",keywd) ; }

  STRFREE(what) ;

  return (cmdOK) ; }



static bool lpctl_active (void)

/*
  This routine processes the 'active' subcommand, which sets values that
  determine the initial size of dylp's copy of the constraint system. The
  values are expressed in terms of fractions of the number of inequalities
  and number of variables in the original system.

  The bnf for the active subcommand is:
    <activeopt> ::= lpcontrol active  <fracspec-LIST(,)> ;
    <fracspec> ::= variables <float> | constraints <float>
  By the time this routine is called, `lpcontrol active' is already parsed.

  Parameters: none

  Returns: TRUE if the remainder of the command parses without error, FALSE
	   otherwise.
*/

{ struct actfrac_struct { int var_seen ;
			  float varfrac ;
			  int con_seen ;
			  float confrac ; } *actfrac ;
  lpopts_struct *opts_lb,*opts_dflt,*opts_ub ;
  parse_any result ;

  const char *rtnnme = "lpctl_active" ;

/*
  BNF for the active command.
*/
  static tdef(zcomma,bnfttD,NULL,",") ;
  static tref(zcomma_ref,zcomma,NULL,NULL) ;
  static tdef(zdnum,bnfttN,10,NULL) ;
  static tref(zactfrac_varfrac,zdnum,bnfstore|bnfflt,
	      mkoff(struct actfrac_struct,varfrac)) ;
  static tref(zactfrac_confrac,zdnum,bnfstore|bnfflt,
	      mkoff(struct actfrac_struct,confrac)) ;
  static tdef(zvar,bnfttID,NULL,"variables") ;
  static tref(zvar_ref,zvar,bnfmin,NULL) ;
  static tdef(zcon,bnfttID,NULL,"constraints") ;
  static tref(zcon_ref,zcon,bnfmin,NULL) ;
  static idef(ziTRUE,TRUE) ;
  static iref(zactfrac_varseen,ziTRUE,mkoff(struct actfrac_struct,var_seen)) ;
  static iref(zactfrac_conseen,ziTRUE,mkoff(struct actfrac_struct,con_seen)) ;

  static comphd(zactfrac_alt1) = { compcnt(3),
    mkcref(zvar_ref),mkcref(zactfrac_varseen),mkcref(zactfrac_varfrac) } ;
  static comphd(zactfrac_alt2) = { compcnt(3),
    mkcref(zcon_ref),mkcref(zactfrac_conseen),mkcref(zactfrac_confrac) } ;
  static althd(zactfrac_alts) = { altcnt(2),
    mkaref(zactfrac_alt1),mkaref(zactfrac_alt2) } ;
  static npdef(zactfrac,zactfrac_alts) ;
  static npref(zactfrac_list,zactfrac,bnflst,zcomma_ref) ;

  static comphd(zactfracs_alt) = { compcnt(1), mkcref(zactfrac_list) } ;
  static gdef(zactfracs_int,sizeof(struct actfrac_struct),NULL,
	      zactfracs_alt) ;
  static gref(zactfracs,zactfracs_int,NULL,NULL,NULLP) ;

  dy_exposeOptDefaults(&opts_lb,&opts_dflt,&opts_ub) ;
/*
  Now to work.  Initialise the reader, attempt to parse the command, and
  check for errors.
*/
  rdrinit() ;
  if (parse(dy_cmdchn,&zactfracs,&result) == FALSE)
  { rdrclear() ;
    errmsg(240,rtnnme,"zactfracs") ;
    return (FALSE) ; }
  actfrac = (struct actfrac_struct *) result.g ;
  rdrclear() ;
/*
  Process the results.
*/
  if (actfrac->var_seen == TRUE)
  { dyio_outfmt(dy_logchn,dy_gtxecho," variables %.2f",actfrac->varfrac) ;
    if (actfrac->con_seen == TRUE)
      dyio_outchr(dy_logchn,dy_gtxecho,',') ;
    if (actfrac->varfrac >= 0)
    { if (actfrac->varfrac > opts_ub->active.vars)
      { warn(244,rtnnme,opts_lb->active.vars,"variables",
	     opts_ub->active.vars,actfrac->varfrac,opts_ub->active.vars) ;
	main_lpopts->active.vars = opts_ub->active.vars ; }
      else
	main_lpopts->active.vars = actfrac->varfrac ; }
    else
    { warn(245,rtnnme,"variables",opts_dflt->active.vars) ; } }

  if (actfrac->con_seen == TRUE)
  { dyio_outfmt(dy_logchn,dy_gtxecho," constraints %.2f",actfrac->confrac) ;
    if (actfrac->confrac >= 0)
    { if (actfrac->confrac > opts_ub->active.cons)
      { warn(244,rtnnme,opts_lb->active.cons,"constraints",
	     opts_ub->active.cons,actfrac->confrac,opts_ub->active.cons) ;
	main_lpopts->active.cons = opts_ub->active.cons ; }
      else
	main_lpopts->active.cons = actfrac->confrac ; }
    else
    { warn(245,rtnnme,"constraints",opts_dflt->active.cons) ; } }

  FREE(actfrac) ;

  return (TRUE) ; }



static bool lpctl_finpurge (void)

/*
  This routine processes the "final" subcommand, which currently controls
  only whether dylp does a final round of variable and/or constraint
  deactivation after finding an optimal solution. Note that this option works
  even when the problem has been solved with the fullsys option. This is
  intended for use in branch&cut, where we usually want to solve the initial
  lp with fullsys, then continue with dynamic lp.

  The bnf for the final purge subcommand is
    <final> ::= lpcontrol final purge <purgespec-LIST(,)> ;
    <purgespec> ::= <what> <sense>
    <what> ::= variables | constraints
    <sense> ::= true | false
  By the time this routine is called, `lpcontrol final' is already parsed.

  Parameters: none

  Returns: TRUE if the remainder of the command parses without error, FALSE
	   otherwise.
*/

{ struct finpurge_struct { int cons ;
			   int vars ; } ;
  struct finpurge_struct *optvals ;
  lpopts_struct *opts_lb,*opts_dflt,*opts_ub ;

  parse_any result ;

  const char *rtnnme = "lpctl_finpurge" ;


/*
  BNF for the final purge command.
*/
  static tdef(znil,bnfttNIL,NULL,NULL) ;
  static tref(znil_ref,znil,NULL,NULL) ;
  static tdef(zcomma,bnfttD,NULL,",") ;
  static tref(zcomma_ref,zcomma,NULL,NULL) ;

  static tdef(zpurge,bnfttID,NULL,"purge") ;
  static tref(zpurge_ref,zpurge,bnfmin,NULL) ;
  static tdef(zvar,bnfttID,NULL,"variables") ;
  static tref(zvar_ref,zvar,bnfmin,NULL) ;
  static tdef(zcon,bnfttID,NULL,"constraints") ;
  static tref(zcon_ref,zcon,bnfmin,NULL) ;

  static idef(ziminus1,-1) ;
  static idef(ziTRUE,TRUE) ;
  static idef(ziFALSE,FALSE) ;
  static tdef(zTRUE,bnfttID,NULL,"TRUE") ;
  static tref(zTRUE_ref,zTRUE,bnfmin,NULL) ;
  static tdef(zFALSE,bnfttID,NULL,"FALSE") ;
  static tref(zFALSE_ref,zFALSE,bnfmin,NULL) ;

  static iref(zwhatvar_ziTRUE,ziTRUE,mkoff(struct finpurge_struct,vars)) ;
  static iref(zwhatvar_ziFALSE,ziFALSE,mkoff(struct finpurge_struct,vars)) ;
  static comphd(zwhat_alt1) = { compcnt(3),
      mkcref(zvar_ref),mkcref(zTRUE_ref),mkcref(zwhatvar_ziTRUE) } ;
  static comphd(zwhat_alt2) = { compcnt(3),
      mkcref(zvar_ref),mkcref(zFALSE_ref),mkcref(zwhatvar_ziFALSE) } ;
  static iref(zwhatcon_ziTRUE,ziTRUE,mkoff(struct finpurge_struct,cons)) ;
  static iref(zwhatcon_ziFALSE,ziFALSE,mkoff(struct finpurge_struct,cons)) ;
  static comphd(zwhat_alt3) = { compcnt(3),
      mkcref(zcon_ref),mkcref(zFALSE_ref),mkcref(zwhatcon_ziFALSE) } ;
  static comphd(zwhat_alt4) = { compcnt(3),
      mkcref(zcon_ref),mkcref(zTRUE_ref),mkcref(zwhatcon_ziTRUE) } ;
  static comphd(zwhat_alt5) = { compcnt(1), mkcref(znil_ref) } ;
  static althd(zwhat_alts) = { altcnt(5),
      mkaref(zwhat_alt1),mkaref(zwhat_alt2),mkaref(zwhat_alt3),
      mkaref(zwhat_alt4),mkaref(zwhat_alt5) } ;
  static npdef(zwhat,zwhat_alts) ;

  static npref(zfinpurge_what,zwhat,bnflst,zcomma_ref) ;

  static iref(zfinpurge_varminus1,ziminus1,
	      mkoff(struct finpurge_struct,vars)) ;
  static iref(zfinpurge_conminus1,ziminus1,
	      mkoff(struct finpurge_struct,cons)) ;
  static comphd(zfinpurge_alt) = { compcnt(4),mkcref(zfinpurge_varminus1),
      mkcref(zfinpurge_conminus1),mkcref(zpurge_ref),mkcref(zfinpurge_what) } ;
  static gdef(zfinpurge_def,
	      sizeof(struct finpurge_struct),NULL,zfinpurge_alt) ;
  static gref(zfinpurge,zfinpurge_def,bnfdebug,NULL,NULLP) ;

  dy_exposeOptDefaults(&opts_lb,&opts_dflt,&opts_ub) ;
/*
  Now to work.  Initialise the reader, attempt to parse the command, and
  check for errors.
*/
  rdrinit() ;
  if (parse(dy_cmdchn,&zfinpurge,&result) == FALSE)
  { rdrclear() ;
    errmsg(240,rtnnme,"zfinspec") ;
    return (FALSE) ; }
  optvals = (struct finpurge_struct *) result.g ;
  rdrclear() ;
/*
  Process the results. An empty command is taken as a request for the default
  value.
*/
  if (optvals->vars < 0 && optvals->cons < 0)
  { warn(246,rtnnme,"final variable deactivation",
	 (opts_dflt->finpurge.vars == TRUE)?"true":"false") ;
    warn(246,rtnnme,"final constraint deactivation",
	 (opts_dflt->finpurge.cons == TRUE)?"true":"false") ;
    return (TRUE) ; }
/*
  Something is specified. Check variables, then constraints.
*/
  if (optvals->vars >= 0)
  { main_lpopts->finpurge.vars = (bool) optvals->vars ;
    dyio_outfmt(dy_logchn,dy_gtxecho,"variables %s",
	        (main_lpopts->finpurge.vars == TRUE)?"true":"false") ; }
  if (optvals->cons >= 0)
  { main_lpopts->finpurge.cons = (bool) optvals->cons ;
    if (optvals->vars >= 0) dyio_outfmt(dy_logchn,dy_gtxecho,", ") ;
    dyio_outfmt(dy_logchn,dy_gtxecho,"constraints %s",
	        (main_lpopts->finpurge.cons == TRUE)?"true":"false") ; }

  FREE(optvals) ;

  return (TRUE) ; }



static bool lpctl_load  (void)

/*
  This routine processes the 'load' subcommand, which sets values that
  control how dylp populates the constraint system in a cold start.

  The bnf for the load subcommand is:
    <loadopt> ::= lpcontrol load [<fraction>] [<interval-LIST(,)>] ;
    <fraction> ::= <float>
    <interval> ::= <open-delim> <ub> <lb> <close-delim>
    <ub> ::= <float>
    <lb> ::= <float>
    <open-delim> ::= [ | )
    <close-delim> ::= [ | )
  By the time this routine is called, `lpcontrol load' is already parsed.

  If only the fraction is specified, the interval defaults are used. Similarly,
  if there's only an interval specification, the fraction default remains. But
  if only one interval is specified, the second interval is marked invalid. If
  you want two, you have to specify two.

  Parameters: none

  Returns: TRUE if the remainder of the command parses without error, FALSE
	   otherwise.
*/

{ struct interval_struct { struct interval_struct *nxt ;
			   char ldelim ;
			   double ub ;
			   double lb ;
			   char rdelim ; } ;
  struct interval_struct *intv,*temp ;
  struct load_struct { int frac_valid ;
		       double frac ;
		       struct interval_struct *intervals ; } ;
  struct load_struct *loadspec ;
  lpopts_struct *opts_lb,*opts_dflt,*opts_ub ;

  char intvstr[50] ;
  int intvndx,intvlen ;

  parse_any result ;

  const char *rtnnme = "lpctl_load" ;


/*
  BNF for the load command.
*/
  static tdef(znil,bnfttNIL,NULL,NULL) ;
  static tref(znil_ref,znil,NULL,NULL) ;
  static tdef(zcomma,bnfttD,NULL,",") ;
  static tref(zcomma_ref,zcomma,NULL,NULL) ;
  static tdef(zlsq,bnfttD,NULL,"[") ;
  static tref(zlsq_ref,zlsq,NULL,NULL) ;
  static tdef(zlpar,bnfttD,NULL,"(") ;
  static tref(zlpar_ref,zlpar,NULL,NULL) ;
  static tdef(zrsq,bnfttD,NULL,"]") ;
  static tref(zrsq_ref,zrsq,NULL,NULL) ;
  static tdef(zrpar,bnfttD,NULL,")") ;
  static tref(zrpar_ref,zrpar,NULL,NULL) ;
  static tdef(zdnum,bnfttN,10,NULL) ;

  static comphd(zldelim_alt1) = { compcnt(1),mkcref(zlsq_ref) } ;
  static comphd(zldelim_alt2) = { compcnt(1),mkcref(zlpar_ref) } ;
  static althd(zldelim_alts) = { altcnt(2),
    mkaref(zldelim_alt1),mkaref(zldelim_alt2) } ;
  static pdef(zldelim,zldelim_alts) ;

  static comphd(zrdelim_alt1) = { compcnt(1),mkcref(zrsq_ref) } ;
  static comphd(zrdelim_alt2) = { compcnt(1),mkcref(zrpar_ref) } ;
  static althd(zrdelim_alts) = { altcnt(2),
    mkaref(zrdelim_alt1),mkaref(zrdelim_alt2) } ;
  static pdef(zrdelim,zrdelim_alts) ;

  static pref(zinterval_ldelim,zldelim,bnfstore|bnfexact,
	      mkoff(struct interval_struct,ldelim),NULLP) ;
  static tref(zinterval_ub,zdnum,bnfstore|bnfdbl,
	      mkoff(struct interval_struct,ub)) ;
  static tref(zinterval_lb,zdnum,bnfstore|bnfdbl,
	      mkoff(struct interval_struct,lb)) ;
  static pref(zinterval_rdelim,zrdelim,bnfstore|bnfexact,
	      mkoff(struct interval_struct,rdelim),NULLP) ;
  static comphd(zinterval_alt) = { compcnt(4),
    mkcref(zinterval_ldelim),
    mkcref(zinterval_ub),mkcref(zinterval_lb),
    mkcref(zinterval_rdelim) } ;

  static gdef(zinterval,sizeof(struct interval_struct),
	      mkoff(struct interval_struct,nxt),zinterval_alt) ;

  static tref(zfraction_dnum,zdnum,bnfstore|bnfdbl,
	      mkoff(struct load_struct,frac)) ;
  static idef(ziTRUE,TRUE) ;
  static iref(zfrac_valid,ziTRUE,mkoff(struct load_struct,frac_valid)) ;
  static comphd(zfraction_alt1) = { compcnt(2),
    mkcref(zfraction_dnum),mkcref(zfrac_valid) } ;
  static comphd(zfraction_alt2) = { compcnt(1),mkcref(znil_ref) } ;
  static althd(zfraction_alts) = { altcnt(2),
    mkaref(zfraction_alt1),mkaref(zfraction_alt2) } ;
  static pdef(zfraction,zfraction_alts) ;

  static pref(zloadbody_fraction,zfraction,NULL,NULL,NULLP) ;
  static gref(zloadbody_interval,zinterval,bnfstore|bnflst,
	      mkoff(struct load_struct,intervals),zcomma_ref) ;

  static comphd(zloadbody_alt1) = { compcnt(2),
    mkcref(zloadbody_fraction),mkcref(zloadbody_interval) } ;
  static comphd(zloadbody_alt2) = { compcnt(1),mkcref(zloadbody_fraction) } ;

  static althd(zloadbody_alts) = { altcnt(2),
    mkaref(zloadbody_alt1),mkaref(zloadbody_alt2) } ;
  static npdef(zloadbody,zloadbody_alts) ;
  static npref(zloadbody_ref,zloadbody,NULL,NULLP) ;

  static comphd(zload_alt) = { compcnt(1),mkcref(zloadbody_ref) } ;
  static gdef(zload_def,sizeof(struct load_struct),NULL,zload_alt) ;
  static gref(zload,zload_def,NULL,NULL,NULLP) ;

  dy_exposeOptDefaults(&opts_lb,&opts_dflt,&opts_ub) ;
/*
  Now to work.  Initialise the reader, attempt to parse the command, and
  check for errors.
*/
  rdrinit() ;
  if (parse(dy_cmdchn,&zload,&result) == FALSE)
  { rdrclear() ;
    errmsg(240,rtnnme,"zload") ;
    return (FALSE) ; }
  loadspec = (struct load_struct *) result.g ;
  rdrclear() ;
/*
  Process the results. An empty command is taken as a request for the default
  value.
*/
  if (loadspec->frac_valid == FALSE && loadspec->intervals == NULL)
  { warn(245,rtnnme,"initial load fraction",opts_dflt->initcons.frac) ;
    intvndx = 0 ;
    intvlen = sizeof(intvstr)-1 ;
    if (opts_dflt->initcons.i1uopen == TRUE)
      intvstr[intvndx] = '(' ;
    else
      intvstr[intvndx] = '[' ;
    intvndx++ ;
    intvndx += dyio_outfxd(&intvstr[intvndx],-(intvlen-intvndx),
			   'l',"%.5f %.5f",
			   opts_dflt->initcons.i1u,opts_dflt->initcons.i1l) ;
    if (opts_dflt->initcons.i1lopen == TRUE)
      intvstr[intvndx] = ')' ;
    else
      intvstr[intvndx] = ']' ;
    intvndx++ ;
    if (opts_dflt->initcons.i2valid == TRUE)
    { if (opts_dflt->initcons.i2uopen == TRUE)
	intvstr[intvndx] = '(' ;
      else
	intvstr[intvndx] = '[' ;
      intvndx++ ;
      intvndx += dyio_outfxd(&intvstr[intvndx],-(intvlen-intvndx),
			     'l',"%.5f %.5f",
			     opts_dflt->initcons.i2u,opts_dflt->initcons.i2l) ;
      if (opts_dflt->initcons.i2lopen == TRUE)
	intvstr[intvndx] = ')' ;
      else
	intvstr[intvndx] = ']' ;
      intvndx++ ; }
    intvstr[intvndx] = '\0' ;
    warn(246,rtnnme,"load interval",intvstr) ;
    return (TRUE) ; }
/*
  The load fraction first.
*/
  if (loadspec->frac_valid == TRUE)
  { dyio_outfmt(dy_logchn,dy_gtxecho," %.2f",loadspec->frac) ;
    if (loadspec->frac < opts_lb->initcons.frac ||
	loadspec->frac > opts_ub->initcons.frac)
    { warn(244,rtnnme,opts_lb->initcons.frac,"initial load fraction",
	   opts_ub->initcons.frac,loadspec->frac,opts_ub->initcons.frac) ;
      main_lpopts->initcons.frac = opts_ub->initcons.frac ; }
    else
    { main_lpopts->initcons.frac = loadspec->frac ; } }
/*
  Now process the intervals, if present. For each interval, check the values
  and load main_lpopts->initcons.
*/
  intv = loadspec->intervals ;
  FREE(loadspec) ;
  if (intv == NULL) return (TRUE) ;

  dyio_outfmt(dy_logchn,dy_gtxecho," %c %.5f %.5f %c",
	      intv->ldelim,intv->ub,intv->lb,intv->rdelim) ;
  if (intv->ldelim == '(')
    main_lpopts->initcons.i1uopen = TRUE ;
  else
    main_lpopts->initcons.i1uopen = FALSE ;
  if (intv->ub > opts_ub->initcons.i1u || intv->ub < opts_lb->initcons.i1u)
  { warn(244,rtnnme,opts_lb->initcons.i1u,"initial load angle bound",
	   opts_ub->initcons.i1u,intv->ub,opts_ub->initcons.i1u) ;
    main_lpopts->initcons.i1u = opts_ub->initcons.i1u ; }
  else
  { main_lpopts->initcons.i1u = intv->ub ; }
  if (intv->lb > opts_ub->initcons.i1l || intv->lb < opts_lb->initcons.i1l)
  { warn(244,rtnnme,opts_lb->initcons.i1l,"initial load angle bound",
	   opts_ub->initcons.i1l,intv->lb,opts_lb->initcons.i1l) ;
    main_lpopts->initcons.i1l = opts_lb->initcons.i1l ; }
  else
  { main_lpopts->initcons.i1l = intv->lb ; }
  if (intv->rdelim == ')')
    main_lpopts->initcons.i1lopen = TRUE ;
  else
    main_lpopts->initcons.i1lopen = FALSE ;
/*
  If we have only one interval, we're done, otherwise repeat the whole thing
  for the second interval.
*/
  temp = intv->nxt ;
  FREE(intv) ;
  if (temp == NULL)
  { main_lpopts->initcons.i2valid = FALSE ;
    return (TRUE) ; }

  intv = temp ;
  dyio_outfmt(dy_logchn,dy_gtxecho,", %c %.5f %.5f %c",
	      intv->ldelim,intv->ub,intv->lb,intv->rdelim) ;
  if (intv->ldelim == '(')
    main_lpopts->initcons.i2uopen = TRUE ;
  else
    main_lpopts->initcons.i2uopen = FALSE ;
  if (intv->ub > opts_ub->initcons.i2u || intv->ub < opts_lb->initcons.i2u)
  { warn(244,rtnnme,opts_lb->initcons.i2u,"initial load angle bound",
	   opts_ub->initcons.i2u,intv->ub,opts_ub->initcons.i2u) ;
    main_lpopts->initcons.i2u = opts_ub->initcons.i2u ; }
  else
  { main_lpopts->initcons.i2u = intv->ub ; }
  if (intv->lb > opts_ub->initcons.i2l || intv->lb < opts_lb->initcons.i2l)
  { warn(244,rtnnme,opts_lb->initcons.i2l,"initial load angle bound",
	   opts_ub->initcons.i2l,intv->lb,opts_lb->initcons.i2l) ;
    main_lpopts->initcons.i2l = opts_lb->initcons.i2l ; }
  else
  { main_lpopts->initcons.i2l = intv->lb ; }
  if (intv->rdelim == ')')
    main_lpopts->initcons.i2lopen = TRUE ;
  else
    main_lpopts->initcons.i2lopen = FALSE ;

  FREE(intv) ;

  return (TRUE) ; }



static bool lpctl_infinity (void)

/*
  This routine handles the `lpcontrol infinity' command. The reason it's broken
  out is so that we can use the strings IEEE and DBL_MAX to specify the most
  common values for infinite and finite infinity. The syntax will also accept
  an arbitrary (double) value.

  The bnf for the infinity subcommand is:
    <loadopt> ::= lpcontrol infinity [<value>]
    <value> ::= IEEE | DBL_MAX | <double>
  By the time this routine is called, `lpcontrol infinity' is already parsed.
  A value of 0 is taken as a request for the current value.

  Parameters: none

  Returns: TRUE if the remainder of the command parses without error, FALSE
	   otherwise.
*/

{ struct infinity_struct { int code ;
			   double val ; } ;
  struct infinity_struct *infinityspec ;
  lptols_struct *tols_dflt ;

  double infinity ;

  parse_any result ;

  const char *rtnnme = "lpctl_infinity" ;


/*
  BNF for the infinity command.
*/
  static tdef(zIEEE,bnfttID,NULL,"IEEE") ;
  static tref(zIEEE_ref,zIEEE,bnfmin,NULL) ;
  static idef(ziOne,1) ;
  static tdef(zDBLMAX,bnfttID,NULL,"DBL_MAX") ;
  static tref(zDBLMAX_ref,zDBLMAX,bnfmin,NULL) ;
  static idef(ziTwo,2) ;
  static tdef(zdnum,bnfttN,10,NULL) ;
  static idef(ziThree,3) ;

  static iref(zinfbody_IEEE,ziOne,mkoff(struct infinity_struct,code)) ;
  static comphd(zinfbody_alt1) = { compcnt(2),
    mkcref(zIEEE_ref),mkcref(zinfbody_IEEE) } ;

  static iref(zinfbody_DBLMAX,ziTwo,mkoff(struct infinity_struct,code)) ;
  static comphd(zinfbody_alt2) = { compcnt(2),
    mkcref(zDBLMAX_ref),mkcref(zinfbody_DBLMAX) } ;

  static iref(zinfbody_dbl,ziThree,mkoff(struct infinity_struct,code)) ;
  static tref(zinfbody_dblval,zdnum,bnfstore|bnfdbl,
	      mkoff(struct infinity_struct,val)) ;
  static comphd(zinfbody_alt3) = { compcnt(2),
    mkcref(zinfbody_dblval),mkcref(zinfbody_dbl) } ;

  static althd(zinfbody_alts) = { altcnt(3),
    mkaref(zinfbody_alt1),mkaref(zinfbody_alt2), mkaref(zinfbody_alt3) } ;
  static npdef(zinfbody,zinfbody_alts) ;
  static npref(zinfbody_ref,zinfbody,NULL,NULLP) ;

  static comphd(zinfinitydef_alt) = { compcnt(1),mkcref(zinfbody_ref) } ;
  static gdef(zinfinity_def,
	      sizeof(struct infinity_struct),NULL,zinfinitydef_alt) ;
  static gref(zinfinity,zinfinity_def,NULL,NULL,NULLP) ;

  dy_exposeTolDefaults(&tols_dflt) ;
/*
  Now to work.  Initialise the reader, attempt to parse the command, and
  check for errors.
*/
  rdrinit() ;
  if (parse(dy_cmdchn,&zinfinity,&result) == FALSE)
  { rdrclear() ;
    errmsg(240,rtnnme,"zinfinity") ;
    return (FALSE) ; }
  infinityspec = (struct infinity_struct *) result.g ;
  rdrclear() ;
/*
  Process the results. A value of 0 is taken as a request for the default.
  Note that we expect HUGE_VAL to be IEEE infinity. We get the default value
  from main_lptols, because the default isn't set until we call dy_defaults.
*/
  switch (infinityspec->code)
  { case 1:
    { dyio_outfmt(dy_logchn,dy_gtxecho," IEEE (%g)",HUGE_VAL) ;
      infinity = HUGE_VAL ;
      if (finite(infinity))
      { warn(314,rtnnme,infinity) ; }
      break ; }
    case 2:
    { dyio_outfmt(dy_logchn,dy_gtxecho," DBL_MAX (%g)",DBL_MAX) ;
      infinity = DBL_MAX ;
      break ; }
    case 3:
    { dyio_outfmt(dy_logchn,dy_gtxecho," %g",infinityspec->val) ;
      infinity = infinityspec->val ;
      if (infinity == 0)
      { infinity = main_lptols->inf ;
	warn(245,rtnnme,"infinity",infinity) ; }
      else
      if (infinity < 0)
      { errmsg(242,rtnnme,infinity,"infinity") ;
	infinity = tols_dflt->inf ; }
      break ; }
    default:
    { FREE(infinityspec) ;
      return (FALSE) ; } }

  main_lptols->inf = infinity ;
  FREE(infinityspec) ;

  return (TRUE) ; }




cmd_retval dy_ctlopt (const char *keywd)

/*
  This routine handles various control options, each taking a single
  parameter (TRUE/FALSE, a number, etc.) If the routine can't understand
  the option, it returns an error.

  The bnf for the control option command is:
    <ctlopt> ::= lpcontrol <what> <value>
    <what> ::= actconlim | actconlvl | active | actvarlim | bogus |
	       check | cold | coldbasis | coldvars | costz |
	       dchk | deactconlvl | degen | degenlite | degenpivs | dfeas |
	       dualacttype | dualmultipiv |
	       factor | final | forcecopy | fullsys |
	       groom | idle | infinity | iters | load |
	       patch | pchk | pfeas | pivot | primmultipiv |
	       purgecon | purgevar | reframe |
	       scaling | scan | swing | usedual | warm | zero
    <value> ::= <bool> | <integer> | <double> | <string>
  
  The tolerances bogus, costz, dchk, pchk, pivot, purgecon, purgevar, reframe,
  swing, and zero expect a double. Infinity is a bit more complicated, and
  further parsing is handled by its own routine.

  The options cold, degen, forcecopy, fullsys, patch, usedual, and warm expect
  a boolean. Coldbasis, deactconlvl, degenlite, and groom expect
  a string from the appropriate list of defined keywords (see below), and the
  others expect an integer.

  Active, final, and load are more complicated, and are handled by private
  parsing routines.

  Parameters:
    keywd:	The command keyword

  Returns: cmdOK or cmdHALTERROR
*/

{ char *what,cmdstr[50] ;
  int code,intdflt,intlb,intub,*intopt,numkwds ;
  double *toler,tolerdflt ;
  double dblopt ;
  bool booldflt,*boolopt ;
  keytab_entry *kwds ;
  cmd_retval retval ;
  lpopts_struct *opts_lb,*opts_dflt,*opts_ub ;
  lptols_struct *tols_dflt ;
  const char *rtnnme = "dy_ctlopt" ;

/*
  A lookup table with the various <what> keywords recognised by the control
  command.
*/

  enum ctlcodes { ctlINV = 0, ctlACTIVESZE, ctlADDVARLIM,
		  ctlBOGUS,
		  ctlCHECK, ctlCOLD, ctlCOLDBASIS, ctlCOLDVARS,
		  ctlCONACTLIM, ctlCONACTLVL, ctlCONDEACTLVL, ctlCONTEXT,
		  ctlCPYORIG, ctlCOSTZ,
		  ctlDCHK, ctlDEGEN, ctlDEGENLITE, ctlDEGENPIVS, ctlDFEAS,
		  ctlDUALADD, ctlDUALMULTIPIV,
		  ctlFACTOR, ctlFINAL, ctlFULLSYS, ctlGROOM,
		  ctlINFINITY, ctlITER, ctlIDLE, ctlLOAD,
		  ctlPATCH, ctlPCHK, ctlPFEAS,
		  ctlPIVOT, ctlPRIMMULTIPIV,
		  ctlPURGE, ctlPURGEVAR, ctlREFRAME,
		  ctlSCALING, ctlSCAN, ctlSWING,
		  ctlUSEDUAL, ctlWARM, ctlZERO } ctlcode ;

  static keytab_entry ctlkwds[] = { { "actconlim", 8, (int) ctlCONACTLIM },
				    { "actconlvl", 8, (int) ctlCONACTLVL },
				    { "active", 4, (int) ctlACTIVESZE },
				    { "actvarlim", 4, (int) ctlADDVARLIM },
				    { "antidegen", 2, (int) ctlDEGEN },
				    { "bogus", 1, (int) ctlBOGUS },
				    { "check", 2, (int) ctlCHECK },
				    { "cold", 5, (int) ctlCOLD },
				    { "coldbasis", 5, (int) ctlCOLDBASIS },
				    { "coldvars", 5, (int) ctlCOLDVARS },
				    { "context", 3, (int) ctlCONTEXT },
				    { "costz", 3, (int) ctlCOSTZ },
				    { "dchk", 2, (int) ctlDCHK },
				    { "deactconlvl", 3, (int) ctlCONDEACTLVL },
				    { "degenlite", 6, (int) ctlDEGENLITE },
				    { "degenpivs", 6, (int) ctlDEGENPIVS },
				    { "dfeas", 2, (int) ctlDFEAS },
				    { "dualacttype", 5, (int) ctlDUALADD },
				    { "dualmultipiv", 5,
				      (int) ctlDUALMULTIPIV },
				    { "factor", 2, (int) ctlFACTOR },
				    { "final", 2, (int) ctlFINAL },
				    { "forcecopy", 2, (int) ctlCPYORIG },
				    { "fullsys", 2, (int) ctlFULLSYS },
				    { "groom", 1, (int) ctlGROOM },
				    { "idle", 2, (int) ctlIDLE },
				    { "infinity", 2, (int) ctlINFINITY },
				    { "iters", 2, (int) ctlITER },
				    { "load", 1, (int) ctlLOAD },
				    { "patch", 2, (int) ctlPATCH },
				    { "pchk", 2, (int) ctlPCHK },
				    { "pfeas", 2, (int) ctlPFEAS },
				    { "pivot", 2, (int) ctlPIVOT },
				    { "primmultipiv", 5,
				      (int) ctlPRIMMULTIPIV },
				    { "purgecon", 6, (int) ctlPURGE },
				    { "purgevar", 6, (int) ctlPURGEVAR },
				    { "reframe", 1, (int) ctlREFRAME },
				    { "scaling", 4, (int) ctlSCALING },
				    { "scan", 4, (int) ctlSCAN },
				    { "swing", 2, (int) ctlSWING },
				    { "usedual", 1, (int) ctlUSEDUAL },
				    { "warm", 1, (int) ctlWARM },
				    { "zero", 1, (int) ctlZERO }
				   } ;

  static int numctlkwds = (sizeof ctlkwds/sizeof (keytab_entry)) ;

/*
  Keywords used by the groom option; the option value is used only by the
  groombasis routine, so there's no real motivation for symbolic codes.
*/

  static keytab_entry groomkwds[] = { { "abort", 1, 2 },
				      { "silent", 1, 0 },
				      { "warn", 1, 1 } } ;

  static int numgroomkwds = (sizeof groomkwds/sizeof (keytab_entry)) ;

/*
  Keywords used by the coldbasis option.
*/

  static keytab_entry basiskwds[] = { { "architectural", 1, ibARCH },
				      { "logical", 1, ibLOGICAL },
				      { "slack", 1, ibSLACK } } ;

  static int numbasiskwds = (sizeof basiskwds/sizeof (keytab_entry)) ;

/*
  Keywords used by the context option
*/

  static keytab_entry contextkwds[] = { { "single", 1, cxSINGLELP },
					{ "initial", 1, cxINITIALLP },
					{ "bandc", 1, cxBANDC } } ;

  static int numcontextkwds = (sizeof contextkwds/sizeof (keytab_entry)) ;

/*
  Keywords used by the degenlite option; the option value is used only by
  the primalout and dualin routines, so there's no real motivation for
  symbolic codes.
*/

  static keytab_entry litekwds[] = { { "alignedge", 6, 3 },
				     { "alignobj", 6, 2 },
				     { "perpedge", 5, 5 },
				     { "perpobj", 5, 4 },
				     { "pivot", 6, 1 },
				     { "pivotabort", 6, 0 } } ;

  static int numlitekwds = (sizeof litekwds/sizeof (keytab_entry)) ;

/*
  Keywords used by the deactconlvl option; the option value is used only by
  the dy_purgecon routine, so there's no real motivation for symbolic codes.
*/
  
  static keytab_entry deactkwds[] = { { "aggressive", 1, 1 },
				      { "fanatic", 1, 2 },
				      { "normal", 1, 0 } } ;

  static int numdeactkwds = (sizeof deactkwds/sizeof (keytab_entry)) ;

/*
  Initialisation, so gcc doesn't complain.
*/
  kwds = NULL ;
  numkwds = 0 ;
  intdflt = -INT_MAX ;
  intlb = INT_MAX ;
  intub = -INT_MAX ;
  intopt = NULL ;
  toler = NULL ;
  tolerdflt = quiet_nan(0) ;
  boolopt = NULL ;
/*
  Acquire the option and tolerance limits and defaults.
*/
  dy_exposeOptDefaults(&opts_lb,&opts_dflt,&opts_ub) ;
  dy_exposeTolDefaults(&tols_dflt) ;
/*
  Now to work. Parse off the <what> keyword and see if we can look it up.
*/
  ctlcode = ctlINV ;
  if (string_opt(&what) == TRUE)
  { code = ambig(what,ctlkwds,numctlkwds) ;
    if (code < 0) 
    { if (code < -1)
	errmsg(233,rtnnme,what) ;
      else
	errmsg(234,rtnnme,what) ; }
    else
      ctlcode = (enum ctlcodes) code ; }
/*
  Set the various variables for each command. There are a few exceptions to the
  pattern, down at the end of the switch, which call command-specific routines
  for more complicated processing.
*/
  dyio_outfxd(cmdstr,-((int) (sizeof(cmdstr)-1)),'l',"%s %s",keywd,what) ;
  switch (ctlcode)
  { case ctlADDVARLIM:
    { intopt = &main_lpopts->addvar ;
      intdflt = opts_dflt->addvar ;
      intlb = opts_lb->addvar ;
      intub = opts_ub->addvar ;
      break ; }
    case ctlBOGUS:
    { toler = &main_lptols->bogus ;
      tolerdflt = tols_dflt->bogus ;
      break ; }
    case ctlCHECK:
    { intopt = &main_lpopts->check ;
      intdflt = opts_dflt->check ;
      intlb = opts_lb->check ;
      intub = opts_ub->check ;
      break ; }
    case ctlCOLD:
    { boolopt = &main_lpopts->forcecold ;
      booldflt = opts_dflt->forcecold ;
      break ; }
    case ctlCOLDBASIS:
    { intopt = (int *) &main_lpopts->coldbasis ;
      intdflt = opts_dflt->coldbasis ;
      intlb = opts_lb->coldbasis ;
      intub = opts_ub->coldbasis ;
      numkwds = numbasiskwds ;
      kwds = basiskwds ;
      break ; }
    case ctlCOLDVARS:
    { intopt = &main_lpopts->coldvars ;
      intdflt = opts_dflt->coldvars ;
      intlb = opts_lb->coldvars ;
      intub = opts_ub->coldvars ;
      break ; }
    case ctlCONACTLIM:
    { intopt = &main_lpopts->con.actlim ;
      intdflt = opts_dflt->con.actlim ;
      intlb = opts_lb->con.actlim ;
      intub = opts_ub->con.actlim ;
      break ; }
    case ctlCONACTLVL:
    { intopt = &main_lpopts->con.actlvl ;
      intdflt = opts_dflt->con.actlvl ;
      intlb = opts_lb->con.actlvl ;
      intub = opts_ub->con.actlvl ;
      break ; }
    case ctlCONTEXT:
    { intopt = (int *) &main_lpopts->context ;
      intdflt = opts_dflt->context ;
      intlb = opts_lb->context ;
      intub = opts_ub->context ;
      numkwds = numcontextkwds ;
      kwds = contextkwds ;
      break ; }
    case ctlCOSTZ:
    { toler = &main_lptols->cost ;
      tolerdflt = tols_dflt->cost ;
      break ; }
    case ctlCPYORIG:
    { boolopt = &main_lpopts->copyorigsys ;
      booldflt = opts_dflt->copyorigsys ;
      break ; }
    case ctlDCHK:
    { toler = &main_lptols->dchk ;
      tolerdflt = tols_dflt->dchk ;
      break ; }
    case ctlDEGEN:
    { boolopt = &main_lpopts->degen ;
      booldflt = opts_dflt->degen ;
      break ; }
    case ctlDEGENLITE:
    { intopt = &main_lpopts->degenlite ;
      intdflt = opts_dflt->degenlite ;
      intlb = opts_lb->degenlite ;
      intub = opts_ub->degenlite ;
      numkwds = numlitekwds ;
      kwds = litekwds ;
      break ; }
    case ctlDEGENPIVS:
    { intopt = &main_lpopts->degenpivlim ;
      intdflt = opts_dflt->degenpivlim ;
      intlb = opts_lb->degenpivlim ;
      intub = opts_ub->degenpivlim ;
      break ; }
    case ctlDFEAS:
    { toler = &main_lptols->dfeas_scale ;
      tolerdflt = tols_dflt->dfeas_scale ;
      break ; }
    case ctlDUALADD:
    { intopt = &main_lpopts->dualadd ;
      intdflt = opts_dflt->dualadd ;
      intlb = opts_lb->dualadd ;
      intub = opts_ub->dualadd ;
      break ; }
    case ctlFACTOR:
    { intopt = &main_lpopts->factor ;
      intdflt = opts_dflt->factor ;
      intlb = opts_lb->factor ;
      intub = opts_ub->factor ;
      break ; }
    case ctlFULLSYS:
    { boolopt = &main_lpopts->fullsys ;
      booldflt = opts_dflt->fullsys ;
      break ; }
    case ctlGROOM:
    { intopt = &main_lpopts->groom ;
      intdflt = opts_dflt->groom ;
      intlb = opts_lb->groom ;
      intub = opts_ub->groom ;
      numkwds = numgroomkwds ;
      kwds = groomkwds ;
      break ; }
    case ctlITER:
    { intopt = &main_lpopts->iterlim ;
      intdflt = opts_dflt->iterlim ;
      intlb = opts_lb->iterlim ;
      intub = opts_ub->iterlim ;
      break ; }
    case ctlIDLE:
    { intopt = &main_lpopts->idlelim ;
      intdflt = opts_dflt->idlelim ;
      intlb = opts_lb->idlelim ;
      intub = opts_ub->idlelim ;
      break ; }
    case ctlDUALMULTIPIV:
    { intopt = &main_lpopts->dpsel.strat ;
      intdflt = opts_dflt->dpsel.strat ;
      intlb = opts_lb->dpsel.strat ;
      intub = opts_ub->dpsel.strat ;
      break ; }
    case ctlPRIMMULTIPIV:
    { intopt = &main_lpopts->ppsel.strat ;
      intdflt = opts_dflt->ppsel.strat ;
      intlb = opts_lb->ppsel.strat ;
      intub = opts_ub->ppsel.strat ;
      break ; }
    case ctlPATCH:
    { boolopt = &main_lpopts->patch ;
      booldflt = opts_dflt->patch ;
      break ; }
    case ctlPCHK:
    { toler = &main_lptols->pchk ;
      tolerdflt = tols_dflt->pchk ;
      break ; }
    case ctlPFEAS:
    { toler = &main_lptols->pfeas_scale ;
      tolerdflt = tols_dflt->pfeas_scale ;
      break ; }
    case ctlPIVOT:
    { toler = &main_lptols->pivot ;
      tolerdflt = tols_dflt->pivot ;
      break ; }
    case ctlPURGE:
    { toler = &main_lptols->purge ;
      tolerdflt = tols_dflt->purge ;
      break ; }
    case ctlCONDEACTLVL:
    { intopt = &main_lpopts->con.deactlvl ;
      intdflt = opts_dflt->con.deactlvl ;
      intlb = opts_lb->con.deactlvl ;
      intub = opts_ub->con.deactlvl ;
      numkwds = numdeactkwds ;
      kwds = deactkwds ;
      break ; }
    case ctlPURGEVAR:
    { toler = &main_lptols->purgevar ;
      tolerdflt = tols_dflt->purgevar ;
      break ; }
    case ctlREFRAME:
    { toler = &main_lptols->reframe ;
      tolerdflt = tols_dflt->reframe ;
      break ; }
    case ctlSCALING:
    { intopt = &main_lpopts->scaling ;
      intdflt = opts_dflt->scaling ;
      intlb = opts_lb->scaling ;
      intub = opts_ub->scaling ;
      break ; }
    case ctlSCAN:
    { intopt = &main_lpopts->scan ;
      intdflt = opts_dflt->scan ;
      intlb = opts_lb->scan ;
      intub = opts_ub->scan ;
      break ; }
    case ctlSWING:
    { toler = &main_lptols->swing ;
      tolerdflt = tols_dflt->swing ;
      break ; }
    case ctlUSEDUAL:
    { boolopt = &main_lpopts->usedual ;
      booldflt = opts_dflt->usedual ;
      break ; }
    case ctlWARM:
    { boolopt = &main_lpopts->forcewarm ;
      booldflt = opts_dflt->forcewarm ;
      break ; }
    case ctlZERO:
    { toler = &main_lptols->zero ;
      tolerdflt = tols_dflt->zero ;
      break ; }
    case ctlACTIVESZE:
    { booldflt = lpctl_active() ;
      if (booldflt == TRUE)
      { return (cmdOK) ; }
      else
      { return (cmdHALTERROR) ; } }
    case ctlINFINITY:
    { booldflt = lpctl_infinity() ;
      if (booldflt == TRUE)
      { return (cmdOK) ; }
      else
      { return (cmdHALTERROR) ; } }
    case ctlLOAD:
    { booldflt = lpctl_load() ;
      if (booldflt == TRUE)
      { return (cmdOK) ; }
      else
      { return (cmdHALTERROR) ; } }
    case ctlFINAL:
    { booldflt = lpctl_finpurge() ;
      if (booldflt == TRUE)
      { return (cmdOK) ; }
      else
      { return (cmdHALTERROR) ; } }
    default:
    { errmsg(236,rtnnme,"<what>","keyword",keywd) ;
      STRFREE(what) ;
      return (cmdHALTERROR) ; } }
/*
  Last but not least, the actual work. For an integer option or a tolerance,
  a negative value is taken as a request from the user to be told the default
  value. Further, for tolerances, 0 is not acceptable. For coldvars and factor,
  gripe about a violation of the upper bound, but don't enforce it. An upper
  bound of -1 is `no upper bound' for the integer-valued options.
  
  Assume the worst, so that we have to explicitly set a successful return code.
*/
  retval = cmdHALTERROR ;
  switch (ctlcode)
  { case ctlADDVARLIM:
    case ctlCHECK:
    case ctlCONACTLIM:
    case ctlCONACTLVL:
    case ctlCONTEXT:
    case ctlDEGENPIVS:
    case ctlDUALADD:
    case ctlDUALMULTIPIV:
    case ctlITER:
    case ctlIDLE:
    case ctlPRIMMULTIPIV:
    case ctlSCALING:
    case ctlSCAN:
    { if (integer_opt(intopt) == TRUE)
      { if (*intopt >= 0)
	{ if (intub > 0 && *intopt > intub)
	  { warn(241,rtnnme,intlb,cmdstr,intub,*intopt,intub) ;
	    *intopt = intub ; } }
	else
	{ warn(243,rtnnme,cmdstr,intdflt) ; }
	retval = cmdOK ; }
      else
      { errmsg(236,rtnnme,"<integer>","parameter",keywd) ; }
      break ; }
    case ctlCOLDVARS:
    case ctlFACTOR:
    { if (integer_opt(intopt) == TRUE)
      { if (*intopt >= 0)
	{ if (intub > 0 && *intopt > intub)
	  { warn(241,rtnnme,intlb,cmdstr,intub,*intopt,intub) ; } }
	else
	{ warn(243,rtnnme,cmdstr,intdflt) ; }
	retval = cmdOK ; }
      else
      { errmsg(236,rtnnme,"<integer>","parameter",keywd) ; }
      break ; }
    case ctlCOLD:
    case ctlCPYORIG:
    case ctlDEGEN:
    case ctlFULLSYS:
    case ctlPATCH:
    case ctlUSEDUAL:
    case ctlWARM:
    { if (bool_opt(boolopt) == TRUE)
      { retval = cmdOK ; }
      else
      { errmsg(236,rtnnme,"<bool>","parameter",keywd) ; }
      break ; }
    case ctlBOGUS:
    case ctlCOSTZ:
    case ctlDCHK:
    case ctlDFEAS:
    case ctlPCHK:
    case ctlPFEAS:
    case ctlPIVOT:
    case ctlPURGE:
    case ctlPURGEVAR:
    case ctlREFRAME:
    case ctlSWING:
    case ctlZERO:
    { if (double_opt(&dblopt) == TRUE)
      { if (dblopt <= 0)
	{ warn(245,rtnnme,cmdstr,tolerdflt) ; }
	else
	{ *toler = dblopt ; }
	retval = cmdOK ; }
      else
      { errmsg(236,rtnnme,"<double>","parameter",keywd) ; }
      break ; }
    case ctlGROOM:
    case ctlCOLDBASIS:
    case ctlDEGENLITE:
    case ctlCONDEACTLVL:
    { if (string_opt(&what) == TRUE)
      { code = ambig(what,kwds,numkwds) ;
	if (code < 0)
	{ if (code < -1)
	    errmsg(233,rtnnme,what) ;
	  else
	    errmsg(234,rtnnme,what) ; }
	else
	{ *intopt = code ;
	  retval = cmdOK ; } }
      else
      { errmsg(236,rtnnme,"<string>","parameter",cmdstr) ; }
      break ; }
    default:
    { errmsg(1,rtnnme,__LINE__) ;
      break ; } }

  STRFREE(what) ;

  return (retval) ; }


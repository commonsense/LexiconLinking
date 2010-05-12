/*
  This file is a portion of the Dylp LP distribution.

        Copyright (C) 2004 Lou Hafer, Stephen Tse
        Copyright (C) 2005 -- 2007 Lou Hafer

        School of Computing Science
        Simon Fraser University
        Burnaby, B.C., V5A 1S6, Canada
        lou@cs.sfu.ca

  This code is licensed under the terms of the Common Public License (CPL).
*/

/*
  This file contains a main program and support routines sufficient to build
  a stand-alone lp solver using the version of dylp as modified for COIN/OSI.
  It's not necessarily pretty, but it should work adequately well.
  
  After processing any command line options, the program establishes
  basic facilities, reads and processes the options file, reads and
  processes the problem specification, calls dylp via a wrapper to solve
  the problem, prints the result and shuts down.

  The command line is expected to be of the form

  osi_dylp [<options>] [<problem-file>]
  
  The options presently in place are:

  -s			Run silent: turns off echo of all generated text to
			  stdout. The default output-file path is changed
			  from stdout to NULL. Silent overpowers terse, in
			  the event both are specified.
  -t			Terse output on stdout. Behaviour is as for silent,
			  but allows an opening title and closing message
			  giving the result of the LP.
  -p <num>		Set overall print level to <num>, [0..5].
  -e <errmsg-file>	Source text for error messages (defaults to
			  dy_errmsgs.txt).
  -E <errlog-file>	A logging file for error messages (default is to
			  direct error messages to stderr and the log file).
  -o <option-file>	Control ('.spc') options for dylp. Default is to not
			look for an options file. Disabled on Windows.
  -m <problem-file>	The problem ('.mps') specification. Defaults to stdin.
  -L <log-file>		A log of dylp's execution (default is no execution
			  logging).
  -O <output-file>	The output file. Defaults to stdout unless the -s
			  option is present, in which case the default is no
			  output.

  -h			Print this help message and exit.
  -v			Print version and exit.

  The -m option is just an alternate way to specify the <problem-file>.

  The error log file is a duplicate of the error messages printed on stderr;
  the execution log file is a duplicate of output on stdout unless -s is
  specified to suppress output to stdout.
*/

#include <string.h>

/*
  We're getting close to the O/S here. Microsoft does it differently, of
  course. There's some disagreement in the unix community about what goes
  where, but I'm hoping that sys/time.h and sys/resource.h together will
  cover it.
*/
#include <time.h>
#if defined(_MSC_VER) || defined(__MSVCRT__)
  /* Nothing to do here. */
#else
# include <sys/time.h>
# include <sys/resource.h>
#endif

#include "dy_cmdint.h"
#include "dylp.h"

static char sccsid[] UNUSED = "@(#)osi_dylp.c	1.10	09/25/04" ;
static char svnid[] UNUSED = "$Id: osi_dylp.c 218 2008-03-28 15:57:22Z lou $" ;

const char *osidylp_version = "1.10" ;			/* sccs! */
const char *osidylp_time ;

/*
  Macro cleverness to specify a default error message file. Depends on ANSI
  C merge of consecutive string constants. DYLP_ERRMSGDIR should have the
  form "path/to/distribution/DyLP/src/Dylp/", including the quotes. See
  DyLP/src/DylpStdLib/DylpConfig.h for further information.
*/

#ifdef DYLP_ERRMSGDIR
# define DYLP_ERRMSGPATH DYLP_ERRMSGDIR "dy_errmsgs.txt"
#else
# define DYLP_ERRMSGPATH "dy_errmsgs.txt"
#endif

/*
  Variables which control i/o operations.

  ttyout		i/o id for output to the user's terminal
  ttyin			i/o id for input from the user's terminal
  dy_cmdchn		i/o id for input from the command file
  dy_logchn		i/o id for the execution log file

  dy_cmdecho		controls echoing of command input to stdout
  dy_gtxecho		controls echoing of generated text to stdout
*/

ioid dy_cmdchn,dy_logchn ;
bool dy_cmdecho, dy_gtxecho ;

/*
  File paths used elsewhere in osi_dylp.

  outpath:	Output file; the optimal solution is written here, as
		well as intermediate incumbent solutions (if requested).
		Intermediate incumbent solutions are solutions which improve
		on the previous solution (hence become the incumbent). Used
		in fathom:check_incumbent
*/

const char *outpath ;

consys_struct *main_sys ;
lpprob_struct *main_lp ;
lpopts_struct *main_lpopts ;
lptols_struct *main_lptols ;



/*
  Routines to produce the version and help messages. The code looks better
  if these are off to the side.
*/

static void print_version (ioid chn, bool echo,
			   const char *cmd, const char *nme, const char *ver)
/*
  Print version, copyright, and free software disclaimer.
*/

{ const char *disc1 = "This is free software; see the source for copying"
		      " conditions. There is NO" ;
  const char *disc2 = "warranty; not even for MERCHANTABILITY or FITNESS FOR A"
		      " PARTICULAR PURPOSE." ;

  dyio_outfmt(dy_logchn,dy_gtxecho,"\n%s (%s) V %s",cmd,nme,ver) ;
  dyio_outfmt(dy_logchn,dy_gtxecho,"\nCopyright (C) 2004 Lou Hafer") ;
  dyio_outfmt(dy_logchn,dy_gtxecho,"\n%s\n%s\n",disc1,disc2) ;

  return ; }

static void print_help (ioid chn, bool echo, const char *name)
/*
  Print help message.
*/
{ dyio_outfmt(chn,echo,"\nusage: %s [<options>] [<problem-file>]",name) ;

  dyio_outfmt(chn,echo,"\n\nThe options presently in place are:\n") ;

  dyio_outfmt(chn,echo,"\n  %s\t\t\t%s","-s",
	      "Run silent: turns off echo of all generated text to") ;
  dyio_outfmt(chn,echo,"\n\t\t\t%s",
	      "stdout. The default output-file path is changed") ;
  dyio_outfmt(chn,echo,"\n\t\t\t%s",
	      "from stdout to NULL. Silent overpowers terse, in") ;
  dyio_outfmt(chn,echo,"\n\t\t\t%s","the event both are specified.") ;

  dyio_outfmt(chn,echo,"\n  %s\t\t\t%s","-t",
	      "Terse output on stdout. Behaviour is as for silent,") ;
  dyio_outfmt(chn,echo,"\n\t\t\t%s",
	      "but allows an opening title and closing message") ;
  dyio_outfmt(chn,echo,"\n\t\t\t%s","giving the result of the LP.") ;

  dyio_outfmt(chn,echo,"\n  %s\t\t%s","-p <num>",
  	      "Set overall print level to <num>, [0..5].") ;

  dyio_outfmt(chn,echo,"\n  %s\t%s","-e <errmsg-file>",
	      "Source text for error messages (defaults to") ;
  dyio_outfmt(chn,echo,"\n\t\t\t%s","dy_errmsgs.txt).") ;

  dyio_outfmt(chn,echo,"\n  %s\t%s","-E <errlog-file>",
	       "A logging file for error messages (default is to") ;
  dyio_outfmt(chn,echo,
	      "\n\t\t\t%s","direct error messages to the log file).") ;

  dyio_outfmt(chn,echo,"\n  %s\t%s","-o <option-file>",
# if defined(_MSC_VER) || defined(__MSVCRT__)
	      "Disabled on Windows.") ;
# else
	      "Control ('.spc') options file (default is no file).") ;
# endif

  dyio_outfmt(chn,echo,"\n  %s\t%s","-m <problem-file>",
	      "The problem ('.mps') specification (no default).") ;

  dyio_outfmt(chn,echo,"\n  %s\t\t%s","-L <log-file>",
	      "A log of execution (default is no execution") ;
  dyio_outfmt(chn,echo,"\n\t\t\t%s","logging).") ;

  dyio_outfmt(chn,echo,"\n  %s\t%s","-O <output-file>",
	      "The output file. Defaults to stdout unless the -s or -t") ;
  dyio_outfmt(chn,echo,"\n\t\t\t%s",
	      "options are present, in which case the default is no") ;
  dyio_outfmt(chn,echo,"\n\t\t\t%s","output.") ;

  dyio_outfmt(chn,echo,"\n  %s\t\t\t%s","-h",
	      "Print this help message and exit.") ;
  dyio_outfmt(chn,echo,"\n  %s\t\t\t%s","-v",
	      "Print version and exit.") ;

  dyio_outfmt(chn,echo,"\n\n%s %s",
	      "The -m option is just an alternate way to specify the",
	      "<problem-file>.") ;

  dyio_outfmt(chn,echo,"\n\n%s\n%s\n%s\n",
  "The error log file is a duplicate of the error messages printed on",
  "stderr; the execution log file is a duplicate of output on stdout unless",
  "-s is specified to suppress output to stdout.") ;

  return ; }



/*
  A few utility routines to produce run time information. We're close to the
  O/S here, everyone has their own idea. `struct timeval' is the unix side of
  life. Microsoft just uses ticks. Microsecond timing is questionable at best,
  so maybe Microsoft has it right here. Easiest way to finesse this (coming
  from the Linux side) is to define a timeval structure. Then the only routine
  we need is get_timeval.
*/

#if defined(_MSC_VER) || defined(__MSVCRT__)

typedef long int __time_t ;
typedef long int __suseconds_t ;

struct timeval
{ __time_t tv_sec ;
  __suseconds_t tv_usec ; } ;

static void get_timeval (struct timeval *tv)
/*
  Call the clock function to get time in ticks, and convert to seconds and
  microseconds.

  Parameter:
    tv:	(i/o) timeval structure, will be filled with the current time

  Returns: undefined
*/

{ double temp ;

  temp = (double) clock() ;
  temp /= CLOCKS_PER_SEC ;
  tv->tv_sec = (__time_t) temp ;
  temp -= tv->tv_sec ;
  tv->tv_usec = (__suseconds_t) (temp*(1e6)) ;

  return ; }

#else	/* end Microsoft (_MSC_VER, __MSVCRT__), start Linux version */
  
static void get_timeval (struct timeval *tv)

/*
  Call getrusage to retrieve the user and system time, sum them, and return
  the lot as a timeval.

  Parameter:
    tv:	(i/o) timeval structure, will be filled with the current time
  
  Returns: undefined
*/

{ struct rusage ruse ;

# ifndef NDEBUG
/*
  Avoid a `read from unitialised' error.
*/
  memset(&ruse,0,sizeof(struct rusage)) ;
# endif

/*
  Retrieve the current time.
*/
  if (getrusage(RUSAGE_SELF,&ruse) < 0)
  { tv->tv_sec = 0 ;
    tv->tv_usec = 0 ;
    return ; }
/*
  Sum up the user and system time and return it. I'm assuming that the number
  of usecs. in the tv_usec field is always less than 1e6 (i.e., less than
  1 sec.) and therefore it's safe to simply add user and system values without
  worrying about integer overflow.
*/
  tv->tv_usec = ruse.ru_utime.tv_usec+ruse.ru_stime.tv_usec ;
  tv->tv_sec = tv->tv_usec/((int) 1e6) ;
  tv->tv_usec = tv->tv_usec%((int) 1e6) ;
  tv->tv_sec += ruse.ru_utime.tv_sec ;
  tv->tv_sec += ruse.ru_stime.tv_sec ;

  return ; }

# endif		/* Microsoft/Linux time function */


static void sub_timevals (struct timeval *diff,
		   struct timeval *end, struct timeval *start)

/*
  Calculate diff = end-start in timeval's.

  Parameters:
    start,end:	start and end timevals.
    diff:	(o) timeval to be filled in with the difference
  
  Returns: undefined
*/

{ 
/*
  Normalise the start and end timevals.
*/
  if (start->tv_usec > 1e6)
  { start->tv_sec += start->tv_usec/((int) 1e6) ;
    start->tv_usec = start->tv_usec%((int) 1e6) ; }
  if (end->tv_usec > 1e6)
  { end->tv_sec += end->tv_usec/((int) 1e6) ;
    end->tv_usec = end->tv_usec%((int) 1e6) ; }
/*
  Do the subtraction. We may have to borrow a few usecs.
*/
  if (end->tv_usec < start->tv_usec)
  { end->tv_usec += (int) 1e6 ;
    end->tv_sec-- ; }
  diff->tv_usec = end->tv_usec-start->tv_usec ;
  diff->tv_sec = end->tv_sec-start->tv_sec ;

  return ; }


static void prt_timeval (ioid chn, bool echo, struct timeval *tv)

/*
  A little local utility to print a timeval as hh:mm:ss.ss.

  Parameters:
    chn:	output channel id
    echo:	controls echo to stdout
    tv:		timeval to be printed
  
  Returns: undefined
*/

{ int hrs, mins, secs, csecs ;

  hrs = tv->tv_sec/3600 ;
  mins = tv->tv_sec%3600 ;
  secs = mins%60 ;
  mins = mins/60 ;
  csecs = tv->tv_usec/1e4 ;
  dyio_outfmt(chn,echo,"%d:%02d:%02d.%02d",hrs,mins,secs,csecs) ;

  return ; }



static void stats_lp (const char *outpath, bool echo, lpprob_struct *lp,
		      struct timeval *lptime, lpstats_struct *lpstats)

/*
  A little shell routine to handle writing detailed statistics on an LP to
  the output file.
  
  Parameters:
    outpath:	the output file path name.
    echo:	TRUE to echo to stdout, FALSE otherwise
    lp:	  	lp problem structure
    lptime:	elapsed time for call to do_lp
    lpstats:	lp statistics structure
 
  Returns : undefined
*/

{ ioid chn ;
  int vndx,bpos ;
  const char *rtnnme = "stats_lp" ;

/*
  Set up the output. Don't echo this to stdout twice.
*/
  if (outpath == NULL)
  { warn(2,rtnnme,"file name") ;
    chn = IOID_NOSTRM ; }
  else
  { chn = dyio_pathtoid(outpath,NULL) ;
    if (chn == IOID_INV) chn = dyio_openfile(outpath,"w") ;
    if (chn == IOID_INV)
    { warn(10,rtnnme,outpath,"w") ;
      chn = IOID_NOSTRM ; }
    if (strcmp(outpath,"stdout") == 0) echo = FALSE ; }
/*
  Print a few items from the lp structure --- name, status, pivot count, and
  lp return code.
*/
  if (lp == NULL)
  { dyio_outfmt(chn,echo,
	        "\n\n<< %s: LP problem structure is NULL! >>\n", rtnnme) ; }
  else
  { dyio_outfmt(chn,echo,
		"\n\nSystem: %s\t\t\tfinal status: %s after %d iterations.",
	        lp->consys->nme,dy_prtlpphase(lp->phase,FALSE),lp->iters) ;
    if (lp->phase == dyDONE)
    { dyio_outfmt(chn,echo,"\n    lp status: %s",dy_prtlpret(lp->lpret)) ;
      switch (lp->lpret)
      { case lpOPTIMAL:
	{ dyio_outfmt(chn,echo,"\t\tobjective: %.9g",lp->obj) ;
	  break ; }
	case lpINFEAS:
	{ dyio_outfmt(chn,echo,"\t\tinfeasibility: %.9g",lp->obj) ;
	  break ; }
	case lpUNBOUNDED:
	{ if (lp->obj != 0)
	  { if (lp->obj < 0)
	    { vndx = abs((int) lp->obj) ;
	      bpos = -1 ; }
	    else
	    { vndx = (int) lp->obj ;
	      bpos = 1 ; }
	    dyio_outfmt(chn,echo,"\t\tunbounded variable %s (%d) (%s)",
		        consys_nme(lp->consys,'v',vndx,FALSE,NULL),vndx,
		        (bpos < 0)?"decreasing":"increasing") ; }
	  break ; }
	default:
	{ break ; } } }
    if (lptime != NULL)
    { dyio_outfmt(chn,echo,"\n    lp time: ") ;
      prt_timeval(chn,echo,lptime) ;
      dyio_outfmt(chn,echo," (%.2f)",lptime->tv_sec+lptime->tv_usec/1e6) ; } }
# ifdef DYLP_STATISTICS
  if (lpstats != NULL) dy_dumpstats(chn,echo,lpstats,lp->consys) ;
# endif
  
  dyio_outfmt(chn,echo,"\n") ;

  dyio_flushio(chn,echo) ;

  return ; }



static lpret_enum do_lp_all (struct timeval *elapsed, int printlvl)

/*
  This routine is a wrapper which makes up the difference between the
  setup completed in the main program and the setup expected by dylp_dolp.
  Structures are created for the lp problem, options, tolerances, and
  statistics, and passed to dylp_dolp. Once the lp completes, the stats are
  printed and the statistics structure is released. Other cleanup is handled
  in main.

  Returns: lp status code; lpINV is used in the event of a fatal error that's
	   not really dylp's fault.
*/

{ int ndx,i,flips ;
  bool *flipped ;

  lpret_enum lpret ;
  lpopts_struct *initial_lpopts ;
  lpstats_struct *initial_lpstats ;
  basis_struct *basis ;

  struct timeval diff,before,after ;

  const char *rtnnme = "do_lp_all" ;

  lpret = lpINV ;

/*
  Create and initialise the main_lp structure to pass the problem in to dylp.
*/
  main_lp = (lpprob_struct *) CALLOC(1,sizeof(lpprob_struct)) ;
  main_lp->phase = dyINV ;
  main_lp->consys = main_sys ;
  main_lp->rowsze = main_sys->rowsze ;
  main_lp->colsze = main_sys->colsze ;
/*
  Step through the constraints and replace ax >= b constraints with (-a)x <=
  -b constraints. consys_mulrow will take care of the necessary
  modifications. Normally this would be handled in mpsin, along with the
  deletion of empty constraints, but the OSI test suite isn't tolerant of
  changing the sense of constraints, let alone removing a few. Arguably a
  good thing.
*/
  flipped = (bool *) CALLOC(main_sys->concnt+1,sizeof(bool)) ;
  flips = 0 ;
  for (ndx = main_sys->concnt ; ndx > 0 ; ndx--)
  { if (main_sys->ctyp[ndx] == contypGE)
    { if (consys_mulrow(main_sys,ndx,-1.0) == FALSE)
      { errmsg(112,rtnnme,main_sys->nme,"scalar multiply","row",
	       consys_nme(main_sys,'c',ndx,FALSE,NULL),ndx) ;
	FREE(flipped) ;
	return (lpFATAL) ; }
      flipped[ndx] = TRUE ;
      flips++ ; } }
/*
  Set up options, statistics, and solve the lp. After we're done, dump the
  statistics and free the dylp statistics structure. The options specified
  here (cold start using a logical basis and the full constraint system) are
  appropriate for the initial solution of an LP.
*/
  initial_lpopts = (lpopts_struct *) MALLOC(sizeof(lpopts_struct)) ;
  memcpy(initial_lpopts,main_lpopts,sizeof(lpopts_struct)) ;

  dy_setprintopts(printlvl,initial_lpopts) ;

# ifdef DYLP_STATISTICS
  initial_lpstats = NULL ;
  dy_initstats(&initial_lpstats,main_sys) ;
# else
  initial_lpstats = NULL ;
# endif

  get_timeval(&before) ;
  lpret = dylp(main_lp,initial_lpopts,main_lptols,initial_lpstats) ;
  get_timeval(&after) ;
  sub_timevals(&diff,&after,&before) ;
  if (elapsed != NULL) (*elapsed) = diff ;

  stats_lp(outpath,FALSE,main_lp,&diff,initial_lpstats) ;
# ifdef DYLP_STATISTICS
  dy_freestats(&initial_lpstats) ;
# endif

/*
  Time to undo any constraint flips. We also have to tweak the corresponding
  duals --- flipping the sign of a row in the basis corresponds to flipping the
  sign of a column in the basis inverse, which means that the sign of the
  corresponding dual is flipped. This requires a little care --- if we're in
  dynamic mode, some constraints may be inactive.
*/
  if (flips > 0)
  { for (ndx = main_sys->concnt ; ndx > 0 ; ndx--)
    { if (flipped[ndx] == TRUE)
      { if (consys_mulrow(main_sys,ndx,-1.0) == FALSE)
	{ errmsg(112,rtnnme,main_sys->nme,"scalar multiply","row",
		 consys_nme(main_sys,'c',ndx,FALSE,NULL),ndx) ;
	  FREE(flipped) ;
	  return (lpFATAL) ; } } }
    if (main_lp->y != NULL)
    { basis = main_lp->basis ;
      for (ndx = 0 ; ndx < basis->len ; ndx++)
      { i = basis->el[ndx].cndx ;
	if (flipped[i] == TRUE) main_lp->y[ndx] = -main_lp->y[ndx] ; } } }
  FREE(flipped) ;

  if (initial_lpopts != NULL) FREE(initial_lpopts) ;

  return (lpret) ;
}



int main (int argc, char *argv[])

{ time_t timeval ;
  struct tm *tm ;
  char runtime[50] ;
  ioid ttyin,ttyout,outchn ;
  int optlett,printlvl ;
  bool silent,terse,swaperrs,errecho,doversion,dohelp ;
  const char *errmsgpath,*errlogpath,*optpath,*mpspath,*logpath ;

  struct timeval lptime ;
  
  const char *rtnnme = argv[0] ;

/*
  getopt() --- we need explicit declarations under strict ANSI compliance
  mode. This function seems to be quite common, however. Present in GCC.
*/
  int getopt(int argc, char*const *argv, const char *optstring) ;
  extern char *optarg ;
  extern int opterr,optind,optopt ;

  const char *optstring = "e:stp:E:o:m:L:O:hv";

  extern struct tm *localtime(const time_t *clock) ;

  /* mpsio.c */

  extern bool dy_mpsin(const char *mpspath,
  		       consys_struct **consys, double infinity);

  /* dy_basis.c */

  extern void dy_initbasis(int concnt, int factor_freq, double zero_tol),
	      dy_freebasis(void) ;

/*
  Set up some defaults, then process the command line options. This is all very
  specific to Unix and SunOS.
*/
  errmsgpath = DYLP_ERRMSGPATH ;
  errlogpath = NULL ;
  optpath = NULL ;
  mpspath = NULL ;
  outpath = NULL ;
  logpath = NULL ;

  ttyout = IOID_INV ;
  ttyin = IOID_INV ;
  dy_logchn = IOID_INV ;
  outchn = IOID_INV ;

  silent = FALSE ;
  terse = FALSE ;
  dy_gtxecho = TRUE ;
  dy_cmdecho = FALSE ;
  doversion = FALSE ;
  dohelp = FALSE ;

  printlvl = -1 ;

  opterr = 0 ;
  for (optlett = getopt(argc,argv,optstring) ;
       optlett != -1 ;
       optlett = getopt(argc,argv,optstring))
    switch (optlett)
    { case 'e':
      { errmsgpath = optarg ;
	break ; }
      case 'E':
      { errlogpath = optarg ;
	break ; }
      case 'o':
      { optpath = optarg ;
	break ; }
      case 'm':
      { mpspath = optarg ;
	break ; }
      case 'L':
      { logpath = optarg ;
	break ; }
      case 'O':
      { outpath = optarg ;
	break ; }
      case 's':
      { silent = TRUE ;
	dy_gtxecho = FALSE ;
	break ; }
      case 't':
      { terse = TRUE ;
	dy_gtxecho = FALSE ;
	break ; }
      case 'p':
      { printlvl = atoi(optarg) ;
	break ; }
      case 'v':
      { doversion = TRUE ;
	break ; }
      case 'h':
      { dohelp = TRUE ;
	break ; }
      case '?':
      { errinit(errmsgpath,errlogpath,TRUE) ;
	dyio_ioinit() ;
	errmsg(3,rtnnme,"command line option",optopt) ;
	exit (1) ; } }
/*
  If there's still a parameter left, it must be the mps file, specified without
  using the -m option. There should not be more than one parameter remaining.
*/
  if (optind < argc)
  { mpspath = argv[optind] ;
    optind++ ; }
  if (optind < argc)
  { dohelp = TRUE ; } 
/*
  What's our output level? If the user has specified a print level, go with it.
  Otherwise, take a cue from any -s or -t flags. Default to print level 2.
*/
  if (printlvl < 0)
  { if (silent == TRUE)
      printlvl = 0 ;
    else
    if (terse == TRUE)
      printlvl = 1 ;
    else
      printlvl = 2 ; }
/*
  Output file name: if the user hasn't specified one, and we're not running
  silent, default to stdout.
*/
  if (outpath == NULL && silent == FALSE) outpath = "stdout" ;
/*
  Grab the time and format it nicely.
*/
  if (time(&timeval) == (time_t)(-1))
  { warn(19,rtnnme) ;
    osidylp_time = "n/a" ; }
  else
  { tm = localtime(&timeval) ;
    if (tm != NULL)
    { strftime(runtime,sizeof(runtime),"%A, %B %d, %Y,  %I:%M:%S %p",tm) ;
      osidylp_time = runtime ; }
    else
    { osidylp_time = "n/a" ; } }
/*
  Figure out the appropriate settings for silent and terse.
  silent is set to (silent || terse), and is passed to process_cmds so that
	 it can properly handle dy_cmdecho and dy_gtxecho.
  terse	 is set for ease of controlling the output specifically mentioned in
	 conjunction with terse mode (which is controlled from this routine).
	 The proper value is (terse || !silent).
  The cryptic little if statement below accomplishes this (try a decision tree
  to convince yourself).
*/
  if (terse == TRUE)
  { if (silent == TRUE)
      terse = FALSE ;
    else
      silent = TRUE ; }
  else
  { if (silent == FALSE) terse = TRUE ; }
/*
  Execute initialization routines for the i/o and error reporting packages.
*/
  if (silent == TRUE)
    errecho = FALSE ;
  else
    errecho = TRUE ;
  errinit(errmsgpath,errlogpath,errecho) ;
  if (dyio_ioinit() != TRUE)
  { errmsg(1,rtnnme,__LINE__) ;
    exit (2) ; }
/*
  Connect ttyout to the standard output. Initialize ttyin, setting the mode to
  line-oriented. Serious internal confusion if we can't manage these. Set the
  initial command input channel to stdin.
*/
  ttyout = dyio_openfile("stdout","w") ;
  if (ttyout == IOID_INV)
  { errmsg(1,rtnnme,__LINE__) ;
    exit(3) ; }
  ttyin = dyio_openfile("stdin","r") ;
  if (ttyin == IOID_INV)
  { errmsg(1,rtnnme,__LINE__) ;
    exit(4) ; }
  (void) dyio_setmode(ttyin,'l') ;
  dy_cmdchn = ttyin ;
/*
  Initialize logging.
*/
  if (logpath != NULL)
  { dy_logchn = dyio_openfile(logpath,"w") ;
    if (dy_logchn == IOID_INV)
    { warn(201,rtnnme,logpath) ;
      dy_logchn = IOID_NOSTRM ; } }
  else
  { dy_logchn = IOID_NOSTRM ; }
/*
  Are we supposed to merge the error messages with the log stream? (Note that
  errors will be echoed to stderr unless we're running silent. If the user's
  turned off logging too, well, it's their business.) swaperrs just tells the
  code that it should reset the error logging channel if it ever resets the
  main logging channel.
*/
  if (errlogpath == NULL && dy_logchn != IOID_NOSTRM)
  { swaperrs = TRUE ;
    errlogpath = logpath ;
    if (dyio_chgerrlog(errlogpath,errecho) == FALSE)
    { warn(18,rtnnme,"<null>",errlogpath) ; } }
/*
  Ok, after all that work, check if we've been asked for the version or usage
  messages. If so, do it and head for the exit. Version preempts help.
*/
  if (doversion == TRUE)
  { print_version(dy_logchn,dy_gtxecho,argv[0],rtnnme,osidylp_version) ;
    goto NOOUTFILE_CLEANUP ; }
  if (dohelp == TRUE)
  { print_help(dy_logchn,dy_gtxecho,argv[0]) ;
    goto NOOUTFILE_CLEANUP ; }
/*
  We're up! Banners to the appropriate places.
*/
  dyio_outfmt(dy_logchn,terse,"\n\t\t    %s\tV %s\n",rtnnme,osidylp_version) ;
  dyio_outfmt(dy_logchn,terse,"\n\t\t%s",runtime) ;
  dyio_outfmt(dy_logchn,terse,"\n\n") ;
  if (outpath != NULL && strcmp(outpath,"stdout") != 0)
  { outchn = dyio_pathtoid(outpath,NULL) ;
    if (outchn == IOID_INV) outchn = dyio_openfile(outpath,"w") ;
    if (outchn == IOID_INV)
    { warn(10,rtnnme,outpath,"w") ; }
    else
    { dyio_outfmt(outchn,FALSE,"\n\t\t    %s\tV %s\n",rtnnme,osidylp_version) ;
      dyio_outfmt(outchn,FALSE,"\n\t\t%s",runtime) ;
      dyio_outfmt(outchn,FALSE,"\n\n") ; } }
/*
  Time to set up the lp options. Establish a set of defaults, then read the
  options file to see what the user has in mind. Because this is a standalone
  shell, doing a one-time solution for an LP, set up a default of cold start
  using the full system and a logical basis. This can be overridden in a .spc
  file if desired.

  For reasons that escape me at the moment, the parser fails on Windows. This
  may get fixed eventually. For now, disabled by the simple expedient of
  forcing optpath to NULL.
*/
  dy_defaults(&main_lpopts,&main_lptols) ;
  main_lpopts->forcecold = TRUE ;
  main_lpopts->fullsys = TRUE ;
  main_lpopts->coldbasis = ibLOGICAL ;
# if defined(_MSC_VER) || defined(__MSVCRT__)
  optpath = NULL ;
# endif
  if (optpath != NULL)
  { dy_cmdchn = dyio_openfile(optpath,"r") ;
    if (dy_cmdchn == IOID_INV) exit (1) ;
    (void) dyio_setmode(dy_cmdchn,'l') ;
    switch (process_cmds(silent))
    { case cmdOK:
      { break ; }
      case cmdHALTERROR:
      { exit (1) ; }
      case cmdHALTNOERROR:
      { exit (0) ; }
      default:
      { exit (1) ; } }
    if (dy_cmdchn != ttyin)
    { (void) dyio_closefile(dy_cmdchn) ;
      dy_cmdchn = IOID_INV ; } }
/*
  Make an attempt to read the mps input file.
*/
  if (dy_mpsin(mpspath,&main_sys,main_lptols->inf) == FALSE)
  { exit (10) ;}
/*
  Check over the option settings, now that we know how big the constraint
  system will be.
*/
  dy_checkdefaults(main_sys,main_lpopts,main_lptols) ;
/*
  Initialise the basis maintenance package. The second parameter controls how
  many basis updates the basis can hold before it requires refactoring.
  Adding 5 to dylp's refactor interval should give a safety margin.
*/
  dy_initbasis(2*main_sys->concnt,main_lpopts->factor+5,0.0) ;
/*
  Run the lp.
*/
  if (do_lp_all(&lptime,printlvl) == FALSE)
  { errmsg(443,rtnnme,main_sys->nme,dy_prtlpphase(main_lp->phase,TRUE),
	   main_lp->iters) ; }
/*
  Should we produce any output? Print to a file, if requested.
*/
  if (outchn != IOID_INV && outchn != ttyout)
  { dy_dumpcompact(outchn,FALSE,main_lp,FALSE) ; }
/*
  Any final terminal output we should do?
*/
  if (printlvl >= 1)
  { if (printlvl >= 2)
    { dy_dumpcompact(dy_logchn,
		     (dy_logchn == IOID_INV)?TRUE:FALSE,main_lp,FALSE) ; }
    dyio_outfmt(dy_logchn,TRUE,"\nReturn code %s",dy_prtlpret(main_lp->lpret)) ;
    if (main_lp->phase == dyDONE)
      dyio_outfmt(dy_logchn,TRUE," after %d pivots",main_lp->iters) ;
    if (main_lp->lpret == lpOPTIMAL)
    { dyio_outfmt(dy_logchn,TRUE,"; objective %.8g",main_lp->obj) ; }
    dyio_outfmt(dy_logchn,TRUE,
		" (%.2f sec.)",lptime.tv_sec+lptime.tv_usec/1e6) ;
    dyio_outfmt(dy_logchn,TRUE,".\n") ;
    dyio_flushio(dy_logchn,dy_gtxecho) ; }
/*
  Final cleanup. Free space used by the remaining main_* structures.
*/
  dy_freebasis() ;

  if (main_lp != NULL)
  { dy_freesoln(main_lp) ;
    if (main_lp->consys != NULL) consys_free(main_lp->consys) ;
    FREE(main_lp) ; }
  if (main_lpopts != NULL) FREE(main_lpopts) ;
  if (main_lptols != NULL) FREE(main_lptols) ;
/*
  Leap to here for shutdown when we opened the output file but didn't solve
  an LP.
*/
  if (outchn != IOID_INV && outchn != ttyout)
  { (void) dyio_closefile(outchn) ; }
/*
  Leap to here for shutdown in cases where we never get as far as opening
  an output file. We still need to close the log file and shut down the i/o
  subsystem.
*/
  NOOUTFILE_CLEANUP:
  if (dy_logchn != IOID_INV && dy_logchn != IOID_NOSTRM)
  { (void) dyio_closefile(dy_logchn) ; }
  dyio_ioterm() ;
  errterm() ;

  exit(0) ;
  
/*
  Just to suppress the silly warning from the compiler, which isn't satisfied
  with the immediately preceding `exit'.
*/

  return (0) ; }

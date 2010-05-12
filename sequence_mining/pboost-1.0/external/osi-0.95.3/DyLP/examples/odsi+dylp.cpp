/*
  This file is a portion of the OsiDylp distribution.

        Copyright (C) 2004 Lou Hafer

        School of Computing Science
        Simon Fraser University
        Burnaby, B.C., V5A 1S6, Canada
        lou@cs.sfu.ca

  This program is free software ; you can redistribute it and/or modify it
  under the terms of the GNU General Public License as published by the Free
  Software Foundation ; either version 2 of the License, or (at your option)
  any later version.

  This program is distributed in the hope that it will be useful, but WITHOUT
  ANY WARRANTY ; without even the implied warranty of MERCHANTABILITY or
  FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for
  more details.

  You should have received a copy of the GNU General Public License along
  with this program ; if not, write to the Free Software Foundation, Inc., 59
  Temple Place, Suite 330, Boston, MA  02111-1307  USA
*/

/*
  This file contains a main program and support routines sufficient to build
  a stand-alone lp solver using dylp working through OsiDylpSolverInterface.
  
  After processing any command line options, the program will deal with the
  .spc file, if specified, then call ODSI::readMps to read the MPS file.
  ODSI::initialSolve is called to solve the LP, and then the results are
  printed.

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
  -o <option-file>	Control ('.spc') options for dylp. Disabled on
  			Windows.
  -m <problem-file>	The problem ('.mps') specification.
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


/*! \file odsi+dylp.cpp
    \brief An alternate main program for running dylp by way of the
    OsiDylpSolverInterface. No frills. Good for testing when you suspect
    the problem lies not in dylp but in the OsiDylpSolverInterface.

    Requires that the COIN libraries be available.
*/

#include <iostream>
#include <assert.h>
#include "dylib_std.h"
#include "OsiDylpSolverInterface.hpp"
#include "CoinTime.hpp"
#include "CoinError.hpp"

/* Begin unnamed local namespace */

namespace {

char sccsid[] UNUSED = "@(#)odsi+dylp.cpp	1.5	11/06/04" ;
char svnid[] UNUSED = "$Id: odsi+dylp.cpp 96 2006-07-02 05:49:51Z lou $" ;

const char *osidylp_version = static_cast<const char *>("1.10") ;

using std::string ;

/*! \brief Print the version message. */

void print_version (const char *cmd, const char *nme, const char *ver)
/*
  Print version, copyright, and free software disclaimer.
*/

{ std::cout << std::endl << cmd << " (" << nme << ") V " << ver ;
  std::cout << std::endl << "Copyright (C) 2004 Lou Hafer" ;
  std::cout << std::endl
	    << "This is free software; see the source for copying"
	    << " conditions. There is NO"
	    << std::endl
	    << "warranty; not even for MERCHANTABILITY or FITNESS FOR A"
	    << " PARTICULAR PURPOSE."
	    << std::endl ;

  return ; }


/*! \brief Print the help message */

void print_help (const char *pgm)

/*
  Print help message.
*/

{ std::cout << "\nusage: " << pgm << " [<options>] [<problem-file>]" ;

  std::cout << "\n\nThe options presently in place are:\n" ;

  std::cout << "\n  -s\t\t\t"
	<< "Run silent: turns off echo of all generated text to"
	<< "\n\t\t\t"
	<< "stdout. The default output-file path is changed"
	<< "\n\t\t\t"
	<< "from stdout to NULL. Silent overpowers terse, in"
	<< "\n\t\t\t"
	<< "the event both are specified." ;

  std::cout << "\n  -t\t\t\t"
	<< "Terse output on stdout. Behaviour is as for silent,"
	<< "\n\t\t\t"
	<< "but allows an opening title and closing message"
	<< "\n\t\t\t"
	<< "giving the result of the LP." ;

  std::cout << "\n  -p <num>\t\t"
  	<< "Set overall print level to <num>, [0..5]." ;

  std::cout << "\n  -o <option-file>\t"
#if defined(_MSC_VER) || defined(__MSVCRT__)
	<< "Disabled on Windows." ;
#else
	<< "Control ('.spc') options file for dylp (default is"
	<< "\n\t\t\t"
	<< "no file)." ;
#endif

  std::cout << "\n  -m <problem-file>\t"
	<< "The problem ('.mps') specification (no default)." ;

  std::cout << "\n  -L <log-file>\t\t"
	<< "A log of dylp's execution (default is no execution"
	<< "\n\t\t\t"
	<< "log)." ;

  std::cout << "\n  -O <output-file>\t"
	<< "The output file. Defaults to stdout unless the -s or"
	<< "\n\t\t\t"
	<< "-t options are present, in which case the default is"
	<< "\n\t\t\t"
	<< "no output." ;

  std::cout << "\n  -h\t\t\t"
	<< "Print this help message and exit." ;
  std::cout << "\n  -v\t\t\t"
	<< "Print version and exit." ;

  std::cout << "\n\n"
	<< "The -m option is just an alternate way to specify the"
	<< " <problem-file>."
	<< std::endl
	<< "The execution log file is a duplicate of output to stdout"
	<< std::endl
	<< "unless -s or -t were specified to suppress output to stdout."
	<< std::endl ;

  return ; }


/*! \brief Routine to dissect a file name.

  Given a path name in the form  [path/prefix/]base[.suffix][.gz], the routine
  will return the prefix, base, and suffix, along with a boolean that tells
  whether the file name indicates compression (.gz suffix).
*/

void dissectPath (const string fullPath, string &prefix, string &base,
		  string &ext, bool &compressed)
/*
  This routine extracts the base file name from a path of the form
      [path/prefix/]base[.suffix][.gz]
  where all parts are optional except the base name.
*/

{ bool gzSeen = false ;
  string::size_type pathpos,sfxpos,gzpos ;

/*
  std::cout << std::endl
	    << "Starting with \"" << fullPath << "\" ("
	    << fullPath.length() << ")\n" ;
*/

/*
  Locate end of leading path prefix, if any. pathpos will point to the start
  of the base name.
*/
  pathpos = fullPath.rfind('/') ;
  if (pathpos == string::npos)
  { pathpos = 0 ;
    prefix = "" ; }
  else
  { prefix = fullPath.substr(0,pathpos) ;
    pathpos++ ; }
/*
  Check for .gz and/or other suffix.
*/
  sfxpos = fullPath.rfind('.') ;
  gzpos = fullPath.length() ;
/*
  We have a valid suffix if sfxpos ends up strictly between pathpos (start of
  base name) and gzpos (end of string). If it's .gz, note we've seen it and
  try again for a standard suffix.
*/
  if (sfxpos > pathpos && sfxpos < gzpos)
  { ext = fullPath.substr(sfxpos) ;
    if (ext == ".gz")
    { gzpos = sfxpos ;
      gzSeen = true ;
      sfxpos = fullPath.rfind('.',gzpos-1) ; } }
  if (sfxpos > pathpos && sfxpos < gzpos)
  { ext = fullPath.substr(sfxpos+1,gzpos-sfxpos-1) ; }
  else
  { ext = "" ;
    sfxpos = gzpos ; }
/*
  The basename should lie between pathpos and sfxpos. Set a last few return
  values and we're done.
*/
  base = fullPath.substr(pathpos,sfxpos-pathpos) ;
  compressed = gzSeen ;

/*
  std::cout << "Parsed \"" << fullPath <<"\" to \""
	    << prefix << "\", \""
	    << base << "\", \"" << ext << "\"" ;
  if (gzSeen == true) std::cout << " (compressed)" ;
  std::cout << "." << std::endl ;
*/

  return ; }


/*! \brief Run an individual problem

  This routine will run a single problem, using an options file if one is
  specified. The routine will create a log file for the run.
*/

void test_user (const char* mpspath, const char* spcpath,
		const char *logpath, const char *outpath,
		bool silent, bool terse, int printlvl)

{ OsiDylpSolverInterface* osi = new OsiDylpSolverInterface ;
  string fullPath ;
  string prefix,base,ext,prefixAndBase ;
  bool compressed ;

/*
  Open thelog file and output file, if requested.
*/
  if (logpath != 0)
  { fullPath = logpath ;
    dissectPath(fullPath,prefix,base,ext,compressed) ;
    if (prefix != "")
      prefixAndBase = prefix+'/' ;
    else
      prefixAndBase = "" ;
    prefixAndBase += base ;
    if (silent == false)
    { std::cout << std::endl
		<< "Logging to " << prefixAndBase << ".log." ; }
    osi->dylp_logfile(prefixAndBase.c_str(),false) ; }

  if (outpath != 0)
  { fullPath = outpath ;
    dissectPath(fullPath,prefix,base,ext,compressed) ;
    if (prefix != "")
      prefixAndBase = prefix+'/' ;
    else
      prefixAndBase = "" ;
    prefixAndBase += base ;
    if (silent == false)
    { std::cout << std::endl
		<< "Output file is " << prefixAndBase << ".out." ; }
    osi->dylp_outfile(prefixAndBase.c_str()) ; }
/*
  Establish output levels. Bits <2:0> of info (info&0x7) are interpreted as
  an integer which is passed to dy_setprintopts. Bit <3> (0x8) will turn on
  CoinMpsIO messages. Bit <4> (0x10) will cause dylp to send output to the
  terminal in addition to any log file. So, for example, info = 0x11 will
  set print level 1 and enable terminal output. At print level 1 and higher,
  you can override individual settings with an options file, processed via
  dylp_controlfile.
*/
  int info = printlvl ;
  if (silent == false && terse == false) info |= 0x10 ;
  osi->setHintParam(OsiDoReducePrint,false,OsiHintTry,&info) ;
/*
  Process the aforementioned options file, if specified.
*/
  if (spcpath != 0)
  { if (silent == false)
    { std::cout << std::endl
		<< "Processing options file " << spcpath << "." ; }
    osi->dylp_controlfile(spcpath,false,true) ; }
/*
  Process the mps file.
*/
  if (silent == false)
  { std::cout << std::endl << "Reading problem file " << mpspath << "." ; }

  fullPath = mpspath ;
  dissectPath(fullPath,prefix,base,ext,compressed) ;
  if (prefix != "")
    prefixAndBase = prefix+'/' ;
  else
    prefixAndBase = "" ;
  prefixAndBase += base ;

  int errs = osi->readMps(prefixAndBase.c_str(),ext.c_str()) ;
  if (errs != 0)
  { std::cout << errs << " errors while reading " << fullPath << "\n" ;
    throw CoinError("MPS input error","test_user","odsi+dylp.cpp") ; }
/*
  And try to solve the LP.
*/
  if (silent == false)
  { std::cout << std::endl << "Starting LP." << std::endl ; }

  try
  { double startTime = CoinCpuTime();
    osi->setHintParam(OsiDoPresolveInInitial,true,OsiHintTry) ;
    osi->initialSolve() ;
    double timeOfSolution = CoinCpuTime()-startTime;
    if (silent == false)
    { std::cout << std::endl << "Finished." ; }
    double val = osi->getObjValue() ;
    int iters = osi->getIterationCount() ;
    std::cout << " obj = " << val
	      << ", iters = " << iters
	      << ", " << timeOfSolution << "secs."
	      << std::endl ;
    osi->dylp_printsoln(true,true) ; }
  catch (CoinError &err)
  { std::cout << std::endl
	      << err.className() << "::" << err.methodName()
	      << " (throw) " << err.message() << std::endl ;
    std::cout.flush() ; }
  
  delete osi ;

  exit(0) ; }

} /* End unnamed local namespace */

/*! \brief An alternate main program for testing

  This is an alternate program for testing dylp and the COIN/OSI layer.
  Given an mps file and (optional) options file on the command line, it will
  run the problem and exit.
*/

int main (int argc, const char* argv[])

{ const char *rtnnme = static_cast<const char *>("dylp") ;

  std::ios::sync_with_stdio() ;
/*
  Process command line options, if any. The user can specify a single,
  unadorned file name as the only argument. The working assumption is that
  if we're looking for the `-' that starts an option, and we don't see it,
  then we're looking at the file name.
*/
  int argNum ;
  bool silent,terse,doversion,dohelp ;
  int printlvl ;
  const char *optpath,*mpspath,*logpath,*outpath ;

  silent = false ;
  terse = false ;
  doversion = false ;
  dohelp = false ;
  printlvl = 1 ;
  optpath = 0 ;
  mpspath = 0 ;
  logpath = 0 ;
  outpath = 0 ;

  for (argNum = 1 ; argNum < argc ; argNum++)
  { char argLett = argv[argNum][0] ;

    if (argLett != '-') break ;
    argLett = argv[argNum][1] ;
/*
  Handle `--' option prefix for the Gnu-ish.
*/
    if (argLett == '-') argLett = argv[argNum][2] ;

    switch (argLett)
    { case 'o':
      { 
        optpath = argv[++argNum] ;
#if defined(_MSC_VER) || defined(__MSVCRT__)
/* Disabled on Windows --- parser fails. */
	optpath = NULL ;
#endif
	break ; }
      case 'm':
      { mpspath = argv[++argNum] ;
	break ; }
      case 'L':
      { logpath = argv[++argNum] ;
	break ; }
      case 'O':
      { outpath = argv[++argNum] ;
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
      { printlvl = atoi(argv[++argNum]) ;
	break ; }
      case 'v':
      { doversion = true ;
	break ; }
      case 'h':
      { dohelp = true ;
	doversion = true ;
	break ; }
      default:
      { std::cout << argv[0] << ": unrecognized option \"-" << argLett
		  << "\"." << std::endl ;
	print_help(argv[0]) ;
	exit (1) ; } }

  if (doversion == true || dohelp == true) break ; }

/*
  If there's exactly one parameter left, it must be the mps file, specified
  without using the -m option. There should be at most one parameter left
  at this point. If we have parameters left, the user is confused. If we don't
  have an mps file, the user is confused.
*/
  if (argNum == argc-1)
  { mpspath = argv[argNum++] ; }
  if (argNum < argc || mpspath == 0)
  { dohelp = TRUE ; } 
/*
  Have we been asked to print help or the version? If so, do it and exit.
*/
  if (doversion == true || dohelp == true)
  { if (doversion == true) print_version(argv[0],rtnnme,osidylp_version) ;
    if (dohelp == true) print_help(argv[0]) ;
    exit (0) ; }
/*
  What's our output level? If the user has specified a print level, go with it.
  Otherwise, take a cue from any -s or -t flags. Default to print level 2.
*/
  if (printlvl < 0)
  { if (silent == true)
      printlvl = 0 ;
    else
    if (terse == true)
      printlvl = 1 ;
    else
      printlvl = 2 ; }
/*
  Call test_user to do all the work.
*/
  try
  { test_user(mpspath,optpath,logpath,outpath,silent,terse,printlvl) ; }
  catch (CoinError &err)
  { std::cout << std::endl
	      << err.className() << "::" << err.methodName()
	      << " (throw) " << err.message() << std::endl ;
    std::cout.flush() ; }

  return (0) ; }


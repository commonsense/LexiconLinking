
! Generic options file specifying copious printing and a few other random
! options. (Comments in an options file start with an `!'.) See the dylp
! printed documentation (Dylp/Doc) for a thorough explanation of options.
! If you're wondering where this comes from, it's an accretion of various
! options that I find handy while debugging  -- lh --

! lpcontrol scaling 0 ;
! lpcontrol dualacttype 3 ;

lpcontrol infinity DBL_MAX ;

lpprint major 1 ;
lpprint phase1 4 ;
lpprint phase2 4 : ;
lpprint dual 4 ;
lpprint pricing 0 ;
lpprint pivoting 1 ;
lpprint pivreject 2 ;
lpprint degen 1 ;
lpprint scaling 2 ;
lpprint basis 5 ;
lpprint crash 1 ;
lpprint setup 2 ;
lpprint conmgmt 1 ;
lpprint varmgmt 1 ;


! lpcontrol antidegen false ;
! lpcontrol iter 550 ;
! lpcontrol idle 100 ;

! lpcontrol final purge variables false, constraints false ;

! lpcontrol load 1.0  [180 90), (90 0] ;
lpcontrol fullsys true ;
lpcontrol coldbasis logical ;
! lpcontrol usedual false ;
! lpcontrol dualmultipiv 0 ;


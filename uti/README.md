This utility takes as input the result of a dynamical matrix
calculation, and produce the spin-modes profiles in silo
format. Then, an interactive visualization program, as VisIt
(https://visit-dav.github.io/visit-website/) can be used to show/plot
the profiles in the desired format.

-------------------------------------------------------------------------------

COMPILATION AND LINKING

Prerequisite:

SILO library (https://wci.llnl.gov/simulation/computer-codes/silo)

=================================

gfortran: 
Use the script compile-and-link-ifor.

=================================

ifort (Intel fortran compiler, not tested): 
Use the script compile-and-link-gfortran.  May be necessary to change 
dmmtosilo.f (replace code for removing/creating directories with the 
PFX... intrinsics of ifort).

-------------------------------------------------------------------------------

USAGE

"dmmtosilo inputfile"

where inputfile is the name of the main input file. Its (fixed) format is
(line by line):


-------------------------------------------------------------------------------

TEST

Use the output produced by micro3d in the test directory,
e.g. test/my-output.dmm.

1. Remove/rename the previous "output" directory; 
2. run "dmmtosilo in-template": this will produce the mode profiles corresponding to ../test/my-output.dmm.  
3. Mode profiles and power spectrum can be found in the "output" directory, in .silo format.

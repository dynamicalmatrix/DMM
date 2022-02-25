# DMM
Dynamical Matrix Method (software for computing spin modes of nanoparticles and magnonic crystals)

This program calculates the spin normal modes of a magnetic particle
(nanodot) according to the "Dynamical Matrix Method".

Those using this code for scientific publications are kindly requested 
to cite [L. Giovannini et al., Phys. Rev. B 70, 172404 (2004)] and
[https://github.com/dynamicalmatrix/DMM.git].

The theory behind this package has been developed by several authors
and published in Refs. [1,2]. This code has been written by me and is
released under the GNU General Public License v3.0 (see LICENSE file).

The code is rather rudimentary and the program only accept input through
files, i.e. lacks a GUI. I cannot answer any question concerning
installation, usage or upgrade of this program for lacking of time,
i.e. this code is unsupported. Anybody wishing to take care of supporting
this package (and hopefully improve it by, e.g., adding a GUI) is welcome.

To contact the developer, send e-mail to <dynamicalmatrix@gmail.com>.

Loris Giovannini

[1] L. Giovannini et al., Phys. Rev. B 70, 172404 (2004).

[2] M. Grimsditch, L. Giovannini, F. Montoncello, F. Nizzoli, G. K. Leaf,
and H. G. Kaper, Phys. Rev. B 70, 054409 (2004).


-------------------------------------------------------------------------------

COMPILATION

gfortran

Prerequisite: LAPACK library (http://www.netlib.org/lapack/). 

The main program (contained in micro3d.f) and all the required 
subroutines/functions can be compiled with the script "compile-gfortran".
Linking is done with the script "link-gfortran".

=================================

ifort (Intel fortran compiler, not tested)

Prerequisite: MKL (Intel) library, modification of inquiry.c (remove the
underscore from routine names called from fortran).  

The main program (contained in micro3d.f) and all the required 
subroutines/functions can be compiled with the script "compile-ifor".  
Linking is done with the script "link-ifor".  
Optimization and auto-parallelization should allow to obtain a faster 
code with ifort.

-------------------------------------------------------------------------------

USAGE

First of all, a file describing the shape, size and composition of the
magnetic particle whose spin modes are to be calculated must be prepared
and provided as input to the "micro3d" program. This is done by using
OOMMF (https://math.nist.gov/oommf/) to find the ground state, usually
with the Oxsii modulus. In turn, this requires the preparation of a
problem-description file (.mif) according to the rules of OOMMF. Important
note: DO NOT NORMALIZE THE MAGNETIZATION! In the .mif file use the option
"normalize_aveM_output 0".
An example is provided in the "test" directory.

=================================

Next, a file describing the material parameters must be prepared. Its
format is (line by line): 
1: comment 
2: description of the first material parameters, as: <material name> <Saturation magnetization, Oe> <Gyromagnetic ratio, 1/(s Oe)> <Anisotropy coefficient k2p, erg/cm^3> <Hard axis direction w.r. to the x-axis, deg> <Hard axis direction w.r. to the z-axis, deg> 
3: description of the second material parameters, as above.  
4: etc. for all materials. The list may include an arbitrary number of material and must at least contain the materials of the problem under investigation.  
n: comment 
n+1: exchange constant between materials 1 and 2 (if they coincides this is the intramaterial exchange constant of that material, otherwise it is and intermaterial exchange coupling) <material name 1> <material name 2> <Exchange stiffness parameter, erg/cm> 
n+2: second exchange constant, as above 
n+3: etc. for all desired material couples. The list may include an arbitrary number of materials and must at least contain the intramaterial exchange for all the materials of the problem under investigation and the intermaterial exchange for all possible combinations.

=================================

Dynamical matrix command usage:

"micro3d inputfile"

where inputfile is the name of the main input file. Its (fixed) format is
(line by line):

 1: comment
 2: comment
 3: comment
 4: <output filename> <.T. or .F. for binary format of the output file>
 5: comment
 6: comment
 7: <0: workspace query; 1: full calculation (approx. normalization); 2: full calculation (exact normalization)> <verbosity level 0...3 (0 is the lowest)>
 8: comment
 9: comment
10: <External magnetic field, Oe> <External field direction w.r. to the x-axis, deg> <External field direction w.r. to the z-axis, deg>
11: comment
12: comment 
13: <omf filename with the description of the magnetic particle and its ground state, see above> <file describing the material parameters, see above> 
14: comment 
15: comment 
16: <Lowest frequency of modes to be calculated, GHz> <Highest frequency of modes to be calculated, GHz> 
17: comment 
18: comment 
19: <Are we calculating a periodic repetition of nanopartiles? 0: no; 1: yes, repeated along x; 2: yes, repeated along y; 3: yes, bidimensional along x and y> <Components of the Bloch wavevector w.r. to the primitive reciprocal vectors, dimensionless >-0.5, <=0.5 (first BZ)>

=================================

The main input file, the .omf file with the ground state and the material
parameter file must be consistent!

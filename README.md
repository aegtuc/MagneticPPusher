Program MagneticPPusher is based on the Boris integration algorithm of relativistic charged particle equation of motion (Boris, 1970; Ripperda et al, 2018) particularized for the case of zero electric field.
  
Three possible magnetic field shapes are considered:
a)	Uniform throughout the domain, with components that can be specified
b)	Helical field in a cylindrical domain of predetermined radius, modulated by Bessel functions (cylindrical magnetic flux rope). Outside the cylinder radius the field is uniform whose components can be specified.
c)	Same as case b) but with a parabolic shock wave (sheath region).
In each of these cases a Gaussian distribution noise level can be specified.
The code is divided into 3 files: 
-	Main.f90
-	constantsModule.f90
-	subsfuncsModule.f90
An elemental batch file is included to compile on Windows with the GNU/gfortran compiler. Although it has also been tested with Intel/Fortran.  
A typical data input file “model.pshr” is also included.
Once compiled, the program can be run from the location of the input file. The results file will be written to the same location.
The result file is a csv file with 7 columns, time, particle possition (x,y,z) and magnetic field at those locations (Bx,By,Bz).
In this version, the format of the input file is not flexible in the sense that the order of the parameters as well as the number of lines between them must be respected. Dummy lines were used as comments describing what each parameter represents.
 
Any questions or comments, please do not hesitate to contact me: eguennam@herrera.unt.edu.ar

References:
Boris, J. P. (1970). Relativistic plasma simulation-optimization of a hybrid code. Fourth conf. numerical simulations of plasmas (Vol. 3).
Ripperda, B., Bacchini, F., Teunissen, J., Xia, C., Porth, O., Sironi, L., Lapenta, G., Keppens, R. (2018, mar). A comprehensive comparison of relativistic particle integrators. The Astrophysical Journal Supplement Series, 235 (1), 21.

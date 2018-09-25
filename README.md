# CmbRecTool

This Library contains Fortran90 subroutines to reconstruct lensing potential, cosmic birefrimgence, and patchy reionization 
from cosmic microwave background anisotropies (CMB). The Library also includes subrotuines for delensing, bi-spectrum calculation, and so on.

# Install

  0) The code assumes "ifort" (intel fortran compiler) as a fortran compiler. 

  1) Go to "pub" directory, and install each pubclic package (FFTW, Healpix, LAPACK and cfitsio). The clibrary utilizes static links to those packages.

  2) Go to the top directory, and type "./MALEALL.sh install". 
  
  3) You will find the Library at "lib" and "mod" directories. 

  4) To use some example codes, you also need "pylib/makefile.py" to create a Makefile if no Makefle exists. 

# Test

You can test the library by using the example codes at "examples" directory. At each subdirectory, you need to compile each code 
by just typying "make" (you do need to edit each Makefile). Then you can run each code by "./exe" or "./exe params.ini" where exe is 
an executable file. 


# References

In flat sky analysis, verification of the source code is shown in 

  - https://arxiv.org/abs/1209.0091 for the temperature lensing reconstruction, 
  - https://arxiv.org/abs/1310.2372 for the polarization lensing reconstruction, 
  - https://arxiv.org/abs/1612.07855 for the cosmic polarization rotation, 
  - https://arxiv.org/abs/1703.00169 for the delensing

For full sky analysis, the library also supports the subroutines for the lensing reconstruction and delensing as shown in 

  - https://arxiv.org/abs/1405.6568

The algorithm written in nldd/xxx.f90 is described in doc/note.pdf


# Acknowledgment

The library software uses the following public codes: FFTW, HEALPix, LAPACK, CFITSIO and also part of CAMB subroutines. 

# Contact

  namikawa at slac.stanford.edu


# CmbRecTool

This Library contains Fortran90 subroutines to reconstruct lensing potential and cosmic polarization rotation 
from cosmic microwave background anisotropies (CMB). 

# Install

  0) The code assumes "ifort" (intel fortran compiler) as a fortran compiler. 

  1) Go to "pub" directory, and install each pubclic package (FFTW, Healpix, LAPACK and cfitsio). 

  2) Go to the top directory, and type "./MALEALL install". 
  
  3) You will find the Library at "lib" and "mod" directories. 


# References

In flat sky analysis, verification of the source code is shown in 

  - https://arxiv.org/abs/1209.0091 for the temperature lensing reconstruction, 
  - https://arxiv.org/abs/1310.2372 for the polarization lensing reconstruction, 
  - https://arxiv.org/abs/1612.07855 for the cosmic polarization rotation, 
  - https://arxiv.org/abs/1703.00169 for the delensing

In full sky analysis, the library also support the subroutines for the lensing reconstruction and delensing as shown in 

  - https://arxiv.org/abs/1405.6568
  

The algorithm written in nldd/xxx.f90 is described in doc/note.pdf


# Acknowledgment

The library software uses the following public codes: FFTW, HEALPix, LAPACK, CFITSIO and also part of CAMB subroutines. 

# Contact

  toshiyan at stanford.edu


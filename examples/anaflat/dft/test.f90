!///////////////////////////////////////////////////////////////////////!
! * DFT test
! * Toshiya Namikawa
! - Last Modified: Tue 16 Aug 2016 06:20:04 PM PDT
!///////////////////////////////////////////////////////////////////////!

program test_dft
  use myutils, only: savetxt, linspace, savearray
  use myconst, only: dlc, pi
  use anaflat, only: elarray_inv
  use myfftw
  implicit none
  integer :: i, el(2), nn(2), npix
  double precision :: D(2)
  complex(dlc), allocatable, dimension(:) :: alm

  !* variables
  nn   = [236,100]
  npix = nn(1)*nn(2)
  D    = [59d0,25d0]*pi/180d0

  allocate(alm(npix))
  write(*,*) 'TEST_dft: Delta(x=y=0)=', dble(nn(1)*nn(2))/D(1)/D(2)
  write(*,*) 'TEST_dft: Delta(lx=ly=0)=', D(1)*D(2)

!  1) alm = 1 -> map = Delta(x,y)
  !* inverse dft
  alm = 1d0
  call dft_1darray(alm,nn,D,-1)
  call savearray('dft_alm.dat',alm,nn)

!  2) map = 1 -> alm = Delta(lx,ly)
  !* dft
  alm = 1d0
  call dft_1darray(alm,nn,D,1)
  call savearray('dft_map.dat',alm,nn)

  call dft_1darray(alm,nn,D,-1)
  call savearray('dft_unity.dat',alm,nn)

! 3) red spectrum
  alm = cmplx(elarray_inv(nn,D))

  !* inverse dft
  call dft_1darray(alm,nn,D,-1)
  call savearray('dft_redalm.dat',alm,nn)

  deallocate(alm)

end program test_dft


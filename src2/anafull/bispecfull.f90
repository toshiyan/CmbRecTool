!////////////////////////////////////////////////////!
! * Bispectrum in fullsky
!////////////////////////////////////////////////////!

module bispec_full
  use alm_tools, only: alm2map, map2alm
  use myconst,   only: dlc, pi

  private alm2map
  private dlc, pi

contains 


subroutine make_quad_gauss(nside,eL,alm)
  implicit none
  !I/O
  integer, intent(in) :: nside, eL(2)
  complex(dlc), intent(inout), dimension(1,0:eL(2),0:eL(2)) :: alm
  !internal
  double precision, allocatable :: kmap(:)

  allocate(kmap(0:12*nside**2-1))
  call alm2map(nside,eL(2),eL(2),alm,kmap)
  kmap = kmap + kmap**2
  call map2alm(nside,eL(2),eL(2),kmap,alm)
 deallocate(kmap)

end subroutine make_quad_gauss


subroutine bispec_equi(eL,bst,alm,bispec)
  implicit none
  !I/O
  integer, intent(in) :: eL(2), bst
  complex(dlc), intent(in), dimension(0:eL(2),0:eL(2)) :: alm
  double precision, intent(out) :: bispec
  !internal
  integer :: l, nside
  double precision, allocatable :: kmap(:)
  complex(dlc), allocatable :: klm(:,:,:)

  nside = eL(2)*bst

  allocate(kmap(0:12*nside**2-1),klm(1,0:eL(2),0:eL(2)))

  klm = 0d0
  do l = eL(1), eL(2) !ell filtering
    klm(1,l,0:l) = alm(l,0:l)
  end do

  call alm2map(nside,eL(2),eL(2),klm,kmap)
  bispec = sum(kmap**3)*(4d0*pi)/dble(12*nside**2)

  deallocate(kmap,klm)

end subroutine bispec_equi


subroutine bispec_fold(eL,bst,alm,bispec)
  implicit none
  !I/O
  integer, intent(in) :: eL(2), bst
  complex(dlc), intent(in), dimension(0:eL(2),0:eL(2)) :: alm
  double precision, intent(out) :: bispec
  !internal
  integer :: l, nside
  double precision, allocatable :: kmap(:,:)
  complex(dlc), allocatable :: klm(:,:,:)

  nside = eL(2)*bst

  allocate(kmap(0:12*nside**2-1,2),klm(2,0:eL(2),0:eL(2)))

  klm = 0d0
  do l = eL(1), eL(2) !ell filtering
    klm(1,l,0:l) = alm(l,0:l)
  end do
  do l = max(2,int(eL(1)/2d0)), int(eL(2)/2d0)
    klm(2,l,0:l) = alm(l,0:l)
  end do

  call alm2map(nside,eL(2),eL(2),klm(1:1,:,:),kmap(:,1))
  call alm2map(nside,eL(2),eL(2),klm(2:2,:,:),kmap(:,2))
  bispec = sum(kmap(:,1)*kmap(:,2)**2) * (4d0*pi)/dble(12*nside**2)

  deallocate(kmap,klm)

end subroutine bispec_fold


subroutine bispec_squi(eL,sL,bst,alm,bispec)
  implicit none
  !I/O
  integer, intent(in) :: eL(2), sL(2), bst
  complex(dlc), intent(in), dimension(0:eL(2),0:eL(2)) :: alm
  double precision, intent(out) :: bispec
  !internal
  integer :: l, nside
  double precision, allocatable :: kmap(:,:)
  complex(dlc), allocatable :: klm(:,:,:)

  nside = eL(2)*bst

  allocate(kmap(0:12*nside**2-1,2),klm(2,0:eL(2),0:eL(2)))

  klm = 0d0
  do l = sL(1), sL(2) !ell filtering
    klm(1,l,0:l) = alm(l,0:l)
  end do
  do l = eL(1), eL(2)
    klm(2,l,0:l) = alm(l,0:l)
  end do

  call alm2map(nside,eL(2),eL(2),klm(1:1,:,:),kmap(:,1))
  call alm2map(nside,eL(2),eL(2),klm(2:2,:,:),kmap(:,2))
  bispec = sum(kmap(:,1)*kmap(:,2)**2) * (4d0*pi)/dble(12*nside**2)

  deallocate(kmap,klm)

end subroutine bispec_squi


subroutine bispec_angl(eL,aL,l1,bst,alm,bispec)
  implicit none
  !I/O
  integer, intent(in) :: eL(2), aL(2), bst, l1
  complex(dlc), intent(in), dimension(0:l1,0:l1) :: alm
  double precision, intent(out) :: bispec
  !internal
  integer :: l, nside
  double precision, allocatable :: kmap(:,:)
  complex(dlc), allocatable :: klm(:,:,:)

  nside = l1*bst

  allocate(kmap(0:12*nside**2-1,2),klm(2,0:l1,0:l1))

  klm = 0d0
  do l = eL(1), eL(2) !ell filtering
    klm(1,l,0:l) = alm(l,0:l)
  end do
  do l = aL(1), aL(2)
    klm(2,l,0:l) = alm(l,0:l)
  end do

  call alm2map(nside,l1,l1,klm(1:1,:,:),kmap(:,1))
  call alm2map(nside,l1,l1,klm(2:2,:,:),kmap(:,2))
  bispec = sum(kmap(:,1)*kmap(:,2)**2) * (4d0*pi)/dble(12*nside**2)

  deallocate(kmap,klm)

end subroutine bispec_angl


end module bispec_full



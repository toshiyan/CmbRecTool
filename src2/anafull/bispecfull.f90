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
  double complex, intent(inout), dimension(1,0:eL(2),0:eL(2)) :: alm
  !internal
  double precision, allocatable :: kmap(:)

  allocate(kmap(0:12*nside**2-1))
  call alm2map(nside,eL(2),eL(2),alm,kmap)
  kmap = kmap + kmap**2
  call map2alm(nside,eL(2),eL(2),kmap,alm)
 deallocate(kmap)

end subroutine make_quad_gauss


subroutine bispec_equi_norm(bp0,bp1,bst,hl)
  implicit none
  integer, intent(in) :: bst
  double precision, intent(in) :: bp0(:), bp1(:)
  double precision, intent(out) :: hl(:)
  integer :: bmax, eL(2), lmax, l, b
  double complex, allocatable :: klm(:,:,:)

  bmax = size(hl)
  lmax = int(bp1(bmax))

  allocate(klm(1,0:lmax,0:lmax)); klm=0d0

  do l = 1, lmax
    klm(1,l,0) = dsqrt(2d0*l+1d0)
  end do

  do b = 1, bmax
    write(*,*) b
    eL = [int(bp0(b)),int(bp1(b))]
    call bispec_equi(eL,bst,klm(1,0:eL(2),0:eL(2)),hl(b))
    hl(b) = hl(b)/dsqrt(4d0*pi)
  end do

  deallocate(klm)

end subroutine bispec_equi_norm


subroutine bispec_equi_x(eL,bst,alm1,alm2,alm3,bispec)
  implicit none
  !I/O
  integer, intent(in) :: eL(2), bst
  double complex, intent(in), dimension(0:eL(2),0:eL(2)) :: alm1, alm2, alm3
  double precision, intent(out) :: bispec
  !internal
  integer :: l, nside
  double precision, allocatable :: kmap(:,:)
  double complex, allocatable :: klm(:,:,:)

  nside = eL(2)*bst

  allocate(kmap(3,0:12*nside**2-1),klm(3,0:eL(2),0:eL(2)))

  klm = 0d0
  do l = eL(1), eL(2) !ell filtering
    klm(1,l,0:l) = alm1(l,0:l)
    klm(2,l,0:l) = alm2(l,0:l)
    klm(3,l,0:l) = alm3(l,0:l)
  end do

  call alm2map(nside,eL(2),eL(2),klm(1:1,:,:),kmap(1,:))
  call alm2map(nside,eL(2),eL(2),klm(2:2,:,:),kmap(2,:))
  call alm2map(nside,eL(2),eL(2),klm(3:3,:,:),kmap(3,:))
  bispec = sum(kmap(1,:)*kmap(2,:)*kmap(3,:))*(4d0*pi)/(12d0*dble(nside)**2)

  deallocate(kmap,klm)

end subroutine bispec_equi_x


subroutine bispec_equi(eL,bst,alm,bispec)
  implicit none
  !I/O
  integer, intent(in) :: eL(2), bst
  double complex, intent(in), dimension(0:eL(2),0:eL(2)) :: alm
  double precision, intent(out) :: bispec
  !internal
  integer :: l, nside
  double precision, allocatable :: kmap(:)
  double complex, allocatable :: klm(:,:,:)

  nside = eL(2)*bst

  allocate(kmap(0:12*nside**2-1),klm(1,0:eL(2),0:eL(2)))

  klm = 0d0
  do l = eL(1), eL(2) !ell filtering
    klm(1,l,0:l) = alm(l,0:l)
  end do

  call alm2map(nside,eL(2),eL(2),klm,kmap)
  bispec = sum(kmap**3)*(4d0*pi)/(12d0*dble(nside)**2)

  deallocate(kmap,klm)

end subroutine bispec_equi


subroutine bispec_fold(eL,bst,alm,bispec)
  implicit none
  !I/O
  integer, intent(in) :: eL(2), bst
  double complex, intent(in), dimension(0:eL(2),0:eL(2)) :: alm
  double precision, intent(out) :: bispec
  !internal
  integer :: l, nside
  double precision, allocatable :: kmap(:,:)
  double complex, allocatable :: klm(:,:,:)

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
  bispec = sum(kmap(:,1)*kmap(:,2)**2) * (4d0*pi)/(12d0*dble(nside)**2)

  deallocate(kmap,klm)

end subroutine bispec_fold


subroutine bispec_sque(eL,sL,bst,alm,bispec)
  implicit none
  !I/O
  integer, intent(in) :: eL(2), sL(2), bst
  double complex, intent(in), dimension(0:eL(2),0:eL(2)) :: alm
  double precision, intent(out) :: bispec
  !internal
  integer :: l, nside, lmax
  double precision, allocatable :: kmap(:,:)
  double complex, allocatable :: klm(:,:,:)

  lmax = max(eL(2),sL(2))

  nside = lmax*bst

  allocate(kmap(0:12*nside**2-1,2),klm(2,0:lmax,0:lmax))

  klm = 0d0
  do l = sL(1), sL(2) !ell filtering
    klm(1,l,0:l) = alm(l,0:l)
  end do
  do l = eL(1), eL(2)
    klm(2,l,0:l) = alm(l,0:l)
  end do

  call alm2map(nside,lmax,lmax,klm(1:1,:,:),kmap(:,1))
  call alm2map(nside,lmax,lmax,klm(2:2,:,:),kmap(:,2))
  bispec = sum(kmap(:,1)*kmap(:,2)**2) * (4d0*pi)/(12d0*dble(nside)**2)

  deallocate(kmap,klm)

end subroutine bispec_sque


subroutine bispec_angl(eL,aL,l1,bst,alm,bispec)
  implicit none
  !I/O
  integer, intent(in) :: eL(2), aL(2), bst, l1
  double complex, intent(in), dimension(0:l1,0:l1) :: alm
  double precision, intent(out) :: bispec
  !internal
  integer :: l, nside
  double precision, allocatable :: kmap(:,:)
  double complex, allocatable :: klm(:,:,:)

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
  bispec = sum(kmap(:,1)*kmap(:,2)**2) * (4d0*pi)/(12d0*dble(nside)**2)

  deallocate(kmap,klm)

end subroutine bispec_angl


end module bispec_full



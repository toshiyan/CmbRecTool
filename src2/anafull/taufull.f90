!////////////////////////////////////////////////////!
! * Inhomogeneous Reion Reconstruction in Fullsky
!////////////////////////////////////////////////////!

module taufull
  use alm_tools, only: alm2map, map2alm
  use myconst, only: dl, dlc, iu

  private alm2map, map2alm
  private dl, dlc, iu

contains 


!* Tlm,Elm,Blm: inverse-variance filtered alms (Cov^-1 (Tlm,Elm,Blm)^t)

subroutine quadtt_tau(nside,Tlm,fC,eL,rL,xlm)
!* fC = ClTT
  implicit none
  !I/O
  integer :: eL(2), rL(2), nside
  real(dl), intent(in) :: fC(:)
  complex(dlc), intent(in), dimension(0:eL(2),0:eL(2)) :: Tlm
  complex(dlc), intent(inout), dimension(0:rL(2),0:rL(2)) :: xlm
  !internal
  integer :: l, npix
  double precision, allocatable :: map(:,:)
  complex(dlc), allocatable :: alm(:,:,:)

  write(*,*) 'calc TT-estimator (tau)'
  npix = 12*nside**2

  !* alm to map 
  allocate(alm(2,0:rL(2),0:rL(2))); alm = 0d0
  do l = rL(1), rL(2)
    alm(1,l,0:l) = Tlm(l,0:l)
    alm(2,l,0:l) = fC(l)*Tlm(l,0:l)
  end do 
  allocate(map(2,0:npix-1))
  call alm2map(nside,rL(2),rL(2),alm(1:1,:,:),map(1,:))
  call alm2map(nside,rL(2),rL(2),alm(2:2,:,:),map(2,:))
  deallocate(alm)

  !* map to alm
  allocate(alm(1,0:rL(2),0:rL(2)))
  call map2alm(nside,rL(2),rL(2),map(1,:)*map(2,:),alm)
  do l = rL(1), rL(2)
    xlm(l,0:l) = alm(1,l,0:l)
  end do
  deallocate(map,alm)

end subroutine quadtt_tau


end module taufull



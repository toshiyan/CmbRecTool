!////////////////////////////////////////////////////!
! * Delensing in Fullsky
!////////////////////////////////////////////////////!

module delensfull
  use alm_tools, only: alm2map, alm2map_spin, map2alm_spin
  use myconst, only: dl, dlc, iu

  private alm2map, alm2map_spin, map2alm_spin
  private dl, dlc, iu

contains 


subroutine LensingB(nside,WElm,Wglm,lBlm,eL,tL)
  implicit none
  !I/O
  integer, intent(in) :: eL(2), tL(2), nside
  complex(dlc), intent(in), dimension(0:tL(2),0:tL(2)) :: WElm, Wglm
  complex(dlc), intent(out), dimension(0:eL(2),0:eL(2)) :: lBlm
  !internal
  integer :: l, m, npix
  real(dl), dimension(:,:), allocatable :: A1,A3,A,map
  complex(dlc), dimension(:,:,:), allocatable :: alm

  npix = 12*nside**2
  allocate(alm(2,0:tL(2),0:tL(2)))

  write(*,*) 'Elm to map (spin-1 transform)'
  alm = 0d0
  do l = tL(1), tL(2)
    alm(1,l,:) = WElm(l,:)*dsqrt(dble((l+2)*(l-1))*0.5)
  end do 
  allocate(A1(0:npix-1,2))
  call alm2map_spin(nside,tL(2),tL(2),1,alm,A1)

  write(*,*) 'Elm to map (spin-3 transform)'
  allocate(A3(0:npix-1,2))
  alm = 0d0
  do l = tL(1), tL(2)
    alm(1,l,:) = WElm(l,:)*dsqrt(dble((l-2)*(l+3))*0.5)
  end do 
  call alm2map_spin(nside,tL(2),tL(2),3,alm,A3)

  write(*,*) 'glm to map (spin-1 transform)'
  allocate(A(0:npix-1,2))
  alm = 0d0
  do l = tL(1), tL(2)
    alm(1,l,:) = Wglm(l,:)*dsqrt(dble(l*(l+1))*0.5)
  end do 
  call alm2map_spin(nside,tL(2),tL(2),1,alm,A)
  deallocate(alm)

  !convolution
  allocate(map(0:npix-1,2))
  !map = A1*A - A3*conjg(A)
  !    = (rA1+iu*iA1)*(rA+iu*iA) - (rA3+iu*iA3)*(rA-iu*iA)
  !    = rA*(rA1-rA3+iu*(iA1-iA3)) + iA*(iu*rA1-iA1+iu*rA3-iA3)
  map(:,1) = A(:,1)*(A1(:,1)-A3(:,1)) - A(:,2)*(A1(:,2)+A3(:,2))
  map(:,2) = A(:,1)*(A1(:,2)-A3(:,2)) + A(:,2)*(A1(:,1)+A3(:,1))
  deallocate(A1,A3,A)

  allocate(alm(2,0:eL(2),0:eL(2)))
  call map2alm_spin(nside,eL(2),eL(2),2,map,alm)
  deallocate(map)

  lBlm = alm(2,:,:)
  deallocate(alm)

end subroutine LensingB


end module delensfull



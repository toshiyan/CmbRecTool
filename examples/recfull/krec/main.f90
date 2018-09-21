! example of fullsky lensing reconstruction

program main
  use readfile,  only: set_params_file, read_prm, read_str, read_int, read_dbl
  use myutils,   only: str, linspace, savetxt, loadtxt, save_average
  use io,        only: loadfile
  use myconst,   only: dlc, pi
  use mycls,     only: alm2bcl, binned_ells
  use nldd_lens, only: AlTT, AlEB
  use anafull,   only: gaussianTEB
  use recfull,   only: quadtt, quadeb
  implicit none
  integer :: i, l, bn, nside, rL(2), eL(2), sn(2)
  double precision, allocatable :: bc(:), ll(:), Fl(:,:,:), cl(:,:), cpp(:,:,:), Al(:,:,:)
  complex(dlc), allocatable :: alm(:,:,:), glm(:,:,:), clm(:,:,:)

  ! read parameters from params.ini
  call set_params_file
  nside = read_int('nside')
  bn    = read_int('bn')
  call read_prm('sn',sn)
  call read_prm('eL',eL)
  call read_prm('rL',rL)

  ! multipole
  allocate(ll(eL(2)),bc(bn))
  ll = linspace(1,eL(2))
  call binned_ells(eL,bc=bc) !binned multipole

  ! read CMB cl (change here to read your cl)
  allocate(cl(4,eL(2))); cl=0d0
  call loadtxt('../../dat/lensedfid_P15.dat',cl(1,2:),cl(2,2:),cl(3,2:),cl(4,2:),rows=[1,eL(2)-1],usecols=[2,3,4,5])
  cl(1,:) = cl(1,:) * 2d0*pi/(ll**2+ll) ! factor out

  ! set filtering function
  allocate(Fl(3,0:eL(2),0:eL(2))); Fl=0d0
  do l = 2, eL(2)
    Fl(1,l,0:l) = 1d0/cl(1,l) !noise not included
  end do

  ! compute normalization
  allocate(Al(2,2,eL(2)))
  call AlTT(rL,eL,Al(1,1,:),Al(1,2,:),cl(1,:),cl(1,:)) !noise not included
  !call AlEB(rL,eL,Al(2,1,:),Al(2,2,:),cl(2,:),cl(2,:),cl(3,:))
  !call savetxt('Al.dat',ll,Al(1,:))

  ! calculate estimator
  allocate(alm(3,0:eL(2),0:eL(2)),glm(2,0:eL(2),0:eL(2)),clm(2,0:eL(2),0:eL(2)),cpp(sn(1):sn(2),2,bn))

  do i = sn(1), sn(2) ! iterate for realizations

    !//// you need some input alm and read it here ////!
    ! 1) read alm
    !call loadfile('alm_r'//str(i)//'.dat',alm) !change this line to read your data appropriately
    ! 2) or generate random gaussian alms
    call gaussianTEB(alm,cl(1,:),cl(2,:),cl(3,:),cl(4,:),eL(2))

    ! inverse variance filtering
    alm = alm*Fl

    ! compute unnormalized quadratic estimator
    call quadtt(nside,alm(1,0:rL(2),0:rL(2)),alm(1,0:rL(2),0:rL(2)),cl(1,:),eL,rL,glm(1,:,:),clm(1,:,:))
    !call quadeb(nside,alm(2,0:rL(2),0:rL(2)),alm(3,0:rL(2),0:rL(2)),cl(2,:),eL,rL,glm(2,:,:),clm(2,:,:))

    ! correct normalization
    do l = 2, eL(2)
      glm(1,l,0:l) = glm(1,l,0:l)*Al(1,1,l)
      clm(1,l,0:l) = clm(1,l,0:l)*Al(1,2,l)
    end do

    ! compute binned power of the phi estimator
    call alm2bcl(bn,eL,glm(1,:,:),cb=cpp(i,1,:),oL=eL)

    ! save to file
    call savetxt('cpp_r'//str(i)//'.dat',bc,cpp(i,1,:),ow=.true.)

  end do

  ! output average of cl over realizations
  call save_average('cpp.dat',cpp,id=[2,3],bc=bc)

  deallocate(alm,glm,clm,Fl,cl,ll,cpp,Al,bc)

end program main


!///////////////////////////////////////////////////////////////////////!
! * CMB Lensing Bispectrum
!///////////////////////////////////////////////////////////////////////!

program main
  use readfile, only: set_params_file, read_prm, read_str, read_int, read_dbl, read_val, read_log
  use myutils,  only: GLdxs, GLpoints, loadtxt, linspace, savetxt, str, gauss_legendre_params, gl_initialize, gl_finalize
  use myfunc
  use mycls,    only: binned_ells
  use cmblbisp
  implicit none
  character(128) :: ltype, fpk, fcl, cpmodel, model
  character(4) :: btypes(4)
  logical :: nonlinear
  integer :: l, l0, zn, kn, oL(2), i, b, bn, eL(2), hL(2), sL(2), aL(2)
  double precision :: h, om, zmax, zmin, lan, kan, del, snr, ns, s0, ncur
  double precision, allocatable :: nef(:), dnq(:), PE(:,:), Ik(:,:), vbl(:,:), D(:), knl(:), bp(:), bc(:), z(:), dz(:), cl(:), nldd(:,:), fac(:), datl(:,:), pl(:,:,:), bl(:,:), wp(:,:), k(:,:), ck(:,:), abc(:,:,:), lambda(:), kappa(:)
  type(gauss_legendre_params) :: GL
  type(cosmoparams) :: cp

  call set_params_file

  ! cosmological parameters (to compute background quantities)
  cpmodel = read_str('cpmodel')
  select case(cpmodel)
  case ('model0')
    cp%Om = 0.31201550908d0
    cp%H0 = 67.51d0
    cp%w0 = -1d0
    cp%wa = 0d0
    cp%nu = 0.06d0/(93.14d0*(cp%H0/100d0)**2*cp%Om)
    kn   = 1380
  case ('modelw')
    cp%Om = 0.279d0
    cp%H0 = 70d0
    cp%w0 = -1d0
    cp%wa = 0d0
    cp%nu = 0d0
    kn   = 1380
  case ('modelp')
    cp%Om = 0.3156d0
    cp%H0 = 67.27d0
    cp%w0 = -1d0
    cp%wa = 0d0
    cp%nu = 0.00443d0
    kn   = 1610
  end select
  h     = cp%H0/100d0
  cp%Ov = 1d0 - cp%Om
  ns    = 0.9645d0

  ! LSS matter bispec fitting model
  model = read_str('model')

  ! MG extension
  lan   = read_dbl('lan')
  kan   = read_dbl('kan')

  ! other parameters
  zn   = read_int('zn')
  call read_prm('oL',oL) !multipole of theoretical bisp
  if (oL(1)<1) stop 'oLmin should be >=1' 
  zmax = read_dbl('zmax')
  zmin = read_dbl('zmin')

  ! filename
  fpk = '/global/homes/t/toshiyan/Work/Ongoing/bisp/data/'//trim(cpmodel)//'/Pk/Pklin.dat'
  fcl = '/global/homes/t/toshiyan/Work/Ongoing/bisp/data/'//trim(cpmodel)//'/cl/fid.dat'

  ! binning
  bn = read_int('bn')
  allocate(bp(bn+1),bc(bn))
  call binned_ells(oL,bp,bc)

  ! load linear Pk at z=0
  allocate(datl(zn+1,kn))
  call loadtxt(fpk,datl(1:2,:),fsize=[2,kn])
  datl(1,:) = datl(1,:)*h    ! k
  datl(2,:) = datl(2,:)/h**3 ! linear Pk at z=0

  ! precomputing interpolation points for z
  allocate(z(zn),dz(zn))
  call gl_initialize(gl,zn,1d-15)
  select case (read_int('zspace')) 
  case(0)
    z  = linspace(zmin,zmax,zn)
    dz = z(2)-z(1)
  case(1)
    z  = glpoints([zmin,zmax],gl%z) ! points
    dz = gldxs([zmin,zmax],gl%w)    ! width
  end select

  ! precompute quantities for bispectrum
  allocate(D(zn),pl(2,zn,oL(2)),k(zn,oL(2)),fac(zn),abc(3,zn,oL(2)),wp(zn,oL(2)),ck(zn,oL(2)),knl(zn));  Pl=0d0; k=0d0; fac=0d0; abc=1d0; wp=0d0; ck=0d0
  call prep_lens_bispectrum(z,dz,read_dbl('zs'),cp,datl(1,:),datl(2,:),model,k,pl,fac,abc,wp,ck,knl=knl,s0=s0,ftype=read_str('nlpk'))
  D = D_z(z,cp)  !growth factor

  allocate(nef(zn),dnq(zn),PE(zn,oL(2)),Ik(zn,oL(2))); nef=1d0; dnq=1d0; PE=1d0; Ik=1d0
  if (model=='RT') then
    do i = 1, zn
      if (z(i)<=9d0) then
        om = omega_m(1d0/(1d0+z(i)),cp) 
        call find_pknl_params(datl(1,:),D(i)**2*datl(2,:),knl(i),nef(i),ncur,nonlinear)
        call RTformula_3h_funcs(k(i,:)/knl(i),h,om,s0*D(i),nef(i),pl(1,i,:),Ik(i,:),PE(i,:),dnq(i))
        !if (i==6)  call savetxt('PE.dat',k(i,:),PE(i,:),Ik(i,:),pl(1,i,:),ow=.true.)
      end if
    end do
  end if

  deallocate(datl)

  ! MG parameter z-evolution
  allocate(lambda(zn),kappa(zn));  lambda=1d0; kappa=1d0
  do i = 1, zn
    if (lan/=0d0) lambda(i) = (omega_m(1d0/(1d0+z(i)),cp))**lan
    if (kan/=0d0) kappa(i)  = (omega_m(1d0/(1d0+z(i)),cp))**kan
  end do
  
  ! bispectrum (density and post-Bron bispectrum)
  write(*,*) 'compute bispectrum'
  allocate(bl(9,oL(2)));  bl=0d0
  ltype = ''
  bl(1,:) = linspace(1,oL(2))
  l0 = int(bc(1))
  if (mod(l0,2)/=0) l0 = l0+1
  btypes = ['equi','fold','sque','angl']
  do i = 1, 4
    write(*,*) btypes(i)
    call bispec_lens(btypes(i),oL,k,pl,fac,abc,wp,ck,bl(i+1,:),l0,h,ns,s0,z,D,knl,nef,dnq,PE,Ik,model,'lss',ltype,lambda,kappa)
    call bispec_lens(btypes(i),oL,k,pl,fac,abc,wp,ck,bl(i+5,:),l0,h,ns,s0,z,D,knl,nef,dnq,PE,Ik,model,'pbn',ltype,lambda,kappa)
  end do

  ! save bispectrum
  call savetxt('bl_'//trim(model)//'.dat',bl(:,3:),ow=.true.)

  !* flat-binned bispectrum
  allocate(vbl(9,bn)); vbl=0d0
  vbl(1,:) = bc
  do b = 1, bn
    eL = int(bp(b:b+1))
    hL = [max(2,int(bp(b)/2)),int(bp(b+1)/2)]
    sL = int(bp(1:2))
    !sL = int(bp(2:3))
    !sL = int(bp(3:4))
    !sL = int(bp(4:5))
    aL = int(bp(bn/2:bn/2+1))
    write(*,*) eL, hL, sL, aL
    call bispec_lens_bin(eL,eL,eL,k,Pl,fac,abc,wp,ck,h,ns,s0,z,D,knl,nef,dnq,PE,Ik,model,'lss',vbl(2,b))
    call bispec_lens_bin(eL,hL,hL,k,Pl,fac,abc,wp,ck,h,ns,s0,z,D,knl,nef,dnq,PE,Ik,model,'lss',vbl(3,b))
    call bispec_lens_bin(sL,eL,eL,k,Pl,fac,abc,wp,ck,h,ns,s0,z,D,knl,nef,dnq,PE,Ik,model,'lss',vbl(4,b))
    call bispec_lens_bin(eL,aL,aL,k,Pl,fac,abc,wp,ck,h,ns,s0,z,D,knl,nef,dnq,PE,Ik,model,'lss',vbl(5,b))
    call bispec_lens_bin(eL,eL,eL,k,Pl,fac,abc,wp,ck,h,ns,s0,z,D,knl,nef,dnq,PE,Ik,model,'pbn',vbl(6,b))
    call bispec_lens_bin(eL,hL,hL,k,Pl,fac,abc,wp,ck,h,ns,s0,z,D,knl,nef,dnq,PE,Ik,model,'pbn',vbl(7,b))
    call bispec_lens_bin(sL,eL,eL,k,Pl,fac,abc,wp,ck,h,ns,s0,z,D,knl,nef,dnq,PE,Ik,model,'pbn',vbl(8,b))
    call bispec_lens_bin(eL,aL,aL,k,Pl,fac,abc,wp,ck,h,ns,s0,z,D,knl,nef,dnq,PE,Ik,model,'pbn',vbl(9,b))
    write(*,*) vbl(:,b)
  end do
  call savetxt('bl_'//trim(model)//'_b'//str(bn)//'.dat',vbl,ow=.true.)
  deallocate(vbl)

end program main


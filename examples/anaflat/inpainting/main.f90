!////////////////////////////////////////////////////!
! Map Making with Inpainting
!////////////////////////////////////////////////////!

program main
  use readFile
  use myconst
  use mycls
  use mycmbexp
  use anaflat
  use myfftw
  implicit none
  character(LEN=50) :: fname
  integer :: simn,el(2),nu,n,npix,nn(2),mm(2),mpix,nc,elt(2)
  real(dl) :: s, mr, DD(2), Wn(4)
  real(dl), allocatable :: Nl(:)
  complex(dlc), dimension(:), allocatable :: P,T,EMAP,BMAP,SimP,Wf
  complex(dlc), dimension(:), allocatable :: readP,readT,nTmap,nPmap
  complex(dlc), dimension(:,:,:), allocatable :: vec

  call SET_PARAMS_FILE

  !realizations
  simn = Read_Int('simnum')
  CALL READ_PARAMS('lrange',el)  
  CALL READ_PARAMS('ltrange',elt)

  !map and mask specification
  CALL READ_PARAMS('nside',nn)
  CALL READ_PARAMS('mside',mm)
  W%sky_frac = Read_Double('sky_frac')
  mr = 0.d0
  if(.not.Read_Int('mask_num')==0) mr = Read_Double('mask_size')

  !#### Initial Settings ####
  s = READ_DOUBLE("pixelsize")*pi/180.d0
  DD = nn*s
  npix = nn(1)*nn(2)
  mpix = mm(1)*mm(2)

  call SET_TYPE_OBS(O)
  allocate(O%LC(O%cln,elt(2)))
  call READCL(O,O%LC,READ_STRING('LCl_infile'),el)

  !Noise map
  nc = READ_INT('nchan')
  allocate(nTmap(npix),nPmap(npix),Nl(elt(2)))
  do nu = 1, nc
    Nl = READ_DOUBLE("sigmaT")*ARCMIN_TO_RAD/Tcmb
    call GAUSSIAN_ALM(nn,DD,elt,nTmap,NL)
    Nl = READ_DOUBLE("sigmaP")*ARCMIN_TO_RAD/Tcmb
    call GAUSSIAN_ALM(nn,DD,elt,nPmap,NL)
  end do
  deallocate(Nl)

  !Mask map
  allocate(Wf(npix))
  call WINDOW_GENERATE(nn,DD,Wf,mr,READ_INT("mask_num"))
  call WINDOW_NORM(nn,1,Wf,Wn(1))

  do n = 1, simn
    allocate(T(npix),P(npix))
    if(Read_Logical('MapRead')) then 
      call SELECT_MAPFILE_RT(n,fname,Read_String("maproot"))
      allocate(readT(mpix),readP(mpix))
      call MAP_READ_CMPLX(mm,fname,readT,readP)
      call MAP_CUT_CMPLX(readT,mm,T,nn)
      call MAP_CUT_CMPLX(readP,mm,P,nn)
      deallocate(readT,readP)
    else 
      CALL GAUSSIAN_MAP(nn,mm,DD,O%LC(O%TT,:),T)
      CALL GAUSSIAN_MAP_POL(nn,mm,DD,O%LC(O%EE,:),O%LC(O%BB,:),P)
    end if
    !add instrumental noise
    T = T + nTmap
    P = P + nPmap
    !construct data vector
    allocate(vec(3,nc,npix))
    do nu = 1, nc 
      T = T*Wf
      P = P*Wf
      call FFTrans(T,nn,DD,1,0)
      call FFTrans(P,nn,DD,1,2)
      vec(1,nu,:) = T
      vec(2,nu,:) = (P+conjg(P))*0.5d0
      vec(3,nu,:) = (P-conjg(P))*0.5d0/iu
    end do
    deallocate(T,P)
    !! need to restore !!
    !call inpaint_interface(DATA,Wf,pmap,Read_Logical('Cinv'))
    deallocate(vec)
  end do

  deallocate(Wf,nTmap,nPmap,O%LC)

end program main


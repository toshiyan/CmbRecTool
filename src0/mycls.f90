!////////////////////////////////////////////////////!
! * Angular power spectrum calculation
!////////////////////////////////////////////////////!

module mycls
  use myconst, only: dlc, pi, TT, TE, EE, BB, dd, Td, Ed, oo
  use myutils, only: splint, spline, FileColumns, FileLines, savetxt, linspace, interp_lin
  implicit none

  interface calccl_flat
    module procedure calccl_flat_alm1d, calccl_flat_alm2d
  end interface calccl_flat

  private dlc, pi, TT, TE, EE, BB, dd, Td, Ed, oo
  private splint, spline, FileColumns, FileLines, savetxt, linspace, interp_lin

contains


subroutine binned_ells(eL,bp,bc,spc)
! * return binned multipole edges and centers
  implicit none
  !I/O
  character(*), intent(in), optional :: spc
  integer, intent(in) :: eL(:)
  double precision, intent(out), optional :: bp(:), bc(:)
  !internal
  character(8) :: spcs=''
  integer :: i, n
  double precision :: d
  double precision, allocatable :: p(:)

  if (present(bc))  n = size(bc)
  if (present(bp))  n = size(bp) - 1
  if (present(spc)) spcs = spc

  allocate(p(n+1))
  p(1) = eL(1)

  select case(spcs)
  case('log')
    d     = dlog(dble(eL(2))/dble(eL(1)))/dble(n)
    p(2:) = [((eL(1)*exp(d*(i-1))),i=2,n+1)]
  case default
    d     = dble(eL(2)-eL(1))/dble(n)
    p(2:) = [((d*(i-1)+p(1)),i=2,n+1)]
  end select
  if (present(bc)) bc = [(((p(i+1)+p(i))*0.5d0),i=1,n)]
  if (present(bp)) bp = p
 
  deallocate(p)

end subroutine binned_ells


subroutine optimal_weight(Cb,Vb,mCb,vCb,f,bc)
! optimally-weighted binned cls
  implicit none
  !I/O
  character(*), intent(in), optional :: f
  double precision, intent(in) :: Cb(:,:)
  double precision, intent(in) :: Vb(:,:)
  double precision, intent(in), optional :: bc(:)
  double precision, intent(out), optional :: mCb(:), vCb(:)
  !internal
  integer :: b, i, bmax
  double precision, allocatable :: W(:), Cl(:)

  bmax = size(Cb,dim=2)
  allocate(W(bmax),Cl(bmax))

  !* optimal weighting
  W = 0d0;  Cl = 0d0
  do b = 1, bmax !loop for bins
    do i = 1, size(Cb,dim=1)
      if(Vb(i,b)<0d0) stop 'error (optimal_weight): variance is negative'
    end do
    if(product(Vb(:,b))>0d0) then
      W(b) = 1d0/sum(1d0/Vb(:,b)**2) !square of variance
      Cl(b) = sum(Cb(:,b)/(Vb(:,b)**2))*W(b)
    end if
  end do

  !* output
  if(present(mCb)) mCb = Cl
  if(present(vCb)) vCb = dsqrt(W)
  if(present(f).and.present(bc))       call savetxt(f,bc,Cl,dsqrt(W))
  if(present(f).and..not.present(bc))  call savetxt(f,linspace(1,bmax,bmax),Cl,dsqrt(W))

  deallocate(W,Cl)

end subroutine optimal_weight


subroutine cl2cb(bc,Cl,Cb,f)
! * compute a simple binned Cl from non-binned Cl
  implicit none
  !I/O
  character(*), intent(in), optional :: f
  double precision, intent(in) :: Cl(:), bc(:)
  double precision, intent(out) :: Cb(:)
  !internal
  integer :: i

  Cb = 0d0
  do i = 1, size(bc)
    Cb(i) = Cl(int(bc(i)))
  end do

  if (present(f)) call savetxt(f,Cb)

end subroutine cl2cb


function cl2cbf(bc,Cl) result(Cb)
! * compute a simple binned Cl from non-binned Cl
  implicit none
  !I/O
  double precision, intent(in) :: Cl(:), bc(:)
  double precision :: Cb(size(bc))
  !internal
  integer :: i

  Cb = 0d0
  do i = 1, size(bc)
    Cb(i) = Cl(int(bc(i)))
  end do

end function cl2cbf


subroutine check_positive(dat,p)
  implicit none
  double precision, intent(in) :: dat(:)
  logical, intent(out) :: p
  integer :: i

  p = .true.
  do i = 1, size(dat)
    if (dat(i)<0d0) p = .false.
  end do

end subroutine check_positive


subroutine cb2cl(bc,Cb,Cl,f,bp,method)
! interpolate non-binned Cl from binned Cl
  implicit none
  !I/O
  double precision, intent(in) :: bc(:), Cb(:)
  double precision, intent(out) :: Cl(:)
  !(optional)
  character(*), intent(in), optional :: f, method
  double precision, intent(in), optional :: bp(:)
  !internal
  logical :: p
  character(8) :: m
  integer :: ln, bn, i, j
  double precision :: y2(size(Cb)), Cb0(size(Cb)+1)

  m  = ''
  bn = size(Cb)
  Cl = 0d0
  ln = size(Cl)
  if (present(method)) m = method

  call check_positive(Cb,p) !check positivity

  !interpolation
  select case (m)
  case ('step')
    if (.not.present(bp)) stop '(cb2cl): bin edge is required'
    do i = 1, bn
      if (ln<int(bp(i))) goto 1
      do j = int(bp(i)), min(ln,int(bp(i+1)))
        Cl(j) = Cb(i)
      end do
    end do
  case ('linear')
    Cl(1:int(bc(1))) = Cb(1)
    do i = 1, bn-1
      if (ln<int(bc(i))) goto 1
      do j = int(bc(i)), min(ln,int(bc(i+1)))
        Cl(j) = interp_lin(dble(j),bc(i),bc(i+1),Cb(i),Cb(i+1))
      end do
    end do
    Cl(int(bc(bn))+1:) = Cb(bn)
  case default
    call spline(bc,Cb,bn,0d0,0d0,y2)
    if (ln<int(bc(1))) goto 1
    do j = int(bc(1)), min(ln,int(bc(bn)))
      Cl(j) = splint(dble(j),bc,Cb,y2)
    end do
1 end select

  !check cl
  do j = 1, size(Cl)
    if (Cl(j)>=0d0.or..not.p) cycle
    call savetxt('mycls_cb2cl.out',Cl,ow=.true.)
    write(*,*) 'cb2cl: interpolated power is negative', j, Cl(j)
    Cl(j) = 0d0
  end do

  !save if you want
  if (present(f)) call savetxt(f,Cl)

end subroutine cb2cl


subroutine cb2c2d(bc,els,eL,Cb,C2d,f,method0,method1,bp)
! Interpolate binned Cl to 2D Cl
  implicit none
  !I/O
  integer, intent(in) :: eL(2)
  double precision, intent(in) :: bc(:), els(:), Cb(:)
  double precision, intent(out) :: C2d(:)
  character(*), intent(in), optional :: f, method0, method1
  double precision, intent(in), optional :: bp(:)
  !internal
  character(8) :: m(2)
  double precision, allocatable :: Cl(:), bp0(:)

  allocate(Cl(eL(2)),bp0(size(bc)+1)); Cl=0d0; bp0=0d0

  m = ['','']
  if(present(bp))      bp0  = bp
  if(present(method0)) m(1) = method0 

  !interpolate Cb -> Cl
  call cb2cl(bc,Cb,Cl,bp=bp0,method=m(1))

  !* output
  if(present(f)) call savetxt(f,linspace(1,eL(2)),Cl,ow=.true.)

  !interpolate Cl -> C2d
  call cl2c2d(els,Cl,eL,c2d)

  deallocate(Cl)

end subroutine cb2c2d


subroutine cl2c2d(els,Cl,eL,c2d)
!* Transform Cl to Cl2D with linear interpolation
  implicit none
  !I/O
  integer, intent(in) :: eL(:)
  double precision, intent(in) :: els(:), Cl(:)
  double precision, intent(out) :: c2d(:)
  !internal
  logical :: p
  integer :: n, l0, l1

  !check positivity
  call check_positive(Cl,p)

  c2d = 0d0
  do n = 1, size(els)
    if(eL(1)>els(n).or.els(n)>eL(2)-1) cycle
    l0 = int(els(n))
    l1 = l0 + 1
    c2d(n) = Cl(l0) + (els(n)-l0)*(Cl(l1)-Cl(l0))
    if (c2d(n)>=0d0.or..not.p) cycle
    write(*,*) Cl(l0), Cl(l1), l0, els(n)
    stop 'error (cl2c2d): interpolated Cl is negative'
  end do

end subroutine cl2c2d


subroutine cl_interp_spline(bls,bCl,tL,Cl,islog)
  implicit none
  !I/O
  logical, intent(in) :: islog 
  integer, intent(in) :: tL(:)
  double precision, intent(in) :: bls(:), bCl(:)
  double precision, intent(out) :: Cl(:)
  !internal
  integer :: l, n, i, nmax
  double precision, allocatable :: tbCl(:), tbls(:), tbCldd(:)

  !choose non-zero bCls
  n = 0
  do i = 1, size(bCl)
    if(bCl(i)>0) n = n + 1
  end do
  nmax = n

  allocate(tbCl(nmax),tbls(nmax),tbCldd(nmax))
  n = 0
  do i = 1, size(bCl)
    if(bCl(i)<=0) cycle
    n = n + 1
    tbCl(n) = bCl(i)
    tbls(n) = bls(i)
  end do

  Cl = 0d0
  call spline(tbls,tbCl,nmax,0d0,0d0,tbCldd)
  do l = tL(1), tL(2)
    if(islog) then 
      Cl(l) = splint(dlog(dble(l)),tbls,tbCl,tbCldd)
    else 
      Cl(l) = splint(dble(l),tbls,tbCl,tbCldd)
    end if
  end do

  deallocate(tbCl,tbls)

end subroutine cl_interp_spline


function alm2cl2d(D,alm1,alm2) result(C)
!* compute power spectrum in 2D Fourier grids: C_ell
  implicit none
  !I/O
  double precision, intent(in) :: D(1:2)
  complex(dlc), intent(in) :: alm1(:)
  complex(dlc), intent(in), optional :: alm2(:)
  double precision :: C(size(alm1))

  ! delta(l=0) = D(1)*D(2)
  if(present(alm2)) then
    C = (alm1*conjg(alm2)+alm2*conjg(alm1))*0.5d0/(D(1)*D(2)) 
  else
    C = alm1*conjg(alm1)/(D(1)*D(2)) 
  end if

end function alm2cl2d


subroutine alm2bcl_flat(bmax,oL,els,D,alm1,alm2,Cb,Vb,norm,f,spc)
!* computing binned cl from alm(s). 
  implicit none
  integer, intent(in) :: bmax, oL(1:2)
  double precision, intent(in) :: D(1:2), els(:)
  complex(dlc), intent(in) :: alm1(:)
!(optional)
  character(*), intent(in), optional :: f, spc
  double precision, intent(in), optional :: norm
  double precision, intent(out), optional :: Cb(:), Vb(:)
  complex(dlc), intent(in), optional :: alm2(:)
  !internal
  character(8) :: spc0, f0
  integer :: npix
  double precision :: norm0
  double precision, allocatable :: C(:), Cb0(:), Vb0(:)
  complex(dlc), allocatable :: alm0(:)

  spc0  = ''
  f0    = ''
  norm0 = 1d0
  npix  = size(alm1)
  if(present(norm))  norm0=norm
  if(present(f))     f0=f
  if(present(spc))   spc0=spc

  !* 2D power spectrum
  allocate(C(npix),alm0(npix))
  alm0 = alm1
  if(present(alm2)) alm0=alm2
  C = (alm1*conjg(alm0)+alm0*conjg(alm1))*0.5d0*norm0/(D(1)*D(2)) 
  deallocate(alm0)

  !* to 1D power spectrum
  allocate(Cb0(bmax),Vb0(bmax))
  call calcbcl_flat(bmax,oL,els,C,Cb0,Vb0,f=f0,spc=spc0)

  if (present(Cb))  Cb = Cb0
  if (present(Vb))  Vb = Vb0

  deallocate(C,Cb0,Vb0)

end subroutine alm2bcl_flat


subroutine calcbcl_flat(bmax,oL,els,Cl,Cb,Vb,tCl,tVl,f,spc)
! * calculate binned cls 
! - [Note]
! - If input alm is zero at oL(1)<=l<=oL(2), the power should be smaller than expectation.
  implicit none
  !I/O
  character(*), intent(in), optional :: f, spc
  integer, intent(in) :: oL(2), bmax
  double precision, intent(in) :: els(:), Cl(:)
  double precision, intent(in), optional :: tCl(:), tVl(:)
  double precision, intent(out), optional :: Cb(:), Vb(:)
  !internal
  character(8) :: spc0
  integer :: i, nmax
  double precision :: bp(bmax+1), b(bmax)
  double precision, allocatable, dimension(:) :: AL, BL, Ab, vAb, iCl, iVl

  spc0 = ''
  if (size(Cl)<1)  stop 'error (calcBcl_flat): size(Cl) is less than 1'
  if (present(Cb).and.size(Cb)/=bmax) stop 'error (calcBcl_flat): size(Cb) is not equal to bmax'
  if (present(spc)) spc0 = spc

  call binned_ells(oL,bp,b,spc0)

  nmax = size(els)
  if (nmax<1) stop 'error (calcBcl_flat): size(els) is less than 1'
  allocate(AL(nmax),BL(nmax),Ab(bmax),vAb(bmax),iCl(nmax),iVl(nmax))
  iCl = 1d0
  iVl = 1d0

  !* interpolate theoretical values at each els(i)
  if(present(tCl)) call cl2c2d(els,tCl,oL,iCl)
  if(present(tVl)) call cl2c2d(els,tVl,oL,iVl)

  !* transfer to the amplitude parameter and its variance
  AL = 0d0
  BL = 0d0
  do i = 1, nmax
    if(els(i)<oL(1).or.els(i)>oL(2)) cycle
    if(iCl(i)<=0d0) stop 'error (calcbcl_flat): iCl'
    if(iVl(i)<=0d0) stop 'error (calcbcl_flat): iVl'
    AL(i) = Cl(i)/iCl(i)
    BL(i) = 1d0
    if(present(tVl)) BL(i) = iCl(i)**2/iVl(i)**2
  end do
  call power_binning_flat(bp,els,AL,BL,Ab,vAb)
  deallocate(AL,BL,iCl,iVl)

  if(present(tCl)) then !transfer to Cls
    Ab  = Ab*cl2cbf(b,tCl)
    vAb = vAb*cl2cbf(b,tCl)
    !Ab = Ab*cl2cb(b,tCl,el)
    !vAb = vAb*cl2cb(b,tCl,el)
  end if

  if(present(Cb)) Cb = Ab
  if(present(Vb)) Vb = vAb
  if(present(f).and.f/='')  call savetxt(f,b,Ab,vAb)

  deallocate(Ab,vAb)

end subroutine calcbcl_flat


subroutine power_binning_flat(bp,els,AL,BL,bA,bV)
!* compute averaging factor: 
!   bA = bV^2 * \sum_{l=bp(i)}^{bp(i+1)} BL(l)*AL(l)
! AL: amplitude
! BL: band-pass filter
!
  implicit none
!
! [input]
!   bp  --- edges of multipole binning
!   els --- multipoles as a function of npix
!   AL  --- amplitude
!   BL  --- band-pass filter
  double precision, intent(in) :: bp(:), els(:), AL(:), BL(:)
!
! [output]
!   bA  --- binned amplitude
!   bV  --- binned variance
  double precision, intent(out) :: bA(:), bV(:)
!
! [internal]
  integer :: b, i, npix, bmax
  double precision, allocatable :: bW(:)

  npix = size(els)
  bmax = size(bp) - 1

  allocate(bW(bmax))
  !* bW = [ sum_i BL(i) ]^-1 = bV^2 : optimal variance at b-th bin
  bW = 0d0;  bA = 0d0;  bV = 0d0
  do b = 1, bmax
    do i = 1, npix
      if (els(i)<bp(b).or.els(i)>=bp(b+1)) cycle
      bW(b) = bW(b) + BL(i) 
      bA(b) = bA(b) + BL(i)*Al(i)
    end do
    if (bW(b) <= 0d0) cycle
    bW(b) = 1d0/bW(b) !square of variance
    bA(b) = bA(b)*bW(b)
    bV(b) = dsqrt(bW(b))
  end do

  deallocate(bW)

end subroutine power_binning_flat


subroutine power_binning_full(bp,eL,Al,vAl,bA,bV)
  implicit none 
  !I/O
  integer, intent(in) :: eL(2)
  double precision, intent(in) :: bp(:), Al(:), vAl(:)
  double precision, intent(out) :: bA(:), bV(:)
  !internal 
  integer :: b, i, bmax
  double precision, allocatable :: bW(:)

  bmax = size(bp) - 1

  allocate(bW(bmax))
  !* compute optimal variance of b-th bin, bW(b)
  bW = 0d0; bA = 0d0; bV = 0d0
  do b = 1, bmax
    do i = el(1), el(2)
      if(i<bp(b).or.i>=bp(b+1)) cycle
      bW(b) = bW(b) + vAl(i)
      bA(b) = bA(b) + vAl(i)*Al(i)
    end do
    if(bW(b)<=0d0) cycle
    bW(b) = 1d0/bW(b)
    bA(b) = bA(b)*BW(b)
    bV(b) = dsqrt(bW(b))
  end do

end subroutine power_binning_full


subroutine power_labeling(el,els,nmax,label)
!* find independent els in 2D grid
  implicit none 
  !I/O
  integer, intent(in) :: el(2)
  double precision, intent(in) :: els(:)
  integer, intent(out) :: nmax, label(:)
  !internal
  integer :: npix, m, n, indep

  npix = size(els)
  nmax = 0

  do n = 1, npix
    indep=-1
    label(n) = 0
    if(els(n)<el(1).or.els(n)>=el(2)) cycle
    do m = 1, n-1
      if(els(n)==els(m)) then
        indep = m
        go to 20
      end if
    end do
20  if(indep==-1) then 
      nmax = nmax + 1
      label(n) = nmax 
    else 
      label(n) = label(indep)
    end if
  end do

end subroutine power_labeling


subroutine cl_elcut(nn,D,Cl) 
!* assign zeros to Cls outsize of the largest annuli
  implicit none
  !I/O
  integer, intent(in) :: nn(1:2)
  double precision, intent(in) :: D(1:2)
  double precision, intent(inout) :: Cl(:)
  !internal
  integer :: i, j, n
  double precision :: x, y, l, lmax

  lmax = 2*pi * min( dble(nn(1)*0.5d0-1)/D(1), dble(nn(2)*0.5d0-1)/D(2))
  n = 1
  do i = 1, nn(1)
    x = 2*pi*dble(i-1-nn(1)*0.5d0)/D(1)
    do j = 1, nn(2)
      y = 2*pi*dble(j-1-nn(2)*0.5d0)/D(2)
      l = dsqrt(x**2+y**2)
      if(l>lmax) Cl(n) = 0d0
      n = n + 1
    end do
  end do

end subroutine cl_elcut


subroutine calccl(alm1,alm2,el,Cl,f)
  implicit none
  !I/O
  character(*), intent(in), optional :: f
  integer, intent(in) :: el(2)
  double precision, intent(out), optional :: Cl(:)
  complex(dlc), intent(in), dimension(0:el(2),0:el(2)) :: alm1, alm2
  !intenral
  integer :: l
  double precision :: tCl(el(2))

  tCl = 0d0
  do l = el(1), el(2)
    tCl(l) = ( real(alm1(l,0)*alm2(l,0)) + 2.*sum(alm1(l,1:l)*conjg(alm2(l,1:l))))/(2.*l+1.)
  end do
  if(present(f))  call savetxt(f,linspace(1,el(2)),tCl)
  if(present(Cl)) Cl=tCl

end subroutine calccl


subroutine calccl_flat_alm1d(D,alm1,alm2,Cl,els,f)
  !compute power spectrum from Fourier grids: C_ell
  implicit none
  !I/O
  character(*), intent(in), optional :: f
  double precision, intent(in) :: D(2)
  double precision, intent(in), optional :: els(:)
  double precision, intent(out), optional :: Cl(:)
  complex(dlc), intent(in) :: alm1(:), alm2(:)
  !internal
  integer :: i, npix
  double precision, allocatable :: C(:)

  npix = size(alm1)
  allocate(C(npix))
  C = (alm1*conjg(alm2)+alm2*conjg(alm1))*0.5d0/(D(1)*D(2))  ! divided by delta(l=0)
  if(present(Cl)) Cl = C
  if(present(f).and.present(els))  call savetxt(f,els,C)
  deallocate(C)

end subroutine calccl_flat_alm1d


subroutine calccl_flat_alm2d(D,alm1,alm2,Cl,els,f)
  !compute power spectrum from Fourier grids: C_ell
  implicit none
  !I/O
  character(*), intent(in), optional :: f
  double precision, intent(in) :: D(2)
  double precision, intent(in), optional :: els(:)
  double precision, intent(out), optional :: Cl(:)
  complex(dlc), intent(in) :: alm1(:,:), alm2(:,:)
  !internal
  integer :: i, j, n, npix
  double precision, allocatable :: C(:)

  npix = size(alm1)
  allocate(C(npix))
  n = 1
  do i = 1, size(alm1,dim=1)
    do j = 1, size(alm1,dim=2)
      C(n) = (alm1(i,j)*conjg(alm2(i,j))+alm2(i,j)*conjg(alm1(i,j)))*0.5d0/(D(1)*D(2))  ! divided by delta(l=0)
      n = n + 1
    end do
  end do
  if(present(Cl)) Cl = C
  if(present(f).and.present(els))  call savetxt(f,els,C)
  deallocate(C)

end subroutine calccl_flat_alm2d


subroutine angle_average(nn,alm,Al,nmax,label,el)
  implicit none
  !subroutine for computing power spectrum from Fourier grids
  !I/O
  integer, intent(in) :: nn(2), nmax, label(:)
  double precision, intent(in) :: el(:)
  double precision, intent(out) :: Al(:)
  complex(dlc), intent(in) :: alm(:)
  !internal
  integer :: i, m, n
  double precision :: num(nmax)

  !generate [ l(i), Al(i), Cl(i) ] table 
  num = 0
  Al = 0.d0
  do i = 1, nmax 
    do n = 1, nn(1)*nn(2)
      if(label(n)==i) then
        num(i) = num(i) + 1
        Al(i) = Al(i) + alm(n) 
      end if
      Al(i) = Al(i)/dble(num(i))
    end do
  end do

end subroutine angle_average


subroutine readcl_camb(Cl,f,eL,HasBB,rowCl)
!* Read cls from a ascii file
!* This routine assumes that the file is obtained by CAMB and contains 
!*   l, TT, EE, TE, dd, Td, (Ed) 
!* if hasbb = false (e.g., scalcls.dat), or 
!*   l, TT, EE, BB, TE, (curl)
!* if hasbb = true (e.g., lensedcls.dat). 
!* For other cases, loadtxt is more convenient. 
!* cls are assumed to be multipiled by l*(l+1)/2pi unless you specify rowCl
  implicit none 
!
![input]
! f : filename
! eL : 
  character(*), intent(in) :: f
  integer, intent(in) :: eL(1:2)
  double precision, intent(out) :: Cl(:,:)
!
!(optional)
! HasBB : file has BB spectrum at 4th column
! rowCl : do not multiply, 2pi/l(l+1), to cls
  logical, intent(in), optional :: HasBB, rowCl
!
![internal]
  double precision :: rC(1:10)
  integer :: i, l, n, m

  open(unit=20,file=trim(f),status='old')
  n = FileColumns(20)
  m = FileLines(20)
  do i = 1, m
    read(20,*) l, rC(1:n-1)
    if (.not.(present(rowCl).and.rowCl))  rC = rC*2d0*pi/dble(l**2+l)
    if (el(1)>l.or.l>el(2)) cycle
    cl(TT,l) = rC(1)
    cl(EE,l) = rC(2)
    if(present(HasBB).and.HasBB) then
      cl(BB,l) = rC(3)
      cl(TE,l) = rC(4)
      if(n>=6) cl(oo,l) = rC(5)
    else
      cl(TE,l) = rC(3)
      cl(dd,l) = rC(4)
      cl(Td,l) = rC(5)
      if(n>=7) cl(Ed,l) = rC(6)
    end if
  end do
  close(20)

end subroutine readcl_camb


subroutine map_vars(sigma,eL,cl)
! * variance of a map and its derivative
  implicit none
  !I/O
  integer, intent(in) :: eL(2)
  double precision, intent(in) :: cl(:)
  double precision, intent(out) :: sigma(0:1)
  !internal
  integer :: l
  double precision :: al
  
  sigma = 0d0
  do l = el(1), el(2)
    al = dble(l)
    sigma(0) = sigma(0) + (2d0*al+1d0)*cl(l)
    sigma(1) = sigma(1) + (2d0*al+1d0)*(al**2+al)*cl(l)
  end do 
  sigma = sigma/(4d0*pi)

end subroutine map_vars


end module mycls


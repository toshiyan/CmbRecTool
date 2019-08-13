!////////////////////////////////////////////////////////////////!
! * CMB lensing bispectrum subroutines for analytic calculation
!////////////////////////////////////////////////////////////////!

module cmblbisp
  use myconst, only: pi
  use myutils, only: neighb, spline, splint, savetxt, gauss_legendre_params, gl_initialize, gl_finalize, glpoints, gldxs, linspace
  use myfunc,  only: C_z, H_z, D_z, NonLinRatios, pk2sigma, cosmoparams, omega_m, find_pknl_params
  implicit none

  type bispecfunc
    !z
    double precision, allocatable :: z(:)
    !basic functions
    double precision, allocatable :: knl(:), n(:), sz(:), zker(:), D(:)
    !interpolated k and Pk^lin, Pk^NL at k=ell/chi
    double precision, allocatable :: kl(:,:), q(:,:), plL(:,:), pl(:,:) 
    !fitting formula coefficients in F2 kernel
    double precision, allocatable :: abc(:,:,:) 
    !MG extension parameters
    double precision, allocatable :: mgp(:,:) 
    !RT
    double precision, allocatable :: dnq(:), Ik(:,:), PE(:,:)
  endtype bispecfunc

  private pi
  private neighb, spline, splint, savetxt, gauss_legendre_params, gl_initialize, gl_finalize, glpoints, gldxs, linspace
  private C_z, H_z, D_z, NonLinRatios, pk2sigma, cosmoparams, omega_m, find_pknl_params

contains


subroutine bispec_lens_lss_init(cp,b,z,dz,zs,k,pkl0,oL,model,ftype,btype)
  implicit none
  !I/O
  type(cosmoparams), intent(in) :: cp
  type(bispecfunc), intent(inout) :: b
  character(*), intent(in) :: model, ftype
  integer, intent(in) :: oL(2)
  double precision, intent(in)  :: z(:), dz(:), zs, k(:), pkl0(:)
  character(*), intent(in), optional :: btype
  !internal
  character(4) :: bisptype = 'kkk'
  logical :: nonlinear
  integer :: i, zn, kn
  double precision :: chis, h, s0, om, ncur
  double precision, allocatable :: chi(:), pk(:,:), pkL(:,:)

  zn = size(z)  !number of z points for z-integral
  kn = size(k)  !number of CAMB output k

  !* lensing weights
  if (present(btype)) bisptype = btype
  allocate(b%zker(zn),b%z(zn))
  b%z = z
  call lens_zweight(cp,z,dz,zs,b%zker,bisptype)

  !* get linear and non-linear Pk at each z
  allocate(pkL(zn,kn),pk(zn,kn))
  call get_pk(cp,z,k,pkL0,pkL,pk,ftype)

  !* sigma_8
  allocate(b%sz(zn),b%D(zn))
  s0   = dsqrt(pk2sigma(8d0/cp%h,k,pkL0))
  b%D  = D_z(z,cp)
  b%sz = s0*b%D
  write(*,*) 'sigma8 = ', s0

  !* interpolate k, Pk at k=l/chi
  allocate(chi(zn),b%kl(zn,oL(2)),b%plL(zn,oL(2)),b%pl(zn,oL(2)))
  chi = C_z(z,cp)  !comoving distance at each z
  if (2d0/chi(zn)<k(1)) write(*,*) 'warning: required minimum k is smaller than input', 2d0/chi(zn), k(1)
  call Limber_k2l(chi,k,pkL,b%kl,b%plL)
  call Limber_k2l(chi,k,pk,b%kl,b%pl)

  !precompute knl
  allocate(b%knl(zn),b%q(zn,oL(2)))
  call get_knl(k,pkl,b%knl)
  do i = 1, zn
    b%q(i,:) = b%kl(i,:)/b%knl(i)
  end do

  !* precompute F2-kernel coefficients a, b, c 
  allocate(b%abc(3,zn,oL(2))); b%abc = 1d0
  if (model=='SC'.or.model=='GM')  call coeff_abc(b,k,pkL,model)

  !* RT
  if (model=='RT') then
    allocate(b%dnq(zn),b%PE(zn,oL(2)),b%Ik(zn,oL(2)),b%n(zn)); b%dnq=1d0; b%PE=1d0; b%Ik=1d0; b%n=1d0
    do i = 1, zn
      if (z(i)<=9d0) then
        om = omega_m(1d0/(1d0+z(i)),cp)
        call find_pknl_params(k,pkL(i,:),b%knl(i),b%n(i),ncur,nonlinear)
        call RTformula_3h_funcs(b%q(i,:),cp%h,om,b%sz(i),b%n(i),b%plL(i,:),b%Ik(i,:),b%PE(i,:),b%dnq(i))
      end if
    end do
  end if

  deallocate(b%knl,b%n,b%dnq,b%Ik,b%PE)

end subroutine bispec_lens_lss_init


subroutine bispec_lens_lss(cp,b,shap,eL,model,ltype,l0,bl)
! lensing bispectrum
  implicit none
  !I/O
  type(cosmoparams), intent(in) :: cp
  type(bispecfunc), intent(in) :: b
  character(*), intent(in) :: shap, model, ltype
  integer, intent(in) :: eL(2), l0
  double precision, intent(out) :: bl(:)
  !internal
  integer :: l, zn, l1, l2, l3

  ! initial set up
  bl = 0d0
  zn = size(b%kl,dim=1)

  ! compute bispectrum
  do l = eL(1), eL(2)

    if (shap=='fold'.and.mod(l,2)/=0) cycle
    if (shap=='sque'.and.2*l<l0)      cycle

    call set_three_ells(shap,l,eL,l0,l1,l2,l3)
    call bispec_lens_lss_kernel(cp,b,l1,l2,l3,bl(l),model)

    if (ltype=='full') bl(l) = bl(l) * W3j_approx(dble(l1),dble(l2),dble(l3)) * dsqrt((2d0*l1+1d0)*(2d0*l2+1d0)*(2d0*l3+1d0)/(4d0*pi))

  end do

end subroutine bispec_lens_lss


subroutine bispec_lens_lss_kernel(cp,b,l1,l2,l3,bl,model)
  implicit none
  !input
  type(cosmoparams), intent(in) :: cp
  type(bispecfunc), intent(in) :: b
  character(*), intent(in) :: model
  integer, intent(in) :: l1, l2, l3
  !output
  double precision, intent(out) :: bl
  !internal
  integer :: i, zn
  double precision :: bk, fh(3)

  zn = size(b%z)

  bl = 0d0
  do i = 1, zn

    write(*,*) i
    select case (model)
    case('SC','GM','TR')
      call bispec_matter_9p(b%q(i,[l1,l2,l3]),b%pl(i,[l1,l2,l3]),reshape(b%abc(:,i,[l1,l2,l3]),[3,3]),bk,b%mgp(:,i))
    case('3B')
      fh = coeff_fih(b%kl(i,l1)+b%kl(i,l2)+b%kl(i,l3),b%D(i),cp%h,b%knl(i))
      call bispec_matter_3B(b%q(i,[l1,l2,l3]),b%plL(i,[l1,l2,l3]),b%pl(i,[l1,l2,l3]),fh,bk)
    case('RT')
      write(*,*) '!'
      fh = b%dnq(i)*b%q(i,[l1,l2,l3])
      write(*,*) '!'
      call bispec_matter_RT(cp,b%q(i,[l1,l2,l3]),b%plL(i,[l1,l2,l3]),b%z(i),b%sz(i),b%n(i),b%Ik(i,l1)*b%Ik(i,l2)*b%Ik(i,l3),b%PE(i,[l1,l2,l3]),fh,bk)
    end select

    bl = bl + b%zker(i) * bk

  end do

end subroutine bispec_lens_lss_kernel


subroutine bispec_matter_9p(q,pl,abc,bk,mgparams)
  implicit none
  !input
  double precision, intent(in) :: q(3), pl(3), abc(3,3)
  double precision, intent(in), optional :: mgparams(2)
  !output
  double precision, intent(out) :: bk
  !internal
  double precision :: F2(3), mgp(2)

  !MG extension
  mgp = 1d0
  if (present(mgparams))  mgp = mgparams

  call F2_Kernel(q([1,2,3]),abc(:,1),abc(:,2),F2(1),mgp)
  call F2_Kernel(q([2,3,1]),abc(:,2),abc(:,3),F2(2),mgp)
  call F2_Kernel(q([3,1,2]),abc(:,3),abc(:,1),F2(3),mgp)
  bk = 2d0*(F2(1)*pl(1)*pl(2) + F2(2)*pl(2)*pl(3) + F2(3)*pl(3)*pl(1))

end subroutine bispec_matter_9p


subroutine bispec_matter_3B(q,plL,pl,fh,bk)
  implicit none
  !input
  double precision, intent(in) :: q(3), plL(3), pl(3), fh(3)
  !output
  double precision, intent(out) :: bk
  !internal
  double precision :: F2(3), p(3)

  p(1) = 1.0011993d0 !Om^{-1/143} ~ 1.008
  p(2) = 1d0
  p(3) = 0.9969955d0
  call F2_Kernel(q([1,2,3]),p,p,F2(1))
  call F2_Kernel(q([2,3,1]),p,p,F2(2))
  call F2_Kernel(q([3,1,2]),p,p,F2(3))
  bk = fh(1) + fh(2)*(plL(1)*plL(2)+plL(2)*plL(3)+plL(3)*plL(1))/3d0 + fh(3)*2d0*(F2(1)*pl(1)*pl(2)+F2(2)*pl(2)*pl(3)+F2(3)*pl(3)*pl(1))

end subroutine bispec_matter_3B


subroutine bispec_matter_RT(cp,q,plL,z,sz,nef,Iks,PE,fh,bk)
  implicit none
  !input
  type(cosmoparams), intent(in) :: cp
  double precision, intent(in) :: q(3), plL(3), z, sz, nef, Iks, PE(3), fh(3)
  !output
  double precision, intent(out) :: bk
  !internal
  double precision :: F2(3), p(3)

  p(1) = 1d0
  p(2) = 1d0
  p(3) = 1d0

  call F2_Kernel(q([1,2,3]),p,p,F2(1))
  call F2_Kernel(q([2,3,1]),p,p,F2(2))
  call F2_Kernel(q([3,1,2]),p,p,F2(3))

  if (z<9d0) then
    call RTformula_1h(q,cp%ns,sz,nef,bk)
    bk = bk/cp%h**6 + 2d0*Iks*((F2(1)+fh(3))*PE(1)*PE(2) + (F2(2)+fh(1))*PE(2)*PE(3) + (F2(3)+fh(2))*PE(3)*PE(1))
  else
    bk = 2d0*(F2(1)*plL(1)*plL(2) + F2(2)*plL(2)*plL(3) + F2(3)*plL(3)*plL(1))
  end if

end subroutine bispec_matter_RT


subroutine F2_Kernel(k,abc1,abc2,F2,mgparams)
! F2 kernel
  implicit none
  double precision, intent(in) :: k(3), abc1(3), abc2(3)
  double precision, intent(out) :: F2
  double precision, intent(in), optional :: mgparams(2)
  !internal
  double precision :: ct, lam, kap

  !MG extension
  lam = 1d0
  kap = 1d0
  if (present(mgparams))  lam = mgparams(1)
  if (present(mgparams))  kap = mgparams(2)

  !Kernel
  ct = (k(3)**2-k(1)**2-k(2)**2)/(2d0*k(1)*k(2)) !cos(theta) of vectors k1 and k2
  F2 = (kap-lam*2d0/7d0)*abc1(1)*abc2(1) + kap*(k(1)**2+k(2)**2)/(2d0*k(1)*k(2))*ct*abc1(2)*abc2(2) + lam*2d0/7d0*ct**2*abc1(3)*abc2(3)

end subroutine F2_Kernel


subroutine coeff_abc(b,k,pkL,model)
! precomputing F2-kernel fitting coefficitents
  implicit none
  !I/O
  type(bispecfunc), intent(inout) :: b
  character(*), intent(in) :: model
  double precision, intent(in) :: k(:), pkL(:,:)
  !internal
  integer :: zn, ln, i, l
  double precision :: Q3f, qan, qbn, qcn, q35, q30, q, n2, a(9)
  double precision, allocatable :: n(:,:)

  !select a model
  select case (model)
  case('SC') 
    a = [0.25d0,3.5d0,2d0,1d0,2d0,-0.2d0,1d0,0d0,0d0]
  case('GM')  
    a = [0.484d0,3.74d0,-0.849d0,0.392d0,1.013d0,-0.575d0,0.128d0,-0.722d0,-0.926d0]
  case default
    stop 'error (precompute_coeff_abc): need to specify model of the nonlinear correction'
  end select

  zn = size(b%kl,dim=1)
  ln = size(b%kl,dim=2)

  allocate(n(zn,ln))
  call get_neff(k,b%kl,pkL,n,model)

  !compute a, b and c
  do i = 1, zn
    do l = 1, ln
      q   = b%kl(i,l)/b%knl(i)
      n2  = 2d0**n(i,l)
      qan = (q*a(1))**(n(i,l)+a(2))
      qbn = (q*a(7))**(n(i,l)+3d0+a(8))
      qcn = (q*a(5))**(n(i,l)+3d0+a(9))
      Q3f = dsqrt(0.7d0*((4d0-n2)/(1d0+2d0*n2)))*b%sz(i)**a(6)
      b%abc(1,i,l) = (1d0 + Q3f*qan )/(1d0+qan)
      b%abc(2,i,l) = (1d0 + 0.2d0*a(3)*(n(i,l)+3d0)*qbn )/(1d0+qbn*dsqrt(q*a(7)))
      b%abc(3,i,l) = (1d0 + (4.5d0*a(4)/(1.5d0+(n(i,l)+3d0)**4d0))*qcn)/(1d0+qcn*dsqrt(q*a(5)))
    end do
  end do

  deallocate(n)

end subroutine coeff_abc


function coeff_fih(Kl,Dz,h,knl)  result(f)
  implicit none
  double precision, intent(in) :: knl, Dz, h, Kl
  double precision :: f(3)

  f(1) = (2.45d6*Dz**8/(0.8d0+0.2d0/Dz**3)) / (1d0+0.054*Dz**2.2d0*Kl**2/h**2)**2
  f(2) = 140*Dz**(-5d0/4d0)/(1d0+1.9d0*Dz**(-1.5d0)*h/Kl)**3
  f(3) = dexp(-Kl/(7.5d0*knl))

end function coeff_fih


subroutine RTformula_1h(q,ns,sigma8,n_eff,bk)
  implicit none
  double precision, intent(in)  :: q(3), ns, sigma8, n_eff
  double precision, intent(out) :: bk
  integer :: i, imin(1), imax(1)
  double precision :: kmin, kmax, kmid, r1, r2
  double precision :: gamman, an, bn, cn, alphan, betan

  !set kmin, kmax and r1, r2
  kmin  = minval(q)
  kmax  = maxval(q)
  if (kmin==kmax) then
    kmid = kmin
  else
    imax = minloc(q)
    imin = maxloc(q)
    do i = 1, 3
      if (i/=imin(1).and.i/=imax(1)) kmid = q(i)
    end do
  end if
  r1 = kmin/kmax
  r2 = (kmid+kmin-kmax)/kmax

  !coefficients
  gamman = 10**(0.182+0.57*n_eff)
  an = -2.167-2.944*log10(sigma8)-1.106*log10(sigma8)**2-2.865*log10(sigma8)**3-0.310*r1**gamman
  bn = -3.428-2.681*log10(sigma8)+1.624*log10(sigma8)**2-0.095*log10(sigma8)**3
  cn = 0.159-1.107*n_eff
  alphan = -4.348 - 3.006*n_eff - 0.5745*n_eff**2 + 10.**(-0.9+0.2*n_eff)*r2**2
  betan  = -1.731 - 2.845*n_eff - 1.4995*n_eff**2 - 0.2811*n_eff**3 + 0.007*r2

  an = 10**an
  bn = 10**bn
  cn = 10**cn
  alphan = 10**alphan
  betan  = 10**betan
  
  if (alphan>1.-(2./3.)*ns)  alphan = 1.-(2./3.)*ns

  bk = 1d0
  do i = 1, 3
    bk = bk * 1d0/(an*q(i)**alphan+bn*q(i)**betan) * 1d0/(1d0+1d0/(cn*q(i)))
  end do

end subroutine RTformula_1h


subroutine RTformula_3h_funcs(q,h,om,sigma8,n_eff,Plin,Ik,PE,dnq)
  implicit none
  double precision, intent(in)  :: q(:), h, om, sigma8, n_eff, Plin(:)
  double precision, intent(out) :: PE(:), IK(:), dnq
  integer :: ki
  double precision :: fn, gn, hn, mn, nn, mun, nun, pn, dn, en

  fn  = -10.533 - 16.838*n_eff - 9.3048*n_eff**2 - 1.8263*n_eff**3
  gn  = 2.787  + 2.405*n_eff + 0.4577*n_eff**2
  hn  = -1.118 - 0.394*n_eff
  mn  = -2.605 - 2.434*log10(sigma8) + 5.71*log10(sigma8)**2
  nn  = -4.468 - 3.08*log10(sigma8) + 1.035*log10(sigma8)**2
  mun = 15.312 + 22.977*n_eff + 10.9579*n_eff**2 + 1.6586*n_eff**3
  nun = 1.347  + 1.246*n_eff  + 0.4525*n_eff**2
  pn  = 0.071  - 0.433*n_eff
  dn  = -0.483 + 0.892*log10(sigma8) - 0.086*om
  en  = -0.632 + 0.646*n_eff

  fn  = 10d0**fn
  gn  = 10d0**gn
  hn  = 10d0**hn
  mn  = 10d0**mn
  nn  = 10d0**nn
  mun = 10d0**mun
  nun = 10d0**nun
  pn  = 10d0**pn
  dn  = 10d0**dn
  en  = 10d0**en

  do ki = 1, size(PE)
    PE(ki) = (1d0+fn*q(ki)**2)/(1d0+gn*q(ki)+hn*q(ki)**2)*Plin(ki) + (1d0/h**3)/(mn*q(ki)**mun+nn*q(ki)**nun) * 1d0/(1d0+1d0/(pn*q(ki))**3)
    Ik(ki) = 1d0/(1d0+en*q(ki))
  end do

  dnq = dn

end subroutine RTformula_3h_funcs


!//// Post-Born Contributions ////!

subroutine bispec_lens_pb_init(cp,kl,pl,z,dz,zs,oL,wp,ck)
  implicit none
  type(cosmoparams), intent(in) :: cp
  integer, intent(in) :: oL(2)
  double precision, intent(in) :: kl(:,:), pl(:,:), z(:), dz(:), zs
  double precision, intent(out) :: wp(:,:), ck(:,:)
  integer :: zn
  double precision :: chis
  double precision, allocatable :: chi(:), Hz(:), wlf(:)

  zn = size(z)

  ! precompute quantities for post-Born bispectrum
  allocate(chi(zn),Hz(zn),wlf(zn))
  chis = C_z(zs,cp) !source comoving distance
  chi  = C_z(z,cp)  !comoving distance at each z
  Hz   = H_z(z,cp)  !expansion rate at z
  wlf  = 1.5d0*cp%Om*(cp%H0/3d5)**2*(1d0+z) !matter -> potential conversion factor (matter dominant)
  call precompute_postborn(dz/Hz,chi,chis,wlf,kl,pl,wp,ck)
  deallocate(chi,Hz,wlf)

end subroutine bispec_lens_pb_init


subroutine precompute_postborn(dchi,chi,chis,wlf,kl,pl,wp,ck)
! precomputing quantities relevant to the post-born correction
  implicit none
  !I/O
  double precision, intent(in) :: chis, dchi(:), chi(:), wlf(:), kl(:,:), pl(:,:)
  double precision, intent(out) :: wp(:,:), ck(:,:)
  !internal
  integer :: i, j, l, zn, ln
  double precision, allocatable :: fchi1(:), fchi2(:,:)

  zn = size(kl,dim=1)
  ln = size(kl,dim=2)

  allocate(fchi1(zn),fchi2(zn,zn)); fchi1=0d0; fchi2=0d0

  fchi1 = (chis-chi)/(chi*chis)

  do i = 1, zn
    do j = 1, zn
      if (chi(i)<=chi(j)) fchi2(i,j) = (chi(j)-chi(i))/(chi(j)*chi(i))
    end do
    wp(i,:) = pl(i,:) * ( wlf(i)/kl(i,:)**2 )**2 * dchi(i) * ((chis-chi(i))/(chis*chi(i)**2))**2
  end do

  do j = 1, zn
    do l = 2, ln
      ck(j,l) = dble(l)**4*sum(pl(:,l)*(wlf/kl(:,l)**2)**2*(dchi/chi**2)*fchi1*fchi2(:,j))
    end do
  end do

  deallocate(fchi1,fchi2)

end subroutine precompute_postborn


subroutine bispec_lens_pb(shap,eL,wp,ck,bl,sql0,ltype)
! lensing bispectrum from post-Born correction
  implicit none
  !I/O
  character(*), intent(in) :: shap
  integer, intent(in) :: eL(2)
  double precision, intent(in) :: wp(:,:), ck(:,:)
  double precision, intent(out) :: bl(:)
  !optional
  character(*), intent(in), optional :: ltype
  integer, intent(in), optional :: sql0
  !internal
  character(8) :: lt=''
  integer :: l, i, zn, l1, l2, l3, l0 = 1

  ! initial set up
  bl = 0d0
  zn = size(wp,dim=1)
  if (present(ltype)) lt = ltype
  if (present(sql0))  l0 = sql0

  ! compute bispectrum
  do l = eL(1), eL(2)

    if (shap=='fold'.and.mod(l,2)/=0) cycle
    if (shap=='sque'.and.2*l<l0)      cycle

    call set_three_ells(shap,l,eL,l0,l1,l2,l3)
    call bispec_lens_pb_kernel(l1,l2,l3,wp,ck,bl(l))

    !fullsky angular bispectrum
    if (lt=='full') bl(l) = bl(l) * W3j_approx(dble(l1),dble(l2),dble(l3)) * dsqrt((2d0*l1+1d0)*(2d0*l2+1d0)*(2d0*l3+1d0)/(4d0*pi))

  end do

end subroutine bispec_lens_pb


subroutine bispec_lens_pb_kernel(l1,l2,l3,wp,ck,bisp)
! compute lensing bispectrum from post-Born correction
  implicit none
  !I/O
  integer, intent(in) :: l1, l2, l3
  double precision, intent(in) :: wp(:,:), ck(:,:)
  double precision, intent(out) :: bisp
  !internal
  double precision :: al1, al2, al3, l1l2, l2l3, l3l1
  
  al1  = dble(l1)
  al2  = dble(l2)
  al3  = dble(l3)
  l1l2 = al3**2-al1**2-al2**2 ! 2L1*L2
  l2l3 = al1**2-al2**2-al3**2
  l3l1 = al2**2-al3**2-al1**2

  !Post-Born correction: Eq.(4.4) of PL16
  !wPp = dchi * W^2(chi,chi_s)/chi^2 * P_psi 
  bisp =        l1l2/(2d0*al1**2*al2**2)*(l3l1*al1**4*sum(wp(:,l1)*ck(:,l2))+l2l3*al2**4*sum(wp(:,l2)*ck(:,l1)))
  bisp = bisp + l2l3/(2d0*al2**2*al3**2)*(l1l2*al2**4*sum(wp(:,l2)*ck(:,l3))+l3l1*al3**4*sum(wp(:,l3)*ck(:,l2)))
  bisp = bisp + l3l1/(2d0*al3**2*al1**2)*(l2l3*al3**4*sum(wp(:,l3)*ck(:,l1))+l1l2*al1**4*sum(wp(:,l1)*ck(:,l3)))

end subroutine bispec_lens_pb_kernel


subroutine set_three_ells(shap,l,eL,l0,l1,l2,l3)
  ! choose a bispectrum configuration
  implicit none
  character(*), intent(in) :: shap
  integer, intent(in)  :: l, eL(2), l0
  integer, intent(out) :: l1, l2, l3

  select case(shap)
    case('equi')
      l1 = l
      l2 = l 
      l3 = l
    case('fold')
      l1 = l
      l2 = l/2
      l3 = l/2
    case('sque')
      l1 = l0
      l2 = l
      l3 = l
    case('angl')
      l1 = l
      l2 = int(eL(2)/2)
      l3 = int(eL(2)/2)
  end select


end subroutine set_three_ells


!//// binned bispectrum ////!

subroutine bispec_lens_bin(cp,b,eL1,eL2,eL3,wp,ck,m,btype,bl)
! reduced bispectrum with flat binning
  !I/O
  implicit none
  type(cosmoparams), intent(in) :: cp
  type(bispecfunc), intent(in) :: b
  character(*), intent(in) :: btype, m
  integer, intent(in) :: eL1(2), eL2(2), eL3(2)
  double precision, intent(in) :: wp(:,:), ck(:,:)
  double precision, intent(out) :: bl(2)
  !internal
  integer :: l1, l2, l3, zn
  double precision :: norm, bisp(2), hlll, tot(2)

  zn = size(b%kl,dim=1)

  tot  = 0d0
  norm = 0d0
  do l1 = eL1(1), eL1(2)
    do l2 = eL2(1), eL2(2)
      do l3 = eL3(1), eL3(2)

        if (l1==0) cycle
        if (l2==0) cycle
        if (l3==0) cycle
        if (l3>l1+l2.or.l3<abs(l1-l2)) cycle
        if (l1>l2+l3.or.l1<abs(l2-l3)) cycle
        if (l2>l3+l1.or.l2<abs(l3-l1)) cycle

        bisp = 0d0
        call bispec_lens_lss_kernel(cp,b,l1,l2,l3,bisp(1),m)
        call bispec_lens_pb_kernel(l1,l2,l3,wp,ck,bisp(2))

        !normalization
        hlll = W3j_approx(dble(l1),dble(l2),dble(l3)) * dsqrt((2d0*l1+1d0)*(2d0*l2+1d0)*(2d0*l3+1d0)/(4d0*pi))
        norm = norm + hlll**2

        !binned bispectrum
        tot  = tot + bisp*hlll**2

      end do
    end do
  end do

  bl = tot/norm

end subroutine bispec_lens_bin


subroutine bispec_gauss_bin(eL1,eL2,eL3,cl,f)
! reduced bispectrum for with flat binning (a=g+g^2)
  implicit none
  integer, intent(in) :: eL1(2), eL2(2), eL3(2)
  double precision, intent(in) :: cl(:)
  double precision, intent(out) :: f
  integer :: l1, l2, l3
  double precision :: norm, bisp, hlll, tot

  tot  = 0d0
  norm = 0d0
  do l1 = eL1(1), eL1(2)
    do l2 = eL2(1), eL2(2)
      do l3 = eL3(1), eL3(2)
        if (l3>l1+l2.or.l3<abs(l1-l2)) cycle
        if (l1>l2+l3.or.l1<abs(l2-l3)) cycle
        if (l2>l3+l1.or.l2<abs(l3-l1)) cycle
        bisp = 2d0*(Cl(l1)*Cl(l2)+Cl(l2)*Cl(l3)+Cl(l3)*Cl(l1))
        !normalization
        hlll = W3j_approx(dble(l1),dble(l2),dble(l3)) * dsqrt((2d0*l1+1d0)*(2d0*l2+1d0)*(2d0*l3+1d0)/(4d0*pi))
        norm = norm + hlll**2
        !binned bispectrum
        tot  = tot + bisp*hlll**2
      end do
    end do
  end do

  f = tot/norm

end subroutine bispec_gauss_bin


subroutine snr_bispec_lens(cp,b,eL,Cl,wp,ck,m,snr)
! SNR sum of lensing bispectrum
  implicit none
  type(cosmoparams), intent(in) :: cp
  type(bispecfunc), intent(in) :: b
  character(*), intent(in) :: m
  integer, intent(in) :: eL(2)
  double precision, intent(in) :: Cl(:), wp(:,:), ck(:,:)
  double precision, intent(out) :: snr
  !internal
  integer :: l1, l2, l3, i, zn
  double precision :: Del, bisp(2), tot

  zn = size(b%kl,dim=1)

  tot = 0d0
  do l1 = eL(1), eL(2)
    do l2 = l1, eL(2)
      do l3 = l2, eL(2)

        if (l3>l1+l2.or.l3<abs(l1-l2)) cycle
        if (l1>l2+l3.or.l1<abs(l2-l3)) cycle
        if (l2>l3+l1.or.l2<abs(l3-l1)) cycle
        if (mod(l1+l2+l3,2)==1) cycle
        Del = 1d0
        if (l1==l2.and.l2/=l3) Del = 2d0
        if (l1/=l2.and.l2==l3) Del = 2d0
        if (l1==l2.and.l2==l3) Del = 6d0

        bisp = 0d0
        call bispec_lens_lss_kernel(cp,b,l1,l2,l3,bisp(1),m)
        call bispec_lens_pb_kernel(l1,l2,l3,wp,ck,bisp(2))

        !flat sky -> full sky
        bisp = bisp * W3j_approx(dble(l1),dble(l2),dble(l3)) * dsqrt((2d0*l1+1d0)*(2d0*l2+1d0)*(2d0*l3+1d0)/(4d0*pi))

        !SNR
        tot = tot + sum(bisp)**2/(Del*Cl(l1)*Cl(l2)*Cl(l3))

      end do
    end do
  end do

  snr = dsqrt(tot)

end subroutine snr_bispec_lens


subroutine snr_bispec_lens_assym(cp,b,eL1,eL2,eL3,Cl,wp,ck,m,snr)
! SNR sum of lensing bispectrum with assymetry in l1,l2,l3
  implicit none
  type(cosmoparams), intent(in) :: cp
  type(bispecfunc), intent(in) :: b
  character(*), intent(in) :: m
  integer, intent(in) :: eL1(2), eL2(2), eL3(2)
  double precision, intent(in) :: Cl(:), wp(:,:), ck(:,:)
  double precision, intent(out) :: snr
  integer :: l1, l2, l3, i, zn
  double precision :: bisp(2), tot

  zn = size(b%kl,dim=1)

  tot = 0d0
  do l1 = eL1(1), eL1(2)
    do l2 = eL2(1), eL2(2)
      do l3 = eL3(1), eL3(2)
        if (l3>l1+l2.or.l3<abs(l1-l2)) cycle
        if (l1>l2+l3.or.l1<abs(l2-l3)) cycle
        if (l2>l3+l1.or.l2<abs(l3-l1)) cycle
        if (mod(l1+l2+l3,2)==1) cycle
        bisp = 0d0
        call bispec_lens_lss_kernel(cp,b,l1,l2,l3,bisp(1),m)
        call bispec_lens_pb_kernel(l1,l2,l3,wp,ck,bisp(2))
        !flat sky -> full sky
        bisp = bisp * W3j_approx(dble(l1),dble(l2),dble(l3)) * dsqrt((2d0*l1+1d0)*(2d0*l2+1d0)*(2d0*l3+1d0)/(4d0*pi))
        !SNR
        tot = tot + sum(bisp)**2/(6d0*Cl(l1)*Cl(l2)*Cl(l3))
      end do
    end do
  end do

  snr = dsqrt(tot)

end subroutine snr_bispec_lens_assym


subroutine snr_xbisp(cp,b,eL,cgg,ckk,btype,m,snr)
! SNR sum of gkk or ggk bispectrum
  implicit none
  type(cosmoparams), intent(in) :: cp
  type(bispecfunc), intent(in) :: b
  character(*), intent(in) :: btype, m
  integer, intent(in) :: eL(2)
  double precision, intent(in) :: cgg(:), ckk(:)
  double precision, intent(out) :: snr
  integer :: l1, l2, l3, i, zn
  double precision :: cov, Del, bisp, tot

  zn = size(b%kl,dim=1)
  tot = 0d0
  do l1 = eL(1), eL(2)
    if (mod(l1,10)==0) write(*,*) l1
    do l2 = eL(1), eL(2)
      do l3 = l2, eL(2)
        if (l3>l1+l2.or.l3<abs(l1-l2)) cycle
        if (l1>l2+l3.or.l1<abs(l2-l3)) cycle
        if (l2>l3+l1.or.l2<abs(l3-l1)) cycle
        if (mod(l1+l2+l3,2)==1) cycle
        Del = 3d0
        if (l2==l3) Del = 6d0
        bisp = 0d0
        call bispec_lens_lss_kernel(cp,b,l1,l2,l3,bisp,m)
        !flat sky -> full sky
        bisp = bisp * W3j_approx(dble(l1),dble(l2),dble(l3)) * dsqrt((2d0*l1+1d0)*(2d0*l2+1d0)*(2d0*l3+1d0)/(4d0*pi))
        !SNR
        if (btype=='gkk') cov = Del*cgg(l1)*ckk(l2)*ckk(l3)
        if (btype=='ggk') cov = Del*ckk(l1)*cgg(l2)*cgg(l3)
        tot = tot + bisp**2/cov
      end do
    end do
  end do

  snr = dsqrt(tot)

end subroutine snr_xbisp


subroutine Limber_k2l(chi,k,Pk,kl,Pl)
! get k and P(k) at k=l/chi
  implicit none
  double precision, intent(in) :: chi(:), k(:), Pk(:,:)
  double precision, intent(out) :: kl(:,:), Pl(:,:)
  integer :: i, l, zn, kn, ln, id
  double precision :: kk
  !double precision, allocatable :: y2a(:)

  zn = size(Pk,dim=1) !num of z points
  kn = size(Pk,dim=2) !num of k points
  ln = size(kl,dim=2) !num of multipole

  do i = 1, zn
    !allocate(y2a(kn))
    !call spline(k,Pk(i,:),kn,0d0,0d0,y2a)
    do l = 1, ln
      kk = dble(l)/chi(i)
      id = neighb(kk,k) !look for neighberest points
      kl(i,l) = kk
      !* linear interpolation
      Pl(i,l) = Pk(i,id) + (Pk(i,id+1)-Pk(i,id))*(kk-k(id))/(k(id+1)-k(id))
      !Pl(i,l) = splint(kk,k(i),Pk(i,:),y2a)
    end do
    !deallocate(y2a)
  end do

end subroutine Limber_k2l


function W3j_approx(l1,l2,l3) result(f)
! approximate W3j symbol
  implicit none
  double precision, intent(in) :: l1,l2,l3
  double precision :: a1,a2,a3,b,f,Lh

  if (mod(int(l1+l2+l3),2)/=0) then 
    f = 0d0
  else
    Lh = dble(l1+l2+l3)*0.5d0
    a1 = ((Lh-l1+0.5d0)/(Lh-l1+1d0))**(Lh-l1+0.25d0)
    a2 = ((Lh-l2+0.5d0)/(Lh-l2+1d0))**(Lh-l2+0.25d0)
    a3 = ((Lh-l3+0.5d0)/(Lh-l3+1d0))**(Lh-l3+0.25d0)
    b = 1d0/((Lh-l1+1d0)*(Lh-l2+1d0)*(Lh-l3+1d0))**(0.25d0)
    f = (-1d0)**Lh/dsqrt(2d0*pi) * exp(1.5d0)* (Lh+1d0)**(-0.25d0) * a1*a2*a3*b
  end if

end function W3j_approx


subroutine zinterp(zmin,zmax,zn,zspace,z,dz)
  ! precomputing interpolation points for z
  implicit none
  integer, intent(in) :: zn, zspace
  double precision, intent(in) :: zmin, zmax
  double precision, intent(out) :: z(:), dz(:)
  type(gauss_legendre_params) :: gl

  call gl_initialize(gl,zn,1d-15)
  select case (zspace)
  case(0)
    z  = linspace(zmin,zmax,zn)
    dz = z(2)-z(1)
  case(1)
    z  = glpoints([zmin,zmax],gl%z) ! points
    dz = gldxs([zmin,zmax],gl%w)    ! width
  end select
  call gl_finalize(gl)


end subroutine zinterp


subroutine lens_zweight(cp,z,dz,zs,zker,btype)
  !lensing weights
  implicit none
  type(cosmoparams), intent(in) :: cp
  double precision, intent(in)  :: z(:), dz(:), zs
  double precision, intent(out) :: zker(:)
  character(*), intent(in), optional :: btype
  character(16) :: bisptype = 'kkk'
  integer :: zn
  double precision :: chis, h
  double precision, allocatable :: Hz(:), chi(:), wlf(:)

  zn = size(z)
 
  allocate(Hz(zn),chi(zn),wlf(zn))
  chis = C_z(zs,cp) !source comoving distance
  chi  = C_z(z,cp)  !comoving distance at each z
  Hz   = H_z(z,cp)  !expansion rate at z
  h    = cp%H0/100d0

  ! lens kernel
  wlf  = 1.5d0*cp%Om*(cp%H0/3d5)**2*(1d0+z) !matter -> potential conversion factor (matter dominant)

  ! check errors
  if (chis<maxval(chi)) stop 'chis is smaller than maximum of chi in the chi integral'

  ! choose kernel for z integral
  if (present(btype)) bisptype = btype
  select case (bisptype)
  case ('kkk')
    zker = dz/Hz * (wlf*(chis-chi)*chi/chis)**3 / chi**4 
  case ('gkk')
    zker = dz/Hz * (wlf*(chis-chi)*chi/chis)**2 / chi**4 * Hz
  case ('ggk')
    zker = dz/Hz * (wlf*(chis-chi)*chi/chis) / chi**4 * Hz**2
  case default
    stop 'need z-kernel type'
  end select

  deallocate(Hz,chi,wlf)

end subroutine lens_zweight


subroutine get_pk(cp,z,k,pkL0,pkL,pk,ftype)
  !get linear and non-linear Pk at each z
  implicit none
  type(cosmoparams), intent(in) :: cp
  double precision, intent(in)  :: z(:), k(:), pkL0(:)
  double precision, intent(out) :: pkL(:,:), pk(:,:)
  character(*), intent(in), optional :: ftype
  character(16) :: fit='T12'
  integer :: i, zn
  double precision, allocatable :: D(:)

  zn = size(z)
  allocate(D(zn))
  D  = D_z(z,cp)  !growth factor
  if (present(ftype)) fit=ftype

  do i = 1, zn
    pkL(i,:) = D(i)**2*pkL0  !linear P(k,z) (i=1 -> z=0)
  end do

  if (fit=='TR')  pk = pkL  !use linear 
  if (fit/='TR')  call NonLinRatios(pkL,z,k,cp,pk,ftype=fit) !nonlinear Pk

  !call savetxt('pklin.dat',ki/h,pklini(pkout,:)*h**3,ow=.true.)
  !call savetxt('pk.dat',ki/h,pki(pkout,:)*h**3,ow=.true.)
  deallocate(D)

end subroutine get_pk


subroutine get_knl(k,PkL,knl)
! get k_NL
  implicit none
  double precision, intent(in) :: k(:), PkL(:,:)
  double precision, intent(out) :: knl(:)
  integer :: kn, zn, i, j
  double precision :: f

  kn = size(k)
  zn = size(PkL,dim=1)

  do i = 1, zn
    do j = 1, kn
      if ( PkL(i,j)*k(j)**3/(2d0*pi**2) > 1d0 ) then
        goto 11
      else
        knl(i)  = k(j)
      end if
    end do
11  continue
  end do

end subroutine get_knl


subroutine get_neff(k,kl,plin,n,model)
! compute neff(k) at k=l/chi
  implicit none
  character(*), intent(in) :: model
  double precision, intent(in) :: k(:), kl(:,:), plin(:,:)
  double precision, intent(out) :: n(:,:)
  integer :: i, l, id, zn, ln, lh, h, l0, l1
  double precision :: v
  double precision, allocatable :: ki(:), ni(:), y(:), x(:), y2(:)

  zn = size(kl,dim=1)
  ln = size(kl,dim=2)

  allocate(ki(ln),ni(ln)); ki=0d0; ni=0d0

  do i = 1, zn
    do l = 2, ln
      ! find nearest index
      id     = neighb(kl(i,l),k)
      ! linear interpolation
      n(i,l) = (plin(i,id+1)-plin(i,id-1))/plin(i,id) * k(id)/(k(id+1)-k(id-1))
    end do
    !if (i==zn/10) call savetxt('test.dat',kl(i,:),n(i,:),ow=.true.)
    if (model=='GM') then !smoothed n_eff
      l0 = neighb(2d-2,kl(i,:)) !ell at BAO kmin
      l1 = neighb(4d-1,kl(i,:)) !ell at BAO kmax
      !find local maxima/minima
      lh = 0
      h  = -1
      v  = n(i,l0)
      do l = l0, l1
        if (h*n(i,l)>=h*v) then
          v  = n(i,l)
        else
          lh = lh + 1
          ki(lh) = kl(i,l)
          ni(lh) = n(i,l)
          h  = -h
        end if
      end do
      if (lh<2)  cycle
      !interpolation points
      allocate(y(lh+1),x(lh+1),y2(lh+1))
      x(1)    = kl(i,l0)
      y(1)    = n(i,l0)
      x(2:lh) = (ki(2:lh)+ki(1:lh-1))*0.5d0
      y(2:lh) = (ni(2:lh)+ni(1:lh-1))*0.5d0
      x(lh+1) = kl(i,l1)
      y(lh+1) = n(i,l1)
      call spline(x,y,lh+1,0d0,0d0,y2)
      !smoothing
      do l = l0, l1
        n(i,l) = splint(kl(i,l),x,y,y2)
      end do
      deallocate(x,y,y2)
      !if (i==zn/10) call savetxt('test2.dat',kl(i,:),n(i,:),ow=.true.)
    end if
  end do
 
  deallocate(ki,ni)

end subroutine get_neff


subroutine precompute_W3j(lmin,W3)
! precomputing W3j symbols
! require large memory for high l
  implicit none
  integer, intent(in) :: lmin
  double precision, intent(out) :: W3(:,:,:)
  integer :: lmax, l1, l2, l3

  lmax = size(W3,dim=1)

  do l1 = lmin, lmax
    do l2 = l1, lmax
      do l3 = l2, lmax
        if (l3>l1+l2.or.l3<abs(l1-l2)) cycle
        if (l1>l2+l3.or.l1<abs(l2-l3)) cycle
        if (l2>l3+l1.or.l2<abs(l3-l1)) cycle
        if (mod(l1+l2+l3,2)==1) cycle
        W3(l1,l2,l3) = W3j_approx(dble(l1),dble(l2),dble(l3)) * dsqrt((2d0*l1+1d0)*(2d0*l2+1d0)*(2d0*l3+1d0)/(4d0*pi))
      end do
    end do
  end do

end subroutine precompute_W3j


subroutine prep_lens_aps(z,dz,zs,cp,ki,pklin0,ck)
  implicit none

  ![input]
  ! z, dz  --- redshift points and thier interval
  ! zs     --- source z
  ! cp     --- cosmological parameters
  ! ki     --- CAMB output k
  ! pklin0 --- CAMB output P(k,z=0)
  type(cosmoparams), intent(in) :: cp
  double precision, intent(in)  :: z(:), dz(:), zs, ki(:), pklin0(:)

  ![output]
  ! ck     --- kappa power spectrum at chis_s
  double precision, intent(out) :: ck(:)

  !internal
  integer :: i, zn, kn, ln, l
  double precision :: chis, h
  double precision, allocatable :: Hz(:), chi(:), D(:), wlf(:), pki(:,:), pklini(:,:), pl(:,:), kl(:,:)

  zn = size(z)  !number of z points for z-integral
  kn = size(ki) !number of CAMB output k
  ln = size(ck)

  !* get distances and lensing weights
  allocate(Hz(zn),chi(zn),D(zn),wlf(zn),pklini(zn,kn),pki(zn,kn))
  D    = D_z(z,cp)  !growth factor
  chis = C_z(zs,cp) !source comoving distance
  chi  = C_z(z,cp)  !comoving distance at each z
  Hz   = H_z(z,cp)  !expansion rate at z
  h    = cp%H0/100d0

  wlf  = 1.5d0*cp%Om*(cp%H0/3d5)**2*(1d0+z) !matter -> potential conversion factor (matter dominant)

  if (2d0/chi(zn)<ki(1)) write(*,*) 'warning: required minimum k is smaller than input', 2d0/chi(zn), ki(1)
  if (chis<maxval(chi)) stop 'chis is smaller than maximum of chi in the chi integral'

  !* get linear and non-linear Pk at each z
  do i = 1, zn
    pklini(i,:) = D(i)**2*pklin0  !linear P(k,z) (i=1 -> z=0)
  end do
  call NonLinRatios(pklini,z,ki,cp,pki) !nonlinear Pk

  ! check nonlinear Pk
  !call savetxt('Pklin.dat',ki/h,pklini(1,:)*h**3,pklini(zn,:)*h**3,ow=.true.)
  !call savetxt('Pk.dat',ki/h,pki(1,:)*h**3,pki(zn,:)*h**3,ow=.true.)

  !* interpolate k, Pk at k=l/chi
  allocate(kl(zn,ln),pl(zn,ln))
  call Limber_k2l(chi,ki,pki,kl,pl)  

  !* kappa aps
  ck = 0d0
  do l = 2, ln
    ck(l) = sum(pl(:,l)*(dz/Hz)*(wlf*(chis-chi)*chi/chis)**2/chi**2)
  end do

  deallocate(Hz,chi,D,wlf,pklini,pki,pl,kl)

end subroutine prep_lens_aps


end module cmblbisp


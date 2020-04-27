!////////////////////////////////////////////////////!
! * Kernel Functions of Normalization
!////////////////////////////////////////////////////!

subroutine Kernels_Tau(rL,WA,WB,X,kernel,gln,gle)
  implicit none
  !I/O
  character(*), intent(in) :: kernel
  integer, intent(in) :: rL(2)
  integer, intent(in), optional :: gln
  double precision, intent(in), optional :: gle
  double precision, intent(in), dimension(rL(1):rL(2)) :: WA, WB
  double precision, intent(out) :: X(:)
  !internal
  type(gauss_legendre_params) :: GL
  integer :: i, l, lmax, gn
  double precision :: mu, II, al, c1_inv, c3, ge
  double precision, dimension(2) :: d00_sup, d00_mid, d00_inf
  double precision, dimension(2) :: ZA00, ZB00, ZA22, ZB22

  !* initialize
  lmax = size(X)

  X = 0d0

  !* GL quadrature
  gn = int((3*max(rL(2),lmax)+1)/2)
  ge = 1d-15
  if (present(gln)) gn = gln
  if (present(gle)) ge = gle
  call gl_initialize(GL,gn,ge)

  do i = 1, GL%n !Gauss-Legendre Integration
    mu = GL%z(i)
    select case(kernel)
    case ('S0','G0')
      call ZETA(0,0,rL,WA,mu,ZA00)
      call ZETA(0,0,rL,WB,mu,ZB00)
      II = 2d0*ZA00(1)*ZB00(1)
    case default
      stop 'error: no kernel'
    end select
    do l = 1, lmax ! loop for l
      al = dble(l)
      if(l==1) then 
        d00_sup = wigd_ini(0,0,mu)
        d00_inf = d00_mid
        d00_mid = d00_sup
        d00_sup(1:2) = mu*d00_mid(1:2)*(al*(2d0*al-1d0))/al**2
      else
        c1_inv = (2d0*al-1d0)/al
        c3 = (al-1d0)**2/((al-1d0)*(2d0*al-1d0))
        d00_sup(1) = (mu*d00_mid(1) - c3*d00_inf(1))*c1_inv
        d00_sup(2) = (mu*d00_mid(2) - c3*d00_inf(2))*c1_inv
      end if
      X(l) = X(l) + II*d00_sup(1)*GL%w(i)*pi
      d00_inf = d00_mid
      d00_mid = d00_sup
    end do
  end do

  call gl_finalize(GL)

end subroutine Kernels_Tau


subroutine qtt(lmax,rlmin,rlmax,fC,OCT,At)
!*  Normalization of reconstructed amplitude modulation from the temperature quadratic estimator
!*
!*  Args:
!*    :lmax (int)        : Maximum multipole of output normalization spectrum
!*    :rlmin/rlmax (int) : Minimum/Maximum multipole of CMB for reconstruction
!*    :fC [l] (double)   : Theory TT spectrum, with bounds (0:rlmax)
!*    :OCT [l] (double)  : Observed TT spectrum, with bounds (0:rlmax)
!*
!*  Returns:
!*    :At [l] (double) : tau normalization, with bounds (0:lmax)
!*
  implicit none
  !I/O
  integer, intent(in) :: lmax, rlmin, rlmax
  double precision, intent(in), dimension(0:rlmax) :: fC, OCT
  double precision, intent(out), dimension(0:lmax) :: At
  !internal
  integer :: rL(2), l
  double precision, dimension(rlmin:rlmax) :: W1, W2
  double precision, dimension(min(2*rlmax,lmax)) :: S0, G0

  write(*,*) 'norm qTT (tau)'
  rL = (/rlmin,rlmax/)

  do l = rlmin, rlmax
    if (OCT(l)==0d0) stop 'error (norm_tau.qtt): observed cltt is zero'
  end do

  !filtering functions
  W1 = 1d0 / OCT(rlmin:rlmax)

  !main calculation
  W2 = W1 * fC(rlmin:rlmax)**2
  S0 = 0d0
  call Kernels_tau(rL,W1,W2,S0,'S0')

  W2 = W1 * fC(rlmin:rlmax)
  G0 = 0d0
  call Kernels_tau(rL,W2,W2,G0,'G0')

  At = 0d0
  do l = 1, lmax
    if (S0(l)+G0(l)/=0d0)  At(l) = 1d0/(S0(l)+G0(l))
  end do

end subroutine qtt




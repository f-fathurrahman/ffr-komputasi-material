!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MODULE EXPO
! ===========
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module expo
  use mesh
  implicit none

  private
  public :: exponential, pot

  integer, parameter :: TAYLOR = 1,  &
                        LANCZOS = 2, &
                        EXPOKIT = 3

  integer, parameter :: mode = LANCZOS
  !integer, parameter :: mode = EXPOKIT

  real(8), allocatable :: pot(:, :)

  contains




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE EXPONENTIAL
! ======================
!
! INPUT:
!   v [real(8), dimension(n, n)] : the potential that defines the Hamiltonian.
!   wf [complex(8), dimension(n, n)] : on input, it contains the wavefunction
!      onto which the exponential of H is applied.
!   dt [real(8)] : the subroutine should return exp[-I*dt*H]|wf>
! --------
! OUTPUT:
!   wf [complex(8), dimension(n, n)] : on output, it should contain the result
!
! This subroutine is merely a driver to the appropriate subroutine, depending
! on the value of the mode parameter, which should be set up above.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine exponential(v, wf, dt)
  implicit none
  real(8), intent(in) :: v(N, N)
  complex(8), intent(inout) :: wf(N, N)
  real(8), intent(in) :: dt

  select case(mode)
    case(TAYLOR);  call exponential_taylor(v, wf, dt)
    case(LANCZOS)
        call exponential_lanczos(v, wf, dt)
    case(EXPOKIT)
        allocate(pot(n, n))
        pot = v
        call exponential_expokit(wf, dt)
        deallocate(pot)
  end select

end subroutine exponential

subroutine exponential_taylor(v, wf, dt)
  implicit none
  real(8), intent(in) :: v(N, N)
  complex(8), intent(inout) :: wf(N, N)
  real(8), intent(in) :: dt

  integer, parameter :: order = 4
  complex(8) :: zfact
  complex(8) :: zpsi(n, n), hzpsi(n, n)
  integer :: i

  zfact = (1.0_8, 0.0_8)
  zpsi = wf
  do i = 1, order
     zfact = zfact*(0.0_8, -1.0_8)*dt/i
     call zhpsi(v, zpsi, hzpsi)
     wf = wf + zfact*hzpsi
     if(i.ne.order) zpsi = hzpsi
  enddo

end subroutine exponential_taylor

subroutine exponential_lanczos(v, wf, dt)
  implicit none
  real(8), intent(in) :: v(N, N)
  complex(8), intent(inout) :: wf(N, N)
  real(8), intent(in) :: dt

  integer :: i
  integer, parameter :: korder = 20
  real(8), parameter :: tol = 1.0e-5_8
  complex(8), allocatable :: x(:, :, :), hm(:, :), expo(:, :), mat(:, :)
  real(8) :: nrm, alpha, beta, res

  complex(8), allocatable :: wsp(:)
  integer :: lwsp, ipiv(korder), iexph, ns, iflag, order

  lwsp = 4*korder*korder+7
  allocate(wsp(lwsp))

  allocate(x(n, n, korder), hm(korder, korder), expo(korder, korder), mat(korder, korder))

  nrm = sqrt(zdotproduct(wf, wf))
  x(:, :, 1) = wf(:, :)/nrm

  call zhpsi(v, x(:, :, 1), wf)
  alpha = zdotproduct(x(:, :, 1), wf)
  wf(:, :)= wf(:, :) - alpha*x(:, :, 1)

  hm = (0.0d0, 0.0d0)
  hm(1, 1) = alpha  
  beta = sqrt(zdotproduct(wf, wf))

  do i = 1, korder - 1
     x(:, :, i + 1) = wf(:, :)/beta
     hm(i+1, i) = beta
     call zhpsi(v, x(:, :, i+1), wf)

     hm(i    , i + 1) = zdotproduct(x(:,:, i)  , wf)
     hm(i + 1, i + 1) = zdotproduct(x(:,:, i+1), wf)
     call zgemv('n', n**2,  2, (-1.0d0,0.0d0), x(1, 1, i), n**2, hm(i:i+1, i+1), 1, (1.0d0, 0.0d0), wf, 1)
     mat = (0.0d0, -1.0d0)*hm
     call ZGPADM(6, i+1, dt, mat, korder, wsp, lwsp, ipiv, iexph, ns, iflag)
     call zcopy((i+1)*(i+1), wsp(iexph), 1, expo(1:i+1, 1:i+1), 1)

     res = abs(beta*abs(expo(1, i+1)))
     beta = sqrt(zdotproduct(wf, wf))

     if(beta < 1.0e-12_8) exit
     if(i>1 .and. res<tol) exit

  enddo

  order = min(korder, i+1)
  if(res> tol .and. beta > 1.0e-12_8) then
     write(*, '(1x,a,2ES18.10)') 'Warning: Lanczos exponential expansion did not converge: ', res, beta
  endif

  call zgemv('n', n**2, korder, cmplx(nrm, 0.0d0, 8), x, n**2, expo, 1, cmplx(0.0d0, 0.0d0, 8), wf, 1)

  deallocate(x, hm, expo, wsp)

end subroutine exponential_lanczos

subroutine exponential_expokit(wf, dt)
  implicit none
  complex(8), intent(inout) :: wf(N, N)
  real(8), intent(in) :: dt

      integer nn, m, lwsp, liwsp
      double precision tol, anorm, s1, s2
      complex*16, allocatable :: x(:), y(:)
      integer, allocatable :: iwsp(:)
      complex*16, allocatable :: wsp(:)
      
      integer i, j, nnz, itrace, iflag, iseed(4)
      complex*16 ZERO, ONE
      parameter( ZERO=(0.0d0,0.0d0), ONE=(1.0d0,0.0d0) )

      external op

      nn = n**2
      m = 8
      lwsp = nn*(m+2)+5*(m+2)**2+7
      liwsp = m+2

      allocate(x(nn), y(nn))
      allocate(iwsp(liwsp), wsp(lwsp))

      call zcopy(nn, wf, 1, x, 1)

      anorm = 1.0d0

!*---  set other input arguments ...
      tol = 1.0e-3_8
      itrace = 0

!*---  compute w = exp(t*A)v with ZGEXPV ...
      call ZGEXPV( nn, m, dt, x, y, tol, anorm, &
                   wsp, lwsp, iwsp, liwsp, op, itrace, iflag )

      call zcopy(nn, y, 1, wf, 1)

end subroutine exponential_expokit

end module expo

subroutine op(u, v)
  use mesh
  use expo
  complex(8), intent(in)  :: u(*)
  complex(8), intent(out)  :: v(*)
  complex(8), allocatable :: a(:, :), b(:, :)
  allocate(a(n, n), b(n, n))
  call zcopy(n*n, u, 1, a, 1)
  call zhpsi(pot, a, b)
  call zscal(n*n, (0.0d0, -1.0d0), b, 1)
  call zcopy(n*n, b, 1, v, 1)
  deallocate(a, b)
end subroutine op

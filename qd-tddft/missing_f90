




!!!!!! MISSING CODE 1
  integer :: ix, iy, k

  integer, parameter :: order = 4
  real(8), allocatable :: c(:)

  allocate(c(-order:order))
  c(-order:order) = (/ &  ! The coefficients of the laplacian...
      -1.785714d-3, 2.539683d-2, -0.2d0, 1.6d0,      &
      -2.847222d0,                                   &
      1.6d0, -0.2d0, 2.539683d-2, -1.785714d-3 /)

  do ix = 1, N
    do iy = 1, N
      lapl(ix, iy) = (0.0_8, 0.0_8)
      do k = -order, order
        if(iy+k>=1 .and. iy+k<=N) lapl(ix, iy) = lapl(ix, iy) + c(k)*f(ix, iy+k)
        if(ix+k>=1 .and. ix+k<=N) lapl(ix, iy) = lapl(ix, iy) + c(k)*f(ix+k, iy)
      end do
    end do
  end do

  lapl = lapl/delta**2
!!!!!! END OF MISSING CODE




!!!!!! MISSING CODE 2
! This defines a harmonic potential:
  integer :: ix, iy
  real(8) :: r2, omega, a, b
  omega = 0.22_8
  do ix = 1, N
  do iy = 1, N
     r2 = x(ix, iy)**2 + y(ix, iy)**2
     v(ix, iy) = 0.5_8*omega**2*r2
  enddo
  enddo
!!!!!! END OF MISSING CODE




!!!!!! MISSING CODE 2
! This defines a quartic potential:
  integer :: ix, iy
  real(8) :: r2, omega, a, b
  a = 0.00008_8
  do ix = 1, N
  do iy = 1, N
     v(ix, iy) = a*x(ix,iy)**4 + a*y(ix, iy)**4
  enddo
  enddo
!!!!!! END OF MISSING CODE






!!!!!! MISSING CODE 3
  do i = 1, max_iter

     call conjugate_gradients(N_wf, wfs, vtot, eigenval, residues)

     call external_pot(vext)
     rhoprev = rho
     call build_rho(wfs, rho)

     diff = sqrt(dotproduct(rho-rhoprev,rho-rhoprev))
     if(diff < tol) exit

     rho = alpha*rho + (1.0 - alpha)*rhoprev
     call interaction_pot(rho, vh, vx, vc, ex, ec)
     vtot = vext + vh + vx + vc

     etot =  energy(eigenval, rho, vh, vc, vx, ec, ex)

     write(*, '(/,a,i6)') 'SCF CYCLE ITER # ',i
     write(*, '(5x,a,es18.4)') 'diff = ',diff
     do j = 1, N_wf
        write(*, '(5x,i4, 2es20.8)') j, eigenval(j), residues(j)
     enddo
     write(*,'(a)')

  enddo
!!!!!! END OF MISSING CODE






!!!!!! MISSING CODE 4
  integer :: mode
  integer, parameter :: INDEPENDENT_PARTICLES = 0, &
                        HARTREE               = 1, &
                        HARTREE_X             = 2, &
                        HARTREE_XC            = 3

  mode = HARTREE

  select case(mode)
    case(INDEPENDENT_PARTICLES)
      vh = 0.0_8; vx = 0.0_8; vc = 0.0_8; ex = 0.0_8; ec = 0.0_8
    case(HARTREE)
      vx = 0.0_8; vc = 0.0_8; ex = 0.0_8; ec = 0.0_8
      call poisson_solve(rho, vh)
    case(HARTREE_X)
      call vxc_lda(rho, vx, vc, Ex, Ec)
      vc = 0.0_8; ec = 0.0_8
      call poisson_solve(rho, vh)
    case(HARTREE_XC)
      call vxc_lda(rho, vx, vc, Ex, Ec)
      call poisson_solve(rho, vh)
  end select
!!!!!! END OF MISSING CODE




!!!!!! MISSING CODE 5
  call laplacian(f, hf)
  hf(:, :) = -0.5_8*hf(:, :) + v(:, :)*f(:, :)
!!!!!! END OF MISSING CODE




!!!!!! MISSING CODE 6
  call zlaplacian(f, hf)
  hf(:, :) = -0.5_8*hf(:, :) + v(:, :)*f(:, :)
!!!!!! END OF MISSING CODE



!!!!!! MISSING CODE 7
  integer  :: ix, iy, jx, jy
  real(8)  :: r1(2), r2(2)
  real(8), parameter :: pi = 3.141592653589793d0

  v = 0.0d0
  do ix = 1, N
    do iy = 1, N
      r1(1) = x(ix, iy); r1(2) = y(ix, iy)
      do jx = 1, N
        do jy = 1, N
           r2(1) = x(jx, jy); r2(2) = y(jx, jy)

          if(ix == jx .and. iy == jy) then
            v(ix, iy) = v(ix, iy) + 2.0d0*sqrt(pi)*rho(ix, iy)/delta
          else
            v(ix, iy) = v(ix, iy) + rho(jx, jy)/sqrt(sum((r1-r2)**2))
          end if
        end do
      end do
      v(ix, iy) = v(ix, iy)*delta**2
    end do
  end do
!!!!!! END OF MISSING CODE

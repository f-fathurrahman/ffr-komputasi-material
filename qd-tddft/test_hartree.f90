subroutine test_hartree
  use mesh
  use fft
  use poisson
  implicit none

  integer :: ix, iy
  real(8), external :: besseli
  real(8), allocatable :: rho(:, :), vhsum(:, :), vhexact(:, :), vhfft(:, :)
  real(8) :: r2, alpha, z, ex
  real(8), parameter :: pi = 3.141592653589793_8

  allocate(rho(n, n), vhsum(n, n), vhexact(n, n), vhfft(n, n))


  alpha = 3.0
  do ix = 1, n
     do iy = 1, n
        r2 = x(ix, iy)**2 + y(ix, iy)**2
        rho(ix, iy) = exp(-r2/alpha**2)
     enddo
  enddo
  rho(:, :) = rho(:, :)/(alpha**2*pi)
  call output(rho,'rho.dat')

  write(*, *) 'Integrated charge: ', sum(rho(:, :))*delta**2, sum(rho(:, :))

  call fft_all_init()

    write(*, *) 'Calling poisson_init...'
  call poisson_init()
   write(*, *) 'Calling poisson_solve...'
  call poisson_sum(rho, vhsum)
   write(*, *) 'Done.'
  call output(vhsum, 'vhsum.dat')
   write(*, *) 'Calling poisson_fft...'
  call poisson_fft(rho, vhfft)
   write(*, *) 'Done.'
  call output(vhfft, 'vhfft.dat')


  do ix = 1, n
     iy = N/2+1
     r2 = x(ix, iy)**2
     z = r2/(2.0_8*alpha**2)
        if(z<200.0_8) then
           ex = (sqrt(pi)/alpha)*exp(-z)*besseli(0, z)
           write(71, '(i6,5es20.12)') ix, x(ix, iy),  rho(ix, iy), vhsum(ix, iy)-ex , vhfft(ix, iy)-ex, ex
        else
           ex = (sqrt(pi)/alpha)*(1.0_8/sqrt(2.0_8*pi*z))*(1.0_8+(1.0/(8.0*z)))
           write(71, '(i6,5es20.12)') ix, x(ix, iy),  rho(ix, iy), vhsum(ix, iy)-ex , vhfft(ix, iy)-ex, ex
        endif     
  enddo

  deallocate(rho, vhsum, vhfft, vhexact)
end subroutine

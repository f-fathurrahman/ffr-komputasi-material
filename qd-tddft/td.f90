subroutine td
  use states
  use mesh
  use fft
  use poisson
  implicit none

  integer :: i


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This is to permit the use of FFTs through the FFTW package.
  write(*, '(a)') 'Initializing FFTs...'
  call fft_all_init()
  write(*, '(a)') 'Done.'


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This initializes the solver for the Hartree term (it needs to be done only
  ! if FFTs are to be used to solve it.
  write(*, '(a)') 'Initializing the Hartree solver...'
  call poisson_init()
  write(*, '(a,/)') 'Done.'


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! The next subroutine gets the wavefunctions from a file called "wfs", that
  ! should contain the wavefunctions obtained from a ground-state calculation.
  write(*, *) 'Reading wavefunctions...'
  call read_wfs('wfs.dat')
  write(*, *) 'Done.'

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! In file wfs, the functions are real wavefunctions. To do the time
  ! propagation, one needs to use complex wavefunctions. So the wavefunctions
  ! are passed to variable zwfs, which is complex.
  allocate(zwfs(n, n, n_wf))
  zwfs = (0.0_8, 0.0_8)
  zwfs = wfs
  deallocate(wfs)

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Subroutine perturbation modifies the KS wavefunctions in some manner
  write(*, *) 'Applying an initial perturbation to the wavefunctions...'
  call perturbation(N_wf, zwfs)
  write(*, *) 'Done.'

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! The real work is now performed by the propagate subroutine.
  call propagate()

  ! clean up and leave.
  deallocate(zwfs)
end subroutine td


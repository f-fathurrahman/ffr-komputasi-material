!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PROGRAM GS
! ==========
!
! This program obtains solves the ground state DFT problem, i.e. it solves the
! Kohn-Sham equations. The result is a file called "wfs" where the ground
! state orbitals are placed, which may then be used by "td" or "excitations",
! to perform TDDFT calculations.
!*/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine gs
  use states
  use poisson
  use mesh
  use fft
  implicit none


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! First of all, specify how many states (occupied and unoccupied) are to be
  ! calculated. The allocates the wfs variable, that holds the wavefunctions
  ! throught all the program.
  N_occ = 3 ! should be an odd number? To avoid getting a metallic system
  N_empty = 2
  N_wf = N_occ + N_empty
  allocate(wfs(N, N, N_wf))


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This is to permit the use of FFTs through the FFTW package.
  !write(*, '(a)') 'Initializing FFTs...'
  !call fft_all_init()
  !write(*, '(a)') 'Done.'


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This initializes the solver for the Hartree term (it needs to be done only
  ! if FFTs are to be used to solve it.
  !write(*, '(a)') 'Initializing the Hartree solver...'
  !call poisson_init()
  !write(*, '(a,/)') 'Done.'


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! Puts random values onto the wavefunctions, to have a starting point for the
  ! following iterations.
  call init_wfs()


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! This is the core of the program, the self-consistent field cycle that
  ! solves the Kohn-Sham system of equations.
  call scf()


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! The obtained wavefunctions are to be used by other programs, such as td.
  ! We write them into a file called 'wfs'.
  call write_wfs('wfs.dat')

  ! clean up
  deallocate(wfs)
  call fft_all_end()
end subroutine gs

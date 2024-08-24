SUBROUTINE band_dft
!  self consistent dft calculation
  use fft_data
  use Gvector
  use PW
  use PW_SMALL
  use pseudopotential
  IMPLICIT NONE 
  INTEGER             :: iteration,ik,i
  double precision   :: rho((N_L(1)+fftinc1),N_L(2),N_L(3))

  write(6,*)'???'
  open(1,file='density.dat')
    read(1,*)store_density
  close(1)
  rho=0.
  store_density(N_L(1)+1:N_L(1)+fftinc1,:,:)=0.d0

  rho=store_density
  call calculate_potential(rho)

  call hamiltonian_small(rho)

!  rho=store_density
!  call calculate_potential(rho)


  do iteration=1,100
    write(6,*)iteration,E_total
!   pseudopotential acting on the wave function
    call nl_pp
!   approximate diagonalization
!    call steepest_descent(rho)
    call sw(rho)
  ENDDO 
    do ik=1,N_k_points
      write(6,'(''eigenvalues at k-point: '',3f6.3)') (k_point(i,ik),i=1,3)
      write(6,'(f7.3,9f8.3)') (eigenvalue(i,ik)*2.d0*rydberg,i=1,n_orbitals)
      write(10,100)ik,(eigenvalue(i,ik)*2.d0*rydberg,i=1,n_orbitals)
    ENDDO 
  100 format(i5,40e12.4)
end








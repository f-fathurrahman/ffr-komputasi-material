SUBROUTINE solve_small
  use fft_data
  USE GVECTOR
  USE PW_SMALL
  USE PW
  USE PSEUDOPOTENTIAL
  IMPLICIT NONE 

  INTEGER           :: iteration,i1,i2,i3
  double precision :: rho(N_L(1)+fftinc1,N_L(2),N_L(3))
  



!      call phase_factor
!      call structure_factor
!      call form_factor
!      call nl_pp_form_factor


!     jellium 
      rho=N_electrons/Volume
      store_density=rho

      call calculate_potential(rho)

      do  iteration=1,N_init
          call hamiltonian_small(rho)
          call calculate_density(rho)
          call calculate_potential(rho)
      ENDDO 

end SUBROUTINE solve_small





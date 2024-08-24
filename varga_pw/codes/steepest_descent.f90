SUBROUTINE steepest_descent(rho)
!   steepest_descent step on the wave function
use fft_data
use Gvector
use PW
use PW_SMALL
use pseudopotential
IMPLICIT NONE 

      
double precision            ::  dt,delt,gamma
double precision,parameter  ::  emass=33.333d0
complex*16                  ::  s, eig
INTEGER                      ::  ik,i,ig,j
double precision            ::  rho((N_L(1)+fftinc1),N_L(2),N_L(3))

  delt=24.d0
  gamma=0.2d0
  dt = delt/emass

  do ik=1,N_k_points 
    do i=1,n_orbitals
!   local part of the hamiltonian acting on the wave function
      call h_psi(i,ik,rho)
!   eigenvalue 
      eig=(0.0, 0.0)
      do ig=1,n_g_vector(ik)
        eig = eig - conjg(wave_function_c(ig,i,ik))*Hpsi(ig)
      ENDDO 
      eigenvalue(i,ik) = real(eig)

      do j=1,n_g_vector(ik)
        s=wave_function_c(j,i,ik)
        wave_function_c(j,i,ik)=s+dt*(Hpsi(j)+eigenvalue(i,ik)*s)
        s=wave_function_c(j,i,ik)-s
      ENDDO 
    ENDDO 
        
    call gram_schmidt(ik,n_orbitals,n_g_vector(ik))        
  ENDDO 

end SUBROUTINE steepest_descent

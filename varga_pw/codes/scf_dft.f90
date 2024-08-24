SUBROUTINE scf_dft
!  self consistent dft calculation
  use fft_data
  use Gvector
  use PW
  use PW_SMALL
  use pseudopotential
  IMPLICIT NONE 
  INTEGER             :: iteration,ik,i,i1,i2,i3
  double precision   :: rho((N_L(1)+fftinc1),N_L(2),N_L(3))


  do iteration=1,40
    if(dens_mix.gt.0.10d0) dens_mix=dens_mix-0.1d0
    write(16,*)'dens_mix',dens_mix
!   pseudopotential acting on the wave function
    call nl_pp
    call calculate_density(rho)
!    if(iteration.eq.1) read(100,*)rho

!    write(300,*)'?',iteration,rho
    call calculate_potential(rho)
!   approximate diagonalization

    call sw(rho)
!    call steepest_descent(rho)
!
    do ik=1,N_k_points
      write(16,'(''eigenvalues at k-point: '',3f6.3)') (k_point(i,ik),i=1,3)
      write(16,'(f7.3,9f8.3)') (eigenvalue(i,ik)*2.d0*rydberg,i=1,n_orbitals)
      write(16,*)'Occupation'
      write(16,'(f7.3,9f8.3)') (occupation(i,ik),i=1,n_orbitals)
    ENDDO 

!c\vk    call fermi

    write(6,*)iteration,E_total
  ENDDO 
  open(1,file='density.dat')
    write(1,*)store_density
  close(1)

end



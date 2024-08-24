SUBROUTINE sw(rho)
use fft_data
use Gvector
use PW
use PW_SMALL
use PSEUDOPOTENTIAL
IMPLICIT NONE 

double precision            ::  dt,delt,gamma,dts
double precision,parameter  ::  emass=33.333d0
complex*16                  ::  s, eig
INTEGER                      ::  ik,i,ig,j,lm,is
double precision            ::  rho((N_L(1)+fftinc1),N_L(2),N_L(3))
double precision            ::  diag,ww,vl,w
double precision            ::  vnl(N_G_K_vector_max),tau(N_G_K_vector_max)



  delt=40.d0
  gamma=0.2d0
  dts = delt
!  <g+k|H|g+k> 
! local potential
  vl=sum(rho)/(N_L(1)*N_L(2)*N_L(3))

  do ik=1,N_K_points

! nonlocal potential
    vnl=0.d0
    do is=1,N_species
       do lm=1,l_max(is)**2+1-2*l_loc(is)
         ww=wnl(is,lm)*N_atom(is)
         do ig=1,N_G_vector(ik)
           vnl(ig)=vnl(ig)+pkg_a(lm,ig,is,ik)**2*ww
         ENDDO 
       ENDDO 
    ENDDO 

    do i=1,N_orbitals
!   local part of the hamiltonian acting on the wave function
      call h_psi(i,ik,rho)
!      write(200,*)hpsi
!   eigenvalue 
      eig=(0.0, 0.0)
      do ig=1,n_g_vector(ik)
        eig = eig - conjg(wave_function_c(ig,i,ik))*Hpsi(ig)
!        write(44,*)ig,i,ik,conjg(wave_function_c(ig,i,ik)),Hpsi(ig)
      ENDDO 
!      write(32,*)i,ik,real(eig)
      eigenvalue(i,ik) = real(eig)
      do ig=1,N_G_vector(ik)
        diag=0.5*GplusK(ig,ik)*tpiba2+vl+vnl(ig)
        tau(ig)=(eigenvalue(i,ik)-diag)/emass
      ENDDO 
      do ig=1,n_g_vector(ik)
        s=wave_function_c(ig,i,ik)
        w=exp(tau(ig)*dts)-1.0
!        write(40,*)dts,w,tau(ig),hpsi(ig),eigenvalue(i,ik),wave_function_c(ig,i,ik)
        hpsi(ig)=w*(hpsi(ig)+eigenvalue(i,ik)*wave_function_c(ig,i,ik))/tau(ig)/emass
        wave_function_c(ig,i,ik)=s+hpsi(ig)
      ENDDO 
    ENDDO 
    call gram_schmidt(ik,n_orbitals,n_g_vector(ik))        
  ENDDO 

end  

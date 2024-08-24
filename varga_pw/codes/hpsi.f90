SUBROUTINE H_Psi(i,ik,v)
! H|Psi>
      use fft_data
      use gvector
      use PSEUDOPOTENTIAL
      use PW

      IMPLICIT NONE 

      INTEGER 		:: i, ik,mode
      real*8    	:: v(N_mesh)

      complex*16	:: t1
      real*8            :: rr((L_pp_max+1)**2), fac
      complex*16	:: t((L_pp_max+1)**2),tt
      INTEGER 		:: ir, ig, igp, is, mlmax, ia, m
      INTEGER 		:: i1,i2,i3, l_mesh
      INTEGER 		:: lm, lm_end
      INTEGER 		:: fft_mesh(3)

      complex*8  	:: c_fft(N_mesh)

      fft_mesh = (/N_L(1)+fftinc1,N_L(2),N_L(3)/)
        c_fft(1:N_mesh) = (0.0, 0.0) 
!     wf--> real space
  do ig=1,n_g_vector(ik)
     c_fft(fft_ind(ig,ik))=wave_function_c(ig,i,ik)
  ENDDO 
  call fft_fast(c_fft,N_L(1),N_L(2),N_L(3),ind3f1(1,ik),ind3f2(1,ik),ind2f1(ik),ind2f2(ik),.true.)

!    multiply local potential and wf in real space         
  c_fft(1:N_mesh) = v(1:N_mesh)*c_fft(1:N_mesh)

 

!  backtransform product of local potential and wavefunction 
  call fft_fast(c_fft,N_L(1),N_L(2),N_L(3),ind3f1(1,ik),ind3f2(1,ik),ind2f1(ik),ind2f2(ik),.false.)
! kinetic energy and potential
  do ig=1,n_g_vector(ik)
     t1=(0.5*gplusk(ig,ik)*tpiba2)*wave_function_c(ig,i,ik)
     HPsi(ig) = -(t1 + c_fft(fft_ind(ig,ik))) 
!     write(55,*)ig,i,t1,hpsi(ig)
  ENDDO 
! pseudopotential

      do is=1,n_species
        lm_end=l_max(is)**2+1-2*l_loc(is)
        rr(1:lm_end) = wnl(is,1:lm_end)
        do ia = 1,n_atom(is)
          do lm=1,lm_end
            t(lm)=-fnl(i,ik,is,ia,lm)*rr(lm)
          ENDDO 
          do ig=1,N_G_vector(ik)
            tt=(0.,0.)
            do lm=1,lm_end
              tt=tt+t(lm)*pkg_a(lm,ig,is,ik)
            ENDDO 
            HPsi(ig) = HPsi(ig) + eigr(ig,ia,is,ik)*tt
          ENDDO 
         ENDDO  
      ENDDO 

end SUBROUTINE H_Psi

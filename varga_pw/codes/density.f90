
      SUBROUTINE calculate_density(electron_density)

      use fft_data
      USE GVECTOR
      USE PW
      USE PSEUDOPOTENTIAL
      IMPLICIT NONE 


      INTEGER       nnr1b, nkpt_run,  ik, ig, igp,ii
      INTEGER       r1, r2, r3,  i_store
      INTEGER       i1, i2, i3, i, j, irot
      real*8       rsum1,dsum,rsum,xkin
      real*8       wkptik, wdotf, wdot2fb, ss, sk, scg
      real*8       electron_density(N_L(1)+fftinc1,N_L(2),N_L(3))
      complex*8    c_t,ssu
      real*4       tt, wdotfbo

      complex*8    c_fft(N_L(1)+fftinc1,N_L(2),N_L(3))

        
      rsum=0.0
      xkin=0.0
      electron_density=0.d0
      do ik=N_k_points,1,-1
            wkptik=w_k_point(ik)
         do i=N_orbitals,1,-1
            wdotf=wkptik*occupation(i,ik)
            wdotfbo=wdotf/volume
            wdotf=wdotf*n_sym
            wdot2fb=wdotf*0.5
            ss=0.0
            sk=0.0
            c_fft=(0.0,0.0)

            do ig=1,N_G_vector(ik)
               c_t=wave_function_c(ig,i,ik)               
               igp = G_index(ig,ik)
               c_fft(G_vector(1,igp),G_vector(2,igp),G_vector(3,igp))=c_t                              
               scg=real(c_t)**2+aimag(c_t)**2
               ss=ss+scg
               sk=sk+gplusk(ig,ik)*scg               
            ENDDO 

            rsum=rsum+ss*wdotf
            xkin=xkin+sk*wdot2fb
            
            call fft_fast(c_fft,N_L(1),N_L(2),N_L(3),ind3f1(1,ik),ind3f2(1,ik),ind2f1(ik),ind2f2(ik),.true.)
                      ssu=(0.,0.)
            ii=0
            do i3=1,N_L(3)
               do i2=1,N_L(2)
                  do i1=1,N_L(1)
                     c_t=c_fft(i1,i2,i3)
                     tt=(real(c_t)**2+aimag(c_t)**2)*wdotfbo
                     electron_density(i1,i2,i3) = electron_density(i1,i2,i3)+tt
                     ssu=ssu+tt
                     ii=ii+1
                  ENDDO 
                ENDDO 
            ENDDO 
         ENDDO 


      ENDDO 

     
      E_kinetic=xkin*tpiba2
           
      call symrho(electron_density)

    

      electron_density=dens_mix*real(store_density)+(1.0-dens_mix)*electron_density
      store_density=electron_density


       end


      SUBROUTINE symrho(electron_density)

      use fft_data
      use gvector
      IMPLICIT NONE 


      INTEGER       r1, r2, r3, sr1, sr2, sr3
      INTEGER       i1, i2, i3, irot
      INTEGER       isr1(48), isr2(48), isr3(48)
      real*8       sum
      logical      check(N_L(1),N_L(2),n_L(3))
      real*8       electron_density(N_L(1)+fftinc1,N_L(2),N_L(3))
  
        check=.false.
     
         do r1=0,N_L(1)-1
            do r2=0,N_L(2)-1
               do r3=0,n_L(3)-1
                  if (.not.check(r1+1,r2+1,r3+1)) then
                     sum=0.0
                     do irot=1,N_sym

                        sr1 =       sym_mat(1,1,irot)*r1
                        sr1 = sr1 + sym_mat(2,1,irot)*r2
                        sr1 = sr1 + sym_mat(3,1,irot)*r3
                        sr1 = mod(sr1,N_L(1))+1
                        if (sr1.lt.1) sr1 = sr1+N_L(1)
!     
                        sr2 =       sym_mat(1,2,irot)*r1
                        sr2 = sr2 + sym_mat(2,2,irot)*r2
                        sr2 = sr2 + sym_mat(3,2,irot)*r3
                        sr2 = mod(sr2,N_L(2))+1
                        if (sr2.lt.1) sr2 = sr2+N_L(2)
!     
                        sr3 =       sym_mat(1,3,irot)*r1
                        sr3 = sr3 + sym_mat(2,3,irot)*r2
                        sr3 = sr3 + sym_mat(3,3,irot)*r3
                        sr3 = mod(sr3,N_L(3))+1
                        if (sr3.lt.1) sr3 = sr3+N_L(3)
!     
                        sum = sum + electron_density(sr1,sr2,sr3)
                        write(401,*)sr1,sr2,sr3,electron_density(sr1,sr2,sr3)
                        isr1(irot) = sr1
                        isr2(irot) = sr2
                        isr3(irot) = sr3
                     ENDDO 
          
                     do irot=1,N_sym
                        electron_density(isr1(irot),isr2(irot),isr3(irot))=sum
                        check(isr1(irot),isr2(irot),isr3(irot))=.true.
                     ENDDO 
                  endif
               ENDDO 
            ENDDO 
         ENDDO 

           


      end

      subroutine fft_fast(c_fft,nr11,nr22,nr33,nr3f1,nr3f2,nr2f1,
     &                    nr2f2,t_inv)
c Computes the Fourier Transform 
c    forward for t_inv=.f.
c    inverse for t_inv=.t.
c of a complex   function c_fft. The Fourier Transform is
c returned in c_fft on output (the input c_fft is overwritten).
c
c Use (single precision) IMSL routines for one-dimensional FFT.
c
        implicit none
        integer nr11, nr22, nr33
        integer nr3f1(*),nr3f2(*)
        integer nr2f1,nr2f2
        complex*8 c_fft(nr11+1,nr22,nr33)
        logical t_inv
c
        integer nmax, i1, i2, i3, nmax2, nnmax
        parameter (nmax=200)
        parameter (nmax2=2*nmax, nnmax=4*nmax+15)
        complex*8 c_tmp(nmax)
        real*4 cpy(nmax2), wfftc(nnmax), scale
c
        if(max(nr11,nr22,nr33).gt.nmax) stop ' nmax in fft'
c
        if(t_inv) then
           scale=1.0
           call cffti(nr33,wfftc)
           do i1=1,nr11
              do i2=1,nr3f1(i1)
c     
                 do i3=1,nr33
                    c_tmp(i3)=c_fft(i1,i2,i3)
                 enddo
c     
                 call cfftf(nr33,c_tmp,wfftc)
c     
                 do i3=1,nr33
                    c_fft(i1,i2,i3)=c_tmp(i3)
                 enddo
              enddo

              do i2=nr3f2(i1),nr22
c     
                 do i3=1,nr33
                    c_tmp(i3)=c_fft(i1,i2,i3)
                 enddo
c     
                 call cfftf(nr33,c_tmp,wfftc)
c     
                 do i3=1,nr33
                    c_fft(i1,i2,i3)=c_tmp(i3)
                 enddo
              enddo
           enddo
c     
c     
           call cffti(nr22,wfftc)
           do i1=1,nr2f1
              do i3=1,nr33
c     
                 do i2=1,nr22
                    c_tmp(i2)=c_fft(i1,i2,i3)
                 enddo
c     
                 call cfftf(nr22,c_tmp,wfftc)
c     
                 do i2=1,nr22
                    c_fft(i1,i2,i3)=c_tmp(i2)
                 enddo
              enddo
           enddo

           do i1=nr2f2,nr11
              do i3=1,nr33
c     
                 do i2=1,nr22
                    c_tmp(i2)=c_fft(i1,i2,i3)
                 enddo
c     
                 call cfftf(nr22,c_tmp,wfftc)
c     
                 do i2=1,nr22
                    c_fft(i1,i2,i3)=c_tmp(i2)
                 enddo
              enddo
           enddo
c     
c     
           call cffti(nr11,wfftc)
           do i2=1,nr22
              do i3=1,nr33
c     
                 do i1=1,nr11
                    c_tmp(i1)=c_fft(i1,i2,i3)
                 enddo
c     
                 call cfftf(nr11,c_tmp,wfftc)
c     
                 do i1=1,nr11
                    c_fft(i1,i2,i3)=c_tmp(i1)
                 enddo
              enddo
           enddo
c     
        else
           scale=1.0/real(nr11*nr22*nr33)
c     
           call cffti(nr11,wfftc)
           do i2=1,nr22
              do i3=1,nr33
c     
                 do i1=1,nr11
                    c_tmp(i1)=c_fft(i1,i2,i3)
                 enddo
c     
                 call cfftb(nr11,c_tmp,wfftc)
c     
                 do i1=1,nr11
                    c_fft(i1,i2,i3)=c_tmp(i1)
                 enddo
              enddo
           enddo
c     
           call cffti(nr22,wfftc)
           do i1=1,nr2f1
              do i3=1,nr33
c     
                 do i2=1,nr22
                    c_tmp(i2)=c_fft(i1,i2,i3)
                 enddo
c     
                 call cfftb(nr22,c_tmp,wfftc)
c     
                 do i2=1,nr22
                    c_fft(i1,i2,i3)=c_tmp(i2)
                 enddo
              enddo
           enddo

           do i1=nr2f2,nr11
              do i3=1,nr33
c     
                 do i2=1,nr22
                    c_tmp(i2)=c_fft(i1,i2,i3)
                 enddo
c     
                 call cfftb(nr22,c_tmp,wfftc)
c     
                 do i2=1,nr22
                    c_fft(i1,i2,i3)=c_tmp(i2)
                 enddo
              enddo
           enddo
c     
c     
           call cffti(nr33,wfftc)
           do i1=1,nr11
              do i2=1,nr3f1(i1)
c     
                 do i3=1,nr33
                    c_tmp(i3)=c_fft(i1,i2,i3)
                 enddo
c     
                 call cfftb(nr33,c_tmp,wfftc)
                 
c     
                 do i3=1,nr33
                    c_fft(i1,i2,i3)=c_tmp(i3)
                 enddo
              enddo

              do i2=nr3f2(i1),nr22
c     
                 do i3=1,nr33
                    c_tmp(i3)=c_fft(i1,i2,i3)
                 enddo
c     
                 call cfftb(nr33,c_tmp,wfftc)
                 
c     
                 do i3=1,nr33
                    c_fft(i1,i2,i3)=c_tmp(i3)
                 enddo
              enddo
           enddo
c     
           do i3=1,nr33
              do i2=1,nr22
                 do i1=1,nr11
                    c_fft(i1,i2,i3)=c_fft(i1,i2,i3)*scale
                 enddo
              enddo
           enddo
        endif
c     
c     
        return
        end

! 
! $Id: fft_fftpack.f90,v 1.2 1999/03/10 17:27:08 noe Exp $
! Last pre-cvs change:  JN    5 Feb 97   10:47 am

      subroutine fft(c_fft,nr11,nr22,nr33,t_inv)

! Computes the Fourier Transform 
!    forward for t_inv=.f.
!    inverse for t_inv=.t.
! of a complex   function c_fft. The Fourier Transform is
! returned in c_fft on output (the input c_fft is overwritten).
!
! Use (double precision) fftpack routines for one-dimensional FFT.

   use param
   use fft_data

   implicit none
   integer	     ::  nr11, nr22, nr33
   complex(c8_kind)  ::  c_fft(nr11+1,nr22,nr33)
   logical           ::  t_inv

   integer           ::  nmax, i1, i2, i3, nnmax

   complex (c_kind),  allocatable	:: c_tmp(:) ! c_tmp(nmax)
   real    (r4_kind), allocatable	:: wfftc(:) ! wfftc(nnmax)
   real    (r4_kind)			:: scale

   nmax = max(nr11,nr22,nr33)        
   nnmax=4*nmax+15

   allocate (c_tmp(nmax))
   allocate (wfftc(nnmax))

   if(t_inv) then
       scale=1.0d0
   else
       scale=1.0d0/real(nr11*nr22*nr33)
   endif

   call cffti(nr33,wfftc)
   do i1=1,nr11
      do i2=1,nr22
         do i3=1,nr33
              c_tmp(i3)=c_fft(i1,i2,i3)
         enddo

         if(.not.t_inv) then
                 call cfftb(nr33,c_tmp,wfftc)
         else
                 call cfftf(nr33,c_tmp,wfftc)
         endif

         do i3=1,nr33
                 c_fft(i1,i2,i3)=c_tmp(i3)
         enddo
      enddo
   enddo

   call cffti(nr22,wfftc)
   do i1=1,nr11
      do i3=1,nr33
         do i2=1,nr22
            c_tmp(i2)=c_fft(i1,i2,i3)
         enddo

         if(.not.t_inv) then
                 call cfftb(nr22,c_tmp,wfftc)
         else
                 call cfftf(nr22,c_tmp,wfftc)
         endif

         do i2=1,nr22
                 c_fft(i1,i2,i3)=c_tmp(i2)
         enddo
      enddo
  enddo

  call cffti(nr11,wfftc)
  do i2=1,nr22
     do i3=1,nr33
        do i1=1,nr11
           c_tmp(i1)=c_fft(i1,i2,i3)
        enddo

        if(.not.t_inv) then
                 call cfftb(nr11,c_tmp,wfftc)
        else
                 call cfftf(nr11,c_tmp,wfftc)
        endif

        do i1=1,nr11
           c_fft(i1,i2,i3)=c_tmp(i1)
        enddo
     enddo
  enddo

  if(.not.t_inv) then
     do i3=1,nr33
        do i2=1,nr22
           do i1=1,nr11
              c_fft(i1,i2,i3)=c_fft(i1,i2,i3)*scale
           enddo
        enddo
     enddo
  endif

  deallocate (c_tmp, wfftc)

  return
end

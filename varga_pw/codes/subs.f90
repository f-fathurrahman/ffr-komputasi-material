
!     Spherical harmonics

      FUNCTION Ylm(x,y,z,lm)

      real(8),intent(IN) :: x,y,z
      INTEGER ,intent(IN) :: lm
      real(8) :: rrrr
      real(8) :: Ylm

      select case( lm )

      case(1)  ; Ylm=1.d0
      case(2)  ; Ylm=y
      case(3)  ; Ylm=z
      case(4)  ; Ylm=x
      case(5)  ; Ylm=sqrt(3.d0)*x*y                     ! lm=5  (2 -2)
      case(6)  ; Ylm=sqrt(3.d0)*y*z                     ! lm=6  (2 -1)
      case(7)  ; Ylm=(2*z*z-x*x-y*y)/2.d0               ! lm=7  (2 0)
      case(8)  ; Ylm=sqrt(3.d0)*x*z                     ! lm=8  (2 1)
      case(9)  ; Ylm=sqrt(3.d0/4.d0)*(x*x-y*y)          ! lm=9  (2 2)
      case(10) ; Ylm=sqrt(5.d0/8.d0)*y*(3*x*x-y*y)      ! lm=10 (3 -3)
      case(11) ; Ylm=sqrt(15.d0)*x*y*z                  ! lm=11 (3 -2)
      case(12) ; Ylm=sqrt(3.d0/8.d0)*y*(4*z*z-x*x-y*y)  ! lm=12 (3 -1)
      case(13) ; Ylm=z*(2*z*z-3*x*x-3*y*y)/2.d0         ! lm=13 (3 0)
      case(14) ; Ylm=sqrt(3.d0/8.d0)*x*(4*z*z-x*x-y*y)  ! lm=14 (3 1)
      case(15) ; Ylm=sqrt(15.d0/4.d0)*z*(x*x-y*y)       ! lm=15 (3 2)
      case(16) ; Ylm=sqrt(5.d0/8.d0)*x*(x*x-3*y*y)      ! lm=16 (3 3)

      end select

      END FUNCTION Ylm

       SUBROUTINE lubksb_c(a,n,np,indx,b)
       INTEGER n,np,indx(n)
       complex*16 a(np,np),b(n)
       INTEGER i,ii,j,ll
       complex*16 sum
       ii=0
       do 12 i=1,n
         ll=indx(i)
         sum=b(ll)
         b(ll)=b(i)
         if (ii.ne.0)then
           do 11 j=ii,i-1
             sum=sum-a(i,j)*b(j)
11         continue
         else if (sum.ne.(0.d0,0.d0)) then
           ii=i
         endif
         b(i)=sum
12     continue
       do 14 i=n,1,-1
         sum=b(i)
         do 13 j=i+1,n
           sum=sum-a(i,j)*b(j)
13       continue
         b(i)=sum/a(i,i)
14     continue
       return
       END SUBROUTINE lubksb_c

       SUBROUTINE ludcmp_c(a,n,np,indx,d)
       INTEGER n,np,indx(n),Ndim
       REAL*8 d,TINY
       COMPLEX*16 a(np,np),sum,du
       PARAMETER (Ndim=5000,TINY=1.0d-20)
       INTEGER i,imax,j,k
       REAL*8 aamax,vv(ndim),dum
       d=1.d0
       do 12 i=1,n
         aamax=0.d0
         do 11 j=1,n
           if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
 11      continue
         if (aamax.eq.0.d0) pause 'singular matrix in ludcmp'
         vv(i)=1.d0/aamax
 12    continue
       do 19 j=1,n
         do 14 i=1,j-1
           sum=a(i,j)
           do 13 k=1,i-1
             sum=sum-a(i,k)*a(k,j)
 13        continue
           a(i,j)=sum
 14      continue
         aamax=0.d0
         do 16 i=j,n
           sum=a(i,j)
           do 15 k=1,j-1
             sum=sum-a(i,k)*a(k,j)
 15        continue
           a(i,j)=sum
           dum=vv(i)*abs(sum)
           if (dum.ge.aamax) then
             imax=i
             aamax=dum
           endif
 16      continue
         if (j.ne.imax)then
           do 17 k=1,n
             du=a(imax,k)
             a(imax,k)=a(j,k)
             a(j,k)=du
 17        continue
           d=-d
           vv(imax)=vv(j)
         endif
         indx(j)=imax
         if(a(j,j).eq.(0.d0,0.d0)) a(j,j)=TINY
         if(j.ne.n)then
           du=1.d0/a(j,j)
           do 18 i=j+1,n
             a(i,j)=a(i,j)*du
 18        continue
         endif
 19    continue
       return
       END SUBROUTINE ludcmp_c

       SUBROUTINE lubksb_r(a,n,np,indx,b)
       INTEGER n,np,indx(n)
       real*8 a(np,np),b(n)
       INTEGER i,ii,j,ll
       real*8 sum
       ii=0
       do 12 i=1,n
         ll=indx(i)
         sum=b(ll)
         b(ll)=b(i)
         if (ii.ne.0)then
           do 11 j=ii,i-1
             sum=sum-a(i,j)*b(j)
11         continue
         else if (sum.ne.0.d0) then
           ii=i
         endif
         b(i)=sum
12     continue
       do 14 i=n,1,-1
         sum=b(i)
         do 13 j=i+1,n
           sum=sum-a(i,j)*b(j)
13       continue
         b(i)=sum/a(i,i)
14     continue
       return
       END SUBROUTINE lubksb_r
      

       SUBROUTINE ludcmp_r(a,n,np,indx,d)
       INTEGER n,np,indx(n),Ndim
       REAL*8 d,TINY
       real*8 a(np,np),sum,du
       PARAMETER (Ndim=5000,TINY=1.0d-20)
       INTEGER i,imax,j,k
       REAL*8 aamax,vv(ndim),dum
       d=1.d0
       do 12 i=1,n
         aamax=0.d0
         do 11 j=1,n
           if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
 11      continue
         if (aamax.eq.0.d0) pause 'singular matrix in ludcmp'
         vv(i)=1.d0/aamax
 12    continue
       do 19 j=1,n
         do 14 i=1,j-1
           sum=a(i,j)
           do 13 k=1,i-1
             sum=sum-a(i,k)*a(k,j)
 13        continue
           a(i,j)=sum
 14      continue
         aamax=0.d0
         do 16 i=j,n
           sum=a(i,j)
           do 15 k=1,j-1
             sum=sum-a(i,k)*a(k,j)
 15        continue
           a(i,j)=sum
           dum=vv(i)*abs(sum)
           if (dum.ge.aamax) then
             imax=i
             aamax=dum
           endif
 16      continue
         if (j.ne.imax)then
           do 17 k=1,n
             du=a(imax,k)
             a(imax,k)=a(j,k)
             a(j,k)=du
 17        continue
           d=-d
           vv(imax)=vv(j)
         endif
         indx(j)=imax
         if(a(j,j).eq.0.d0) a(j,j)=TINY
         if(j.ne.n)then
           du=1.d0/a(j,j)
           do 18 i=j+1,n
             a(i,j)=a(i,j)*du
 18        continue
         endif
 19    continue
       return
       END SUBROUTINE ludcmp_r

SUBROUTINE inv_c(a,n,ai)
IMPLICIT NONE 
  INTEGER       :: n,i
  complex*16   :: a(n,n),ai(n,n)
  INTEGER       :: indx(n)
  real*8       :: d

  ai=0.d0
  do i=1,n
    ai(i,i)=1.d0
  ENDDO 
  call ludcmp_c(a,n,n,indx,d)
  do i=1,n
    call lubksb_c(a,n,n,indx,ai(1,i))
  ENDDO 

end SUBROUTINE inv_c


SUBROUTINE inv_r(a,n,ai)
IMPLICIT NONE 
  INTEGER       :: n,i
  real*8       :: a(n,n),ai(n,n)
  INTEGER       :: indx(n)
  real*8       :: d

  ai=0.d0
  do i=1,n
    ai(i,i)=1.d0
  ENDDO 
  call ludcmp_r(a,n,n,indx,d)
  do i=1,n
    call lubksb_r(a,n,n,indx,ai(1,i))
  ENDDO 

end SUBROUTINE inv_r




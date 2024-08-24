PROGRAM per_fd1d

  USE linalg
  implicit none
  INTEGER, PARAMETER :: n = 401
  REAL(8), PARAMETER :: h2m = 0.5d0
  REAL(8) :: h(n,n),x(n),e(n),v(n,n),vv(n)
  REAL(8) :: t,p,ss,hh,omega,a,b,w,v0,pi,l
  INTEGER :: i,j

  a=0.d0
  b=5.d0
  l=b-a
  v0=1.5d0
  pi=4.d0*atan(1.d0)
  hh=(b-a)/(n-1)
  do i=1,n
    x(i)=(i-1)*hh
  end do
  w=h2m/hh**2
  h=0.d0
  do i=1,n
    h(i,i)=2.d0*w
    if(i.ne.n) then
      h(i,i+1)=-w
      h(i+1,i)=-w
    endif
  end do
  h(1,n)=-w
  h(n,1)=-w

  vv=0.d0
  do i=1,n
! Matheiu
!    vv(i)=v0*(1.d0+cos(2.d0*pi*x(i)/l))
!    h(i,i)=h(i,i)+vv(i)

!  Kronig-Penney    
    if(x(i).lt.l/4.d0) h(i,i)=h(i,i)+v0
    if(x(i).gt.3.d0*l/4.d0) h(i,i)=h(i,i)+v0
    if(x(i).lt.l/4.d0.or.x(i).gt.3.d0*l/4.d0) vv(i)=v0
  end do  

  call diag(h,n,e,v)

  do i=1,n
    write(6,*)e(i)
  end do

  do i=1,5
    do j=1,n
      write(10,*)x(j),4*v(j,i)+i*1.d0
    end do
    do j=1,n
      write(10,*)l+x(j),4*v(j,i)+i*1.d0
    end do
    write(10,*)
  end do
  do j=1,n
    write(10,*)x(j),vv(j)
  end do
  do j=1,n
    write(10,*)l+x(j),vv(j)
  end do
  
END PROGRAM 

  




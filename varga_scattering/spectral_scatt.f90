!Use of spectral projection to obtain eigenstates of a Gaussian function from a step potential 

USE linalg
implicit none
  complex*16,parameter                    :: zi=(0.d0,1.d0)
!Constants
  real*8,parameter                        :: h2m=0.5d0,pi=3.1415926,alpha_w=2.d0
!Left and right boundary values
  real*8                                  :: vpl,vpr,kl,kr,al,l0,v0,k1,k2
!Values to be used in do-loops
  integer                                 :: N,i,j,k,l,Nd,i1,i2,en,m
!Spectral operator/green's function/hamiltonians
  complex*16,dimension(:,:),allocatable   :: am,gm,hm,hm0,gmt
!Wave functions/gaussian function
  complex*16,dimension(:),allocatable     :: capl,capr,wf,wfa,wfb,psi_a,psi_r
!Various values
  real*8,dimension(:),allocatable         :: pot,xgr
  real*8,dimension(:),allocatable         :: eva
  real*8,dimension(:,:),allocatable       :: hmat,eve
  real*8                                  :: h,t,tt,a,b,c1,c2,x,cf,ee,dd,yl,yr,dr,c,xl,xr,k0,ke,x0
  complex*16                              :: e,tc,su,af,bf,s11,s21
  
  c2=10.d0
  Nd=600
  a=-30.d0
  b=+30.d0
  c1=-c2
  dr=b-c2+0.01d0
  c=2.62
  h=(b-a)/(Nd+1)
  vpl=-1.d0
  vpr=1.d0
  t=h2m/h**2
  k0=0.1d0
  x0=0.d0
  
  allocate(hm(Nd,Nd),hm0(Nd,Nd),am(Nd,Nd),capl(Nd),capr(Nd),gm(Nd,Nd),pot(Nd), &
&   hmat(Nd,Nd),eva(Nd),Eve(Nd,Nd),wf(Nd),wfa(nd),wfb(nd),xgr(nd),gmt(nd,nd),psi_a(nd),psi_r(nd))

  hm=(0.d0,0.d0)
  hm0=(0.d0,0.d0)
  do i=1,Nd
    x=a+(i-1)*h
    xgr(i)=x
!Setting up gaussian function
    wfa(i)=1.d0/(pi*alpha_w**2)**0.25d0*exp(zi*k0*x-0.5d0*(x-x0)**2/alpha_w**2)
!Setting up boundaries
    capl(i)=0.d0
    capr(i)=0.d0
    xl=c*(x-c1)/dr
    xr=c*(x-c2)/dr
    yr=4.d0/(c-xr)**2+4.d0/(c+xr)**2-8.d0/c**2   
    yl=4.d0/(c-xl)**2+4.d0/(c+xl)**2-8.d0/c**2   
    if(x.lt.c1) capl(i)=-zi*h2m*(2*pi/dr)**2*yl
    if(x.gt.c2) capr(i)=-zi*h2m*(2*pi/dr)**2*yr
    write(20,*)x,imag(capl(i)),imag(capr(i))
    hm0(i,i)=205.d0/72.d0/h**2*h2m
    if(i.ne.1) then
      hm0(i,i-1)=-8.d0/5.d0/h**2*h2m
      hm0(i-1,i)=-8.d0/5.d0/h**2*h2m
    endif
    if(i.gt.2) then
      hm0(i,i-2)=1.d0/5.d0/h**2*h2m
      hm0(i-2,i)=1.d0/5.d0/h**2*h2m
    endif
    if(i.gt.3) then
      hm0(i,i-3)=-8.d0/315.d0/h**2*h2m
      hm0(i-3,i)=-8.d0/315.d0/h**2*h2m
    endif
    if(i.gt.4) then
      hm0(i,i-4)=1.d0/560.d0/h**2*h2m
      hm0(i-4,i)=1.d0/560.d0/h**2*h2m
    endif
  end do
!Setting up step potential
  do i=1,Nd
    pot(i)=0.d0
    if(i.le.Nd/2) then
      pot(i)=vpl
    else
      pot(i)=vpr
    endif       
    write(20,*)xgr(i),pot(i)
  end do

  do i=1,Nd
    hm0(i,i)=hm0(i,i)+pot(i)+capl(i)/2.d0+capr(i)/2.d0
  end do

!Using only 1 energy (E=1.5) 
do en=5,5
  e=1.d0+en*0.1d0+0.0001d0*zi
   hm=hm0

   am=-hm
   do i=1,Nd
     am(i,i)=e-hm(i,i)
   end do
 
   call inv(am,nd,gm)
   do i=1,nd
     do j=1,nd
       gmt(j,i)=zi*(Conjg(gm(i,j))-gm(j,i))
     end do
   end do
   
   wf=matmul(gmt,wfa)

   do i=1,Nd-1
     x=a+(i-1)*h
     write(40,*)x,real(wf(i)),imag(wf(i))
     write(41,*)x,imag(wf(i))
   end do

end do   


  end

USE LINALG
implicit none
  real*8,parameter                        :: h2m=0.5d0,pi=3.1415926,p=5.d0
  complex*16,parameter                    :: zi=(0.d0,1.d0)
  integer                                 :: i,j,Nd,ie,kk
  real*8,dimension(:),allocatable         :: pot
  real*8,dimension(:),allocatable         :: eva
  real*8,dimension(:,:),allocatable       :: hmat,eve,hm
  real*8                                  :: h,x,a,b,vpl,vpr,t
  real*8                                  :: ka,kb,Energy,aa,bb,cc
  real*8                                  :: fa,fb,Raa,Rab,Rbb,fx,Rax,Rbx
  complex*16                              :: a11,a12,a21,a22,b1,b2,d,tca,rca,ea,eb,w,w1,w2,ww,xw
 
  Nd=500
  h=0.01d0
  vpl=0.d0
  vpr=0.d0
  t=h2m/h**2
  a=0.d0
  b=p
  allocate(hm(Nd,Nd),pot(Nd),hmat(Nd,Nd),eve(Nd,Nd),eva(Nd))
  
  hm=(0.d0,0.d0)
  do i=1,Nd
    hm(i,i)=2.d0*t
    if(i.gt.1) then
      hm(i,i-1)=-t
      hm(i-1,i)=-t
    endif
  end do
  hm(1,1)=t
  hm(Nd,Nd)=t

  do i=1,Nd
    x=i*h
    hm(i,i)=hm(i,i)+pot(i)+1.5d0*(1.d0+cos(2.d0*pi*x/P))
  end do
  
  hmat=hm
  call diag(hmat,Nd,eva,eve)


do kk=0,800
  Energy=kk*0.01d0
  ka=sqrt((energy-vpl)/h2m)
  kb=sqrt((energy-vpr)/h2m)

  Raa=0.d0
  Rab=0.d0
  Rbb=0.d0

  do i=1,Nd
    fa=eve(1,i)/sqrt(h)
    fb=eve(Nd,i)/sqrt(h)
    Raa=Raa-0.5d0*fa*fa/(Energy-eva(i))
    Rab=Rab-0.5d0*fa*fb/(Energy-eva(i))
    Rbb=Rbb-0.5d0*fb*fb/(Energy-eva(i))
  end do

  ea=exp(zi*ka*a)
  eb=exp(zi*kb*b)
  xw=exp(zi*ka*(b-a))

  aa=Rab
  bb=-Raa-Rbb
  cc=Rab

  ww=bb**2-4.d0*aa*cc
  w=sqrt(ww)
  w1=(-bb+w)/aa*0.5d0
  w2=(-bb-w)/aa*0.5d0
  xw=w1
!  write(6,*)xw
  write(6,*)'----------',xw
  write(6,*)energy,real(log(xw)/(b-a)/zi)
  write(10,*)abs(real(log(xw)/(b-a)/zi))*P,energy
  write(11,*)-abs(imag(log(xw)/(b-a)/zi))*2.d0,energy
!  xw=w2
!  write(6,*)xw
!  write(6,*)'----------'
!  write(6,*)energy,real(log(xw)/(b-a)/zi)
!  write(11,*)abs(real(log(xw)/(b-a)/zi)),energy
end do
end






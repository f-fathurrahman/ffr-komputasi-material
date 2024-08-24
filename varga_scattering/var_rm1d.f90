implicit none                  
real*8,parameter                  :: h2m=0.5d0,a=1.d0,b=5.d0,V0=2.5d0,x0=3.d0,&
&                                    mu=0.2d0,d=1.d0
complex*16,parameter              :: zi=(0.d0,1.d0)
integer                           :: i,j,k,n,m,Ndim
real*8                            :: st,sv,f1,f2,df1,df2,sh,so,xxa
real*8, external                  :: potential

real*8,dimension(:),allocatable   :: pot,evec

integer                           :: ifail,isec,ii,ien
integer                           :: Nstep

real*8,dimension(:,:),allocatable :: z,T,S,V,da,db,hm,om,am,ai,hm0
real*8,dimension(:),allocatable   :: alfi,alfr,beta,iter,kappa,xx,xv

real*8                             :: w,lam_a,lam_b,eps1,k_L,k_R,v_L,v_R,energy,h,x,ss
complex*16                         :: A_L1,A_R1,B_L1,B_R1,A_L2,A_R2,B_L2,B_R2,det,rc,tc
logical                            :: matz
real*8,external                    :: x02ajf


Nstep=100000

open(1,file='basis.inp')
read(1,*)Ndim

allocate(T(Ndim,Ndim),S(Ndim,Ndim),V(Ndim,Ndim),pot(0:Nstep),da(Ndim,Ndim),db(Ndim,Ndim),kappa(ndim))
allocate(hm(Ndim,Ndim),om(Ndim,Ndim),evec(Ndim),xx(Ndim),ai(Ndim,ndim),am(Ndim,Ndim))

do i=1,Ndim
  read(1,*)kappa(i)
end do
close(1)

h=(b-a)/Nstep
do k=0,Nstep
  x=a+k*h      
  w=(x-x0)/d
  pot(k)=-V0*sinh(w)**2/cosh(w-mu)**2
  write(10,*)x,pot(k),potential(x)
end do
  x=-10.d0      
  w=(x-x0)/d
  v_L=-V0*sinh(w)**2/cosh(w-mu)**2
  x=+10.d0      
  w=(x-x0)/d
  v_R=-V0*sinh(w)**2/cosh(w-mu)**2

do i=1,Ndim
  xx(i)=cos(kappa(i)*a)
  do j=1,Ndim
    sv=0.d0
    st=0.d0
    ss=0.d0
    do k=0,Nstep
      x=a+k*h      
      f1=cos(kappa(i)*x)
      f2=cos(kappa(j)*x)
      df1=-kappa(i)*sin(kappa(i)*x)
      df2=-kappa(j)*sin(kappa(j)*x)
      ss=ss+h*f1*f2
      sv=sv+h*f1*f2*pot(k)
      st=st+h*df1*df2
    end do
    T(i,j)=st
    S(i,j)=ss
    V(i,j)=sv
    f1=cos(kappa(i)*a)
    f2=cos(kappa(j)*a)
    da(i,j)=f1*f2
    f1=cos(kappa(i)*b)
    f2=cos(kappa(j)*b)
    db(i,j)=f1*f2
  end do
end do

Nstep=1000
h=(b-a)/Nstep
do ien=0,0
  energy=-0.3d0
  k_L=sqrt((energy-v_L)/h2m)
  k_R=sqrt((energy-v_R)/h2m)
  do isec=1,2
  if(isec.eq.1) lam_b=+2.d0
  if(isec.eq.2) lam_b=-2.d0
  do i=1,Ndim
    do j=1,Ndim
      hm(i,j)=-T(i,j)+energy*S(i,j)/h2m-V(i,j)/h2m+lam_b*db(i,j)
      om(i,j)=da(i,j)
    end do
  end do
  am=hm

  call inv_r(am,ndim,ai)
  xxa=0.d0
  do i=1,Ndim
    do j=1,Ndim
      xxa=xxa+xx(i)*ai(i,j)*xx(j)
    end do
  end do
  write(6,*)'sol1',1.d0/xxa

  am=hm-1.d0/xxa*om
  call inv_r(am,ndim,ai)
  evec=0.d0
  evec(1)=1.d0
  evec=matmul(ai,evec)
  xxa=sqrt(sum(evec(:)**2))
  write(6,*)'norm',xxa
  evec=evec/xxa
  xxa=sqrt(sum(evec(:)**2))
  write(6,*)'norm',xxa
  
  if(isec.eq.1)  call wave_function_match(a,b,Nstep,Ndim,kappa,k_L,k_R,evec,A_L1,B_L1,A_R1,B_R1)
  if(isec.eq.2)  call wave_function_match(a,b,Nstep,Ndim,kappa,k_L,k_R,evec,A_L2,B_L2,A_R2,B_R2)

end do
det=A_L2*B_R1-B_R2*A_L1
rc=(B_L2*B_R1-B_R2*B_L1)/det
tc=(A_R2*B_R1-B_R2*A_R1)/det

write(6,*)'R',real(Conjg(rc)*rc)
write(6,*)'T',real(Conjg(tc)*tc*k_R/k_L)
write(20,*)energy,real(Conjg(rc)*rc),real(Conjg(tc)*tc*k_R/k_L)
write(6,*)'R+T',real(Conjg(rc)*rc+Conjg(tc)*tc*k_R/k_L)

end do

end


function potential(x)
  real*8,parameter            :: V0=2.5d0,x0=3.d0,mu=0.2d0,d=1.d0
  real*8                      :: x,potential,w
  w=(x-x0)/d
  potential=-V0*sinh(w)**2/cosh(w-mu)**2
end


subroutine wave_function_match(a,b,Nstep,Ndim,kappa,ka,kb,z,a1,b1,a2,b2)
  implicit none
  integer                       :: Nstep,Ndim
  double precision              :: wf(0:Nstep),kappa(Ndim),z(Ndim)
  double precision              :: a,b,h,ka,kb
  complex*16                    :: a1,a2,b1,b2
  double precision              :: fa1,fb1,dfa1,dfb1,alp1,bet1,alp2,bet2,x,ss,f
  complex*16,parameter          :: zi=(0.d0,1.d0)
  integer                       :: i,j,k
real*8, external                  :: potential

  h=(b-a)/Nstep
  do k=0,Nstep
    x=a+k*h      
    ss=0.d0
    do j=1,Ndim
      f=cos(kappa(j)*x)
      ss=ss+f*z(j)  
    end do
    wf(k)=ss
  end do

  fa1=wf(0)
  dfa1=(wf(1)-wf(0))/h
  alp1=fa1*cos(ka*a)-dfa1*sin(ka*a)/ka
  bet1=fa1*sin(ka*a)+dfa1*cos(ka*a)/ka        
  a1=0.5d0*(alp1-zi*bet1)
  b1=0.5d0*(alp1+zi*bet1)
  fb1=wf(Nstep)
  dfb1=(wf(Nstep)-wf(Nstep-1))/h
  alp2=fb1*cos(kb*b)-dfb1*sin(kb*b)/kb
  bet2=fb1*sin(kb*b)+dfb1*cos(kb*b)/kb
  a2=0.5d0*(alp2-zi*bet2)
  b2=0.5d0*(alp2+zi*bet2)

! wave function
  do i=-500,1100
    x=0.01d0*i
    if(x.lt.a) ss=alp1*cos(ka*x)+bet1*sin(ka*x)            
    if(x.gt.b) ss=alp2*cos(kb*x)+bet2*sin(kb*x)            
    if(x.ge.a.and.x.le.b) then
      ss=0.d0
      do j=1,Ndim
        f=cos(kappa(j)*x)
        ss=ss+f*z(j)  
      end do
    endif
    write(12,*)x,ss
    write(13,*)x,potential(x)
  end do
end subroutine wave_function_match




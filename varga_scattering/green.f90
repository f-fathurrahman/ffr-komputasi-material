implicit none
  real*8,parameter                      :: h2m=0.5d0,h=0.006d0,pi=3.141592653589793d0
  complex*16,parameter                  :: zi=(0.d0,1.d0)
  integer                               :: i,j,k,l,n,m,i1,i2,energy
  real*8,dimension(:),allocatable       :: xr
  real*8,dimension(:,:),allocatable     :: hlc,hrc,hcc,hl00,hl10,hr00,hr10,cm
  complex*16,dimension(:,:),allocatable :: GCC_R,GCC_A,SigmaL,SigmaR,GammaL,GammaR,T,GL00,GR00
  complex*16                            :: e
  real*8                                :: ee,ww1,ww2,x,tc,v0,l0,al,vl,vr,tk,tcc
  complex*16                            :: kl,kr,wf,su1,su2,ds,ts


  tk=h2m/h**2
  vl=-1.d0
  vr=1.d0
  n=20
  m=20

  allocate(GL00(m,m),GR00(m,m))
  allocate(hlc(m,n),hrc(m,n),hcc(n,n),hl00(m,m),hl10(m,m),hr00(m,m),hr10(m,m),cm(n,n))
  allocate(GammaR(n,n),GammaL(n,n),SigmaR(n,n),SigmaL(n,n),GCC_A(n,n),GCC_R(n,n),T(n,n))

    
  hlc=0.d0
  hrc=0.d0
  hlc(m,1)=-tk
  hrc(1,n)=-tk
  hcc=0.d0
  do i=1,n
    if(i.ne.1) hcc(i,i-1)=-tk
    if(i.ne.1) hcc(i-1,i)=-tk
  end do            
  do i=1,n/2
    hcc(i,i)=2.d0*tk+vl
  end do
  do i=n/2+1,n
    hcc(i,i)=2.d0*tk+vr
  end do

  hl00=0.d0
  do i=1,m
    hl00(i,i)=2.d0*tk+vl
    if(i.ne.1) hl00(i,i-1)=-tk
    if(i.ne.1) hl00(i-1,i)=-tk
  end do
  hl10=0.d0
  hl10(m,1)=-tk
    hr00=0.d0
    do i=1,m
      hr00(i,i)=2.d0*tk+vr
      if(i.ne.1) hr00(i,i-1)=-tk
      if(i.ne.1) hr00(i-1,i)=-tk
    end do
    hr10=0.d0
    hr10(1,m)=-tk
    cm=0.d0
    do i=1,n
      cm(i,i)=1.d0
    end do

open(1,file="t_vs_e.dat")
do energy=1,51

  e=1.d0+(energy-1)*0.1d0+zi*0.0001d0
  call green(e,hcc,hlc,hrc,hl00,hl10,hr00,hr10,n,m,GCC_R,SigmaL,SigmaR,GL00,GR00)

  write(10,*)real(e-2.d0*tk)/tk,real(gl00(m,m)),imag(gl00(m,m))
  su1=(0.d0,0.d0)
  do i=1,m
    su1=su1+gl00(i,i)
  end do
  write(12,*)real(e-2.d0*tk)/tk,real(su1),imag(su1)
  



  do i=1,n
    do j=1,n        
      GammaL(i,j)=zi*(SigmaL(i,j)-Conjg(SigmaL(j,i)))
      GammaR(i,j)=zi*(SigmaR(i,j)-Conjg(SigmaR(j,i)))
      GCC_A(i,j)=Conjg(GCC_R(j,i))
    end do
  end do

    T=matmul(GCC_R,matmul(GammaL,matmul(GCC_A,GammaR)))
      
    tcc=0.d0  
    do i=1,n
      tcc=tcc+T(i,i)
    enddo
   
    ee=real(e)

    write(1,'(2ES16.8)')ee,tcc
    
    kl=sqrt(ee/h2m-vl/h2m)
    kr=sqrt(ee/h2m-vr/h2m)

    ds=(0.d0,0.d0)
    ts=(0.d0,0.d0)
    do i=1,n
      su1=(0.d0,0.d0)
      su2=(0.d0,0.d0)
      do i1=1,n
        do i2=1,n            
          su1=su1+cm(i,i1)*cm(i,i2)*GCC_R(i1,i2)
          su2=su2+cm(i,i1)*cm(i,i2)*T(i1,i2)
        end do
      end do
      write(11,*)x,real(su1),-imag(su1)
      ds=ds+su1
      ts=ts+su2
    end do
!   potential step
    tc=real(4.d0*kr*kl/(kl+kr)**2)
!   pot well
!    v0=1.d0
!    l0=8.d0
!    al=sqrt((ee+v0)/h2m)
!    tc=1.d0/(1.d0+v0**2*sin(al*l0)**2/(4.d0*ee*(ee+v0)))
    write(6,*)real(e),real(ts),tc
    write(33,*)real(e),real(ts),tc
    write(30,*)real(e),real(ds),-imag(ds)/pi


end do
close(1)
end

subroutine green(e,hcc,hlc,hrc,hl00,hl10,hr00,hr10,n,m,gcc,SigmaL,SigmaR,GL00,GR00)
USE linalg
implicit none
  complex*16,parameter              :: zi=(0.d0,1.d0)
  integer                           :: i,j,k,l,n,m
  real*8                            :: hlc(m,n),hrc(m,n),hcc(n,n)
  real*8                            :: hl00(m,m),hl10(m,m),hr00(m,m),hr10(m,m)
  complex*16                        :: e,SigmaL(n,n),SigmaR(n,n),GCC(n,n), &
&                                      um(n,n),a(n,n),GL00(m,m),GR00(m,m)

  um=(0.d0,0.d0)
  do i=1,n
    um(i,i)=(1.d0,0.d0)
  end do            

  call  green_decimation(e,hl00,hl10,m,gL00)
  call  green_decimation(e,hr00,hr10,m,gR00)
!  call  green_iteration(e,hl00,hl10,m,gL00)
!  call  green_iteration(e,hr00,hr10,m,gR00)



  SigmaL=(0.d0,0.d0)
  SigmaR=(0.d0,0.d0)
  do i=1,n
    do j=1,n        
      do k=1,m
        do l=1,m
          SigmaL(i,j)=SigmaL(i,j)+hlc(k,i)*GL00(k,l)*hlc(l,j)
          SigmaR(i,j)=SigmaR(i,j)+hrc(k,i)*GR00(k,l)*hrc(l,j)
        end do
      end do
    end do
  end do

  do i=1,n
    do j=1,n        
      a(i,j)=e*um(i,j)-hcc(i,j)-SigmaL(i,j)-SigmaR(i,j)
    end do
  end do
  call inv(a,n,gcc)
     
end subroutine green


subroutine green_iteration(e,h00,h10,n,g00)
USE linalg
implicit none
  integer                           :: i,j,k,n
  real*8                            :: h00(n,n),h10(n,n),h10t(n,n)
  real*8                            :: a,b,su,t
  complex*16                        :: e,gs,gb
  complex*16                        :: am(n,n),g00(n,n),ami(n,n),vv(n,n),u0(n,n),um(n,n)


  um=(0.d0,0.d0)
  do i=1,n
    um(i,i)=(1.d0,0.d0)
    do j=1,n
      h10t(i,j)=h10(j,i)
    end do
  end do

  am=e*um-h00
  vv=am
  do k=1,10000
    call inv(vv,n,g00)
    u0=matmul(g00,h10)
    vv=am-matmul(h10t,u0)
  end do
  call inv(vv,n,g00)
    
end subroutine green_iteration


subroutine green_decimation(e,h00,h10,n,g00)
USE linalg
implicit none
  integer                           :: i,j,k,n
  real*8                            :: h00(n,n),h10(n,n),h10t(n,n)
  real*8                            :: a,b,su,t
  complex*16                        :: e,gs,gb
  complex*16                        :: am(n,n),g00(n,n),u0(n,n),v0(n,n),tv(n,n),tu(n,n)
  complex*16                        :: ami(n,n),fv(n,n),fu(n,n),u1(n,n),v1(n,n),g1(n,n)
  complex*16                        :: uu(n,n),vv(n,n),um(n,n),fu1(n,n),fv1(n,n)
  

  um=(0.d0,0.d0)
  do i=1,n
    um(i,i)=(1.d0,0.d0)
    do j=1,n
      h10t(i,j)=h10(j,i)
    end do
  end do
  
  am=e*um-h00
  call inv(am,n,g00)

  u0=matmul(g00,h10)
  v0=matmul(g00,h10t)
  tv=v0
  tu=u0
  fv=v0
  fu=u0
  do j=1,20
    am=um-matmul(u0,v0)-matmul(v0,u0)
    call inv(am,n,ami)
    uu=matmul(u0,u0)
    vv=matmul(v0,v0)
    u1=matmul(ami,uu)
    v1=matmul(ami,vv)
    fu1=matmul(fv,u1)
    fv1=matmul(fu,v1)
    fu=matmul(fu,u1)
    fv=matmul(fv,v1)
    tu=tu+fu1
    tv=tv+fv1
    u0=u1
    v0=v1
  end do
  am=e*um-h00-matmul(h10t,tu)
  call inv(am,n,g00) 
  
end subroutine green_decimation


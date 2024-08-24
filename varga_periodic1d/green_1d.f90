MODULE GREEN_1D 
use linalg
implicit none
! Number of discretization points per cell
integer,parameter             :: N=80
! Number of cells
integer,parameter             :: M=8000
! Length of cell
double precision              :: L_cell=5.d0,b=2.d0
! hbar**2/2m
double precision              :: h2m=3.81d0
double precision,parameter    :: Pi=3.141592653589793d0,epsilon=0.02d0
! energy
integer,parameter             :: N_energy_points=120
double precision              :: E_energy_points(N_energy_points),Energy(N)
double precision,parameter    :: de=0.1d0
!density of states
double precision              :: D0(N_Energy_points)
!wave function
complex*16                    :: phi(0:N+1,N)
complex*16,parameter          :: zi=(0.d0,1.d0)
! potential
double precision              :: V_0(N),V_B(N)
! Green's function
complex*16                    :: Green(N,N,N_Energy_points)


CONTAINS

subroutine calculate_bloch_states(k)
  implicit none
  double precision         :: k
  integer                  :: i,j
  double precision         :: a
  double precision         :: hm(2*N,2*N),evec(2*N,2*N),eval(2*N)

  a=h2m/(L_cell/N)**2
  hm=0.d0
  do i=1,N
    hm(i,i)=2.d0*a+V_0(i); hm(N+i,N+i)=2.d0*a+V_0(i)
    if(i/=N) then
      hm(i,i+1)=-a;       hm(i+1,i)=-a
      hm(N+i,N+i+1)=-a ;  hm(N+i+1,N+i)=-a
    endif
  end do
  hm(1,N)=-a*cos(k*L_cell);     hm(N,1)=-a*cos(k*L_cell)
  hm(N+1,2*N)=-a*cos(k*L_cell); hm(2*N,N+1)=-a*cos(k*L_cell)
  hm(N,N+1)=-a*sin(k*L_cell);   hm(1,2*N)=a*sin(k*L_cell)
  hm(N+1,N)=-a*sin(k*L_cell);   hm(2*N,1)=a*sin(k*L_cell)

  call diag(hm,2*n,2*n,eval,evec)
  write(77,*)k,eval(1),eval(2),eval(3),eval(4)
  
  do i=1,N
    Energy(i)=eval((i-1)*2+1)
    do j=1,N
      phi(j,i)=evec(j,(i-1)*2+1)+zi*evec(N+j,(i-1)*2+1)
    end do
    phi(0,i)=exp(zi*k*L_cell)*phi(N,i)
    phi(N+1,i)=exp(-zi*k*L_cell)*phi(1,i)
     do j=1,N_energy_points
       if(energy(i)<30.d0) then
         D0(j)=D0(j)-imag(1.d0/(E_energy_points(j)+zi*epsilon-energy(i)))/pi
       endif
     end do
  end do
end subroutine calculate_bloch_states

subroutine calculate_Green_spectral
  implicit none
  integer         :: i,j,k,l
  complex*16      :: su,z
  do j=1,N_energy_points
    do l=1,N
      su=1.d0/(E_energy_points(j)+zi*epsilon-energy(l))
      do i=1,N
        z=phi(i,l)*su
        do k=1,N
          Green(k,i,j)=Green(k,i,j)+Conjg(phi(k,l))*z
        end do
      end do
    end do
  end do
end subroutine calculate_Green_spectral

subroutine calculate_Green_wf
  implicit none
  complex*16       :: W
  integer          :: i,j,k,l
  complex*16       :: phi_p(N),deter
  double precision :: d1,d2
  complex*16       :: G0(0:N+1,0:N+1),hp(N),XX(N,N),XI(N,N),r,t,su

  do l=1,N
    i=N/2
    W=h2m/(L_cell/N)**2*(phi(i+1,l)*Conjg(phi(i,l))- &
&     Conjg(phi(i+1,l))*phi(i,l)- &
&     phi(i-1,l)*Conjg(phi(i,l))+ &
&     Conjg(phi(i-1,l))*phi(i,l))/2.d0
    do i=0,N+1
      do j=0,N+1
        if(i.le.j) then
        G0(i,j)=Conjg(phi(i,l))*phi(j,l)/W
        else
        G0(i,j)=phi(i,l)*Conjg(phi(j,l))/W
        endif          
      end do
    end do    
    XX=(0.d0,0.d0)
    do i=1,N
      XX(i,i)=(1.d0,0.d0)
    end do
    do i=1,N
      do j=1,N
        XX(i,j)=XX(i,j)-G0(i,j)*V_B(j)
      end do
    end do
    call inv(XX,N,XI)
    call det(XX,N,deter)
    d2=-imag(log(deter))/pi
    do i=1,N
      phi_p(i)=dot_product(Xi(i,:),phi(1:N,l))  
    end do
    R=(0.d0,0.d0); T=(1.d0,0.d0)
    do i=1,N
      R=R+phi(i,l)*V_B(i)*phi_p(i)/W
      T=T-Conjg(phi(i,l))*V_B(i)*phi_p(i)/W
    end do
    if(energy(l)<10.d0) write(20,*)energy(l),1.d0-real(R*Conjg(R))
    if(energy(l)<10.d0) write(21,*)energy(l),real(T*Conjg(T))
    su=(0.d0,0.d0)
    do i=1,N
      su=su+G0(i,i)
    end do
    if(energy(l).lt.10.d0.and.-imag(su)>0)write(13,*)energy(l),-imag(su)

    if(energy(l).lt.10.d0.and.l.gt.1) then
      write(41,*)energy(l),(d2-d1)/(energy(l)-energy(l-1))
    endif
    d1=d2
  end do

end subroutine calculate_Green_wf


subroutine calculate_perturbed_dos
  implicit none
  complex*16       :: W
  integer          :: i,j,k,l
  complex*16       :: deter
  double precision :: dos(N_energy_points)
  complex*16       :: XX(N,N),XI(N,N)

  do l=1,N_energy_points
    XX=(0.d0,0.d0)
    do i=1,N
      XX(i,i)=(1.d0,0.d0)
    end do
    do i=1,N
      do j=1,N
        XX(i,j)=XX(i,j)-Green(i,j,l)*V_B(j)
      end do
    end do
    call inv(XX,N,XI)
    call det(XX,N,deter)
    dos(l)=-imag(log(deter))/pi
  end do

  do l=1,N_energy_points-1
    write(61,*)e_energy_points(l),(dos(l+1)-dos(l))/de
  end do
  
end subroutine calculate_perturbed_dos

END MODULE GREEN_1D

USE GREEN_1D
implicit none
double precision              :: dx,k,x,a,V_tot(n)
integer                       :: ii,i,j,ii1,ii2,kk,jj,mm,beta
complex*16                    :: su,z,w

do i=1,N_energy_points
  E_energy_points(i)=-2.d0+(i-1)*de
end do

dx=L_cell/N
a=h2m/dx**2
V_tot=-2.d0
do i=1,N
  V_0(i)=0.d0
  x=i*dx
  if(x.ge.(L_cell-b)/2.d0.and.x.le.b+(L_cell-b)/2.d0) V_0(i)=1.d0
end do
V_B=V_tot-V_0

do j=-5,5
  do i=1,N
    x=j*L_cell+i*dx-0.5d0*L_cell
    if(j.ne.0) then
      write(30,*)x,V_0(i)
    else
      write(30,*)x,V_tot(i)
    endif
  end do
end do

! do loop over k vector
do beta=1,M
  write(6,*)beta
  k=-pi/L_cell+(beta-0.5d0)*2.d0*pi/(M*L_cell)
  call calculate_bloch_states(k)
  write(1,'(5d16.8)')k/(2.d0*pi/L_cell),Energy(1),Energy(2),Energy(3),Energy(4)
  call calculate_Green_spectral
  call calculate_Green_wf
end do

  call calculate_perturbed_dos

  do kk=1,N_energy_points
    write(12,*)E_energy_points(kk),D0(kk)
    su=(0.d0,0.d0)
    do i=1,N
      su=su+Green(i,i,kk)
    end do
    write(11,*)E_energy_points(kk),-imag(su)    
  end do

end




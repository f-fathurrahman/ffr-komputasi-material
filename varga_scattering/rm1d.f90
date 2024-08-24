USE linalg
implicit none
real*8,dimension(:),allocatable    :: e,xr,wr
real*8,dimension(:,:),allocatable  :: hm,v,t

real*8,external                    :: potential
real*8                             :: xmu,sc,x,y,dy,energy,Raa,Rab,Rbb
real*8 :: energy_start,energy_step_size,fa,fb,a,b,ka,kb,va,vb,su,h
integer                            :: N,i,j,k,k1,k2,nrot,Nrm,num_energies,ii
real*8,parameter                   :: Pi=3.141592653589793d0,h2m=0.5d0
complex*16,parameter               :: zi=(0.d0,1.d0)
complex*16                         :: a11,a12,a21,a22,b1,b2,x1,x2,d,tc,rc,&
&                                     exa,exb

N=400
a=-5.d0
b=5.d0
num_energies=100
energy_start=-1.65d0
energy_step_size=0.05d0

allocate(wr(N),xr(N),t(N,N),v(N,N),e(N),hm(N,N))

  h=(b-a)/(1.d0*(N+1))
  do i=1,N
    x=a+i*h
    xr(i)=x
    wr(i)=h
  end do
  t=0.d0
  do i=1,N
    t(i,i)=-2.d0/h**2
    if(i.ne.1) then
      t(i,i-1)=1.d0/h**2
      t(i-1,i)=1.d0/h**2
    endif
  end do
!   BC!
  t(1,1)=-1.d0/h**2
  t(N,N)=-1.d0/h**2

  hm=-h2m*t

do i=1,n
  hm(i,i)=hm(i,i)+potential(xr(i))
end do

call diag(hm,n,e,v)

open(1,file="T_vs_E.dat")
  do ii=1,num_energies
    energy=energy_start+(ii-1)*energy_step_size

		va=potential(a)
		vb=potential(b)
		ka=sqrt((energy-va)/h2m)
		kb=sqrt((energy-vb)/h2m)
		Raa=0.d0
		Rab=0.d0
		Rbb=0.d0

		do i=1,N
			write(6,*)e(i)
			fa=v(1,i)
			fb=v(N,i)
			Raa=Raa-0.5d0*fa*fa/(Energy-e(i))/h
			Rab=Rab-0.5d0*fa*fb/(Energy-e(i))/h
			Rbb=Rbb-0.5d0*fb*fb/(Energy-e(i))/h
		end do
		a11=1.d0-zi*ka*Raa
		a12=-zi*kb*Rab
		a21=-zi*ka*Rab
		a22=1.d0-zi*kb*Rbb
		b1=-2.d0*zi*ka*Raa*exp(zi*ka*a)
		b2=-2.d0*zi*ka*Rab*exp(zi*ka*a)
		d=a11*a22-a12*a21
		x1=(a22*b1-a12*b2)/d
		x2=(a11*b2-a21*b1)/d
		rc=(x1-exp(zi*ka*a))/exp(-zi*ka*a)
		tc=x2/exp(zi*kb*b)

		write(6,*)'rc',rc*Conjg(rc)
		write(6,*)'tc',tc*Conjg(tc)*kb/ka

		write(6,*)'r',(ka-kb)**2/(ka+kb)**2
		write(6,*)'t',(2*ka)**2/(ka+kb)**2*kb/ka



!		exa=exp(-zi*ka*a)
!		exb=exp(zi*kb*b)

!		a11=(1.d0-zi*ka*Raa)*exa
!		a12=-zi*kb*Rab*exb
!		a21=-zi*ka*Rab*exa
!		a22=(1.d0-zi*kb*Rbb)*exb
!		b1=-(1.d0+zi*ka*Raa)/exa
!		b2=-zi*ka*Rab/exa
!		d=a11*a22-a12*a21
!		x1=(a22*b1-a12*b2)/d
!		x2=(a11*b2-a21*b1)/d
!		rc=x1
!		tc=x2
		write(6,*)'rc',rc*Conjg(rc)
		write(6,*)'tc',tc*Conjg(tc)*kb/ka
		write(1,'(2ES16.8)')energy,real(tc*Conjg(tc)*kb/ka)
enddo
close(1)


end


function potential(x)
  implicit none
  real*8,parameter              :: V0=2.5d0,x0=3.d0,mu=0.2d0,d=1.d0
  real*8                        :: w,x,potential
!  w=(x-x0)/d
!  potential=-V0*sinh(w)**2/cosh(w-mu)**2
  potential=1.d0
  if(x.le.0.d0) potential=-1.d0
end function potential



subroutine finite_diff_1(N,a,b,xr,wr,t)
  implicit none
  real*8                             :: x,a,b,wt,h
  integer                            :: i,j,N
  real*8                             :: xr(N),wr(N),t(N,N)

end subroutine finite_diff_1

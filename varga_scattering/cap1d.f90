PROGRAM cap1d
	implicit none
  complex*16,parameter                    :: zi=(0.d0,1.d0)
  real*8,parameter                        :: h2m=0.5d0,pi=3.141592653589793d0,hbar=1.d0
  real*8                                  :: vpl,vpr,kl,kr,al,l0,v0,x0,d,mu,v1,v3,gamm,denominator,width
  integer                                 :: N,i,j,k,l,N_Lattice,i1,i2,en,m,potential_type
  complex*16,dimension(:,:),allocatable   :: am,gm,hm
  complex*16,dimension(:),allocatable     :: phi,psi,wf_exact,w_l,w_r
  real*8,dimension(:),allocatable         :: eva,pot
  real*8,dimension(:,:),allocatable       :: hmat,eve
  real*8                                  :: h,t,tt,a,b,c1,c2,x,cf,ee,dd,yl,yr,dr,c,xl,xr,e_start,de
  complex*16                              :: e,tc,su,cur,d1,transmission
  
  potential_type=3 ! 1=barrier, 2=step, 3=Morse-Feshbach
  N_Lattice=500
  a=-25.d0
  b=+25.d0
  c2=10.d0
  h=(b-a)/(N_Lattice+1)
  e_start=0.d0
  de=0.1d0
  
  allocate(hm(N_Lattice,N_Lattice),am(N_Lattice,N_Lattice),w_l(N_Lattice),w_r(N_Lattice))
  allocate(gm(N_Lattice,N_Lattice),pot(N_Lattice),phi(N_Lattice),psi(N_Lattice),wf_exact(N_Lattice))
  
  ! Parameters for the Morse-Feshbach potential
  m=1.0d0
  V0=2.5d0
  x0=0.0d0
  d=1.0d0
  mu=0.2d0
  V1=-V0*exp(-2.0d0*mu)
	V3=-V0*exp(2.0d0*mu)
  gamm=sqrt((8.0d0*m*d**2.0d0/hbar**2.0d0)*V0*cosh(mu)**2.0d0-1.0d0)  
    
  ! Set up the scattering potential
  select case(potential_type)
    case(1) ! barrier potential     
		  vpl=0.d0
  		vpr=1.d0
			width=2.5d0     
      do i=1,N_Lattice
        x=a+i*h
        if((x<0.d0).OR.(x>width)) then
          pot(i)=vpl       
        else
          pot(i)=vpr
        endif       
      enddo
    case(2) ! step potential
      vpl=-1.d0
      vpr=1.d0
      do i=1,N_Lattice
        x=a+i*h      
        if(x<0.d0) then
          pot(i)=vpl
        else
          pot(i)=vpr
        endif       
      enddo
    case(3) ! Morse-Feshbach potential 
		  do i=1,N_Lattice
  		  x=a+i*h
			  pot(i)=-V0*sinh((x-x0)/d)**2/cosh((x-x0)/d-mu)**2
			enddo
      e_start=-1.55d0
  end select

  ! Output the scattering potential vs. x
  open(1,file="potential.dat")
  do i=1,N_Lattice
    write(1,'(2ES16.8)')a+i*h,pot(i)
  enddo
  write(1,*)
  
  ! Set up the complex absorbing potential and the kinetic energy part
  ! of the Hamiltonian
  c1=-c2
  dr=b-c2
  c=2.62d0
  hm=(0.d0,0.d0)
  t=h2m/h**2
  do i=1,N_Lattice
    x=a+i*h
    w_l(i)=0.d0
    w_r(i)=0.d0
    xl=c*(x-c1)/dr
    xr=c*(x-c2)/dr
    yr=4.d0/(c-xr)**2+4.d0/(c+xr)**2-8.d0/c**2   
    yl=4.d0/(c-xl)**2+4.d0/(c+xl)**2-8.d0/c**2   
    if(x<c1) w_l(i)=-zi*h2m*(2*pi/dr)**2*yl
    if(x>c2) w_r(i)=-zi*h2m*(2*pi/dr)**2*yr
    write(20,*)x,imag(w_l(i)+w_r(i))
    hm(i,i)=2.d0*t
    if(i>1) then
      hm(i,i-1)=-1.d0*t
      hm(i-1,i)=-1.d0*t
    endif
  enddo
  
  ! Add the scattering and complex absorbing potentials to the Hamiltonian
  do i=1,N_Lattice
    hm(i,i)=hm(i,i)+pot(i)+w_l(i)+w_r(i)
  enddo
  
  ! Calculate the transmission for several energies
	do en=1,100
		e=cmplx(e_start+en*de)
		am=-hm
		do i=1,N_Lattice
		  am(i,i)=e-hm(i,i)
		enddo 
		call inv_complex(am,N_Lattice,gm)
		transmission=(0.d0,0.d0)
		do i=1,N_Lattice
	    x=a+(i-1)*h
		  write(30,*)x,real(gm(i,i)),imag(gm(i,i))
		  do j=1,N_Lattice
		    transmission=transmission+abs(gm(i,j))**2*w_r(i)*w_l(j)*4.d0
		  enddo
		enddo

    ! Calculate the analytical transmission coefficient
		select case(potential_type)
		  case(1) ! barrier potential      
				ee=real(e)
				kl=sqrt(ee/h2m-vpl/h2m)
				kr=sqrt(ee/h2m-vpr/h2m)
				
				if(ee<vpr) then
				  kr=sqrt((vpr-ee)/h2m)
				  tc=1.d0/(1.d0+(vpr*sinh(kr*width))**2/(4.d0*ee*(vpr-ee)))
				elseif(ee>vpr) then
				  tc=1.d0/(1.d0+(vpr*sin(kr*width))**2/(4.d0*ee*(ee-vpr)))
				else
				  tc=1.d0/(1.d0+vpr*width**2/(4.d0*h2m))
				endif
				
!				tc=((4.d0*kr*kl*exp(-zi*width*(kl-kr)))/((kl+kr)**2-exp(2.d0*width*zi*kr)*(kl-kr)**2))**2
		  case(2) ! step potential
				ee=real(e)
				kl=sqrt(ee/h2m-vpl/h2m)
				kr=sqrt(ee/h2m-vpr/h2m)
				tc=real(4.d0*kr*kl/(kl+kr)**2)
		  case(3) ! Morse-Feshbach potential
 				ee=real(e)
		    kl=sqrt(2.0d0*(ee-V1))/hbar
    		kr=sqrt(2.0d0*(ee-V3))/hbar
    		denominator=cosh(pi*d*(kl+kr))+cosh(pi*gamm)
		    tc=2.0d0*sinh(pi*d*kl)*sinh(pi*d*kr)/denominator	
		end select
  
		write(6,*)real(e),-real(transmission),real(tc)    
		write(16,*)real(e),-real(transmission),real(tc)

!		do i=1,N_Lattice
!		  x=a+i*h
!		  phi(i)=exp(zi*kl*x)/sqrt(kl)
!		  if(x<0.d0) then
!		    wf_exact(i)=phi(i)+sqrt(rc)*exp(-zi*kl*x)/sqrt(kl)
!		  else
!		    wf_exact(i)=exp(zi*kr*x)*sqrt(tc)
!		  endif
!		enddo

!		do i=1,N_Lattice
!		  x=a+i*h
!		  phi(i)=exp(zi*kl*x)/sqrt(kl)
!		  su=(0.d0,0.d0)
!		  do j=1,N_Lattice
!		    su=su-gm(i,j)*w_l(j)*phi(j)
!		  enddo
!		  psi(i)=su
!		  write(26,*)x,real(su),imag(su)
!		  write(27,*)x,real(wf_exact(i)),imag(wf_exact(i))
!		enddo

!		do i=5,N_Lattice-4
!		  x=a+i*h
!		  d1=1.d0/2.d0*(psi(i+1)-psi(i-1))
!		  cur=conjg(psi(i))*d1/h
!		  write(40,*)x,real(cur),imag(cur)
!		enddo

!		do i=5,N_Lattice-4
!		  x=a+i*h
!			d1=4.d0/5.d0*(psi(i+1)-psi(i-1))-1.d0/5.d0*(psi(i+2)-psi(i-2)) &
!				+4.d0/105.d0*(psi(i+3)-psi(i-3))-1.d0/280.d0*(psi(i+4)-psi(i-4))
!      cur=conjg(psi(i))*d1/h
!	    write(41,*)x,real(cur),imag(cur)
!		enddo
	enddo
END PROGRAM cap1d 

SUBROUTINE inv_complex(A,N,AI)
  implicit none
  integer                :: N,retval,lwork
  integer,allocatable    :: ipiv(:)
  complex*16             :: A(N,N),AI(N,N),temp(1)
  complex*16,allocatable :: work(:)

  AI=A
  allocate(ipiv(N))
  call zgetrf(N,N,AI,N,ipiv,retval)
  call zgetri(N,AI,N,ipiv,temp,-1,retval)
  lwork=int(temp(1))
  allocate(work(lwork))
  call zgetri(N,AI,N,ipiv,work,lwork,retval)
  deallocate(ipiv,work)
END SUBROUTINE inv_complex


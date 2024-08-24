PROGRAM RMATRIX_1D
  implicit none
  real*8,parameter                  :: h2m=0.5d0,pi=3.141592653589793d0
  complex*16,parameter              :: zi=(0.d0,1.d0)
  integer                           :: N,i,j,num_basis_states,num_energies,num_lattice_points_in_basis
  integer                           :: swapped,retval,LWORK,fd_radius,num_lattice_points
  real*8                            :: ka,kb,va,vb,energy,energy_start,energy_step_size,Raa,Rab,Rbb,t2
  real*8                            :: r2,lattice_step_size,temp_eigenvalue,start_coordinate,end_coordinate
  real*8                            :: phi_a,phi_b,t2_analytical,r2_analytical
  real*8,dimension(:),allocatable   :: eigenvalues,xr,wr,work,ev_real,ev_imag,temp_eigenvector
  real*8,dimension(:,:),allocatable :: eigenvectors,t,h,temp
  real*8,external                   :: potential

  ! User-defined values
  fd_radius=1
  num_energies=100
  energy_start=-1.65d0
  energy_step_size=0.05d0
  num_lattice_points=400
  start_coordinate=-5.d0
  end_coordinate=5.d0

  ! Setup the kinetic energy part of the Hamiltonian
  N=num_lattice_points
  allocate(wr(N),xr(N),h(N,N),t(N,N),eigenvectors(N,N),eigenvalues(N))
  allocate(temp(1,N),ev_real(N),ev_imag(N),temp_eigenvector(N))
  if(fd_radius==1) then
    call finite_diff_1(N,start_coordinate,end_coordinate,xr,wr,t)
  elseif(fd_radius==2) then
    call finite_diff_2(N,start_coordinate,end_coordinate,xr,wr,t)
  elseif(fd_radius==3) then
    call finite_diff_3(N,start_coordinate,end_coordinate,xr,wr,t)
  else
    call finite_diff_4(N,start_coordinate,end_coordinate,xr,wr,t)
  endif
  h=-h2m*t

  ! Add the potential to the Hamiltonian
  do i=1,N
    h(i,i)=h(i,i)+potential(xr(i))
  end do

  ! Diagonalize the Hamiltonian
  call dgeev('N','V',N,h,N,ev_real,ev_imag,temp,1,eigenvectors,N,ev_real,-1,retval)
  LWORK=ev_real(1)
  allocate(work(LWORK))
  call dgeev('N','V',N,h,N,ev_real,ev_imag,temp,1,eigenvectors,N,work,LWORK,retval)
  deallocate(work)
  eigenvalues=ev_real+zi*ev_imag

  ! Sort the results by eigenvalue
  swapped=1
  do while(swapped==1)
    swapped=0
    do i=1,(N-1)
      if(real(eigenvalues(i))>real(eigenvalues(i+1))) then
        temp_eigenvalue=eigenvalues(i)
        temp_eigenvector=eigenvectors(:,i)
        eigenvalues(i)=eigenvalues(i+1)
        eigenvectors(:,i)=eigenvectors(:,i+1)
        eigenvalues(i+1)=temp_eigenvalue
        eigenvectors(:,i+1)=temp_eigenvector
        swapped=1
      endif
    enddo
  enddo

  ! For each energy value, calculate the R-matrix elements and the transmission/reflection probabilities
  lattice_step_size=xr(2)-xr(1)
  va=potential(start_coordinate)
  vb=potential(end_coordinate)
  open(11,file="T_vs_E.dat")
  do i=1,num_energies
    energy=energy_start+(i-1)*energy_step_size
    Raa=0.d0
    Rab=0.d0
    Rbb=0.d0
    do j=1,N
      if(fd_radius==1) then
        phi_a=eigenvectors(1,j)
        phi_b=eigenvectors(N,j)
      elseif(fd_radius==2) then
        phi_a=(-4.d0*eigenvectors(1,j)+eigenvectors(2,j))/3.d0
        phi_b=(-4.d0*eigenvectors(N,j)+eigenvectors(N-1,j))/3.d0
      elseif(fd_radius==3) then
        phi_a=(18.d0*eigenvectors(1,j)-9.d0*eigenvectors(2,j)+2.d0*eigenvectors(3,j))/11.d0
        phi_b=(18.d0*eigenvectors(N,j)-9.d0*eigenvectors(N-1,j)+2.d0*eigenvectors(N-2,j))/11.d0
      else
        phi_a=(48.d0*eigenvectors(1,j)-36.d0*eigenvectors(2,j)+16.d0*eigenvectors(3,j)-3.d0*eigenvectors(4,j))/25.d0
        phi_b=(48.d0*eigenvectors(N,j)-36.d0*eigenvectors(N-1,j)+16.d0*eigenvectors(N-2,j)-3.d0*eigenvectors(N-3,j))/25.d0
      endif
      Raa=Raa+phi_a*phi_a/(energy-eigenvalues(j))/lattice_step_size
      Rab=Rab+phi_a*phi_b/(energy-eigenvalues(j))/lattice_step_size
      Rbb=Rbb+phi_b*phi_b/(energy-eigenvalues(j))/lattice_step_size
    enddo
    Raa=-Raa*h2m
    Rab=-Rab*h2m
    Rbb=-Rbb*h2m

    ka=sqrt((energy-va)/h2m)
    kb=sqrt((energy-vb)/h2m)

    t2_analytical=(2*ka)**2/(ka+kb)**2*kb/ka
    r2_analytical=(ka-kb)**2/(ka+kb)**2

    t2=4.d0*ka*kb*(Rab**2)/((1+ka*kb*(Rab**2-Raa*Rbb))**2+(ka*Raa+kb*Rbb)**2)
    r2=((ka*kb*Rab**2)**2-2.d0*Raa*Rbb*(ka*kb*Rab)**2+(ka*kb*Raa*Rbb)**2+1-2.d0*ka*kb*Rab**2+(ka*Raa)**2+(kb*Rbb)**2) &
       /((1+ka*kb*(Rab**2-Raa*Rbb))**2+(ka*Raa+kb*Rbb)**2)

    write(*,*)"Energy=",energy
    write(*,*)"T: ",t2,t2_analytical
    write(*,*)"R: ",r2,r2_analytical
    write(1,*)energy,t2,t2_analytical
    write(11,'(2ES16.8)')energy,t2
  enddo
  close(11)
  deallocate(wr,xr,h,t,eigenvectors,eigenvalues,temp,ev_real,ev_imag,temp_eigenvector)
END PROGRAM RMATRIX_1D

FUNCTION potential(x)
  real*8,parameter              :: V0=2.5d0,x0=3.d0,mu=0.2d0,d=1.d0
  real*8 :: potential,x,w
  if(x>=0.d0) then
    potential=1.d0
  else
    potential=-1.d0
  endif
  
!  w=(x-x0)/d
!  potential=-V0*sinh(w)**2/cosh(w-mu)**2
END FUNCTION potential

SUBROUTINE finite_diff_1(N,a,b,xr,wr,t)
  implicit none
  integer                            :: i,j,N
  real*8                             :: x,a,b,wt,h,xr(N),wr(N),t(N,N)

  h=(b-a)/(1.d0*(N-1))
  do i=1,N
    x=a+(i-1)*h
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
END SUBROUTINE finite_diff_1

SUBROUTINE finite_diff_2(N,a,b,xr,wr,t)
  implicit none
  integer                            :: i,j,N
  real*8                             :: x,a,b,wt,h,xr(N),wr(N),t(N,N)

  h=(b-a)/(1.d0*(N-1))
  do i=1,N
    x=a+(i-1)*h
    xr(i)=x
    wr(i)=h
  end do
    t=0.d0
  do i=1,n
    x=a+(i-1)*h
    t(i,i)=-26.d0/4.d0/h**2
    if(i.gt.1) then
      t(i,i-1)=4.d0/h**2
      t(i-1,i)=4.d0/h**2
    endif
    if(i.gt.2) then
      t(i,i-2)=-3.d0/4.d0/h**2
      t(i-2,i)=-3.d0/4.d0/h**2
    endif
  end do
  t(1,1)=-9.d0/4.d0/h**2
  t(1,2)=12.d0/4.d0/h**2
  t(2,1)=12.d0/4.d0/h**2
  t(2,2)=-25.d0/4.d0/h**2
  t(n,n)=-9.d0/4.d0/h**2
  t(n,n-1)=12.d0/4.d0/h**2
  t(n-1,n)=12.d0/4.d0/h**2
  t(n-1,n-1)=-25.d0/4.d0/h**2
END SUBROUTINE finite_diff_2

SUBROUTINE finite_diff_3(N,a,b,xr,wr,t)
  implicit none
  integer                            :: i,j,N
  real*8                             :: x,a,b,wt,h,xr(N),wr(N),t(N,N)

  h=(b-a)/(1.d0*(N-1))
  do i=1,N
    x=a+(i-1)*h
    xr(i)=x
    wr(i)=h
  end do
    t=0.d0
    do i=1,N
      t(i,i)=-530.d0/36.d0/h**2
      if(i.gt.1) then
        t(i,i-1)=378.d0/36.d0/h**2
        t(i-1,i)=378.d0/36.d0/h**2
      endif
      if(i.gt.2) then
        t(i,i-2)=-135.d0/36.d0/h**2
        t(i-2,i)=-135.d0/36.d0/h**2
      endif
      if(i.gt.3) then
        t(i,i-3)=22.d0/36.d0/h**2
        t(i-3,i)=22.d0/36.d0/h**2
      endif
    end do
!   BC!
    t(1,1)=-121.d0/36.d0/h**2
    t(1,2)=198.d0/36.d0/h**2
    t(2,1)=198.d0/36.d0/h**2
    t(2,2)=-445.d0/36.d0/h**2
    t(1,3)=-99.d0/36.d0/h**2
    t(3,1)=-99.d0/36.d0/h**2
    t(2,3)=360.d0/36.d0/h**2
    t(3,2)=360.d0/36.d0/h**2
    t(3,3)=-526.d0/36.d0/h**2

    t(N,N)=-121.d0/36.d0/h**2
    t(N,N-1)=198.d0/36.d0/h**2
    t(N-1,N)=198.d0/36.d0/h**2
    t(N-1,N-1)=-445.d0/36.d0/h**2
    t(N,N-2)=-99.d0/36.d0/h**2
    t(N-2,N)=-99.d0/36.d0/h**2
    t(N-1,N-2)=360.d0/36.d0/h**2
    t(N-2,N-1)=360.d0/36.d0/h**2
    t(N-2,N-2)=-526.d0/36.d0/h**2
END SUBROUTINE finite_diff_3

SUBROUTINE finite_diff_4(N,a,b,xr,wr,t)
  implicit none
  integer                            :: i,j,N
  real*8                             :: x,a,b,wt,h,xr(N),wr(N),t(N,N)

  h=(b-a)/(1.d0*(N-1))
  do i=1,N
    x=a+(i-1)*h
    xr(i)=x
    wr(i)=h
  end do
    t=0.d0
    do i=1,N
      t(i,i)=-4490.d0/144.d0/h**2
      if(i.gt.1) then
        t(i,i-1)=3552.d0/144.d0/h**2
        t(i-1,i)=3552.d0/144.d0/h**2
      endif
      if(i.gt.2) then
        t(i,i-2)=-1776.d0/144.d0/h**2
        t(i-2,i)=-1776.d0/144.d0/h**2
      endif
      if(i.gt.3) then
        t(i,i-3)=544.d0/144.d0/h**2
        t(i-3,i)=544.d0/144.d0/h**2
      endif
      if(i.gt.4) then
        t(i,i-4)=-75.d0/144.d0/h**2
        t(i-4,i)=-75.d0/144.d0/h**2
      endif
    end do
!   BC!
    t(1,1)=-625.d0/144.d0/h**2
    t(1,2)=1200.d0/144.d0/h**2
    t(2,1)=1200.d0/144.d0/h**2
    t(2,2)=-2929.d0/144.d0/h**2
    t(3,1)=-900.d0/144.d0/h**2
    t(1,3)=-900.d0/144.d0/h**2
    t(3,2)=2928.d0/144.d0/h**2
    t(2,3)=2928.d0/144.d0/h**2
    t(3,3)=-4225.d0/144.d0/h**2
    t(1,4)=400.d0/144.d0/h**2
    t(4,1)=400.d0/144.d0/h**2
    t(2,4)=-1668.d0/144.d0/h**2
    t(4,2)=-1668.d0/144.d0/h**2
    t(3,4)=3504.d0/144.d0/h**2
    t(4,3)=3504.d0/144.d0/h**2
    t(4,4)=-4481.d0/144.d0/h**2

    t(N,N)=-625.d0/144.d0/h**2
    t(N,N-1)=1200.d0/144.d0/h**2
    t(N-1,N)=1200.d0/144.d0/h**2
    t(N-1,N-1)=-2929.d0/144.d0/h**2
    t(N-2,N)=-900.d0/144.d0/h**2
    t(N,N-2)=-900.d0/144.d0/h**2
    t(N-2,N-1)=2928.d0/144.d0/h**2
    t(N-1,N-2)=2928.d0/144.d0/h**2
    t(N-2,N-2)=-4225.d0/144.d0/h**2
    t(N,N-3)=400.d0/144.d0/h**2
    t(N-3,N)=400.d0/144.d0/h**2
    t(N-1,N-3)=-1668.d0/144.d0/h**2
    t(N-3,N-1)=-1668.d0/144.d0/h**2
    t(N-2,N-3)=3504.d0/144.d0/h**2
    t(N-3,N-2)=3504.d0/144.d0/h**2
    t(N-3,N-3)=-4481.d0/144.d0/h**2
END SUBROUTINE finite_diff_4

SUBROUTINE finite_diff_1_alt(N,a,b,xr,wr,t)
  implicit none
  real*8                             :: a,b,wt,h
  integer                            :: i,j,N
  real*8                             :: xr(N),wr(N),t(N,N)

  h=(b-a)/(1.d0*(N-1))
  wr=h
  do i=1,N
    xr(i)=a+(i-1)*h
  enddo

  t=0.d0
  t(1,1)=-1.d0
  t(1,2)=1.d0
  t(1,3)=0.d0
  do i=2,N-1
    t(i,i-1)=1.d0
    t(i,i)=-2.d0
    t(i,i+1)=1.d0
  enddo

  t(N,(N-2))=0.d0
  t(N,(N-1))=1.d0
  t(N,N)=-1.d0

  t=t/h**2
END SUBROUTINE finite_diff_1_alt

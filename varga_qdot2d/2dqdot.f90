!--------------------------------------------------
SUBROUTINE poisson_solve_2d()
!--------------------------------------------------
  USE m_qd2d
  IMPLICIT NONE
  INTEGER :: i1,i2
  rho = -rho
  ! rho=-4.d0*pi*rho
  DO i1=1,N_L(1)
    DO i2=1,N_L(2)
      rho2(i1-1,i2-1) = rho(Lattice_inv(i1,i2))
    ENDDO
  ENDDO 

  CALL hartree_direct(N_L(1),N_L(2),N_L(1)*grid_step(1),N_L(2)*grid_step(2),rho2)
  DO i1 = 1,N_L(1)
    DO i2 = 1,N_L(2)
      v_pot(Lattice_inv(i1,i2)) = E2*rho2(i1-1,i2-1)
    ENDDO
  ENDDO 

END SUBROUTINE 

!------------------------
SUBROUTINE init_lattice()
!------------------------
  USE m_qd2d
  IMPLICIT NONE
  integer :: k1,k2,k3,num,i,k
  num = 0
  DO k1=1,N_L(1)
    DO k2=1,N_L(2)
      num = num + 1 
      Lattice(1,num) = k1
      Lattice(2,num) = k2
      Lattice_inv(k1,k2) = num
      grid_point(1,num) = -0.5d0*(N_L(1)-1)*grid_step(1)+(k1-1)*grid_step(1)
      grid_point(2,num) = -0.5d0*(N_L(2)-1)*grid_step(2)+(k2-1)*grid_step(2)
    ENDDO 
  ENDDO 
END SUBROUTINE 

!------------------------------------
SUBROUTINE init_confining_potential()
!------------------------------------
  USE m_qd2d
  IMPLICIT NONE
  INTEGER :: i
  DO i=1,N_L_points
    V_ext(i) = 0.5d0*omega0**2*( grid_point(1,i)**2 + grid_point(2,i)**2 )
  ENDDO 
END SUBROUTINE

!-----------------------------
SUBROUTINE laplace_operator1()
!-----------------------------
  USE m_qd2d
  IMPLICIT NONE  
  integer          :: i1,i2,kpoint,i
  REAL(8) :: k_x,K_y
  wf=0.d0    
  do i=1,N_L_points
    i1=Lattice(1,i); i2=Lattice(2,i)
    wf(i1,i2)=phi(i)
  end do
  do i=1,N_L_points
    i1=Lattice(1,i); i2=Lattice(2,i)
    K_x=(-2.d0*wf(i1,i2)+wf(i1+1,i2)+wf(i1-1,i2))/grid_step(1)**2
    K_y=(-2.d0*wf(i1,i2)+wf(i1,i2+1)+wf(i1,i2-1))/grid_step(2)**2
    L_phi(i)=K_x+K_y
  end do
end subroutine laplace_operator1

!----------------------------
SUBROUTINE laplace_operator()
!----------------------------
  USE m_qd2d
  IMPLICIT NONE

  REAL(8), parameter :: cN0=-205.d0/72.d0, cN1=8.d0/5.d0, &
                        cN2=-1.d0/5.d0,    cN3=8.d0/315.d0, cN4=-1.d0/560.d0
  integer :: i,i1,i2,i3
  REAL(8) :: K_x,K_y

  do i = 1,N_L_points
    i1 = Lattice(1,i)
    i2 = Lattice(2,i)
    wf(i1,i2)=Phi(i)
  enddo

  do i=1,N_L_points
    i1 = Lattice(1,i)
    i2 = Lattice(2,i)
    K_x = ( cN0*wf(i1,i2) + & 
            cN1*(wf(i1+1,i2) + wf(i1-1,i2)) + & 
            cN2*(wf(i1+2,i2) + wf(i1-2,i2)) + & 
            cN3*(wf(i1+3,i2) + wf(i1-3,i2)) + & 
            cN4*(wf(i1+4,i2) + wf(i1-4,i2)))/Grid_Step(1)**2
    K_y = ( cN0* wf(i1,i2) + & 
            cN1*(wf(i1,i2+1)+wf(i1,i2-1)) + & 
            cN2*(wf(i1,i2+2)+wf(i1,i2-2)) + & 
            cN3*(wf(i1,i2+3)+wf(i1,i2-3)) + & 
            cN4*(wf(i1,i2+4)+wf(i1,i2-4)))/Grid_Step(2)**2
    L_Phi(i) = K_x + K_y
  ENDDO
end subroutine

!--------------------
subroutine CoGr(b,NL)
!--------------------
  USE m_qd2d
  IMPLICIT NONE
	integer :: NL
	REAL(8) :: b(NL)
	REAL(8) :: c0,c1,alfa,beta,rr,bn,con
	integer :: iteration,i,j
	integer, parameter :: N_iter=1500
	REAL(8), parameter :: eps=1.d-10
	REAL(8),  dimension(:), allocatable  :: g,d,h,x

	allocate( g(NL), d(NL), h(NL), x(NL) )

	bn=sqrt(dot_product(b,b))
	x = 0.d0
	phi = x
  call laplace_operator()
  g = b + L_phi
	d = g
	c0 = dot_product(g,g)

	do iteration=1,N_iter
		con=abs(sqrt(dot_product(g,g))/bn)
		write(1,*)iteration,con
		if(con.gt.eps) then
		  phi=d; call laplace_operator; h=L_phi
		  alfa=-c0/dot_product(d,h)
		  x=x+alfa*d
		  g=g+alfa*h
		  c1=dot_product(g,g)
		  beta=c1/c0; c0=c1
		  d=g+beta*d
		endif
	end do
	b=x
	if(con.gt.eps) then
		write(6,*)'Poisson is not converged!'
	endif
	deallocate(g,d,x,h)
end subroutine

!-----------------------------------
SUBROUTINE Hamiltonian_wavefn(ispin)
!-----------------------------------
  USE m_qd2d
  IMPLICIT NONE
  ! calculate H|Phi>
  integer :: ispin
  call laplace_operator()
  if(ispin==1) H_Phi=(V_ext+VH+V_exchange_dw)*Phi-h2m*L_Phi
  if(ispin==2) H_Phi=(V_ext+VH+V_exchange_up)*Phi-h2m*L_Phi

end subroutine

!-----------------------------------
subroutine conjugate_gradient(ispin)
!-----------------------------------
  USE m_qd2d
  IMPLICIT NONE
  
  !  conjugate gradient minimization to diagonalize the Hamiltonian
  
  integer              :: orbital,iteration,i,ispin
  REAL(8)     :: Phi_H_Phi,gamma,beta_Phi,&
         &                beta_beta,beta_H_Phi,&
         &                Phi0_Phi0,Phi0_H_Phi0,&
         &                delta,A,B,C,omega,overlap
  REAL(8), allocatable :: alpha(:),beta(:),Phi0(:),H_Phi0(:)

  allocate(alpha(N_L_points),beta(N_L_points),Phi0(N_L_points),H_Phi0(N_L_points))
    do orbital=1,N_orbitals(ispin)
       Phi=Psi(:,orbital,ispin)
       Phi0=Phi
       call Hamiltonian_wavefn(ispin)
       H_Phi0=H_Phi
       Phi0_H_Phi0=sum(Phi*H_Phi)*Grid_Volume
       Phi0_Phi0=1
       delta=Phi0_H_Phi0
       do iteration=1,N_iteration
          alpha=2*(H_Phi0-delta*Phi0)
          do i=1,orbital-1
             alpha=alpha-Psi(:,i,ispin)*sum(Psi(:,i,ispin)*alpha)*Grid_Volume
          end do
          gamma=sum(alpha*alpha)*Grid_Volume
          if(iteration.eq.1) beta=-alpha
          if(iteration.gt.1) beta=-alpha+gamma/overlap*beta
          overlap=gamma
          beta_Phi=sum(Phi0*beta)*Grid_Volume
          beta_beta=sum(beta*beta)*Grid_Volume
          beta_H_Phi=sum(H_Phi0*beta)*Grid_Volume

          Phi=beta
          call Hamiltonian_wavefn(ispin)
          alpha=H_Phi
          Phi_H_Phi=sum(Phi*H_Phi)*Grid_Volume
          A = Phi_H_Phi*beta_Phi  - beta_H_Phi*beta_beta
          B = Phi_H_Phi*Phi0_Phi0   - Phi0_H_Phi0*beta_beta
          C = beta_H_Phi*Phi0_Phi0 - Phi0_H_Phi0*beta_Phi
          omega=(-B+sqrt(B*B-4*A*C))/(2*A)
          Phi0   = Phi0   +omega*Phi
          H_Phi0 = H_Phi0 +omega*H_Phi
          Phi0_Phi0   =sum(Phi0*Phi0)*Grid_Volume
          Phi0_H_Phi0 =sum(Phi0*H_Phi0)*Grid_Volume
          delta=Phi0_H_Phi0/Phi0_Phi0
       end do
       Phi0_Phi0=sum(Phi0*Phi0)*Grid_Volume
       Psi(:,orbital,ispin)=Phi0/sqrt(Phi0_Phi0)
    end do
    deallocate(alpha,beta,Phi0,H_Phi0)
end subroutine

!----------------------------------
SUBROUTINE Orthogonalization(ispin)
!----------------------------------
  USE m_qd2d
  IMPLICIT NONE
  
  ! Gram-Schmidt orthogonalization
  integer :: orbital,i,ispin
  REAL(8) :: s
  complex(8) :: overlap

  do orbital=1,N_orbitals(ispin)
    do i=1,orbital-1
      overlap=sum(Psi(:,i,ispin)*Psi(:,orbital,ispin))*Grid_Volume
      Psi(:,orbital,ispin)=Psi(:,orbital,ispin)-overlap*Psi(:,i,ispin)
    end do
    s=sum(abs(Psi(:,orbital,ispin))**2)*Grid_Volume
    Psi(:,orbital,ispin)=Psi(:,orbital,ispin)/sqrt(s)
  end do
end subroutine 

!-------------------------------
SUBROUTINE Hamiltonian_density()
!-------------------------------
  USE m_qd2d
  IMPLICIT NONE
  ! Calculate the density dependent part of the Hamiltonian
  rho = density
  CALL poisson_solve_2d()
  vh = v_pot
  call Exchange_Correlation()
END SUBROUTINE

!------------------------
SUBROUTINE Total_Energy()
!------------------------
  USE m_qd2d, ONLY: psi, E_exchange, E_Hartree, Grid_Volume, phi, N_orbitals, h_phi, Etot
  IMPLICIT NONE
  ! Calculate the  Total Energy
  integer :: i,orbital,ispin
  REAL(8) :: E_sp
  
  Etot = 0.d0

  E_sp = 0.d0
  do ispin = 1,2
    do orbital = 1,N_orbitals(ispin)    
      Phi = Psi(:,orbital,ispin)
      call Hamiltonian_wavefn(ispin)
      E_sp = E_sp + sum(Phi*H_Phi)*Grid_Volume
    enddo 
  enddo 

  call Exchange_Correlation()
  
  call Hartree_Energy()
  
  !write(*,*)'E_sp    ', Energy
  !write(*,*)'E_ex    ', E_exchange
  !write(*,*)'E_H     ', E_hartree
  
  Etot = E_sp + E_Exchange + E_Hartree
  
  !write(*,*) 'Etot ', Etot
end subroutine

!----------------------------------
SUBROUTINE Calculate_Density(c1,c2)
!----------------------------------
  USE m_qd2d
  IMPLICIT NONE
  !  Calculate Density  (linear mixing)
  integer :: orbital,i,k,grid,ispin
  REAL(8) :: c1,c2
      
  Density_old = Density
  Density = 0.d0
  Density_up_old = Density_up
  Density_dw_old = Density_dw
  Density_up = 0.d0
  Density_dw = 0.d0 

do ispin=1,2
 do orbital=1,N_orbitals(ispin)
   Phi(:)=Psi(:,orbital,ispin)
   do i=1,N_L_points
    if(ispin==1) Density_dw(i)=Density_dw(i)+Phi(i)*Phi(i)
    if(ispin==2) Density_up(i)=Density_up(i)+Phi(i)*Phi(i)
  end do
  end do
end do     
Density_up=c1*Density_up_old+c2*Density_up
Density_dw=c1*Density_dw_old+c2*Density_dw
Density=Density_up+Density_dw
end subroutine


!--------------------------------
SUBROUTINE Exchange_Correlation()
!--------------------------------
  USE m_qd2d
  USE XC, ONLY: xc_pot
  IMPLICIT NONE
  ! calculate the EXC potential
  integer :: i
  REAL(8), parameter :: aB=0.529177d0,Ryd=13.6058d0 
  REAL(8) :: V_xc,su,eps_c,eps_x,v_x,v_c, &
             rho_up,rho_down,v_x_up,v_x_down,v_c_up,v_c_down

  su=0.d0
  do i=1,N_L_points
    rho_up = density_up(i) + 1.d-20
    rho_down = density_dw(i) + 1.d-20
    call xc_pot(rho_up,rho_down,v_x_up,v_x_down,v_c_up,v_c_down,eps_x,eps_c)
    V_Exchange_up(i) = 2*Ryd*(v_x_up+v_c_up)
    V_Exchange_dw(i) = 2*Ryd*(v_x_down+v_c_down)
    su = su + Density(i)*(eps_x+eps_c) - &
              Density_up(i)*V_exchange_up(i)- &
              Density_dw(i)*V_exchange_dw(i)
  end do
  E_exchange = su*Grid_Volume
end subroutine

!--------------------------
SUBROUTINE Hartree_Energy()
!--------------------------
  USE m_qd2d
  IMPLICIT NONE
  ! calculate the Hartree-Energy
  E_Hartree = -0.5d0*dot_product(density,VH)*grid_volume
end subroutine

!---------------------------------------
subroutine hartree_direct(n1,n2,l1,l2,u)
!---------------------------------------
  IMPLICIT NONE
  integer :: N1,N2
  REAL(8) :: L1,L2
  REAL(8) :: u(0:N1-1,0:N2-1)
  REAL(8) :: rho(0:N1-1,0:N2-1)
  integer :: i,j,ip,jp,ka,kb
  REAL(8)  :: the_sum,g00,d,p,x1,x2,y1,y2,dx1,dx2

  dx1 = l1/dble(n1)
  dx2 = l2/dble(n2)
  g00 = 2.d0*dx1*log( (sqrt(2.d0)+1.d0) / (sqrt(2.d0)-1.d0) )
  rho = u
  do i = 0,N1-1
    do j=0,N2-1
      x1 = i*dx1
      y1 = j*dx2
      the_sum=0.d0
      do ip=0,N1-1
        do jp=0,N2-1
          p=rho(ip,jp)
          if((i==ip).AND.(j==jp)) then
            the_sum=the_sum+p*g00
          else
            x2=ip*dx1
            y2=jp*dx2
            d=sqrt((x1-x2)**2+(y1-y2)**2)
            the_sum=the_sum+p*dx1*dx2/d
          endif
        enddo
      enddo
      u(i,j)=-the_sum
    enddo
  enddo
end subroutine


!------------------------------------------------------------------------------
PROGRAM qdot2d
!------------------------------------------------------------------------------
  USE m_qd2d
  IMPLICIT NONE
  integer :: i,k,i1,i2,i3,orbital,ispin,kk
  REAL(8) :: x,y,z,x0,y0,z0
  REAL(8) :: dx_up(N_L(1))
  REAL(8) :: dy_up(N_L(2))
  REAL(8) :: dx_dw(N_L(1))
  REAL(8) :: dy_dw(N_L(2))
  character(255) :: total_density_filename,energy_vs_iteration_filename
  character(255) :: confining_potential_filename
  
  total_density_filename="total_3d_density.3D.dat"
  confining_potential_filename="confining_potential.3D.dat"
  energy_vs_iteration_filename="energy_vs_iteration.dat"
  
  ! Setup the lattice and initial guess for wavefunctions 
  call init_lattice()

  do ispin=1,2
    do k=1,N_orbitals(ispin)
      call random_number(x0)
      call random_number(y0)
      do i=1,N_L_points
        x = grid_point(1,i) - x0
        y = grid_point(2,i) - y0
        Psi(i,k,ispin) = exp(-0.5d0*omega0*(x**2 + y**2))
      enddo
    enddo
  enddo
  
  call init_confining_potential()
  
  ! Orthogonalize the initial wavefunctions, and use them to calculate the initial density and energy
  call orthogonalization(1) ! spin up
  call orthogonalization(2) ! spin dn
  call calculate_density(0.d0,1.d0)
  call total_energy()

  ! Use the conjugate gradient method to diagonalize the Hamiltonian
  do k = 1,N_scf_iter
    ! Spin up
    call orthogonalization(1)
    call conjugate_gradient(1)
    ! Spin dn
    call orthogonalization(2)
    call conjugate_gradient(2)
    ! Mixing
    call calculate_density(0.5d0, 0.5d0)
    ! New density and potential
    call Hamiltonian_density()
    !
    call total_energy()
    !
    WRITE(*,'(1x,I5,F18.10)') k, Etot
  enddo  

  ! Output the final 3D total density to a Point3D format file, suitable for VisIt
  open(1,file=trim(adjustl(total_density_filename)))
  open(2,file=trim(adjustl(confining_potential_filename)))
  write(1,*)"x y z density"
  write(2,*)"x y z V_ext"
  do i=1,N_L_points
    write(1,'(4ES16.8)')grid_point(1,i),grid_point(2,i),0.d0,density(i)*grid_volume
    write(2,'(4ES16.8)')grid_point(1,i),grid_point(2,i),0.d0,V_ext(i)   
  enddo
  close(1)
  close(2)
  
  write(30,*)N_L
  do i=1,N_L_points
    write(30,*)density(i)
  enddo
  do i=1,N_L_points
    write(30,*)density_up(i)
  enddo
  do i=1,N_L_points
    write(30,*)density_dw(i)
  enddo

  do i1=1,N_L(1)
    do i2=1,N_L(2)
      write(40,*)i1,i2,density(Lattice_inv(i1,i2))
    enddo
    write(40,*)
  enddo
  
  dx_up=0.d0
  dx_dw=0.d0
  dy_up=0.d0
  dy_dw=0.d0

  do i1=1,N_L(1)
    do i2=1,N_L(2)
      dx_up(i1)=dx_up(i1)+density_up(Lattice_inv(i1,i2))
      dx_dw(i1)=dx_dw(i1)+density_dw(Lattice_inv(i1,i2))
      dy_up(i2)=dy_up(i2)+density_up(Lattice_inv(i1,i2))
      dy_dw(i2)=dy_dw(i2)+density_dw(Lattice_inv(i1,i2))
    enddo
  enddo

  do i1=1,N_L(1)
    x=-0.5d0*(N_L(1)-1)*grid_step(1)+(i1-1)*grid_step(1)
    write(10,*)x,dx_up(i1)
    write(11,*)x,dx_dw(i1)
  enddo

  do i2=1,N_L(2)
    y=-0.5d0*(N_L(2)-1)*grid_step(2)+(i2-1)*grid_step(2)
    write(12,*)y,dy_up(i2)
    write(13,*)y,dy_dw(i2)
  enddo

	kk=0
	do ispin=1,2
		do orbital=1,N_orbitals(ispin)
			kk=kk+1
			density=0.d0 
			Phi(:)=Psi(:,orbital,ispin)
			do i=1,N_L_points
				Density(i)=Density(i)+Phi(i)*Phi(i)
			enddo
			do i1=1,N_L(1)
				do i2=1,N_L(2)
				  write(40+kk,*)i1,i2,density(Lattice_inv(i1,i2))
				enddo
				write(40+kk,*)
			enddo
		enddo
	enddo
END PROGRAM


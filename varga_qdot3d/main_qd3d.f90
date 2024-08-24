!------------------------------------------------------------------------------
PROGRAM qd3d
!------------------------------------------------------------------------------
  USE m_qd3d
  implicit none
  integer :: i, k, ispin
  REAL(8) :: x, y, z, x0, y0, z0
  REAL(8), allocatable :: dx_up(:),dy_up(:),dz_up(:),dx_dw(:),dy_dw(:),dz_dw(:)
  character(255) :: total_density_filename
  character(255) :: confining_potential_filename
  real(8) :: Etot, Etot_old, dEtot
  
  CALL init_random_seed()

  total_density_filename       = "total_3d_density.3D.dat"
  confining_potential_filename = "confining_potential.3D.dat"
  
  ! INPUT-----------------------------------
  N_orbitals(1) = 1
  N_orbitals(2) = 1
  omega0 = 2.d0
  N_L = (/ 20, 20, 20 /)
  grid_step = (/ 0.3d0, 0.3d0, 0.3d0 /)
  !-----------------------------------------

  N_L_points = product(N_L)
  dVol = product(grid_step)
  
  ALLOCATE( sp_energy((N_orbitals(1)+N_orbitals(2))) )

  ALLOCATE(Psi(N_L_points,(N_orbitals(1)+N_orbitals(2)),2))
  
  ALLOCATE( Lattice( 3,N_L_points ) )
  ALLOCATE( Lattice_inv( N_L(1), N_L(2), N_L(3) ) )
  ALLOCATE( grid_point( 3, N_L_Points ) )

  ALLOCATE( V_X0( N_L(2), N_L(3) ) )
  ALLOCATE( V_XN( N_L(2), N_L(3) ) )
  ALLOCATE( V_Y0( N_L(1), N_L(3) ) )
  ALLOCATE( V_YN( N_L(1), N_L(3) ) )
  ALLOCATE( V_Z0( N_L(1), N_L(2) ) )
  ALLOCATE( V_ZN( N_L(1), N_L(2) ) )

  ALLOCATE( rho(N_L_points) )

  ALLOCATE( V_POT(N_L_points) )
  
  ALLOCATE( density(N_L_points) )
  ALLOCATE( density_old(N_L_points) )
  ALLOCATE( density_up(N_L_points) )
  
  ALLOCATE( density_up_old(N_L_points) )
  ALLOCATE( density_dw(N_L_points) )
  ALLOCATE( density_dw_old(N_L_points) )
  
  ALLOCATE( phi(N_L_points) )
  ALLOCATE( L_phi(N_L_points) )
  ALLOCATE( VH(N_L_points) )
  ALLOCATE( V_ext(N_L_points) )
  
  ALLOCATE( V_exchange_up(N_L_points) )
  ALLOCATE( V_exchange_dw(N_L_points) )

  ALLOCATE( H_Phi(N_L_Points) )
  
  ALLOCATE( wf(-N_d:N_L(1)+N_d, -N_d:N_L(2)+N_d, -N_d:N_L(3)+N_d) )  
  
  ALLOCATE( dx_up(N_L(1)) )
  ALLOCATE( dy_up(N_L(2)) )
  ALLOCATE( dz_up(N_L(3)) )
  ALLOCATE( dx_dw(N_L(1)) )
  ALLOCATE( dy_dw(N_L(2)) )
  ALLOCATE( dz_dw(N_L(3)) )
    
  ! Setup the lattice and initial guess for wavefunctions  
  call init_lattice()
  do ispin=1,2
    do k=1,N_orbitals(ispin)
      call random_number(x0)
      call random_number(y0)
      call random_number(z0)
      do i=1,N_L_points
        x=grid_point(1,i)-x0
        y=grid_point(2,i)-y0
        z=grid_point(3,i)-z0
        Psi(i,k,ispin)=exp(-0.5d0*omega0*(x**2+y**2+z**2))
      enddo
    enddo
  enddo
  
  call init_confining_potential()
  
  ! Orthogonalize the initial wavefunctions, and use them to calculate the initial density and energy
  call orthogonalization(1)
  call orthogonalization(2)
  call calculate_density(0.d0,1.d0)
  call total_energy(Etot)

  Etot_old = Etot
  dEtot = 0.d0
  ! Use the conjugate gradient method to diagonalize the Hamiltonian
  do k=1,N_scf_iter
    call orthogonalization(1)
    call conjugate_gradient(1)
    call orthogonalization(2)
    call conjugate_gradient(2)
    call calculate_density(0.5d0,0.5d0)
    call total_energy(Etot)
    dEtot = abs(Etot - Etot_old)
    write(*,'(1x,I5,F18.10,ES18.10)') k, Etot, dEtot
    if(dEtot < 1e-7) then
      WRITE(*,*) 'Convergence in Etot is achieved'
      exit
    endif
    Etot_old = Etot
  enddo  
    
  ! Output the final 3D total density to a Point3D format file, suitable for VisIt
  open(1,file=trim(adjustl(total_density_filename)))
  open(2,file=trim(adjustl(confining_potential_filename)))
  write(1,*)"x y z density"
  write(2,*)"x y z V_ext"
  do i=1,N_L_points
    write(1,'(4ES16.8)')grid_point(1,i),grid_point(2,i),grid_point(3,i),density(i)*dVol
    write(2,'(4ES16.8)')grid_point(1,i),grid_point(2,i),grid_point(3,i),V_ext(i)   
  enddo
  close(1)
  close(2)
      
  deallocate(V_X0,V_XN,V_Y0,V_YN,V_Z0,V_ZN,rho,V_POT,density,density_old,density_up)
  deallocate(V_exchange_up,V_exchange_dw,H_Phi,dx_up,dy_up,dz_up,dx_dw,dy_dw,dz_dw)
  deallocate(density_up_old,density_dw,density_dw_old,phi,L_phi,VH,V_ext,wf)
  deallocate(sp_energy,Psi,lattice,lattice_inv,grid_point)
END PROGRAM 


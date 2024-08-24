MODULE m_qd2d
  IMPLICIT NONE

  ! parameters: atomic units
  REAL(8), PARAMETER :: E2=1.0, H2M=0.5d0, a_B=1.d0
  REAL(8), PARAMETER :: Ry=0.5d0, Pi=3.14159265358979d0
  
  ! Number of lattice points 
  INTEGER, PARAMETER :: N_L(2)=(/81,81/)
  INTEGER, PARAMETER :: N_L_points=N_L(1)*N_L(2)
  
  ! Lattice index
  INTEGER :: Lattice(2,N_L_points), Lattice_inv(N_L(1),N_L(2))
  
  ! grid spacing
  REAL(8), PARAMETER :: grid_step(2)=(/0.2d0,0.2d0/)
  REAL(8), PARAMETER :: grid_volume=grid_step(1)*grid_step(2)
  
  ! grid points
  REAL(8) :: grid_point(2,N_L_Points)  
  
  ! boundary conditions 
  REAL(8) :: V_X0(N_L(2)),V_XN(N_L(2))
  REAL(8):: V_Y0(N_L(1)),V_YN(N_L(1))
  
  ! order of finite  difference
  INTEGER, PARAMETER :: N_d=4
  
  ! charge density
  REAL(8) :: rho(N_L_points)
  REAL(8) :: rho2(0:N_L(1)-1,0:N_L(2)-1)
  
  ! potential
  REAL(8) :: V_POT(N_L_points)
  
  ! auxiliary arrays
  REAL(8) :: wf(-N_d:N_L(1) + N_d, -N_d:N_L(2) + N_d)
  REAL(8) :: phi(N_L_points), L_phi(N_L_points), VH(N_L_points), V_ext(N_L_points)

  ! max L in multipole expansion
  INTEGER, PARAMETER :: L_max=4
  REAL(8), PARAMETER :: small=1.d-50

  ! density    
  REAL(8) :: density(N_L_points),density_old(N_L_points)
  REAL(8) :: density_up(N_L_points),density_up_old(N_L_points)
  REAL(8) :: density_dw(N_L_points),density_dw_old(N_L_points)
 
  ! exchange correlation potential  
  REAL(8) :: V_exchange_up(N_L_points),V_exchange_dw(N_L_points)

  ! Number of orbitals
  INTEGER, PARAMETER :: N_orbitals_max=100
  INTEGER, PARAMETER :: N_orbitals_up=1,N_orbitals_dw=1
  INTEGER, PARAMETER :: N_orbitals(2)=(/N_orbitals_up,N_orbitals_dw/)

  ! wave function  
  REAL(8) :: Psi(N_L_Points,N_orbitals_max,2)
  REAL(8) :: H_Phi(N_L_Points)
  INTEGER :: N_iteration=4

  ! Energies 
  REAL(8) :: E_hartree, E_exchange, Etot
  REAL(8), PARAMETER :: omega0=0.5d0
  INTEGER, PARAMETER :: N_scf_iter=80

END MODULE
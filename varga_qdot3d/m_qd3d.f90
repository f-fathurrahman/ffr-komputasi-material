MODULE m_qd3d
  implicit none
  ! parameters: atomic units
  REAL(8), PARAMETER :: E2 = 1.d0
  REAL(8), PARAMETER :: h2m = 0.5d0
  REAL(8), PARAMETER :: a_B = 1.d0
  REAL(8), PARAMETER :: Ry = 0.5d0
  REAL(8), PARAMETER :: pi=3.141592653589793d0
  
  ! Number of lattice points 
  integer :: N_L(3), N_L_points
   
! Lattice index
  integer, allocatable :: Lattice(:,:), Lattice_inv(:,:,:)
! grid points  
  REAL(8), allocatable :: grid_point(:,:)  
! grid spacing
  REAL(8) :: grid_step(3), dVol
 
! boundary conditions 
  REAL(8), allocatable :: V_X0(:,:), V_XN(:,:), V_Y0(:,:), V_YN(:,:), V_Z0(:,:), V_ZN(:,:)
    
  REAL(8), allocatable :: rho(:),V_POT(:),wf(:,:,:),V_exchange_up(:)
  REAL(8), allocatable :: phi(:),L_phi(:),VH(:),V_ext(:),V_exchange_dw(:)
  REAL(8), allocatable :: density(:),density_old(:),density_up(:),H_Phi(:)
  REAL(8), allocatable :: density_up_old(:),density_dw(:),density_dw_old(:)
  
  ! order of finite  difference
  integer, parameter           :: N_d=4
  
  ! max L in multipole expansion
  integer, parameter :: L_max=4
  REAL(8), parameter :: small=1.d-50
    
  REAL(8) :: omega0

  ! Number of orbitals
  integer :: N_orbitals(2)

  ! single_particle energy
  REAL(8), allocatable :: sp_energy(:),Psi(:,:,:)

  integer :: N_iteration=4
  
  ! Energies 
  REAL(8)            :: E_hartree,E_exchange
  integer, parameter :: N_scf_iter = 50

end module
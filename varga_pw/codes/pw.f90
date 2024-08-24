MODULE PW

  IMPLICIT NONE
  ! density
  COMPLEX(8), ALLOCATABLE :: store_density(:,:,:)
  
  !  REAL(8),allocatable    :: density(:,:,:)
  
  COMPLEX(8), ALLOCATABLE :: HPsi(:)

  COMPLEX(4), ALLOCATABLE :: wave_function_c(:,:,:)
  REAL(8), ALLOCATABLE :: Occupation(:,:)
  REAL(8), ALLOCATABLE :: eigenvalue(:,:)
  !
  INTEGER :: N_mesh

  REAL(8)                :: E_total
  REAL(8)                :: E_kinetic
  REAL(8)                :: E_Hartree
  REAL(8)                :: E_exchange
  REAL(8)                :: E_non_local
  REAL(8)                :: E_self
  REAL(8)                :: E_eigen
  REAL(8)                :: E_pseudo
  REAL(8)                :: E_es
  REAL(8)                :: E_Gauss
  REAL(8)                :: E_Fermi
  REAL(8)                :: dens_mix

  REAL(8), PARAMETER      :: rydberg=13.6058d0,pi=3.1415926535897d0
  REAL(8), PARAMETER      :: kt=0.1d0

  INTEGER, PARAMETER               :: N_init=1

  LOGICAL                          :: Band_Structure

  CONTAINS 


!======================
SUBROUTINE initialize()
!======================

  USE GVECTOR
  USE PW_SMALL
  USE PSEUDOPOTENTIAL
  USE FFT_DATA
  
  IMPLICIT NONE 
  
  INTEGER :: i,ik

  CALL input()
  CALL input_atoms()

  CALL initialize_lattice()
  CALL initialize_symmetry()

  CALL generate_G_vectors()
  CALL pw_small_init()
  CALL initialize_pp()

  N_mesh = (N_L(1) + fftinc1) * N_L(2) * N_L(3)
  
  ALLOCATE( store_density( (N_L(1) + fftinc1), N_L(2), N_L(3) ) )
  ! ALLOCATE(density((N_L(1)+fftinc1),N_L(2),N_L(3)))
  
  ALLOCATE( HPsi(N_G_K_vector_max) )
  
  ALLOCATE( wave_function_c(N_G_K_vector_max,N_orbitals,N_K_points) )
  
  ALLOCATE( eigenvalue(N_orbitals,N_K_points) )
  
  ALLOCATE( occupation(N_orbitals,N_K_points) )

  write(16,*) N_orbitals
  
  occupation(:,:) = 0.d0
  DO i=1,N_electrons/2
    DO ik=1,N_K_points
      occupation(i,ik) = 2.0d0
    ENDDO 
  ENDDO 
  IF( (-1)**(N_electrons/2) < 0) THEN 
    DO ik=1,N_K_points
      occupation( N_electrons/2+1, ik ) = 1.0d0
    ENDDO
  ENDIF

  wave_function_c = (0.0, 0.0)
  
  E_kinetic = 0.d0
  
  dens_mix = 0.95d0

END SUBROUTINE 



!=================
SUBROUTINE input()
!=================

  USE GVECTOR
  USE PW_SMALL
  IMPLICIT NONE 
  OPEN(1, file='pw.inp')
  READ(1,*) E_cut,E_cut_small
  READ(1,*) BAND_STRUCTURE
  CLOSE(1)

  WRITE(*,*)
  WRITE(*,*) 'Finished reading pw.inp file in SUBROUTINE input()'
  WRITE(*,*)

END SUBROUTINE


END MODULE



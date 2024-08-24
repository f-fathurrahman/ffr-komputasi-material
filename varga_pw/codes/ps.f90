MODULE PSEUDOPOTENTIAL

  IMPLICIT NONE 
  
  ! number of species
  INTEGER :: N_species
  
  ! number of atoms per species
  INTEGER, ALLOCATABLE :: N_atom(:), z_atom(:)
  INTEGER :: N_atom_max
  
  ! positions
  REAL(8), ALLOCATABLE :: atpos(:,:,:)
  INTEGER :: N_orbitals
  INTEGER :: N_electrons,N_conduction
  
  ! max L in PS
  INTEGER, PARAMETER :: L_pp_max=4
  
  ! max number of radial points
  INTEGER, PARAMETER :: N_pp_max=1000
  
  ! pseudocore charge
  REAL(8), ALLOCATABLE :: charge_pp(:)
  
  ! gaussian
  REAL(8), PARAMETER :: beta=1.d0
  
  ! PP
  REAL(8), ALLOCATABLE :: vion(:,:,:),p(:,:,:),r(:,:), cclog(:)
  INTEGER, ALLOCATABLE :: l_loc(:),l_max(:),N_pp(:)
  REAL(8), ALLOCATABLE :: rhops(:,:),vps(:,:)
  COMPLEX(8), ALLOCATABLE :: fnl(:,:,:,:,:),ei1(:,:,:),ei2(:,:,:),ei3(:,:,:),eigr(:,:,:,:),sfac(:,:)
  REAL(8), ALLOCATABLE :: wnl(:,:)
  REAL(8),ALLOCATABLE :: pkg_a(:,:,:,:)



CONTAINS

!=======================
SUBROUTINE input_atoms()
!=======================
  IMPLICIT NONE 
  INTEGER :: i,j,k

  WRITE(6,*) '?1'
  
  N_atom_max = 0
  
  OPEN(1,file='atom.inp')
  READ(1,*) N_electrons,N_conduction
  READ(1,*) n_species
  
  ALLOCATE( N_atom(N_species), charge_pp(N_species), z_atom(N_species) )
  
  DO i=1,n_species
    READ(1,*) n_atom(i), z_atom(i)
    IF( n_atom(i) > N_atom_max ) N_atom_max = n_atom(i)
  ENDDO 
  
  
  ALLOCATE( atpos(3,N_atom_max,N_species) )
  DO i = 1,N_species
    DO j = 1,n_atom(i)
      READ(1,*) ( atpos(k,j,i), k=1,3 )
    ENDDO 
  ENDDO 

  N_orbitals = N_electrons/2 + N_conduction
  
  IF( (-1)**(N_electrons/2) < 0) N_orbitals = N_electrons/2 + 1 + N_conduction

END SUBROUTINE 


!===================
SUBROUTINE read_pp()
!===================
  IMPLICIT NONE 
  INTEGER :: total_charge,i,l,j,is

  ALLOCATE( l_loc(N_species), l_max(N_species) )
  ALLOCATE( vion(N_pp_max,N_species,L_pp_max) )
  ALLOCATE( N_pp(N_species), cclog(N_species) )
  ALLOCATE( p(N_pp_max,N_species,3), r(N_pp_max,N_species) )
  
  total_charge = 0
  i = 0
  DO is = 1,n_species
    IF(z_atom(is) == 1)  OPEN(1,file='hpp.dat')
    IF(z_atom(is) == 3)  OPEN(1,file='lipp.dat')
    IF(z_atom(is) == 4)  OPEN(1,file='bepp.dat')
    IF(z_atom(is) == 5)  OPEN(1,file='bpp.dat')
    IF(z_atom(is) == 6)  OPEN(1,file='cpp.dat')
    IF(z_atom(is) == 7)  OPEN(1,file='npp.dat')
    IF(z_atom(is) == 8)  OPEN(1,file='opp.dat')
    IF(z_atom(is) == 11) OPEN(1,file='napp.dat')
    IF(z_atom(is) == 13) OPEN(1,file='alpp.dat')
    IF(z_atom(is) == 14) OPEN(1,file='sipp.dat')
    IF(z_atom(is) == 16) OPEN(1,file='spp.dat')
    IF(z_atom(is) == 28) OPEN(1,file='nipp.dat')
    IF(z_atom(is) == 31) OPEN(1,file='gapp.dat')
    IF(z_atom(is) == 33) OPEN(1,file='aspp.dat')

    READ(1,*) charge_pp(is), l_max(is), l_loc(is)

    total_charge = total_charge + n_atom(is) * charge_pp(is)

    DO l=1,l_max(is)
      READ(1,*) N_pp(is),cclog(is),(i,r(j,is),p(j,is,l),vion(j,is,l),j=1,N_pp(is))
    ENDDO
    cclog(is) = log( cclog(is) )
    CLOSE(1)
  ENDDO

END SUBROUTINE 

!=========================
SUBROUTINE initialize_pp()
!=========================
  
  USE GVECTOR
  USE PW_SMALL

  IMPLICIT NONE 

  ALLOCATE( rhops(N_species,N_G_vector_max) )
  ALLOCATE( vps(N_species,N_G_vector_max) )
  ALLOCATE( fnl(N_G_wf_max,N_K_points,N_species,N_atom_max,(L_PP_max+1)**2) )
  ALLOCATE( wnl(N_species,(L_PP_max+1)**2) )
  ALLOCATE( pkg_a((L_PP_max+1)**2,N_G_K_vector_max,N_species,N_K_points) )
  ALLOCATE( ei1(-(N_L(1)+1)/2:(N_L(1)+1)/2,N_atom_max,N_species) )
  ALLOCATE( ei2(-(N_L(2)+1)/2:(N_L(2)+1)/2,N_atom_max,N_species) )
  ALLOCATE( ei3(-(N_L(3)+1)/2:(N_L(3)+1)/2,N_atom_max,N_species) )
  ALLOCATE( eigr(N_G_K_vector_max,N_atom_max,N_species,N_k_points) )
  ALLOCATE( sfac(N_species,N_G_vector_max) )

  CALL read_pp()

END SUBROUTINE 


END MODULE 


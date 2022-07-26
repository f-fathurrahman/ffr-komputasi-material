!-------------
MODULE globals
!-------------

IMPLICIT NONE

REAL(8):: PI 
COMPLEX(8):: Im

!#################################################
! Geometry and basis set size. Values are set during initialisation
!#################################################
REAL(8):: BoxL, Omega,& ! Linear Box Size
          TimeStepOrt, TimeStepCP, RedFac, mu, Ediff
INTEGER :: MaxBas, GridSize, &! Maximum number of waves along 
           DiagSize, NoOfPW,   &! LINEAR direction, Number of Plane Waves
           MaxIter
LOGICAL :: PrintOut
LOGICAL :: OnlyStatic

!#################################################
! Physical fields
!#################################################
COMPLEX(8), ALLOCATABLE :: Density_K(:,:,:)
COMPLEX(8), ALLOCATABLE :: Density_R(:,:,:)

!#################################################
! Physical system: nr of ions, electrons and orbitals
!#################################################
INTEGER :: N_ion, N_electron, N_orbitals


!#################################################
! Wavefunction coefficients
!#################################################
! complex(8), ALLOCATABLE :: Coeffs_K(:,:,:,:), Coeffs_R(:,:,:,:)

!#################################################
! Stored values of pseudopot.
! core charges and short range part of local pseudopot
!#################################################
COMPLEX(8), ALLOCATABLE :: CoreCharge(:,:,:,:),  &
                           totCoreCharge(:,:,:), &
                           ShortLocal(:,:,:,:), &
                           NonLocal(:,:,:), totShortLocal(:,:,:), &
                           PseudoGrid(:,:,:)

!#################################################
! Data of ions, stored Coulomb potential and kinetic
! term.
!#################################################
TYPE Type_Ions
  INTEGER :: AtomNum
  real(8):: Mass
  real(8):: R_I(3)
END TYPE Type_Ions
TYPE(Type_Ions), ALLOCATABLE :: Ions(:)

REAL(8), ALLOCATABLE :: Gmin2Grid(:,:,:), G2Grid (:), FillFac(:)


!#################################################
! Connection between linear indices of H_KS and 
! grid positions 
!#################################################
INTEGER, ALLOCATABLE :: GridIndex(:,:,:), GridPos(:,:), GGrid(:,:,:,:)

!#################################################
! Cut-off in reciprocal space
!#################################################
REAL(8):: GMax


!#################################################
! Pseudopotential parameters
!#################################################
TYPE Type_PP_Params
  INTEGER :: AtomNum
  INTEGER :: Zion
  INTEGER :: N_nonzero_C_i
  REAL(8):: Xi
  REAL(8):: C(4)
  REAL(8):: MaxL
  REAL(8):: r_s
  REAL(8):: h_1s, h_2s, r_p, h_1p
END TYPE Type_PP_Params

INTEGER :: No_OF_DIFF_IONS !Has to be assigned in main/InCP

TYPE(Type_PP_Params), ALLOCATABLE :: PP_params(:) 


END MODULE Globals



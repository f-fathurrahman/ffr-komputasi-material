MODULE globals

 IMPLICIT NONE

 DOUBLE PRECISION :: PI 
 DOUBLE COMPLEX   :: Im

 !#################################################
 ! Geometry and basis set size. Values are set during initialisation
 !#################################################

 DOUBLE PRECISION :: BoxL, Omega,& ! Linear Box Size
                    TimeStepOrt, TimeStepCP, RedFac, mu, Ediff
 INTEGER :: MaxBas, GridSize, &! Maximum number of waves along 
           DiagSize, NoOfPW,   &! LINEAR direction, Number of Plane Waves
           MaxIter
 LOGICAL :: PrintOut
 LOGICAL :: OnlyStatic
 !#################################################
 ! Physical fields
 !#################################################

 DOUBLE COMPLEX, ALLOCATABLE :: Density_K(:,:,:)
 DOUBLE COMPLEX, ALLOCATABLE :: Density_R(:,:,:)

 !#################################################
 ! Physical system: nr of ions, electrons and orbitals
 !#################################################

 INTEGER :: N_ion, N_electron, N_orbitals

 !#################################################
 ! Wavefunction coefficients
 !#################################################

 ! DOUBLE COMPLEX, ALLOCATABLE :: Coeffs_K(:,:,:,:), Coeffs_R(:,:,:,:)

 !#################################################
 ! Stored values of pseudopot.
 ! core charges and short range part of local pseudopot
 !#################################################
 DOUBLE COMPLEX, ALLOCATABLE :: CoreCharge(:,:,:,:),  &
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
   DOUBLE PRECISION :: Mass
   DOUBLE PRECISION :: R_I(3)
 END TYPE Type_Ions
 TYPE(Type_Ions), ALLOCATABLE :: Ions(:)
 
 DOUBLE PRECISION, ALLOCATABLE :: Gmin2Grid(:,:,:), &
                                 G2Grid (:), FillFac(:)
 !#################################################
 ! Connection between linear indices of H_KS and 
 ! grid positions 
 !#################################################


 INTEGER, ALLOCATABLE :: GridIndex(:,:,:), GridPos(:,:), GGrid(:,:,:,:)

 !#################################################
 ! Cut-off in reciprocal space
 !#################################################

 DOUBLE PRECISION :: GMax


 !#################################################
 ! Pseudopotential parameters
 !#################################################

 TYPE Type_PP_Params
   INTEGER :: AtomNum
   INTEGER :: Zion
   INTEGER :: N_nonzero_C_i
   DOUBLE PRECISION :: Xi
   DOUBLE PRECISION :: C(4)
   DOUBLE PRECISION :: MaxL
   DOUBLE PRECISION :: r_s
   DOUBLE PRECISION :: h_1s, h_2s, r_p, h_1p
 END TYPE Type_PP_Params

 INTEGER :: No_OF_DIFF_IONS !Has to be assigned in main/InCP

 TYPE(Type_PP_Params), ALLOCATABLE :: PP_params(:) 


END MODULE Globals

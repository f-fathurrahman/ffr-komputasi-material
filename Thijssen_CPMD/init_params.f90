!
! Read some global parameters from file InCP and initialize
! some other variables
!
!-----------------------
SUBROUTINE init_params()
!-----------------------
  use globals
  IMPLICIT NONE
  REAL(8) :: Energy_CutOff
  INTEGER :: I, I2, I3, I5, N

  ! Some constants
  PI = 4.D0*ATAN(1.D0)
  Im = CMPLX(0.D0, 1.D0, kind=8)

  OPEN(8,File="InCP")
  READ(8,*) OnlyStatic
  READ(8,*) BoxL
  Omega = BoxL**3
  READ(8,*) Energy_CutOff
  GMax = SQRT(2.D0*Energy_Cutoff)
  GMax = GMax*BoxL*0.5D0/PI
  ! Determine the FFT grid size here
  I2 = ceiling(log(4*Gmax)/Log(2.D0))
  I2 = 2**I2
  I3 = ceiling(log(4*Gmax)/Log(3.D0))
  I3 = 3**I3
  I5 = ceiling(log(4*Gmax)/Log(5.D0))
  I5 = 5**I5
  GridSize = MIN(I2,I3,I5)
  !
  WRITE(*,*) 'GridSize = ', GridSize
  !
  MaxBas = INT(GMax)
  ALLOCATE(Density_K(0:GridSize-1,0:GridSize-1,0:GridSize-1))
  ALLOCATE(Density_R(0:GridSize-1,0:GridSize-1,0:GridSize-1))
  !
  CALL FFT_Plan(GridSize)
  !
  WRITE(*,'(1x,A,I8)') 'Linear index of PWs runs from 0 to ', GridSize-1 
  
  READ(8,*) mu
  READ(8,*) TimeStepOrt      
  READ(8,*) TimeStepCP
  READ(8,*) Ediff
  READ(8,*) MaxIter
  READ(8,*) RedFac
  READ(8,*) No_OF_DIFF_IONS
  
  WRITE(*,*) 'There are ',  No_OF_DIFF_IONS, ' different ions'
  
  ! Allocate memory for pseudopotential data
  ALLOCATE( PP_Params(No_OF_DIFF_IONS) )
  PP_Params(:)%AtomNum = 0

  READ(8, *) N_ion
  ALLOCATE(Ions(N_ion))
  DO I = 1, N_ion
    READ(8,*) Ions(I)%AtomNum, &  ! Atomic Number
               Ions(I)%Mass, &     ! Mass of Ion
               Ions(I)%R_I(1), Ions(I)%R_I(2), Ions(I)%R_I(3) ! Positions
    Ions(I)%Mass = Ions(I)%Mass*1836.D0 ! convert to amu
    CALL init_pp_params(Ions(I)%AtomNum)
  ENDDO
  WRITE(*,*) 'No of ions', N_ion
  DO i = 1, N_ion
    WRITE(*,'(1x,I4,4F18.10)') ions(i)%AtomNum, ions(i)%Mass, &
                            &  ions(i)%R_I(1), ions(i)%R_I(2), ions(i)%R_I(3)
  ENDDO 
  
  CALL init_grids()

  ! Specify number of electrons and also the occupation numbers
  READ(8,*) N_electron
  READ(8,*) N_Orbitals
  ALLOCATE(FillFac(N_orbitals))
  DO N = 1, N_orbitals
    READ(8,*) FillFac(N)
  ENDDO
  CLOSE(8)

  WRITE(*,*) 'Finished reading input file'

END SUBROUTINE

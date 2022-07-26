!----------------------
SUBROUTINE init_grids()
!----------------------
  use globals
  ! The various grids are allocated, initialized and subsequently the grids that
  ! are independent of the position of the Ions are filled, FillStaticGrids,
  ! then the grids that are dependent on the position are filled, FillDynGrids
  IMPLICIT NONE
  INTEGER :: I, J, K
  REAL(8) :: G2
  ! Functions
  REAL(8) :: G2_Short

  NoOfPW = 0
  DO I=0, GridSize-1
  DO J=0, GridSize-1
  DO K=0, GridSize-1
    G2 = G2_Short(I,J,K)
    IF( G2 < Gmax**2 ) NoOfPW = NoOfPW + 1
  ENDDO
  ENDDO
  ENDDO

  ALLOCATE(PseudoGrid(0:GridSize-1, 0:GridSize-1, 0:GridSize-1))
  ALLOCATE(CoreCharge(N_ion,0:GridSize-1, 0:GridSize-1, 0:GridSize-1))
  ALLOCATE(totCoreCharge(0:GridSize-1, 0:GridSize-1, 0:GridSize-1))
  ALLOCATE(ShortLocal(N_ion,0:GridSize-1, 0:GridSize-1, 0:GridSize-1))
  ALLOCATE(totShortLocal(0:GridSize-1, 0:GridSize-1, 0:GridSize-1))
  ALLOCATE(NonLocal(NoOfPW,N_ion,5))
  ALLOCATE(GridIndex(0:GridSize-1, 0:GridSize-1, 0:GridSize-1))
  ALLOCATE(GridPos(NoOfPW,3))
  ALLOCATE(G2Grid(NoOfPW))
  ALLOCATE(GGrid(3,0:GridSize-1, 0:GridSize-1, 0:GridSize-1))
  ALLOCATE(Gmin2Grid(0:GridSize-1, 0:GridSize-1, 0:GridSize-1))

  ! Initialize G-vectors and related quantities
  CALL fill_static_grids()

  ! Initialize structure factors,
  ! local potentials, nonlocal potentials due to ions
  CALL fill_dyn_grids()   
    
END SUBROUTINE 


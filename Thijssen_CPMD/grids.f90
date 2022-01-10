MODULE grids

USE globals
USE pseudo
USE utilities

CONTAINS

  SUBROUTINE InitGrids()
  ! The various grids are allocated, initialized and subsequently the grids that
  ! are independent of the position of the Ions are filled, FillStaticGrids,
  ! then the grids that are dependent on the position are filled, FillDynGrids
  IMPLICIT NONE
    INTEGER :: I, J, K
    DOUBLE PRECISION :: G2
    
    NoOfPW = 0
    DO I=0, GridSize-1
      DO J=0, GridSize-1
        DO K=0, GridSize-1
          G2 = G2_Short(I,J,K)
          IF (G2<Gmax**2) NoOfPW = NoOfPW + 1
        END DO
      END DO
    END DO
    
    ALLOCATE (PseudoGrid(0:GridSize-1, 0:GridSize-1, 0:GridSize-1))
    ALLOCATE (CoreCharge(N_ion,0:GridSize-1, 0:GridSize-1, 0:GridSize-1))
    ALLOCATE (totCoreCharge(0:GridSize-1, 0:GridSize-1, 0:GridSize-1))
    ALLOCATE (ShortLocal(N_ion,0:GridSize-1, 0:GridSize-1, 0:GridSize-1))
    ALLOCATE (totShortLocal(0:GridSize-1, 0:GridSize-1, 0:GridSize-1))
    ALLOCATE (NonLocal(NoOfPW,N_ion,5))
    ALLOCATE (GridIndex(0:GridSize-1, 0:GridSize-1, 0:GridSize-1))
    ALLOCATE (GridPos(NoOfPW,3))
    ALLOCATE (G2Grid(NoOfPW))
    ALLOCATE (GGrid(3,0:GridSize-1, 0:GridSize-1, 0:GridSize-1))
    ALLOCATE (Gmin2Grid(0:GridSize-1, 0:GridSize-1, 0:GridSize-1))
    
    CALL FillStaticGrids()
    CALL FillDynGrids()   
    
  END SUBROUTINE InitGrids


  SUBROUTINE FillDynGrids()
  ! In this subroutine the grids are filled that are dependent on the position
  ! of the Ions
    IMPLICIT NONE
  
    INTEGER :: I, J, K, N, G2, Ind, I1, J1, K1, PP_Index, AtomNum
    DOUBLE COMPLEX :: StructFac
    totShortLocal = CMPLX(0.D0)
    PseudoGrid = CMPLX(0.D0)
    ShortLocal = CMPLX(0.D0)
    totCoreCharge = CMPLX(0.D0)
    NonLocal = CMPLX(0.D0)

    Ind = 0
    DO I=0, GridSize-1
      DO J=0, GridSize-1
        DO K=0, GridSize-1
          G2 = G2_Short(I,J,K)
          IF (G2<Gmax**2) ind = ind + 1
            DO N = 1, N_ion
              AtomNum = Ions(N)%AtomNum
              PP_Index = GetIndexPP(AtomNum)
              I1 = I-INT((2.D0*I)/GridSize)*GridSize
              J1 = J-INT((2.D0*J)/GridSize)*GridSize
              K1 = K-INT((2.D0*K)/GridSize)*GridSize
              CALL CalcFacs(I, J, K, N, G2, StructFac)
              PseudoGrid(I,J,K)=PseudoGrid(I,J,K)+ & 
                  Local_PP(G2,Ions(N)%AtomNum)*StructFac
              ShortLocal(N,I,J,K) = ShortLocal(N,I,J,K) +&
                                Short_PP(G2,Ions(N)%AtomNum)*StructFac
              totShortLocal(I,J,K) = totShortLocal(I,J,K) + &
                                ShortLocal(N,I,J,K)
              CoreCharge(N,I,J,K) = CoreDens(G2,Ions(N)%AtomNum)*StructFac
              totCoreCharge(I,J,K) = totCoreCharge(I,J,K) + CoreCharge(N,I,J,K)
              IF (G2<Gmax**2) THEN
                NonLocal(Ind,N,1) = NonLoc(I1,J1,K1,AtomNum,0,0,1)*StructFac
                IF (PP_Params(PP_Index)%MaxL>0) THEN
                  NonLocal(Ind,N,2) = NonLoc(I1,J1,K1,AtomNum,0,0,2)*StructFac
                  NonLocal(Ind,N,3) = NonLoc(I1,J1,K1,AtomNum,1,1,1)*StructFac
                  NonLocal(Ind,N,4) = NonLoc(I1,J1,K1,AtomNum,1,0,1)*StructFac
                  NonLocal(Ind,N,5) = NonLoc(I1,J1,K1,AtomNum,1,-1,1)*StructFac
                END IF
              END IF !G2<Gmax**2
          END DO !N
        END DO !K
      END DO !J
    END DO !I
    
  END SUBROUTINE FillDynGrids
  
  SUBROUTINE FillStaticGrids()
  ! In this subroutine the grids are defined that are independent on the
  ! position of the Ions
  ! G2 stands for G squared, Gmin2 for 1/G Squared

    INTEGER          :: I, J, K, Ind
    INTEGER          :: II, JJ, KK
    DOUBLE PRECISION :: G2
    Ind = 0
    
    GGrid = 0
    GridIndex = 0
    GridPos = 0
    DO I=0, GridSize-1
      DO J=0, GridSize-1
        DO K=0, GridSize-1
          G2 = G2_Short(I,J,K)
          II = I-INT((2.D0*I)/GridSize)*GridSize
          JJ = J-INT((2.D0*J)/GridSize)*GridSize
          KK = K-INT((2.D0*K)/GridSize)*GridSize
          IF (G2<4*GMax**2) THEN
            GGrid(:,I,J,K)=(/II,JJ,KK/)
          ELSE
            GGrid(:,I,J,K) = (/0, 0, 0/)
          END IF
          IF (G2<Gmax**2) THEN
            Ind = Ind + 1
            GridIndex(I, J, K) = Ind
            GridPos(Ind,:) = (/I,J,K/)
            G2Grid(Ind) = G2
          END IF
          IF (G2/=0) THEN
            Gmin2Grid(I,J,K) = 1.D0/G2
          ELSE
            Gmin2Grid(I,J,K) = 0.D0
          END IF
        END DO !K
      END DO !J
    END DO !I
  
  END SUBROUTINE FillStaticGrids


END MODULE grids

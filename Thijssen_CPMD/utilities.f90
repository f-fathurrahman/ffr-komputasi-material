MODULE utilities

USE globals

CONTAINS
  
  SUBROUTINE Gram_Schmidt(Vectors, Number, Dimen)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: Number, Dimen
  DOUBLE COMPLEX, INTENT(INOUT) :: Vectors(Number, Dimen)
  INTEGER Iorb1, Iorb2
  DOUBLE COMPLEX :: IP

  DO Iorb1 = 1, Number
    DO Iorb2 = 1, Iorb1-1
      IP = DOT_PRODUCT(Vectors(Iorb2,:),Vectors(Iorb1,:))
      Vectors(Iorb1,:) = Vectors(IOrb1,:)-IP*Vectors(Iorb2,:)
    END DO
    IP = DOT_PRODUCT(Vectors(Iorb1,:),Vectors(Iorb1,:))
    IP = 1/SQRT(IP)
    Vectors(Iorb1,:) = Vectors(Iorb1,:)*IP
  END DO

  END SUBROUTINE Gram_Schmidt

  

  DOUBLE COMPLEX FUNCTION InnerProd(Arr1, Arr2)

  IMPLICIT NONE

  DOUBLE COMPLEX :: Arr1(0:GridSize-1,0:GridSize-1,0:GridSize-1), &
                    Arr2(0:GridSize-1,0:GridSize-1,0:GridSize-1)
  InnerProd = SUM(CONJG(Arr1)*Arr2)
  END FUNCTION InnerProd  




  SUBROUTINE CalcDensAndCoeffs_R(NZCoeffs, Coeffs_R)
! Calculates the density_K from Coeffs_K, and normalize all
  IMPLICIT NONE
  DOUBLE COMPLEX, INTENT(IN) :: NZCoeffs(N_orbitals, NoOfPW)
  DOUBLE COMPLEX, INTENT(INOUT) :: Coeffs_R(N_orbitals,0:GridSize-1,0:GridSize-1,0:GridSize-1)
  DOUBLE COMPLEX, ALLOCATABLE :: Coeffs_K(:,:,:,:)
  
  INTEGER :: ElecCnt, IIndex, N
  DOUBLE PRECISION :: Norm
  
  ALLOCATE(Coeffs_K(N_orbitals,0:GridSize-1,0:GridSize-1,0:GridSize-1)) 
  Coeffs_K = CMPLX(0.D0)
  DO IIndex = 1, NoOfPW
    Coeffs_K(:, GridPos(IIndex,1), GridPos(IIndex,2), GridPos(IIndex,3)) = &
         NZCoeffs(:,IIndex)
  END DO 
  DO N=1, N_orbitals
    CALL Forward_FFT(GridSize, Coeffs_K(N,:,:,:),Coeffs_R(N,:,:,:))
  END DO
  Coeffs_R = Coeffs_R/SQRT(Omega)
  Density_R = CMPLX(0.D0)
  DO N = 1, N_orbitals
    Density_R = Density_R + FillFac(N)*Coeffs_R(N,:,:,:)*CONJG(Coeffs_R(N,:,:,:))
  END DO
  CALL Backward_FFT(GridSize, Density_R, Density_K)
  END SUBROUTINE  CalcDensAndCoeffs_R

  



  SUBROUTINE Rattle (NZCoeffs, OldNZCoeffs, NZCoeffsdot)
  IMPLICIT NONE
  DOUBLE COMPLEX, INTENT(INOUT) :: NZCoeffs(N_Orbitals, NoOfPW), &
                                   NZCoeffsdot(N_Orbitals, NoOfPW), &
                                   OldNZCoeffs(N_Orbitals, NoOfPW)

  INTEGER :: N, J, K
  DOUBLE PRECISION :: Norm
  DOUBLE COMPLEX, ALLOCATABLE :: I(:,:), A(:,:), B(:,:), &
                                 X(:,:), Correction(:,:)
  IF (.NOT.ALLOCATED(A)) THEN
    ALLOCATE (A(N_Orbitals, N_Orbitals))
    ALLOCATE (B(N_Orbitals, N_Orbitals))
    ALLOCATE (X(N_Orbitals, N_Orbitals))
    ALLOCATE (Correction(N_Orbitals, N_Orbitals))
    ALLOCATE (I(N_Orbitals, N_Orbitals))
  END IF

  I = CMPLX(0.D0)
  DO N = 1, N_Orbitals 
    I(N,N) = CMPLX(1.D0)
  END DO

  A = MATMUL(CONJG(NZCoeffs),TRANSPOSE(NZCoeffs))
  B = MATMUL(CONJG(OldNZCoeffs),TRANSPOSE(NZCoeffs))

  X = 0.5D0*(I-A)

  DO 
    Correction = I-A - MATMUL(TRANSPOSE(CONJG(B)), X) - MATMUL(X, B) - &
               MATMUL(X,X)
    X = X + 0.5*Correction
    Norm = SUM(CONJG(Correction)*Correction)
    IF (Norm<1.D-10) EXIT  
  END DO 
  NZCoeffs = NZCoeffs + MATMUL(CONJG(X),OldNZCoeffs)
  

  END SUBROUTINE Rattle





  DOUBLE PRECISION FUNCTION G2_Short (I,J,K)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: I, J, K
  INTEGER :: II, JJ, KK

  II = I-INT((2.D0*I)/GridSize)*GridSize
  JJ = J-INT((2.D0*J)/GridSize)*GridSize
  KK = K-INT((2.D0*K)/GridSize)*GridSize

  G2_Short = II*II+JJ*JJ+KK*KK
  END FUNCTION G2_Short




  SUBROUTINE CalcFacs(I,J,K,N,G2,StructFac)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: I, J, K, N
  INTEGER, INTENT(OUT) :: G2
  DOUBLE COMPLEX, INTENT(OUT) :: StructFac
  INTEGER :: II, JJ, KK
  INTEGER :: pos(3)

  II = I-INT((2.D0*I)/GridSize)*GridSize
  JJ = J-INT((2.D0*J)/GridSize)*GridSize
  KK = K-INT((2.D0*K)/GridSize)*GridSize

  G2 = II*II+JJ*JJ+KK*KK
  pos = (/II, JJ, KK/)
  StructFac = EXP(-Im*2*PI*DOT_PRODUCT(Ions(N)%R_I,pos)/BoxL)
  END SUBROUTINE CalcFacs

  SUBROUTINE PrintIons()
  IMPLICIT NONE
  INTEGER I
  OPEN (UNIT=11, FILE='Ions.dat', POSITION='APPEND')
  DO I= 1, N_Ion
        print '(A23 I3)', 'Ion Nr:', I
        print '(A23 F15.8 F15.8 F15.8)', 'Ion R_I:', Ions(I)%R_I
        print *
        WRITE (11,'(3F15.8 $)') Ions(I)%R_I
  END DO
  WRITE (11,*)
  CLOSE (11) 
  END SUBROUTINE


  SUBROUTINE StoreOptimal (NZCoeffs,NZCoeffsDot)
  IMPLICIT NONE
  DOUBLE COMPLEX, INTENT(IN) :: NZCoeffs(N_orbitals, NoOfPW)
  DOUBLE COMPLEX, INTENT(IN) :: NZCoeffsDot(N_orbitals, NoOfPW)

  OPEN ( 10, FILE='opt.dat')
  WRITE (10,*) NZCoeffs
  WRITE (10,*) NZCoeffsDot
  CLOSE( 10 )
  END SUBROUTINE

  SUBROUTINE GetOptimal (NZCoeffs,NZCoeffsDot)
  IMPLICIT NONE
  LOGICAL       :: Exists
  DOUBLE COMPLEX, INTENT(OUT) :: NZCoeffs(N_orbitals, NoOfPW)
  DOUBLE COMPLEX, INTENT(OUT) :: NZCoeffsDot(N_orbitals, NoOfPW)
  
  INQUIRE (FILE='opt.dat', EXIST=Exists)
  IF (Exists) THEN
    OPEN ( 10, FILE='opt.dat')
    READ (10,*) NZCoeffs
    READ (10,*) NZCoeffsDot
    CLOSE(10)
  END IF
  END SUBROUTINE

  
  SUBROUTINE StoreLast (NZCoeffs,NZCoeffsDot, R_IonDot, Time)
  IMPLICIT NONE
  DOUBLE COMPLEX, INTENT(IN) :: NZCoeffs(N_orbitals, NoOfPW)
  DOUBLE COMPLEX, INTENT(IN) :: NZCoeffsDot(N_orbitals, NoOfPW)
  DOUBLE COMPLEX, INTENT(IN) :: R_IonDot(N_ion, 3)
  DOUBLE PRECISION, INTENT(IN) :: Time
  
  OPEN ( 10, FILE='last.dat')
  WRITE (10,*) Time
  WRITE (10,*) NZCoeffs
  WRITE (10,*) NZCoeffsDot
  WRITE (10,*) Ions
  WRITE (10,*) R_IonDot
  CLOSE( 10 )
  END SUBROUTINE

  SUBROUTINE GetLast (NZCoeffs,NZCoeffsDot, R_IonDot, Time)
  IMPLICIT NONE
  DOUBLE COMPLEX, INTENT(OUT) :: NZCoeffs(N_orbitals, NoOfPW)
  DOUBLE COMPLEX, INTENT(OUT) :: NZCoeffsDot(N_orbitals, NoOfPW)
  DOUBLE COMPLEX, INTENT(OUT) :: R_IonDot(N_ion, 3)
  DOUBLE PRECISION, INTENT(OUT) :: Time
         
  LOGICAL       :: Exists
  INQUIRE(FILE='last.dat', EXIST=Exists)
  IF (Exists) THEN
    OPEN ( 10, FILE='last.dat')
    READ (10,*) Time
    READ (10,*) NZCoeffs
    READ (10,*) NZCoeffsDot
    READ (10,*) Ions
    READ (10,*) R_IonDot
    CLOSE(10)
  END IF 
  END SUBROUTINE GetLast  
  
  SUBROUTINE PeriodicBoundary()
  IMPLICIT NONE
  INTEGER       :: N, I 
     DO N=1, N_Ion
        DO I=1, 3
          Ions(N)%R_I(I) =  Ions(N)%R_I(I) + BoxL/2 - & 
              BoxL*DNINT(Ions(N)%R_I(I)/BoxL) 
        END DO
     END DO
  END SUBROUTINE

END MODULE Utilities

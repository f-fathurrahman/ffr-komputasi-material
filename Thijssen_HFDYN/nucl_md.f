      PROGRAM Twoatom
C*****************************************************************
C* This program performs a Quantum Molecular Dynamics (MD)       *
C* calculation of the ground state of the hydrogen molecule with *
C* the nuclear motion included in the MD algorithm.              *
C*                                                               *
C* Program described in "Computational Physics", J. M. Thijssen  *
C* Cambridge University Press, 1999                              *
C* Section 9.3.2                                                 *
C*                                                               *
C* Program written by Jos Thijssen 1997-1999                     *
C*****************************************************************

      include "nucl.glob"

      CALL Initialise 
      CALL Verlet

      END

      SUBROUTINE Initialise
      include "nucl.glob"
      DATA Alpha/13.00773D0,1.962079D0,0.444529D0,0.1219492D0,
     .           13.00773D0,1.962079D0,0.444529D0,0.1219492D0/

      WRITE (*,*) 'Geef Dist'
      READ(5,*) Dist

      OPEN (3, File='distances')

      PI = 4.D0*atan(1.D0)

      CALL FillExpArray
      CALL Overlap
      CALL DerivOverlap
C This factor is stored for efficiency:
      FFac = 0.5D0*SQRT(PI)
      CALL CalcHamilton
      CALL BuildSuper

      END



      SUBROUTINE Verlet

      include "nucl.glob"

      DOUBLE PRECISION Force(NSize), Coeffs(NSize), OldCoeffs(NSize),
     .       h, Ener, Gamma, NewCoeffs(NSize),
     .       NucForce, TotNucForce, NucTimeStep, Lambda, 
     .       NucMass, OldDist, NewDist, GHDiv2
      PARAMETER (NucMass = 1837.D0)

      INTEGER R, Count, StepNum, NucStepNum

      PRINT *, 'Give time step h'
      READ *, h
      PRINT *, 'Give friction coefficient Gamma'
      READ *, Gamma
      PRINT *, 'Give number of time steps'
      READ *, StepNum

C The factor GHDiv2 occurs in the verlet algorithm with friction
      GHDiv2 = Gamma*h*0.5D0
C NucStepNum is the number of electronic integration steps between
C two nuclear integration steps
      NucStepNum = NINT(SQRT(NucMass))
      NucTimeStep = NucStepNum*h
C Initial values for the coefficients
      DO R=1, NSize
        Coeffs(R) = 1.D0
        OldCoeffs(R) = 1.D0
      END DO
C Initial value for the Lgrange parameter Lambda
      Lambda = 0.0D0

C Standard normalisation of the coefficients C(R)
      CALL StdNormalise (OldCoeffs)
      CALL StdNormalise (Coeffs)
      CALL UpdateNewC(Coeffs)

C For nuclear motion:
      TotNucForce = 0.D0
      OldDist = Dist
    
      DO Count=1, StepNum
        CALL CalcNucForce(NucForce, Coeffs, Lambda)
C The average of the nuclear force during the (long) nuclear
C time steps is calculated
        TotNucForce = TotNucForce+NucForce*h

        CALL CalcCForce(Force, Coeffs)
C Verlet step:
        DO R=1, NSize
C Friction Verlet algorithm
          NewCoeffs(R) = (2*Coeffs(R)-(1-GHDiv2)*OldCoeffs(R)
     .                  +h*h*(-Force(R)))/(1+GHDiv2)
          OldCoeffs(R) = Coeffs(R)
        END DO
        CALL Normalise (NewCoeffs, Coeffs, Lambda)
C The Lagrange parameter from the normalisation routine must
C be corrected because of the frictional Verlet algorithm and 
C because an extra factor h^2 occurs with the force. 
        Lambda = Lambda*(1+GHDiv2)/(h*h)
C Nuclear Motion, only every NucStepNum steps:
        IF (MOD(Count,NucStepNum) .EQ. 0) THEN
          NewDist = 2*Dist-OldDist-2.D0*NucTimeStep*TotNucForce/NucMass
C Note that a factor NucTimeStep is already included in TotNucForce
          OldDist = Dist
          Dist = NewDist
          CALL CalcCEner(Ener, Coeffs)
          CALL UpdateNewDist
C Some output
          print *, Ener+1/Dist, TotNucForce/NucTimeStep, Dist
          write(3, '(F10.5)') dist
          TotNucForce = 0.D0
        END IF
C Bookkeeping and results for electronic coefficients:
        CALL CalcCEner(Ener, Coeffs)
        CALL UpdateNewC(NewCoeffs)
        DO R=1, NSize
          Coeffs(R) = NewCoeffs(R)
        END DO
      END DO
      END


      SUBROUTINE UpdateNewDist

C After the distance has changed, calculate all electronic 
C matrix elements depending on the distance

      CALL FillExpArray
      CALL Overlap
      CALL DerivOverlap
      CALL CalcHamilton
      CALL BuildSuper
      END


      SUBROUTINE UpdateNewC(Coeffs)

      include "nucl.glob"

      DOUBLE PRECISION Coeffs(NSize)

C After the Coeffs have changed, recalculate matrices depending 
C on those Coeffs

      CALL BuildDensMat(Coeffs)
      CALL BuildG
      CALL CalcFock
      END

      




      SUBROUTINE Normalise (NewCoeffs, Coeffs, Lambda)

C Normalisation of Coeffs according to Eq. (9.32)

      include "nucl.glob"

      DOUBLE PRECISION Coeffs(NSize), NewCoeffs(NSize), Norm, 
     .       SCTild(NSize), SC(NSize),
     .       SSC(NSize), A, B, CC, Discr, Lambda2, Lambda

      INTEGER R, S

      DO R=1, NSize
        SCTild(R) = 0.D0
        SC(R) = 0.D0
        DO S=1, NSize
          SCTild(R) = SCTild(R)+SMatrix(R, S)*NewCoeffs(S)
          SC(R) = SC(R)+SMatrix(R, S)*Coeffs(S)
        END DO
      END DO
      DO R=1, NSize
        SSC(R) = 0.D0
        DO S=1, NSize
          SSC(R) = SSC(R)+SMatrix(R, S)*SC(S)
        END DO
      END DO
      A = 0.D0
      B = 0.D0
      CC = 0.D0
      DO R=1, NSize
        CC = CC + NewCoeffs(R)*SCTild(R)
        A = A + SSC(R)*SC(R)
        B = B - 2*SC(R)*SCTild(R)
      END DO
      CC = CC-1.D0
      IF (B*B.LT.4*A*CC) print *, 'discr neg'
      Discr = SQRT(B*B-4*A*CC)
      Lambda2 = 0.5D0*(-B-Discr)/A
      DO R=1, NSize
        NewCoeffs(R) = NewCoeffs(R)-Lambda2*SC(R)
      END DO
      Norm = 0.D0
      DO R=1, NSize
        DO S=1, NSize
          Norm = Norm + NewCoeffs(R)*SMatrix(R,S)*NewCoeffs(S)
        END DO
      END DO
      Lambda = Lambda2
      END




      SUBROUTINE CalcCForce(Force, Coeffs)
C Calculation of "forces" acting on the coefficients,
C see Eq. (9.29)

      include "nucl.glob"

      DOUBLE PRECISION Force(NSize), Coeffs(NSize)

      INTEGER R, S

      DO R=1, NSize
        Force(R) = 0.D0
      END DO
      DO R=1, NSize
        DO S = 1, R-1
          FMatrix(S,R) = FMatrix(R, S)
        END DO
      END DO
      DO R=1, NSize
        DO S = 1, NSize
          Force(R) = Force(R) + 4*FMatrix(R,S)*Coeffs(S)
        END DO
      END DO
      END



      SUBROUTINE CalcNucForce(NucForce, Coeffs, Lambda)

C Calculation of the nuclear forces using Eq. (4.40) and those
C given in section 9.3.2

      include "nucl.glob"

      DOUBLE PRECISION NucForce, Coeffs(NSize), F1, Lambda

      INTEGER R, S, T, U

      NucForce = 0.D0
C Contribution due to the 1-electron Hamiltonian
      DO R=1, NSize
        DO S = 1, NSize
          NucForce = NucForce + DHamilton(R, S)*Coeffs(R)*Coeffs(S)
        ENDDO
      ENDDO
C Contribution due to the two-electron matrix elements
      DO R=1, NSize
        DO S = 1, NSize
          DO T = 1, NSize
            DO U = 1, NSize
              NucForce = NucForce + DQMat(R,S,T,U)*
     .                 Coeffs(R)*Coeffs(S)*Coeffs(T)*Coeffs(U)*0.5D0
            ENDDO
          ENDDO
        ENDDO
      ENDDO
C Contribution of the electrostatic repulsion force between the nuclei
      NucForce = 2*NucForce-1/Dist/Dist
C Contribution due to the Lagrange parameter (2nd term in Eq. 9.21):
      F1 = 0.D0
      DO R=1, N
        DO S = 1, N
          F1 = F1 - Lambda*DSMat(R+N,S)*Coeffs(R)*Coeffs(S)
          F1 = F1 - Lambda*DSMat(R,S+N)*Coeffs(R)*Coeffs(S)
        ENDDO
      ENDDO
      NucForce = NucForce - 0.5D0*F1
 
      END



      SUBROUTINE CalcCEner(Ener, Coeffs)
C Calculation of the energy, see Eqs. (4.40) and (4.46)

      include "nucl.glob"

      DOUBLE PRECISION Ener, Coeffs(NSize)

      INTEGER R, S

      Ener = 0.d0
      DO R=1, NSize
        DO S = 1, R-1
          Ener = Ener+2*Coeffs(R)*Coeffs(S)*(FMatrix(R,S)+
     .                             Hamilton(R,S))
        END DO
        Ener = Ener+Coeffs(R)*Coeffs(R)*
     .          (FMatrix(R,R)+Hamilton(R,R))
      END DO

      END






      SUBROUTINE FillExpArray
C The array ExpArray is used in the calculations of the overlap
C matrix etcetera

      include "nucl.glob"

      INTEGER R, S

      DO R=1,N
         DO S=1,N
            ExpArray (R,S+N) =
     .       DEXP(-Alpha(R)*Alpha(S+N)/(Alpha(R)+Alpha(S+N))*Dist*Dist)
            ExpArray(S+N,R) = ExpArray(R,S+N)
         END DO
      END DO
      END




      SUBROUTINE Overlap
C Overlap matrix

      include "nucl.glob"

      INTEGER R, S

      DOUBLE PRECISION TussenResult, Factor

      TussenResult = .0D0
      DO R=1,N
         DO S=1,R-1
            Factor = Pi/(Alpha(R)+Alpha(S))
            SMatrix(R,S) = Factor*DSQRT(Factor)
            SMatrix(S,R) = SMatrix(R,S)

            Factor     = Pi/(Alpha(R+N)+Alpha(S+N))
            SMatrix(R+N,S+N) = Factor*DSQRT(Factor)
            SMatrix(S+N,R+N) = SMatrix(R+N,S+N)
         END DO
         Factor    = Pi/(Alpha(R)+Alpha(R))
         SMatrix(R,R)    =  Factor*DSQRT(Factor)

         Factor      = Pi/(Alpha(R+N)+Alpha(R+N))
         SMatrix(R+N,R+N)  =  Factor*DSQRT(Factor)
      END DO
      DO R=1,N
         DO S=1,N
            Factor    = Pi/(Alpha(R)+Alpha(S+N))
            SMatrix(R, S+N) = Factor*DSQRT(Factor)*ExpArray(R,S+N)
            SMatrix(S+N, R) = SMatrix(R,S+N)
         END DO
      END DO
      END


      SUBROUTINE DerivOverlap

      include "nucl.glob"

      INTEGER R, S

      DOUBLE PRECISION Alph

      DO R=1,N
        DO S=1,R-1
          DSMat(R,S) = 0.D0
          DSMat(S,R) = 0.D0

          DSMat(R+N,S+N) = 0.D0
          DSMat(S+N,R+N) = 0.D0
        ENDDO
        DSMat(R,R)    =  0.D0
      ENDDO
      DO R=1,N
        DO S=1,N
          Alph = Alpha(R)+Alpha(S)
          DSMat(R+N, S) = -2*Alpha(R)*Alpha(S)/Alph*
     .                     Dist*SMatrix(R+N,S)
          DSMat(R, S+N) = DSMat(R+N, S)
        ENDDO
      ENDDO
      END



      SUBROUTINE CalcHamilton
C The Hamiltonian and its derivative DHamilton w.r.t. the nuclear separation 
C "Dist" is calculated

      include "nucl.glob"

      INTEGER R, S

      DOUBLE PRECISION Kinet(NSize,NSize), Coul(NSize,NSize),
     .                 DKinet(NSize,NSize), DCoul(NSize,NSize)

      CALL CalcKinet(Kinet, DKinet)
      CALL CalcCoul(Coul)
      CALL CalcDCoul(DCoul)
      DO R=1, N
        DO S=1, R
          Hamilton(R,S) = Kinet(R,S) - Coul(R, S)
          Hamilton(R+N, S+N) = Hamilton(R, S)
          Hamilton(R+N, S) = Kinet(R+N, S) - Coul(R+N, S)
          Hamilton(S+N, R) = Hamilton(R+N, S)
        END DO
      END DO

      DO R=1, N
        DO S=1, N
          DHamilton(R,S) = - DCoul(R, S)
          DHamilton(R+N, S+N) = DHamilton(R, S)
          DHamilton(R+N, S) = DKinet(R+N, S) - DCoul(R, S+N)
          DHamilton(R, S+N) = DHamilton(R+N, S)
        ENDDO
      ENDDO

      END




      SUBROUTINE CalcKinet (Kinet, DKinet)
C Calculates the kinetic energy matrix element and its derivative w.r.t.
C the separation "Dist" between the nuclei

      include "nucl.glob"

      DOUBLE PRECISION Kinet(NSize,NSize), K0, K1,
     .       A, B, DKinet(NSize,NSize), ApB

      INTEGER R, S
      
      DO R=1, N
        DO S=1, R
          A = Alpha(R)
          B = Alpha(S)
          ApB = A + B
          K0 = 3*A*B/ApB
          K1 = K0-2*(A*B/ApB*Dist)**2
          Kinet(R,S) = K0*SMatrix(R,S)
          Kinet(R+N,S) = K1*SMatrix(R+N,S)
          DKinet(R,S) = 0.D0
          DKinet(S,R) = 0.D0
          DKinet(R+N,S+N) = 0.D0
          DKinet(S+N,R+N) = 0.D0
          DKinet(R+N, S) = K1*DSmat(R+N,S) - 
     .                       4*SMatrix(R+N,S)*(A*B/ApB)**2*Dist
          DKinet(S+N, R) = DKinet(R+N,S)
          DKinet(R, S+N) = DKinet(R+N,S)
          DKinet(S, R+N) = DKinet(R+N,S)
        ENDDO
      ENDDO
      END      




      DOUBLE PRECISION FUNCTION F (X)
C The function F_0(x) defined  in Eq. (4.112)

      include "nucl.glob"
      DOUBLE PRECISION X

      IF (X.EQ.0.D0) THEN
        F = 1.d0
      ELSE
        F = FFac/SQRT(X)*erf(SQRT(x))
      END IF
      END


      DOUBLE PRECISION FUNCTION F1 (X)

      include "nucl.glob"
      DOUBLE PRECISION X, F

      IF (X.EQ.0.D0) THEN
        F1 = 0.3333333333333333333d0
      ELSE
        F1 = 0.5D0*(F(X)-exp(-X))/X
      ENDIF
      END



      SUBROUTINE CalcCoul (Coul)
C Matrix element of the Coulomb energy
      include "nucl.glob"

      DOUBLE PRECISION Coul(NSize,NSize), L0, L1, 
     .       A, B, Theta, F, T

      INTEGER R, S
      
      DO R=1, N
        DO S=1, R
          A = Alpha(R)
          B = Alpha(S)
          Theta = 2*SQRT((A+B)/PI)
          t = (A+B)*Dist**2
          L0 = F(t) + 1.d0
          t = Dist**2/(A+B)
          L1 = F(A**2*T) + F(B**2*T)
          Coul(R,S) = Theta*L0*SMatrix(R,S)
          Coul(R+N,S) = Theta*L1*SMatrix(R+N,S)
        END DO
      END DO
      END      



      SUBROUTINE CalcDCoul (DCoul)
C Matrix element of the derivative of the Coulomb energy w.r.t. nuclear
C separation "Dist"
      include "nucl.glob"

      DOUBLE PRECISION DCoul(NSize,NSize), L0, LI, 
     .       A, B, Theta, F, T, F1, F1T, F1A, F1B

      INTEGER R, S
      
      DO R=1, N
        DO S=1, N
          A = Alpha(R)
          B = Alpha(S)
          Theta = 2*SQRT((A+B)/PI)
          t = (A+B)*Dist**2
          L0 = F(t) + 1.d0
          F1T = F1(T)
          LI = Dist*F1T
          DCoul(R, S) = Theta*(-SMatrix(R,S)*2*(A+B)*Dist*F1T)
          DCoul(R+N, S+N) = DCoul(R, S)
          t = Dist**2/(A+B)
          F1A = F1(A**2*T)
          F1B = F1(B**2*T)
          L0 = F(A**2*T) + F(B**2*T)
          LI = 2*Dist*(A*A*F1A + B*B*F1B)/(A+B)
          DCoul(R+N,S) = Theta*(L0*DSMat(R+N,S)-
     .                          LI*SMatrix(R+N,S))
          DCoul(R, N+S) = DCoul(R+N, S)
        ENDDO
      ENDDO
      END      




      SUBROUTINE BuildDensMat(Coeffs)
C The density matrix is defined in Eq. (4.66)
C It occurs in the construction of the Fock matrix

      include "nucl.glob"

      DOUBLE PRECISION Coeffs(NSize)

      INTEGER T, U

      DO T=1,NSize
         DO U=1,T-1
            DensMat(T, U) = 2.D0*Coeffs(T)*Coeffs(U)
            DensMat(U, T) = DensMat(T,U)
         END DO
         DensMat(T,T) = 2.D0*Coeffs(T)*Coeffs(T)
      END DO
      END




      SUBROUTINE StdNormalise(Coeffs)
C This is a "standard" normalisation: it does not 
C solve for a Lagrange parameter Lambda which must be
C chosen such as to ensure proper normalisation, as 
C is done in the routine "Normalise" above.
      include "nucl.glob"

      DOUBLE PRECISION Coeffs(NSize), Norm

      INTEGER R, S

      Norm = 0.D0
      DO R=1, NSize
        DO S=1, NSize
          Norm = Norm + Coeffs(R)*Smatrix(R,S)*Coeffs(S)
        END DO
      END DO

      Norm = SQRT(Norm)

      DO R=1, NSize
        Coeffs(R) = Coeffs(R)/Norm
      END DO
      END





      SUBROUTINE BuildG
C Builds the matrix G which is a contraction of the two-electron
C matrix elements with the density matrix, see Eq. (4.86)
C Full use is made of symmetries

      include "nucl.glob"

      DOUBLE PRECISION SuperElem


      INTEGER R, S, T, U, MaxU, Count

      DO R=1,NSize
         DO S=1,R
            GMatrix(R,S) = 0.D0
         END DO
      END DO

      Count = 0
      DO R=1,NSize
        DO S=1,R
          DO T=1,R
            IF (T .LT. R) THEN
              MaxU = T
            ELSE
              MaxU = S
            END IF
            DO U=1,MaxU
              Count = Count + 1
              SuperElem = QMatrix(R,S,T,U)
              IF (TypeMatrix(R,S,T,U).EQ.1) THEN
                GMatrix(R,R) = GMatrix(R,R)+0.5*SuperElem*
     .                                              DensMat(R,R)
              ELSE IF (TypeMatrix(R,S,T,U).EQ.2) THEN
                GMatrix(R,R) = GMatrix(R,R)+0.5D0*SuperElem*
     .                                              DensMat(T,T)
                GMatrix(T,T) = GMatrix(T,T)+0.5D0*SuperElem*
     .                                              DensMat(R,R)
              ELSE IF (TypeMatrix(R,S,T,U).EQ.3) THEN
                GMatrix(R,R) = GMatrix(R,R)+SuperElem*
     .                                              DensMat(T,U)
                GMatrix(T,U) = GMatrix(T,U)+0.5D0*SuperElem*
     .                                              DensMat(R,R)
              ELSE IF (TypeMatrix(R,S,T,U).EQ.4) THEN
                GMatrix(T,T) = GMatrix(T,T)+SuperElem*
     .                                              DensMat(R,S)
                GMatrix(R,S) = GMatrix(R,S)+0.5D0*SuperElem*
     .                                              DensMat(T,T)
              ELSE IF (TypeMatrix(R,S,T,U).EQ.5) THEN
                GMatrix(R,S) = GMatrix(R,S)+SuperElem*
     .                                              DensMat(R,S)
              ELSE IF (TypeMatrix(R,S,T,U).EQ.6) THEN
                GMatrix(R,S) = GMatrix(R,S)+SuperElem*
     .                                              DensMat(T,U)
                GMatrix(T,U) = GMatrix(T,U)+SuperElem*
     .                                              DensMat(R,S)
              END IF
            END DO
          END DO
        END DO
      END DO   
      END



      SUBROUTINE CalcFock
C The fock matrix is evaluated, see Eq. (4.88)

      include "nucl.glob"

      INTEGER  R, S

      DO R = 1, NSize
        DO S = 1, R
           FMatrix(R,S) = Hamilton(R,S) + GMatrix(R,S)
        END DO
      END DO
      END





      SUBROUTINE BuildSuper
C The Supermatrix Q(r,s,t,u) contains the two-electron
C matrix elements. From a particular set of indices, the 
C full matrix is built. 

      include "nucl.glob"

      INTEGER DetermineType, R, S, T, U, MaxU
      DOUBLE PRECISION A

      CALL FillQMatrix
      CALL FillDQ
      DO R=1,NSize
        DO S=1,R
          DO T=1,R
            IF (T .LT. R) THEN
              MaxU = T
            ELSE
              MaxU = S
            END IF
            DO U=1,MaxU
              TypeMatrix (R,S,T,U) = DetermineType (R,S,T,U)
              A = QMatrix(R, S, T, U)
              QMatrix(S,R,T,U) = A
              QMatrix(R,S,U,T) = A
              QMatrix(S,R,U,T) = A
              QMatrix(T,U,R,S) = A
              QMatrix(U,T,R,S) = A
              QMatrix(T,U,S,R) = A
              QMatrix(U,T,S,R) = A
            END DO
          END DO
        END DO
      END DO
      END



      INTEGER FUNCTION DetermineType(R,S,T,U)
C Determines relative values of coeffs rstu. See the
C discussion on page 75.

      IMPLICIT NONE 

      INTEGER R, S, T, U, Type

      TYPE = 0
      IF ((R .EQ. S) .AND. (T .EQ. U) .AND. (R.EQ.T)) THEN
        Type = 1
      ELSEIF ((R .EQ. S) .AND. (T .EQ. U)) THEN
        Type = 2
      ELSEIF ((R .EQ. S) .AND. (T .GT. U)) THEN
        Type = 3
      ELSEIF ((R .GT. S) .AND. (T .EQ. U)) THEN
        Type = 4
      ELSEIF ((R .GT. S) .AND. (T .GT. U) .AND. 
     .        (R.EQ.T) .AND. (S.EQ.U))THEN
        Type = 5
      ELSE
        Type = 6
      END IF

      DetermineType = Type

      END



      SUBROUTINE FillQMatrix
C The actual evaluation of the two-electron matrix elements, 
C see Eq. (4.121)
      include "nucl.glob"
      DOUBLE PRECISION A, B, CC, D, P1, P2, P3, G, Lambda, F, tt, P, Q
      INTEGER R, S, T, U, MaxU

      DO R=1, NSize
        DO S=1, R
          DO T=1, R
            IF (T .LT. R) THEN
              MaxU = T
            ELSE
              MaxU = S
            END IF
            DO U=1,MaxU
              A = Alpha(R) 
              B = Alpha(S) 
              CC = Alpha(T) 
              D = Alpha(U) 
              P1 = A+B
              P2 = CC+D
              P3 = P1+P2
              P = ((2*((R-1)/4)-1)*A+(2*((S-1)/4)-1)*B)/(A+B)*Dist/2
              Q = ((2*((T-1)/4)-1)*CC+(2*((U-1)/4)-1)*D)/(CC+D)*Dist/2
              tt = P1*P2/P3*(P-Q)**2
              Lambda = 2*SQRT(P1*P2/PI/P3)
              G = F(tt)
              QMatrix(R, S, T, U) = SMatrix(R,S)*SMatrix(T,U)*Lambda*G
            END DO
          END DO
        END DO
      END DO
      END



      SUBROUTINE FillDQ
      include "nucl.glob"
      DOUBLE PRECISION A, B, CC, D, P1, P2, P3, G, Lambda, F, tt, 
     .                 P, Q, F1, GI
      INTEGER R, S, T, U

      DO R=1, NSize
        DO S=1, NSize
          DO T=1, NSize
            DO U=1, NSize
              A = Alpha(R) 
              B = Alpha(S) 
              CC = Alpha(T) 
              D = Alpha(U) 
              P1 = A+B
              P2 = CC+D
              P3 = P1+P2
              P = ((2*((R-1)/4)-1)*A+(2*((S-1)/4)-1)*B)/(A+B)
              Q = ((2*((T-1)/4)-1)*CC+(2*((U-1)/4)-1)*D)/(CC+D)
              tt = P1*P2/P3*((P-Q)*Dist/2)**2
              Lambda = 2*SQRT(P1*P2/PI/P3)
              G = F(tt)
              GI = P1*P2/P3*F1(tt)*(P-Q)**2*Dist/2
              DQMat(R, S, T, U) = Lambda*(DSMat(R,S)*SMatrix(T,U)*G+
     .                 DSMat(T,U)*SMatrix(R,S)*G-
     .                 SMatrix(R,S)*SMatrix(T,U)*GI)
            ENDDO
          ENDDO
        ENDDO
      ENDDO
      END











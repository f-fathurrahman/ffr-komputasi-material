      PROGRAM TwoAtom
C*****************************************************************
C* This program performs a Quantum Molecular Dynamics (MD)       *
C* calculation of the ground state of the hydrogen molecule with *
C* fixed distance between the nuclei.                            *
C*                                                               *
C* Program described in "Computational Physics", J. M. Thijssen  *
C* Cambridge University Press, 1999                              *
C* Second edition: 2007                                          *
C* Section 9.3.1                                                 *
C*                                                               *
C* Program written by Jos Thijssen 1997-1999                     *
C* Small additions made in 2007                                  *
C*****************************************************************

      include "elec.glob"

      CALL Initialise 
      CALL Verlet

      END

      SUBROUTINE Initialise
      include "elec.glob"
      DATA Alpha/13.00773D0,1.962079D0,0.444529D0,0.1219492D0,
     .           13.00773D0,1.962079D0,0.444529D0,0.1219492D0/

      WRITE (*,*) 'Geef Dist'
      READ(5,*) Dist

      PI = 4.D0*atan(1.D0)

      Call FillExpArray
      Call Overlap
C This factor is stored for efficiency:
      FFac = 0.5D0*SQRT(PI)
      Call CalcHamilton
      Call BuildSuper

      END



      SUBROUTINE Verlet

      include "elec.glob"

      DOUBLE PRECISION Force(NSize), Coeffs(NSize), OldCoeffs(NSize),
     .       h, Ener, Gamma, NewCoeffs(NSize), GHDiv2

      INTEGER R, Count, StepNum

      PRINT *, 'Give time step h'
      READ *, h
      PRINT *, 'Give friction coefficient Gamma'
      READ *, Gamma
      PRINT *, 'Give number of time steps'
      READ *, StepNum

      GHDiv2 = Gamma*h*0.5D0
C Initial values for the coefficients
      DO R=1, NSize
        Coeffs(R) = 1.D0
        OldCoeffs(R) = 1.D0
      END DO

      CALL StdNormalise (OldCoeffs)
      CALL StdNormalise (Coeffs)
      CALL UpdateNewC(Coeffs)

      DO Count=1, StepNum
        CALL CalcCForce(Force, Coeffs)
C Verlet step:
        DO R=1, NSize
          NewCoeffs(R) = (2*Coeffs(R)-(1-Gamma)*OldCoeffs(R)
     .                  +h*h*(-Force(R)))/(1+Gamma)
          OldCoeffs(R) = Coeffs(R)
        END DO
C Bookkeeping and results:
        CALL CalcCEner(Ener, Coeffs)
        print *, 'Energy  ', Ener+1.D0/Dist
        CALL Normalise (NewCoeffs, Coeffs, h)
        CALL UpdateNewC(NewCoeffs)
        DO R=1, NSize
          Coeffs(R) = NewCoeffs(R)
        END DO
      END DO
      END



      SUBROUTINE UpdateNewC(Coeffs)

      include "elec.glob"

      DOUBLE PRECISION Coeffs(NSize)

C After the Coeffs have changed, recalculate matrices depending 
C on those Coeffs

      CALL BuildDensMat(Coeffs)
      CALL BuildG
      CALL CalcFock
      END

      




      SUBROUTINE Normalise (NewCoeffs, Coeffs, h)

C Normalisation of Coeffs according to Eq. (9.32)

      include "elec.glob"

      DOUBLE PRECISION Coeffs(NSize), NewCoeffs(NSize), Norm, 
     .       h, SCTild(NSize), SC(NSize),
     .       SSC(NSize), A, B, CC, Discr, Lambda2

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
        A = A + h*h*SSC(R)*SC(R)
        B = B - 2*h*SC(R)*SCTild(R)
      END DO
      CC = CC-1.D0
      IF (B*B.LT.4*A*CC) print *, 'discr neg'
      Discr = SQRT(B*B-4*A*CC)
      Lambda2 = 0.5D0*(-B-Discr)/A
      DO R=1, NSize
        NewCoeffs(R) = NewCoeffs(R)-h*Lambda2*SC(R)
      END DO
      Norm = 0.D0
      DO R=1, NSize
        DO S=1, NSize
          Norm = Norm + NewCoeffs(R)*SMatrix(R,S)*NewCoeffs(S)
        END DO
      END DO
      END




      SUBROUTINE CalcCForce(Force, Coeffs)
C Calculation of "forces" acting on the coefficients,
C see Eq. (9.29)

      include "elec.glob"

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





      SUBROUTINE CalcCEner(Ener, Coeffs)
C Calculation of the energy, see Eqs. (4.40) and (4.46)

      include "elec.glob"

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

      include "elec.glob"

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

      include "elec.glob"

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





      SUBROUTINE CalcHamilton
C The Hamiltonian is calculated

      include "elec.glob"

      INTEGER R, S

      DOUBLE PRECISION Kinet(NSize,NSize), Coul(NSize,NSize)

      CALL CalcKinet(Kinet)
      CALL CalcCoul(Coul)
      DO R=1, N
        DO S=1, R
          Hamilton(R,S) = Kinet(R,S) - Coul(R, S)
          Hamilton(R+N, S+N) = Hamilton(R, S)
          Hamilton(R+N, S) = Kinet(R+N, S) - Coul(R+N, S)
          Hamilton(S+N, R) = Hamilton(R+N, S)
        END DO
      END DO
      END



      SUBROUTINE CalcKinet (Kinet)
C Kinetic energy matrix element
      include "elec.glob"

      DOUBLE PRECISION Kinet(NSize,NSize), K0, K1,
     .                 A, B

      INTEGER R, S
      
      DO R=1, N
        DO S=1, R
          A = Alpha(R)
          B = Alpha(S)
          K0 = 3*A*B/(A+B)
          K1 = K0-2*((A*B)/(A+B)*Dist)**2
          Kinet(R,S) = K0*SMatrix(R,S)
          Kinet(R+N,S) = K1*SMatrix(R+N,S)
        END DO
      END DO
      END




      DOUBLE PRECISION FUNCTION F (X)
C The function F_0(x) defined  in Eq. (4.112)

      include "elec.glob"
      DOUBLE PRECISION X, erf

      IF (X.EQ.0.D0) THEN
        F = 1.d0
      ELSE
        F = FFac/SQRT(X)*erf(SQRT(x))
      END IF
      END





      SUBROUTINE CalcCoul (Coul)
C Matrix element of the Coulomb energy
      include "elec.glob"

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






      SUBROUTINE BuildDensMat(Coeffs)
C The density matrix is defined in Eq. (4.66)
C It occurs in the construction of the Fock matrix

      include "elec.glob"

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
      include "elec.glob"

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

      include "elec.glob"

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
C The fock matrix is evaluated, see (4.88)

      include "elec.glob"

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

      include "elec.glob"

      INTEGER DetermineType, R, S, T, U, N2, MaxU
      DOUBLE PRECISION A

      N2 = 2*N

      CALL FillQMatrix
      DO R=1,N2
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

      IF (Type .EQ. 0) THEN
         WRITE (*,*) R, S, T, U
      END IF

      DetermineType = Type

      END



      SUBROUTINE FillQMatrix
C The actual evaluation of the two-electron matrix elements, 
C see Eq. (4.121)
      include "elec.glob"
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









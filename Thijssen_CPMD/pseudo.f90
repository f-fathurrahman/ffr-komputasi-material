MODULE pseudo
! This module contains the routines for evaluating the different
! terms of the pseudopotential. The parametrisation of
! Goedecker, Teter and Hutter is used (PRB 54, p1703 (1996))

USE globals
USE utilities
USE excorr

CONTAINS

CHARACTER(LEN=8) FUNCTION GetPPFilename(AtomNum)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: AtomNum
  CHARACTER(LEN=3) StrAtomNum

  IF (atomnum<10) THEN
     write(StrAtomNum,'(I1)') atomnum
  ELSE IF (atomnum<100) THEN
     write(StrAtomNum,'(I2)') atomnum
  ELSE
     write(StrAtomNum,'(I3)') atomnum
  END IF

  GetPPFilename = "PP//PP"//StrAtomNum

END FUNCTION GetPPFilename

INTEGER FUNCTION GetIndexPP(AtomNum)
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: AtomNum
  INTEGER :: I
  GetIndexPP=-1
  DO I=1, No_OF_DIFF_IONS
    IF (PP_Params(I)%AtomNum.EQ.AtomNum) THEN
        GetIndexPP=I
        EXIT
    END IF
  END DO
END FUNCTION GetIndexPP

SUBROUTINE Init_PP_Params (AtomNum)

  IMPLICIT NONE
  INTEGER, SAVE ::  NoPP
  INTEGER, INTENT(IN) :: AtomNum
  INTEGER :: ind

  CHARACTER(LEN=8) :: filename
  ind = GetIndexPP(AtomNum)
  IF (ind.EQ.-1) THEN
    NoPP=NoPP+1
    filename = GetPPFilename(AtomNum)
    print *, "Reading PP_Params from file "//filename
    PP_Params(NoPP)%AtomNum = AtomNum
    OPEN (9, File=filename)
    READ (9, *) PP_Params(NoPP)%Zion
    READ (9, *) PP_Params(NoPP)%N_nonzero_C_i
    READ (9, *) PP_Params(NoPP)%Xi
    READ (9, *) PP_Params(NoPP)%C(1)
    READ (9, *) PP_Params(NoPP)%C(2)
    READ (9, *) PP_Params(NoPP)%C(3)
    READ (9, *) PP_Params(NoPP)%C(4)
    READ (9, *) PP_Params(NoPP)%MaxL
    READ (9, *) PP_Params(NoPP)%r_s
    READ (9, *) PP_Params(NoPP)%h_1s
    READ (9, *) PP_Params(NoPP)%h_2s
    READ (9, *) PP_Params(NoPP)%r_p
    READ (9, *) PP_Params(NoPP)%h_1p
    CLOSE(9)
  END IF
END SUBROUTINE Init_PP_Params

DOUBLE PRECISION FUNCTION Local_PP(G2,AtomNum)
  ! Returns the long-and short range part of the FT of the local pseudopotential
  ! This is useful for calculating the Kohn-Sham Hamiltonian, but not for the 
  ! total energy, as the long range part should be treated differently in that case
  ! (Ewald sum)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: G2, AtomNum
  DOUBLE PRECISION :: Xi, C1, C2, Zion, CNum
  INTEGER :: ind

  ind = GetIndexPP(AtomNum)
  Zion = DBLE(PP_Params(ind)%Zion)
  Xi = PP_Params(ind)%Xi/BoxL

  IF (G2/=0) THEN
    Local_PP = -Zion/(BoxL*G2*PI)*EXP(-2*PI*PI*G2*Xi*Xi)+Short_PP(G2,AtomNum)
  ELSE
    Local_PP = 0.D0
    RETURN
  END IF
END FUNCTION Local_PP


DOUBLE PRECISION FUNCTION Short_PP(G2, AtomNum)
  ! Returns the long-and short range part of the FT of the local pseudopotential
  ! This is useful for calculating the Kohn-Sham Hamiltonian, but not for the 
  ! total energy, as the long range part should be treated differently in that case
  ! (Ewald sum)
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: G2, AtomNum
  DOUBLE PRECISION :: Xi, C1, C2
  INTEGER :: ind, CNum

  ind = GetIndexPP(AtomNum)
  Cnum = PP_Params(ind)%N_nonzero_C_i 
  Xi = PP_Params(ind)%Xi
  C1 = PP_Params(ind)%C(1)
  Xi = Xi/BoxL
  ! First, all prefactors of the common Gaussian exponent are calculated.

  Short_PP = (2.D0*PI)**1.5D0*Xi**3*C1
  IF (CNum >= 2) THEN
    C2 =  PP_Params(ind)%C(2)
    Short_PP =  Short_PP+(2.D0*PI)**1.5D0*Xi**3*C2*(3-4*PI*PI*G2*Xi*Xi)
  END IF
  Short_PP = Short_PP*EXP(-2*PI*PI*G2*Xi*Xi)
END FUNCTION Short_PP



DOUBLE PRECISION FUNCTION CoreDens(G2,AtomNum)
  ! Returns the long-and short range part of the FT of the local pseudopotential
  ! This is useful for calculating the Kohn-Sham Hamiltonian, but not for the 
  ! total energy, as the long range part should be treated differently in that case
  ! (Ewald sum)
  IMPLICIT NONE

  ! INTEGER, INTENT(IN) :: I, J, K, AtomNum
  INTEGER, INTENT(IN) :: G2, AtomNum
  DOUBLE PRECISION :: Xi, C1, C2, Zion
  INTEGER :: ind,  CNum

  ind = GetIndexPP(AtomNum)
  Zion = DBLE(PP_Params(ind)%Zion)
  Xi = PP_Params(ind)%Xi
  Xi = Xi/BoxL

  ! First, all prefactors of the common Gaussian exponent are calculated.

  CoreDens = -Zion*EXP(-2*PI*PI*G2*Xi*Xi)/Omega

END FUNCTION CoreDens



DOUBLE COMPLEX FUNCTION NonLoc (I, J, K, AtomNum, L, M, II)
! Returns the term P_\alpha^I from Eq. (178)
  ! Parameters: I, J, K: G-vector
  !             L, M : Angular Momentum
  !             II : Index of projector (for s-waves)

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: I, J, K, L, M, II, AtomNum;
  INTEGER :: G2
  DOUBLE PRECISION :: Pl, Xi, RInv
  DOUBLE COMPLEX :: Ylm
  INTEGER :: ind

  ind = GetIndexPP(AtomNum)
  Ylm = CMPLX(1/SQRT(4*PI))

  G2 = G2_Short(I,J,K)
  IF (L==0) THEN
    Xi = PP_Params(ind)%r_s/BoxL
    IF (II==1) THEN
      Pl = Xi**1.5D0*PI**1.25D0*4*SQRT(2.D0)*EXP(-2*PI*PI*G2*Xi*Xi)
    ELSE
      Pl = Xi**1.5D0*PI**1.25D0*8/SQRT(7.5D0)*EXP(-2*PI*PI*G2*Xi*Xi)*&
                (3-4*PI*PI*G2*Xi*Xi)
    END IF
  ELSE
    Xi = PP_Params(ind)%r_p/BoxL
    Pl = Xi**2.5D0*PI**1.25D0*8*2*PI*SQRT(G2/3.D0)*EXP(-2*PI*PI*G2*Xi*Xi)
    IF (G2 /=0) THEN
      RInv = 1/SQRT(DBLE(G2))
      SELECT CASE (M)
        CASE( 1); Ylm = SQRT(0.5D0)*Ylm*CMPLX(-DBLE(I), -DBLE(J))
        CASE( 0); Ylm = Ylm*CMPLX(DBLE(K), 0)
        CASE(-1); Ylm = SQRT(0.5D0)*Ylm*CMPLX(DBLE(I), -DBLE(J))
      END SELECT
      Ylm = SQRT(3.D0)*Ylm*RInv
    ELSE
      Ylm = CMPLX(0.D0)
    END IF   
  END IF
  NonLoc = Pl*Ylm
END FUNCTION NonLoc

END MODULE pseudo

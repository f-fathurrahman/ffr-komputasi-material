!-------------------------------------
real(8) FUNCTION short_pp(G2, AtomNum)
!-------------------------------------
  ! Returns the long-and short range part of the FT of the local pseudopotential
  ! This is useful for calculating the Kohn-Sham Hamiltonian, but not for the 
  ! total energy, as the long range part should be treated differently in that case
  ! (Ewald sum)
  
  use globals
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: G2, AtomNum
  real(8) :: Xi, C1, C2
  INTEGER :: ind, CNum
  INTEGER :: get_index_pp

  ind = get_index_pp(AtomNum)
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
END FUNCTION
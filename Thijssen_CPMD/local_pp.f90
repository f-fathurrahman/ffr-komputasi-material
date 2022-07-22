!------------------------------------
real(8) FUNCTION local_pp(G2,AtomNum)
!------------------------------------
  ! Returns the long-and short range part of the FT of the local pseudopotential
  ! This is useful for calculating the Kohn-Sham Hamiltonian, but not for the 
  ! total energy, as the long range part should be treated differently in that case
  ! (Ewald sum)
  use globals
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: G2, AtomNum
  real(8) :: Xi, C1, C2, Zion, CNum
  INTEGER :: ind
  INTEGER :: get_index_pp
  real(8) :: Short_PP

  ind = Get_Index_pp(AtomNum)
  Zion = DBLE(PP_Params(ind)%Zion)
  Xi = PP_Params(ind)%Xi/BoxL

  IF (G2/=0) THEN
    Local_PP = -Zion/(BoxL*G2*PI)*EXP(-2*PI*PI*G2*Xi*Xi)+Short_PP(G2,AtomNum)
  ELSE
    Local_PP = 0.D0
    RETURN
  END IF
END FUNCTION
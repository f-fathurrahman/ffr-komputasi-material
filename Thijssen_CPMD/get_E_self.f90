real(8) FUNCTION get_E_self(AtomNum)
  ! Returns the long-and short range part of the FT of the local pseudopotential
  ! This is useful for calculating the Kohn-Sham Hamiltonian, but not for the 
  ! total energy, as the long range part should be treated differently in that case
  ! (Ewald sum)
  use globals
  IMPLICIT NONE

  INTEGER, INTENT(IN) :: AtomNum
  real(8) :: Xi, Zion
  INTEGER :: ind
  integer :: get_index_pp

  ind = get_index_pp(AtomNum)
  Zion = DBLE(PP_Params(ind)%Zion)
  Xi = PP_Params(ind)%Xi

  Get_E_Self = Zion*Zion/(2*SQRT(PI)*Xi)

END FUNCTION 

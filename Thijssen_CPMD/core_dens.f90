!-------------------------------------
real(8) FUNCTION core_dens(G2,AtomNum)
!-------------------------------------
  use globals
  ! Returns the long-and short range part of the FT of the local pseudopotential
  ! This is useful for calculating the Kohn-Sham Hamiltonian, but not for the 
  ! total energy, as the long range part should be treated differently in that case
  ! (Ewald sum)
  IMPLICIT NONE

  ! INTEGER, INTENT(IN) :: I, J, K, AtomNum
  INTEGER, INTENT(IN) :: G2, AtomNum
  real(8) :: Xi, C1, C2, Zion
  INTEGER :: ind,  CNum
  ! Function
  integer :: get_index_pp


  ind = get_index_pp(AtomNum)
  Zion = DBLE(PP_Params(ind)%Zion)
  Xi = PP_Params(ind)%Xi
  Xi = Xi/BoxL

  ! First, all prefactors of the common Gaussian exponent are calculated.
  core_dens = -Zion*EXP(-2*PI*PI*G2*Xi*Xi)/Omega

END FUNCTION
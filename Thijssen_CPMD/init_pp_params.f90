!---------------------------------
SUBROUTINE init_pp_params(AtomNum)
!---------------------------------
  use globals
  IMPLICIT NONE
  INTEGER, SAVE :: NoPP
  INTEGER, INTENT(IN) :: AtomNum
  INTEGER :: ind
  ! Functions
  integer :: get_index_pp
  CHARACTER(LEN=8) :: get_pp_filename

  CHARACTER(LEN=8) :: filename

  ind = get_index_pp(AtomNum)
  IF( ind == -1) THEN
    NoPP = NoPP + 1
    filename = get_pp_filename(AtomNum)
    WRITE(*,*) 'Reading PP_Params from file '//filename
    PP_Params(NoPP)%AtomNum = AtomNum
    OPEN(9, File=filename)
    READ(9, *) PP_Params(NoPP)%Zion
    READ(9, *) PP_Params(NoPP)%N_nonzero_C_i
    READ(9, *) PP_Params(NoPP)%Xi
    READ(9, *) PP_Params(NoPP)%C(1)
    READ(9, *) PP_Params(NoPP)%C(2)
    READ(9, *) PP_Params(NoPP)%C(3)
    READ(9, *) PP_Params(NoPP)%C(4)
    READ(9, *) PP_Params(NoPP)%MaxL
    READ(9, *) PP_Params(NoPP)%r_s
    READ(9, *) PP_Params(NoPP)%h_1s
    READ(9, *) PP_Params(NoPP)%h_2s
    READ(9, *) PP_Params(NoPP)%r_p
    READ(9, *) PP_Params(NoPP)%h_1p
    CLOSE(9)
  ENDIF
END SUBROUTINE

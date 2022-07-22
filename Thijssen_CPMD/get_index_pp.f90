!-------------------------------------
INTEGER FUNCTION get_index_pp(AtomNum)
!-------------------------------------
  use globals, only: No_OF_DIFF_IONS, PP_Params
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: AtomNum
  INTEGER :: I
  get_index_pp=-1
  DO I=1, No_OF_DIFF_IONS
    IF (PP_Params(I)%AtomNum.EQ.AtomNum) THEN
        get_index_pp = I
        EXIT
    END IF
  END DO
END FUNCTION
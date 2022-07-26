!-----------------------------------------
FUNCTION get_index_pp(AtomNum) RESULT(res)
!-----------------------------------------
  USE globals, ONLY: No_OF_DIFF_IONS, PP_Params
  !
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: AtomNum
  INTEGER :: I
  INTEGER :: res
  
  res = -1
  DO i = 1,No_OF_DIFF_IONS
    IF( PP_Params(i)%AtomNum == AtomNum ) THEN
      res = i
      EXIT
    ENDIF
  ENDDO
END FUNCTION

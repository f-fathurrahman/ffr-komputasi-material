!-------------------------------------------------
CHARACTER(LEN=8) FUNCTION get_pp_filename(AtomNum)
!-------------------------------------------------
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

  Get_PP_Filename = "PP//PP"//StrAtomNum

END FUNCTION
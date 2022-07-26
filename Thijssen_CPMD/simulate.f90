!*************************************************************
!****                  Simulate                           ****
!*************************************************************
SUBROUTINE simulate()
  USE globals, only: OnlyStatic
  IMPLICIT NONE
  CALL calc_orbitals()
  IF( .NOT. OnlyStatic ) CALL run_carpar()
END SUBROUTINE


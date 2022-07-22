!*************************************************************
!****                  Simulate                           ****
!*************************************************************
SUBROUTINE Simulate()
  use globals, only: OnlyStatic
  IMPLICIT NONE
  CALL calc_orbitals()
  IF( .NOT. OnlyStatic ) CALL run_carpar()
END SUBROUTINE Simulate
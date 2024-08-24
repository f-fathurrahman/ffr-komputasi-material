PROGRAM pw_dft
  USE PW
  
  IMPLICIT NONE
  
  ! input and initialize
  CALL initialize()
  
  ! precalculate matrix elements
  CALL phase_factor()
  CALL structure_factor()
  CALL form_factor()
  CALL nl_pp_form_factor()

  IF( BAND_STRUCTURE .eqv. .true. ) THEN 
    ! band structure calculation
    CALL band_dft()
  ELSE 
    ! scf calculation
    ! initialize the density and wave function by diagonalizing on a small basis
    dens_mix = 0.95
    CALL solve_small()
    dens_mix = 0.95
    ! scf 
    CALL scf_dft()
  ENDIF 
  
  CALL calculate_dos()   

  WRITE(*,*)
  WRITE(*,*) 'Program pw_dft ended normally'
  WRITE(*,*)

END PROGRAM 



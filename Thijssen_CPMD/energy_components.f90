module energy_components

implicit none

complex(8) :: E_kin, E_locPP, E_locPPsr, E_hartree, E_xc, E_nonlocPP
complex(8) :: E_ovrl, E_self, E_totSelf, E_core

contains

subroutine print_energy_components()

  implicit none
  complex(8) :: E_KS

  E_KS = E_kin + E_locPPsr + E_xc + E_nonlocPP + E_self + E_ovrl - E_totSelf

  print '(A23 F15.8)', 'kin:', DBLE(E_kin)
  print '(A23 F15.8)', 'pp_sr:', DBLE(E_locPPsr)
  print '(A23 F15.8)', 'loc pp:', DBLE(E_locPP)
  print '(A23 F15.8)', 'xc:', DBLE(E_xc)
  print '(A23 F15.8)', 'Hartree energy:', DBLE(E_hartree)
  print '(A23 F15.8)', 'Nonlocal Psp:', DBLE(E_nonlocPP)
  print '(A23 F15.8)', 'Local core energy:', DBLE(E_core)
  print '(A23 F15.8)', 'self-energy:', DBLE(E_self)
  print '(A23 F15.8)', 'ovrl:', DBLE(E_ovrl)  
  print '(A23 F15.8)', 'Total energy:', DBLE(E_KS) 
  print *
end

end
subroutine fortranmain(mode)
  integer, intent(in) :: mode

  write(*, '(a)') 

  select case(mode)
    case(1); call coeff
    case(2); call test_hartree
    case(3); call test_laplacian
    case(4); call test_exp
    case(5); call gs
    case(6); call td
    case(7); call strength_function
    case(8); call excitations
  end select

end subroutine fortranmain

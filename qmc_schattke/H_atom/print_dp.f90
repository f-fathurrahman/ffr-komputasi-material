program main
  integer, parameter :: dp=selected_real_kind(2*precision(1.0))
  write(*,*) 'dp = ', dp
end program
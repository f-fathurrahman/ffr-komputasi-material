module m_constants

implicit none
integer, parameter, public :: dp=selected_real_kind(2*precision(1.0))
real(8), parameter, public :: PI=3.1415926535897932d0

end module
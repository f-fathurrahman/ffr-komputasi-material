!------------------
program test_random
!------------------
  use m_random
  implicit none
  integer :: i
  real(8) :: r
  !
  call initran()
  do i = 1,5
    call genran(r)
    write(*,'(1x,I4,F18.10)') i, r
  enddo

  write(*,*)
  write(*,*) 'Using gfortran random_number'
  !call set_seed()
  !call test_seed()
  do i = 1,5
    call random_number(r)
    write(*,'(1x,I4,F18.10)') i, r
  enddo

end program


subroutine test_seed()
  implicit none
  integer, allocatable :: seed(:)
  integer :: n
  !
  call random_seed(size = n)
  allocate(seed(n))
  call random_seed(get=seed)
  write(*,*) seed
end subroutine

!--------------------
subroutine set_seed()
!--------------------
  implicit none
  integer, allocatable :: seed(:)
  integer :: n
  !
  call random_seed(size = n)
  allocate(seed(n))
  seed(:) = 123456
  call random_seed(put=seed)
end
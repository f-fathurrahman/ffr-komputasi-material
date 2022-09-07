!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PROGRAM COEFF
! =============
!
! It outputs the coefficients of the nine-point formula for the second
! derivative in a regular one-dimensional grid (a Laplacian in two or three
! dimensions may be easily worked out by summing the terms).
!
! By changing the variables "m" (derivative order) and "n" (order of the
! approximation, e.g 4 for a 2*4+1=9 points formula), you may obtain the
! coefficients for other cases.
!*/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine coeff

  integer :: m, n, ierr
  real(8), allocatable :: c(:), x(:)
  real(8) :: delta

  delta = 1.0

  m = 2
  n = 4

  allocate(c(2*n+1), x(2*n))

  x(1) = delta
  x(2) = 2*delta
  x(3) = 3*delta
  x(4) = 4*delta
  x(5) = -delta
  x(6) = -2*delta
  x(7) = -3*delta
  x(8) = -4*delta

  call coefficients(m, n, x, c, ierr)

  ! Note that the indexes are changed... just notation.
  write(*,'(a,e16.7)')  'c(0)     = ', c(1)
  write(*, *)
  write(*,'(a,8e16.7)') 'c( 1: n) = ', c(2:n+1)
  write(*, *)
  write(*,'(a,8e16.7)') 'c(-1:-n) = ', c(n+2:2*n+1)

end subroutine coeff




!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! FUNCTION COEFF
! ==============
!
! INPUT:
!   m [integer] : the order of the derivatives representation.
!   n [integer] : The number of points given is 2*n
!   x [real(8), dimension(2*n)] : positions of the points. The "problem" point position is not given,
!     and assumed to be zero. 
! ---------
! OUTPUT:
!   c [real(8), dimension(2*n+1)] : the coefficients of the points. The first one corresponds to the
!     the coefficient at the problem points (which is always minus the sum of all the others), whereas 
!     the rest are ordered in the same manner that were given in array x.
!   coeff [integer] : error code. It is the error code of the LAPACK subroutine dgels
!
! Calculates the coefficients for the representation of the m-th order
! derivative of a function at a given point, given that we will have access 
! to the values of this function at 2*n points around it (besides the value
! of this function at the problem point). Typically this means n points to the
! left, and n points to the right, but that is not mandatory.
!
! NOTES:
! ------
! It requires BLAS dgesv subroutine.
!*/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine coefficients(m, n, x, c, ierr)
  integer, intent(in) :: m, n
  real(8), intent(in) :: x(2*n)
  real(8), intent(out) :: c(2*n+1)
  integer :: ierr
  integer :: i, j, k, lwork, info
  real(8), allocatable :: a(:, :), e(:), work(:)
  
  allocate(a(2*n, 2*n), e(2*n))

  do i = 1, 2*n
     do j = 1, 2*n
        a(i, j) = x(j)**i
     enddo
  enddo

  k = 1
  e = 0.0
  do i = 1, 2*n
     k = k*i
     if(m==i) then
       e(i) = k
       exit
     endif
  enddo

  lwork = -1
  allocate(work(1))
  call dgels('n', 2*n, 2*n, 1, a, 2*n, e, 2*n, work, lwork, info)
  lwork = int(work(1))
  deallocate(work); allocate(work(lwork))
  call dgels('n', 2*n, 2*n, 1, a, 2*n, e, 2*n, work, lwork, info)

  c(1) = - sum(e(1:2*n))
  do j = 1, 2*n
     c(j+1) = e(j)
  enddo

  deallocate(work, a, e)
  ierr = int(info)
end subroutine 







module expodata

  complex*16, allocatable :: a(:, :)
  integer :: n

end module expodata


subroutine matvec(x, y)
  use expodata
  implicit none
  !integer :: n
  complex*16 :: x(n), y(n)
  !complex*16 :: a(n, n)
  integer :: i, j
  do i = 1, n
     y(i) = (0.0_8, 0.0_8)
     do j = 1, n
        y(i) = y(i) + a(i, j)*x(j)
     enddo
  enddo
end subroutine matvec



subroutine test_exp
  use expodata
  implicit none

  complex*16, allocatable :: v(:), w(:), wsp(:)
  integer :: m, lwsp, liwsp, itrace, i, j, iflag, ns, iexph
  integer, allocatable :: iwsp(:), ipiv(:)
  double precision :: t, tol, anorm
  external :: matvec
  complex*16, parameter :: zi = (0.0_8, 1.0_8)

  n = 3 ! Dimension of the matrix
  m = min(n-1, 15) ! maximum size of the Krylov basis.
  t = 1.00_8
  tol = 1.0e-6_8 ! Error tolerance.
  lwsp = n*(m+2)+(m+2)**2+5*(m+2)**2+6+1
!!$  lwsp = 4*n*n+6+1.
  liwsp = m + 5
  itrace = 0
  iflag = 0

  allocate(a(n, n), v(n), w(n), wsp(lwsp), iwsp(liwsp), ipiv(n))

  a(1, 1) = 3.0_8
  a(1, 2) = 5.0_8
  a(1, 3) = 2.0_8
  a(2, 1) =-1.0_8
  a(2, 2) =-2.0_8
  a(2, 3) = 0.0_8
  a(3, 1) = 0.0_8
  a(3, 2) =-5.0_8
  a(3, 3) = 2.0_8

  anorm = sqrt(10.0_8)

  v(1) = (1.0_8, 0.0_8)
  v(2) = (0.0_8, 0.0_8)
  v(3) = (0.0_8, 0.0_8)

!!$  call ZGPADM(6,n,t,a,n,wsp,lwsp,ipiv,iexph,ns,iflag)

!!$  do i = 1, n**2
!!$     write(*, *) wsp(iexph+i-1)
!!$     write(*, *)
!!$  enddo

  write(*, *) 'step1'
  call ZGEXPV( n, m, t, v, w, tol, anorm, wsp, lwsp, iwsp, liwsp, matvec, itrace, iflag)
  write(*, *) 'step2'


  write(*, *) 'iflag = ', iflag
  do i = 1, n
     write(*, '(4es18.8)') w(i)
  enddo

  deallocate(a, v, w)



end subroutine test_exp

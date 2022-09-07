!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! MODULE MESH
! ===========
!
! This module should defines the mesh, that is, the simulation arena in a 
! real-space representation based scheme. It should contain:
! (i) the necessary parameters that define the grid (i.e. number of points, 
!     grid spacing...),
! (ii) functions to relate the "indexes" of the arrays that define each
!     function to real coordinates,
! (iii) the definition of the Hilbert space, which amounts to defining a
!     dot product between functions defined on the grid, and
! (iv) the needed differential operators, which in this case is only the
!     laplacian.
!*/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
MODULE mesh
  IMPLICIT NONE

!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! First of all, we need to define the dimensions of the simulation cell, and
! the grid spacing. The next variables are suitable for the examples that will
! be done later.
!
! We will use a regular rectangular mesh, of equal sides (that is, a square).
! This can be easily extended to more complicated geometries, but let us
! keep it simple for the moment.
!*/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  REAL(8), PARAMETER :: mesh_length = 50.d0 ! the length L of the cube
  INTEGER, PARAMETER :: n = 81       ! total number of points in each direction
  REAL(8), PARAMETER :: delta = mesh_length/(N-1)  ! spacing between the points

CONTAINS 

!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Next functions assign the real "x" and "y" coordinate to a set of integer
! indexes. The following definitions maps the indexes onto the square
! [-L/2, L/2]. It is better if we assume that N is an odd number, and in this
! way (0, 0) belongs to the mesh, and we have an equal number points in each
! direction.
!
! These definitions are once again just simple-minded ones; one could define
! more elaborate, non-regular meshes. Note also that we make "x" depend on
! both ix and iy, which is not necessary for this simple example.
!*/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
REAL(8) FUNCTION x(ix, iy)
   INTEGER, INTENT(in) :: ix, iy
   !x = (ix-1)*delta - (N/2)*delta
   x = (ix-1)*delta - 0.5d0*N*delta
END FUNCTION 

REAL(8) FUNCTION y(ix, iy)
   INTEGER, INTENT(in) :: ix, iy
   y = (iy-1)*delta - 0.5d0*N*delta
END FUNCTION 


!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! To define the Hilbert space necessary for any Quantum-Mechanical code, we
! need to define a dot product. This defines the norm and the distance.
!
! Note that we have two dot products, one for real functions and one for
! complex functions. If you know some Fortran 90, you can bundle the two
! definitions together by means of an module interface; in this way you do not
! have to worry about how the functions are in the rest of the code.
!*/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
REAL(8) FUNCTION dotproduct(a, b) RESULT(r)
  REAL(8), INTENT(in) :: a(n, n), b(n, n)
  r = sum(a(:,:)*b(:, :))*delta**2
END FUNCTION 

COMPLEX(8) FUNCTION zdotproduct(a, b) RESULT(r)
  COMPLEX(8), INTENT(in) :: a(n, n), b(n, n)
  r = sum(conjg(a(:,:))*b(:, :))*delta**2
END FUNCTION 


!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE LAPLACIAN.
! =====================
!
! INPUT:  
!   f [real(8), dimension(n, n)] : the function whose Laplacian is calculated
! ---------
! OUTPUT: 
!   lapl [real(8), dimension(n, n)] : the Laplacian of f.
!
! The kinetic operator is a Laplacian in real space, so we need a procedure
! that calculates the Laplacian of a function.
!*/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE laplacian(f, lapl)

  IMPLICIT NONE

  REAL(8), INTENT(in)  :: f(N, N)    ! function whose laplacian is to be taken
  REAL(8), INTENT(out) :: lapl(N, N) ! the laplacian

  INTEGER :: ix, iy, k

  INTEGER, PARAMETER :: order = 4
  REAL(8), ALLOCATABLE :: c(:)

  ALLOCATE(c(-order:order))
  c(-order:order) = (/ &  ! The coefficients of the laplacian...
      -1.785714d-3, 2.539683d-2, -0.2d0, 1.6d0,      &
      -2.847222d0,                                   &
      1.6d0, -0.2d0, 2.539683d-2, -1.785714d-3 /)

  DO ix = 1, N
    DO iy = 1, N
      lapl(ix, iy) = 0.d0
      DO k = -order, order
        IF(iy+k>=1 .and. iy+k<=N) lapl(ix, iy) = lapl(ix, iy) + c(k)*f(ix, iy+k)
        IF(ix+k>=1 .and. ix+k<=N) lapl(ix, iy) = lapl(ix, iy) + c(k)*f(ix+k, iy)
      ENDDO 
    ENDDO 
  ENDDO 

  lapl(:,:) = lapl(:,:)/delta**2

END SUBROUTINE laplacian


!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE ZLAPLACIAN.
! =====================
!
! INPUT:  
!   f [complex(8), dimension(n, n)] : the function whose Laplacian is calculated
! ---------
! OUTPUT: 
!   lapl [complex(8), dimension(n, n)] : the Laplacian of f.
!
! This is exactly the same that laplacian, but for complex functions. The
! missing code is identical.
!*/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
SUBROUTINE zlaplacian(f, lapl)

  IMPLICIT NONE

  COMPLEX(8), INTENT(in)  :: f(N, N)    ! function whose laplacian is to be taken
  COMPLEX(8), INTENT(out) :: lapl(N, N) ! the laplacian

  INTEGER :: ix, iy, k

  INTEGER, PARAMETER :: order = 4
  REAL(8), ALLOCATABLE :: c(:)

  ALLOCATE(c(-order:order))
  c(-order:order) = (/ &  ! The coefficients of the laplacian...
      -1.785714d-3, 2.539683d-2, -0.2d0, 1.6d0,      &
      -2.847222d0,                                   &
      1.6d0, -0.2d0, 2.539683d-2, -1.785714d-3 /)

  DO ix = 1, N
    DO iy = 1, N
      lapl(ix, iy) = (0.d0, 0.d0)
      DO k = -order, order
        IF(iy+k>=1 .and. iy+k<=N) lapl(ix, iy) = lapl(ix, iy) + c(k)*f(ix, iy+k)
        IF(ix+k>=1 .and. ix+k<=N) lapl(ix, iy) = lapl(ix, iy) + c(k)*f(ix+k, iy)
      ENDDO 
    ENDDO 
  ENDDO 

  lapl(:,:) = lapl(:,:)/delta**2

END SUBROUTINE 




end module mesh





!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! PROGRAM TEST_LAPLACIAN
! ======================
!
! This program is built with the purpose of testing the implementation of the
! Laplacian subroutines in the mesh module. It builds a Gaussian distribution
! (normalized so that it integrates to one), and calculates its Laplacian.
! Then it compares it with the exact analytical vales, and outputs the error,
! defined to be <(f_approx - f_exact)|(f_aaprox - f_exact)>
!
! You may try different function, or vary the "hardness" of the function by
! changing the alpha parameter. You may also check how the error depends
! on the order of the discretization of the Laplacian. 
!*/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine test_laplacian
  USE mesh
  IMPLICIT NONE 
  INTEGER :: ix, iy
  REAL(8), ALLOCATABLE :: rho(:,:), exl(:,:), apl(:,:)
  REAL(8) :: r2, alpha
  REAL(8), PARAMETER :: pi = 3.141592653589793d0

  ALLOCATE(rho(n, n), exl(n, n), apl(n, n))

  ! Define the problem density, which is just a Gaussian.
  alpha = 3.0d0
  DO ix = 1, n
     DO iy = 1, n
        r2 = x(ix,iy)**2 + y(ix,iy)**2
        rho(ix,iy) = exp(-r2/alpha**2)
     ENDDO 
  ENDDO 
  rho(:,:) = rho(:,:)/(alpha**2*pi)

  ! The exact value of the Laplacian of the problem density is put in exl variable
  DO ix = 1, n
  DO iy = 1, n
    r2 = x(ix, iy)**2 + y(ix, iy)**2
    exl(ix,iy) = (4.0d0/alpha**2)*(r2/alpha**2-1.0d0)*rho(ix, iy)
  ENDDO 
  ENDDO 

  ! Calculate the Laplacian of the problem density throught the laplacian
  ! subroutine, and put the result into apl variable.
  CALL laplacian(rho, apl)

  ! For visualization purposes:
  CALL output(rho,'rho.dat')
  CALL output(exl,'exact_laplacian.dat')
  CALL output(apl,'approximated_laplacian.dat')

  ! Outputs an estimation of the error of the laplacian subroutine.
  WRITE(*,'(A,ES20.8)') 'Error: ', dotproduct(apl-exl,apl-exl)

  DEALLOCATE(rho, exl, apl)
END subroutine


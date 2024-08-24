PROGRAM transfer_matrix
  IMPLICIT NONE 
  COMPLEX(8), PARAMETER :: zi=(0.d0,1.d0)
  INTEGER, PARAMETER :: N=5000, N_e=200
  REAL(8), PARAMETER :: xl=-15.d0, xu=15.d0, h2m=0.5d0
  REAL(8) :: x,xp,e,h
  INTEGER :: j,i
  COMPLEX(8) :: a,b,ap,bp,k,kp,k0,ec,v,vp

  ! function
  COMPLEX(8) :: potential

  h = (xu - xl)/dble(N+1)

  DO i=1,N_e
    V = potential(xu)
    e = 0.d0 + 0.5d0*i
    ec = (E-V)/h2m
    k0 = sqrt(ec)
    ap = (1.d0,0.d0)
    bp = (0.d0,0.d0)
    DO j = N,0,-1
      x = xl+j*h
      V = potential(x)
      xp = x+h
      Vp = potential(xp)
      
      ec = (E-V)/h2m
      k = sqrt(ec)
      
      ec = (E-Vp)/h2m
      kp = sqrt(ec)
      
      a = 0.5d0*(ap*(1.d0 + kp/k)*exp(zi*kp*x) + &
          bp*(1.d0-kp/k)*exp(-zi*kp*x))*exp(-zi*k*x)
      
      b = 0.5d0*(ap*(1.d0 - kp/k)*exp(zi*kp*x) + &
          bp*(1.d0 + kp/k)*exp(-zi*kp*x))*exp(zi*k*x)
      ap = a
      bp = b
    ENDDO 
    WRITE(1,*) e, abs(b/a)**2, abs(abs(1.d0/a)**2/k*k0)
  ENDDO

END PROGRAM 


FUNCTION potential(x) result(v)
  IMPLICIT NONE 
  COMPLEX(8) :: v
  REAL(8) :: x
  !
  REAL(8), PARAMETER :: V0 = 2.5d0
  REAL(8), PARAMETER :: d = 1.d0
  REAL(8), PARAMETER :: x0 = 3.d0
  REAL(8), PARAMETER :: mu = 0.2d0

  v = -V0*sinh((x-x0)/d)**2/cosh((x-x0)/d-mu)**2

end function 

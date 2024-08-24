PROGRAM fd_td
  USE linalg
  IMPLICIT NONE
  INTEGER, PARAMETER :: n=200, nt=500
  REAL(8), PARAMETER :: h2m=0.5d0, a=-25.d0, b=25.d0
  REAL(8) :: h(n,n)
  REAL(8) :: u(n,n)
  REAL(8) :: dx, x, dt, norm, p0, a0, t
  COMPLEX(8) :: zo(n,n), am(n,n), bm(n,n), cm(n,n), c(n)
  COMPLEX(8) :: f, csu, hc(n,n)
  COMPLEX(8), PARAMETER :: zi=(0.d0,1.d0)
  INTEGER :: i, k
  REAL(8), PARAMETER :: pi=3.1415926535897932384d0
  COMPLEX(8) :: Gauss_wave_packet_1d

  dx = (b-a)/(n+1)
  u(:,:) = 0.d0
  DO i = 1,n
    u(i,i) = 1.d0
    x = a + i*dx
  ENDDO

  t = h2m/dx**2

  DO i=1,N
    h(i,i) = 205.d0/72.d0*t
    IF( i > 1 ) THEN 
      h(i,i-1) = -8.d0/5.d0*t
      h(i-1,i) = -8.d0/5.d0*t
    ENDIF 
    IF( i > 2 ) THEN 
      h(i,i-2) = 1.d0/5.d0*t
      h(i-2,i) = 1.d0/5.d0*t
    ENDIF 
    IF( i > 3 ) THEN 
      h(i,i-3) = -8.d0/315.d0*t
      h(i-3,i) = -8.d0/315.d0*t
    ENDIF 
    IF( i > 4 ) THEN 
      h(i,i-4) = 1.d0/560.d0*t
      h(i-4,i) = 1.d0/560.d0*t
    ENDIF 
  ENDDO 

  hc = h
   
  dt = 0.001d0
  am = u - 0.5d0*zi*hc*dt
  zo = u + 0.5d0*zi*hc*dt

  CALL inv(zo, n, bm)
  cm = matmul(am,bm)

  ! initial wave function
  a0 = 1.d0
  p0 = 1.d0
  t  = 0.d0
  DO i = 1,N
    x = a + i*dx
    f=Gauss_wave_packet_1d(x,t,a0,p0)
    c(i)=f
  ENDDO 

  norm = 0.d0
  do k=1,n
    x = a+k*dx
    csu = c(k)
    norm = norm + real(conjg(csu)*csu*dx)
    WRITE(1,*) x, real(csu), imag(csu)
    WRITE(2,*) x, real(conjg(csu)*csu)
  ENDDO 
  WRITE(6,*) norm

  DO i=1,nt
    c = matmul(cm,c)
  ENDDO 
 
  ! Numerical
  t = nt*dt
  norm = 0.d0
  DO k=1,n
    x = a + k*dx
    csu = c(k)
    norm = norm + real(conjg(csu)*csu*dx)
        
    WRITE(11,*) x, real(csu), imag(csu)
    WRITE(12,*) x, real(conjg(csu)*csu)
    
  ENDDO 
  WRITE(6,*) norm


  !ANALYTICAL
  t = nt*dt
  DO k = 1,n
    x = a + k*dx
    f = Gauss_wave_packet_1d(x,t,a0,p0)
    WRITE(21,*) x, real(f), imag(f)
    WRITE(22,*) x, real(Conjg(f)*f)
  ENDDO 
  WRITE(6,*) norm

END PROGRAM 


FUNCTION Gauss_wave_packet_1d(x,t,a0,p0)
!This function computes the value of a 1D Gausian wave packet propagating in free space.
!Parameters:
!  x [in Angstrom] - spatial coordinates of the point where the packet needs to be computed
!  t [in femtosec]    - time
!  a [in Angstrom]    - spatial width of the packet at t=0
!  p0x [in eV*fs/Angstrom] - the average momentum of the electron in the wave packet
!
!It is assumed that the following global constants have been defined
!  hbar=0.658211899d0 [eV*fs]; H2M=3.80998174348d0 [eV*Angstom^2]; ElectronMass=hbar*hbar/(2*H2M) [eV*Angstom^2]
!
!It is also possible to use atomic units for input/output of this function. In this case one needs
!to pass all input values in atomic units while the global constants must be set as follows:
!  hbar=1.0;  H2M=0.5d0;  ElectronMass=1.0d0;
  IMPLICIT NONE 
  REAL(8) :: x,t,a0,p0
  COMPLEX(8) :: Gauss_wave_packet_1d,f,dc
  COMPLEX(8), PARAMETER :: zi=(0.d0,1.d0)
  REAL(8), PARAMETER :: pi=3.1415926535897932384d0

  dc = 1.d0+2.d0*zi*t/a0**2
  f = (2.d0/(pi*a0**2*dc**2))**0.25d0*exp(-(x-p0*t)**2/(a0**2*dc))*exp(zi*p0*(x-0.5d0*p0**2*t))
  Gauss_wave_packet_1d = f

END FUNCTION 



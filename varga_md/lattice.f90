PROGRAM lattice
  implicit none
  call initialize_fcc
END PROGRAM lattice

!------------------------
SUBROUTINE initialize_fcc
!------------------------
  implicit none
  integer            :: ncx,ncy,ncz,N_particles,i,j,ix,iy,iz,il,ic
  real(8)             :: a0, b0, c0, fac, t, v_cm(3), T_0, K_E, L(3)
  real(8),allocatable :: r(:,:), v(:,:)
  real(8),external    :: random_gaussian

  ! Read in the cell size, number of cells, and desired temperature
  open(1, file='fcc.inp')
  read(1,*) a0, b0, c0     ! x,y,z size of cell
  read(1,*) ncx, ncy, ncz  ! number of cells in x,y,z directions
  read(1,*) T_0          ! desired temperature  
  close(1)

  ! Allocate memory
  N_particles = 8*ncx*ncy*ncz ! 8 atoms per cell
  allocate(r(3,N_particles), v(3,N_particles))

  ! Define the coordinates in the first cell (FCC lattice)
  r(1,1) = 0.0      ;  r(2,1) = 0.0      ;  r(3,1) = 0.0
  r(1,2) = 0.5d0*a0 ;  r(2,2) = 0.5d0*b0 ;  r(3,2) = 0.0
  r(1,3) = 0.0      ;  r(2,3) = 0.5d0*b0 ;  r(3,3) = 0.5d0*c0
  r(1,4) = 0.5d0*a0 ;  r(2,4) = 0.0      ;  r(3,4) = 0.5d0*c0
  r(1,5) = 0.25d0*a0;  r(2,5) = 0.25d0*b0;  r(3,5) = 0.25d0*c0
  r(1,6) = 0.75d0*a0;  r(2,6) = 0.75d0*b0;  r(3,6) = 0.25d0*c0
  r(1,7) = 0.25d0*a0;  r(2,7) = 0.75d0*b0;  r(3,7) = 0.75d0*c0
  r(1,8) = 0.75d0*a0;  r(2,8) = 0.25d0*b0;  r(3,8) = 0.75d0*c0

  ! Calculate the total size occupied by the cells
  L(1) = ncx*a0
  L(2) = ncy*b0
  L(3) = ncz*c0

  ! Copy the coordinates into the other cells
  il = 0
  do iz = 0,ncz-1
    do iy = 0,ncy-1
      do ix = 0,ncx-1
        do ic = 1,8
          r(1,ic+il) = r(1,ic) + a0*ix
          r(2,ic+il) = r(2,ic) + b0*iy
          r(3,ic+il) = r(3,ic) + c0*iz
        enddo
        il = il + 8 ! increment by 8 (number of atoms for each cells)
      enddo
    enddo
  enddo

  ! Randomly set initial velocities
  do i = 1,N_particles
    v(1,i) = random_gaussian(1.d0)
    v(2,i) = random_gaussian(1.d0)
    v(3,i) = random_gaussian(1.d0)
  enddo

  ! Calculate the center of mass velocity
  v_cm = 0.d0
  do i = 1,N_particles
    v_cm(:) = v_cm(:) + v(:,i)
  enddo

  ! Make the center of mass velocity zero and calculate the kinetic energy
  K_E = 0.d0
  do i=1,N_particles
    v(:,i) = v(:,i) - v_cm(:)/N_particles
    K_E = K_E + 0.5d0*(v(1,i)**2 + v(2,i)**2 + v(3,i)**2)
  enddo

  ! Scale velocities to obtain desired temperature
  if(T_0 > 0.d0) then
    T = (2.d0/3.d0)*(K_E/N_particles)
    fac = sqrt(T_0/T)
    K_E = 0.d0
    v = v*fac
    do i = 1,N_particles
      K_E = K_E + 0.5d0*(v(1,i)**2+v(2,i)**2+v(3,i)**2)
    enddo
  endif

  ! Output to two files: first a file with positions and velocities, and second
  ! a file in XYZ format with just the positions
  open(1, file='coord.inp')
  open(2, file='coord.xyz')
  write(2,*) N_Particles; write(2,*)
  write(1,'(3ES16.8)') L(1), L(2), L(3)
  write(1,*) N_Particles
  do i=1,N_Particles
    write(1,'(6ES16.8)')(r(j,i),j=1,3),(v(j,i),j=1,3)
    write(2,'(a,3ES16.8)')"Si",(r(j,i),j=1,3)
  enddo
  close(1)
  close(2)

  ! Deallocate memory
  deallocate(r,v)
END SUBROUTINE initialize_fcc


!------------------------------
FUNCTION random_gaussian(sigma)
!------------------------------
  ! Calculate a random number from a Gaussian distribution
  implicit none
  real(8)           :: random_gaussian,sigma,x,y
  real(8),parameter :: pi2=2.d0*3.14159265358979d0
  call random_number(x)
  call random_number(y)
  random_gaussian = sqrt(sigma)*sqrt(-2.*log(x))*cos(pi2*y)
END FUNCTION random_gaussian

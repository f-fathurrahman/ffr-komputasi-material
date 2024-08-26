!------------------
MODULE m_hubbard_2d
!------------------

IMPLICIT NONE

! Lattice size
INTEGER :: N_x,N_y

! Lattice dimension
INTEGER :: Ndim

! Boundary condition
LOGICAL :: PBC

! Hopping matrix
INTEGER, ALLOCATABLE :: hopping_matrix(:,:)

! Hubbard U (must be positive)
real(8) :: U

! number of slices
integer :: N_S

! chemical potential
real(8) :: Mu

! delta_t for Suzuki-Trotter
real(8) :: delta_t

! Automatic arrays
real(8) :: exp_up(-1:1), exp_down(-1:1)

! sigma
integer, allocatable :: Sigma(:,:)

! random seed
integer  :: idum

real(8), allocatable :: e_up(:,:), e_down(:,:), t_m(:,:)

!
integer :: N_spin_flip,N_steps,N_thermolization

CONTAINS

!----------------------
subroutine initialize()
!----------------------
  implicit none
  integer                  :: i,j,m,mx
  real(8)                   :: lambda,v_p,v_m,x
  real(8),external          :: ran2
  real(8),allocatable       :: am(:,:),eval(:),evec(:,:)

  open(1,file='hub2d.inp')
  read(1,*) N_x, N_y
  read(1,*) U
  read(1,*) mu
  read(1,*) N_S
  read(1,*) N_spin_flip
  read(1,*) N_steps
  read(1,*) N_thermolization
  read(1,*) PBC
  close(1)
    
  Ndim = N_x*N_y
  allocate(Hopping_Matrix(Ndim,Ndim), Sigma(Ndim,N_S))
  allocate(e_up(Ndim,N_S), e_down(Ndim,N_S), T_m(Ndim,Ndim))

  idum = 1231
  ! initialize the Hopping Matrix
  Hopping_matrix=0
  do i = 1,N_x-1
    do j = 1,N_y
      m = 1 + (i-1) + N_x*(j-1)
      Hopping_matrix(m,m+1) = 1
      Hopping_matrix(m+1,m) = 1
    enddo
  enddo

  do i=1,N_y-1
    do j=1,N_x
      m=j+N_x*(i-1)
      mx=m+N_x   
      Hopping_matrix(m,mx)=1
      Hopping_matrix(mx,m)=1
    enddo
  enddo

  ! boundary-condition: periodic hopping or not
  if( PBC ) then
    do i=1,N_y
      m = N_x + N_x*(i-1)
      Hopping_matrix(m,m-(N_x-1)) = 1
      Hopping_matrix(m-(N_x-1),m) = 1
    end do
  endif
  
  do i=2,N_y-1,2
    do j=1,N_x
      m=j+N_x*(i-1)
      Hopping_matrix(m,m+N_x)=1
      Hopping_matrix(m+N_x,m)=1
    end do                  
  enddo

  if (PBC .and. (N_y > 1)) then
    do i=1,N_x
      m = i + N_x*(N_y-1)
      Hopping_matrix(m,m-(Ndim-N_x)) = 1
      Hopping_matrix(m-(Ndim-N_x),m) = 1
    end do
  endif

  do j=1,Ndim
    do i=1,Ndim
      write(16,*)j,i,Hopping_matrix(j,i)
    end do
  end do

  lambda = 2.d0*atanh(sqrt(tanh(U*delta_t/4.d0)))
  v_p = Exp(+lambda-delta_t*Mu)
  v_m = Exp(-lambda-delta_t*Mu)

  exp_up(-1)=v_m
  exp_up(1)=v_p

  exp_down(-1)=v_p
  exp_down(1)=v_m
  ! initialize the HS fields
  do i=1,N_S
    do j=1,Ndim
      x=ran2(idum)
      if(x < 0.5d0) then
        sigma(j,i) = -1
      else
        sigma(j,i) = 1
      endif
    enddo
  enddo
  call calculate_fields()
  
  allocate( am(Ndim,Ndim), evec(Ndim,Ndim), eval(Ndim))
  
  am = Hopping_matrix*delta_t

  call diag_real(am,Ndim,eval,evec)
  
  T_M(:,:) = 0.d0
  do i=1,Ndim
    do j=1,Ndim
      do m=1,Ndim
        T_m(j,i) = T_m(j,i) + evec(m,i)*exp(eval(m))*evec(m,j)
      enddo
    enddo
  enddo

  write(*,*) 'Finished initialize'
 
end subroutine


!---------------------
subroutine spin_flip()
!---------------------
  implicit none
  integer       :: i,j,k,ii
  real(8)        :: x
  real(8),external          :: ran2

  ii=N_Spin_flip*Ndim*N_s
  do i=1,ii
    x=ran2(idum)
    k=idint(x*N_S)+1
    x=ran2(idum)
    j=idint(x*Ndim)+1
    if( sigma(j,k) == -1) then
      sigma(j,k)=1
    else
      sigma(j,k)=1
    endif
  enddo

end subroutine


!----------------------------
subroutine calculate_fields()
!----------------------------
  implicit none
  integer       :: i,j

  do i=1,N_S
    do j=1,Ndim
        e_up(j,i)=exp_up(sigma(j,i))
        e_down(j,i)=exp_down(sigma(j,i))
    end do
  end do
end subroutine

!------------------------------
subroutine Greens_function(A,G)
!------------------------------
  implicit none
  integer :: i,j,k
  real(8) :: A(Ndim,N_S),B(Ndim,Ndim),C(Ndim,Ndim),G(Ndim,Ndim)

  do i=1,N_S
    B=0.d0
    if(i==1) then
      do j=1,Ndim
        B(j,j)=A(j,i)
      end do
    else
      do j=1,Ndim
        do k=1,Ndim
          B(k,j)=C(k,j)*A(j,i)
        end do
      end do
    endif
    C=matmul(B,T_m)
  end do

  do i=1,Ndim
    C(i,i)=1.d0+C(i,i)
  end do
  
  call inv(C,Ndim,G)

end subroutine


!-----------------------
subroutine Monte_Carlo()
!-----------------------
  implicit none
  integer                :: i,j
  real(8), allocatable :: G_up(:,:),G_down(:,:),G_up_new(:,:),G_down_new(:,:),A(:,:),A_new(:,:),B(:,:)
  real(8) :: d,d_new,x,w
  real(8), external :: ran2

  allocate(G_down(Ndim,Ndim),G_up(Ndim,Ndim),G_down_new(Ndim,Ndim),G_up_new(Ndim,Ndim), &
&   A(Ndim,Ndim),A_new(Ndim,Ndim),B(Ndim,N_S))

  B=E_up
  call Greens_function(A,G_up)
  B=E_down
  call Greens_function(A,G_down)

  do i=1,N_steps
    if( mod(i, 100) == 0 ) THEN
      write(*,*) 'Step: ', i
    endif
    do j=1,Ndim
      call spin_flip()
      call calculate_fields()
      B = E_up
      call Greens_function(A,G_up_new)
      B = E_down
      call Greens_function(A,G_down_new)
      A = matmul(G_up,G_down)
      call det(A,Ndim,d)
      A_new = matmul(G_up_new,G_down_new)
      call det(A,Ndim,d_new)
      w = d/d_new
      ! heat bath
      w = w/(1.d0+w)
      x = ran2(idum)
      if(abs(w) > x) then
        G_up=G_up_new
        G_down=G_down_new
      endif        
    enddo 
    ! measuerement
    if(i > N_thermolization) then
      ! do something?
    endif
  enddo
end subroutine


subroutine diag_real(A,n,e,v)
  !Arguments(input):
  real(8)   A(n,n)          !Symmetric matrix
  integer   n               !Eigenvalue problem size
  !Arguments(output):
  real(8)   e(n)   !Eigenvalues
  real(8)   v(n,n) !Corresponding eigenvectors

  !local variables
  integer,allocatable,dimension(:)   :: IWORK
  character(1)                       :: JOBZ, RANGE, UPLO
  real(8)                            :: VL, VU, ABSTOL
  real(8),allocatable,dimension(:)   :: WORK
  integer                            :: LDA, IL, IU, M, LDZ, NB, LWORK,INFO
  integer,allocatable,dimension(:)   :: IFAIL
  integer                            :: ILAENV

  JOBZ='V'       !'V' means compute both eigenvalues and eigenvectors
  RANGE='I'      !'I' means only selected eigenvalues/eigenvectors will be found
  UPLO='U'       !'U' means the upper triangle of A is provided.
  LDA=n          !Leading dimension of A is N
  VL=0.0D0       !VL can be set to anything because it is not used when
  RANGE='I'
  VU=0.0D0       !VU can be set to anything because it is not used when
  RANGE='I'
  IL=1           !The smallest eigenvalue number
  IU=n        !The largest eigenvalue number.
  ABSTOL=2*tiny(ABSTOL) !Set tolerance that yields highest possible accuracy in the
  !calculations; equivalent to  ABSTOL=2*DLAMCH('S')
  LDZ=n          !Leading dimension of v is N
  NB=ILAENV(1,'DSYTRD','VIU',n,n,n,n) !Determine the optimal block size
  LWORK=(NB+3)*N !Set the size of array WORK
  allocate(WORK(LWORK)) !Allocate array WORK, which is a real(8) work array used by DSYGVX
  allocate(IWORK(5*N))  !Allocate array IWORK, which is an integer work array used by DSYGVX
  allocate(IFAIL(N)) !Allocate IFAIL, an integer flag array used by DSYGVX

  call DSYEVX(JOBZ, RANGE, UPLO, n, A, LDA, VL, VU, IL, IU, &
&   ABSTOL, M, e, v, LDZ, WORK, LWORK, IWORK, &
&   IFAIL, INFO )

  if (INFO/=0) then
    write(6,*) 'Error in diag4: subroutine DSYEVX failed with INFO=',INFO
    stop
  endif

  deallocate(IFAIL)
  deallocate(IWORK)
  deallocate(WORK)
end subroutine diag_real

END MODULE


       SUBROUTINE lubksb_r(a,n,np,indx,b)
       INTEGER n,np,indx(n)
       real(8) a(np,np),b(n)
       INTEGER i,ii,j,ll
       real(8) sum
       ii=0
       do 12 i=1,n
         ll=indx(i)
         sum=b(ll)
         b(ll)=b(i)
         if (ii.ne.0)then
           do 11 j=ii,i-1
             sum=sum-a(i,j)*b(j)
11         continue
         else if (sum.ne.0.d0) then
           ii=i
         endif
         b(i)=sum
12     continue
       do 14 i=n,1,-1
         sum=b(i)
         do 13 j=i+1,n
           sum=sum-a(i,j)*b(j)
13       continue
         b(i)=sum/a(i,i)
14     continue
       return
       END SUBROUTINE lubksb_r
      

       SUBROUTINE ludcmp_r(a,n,np,indx,d)
       INTEGER n,np,indx(n),Ndim
       real(8) d,TINY
       real(8) a(np,np),sum,du
       PARAMETER (Ndim=5000,TINY=1.0d-20)
       INTEGER i,imax,j,k
       real(8) aamax,vv(ndim),dum
       d=1.d0
       do 12 i=1,n
         aamax=0.d0
         do 11 j=1,n
           if (abs(a(i,j)).gt.aamax) aamax=abs(a(i,j))
 11      continue
         if (aamax.eq.0.d0) stop 'singular matrix in ludcmp'
         vv(i)=1.d0/aamax
 12    continue
       do 19 j=1,n
         do 14 i=1,j-1
           sum=a(i,j)
           do 13 k=1,i-1
             sum=sum-a(i,k)*a(k,j)
 13        continue
           a(i,j)=sum
 14      continue
         aamax=0.d0
         do 16 i=j,n
           sum=a(i,j)
           do 15 k=1,j-1
             sum=sum-a(i,k)*a(k,j)
 15        continue
           a(i,j)=sum
           dum=vv(i)*abs(sum)
           if (dum.ge.aamax) then
             imax=i
             aamax=dum
           endif
 16      continue
         if (j.ne.imax)then
           do 17 k=1,n
             du=a(imax,k)
             a(imax,k)=a(j,k)
             a(j,k)=du
 17        continue
           d=-d
           vv(imax)=vv(j)
         endif
         indx(j)=imax
         if(a(j,j).eq.0.d0) a(j,j)=TINY
         if(j.ne.n)then
           du=1.d0/a(j,j)
           do 18 i=j+1,n
             a(i,j)=a(i,j)*du
 18        continue
         endif
 19    continue
       return
       END SUBROUTINE ludcmp_r

subroutine inv(a,n,ai)
implicit none
  integer      :: n,i
  real(8)       :: a(n,n),ai(n,n)
  integer      :: indx(n)
  real(8)       :: d

  ai=0.d0
  do i=1,n
    ai(i,i)=1.d0
  end do
  call ludcmp_r(a,n,n,indx,d)
  do i=1,n
    call lubksb_r(a,n,n,indx,ai(1,i))
  end do

end subroutine inv


subroutine det(a,n,d)
implicit none
  integer      :: n,i
  real(8)       :: a(n,n)
  integer      :: indx(n)
  real(8)       :: d

  call ludcmp_r(a,n,n,indx,d)
  do i=1,n
    d=d*a(i,i)
  end do

end subroutine det



      function ran2(idum)
      implicit real(8)(a-h,o-z)
      parameter (m=714025,ia=1366,ic=150889,rm=1.4005112d-6)
      common /com1 /iy,ir(97)
      data iff /0/
      if(idum.lt.0.or.iff.eq.0)then
        iff=1
        idum=mod(ic-idum,m)
        do 11 j=1,97
          idum=mod(ia*idum+ic,m)
          ir(j)=idum
11      continue
        idum=mod(ia*idum+ic,m)
        iy=idum
      endif
      j=1+(97*iy)/m
      if(j.gt.97.or.j.lt.1) stop 'Something wrong'
      iy=ir(j)
      ran2=iy*rm
      idum=mod(ia*idum+ic,m)
      ir(j)=idum
      return
      end


!-----------
program main
!-----------
  USE m_hubbard_2d
  call initialize()
  call Monte_Carlo()
end program



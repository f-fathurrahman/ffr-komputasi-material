MODULE linalg
  implicit none
  
  interface diag
    module procedure diag_gen_complex,diag_gen_complex_partial,diag_gen_real, &
      diag_gen_real_partial,diag_real,diag_real_partial,diag_complex
  end interface
    
  interface inv
    module procedure inv_real,inv_complex
  end interface

  interface det
    module procedure det_real,det_complex
  end interface  

CONTAINS

SUBROUTINE diag_complex(hamiltonian,N,eigenvalues,eigenvectors)
  implicit none
  integer    :: N,swapped,i,LWORK,retval
  real*8     :: rwork(2*N)
  complex*16 :: hamiltonian(N,N),eigenvectors(N,N),eigenvalues(N),temp_eigenvalue
  complex*16,allocatable :: temp(:,:),temp_eigenvector(:),work(:)

  ! This first call just finds the optimal size for "work"
  allocate(temp(1,N),temp_eigenvector(N))
  call zgeev('N','V',N,hamiltonian,N,eigenvalues,temp,1,eigenvectors,N,temp_eigenvector,-1,rwork,retval)
  LWORK=int(real(temp_eigenvector(1)))
  allocate(work(LWORK))

  ! Now call to actually diagonalize
  call zgeev('N','V',N,hamiltonian,N,eigenvalues,temp,1,eigenvectors,N,work,LWORK,rwork,retval)
  deallocate(temp,work)

  ! Sort the results by the real part of the eigenvalue
  swapped=1
  do while(swapped==1)
    swapped=0
    do i=1,N-1
      if(real(eigenvalues(i))>real(eigenvalues(i+1))) then
        temp_eigenvalue=eigenvalues(i)
        temp_eigenvector=eigenvectors(:,i)
        eigenvalues(i)=eigenvalues(i+1)
        eigenvectors(:,i)=eigenvectors(:,i+1)
        eigenvalues(i+1)=temp_eigenvalue
        eigenvectors(:,i+1)=temp_eigenvector
        swapped=1
      endif
    enddo
  enddo
  deallocate(temp_eigenvector)
END SUBROUTINE diag_complex

subroutine diag_gen_real(hm,om,n,e,v)
  real(8)   hm(n,n)         !Hamiltonian matrix
  real(8)   om(n,n)         !Overlap matrix
  integer   n               !Eigenvalue problem size
  real(8)   e(n)   !Eigenvalues
  real(8)   v(n,n) !Corresponding eigenvectors
  call diag_gen_real_partial(hm,om,n,n,e,v)
end subroutine diag_gen_real

subroutine diag_gen_complex(hm,om,n,e,v)
  complex(8)   hm(n,n)     !Hamiltonian matrix
  complex(8)   om(n,n)     !Overlap matrix
  integer   n               !Eigenvalue problem size
  real(8)   e(n)   !Eigenvalues
  complex(8)   v(n,n) !Corresponding eigenvectors
  call diag_gen_complex_partial(hm,om,n,n,e,v)
end subroutine diag_gen_complex

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

subroutine diag_real_partial(A,n,neig,e,v)
  !Arguments(input):
  real(8)   A(n,n)          !Symmetric matrix
  integer   n               !Eigenvalue problem size
  integer   neig            !Number of lowest eigenvalues/eigenvectors to compute
  !Arguments(output):
  real(8)   e(neig)   !Eigenvalues
  real(8)   v(n,neig) !Corresponding eigenvectors

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
  IU=neig        !The largest eigenvalue number.
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
end subroutine diag_real_partial
  
subroutine diag_gen_complex_partial(hm,om,n,neig,e,v)
!Arguments(input):
complex(8)   hm(n,n)     !Hamiltonian matrix
complex(8)   om(n,n)     !Overlap matrix
integer   n               !Eigenvalue problem size
integer   neig            !Number of lowest eigenvalues/eigenvectors to compute
!Arguments(output):
real(8)   e(neig)   !Eigenvalues
complex(8)   v(n,neig) !Corresponding eigenvectors

!local variables
integer,allocatable,dimension(:)   :: IWORK
real(8)                            :: ABSTOL
real(8),allocatable,dimension(:)  :: RWORK
complex(8),allocatable,dimension(:)   :: WORK
integer                            :: LDA, LDB, IU, M, LDZ, NB, LWORK,INFO
integer,allocatable,dimension(:)   :: IFAIL
integer                            :: ILAENV
character(1), parameter :: JOBZ='V',RANGE='I',UPLO='U'
real(8), parameter :: VL=0.0d0, VU=0.0d0
integer, parameter :: IL=1
integer, parameter :: ITYPE = 1




      !The problem type is A*x = lambda*B*x
      !'V' means compute both eigenvalues and eigenvectors
      !'I' means only selected eigenvalues/eigenvectors will be found
      !'U' means the upper triangles of A and B are provided.
LDA=N          !Leading dimension of A is N
LDB=N          !Leading dimension of B is N
       !VL can be set to anything because it is not used when RANGE='I'
       !VU can be set to anything because it is not used when RANGE='I'
       !The smallest eigenvalue number
IU=neig        !The largest eigenvalue number.
ABSTOL=2*tiny(ABSTOL) !Set tolerance that yields highest possible accuracy in the
                      !calculations; equivalent to  ABSTOL=2*DLAMCH('S')
LDZ=N          !Leading dimension of Z is N
NB=ILAENV(1,'DSYTRD','VIU',N,N,N,N) !Determine the optimal block size
LWORK=(NB+3)*N !Set the size of array WORK
allocate(WORK(LWORK)) !Allocate array WORK, which is a real(8) work array used by DSYGVX
allocate(IWORK(5*N))  !Allocate array IWORK, which is an integer work array used by DSYGVX
allocate(IFAIL(N)) !Allocate IFAIL, an integer flag array used by DSYGVX
allocate(RWORK(7*N))

call ZHEGVX( ITYPE, JOBZ, RANGE, UPLO, N, hm, LDA, om, LDB, VL, VU, IL, IU, &
     &  ABSTOL, M, e, v, LDZ, WORK, LWORK, RWORK, &
     &  IWORK, IFAIL, INFO )


if (INFO/=0) then
  write(6,*) 'Error in diag3: subroutine DSYGVX failed with INFO=',INFO
  stop
endif

deallocate(IFAIL)
deallocate(IWORK)
deallocate(WORK)
deallocate(RWORK)

end subroutine diag_gen_complex_partial

subroutine diag_gen_real_partial(hm,om,n,neig,e,v)
!Arguments(input):
real(8)   hm(n,n)         !Hamiltonian matrix
real(8)   om(n,n)         !Overlap matrix
integer   n               !Eigenvalue problem size
integer   neig            !Number of lowest eigenvalues/eigenvectors to compute
!Arguments(output):
real(8)   e(neig)   !Eigenvalues
real(8)   v(n,neig) !Corresponding eigenvectors

!local variables
integer,allocatable,dimension(:)   :: IWORK
character(1)                       :: JOBZ, RANGE, UPLO
real(8)                            :: VL, VU, ABSTOL
real(8),allocatable,dimension(:)   :: WORK
integer                            :: ITYPE, LDA, LDB, IL, IU, M, LDZ, NB, LWORK,INFO
integer,allocatable,dimension(:)   :: IFAIL
integer                            :: ILAENV
ITYPE=1        !The problem type is A*x = lambda*B*x
JOBZ='V'       !'V' means compute both eigenvalues and eigenvectors
RANGE='I'      !'I' means only selected eigenvalues/eigenvectors will be found
UPLO='U'       !'U' means the upper triangles of A and B are provided.
LDA=N          !Leading dimension of A is N
LDB=N          !Leading dimension of B is N
VL=0.0D0       !VL can be set to anything because it is not used when RANGE='I'
VU=0.0D0       !VU can be set to anything because it is not used when RANGE='I'
IL=1           !The smallest eigenvalue number
IU=neig        !The largest eigenvalue number.
ABSTOL=2*tiny(ABSTOL) !Set tolerance that yields highest possible accuracy in the
                      !calculations; equivalent to  ABSTOL=2*DLAMCH('S')
LDZ=N          !Leading dimension of Z is N
NB=ILAENV(1,'DSYTRD','VIU',N,N,N,N) !Determine the optimal block size
LWORK=(NB+3)*N !Set the size of array WORK
allocate(WORK(LWORK)) !Allocate array WORK, which is a real(8) work array used by DSYGVX
allocate(IWORK(5*N))  !Allocate array IWORK, which is an integer work array used by DSYGVX
allocate(IFAIL(N)) !Allocate IFAIL, an integer flag array used by DSYGVX

!Call DSYGVX
call DSYGVX( ITYPE, JOBZ, RANGE, UPLO, N, hm, LDA, om, LDB,    &
             VL, VU, IL, IU, ABSTOL, M, e, v, LDZ, WORK,     &
             LWORK, IWORK, IFAIL, INFO )
if (INFO/=0) then
  write(6,*) 'Error in diag3: subroutine DSYGVX failed with INFO=',INFO
  stop
endif

deallocate(IFAIL)
deallocate(IWORK)
deallocate(WORK)

end subroutine diag_gen_real_partial

SUBROUTINE inv_real(A,N,AI)
  integer             :: N,retval,lwork
  integer,allocatable :: ipiv(:)
  real*8              :: A(N,N),AI(N,N),temp(1)
  real*8,allocatable  :: work(:)

  AI=A
  allocate(ipiv(N))
  call dgetrf(N,N,AI,N,ipiv,retval)
  call dgetri(N,AI,N,ipiv,temp,-1,retval)
  lwork=int(temp(1))
  allocate(work(lwork))
  call dgetri(N,AI,N,ipiv,work,lwork,retval)
  deallocate(ipiv,work)
END SUBROUTINE inv_real

SUBROUTINE inv_complex(A,N,AI)
  integer                :: N,retval,lwork
  integer,allocatable    :: ipiv(:)
  complex*16             :: A(N,N),AI(N,N),temp(1)
  complex*16,allocatable :: work(:)

  AI=A
  allocate(ipiv(N))
  call zgetrf(N,N,AI,N,ipiv,retval)
  call zgetri(N,AI,N,ipiv,temp,-1,retval)
  lwork=int(temp(1))
  allocate(work(lwork))
  call zgetri(N,AI,N,ipiv,work,lwork,retval)
  deallocate(ipiv,work)
END SUBROUTINE inv_complex

SUBROUTINE det_real(A,N,retval)
  integer             :: i,N,info
  integer,allocatable :: ipiv(:)
  real*8              :: A(N,N),retval

  allocate(ipiv(N))
  call dgetrf(N,N,A,N,ipiv,info)
  retval=1.d0
  do i=1,N
    retval=retval*A(i,i)
    if(ipiv(i)/=i) retval=retval*(-1.d0)
  enddo
  deallocate(ipiv)
END SUBROUTINE det_real

SUBROUTINE det_complex(A,N,retval)
  integer             :: i,N,info
  integer,allocatable :: ipiv(:)
  complex*16          :: A(N,N),retval

  allocate(ipiv(N))
  call zgetrf(N,N,A,N,ipiv,info)
  retval=(1.d0,0.d0)
  do i=1,N
    retval=retval*A(i,i)
    if(ipiv(i)/=i) retval=retval*(-1.d0,0.d0)
  enddo
  deallocate(ipiv)
END SUBROUTINE det_complex

SUBROUTINE indexx(n,arr,indx)
      INTEGER :: n,indx(n),M,NSTACK
      REAL*8 :: arr(n)
      PARAMETER (M=7,NSTACK=50)
      INTEGER :: i,indxt,ir,itemp,j,jstack,k,l,istack(NSTACK)
      REAL*8 a
      do 11 j=1,n
        indx(j)=j
11    continue
      jstack=0
      l=1
      ir=n
1     if(ir-l.lt.M)then
        do 13 j=l+1,ir
          indxt=indx(j)
          a=arr(indxt)
          do 12 i=j-1,1,-1
            if(arr(indx(i)).le.a)goto 2
            indx(i+1)=indx(i)
12        continue
          i=0
2         indx(i+1)=indxt
13      continue
        if(jstack.eq.0)return
        ir=istack(jstack)
        l=istack(jstack-1)
        jstack=jstack-2
      else
        k=(l+ir)/2
        itemp=indx(k)
        indx(k)=indx(l+1)
        indx(l+1)=itemp
        if(arr(indx(l+1)).gt.arr(indx(ir)))then
          itemp=indx(l+1)
          indx(l+1)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l)).gt.arr(indx(ir)))then
          itemp=indx(l)
          indx(l)=indx(ir)
          indx(ir)=itemp
        endif
        if(arr(indx(l+1)).gt.arr(indx(l)))then
          itemp=indx(l+1)
          indx(l+1)=indx(l)
          indx(l)=itemp
        endif
        i=l+1
        j=ir
        indxt=indx(l)
        a=arr(indxt)
3       continue
          i=i+1
        if(arr(indx(i)).lt.a)goto 3
4       continue
          j=j-1
        if(arr(indx(j)).gt.a)goto 4
        if(j.lt.i)goto 5
        itemp=indx(i)
        indx(i)=indx(j)
        indx(j)=itemp
        goto 3
5       indx(l)=indx(j)
        indx(j)=indxt
        jstack=jstack+2
        if(jstack.gt.NSTACK) then
          write(6,*)'NSTACK too small in indexx'
          stop
        endif
        if(ir-i+1.ge.j-l)then
          istack(jstack)=ir
          istack(jstack-1)=i
          ir=j-1
        else
          istack(jstack)=j-1
          istack(jstack-1)=l
          l=i
        endif
      endif
      goto 1
      END subroutine indexx
END MODULE linalg


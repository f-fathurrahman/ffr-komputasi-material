!PROGRAM LAPACK_EIGENSOLVER_TEST
!  IMPLICIT NONE 

!  write(*,*);write(*,'(a)')"*** Testing real SUBROUTINEs ***"
!  call test_real

!  write(*,'(a)')"*** Testing complex SUBROUTINEs ***"
!  call test_complex
!END PROGRAM LAPACK_EIGENSOLVER_TEST

SUBROUTINE test_real
  IMPLICIT NONE 
  INTEGER  :: N,i
  real*8,allocatable :: A(:,:),B(:,:),eigenvalues(:),eigenvectors(:,:)

  N=2 ! A and B are NxN matrices

  allocate(A(N,N),B(N,N),eigenvalues(N),eigenvectors(N,N))

  ! Set the values for the A and B matrices
  A(1,1)=1.2d0; A(1,2)=3.4d0
  A(2,1)=5.6d0; A(2,2)=7.8d0

  B(1,1)=9.0d0; B(1,2)=1.2d0
  B(2,1)=3.4d0; B(2,2)=5.6d0

  ! Calculate the eigenvalues of A
  write(*,*)"Eigenvalues of real matrix A:"
  call diag_lapack_real(A,N,eigenvalues,eigenvectors)
  do i=1,N
    write(*,'(a,i8,a,F16.8)')"Eigenvalue ",i,": ",eigenvalues(i)
  ENDDO 
  write(*,*)

  ! Calculate the eigenvalues of B
  write(*,*)"Eigenvalues of real matrix B:"
  call diag_lapack_real(B,N,eigenvalues,eigenvectors)
  do i=1,N
    write(*,'(a,i8,a,F16.8)')"Eigenvalue ",i,": ",eigenvalues(i)
  ENDDO 
  write(*,*)

  ! The LAPACK calls overwrite A and B, so reset them before the next call
  A(1,1)=1.2d0; A(1,2)=3.4d0
  A(2,1)=5.6d0; A(2,2)=7.8d0

  B(1,1)=9.0d0; B(1,2)=1.2d0
  B(2,1)=3.4d0; B(2,2)=5.6d0

  ! Calculate the generalized eigenvalues of A with respect to B
  write(*,*)"Generalized eigenvalues of real matrix A with respect to real matrix B"
  call diag_lapack_real_generalized(A,B,N,eigenvalues,eigenvectors)
  do i=1,N
    write(*,'(a,i8,a,F16.8)')"Eigenvalue ",i,": ",eigenvalues(i)
  ENDDO 
  write(*,*)

  deallocate(A,B,eigenvalues,eigenvectors)
END SUBROUTINE test_real

SUBROUTINE test_complex
  IMPLICIT NONE 
  INTEGER  :: N,i
  complex*16,allocatable :: A(:,:),B(:,:),eigenvalues(:),eigenvectors(:,:)
  complex*16,parameter :: zi=(0.d0,1.d0)

  N=2 ! A and B are NxN matrices

  allocate(A(N,N),B(N,N),eigenvalues(N),eigenvectors(N,N))

  ! Set the values for the A and B matrices
  A(1,1)=1.2d0; A(1,2)=3.4d0
  A(2,1)=zi*5.6d0; A(2,2)=7.8d0

  B(1,1)=9.0d0; B(1,2)=1.2d0
  B(2,1)=3.4d0; B(2,2)=zi*5.6d0

  ! Calculate the eigenvalues of A
  write(*,*)"Eigenvalues of complex matrix A:"
  call diag_lapack_complex(A,N,eigenvalues,eigenvectors)
  do i=1,N
    write(*,'(a,i8,a,2F16.8)')"Eigenvalue ",i,": ",eigenvalues(i)
  ENDDO 
  write(*,*)

  ! Calculate the eigenvalues of B
  write(*,*)"Eigenvalues of complex matrix B:"
  call diag_lapack_complex(B,N,eigenvalues,eigenvectors)
  do i=1,N
    write(*,'(a,i8,a,2F16.8)')"Eigenvalue ",i,": ",eigenvalues(i)
  ENDDO 
  write(*,*)

  ! The LAPACK calls overwrite A and B, so reset them before the next call
  A(1,1)=1.2d0; A(1,2)=3.4d0
  A(2,1)=zi*5.6d0; A(2,2)=7.8d0

  B(1,1)=9.0d0; B(1,2)=1.2d0
  B(2,1)=3.4d0; B(2,2)=zi*5.6d0

  ! Calculate the generalized eigenvalues of A with respect to B
  write(*,*)"Generalized eigenvalues of complex matrix A with respect to complex matrix B"
  call diag_lapack_complex_generalized(A,B,N,eigenvalues,eigenvectors)
  do i=1,N
    write(*,'(a,i8,a,2F16.8)')"Eigenvalue ",i,": ",eigenvalues(i)
  ENDDO 
  write(*,*)

  deallocate(A,B,eigenvalues,eigenvectors)
END SUBROUTINE test_complex

SUBROUTINE diag_lapack_real(hamiltonian,N,eigenvalues,eigenvectors)
  IMPLICIT NONE 
  INTEGER  :: N,swapped,i,LWORK,retval
  real*8  :: hamiltonian(N,N),eigenvectors(N,N),eigenvalues(N),temp_eigenvalue
  real*8,allocatable :: ev_real(:),ev_imag(:),temp(:,:),temp_eigenvector(:),work(:)

  ! This first call just finds the optimal size for "work"
  allocate(ev_real(N),ev_imag(N),temp(1,N),temp_eigenvector(N))
  call dgeev('N','V',N,hamiltonian,N,ev_real,ev_imag,temp,1,eigenvectors,N,ev_real,-1,retval)
  LWORK=ev_real(1)
  allocate(work(LWORK))

  ! Now call to actually diagonalize
  call dgeev('N','V',N,hamiltonian,N,ev_real,ev_imag,temp,1,eigenvectors,N,work,LWORK,retval)
  eigenvalues=ev_real!+zi*ev_imag
  deallocate(ev_real,ev_imag,temp,work)

  ! Sort the results by eigenvalue
  swapped=1
  do while(swapped==1)
    swapped=0
    do i=1,N-1
      if(eigenvalues(i)>eigenvalues(i+1)) then
        temp_eigenvalue=eigenvalues(i)
        temp_eigenvector=eigenvectors(:,i)
        eigenvalues(i)=eigenvalues(i+1)
        eigenvectors(:,i)=eigenvectors(:,i+1)
        eigenvalues(i+1)=temp_eigenvalue
        eigenvectors(:,i+1)=temp_eigenvector
        swapped=1
      endif
    ENDDO 
  ENDDO 
  deallocate(temp_eigenvector)
END SUBROUTINE diag_lapack_real

SUBROUTINE diag_lapack_complex(hamiltonian,N,eigenvalues,eigenvectors)
  IMPLICIT NONE 
  INTEGER     :: N,swapped,i,LWORK,retval
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
    ENDDO 
  ENDDO 
  deallocate(temp_eigenvector)
END SUBROUTINE diag_lapack_complex

SUBROUTINE diag_lapack_real_generalized(A,B,N,eigenvalues,eigenvectors)
  IMPLICIT NONE 
  INTEGER  :: N,swapped,i,LWORK,retval
  real*8  :: A(N,N),B(N,N),eigenvectors(N,N),eigenvalues(N),temp_eigenvalue
  real*8,allocatable    :: ev_real(:),ev_imag(:),temp(:,:),temp_eigenvector(:),work(:),beta(:)
  complex*16, parameter :: zi=(0.d0,1.d0)

  ! This first call just finds the optimal size for "work"
  allocate(ev_real(N),ev_imag(N),temp(1,N),temp_eigenvector(N),beta(N))
  call dggev('N','V',N,A,N,B,N,ev_real,ev_imag,beta,temp,1,eigenvectors,N,temp_eigenvector,-1,retval)
  LWORK=temp_eigenvector(1)
  allocate(work(LWORK))

  ! Now call to actually diagonalize
  call dggev('N','V',N,A,N,B,N,ev_real,ev_imag,beta,temp,1,eigenvectors,N,work,LWORK,retval)
  eigenvalues=real((ev_real+zi*ev_imag)/beta)
  deallocate(ev_real,ev_imag,temp,work,beta)

  ! Sort the results by eigenvalue
  swapped=1
  do while(swapped==1)
    swapped=0
    do i=1,N-1
      if(eigenvalues(i)>eigenvalues(i+1)) then
        temp_eigenvalue=eigenvalues(i)
        temp_eigenvector=eigenvectors(:,i)
        eigenvalues(i)=eigenvalues(i+1)
        eigenvectors(:,i)=eigenvectors(:,i+1)
        eigenvalues(i+1)=temp_eigenvalue
        eigenvectors(:,i+1)=temp_eigenvector
        swapped=1
      endif
    ENDDO 
  ENDDO 
  deallocate(temp_eigenvector)
END SUBROUTINE diag_lapack_real_generalized

SUBROUTINE diag_lapack_complex_generalized(A,B,N,eigenvalues,eigenvectors)
  IMPLICIT NONE 
  INTEGER     :: N,swapped,i,LWORK,retval
  real*8     :: rwork(8*N)
  complex*16 :: A(N,N),B(N,N),eigenvectors(N,N),eigenvalues(N),temp_eigenvalue
  complex*16,allocatable :: temp(:,:),temp_eigenvector(:),work(:),beta(:)

  ! This first call just finds the optimal size for "work"
  allocate(temp(1,N),temp_eigenvector(N),beta(N))
  call zggev('N','V',N,A,N,B,N,eigenvalues,beta,temp,1,eigenvectors,N,temp_eigenvector,-1,rwork,retval)
  LWORK=int(real(temp_eigenvector(1)))
  allocate(work(LWORK))

  ! Now call to actually diagonalize
  call zggev('N','V',N,A,N,B,N,eigenvalues,beta,temp,1,eigenvectors,N,work,LWORK,rwork,retval)
  eigenvalues=eigenvalues/beta
  deallocate(temp,work,beta)

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
    ENDDO 
  ENDDO 
  deallocate(temp_eigenvector)
END SUBROUTINE diag_lapack_complex_generalized

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
        if(jstack.gt.NSTACK)pause 'NSTACK too small in indexx'
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
      END SUBROUTINE indexx
      

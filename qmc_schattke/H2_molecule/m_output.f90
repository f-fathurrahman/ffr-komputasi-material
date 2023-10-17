      module output 
       use highlevel
       use midlevel
C  output for 1-,2-,3-dimensional arrays on files termed PRONAME
C  with some additional naming
       implicit none
       public :: OUT1D,OUT2D,OUT3D
       integer,parameter,public :: FE1MAX=100,FE2MAX=100,FE3MAX=100
       integer :: ios
       character(10),public :: PRONAME 
      contains 
C
       subroutine OUT1D(NDIV1,FELD1D,FELDNAM)
        integer :: NDIV1
        real(kind=dp),intent(in),dimension(FE1MAX):: FELD1D
        character(7),intent(in) :: FELDNAM
        integer :: nx
        open(unit=36,file=PRONAME//FELDNAM//".dat",status="unknown")
        do nx=1,NDIV1
         write(36,*) nx,FELD1D(nx)
        end do
        close(36,iostat=ios)
        if (ios .gt. 0) then
          write(35,*) 'output error from OUT3D'
        end if
       end subroutine OUT1D
C
       subroutine OUT2D(FELD2D)
        real(kind=dp),intent(in),dimension(FE2MAX,FE2MAX):: FELD2D
        open(unit=36,file=PRONAME//"_erw_xfarbe.dat",status="unknown")
        open(unit=37,file=PRONAME//"_var_xfarbe.dat",status="unknown")
C        write(36,*) 'Energy expectation value on (size,alpha) - plane'
C        write(36,*) NXFA2+1,'  ',NXFA1+1
C        write(37,*) 'Variance on (size,alpha) - plane'
C        write(37,*) NXFA2+1,'  ',NXFA1+1
C
C        write(36,*) LENGTH0,LENGTH1
C        write(36,*) ALPHA0,ALPHA1
C        write(37,*) LENGTH0,LENGTH1
C        write(37,*) ALPHA0,ALPHA1
        close(36,iostat=ios)
        close(37,iostat=ios)
       end subroutine OUT2D
C
       subroutine OUT3D(NDIV1,NDIV2,NDIV3,FELD3D)
        integer :: NDIV1,NDIV2,NDIV3
        real(kind=dp),intent(in),dimension(FE3MAX,FE3MAX,FE3MAX):: 
     &          FELD3D
        integer :: nx,ny,nz
        open(unit=36,file=PRONAME//"_erw.dat",position="append",
     &               status="unknown")
        do nz=1,NDIV3
         do ny=1,NDIV2
          do nx=1,NDIV1
           write(36,55) nx,ny,nz,FELD3D(nx,ny,nz)
          end do
         end do
         write(36,*)'  '
         write(36,*)'  '
        end do
        close(36,iostat=ios)
        if (ios .gt. 0) then
          write(35,*) 'output error from OUT3D'
        end if
55     format(t3,i3,i3,i3,e12.3)
       end subroutine OUT3D
      end module output
C
C
C-----------------------------------------------------------------------       

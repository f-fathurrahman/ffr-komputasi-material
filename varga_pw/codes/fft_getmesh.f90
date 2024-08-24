  !
  ! $Id: fft_getmesh.f90,v 1.3 1999/03/25 13:31:41 noe Exp $
  !
  ! gets a list of allowed mesh points and dimension paddings, either from
  ! file or by using the default settings
  !
  SUBROUTINE getmesh

    INTEGER  :: nm
    logical :: called = .false.
    save    :: called

    if( .not. called ) then
       called = .true.

       open(39,file='fftmesh.dat',status='old',err=999)
       read(39,*) fftinc1, fftinc2, fftinc3       ! read mesh pads
       read(39,*) nm                              ! number of meshpoints
       allocate( fftmesh(nm) )                  ! alloc storage for meshpts
       do i=1,nm
          read(39,*) fftmesh(i)
       ENDDO 
       read(39,fmt=*,err=998,end=998) fft_conf
 998   close(39)
       write(6,*) 'read fft mesh data from ''fftmesh.dat'''
       return

 999   fftmesh => fft_defaultmesh
       return
    endif
  end SUBROUTINE



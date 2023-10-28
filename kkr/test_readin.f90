!------------------------------------------------
program test
!------------------------------------------------

implicit none

integer, parameter :: ngmx=15, msr=400 ,nbuf=1000

real(8) :: a, alpha, beta, boa
real(8) :: coa, edelt
real(8) :: ewidth, gamma

! i, j, k, l, m, n <- integers
integer :: maxitr
integer :: natm, ntyp

real(8) :: r(3,3),angl(3)
real(8), allocatable :: xrmt(:), xfield(:), xconc(:), xanclr(:)

integer, allocatable :: jlmxtyp(:), jncmp(:)
character  go*6,file*256,brvtyp*6,reltyp*6,sdftyp*12,magtyp*4 &
   &          ,outtyp*6,bzqlty*8,record*4, pmxtyp*16
character, allocatable:: xtype(:)*8, xatmtyp(:)*8, xatmicv(:)*24


allocate( xanclr(nbuf) )
allocate( xrmt(nbuf) )
allocate( xfield(nbuf) )
allocate( xconc(nbuf) )
allocate( jlmxtyp(nbuf) )
allocate( jncmp(nbuf) )
allocate( xtype(nbuf) )
allocate( xatmtyp(nbuf) )
allocate( xatmicv(3*nbuf) )

call readin(go,file,brvtyp,a,coa,boa,alpha,beta,gamma &
  &  ,edelt,ewidth,reltyp,sdftyp,magtyp,record,outtyp &
  &  ,bzqlty,maxitr,pmxtyp,ntyp,xtype,jncmp,xrmt,xfield &
  &  ,jlmxtyp,xanclr,xconc,natm,xatmicv,xatmtyp,r,angl,*20)



deallocate( xanclr )
deallocate( xrmt )
deallocate( xfield )
deallocate( xconc )
deallocate( jlmxtyp )
deallocate( jncmp )
deallocate( xtype )
deallocate( xatmtyp )
deallocate( xatmicv )


20 continue
write(*,*)
write(*,*) 'Program finished normally'
write(*,*) '-------------------------'
write(*,*)

end program


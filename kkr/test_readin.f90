!------------------------------------------------
program test
!------------------------------------------------
!
implicit none
!
integer, parameter :: ngmx=15, msr=400 ,nbuf=1000
!
real(8) :: a, alpha, beta, boa
real(8) :: coa, edelt
real(8) :: ewidth, gamma
!
! i, j, k, l, m, n <- integers
integer :: maxitr
integer :: natm, ntyp
!
real(8) :: r(3,3),angl(3)
real(8), allocatable :: xrmt(:), xfield(:), xconc(:), xanclr(:)
!
integer, allocatable :: jlmxtyp(:), jncmp(:)
character :: go*6, file*256, brvtyp*6, reltyp*6, sdftyp*12, magtyp*4
character :: outtyp*6, bzqlty*8, record*4, pmxtyp*16
character, allocatable:: xtype(:)*8, xatmtyp(:)*8, xatmicv(:)*24
!
integer :: i, nf, mse, mse0, lmpmx, msex, ids, ibrav
integer :: ncmpx, ncmax, ndmx, mxl, lastmx
character :: dmy*6, trmkey*16
! Functions
integer :: nfqlty
integer :: ibrava


allocate( xanclr(nbuf) )
allocate( xrmt(nbuf) )
allocate( xfield(nbuf) )
allocate( xconc(nbuf) )
allocate( jlmxtyp(nbuf) )
allocate( jncmp(nbuf) )
allocate( xtype(nbuf) )
allocate( xatmtyp(nbuf) )
allocate( xatmicv(3*nbuf) )

call readin(go, file, &
    brvtyp, a, coa, boa, alpha, beta, gamma, &
    edelt, ewidth, reltyp, sdftyp, magtyp, record, outtyp, &
    bzqlty, maxitr, pmxtyp, ntyp, xtype, jncmp, xrmt, xfield, &
    jlmxtyp, xanclr, xconc, natm, xatmicv, xatmtyp, r, angl, *20)

write(*,*) 'ntyp = ', ntyp
write(*,*) 'natm = ', natm
write(*,*) 'xconc = ', xconc(1:ntyp)

! These are strings
write(*,*) 'shape xatmtyp = ', shape(xatmtyp)
write(*,*) 'xatmtyp(1) = ', xatmtyp(1)
write(*,*)
write(*,*) 'shape xtype = ', shape(xtype)
write(*,*) 'xtype(1) = ', xtype(1)
write(*,*)
write(*,*) 'shape xatmicv = ', shape(xatmicv)
write(*,*) 'xatmicv(1) = ', xatmicv(1)
write(*,*)
write(*,*) 'shape jncmp = ', shape(jncmp)
do i = 1,ntyp
  write(*,'(1x,A,I4,I4)') 'i, jncmp(i) = ', i, jncmp(i)
enddo
write(*,*)
write(*,*) 'shape jlmxtyp = ', shape(jlmxtyp)
do i = 1,ntyp
  write(*,'(1x,A,I4,I4)') 'i, jlmxtyp(i) = ', i, jlmxtyp(i)
enddo

ncmpx = 0
ncmax = 0
mxl = 1
do i = 1,ntyp
  mxl = max(mxl, jlmxtyp(i) + 1)
  ncmax = max(ncmax, jncmp(i))
  ncmpx = ncmpx + jncmp(i)
enddo

write(*,*)
write(*,*) 'mxl = ', mxl
write(*,*) 'ncmax = ', ncmax
write(*,*) 'ncmpx = ', ncmpx

ndmx = natm*(natm - 1) + 1
write(*,*) 'ndmx = ', ndmx
      
lastmx = int(2900d0/3d0**(5-mxl))
write(*,*) 'lastmx = ', lastmx

! data
mse0 = 35
msex = 201

! what's this? For LMD calculation?
lmpmx = ncmpx*(ncmpx - 1)/2
write(*,*) 'lmpmx = ', lmpmx

! DOS calculation flag, used in drvmsh
ids = 0
! ids will be set according to value of go (type of calculation)

! mse is no. of energy mesh
! it will be set to msex for DOS calculation
mse = min(mse0, msex)
write(*,*) 'mse = ', mse
! Many arrays have mse in their dimensioning

dmy = brvtyp
ibrav = ibrava(trmkey('tlt', dmy))
write(*,*) 'dmy = ', dmy
write(*,*) 'ibrav = ', ibrav

nf = nfqlty(bzqlty, ibrav)
write(*,*) 'nf = ', nf

!
! Do some conversion (?)
!
allocate( anclr(ncmpx), rmt(ntyp), field(ntyp), conc(ncmpx), &
          lmxtyp(ntyp), ncmp(ntyp), type(ntyp), atmtyp(natm), &
          atmicv(3*natm) )

call equarr(xrmt,rmt,ntyp)
call equarr(xfield,field,ntyp)
call equarr(xconc,conc,ncmpx)
      
    !call equarr(xtype,type,ntyp)
    type(1:ntyp) = xtype(1:ntyp)

    call equari(jncmp,ncmp,ntyp)
    call equarr(xanclr,anclr,ncmpx)
    call equari(jlmxtyp,lmxtyp,ntyp)
      
    !call equarr(xatmtyp,atmtyp,natm)
    atmtyp(1:natm) = xatmtyp(1:natm)

    !write(*,*) 'xatmicv = ', xatmicv
    !write(*,*) 'atmicv = ', atmicv

    call equarr(xatmicv,atmicv,9*natm)
    !atmicv(1:9*natm) = xatmicv(1:9*natm)
    !atmicv(:) = xatmicv(:)




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


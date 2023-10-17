C-----------------------------------------------------------------------
      module determinant
       use highlevel
       use midlevel
       use orbital
C   Lithium version
C   4.10.2009: Uses new module strategy
C   Here: update of determinant and its derivatives
C   QD(2)= determinantal acceptance ratio for both spins
C   DET= actual determinant
C   ANEW= matrix inverse referring to the new determinant for each spin
C   AOLD= matrix inverse referring to the old determinant for each spin
       implicit none
       real(dp),public,dimension(2)       :: QD,QDOPT,DET,DETOPT
       real(dp),public,dimension(NE,2)    :: LAPLDET,LAPLDETOLD
       real(dp),public,dimension(3,NE,2)  :: GRADDET,GRADDETOLD
       real(dp),dimension(NE,NE,2),public :: ANEW,AOLD,ANEWOPT,AOLDOPT
       public :: SLAQUOT,SLASM,SLAKIN,DETARR

      contains
C-----------------------------------------------------------------------
      subroutine INITDET
       integer  :: j,hs
       hs=NES(1)
       AOLD=0.0_dp
       ANEW=0.0_dp
       AOLDOPT=0.0_dp
       ANEWOPT=0.0_dp
       call DETARR(hs,PSIMAT(1:hs,1:hs,1),DET(1))
       call DETARR(hs,PSIMATOPT(1:hs,1:hs,1),DETOPT(1))
       do j = 1,NE
        AOLD(j,j,1:2) = 1.0_dp
        AOLD(1,2,1)=0.11111_dp
        ANEW(j,j,1:2) = 1.0_dp
        AOLDOPT(j,j,1:2) = 1.0_dp
        AOLDOPT(1,2,1)=0.11111_dp
        ANEWOPT(j,j,1:2) = 1.0_dp
       end do
       AOLD(1,1,2)=1.0_dp/PSIMAT(1,1,2)
       AOLDOPT(1,1,2)=1.0_dp/PSIMATOPT(1,1,2)
      end subroutine INITDET
C-----------------------------------------------------------------------
      subroutine SLAQUOT(aold,psinew,q)
C   Calculates quotient q, the ratio between new and old determinant,
C   as determinantal acceptance ratio
C   for an offered move of electron IEES of spin IES, calculated by
C   decomposition with respect to the cofactors
       real(dp),intent(in),dimension(NES(IES),NES(IES))   :: aold
       real(dp),intent(in),dimension(NES(IES),NES(IES))   :: psinew
       real(dp),intent(out)                         :: q
       integer                                      :: ns
C
       ns = NES(IES)
       q=dot_product (aold(IEES,1:ns),psinew(1:ns,IEES))
      end subroutine slaquot

C-----------------------------------------------------------------------
      subroutine SLASM(q,aold,psinew,anew)
C   Calculates according to the Sherman-Morrison formula the
C   change of the inverse of a matrix, if the global electron
C   variable IE has been moved.
C   q=  ratio new-determinant/old-determinant
C   aold= matrix inverse of the old-determinant matrix before move
C   psinew(i,k)=  new wavefunction matrix after move
C   first index refers to orbital, second index to electron number
C   anew=  matrix inverse of the new-determinant matrix after move
       real(dp),intent(in)                          :: q
       real(dp),intent(in),dimension(NES(IES),NES(IES))   :: aold
       real(dp),intent(in),dimension(NES(IES),NES(IES))  :: psinew
       real(dp),intent(out),dimension(NES(IES),NES(IES))  :: anew
C
       integer                                      :: j,ns
       real(dp),dimension(NES(IES))                    :: vsum
C   The IEES'th electron of spin IES is actually considered
       ns = NES(IES)
       if (q == 0.0_dp) then
        anew = 0.0_dp
        write (*,*) 'QD=0, node of determinant'
       else
        do j=1,ns
         vsum(j) = dot_product (aold(j,1:ns),psinew(1:ns,IEES))
        end do
        vsum =-vsum/q
        do j=1,ns
         anew(1:ns,j) = aold(1:ns,j) +
     &                 vsum(1:ns)*aold(IEES,j)
        end do
        anew(IEES,1:ns) = aold(IEES,1:ns)/q
       end if
      end subroutine SLASM
C-----------------------------------------------------------------------
      subroutine SLAKIN(aold,pgrnew,planew,graddet,lapldet)
C  Calculates the contributions of the Slater determinant to the local
C  kinetic energy, i.e. gradient and laplacian, of electron IE
C  with spin  IES in one-eletron wavefunction psi(orbital). ORBDER
C  has to be called before in order to obtain the input, pgrnew and
C  planew, at the new positions of electron IE
C  aold= inverse of determinant matrix
C  pgrnew= gradient of one-electron wavefunction at new position
C  planew= laplacian of one-electron wavefunction at new position
       real(dp),intent(in),dimension(NES(IES),NES(IES)) :: aold
       real(dp),intent(in),dimension(3,NORB)            :: pgrnew
       real(dp),intent(in),dimension(NORB)              :: planew
       real(dp),intent(out),dimension(3,NES(IES))       :: graddet
       real(dp),intent(out),dimension(NES(IES))         :: lapldet
C
       integer                     :: js,ns,ns0
       real(dp)                    :: hpla
       real(dp),dimension(3)       :: hpgr
       logical                     :: toolarge
       ns = NES(IES)
       ns0 = (IES-1)*NES(1)
       hpgr = 0.0_dp
       hpla = 0.0_dp
       do js=1,ns
        hpgr(1:3) = hpgr(1:3) +
     &              aold(IEES,js)*pgrnew(1:3,NELORB(ns0+js))
        hpla = hpla + aold(IEES,js)*planew(NELORB(ns0+js))
       end do
       graddet(1:3,IEES) = hpgr(1:3)
       lapldet(IEES) = hpla
       toolarge = .false.
       if ((sum(dabs(hpgr)) > 1000.0_dp) .OR.
     &         (dabs(hpla) > 1000.0_dp)) toolarge = .true.
       if (toolarge) write(39,*) 'kinetic energy large grad,lapl=',
     &              abs(sum(hpgr)),hpla
      end subroutine SLAKIN
C-----------------------------------------------------------------------
      subroutine DETARR(nies,psiarr,det)
C  Calculates the determinant
C  This routine is mainly for test purposes, as the determinant is
C  generated at start by random filling and always updated with help of
C  the inverse. Here, we start with unit matrix c(i,k) and determinant,
C  and update similarly by introducing at each step a column of the
C  matrix psiarr(i,k) via the inverse aold (i,k) of updated c(i,k).
       integer,intent(in)  :: nies
       real(dp),intent(in),dimension(nies,nies) :: psiarr
       real(dp),intent(out)  :: det
C
       integer   :: m,n,j,k,hnes
       real(dp)  :: q,dold,dnew
       real(dp),dimension(nies)          :: ha
       real(dp),dimension(nies,nies) :: aold,anew,c
       hnes = nies
       c = 0.0_dp
       do j=1,hnes
        c(j,j)=1.0_dp
       end do
       dold=1.0_dp
       aold=c
       q=0.0_dp
       do k=1,hnes
        q = dot_product (aold(k,1:hnes),psiarr(1:hnes,k))
        do m=1,hnes
         ha(m)=dot_product (aold(m,1:hnes),psiarr(1:hnes,k))
         if (m == k) then
          anew(m,1:hnes)=aold(m,1:hnes)/q
          else
          anew(m,1:hnes)=aold(m,1:hnes)-ha(m)*aold(k,1:hnes)/q
         end if
        end do
        dnew=q*dold
        aold=anew
        dold=dnew
       end do
       det = dnew
      end subroutine DETARR
C-----------------------------------------------------------------------
      end module determinant
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

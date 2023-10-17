C-----------------------------------------------------------------------
      module jastrow
       use highlevel
       use midlevel
C  Calculates the change of the Jastrow two-particle potential and
C  updates the derivatives. The Jastrow factor J is written as
C  J=exp(-\sum_k u) with u as in the subroutines.
C  GRJAS(3,NE,IES)=  Gradient of J divided by J for spin IES
C  LAPJAS(NE,IES)= Laplace of J divided by J for spin IES
C  QJ(IES)= J^new/J^old, Jastrow part of acceptance for spin IES,
C  to be squared to yield the rates
C  calculated for actual electron, IE.
C  CJAS= only Jastrow constant for exponential like Jastrow factor
C  BETA1,BETA2= Additional constants for atomic like Jastrow factor,
C  which introduce a short range repulsion (attraction) for
C  positive (negative) values. Set BETA1=0.0 to cancel the latter part.
       implicit none
       public :: INITJAS,JASEXP,JASEXPATOM,JASEXPATOM1
       real(kind=dp),public                   :: CJAS,BETA1,BETA2
       real(kind=dp),dimension(2),public      :: QJ,QJOPT
       real(kind=dp),public              :: CJAS1,CJASOPT,WJAS,WJASOPT
       real(kind=dp),dimension(NE),public     :: JPSI,JPSIOPT
       real(kind=dp),dimension(NE,2),public   :: LAPJAS,LAPJASOLD
       real(kind=dp),dimension(3,NE,2),public :: GRJAS,GRJASOLD
      contains
C-----------------------------------------------------------------------
      subroutine INITJAS
       BETA1=0.0_dp
       BETA2=0.0_dp
       if (CJAS < EMACH) CJAS=EMACH
      end subroutine INITJAS
C-----------------------------------------------------------------------
       subroutine JASEXP
C  Determines Jastrow exponent and derivatives for a form
C  u=as*F/r*(1-e^(-r/F)), as=0.5 for parallel spin and
C  as=1.0 for antiparallel spins. From the derivatives
C  the kinetic energy contribution of the Jastrow factor
C  is obtained.
        integer :: i,k,n
        real(dp) :: as,jasu,jasdif,u2d,u2do,woo,won
        real(dp),dimension(3) :: u1d,u1do
        jasu=0.0_dp
        jasdif=0.0_dp
        u1d=0.0_dp
        u1do=0.0_dp
        u2d=0.0_dp
        u2do=0.0_dp
        call RDIST(RE,RENEW)
        ielek:do k=1,NE
         if (k .eq. IE) then
          cycle ielek
         end if
         as=0.5_dp     ! for equal spins
         if (((IES == 1) .and. (k > NES(1))) .or.
     &       ((IES == 2) .and. (k <= NES(1)))) as =1.0_dp
         woo = DIST(4,IE,k)
         won = DISTNEW(4,IE,k)
         jasdif=jasdif+(as/won*(1.0_dp-dexp(-won/CJAS))-
     &         as/woo*(1.0_dp-dexp(-woo/CJAS)))*CJAS**2
         jasu=jasu+as/won*(1.0_dp-dexp(-won/CJAS))*CJAS**2
         u1d(1:3)=u1d(1:3)-as*DISTNEW(1:3,IE,k)/won**3*(1.0_dp-
     &       dexp(-won/CJAS)*(1.0_dp+won/CJAS))*CJAS**2
         u1do(1:3)=u1do(1:3)-as*DIST(1:3,IE,k)/woo**3*(1.0_dp-
     &       dexp(-woo/CJAS)*(1.0_dp+woo/CJAS))*CJAS**2
         u2d=u2d-as/won*dexp(-won/CJAS)
         u2do=u2do-as/woo*dexp(-woo/CJAS)
        end do ielek
        QJ(IES) = dexp(-jasdif)
        GRJAS(1:3,IEES,IES)=-u1d(1:3)
        GRJASOLD(1:3,IEES,IES)=-u1do(1:3)
        LAPJAS(IEES,IES)=-u2d+(u1d(1)**2+u1d(2)**2+u1d(3)**2)
        LAPJASOLD(IEES,IES)=-u2do+(u1do(1)**2+u1do(2)**2+u1do(3)**2)
      end subroutine JASEXP
C-----------------------------------------------------------------------
      subroutine JASEXPATOM
C  Determines Jastrow exponent and derivatives for a form
C  u=F/2/(1+r/F)*(delta_{s,-s'}+1/2*delta_{s,s'})+
C    BETA_1exp(-BETA2*r^2). From the derivatives
C  the kinetic energy contribution of the Jastrow factor
C  is obtained.
        integer :: i,k,n
        real(dp) :: as,jasu,jasdif,u2d,u2do,woo,won
        real(dp),dimension(3) :: u1d,u1do
        jasu=0.0_dp
        jasdif=0.0_dp
        u1d=0.0_dp
        u1do=0.0_dp
        u2d=0.0_dp
        u2do=0.0_dp
        call RDIST(RE,RENEW)
        ielekat:do k=1,NE
         if (k .eq. IE) then
          cycle ielekat
         end if
         as=0.5_dp     ! for equal spins
         if (((IES == 1) .and. (k > NES(1))) .or.
     &       ((IES == 2) .and. (k <= NES(1)))) as =1.0_dp
         woo = DIST(4,IE,k)
         won = DISTNEW(4,IE,k)
         jasdif=jasdif+(as/2.0_dp/(1.0_dp+won/CJAS)-
     &       as/2.0_dp/(1.0_dp+woo/CJAS))*CJAS+
     &       BETA1*(dexp(-BETA2*won**2)-dexp(-BETA2*woo**2))
         jasu=jasu+as/2.0_dp/(1.0_dp+won/CJAS)*CJAS+
     &       BETA1*dexp(-BETA2*won**2)
         u1d(1:3)=u1d(1:3)-as*DISTNEW(1:3,IE,k)/won/
     &            (1.0_dp+won/CJAS)**2/2.0_dp-
     &       2.0_dp*DISTNEW(1:3,IE,k)*BETA1*BETA2*dexp(-BETA2*won**2)
         u1do(1:3)=u1do(1:3)-as*DIST(1:3,IE,k)/woo/
     &            (1.0_dp+woo/CJAS)**2/2.0_dp-
     &       2.0_dp*DIST(1:3,IE,k)*BETA1*BETA2*dexp(-BETA2*woo**2)
         u2d=u2d-as/won/(1.0_dp+won/CJAS)**3+
     &    BETA1*BETA2*dexp(-BETA2*won**2)*(4.0_dp*BETA2*won**2-6.0_dp)
         u2do=u2do-as/woo/(1.0_dp+woo/CJAS)**3+
     &    BETA1*BETA2*dexp(-BETA2*woo**2)*(4.0_dp*BETA2*woo**2-6.0_dp)
        end do ielekat
        JPSI(IE) = dexp(-jasu)
        QJ(IES) = dexp(-jasdif)
        GRJAS(1:3,IEES,IES)=-u1d(1:3)
        GRJASOLD(1:3,IEES,IES)=-u1do(1:3)
        LAPJAS(IEES,IES)=-u2d+(u1d(1)**2+u1d(2)**2+u1d(3)**2)
        LAPJASOLD(IEES,IES)=-u2do+(u1do(1)**2+u1do(2)**2+u1do(3)**2)
       end subroutine JASEXPATOM
C-----------------------------------------------------------------------
      subroutine JASEXPATOM1
C  Determines Jastrow exponent and derivatives for a form
C  u=F/2/(1+r/F)*(delta_{s,-s'}+1/2*delta_{s,s'})+
C    BETA_1exp(-BETA2*r^2). From the derivatives
C  the kinetic energy contribution of the Jastrow factor
C  is obtained.
        integer :: i,k,n
        real(dp) :: as,jasu,jasdif,u2d,u2do,woo,won
        real(dp),dimension(3) :: u1d,u1do
        jasu=0.0_dp
        jasdif=0.0_dp
        u1d=0.0_dp
        u1do=0.0_dp
        u2d=0.0_dp
        u2do=0.0_dp
        call RDIST(RE,RENEW)
        CJAS1 = CJASOPT
        ielekat:do k=1,NE
         if (k .eq. IE) then
          cycle ielekat
         end if
         as=0.5_dp     ! for equal spins
         if (((IES == 1) .and. (k > NES(1))) .or.
     &       ((IES == 2) .and. (k <= NES(1)))) as =1.0_dp
         woo = DIST(4,IE,k)
         won = DISTNEW(4,IE,k)
         jasdif=jasdif+(as/2.0_dp/(1.0_dp+won/CJAS1)-
     &       as/2.0_dp/(1.0_dp+woo/CJAS1))*CJAS1+
     &       BETA1*(dexp(-BETA2*won**2)-dexp(-BETA2*woo**2))
         jasu=jasu+as/2.0_dp/(1.0_dp+won/CJAS1)*CJAS1+
     &       BETA1*dexp(-BETA2*won**2)
         u1d(1:3)=u1d(1:3)-as*DISTNEW(1:3,IE,k)/won/
     &            (1.0_dp+won/CJAS1)**2/2.0_dp-
     &       2.0_dp*DISTNEW(1:3,IE,k)*BETA1*BETA2*dexp(-BETA2*won**2)
         u1do(1:3)=u1do(1:3)-as*DIST(1:3,IE,k)/woo/
     &            (1.0_dp+woo/CJAS1)**2/2.0_dp-
     &       2.0_dp*DIST(1:3,IE,k)*BETA1*BETA2*dexp(-BETA2*woo**2)
         u2d=u2d-as/won/(1.0_dp+won/CJAS1)**3+
     &    BETA1*BETA2*dexp(-BETA2*won**2)*(4.0_dp*BETA2*won**2-6.0_dp)
         u2do=u2do-as/woo/(1.0_dp+woo/CJAS1)**3+
     &    BETA1*BETA2*dexp(-BETA2*woo**2)*(4.0_dp*BETA2*woo**2-6.0_dp)
        end do ielekat
        JPSIOPT(IE) = dexp(-jasu)
        QJOPT(IES) = dexp(-jasdif)
C        GRJAS(1:3,IEES,IES)=-u1d(1:3)
C        GRJASOLD(1:3,IEES,IES)=-u1do(1:3)
C        LAPJAS(IEES,IES)=-u2d+(u1d(1)**2+u1d(2)**2+u1d(3)**2)
C        LAPJASOLD(IEES,IES)=-u2do+(u1do(1)**2+u1do(2)**2+u1do(3)**2)
       end subroutine JASEXPATOM1
C-----------------------------------------------------------------------
      end module jastrow
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------

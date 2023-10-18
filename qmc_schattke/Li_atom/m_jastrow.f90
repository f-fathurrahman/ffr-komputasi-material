module m_jastrow
  use m_highlevel
  use m_midlevel

!  Calculates the change of the Jastrow two-particle potential and
!  updates the derivatives. The Jastrow factor J is written as
!  J=exp(-\sum_k u) with u as in the subroutines.
!  GRJAS(3,NE,IES)=  Gradient of J divided by J for spin IES
!  LAPJAS(NE,IES)= Laplace of J divided by J for spin IES
!  QJ(IES)= J^new/J^old, Jastrow part of acceptance for spin IES,
!  to be squared to yield the rates
!  calculated for actual electron, IE.
!  CJAS= only Jastrow constant for exponential like Jastrow factor
!  BETA1,BETA2= Additional constants for atomic like Jastrow factor,
!  which introduce a short range repulsion (attraction) for
!  positive (negative) values. Set BETA1=0.0 to cancel the latter part.
  
  implicit none
  public :: INITJAS,JASEXP,JASEXPATOM
  real(kind=dp),public:: CJAS,BETA1,BETA2,QJ(2)
  real(kind=dp),dimension(NE,2),public   :: LAPJAS,LAPJASOLD
  real(kind=dp),dimension(3,NE,2),public :: GRJAS,GRJASOLD

contains

subroutine INITJAS
  BETA1=0.0_dp
  BETA2=0.0_dp
  if (CJAS < EMACH) CJAS=EMACH
end subroutine INITJAS


subroutine JASEXP()
!  Determines Jastrow exponent and derivatives for a form
!  u=as*F/r*(1-e^(-r/F)), as=0.5 for parallel spin and
!  as=1.0 for antiparallel spins. From the derivatives
!  the kinetic energy contribution of the Jastrow factor
!  is obtained.
  integer :: k
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
         if (((IES == 1) .and. (k > NES(1))) .or. &
     &       ((IES == 2) .and. (k <= NES(1)))) as =1.0_dp
         woo = DIST(4,IE,k)
         won = DISTNEW(4,IE,k)

         jasdif=jasdif+(as/won*(1.0_dp-dexp(-won/CJAS))- &
     &         as/woo*(1.0_dp-dexp(-woo/CJAS)))*CJAS**2
         
         jasu=jasu+as/won*(1.0_dp-dexp(-won/CJAS))*CJAS**2
         
         u1d(1:3)=u1d(1:3)-as*DISTNEW(1:3,IE,k)/won**3*(1.0_dp - &
     &       dexp(-won/CJAS)*(1.0_dp+won/CJAS))*CJAS**2

         u1do(1:3)=u1do(1:3)-as*DIST(1:3,IE,k)/woo**3*(1.0_dp - &
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



subroutine JASEXPATOM()
!  Determines Jastrow exponent and derivatives for a form
!  u=F/2/(1+r/F)*(delta_{s,-s'}+1/2*delta_{s,s'})+
!    BETA_1exp(-BETA2*r^2). From the derivatives
!  the kinetic energy contribution of the Jastrow factor
!  is obtained.
        integer :: k
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
         
         as = 0.5_dp     ! for equal spins
         
         if (((IES == 1) .and. (k > NES(1))) .or. &
     &       ((IES == 2) .and. (k <= NES(1)))) as =1.0_dp
         
         woo = DIST(4,IE,k)
         won = DISTNEW(4,IE,k)

         jasdif=jasdif+(as/2.0_dp/(1.0_dp+won/CJAS)- &
     &       as/2.0_dp/(1.0_dp+woo/CJAS))*CJAS+ &
     &       BETA1*(dexp(-BETA2*won**2)-dexp(-BETA2*woo**2))
         
         jasu=jasu+as/2.0_dp/(1.0_dp+won/CJAS)*CJAS + &
     &       BETA1*dexp(-BETA2*won**2)
         
         u1d(1:3)=u1d(1:3)-as*DISTNEW(1:3,IE,k) / won / &
     &            (1.0_dp+won/CJAS)**2/2.0_dp - &
     &       2.0_dp*DISTNEW(1:3,IE,k)*BETA1*BETA2*dexp(-BETA2*won**2)
         
         u1do(1:3)=u1do(1:3)-as*DIST(1:3,IE,k) / woo / &
     &            (1.0_dp+woo/CJAS)**2/2.0_dp- &
     &       2.0_dp*DIST(1:3,IE,k)*BETA1*BETA2*dexp(-BETA2*woo**2)
         
         u2d=u2d-as/won/(1.0_dp+won/CJAS)**3 + &
     &    BETA1*BETA2*dexp(-BETA2*won**2)*(4.0_dp*BETA2*won**2-6.0_dp)
         
         u2do=u2do-as/woo/(1.0_dp+woo/CJAS)**3 + &
     &    BETA1*BETA2*dexp(-BETA2*woo**2)*(4.0_dp*BETA2*woo**2-6.0_dp)
        end do ielekat
        QJ(IES) = dexp(-jasdif)
        GRJAS(1:3,IEES,IES)=-u1d(1:3)
        GRJASOLD(1:3,IEES,IES)=-u1do(1:3)
        LAPJAS(IEES,IES)=-u2d+(u1d(1)**2+u1d(2)**2+u1d(3)**2)
        LAPJASOLD(IEES,IES)=-u2do+(u1do(1)**2+u1do(2)**2+u1do(3)**2)
       end subroutine JASEXPATOM



end module

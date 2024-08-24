MODULE XC
  implicit none
  real*8,parameter :: pi=3.141592653589793d0
  real*8,parameter :: Ry=13.60569193d0    !Rydberg constant in eV
  real*8,parameter :: a_B=0.52917720859d0 !Bohr radius
  real*8,parameter :: gammaU=-0.1423d0,beta1U=1.0529d0,beta2U=0.3334d0,AU=0.0311d0,BU=-0.048d0
  real*8,parameter :: CU=0.002d0,DU=-0.0116d0  
 
  CONTAINS


      subroutine xc_pot(rho_up,rho_down,v_x_up,v_x_down,&
&                                   v_c_up,v_c_down,eps_x,eps_c)
      
      implicit none
      real*8                           :: eps_x,eps_c,rho,rho_up,rho_down,&
&                                         v_x_up,v_x_down,v_c_up,v_c_down,&
&                                         v_x,v_c

!      rho=rho_up+rho_down
!      call  xc_lda(rho,v_x,v_c,eps_x,eps_c)
!      v_x_up=v_x
!      v_x_down=v_x
!      v_c_up=v_c
!      v_c_down=v_c

      call xc_lsda(rho_up,rho_down,&
&                         v_x_up,v_x_down,v_c_up,v_c_down,eps_x,eps_c)

      end  subroutine xc_pot

      subroutine xc_lda(rho,v_x,v_c,eps_x,eps_c)
      
      implicit none
      integer                                    :: i
      real*8                                     :: rho,c1,c2,rs,&
     &                                              rssq,rsln,v_x,v_c,&
     &                                              sum,eps_x,eps_c

      c1=3/(4*Pi)
      c2=4.d0/3.d0*0.4582d0
      sum=0.d0

        rs=(c1/rho)**(1.d0/3.d0)
        rs=rs/a_B
        v_x=-c2/rs
        if(rs.gt.1.d0) then
          rssq=sqrt(rs)
          v_c=gammaU*(1.d0+7.d0/6.d0*beta1U*rssq+4.d0/3.d0*beta2U*rs)&
     &             /(1.d0+beta1U*rssq+beta2U*rs)**2
        else
          rsln=log(rs)
          v_c=AU*rsln+(BU-AU/3.d0)+2.d0/3.d0*CU*rs*rsln+(2.d0*DU-CU)/3.d0*rs
        endif


        eps_x=-.4582d0/rs*2.d0*Ry
        if(rs.gt.1.d0) then
          rssq=sqrt(rs)
          eps_c=gammaU/(1.d0+beta1U*rssq+beta2U*rs)
        else
          rsln=log(rs)
          eps_c=AU*rsln+BU+CU*rs*rsln+DU*rs
        endif
        eps_c=2.d0*Ry*eps_c

      end subroutine xc_lda

      subroutine xc_lsda(rho_up,rho_down,&
&                         v_x_up,v_x_down,v_c_up,v_c_down,eps_x,eps_c)
! ======================================================================
!      .......
!      exc_tot=0.0
!      do (over grid)
!         ......
!         call exco_lda(...)
!         exc_tot = exc_tot + rho*(eps_x+eps_c)*weight
!         ......
!      enddo
!      ......
! USAGE :  For potentials
!      .......
!         call exco_lda(...)
!      FOR EACH INPUT THE SCALAR POTENTIALS ARE
!         v_x_up, v_x_down, v_c_up, v_c_down
!      .......
! **********************************************************************
! INPUT :
!     real up,dn    :  UP/DOWN DENSITY
! OUTPUT:
!     real eps_x   : LSDA EXCHANGE ENERGY PER PTL. 
!     real v_x_up : UP LSDA EXCHANGE POTENTIAL
!     real v_x_down : DOWN LSDA EXCHANGE POTENTIAL 
!     real eps_c   : LSDA CORRELATION ENERGY PER PTL.
!     real v_c_up : UP LSDA CORRELATION POTENTIAL
!     real v_c_down : DOWN LSDA CORRELATION POTENTIAL 
! LOCAL :
!     real rhomin   : MIN. CUTOFF OF DENSITY
!     real pi
!     real thrd,thrd2,thrd4
!     real conf,alpha,ax : rs=alpha/fk, e_x[LDA]=ax*rho^(4/3)  
!     real fk       : FERMI WAVEVECTOR  = conf*(rho)^(1/3) 
!     real rho2,rho
! EXCHANGE PART:
!     real exuplsda,exdnlsda
! CORRELATION PART:
!     real zet      : (up-dn)/rho
!     real rs       : LOCAL SEITZ RADIUS=(3/(4pi*rho))^(1/3)=alpha/fk
! **********************************************************************
      implicit none
      real*8 rho_up,rho_down,eps_x,v_x_up,v_x_down,eps_c,v_c_up,v_c_down,&
          rhomin,pi,thrd,thrd2,thrd4,conf,alpha,ax,fk,rho2,rho,&
          exuplsda,exdnlsda,zet,rs
!
      parameter(rhomin=1.e-20)  ! MACHINE DEPENDENT. CHECK!
      parameter(thrd=1./3.,thrd2=2./3.,thrd4=4./3.)
! ----------------------------------------------------------------------
      pi=4.*atan(1.0)
      conf=(3.*pi**2)**thrd     ! fk=conf*d^(1/3)
      alpha=(9.*pi/4.)**thrd    ! rs=alpha/fk
!
! --LSDA EXCHANGE--
!
! USE SPIN-SCALING RELATION: 
!     Ex[up,dn]=0.5*(Ex[2*up]+Ex[2*dn]) 
!   OR
!     ex[up,dn]=Ex/N=(ex[2*up]*up+ex[2*dn]*dn)/(up+dn)
! FORMULAS:
!     e_x[LDA]=ax*rho^(4/3)  
!   WHERE
!     ax = -3/4*(3/pi)^(1/3)
!
      ax=-0.75*(3./pi)**thrd
! DO UP EXCHANGE
      rho2=2.*rho_up
      if(rho2.gt.rhomin)then     ! AVOID OVERFLOW
        exuplsda=ax*rho2**thrd
        v_x_up=exuplsda*thrd4
      else
        exuplsda=0.0
        v_x_up=0.0
      endif
! REPEAT FOR DOWN
      rho2=2.*rho_down
      if(rho2.gt.rhomin)then
        exdnlsda=ax*rho2**thrd
        v_x_down=exdnlsda*thrd4
      else
	exdnlsda=0.0
	v_x_down=0.0
      endif
! CONSTRUCT TOTAL DENSITY AND CONTRIBUTION TO ex
      rho=rho_up+rho_down
      eps_x=(exuplsda*rho_up+exdnlsda*rho_down)/rho
!
! --CORRELATION--
!
      if(rho.gt.rhomin) then
         fk=conf*rho**thrd
         rs=alpha/fk
         zet=(rho_up-rho_down)/rho
         call corlsd(rs,zet,eps_c,v_c_up,v_c_down)
      else
         eps_c=0.0
         v_c_up=0.0
         v_c_down=0.0
      endif
      return
      end subroutine xc_lsda
! =====================================================================
      subroutine corlsd(rs,zet,ec,vcup,vcdn)
! =====================================================================
! *********************************************************************
! REFERENCES:
! [a] J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
! *********************************************************************
! INPUT:
!      real rs   : SEITZ RADIUS=(3/4pi rho)^(1/3)
!      real zet  : RELATIVE SPIN POLARIZATION = (rhoup-rhodn)/rho
! OUTPUT: 
!      real ec   : LSD CORRELATION ENERGY FROM [a]
!      real vcup : LSD UP CORRELATION POTENTIAL
!      real vcdn : LSD DN CORRELATION POTENTIAL
! LOCAL:
!      real thrd,thrdm,thrd2,thrd4,sixthm
! NUMBERS FOR USE IN LSD ENERGY SPIN-INTERPOLATION FORMULA, [a](9)
!      real gam     : 2^(4/3)-2
!      real fzz     : f''(0)= 8/(9*GAM)
! CORRELATION ENERGY PART
!      real rtrs    : sqrt(rs)
!      real z4      : zet^4
!      real eu      : UNPOLARIZED LSD CORRELATION ENERGY
!      real eurs    : d(eu)/d(rs)
!      real ep      : FULLY POLARIZED LSD CORRELATION ENERGY
!      real eprs    : d(ep)/d(rs)
!      real alfc,alfm    : +/- SPIN STIFFNESS, [a](3).
!      real alfrsm  : -d(alpha)/d(rs)
!      real f       : spin-scaling factor from [a](9).
!      real tmp1,tmp2,tmp3,tmp4,tmp5,tmp6
! CORRELATION POTENTIAL PART
!      real ecrs    : d(ec)/d(rs) [a](A2)
!      real eczet   : d(ec)/d(zeta) [a](A3)
!      real fz      : d(f)/d(zeta) [a](A4)
!      real comm
! **********************************************************************    
      implicit none
      real*8 rs,zet,ec,vcup,vcdn,thrd,thrdm,thrd2,thrd4,sixthm,&
          gam,fzz,rtrs,z4,eu,eurs,ep,eprs,alfc,alfm,alfrsm,f,&
          tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,ecrs,eczet,fz,comm
!
      parameter(thrd=1./3.,thrdm=-thrd,thrd2=2.*thrd,thrd4=4.*thrd)
      parameter(sixthm=thrdm/2.0)
! ----------------------------------------------------------------------
      gam = 2.**thrd4-2.0
      fzz=8./(9.*gam)
! LSD ENERGY CONTRIBUTIONS, USING [a](10) and Table I[a].
! CONSTRUCT EC, USING [a](8)
      rtrs=sqrt(rs)
      tmp1= 0.03109070
      tmp2=0.213700
      tmp3=7.59570
      tmp4=3.58760
      tmp5=1.63820
      tmp6=0.492940
      CALL gcor(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,rtrs,eu,eurs)
      tmp1= 0.015545350
      tmp2=0.205480
      tmp3=14.11890
      tmp4=6.19770
      tmp5=3.36620
      tmp6=0.625170
      CALL gcor(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,rtrs,ep,eprs)
      tmp1= 0.01688690
      tmp2=0.111250
      tmp3=10.3570
      tmp4=3.62310
      tmp5=0.880260
      tmp6=0.496710
      CALL gcor(tmp1,tmp2,tmp3,tmp4,tmp5,tmp6,rtrs,alfm,alfrsm)
      alfc = -alfm
      z4 = zet**4
      f=((1.0+zet)**thrd4+(1.0-zet)**thrd4-2.0)/gam
      ec = eu*(1.0-f*z4)+ep*f*z4-alfm*f*(1.0-z4)/fzz
!
!
! LSD POTENTIAL FROM [a](A1)
!
      ecrs = eurs*(1.0-f*z4)+eprs*f*z4-alfrsm*f*(1.0-z4)/fzz
      fz = thrd4*((1.0+zet)**thrd-(1.0-zet)**thrd)/gam
      eczet = 4.0*(zet**3)*f*(ep-eu+alfm/fzz)+fz*(z4*ep-z4*eu&
             -(1.0-z4)*alfm/fzz)
      comm = ec -rs*ecrs/3.0-zet*eczet
      vcup = comm + eczet
      vcdn = comm - eczet
!
      return
      end subroutine corlsd
! =====================================================================
      subroutine gcor(a,a1,b1,b2,b3,b4,rtrs,gg,ggrs)
! =====================================================================
! *********************************************************************
! SLIMMED DOWN VERSION OF gcor USED IN PW91 ROUTINES, TO INTERPOLATE
! LSD CORRELATION ENERGY, AS GIVEN BY (10) OF
! J. P. Perdew and Y. Wang, Phys. Rev. B {\bf 45}, 13244 (1992).
! *********************************************************************
! INPUT :
!      real a,a1,b1,b2,b3,b4,rtrs
! OUTPUT:
!      real gg,ggrs
! LOCAL :
!      real q0,q1,q2,q3
! *********************************************************************
      implicit none
      real*8 a,a1,b1,b2,b3,b4,rtrs,gg,ggrs,q0,q1,q2,q3
! ---------------------------------------------------------------------
      q0 = -2.0*a*(1.0+a1*rtrs*rtrs)
      q1 = 2.0*a*rtrs*(b1+rtrs*(b2+rtrs*(b3+b4*rtrs)))
      q2 = log(1.0+1.0/q1)
      gg = q0*q2
!
      q3 = a*(b1/rtrs+2.0*b2+rtrs*(3.0*b3+4.0*b4*rtrs))
      ggrs = -2.0*a*a1*q2-q0*q3/(q1*(1.0+q1))
!
      return
      end subroutine gcor




END MODULE XC

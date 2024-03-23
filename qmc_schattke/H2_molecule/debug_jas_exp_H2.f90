!------------------------------------
subroutine JAS_EXP_H2(vj, vjd, v, vd)
!------------------------------------
  use m_highlevel, only: DP, Nelectrons, NES1, EMACH
  use m_midlevel, only: IE, RE, RNEU, RK
  use m_jastrow, only: gam, QJC, DIST, DISTNEU, LAPJASOLD, LAPJAS
  use m_jastrow, only: DKX, CJAS, BETA1, BETA2, GRJAS, GRJASOLD
  implicit none
! Updates the distances from the active electron to all others
! and determines Jastrow exponent and derivatives for atoms:
! u=F/2/(1+r/F)*(delta_{s,-s'}+1/2*delta_{s,s'})
! In addition, the wave function of Kolos and Roothaan
! (James and Coolidge)
! (Rev.Mod..Phys.32,219(1960)) is programmed as part of the Jastrow
! factor. It has to be combined with product ansatz of orbital wave.
  real(dp), intent(out), dimension(Nelectrons):: vj,vjd,v,vd
  integer :: k,kie,ii
  real(dp) :: as,jasu,jasdif,u2d,woo,won,u3,u4,u2do
  real(dp) :: jasup,jasdifp
  real(dp) :: g1mrn,g1mro,g1lmn,g1lmo,g2ln,g2lo,g2mn,g2mo,lapln, laplo,lapmn,lapmo
  real(dp), dimension(3) :: u1d,u1do,jc1dn,jc1do,u1dp,u1dpo
  real(dp), dimension(3) :: rran,rrbn,rrao,rrbo
  real(dp) :: rhn,rho,jasjcn,jasjco,jasujcn,jasujco,jc2dn,jc2do
  real(dp), dimension(2) :: ran,rbn,rao,rbo,lan,lao,mun,muo
  real(dp), dimension(12) :: a,xn,xo,zn,zo
  real(dp), dimension(3,12) :: yn,yo
        
  a(:) = 0.d0
  a(1) = 1.d0
  ! Without series expansion comment 12 lines of constants from KR
  ! Constants from Kolos and Roothaan
  a(1) = +2.192089d0
  a(2) = +1.098975d0
  a(3) = -0.377500d0
  a(4) = -0.139338d0
  a(5) = +0.859247d0
  a(6) = -0.058316d0
  a(7) = +0.078257d0
  a(8) = +0.150633d0
  a(9) = -0.052156d0
  a(10) = -0.126629d0
  a(11) = +0.132561d0
  a(12) = +0.248411d0
  ! For the 5 term series of JC comment the last 12 lines and
  ! uncomment the following 5 constants from James and Coolidge
  ! a(1)=2.23779_dp
  ! a(2)=0.80483_dp
  ! a(3)=-0.27997_dp
  ! a(4)=-0.60985_dp
  ! a(5)=0.19917_dp
  jasu = 0.d0
  jasup = 0.d0
  jasdif = 0.d0
  jasdifp = 0.d0
  u1d(1:3) = 0.d0
  u1dp(1:3) = 0.d0
  u1do(1:3) = 0.d0
  u1dpo(1:3) = 0.d0
  u2d = 0.d0
  u2do = 0.d0
  u3 = 0.d0
  u4 = 0.d0

  ielekat: do k=1,Nelectrons
    !
    if( k .eq. IE ) then
      cycle ielekat
    endif
    !
    as = 0.5d0     ! for equal spins
    !
    if( ( (IE.le.NES1).and.(k.gt.NES1)) .or. &
     &       ((IE.gt.NES1).and.(k.le.NES1))) as=1.0d0
  
    DIST(1:3,IE,k) = RE(1:3,IE) - RE(1:3,k)
    DISTNEU(1:3,IE,k) = RNEU(1:3,IE) - RE(1:3,k)
    woo = max(2.d0/3.d0*EMACH, sqrt(sum(DIST(1:3,IE,k)**2)))
    won = max(2.d0/3.d0*EMACH, sqrt(sum(DISTNEU(1:3,IE,k)**2)))
    DIST(4,IE,k) = woo
    DISTNEU(4,IE,k) = won
    !
    jasdif = jasdif + CJAS*(as/2.0d0/(1.0d0 + won/CJAS) - &
          &  as/2.0d0/(1.0d0 + woo/CJAS)) + &
          &  BETA1*(exp(-BETA2*won**2) - exp(-BETA2*woo**2))
    !
    jasdifp = jasdifp + &
          &   sum( (RNEU(1:3,IE) + RE(1:3,k) - RK(1:3,2) )**2 ) - &
          &   sum( (RE(1:3,IE) + RE(1:3,k) - RK(1:3,2) )**2 )
    !
    jasu = jasu + as*CJAS/2.0d0/(1.0d0+won/CJAS) + BETA1*exp(-BETA2*won**2)
    !
    jasup = jasup + sum( (RNEU(1:3,IE) + RE(1:3,k) - RK(1:3,2))**2 )
    !
    u1d(1:3) = u1d(1:3)-as*DISTNEU(1:3,IE,k) / won / &
            & (1.0d0+won/CJAS)**2/2.0d0 - &
            &  2.0d0*DISTNEU(1:3,IE,k)*BETA1*BETA2*exp(-BETA2*won**2)
    !     
    u1dp(1:3) = RNEU(1:3,IE) + RE(1:3,k) - RK(1:3,2)
    !     
    u1do(1:3) = u1do(1:3)-as*DIST(1:3,IE,k) / woo / &
            & (1.0d0+woo/CJAS)**2/2.0d0 - &
            &  2.0d0*DIST(1:3,IE,k)*BETA1*BETA2*exp(-BETA2*woo**2)
    !
    u1dpo(1:3) = RE(1:3,IE) + RE(1:3,k) - RK(1:3,2)
    !
    u2d = u2d - as/won/(1.0d0 + won/CJAS)**3 + &
            & BETA1*BETA2*exp(-BETA2*won**2)*(4.0d0*BETA2*won**2 - 6.0d0)
  
    u2do = u2do - as/woo/(1.0d0+woo/CJAS)**3 + &
            & BETA1*BETA2*exp(-BETA2*woo**2)*(4.0d0*BETA2*woo**2 - 6.0d0)
    !
    u3 = u3 + 1.0d0/won
    u4 = u4 + 1.0d0/won - 1.0d0/woo
  enddo ielekat
  !
  jasu = jasu + jasup*GAM/2.0d0
  jasdif = jasdif + jasdifp*GAM/2.0d0
  u1d = u1d + u1dp*GAM
  u1do = u1do + u1dpo*GAM
  u2d = u2d + 3.d0*GAM
  u2do = u2do + 3.d0*GAM
  ! For an additional Jastrow factor: psi_J=J*sum_{ii=1}^5 a(ii)*x(ii):
  ! The basis functions x(ii) are formulated in terms of
  ! lambda1,lambda2,mu1,mu2,rho from James and Coolidge, which are abbreviated
  ! by the first two letters; the letters a and b denote the two protons and
  ! the letters n and o are appended to refer to new and old;
  ! the suffixes 1 and 2 are displayed by the
  ! index IE and kie of the actual electron which might have been moved
  ! and the other second electron which keeps its old position, resp..
  kie = 2
  if( IE == 2 ) kie=1
  rhn = 2.d0*DISTNEU(4,IE,kie)/DKX
  rho = 2.d0*DIST(4,IE,kie)/DKX
  ran(IE) = max(EMACH, sqrt(sum((RNEU(1:3,IE)-RK(1:3,1))**2)))
  rao(IE) = max(EMACH, sqrt(sum((RE(1:3,IE)-RK(1:3,1))**2)))
  rbn(IE) = max(EMACH, sqrt(sum((RNEU(1:3,IE)-RK(1:3,2))**2)))
  rbo(IE) = max(EMACH, sqrt(sum((RE(1:3,IE)-RK(1:3,2))**2)))
  lan(IE) = (ran(IE) + rbn(IE))/DKX
  lao(IE) = (rao(IE) + rbo(IE))/DKX
  mun(IE) = (ran(IE) - rbn(IE))/DKX
  muo(IE) = (rao(IE) - rbo(IE))/DKX
  rao(kie) = max(EMACH, sqrt(sum((RE(1:3,kie) - RK(1:3,1))**2)))
  rbo(kie) = max(EMACH, sqrt(sum((RE(1:3,kie) - RK(1:3,2))**2)))
  lao(kie) = (rao(kie) + rbo(kie))/DKX
  muo(kie) = (rao(kie) - rbo(kie))/DKX
  ! Accepted step: new coordinates for IE and old for kie
  xn(1) = 2.d0
  xn(2) = mun(IE)**2+muo(kie)**2
  xn(3) = 2.d0*mun(IE)*muo(kie)
  xn(4) = lan(IE)+lao(kie)
  xn(5) = 2.d0*rhn
  xn(6) = (lan(IE)+lao(kie))*mun(IE)*muo(kie)
  xn(7) = lan(IE)*muo(kie)**2+lao(kie)*mun(IE)**2
  xn(8) = lao(kie)**2+lan(IE)**2
  xn(9) = 2.d0*rhn**2
  xn(10) = 2.0d0*lan(IE)*lao(kie)
  xn(11) = 2.0d0*mun(IE)**2*muo(kie)**2
  xn(12) = (muo(kie)**2 + mun(IE)**2)*rhn
  ! Not accepted step: old coordinates for both IE and kie
  xo(1) = 2.0d0
  xo(2) = muo(IE)**2 + muo(kie)**2
  xo(3) = 2.0d0*muo(IE)*muo(kie)
  xo(4) = lao(IE) + lao(kie)
  xo(5) = 2.0d0*rho
  xo(6) = (lao(IE) + lao(kie))*muo(IE)*muo(kie)
  xo(7) = lao(IE)*muo(kie)**2 + lao(kie)*muo(IE)**2
  xo(8) = lao(kie)**2 + lao(IE)**2
  xo(9) = 2.0d0*rho**2
  xo(10) = 2.0d0*lao(IE)*lao(kie)
  xo(11) = 2.0d0*muo(IE)**2*muo(kie)**2
  xo(12) = (muo(kie)**2 + muo(IE)**2)*rho
  ! The 1st derivative (new and old) is denoted by yn(3,12) and yo(3,12)
  rran(1:3)=(RNEU(1:3,IE)-RK(1:3,1))/DKX/ran(IE)
  rrbn(1:3)=(RNEU(1:3,IE)-RK(1:3,2))/DKX/rbn(IE)
  rrao(1:3)=(RE(1:3,IE)-RK(1:3,1))/DKX/rao(IE)
  rrbo(1:3)=(RE(1:3,IE)-RK(1:3,2))/DKX/rbo(IE)
  ! Accepted step
  yn(1:3,1) = 0.0d0
  yn(1:3,2) = 2.0d0*mun(IE)*(rran-rrbn)
  yn(1:3,3) = 2.0d0*muo(kie)*(rran-rrbn)
  yn(1:3,4) = rran + rrbn
  yn(1:3,5) = 4.0d0*DISTNEU(1:3,IE,kie)/DISTNEU(4,IE,kie)/DKX
  yn(1:3,6) = (rran + rrbn)*mun(IE)*muo(kie) + (lan(IE)+lao(kie))*(rran-rrbn)*muo(kie)
  yn(1:3,7) = (rran + rrbn)*muo(kie)**2 + 2.0d0*lao(kie)*mun(IE)*(rran-rrbn)
  yn(1:3,8) = 2.0d0*lan(IE)*(rran+rrbn)
  yn(1:3,9) = 16.0d0*DISTNEU(1:3,IE,kie)/DKX**2
  yn(1:3,10) = 2.0d0*(rran + rrbn)*lao(kie)
  yn(1:3,11) = 4.0d0*(rran - rrbn)*mun(IE)*muo(kie)**2
  !
  yn(1:3,12) = 2.0d0*mun(IE)*(rran - rrbn)*rhn + &
            & (muo(kie)**2+mun(IE)**2)*2.0d0*DISTNEU(1:3,IE,kie)/ &
            &  DISTNEU(4,IE,kie)/DKX
  !
  ! Not accepted step
  yo(1:3,1) = 0.0d0
  yo(1:3,2) = 2.0d0*muo(IE)*(rrao-rrbo)
  yo(1:3,3) = 2.0d0*muo(kie)*(rrao-rrbo)
  yo(1:3,4) = rrao + rrbo
  yo(1:3,5) = 4.0d0*DIST(1:3,IE,kie)/DIST(4,IE,kie)/DKX     
  yo(1:3,6) = (rrao + rrbo)*muo(IE)*muo(kie) + &
           &  (lao(IE) + lao(kie))*(rrao - rrbo)*muo(kie)
  yo(1:3,7) = (rrao+rrbo)*muo(kie)**2 + 2.0d0*lao(kie)*muo(IE)*(rrao-rrbo)
  yo(1:3,8) = 2.0d0*lao(IE)*(rrao+rrbo)
  yo(1:3,9) = 16.0d0*DIST(1:3,IE,kie)/DKX**2
  yo(1:3,10) = 2.0d0*(rrao + rrbo)*lao(kie)
  yo(1:3,11) = 4.0d0*(rrao - rrbo)*muo(IE)*muo(kie)**2
  !
  yo(1:3,12) = 2.0d0*muo(IE)*(rrao-rrbo)*rho+ &
           & (muo(kie)**2+muo(IE)**2)*2.0d0*DIST(1:3,IE,kie)/DIST(4,IE,kie)/DKX
  ! The 2nd derivative (new and old) is denoted by zn(5) and zo(5)
  g1mrn = dot_product(DISTNEU(1:3,IE,kie), (rran(1:3)-rrbn(1:3))) * 4.0d0/DKX**2/rhn
  !
  g1mro = dot_product(DIST(1:3,IE,kie), (rrao(1:3)-rrbo(1:3))) * 4.0d0/DKX**2/rho
  !
  g1lmn = dot_product(rran,rrbn)*DKX**2
  g1lmo = dot_product(rrao,rrbo)*DKX**2
  g2ln = 2.0d0*(1.0d0 + g1lmn)/DKX**2
  g2lo = 2.0d0*(1.0d0 + g1lmo)/DKX**2
  g2mn = 2.0d0*(1.0d0 - g1lmn)/DKX**2
  g2mo = 2.0d0*(1.0d0 - g1lmo)/DKX**2
  lapln = 2.0d0*(1.0d0/ran(IE) + 1.0d0/rbn(IE))/DKX
  laplo = 2.0d0*(1.0d0/rao(IE) + 1.0d0/rbo(IE))/DKX
  lapmn = 2.0d0*(1.0d0/ran(IE) - 1.0d0/rbn(IE))/DKX
  lapmo = 2.0d0*(1.0d0/rao(IE) - 1.0d0/rbo(IE))/DKX
  !
  ! Accepted step
  zn(1) = 0.0d0
  zn(2) = 2.0d0*(g2mn+mun(IE)*lapmn)
  zn(3) = 2.0d0*muo(kie)*lapmn
  zn(4) = lapln
  zn(5) = 16.0d0/rhn/DKX**2
  zn(6) = muo(kie)*(lapln*mun(IE)+lapmn*(lan(IE)+lao(kie)))
  zn(7) = lapln*muo(kie)**2+2.0d0*lao(kie)*(g2mn+mun(IE)*lapmn)
  zn(8) = 2.0d0*(g2ln+lan(IE)*lapln)
  zn(9) = 48.0d0/DKX**2
  zn(10) = 2.0d0*lao(kie)*lapln
  zn(11) = 4.0d0*muo(kie)**2*(g2mn+mun(IE)*lapmn)
  zn(12) = 2.0d0*rhn*(g2mn+mun(IE)*lapmn) + 4.0d0*mun(IE)*g1mrn + &
     &     8.0d0*(mun(IE)**2+muo(kie)**2)/rhn/DKX**2

  ! Not accepted step
  zo(1) = 0.0d0
  zo(2) = 2.0d0*(g2mo+muo(IE)*lapmo)
  zo(3) = 2.0d0*muo(kie)*lapmo
  zo(4) = laplo
  zo(5) = 16.0d0/rho/DKX**2
  zo(6) = muo(kie)*(laplo*muo(IE)+lapmo*(lao(IE)+lao(kie)))
  zo(7) = laplo*muo(kie)**2+2.0d0*lao(kie)*(g2mo+muo(IE)*lapmo)
  zo(8) = 2.0d0*(g2lo+lao(IE)*laplo)
  zo(9) = 48.0d0/DKX**2
  zo(10) = 2.0d0*lao(kie)*laplo
  zo(11) = 4.0d0*muo(kie)**2*(g2mo+muo(IE)*lapmo)
  zo(12) = 2.0d0*rho*(g2mo+muo(IE)*lapmo)+4.0d0*muo(IE)*g1mro + &
     &     8.0d0*(muo(IE)**2+muo(kie)**2)/rho/DKX**2
  jasjcn = 0.0d0
  jasjco = 0.0d0
  jc1dn(1:3) = 0.0d0
  jc1do(1:3) = 0.0d0
  jc2dn = 0.0d0
  jc2do = 0.0d0
  !
  ljc: do ii = 1,12
    jasjcn = jasjcn + a(ii)*xn(ii)
    jasjco = jasjco + a(ii)*xo(ii)
    jc1dn(1:3) = jc1dn(1:3)+a(ii)*yn(1:3,ii)
    jc1do(1:3) = jc1do(1:3)+a(ii)*yo(1:3,ii)
    jc2dn = jc2dn+a(ii)*zn(ii)
    jc2do = jc2do+a(ii)*zo(ii)
  enddo ljc
  
  ! if((jasjcn .le. 0.0d0) .or. (jasjco .le. 0.0d0)) then
  !   write(*,*)'non-positive argument of log'
  !   stop
  ! endif

  ! jasjc or jasjco might be negative, but only the square is needed
  ! for use in probability measure. So, take modulus.
  jasujcn=-dlog(dabs(jasjcn))
  jasujco=-dlog(dabs(jasjco))
  ! Instead calculate directly the acceptance ratio qjc
  QJC = jasjcn/jasjco 
  ! For derivatives wave function is just a factor, no exponentiation
  jc1dn(1:3) = jc1dn(1:3)/jasjcn
  jc1do(1:3) = jc1do(1:3)/jasjco
  jc2dn = jc2dn/jasjcn
  jc2do = jc2do/jasjco
  vjd(IE) = jasdif + jasujcn - jasujco
  vj(IE) = jasu + jasujcn
  GRJAS(1:3,IE) = - u1d(1:3) + jc1dn(1:3)
  
  LAPJAS(IE) = -u2d + sum(u1d(1:3)**2) + jc2dn - 2.0d0*dot_product(jc1dn(1:3),u1d(1:3))
  GRJASOLD(1:3,IE) = - u1do(1:3) + jc1do(1:3)

  LAPJASOLD(IE) = - u2do + sum (u1do(1:3)**2) + jc2do - 2.0d0*dot_product(jc1do(1:3),u1do(1:3))
  v(IE) = u3
  vd(IE) = u4
end subroutine 


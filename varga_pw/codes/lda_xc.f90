SUBROUTINE lda_xc(rr,vxc,exc,ld)   
!   
!   LDA exchange correlation potential Phys.Rev. B 23,5048(1981)                            
!   (atomic units used)
!
! Perdew-Zunger LDA
IMPLICIT NONE 
  real*8, parameter :: x0 = 0.0d0, x6 = 6.0d0, x1 = 1.0d0, x7 = 7.0d0, &
& x2 = 2.0d0, x3 = 3.0d0, x4 = 4.0d0, x9 = 9.0d0, eps0= 1.0d-1, &
& pi=3.14159265358979323d0,au=0.0311d0,bu=-0.048d0,cu=0.002d0,du=-0.0116d0,aph=-0.1423d0, &
& bt1=1.0529d0,bt2=0.3334d0,aa=-x3/(x4*pi),bu1=bu-au/x3,cu1=x2*cu/x3, &
& du1= (x2*du-cu)/x3,p13=x1/x3,pi4=x4*pi,c76=x7/x6,c43=x4/x3
  INTEGER  :: ld
  real*8  :: exc,rr(ld), vxc(ld)
  real*8  :: eps,bpr,ec,ex,r,rrs,rs,rsl,rss,uc,ux,ax
  INTEGER  :: i,ip,j,l
      eps = 1.d-13
      exc = 0.0
      ax = aa * (x9 * pi / x4)**p13 
      do  i=1,ld
        r = rr(i)
        if (r .lt. -eps0) then                  
          write(6,*)'negative electron density'
          stop
          return
        else if (r.le.eps) then
          ux = x0
          uc = x0
          ec = x0
          ex = x0
        else
          rs = (x3/(pi4*r))**p13
          rrs = sqrt(rs)
          rss = rs*rs
          ex = ax/rs
          ux = c43*ex
          if (rs.ge.x1) then
            bpr = bt1*rrs + x1 + bt2*rs
            ec = aph/bpr
            uc = (c76*bt1*rrs+x1+c43*bt2*rs)*ec/bpr
          else
            rsl = log(rs)
            ec = du*rs + bu + au*rsl + cu*rs*rsl
            uc = du1*rs + bu1 + au*rsl + cu1*rs*rsl
          endif
        endif
        vxc(i) = ux + uc
        exc = exc + r * (ex + ec)
      ENDDO 
  return
end SUBROUTINE lda_xc

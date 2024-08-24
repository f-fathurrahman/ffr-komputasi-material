implicit none
integer                      :: i,j
double precision,parameter   :: del=1.d-14,small=1.d-28
double precision,parameter   :: h2m=0.5d0,dx=0.1d0,epsilon=0.05d0
integer,parameter            :: nmax=10000,n_e=3000
complex*16                   :: a(0:nmax),b(0:nmax),f,c,d
double precision             :: alpha,beta,x
complex*16                   :: energy,g0,g
complex*16,parameter         :: zi=(0.d0,1.d0)

alpha=2.d0*h2m/dx**2; beta=-h2m/dx**2
b(0)=(0.d0,0.d0); a(1)=(1.d0,0.d0)
do i=1,n_e
  energy=-3.d0*beta+alpha+6.d0*beta/dfloat(n_e)*i+zi*epsilon
  do j=1,nmax
    b(j)=energy-alpha
    if(j.ne.1) a(j)=-beta**2
  end do
  if(b(0)==0.d0) then
    f=small
  else
    f=b(0)
  endif
  c=f
  d=0.d0
  do j=1,nmax
    d=b(j)+a(j)*d
    if(d==0.d0) d=small
    c=b(j)+a(j)/c
    if(c==0.d0) c=small
    d=1.d0/d
    f=f*c*d
    if(abs(1.d0-c*d).lt.del) exit
  end do
! surface green's function
  g0=f
! green's function
  g=1.d0/(energy-alpha-2.d0*beta**2*g0)
end do
end
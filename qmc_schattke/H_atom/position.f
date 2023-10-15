C-----------------------------------------------------------------------
      module position
       implicit none
       public :: DENSITY1D
       integer,parameter,public          :: NRHO=1000
       real(kind=8),parameter,public    :: EMACH=1.0e-8_dp
       integer,public                    :: IMC
       real(kind=8),public              :: LENGTH,DRHO
       real(kind=8),dimension(3),public :: RE,RNEU
       real(8),public,dimension(NRHO)   :: RHORAD,AVRHORAD
      contains
C-----------------------------------------------------------------------
      subroutine DENSITY1D
C   Radial density RHORAD is a function of distance s from nucleus and
C   normalized to 1=sum_s 4*PI*s**2 ds RHORAD(s). It is discretized
C   in units of DRHO with NRHO sections. Values below DRHO are added
C   to first unit and those above NRHO*DRHO added to last unit.
       integer    :: j
       real(8)   :: s,h
       RHORAD=0.d0
       h=4.d0*PI*DRHO
       s=max(sqrt(sum (RNEU(1:3)**2)),DRHO+EMACH)
       s=min(s,NRHO*DRHO)
       j=int(s/DRHO)
       RHORAD(j)=1/h/s**2
      end subroutine DENSITY1D
C-----------------------------------------------------------------------
      end module position


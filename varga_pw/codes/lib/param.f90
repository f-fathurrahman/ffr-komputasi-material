!
! $Id: param.f90,v 1.3 1999/03/22 11:52:48 noe Exp $
! Last pre-CVS change:  GUE   3 Sep 96    7:52 pm
!
module param

   implicit none

      INTEGER, parameter :: r4_kind=KIND(0.0)
      INTEGER, parameter :: r8_kind=KIND(0.d0)

      INTEGER, parameter :: c_kind=KIND((0.0,0.0))
      INTEGER, parameter :: c4_kind=c_kind
      INTEGER, parameter :: c8_kind=KIND((0.d0,0.d0))

      integer nsx,nax,nx,ngwx,ngx
      integer ngwix,nx_init,nr1x,nr2x,nr3x,nschltz
      integer nx_basis,max_basis_n,nlmax_init
      integer nnrx,nkptx,nlmax,mmaxx,n_fft_store
      INTEGER nstatx !=nkptx*nx


end module


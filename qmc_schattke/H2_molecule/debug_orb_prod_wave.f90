!---------------------------
subroutine orb_prod_wav(r, psi)
!---------------------------
  use m_highlevel, only: DP, Nelectrons, Natoms, NelectronsPerSpin, EMACH
  use m_midlevel, only: RK
  use m_orbital, only: WAVEC, NORB, ALPHA
  implicit none
! Calculate orbital phi=exp(-alpha*r) and wavefunction psi as 
! product of phi.Only one orbital.
! Up spin runs from IE=1,NelectronsPerSpin, down spin from IE=NelectronsPerSpin+1,Nelectrons
! Should be tabulated later when we use determinants!
  real(dp),intent(in), dimension(3) :: r
  real(dp),intent(out), dimension(NORB,Natoms) :: psi
  integer :: k
  real(dp) :: s
  real(dp),dimension(3) :: rad
  real(dp),dimension(NORB,Natoms) :: phi
  do k=1,Natoms
    rad(1:3)=r(1:3)-RK(1:3,k)
    s = sqrt( sum(rad(1:3)**2) ) 
    ! Single orbitals look like (note: orbitals differ only by shifting)
    ! The extension by WAVEC .NEQ. 0.0 seems useless as covered by change 
    ! of ALPHA, for small values at least
    phi(1:NORB,k)=(1.0_dp+WAVEC*s)*dexp(-ALPHA*s)
  enddo
  ! Associate single particle wavefunction with electron, here 1 to 1
  psi(1:NORB,1) = phi(1:NORB,1)*phi(1:NORB,2) ! product ansatz
  psi(1:NORB,2) = phi(1:NORB,1)*phi(1:NORB,2) ! product ansatz
end subroutine



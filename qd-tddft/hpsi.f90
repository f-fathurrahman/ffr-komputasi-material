!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE HPSI
! ===============
!
! INPUT:
!   v [real(8), dimension(n, n)] : Kohn-Sham potential.
!   f [real(8), dimension(n, n)] : wave function on which the KS Hamiltonian is applied.
! --------
! OUTPUT:
!   hf [real(8), dimension(n, n)] : the resulting wavefunction: |hf> = H|f>
!
! This subroutine operates the Kohn-Sham Hamiltonian (the kinetic term, plus
! the local potential v that should contain the sum of all the terms -- external,
! Hartree and exchange and correlatioin) on a given wavefunction f, and puts
! the result of hf.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine hpsi(v, f, hf)
  use mesh

  implicit none

  real(8), intent(in)  :: v(N, N)
  real(8), intent(in)  :: f(N, N)    ! function whose hamiltonian is going to be calculated
  real(8), intent(out) :: hf(N, N) ! the results

  hf = 0.0
  call laplacian(f, hf)
  hf(:, :) = -0.5_8*hf(:, :) + v(:, :)*f(:, :)

end subroutine hpsi




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE ZHPSI
! ===============
!
! INPUT:
!   v [real(8), dimension(n, n)] : Kohn-Sham potential.
!   f [complex(8), dimension(n, n)] : wave function on which the KS Hamiltonian is applied.
! --------
! OUTPUT:
!   hf [complex(8), dimension(n, n)] : the resulting wavefunction: |hf> = H|f>
!
! This subroutine operates the Kohn-Sham Hamiltonian (the kinetic term, plus
! the local potential v that should contain the sum of all the terms -- external,
! Hartree and exchange and correlatioin) on a given wavefunction f, and puts
! the result of hf.
!
! NOTES:
! ------
! It is a repetition of previous subroutine hpsi, but that operates on
! complex wavefunctions.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine zhpsi(v, f, hf)
  use mesh

  implicit none

  real(8), intent(in)  :: v(N, N)
  complex(8), intent(in)  :: f(N, N)    ! function whose hamiltonian is going to be calculated
  complex(8), intent(inout) :: hf(N, N) ! the results

  hf = (0.0, 0.0)
  call zlaplacian(f, hf)
  hf(:, :) = -0.5_8*hf(:, :) + v(:, :)*f(:, :)

end subroutine zhpsi

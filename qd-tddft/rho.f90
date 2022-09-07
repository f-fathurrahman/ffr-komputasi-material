!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Calculates the density given the Kohn-Sham orbitals.
subroutine build_rho(w, rho)
  use states
  use mesh

  implicit none
  real(8), intent(in)  :: w(N, N, N_wf)
  real(8), intent(out) :: rho(N, N)

  integer :: ix, iy, i

  rho = 0.0_8
  do i = 1, N_occ
     do ix = 1, N
        do iy = 1, N
           rho(ix, iy) = rho(ix, iy) + 2.0_8*w(ix, iy, i)**2
        enddo
     enddo
  enddo

end subroutine build_rho

subroutine zbuild_rho(w, rho)
  use states
  use mesh

  implicit none
  complex(8), intent(in) :: w(N, N, N_wf)
  real(8), intent(out) :: rho(N, N)

  integer :: ix, iy, i

  rho = 0.0_8
  do i = 1, N_occ
     do ix = 1, N
        do iy = 1, N
           rho(ix, iy) = rho(ix, iy) + 2.0_8*conjg(w(ix, iy, i))*w(ix, iy, i)
        enddo
     enddo
  enddo

end subroutine zbuild_rho

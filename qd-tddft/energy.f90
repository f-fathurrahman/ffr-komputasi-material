real(8) function energy(eigenval, rho, vh, vc, vx, ec, ex) result(e)
  use states
  use mesh

  implicit none
  real(8), intent(in) :: eigenval(N_wf)
  real(8), intent(in) :: ec, ex, rho(N, N), vh(N, N), vc(N, N), vx(N, N)

  e = 2.0*sum(eigenval(1:N_occ)) + ec + ex
  e = e - sum(rho(:, :)*(vx(:, :) + vc(:, :) + vh(:, :)/2.0))*delta**2

end function energy

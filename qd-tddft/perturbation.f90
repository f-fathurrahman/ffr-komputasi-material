!/*!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! SUBROUTINE PERTURBATION
!
!
!*/!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine perturbation(N_wf, wfs)
  use mesh
  implicit none
  integer, intent(in) :: N_wf
  complex(8), intent(inout) :: wfs(N, N, N_wf)

  integer :: j, ix, iy

  do j = 1, N_wf
     do ix = 1, n
        do iy = 1, n
           wfs(ix, iy, j) = exp((0._8,1._8)*0.01_8*x(ix,iy))*wfs(ix, iy, j)
        enddo
     enddo
  enddo

end subroutine perturbation

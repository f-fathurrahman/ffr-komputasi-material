!! This program is free software; you can redistribute it and/or modify
!! it under the terms of the GNU General Public License as published by
!! the Free Software Foundation; either version 2, or (at your option)
!! any later version.
!!
!! This program is distributed in the hope that it will be useful,
!! but WITHOUT ANY WARRANTY; without even the implied warranty of
!! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!! GNU General Public License for more details.
!!
!! You should have received a copy of the GNU General Public License
!! along with this program; if not, write to the Free Software
!! Foundation, Inc., 59 Temple Place - Suite 330, Boston, MA


subroutine init_wfs()
  use states
  use mesh

  implicit none
  integer :: i, ix, iy
  do i = 1, N_wf
    do ix = 1, N
       do iy = 1, N
          call random_number(wfs(ix, iy, i))
       enddo
    enddo
  enddo
end subroutine init_wfs

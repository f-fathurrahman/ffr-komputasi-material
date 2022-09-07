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

subroutine write_wfs(filename)
  use states
  use mesh
  implicit none

  character(len=*) :: filename

  integer :: i

  open(unit=11, file=trim(filename), form='unformatted')
  write(unit=11) N_occ, N_empty
  do i = 1, n_wf
     write(unit=11) wfs(:, :, i)
  enddo

  close(unit=11)

end subroutine write_wfs

! Declerations for EM force routine
Module Defs
save

real, parameter :: PI = 3.14159265358979


type elements
  integer :: atno
  character(2) :: name
  real :: a0, a1, a2, rm, om, m, rc, alf, ec, c
  integer :: z
end type


type(elements) :: sp(10)
type(elements) :: el(56)

real :: rr, rrdir(3), omave, zave, mave, x, y, z, box_x, box_y, box_z
real :: cutoff, xr, yr, zr
real :: kf, ef, vf, vd, vddir(3), rs, ksi
integer :: i1, i2, force_j, force_i



end module Defs


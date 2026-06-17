! =============================================================================
! MODULE: commondefs
! =============================================================================
module commondefs
  save
  type atom
    real x0, y0, z0, x, y, z, fx, fy, fz, elfx, elfy, elfz, vx, vy, vz
    real ax, ay, az, bx, by, bz, cx, cy, cz, rox, roy, roz, roxl, royl, roz1, &
         xp, yp, zp, xb, yb, zb
    integer kind
  end type atom

  type ionpos
    real x, y, z
    integer kind
  end type ionpos

  type(atom), allocatable :: atoms(:)
  type(ionpos), allocatable :: ions(:)

  integer n, ntype, nvac, nwall
  real boxsize_x, boxsize_y, boxsize_z, rcut, a
  integer ncellx, ncelly, ncellz
  real cell_factorx, cell_factory, cell_factorz, bxr, byr, bzr
  integer head_chain(0:1000), linked_list(40010), n_v_list(40000)
  real atno(10), conc(10), vdd(3)
  real, parameter :: avo = 6.022136736d23
end module commondefs




! =============================================================================
! MAIN PROGRAM
! =============================================================================
program md
  use commondefs
  use vardecs
  implicit none

  integer i, counter
  real sump, pavr, totei, totkini, msd, rp, r0
  real elapsedtime
  character(len=80) :: controlfile

  ! ---------------------------------------------------------------------------
  ! Initialisation
  ! ---------------------------------------------------------------------------
  call read_control(controlfile)
  call open_output(controlfile, iconf, ivac, iforce, imsd, imovie)

  call initvel(temp)

  ! Compute initial forces
  call forceeval(np, rcut, boxsize_x, boxsize_y, boxsize_z)

  ! Initialize linked list
  call new_nlist()

  counter = 0
  sump = 0.0
  tottime = 0.0

  ! ---------------------------------------------------------------------------
  ! Equilibration loop (steps 1 .. scalestep)
  ! ---------------------------------------------------------------------------
  do tstep = 1, scalestep
    ! Integrator (uncomment the desired one)
    call nh2(tstp)                     ! reversible Nose-Hoover Thermostat
    ! call anders2()                   ! Andersen pressure control
    ! call verlet2(tstp)               ! non-reversible Nose-Hoover Thermostat
    ! call correct(tstp)               ! Gear 5 value
    ! call berendsen()                 ! Berendsen pressure control
    ! call ver2(tstp)                  ! traditional velocity Verlet

    tottime = counter * tstp

    ! Temperature scaling (uncomment if needed)
    ! call velscale()

    ! Compute kinetic energy
    sumvel = 0.0
    do i = 1, n
      sumvel = sumvel + (atoms(i)%vx**2 + atoms(i)%vy**2 + atoms(i)%vz**2) * molwt(atoms(i)%kind)
    end do

    ! Properties
    totkin = 0.5 * sumvel
    counter = counter + 1
    tempdum = sumvel / (3.0 * real(n) * bolz)
    pressure = ( (n * bolz * tempdum) + (1.0/3.0) * tow ) / volume
    sump = sump + pressure
    pavr = sump / real(counter)
    totei = tote / real(n + nwall)
    totkini = totkin / real(n + nwall)
    hamilton = tote + totkin + (ps**2 * q) / 2.0 + real(g) * tempdum * bolz * s

    write(*, 2) tstep, tempdum, tote + totkin, hamilton
    2 format(i8, 3es13.5)

    ! Update linked list
    call new_nlist()

    ! Write properties to file
    if (enunit == 1) then
      write(5, 15) tottime, totei, totkini, totei + totkini + encorr, tempdum, pressure, hamilton, volume
      15 format(1x, 8es25.17)
    else if (enunit == 0) then
      write(5, 16) tottime, totei, totkini, totei + totkini + encorr, tempdum, pressure, hamilton, volume
      16 format(1x, 8es25.17)
    end if

    ! Write vacancy positions
    if (ivac > 0) then
      if (mod(tstep, ivac) == 0) then
        write(8, *) nvac
        write(8, *) " vacancy positions: ", " at step= ", tstep, " time= ", tottime, " fs "
        do i = n + nwall + 1, n + nwall + nvac
          write(8, 14) atoms(i)%x, atoms(i)%y, atoms(i)%z
          14 format(1x, 3f20.10)
        end do
      end if
    end if

    ! Write forces
    if (iforce > 0) then
      if (mod(tstep, iforce) == 0) then
        write(7, *) " forces ", " at step= ", tstep, " time= ", tottime, " fs "
        do i = 1, n + nwall + nvac
          write(7, 13) atoms(i)%fx, atoms(i)%fy, atoms(i)%fz, atoms(i)%elfx, atoms(i)%elfy, atoms(i)%elfz
          13 format(1x, 6es15.7)
        end do
      end if
    end if

    ! Write trajectory to .xyz file
    if (imovie > 0) then
      if (mod(tstep, imovie) == 0) then
        write(15, *) n + nwall
        write(15, *) 'md step ', tstep
        do i = 1, n + nwall
          write(15, '(i3, 3f10.4)') atoms(i)%kind, atoms(i)%x, atoms(i)%y, atoms(i)%z
        end do
      end if
    end if

  end do   ! end of equilibration loop

  ! ---------------------------------------------------------------------------
  ! Store initial positions for MSD
  ! ---------------------------------------------------------------------------
  do i = 1, n
    atoms(i)%x0 = atoms(i)%x
    atoms(i)%y0 = atoms(i)%y
    atoms(i)%z0 = atoms(i)%z
  end do

  ! ---------------------------------------------------------------------------
  ! Production loop (steps scalestep+1 .. totstep)
  ! ---------------------------------------------------------------------------
  do tstep = scalestep + 1, totstep
    ! Store first positions for MSD
    do i = 1, n
      atoms(i)%xb = atoms(i)%x
      atoms(i)%yb = atoms(i)%y
      atoms(i)%zb = atoms(i)%z
    end do

    ! Integrator predictor (uncomment the desired one)
    call nh1(tstp)                     ! reversible Nose-Hoover Thermostat
    ! call predict(tstp)               ! Gear 5 value
    ! call verlet1(tstp)               ! non-reversible Nose-Hoover Thermostat
    ! call ver1(tstp)                  ! traditional velocity Verlet

    ! Periodic boundary conditions (x and z only, y commented)
    do i = 1, n
      if (atoms(i)%x < 0.0) then
        atoms(i)%x = atoms(i)%x + boxsize_x
      else if (atoms(i)%x > boxsize_x) then
        atoms(i)%x = atoms(i)%x - boxsize_x
      end if
      ! Uncomment for PBC in y
      ! if (atoms(i)%y < 0.0) then
      !   atoms(i)%y = atoms(i)%y + boxsize_y
      ! else if (atoms(i)%y > boxsize_y) then
      !   atoms(i)%y = atoms(i)%y - boxsize_y
      ! end if
      if (atoms(i)%z < 0.0) then
        atoms(i)%z = atoms(i)%z + boxsize_z
      else if (atoms(i)%z > boxsize_z) then
        atoms(i)%z = atoms(i)%z - boxsize_z
      end if
    end do

    ! Evaluate forces and energy
    call forceeval(np, rcut, boxsize_x, boxsize_y, boxsize_z)

    ! Electromigration force (if step matches)
    if (mod(tstep, ielmigstep) == 0) then
      call findvac()
      call elmig()
    end if

    ! Add EM force to MD force
    do i = 1, n
      atoms(i)%fx = atoms(i)%fx + atoms(i)%elfx
      atoms(i)%fy = atoms(i)%fy + atoms(i)%elfy
      atoms(i)%fz = atoms(i)%fz + atoms(i)%elfz
    end do

    ! Integrator corrector (uncomment the one matching predictor)
    call nh2(tstp)                     ! reversible Nose-Hoover Thermostat
    ! call verlet2(tstp)               ! non-reversible Nose-Hoover Thermostat
    ! call correct(tstp)               ! Gear 5 value
    ! call ver2(tstp)                  ! traditional velocity Verlet

    tottime = counter * tstp

    ! Mean square displacement
    if (imsd /= 0) then
      if (mod(tstep, imsd) == 0) then
        msd = 0.0
        ! Rewind PBC for MSD
        do i = 1, n
          atoms(i)%xp = atoms(i)%x + dnint((atoms(i)%xb - atoms(i)%x) / boxsize_x) * boxsize_x
          atoms(i)%yp = atoms(i)%y + dnint((atoms(i)%yb - atoms(i)%y) / boxsize_y) * boxsize_y
          atoms(i)%zp = atoms(i)%z + dnint((atoms(i)%zb - atoms(i)%z) / boxsize_z) * boxsize_z
        end do
        do i = msdstart, msdend
          rp = dsqrt(atoms(i)%xp**2 + atoms(i)%yp**2 + atoms(i)%zp**2)
          r0 = dsqrt(atoms(i)%x0**2 + atoms(i)%y0**2 + atoms(i)%z0**2)
          msd = msd + (rp - r0)**2
        end do
        msd = msd / real(msdend - msdstart + 1)
        write(4, '(f20.13,1x,es25.17)') real(tstep - (scalestep + imsd)) * tstp, msd
      end if
    end if

    ! Kinetic energy
    sumvel = 0.0
    do i = 1, n
      sumvel = sumvel + (atoms(i)%vx**2 + atoms(i)%vy**2 + atoms(i)%vz**2) * molwt(atoms(i)%kind)
    end do

    ! Properties
    totkin = 0.5 * sumvel
    counter = counter + 1
    tempdum = sumvel / (3.0 * real(n) * bolz)
    pressure = ( (n * bolz * tempdum) + (1.0/3.0) * tow ) / volume
    sump = sump + pressure
    pavr = sump / real(counter)
    totei = tote / real(n + nwall)
    totkini = totkin / real(n + nwall)
    hamilton = tote + totkin + (ps**2 * q) / 2.0 + real(g) * tempdum * bolz * s

    write(*, 3) tstep, tempdum, tote, hamilton
    3 format(i8, 3es13.5)

    ! Update neighbor list
    call new_nlist()

    ! Write properties
    if (enunit == 1) then
      write(5, 18) tottime, totei, totkini, totei + totkini + encorr, tempdum, pressure, hamilton, volume, msd
      18 format(1x, 9es25.17)
    else if (enunit == 0) then
      write(5, 19) tottime, totei, totkini, totei + totkini + encorr, tempdum, pressure, hamilton, volume, msd
      19 format(1x, 9es25.17)
    end if

    ! Write configuration
    if (iconf > 0) then
      if (mod(tstep, iconf) == 0) then
        write(9, *) tstep
        do i = 1, n
          write(9, 21) atoms(i)%x, atoms(i)%y, atoms(i)%z, atoms(i)%vx, atoms(i)%vy, &
                        atoms(i)%vz, atoms(i)%ax, atoms(i)%ay, atoms(i)%az
          21 format(1x, 9es15.7)
        end do
      end if
    end if

    ! Write forces
    if (iforce > 0) then
      if (mod(tstep, iforce) == 0) then
        write(7, *) " forces ", " at step= ", tstep, " time= ", tottime, " fs "
        do i = 1, n + nvac + nwall
          write(7, 24) atoms(i)%fx, atoms(i)%fy, atoms(i)%fz, atoms(i)%elfx, atoms(i)%elfy, atoms(i)%elfz
          24 format(1x, 6es15.7)
        end do
      end if
    end if

    ! Write vacancy positions
    if (ivac > 0) then
      if (mod(tstep, ivac) == 0) then
        write(8, *) nvac
        write(8, *) " vacancy positions: ", " at step= ", tstep, " time= ", tottime, " fs "
        do i = n + nwall + 1, n + nvac + nwall
          write(8, 22) atoms(i)%x, atoms(i)%y, atoms(i)%z
          22 format(1x, 3f20.10)
        end do
      end if
    end if

    ! Write trajectory to .xyz
    if (imovie > 0) then
      if (mod(tstep, imovie) == 0) then
        write(15, *) n + nwall
        write(15, *) 'md step ', tstep
        do i = 1, n + nwall
          write(15, "(i3,3f10.4)") atoms(i)%kind, atoms(i)%x, atoms(i)%y, atoms(i)%z
        end do
      end if
    end if

  end do   ! end of production loop

  ! ---------------------------------------------------------------------------
  ! Finalisation
  ! ---------------------------------------------------------------------------
  if (iconf > 0) close(9)
  if (imsd > 0) close(4)
  if (iforce > 0) close(7)
  if (ivac > 0) close(8)
  if (imovie > 0) close(15)

  ! Velocity autocorrelation function (requires external routine)
  ! call vacf(controlfile, totstep, scalestep, tstp, n, iconf)

  elapsedtime = timef()   ! external function
  write(*, *) "total_elapsed_time_is: ", elapsedtime
  write(5, *) "total_elapsed_time_is: ", elapsedtime, "secs."
  close(5)

end program md


! =============================================================================
! SUBROUTINE: open_output
! =============================================================================
subroutine open_output(name1, ic, iv, ifrc, ims, imov)
  integer ic, iv, ifrc, ims, imov
  character name1*(*)
  character(len=len_trim(name1)) name2

  name2 = trim(name1)
  open(unit=5, file=name2//"out")
  if (ic > 0) open(unit=9, file=name2//"con")
  if (iv > 0) open(unit=8, file=name2//"vac")
  if (ifrc > 0) open(unit=7, file=name2//"elf")
  if (ims > 0) open(unit=4, file=name2//"msd")
  if (imov > 0) open(unit=15, file=name2//"xyz")
end subroutine open_output

! =============================================================================
! SUBROUTINE: read_control
! =============================================================================
subroutine read_control(name1)
  use vardecs
  use commondefs
  implicit none

  integer i, npot, pottype, at1, at2, bit
  character name1*(*), name*22
  character(len=len_trim(name1)) :: name2

  name2 = trim(name1)
  open(unit=10, file=name2)
  open(unit=11, file='fcc.txt')

  read(10, *) n, nvac, nwall
  read(10, *) ntype
  read(10, *) totstep, scalestep, scalelint
  read(10, *) deltat
  read(10, *) tol
  read(10, *) msdstart, msdend
  read(10, *) temp
  read(10, *) q
  read(10, *) iconf, ivac, iforce, imsd, imovie, ielmigstep
  read(10, *) enunit
  read(10, *) a, boxsize_x, boxsize_y, boxsize_z, rcut
  read(10, *) (molwt(i), i = 1, ntype)
  do i = 1, ntype
    read(10, *) atno(i), conc(i)
  end do
  read(10, *) vdd(1), vdd(2), vdd(3)

  allocate(atoms(n + nvac + nwall + 150))
  allocate(ions(n + nvac + nwall + 150))

  do i = 1, 50
    pairtype(i)%pot = 0
  end do

  read(10, *) npot
  do i = 1, npot
    read(10, *) pottype
    select case (pottype)
    case (1)   ! Lennard-Jones
      read(10, *) at1, at2, ljpot(i)%eps, ljpot(i)%sigma
      bit = 2**at1 + 2**at2
      pairtype(bit)%index = i
      pairtype(bit)%pot = 1
    case (2)   ! Fourier potential
      read(10, *) at1, at2, fourierpot(i)%a, fourierpot(i)%b, fourierpot(i)%c, &
                  fourierpot(i)%d, fourierpot(i)%e, fourierpot(i)%f, fourierpot(i)%g, &
                  fourierpot(i)%h, fourierpot(i)%ix, fourierpot(i)%jx, fourierpot(i)%kx, &
                  fourierpot(i)%lx, fourierpot(i)%mx, fourierpot(i)%rmin, fourierpot(i)%rmax, &
                  fourierpot(i)%eps, fourierpot(i)%sigma
      bit = 2**at1 + 2**at2
      pairtype(bit)%index = i
      pairtype(bit)%pot = 2
    end select
  end do

  nbox_x = int(boxsize_x / a)
  nbox_y = int(boxsize_y / a)
  nbox_z = int(boxsize_z / a)
  nperf = nbox_x * nbox_y * nbox_z * 4

  ! Read atomic coordinates
  do i = 1, n + nvac + nwall
    read(10, 10) atoms(i)%kind, atoms(i)%x, atoms(i)%y, atoms(i)%z
  end do
10 format(i2, 1x, f20.17, 4x, f20.17, 4x, f20.17)

  allocate(xini(nperf), yini(nperf), zini(nperf))
  close(11)

  ncellx = int(boxsize_x / rcut / 1.01)
  ncelly = int(boxsize_y / rcut / 1.01)
  ncellz = int(boxsize_z / rcut / 1.01)
  cell_factorx = real(ncellx) / boxsize_x
  cell_factory = real(ncelly) / boxsize_y
  cell_factorz = real(ncellz) / boxsize_z

  call new_nlist()
  call neighbour_cells()

end subroutine read_control

! =============================================================================
! SUBROUTINE: initvel
! =============================================================================
subroutine initvel(temper)
  use commondefs
  use vardecs
  implicit none

  real temper, rtemp, sumx, sumy, sumz, scale, sumvel2, aheat2
  real dummy, gauss
  integer i

  dummy = 0.5
  rtemp = sqrt(temper)
  aheat2 = real(n) * 3.0 * bolz * temper

  do i = 1, n
    atoms(i)%vx = rtemp * gauss(dummy)
    atoms(i)%vy = rtemp * gauss(dummy)
    atoms(i)%vz = rtemp * gauss(dummy)
  end do

  sumx = 0.0; sumy = 0.0; sumz = 0.0
  do i = 1, n
    sumx = sumx + atoms(i)%vx
    sumy = sumy + atoms(i)%vy
    sumz = sumz + atoms(i)%vz
  end do
  sumx = sumx / real(n)
  sumy = sumy / real(n)
  sumz = sumz / real(n)

  do i = 1, n
    atoms(i)%vx = atoms(i)%vx - sumx
    atoms(i)%vy = atoms(i)%vy - sumy
    atoms(i)%vz = atoms(i)%vz - sumz
  end do

  sumvel2 = 0.0
  do i = 1, n
    sumvel2 = sumvel2 + (atoms(i)%vx**2 + atoms(i)%vy**2 + atoms(i)%vz**2) * molwt(atoms(i)%kind)
  end do

  scale = sqrt(aheat2 / sumvel2)
  sumvel = 0.0
  do i = 1, n
    atoms(i)%vx = atoms(i)%vx * scale
    atoms(i)%vy = atoms(i)%vy * scale
    atoms(i)%vz = atoms(i)%vz * scale
    sumvel = sumvel + (atoms(i)%vx**2 + atoms(i)%vy**2 + atoms(i)%vz**2) * molwt(atoms(i)%kind)
  end do

end subroutine initvel

! =============================================================================
! FUNCTION: gauss (returns Gaussian random number)
! =============================================================================
real function gauss(dummy)
  real a1, a3, a5, a7, a9
  parameter (a1 = 3.949846138, a3 = 0.252408784)
  parameter (a5 = 0.076542912, a7 = 0.008355968)
  parameter (a9 = 0.029899776)
  real sum, r, r2
  real ranf, dummy
  integer i

  sum = 0.0
  do i = 1, 12
    sum = sum + ranf(dummy)
  end do
  r = (sum - 6.0) / 4.0
  r2 = r * r
  gauss = ((((a9 * r2 + a7) * r2 + a5) * r2 + a3) * r2 + a1) * r
  return
end function gauss

! -----------------------------------------------------------------------------
! Note: the function ranf() is not defined; it is expected to return a random
!       number in [0,1]. The user must provide an implementation.
! -----------------------------------------------------------------------------

! =============================================================================
! SUBROUTINE: forceeval
! =============================================================================
subroutine forceeval(npart, rcut, box_x, box_y, box_z)
  use commondefs
  use vardecs
  implicit none

  integer npart
  real rcut, box_x, box_y, box_z
  real xi, yi, zi, rxij, ryij, rzij, xr, yr, zr, rijsq, rijsqrt, rij, rpl
  real fxij, fyij, fzij, fij
  real rcutsq2, box_x2, box_y2, box_z2
  integer i, j, ic, inn, in, bit

  integer, external :: find_cell

  tote = 0.0
  tow = 0.0
  rcutsq2 = rcut**2
  box_x2 = box_x / 2.0
  box_y2 = box_y / 2.0
  box_z2 = box_z / 2.0

  do i = 1, npart + nwall
    atoms(i)%fx = 0.0
    atoms(i)%fy = 0.0
    atoms(i)%fz = 0.0
  end do

  do i = 1, npart + nwall - 1
    xi = atoms(i)%x
    yi = atoms(i)%y
    zi = atoms(i)%z
    ic = find_cell(xi, yi, zi)

    do inn = int(ic * 27 + 1), int(ic * 27 + 27)
      in = n_v_list(inn)
      j = head_chain(in)
      do while (j /= 0)
        if ((j /= i) .and. (j >= i + 1)) then
          bit = 2**atoms(i)%kind + 2**atoms(j)%kind
          if (pairtype(bit)%pot == 0) then
            j = linked_list(j)
            cycle
          end if

          rxij = xi - atoms(j)%x
          ryij = yi - atoms(j)%y
          rzij = zi - atoms(j)%z

          ! Minimum image convention
          xr = abs(rxij)
          yr = abs(ryij)
          zr = abs(rzij)
          if (xr > box_x2) rxij = (xr - box_x) * sign(1.0, rxij)
          if (yr > box_y2) ryij = (yr - box_y) * sign(1.0, ryij)
          if (zr > box_z2) rzij = (zr - box_z) * sign(1.0, rzij)

          rijsq = rxij*rxij + ryij*ryij + rzij*rzij
          rijsqrt = sqrt(rijsq)

          if (rijsq < rcutsq2) then
            select case (pairtype(bit)%pot)
            case (1)   ! Lennard-Jones
              fij = lennardforce(rijsqrt, bit)
              tote = tote + lennardenergy(rijsqrt, bit)
            case (2)   ! Fourier
              fij = forcefunc(rijsqrt, bit)
              tote = tote + potfunc(rijsqrt, bit)
            end select

            rpl = fij * rijsqrt
            fxij = fij * rxij / rijsqrt
            fyij = fij * ryij / rijsqrt
            fzij = fij * rzij / rijsqrt

            atoms(i)%fx = atoms(i)%fx + fxij
            atoms(j)%fx = atoms(j)%fx - fxij
            atoms(i)%fy = atoms(i)%fy + fyij
            atoms(j)%fy = atoms(j)%fy - fyij
            atoms(i)%fz = atoms(i)%fz + fzij
            atoms(j)%fz = atoms(j)%fz - fzij
            tow = tow + rpl
          end if
        end if
        j = linked_list(j)
      end do
    end do
  end do

end subroutine forceeval

! =============================================================================
! FUNCTION: potfunc (Fourier potential energy)
! =============================================================================
function potfunc(rij, bit)
  use vardecs
  implicit none
  real potfunc, rij, rx
  integer bit
  real pi
  parameter (pi = 3.141592653589793)

  rx = pi * (rij - fourierpot(pairtype(bit)%index)%rmin) / &
       (fourierpot(pairtype(bit)%index)%rmax - fourierpot(pairtype(bit)%index)%rmin)

  potfunc = fourierpot(pairtype(bit)%index)%a + &
            fourierpot(pairtype(bit)%index)%b * cos(rx) + &
            fourierpot(pairtype(bit)%index)%c * sin(rx) + &
            fourierpot(pairtype(bit)%index)%d * cos(2.0*rx) + &
            fourierpot(pairtype(bit)%index)%e * sin(2.0*rx) + &
            fourierpot(pairtype(bit)%index)%f * cos(3.0*rx) + &
            fourierpot(pairtype(bit)%index)%g * sin(3.0*rx) + &
            fourierpot(pairtype(bit)%index)%h * cos(4.0*rx) + &
            fourierpot(pairtype(bit)%index)%ix * sin(4.0*rx) + &
            fourierpot(pairtype(bit)%index)%jx * cos(5.0*rx) + &
            fourierpot(pairtype(bit)%index)%kx * sin(5.0*rx) + &
            fourierpot(pairtype(bit)%index)%lx * cos(6.0*rx) + &
            fourierpot(pairtype(bit)%index)%mx * sin(6.0*rx)
  return
end function potfunc

! =============================================================================
! FUNCTION: forcefunc (Fourier force)
! =============================================================================
function forcefunc(rij2, bit)
  use vardecs
  implicit none
  real forcefunc, rij2, y, ry
  integer bit
  real pi
  parameter (pi = 3.141592653589793)

  ry = pi * (rij2 - fourierpot(pairtype(bit)%index)%rmin) / &
       (fourierpot(pairtype(bit)%index)%rmax - fourierpot(pairtype(bit)%index)%rmin)
  y = pi / (fourierpot(pairtype(bit)%index)%rmax - fourierpot(pairtype(bit)%index)%rmin)

  forcefunc = fourierpot(pairtype(bit)%index)%b * sin(ry) * y - &
              fourierpot(pairtype(bit)%index)%c * cos(ry) * y + &
              2.0 * fourierpot(pairtype(bit)%index)%d * sin(2.0*ry) * y - &
              2.0 * fourierpot(pairtype(bit)%index)%e * cos(2.0*ry) * y + &
              3.0 * fourierpot(pairtype(bit)%index)%f * sin(3.0*ry) * y - &
              3.0 * fourierpot(pairtype(bit)%index)%g * cos(3.0*ry) * y + &
              4.0 * fourierpot(pairtype(bit)%index)%h * sin(4.0*ry) * y - &
              4.0 * fourierpot(pairtype(bit)%index)%ix * cos(4.0*ry) * y + &
              5.0 * fourierpot(pairtype(bit)%index)%jx * sin(5.0*ry) * y - &
              5.0 * fourierpot(pairtype(bit)%index)%kx * cos(5.0*ry) * y + &
              6.0 * fourierpot(pairtype(bit)%index)%lx * sin(6.0*ry) * y - &
              6.0 * fourierpot(pairtype(bit)%index)%mx * cos(6.0*ry) * y
  return
end function forcefunc

! =============================================================================
! FUNCTION: lennardforce
! =============================================================================
function lennardforce(rij3, bit)
  use vardecs
  implicit none
  real lennardforce, rij3
  integer bit

  lennardforce = 24.0 * ljpot(pairtype(bit)%index)%eps * &
                 (2.0 * (ljpot(pairtype(bit)%index)%sigma**12 / rij3**13) - &
                        (ljpot(pairtype(bit)%index)%sigma**6  / rij3**7))
  return
end function lennardforce

! =============================================================================
! FUNCTION: lennardenergy
! =============================================================================
function lennardenergy(rij, bit)
  use vardecs
  implicit none
  real lennardenergy, rij
  integer bit

  lennardenergy = 4.0 * ljpot(pairtype(bit)%index)%eps * &
                  ((ljpot(pairtype(bit)%index)%sigma / rij)**12 - &
                   (ljpot(pairtype(bit)%index)%sigma / rij)**6)
  return
end function lennardenergy

! =============================================================================
! FUNCTION: find_cell
! =============================================================================
integer function find_cell(ionx, iony, ionz)
  use commondefs
  implicit none
  real ionx, iony, ionz
  integer ix, iy, iz

  ix = int(ionx * cell_factorx)
  iy = int(iony * cell_factory)
  iz = int(ionz * cell_factorz)
  find_cell = ix + iy * ncellx + iz * ncellx * ncelly
end function find_cell

! =============================================================================
! SUBROUTINE: neighbour_cells
! =============================================================================
subroutine neighbour_cells()
  use commondefs
  implicit none
  integer ix, iy, iz, ic, in, itel, icx, icy, icz, iccx, iccy, iccz

  do ic = 0, ncellx * ncelly * ncellz - 1
    itel = 0
    icz = ic / (ncellx * ncelly)
    icy = (ic - icz * ncellx * ncelly) / ncellx
    icx = ic - icy * ncellx - icz * ncellx * ncelly

    do iz = -1, 1
      iccz = icz + iz
      if (iccz < 0) then
        iccz = iccz + ncellz
      else if (iccz >= ncellz) then
        iccz = iccz - ncellz
      end if

      do iy = -1, 1
        iccy = icy + iy
        if (iccy < 0) then
          iccy = iccy + ncelly
        else if (iccy >= ncelly) then
          iccy = iccy - ncelly
        end if

        do ix = -1, 1
          iccx = icx + ix
          if (iccx < 0) then
            iccx = iccx + ncellx
          else if (iccx >= ncellx) then
            iccx = iccx - ncellx
          end if

          in = iccx + iccy * ncellx + iccz * ncellx * ncelly
          itel = itel + 1
          n_v_list(int(ic * 27 + itel)) = in
        end do
      end do
    end do
  end do
end subroutine neighbour_cells

! =============================================================================
! SUBROUTINE: new_nlist
! =============================================================================
subroutine new_nlist()
  use commondefs
  implicit none
  integer i, ic
  integer, external :: find_cell

  do ic = 0, ncellx * ncelly * ncellz - 1
    head_chain(ic) = 0
  end do

  do i = 1, n + nwall
    ic = find_cell(atoms(i)%x, atoms(i)%y, atoms(i)%z)
    linked_list(i) = head_chain(ic)
    head_chain(ic) = i
  end do
end subroutine new_nlist

! =============================================================================
! The following routines are related to electromigration.
! They require a module 'defs' which is not provided in the document.
! The code is kept as is, but the user must supply the missing module.
! =============================================================================

! -----------------------------------------------------------------------------
! SUBROUTINE: init (electromigration initialisation)
! -----------------------------------------------------------------------------
subroutine init()
  use defs
  use commondefs
  implicit none
  integer i, j

  vd = sqrt(vdd(1)**2 + vdd(2)**2 + vdd(3)**2)
  vddir(1) = vdd(1) / vd
  vddir(2) = vdd(2) / vd
  vddir(3) = vdd(3) / vd

  do i = 1, ntype
    do j = 1, 56
      if (atno(i) == el(j)%atno) sp(i) = el(j)
    end do
  end do

  box_x = boxsize_x / 0.529177
  box_y = boxsize_y / 0.529177
  box_z = boxsize_z / 0.529177
  cutoff = rcut / 0.529177

  do i = 1, n + nvac + nwall
    ions(i)%x = atoms(i)%x / 0.529177249
    ions(i)%y = atoms(i)%y / 0.529177249
    ions(i)%z = atoms(i)%z / 0.529177249
    ions(i)%kind = atoms(i)%kind
  end do

  close(6)
  omave = 0.0
  zave = 0.0
  mave = 0.0
  do i = 1, ntype
    omave = omave + conc(i) * sp(i)%om
    zave = zave + conc(i) * sp(i)%z
    mave = mave + conc(i) * sp(i)%m
  end do
  mave = 1.0

  kf = (3.0 * pi * pi * zave / omave)**(1.0/3.0)
  vf = kf / mave
  ef = kf * kf / 2.0 / mave
  rs = mave * (3.0 * omave / 4.0 / pi / zave)**(1.0/3.0)
  ksi = 2.0 / (1.0 + 0.02544716775 * rs)

  do i = 1, n + nwall
    atoms(i)%elfx = 0.0
    atoms(i)%elfy = 0.0
    atoms(i)%elfz = 0.0
  end do
end subroutine init

! -----------------------------------------------------------------------------
! FUNCTION: cross (cosine of angle between vd and Rj-Rm)
! -----------------------------------------------------------------------------
real function cross(rx, ry, rz)
  use defs
  implicit none
  real rx, ry, rz, rr
  rr = sqrt(rx*rx + ry*ry + rz*rz)
  cross = (rx/rr)*vddir(1) + (ry/rr)*vddir(2) + (rz/rr)*vddir(3)
end function cross

! -----------------------------------------------------------------------------
! FUNCTION: dotproduct
! -----------------------------------------------------------------------------
real function dotproduct(ii, var)
  use defs
  use commondefs
  implicit none
  real var
  integer ii

  if (ii == 1) then
    dotproduct = vdd(1)*var + vdd(2)*y + vdd(3)*z
  else if (ii == 2) then
    dotproduct = vdd(1)*x + vdd(2)*var + vdd(3)*z
  else if (ii == 3) then
    dotproduct = vdd(1)*x + vdd(2)*y + vdd(3)*var
  end if
end function dotproduct

! -----------------------------------------------------------------------------
! FUNCTION: bessel (spherical Bessel function)
! -----------------------------------------------------------------------------
real function bessel(q)
  use defs
  implicit none
  real q, xx
  xx = q * rr
  bessel = sin(xx) / xx / xx - cos(xx) / xx
end function bessel

! -----------------------------------------------------------------------------
! FUNCTION: integrant1 (used in integration)
! -----------------------------------------------------------------------------
real function integrant1(q)
  use defs
  use commondefs
  implicit none
  real q, wi, wj
  real, external :: w, bessel

  wj = w(sp(ions(force_j)%kind), q)
  if (ions(force_i)%kind < 0) then
    wi = w(sp(1), q)
  else
    wi = w(sp(ions(force_i)%kind), q)
  endif
  integrant1 = wi * wj * bessel(q)
end function integrant1

! -----------------------------------------------------------------------------
! FUNCTION: uu_s (part of EM force)
! -----------------------------------------------------------------------------
function uu_s(ii, var)
  use defs
  use commondefs
  implicit none
  real uu_s, var, uu, term1, term2, term3, t1
  integer ii, i
  real, external :: integrate, integrant1, integrant2, integrant3, cross, dotproduct

  uu = 0.0
  do i = 1, n
    if (i == force_j) cycle
    t1 = 1.0 / (sp(ions(i)%kind)%z)
    if (sp(ions(force_j)%kind)%om) then
      if (ions(force_i)%kind < 0) then
        term1 = -t1 * vd / vf * cross(rx,ry,rz) * integrate(integrant1)
        uu = uu + term1
      else
        term1 = t1 * vd / vf * cross(rx,ry,rz) * integrate(integrant1)
        uu = uu + term1
      endif
    else
      uu = uu + term1
    endif
    if (ions(force_i)%kind > 1) then
      term3 = t1 * vd / vf * cross(-rx, -ry, -rz) * integrate(integrant3)
      uu = uu - term3
    endif
    if (i == force_i) then
      term2 = t1 / 3.0 * dotproduct(ii, var) / vf * integrate(integrant2)
      uu = uu + term2
    endif
  enddo
  uu_s = uu
end function uu_s

! -----------------------------------------------------------------------------
! SUBROUTINE: elmig (main EM force routine)
! -----------------------------------------------------------------------------
subroutine elmig()
  use defs
  use commondefs
  implicit none
  real error
  real, external :: uu_s, dfridr, init
  !include 'data.h'

  call init()
  do force_j = 1, n
    x = ions(force_j)%x
    y = ions(force_j)%y
    z = ions(force_j)%z
    ! Convert to eV/A
    atoms(force_j)%elfx = -dfridr(1, x, 0.2, error) * 51.42210754
    atoms(force_j)%elfy = -dfridr(2, y, 0.2, error) * 51.42210754
    atoms(force_j)%elfz = -dfridr(3, z, 0.2, error) * 51.42210754
  end do
end subroutine elmig

! -----------------------------------------------------------------------------
! FUNCTION: bofq (Heine-Abarenkov-Animalu psp N(q))
! -----------------------------------------------------------------------------
real function bofq(ee, q)
  use defs
  implicit none
  real q, qrm, qrc, t1, t2, t3, t4, t5, t6, t7
  type(elements) :: ee

  qrm = q * ee%rm
  qrc = q * ee%rc
  t1 = pi / omave / q / q
  t2 = 8.0 * t1 * ee%c / q
  t3 = sin(qrm) - qrm * cos(qrm)
  t4 = 8.0 * t1 * ee%z * cos(qrm)
  t5 = 4.0 * t1 * ee%c / q
  t6 = 24.0 * t1 * ee%z * ee%alf / qrc / qrc / qrc
  t7 = sin(qrc) - qrc * cos(qrc)
  bofq = t2 * t3 - t4 + (t5 - t6) * t7
end function bofq

! -----------------------------------------------------------------------------
! FUNCTION: nofq (Heine-Abarenkov-Animalu psp N(q) – second part)
! -----------------------------------------------------------------------------
real function nofq(ee, q)
  use defs
  implicit none
  real q, t1, t2, t3, t4, xx, yy, cost, costp
  real j0x, j1x, j2x, j3x, j0y, j1y, j2y, j3y
  type(elements) :: ee

  t1 = pi * ee%rm * ee%rm * ee%rm / omave
  xx = kf * ee%rm
  j0x = sin(xx) / xx
  j1x = sin(xx) / xx / xx - cos(xx) / xx
  j2x = 3.0 / xx * j1x - j0x
  j3x = 5.0 / xx * j2x - j1x

  if ((q - 2.0 * kf) <= 0.0) then
    cost = (1.0 - q*q / 2.0 / kf / kf)
    t2 = -4.0 * t1 * (ee%a0 - ee%c) * (j0x*j0x - j1x*cos(xx)/xx)
    t3 = -12.0 * t1 * (ee%a1 - ee%c) * (j1x*j1x - j0x*j2x) * cost
    t4 = -20.0 * t1 * (ee%a2 - ee%c) * (j2x*j2x - j1x*j3x) * (3.0 * cost*cost - 1.0) / 2.0
    nofq = t2 + t3 + t4
  else
    nofq = 0.0
  endif
end function nofq

! =============================================================================
! SUBROUTINE: findvac (find vacancies)
! =============================================================================
subroutine findvac()
  use commondefs
  implicit none
  integer i, j, dum
  logical flag1
  real rxij, ryij, rzij
  real, allocatable :: xini(:), yini(:), zini(:)

  allocate(xini(nperf), yini(nperf), zini(nperf))
  call fcc2(a, int(nbox_x), int(nbox_y), int(nbox_z))

  rewind(11)
  do j = 1, nperf
    read(11, *) xini(j), yini(j), zini(j)
  end do
  close(11)

  dum = 1
  flag1 = .false.
  do j = 1, nperf
    do i = 1, n + nwall
      rxij = xini(j) - atoms(i)%x
      ryij = yini(j) - atoms(i)%y
      rzij = zini(j) - atoms(i)%z
      if (sqrt(rxij**2 + ryij**2 + rzij**2) < tol * sqrt(3.0)) flag1 = .true.
    end do
    if (.not. flag1) then
      atoms(n + nwall + dum)%x = xini(j)
      atoms(n + nwall + dum)%y = yini(j)
      atoms(n + nwall + dum)%z = zini(j)
      atoms(n + nwall + dum)%kind = -1
      dum = dum + 1
    end if
    flag1 = .false.
  end do
  nvac = dum - 1
end subroutine findvac

! =============================================================================
! SUBROUTINE: fcc (generate FCC lattice)
! =============================================================================
subroutine fcc(a, nx, ny, nz)
  real a
  real, dimension(4) :: x0, y0, z0
  integer ix, iy, iz, i, nx, ny, nz
  real, allocatable :: x(:), y(:), z(:)

  allocate(x(nx*ny*nz*4), y(nx*ny*nz*4), z(nx*ny*nz*4))
  open(11, file="fcc.txt", position="rewind", status="old")

  x0(1) = 0.0; y0(1) = 0.0; z0(1) = 0.0
  x0(2) = 0.0; y0(2) = 0.5*a; z0(2) = 0.5*a
  x0(3) = 0.5*a; y0(3) = 0.0; z0(3) = 0.5*a
  x0(4) = 0.5*a; y0(4) = 0.5*a; z0(4) = 0.0

  do iy = -1, ny - 2
    do ix = 0, nx - 1
      do iz = 0, nz - 1
        do i = 1, 4
          x(i) = ix * a + x0(i)
          y(i) = iy * a + y0(i)
          z(i) = iz * a + z0(i)
          write(11, 10) x(i), y(i), z(i)
        end do
      end do
    end do
  end do
10 format(f20.17, 4x, f20.17, 4x, f20.17)
  close(11)
end subroutine fcc

! =============================================================================
! SUBROUTINE: fcc2 (alternative FCC generation)
! =============================================================================
subroutine fcc2(a, nx, ny, nz)
  real a
  real, dimension(4) :: x0, y0, z0
  integer ix, iy, iz, i, nx, ny, nz
  real, allocatable :: x(:), y(:), z(:)

  allocate(x(nx*ny*nz*4), y(nx*ny*nz*4), z(nx*ny*nz*4))
  open(11, file="fcc.txt", position="rewind", status="old")

  x0(1) = 0.0; y0(1) = 0.0; z0(1) = 0.0
  x0(2) = 0.0; y0(2) = 0.5*a; z0(2) = 0.5*a
  x0(3) = 0.5*a; y0(3) = 0.0; z0(3) = 0.5*a
  x0(4) = 0.5*a; y0(4) = 0.5*a; z0(4) = 0.0

  do iy = -1, ny - 2
    do ix = 0, nx - 1
      do iz = 0, nz - 1
        do i = 1, 4
          x(i) = ix * a + x0(i)
          y(i) = iy * a + y0(i)
          z(i) = iz * a + z0(i)
          write(11, 10) x(i), y(i), z(i)
        end do
      end do
    end do
  end do
10 format(f20.17, 4x, f20.17, 4x, f20.17)
  close(11)
end subroutine fcc2

! =============================================================================
! SUBROUTINE: velscale (velocity scaling for temperature control)
! =============================================================================
subroutine velscale()
  use vardecs
  use commondefs
  implicit none
  real scalefactor, scalecount
  integer i

  sumvel = 0.0
  do i = 1, n
    sumvel = sumvel + (atoms(i)%vx**2 + atoms(i)%vy**2 + atoms(i)%vz**2) * molwt(atoms(i)%kind)
  end do

  if (tstep <= scalestep) then
    scalecount = mod(real(tstep), real(scalelint))
    if (scalecount == 0.0) then
      scalefactor = sqrt(3.0 * real(n) * bolz * temp / sumvel)
      do i = 1, n
        atoms(i)%vx = atoms(i)%vx * scalefactor
        atoms(i)%vy = atoms(i)%vy * scalefactor
        atoms(i)%vz = atoms(i)%vz * scalefactor
      end do
    end if
  end if
end subroutine velscale

! =============================================================================
! SUBROUTINE: anders1 (Andersen pressure control – predictor)
! =============================================================================
subroutine anders1()
  use vardecs
  use commondefs
  implicit none
  real sumro
  integer i

  sumro = 0.0
  do i = 1, n
    atoms(i)%rox  = atoms(i)%x / boxsize_x
    atoms(i)%roy  = atoms(i)%y / boxsize_y
    atoms(i)%roz  = atoms(i)%z / boxsize_z
    atoms(i)%rox1 = atoms(i)%vx / boxsize_x
    atoms(i)%roy1 = atoms(i)%vy / boxsize_y
    atoms(i)%roz1 = atoms(i)%vz / boxsize_z
    sumro = sumro + (atoms(i)%rox1**2 + atoms(i)%roy1**2 + atoms(i)%roz1**2) * molwt(atoms(i)%kind)
  end do

  do i = 1, n
    atoms(i)%rox = atoms(i)%rox + tstp * atoms(i)%rox1 + 0.5 * tstp**2 * &
                   ( (atoms(i)%fx / molwt(atoms(i)%kind)) / boxsize_x - 2.0/3.0 * atoms(i)%rox1 * volume1 / volume )
    atoms(i)%roy = atoms(i)%roy + tstp * atoms(i)%roy1 + 0.5 * tstp**2 * &
                   ( (atoms(i)%fy / molwt(atoms(i)%kind)) / boxsize_y - 2.0/3.0 * atoms(i)%roy1 * volume1 / volume )
    atoms(i)%roz = atoms(i)%roz + tstp * atoms(i)%roz1 + 0.5 * tstp**2 * &
                   ( (atoms(i)%fz / molwt(atoms(i)%kind)) / boxsize_z - 2.0/3.0 * atoms(i)%roz1 * volume1 / volume )
    atoms(i)%rox1 = atoms(i)%rox1 + 0.5 * tstp * &
                    ( (atoms(i)%fx / molwt(atoms(i)%kind)) / boxsize_x - 2.0/3.0 * atoms(i)%rox1 * volume1 / volume )
    atoms(i)%roy1 = atoms(i)%roy1 + 0.5 * tstp * &
                    ( (atoms(i)%fy / molwt(atoms(i)%kind)) / boxsize_y - 2.0/3.0 * atoms(i)%roy1 * volume1 / volume )
    atoms(i)%roz1 = atoms(i)%roz1 + 0.5 * tstp * &
                    ( (atoms(i)%fz / molwt(atoms(i)%kind)) / boxsize_z - 2.0/3.0 * atoms(i)%roz1 * volume1 / volume )
  end do
end subroutine anders1

! =============================================================================
! SUBROUTINE: anders2 (Andersen pressure control – corrector)
! =============================================================================
subroutine anders2()
  use vardecs
  use commondefs
  implicit none
  integer i
  real sumro

  sumro = 0.0
  do i = 1, n
    atoms(i)%x = atoms(i)%rox * boxsize_x
    atoms(i)%y = atoms(i)%roy * boxsize_y
    atoms(i)%z = atoms(i)%roz * boxsize_z
    atoms(i)%vx = atoms(i)%rox1 * boxsize_x
    atoms(i)%vy = atoms(i)%roy1 * boxsize_y
    atoms(i)%vz = atoms(i)%roz1 * boxsize_z
    sumro = sumro + (atoms(i)%rox1**2 + atoms(i)%roy1**2 + atoms(i)%roz1**2) * molwt(atoms(i)%kind)
  end do

  totkin = 0.5 * sumro
  tempdum = (sumro) / (3.0 * real(n) * bolz)
  pressure = ((n * bolz * tempdum) + (1.0/3.0) * tow) / volume
  volume1 = volume1 + 0.5 * tstp * (pressure - pset) / m
end subroutine anders2

! =============================================================================
! SUBROUTINE: berendsen (Berendsen pressure and temperature control)
! =============================================================================
subroutine berendsen()
  use vardecs
  use commondefs
  implicit none
  integer i

  sumvel = 0.0
  do i = 1, n
    sumvel = sumvel + (atoms(i)%vx**2 + atoms(i)%vy**2 + atoms(i)%vz**2) * molwt(atoms(i)%kind)
  end do

  tempdum = sumvel / (3.0 * real(n) * bolz)
  pressure = ((n * bolz * tempdum) + (1.0/3.0) * tow) / volume
  mu = (1.0 + (tstp / tau_p) * beta * (pressure - pset))**(1.0/3.0)
  lambda = (1.0 + (tstp / tau_t) * (temp / tempdum - 1.0)) * 0.5

  do i = 1, n
    atoms(i)%vx = atoms(i)%vx + tstp * atoms(i)%fx / molwt(atoms(i)%kind)
    atoms(i)%vy = atoms(i)%vy + tstp * atoms(i)%fy / molwt(atoms(i)%kind)
    atoms(i)%vz = atoms(i)%vz + tstp * atoms(i)%fz / molwt(atoms(i)%kind)
    atoms(i)%vx = atoms(i)%vx * lambda
    atoms(i)%vy = atoms(i)%vy * lambda
    atoms(i)%vz = atoms(i)%vz * lambda
    atoms(i)%x = atoms(i)%x + tstp * atoms(i)%vx
    atoms(i)%y = atoms(i)%y + tstp * atoms(i)%vy
    atoms(i)%z = atoms(i)%z + tstp * atoms(i)%vz
    atoms(i)%x = atoms(i)%x * mu
    atoms(i)%y = atoms(i)%y * mu
    atoms(i)%z = atoms(i)%z * mu
  end do

  boxsize_x = boxsize_x * mu
  boxsize_y = boxsize_y * mu
  boxsize_z = boxsize_z * mu
  volume = volume * (mu**3)
end subroutine berendsen

! =============================================================================
! SUBROUTINE: predict (Gear 5th order predictor)
! =============================================================================
subroutine predict(dt)
  use vardecs
  use commondefs
  implicit none
  real dt, c1, c2, c3, c4
  integer i

  c1 = dt
  c2 = c1 * dt / 2.0
  c3 = c2 * dt / 3.0
  c4 = c3 * dt / 4.0

  do i = 1, n
    atoms(i)%x = atoms(i)%x + c1 * atoms(i)%vx + c2 * atoms(i)%ax + c3 * atoms(i)%bx + c4 * atoms(i)%cx
    atoms(i)%y = atoms(i)%y + c1 * atoms(i)%vy + c2 * atoms(i)%ay + c3 * atoms(i)%by + c4 * atoms(i)%cy
    atoms(i)%z = atoms(i)%z + c1 * atoms(i)%vz + c2 * atoms(i)%az + c3 * atoms(i)%bz + c4 * atoms(i)%cz
    atoms(i)%vx = atoms(i)%vx + c1 * atoms(i)%ax + c2 * atoms(i)%bx + c3 * atoms(i)%cx
    atoms(i)%vy = atoms(i)%vy + c1 * atoms(i)%ay + c2 * atoms(i)%by + c3 * atoms(i)%cy
    atoms(i)%vz = atoms(i)%vz + c1 * atoms(i)%az + c2 * atoms(i)%bz + c3 * atoms(i)%cz
    atoms(i)%ax = atoms(i)%ax + c1 * atoms(i)%bx + c2 * atoms(i)%cx
    atoms(i)%ay = atoms(i)%ay + c1 * atoms(i)%by + c2 * atoms(i)%cy
    atoms(i)%az = atoms(i)%az + c1 * atoms(i)%bz + c2 * atoms(i)%cz
    atoms(i)%bx = atoms(i)%bx + c1 * atoms(i)%cx
    atoms(i)%by = atoms(i)%by + c1 * atoms(i)%cy
    atoms(i)%bz = atoms(i)%bz + c1 * atoms(i)%cz
  end do
end subroutine predict

! =============================================================================
! SUBROUTINE: correct (Gear 5th order corrector)
! =============================================================================
subroutine correct(dt)
  use vardecs
  use commondefs
  implicit none
  real dt, axi, ayi, azi, corrx, corry, corrz
  real c1, c2, c3, c4
  real gear0, gear1, gear3, gear4
  real cr, cv, cb, cc
  integer i

  gear0 = 19.0 / 120.0
  gear1 = 3.0 / 4.0
  gear3 = 1.0 / 2.0
  gear4 = 1.0 / 12.0

  c1 = dt
  c2 = c1 * dt / 2.0
  c3 = c2 * dt / 3.0
  c4 = c3 * dt / 4.0
  cr = gear0 * c2
  cv = gear1 * c2 / c1
  cb = gear3 * c2 / c3
  cc = gear4 * c2 / c4

  do i = 1, n
    axi = atoms(i)%fx / molwt(atoms(i)%kind)
    ayi = atoms(i)%fy / molwt(atoms(i)%kind)
    azi = atoms(i)%fz / molwt(atoms(i)%kind)
    corrx = axi - atoms(i)%ax
    corry = ayi - atoms(i)%ay
    corrz = azi - atoms(i)%az

    atoms(i)%x = atoms(i)%x + cr * corrx
    atoms(i)%y = atoms(i)%y + cr * corry
    atoms(i)%z = atoms(i)%z + cr * corrz
    atoms(i)%vx = atoms(i)%vx + cv * corrx
    atoms(i)%vy = atoms(i)%vy + cv * corry
    atoms(i)%vz = atoms(i)%vz + cv * corrz
    atoms(i)%ax = axi
    atoms(i)%ay = ayi
    atoms(i)%az = azi
    atoms(i)%bx = atoms(i)%bx + cb * corrx
    atoms(i)%by = atoms(i)%by + cb * corry
    atoms(i)%bz = atoms(i)%bz + cb * corrz
    atoms(i)%cx = atoms(i)%cx + cc * corrx
    atoms(i)%cy = atoms(i)%cy + cc * corry
    atoms(i)%cz = atoms(i)%cz + cc * corrz
  end do
end subroutine correct

! =============================================================================
! SUBROUTINE: verll (velocity Verlet predictor)
! =============================================================================
subroutine verll(dt)
  use vardecs
  use commondefs
  implicit none
  real dt, dt2, dtsq2
  integer i

  dt2 = dt / 2.0
  dtsq2 = dt * dt2

  do i = 1, n
    atoms(i)%x = atoms(i)%x + dt * atoms(i)%vx + dtsq2 * (atoms(i)%fx / molwt(atoms(i)%kind))
    atoms(i)%y = atoms(i)%y + dt * atoms(i)%vy + dtsq2 * (atoms(i)%fy / molwt(atoms(i)%kind))
    atoms(i)%z = atoms(i)%z + dt * atoms(i)%vz + dtsq2 * (atoms(i)%fz / molwt(atoms(i)%kind))
  end do
end subroutine verll

! =============================================================================
! NOTE: The following subroutines/functions are referenced but not provided in
!       the document. The user must supply implementations:
!       - nh1, nh2 (Nose-Hoover thermostat)
!       - verlet1, verlet2 (Nose-Hoover thermostat)
!       - ver1, ver2 (velocity Verlet)
!       - vacf (velocity autocorrelation function)
!       - timef (elapsed time function)
!       - dfridr (numerical derivative)
!       - integrate (numerical integration)
!       - w, ranf, and other functions from module 'defs'
!       - Module 'defs' containing types: elements, sp, el, etc.
! =============================================================================

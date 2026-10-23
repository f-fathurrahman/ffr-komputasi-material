! APPENDIX A
! LIST OF COMPUTER PROGRAM


!! MD program
program md
    use vardecs
    use commondefs
    implicit none

    real :: rp, r0, msd
    real :: deltatsq, hdeltatsq, boxsize_x2, boxsize_y2, boxsize_z2
    real :: elapsedtime
    real :: pavr, sump
    real :: totei, totkini
    integer :: i, counter
    real, external :: timef, read_control, open_output, fcc2, initvel
    real, external :: forceeval, nh1, nh2, elmig, new_nlist, findvac, vacf
    character(len=256) :: controlfile
    integer(2) :: status

    elapsedtime = timef()
    call getarg(1, controlfile, status)
    if (status == -1) then
        write(*,*) "enter controlfilename:"
        read(*,*) controlfile
        if (controlfile == "") stop
    endif

    call read_control(controlfile)
    call open_output(controlfile, iconf, ivac, iforce, imsd, imovie)
    call fcc2(a, int(nbox_x), int(nbox_y), int(nbox_z))

    if (enunit == 0) then
        bolz = 1.380658d-23
        tunit = 1.0d-15
        do i = 1, ntype
            molwt(i) = molwt(i) / avo
        enddo
    endif

    if (enunit == 1) then
        bolz = 8.617317547d-5
        tunit = 0.09822732370761d0
    endif

    ! parameters in cut-off corrections
    eps = 1.507212225d-20
    sigma = 2.5748426d-10
    pset = 2.5d-3
    tau_p = 0.01
    tau_t = 1.0
    beta = 0.744
    m = 1.0d-5
    volumel = 0.0
    counter = 0
    seed = 0.85
    np = n
    deltatsq = (deltat**2) * (tunit**2)
    hdeltatsq = deltatsq * 0.5
    tstp = deltat * tunit

    boxsize_x2 = boxsize_x / 2.0
    boxsize_y2 = boxsize_y / 2.0
    boxsize_z2 = boxsize_z / 2.0

    volume = boxsize_x * boxsize_y * boxsize_z
    bxr = boxsize_x / (volume**0.3333333333)
    byr = boxsize_y / (volume**0.3333333333)
    bzr = boxsize_z / (volume**0.3333333333)

    ! cut-off corrections
    encorr = (2.0 * pi * np * volume * 4.0 * eps) * &
             ((-1.0/9.0) * (sigma**12 / rcut**9) - (-1.0/3.0) * (sigma**6 / rcut**3))

    ! initialize property values
    msd = 0.0
    sumvel = 0.0
    pavr = 0.0

    !! Give initial velocities to atoms
    call initvel(temp)
    call forceeval(np, rcut, boxsize_x, boxsize_y, boxsize_z)

    !! :: nose-hoover thermostat ::
    s = 0.0
    ! initial accelerations, second and third derivatives
    g = 3 * n
    ps = (sumvel - real(g) * bolz * temp) / q

    do i = 1, n
        atoms(i)%ax = atoms(i)%fx / molwt(atoms(i)%kind)
        atoms(i)%ay = atoms(i)%fy / molwt(atoms(i)%kind)
        atoms(i)%az = atoms(i)%fz / molwt(atoms(i)%kind)
        atoms(i)%bx = 0.0
        atoms(i)%by = 0.0
        atoms(i)%bz = 0.0
        atoms(i)%cx = 0.0
        atoms(i)%cy = 0.0
        atoms(i)%cz = 0.0
    end do

    write(5,1) "tottime", "poten", "kinetic", "toten", "temp", "pressure", &
               "hamilton", "volume", "diffusion"
1   format(10x, a7, 8a25)

    !! MD Loop
    ! equilibration
    do tstep = 1, scalestep
        !! Select the appropriate integration routine
        ! call predict(tstp)    !! Uncomment this for Gear 5 value
        ! call verlet1(tstp)    !! Uncomment this for non-reversible Nose-Hoover Thermostat
        ! call verl1(tstp)      !! Uncomment this for traditional velocity Verlet
        ! call anders1()        !! Uncomment this for Andersen pressure control
        call nh1(tstp)          !! Uncomment this for reversible Nose-Hoover Thermostat

        ! pbc
        do i = 1, n
            if (atoms(i)%x < 0.0) then
                atoms(i)%x = atoms(i)%x + boxsize_x
            else if (atoms(i)%x > boxsize_x) then
                atoms(i)%x = atoms(i)%x - boxsize_x
            end if

            !! Uncomment the following lines for PBC in y-direction
            ! if (atoms(i)%y < 0.0) then
            !    atoms(i)%y = atoms(i)%y + boxsize_y
            ! else if (atoms(i)%y > boxsize_y) then
            !    atoms(i)%y = atoms(i)%y - boxsize_y
            ! end if

            if (atoms(i)%z < 0.0) then
                atoms(i)%z = atoms(i)%z + boxsize_z
            else if (atoms(i)%z > boxsize_z) then
                atoms(i)%z = atoms(i)%z - boxsize_z
            end if
        end do

        call forceeval(np, rcut, boxsize_x, boxsize_y, boxsize_z)

        !! Select the appropriate integration routine (corrector step)
        call nh2(tstp)           !! Reversible Nose-Hoover
        ! call anders2()         !! Andersen pressure control
        ! call verlet2(tstp)     !! Non-reversible Nose-Hoover Thermostat
        ! call correct(tstp)     !! Gear 5 value
        ! call berendsen()       !! Berendsen pressure control
        ! call verl2(tstp)       !! Traditional velocity Verlet

        tottime = counter * tstp
        ! call velscale()        !! Temperature control by velocity scaling

        sumvel = 0.0
        do i = 1, n
            sumvel = sumvel + (atoms(i)%vx**2 + atoms(i)%vy**2 + atoms(i)%vz**2) * &
                     molwt(atoms(i)%kind)
        end do

        !! Uncomment the following in a pressure controlled simulation
        ! boxsize_x = volume * bxr
        ! boxsize_y = volume * byr
        ! boxsize_z = volume * bzr

        ! calculate properties
        totkin = 0.5 * sumvel
        counter = counter + 1
        tempdum = (sumvel) / (3.0 * real(n) * bolz)
        pressure = ((n * bolz * tempdum) + (1.0/3.0) * totw) / volume
        sump = sump + pressure
        pavr = sump / counter
        totei = tote / (n + nwall)
        totkini = totkin / (n + nwall)
        hamilton = tote + totkin + (ps**2 * q) / (2.0) + real(g) * tempdum * bolz * s

        write(*,2) tstep, tempdum, tote, totkin, hamilton
2       format(i8, 3es13.5)

        !! update Linked list
        call new_nlist()

        !! Write the properties to a file
        if (enunit == 1) then
            write(5,15) tottime, totei, totkini, totei + totkini + encorr, &
                        tempdum, pressure, hamilton, volume
15          format(1x, 8es25.17)
        endif

        if (enunit == 0) then
            write(5,16) tottime, totei, totkini, totei + totkini + encorr, &
                        tempdum, pressure, hamilton, volume
16          format(1x, 8es25.17)
        endif

        !! write vacancy positions at intervals
        if (ivac > 0) then
            if (mod(tstep, ivac) == 0) then
                write(8,*) nvac
                write(8,*) "vacancy positions:", "at step=", tstep, " time=", tottime, "fs"
                do i = n + nwall + 1, n + nwall + nvac
                    write(8,14) atoms(i)%x, atoms(i)%y, atoms(i)%z
14                  format(1x, 3f20.10)
                enddo
            endif
        endif

        !! Write the force on the atoms at intervals
        if (iforce > 0) then
            if (mod(tstep, iforce) == 0) then
                write(7,*) "forces", "at step=", tstep, " time=", tottime, "fs"
                do i = 1, n + nwall + nvac
                    write(7,13) atoms(i)%fx, atoms(i)%fy, atoms(i)%fz, &
                                atoms(i)%elfx, atoms(i)%elfy, atoms(i)%elfz
13                  format(1x, 6es15.7)
                enddo
            endif
        endif

        !! write the trajectory of atoms to .xyz file
        if (imovie > 0) then
            if (mod(tstep, imovie) == 0) then
                write(15,*) n + nwall
                write(15,*) 'md step', tstep
                do i = 1, n + nwall
                    write(15, '(i3, 3f10.4)') atoms(i)%kind, atoms(i)%x, atoms(i)%y, atoms(i)%z
                enddo
            endif
        endif
    end do

    ! production
    ! store initial positions for MSD calculations
    do i = 1, n
        atoms(i)%x0 = atoms(i)%x
        atoms(i)%y0 = atoms(i)%y
        atoms(i)%z0 = atoms(i)%z
    end do

    do tstep = scalestep + 1, totstep
        ! store first positions for msd
        do i = 1, n
            atoms(i)%xb = atoms(i)%x
            atoms(i)%yb = atoms(i)%y
            atoms(i)%zb = atoms(i)%z
        enddo

        call nh1(tstp)
        ! call predict(tstp)
        ! call verlet1(tstp)
        ! call verl1(tstp)

        ! pbc
        do i = 1, n
            if (atoms(i)%x < 0.0) then
                atoms(i)%x = atoms(i)%x + boxsize_x
            else if (atoms(i)%x > boxsize_x) then
                atoms(i)%x = atoms(i)%x - boxsize_x
            end if

            ! if (atoms(i)%y < 0.0) then
            !    atoms(i)%y = atoms(i)%y + boxsize_y
            ! else if (atoms(i)%y > boxsize_y) then
            !    atoms(i)%y = atoms(i)%y - boxsize_y
            ! end if

            if (atoms(i)%z < 0.0) then
                atoms(i)%z = atoms(i)%z + boxsize_z
            else if (atoms(i)%z > boxsize_z) then
                atoms(i)%z = atoms(i)%z - boxsize_z
            end if
        end do

        !! Evaluate Forces and calculate energy
        call forceeval(np, rcut, boxsize_x, boxsize_y, boxsize_z)

        if (mod(tstep, ielmigstep) == 0) then
            ! Find vacancies
            call findvac()
            ! add electromigration force
            call elmig()
            ! EM force to MD force
        end if

        do i = 1, n
            atoms(i)%fx = atoms(i)%fx + atoms(i)%elfx
            atoms(i)%fy = atoms(i)%fy + atoms(i)%elfy
            atoms(i)%fz = atoms(i)%fz + atoms(i)%elfz
        end do

        call nh2(tstp)
        ! call verlet2(tstp)
        ! call correct(tstp)
        ! call verl2(tstp)

        tottime = counter * tstp

        if (imsd /= 0) then
            if (mod(tstep, imsd) == 0) then
                msd = 0.0
                ! rewind pbc for msd calculation
                do i = 1, n
                    atoms(i)%xp = atoms(i)%x + nint((atoms(i)%xb - atoms(i)%x) / boxsize_x) * boxsize_x
                    atoms(i)%yp = atoms(i)%y + nint((atoms(i)%yb - atoms(i)%y) / boxsize_y) * boxsize_y
                    atoms(i)%zp = atoms(i)%z + nint((atoms(i)%zb - atoms(i)%z) / boxsize_z) * boxsize_z
                enddo

                ! MSD calculation
                do i = msdstart, msdend
                    rp = sqrt(atoms(i)%xp**2 + atoms(i)%yp**2 + atoms(i)%zp**2)
                    r0 = sqrt(atoms(i)%x0**2 + atoms(i)%y0**2 + atoms(i)%z0**2)
                    msd = msd + (rp - r0)**2
                enddo
                msd = msd / real(msdend - msdstart + 1)
                write(4, '(f20.13, 1x, es25.17)') real((tstep - (scalestep + imsd))) * tstp, msd
            endif
        endif

        sumvel = 0.0
        do i = 1, n
            sumvel = sumvel + (atoms(i)%vx**2 + atoms(i)%vy**2 + atoms(i)%vz**2) * molwt(atoms(i)%kind)
        end do

        ! calculate properties
        totkin = 0.5 * sumvel
        counter = counter + 1
        tempdum = (sumvel) / (3.0 * real(n) * bolz)
        pressure = ((n * bolz * tempdum) + (1.0/3.0) * totw) / volume
        sump = sump + pressure
        pavr = sump / counter
        totei = tote / (n + nwall)
        totkini = totkin / (n + nwall)
        hamilton = tote + totkin + (ps**2 * q) / (2.0) + real(g) * tempdum * bolz * s

        write(*,3) tstep, tempdum, tote, hamilton
3       format(i8, 3es13.5)

        ! Update neighbor list
        call new_nlist()

        if (enunit == 1) then
            write(5,18) tottime, totei, totkini, totei + totkini + encorr, &
                        tempdum, pressure, hamilton, volume, msd
18          format(1x, 9es25.17)
        endif

        if (enunit == 0) then
            write(5,19) tottime, totei, totkini, totei + totkini + encorr, &
                        tempdum, pressure, hamilton, volume, msd
19          format(1x, 9es25.17)
        endif

        if (iconf > 0) then
            if (mod(tstep, iconf) == 0) then
                write(9,*) tstep
                do i = 1, n
                    write(9,21) atoms(i)%x, atoms(i)%y, atoms(i)%z, &
                                atoms(i)%vx, atoms(i)%vy, atoms(i)%vz
21                  format(1x, 6es15.7)
                end do
            endif
        endif

        if (iforce > 0) then
            if (mod(tstep, iforce) == 0) then
                write(7,*) "forces", "at step=", tstep, "time", tottime, "fs"
                do i = 1, n + nvac + nwall
                    write(7,24) atoms(i)%fx, atoms(i)%fy, atoms(i)%fz, &
                                atoms(i)%elfx, atoms(i)%elfy, atoms(i)%elfz
24                  format(1x, 6es15.7)
                enddo
            endif
        endif

        if (ivac > 0) then
            if (mod(tstep, ivac) == 0) then
                write(8,*) nvac
                write(8,*) "vacancy positions:", "at step=", tstep, " time=", tottime, "fs"
                do i = n + nwall + 1, n + nvac + nwall
                    write(8,22) atoms(i)%x, atoms(i)%y, atoms(i)%z
22                  format(1x, 3f20.10)
                enddo
            endif
        endif

        if (imovie > 0) then
            if (mod(tstep, imovie) == 0) then
                write(15,*) n + nwall
                write(15,*) 'md step', tstep
                do i = 1, n + nwall
                    write(15, '(i3, 3f10.4)') atoms(i)%kind, atoms(i)%x, atoms(i)%y, atoms(i)%z
                enddo
            endif
        endif
    enddo

    if (iconf > 0) close(9)
    if (imsd > 0) close(4)
    if (iforce > 0) close(7)
    if (ivac > 0) close(8)
    if (imovie > 0) close(15)

    !! Calculate velocity autocorrelation function
    call vacf(controlfile, totstep, scalestep, tstp, n, iconf)

    elapsedtime = timef()
    write(*,*) "total elapsed time is:", elapsedtime
    write(5,*) "total elapsed time is:", elapsedtime, "secs."
    close(5)

end program md

!! This routine opens output files
subroutine open_output(name1, ic, iv, ifrc, ims, imov)
    implicit none
    integer :: ic, iv, ifrc, ims, imov
    character(*) :: name1
    character(len=len_trim(name1)) :: name2

    name2 = trim(name1)
    open(unit=5, file=name2//".out")
    if (ic > 0) open(unit=9, file=name2//".con")
    if (iv > 0) open(unit=8, file=name2//".vac")
    if (ifrc > 0) open(unit=7, file=name2//".elf")
    if (ims > 0) open(unit=4, file=name2//".msd")
    if (imov > 0) open(unit=15, file=name2//".xyz")
end subroutine open_output

!! Read input data
subroutine read_control(name1)
    use vardecs
    use commondefs
    implicit none

    integer :: i
    character(*) :: name1
    character(len=len_trim(name1)) :: name2
    integer :: npot, pottype, at1, at2, bit

    name2 = trim(name1)
    open(unit=10, file=name2)
    open(unit=11, file='fcc.txt')

    read(10,*) n, nvac, nwall
    read(10,*) ntype
    read(10,*) totstep, scalestep, scaleint
    read(10,*) deltat
    read(10,*) tol
    read(10,*) msdstart, msdend
    read(10,*) temp
    read(10,*) q
    read(10,*) enunit
    read(10,*) a, boxsize_x, boxsize_y, boxsize_z, rcut
    read(10,*) iconf, ivac, iforce, imsd, imovie, ielmigstep
    read(10,*) (molwt(i), i=1, ntype)

    do i = 1, ntype
        read(10,*) atno(i), conc(i)
    enddo

    read(10,*) vdd(1), vdd(2), vdd(3)

    allocate(atoms(n + nvac + nwall + 150))
    allocate(ions(n + nvac + nwall + 150))

    do i = 1, 50
        pairtype(i)%pot = 0
    enddo

    read(10,*) npot
    do i = 1, npot
        read(10,*) pottype
        select case (pottype)
        case (1) ! LJ potential
            read(10,*) at1, at2, ljpot(i)%eps, ljpot(i)%sigma
            bit = 2**at1 + 2**at2
            pairtype(bit)%index = i
            pairtype(bit)%pot = 1
        case (2) ! fourier potential
            read(10,*) at1, at2, fourierpot(i)%a, fourierpot(i)%b, &
                       fourierpot(i)%c, fourierpot(i)%d, fourierpot(i)%e, fourierpot(i)%f, &
                       fourierpot(i)%g, fourierpot(i)%h, fourierpot(i)%ix, fourierpot(i)%jx, &
                       fourierpot(i)%kx, fourierpot(i)%lx, fourierpot(i)%mx, &
                       fourierpot(i)%rmin, fourierpot(i)%rmax, fourierpot(i)%eps, fourierpot(i)%sigma
            bit = 2**at1 + 2**at2
            pairtype(bit)%index = i
            pairtype(bit)%pot = 2
        end select
    enddo

    nbox_x = int(boxsize_x / a)
    nbox_y = int(boxsize_y / a)
    nbox_z = int(boxsize_z / a)
    nperf = nbox_x * nbox_y * nbox_z * 4

    ! read atomic coordinates
    do i = 1, n + nvac + nwall
        read(10,10) atoms(i)%kind, atoms(i)%x, atoms(i)%y, atoms(i)%z
10      format(i2, 1x, f20.17, 4x, f20.17, 4x, f20.17)
    end do

    allocate(xini(nperf), yini(nperf), zini(nperf))
    close(11)

    ncellx = int(boxsize_x / rcut / 1.01)
    ncelly = int(boxsize_y / rcut / 1.01)
    ncellz = int(boxsize_z / rcut / 1.01)
    cell_factorx = ncellx / boxsize_x
    cell_factory = ncelly / boxsize_y
    cell_factorz = ncellz / boxsize_z

    call new_nlist()
    call neighbour_cells()
end subroutine read_control

!! This routine gives initial velocities to atoms
subroutine initvel(temper)
    use commondefs
    use vardecs
    implicit none

    real :: rtemp, sumx, sumy, sumz
    real :: gauss, dummy, temper, scale, sumvel2
    real :: aheat2
    integer :: i

    dummy = 0.5
    rtemp = sqrt(temper)
    aheat2 = real(n) * 3.0 * bolz * temper

    do i = 1, n
        atoms(i)%vx = rtemp * gauss(dummy)
        atoms(i)%vy = rtemp * gauss(dummy)
        atoms(i)%vz = rtemp * gauss(dummy)
    end do

    sumx = 0.0
    sumy = 0.0
    sumz = 0.0
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
    do i = 1, np
        sumvel2 = sumvel2 + (atoms(i)%vx**2 + atoms(i)%vy**2 + atoms(i)%vz**2) * molwt(atoms(i)%kind)
    end do

    scale = sqrt(aheat2 / sumvel2)

    sumvel = 0.0
    do i = 1, np
        atoms(i)%vx = atoms(i)%vx * scale
        atoms(i)%vy = atoms(i)%vy * scale
        atoms(i)%vz = atoms(i)%vz * scale
        sumvel = sumvel + (atoms(i)%vx**2 + atoms(i)%vy**2 + atoms(i)%vz**2) * molwt(atoms(i)%kind)
    end do

    return
end subroutine initvel

!! This function gives a random number in Gaussian distribution
real function gauss(dummy)
    implicit none
    real :: a1, a3, a5, a7, a9
    parameter (a1 = 3.949846138, a3 = 0.252408784)
    parameter (a5 = 0.076542912, a7 = 0.008355968)
    parameter (a9 = 0.029899776)
    real :: sum, r, r2
    real :: ranf, dummy
    integer :: i

    sum = 0.0
    do i = 1, 12
        sum = sum + ranf(dummy)
    enddo

    r = (sum - 6.0) / 4.0
    r2 = r * r
    gauss = ((((a9 * r2 + a7) * r2 + a5) * r2 + a3) * r2 + a1) * r
    return
end function gauss

!! This function generates a random number
real function ranf(dummy)
    implicit none
    integer :: l, c, m
    parameter (l = 1029, c = 221591, m = 1048576)
    integer, save :: seed = 0
    real :: dummy

    seed = mod(seed * l + c, m)
    ranf = real(seed) / m
    return
end function ranf

!! Reversible Nose-Hoover thermostat predictor step
subroutine nh1(dt)
    use vardecs
    use commondefs
    implicit none

    real :: dt, dtsq2, dt2, sumvel2
    integer :: i

    dt2 = dt / 2.0
    dtsq2 = dt * dt2

    sumvel2 = 0.0
    do i = 1, n
        sumvel2 = sumvel2 + (atoms(i)%vx**2 + atoms(i)%vy**2 + atoms(i)%vz**2) * molwt(atoms(i)%kind)
    end do

    ps = ps + (sumvel2 - real(g) * bolz * temp) * dt2 / q
    s = s + ps * dt

    do i = 1, n
        atoms(i)%vx = atoms(i)%vx * dexp(-ps * dt2) + (atoms(i)%fx / molwt(atoms(i)%kind)) * dt2
        atoms(i)%vy = atoms(i)%vy * dexp(-ps * dt2) + (atoms(i)%fy / molwt(atoms(i)%kind)) * dt2
        atoms(i)%vz = atoms(i)%vz * dexp(-ps * dt2) + (atoms(i)%fz / molwt(atoms(i)%kind)) * dt2

        atoms(i)%x = atoms(i)%x + dt * atoms(i)%vx
        atoms(i)%y = atoms(i)%y + dt * atoms(i)%vy
        atoms(i)%z = atoms(i)%z + dt * atoms(i)%vz
    end do
end subroutine nh1

!! Reversible Nose-Hoover thermostat corrector step
subroutine nh2(dt)
    use vardecs
    use commondefs
    implicit none

    real :: dt, dt2, sumvel2
    integer :: i

    dt2 = dt / 2.0
    sumvel2 = 0.0
    totkin = 0.0

    do i = 1, n
        atoms(i)%vx = (atoms(i)%vx + (atoms(i)%fx / molwt(atoms(i)%kind)) * dt2) * dexp(-ps * dt2)
        atoms(i)%vy = (atoms(i)%vy + (atoms(i)%fy / molwt(atoms(i)%kind)) * dt2) * dexp(-ps * dt2)
        atoms(i)%vz = (atoms(i)%vz + (atoms(i)%fz / molwt(atoms(i)%kind)) * dt2) * dexp(-ps * dt2)

        sumvel2 = sumvel2 + (atoms(i)%vx**2 + atoms(i)%vy**2 + atoms(i)%vz**2) * molwt(atoms(i)%kind)
    end do

    ps = ps + (sumvel2 - real(g) * bolz * temp) * dt2 / q
end subroutine nh2

!! This routine calculates the forces from the interatomic potential defined on atoms in MD
subroutine forceeval(npart, rcut2, box_x, box_y, box_z)
    use vardecs
    use commondefs
    implicit none

    real :: rcut2, box_x, box_y, box_z, rcutsq2, box_x2, box_y2, box_z2, xr, yr, zr, xi, yi, zi
    real, external :: forcefunc, potfunc, lennardforce, lennardenergy
    integer, external :: find_cell
    real :: rxij, ryij, rzij, rijsq, rijsqrt
    real :: fij, rpl, fxij, fyij, fzij
    integer :: npart, bit, i, j, ic, inn, in

    tote = 0.0
    totw = 0.0
    rcutsq2 = rcut2**2
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

        ic = find_cell(atoms(i)%x, atoms(i)%y, atoms(i)%z)

        do inn = int(ic * 27 + 1), int(ic * 27 + 27)
            in = n_v_list(inn)
            j = head_chain(in)

            do while (j /= 0)
                if ((j /= i) .and. (j >= i + 1)) then
                    bit = 2**atoms(i)%kind + 2**atoms(j)%kind
                    if (pairtype(bit)%pot == 0) then
                        j = linked_list(j)
                        cycle
                    endif

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

                    rijsq = (rxij * rxij) + (ryij * ryij) + (rzij * rzij)
                    rijsqrt = sqrt(rijsq)

                    if (rijsq < rcutsq2) then
                        select case (pairtype(bit)%pot)
                        case (1)
                            fij = lennardforce(rijsqrt, bit)
                            tote = tote + lennardenergy(rijsqrt, bit)
                        case (2)
                            fij = forcefunc(rijsqrt, bit)
                            tote = tote + potfunc(rijsqrt, bit)
                        end select

                        rpl = fij
                        fxij = fij * rxij / rijsqrt
                        fyij = fij * ryij / rijsqrt
                        fzij = fij * rzij / rijsqrt

                        atoms(i)%fx = atoms(i)%fx + fxij
                        atoms(i)%fy = atoms(i)%fy + fyij
                        atoms(i)%fz = atoms(i)%fz + fzij

                        atoms(j)%fx = atoms(j)%fx - fxij
                        atoms(j)%fy = atoms(j)%fy - fyij
                        atoms(j)%fz = atoms(j)%fz - fzij

                        totw = totw + rpl
                    end if
                endif
                j = linked_list(j)
            enddo
        end do
    end do
    return
end subroutine forceeval

!! This function calculates potential energy by Fourier series
function potfunc(rij, bit)
    use vardecs
    implicit none

    real :: potfunc, rij, rx
    integer :: bit

    rx = pi * (rij - fourierpot(pairtype(bit)%index)%rmin) / &
         (fourierpot(pairtype(bit)%index)%rmax - fourierpot(pairtype(bit)%index)%rmin)

    potfunc = fourierpot(pairtype(bit)%index)%a + &
              fourierpot(pairtype(bit)%index)%b * cos(rx) + &
              fourierpot(pairtype(bit)%index)%c * sin(rx) + &
              fourierpot(pairtype(bit)%index)%d * cos(2.0 * rx) + &
              fourierpot(pairtype(bit)%index)%e * sin(2.0 * rx) + &
              fourierpot(pairtype(bit)%index)%f * cos(3.0 * rx) + &
              fourierpot(pairtype(bit)%index)%g * sin(3.0 * rx) + &
              fourierpot(pairtype(bit)%index)%h * cos(4.0 * rx) + &
              fourierpot(pairtype(bit)%index)%ix * sin(4.0 * rx) + &
              fourierpot(pairtype(bit)%index)%jx * cos(5.0 * rx) + &
              fourierpot(pairtype(bit)%index)%kx * sin(5.0 * rx) + &
              fourierpot(pairtype(bit)%index)%lx * cos(6.0 * rx) + &
              fourierpot(pairtype(bit)%index)%mx * sin(6.0 * rx)
    return
end function potfunc

!! This function calculates the force by Fourier series
function forcefunc(rij2, bit)
    use vardecs
    implicit none

    real :: forcefunc, rij2, ry, y
    integer :: bit

    ry = pi * (rij2 - fourierpot(pairtype(bit)%index)%rmin) / &
         (fourierpot(pairtype(bit)%index)%rmax - fourierpot(pairtype(bit)%index)%rmin)

    y = pi / (fourierpot(pairtype(bit)%index)%rmax - fourierpot(pairtype(bit)%index)%rmin)

    forcefunc = fourierpot(pairtype(bit)%index)%b * sin(ry) * y - &
                fourierpot(pairtype(bit)%index)%c * cos(ry) * y - &
                2.0 * fourierpot(pairtype(bit)%index)%d * sin(2.0 * ry) * y - &
                2.0 * fourierpot(pairtype(bit)%index)%e * cos(2.0 * ry) * y - &
                3.0 * fourierpot(pairtype(bit)%index)%f * sin(3.0 * ry) * y - &
                3.0 * fourierpot(pairtype(bit)%index)%g * cos(3.0 * ry) * y - &
                4.0 * fourierpot(pairtype(bit)%index)%h * sin(4.0 * ry) * y - &
                4.0 * fourierpot(pairtype(bit)%index)%ix * cos(4.0 * ry) * y + &
                5.0 * fourierpot(pairtype(bit)%index)%jx * sin(5.0 * ry) * y - &
                5.0 * fourierpot(pairtype(bit)%index)%kx * cos(5.0 * ry) * y + &
                6.0 * fourierpot(pairtype(bit)%index)%lx * sin(6.0 * ry) * y - &
                6.0 * fourierpot(pairtype(bit)%index)%mx * cos(6.0 * ry) * y
    return
end function forcefunc

!! Calculate forces if Lennard-Jones potential is used
function lennardforce(rij3, bit)
    use vardecs
    implicit none

    integer :: bit
    real :: lennardforce, rij3

    lennardforce = 24.0 * ljpot(pairtype(bit)%index)%eps * &
                   (2.0 * (ljpot(pairtype(bit)%index)%sigma**12 / rij3**13) - &
                          (ljpot(pairtype(bit)%index)%sigma**6 / rij3**7))
    return
end function lennardforce

!! Potential energy if Lennard-Jones potential is used
function lennardenergy(rij, bit)
    use vardecs
    implicit none

    integer :: bit
    real :: lennardenergy, rij

    lennardenergy = 4.0 * ljpot(pairtype(bit)%index)%eps * &
                    ((ljpot(pairtype(bit)%index)%sigma / rij)**12 - &
                     (ljpot(pairtype(bit)%index)%sigma / rij)**6)
    return
end function lennardenergy

!! This routine finds the cell of atom in linked list
integer function find_cell(ionx, iony, ionz)
    use commondefs
    implicit none

    real :: ionx, iony, ionz
    integer :: ix, iy, iz

    ix = int(ionx * cell_factorx)
    iy = int(iony * cell_factory)
    iz = int(ionz * cell_factorz)

    find_cell = ix + iy * ncellx + iz * ncellx * ncelly
end function find_cell

!! Linked list cell generator
subroutine neighbour_cells()
    use commondefs
    implicit none

    integer :: ix, iy, iz, ic, in, itel, icx, icy, icz, iccx, iccy, iccz

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
            endif

            do iy = -1, 1
                iccy = icy + iy
                if (iccy < 0) then
                    iccy = iccy + ncelly
                else if (iccy >= ncelly) then
                    iccy = iccy - ncelly
                endif

                do ix = -1, 1
                    iccx = icx + ix
                    if (iccx < 0) then
                        iccx = iccx + ncellx
                    else if (iccx >= ncellx) then
                        iccx = iccx - ncellx
                    endif

                    in = iccx + iccy * ncellx + iccz * ncellx * ncelly
                    itel = itel + 1
                    n_v_list(int(ic * 27 + itel)) = in
                enddo
            enddo
        enddo
    enddo
end subroutine neighbour_cells

!! Construct linked list
subroutine new_nlist()
    use commondefs
    implicit none

    integer :: i, ic
    integer, external :: find_cell

    do ic = 0, ncellx * ncelly * ncellz - 1
        head_chain(ic) = 0
    enddo

    do i = 1, n + nwall
        ic = find_cell(atoms(i)%x, atoms(i)%y, atoms(i)%z)
        linked_list(i) = head_chain(ic)
        head_chain(ic) = i
    enddo
end subroutine new_nlist

!! Initialization part of the electromigration program
subroutine init()
    use defs
    use commondefs
    implicit none

    integer :: i, j

    vd = sqrt(vdd(1)**2 + vdd(2)**2 + vdd(3)**2)
    vddir(1) = vdd(1) / vd
    vddir(2) = vdd(2) / vd
    vddir(3) = vdd(3) / vd

    do i = 1, ntype
        do j = 1, 56
            if (atno(i) == el(j)%atno) sp(i) = el(j)
        enddo
    enddo

    box_x = boxsize_x / 0.529177
    box_y = boxsize_y / 0.529177
    box_z = boxsize_z / 0.529177
    cutoff = rcut / 0.529177

    do i = 1, n + nvac + nwall
        ions(i)%x = atoms(i)%x / 0.529177249
        ions(i)%y = atoms(i)%y / 0.529177249
        ions(i)%z = atoms(i)%z / 0.529177249
        ions(i)%kind = atoms(i)%kind
    enddo

    close(6)

    omave = 0.0
    zave = 0.0
    mave = 0.0

    do i = 1, ntype
        omave = omave + conc(i) * sp(i)%om
        zave = zave + conc(i) * sp(i)%z
        mave = mave + conc(i) * sp(i)%m
    enddo

    mave = 1.0
    kf = (3.0 * pi * pi * zave / omave)**0.33333333333333
    vf = kf / mave
    ef = kf * kf / 2.0 / mave
    rs = mave * (3.0 / 4.0 / pi / zave * omave)**0.33333333333333
    ksi = 2.0 * (1.0 + 0.02544716775 * rs)

    do i = 1, n + nwall
        atoms(i)%elfx = 0.0
        atoms(i)%elfy = 0.0
        atoms(i)%elfz = 0.0
    enddo
end subroutine init

!! Cosine of angle between vd and Rj-Rm
real function coss(rx, ry, rz)
    use defs
    implicit none

    real :: rx, ry, rz
    coss = rx * rr * vddir(1) + ry * rr * vddir(2) + rz * rr * vddir(3)
end function coss

!! Dot product of vd and Rj
real function dotproduct(ii, var)
    use defs
    use commondefs
    implicit none

    real :: var
    integer :: ii

    if (ii == 1) dotproduct = vdd(1) * var + vdd(2) * y + vdd(3) * z
    if (ii == 2) dotproduct = vdd(1) * x + vdd(2) * var + vdd(3) * z
    if (ii == 3) dotproduct = vdd(1) * x + vdd(2) * y + vdd(3) * var
end function dotproduct

!! Spherical Bessel function
real function bessel(q)
    use defs
    implicit none

    real :: q, xx
    xx = q * rr
    bessel = sin(xx) / (xx * xx) - cos(xx) / xx
end function bessel

!! Functions for integration in force calculation
real function integrant1(q)
    use defs
    use commondefs
    implicit none

    real :: wi, wj, q
    real, external :: w, bessel

    wj = w(sp(ions(force_j)%kind), q)
    if (ions(force_i)%kind < 0) then
        wi = w(sp(1), q)
    else
        wi = w(sp(ions(force_i)%kind), q)
    endif
    integrant1 = wj * wi * q * q * q
end function integrant1

real function integrant2(q)
    use defs
    use commondefs
    implicit none

    real :: wi, wj, q
    real, external :: w, bessel

    wj = w(sp(ions(force_j)%kind), q)
    wi = w(sp(1), q)
    integrant2 = wj * (wj - wi) * q * q * q * q * bessel(q)
end function integrant2

real function integrant3(q)
    use defs
    use commondefs
    implicit none

    real :: wi, wj, q
    real, external :: w, bessel

    wi = w(sp(1), q)
    wj = w(sp(ions(force_j)%kind), q)
    integrant3 = wi * wj * q * q * q * bessel(q)
end function integrant3

!! Gauss quadrature integration
real function integrate(integrant)
    use defs
    implicit none

    real :: a, b, epsabs, epsrel, result, abserr
    integer :: key, neval, ier
    real, external :: integrant
    external :: qag

    a = 0.00001
    b = 2.0 * kf
    epsabs = 0.0
    epsrel = 0.001
    key = 3

    call qag(integrant, a, b, epsabs, epsrel, key, result, abserr, neval, ier)
    integrate = result
end function integrate

!! Calculate potential used in EM force calculation
real function uu_s(ii, var)
    use defs
    use commondefs
    implicit none

    real, external :: coss, dotproduct, integrant1, integrant2, integrant3, integrate
    real :: t1, term1, term2, term3, var, rx, ry, rz
    integer :: ii

    t1 = omave**2 * mave * kf / 4.0 / pi**3
    uu_s = 0.0

    do force_i = 1, n + nvac + nwall
        if (ions(force_i)%kind /= 1) then
            rx = ions(force_i)%x - x
            ry = ions(force_i)%y - y
            rz = ions(force_i)%z - z

            if (ii == 1) rx = ions(force_i)%x - var
            if (ii == 2) ry = ions(force_i)%y - var
            if (ii == 3) rz = ions(force_i)%z - var

            xr = abs(rx)
            yr = abs(ry)
            zr = abs(rz)

            if (xr > box_x / 2.0) rx = (xr - box_x) * sign(1.0, rx)
            if (yr > box_y / 2.0) ry = (yr - box_y) * sign(1.0, ry)
            if (zr > box_z / 2.0) rz = (zr - box_z) * sign(1.0, rz)

            rr = sqrt(rx*rx + ry*ry + rz*rz)

            if (rr < cutoff) then
                if (force_i /= force_j) then
                    term1 = t1 * vd / vf * coss(-rx, ry, rz) * integrate(integrant1)
                    if (ions(force_i)%kind < 0) then
                        uu_s = uu_s - term1
                    else
                        uu_s = uu_s + term1
                    endif

                    if (ions(force_i)%kind > 1) then
                        term3 = t1 * vd / vf * coss(-rx, ry, rz) * integrate(integrant3)
                        uu_s = uu_s - term3
                    endif
                else
                    term2 = t1 / 3.0 * dotproduct(ii, var) / vf * integrate(integrant2)
                    uu_s = uu_s + term2
                endif
            endif
        endif
    enddo
end function uu_s

!! Main EM force routine
subroutine elmig()
    use defs
    use commondefs
    implicit none

    real :: error
    real, external :: uu_s, dfridr
    include 'data.h'

    !! All units in EM force are in au.
    call init()

    do force_j = 1, n
        x = atoms(force_j)%x
        y = atoms(force_j)%y
        z = atoms(force_j)%z

        ! convert to eV/A
        atoms(force_j)%elfx = -dfridr(1, x, 0.2, error) * 51.42210754
        atoms(force_j)%elfy = -dfridr(2, y, 0.2, error) * 51.42210754
        atoms(force_j)%elfz = -dfridr(3, z, 0.2, error) * 51.42210754
    enddo
end subroutine elmig

!! Heine-Abarenkov-Aniamalu psp N(q)
real function bofq(ele, q)
    use defs
    implicit none

    real :: q, qrm, qrc, t1, t2, t3, t4, t5, t6, t7
    type(elements) :: ele

    qrm = q * ele%rm
    qrc = q * ele%rc

    t1 = pi / omave / q / q
    t2 = 8.0 * t1 * ele%c / q
    t3 = sin(qrm) - qrm * cos(qrm)
    t4 = 8.0 * t1 * ele%z * cos(qrm)
    t5 = 4.0 * t1 * ele%ec / q
    t6 = 24.0 * t1 * ele%z * ele%alf * qrc / qrc / qrc
    t7 = sin(qrc) - qrc * cos(qrc)

    bofq = -t2 * t3 - t4 + (t5 - t6) * t7
end function bofq

!! Heine-Abarenkov-Aniamalu psp N(q)
real function nofq(ele, q)
    use defs
    implicit none

    real :: q, t1, t2, t3, t4, xx, yy, cost, costp
    real :: j0x, j1x, j2x, j3x, j0y, j1y, j2y, j3y
    type(elements) :: ele

    t1 = pi * ele%m * ele%m / omave
    xx = kf * ele%rm

    j0x = sin(xx) / xx
    j1x = sin(xx) / xx / xx - cos(xx) / xx
    j2x = 3.0 / xx * j1x - j0x
    j3x = 5.0 / xx * j2x - j1x

    if ((q - 2.0 * kf) <= 0.0) then
        cost = (1.0 - q * q / 2.0 / kf / kf)
        t2 = -4.0 * t1 * (ele%a0 - ele%c) * (j0x * j0x - j1x * cos(xx) / xx)
        t3 = -12.0 * t1 * (ele%a1 - ele%c) * (j1x * j1x - j0x * j2x) * cost
        t4 = -20.0 * t1 * (ele%a2 - ele%c) * (j2x * j2x - j1x * j3x) * (3.0 * cost * cost - 1.0) / 2.0
    else
        yy = (q - kf) * ele%rm
        j0y = sin(yy) / yy
        j1y = sin(yy) / yy / yy - cos(yy) / yy
        j2y = 3.0 / yy * j1y - j0y
        j3y = 5.0 / yy * j2y - j1y

        costp = (xx * xx + yy * yy - (q * ele%rm)**2) / 2.0 / xx / yy
        t2 = -8.0 * t1 * (ele%a0 - ele%c) / (xx * xx - yy * yy) * (xx * j1x * j0y - yy * j1y * j0x)
        t3 = -24.0 * t1 * (ele%a1 - ele%c) / (xx * xx - yy * yy) * (xx * j2x * j1y - yy * j2y * j1x) * costp
        t4 = -40.0 * t1 * (ele%a2 - ele%c) / (xx * xx - yy * yy) * (xx * j3x * j2y - yy * j3y * j2x) * (3.0 * costp * costp - 1.0) / 2.0
    end if

    nofq = t2 + t3 + t4
    return
end function nofq

!! Unscreened HAA form factor
real function w0(ele, q)
    use defs
    implicit none

    real :: q
    type(elements) :: ele
    real, external :: nofq, bofq

    w0 = (nofq(ele, q) + bofq(ele, q)) / 2.0
end function w0

!! Screened HAA form factor
real function w(ele, q)
    use defs
    implicit none

    real :: q, qq, ggvosko, lindhard, vc, hartree
    real, external :: w0
    type(elements) :: ele

    qq = q / 2.0 / kf
    ggvosko = 0.5 * (q**2 / (q**2 + ksi * kf**2))
    lindhard = mave * kf / pi**2 * (0.5 + (1.0 - qq**2) / 4.0 / qq * log(abs((1.0 + qq) / (1.0 - qq))))
    vc = 4.0 * pi / q**2
    hartree = 1.0 - vc * lindhard * (1.0 + vc * lindhard * ggvosko)

    w = w0(ele, q) / hartree
end function w

!! Find vacancies in the md cell
subroutine findvac()
    use commondefs
    use vardecs
    implicit none

    integer :: dum, i, j
    real :: rxij, ryij, rzij
    real, external :: fcc
    logical :: flag1

    call fcc(a, int(nbox_x), int(nbox_y), int(nbox_z))
    open(11, file="fcc.txt")

    do i = 1, nperf
        read(11,11) xini(i), yini(i), zini(i)
11      format(f20.17, 4x, f20.17, 4x, f20.17)
    end do
    close(11)

    dum = 1
    flag1 = .false.

    do j = 1, nperf
        do i = 1, n + nwall
            rxij = xini(j) - atoms(i)%x
            ryij = yini(j) - atoms(i)%y
            rzij = zini(j) - atoms(i)%z

            if (abs(rxij) > tol) rxij = (abs(rxij) - boxsize_x) * sign(1.0, rxij)
            if (abs(ryij) > tol) ryij = (abs(ryij) - boxsize_y) * sign(1.0, ryij)
            if (abs(rzij) > tol) rzij = (abs(rzij) - boxsize_z) * sign(1.0, rzij)

            if (sqrt(rxij**2 + ryij**2 + rzij**2) < (tol * sqrt(3.0))) flag1 = .true.
        enddo

        if (.not. flag1) then
            atoms(n + nwall + dum)%x = xini(j)
            atoms(n + nwall + dum)%y = yini(j)
            atoms(n + nwall + dum)%z = zini(j)
            atoms(n + nwall + dum)%kind = -1
            dum = dum + 1
        endif
        flag1 = .false.
    enddo
    nvac = dum - 1
end subroutine findvac

!! Construct FCC cell from old data used in findvac
subroutine fcc(a, nx, ny, nz)
    implicit none
    real :: a
    real, dimension(4) :: x0, y0, z0
    integer :: ix, iy, iz, i, nx, ny, nz
    real, allocatable :: x(:), y(:), z(:)

    allocate(x(nx * ny * nz * 4), y(nx * ny * nz * 4), z(nx * ny * nz * 4))
    open(11, file="fcc.txt", position="rewind", status="old")

    x0(1) = 0.0; y0(1) = 0.0; z0(1) = 0.0
    x0(2) = 0.0; y0(2) = 0.5 * a; z0(2) = 0.5 * a
    x0(3) = 0.5 * a; y0(3) = 0.0; z0(3) = 0.5 * a
    x0(4) = 0.5 * a; y0(4) = 0.5 * a; z0(4) = 0.0

    do iy = 0, ny - 1
        do ix = 0, nx - 1
            do iz = 0, nz - 1
                do i = 1, 4
                    x(i) = ix * a + x0(i)
                    y(i) = iy * a + y0(i)
                    z(i) = iz * a + z0(i)
                    write(11,10) x(i), y(i), z(i)
10                  format(f20.17, 4x, f20.17, 4x, f20.17)
                end do
            end do
        end do
    end do
    close(11)
end subroutine fcc

subroutine fcc2(a, nx, ny, nz)
    implicit none
    real :: a
    real, dimension(4) :: x0, y0, z0
    integer :: ix, iy, iz, i, nx, ny, nz
    real, allocatable :: x(:), y(:), z(:)

    allocate(x(nx * ny * nz * 4), y(nx * ny * nz * 4), z(nx * ny * nz * 4))
    open(11, file="fcc.txt", position="rewind", status="old")

    x0(1) = 0.0; y0(1) = 0.0; z0(1) = 0.0
    x0(2) = 0.0; y0(2) = 0.5 * a; z0(2) = 0.5 * a
    x0(3) = 0.5 * a; y0(3) = 0.0; z0(3) = 0.5 * a
    x0(4) = 0.5 * a; y0(4) = 0.5 * a; z0(4) = 0.0

    do iy = 0, ny - 1
        do ix = 0, nx - 1
            do iz = 0, nz - 1
                do i = 1, 4
                    x(i) = ix * a + x0(i)
                    y(i) = iy * a + y0(i)
                    z(i) = iz * a + z0(i)
                    write(11,10) x(i), y(i), z(i)
10                  format(f20.17, 4x, f20.17, 4x, f20.17)
                end do
            end do
        end do
    end do
    close(11)
end subroutine fcc2

!! Temperature control with velocity scaling
subroutine velscale()
    use vardecs
    use commondefs
    implicit none

    real :: scalefactor
    integer :: i, scalecount

    sumvel = 0.0
    do i = 1, n
        sumvel = sumvel + (atoms(i)%vx**2 + atoms(i)%vy**2 + atoms(i)%vz**2) * molwt(atoms(i)%kind)
    end do

    ! scale velocities at scalestep intervals
    if (tstep <= scalestep) then
        scalecount = mod(tstep, scaleint)
        if (scalecount == 0) then
            scalefactor = sqrt(3.0 * real(n) * bolz * temp / sumvel)
            do i = 1, np
                atoms(i)%vx = atoms(i)%vx * scalefactor
                atoms(i)%vy = atoms(i)%vy * scalefactor
                atoms(i)%vz = atoms(i)%vz * scalefactor
            end do
        end if
    end if
end subroutine velscale

!! Andersen pressure control with velocity Verlet predictor step
subroutine anders1()
    use vardecs
    use commondefs
    implicit none

    real :: sumro
    integer :: i

    sumro = 0.0
    do i = 1, n
        atoms(i)%rox = atoms(i)%x / boxsize_x
        atoms(i)%roy = atoms(i)%y / boxsize_y
        atoms(i)%roz = atoms(i)%z / boxsize_z

        atoms(i)%rox1 = atoms(i)%vx / boxsize_x
        atoms(i)%royl = atoms(i)%vy / boxsize_y
        atoms(i)%rozl = atoms(i)%vz / boxsize_z

        sumro = sumro + (atoms(i)%rox1**2 + atoms(i)%royl**2 + atoms(i)%rozl**2) * molwt(atoms(i)%kind)
    enddo

    do i = 1, n
        atoms(i)%rox = atoms(i)%rox + tstp * atoms(i)%rox1 + 0.5 * tstp**2 * &
                       ((atoms(i)%fx / molwt(atoms(i)%kind)) / boxsize_x - &
                        2.0/3.0 * atoms(i)%rox1 * volumel / volume)

        atoms(i)%roy = atoms(i)%roy + tstp * atoms(i)%royl + 0.5 * tstp**2 * &
                       ((atoms(i)%fy / molwt(atoms(i)%kind)) / boxsize_y - &
                        2.0/3.0 * atoms(i)%royl * volumel / volume)

        atoms(i)%roz = atoms(i)%roz + tstp * atoms(i)%rozl + 0.5 * tstp**2 * &
                       ((atoms(i)%fz / molwt(atoms(i)%kind)) / boxsize_z - &
                        2.0/3.0 * atoms(i)%rozl * volumel / volume)

        atoms(i)%rox1 = atoms(i)%rox1 + 0.5 * tstp * &
                        ((atoms(i)%fx / molwt(atoms(i)%kind)) / boxsize_x - &
                         2.0/3.0 * atoms(i)%rox1 * volumel / volume)

        atoms(i)%royl = atoms(i)%royl + 0.5 * tstp * &
                        ((atoms(i)%fy / molwt(atoms(i)%kind)) / boxsize_y - &
                         2.0/3.0 * atoms(i)%royl * volumel / volume)

        atoms(i)%rozl = atoms(i)%rozl + 0.5 * tstp * &
                        ((atoms(i)%fz / molwt(atoms(i)%kind)) / boxsize_z - &
                         2.0/3.0 * atoms(i)%rozl * volumel / volume)
    end do

    sumvel = 0.0
    do i = 1, n
        atoms(i)%x = atoms(i)%rox * boxsize_x
        atoms(i)%y = atoms(i)%roy * boxsize_y
        atoms(i)%z = atoms(i)%roz * boxsize_z

        atoms(i)%vx = atoms(i)%rox1 * boxsize_x
        atoms(i)%vy = atoms(i)%royl * boxsize_y
        atoms(i)%vz = atoms(i)%rozl * boxsize_z

        sumvel = sumvel + (atoms(i)%vx**2 + atoms(i)%vy**2 + atoms(i)%vz**2) * molwt(atoms(i)%kind)
    enddo

    tempdum = (sumvel) / (3.0 * real(n) * bolz)
    pressure = ((n * bolz * tempdum) + (1.0/3.0) * totw) / volume

    volumel = volumel + 0.5 * tstp**2 * (pressure - pset) / m
    volume = volume + tstp * volumel + 0.5 * tstp**2 * (pressure - pset) / m

    boxsize_x = (volume**0.33333333) * bxr
    boxsize_y = (volume**0.33333333) * byr
    boxsize_z = (volume**0.33333333) * bzr

    sumvel = 0.0
    do i = 1, n
        atoms(i)%x = atoms(i)%rox * boxsize_x
        atoms(i)%y = atoms(i)%roy * boxsize_y
        atoms(i)%z = atoms(i)%roz * boxsize_z

        atoms(i)%vx = atoms(i)%rox1 * boxsize_x
        atoms(i)%vy = atoms(i)%royl * boxsize_y
        atoms(i)%vz = atoms(i)%rozl * boxsize_z

        sumvel = sumvel + (atoms(i)%vx**2 + atoms(i)%vy**2 + atoms(i)%vz**2) * molwt(atoms(i)%kind)
    enddo
end subroutine anders1

!! Andersen pressure control with velocity Verlet corrector step
subroutine anders2()
    use vardecs
    use commondefs
    implicit none

    integer :: i

    do i = 1, n
        atoms(i)%rox1 = atoms(i)%rox1 + 0.5 * tstp * &
                        ((atoms(i)%fx / molwt(atoms(i)%kind)) / boxsize_x - &
                         2.0/3.0 * atoms(i)%rox1 * volumel / volume)

        atoms(i)%royl = atoms(i)%royl + 0.5 * tstp * &
                        ((atoms(i)%fy / molwt(atoms(i)%kind)) / boxsize_y - &
                         2.0/3.0 * atoms(i)%royl * volumel / volume)

        atoms(i)%rozl = atoms(i)%rozl + 0.5 * tstp * &
                        ((atoms(i)%fz / molwt(atoms(i)%kind)) / boxsize_z - &
                         2.0/3.0 * atoms(i)%rozl * volumel / volume)
    end do

    sumvel = 0.0
    do i = 1, n
        atoms(i)%x = atoms(i)%rox * boxsize_x
        atoms(i)%y = atoms(i)%roy * boxsize_y
        atoms(i)%z = atoms(i)%roz * boxsize_z

        atoms(i)%vx = atoms(i)%rox1 * boxsize_x
        atoms(i)%vy = atoms(i)%royl * boxsize_y
        atoms(i)%vz = atoms(i)%rozl * boxsize_z

        sumvel = sumvel + (atoms(i)%vx**2 + atoms(i)%vy**2 + atoms(i)%vz**2) * molwt(atoms(i)%kind)
    enddo

    tempdum = (sumvel) / (3.0 * real(n) * bolz)
    pressure = ((n * bolz * tempdum) + (1.0/3.0) * totw) / volume
    volumel = volumel + 0.5 * tstp * (pressure - pset) / m
end subroutine anders2

!! Berendsen pressure and temperature control
subroutine berendsen()
    use vardecs
    use commondefs
    implicit none

    integer :: i

    sumvel = 0.0
    do i = 1, n
        sumvel = sumvel + (atoms(i)%vx**2 + atoms(i)%vy**2 + atoms(i)%vz**2) * molwt(atoms(i)%kind)
    end do

    tempdum = (sumvel) / (3.0 * real(n) * bolz)
    pressure = ((n * bolz * tempdum) + (1.0/3.0) * totw) / volume

    mu = (1.0 + (tstp / tau_p) * beta * (pressure - pset))**(1.0/3.0)
    lambda = (1.0 + (tstp / tau_t) * (temp / tempdum - 1.0))**0.5

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

!! Gear 5th order integration predictor step
subroutine predict(dt)
    use vardecs
    use commondefs
    implicit none

    integer :: i
    real :: c1, c2, c3, c4
    real :: dt

    c1 = dt
    c2 = c1 * dt / 2.0
    c3 = c2 * dt / 3.0
    c4 = c3 * dt / 4.0

    do i = 1, n
        atoms(i)%x = atoms(i)%x + c1 * atoms(i)%vx + c2 * atoms(i)%ax + &
                     c3 * atoms(i)%bx + c4 * atoms(i)%cx
        atoms(i)%y = atoms(i)%y + c1 * atoms(i)%vy + c2 * atoms(i)%ay + &
                     c3 * atoms(i)%by + c4 * atoms(i)%cy
        atoms(i)%z = atoms(i)%z + c1 * atoms(i)%vz + c2 * atoms(i)%az + &
                     c3 * atoms(i)%bz + c4 * atoms(i)%cz

        atoms(i)%vx = atoms(i)%vx + c1 * atoms(i)%ax + c2 * atoms(i)%bx + &
                      c3 * atoms(i)%cx
        atoms(i)%vy = atoms(i)%vy + c1 * atoms(i)%ay + c2 * atoms(i)%by + &
                      c3 * atoms(i)%cy
        atoms(i)%vz = atoms(i)%vz + c1 * atoms(i)%az + c2 * atoms(i)%bz + &
                      c3 * atoms(i)%cz

        atoms(i)%ax = atoms(i)%ax + c1 * atoms(i)%bx + c2 * atoms(i)%cx
        atoms(i)%ay = atoms(i)%ay + c1 * atoms(i)%by + c2 * atoms(i)%cy
        atoms(i)%az = atoms(i)%az + c1 * atoms(i)%bz + c2 * atoms(i)%cz

        atoms(i)%bx = atoms(i)%bx + c1 * atoms(i)%cx
        atoms(i)%by = atoms(i)%by + c1 * atoms(i)%cy
        atoms(i)%bz = atoms(i)%bz + c1 * atoms(i)%cz
    end do
    return
end subroutine predict

!! Corrector stage of Gear 5th order integration
subroutine correct(dt)
    use vardecs
    use commondefs
    implicit none

    real :: axi, ayi, azi, corrx, corry, corrz
    real :: c1, c2, c3, c4
    real :: gear0, gear1, gear3, gear4
    real :: cr, cv, cb, cc
    real :: dt
    integer :: i

    gear0 = 19.0 / 120.0
    gear1 = 3.0 / 4.0
    gear3 = 1.0 / 2.0
    gear4 = 1.0 / 12.0

    c1 = dt
    c2 = c1 * dt / 2.0
    c3 = c2 / 3.0
    c4 = c3 / 4.0

    cr = gear0 * c2
    cv = gear1 * c2 / c1
    cb = gear3 / c3
    cc = gear4 / c4

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
    return
end subroutine correct

!! velocity Verlet predictor step
subroutine verl1(dt)
    use vardecs
    use commondefs
    implicit none

    real :: dt, dtsq2, dt2
    integer :: i

    dt2 = dt / 2.0
    dtsq2 = dt * dt2

    do i = 1, np
        atoms(i)%x = atoms(i)%x + dt * atoms(i)%vx + dtsq2 * (atoms(i)%fx / molwt(atoms(i)%kind))
        atoms(i)%y = atoms(i)%y + dt * atoms(i)%vy + dtsq2 * (atoms(i)%fy / molwt(atoms(i)%kind))
        atoms(i)%z = atoms(i)%z + dt * atoms(i)%vz + dtsq2 * (atoms(i)%fz / molwt(atoms(i)%kind))

        atoms(i)%vx = atoms(i)%vx + dt2 * (atoms(i)%fx / molwt(atoms(i)%kind))
        atoms(i)%vy = atoms(i)%vy + dt2 * (atoms(i)%fy / molwt(atoms(i)%kind))
        atoms(i)%vz = atoms(i)%vz + dt2 * (atoms(i)%fz / molwt(atoms(i)%kind))
    end do
    return
end subroutine verl1

!! velocity Verlet corrector step
subroutine verl2(dt)
    use vardecs
    use commondefs
    implicit none

    real :: dt, dt2
    integer :: i

    dt2 = dt / 2.0
    totkin = 0.0

    do i = 1, np
        atoms(i)%vx = atoms(i)%vx + dt2 * (atoms(i)%fx / molwt(atoms(i)%kind))
        atoms(i)%vy = atoms(i)%vy + dt2 * (atoms(i)%fy / molwt(atoms(i)%kind))
        atoms(i)%vz = atoms(i)%vz + dt2 * (atoms(i)%fz / molwt(atoms(i)%kind))

        totkin = totkin + (atoms(i)%vx**2 + atoms(i)%vy**2 + atoms(i)%vz**2) * molwt(atoms(i)%kind)
    end do
    totkin = totkin * 0.5
    return
end subroutine verl2

!! Velocity autocorrelation function calculation
subroutine vacf(name1, tstep, sstep, dtime, n, icf)
    implicit none

    real :: vaf
    integer :: i, j, n, step, nsample, tstep, sstep, icf
    real :: dtime, time
    real :: vx0(n), vy0(n), vz0(n), x0(n), y0(n), z0(n)
    real :: vx(n), vy(n), vz(n), x(n), y(n), z(n)
    real :: v0(n), r0(n)
    character(*) :: name1
    character(len=len_trim(name1)) :: name2

    name2 = trim(name1)
    vaf = 0.0

    open(unit=3, file=name2//".con")
    open(unit=2, file=name2//".vaf")

    read(3,*) step
    do i = 1, n
        read(3,31) x0(i), y0(i), z0(i), vx0(i), vy0(i), vz0(i)
31      format(1x, 6es15.7)
        v0(i) = sqrt(vx0(i)*vx0(i) + vy0(i)*vy0(i) + vz0(i)*vz0(i))
        r0(i) = sqrt(x0(i)**2 + y0(i)**2 + z0(i)**2)
        vaf = vaf + (v0(i) * v0(i)) / (v0(i) * v0(i))
    end do

    vaf = vaf / real(n)
    time = real((step - (sstep + icf))) * dtime
    nsample = (tstep - sstep) / icf

    write(2,*) time, vaf

    do j = 1, nsample
        vaf = 0.0
        read(3,*) step
        do i = 1, n
            read(3,32) x(i), y(i), z(i), vx(i), vy(i), vz(i)
32          format(1x, 6es15.7)
            vaf = vaf + (vx0(i)*vx(i) + vy0(i)*vy(i) + vz0(i)*vz(i)) / (v0(i) * v0(i))
        end do
        vaf = vaf / real(n)
        time = real((step - (sstep + icf))) * dtime
        write(2,*) time, vaf
    end do

    close(2)
    close(3)
end subroutine vacf


!! Declarations for MD
module vardecs
    implicit none

    type lj
        real :: eps, sigma
    end type lj

    type fourier
        real :: a, b, c, d, e, f, g, h, ix, jx, kx, lx, mx, rmin, rmax, eps, sigma
    end type fourier

    type poten
        integer :: index, pot
    end type poten

    type(poten) :: pairtype(50)
    type(lj) :: ljpot(20)
    type(fourier) :: fourierpot(20)
    integer :: totstep, np, enunit, scaleint, scalestep, nperf, g
    real :: tote, totkin, hamilton, totw, volume, volumel, boxsize_xl, deltat
    real :: temp, tunit, tstp, tottime
    real :: seed, tempdum, sumvel, tol
    real :: nbox_x, nbox_y, nbox_z
    real :: bolz
    integer :: tstep, iconf, imovie, ivac, iforce, imsd, ielmigstep, msdstart, msdend
    real :: q, ps, s
    real :: eps, sigma, encorr
    real :: pressure, pset, m
    real :: mu, tau_t, tau_p, lambda, beta
    real, dimension(10) :: molwt
    real, allocatable :: xini(:), yini(:), zini(:)
    real, parameter :: pi = 3.14159265358979
end module vardecs

! data.h
! Wrapper functions enabling dfridr to differentiate uu_s along x, y, and z directions

real function uux(var)
    use defs
    use commondefs
    implicit none
    real :: var
    real, external :: uu_s
    uux = uu_s(1, var)
end function uux

real function uuy(var)
    use defs
    use commondefs
    implicit none
    real :: var
    real, external :: uu_s
    uuy = uu_s(2, var)
end function uuy

real function uuz(var)
    use defs
    use commondefs
    implicit none
    real :: var
    real, external :: uu_s
    uuz = uu_s(3, var)
end function uuz
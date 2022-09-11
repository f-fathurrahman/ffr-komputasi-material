subroutine excitations
  use states
  use mesh
  use fft
  use poisson
  implicit none


  !variables global to the subroutine
  integer, allocatable :: pair_i(:), pair_a(:)  ! holds separated indices of compund index ia
  real(8), allocatable :: mat(:,:)              ! general purpose matrix
  real(8), allocatable :: energies(:,:)
  real(8), allocatable :: eigenval(:)         ! excitation energies and intensities
  real(8), allocatable :: os(:, :)
  real(8), allocatable :: vext(:, :), vx(:, :), vc(:, :), vh(:, :), vtot(:, :), rho(:, :), hwf(:, :)
  !logical, allocatable :: saved_K(:, :)         ! which matrix elements have been loaded
  integer :: n_pairs, ia, jb, a, i, b, j
  real(8) :: temp, ex, ec

  write(*, '(a)') 'Initializing FFTs...'
  call fft_all_init()
  call poisson_init()
  write(*, '(a,/)') 'Done.'

  write(*, *) 'Reading wavefunctions...'
  call read_wfs('wfs.dat')
  write(*, *) 'Done.'

  write(*,*) 'n_wf = ', n_wf
  allocate(eigenval(n_wf))

  ! Now we should get the potentials, eigenvalues, etc...
  allocate(vext(n, n), vh(n, n), vc(n, n), vx(n, n), vtot(n, n), hwf(n, n), rho(n, n))
  call external_pot(vext)
  call build_rho(wfs, rho)
  call interaction_pot(rho, vh, vx, vc, ex, ec)
  vtot = vext + vh + vx + vc
  !
  write(*,*) 'Eigenvalues: '
  do j = 1, N_wf
    call hpsi(vtot, wfs(:, :, j), hwf)
    eigenval(j) = dotproduct(wfs(:, :, j), hwf)
    write(*,*) j, eigenval(j)
  enddo

  n_pairs = N_occ*N_empty
  allocate(pair_i(n_pairs), pair_a(n_pairs))

  j = 1
  do a = n_occ+1, n_wf
    do i = 1, n_occ
       pair_i(j) = i
       pair_a(j) = a
       j = j + 1
    enddo
  enddo

  allocate(mat(n_pairs, n_pairs), energies(n_pairs, 3))
  mat = 0.0d0
  energies = 0.0d0

  ! calculate the matrix elements of (v + fxc)
  write(*,*) 'Calculating matrix elements of (v + fxc): '
  do ia = 1, n_pairs
    i = pair_i(ia)
    a = pair_a(ia)
    do jb = 1, n_pairs
       j = pair_i(jb)
       b = pair_a(jb)
       mat(ia, jb) = k_term(i, a, j, b)
    enddo
  enddo

  ! diagonal
  do ia = 1, n_pairs
    i = pair_i(ia)
    a = pair_a(ia)
    temp = eigenval(a) - eigenval(i)
    do jb = 1, n_pairs
      j = pair_i(jb)
      b = pair_a(jb)
      mat(ia, jb)  = 4.0d0 * sqrt(temp) * mat(ia, jb) * sqrt(eigenval(b) - eigenval(j))
      if(jb /= ia) mat(jb, ia) = mat(ia, jb) ! the matrix is symmetric
    enddo
    mat(ia, ia) = temp**2 + mat(ia, ia)
  enddo

  !write(*, *)
  !write(*, *) mat


  ! Diagonalization of the matrix mat. Eigenvalues go to energies(:, 1)
  call eigensolve(n_pairs, mat, energies)
  energies(:, 1) = sqrt(energies(:, 1))

  allocate(os(n_pairs, 2))
  do ia = 1, n_pairs
     i = pair_i(ia)
     a = pair_a(ia)
     call dipole_matrix_elem(i, a, os(ia,:))
  enddo

  do ia = 1, n_pairs
    do j = 1, 2
        energies(ia,1+j) = 2.d0 * (sum(os(:,j)*mat(:,ia) * sqrt(eigenval(pair_a(:)) - eigenval(pair_i(:))) ))**2
    enddo
  enddo

  ! And write down the results
  write(*, '(4(a15,1x))') 'E' , '<x>', '<y>', '<f>'
  do ia = 1, n_pairs
     write(*, '(5(e15.8,1x))') energies(ia,1), &
     energies(ia, 2:3), (2.0/3.0)*sum(energies(ia, 2:3))
  end do
  write(*, *)
!!$  write(iunit, '(a14)', advance = 'no') ' value '
!!$  do ia = 1, n_pairs
!!$     write(iunit, '(3x,i4,a1,i4,2x)', advance='no') pair_i(ia), ' - ', pair_a(ia)
!!$  end do
!!$  write(iunit, '(1x)')
!!$      
!!$  do ia = 1, n_pairs
!!$    write(iunit, '(es14.6)', advance='no') energies(ia, 1) / units_out%energy%factor
!!$    temp = M_ONE
!!$    if(maxval(mat(:, ia)) < abs(minval(mat(:, ia)))) temp = -temp
!!$    do jb = 1, n_pairs
!!$       write(iunit, '(es14.6)', advance='no') temp*mat(jb, ia)
!!$    end do
!!$    write(iunit, '(1x)')
!!$  end do
      
  contains

  !-----------------------------
  subroutine eigensolve(k, a, e)
  !-----------------------------
    integer, intent(in)     :: k
    real(8), intent(inout)  :: a(k, k)
    real(8), intent(out)    :: e(k)

    integer :: info, lwork
    real(8), allocatable :: work(:)

    lwork = 6*k
    allocate(work(lwork))
    call dsyev('V', 'U', k, a(1,1), k, e(1), work(1), lwork, info)
    deallocate(work)
  !------------------------
  end subroutine eigensolve
  !------------------------


  !----------------------------------
  real(8) function k_term(i, a, j, b)
  !----------------------------------
    integer, intent(in) :: i, a, j, b
    real(8), allocatable :: rho_i(:, :), rho_j(:, :), pot(:, :)
    integer :: ix, iy
    real(8) :: fxc

    allocate(rho_i(n, n), rho_j(n, n), pot(n, n))

    rho_i(:, :) = wfs(:, :, i)*wfs(:, :, a)
    rho_j(:, :) = wfs(:, :, j)*wfs(:, :, b)

    pot = 0.0d0
    call poisson_solve(rho_j, pot)
    k_term = dotproduct(rho_i, pot)

    ! now we have fxc
    do ix = 1, n
       do iy = 1, n
          call fxc_LDA(rho(ix, iy), fxc)
          K_term = K_term + rho_i(ix, iy)*rho_j(ix, iy)*fxc*delta**2
       enddo
    enddo

    deallocate(rho_i, rho_j, pot)
  !------------------
  end function k_term
  !------------------


  !-------------------------------------
  subroutine dipole_matrix_elem(i, j, s)
  !-------------------------------------
    integer, intent(in) :: i, j
    real(8), intent(out) :: s(2)

    integer :: ix, iy

    s = 0.0d0
    do ix = 1, n
       do iy = 1, n
          s(1) = s(1) + x(ix, iy)*wfs(ix, iy, i)*wfs(ix, iy, j)
          s(2) = s(2) + y(ix, iy)*wfs(ix, iy, i)*wfs(ix, iy, j)
       enddo
    enddo
    s = s*delta**2
  !--------------------------------
  end subroutine dipole_matrix_elem
  !--------------------------------



end subroutine excitations

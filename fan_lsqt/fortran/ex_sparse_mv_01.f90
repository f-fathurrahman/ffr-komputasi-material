subroutine ccs_matvec(n, col_ptr, row_idx, values, x, y)
    implicit none
    integer, intent(in) :: n                 ! Matrix size (n x n)
    integer, intent(in) :: col_ptr(n+1)      ! Column pointers (size n+1)
    integer, intent(in) :: row_idx(*)        ! Row indices of nonzeros
    real(8), intent(in) :: values(*)         ! Nonzero values
    real(8), intent(in) :: x(n)              ! Input vector
    real(8), intent(out) :: y(n)             ! Output vector
    
    integer :: j, i
    
    ! Initialize output vector to zero
    y = 0.0d0
    
    ! Loop over columns
    do j = 1, n
        do i = col_ptr(j), col_ptr(j+1) - 1
            y(row_idx(i)) = y(row_idx(i)) + values(i) * x(j)
        end do
    end do
end subroutine ccs_matvec

program test_ccs_matvec
    implicit none
    integer, parameter :: n = 3
    integer :: col_ptr(n+1), row_idx(5)
    real(8) :: values(5), x(n), y(n)
    
    ! Corrected CCS format for:
    !  | 1 0 2 |
    !  | 0 3 0 |
    !  | 4 0 5 |
    col_ptr = (/1, 3, 4, 6/)
    row_idx = (/1, 3, 2, 1, 3/)
    values  = (/1.0d0, 4.0d0, 3.0d0, 2.0d0, 5.0d0/)
    x = (/1.0d0, 2.0d0, 3.0d0/)
    
    call ccs_matvec(n, col_ptr, row_idx, values, x, y)
    
    print *, 'Resulting vector y: ', y
end program test_ccs_matvec

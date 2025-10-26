module grid_mod
  use params
  implicit none
  private
  public :: create_grid, validate_grid, report_grid_info
  public :: Nx, Ny, dx, dy, x_coords, y_coords
  public :: inside_mask, boundary_mask, idx2i, idx2j

  !-------------------------------------------------------------
  ! Grid size (global parameters)
  !-------------------------------------------------------------
  integer, parameter :: Nx = Nx_glob
  integer, parameter :: Ny = Ny_glob

  real(dp) :: dx, dy
  real(dp), allocatable :: x_coords(:), y_coords(:)
  logical, allocatable :: inside_mask(:,:), boundary_mask(:,:)

contains

  subroutine create_grid()
    !-------------------------------------------------------------
    ! Allocate and initialize grid coordinates and masks
    !-------------------------------------------------------------
    integer :: i, j
    real(dp) :: x_start, y_start

    ! Safe deallocation
    if (allocated(x_coords))      deallocate(x_coords)
    if (allocated(y_coords))      deallocate(y_coords)
    if (allocated(inside_mask))   deallocate(inside_mask)
    if (allocated(boundary_mask)) deallocate(boundary_mask)

    ! Validate grid size
    if (Nx < 2 .or. Ny < 2) then
      error stop "create_grid: Nx and Ny must be >= 2"
    end if

    ! Domain limits
    x_start = Xmin
    y_start = Ymin

    dx = Lx_glob / real(Nx-1, dp)
    dy = Ly_glob / real(Ny-1, dp)

    ! Allocate arrays
    allocate(x_coords(Nx))
    allocate(y_coords(Ny))
    allocate(inside_mask(Nx,Ny))
    allocate(boundary_mask(Nx,Ny))

    ! Fill coordinate arrays
    do i = 1, Nx
      x_coords(i) = x_start + real(i-1, dp) * dx
    end do
    do j = 1, Ny
      y_coords(j) = y_start + real(j-1, dp) * dy
    end do

    ! Create masks for boundary and interior points
    !$omp parallel do collapse(2) private(i,j) schedule(static)
    do j = 1, Ny
      do i = 1, Nx
        boundary_mask(i,j) = (i == 1 .or. i == Nx .or. j == 1 .or. j == Ny)
        inside_mask(i,j)   = .not. boundary_mask(i,j)
      end do
    end do

    ! Validate grid quality
    call validate_grid()
    call report_grid_info()
  end subroutine create_grid

  subroutine validate_grid()
    !-------------------------------------------------------------
    ! Checks uniformity, spacing, and aspect ratio
    !-------------------------------------------------------------
    real(dp) :: dx_ref, dy_ref, tol, aspect_ratio
    integer :: i

    tol = 1.0e-12_dp

    if (.not. allocated(x_coords) .or. .not. allocated(y_coords)) then
      error stop "validate_grid: coordinate arrays not allocated"
    end if

    ! Reference spacing
    dx_ref = (x_coords(Nx) - x_coords(1)) / real(Nx-1, dp)
    dy_ref = (y_coords(Ny) - y_coords(1)) / real(Ny-1, dp)

    ! Uniformity check
    do i = 2, Nx
      if (abs((x_coords(i)-x_coords(i-1)) - dx_ref) > tol) then
        write(*,*) "validate_grid: Non-uniform spacing in x at index=", i
      end if
    end do
    do i = 2, Ny
      if (abs((y_coords(i)-y_coords(i-1)) - dy_ref) > tol) then
        write(*,*) "validate_grid: Non-uniform spacing in y at index=", i
      end if
    end do

    ! Aspect ratio metric
    aspect_ratio = dx_ref / dy_ref
    if (aspect_ratio < 0.2_dp .or. aspect_ratio > 5.0_dp) then
      write(*,'(A,F8.3)') "Warning: poor aspect ratio detected: ", aspect_ratio
    end if
  end subroutine validate_grid

  subroutine report_grid_info()
    !-------------------------------------------------------------
    ! Prints grid metrics and spacing
    !-------------------------------------------------------------
    real(dp) :: aspect_ratio

    if (.not. allocated(x_coords) .or. .not. allocated(y_coords)) then
      write(*,*) "report_grid_info: grid not allocated"
      return
    end if

    aspect_ratio = dx / dy

    write(*,'(A)') "-------------------------------------------"
    write(*,'(A)') " Grid Information"
    write(*,'(A)') "-------------------------------------------"
    write(*,'(A,I6)') " Nx = ", Nx
    write(*,'(A,I6)') " Ny = ", Ny
    write(*,'(A,ES12.5)') " dx = ", dx
    write(*,'(A,ES12.5)') " dy = ", dy
    write(*,'(A,ES12.5)') " Aspect ratio (dx/dy) = ", aspect_ratio
    write(*,'(A,ES12.5,1X,ES12.5)') " Domain extents Xmin,Xmax = ", Xmin, Xmax
    write(*,'(A,ES12.5,1X,ES12.5)') " Domain extents Ymin,Ymax = ", Ymin, Ymax
    write(*,'(A)') "-------------------------------------------"
  end subroutine report_grid_info

  pure integer function idx2i(idx) result(i)
    ! Convert linear index to i-coordinate
    integer, intent(in) :: idx
    i = ((idx-1) / Ny) + 1
  end function idx2i

  pure integer function idx2j(idx) result(j)
    ! Convert linear index to j-coordinate
    integer, intent(in) :: idx
    j = mod(idx-1, Ny) + 1
  end function idx2j

end module grid_mod

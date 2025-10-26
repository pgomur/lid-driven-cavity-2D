module numerics_mod
  use params
  use grid_mod
  implicit none
  private

  ! Public routines
  public :: gradient_x, gradient_y, gradient_x_quick, gradient_y_quick
  public :: laplacian_2d, apply_dirichlet_bc_scalar, apply_neumann_bc_scalar
  public :: check_field_validity, compute_vorticity, smooth_field
  public :: add, subtract, multiply, divide
  public :: power, sqrt_field, abs_field, min_field, max_field
  public :: sin_field, cos_field, tan_field
  public :: exp_field, log_field, log10_field
  public :: floor_field, ceil_field, sign_field
  public :: mean_field, sum_field

contains

  subroutine gradient_x(f, gx, dx)
    ! Central difference in x (interior), one-sided at boundaries
    real(dp), intent(in) :: f(:,:)
    real(dp), intent(out) :: gx(size(f,1), size(f,2))
    real(dp), intent(in) :: dx
    integer :: i, j, nx, ny
    nx = size(f,1)
    ny = size(f,2)
    gx = 0.0_dp

    !$omp parallel do collapse(2) private(i,j) schedule(static)
    do j = 1, ny
      do i = 2, nx-1
        gx(i,j) = ( f(i+1,j) - f(i-1,j) ) / (2.0_dp*dx)
      end do
    end do

    ! Boundaries (forward/backward)
    do j = 1, ny
      gx(1,j) = ( f(2,j) - f(1,j) ) / dx
      gx(nx,j) = ( f(nx,j) - f(nx-1,j) ) / dx
    end do
  end subroutine gradient_x

  subroutine gradient_y(f, gy, dy)
    ! Central difference in y
    real(dp), intent(in) :: f(:,:)
    real(dp), intent(out) :: gy(size(f,1), size(f,2))
    real(dp), intent(in) :: dy
    integer :: i, j, nx, ny
    nx = size(f,1)
    ny = size(f,2)
    gy = 0.0_dp

    !$omp parallel do collapse(2) private(i,j) schedule(static)
    do i = 1, nx
      do j = 2, ny-1
        gy(i,j) = ( f(i,j+1) - f(i,j-1) ) / (2.0_dp*dy)
      end do
    end do

    ! Boundaries
    do i = 1, nx
      gy(i,1) = ( f(i,2) - f(i,1) ) / dy
      gy(i,ny) = ( f(i,ny) - f(i,ny-1) ) / dy
    end do
  end subroutine gradient_y

  subroutine gradient_x_quick(f, gx, dx)
    ! QUICK scheme 3er orden en x
    real(dp), intent(in) :: f(:,:)
    real(dp), intent(out) :: gx(size(f,1), size(f,2))
    real(dp), intent(in) :: dx
    integer :: i, j, nx, ny
    nx = size(f,1)
    ny = size(f,2)
    gx = 0.0_dp

    !$omp parallel do collapse(2) private(i,j) schedule(static)
    do j = 1, ny
      do i = 3, nx-2
        gx(i,j) = (-f(i+2,j) + 6.0_dp*f(i+1,j) - 3.0_dp*f(i,j) - 2.0_dp*f(i-1,j)) / (6.0_dp*dx)
      end do
    end do

    ! Bordes: fallback a central
    call gradient_x(f, gx, dx)
  end subroutine gradient_x_quick

  subroutine gradient_y_quick(f, gy, dy)
    ! QUICK scheme 3er orden en y
    real(dp), intent(in) :: f(:,:)
    real(dp), intent(out) :: gy(size(f,1), size(f,2))
    real(dp), intent(in) :: dy
    integer :: i, j, nx, ny
    nx = size(f,1)
    ny = size(f,2)
    gy = 0.0_dp

    !$omp parallel do collapse(2) private(i,j) schedule(static)
    do i = 1, nx
      do j = 3, ny-2
        gy(i,j) = (-f(i,j+2) + 6.0_dp*f(i,j+1) - 3.0_dp*f(i,j) - 2.0_dp*f(i,j-1)) / (6.0_dp*dy)
      end do
    end do

    ! Bordes: fallback a central
    call gradient_y(f, gy, dy)
  end subroutine gradient_y_quick

  subroutine laplacian_2d(f, dx_in, dy_in, lap)
    ! 5-point Laplacian
    real(dp), intent(in) :: f(:,:)
    real(dp), intent(in) :: dx_in, dy_in
    real(dp), intent(out) :: lap(size(f,1), size(f,2))
    integer :: i, j, nx, ny
    nx = size(f,1)
    ny = size(f,2)
    lap = 0.0_dp

    !$omp parallel do collapse(2) private(i,j) schedule(static)
    do j = 2, ny-1
      do i = 2, nx-1
        lap(i,j) = ( f(i+1,j) - 2.0_dp*f(i,j) + f(i-1,j) ) / (dx_in*dx_in) + &
                   ( f(i,j+1) - 2.0_dp*f(i,j) + f(i,j-1) ) / (dy_in*dy_in)
      end do
    end do
  end subroutine laplacian_2d

  subroutine apply_dirichlet_bc_scalar(f, top, bottom, left, right)
    ! Aplica Dirichlet
    real(dp), intent(inout) :: f(:,:)
    real(dp), intent(in) :: top, bottom, left, right
    integer :: i, j, nx, ny
    nx = size(f,1)
    ny = size(f,2)

    do i = 1, nx
      f(i,ny) = top
      f(i,1) = bottom
    end do
    do j = 1, ny
      f(1,j) = left
      f(nx,j) = right
    end do
  end subroutine apply_dirichlet_bc_scalar

  subroutine apply_neumann_bc_scalar(f, dtop, dbottom, dleft, dright)
    ! Aplica Neumann: ∂f/∂n = valor
    real(dp), intent(inout) :: f(:,:)
    real(dp), intent(in) :: dtop, dbottom, dleft, dright
    integer :: i, j, nx, ny
    nx = size(f,1)
    ny = size(f,2)

    do i = 2, nx-1
      f(i,ny) = f(i,ny-1) + dtop
      f(i,1)  = f(i,2)  + dbottom
    end do
    do j = 2, ny-1
      f(1,j)  = f(2,j)  + dleft
      f(nx,j) = f(nx-1,j) + dright
    end do
  end subroutine apply_neumann_bc_scalar

  subroutine check_field_validity(f)
    ! Comprueba NaN e Inf
    real(dp), intent(in) :: f(:,:)
    integer :: i,j
    logical :: has_invalid
    has_invalid = .false.
    !$omp parallel do collapse(2) private(i,j) reduction(.or.:has_invalid)
    do j = 1, size(f,2)
      do i = 1, size(f,1)
        if (.not.(f(i,j) == f(i,j))) has_invalid = .true. ! NaN
        if (abs(f(i,j)) > huge(1.0_dp)) has_invalid = .true. ! Inf
      end do
    end do
    !$omp end parallel do
    if (has_invalid) stop "ERROR: Field contains NaN or Inf"
  end subroutine check_field_validity

  subroutine compute_vorticity(u, v, omega, dx_in, dy_in)
    ! Vorticidad ω = ∂v/∂x - ∂u/∂y
    real(dp), intent(in) :: u(:,:), v(:,:)
    real(dp), intent(out) :: omega(size(u,1), size(u,2))
    real(dp), intent(in) :: dx_in, dy_in
    real(dp), allocatable :: gx(:,:), gy(:,:)
    allocate(gx(size(u,1),size(u,2)))
    allocate(gy(size(u,1),size(u,2)))
    call gradient_x(v, gx, dx_in)
    call gradient_y(u, gy, dy_in)
    omega = gx - gy
    deallocate(gx, gy)
  end subroutine compute_vorticity

  subroutine smooth_field(f, f_smooth)
    ! Filtro simple promedio 3x3
    real(dp), intent(in) :: f(:,:)
    real(dp), intent(out) :: f_smooth(size(f,1), size(f,2))
    integer :: i,j,nx,ny
    nx = size(f,1)
    ny = size(f,2)
    f_smooth = f
    !$omp parallel do collapse(2) private(i,j) schedule(static)
    do j = 2, ny-1
      do i = 2, nx-1
        f_smooth(i,j) = (f(i,j)+f(i+1,j)+f(i-1,j)+f(i,j+1)+f(i,j-1))/5.0_dp
      end do
    end do
  end subroutine smooth_field

  pure function add(a,b) result(c)
    real(dp), intent(in) :: a(:,:), b(:,:)
    real(dp) :: c(size(a,1), size(a,2))
    c = a + b
  end function add

  pure function subtract(a,b) result(c)
    real(dp), intent(in) :: a(:,:), b(:,:)
    real(dp) :: c(size(a,1), size(a,2))
    c = a - b
  end function subtract

  pure function multiply(a,b) result(c)
    real(dp), intent(in) :: a(:,:), b(:,:)
    real(dp) :: c(size(a,1), size(a,2))
    c = a * b
  end function multiply

  pure function divide(a,b) result(c)
    real(dp), intent(in) :: a(:,:), b(:,:)
    real(dp) :: c(size(a,1), size(a,2))
    c = a / b
  end function divide

  pure function power(a,p) result(c)
    real(dp), intent(in) :: a(:,:), p
    real(dp) :: c(size(a,1), size(a,2))
    c = a**p
  end function power

  pure function sqrt_field(a) result(c)
    real(dp), intent(in) :: a(:,:)
    real(dp) :: c(size(a,1), size(a,2))
    c = sqrt(a)
  end function sqrt_field

  pure function abs_field(a) result(c)
    real(dp), intent(in) :: a(:,:)
    real(dp) :: c(size(a,1), size(a,2))
    c = abs(a)
  end function abs_field

  pure function min_field(a,b) result(c)
    real(dp), intent(in) :: a(:,:), b(:,:)
    real(dp) :: c(size(a,1), size(a,2))
    c = min(a,b)
  end function min_field

  pure function max_field(a,b) result(c)
    real(dp), intent(in) :: a(:,:), b(:,:)
    real(dp) :: c(size(a,1), size(a,2))
    c = max(a,b)
  end function max_field

  pure function sin_field(a) result(c)
    real(dp), intent(in) :: a(:,:)
    real(dp) :: c(size(a,1), size(a,2))
    c = sin(a)
  end function sin_field

  pure function cos_field(a) result(c)
    real(dp), intent(in) :: a(:,:)
    real(dp) :: c(size(a,1), size(a,2))
    c = cos(a)
  end function cos_field

  pure function tan_field(a) result(c)
    real(dp), intent(in) :: a(:,:)
    real(dp) :: c(size(a,1), size(a,2))
    c = tan(a)
  end function tan_field

  pure function exp_field(a) result(c)
    real(dp), intent(in) :: a(:,:)
    real(dp) :: c(size(a,1), size(a,2))
    c = exp(a)
  end function exp_field

  pure function log_field(a) result(c)
    real(dp), intent(in) :: a(:,:)
    real(dp) :: c(size(a,1), size(a,2))
    c = log(a)
  end function log_field

  pure function log10_field(a) result(c)
    real(dp), intent(in) :: a(:,:)
    real(dp) :: c(size(a,1), size(a,2))
    c = log10(a)
  end function log10_field

  pure function floor_field(a) result(c)
      real(dp), intent(in) :: a(:,:)
      real(dp) :: c(size(a,1), size(a,2))
      integer :: i,j
      do j = 1, size(a,2)
          do i = 1, size(a,1)
              c(i,j) = real(floor(a(i,j)), dp)
          end do
      end do
  end function floor_field

  pure function ceil_field(a) result(c)
      real(dp), intent(in) :: a(:,:)
      real(dp) :: c(size(a,1), size(a,2))
      integer :: i,j
      do j = 1, size(a,2)
          do i = 1, size(a,1)
              c(i,j) = real(ceiling(a(i,j)), dp)
          end do
      end do
  end function ceil_field

  pure function sign_field(a) result(c)
    real(dp), intent(in) :: a(:,:)
    real(dp) :: c(size(a,1), size(a,2))
    c = sign(1.0_dp, a)
  end function sign_field

  pure function mean_field(a) result(c)
    real(dp), intent(in) :: a(:,:)
    real(dp) :: c
    c = sum(a)/real(size(a,1)*size(a,2), dp)
  end function mean_field

  pure function sum_field(a) result(c)
    real(dp), intent(in) :: a(:,:)
    real(dp) :: c
    c = sum(a)
  end function sum_field

end module numerics_mod

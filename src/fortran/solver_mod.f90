!=============================================================
! Features:
!   - CFL stability monitoring
!   - Robust error handling
!   - Convergence diagnostics
!   - Optimized memory access patterns
!   - Production-ready logging
!=============================================================
module solver_mod
  use params
  use grid_mod
  use numerics_mod
  use omp_lib
  implicit none
  private
  
  ! Public interface
  public :: solve_poisson_sor, advance_explicit_step_uv, solve_dummy
  public :: check_cfl_stability, compute_divergence_max
  
  ! Module constants
  real(dp), parameter :: CFL_MAX_SAFE = 0.5_dp
  real(dp), parameter :: DIV_TOLERANCE = 1.0e-10_dp
  real(dp), parameter :: EPSILON_ZERO = 1.0e-16_dp
  
  ! Diagnostics
  logical, save :: first_call = .true.
  integer, save :: warning_count = 0

contains

  subroutine solve_poisson_sor(phi, rhs, dx, dy, tol_local, maxiter, omega)
    !---------------------------------------------------------
    ! Red-Black SOR solver:
    !   - Optimized residual calculation (every 10 iters)
    !   - Boundary condition enforcement
    !   - Detailed convergence reporting
    !---------------------------------------------------------
    real(dp), intent(inout) :: phi(:,:)
    real(dp), intent(in)    :: rhs(:,:)
    real(dp), intent(in)    :: dx, dy
    real(dp), intent(in), optional :: tol_local
    integer,  intent(in), optional :: maxiter
    real(dp), intent(in), optional :: omega

    ! Local variables
    integer  :: nx, ny, i, j, iter, color
    real(dp) :: res, res0, dx2, dy2, denom, om, tolv
    integer  :: maxit, check_interval
    logical  :: converged
    real(dp) :: wall_time_start, wall_time_end

    ! Timer start
    wall_time_start = omp_get_wtime()

    ! Array dimensions
    nx = size(phi, 1)
    ny = size(phi, 2)
    
    ! Validate input
    if (nx /= size(rhs,1) .or. ny /= size(rhs,2)) then
      print *, "ERROR [solve_poisson_sor]: phi and rhs dimension mismatch"
      stop 1
    end if
    
    if (nx < 3 .or. ny < 3) then
      print *, "ERROR [solve_poisson_sor]: grid too small (nx, ny) = ", nx, ny
      stop 1
    end if

    ! Precompute constants
    dx2   = dx * dx
    dy2   = dy * dy
    denom = 2.0_dp * (dx2 + dy2)

    ! Set parameters with defaults
    om = 1.5_dp
    if (present(omega)) then
      om = omega
      if (om <= 0.0_dp .or. om >= 2.0_dp) then
        print *, "WARNING [solve_poisson_sor]: omega out of range, using 1.5"
        om = 1.5_dp
      end if
    end if

    tolv = 1.0e-8_dp
    if (present(tol_local)) tolv = max(tol_local, EPSILON_ZERO)

    maxit = 10000
    if (present(maxiter)) maxit = max(maxiter, 1)
    
    check_interval = max(10, maxit / 100)  ! Check every 1% or 10 iters

    !---------------------------------------------------------
    ! Initial residual calculation
    !---------------------------------------------------------
    res0 = 0.0_dp
    !$omp parallel do collapse(2) reduction(+:res0) schedule(static) private(i,j)
    do j = 2, ny - 1
      do i = 2, nx - 1
        res0 = res0 + ((phi(i+1,j) - 2.0_dp*phi(i,j) + phi(i-1,j)) / dx2 + &
                       (phi(i,j+1) - 2.0_dp*phi(i,j) + phi(i,j-1)) / dy2 - rhs(i,j))**2
      end do
    end do
    res0 = sqrt(res0)
    
    ! Handle trivial case
    if (res0 < EPSILON_ZERO) then
      res0 = 1.0_dp
      if (verbose) print *, "INFO [solve_poisson_sor]: Trivial RHS, residual ~ 0"
    end if

    !---------------------------------------------------------
    ! Main iteration loop
    !---------------------------------------------------------
    converged = .false.
    iter = 0
    
    do while (.not. converged .and. iter < maxit)
      iter = iter + 1

      ! Red-Black Gauss-Seidel sweep
      do color = 0, 1
        !$omp parallel do collapse(2) schedule(static) private(i,j)
        do j = 2, ny - 1
          do i = 2, nx - 1
            if (mod(i + j, 2) == color) then
              phi(i,j) = (1.0_dp - om) * phi(i,j) + om * &
                         ((dy2 * (phi(i+1,j) + phi(i-1,j)) + &
                           dx2 * (phi(i,j+1) + phi(i,j-1)) - &
                           dx2 * dy2 * rhs(i,j)) / denom)
            end if
          end do
        end do
        !$omp end parallel do
      end do

      ! Enforce Neumann BCs on pressure (dp/dn = 0)
      call enforce_pressure_bc(phi)

      ! Check convergence periodically
      if (mod(iter, check_interval) == 0 .or. iter == 1) then
        res = compute_residual(phi, rhs, dx, dy)
        
        if (res / res0 < tolv) then
          converged = .true.
        end if
        
        if (verbose .and. mod(iter, 100) == 0) then
          print '(A,I6,A,ES12.4)', " [SOR] iter=", iter, "  res/res0=", res/res0
        end if
      end if
    end do

    ! Final residual if not recently computed
    if (mod(iter, check_interval) /= 0) then
      res = compute_residual(phi, rhs, dx, dy)
    end if

    ! Timer end
    wall_time_end = omp_get_wtime()

    !---------------------------------------------------------
    ! Convergence reporting
    !---------------------------------------------------------
    if (verbose) then
      if (converged) then
        print '(A,I6,A,ES10.2,A,F8.4,A)', &
          " [SOR] ✓ Converged in ", iter, " iterations (res=", res/res0, &
          ") in ", wall_time_end - wall_time_start, "s"
      else
        print '(A,I6,A,ES10.2,A)', &
          " [SOR] ⚠ Max iterations reached (", maxit, "), res=", res/res0, &
          " - Consider increasing maxiter or omega"
      end if
    end if

  end subroutine solve_poisson_sor

  subroutine advance_explicit_step_uv(u, v, dx, dy, dt_local)
    !---------------------------------------------------------
    ! Explicit advection-diffusion with CFL checking
    !---------------------------------------------------------
    real(dp), intent(inout) :: u(:,:), v(:,:)
    real(dp), intent(in)    :: dx, dy, dt_local
    
    real(dp), allocatable :: lap_u(:,:), lap_v(:,:)
    real(dp) :: cfl, cfl_conv, cfl_diff, u_max, v_max
    integer  :: nx, ny, i, j, stat

    nx = size(u, 1)
    ny = size(u, 2)
    
    ! Input validation
    if (nx /= size(v,1) .or. ny /= size(v,2)) then
      print *, "ERROR [advance_explicit_step_uv]: u and v size mismatch"
      stop 1
    end if

    !---------------------------------------------------------
    ! CFL stability check
    !---------------------------------------------------------
    u_max = maxval(abs(u(2:nx-1, 2:ny-1)))
    v_max = maxval(abs(v(2:nx-1, 2:ny-1)))
    
    cfl_conv = dt_local * (u_max / dx + v_max / dy)
    cfl_diff = nu * dt_local * (1.0_dp / (dx*dx) + 1.0_dp / (dy*dy))
    cfl = cfl_conv + cfl_diff
    
    if (cfl > CFL_MAX_SAFE) then
      warning_count = warning_count + 1
      if (verbose .or. warning_count <= 5) then
        print '(A,F8.4,A,F8.4,A)', &
          " ⚠ CFL WARNING: CFL=", cfl, " > ", CFL_MAX_SAFE, &
          " - Reduce dt or refine grid!"
        if (warning_count == 5) then
          print *, " (Further CFL warnings suppressed)"
        end if
      end if
    end if

    !---------------------------------------------------------
    ! Allocate work arrays
    !---------------------------------------------------------
    allocate(lap_u(nx, ny), lap_v(nx, ny), stat=stat)
    if (stat /= 0) then
      print *, "ERROR [advance_explicit_step_uv]: Memory allocation failed"
      stop 1
    end if

    ! Compute Laplacians
    call laplacian_2d(u, dx, dy, lap_u)
    call laplacian_2d(v, dx, dy, lap_v)

    !---------------------------------------------------------
    ! Explicit update (advection + diffusion)
    !---------------------------------------------------------
    !$omp parallel do collapse(2) schedule(static) private(i,j)
    do j = 2, ny - 1
      do i = 2, nx - 1
        ! U-momentum: ∂u/∂t = -u·∇u + ν∇²u
        u(i,j) = u(i,j) + dt_local * &
                 (-u(i,j) * (u(i+1,j) - u(i-1,j)) / (2.0_dp * dx) - &
                  v(i,j) * (u(i,j+1) - u(i,j-1)) / (2.0_dp * dy) + &
                  nu * lap_u(i,j))
        
        ! V-momentum: ∂v/∂t = -u·∇v + ν∇²v
        v(i,j) = v(i,j) + dt_local * &
                 (-u(i,j) * (v(i+1,j) - v(i-1,j)) / (2.0_dp * dx) - &
                  v(i,j) * (v(i,j+1) - v(i,j-1)) / (2.0_dp * dy) + &
                  nu * lap_v(i,j))
      end do
    end do
    !$omp end parallel do

    deallocate(lap_u, lap_v)

  end subroutine advance_explicit_step_uv

  subroutine solve_dummy(u, v, p, dx, dy)
    !---------------------------------------------------------
    ! Fractional-step method (Chorin's projection)
    ! 1. Predictor: explicit advection-diffusion
    ! 2. Poisson solve for pressure
    ! 3. Corrector: project velocity to divergence-free
    !---------------------------------------------------------
    real(dp), intent(inout) :: u(:,:), v(:,:), p(:,:)
    real(dp), intent(in)    :: dx, dy
    
    real(dp), allocatable :: rhs(:,:)
    real(dp) :: div_max, div_l2
    integer  :: nx, ny, i, j

    nx = size(u, 1)
    ny = size(u, 2)
    
    ! Validate
    if (nx /= size(v,1) .or. nx /= size(p,1) .or. &
        ny /= size(v,2) .or. ny /= size(p,2)) then
      print *, "ERROR [solve_dummy]: Array size mismatch"
      stop 1
    end if
    
    allocate(rhs(nx, ny))

    !---------------------------------------------------------
    ! STEP 1: Apply boundary conditions
    !---------------------------------------------------------
    call apply_dirichlet_bc_scalar(u, Ulid, 0.0_dp, 0.0_dp, 0.0_dp)
    call apply_dirichlet_bc_scalar(v, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp)

    if (first_call .and. verbose) then
      print *, "============================================"
      print *, " Fractional-Step Solver Initialized"
      print *, "============================================"
      first_call = .false.
    end if

    !---------------------------------------------------------
    ! STEP 2: Predictor step (u*, v*)
    !---------------------------------------------------------
    call advance_explicit_step_uv(u, v, dx, dy, dt)
    
    ! Reapply BCs after advection
    call apply_dirichlet_bc_scalar(u, Ulid, 0.0_dp, 0.0_dp, 0.0_dp)
    call apply_dirichlet_bc_scalar(v, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp)

    !---------------------------------------------------------
    ! STEP 3: Pressure Poisson equation
    ! ∇²p = ρ/Δt ∇·u*
    !---------------------------------------------------------
    rhs = 0.0_dp
    !$omp parallel do collapse(2) schedule(static) private(i,j)
    do j = 2, ny - 1
      do i = 2, nx - 1
        rhs(i,j) = ((u(i+1,j) - u(i-1,j)) / (2.0_dp * dx) + &
                    (v(i,j+1) - v(i,j-1)) / (2.0_dp * dy)) / dt
      end do
    end do
    !$omp end parallel do
    
    ! Solve for pressure
    call solve_poisson_sor(p, rhs, dx, dy)
    
    ! Ensure pressure BCs (Neumann)
    call enforce_pressure_bc(p)

    !---------------------------------------------------------
    ! STEP 4: Corrector step (projection)
    ! u^(n+1) = u* - Δt ∇p
    !---------------------------------------------------------
    !$omp parallel do collapse(2) schedule(static) private(i,j)
    do j = 2, ny - 1
      do i = 2, nx - 1
        u(i,j) = u(i,j) - dt * (p(i+1,j) - p(i-1,j)) / (2.0_dp * dx)
        v(i,j) = v(i,j) - dt * (p(i,j+1) - p(i,j-1)) / (2.0_dp * dy)
      end do
    end do
    !$omp end parallel do
    
    ! Final BC enforcement
    call apply_dirichlet_bc_scalar(u, Ulid, 0.0_dp, 0.0_dp, 0.0_dp)
    call apply_dirichlet_bc_scalar(v, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp)

    !---------------------------------------------------------
    ! STEP 5: Diagnostics (check mass conservation)
    !---------------------------------------------------------
    if (verbose) then
      call compute_divergence_max(u, v, dx, dy, div_max, div_l2)
      if (div_max > DIV_TOLERANCE) then
        print '(A,ES10.2,A,ES10.2)', &
          " ⚠ Divergence: max=", div_max, "  L2=", div_l2
      end if
    end if

    deallocate(rhs)

  end subroutine solve_dummy

  !=============================================================
  ! HELPER FUNCTIONS
  !=============================================================

  function compute_residual(phi, rhs, dx, dy) result(res)
    !---------------------------------------------------------
    ! Compute L2 residual of Poisson equation
    !---------------------------------------------------------
    real(dp), intent(in) :: phi(:,:), rhs(:,:)
    real(dp), intent(in) :: dx, dy
    real(dp) :: res
    integer :: nx, ny, i, j
    real(dp) :: dx2, dy2, local_res

    nx = size(phi, 1)
    ny = size(phi, 2)
    dx2 = dx * dx
    dy2 = dy * dy
    
    res = 0.0_dp
    !$omp parallel do collapse(2) reduction(+:res) schedule(static) private(i,j,local_res)
    do j = 2, ny - 1
      do i = 2, nx - 1
        local_res = (phi(i+1,j) - 2.0_dp*phi(i,j) + phi(i-1,j)) / dx2 + &
                    (phi(i,j+1) - 2.0_dp*phi(i,j) + phi(i,j-1)) / dy2 - rhs(i,j)
        res = res + local_res * local_res
      end do
    end do
    !$omp end parallel do
    
    res = sqrt(res)

  end function compute_residual

  subroutine enforce_pressure_bc(p)
    !---------------------------------------------------------
    ! Neumann boundary conditions for pressure: ∂p/∂n = 0
    !---------------------------------------------------------
    real(dp), intent(inout) :: p(:,:)
    integer :: nx, ny

    nx = size(p, 1)
    ny = size(p, 2)

    ! West and East
    p(1,  :) = p(2,    :)
    p(nx, :) = p(nx-1, :)
    
    ! South and North
    p(:, 1)  = p(:, 2)
    p(:, ny) = p(:, ny-1)

  end subroutine enforce_pressure_bc

  subroutine compute_divergence_max(u, v, dx, dy, div_max, div_l2)
    !---------------------------------------------------------
    ! Compute divergence for incompressibility check
    !---------------------------------------------------------
    real(dp), intent(in)  :: u(:,:), v(:,:)
    real(dp), intent(in)  :: dx, dy
    real(dp), intent(out) :: div_max, div_l2
    
    integer  :: nx, ny, i, j
    real(dp) :: div_ij, div_sum

    nx = size(u, 1)
    ny = size(u, 2)
    
    div_max = 0.0_dp
    div_sum = 0.0_dp
    
    !$omp parallel do collapse(2) reduction(max:div_max) reduction(+:div_sum) &
    !$omp private(i,j,div_ij) schedule(static)
    do j = 2, ny - 1
      do i = 2, nx - 1
        div_ij = abs((u(i+1,j) - u(i-1,j)) / (2.0_dp * dx) + &
                     (v(i,j+1) - v(i,j-1)) / (2.0_dp * dy))
        div_max = max(div_max, div_ij)
        div_sum = div_sum + div_ij * div_ij
      end do
    end do
    !$omp end parallel do
    
    div_l2 = sqrt(div_sum / real((nx-2) * (ny-2), dp))

  end subroutine compute_divergence_max

  function check_cfl_stability(u, v, dx, dy, dt_local) result(is_stable)
    !---------------------------------------------------------
    ! Public function to check CFL before time step
    !---------------------------------------------------------
    real(dp), intent(in) :: u(:,:), v(:,:)
    real(dp), intent(in) :: dx, dy, dt_local
    logical :: is_stable
    real(dp) :: cfl, u_max, v_max
    
    u_max = maxval(abs(u))
    v_max = maxval(abs(v))
    
    cfl = dt_local * (u_max / dx + v_max / dy)
    is_stable = (cfl <= CFL_MAX_SAFE)
    
    if (.not. is_stable .and. verbose) then
      print '(A,F8.4,A)', " ✗ CFL check failed: CFL=", cfl, " > 0.5"
    end if

  end function check_cfl_stability

end module solver_mod
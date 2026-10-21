! Initial conditions s0 of the classical kicked top on the unit sphere.
! Lecture notes: Sec. 10 (initial conditions), eqs. (13) [eq:unit] and (9) [eq:area].
module init_mod
  use kinds_mod, only: wp, pi, two_pi
  use config_mod, only: config_t
  use sphere_mod, only: spherical_to_cartesian, equal_area_grid
  implicit none
  private
  public :: initial_conditions

contains

  ! Returns the initial points s0(3, n_ic) selected by cfg%mode:
  !   'single': one point (theta0, phi0); 'grid': equal-area grid; 'random': uniform random points.
  subroutine initial_conditions(cfg, s0)
    type(config_t), intent(in) :: cfg
    real(wp), allocatable, intent(out) :: s0(:,:)

    select case (trim(cfg%mode))
    case ('single')
      allocate(s0(3, 1))
      s0(:, 1) = spherical_to_cartesian(cfg%theta0, cfg%phi0)
    case ('grid')
      call grid_points(cfg%n_theta, cfg%n_phi, s0)
    case ('random')
      call random_points(cfg%n_random, cfg%seed, s0)
    end select
  end subroutine initial_conditions

  ! The cell centres of the equal-area grid of sphere_mod, as unit vectors.
  subroutine grid_points(n_theta, n_phi, s0)
    integer, intent(in) :: n_theta
    integer, intent(in) :: n_phi
    real(wp), allocatable, intent(out) :: s0(:,:)
    real(wp) :: theta(n_theta*n_phi), phi(n_theta*n_phi)
    integer :: ic

    call equal_area_grid(n_theta, n_phi, theta, phi)
    allocate(s0(3, n_theta*n_phi))
    do ic = 1, n_theta*n_phi
      s0(:, ic) = spherical_to_cartesian(theta(ic), phi(ic))
    end do
  end subroutine grid_points

  ! n_random points uniformly distributed on the sphere: Z and phi uniform,
  ! which by eq. (9) [eq:area] is the uniform distribution.
  subroutine random_points(n_random, seed, s0)
    integer, intent(in) :: n_random
    integer, intent(in) :: seed
    real(wp), allocatable, intent(out) :: s0(:,:)
    integer :: ic, n_seed, i
    integer, allocatable :: seed_arr(:)
    real(wp) :: r(2)

    call random_seed(size=n_seed)
    allocate(seed_arr(n_seed))
    seed_arr = [(seed + 37*i, i = 1, n_seed)]
    call random_seed(put=seed_arr)

    allocate(s0(3, n_random))
    do ic = 1, n_random
      call random_number(r)
      s0(:, ic) = spherical_to_cartesian(acos(2.0_wp*r(1) - 1.0_wp), -pi + two_pi*r(2))
    end do
  end subroutine random_points

end module init_mod

! Largest Lyapunov exponent of each initial condition.
! Lecture notes: Sec. 8 [sec:lyapunov], eq. (33) [eq:lyapunov] (definition) and
! Sec. 8.1 [sec:benettin] (Benettin algorithm, tangent map eq. (34) [eq:tangent]).
module lyapunov_mod
  use kinds_mod, only: wp
  use config_mod, only: config_t
  use map_mod, only: kick, kick_tangent
  implicit none
  private
  public :: lyapunov

contains

  ! Finite-time Lyapunov exponent Lambda_N, N = n_kicks, of every initial
  ! condition (Benettin algorithm, Sec. 8.1; OpenMP over initial conditions).
  ! After n_transient kicks, the state s and a tangent vector ds are evolved
  ! together for n_kicks kicks. Every n_renorm kicks ds is projected back onto
  ! the tangent plane, its growth ln|ds| is accumulated, and it is rescaled to
  ! unit length. lam(ic) = (sum of ln|ds|)/n_kicks, which equals eq. (33)
  ! [eq:lyapunov] evaluated at the finite time N.
  subroutine lyapunov(cfg, s0, lam)
    type(config_t), intent(in) :: cfg
    real(wp), intent(in) :: s0(:,:)
    real(wp), allocatable, intent(out) :: lam(:)
    integer :: ic, n
    real(wp) :: s(3), ds(3), log_sum

    allocate(lam(size(s0, 2)))

    !$omp parallel do private(s, ds, n, log_sum) schedule(dynamic)
    do ic = 1, size(s0, 2)
      s = s0(:, ic)
      do n = 1, cfg%n_transient
        s = kick(s, cfg%p, cfg%k)
      end do

      ds = initial_tangent(s)
      log_sum = 0.0_wp
      do n = 1, cfg%n_kicks
        ds = kick_tangent(s, ds, cfg%p, cfg%k)
        s = kick(s, cfg%p, cfg%k)
        if (mod(n, cfg%n_renorm) == 0 .or. n == cfg%n_kicks) call renormalise(s, ds, log_sum)
      end do
      lam(ic) = log_sum/real(cfg%n_kicks, wp)
    end do
    !$omp end parallel do
  end subroutine lyapunov

  ! A unit vector tangent to the sphere at s: z_hat x s normalised, or
  ! x_hat x s near the poles, where z_hat x s vanishes.
  pure function initial_tangent(s) result(ds)
    real(wp), intent(in) :: s(3)
    real(wp) :: ds(3)

    if (s(1)**2 + s(2)**2 > 1.0e-2_wp) then
      ds = [-s(2), s(1), 0.0_wp]
    else
      ds = [0.0_wp, -s(3), s(2)]
    end if
    ds = ds/norm2(ds)
  end function initial_tangent

  ! Removes the round-off component of ds along s, adds ln|ds| to log_sum
  ! and rescales ds to unit length.
  pure subroutine renormalise(s, ds, log_sum)
    real(wp), intent(in) :: s(3)
    real(wp), intent(inout) :: ds(3)
    real(wp), intent(inout) :: log_sum
    real(wp) :: d

    ds = ds - dot_product(ds, s)*s
    d = norm2(ds)
    log_sum = log_sum + log(d)
    ds = ds/d
  end subroutine renormalise

end module lyapunov_mod

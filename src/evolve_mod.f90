! Iteration of the stroboscopic map, s_n = F^n(s_0), for all initial conditions.
! Lecture notes: Sec. 3 [sec:map] (the map) and Sec. 10 (recording, OpenMP).
module evolve_mod
  use kinds_mod, only: wp
  use config_mod, only: config_t
  use map_mod, only: kick
  implicit none
  private
  public :: evolve

contains

  ! Iterates kick, eq. (20) [eq:map], for every initial condition (in parallel
  ! over initial conditions with OpenMP). After n_transient discarded
  ! kicks, the state is stored at kicks 0, save_every, 2*save_every, ... <= n_kicks.
  ! orbit has shape (3, n_saved, n_ic).
  subroutine evolve(cfg, s0, orbit)
    type(config_t), intent(in) :: cfg
    real(wp), intent(in) :: s0(:,:)
    real(wp), allocatable, intent(out) :: orbit(:,:,:)
    integer :: ic, n, isave, n_saved
    real(wp) :: s(3)

    n_saved = cfg%n_kicks/cfg%save_every + 1
    allocate(orbit(3, n_saved, size(s0, 2)))

    !$omp parallel do private(s, n, isave) schedule(static)
    do ic = 1, size(s0, 2)
      s = s0(:, ic)
      do n = 1, cfg%n_transient
        s = kick(s, cfg%p, cfg%k)
      end do

      isave = 1
      orbit(:, isave, ic) = s
      do n = 1, cfg%n_kicks
        s = kick(s, cfg%p, cfg%k)
        if (mod(n, cfg%save_every) == 0) then
          isave = isave + 1
          orbit(:, isave, ic) = s
        end if
      end do
    end do
    !$omp end parallel do
  end subroutine evolve

end module evolve_mod

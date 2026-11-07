! Driver of the three-regime study: for every rotation angle h = p in the
! config it runs
!   S1  the classical kicked top (Poincare section on a grid and the orbit of
!       the initial state), with the torsion k = 2 n_q eps of the chain;
!   S3  the qubit chain, n_q spins 1/2 with all-to-all zz couplings, evolved by
!       its Floquet step on the state vector only, with Husimi snapshots (S2);
!   S4  time reversal of the chain: t_r forward steps with T = 2 pi + eps,
!       then t_r backward steps with T' = 2 pi - eps and -h.
! Usage: kicked_top_regimes [namelist file]   (default config_regimes.nml)
program kicked_top_regimes
  use kinds_mod, only: wp
  use io_utils_mod, only: make_dir
  use sphere_mod, only: equal_area_grid
  use config_mod, only: config_t
  use init_mod, only: initial_conditions
  use evolve_mod, only: evolve
  use io_mod, only: write_results
  use husimi_mod, only: husimi, husimi_norm
  use regimes_config_mod, only: rconfig_t, read_regimes_config, str_len
  use regimes_io_mod, only: run_dir, write_regimes, write_chain, write_husimi, write_summary
  use chain_state_mod, only: ferromagnetic_direction, product_coherent_state, dicke_amplitudes
  use chain_floquet_mod, only: forward_step, backward_step
  use chain_observables_mod, only: fidelity, site_spins, total_spin
  implicit none

  type(rconfig_t) :: cfg
  character(len=str_len) :: cfg_path
  real(wp) :: theta0, phi0
  real(wp), allocatable :: hus_theta(:), hus_phi(:)
  integer :: i

  cfg_path = 'config_regimes.nml'
  if (command_argument_count() >= 1) call get_command_argument(1, cfg_path)

  call read_regimes_config(trim(cfg_path), cfg)
  call ferromagnetic_direction(cfg%s0, cfg%n_q, cfg%pol, theta0, phi0)
  call make_dir(trim(cfg%out_dir))
  call write_regimes(cfg)

  allocate(hus_theta(cfg%hus_n_theta*cfg%hus_n_phi), hus_phi(cfg%hus_n_theta*cfg%hus_n_phi))
  call equal_area_grid(cfg%hus_n_theta, cfg%hus_n_phi, hus_theta, hus_phi)

  do i = 1, size(cfg%h)
    call run_classical(cfg%h(i), run_dir(cfg, i)//'/classical')
    call run_chain(cfg%h(i), run_dir(cfg, i)//'/chain')
  end do

contains

  ! S1: classical map, lecture notes eq. (20) [eq:map], with p = h and k = 2 n_q eps.
  subroutine run_classical(h, dir)
    real(wp), intent(in) :: h
    character(len=*), intent(in) :: dir
    type(config_t) :: c
    real(wp), allocatable :: s0(:,:), orbit(:,:,:)

    c%p = h
    c%k = cfg%k
    c%n_kicks = cfg%n_kicks
    c%n_transient = 0
    c%save_every = 1
    c%theta0 = theta0
    c%phi0 = phi0
    c%n_theta = cfg%n_theta
    c%n_phi = cfg%n_phi
    c%n_random = 1
    c%seed = 1
    c%lyapunov = .false.
    c%n_renorm = 1

    c%mode = 'grid'
    c%out_dir = dir//'/section'
    call initial_conditions(c, s0)
    call evolve(c, s0, orbit)
    call write_results(c, orbit)

    c%mode = 'single'
    c%out_dir = dir//'/orbit'
    call initial_conditions(c, s0)
    call evolve(c, s0, orbit)
    call write_results(c, orbit)
  end subroutine run_classical

  ! S3/S4: forward evolution for t_r steps, then backward evolution for t_r steps.
  subroutine run_chain(h, dir)
    real(wp), intent(in) :: h
    character(len=*), intent(in) :: dir
    complex(wp), allocatable :: psi0(:), psi(:)
    real(wp), allocatable :: obs(:,:), q(:,:)
    real(wp) :: sig(3, cfg%n_q), j, sym_dev
    integer :: steps(3), t, l

    j = 0.5_wp*cfg%n_q
    steps = [0, cfg%t_r, 2*cfg%t_r]
    allocate(obs(9 + 3*cfg%n_q, 0:2*cfg%t_r), q(size(hus_theta), size(steps)))
    psi0 = product_coherent_state(cfg%n_q, theta0, phi0)
    psi = psi0
    sym_dev = 0.0_wp

    do t = 0, 2*cfg%t_r
      if (t >= 1 .and. t <= cfg%t_r) call forward_step(psi, cfg%n_q, h, cfg%eps)
      if (t > cfg%t_r) call backward_step(psi, cfg%n_q, h, cfg%eps)

      sig = site_spins(psi, cfg%n_q)
      obs(1, t) = fidelity(psi0, psi)
      obs(2, t) = sum(abs(psi)**2)
      obs(3:9, t) = total_spin(psi, cfg%n_q)
      obs(10:, t) = [sig(1, :), sig(2, :), sig(3, :)]

      do l = 1, size(steps)
        if (t == steps(l)) then
          sym_dev = max(sym_dev, abs(1.0_wp - sum(abs(dicke_amplitudes(psi, cfg%n_q))**2)))
          call husimi(dicke_amplitudes(psi, cfg%n_q), j, hus_theta, hus_phi, q(:, l))
        end if
      end do
    end do

    call make_dir(dir)
    call write_chain(cfg, h, dir, obs)
    call write_husimi(cfg, h, dir, steps, hus_theta, hus_phi, q)
    call write_summary(cfg, h, dir, &
      [character(len=24) :: 'F(t_r)', 'F(2t_r)', 'max|norm-1|', 'max|S^2-j(j+1)|', &
       'max(1-symmetric_weight)', 'husimi_norm(0)', 'husimi_norm(t_r)', 'husimi_norm(2t_r)'], &
      [obs(1, cfg%t_r), obs(1, 2*cfg%t_r), maxval(abs(obs(2, :) - 1.0_wp)), &
       maxval(abs(obs(9, :) - j*(j + 1.0_wp))), sym_dev, &
       husimi_norm(q(:, 1), j), husimi_norm(q(:, 2), j), husimi_norm(q(:, 3), j)])

    write(*, '(a,f0.4,a,f0.4,a,es10.3,a)') 'kicked_top_regimes: h = ', h, ', k = ', cfg%k, &
      ', F(2t_r) - 1 = ', obs(1, 2*cfg%t_r) - 1.0_wp, ' -> '//dir
  end subroutine run_chain

end program kicked_top_regimes

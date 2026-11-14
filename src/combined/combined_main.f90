! Driver of the three-regime study: for every rotation angle p in the
! config it runs
! (make combined; make spin_chain runs the chain part alone, kicked_top_spin_chain)
!   S1  the classical kicked top (Poincare section on a grid and the orbit of
!       the initial state), with the torsion k of the chain;
!   S3  the qubit chain, n_q spins 1/2 with all-to-all zz couplings V_il = 1/4,
!       evolved by its Floquet step on the state vector only, with Husimi
!       snapshots (S2);
!   S4  time reversal of the chain: t_r forward steps U, then t_r backward
!       steps U^dagger (coupling 2 pi - k/(2j), kick -p); with butterfly_effect
!       qubit i_flip is rotated about a = flip_axis by flip_angle
!       (exp(-i flip_angle/2 sigma^a), pi = flip sigma^a) at t_r, before the backward leg.
! Usage: kicked_top_combined [namelist file]   (default config_combined.nml)
program kicked_top_combined
  use kinds_mod, only: wp
  use io_utils_mod, only: make_dir
  use sphere_mod, only: equal_area_grid
  use config_mod, only: config_t
  use init_mod, only: initial_conditions
  use evolve_mod, only: evolve
  use io_mod, only: write_results
  use combined_config_mod, only: rconfig_t, read_combined_config, str_len
  use combined_io_mod, only: run_dir, write_regimes
  use chain_state_mod, only: ferromagnetic_direction
  use chain_floquet_mod, only: uniform_coupling, zz_energies
  use chain_run_mod, only: run_chain
  implicit none

  type(rconfig_t) :: cfg
  character(len=str_len) :: cfg_path
  real(wp) :: theta0, phi0
  real(wp), allocatable :: hus_theta(:), hus_phi(:), e_zz(:)
  integer :: i

  cfg_path = 'config_combined.nml'
  if (command_argument_count() >= 1) call get_command_argument(1, cfg_path)

  call read_combined_config(trim(cfg_path), cfg)
  call ferromagnetic_direction(cfg%s0, cfg%n_q, cfg%pol, theta0, phi0)
  call make_dir(trim(cfg%out_dir))
  call write_regimes(cfg)

  allocate(hus_theta(cfg%hus_n_theta*cfg%hus_n_phi), hus_phi(cfg%hus_n_theta*cfg%hus_n_phi))
  call equal_area_grid(cfg%hus_n_theta, cfg%hus_n_phi, hus_theta, hus_phi)

  e_zz = zz_energies(cfg%n_q, uniform_coupling(cfg%n_q))

  do i = 1, size(cfg%p)
    call run_classical(cfg%p(i), run_dir(cfg, i)//'/classical')
    call run_chain(cfg, cfg%p(i), theta0, phi0, e_zz, hus_theta, hus_phi, run_dir(cfg, i)//'/chain', 'kicked_top_combined')
  end do

contains

  ! S1: classical map, lecture notes eq. (20) [eq:map], with the p and k of the chain.
  subroutine run_classical(p, dir)
    real(wp), intent(in) :: p
    character(len=*), intent(in) :: dir
    type(config_t) :: c
    real(wp), allocatable :: s0(:,:), orbit(:,:,:)

    c%p = p
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

end program kicked_top_combined

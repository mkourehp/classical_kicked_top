! Driver of the qubit chain alone (make spin_chain): for every rotation angle p
! in the config, n_q spins 1/2 with all-to-all zz couplings V_il = 1/4 are
! evolved t_r steps forward, then (with butterfly_effect) sigma^a, a = flip_axis,
! acts on qubit i_flip, then t_r steps backward (coupling 2 pi - k/(2j), kick -p).
! Same namelist and chain output as kicked_top_combined, without the classical part.
! Usage: kicked_top_spin_chain [namelist file]   (default config_spin_chain.nml)
program kicked_top_spin_chain
  use kinds_mod, only: wp
  use io_utils_mod, only: make_dir
  use sphere_mod, only: equal_area_grid
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

  cfg_path = 'config_spin_chain.nml'
  if (command_argument_count() >= 1) call get_command_argument(1, cfg_path)

  call read_combined_config(trim(cfg_path), cfg)
  call ferromagnetic_direction(cfg%s0, cfg%n_q, cfg%pol, theta0, phi0)
  call make_dir(trim(cfg%out_dir))
  call write_regimes(cfg)

  allocate(hus_theta(cfg%hus_n_theta*cfg%hus_n_phi), hus_phi(cfg%hus_n_theta*cfg%hus_n_phi))
  call equal_area_grid(cfg%hus_n_theta, cfg%hus_n_phi, hus_theta, hus_phi)

  e_zz = zz_energies(cfg%n_q, uniform_coupling(cfg%n_q))

  do i = 1, size(cfg%p)
    call run_chain(cfg, cfg%p(i), theta0, phi0, e_zz, hus_theta, hus_phi, run_dir(cfg, i)//'/chain', &
      'kicked_top_spin_chain')
  end do

end program kicked_top_spin_chain

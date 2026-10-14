program kicked_top
  use kinds_mod, only: wp
  use config_mod, only: config_t, read_config, str_len
  use init_mod, only: initial_conditions
  use evolve_mod, only: evolve
  use io_mod, only: write_results
  implicit none

  type(config_t) :: cfg
  character(len=str_len) :: cfg_path
  real(wp), allocatable :: s0(:,:)
  real(wp), allocatable :: orbit(:,:,:)

  cfg_path = 'config.nml'
  if (command_argument_count() >= 1) call get_command_argument(1, cfg_path)

  call read_config(trim(cfg_path), cfg)
  call initial_conditions(cfg, s0)
  call evolve(cfg, s0, orbit)
  call write_results(cfg, orbit)

  write(*, '(a,i0,a,i0,a)') 'kicked_top: ', size(orbit, 3), ' orbits, ', size(orbit, 2), &
    ' points each -> '//trim(cfg%out_dir)
end program kicked_top

! Driver of the classical kicked top.
! Usage: kicked_top [namelist file]   (default config.nml)
! read_config -> initial_conditions -> evolve -> write_results, and, if
! &lyapunov is enabled, lyapunov -> write_lyapunov.
! Lecture notes: Sec. 11 [sec:code], Fig. 10 [fig:flow].
program kicked_top
  use kinds_mod, only: wp
  use config_mod, only: config_t, read_config, str_len
  use init_mod, only: initial_conditions
  use evolve_mod, only: evolve
  use lyapunov_mod, only: lyapunov
  use io_mod, only: write_results, write_lyapunov
  implicit none

  type(config_t) :: cfg
  character(len=str_len) :: cfg_path
  real(wp), allocatable :: s0(:,:)
  real(wp), allocatable :: orbit(:,:,:)
  real(wp), allocatable :: lam(:)

  cfg_path = 'config.nml'
  if (command_argument_count() >= 1) call get_command_argument(1, cfg_path)

  call read_config(trim(cfg_path), cfg)
  call initial_conditions(cfg, s0)
  call evolve(cfg, s0, orbit)
  call write_results(cfg, orbit)
  if (cfg%lyapunov) then
    call lyapunov(cfg, s0, lam)
    call write_lyapunov(cfg, s0, lam)
  end if

  write(*, '(a,i0,a,i0,a)') 'kicked_top: ', size(orbit, 3), ' orbits, ', size(orbit, 2), &
    ' points each -> '//trim(cfg%out_dir)
end program kicked_top

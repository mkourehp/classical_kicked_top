module config_mod
  use kinds_mod, only: wp
  implicit none
  private
  public :: config_t, read_config, str_len

  integer, parameter :: str_len = 256

  type :: config_t
    ! &system
    real(wp) :: p
    real(wp) :: k
    ! &simulation
    integer :: n_kicks
    integer :: n_transient
    integer :: save_every
    ! &initial
    character(len=str_len) :: mode
    real(wp) :: theta0
    real(wp) :: phi0
    integer :: n_theta
    integer :: n_phi
    integer :: n_random
    integer :: seed
    ! &output
    character(len=str_len) :: out_dir
  end type config_t

contains

  subroutine read_config(path, cfg)
    character(len=*), intent(in) :: path
    type(config_t), intent(out) :: cfg

    real(wp) :: p, k, theta0, phi0
    integer :: n_kicks, n_transient, save_every
    integer :: n_theta, n_phi, n_random, seed
    character(len=str_len) :: mode, out_dir
    integer :: u, ios

    namelist /system/ p, k
    namelist /simulation/ n_kicks, n_transient, save_every
    namelist /initial/ mode, theta0, phi0, n_theta, n_phi, n_random, seed
    namelist /output/ out_dir

    p = acos(-1.0_wp)/2.0_wp
    k = 3.0_wp
    n_kicks = 500
    n_transient = 0
    save_every = 1
    mode = 'grid'
    theta0 = 1.0_wp
    phi0 = 0.5_wp
    n_theta = 20
    n_phi = 20
    n_random = 100
    seed = 12345
    out_dir = 'output'

    open(newunit=u, file=path, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'config_mod: cannot open namelist file '//trim(path)

    read(u, nml=system, iostat=ios)
    call check_read(ios, 'system')
    rewind(u)
    read(u, nml=simulation, iostat=ios)
    call check_read(ios, 'simulation')
    rewind(u)
    read(u, nml=initial, iostat=ios)
    call check_read(ios, 'initial')
    rewind(u)
    read(u, nml=output, iostat=ios)
    call check_read(ios, 'output')
    close(u)

    cfg%p = p
    cfg%k = k
    cfg%n_kicks = n_kicks
    cfg%n_transient = n_transient
    cfg%save_every = save_every
    cfg%mode = mode
    cfg%theta0 = theta0
    cfg%phi0 = phi0
    cfg%n_theta = n_theta
    cfg%n_phi = n_phi
    cfg%n_random = n_random
    cfg%seed = seed
    cfg%out_dir = out_dir

    call validate(cfg)
  end subroutine read_config

  subroutine check_read(ios, group)
    integer, intent(in) :: ios
    character(len=*), intent(in) :: group

    ! A missing group (end of file) keeps the defaults.
    if (ios /= 0 .and. .not. is_iostat_end(ios)) then
      error stop 'config_mod: error reading namelist group &'//group
    end if
  end subroutine check_read

  subroutine validate(cfg)
    type(config_t), intent(in) :: cfg

    if (cfg%n_kicks < 1) error stop 'config_mod: n_kicks must be >= 1'
    if (cfg%n_transient < 0) error stop 'config_mod: n_transient must be >= 0'
    if (cfg%save_every < 1) error stop 'config_mod: save_every must be >= 1'

    select case (trim(cfg%mode))
    case ('single')
    case ('grid')
      if (cfg%n_theta < 1 .or. cfg%n_phi < 1) error stop 'config_mod: n_theta and n_phi must be >= 1'
    case ('random')
      if (cfg%n_random < 1) error stop 'config_mod: n_random must be >= 1'
    case default
      error stop "config_mod: mode must be 'single', 'grid' or 'random'"
    end select
  end subroutine validate

end module config_mod

! Run configuration: the derived type config_t and the reader of the namelist
! file (default config.nml). Each namelist group has its own block of fields;
! a missing group or variable keeps its default value.
! Lecture notes: Sec. 11.1 (table of parameters); p and k are the parameters
! of the Hamiltonian, eq. (12) [eq:H].
module config_mod
  use kinds_mod, only: wp
  implicit none
  private
  public :: config_t, read_config, str_len

  integer, parameter :: str_len = 256   ! length of string parameters

  type :: config_t
    ! &system
    real(wp) :: p              ! rotation angle per period, eq. (12)
    real(wp) :: k              ! torsion (kick) strength, eq. (12)
    ! &simulation
    integer :: n_kicks         ! kicks recorded after the transient
    integer :: n_transient     ! kicks discarded before recording
    integer :: save_every      ! store every save_every-th kick
    ! &initial
    character(len=str_len) :: mode   ! 'single' | 'grid' | 'random'
    real(wp) :: theta0         ! single: polar angle of the initial point
    real(wp) :: phi0           ! single: azimuth of the initial point
    integer :: n_theta         ! grid: number of points in Z = cos(theta)
    integer :: n_phi           ! grid: number of points in phi
    integer :: n_random        ! random: number of points
    integer :: seed            ! random: seed of the random number generator
    ! &output
    character(len=str_len) :: out_dir   ! directory of the output files
    ! &lyapunov
    logical :: lyapunov        ! compute the Lyapunov exponent, Sec. 8.1 [sec:benettin]
    integer :: n_renorm        ! rescale the tangent vector every n_renorm kicks
  end type config_t

contains

  ! Reads all namelist groups from the file at path into cfg and validates them.
  subroutine read_config(path, cfg)
    character(len=*), intent(in) :: path
    type(config_t), intent(out) :: cfg

    real(wp) :: p, k, theta0, phi0
    integer :: n_kicks, n_transient, save_every
    integer :: n_theta, n_phi, n_random, seed
    character(len=str_len) :: mode, out_dir
    logical :: enabled
    integer :: n_renorm
    integer :: u, ios

    namelist /system/ p, k
    namelist /simulation/ n_kicks, n_transient, save_every
    namelist /initial/ mode, theta0, phi0, n_theta, n_phi, n_random, seed
    namelist /output/ out_dir
    namelist /lyapunov/ enabled, n_renorm

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
    enabled = .false.
    n_renorm = 1

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
    rewind(u)
    read(u, nml=lyapunov, iostat=ios)
    call check_read(ios, 'lyapunov')
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
    cfg%lyapunov = enabled
    cfg%n_renorm = n_renorm

    call validate(cfg)
  end subroutine read_config

  ! Stops on a malformed namelist group.
  subroutine check_read(ios, group)
    integer, intent(in) :: ios
    character(len=*), intent(in) :: group

    ! A missing group (end of file) keeps the defaults.
    if (ios /= 0 .and. .not. is_iostat_end(ios)) then
      error stop 'config_mod: error reading namelist group &'//group
    end if
  end subroutine check_read

  ! Stops with a message if a parameter is out of range.
  subroutine validate(cfg)
    type(config_t), intent(in) :: cfg

    if (cfg%n_kicks < 1) error stop 'config_mod: n_kicks must be >= 1'
    if (cfg%n_transient < 0) error stop 'config_mod: n_transient must be >= 0'
    if (cfg%save_every < 1) error stop 'config_mod: save_every must be >= 1'
    if (cfg%n_renorm < 1) error stop 'config_mod: n_renorm must be >= 1'

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

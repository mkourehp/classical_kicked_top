! Configuration of kicked_top_regimes, read from a namelist file (default
! config_regimes.nml). One parameter set drives the classical kicked top and
! the qubit chain for every rotation angle h = p in the list h(:).
module regimes_config_mod
  use kinds_mod, only: wp
  use io_utils_mod, only: check_nml_read
  implicit none
  private
  public :: rconfig_t, read_regimes_config, str_len

  integer, parameter :: str_len = 256   ! length of string parameters
  integer, parameter :: max_h = 20      ! maximum number of h values

  type :: rconfig_t
    integer :: n_q                 ! number of qubits, j = n_q/2
    real(wp) :: eps                ! coupling time offset: T = 2 pi + eps, T' = 2 pi - eps
    real(wp) :: k                  ! torsion strength of the kicked top, k = 2 n_q eps
    real(wp), allocatable :: h(:)  ! rotation angles h = p, one run per value
    integer :: s0                  ! initial state: -1 all up, 0 all down (along pol)
    character(len=1) :: pol        ! axis of the initial polarisation: 'z', 'x', 'y'
    integer :: n_kicks             ! classical: kicks per orbit
    integer :: n_theta             ! classical: grid points in Z = cos(theta)
    integer :: n_phi               ! classical: grid points in phi
    integer :: t_r                 ! chain: forward steps before the reversal
    integer :: hus_n_theta         ! Husimi grid points in Z
    integer :: hus_n_phi           ! Husimi grid points in phi
    character(len=str_len) :: out_dir   ! output directory
  end type rconfig_t

contains

  subroutine read_regimes_config(path, cfg)
    character(len=*), intent(in) :: path
    type(rconfig_t), intent(out) :: cfg

    integer :: n_q, s0, n_kicks, n_theta, n_phi, t_r
    real(wp) :: eps, h(max_h)
    character(len=1) :: pol
    character(len=str_len) :: out_dir
    integer :: u, ios

    namelist /system/ n_q, eps, h
    namelist /initial/ s0, pol
    namelist /classical/ n_kicks, n_theta, n_phi
    namelist /chain/ t_r
    namelist /husimi/ n_theta, n_phi
    namelist /output/ out_dir

    n_q = 10
    eps = 1.0_wp
    h = huge(1.0_wp)
    h(1:3) = [0.02_wp, 0.1_wp, 1.0_wp]
    s0 = -1
    pol = 'z'
    t_r = 30
    out_dir = 'output_regimes'

    open(newunit=u, file=path, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'regimes_config_mod: cannot open namelist file '//trim(path)

    read(u, nml=system, iostat=ios)
    call check_nml_read(ios, 'system', 'regimes_config_mod')
    rewind(u)
    read(u, nml=initial, iostat=ios)
    call check_nml_read(ios, 'initial', 'regimes_config_mod')
    rewind(u)
    n_kicks = 500
    n_theta = 20
    n_phi = 20
    read(u, nml=classical, iostat=ios)
    call check_nml_read(ios, 'classical', 'regimes_config_mod')
    cfg%n_kicks = n_kicks
    cfg%n_theta = n_theta
    cfg%n_phi = n_phi
    rewind(u)
    read(u, nml=chain, iostat=ios)
    call check_nml_read(ios, 'chain', 'regimes_config_mod')
    rewind(u)
    n_theta = 100
    n_phi = 200
    read(u, nml=husimi, iostat=ios)
    call check_nml_read(ios, 'husimi', 'regimes_config_mod')
    cfg%hus_n_theta = n_theta
    cfg%hus_n_phi = n_phi
    rewind(u)
    read(u, nml=output, iostat=ios)
    call check_nml_read(ios, 'output', 'regimes_config_mod')
    close(u)

    cfg%n_q = n_q
    cfg%eps = eps
    cfg%k = 2.0_wp*n_q*eps
    cfg%h = pack(h, h < huge(1.0_wp))
    cfg%s0 = s0
    cfg%pol = pol
    cfg%t_r = t_r
    cfg%out_dir = out_dir

    call validate(cfg)
  end subroutine read_regimes_config

  subroutine validate(cfg)
    type(rconfig_t), intent(in) :: cfg

    if (cfg%n_q < 1 .or. cfg%n_q > 26) error stop 'regimes_config_mod: n_q must be in 1 .. 26'
    if (size(cfg%h) < 1) error stop 'regimes_config_mod: give at least one value of h'
    if (cfg%n_kicks < 1) error stop 'regimes_config_mod: n_kicks must be >= 1'
    if (cfg%n_theta < 1 .or. cfg%n_phi < 1) error stop 'regimes_config_mod: n_theta and n_phi must be >= 1'
    if (cfg%t_r < 1) error stop 'regimes_config_mod: t_r must be >= 1'
    if (cfg%hus_n_theta < 1 .or. cfg%hus_n_phi < 1) error stop 'regimes_config_mod: Husimi grid must be >= 1'
  end subroutine validate

end module regimes_config_mod

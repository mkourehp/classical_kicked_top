! Run configuration of the quantum kicked top: the derived type qconfig_t and the
! reader of its namelist file (default config_quantum.nml). A missing group or
! variable keeps its default value.
module quantum_config_mod
  use kinds_mod, only: wp, pi
  use io_utils_mod, only: check_nml_read
  implicit none
  private
  public :: qconfig_t, read_qconfig

  integer, parameter :: str_len = 256
  integer, parameter :: max_snapshots = 100

  type :: qconfig_t
    ! &system
    real(wp) :: p                  ! rotation angle per period
    real(wp) :: k                  ! torsion (kick) strength
    real(wp) :: j                  ! spin quantum number (integer or half-integer)
    integer :: dim                 ! dimension 2j+1 of the Hilbert space
    ! &simulation
    integer :: n_kicks             ! number of kicks
    ! &initial
    real(wp) :: theta0             ! polar angle of the initial spin coherent state
    real(wp) :: phi0               ! azimuth of the initial spin coherent state
    ! &husimi
    logical :: husimi              ! write Husimi snapshots
    integer :: n_theta             ! grid points in Z = cos(theta)
    integer :: n_phi               ! grid points in phi
    integer, allocatable :: husimi_kicks(:)   ! kicks at which Q is written
    ! &spectrum
    logical :: spectrum            ! compute quasienergies and spacing statistics
    integer :: n_bins              ! histogram bins on [0, s_max]
    real(wp) :: s_max              ! largest unfolded spacing in the histogram
    ! &output
    character(len=str_len) :: out_dir   ! directory of the output files
  end type qconfig_t

contains

  ! Reads all namelist groups from the file at path into cfg and validates them.
  subroutine read_qconfig(path, cfg)
    character(len=*), intent(in) :: path
    type(qconfig_t), intent(out) :: cfg

    real(wp) :: p, k, j, theta0, phi0, s_max
    integer :: n_kicks, n_theta, n_phi, n_bins, kicks(max_snapshots)
    logical :: enabled
    character(len=str_len) :: out_dir
    integer :: u, ios, two_j

    namelist /system/ p, k, j
    namelist /simulation/ n_kicks
    namelist /initial/ theta0, phi0
    namelist /husimi/ enabled, n_theta, n_phi, kicks
    namelist /spectrum/ enabled, n_bins, s_max
    namelist /output/ out_dir

    p = pi/2.0_wp
    k = 3.0_wp
    j = 50.0_wp
    n_kicks = 100
    theta0 = 2.25_wp
    phi0 = 0.63_wp
    n_theta = 100
    n_phi = 200
    kicks = -1
    kicks(1:4) = [0, 10, 50, 100]
    n_bins = 30
    s_max = 3.0_wp
    out_dir = 'output_quantum'

    open(newunit=u, file=path, status='old', action='read', iostat=ios)
    if (ios /= 0) error stop 'quantum_config_mod: cannot open namelist file '//trim(path)

    read(u, nml=system, iostat=ios)
    call check_nml_read(ios, 'system', 'quantum_config_mod')
    rewind(u)
    read(u, nml=simulation, iostat=ios)
    call check_nml_read(ios, 'simulation', 'quantum_config_mod')
    rewind(u)
    read(u, nml=initial, iostat=ios)
    call check_nml_read(ios, 'initial', 'quantum_config_mod')
    rewind(u)
    enabled = .true.
    read(u, nml=husimi, iostat=ios)
    call check_nml_read(ios, 'husimi', 'quantum_config_mod')
    cfg%husimi = enabled
    rewind(u)
    enabled = .true.
    read(u, nml=spectrum, iostat=ios)
    call check_nml_read(ios, 'spectrum', 'quantum_config_mod')
    cfg%spectrum = enabled
    rewind(u)
    read(u, nml=output, iostat=ios)
    call check_nml_read(ios, 'output', 'quantum_config_mod')
    close(u)

    two_j = nint(2.0_wp*j)
    if (two_j < 1 .or. abs(2.0_wp*j - real(two_j, wp)) > 1.0e-9_wp) &
      error stop 'quantum_config_mod: j must be a positive integer or half-integer'

    cfg%p = p
    cfg%k = k
    cfg%j = 0.5_wp*real(two_j, wp)
    cfg%dim = two_j + 1
    cfg%n_kicks = n_kicks
    cfg%theta0 = theta0
    cfg%phi0 = phi0
    cfg%n_theta = n_theta
    cfg%n_phi = n_phi
    cfg%husimi_kicks = pack(kicks, kicks >= 0 .and. kicks <= n_kicks)
    cfg%n_bins = n_bins
    cfg%s_max = s_max
    cfg%out_dir = out_dir

    call validate(cfg)
  end subroutine read_qconfig

  ! Stops with a message if a parameter is out of range.
  subroutine validate(cfg)
    type(qconfig_t), intent(in) :: cfg

    if (cfg%n_kicks < 0) error stop 'quantum_config_mod: n_kicks must be >= 0'
    if (cfg%theta0 < 0.0_wp .or. cfg%theta0 > pi) error stop 'quantum_config_mod: theta0 must be in [0, pi]'
    if (cfg%husimi .and. (cfg%n_theta < 1 .or. cfg%n_phi < 1)) &
      error stop 'quantum_config_mod: n_theta and n_phi must be >= 1'
    if (cfg%spectrum .and. (cfg%n_bins < 1 .or. cfg%s_max <= 0.0_wp)) &
      error stop 'quantum_config_mod: n_bins must be >= 1 and s_max > 0'
  end subroutine validate

end module quantum_config_mod

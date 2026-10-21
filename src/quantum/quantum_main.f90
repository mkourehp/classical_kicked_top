! Driver of the quantum kicked top.
! Usage: kicked_top_quantum [namelist file]   (default config_quantum.nml)
! read_qconfig -> Floquet operator U -> evolve a spin coherent state kick by kick
! (expectation values, Husimi snapshots) -> quasienergies and spacing statistics
! -> summary of the correctness checks.
program kicked_top_quantum
  use kinds_mod, only: wp, pi
  use io_utils_mod, only: make_dir
  use quantum_config_mod, only: qconfig_t, read_qconfig
  use spin_mod, only: coherent_state, expectation_j
  use floquet_mod, only: diagonalise_jy, exp_jy, floquet_operator
  use sphere_mod, only: equal_area_grid
  use husimi_mod, only: husimi, husimi_norm
  use spectrum_mod, only: quasienergies, spacing_statistics, commutator_error, unitarity_error
  use quantum_io_mod, only: write_expectation, write_husimi, write_quasienergies, write_spacings, &
                            write_summary
  implicit none

  type(qconfig_t) :: cfg
  character(len=256) :: cfg_path
  complex(wp), allocatable :: v(:,:), u(:,:), psi(:)
  real(wp), allocatable :: mu(:), obs(:,:), theta(:), phi(:), q(:,:)
  real(wp), allocatable :: eps(:), s_centre(:), density(:)
  integer, allocatable :: sector(:)
  character(len=40), allocatable :: names(:)
  real(wp), allocatable :: values(:)
  character(len=40) :: name
  real(wp) :: mean_r, max_dev
  integer :: n, i

  cfg_path = 'config_quantum.nml'
  if (command_argument_count() >= 1) call get_command_argument(1, cfg_path)
  call read_qconfig(trim(cfg_path), cfg)
  call make_dir(cfg%out_dir)

  allocate(v(cfg%dim, cfg%dim), mu(cfg%dim))
  call diagonalise_jy(cfg%j, cfg%dim, v, mu)
  u = floquet_operator(cfg%j, cfg%p, cfg%k, v, mu)

  allocate(names(0), values(0))
  call add('unitarity_error', unitarity_error(u))
  call add('parity_commutator', commutator_error(u, exp_jy(v, mu, pi)))

  if (cfg%husimi) then
    allocate(theta(cfg%n_theta*cfg%n_phi), phi(cfg%n_theta*cfg%n_phi))
    allocate(q(size(theta), size(cfg%husimi_kicks)))
    call equal_area_grid(cfg%n_theta, cfg%n_phi, theta, phi)
  end if

  psi = coherent_state(cfg%j, cfg%dim, cfg%theta0, cfg%phi0)
  allocate(obs(5, 0:cfg%n_kicks))
  do n = 0, cfg%n_kicks
    if (n > 0) psi = matmul(u, psi)
    obs(1:3, n) = expectation_j(psi, cfg%j)/cfg%j
    obs(4, n) = sum(abs(psi)**2)
    obs(5, n) = (cfg%j*(cfg%j + 1.0_wp) - cfg%j**2*sum(obs(1:3, n)**2))/cfg%j**2
    if (cfg%husimi) then
      do i = 1, size(cfg%husimi_kicks)
        if (cfg%husimi_kicks(i) == n) call husimi(psi, cfg%j, theta, phi, q(:, i))
      end do
    end if
  end do
  call write_expectation(cfg, obs)
  call add('max_norm_error', maxval(abs(obs(4, :) - 1.0_wp)))

  if (cfg%husimi) then
    call write_husimi(cfg, theta, phi, q)
    do i = 1, size(cfg%husimi_kicks)
      write(name, '(a,i0)') 'husimi_norm_n', cfg%husimi_kicks(i)
      call add(name, husimi_norm(q(:, i), cfg%j))
    end do
  end if

  if (cfg%spectrum) then
    call quasienergies(u, v, eps, sector, max_dev)
    allocate(s_centre(cfg%n_bins), density(cfg%n_bins))
    call spacing_statistics(eps, sector, cfg%n_bins, cfg%s_max, s_centre, density, mean_r)
    call write_quasienergies(cfg, eps, sector)
    call write_spacings(cfg, s_centre, density)
    call add('max_eigenvalue_modulus_error', max_dev)
    call add('mean_gap_ratio', mean_r)
    call add('mean_gap_ratio_poisson', 2.0_wp*log(2.0_wp) - 1.0_wp)
    call add('mean_gap_ratio_coe', 0.5307_wp)
  end if

  write(*, '(a,i0,a)') 'kicked_top_quantum: dim = ', cfg%dim, ' -> '//trim(cfg%out_dir)
  call write_summary(cfg, names, values)

contains

  ! Appends one named check to the summary.
  subroutine add(label, value)
    character(len=*), intent(in) :: label
    real(wp), intent(in) :: value

    names = [character(len=40) :: names, label]
    values = [values, value]
  end subroutine add

end program kicked_top_quantum

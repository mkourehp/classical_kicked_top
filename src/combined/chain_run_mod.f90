! Time evolution of the qubit chain for one rotation angle p, shared by
! kicked_top_combined (make combined) and kicked_top_spin_chain (make spin_chain):
! t_r forward steps U, then (with butterfly_effect) sigma^a on qubit i_flip,
! a = flip_axis, then t_r backward steps (coupling 2 pi - k/(2j), kick -p).
! Writes the files of combined_io_mod in dir.
module chain_run_mod
  use kinds_mod, only: wp
  use io_utils_mod, only: make_dir
  use husimi_mod, only: husimi, husimi_norm
  use combined_config_mod, only: rconfig_t
  use combined_io_mod, only: write_chain, write_husimi, write_summary
  use chain_state_mod, only: product_coherent_state, dicke_amplitudes
  use chain_floquet_mod, only: flip_qubit, forward_step, backward_step
  use chain_observables_mod, only: fidelity, site_spins, total_spin
  implicit none
  private
  public :: run_chain

contains

  ! Forward evolution for t_r steps, then (after the optional flip) backward evolution for t_r steps.
  ! obs(:, t_r) and the Husimi snapshot at t_r are taken before the flip.
  subroutine run_chain(cfg, p, theta0, phi0, e_zz, hus_theta, hus_phi, dir, prog)
    type(rconfig_t), intent(in) :: cfg
    real(wp), intent(in) :: p
    real(wp), intent(in) :: theta0, phi0      ! direction of the initial polarisation
    real(wp), intent(in) :: e_zz(0:)          ! diagonal of H_zz
    real(wp), intent(in) :: hus_theta(:), hus_phi(:)
    character(len=*), intent(in) :: dir
    character(len=*), intent(in) :: prog      ! program name for the progress line
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
      if (t >= 1 .and. t <= cfg%t_r) call forward_step(psi, cfg%n_q, p, cfg%k, e_zz)
      if (t == cfg%t_r + 1 .and. cfg%butterfly_effect) call flip_qubit(psi, cfg%i_flip, cfg%flip_axis)
      if (t > cfg%t_r) call backward_step(psi, cfg%n_q, p, cfg%k, e_zz)

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
    call write_chain(cfg, p, dir, obs)
    call write_husimi(cfg, p, dir, steps, hus_theta, hus_phi, q)
    call write_summary(cfg, p, dir, &
      [character(len=24) :: 'F(t_r)', 'F(2t_r)', 'max|norm-1|', 'max|S^2-j(j+1)|', &
       'max(1-symmetric_weight)', 'husimi_norm(0)', 'husimi_norm(t_r)', 'husimi_norm(2t_r)'], &
      [obs(1, cfg%t_r), obs(1, 2*cfg%t_r), maxval(abs(obs(2, :) - 1.0_wp)), &
       maxval(abs(obs(9, :) - j*(j + 1.0_wp))), sym_dev, &
       husimi_norm(q(:, 1), j), husimi_norm(q(:, 2), j), husimi_norm(q(:, 3), j)])

    write(*, '(a,f0.4,a,f0.4,a,es10.3,a)') prog//': p = ', p, ', k = ', cfg%k, &
      ', F(2t_r) - 1 = ', obs(1, 2*cfg%t_r) - 1.0_wp, ' -> '//dir
  end subroutine run_chain

end module chain_run_mod

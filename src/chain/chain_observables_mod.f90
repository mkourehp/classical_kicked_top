! Observables of the qubit chain, computed from the state vector (basis as in
! chain_state_mod): fidelity, single-qubit spins and the total spin S = sum_i sigma_i/2.
module chain_observables_mod
  use kinds_mod, only: wp
  implicit none
  private
  public :: fidelity, site_spins, total_spin

contains

  ! Fidelity F = |<psi0|psi>|**2.
  pure function fidelity(psi0, psi) result(f)
    complex(wp), intent(in) :: psi0(0:)
    complex(wp), intent(in) :: psi(0:)
    real(wp) :: f

    f = abs(dot_product(psi0, psi))**2
  end function fidelity

  ! sig(:, i) = (<sigma^x_i>, <sigma^y_i>, <sigma^z_i>) of qubit i = 1 .. n_q.
  ! For a pair (up, down) differing in bit i: <sigma^x> += 2 Re(conj(up) down),
  ! <sigma^y> += 2 Im(conj(up) down); <sigma^z> = sum |psi_b|**2 (+1 up, -1 down).
  pure function site_spins(psi, n_q) result(sig)
    complex(wp), intent(in) :: psi(0:)
    integer, intent(in) :: n_q
    real(wp) :: sig(3, n_q)
    complex(wp) :: z
    integer :: i, b, mask

    sig = 0.0_wp
    do i = 1, n_q
      mask = ishft(1, i - 1)
      do b = 0, size(psi) - 1
        if (iand(b, mask) /= 0) then
          sig(3, i) = sig(3, i) + abs(psi(b))**2
        else
          sig(3, i) = sig(3, i) - abs(psi(b))**2
          z = conjg(psi(ior(b, mask)))*psi(b)
          sig(1, i) = sig(1, i) + 2.0_wp*real(z, wp)
          sig(2, i) = sig(2, i) + 2.0_wp*aimag(z)
        end if
      end do
    end do
  end function site_spins

  ! Total spin: [<S_x>, <S_y>, <S_z>, <S_x**2>, <S_y**2>, <S_z**2>, <S**2>],
  ! with <S_a**2> = ||S_a psi||**2. S_a psi is built in a work vector of the
  ! size of psi; S**2 = j(j+1) with j = n_q/2 in the symmetric subspace.
  function total_spin(psi, n_q) result(st)
    complex(wp), intent(in) :: psi(0:)
    integer, intent(in) :: n_q
    real(wp) :: st(7)
    complex(wp), allocatable :: w(:)
    real(wp) :: sig(3, n_q)
    integer :: i, b, mask

    sig = site_spins(psi, n_q)
    st(1:3) = 0.5_wp*sum(sig, dim=2)

    allocate(w(0:size(psi) - 1))
    ! S_x psi
    !$omp parallel do private(i) schedule(static)
    do b = 0, size(psi) - 1
      w(b) = (0.0_wp, 0.0_wp)
      do i = 0, n_q - 1
        w(b) = w(b) + 0.5_wp*psi(ieor(b, ishft(1, i)))
      end do
    end do
    !$omp end parallel do
    st(4) = sum(abs(w)**2)
    ! S_y psi: sigma^y |down> = i |up>, sigma^y |up> = -i |down>
    !$omp parallel do private(i, mask) schedule(static)
    do b = 0, size(psi) - 1
      w(b) = (0.0_wp, 0.0_wp)
      do i = 0, n_q - 1
        mask = ishft(1, i)
        if (iand(b, mask) /= 0) then
          w(b) = w(b) + cmplx(0.0_wp, 0.5_wp, wp)*psi(ieor(b, mask))
        else
          w(b) = w(b) - cmplx(0.0_wp, 0.5_wp, wp)*psi(ieor(b, mask))
        end if
      end do
    end do
    !$omp end parallel do
    st(5) = sum(abs(w)**2)
    ! S_z psi: diagonal, S_z = u - n_q/2
    st(6) = 0.0_wp
    do b = 0, size(psi) - 1
      st(6) = st(6) + (popcnt(b) - 0.5_wp*n_q)**2*abs(psi(b))**2
    end do
    st(7) = st(4) + st(5) + st(6)
  end function total_spin

end module chain_observables_mod

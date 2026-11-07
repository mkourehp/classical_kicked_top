! Floquet step of the qubit chain, applied to the state vector only (no
! Hamiltonian or 2**n_q x 2**n_q matrix is built). Basis as in chain_state_mod.
!
! The chain realises the kicked top with J = (1/2) sum_i sigma_i, j = n_q/2:
!   H(t) = (k/(2j)) H_zz + p J_y sum_n delta(t-n),   H_zz = sum_{i,l} V_il sigma^z_i sigma^z_l,
! with i and l running over all qubits (i = l included). The couplings act
! continuously; at integer times a pulse rotates every qubit about y by p.
! For V_il = 1/4, H_zz = (1/2 sum_i sigma^z_i)**2 = J_z**2. Sampled just before
! each kick, one period is
!   U = exp(-i (k/(2j)) H_zz) exp(-i p J_y),
! the kicked-top Floquet operator of floquet_mod, lecture notes eq. (39) [eq:floquet].
module chain_floquet_mod
  use kinds_mod, only: wp
  implicit none
  private
  public :: uniform_coupling, zz_energies, zz_coupling, rotate_y, forward_step, backward_step

contains

  ! V_il = 1/4 for all i, l: H_zz = J_z**2 (kicked-top limit).
  pure function uniform_coupling(n_q) result(v)
    integer, intent(in) :: n_q
    real(wp) :: v(n_q, n_q)

    v = 0.25_wp
  end function uniform_coupling

  ! Diagonal of H_zz: e(b) = sum_{i,l} V_il z_i z_l, z_i = +1 if bit i-1 of b
  ! is 1 (qubit i up) and -1 otherwise. Computed once per run.
  function zz_energies(n_q, v) result(e)
    integer, intent(in) :: n_q
    real(wp), intent(in) :: v(n_q, n_q)
    real(wp), allocatable :: e(:)
    real(wp) :: z(n_q)
    integer :: b, i

    allocate(e(0:2**n_q - 1))
    !$omp parallel do private(z, i) schedule(static)
    do b = 0, 2**n_q - 1
      do i = 1, n_q
        z(i) = merge(1.0_wp, -1.0_wp, btest(b, i - 1))
      end do
      e(b) = dot_product(z, matmul(v, z))
    end do
    !$omp end parallel do
  end function zz_energies

  ! psi -> exp(-i t H_zz) psi, with e the diagonal of H_zz from zz_energies.
  subroutine zz_coupling(psi, e, t)
    complex(wp), intent(inout) :: psi(0:)
    real(wp), intent(in) :: e(0:)
    real(wp), intent(in) :: t
    integer :: b

    !$omp parallel do schedule(static)
    do b = 0, size(psi) - 1
      psi(b) = psi(b)*cmplx(cos(t*e(b)), -sin(t*e(b)), wp)
    end do
    !$omp end parallel do
  end subroutine zz_coupling

  ! psi -> exp(-i alpha J_y) psi = prod_i exp(-i alpha/2 sigma^y_i) psi.
  ! On qubit i, in the basis (up, down): [[cos(alpha/2), -sin(alpha/2)], [sin(alpha/2), cos(alpha/2)]],
  ! applied to every pair of basis states that differ only in bit i.
  ! Rotates the spin about y by alpha, lecture notes eq. (17) [eq:rotation-y].
  subroutine rotate_y(psi, n_q, alpha)
    complex(wp), intent(inout) :: psi(0:)
    integer, intent(in) :: n_q
    real(wp), intent(in) :: alpha
    real(wp) :: c, s
    complex(wp) :: a_up, a_dn
    integer :: i, b, mask

    c = cos(alpha/2.0_wp)
    s = sin(alpha/2.0_wp)
    do i = 0, n_q - 1
      mask = ishft(1, i)
      !$omp parallel do private(a_up, a_dn) schedule(static)
      do b = 0, size(psi) - 1
        if (iand(b, mask) /= 0) cycle
        a_dn = psi(b)
        a_up = psi(ior(b, mask))
        psi(ior(b, mask)) = c*a_up - s*a_dn
        psi(b) = s*a_up + c*a_dn
      end do
      !$omp end parallel do
    end do
  end subroutine rotate_y

  ! One period: U = exp(-i (k/(2j)) H_zz) exp(-i p J_y), j = n_q/2
  ! (first the kick, then the couplings for one period).
  subroutine forward_step(psi, n_q, p, k, e)
    complex(wp), intent(inout) :: psi(0:)
    integer, intent(in) :: n_q
    real(wp), intent(in) :: p
    real(wp), intent(in) :: k
    real(wp), intent(in) :: e(0:)

    call rotate_y(psi, n_q, p)
    call zz_coupling(psi, e, k/n_q)
  end subroutine forward_step

  ! Inverse period: U^dagger = exp(+i p J_y) exp(+i (k/(2j)) H_zz).
  subroutine backward_step(psi, n_q, p, k, e)
    complex(wp), intent(inout) :: psi(0:)
    integer, intent(in) :: n_q
    real(wp), intent(in) :: p
    real(wp), intent(in) :: k
    real(wp), intent(in) :: e(0:)

    call zz_coupling(psi, e, -k/n_q)
    call rotate_y(psi, n_q, -p)
  end subroutine backward_step

end module chain_floquet_mod

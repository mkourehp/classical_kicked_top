! Floquet steps of the qubit chain, applied to the state vector only (no
! Hamiltonian or 2**n_q x 2**n_q matrix is built). Basis as in chain_state_mod.
!
! J_z**2 is realised by static all-to-all couplings,
!   H_zz = sum_{i<l} sigma^z_i sigma^z_l = 2 S_z**2 - n_q/2 = (M**2 - n_q)/2,   M = 2u - n_q,
! with integer eigenvalues, so exp(-i (2 pi + eps) H_zz) = exp(-i eps H_zz):
! the coupling acting for T = 2 pi + eps is the kicked-top torsion with
! k/(2j) = 2 eps, i.e. k = 2 n_q eps, lecture notes eq. (39) [eq:floquet].
! The rotation exp(-i h S_y) is a product of single-qubit rotations; h = p.
!
! Time reversal (as in spin_chain): T -> T' = 2 pi - eps, h -> -h and the
! order of the two factors swapped, so that backward_step = forward_step^dagger.
module chain_floquet_mod
  use kinds_mod, only: wp, two_pi
  implicit none
  private
  public :: zz_coupling, rotate_y, forward_step, backward_step

contains

  ! psi -> exp(-i t H_zz) psi. H_zz is diagonal and depends only on the number
  ! u of up qubits, so one phase per u.
  subroutine zz_coupling(psi, n_q, t)
    complex(wp), intent(inout) :: psi(0:)
    integer, intent(in) :: n_q
    real(wp), intent(in) :: t
    complex(wp) :: phase(0:n_q)
    integer :: u, b

    do u = 0, n_q
      phase(u) = exp(cmplx(0.0_wp, -t*real(((2*u - n_q)**2 - n_q)/2, wp), wp))
    end do
    !$omp parallel do schedule(static)
    do b = 0, size(psi) - 1
      psi(b) = psi(b)*phase(popcnt(b))
    end do
    !$omp end parallel do
  end subroutine zz_coupling

  ! psi -> exp(-i alpha S_y) psi = prod_i exp(-i alpha/2 sigma^y_i) psi.
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

  ! One forward period: U_f = exp(-i (2 pi + eps) H_zz) exp(-i h S_y)
  ! (rotation first, then the couplings, as in the kicked-top map).
  subroutine forward_step(psi, n_q, h, eps)
    complex(wp), intent(inout) :: psi(0:)
    integer, intent(in) :: n_q
    real(wp), intent(in) :: h
    real(wp), intent(in) :: eps

    call rotate_y(psi, n_q, h)
    call zz_coupling(psi, n_q, two_pi + eps)
  end subroutine forward_step

  ! One backward period: U_b = exp(+i h S_y) exp(-i (2 pi - eps) H_zz) = U_f^dagger.
  subroutine backward_step(psi, n_q, h, eps)
    complex(wp), intent(inout) :: psi(0:)
    integer, intent(in) :: n_q
    real(wp), intent(in) :: h
    real(wp), intent(in) :: eps

    call zz_coupling(psi, n_q, two_pi - eps)
    call rotate_y(psi, n_q, -h)
  end subroutine backward_step

end module chain_floquet_mod

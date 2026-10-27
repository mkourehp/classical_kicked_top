! States of the qubit chain: n_q spins 1/2 whose total spin S = sum_i sigma_i/2
! realises the spin J of the kicked top, j = n_q/2 for a ferromagnetic state.
! Computational basis |b>, b = 0 .. 2**n_q - 1: bit i of b is 1 if qubit i is
! up (sigma^z_i = +1) and 0 if it is down, as in the spin_chain repository.
! Only the 2**n_q state vector is stored; no operator matrix is built.
module chain_state_mod
  use kinds_mod, only: wp, pi
  implicit none
  private
  public :: ferromagnetic_direction, product_coherent_state, dicke_amplitudes

contains

  ! Direction (theta0, phi0) of the ferromagnetic initial state given by s0 and
  ! pol, as in spin_chain: s0 = -1 (or 2**n_q - 1) puts all qubits up along the
  ! axis pol = 'z', 'x' or 'y', s0 = 0 puts all of them down. Other s0 are not
  ! ferromagnetic and are rejected.
  subroutine ferromagnetic_direction(s0, n_q, pol, theta0, phi0)
    integer, intent(in) :: s0
    integer, intent(in) :: n_q
    character(len=*), intent(in) :: pol
    real(wp), intent(out) :: theta0
    real(wp), intent(out) :: phi0
    logical :: up

    if (s0 == -1 .or. s0 == 2**n_q - 1) then
      up = .true.
    else if (s0 == 0) then
      up = .false.
    else
      error stop 'chain_state_mod: only ferromagnetic initial states are supported (s0 = -1 or 0)'
    end if

    select case (pol)
    case ('z')
      theta0 = 0.0_wp
      phi0 = 0.0_wp
    case ('x')
      theta0 = pi/2.0_wp
      phi0 = 0.0_wp
    case ('y')
      theta0 = pi/2.0_wp
      phi0 = pi/2.0_wp
    case default
      error stop "chain_state_mod: pol must be 'z', 'x' or 'y'"
    end select

    if (.not. up) then
      theta0 = pi - theta0
      phi0 = phi0 + pi
    end if
  end subroutine ferromagnetic_direction

  ! Product of n_q single-qubit coherent states,
  !   |theta,phi>^(x n_q),  |theta,phi>_1 = e^{-i phi sigma^z/2} e^{-i theta sigma^y/2} |up>
  !                                       = cos(theta/2) e^{-i phi/2} |up> + sin(theta/2) e^{i phi/2} |down>.
  ! The amplitude of |b> depends only on the number u of up qubits in b. Its
  ! projection on the symmetric subspace is the spin coherent state of spin_mod,
  ! lecture notes eq. (41) [eq:coherent].
  pure function product_coherent_state(n_q, theta, phi) result(psi)
    integer, intent(in) :: n_q
    real(wp), intent(in) :: theta
    real(wp), intent(in) :: phi
    complex(wp) :: psi(0:2**n_q - 1)
    complex(wp) :: a_up, a_dn, amp(0:n_q)
    integer :: u, b

    a_up = cos(theta/2.0_wp)*exp(cmplx(0.0_wp, -phi/2.0_wp, wp))
    a_dn = sin(theta/2.0_wp)*exp(cmplx(0.0_wp, phi/2.0_wp, wp))
    do u = 0, n_q
      amp(u) = a_up**u*a_dn**(n_q - u)
    end do
    do b = 0, size(psi) - 1
      psi(b) = amp(popcnt(b))
    end do
  end function product_coherent_state

  ! Amplitudes c_m = <j,m|psi> of psi on the Dicke states (symmetric subspace,
  ! j = n_q/2), |j,m> = binom(n_q,u)**(-1/2) * sum of |b> with u = j + m up qubits.
  ! Ordered as in spin_mod: index i = 1 .. n_q+1 holds m = j - (i-1), so that
  ! spin_mod and husimi_mod can be applied to c. sum |c_m|**2 = 1 if psi is
  ! symmetric under permutations of the qubits.
  pure function dicke_amplitudes(psi, n_q) result(c)
    complex(wp), intent(in) :: psi(0:)
    integer, intent(in) :: n_q
    complex(wp) :: c(n_q + 1)
    complex(wp) :: a(0:n_q)
    real(wp) :: log_binom
    integer :: u, b

    a = (0.0_wp, 0.0_wp)
    do b = 0, size(psi) - 1
      a(popcnt(b)) = a(popcnt(b)) + psi(b)
    end do
    do u = 0, n_q
      log_binom = log_gamma(real(n_q + 1, wp)) - log_gamma(real(u + 1, wp)) - log_gamma(real(n_q - u + 1, wp))
      c(n_q - u + 1) = a(u)*exp(-0.5_wp*log_binom)
    end do
  end function dicke_amplitudes

end module chain_state_mod

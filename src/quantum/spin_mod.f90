! Spin-j operators and spin coherent states (units hbar = 1).
! Basis |j,m>, m = j, j-1, ..., -j, stored at index i = 1, ..., 2j+1 with m = j - (i-1).
module spin_mod
  use kinds_mod, only: wp
  implicit none
  private
  public :: m_values, jy_matrix, expectation_j, coherent_state

contains

  ! The values m_i = j - (i-1), i = 1, ..., n = 2j+1.
  pure function m_values(j, n) result(m)
    real(wp), intent(in) :: j
    integer, intent(in) :: n
    real(wp) :: m(n)
    integer :: i

    m = [(j - real(i - 1, wp), i = 1, n)]
  end function m_values

  ! Ladder coefficient a(m) = sqrt(j(j+1) - m(m+1)) of J_+|j,m> = a(m)|j,m+1>.
  elemental function ladder(j, m) result(a)
    real(wp), intent(in) :: j
    real(wp), intent(in) :: m
    real(wp) :: a

    a = sqrt(max(j*(j + 1.0_wp) - m*(m + 1.0_wp), 0.0_wp))
  end function ladder

  ! J_y = (J_+ - J_-)/(2i) as a dense Hermitian matrix:
  ! <j,m+1|J_y|j,m> = -i a(m)/2, the only non-zero elements besides their conjugates.
  pure function jy_matrix(j, n) result(jy)
    real(wp), intent(in) :: j
    integer, intent(in) :: n
    complex(wp) :: jy(n, n)
    real(wp) :: m(n)
    integer :: i

    m = m_values(j, n)
    jy = (0.0_wp, 0.0_wp)
    do i = 2, n
      ! index i holds m_i, index i-1 holds m_i + 1
      jy(i - 1, i) = cmplx(0.0_wp, -0.5_wp*ladder(j, m(i)), wp)
      jy(i, i - 1) = conjg(jy(i - 1, i))
    end do
  end function jy_matrix

  ! [<J_x>, <J_y>, <J_z>] in the normalised state psi, using
  ! <J_+> = <J_x> + i<J_y> = sum_m conj(psi_{m+1}) a(m) psi_m.
  pure function expectation_j(psi, j) result(ev)
    complex(wp), intent(in) :: psi(:)
    real(wp), intent(in) :: j
    real(wp) :: ev(3)
    real(wp) :: m(size(psi))
    complex(wp) :: jplus
    integer :: i

    m = m_values(j, size(psi))
    jplus = (0.0_wp, 0.0_wp)
    do i = 2, size(psi)
      jplus = jplus + conjg(psi(i - 1))*ladder(j, m(i))*psi(i)
    end do
    ev = [real(jplus, wp), aimag(jplus), sum(m*abs(psi)**2)]
  end function expectation_j

  ! Spin coherent state |theta,phi> = exp(-i phi J_z) exp(-i theta J_y) |j,j>, theta in [0, pi]:
  !   c_m = sqrt(binom(2j, j-m)) cos(theta/2)**(j+m) sin(theta/2)**(j-m) exp(-i m phi).
  ! It is normalised and <J>/j = (sin(theta)cos(phi), sin(theta)sin(phi), cos(theta)).
  ! The amplitudes are computed from logarithms, so large j does not overflow.
  pure function coherent_state(j, n, theta, phi) result(c)
    real(wp), intent(in) :: j
    integer, intent(in) :: n
    real(wp), intent(in) :: theta
    real(wp), intent(in) :: phi
    complex(wp) :: c(n)
    real(wp) :: m(n), log_binom, log_amp
    integer :: i, r

    m = m_values(j, n)
    log_binom = 0.0_wp                       ! log binom(2j, r), r = j - m = i - 1
    do i = 1, n
      r = i - 1
      if (r > 0) log_binom = log_binom + log((2.0_wp*j - real(r - 1, wp))/real(r, wp))
      log_amp = 0.5_wp*log_binom + log_power(cos(0.5_wp*theta), j + m(i)) &
                + log_power(sin(0.5_wp*theta), j - m(i))
      c(i) = exp(log_amp)*cmplx(cos(m(i)*phi), -sin(m(i)*phi), wp)
    end do
  end function coherent_state

  ! e*log(x) for x >= 0, with 0**0 = 1 and a large negative value for 0**e, e > 0.
  pure function log_power(x, e) result(r)
    real(wp), intent(in) :: x
    real(wp), intent(in) :: e
    real(wp) :: r

    if (e < 0.5_wp) then
      r = 0.0_wp
    else if (x <= 0.0_wp) then
      r = -1.0e300_wp
    else
      r = e*log(x)
    end if
  end function log_power

end module spin_mod

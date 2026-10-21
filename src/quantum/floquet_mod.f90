! One-period (Floquet) operator of the quantum kicked top,
!   U = exp(-i k/(2j) J_z**2) exp(-i p J_y)      (first the rotation, then the kick),
! the quantum counterpart of the classical map, Sec. 9 [sec:quantum].
module floquet_mod
  use kinds_mod, only: wp
  use spin_mod, only: m_values, jy_matrix
  implicit none
  private
  public :: diagonalise_jy, exp_jy, floquet_operator

contains

  ! Eigen-decomposition J_y = V diag(mu) V^dagger with LAPACK zheev.
  ! The eigenvalues mu are returned in ascending order (-j, -j+1, ..., j).
  subroutine diagonalise_jy(j, n, v, mu)
    real(wp), intent(in) :: j
    integer, intent(in) :: n
    complex(wp), intent(out) :: v(n, n)
    real(wp), intent(out) :: mu(n)
    complex(wp), allocatable :: work(:)
    complex(wp) :: work_query(1)
    real(wp) :: rwork(max(1, 3*n - 2))
    integer :: info
    external :: zheev

    v = jy_matrix(j, n)
    call zheev('V', 'U', n, v, n, mu, work_query, -1, rwork, info)
    allocate(work(max(1, int(real(work_query(1), wp)))))
    call zheev('V', 'U', n, v, n, mu, work, size(work), rwork, info)
    if (info /= 0) error stop 'floquet_mod: zheev failed'
  end subroutine diagonalise_jy

  ! exp(-i alpha J_y) = V diag(exp(-i alpha mu)) V^dagger, a rotation about y by alpha.
  pure function exp_jy(v, mu, alpha) result(r)
    complex(wp), intent(in) :: v(:,:)
    real(wp), intent(in) :: mu(:)
    real(wp), intent(in) :: alpha
    complex(wp) :: r(size(mu), size(mu))
    complex(wp) :: w(size(mu), size(mu))
    integer :: l

    do l = 1, size(mu)
      w(:, l) = v(:, l)*cmplx(cos(alpha*mu(l)), -sin(alpha*mu(l)), wp)
    end do
    r = matmul(w, conjg(transpose(v)))
  end function exp_jy

  ! U = K R with R = exp(-i p J_y) and the diagonal kick
  ! K = diag(exp(-i k m**2/(2j))). v, mu: eigen-decomposition of J_y.
  pure function floquet_operator(j, p, k, v, mu) result(u)
    real(wp), intent(in) :: j
    real(wp), intent(in) :: p
    real(wp), intent(in) :: k
    complex(wp), intent(in) :: v(:,:)
    real(wp), intent(in) :: mu(:)
    complex(wp) :: u(size(mu), size(mu))
    real(wp) :: m(size(mu)), phase
    integer :: i

    m = m_values(j, size(mu))
    u = exp_jy(v, mu, p)
    do i = 1, size(mu)
      phase = 0.5_wp*k*m(i)**2/j
      u(i, :) = u(i, :)*cmplx(cos(phase), -sin(phase), wp)
    end do
  end function floquet_operator

end module floquet_mod

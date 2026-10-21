! Husimi distribution Q(theta, phi) = |<theta,phi|psi>|**2 of a spin-j state psi,
! with |theta,phi> the spin coherent state of spin_mod. Q >= 0, and
! (2j+1)/(4 pi) * integral of Q over the sphere = 1 for a normalised psi.
module husimi_mod
  use kinds_mod, only: wp
  use spin_mod, only: coherent_state
  implicit none
  private
  public :: husimi, husimi_norm

contains

  ! Q at the points (theta(i), phi(i)), in parallel over the points (OpenMP).
  subroutine husimi(psi, j, theta, phi, q)
    complex(wp), intent(in) :: psi(:)
    real(wp), intent(in) :: j
    real(wp), intent(in) :: theta(:)
    real(wp), intent(in) :: phi(:)
    real(wp), intent(out) :: q(:)
    integer :: i

    !$omp parallel do schedule(static)
    do i = 1, size(theta)
      q(i) = abs(dot_product(coherent_state(j, size(psi), theta(i), phi(i)), psi))**2
    end do
    !$omp end parallel do
  end subroutine husimi

  ! (2j+1)/(4 pi) * integral of Q, evaluated on an equal-area grid of n points
  ! (each of area 4 pi/n): (2j+1) * mean(Q). Should be 1.
  pure function husimi_norm(q, j) result(norm)
    real(wp), intent(in) :: q(:)
    real(wp), intent(in) :: j
    real(wp) :: norm

    norm = (2.0_wp*j + 1.0_wp)*sum(q)/real(size(q), wp)
  end function husimi_norm

end module husimi_mod

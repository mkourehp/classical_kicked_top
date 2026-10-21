! Geometry of the unit sphere, shared by the classical and the quantum code.
! Lecture notes: eq. (13) [eq:unit] and eq. (9) [eq:area].
module sphere_mod
  use kinds_mod, only: wp, pi, two_pi
  implicit none
  private
  public :: spherical_to_cartesian, equal_area_grid

contains

  ! The unit vector s = (sin(theta)cos(phi), sin(theta)sin(phi), cos(theta)), eq. (13) [eq:unit].
  pure function spherical_to_cartesian(theta, phi) result(s)
    real(wp), intent(in) :: theta
    real(wp), intent(in) :: phi
    real(wp) :: s(3)

    s = [sin(theta)*cos(phi), sin(theta)*sin(phi), cos(theta)]
  end function spherical_to_cartesian

  ! Cell centres (theta(i), phi(i)) of an n_theta x n_phi grid that is uniform
  ! in the canonical coordinates (phi, Z=cos(theta)) of Sec. 1.3 [sec:canonical].
  ! By eq. (9) [eq:area] every cell has the same area 4*pi/(n_theta*n_phi).
  ! Order: Z in the outer loop (from -1 to 1), phi in the inner loop (from -pi to pi).
  pure subroutine equal_area_grid(n_theta, n_phi, theta, phi)
    integer, intent(in) :: n_theta
    integer, intent(in) :: n_phi
    real(wp), intent(out) :: theta(n_theta*n_phi)
    real(wp), intent(out) :: phi(n_theta*n_phi)
    integer :: i, l, n
    real(wp) :: z

    n = 0
    do i = 1, n_theta
      z = -1.0_wp + (real(i, wp) - 0.5_wp)*2.0_wp/real(n_theta, wp)
      do l = 1, n_phi
        n = n + 1
        theta(n) = acos(z)
        phi(n) = -pi + (real(l, wp) - 0.5_wp)*two_pi/real(n_phi, wp)
      end do
    end do
  end subroutine equal_area_grid

end module sphere_mod

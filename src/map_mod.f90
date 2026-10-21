! The stroboscopic map of the classical kicked top and its tangent map.
! The state is the unit vector s = J/j = (X, Y, Z), eq. (13) [eq:unit].
! Lecture notes: Sec. 3 [sec:map] (derivation of the map) and
! Sec. 8.1 [sec:benettin] (tangent map).
module map_mod
  use kinds_mod, only: wp
  implicit none
  private
  public :: kick, kick_tangent

contains

  ! One period of the kicked top, s_{n+1} = F(s_n), eq. (20) [eq:map]:
  ! free precession about y by the angle p, eq. (17) [eq:rotation-y],
  ! followed by the torsion about z by the angle k*Z', eq. (19) [eq:torsion].
  pure function kick(s, p, k) result(s_new)
    real(wp), intent(in) :: s(3)
    real(wp), intent(in) :: p
    real(wp), intent(in) :: k
    real(wp) :: s_new(3)
    real(wp) :: s1(3)

    s1 = rotate_y(s, p)
    s_new = rotate_z(s1, k*s1(3))
  end function kick

  ! Tangent map: the image DF(s)*ds of a small displacement ds of s under one
  ! kick, eq. (34) [eq:tangent]. The rotation is linear; the torsion
  ! R_z(k*Z')s' also turns with the change k*dZ' of its angle, which adds
  ! k*dZ' * (z_hat x s_{n+1}). If ds is tangent to the sphere at s, the result
  ! is tangent at F(s).
  pure function kick_tangent(s, ds, p, k) result(ds_new)
    real(wp), intent(in) :: s(3)
    real(wp), intent(in) :: ds(3)
    real(wp), intent(in) :: p
    real(wp), intent(in) :: k
    real(wp) :: ds_new(3)
    real(wp) :: s1(3), ds1(3), s2(3)

    s1 = rotate_y(s, p)
    ds1 = rotate_y(ds, p)
    s2 = rotate_z(s1, k*s1(3))
    ds_new = rotate_z(ds1, k*s1(3)) + k*ds1(3)*[-s2(2), s2(1), 0.0_wp]
  end function kick_tangent

  ! R_y(alpha)*v, rotation about the y axis by the angle alpha, eq. (11) [eq:rotations].
  pure function rotate_y(v, alpha) result(w)
    real(wp), intent(in) :: v(3)
    real(wp), intent(in) :: alpha
    real(wp) :: w(3)
    real(wp) :: c, sn

    c = cos(alpha)
    sn = sin(alpha)
    w = [v(1)*c + v(3)*sn, v(2), -v(1)*sn + v(3)*c]
  end function rotate_y

  ! R_z(alpha)*v, rotation about the z axis by the angle alpha, eq. (11) [eq:rotations].
  pure function rotate_z(v, alpha) result(w)
    real(wp), intent(in) :: v(3)
    real(wp), intent(in) :: alpha
    real(wp) :: w(3)
    real(wp) :: c, sn

    c = cos(alpha)
    sn = sin(alpha)
    w = [v(1)*c - v(2)*sn, v(1)*sn + v(2)*c, v(3)]
  end function rotate_z

end module map_mod

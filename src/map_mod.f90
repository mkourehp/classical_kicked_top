module map_mod
  use kinds_mod, only: wp
  implicit none
  private
  public :: kick

contains

  ! One period of the kicked top on the unit sphere (Haake convention):
  ! rotation about y by p, then torsion about z by the angle k*Z.
  pure function kick(s, p, k) result(s_new)
    real(wp), intent(in) :: s(3)
    real(wp), intent(in) :: p
    real(wp), intent(in) :: k
    real(wp) :: s_new(3)
    real(wp) :: x, y, z, c, sn

    c = cos(p)
    sn = sin(p)
    x = s(1)*c + s(3)*sn
    y = s(2)
    z = -s(1)*sn + s(3)*c

    c = cos(k*z)
    sn = sin(k*z)
    s_new(1) = x*c - y*sn
    s_new(2) = x*sn + y*c
    s_new(3) = z
  end function kick

end module map_mod

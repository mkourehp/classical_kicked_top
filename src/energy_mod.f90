module energy_mod
  use kinds_mod, only: wp
  implicit none
  private
  public :: energy

contains

  ! Energy terms of H = p*Y + (k/2)*Z**2 at a point s on the unit sphere:
  ! [rotation term p*Y, kick term (k/2)*Z**2, total].
  pure function energy(s, p, k) result(e)
    real(wp), intent(in) :: s(3)
    real(wp), intent(in) :: p
    real(wp), intent(in) :: k
    real(wp) :: e(3)

    e(1) = p*s(2)
    e(2) = 0.5_wp*k*s(3)**2
    e(3) = e(1) + e(2)
  end function energy

end module energy_mod

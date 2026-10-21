! Energy terms of the scaled kicked-top Hamiltonian.
! Lecture notes: Sec. 4 [sec:energy], eq. (23) [eq:energy].
module energy_mod
  use kinds_mod, only: wp
  implicit none
  private
  public :: energy

contains

  ! Energy terms at a point s = (X, Y, Z) of the unit sphere, i.e. the scaled
  ! Hamiltonian h of eq. (14) [eq:h] without the delta comb:
  !   e(1) = E_rot  = p*Y          (free precession about y)
  !   e(2) = E_kick = (k/2)*Z**2   (torsion)
  !   e(3) = E      = e(1) + e(2)
  ! E is not conserved by the map (Sec. 4); it is nearly conserved only when p
  ! and k are both small, where it is the effective Hamiltonian of eq. (24) [eq:heff].
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

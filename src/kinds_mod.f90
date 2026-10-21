! Working precision and mathematical constants shared by all modules.
module kinds_mod
  use, intrinsic :: iso_fortran_env, only: real64
  implicit none
  private
  public :: wp, pi, two_pi

  integer, parameter :: wp = real64          ! double precision for all reals
  real(wp), parameter :: pi = acos(-1.0_wp)
  real(wp), parameter :: two_pi = 2.0_wp*pi
end module kinds_mod

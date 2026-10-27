! Quasienergy spectrum of the Floquet operator U and its level-spacing statistics.
! U commutes with the parity P = exp(-i pi J_y), the quantum version of the
! classical symmetry S = R_y(pi), eq. (25) [eq:symmetry]. Spacing statistics are
! only meaningful within one symmetry sector, so U is diagonalised separately in
! the two parity sectors.
module spectrum_mod
  use kinds_mod, only: wp, pi, two_pi
  implicit none
  private
  public :: quasienergies, spacing_statistics, theory_spacing, commutator_error, unitarity_error

contains

  ! Quasienergies eps in [-pi, pi), U|n> = exp(-i eps_n)|n>, sorted within each sector.
  ! v holds the eigenvectors of J_y for the ascending eigenvalues mu = -j, ..., j.
  ! P = exp(-i pi mu) changes sign from one mu to the next, so the odd and the even
  ! columns of v span the two parity sectors; sector(i) = 1 or 2 labels eps(i).
  ! max_dev = max | |lambda| - 1 | over the eigenvalues lambda of U (0 if U is unitary).
  subroutine quasienergies(u, v, eps, sector, max_dev)
    complex(wp), intent(in) :: u(:,:)
    complex(wp), intent(in) :: v(:,:)
    real(wp), allocatable, intent(out) :: eps(:)
    integer, allocatable, intent(out) :: sector(:)
    real(wp), intent(out) :: max_dev
    complex(wp), allocatable :: vs(:,:), us(:,:), lam(:), work(:)
    real(wp), allocatable :: rwork(:), e(:)
    complex(wp) :: dummy(1, 1), work_query(1)
    integer :: n, ns, s, info, i0
    external :: zgeev

    n = size(u, 1)
    allocate(eps(n), sector(n))
    max_dev = 0.0_wp
    i0 = 0
    do s = 1, 2
      ns = (n - s)/2 + 1                     ! columns s, s+2, ..., <= n
      allocate(vs(n, ns), us(ns, ns), lam(ns), rwork(2*ns), e(ns))
      vs = v(:, s:n:2)
      us = matmul(conjg(transpose(vs)), matmul(u, vs))
      call zgeev('N', 'N', ns, us, ns, lam, dummy, 1, dummy, 1, work_query, -1, rwork, info)
      allocate(work(max(1, int(real(work_query(1), wp)))))
      call zgeev('N', 'N', ns, us, ns, lam, dummy, 1, dummy, 1, work, size(work), rwork, info)
      if (info /= 0) error stop 'spectrum_mod: zgeev failed'
      max_dev = max(max_dev, maxval(abs(abs(lam) - 1.0_wp)))
      e = -atan2(aimag(lam), real(lam, wp))
      call sort(e)
      eps(i0 + 1:i0 + ns) = e
      sector(i0 + 1:i0 + ns) = s
      i0 = i0 + ns
      deallocate(vs, us, lam, rwork, work, e)
    end do
  end subroutine quasienergies

  ! Unfolded nearest-neighbour spacings and their gap ratios, per sector.
  ! A sector of ns levels on the circle has mean spacing 2 pi/ns, so
  ! s_i = (eps_{i+1} - eps_i) ns/(2 pi) has mean 1 (the last spacing wraps around the circle).
  ! density(b): histogram of the s_i on n_bins bins of [0, s_max] with centres s_centre(b),
  ! normalised as a probability density (sum(density)*ds = fraction of s_i below s_max).
  ! mean_r: average of r_i = min(d_i, d_{i+1})/max(d_i, d_{i+1}) over neighbouring spacings,
  ! which needs no unfolding.
  subroutine spacing_statistics(eps, sector, n_bins, s_max, s_centre, density, mean_r)
    real(wp), intent(in) :: eps(:)
    integer, intent(in) :: sector(:)
    integer, intent(in) :: n_bins
    real(wp), intent(in) :: s_max
    real(wp), intent(out) :: s_centre(n_bins)
    real(wp), intent(out) :: density(n_bins)
    real(wp), intent(out) :: mean_r
    real(wp), allocatable :: e(:), d(:)
    real(wp) :: ds, r_sum, d_next
    integer :: s, ns, i, b, n_total

    ds = s_max/real(n_bins, wp)
    s_centre = [((real(b, wp) - 0.5_wp)*ds, b = 1, n_bins)]
    density = 0.0_wp
    r_sum = 0.0_wp
    n_total = 0
    do s = 1, 2
      ns = count(sector == s)
      if (ns < 3) cycle
      allocate(e(ns), d(ns))
      e = pack(eps, sector == s)
      d(1:ns - 1) = e(2:ns) - e(1:ns - 1)
      d(ns) = e(1) + two_pi - e(ns)
      do i = 1, ns
        b = int(d(i)*real(ns, wp)/(two_pi*ds)) + 1
        if (b <= n_bins) density(b) = density(b) + 1.0_wp
        d_next = d(mod(i, ns) + 1)
        if (max(d(i), d_next) > 0.0_wp) r_sum = r_sum + min(d(i), d_next)/max(d(i), d_next)
      end do
      n_total = n_total + ns
      deallocate(e, d)
    end do
    if (n_total == 0) error stop 'spectrum_mod: need at least 3 levels per sector'
    density = density/(real(n_total, wp)*ds)
    mean_r = r_sum/real(n_total, wp)
  end subroutine spacing_statistics

  ! Reference spacing distributions at s: [Poisson exp(-s) (regular dynamics),
  ! COE Wigner surmise (pi/2) s exp(-pi s**2/4) (chaotic, time-reversal invariant),
  ! CUE (32/pi**2) s**2 exp(-4 s**2/pi) (chaotic, no time-reversal symmetry)].
  elemental function theory_spacing(s, ensemble) result(ps)
    real(wp), intent(in) :: s
    integer, intent(in) :: ensemble      ! 1 Poisson, 2 COE, 3 CUE
    real(wp) :: ps

    select case (ensemble)
    case (1)
      ps = exp(-s)
    case (2)
      ps = 0.5_wp*pi*s*exp(-0.25_wp*pi*s**2)
    case default
      ps = 32.0_wp/pi**2*s**2*exp(-4.0_wp*s**2/pi)
    end select
  end function theory_spacing

  ! max |(AB - BA)_{ab}|; zero if A and B commute.
  pure function commutator_error(a, b) result(err)
    complex(wp), intent(in) :: a(:,:)
    complex(wp), intent(in) :: b(:,:)
    real(wp) :: err

    err = maxval(abs(matmul(a, b) - matmul(b, a)))
  end function commutator_error

  ! max |(U^dagger U - 1)_{ab}|; zero if U is unitary.
  pure function unitarity_error(u) result(err)
    complex(wp), intent(in) :: u(:,:)
    real(wp) :: err
    complex(wp) :: w(size(u, 1), size(u, 2))
    integer :: i

    w = matmul(conjg(transpose(u)), u)
    do i = 1, size(u, 1)
      w(i, i) = w(i, i) - 1.0_wp
    end do
    err = maxval(abs(w))
  end function unitarity_error

  ! Ascending insertion sort.
  pure subroutine sort(x)
    real(wp), intent(inout) :: x(:)
    real(wp) :: t
    integer :: i, l

    do i = 2, size(x)
      t = x(i)
      l = i - 1
      do while (l >= 1)
        if (x(l) <= t) exit
        x(l + 1) = x(l)
        l = l - 1
      end do
      x(l + 1) = t
    end do
  end subroutine sort

end module spectrum_mod

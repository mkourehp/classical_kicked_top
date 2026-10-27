! Plot-ready ASCII output of the quantum kicked top. Python only reads these files.
! Every file starts with a header line '# p = ..., k = ..., j = ...' and a column line.
module quantum_io_mod
  use kinds_mod, only: wp
  use quantum_config_mod, only: qconfig_t
  use sphere_mod, only: spherical_to_cartesian
  use spectrum_mod, only: theory_spacing
  implicit none
  private
  public :: write_expectation, write_husimi, write_quasienergies, write_spacings, write_summary

contains

  ! expectation.dat, one line per kick n:
  !   n  <J_x>/j  <J_y>/j  <J_z>/j  norm  spread
  ! obs(:, n) = [<J>/j (3), <psi|psi>, (<J^2> - |<J>|^2)/j^2].
  subroutine write_expectation(cfg, obs)
    type(qconfig_t), intent(in) :: cfg
    real(wp), intent(in) :: obs(:, 0:)
    integer :: u, n

    call open_file(cfg, 'expectation.dat', 'n Jx/j Jy/j Jz/j norm spread', u)
    do n = 0, ubound(obs, 2)
      write(u, '(i0,5(1x,es23.15e3))') n, obs(:, n)
    end do
    close(u)
  end subroutine write_expectation

  ! husimi.dat, one block per snapshot (kick n), blocks separated by a blank line:
  !   n  phi  Z  X  Y  Q
  ! on the equal-area grid (theta, phi); q(:, i) is the snapshot at cfg%husimi_kicks(i).
  subroutine write_husimi(cfg, theta, phi, q)
    type(qconfig_t), intent(in) :: cfg
    real(wp), intent(in) :: theta(:)
    real(wp), intent(in) :: phi(:)
    real(wp), intent(in) :: q(:,:)
    real(wp) :: s(3)
    integer :: u, i, l

    call open_file(cfg, 'husimi.dat', 'n phi Z X Y Q', u)
    do i = 1, size(q, 2)
      do l = 1, size(theta)
        s = spherical_to_cartesian(theta(l), phi(l))
        write(u, '(i0,5(1x,es23.15e3))') cfg%husimi_kicks(i), phi(l), s(3), s(1), s(2), q(l, i)
      end do
      write(u, '(a)') ''
    end do
    close(u)
  end subroutine write_husimi

  ! quasienergies.dat, one line per eigenvalue exp(-i eps) of U:  sector  eps
  subroutine write_quasienergies(cfg, eps, sector)
    type(qconfig_t), intent(in) :: cfg
    real(wp), intent(in) :: eps(:)
    integer, intent(in) :: sector(:)
    integer :: u, i

    call open_file(cfg, 'quasienergies.dat', 'sector eps', u)
    do i = 1, size(eps)
      write(u, '(i0,1x,es23.15e3)') sector(i), eps(i)
    end do
    close(u)
  end subroutine write_quasienergies

  ! spacings.dat, one line per histogram bin:  s  P(s)  Poisson  COE  CUE
  subroutine write_spacings(cfg, s_centre, density)
    type(qconfig_t), intent(in) :: cfg
    real(wp), intent(in) :: s_centre(:)
    real(wp), intent(in) :: density(:)
    integer :: u, b

    call open_file(cfg, 'spacings.dat', 's P(s) Poisson COE CUE', u)
    do b = 1, size(s_centre)
      write(u, '(5(es23.15e3,1x))') s_centre(b), density(b), theory_spacing(s_centre(b), [1, 2, 3])
    end do
    close(u)
  end subroutine write_spacings

  ! summary.dat: one 'name value' line per check, also printed to standard output.
  subroutine write_summary(cfg, names, values)
    type(qconfig_t), intent(in) :: cfg
    character(len=*), intent(in) :: names(:)
    real(wp), intent(in) :: values(:)
    integer :: u, i

    call open_file(cfg, 'summary.dat', 'name value', u)
    do i = 1, size(names)
      write(u, '(a,1x,es23.15e3)') trim(names(i)), values(i)
      write(*, '(2x,a,t36,es12.4)') trim(names(i)), values(i)
    end do
    close(u)
  end subroutine write_summary

  ! Opens cfg%out_dir/name for writing and writes the two header lines.
  subroutine open_file(cfg, name, columns, u)
    type(qconfig_t), intent(in) :: cfg
    character(len=*), intent(in) :: name
    character(len=*), intent(in) :: columns
    integer, intent(out) :: u
    character(len=96) :: buf

    open(newunit=u, file=trim(cfg%out_dir)//'/'//name, status='replace', action='write')
    write(buf, '(a,f0.6,a,f0.6,a,f0.1)') 'p = ', cfg%p, ', k = ', cfg%k, ', j = ', cfg%j
    write(u, '(a)') '# '//trim(buf)
    write(u, '(a)') '# '//columns
  end subroutine open_file

end module quantum_io_mod

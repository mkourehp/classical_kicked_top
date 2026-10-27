! ASCII output of kicked_top_regimes. Layout of cfg%out_dir:
!   regimes.dat                       one line per run: i, h, k, eps, n_q
!   h_<i>/classical/section/          classical Poincare section on a grid (files of io_mod)
!   h_<i>/classical/orbit/            classical orbit from the initial state (files of io_mod)
!   h_<i>/chain/fidelity.dat          t  F  norm
!   h_<i>/chain/total_spin.dat        t  <Sx> <Sy> <Sz>  <Sx^2> <Sy^2> <Sz^2>  <S^2>
!   h_<i>/chain/site_spins.dat        t  <sx_1..n_q>  <sy_1..n_q>  <sz_1..n_q>
!   h_<i>/chain/husimi.dat            t  phi Z X Y Q   (blank line between snapshots)
!   h_<i>/chain/summary.dat           name value (correctness checks)
! Steps t = 0 .. t_r are forward, t_r+1 .. 2 t_r backward.
module regimes_io_mod
  use kinds_mod, only: wp
  use regimes_config_mod, only: rconfig_t
  use sphere_mod, only: spherical_to_cartesian
  implicit none
  private
  public :: run_dir, write_regimes, write_chain, write_husimi, write_summary

contains

  ! Output directory of run i.
  function run_dir(cfg, i) result(dir)
    type(rconfig_t), intent(in) :: cfg
    integer, intent(in) :: i
    character(len=:), allocatable :: dir
    character(len=16) :: buf

    write(buf, '(i0)') i
    dir = trim(cfg%out_dir)//'/h_'//trim(buf)
  end function run_dir

  subroutine write_regimes(cfg)
    type(rconfig_t), intent(in) :: cfg
    integer :: u, i

    open(newunit=u, file=trim(cfg%out_dir)//'/regimes.dat', status='replace', action='write')
    write(u, '(a)') '# kicked_top_regimes: k = 2 n_q eps, T = 2 pi + eps, T'' = 2 pi - eps'
    write(u, '(a)') '# i h k eps n_q'
    do i = 1, size(cfg%h)
      write(u, '(i0,3(1x,es23.15e3),1x,i0)') i, cfg%h(i), cfg%k, cfg%eps, cfg%n_q
    end do
    close(u)
  end subroutine write_regimes

  subroutine open_file(dir, name, title, columns, u)
    character(len=*), intent(in) :: dir, name, title, columns
    integer, intent(out) :: u

    open(newunit=u, file=dir//'/'//name, status='replace', action='write')
    write(u, '(a)') '# '//title
    write(u, '(a)') '# '//columns
  end subroutine open_file

  ! Header line of the chain files.
  function title(cfg, h) result(t)
    type(rconfig_t), intent(in) :: cfg
    real(wp), intent(in) :: h
    character(len=128) :: t

    write(t, '(a,i0,a,f0.4,a,f0.4,a,f0.4,a,i0)') 'n_q = ', cfg%n_q, ', h = ', h, ', eps = ', cfg%eps, &
      ', k = ', cfg%k, ', t_r = ', cfg%t_r
  end function title

  ! obs(:, t): F, norm, <S_a>, <S_a^2>, <S^2>, site spins (3*n_q), t = 0 .. 2 t_r.
  subroutine write_chain(cfg, h, dir, obs)
    type(rconfig_t), intent(in) :: cfg
    real(wp), intent(in) :: h
    character(len=*), intent(in) :: dir
    real(wp), intent(in) :: obs(:, 0:)
    integer :: uf, us, uq, t
    character(len=32) :: fmt

    call open_file(dir, 'fidelity.dat', trim(title(cfg, h)), 't F norm', uf)
    call open_file(dir, 'total_spin.dat', trim(title(cfg, h)), &
      't <Sx> <Sy> <Sz> <Sx^2> <Sy^2> <Sz^2> <S^2>', us)
    call open_file(dir, 'site_spins.dat', trim(title(cfg, h)), &
      't <sigma_x>(1..n_q) <sigma_y>(1..n_q) <sigma_z>(1..n_q)', uq)
    write(fmt, '(a,i0,a)') '(i0,', 3*cfg%n_q, '(1x,es23.15e3))'
    do t = 0, ubound(obs, 2)
      write(uf, '(i0,2(1x,es23.15e3))') t, obs(1:2, t)
      write(us, '(i0,7(1x,es23.15e3))') t, obs(3:9, t)
      write(uq, fmt) t, obs(10:, t)
    end do
    close(uf)
    close(us)
    close(uq)
  end subroutine write_chain

  ! q(:, i) is the Husimi distribution at step steps(i) on the points (theta, phi).
  subroutine write_husimi(cfg, h, dir, steps, theta, phi, q)
    type(rconfig_t), intent(in) :: cfg
    real(wp), intent(in) :: h
    character(len=*), intent(in) :: dir
    integer, intent(in) :: steps(:)
    real(wp), intent(in) :: theta(:), phi(:), q(:,:)
    real(wp) :: s(3)
    integer :: u, i, l

    call open_file(dir, 'husimi.dat', trim(title(cfg, h)), 't phi Z X Y Q', u)
    do i = 1, size(steps)
      do l = 1, size(theta)
        s = spherical_to_cartesian(theta(l), phi(l))
        write(u, '(i0,5(1x,es23.15e3))') steps(i), phi(l), s(3), s(1), s(2), q(l, i)
      end do
      write(u, '(a)') ''
    end do
    close(u)
  end subroutine write_husimi

  subroutine write_summary(cfg, h, dir, names, values)
    type(rconfig_t), intent(in) :: cfg
    real(wp), intent(in) :: h
    character(len=*), intent(in) :: dir
    character(len=*), intent(in) :: names(:)
    real(wp), intent(in) :: values(:)
    integer :: u, i

    call open_file(dir, 'summary.dat', trim(title(cfg, h)), 'name value', u)
    do i = 1, size(names)
      write(u, '(a,1x,es23.15e3)') trim(names(i)), values(i)
    end do
    close(u)
  end subroutine write_summary

end module regimes_io_mod

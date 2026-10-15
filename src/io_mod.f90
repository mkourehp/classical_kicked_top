module io_mod
  use kinds_mod, only: wp
  use config_mod, only: config_t
  use energy_mod, only: energy
  implicit none
  private
  public :: write_results

  character(len=*), parameter :: trajectory_file = 'trajectory.dat'
  character(len=*), parameter :: phase_space_file = 'phase_space.dat'
  character(len=*), parameter :: energy_file = 'energy.dat'

contains

  ! Writes plot-ready ASCII files to cfg%out_dir. Orbits are separated by a blank line.
  !   trajectory.dat  : ic  n  X  Y  Z
  !   phase_space.dat : ic  n  phi  Z      (phi in (-pi, pi])
  !   energy.dat      : pY  k/2*Z**2  total (same rows as the files above)
  subroutine write_results(cfg, orbit)
    type(config_t), intent(in) :: cfg
    real(wp), intent(in) :: orbit(:,:,:)
    character(len=:), allocatable :: title
    integer :: u_traj, u_ps, u_en, ic, isave, n
    character(len=64) :: buf

    call execute_command_line('mkdir -p "'//trim(cfg%out_dir)//'"')

    write(buf, '(a,f0.6,a,f0.6)') 'p = ', cfg%p, ', k = ', cfg%k
    title = trim(buf)

    open(newunit=u_traj, file=trim(cfg%out_dir)//'/'//trajectory_file, status='replace', action='write')
    open(newunit=u_ps, file=trim(cfg%out_dir)//'/'//phase_space_file, status='replace', action='write')
    open(newunit=u_en, file=trim(cfg%out_dir)//'/'//energy_file, status='replace', action='write')

    write(u_traj, '(a)') '# '//title
    write(u_traj, '(a)') '# ic n X Y Z'
    write(u_ps, '(a)') '# '//title
    write(u_ps, '(a)') '# ic n phi Z'
    write(u_en, '(a)') '# '//title
    write(u_en, '(a)') '# pY k/2*Z^2 total'

    do ic = 1, size(orbit, 3)
      do isave = 1, size(orbit, 2)
        n = (isave - 1)*cfg%save_every
        write(u_traj, '(i0,1x,i0,3(1x,es23.15e3))') ic, n, orbit(:, isave, ic)
        write(u_ps, '(i0,1x,i0,2(1x,es23.15e3))') ic, n, &
          atan2(orbit(2, isave, ic), orbit(1, isave, ic)), orbit(3, isave, ic)
        write(u_en, '(3(es23.15e3,1x))') energy(orbit(:, isave, ic), cfg%p, cfg%k)
      end do
      write(u_traj, '(a)') ''
      write(u_ps, '(a)') ''
      write(u_en, '(a)') ''
    end do

    close(u_traj)
    close(u_ps)
    close(u_en)
  end subroutine write_results

end module io_mod

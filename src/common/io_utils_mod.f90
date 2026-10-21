! Helpers for reading namelist files and writing output, shared by the
! classical and the quantum code.
module io_utils_mod
  implicit none
  private
  public :: check_nml_read, make_dir

contains

  ! Stops on a malformed namelist group. A missing group (end of file) is
  ! accepted, so its variables keep their defaults. caller names the module.
  subroutine check_nml_read(ios, group, caller)
    integer, intent(in) :: ios
    character(len=*), intent(in) :: group
    character(len=*), intent(in) :: caller

    if (ios /= 0 .and. .not. is_iostat_end(ios)) then
      error stop caller//': error reading namelist group &'//group
    end if
  end subroutine check_nml_read

  ! Creates the directory path (with parents) if it does not exist.
  subroutine make_dir(path)
    character(len=*), intent(in) :: path

    call execute_command_line('mkdir -p "'//trim(path)//'"')
  end subroutine make_dir

end module io_utils_mod

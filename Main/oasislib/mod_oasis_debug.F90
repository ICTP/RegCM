!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#ifdef OASIS

module mod_oasis_debug

  use mod_intkinds
  use mod_message
  use mod_stdio
  use mod_dynparam
  use mod_mppparam

  implicit none

  private

  integer(ik4) , public :: proc_unit
#ifdef DEBUG
  integer(ik4) :: tab_level
#endif

  public :: oasisxregcm_open_log , oasisxregcm_close_log
#ifdef DEBUG
  public :: checkpoint_main , checkpoint_proc
  !public :: checkpoint_time
  public :: checkpoint_enter , checkpoint_exit
#endif

  contains

    ! create oasisxregcm information log files
    subroutine oasisxregcm_open_log(component_name,component_id)
      implicit none
      character(len=6) , intent(in) :: component_name
      integer(ik4) , intent(in) :: component_id
      character(len=128) :: comp_log
      character(len=3) :: suffix
      integer(ik4) :: iretval
#ifdef DEBUG
      character(len=*) , parameter :: sub_name = 'oasisxregcm_open_log'
#endif
      !--------------------------------------------------------------------------
      write(suffix,'(i3.3)') myid
      comp_log = 'oaxre_'//suffix
      open(newunit=proc_unit,file=trim(comp_log),status='replace', &
              action='write',form='formatted',iostat=iretval)
      if ( iretval /= 0 ) then
        call fatal(__FILE__,__LINE__,'Cannot open OASISxREGCM log files!')
      end if
      write(proc_unit,*) '--------------------------------------------------'
      write(proc_unit,*) 'This file is an output file generated if the model'
      write(proc_unit,*) ' was compiled with -DOASIS flag (--enable-oasis'
      write(proc_unit,*) ' with configure)'
      write(proc_unit,*) ''
      write(proc_unit,*) '     Component: ', trim(component_name)
      write(proc_unit,"(A,I2)") '        with id: ', component_id 
      write(proc_unit,*) '          Rank: ', myid
      write(proc_unit,*) '--------------------------------------------------'
      call flush(proc_unit)
#ifdef DEBUG
      tab_level = 0
      call checkpoint_enter(sub_name)
#endif
      write(proc_unit,*) "All processes output files setup"
      call flush(proc_unit)
#ifdef DEBUG
      call checkpoint_exit(sub_name)
#endif
    end subroutine oasisxregcm_open_log

    ! close the opened log files
    subroutine oasisxregcm_close_log
      implicit none
#ifdef DEBUG
      character(len=*) , parameter :: sub_name = 'oasisxregcm_close_log'
#endif
      !--------------------------------------------------------------------------
#ifdef DEBUG
      call checkpoint_enter(sub_name)
#endif
      write(proc_unit,*) '--------------------------------------------------'
      write(proc_unit,*) 'This file is being closed'
      write(proc_unit,*) '--------------------------------------------------'
      call flush(proc_unit)
      close(proc_unit)
    end subroutine oasisxregcm_close_log

#ifdef DEBUG
    ! write a debug checkpoint statement in stdout
    ! (mixed with other components' outputs)
    subroutine checkpoint_main(sub,comment)
      implicit none
      ! the subroutine name where this subroutine is called in
      character(len=*) , intent(in) :: sub
      ! comment to add to the debug statement
      character(len=*) , intent(in) :: comment
      !--------------------------------------------------------------------------
      if ( myid == italk ) write(stdout,*) '(',sub,') ',comment
    end subroutine checkpoint_main

    ! write a debug checkpoint statement in each proc log file
    subroutine checkpoint_proc(sub,comment)
      implicit none
      ! the subroutine name where this subroutine is called in
      character(len=*) , intent(in) :: sub
      ! comment to add to the debug statement
      character(len=*) , intent(in) :: comment
      !--------------------------------------------------------------------------
      write(proc_unit,*) '(',sub,') ',comment 
      call flush(proc_unit)
    end subroutine checkpoint_proc

    ! write the a time indication statement
    subroutine checkpoint_time(sub,time)
      implicit none
      ! the subroutine name where this subroutine is called in
      character(len=*) , intent(in) :: sub
      ! time to be printed
      integer(ik4) , intent(in) :: time
      character(len=128) :: comment
      !--------------------------------------------------------------------------
      write(comment,*) 'execution time: ', time
      call checkpoint_proc(sub,trim(comment))
    end subroutine checkpoint_time

    ! write an 'ENTER' statement and change of indent level appropriately
    subroutine checkpoint_enter(sub)
      implicit none
      character(len=*) , intent(in) :: sub
      character(len=:) , allocatable :: indent
      !--------------------------------------------------------------------------
      call tab_string(indent)
      write(proc_unit,*) indent,' ENTER (',sub,')'
      call flush(proc_unit)
      tab_level=tab_level+1
    end subroutine checkpoint_enter

    ! write an 'EXIT' statement and change of indent level appropriately
    subroutine checkpoint_exit(sub)
      implicit none
      character(len=*) , intent(in) :: sub
      character(len=:) , allocatable :: indent
      !--------------------------------------------------------------------------
      tab_level=tab_level-1
      call tab_string(indent)
      write(proc_unit,*) indent,' EXIT  (',sub,')'
      call flush(proc_unit)
    end subroutine checkpoint_exit

    ! return the indentation string corresponding on tab_level
    subroutine tab_string(string)
      implicit none
      character(len=:) , allocatable , intent(out) :: string
      !--------------------------------------------------------------------------
      select case ( tab_level )
      case (0)
        string = '-'
      case (1)
        string = '---'
      case (2)
        string = '-----'
      case (3)
        string = '-------'
      case (4)
        string = '---------'
      case (5)
        string = '-----------'
      case (6)
        string = '-------------'
      case default
        string = '-----sup-----'
      end select
    end subroutine tab_string
#endif

end module mod_oasis_debug
!
#endif
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

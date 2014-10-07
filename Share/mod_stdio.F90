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

module mod_stdio
!
  use mod_intkinds
  use iso_fortran_env

  implicit none

  private

  integer(ik4) , public , parameter :: stdin = input_unit
  integer(ik4) , public , parameter :: stdout = output_unit
  integer(ik4) , public , parameter :: stderr = error_unit

  integer(ik4) , parameter :: file_maxunit = 99
  integer(ik4) , parameter :: file_minunit = 10

  logical , dimension(file_minunit:file_maxunit) :: unit_tag = .false.

  public :: file_getunit , file_freeunit

  contains

    integer(ik4) function file_getunit(unitn)
      implicit none
      integer(ik4) , optional , intent(in) :: unitn
      integer(ik4) :: n
      logical :: isopened
      ! The "I want that unit" case
      file_getunit = -1
      if (present (unitn)) then
        ! Is a valid unit number ? Has been already given out ?
        if (unitn < file_minunit .or. unitn > file_maxunit .or. &
            unit_tag(unitn)) then
          return
        end if
        ! Is it opened by some other part of the program ?
        inquire(unitn, opened=isopened)
        if (isopened) then
          file_getunit = -1
          return
        else
          ! Nulla osta, grant it.
          unit_tag (unitn) = .true.
          file_getunit = unitn
          return
        end if
      else
        do n = file_maxunit, file_minunit, -1
          inquire(n,opened=isopened)
          if ( isopened ) then
            ! First loop get all opened files and flags then
            if ( .not. unit_tag(n) ) unit_tag(n) = .true.
            cycle
          end if
          ! This one is not opened and is not granted.
          if (.not. unit_tag(n)) then
            ! Grant it
            unit_tag(n) = .true.
            file_getunit = n
            return
          end if
        end do
      end if
      ! No avail...
      return
    end function file_getunit

    subroutine file_freeunit(unitn)
      implicit none
      integer(ik4) , intent(in) :: unitn
      if ( unit_tag(unitn) ) then
        close(unitn)
        unit_tag(unitn) = .false.
      end if
    end subroutine file_freeunit

end module mod_stdio
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

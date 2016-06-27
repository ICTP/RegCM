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

module mod_ch_fnest
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_memutil
  use mod_grid
  use mod_wrtoxd
  use mod_interp
  use mod_date
  use mod_nchelper
  use netcdf

  implicit none

  private

  integer(ik4) :: ncid

  data ncid /-1/

  public :: header_fnest , get_fnest , close_fnest

  contains

  subroutine header_fnest(idate,cdir,cname)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    character(len=*) , intent(in) :: cdir , cname
  end subroutine header_fnest

  subroutine get_fnest(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
  end subroutine get_fnest

  subroutine close_fnest
    use netcdf
    implicit none
    integer(ik4) :: istatus
    if ( ncid > 0 ) then
      istatus = nf90_close(ncid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error close fnest file')
    end if
  end subroutine close_fnest

end module mod_ch_fnest
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

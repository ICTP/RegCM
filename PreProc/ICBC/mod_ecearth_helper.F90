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

module mod_ecearth_helper

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_date

  private

  public :: echvars
  public :: find_ecearth_sst
  public :: find_ecearth_dim , find_ecearth_file

  integer(ik4) , parameter :: nvars = 6
  character(len=3) , target , dimension(nvars) :: echvars = &
            ['t  ' , 'z  ' , 'q  ' , 'u  ' , 'v  ', 'XXX']

  contains

  subroutine find_ecearth_sst(fname,idate)
    implicit none
    character(len=256) , intent(out) :: fname
    type(rcm_time_and_date) , intent(in) :: idate
    if ( .not. date_in_scenario(idate,5,.true.) ) then
      fname = trim(inpglob)//'/EC-EARTH/SST/RF/ich1_sst_1950-2009.nc'
    else
      fname = trim(inpglob)//'/EC-EARTH/SST/RCP'//ssttyp(4:5)//&
                '/ic'//ssttyp(4:4)//'1_sst_2006-2100.nc'
    end if
  end subroutine find_ecearth_sst

  subroutine find_ecearth_dim(dim_filename)
    implicit none
    character(len=256) , intent(out) :: dim_filename
    ! Just return the name of one file in the historical dataset
    ! we hope is there.
    dim_filename = trim(inpglob)//'/EC-EARTH/fixed/ecearth.nc'
  end subroutine find_ecearth_dim

  subroutine find_ecearth_file(ecearth_filename,var,idate)
    implicit none
    character(len=256) , intent(out) :: ecearth_filename
    character(len=*) , intent(in) :: var
    type(rcm_time_and_date) , intent(in) :: idate
    character(len=256) :: inname
    integer(ik4) :: y , m , d , h
    call split_idate(idate,y,m,d,h)
    if ( m == 1 .and. d == 1 .and. h == 0 ) y = y - 1
    if ( .not. date_in_scenario(idate,5,.true.) ) then
      write (inname,'(a,a,i0.4,a,a,a,i0.4,a)') 'RF', pthsep, y, pthsep, &
                 'ich1_', trim(var)//'_', y, '.nc'
    else
      write (inname,'(a,a,i0.4,a,a,a,i0.4,a)') 'RCP'//dattyp(4:5), pthsep, y, &
                  pthsep, 'ich1_', trim(var)//'_', y, '.nc'
    end if
    ecearth_filename = trim(inpglob)//'/EC-EARTH/'//inname
  end subroutine find_ecearth_file

end module mod_ecearth_helper
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

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

module mod_canesm_helper

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_date

  private

  public :: cavars
  public :: find_canesm_sst
  public :: find_canesm_dim , find_canesm_topo , find_canesm_file

  integer(ik4) , parameter :: nvars = 6
  character(len=3) , target , dimension(nvars) :: cavars = &
            (/'ta ' , 'XXX' , 'hus' , 'ua ' , 'va ' , 'ps '/)

  character(len=64) :: cabase1  = '_6hrLev_CanESM2_historical'
  character(len=64) :: cabase2  = '_6hrLev_CanESM2_rcp'
  character(len=64) :: cabase3 = '_r1i1p1_'

  contains

  subroutine find_canesm_sst(fname,idate)
    implicit none
    character(len=256) , intent(out) :: fname
    type(rcm_time_and_date) , intent(in) :: idate
    if ( .not. date_in_scenario(idate,5,.true.) ) then
      fname = trim(inpglob)//pthsep//'CanESM2'//pthsep//'SST'// &
              pthsep//'ts_Amon_CanESM2_historical'// &
              '_r1i1p1_185001-200512.nc'
    else
      fname = trim(inpglob)//pthsep//'CanESM2'//pthsep//'SST'// &
              pthsep//'ts_Amon_CanESM2_rcp'//ssttyp(4:5)//  &
              '_r1i1p1_200601-210012.nc'
    end if
  end subroutine find_canesm_sst

  subroutine assemble_path(fname,scen,var,d1,d2)
    implicit none
    character(len=256) , intent(out) :: fname
    character(len=*) , intent(in) :: scen
    character(len=*) , intent(in) :: var
    character(len=*) , intent(in) :: d1
    character(len=*) , intent(in) :: d2
    if ( scen == 'RF' ) then
      fname = trim(inpglob)//pthsep//'CanESM2'//pthsep//trim(scen)// &
              pthsep//trim(var)//pthsep//trim(var)//trim(cabase1)//  &
              trim(cabase3)//trim(d1)//'00-'//trim(d2)//'00.nc'
    else
      fname = trim(inpglob)//pthsep//'CanESM2'//pthsep//trim(scen)// &
              pthsep//trim(var)//pthsep//trim(var)//trim(cabase2)//  &
              scen(4:5)//trim(cabase3)//trim(d1)//'00-'//trim(d2)//'00.nc'
    end if
  end subroutine assemble_path

  subroutine find_canesm_dim(dim_filename)
    implicit none
    character(len=256) , intent(out) :: dim_filename
    ! Just return the name of one file in the historical dataset
    ! we hope is there.
    call assemble_path(dim_filename,'RF','ta','1950010100','1950123118')
  end subroutine find_canesm_dim

  subroutine find_canesm_topo(topo_filename)
    implicit none
    character(len=256) , intent(out) :: topo_filename
    topo_filename = trim(inpglob)//pthsep//'CanESM2'//pthsep//'fixed'// &
              pthsep//'orog_fx_CanESM2_historical_r0i0p0.nc'
  end subroutine find_canesm_topo

  subroutine find_canesm_file(canesm_filename,var,idate)
    implicit none
    character(len=256) , intent(out) :: canesm_filename
    character(len=*) , intent(in) :: var
    type(rcm_time_and_date) , intent(in) :: idate
    character(len=10) :: d1 , d2
    integer(ik4) :: y , m , d , h
    call split_idate(idate,y,m,d,h)
    write(d1,'(i0.4,i0.2,i0.2,i0.2)') y, 1, 1, 0
    write(d2,'(i0.4,i0.2,i0.2,i0.2)') y, 12, 31, 18
    if ( .not. date_in_scenario(idate,5,.true.) ) then
      call assemble_path(canesm_filename,'RF',var,d1,d2)
    else
      call assemble_path(canesm_filename,'RCP'//dattyp(4:5),var,d1,d2)
    end if
  end subroutine find_canesm_file

end module mod_canesm_helper
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

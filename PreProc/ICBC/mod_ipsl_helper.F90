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

module mod_ipsl_helper

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_date

  private

  public :: ipvars
  public :: find_ipsl_sst
  public :: find_ipsl_dim , find_ipsl_topo , find_ipsl_file

  integer(ik4) , parameter :: nvars = 6
  character(len=3) , target , dimension(nvars) :: ipvars = &
            (/'ta ' , 'XXX' , 'hus' , 'ua ' , 'va ' , 'ps '/)

  character(len=64) :: ipbase1  = '_6hrLev_IPSL-CM5A-LR_historical'
  character(len=64) :: ipbase2  = '_6hrLev_IPSL-CM5A-LR_rcp'
  character(len=64) :: ipbase3 = '_r1i1p1_'

  contains

  subroutine find_ipsl_sst(fname,idate)
    implicit none
    character(len=256) , intent(out) :: fname
    type(rcm_time_and_date) , intent(in) :: idate
    if ( .not. date_in_scenario(idate,5,.true.) ) then
      fname = trim(inpglob)//pthsep//'IPSL-CM5A-LR'//pthsep//'SST'// &
              pthsep//'tos_Omon_IPSL-CM5A-LR_historical'// &
              '_r1i1p1_185001-200512.nc'
    else
      fname = trim(inpglob)//pthsep//'IPSL-CM5A-LR'//pthsep//'SST'// &
              pthsep//'tos_Omon_IPSL-CM5A-LR_rcp'//ssttyp(4:5)//  &
              '_r1i1p1_200601-230012.nc'
    end if
  end subroutine find_ipsl_sst

  subroutine assemble_path(fname,scen,var,d1,d2)
    implicit none
    character(len=256) , intent(out) :: fname
    character(len=*) , intent(in) :: scen
    character(len=*) , intent(in) :: var
    character(len=*) , intent(in) :: d1
    character(len=*) , intent(in) :: d2
    if ( scen == 'RF' ) then
      fname = trim(inpglob)//pthsep//'IPSL-CM5A-LR'//pthsep//trim(scen)// &
              pthsep//trim(var)//pthsep//trim(var)//trim(ipbase1)//  &
              trim(ipbase3)//trim(d1)//'00-'//trim(d2)//'00.nc'
    else
      fname = trim(inpglob)//pthsep//'IPSL-CM5A-LR'//pthsep//trim(scen)// &
              pthsep//trim(var)//pthsep//trim(var)//trim(ipbase2)//  &
              scen(4:5)//trim(ipbase3)//trim(d1)//'00-'//trim(d2)//'00.nc'
    end if
  end subroutine assemble_path

  subroutine find_ipsl_dim(dim_filename)
    implicit none
    character(len=256) , intent(out) :: dim_filename
    ! Just return the name of one file in the historical dataset
    ! we hope is there.
    call assemble_path(dim_filename,'RF','ta','1950010103','1959123121')
  end subroutine find_ipsl_dim

  subroutine find_ipsl_topo(topo_filename)
    implicit none
    character(len=256) , intent(out) :: topo_filename
    topo_filename = trim(inpglob)//pthsep//'IPSL-CM5A-LR'//pthsep//'fixed'// &
              pthsep//'orog_fx_IPSL-CM5A-LR_historical_r0i0p0.nc'
  end subroutine find_ipsl_topo

  subroutine find_ipsl_file(ipsl_filename,var,idate)
    implicit none
    character(len=256) , intent(out) :: ipsl_filename
    character(len=*) , intent(in) :: var
    type(rcm_time_and_date) , intent(in) :: idate
    character(len=10) :: d1 , d2
    integer(ik4) :: y , m , d , h
    integer(ik4) :: y1 , y2
    call split_idate(idate,y,m,d,h)
    select case (var)
      case ('ps')
        if ( y < 1950 ) then
          y1 = 1900
          y2 = 1949
        else if ( y < 2000 ) then
          y1 = 1950
          y2 = 1999
        else if ( y > 2000 .and. y < 2006 ) then
          y1 = 2000
          y2 = 2005
        else
          y1 = y/50*50+6
          y2 = y1 + 49
        end if
      case default
        if ( y > 2000 .and. y < 2006 ) then
          y1 = 2000
          y2 = 2005
        else if ( y < 2000 ) then
          y1 = y/10*10
          y2 = y1+9
        else
          y1 = y/10*10+6
          y2 = y1 + 9
        end if
    end select
    write(d1,'(i0.4i0.2i0.2i0.2)') y1, 1, 1, 3
    write(d2,'(i0.4i0.2i0.2i0.2)') y2, 12, 31, 21
    if ( .not. date_in_scenario(idate,5,.true.) ) then
      call assemble_path(ipsl_filename,'RF',var,d1,d2)
    else
      call assemble_path(ipsl_filename,'RCP'//dattyp(4:5),var,d1,d2)
    end if
  end subroutine find_ipsl_file

end module mod_ipsl_helper
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

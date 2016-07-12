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

module mod_csiro_helper

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_date

  private

  public :: csirvars
  public :: find_csiro_sst
  public :: find_csiro_dim , find_csiro_topo , find_csiro_file

  integer(ik4) , parameter :: nvars = 6
  character(len=3) , target , dimension(nvars) :: csirvars = &
            (/'ta ' , 'XXX' , 'hus' , 'ua ' , 'va ' , 'ps '/)

  character(len=64) :: csirbase1  = '_6hrLev_CSIRO-Mk3-6-0_historical'
  character(len=64) :: csirbase2  = '_6hrLev_CSIRO-Mk3-6-0_rcp'
  character(len=64) :: csirbase3 = '_r1i1p1_'

  contains

  subroutine find_csiro_sst(fname,idate)
    implicit none
    character(len=256) , intent(out) :: fname
    type(rcm_time_and_date) , intent(in) :: idate
    if ( .not. date_in_scenario(idate,5) ) then
      fname = trim(inpglob)//pthsep//'CSIRO-MK36'//pthsep//'SST'// &
              pthsep//'tos_Omon_CSIRO-Mk3-6-0_historical'// &
              '_r1i1p1_185001-200512.nc'
    else
      fname = trim(inpglob)//pthsep//'CSIRO-MK36'//pthsep//'SST'// &
              pthsep//'tos_Omon_CSIRO-Mk3-6-0_rcp'//ssttyp(4:5)//  &
              '_r1i1p1_200601-210012.nc'
    end if
  end subroutine find_csiro_sst

  subroutine assemble_path(fname,scen,var,d1,d2)
    implicit none
    character(len=256) , intent(out) :: fname
    character(len=*) , intent(in) :: scen
    character(len=*) , intent(in) :: var
    character(len=*) , intent(in) :: d1
    character(len=*) , intent(in) :: d2
    if ( scen == 'RF' ) then
      fname = trim(inpglob)//pthsep//'CSIRO-MK36'//pthsep//trim(scen)// &
              pthsep//trim(var)//pthsep//trim(var)//trim(csirbase1)//  &
              trim(csirbase3)//trim(d1)//'00-'//trim(d2)//'00.nc'
    else
      fname = trim(inpglob)//pthsep//'CSIRO-MK36'//pthsep//trim(scen)// &
              pthsep//trim(var)//pthsep//trim(var)//trim(csirbase2)//  &
              scen(4:5)//trim(csirbase3)//trim(d1)//'00-'//trim(d2)//'00.nc'
    end if
  end subroutine assemble_path

  subroutine find_csiro_dim(dim_filename)
    implicit none
    character(len=256) , intent(out) :: dim_filename
    ! Just return the name of one file in the historical dataset
    ! we hope is there.
    call assemble_path(dim_filename,'RF','ta','1950010106','1951010100')
  end subroutine find_csiro_dim

  subroutine find_csiro_topo(topo_filename)
    implicit none
    character(len=256) , intent(out) :: topo_filename
    topo_filename = trim(inpglob)//pthsep//'CSIRO-MK36'//pthsep//'fixed'// &
              pthsep//'orog_fx_CSIRO-Mk3-6-0_historical_r0i0p0.nc'
  end subroutine find_csiro_topo

  subroutine find_csiro_file(csiro_filename,var,idate)
    implicit none
    character(len=256) , intent(out) :: csiro_filename
    character(len=*) , intent(in) :: var
    type(rcm_time_and_date) , intent(in) :: idate
    character(len=10) :: d1 , d2
    integer(ik4) :: y , m , d , h
    integer(ik4) :: iyear1 , iyear2
    call split_idate(idate,y,m,d,h)
    select case (var)
      case ('ps')
        if ( .not. date_in_scenario(idate,5,.true.) ) then
          if ( y == 2005 ) then
            if ( m == 1 .and. d == 1 .and. h == 0 ) then
              iyear1 = y/5*5
              if ( mod(y,5) == 0 .and. m == 1 .and. &
                   d == 1 .and. h == 0 ) then
                iyear1 = iyear1 - 5
              end if
              iyear2 = iyear1 + 5
            else
              iyear1 = 2005
              iyear2 = 2006
            end if
          else
            iyear1 = y/5*5
            if ( mod(y,5) == 0 .and. m == 1 .and. &
                 d == 1 .and. h == 0 ) then
              iyear1 = iyear1 - 5
            end if
            iyear2 = iyear1 + 5
          end if
        else
          iyear1 = (y-2006)/5*5+2006
          iyear2 = iyear1 + 5
        end if
        write(d1,'(i0.4i0.2i0.2i0.2)') iyear1, 1, 1, 6
        write(d2,'(i0.4i0.2i0.2i0.2)') iyear2, 1, 1, 0
        if ( .not. date_in_scenario(idate,5,.true.) ) then
          call assemble_path(csiro_filename,'RF',csirvars(6),d1,d2)
        else
          call assemble_path(csiro_filename,'RCP'//dattyp(4:5), &
                             csirvars(6),d1,d2)
        end if
      case default
        iyear1 = y
        if ( m == 1 .and. d == 1 .and. h == 0 ) then
          iyear1 = y-1
        end if
        iyear2 = iyear1+1
        write(d1,'(i0.4i0.2i0.2i0.2)') iyear1, 1, 1, 6
        write(d2,'(i0.4i0.2i0.2i0.2)') iyear2, 1, 1, 0
        if ( .not. date_in_scenario(idate,5,.true.) ) then
          call assemble_path(csiro_filename,'RF',var,d1,d2)
        else
          call assemble_path(csiro_filename,'RCP'//dattyp(4:5),var,d1,d2)
        end if
    end select
  end subroutine find_csiro_file

end module mod_csiro_helper
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

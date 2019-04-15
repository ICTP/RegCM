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

module mod_ccsm4_helper

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_date

  private

  public :: ccsm4vars
  public :: find_ccsm4_sst
  public :: find_ccsm4_dim , find_ccsm4_topo , find_ccsm4_file

  integer(ik4) , parameter :: nvars = 6
  character(len=3) , target , dimension(nvars) :: ccsm4vars = &
            ['ta ' , 'XXX' , 'hus' , 'ua ' , 'va ' , 'ps ']

  character(len=64) :: cabase1  = '_6hrLev_CCSM4_historical'
  character(len=64) :: cabase2  = '_6hrLev_CCSM4_rcp'
  character(len=64) :: cabase3 = '_r6i1p1_'

  contains

  subroutine find_ccsm4_sst(fname,idate)
    implicit none
    character(len=256) , intent(out) :: fname
    type(rcm_time_and_date) , intent(in) :: idate
    if ( .not. date_in_scenario(idate,5,.true.) ) then
      fname = trim(inpglob)//pthsep//'CCSM4'//pthsep//'SST'// &
              pthsep//'tos_Omon_CCSM4_historical'// &
              '_r6i1p1_185001-200512.nc'
    else
      fname = trim(inpglob)//pthsep//'CCSM4'//pthsep//'SST'// &
              pthsep//'tos_Omon_CCSM4_rcp'//ssttyp(4:5)//  &
              '_r6i1p1_200601-210012.nc'
    end if
  end subroutine find_ccsm4_sst

  subroutine assemble_path(fname,scen,var,d1,d2)
    implicit none
    character(len=256) , intent(out) :: fname
    character(len=*) , intent(in) :: scen
    character(len=*) , intent(in) :: var
    character(len=*) , intent(in) :: d1
    character(len=*) , intent(in) :: d2
    if ( scen == 'RF' ) then
      fname = trim(inpglob)//pthsep//'CCSM4'//pthsep//trim(scen)// &
              pthsep//trim(var)//pthsep//trim(var)//trim(cabase1)//  &
              trim(cabase3)//trim(d1)//'-'//trim(d2)//'.nc'
    else
      fname = trim(inpglob)//pthsep//'CCSM4'//pthsep//trim(scen)// &
              pthsep//trim(var)//pthsep//trim(var)//trim(cabase2)//  &
              scen(4:5)//trim(cabase3)//trim(d1)//'-'//trim(d2)//'.nc'
    end if
  end subroutine assemble_path

  subroutine find_ccsm4_dim(dim_filename)
    implicit none
    character(len=256) , intent(out) :: dim_filename
    ! Just return the name of one file in the historical dataset
    ! we hope is there.
    call assemble_path(dim_filename,'RF','ta','1950010106','1950033118')
  end subroutine find_ccsm4_dim

  subroutine find_ccsm4_topo(topo_filename)
    implicit none
    character(len=256) , intent(out) :: topo_filename
    topo_filename = trim(inpglob)//pthsep//'CCSM4'//pthsep//'fixed'// &
              pthsep//'orog_fx_CCSM4_historical_r0i0p0.nc'
  end subroutine find_ccsm4_topo

  subroutine find_ccsm4_file(ccsm4_filename,var,idate)
    implicit none
    character(len=256) , intent(out) :: ccsm4_filename
    character(len=*) , intent(in) :: var
    type(rcm_time_and_date) , intent(in) :: idate
    character(len=10) :: d1 , d2
    integer(ik4) :: y , m , d , h
    call split_idate(idate,y,m,d,h)
    if ( var == 'ps' ) then
      if ( y == 1950 ) then
        write(d1,'(i0.4,i0.2,i0.2,i0.2)') y, 1, 1, 6
      else
        write(d1,'(i0.4,i0.2,i0.2,i0.2)') y, 1, 1, 0
      end if
      write(d2,'(i0.4,i0.2,i0.2,i0.2)') y, 12, 31, 18
    else
      if ( y == 1950 .and. m < 4 ) then
        write(d1,'(i0.4,i0.2,i0.2,i0.2)') y, 1, 1, 6
      else
        m = (m-1)/3*3+1
        write(d1,'(i0.4,i0.2,i0.2,i0.2)') y, (m-1)/3*3+1, 1, 0
      end if
      m = m + 2
      write(d2,'(i0.4,i0.2,i0.2,i0.2)') y, m, ndaypm(y,m,noleap), 18
    end if
    if ( .not. date_in_scenario(idate,5,.true.) ) then
      call assemble_path(ccsm4_filename,'RF',var,d1,d2)
    else
      call assemble_path(ccsm4_filename,'RCP'//dattyp(4:5),var,d1,d2)
    end if
  end subroutine find_ccsm4_file

end module mod_ccsm4_helper
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

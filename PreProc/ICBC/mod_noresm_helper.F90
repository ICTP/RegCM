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

module mod_noresm_helper

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_date

  private

  public :: novars
  public :: find_noresm_sst
  public :: find_noresm_dim , find_noresm_topo , find_noresm_file

  integer(ik4) , parameter :: nvars = 6
  character(len=3) , target , dimension(nvars) :: novars = &
            ['ta ' , 'XXX' , 'hus' , 'ua ' , 'va ' , 'ps ']

  character(len=64) :: nobase1  = '_6hrLev_NorESM1-M_historical'
  character(len=64) :: nobase2  = '_6hrLev_NorESM1-M_rcp'
  character(len=64) :: nobase3 = '_r1i1p1_'

  contains

  subroutine find_noresm_sst(fname,idate)
    implicit none
    character(len=256) , intent(out) :: fname
    type(rcm_time_and_date) , intent(in) :: idate
    character(len=8) :: d1 , d2
    integer(ik4) :: y , m , d , h
    call split_idate(idate,y,m,d,h)
    if ( y == 2006 ) then
      fname = trim(inpglob)//pthsep//'NorESM1-M'//pthsep//'SST'// &
              pthsep//'tos_day_NorESM1-M_historicalExt'// &
              '_r1i1p1_20060101-20061231.nc'
    else
      write(d1,'(i0.4,i0.2,i0.2)') y, 1, 1
      write(d2,'(i0.4,i0.2,i0.2)') y, 12, 31
      if ( .not. date_in_scenario(idate,5,.true.) ) then
        fname = trim(inpglob)//pthsep//'NorESM1-M'//pthsep//'SST'// &
                pthsep//'tos_day_NorESM1-M_historical'// &
                '_r1i1p1_'//d1//'-'//d2//'.nc'
      else
        fname = trim(inpglob)//pthsep//'NorESM1-M'//pthsep//'SST'// &
                pthsep//'tos_day_NorESM1-M_rcp'//ssttyp(4:5)//  &
                '_r1i1p1_'//d1//'-'//d2//'.nc'
      end if
    end if
  end subroutine find_noresm_sst

  subroutine assemble_path(fname,scen,var,d1,d2)
    implicit none
    character(len=256) , intent(out) :: fname
    character(len=*) , intent(in) :: scen
    character(len=*) , intent(in) :: var
    character(len=*) , intent(in) :: d1
    character(len=*) , intent(in) :: d2
    if ( scen == 'RF' ) then
      fname = trim(inpglob)//pthsep//'NorESM1-M'//pthsep//trim(scen)// &
              pthsep//trim(var)//pthsep//trim(var)//trim(nobase1)//    &
              trim(nobase3)//trim(d1)//'-'//trim(d2)//'.nc'
    else
      fname = trim(inpglob)//pthsep//'NorESM1-M'//pthsep//trim(scen)// &
              pthsep//trim(var)//pthsep//trim(var)//trim(nobase2)//    &
              scen(4:5)//trim(nobase3)//trim(d1)//'-'//trim(d2)//'.nc'
    end if
  end subroutine assemble_path

  subroutine find_noresm_dim(dim_filename)
    implicit none
    character(len=256) , intent(out) :: dim_filename
    ! Just return the name of one file in the historical dataset
    ! we hope is there.
    call assemble_path(dim_filename,'RF','ta','1950010100','1950063018')
  end subroutine find_noresm_dim

  subroutine find_noresm_topo(topo_filename)
    implicit none
    character(len=256) , intent(out) :: topo_filename
    topo_filename = trim(inpglob)//pthsep//'NorESM1-M'//pthsep//'fixed'// &
              pthsep//'orog_fx_NorESM1-M_historical_r0i0p0.nc'
  end subroutine find_noresm_topo

  subroutine find_noresm_file(noresm_filename,var,idate)
    implicit none
    character(len=256) , intent(out) :: noresm_filename
    character(len=*) , intent(in) :: var
    type(rcm_time_and_date) , intent(in) :: idate
    character(len=10) :: d1 , d2
    integer(ik4) :: y , m , d , h , y1
    call split_idate(idate,y,m,d,h)
    if ( var == 'ps' ) then
      if ( y < 2000 ) then
        y = y/10*10
        write(d1,'(i0.4,i0.2,i0.2,i0.2)') y, 1, 1, 0
        write(d2,'(i0.4,i0.2,i0.2,i0.2)') y+9, 12, 31, 18
      else if ( y >= 2000 .and. y < 2006 ) then
        write(d1,'(i0.4,i0.2,i0.2,i0.2)') 2000, 1, 1, 0
        write(d2,'(i0.4,i0.2,i0.2,i0.2)') 2005, 12, 31, 18
      else if ( y >= 2006 .and. y < 2100 ) then
        y1 = y/10*10+6
        if ( y < y1 ) y1 = y1 - 10
        write(d1,'(i0.4,i0.2,i0.2,i0.2)') y1, 1, 1, 0
        if ( y >= 2096 ) then
          write(d2,'(i0.4,i0.2,i0.2,i0.2)') 2100, 12, 31, 18
        else
          write(d2,'(i0.4,i0.2,i0.2,i0.2)') y1+9, 12, 31, 18
        end if
      else
        y = y/10*10+1
        write(d1,'(i0.4,i0.2,i0.2,i0.2)') y, 1, 1, 0
        write(d2,'(i0.4,i0.2,i0.2,i0.2)') y+9, 12, 31, 18
      end if
    else
      if ( m <= 6 ) then
        write(d1,'(i0.4,i0.2,i0.2,i0.2)') y, 1, 1, 0
        write(d2,'(i0.4,i0.2,i0.2,i0.2)') y, 6, 30, 18
      else
        write(d1,'(i0.4,i0.2,i0.2,i0.2)') y, 7, 1, 0
        write(d2,'(i0.4,i0.2,i0.2,i0.2)') y, 12, 31, 18
      end if
    end if
    if ( .not. date_in_scenario(idate,5,.true.) ) then
      call assemble_path(noresm_filename,'RF',var,d1,d2)
    else
      call assemble_path(noresm_filename,'RCP'//dattyp(4:5),var,d1,d2)
    end if
  end subroutine find_noresm_file

end module mod_noresm_helper
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

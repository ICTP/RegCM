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

module mod_cnrm_helper

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_date

  private

  public :: cnrmvars
  public :: find_cnrm_sst
  public :: find_cnrm_dim , find_cnrm_topo , find_cnrm_file

  integer(ik4) , parameter :: nvars = 6
  character(len=3) , target , dimension(nvars) :: cnrmvars = &
            (/'ta ' , 'XXX' , 'hus' , 'ua ' , 'va ' , 'ps '/)

  character(len=64) :: cnrmbase1  = '_6hrLev_CNRM-CM5_historical'
  character(len=64) :: cnrmbase2  = '_6hrLev_CNRM-CM5_rcp'
  character(len=64) :: cnrmbase3 = '_r1i1p1_'

  contains

  subroutine find_cnrm_sst(fname,idate)
    implicit none
    character(len=256) , intent(out) :: fname
    type(rcm_time_and_date) , intent(in) :: idate
    integer(ik4) :: y1 , y2
    character(len=6) :: d1 , d2
    integer(ik4) :: y , m , d , h
    call split_idate(idate,y,m,d,h)
    if ( y < 2006 ) then
      y1 = (y)/10*10
      if ( y1 == 2000 ) then
        y2 = 2005
      else
        y2 = y1 + 9
      end if
    else
      y1 = (y-6)/10*10+6
      y2 = y1 + 9
    end if
    write(d1,'(i0.4,i0.2)') y1, 1
    write(d2,'(i0.4,i0.2)') y2, 12
    if ( y1 < 2005 ) then
      fname = trim(inpglob)//pthsep//'CNRM-CM5'//pthsep//'SST'// &
              pthsep//'tos_Omon_CNRM-CM5_historical'// &
              '_r1i1p1_'//d1//'-'//d2//'.nc'
    else
      fname = trim(inpglob)//pthsep//'CNRM-CM5'//pthsep//'SST'// &
              pthsep//'tos_Omon_CNRM-CM5_rcp'//ssttyp(4:5)//  &
              '_r1i1p1_'//d1//'-'//d2//'.nc'
    end if
  end subroutine find_cnrm_sst

  subroutine assemble_path(fname,scen,var,d1,d2)
    implicit none
    character(len=256) , intent(out) :: fname
    character(len=*) , intent(in) :: scen
    character(len=*) , intent(in) :: var
    character(len=*) , intent(in) :: d1
    character(len=*) , intent(in) :: d2
    if ( scen == 'RF' ) then
      fname = trim(inpglob)//pthsep//'CNRM-CM5'//pthsep//trim(scen)// &
              pthsep//trim(var)//pthsep//trim(var)//trim(cnrmbase1)//  &
              trim(cnrmbase3)//trim(d1)//'00-'//trim(d2)//'00.nc'
    else
      fname = trim(inpglob)//pthsep//'CNRM-CM5'//pthsep//trim(scen)// &
              pthsep//trim(var)//pthsep//trim(var)//trim(cnrmbase2)//  &
              scen(4:5)//trim(cnrmbase3)//trim(d1)//'00-'//trim(d2)//'00.nc'
    end if
  end subroutine assemble_path

  subroutine find_cnrm_dim(dim_filename)
    implicit none
    character(len=256) , intent(out) :: dim_filename
    ! Just return the name of one file in the historical dataset
    ! we hope is there.
    call assemble_path(dim_filename,'RF','ta','1950010106','1950020100')
  end subroutine find_cnrm_dim

  subroutine find_cnrm_topo(topo_filename)
    implicit none
    character(len=256) , intent(out) :: topo_filename
    topo_filename = trim(inpglob)//pthsep//'CNRM-CM5'//pthsep//'fixed'// &
              pthsep//'orog_fx_CNRM-CM5_historical_r0i0p0.nc'
  end subroutine find_cnrm_topo

  subroutine find_cnrm_file(cnrm_filename,var,idate)
    implicit none
    character(len=256) , intent(out) :: cnrm_filename
    character(len=*) , intent(in) :: var
    type(rcm_time_and_date) , intent(in) :: idate
    character(len=10) :: d1 , d2
    integer(ik4) :: y , m , d , h
    integer(ik4) :: y1 , y2 , m1 , m2
    call split_idate(idate,y,m,d,h)
    y1 = y
    y2 = y
    if ( d == 1 .and. h == 0 ) then
      m1 = m - 1
      if ( m1 < 1 ) then
        m1 = 12
        y1 = y1 - 1
      end if
    else
      m1 = m
    end if
    m2 = m1+1
    if ( m2 > 12 ) then
      m2 = 1
      y2 = y2+1
    end if
    write(d1,'(i0.4,i0.2,i0.2,i0.2)') y1, m1, 1, 6
    write(d2,'(i0.4,i0.2,i0.2,i0.2)') y2, m2, 1, 0
    if ( .not. date_in_scenario(idate,5,.true.) ) then
      call assemble_path(cnrm_filename,'RF',var,d1,d2)
    else
      call assemble_path(cnrm_filename,'RCP'//dattyp(4:5),var,d1,d2)
    end if
  end subroutine find_cnrm_file

end module mod_cnrm_helper
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

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

module mod_gfdl_helper

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_date

  private

  public :: gfdlvars
  public :: find_gfdl_sst
  public :: find_gfdl_dim , find_gfdl_topo , find_gfdl_file

  integer(ik4) , parameter :: nvars = 6
  character(len=3) , target , dimension(nvars) :: gfdlvars = &
            ['ta ' , 'XXX' , 'hus' , 'ua ' , 'va ' , 'ps ']

  character(len=64) :: gfdlbase1  = '_6hrLev_GFDL-ESM2M_historical'
  character(len=64) :: gfdlbase2  = '_6hrLev_GFDL-ESM2M_rcp'
  character(len=64) :: gfdlbase3 = '_r1i1p1_'

  contains

  subroutine find_gfdl_sst(fname,idate)
    implicit none
    character(len=256) , intent(out) :: fname
    type(rcm_time_and_date) , intent(in) :: idate
    integer(ik4) :: y1 , y2
    character(len=6) :: d1 , d2
    integer(ik4) :: y , m , d , h
    call split_idate(idate,y,m,d,h)
    y1 = (y-1)/5*5+1
    y2 = y1+4
    write(d1,'(i0.4,i0.2)') y1, 1
    write(d2,'(i0.4,i0.2)') y2, 12
    if ( y1 < 2005 ) then
      fname = trim(inpglob)//pthsep//'GFDL-ESM2M'//pthsep//'SST'// &
              pthsep//'ts_Amon_GFDL-ESM2M_historical'// &
              '_r1i1p1_'//d1//'-'//d2//'.nc'
    else
      fname = trim(inpglob)//pthsep//'GFDL-ESM2M'//pthsep//'SST'// &
              pthsep//'ts_Amon_GFDL-ESM2M_rcp'//ssttyp(4:5)//  &
              '_r1i1p1_'//d1//'-'//d2//'.nc'
    end if
  end subroutine find_gfdl_sst

  subroutine assemble_path(fname,scen,var,d1,d2)
    implicit none
    character(len=256) , intent(out) :: fname
    character(len=*) , intent(in) :: scen
    character(len=*) , intent(in) :: var
    character(len=*) , intent(in) :: d1
    character(len=*) , intent(in) :: d2
    if ( scen == 'RF' ) then
      fname = trim(inpglob)//pthsep//'GFDL-ESM2M'//pthsep//trim(scen)// &
              pthsep//trim(var)//pthsep//trim(var)//trim(gfdlbase1)//  &
              trim(gfdlbase3)//trim(d1)//'-'//trim(d2)//'.nc'
    else
      fname = trim(inpglob)//pthsep//'GFDL-ESM2M'//pthsep//trim(scen)// &
              pthsep//trim(var)//pthsep//trim(var)//trim(gfdlbase2)//  &
              scen(4:5)//trim(gfdlbase3)//trim(d1)//'-'//trim(d2)//'.nc'
    end if
  end subroutine assemble_path

  subroutine find_gfdl_dim(dim_filename)
    implicit none
    character(len=256) , intent(out) :: dim_filename
    ! Just return the name of one file in the historical dataset
    ! we hope is there.
    call assemble_path(dim_filename,'RF','ta','1951010100','1955123123')
  end subroutine find_gfdl_dim

  subroutine find_gfdl_topo(topo_filename)
    implicit none
    character(len=256) , intent(out) :: topo_filename
    topo_filename = trim(inpglob)//pthsep//'GFDL-ESM2M'//pthsep//'fixed'// &
              pthsep//'orog_fx_GFDL-ESM2M_historical_r0i0p0.nc'
  end subroutine find_gfdl_topo

  subroutine find_gfdl_file(gfdl_filename,var,idate)
    implicit none
    character(len=256) , intent(out) :: gfdl_filename
    character(len=*) , intent(in) :: var
    type(rcm_time_and_date) , intent(in) :: idate
    character(len=10) :: d1 , d2
    integer(ik4) :: y , m , d , h
    integer(ik4) :: y1 , y2
    call split_idate(idate,y,m,d,h)
    y1 = (y-1)/5*5+1
    y2 = y1+4
    if ( y == y1 .and. m == 1 .and. d == 1 .and. h == 0 ) then
      y1 = y1-5
      y2 = y2-5
    end if
    write(d1,'(i0.4,i0.2,i0.2,i0.2)') y1, 1, 1, 0
    write(d2,'(i0.4,i0.2,i0.2,i0.2)') y2, 12, 31, 23
    if ( .not. date_in_scenario(idate,5,.true.) .or. &
      (y == 2006 .and.  m == 1 .and. d == 1 .and. h == 0) ) then
      call assemble_path(gfdl_filename,'RF',var,d1,d2)
    else
      call assemble_path(gfdl_filename,'RCP'//dattyp(4:5),var,d1,d2)
    end if
  end subroutine find_gfdl_file

end module mod_gfdl_helper
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

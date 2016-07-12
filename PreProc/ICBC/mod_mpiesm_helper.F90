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

module mod_mpiesm_helper

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_date

  private

  public :: mpievars
  public :: find_mpiesm_sst
  public :: find_mpiesm_dim , find_mpiesm_topo , find_mpiesm_file

  integer(ik4) , parameter :: nvars = 6
  character(len=3) , target , dimension(nvars) :: mpievars = &
            (/'ta ' , 'XXX' , 'hus' , 'ua ' , 'va ' , 'aps'/)

  character(len=64) :: mpiebase1  = '_6hrLev_MPI-ESM-MR_historical'
  character(len=64) :: mpiebase2  = '_6hrLev_MPI-ESM-MR_rcp'
  character(len=64) :: mpiebase3 = '_r1i1p1_'

  contains

  subroutine find_mpiesm_sst(fname,idate)
    implicit none
    character(len=256) , intent(out) :: fname
    type(rcm_time_and_date) , intent(in) :: idate
    character(len=10) :: d1 , d2
    integer(ik4) :: y , m , d , h
    integer(ik4) :: y1 , y2 , m1 , m2
    call split_idate(idate,y,m,d,h)
    y1 = y
    m1 = m
    y2 = y1
    m2 = m1 + 1
    if ( m2 > 12 ) then
      m2 = 1
      y2 = y2 + 1
    end if
    write(d1,'(i0.4,i0.2,i0.2,i0.2)') y1, m1, 1, 0
    write(d2,'(i0.4,i0.2,i0.2,i0.2)') y2, m2, 1, 0
    if ( idate < 2005120100 ) then
      fname = trim(inpglob)//pthsep//'MPI-ESM-MR'//pthsep//'SST'// &
              pthsep//'tos_6hrLev_MPI-ESM-MR_historical'// &
              '_r1i1p1_'//d1//'00-'//d2//'00.nc'
    else
      fname = trim(inpglob)//pthsep//'MPI-ESM-MR'//pthsep//'SST'// &
              pthsep//'tos_6hrLev_MPI-ESM-MR_rcp'//ssttyp(4:5)//  &
              '_r1i1p1_'//d1//'00-'//d2//'00.nc'
    end if
  end subroutine find_mpiesm_sst

  subroutine assemble_path(fname,scen,var,d1,d2)
    implicit none
    character(len=256) , intent(out) :: fname
    character(len=*) , intent(in) :: scen
    character(len=*) , intent(in) :: var
    character(len=*) , intent(in) :: d1
    character(len=*) , intent(in) :: d2
    if ( scen == 'RF' ) then
      fname = trim(inpglob)//pthsep//'MPI-ESM-MR'//pthsep//trim(scen)// &
              pthsep//trim(var)//pthsep//trim(var)//trim(mpiebase1)//  &
              trim(mpiebase3)//trim(d1)//'00-'//trim(d2)//'00.nc'
    else
      fname = trim(inpglob)//pthsep//'MPI-ESM-MR'//pthsep//trim(scen)// &
              pthsep//trim(var)//pthsep//trim(var)//trim(mpiebase2)//  &
              scen(4:5)//trim(mpiebase3)//trim(d1)//'00-'//trim(d2)//'00.nc'
    end if
  end subroutine assemble_path

  subroutine find_mpiesm_dim(dim_filename)
    implicit none
    character(len=256) , intent(out) :: dim_filename
    ! Just return the name of one file in the historical dataset
    ! we hope is there.
    call assemble_path(dim_filename,'RF','ta','1970010100','1970020100')
  end subroutine find_mpiesm_dim

  subroutine find_mpiesm_topo(topo_filename)
    implicit none
    character(len=256) , intent(out) :: topo_filename
    topo_filename = trim(inpglob)//pthsep//'MPI-ESM-MR'//pthsep//'fixed'// &
              pthsep//'geosp_fx_MPI-ESM-MR_historical_r1i1p1.nc'
  end subroutine find_mpiesm_topo

  subroutine find_mpiesm_file(mpiesm_filename,var,idate)
    implicit none
    character(len=256) , intent(out) :: mpiesm_filename
    character(len=*) , intent(in) :: var
    type(rcm_time_and_date) , intent(in) :: idate
    character(len=10) :: d1 , d2
    integer(ik4) :: y , m , d , h
    integer(ik4) :: y1 , y2 , m1 , m2
    call split_idate(idate,y,m,d,h)
    y1 = y
    m1 = m
    y2 = y
    m2 = m1+1
    if ( m2 > 12 ) then
      m2 = 1
      y2 = y2+1
    end if
    write(d1,'(i0.4,i0.2,i0.2,i0.2)') y1, m1, 1, 0
    write(d2,'(i0.4,i0.2,i0.2,i0.2)') y2, m2, 1, 0
    if ( .not. date_in_scenario(idate,5,.true.) ) then
      call assemble_path(mpiesm_filename,'RF',var,d1,d2)
    else
      call assemble_path(mpiesm_filename,'RCP'//dattyp(4:5),var,d1,d2)
    end if
  end subroutine find_mpiesm_file

end module mod_mpiesm_helper
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

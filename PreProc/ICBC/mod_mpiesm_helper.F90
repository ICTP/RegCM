!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_mpiesm_helper

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_message
  use mod_date
  use mod_stdio

  private

  public :: mpievars
  public :: find_mpiesm_sst
  public :: find_mpiesm_dim, find_mpiesm_topo, find_mpiesm_file

  integer(ik4), parameter :: nvars = 6
  character(len=3), target, dimension(nvars) :: mpievars = &
            ['ta ', 'XXX', 'hus', 'ua ', 'va ', 'aps']

  character(len=64) :: mpiebase1  = '_6hrLev_MPI-ESM-MR_historical'
  character(len=64) :: mpiebase2  = '_6hrLev_MPI-ESM-MR_rcp'
  character(len=64) :: mpiebase3 = '_r1i1p1_'
  character(len=64) :: mpilbase1  = '_6hrLev_MPI-ESM-LR_historical'
  character(len=64) :: mpilbase2  = '_6hrLev_MPI-ESM-LR_rcp'
  character(len=64) :: mpilbase3 = '_r1i1p1_'

  contains

  subroutine find_mpiesm_sst(fname,idate,res)
    implicit none
    character(len=256), intent(out) :: fname
    type(rcm_time_and_date), intent(in) :: idate
    character(len=1), intent(in) :: res
    character(len=10) :: d1, d2
    integer(ik4) :: y, m, d, h
    integer(ik4) :: y1, y2, m1, m2
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
    if ( res == 'M' ) then
      if ( idate < 2005120100 ) then
        fname = trim(inpglob)//pthsep//'MPI-ESM-MR'//pthsep//'SST'// &
                pthsep//'tos_6hrLev_MPI-ESM-MR_historical'// &
                '_r1i1p1_'//d1//'00-'//d2//'00.nc'
      else
        fname = trim(inpglob)//pthsep//'MPI-ESM-MR'//pthsep//'SST'// &
                pthsep//'tos_6hrLev_MPI-ESM-MR_rcp'//ssttyp(4:5)//  &
                '_r1i1p1_'//d1//'00-'//d2//'00.nc'
      end if
    else if ( res == 'L' ) then
      if ( idate < 2005120100 ) then
        fname = trim(inpglob)//pthsep//'MPI-ESM-LR'//pthsep//'SST'// &
                pthsep//'tos_6hrLev_MPI-ESM-LR_historical'// &
                '_r1i1p1_'//d1//'00-'//d2//'00.nc'
      else
        fname = trim(inpglob)//pthsep//'MPI-ESM-LR'//pthsep//'SST'// &
                pthsep//'tos_6hrLev_MPI-ESM-LR_rcp'//ssttyp(4:5)//  &
                '_r1i1p1_'//d1//'00-'//d2//'00.nc'
      end if
    else
      write(stderr,*) 'Resolution requested: ', res
      write(stderr,*) 'Resolution supported: L,M'
      call die('sst','Unknown resolution for MPI-ESM',1)
    end if
  end subroutine find_mpiesm_sst

  subroutine assemble_path(fname,scen,var,d1,d2,res)
    implicit none
    character(len=256), intent(out) :: fname
    character(len=*), intent(in) :: scen
    character(len=*), intent(in) :: var
    character(len=*), intent(in) :: d1
    character(len=*), intent(in) :: d2
    character(len=1), intent(in) :: res
    if ( res == 'M' ) then
      if ( scen == 'RF' ) then
        fname = trim(inpglob)//pthsep//'MPI-ESM-MR'//pthsep//trim(scen)// &
                pthsep//trim(var)//pthsep//trim(var)//trim(mpiebase1)//  &
                trim(mpiebase3)//trim(d1)//'00-'//trim(d2)//'00.nc'
      else
        fname = trim(inpglob)//pthsep//'MPI-ESM-MR'//pthsep//trim(scen)// &
                pthsep//trim(var)//pthsep//trim(var)//trim(mpiebase2)//  &
                scen(4:5)//trim(mpiebase3)//trim(d1)//'00-'//trim(d2)//'00.nc'
      end if
    else if ( res == 'L' ) then
      if ( scen == 'RF' ) then
        fname = trim(inpglob)//pthsep//'MPI-ESM-LR'//pthsep//trim(scen)// &
                pthsep//trim(var)//pthsep//trim(var)//trim(mpilbase1)//  &
                trim(mpilbase3)//trim(d1)//'00-'//trim(d2)//'00.nc'
      else
        fname = trim(inpglob)//pthsep//'MPI-ESM-LR'//pthsep//trim(scen)// &
                pthsep//trim(var)//pthsep//trim(var)//trim(mpilbase2)//  &
                scen(4:5)//trim(mpilbase3)//trim(d1)//'00-'//trim(d2)//'00.nc'
      end if
    else
      write(stderr,*) 'Resolution requested: ', res
      write(stderr,*) 'Resolution supported: L,M'
      call die('sst','Unknown resolution for MPI-ESM',1)
    end if
  end subroutine assemble_path

  subroutine find_mpiesm_dim(dim_filename,res)
    implicit none
    character(len=256), intent(out) :: dim_filename
    character(len=1), intent(in) :: res
    ! Just return the name of one file in the historical dataset
    ! we hope is there.
    call assemble_path(dim_filename,'RF','ta','1970010100','1970020100',res)
  end subroutine find_mpiesm_dim

  subroutine find_mpiesm_topo(topo_filename,res)
    implicit none
    character(len=256), intent(out) :: topo_filename
    character(len=1), intent(in) :: res
    if ( res == 'M' ) then
      topo_filename = trim(inpglob)//pthsep//'MPI-ESM-MR'//pthsep//'fixed'// &
              pthsep//'geosp_fx_MPI-ESM-MR_historical_r1i1p1.nc'
    else if ( res == 'L' ) then
      topo_filename = trim(inpglob)//pthsep//'MPI-ESM-LR'//pthsep//'fixed'// &
              pthsep//'geosp_fx_MPI-ESM-LR_historical_r1i1p1.nc'
    else
      write(stderr,*) 'Resolution requested: ', res
      write(stderr,*) 'Resolution supported: L,M'
      call die('sst','Unknown resolution for MPI-ESM',1)
    end if
  end subroutine find_mpiesm_topo

  subroutine find_mpiesm_file(mpiesm_filename,var,idate,res)
    implicit none
    character(len=256), intent(out) :: mpiesm_filename
    character(len=*), intent(in) :: var
    type(rcm_time_and_date), intent(in) :: idate
    character(len=1), intent(in) :: res
    character(len=10) :: d1, d2
    integer(ik4) :: y, m, d, h
    integer(ik4) :: y1, y2, m1, m2
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
      call assemble_path(mpiesm_filename,'RF',var,d1,d2,res)
    else
      call assemble_path(mpiesm_filename,'RCP'//dattyp(4:5),var,d1,d2,res)
    end if
  end subroutine find_mpiesm_file

end module mod_mpiesm_helper
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

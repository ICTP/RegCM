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

module mod_miroc_helper

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_date

  private

  public :: mirocvars
  public :: find_miroc_sst
  public :: find_miroc_dim, find_miroc_topo, find_miroc_file

  integer(ik4), parameter :: nvars = 6
  character(len=3), target, dimension(nvars) :: mirocvars = &
            ['ta ', 'XXX', 'hus', 'ua ', 'va ', 'ps ']

  character(len=64) :: mirocbase1  = '_6hrLev_MIROC5_historical'
  character(len=64) :: mirocbase2  = '_6hrLev_MIROC5_rcp'
  character(len=64) :: mirocbase3 = '_r1i1p1_'

  integer(ik4), dimension(12) :: dm

  data dm /31, 28, 31, 30, 31, 30, 31, 31, 30, 31, 30, 31/

  contains

  subroutine find_miroc_sst(fname,idate)
    implicit none
    character(len=256), intent(out) :: fname
    type(rcm_time_and_date), intent(in) :: idate
    if ( .not. date_in_scenario(idate,5) ) then
      fname = trim(inpglob)//pthsep//'MIROC5'//pthsep//'SST'// &
              pthsep//'tos_Omon_MIROC5_historical'// &
              '_r1i1p1_185001-201212.nc'
    else
      fname = trim(inpglob)//pthsep//'MIROC5'//pthsep//'SST'// &
              pthsep//'tos_Omon_MIROC5_rcp'//ssttyp(4:5)//  &
              '_r1i1p1_200601-210012.nc'
    end if
  end subroutine find_miroc_sst

  subroutine assemble_path(fname,scen,var,d1,d2)
    implicit none
    character(len=256), intent(out) :: fname
    character(len=*), intent(in) :: scen
    character(len=*), intent(in) :: var
    character(len=*), intent(in) :: d1
    character(len=*), intent(in) :: d2
    if ( scen == 'RF' ) then
      fname = trim(inpglob)//pthsep//'MIROC5'//pthsep//trim(scen)// &
              pthsep//trim(var)//pthsep//trim(var)//trim(mirocbase1)//  &
              trim(mirocbase3)//trim(d1)//'-'//trim(d2)//'.nc'
    else
      fname = trim(inpglob)//pthsep//'MIROC5'//pthsep//trim(scen)// &
              pthsep//trim(var)//pthsep//trim(var)//trim(mirocbase2)//  &
              scen(4:5)//trim(mirocbase3)//trim(d1)//'-'//trim(d2)//'.nc'
    end if
  end subroutine assemble_path

  subroutine find_miroc_dim(dim_filename)
    implicit none
    character(len=256), intent(out) :: dim_filename
    ! Just return the name of one file in the historical dataset
    ! we hope is there.
    call assemble_path(dim_filename,'RF','ta','2000120100','2000123118')
  end subroutine find_miroc_dim

  subroutine find_miroc_topo(topo_filename)
    implicit none
    character(len=256), intent(out) :: topo_filename
    topo_filename = trim(inpglob)//pthsep//'MIROC5'//pthsep//'fixed'// &
              pthsep//'orog_fx_MIROC5_historical_r0i0p0.nc'
  end subroutine find_miroc_topo

  subroutine find_miroc_file(miroc_filename,var,idate)
    implicit none
    character(len=256), intent(out) :: miroc_filename
    character(len=*), intent(in) :: var
    type(rcm_time_and_date), intent(in) :: idate
    character(len=10) :: d1, d2
    integer(ik4) :: y, m, d, h
    call split_idate(idate,y,m,d,h)
    select case (var)
      case ('ps')
        write(d1,'(i0.4,i0.2,i0.2,i0.2)') y, 1, 1, 0
        write(d2,'(i0.4,i0.2,i0.2,i0.2)') y, 12, 31, 18
        if ( .not. date_in_scenario(idate,5,.true.) ) then
          call assemble_path(miroc_filename,'RF',var,d1,d2)
        else
          call assemble_path(miroc_filename,'RCP'//dattyp(4:5),var,d1,d2)
        end if
      case default
        write(d1,'(i0.4,i0.2,i0.2,i0.2)') y, m, 1, 0
        write(d2,'(i0.4,i0.2,i0.2,i0.2)') y, m, dm(m), 18
        if ( .not. date_in_scenario(idate,5,.true.) ) then
          call assemble_path(miroc_filename,'RF',var,d1,d2)
        else
          call assemble_path(miroc_filename,'RCP'//dattyp(4:5),var,d1,d2)
        end if
    end select
  end subroutine find_miroc_file

end module mod_miroc_helper
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

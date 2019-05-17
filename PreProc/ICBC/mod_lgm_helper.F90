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

module mod_lgm_helper

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_message
  use mod_date

  private

  public :: lgmvars
  public :: find_lgm_sst
  public :: find_lgm_dim , find_lgm_topo , find_lgm_file

  integer(ik4) , parameter :: nvars = 6
  character(len=3) , target , dimension(nvars) :: lgmvars = &
            ['t  ' , 'XXX' , 'sq ' , 'u  ' , 'v  ' , 'aps']

  contains

  subroutine find_lgm_sst(fname,idate,res)
    implicit none
    character(len=256) , intent(out) :: fname
    type(rcm_time_and_date) , intent(in) :: idate
    character(len=1) , intent(in) :: res
    character(len=4) :: cy
    integer(ik4) :: y , m , d , h
    call split_idate(idate,y,m,d,h)
    write(cy,'(i0.4)') y
    if ( res == 'P' ) then
      fname = trim(inpglob)//pthsep//'LGM'//pthsep//'GCM_LGM_mod'// &
                  pthsep//'lgm_r1i1p1-P_echam6_echam_'// &
                  cy//'_mod.nc'
    else
      if ( y > 1939 ) then
        fname = trim(inpglob)//pthsep//'LGM'//pthsep//'GCM_PiControl_mod'// &
                    pthsep//'piControl_r1i1p1-LR_echam6_echam_'// &
                    cy//'_mod.nc'
      else
        fname = trim(inpglob)//pthsep//'LGM'//pthsep//'GCM_PiControl_mod'// &
                    pthsep//'piControl_r1i1p1-P_echam6_echam_'// &
                    cy//'_mod.nc'
      end if
    end if
  end subroutine find_lgm_sst

  subroutine find_lgm_dim(fname,res)
    implicit none
    character(len=256) , intent(out) :: fname
    character(len=1) , intent(in) :: res
    character(len=4) , parameter :: cy = '1930'
    ! Just return the name of one file in the historical dataset
    ! we hope is there.
    if ( res == 'P' ) then
      fname = trim(inpglob)//pthsep//'LGM'//pthsep//'GCM_LGM_mod'// &
                  pthsep//'lgm_r1i1p1-P_echam6_echam_'// &
                  cy//'_mod.nc'
    else
      fname = trim(inpglob)//pthsep//'LGM'//pthsep//'GCM_PiControl_mod'// &
                  pthsep//'piControl_r1i1p1-P_echam6_echam_'// &
                  cy//'_mod.nc'
    end if
  end subroutine find_lgm_dim

  subroutine find_lgm_topo(topo_filename,res)
    implicit none
    character(len=256) , intent(out) :: topo_filename
    character(len=1) , intent(in) :: res
    if ( res == 'P' ) then
      topo_filename = trim(inpglob)//pthsep//'LGM'//pthsep//'fixed'// &
              pthsep//'lgm_r1i1p1-P_echam6_echam_fx_geosp.nc'
    else
      topo_filename = trim(inpglob)//pthsep//'LGM'//pthsep//'fixed'// &
              pthsep//'piControl_r1i1p1-LR_echam6_echam_fx_geosp.nc'
    end if
  end subroutine find_lgm_topo

  subroutine find_lgm_file(fname,idate,res)
    implicit none
    character(len=256) , intent(out) :: fname
    type(rcm_time_and_date) , intent(in) :: idate
    character(len=1) , intent(in) :: res
    character(len=4) :: cy
    integer(ik4) :: y , m , d , h
    call split_idate(idate,y,m,d,h)
    write(cy,'(i0.4)') y
    if ( res == 'P' ) then
      fname = trim(inpglob)//pthsep//'LGM'//pthsep//'GCM_LGM_mod'// &
                  pthsep//'lgm_r1i1p1-P_echam6_echam_'// &
                  cy//'_mod.nc'
    else
      if ( y > 1939 ) then
        fname = trim(inpglob)//pthsep//'LGM'//pthsep//'GCM_PiControl_mod'// &
                    pthsep//'piControl_r1i1p1-LR_echam6_echam_'// &
                    cy//'_mod.nc'
      else
        fname = trim(inpglob)//pthsep//'LGM'//pthsep//'GCM_PiControl_mod'// &
                    pthsep//'piControl_r1i1p1-P_echam6_echam_'// &
                    cy//'_mod.nc'
      end if
    end if
  end subroutine find_lgm_file

end module mod_lgm_helper
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

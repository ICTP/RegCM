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
module mod_mkwetland
  use mod_realkinds
  use mod_intkinds
  use mod_dynparam
  use mod_constants
  use mod_grid
  use mod_rdldtr
  use mod_intldtr
  use mod_message
  use mod_memutil
  use netcdf

  implicit none

  private

  public :: mkwetland

  character(len=16) , parameter :: varname1 = 'PCT_WETLAND'
  character(len=16) , parameter :: varname2 = 'PCT_LAKE'

  real(rkx) , parameter :: vcutoff = 1.0_rkx

  contains

  subroutine mkwetland(wetfile,mask,wetland,lake)
    implicit none
    character(len=*) , intent(in) :: wetfile
    real(rkx) , dimension(:,:) , intent(in) :: mask
    real(rkx) , dimension(:,:) , intent(out) :: wetland , lake
    integer(ik4) :: i , j
    type(globalfile) :: gfile
    character(len=256) :: inpfile

    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//wetfile
    call gfopen(gfile,inpfile,xlat,xlon,ds*nsg,roidem,i_band)
    call gfread(gfile,varname1,wetland,d_zero)
    call gfread(gfile,varname2,lake,d_zero)
    call gfclose(gfile)

    do i = 1 , iysg
      do j = 1 , jxsg
        if ( mask(j,i) < 0.5_rkx ) then
          wetland(j,i) = h_missing_value
        else
          if ( wetland(j,i) < d_zero ) then
            wetland(j,i) = d_zero
          else
            if ( wetland(j,i) < vcutoff ) then
              wetland(j,i) = d_zero
            else
              wetland(j,i) = min(max(0,nint(wetland(j,i))),100)
            end if
          end if
        end if
      end do
    end do
    do i = 1 , iysg
      do j = 1 , jxsg
        if ( mask(j,i) < 0.5_rkx ) then
          lake(j,i) = h_missing_value
        else
          if ( lake(j,i) < d_zero ) then
            lake(j,i) = d_zero
          else
            if ( lake(j,i) < vcutoff ) then
              lake(j,i) = d_zero
            else
              lake(j,i) = min(max(0,nint(lake(j,i))),100)
            end if
          end if
        end if
      end do
    end do
  end subroutine mkwetland

end module mod_mkwetland
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

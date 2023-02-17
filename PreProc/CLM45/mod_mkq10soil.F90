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
module mod_mkq10soil
#if defined(CN)
  use mod_realkinds
  use mod_intkinds
  use mod_constants
  use mod_dynparam
  use mod_grid
  use mod_rdldtr

  implicit none

  private

  public :: mkq10soil

  character(len=16) , parameter :: varname = 'q10'

  real(rkx) :: vmin = 0.0_rkx

  contains

  subroutine mkq10soil(q10soilfile,mask,q10soil)
    implicit none
    character(len=*) , intent(in) :: q10soilfile
    real(rkx) , dimension(:,:) , intent(in) :: mask
    real(rkx) , dimension(:,:) , intent(out) :: q10soil
    integer(ik4) :: i , j
    type(globalfile) :: gfile
    character(len=256) :: inpfile

    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//q10soilfile
    call gfopen(gfile,inpfile,xlat,xlon,ds*nsg,roidem,i_band)
    call gfread(gfile,varname,q10soil,h_missing_value)
    call gfclose(gfile)

    call bestaround(q10soil,h_missing_value)
    do i = 1 , iysg
      do j = 1 , jxsg
        if ( mask(j,i) < 0.5_rkx ) then
          q10soil(j,i) = h_missing_value
        else
          q10soil(j,i) = max(d_zero,q10soil(j,i))
        end if
      end do
    end do
  end subroutine mkq10soil
#endif
end module mod_mkq10soil
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

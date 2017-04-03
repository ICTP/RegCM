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
module mod_mksoilcol
  use mod_realkinds
  use mod_intkinds
  use mod_dynparam
  use mod_grid
  use mod_rdldtr

  implicit none

  private

  public :: mksoilcol

  character(len=16) , parameter :: varname = 'SOIL_COLOR'

  contains

  subroutine mksoilcol(soilcolfile,mask,soilcol)
    implicit none
    character(len=*) , intent(in) :: soilcolfile
    real(rkx) , dimension(:,:) , intent(in) :: mask
    real(rkx) , dimension(:,:) , intent(out) :: soilcol
    integer(ik4) :: i , j
    type(globalfile) :: gfile
    character(len=256) :: inpfile

    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//soilcolfile
    call gfopen(gfile,inpfile,xlat,xlon,ds*nsg,i_band)
    call gfread(gfile,varname,soilcol,h_missing_value)
    call gfclose(gfile)

    do i = 1 , iysg
      do j = 1 , jxsg
        if ( mask(j,i) < 0.5_rkx ) then
          soilcol(j,i) = h_missing_value
        else
          if ( soilcol(j,i) > h_missing_value ) then
            soilcol(j,i) = max(d_zero,soilcol(j,i))
          else
            call bestaround(soilcol,i,j)
            soilcol(j,i) = max(d_zero,soilcol(j,i))
          end if
        end if
      end do
    end do
  end subroutine mksoilcol

end module mod_mksoilcol
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

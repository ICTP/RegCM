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
module mod_mkglacier
  use mod_realkinds
  use mod_intkinds
  use mod_dynparam
  use mod_constants
  use mod_grid
  use mod_rdldtr

  implicit none

  private

  public :: mkglacier

  character(len=16) , parameter :: varname = 'PCT_GLACIER'

  real(rkx) , parameter :: vcutoff = 1.0_rkx

  contains

  subroutine mkglacier(glcfile,mask,glc)
    implicit none
    character(len=*) , intent(in) :: glcfile
    real(rkx) , dimension(:,:) , intent(in) :: mask
    real(rkx) , dimension(:,:) , intent(out) :: glc
    integer(ik4) :: i , j
    type(globalfile) :: gfile
    character(len=256) :: inpfile

    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//glcfile
    call gfopen(gfile,inpfile,xlat,xlon,ds*nsg,roidem,i_band)
    call gfread(gfile,varname,glc,d_zero)
    call gfclose(gfile)

    do i = 1 , iysg
      do j = 1 , jxsg
        if ( mask(j,i) < 0.5_rkx ) then
          glc(j,i) = h_missing_value
        else
          if ( glc(j,i) < d_zero ) then
            glc(j,i) = d_zero
          else
            if ( glc(j,i) < vcutoff ) then
              glc(j,i) = d_zero
            else
              glc(j,i) = min(max(0,nint(glc(j,i))),100)
            end if
          end if
        end if
      end do
    end do
  end subroutine mkglacier

end module mod_mkglacier
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

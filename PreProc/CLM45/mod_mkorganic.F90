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
module mod_mkorganic
  use mod_realkinds
  use mod_intkinds
  use mod_constants
  use mod_dynparam
  use mod_grid
  use mod_rdldtr

  implicit none

  private

  public :: mkorganic

  character(len=16) , parameter :: varname = 'ORGANIC'

  contains

  subroutine mkorganic(orgfile,mask,organic)
    implicit none
    character(len=*) , intent(in) :: orgfile
    real(rkx) , dimension(:,:) , intent(in) :: mask
    real(rkx) , dimension(:,:,:) , intent(out) :: organic
    integer(ik4) :: n , j , i , norg
    type(globalfile) :: gfile
    character(len=256) :: inpfile

    norg = size(organic,3)
    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//orgfile
    call gfopen(gfile,inpfile,xlat,xlon,ds*nsg,roidem,i_band)
    call gfread(gfile,varname,organic,h_missing_value)
    call gfclose(gfile)
    do n = 1 , norg
      call bestaround(organic(:,:,n),h_missing_value)
      do i = 1 , iysg
        do j = 1 , jxsg
          if ( mask(j,i) < 0.5_rkx ) then
            organic(j,i,n) = h_missing_value
          else
            organic(j,i,n) = min(max(d_zero,organic(j,i,n)),130.0_rkx)
          end if
        end do
      end do
    end do
  end subroutine mkorganic

end module mod_mkorganic
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

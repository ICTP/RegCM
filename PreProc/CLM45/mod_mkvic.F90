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
module mod_mkvic
#ifdef VICHYDRO
  use mod_realkinds
  use mod_intkinds
  use mod_dynparam
  use mod_constants
  use mod_grid
  use mod_rdldtr

  implicit none

  private

  public :: mkvic

  integer , parameter :: nvic = 4

  character(len=16) , parameter , dimension(nvic):: varname = &
          [ 'binfl' , 'Ds   ' , 'Dsmax' , 'Ws   ']

  contains

  subroutine mkvic(vicfile,mask,vic)
    implicit none
    character(len=*) , intent(in) :: vicfile
    real(rkx) , dimension(:,:) , intent(in) :: mask
    real(rkx) , dimension(:,:,:) , intent(out) :: vic
    integer(ik4) :: n , j , i
    type(globalfile) :: gfile
    character(len=256) :: inpfile

    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//vicfile

    call gfopen(gfile,inpfile,xlat,xlon,ds*nsg,roidem,i_band)
    do n = 1 , nvic
      call gfread(gfile,varname(n),vic(:,:,n),h_missing_value)
      call bestaround(vic(:,:,n),h_missing_value)
      do i = 1 , iysg
        do j = 1 , jxsg
          if ( mask(j,i) < 0.5_rkx ) then
            vic(j,i,n) = h_missing_value
          else
            vic(j,i,n) = max(d_zero,vic(j,i,n))
          end if
        end do
      end do
    end do
    call gfclose(gfile)
    do i = 1 , iysg
      do j = 1 , jxsg
        if ( mask(j,i) > 0.5_rkx ) then
          vic(j,i,1) = max(0.005_rkx,vic(j,i,1))
          vic(j,i,2) = max(0.1_rkx,vic(j,i,2))
          vic(j,i,3) = max(0.5_rkx,vic(j,i,3))
          vic(j,i,4) = max(0.5_rkx,vic(j,i,4))
        end if
      end do
    end do
  end subroutine mkvic
#endif

end module mod_mkvic
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

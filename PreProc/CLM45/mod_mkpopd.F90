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
module mod_mkpopd
#ifdef CN
  use mod_realkinds
  use mod_intkinds
  use mod_dynparam
  use mod_grid
  use mod_rdldtr

  implicit none

  private

  public :: mkpopd

  character(len=16) , parameter :: maskname = 'LANDMASK'
  character(len=16) , parameter :: varname = 'PDENS'

  real(rkx) :: vmisdat = -9999.0_rkx

  contains

  subroutine mkpopd(popdfile,popd,it)
    implicit none
    character(len=*) , intent(in) :: popdfile
    real(rkx) , dimension(:,:) , intent(out) :: popd
    integer(ik4) , intent(in) :: it
    integer(ik4) :: i , j
    real(rkx) , pointer , dimension(:,:) :: mask
    type(globalfile) :: gfile

    character(len=256) :: inpfile

    allocate(mask(jxsg,iysg))
    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//popdfile
    call gfopen(gfile,inpfile,xlat,xlon,ds*nsg,i_band)
    call gfread(gfile,maskname,mask)
    call gfread(gfile,varname,popd,it)

    do i = 1 , iysg
      do j = 1 , jxsg
        if ( mask(j,i) < 1.0_rkx ) then
          popd(j,i) = vmisdat
        else
          if ( popd(j,i) < 0.1_rkx ) popd(j,i) = 0.1_rkx
        end if
      end do
    end do
    deallocate(mask)
    call gfclose(gfile)
  end subroutine mkpopd
#endif

end module mod_mkpopd
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

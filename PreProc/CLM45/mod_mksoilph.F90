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
module mod_mksoilph
#if defined(CN) && defined(LCH4)
  use mod_realkinds
  use mod_intkinds
  use mod_dynparam
  use mod_grid
  use mod_rdldtr

  implicit none

  private

  public :: mksoilph

  character(len=16) , parameter :: varname = 'soilph'
  character(len=16) , parameter :: maskname = 'landmask'

  real(rkx) :: vmin = 4.4855_rkx
  real(rkx) :: vmisdat = 0.0_rkx

  contains

  subroutine mksoilph(soilphfile,soilph)
    implicit none
    character(len=*) , intent(in) :: soilphfile
    real(rkx) , dimension(:,:) , intent(out) :: soilph
    integer(ik4) :: i , j
    real(rkx) , pointer , dimension(:,:) :: mask
    type(globalfile) :: gfile

    character(len=256) :: inpfile

    allocate(mask(jxsg,iysg))
    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//soilphfile
    call gfopen(gfile,inpfile,xlat,xlon,ds*nsg,i_band)
    call gfread(gfile,maskname,mask)
    call gfread(gfile,varname,soilph)

    do i = 1 , iysg
      do j = 1 , jxsg
        if ( mask(j,i) < 1.0_rkx ) then
          soilph(j,i) = vmisdat
        else
          if ( soilph(j,i) < vmin ) soilph(j,i) = vmin
        end if
      end do
    end do
    deallocate(mask)
    call gfclose(gfile)
  end subroutine mksoilph
#endif

end module mod_mksoilph
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

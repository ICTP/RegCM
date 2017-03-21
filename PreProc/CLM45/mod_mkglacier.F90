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
  use mod_grid
  use mod_rdldtr
  use mod_intldtr
  use mod_message
  use mod_memutil

  implicit none

  private

  public :: mkglacier

  character(len=16) , parameter :: varname = 'PCT_GLACIER'
  character(len=16) , parameter :: maskname = 'LANDMASK'

  real(rkx) :: vmisdat = -9999.0_rkx
  real(rkx) :: vcutoff = 25.0_rkx

  contains

  subroutine mkglacier(glcfile,glc)
    implicit none
    character(len=*) , intent(in) :: glcfile
    real(rkx) , dimension(:,:) , intent(out) :: glc
    integer(ik4) :: i , j
    real(rkx) , pointer , dimension(:,:) :: rvar
    real(rkx) , pointer , dimension(:,:) :: rmask , mask

    character(len=256) :: inpfile

    allocate(mask(jxsg,iysg))
    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//glcfile
    call read_ncglob(inpfile,varname,0,4,i_band,xlat,xlon, &
                     grdlnma,grdlnmn,grdltma,grdltmn,nlatin,nlonin,rvar)
    call read_ncglob(inpfile,maskname,0,4,i_band,xlat,xlon, &
                     grdlnma,grdlnmn,grdltma,grdltmn,nlatin,nlonin,rmask)

    glc = 0.0_rkx
    call interp(ds*nsg,jxsg,iysg,xlat,xlon,mask,rmask,3)
    call interp(ds*nsg,jxsg,iysg,xlat,xlon,glc,rvar,6,rdem=roidem)
    do i = 1 , iysg
      do j = 1 , jxsg
        if ( mask(j,i) < 1.0_rkx ) then
          glc(j,i) = vmisdat
        else
          if ( glc(j,i) < vcutoff ) glc(j,i) = d_zero
        end if
      end do
    end do

    call relmem2d(rvar)
    call relmem2d(rmask)
  end subroutine mkglacier

end module mod_mkglacier
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

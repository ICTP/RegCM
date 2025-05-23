!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_mklake
  use mod_realkinds
  use mod_intkinds
  use mod_constants
  use mod_dynparam
  use mod_grid
  use mod_rdldtr

  implicit none

  private

  public :: mklake

  character(len=16), parameter :: varname = 'LAKEDEPTH'

  contains

  subroutine mklake(lakefile,mask,lake)
    implicit none
    character(len=*), intent(in) :: lakefile
    real(rkx), dimension(:,:), intent(in) :: mask
    real(rkx), dimension(:,:), intent(out) :: lake
    integer(ik4) :: i, j
    type(globalfile) :: gfile
    character(len=256) :: inpfile

    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//lakefile
    call gfopen(gfile,inpfile,xlat,xlon,ds*nsg,roidem,i_band)
    call gfread(gfile,varname,lake,h_missing_value)
    call gfclose(gfile)
    call bestaround(lake,h_missing_value)
    do i = 1, iysg
      do j = 1, jxsg
        if ( mask(j,i) < 0.5_rkx ) then
          lake(j,i) = h_missing_value
        else
          lake(j,i) = max(d_zero,lake(j,i))
        end if
      end do
    end do
  end subroutine mklake

end module mod_mklake
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

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

module mod_mkrough
  use mod_realkinds
  use mod_intkinds
  use mod_dynparam
  use mod_constants
  use mod_grid
  use mod_rdldtr

  implicit none

  private

  public :: do_roughness, mkrough

  character(len=16), parameter :: varname = 'z0'

  contains

  logical function do_roughness(roughfile) result(doit)
    implicit none
    character(len=*), intent(in) :: roughfile
    character(len=256) :: inpfile
    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//roughfile
    inquire(file=inpfile,exist=doit)
  end function do_roughness

  subroutine mkrough(roughfile,mask,rough)
    implicit none
    character(len=*), intent(in) :: roughfile
    real(rkx), dimension(:,:), intent(in) :: mask
    real(rkx), dimension(:,:), intent(out) :: rough
    integer(ik4) :: i, j
    type(globalfile) :: gfile
    character(len=256) :: inpfile

    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//roughfile
    call gfopen(gfile,inpfile,xlat,xlon,ds/nsg,roidem,i_band)
    call gfread(gfile,varname,rough,h_missing_value)
    call gfclose(gfile)

    call bestaround(rough,h_missing_value)

    do i = 1, iysg
      do j = 1, jxsg
        if ( mask(j,i) < 0.5_rkx ) then
          rough(j,i) = h_missing_value
        else
          rough(j,i) = max(d_zero,rough(j,i))
        end if
      end do
    end do
  end subroutine mkrough

end module mod_mkrough
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

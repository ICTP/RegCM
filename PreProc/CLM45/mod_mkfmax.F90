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

module mod_mkfmax
  use mod_realkinds
  use mod_intkinds
  use mod_dynparam
  use mod_constants
  use mod_grid
  use mod_rdldtr

  implicit none

  private

  public :: mkfmax

  character(len=16), parameter :: varname = 'FMAX'

  contains

  subroutine mkfmax(fmaxfile,mask,fmax)
    implicit none
    character(len=*), intent(in) :: fmaxfile
    real(rkx), dimension(:,:), intent(in) :: mask
    real(rkx), dimension(:,:), intent(out) :: fmax
    integer(ik4) :: j, i
    type(globalfile) :: gfile

    character(len=256) :: inpfile

    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//fmaxfile
    call gfopen(gfile,inpfile,xlat,xlon,ds*nsg,roidem,i_band)
    call gfread(gfile,varname,fmax,h_missing_value)
    call gfclose(gfile)

    call bestaround(fmax,h_missing_value)

    do i = 1, iysg
      do j = 1, jxsg
        if ( mask(j,i) < 0.5_rkx ) then
          fmax(j,i) = h_missing_value
        else
          fmax(j,i) = max(d_zero,min(d_one,fmax(j,i)))
        end if
      end do
    end do
  end subroutine mkfmax

end module mod_mkfmax
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

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

module mod_mklightning
#ifdef CN
  use mod_realkinds
  use mod_intkinds
  use mod_constants
  use mod_dynparam
  use mod_grid
  use mod_rdldtr

  implicit none

  private

  public :: mklightning, mklightning_init, mklightning_close

  character(len=16), parameter :: varname = 'lnfm'

  type(globalfile) :: gfile

  contains

  subroutine mklightning_init(lnfmfile)
    implicit none
    character(len=*), intent(in) :: lnfmfile
    character(len=256) :: inpfile
    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//lnfmfile
    call gfopen(gfile,inpfile,xlat,xlon,ds*nsg,roidem,i_band)
  end subroutine mklightning_init

  subroutine mklightning(lightning,mask,it)
    implicit none
    real(rkx), dimension(:,:), intent(in) :: mask
    real(rkx), dimension(:,:), intent(out) :: lightning
    integer(ik4), intent(in) :: it
    integer(ik4) :: i, j
    call gfread(gfile,varname,lightning,it,h_missing_value)
    call bestaround(lightning,h_missing_value)
    do i = 1, iysg
      do j = 1, jxsg
        if ( mask(j,i) < 0.5_rkx ) then
          lightning(j,i) = h_missing_value
        else
          lightning(j,i) = max(d_zero,lightning(j,i))
        end if
      end do
    end do
  end subroutine mklightning

  subroutine mklightning_close
    implicit none
    call gfclose(gfile)
  end subroutine mklightning_close
#endif

end module mod_mklightning
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

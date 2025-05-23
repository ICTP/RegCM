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

module mod_mkpopd
#ifdef CN
  use mod_realkinds
  use mod_intkinds
  use mod_constants
  use mod_dynparam
  use mod_grid
  use mod_rdldtr

  implicit none

  private

  public :: mkpopd, mkpopd_init, mkpopd_close

  character(len=16), parameter :: varname = 'PDENS'

  type(globalfile) :: gfile

  contains

  subroutine mkpopd_init(popdfile)
    implicit none
    character(len=*), intent(in) :: popdfile
    character(len=256) :: inpfile

    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//popdfile
    call gfopen(gfile,inpfile,xlat,xlon,ds*nsg,roidem,i_band)
  end subroutine mkpopd_init

  subroutine mkpopd(popd,mask,it)
    implicit none
    real(rkx), dimension(:,:), intent(in) :: mask
    real(rkx), dimension(:,:), intent(out) :: popd
    integer(ik4), intent(in) :: it
    integer(ik4) :: i, j
    call gfread(gfile,varname,popd,it,h_missing_value)
    call bestaround(popd,h_missing_value)
    do i = 1, iysg
      do j = 1, jxsg
        if ( mask(j,i) < 0.5_rkx ) then
          popd(j,i) = h_missing_value
        else
          popd(j,i) = max(d_zero,popd(j,i))
        end if
      end do
    end do
  end subroutine mkpopd

  subroutine mkpopd_close
    implicit none
    call gfclose(gfile)
  end subroutine mkpopd_close
#endif

end module mod_mkpopd
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

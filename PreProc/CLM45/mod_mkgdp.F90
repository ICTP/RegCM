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

module mod_mkgdp
  use mod_realkinds
  use mod_intkinds
  use mod_dynparam
  use mod_constants
  use mod_grid
  use mod_rdldtr

  implicit none

  private

  public :: mkgdp

  character(len=16), parameter :: varname = 'gdp'

  contains

  subroutine mkgdp(gdpfile,mask,gdp)
    implicit none
    character(len=*), intent(in) :: gdpfile
    real(rkx), dimension(:,:), intent(in) :: mask
    real(rkx), dimension(:,:), intent(out) :: gdp
    integer(ik4) :: i, j
    type(globalfile) :: gfile
    character(len=256) :: inpfile

    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//gdpfile
    call gfopen(gfile,inpfile,xlat,xlon,ds*nsg,roidem,i_band)
    call gfread(gfile,varname,gdp,h_missing_value)
    call gfclose(gfile)

    call bestaround(gdp,h_missing_value)

    do i = 1, iysg
      do j = 1, jxsg
        if ( mask(j,i) < 0.5_rkx ) then
          gdp(j,i) = h_missing_value
        else
          gdp(j,i) = max(d_zero,gdp(j,i))
        end if
      end do
    end do
  end subroutine mkgdp

end module mod_mkgdp
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

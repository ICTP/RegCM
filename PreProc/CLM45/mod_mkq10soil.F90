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

#if defined(CN)
module mod_mkq10soil
  use mod_realkinds, only : rkx
  use mod_intkinds, only : ik4
  use mod_constants, only : d_zero
  use mod_dynparam, only : ds, i_band, inpglob, iysg, jxsg, nsg, pthsep, &
    roidem
  use mod_grid, only : h_missing_value, xlat, xlon
  use mod_rdldtr, only : bestaround, gfclose, gfopen, gfread, globalfile

  implicit none

  private

  public :: mkq10soil

  character(len=16), parameter :: varname = 'q10'

  real(rkx) :: vmin = 0.0_rkx

  contains

  subroutine mkq10soil(q10soilfile,mask,q10soil)
    implicit none
    character(len=*), intent(in) :: q10soilfile
    real(rkx), dimension(:,:), intent(in) :: mask
    real(rkx), dimension(:,:), intent(out) :: q10soil
    integer(ik4) :: i, j
    type(globalfile) :: gfile
    character(len=256) :: inpfile

    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//q10soilfile
    call gfopen(gfile,inpfile,xlat,xlon,ds/nsg,roidem,i_band)
    call gfread(gfile,varname,q10soil,h_missing_value)
    call gfclose(gfile)

    call bestaround(q10soil,h_missing_value)
    do i = 1, iysg
      do j = 1, jxsg
        if ( mask(j,i) < 0.5_rkx ) then
          q10soil(j,i) = h_missing_value
        else
          q10soil(j,i) = max(d_zero,q10soil(j,i))
        end if
      end do
    end do
  end subroutine mkq10soil
end module mod_mkq10soil
#else
module mod_mkq10soil
  implicit none
  private
end module mod_mkq10soil
#endif
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

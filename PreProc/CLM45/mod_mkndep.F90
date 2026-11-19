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
module mod_mkndep
  use mod_realkinds, only : rkx
  use mod_intkinds, only : ik4
  use mod_constants, only : d_zero
  use mod_dynparam, only : ds, i_band, inpglob, iysg, jxsg, nsg, pthsep, &
    roidem
  use mod_grid, only : h_missing_value, xlat, xlon
  use mod_rdldtr, only : bestaround, gfclose, gfopen, gfread, globalfile

  implicit none

  private

  public :: mkndep

  character(len=16), parameter :: varname = 'NDEP_year'

  real(rkx) :: vmin = 0.0_rkx

  contains

  subroutine mkndep(ndepfile,mask,ndep)
    implicit none
    character(len=*), intent(in) :: ndepfile
    real(rkx), dimension(:,:), intent(in) :: mask
    real(rkx), dimension(:,:), intent(out) :: ndep
    integer(ik4) :: i, j
    type(globalfile) :: gfile
    character(len=256) :: inpfile

    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//ndepfile
    call gfopen(gfile,inpfile,xlat,xlon,ds/nsg,roidem,i_band)
    call gfread(gfile,varname,ndep,h_missing_value)
    call gfclose(gfile)

    call bestaround(ndep,h_missing_value)

    do i = 1, iysg
      do j = 1, jxsg
        if ( mask(j,i) < 0.5_rkx ) then
          ndep(j,i) = h_missing_value
        else
          ndep(j,i) = max(d_zero,ndep(j,i))
        end if
      end do
    end do
  end subroutine mkndep
end module mod_mkndep
#else
module mod_mkndep
  implicit none
  private
end module mod_mkndep
#endif
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

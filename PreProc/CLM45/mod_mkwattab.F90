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

module mod_mkwattab
  !
  ! Initial soil water state for a CLM cold start: volumetric water of
  ! the hydrologically active soil layers and water table depth, from a
  ! spun up CLM4.5 global state (mksrf_water.nc).
  !
  use mod_realkinds
  use mod_intkinds
  use mod_constants
  use mod_dynparam
  use mod_grid
  use mod_rdldtr

  implicit none

  private

  public :: mkwattab

  character(len=16), parameter :: h2oname = 'H2OSOI'
  character(len=16), parameter :: zwtname = 'ZWT'

  contains

  subroutine mkwattab(watfile,mask,h2osoi,zwt)
    implicit none
    character(len=*), intent(in) :: watfile
    real(rkx), dimension(:,:), intent(in) :: mask
    real(rkx), dimension(:,:,:), intent(out) :: h2osoi
    real(rkx), dimension(:,:), intent(out) :: zwt
    integer(ik4) :: n, j, i, nlev
    type(globalfile) :: gfile
    character(len=256) :: inpfile

    nlev = size(h2osoi,3)
    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//watfile
    call gfopen(gfile,inpfile,xlat,xlon,ds/nsg,roidem,i_band)
    call gfread(gfile,h2oname,h2osoi,h_missing_value)
    call gfread(gfile,zwtname,zwt,h_missing_value)
    call gfclose(gfile)
    do n = 1, nlev
      call bestaround(h2osoi(:,:,n),h_missing_value)
      do i = 1, iysg
        do j = 1, jxsg
          if ( mask(j,i) < 0.5_rkx ) then
            h2osoi(j,i,n) = h_missing_value
          else
            h2osoi(j,i,n) = min(max(d_zero,h2osoi(j,i,n)),d_one)
          end if
        end do
      end do
    end do
    call bestaround(zwt,h_missing_value)
    do i = 1, iysg
      do j = 1, jxsg
        if ( mask(j,i) < 0.5_rkx ) then
          zwt(j,i) = h_missing_value
        else
          zwt(j,i) = max(d_zero,zwt(j,i))
        end if
      end do
    end do
  end subroutine mkwattab

end module mod_mkwattab
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

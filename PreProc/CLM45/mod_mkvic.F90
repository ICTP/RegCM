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

#ifdef VICHYDRO
module mod_mkvic
  use mod_realkinds, only : rkx
  use mod_intkinds, only : ik4
  use mod_dynparam, only : ds, i_band, inpglob, iysg, jxsg, nsg, pthsep, &
    roidem
  use mod_constants, only : d_zero
  use mod_grid, only : h_missing_value, xlat, xlon
  use mod_rdldtr, only : bestaround, gfclose, gfopen, gfread, globalfile

  implicit none

  private

  public :: mkvic

  integer, parameter :: nvic = 4

  character(len=16), parameter, dimension(nvic):: varname = &
          [ 'binfl', 'Ds   ', 'Dsmax', 'Ws   ']

  contains

  subroutine mkvic(vicfile,mask,vic)
    implicit none
    character(len=*), intent(in) :: vicfile
    real(rkx), dimension(:,:), intent(in) :: mask
    real(rkx), dimension(:,:,:), intent(out) :: vic
    integer(ik4) :: n, j, i
    type(globalfile) :: gfile
    character(len=256) :: inpfile

    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//vicfile

    call gfopen(gfile,inpfile,xlat,xlon,ds/nsg,roidem,i_band)
    do n = 1, nvic
      call gfread(gfile,varname(n),vic(:,:,n),h_missing_value)
      call bestaround(vic(:,:,n),h_missing_value)
      do i = 1, iysg
        do j = 1, jxsg
          if ( mask(j,i) < 0.5_rkx ) then
            vic(j,i,n) = h_missing_value
          else
            vic(j,i,n) = max(d_zero,vic(j,i,n))
          end if
        end do
      end do
    end do
    call gfclose(gfile)
    do i = 1, iysg
      do j = 1, jxsg
        if ( mask(j,i) > 0.5_rkx ) then
          vic(j,i,1) = max(0.005_rkx,vic(j,i,1))
          vic(j,i,2) = max(0.1_rkx,vic(j,i,2))
          vic(j,i,3) = max(0.5_rkx,vic(j,i,3))
          vic(j,i,4) = max(0.5_rkx,vic(j,i,4))
        end if
      end do
    end do
  end subroutine mkvic

end module mod_mkvic
#else
module mod_mkvic
  implicit none
  private
end module mod_mkvic
#endif
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

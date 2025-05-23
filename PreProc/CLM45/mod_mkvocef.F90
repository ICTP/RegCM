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

module mod_mkvocef
  use mod_realkinds
  use mod_intkinds
  use mod_constants
  use mod_dynparam
  use mod_grid
  use mod_rdldtr

  implicit none

  private

  public :: mkvocef

  integer, parameter :: nvocs = 6

  character(len=16), parameter, dimension(nvocs):: varname = &
          [ 'ef_btr', 'ef_crp', 'ef_fdt', &
             'ef_fet', 'ef_grs', 'ef_shr']

  contains

  subroutine mkvocef(vocfile,mask,vocef)
    implicit none
    character(len=*), intent(in) :: vocfile
    real(rkx), dimension(:,:), intent(in) :: mask
    real(rkx), dimension(:,:,:), intent(out) :: vocef
    integer(ik4) :: n, i, j
    type(globalfile) :: gfile
    character(len=256) :: inpfile

    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//vocfile
    call gfopen(gfile,inpfile,xlat,xlon,ds*nsg,roidem,i_band)
    do n = 1, nvocs
      call gfread(gfile,varname(n),vocef(:,:,n),h_missing_value)
    end do
    call gfclose(gfile)

    do n = 1, nvocs
      call bestaround(vocef(:,:,n),h_missing_value)
      do i = 1, iysg
        do j = 1, jxsg
          if ( mask(j,i) < 0.5_rkx ) then
            vocef(j,i,n) = h_missing_value
          else
            vocef(j,i,n) = max(d_zero,vocef(j,i,n))
          end if
        end do
      end do
    end do
  end subroutine mkvocef

end module mod_mkvocef
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

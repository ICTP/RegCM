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

module mod_mkalbedo
  use mod_realkinds
  use mod_intkinds
  use mod_dynparam
  use mod_constants
  use mod_grid
  use mod_rdldtr

  implicit none (type, external)

  private

  public :: mkalbedo

  character(len=16), parameter :: varname = 'albedo'

  contains

  subroutine mkalbedo(albedo_file,mask,albedo)
    implicit none (type, external)
    character(len=*), intent(in) :: albedo_file
    real(rkx), dimension(:,:), intent(in) :: mask
    real(rkx), dimension(:,:,:), intent(out) :: albedo
    type(globalfile) :: gfile
    character(len=256) :: inpfile

    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//albedo_file
    call gfopen(gfile,inpfile,xlat,xlon,ds/nsg,roidem,i_band)
    call gfread(gfile,varname,albedo,h_missing_value)
    call gfclose(gfile)

    call inrange(albedo,0.0_rkx,1.0_rkx)

    contains

      subroutine inrange(f,mi,ma)
        implicit none (type, external)
        real(rkx), dimension(:,:,:), intent(inout) :: f
        real(rkx), intent(in) :: mi, ma
        integer :: i, j, n
        do n = 1, 2
          call bestaround(f(:,:,n),h_missing_value)
          do i = 1, iysg
            do j = 1, jxsg
              if ( mask(j,i) < 0.5_rkx ) then
                f(j,i,n) = h_missing_value
              else
                f(j,i,n) = min(ma,max(mi,f(j,i,n)))
              end if
            end do
          end do
        end do
      end subroutine inrange

  end subroutine mkalbedo

end module mod_mkalbedo
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

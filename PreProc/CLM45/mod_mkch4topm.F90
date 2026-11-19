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

#if defined(CN) && defined(LCH4)
module mod_mkch4topm
  use mod_realkinds, only : rkx
  use mod_intkinds, only : ik4
  use mod_dynparam, only : ds, i_band, inpglob, iysg, jxsg, nsg, pthsep, &
    roidem
  use mod_constants, only : d_zero
  use mod_grid, only : h_missing_value, xlat, xlon
  use mod_rdldtr, only : bestaround, gfclose, gfopen, gfread, globalfile

  implicit none

  private

  public :: mkch4topm

  integer, parameter :: nlch4 = 4

  character(len=16), parameter, dimension(nlch4):: varname = &
          [ 'K_PAR ', 'XM_PAR', 'V_PAR ', 'MAXF  ' ]

  contains

  subroutine mkch4topm(ch4topmfile,mask,lch4)
    implicit none
    character(len=*), intent(in) :: ch4topmfile
    real(rkx), dimension(:,:), intent(in) :: mask
    real(rkx), dimension(:,:,:), intent(out) :: lch4
    integer(ik4) :: n, i, j
    type(globalfile) :: gfile

    character(len=256) :: inpfile

    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//ch4topmfile

    call gfopen(gfile,inpfile,xlat,xlon,ds/nsg,roidem,i_band)
    do n = 1, nlch4
      call gfread(gfile,varname(n),lch4(:,:,n),h_missing_value)
      call bestaround(lch4(:,:,n),h_missing_value)
      do i = 1, iysg
        do j = 1, jxsg
          if ( mask(j,i) < 0.5_rkx ) then
            lch4(j,i,n) = h_missing_value
          else
            lch4(j,i,n) = max(d_zero,lch4(j,i,n))
          end if
        end do
      end do
    end do
    call gfclose(gfile)
  end subroutine mkch4topm
end module mod_mkch4topm
#else
module mod_mkch4topm
  implicit none
  private
end module mod_mkch4topm
#endif

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

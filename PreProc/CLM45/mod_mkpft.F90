!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
module mod_mkpft
  use mod_realkinds
  use mod_intkinds
  use mod_dynparam
  use mod_grid
  use mod_rdldtr
  use mod_intldtr
  use mod_message
  use mod_memutil

  implicit none

  private

  public :: mkpft

  character(len=16) , parameter :: varname = 'PCT_PFT'
  character(len=16) , parameter :: maskname = 'LANDMASK'

  real(rkx) :: vmisdat = -9999.0_rkx
  real(rkx) :: vcutoff = 19.0_rkx

  contains

  subroutine mkpft(pftfile,pft)
    implicit none
    character(len=*) , intent(in) :: pftfile
    real(rkx) , dimension(:,:,:) , intent(out) :: pft
    integer(ik4) :: i , j , n , npft
    integer(ik4) , dimension(1) :: il
    real(rkx) , pointer , dimension(:,:) :: mask
    type(globalfile) :: gfile

    character(len=256) :: inpfile

    npft = size(pft,3)
    allocate(mask(jxsg,iysg))
    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//pftfile
    call gfopen(gfile,inpfile,xlat,xlon,ds*nsg,i_band)
    call gfread(gfile,maskname,mask)
    call gfread(gfile,varname,pft)

    do i = 1 , iysg
      do j = 1 , jxsg
        if ( mask(j,i) < 1.0_rkx ) then
          pft(j,i,:) = vmisdat
        else
          il = maxloc(pft(j,i,:))
          do n = 1 , npft
            if ( n == il(1) ) cycle
            if ( pft(j,i,n) < vcutoff ) then
              pft(j,i,il(1)) = pft(j,i,il(1)) + pft(j,i,n)
              pft(j,i,n) = d_zero
            end if
          end do
          do n = 1 , npft
            pft(j,i,n) = min(pft(j,i,n),100.0_rkx)
          end do
        end if
      end do
    end do
    deallocate(mask)
    call gfclose(gfile)
  end subroutine mkpft

end module mod_mkpft
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

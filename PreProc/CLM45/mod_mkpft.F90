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
    real(rkx) , pointer , dimension(:,:,:) :: rvar
    real(rkx) , pointer , dimension(:,:) :: rmask , mask

    character(len=256) :: inpfile

    allocate(mask(jxsg,iysg))
    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//pftfile
    call read_ncglob(inpfile,varname,0,4,i_band,xlat,xlon, &
                     grdlnma,grdlnmn,grdltma,grdltmn,nlatin,nlonin,rvar)
    call read_ncglob(inpfile,maskname,0,4,i_band,xlat,xlon, &
                     grdlnma,grdlnmn,grdltma,grdltmn,nlatin,nlonin,rmask)

    npft = size(pft,3)
    if ( size(rvar,3) /= npft ) then
      call die('mkpft','Requested pft differ from pft in file '// &
               trim(pftfile),1)
    end if

    call interp(ds*nsg,jxsg,iysg,xlat,xlon,mask,rmask,3)

    pft = 0.0_rkx
    do n = 1 , npft
      call interp(ds*nsg,jxsg,iysg,xlat,xlon,pft(:,:,n), &
                  rvar(:,:,n),6,rdem=roidem)
    end do
    do i = 1 , iysg
      do j = 1 , jxsg
        if ( mask(j,i) < 1.0_rkx ) then
          pft(j,i,:) = vmisdat
        else
          il = maxloc(pft(j,i,:))
          do n = 1 , npft
            if ( n == il(1) ) cycle
            if ( pft(j,i,n) < vcutoff ) then
              pft(j,i,n) = d_zero
              pft(j,i,il(1)) = pft(j,i,il(1)) + pft(j,i,n)
            end if
          end do
        end if
      end do
    end do

    call relmem3d(rvar)
    call relmem2d(rmask)

    deallocate(mask)
  end subroutine mkpft

end module mod_mkpft
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

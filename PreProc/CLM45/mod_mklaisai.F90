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
module mod_mklaisai
  use mod_realkinds
  use mod_intkinds
  use mod_dynparam
  use mod_constants
  use mod_grid
  use mod_rdldtr

  implicit none

  private

  public :: mklaisai

  character(len=16) , parameter :: timedim = 'time'
  character(len=16) , parameter :: varname1 = 'MONTHLY_LAI'
  character(len=16) , parameter :: varname2 = 'MONTHLY_SAI'
  character(len=32) , parameter :: varname3 = 'MONTHLY_HEIGHT_TOP'
  character(len=32) , parameter :: varname4 = 'MONTHLY_HEIGHT_BOT'

  contains

  subroutine mklaisai(laisaifile,mask, &
                  monthly_lai,monthly_sai,monthly_top,monthly_bot)
    implicit none
    character(len=*) , intent(in) :: laisaifile
    real(rkx) , dimension(:,:) , intent(in) :: mask
    real(rkx) , dimension(:,:,:,:) , intent(out) :: monthly_sai , monthly_lai
    real(rkx) , dimension(:,:,:,:) , intent(out) :: monthly_top , monthly_bot
    integer(ik4) :: nm , np
    type(globalfile) :: gfile
    character(len=256) :: inpfile

    np = size(monthly_lai,3)
    nm = size(monthly_lai,4)
    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//laisaifile
    call gfopen(gfile,inpfile,xlat,xlon,ds*nsg,roidem,i_band)
    call gfread(gfile,varname1,monthly_lai,h_missing_value)
    call gfread(gfile,varname2,monthly_sai,h_missing_value)
    call gfread(gfile,varname3,monthly_top,h_missing_value)
    call gfread(gfile,varname4,monthly_bot,h_missing_value)
    call gfclose(gfile)

    call inrange(monthly_lai,0.0_rkx,7.0_rkx)
    call inrange(monthly_sai,0.0_rkx,6.5_rkx)
    call inrange(monthly_top,0.0_rkx,35.0_rkx)
    call inrange(monthly_bot,0.0_rkx,11.5_rkx)

    contains

      subroutine inrange(f,mi,ma)
        implicit none
        real(rkx) , dimension(:,:,:,:) , intent(inout) :: f
        real(rkx) , intent(in) :: mi , ma
        integer :: i , j , p , n
        do n = 1 , nm
          do p = 1, np
            call bestaround(f(:,:,p,n),h_missing_value)
            do i = 1 , iysg
              do j = 1 , jxsg
                if ( mask(j,i) < 0.5_rkx ) then
                  f(j,i,p,n) = h_missing_value
                else
                  f(j,i,p,n) = real(int(f(j,i,p,n)*100),rkx)/d_100
                  f(j,i,p,n) = min(ma,max(mi,f(j,i,p,n)))
                end if
              end do
            end do
          end do
        end do
      end subroutine inrange

  end subroutine mklaisai

end module mod_mklaisai
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

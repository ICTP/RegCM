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
module mod_mklch4
#if defined(CN) && defined(LCH4)
  use mod_realkinds
  use mod_intkinds
  use mod_dynparam
  use mod_grid
  use mod_rdldtr

  implicit none

  private

  public :: mklch4

  integer , parameter :: nlch4 = 3

  character(len=16) , parameter , dimension(nlch4):: varname = &
          [ 'F0  ' , 'P3  ' , 'ZWT0']

  contains

  subroutine mklch4(lch4file,mask,lch4)
    implicit none
    character(len=*) , intent(in) :: lch4file
    real(rkx) , dimension(:,:) , intent(in) :: mask
    real(rkx) , dimension(:,:,:) , intent(out) :: lch4
    integer(ik4) :: n , i , j
    type(globalfile) :: gfile

    character(len=256) :: inpfile

    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//lch4file

    call gfopen(gfile,inpfile,xlat,xlon,ds*nsg,roidem,i_band)
    do n = 1 , nlch4
      call gfread(gfile,varname(n),lch4(:,:,n),h_missing_value)
      do i = 1 , iysg
        do j = 1 , jxsg
          if ( mask(j,i) < 0.5_rkx ) then
            lch4(j,i,n) = h_missing_value
          else
            if ( lch4(j,i,n) > h_missing_value ) then
              lch4(j,i,n) = max(d_zero,lch4(j,i,n))
            else
              call bestaround(lch4(:,:,n),i,j)
              lch4(j,i,n) = max(d_zero,lch4(j,i,n))
            end if
          end if
        end do
      end do
    end do
    call gfclose(gfile)
  end subroutine mklch4
#endif

end module mod_mklch4
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

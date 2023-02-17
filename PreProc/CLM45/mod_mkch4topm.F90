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
module mod_mkch4topm
#if defined(CN) && defined(LCH4)
  use mod_realkinds
  use mod_intkinds
  use mod_dynparam
  use mod_constants
  use mod_grid
  use mod_rdldtr

  implicit none

  private

  public :: mkch4topm

  integer , parameter :: nlch4 = 4

  character(len=16) , parameter , dimension(nlch4):: varname = &
          [ 'K_PAR ' , 'XM_PAR' , 'V_PAR ', 'MAXF  ' ]

  contains

  subroutine mkch4topm(ch4topmfile,mask,lch4)
    implicit none
    character(len=*) , intent(in) :: ch4topmfile
    real(rkx) , dimension(:,:) , intent(in) :: mask
    real(rkx) , dimension(:,:,:) , intent(out) :: lch4
    integer(ik4) :: n , i , j
    type(globalfile) :: gfile

    character(len=256) :: inpfile

    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//ch4topmfile

    call gfopen(gfile,inpfile,xlat,xlon,ds*nsg,roidem,i_band)
    do n = 1 , nlch4
      call gfread(gfile,varname(n),lch4(:,:,n),h_missing_value)
      call bestaround(lch4(:,:,n),h_missing_value)
      do i = 1 , iysg
        do j = 1 , jxsg
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
#endif

end module mod_mkch4topm
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

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
module mod_mksoitex
  use mod_realkinds
  use mod_intkinds
  use mod_constants
  use mod_dynparam
  use mod_grid
  use mod_rdldtr

  implicit none

  private

  public :: mksoitex

  character(len=16) , parameter :: mapname = 'MAPUNITS'
  character(len=24) , parameter :: mapdim = 'max_value_mapunit'
  character(len=16) , parameter :: varname1 = 'PCT_SAND'
  character(len=16) , parameter :: varname2 = 'PCT_CLAY'

  contains

  subroutine mksoitex(soitexfile,mask,sand,clay)
    implicit none
    character(len=*) , intent(in) :: soitexfile
    real(rkx) , dimension(:,:) , intent(in) :: mask
    real(rkx) , dimension(:,:,:) , intent(out) :: sand , clay
    integer(ik4) :: i , j , nc , n
    type(globalfile) :: gfile
    real(rkx) :: tsum
    character(len=256) :: inpfile

    nc = size(sand,3)

    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//soitexfile
    call gfopen(gfile,inpfile,xlat,xlon,ds*nsg,roidem,i_band)
    call gfread(gfile,varname1,mapdim,mapname,sand,.false.,h_missing_value)
    call gfread(gfile,varname2,mapdim,mapname,clay,.false.,h_missing_value)

    do n = 1 , nc
      do i = 1 , iysg
        do j = 1 , jxsg
          if ( mask(j,i) < 0.5_rkx ) then
            sand(j,i,n) = h_missing_value
            clay(j,i,n) = h_missing_value
          else
            sand(j,i,n) = int(sand(j,i,n))
            clay(j,i,n) = int(clay(j,i,n))
            if ( sand(j,i,n) < 1.0 .and. clay(j,i,n) < 100.0 ) then
              sand(j,i,n) = 1.0
            end if
            if ( clay(j,i,n) < 1.0 .and. sand(j,i,n) < 100.0 ) then
              clay(j,i,n) = 1.0
            end if
          end if
        end do
      end do
    end do

    call gfclose(gfile)
  end subroutine mksoitex

end module mod_mksoitex
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

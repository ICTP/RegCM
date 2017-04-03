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
  use mod_dynparam
  use mod_grid
  use mod_rdldtr

  implicit none

  private

  public :: mksoitex

  character(len=16) , parameter :: mapname = 'MAPUNITS'
  character(len=24) , parameter :: mapdim = 'number_of_mapunits'
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
    character(len=256) :: inpfile

    nc = size(sand,3)

    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//soitexfile
    call gfopen(gfile,inpfile,xlat,xlon,ds*nsg,i_band)
    call gfread(gfile,varname1,mapdim,mapname,sand,.false.,h_missing_value)
    call gfread(gfile,varname2,mapdim,mapname,clay,.false.,h_missing_value)

    do i = 1 , iysg
      do j = 1 , jxsg
        if ( mask(j,i) > 0.5_rkx ) then
          if ( clay(j,i,1) < 0.0_rkx .or. sand(j,i,1) < 0.0_rkx ) then
            do n = 1 , nc
              call bestaround(sand(:,:,n),i,j)
              call bestaround(clay(:,:,n),i,j)
            end do
          end if
          do n = 1 , nc
            if ( clay(j,i,n) + sand(j,i,n) /= 100.0_rkx ) then
              if ( clay(j,i,n) > sand(j,i,n) ) then
                sand(j,i,n) = 100.0_rkx - clay(j,i,n)
              else
                clay(j,i,n) = 100.0_rkx - sand(j,i,n)
              end if
            end if
          end do
        else
          sand(j,i,:) = h_missing_value
          clay(j,i,:) = h_missing_value
        end if
      end do
    end do

    call gfclose(gfile)
  end subroutine mksoitex

end module mod_mksoitex
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

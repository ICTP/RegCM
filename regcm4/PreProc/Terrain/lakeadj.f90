!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      subroutine lakeadj(lsmtyp,lnduse,htgrid,xlat,xlon,imx,jmx)
 
      implicit none
!
! PARAMETER definitions
!
      real(4) , parameter :: zerie = 174. , zhuron = 177. ,             &
                           & zontar = 75. , zsup = 183. , zmich = 177.
!
! Dummy arguments
!
      integer :: imx , jmx
      character(4) :: lsmtyp
      real(4) , dimension(imx,jmx) :: htgrid , xlat , xlon
      integer , dimension(imx,jmx) :: lnduse
      intent (in) imx , jmx , lnduse , lsmtyp , xlat , xlon
      intent (inout) htgrid
!
! Local variables
!
      integer :: i , j
      real(4) :: xx , yy
!
!     ****  ADJUST GREAT LAKE ELEVATION **** C
!
      do i = 1 , imx
        do j = 1 , jmx
          if ( (lsmtyp=='BATS' .and. lnduse(i,j)==14) .or.              &
              &(lsmtyp=='USGS' .and. lnduse(i,j)==16) ) then
            xx = xlon(i,j)
            yy = xlat(i,j)
            if ( yy<=43.2 .and. yy>=41.0 .and. xx<=-78.0 .and.          &
               & xx>=-84.0 ) then                       ! LAKE ERIE
              print * , '**** ADUJUSTING LAKE ERIE LEVEL ****'
              print * , '     NEW:' , zerie , '    OLD:' , htgrid(i,j) ,&
                  & i , j
              htgrid(i,j) = zerie
            else if ( yy<=46.4 .and. yy>=43.0 .and. xx<=-79.9 .and.     &
                    & yy>=-85.0 ) then                  ! LAKE HURON
              print * , '**** ADUJUSTING LAKE HURON LEVEL ****'
              print * , '     NEW:' , zhuron , '    OLD:' , htgrid(i,j) &
                  & , i , j
              htgrid(i,j) = zhuron
            else if ( yy<=44.5 .and. yy>=43.2 .and. xx<=-75.0 .and.     &
                    & yy>=-79.9 ) then                  ! LAKE ONTARIO
              print * , '**** ADUJUSTING LAKE ONTARIO LEVEL ****'
              print * , '     NEW:' , zontar , '    OLD:' , htgrid(i,j) &
                  & , i , j
              htgrid(i,j) = zontar
            else if ( yy<=49.4 .and. yy>=46.2 .and. xx<=-84.2 .and.     &
                    & xx>=-93.0 ) then                  ! LAKE SUPERIOR
              print * , '**** ADUJUSTING LAKE SUPERIOR LEVEL ****'
              print * , '     NEW:' , zsup , '    OLD:' , htgrid(i,j) , &
                  & i , j
              htgrid(i,j) = zsup
            else if ( yy<=46.2 .and. yy>=41.0 .and. xx<=-84.8 .and.     &
                    & xx>=-89.0 ) then                  ! LAKE MICHIGAN
              print * , '**** ADUJUSTING LAKE MICHIGAN LEVEL ****'
              print * , '     NEW:' , zmich , '    OLD:' , htgrid(i,j) ,&
                  & i , j
              htgrid(i,j) = zmich
            else
            end if
          end if
        end do
      end do
 
      end subroutine lakeadj

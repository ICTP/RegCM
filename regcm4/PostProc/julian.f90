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

      subroutine julian(idate,julnc,iyr,imo,idy,ihr)
 
      implicit none
!
! Dummy arguments
!
      integer :: idate , idy , ihr , imo , iyr , julnc
      intent (out) julnc
      intent (inout) idate , idy , ihr , imo , iyr
!
! Local variables
!
      integer :: ileap , iyrm1 , j , julday
      integer , dimension(12) :: jprev , lenmon
!
      data lenmon/31 , 28 , 31 , 30 , 31 , 30 , 31 , 31 , 30 , 31 , 30 ,&
         & 31/
 
      iyr = idate/1000000
      imo = idate/10000 - iyr*100
      idy = idate/100 - iyr*10000 - imo*100
      ihr = idate - idate/100*100
      ileap = mod(iyr,4)
      if ( ileap==0 ) then
        lenmon(2) = 29
      else
        lenmon(2) = 28
      end if
 
      if ( ihr>23 ) then
        ihr = ihr - 24
        idy = idy + 1
      end if
      if ( idy>lenmon(imo) ) then
        idy = 1
        imo = imo + 1
      end if
      if ( imo>12 ) then
        imo = 1
        iyr = iyr + 1
      end if
      idate = iyr*1000000 + imo*10000 + idy*100 + ihr
 
      iyrm1 = iyr - 1
 
 
      jprev(1) = 0
      do j = 2 , 12
        jprev(j) = jprev(j-1) + lenmon(j-1)
      end do
 
      julday = idy + jprev(imo) - 1
 
      julnc = ((iyr-1900)*365+julday+int((iyrm1-1900)/4))*24 + ihr
 
      end subroutine julian

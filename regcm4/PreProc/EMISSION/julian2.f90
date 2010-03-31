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

      subroutine julian2(mdate,dt,julday,julncep,jul1900,iyr,imo,idy,   &
                       & ihr,ayr)
 
!----------------------------------------------------------------
!   subroutine JULIAN calculate the year, month, day, hour
!   from the format yyyymmddhh so
!   if your input (mdate) in the above format
!   you will get yyyy mm dd hh as a sperate variabls
!   input :
!   mdate in the formate yyyymmddhh
!   dt = 1
!   output :
!   iyr  for year, imo for month, idy for day, ihr for hours
!----------------------------------------------------------------
      implicit none
!
! Dummy arguments
!
      character(4) :: ayr
      integer :: dt , idy , ihr , imo , iyr , jul1900 , julday ,        &
               & julncep , mdate
      intent (in) dt
      intent (out) ayr , jul1900 , julncep
      intent (inout) idy , ihr , imo , iyr , julday , mdate
!
! Local variables
!
      integer :: idate , idec , ihun , ileap , ione , iths , iyrm1 , j
      integer , dimension(12) :: jprev , lenmon
!
      data lenmon/31 , 28 , 31 , 30 , 31 , 30 , 31 , 31 , 30 , 31 , 30 ,&
         & 31/
 
      idate = mdate/100
      iyr = idate/10000
      iyrm1 = iyr - 1
      imo = (idate-iyr*10000)/100
      idy = mod(idate,100)
      ihr = mdate - idate*100
 
      ileap = mod(iyr,4)
      if ( ileap==0 ) lenmon(2) = 29
 
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
      mdate = iyr*1000000 + imo*10000 + idy*100 + ihr
 
      jprev(1) = 0
      do j = 2 , 12
        jprev(j) = jprev(j-1) + lenmon(j-1)
      end do
 
      julday = idy + jprev(imo) - 1
      julncep = (julday+1)*24/dt - (24/dt-1) + ihr/dt
      jul1900 = ((iyr-1900)*365+julday+int((iyrm1-1900)/4))*24 + ihr
 
      iths = iyr/1000
      ihun = (iyr-iths*1000)/100
      idec = (iyr-iths*1000-ihun*100)/10
      ione = (iyr-iths*1000-ihun*100-idec*10)
      ayr = char(iths+48)//char(ihun+48)//char(idec+48)//char(ione+48)
 
      end subroutine julian2

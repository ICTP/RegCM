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

      subroutine julian(mdate,nyrp,nmop,wt)
      implicit none
!
! Dummy arguments
!
      integer :: mdate , nmop , nyrp
      real :: wt
      intent (in) mdate
      intent (out) nyrp , wt
      intent (inout) nmop
!
! Local variables
!
      real :: fdenom , fnumer
      integer :: idate , iday , ileap , imo , iyr , j , julday , nmo ,  &
               & nyr
      integer , dimension(12) :: jprev , julmid , lenmon , midmon
! 
      data lenmon/31 , 28 , 31 , 30 , 31 , 30 , 31 , 31 , 30 , 31 , 30 ,&
         & 31/
      data midmon/16 , 14 , 16 , 15 , 16 , 15 , 16 , 16 , 15 , 16 , 15 ,&
         & 16/
 
!     ******           INITIALIZE NMOP, NYRP

      nmop = 1
      nyrp = 0
 
      idate = mdate/100
      iyr = idate/10000
      imo = (idate-iyr*10000)/100
      iday = mod(idate,100)
 
      ileap = mod(iyr,4)
      lenmon(2) = 28
      if ( ileap==0 ) lenmon(2) = 29
 
      jprev(1) = 0
      do j = 2 , 12
        jprev(j) = jprev(j-1) + lenmon(j-1)
      end do
      do j = 1 , 12
        julmid(j) = jprev(j) + midmon(j)
      end do
      julday = iday + jprev(imo)
 
!     PRINT *, 'MDATE, IYR, IMO, IDAY, JULDAY = '
!     A       ,  MDATE, IYR, IMO, IDAY, JULDAY
 
      do nyr = 1948 , 2145  !94
        do nmo = 1 , 12
 
          if ( (nyr==iyr) .and. (julmid(nmo)>julday) ) go to 100
          if ( nyr>iyr ) go to 100
 
          nmop = nmo
          nyrp = nyr
 
        end do
      end do
 
 100  continue
      fnumer = float(julday-julmid(nmop))
      if ( fnumer<0. ) fnumer = fnumer + 365.
      fdenom = float(julmid(nmo)-julmid(nmop))
      if ( fdenom<=0. ) fdenom = fdenom + 365.
      wt = fnumer/fdenom
 
!     PRINT *, 'JULMID(NMOP), JULDAY, JULMID(NMO), WT ='
!     A       ,  JULMID(NMOP), JULDAY, JULMID(NMO), WT
 
      end subroutine julian

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
 
      subroutine for_next
      use mod_dynparam
      use mod_param1 , only : nslice
      use mod_date
! 
      open(99, file='restparam.nl')
      write (99,99001) '&restartparam'
      if ( idate1.lt.globidate2 ) then
        write (99,99001) 'ifrest  = .true. '
      else
        write (99,99001) 'ifrest  = .false.'
      end if
      write (99,99002) 'idate0  = ' , idate0
      write (99,99002) 'idate1  = ' , idate1
      write (99,99002) 'idate2  = ' , globidate2
      write (99,99002) 'nslice  = ' , nslice
      write (99,99002) '/'
      close (99)

99001 format (1x,a)
99002 format (1x,a,i10,',')

      end subroutine for_next

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
 
      subroutine solar1(xtime)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine computes the solar declination angle from the   c
!     julian date.                                                    c
!                                                                     c
!     xtime  : forecast time in minutes.                              c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use message
      use param3
      use date
      implicit none
!
! Dummy arguments
!
      real(8) :: xtime
      intent (in) xtime
!
! Local variables
!
      real(8) :: calday , dayspy , decdeg , delta , pie , theta
!
!----------------------------------------------------------------------
!
!-----obecl : obliquity = 23.5 degree.
!
!KN   obecl=23.5*degrad
!KN   sinob=dsin(obecl)
!
!-----calculate longitude of the sun from vernal equinox:
!
!KN   julian=julday+nint((xtime/60.+gmt)/24.)
!KN   if (julian .ge. 81) xlong=dpd*dble(julian-81)
!KN   if (julian .lt. 81) xlong=dpd*dble(julian+284)
!KN   xlong=xlong*degrad
!KN   arg=sinob*dsin(xlong)
!KN   declin=dasin(arg)
 
!KN   added below
!KN   to take the earth's orbital eccentricity into consideration
!
      pie = 3.141592
!     dayspy = 365.
      dayspy = 365.24
      calday = dble(julday) + (nnnnnn-nstrt0)/4. + (xtime/60.+gmt)/24.
      theta = 2.*pie*calday/dayspy
!
!     Solar declination in radians:
!
      delta = .006918 - .399912*dcos(theta) + .070257*dsin(theta)       &
            & - .006758*dcos(2.*theta) + .000907*dsin(2.*theta)         &
            & - .002697*dcos(3.*theta) + .001480*dsin(3.*theta)
!
      declin = delta
!
!KN   added above
 
      decdeg = declin/degrad
!
      write (aline, 99001) decdeg
      call say
99001 format (11x,'*** solar declination angle = ',f6.2,' degrees.')
!
      end subroutine solar1

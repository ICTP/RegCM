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
 
      subroutine solar1(xtime)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine computes the solar declination angle from the   c
!     julian date.                                                    c
!                                                                     c
!     xtime  : forecast time in minutes.                              c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use mod_message
      use mod_date , only : declin , ldatez , yeardayfrac
      use mod_constants , only : mathpi , degrad , dayspy
      implicit none
!
! Dummy arguments
!
      real(8) :: xtime
      intent (in) xtime
!
! Local variables
!
      real(8) :: calday , decdeg , delta , theta
!
!----------------------------------------------------------------------
!
      calday = yeardayfrac(ldatez)
      theta = 2.*mathpi*calday/dayspy
!
!     Solar declination in radians:
!
      delta = .006918 - .399912*dcos(theta) + .070257*dsin(theta)       &
            & - .006758*dcos(2.*theta) + .000907*dsin(2.*theta)         &
            & - .002697*dcos(3.*theta) + .001480*dsin(3.*theta)
!
      declin = delta
      decdeg = declin/degrad
!
      write (aline, 99001) calday, decdeg
      call say
99001 format (11x,'*** Day ',f6.2,' solar declination angle = ',f6.2,   &
          &   ' degrees.')
!
      end subroutine solar1

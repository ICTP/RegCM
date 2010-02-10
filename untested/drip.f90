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
 
      subroutine drip
!
!  excess leaf water is put into rain or snow;
!  leaf water is reset to its maximum value.
!
      use regcm_param
      use mod_bats , only : etrrun , ldew1d , vegt , sdrop , xrun ,     &
                    & dewmx , npts , ldoc1d , sigf , c , tm
      implicit none
!
! Local variables
!
      integer :: n , np
!
      do np = np1 , npts
        do n = 1 , nnsg
          if ( ldoc1d(n,np).gt.0.5 ) then
            if ( sigf(n,np).gt.0.001 ) then
!             ***********         xrun = leaf drip ; sdrop = snow drop
!             off foliage
              etrrun(n,np) = 0.
              xrun(n,np) = ldew1d(n,np) - dewmx*vegt(n,np)
              sdrop(n,np) = 0.
 
!             ***********         test on maximum value of dew
              if ( xrun(n,np).gt.0. ) then
                etrrun(n,np) = xrun(n,np)
                ldew1d(n,np) = dewmx*vegt(n,np)
              end if
 
!             ***********         below freezing excess leaf water
!             falls as snow
              if ( (xrun(n,np).gt.0.) .and. (tm(n,np).lt.c(67)) ) then
                etrrun(n,np) = 0.
                sdrop(n,np) = xrun(n,np)
              end if
            end if
          end if
        end do
      end do
 
      end subroutine drip

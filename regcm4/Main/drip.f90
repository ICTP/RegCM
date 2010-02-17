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
      use mod_regcm_param
      use mod_bats , only : etrrun , ldew1d , vegt , sdrop , xrun ,     &
                    & ldoc1d , sigf , tm
      use mod_constants , only : dewmx , tmelt
      implicit none
!
! Local variables
!
      integer :: n , i
!
      do i = 2 , ixm1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 ) then
            if ( sigf(n,i).gt.0.001 ) then
!             ***********         xrun = leaf drip ; sdrop = snow drop
!             off foliage
              etrrun(n,i) = 0.
              xrun(n,i) = ldew1d(n,i) - dewmx*vegt(n,i)
              sdrop(n,i) = 0.
 
!             ***********         test on maximum value of dew
              if ( xrun(n,i).gt.0. ) then
                etrrun(n,i) = xrun(n,i)
                ldew1d(n,i) = dewmx*vegt(n,i)
              end if
 
!             ***********         below freezing excess leaf water
!             falls as snow
              if ( (xrun(n,i).gt.0.) .and. (tm(n,i).lt.tmelt) ) then
                etrrun(n,i) = 0.
                sdrop(n,i) = xrun(n,i)
              end if
            end if
          end if
        end do
      end do
 
      end subroutine drip

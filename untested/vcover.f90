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
 
      subroutine vcover
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
!     provides leaf and stem area parameters;
!     depends on climate through subsoil temperatures.
 
      use regcm_param
      use bats
      implicit none
!
! Local variables
!
      integer :: n , np
!
      do np = np1 , npts
        do n = 1 , nnsg
          if ( ldoc1d(n,np).gt.0.5 ) then
            if ( sigf(n,np).gt.0.001 ) seasb(n,np)                      &
               & = dmax1(0.D0,1.-0.0016*dmax1(298.-tgb1d(n,np),0.D0)**2)
          end if
        end do
      end do
 
      do np = np1 , npts
        do n = 1 , nnsg
          if ( ldoc1d(n,np).gt.0.5 ) then
            if ( sigf(n,np).gt.0.001 ) then
              xlai(n,np) = xla(lveg(n,np))
              xlai(n,np) = xlai(n,np) + (xlai0(lveg(n,np))-xlai(n,np))  &
                         & *(1.-seasb(n,np))
              rlai(n,np) = xlai(n,np) + sai(lveg(n,np))
              xlsai(n,np) = xlai(n,np) + sai(lveg(n,np))
              vegt(n,np) = sigf(n,np)*xlsai(n,np)
            end if
          end if
        end do
      end do
!
      end subroutine vcover

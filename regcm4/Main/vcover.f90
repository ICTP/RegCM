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
 
      subroutine vcover
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
!     provides leaf and stem area parameters;
!     depends on climate through subsoil temperatures.
 
      use mod_regcm_param
      use mod_bats
      implicit none
!
! Local variables
!
      integer :: n , i
!
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 ) then
            if ( sigf(n,i).gt.0.001 ) seasb(n,i)                        &
               & = dmax1(0.D0,1.-0.0016*dmax1(298.-tgb1d(n,i),0.D0)**2)
          end if
        end do
      end do
 
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 ) then
            if ( sigf(n,i).gt.0.001 ) then
              xlai(n,i) = xla(lveg(n,i))
              xlai(n,i) = xlai(n,i) + (xlai0(lveg(n,i))-xlai(n,i))      &
                         & *(1.-seasb(n,i))
              rlai(n,i) = xlai(n,i) + sai(lveg(n,i))
              xlsai(n,i) = xlai(n,i) + sai(lveg(n,i))
              vegt(n,i) = sigf(n,i)*xlsai(n,i)
            end if
          end if
        end do
      end do
!
      end subroutine vcover

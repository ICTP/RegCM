C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C    This file is part of RegCM model.
C
C    RegCM model is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    RegCM model is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
C
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      SUBROUTINE MIXER ( kmin, surf, dz, T, dnsty, depmax, depth )
c
c Simulates convective mixing
c
      implicit none

c     Arguments
      integer kmin, depmax, depth
      real(kind=8)  T(depmax,2), dnsty(depmax)
      real(kind=8)  surf, dz

c     Local variables
      real(kind=8)  avet,avev,vol,tav,tdiff
      integer k,m,k2

      do k = kmin,depth-1
        avet=0.0
        avev=0.0

        if (dnsty(k).gt.dnsty(k+1)) then

          do m = kmin, k+1
            if (m .eq. 1) then
              vol = surf
            else
              vol = dz
            end if
            avet=avet+T(m,2)*vol
            avev=avev+vol
          end do

          Tav=avet/avev
          do k2 = kmin, k+1
            T(k2,2)=Tav
            dnsty(k2)=1000.0*(1.0-1.9549E-05*(dabs((T(k2,2)+273.15)
     +               -277.0))**1.68)
          end do
        end if

      end do ! K loop 

      tdiff = 0.0
      do k = kmin, depth
        tdiff = tdiff + t(k,2) - t(k,1)
        if (k.eq.1) tdiff = tdiff * surf
        T(k,1) = T(k,2)
      end do

      return
      end subroutine MIXER

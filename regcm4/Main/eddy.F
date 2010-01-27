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

      subroutine eddy (dt,surf,dz,kv,u2,T,dnsty,de,depth)

C Computes density, eddy diffusivity and variable time step

      implicit none

      real(kind=8)  dm
      parameter(dm=1.38889E-07)
Ccc   parameter(dm=5.148e-04,       ! value used in prev simulations

      real(kind=8)  dt, surf, dz, kv, u2
      real(kind=8)  ks, N2
      integer depth
      real(kind=8)  T(depth,2), de(depth), dnsty(depth)
      real(kind=8)  demax,rimax,ws,po,dpdz,z,rad,ri
      integer k,kmax

      demax = .5*dz**2 / dt
      demax = .99 * demax
      Rimax=0.0
      do k=1,depth
       dnsty(k)=1000.0*(1.0-1.9549E-05*(dabs((T(k,1)+273.15)
     +          -277.0))**1.68)
      end do

      if (u2.lt.0.5) u2 = 0.5

C******     compute eddy diffusion profile
C N2 Brunt-Vaisala number
C Ri gradient Richardson number
C dm molecular diffusion of water
C******     compute eddy diffusion profile

      ks = 0.745 * u2**(-1.84)
      ws = 0.0012 * u2
      Po = 1.0

      do k = 1, depth-1
        dpdz = (dnsty(k+1)-dnsty(k))/dz   ! gtb removed /2.0
        N2   = dpdz/dnsty(k)*9.8          ! gtb removed minus
        z    = surf + dble(k-1)           ! gtb: k was k-1
        rad = 1. + 40. * N2 * (kv*z*dexp(ks*z)/ws)**2
        if (rad .lt. 0) then
          rad = 0.0
        end if
        Ri    = (-1.0+dsqrt(rad))/20.0
        de(k) = dm+kv*ws*z*Po*dexp(-ks*z)/(1.0+37.0*Ri**2)
        if (de(k) .gt. demax) then
          de(k) = demax
        end if
        if (dabs(Ri) .gt. Rimax) then
          kmax  = k
          Rimax = dabs(Ri)
        end if
      end do
      de(depth) = 0.0

Csb   print *, 'kmax, Rimax = ', kmax, Rimax

      return
      end subroutine eddy

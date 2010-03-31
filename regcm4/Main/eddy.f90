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
 
      subroutine eddy(dt,surf,dz,kv,u2,t,dnsty,de,depth)
 
! Computes density, eddy diffusivity and variable time step
 
      use mod_constants , only : gti , tzero
      implicit none
!
! PARAMETER definitions
!
!cc   dm=5.148e-04       ! value used in prev simulations
      real(8) , parameter :: dm = 1.38889E-07
!
! Dummy arguments
!
      integer :: depth
      real(8) :: dt , dz , kv , surf , u2
      real(8) , dimension(depth) :: de , dnsty
      real(8) , dimension(depth,2) :: t
      intent (in) depth , dt , dz , kv , surf , t
      intent (inout) de , dnsty , u2
!
! Local variables
!
      real(8) :: demax , dpdz , ks , n2 , po , rad , ri , rimax , ws , z
      integer :: k
!
      demax = .5*dz**2/dt
      demax = .99*demax
      rimax = 0.0
      do k = 1 , depth
        dnsty(k) = 1000.0*(1.0-1.9549E-05*(dabs((t(k,1)+tzero)-277.0))  &
                 & **1.68)
      end do
 
      if ( u2.lt.0.5 ) u2 = 0.5
 
!******     compute eddy diffusion profile
!     N2 Brunt-Vaisala number
!     Ri gradient Richardson number
!     dm molecular diffusion of water
!******     compute eddy diffusion profile
 
      ks = 0.745*u2**(-1.84)
      ws = 0.0012*u2
      po = 1.0
 
      do k = 1 , depth - 1
        dpdz = (dnsty(k+1)-dnsty(k))/dz   ! gtb removed /2.0
        n2 = dpdz/dnsty(k)*gti            ! gtb removed minus
        z = surf + dble(k-1)              ! gtb: k was k-1
        rad = 1. + 40.*n2*(kv*z*dexp(ks*z)/ws)**2
        if ( rad.lt.0 ) rad = 0.0
        ri = (-1.0+dsqrt(rad))/20.0
        de(k) = dm + kv*ws*z*po*dexp(-ks*z)/(1.0+37.0*ri**2)
        if ( de(k).gt.demax ) de(k) = demax
        if ( dabs(ri).gt.rimax ) rimax = dabs(ri)
      end do
      de(depth) = 0.0
 
      end subroutine eddy

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
 
      subroutine temp(dt,surf,dz,t,sw,lnet,qe,qh,dnsty,de,eta,depmax,   &
                    & depth)
!*****************BEGIN SUBROUTINE TEMP********************
!             COMPUTES TEMPERATURE PROFILE                *
!**********************************************************
      implicit none
!
! PARAMETER definitions
!
      real(8) , parameter :: cp = 4179.98
!
! Dummy arguments
!
      integer :: depmax , depth
      real(8) :: dt , dz , eta , lnet , qe , qh , surf , sw
      real(8) , dimension(depmax) :: de , dnsty
      real(8) , dimension(depmax,2) :: t
      intent (in) de , depmax , depth , dt , dz , eta , lnet , qe , qh ,&
                & sw
      intent (inout) dnsty , surf , t
!
! Local variables
!
      real(8) :: bot , t1 , t2 , tdiff , top
      integer :: k
 
      surf = 1.0
 
!******    solve differential equations of heat transfer
      do k = 1 , depth
        t(k,2) = t(k,1)
      end do
 
      k = 1
      t1 = sw*(1.-dexp(-eta*surf))/(surf*dnsty(k)*cp) + (lnet+qe+qh)    &
         & /(surf*dnsty(k)*cp)
      t2 = -de(k)*(t(k,1)-t(k+1,1))/surf
      t(k,2) = t(k,2) + (t1+t2)*dt
 
      do k = 2 , depth - 1
        top = (surf+(k-2)*dz)
        bot = (surf+(k-1)*dz)
        t1 = sw*(dexp(-eta*top)-dexp(-eta*bot))/(dz*dnsty(k)*cp)
        t2 = (de(k-1)*(t(k-1,1)-t(k,1))-de(k)*(t(k,1)-t(k+1,1)))/dz
        t(k,2) = t(k,2) + (t1+t2)*dt
      end do
 
      k = depth
      top = (surf+(k-2)*dz)
      t1 = sw*dexp(-eta*top)/(dz*dnsty(k)*cp)
      t2 = de(k-1)*(t(depth-1,1)-t(depth,1))/dz
      t(k,2) = t(k,2) + (t1+t2)*dt
 
      tdiff = 0.
      do k = 1 , depth
        tdiff = tdiff + t(k,2) - t(k,1)
        if ( k.eq.1 ) tdiff = tdiff*surf
        t(k,1) = t(k,2)
        dnsty(k) = 1000.0*(1.0-1.9549E-05*(dabs((t(k,2)+273.15)-277.0)) &
                 & **1.68)
      end do
!sb   print *, 'TEMP: Total temp change = ', tdiff
 
      end subroutine temp

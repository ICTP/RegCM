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
 
      subroutine lake(iutlak,day,dt,ta,ua,za,q,sw,lnet,hsen,hlat,ts,    &
                    & evap,prec,hice,hsnow)
 
      use mod_lake
      use mod_constants , only : ep2
      implicit none
!
! PARAMETER definitions
!
!***  delta: Stephan-Boltzmann const
!***  Le:    latent heat of evaporation
      real(8) , parameter :: delta = 5.67E-08 , le = 2.47E06
      integer , parameter :: depmax = 400
!
! Dummy arguments
!
      real(8) :: day , dt , evap , hice , hlat , hsen , hsnow , lnet ,  &
               & prec , q , sw , ta , ts , ua , za
      integer :: iutlak
      intent (in) day , hlat , hsen , iutlak , q , ta , ua , za
      intent (out) ts
      intent (inout) evap , hice , hsnow
!
! Local variables
!
      real(8) , dimension(depmax) :: de , dnsty
      integer :: depth , freeze , ilake , j , jlake , k , kmin
      real(8) :: dz , ea , eta , hi , hs , kv , ld , lu , qe , qh ,     &
               & surf , tac , tcutoff , tk , u2 , z2 , zo
      real(8) , dimension(depmax,2) :: t
!
!***  dt:  time step in seconds
!***  surf:surface thickness
!***  dz:  vertical grid spacing in m
!***  kv:  von Karman const
!***  zo:  surface roughness length
      surf = 1.0
      dz = 1.0
      kv = 0.4
      zo = 0.001
      z2 = 2.0
 
!     interpolate winds at z1 m to 2m via log wind profile
      u2 = ua*dlog(z2/zo)/dlog(za/zo)
 
!******    depth: 1-m slices of lake depth
      read (iin) ilake , jlake , depth , freeze , hi , hice , hsnow ,   &
               & eta , (t(j,1),j=1,depth)
      do k = 1 , depth
        t(k,2) = t(k,1)
      end do
 
      tac = ta - 273.15
      tk = 273.15 + t(1,1)
      lu = -0.97*delta*tk**4
      ld = lnet - lu
      qe = hlat*le
      qh = hsen
 
!     convert mixing ratio to air vapor pressure
      ea = q*88.0/(ep2+0.378*q)
 
!     ******    Check if conditions exist for lake ice
      tcutoff = -0.001
      if ( (hice.eq.0.0) .and. (t(1,1).gt.tcutoff) ) then
 
!       ******    Calculate eddy diffusivities
        call eddy(dt,surf,dz,kv,u2,t,dnsty,de,depth)
 
!       ******    Lake temperature calc using BATS sensible and latent
!       heats
        call temp(dt,surf,dz,t,sw,lnet,qe,qh,dnsty,de,eta,depmax,depth)
 
!       ******    Convective mixer
        kmin = 1
        call mixer(kmin,surf,dz,t,dnsty,depmax,depth)
 
      else
 
        call ice(sw,ld,tac,u2,ea,hs,hi,hice,evap,t,depth,prec)
        if ( freeze.eq.0 ) t(1,1) = t(1,2)
 
      end if
 
      write (iout) ilake , jlake , depth , freeze , hi , hice , hsnow , &
                 & eta , (t(j,1),j=1,depth)
 
      write (iutlak) day , ilake , jlake , depth , evap , hi , hice ,   &
                   & hsnow , (t(j,1),j=1,depth)
 
      ts = t(1,1) + 273.15
      evap = evap/3600.          !  convert evap from mm/hr to mm/sec
      hice = hice*1000.          !  convert ice  from m to mm
      hsnow = hsnow*100.         !  convert snow from m depth to mm h20
 
      lcount = lcount + 1
      if ( lcount.eq.numpts ) then
        lcount = 0
        iin = 83 - iin
        iout = 83 - iout
        rewind (iin)
        rewind (iout)
      end if
 
      end subroutine lake

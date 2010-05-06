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
 
      subroutine tseice

      use mod_dynparam
      use mod_param1 , only : dtbat
      use mod_bats
      use mod_constants , only : ch2o , cice , csnw , tzero , stdp ,    &
                               & wlhf , wlhs , cpd
      implicit none
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
! 
!     routine provides sensible and latent fluxes
!             and snow melt over ice
!
!     fss = conductive heat flow through ice
!     hrl = latent heat through leads
!     hsl = sensible heat through leads
!     hs  = heat energy balance at surface of ice
!
!         sea-ice mask could be reset in here with "imelt" - but not
!                  done at present
!
!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
! Local variables
!
      real(8) :: bb , fact , fss , hrl , hs , hsl , qgrnd , ratsi ,     &
               & rhosw3 , rsd1 , rss , smc4 , smt , tg , tgrnd , wss ,  &
               & wtt
      integer :: n , i
!
      do i = 2 , iym1
        do n = 1 , nnsg
 
          if ( lveg(n,i).eq.14 ) exit   ! lake model handles this case
 
          if ( ldoc1d(n,i).gt.1.5 ) then
#ifdef SEAICE
            if ( sice1d(n,i).lt.1000. ) sice1d(n,i) = 1000.0
#endif
 
! ******                rhosw = density of snow relative to water
            rhosw3 = rhosw(n,i)**3
            imelt(n,i) = 0
! ******                cice = specific heat of sea-ice per unit volume
            rsd1 = cice*sice1d(n,i)/1000.0
            if ( scv1d(n,i).gt.0.0 ) then
              rss = csnw*scv1d(n,i)/1000.0
              ratsi = scv1d(n,i)/(1.4*rhosw3*sice1d(n,i))
              wtt = 1./(1.+ratsi)
              wss = (scv1d(n,i)+2.8*rhosw3*sice1d(n,i))                 &
                  & /(scv1d(n,i)+1.4*rhosw3*sice1d(n,i))
! ******                include snow heat capacity
              rsd1 = 0.5*(wss*rss+wtt*rsd1)
            end if
            tgb1d(n,i) = -2.0 + tzero
! ******                subsurface heat flux through ice
            fss = 7.E-4*(tgb1d(n,i)-tg1d(n,i))                          &
                & *ch2o*rhosw3/(scv1d(n,i)+1.4*rhosw3*sice1d(n,i))
            sice1d(n,i) = sice1d(n,i) + fss*dtbat/wlhf*1.087
 
! ******                set sea ice parameter for melting and return
            if ( sice1d(n,i).le.0.0 ) then
              imelt(n,i) = 1
              exit
            end if
! ******                assume lead ocean temp is -1.8c
! ******                flux of heat and moisture through leads
! ******                sat. mixing ratio at t=-1.8c is 3.3e-3
            qice(n,i) = 3.3E-3*stdp/p1d(n,i)
!
!  determine effective surface fluxes over ice, allowing for leads;
!  aarea(n,i) is set in subroutine drag.
!
            tlef1d(n,i) = ts1d(n,i)
            qgrnd = ((1.-aarea(n,i))*cdr(n,i)*qg1d(n,i)+aarea(n,i)      &
                  & *clead(n,i)*qice(n,i))/cdrx(n,i)
            tgrnd = ((1.-aarea(n,i))*cdr(n,i)*tg1d(n,i)+aarea(n,i)      &
                  & *clead(n,i)*(tzero-1.8))/cdrx(n,i)
            fact = -rhs1d(n,i)*cdrx(n,i)*vspda(n,i)
            delq1d(n,i) = (qs1d(n,i)-qgrnd)*gwet1d(n,i)
            delt1d(n,i) = ts1d(n,i) - tgrnd
! ******           output fluxes, averaged over leads and ice
            evpr1d(n,i) = fact*delq1d(n,i)
            sent1d(n,i) = fact*cpd*delt1d(n,i)
            hrl = rhs1d(n,i)*vspda(n,i)*clead(n,i)*(qice(n,i)-qs1d(n,i))
            hsl = rhs1d(n,i)*vspda(n,i)*clead(n,i)                      &
                & *(tzero-1.8-ts1d(n,i))*cpd
! ******           get fluxes over ice for sublimation (subrout snow)
! ******              and melt (below) calculation
            fseng(n,i) = (sent1d(n,i)-aarea(n,i)*hsl)/(1.-aarea(n,i))
            fevpg(n,i) = (evpr1d(n,i)-aarea(n,i)*hrl)/(1.-aarea(n,i))
            hs = fsw1d(i) - flw1d(i) - fseng(n,i) - wlhs*fevpg(n,i)
            bb = dtbat*(hs+fss)/rsd1
! ******           snow melt
            sm(n,i) = 0.
            if ( tg1d(n,i).ge.tzero ) sm(n,i) = (hs+fss)/wlhf
            if ( sm(n,i).le.0. ) sm(n,i) = 0.
            smc4 = sm(n,i)*dtbat
            if ( scv1d(n,i).lt.smc4 ) then
! ******                all snow removed, melt ice
              smt = (scv1d(n,i)/dtbat)
! ******                rho(h2o)/rho(ice) = 1.087
              sice1d(n,i) = sice1d(n,i) + dtbat*(smt-sm(n,i))*1.087
              sm(n,i) = smt
              tg1d(n,i) = tzero
! ******                set sea ice parameter for melting and return
              if ( sice1d(n,i).le.0.0 ) then
                imelt(n,i) = 1
                exit
              end if
            else
!  **********             snow or ice with no snow melting
              tg = tg1d(n,i) + bb
              if ( tg.ge.tzero ) tg1d(n,i) = tzero
              if ( tg.lt.tzero ) tg1d(n,i) = tg
            end if
          end if
        end do
      end do
 
      end subroutine tseice

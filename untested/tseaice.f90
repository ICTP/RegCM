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
 
      subroutine tseice

      use mod_regcm_param
      use mod_bats
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
      integer :: n , np
!
      do np = np1 , npts
        do n = 1 , nnsg
 
          if ( lveg(n,np).eq.14 ) exit   ! lake model handles this case
 
          if ( ldoc1d(n,np).gt.1.5 ) then
 
! ******                rhosw = density of snow relative to water
            rhosw3 = rhosw(n,np)**3
            imelt(n,np) = 0
! ******                cice = specific heat of sea-ice per unit volume
            rsd1 = cice*sice1d(n,np)/1000.0
            if ( scv1d(n,np).gt.0.0 ) then
              rss = csnw*scv1d(n,np)/1000.0
              ratsi = scv1d(n,np)/(1.4*rhosw3*sice1d(n,np))
              wtt = 1./(1.+ratsi)
              wss = (scv1d(n,np)+2.8*rhosw3*sice1d(n,np))               &
                  & /(scv1d(n,np)+1.4*rhosw3*sice1d(n,np))
! ******                include snow heat capacity
              rsd1 = 0.5*(wss*rss+wtt*rsd1)
            end if
            tgb1d(n,np) = -2.0 + c(67)
! ******                subsurface heat flux through ice
            fss = 7.E-4*(tgb1d(n,np)-tg1d(n,np))                        &
                & *ch2o*rhosw3/(scv1d(n,np)+1.4*rhosw3*sice1d(n,np))
            sice1d(n,np) = sice1d(n,np) + fss*c(4)/c(127)*1.087
 
! ******                set sea ice parameter for melting and return
            if ( sice1d(n,np).le.0.0 ) then
              imelt(n,np) = 1
              exit
            end if
! ******                assume lead ocean temp is -1.8c
! ******                flux of heat and moisture through leads
! ******                sat. mixing ratio at t=-1.8c is 3.3e-3
            qice(n,np) = 3.3E-3*c(81)/p1d(n,np)
!
!  determine effective surface fluxes over ice, allowing for leads;
!  aarea(n,np) is set in subroutine drag.
!
            tlef1d(n,np) = ts1d(n,np)
            qgrnd = ((1.-aarea(n,np))*cdr(n,np)*qg1d(n,np)+aarea(n,np)  &
                  & *clead(n,np)*qice(n,np))/cdrx(n,np)
            tgrnd = ((1.-aarea(n,np))*cdr(n,np)*tg1d(n,np)+aarea(n,np)  &
                  & *clead(n,np)*(c(67)-1.8))/cdrx(n,np)
            fact = -rhs1d(n,np)*cdrx(n,np)*vspda(n,np)
            delq1d(n,np) = (qs1d(n,np)-qgrnd)*gwet1d(n,np)
            delt1d(n,np) = ts1d(n,np) - tgrnd
! ******           output fluxes, averaged over leads and ice
            evpr1d(n,np) = fact*delq1d(n,np)
            sent1d(n,np) = fact*c(58)*delt1d(n,np)
            hrl = rhs1d(n,np)*vspda(n,np)*clead(n,np)                   &
                & *(qice(n,np)-qs1d(n,np))
            hsl = rhs1d(n,np)*vspda(n,np)*clead(n,np)                   &
                & *(c(67)-1.8-ts1d(n,np))*c(58)
! ******           get fluxes over ice for sublimation (subrout snow)
! ******              and melt (below) calculation
            fseng(n,np) = (sent1d(n,np)-aarea(n,np)*hsl)                &
                        & /(1.-aarea(n,np))
            fevpg(n,np) = (evpr1d(n,np)-aarea(n,np)*hrl)                &
                        & /(1.-aarea(n,np))
            hs = fsw1d(np) - flw1d(np) - fseng(n,np) - c(126)           &
               & *fevpg(n,np)
            bb = c(4)*(hs+fss)/rsd1
! ******           snow melt
            sm(n,np) = 0.
            if ( tg1d(n,np).ge.c(67) ) sm(n,np) = (hs+fss)/c(127)
            if ( sm(n,np).le.0. ) sm(n,np) = 0.
            smc4 = sm(n,np)*c(4)
            if ( scv1d(n,np).lt.smc4 ) then
! ******                all snow removed, melt ice
              smt = (scv1d(n,np)/c(4))
! ******                rho(h2o)/rho(ice) = 1.087
              sice1d(n,np) = sice1d(n,np) + c(4)*(smt-sm(n,np))*1.087
              sm(n,np) = smt
              tg1d(n,np) = c(67)
! ******                set sea ice parameter for melting and return
              if ( sice1d(n,np).le.0.0 ) then
                imelt(n,np) = 1
                exit
              end if
            else
!  **********             snow or ice with no snow melting
              tg = tg1d(n,np) + bb
              if ( tg.ge.c(67) ) tg1d(n,np) = c(67)
              if ( tg.lt.c(67) ) tg1d(n,np) = tg
            end if
          end if
        end do
      end do
 
      end subroutine tseice

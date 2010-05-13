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
 
      subroutine radctl(jslc,alat,coslat,ts,pmid,pint,pmln,piln,t,      &
                      & h2ommr,cld,effcld,clwp,albs,albsd,albl,albld,   &
                      & fsns,qrs,qrl,flwds,rel,rei,fice,sols,soll,solsd,&
                      & solld,emiss1d,fsnt,fsntc,fsnsc,flnt,flns,flntc, &
                      & flnsc,solin,alb,albc,fsds,fsnirt,fsnrtc,        &
                      & fsnirtsq,eccf,o3vmr)
!
!-----------------------------------------------------------------------
!
! Driver for radiation computation.
!
! Radiation uses cgs units, so conversions must be done from
! model fields to radiation fields.
!
! Calling sequence:
!
!     radinp      Converts units of model fields and computes ozone
!                 mixing ratio for solar scheme
!
!     radcsw      Performs solar computation
!       radalb    Computes surface albedos
!       radded    Computes delta-Eddington solution
!       radclr    Computes diagnostic clear sky fluxes
!
!     radclw      Performs longwave computation
!
!       radtpl    Computes path quantities
!       radems    Computes emissivity
!       radabs    Computes absorptivity
!
!---------------------------Code history--------------------------------
!
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          J. Kiehl, B. Briegleb, August 1992
!
! Modified:          B. Briegleb, March 1995 to add aerosol
!                    to shortwave code
!
!-----------------------------------------------------------------------
!
      use mod_dynparam
      use mod_date
      use mod_aerosol , only : aermix , aeroppt , aerout
      implicit none
!
!     Input arguments
!
! ts      - Surface (skin) temperature
! pmid    - Model level pressures
! pint    - Model interface pressures
! pmln    - Natural log of pmid
! rel     - liquid cloud particle effective radius
! rei     - ice effective drop size (microns)
! fice    - fractional ice content within cloud
! piln    - Natural log of pint
! t       - Model level temperatures
! h2ommr  - Model level specific humidity
! cld     - Fractional cloud cover
! effcld  - Effective fractional cloud cover
! clwp    - Cloud liquid water path
! alat    - current latitude(radians)
! coslat  - cosine latitude
!
!     Output solar arguments
!
! fsns    - Surface absorbed solar flux
! sols    - Downward solar rad onto surface (sw direct)
! soll    - Downward solar rad onto surface (lw direct)
! solsd   - Downward solar rad onto surface (sw diffuse)
! solld   - Downward solar rad onto surface (lw diffuse)
! qrs     - Solar heating rate
!
!     Output longwave arguments
!
! qrl     - Longwave cooling rate
! flwds   - Surface down longwave flux
!
!
! Dummy arguments
!
      real(8) :: eccf
      integer :: jslc
      real(8) , dimension(iym1) :: alb , albc , albl , albld , albs ,  &
                                  & albsd , alat , coslat , emiss1d ,   &
                                  & flns , flnsc , flnt , flntc ,       &
                                  & flwds , fsds , fsnirt , fsnirtsq ,  &
                                  & fsnrtc , fsns , fsnsc , fsnt ,      &
                                  & fsntc , solin , soll , solld ,      &
                                  & sols , solsd , ts
      real(8) , dimension(iym1,kzp1) :: cld , effcld , piln , pint
      real(8) , dimension(iym1,kz) :: clwp , fice , h2ommr , pmid ,  &
           & pmln , qrl , qrs , rei , rel , t
      real(8) , dimension(iym1,kz) :: o3vmr
      intent (out) alb , albc
      intent (inout) flns , flnsc , flnt , flntc , flwds , fsds ,       &
                   & fsnirt , fsnirtsq , fsnrtc , fsns , fsnsc , fsnt , &
                   & fsntc , solin
!
!---------------------------Local variables-----------------------------
!
! i        - index
! solin    - Solar incident flux
! fsnt     - Net column abs solar flux at model top
! fsntc    - Clear sky total column abs solar flux
! fsnsc    - Clear sky surface abs solar flux
! fsnirt   - Near-IR flux absorbed at toa
! fsnrtc   - Clear sky near-IR flux absorbed at toa
! fsnirtsq - Near-IR flux absorbed at toa >= 0.7 microns
! fsds     - Flux Shortwave Downwelling Surface
! flnt     - Net outgoing lw flux at model top
! flns     - Srf longwave cooling (up-down) flux
! flntc    - Clear sky lw flux at model top
! flnsc    - Clear sky lw flux at srf (up-down)
! pbr      - Model mid-level pressures (dynes/cm2)
! pnm      - Model interface pressures (dynes/cm2)
! o3vmr    - Ozone volume mixing ratio
! o3mmr    - Ozone mass mixing ratio
! plco2    - Prs weighted CO2 path
! plh2o    - Prs weighted H2O path
! tclrsf   - Total clear sky fraction, level to space
! eccf     - Earth/sun distance factor
! n2o      - nitrous oxide mass mixing ratio
! ch4      - methane mass mixing ratio
! cfc11    - cfc11 mass mixing ratio
! cfc12    - cfc12 mass mixing ratio
! rh       - level relative humidity (fraction)
!
! Local variables
!
      real(8) , dimension(iym1) :: aeradfo , aeradfos
      real(8) , dimension(iym1,kz) :: cfc11 , cfc12 , ch4 , n2o
      integer :: i
      real(8) , dimension(iym1,kz) :: o3mmr , pbr , rh
      real(8) , dimension(iym1,kzp1) :: plco2 , plh2o , pnm , tclrsf
!
!     Instead of interpolating the o3vmr from the time-interpolated
!     values, we pass compute o3vmr in getdat() and pass it directly
!     into radctl(). o3mmr will be computed in radinp().
!
!     Set latitude dependent radiation input
!
      call radinp(pmid,pint,h2ommr,cld,o3vmr,pbr,pnm,plco2,plh2o,tclrsf,&
                & eccf,o3mmr)
!
!     Solar radiation computation
!
      if ( dosw ) then
!
!       Specify aerosol mass mixing ratio
!
        call aermix(pnm,rh,jslc,1,iym1,iym1,kz,ntr)
 
        call aeroppt(rh,pint)
!
        call radcsw(pnm,h2ommr,o3mmr,cld,clwp,rel,rei,fice,eccf,albs,   &
                  & albsd,albl,albld,solin,qrs,fsns,fsnt,fsds,fsnsc,    &
                  & fsntc,sols,soll,solsd,solld,fsnirt,fsnrtc,fsnirtsq, &
                  & aeradfo,aeradfos)
!
        call aerout(jslc,aeradfo,aeradfos)
!
!       Convert units of shortwave fields needed by rest of model from
!       CGS to MKS
        do i = 1 , iym1
          solin(i) = solin(i)*1.E-3
          fsnt(i) = fsnt(i)*1.E-3
          fsns(i) = fsns(i)*1.E-3
          fsntc(i) = fsntc(i)*1.E-3
          fsnsc(i) = fsnsc(i)*1.E-3
 
          fsds(i) = fsds(i)*1.E-3
          fsnirt(i) = fsnirt(i)*1.E-3
          fsnrtc(i) = fsnrtc(i)*1.E-3
          fsnirtsq(i) = fsnirtsq(i)*1.E-3
        end do
!
!       Calculate/outfld albedo and clear sky albedo
!
        do i = 1 , iym1
          if ( solin(i).gt.0. ) then
            alb(i) = (solin(i)-fsnt(i))/solin(i)
          else
            alb(i) = 0.
          end if
        end do
!
        do i = 1 , iym1
          if ( solin(i).gt.0. ) then
            albc(i) = (solin(i)-fsntc(i))/solin(i)
          else
            albc(i) = 0.
          end if
        end do
      end if
!
!     Longwave radiation computation
!
      if ( dolw ) then
!
!       Specify trace gas mixing ratios
!
        call trcmix(pmid,alat,coslat,n2o,ch4,cfc11,cfc12)
!
        call radclw(jslc,ts,t,h2ommr,o3vmr,pbr,pnm,pmln,piln,plco2,     &
                  & plh2o,n2o,ch4,cfc11,cfc12,effcld,tclrsf,qrl,flns,   &
                  & flnt,flnsc,flntc,flwds,emiss1d)
!
!       Convert units of longwave fields needed by rest of model from
!       CGS to MKS
        do i = 1 , iym1
          flnt(i) = flnt(i)*1.E-3
          flns(i) = flns(i)*1.E-3
          flntc(i) = flntc(i)*1.E-3
          flnsc(i) = flnsc(i)*1.E-3
          flwds(i) = flwds(i)*1.E-3
        end do
      end if
 
      end subroutine radctl

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
 
      module mod_radiation

      use mod_constants
      use mod_dynparam
      use mod_runparams
      use mod_bats
      use mod_aerosol
      use mod_date
      use mod_message
#ifdef CLM
      use mod_clm
#endif

!     Used by this module only

      use mod_tracer
      use mod_scenarios

      private

      public :: allocate_mod_radiation , radini , radctl
      public :: absnxt , abstot , emstot

      ! absnxt  - Nearest layer absorptivities
      ! abstot  - Non-adjacent layer absorptivites
      ! emstot  - Total emissivity

      real(8) , allocatable, dimension(:,:,:,:)  :: absnxt , absnxt0
      real(8) , allocatable, dimension(:,:,:,:)  :: abstot , abstot0
      real(8) , allocatable, dimension(:,:,:) :: emstot , emstot0
      real(8) , allocatable, dimension(:,:,:,:):: xuinpl
!
!     Ozone
!
      integer , parameter :: pnoz = 100
      integer , parameter :: pozlon = 1
      real(8) :: cplol , cplos , ldoyoz , ndoyoz
      real(8) , dimension(1,pnoz) :: ozmix
      real(8) , dimension(pozlon,pnoz,1,2) :: ozmixm
      real(8) , dimension(pnoz) :: pin
      integer :: koz , nyroz
!
      real(8) , dimension(2) :: a1 , a2 , b1 , b2 , realk , st
      real(8) , dimension(4) :: c1 , c2 , c3 , c4 , c5 , c6 , c7
      real(8) :: c10 , c11 , c12 , c13 , c14 , c15 , c16 , c17 , c18 ,  &
               & c19 , c20 , c21 , c22 , c23 , c24 , c25 , c26 , c27 ,  &
               & c28 , c29 , c30 , c31 , c8 , c9 , cfa1 , fc1 , fwc1 ,  &
               & fwc2 , fwcoef
      real(8) :: co2vmr
      real(8) , dimension(3,4) :: coefa , coefc , coefe
      real(8) , dimension(4,4) :: coefb , coefd
      real(8) , dimension(6,2) :: coeff , coefi
      real(8) , dimension(2,4) :: coefg , coefh
      real(8) , dimension(3,2) :: coefj , coefk
!
!     H2O DMISSIVITY AND ABSORTIVITY CODFFICIDNTS
!
      data coefa/1.01400D+00 , 6.41695D-03 , 2.85787D-05 , 1.01320D+00 ,&
         & 6.86400D-03 , 2.96961D-05 , 1.02920D+00 , 1.01680D-02 ,      &
         & 5.30226D-05 , 1.02743D+00 , 9.85113D-03 , 5.00233D-05/
!
      data coefb/8.85675D+00 , -3.51620D-02 , 2.38653D-04 ,             &
         & -1.71439D-06 , 5.73841D+00 , -1.91919D-02 , 1.65993D-04 ,    &
         & -1.54665D-06 , 6.64034D+00 , 1.56651D-02 , -9.73357D-05 ,    &
         & 0.0 , 7.09281D+00 , 1.40056D-02 , -1.15774D-04 , 0.0/
!
      data coefc/9.90127D-01 , 1.22475D-03 , 4.90135D-06 , 9.89753D-01 ,&
         & 1.97081D-03 , 3.42046D-06 , 9.75230D-01 , 1.03341D-03 , 0.0 ,&
         & 9.77366D-01 , 8.60014D-04 , 0.0/
!
      data coefd/7.03047D-01 , -2.63501D-03 , -1.57023D-06 , 0.0 ,      &
         & 5.29269D-01 , -3.14754D-03 , 4.39595D-06 , 0.0 ,             &
         & 7.88193D-02 , 1.31290D-03 , 4.25827D-06 , -1.23982D-08 ,     &
         & 1.62744D-01 , 2.22847D-03 , 2.60102D-06 , -4.30133D-08/
!
      data coefe/3.93137D-02 , -4.34341D-05 , 3.74545D-07 ,             &
         & 3.67785D-02 , -3.10794D-05 , 2.94436D-07 , 7.42500D-02 ,     &
         & 3.97397D-05 , 0.0 , 7.52859D-02 , 4.18073D-05 , 0.0/
!
      data coeff/2.2037D-01 , 1.39719D-03 , -7.32011D-06 ,              &
         & -1.40262D-08 , 2.13638D-10 , -2.35955D-13 , 3.07431D-01 ,    &
         & 8.27225D-04 , -1.30067D-05 , 3.49847D-08 , 2.07835D-10 ,     &
         & -1.98937D-12/
!
      data coefg/9.04489D+00 , -9.56499D-03 , 1.80898D+01 ,             &
         & -1.91300D-02 , 8.72239D+00 , -9.53359D-03 , 1.74448D+01 ,    &
         & -1.90672D-02/
!
      data coefh/5.46557D+01 , -7.30387D-02 , 1.09311D+02 ,             &
         & -1.46077D-01 , 5.11479D+01 , -6.82615D-02 , 1.02296D+02 ,    &
         & -1.36523D-01/
!
      data coefi/3.31654D-01 , -2.86103D-04 , -7.87860D-06 ,            &
         & 5.88187D-08 , -1.25340D-10 , -1.37731D-12 , 3.14365D-01 ,    &
         & -1.33872D-03 , -2.15585D-06 , 6.07798D-08 , -3.45612D-10 ,   &
         & -9.34139D-15/
!
      data coefj/2.82096D-02 , 2.47836D-04 , 1.16904D-06 , 9.27379D-02 ,&
         & 8.04454D-04 , 6.88844D-06/
!
      data coefk/2.48852D-01 , 2.09667D-03 , 2.60377D-06 , 1.03594D+00 ,&
         & 6.58620D-03 , 4.04456D-06/
!
!     Narrow band data for H2O
!     200CM data for 800-1000 CM-1 and 1000-1200 CM-1.
!
      data realk/0.18967069430426D-04 , 0.70172244841851D-04/
      data st/0.31930234492350D-03 , 0.97907319939060D-03/
      data a1/0.28775403075736D-01 , 0.23236701470511D-01/
      data a2/ -0.57966222388131D-04 , -0.95105504388411D-04/
      data b1/0.29927771523756D-01 , 0.21737073577293D-01/
      data b2/ -0.86322071248593D-04 , -0.78543550629536D-04/
!
!
      contains 

      subroutine allocate_mod_radiation 
        implicit none        
#ifdef MPP1
        allocate(absnxt(iym1,kz,4,jxp))
        allocate(abstot(iym1,kzp1,kzp1,jxp))
        allocate(emstot(iym1,kzp1,jxp))
        allocate(absnxt0(iym1,kz,4,jxp))
        allocate(abstot0(iym1,kzp1,kzp1,jxp))
        allocate(emstot0(iym1,kzp1,jxp))
        allocate(xuinpl(iym1,kzp1,4,jxp))
#else
#ifdef BAND
        allocate(absnxt(iym1,kz,4,jx))
        allocate(abstot(iym1,kzp1,kzp1,jx))
        allocate(emstot(iym1,kzp1,jx))
        allocate(absnxt0(iym1,kz,4,jx))
        allocate(abstot0(iym1,kzp1,kzp1,jx))
        allocate(emstot0(iym1,kzp1,jx))
        allocate(xuinpl(iym1,kzp1,4,jx))       
#else
        allocate(absnxt(iym1,kz,4,jxm1))
        allocate(abstot(iym1,kzp1,kzp1,jxm1))
        allocate(emstot(iym1,kzp1,jxm1))
        allocate(absnxt0(iym1,kz,4,jxm1))
        allocate(abstot0(iym1,kzp1,kzp1,jxm1))
        allocate(emstot0(iym1,kzp1,jxm1))
        allocate(xuinpl(iym1,kzp1,4,jxm1))       
#endif 

#endif 
      end subroutine allocate_mod_radiation 

      subroutine radini

!-----------------------------------------------------------------------
!
! Initialize various constants for radiation scheme; note that
! the radiation scheme uses cgs units.
!
!---------------------------Code history--------------------------------
!
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          J. Kiehl, B. Briegleb, August 1992
!
!-----------------------------------------------------------------------
!
      implicit none
!
!---------------------------Local variables-----------------------------
!
! iband  - H2O band index
! v0     - Volume of a gas at stp (m**3/kmol)
! p0     - Standard pressure (pascals)
! amd    - Effective molecular weight of dry air (kg/kmol)
!
! Local variables
!
      real(8) :: amd , p0 , v0
      integer :: iband
!
!     Set general radiation consts; convert to cgs units where
!     appropriate:
!
!IPCC
!1991-1995
!     co2vmr  =  3.55e-4
!     ch40 = 0.55241 * 1.714e-6
!     n2o0 = 1.51913 * 0.311e-6
!1961-1965
!     co2vmr  =  3.10e-4
!     ch40 = 0.55241 * 1.414e-6
!     n2o0 = 1.51913 * 0.287e-6
 
!     cfc110 = 4.69548 * 0.280e-9
!     cfc120 = 4.14307 * 0.503e-9
!     co2mmr = 1.51913 * co2vmr
      if ( lyear.ge.1750 .and. lyear.le.2100 ) then
        co2vmr = cgas(2,lyear)*1.E-6
        co2mmr = co2vmr*44.0/28.9644
        ch40 = cgas(3,lyear)*1.E-9*0.55241
        n2o0 = cgas(4,lyear)*1.E-9*1.51913
        cfc110 = cgas(5,lyear)*1.E-12*4.69548
        cfc120 = cgas(6,lyear)*1.E-12*4.14307
      else
        write (aline,*) '  Simulation date:  ' , lyear
        call say
        call fatal(__FILE__,__LINE__,                                   &
               &'CONCENTRATION VALUES OUTSIDE OF DATE RANGE (1750-2100)'&
              & )
      end if
!     print*,'IN RADINI (TOP)'
!     print*,'  co2vmr= ',co2vmr
!     print*,'  co2mmr= ',co2mmr
!     print*,'  ch40  = ',ch40
!     print*,'  n2o0  = ',n2o0
!     print*,'  cfc110= ',cfc110
!     print*,'  cfc120= ',cfc120
!
!     Coefficients for h2o emissivity and absorptivity.
!
      do iband = 1 , 4
        c1(iband) = coefe(3,iband)/coefe(2,iband)
        c2(iband) = coefb(3,iband)/coefb(2,iband)
        c3(iband) = coefb(4,iband)/coefb(3,iband)
        c4(iband) = coefd(3,iband)/coefd(2,iband)
        c5(iband) = coefd(4,iband)/coefd(3,iband)
        c6(iband) = coefa(3,iband)/coefa(2,iband)
        c7(iband) = coefc(3,iband)/coefc(2,iband)
      end do
      c8 = coeff(3,1)/coeff(2,1)
      c9 = coeff(3,2)/coeff(2,2)
      c10 = coeff(4,1)/coeff(3,1)
      c11 = coeff(4,2)/coeff(3,2)
      c12 = coeff(5,1)/coeff(4,1)
      c13 = coeff(5,2)/coeff(4,2)
      c14 = coeff(6,1)/coeff(5,1)
      c15 = coeff(6,2)/coeff(5,2)
      c16 = coefj(3,1)/coefj(2,1)
      c17 = coefk(3,1)/coefk(2,1)
      c18 = coefi(3,1)/coefi(2,1)
      c19 = coefi(3,2)/coefi(2,2)
      c20 = coefi(4,1)/coefi(3,1)
      c21 = coefi(4,2)/coefi(3,2)
      c22 = coefi(5,1)/coefi(4,1)
      c23 = coefi(5,2)/coefi(4,2)
      c24 = coefi(6,1)/coefi(5,1)
      c25 = coefi(6,2)/coefi(5,2)
      c26 = coefj(3,2)/coefj(2,2)
      c27 = coefk(3,2)/coefk(2,2)
      c28 = .5
      c29 = .002053
      c30 = .1
      c31 = 3.0E-5
      cfa1 = .61
!
!     Initialize further longwave constants referring to far wing
!     correction; R&D refers to:
!
!     Ramanathan, V. and  P.Downey, 1986: A Nonisothermal
!     Emissivity and Absorptivity Formulation for Water Vapor
!     Journal of Geophysical Research, vol. 91., D8, pp 8649-8666
!
      fwcoef = .1           ! See eq(33) R&D
      fwc1 = .30            ! See eq(33) R&D
      fwc2 = 4.5            ! See eq(33) and eq(34) in R&D
      fc1 = 2.6             ! See eq(34) R&D
!
!     Initialize ozone data.
!
      v0 = 22.4136          ! Volume of a gas at stp (m**3/kmol)
      p0 = 0.1*sslp         ! Standard pressure (pascals)
      amd = 28.9644         ! Molecular weight of dry air (kg/kmol)
!
!     Constants for ozone path integrals (multiplication by 100 for unit
!     conversion to cgs from mks):
!
      cplos = v0/(amd*gti)*100.0
      cplol = v0/(amd*gti*p0)*0.5*100.0
!
      end subroutine radini
!
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
      real(8),  dimension(iym1)::  aerlwfo , aerlwfos
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
!       call aerout(jslc,aeradfo,aeradfos)
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
                  & flnt,flnsc,flntc,flwds,emiss1d,aerlwfo,aerlwfos)
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

!     aersols diagnostics
      if ( ichem==1 ) then
        call aerout(jslc,aeradfo,aeradfos,aerlwfo,aerlwfos)
      end if   

      end subroutine radctl
!
      subroutine radcsw(pint,h2ommr,o3mmr,cld,clwp,rel,rei,fice,eccf,   &
                      & asdir,asdif,aldir,aldif,solin,qrs,fsns,fsnt,    &
                      & fsds,fsnsc,fsntc,sols,soll,solsd,solld,fsnirt,  &
                      & fsnrtc,fsnirtsq,aeradfo,aeradfos)
 
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!qian 30/06/99,  csm new scheme: hygroscopic growth effect of
!qian            sulfate has been included
!qian            main changed codes: radcsw,
 
!-----------------------------------------------------------------------
!
! Solar radiation code
!
! Basic method is Delta-Eddington as described in:
!
!    Briegleb, Bruce P., 1992: Delta-Eddington
!    Approximation for Solar Radiation in the NCAR Community Climate Model,
!    Journal of Geophysical Research, Vol 97, D7, pp7603-7612).
!
! Two changes to the basic method described above are: (1) the distinction
! between liquid and ice particle clouds, and (2) the addition of an
! aerosol with sulfate radiative properties.
!
! Divides solar spectrum into 18 intervals from 0.2-5.0 micro-meters.
! solar flux fractions specified for each interval. allows for
! seasonally and diurnally varying solar input.  Includes molecular,
! cloud, aerosol, and surface scattering, along with h2o,o3,co2,o2,cloud,
! and surface absorption. Computes delta-eddington reflections and
! transmissions assuming homogeneously mixed layers. Adds the layers
! assuming scattering between layers to be isotropic, and distinguishes
! direct solar beam from scattered radiation.
!
! Longitude loops are broken into 1 or 2 sections, so that only daylight
! (i.e. coszrs > 0) computations are done.
!
! Note that an extra layer above the model top layer is added.
!
! cgs units are used.
!
! Special diagnostic calculation of the clear sky surface and total column
! absorbed flux is also done for cloud forcing diagnostics.
!
!
!---------------------------Code history--------------------------------
!
! Modified March 1995 to add aerosols
! Original version:  B. Briegleb
! Standardized:      J. Rosinski, June 1992
! Reviewed:          J. Kiehl, B. Briegleb, August 1992
! Reviewed:          J. Kiehl, April 1996
! Reviewed:          B. Briegleb, May 1996
! 19th Band Added:   W. Collins March 1997
! Merge Wong optics: T. Schneider, Mar 1998
!
!-----------------------------------------------------------------------
!
      implicit none
!
! PARAMETER definitions
!
! v_raytau_xx - Constants for new bands
! v_abo3_xx   - Constants for new bands
!
      real(8) , parameter :: v_raytau_35 = 0.155208 ,                   &
                           & v_raytau_64 = 0.0392 ,                     &
                           & v_abo3_35 = 2.4058030E+01 ,                &
                           & v_abo3_64 = 2.210E+01
!
!     Input arguments
!
! pint    - Interface pressure
! h2ommr  - Specific humidity (h2o mass mix ratio)
! o3mmr   - Ozone mass mixing ratio
! cld     - Fractional cloud cover
! clwp    - Layer liquid water path
! rel     - Liquid effective drop size (microns)
! rei     - Ice effective drop size (microns)
! fice    - Fractional ice content within cloud
! eccf    - Eccentricity factor (1./earth-sun dist ** 2)
! asdir   - 0.2-0.7 micro-meter srfc alb to direct rad
! aldir   - 0.7-5.0 micro-meter srfc alb to direct rad
! asdif   - 0.2-0.7 micro-meter srfc alb to diffuse mod_rad
! aldif   - 0.7-5.0 micro-meter srfc alb to diffuse mod_rad
!
!     Output arguments
!
! solin    - Incident solar flux
! qrs      - Solar heating rate
! fsns     - Surface absorbed solar flux
! fsnt     - Total column absorbed solar flux
! fsds     - Flux Shortwave Downwelling Surface
! fsnsc    - Clear sky surface absorbed solar flux
! fsntc    - Clear sky total column absorbed solar flx
! sols     - Direct solar rad incident on surface (< 0.7)
! soll     - Direct solar rad incident on surface (>= 0.7)
! solsd    - Diffuse mod_solar rad incident on surface (< 0.7)
! solld    - Diffuse mod_solar rad incident on surface (>= 0.7)
! fsnirt   - Near-IR flux absorbed at toa
! fsnrtc   - Clear sky near-IR flux absorbed at toa
! fsnirtsq - Near-IR flux absorbed at toa >= 0.7 microns
!
!
! Dummy arguments
!
      real(8) :: eccf
      real(8) , dimension(iym1) :: aeradfo , aeradfos , aldif , aldir ,&
                                  & asdif , asdir , fsds , fsnirt ,     &
                                  & fsnirtsq , fsnrtc , fsns , fsnsc ,  &
                                  & fsnt , fsntc , solin , soll ,       &
                                  & solld , sols , solsd
      real(8) , dimension(iym1,kzp1) :: cld , pint
      real(8) , dimension(iym1,kz) :: clwp , fice , h2ommr , o3mmr , &
           & qrs , rei , rel
      intent (in) aldif , aldir , asdif , asdir , cld , clwp , eccf ,   &
                & fice , h2ommr , o3mmr , pint , rei , rel
      intent (out) aeradfo , aeradfos , fsds , qrs
      intent (inout) fsnirt , fsnirtsq , fsnrtc , fsns , fsnsc , fsnt , &
                   & fsntc , solin , soll , solld , sols , solsd
!
!---------------------------Local variables-----------------------------
!
! ns       - Spectral loop index
! i        - Longitude loop index
! k        - Level loop index
! n        - Loop index for daylight
! nloop    - Number of daylight loops
! is       - Daytime start indices
! ie       - Daytime end indices
! indxsl   - Index for cloud particle properties
!
!     A. Slingo's data for cloud particle radiative properties (from 'A
!     GCM Parameterization for the Shortwave Properties of Water
!     Clouds' JAS vol. 46 may 1989 pp 1419-1427)
!
! abarl    - A coefficient for extinction optical depth
! bbarl    - B coefficient for extinction optical depth
! cbarl    - C coefficient for single particle scat albedo
! dbarl    - D coefficient for single particle scat albedo
! ebarl    - E coefficient for asymmetry parameter
! fbarl    - F coefficient for asymmetry parameter
! abarli   - A coefficient for current spectral interval
! bbarli   - B coefficient for current spectral interval
! cbarli   - C coefficient for current spectral interval
! dbarli   - D coefficient for current spectral interval
! ebarli   - E coefficient for current spectral interval
! fbarli   - F coefficient for current spectral interval
!
!     Caution... A. Slingo recommends no less than 4.0 micro-meters nor
!     greater than 20 micro-meters
!
!     ice water coefficients (Ebert and Curry,1992, JGR, 97, 3831-3836)
!
! abari    - a coefficient for extinction optical depth
! bbari    - b coefficient for extinction optical depth
! cbari    - c coefficient for single particle scat albedo
! dbari    - d coefficient for single particle scat albedo
! ebari    - e coefficient for asymmetry parameter
! fbari    - f coefficient for asymmetry parameter
! abarii   - A coefficient for current spectral interval
! bbarii   - B coefficient for current spectral interval
! cbarii   - C coefficient for current spectral interval
! dbarii   - D coefficient for current spectral interval
! ebarii   - E coefficient for current spectral interval
! fbarii   - F coefficient for current spectral interval
! delta    - Pressure (atmospheres) for stratos. h2o limit
! o2mmr    - O2 mass mixing ratio
!
!     Next series depends on spectral interval
!
! frcsol   - Fraction of solar flux in each spectral interval
! wavmin   - Min wavelength (micro-meters) of interval
! wavmax   - Max wavelength (micro-meters) of interval
! raytau   - Rayleigh scattering optical depth
! abh2o    - Absorption coefficiant for h2o (cm2/g)
! abo3     - Absorption coefficiant for o3  (cm2/g)
! abco2    - Absorption coefficiant for co2 (cm2/g)
! abo2     - Absorption coefficiant for o2  (cm2/g)
! ph2o     - Weight of h2o in spectral interval
! pco2     - Weight of co2 in spectral interval
! po2      - Weight of o2  in spectral interval
! nirwgt   - Weight for intervals to simulate satellite filter
! wgtint   - Weight for specific spectral interval
!
!     Diagnostic and accumulation arrays; note that sfltot, fswup, and
!     fswdn are not used in the computation,but are retained for future
!     use.
!
! solflx   - Solar flux in current interval
! sfltot   - Spectrally summed total solar flux
! totfld   - Spectrally summed flux divergence
! fswup    - Spectrally summed up flux
! fswdn    - Spectrally summed down flux
!
!     Cloud radiative property arrays
!
! tauxcl   - water cloud extinction optical depth
! tauxci   - ice cloud extinction optical depth
! wcl      - liquid cloud single scattering albedo
! gcl      - liquid cloud asymmetry parameter
! fcl      - liquid cloud forward scattered fraction
! wci      - ice cloud single scattering albedo
! gci      - ice cloud asymmetry parameter
! fci      - ice cloud forward scattered fraction
!
!     Various arrays and other constants:
!
! pflx     - Interface press, including extra layer
! zenfac   - Square root of cos solar zenith angle
! sqrco2   - Square root of the co2 mass mixg ratio
! tmp1     - Temporary constant array
! tmp2     - Temporary constant array
! pdel     - Pressure difference across layer
! path     - Mass path of layer
! xptop    - Lower interface pressure of extra layer
! ptho2    - Used to compute mass path of o2
! ptho3    - Used to compute mass path of o3
! pthco2   - Used to compute mass path of co2
! pthh2o   - Used to compute mass path of h2o
! h2ostr   - Inverse square root h2o mass mixing ratio
! wavmid   - Spectral interval middle wavelength
! trayoslp - Rayleigh optical depth/standard pressure
! tmp1l    - Temporary constant array
! tmp2l    - Temporary constant array
! tmp3l    - Temporary constant array
! tmp1i    - Temporary constant array
! tmp2i    - Temporary constant array
! tmp3i    - Temporary constant array
! rdenom   - Multiple scattering term
! psf      - Frac of solar flux in spect interval
!
!     Layer absorber amounts; note that 0 refers to the extra layer
!     added above the top model layer
!
! uh2o     - Layer absorber amount of h2o
! uo3      - Layer absorber amount of  o3
! uco2     - Layer absorber amount of co2
! uo2      - Layer absorber amount of  o2
!
!     Total column absorber amounts:
!
! uth2o    - Total column  absorber amount of h2o
! uto3     - Total column  absorber amount of  o3
! utco2    - Total column  absorber amount of co2
! uto2     - Total column  absorber amount of  o2
!
!     These arrays are defined for kz model layers; 0 refers to the
!     extra layer on top:
!
! rdir     - Layer reflectivity to direct rad
! rdif     - Layer reflectivity to diffuse mod_rad
! tdir     - Layer transmission to direct rad
! tdif     - Layer transmission to diffuse mod_rad
! explay   - Solar beam exp transmission for layer
! flxdiv   - Flux divergence for layer
!
!     These arrays are defined at model interfaces; 0 is the top of the
!     extra layer above the model top; kzp1 is the earth surface:
!
! rupdir   - Ref to dir rad for layers below
! rupdif   - Ref to dif rad for layers below
! rdndif   - Ref to dif rad for layers above
! exptdn   - Solar beam exp down transm from top
! tottrn   - Total transmission for layers above
! fluxup   - Up   flux at model interface
! fluxdn   - Down flux at model interface
!
! wkaer    - works table
! aeradfo  - spectrally integrated aerosol radiative forcing ( TOA)
!-----------------------------------------------------------------------
!
! Local variables
!
      real(8) , dimension(4) :: abari , abarl , bbari , bbarl , cbari , &
                              & cbarl , dbari , dbarl , ebari , ebarl , &
                              & fbari , fbarl
      real(8) :: abarii , abarli , bbarii , bbarli , cbarii , cbarli ,  &
               & dbarii , dbarli , delta , ebarii , ebarli , fbarii ,   &
               & fbarli , h2ostr , o2mmr , path , pdel , psf ,          &
               & pthco2 , pthh2o , ptho2 , ptho3 , xptop , rdenom ,     &
               & sqrco2 , tmp1 , tmp1i , tmp1l , tmp2 , tmp2i , tmp2l , &
               & tmp3i , tmp3l , trayoslp , wavmid , wgtint
      real(8) , dimension(nspi) :: abco2 , abh2o , abo2 , abo3 ,        &
                                   & frcsol , nirwgt , pco2 , ph2o ,    &
                                   & po2 , raytau , wavmax , wavmin
      real(8) , dimension(iym1,0:kz) :: explay , fci , fcl , flxdiv ,&
           & gci , gcl , rdif , rdir , tauxci , tauxcl , tdif , tdir ,  &
           & totfld , uco2 , uh2o , uo2 , uo3 , wci , wcl
      real(8) , dimension(iym1,0:kzp1) :: exptdn , fluxdn , fluxup ,  &
           & fswdn , fswup , pflx , rdndif , rupdif , rupdir , tottrn
      integer :: i , indxsl , k , n , nloop , ns
      integer , dimension(2) :: ie , is
      real(8) , dimension(iym1) :: sfltot , solflx , utco2 , uth2o ,   &
                                  & uto2 , uto3 , x0fsnrtc , x0fsnsc ,  &
                                  & x0fsntc , zenfac
      real(8) , dimension(iym1,0:kz,4) :: wkaer
      real(8) , dimension(iym1,4) :: zero
!
      data abarl/2.817E-02 , 2.682E-02 , 2.264E-02 , 1.281E-02/
      data bbarl/1.305 , 1.346 , 1.454 , 1.641/
      data cbarl/ - 5.62E-08 , -6.94E-06 , 4.64E-04 , 0.201/
      data dbarl/1.63E-07 , 2.35E-05 , 1.24E-03 , 7.56E-03/
      data ebarl/0.829 , 0.794 , 0.754 , 0.826/
      data fbarl/2.482E-03 , 4.226E-03 , 6.560E-03 , 4.353E-03/
 
      data abari/3.448E-03 , 3.448E-03 , 3.448E-03 , 3.448E-03/
      data bbari/2.431 , 2.431 , 2.431 , 2.431/
      data cbari/1.00E-05 , 1.10E-04 , 1.861E-02 , .46658/
      data dbari/0.0 , 1.405E-05 , 8.328E-04 , 2.05E-05/
      data ebari/0.7661 , 0.7730 , 0.794 , 0.9595/
      data fbari/5.851E-04 , 5.665E-04 , 7.267E-04 , 1.076E-04/
      data delta/1.70E-3/
      data o2mmr/.23143/
 
      data frcsol/.001488 , .001389 , .001290 , .001686 , .002877 ,     &
         & .003869 , .026336 , .360739 , .065392 , .526861 , .526861 ,  &
         & .526861 , .526861 , .526861 , .526861 , .526861 , .006239 ,  &
         & .001834 , .001834/

! weight for 0.64 - 0.7 microns  appropriate to clear skies over oceans

      data nirwgt/0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 , 0.0 ,       &
         & 0.320518 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 , 1.0 ,   &
         & 1.0 , 1.0/
 
      data wavmin/.200 , .245 , .265 , .275 , .285 , .295 , .305 ,      &
         & .350 , .640 , .700 , .701 , .701 , .701 , .701 , .702 ,      &
         & .702 , 2.630 , 4.160 , 4.160/
 
      data wavmax/.245 , .265 , .275 , .285 , .295 , .305 , .350 ,      &
         & .640 , .700 , 5.000 , 5.000 , 5.000 , 5.000 , 5.000 , 5.000 ,&
         & 5.000 , 2.860 , 4.550 , 4.550/
 
      data raytau/4.020 , 2.180 , 1.700 , 1.450 , 1.250 , 1.085 ,       &
         & 0.730 , v_raytau_35 , v_raytau_64 , 0.020 , .0001 , .0001 ,  &
         & .0001 , .0001 , .0001 , .0001 , .0001 , .0001 , .0001/
!
!     Absorption coefficients
!
      data abh2o/.000 , .000 , .000 , .000 , .000 , .000 , .000 , .000 ,&
         & .000 , .002 , .035 , .377 , 1.950 , 9.400 , 44.600 ,         &
         & 190.000 , .000 , .000 , .000/
 
      data abo3/5.370E+04 , 13.080E+04 , 9.292E+04 , 4.530E+04 ,        &
         & 1.616E+04 , 4.441E+03 , 1.775E+02 , v_abo3_35 , v_abo3_64 ,  &
         & .000 , .000 , .000 , .000 , .000 , .000 , .000 , .000 ,      &
         & .000 , .000/
 
      data abco2/.000 , .000 , .000 , .000 , .000 , .000 , .000 , .000 ,&
         & .000 , .000 , .000 , .000 , .000 , .000 , .000 , .000 ,      &
         & .094 , .196 , 1.963/
 
      data abo2/.000 , .000 , .000 , .000 , .000 , .000 , .000 , .000 , &
         & 1.11E-05 , 6.69E-05 , .000 , .000 , .000 , .000 , .000 ,     &
         & .000 , .000 , .000 , .000/
!
!     Spectral interval weights
!
      data ph2o/.000 , .000 , .000 , .000 , .000 , .000 , .000 , .000 , &
         & .000 , .505 , .210 , .120 , .070 , .048 , .029 , .018 ,      &
         & .000 , .000 , .000/
 
      data pco2/.000 , .000 , .000 , .000 , .000 , .000 , .000 , .000 , &
         & .000 , .000 , .000 , .000 , .000 , .000 , .000 , .000 ,      &
         & 1.000 , .640 , .360/
 
      data po2/.000 , .000 , .000 , .000 , .000 , .000 , .000 , .000 ,  &
         & 1.000 , 1.000 , .000 , .000 , .000 , .000 , .000 , .000 ,    &
         & .000 , .000 , .000/
!
!     Initialize output fields:
!
      fsds(:) = 0.0
      fsnirt(:) = 0.0
      fsnrtc(:) = 0.0
      fsnirtsq(:) = 0.0
      fsnt(:) = 0.0
      fsns(:) = 0.0
      solin(:) = 0.0
      fsnsc(:) = 0.0
      fsntc(:) = 0.0
      sols(:) = 0.0
      soll(:) = 0.0
      solsd(:) = 0.0
      solld(:) = 0.0
      sabveg(:) = 0.0
      solis(:) = 0.0
      solvs(:) = 0.0
      solvd(:) = 0.0
!
      aeradfo(:) = 0.0
      aeradfos(:) = 0.0
      x0fsntc(:) = 0.0
      x0fsnsc(:) = 0.0
      x0fsnrtc(:) = 0.0
!
      qrs(:,:) = 0.0
!
      do k = 1 , kz
        do i = 1 , iym1
          pdel = pint(i,k+1) - pint(i,k)
          path = pdel*rga
        end do
      end do
!
!     Compute starting, ending daytime loop indices:
!
      nloop = 0
      is = 0
      ie = 0
      is(1) = isrchfgt(iym1,coszrs,1,0.D0)
!
!     If night everywhere, return:
!
      if ( is(1).gt.iym1 ) return
      ie(1) = isrchfle(iym1-is(1),coszrs(is(1)+1),1,0.D0) + is(1) - 1
      nloop = 1
!
!     Possibly 2 daytime loops needed:
!
      if ( ie(1).ne.iym1 ) then
        is(2) = isrchfgt(iym1-ie(1),coszrs(ie(1)+1),1,0.D0) + ie(1)
        if ( is(2).le.iym1 ) then
          nloop = 2
          ie(2) = iym1
        end if
      end if
!
!     Define solar incident radiation and interface pressures:
!
      do n = 1 , nloop
        do i = is(n) , ie(n)
          solin(i) = scon*eccf*coszrs(i)
          pflx(i,0) = 0.
        end do
      end do
      do k = 1 , kzp1
        do n = 1 , nloop
          do i = is(n) , ie(n)
            pflx(i,k) = pint(i,k)
          end do
        end do
      end do
!
!     Compute optical paths:
!     CO2, use old scheme(as constant)
!
      tmp1 = 0.5/(gtigts*sslp)
!     co2mmr = co2vmr*(mmwco2/mmwair)
 
      sqrco2 = sqrt(co2mmr)
      do n = 1 , nloop
        do i = is(n) , ie(n)
          xptop = pflx(i,1)
          ptho2 = o2mmr*xptop*rga
          ptho3 = o3mmr(i,1)*xptop*rga
          pthco2 = sqrco2*(xptop*rga)
          h2ostr = sqrt(1./h2ommr(i,1))
          zenfac(i) = sqrt(coszrs(i))
          pthh2o = xptop**2*tmp1 + (xptop*rga)*(h2ostr*zenfac(i)*delta)
          uh2o(i,0) = h2ommr(i,1)*pthh2o
          uco2(i,0) = zenfac(i)*pthco2
          uo2(i,0) = zenfac(i)*ptho2
          uo3(i,0) = ptho3
        end do
      end do
!
      tmp2 = delta*rga
      do k = 1 , kz
        do n = 1 , nloop
          do i = is(n) , ie(n)
            pdel = pflx(i,k+1) - pflx(i,k)
            path = pdel*rga
            ptho2 = o2mmr*path
            ptho3 = o3mmr(i,k)*path
            pthco2 = sqrco2*path
            h2ostr = sqrt(1.0/h2ommr(i,k))
            pthh2o = (pflx(i,k+1)**2-pflx(i,k)**2)                      &
                   & *tmp1 + pdel*h2ostr*zenfac(i)*tmp2
            uh2o(i,k) = h2ommr(i,k)*pthh2o
            uco2(i,k) = zenfac(i)*pthco2
            uo2(i,k) = zenfac(i)*ptho2
            uo3(i,k) = ptho3
          end do
        end do
      end do
!
!     Compute column absorber amounts for the clear sky computation:
!
      do n = 1 , nloop
        do i = is(n) , ie(n)
          uth2o(i) = 0.0
          uto3(i) = 0.0
          utco2(i) = 0.0
          uto2(i) = 0.0
        end do
      end do
      do k = 1 , kz
        do n = 1 , nloop
          do i = is(n) , ie(n)
            uth2o(i) = uth2o(i) + uh2o(i,k)
            uto3(i) = uto3(i) + uo3(i,k)
            utco2(i) = utco2(i) + uco2(i,k)
            uto2(i) = uto2(i) + uo2(i,k)
          end do
        end do
      end do
!
!     Initialize spectrally integrated totals:
!
      do k = 0 , kz
        do i = 1 , iym1
          totfld(i,k) = 0.0
          fswup(i,k) = 0.0
          fswdn(i,k) = 0.0
        end do
      end do
      do i = 1 , iym1
        sfltot(i) = 0.0
        fswup(i,kzp1) = 0.0
        fswdn(i,kzp1) = 0.0
      end do
!
!     Set cloud properties for top (0) layer; so long as tauxcl is zero,
!     there is no cloud above top of model; the other cloud properties
!     are arbitrary:
!
      do n = 1 , nloop
        do i = is(n) , ie(n)
          tauxcl(i,0) = 0.
          wcl(i,0) = 0.999999
          gcl(i,0) = 0.85
          fcl(i,0) = 0.725
          tauxci(i,0) = 0.
          wci(i,0) = 0.999999
          gci(i,0) = 0.85
          fci(i,0) = 0.725
        end do
      end do
!
!     Begin spectral loop
!
      do ns = 1 , nspi
        wgtint = nirwgt(ns)
!
!       Set index for cloud particle properties based on the wavelength,
!       according to A. Slingo (1989) equations 1-3:
!       Use index 1 (0.25 to 0.69 micrometers) for visible
!       Use index 2 (0.69 - 1.19 micrometers) for near-infrared
!       Use index 3 (1.19 to 2.38 micrometers) for near-infrared
!       Use index 4 (2.38 to 4.00 micrometers) for near-infrared
!
!       Note that the minimum wavelength is encoded (with .001, .002,
!       .003) in order to specify the index appropriate for the
!       near-infrared cloud absorption properties
!
        indxsl = 0
        if ( wavmax(ns).le.0.7 ) then
          indxsl = 1
        else if ( wavmin(ns).eq.0.700 ) then
          indxsl = 2
        else if ( wavmin(ns).eq.0.701 ) then
          indxsl = 3
        else if ( wavmin(ns).eq.0.702 .or. wavmin(ns).gt.2.38 ) then
          indxsl = 4
        else
        end if
!
!       Set cloud extinction optical depth, single scatter albedo,
!       asymmetry parameter, and forward scattered fraction:
!
        abarli = abarl(indxsl)
        bbarli = bbarl(indxsl)
        cbarli = cbarl(indxsl)
        dbarli = dbarl(indxsl)
        ebarli = ebarl(indxsl)
        fbarli = fbarl(indxsl)
!
        abarii = abari(indxsl)
        bbarii = bbari(indxsl)
        cbarii = cbari(indxsl)
        dbarii = dbari(indxsl)
        ebarii = ebari(indxsl)
        fbarii = fbari(indxsl)
!
        do k = 1 , kz
          do n = 1 , nloop
            do i = is(n) , ie(n)
!
!             liquid
!
              tmp1l = abarli + bbarli/rel(i,k)
              tmp2l = 1. - cbarli - dbarli*rel(i,k)
              tmp3l = fbarli*rel(i,k)
!
!             ice
!
              tmp1i = abarii + bbarii/rei(i,k)
              tmp2i = 1. - cbarii - dbarii*rei(i,k)
              tmp3i = fbarii*rei(i,k)
!
!             Cloud fraction incorporated into cloud extinction optical
!             depth
!found
!             April 12 2000, Filippo found the different scheme here:
 
!scheme       1
!ccm3.6.6
!             tauxcl(i,k) = clwp(i,k)*tmp1l*(1.-fice(i,k))
!             $                     *cld(i,k)*sqrt(cld(i,k))
!             tauxci(i,k) = clwp(i,k)*tmp1i*fice(i,k)
!             $                     *cld(i,k)*sqrt(cld(i,k))
!
!scheme       2
!KN
              tauxcl(i,k) = clwp(i,k)*tmp1l*(1.-fice(i,k))*cld(i,k)     &
                          & /(1.+(1.-0.85)*(1.-cld(i,k))*clwp(i,k)      &
                          & *tmp1l*(1.-fice(i,k)))
              tauxci(i,k) = clwp(i,k)*tmp1i*fice(i,k)*cld(i,k)          &
                          & /(1.+(1.-0.78)*(1.-cld(i,k))*clwp(i,k)      &
                          & *tmp1i*fice(i,k))
 
!scheme       3
!EES          below replaced
!             tauxcl(i,k) = clwp(i,k)*tmp1l*(1.-fice(i,k))
!             $                     *cld(i,k)**0.85
!             tauxci(i,k) = clwp(i,k)*tmp1i*fice(i,k)
!             $                     *cld(i,k)**0.85
!found_
 
!
!             Do not let single scatter albedo be 1; delta-eddington
!             solution for non-conservative case:
!
!qian         30/06/99        wcl(i,k) = dmin1(tmp2l,.999999)
              wcl(i,k) = dmin1(tmp2l,.999999D0)
              gcl(i,k) = ebarli + tmp3l
              fcl(i,k) = gcl(i,k)*gcl(i,k)
!
              wci(i,k) = dmin1(tmp2i,.999999D0)
              gci(i,k) = ebarii + tmp3i
              fci(i,k) = gci(i,k)*gci(i,k)
!
            end do
          end do
        end do
!
!       Set reflectivities for surface based on mid-point wavelength
!
        wavmid = 0.5*(wavmin(ns)+wavmax(ns))
!
!       Wavelength less  than 0.7 micro-meter
!
        if ( wavmid.lt.0.7 ) then
          do n = 1 , nloop
            do i = is(n) , ie(n)
              albdir(i) = asdir(i)
              albdif(i) = asdif(i)
            end do
          end do
!
!         Wavelength greater than 0.7 micro-meter
!
        else
          do n = 1 , nloop
            do i = is(n) , ie(n)
              albdir(i) = aldir(i)
              albdif(i) = aldif(i)
            end do
          end do
        end if
        trayoslp = raytau(ns)/sslp
!
!       Layer input properties now completely specified; compute the
!       delta-Eddington solution reflectivities and transmissivities
!       for each layer, starting from the top and working downwards:
 
!       options for aerosol: no climatic feedback if idirect .eq. 1
!       should be consistent with aeroppt routine

        if (ichem==1 ) then
          if ( idirect.eq.2 ) then
            do k = 0 , kz
              do i = 1 , iym1
                wkaer(i,k,1) = tauxar_mix(i,k,ns)
                wkaer(i,k,2) = tauasc_mix(i,k,ns)
                wkaer(i,k,3) = gtota_mix(i,k,ns)
                wkaer(i,k,4) = ftota_mix(i,k,ns)
              end do
            end do
          else if ( idirect.eq.1 ) then
            do k = 0 , kz
              do i = 1 , iym1
                wkaer(i,k,1) = 0.
                wkaer(i,k,2) = 0.
                wkaer(i,k,3) = 0.
                wkaer(i,k,4) = 0.
              end do
            end do
          end if
        else
          do k = 0 , kz
            do i = 1 , iym1
              wkaer(i,k,1) = 0.
              wkaer(i,k,2) = 0.
              wkaer(i,k,3) = 0.
              wkaer(i,k,4) = 0.
            end do
          end do
        end if
 
        call radded(coszrs,trayoslp,pflx,abh2o(ns),abo3(ns),abco2(ns),  &
                  & abo2(ns),uh2o,uo3,uco2,uo2,tauxcl,wcl,gcl,fcl,      &
                  & tauxci,wci,gci,fci,wkaer(1,0,1),wkaer(1,0,2),       &
                  & wkaer(1,0,3),wkaer(1,0,4),nloop,is,ie,rdir,rdif,    &
                  & tdir,tdif,explay,exptdn,rdndif,tottrn)
!
!       Compute reflectivity to direct and diffuse mod_radiation for layers
!       below by adding succesive layers starting from the surface and
!       working upwards:
!
        do n = 1 , nloop
          do i = is(n) , ie(n)
            rupdir(i,kzp1) = albdir(i)
            rupdif(i,kzp1) = albdif(i)
          end do
        end do
        do k = kz , 0 , -1
          do n = 1 , nloop
            do i = is(n) , ie(n)
              rdenom = 1./(1.-rdif(i,k)*rupdif(i,k+1))
              rupdir(i,k) = rdir(i,k) + tdif(i,k)                       &
                          & *(rupdir(i,k+1)*explay(i,k)+rupdif(i,k+1)   &
                          & *(tdir(i,k)-explay(i,k)))*rdenom
              rupdif(i,k) = rdif(i,k) + rupdif(i,k+1)*tdif(i,k)         &
                          & **2*rdenom
            end do
          end do
        end do
!
!       Compute up and down fluxes for each interface, using the added
!       atmospheric layer properties at each interface:
!
        do k = 0 , kzp1
          do n = 1 , nloop
            do i = is(n) , ie(n)
              rdenom = 1./(1.-rdndif(i,k)*rupdif(i,k))
              fluxup(i,k) = (exptdn(i,k)*rupdir(i,k)+                   &
                          & (tottrn(i,k)-exptdn(i,k))*rupdif(i,k))*     &
                          & rdenom
              fluxdn(i,k) = exptdn(i,k)                                 &
                          & + (tottrn(i,k)-exptdn(i,k)+exptdn(i,k)      &
                          & *rupdir(i,k)*rdndif(i,k))*rdenom
            end do
          end do
        end do
!
!       Compute flux divergence in each layer using the interface up
!       and down fluxes:
!
        do k = 0 , kz
          do n = 1 , nloop
            do i = is(n) , ie(n)
              flxdiv(i,k) = (fluxdn(i,k)-fluxdn(i,k+1))                 &
                          & + (fluxup(i,k+1)-fluxup(i,k))
            end do
          end do
        end do
!
!       Monochromatic computation completed; accumulate in totals;
!       adjust fraction within spectral interval to allow for the
!       possibility of sub-divisions within a particular interval:
!
        psf = 1.0
        if ( ph2o(ns).ne.0. ) psf = psf*ph2o(ns)
        if ( pco2(ns).ne.0. ) psf = psf*pco2(ns)
        if ( po2(ns).ne.0. ) psf = psf*po2(ns)
        do n = 1 , nloop
          do i = is(n) , ie(n)
            solflx(i) = solin(i)*frcsol(ns)*psf
            fsnt(i) = fsnt(i) + solflx(i)*(fluxdn(i,1)-fluxup(i,1))
 
            fsns(i) = fsns(i) + solflx(i)                               &
                    & *(fluxdn(i,kzp1)-fluxup(i,kz + 1))
 
            sfltot(i) = sfltot(i) + solflx(i)
            fswup(i,0) = fswup(i,0) + solflx(i)*fluxup(i,0)
            fswdn(i,0) = fswdn(i,0) + solflx(i)*fluxdn(i,0)
!
!           Down spectral fluxes need to be in mks; thus the .001
!           conversion factors
            if ( wavmid.lt.0.7 ) then
              sols(i) = sols(i) + exptdn(i,kzp1)*solflx(i)*0.001
              solsd(i) = solsd(i) + (fluxdn(i,kzp1)-exptdn(i,kz + 1))   &
                       & *solflx(i)*0.001
!KN           added below
              sabveg(i) = sabveg(i)                                     &
                        & + (solflx(i)*(fluxdn(i,kzp1)-fluxup(i,kz + 1))&
                        & )*(1.-albvs(i))/(1.-albdir(i))*0.001
!KN           added above
            else
              soll(i) = soll(i) + exptdn(i,kzp1)*solflx(i)*0.001
              solld(i) = solld(i) + (fluxdn(i,kzp1)-exptdn(i,kz + 1))   &
                       & *solflx(i)*0.001
              fsnirtsq(i) = fsnirtsq(i) + solflx(i)                     &
                          & *(fluxdn(i,0)-fluxup(i,0))
!KN           added below
              sabveg(i) = sabveg(i)                                     &
                        & + (solflx(i)*(fluxdn(i,kzp1)-fluxup(i,kz + 1))&
                        & )*(1.-albvl(i))/(1.-albdir(i))*0.001
!KN           added above
            end if
            fsnirt(i) = fsnirt(i) + wgtint*solflx(i)                    &
                      & *(fluxdn(i,0)-fluxup(i,0))
 
!
          end do
        end do
        do k = 0 , kz
          do n = 1 , nloop
            do i = is(n) , ie(n)
              totfld(i,k) = totfld(i,k) + solflx(i)*flxdiv(i,k)
              fswup(i,k+1) = fswup(i,k+1) + solflx(i)*fluxup(i,k+1)
              fswdn(i,k+1) = fswdn(i,k+1) + solflx(i)*fluxdn(i,k+1)
            end do
          end do
        end do
 
!       solis is incident visible solar radiation
        if ( ns.eq.8 ) then
!         -trapuv
!         do i=1,iym1
!         solis(i)=solflx(i)*0.001*fluxdn(i,kzp1)
!         end do
!         -trapuv_
          do n = 1 , nloop
            do i = is(n) , ie(n)
              solvs(i) = exptdn(i,kzp1)*solflx(i)*0.001
              solvd(i) = (fluxdn(i,kzp1)-exptdn(i,kz + 1))*solflx(i)    &
                       & *0.001
              solis(i) = solflx(i)*0.001*fluxdn(i,kzp1)
            end do
          end do
        end if
!EES    apr 20
 
!FAB
!       CLEAR SKY CALCULATION PLUS AEROSOL
!       FORCING RAD CLR is called 2 times , one with O aerosol OP , and
!       one with actual aerosol. DIFFERENCE  in net TOA SW for the two
!       case is saved as one more variable in the rad file. The
!       outputed TOASW ( fsntc, clrst) is accounting for aerosol.
        if ( ichem==1 .and. idirect.ge.1 ) then
 
          do i = 1 , iym1
            zero(i,1) = 0.
            zero(i,2) = 0.
            zero(i,3) = 0.
            zero(i,4) = 0.
          end do
 
!         Following code is the diagnostic clear sky computation:
!
!         Compute delta-Eddington solution reflectivities and
!         transmissivities for the entire column; note, for
!         convenience, we use mod_the same reflectivity and transmissivity
!         arrays as for the full calculation above, where 0 for layer
!         quantities refers to the entire atmospheric column, and where
!         0 for interface quantities refers to top of atmos- phere,
!         while 1 refers to the surface:
          call radclr(coszrs,trayoslp,pflx,abh2o(ns),abo3(ns),abco2(ns),&
                    & abo2(ns),uth2o,uto3,utco2,uto2,zero(1,1),zero(1,2)&
                    & ,zero(1,3),zero(1,4),nloop,is,ie,rdir,rdif,tdir,  &
                    & tdif,explay,exptdn,rdndif,tottrn)
 
 
!
!         Compute reflectivity to direct and diffuse mod_radiation for
!         entire column; 0,1 on layer quantities refers to two
!         effective layers overlying surface; 0 on interface quantities
!         refers to top of column; 2 on interface quantities refers to
!         the surface:
          do n = 1 , nloop
            do i = is(n) , ie(n)
              rupdir(i,2) = albdir(i)
              rupdif(i,2) = albdif(i)
            end do
          end do
!
          do k = 1 , 0 , -1
            do n = 1 , nloop
              do i = is(n) , ie(n)
                rdenom = 1./(1.-rdif(i,k)*rupdif(i,k+1))
                rupdir(i,k) = rdir(i,k) + tdif(i,k)                     &
                            & *(rupdir(i,k+1)*explay(i,k)+rupdif(i,k+1) &
                            & *(tdir(i,k)-explay(i,k)))*rdenom
                rupdif(i,k) = rdif(i,k) + rupdif(i,k+1)*tdif(i,k)       &
                            & **2*rdenom
              end do
            end do
          end do
!
!         Compute up and down fluxes for each interface, using the added
!         atmospheric layer properties at each interface:
!
          do k = 0 , 2
            do n = 1 , nloop
              do i = is(n) , ie(n)
                rdenom = 1./(1.-rdndif(i,k)*rupdif(i,k))
                fluxup(i,k) = (exptdn(i,k)*rupdir(i,k)+(tottrn(i,k)-    &
                            & exptdn(i,k))*rupdif(i,k))*rdenom
                fluxdn(i,k) = exptdn(i,k)                               &
                            & + (tottrn(i,k)-exptdn(i,k)+exptdn(i,k)    &
                            & *rupdir(i,k)*rdndif(i,k))*rdenom
              end do
            end do
          end do
!
          do n = 1 , nloop
            do i = is(n) , ie(n)
              x0fsntc(i) = x0fsntc(i) + solflx(i)                       &
                         & *(fluxdn(i,0)-fluxup(i,0))
              x0fsnsc(i) = x0fsnsc(i) + solflx(i)                       &
                         & *(fluxdn(i,2)-fluxup(i,2))
              x0fsnrtc(i) = x0fsnrtc(i) + wgtint*solflx(i)              &
                          & *(fluxdn(i,0)-fluxup(i,0))
 
!             SAVE the ref net TOA flux ( and put back the cumul
 
!             variables to 0.)
            end do
          end do
!
!         End of clear sky calculation with O aerosol OP
!
        end if
!
!       Following code is the diagnostic clear sky computation:
!
!       Compute delta-Eddington solution reflectivities and
!       transmissivities for the entire column; note, for convenience,
!       we use mod_the same reflectivity and transmissivity arrays as for
!       the full calculation above, where 0 for layer quantities refers
!       to the entire atmospheric column, and where 0 for interface
!       quantities refers to top of atmos- phere, while 1 refers to the
!       surface:
        call radclr(coszrs,trayoslp,pflx,abh2o(ns),abo3(ns),abco2(ns),  &
                  & abo2(ns),uth2o,uto3,utco2,uto2,tauxar_mix_cs(1,ns), &
                  & tauasc_mix_cs(1,ns),gtota_mix_cs(1,ns),             &
                  & ftota_mix_cs(1,ns),nloop,is,ie,rdir,rdif,tdir,tdif, &
                  & explay,exptdn,rdndif,tottrn)
 
!
!       Compute reflectivity to direct and diffuse mod_radiation for entire
!       column; 0,1 on layer quantities refers to two effective layers
!       overlying surface; 0 on interface quantities refers to top of
!       column; 2 on interface quantities refers to the surface:
!
 
        do n = 1 , nloop
          do i = is(n) , ie(n)
 
            rupdir(i,2) = albdir(i)
            rupdif(i,2) = albdif(i)
 
          end do
        end do
!
        do k = 1 , 0 , -1
          do n = 1 , nloop
            do i = is(n) , ie(n)
              rdenom = 1./(1.-rdif(i,k)*rupdif(i,k+1))
              rupdir(i,k) = rdir(i,k) + tdif(i,k)                       &
                          & *(rupdir(i,k+1)*explay(i,k)+rupdif(i,k+1)   &
                          & *(tdir(i,k)-explay(i,k)))*rdenom
              rupdif(i,k) = rdif(i,k) + rupdif(i,k+1)*tdif(i,k)         &
                          & **2*rdenom
            end do
          end do
        end do
!
!       Compute up and down fluxes for each interface, using the added
!       atmospheric layer properties at each interface:
!
        do k = 0 , 2
          do n = 1 , nloop
            do i = is(n) , ie(n)
              rdenom = 1./(1.-rdndif(i,k)*rupdif(i,k))
              fluxup(i,k) = (exptdn(i,k)*rupdir(i,k)+(tottrn(i,k)-exptdn&
                          & (i,k))*rupdif(i,k))*rdenom
              fluxdn(i,k) = exptdn(i,k)                                 &
                          & + (tottrn(i,k)-exptdn(i,k)+exptdn(i,k)      &
                          & *rupdir(i,k)*rdndif(i,k))*rdenom
            end do
          end do
        end do
!
        do n = 1 , nloop
          do i = is(n) , ie(n)
            fsntc(i) = fsntc(i) + solflx(i)*(fluxdn(i,0)-fluxup(i,0))
            fsnsc(i) = fsnsc(i) + solflx(i)*(fluxdn(i,2)-fluxup(i,2))
            fsnrtc(i) = fsnrtc(i) + wgtint*solflx(i)                    &
                      & *(fluxdn(i,0)-fluxup(i,0))
 
          end do
        end do
!
!       End of clear sky calculation
!
      end do                    ! End of spectral interval loop
 
!     FAB calculation of TOA aerosol radiative forcing
      if ( ichem==1 .and. idirect.ge.1 ) then
        do n = 1 , nloop
          do i = is(n) , ie(n)
            aeradfo(i) = -(x0fsntc(i)-fsntc(i))
            aeradfos(i) = -(x0fsnsc(i)-fsnsc(i))
          end do
        end do
      end if
!
!     Compute solar heating rate (k/s)
!
      do k = 1 , kz
        do n = 1 , nloop
          do i = is(n) , ie(n)
            qrs(i,k) = -gocp*totfld(i,k)/(pint(i,k)-pint(i,k+1))
          end do
        end do
      end do
!
!     Set the downwelling flux at the surface
!
      do i = 1 , iym1
        fsds(i) = fswdn(i,kzp1)
      end do
!
      end subroutine radcsw
!
      subroutine radclw(jslc,ts,tnm,qnm,o3vmr,pmid,pint,pmln,piln,plco2,&
                      & plh2o,n2o,ch4,cfc11,cfc12,cld,tclrsf,qrl,flns,  &
                      & flnt,flnsc,flntc,flwds,emiss1d,aerlwfo,aerlwfos)

!-----------------------------------------------------------------------
!
! Compute longwave radiation heating rates and boundary fluxes
!
! Uses broad band absorptivity/emissivity method to compute clear sky;
! assumes randomly overlapped clouds with variable cloud emissivity to
! include effects of clouds.
!
! Computes clear sky absorptivity/emissivity at lower frequency (in
! general) than the model radiation frequency; uses previously computed
! and stored values for efficiency
!
! Note: This subroutine contains vertical indexing which proceeds
!       from bottom to top rather than the top to bottom indexing
!       used in the rest of the model.
!
!---------------------------Code history--------------------------------
!
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          J. Kiehl, B. Briegleb, August 1992
!
!-----------------------------------------------------------------------
!
      implicit none
!
!     Input arguments
!
! ts      - Ground (skin) temperature
! emiss1d - Emissivity of surface
!
!     Input arguments which are only passed to other routines
!
! tnm     - Level temperature
! qnm     - Level moisture field
! o3vmr   - ozone volume mixing ratio
! pmid    - Level pressure
! pint    - Model interface pressure
! pmln    - Ln(pmid)
! piln    - Ln(pint)
! plco2   - Path length co2
! plh2o   - Path length h2o
! n2o     - nitrous oxide mass mixing ratio
! ch4     - methane mass mixing ratio
! cfc11   - cfc11 mass mixing ratio
! cfc12   - cfc12 mass mixing ratio
!
!     Input/Output arguments
!
! cld     - Cloud cover
! tclrsf  - Clear sky fraction
!
!     Output arguments
!
! qrl     - Longwave heating rate
! flns    - Surface cooling flux
! flnt    - Net outgoing flux
! flnsc   - Clear sky surface cooing
! flntc   - Net clear sky outgoing flux
! flwds   - Down longwave flux at surface
!
! Dummy arguments
!
      integer :: jslc
      real(8) , dimension(iym1,kz) :: cfc11 , cfc12 , ch4 , n2o ,    &
           & o3vmr , pmid , pmln , qnm , qrl , tnm
      real(8) , dimension(iym1,kzp1) :: cld , piln , pint , plco2 ,   &
           & plh2o , tclrsf
      real(8) , dimension(iym1) :: emiss1d , flns , flnsc , flnt ,     &
                                  & flntc , flwds , ts
      real(8), dimension(iym1) :: aerlwfo , aerlwfos


      intent (in) cld , emiss1d
      intent (out) flns , flnsc , flnt , flntc , flwds , qrl , aerlwfo, &
                 & aerlwfos

      intent (inout) tclrsf
!
!---------------------------Local variables-----------------------------
!
! i       - Longitude index
! k       - Level index
! k1      - Level index
! k2      - Level index
! k3      - Level index
! km      - Level index
! km1     - Level index
! km2     - Level index
! km3     - Level index
! km4     - Level index
! tmp     - Temporary
! tmp1    - Temporary 1
! absbt   - Downward emission at model top
! plol    - O3 pressure wghted path length
! plos    - O3 path length
! co2em   - Layer co2 normalized planck funct. derivative
! co2eml  - Interface co2 normalized planck funct. deriv.
! delt    - Diff t**4 mid layer to top interface
! delt1   - Diff t**4 lower intrfc to mid layer
! bk1     - Absrptvty for vertical quadrature
! bk2     - Absrptvty for vertical quadrature
! ful     - Total upwards longwave flux
! fsul    - Clear sky upwards longwave flux
! fdl     - Total downwards longwave flux
! fsdl    - Clear sky downwards longwv flux
! fclb4   - Sig t**4 for cld bottom interfc
! fclt4   - Sig t**4 for cloud top interfc
! s       - Flx integral sum
! tplnka  - Planck fnctn temperature
! s2c     - H2o cont amount
! s2t     - H2o cont temperature
! w       - H2o path
! tplnke  - Planck fnctn temperature
! h2otr   - H2o trnmsn for o3 overlap
! co2t    - Prs wghted temperature path
! tint    - Interface temperature
! tint4   - Interface temperature**4
! tlayr   - Level temperature
! tlayr4  - Level temperature**4
! rtclrsf - 1./tclrsf(i,k)
! klov    - Cloud lowest level index
! khiv    - Cloud highest level index
! khivm   - khiv(i) - 1
!
!     Trace gas variables
!
! ucfc11  - CFC11 path length
! ucfc12  - CFC12 path length
! un2o0   - N2O path length
! un2o1   - N2O path length (hot band)
! uch4    - CH4 path length
! uco211  - CO2 9.4 micron band path length
! uco212  - CO2 9.4 micron band path length
! uco213  - CO2 9.4 micron band path length
! uco221  - CO2 10.4 micron band path length
! uco222  - CO2 10.4 micron band path length
! uco223  - CO2 10.4 micron band path length
! bn2o0   - pressure factor for n2o
! bn2o1   - pressure factor for n2o
! bch4    - pressure factor for ch4
! uptype  - p-type continuum path length
! abplnk1 - non-nearest layer Plack factor
! abplnk2 - nearest layer factor
!
! Local variables
!
      real(8) , dimension(14,iym1,kzp1) :: abplnk1 , abplnk2
      real(8) , dimension(iym1) :: absbt , bk1 , bk2 , delt , delt1 ,  &
                                  & tmp , tplnke
      real(8) , dimension(iym1,kzp1) :: bch4 , bn2o0 , bn2o1 , co2em ,&
           & co2t , fdl , fsdl , fsul , ful , h2otr , plol , plos ,     &
           & rtclrsf , s2c , s2t , tint , tint4 , tlayr , tlayr4 ,      &
           & tplnka , ucfc11 , ucfc12 , uch4 , uco211 , uco212 ,        &
           & uco213 , uco221 , uco222 , uco223 , un2o0 , un2o1 ,        &
           & uptype , w , fsul0 , fsdl0 , ful0 , fdl0
      real(8) , dimension(iym1,kz) :: co2eml , fclb4 , fclt4
      logical , dimension(iym1) :: done , start
      real(8) :: tmp1
      integer :: i , ii , k , k1 , k2 , k3 , khighest , km , km1 , km2 ,&
               & km3 , km4 , iym1c , rad , n , nradaer
      integer , dimension(iym1) :: indx , khiv , khivm , klov
      real(8) , dimension(iym1,kzp1,kzp1) :: s , s0
      real(8) , dimension(iym1,kzp1,kzp1) :: tone
!
      tone(:,:,:)=1. 
!
      do i = 1 , iym1
        rtclrsf(i,1) = 1.0/tclrsf(i,1)
      end do
!
      do k = 1 , kz
        do i = 1 , iym1
          fclb4(i,k) = 0.
          fclt4(i,k) = 0.
          tclrsf(i,k+1) = tclrsf(i,k)*(1.-cld(i,k+1))
          rtclrsf(i,k+1) = 1./tclrsf(i,k+1)
        end do
      end do
!
!     Calculate some temperatures needed to derive absorptivity and
!     emissivity, as well as some h2o path lengths
!
      call radtpl(tnm,ts,qnm,pint,plh2o,tplnka,s2c,s2t,w,tplnke,tint,   &
                & tint4,tlayr,tlayr4,pmln,piln)
 
!     do emissivity and absorptivity calculations
!     only if abs/ems computation
!
      if ( (jyear.eq.jyear0 .and. ktau.eq.0) .or.                       &
         & (mod(ktau+1,ifrabe).eq.0) ) then
 
!
!       Compute ozone path lengths at frequency of a/e calculation.
!
        call radoz2(o3vmr,pint,plol,plos)
!
!       Compute trace gas path lengths
!
        call trcpth(tnm,pint,cfc11,cfc12,n2o,ch4,qnm,ucfc11,ucfc12,     &
                  & un2o0,un2o1,uch4,uco211,uco212,uco213,uco221,uco222,&
                  & uco223,bn2o0,bn2o1,bch4,uptype)
!
!
!       Compute total emissivity:
!
        call radems(s2c,s2t,w,tplnke,plh2o,pint,plco2,tint,tint4,tlayr, &
                  & tlayr4,plol,plos,ucfc11,ucfc12,un2o0,un2o1,uch4,    &
                  & uco211,uco212,uco213,uco221,uco222,uco223,uptype,   &
                  & bn2o0,bn2o1,bch4,co2em,co2eml,co2t,h2otr,abplnk1,   &
                  & abplnk2,jslc)
!
!       Compute total absorptivity:
!
        call radabs(pmid,pint,co2em,co2eml,tplnka,s2c,s2t,w,h2otr,plco2,&
                  & plh2o,co2t,tint,tlayr,plol,plos,pmln,piln,ucfc11,   &
                  & ucfc12,un2o0,un2o1,uch4,uco211,uco212,uco213,uco221,&
                  & uco222,uco223,uptype,bn2o0,bn2o1,bch4,abplnk1,      &
                  & abplnk2,jslc)

         if (ichem .ne. 0 .and. idirect > 0) then
           abstot0(:,:,:,jslc) = abstot(:,:,:,jslc)
           emstot0(:,:,jslc) = emstot(:,:,jslc)
           absnxt0(:,:,:,jslc) = absnxt(:,:,:,jslc)   
        end if
      end if
!
!     Find the lowest and highest level cloud for each grid point
!     Note: Vertical indexing here proceeds from bottom to top
!
      do i = 1 , iym1
        klov(i) = 0
        done(i) = .false.
      end do
      do k = 1 , kz
        do i = 1 , iym1
          if ( .not.done(i) .and. cld(i,kzp2-k).gt.0.0 ) then
            done(i) = .true.
            klov(i) = k
          end if
        end do
      end do
      call whenne(iym1,klov,1,0,indx,iym1c)
      do i = 1 , iym1
        khiv(i) = klov(i)
        done(i) = .false.
      end do
      do k = kz , 1 , -1
        do ii = 1 , iym1c
          i = indx(ii)
          if ( .not.done(i) .and. cld(i,kzp2-k).gt.0.0 ) then
            done(i) = .true.
            khiv(i) = k
          end if
        end do
      end do
      do i = 1 , iym1
        khivm(i) = khiv(i) - 1
      end do
!
!     Note: Vertical indexing here proceeds from bottom to top
!
      do ii = 1 , iym1c
        i = indx(ii)
        do k = klov(i) , khiv(i)
          fclt4(i,kzp1-k) = stebol*tint4(i,kzp2-k)
          fclb4(i,kzp1-k) = stebol*tint4(i,kzp3-k)
        end do
      end do

!
! option to calculate LW aerosol radiative forcing

! FAB LW radiative forcing ( rad=1 : avec dust)
      if (ichem .ne. 0 .and. idirect > 0) then
        nradaer = 2
        fsul0(:,:) = 0.
        fsdl0(:,:) = 0.
        abstot(:,:,:,jslc) = abstot0(:,:,:,jslc)
        emstot(:,:,jslc) = emstot0(:,:,jslc)
        absnxt(:,:,:,jslc) = absnxt0(:,:,:,jslc)
      else
        nradaer = 1
      end if

      do rad = 1 , nradaer

        if (ichem==1 .and. idirect > 0 .and. rad==2 ) then
          abstot(:,:,:,jslc) = 1-(1-abstot0(:,:,:,jslc))*aerlwtr(:,:,:)
          emstot(:,:,jslc) = 1-(1-emstot0(:,:,jslc))*aerlwtr(:,:,1)
          do k = 1 , kz  ! aerlwtr defined on plev levels
            do n = 1 , 4
              absnxt(:,k,n,jslc) = 1-(1-absnxt0(:,k,n,jslc))*           &
                                & (aerlwtr(:,k,k+1)**xuinpl(:,k,n,jslc))
            end do
          end do
        end if
!
!     Compute sums used in integrals (all longitude points)
!
!     Definition of bk1 & bk2 depends on finite differencing.  for
!     trapezoidal rule bk1=bk2. trapezoidal rule applied for nonadjacent
!     layers only.
!
!     delt=t**4 in layer above current sigma level km.
!     delt1=t**4 in layer below current sigma level km.
!
        do i = 1 , iym1
          delt(i) = tint4(i,kz) - tlayr4(i,kzp1)
          delt1(i) = tlayr4(i,kzp1) - tint4(i,kzp1)
          s(i,kzp1,kzp1) = stebol*(delt1(i)*absnxt(i,kz,1,jslc)         &
                         & +delt(i)*absnxt(i,kz,4,jslc))
          s(i,kz,kzp1) = stebol*(delt(i)*absnxt(i,kz,2,jslc)+delt1(i)   &
                        & *absnxt(i,kz,3,jslc))
        end do
        do k = 1 , kz - 1
          do i = 1 , iym1
            bk2(i) = (abstot(i,k,kz,jslc)+abstot(i,k,kzp1,jslc))*0.5
            bk1(i) = bk2(i)
            s(i,k,kzp1) = stebol*(bk2(i)*delt(i)+bk1(i)*delt1(i))
          end do
        end do
!
!       All k, km>1
!
        do km = kz , 2 , -1
          do i = 1 , iym1
            delt(i) = tint4(i,km-1) - tlayr4(i,km)
            delt1(i) = tlayr4(i,km) - tint4(i,km)
          end do
          do k = kzp1 , 1 , -1
            if ( k.eq.km ) then
              do i = 1 , iym1
                bk2(i) = absnxt(i,km-1,4,jslc)
                bk1(i) = absnxt(i,km-1,1,jslc)
              end do
            else if ( k.eq.km-1 ) then
              do i = 1 , iym1
                bk2(i) = absnxt(i,km-1,2,jslc)
                bk1(i) = absnxt(i,km-1,3,jslc)
              end do
            else
              do i = 1 , iym1
                bk2(i) = (abstot(i,k,km-1,jslc)+abstot(i,k,km,jslc))*0.5
                bk1(i) = bk2(i)
              end do
            end if
            do i = 1 , iym1
              s(i,k,km) = s(i,k,km+1)                                   &
                      & + stebol*(bk2(i)*delt(i)+bk1(i)*delt1(i))
            end do
          end do
        end do
!
!       Computation of clear sky fluxes always set first level of fsul
!
        do i = 1 , iym1
          if ( iemiss.eq.1 ) then
            fsul(i,kzp1) = emiss1d(i)*(stebol*(ts(i)**4))
          else
            fsul(i,kzp1) = stebol*(ts(i)**4)
          end if
        end do
!
!       Downward clear sky fluxes store intermediate quantities in down
!       flux Initialize fluxes to clear sky values.
!
        do i = 1 , iym1
          tmp(i) = fsul(i,kzp1) - stebol*tint4(i,kzp1)
          fsul(i,1) = fsul(i,kzp1) - abstot(i,1,kzp1,jslc)*tmp(i)       &
                  & + s(i,1,2)
          fsdl(i,1) = stebol*(tplnke(i)**4)*emstot(i,1,jslc)
          ful(i,1) = fsul(i,1)
          fdl(i,1) = fsdl(i,1)
        end do
!
!       fsdl(i,kzp1) assumes isothermal layer
!
        do k = 2 , kz
          do i = 1 , iym1
            fsul(i,k) = fsul(i,kzp1) - abstot(i,k,kzp1,jslc)*tmp(i)     &
                    & + s(i,k,k+1)
            ful(i,k) = fsul(i,k)
            fsdl(i,k) = stebol*(tplnke(i)**4)*emstot(i,k,jslc)          &
                    & - (s(i,k,2)-s(i,k,k+1))
            fdl(i,k) = fsdl(i,k)
          end do
        end do
!
!       Store the downward emission from level 1 = total gas emission *
!       sigma t**4.  fsdl does not yet include all terms
!
        do i = 1 , iym1
          ful(i,kzp1) = fsul(i,kzp1)
          absbt(i) = stebol*(tplnke(i)**4)*emstot(i,kzp1,jslc)
          fsdl(i,kzp1) = absbt(i) - s(i,kzp1,2)
          fdl(i,kzp1) = fsdl(i,kzp1)
        end do

!     FAB radiative forcing sur fsul

        if (ichem==1 .and. idirect > 0 .and. rad==1 ) then
          fsul0(:,:) = fsul(:,:)! save fsul0 = no dust
          fsdl0(:,:) = fsdl(:,:)!
          ful0(:,:) = ful(:,:)
          fdl0(:,:) = fdl(:,:)
          s0(:,:,:) = s(:,:,:)
        end if

      end do ! end rad loop

!     FAB after this DO loop fsul account for dust LW effect
!     which is OK in case of idirect=2

      if (ichem ==1 .and. idirect > 0 ) then

        aerlwfo(:) = fsul0(:,1) - fsul(:,1)

!       surface lw net ! fsul(i,plevp) - fsdl(i,plevp)
!       aerlwfos(:)= fsdl0(:,kz)-fsdl(:,kz)
        aerlwfos(:) = (fsul0(:,kzp1)-fsdl0(:,kzp1))-                    &
                 &    (fsul(:,kzp1) - fsdl(:,kzp1))
         
!       return to no aerosol LW effect  situation if idirect ==1
        if ( idirect==1 ) then
          fsul(:,:) = fsul0(:,:)
          fsdl(:,:) = fsdl0(:,:)
          ful(:,:) = ful0(:,:)
          fdl(:,:) = fdl0(:,:)
          s(:,:,:) = s0(:,:,:)
        end if 

      end if ! end aersol rad diagnostic

!
!     Modifications for clouds
!
!     Further qualify longitude subset for computations.  Select only
!     those locations where there are clouds (total cloud fraction <=
!     1.e-3 treated as clear)
!
      call whenflt(iym1,tclrsf(1,kzp1),1,0.999D0,indx,iym1c)
!
!     Compute downflux at level 1 for cloudy sky
!
      do ii = 1 , iym1c
        i = indx(ii)
!
!       First clear sky flux plus flux from cloud at level 1
!
        fdl(i,kzp1) = fsdl(i,kzp1)*tclrsf(i,kz)                     &
                     & *rtclrsf(i,kzp1-khiv(i)) + fclb4(i,kz-1)      &
                     & *cld(i,kz)
      end do
!
!     Flux emitted by other layers
!     Note: Vertical indexing here proceeds from bottom to top
!
      khighest = khiv(intmax(iym1,khiv,1))
      do km = 3 , khighest
        km1 = kzp1 - km
        km2 = kzp2 - km
        km4 = kzp4 - km
        do ii = 1 , iym1c
          i = indx(ii)
          if ( km.le.khiv(i) ) then
            tmp1 = cld(i,km2)*tclrsf(i,kz)*rtclrsf(i,km2)
            fdl(i,kzp1) = fdl(i,kzp1) + (fclb4(i,km1)-s(i,kzp1,km4)) &
                         & *tmp1
          end if
        end do
      end do
!
!     Note: Vertical indexing here proceeds from bottom to top
!
      do k = 1 , khighest - 1
        k1 = kzp1 - k
        k2 = kzp2 - k
        k3 = kzp3 - k
        do ii = 1 , iym1c
          i = indx(ii)
          if ( k.ge.klov(i) .and. k.le.khivm(i) ) ful(i,k2) = fsul(i,k2)&
             & *(tclrsf(i,kzp1)*rtclrsf(i,k1))
        end do
        do km = 1 , k
          km1 = kzp1 - km
          km2 = kzp2 - km
          km3 = kzp3 - km
          do ii = 1 , iym1c
            i = indx(ii)
!
            if ( k.le.khivm(i) .and. km.ge.klov(i) .and. km.le.khivm(i) &
               & ) ful(i,k2) = ful(i,k2)                                &
                             & + (fclt4(i,km1)+s(i,k2,k3)-s(i,k2,km3))  &
                             & *cld(i,km2)*(tclrsf(i,km1)*rtclrsf(i,k1))
          end do
        end do             ! km=1,k
      end do               ! k=1,khighest-1
!
      do k = 1 , kzp1
        k2 = kzp2 - k
        k3 = kzp3 - k
        do i = 1 , iym1
          start(i) = .false.
        end do
        do ii = 1 , iym1c
          i = indx(ii)
          if ( k.ge.khiv(i) ) then
            start(i) = .true.
            ful(i,k2) = fsul(i,k2)*tclrsf(i,kzp1)                      &
                      & *rtclrsf(i,kzp1-khiv(i))
          end if
        end do
        do km = 1 , khighest
          km1 = kzp1 - km
          km2 = kzp2 - km
          km3 = kzp3 - km
          do ii = 1 , iym1c
            i = indx(ii)
            if ( start(i) .and. km.ge.klov(i) .and. km.le.khiv(i) )     &
               & ful(i,k2) = ful(i,k2)                                  &
                           & + (cld(i,km2)*tclrsf(i,km1)*rtclrsf(i,     &
                           & kzp1-khiv(i)))                            &
                           & *(fclt4(i,km1)+s(i,k2,k3)-s(i,k2,km3))
          end do
        end do          ! km=1,khighest
      end do            ! k=1,kzp1
!
!     Computation of the downward fluxes
!
      do k = 2 , khighest - 1
        k1 = kzp1 - k
        k2 = kzp2 - k
        k3 = kzp3 - k
        do ii = 1 , iym1c
          i = indx(ii)
          if ( k.le.khivm(i) ) fdl(i,k2) = 0.
        end do
        do km = k + 1 , khighest
          km1 = kzp1 - km
          km2 = kzp2 - km
          km4 = kzp4 - km
          do ii = 1 , iym1c
            i = indx(ii)
!
            if ( k.le.khiv(i) .and. km.ge.max0(k+1,klov(i)) .and.       &
               & km.le.khiv(i) ) fdl(i,k2) = fdl(i,k2)                  &
               & + (cld(i,km2)*tclrsf(i,k1)*rtclrsf(i,km2))             &
               & *(fclb4(i,km1)-s(i,k2,km4)+s(i,k2,k3))
          end do
        end do             ! km=k+1,khighest
        do ii = 1 , iym1c
          i = indx(ii)
          if ( k.le.khivm(i) ) fdl(i,k2) = fdl(i,k2) + fsdl(i,k2)       &
             & *(tclrsf(i,k1)*rtclrsf(i,kzp1-khiv(i)))
        end do
      end do               ! k=1,khighest-1
!
!     End cloud modification loops
!
!     All longitudes: store history tape quantities
!
      do i = 1 , iym1
!
!       Downward longwave flux
!
        flwds(i) = fdl(i,kzp1)
!
!       Net flux
!
        flns(i) = ful(i,kzp1) - fdl(i,kzp1)
!
!       Clear sky flux at top of atmosphere
!
        flntc(i) = fsul(i,1)
        flnsc(i) = fsul(i,kzp1) - fsdl(i,kzp1)
!
!       Outgoing ir
!
        flnt(i) = ful(i,1) - fdl(i,1)
      end do
!
!     Computation of longwave heating (k per sec)
!
      do k = 1 , kz
        do i = 1 , iym1
          qrl(i,k) = (ful(i,k)-fdl(i,k)-ful(i,k+1)+fdl(i,k+1))          &
                   & *gocp/((pint(i,k)-pint(i,k+1)))
        end do
      end do
!
      end subroutine radclw
!
      subroutine radclr(coszrs,trayoslp,pflx,abh2o,abo3,abco2,abo2,     &
                      & uth2o,uto3,utco2,uto2,tauxar_mix_css,           &
                      & tauasc_mix_css,gtota_mix_css,ftota_mix_css,     &
                      & nloop,is,ie,rdir,rdif,tdir,tdif,explay,exptdn,  &
                      & rdndif,tottrn)
!
!-----------------------------------------------------------------------
!
! Delta-Eddington solution for special clear sky computation
!
! Computes total reflectivities and transmissivities for two atmospheric
! layers: an overlying purely ozone absorbing layer, and the rest of the
! column below.
!
! For more details , see Briegleb, Bruce P., 1992: Delta-Eddington
! Approximation for Solar Radiation in the NCAR Community Climate Model,
! Journal of Geophysical Research, Vol 97, D7, pp7603-7612).
!
!---------------------------Code history--------------------------------
!
! Original version:  B. Briegleb
! Standardized:      J. Rosinski, June 1992
! Reviewed:          J. Kiehl, B. Briegleb, August 1992
!
!-----------------------------------------------------------------------
!
      implicit none
!
!     Minimum total transmission below which no layer computation are
!     done:
!
! trmin   - Minimum total transmission allowed
! wray    - Rayleigh single scatter albedo
! gray    - Rayleigh asymetry parameter
! fray    - Rayleigh forward scattered fraction
!
!
! PARAMETER definitions
!
      real(8) , parameter :: trmin = 1.E-3 , wray = 0.999999 ,          &
                           & gray = 0.0 , fray = 0.1
!
!------------------------------Arguments--------------------------------
!
!     Input arguments
!
! coszrs   - Cosine zenith angle
! trayoslp - Tray/sslp
! pflx     - Interface pressure
! abh2o    - Absorption coefficiant for h2o
! abo3     - Absorption coefficiant for o3
! abco2    - Absorption coefficiant for co2
! abo2     - Absorption coefficiant for o2
! uth2o    - Total column absorber amount of h2o
! uto3     - Total column absorber amount of  o3
! utco2    - Total column absorber amount of co2
! uto2     - Total column absorber amount of  o2
! nloop    - Number of loops (1 or 2)
! is       - Starting index for 1 or 2 loops
! ie       - Ending index for 1 or 2 loops
!
!     Input/Output arguments
!
!     Following variables are defined for each layer; note, we use
!     layer 0 to refer to the entire atmospheric column:
!
! rdir     - Layer reflectivity to direct rad
! rdif     - Layer refflectivity to diffuse mod_rad
! tdir     - Layer transmission to direct rad
! tdif     - Layer transmission to diffuse mod_rad
! explay   - Solar beam exp transmn for layer
!
!     Note that the following variables are defined on interfaces, with
!     the index k referring to the top interface of the kth layer:
!     exptdn,rdndif,tottrn; for example, tottrn(k=5) refers to the total
!     transmission to the top interface of the 5th layer.
!
! exptdn   - Solar beam exp down transmn from top
! rdndif   - Added dif ref for layers above
! tottrn   - Total transmission for layers above
!
!
! Dummy arguments
!
      real(8) :: abco2 , abh2o , abo2 , abo3 , trayoslp
      integer :: nloop
      real(8) , dimension(iym1) :: coszrs , ftota_mix_css ,            &
                                  & gtota_mix_css , tauasc_mix_css ,    &
                                  & tauxar_mix_css , utco2 , uth2o ,    &
                                  & uto2 , uto3
      real(8) , dimension(iym1,0:kz) :: explay , rdif , rdir , tdif ,&
           & tdir
      real(8) , dimension(iym1,0:kzp1) :: exptdn , pflx , rdndif ,    &
           & tottrn
      integer , dimension(2) :: ie , is
      intent (in) abco2 , abh2o , abo2 , abo3 , coszrs , ftota_mix_css ,&
                & gtota_mix_css , ie , is , nloop , pflx ,              &
                & tauasc_mix_css , tauxar_mix_css , trayoslp , utco2 ,  &
                & uth2o , uto2 , uto3
      intent (inout) explay , exptdn , rdif , rdir , rdndif , tdif ,    &
                   & tdir , tottrn
!
!---------------------------Local variables-----------------------------
!
! i        - Longitude index
! k        - Level index
! nn       - Index of longitude loops (max=nloop)
! ii       - Longitude index
! nval     - Number of long values satisfying criteria
! indx     - Array of longitude indices
! taugab   - Total column gas absorption optical depth
! tauray   - Column rayleigh optical depth
! tautot   - Total column optical depth
! wtot     - Total column single scatter albedo
! gtot     - Total column asymmetry parameter
! ftot     - Total column forward scatter fraction
! ts       - Column scaled extinction optical depth
! ws       - Column scaled single scattering albedo
! gs       - Column scaled asymmetry parameter
! rdenom   - Mulitiple scattering term
! rdirexp  - Layer direct ref times exp transmission
! tdnmexp  - Total transmission minus exp transmission
!
!---------------------------Statement functions-------------------------
!
!     Statement functions for delta-Eddington solution; for detailed
!     explanation of individual terms, see the routine 'radded'.
!
!     Intermediate terms for delta-Eddington solution
!
!
! Local variables
!
      real(8) :: alp , amg , apg , arg , e , et , extins , f , ftot ,   &
               & g , gam , gs , gtot , lm , ne , rdenom , rdirexp , t , &
               & tautot , tdnmexp , ts , ue , uu , w , ws , wtot
      real(8) :: alpha , asys , el , xgamma , n , omgs , taus , u
      integer :: i , ii , k , nn , nval
      integer , dimension(iym1) :: indx
      real(8) , dimension(iym1) :: taugab , tauray
!
      alpha(w,uu,g,e) = .75*w*uu*((1.+g*(1-w))/(1.-e*e*uu*uu))
      xgamma(w,uu,g,e) = .50*w*((3.*g*(1.-w)*uu*uu+1.)/(1.-e*e*uu*uu))
      el(w,g) = dsqrt(3.*(1-w)*(1.-w*g))
      taus(w,f,t) = (1.-w*f)*t
      omgs(w,f) = (1.-f)*w/(1.-w*f)
      asys(g,f) = (g-f)/(1.-f)
      u(w,g,e) = 1.5*(1.-w*g)/e
      n(uu,et) = ((uu+1.)*(uu+1.)/et) - ((uu-1.)*(uu-1.)*et)
!
!-----------------------------------------------------------------------
!
!     Initialize all total transmimission values to 0, so that nighttime
!     values from previous computations are not used:
!
      tottrn = 0.0
!     print*,'dans radclr', maxval(tottrn)
      do k = 0 , kzp1
        do i = 1 , iym1
          tottrn(i,k) = 0.
        end do
      end do
!
!     Compute total direct beam transmission, total transmission, and
!     reflectivity for diffuse mod_radiation (from below) for all layers
!     above each interface by starting from the top and adding layers
!     down:
!
!     The top layer is assumed to be a purely absorbing ozone layer, and
!     that the mean diffusivity for diffuse mod_transmission is 1.66:
!
      do nn = 1 , nloop
        do i = is(nn) , ie(nn)
!
          taugab(i) = abo3*uto3(i)
!
!         Limit argument of exponential to 25, in case coszrs is very
!         small:
          arg = dmin1(taugab(i)/coszrs(i),25.D0)
          explay(i,0) = dexp(-arg)
          tdir(i,0) = explay(i,0)
!
!         Same limit for diffuse mod_transmission:
!
          arg = dmin1(1.66*taugab(i),25.D0)
          tdif(i,0) = dexp(-arg)
!
          rdir(i,0) = 0.0
          rdif(i,0) = 0.0
!
!         Initialize top interface of extra layer:
!
          exptdn(i,0) = 1.0
          rdndif(i,0) = 0.0
          tottrn(i,0) = 1.0
!
          rdndif(i,1) = rdif(i,0)
          tottrn(i,1) = tdir(i,0)
!
        end do
      end do
!
!     Now, complete the rest of the column; if the total transmission
!     through the top ozone layer is less than trmin, then no
!     delta-Eddington computation for the underlying column is done:
!
      do k = 1 , 1
!
!       Initialize current layer properties to zero;only if total
!       transmission to the top interface of the current layer exceeds
!       the minimum, will these values be computed below:
!
        do nn = 1 , nloop
          do i = is(nn) , ie(nn)
!
            rdir(i,k) = 0.0
            rdif(i,k) = 0.0
            tdir(i,k) = 0.0
            tdif(i,k) = 0.0
            explay(i,k) = 0.0
!
!           Calculates the solar beam transmission, total transmission,
!           and reflectivity for diffuse mod_radiation from below at the
!           top of the current layer:
!
            exptdn(i,k) = exptdn(i,k-1)*explay(i,k-1)
            rdenom = 1./(1.-rdif(i,k-1)*rdndif(i,k-1))
            rdirexp = rdir(i,k-1)*exptdn(i,k-1)
            tdnmexp = tottrn(i,k-1) - exptdn(i,k-1)
            tottrn(i,k) = exptdn(i,k-1)*tdir(i,k-1) + tdif(i,k-1)       &
                        & *(tdnmexp+rdndif(i,k-1)*rdirexp)*rdenom
            rdndif(i,k) = rdif(i,k-1) + (rdndif(i,k-1)*tdif(i,k-1))     &
                        & *(tdif(i,k-1)*rdenom)
!
          end do
        end do
!
!       Compute next layer delta-Eddington solution only if total
!       transmission of radiation to the interface just above the layer
!       exceeds trmin.
        call whenfgt(iym1,tottrn(1,k),1,trmin,indx,nval)
        if ( nval.gt.0 ) then
!CDIR$    IVDEP
          do ii = 1 , nval
            i = indx(ii)
!
!           Remember, no ozone absorption in this layer:
!
            tauray(i) = trayoslp*pflx(i,kzp1)
            taugab(i) = abh2o*uth2o(i) + abco2*utco2(i) + abo2*uto2(i)
!
            tautot = tauray(i) + taugab(i) + tauxar_mix_css(i)
!
            wtot = (wray*tauray(i)+tauasc_mix_css(i))/tautot
!
            gtot = (gray*wray*tauray(i)+gtota_mix_css(i))/(wtot*tautot)
!
            ftot = (fray*wray*tauray(i)+ftota_mix_css(i))/(wtot*tautot)
!
            ts = taus(wtot,ftot,tautot)
            ws = omgs(wtot,ftot)
            gs = asys(gtot,ftot)
            lm = el(ws,gs)
            alp = alpha(ws,coszrs(i),gs,lm)
            gam = xgamma(ws,coszrs(i),gs,lm)
            ue = u(ws,gs,lm)
!
!           Limit argument of exponential to 25, in case lm very large:
!
            arg = dmin1(lm*ts,25.D0)
            extins = dexp(-arg)
            ne = n(ue,extins)
!
            rdif(i,k) = (ue+1.)*(ue-1.)*(1./extins-extins)/ne
            tdif(i,k) = 4.*ue/ne
!
!           Limit argument of exponential to 25, in case coszrs is very
!           small:
            arg = dmin1(ts/coszrs(i),25.D0)
            explay(i,k) = dexp(-arg)
!
            apg = alp + gam
            amg = alp - gam
            rdir(i,k) = amg*(tdif(i,k)*explay(i,k)-1.) + apg*rdif(i,k)
            tdir(i,k) = apg*tdif(i,k) + (amg*rdif(i,k)-(apg-1.))        &
                      & *explay(i,k)
!
!           Under rare conditions, reflectivies and transmissivities
!           can be negative; zero out any negative values
!
            rdir(i,k) = dmax1(rdir(i,k),0.D0)
            tdir(i,k) = dmax1(tdir(i,k),0.D0)
            rdif(i,k) = dmax1(rdif(i,k),0.D0)
            tdif(i,k) = dmax1(tdif(i,k),0.D0)
!
          end do
        end if
!
      end do
!
!     Compute total direct beam transmission, total transmission, and
!     reflectivity for diffuse mod_radiation (from below) for both layers
!     above the surface:
!
      k = 2
      do nn = 1 , nloop
        do i = is(nn) , ie(nn)
          exptdn(i,k) = exptdn(i,k-1)*explay(i,k-1)
          rdenom = 1./(1.-rdif(i,k-1)*rdndif(i,k-1))
          rdirexp = rdir(i,k-1)*exptdn(i,k-1)
          tdnmexp = tottrn(i,k-1) - exptdn(i,k-1)
          tottrn(i,k) = exptdn(i,k-1)*tdir(i,k-1) + tdif(i,k-1)         &
                      & *(tdnmexp+rdndif(i,k-1)*rdirexp)*rdenom
          rdndif(i,k) = rdif(i,k-1) + (rdndif(i,k-1)*tdif(i,k-1))       &
                      & *(tdif(i,k-1)*rdenom)
        end do
      end do
!
      end subroutine radclr
!
      subroutine radded(coszrs,trayoslp,pflx,abh2o,abo3,abco2,abo2,uh2o,&
                      & uo3,uco2,uo2,tauxcl,wcl,gcl,fcl,tauxci,wci,gci, &
                      & fci,tauxar_mixs,tauasc_mixs,gtota_mixs,         &
                      & ftota_mixs,nloop,is,ie,rdir,rdif,tdir,tdif,     &
                      & explay,exptdn,rdndif,tottrn)

!-----------------------------------------------------------------------
!
! Computes layer reflectivities and transmissivities, from the top down
! to the surface using the delta-Eddington solutions for each layer;
! adds layers from top down to surface as well.
!
! If total transmission to the interface above a particular layer is
! less than trmin, then no further delta-Eddington solutions are
! evaluated for layers below
!
! For more details , see Briegleb, Bruce P., 1992: Delta-Eddington
! Approximation for Solar Radiation in the NCAR Community Climate Model,
! Journal of Geophysical Research, Vol 97, D7, pp7603-7612).
!
!---------------------------Code history--------------------------------
!
! Original version:  B. Briegleb
! Standardized:      J. Rosinski, June 1992
! Reviewed:          J. Kiehl, B. Briegleb, August 1992
!
!-----------------------------------------------------------------------
!
      implicit none
!
! PARAMETER definitions
!
!     Minimum total transmission below which no layer computation are
!     done:
!
! trmin  - Minimum total transmission allowed
! wray   - Rayleigh single scatter albedo
! gray   - Rayleigh asymetry parameter
! fray   - Rayleigh forward scattered fraction
!
      real(8) , parameter :: trmin = 1.E-3 , wray = 0.999999 ,          &
                           & gray = 0.0 , fray = 0.1
!
!------------------------------Arguments--------------------------------
!
!     Input arguments
!
! coszrs   - Cosine zenith angle
! trayoslp - Tray/sslp
! pflx     - Interface pressure
! abh2o    - Absorption coefficiant for h2o
! abo3     - Absorption coefficiant for o3
! abco2    - Absorption coefficiant for co2
! abo2     - Absorption coefficiant for o2
! uh2o     - Layer absorber amount of h2o
! uo3      - Layer absorber amount of  o3
! uco2     - Layer absorber amount of co2
! uo2      - Layer absorber amount of  o2
! tauxcl   - Cloud extinction optical depth
! wcl      - Cloud single scattering albedo
! gcl      - Cloud assymetry parameter
! fcl      - Cloud forward scattered fraction
! tauxci   - Cloud extinction optical depth
! wci      - Cloud single scattering albedo
! gci      - Cloud assymetry parameter
! fci      - Cloud forward scattered fraction
! nloop    - Number of loops (1 or 2)
! is       - Starting index for 1 or 2 loops
! ie       - Ending index for 1 or 2 loops
!
!     Input/Output arguments
!
!     Following variables are defined for each layer; 0 refers to extra
!     layer above top of model:
!
! rdir     - Layer reflectivity to direct rad
! rdif     - Layer refflectivity to diffuse mod_rad
! tdir     - Layer transmission to direct rad
! tdif     - Layer transmission to diffuse mod_rad
! explay   - Solar beam exp transm for layer
!
!     (Note that the following variables are defined on interfaces,
!     with the index k referring to the top interface of the kth layer:
!     exptdn,rdndif,tottrn; for example, tottrn(k=5) refers to the total
!     transmission to the top interface of the 5th layer; kzp1 refers
!     to the earth surface)
!
! rdndif   - Added dif ref for layers above
! exptdn   - Solar beam exp down transm from top
! tottrn   - Total transmission for layers above
!
!
! Dummy arguments
!
      real(8) :: abco2 , abh2o , abo2 , abo3 , trayoslp
      integer :: nloop
      real(8) , dimension(iym1) :: coszrs
      real(8) , dimension(iym1,0:kz) :: explay , fci , fcl ,         &
           & ftota_mixs , gci , gcl , gtota_mixs , rdif , rdir ,        &
           & tauasc_mixs , tauxar_mixs , tauxci , tauxcl , tdif , tdir ,&
           & uco2 , uh2o , uo2 , uo3 , wci , wcl
      real(8) , dimension(iym1,0:kzp1) :: exptdn , pflx , rdndif ,    &
           & tottrn
      integer , dimension(2) :: ie , is
      intent (in) abco2 , abh2o , abo2 , abo3 , coszrs , fci , fcl ,    &
                & ftota_mixs , gci , gcl , gtota_mixs , ie , is ,       &
                & nloop , pflx , tauasc_mixs , tauxar_mixs , tauxci ,   &
                & tauxcl , trayoslp , uco2 , uh2o , uo2 , uo3 , wci ,   &
                & wcl
      intent (inout) explay , exptdn , rdif , rdir , rdndif , tdif ,    &
                   & tdir , tottrn
!
!---------------------------Local variables-----------------------------
!
! i        - Longitude index
! k        - Level index
! nn       - Index of longitude loops (max=nloop)
! ii       - Longitude index
! nval     - Number of long values satisfying criteria
! indx     - Array of longitude indices
! taugab   - Layer total gas absorption optical depth
! tauray   - Layer rayleigh optical depth
! taucsc   - Layer cloud scattering optical depth
! tautot   - Total layer optical depth
! wtot     - Total layer single scatter albedo
! gtot     - Total layer asymmetry parameter
! ftot     - Total layer forward scatter fraction
! wtau     - rayleigh layer scattering optical depth
! wt       - layer total single scattering albedo
! ts       - layer scaled extinction optical depth
! ws       - layer scaled single scattering albedo
! gs       - layer scaled asymmetry parameter
! rdenom   - mulitiple scattering term
! rdirexp  - layer direct ref times exp transmission
! tdnmexp  - total transmission minus exp transmission
!
!---------------------------Statement functions-------------------------
!
!     Statement functions and other local variables
!
! alpha    - Term in direct reflect and transmissivity
! xgamm    - Term in direct reflect and transmissivity
! el       - Term in alpha,xgamm,n,u
! taus     - Scaled extinction optical depth
! omgs     - Scaled single particle scattering albedo
! asys     - Scaled asymmetry parameter
! u        - Term in diffuse reflect and transmissivity
! n        - Term in diffuse reflect and transmissivity
! lm       - Temporary for el
! ne       - Temporary for n
! w        - Dummy argument for statement function
! uu       - Dummy argument for statement function
! g        - Dummy argument for statement function
! e        - Dummy argument for statement function
! f        - Dummy argument for statement function
! t        - Dummy argument for statement function
! et       - Dummy argument for statement function
!
!     Intermediate terms for delta-eddington solution
!
! alp      - Temporary for alpha
! gam      - Temporary for xgamm
! ue       - Temporary for u
! arg      - Exponential argument
! extins   - Extinction
! amg      - Alp - gam
! apg      - Alp + gam
!
! Local variables
!
      real(8) :: alp , amg , apg , arg , e , et , extins , f , ftot ,   &
               & g , gam , gs , gtot , lm , ne , rdenom , rdirexp , t , &
               & taucsc , tautot , tdnmexp , ts , ue , uu , w , ws ,    &
               & wt , wtau , wtot
      real(8) :: alpha , asys , el , xgamm , n , omgs , taus , u
      integer :: i , ii , k , nn , nval
      integer , dimension(iym1) :: indx
      real(8) , dimension(iym1) :: taugab , tauray
!
      alpha(w,uu,g,e) = .75*w*uu*((1.+g*(1-w))/(1.-e*e*uu*uu))
      xgamm(w,uu,g,e) = .50*w*((3.*g*(1.-w)*uu*uu+1.)/(1.-e*e*uu*uu))
      el(w,g) = dsqrt(3.*(1-w)*(1.-w*g))
      taus(w,f,t) = (1.-w*f)*t
      omgs(w,f) = (1.-f)*w/(1.-w*f)
      asys(g,f) = (g-f)/(1.-f)
      u(w,g,e) = 1.5*(1.-w*g)/e
      n(uu,et) = ((uu+1.)*(uu+1.)/et) - ((uu-1.)*(uu-1.)*et)
!
!-----------------------------------------------------------------------
!
!     Initialize all total transmission values to 0, so that nighttime
!     values from previous computations are not used:
!
      tottrn = 0.D0
!
!     Compute total direct beam transmission, total transmission, and
!     reflectivity for diffuse mod_radiation (from below) for all layers
!     above each interface by starting from the top and adding layers
!     down:
!     For the extra layer above model top:
!
      do nn = 1 , nloop
        do i = is(nn) , ie(nn)
!
          tauray(i) = trayoslp*(pflx(i,1)-pflx(i,0))
          taugab(i) = abh2o*uh2o(i,0) + abo3*uo3(i,0) + abco2*uco2(i,0) &
                    & + abo2*uo2(i,0)
!
          tautot = tauxcl(i,0) + tauxci(i,0) + tauray(i) + taugab(i)    &
                 & + tauxar_mixs(i,0)
          taucsc = tauxcl(i,0)*wcl(i,0) + tauxci(i,0)*wci(i,0)          &
                 & + tauasc_mixs(i,0)
          wtau = wray*tauray(i)
          wt = wtau + taucsc
          wtot = wt/tautot
          gtot = (wtau*gray+gcl(i,0)*tauxcl(i,0)*wcl(i,0)+gci(i,0)      &
               & *tauxci(i,0)*wci(i,0)+gtota_mixs(i,0))/wt
          ftot = (wtau*fray+fcl(i,0)*tauxcl(i,0)*wcl(i,0)+fci(i,0)      &
               & *tauxci(i,0)*wci(i,0)+ftota_mixs(i,0))/wt
!
          ts = taus(wtot,ftot,tautot)
          ws = omgs(wtot,ftot)
          gs = asys(gtot,ftot)
          lm = el(ws,gs)
          alp = alpha(ws,coszrs(i),gs,lm)
          gam = xgamm(ws,coszrs(i),gs,lm)
          ue = u(ws,gs,lm)
!
!         Limit argument of exponential to 25, in case lm*ts very large:
!
          arg = dmin1(lm*ts,25.D0)
          extins = dexp(-arg)
          ne = n(ue,extins)
!
          rdif(i,0) = (ue+1.)*(ue-1.)*(1./extins-extins)/ne
          tdif(i,0) = 4.*ue/ne
!
!         Limit argument of exponential to 25, in case coszrs is very
!         small:
          arg = dmin1(ts/coszrs(i),25.D0)
          explay(i,0) = dexp(-arg)
!
          apg = alp + gam
          amg = alp - gam
          rdir(i,0) = amg*(tdif(i,0)*explay(i,0)-1.) + apg*rdif(i,0)
          tdir(i,0) = apg*tdif(i,0) + (amg*rdif(i,0)-(apg-1.))          &
                    & *explay(i,0)
!
!         Under rare conditions, reflectivies and transmissivities can
!         be negative; zero out any negative values
!
          rdir(i,0) = dmax1(rdir(i,0),0.D0)
          tdir(i,0) = dmax1(tdir(i,0),0.D0)
          rdif(i,0) = dmax1(rdif(i,0),0.D0)
          tdif(i,0) = dmax1(tdif(i,0),0.D0)
!
!         Initialize top interface of extra layer:
!
          exptdn(i,0) = 1.0
          rdndif(i,0) = 0.0
          tottrn(i,0) = 1.0
!
          rdndif(i,1) = rdif(i,0)
          tottrn(i,1) = tdir(i,0)
!
        end do
      end do
!
!     Now, continue down one layer at a time; if the total transmission
!     to the interface just above a given layer is less than trmin,
!     then no delta-eddington computation for that layer is done:
!
      do k = 1 , kz
!
!       Initialize current layer properties to zero; only if total
!       transmission to the top interface of the current layer exceeds
!       the minimum, will these values be computed below:
!
        do nn = 1 , nloop
          do i = is(nn) , ie(nn)
!
            rdir(i,k) = 0.0
            rdif(i,k) = 0.0
            tdir(i,k) = 0.0
            tdif(i,k) = 0.0
            explay(i,k) = 0.0
!
!           Calculates the solar beam transmission, total transmission,
!           and reflectivity for diffuse mod_radiation from below at the
!           top of the current layer:
!
            exptdn(i,k) = exptdn(i,k-1)*explay(i,k-1)
!KN         modified below (for computational stability)
!           rdenom      = 1./(1. - rdif(i,k-1)*rdndif(i,k-1))
            rdenom = 1./(1.-dmin1(rdif(i,k-1)*rdndif(i,k-1),0.999999D0))
!KN         modified above
            rdirexp = rdir(i,k-1)*exptdn(i,k-1)
            tdnmexp = tottrn(i,k-1) - exptdn(i,k-1)
            tottrn(i,k) = exptdn(i,k-1)*tdir(i,k-1) + tdif(i,k-1)       &
                        & *(tdnmexp+rdndif(i,k-1)*rdirexp)*rdenom
            rdndif(i,k) = rdif(i,k-1) + (rdndif(i,k-1)*tdif(i,k-1))     &
                        & *(tdif(i,k-1)*rdenom)
!
          end do
        end do
!
!       Compute next layer delta-eddington solution only if total
!       transmission of radiation to the interface just above the layer
!       exceeds trmin.
        call whenfgt(iym1,tottrn(1,k),1,trmin,indx,nval)
        if ( nval.gt.0 ) then
!CDIR$    IVDEP
          do ii = 1 , nval
            i = indx(ii)
!
            tauray(i) = trayoslp*(pflx(i,k+1)-pflx(i,k))
            taugab(i) = abh2o*uh2o(i,k) + abo3*uo3(i,k)                 &
                      & + abco2*uco2(i,k) + abo2*uo2(i,k)
!
            tautot = tauxcl(i,k) + tauxci(i,k) + tauray(i) + taugab(i)  &
                   & + tauxar_mixs(i,k)
            taucsc = tauxcl(i,k)*wcl(i,k) + tauxci(i,k)*wci(i,k)        &
                   & + tauasc_mixs(i,k)
            wtau = wray*tauray(i)
            wt = wtau + taucsc
            wtot = wt/tautot
            gtot = (wtau*gray+gcl(i,k)*wcl(i,k)*tauxcl(i,k)+gci(i,k)    &
                 & *wci(i,k)*tauxci(i,k)+gtota_mixs(i,k))/wt
            ftot = (wtau*fray+fcl(i,k)*wcl(i,k)*tauxcl(i,k)+fci(i,k)    &
                 & *wci(i,k)*tauxci(i,k)+ftota_mixs(i,k))/wt
!
            ts = taus(wtot,ftot,tautot)
            ws = omgs(wtot,ftot)
            gs = asys(gtot,ftot)
            lm = el(ws,gs)
            alp = alpha(ws,coszrs(i),gs,lm)
            gam = xgamm(ws,coszrs(i),gs,lm)
            ue = u(ws,gs,lm)
!
!           Limit argument of exponential to 25, in case lm very large:
!
            arg = dmin1(lm*ts,25.D0)
            extins = dexp(-arg)
            ne = n(ue,extins)
!
            rdif(i,k) = (ue+1.)*(ue-1.)*(1./extins-extins)/ne
            tdif(i,k) = 4.*ue/ne
!
!           Limit argument of exponential to 25, in case coszrs is very
!           small:
            arg = dmin1(ts/coszrs(i),25.D0)
            explay(i,k) = dexp(-arg)
!
            apg = alp + gam
            amg = alp - gam
            rdir(i,k) = amg*(tdif(i,k)*explay(i,k)-1.) + apg*rdif(i,k)
            tdir(i,k) = apg*tdif(i,k) + (amg*rdif(i,k)-(apg-1.))        &
                      & *explay(i,k)
!
!           Under rare conditions, reflectivies and transmissivities
!           can be negative; zero out any negative values
!
            rdir(i,k) = dmax1(rdir(i,k),0.D0)
            tdir(i,k) = dmax1(tdir(i,k),0.D0)
            rdif(i,k) = dmax1(rdif(i,k),0.D0)
            tdif(i,k) = dmax1(tdif(i,k),0.D0)
          end do
        end if
!
      end do
!
!     Compute total direct beam transmission, total transmission, and
!     reflectivity for diffuse mod_radiation (from below) for all layers
!     above the surface:
!
      k = kzp1
      do nn = 1 , nloop
        do i = is(nn) , ie(nn)
          exptdn(i,k) = exptdn(i,k-1)*explay(i,k-1)
!KN       modified below (for computational stability)
!         rdenom = 1./(1. - rdif(i,k-1)*rdndif(i,k-1))
          rdenom = 1./(1.-dmin1(rdif(i,k-1)*rdndif(i,k-1),0.999999D0))
!KN       modified above
          rdirexp = rdir(i,k-1)*exptdn(i,k-1)
          tdnmexp = tottrn(i,k-1) - exptdn(i,k-1)
          tottrn(i,k) = exptdn(i,k-1)*tdir(i,k-1) + tdif(i,k-1)         &
                      & *(tdnmexp+rdndif(i,k-1)*rdirexp)*rdenom
          rdndif(i,k) = rdif(i,k-1) + (rdndif(i,k-1)*tdif(i,k-1))       &
                      & *(tdif(i,k-1)*rdenom)
        end do
      end do
!
      end subroutine radded
!
      subroutine radabs(pbr,pnm,co2em,co2eml,tplnka,s2c,s2t,w,h2otr,    &
                      & plco2,plh2o,co2t,tint,tlayr,plol,plos,pmln,piln,&
                      & ucfc11,ucfc12,un2o0,un2o1,uch4,uco211,uco212,   &
                      & uco213,uco221,uco222,uco223,uptype,bn2o0,bn2o1, &
                      & bch4,abplnk1,abplnk2,jslc)

!-----------------------------------------------------------------------
!
! Compute absorptivities for h2o, co2, o3, ch4, n2o, cfc11 and cfc12
!
! h2o  ....  Uses nonisothermal emissivity for water vapor from
!            Ramanathan, V. and  P.Downey, 1986: A Nonisothermal
!            Emissivity and Absorptivity Formulation for Water Vapor
!            Journal of Geophysical Research, vol. 91., D8, pp 8649-8666
!
! co2  ....  Uses absorptance parameterization of the 15 micro-meter
!            (500 - 800 cm-1) band system of Carbon Dioxide, from
!            Kiehl, J.T. and B.P.Briegleb, 1991: A New Parameterization
!            of the Absorptance Due to the 15 micro-meter Band System
!            of Carbon Dioxide Jouranl of Geophysical Research,
!            vol. 96., D5, pp 9013-9019.
!            Parameterizations for the 9.4 and 10.4 mircon bands of CO2
!            are also included.
!
! o3   ....  Uses absorptance parameterization of the 9.6 micro-meter
!            band system of ozone, from Ramanathan, V. and R.Dickinson,
!            1979: The Role of stratospheric ozone in the zonal and
!            seasonal radiative energy balance of the earth-troposphere
!            system. Journal of the Atmospheric Sciences, Vol. 36,
!            pp 1084-1104
!
! ch4  ....  Uses a broad band model for the 7.7 micron band of methane.
!
! n20  ....  Uses a broad band model for the 7.8, 8.6 and 17.0 micron
!            bands of nitrous oxide
!
! cfc11 ...  Uses a quasi-linear model for the 9.2, 10.7, 11.8 and 12.5
!            micron bands of CFC11
!
! cfc12 ...  Uses a quasi-linear model for the 8.6, 9.1, 10.8 and 11.2
!            micron bands of CFC12
!
!
! Computes individual absorptivities for non-adjacent layers, accounting
! for band overlap, and sums to obtain the total; then, computes the
! nearest layer contribution.
!
!---------------------------Code history--------------------------------
!
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          J. Kiehl, B. Briegleb, August 1992
! Reviewed:          J. Kiehl, April 1996
! Reviewed:          B. Briegleb, May 1996
!
!-----------------------------------------------------------------------
!
      implicit none
!
!------------------------------Arguments--------------------------------
!
!     Input arguments
!
! pbr      - Prssr at mid-levels (dynes/cm2)
! pnm      - Prssr at interfaces (dynes/cm2)
! co2em    - Co2 emissivity function
! co2eml   - Co2 emissivity function
! tplnka   - Planck fnctn level temperature
! s2c      - H2o continuum path length
! s2t      - H2o tmp and prs wghted path
! w        - H2o prs wghted path
! h2otr    - H2o trnsmssn fnct for o3 overlap
! plco2    - Co2 prs wghted path length
! plh2o    - H2o prs wfhted path length
! co2t     - Tmp and prs wghted path length
! tint     - Interface temperatures
! tlayr    - K-1 level temperatures
! plol     - Ozone prs wghted path length
! plos     - Ozone path length
! pmln     - Ln(pmidm1)
! piln     - Ln(pintm1)
!
!     Trace gas variables
!
! ucfc11   - CFC11 path length
! ucfc12   - CFC12 path length
! un2o0    - N2O path length
! un2o1    - N2O path length (hot band)
! uch4     - CH4 path length
! uco211   - CO2 9.4 micron band path length
! uco212   - CO2 9.4 micron band path length
! uco213   - CO2 9.4 micron band path length
! uco221   - CO2 10.4 micron band path length
! uco222   - CO2 10.4 micron band path length
! uco223   - CO2 10.4 micron band path length
! uptype   - continuum path length
! bn2o0    - pressure factor for n2o
! bn2o1    - pressure factor for n2o
! bch4     - pressure factor for ch4
! abplnk1   - non-nearest layer Plack factor
! abplnk2   - nearest layer factor
! abstrc    - total trace gas absorptivity
! bplnk     - Planck functions for sub-divided layers
!
!     Output arguments (radbuf)
!
! Nearest layer absorptivities
! Non-adjacent layer absorptivites
! Total emissivity
!
! Dummy arguments
!
      integer :: jslc
      real(8) , dimension(14,iym1,kzp1) :: abplnk1 , abplnk2
      real(8) , dimension(iym1,kzp1) :: bch4 , bn2o0 , bn2o1 , co2em ,&
           & co2t , h2otr , piln , plco2 , plh2o , plol , plos , pnm ,  &
           & s2c , s2t , tint , tlayr , tplnka , ucfc11 , ucfc12 ,      &
           & uch4 , uco211 , uco212 , uco213 , uco221 , uco222 ,        &
           & uco223 , un2o0 , un2o1 , uptype , w
      real(8) , dimension(iym1,kz) :: co2eml , pbr , pmln
      intent (in) abplnk2 , co2em , co2eml , co2t , h2otr , jslc , pbr ,&
                & piln , plco2 , plh2o , plol , plos , pmln , s2t ,     &
                & tint , tlayr , tplnka , w
!
!---------------------------Local variables-----------------------------
!
! i        - Longitude index
! k        - Level index
! k1       - Level index
! k2       - Level index
! kn       - Nearest level index
! iband    - Band  index
! pnew     - Effective pressure for H2O vapor linewidth
! trline   - Transmission due to H2O lines in window
! u        - Pressure weighted H2O path length
! tbar     - Mean layer temperature
! emm      - Mean co2 emissivity
! o3emm    - Mean o3 emissivity
! o3bndi   - Ozone band parameter
! temh2o   - Mean layer temperature equivalent to tbar
! k21      - Exponential coefficient used to calculate rotation band
!            transmissvty in the 650-800 cm-1 region (tr1)
! k22      - Exponential coefficient used to calculate rotation band
!            transmissvty in the 500-650 cm-1 region (tr2)
! uc1      - H2o continuum pathlength in 500-800 cm-1
! to3h2o   - H2o trnsmsn for overlap with o3
! pi       - For co2 absorptivity computation
! sqti     - Used to store sqrt of mean temperature
! et       - Co2 hot band factor
! et2      - Co2 hot band factor squared
! et4      - Co2 hot band factor to fourth power
! omet     - Co2 stimulated emission term
! f1co2    - Co2 central band factor
! f2co2    - Co2 weak band factor
! f3co2    - Co2 weak band factor
! t1co2    - Overlap factr weak bands on strong band
! sqwp     - Sqrt of co2 pathlength
! f1sqwp   - Main co2 band factor
! oneme    - Co2 stimulated emission term
! alphat   - Part of the co2 stimulated emission term
! wco2     - Constants used to define co2 pathlength
! posqt    - Effective pressure for co2 line width
! u7       - Co2 hot band path length
! u8       - Co2 hot band path length
! u9       - Co2 hot band path length
! u13      - Co2 hot band path length
! rbeta7   - Inverse of co2 hot band line width par
! rbeta8   - Inverse of co2 hot band line width par
! rbeta9   - Inverse of co2 hot band line width par
! rbeta13  - Inverse of co2 hot band line width par
! tpatha   - For absorptivity computation
! a        - Eq(2) in table A3a of R&D
! abso     - Absorptivity for various gases/bands
! dtp      - Path temp minus 300 K used in h2o rotation band absorptivity
! dtx      - Planck temperature minus 250 K
! dty      - Path temperature minus 250 K
! dtz      - Planck temperature minus 300 K
! term1    - Equation(5) in table A3a of R&D(1986)
! term2    - Delta a(Te) in table A3a of R&D(1986)
! term3    - DB/dT function for rotation and
! term4    - Equation(6) in table A3a of R&D(1986)
! term5    - Delta a(Tp) in table A3a of R&D(1986)
! term6    - DB/dT function for window region
! term7    - Kl_inf(i) in eq(8) of table A3a of R&D
! term8    - Delta kl_inf(i) in eq(8)
! term9    - DB/dT function for 500-800 cm-1 region
! tr1      - Eqn(6) in table A2 of R&D for 650-800
! tr10     - Eqn (6) times eq(4) in table A2 of R&D for 500-650 cm-1 region
! tr2      - Eqn(6) in table A2 of R&D for 500-650
! tr5      - Eqn(4) in table A2 of R&D for 650-800
! tr6      - Eqn(4) in table A2 of R&D for 500-650
! tr9      - Equation (6) times eq(4) in table A2 of R&D for 650-800 cm-1 region
! uc       - Y + 0.002U in eq(8) of table A2 of R&D
! sqrtu    - Sqrt of pressure weighted h20 pathlength
! fwk      - Equation(33) in R&D far wing correction
! fwku     - GU term in eqs(1) and (6) in table A2
! r2st     - 1/(2*beta) in eq(10) in table A2
! dtyp15   - DeltaTp in eqs(11) & (12) in table A3a
! dtyp15sq - (DeltaTp)^2 in eqs(11) & (12) table A3a
! to3co2   - P weighted temp in ozone band model
! dpnm     - Pressure difference between two levels
! pnmsq    - Pressure squared
! dw       - Amount of h2o between two levels
! uinpl    - Nearest layer subdivision factor
! winpl    - Nearest layer subdivision factor
! zinpl    - Nearest layer subdivision factor
! pinpl    - Nearest layer subdivision factor
! dplh2o   - Difference in press weighted h2o amount
! r80257   - Conversion factor for h2o pathlength
! r293     - 1/293
! r250     - 1/250
! r3205    - Line width factor for o3 (see R&Di)
! r300     - 1/300
! rsslp    - Reciprocal of sea level pressure
! r2sslp   - 1/2 of rsslp
! ds2c     - Y in eq(7) in table A2 of R&D
! a11      - A1 in table A3b for rotation band absorptivity
! a31      - A3 in table A3b for rotation band absorptivity
! a21      - First part in numerator of A2 in table A3b
! a22      - Second part in numerator of A2 in table A3b
! a23      - Denominator of A2 in table A3b (rotation band)
! t1t4     - Eq(3) in table A3a of R&D
! t2t5     - Eq(4) in table A3a of R&D
! rsum     - Eq(1) in table A2 of R&D
! a41      - Numerator in A2 in Vib-rot abstivity(table A3b)
! a51      - Denominator in A2 in Vib-rot (table A3b)
! a61      - A3 factor for Vib-rot band in table A3b
! phi      - Eq(11) in table A3a of R&D
! psi      - Eq(12) in table A3a of R&D
! cf812    - Eq(11) in table A2 of R&D
! ubar     - H2o scaled path see comment for eq(10) table A2
! pbar     - H2o scaled pres see comment for eq(10) table A2
! g4       - Arguement in exp() in eq(10) table A2
! dplos    - Ozone pathlength eq(A2) in R&Di
! dplol    - Presure weighted ozone pathlength
! beta     - Local interface temperature (includes Voigt line correction factor)
! rphat    - Effective pressure for ozone beta
! tcrfac   - Ozone temperature factor table 1 R&Di
! tmp1     - Ozone band factor see eq(A1) in R&Di
! u1       - Effective ozone pathlength eq(A2) in R&Di
! realnu   - 1/beta factor in ozone band model eq(A1)
! tmp2     - Ozone band factor see eq(A1) in R&Di
! u2       - Effective ozone pathlength eq(A2) in R&Di
! rsqti    - Reciprocal of sqrt of path temperature
! tpath    - Path temperature used in co2 band model
! tmp3     - Weak band factor see K&B
! rdpnmsq  - Reciprocal of difference in press^2
! rdpnm    - Reciprocal of difference in press
! p1       - Mean pressure factor
! p2       - Mean pressure factor
! dtym10   - T - 260 used in eq(9) and (10) table A3a
! dplco2   - Co2 pathlength
! corfac   - Correction factors in table A3b
! g2       - Part of arguement in eq(10) in table A2
! te       - A_0 T factor in ozone model table 1 of R&Di
! denom    - Denominator in eq(8) of table A3a of R&D
! trab2    - Transmission terms for H2o  500 -  800 cm-1
! trab4    - Transmission terms for H2o  800 - 1000 cm-1
! trab6    - Transmission terms for H2o 1000 - 1200 cm-1
! absbnd   - Proportional to co2 band absorptance
! dbvtit   - Intrfc drvtv plnck fnctn for o3
! dbvtly   - Level drvtv plnck fnctn for o3
! dbvt     - Planck fnctn tmp derivative for o3
!
!
! Local variables
!
      real(8) :: a , a11 , a21 , a22 , a23 , a31 , a41 , a51 , a61 ,    &
               & absbnd , alphat , beta , cf812 , corfac , denom ,      &
               & dplco2 , dplol , dplos , ds2c , dtym10 , et , et2 ,    &
               & et4 , f1co2 , g2 , g4 , k21 , k22 , o3bndi , omet ,    &
               & oneme , p1 , p2 , pbar , phi , pi , posqt , psi ,      &
               & r250 , r293 , r2sslp , r300 , r3205 , r80257 ,         &
               & rbeta13 , rbeta8 , rbeta9 , rdpnm , rdpnmsq , realnu , &
               & rphat , rsqti , rsslp , rsum , sqwp , t , t1t4 , t2t5 ,&
               & tcrfac , te , tlocal , tmp1 , tmp2 , tmp3 , tpath ,    &
               & tr1 , tr2 , tr5 , tr6 , u1 , u13 , u2 , u8 , u9 ,      &
               & ubar , wco2
      real(8) , dimension(iym1,6) :: abso
      real(8) , dimension(iym1) :: abstrc , dplh2o , dpnm , dtp , dtx ,&
                                  & dty , dtyp15 , dtyp15sq , dtz , dw ,&
                                  & f1sqwp , f2co2 , f3co2 , fwk ,      &
                                  & fwku , pnew , rbeta7 , sqrtu ,      &
                                  & sqti , t1co2 , tco2 , th2o , to3 ,  &
                                  & to3co2 , to3h2o , tpatha , tr10 ,   &
                                  & tr9 , trab2 , trab4 , trab6 , u ,   &
                                  & u7 , uc , uc1
      real(8) , dimension(14,iym1,4) :: bplnk
      real(8) :: dbvt
      real(8) , dimension(iym1,kzp1) :: dbvtit , pnmsq , term6 , term9
      real(8) , dimension(iym1,kz) :: dbvtly
      real(8) , dimension(iym1,4) :: emm , o3emm , pinpl , tbar ,      &
                                    & temh2o , term1 , term2 , term3 ,  &
                                    & term4 , term5 , uinpl , winpl ,   &
                                    & zinpl
      integer :: i , iband , k , k1 , k2 , kn , wvl
      real(8) , dimension(2) :: r2st
      real(8) , dimension(iym1,2) :: term7 , term8 , trline
!
!--------------------------Statement function---------------------------
!
      dbvt(t) = (-2.8911366682E-4+(2.3771251896E-6+1.1305188929E-10*t)  &
              & *t)/(1.0+(-6.1364820707E-3+1.5550319767E-5*t)*t)
!
!-----------------------------------------------------------------------
!
!     Initialize
!
      do k = 1 , kz
        do i = 1 , iym1
          dbvtly(i,k) = dbvt(tlayr(i,k+1))
          dbvtit(i,k) = dbvt(tint(i,k))
        end do
      end do
      do i = 1 , iym1
        dbvtit(i,kzp1) = dbvt(tint(i,kz + 1))
      end do
!
      r80257 = 1./8.0257D-04
      r293 = 1./293.D0
      r250 = 1./250.D0
      r3205 = 1./.3205D0
      r300 = 1./300.D0
      rsslp = 1./sslp
      r2sslp = 1./(2.*sslp)
      r2st(1) = 1./(2.*st(1))
      r2st(2) = 1./(2.*st(2))
!     bndfct  = 2.0*22.18d0/(dsqrt(196.d0)*300.)
!
!     Non-adjacent layer absorptivity:
!
!     abso(i,1)     0 -  800 cm-1   h2o rotation band
!     abso(i,2)  1200 - 2200 cm-1   h2o vibration-rotation band
!     abso(i,3)   800 - 1200 cm-1   h2o window
!     abso(i,4)   500 -  800 cm-1   h2o rotation band overlap with co2
!     abso(i,5)   o3  9.6 micrometer band (nu3 and nu1 bands)
!     abso(i,6)   co2 15  micrometer band system
!
      do k = 1 , kzp1
        do i = 1 , iym1
          pnmsq(i,k) = pnm(i,k)**2
          dtx(i) = tplnka(i,k) - 250.
          term6(i,k) = coeff(1,2) + coeff(2,2)*dtx(i)                   &
                     & *(1.+c9*dtx(i)*(1.+c11*dtx(i)                    &
                     & *(1.+c13*dtx(i)*(1.+c15*dtx(i)))))
          term9(i,k) = coefi(1,2) + coefi(2,2)*dtx(i)                   &
                     & *(1.+c19*dtx(i)*(1.+c21*dtx(i)                   &
                     & *(1.+c23*dtx(i)*(1.+c25*dtx(i)))))
        end do
      end do
!
!     Non-nearest layer level loops
!
      do k1 = kzp1 , 1 , -1
        do k2 = kzp1 , 1 , -1
          if ( k1.ne.k2 ) then
            do i = 1 , iym1
              dplh2o(i) = plh2o(i,k1) - plh2o(i,k2)
              u(i) = dabs(dplh2o(i))
              sqrtu(i) = dsqrt(u(i))
              ds2c = dabs(s2c(i,k1)-s2c(i,k2))
              dw(i) = dabs(w(i,k1)-w(i,k2))
              uc1(i) = (ds2c+1.7E-3*u(i))*(1.+2.*ds2c)/(1.+15.*ds2c)
              uc(i) = ds2c + 2.E-3*u(i)
            end do
            do i = 1 , iym1
              pnew(i) = u(i)/dw(i)
              tpatha(i) = (s2t(i,k1)-s2t(i,k2))/dplh2o(i)
              dtx(i) = tplnka(i,k2) - 250.
              dty(i) = tpatha(i) - 250.
              dtyp15(i) = dty(i) + 15.
              dtyp15sq(i) = dtyp15(i)**2
              dtz(i) = dtx(i) - 50.
              dtp(i) = dty(i) - 50.
            end do
            do iband = 2 , 4 , 2
              do i = 1 , iym1
                term1(i,iband) = coefe(1,iband) + coefe(2,iband)*dtx(i) &
                               & *(1.+c1(iband)*dtx(i))
                term2(i,iband) = coefb(1,iband) + coefb(2,iband)*dtx(i) &
                               & *(1.+c2(iband)*dtx(i)                  &
                               & *(1.+c3(iband)*dtx(i)))
                term3(i,iband) = coefd(1,iband) + coefd(2,iband)*dtx(i) &
                               & *(1.+c4(iband)*dtx(i)                  &
                               & *(1.+c5(iband)*dtx(i)))
                term4(i,iband) = coefa(1,iband) + coefa(2,iband)*dty(i) &
                               & *(1.+c6(iband)*dty(i))
                term5(i,iband) = coefc(1,iband) + coefc(2,iband)*dty(i) &
                               & *(1.+c7(iband)*dty(i))
              end do
            end do
!
!           abso(i,1)     0 -  800 cm-1   h2o rotation band
!
            do i = 1 , iym1
              a11 = 0.44 + 3.380E-4*dtz(i) - 1.520E-6*dtz(i)*dtz(i)
              a31 = 1.05 - 6.000E-3*dtp(i) + 3.000E-6*dtp(i)*dtp(i)
              a21 = 1.00 + 1.717E-3*dtz(i) - 1.133E-5*dtz(i)*dtz(i)
              a22 = 1.00 + 4.443E-3*dtp(i) + 2.750E-5*dtp(i)*dtp(i)
              a23 = 1.00 + 3.600*sqrtu(i)
              corfac = a31*(a11+((2.*a21*a22)/a23))
              t1t4 = term1(i,2)*term4(i,2)
              t2t5 = term2(i,2)*term5(i,2)
              a = t1t4 + t2t5/(1.+t2t5*sqrtu(i)*corfac)
              fwk(i) = fwcoef + fwc1/(1.+fwc2*u(i))
              fwku(i) = fwk(i)*u(i)
              rsum = dexp(-a*(sqrtu(i)+fwku(i)))
              abso(i,1) = (1.-rsum)*term3(i,2)
!             trab1(i)  = rsum
            end do
!
!           abso(i,2)  1200 - 2200 cm-1   h2o vibration-rotation band
!
            do i = 1 , iym1
              a41 = 1.75 - 3.960E-03*dtz(i)
              a51 = 1.00 + 1.3*sqrtu(i)
              a61 = 1.00 + 1.250E-03*dtp(i) + 6.250E-05*dtp(i)*dtp(i)
              corfac = .29*(1.+a41/a51)*a61
              t1t4 = term1(i,4)*term4(i,4)
              t2t5 = term2(i,4)*term5(i,4)
              a = t1t4 + t2t5/(1.+t2t5*sqrtu(i)*corfac)
              rsum = dexp(-a*(sqrtu(i)+fwku(i)))
              abso(i,2) = (1.-rsum)*term3(i,4)
!             trab7(i)  = rsum
            end do
!
!           Line transmission in 800-1000 and 1000-1200 cm-1 intervals
!
            do k = 1 , 2
              do i = 1 , iym1
                phi = dexp(a1(k)*dtyp15(i)+a2(k)*dtyp15sq(i))
                psi = dexp(b1(k)*dtyp15(i)+b2(k)*dtyp15sq(i))
                ubar = dw(i)*phi*1.66*r80257
                pbar = pnew(i)*(psi/phi)
                cf812 = cfa1 + (1.-cfa1)/(1.+ubar*pbar*10.)
                g2 = 1. + ubar*4.0*st(k)*cf812/pbar
                g4 = realk(k)*pbar*r2st(k)*(dsqrt(g2)-1.)
                trline(i,k) = dexp(-g4)
              end do
            end do
            do i = 1 , iym1
              term7(i,1) = coefj(1,1) + coefj(2,1)*dty(i)               &
                         & *(1.+c16*dty(i))
              term8(i,1) = coefk(1,1) + coefk(2,1)*dty(i)               &
                         & *(1.+c17*dty(i))
              term7(i,2) = coefj(1,2) + coefj(2,2)*dty(i)               &
                         & *(1.+c26*dty(i))
              term8(i,2) = coefk(1,2) + coefk(2,2)*dty(i)               &
                         & *(1.+c27*dty(i))
            end do
!
!           abso(i,3)   800 - 1200 cm-1   h2o window
!           abso(i,4)   500 -  800 cm-1   h2o rotation band overlap
!           with co2
            do i = 1 , iym1
              k21 = term7(i,1) + term8(i,1)                             &
                  & /(1.+(c30+c31*(dty(i)-10.)*(dty(i)-10.))*sqrtu(i))
              k22 = term7(i,2) + term8(i,2)                             &
                  & /(1.+(c28+c29*(dty(i)-10.))*sqrtu(i))
              tr1 = dexp(-(k21*(sqrtu(i)+fc1*fwku(i))))
              tr2 = dexp(-(k22*(sqrtu(i)+fc1*fwku(i))))
              tr5 = dexp(-((coefh(1,3)+coefh(2,3)*dtx(i))*uc1(i)))
              tr6 = dexp(-((coefh(1,4)+coefh(2,4)*dtx(i))*uc1(i)))
              tr9(i) = tr1*tr5
              tr10(i) = tr2*tr6
              th2o(i) = tr10(i)
              trab2(i) = 0.65*tr9(i) + 0.35*tr10(i)
              trab4(i) = dexp(-(coefg(1,3)+coefg(2,3)*dtx(i))*uc(i))
              trab6(i) = dexp(-(coefg(1,4)+coefg(2,4)*dtx(i))*uc(i))
              abso(i,3) = term6(i,k2)                                   &
                        & *(1.-.5*trab4(i)*trline(i,2)-.5*trab6(i)      &
                        & *trline(i,1))
              abso(i,4) = term9(i,k2)*.5*(tr1-tr9(i)+tr2-tr10(i))
            end do
            if ( k2.lt.k1 ) then
              do i = 1 , iym1
                to3h2o(i) = h2otr(i,k1)/h2otr(i,k2)
              end do
            else
              do i = 1 , iym1
                to3h2o(i) = h2otr(i,k2)/h2otr(i,k1)
              end do
            end if
!
!           abso(i,5)   o3  9.6 micrometer band (nu3 and nu1 bands)
!
            do i = 1 , iym1
              dpnm(i) = pnm(i,k1) - pnm(i,k2)
              to3co2(i) = (pnm(i,k1)*co2t(i,k1)-pnm(i,k2)*co2t(i,k2))   &
                        & /dpnm(i)
              te = (to3co2(i)*r293)**.7
              dplos = plos(i,k1) - plos(i,k2)
              dplol = plol(i,k1) - plol(i,k2)
              u1 = 18.29*dabs(dplos)/te
              u2 = .5649*dabs(dplos)/te
              rphat = dplol/dplos
              tlocal = tint(i,k2)
              tcrfac = dsqrt(tlocal*r250)*te
              beta = r3205*(rphat+dpfo3*tcrfac)
              realnu = te/beta
              tmp1 = u1/dsqrt(4.+u1*(1.+realnu))
              tmp2 = u2/dsqrt(4.+u2*(1.+realnu))
              o3bndi = 74.*te*dlog(1.+tmp1+tmp2)
              abso(i,5) = o3bndi*to3h2o(i)*dbvtit(i,k2)
              to3(i) = 1.0/(1.+0.1*tmp1+0.1*tmp2)
!             trab5(i)  = 1.-(o3bndi/(1060-980.))
            end do
!
!           abso(i,6)      co2 15  micrometer band system
!
            do i = 1 , iym1
              sqwp = dsqrt(dabs(plco2(i,k1)-plco2(i,k2)))
              et = dexp(-480./to3co2(i))
              sqti(i) = dsqrt(to3co2(i))
              rsqti = 1./sqti(i)
              et2 = et*et
              et4 = et2*et2
              omet = 1. - 1.5*et2
              f1co2 = 899.70*omet*(1.+1.94774*et+4.73486*et2)*rsqti
              f1sqwp(i) = f1co2*sqwp
              t1co2(i) = 1./(1.+(245.18*omet*sqwp*rsqti))
              oneme = 1. - et2
              alphat = oneme**3*rsqti
              pi = dabs(dpnm(i))
              wco2 = 2.5221*co2vmr*pi*rga
              u7(i) = 4.9411E4*alphat*et2*wco2
              u8 = 3.9744E4*alphat*et4*wco2
              u9 = 1.0447E5*alphat*et4*et2*wco2
              u13 = 2.8388E3*alphat*et4*wco2
              tpath = to3co2(i)
              tlocal = tint(i,k2)
              tcrfac = dsqrt(tlocal*r250*tpath*r300)
              posqt = ((pnm(i,k2)+pnm(i,k1))*r2sslp+dpfco2*tcrfac)*rsqti
              rbeta7(i) = 1./(5.3228*posqt)
              rbeta8 = 1./(10.6576*posqt)
              rbeta9 = rbeta7(i)
              rbeta13 = rbeta9
              f2co2(i) = (u7(i)/dsqrt(4.+u7(i)*(1.+rbeta7(i))))         &
                       & + (u8/dsqrt(4.+u8*(1.+rbeta8)))                &
                       & + (u9/dsqrt(4.+u9*(1.+rbeta9)))
              f3co2(i) = u13/dsqrt(4.+u13*(1.+rbeta13))
            end do
            if ( k2.ge.k1 ) then
              do i = 1 , iym1
                sqti(i) = dsqrt(tlayr(i,k2))
              end do
            end if
!
            do i = 1 , iym1
              tmp1 = dlog(1.+f1sqwp(i))
              tmp2 = dlog(1.+f2co2(i))
              tmp3 = dlog(1.+f3co2(i))
              absbnd = (tmp1+2.*t1co2(i)*tmp2+2.*tmp3)*sqti(i)
              abso(i,6) = trab2(i)*co2em(i,k2)*absbnd
              tco2(i) = 1./(1.0+10.0*(u7(i)/dsqrt(4.+u7(i)*(1.+rbeta7(i)&
                      & ))))
!             trab3(i)  = 1. - bndfct*absbnd
            end do
!
!           Calculate absorptivity due to trace gases
!
            call trcab(k1,k2,ucfc11,ucfc12,un2o0,un2o1,uch4,uco211,     &
                     & uco212,uco213,uco221,uco222,uco223,bn2o0,bn2o1,  &
                     & bch4,to3co2,pnm,dw,pnew,s2c,uptype,u,abplnk1,    &
                     & tco2,th2o,to3,abstrc)
!
!           Sum total absorptivity
!
            do i = 1 , iym1
              abstot(i,k1,k2,jslc) = abso(i,1) + abso(i,2) + abso(i,3)  &
                                   & + abso(i,4) + abso(i,5) + abso(i,6)&
                                   & + abstrc(i)
            end do
          end if
        end do
      end do          ! End of non-nearest layer level loops
!
!     Non-adjacent layer absorptivity:
!
!     abso(i,1)     0 -  800 cm-1   h2o rotation band
!     abso(i,2)  1200 - 2200 cm-1   h2o vibration-rotation band
!     abso(i,3)   800 - 1200 cm-1   h2o window
!     abso(i,4)   500 -  800 cm-1   h2o rotation band overlap with co2
!     abso(i,5)   o3  9.6 micrometer band (nu3 and nu1 bands)
!     abso(i,6)   co2 15  micrometer band system
!
!     Nearest layer level loop
!
      do k2 = kz , 1 , -1
        do i = 1 , iym1
          tbar(i,1) = 0.5*(tint(i,k2+1)+tlayr(i,k2+1))
          emm(i,1) = 0.5*(co2em(i,k2+1)+co2eml(i,k2))
          tbar(i,2) = 0.5*(tlayr(i,k2+1)+tint(i,k2))
          emm(i,2) = 0.5*(co2em(i,k2)+co2eml(i,k2))
          tbar(i,3) = 0.5*(tbar(i,2)+tbar(i,1))
          emm(i,3) = emm(i,1)
          tbar(i,4) = tbar(i,3)
          emm(i,4) = emm(i,2)
          o3emm(i,1) = 0.5*(dbvtit(i,k2+1)+dbvtly(i,k2))
          o3emm(i,2) = 0.5*(dbvtit(i,k2)+dbvtly(i,k2))
          o3emm(i,3) = o3emm(i,1)
          o3emm(i,4) = o3emm(i,2)
          temh2o(i,1) = tbar(i,1)
          temh2o(i,2) = tbar(i,2)
          temh2o(i,3) = tbar(i,1)
          temh2o(i,4) = tbar(i,2)
          dpnm(i) = pnm(i,k2+1) - pnm(i,k2)
        end do
!----------------------------------------------------------
!       Weighted Planck functions for trace gases
!
        do wvl = 1 , 14
          do i = 1 , iym1
            bplnk(wvl,i,1) = 0.5*(abplnk1(wvl,i,k2+1)+abplnk2(wvl,i,k2))
            bplnk(wvl,i,2) = 0.5*(abplnk1(wvl,i,k2)+abplnk2(wvl,i,k2))
            bplnk(wvl,i,3) = bplnk(wvl,i,1)
            bplnk(wvl,i,4) = bplnk(wvl,i,2)
          end do
        end do
!---------------------------------------------------------
        do i = 1 , iym1
          rdpnmsq = 1./(pnmsq(i,k2+1)-pnmsq(i,k2))
          rdpnm = 1./dpnm(i)
          p1 = .5*(pbr(i,k2)+pnm(i,k2+1))
          p2 = .5*(pbr(i,k2)+pnm(i,k2))
          uinpl(i,1) = (pnmsq(i,k2+1)-p1**2)*rdpnmsq
          uinpl(i,2) = -(pnmsq(i,k2)-p2**2)*rdpnmsq
          uinpl(i,3) = -(pnmsq(i,k2)-p1**2)*rdpnmsq
          uinpl(i,4) = (pnmsq(i,k2+1)-p2**2)*rdpnmsq
          winpl(i,1) = (.5*(pnm(i,k2+1)-pbr(i,k2)))*rdpnm
          winpl(i,2) = (.5*(-pnm(i,k2)+pbr(i,k2)))*rdpnm
          winpl(i,3) = (.5*(pnm(i,k2+1)+pbr(i,k2))-pnm(i,k2))*rdpnm
          winpl(i,4) = (.5*(-pnm(i,k2)-pbr(i,k2))+pnm(i,k2+1))*rdpnm
          tmp1 = 1./(piln(i,k2+1)-piln(i,k2))
          tmp2 = piln(i,k2+1) - pmln(i,k2)
          tmp3 = piln(i,k2) - pmln(i,k2)
          zinpl(i,1) = (.5*tmp2)*tmp1
          zinpl(i,2) = (-.5*tmp3)*tmp1
          zinpl(i,3) = (.5*tmp2-tmp3)*tmp1
          zinpl(i,4) = (tmp2-.5*tmp3)*tmp1
          pinpl(i,1) = 0.5*(p1+pnm(i,k2+1))
          pinpl(i,2) = 0.5*(p2+pnm(i,k2))
          pinpl(i,3) = 0.5*(p1+pnm(i,k2))
          pinpl(i,4) = 0.5*(p2+pnm(i,k2+1))
        end do

!       FAB AER SAVE uinpl  for aerosl LW forcing calculation
        xuinpl (:,k2,:,jslc) =  uinpl(:,:)
!       FAB AER SAVE uinpl  for aerosl LW forcing calculation

        do kn = 1 , 4
          do i = 1 , iym1
            u(i) = uinpl(i,kn)*dabs(plh2o(i,k2)-plh2o(i,k2+1))
            sqrtu(i) = dsqrt(u(i))
            dw(i) = dabs(w(i,k2)-w(i,k2+1))
            pnew(i) = u(i)/(winpl(i,kn)*dw(i))
            ds2c = dabs(s2c(i,k2)-s2c(i,k2+1))
            uc1(i) = uinpl(i,kn)*ds2c
            uc1(i) = (uc1(i)+1.7E-3*u(i))*(1.+2.*uc1(i))/(1.+15.*uc1(i))
            uc(i) = uinpl(i,kn)*ds2c + 2.E-3*u(i)
          end do
          do i = 1 , iym1
            dtx(i) = temh2o(i,kn) - 250.
            dty(i) = tbar(i,kn) - 250.
            dtyp15(i) = dty(i) + 15.
            dtyp15sq(i) = dtyp15(i)**2
            dtz(i) = dtx(i) - 50.
            dtp(i) = dty(i) - 50.
          end do
          do iband = 2 , 4 , 2
            do i = 1 , iym1
              term1(i,iband) = coefe(1,iband) + coefe(2,iband)*dtx(i)   &
                             & *(1.+c1(iband)*dtx(i))
              term2(i,iband) = coefb(1,iband) + coefb(2,iband)*dtx(i)   &
                             & *(1.+c2(iband)*dtx(i)                    &
                             & *(1.+c3(iband)*dtx(i)))
              term3(i,iband) = coefd(1,iband) + coefd(2,iband)*dtx(i)   &
                             & *(1.+c4(iband)*dtx(i)                    &
                             & *(1.+c5(iband)*dtx(i)))
              term4(i,iband) = coefa(1,iband) + coefa(2,iband)*dty(i)   &
                             & *(1.+c6(iband)*dty(i))
              term5(i,iband) = coefc(1,iband) + coefc(2,iband)*dty(i)   &
                             & *(1.+c7(iband)*dty(i))
            end do
          end do
!
!         abso(i,1)     0 -  800 cm-1   h2o rotation band
!
          do i = 1 , iym1
            a11 = 0.44 + 3.380E-4*dtz(i) - 1.520E-6*dtz(i)*dtz(i)
            a31 = 1.05 - 6.000E-3*dtp(i) + 3.000E-6*dtp(i)*dtp(i)
            a21 = 1.00 + 1.717E-3*dtz(i) - 1.133E-5*dtz(i)*dtz(i)
            a22 = 1.00 + 4.443E-3*dtp(i) + 2.750E-5*dtp(i)*dtp(i)
            a23 = 1.00 + 3.600*sqrtu(i)
            corfac = a31*(a11+((2.*a21*a22)/a23))
            t1t4 = term1(i,2)*term4(i,2)
            t2t5 = term2(i,2)*term5(i,2)
            a = t1t4 + t2t5/(1.+t2t5*sqrtu(i)*corfac)
            fwk(i) = fwcoef + fwc1/(1.+fwc2*u(i))
            fwku(i) = fwk(i)*u(i)
            rsum = dexp(-a*(sqrtu(i)+fwku(i)))
            abso(i,1) = (1.-rsum)*term3(i,2)
!           trab1(i) = rsum
          end do
!
!         abso(i,2)  1200 - 2200 cm-1   h2o vibration-rotation band
!
          do i = 1 , iym1
            a41 = 1.75 - 3.960E-03*dtz(i)
            a51 = 1.00 + 1.3*sqrtu(i)
            a61 = 1.00 + 1.250E-03*dtp(i) + 6.250E-05*dtp(i)*dtp(i)
            corfac = .29*(1.+a41/a51)*a61
            t1t4 = term1(i,4)*term4(i,4)
            t2t5 = term2(i,4)*term5(i,4)
            a = t1t4 + t2t5/(1.+t2t5*sqrtu(i)*corfac)
            rsum = dexp(-a*(sqrtu(i)+fwku(i)))
            abso(i,2) = (1.-rsum)*term3(i,4)
!           trab7(i) = rsum
          end do
!
!         Line transmission in 800-1000 and 1000-1200 cm-1 intervals
!
          do k = 1 , 2
            do i = 1 , iym1
              phi = dexp(a1(k)*dtyp15(i)+a2(k)*dtyp15sq(i))
              psi = dexp(b1(k)*dtyp15(i)+b2(k)*dtyp15sq(i))
              ubar = dw(i)*phi*winpl(i,kn)*1.66*r80257
              pbar = pnew(i)*(psi/phi)
              cf812 = cfa1 + (1.-cfa1)/(1.+ubar*pbar*10.)
              g2 = 1. + ubar*4.0*st(k)*cf812/pbar
              g4 = realk(k)*pbar*r2st(k)*(dsqrt(g2)-1.)
              trline(i,k) = dexp(-g4)
            end do
          end do
          do i = 1 , iym1
            term7(i,1) = coefj(1,1) + coefj(2,1)*dty(i)*(1.+c16*dty(i))
            term8(i,1) = coefk(1,1) + coefk(2,1)*dty(i)*(1.+c17*dty(i))
            term7(i,2) = coefj(1,2) + coefj(2,2)*dty(i)*(1.+c26*dty(i))
            term8(i,2) = coefk(1,2) + coefk(2,2)*dty(i)*(1.+c27*dty(i))
          end do
!
!         abso(i,3)   800 - 1200 cm-1   h2o window
!         abso(i,4)   500 -  800 cm-1   h2o rotation band overlap with
!         co2
          do i = 1 , iym1
            dtym10 = dty(i) - 10.
            denom = 1. + (c30+c31*dtym10*dtym10)*sqrtu(i)
            k21 = term7(i,1) + term8(i,1)/denom
            denom = 1. + (c28+c29*dtym10)*sqrtu(i)
            k22 = term7(i,2) + term8(i,2)/denom
            term9(i,2) = coefi(1,2) + coefi(2,2)*dtx(i)                 &
                       & *(1.+c19*dtx(i)*(1.+c21*dtx(i)                 &
                       & *(1.+c23*dtx(i)*(1.+c25*dtx(i)))))
            tr1 = dexp(-(k21*(sqrtu(i)+fc1*fwku(i))))
            tr2 = dexp(-(k22*(sqrtu(i)+fc1*fwku(i))))
            tr5 = dexp(-((coefh(1,3)+coefh(2,3)*dtx(i))*uc1(i)))
            tr6 = dexp(-((coefh(1,4)+coefh(2,4)*dtx(i))*uc1(i)))
            tr9(i) = tr1*tr5
            tr10(i) = tr2*tr6
            trab2(i) = 0.65*tr9(i) + 0.35*tr10(i)
            th2o(i) = tr10(i)
            trab4(i) = dexp(-(coefg(1,3)+coefg(2,3)*dtx(i))*uc(i))
            trab6(i) = dexp(-(coefg(1,4)+coefg(2,4)*dtx(i))*uc(i))
            term6(i,2) = coeff(1,2) + coeff(2,2)*dtx(i)                 &
                       & *(1.+c9*dtx(i)*(1.+c11*dtx(i)                  &
                       & *(1.+c13*dtx(i)*(1.+c15*dtx(i)))))
            abso(i,3) = term6(i,2)                                      &
                      & *(1.-.5*trab4(i)*trline(i,2)-.5*trab6(i)        &
                      & *trline(i,1))
            abso(i,4) = term9(i,2)*.5*(tr1-tr9(i)+tr2-tr10(i))
          end do
!
!         abso(i,5)  o3  9.6 micrometer (nu3 and nu1 bands)
!
          do i = 1 , iym1
            te = (tbar(i,kn)*r293)**.7
            dplos = dabs(plos(i,k2+1)-plos(i,k2))
            u1 = zinpl(i,kn)*18.29*dplos/te
            u2 = zinpl(i,kn)*.5649*dplos/te
            tlocal = tbar(i,kn)
            tcrfac = dsqrt(tlocal*r250)*te
            beta = r3205*(pinpl(i,kn)*rsslp+dpfo3*tcrfac)
            realnu = te/beta
            tmp1 = u1/dsqrt(4.+u1*(1.+realnu))
            tmp2 = u2/dsqrt(4.+u2*(1.+realnu))
            o3bndi = 74.*te*dlog(1.+tmp1+tmp2)
            abso(i,5) = o3bndi*o3emm(i,kn)*(h2otr(i,k2+1)/h2otr(i,k2))
            to3(i) = 1.0/(1.+0.1*tmp1+0.1*tmp2)
!           trab5(i) = 1.-(o3bndi/(1060-980.))
          end do
!
!         abso(i,6)   co2 15  micrometer band system
!
          do i = 1 , iym1
            dplco2 = plco2(i,k2+1) - plco2(i,k2)
            sqwp = dsqrt(uinpl(i,kn)*dplco2)
            et = dexp(-480./tbar(i,kn))
            sqti(i) = dsqrt(tbar(i,kn))
            rsqti = 1./sqti(i)
            et2 = et*et
            et4 = et2*et2
            omet = (1.-1.5*et2)
            f1co2 = 899.70*omet*(1.+1.94774*et+4.73486*et2)*rsqti
            f1sqwp(i) = f1co2*sqwp
            t1co2(i) = 1./(1.+(245.18*omet*sqwp*rsqti))
            oneme = 1. - et2
            alphat = oneme**3*rsqti
            pi = dabs(dpnm(i))*winpl(i,kn)
            wco2 = 2.5221*co2vmr*pi*rga
            u7(i) = 4.9411E4*alphat*et2*wco2
            u8 = 3.9744E4*alphat*et4*wco2
            u9 = 1.0447E5*alphat*et4*et2*wco2
            u13 = 2.8388E3*alphat*et4*wco2
            tpath = tbar(i,kn)
            tlocal = tbar(i,kn)
            tcrfac = dsqrt((tlocal*r250)*(tpath*r300))
            posqt = (pinpl(i,kn)*rsslp+dpfco2*tcrfac)*rsqti
            rbeta7(i) = 1./(5.3228*posqt)
            rbeta8 = 1./(10.6576*posqt)
            rbeta9 = rbeta7(i)
            rbeta13 = rbeta9
            f2co2(i) = u7(i)/dsqrt(4.+u7(i)*(1.+rbeta7(i)))             &
                     & + u8/dsqrt(4.+u8*(1.+rbeta8))                    &
                     & + u9/dsqrt(4.+u9*(1.+rbeta9))
            f3co2(i) = u13/dsqrt(4.+u13*(1.+rbeta13))
            tmp1 = dlog(1.+f1sqwp(i))
            tmp2 = dlog(1.+f2co2(i))
            tmp3 = dlog(1.+f3co2(i))
            absbnd = (tmp1+2.*t1co2(i)*tmp2+2.*tmp3)*sqti(i)
            abso(i,6) = trab2(i)*emm(i,kn)*absbnd
            tco2(i) = 1.0/(1.0+10.0*u7(i)/dsqrt(4.+u7(i)*(1.+rbeta7(i)))&
                    & )
!           trab3(i) = 1. - bndfct*absbnd
          end do
!
!         Calculate trace gas absorptivity for nearest layer
!
          call trcabn(k2,kn,ucfc11,ucfc12,un2o0,un2o1,uch4,uco211,      &
                    & uco212,uco213,uco221,uco222,uco223,tbar,bplnk,    &
                    & winpl,pinpl,tco2,th2o,to3,uptype,dw,s2c,u,pnew,   &
                    & abstrc,uinpl)
!
!         Total next layer absorptivity:
!
          do i = 1 , iym1
            absnxt(i,k2,kn,jslc) = abso(i,1) + abso(i,2) + abso(i,3)    &
                                 & + abso(i,4) + abso(i,5) + abso(i,6)  &
                                 & + abstrc(i)
          end do
        end do
      end do                    !  end of nearest layer level loop
!
      end subroutine radabs
!
      subroutine radems(s2c,s2t,w,tplnke,plh2o,pnm,plco2,tint,tint4,    &
                      & tlayr,tlayr4,plol,plos,ucfc11,ucfc12,un2o0,     &
                      & un2o1,uch4,uco211,uco212,uco213,uco221,uco222,  &
                      & uco223,uptype,bn2o0,bn2o1,bch4,co2em,co2eml,    &
                      & co2t,h2otr,abplnk1,abplnk2,jslc)
!
!-----------------------------------------------------------------------
!
! Compute emissivity for H2O, CO2, O3
!
! H2O  ....  Uses nonisothermal emissivity for water vapor from
!            Ramanathan, V. and  P.Downey, 1986: A Nonisothermal
!            Emissivity and Absorptivity Formulation for Water Vapor
!            Jouranl of Geophysical Research, vol. 91., D8, pp 8649-8666
!
!
! CO2  ....  Uses absorptance parameterization of the 15 micro-meter
!            (500 - 800 cm-1) band system of Carbon Dioxide, from
!            Kiehl, J.T. and B.P.Briegleb, 1991: A New Parameterization
!            of the Absorptance Due to the 15 micro-meter Band System
!            of Carbon Dioxide Jouranl of Geophysical Research,
!            vol. 96., D5, pp 9013-9019
!
! O3   ....  Uses absorptance parameterization of the 9.6 micro-meter
!            band system of ozone, from Ramanathan, V. and R. Dickinson,
!            1979: The Role of stratospheric ozone in the zonal and
!            seasonal radiative energy balance of the earth-troposphere
!            system. Journal of the Atmospheric Sciences, Vol. 36,
!            pp 1084-1104
!
! Computes individual emissivities, accounting for band overlap, and
! sums to obtain the total.
!
!---------------------------Code history--------------------------------
!
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          J. Kiehl, B. Briegleb, August 1992
!
!-----------------------------------------------------------------------
!
      implicit none
!
!------------------------------Arguments--------------------------------
!
!     Input arguments
!
! s2c     - H2o continuum path length
! s2t     - Tmp and prs wghted h2o path length
! w       - H2o path length
! tplnke  - Layer planck temperature
! plh2o   - H2o prs wghted path length
! pnm     - Model interface pressure
! plco2   - Prs wghted path of co2
! tint    - Model interface temperatures
! tint4   - Tint to the 4th power
! tlayr   - K-1 model layer temperature
! tlayr4  - Tlayr to the 4th power
! plol    - Pressure wghtd ozone path
! plos    - Ozone path
!
!     Trace gas variables
!
! ucfc11  - CFC11 path length
! ucfc12  - CFC12 path length
! un2o0   - N2O path length
! un2o1   - N2O path length (hot band)
! uch4    - CH4 path length
! uco211  - CO2 9.4 micron band path length
! uco212  - CO2 9.4 micron band path length
! uco213  - CO2 9.4 micron band path length
! uco221  - CO2 10.4 micron band path length
! uco222  - CO2 10.4 micron band path length
! uco223  - CO2 10.4 micron band path length
! bn2o0   - pressure factor for n2o
! bn2o1   - pressure factor for n2o
! bch4    - pressure factor for ch4
! uptype  - p-type continuum path length
!
!     Output arguments
!
! co2em   - Layer co2 normalzd plnck funct drvtv
! co2eml  - Intrfc co2 normalzd plnck func drvtv
! co2t    - Tmp and prs weighted path length
! h2otr   - H2o transmission over o3 band
! emplnk  - emissivity Planck factor
! abplnk1 - non-nearest layer Plack factor
! abplnk2 - nearest layer factor
! emstrc  - total trace gas emissivity
!
! Dummy arguments
!
      integer :: jslc
      real(8) , dimension(14,iym1,kzp1) :: abplnk1 , abplnk2
      real(8) , dimension(iym1,kzp1) :: bch4 , bn2o0 , bn2o1 , co2em ,&
           & co2t , h2otr , plco2 , plh2o , plol , plos , pnm , s2c ,   &
           & s2t , tint , tint4 , tlayr , tlayr4 , ucfc11 , ucfc12 ,    &
           & uch4 , uco211 , uco212 , uco213 , uco221 , uco222 ,        &
           & uco223 , un2o0 , un2o1 , uptype , w
      real(8) , dimension(iym1,kz) :: co2eml
      real(8) , dimension(iym1) :: tplnke
      intent (in) jslc , plco2 , plh2o , plol , plos , s2t , tint4 ,    &
                & tlayr4
      intent (out) co2em , co2eml
      intent (inout) co2t , h2otr
!
!---------------------------Local variables-----------------------------
!
! i       - Longitude index
! k       - Level index
! k1      - Level index
! iband   - H2o band index
!
!     Local variables for H2O:
!
! h2oems  - H2o emissivity
! tpathe  - Used to compute h2o emissivity
! a       - Eq(2) in table A3a of R&D
! corfac  - Correction factors in table A3b rotation band absorptivity
! dtp     - Path temperature minus 300 K used in
! dtx     - Planck temperature minus 250 K
! dty     - Path temperature minus 250 K
! dtz     - Planck temperature minus 300 K
! emis    - Total emissivity (h2o+co2+o3)
! rsum    - Eq(1) in table A2 of R&D
! term1   - Equation(5) in table A3a of R&D(1986)
! term2   - Delta a(Te) in table A3a of R&D(1986)
! term3   - B(T) function for rotation and vibration-rotation band emissivity
! term4   - Equation(6) in table A3a of R&D(1986)
! term5   - Delta a(Tp) in table A3a of R&D(1986)
! term6   - B(T) function for window region
! term7   - Kl_inf(i) in eq(8) of table A3a of R&D
! term8   - Delta kl_inf(i) in eq(8)
! term9   - B(T) function for 500-800 cm-1 region
! tr1     - Equation(6) in table A2 for 650-800
! tr2     - Equation(6) in table A2 for 500-650
! tr3     - Equation(4) in table A2 for 650-800
! tr4     - Equation(4),table A2 of R&D for 500-650 
! tr7     - Equation (6) times eq(4) in table A2 of R&D for 650-800 cm-1 region
! tr8     - Equation (6) times eq(4) in table A2 of R&D for 500-650 cm-1 region
! uc      - Y + 0.002U in eq(8) of table A2 of R&D
! pnew    - Effective pressure for h2o linewidth
! trline  - Transmission due to H2O lines in window
! k21     - Exponential coefficient used to calc rot band transmissivity
!           in the 650-800 cm-1 region (tr1)
! k22     - Exponential coefficient used to calc rot band transmissivity
!           in the 500-650 cm-1 region (tr2)
! u       - Pressure weighted H2O path length
! uc1     - H2o continuum pathlength 500-800 cm-1
! r80257  - Conversion factor for h2o pathlength
! a11     - A1 in table A3b for rotation band emiss
! a31     - A3 in table A3b for rotation band emiss
! a21     - First part in numerator of A2 table A3b
! a22     - Second part in numertor of A2 table A3b
! a23     - Denominator of A2 table A3b (rot band)
! t1t4    - Eq(3) in table A3a of R&D
! t2t5    - Eq(4) in table A3a of R&D
! fwk     - Equation(33) in R&D far wing correction
! a41     - Numerator in A2 in Vib-rot (table A3b)
! a51     - Denominator in A2 in Vib-rot(table A3b)
! a61     - A3 factor for Vib-rot band in table A3b
! phi     - Eq(11) in table A3a of R&D
! psi     - Eq(12) in table A3a of R&D
! ubar    - H2o scaled path comment eq(10) table A2
! g1      - Part of eq(10) table A2
! pbar    - H2o scaled pres comment eq(10) table A2
! g3      - Part of eq(10) table A2
! g2      - Part of arguement in eq(10) in table A2
! g4      - Arguement in exp() in eq(10) table A2
! cf812   - Eq(11) in table A2 of R&D
! troco2  - H2o overlap factor for co2 absorption
!
!     Local variables for CO2:
!
! co2ems  - Co2 emissivity
! co2plk  - Used to compute co2 emissivity
! xsum    - Used to calculate path temperature
! t1i     - Co2 hot band temperature factor
! sqti    - Sqrt of temperature
! pi      - Pressure used in co2 mean line width
! et      - Co2 hot band factor
! et2     - Co2 hot band factor
! et4     - Co2 hot band factor
! omet    - Co2 stimulated emission term
! ex      - Part of co2 planck function
! f1co2   - Co2 weak band factor
! f2co2   - Co2 weak band factor
! f3co2   - Co2 weak band factor
! t1co2   - Overlap factor weak bands strong band
! sqwp    - Sqrt of co2 pathlength
! f1sqwp  - Main co2 band factor
! oneme   - Co2 stimulated emission term
! alphat  - Part of the co2 stimulated emiss term
! wco2    - Consts used to define co2 pathlength
! posqt   - Effective pressure for co2 line width
! rbeta7  - Inverse of co2 hot band line width par
! rbeta8  - Inverse of co2 hot band line width par
! rbeta9  - Inverse of co2 hot band line width par
! rbeta13 - Inverse of co2 hot band line width par
! tpath   - Path temp used in co2 band model
! tmp1    - Co2 band factor
! tmp2    - Co2 band factor
! tmp3    - Co2 band factor
! tlayr5  - Temperature factor in co2 Planck func
! rsqti   - Reciprocal of sqrt of temperature
! exm1sq  - Part of co2 Planck function
! u7      - Absorber amount for various co2 band systems
! u8      - Absorber amount for various co2 band systems
! u9      - Absorber amount for various co2 band systems
! u13     - Absorber amount for various co2 band systems
! r250    - Inverse 250K
! r300    - Inverse 300K
! rsslp   - Inverse standard sea-level pressure
!
!     Local variables for O3:
!
! o3ems   - Ozone emissivity
! dbvtt   - Tmp drvtv of planck fctn for tplnke
! te      - Temperature factor
! u1      - Path length factor
! u2      - Path length factor
! phat    - Effecitive path length pressure
! tlocal  - Local planck function temperature
! tcrfac  - Scaled temperature factor
! beta    - Absorption funct factor voigt effect
! realnu  - Absorption function factor
! o3bndi  - Band absorption factor
!
!     Transmission terms for various spectral intervals:
!
! trem4   - H2o   800 - 1000 cm-1
! trem6   - H2o  1000 - 1200 cm-1
! absbnd  - Proportional to co2 band absorptance
! tco2    - co2 overlap factor
! th2o    - h2o overlap factor
! to3     - o3 overlap factor
!
!
! Local variables
!
      real(8) , dimension(iym1) :: a , co2plk , corfac , dbvtt , dtp ,  &
                                  & dtx , dty , dtz , k21 , k22 , pnew ,&
                                  & rsum , tco2 , th2o , to3 , tpathe , &
                                  & tr1 , tr2 , tr3 , tr4 , tr7 , tr8 , &
                                  & trem4 , trem6 , u , uc , uc1 , xsum
      real(8) :: a11 , a21 , a22 , a23 , a31 , a41 , a51 , a61 ,        &
               & absbnd , alphat , beta , cf812 , et , et2 , et4 , ex , &
               & exm1sq , f1co2 , f1sqwp , f2co2 , f3co2 , fwk , g1 ,   &
               & g2 , g3 , g4 , o3bndi , omet , oneme , pbar , phat ,   &
               & phi , pi , posqt , psi , r250 , r300 , r80257 ,        &
               & rbeta13 , rbeta7 , rbeta8 , rbeta9 , realnu , rsqti ,  &
               & sqti , sqwp , t , t1co2 , t1i , t1t4 , t2t5 ,          &
               & tcrfac , te , tlayr5 , tlocal , tmp1 , tmp2 , tmp3 ,   &
               & tpath , u1 , u13 , u2 , u7 , u8 , u9 , ubar , ux , vx ,&
               & wco2
      real(8) , dimension(iym1,kzp1) :: co2ems , emstrc , h2oems ,    &
           & o3ems , troco2
      real(8) :: dbvt , fo3
      real(8) , dimension(iym1,4) :: emis , term1 , term2 , term3 ,    &
                                    & term4 , term5
      real(8) , dimension(14,iym1) :: emplnk
      integer :: i , iband , k , k1
      real(8) , dimension(iym1,2) :: term6 , term7 , term8 , term9 ,   &
                                    & trline
!
!---------------------------Statement functions-------------------------
!
!     Statement functions
!     Derivative of planck function at 9.6 micro-meter wavelength, and
!     an absorption function factor:
!
!
      dbvt(t) = (-2.8911366682D-4+(2.3771251896D-6+1.1305188929D-10*t)  &
              & *t)/(1.D0+(-6.1364820707D-3+1.5550319767D-5*t)*t)
!
      fo3(ux,vx) = ux/dsqrt(4.D0+ux*(1.D0+vx))
!
!-----------------------------------------------------------------------
!
!     Initialize
!
      r80257 = 1.D0/8.0257D-04
!
      r250 = 1.D0/250.D0
      r300 = 1.D0/300.D0
!
!     Planck function for co2
!
      do i = 1 , iym1
        ex = dexp(960.D0/tplnke(i))
        co2plk(i) = 5.D8/((tplnke(i)**4)*(ex-1.))
        co2t(i,1) = tplnke(i)
        xsum(i) = co2t(i,1)*pnm(i,1)
      end do
      k = 1
      do k1 = kzp1 , 2 , -1
        k = k + 1
        do i = 1 , iym1
          xsum(i) = xsum(i) + tlayr(i,k)*(pnm(i,k)-pnm(i,k-1))
          ex = dexp(960./tlayr(i,k1))
          tlayr5 = tlayr(i,k1)*tlayr4(i,k1)
          co2eml(i,k1-1) = 1.2E11*ex/(tlayr5*(ex-1.)**2)
          co2t(i,k) = xsum(i)/pnm(i,k)
        end do
      end do
!     bndfct = 2.d0*22.18/(dsqrt(196.d0)*300.)
!
!     Initialize planck function derivative for O3
!
      do i = 1 , iym1
        dbvtt(i) = dbvt(tplnke(i))
      end do
!
!     Calculate trace gas Planck functions
!
      call trcplk(tint,tlayr,tplnke,emplnk,abplnk1,abplnk2)
!
!     Interface loop
!
      do k1 = 1 , kzp1
!
!       H2O emissivity
!
!       emis(i,1)     0 -  800 cm-1   rotation band
!       emis(i,2)  1200 - 2200 cm-1   vibration-rotation band
!       emis(i,3)   800 - 1200 cm-1   window
!       emis(i,4)   500 -  800 cm-1   rotation band overlap with co2
!
!       For the p type continuum
!
        do i = 1 , iym1
          uc(i) = s2c(i,k1) + 2.E-3*plh2o(i,k1)
          u(i) = plh2o(i,k1)
          pnew(i) = u(i)/w(i,k1)
!
!         Apply scaling factor for 500-800 continuum
!
          uc1(i) = (s2c(i,k1)+1.7E-3*plh2o(i,k1))*(1.+2.*s2c(i,k1))     &
                 & /(1.+15.*s2c(i,k1))
          tpathe(i) = s2t(i,k1)/plh2o(i,k1)
        end do
        do i = 1 , iym1
          dtx(i) = tplnke(i) - 250.
          dty(i) = tpathe(i) - 250.
          dtz(i) = dtx(i) - 50.
          dtp(i) = dty(i) - 50.
        end do
        do iband = 1 , 3 , 2
          do i = 1 , iym1
            term1(i,iband) = coefe(1,iband) + coefe(2,iband)*dtx(i)     &
                           & *(1.+c1(iband)*dtx(i))
            term2(i,iband) = coefb(1,iband) + coefb(2,iband)*dtx(i)     &
                           & *(1.+c2(iband)*dtx(i)*(1.+c3(iband)*dtx(i))&
                           & )
            term3(i,iband) = coefd(1,iband) + coefd(2,iband)*dtx(i)     &
                           & *(1.+c4(iband)*dtx(i)*(1.+c5(iband)*dtx(i))&
                           & )
            term4(i,iband) = coefa(1,iband) + coefa(2,iband)*dty(i)     &
                           & *(1.+c6(iband)*dty(i))
            term5(i,iband) = coefc(1,iband) + coefc(2,iband)*dty(i)     &
                           & *(1.+c7(iband)*dty(i))
          end do
        end do
!
!       emis(i,1)     0 -  800 cm-1   rotation band
!
        do i = 1 , iym1
          a11 = .37 - 3.33E-5*dtz(i) + 3.33E-6*dtz(i)*dtz(i)
          a31 = 1.07 - 1.00E-3*dtp(i) + 1.475E-5*dtp(i)*dtp(i)
          a21 = 1.3870 + 3.80E-3*dtz(i) - 7.8E-6*dtz(i)*dtz(i)
          a22 = 1.0 - 1.21E-3*dtp(i) - 5.33E-6*dtp(i)*dtp(i)
          a23 = 0.9 + 2.62*dsqrt(u(i))
          corfac(i) = a31*(a11+((a21*a22)/a23))
          t1t4 = term1(i,1)*term4(i,1)
          t2t5 = term2(i,1)*term5(i,1)
          a(i) = t1t4 + t2t5/(1.+t2t5*dsqrt(u(i))*corfac(i))
          fwk = fwcoef + fwc1/(1.+fwc2*u(i))
          rsum(i) = dexp(-a(i)*(dsqrt(u(i))+fwk*u(i)))
          emis(i,1) = (1.-rsum(i))*term3(i,1)
!         trem1(i)  = rsum(i)
!
!         emis(i,2)  1200 - 2200 cm-1   vibration-rotation band
!
          a41 = 1.75 - 3.96E-3*dtz(i)
          a51 = 1.00 + 1.3*dsqrt(u(i))
          a61 = 1.00 + 1.25E-3*dtp(i) + 6.25E-5*dtp(i)*dtp(i)
          corfac(i) = .3*(1.+(a41)/(a51))*a61
          t1t4 = term1(i,3)*term4(i,3)
          t2t5 = term2(i,3)*term5(i,3)
          a(i) = t1t4 + t2t5/(1.+t2t5*dsqrt(u(i))*corfac(i))
          fwk = fwcoef + fwc1/(1.+fwc2*u(i))
          rsum(i) = dexp(-a(i)*(dsqrt(u(i))+fwk*u(i)))
          emis(i,2) = (1.-rsum(i))*term3(i,3)
!         trem7(i) = rsum(i)
        end do
!
!       Line transmission in 800-1000 and 1000-1200 cm-1 intervals
!
        do k = 1 , 2
          do i = 1 , iym1
            phi = a1(k)*(dty(i)+15.) + a2(k)*(dty(i)+15.)**2
            psi = b1(k)*(dty(i)+15.) + b2(k)*(dty(i)+15.)**2
            phi = dexp(phi)
            psi = dexp(psi)
            ubar = w(i,k1)*phi
            ubar = (ubar*1.66)*r80257
            pbar = pnew(i)*(psi/phi)
            cf812 = cfa1 + ((1.-cfa1)/(1.+ubar*pbar*10.))
            g1 = (realk(k)*pbar)/(2.*st(k))
            g2 = 1. + (ubar*4.0*st(k)*cf812)/pbar
            g3 = dsqrt(g2) - 1.
            g4 = g1*g3
            trline(i,k) = dexp(-g4)
          end do
        end do
        do i = 1 , iym1
          term7(i,1) = coefj(1,1) + coefj(2,1)*dty(i)*(1.+c16*dty(i))
          term8(i,1) = coefk(1,1) + coefk(2,1)*dty(i)*(1.+c17*dty(i))
          term7(i,2) = coefj(1,2) + coefj(2,2)*dty(i)*(1.+c26*dty(i))
          term8(i,2) = coefk(1,2) + coefk(2,2)*dty(i)*(1.+c27*dty(i))
        end do
!
!       emis(i,3)   800 - 1200 cm-1   window
!
        do i = 1 , iym1
          term6(i,1) = coeff(1,1) + coeff(2,1)*dtx(i)                   &
                     & *(1.+c8*dtx(i)*(1.+c10*dtx(i)                    &
                     & *(1.+c12*dtx(i)*(1.+c14*dtx(i)))))
          trem4(i) = dexp(-(coefg(1,1)+coefg(2,1)*dtx(i))*uc(i))        &
                   & *trline(i,2)
          trem6(i) = dexp(-(coefg(1,2)+coefg(2,2)*dtx(i))*uc(i))        &
                   & *trline(i,1)
          emis(i,3) = term6(i,1)*(1.-.5*trem4(i)-.5*trem6(i))
!
!         emis(i,4)   500 -  800 cm-1   rotation band overlap with co2
!
          k21(i) = term7(i,1) + term8(i,1)                              &
                 & /(1.+(c30+c31*(dty(i)-10.)*(dty(i)-10.))*dsqrt(u(i)))
          k22(i) = term7(i,2) + term8(i,2)                              &
                 & /(1.+(c28+c29*(dty(i)-10.))*dsqrt(u(i)))
          term9(i,1) = coefi(1,1) + coefi(2,1)*dtx(i)                   &
                     & *(1.+c18*dtx(i)*(1.+c20*dtx(i)                   &
                     & *(1.+c22*dtx(i)*(1.+c24*dtx(i)))))
          fwk = fwcoef + fwc1/(1.+fwc2*u(i))
          tr1(i) = dexp(-(k21(i)*(dsqrt(u(i))+fc1*fwk*u(i))))
          tr2(i) = dexp(-(k22(i)*(dsqrt(u(i))+fc1*fwk*u(i))))
          tr3(i) = dexp(-((coefh(1,1)+coefh(2,1)*dtx(i))*uc1(i)))
          tr4(i) = dexp(-((coefh(1,2)+coefh(2,2)*dtx(i))*uc1(i)))
          tr7(i) = tr1(i)*tr3(i)
          tr8(i) = tr2(i)*tr4(i)
          emis(i,4) = term9(i,1)*.5*(tr1(i)-tr7(i)+tr2(i)-tr8(i))
          h2oems(i,k1) = emis(i,1) + emis(i,2) + emis(i,3) + emis(i,4)
          troco2(i,k1) = 0.65*tr7(i) + 0.35*tr8(i)
          th2o(i) = tr8(i)
!         trem2(i)     = troco2(i,k1)
        end do
!
!       CO2 emissivity for 15 micron band system
!
        do i = 1 , iym1
          t1i = dexp(-480./co2t(i,k1))
          sqti = dsqrt(co2t(i,k1))
          rsqti = 1./sqti
          et = t1i
          et2 = et*et
          et4 = et2*et2
          omet = 1. - 1.5*et2
          f1co2 = 899.70*omet*(1.+1.94774*et+4.73486*et2)*rsqti
          sqwp = dsqrt(plco2(i,k1))
          f1sqwp = f1co2*sqwp
          t1co2 = 1./(1.+245.18*omet*sqwp*rsqti)
          oneme = 1. - et2
          alphat = oneme**3*rsqti
          wco2 = 2.5221*co2vmr*pnm(i,k1)*rga
          u7 = 4.9411E4*alphat*et2*wco2
          u8 = 3.9744E4*alphat*et4*wco2
          u9 = 1.0447E5*alphat*et4*et2*wco2
          u13 = 2.8388E3*alphat*et4*wco2
!
          tpath = co2t(i,k1)
          tlocal = tplnke(i)
          tcrfac = dsqrt((tlocal*r250)*(tpath*r300))
          pi = pnm(i,k1)*rsslp + 2.*dpfco2*tcrfac
          posqt = pi/(2.*sqti)
          rbeta7 = 1./(5.3288*posqt)
          rbeta8 = 1./(10.6576*posqt)
          rbeta9 = rbeta7
          rbeta13 = rbeta9
          f2co2 = (u7/dsqrt(4.+u7*(1.+rbeta7)))                         &
                & + (u8/dsqrt(4.+u8*(1.+rbeta8)))                       &
                & + (u9/dsqrt(4.+u9*(1.+rbeta9)))
          f3co2 = u13/dsqrt(4.+u13*(1.+rbeta13))
          tmp1 = dlog(1.+f1sqwp)
          tmp2 = dlog(1.+f2co2)
          tmp3 = dlog(1.+f3co2)
          absbnd = (tmp1+2.*t1co2*tmp2+2.*tmp3)*sqti
          tco2(i) = 1.0/(1.0+10.0*(u7/dsqrt(4.+u7*(1.+rbeta7))))
          co2ems(i,k1) = troco2(i,k1)*absbnd*co2plk(i)
          ex = dexp(960./tint(i,k1))
          exm1sq = (ex-1.)**2
          co2em(i,k1) = 1.2E11*ex/(tint(i,k1)*tint4(i,k1)*exm1sq)
!         trem3(i) = 1. - bndfct*absbnd
        end do
!
!       O3 emissivity
!
        do i = 1 , iym1
          h2otr(i,k1) = dexp(-12.*s2c(i,k1))
          te = (co2t(i,k1)/293.)**.7
          u1 = 18.29*plos(i,k1)/te
          u2 = .5649*plos(i,k1)/te
          phat = plos(i,k1)/plol(i,k1)
          tlocal = tplnke(i)
          tcrfac = dsqrt(tlocal*r250)*te
          beta = (1./.3205D0)*((1./phat)+(dpfo3*tcrfac))
          realnu = (1./beta)*te
          o3bndi = 74.*te*(tplnke(i)/375.)                              &
                 & *dlog(1.+fo3(u1,realnu)+fo3(u2,realnu))
          o3ems(i,k1) = dbvtt(i)*h2otr(i,k1)*o3bndi
          to3(i) = 1.0/(1.+0.1*fo3(u1,realnu)+0.1*fo3(u2,realnu))
!         trem5(i)    = 1.-(o3bndi/(1060-980.))
        end do
!
!       Calculate trace gas emissivities
!
        call trcems(k1,co2t,pnm,ucfc11,ucfc12,un2o0,un2o1,bn2o0,bn2o1,  &
                  & uch4,bch4,uco211,uco212,uco213,uco221,uco222,uco223,&
                  & uptype,w,s2c,u,emplnk,th2o,tco2,to3,emstrc)
!
!       Total emissivity:
!
        do i = 1 , iym1
          emstot(i,k1,jslc) = h2oems(i,k1) + co2ems(i,k1) + o3ems(i,k1) &
                            & + emstrc(i,k1)
        end do
      end do            ! End of interface loop
!
      end subroutine radems
!
      subroutine radoz2(o3vmr,pint,plol,plos)

!-----------------------------------------------------------------------
!
! Computes the path length integrals to the model interfaces given the
! ozone volume mixing ratio
!
!---------------------------Code history--------------------------------
!
! Original version:     CCM1
! Standardized:         J. Rosinski, June 1992
! Reviewed:             J. Kiehl, B. Briegleb, August 1992
! Mixing ratio version: Bruce Biegleb, September 1992
!
!-----------------------------------------------------------------------
!
      implicit none
!
!------------------------------Input arguments--------------------------
!
! o3vmr   - ozone volume mixing ratio
! pint    - Model interface pressures
!
!----------------------------Output arguments---------------------------
!
! plol    - Ozone prs weighted path length (cm)
! plos    - Ozone path length (cm)
!
!
! Dummy arguments
!
      real(8) , dimension(iym1,kz) :: o3vmr
      real(8) , dimension(iym1,kzp1) :: pint , plol , plos
      intent (in) o3vmr , pint
      intent (inout) plol , plos
!
!---------------------------Local workspace-----------------------------
!
! i       - longitude index
! k       - level index
!
!-----------------------------------------------------------------------
!
! Local variables
!
      integer :: i , k
!
!     Evaluate the ozone path length integrals to interfaces;
!     factors of .1 and .01 to convert pressures from cgs to mks:
!
!     Bug fix, 24 May 1996:  the 0.5 and 0.25 factors removed.
!
      do i = 1 , iym1
        plos(i,1) = 0.1*cplos*o3vmr(i,1)*pint(i,1)
        plol(i,1) = 0.01*cplol*o3vmr(i,1)*pint(i,1)*pint(i,1)
      end do
      do k = 2 , kzp1
        do i = 1 , iym1
          plos(i,k) = plos(i,k-1) + 0.1*cplos*o3vmr(i,k-1)              &
                    & *(pint(i,k)-pint(i,k-1))
          plol(i,k) = plol(i,k-1) + 0.01*cplol*o3vmr(i,k-1)             &
                    & *(pint(i,k)*pint(i,k)-pint(i,k-1)*pint(i,k-1))
        end do
      end do
!
      end subroutine radoz2
!
      subroutine radtpl(tnm,ts,qnm,pnm,plh2o,tplnka,s2c,s2t,w,tplnke,   &
                      & tint,tint4,tlayr,tlayr4,pmln,piln)
!-----------------------------------------------------------------------
!
! Compute temperatures and path lengths for longwave radiation
!
!---------------------------Code history--------------------------------
!
! Original version:  CCM1
! Standardized:      L. Buja, June 1992
! Reviewed:          J. Kiehl, B. Briegleb, August 1992
!
!-----------------------------------------------------------------------
!
      implicit none
!
!     Input arguments
!
! tnm    - Model level temperatures
! ts     - Surface skin temperature
! qnm    - Model level specific humidity
! pnm    - Pressure at model interfaces (dynes/cm2)
! plh2o  - Pressure weighted h2o path
!
!     Output arguments
!
! tplnka - Level temperature from interface temperatures
! s2c    - H2o continuum path length
! s2t    - H2o tmp and prs wghtd path length
! w      - H2o path length
! tplnke - Equal to tplnka
! tint   - Layer interface temperature
! tint4  - Tint to the 4th power
! tlayr  - K-1 level temperature
! tlayr4 - Tlayr to the 4th power
! pmln   - Ln(pmidm1)
! piln   - Ln(pintm1)
!
!
! Dummy arguments
!
      real(8) , dimension(iym1,kzp1) :: piln , plh2o , pnm , s2c ,    &
           & s2t , tint , tint4 , tlayr , tlayr4 , tplnka , w
      real(8) , dimension(iym1,kz) :: pmln , qnm , tnm
      real(8) , dimension(iym1) :: tplnke , ts
      intent (in) piln , plh2o , pmln , pnm , qnm , tnm , ts
      intent (out) tint4 , tplnke
      intent (inout) s2c , s2t , tint , tlayr , tlayr4 , tplnka , w
!
!---------------------------Local variables-----------------------------
!
! i      - Longitude index
! k      - Level index
! r296   - Inverse stand temp for h2o continuum
! repsil - Inver ratio mol weight h2o to dry air
! dy     - Thickness of layer for tmp interp
! dpnm   - Pressure thickness of layer
! dpnmsq - Prs squared difference across layer
! rtnm   - Inverse level temperature
!
!-----------------------------------------------------------------------
!
! Local variables
!
      real(8) :: dpnm , dpnmsq , dy , r296 , repsil , rtnm
      integer :: i , k
!
      r296 = 1./296.D0
      repsil = 1./ep2
!
!     Set the top and bottom intermediate level temperatures,
!     top level planck temperature and top layer temp**4.
!
!     Tint is lower interface temperature
!     (not available for bottom layer, so use ground temperature)
!
      do i = 1 , iym1
        tint(i,kzp1) = ts(i)
        tint4(i,kzp1) = tint(i,kz + 1)**4
        tplnka(i,1) = tnm(i,1)
        tint(i,1) = tplnka(i,1)
        tlayr4(i,1) = tplnka(i,1)**4
        tint4(i,1) = tlayr4(i,1)
      end do
!
!     Intermediate level temperatures are computed using temperature
!     at the full level below less dy*delta t,between the full level
!
      do k = 2 , kz
        do i = 1 , iym1
          dy = (piln(i,k)-pmln(i,k))/(pmln(i,k-1)-pmln(i,k))
          tint(i,k) = tnm(i,k) - dy*(tnm(i,k)-tnm(i,k-1))
          tint4(i,k) = tint(i,k)**4
        end do
      end do
!
!     Now set the layer temp=full level temperatures and establish a
!     planck temperature for absorption (tplnka) which is the average
!     the intermediate level temperatures.  Note that tplnka is not
!     equal to the full level temperatures.
!
      do k = 2 , kzp1
        do i = 1 , iym1
          tlayr(i,k) = tnm(i,k-1)
          tlayr4(i,k) = tlayr(i,k)**4
          tplnka(i,k) = .5*(tint(i,k)+tint(i,k-1))
        end do
      end do
!
!     Calculate tplank for emissivity calculation.
!     Assume isothermal tplnke i.e. all levels=ttop.
!
      do i = 1 , iym1
        tplnke(i) = tplnka(i,1)
        tlayr(i,1) = tint(i,1)
      end do
!
!     Now compute h2o path fields:
!
      do i = 1 , iym1
        s2t(i,1) = plh2o(i,1)*tnm(i,1)
!       ccm3.2
!       w(i,1)   = (plh2o(i,1)*2.) / pnm(i,1)
!       s2c(i,1) = plh2o(i,1) * qnm(i,1) * repsil
 
!       ccm3.6.6
        w(i,1) = sslp*(plh2o(i,1)*2.)/pnm(i,1)
        rtnm = 1./tnm(i,1)
        s2c(i,1) = plh2o(i,1)*exp(1800.*(rtnm-r296))*qnm(i,1)*repsil
      end do
      do k = 1 , kz
        do i = 1 , iym1
          dpnm = pnm(i,k+1) - pnm(i,k)
          dpnmsq = pnm(i,k+1)**2 - pnm(i,k)**2
          rtnm = 1./tnm(i,k)
          s2t(i,k+1) = s2t(i,k) + rgsslp*dpnmsq*qnm(i,k)*tnm(i,k)
          w(i,k+1) = w(i,k) + rga*qnm(i,k)*dpnm
          s2c(i,k+1) = s2c(i,k) + rgsslp*dpnmsq*qnm(i,k)                &
                     & *dexp(1800.*(rtnm-r296))*qnm(i,k)*repsil
        end do
      end do
!
      end subroutine radtpl
!
      subroutine radinp(pmid,pint,h2ommr,cld,o3vmr,pmidrd,pintrd,plco2, &
                      & plh2o,tclrsf,eccf,o3mmr)

!-----------------------------------------------------------------------
!
! Set latitude and time dependent arrays for input to solar
! and longwave radiation.
!
! Convert model pressures to cgs, compute path length arrays needed for the
! longwave radiation, and compute ozone mixing ratio, needed for the solar
! radiation.
!
!---------------------------Code history--------------------------------
!
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          J. Kiehl, B. Briegleb, August 1992
!
!-----------------------------------------------------------------------
!
      implicit none
!
! Dummy arguments
!
      real(8) :: eccf
      real(8) , dimension(iym1,kzp1) :: cld , pint , pintrd , plco2 ,   &
           & plh2o , tclrsf
      real(8) , dimension(iym1,kz) :: h2ommr , o3mmr , o3vmr , pmid ,   &
           & pmidrd
      intent (in) cld , h2ommr , o3vmr , pint , pmid
      intent (out) eccf , o3mmr , plco2 , pmidrd
      intent (inout) pintrd , plh2o , tclrsf
!
! Local variables
!
#ifndef CLM
      real(8) :: theta
#endif
      real(8) :: cpwpl , vmmr
      integer :: i , k
!
!------------------------------Arguments--------------------------------
!
!     Input arguments
!
! pmid    - Pressure at model mid-levels (pascals)
! pint    - Pressure at model interfaces (pascals)
! h2ommr  - H2o mass mixing ratio
! cld     - Fractional cloud cover
! o3vmr   - ozone volume mixing ratio
!
!     Output arguments
!
! pmidrd  - Pressure at mid-levels (dynes/cm*2)
! pintrd  - Pressure at interfaces (dynes/cm*2)
! plco2   - Vert. pth lngth of co2 (prs-weighted)
! plh2o   - Vert. pth lngth h2o vap.(prs-weighted)
! tclrsf  - Product of clr-sky fractions from top of atmosphere to level.
! eccf    - Earth-sun distance factor
! o3mmr   - Ozone mass mixing ratio
!
!---------------------------Local variables-----------------------------
!
! i       - Longitude loop index
! k       - Vertical loop index
! theta   - Earth orbit seasonal angle in radians
! cpwpl   - Const in co2 mixing ratio to path length conversn
! vmmr    - Ozone volume mixing ratio
!
!-----------------------------------------------------------------------
!
!     Compute solar distance factor and cosine solar zenith angle usi
!     day value where a round day (such as 213.0) refers to 0z at
!     Greenwich longitude.
!
!     Use formulas from Paltridge, G.W. and C.M.R. Platt 1976: Radiative
!     Processes in Meterology and Climatology, Elsevier Scientific
!     Publishing Company, New York  p. 57, p. 62,63.
!
!     Compute eccentricity factor (sun-earth distance factor)
!
#ifdef CLM
      eccf  = r2ceccf
#else
      theta = twopi*calday/dayspy
      eccf = 1.000110 + .034221*dcos(theta) + .001280*dsin(theta)       &
           & + .000719*dcos(2.*theta) + .000077*dsin(2.*theta)
#endif
!
!     Convert pressure from pascals to dynes/cm2
!
      do k = 1 , kz
        do i = 1 , iym1
          pmidrd(i,k) = pmid(i,k)*10.0
          pintrd(i,k) = pint(i,k)*10.0
        end do
      end do
      do i = 1 , iym1
        pintrd(i,kzp1) = pint(i,kz + 1)*10.0
      end do
!
!     Compute path quantities used in the longwave radiation:
!
      vmmr = amco2/amd
      cpwpl = vmmr*0.5/(gtigts*sslp)
      do i = 1 , iym1
        plh2o(i,1) = rgsslp*h2ommr(i,1)*pintrd(i,1)*pintrd(i,1)
        plco2(i,1) = co2vmr*cpwpl*pintrd(i,1)*pintrd(i,1)
        tclrsf(i,1) = 1.
      end do
      do k = 1 , kz
        do i = 1 , iym1
          plh2o(i,k+1) = plh2o(i,k)                                     &
                       & + rgsslp*(pintrd(i,k+1)**2-pintrd(i,k)**2)     &
                       & *h2ommr(i,k)
          plco2(i,k+1) = co2vmr*cpwpl*pintrd(i,k+1)**2
          tclrsf(i,k+1) = tclrsf(i,k)*(1.-cld(i,k+1))
        end do
      end do
!
!     Convert ozone volume mixing ratio to mass mixing ratio:
!
      vmmr = amo/amd
      do k = 1 , kz
        do i = 1 , iym1
          o3mmr(i,k) = vmmr*o3vmr(i,k)
        end do
      end do
!
      end subroutine radinp
!
      function isrchfgt(n,array,inc,rtarg)
      implicit none
!
! Dummy arguments
!
      integer :: inc , n
      real(8) :: rtarg
      real(8) , dimension(*) :: array
      integer :: isrchfgt
      intent (in) array , inc , n , rtarg
!
! Local variables
!
      integer :: i , ind
!
      if ( n.le.0 ) then
        isrchfgt = 0
        return
      end if
      ind = 1
      do i = 1 , n
        if ( array(ind).gt.rtarg ) then
          isrchfgt = i
          return
        else
          ind = ind + inc
        end if
      end do
      isrchfgt = ind
      end function isrchfgt
 
      function isrchfle(n,array,inc,rtarg)
      implicit none
!
! Dummy arguments
!
      integer :: inc , n
      real(8) :: rtarg
      real(8) , dimension(*) :: array
      integer :: isrchfle
      intent (in) array , inc , n , rtarg
!
! Local variables
!
      integer :: i , ind
!
      if ( n.le.0 ) then
        isrchfle = 0
        return
      end if
      ind = 1
      do i = 1 , n
        if ( array(ind).le.rtarg ) then
          isrchfle = i
          return
        else
          ind = ind + inc
        end if
      end do
      isrchfle = ind
      end function isrchfle
!
      subroutine wheneq(n,array,inc,itarg,indx,nval)
 
      implicit none
!
! Dummy arguments
!
      integer :: inc , itarg , n , nval
      integer , dimension(*) :: array , indx
      intent (in) array , inc , itarg , n
      intent (out) indx
      intent (inout) nval
!
! Local variables
!
      integer :: i , ina
!
      ina = 1
      nval = 0
      if ( inc.lt.0 ) ina = (-inc)*(n-1) + 1
      do i = 1 , n
        if ( array(ina).eq.itarg ) then
          nval = nval + 1
          indx(nval) = i
        end if
        ina = ina + inc
      end do
 
      end subroutine wheneq
!
      subroutine whenne(n,array,inc,itarg,indx,nval)
 
      implicit none
!
! Dummy arguments
!
      integer :: inc , itarg , n , nval
      integer , dimension(*) :: array , indx
      intent (in) array , inc , itarg , n
      intent (out) indx
      intent (inout) nval
!
! Local variables
!
      integer :: i , ina
!
      ina = 1
      nval = 0
      if ( inc.lt.0 ) ina = (-inc)*(n-1) + 1
      do i = 1 , n
        if ( array(ina).ne.itarg ) then
          nval = nval + 1
          indx(nval) = i
        end if
        ina = ina + inc
      end do
      end subroutine whenne
!
      subroutine whenfgt(n,array,inc,rtarg,indx,nval)
 
      implicit none
!
! Dummy arguments
!
      integer :: inc , n , nval
      real(8) :: rtarg
      real(8) , dimension(*) :: array
      integer , dimension(*) :: indx
      intent (in) array , inc , n , rtarg
      intent (out) indx
      intent (inout) nval
!
! Local variables
!
      integer :: i , ina
!
      ina = 1
      nval = 0
      if ( inc.lt.0 ) ina = (-inc)*(n-1) + 1
      do i = 1 , n
        if ( array(ina).gt.rtarg ) then
          nval = nval + 1
          indx(nval) = i
        end if
        ina = ina + inc
      end do
 
      end subroutine whenfgt
!
      subroutine whenflt(n,array,inc,rtarg,indx,nval)
 
      implicit none
!
! Dummy arguments
!
      integer :: inc , n , nval
      real(8) :: rtarg
      real(8) , dimension(*) :: array
      integer , dimension(*) :: indx
      intent (in) array , inc , n , rtarg
      intent (out) indx
      intent (inout) nval
!
! Local variables
!
      integer :: i , ina
!
      ina = 1
      nval = 0
      if ( inc.lt.0 ) ina = (-inc)*(n-1) + 1
      do i = 1 , n
        if ( array(ina).lt.rtarg ) then
          nval = nval + 1
          indx(nval) = i
        end if
        ina = ina + inc
      end do
      end subroutine whenflt
!
      function intmax(n,iy,inc)
!
      implicit none
!
! Dummy arguments
!
      integer :: inc , n
      integer :: intmax
      integer , dimension(*) :: iy
      intent (in) inc , iy , n
!
! Local variables
!
      integer :: i , mx
!
      mx = iy(1)
      intmax = 1
      do i = 1 + inc , inc*n , inc
        if ( iy(i).gt.mx ) then
          mx = iy(i)
          intmax = i
        end if
      end do
      end function intmax
!
      end module mod_radiation

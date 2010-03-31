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
 
      subroutine radcsw(pint,h2ommr,o3mmr,cld,clwp,rel,rei,fice,eccf,   &
                      & asdir,asdif,aldir,aldif,solin,qrs,fsns,fsnt,    &
                      & fsds,fsnsc,fsntc,sols,soll,solsd,solld,fsnirt,  &
                      & fsnrtc,fsnirtsq,tauxar_mix,tauasc_mix,gtota_mix,&
                      & ftota_mix,tauxar_mix_cs,tauasc_mix_cs,          &
                      & gtota_mix_cs,ftota_mix_cs,aeradfo,aeradfos)
 
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
      use mod_regcm_param
      use mod_param2
      use mod_bats
      use mod_tracer
      use mod_aerosol , only : nspi
      use mod_constants , only : gocp , gtigts , sslp , rga
      implicit none
!
! PARAMETER definitions
!
! scon        - Solar constant (in erg/cm**2/sec)
! nspint      - Num of spctrl intervals across solar spectrum
! v_raytau_xx - Constants for new bands
! v_abo3_xx   - Constants for new bands
!
      real(8) , parameter :: scon = 1.367E6
      integer , parameter :: nspint = 19
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
      real(8) , dimension(iym1,0:kz,nspi) :: ftota_mix , gtota_mix , &
           & tauasc_mix , tauxar_mix
      real(8) , dimension(iym1,nspi) :: ftota_mix_cs , gtota_mix_cs ,  &
           & tauasc_mix_cs , tauxar_mix_cs
      intent (in) aldif , aldir , asdif , asdir , cld , clwp , eccf ,   &
                & fice , ftota_mix , gtota_mix , h2ommr , o3mmr , pint ,&
                & rei , rel , tauasc_mix , tauxar_mix
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
! ptop     - Lower interface pressure of extra layer
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
               & pthco2 , pthh2o , ptho2 , ptho3 , ptop , rdenom ,      &
               & sqrco2 , tmp1 , tmp1i , tmp1l , tmp2 , tmp2i , tmp2l , &
               & tmp3i , tmp3l , trayoslp , wavmid , wgtint
      real(8) , dimension(nspint) :: abco2 , abh2o , abo2 , abo3 ,      &
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
      integer, external :: isrchfgt , isrchfle
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
      do i = 1 , iym1
        fsds(i) = 0.0
        fsnirt(i) = 0.0
        fsnrtc(i) = 0.0
        fsnirtsq(i) = 0.0
        fsnt(i) = 0.0
        fsns(i) = 0.0
        solin(i) = 0.0
        fsnsc(i) = 0.0
        fsntc(i) = 0.0
        sols(i) = 0.0
        soll(i) = 0.0
        solsd(i) = 0.0
        solld(i) = 0.0
        sabveg(i) = 0.0
        solis(i) = 0.0
        solvs(i) = 0.0
        solvd(i) = 0.0
!
        aeradfo(i) = 0.0
        aeradfos(i) = 0.0
        x0fsntc(i) = 0.0
        x0fsnsc(i) = 0.0
        x0fsnrtc(i) = 0.0
      end do
      do k = 1 , kz
        do i = 1 , iym1
          qrs(i,k) = 0.0
        end do
      end do
      do k = 1 , kz
        do i = 1 , iym1
          pdel = pint(i,k+1) - pint(i,k)
          path = pdel/gtigts
        end do
      end do
!
!     Compute starting, ending daytime loop indices:
!
      nloop = 0
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
          ptop = pflx(i,1)
          ptho2 = o2mmr*ptop/gtigts
          ptho3 = o3mmr(i,1)*ptop/gtigts
          pthco2 = sqrco2*(ptop/gtigts)
          h2ostr = sqrt(1./h2ommr(i,1))
          zenfac(i) = sqrt(coszrs(i))
          pthh2o = ptop**2*tmp1 + (ptop*rga)*(h2ostr*zenfac(i)*delta)
          uh2o(i,0) = h2ommr(i,1)*pthh2o
          uco2(i,0) = zenfac(i)*pthco2
          uo2(i,0) = zenfac(i)*ptho2
          uo3(i,0) = ptho3
        end do
      end do
!
      tmp2 = delta/gtigts
      do k = 1 , kz
        do n = 1 , nloop
          do i = is(n) , ie(n)
            pdel = pflx(i,k+1) - pflx(i,k)
            path = pdel/gtigts
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
      do ns = 1 , nspint
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
        else
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
              fluxup(i,k) = (exptdn(i,k)*rupdir(i,k)+(tottrn(i,k)-exptdn&
                          & (i,k))*rupdif(i,k))*rdenom
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
!       CLEAR SKY CALCULATION PLUS AUTOMATCI CALCULATEION OF AEROSOL
!       FORCING RAD CLR is called 2 times , one with O aerosol OP , and
!       one with actual aerosol. DIFFERENCE  in net TOA SW for the two
!       case is saved as one more variable in the rad file. The
!       outputed TOASW ( fsntc, clrst) is accounting for aerosol.
        if ( idirect.ge.1 ) then
 
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
!         if(ns==8)  print *,'DANS radscw FUPav',
!         &    maxval(fluxup(:,0)* solflx(:)*1.E-3)
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
 
!       if(ns==8)  print *,'DANS radscw FUPav',
!       &    maxval(fluxup(:,0)* solflx(:)*1.E-3)
 
!
!       End of clear sky calculation
!
      end do                    ! End of spectral interval loop
 
!     FAB calculation of TOA aerosol radiative forcing
      if ( idirect.ge.1 ) then
        do n = 1 , nloop
          do i = is(n) , ie(n)
!           test
            aeradfo(i) = -(x0fsntc(i)-fsntc(i))
            aeradfos(i) = -(x0fsnsc(i)-fsnsc(i))
          end do
        end do
!TEST
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

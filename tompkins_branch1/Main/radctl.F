cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine radctl(jslc    ,clat    ,coslat  ,lat  ,ts   ,
     $                  pmid    ,pint    ,pmln    ,piln    ,t       ,
     $                  h2ommr  ,cld     ,effcld  ,clwp    ,coszrs  ,
     $                  albs    ,albsd   ,albl    ,albld   ,fsns    ,
     $                  qrs     ,qrl     ,flwds   ,lwup    ,rel     ,
     $                  rei     ,fice    ,sols    ,soll    ,solsd   ,
     $                  solld   ,emiss1d ,
c++csz
     $     fsnt,fsntc,fsnsc,flnt,flns,flntc,flnsc,solin, ! output
Cadd by bixq
     &     alb,albc,fsds,fsnirt,fsnrtc,fsnirtsq,
Cadd_
     $     eccf, ! output
     $     o3vmr) ! input
C    

c
C-----------------------------------------------------------------------
C
C Driver for radiation computation.
C
C Radiation uses cgs units, so conversions must be done from
C model fields to radiation fields.
C
C Calling sequence:
C
C     radinp      Converts units of model fields and computes ozone
C                 mixing ratio for solar scheme
C
C     radcsw      Performs solar computation
C       radalb    Computes surface albedos
C       radded    Computes delta-Eddington solution
C       radclr    Computes diagnostic clear sky fluxes
C
C     radclw      Performs longwave computation
C
C       radtpl    Computes path quantities
C       radems    Computes emissivity
C       radabs    Computes absorptivity
C
C---------------------------Code history--------------------------------
C
C Original version:  CCM1
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Kiehl, B. Briegleb, August 1992
C
C Modified:          B. Briegleb, March 1995 to add aerosol
C                    to shortwave code
C
C-----------------------------------------------------------------------
c
c $Id: crm.F,v 1.9 1995/11/10 17:39:12 zender Exp $
c $Author: zender $
c
c
c $Id: implicit.h,v 1.1.1.1 1995/02/09 23:26:52 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
      integer jslc
c
C------------------------------Parameters-------------------------------
c
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
c
      include 'Commons/parrad.cb'

C------------------------------Commons----------------------------------
      include 'Commons/comctl.cb'
      include 'Commons/comhst.cb'
      include 'Commons/comtim.cb'
      include 'Commons/aeroppt.cb'
C-----------------------------------------------------------------------
C
      real*8  aermmr        ! Radiation level aerosol mass mixing ratio
      common/aermr/ aermmr(plond,plevr)
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      integer lat               ! Latitude row index
      real*8  ts(plond),           ! Surface (skin) temperature
     $     pmid(plond,plev),    ! Model level pressures
     $     pint(plond,plevp),   ! Model interface pressures
     $     pmln(plond,plev),    ! Natural log of pmid
     $     rel(plond,plev),     ! liquid cloud particle effective radius
     $     rei(plond,plev),     ! ice effective drop size (microns)
     $     fice(plond,plev),    ! fractional ice content within cloud
     $     piln(plond,plevp),   ! Natural log of pint
     $     t(plond,plev),       ! Model level temperatures
     $     h2ommr(plond,plev),  ! Model level specific humidity
     $     cld(plond,plevp),    ! Fractional cloud cover
     $     effcld(plond,plevp), ! Effective fractional cloud cover
     $     clwp(plond,plev)     ! Cloud liquid water path
      real*8  coszrs(plond),       ! Cosine solar zenith angle
     $     albs(plond),
     $     albsd(plond),
     $     albl(plond),
     $     albld(plond)
      real*8  emiss1d(plond)
      real*8  clat(plond),                  ! current latitude(radians)
     $     coslat(plond)                 ! cosine latitude
C
C Output solar arguments
C
      real*8  fsns(plond),         ! Surface absorbed solar flux
     $     sols(plond),         ! Downward solar rad onto surface (sw direct)
     $     soll(plond),         ! Downward solar rad onto surface (lw direct)
     $     solsd(plond),        ! Downward solar rad onto surface (sw diffuse)
     $     solld(plond),        ! Downward solar rad onto surface (lw diffuse)
     $     qrs(plond,plev)      ! Solar heating rate
C
C Output longwave arguments
C
      real*8  qrl(plond,plev),     ! Longwave cooling rate
     $     flwds(plond),        ! Surface down longwave flux
     $     lwup(plond)          ! Surface up longwave flux from coupler
C
C---------------------------Local variables-----------------------------
C
      integer i                 ! index
      real*8  solin(plond),        ! Solar incident flux
     $     fsnt(plond),         ! Net column abs solar flux at model top
     $     fsntc(plond),        ! Clear sky total column abs solar flux
     $     fsnsc(plond)         ! Clear sky surface abs solar flux

CEES next 3 added, they are calculated in radcsw
      real*8 fsnirt(plond)      ! Near-IR flux absorbed at toa
      real*8 fsnrtc(plond)      ! Clear sky near-IR flux absorbed at toa
      real*8 fsnirtsq(plond)    ! Near-IR flux absorbed at toa >= 0.7 microns
      real*8 fsds(plond)        ! Flux Shortwave Downwelling Surface
      real*8 alb(plond)
      real*8 albc(plond)

      real*8  flnt(plond),         ! Net outgoing lw flux at model top
     $     flns(plond),         ! Srf longwave cooling (up-down) flux
     $     flntc(plond),        ! Clear sky lw flux at model top
     $     flnsc(plond)         ! Clear sky lw flux at srf (up-down)
      real*8  pbr(plond,plevr),    ! Model mid-level pressures (dynes/cm2)
     $     pnm(plond,plevrp),   ! Model interface pressures (dynes/cm2)
     $     o3vmr(plond,plevr),  ! Ozone volume mixing ratio
     $     o3mmr(plond,plevr),  ! Ozone mass mixing ratio
     $     plco2(plond,plevrp), ! Prs weighted CO2 path
     $     plh2o(plond,plevrp), ! Prs weighted H2O path
     $     tclrsf(plond,plevrp),! Total clear sky fraction, level to space
     $     eccf                 ! Earth/sun distance factor
      real*8  n2o(plond,plev),        ! nitrous oxide mass mixing ratio
     $     ch4(plond,plev),        ! methane mass mixing ratio
     $     cfc11(plond,plev),      ! cfc11 mass mixing ratio
     $     cfc12(plond,plev)       ! cfc12 mass mixing ratio
      real*8  rh(plond,plevr)      ! level relative humidity (fraction)
C


C  aerosol mixing optical properties output 
  
      real*8  tauxar_mix(plond,0:plev,nspi)  ! aerosol extinction optical depth
      real*8  tauasc_mix(plond,0:plev,nspi)
      real*8  gtota_mix(plond,0:plev,nspi)
      real*8  ftota_mix(plond,0:plev,nspi)
      real*8  tauxar_mix_cs(plond,nspi)  ! aerosol extinction optical depth
      real*8  tauasc_mix_cs(plond,nspi)
      real*8  gtota_mix_cs(plond,nspi)
      real*8  ftota_mix_cs(plond,nspi)
      real*8  aeradfo (plond),  aeradfos (plond) 


C Declare local arrays to which model input arrays are interpolated here.
C Current default is none since radiation grid = model grid.
C
C Externals.
C
c      external radinp,          ! Computes latitude dependent radiation input
c     $         aermix,          ! Specifies aerosol mass mixing ratio
c     $         radcsw,          ! Computes solar radiation
c     $         radozn,          ! Computes ozone volume mixing ratio
c     $         radclw,          ! Computes longwave radiation
c     $         torgrid,         ! Interpolate model variables to radiation grid
c     $         fmrgrid          ! Interpolate radiation variables to model grid
C--------------------------------------------------------------------------
C
C Interpolate ozone volume mixing ratio to model levels
C
c++csz
c
c     This is a do nothing routine in the CRM.
c     Instead of interpolating the o3vmr from the time-interpolated
c     values, we pass compute o3vmr in getdat() and pass it directly
c     into radctl(). o3mmr will be computed in radinp().
c
c      call radozn(lat     ,pmid    ,o3vmr   )
c--csz
C
C Set latitude dependent radiation input
C
      call radinp(pmid    ,pint    ,h2ommr  ,cld     ,o3vmr   ,
     $            pbr     ,pnm     ,plco2   ,plh2o   ,tclrsf  ,
     $            eccf    ,o3mmr   )
C
C Solar radiation computation
C
      if (dosw) then
C
C Specify aerosol mass mixing ratio
C
         call aermix(pnm     ,h2ommr  ,rh ,jslc)

         call aeroppt( rh,pint,
     &        tauxar_mix, tauasc_mix, gtota_mix, ftota_mix,
     &        tauxar_mix_cs, tauasc_mix_cs, gtota_mix_cs, ftota_mix_cs)  

C



         call radcsw(pnm     ,h2ommr  ,o3mmr   ,rh      ,
     $               cld     ,clwp    ,rel     ,rei     ,fice    ,
     $               eccf             ,albs    ,albsd   ,albl    ,
     $               albld   ,solin   ,qrs     ,fsns    ,fsnt    ,
     $               fsds    ,fsnsc   ,fsntc   ,sols    ,soll    ,
     $               solsd   ,solld   ,fsnirt  ,fsnrtc  ,fsnirtsq,
     $               tauxar_mix,  tauasc_mix, gtota_mix,ftota_mix,
     $         tauxar_mix_cs, tauasc_mix_cs, gtota_mix_cs,ftota_mix_cs,
     $                aeradfo, aeradfos)


         call aerout(jslc,tauxar_mix,tauasc_mix,gtota_mix,aeradfo,
     &               aeradfos)  


C Convert units of shortwave fields needed by rest of model from CGS to MKS
C
         do i=1,plon
            solin(i) = solin(i)*1.e-3
            fsnt(i)  = fsnt(i) *1.e-3
            fsns(i)  = fsns(i) *1.e-3
            fsntc(i) = fsntc(i)*1.e-3
            fsnsc(i) = fsnsc(i)*1.e-3

            fsds(i)  = fsds(i)*1.e-3
            fsnirt(i)= fsnirt(i)*1.e-3
            fsnrtc(i)= fsnrtc(i)*1.e-3
            fsnirtsq(i)= fsnirtsq(i)*1.e-3
         end do
C
C Calculate/outfld albedo and clear sky albedo
C
          do i=1,plon
             if (solin(i).gt.0.) then
                alb(i) = (solin(i) - fsnt(i)) / solin(i)
             else
                alb(i) = 0.
             end if
          end do
C
          do i=1,plon
             if (solin(i).gt.0.) then
                albc(i) = (solin(i) - fsntc(i)) / solin(i)
             else
                albc(i) = 0.
             end if
          end do
      end if
C
C Longwave radiation computation
C
      if (dolw) then
c
c Specify trace gas mixing ratios
c
         call trcmix(pmid, clat, coslat, n2o, ch4, cfc11, cfc12)
c
         call radclw(jslc, lat  ,ts   ,t    ,h2ommr  ,o3vmr   ,
     $               pbr     ,pnm     ,pmln    ,piln    ,plco2   ,
     $               plh2o   ,n2o     ,ch4     ,cfc11   ,cfc12   ,
     $               effcld  ,tclrsf  ,qrl     ,flns    ,flnt    ,
     $               flnsc   ,flntc   ,flwds   ,lwup    ,emiss1d )
C
C Convert units of longwave fields needed by rest of model from CGS to MKS
C
         do i=1,plon
            flnt(i)  = flnt(i)*1.e-3
            flns(i)  = flns(i)*1.e-3
            flntc(i) = flntc(i)*1.e-3
            flnsc(i) = flnsc(i)*1.e-3
            flwds(i) = flwds(i)*1.e-3
         end do
      end if

      return
      end

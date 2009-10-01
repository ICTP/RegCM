
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine radinp(pmid    ,pint    ,h2ommr  ,cld     ,o3vmr   ,
     $                  pmidrd  ,pintrd  ,plco2   ,plh2o   ,tclrsf  ,
     $                  eccf    ,o3mmr   )
C-----------------------------------------------------------------------
C
C Set latitude and time dependent arrays for input to solar
C and longwave radiation.
C
C Convert model pressures to cgs, compute path length arrays needed for the
C longwave radiation, and compute ozone mixing ratio, needed for the solar
C radiation.
C
C---------------------------Code history--------------------------------
C
C Original version:  CCM1
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Kiehl, B. Briegleb, August 1992
C
C-----------------------------------------------------------------------
c
c $Id: radinp.F,v 1.1.1.1 1995/02/09 23:27:02 ccm2 Exp $
c $Author: ccm2 $
c
c
c $Id: implicit.h,v 1.1.1.1 1995/02/09 23:26:52 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
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
C
C------------------------------Commons----------------------------------
      include 'Commons/comtim.cb'
      include 'Commons/crdcon.cb'
C
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real*8  pmid(plond,plev),   ! Pressure at model mid-levels (pascals)
     $     pint(plond,plevp),  ! Pressure at model interfaces (pascals)
     $     h2ommr(plond,plev), ! H2o mass mixing ratio
     $     cld(plond,plevp),   ! Fractional cloud cover
     $     o3vmr(plond,plev)   ! ozone volume mixing ratio
C
C Output arguments
C
      real*8  pmidrd(plond,plev), ! Pressure at mid-levels (dynes/cm*2)
     $     pintrd(plond,plevp),! Pressure at interfaces (dynes/cm*2)
     $     plco2(plond,plevp), ! Vert. pth lngth of co2 (prs-weighted)
     $     plh2o(plond,plevp), ! Vert. pth lngth h2o vap.(prs-weighted)
     $     tclrsf(plond,plevp) ! Product of clr-sky fractions from top
C                              ! of atmosphere to level.
      real*8  eccf,               ! Earth-sun distance factor
     $     o3mmr(plond,plev)   ! Ozone mass mixing ratio
C
C---------------------------Local variables-----------------------------
C
      integer i,    ! Longitude loop index
     $        k     ! Vertical loop index
      real*8  theta,   ! Earth orbit seasonal angle in radians
     $     p0 ,     ! Standard pressure (dynes/cm**2)
     $     amd,     ! Effective molecular weight of dry air (g/mol)
     $     amo,     ! Molecular weight of ozone (g/mol)
     $     amco2,   ! Molecular weight of co2   (g/mol)
     $     cpwpl,   ! Const in co2 mixing ratio to path length conversn
     $     vmmr     ! Ozone volume mixing ratio
      save     p0   ,amd   ,amo  ,amco2
C
      data p0    /  1.01325e6 /
      data amd   /  28.9644   /
      data amo   /  48.0000   /
      data amco2 /  44.0000   /
C
C-----------------------------------------------------------------------
C
C Compute solar distance factor and cosine solar zenith angle usi
C day value where a round day (such as 213.0) refers to 0z at
C Greenwich longitude.
C
C Use formulas from Paltridge, G.W. and C.M.R. Platt 1976: Radiative
C Processes in Meterology and Climatology, Elsevier Scientific
C Publishing Company, New York  p. 57, p. 62,63.
C
C Compute eccentricity factor (sun-earth distance factor)
C
      theta = 2.*pie*calday/dayspy
      eccf  = 1.000110 + .034221*dcos(theta) + .001280*dsin(theta) +
     $         .000719*dcos(2.*theta) + .000077*dsin(2.*theta)
C
C Convert pressure from pascals to dynes/cm2
C
      do k=1,plev
         do i=1,plon
            pmidrd(i,k) = pmid(i,k)*10.0
            pintrd(i,k) = pint(i,k)*10.0
         end do
      end do
      do i=1,plon
         pintrd(i,plevp) = pint(i,plevp)*10.0
      end do
C
C Compute path quantities used in the longwave radiation:
C
      vmmr  = amco2/amd
      cpwpl = vmmr*0.5/(gravit*p0)
      do i=1,plon
         plh2o(i,1)  = rgsslp*h2ommr(i,1)*pintrd(i,1)*pintrd(i,1)
         plco2(i,1)  = co2vmr*cpwpl*pintrd(i,1)*pintrd(i,1)
         tclrsf(i,1) = 1.
      end do
      do k=1,plev
         do i=1,plon
            plh2o(i,k+1)  = plh2o(i,k) + rgsslp*
     $             (pintrd(i,k+1)**2 - pintrd(i,k)**2)*h2ommr(i,k)
            plco2(i,k+1)  = co2vmr*cpwpl*pintrd(i,k+1)**2
            tclrsf(i,k+1) = tclrsf(i,k)*(1.-cld(i,k+1))
         end do
      end do
C
C Convert ozone volume mixing ratio to mass mixing ratio:
C
      vmmr = amo/amd
      do k=1,plev
         do i=1,plon
            o3mmr(i,k) = vmmr*o3vmr(i,k)
         end do
      end do
C
      return
      end

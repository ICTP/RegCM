C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C    This file is part of RegCM model.
C
C    RegCM model is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    RegCM model is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
C
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      subroutine radded(coszrs  ,trayoslp,pflx    ,abh2o   ,abo3    ,
     $                  abco2   ,abo2    ,uh2o    ,uo3     ,uco2    ,
     $                  uo2     ,tauxcl  ,wcl     ,gcl     ,fcl     ,
     $                  tauxci  ,wci     ,gci     ,fci     ,
     $                  tauxar_mixs,tauasc_mixs ,gtota_mixs,ftota_mixs,
     $                  nloop   ,is      ,
     $                  ie      ,rdir    ,rdif    ,tdir    ,tdif    ,
     $                  explay  ,exptdn  ,rdndif  ,tottrn  )
C-----------------------------------------------------------------------
C
C Computes layer reflectivities and transmissivities, from the top down
C to the surface using the delta-Eddington solutions for each layer;
C adds layers from top down to surface as well.
C
C If total transmission to the interface above a particular layer is
C less than trmin, then no further delta-Eddington solutions are
C evaluated for layers below
C
C For more details , see Briegleb, Bruce P., 1992: Delta-Eddington
C Approximation for Solar Radiation in the NCAR Community Climate Model,
C Journal of Geophysical Research, Vol 97, D7, pp7603-7612).
C
C---------------------------Code history--------------------------------
C
C Original version:  B. Briegleb
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Kiehl, B. Briegleb, August 1992
C
C-----------------------------------------------------------------------
c
c $Id: radded.F,v 1.3 1995/03/17 18:54:11 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
c
C------------------------------Parameters-------------------------------
c
#     include "parameter.inc"
c
#     include "include/parrad.cb"
C
C-----------------------------------------------------------------------
C
C Minimum total transmission below which no layer computation are done:
C
      real(kind=8)  trmin,          ! Minimum total transmission allowed
     $      wray,           ! Rayleigh single scatter albedo
     $      gray,           ! Rayleigh asymetry parameter
     $      fray            ! Rayleigh forward scattered fraction
      parameter (trmin = 1.e-3,
     $           wray = 0.999999,
     $           gray = 0.0,
     $           fray = 0.1)
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real(kind=8)  coszrs(plond),         ! Cosine zenith angle
     $     trayoslp,              ! Tray/sslp
     $     pflx(plond,0:plevp),   ! Interface pressure
     $     abh2o,                 ! Absorption coefficiant for h2o
     $     abo3 ,                 ! Absorption coefficiant for o3
     $     abco2,                 ! Absorption coefficiant for co2
     $     abo2 ,                 ! Absorption coefficiant for o2
     $     uh2o(plond,0:plev),    ! Layer absorber amount of h2o
     $     uo3(plond,0:plev),     ! Layer absorber amount of  o3
     $     uco2(plond,0:plev),    ! Layer absorber amount of co2
     $     uo2(plond,0:plev)      ! Layer absorber amount of  o2
      real(kind=8)  tauxcl(plond,0:plev),  ! Cloud extinction optical depth
     $     wcl(plond,0:plev),      ! Cloud single scattering albedo
     $     gcl(plond,0:plev),      ! Cloud assymetry parameter
     $     fcl(plond,0:plev),       ! Cloud forward scattered fraction
     $     tauxci(plond,0:plev),  ! Cloud extinction optical depth
     $     wci(plond,0:plev),      ! Cloud single scattering albedo
     $     gci(plond,0:plev),      ! Cloud assymetry parameter
     $     fci(plond,0:plev)
      real(kind=8)  tauxar_mixs(plond,0:plev), ! Aerosol extinction optical depth
     $        tauasc_mixs(plond,0:plev), ! Aerosol single scattering albedo
     $     gtota_mixs(plond,0:plev),     ! Aerosol assymetry parameter
     $     ftota_mixs(plond,0:plev)      ! Aerosol forward scattered fraction
      integer nloop,              ! Number of loops (1 or 2)
     $        is(2),              ! Starting index for 1 or 2 loops
     $        ie(2)               ! Ending index for 1 or 2 loops
C
C Input/Output arguments
C
C Following variables are defined for each layer; 0 refers to extra
C layer above top of model:
C
      real(kind=8)  rdir(plond,0:plev),    ! Layer reflectivity to direct rad
     $     rdif(plond,0:plev),    ! Layer refflectivity to diffuse rad
     $     tdir(plond,0:plev),    ! Layer transmission to direct rad
     $     tdif(plond,0:plev),    ! Layer transmission to diffuse rad
     $     explay(plond,0:plev)   ! Solar beam exp transm for layer
C
C (Note that the following variables are defined on interfaces, with the
C  index k referring to the top interface of the kth layer:
C  exptdn,rdndif,tottrn; for example, tottrn(k=5) refers to the total
C  transmission to the top interface of the 5th layer; plevp refers to
C  the earth surface)
C
      real(kind=8)  rdndif(plond,0:plevp), ! Added dif ref for layers above
     $     exptdn(plond,0:plevp), ! Solar beam exp down transm from top
     $     tottrn(plond,0:plevp)  ! Total transmission for layers above
C
C---------------------------Local variables-----------------------------
C
      integer i,            ! Longitude index
     $        k,            ! Level index
     $        nn,           ! Index of longitude loops (max=nloop)
     $        ii,           ! Longitude index
     $        nval,         ! Number of long values satisfying criteria
     $        index(plond)  ! Array of longitude indices
C
      real(kind=8)  taugab(plond),   ! Layer total gas absorption optical depth
     $     tauray(plond),   ! Layer rayleigh optical depth
     $     taucsc,          ! Layer cloud scattering optical depth
     $     tautot,          ! Total layer optical depth
     $     wtot,            ! Total layer single scatter albedo
     $     gtot,            ! Total layer asymmetry parameter
     $     ftot             ! Total layer forward scatter fraction
C
      real(kind=8)  wtau,            !  rayleigh layer scattering optical depth
     $     wt,              !  layer total single scattering albedo
     $     ts,              !  layer scaled extinction optical depth
     $     ws,              !  layer scaled single scattering albedo
     $     gs               !  layer scaled asymmetry parameter
C
      real(kind=8)  rdenom,          !  mulitiple scattering term
     $     rdirexp,         !  layer direct ref times exp transmission
     $     tdnmexp          !  total transmission minus exp transmission
C
C---------------------------Statement functions-------------------------
C
C Statement functions and other local variables
C
      real(kind=8)  alpha,           ! Term in direct reflect and transmissivity
     $     gamma,           ! Term in direct reflect and transmissivity
     $     el,              ! Term in alpha,gamma,n,u
     $     taus,            ! Scaled extinction optical depth
     $     omgs,            ! Scaled single particle scattering albedo
     $     asys,            ! Scaled asymmetry parameter
     $     u,               ! Term in diffuse reflect and transmissivity
     $     n,               ! Term in diffuse reflect and transmissivity
     $     lm,              ! Temporary for el
     $     ne               ! Temporary for n
      real(kind=8)  w,               ! Dummy argument for statement function
     $     uu,              ! Dummy argument for statement function
     $     g,               ! Dummy argument for statement function
     $     e,               ! Dummy argument for statement function
     $     f,               ! Dummy argument for statement function
     $     t,               ! Dummy argument for statement function
     $     et               ! Dummy argument for statement function
C
C Intermediate terms for delta-eddington solution
C
      real(kind=8)  alp,             ! Temporary for alpha
     $     gam,             ! Temporary for gamma
     $     ue,              ! Temporary for u
     $     arg,             ! Exponential argument
     $     extins,          ! Extinction
     $     amg,             ! Alp - gam
     $     apg              ! Alp + gam
C
      alpha(w,uu,g,e) = .75*w*uu*((1. + g*(1-w))/(1. - e*e*uu*uu))
      gamma(w,uu,g,e) = .50*w*((3.*g*(1.-w)*uu*uu + 1.)/(1.-e*e*uu*uu))
      el(w,g)         = dsqrt(3.*(1-w)*(1. - w*g))
      taus(w,f,t)     = (1. - w*f)*t
      omgs(w,f)       = (1. - f)*w/(1. - w*f)
      asys(g,f)       = (g - f)/(1. - f)
      u(w,g,e)        = 1.5*(1. - w*g)/e
      n(uu,et)        = ((uu+1.)*(uu+1.)/et ) - ((uu-1.)*(uu-1.)*et)
C
C-----------------------------------------------------------------------
C
C Initialize all total transmission values to 0, so that nighttime
C values from previous computations are not used:
C
      tottrn = 0.d0
C
C Compute total direct beam transmission, total transmission, and
C reflectivity for diffuse radiation (from below) for all layers above
C each interface by starting from the top and adding layers down:
C
C For the extra layer above model top:
C
      do 200 nn=1,nloop
         do 100 i=is(nn),ie(nn)
C
            tauray(i) = trayoslp*(pflx(i,1)-pflx(i,0))
            taugab(i) = abh2o*uh2o(i,0) + abo3*uo3(i,0) +
     $                  abco2*uco2(i,0) + abo2*uo2(i,0)
C
            tautot  = tauxcl(i,0) + tauxci(i,0) + tauray(i) + taugab(i)
     $                            + tauxar_mixs(i,0)
            taucsc  = tauxcl(i,0)*wcl(i,0)+tauxci(i,0)*wci(i,0)
     $                            + tauasc_mixs (i,0)
            wtau    = wray*tauray(i)
            wt      = wtau + taucsc
            wtot = wt/tautot
            gtot = (wtau*gray + gcl(i,0)*tauxcl(i,0)*wcl(i,0) +
     $                          gci(i,0)*tauxci(i,0)*wci(i,0) +
     $                          gtota_mixs(i,0)
     $                                                         )/wt
            ftot = (wtau*fray + fcl(i,0)*tauxcl(i,0)*wcl(i,0) +
     $                          fci(i,0)*tauxci(i,0)*wci(i,0) +
     $                          ftota_mixs(i,0)
     $                                                         )/wt
C
            ts   = taus(wtot,ftot,tautot)
            ws   = omgs(wtot,ftot)
            gs   = asys(gtot,ftot)
            lm   = el(ws,gs)
            alp  = alpha(ws,coszrs(i),gs,lm)
            gam  = gamma(ws,coszrs(i),gs,lm)
            ue   = u(ws,gs,lm)
C
C Limit argument of exponential to 25, in case lm*ts very large:
C
            arg  = dmin1(lm*ts,25.d0)
            extins = dexp(-arg)
            ne = n(ue,extins)
C
            rdif(i,0) = (ue+1.)*(ue-1.)*(1./extins - extins)/ne
            tdif(i,0) = 4.*ue/ne
C
C Limit argument of exponential to 25, in case coszrs is very small:
C
            arg       = dmin1(ts/coszrs(i),25.d0)
            explay(i,0) = dexp(-arg)
C
            apg = alp + gam
            amg = alp - gam
            rdir(i,0) = amg*(tdif(i,0)*explay(i,0) - 1.) + apg*rdif(i,0)
            tdir(i,0) = apg*tdif(i,0) +
     $                  (amg*rdif(i,0) - (apg-1.))*explay(i,0)
C
C Under rare conditions, reflectivies and transmissivities can be
C negative; zero out any negative values
C
            rdir(i,0) = dmax1(rdir(i,0),0.d0)
            tdir(i,0) = dmax1(tdir(i,0),0.d0)
            rdif(i,0) = dmax1(rdif(i,0),0.d0)
            tdif(i,0) = dmax1(tdif(i,0),0.d0)
C
C Initialize top interface of extra layer:
C
            exptdn(i,0) =   1.0
            rdndif(i,0) =   0.0
            tottrn(i,0) =   1.0
C
            rdndif(i,1) = rdif(i,0)
            tottrn(i,1) = tdir(i,0)
C
  100    continue
  200 continue
C
C Now, continue down one layer at a time; if the total transmission to
C the interface just above a given layer is less than trmin, then no
C delta-eddington computation for that layer is done:
C
      do 400 k=1,plev
C
C Initialize current layer properties to zero; only if total
C transmission to the top interface of the current layer exceeds the
C minimum, will these values be computed below:
C
         do nn=1,nloop
            do i=is(nn),ie(nn)
C
               rdir(i,k)   =  0.0
               rdif(i,k)   =  0.0
               tdir(i,k)   =  0.0
               tdif(i,k)   =  0.0
               explay(i,k) =  0.0
C
C Calculates the solar beam transmission, total transmission, and
C reflectivity for diffuse radiation from below at the top of the
C current layer:
C
               exptdn(i,k) = exptdn(i,k-1)*explay(i,k-1)
cKN   modified below (for computational stability)
c               rdenom      = 1./(1. - rdif(i,k-1)*rdndif(i,k-1))
      rdenom  = 1./(1.-dmin1(rdif(i,k-1)*rdndif(i,k-1),0.999999d0))
cKN   modified above
               rdirexp     = rdir(i,k-1)*exptdn(i,k-1)
               tdnmexp     = tottrn(i,k-1) - exptdn(i,k-1)
               tottrn(i,k) = exptdn(i,k-1)*tdir(i,k-1) + tdif(i,k-1)*
     $                      (tdnmexp + rdndif(i,k-1)*rdirexp)*rdenom
               rdndif(i,k) = rdif(i,k-1)  +
     $                (rdndif(i,k-1)*tdif(i,k-1))*(tdif(i,k-1)*rdenom)
C
            end do
         end do
C
C Compute next layer delta-eddington solution only if total transmission
C of radiation to the interface just above the layer exceeds trmin.
C
         call whenfgt(plon,tottrn(1,k),1,trmin,index,nval)
         if(nval.gt.0) then
CCDIR$ IVDEP
            do 300 ii=1,nval
               i=index(ii)
C
               tauray(i) = trayoslp*(pflx(i,k+1)-pflx(i,k))
               taugab(i) = abh2o*uh2o(i,k) + abo3*uo3(i,k) +
     $                     abco2*uco2(i,k) + abo2*uo2(i,k)
C
               tautot = tauxcl(i,k) + tauxci(i,k) +
     $                  tauray(i) + taugab(i) + tauxar_mixs(i,k)
               taucsc    = tauxcl(i,k)*wcl(i,k) + tauxci(i,k)*wci(i,k)
     $                            + tauasc_mixs(i,k)
               wtau      = wray*tauray(i)
               wt        = wtau + taucsc
               wtot   = wt/tautot
               gtot   = (wtau*gray + gcl(i,k)*wcl(i,k)*tauxcl(i,k)
     $                   + gci(i,k)*wci(i,k)*tauxci(i,k)
     $                   + gtota_mixs(i,k)
     $                                                  )/wt
               ftot   = (wtau*fray + fcl(i,k)*wcl(i,k)*tauxcl(i,k)
     $                   + fci(i,k)*wci(i,k)*tauxci(i,k)
     $                   + ftota_mixs(i,k)
     $                                                  )/wt
C
               ts   = taus(wtot,ftot,tautot)
               ws   = omgs(wtot,ftot)
               gs   = asys(gtot,ftot)
               lm   = el(ws,gs)
               alp  = alpha(ws,coszrs(i),gs,lm)
               gam  = gamma(ws,coszrs(i),gs,lm)
               ue   = u(ws,gs,lm)
C
C Limit argument of exponential to 25, in case lm very large:
C
               arg  = dmin1(lm*ts,25.d0)
               extins = dexp(-arg)
               ne = n(ue,extins)
C
               rdif(i,k) = (ue+1.)*(ue-1.)*(1./extins - extins)/ne
               tdif(i,k)   =   4.*ue/ne
C
C Limit argument of exponential to 25, in case coszrs is very small:
C
               arg       = dmin1(ts/coszrs(i),25.d0)
               explay(i,k) = dexp(-arg)
C
               apg = alp + gam
               amg = alp - gam
               rdir(i,k) = amg*(tdif(i,k)*explay(i,k) - 1.) +
     $                     apg*rdif(i,k)
               tdir(i,k) = apg*tdif(i,k) +
     $                     (amg*rdif(i,k) - (apg-1.))*explay(i,k)
C
C Under rare conditions, reflectivies and transmissivities can be
C negative; zero out any negative values
C
               rdir(i,k) = dmax1(rdir(i,k),0.d0)
               tdir(i,k) = dmax1(tdir(i,k),0.d0)
               rdif(i,k) = dmax1(rdif(i,k),0.d0)
               tdif(i,k) = dmax1(tdif(i,k),0.d0)
  300       continue
         end if
C
  400 continue
C
C Compute total direct beam transmission, total transmission, and
C reflectivity for diffuse radiation (from below) for all layers
C above the surface:
C
      k = plevp
      do nn=1,nloop
         do i=is(nn),ie(nn)
            exptdn(i,k) = exptdn(i,k-1)*explay(i,k-1)
cKN   modified below (for computational stability)
c            rdenom = 1./(1. - rdif(i,k-1)*rdndif(i,k-1))
      rdenom = 1./(1.-dmin1(rdif(i,k-1)*rdndif(i,k-1),0.999999d0))
cKN   modified above
            rdirexp = rdir(i,k-1)*exptdn(i,k-1)
            tdnmexp = tottrn(i,k-1) - exptdn(i,k-1)
            tottrn(i,k) = exptdn(i,k-1)*tdir(i,k-1) + tdif(i,k-1)*
     $                   (tdnmexp + rdndif(i,k-1)*rdirexp)*rdenom
            rdndif(i,k) = rdif(i,k-1)  +
     $               (rdndif(i,k-1)*tdif(i,k-1))*(tdif(i,k-1)*rdenom)
         end do
      end do
C
      return
      end subroutine radded

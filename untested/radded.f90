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
      use mod_regcm_param
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
!     transmission to the top interface of the 5th layer; kx + 1 refers
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
      real(8) , dimension(ix - 1) :: coszrs
      real(8) , dimension(ix - 1,0:kx) :: explay , fci , fcl ,         &
           & ftota_mixs , gci , gcl , gtota_mixs , rdif , rdir ,        &
           & tauasc_mixs , tauxar_mixs , tauxci , tauxcl , tdif , tdir ,&
           & uco2 , uh2o , uo2 , uo3 , wci , wcl
      real(8) , dimension(ix - 1,0:kx + 1) :: exptdn , pflx , rdndif ,    &
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
      integer , dimension(ix - 1) :: indx
      real(8) , dimension(ix - 1) :: taugab , tauray
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
      do k = 1 , kx
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
        call whenfgt(ix - 1,tottrn(1,k),1,trmin,indx,nval)
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
      k = kx + 1
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

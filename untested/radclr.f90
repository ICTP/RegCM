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
      use regcm_param
      use parrad
      use aerm
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
! rdif     - Layer refflectivity to diffuse rad
! tdir     - Layer transmission to direct rad
! tdif     - Layer transmission to diffuse rad
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
      real(8) , dimension(plond) :: coszrs , ftota_mix_css ,            &
                                  & gtota_mix_css , tauasc_mix_css ,    &
                                  & tauxar_mix_css , utco2 , uth2o ,    &
                                  & uto2 , uto3
      real(8) , dimension(plond,0:plev) :: explay , rdif , rdir , tdif ,&
           & tdir
      real(8) , dimension(plond,0:plevp) :: exptdn , pflx , rdndif ,    &
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
      integer , dimension(plond) :: indx
      real(8) , dimension(plond) :: taugab , tauray
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
      do k = 0 , plevp
        do i = 1 , plond
          tottrn(i,k) = 0.
        end do
      end do
!
!     Compute total direct beam transmission, total transmission, and
!     reflectivity for diffuse radiation (from below) for all layers
!     above each interface by starting from the top and adding layers
!     down:
!
!     The top layer is assumed to be a purely absorbing ozone layer, and
!     that the mean diffusivity for diffuse transmission is 1.66:
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
!         Same limit for diffuse transmission:
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
!           and reflectivity for diffuse radiation from below at the
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
        call whenfgt(plon,tottrn(1,k),1,trmin,indx,nval)
        if ( nval.gt.0 ) then
!CDIR$    IVDEP
          do ii = 1 , nval
            i = indx(ii)
!
!           Remember, no ozone absorption in this layer:
!
            tauray(i) = trayoslp*pflx(i,plevp)
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
!     reflectivity for diffuse radiation (from below) for both layers
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

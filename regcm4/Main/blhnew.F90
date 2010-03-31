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
 
      subroutine blhnew
!
! ------------------------------------------------------------
! this routine computes the boundary layer eddy diffusivities
! for momentum, heat and moisture and the counter-gradient
! terms for heat and moisture.
!
! reference : holtslag, de bruijn and pan - mwr - 8/90
!
! input arguments :  j       longitudinal position index
!                    ubx3d   u wind component
!                    vbx3d   v wind component
!                    thx3d   potential temperature
!                    thvx    virtual potential temperature
!                    za      height of half sigma levels
!                    f       coriolis parameter
!                    shum    specific humidity
!                    xhfx    sensible heat flux
!                    xqfx    sfc kinematic moisture flux
!                    th10    virt. pot. temp. at 10m
!                    hfxv    surface virtual heat flux
!                    obklen  monin obukov length
!                    ustr    friction velocity
!                    kzo     minimum eddy diffusivity
!
! input/output
! arguments :        therm   thermal temperature excess
!
! output arguments : cgh     counter-gradient term for heat
!                    cgq     counter-gradient term for moisture
!                    kvm     eddy diffusivity for momentum
!                    kvh     eddy diffusivity for heat
!                    kvq     eddy diffusivity for moisture
!                    zpbl     boundary layer height
!
      use mod_regcm_param
      use mod_param2
      use mod_param3 , only : kt
      use mod_main
      use mod_slice
      use mod_pbldim
      use mod_blh_tmp
      use mod_constants , only : gti , vonkar
      implicit none
!
! Local variables
!
      real(8) :: betah , betam , betas , binh , binm , ccon , fak ,     &
               & fak1 , fak2 , fht , xfmt , kzo , onet , pblk , pblk1 , &
               & pblk2 , pfcor , phpblm , pink , pr , ricr , sffrac ,   &
               & therm2 , tkv , tlv , ttkl , vv , vvl , wsc , z , zh ,  &
               & zl , zm , zp , zzh , zzhnew , zzhnew2
      integer :: i , j , k , k2
      real(8) , dimension(iy,kz) :: ri
      real(8) , dimension(iy) :: therm
!
      data kzo/1./
!
!
!     ------------------------------------------------------------
!
!     real(kind=8)  cgq(iy,kz)
!     -----------------------------------------------------------
 
!     gravity
!     coef. of proportionality and lower % of bl in sfc layer
      data fak , sffrac/8.5 , 0.1/
!     beta coefs. for momentum, stable conditions and heat
      data betam , betas , betah/15.0 , 5.0 , 15.0/
!     power in formula for k and critical ri for judging stability
      data pink , ricr/2.0 , 0.25/
 
!     exponent : one third
      onet = 1./3.
!     set constants
      ccon = fak*sffrac*vonkar
      binm = betam*sffrac
      binh = betah*sffrac
#ifdef MPP1
      do j = jbegin , jendx
#else
      do j = 2 , jxm1
#endif
 
!       ****note: kt, max no. of pbl levels, calculated in param
!       ******   compute richardson number
        do i = 2 , iym1
          therm(i) = 0.0
        end do
 
        do k = kz , kt , -1
          do i = 2 , iym1
            vv = ubx3d(i,k,j)*ubx3d(i,k,j) + vbx3d(i,k,j)*vbx3d(i,k,j)
            ri(i,k) = gti*(thvx(i,k,j)-th10(i,j))*za(i,k,j)/            &
                    & (th10(i,j)*vv)
          end do
        end do
 
!       ******   first, set bl height to height of lowest model level
        do i = 2 , iym1
          zpbl(i,j) = za(i,kz,j)
        end do
 
!       ******   looking for bl top
        do k = kz , kt + 1 , -1
          k2 = k - 1
          do i = 2 , iym1
!     ******   bl height lies between this level and the last
!     ******   use linear interp. of rich. no. to height of ri=ricr
            if ( (ri(i,k).lt.ricr) .and. (ri(i,k2).ge.ricr) ) zpbl(i,j) &
               & = za(i,k,j) + (za(i,k2,j)-za(i,k,j))                   &
                 & *((ricr-ri(i,k))/(ri(i,k2)-ri(i,k)))
          end do
        end do
 
        do i = 2 , iym1
!     ******   set bl top to highest allowable model layer
          if ( ri(i,kt).lt.ricr ) zpbl(i,j) = za(i,kt,j)
        end do
 
!       ******   recompute richardson no. at lowest model level
        do i = 2 , iym1
          if ( hfxv(i,j).gt.0. ) then
!           ******   estimate of convective velocity scale
            xfmt = (1.0-(binm*zpbl(i,j)/obklen(i,j)))**onet
            wsc = ustr(i,j)*xfmt
!           ******   thermal temperature excess
            therm(i) = (xhfx(i,j)+0.61*thx3d(i,kz,j)*xqfx(i,j))*fak/wsc
            vvl = ubx3d(i,kz,j)*ubx3d(i,kz,j) + vbx3d(i,kz,j)           &
                & *vbx3d(i,kz,j)
            ri(i,kz) = -gti*therm(i)*za(i,kz,j)/(th10(i,j)*vvl)
          end if
        end do
 
!       ******   recompute richardson no. at other model levels
        do k = kz - 1 , kt , -1
          do i = 2 , iym1
            if ( hfxv(i,j).gt.0. ) then
              tlv = th10(i,j) + therm(i)
              tkv = thx3d(i,k,j)                                        &
                  & *(1.0+0.61*(qvb3d(i,k,j)/(qvb3d(i,k,j)+1)))
              ttkl = tkv - tlv
              vv = ubx3d(i,k,j)*ubx3d(i,k,j) + vbx3d(i,k,j)*vbx3d(i,k,j)
              ri(i,k) = gti*ttkl*za(i,k,j)/(th10(i,j)*vv)
            end if
          end do
        end do
 
!       ******   improve estimate of bl height under convective
!       conditions ******   using convective temperature excess (therm)
        do k = kz , kt + 1 , -1
          k2 = k - 1
          do i = 2 , iym1
            if ( hfxv(i,j).gt.0. ) then
!     ******   bl height lies between this level and the last
!     ******   use linear interp. of rich. no. to height of ri=ricr
              if ( (ri(i,k).lt.ricr) .and. (ri(i,k2).ge.ricr) )         &
                 & zpbl(i,j) = za(i,k,j) + (za(i,k2,j)-za(i,k,j))       &
                             & *((ricr-ri(i,k))/(ri(i,k2)-ri(i,k)))
            end if
          end do
        end do
 
        do i = 2 , iym1
          if ( hfxv(i,j).gt.0. ) then
!     ******   set bl top to highest allowable model layer
            if ( ri(i,kt).lt.ricr ) zpbl(i,j) = za(i,kt,j)
          end if
        end do
 
!       ******   limit bl height to be at least mech. mixing depth
        do i = 2 , iym1
!         ******   limit coriolis parameter to value at 10 deg. latitude
          pfcor = dmax1(dabs(f(i,j)),2.546D-5)
!         ******   compute mechanical mixing depth,
!         ******   set to lowest model level if lower
          phpblm = 0.07*ustr(i,j)/pfcor
          phpblm = dmax1(phpblm,za(i,kz,j))
          zpbl(i,j) = dmax1(zpbl(i,j),phpblm)
        end do
 
        do k = kz , kt + 1 , -1
          k2 = k - 1
          do i = 2 , iym1
            pblk = 0.0
            zm = za(i,k,j)
            zp = za(i,k2,j)
            if ( zm.lt.zpbl(i,j) ) then
              zp = dmin1(zp,zpbl(i,j))
              z = 0.5*(zm+zp)
              zh = z/zpbl(i,j)
              zl = z/obklen(i,j)
              if ( zh.le.1. ) then
                zzh = 1. - zh
                zzh = zzh**pink
!xexp4          zzhnew = zpbl(i,j)*(1.-zh)*zh**1.5
!xexp5          zzhnew = 0.5*zpbl(i,j)*(1.-zh)*zh**1.5
!xexp6          zzhnew = 1. - zh
!xexp7          zzhnew =0.5* (1. - zh)
!Sara
!               zzhnew =0.25* (1. - zh)
!               zzhnew =0.75* (1. - zh)
!Sara_
                zzhnew = 0.25*(1.-zh)
!xexp10         zzhnew =zh * (1. - zh)**2
!chem
                if ( ichem.eq.1 ) zzhnew2 = (1.-zh)**2
!chem_
              else
                zzh = 0.
                zzhnew = 0.
!chem
                zzhnew2 = 0.
!chem_
              end if
              fak1 = ustr(i,j)*zpbl(i,j)*vonkar
              if ( hfxv(i,j).le.0. ) then
!**             stable and neutral conditions
!**             igroup = 1
 
!**             prevent pblk from becoming too small in very stable
!               conditions
                if ( zl.le.1. ) then
                  pblk = fak1*zh*zzh/(1.+betas*zl)
!xexp5            pblk1 = vonkar * ustr(i,j) / (1.+betas*zl) * zzhnew
                  pblk1 = fak1*zh*zzhnew/(1.+betas*zl)
!chem
                  if ( ichem.eq.1 )                                     &
                     & pblk2 = fak1*zh*zzhnew2/(1.+betas*zl)
!chem_
                else
                  pblk = fak1*zh*zzh/(betas+zl)
!xexp5            pblk1 = vonkar * ustr(i,j) / (betas+zl) * zzhnew
                  pblk1 = fak1*zh*zzhnew/(betas+zl)
!chem
                  if ( ichem.eq.1 ) pblk2 = fak1*zh*zzhnew2/(betas+zl)
!chem_
                end if
!**             compute eddy diffusivities
                kvm(i,k,j) = dmax1(pblk,kzo)
                kvh(i,k,j) = kvm(i,k,j)
                kvq(i,k,j) = dmax1(pblk1,kzo)
!chem
                if ( ichem.eq.1 ) kvc(i,k,j) = dmax1(pblk2,kzo)
!chem_
!**             compute counter-gradient term
                cgh(i,k,j) = 0.0
!               cgq(i,k) = 0.0
              else
!**             unstable conditions
 
!**             compute counter gradient term
                if ( zh.ge.sffrac ) then
!**               igroup = 2
                  xfmt = (1.-binm*zpbl(i,j)/obklen(i,j))**onet
                  fht = dsqrt(1.-binh*zpbl(i,j)/obklen(i,j))
                  wsc = ustr(i,j)*xfmt
                  pr = (xfmt/fht) + ccon
                  fak2 = wsc*zpbl(i,j)*vonkar
                  pblk = fak2*zh*zzh
!xexp5            pblk1 = vonkar * wsc * zzhnew
                  pblk1 = fak2*zh*zzhnew
!chem
                  if ( ichem.eq.1 ) pblk2 = fak2*zh*zzhnew2
!chem_
                  therm2 = fak/(zpbl(i,j)*wsc)
                  cgh(i,k,j) = hfxv(i,j)*therm2
!                 cgq(i,k) = xqfx(i,j)*therm2
!                 cgq(i,k) = 0.0
                else
!**               igroup = 3
                  pblk = fak1*zh*zzh*(1.-betam*zl)**onet
!xexp5            pblk1 = vonkar * ustr(i,j) * zzhnew *
!                 (1.-betam*zl)**onet
                  pblk1 = fak1*zh*zzhnew*(1.-betam*zl)**onet
!chem
                  if ( ichem.eq.1 )                                     &
                     & pblk2 = fak1*zh*zzhnew2*(1.-betam*zl)**onet
!chem_
                  pr = ((1.-betam*zl)**onet)/dsqrt(1.-betah*zl)
                  cgh(i,k,j) = 0.0
!                 cgq(i,k) = 0.0
                end if
 
!**             compute eddy diffusivities
                kvm(i,k,j) = dmax1(pblk,kzo)
                kvh(i,k,j) = dmax1((pblk/pr),kzo)
!               kvq(i,k,j) = kvh(i,k,j)
                kvq(i,k,j) = dmax1(pblk1,kzo)
!chem
                if ( ichem.eq.1 ) kvc(i,k,j) = dmax1(pblk2,kzo)
!chem_
 
              end if
            end if
          end do
        end do
      end do
 
      end subroutine blhnew

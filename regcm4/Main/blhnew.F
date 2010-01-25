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

c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine blhnew
c
c ------------------------------------------------------------
c this routine computes the boundary layer eddy diffusivities
c for momentum, heat and moisture and the counter-gradient
c terms for heat and moisture.
c
c reference : holtslag, de bruijn and pan - mwr - 8/90
c
c input arguments :  j       longitudinal position index
c                    ubx3d   u wind component
c                    vbx3d   v wind component
c                    thx3d   potential temperature
c                    thvx    virtual potential temperature
c                    za      height of half sigma levels
c                    f       coriolis parameter
c                    shum    specific humidity
c                    xhfx    sensible heat flux
c                    xqfx    sfc kinematic moisture flux
c                    th10    virt. pot. temp. at 10m
c                    hfxv    surface virtual heat flux
c                    obklen  monin obukov length
c                    ustr    friction velocity
c                    kzo     minimum eddy diffusivity
c
c input/output
c arguments :        therm   thermal temperature excess
c
c output arguments : cgh     counter-gradient term for heat
c                    cgq     counter-gradient term for moisture
c                    kvm     eddy diffusivity for momentum
c                    kvh     eddy diffusivity for heat
c                    kvq     eddy diffusivity for moisture
c                    zpbl     boundary layer height
c
      implicit none
      integer j

#     include "parameter.inc"

#     include "include/param1.cb"
#     include "include/param2.cb"
#     include "include/param3.cb"
#     include "include/main.cb"
#     include "include/slice.cb"
#     include "include/pbldim.cb"

#ifdef MPP1
#     include "mpiregcm.inc"
#endif

      real(kind=8)  kzo
      data    kzo/1./
c
#ifdef MPP1
      real(kind=8)  kvm(ix,kx,jxp),kvh(ix,kx,jxp)
      real(kind=8)  kvq(ix,kx,jxp),kvc(ix,kx,jxp),cgh(ix,kx,jxp)
      real(kind=8)  xhfx(ix,jxp),xqfx(ix,jxp),hfxv(ix,jxp)
      real(kind=8)  obklen(ix,jxp),th10(ix,jxp),ustr(ix,jxp)
#else
      real(kind=8)  kvm(ix,kx,jlx),kvh(ix,kx,jlx)
      real(kind=8)  kvq(ix,kx,jlx),kvc(ix,kx,jlx),cgh(ix,kx,jlx)
      real(kind=8)  xhfx(ix,jx),xqfx(ix,jx),hfxv(ix,jx)
      real(kind=8)  obklen(ix,jx),th10(ix,jx),ustr(ix,jx)
#endif
      COMMON /BLH_TMP/ kvm,kvh,kvq,kvc,cgh
     &       , xhfx,xqfx,obklen,hfxv,th10,ustr
C
c ------------------------------------------------------------
c
c --- input variables
c     real(kind=8)  shum(ix,kx)
c#ifdef MPP1
c      real(kind=8)  xhfx(ix,jxp),xqfx(ix,jxp)
c#else
c      real(kind=8)  xhfx(ix,jx),xqfx(ix,jx)
c#endif
c      real(kind=8)  th10(ix),     hfxv(ix),  obklen(ix),   ustr(ix)
c
c --- input/ouput variables
c      real(kind=8)  therm(ix)
       real(kind=8)  ri(ix,kx),cgq(ix,kx),therm(ix)
cc
c --- output variables
c      real(kind=8)  cgh(ix,kx), cgq(ix,kx)
c      real(kind=8)  kvm(ix,kx), kvh(ix,kx), kvq(ix,kx), kvc(ix,kx)
c
c -----------------------------------------------------------

c     gravity
c     coef. of proportionality and lower % of bl in sfc layer
      real(kind=8)  fak,sffrac
      data    fak,sffrac         /8.5,0.1/
c     beta coefs. for momentum, stable conditions and heat
      real(kind=8)  betam,betas,betah
      data    betam,betas,betah  /15.0,5.0,15.0/
c     power in formula for k and critical ri for judging stability
      real(kind=8)  pink,ricr
      data    pink,ricr          /2.0,0.25/
      real(kind=8)  onet,ccon,binm,binh,vv,fmt,wsc,vvl,tlv,tkv,ttkl
      integer i,k,k2
      real(kind=8)  pfcor,phpblm,pblk,zm,zp,z,zh,zl,zzh,zzhnew,zzhnew2
      real(kind=8)  pblk1,fht,pr,fak2,therm2,pblk2,fak1

c     exponent : one third
      onet = 1./3.
c     set constants
      ccon = fak * sffrac * karman
      binm = betam * sffrac
      binh = betah * sffrac
#ifdef MPP1
      do 1079 j=jbegin,jendx
#else
      do 1079 j=2,jlx
#endif

c ****note: kt, max no. of pbl levels, calculated in param
c ******   compute richardson number
      do i=2,ilx
        therm(i) = 0.0
      enddo

      do 50 k=kx,kt,-1
      do 50 i=2,ilx
       vv = ubx3d(i,k,j)*ubx3d(i,k,j) + vbx3d(i,k,j)*vbx3d(i,k,j)
       ri(i,k) = g * (thvx(i,k,j)-th10(i,j)) * za(i,k,j) /
     1                           (th10(i,j)  * vv)
50    continue

c ******   first, set bl height to height of lowest model level
      do 100 i=2,ilx
       zpbl(i,j) = za(i,kx,j)
100   continue

c ******   looking for bl top
      do 200 k=kx,kt+1,-1
       k2 = k-1
      do 200 i=2,ilx
       if( (ri(i,k).lt.ricr) .and. (ri(i,k2).ge.ricr) ) then
c ******   bl height lies between this level and the last
c ******   use linear interp. of rich. no. to height of ri=ricr
        zpbl(i,j) = za(i,k,j) + (za(i,k2,j) - za(i,k,j)) *
     1          ( (ricr-ri(i,k))/(ri(i,k2)-ri(i,k)) )
       end if
200   continue

      do 300 i=2,ilx
       if (ri(i,kt).lt.ricr) then
c ******   set bl top to highest allowable model layer
        zpbl(i,j) = za(i,kt,j)
       end if
300   continue

c ******   recompute richardson no. at lowest model level
      do 400 i=2,ilx
       if (hfxv(i,j).gt.0.) then
c ******   estimate of convective velocity scale
        fmt = (1.0 - (binm*zpbl(i,j)/obklen(i,j))) ** onet
        wsc = ustr(i,j) * fmt
c ******   thermal temperature excess
        therm(i) = (xhfx(i,j)+0.61*thx3d(i,kx,j)*xqfx(i,j))*fak/wsc
        vvl = ubx3d(i,kx,j)*ubx3d(i,kx,j) + vbx3d(i,kx,j)*vbx3d(i,kx,j)
        ri(i,kx) = -g * therm(i) * za(i,kx,j) / (th10(i,j)*vvl)
       end if
400   continue

c ******   recompute richardson no. at other model levels
      do 420 k=kx-1,kt,-1
      do 420 i=2,ilx
       if (hfxv(i,j).gt.0.) then
        tlv = th10(i,j) + therm(i)
        tkv = thx3d(i,k,j)*(1.0+0.61*(qvb3d(i,k,j)/(qvb3d(i,k,j)+1)))
        ttkl = tkv - tlv
        vv = ubx3d(i,k,j)*ubx3d(i,k,j) + vbx3d(i,k,j)*vbx3d(i,k,j)
        ri(i,k) = g * ttkl * za(i,k,j) / (th10(i,j) * vv)
       end if
420   continue

c ******   improve estimate of bl height under convective conditions
c ******   using convective temperature excess (therm)
      do 500 k=kx,kt+1,-1
       k2 = k-1
      do 500 i=2,ilx
       if (hfxv(i,j).gt.0.) then
        if( (ri(i,k).lt.ricr) .and. (ri(i,k2).ge.ricr) ) then
c ******   bl height lies between this level and the last
c ******   use linear interp. of rich. no. to height of ri=ricr
         zpbl(i,j) = za(i,k,j) + (za(i,k2,j) - za(i,k,j)) *
     1           ( (ricr-ri(i,k))/(ri(i,k2)-ri(i,k)) )
        end if
       end if
500   continue

      do 600 i=2,ilx
       if (hfxv(i,j).gt.0.) then
        if (ri(i,kt).lt.ricr) then
c ******   set bl top to highest allowable model layer
         zpbl(i,j) = za(i,kt,j)
        end if
       end if
600   continue

c ******   limit bl height to be at least mech. mixing depth
      do 700 i=2,ilx
c ******   limit coriolis parameter to value at 10 deg. latitude
       pfcor = dmax1(dabs(f(i,j)),2.546d-5)
c ******   compute mechanical mixing depth,
c ******   set to lowest model level if lower
       phpblm = 0.07 * ustr(i,j) / pfcor
       phpblm = dmax1( phpblm, za(i,kx,j) )
       zpbl(i,j) = dmax1( zpbl(i,j), phpblm )
700   continue

      do 910 k=kx,kt+1,-1
       k2 = k - 1
      do 900 i=2,ilx
       pblk = 0.0
       zm = za(i,k,j)
       zp = za(i,k2,j)
       if (zm.lt.zpbl(i,j)) then
        zp = dmin1( zp, zpbl(i,j) )
        z = 0.5*(zm+zp)
        zh = z/zpbl(i,j)
        zl = z/obklen(i,j)
        if (zh.le.1.) then
         zzh = 1. - zh
         zzh = zzh ** pink
cxexp4         zzhnew = zpbl(i,j)*(1.-zh)*zh**1.5
cxexp5         zzhnew = 0.5*zpbl(i,j)*(1.-zh)*zh**1.5
cxexp6         zzhnew = 1. - zh
cxexp7         zzhnew =0.5* (1. - zh)
!Sara
!        zzhnew =0.25* (1. - zh)
!        zzhnew =0.75* (1. - zh)
!Sara_
         zzhnew =0.25* (1. - zh)
cxexp10         zzhnew =zh * (1. - zh)**2
Cchem
         if(ichem .eq.1) zzhnew2= (1. - zh)**2
Cchem_
        else
         zzh = 0.
         zzhnew = 0.
Cchem
         zzhnew2= 0.
Cchem_
        end if
        fak1 = ustr(i,j) * zpbl(i,j) * karman
        if (hfxv(i,j).le.0.) then
c** stable and neutral conditions
c** igroup = 1

c** prevent pblk from becoming too small in very stable conditions
         if (zl.le.1.) then
          pblk = fak1*zh*zzh/(1.+betas*zl)
cxexp5          pblk1 = karman * ustr(i,j) / (1.+betas*zl) * zzhnew
          pblk1 = fak1*zh*zzhnew/(1.+betas*zl)
Cchem
          if(ichem.eq.1) pblk2 = fak1*zh*zzhnew2/(1.+betas*zl)
Cchem_
         else
          pblk = fak1*zh*zzh/(betas+zl)
cxexp5          pblk1 = karman * ustr(i,j) / (betas+zl) * zzhnew
          pblk1 = fak1*zh*zzhnew/(betas+zl)
Cchem
          if(ichem.eq.1) pblk2 = fak1*zh*zzhnew2/(betas+zl)
Cchem_
         end if
c** compute eddy diffusivities
         kvm(i,k,j) = dmax1(pblk,kzo)
         kvh(i,k,j) = kvm(i,k,j)
         kvq(i,k,j) = dmax1(pblk1,kzo)
Cchem
         if(ichem .eq.1) kvc(i,k,j) = dmax1(pblk2,kzo)
Cchem_
c** compute counter-gradient term
         cgh(i,k,j) = 0.0
         cgq(i,k) = 0.0
        else
c** unstable conditions

c** compute counter gradient term
         if (zh.ge.sffrac) then
c** igroup = 2
          fmt = (1.-binm*zpbl(i,j)/obklen(i,j))**onet
          fht = dsqrt(1.-binh*zpbl(i,j)/obklen(i,j))
          wsc = ustr(i,j)*fmt
          pr = (fmt/fht) + ccon
          fak2 = wsc*zpbl(i,j)*karman
          pblk = fak2*zh*zzh
cxexp5          pblk1 = karman * wsc * zzhnew
          pblk1 = fak2*zh*zzhnew
Cchem
          if(ichem.eq.1) pblk2 = fak2*zh*zzhnew2
Cchem_
          therm2 = fak/(zpbl(i,j)*wsc)
          cgh(i,k,j) = hfxv(i,j)*therm2
c         cgq(i,k) = xqfx(i,j)*therm2
          cgq(i,k) = 0.0
         else
c** igroup = 3
          pblk = fak1*zh*zzh*(1.-betam*zl)**onet
cxexp5          pblk1 = karman * ustr(i,j) * zzhnew * (1.-betam*zl)**onet
          pblk1 = fak1*zh*zzhnew*(1.-betam*zl)**onet
Cchem
          if(ichem.eq.1) pblk2 = fak1*zh*zzhnew2*(1.-betam*zl)**onet
Cchem_
          pr = ((1.-betam*zl)**onet)/dsqrt(1.-betah*zl)
          cgh(i,k,j) = 0.0
          cgq(i,k) = 0.0
         end if

c** compute eddy diffusivities
         kvm(i,k,j) = dmax1(pblk,kzo)
         kvh(i,k,j) = dmax1((pblk/pr),kzo)
c        kvq(i,k,j) = kvh(i,k,j)
         kvq(i,k,j) = dmax1(pblk1,kzo)
Cchem
         if(ichem .eq.1) kvc(i,k,j) = dmax1(pblk2,kzo)
Cchem_

        end if
       end if
900   continue
910   continue
 1079 continue

      return
      end

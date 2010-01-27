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

      subroutine initb
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c     ***  provides constant initial fields to boundary subroutine
c     ***    soil textures are set in soilbc
c     ***    soil   colors are set in albedo
c
c     ***   units are si
c
      implicit none
#     include "parameter.inc"
#     include "include/param1.cb"
#     include "include/param2.cb"
#     include "include/iunits.cb"
#     include "include/bats.cb"
#     include "include/main.cb"
#     include "include/bats.cb2"
#ifdef MPP1
#     include "mpiregcm.inc"
#endif
c
c----------------------------------------------------------------------
c     slmo   : surface moisture availability in fraction of one.
      real(kind=8)  slmo(20)                !BATS land types
c     data slmo/0.75,0.60,0.75,0.75,0.75,
c    &          0.80,0.60,0.10,0.90,0.90,
c    &          0.10,0.95,0.90,1.00,1.00,
c    &          0.60,0.75,0.85,0.85,0.85 /
      data slmo/0.65,0.45,0.60,0.60,0.65,
     &          0.65,0.55,0.10,0.90,0.80,
     &          0.20,0.90,0.90,1.00,1.00,
     &          0.50,0.50,0.65,0.60,0.60/

c ****** vgtran transforms mm4 surface types into bats surface types
cMM4  real(kind=8)  vgtran(13)
cMM4  data vgtran /5.,1.,2.,5.,3.,18.,0.,13.,11.,9.,12.,6.,7./
cMM4  integer ist
      integer is,ill,jll,nlveg,itex,k

c ****** typically resp=1.0 kg/m**2/s and changes by<10% in 10 days
c
#ifdef MPP1
      do 18 jll =1,jendx
#else
      do 18 jll =1,jx-1
#endif
      do 18 ill =1,ix-1
       pptnc(ill,jll) = 0.
       pptc(ill,jll) = 0.
cMM4   ist=nint(satbrt(ill,jll))
cMM4   if(ist.le.13)then
cMM4    ist=ist
cMM4   else
cMM4    ist=7
cMM4   end if
cMM4   veg2d(ill,jll)=vgtran(ist)
ceros  note:  when using bats dataset comment line above
       veg2d(ill,jll)=satbrt(ill,jll)
       if (satbrt(ill,jll).gt.13.9.and.satbrt(ill,jll).lt.15.1)
     &        veg2d(ill,jll) = 0.
18    continue
#ifdef MPP1
      do jll =1,jendx
#else
      do jll =1,jx-1
#endif
      do ill=1,ix-1
      do k=1,NNSG
         veg2d1(k,ill,jll)=satbrt1(k,ill,jll)
         if(satbrt1(k,ill,jll).gt.13.9.and.satbrt1(k,ill,jll).lt.15.1)
     &       veg2d1(k,ill,jll) = 0.
      end do
      end do
      end do

c ******  initialize hostetler lake model
      if (lakemod.eq.1) then
#ifdef MPP1
       call initlk(veg2d1, (ix-1)*NSG, jxp*NSG)
#else
       call initlk(veg2d1, (ix-1)*NSG, (jx-1)*NSG)
#endif
      end if

#ifdef MPP1
      do 20 jll =1,jendx
#else
      do 20 jll =1,jx-1
#endif
      do 20 ill =1,ix-1
      do 20 k=1,NNSG
       if(veg2d1(k,ill,jll).lt.0.5)then
        ocld2d(k,ill,jll)=0.
       else
        ocld2d(k,ill,jll)=1.
       end if
       nlveg=nint(veg2d1(k,ill,jll))
       if(nlveg.eq.0)then
        nlveg=15
       else
        nlveg=nlveg
       end if
csol   itex=int(text2d(k,ill,jll))
       itex=iexsol(nlveg)
       tg2d(k,ill,jll)=tgb(ill,jll)
       tgb2d(k,ill,jll)=tgb(ill,jll)
       taf2d(k,ill,jll)=tgb(ill,jll)
       tlef2d(k,ill,jll)=tgb(ill,jll)

c ******  initialize soil moisture in the 3 layers
       is=nint(satbrt1(k,ill,jll))
       swt2d(k,ill,jll)=deptv(nlveg)*xmopor(itex)*slmo(is)
       srw2d(k,ill,jll)=deprv(nlveg)*xmopor(itex)*slmo(is)
       ssw2d(k,ill,jll)=depuv(nlveg)*xmopor(itex)*slmo(is)

       dew2d(k,ill,jll)=0.
       sag2d(k,ill,jll)=0.
       scv2d(k,ill,jll)=dmax1(snowc(k,ill,jll),0.d0)
       sice2d(k,ill,jll)=0.
       gwet2d(k,ill,jll)=0.5
       sena2d(k,ill,jll)=0.
       evpa2d(k,ill,jll)=0.
       rnos2d(k,ill,jll)=0.
       rno2d(k,ill,jll)=0.
       if(sice2d(k,ill,jll).gt.0.)then
        ocld2d(k,ill,jll)=2.
       end if
       ircp2d(k,ill,jll)=0.
20    continue

#ifdef MPP1
      do 30 jll =1,jendx
#else
      do 30 jll =1,jx-1
#endif
      do 30 ill =1,ix-1
       fsw2d(ill,jll)=0.
       flw2d(ill,jll)=0.
       sabv2d(ill,jll)=0.
       sol2d(ill,jll)=0.
       fswa2d(ill,jll)=0.
       flwa2d(ill,jll)=0.
       prca2d(ill,jll)=0.
       prnca2d(ill,jll)=0.
       svga2d(ill,jll)=0.
       sina2d(ill,jll)=0.
30    continue

      return
      end subroutine initb

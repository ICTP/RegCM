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

      SUBROUTINE LAKEDRV(jslc)
      implicit none
#     include "parameter.inc"
#     include "include/param1.cb"
#     include "include/iunits.cb"
#     include "include/bats.cb"
#     include "include/bats.cb2"
#     include "include/date.cb"

      real(kind=8)  dayl,tl,vl,zl,ql,fsw,flw,prec,hsen
      real(kind=8)  hlat,tgl,evl,aveice,snow
      integer ilake,jslc,n

      do 7 ilake=1,ix-1
      do n=1,NNSG
        if( veg2d1(n,ilake,jslc).eq.14. ) then
          dayl = (NNNNNN-NSTRT0)/4.+(xtime + dtmin) / 1440.
          tl  =       ts1d(n,ilake)
          vl  = dsqrt( us1d(ilake)**2 + vs1d(ilake)**2 )
          zl  =       z1d(n,ilake)
          ql  =       qs1d(n,ilake)
          fsw =       fsw1d(ilake)
          flw = -1. * flw1d(ilake)
          prec=       prca2d(ilake,jslc)   !  units of prec = mm
          hsen= -1. * sent1d(n,ilake)
          hlat= -1. * evpr1d(n,ilake)
          call lake( iutlak, dayl, dtlake, tl, vl, zl, ql, fsw, flw
     a             , hsen, hlat, tgl, evl, prec, aveice, snow)

          tg1d(n,ilake) = tgl
          tgb1d(n,ilake) = tgl
          if(aveice.le.10.)then
            ldoc1d(n,ilake) = 0.
            sice1d(n,ilake) = 0.
            scv1d(n,ilake) = 0.
            sag1d(n,ilake) = 0.
          else
            ldoc1d(n,ilake) = 2.
            sice1d(n,ilake) = aveice    !  units of ice = mm
            scv1d(n,ilake)  = snow    !  units of snow = mm h2o
          end if

        end if
      enddo
7     continue

      return
      end


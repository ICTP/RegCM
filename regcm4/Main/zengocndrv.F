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

      subroutine zengocndrv(j)

      implicit none
#     include "parameter.inc"
#     include "include/param1.cb"
#     include "include/param2.cb"
#     include "include/param3.cb"
#     include "include/main.cb"
#     include "include/pbldim.cb"
#     include "include/pmoist.cb"
#     include "include/slice.cb"
#     include "include/rad.cb"
#     include "include/bats.cb2"
#     include "include/bats.cb"

c ***** Input Variables ***** c
      real(kind=8)  uv995, tsurf, t995, q995, z995, zi, psurf
      integer n,i,j
c ***** Output Variables ***** c
      real(kind=8)  uv10, tau, lh, sh, dth, dqh, ustar, zo, qs

      real(kind=8)  facttq
      real(kind=8)  mmpd,wpm2

      do i=np1,nbmax
      do n=1,NNSG
        if (ocld2d(n,i,j).lt.0.5) then
          uv995 = sqrt(ubx3d(i,kx,j)**2+vbx3d(i,kx,j)**2)
          tsurf = tgb(i,j)-273.16
          t995 = tb3d(i,kx,j)-273.16
          q995 = qvb3d(i,kx,j)/(1.+qvb3d(i,kx,j))
          z995 = za(i,kx,j)
          zi = zpbl(i,j)
          psurf = (psb(i,j)+ptop)*10.
          CALL ZENGOCN(uv995,tsurf,t995,q995,z995,zi,psurf,qs
     &       , karman,g,r,cp,uv10,tau,lh,sh,dth,dqh,ustar,zo)
          tg1d(n,i) = tgb(i,j)
          tgb1d(n,i)= tgb(i,j)
          sent1d(n,i) = sh
          evpr1d(n,i) = lh/xlv
c Back out Drag Coefficient
          drag1d(n,i) = ustar**2*rhox2d(i,j)/uv995
          facttq = dlog(z995/2.)/dlog(z995/zo)
          u10m1d(n,i) = ubx3d(i,kx,j)*uv10/uv995
          v10m1d(n,i) = vbx3d(i,kx,j)*uv10/uv995
          t2m_1d(n,i) = t995 + 273.15 - dth*facttq

          if (mod(ntime+nint(dtmin*60.),kbats).eq.0
     &        .or. (jyear.eq.jyearr.and.ktau.eq.ktaur)) then
            if (jyear.eq.jyear0.and.ktau.le.1) then
              mmpd = 86400./dtbat
              wpm2 = 1./dtbat
            else if (jyear.eq.jyear0.and.
     &               dble(ktau*dtmin).le.batfrq*60.+0.01) then
              mmpd = 24./(batfrq-dtmin/60.)
              wpm2 = 1./((batfrq-dtmin/60.)*3600.)
            else
              mmpd = 24./batfrq
              wpm2 = 1./(batfrq*3600.)
            end if
            facttq = dlog(z995/2.)/dlog(z995/zo)
            q2m_1d(n,i) = q995 - dqh*facttq
            tgb2d(n,i,j)= tgb(i,j)
          end if
        end if
      end do
      end do

      return
      end subroutine zengocndrv

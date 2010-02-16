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
C
      subroutine zengocndrv(j)
C
      implicit none
C
#     include "parameter.inc"
#     include "include/param1.cb"
#     include "include/param2.cb"
#     include "include/param3.cb"
#     include "include/main.cb"
#     include "include/pbldim.cb"
#     include "include/pmoist.cb"
#     include "include/slice.cb"
#     include "include/bats2.cb"
#     include "include/bats.cb"
C
C     ***** Input Variables ***** c
C
      real(kind=8) uv995 , tsurf , t995 , q995 , z995 , zi , psurf
      integer n , i , j
C
C     ***** Output Variables ***** c
C
      real(kind=8) uv10 , tau , lh , sh , dth , dqh , ustar , zo , qs
C
      real(kind=8) facttq
C
      do i = np1 , nbmax
        do n = 1 , nnsg
          if ( ocld2d(n,i,j)<0.5 ) then
            uv995 = sqrt(ubx3d(i,kx,j)**2+vbx3d(i,kx,j)**2)
            tsurf = tgb(i,j) - 273.16
            t995 = tb3d(i,kx,j) - 273.16
            q995 = qvb3d(i,kx,j)/(1.+qvb3d(i,kx,j))
            z995 = za(i,kx,j)
            zi = zpbl(i,j)
            psurf = (psb(i,j)+ptop)*10.
            call zengocn(uv995,tsurf,t995,q995,z995,zi,psurf,qs,karman,
     &                   g,r,cp,uv10,tau,lh,sh,dth,dqh,ustar,zo)
            tg1d(n,i) = tgb(i,j)
            tgb1d(n,i) = tgb(i,j)
            sent1d(n,i) = sh
            evpr1d(n,i) = lh/xlv
C           Back out Drag Coefficient
            drag1d(n,i) = ustar**2*rhox2d(i,j)/uv995
            facttq = dlog(z995/2.)/dlog(z995/zo)
            u10m1d(n,i) = ubx3d(i,kx,j)*uv10/uv995
            v10m1d(n,i) = vbx3d(i,kx,j)*uv10/uv995
            t2m_1d(n,i) = t995 + 273.15 - dth*facttq
C
            if ( mod(ntime+nint(dtmin*60.),kbats)==0 .or. 
     &           (jyear==jyearr .and. ktau==ktaur) ) then
              facttq = dlog(z995/2.)/dlog(z995/zo)
              q2m_1d(n,i) = q995 - dqh*facttq
              tgb2d(n,i,j) = tgb(i,j)
            end if
          end if
        end do
      end do
C
      end subroutine zengocndrv

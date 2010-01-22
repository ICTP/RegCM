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

      subroutine vcover
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

c     provides leaf and stem area parameters;
c     depends on climate through subsoil temperatures.

      implicit none
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#     include "include/bats.cb"

      integer n,np

      do 90 np=np1,npts
      do n=1,NNSG
       if (ldoc1d(n,np).gt.0.5) then
        if (sigf(n,np).gt.0.001) then
         seasb(n,np)=dmax1(0.d0
     &                    ,1.-0.0016*dmax1(298.-tgb1d(n,np),0.d0)**2)
        endif
       endif
      enddo
90    continue

      do 91 np=np1,npts
      do n=1,NNSG
       if (ldoc1d(n,np).gt.0.5) then
        if (sigf(n,np).gt.0.001) then
         xlai(n,np)=xla(lveg(n,np))
         xlai(n,np)=xlai(n,np)
     &             +(xlai0(lveg(n,np))-xlai(n,np))*(1.-seasb(n,np))
         rlai(n,np)  = xlai(n,np) + sai(lveg(n,np))
         xlsai(n,np) = xlai(n,np) + sai(lveg(n,np))
         vegt(n,np)  = sigf(n,np) * xlsai(n,np)
        endif
       endif
      enddo
91    continue
c
      return
      end

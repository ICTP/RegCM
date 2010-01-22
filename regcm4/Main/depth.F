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

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine depth
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c          wt = fraction of vegetation covered by snow
c        sigf = fraction of veg. cover, excluding snow-covered veg.
c        scvk = fraction of soil covered by snow
c
c     scrat = snow depth (m) =  .001 snow depth (mm) / snow density
c     rhosw = ratio of snow density to density of h2o
c
c     height of vegetation assumed to be 10 x vegetation roughness ht
c     densi is defined the same here as in subr. albedo
c
c     wt now scaled so that runs betw. 0 & 1 and wt=0.5 when depth
c     of snow equals height of vegetation
c

      implicit none
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#     include "include/bats.cb"

      real(kind=8)  age
      integer n,np

      do 90 np=np1,npts
      do n=1,NNSG

       if( ldoc1d(n,np).gt.0.5 )then
        age       = (1.-1./(1.+sag1d(n,np)) )
        rhosw(n,np) = .10 * (1.+3.*age)
        densi(n,np) = .01 / (1.+3.*age)
        scrat(n,np) = scv1d(n,np)*densi(n,np)
        wt(n,np) = 1.0
        if (lveg(n,np) .gt. 0) then
           wt(n,np) = 0.1*scrat(n,np) / rough(lveg(n,np))
           wt(n,np) = wt(n,np) / (1.+wt(n,np))
        endif
        sigf(n,np) = (1.-wt(n,np)) * veg1d(n,np)
        scvk(n,np) = scrat(n,np) / (0.1+scrat(n,np))
       endif

      enddo
90    continue

      return
      end

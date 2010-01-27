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

      subroutine slice(j)

      implicit none

#     include "parameter.inc"
#     include "include/param1.cb"
#     include "include/param2.cb"
#     include "include/param3.cb"
#     include "include/main.cb"
#     include "include/cvaria.cb"
#     include "include/pbldim.cb"
#     include "include/pmoist.cb"
#     include "include/slice.cb"
#     include "include/bats2.cb"
#     include "include/bats.cb"

      integer j,n,np
      real(kind=8)  amxtem,sfac

c ********* For albedov
c **********            This is taken from subroutine interf so that
c **********            radiation can be called in tend (not vecbats).
      npts=nbmax
      do np=np1,npts
      do n=1,NNSG
        ldoc1d(n,np)=ocld2d(n,np,j)
        sice1d(n,np)=sice2d(n,np,j)
        tgb1d(n,np)=tgb2d(n,np,j)
        ssw1d(n,np)=ssw2d(n,np,j)
        lveg(n,np) = nint(veg2d1(n,np,j))
        amxtem=dmax1(298.-tgb1d(n,np),0.d0)
        sfac=1.-dmax1(0.d0,1.-0.0016*amxtem**2)
        if(lveg(n,np).eq.0)then
         veg1d(n,np)=0.
        else
         veg1d(n,np)=vegc(lveg(n,np))-seasf(lveg(n,np))*sfac
        end if
        ts1d(n,np)=thx3d(np,kx,j)-6.5e-3/g*(ht1(n,np,j)-ht(np,j))
        scv1d(n,np)=scv2d(n,np,j)
        sag1d(n,np)=sag2d(n,np,j)
      end do
      end do

      return
      end subroutine slice

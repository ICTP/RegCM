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

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
cqhy15/07/99 subroutine cumtran(il,jl,kx,ntr,dsigma,icumtop,chia,chib)
      subroutine cumtran
      implicit none
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#     include "include/trachem.cb"
#     include "include/main.cb"
#     include "include/mainchem.cb"
#ifdef MPP1
      include 'mpiregcm.h'
#endif
#     include "include/param3.cb"

      real(kind=8)  cumfrc,deltas,chiabar,chibbar
      integer i,j,k,n,kcumtop

chy the well-mixing only over cumulus cloud fraction (assuming 0.3)
c add the well-mixing fraction in dt step (5%) - only updraft
chy     cumfrc = 0.3
        cumfrc = 0.015

      do 100 n=1,ntr
#ifdef MPP1
      do 100 j=jbegin,jendx
#else
      do 100 j=2,jx-1
#endif
      do 100 i=2,ix-1
      if (icumtop(i,j).gt.0) then
         deltas=0.
         chiabar=0.
         chibbar=0.
         kcumtop=max0(icumtop(i,j),4)
         do 90 k=kcumtop,kx
         deltas=deltas+dsigma(k)
         chiabar=chiabar+chia(i,k,j,n)*dsigma(k)
         chibbar=chibbar+chib(i,k,j,n)*dsigma(k)
  90     continue
c?         do 95 k=icumtop(i,j),kx      ! yhuang, 12/98
cqhy
         do 95 k=kcumtop,kx
         chia(i,k,j,n)=chia(i,k,j,n)*(1.-cumfrc)+cumfrc*chiabar/deltas
         chib(i,k,j,n)=chib(i,k,j,n)*(1.-cumfrc)+cumfrc*chibbar/deltas
  95     continue
      endif
  100 continue
      return
      end

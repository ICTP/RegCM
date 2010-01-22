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

c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      subroutine frawat
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c  **** determines fraction of foliage covered by water fwet, and
c  **** the fraction of foliage that is dry transpiring leaf fdry.
c  note: their defns differ - fwet is the fraction of all veg surfaces
c  which are wet because stems can evaporate, fdry is the fraction
c  of lai which is dry because only leaves can transpire
c
c  ldew1d(np) is in kg/m**2/s
c  fwet   = ratio of dew to max value to 2/3 power
c           ( 2/3 power comes from deardorff (1978) )
c              ** keep fwet le 1.0 **
c  dewmxi = inverse of max allowed dew depth on leaf in mm
c
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
         fwet(n,np) = 0.
         if(ldew1d(n,np).gt.0.) then
          fwet(n,np) = ((dewmxi/vegt(n,np))*ldew1d(n,np))**.666666666666
          fwet(n,np) = dmin1(fwet(n,np),1.d0)
         endif
         fdry(n,np) = (1.-fwet(n,np))*xlai(n,np) / xlsai(n,np)
        endif
       endif
      enddo
90    continue

      return
      end

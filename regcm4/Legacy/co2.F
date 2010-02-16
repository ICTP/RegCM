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

      subroutine co2
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c     model for plant carbon uptake and dead carbon decomposition
c
c            c(4) = number of seconds in half-hour
c            cari = rate of carbon uptake by plants
c       pbp1d(n,np) = accumulated primary biomass
c      resp1d(n,np) = carbon in soil and in veg (kg c / m**2 / sec)
c              ra = leaf aerodynamic resistance factor
c             rap = leaf boundary layer resistance to co2
c           resps = soil respiration
c             rmp = physical mesophyllic resistance
c              rs = stomatal resistance
c             rsp = stomatal resistance to co2
c              rt = total mechanical resistance to co2
c
c     0.33 pptv co2 = 0.30 pptm dry matter in plants
c
c     decomposable soil carbon assumed incremented by npp-soil resp.
c
c     resistances for co2 are larger than those for h2o due to
c                    difference in molecular weight
c
      implicit none
#     include "parameter.inc"
#     include "include/bats.cb"
#     include "include/bats2.cb"

      real(kind=8)  rs
      real(kind=8)  ra
      common /ictp01/ rs(NNSG,nbmax), ra(NNSG,nbmax)
      real(kind=8)  cari(NNSG,nbmax)
      integer n,np
      real(kind=8)  rmp,rsp,rap,rt,carbon,resps

      rmp = 800.0
      do 90 np=np1,npts
      do n=1,NNSG
       if (ldoc1d(n,np).gt.0.5) then
        if (sigf(n,np).gt.0.001) then
         rsp = rs(n,np)*1.7
         rap = ra(n,np)*1.5
         rt=rsp+rap+rmp
         cari(n,np) = sigf(n,np) * xlsai(n,np) * fdry(n,np)
     a         * carbon( solis(np)*rlai(n,np), tlef1d(n,np), rt
     b                 , tg1d(n,np), xlai(n,np), xlsai(n,np) )
         pbp1d(n,np) = pbp1d(n,np) + cari(n,np)*c(4)
        end if
       end if
      end do
90    continue

      do 91 np=np1,npts
      do n=1,NNSG
       if (ldoc1d(n,np).gt.0.5) then
        if (sigf(n,np).gt.0.001) then
        if(pbp1d(n,np).lt.0) pbp1d(n,np) = 0.
         resps = 0.7e-7* resp1d(n,np)*dexp(0.1*(tg1d(n,np)-300.))
     a           *dmin1(1.d0,ssw1d(n,np)/(0.6*gwmx0(n,np)))
         resp1d(n,np) = resp1d(n,np) + (cari(n,np)-resps)*c(4)
         if (resp1d(n,np).lt.0.0) resp1d(n,np) = 0.
        end if
       end if
      end do
91    continue
      return
      end subroutine co2

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

      subroutine stomat
c:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c     gives leaf stomatal resistance from environmental parameters
c             under conditions of no moisture stress
c
c     standard lai from xla=max & xlai0=min lai (set in bconst)
c     seasb = fseas(tgb1d(n,np) (set in bndry) is a seasonal
c             factor for reduced winter lai and root water uptake
c        fc = light sensitivity for crops and grasses and has inverse
c             radiation units (m**2/watt)
c      rlai = sum of leaf and stem area indices
c     rmax0 = 5000. s/m (maximum resistance)
c     radu & radl = visible light intensity in upper & lower canopy
c     ft & fb = the fractional intercepted photo-active radiation
c             per unit (leaf & stem) area in the top (upper) and
c             bottom (lower) canopies, respectively
c     radfi = average of upper and lower canopy light factors
c        rs = stomatal resistance = min.res. * rad.factor * leaf factor
c      trup = transmission of the upper canopy, assumed to be the same
c             for the lower canopy,i.e., trup=dexp(-0.5*g*rlai/czen),
c             where g = attenuation factor
c
c     documented in ncar tech note, dickinson et al., 1986
c     improved stomatal shading, dickinson, nov 88.
c
      implicit none
#     include "parameter.inc"
#     include "include/bats.cb"
#     include "include/bats.cb2"

      real(kind=8)  rad(10), radd(10)

      real(kind=8)  vpdc
      Common/abb/ vpdc(NNSG,nbmax)
      real(kind=8)  rs
      real(kind=8)  ra
      common /ictp01/ rs(NNSG,nbmax), ra(NNSG,nbmax)

      real(kind=8)  trup(NNSG,nbmax), trupd(NNSG,nbmax)
      real(kind=8)  fsold(NNSG,nbmax), radf(NNSG,nbmax)
      real(kind=8)  fsol0(NNSG,nbmax),rmini(NNSG,nbmax)
      real(kind=8)  g,difzen,radfi,vpdf,seas
      integer n,il,np,ilmax
c
c ***** seasonal temperature factor
      real(kind=8)  x,fseas
      fseas(x)=dmax1(0.0d0,1.-0.0016*dmax1(298.-x,0.d0)**2)

c ***** g is average leaf crosssection per unit lai
c ***** difzen is ave of inverse of cos of angle of diffuse vis light
c ***** ilmax is number of canopy layers
c ***** czen is cosine solar zenith angle for incident light
c *****   (to spec from input data need a good treatment of diffuse rad)
c ***** trup is transmission of direct beam light in one canopy layer
c ***** trupd is transmission of diffuse light in one canopy layer
      g = 0.5
      difzen = 2.0
      ilmax = 4
c* delete fracd here to put in diffuse radiation from ccm
ccc   fracd = difrat         !  from shuttleworth mods #2

      do 91 np=np1,npts
      do n=1,NNSG
       if (ldoc1d(n,np).gt.0.5) then
        if (sigf(n,np).gt.0.001) then
c **********            zenith angle set in zenitm
         if(czen(np).gt.0.001) then
          trup(n,np) = dexp(-g*rlai(n,np)/(ilmax*czen(np)))
          trupd(n,np)= dexp(-difzen*g*rlai(n,np)/(ilmax))
          fsold(n,np) =     fracd(np) *solis(np)*fc(lveg(n,np))
          fsol0(n,np) = (1.-fracd(np))*solis(np)*fc(lveg(n,np))
          rmini(n,np) = rsmin(lveg(n,np))/rmax0
         end if
        end if
       end if
      end do
91    continue

      do 92 np=np1,npts
      do n=1,NNSG
       if (ldoc1d(n,np).gt.0.5) then
        if (sigf(n,np).gt.0.001) then
         if(czen(np).gt.0.001) then
          rad(1)  = (1.-trup(n,np) )*fsol0(n,np)*ilmax / rlai(n,np)
          radd(1) = (1.-trupd(n,np))*fsold(n,np)*ilmax / rlai(n,np)
          do 1 il = 2, ilmax
           rad(il)  = trup(n,np) *rad(il-1)
           radd(il) = trupd(n,np)*radd(il-1)
1         continue
          radfi = 0.
          do 2 il = 1, ilmax
           radfi=radfi+(rad(il)+radd(il)+rmini(n,np))/
     1     (1.+rad(il)+radd(il))
2         continue
          radf(n,np) = ilmax / radfi
         end if
        end if
       end if
      end do
92    continue

      do 93 np=np1,npts
      do n=1,NNSG
       if (ldoc1d(n,np).gt.0.5) then
        if (sigf(n,np).gt.0.001) then
         if(czen(np).gt.0.001) then
          vpdf = 1. / dmax1(0.3d0,1.d0- vpdc(n,np)*0.025)
          seas = 1. / (rmini(n,np)+fseas(tlef1d(n,np)))
          rs(n,np) = rsmin(lveg(n,np)) * radf(n,np) * seas * vpdf
          rs(n,np) = dmin1( rs(n,np), rmax0 )
         else
          rs(n,np) = rmax0
         end if
        end if
       end if
      end do
93    continue
      return
      end subroutine stomat

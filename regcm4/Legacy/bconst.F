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

      subroutine bconst
c::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
      implicit none

c     ***   defines constants for boundary subroutine
#     include "parameter.inc"

#     include "include/bats.cb"

      pi=3.141592

c     ****    specific heats per m**3            (joules/m**3/k)
      ch2o=4.186e6
      cice=0.45*ch2o
      cwi=1./0.45
c     ****  specific heat for snow after multiplied by snow density
      csnw=0.49*ch2o
      cws=1./0.49
c     ****
c     **** data for vegetation calculation
c     ***         maximum stomatl resistance    (s/m)
      rmax0 = 2.e4
c     ****                 maximum allowed dew(mm) and inverse (dewmxi)
      dewmx=0.1
      dewmxi=1.0/dewmx
c     ****
c     ****   trasmx0: maximum rate of transpiration with saturated soil
c     ****                    (kg/m**2/s)
      trsmx0=2.e-4
c     ****    drain is drainage out of 10m layer bottom (mm/s);
c     ****    drain is set fairly large to prevent swamping the soil
c     drain = 4.e-4
      drain = 1.e-4
c     ****        length of day in seconds
      tau1=8.64e4

c     ****    time steps c(6) and c(7) set in vecbats.

c     **** gravity (m/sec**2)
      c(54)=9.80616
c     **** gas constant for dry air (joules/kg/k)
      c(57)=2.8704e2
c     **** specific heat at constant pressure for dry air (cpd)
      c(58)=3.5*c(57)
c     **** ratio of rd/cpd is kappa  - comes into defn of potential temp
      c(68)=c(57)/c(58)
c     **** freezing point value
      c(67) = 273.16
c     **** si/non-d constants for tetens formula  - used in sbrt
c      'satur' which computes specific humidity qsat which is non-dim
c
      c(70)=21.874
      c(71)=7.66
      c(72)=17.269
      c(73)=35.86
      c(74)=6.11e2
      c(75)=0.622
      c(76)=0.378
c     ****  standard pressure  (pa)
      c(81)=1.013250e5
c     *****  stefans constant  (watts/m**2/k**4)
      c(83)=5.67e-8
c     ****  turbulent wind for stable conditions  (m/sec)
      c(90)=0.1
c     **** latent heats   (joules/kg)
      c(125)=2.51040e6
      c(127)=0.3336e6
      c(126)=c(125)+c(127)
c     ****    zlnd = soil roughness length
c     ****    zoce = ocean roughness length
c     ****    zsno = snow roughness length
c     ****    vonkar = von karman constant
      zlnd=0.01
      zoce=0.00040          ! changed to agree with pat's mods (4/92)
      zsno=0.00040
      vonkar=0.378
c     ****    csoilc = drag coefficient for soil under canopy
      csoilc = 4.e-3
      return
      end subroutine bconst

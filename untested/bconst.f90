!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
      subroutine bconst

      use regcm_param
      use mod_bats , only : c , ch2o , cice , csnw , csoilc , cwi ,     &
                    & cws , dewmx , dewmxi , drain , pi , rmax0 , tau1 ,&
                    & trsmx0 , vonkar , zlnd , zoce , zsno
      implicit none
!     ***   defines constants for boundary subroutine
 
      pi = 3.14159265358979323846
 
!     ****    specific heats per m**3            (joules/m**3/k)
      ch2o = 4.186E6
      cice = 0.45*ch2o
      cwi = 1./0.45
!     ****  specific heat for snow after multiplied by snow density
      csnw = 0.49*ch2o
      cws = 1./0.49
!     ****
!     **** data for vegetation calculation
!     ***         maximum stomatl resistance    (s/m)
      rmax0 = 2.E4
!     ****                 maximum allowed dew(mm) and inverse (dewmxi)
      dewmx = 0.1
      dewmxi = 1.0/dewmx
!     ****
!     ****   trasmx0: maximum rate of transpiration with saturated soil
!     ****                    (kg/m**2/s)
      trsmx0 = 2.E-4
!     ****    drain is drainage out of 10m layer bottom (mm/s);
!     ****    drain is set fairly large to prevent swamping the soil
!     drain = 4.e-4
      drain = 1.E-4
!     ****        length of day in seconds
      tau1 = 8.64E4
 
!     ****    time steps c(6) and c(7) set in vecbats.
 
!     **** gravity (m/sec**2)
      c(54) = 9.80616
!     **** gas constant for dry air (joules/kg/k)
      c(57) = 2.8704E2
!     **** specific heat at constant pressure for dry air (cpd)
      c(58) = 3.5*c(57)
!     **** ratio of rd/cpd is kappa  - comes into defn of potential temp
      c(68) = c(57)/c(58)
!     **** freezing point value
      c(67) = 273.16
!     **** si/non-d constants for tetens formula  - used in sbrt
!     'satur' which computes specific humidity qsat which is non-dim
!
      c(70) = 21.874
      c(71) = 7.66
      c(72) = 17.269
      c(73) = 35.86
      c(74) = 6.11E2
      c(75) = 0.622
      c(76) = 0.378
!     ****  standard pressure  (pa)
      c(81) = 1.013250E5
!     *****  stefans constant  (watts/m**2/k**4)
      c(83) = 5.67E-8
!     ****  turbulent wind for stable conditions  (m/sec)
      c(90) = 0.1
!     **** latent heats   (joules/kg)
      c(125) = 2.51040E6
      c(127) = 0.3336E6
      c(126) = c(125) + c(127)
!     ****    zlnd = soil roughness length
!     ****    zoce = ocean roughness length
!     ****    zsno = snow roughness length
!     ****    vonkar = von karman constant
      zlnd = 0.01
      zoce = 0.00040        ! changed to agree with pat's mods (4/92)
      zsno = 0.00040
      vonkar = 0.378
!     ****    csoilc = drag coefficient for soil under canopy
      csoilc = 4.E-3

      end subroutine bconst

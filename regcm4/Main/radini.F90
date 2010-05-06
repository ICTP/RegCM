!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
      subroutine radini

!-----------------------------------------------------------------------
!
! Initialize various constants for radiation scheme; note that
! the radiation scheme uses cgs units.
!
!---------------------------Code history--------------------------------
!
! Original version:  CCM1
! Standardized:      J. Rosinski, June 1992
! Reviewed:          J. Kiehl, B. Briegleb, August 1992
!
!-----------------------------------------------------------------------
!
      use mod_dynparam
      use mod_comozp
      use mod_crdcae
      use mod_tracer0
      use mod_tracer
      use mod_date
      use mod_message
      use mod_constants , only : gti , gtigts , cpd , stebol , sslp
      implicit none
!
!---------------------------Local variables-----------------------------
!
! iband  - H2O band index
! v0     - Volume of a gas at stp (m**3/kmol)
! p0     - Standard pressure (pascals)
! amd    - Effective molecular weight of dry air (kg/kmol)
!
! Local variables
!
      real(8) :: amd , p0 , v0
      integer :: iband
!
!     Set general radiation consts; convert to cgs units where
!     appropriate:
!
!IPCC
!1991-1995
!     co2vmr  =  3.55e-4
!     ch40 = 0.55241 * 1.714e-6
!     n2o0 = 1.51913 * 0.311e-6
!1961-1965
!     co2vmr  =  3.10e-4
!     ch40 = 0.55241 * 1.414e-6
!     n2o0 = 1.51913 * 0.287e-6
 
!     cfc110 = 4.69548 * 0.280e-9
!     cfc120 = 4.14307 * 0.503e-9
!     co2mmr = 1.51913 * co2vmr
      if ( lyear.ge.1750 .and. lyear.le.2100 ) then
        co2vmr = cgas(2,lyear)*1.E-6
        co2mmr = co2vmr*44.0/28.9644
        ch40 = cgas(3,lyear)*1.E-9*0.55241
        n2o0 = cgas(4,lyear)*1.E-9*1.51913
        cfc110 = cgas(5,lyear)*1.E-12*4.69548
        cfc120 = cgas(6,lyear)*1.E-12*4.14307
      else
        write (aline,*) '  Simulation date:  ' , lyear
        call say
        call fatal(__FILE__,__LINE__,                                   &
               &'CONCENTRATION VALUES OUTSIDE OF DATE RANGE (1750-2100)'&
              & )
      end if
!     print*,'IN RADINI (TOP)'
!     print*,'  co2vmr= ',co2vmr
!     print*,'  co2mmr= ',co2mmr
!     print*,'  ch40  = ',ch40
!     print*,'  n2o0  = ',n2o0
!     print*,'  cfc110= ',cfc110
!     print*,'  cfc120= ',cfc120
!
!     Coefficients for h2o emissivity and absorptivity.
!
      do iband = 1 , 4
        c1(iband) = coefe(3,iband)/coefe(2,iband)
        c2(iband) = coefb(3,iband)/coefb(2,iband)
        c3(iband) = coefb(4,iband)/coefb(3,iband)
        c4(iband) = coefd(3,iband)/coefd(2,iband)
        c5(iband) = coefd(4,iband)/coefd(3,iband)
        c6(iband) = coefa(3,iband)/coefa(2,iband)
        c7(iband) = coefc(3,iband)/coefc(2,iband)
      end do
      c8 = coeff(3,1)/coeff(2,1)
      c9 = coeff(3,2)/coeff(2,2)
      c10 = coeff(4,1)/coeff(3,1)
      c11 = coeff(4,2)/coeff(3,2)
      c12 = coeff(5,1)/coeff(4,1)
      c13 = coeff(5,2)/coeff(4,2)
      c14 = coeff(6,1)/coeff(5,1)
      c15 = coeff(6,2)/coeff(5,2)
      c16 = coefj(3,1)/coefj(2,1)
      c17 = coefk(3,1)/coefk(2,1)
      c18 = coefi(3,1)/coefi(2,1)
      c19 = coefi(3,2)/coefi(2,2)
      c20 = coefi(4,1)/coefi(3,1)
      c21 = coefi(4,2)/coefi(3,2)
      c22 = coefi(5,1)/coefi(4,1)
      c23 = coefi(5,2)/coefi(4,2)
      c24 = coefi(6,1)/coefi(5,1)
      c25 = coefi(6,2)/coefi(5,2)
      c26 = coefj(3,2)/coefj(2,2)
      c27 = coefk(3,2)/coefk(2,2)
      c28 = .5
      c29 = .002053
      c30 = .1
      c31 = 3.0E-5
      cfa1 = .61
!
!     Initialize further longwave constants referring to far wing
!     correction; R&D refers to:
!
!     Ramanathan, V. and  P.Downey, 1986: A Nonisothermal
!     Emissivity and Absorptivity Formulation for Water Vapor
!     Journal of Geophysical Research, vol. 91., D8, pp 8649-8666
!
      fwcoef = .1           ! See eq(33) R&D
      fwc1 = .30            ! See eq(33) R&D
      fwc2 = 4.5            ! See eq(33) and eq(34) in R&D
      fc1 = 2.6             ! See eq(34) R&D
!
!     Initialize ozone data.
!
      v0 = 22.4136          ! Volume of a gas at stp (m**3/kmol)
      p0 = 0.1*sslp         ! Standard pressure (pascals)
      amd = 28.9644         ! Molecular weight of dry air (kg/kmol)
!
!     Constants for ozone path integrals (multiplication by 100 for unit
!     conversion to cgs from mks):
!
      cplos = v0/(amd*gti)*100.0
      cplol = v0/(amd*gti*p0)*0.5*100.0
!
      end subroutine radini

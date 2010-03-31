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
!
      module mod_aero_param

      use mod_regcm_param , only : nveg

      implicit none
!
      integer , parameter :: isize = 12

!     Basic aerosol density (ACE-2 ?) in kg/m^3
      real(8) , parameter :: rhop = 2650.0

!     Dynamic Viscosity Parameters
      real(8) , parameter :: a1 = 145.8
      real(8) , parameter :: a2 = 1.5
      real(8) , parameter :: a3 = 110.4

!     Molecular Free Path calculation parameters
      real(8) , parameter :: c1 = 6.54E-8
      real(8) , parameter :: c2 = 1.818E-5
      real(8) , parameter :: c3 = 1.013E5
      real(8) , parameter :: c4 = 293.15

!     Cunningham slip correction factor parameters
      real(8) , parameter :: aa1 = 1.257
      real(8) , parameter :: aa2 = 0.4
      real(8) , parameter :: aa3 = 1.1

!     Aerosol Dry radius
      real(8) , dimension(2,isize) :: aerosize
      real(8) , dimension(isize) :: avesize
!     Stokes parameter
      real(8) , dimension(nveg) :: aest
      real(8) , dimension(nveg) :: arye
!
! DATA SECTION
!
      data aerosize/1.0E-08 , 2.0E-08 , 2.0E-08 , 4.0E-08 , 4.0E-08 ,   &
         & 8.0E-08 , 8.0E-08 , 1.6E-07 , 1.6E-07 , 3.2E-07 , 3.2E-07 ,  &
         & 6.4E-07 , 6.4E-07 , 1.28E-06 , 1.28E-06 , 2.56E-06 ,         &
         & 2.56E-06 , 5.12E-06 , 5.12E-06 , 10.4E-06 , 10.24E-06 ,      &
         & 20.48E-06 , 20.48E-06 , 40.6E-06/

      data aest/0.80 , 0.80 , 0.8 , 0.8 , 1.2 , 1.20 , 2.00 , 1.5 ,     &
         & 1.5 , 2.0 , 15.0 , 15.0 , 1.5 , 1.5 , 1.5 , 15.0 , 1.20 ,    &
         & 1.2 , 1.2 , 1.2/

      data arye/0.5 , 5.0 , 0.5 , 5.0 , 1.0 , 1.0 , 0.0001 , 5.0 ,      &
         & 10.0 , 10.0 , 0.0001 , 0.0001 , 0.56 , 0.56 , 0.56 , 0.56 ,  &
         & 0.56 , 0.56 , 0.56 , 0.56/

      end module mod_aero_param

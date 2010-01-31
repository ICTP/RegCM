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
 
      subroutine inicon
      use constants
      implicit none

!     Preset constants in mo_constants.
!
!     Method:
!
!     *inicon* is called from *setdyn*.
!
!     Authors:
!
!     M. Jarraud, ECMWF, December 1982, original source
!     L. Kornblueh, MPI, May 1998, f90 rewrite
!     U. Schulzweida, MPI, May 1998, f90 rewrite
!     H.-S. Bauer, MPI, Jul 1998, changed
!     A. Rhodin, MPI, Jan 1999, subroutine inicon put into module
!     mo_constants
!     for more details see file AUTHORS
!
 
!--   1. Preset constants
 
      gti = 9.80665
      cpd = 1005.46
      cpv = 1869.46
      rdti = 287.05
      rvti = 461.51
 
      rcpd = 1./cpd
      vtmpc1 = rvti/rdti - 1.
      vtmpc2 = cpv/cpd - 1.
 
      rhoh2o = 1000.
      alv = 2.5008E6
      als = 2.8345E6
      alf = als - alv
 
      tmelt = 273.16
 
      c1es = 610.78
      c2es = c1es*rdti/rvti
      c3les = 17.269
      c3ies = 21.875
      c4les = 35.86
      c4ies = 7.66
      c5les = c3les*(tmelt-c4les)
      c5ies = c3ies*(tmelt-c4ies)
      c5alvcp = c5les*alv/cpd
      c5alscp = c5ies*als/cpd
      alvdcp = alv/cpd
      alsdcp = als/cpd
 
      end subroutine inicon

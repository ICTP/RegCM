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

      subroutine rot2nrot(lams,phis,pollon,pollat,lam,phi)
      implicit none
!
! Dummy arguments
!
      real(4) :: lam , lams , phi , phis , pollat , pollon
      intent (in) lams , phis , pollat , pollon
      intent (out) lam , phi
!
! Local variables
!
      real(4) :: arg , d2r , plam , pphi , r2d , zarg1 , zarg2 ,        &
               & zcospol , zlampol , zlams , zphis , zsinpol
!
!----------------------------------------------------------------------------
!     Purpose:
!     Adaption of the DWD-Functions to convert rotated pole coordinates
!     (PHIS,LAMS) into geogrphic coordinates (PHI,LAM). The location of
!     the rotated pole is passed trough POLLON and POLLAT. POLLON and
!     POLLAT give the origin of the rotated grid. The
!     first four arguments are input, the last two are output. All
!     angles are in degrees (north>0, east>0)
!     History:
!     05/90   D.MAJEWSKI (DWD)
!     03/93   D.BRESCH (ETHZ)
!     09/96   D.LUETHI (ETHZ)
 
      r2d = 45./atan(1.)
      d2r = atan(1.)/45.
 
      plam = pollon + 180.
      pphi = 90. - pollat
 
      if ( plam>180. ) plam = plam - 360.
      zsinpol = sin(d2r*pphi)
      zcospol = cos(d2r*pphi)
      zlampol = d2r*plam
 
      zphis = 2*atan(exp(d2r*phis)) - atan(1.)*2.
      zlams = lams
      if ( zlams>180.0 ) zlams = zlams - 360.0
      zlams = d2r*zlams
 
!     first, the conversion of PHIS to PHI:
      arg = zcospol*cos(zphis)*cos(zlams) + zsinpol*sin(zphis)
      phi = r2d*asin(arg)
 
!     follows conversion of LAMS to LAM:
      zarg1 = sin(zlampol)                                              &
            & *(-zsinpol*cos(zlams)*cos(zphis)+zcospol*sin(zphis))      &
            & - cos(zlampol)*sin(zlams)*cos(zphis)
      zarg2 = cos(zlampol)                                              &
            & *(-zsinpol*cos(zlams)*cos(zphis)+zcospol*sin(zphis))      &
            & + sin(zlampol)*sin(zlams)*cos(zphis)
      if ( abs(zarg2)>=1.E-30 ) then
        lam = r2d*atan2(zarg1,zarg2)
      else if ( abs(zarg1)<1.E-30 ) then
        lam = 0.0
      else if ( zarg1>0. ) then
        lam = 90.0
      else
        lam = -90.0
      end if
 
      end subroutine rot2nrot

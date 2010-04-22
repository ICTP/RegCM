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

      subroutine nrot2rot(lam,phi,pollon,pollat,lams,phis)
      use mod_constants , only : r2d => raddeg
      use mod_constants , only : d2r => degrad
      implicit none
!
! Dummy arguments
!
      real(8) :: lam , lams , phi , phis , pollat , pollon
      intent (in) lam , phi , pollat , pollon
      intent (out) lams
      intent (inout) phis
!
! Local variables
!
      real(8) :: plam , pphi , zarg , zarg1 , zarg2 , zcospol , zlam ,  &
             &   zlampol , zphi , zsinpol
!
! ----------------------------------------------------------------------
!     Purpose:
!     Adaption of the DWD-Functions to convert real geographical
!     coordinates (PHI,LAM) into coordinates in the rotated system
!     (PHIS,LAMS). The rotated pole is passed trough POLLON and POLLAT.
!     POLLON and POLLAT give the origin of the new rotated grid. The
!     first four arguments are input, the second two are output. All
!     angles are in degrees (north>0, east>0)
!     History:
!     05/90   D.MAJEWSKI (DWD), G. DE MORSIER (SMA)
!     03/93   D.BRESCH (ETHZ)
!     11/97   D.LUETHI (ETHZ)
 
      plam = pollon + 180.
      pphi = 90. - pollat
 
      if ( plam>180. ) plam = plam - 360.
 
      zsinpol = sin(d2r*pphi)
      zcospol = cos(d2r*pphi)
      zlampol = d2r*plam
 
!     first, the conversion of PHI to PHIS:
      zphi = d2r*phi
      zlam = lam
      if ( zlam>180.0 ) zlam = zlam - 360.0
      zlam = d2r*zlam
      zarg = zcospol*cos(zphi)*cos(zlam-zlampol) + zsinpol*sin(zphi)
      phis = asin(zarg)
      phis = log(tan(phis/2.+atan(1.)))*r2d
 
!     now, the conversion for LAMS follws:
      zphi = d2r*phi
      zlam = lam
      if ( zlam>180.0 ) zlam = zlam - 360.0
      zlam = d2r*zlam
      zarg1 = -sin(zlam-zlampol)*cos(zphi)
      zarg2 = -zsinpol*cos(zphi)*cos(zlam-zlampol) + zcospol*sin(zphi)
      if ( abs(zarg2)>=1.E-30 ) then
        lams = r2d*atan2(zarg1,zarg2)
      else if ( abs(zarg1)<1.E-30 ) then
        lams = 0.0
      else if ( zarg1>0. ) then
        lams = 90.0
      else
        lams = -90.0
      end if
 
      end subroutine nrot2rot

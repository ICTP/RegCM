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

      module mod_projections

      use mod_constants , only : earthrad , eomeg2
      use mod_constants , only : mathpi , degrad , raddeg

      implicit none

      contains

      subroutine lambrt(xlon,xlat,smap,coriol,iy,jx,clon,clat,ds,idot,  &
                      & xn,truelatl,truelath)

      implicit none
!
! Dummy arguments
!
      real(4) :: clat , clon , ds , truelath , truelatl , xn
      integer :: idot , iy , jx
      real(4) , dimension(iy,jx) :: coriol , smap , xlat , xlon
      intent (in) clat , clon , ds , idot , iy , jx , truelath ,        &
                & truelatl
      intent (out) coriol , smap , xlat , xlon , xn
!
! Local variables
!
      real(8) :: cell , cell1 , cell2 , cntri , cntrj , flp , flpp ,    &
               & pole , psi1 , psix , psx , r , xsign , truelat1 ,      &
               & truelat2 , x , xcntr , y , ycntr
      integer :: i , j
!
!     CLAT IS THE CENTRAL LATITUDE OF THE COARSE DOMAIN.
!     CLON IS THE CENTRAL LONGITUDE OF THE COARSE DOMAIN.
!
!     THIS ROUTINE CALCULATES MESO MAP(LAT,LONG,CORIOLIS,MAP SCALE)
!     FOR LAMBERT CONFORMAL PROJECTION
!
!     IY IS THE I DIMENSION FOR THIS DOMAIN.
!     JX IS THE J DIMENSION FOR THIS DOMAIN.
!     IDOT IS ICROSS ( = 1)  OR IDOT ( = 0).
!
!---------------------------------------------------------------------
!
!
!     XN IS CONE FACTOR FOR THE PROJECTION (FROM PROGRAM TERRAIN).
!     PSI1 IS THE COLATITUDE OF TRUELAT 1, IN RADIANS.
!     PI IS PI.
!
!---------------------------------------------------------------------
!
      if ( clat<0. ) then
        xsign = -1.       ! SOUTH HEMESPHERE
      else
        xsign = 1.        ! NORTH HEMESPHERE
      end if
      pole = xsign*90.0
 
      truelat1 = truelath
      truelat2 = truelatl
      if ( abs(truelat1-truelat2)>1.E-1 ) then
        xn = (dlog10(cos(truelat1*degrad))-dlog10(cos(truelat2*degrad)))&
           & /(dlog10(tan((45.0-xsign*truelat1/2.0)*degrad))            &
           & -dlog10(tan((45.0-xsign*truelat2/2.0)*degrad)))
      else
        xn = xsign*sin(truelat1*degrad)
      end if
!     XN=0.716
 
      psi1 = 90. - xsign*truelat1
      if ( clat<0. ) psi1 = -psi1
!
      psi1 = psi1*degrad
      cntrj = (jx+idot)/2.
      cntri = (iy+idot)/2.
!
      psx = (pole-clat)*degrad
      cell = earthrad*sin(psi1)/xn
      write (*,*) 'PSX,PSI1 = ' , psx , psi1
      cell2 = (tan(psx/2.))/(tan(psi1/2.))
      r = cell*(cell2)**xn
      xcntr = 0.
      ycntr = -r
!
      do j = 1 , jx
        x = xcntr + (j-cntrj)*ds
        do i = 1 , iy
          y = ycntr + (i-cntri)*ds
          r = sqrt(x*x+y*y)
          if ( y==0. ) then
            if ( x>=0. ) then
              flp = 90.*degrad
            else
              flp = -90.*degrad
            end if
          else if ( clat<0.0 ) then
            flp = atan2(x,y)
          else
            flp = atan2(x,-y)
          end if
          flpp = flp/xn/degrad + clon
!         IF (FLPP.GT.180.0) FLPP = FLPP-360.0
!         IF (FLPP.LT.-180.0) FLPP = FLPP+360.0
          xlon(i,j) = flpp
          if ( clat<0.0 ) r = -r
          cell = r*xn/(earthrad*sin(psi1))
          cell1 = tan(psi1/2.)*cell**(1./xn)
          cell2 = atan(cell1)
          psx = 2.*cell2/degrad
          xlat(i,j) = pole - psx
          psix = psx*degrad
          smap(i,j) = (sin(psi1)/sin(psix))                             &
                    & *((tan(psix/2.)/tan(psi1/2.))**xn)
        end do
      end do
      if ( idot==1 ) then
        do i = 1 , iy
          do j = 1 , jx
            coriol(i,j) = eomeg2*sin(xlat(i,j)*degrad)
          end do
        end do
      end if
      end subroutine lambrt

      subroutine mappol(xlon,xlat,xmap,coriol,iy,jx,clon,clat,delx,idot)

      implicit none
!
! Dummy arguments
!
      real(4) :: clat , clon , delx
      integer :: idot , iy , jx
      real(4) , dimension(iy,jx) :: coriol , xlat , xlon , xmap
      intent (in) clat , clon , delx , idot , iy , jx
      intent (out) coriol , xlat , xlon , xmap
!
! Local variables
!
      real(4) :: cell , cell1 , cell2 , cntri , cntrj , flp , flpp ,    &
               & pole , psi1 , psix , psx , r , x , xcntr , xn , y ,    &
               & ycntr
      integer :: i , ii1 , j , jj1
!
!     XN IS CONE FACTOR FOR THE PROJECTION (FROM PROGRAM TERRAIN).
!     PSI1 IS THE COLATITUDE OF TRUELAT 1, IN RADIANS.
!     PI IS PI.
!
!---------------------------------------------------------------------
 
!     COMPUTE LATS, LONS, AND MAP-SCALE FACTORS FOR
!     LAMBERT CONFORMAL MAP CENTERED AT CLON,CLAT. TRUE AT 30.N AND
 
!     60.N. IY IS NUMBER OF N-S POINTS.  JX IS NUMBER OF E-W POINTS.
!     CLON, CLAT IS LAT, LON OF CENTER OF GRID (DEGREES EAST, NORTH).
!     DELX IS GRID SPACING IN METERS.
!     ALWAYS FOR CROSS GRID.
 
!
      xn = 1.0
!
      pole = 90.
      psi1 = 30.
      psi1 = psi1*degrad
      if ( clat<0.0 ) then
        pole = -90.0
        psi1 = -30.
        psi1 = psi1*degrad
      end if
      cntrj = float(jx+idot)/2.
      cntri = float(iy+idot)/2.
!
      psx = (pole-clat)*degrad
      cell = earthrad*sin(psx)/xn
      cell2 = (1.+cos(psi1))/(1.+cos(psx))
      r = cell*(cell2)**xn
      xcntr = 0.
      ycntr = -r
!
      ii1 = iy
      jj1 = jx
      do i = 1 , ii1
        y = ycntr + (i-cntri)*delx
        do j = 1 , jj1
          x = xcntr + (j-cntrj)*delx
          r = sqrt(x*x+y*y)
          if ( y==0. ) then
            if ( x>=0. ) then
              flp = 90.*degrad
            else
              flp = -90.*degrad
            end if
          else if ( clat<0.0 ) then
            flp = atan2(x,y)
          else
            flp = atan2(x,-y)
          end if
          flpp = flp/xn/degrad + clon
!         IF (FLPP.GT.180.0) FLPP = FLPP-360.0
!         IF (FLPP.LT.-180.0) FLPP = FLPP+360.0
          xlon(i,j) = flpp
          if ( clat<0.0 ) r = -r
          cell = r/earthrad
          cell1 = cell/(1.0+cos(psi1))
          cell2 = atan(cell1)
          psx = 2.*cell2/degrad
          xlat(i,j) = pole - psx
          psix = psx*degrad
          xmap(i,j) = ((1.0+cos(psi1))/(1.0+cos(psix)))**xn
        end do
      end do
 
      if ( idot==1 ) then
        do i = 1 , iy
          do j = 1 , jx
            coriol(i,j) = eomeg2*sin(xlat(i,j)*degrad)
          end do
        end do
      end if
      end subroutine mappol

      subroutine normer(xlon,xlat,xmap,coriol,iy,jx,clon,clat,delx,idot)

      implicit none
!
! Dummy arguments
!
      real(4) :: clat , clon , delx
      integer :: idot , iy , jx
      real(4) , dimension(iy,jx) :: coriol , xlat , xlon , xmap
      intent (in) clat , clon , delx , idot , iy , jx
      intent (out) coriol , xlat , xlon , xmap
!
! Local variables
!
      real(4) :: c2 , cell , cntri , cntrj , deglat , phi1 , phictr ,   &
               & pole , x , xcntr , y , ycntr
      integer :: i , ii1 , j , jj1
!
!     COMPUTE LATS, LONS, AND MAP-SCALE FACTORS FOR
!     LAMBERT CONFORMAL MAP CENTERED AT CLON,CLAT. TRUE AT 30.N AND
!     60.N. IY IS NUMBER OF N-S POINTS.  JX IS NUMBER OF E-W POINTS.
!     CLON, CLAT IS LAT, LON OF CENTER OF GRID (DEGREES EAST, NORTH).
!     DELX IS GRID SPACING IN METERS.
!     ALWAYS FOR CROSS GRID.
 
      pole = 90.
      cntrj = (jx+idot)/2.
      cntri = (iy+idot)/2.
      if ( clat<0.0 ) pole = -90.0
!
!     FOR MERCATOR PROJECTION TRUE AT PHI1
!
      phi1 = 0.
      phi1 = phi1*degrad
      c2 = earthrad*cos(phi1)
      xcntr = 0.
      phictr = clat*degrad
      cell = cos(phictr)/(1.+sin(phictr))
      ycntr = -c2*log(cell)
!
      ii1 = iy
      jj1 = jx
      do i = 1 , ii1
        y = ycntr + (i-cntri)*delx
        do j = 1 , jj1
          x = xcntr + (j-cntrj)*delx
!
!         CALCULATIONS FOR MERCATOR
!
          xlon(i,j) = clon + ((x-xcntr)/c2)/degrad
          cell = exp(y/c2)
          xlat(i,j) = 2.*(atan(cell)/degrad) - 90.
          deglat = xlat(i,j)*degrad
          xmap(i,j) = cos(phi1)/cos(deglat)
        end do
      end do
 
      if ( idot==1 ) then
        do i = 1 , iy
          do j = 1 , jx
            coriol(i,j) = eomeg2*sin(xlat(i,j)*degrad)
          end do
        end do
      end if
 
      end subroutine normer

      subroutine rotmer(xlon,xlat,xmap,coriol,iy,jx,clon,clat,pollon,   &
                      & pollat,ds,idot)

      implicit none
!
! Dummy arguments
!
      real(4) :: clat , clon , ds , pollat , pollon
      integer :: idot , iy , jx
      real(4) , dimension(iy,jx) :: coriol , xlat , xlon , xmap
      intent (in) clat , clon , ds , idot , iy , jx
      intent (out) coriol , xlat , xlon , xmap
!
! Local variables
!
      real(8) :: cntri , cntrj , ddeg , fai , xoff , yoff
      real(8) :: xr , yr , x , y , plat , plon
      integer :: i , j
!
!---------------------------------------------------------------------
!     COMPUTE LATS, LONS, MAP-SCALE FACTORS, AND CORIOLIS PARAMETER FOR
!     ROTATED POLE MAP CENTERED AT CLON,CLAT. ORIGIN OF ROTATED GRID IS
!     GIVEN BY POLLON AND POLLAT.
!     IMX,JMX,KSV,KSH AND LSIG2 MUST CORRESPOND TO IY,JMAX,KSTRPV,KSTRPH
!     AND NVERT2 IN THE MASTER INPUT FILE; OTHERWISE THE PROGRAM WILL
!     ABORT:LMX MUST BE THE MAXIMUM NUMBER OF LEVELS
!     (PRESSURE--OR--SIGMA) IMAXN AND IMXC ARE NESTED AND NON-EXPANDED
!     GRID DIMENSIONS. IMXC IS EQUAL TO IMX IF YOU ARE NOT USING THE
!     EXPANDED GRID. SAME FOR J.
!
!-----CENTER OF GRID
!
      cntrj = (jx+idot)/2.
      cntri = (iy+idot)/2.

!     GRID SPACING IN DEGREES

      ddeg = ds*raddeg/earthrad
      plat = pollat
      plon = pollon
!
      xoff = clon - pollon
      yoff = clat - pollat
!
!-----CALCULATE X AND Y POSITIONS OF GRID
!
      do i = 1 , iy
        do j = 1 , jx
          xr = xoff + (j-cntrj)*ddeg
          yr = yoff + (i-cntri)*ddeg
!
!-----NOW CALCULATE LAT AND LON OF THIS POINT
!-----    ROTATE COORDINATES BACK TO NONRATED POLE
!
          call rot2nrot(xr,yr,plon,plat,x,y)
!
          xlon(i,j) = x
          xlat(i,j) = y
          fai = degrad*yr
          xmap(i,j) = 1.0/cos(fai)
        end do
      end do
 
      if ( idot==1 ) then
        do i = 1 , iy
          do j = 1 , jx
            coriol(i,j) = eomeg2*sin(xlat(i,j)*degrad)
          end do
        end do
      end if
 
      end subroutine rotmer

      subroutine rot2nrot(lams,phis,pollon,pollat,lam,phi)

      implicit none
!
! Dummy arguments
!
      real(8) :: lam , lams , phi , phis , pollat , pollon
      intent (in) lams , phis , pollat , pollon
      intent (out) lam , phi
!
! Local variables
!
      real(8) :: arg , plam , pphi , zarg1 , zarg2 , zcospol , zlampol ,&
              &  zlams , zphis , zsinpol
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
 
      plam = pollon + 180.
      pphi = 90. - pollat
 
      if ( plam>180. ) plam = plam - 360.
      zsinpol = sin(degrad*pphi)
      zcospol = cos(degrad*pphi)
      zlampol = degrad*plam
 
      zphis = 2*atan(exp(degrad*phis)) - atan(1.)*2.
      zlams = lams
      if ( zlams>180.0 ) zlams = zlams - 360.0
      zlams = degrad*zlams
 
!     first, the conversion of PHIS to PHI:
      arg = zcospol*cos(zphis)*cos(zlams) + zsinpol*sin(zphis)
      phi = raddeg*asin(arg)
 
!     follows conversion of LAMS to LAM:
      zarg1 = sin(zlampol)                                              &
            & *(-zsinpol*cos(zlams)*cos(zphis)+zcospol*sin(zphis))      &
            & - cos(zlampol)*sin(zlams)*cos(zphis)
      zarg2 = cos(zlampol)                                              &
            & *(-zsinpol*cos(zlams)*cos(zphis)+zcospol*sin(zphis))      &
            & + sin(zlampol)*sin(zlams)*cos(zphis)
      if ( abs(zarg2)>=1.E-37 ) then
        lam = raddeg*atan2(zarg1,zarg2)
      else if ( abs(zarg1)<1.E-37 ) then
        lam = 0.0
      else if ( zarg1>0. ) then
        lam = 90.0
      else
        lam = -90.0
      end if

      end subroutine rot2nrot

      subroutine nrot2rot(lam,phi,pollon,pollat,lams,phis)

      implicit none
!
! Dummy arguments
!
      real(8) :: lam , lams , phi , phis , pollat , pollon
      intent (in) lam , phi , pollat , pollon
      intent (out) lams , phis
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
 
      zsinpol = sin(degrad*pphi)
      zcospol = cos(degrad*pphi)
      zlampol = degrad*plam
 
!     first, the conversion of PHI to PHIS:
      zphi = degrad*phi
      zlam = lam
      if ( zlam>180.0 ) zlam = zlam - 360.0
      zlam = degrad*zlam
      zarg = zcospol*cos(zphi)*cos(zlam-zlampol) + zsinpol*sin(zphi)
      phis = asin(zarg)
      phis = log(tan(phis/2.+atan(1.)))*raddeg
 
!     now, the conversion for LAMS follws:
      zphi = degrad*phi
      zlam = lam
      if ( zlam>180.0 ) zlam = zlam - 360.0
      zlam = degrad*zlam
      zarg1 = -sin(zlam-zlampol)*cos(zphi)
      zarg2 = -zsinpol*cos(zphi)*cos(zlam-zlampol) + zcospol*sin(zphi)
      if ( abs(zarg2)>=1.E-37 ) then
        lams = raddeg*atan2(zarg1,zarg2)
      else if ( abs(zarg1)<1.E-37 ) then
        lams = 0.0
      else if ( zarg1>0. ) then
        lams = 90.0
      else
        lams = -90.0
      end if
 
      end subroutine nrot2rot

      end module mod_projections

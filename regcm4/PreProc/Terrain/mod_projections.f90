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
      real(8) :: ala1 , alo1 , x , y , hemi , rebydx , reflon
      real(8) :: scale_top , rsw , polei , polej
      real(8) :: r2 , gi2 , arcc , cntrj , cntri
      integer :: i , j
!
      cntrj = float(jx+idot)/2.
      cntri = float(iy+idot)/2.
      rebydx = earthrad / delx
      if (clat > 0.0) then
        hemi = 1.0
      else
        hemi = -1.0
      end if
      reflon = clon + 90.0
      ala1 = clat*degrad
      alo1 = (clon-reflon)*degrad
      scale_top = 1. + hemi * sin(ala1)
      rsw = rebydx*cos(ala1)*scale_top/(1.0+hemi*sin(ala1))
      polei = cntri - hemi * rsw * sin(alo1)
      polej = cntrj - rsw * cos(alo1)
!
      do i = 1 , iy
        y = (i - polei) * hemi
        do j = 1 , jx
          x = j - polej
          r2 = x**2 + y**2
          if (abs(r2) < 1e-30) then
            xlat(i,j) = hemi*90.0
            xlon(i,j) = reflon
          else
            gi2 = (rebydx * scale_top)**2.
            xlat(i,j) = raddeg * hemi * asin((gi2-r2)/(gi2+r2))
            arcc = acos(x/sqrt(r2))
            if (y > 0) then
              xlon(i,j) = reflon + raddeg * arcc
            else
              xlon(i,j) = reflon - raddeg * arcc
            end if
          end if
          xmap(i,j) = scale_top/(1. + hemi * sin(xlat(i,j)*degrad))
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

      subroutine xyobsll(iy,jx,iproj,clat,clon,plat,plon,truelath)
      use mod_maps
      use mod_block
      use mod_constants , only : raddeg , degrad , erkm
      implicit none
!
! Dummy arguments
!
      real(4) :: clat , clon , plat , plon , truelath
      character(6) :: iproj
      integer :: iy , jx
      intent (in) clat , clon , iproj , iy , jx , truelath
!
! Local variables
!
      real(8) :: c2 , cell , cell2 , cntri , cntrj , flp , flpp ,       &
               & phi1 , phic , phir , phix , pole , psi1 , psx , r ,    &
               & xcntr , xlonx ,  xrot , ycntr
      real(8) :: xnr , ynr , xr , yr , pla , plo, cla , clo
      integer :: ie , je , ii , im , ilen
!
      pla = plat
      plo = plon
      cla = clat
      clo = clon
      ilen = iy*jx + 2
      ie = iy - 1
      je = jx - 1
      psi1 = 1.0D+36
      pole = 90.0D+00
      if ( iproj=='LAMCON' ) psi1 = 90. - truelath
      if ( iproj=='POLSTR' ) psi1 = 30.0
      psi1 = psi1/raddeg

!-----psi1 is colatitude of lat where cone or plane intersects earth

      cell = 0.0
      cell2 = 0.0
      ycntr = 0.0

      if ( cla<0. ) then
        if ( truelath>0. ) then
          psi1 = -(90.-truelath)
        else
          psi1 = -(90.+truelath)
        end if
        pole = -90.
        psi1 = psi1/raddeg
      end if
      if ( iproj=='LAMCON' .or. iproj=='POLSTR' ) then
        psx = (pole-cla)/raddeg
        if ( iproj=='LAMCON' ) then
          cell = erkm*sin(psi1)/xn
          cell2 = (tan(psx/2.))/(tan(psi1/2.))
        else if ( iproj=='POLSTR' ) then
          cell = erkm*sin(psx)/xn
          cell2 = (1.+cos(psi1))/(1.+cos(psx))
        else
        end if
        r = cell*(cell2)**xn
        xcntr = 0.0
        ycntr = -r
      end if
      cntrj = float(je)/2.
      cntri = float(ie)/2.
!
!-----grid incoming data.  grdltmn=minimum latitude of incoming data.
!-----grdlnmn=minimum longitude of incoming data.
!
      do ii = 1 , nobs
        im = ii - 1
        if ( iproj=='LAMCON' .or. iproj=='POLSTR' ) then
          xrot = clo + 90./xn
          phix = yobs(ii)
          xlonx = xobs(ii)
          flpp = (xlonx-xrot)/raddeg
          flp = xn*flpp
          psx = (pole-phix)/raddeg
          if ( iproj=='LAMCON' ) then
            cell = erkm*sin(psi1)/xn
            cell2 = (tan(psx/2.))/(tan(psi1/2.))
          else if ( iproj=='POLSTR' ) then
            cell = erkm*sin(psx)/xn
            cell2 = (1.+cos(psi1))/(1.+cos(psx))
          else
          end if
          r = cell*(cell2)**xn
          xobs(ii) = (r*cos(flp)-xcntr)*1000.
          yobs(ii) = (r*sin(flp)-ycntr)*1000.
          if ( cla<0.0 ) xobs(ii) = -xobs(ii)
        end if
        if ( iproj=='NORMER' ) then
          phi1 = 0.0   ! plat/raddeg
          phir = yobs(ii)/raddeg
          phic = cla/raddeg
          c2 = erkm*cos(phi1)
          cell = cos(phir)/(1.0+sin(phir))
          cell2 = cos(phic)/(1.0+sin(phic))
          ycntr = -c2*log(cell2)
          xobs(ii) = (c2*(xobs(ii)-clo)/raddeg)*1000.
          yobs(ii) = (-c2*log(cell)-ycntr)*1000.
        end if
        if ( iproj=='ROTMER' ) then
          xcntr = plo - clo
          ycntr = pla - cla
          xnr = xobs(ii)
          ynr = yobs(ii)
          call nrot2rot(xnr,ynr,plo,pla,xr,yr)
          xobs(ii) = erkm*degrad*(xcntr+xr)*1000.
          yobs(ii) = erkm*degrad*(ycntr+yr)*1000.
        end if
        ht(ii) = ht(ii)/100.
        ht2(ii) = ht2(ii)/100000.
      end do

      end subroutine xyobsll
!
      end module mod_projections

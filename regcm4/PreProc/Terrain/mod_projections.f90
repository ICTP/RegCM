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
      use mod_constants , only : degrad , raddeg

      real(8) , private :: truelat1 , truelat2 , stdlon , tl1r , ctl1r
      real(8) , private :: conefac , rsw , rebydx , hemi
      real(8) , private :: reflon , dlon , scale_top
      real(8) , private :: polei , polej
      real(8) , private :: polelon , polelat , xoff , yoff
      real(8) , private :: zsinpol , zcospol , zlampol

      contains

      subroutine setup_lcc(clat,clon,ci,cj,ds,slon,trlat1,trlat2)
        implicit none
        real(4) , intent(in) :: ci , cj , slon , clat , clon , ds ,     &
                      & trlat1 , trlat2
        real(8) :: arg , deltalon1 , tl2r
!
        stdlon = slon
        truelat1 = trlat1
        truelat2 = trlat2

        tl1r = truelat1*degrad
        tl2r = truelat2*degrad

        if (truelat1 > 0.0) then
          hemi = 1.0
        else
          hemi = -1.0
        end if
        rebydx = earthrad / ds

        if ( abs(truelat1-truelat2) > 0.1) then
          conefac = log10(cos(tl1r)) - log10(cos(tl2r))
          conefac = conefac/(log10(tan((45.0-abs(truelat1)/2.0)*degrad))&
                    -log10(tan((45.0 - abs(truelat2)/2.0) * degrad)))   
        else
          conefac = sin(abs(tl1r))
        end if
        deltalon1 = clon - stdlon
        if (deltalon1 >  180.0) deltalon1 = deltalon1 - 360.
        if (deltalon1 < -180.0) deltalon1 = deltalon1 + 360.

        ctl1r = cos(tl1r)
        rsw = rebydx * ctl1r/conefac * (tan((90.*hemi-clat)*degrad/2.) /&
                  tan((90.*hemi-truelat1)*degrad/2.))**conefac
        arg = conefac*(deltalon1*degrad)
        polei = hemi*ci - hemi * rsw * sin(arg)
        polej = hemi*cj + rsw * cos(arg)

      end subroutine setup_lcc

      subroutine ijll_lc(i,j,lat,lon)
        real(4) , intent(in) :: i , j
        real(4) , intent(out) :: lat , lon
        real(8) :: chi1 , chi2 , chi
        real(8) :: inew , jnew , xx , yy , r2 , r

        chi1 = (90. - hemi*truelat1)*degrad
        chi2 = (90. - hemi*truelat2)*degrad
        inew = hemi * i
        jnew = hemi * j
        xx = inew - polei
        yy = polej - jnew
        r2 = (xx*xx + yy*yy)
        r = sqrt(r2)/rebydx
        if (abs(r2) < 1e-30) then
          lat = hemi * 90.
          lon = stdlon
        else
          lon = stdlon + raddeg * atan2(hemi*xx,yy)/conefac
          lon = mod(lon+360.0, 360.0)
          if (abs(chi1-chi2) < 1e-30) then
            chi = 2.0*atan((r/tan(chi1) )**(1./conefac) * tan(chi1*0.5))
          else
            chi = 2.0*atan((r*conefac/sin(chi1))**(1./conefac) *        &
                & tan(chi1*0.5))
          end if
          lat = (90.0-chi*raddeg)*hemi
        end if
        if (lon >  180.0) lon = lon - 360.0
        if (lon < -180.0) lon = lon + 360.0
      end subroutine ijll_lc

      subroutine llij_lc(lat,lon,i,j)
        implicit none
        real(4) , intent(in) :: lat , lon
        real(4) , intent(out) :: i , j
        real(8) :: arg , deltalon , rm

        deltalon = lon - stdlon
        if (deltalon > +180.) deltalon = deltalon - 360.0
        if (deltalon < -180.) deltalon = deltalon + 360.0

        rm = rebydx * ctl1r/conefac * (tan((90.*hemi-lat)*degrad/2.) /  &
                tan((90.*hemi-truelat1)*degrad/2.))**conefac
        arg = conefac*(deltalon*degrad)
        i = polei + hemi * rm * sin(arg)
        j = polej - rm * cos(arg)
        i = hemi * i
        j = hemi * j
      end subroutine llij_lc

      subroutine uvrot_lc(lon, alpha)
        implicit none
        real(4) , intent(in) :: lon
        real(4) , intent(out) :: alpha
        real(4) :: deltalon

        deltalon = lon - stdlon
        if (deltalon > +180.) deltalon = deltalon - 360.0
        if (deltalon < -180.) deltalon = deltalon + 360.0
        alpha = deltalon*degrad*conefac

      end subroutine uvrot_lc

      subroutine setup_plr(clat,clon,ci,cj,ds,slon)
        implicit none
        real(4) , intent(in) :: clat , clon , cj , ci , ds , slon
        real(8) :: ala1 , alo1

        stdlon = slon
        if (clat > 0.0) then
          hemi = 1.0
        else
          hemi = -1.0
        end if
        rebydx = earthrad / ds
        reflon = stdlon + 90.0
        ala1 = clat*degrad
        alo1 = (clon-reflon)*degrad
        scale_top = 1. + hemi * sin(ala1)
        rsw = rebydx*cos(ala1)*scale_top/(1.0+hemi*sin(ala1))
        polei = ci - rsw * cos(alo1)
        polej = cj - hemi * rsw * sin(alo1)
      end subroutine setup_plr

      subroutine llij_ps(lat,lon,i,j)
        implicit none
        real(4) , intent(in) :: lat , lon
        real(4) , intent(out) :: i , j
        real(8) :: ala , alo , rm

        ala = lat * degrad
        rm = rebydx * cos(ala) * scale_top/(1.0 + hemi * sin(ala))
        alo = (lon - reflon) * degrad
        i = polei + rm * cos(alo)
        j = polej + hemi * rm * sin(alo)
      end subroutine llij_ps

      subroutine ijll_ps(i,j,lat,lon)
        implicit none
        real(4) , intent(in) :: i , j
        real(4) , intent(out) :: lat , lon
        real(8) :: xx , yy , r2 , gi2 , arcc

        xx = i - polei
        yy = (j - polej) * hemi
        r2 = xx**2 + yy**2
        if (abs(r2) < 1e-30) then
          lat = hemi*90.0
          lon = reflon
        else
          gi2 = (rebydx * scale_top)**2.0
          lat = raddeg * hemi * asin((gi2-r2)/(gi2+r2))
          arcc = acos(xx/sqrt(r2))
          if (yy > 0) then
            lon = reflon + raddeg * arcc
          else
            lon = reflon - raddeg * arcc
          end if
        end if
        if (lon >  180.) lon = lon - 360.0
        if (lon < -180.) lon = lon + 360.0
      end subroutine ijll_ps
     
      subroutine uvrot_ps(lon, alpha)
        implicit none
        real(4) , intent(in) :: lon
        real(4) , intent(out) :: alpha
        real(4) :: deltalon

        deltalon = lon - stdlon
        if (deltalon > +180.) deltalon = deltalon - 360.0
        if (deltalon < -180.) deltalon = deltalon + 360.0
        alpha = deltalon*degrad*hemi

      end subroutine uvrot_ps

      subroutine setup_mrc(clat,clon,ci,cj,ds)
        implicit none
        real(4) , intent(in) :: clat , clon , cj , ci , ds
        real(8) :: clain

        stdlon = clon
        clain = cos(clat*degrad)
        dlon = ds / (earthrad * clain)
        rsw = 0.0
        if (abs(clat) > 1e-30) then
          rsw = (log(tan(0.5*((clat+90.)*degrad))))/dlon
        end if
        polei = ci
        polej = cj
      end subroutine setup_mrc

      subroutine llij_mc(lat,lon,i,j)
        implicit none
        real(4) , intent(in) :: lat , lon
        real(4) , intent(out) :: i , j
        real(8) :: deltalon

        deltalon = lon - stdlon
        if (deltalon > +180.) deltalon = deltalon - 360.0
        if (deltalon < -180.) deltalon = deltalon + 360.0
        i = polei + (deltalon/(dlon*raddeg))
        j = polej + (log(tan(0.5*((lat + 90.)*degrad)))) / dlon - rsw
      end subroutine llij_mc

      subroutine ijll_mc(i,j,lat,lon)
        implicit none
        real(4) , intent(in) :: i , j
        real(4) , intent(out) :: lat , lon

        lat = 2.0*atan(exp(dlon*(rsw + j-polej)))*raddeg - 90.
        lon = (i-polei)*dlon*raddeg + stdlon
        if (lon >  180.) lon = lon - 360.0
        if (lon < -180.) lon = lon + 360.0
      end subroutine ijll_mc
     
      subroutine setup_rmc(clat,clon,ci,cj,ds,plon,plat)
        implicit none
        real(4) , intent(in) :: clat , clon , cj , ci , ds , plon , plat
        real(8) :: plam , pphi

        polelon = plon
        polelat = plat
        dlon = ds*raddeg/earthrad
        xoff = clon - plon
        yoff = clat - plat
        polei = ci
        polej = cj
        pphi = 90. - plat
        plam = plon + 180.
        if ( plam>180. ) plam = plam - 360.
        zsinpol = sin(degrad*pphi)
        zcospol = cos(degrad*pphi)
        zlampol = degrad*plam
      end subroutine setup_rmc

      subroutine llij_rc(lat,lon,i,j)
        implicit none
        real(4) , intent(in) :: lat , lon
        real(4) , intent(out) :: i , j
        real(8) :: zarg , zarg1 , zarg2 , zlam , zphi
        real(8) :: lams , phis
     
        zphi = degrad*lat
        zlam = lon
        if ( zlam>180.0 ) zlam = zlam - 360.0
        zlam = degrad*zlam

        zarg = zcospol*cos(zphi)*cos(zlam-zlampol) + zsinpol*sin(zphi)
        phis = asin(zarg)
        phis = log(tan(phis/2.+atan(1.)))*raddeg
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
        i = polei + (lams-xoff)/dlon
        j = polej + (phis-yoff)/dlon
      end subroutine llij_rc

      subroutine ijll_rc(i,j,lat,lon)
        implicit none
        real(4) , intent(in) :: i , j
        real(4) , intent(out) :: lat , lon
        real(8) :: xr , yr , arg , zarg1 , zarg2

        xr = xoff + (j-polej)*dlon
        if ( xr>180.0 ) xr = xr - 360.0
        xr = degrad*xr
        yr = yoff + (i-polej)*dlon
        yr = 2*atan(exp(degrad*yr)) - atan(1.)*2.
     
        arg = zcospol*cos(yr)*cos(xr) + zsinpol*sin(yr)
        lat = raddeg*asin(arg)
        zarg1 = sin(zlampol)*(-zsinpol*cos(xr)*cos(yr)+zcospol*sin(yr))-&
                cos(zlampol)*sin(xr)*cos(yr)
        zarg2 = cos(zlampol)*(-zsinpol*cos(xr)*cos(yr)+zcospol*sin(yr))+&
                sin(zlampol)*sin(xr)*cos(yr)
        if ( abs(zarg2)>=1.E-37 ) then
          lon = raddeg*atan2(zarg1,zarg2)
        else if ( abs(zarg1)<1.E-37 ) then
          lon = 0.0
        else if ( zarg1>0. ) then
          lon = 90.0
        else
          lon = -90.0
        end if
        if (lon >  180.) lon = lon - 360.0
        if (lon < -180.) lon = lon + 360.0
      end subroutine ijll_rc
     
      function rounder(value,ltop)
        implicit none
        real(4) , intent(in) :: value
        logical, intent(in) :: ltop
        real(4) :: rounder
        integer :: tmpval
        if (ltop) then
          tmpval = ceiling(value*100.0)
        else
          tmpval = floor(value*100.0)
        end if
        rounder = real(tmpval)/100.0
      end function rounder

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
      real(8) :: clain , cntri , cntrj , deglat , dlon , rsw , x , y
      integer :: i , j
!
!     COMPUTE LATS, LONS, AND MAP-SCALE FACTORS FOR MERCATOR
 
      cntrj = (jx+idot)/2.
      cntri = (iy+idot)/2.
!
      clain = cos(clat*degrad)
      dlon = delx / (earthrad * clain)
      rsw = 0.0
      if (abs(clat) > 1e-30) then
        rsw = (log(tan(0.5*((clat+90.)*degrad))))/dlon
      end if
!
      do i = 1 , iy
        y = i
        do j = 1 , jx
          x = j
          xlon(i,j) = (j-cntrj)*dlon*raddeg + clon
          xlat(i,j) = 2.0*atan(exp(dlon*(rsw + i-cntri)))*raddeg - 90.
          deglat = xlat(i,j)*degrad
          xmap(i,j) = 1.0/cos(deglat)
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

      subroutine xyobsll(iy,jx,iproj,clat,clon,plat,plon,ds,truelatl,   &
                       & truelath)
      use mod_maps
      use mod_block
      implicit none
!
! Dummy arguments
!
      real(4) :: clat , clon , plat , plon , truelath , truelatl , ds
      character(6) :: iproj
      integer :: iy , jx
      intent (in) iy , jx
      intent (in) clat , clon , iproj , truelath , truelatl , ds
      real(4) :: cntrj , cntri , alon , alat
      integer :: ii
!
! Local variables
!
      cntri = float(iy)/2.
      cntrj = float(jx)/2.
!
      if ( iproj=='LAMCON' ) then
       call setup_lcc(clat,clon,cntrj,cntri,dsinm,clon,truelatl,       &
                     & truelath)
      else if (iproj == 'POLSTR') then
        call setup_plr(clat,clon,cntrj,cntri,dsinm,clon)
      else if (iproj == 'NORMER') then
        call setup_mrc(clat,clon,cntrj,cntri,dsinm)
      else if (iproj == 'ROTMER') then
        call setup_rmc(clat,clon,cntrj,cntri,dsinm,plon,plat)
      else
        write (6,*) 'Unknown projection in xyobsll'
        stop
      end if

      do ii = 1 , nobs
        alon = xobs(ii)
        alat = yobs(ii)
        if ( iproj=='LAMCON' ) then
          call llij_lc(alat,alon,xobs(ii),yobs(ii))
        else if ( iproj=='POLSTR' ) then
          call llij_ps(alat,alon,xobs(ii),yobs(ii))
        else if ( iproj=='NORMER' ) then
          call llij_mc(alat,alon,xobs(ii),yobs(ii))
        else if ( iproj=='ROTMER' ) then
          call llij_rc(alat,alon,xobs(ii),yobs(ii))
        end if
        xobs(ii) = (xobs(ii)-cntrj) * dsinm
        yobs(ii) = (yobs(ii)-cntri) * dsinm
        ht(ii) = ht(ii)/100.
        ht2(ii) = ht2(ii)/100000.
      end do
      print *, minval(xobs), ' : ' , maxval(xobs)
      print *, minval(yobs), ' : ' , maxval(yobs)

      end subroutine xyobsll
!
      end module mod_projections

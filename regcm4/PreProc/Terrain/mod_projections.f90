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

      use mod_constants , only : earthrad
      use mod_constants , only : degrad , raddeg

      real(8) :: conefac
      real(8) , private :: stdlon
      real(8) , private :: truelat1 , truelat2 , tl1r , tl2r , ctl1r
      real(8) , private :: colat1 , colat2 , nfac
      real(8) , private :: rsw , rebydx , hemi
      real(8) , private :: reflon , dlon , scale_top
      real(8) , private :: polei , polej
      real(8) , private :: polelon , polelat , xoff , yoff
      real(8) , private :: zsinpol , zcospol , zlampol , zphipol
      logical , private :: lamtan

      contains

      subroutine setup_lcc(clat,clon,ci,cj,ds,slon,trlat1,trlat2)
        implicit none
        real(4) , intent(in) :: ci , cj , slon , clat , clon , ds , &
                             &  trlat1 , trlat2
        real(8) :: arg , deltalon1
!
        stdlon = slon
        truelat1 = trlat1
        truelat2 = trlat2

        tl1r = truelat1*degrad
        tl2r = truelat2*degrad
        colat1  = degrad*(90.0 - truelat1)
        colat2  = degrad*(90.0 - truelat2)

        nfac = (log(sin(colat1)) - log(sin(colat2))) &
                / (log(tan(colat1/2.0)) - log(tan(colat2/2.0)))

        if (truelat1 > 0.0) then
          hemi = 1.0
        else
          hemi = -1.0
        end if
        rebydx = earthrad / ds

        if ( abs(truelat1-truelat2) > 0.1) then
          conefac = log10(cos(tl1r)) - log10(cos(tl2r))
          conefac = conefac / &
                  & (log10(tan((45.0-abs(truelat1)/2.0)*degrad)) - &
                  & log10(tan((45.0 - abs(truelat2)/2.0) * degrad)))   
          lamtan = .false.
        else
          conefac = sin(abs(tl1r))
          lamtan = .true.
        end if
        deltalon1 = clon - stdlon
        if (deltalon1 >  180.0) deltalon1 = deltalon1 - 360.
        if (deltalon1 < -180.0) deltalon1 = deltalon1 + 360.

        ctl1r = cos(tl1r)
        rsw = rebydx * ctl1r/conefac * (tan((90.*hemi-clat)*degrad/2.) /&
            &     tan((90.*hemi-truelat1)*degrad/2.))**conefac
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
            chi = 2.0*atan((r/tan(chi1))**(1./conefac)*tan(chi1*0.5))
          else
            chi = 2.0*atan((r*conefac/sin(chi1))**(1./conefac) * &
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

        rm = rebydx * ctl1r/conefac * (tan((90.*hemi-lat)*degrad/2.) / &
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

        deltalon = stdlon - lon
        if (deltalon > +180.) deltalon = deltalon - 360.0
        if (deltalon < -180.) deltalon = deltalon + 360.0
        alpha = deltalon*degrad*conefac

      end subroutine uvrot_lc

      subroutine mapfac_lc(lat, xmap)
        implicit none
        real(4) , intent(in) :: lat
        real(4) , intent(out) :: xmap
        real(8) :: colat

        colat = degrad*(90.0-lat)
        if (.not. lamtan) then
          xmap = sin(colat2)/sin(colat) * &
               & (tan(colat/2.0)/tan(colat2/2.0))**nfac
        else
          xmap = sin(colat1)/sin(colat) * &
               & (tan(colat/2.0)/tan(colat1/2.0))**cos(colat1)
        endif
      end subroutine mapfac_lc

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
        real(8) :: ala , alo , rm , deltalon

        deltalon = lon - reflon
        if (deltalon > +180.) deltalon = deltalon - 360.0
        if (deltalon < -180.) deltalon = deltalon + 360.0
        alo = deltalon * degrad
        ala = lat * degrad

        rm = rebydx * cos(ala) * scale_top/(1.0 + hemi * sin(ala))
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
 
      subroutine mapfac_ps(lat, xmap)
        implicit none
        real(4) , intent(in) :: lat
        real(4) , intent(out) :: xmap
        xmap = scale_top/(1. + hemi * sin(lat*degrad))
      end subroutine mapfac_ps

      subroutine uvrot_ps(lon, alpha)
        implicit none
        real(4) , intent(in) :: lon
        real(4) , intent(out) :: alpha
        real(4) :: deltalon

        deltalon = stdlon - lon
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
 
      subroutine mapfac_mc(lat, xmap)
        implicit none
        real(4) , intent(in) :: lat
        real(4) , intent(out) :: xmap
        xmap = 1.0/cos(lat*degrad)
      end subroutine mapfac_mc

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
        if (polelat > 0.0) then
          pphi = 90. - plat
          plam = plon + 180.
        else
          pphi = 90. + plat
          plam = plon
        end if
        if ( plam>180. ) plam = plam - 360.
        zlampol = degrad*plam
        zphipol = degrad*pphi
        zsinpol = sin(zphipol)
        zcospol = cos(zphipol)
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

        xr = xoff + (i-polei)*dlon
        if ( xr>180.0 ) xr = xr - 360.0
        xr = degrad*xr
        yr = yoff + (j-polej)*dlon
        yr = 2*atan(exp(degrad*yr)) - atan(1.)*2.
 
        arg = zcospol*cos(yr)*cos(xr) + zsinpol*sin(yr)
        lat = raddeg*asin(arg)
        zarg1 = sin(zlampol)*(-zsinpol*cos(xr)*cos(yr)+ &
              & zcospol*sin(yr))-cos(zlampol)*sin(xr)*cos(yr)
        zarg2 = cos(zlampol)*(-zsinpol*cos(xr)*cos(yr)+ &
              & zcospol*sin(yr))+sin(zlampol)*sin(xr)*cos(yr)
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
 
      subroutine uvrot_rc(lat, lon, alpha)
        implicit none
        real(4) , intent(in) :: lon , lat
        real(4) , intent(out) :: alpha
        real(8) :: zphi , zrla , zrlap , zarg1 , zarg2 , znorm
        zphi = lat*degrad
        zrla = lon*degrad
        if (lat > 89.999999) zrla = 0.0
        zrlap = zlampol - zrla
        zarg1 = zcospol*sin(zrlap)
        zarg2 = zsinpol*cos(zphi) - zcospol*sin(zphi)*cos(zrlap)
        znorm = 1.0/sqrt(zarg1**2+zarg2**2)
        alpha = acos(zarg2*znorm)
      end subroutine uvrot_rc

      subroutine mapfac_rc(ir, xmap)
        implicit none
        real(4) , intent(in) :: ir
        real(4) , intent(out) :: xmap
        real(8) :: yr
        yr = yoff + (ir-polej)*dlon
        xmap = 1.0/cos(yr*degrad)
      end subroutine mapfac_rc

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
!
      end module mod_projections

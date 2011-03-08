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
        real(8) , intent(in) :: ci , cj , slon , clat , clon , ds , &
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
          hemi = 1.0D0
        else
          hemi = -1.0D0
        end if
        rebydx = earthrad / ds

        if ( dabs(truelat1-truelat2) > 0.1D0) then
          conefac = dlog10(dcos(tl1r)) - dlog10(dcos(tl2r))
          conefac = conefac / &
             & (dlog10(dtan((45.0D0-dabs(truelat1)/2.0D0)*degrad)) - &
             &  dlog10(dtan((45.0D0 - dabs(truelat2)/2.0D0) * degrad)))
          lamtan = .false.
        else
          conefac = dsin(dabs(tl1r))
          lamtan = .true.
        end if
        deltalon1 = clon - stdlon
        if (deltalon1 >  180.0D0) deltalon1 = deltalon1 - 360.0D0
        if (deltalon1 < -180.0D0) deltalon1 = deltalon1 + 360.0D0

        ctl1r = dcos(tl1r)
        rsw = rebydx * ctl1r/conefac * &
                (dtan((90.0D0*hemi-clat)*degrad/2.0D0) /&
            &    dtan((90.0D0*hemi-truelat1)*degrad/2.0D0))**conefac
        arg = conefac*(deltalon1*degrad)
        polei = hemi*ci - hemi * rsw * dsin(arg)
        polej = hemi*cj + rsw * dcos(arg)
      end subroutine setup_lcc

      subroutine ijll_lc(i,j,lat,lon)
        real(4) , intent(in) :: i , j
        real(4) , intent(out) :: lat , lon
        real(8) :: chi1 , chi2 , chi
        real(8) :: inew , jnew , xx , yy , r2 , r

        chi1 = (90.0D0 - hemi*truelat1)*degrad
        chi2 = (90.0D0 - hemi*truelat2)*degrad
        inew = hemi * dble(i)
        jnew = hemi * dble(j)
        xx = inew - polei
        yy = polej - jnew
        r2 = (xx*xx + yy*yy)
        r = dsqrt(r2)/rebydx
        if (dabs(r2) < 1D-30) then
          lat = real(hemi * 90.0D0)
          lon = real(stdlon)
        else
          lon = real(stdlon + raddeg * datan2(hemi*xx,yy)/conefac)
          lon = mod(lon+360.0, 360.0)
          if (dabs(chi1-chi2) < 1D-30) then
            chi = 2.0D0*datan((r/dtan(chi1))** &
                          (1.0D0/conefac)*tan(chi1*0.5D0))
          else
            chi = 2.0D0*datan((r*conefac/dsin(chi1))** &
                     (1.0D0/conefac)*dtan(chi1*0.5D0))
          end if
          lat = real((90.0D0-chi*raddeg)*hemi)
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
        i = real(hemi*(polei + hemi * rm * dsin(arg)))
        j = real(hemi*(polej - rm * cos(arg)))
      end subroutine llij_lc

      subroutine uvrot_lc(lon, alpha)
        implicit none
        real(4) , intent(in) :: lon
        real(4) , intent(out) :: alpha
        real(8) :: deltalon

        deltalon = stdlon - dble(lon)
        if (deltalon > +180.D0) deltalon = deltalon - 360.0D0
        if (deltalon < -180.D0) deltalon = deltalon + 360.0D0
        alpha = real(deltalon*degrad*conefac)
      end subroutine uvrot_lc

      subroutine mapfac_lc(lat, xmap)
        implicit none
        real(4) , intent(in) :: lat
        real(4) , intent(out) :: xmap
        real(8) :: colat

        colat = degrad*(90.0D0-dble(lat))
        if (.not. lamtan) then
          xmap = real(dsin(colat2)/dsin(colat) * &
               & (dtan(colat/2.0D0)/dtan(colat2/2.0D0))**nfac)
        else
          xmap = real(dsin(colat1)/dsin(colat) * &
               & (dtan(colat/2.0D0)/dtan(colat1/2.0D0))**dcos(colat1))
        endif
      end subroutine mapfac_lc

      subroutine setup_plr(clat,clon,ci,cj,ds,slon)
        implicit none
        real(8) , intent(in) :: clat , clon , cj , ci , ds , slon
        real(8) :: ala1 , alo1

        stdlon = slon
        if (clat > 0.0D0) then
          hemi = 1.0D0
        else
          hemi = -1.0D0
        end if
        rebydx = earthrad / ds
        reflon = stdlon + 90.0D0
        ala1 = clat*degrad
        alo1 = (clon-reflon)*degrad
        scale_top = 1.0D0 + hemi * dsin(ala1)
        rsw = rebydx*dcos(ala1)*scale_top/(1.0D0+hemi*dsin(ala1))
        polei = ci - rsw * dcos(alo1)
        polej = cj - hemi * rsw * dsin(alo1)
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

        rm = rebydx * dcos(ala) * scale_top/(1.0D0 + hemi * dsin(ala))
        i = real(polei + rm * dcos(alo))
        j = real(polej + hemi * rm * dsin(alo))
      end subroutine llij_ps

      subroutine ijll_ps(i,j,lat,lon)
        implicit none
        real(4) , intent(in) :: i , j
        real(4) , intent(out) :: lat , lon
        real(8) :: xx , yy , r2 , gi2 , arcc

        xx = i - polei
        yy = (j - polej) * hemi
        r2 = xx**2 + yy**2
        if (abs(r2) < 1D-30) then
          lat = real(hemi*90.0D0)
          lon = real(reflon)
        else
          gi2 = (rebydx * scale_top)**2.0D0
          lat = real(raddeg * hemi * dasin((gi2-r2)/(gi2+r2)))
          arcc = dacos(xx/dsqrt(r2))
          if (yy > 0) then
            lon = real(reflon + raddeg * arcc)
          else
            lon = real(reflon - raddeg * arcc)
          end if
        end if
        if (lon >  180.) lon = lon - 360.0
        if (lon < -180.) lon = lon + 360.0
      end subroutine ijll_ps
 
      subroutine mapfac_ps(lat, xmap)
        implicit none
        real(4) , intent(in) :: lat
        real(4) , intent(out) :: xmap
        xmap = real(scale_top/(1.0D0 + hemi * dsin(dble(lat)*degrad)))
      end subroutine mapfac_ps

      subroutine uvrot_ps(lon, alpha)
        implicit none
        real(4) , intent(in) :: lon
        real(4) , intent(out) :: alpha
        real(8) :: deltalon

        deltalon = stdlon - dble(lon)
        if (deltalon > +180.0D0) deltalon = deltalon - 360.0D0
        if (deltalon < -180.0D0) deltalon = deltalon + 360.0D0
        alpha = real(deltalon*degrad*hemi)
      end subroutine uvrot_ps

      subroutine setup_mrc(clat,clon,ci,cj,ds)
        implicit none
        real(8) , intent(in) :: clat , clon , cj , ci , ds
        real(8) :: clain

        stdlon = clon
        clain = cos(clat*degrad)
        dlon = ds / (earthrad * clain)
        rsw = 0.0D0
        if (dabs(clat) > 1D-30) then
          rsw = (dlog(dtan(0.5D0*((clat+90.0D0)*degrad))))/dlon
        end if
        polei = ci
        polej = cj
      end subroutine setup_mrc

      subroutine llij_mc(lat,lon,i,j)
        implicit none
        real(4) , intent(in) :: lat , lon
        real(4) , intent(out) :: i , j
        real(8) :: deltalon

        deltalon = dble(lon) - stdlon
        if (deltalon > +180.D0) deltalon = deltalon - 360.0D0
        if (deltalon < -180.D0) deltalon = deltalon + 360.0D0
        i = real(polei + (deltalon/(dlon*raddeg)))
        j = real(polej + (dlog(dtan(0.5D0* &
                  ((dble(lat) + 90.0D0)*degrad)))) / dlon - rsw)
      end subroutine llij_mc

      subroutine ijll_mc(i,j,lat,lon)
        implicit none
        real(4) , intent(in) :: i , j
        real(4) , intent(out) :: lat , lon

        lat = real(2.0D0*datan(dexp(dlon*(rsw + &
                                          j-polej)))*raddeg - 90.0D0)
        lon = real((dble(i)-polei)*dlon*raddeg + stdlon)
        if (lon >  180.) lon = lon - 360.0
        if (lon < -180.) lon = lon + 360.0
      end subroutine ijll_mc
 
      subroutine mapfac_mc(lat, xmap)
        implicit none
        real(4) , intent(in) :: lat
        real(4) , intent(out) :: xmap
        xmap = real(1.0D0/dcos(dble(lat)*degrad))
      end subroutine mapfac_mc

      subroutine setup_rmc(clat,clon,ci,cj,ds,plon,plat)
        implicit none
        real(8) , intent(in) :: clat , clon , cj , ci , ds , plon , plat
        real(8) :: plam , pphi

        polelon = plon
        polelat = plat
        dlon = ds*raddeg/earthrad
        xoff = clon - plon
        yoff = clat - plat
        polei = ci
        polej = cj
        pphi = 90.0D0 - plat
        plam = plon + 180.0D0
        if ( plam>180.0D0 ) plam = plam - 360.0D0
        zlampol = degrad*plam
        zphipol = degrad*pphi
        zsinpol = dsin(zphipol)
        zcospol = dcos(zphipol)
      end subroutine setup_rmc

      subroutine llij_rc(lat,lon,i,j)
        implicit none
        real(4) , intent(in) :: lat , lon
        real(4) , intent(out) :: i , j
        real(8) :: zarg , zarg1 , zarg2 , zlam , zphi
        real(8) :: lams , phis
 
        zphi = degrad*dble(lat)
        zlam = dble(lon)
        if ( zlam>180.0D0 ) zlam = zlam - 360.0D0
        zlam = degrad*zlam

        zarg = zcospol*dcos(zphi)*dcos(zlam-zlampol) + zsinpol*dsin(zphi)
        phis = dasin(zarg)
        phis = dlog(dtan(phis/2.0D0+datan(1.0D0)))*raddeg
        zarg1 = -dsin(zlam-zlampol)*dcos(zphi)
        zarg2 = -zsinpol*dcos(zphi)*dcos(zlam-zlampol) + zcospol*dsin(zphi)
        if ( dabs(zarg2)>=1.D-37 ) then
          lams = raddeg*datan2(zarg1,zarg2)
        else if ( dabs(zarg1)<1.D-37 ) then
          lams = 0.0D0
        else if ( zarg1>0.0D0 ) then
          lams = 90.0D0
        else
          lams = -90.0D0
        end if
        i = real(polei + (lams-xoff)/dlon)
        j = real(polej + (phis-yoff)/dlon)
      end subroutine llij_rc

      subroutine ijll_rc(i,j,lat,lon)
        implicit none
        real(4) , intent(in) :: i , j
        real(4) , intent(out) :: lat , lon
        real(8) :: xr , yr , arg , zarg1 , zarg2

        xr = xoff + (i-polei)*dlon
        if ( xr>180.0D0 ) xr = xr - 360.0D0
        xr = degrad*xr
        yr = yoff + dble((j-polej))*dlon
        yr = 2.0D0*datan(dexp(degrad*yr)) - datan(1.0D0)*2.0D0
 
        arg = zcospol*dcos(yr)*dcos(xr) + zsinpol*dsin(yr)
        lat = real(raddeg*dasin(arg))
        zarg1 = dsin(zlampol)*(-zsinpol*dcos(xr)*dcos(yr)+ &
              & zcospol*dsin(yr))-dcos(zlampol)*dsin(xr)*dcos(yr)
        zarg2 = dcos(zlampol)*(-zsinpol*dcos(xr)*dcos(yr)+ &
              & zcospol*dsin(yr))+dsin(zlampol)*dsin(xr)*dcos(yr)
        if ( dabs(zarg2)>=1.D-37 ) then
          lon = real(raddeg*datan2(zarg1,zarg2))
        else if ( dabs(zarg1)<1.D-37 ) then
          lon = 0.0
        else if ( zarg1>0.0D0 ) then
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
        if (lat > 89.999999) zrla = 0.0D0
        zrlap = zlampol - zrla
        zarg1 = zcospol*dsin(zrlap)
        zarg2 = zsinpol*dcos(zphi) - zcospol*dsin(zphi)*dcos(zrlap)
        znorm = 1.0D0/dsqrt(zarg1**2.0D0+zarg2**2.0D0)
        alpha = real(dacos(zarg2*znorm))
      end subroutine uvrot_rc

      subroutine mapfac_rc(ir, xmap)
        implicit none
        real(4) , intent(in) :: ir
        real(4) , intent(out) :: xmap
        real(8) :: yr
        yr = yoff + (ir-polej)*dlon
        xmap = real(1.0D0/dcos(yr*degrad))
      end subroutine mapfac_rc

      function rounder(rval,ltop)
        implicit none
        real(4) , intent(in) :: rval
        logical, intent(in) :: ltop
        real(4) :: rounder
        integer :: tmpval
        if (ltop) then
          tmpval = ceiling(rval*100.0)
        else
          tmpval = floor(rval*100.0)
        end if
        rounder = real(tmpval)/100.0
      end function rounder
!
      end module mod_projections

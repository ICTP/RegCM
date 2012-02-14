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

  use mod_realkinds
  use mod_constants

  real(dp) :: conefac
  real(dp) , private :: stdlon
  real(dp) , private :: truelat1 , truelat2 , tl1r , tl2r , ctl1r
  real(dp) , private :: colat1 , colat2 , nfac
  real(dp) , private :: rsw , rebydx , hemi
  real(dp) , private :: reflon , dlon , scale_top
  real(dp) , private :: polei , polej
  real(dp) , private :: polelon , polelat , xoff , yoff
  real(dp) , private :: zsinpol , zcospol , zlampol , zphipol
  logical , private :: lamtan

  contains

  subroutine setup_lcc(clat,clon,ci,cj,ds,slon,trlat1,trlat2)
    implicit none
    real(dp) , intent(in) :: ci , cj , slon , clat , clon , ds , &
                            trlat1 , trlat2
    real(dp) :: arg , deltalon1
!
    stdlon = slon
    truelat1 = trlat1
    truelat2 = trlat2

    tl1r = truelat1*degrad
    tl2r = truelat2*degrad
    colat1  = degrad*(deg90 - truelat1)
    colat2  = degrad*(deg90 - truelat2)

    nfac = (log(sin(colat1))        - log(sin(colat2))) / &
           (log(tan(colat1*d_half)) - log(tan(colat2*d_half)))

    if (truelat1 > d_zero) then
      hemi =  d_one
    else
      hemi = -d_one
    end if
    rebydx = earthrad / ds

    if ( abs(truelat1-truelat2) > 0.1D0) then
      conefac = log10(cos(tl1r)) - log10(cos(tl2r))
      conefac = conefac / &
                (log10(tan((deg45-abs(truelat1)/d_two)*degrad)) - &
                 log10(tan((deg45-abs(truelat2)/d_two)*degrad)))
      lamtan = .false.
    else
      conefac = sin(abs(tl1r))
      lamtan = .true.
    end if
    deltalon1 = clon - stdlon
    if (deltalon1 >  deg180) deltalon1 = deltalon1 - deg360
    if (deltalon1 < -deg180) deltalon1 = deltalon1 + deg360

    ctl1r = cos(tl1r)
    rsw = rebydx * ctl1r/conefac * &
             (tan((deg90*hemi-clat)    *degrad*d_half) / &
              tan((deg90*hemi-truelat1)*degrad*d_half))**conefac
    arg = conefac*(deltalon1*degrad)
    polei = hemi*ci - hemi * rsw * sin(arg)
    polej = hemi*cj + rsw * cos(arg)

  end subroutine setup_lcc

  subroutine ijll_lc(i,j,lat,lon)
    real(sp) , intent(in) :: i , j
    real(sp) , intent(out) :: lat , lon
    real(dp) :: chi1 , chi2 , chi
    real(dp) :: inew , jnew , xx , yy , r2 , r

    chi1 = (deg90 - hemi*truelat1)*degrad
    chi2 = (deg90 - hemi*truelat2)*degrad
    inew = hemi * dble(i)
    jnew = hemi * dble(j)
    xx = inew - polei
    yy = polej - jnew
    r2 = (xx*xx + yy*yy)
    r = sqrt(r2)/rebydx
    if (abs(r2) < dlowval) then
      lat = real(hemi * deg90)
      lon = real(stdlon)
    else
      lon = real(stdlon + raddeg * atan2(hemi*xx,yy)/conefac)
      lon = mod(lon+360.0, 360.0)
      if (abs(chi1-chi2) < dlowval) then
        chi = d_two*atan((r/tan(chi1))**(d_one/conefac)*tan(chi1*d_half))
      else
        chi = d_two*atan((r*conefac/sin(chi1))**(d_one/conefac) * &
              tan(chi1*d_half))
      end if
      lat = real((deg90-chi*raddeg)*hemi)
    end if
    if (lon >  180.0) lon = lon - 360.0
    if (lon < -180.0) lon = lon + 360.0
  end subroutine ijll_lc

  subroutine llij_lc(lat,lon,i,j)
    implicit none
    real(sp) , intent(in) :: lat , lon
    real(sp) , intent(out) :: i , j
    real(dp) :: arg , deltalon , rm

    deltalon = lon - stdlon
    if (deltalon > +deg180) deltalon = deltalon - deg360
    if (deltalon < -deg180) deltalon = deltalon + deg360

    rm = rebydx * ctl1r/conefac * &
           (tan((deg90*hemi-dble(lat))*degrad*d_half) / &
            tan((deg90*hemi-truelat1)*degrad*d_half))**conefac
    arg = conefac*(deltalon*degrad)
    i = real(hemi*(polei + hemi * rm * sin(arg)))
    j = real(hemi*(polej - rm * cos(arg)))
  end subroutine llij_lc

  subroutine uvrot_lc(lon, alpha)
    implicit none
    real(sp) , intent(in) :: lon
    real(sp) , intent(out) :: alpha
    real(dp) :: deltalon

    deltalon = stdlon - dble(lon)
    if (deltalon > +deg180) deltalon = deltalon - deg360
    if (deltalon < -deg180) deltalon = deltalon + deg360
    alpha = real(deltalon*degrad*conefac)
  end subroutine uvrot_lc

  subroutine mapfac_lc(lat, xmap)
    implicit none
    real(sp) , intent(in) :: lat
    real(sp) , intent(out) :: xmap
    real(dp) :: colat

    colat = degrad*(deg90-lat)
    if (.not. lamtan) then
      xmap = real(sin(colat2)/sin(colat) * &
             (tan(colat*d_half)/tan(colat2*d_half))**nfac)
    else
      xmap = real(sin(colat1)/sin(colat) * &
             (tan(colat*d_half)/tan(colat1*d_half))**cos(colat1))
    endif
  end subroutine mapfac_lc

  subroutine setup_plr(clat,clon,ci,cj,ds,slon)
    implicit none
    real(dp) , intent(in) :: clat , clon , cj , ci , ds , slon
    real(dp) :: ala1 , alo1

    stdlon = slon
    if (clat > 0.0) then
      hemi = d_one
    else
      hemi = -d_one
    end if
    rebydx = earthrad / ds
    reflon = stdlon + deg90
    ala1 = clat*degrad
    alo1 = (clon-reflon)*degrad
    scale_top = d_one + hemi * sin(ala1)
    rsw = rebydx*cos(ala1)*scale_top/(d_one+hemi*sin(ala1))
    polei = ci - rsw * cos(alo1)
    polej = cj - hemi * rsw * sin(alo1)
  end subroutine setup_plr

  subroutine llij_ps(lat,lon,i,j)
    implicit none
    real(sp) , intent(in) :: lat , lon
    real(sp) , intent(out) :: i , j
    real(dp) :: ala , alo , rm , deltalon

    deltalon = dble(lon) - reflon
    if (deltalon > +deg180) deltalon = deltalon - deg360
    if (deltalon < -deg180) deltalon = deltalon + deg360
    alo = deltalon * degrad
    ala = lat * degrad

    rm = rebydx * cos(ala) * scale_top/(d_one + hemi * sin(ala))
    i = real(polei + rm * cos(alo))
    j = real(polej + hemi * rm * sin(alo))

  end subroutine llij_ps

  subroutine ijll_ps(i,j,lat,lon)
    implicit none
    real(sp) , intent(in) :: i , j
    real(sp) , intent(out) :: lat , lon
    real(dp) :: xx , yy , r2 , gi2 , arcc

    xx = dble(i) - polei
    yy = (dble(j) - polej) * hemi
    r2 = xx**d_two + yy**d_two
    if (abs(r2) < dlowval) then
      lat = real(hemi*deg90)
      lon = real(reflon)
    else
      gi2 = (rebydx * scale_top)**d_two
      lat = real(raddeg * hemi * asin((gi2-r2)/(gi2+r2)))
      arcc = acos(xx/sqrt(r2))
      if (yy > d_zero) then
        lon = real(reflon + raddeg * arcc)
      else
        lon = real(reflon - raddeg * arcc)
      end if
    end if
    if (lon >  180.0) lon = lon - 360.0
    if (lon < -180.0) lon = lon + 360.0
  end subroutine ijll_ps
 
  subroutine mapfac_ps(lat, xmap)
    implicit none
    real(sp) , intent(in) :: lat
    real(sp) , intent(out) :: xmap
    xmap = real(scale_top/(d_one + hemi * sin(lat*degrad)))
  end subroutine mapfac_ps

  subroutine uvrot_ps(lon, alpha)
    implicit none
    real(sp) , intent(in) :: lon
    real(sp) , intent(out) :: alpha
    real(dp) :: deltalon

    deltalon = stdlon - dble(lon)
    if (deltalon > +deg180) deltalon = deltalon - deg360
    if (deltalon < -deg180) deltalon = deltalon + deg360
    alpha = real(deltalon*degrad*hemi)
  end subroutine uvrot_ps

  subroutine setup_mrc(clat,clon,ci,cj,ds)
    implicit none
    real(dp) , intent(in) :: clat , clon , cj , ci , ds
    real(dp) :: clain

    stdlon = clon
    clain = cos(clat*degrad)
    dlon = ds / (earthrad * clain)
    rsw = d_zero
    if (abs(clat) > dlowval) then
      rsw = (log(tan(d_half*((clat+deg90)*degrad))))/dlon
    end if
    polei = ci
    polej = cj
  end subroutine setup_mrc

  subroutine llij_mc(lat,lon,i,j)
    implicit none
    real(sp) , intent(in) :: lat , lon
    real(sp) , intent(out) :: i , j
    real(dp) :: deltalon

    deltalon = lon - stdlon
    if (deltalon > +deg180) deltalon = deltalon - deg360
    if (deltalon < -deg180) deltalon = deltalon + deg360
    i = real(polei + (deltalon/(dlon*raddeg)))
    j = real(polej + &
         (log(tan(d_half*((dble(lat)+deg90)*degrad)))) / dlon - rsw)
  end subroutine llij_mc

  subroutine ijll_mc(i,j,lat,lon)
    implicit none
    real(sp) , intent(in) :: i , j
    real(sp) , intent(out) :: lat , lon

    lat = real(d_two*atan(exp(dlon*(rsw + dble(j)-polej)))*raddeg-deg90)
    lon = real((dble(i)-polei)*dlon*raddeg + stdlon)
    if (lon >  180.0) lon = lon - 360.0
    if (lon < -180.0) lon = lon + 360.0
  end subroutine ijll_mc
 
  subroutine mapfac_mc(lat, xmap)
    implicit none
    real(sp) , intent(in) :: lat
    real(sp) , intent(out) :: xmap
    xmap = real(d_one/cos(dble(lat)*degrad))
  end subroutine mapfac_mc

  subroutine setup_rmc(clat,clon,ci,cj,ds,plon,plat)
    implicit none
    real(dp) , intent(in) :: clat , clon , cj , ci , ds , plon , plat
    real(dp) :: plam , pphi

    polelon = plon
    polelat = plat
    dlon = ds*raddeg/earthrad
    xoff = clon - plon
    yoff = clat - plat
    polei = ci
    polej = cj
    pphi = deg90 - plat
    plam = plon + deg180
    if ( plam>deg180 ) plam = plam - deg360
    zlampol = degrad*plam
    zphipol = degrad*pphi
    zsinpol = sin(zphipol)
    zcospol = cos(zphipol)
  end subroutine setup_rmc

  subroutine llij_rc(lat,lon,i,j)
    implicit none
    real(sp) , intent(in) :: lat , lon
    real(sp) , intent(out) :: i , j
    real(dp) :: zarg , zarg1 , zarg2 , zlam , zphi
    real(dp) :: lams , phis
 
    zphi = degrad*dble(lat)
    zlam = dble(lon)
    if ( zlam>deg180 ) zlam = zlam - deg360
    zlam = degrad*zlam

    zarg = zcospol*cos(zphi)*cos(zlam-zlampol) + zsinpol*sin(zphi)
    phis = asin(zarg)
    phis = log(tan(phis*d_half+atan(d_one)))*raddeg
    zarg1 = -sin(zlam-zlampol)*cos(zphi)
    zarg2 = -zsinpol*cos(zphi)*cos(zlam-zlampol) + zcospol*sin(zphi)
    if ( abs(zarg2)>=dlowval ) then
      lams = raddeg*atan2(zarg1,zarg2)
    else if ( abs(zarg1)<dlowval ) then
      lams = deg00
    else if ( zarg1>d_zero ) then
      lams = deg90
    else
      lams = -deg90
    end if
    i = real(polei + (lams-xoff)/dlon)
    j = real(polej + (phis-yoff)/dlon)
  end subroutine llij_rc

  subroutine ijll_rc(i,j,lat,lon)
    implicit none
    real(sp) , intent(in) :: i , j
    real(sp) , intent(out) :: lat , lon
    real(dp) :: xr , yr , arg , zarg1 , zarg2

    xr = xoff + (i-polei)*dlon
    if ( xr>deg180 ) xr = xr - deg360
    xr = degrad*xr
    yr = yoff + (j-polej)*dlon
    yr = d_two*atan(exp(degrad*yr)) - atan(d_one)*d_two
 
    arg = zcospol*cos(yr)*cos(xr) + zsinpol*sin(yr)
    lat = real(raddeg*asin(arg))
    zarg1 = sin(zlampol)*(-zsinpol*cos(xr)*cos(yr)+ &
            zcospol*sin(yr))-cos(zlampol)*sin(xr)*cos(yr)
    zarg2 = cos(zlampol)*(-zsinpol*cos(xr)*cos(yr)+ &
            zcospol*sin(yr))+sin(zlampol)*sin(xr)*cos(yr)
    if ( abs(zarg2)>=dlowval ) then
      lon = real(raddeg*atan2(zarg1,zarg2))
    else if ( abs(zarg1)<dlowval ) then
      lon = real(deg00)
    else if ( zarg1>d_zero ) then
      lon = real(deg90)
    else
      lon = real(-deg90)
    end if
    if (lon >  180.0) lon = lon - 360.0
    if (lon < -180.0) lon = lon + 360.0
  end subroutine ijll_rc
 
  subroutine uvrot_rc(lat, lon, alpha)
    implicit none
    real(sp) , intent(in) :: lon , lat
    real(sp) , intent(out) :: alpha
    real(dp) :: zphi , zrla , zrlap , zarg1 , zarg2 , znorm
    zphi = lat*degrad
    zrla = lon*degrad
    if (lat > deg90-0.00001D0) zrla = d_zero
    zrlap = zlampol - zrla
    zarg1 = zcospol*sin(zrlap)
    zarg2 = zsinpol*cos(zphi) - zcospol*sin(zphi)*cos(zrlap)
    znorm = d_one/dsqrt(zarg1**d_two+zarg2**d_two)
    alpha = real(dacos(zarg2*znorm))
  end subroutine uvrot_rc

  subroutine mapfac_rc(ir, xmap)
    implicit none
    real(sp) , intent(in) :: ir
    real(sp) , intent(out) :: xmap
    real(dp) :: yr
    yr = yoff + (ir-polej)*dlon
    xmap = real(d_one/dcos(yr*degrad))
  end subroutine mapfac_rc

  function rounder(xval,ltop)
    implicit none
    real(sp) , intent(in) :: xval
    logical, intent(in) :: ltop
    real(sp) :: rounder
    integer :: tmpval
    if (ltop) then
      tmpval = ceiling(xval*100.0)
    else
      tmpval = floor(xval*100.0)
    end if
    rounder = real(tmpval)/100.0
  end function rounder
!
end module mod_projections

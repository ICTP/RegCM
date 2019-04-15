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

  use mod_intkinds
  use mod_realkinds
  use mod_constants

  implicit none

  private

  real(rkx) :: stdlon
  real(rkx) :: truelat1 , truelat2 , tl1r , tl2r , ctl1r
  real(rkx) :: colat1 , colat2 , nfac
  real(rkx) :: rsw , rebydx , hemi
  real(rkx) :: reflon , dlon , scale_top
  real(rkx) :: polei , polej
  real(rkx) :: polelon , polelat , xoff , yoff
  real(rkx) :: zsinpol , zcospol , zlampol , zphipol
  logical :: lamtan

  real(rkx) , public :: conefac
  public :: setup_lcc , llij_lc , ijll_lc , uvrot_lc , mapfac_lc
  public :: setup_plr , llij_ps , ijll_ps , uvrot_ps , mapfac_ps
  public :: setup_mrc , llij_mc , ijll_mc , mapfac_mc
  public :: setup_rmc , llij_rc , ijll_rc , uvrot_rc , mapfac_rc
  public :: rounder

  contains

  subroutine setup_lcc(clat,clon,ci,cj,ds,slon,trlat1,trlat2)
    implicit none
    real(rkx) , intent(in) :: ci , cj , slon , clat , clon , ds , &
                            trlat1 , trlat2
    real(rkx) :: arg , deltalon1

    stdlon = slon
    truelat1 = trlat1
    truelat2 = trlat2
    tl1r = truelat1*degrad
    tl2r = truelat2*degrad
    colat1  = degrad*(deg90 - truelat1)
    colat2  = degrad*(deg90 - truelat2)
    nfac = (log(sin(colat1))        - log(sin(colat2))) / &
           (log(tan(colat1*d_half)) - log(tan(colat2*d_half)))
    if ( truelat1 > d_zero ) then
      hemi =  d_one
    else
      hemi = -d_one
    end if
    rebydx = earthrad / ds
    if ( abs(truelat1-truelat2) > 0.1_rkx ) then
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
    if ( deltalon1 >  deg180 ) deltalon1 = deltalon1 - deg360
    if ( deltalon1 < -deg180 ) deltalon1 = deltalon1 + deg360
    ctl1r = cos(tl1r)
    rsw = rebydx * ctl1r/conefac * &
             (tan((deg90*hemi-clat)    *degrad*d_half) / &
              tan((deg90*hemi-truelat1)*degrad*d_half))**conefac
    arg = conefac*(deltalon1*degrad)
    polei = hemi*ci - hemi * rsw * sin(arg)
    polej = hemi*cj + rsw * cos(arg)
  end subroutine setup_lcc

  subroutine ijll_lc(i,j,lat,lon)
    real(rkx) , intent(in) :: i , j
    real(rkx) , intent(out) :: lat , lon
    real(rkx) :: chi1 , chi2 , chi
    real(rkx) :: inew , jnew , xx , yy , r2 , r

    chi1 = (deg90 - hemi*truelat1)*degrad
    chi2 = (deg90 - hemi*truelat2)*degrad
    inew = hemi * i
    jnew = hemi * j
    xx = inew - polei
    yy = polej - jnew
    r2 = (xx*xx + yy*yy)
    r = sqrt(r2)/rebydx
    if ( abs(r2) < dlowval ) then
      lat = hemi * deg90
      lon = stdlon
    else
      lon = stdlon + raddeg * atan2(hemi*xx,yy)/conefac
      lon = mod(lon+deg360, deg360)
      if ( abs(chi1-chi2) < dlowval ) then
        chi = d_two*atan((r/tan(chi1))**(d_one/conefac)*tan(chi1*d_half))
      else
        chi = d_two*atan((r*conefac/sin(chi1))**(d_one/conefac) * &
              tan(chi1*d_half))
      end if
      lat = (deg90-chi*raddeg)*hemi
    end if
    if ( lon >  deg180 ) lon = lon - deg360
    if ( lon < -deg180 ) lon = lon + deg360
  end subroutine ijll_lc

  subroutine llij_lc(lat,lon,i,j)
    implicit none
    real(rkx) , intent(in) :: lat , lon
    real(rkx) , intent(out) :: i , j
    real(rkx) :: arg , deltalon , rm

    deltalon = lon - stdlon
    if ( deltalon > +deg180 ) deltalon = deltalon - deg360
    if ( deltalon < -deg180 ) deltalon = deltalon + deg360
    rm = rebydx * ctl1r/conefac * &
           (tan((deg90*hemi-lat)*degrad*d_half) / &
            tan((deg90*hemi-truelat1)*degrad*d_half))**conefac
    arg = conefac*(deltalon*degrad)
    i = hemi*(polei + hemi * rm * sin(arg))
    j = hemi*(polej - rm * cos(arg))
  end subroutine llij_lc

  subroutine uvrot_lc(lon, alpha)
    implicit none
    real(rkx) , intent(in) :: lon
    real(rkx) , intent(out) :: alpha
    real(rkx) :: deltalon

    deltalon = stdlon - lon
    if ( deltalon > +deg180 ) deltalon = deltalon - deg360
    if ( deltalon < -deg180 ) deltalon = deltalon + deg360
    alpha = deltalon*degrad*conefac
  end subroutine uvrot_lc

  subroutine mapfac_lc(lat, xmap)
    implicit none
    real(rkx) , intent(in) :: lat
    real(rkx) , intent(out) :: xmap
    real(rkx) :: colat

    colat = degrad*(deg90-lat)
    if ( .not. lamtan ) then
      xmap = sin(colat2)/sin(colat) * &
             (tan(colat*d_half)/tan(colat2*d_half))**nfac
    else
      xmap = sin(colat1)/sin(colat) * &
             (tan(colat*d_half)/tan(colat1*d_half))**cos(colat1)
    endif
  end subroutine mapfac_lc

  subroutine setup_plr(clat,clon,ci,cj,ds,slon)
    implicit none
    real(rkx) , intent(in) :: clat , clon , cj , ci , ds , slon
    real(rkx) :: ala1 , alo1

    stdlon = slon
    if ( clat > d_zero ) then
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
    real(rkx) , intent(in) :: lat , lon
    real(rkx) , intent(out) :: i , j
    real(rkx) :: ala , alo , rm , deltalon

    deltalon = lon - reflon
    if ( deltalon > +deg180 ) deltalon = deltalon - deg360
    if ( deltalon < -deg180 ) deltalon = deltalon + deg360
    alo = deltalon * degrad
    ala = lat * degrad
    rm = rebydx * cos(ala) * scale_top/(d_one + hemi * sin(ala))
    i = polei + rm * cos(alo)
    j = polej + hemi * rm * sin(alo)
  end subroutine llij_ps

  subroutine ijll_ps(i,j,lat,lon)
    implicit none
    real(rkx) , intent(in) :: i , j
    real(rkx) , intent(out) :: lat , lon
    real(rkx) :: xx , yy , r2 , gi2 , arcc

    xx = i - polei
    yy = (j - polej) * hemi
    r2 = xx**d_two + yy**d_two
    if ( abs(r2) < dlowval ) then
      lat = hemi*deg90
      lon = reflon
    else
      gi2 = (rebydx * scale_top)**d_two
      lat = raddeg * hemi * asin((gi2-r2)/(gi2+r2))
      arcc = acos(xx/sqrt(r2))
      if ( yy > d_zero ) then
        lon = reflon + raddeg * arcc
      else
        lon = reflon - raddeg * arcc
      end if
    end if
    if ( lon >  deg180 ) lon = lon - deg360
    if ( lon < -deg180 ) lon = lon + deg360
  end subroutine ijll_ps

  subroutine mapfac_ps(lat, xmap)
    implicit none
    real(rkx) , intent(in) :: lat
    real(rkx) , intent(out) :: xmap
    xmap = scale_top/(d_one + hemi * sin(lat*degrad))
  end subroutine mapfac_ps

  subroutine uvrot_ps(lon, alpha)
    implicit none
    real(rkx) , intent(in) :: lon
    real(rkx) , intent(out) :: alpha
    real(rkx) :: deltalon

    deltalon = stdlon - lon
    if ( deltalon > +deg180 ) deltalon = deltalon - deg360
    if ( deltalon < -deg180 ) deltalon = deltalon + deg360
    alpha = deltalon*degrad*hemi
  end subroutine uvrot_ps

  subroutine setup_mrc(clat,clon,ci,cj,ds)
    implicit none
    real(rkx) , intent(in) :: clat , clon , cj , ci , ds
    real(rkx) :: clain

    stdlon = clon
    clain = cos(clat*degrad)
    dlon = ds / (earthrad * clain)
    rsw = d_zero
    if ( abs(clat) > dlowval ) then
      rsw = (log(tan(d_half*((clat+deg90)*degrad))))/dlon
    end if
    polei = ci
    polej = cj
  end subroutine setup_mrc

  subroutine llij_mc(lat,lon,i,j)
    implicit none
    real(rkx) , intent(in) :: lat , lon
    real(rkx) , intent(out) :: i , j
    real(rkx) :: deltalon

    deltalon = lon - stdlon
    if ( deltalon > +deg180 ) deltalon = deltalon - deg360
    if ( deltalon < -deg180 ) deltalon = deltalon + deg360
    i = polei + (deltalon/(dlon*raddeg))
    j = polej + (log(tan(d_half*((lat+deg90)*degrad)))) / dlon - rsw
  end subroutine llij_mc

  subroutine ijll_mc(i,j,lat,lon)
    implicit none
    real(rkx) , intent(in) :: i , j
    real(rkx) , intent(out) :: lat , lon

    lat = d_two*atan(exp(dlon*(rsw + j-polej)))*raddeg-deg90
    lon = (i-polei)*dlon*raddeg + stdlon
    if ( lon >  deg180 ) lon = lon - deg360
    if ( lon < -deg180 ) lon = lon + deg360
  end subroutine ijll_mc

  subroutine mapfac_mc(lat, xmap)
    implicit none
    real(rkx) , intent(in) :: lat
    real(rkx) , intent(out) :: xmap
    xmap = d_one/cos(lat*degrad)
  end subroutine mapfac_mc

  subroutine setup_rmc(clat,clon,ci,cj,ds,plon,plat)
    implicit none
    real(rkx) , intent(in) :: clat , clon , cj , ci , ds , plon , plat
    real(rkx) :: plam , pphi

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
    real(rkx) , intent(in) :: lat , lon
    real(rkx) , intent(out) :: i , j
    real(rkx) :: zarg , zarg1 , zarg2 , zlam , zphi
    real(rkx) :: lams , phis

    zphi = degrad*lat
    zlam = lon
    if ( zlam>deg180 ) zlam = zlam - deg360
    zlam = degrad*zlam
    zarg = zcospol*cos(zphi)*cos(zlam-zlampol) + zsinpol*sin(zphi)
    phis = asin(zarg)
    phis = log(tan(phis*d_half+atan(d_one)))*raddeg
    zarg1 = -sin(zlam-zlampol)*cos(zphi)
    zarg2 = -zsinpol*cos(zphi)*cos(zlam-zlampol) + zcospol*sin(zphi)
    if ( abs(zarg2) >= dlowval ) then
      lams = raddeg*atan2(zarg1,zarg2)
    else if ( abs(zarg1) < dlowval ) then
      lams = deg00
    else if ( zarg1 > d_zero ) then
      lams = deg90
    else
      lams = -deg90
    end if
    i = polei + (lams-xoff)/dlon
    j = polej + (phis-yoff)/dlon
  end subroutine llij_rc

  subroutine ijll_rc(i,j,lat,lon)
    implicit none
    real(rkx) , intent(in) :: i , j
    real(rkx) , intent(out) :: lat , lon
    real(rkx) :: xr , yr , arg , zarg1 , zarg2

    xr = xoff + (i-polei)*dlon
    if ( xr > deg180 ) xr = xr - deg360
    xr = degrad*xr
    yr = yoff + (j-polej)*dlon
    yr = d_two*atan(exp(degrad*yr)) - atan(d_one)*d_two
    arg = zcospol*cos(yr)*cos(xr) + zsinpol*sin(yr)
    lat = raddeg*asin(arg)
    zarg1 = sin(zlampol)*(-zsinpol*cos(xr)*cos(yr)+ &
            zcospol*sin(yr))-cos(zlampol)*sin(xr)*cos(yr)
    zarg2 = cos(zlampol)*(-zsinpol*cos(xr)*cos(yr)+ &
            zcospol*sin(yr))+sin(zlampol)*sin(xr)*cos(yr)
    if ( abs(zarg2) >= dlowval ) then
      lon = raddeg*atan2(zarg1,zarg2)
    else if ( abs(zarg1) < dlowval ) then
      lon = deg00
    else if ( zarg1 > d_zero ) then
      lon = deg90
    else
      lon = -deg90
    end if
    if ( lon >  deg180 ) lon = lon - deg360
    if ( lon < -deg180 ) lon = lon + deg360
  end subroutine ijll_rc

  subroutine uvrot_rc(lat, lon, alpha)
    implicit none
    real(rkx) , intent(in) :: lon , lat
    real(rkx) , intent(out) :: alpha
    real(rkx) :: zphi , zrla , zrlap , zarg1 , zarg2 , znorm
    zphi = lat*degrad
    zrla = lon*degrad
    if ( lat > deg90-0.00001_rkx ) zrla = d_zero
    zrlap = zlampol - zrla
    zarg1 = zcospol*sin(zrlap)
    zarg2 = zsinpol*cos(zphi) - zcospol*sin(zphi)*cos(zrlap)
    znorm = d_one/sqrt(zarg1**d_two+zarg2**d_two)
    alpha = acos(zarg2*znorm)
  end subroutine uvrot_rc

  subroutine mapfac_rc(ir, xmap)
    implicit none
    real(rkx) , intent(in) :: ir
    real(rkx) , intent(out) :: xmap
    real(rkx) :: yr
    yr = yoff + (ir-polej)*dlon
    xmap = d_one/cos(yr*degrad)
  end subroutine mapfac_rc

  function rounder(xval,ltop)
    implicit none
    real(rkx) , intent(in) :: xval
    logical , intent(in) :: ltop
    real(rkx) :: rounder
    integer(ik4) :: tmpval
    if ( ltop ) then
      tmpval = ceiling(xval*100.0_rkx)
    else
      tmpval = floor(xval*100.0_rkx)
    end if
    rounder = tmpval/100.0_rkx
  end function rounder

end module mod_projections
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

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

module mod_param_clm

  use mod_intkinds
  use mod_realkinds
  use mod_message
  use mod_stdio

  implicit none
!
  integer, private :: k

  integer(ik4) , parameter :: nfld = 18 , npft = 17 , nsoi = 10 ,        &
                         ipft = 1 , ilai = 2 , isai = 3 , itop = 4 ,&
                         ibot = 5 , ilak = 6 , iwtl = 7 , iglc = 8 ,&
                         iurb = 9 , isnd = 10 , icly = 11 ,         &
                         icol = 12 , ioro = 13 , iiso = 15 ,        &
                         ifma = 14 , iapin = 16 , ibpin = 17 ,      &
                         imbo = 18
!
  real(rk4) , dimension(nfld) :: glat1 , glat2 , glon1 , glon2 , vmin
  integer(ik4) , dimension(nfld) :: nlat , nlev , nlon , ntim
!     ** glev_st = soil level depths
  real(rk4) , dimension(nsoi) :: glev_st
  character(len=256) , dimension(nfld) :: infil
  character(len=64) , dimension(nfld) :: vnam
  character(len=64) :: vnam_lm , vnam_st
!
  data (glev_st(k),k=1,nsoi)/0.0175 , 0.0451 , 0.0906 , 0.1656 ,    &
        0.2892 , 0.4930 , 0.8290 , 1.3829 , 2.2962 , 3.4332/
 
!     ** Landmask information in each file
!     ** vnam_lm = name of landmask variable
  data vnam_lm/'LANDMASK'/
 
!     ***** INFORMATION ON EACH CLM3 VARIABLE
 
!     ** Plant functional types
  data vnam(ipft) , vmin(ipft)/'PCT_PFT' , -99.0/
  data infil(ipft)/'/CLM/mksrf_pft.nc'/
  data nlon(ipft) , nlat(ipft) , nlev(ipft) , ntim(ipft)/720 , 360 ,&
       npft , 1/
  data glon1(ipft) , glon2(ipft) , glat1(ipft) , glat2(ipft)        &
       / - 179.75 , 179.75 , -89.75 , 89.75/
 
!     ** Vegetation parameters (LAI, SAI, etc)
  data vnam(ilai) , vmin(ilai)/'MONTHLY_LAI' , -99.0/
  data vnam(isai) , vmin(isai)/'MONTHLY_SAI' , -99.0/
  data vnam(itop) , vmin(itop)/'MONTHLY_HEIGHT_TOP' , -99.0/
  data vnam(ibot) , vmin(ibot)/'MONTHLY_HEIGHT_BOT' , -99.0/
  data infil(ilai)/'/CLM/mksrf_lai.nc'/
  data infil(isai)/'/CLM/mksrf_lai.nc'/
  data infil(itop)/'/CLM/mksrf_lai.nc'/
  data infil(ibot)/'/CLM/mksrf_lai.nc'/
  data nlon(ilai) , nlat(ilai) , nlev(ilai) , ntim(ilai)/720 , 360 ,&
       npft , 12/
  data nlon(isai) , nlat(isai) , nlev(isai) , ntim(isai)/720 , 360 ,&
       npft , 12/
  data nlon(itop) , nlat(itop) , nlev(itop) , ntim(itop)/720 , 360 ,&
       npft , 12/
  data nlon(ibot) , nlat(ibot) , nlev(ibot) , ntim(ibot)/720 , 360 ,&
       npft , 12/
  data glon1(ilai) , glon2(ilai) , glat1(ilai) , glat2(ilai)        &
       / - 179.75 , 179.75 , -89.75 , 89.75/
  data glon1(isai) , glon2(isai) , glat1(isai) , glat2(isai)        &
       / - 179.75 , 179.75 , -89.75 , 89.75/
  data glon1(itop) , glon2(itop) , glat1(itop) , glat2(itop)        &
       / - 179.75 , 179.75 , -89.75 , 89.75/
  data glon1(ibot) , glon2(ibot) , glat1(ibot) , glat2(ibot)        &
       / - 179.75 , 179.75 , -89.75 , 89.75/
 
!     ** Land water (lake and wetland)
!     data vnam(ilak), vmin(ilak)  / 'PCT_LAKE', -99.0 /
!     data vnam(iwtl), vmin(iwtl)  / 'PCT_WETLAND', -99.0 /
!     data infil(ilak) / '/CLM/mksrf_lanwat.nc' /
!     data infil(iwtl) / '/CLM/mksrf_lanwat.nc' /
!     data nlon(ilak),nlat(ilak),nlev(ilak),ntim(ilak)
!     &   /        360,       180,         1,        1 /
!     data nlon(iwtl),nlat(iwtl),nlev(iwtl),ntim(iwtl)
!     &   /        360,       180,         1,        1 /
!     data glon1(ilak),glon2(ilak),glat1(ilak),glat2(ilak)
!     &   /         0.5,      359.5,      -89.5,       89.5 /
!     data glon1(iwtl),glon2(iwtl),glat1(iwtl),glat2(iwtl)
!     &   /         0.5,      359.5,      -89.5,       89.5 /
  data vnam(ilak) , vmin(ilak)/'PCT_LAKE' , -99.0/
  data vnam(iwtl) , vmin(iwtl)/'PCT_WETLAND' , -99.0/
  data infil(ilak)/'/CLM/mksrf_lanwat.nc'/
  data infil(iwtl)/'/CLM/mksrf_lanwat.nc'/
  data nlon(ilak) , nlat(ilak) , nlev(ilak) , ntim(ilak)/7200 ,     &
       3600 , 1 , 1/
  data nlon(iwtl) , nlat(iwtl) , nlev(iwtl) , ntim(iwtl)/7200 ,     &
       3600 , 1 , 1/
  data glon1(ilak) , glon2(ilak) , glat1(ilak) , glat2(ilak)        &
       / - 179.975 , 179.975 , -89.975 , 89.975/
  data glon1(iwtl) , glon2(iwtl) , glat1(iwtl) , glat2(iwtl)        &
       / - 179.975 , 179.975 , -89.975 , 89.975/
 
!     ** Glacier
  data vnam(iglc) , vmin(iglc)/'PCT_GLACIER' , -99.0/
  data infil(iglc)/'/CLM/mksrf_glacier.nc'/
  data nlon(iglc) , nlat(iglc) , nlev(iglc) , ntim(iglc)/7200 ,     &
       3600 , 1 , 1/
  data glon1(iglc) , glon2(iglc) , glat1(iglc) , glat2(iglc)        &
       / - 179.975 , 179.975 , -89.975 , 89.975/
 
!     ** Urban
  data vnam(iurb) , vmin(iurb)/'PCT_URBAN' , -99.0/
  data infil(iurb)/'/CLM/mksrf_urban.nc'/
  data nlon(iurb) , nlat(iurb) , nlev(iurb) , ntim(iurb)/7200 ,     &
       3600 , 1 , 1/
  data glon1(iurb) , glon2(iurb) , glat1(iurb) , glat2(iurb)        &
       / - 179.975 , 179.975 , -89.975 , 89.975/
 
!     ** Soil texture
  data vnam_st/'MAPUNITS'/
  data vnam(isnd) , vmin(isnd)/'PCT_SAND' , -99.0/
  data vnam(icly) , vmin(icly)/'PCT_CLAY' , -99.0/
  data infil(isnd)/'/CLM/mksrf_soitex.10level.nc'/
  data infil(icly)/'/CLM/mksrf_soitex.10level.nc'/
  data nlon(isnd) , nlat(isnd) , nlev(isnd) , ntim(isnd)/4320 ,     &
       2160 , 10 , 6998/
  data nlon(icly) , nlat(icly) , nlev(icly) , ntim(icly)/4320 ,     &
       2160 , 10 , 6998/
  data glon1(isnd) , glon2(isnd) , glat1(isnd) , glat2(isnd)        &
       / - 179.9583 , 179.9583 , -89.95834 , 89.95834/
  data glon1(icly) , glon2(icly) , glat1(icly) , glat2(icly)        &
       / - 179.9583 , 179.9583 , -89.95834 , 89.95834/
 
!     ** Soil color
!     data vnam(icol), vmin(icol)  / 'SOIL_COLOR', 0./
!     data infil(icol) / '/CLM/mksrf_soicol_clm2.nc' /
!     data nlon(icol),nlat(icol),nlev(icol),ntim(icol)
!     &   /       128,       64,        1,       1./
!     data glon1(icol),glon2(icol),glat1(icol),glat2(icol)
!     &   /          0.,   357.1875,   -87.8638,    87.8638 /
  data vnam(icol) , vmin(icol)/'SOIL_COLOR' , 0.0/
  data infil(icol)/'/CLM/mksrf_soicol_clm2.nc'/
  data nlon(icol) , nlat(icol) , nlev(icol) , ntim(icol)/7200 ,     &
       3600 , 1 , 1/
  data glon1(icol) , glon2(icol) , glat1(icol) , glat2(icol)        &
       / - 179.975 , 179.975 , -89.975 , 89.975/
 
!     ** Orography
  data vnam(ioro) , vmin(ioro)/'LANDFRAC' , -99./
  data infil(ioro)/'/CLM/mksrf_navyoro_20min.nc'/
  data nlon(ioro) , nlat(ioro) , nlev(ioro) , ntim(ioro)/7200 ,     &
       3600 , 1 , 1/
  data glon1(ioro) , glon2(ioro) , glat1(ioro) , glat2(ioro)        &
       / - 179.975 , 179.975 , -89.975 , 89.975/
!     data vnam(ioro), vmin(ioro)  / 'lufrac', -99. /
!     data infil(ioro) / '/CLM/mksrf_navyoro_20min.nc' /
!     data nlon(ioro),nlat(ioro),nlev(ioro),ntim(ioro)
!     &   /      1080,      540,        20,       1 /
!     data glon1(ioro),glon2(ioro),glat1(ioro),glat2(ioro)
!     &   /   0.1666667,   359.8333,  -89.83334,   89.83334 /
 
!     ***** ADDITION ISOPRENE *****
  data vnam(iiso) , vmin(iiso)/'ISOP' , -99./
  data infil(iiso)/'mksrf_iso.nc'/
  data nlon(iiso) , nlat(iiso) , nlev(iiso) , ntim(iiso)/8571 ,     &
       3333 , 1 , 1/
  data glon1(iiso) , glon2(iiso) , glat1(iiso) , glat2(iiso)        &
       / - 179.97886665 , 179.961133350012 , -56.97857335 ,         &
       82.965426650004/
 
!     ***** ADDITION B-PINENE *****
  data vnam(ibpin) , vmin(ibpin)/'BPINE' , -99./
  data infil(ibpin)/'mksrf_pinb.nc'/
  data nlon(ibpin) , nlat(ibpin) , nlev(ibpin) , ntim(ibpin)/8571 , &
       3333 , 1 , 1/
  data glon1(ibpin) , glon2(ibpin) , glat1(ibpin) , glat2(ibpin)    &
       / - 179.97886665 , 179.961133350012 , -56.97857335 ,         &
       82.965426650004/
 
!     ***** ADDITION A-PINENE *****
  data vnam(iapin) , vmin(iapin)/'APINE' , -99./
  data infil(iapin)/'mksrf_pina.nc'/
  data nlon(iapin) , nlat(iapin) , nlev(iapin) , ntim(iapin)/8571 , &
       3333 , 1 , 1/
  data glon1(iapin) , glon2(iapin) , glat1(iapin) , glat2(iapin)    &
       / - 179.97886665 , 179.961133350012 , -56.97857335 ,         &
       82.965426650004/
 
!     ***** ADDITION METHYLBUTENOL *****
  data vnam(imbo) , vmin(imbo)/'MBO' , -99./
  data infil(imbo)/'mksrf_mbo.nc'/
  data nlon(imbo) , nlat(imbo) , nlev(imbo) , ntim(imbo)/8571 ,     &
       4286 , 1 , 1/
  data glon1(imbo) , glon2(imbo) , glat1(imbo) , glat2(imbo)        &
       / - 179.979 , 179.961000000012 , -89.979 , 89.9910000000055/
 
!     **** ADDITION Maximum Fractional Saturated Area ***
  data vnam(ifma) , vmin(ifma)/'FMAX' , -99.0/
  data infil(ifma)/'/CLM/mksrf_fmax.nc'/
  data nlon(ifma) , nlat(ifma) , nlev(ifma) , ntim(ifma)/720 , 360 ,&
       1 , 1/
  data glon1(ifma) , glon2(ifma) , glat1(ifma) , glat2(ifma)        &
       / - 179.75 , 179.75 , -89.75 , 89.75/

  contains

  subroutine param(nx,ny,kz,xlat,xlon,varmin,varmax,xlat1d,xlon1d,  &
                   xlonmin,xlonmax,xlatmin,xlatmax,iadim,ndim)
 
  implicit none
!
  real(rk4) :: xlatmax , xlatmin , xlonmax , xlonmin
  integer(ik4) :: kz , ndim , nx , ny
  integer(ik4) , dimension(ndim) :: iadim
  real(rk4) , dimension(ndim) :: varmax , varmin
  real(rk4) , dimension(ny,nx) :: xlat , xlon
  real(rk4) , dimension(ny) :: xlat1d
  real(rk4) , dimension(nx) :: xlon1d
  intent (in) kz , ndim , nx , ny , xlat , xlon
  intent (out) iadim , varmax , varmin , xlat1d , xlatmax , xlatmin ,&
               xlon1d , xlonmax , xlonmin
!
  integer(ik4) :: i , j
!
  varmin(1) = minval(xlon)
  varmin(2) = minval(xlat)
  varmin(3) = 1.0 !1050.
  varmax(1) = maxval(xlon)
  varmax(2) = maxval(xlat)
  varmax(3) = real(kz) !1050.
  iadim(1) = nx
  iadim(2) = ny
  iadim(3) = kz
  do i = 1 , nx
    xlon1d(i) = xlon(ny/2,i)
  end do
  do j = 1 , ny
    xlat1d(j) = xlat(j,nx/2)
  end do
  xlonmin = minval(xlon)
  xlonmax = maxval(xlon)
  xlatmin = minval(xlat)
  xlatmax = maxval(xlat)
 
  end subroutine param

  subroutine comp(fields,bvoc)
  implicit none
!
  integer(ik4) :: fields
  logical :: bvoc
  intent (in) bvoc
  intent (out) fields
!
  integer(ik4) :: numcompounds
 
  if ( bvoc ) then
 50     continue
    write(stdout,*) ' '
    write(stdout,*) ' '
    write(stdout,*) '********************************************'
    write(stdout,*) 'Creating biogenic emissions files'
    write(stdout,*) 'ENTER NUMBER OF SPECIES'
    read (*,*) numcompounds
    fields = 14 + numcompounds
    if ( fields>50 ) then
      call die('comp','Field maximum number is 50',1)
    else if ( fields<14 ) then
      go to 50
    else
    end if
  else
    fields = 14
  end if
 
  write(stdout,*) 'producing ' , fields , ' files'
  end subroutine comp

end module mod_param_clm

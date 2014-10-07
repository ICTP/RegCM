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

  integer(ik4) , parameter :: nfld = 19
  integer(ik4) , parameter :: npft = 17
  integer(ik4) , parameter :: nsoi = 10
  integer(ik4) , parameter :: iglc = 1 ,  ilai = 2 ,  isai = 3 ,  &
      itop = 4 ,  ibot = 5 ,  ilak = 6 ,  iwtl = 7 ,  ifma = 8 ,  &
      iurb = 9 ,  isnd = 10 , icly = 11 , icol = 12 , ioro = 13 , &
      ipft = 14 , iiso = 15 , iapin = 16 , ilimo = 17 , ibpin = 18 , imbo = 19
!
  real(rk4) , dimension(nfld) :: glat1 , glat2 , glon1 , glon2 , &
                                 dlat , dlon , vmin
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
  data nlev(ipft) , ntim(ipft)/npft , 1/

!     ** Vegetation parameters (LAI, SAI, etc)
  data vnam(ilai) , vmin(ilai)/'MONTHLY_LAI' , -99.0/
  data vnam(isai) , vmin(isai)/'MONTHLY_SAI' , -99.0/
  data vnam(itop) , vmin(itop)/'MONTHLY_HEIGHT_TOP' , -99.0/
  data vnam(ibot) , vmin(ibot)/'MONTHLY_HEIGHT_BOT' , -99.0/
  data infil(ilai)/'/CLM/mksrf_lai.nc'/
  data infil(isai)/'/CLM/mksrf_lai.nc'/
  data infil(itop)/'/CLM/mksrf_lai.nc'/
  data infil(ibot)/'/CLM/mksrf_lai.nc'/
  data nlev(ilai) , ntim(ilai)/ npft , 12/
  data nlev(isai) , ntim(isai)/ npft , 12/
  data nlev(itop) , ntim(itop)/ npft , 12/
  data nlev(ibot) , ntim(ibot)/ npft , 12/

!     ** Land water (lake and wetland)
!     data vnam(ilak), vmin(ilak)  / 'PCT_LAKE', -99.0 /
!     data vnam(iwtl), vmin(iwtl)  / 'PCT_WETLAND', -99.0 /
!     data infil(ilak) / '/CLM/mksrf_lanwat.nc' /
!     data infil(iwtl) / '/CLM/mksrf_lanwat.nc' /
!     data nlev(ilak),ntim(ilak) / 1,        1 /
!     data nlev(iwtl),ntim(iwtl) / 1,        1 /
  data vnam(ilak) , vmin(ilak)/'PCT_LAKE' , -99.0/
  data vnam(iwtl) , vmin(iwtl)/'PCT_WETLAND' , -99.0/
  data infil(ilak)/'/CLM/mksrf_lanwat.nc'/
  data infil(iwtl)/'/CLM/mksrf_lanwat.nc'/
  data nlev(ilak) , ntim(ilak)/ 1 , 1/
  data nlev(iwtl) , ntim(iwtl)/ 1 , 1/

!     ** Glacier
  data vnam(iglc) , vmin(iglc)/'PCT_GLACIER' , -99.0/
  data infil(iglc)/'/CLM/mksrf_glacier.nc'/
  data nlev(iglc) , ntim(iglc)/ 1 , 1/

!     ** Urban
  data vnam(iurb) , vmin(iurb)/'PCT_URBAN' , -99.0/
  data infil(iurb)/'/CLM/mksrf_urban.nc'/
  data nlev(iurb) , ntim(iurb)/ 1 , 1/

!     ** Soil texture
  data vnam_st/'MAPUNITS'/
  data vnam(isnd) , vmin(isnd)/'PCT_SAND' , -99.0/
  data vnam(icly) , vmin(icly)/'PCT_CLAY' , -99.0/
  data infil(isnd)/'/CLM/mksrf_soitex.10level.nc'/
  data infil(icly)/'/CLM/mksrf_soitex.10level.nc'/
  data nlev(isnd) , ntim(isnd)/ 10 , 6998/
  data nlev(icly) , ntim(icly)/ 10 , 6998/

!     ** Soil color
!     data vnam(icol), vmin(icol)  / 'SOIL_COLOR', 0./
!     data infil(icol) / '/CLM/mksrf_soicol_clm2.nc' /
!     data nlev(icol),ntim(icol) / 1,       1./
  data vnam(icol) , vmin(icol)/'SOIL_COLOR' , 0.0/
  data infil(icol)/'/CLM/mksrf_soicol_clm2.nc'/
  data nlev(icol) , ntim(icol)/ 1 , 1/

!     ** Orography
  data vnam(ioro) , vmin(ioro)/'LANDFRAC' , -99./
  data infil(ioro)/'/CLM/mksrf_navyoro_20min.nc'/
  data nlev(ioro) , ntim(ioro)/ 1 , 1/
!     data vnam(ioro), vmin(ioro)  / 'lufrac', -99. /
!     data infil(ioro) / '/CLM/mksrf_navyoro_20min.nc' /
!     data nlev(ioro),ntim(ioro) / 20,       1 /

!     ***** ADDITION ISOPRENE *****
  data vnam(iiso) , vmin(iiso)/'ISOP' , -99./
  data infil(iiso)/'CLM/mksrf_iso.nc'/
  data nlev(iiso) , ntim(iiso)/ 1 , 1/

!     ***** ADDITION B-PINENE *****
  data vnam(ibpin) , vmin(ibpin)/'BPINE' , -99./
  data infil(ibpin)/'CLM/mksrf_pinb.nc'/
  data nlev(ibpin) , ntim(ibpin)/ 1 , 1/

!     ***** ADDITION A-PINENE *****
  data vnam(iapin) , vmin(iapin)/'APIN' , -99./
  data infil(iapin)/'CLM/mksrf_pina.nc'/
  data nlev(iapin) , ntim(iapin)/ 1 , 1/

!     ***** ADDITION LIMONENE *****
  data vnam(ilimo) , vmin(ilimo)/'LIMO' , -99./
  data infil(ilimo)/'CLM/mksrf_limo.nc'/
  data nlev(ilimo) , ntim(ilimo)/ 1 , 1/

!     ***** ADDITION METHYLBUTENOL *****
  data vnam(imbo) , vmin(imbo)/'MBO' , -99./
  data infil(imbo)/'CLM/mksrf_mbo.nc'/
  data nlev(imbo) , ntim(imbo)/ 1 , 1/

!     **** ADDITION Maximum Fractional Saturated Area ***
  data vnam(ifma) , vmin(ifma)/'FMAX' , -99.0/
  data infil(ifma)/'/CLM/mksrf_fmax.nc'/
  data nlev(ifma) , ntim(ifma)/ 1 , 1/

  contains

  subroutine param(nx,ny,kz,xlat,xlon,varmin,varmax,xlat1d,xlon1d,  &
                   xlonmin,xlonmax,xlatmin,xlatmax,iadim,ndim)

  implicit none
!
  real(rk4) :: xlatmax , xlatmin , xlonmax , xlonmin
  integer(ik4) :: kz , ndim , nx , ny
  integer(ik4) , dimension(ndim) :: iadim
  real(rk4) , dimension(ndim) :: varmax , varmin
  real(rk4) , dimension(nx,ny) :: xlat , xlon
  real(rk4) , dimension(nx) :: xlon1d
  real(rk4) , dimension(ny) :: xlat1d
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
    xlon1d(i) = xlon(i,ny/2)
  end do
  do j = 1 , ny
    xlat1d(j) = xlat(nx/2,j)
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
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

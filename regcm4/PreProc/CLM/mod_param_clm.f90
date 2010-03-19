      module mod_param_clm
      implicit none
!
! PARAMETER definitions
!
!     ** Define number of CLM3 fields to interpolate to RegCM3 grid
!     ** nlon   = number of longitudes
!     ** nlat   = number of latitudes
!     ** nlev   = number of vertical levels or PFTs
!     ** ntim   = number of times
!     ** nfil   = number of fields
!     ** glon1  = first longitude in CLM3 file
!     ** glon2  = last longitude in CLM3 file
!     ** glat1  = first latitude in CLM3 file
!     ** glat2  = last latitude in CLM3 file
!     ** vnam   = variable name (short name)
!     ** vmin   = minimum possible data value (for undefined values)
!     ** outdir = name of directory where to put the interpolated files
!     ** outfil = name of interpolated CLM3/RegCM3 file
!     ** infil  = names of global CLM3 input files (NetCDF is defined
!     elsewhere)
!     ** npft = number of plant functional types
!     ** nsoi = number of soil levels

      integer , parameter :: nfld = 18 , npft = 17 , nsoi = 10 ,        &
                           & ipft = 1 , ilai = 2 , isai = 3 , itop = 4 ,&
                           & ibot = 5 , ilak = 6 , iwtl = 7 , iglc = 8 ,&
                           & iurb = 9 , isnd = 10 , icly = 11 ,         &
                           & icol = 12 , ioro = 13 , iiso = 15 ,        &
                           & ifma = 14 , iapin = 16 , ibpin = 17 ,      &
                           & imbo = 18
!
! Local variables
!
      real , dimension(nfld) :: glat1 , glat2 , glon1 , glon2 , vmin
!     ** glev_st = soil level depths
      real , dimension(nsoi) :: glev_st
      character(50) , dimension(nfld) :: infil
      integer :: k
      integer , dimension(nfld) :: nlat , nlev , nlon , ntim
      character(30) :: outdir , outfil
      character(20) , dimension(nfld) :: vnam
      character(20) :: vnam_lm , vnam_st
!
      data outdir/'../../Input/'/
      data outfil/'CLM3.INFO'/
 
      data (glev_st(k),k=1,nsoi)/0.0175 , 0.0451 , 0.0906 , 0.1656 ,    &
          & 0.2892 , 0.4930 , 0.8290 , 1.3829 , 2.2962 , 3.4332/
 
!     ** Landmask information in each file
!     ** vnam_lm = name of landmask variable
      data vnam_lm/'LANDMASK'/
 
!     ***** INFORMATION ON EACH CLM3 VARIABLE
 
!     ** Plant functional types
      data vnam(ipft) , vmin(ipft)/'PCT_PFT' , -99/
      data infil(ipft)/'mksrf_pft.nc'/
      data nlon(ipft) , nlat(ipft) , nlev(ipft) , ntim(ipft)/720 , 360 ,&
         & npft , 1/
      data glon1(ipft) , glon2(ipft) , glat1(ipft) , glat2(ipft)        &
         & / - 179.75 , 179.75 , -89.75 , 89.75/
 
!     ** Vegetation parameters (LAI, SAI, etc)
      data vnam(ilai) , vmin(ilai)/'MONTHLY_LAI' , -99/
      data vnam(isai) , vmin(isai)/'MONTHLY_SAI' , -99/
      data vnam(itop) , vmin(itop)/'MONTHLY_HEIGHT_TOP' , -99/
      data vnam(ibot) , vmin(ibot)/'MONTHLY_HEIGHT_BOT' , -99/
      data infil(ilai)/'mksrf_lai.nc'/
      data infil(isai)/'mksrf_lai.nc'/
      data infil(itop)/'mksrf_lai.nc'/
      data infil(ibot)/'mksrf_lai.nc'/
      data nlon(ilai) , nlat(ilai) , nlev(ilai) , ntim(ilai)/720 , 360 ,&
         & npft , 12/
      data nlon(isai) , nlat(isai) , nlev(isai) , ntim(isai)/720 , 360 ,&
         & npft , 12/
      data nlon(itop) , nlat(itop) , nlev(itop) , ntim(itop)/720 , 360 ,&
         & npft , 12/
      data nlon(ibot) , nlat(ibot) , nlev(ibot) , ntim(ibot)/720 , 360 ,&
         & npft , 12/
      data glon1(ilai) , glon2(ilai) , glat1(ilai) , glat2(ilai)        &
         & / - 179.75 , 179.75 , -89.75 , 89.75/
      data glon1(isai) , glon2(isai) , glat1(isai) , glat2(isai)        &
         & / - 179.75 , 179.75 , -89.75 , 89.75/
      data glon1(itop) , glon2(itop) , glat1(itop) , glat2(itop)        &
         & / - 179.75 , 179.75 , -89.75 , 89.75/
      data glon1(ibot) , glon2(ibot) , glat1(ibot) , glat2(ibot)        &
         & / - 179.75 , 179.75 , -89.75 , 89.75/
 
!     ** Land water (lake and wetland)
!     data vnam(ilak), vmin(ilak)  / 'PCT_LAKE', -99 /
!     data vnam(iwtl), vmin(iwtl)  / 'PCT_WETLAND', -99 /
!     data infil(ilak) / 'mksrf_lanwat.nc' /
!     data infil(iwtl) / 'mksrf_lanwat.nc' /
!     data nlon(ilak),nlat(ilak),nlev(ilak),ntim(ilak)
!     &   /        360,       180,         1,        1 /
!     data nlon(iwtl),nlat(iwtl),nlev(iwtl),ntim(iwtl)
!     &   /        360,       180,         1,        1 /
!     data glon1(ilak),glon2(ilak),glat1(ilak),glat2(ilak)
!     &   /         0.5,      359.5,      -89.5,       89.5 /
!     data glon1(iwtl),glon2(iwtl),glat1(iwtl),glat2(iwtl)
!     &   /         0.5,      359.5,      -89.5,       89.5 /
      data vnam(ilak) , vmin(ilak)/'PCT_LAKE' , -99/
      data vnam(iwtl) , vmin(iwtl)/'PCT_WETLAND' , -99/
      data infil(ilak)/'mksrf_lanwat.nc'/
      data infil(iwtl)/'mksrf_lanwat.nc'/
      data nlon(ilak) , nlat(ilak) , nlev(ilak) , ntim(ilak)/7200 ,     &
         & 3600 , 1 , 1/
      data nlon(iwtl) , nlat(iwtl) , nlev(iwtl) , ntim(iwtl)/7200 ,     &
         & 3600 , 1 , 1/
      data glon1(ilak) , glon2(ilak) , glat1(ilak) , glat2(ilak)        &
         & / - 179.975 , 179.975 , -89.975 , 89.975/
      data glon1(iwtl) , glon2(iwtl) , glat1(iwtl) , glat2(iwtl)        &
         & / - 179.975 , 179.975 , -89.975 , 89.975/
 
!     ** Glacier
      data vnam(iglc) , vmin(iglc)/'PCT_GLACIER' , -99/
      data infil(iglc)/'mksrf_glacier.nc'/
      data nlon(iglc) , nlat(iglc) , nlev(iglc) , ntim(iglc)/7200 ,     &
         & 3600 , 1 , 1/
      data glon1(iglc) , glon2(iglc) , glat1(iglc) , glat2(iglc)        &
         & / - 179.975 , 179.975 , -89.975 , 89.975/
 
!     ** Urban
      data vnam(iurb) , vmin(iurb)/'PCT_URBAN' , -99/
      data infil(iurb)/'mksrf_urban.nc'/
      data nlon(iurb) , nlat(iurb) , nlev(iurb) , ntim(iurb)/7200 ,     &
         & 3600 , 1 , 1/
      data glon1(iurb) , glon2(iurb) , glat1(iurb) , glat2(iurb)        &
         & / - 179.975 , 179.975 , -89.975 , 89.975/
 
!     ** Soil texture
      data vnam_st/'MAPUNITS'/
      data vnam(isnd) , vmin(isnd)/'PCT_SAND' , -99/
      data vnam(icly) , vmin(icly)/'PCT_CLAY' , -99/
      data infil(isnd)/'mksrf_soitex.10level.nc'/
      data infil(icly)/'mksrf_soitex.10level.nc'/
      data nlon(isnd) , nlat(isnd) , nlev(isnd) , ntim(isnd)/4320 ,     &
         & 2160 , 10 , 6998/
      data nlon(icly) , nlat(icly) , nlev(icly) , ntim(icly)/4320 ,     &
         & 2160 , 10 , 6998/
      data glon1(isnd) , glon2(isnd) , glat1(isnd) , glat2(isnd)        &
         & / - 179.9583 , 179.9583 , -89.95834 , 89.95834/
      data glon1(icly) , glon2(icly) , glat1(icly) , glat2(icly)        &
         & / - 179.9583 , 179.9583 , -89.95834 , 89.95834/
 
!     ** Soil color
!     data vnam(icol), vmin(icol)  / 'SOIL_COLOR', 0 /
!     data infil(icol) / 'mksrf_soicol_clm2.nc' /
!     data nlon(icol),nlat(icol),nlev(icol),ntim(icol)
!     &   /       128,       64,        1,       1 /
!     data glon1(icol),glon2(icol),glat1(icol),glat2(icol)
!     &   /          0.,   357.1875,   -87.8638,    87.8638 /
      data vnam(icol) , vmin(icol)/'SOIL_COLOR' , 0/
      data infil(icol)/'mksrf_soicol_clm2.nc'/
      data nlon(icol) , nlat(icol) , nlev(icol) , ntim(icol)/7200 ,     &
         & 3600 , 1 , 1/
      data glon1(icol) , glon2(icol) , glat1(icol) , glat2(icol)        &
         & / - 179.975 , 179.975 , -89.975 , 89.975/
 
!     ** Orography
      data vnam(ioro) , vmin(ioro)/'LANDFRAC' , -99./
      data infil(ioro)/'mksrf_navyoro_20min.nc'/
      data nlon(ioro) , nlat(ioro) , nlev(ioro) , ntim(ioro)/7200 ,     &
         & 3600 , 1 , 1/
      data glon1(ioro) , glon2(ioro) , glat1(ioro) , glat2(ioro)        &
         & / - 179.975 , 179.975 , -89.975 , 89.975/
!     data vnam(ioro), vmin(ioro)  / 'lufrac', -99. /
!     data infil(ioro) / 'mksrf_navyoro_20min.nc' /
!     data nlon(ioro),nlat(ioro),nlev(ioro),ntim(ioro)
!     &   /      1080,      540,        20,       1 /
!     data glon1(ioro),glon2(ioro),glat1(ioro),glat2(ioro)
!     &   /   0.1666667,   359.8333,  -89.83334,   89.83334 /
 
!     ***** ADDITION ISOPRENE *****
      data vnam(iiso) , vmin(iiso)/'ISOP' , -99./
      data infil(iiso)/'mksrf_iso.nc'/
      data nlon(iiso) , nlat(iiso) , nlev(iiso) , ntim(iiso)/8571 ,     &
         & 3333 , 1 , 1/
      data glon1(iiso) , glon2(iiso) , glat1(iiso) , glat2(iiso)        &
         & / - 179.97886665 , 179.961133350012 , -56.97857335 ,         &
         & 82.965426650004/
 
!     ***** ADDITION B-PINENE *****
      data vnam(ibpin) , vmin(ibpin)/'BPINE' , -99./
      data infil(ibpin)/'mksrf_pinb.nc'/
      data nlon(ibpin) , nlat(ibpin) , nlev(ibpin) , ntim(ibpin)/8571 , &
         & 3333 , 1 , 1/
      data glon1(ibpin) , glon2(ibpin) , glat1(ibpin) , glat2(ibpin)    &
         & / - 179.97886665 , 179.961133350012 , -56.97857335 ,         &
         & 82.965426650004/
 
!     ***** ADDITION A-PINENE *****
      data vnam(iapin) , vmin(iapin)/'APINE' , -99./
      data infil(iapin)/'mksrf_pina.nc'/
      data nlon(iapin) , nlat(iapin) , nlev(iapin) , ntim(iapin)/8571 , &
         & 3333 , 1 , 1/
      data glon1(iapin) , glon2(iapin) , glat1(iapin) , glat2(iapin)    &
         & / - 179.97886665 , 179.961133350012 , -56.97857335 ,         &
         & 82.965426650004/
 
!     ***** ADDITION METHYLBUTENOL *****
      data vnam(imbo) , vmin(imbo)/'MBO' , -99./
      data infil(imbo)/'mksrf_mbo.nc'/
      data nlon(imbo) , nlat(imbo) , nlev(imbo) , ntim(imbo)/8571 ,     &
         & 4286 , 1 , 1/
      data glon1(imbo) , glon2(imbo) , glat1(imbo) , glat2(imbo)        &
         & / - 179.979 , 179.961000000012 , -89.979 , 89.9910000000055/
 
!     **** ADDITION Maximum Fractional Saturated Area ***
      data vnam(ifma) , vmin(ifma)/'FMAX' , -99/
      data infil(ifma)/'mksrf_fmax.nc'/
      data nlon(ifma) , nlat(ifma) , nlev(ifma) , ntim(ifma)/720 , 360 ,&
         & 1 , 1/
      data glon1(ifma) , glon2(ifma) , glat1(ifma) , glat2(ifma)        &
         & / - 179.75 , 179.75 , -89.75 , 89.75/

      end module mod_param_clm

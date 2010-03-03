      module mod_ccsm

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx!
!
! This is a package of subroutines to read CCSM T85 and T42 L26 data in
! NETCDF format and to prepare Initial and boundary conditions for RegCM3.
! Both Global and Window of CCSM data are acceptable.
! Written By Moetasim Ashfaq Dec-2005 @ PURDUE.EDU
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx!
!
!  SUBROUTINE CAM85
!  Read unpacked CCSM NETCDF T85 L26 (six hourly) data and save into data
!    arrays. Each varaible is read in seperate monthly data files.
!
!  SUBROUTINE CAM42
!  Read unpacked CCSM NETCDF T42 L26 (six hourly) data and save into data
!     arrays. Each varaible is read in seperate yearly data files.
!
!  SUBROUTINE HEADER_CAM85 & HEADER_CAM42
!  Define pressure levels, Ak and Bk coeffcients and global grid dimensions
!     In CCSM, the vertical coordinate is a hybrid sigma-pressure system
!     The pressure is defined as:
!     P=Ak*Po+Bk*PS(i,j)
!     All 3D fields required for ICBC are at mid-points, so
!     Ak refers to hyam, the hybrid A coefficient at midpoints, and
!     Bk refers to hybm, the hybrid B coefficient at midpoints
!     Po=1000mb
!
!  SUBROUTINE GET_CAM85
!  Main subroutine to read data arrays from CAM85 and prepare ICBCs at RCM Grid
!
!  SUBROUTINE GET_CAM42
!  Main subroutine to read data arrays from CAM42 and prepare ICBCs at RCM Grid
!
!  SUBROUTINE INITDATE3
!  Initialize CCSM 365 days calendar (No leap years)
!
!  SUBROUTINE CAMCLNDR
!  Subroutine for SST preparation with No leap year
!
!  SUBROUTINE HANDLE_ERR
!  Handle error for NETCDF calls
!
!  SUBROUTINES CALLED FROM ICBC.f
!     1) INTLIN
!     2) INTLOG
!     3) HUMIF1FV
!     3) BILINX2CR
!     4) BILINX2DT
!     5) TOP2BTM
!     6) UVROT4
!     7) INTGTB
!     8) P1P2
!     8) INTPSN
!     9) INTV3
!     10) MKSST
!     10) HUMID2FV
!     10) INTV1
!     11) WRITEF
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx!
!  DATA PREPARATION
!  Dataset required to use this code can be preapred using NCO utilitiles
!     such as NCKS, NCRCAT etc.
!  Prepare:
!     Monthly data files for CAM85
!     Yearly data files for CAM42
!  For example:
!     To extract global data of CAM42 for Specific Humidity
!     ncks -v Q input.nc ccsm.shum.nyear.nc    ,and
!     to extract a subset (window) of CAM85 data for Specific Humidity
!     ncks -d lat,min,max -d lon,min,max -v Q input.nc ccsm.shumJAN.nyear.nc
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx!
!  NAMING CONVENTION (Global Data Files) CAM85
!  (MONTH) = JAN/FEB/MAR/APR/MAY/JUN/JUL/AUG/SEP/OCT/NOV/DEC
!
!  ccsm.air(MONTH).nyear.nc     for 'T'     (Temperature)
!  ccsm.hgt(MONTH).nyear.nc     for 'Z3'    (Geopotential Height)
!  ccsm.shum(MONTH).nyear.nc    for 'Q'     (Specific Humidity)
!  ccsm.uwnd(MONTH).nyear.nc    for 'U'     (Zonal Wind)
!  ccsm.vwnd(MONTH).nyear.nc    for 'V'     (Meridonial Wind)
!  ccsm.pres(MONTH).nyear.nc    for 'PS'    (Surface Pressure)
!
!  PATH /DATA/CAM85/NYEAR/
!  ccsm_ht.nc      for 'PHIS'  (Surface Geopotential-static field)
!
!  PATH /DATA/CAM85/
!
!  NAMING CONVENTION (Global Data Files) CAM42
!  ccsm.air.nyear        for 'T'   (Temperature)
!  ccsm.hgt.nyear.nc     for 'Z3'  (Geopotential Height)
!  ccsm.shum.nyear.nc    for 'Q'   (Specific Humidity)
!  ccsm.uwnd.nyear.nc    for 'U'   (Zonal Wind)
!  ccsm.vwnd.nyear.nc    for 'V'   (Meridonial Wind)
!  ccsm.pres.nyear.nc    for 'PS'  (Surface Pressure)
!
!  PATH /DATA/CAM42/NYEAR/
!  ccsm_ht.nc      for 'PHIS'  (Surface Geopotential-static field)
!
!  PATH /DATA/CAM42/
!
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx!

      use mod_domain

      implicit none

      ! CAM58
      integer , parameter :: ilon = 256 , jlat = 128
      ! CAM42
      integer , parameter :: ilonh = 128 , jlath = 64

      integer , parameter :: klev = 26 , npl = 18

      integer , dimension(10) :: icount , istart
      integer :: inet1 , ivar1
      integer , dimension(6) :: inet6 , ivar6

      equivalence (b2(1,1,1),tp(1,1,1))
      equivalence (d2(1,1,1),up(1,1,1))
      real , dimension(jx,iy,npl*3) :: b3
      real , dimension(jx,iy) :: b3pd
      real , dimension(jx,iy,npl*2) :: d3
      real , dimension(jx,iy,npl) :: h3 , q3 , t3 , u3 , v3
      equivalence (b3(1,1,1),t3(1,1,1))
      equivalence (d3(1,1,1),u3(1,1,1))
      real , dimension(jx,iy,kz) :: h4 , q4 , t4 , u4 , v4
      real , dimension(jx,iy) :: ps4 , ts4
      real , dimension(ilon,jlat,npl*3) :: b2
      real , dimension(ilon,jlat,npl*2) :: d2
      real , dimension(ilon,jlat,npl) :: hp , qp , tp , up , vp
      real , dimension(ilon,jlat,klev) :: pp3d
      real , dimension(ilon,jlat,klev) :: hvar , qvar , tvar , uvar ,   &
                                        & vvar
      real , dimension(ilon,jlat) :: psvar , zsvar
      real , dimension(klev+1) :: ak , bk
      real , dimension(jlat) :: glat
      real , dimension(ilon) :: glon
      real , dimension(npl) :: pplev , sigma1 , sigmar
      real :: psref

      real , dimension(jx,iy) :: coriol , dlat , dlon , msfx , snowcv , &
                               & topogm , toposdgm , xlandu , xlat ,    &
                               & xlon
      real , dimension(kz) :: dsigma , sigma2
      real , dimension(kz+1) :: sigmaf
      real :: delx , grdfac
      character(6) :: lgtype

      contains

      subroutine get_cam42(idate)

        use netcdf

        implicit none
!
! Dummy arguments
!
        integer :: idate
!
! Local variables
!
        real , dimension(npl) :: c1
        real , dimension(kz) :: c2
        integer :: checklat , checklon , i , i0 , i1 , ii ,             &
                &  imax , imin , j , j0 , j1 , jj , k , latid ,         &
                &  latlen , lonid , lonlen , nmop , nyrp , istatus
        integer :: iomega
        real , dimension(jx,iy,npl) :: dum1
        real , dimension(jx,iy,npl,2) :: dum2
        real :: pp , w , wt
        real , dimension(jx,iy) :: pa , sst1 , sst2 , tlayer , za
        character(38) :: pathaddname
        logical :: there
        character(4) :: varname
        real , allocatable , dimension(:,:) :: work

        data varname/'PHIS'/
 
        call cam42(idate,idate1,xlon,xlat,glat,jx,iy,i0,i1,j0,j1)

        if ( idate==idate1 ) then
          pathaddname = '../DATA/CAM42/ccsm_ht.nc'
          inquire (file=pathaddname,exist=there)
          if ( .not.there ) then
            print * , pathaddname , 'is not available'
            stop
          end if
          istatus = nf90_open(pathaddname,nf90_nowrite,inet1)
          if ( istatus/=nf90_noerr ) call handle_err(istatus)
          istatus = nf90_inq_varid(inet1,varname,ivar1)
          if ( istatus/=nf90_noerr ) call handle_err(istatus)
          write (*,*) inet1 , pathaddname , ivar1 , varname
          istatus = nf90_inq_dimid(inet1,'lon',lonid)
          if ( istatus/=nf90_noerr ) call handle_err(istatus)
          istatus = nf90_inq_dimid(inet1,'lat',latid)
          if ( istatus/=nf90_noerr ) call handle_err(istatus)
          istatus = nf90_inquire_dimension(inet1,lonid,len=lonlen)
          if ( istatus/=nf90_noerr ) call handle_err(istatus)
          istatus = nf90_inquire_dimension(inet1,latid,len=latlen)
          if ( istatus/=nf90_noerr ) call handle_err(istatus)
          imin = i0
          imax = i1
          if ( i0>i1 ) imax = i1 + ilonh
          checklon = (imax-imin) + 1
          checklat = (j1-j0) + 1
          print * , 'i0,i1' , i0 , i1
          if ( checklon/=lonlen .or. checklat/=latlen ) then
            print * , 'DOMAIN DIMENSIONS DO NOT MATCH'
            print * , 'LAT for 3D Variables = ' , checklat
            print * , 'LAT for' , varname , '= ' , latlen
            print * , 'LON for 3D Variables = ' , checklon
            print * , 'LON for' , varname , '= ' , lonlen
            stop 'Check Dimensions in CCSM Data Files'
          end if
          allocate (work(lonlen,latlen))
          icount(1) = lonlen
          icount(2) = latlen
          icount(3) = 1
          istart(1) = 1
          istart(2) = 1
          istart(3) = 1
          istatus = nf90_get_var(inet1,ivar1,work,istart,icount)
          if ( istatus/=nf90_noerr ) call handle_err(istatus)
          istatus = nf90_close(inet1)
          do jj = j0 , j1
            j = jj - j0 + 1
            if ( i0>i1 ) then
              do ii = i0 , ilonh
                i = ii - i0 + 1
                zsvar(ii,jj) = work(i,j)/9.80665
              end do
              do ii = 1 , i1
                i = ii + (ilonh-i0)
                zsvar(ii,jj) = work(i,j)/9.80665
              end do
            else
              do ii = i0 , i1
                i = ii - i0 + 1
                zsvar(ii,jj) = work(i,j)/9.80665
              end do
            end if
          end do
        end if
        write (*,*) 'Read in fields at Date: ' , idate
 
        do k = 1 , klev
          do j = 1 , jlath
            do i = 1 , ilonh
              if ( psvar(i,j)>-9995. ) then
                pp3d(i,j,k) = psvar(i,j)*0.5*(bk(k)+bk(k+1))            &
                          & + 0.5*(ak(k)+ak(k+1))
              else
                pp3d(i,j,k) = -9999.0
              end if
            end do
          end do
        end do
 
        call height(hp,hvar,tvar,psvar,pp3d,zsvar,                      &
                  & ilonh,jlath,klev,pplev,npl)
 
        call intlin(up,uvar,psvar,pp3d,ilonh,jlath,klev,pplev,npl)
        call intlin(vp,vvar,psvar,pp3d,ilonh,jlath,klev,pplev,npl)
 
        call intlog(tp,tvar,psvar,pp3d,ilonh,jlath,klev,pplev,npl)
 
        call humid1fv(tvar,qvar,pp3d,ilonh,jlath,klev)
        call intlin(qp,qvar,psvar,pp3d,ilonh,jlath,klev,pplev,npl)
 
        call bilinxcx(b3,b2,xlon,xlat,glon,glat,jx,iy,ilonh,jlath,npl)
        call bilinxdt(d3,d2,dlon,dlat,glon,glat,jx,iy,ilonh,jlath,npl)
        call uvrot4(u3,v3,dlon,dlat,clon,clat,grdfac,jx,iy,             &
                &   npl,plon,plat,lgtype)
 
        call top2btm(t3,jx,iy,npl)
        call top2btm(q3,jx,iy,npl)
        call top2btm(h3,jx,iy,npl)
        call top2btm(u3,jx,iy,npl)
        call top2btm(v3,jx,iy,npl)
 
        call intgtb(pa,za,tlayer,topogm,t3,h3,sigmar,                   &
                  & jx,iy,npl,dum1,dum2)
 
        call intpsn(ps4,topogm,pa,za,tlayer,ptop,jx,iy)
        call p1p2(b3pd,ps4,jx,iy)
 
        call intv3(ts4,t3,ps4,sigmar,ptop,jx,iy,npl,dum1)
        call camclndr(idate,nyrp,nmop,wt)
        call mksst(ts4,sst1,sst2,topogm,xlandu,jx,iy,nyrp,nmop,wt)
        call intv1(u4,u3,b3pd,sigma2,sigmar,ptop,                       &
                &  jx,iy,kz,npl,1,1,dum2,c1,c2)
        call intv1(v4,v3,b3pd,sigma2,sigmar,ptop,                       &
               &   jx,iy,kz,npl,1,1,dum2,c1,c2)
        call intv2(t4,t3,ps4,sigma2,sigmar,ptop,                        &
               &   jx,iy,kz,npl,1,dum2,c1,c2)
 
        call intv1(q4,q3,ps4,sigma2,sigmar,ptop,                        &
               &   jx,iy,kz,npl,1,0,dum2,c1,c2)
        call humid2fv(t4,q4,ps4,ptop,sigma2,jx,iy,kz)
 
        call hydrost(h4,t4,topogm,ps4,ptop,sigmaf,sigma2,               &
               &     dsigma,jx,iy,kz)
 
        call writef(u4,v4,t4,q4,pp,w,ps4,ts4,ptop,jx,iy,kz,idate,iomega)

        deallocate(work)
      end subroutine get_cam42
!
!-----------------------------------------------------------------------
! 
      subroutine head_cam42
        implicit none
!
! Local variables
!
        integer :: i , k , kr
        real, dimension(jlath) :: fixlat
!
        data fixlat/ - 87.8638 , -85.0965 , -82.3129 , -79.5256 ,       &
         & -76.7369 ,                                                   &
         & -73.9475 , -71.1578 , -68.3678 , -65.5776 , -62.7874 ,       &
         & -59.9970 , -57.2066 , -54.4162 , -51.6257 , -48.8352 ,       &
         & -46.0447 , -43.2542 , -40.4636 , -37.6731 , -34.8825 ,       &
         & -32.0919 , -29.3014 , -26.5108 , -23.7202 , -20.9296 ,       &
         & -18.1390 , -15.3484 , -12.5578 , -9.76715 , -6.97653 ,       &
         & -4.18592 , -1.39531 , 1.39531 , 4.18592 , 6.97653 , 9.76715 ,&
         & 12.5578 , 15.3484 , 18.1390 , 20.9296 , 23.7202 , 26.5108 ,  &
         & 29.3014 , 32.0919 , 34.8825 , 37.6731 , 40.4636 , 43.2542 ,  &
         & 46.0447 , 48.8352 , 51.6257 , 54.4162 , 57.2066 , 59.9970 ,  &
         & 62.7874 , 65.5776 , 68.3678 , 71.1578 , 73.9475 , 76.7369 ,  &
         & 79.5256 , 82.3129 , 85.0965 , 87.8638/

        do i = 1 , ilonh
          glon(i) = float(i-1)*2.8125
        end do
        do i = 1 , jlath
          glat(i) = fixlat(i)
        end do
 
        pplev(1) = 30.
        pplev(2) = 50.
        pplev(3) = 70.
        pplev(4) = 100.
        pplev(5) = 150.
        pplev(6) = 200.
        pplev(7) = 250.
        pplev(8) = 300.
        pplev(9) = 350.
        pplev(10) = 420.
        pplev(11) = 500.
        pplev(12) = 600.
        pplev(13) = 700.
        pplev(14) = 780.
        pplev(15) = 850.
        pplev(16) = 920.
        pplev(17) = 960.
        pplev(18) = 1000.
 
!
        do k = 1 , npl
          sigmar(k) = pplev(k)*0.001
        end do
 
        do k = 1 , npl
          kr = npl - k + 1
          sigma1(k) = sigmar(kr)
        end do
        psref = 1000.
 
!       Ak is acatually Ak*Po
 
        ak(1) = 2.194067
        ak(2) = 4.895209
        ak(3) = 9.882418
        ak(4) = 18.05201
        ak(5) = 29.83724
        ak(6) = 44.62334
        ak(7) = 61.60587
        ak(8) = 78.51243
        ak(9) = 77.31271
        ak(10) = 75.90131
        ak(11) = 74.24086
        ak(12) = 72.28744
        ak(13) = 69.98932
        ak(14) = 67.28574
        ak(15) = 64.10509
        ak(16) = 60.36322
        ak(17) = 55.96111
        ak(18) = 50.78225
        ak(19) = 44.68960
        ak(20) = 37.52190
        ak(21) = 29.08949
        ak(22) = 20.84739
        ak(23) = 13.34443
        ak(24) = 7.084990
        ak(25) = 2.521360
        ak(26) = 0.0
        ak(27) = 0.0
 
        bk(1) = 0.0
        bk(2) = 0.0
        bk(3) = 0.0
        bk(4) = 0.0
        bk(5) = 0.0
        bk(6) = 0.0
        bk(7) = 0.0
        bk(8) = 0.0
        bk(9) = 0.0150530
        bk(10) = 0.0327622
        bk(11) = 0.0535962
        bk(12) = 0.0781062
        bk(13) = 0.1069411
        bk(14) = 0.1408637
        bk(15) = 0.1807720
        bk(16) = 0.2277220
        bk(17) = 0.2829562
        bk(18) = 0.3479364
        bk(19) = 0.4243822
        bk(20) = 0.5143168
        bk(21) = 0.6201202
        bk(22) = 0.7235355
        bk(23) = 0.8176768
        bk(24) = 0.8962153
        bk(25) = 0.9534761
        bk(26) = 0.9851122
        bk(27) = 1.
 
      end subroutine head_cam42
!
!-----------------------------------------------------------------------
! 
      subroutine cam42(idate,idate0,xlon,xlat,glat,jx,iy,i0,i1,j0,j1)
      use netcdf
      implicit none
!
! Dummy arguments
!
      integer :: i0 , i1 , idate , idate0 , iy , j0 , j1 , jx
      real , dimension(jlath) :: glat
      real , dimension(jx,iy) :: xlat , xlon
      intent (in) glat , idate , idate0 , iy , jx , xlat , xlon
      intent (inout) i0 , i1 , j0 , j1
!
! Local variables
!
      integer :: i , ii , ilev , inet , it , ivar , j , jj ,            &
               & jmax , jmin , kkrec , latid , latlen , lonid ,         &
               & lonlen , month , nday , nhour , nlev , nyear ,         &
               & istatus , timid , timlen
      character(25) :: inname
      real :: nlat0 , nlat1 , nlon0 , nlon1
      character(39) :: pathaddname
      logical :: there
      real , dimension(6) :: checklat , checklon , checktim
      character(2) , dimension(6) :: varname
      real , allocatable , dimension(:,:,:) :: work
      real , allocatable , dimension(:) :: work1 , work2
!
      data varname/'T' , 'Z3' , 'Q' , 'U' , 'V' , 'PS'/
!
      nyear = idate/1000000
      month = idate/10000 - nyear*100
      nday = idate/100 - nyear*10000 - month*100
      nhour = idate - nyear*1000000 - month*10000 - nday*100
      do kkrec = 1 , 6
        if ( kkrec<=5 ) nlev = klev
        if ( kkrec==6 ) nlev = 0
        if ( idate==idate0 .or.                                         &
           & ((mod(idate,100000)==10100) .and. mod(idate,1000000)       &
           & /=110100) ) then
          if ( kkrec==1 ) then
            write (inname,99001) nyear , 'air.' , nyear
          else if ( kkrec==2 ) then
            write (inname,99001) nyear , 'hgt.' , nyear
          else if ( kkrec==3 ) then
            write (inname,99002) nyear , 'shum.' , nyear
          else if ( kkrec==4 ) then
            write (inname,99002) nyear , 'uwnd.' , nyear
          else if ( kkrec==5 ) then
            write (inname,99002) nyear , 'vwnd.' , nyear
          else if ( kkrec==6 ) then
            write (inname,99002) nyear , 'pres.' , nyear
          else
          end if
 
          pathaddname = '../DATA/CAM42/'//inname
          inquire (file=pathaddname,exist=there)
          if ( .not.there ) then
            print * , pathaddname , ' is not available'
            stop
          end if
 
          istatus = nf90_open(pathaddname,nf90_nowrite,inet6(kkrec))
          if ( istatus/=nf90_noerr ) call handle_err(istatus)
          istatus = nf90_inq_varid(inet6(kkrec),varname(kkrec),         &
                  &                ivar6(kkrec))
          if ( istatus/=nf90_noerr ) call handle_err(istatus)
          write (*,*) inet6(kkrec) , pathaddname , ivar6(kkrec) ,       &
                    & varname(kkrec)
        end if
 
        istatus = nf90_inq_dimid(inet6(kkrec),'lon',lonid)
        if ( istatus/=nf90_noerr ) call handle_err(istatus)
        istatus = nf90_inq_dimid(inet6(kkrec),'lat',latid)
        if ( istatus/=nf90_noerr ) call handle_err(istatus)
        istatus = nf90_inq_dimid(inet6(kkrec),'time',timid)
        if ( istatus/=nf90_noerr ) call handle_err(istatus)
        istatus = nf90_inquire_dimension(inet6(kkrec),lonid,len=lonlen)
        if ( istatus/=nf90_noerr ) call handle_err(istatus)
        istatus = nf90_inquire_dimension(inet6(kkrec),latid,len=latlen)
        if ( istatus/=nf90_noerr ) call handle_err(istatus)
        istatus = nf90_inquire_dimension(inet6(kkrec),timid,len=timlen)
        if ( istatus/=nf90_noerr ) call handle_err(istatus)
 
 
        checklon(kkrec) = lonlen
        checklat(kkrec) = latlen
        checktim(kkrec) = timlen
        if ( kkrec>1 ) then
          if ( checklat(kkrec)/=checklat(kkrec-1) .or. checklon(kkrec)  &
             & /=checklon(kkrec-1) .or. checktim(kkrec)                 &
             & /=checktim(kkrec-1) ) then
            print * , 'DOMAIN DIMENSIONS DO NOT MATCH'
            print * , 'LAT for' , varname(kkrec+1) , '=' ,              &
                & checklat(kkrec)
            print * , 'LAT for' , varname(kkrec) , '=' , checklat(kkrec)
            print * , 'LON for' , varname(kkrec+1) , '=' ,              &
                & checklon(kkrec)
            print * , 'LON for' , varname(kkrec) , '=' , checklon(kkrec)
            print * , 'TIME for' , varname(kkrec+1) , '=' ,             &
                & checktim(kkrec)
            print * , 'TIME for' , varname(kkrec) , '=' ,               &
                & checktim(kkrec)
            stop 'Check CCSM Data Files'
          end if
        end if
        if ( kkrec==1 ) then
            allocate (work(lonlen,latlen,klev))
            allocate (work1(lonlen))
            allocate (work2(latlen))
        end if
        if ( xlon(1,1)<0.0 ) then
          nlon0 = xlon(1,1) + 360.
        else
          nlon0 = xlon(1,1)
        end if
        if ( xlon(jx,1)<0.0 ) then
          nlon1 = xlon(jx,1) + 360.
        else
          nlon1 = xlon(jx,1)
        end if
        nlat0 = xlat(1,1)
        nlat1 = xlat(1,iy)
        istart(1) = 1
        icount(1) = lonlen
        istatus = nf90_inq_varid(inet6(kkrec),'lon',lonid)
        if ( istatus/=nf90_noerr ) call handle_err(istatus)
        istatus = nf90_inq_varid(inet6(kkrec),'lat',latid)
        if ( istatus/=nf90_noerr ) call handle_err(istatus)
        istatus = nf90_get_var(inet6(kkrec),lonid,work1,istart,icount)
        icount(1) = latlen
        istatus = nf90_get_var(inet6(kkrec),latid,work2,istart,icount)
 
        if ( nlon0<work1(1) .or. nlon1>work1(lonlen) .or. nlat0<work2(1)&
           & .or. nlat1>work2(latlen) ) then
          print * , 'DOMAIN DIMENSIONS DO NOT MATCH'
          print * , 'CCSM Window LON min=' , work1(1) , 'max=' ,        &
              & work1(lonlen)
          print * , 'RCM  Domain LON min=' , nlon0 , 'max=' , nlon1
          print * , 'CCSM Window LAT min=' , work2(1) , 'max=' ,        &
              & work2(latlen)
          print * , 'RCM  Domain LAT min=' , nlat0 , 'max=' , nlat1
          stop 'Correct Domain Parameters in domain.param'
        end if
 
        if ( idate==idate0 ) then
          i0 = work1(1)/2.8125 + 1
          i1 = work1(lonlen)/2.8125 + 1
          do i = 1 , jlath
            jmin = nint(glat(i)-work2(1))
            jmax = nint(glat(i)-work2(latlen))
            if ( jmin==0 ) j0 = i
            if ( jmax==0 ) j1 = i
          end do
        end if
 
        it = (nday-1)*4 + nhour/6 + 1
        if ( month==2 ) it = it + 31*4
        if ( month==3 ) it = it + 59*4
        if ( month==4 ) it = it + 90*4
        if ( month==5 ) it = it + 120*4
        if ( month==6 ) it = it + 151*4
        if ( month==7 ) it = it + 181*4
        if ( month==8 ) it = it + 212*4
        if ( month==9 ) it = it + 243*4
        if ( month==10 ) it = it + 273*4
        if ( month==11 ) it = it + 304*4
        if ( month==12 ) it = it + 334*4
 
        icount(1) = lonlen
        icount(2) = latlen
        icount(3) = klev
        icount(4) = 1
        istart(1) = 1
        istart(2) = 1
        istart(3) = 1
        istart(4) = it
        inet = inet6(kkrec)
        ivar = ivar6(kkrec)
        if ( nlev>0 ) then
          icount(3) = nlev
          istatus = nf90_get_var(inet,ivar,work,istart,icount)
          do ilev = 1 , nlev
            if ( kkrec==1 ) then
              do jj = j0 , j1
                j = jj - j0 + 1
                if ( i0>i1 ) then
                  do ii = i0 , ilonh
                    i = ii - i0 + 1
                    tvar(ii,jj,ilev) = work(i,j,ilev)
                  end do
                  do ii = 1 , i1
                    i = ii + (ilonh-i0) + 1
                    tvar(ii,jj,ilev) = work(i,j,ilev)
                  end do
                else
                  do ii = i0 , i1
                    i = ii - i0 + 1
                    tvar(ii,jj,ilev) = work(i,j,ilev)
                  end do
                end if
              end do
!
            else if ( kkrec==2 ) then
              do jj = j0 , j1
                j = jj - j0 + 1
                if ( i0>i1 ) then
                  do ii = i0 , ilonh
                    i = ii - i0 + 1
                    hvar(ii,jj,ilev) = work(i,j,ilev)
                  end do
                  do ii = 1 , i1
                    i = ii + (ilonh-i0) + 1
                    hvar(ii,jj,ilev) = work(i,j,ilev)
                  end do
                else
                  do ii = i0 , i1
                    i = ii - i0 + 1
                    hvar(ii,jj,ilev) = work(i,j,ilev)
                  end do
                end if
              end do
            else if ( kkrec==3 ) then
              do jj = j0 , j1
                j = jj - j0 + 1
                if ( i0>i1 ) then
                  do ii = i0 , ilonh
                    i = ii - i0 + 1
                    qvar(ii,jj,ilev) = work(i,j,ilev)
                  end do
                  do ii = 1 , i1
                    i = ii + (ilonh-i0) + 1
                    qvar(ii,jj,ilev) = work(i,j,ilev)
                  end do
                else
                  do ii = i0 , i1
                    i = ii - i0 + 1
                    qvar(ii,jj,ilev) = work(i,j,ilev)
                  end do
                end if
              end do
            else if ( kkrec==4 ) then
              do jj = j0 , j1
                j = jj - j0 + 1
                if ( i0>i1 ) then
                  do ii = i0 , ilonh
                    i = ii - i0 + 1
                    uvar(ii,jj,ilev) = work(i,j,ilev)
                  end do
                  do ii = 1 , i1
                    i = ii + (ilonh-i0) + 1
                    uvar(ii,jj,ilev) = work(i,j,ilev)
                  end do
                else
                  do ii = i0 , i1
                    i = ii - i0 + 1
                    uvar(ii,jj,ilev) = work(i,j,ilev)
                  end do
                end if
              end do
            else if ( kkrec==5 ) then
              do jj = j0 , j1
                j = jj - j0 + 1
                if ( i0>i1 ) then
                  do ii = i0 , ilonh
                    i = ii - i0 + 1
                    vvar(ii,jj,ilev) = work(i,j,ilev)
                  end do
                  do ii = 1 , i1
                    i = ii + (ilonh-i0) + 1
                    vvar(ii,jj,ilev) = work(i,j,ilev)
                  end do
                else
                  do ii = i0 , i1
                    i = ii - i0 + 1
                    vvar(ii,jj,ilev) = work(i,j,ilev)
                  end do
                end if
              end do
            else
            end if
          end do
        else if ( nlev==0 ) then
          icount(3) = 1
          istart(3) = it
          istatus = nf90_get_var(inet,ivar,work,istart,icount)
          if ( kkrec==6 ) then
            do jj = j0 , j1
              j = jj - j0 + 1
              if ( i0>i1 ) then
                do ii = i0 , ilonh
                  i = ii - i0 + 1
                  psvar(ii,jj) = work(i,j,1)*0.01
                end do
                do ii = 1 , i1
                  i = ii + (ilonh-i0) + 1
                  psvar(ii,jj) = work(i,j,1)*0.01
                end do
              else
                do ii = i0 , i1
                  i = ii - i0 + 1
                  psvar(ii,jj) = work(i,j,1)*0.01
                end do
              end if
            end do
          end if
        else
        end if
      end do
      do kkrec = 1 , 6
        istatus = nf90_close(inet6(kkrec))
      end do
      deallocate(work)
      deallocate(work1)
      deallocate(work2)
99001 format (i4,'/','ccsm.',a4,i4,'.nc')
99002 format (i4,'/','ccsm.',a5,i4,'.nc')
      end subroutine cam42
!
!-----------------------------------------------------------------------
! 
      subroutine get_cam85(idate)

      use netcdf

      implicit none
!
! Dummy arguments
!
      integer :: idate
!
! Local variables
!
      real , dimension(npl) :: c1
      real , dimension(kz) :: c2
      integer :: checklat , checklon , i , i0 , i1 , ii , imax ,        &
               & imin , j , j0 , j1 , jj , k , latid ,  latlen , lonid ,&
               & lonlen , nmop , nyrp , istatus
      integer :: iomega
      real , dimension(jx,iy,npl) :: dum1
      real , dimension(jx,iy,npl,2) :: dum2
      real :: pp , w , wt
      real , dimension(jx,iy) :: pa , sst1 , sst2 , tlayer , za
      character(38) :: pathaddname
      logical :: there
      character(4) :: varname
      real , allocatable , dimension(:,:) :: work
!
      data varname/'PHIS'/
 
      call cam85(idate,idate1,xlon,xlat,glat,jx,iy,i0,i1,j0,j1)
      if ( idate==idate1 ) then
        pathaddname = '../DATA/CAM85/ccsm_ht.nc'
        inquire (file=pathaddname,exist=there)
        if ( .not.there ) then
          print * , pathaddname , 'is not available'
          stop
        end if
        istatus = nf90_open(pathaddname,nf90_nowrite,inet1)
        if ( istatus/=nf90_noerr ) call handle_err(istatus)
        istatus = nf90_inq_varid(inet1,varname,ivar1)
        if ( istatus/=nf90_noerr ) call handle_err(istatus)
        write (*,*) inet1 , ivar1 , pathaddname , ivar1 , varname
        istatus = nf90_inq_dimid(inet1,'lon',lonid)
        if ( istatus/=nf90_noerr ) call handle_err(istatus)
        istatus = nf90_inq_dimid(inet1,'lat',latid)
        if ( istatus/=nf90_noerr ) call handle_err(istatus)
        istatus = nf90_inquire_dimension(inet1,lonid,len=lonlen)
        if ( istatus/=nf90_noerr ) call handle_err(istatus)
        istatus = nf90_inquire_dimension(inet1,latid,len=latlen)
        if ( istatus/=nf90_noerr ) call handle_err(istatus)
        imin = i0
        imax = i1
        if ( i0>i1 ) imax = i1 + ilon
        checklon = (imax-imin) + 1
        checklat = (j1-j0) + 1
        if ( checklon/=lonlen .or. checklat/=latlen ) then
          print * , 'DOMAIN DIMENSIONS DO NOT MATCH'
          print * , 'LAT for 3D Variables = ' , checklat
          print * , 'LAT for' , varname , '= ' , latlen
          print * , 'LON for 3D Variables = ' , checklon
          print * , 'LON for' , varname , '= ' , lonlen
          stop 'Check Dimensions in CCSM Data Files'
        end if
        allocate (work(lonlen,latlen))
        icount(1) = lonlen
        icount(2) = latlen
        icount(3) = 1
        istart(1) = 1
        istart(2) = 1
        istart(3) = 1
        istatus = nf90_get_var(inet1,ivar1,work,istart,icount)
        if ( istatus/=nf90_noerr ) call handle_err(istatus)
        istatus = nf90_close(inet1)
        do jj = j0 , j1
          j = jj - j0 + 1
          if ( i0>i1 ) then
            do ii = i0 , ilon
              i = ii - i0 + 1
              zsvar(ii,jj) = work(i,j)/9.80665
            end do
            do ii = 1 , i1
              i = ii + (ilon-i0)
              zsvar(ii,jj) = work(i,j)/9.80665
            end do
          else
            do ii = i0 , i1
              i = ii - i0 + 1
              zsvar(ii,jj) = work(i,j)/9.80665
            end do
          end if
        end do
      end if
      write (*,*) 'Read in fields at Date: ' , idate
 
      do k = 1 , klev
        do j = 1 , jlat
          do i = 1 , ilon
            if ( psvar(i,j)>-9995. ) then
              pp3d(i,j,k) = psvar(i,j)*0.5*(ak(k)+bk(k+1))              &
                          & + 0.5*(ak(k)+bk(k+1))
            else
              pp3d(i,j,k) = -9999.0
            end if
          end do
        end do
      end do
 
      call height(hp,hvar,tvar,psvar,pp3d,zsvar,ilon,jlat,klev,pplev,   &
                & npl)
 
      call intlin(up,uvar,psvar,pp3d,ilon,jlat,klev,pplev,npl)
      call intlin(vp,vvar,psvar,pp3d,ilon,jlat,klev,pplev,npl)
 
      call intlog(tp,tvar,psvar,pp3d,ilon,jlat,klev,pplev,npl)
 
      call humid1fv(tvar,qvar,pp3d,ilon,jlat,klev)
      call intlin(qp,qvar,psvar,pp3d,ilon,jlat,klev,pplev,npl)
 
      call bilinxcx(b3,b2,xlon,xlat,glon,glat,jx,iy,ilon,jlat,npl)
      call bilinxdt(d3,d2,dlon,dlat,glon,glat,jx,iy,ilon,jlat,npl)
      call uvrot4(u3,v3,dlon,dlat,clon,clat,grdfac,jx,iy,npl,plon,plat, &
                & lgtype)
 
      call top2btm(t3,jx,iy,npl)
      call top2btm(q3,jx,iy,npl)
      call top2btm(h3,jx,iy,npl)
      call top2btm(u3,jx,iy,npl)
      call top2btm(v3,jx,iy,npl)
 
      call intgtb(pa,za,tlayer,topogm,t3,h3,sigmar,jx,iy,npl,dum1,dum2)
 
      call intpsn(ps4,topogm,pa,za,tlayer,ptop,jx,iy)
      call p1p2(b3pd,ps4,jx,iy)
 
      call intv3(ts4,t3,ps4,sigmar,ptop,jx,iy,npl,dum1)
      call camclndr(idate,nyrp,nmop,wt)
      call mksst(ts4,sst1,sst2,topogm,xlandu,jx,iy,nyrp,nmop,wt)
      call intv1(u4,u3,b3pd,sigma2,sigmar,ptop,jx,iy,kz,npl,1,1,dum2,c1,&
               & c2)
      call intv1(v4,v3,b3pd,sigma2,sigmar,ptop,jx,iy,kz,npl,1,1,dum2,c1,&
               & c2)
      call intv2(t4,t3,ps4,sigma2,sigmar,ptop,jx,iy,kz,npl,1,dum2,c1,c2)
 
      call intv1(q4,q3,ps4,sigma2,sigmar,ptop,jx,iy,kz,npl,1,0,dum2,c1, &
               & c2)
      call humid2fv(t4,q4,ps4,ptop,sigma2,jx,iy,kz)
 
      call hydrost(h4,t4,topogm,ps4,ptop,sigmaf,sigma2,dsigma,jx,iy,kz)
 
      call writef(u4,v4,t4,q4,pp,w,ps4,ts4,ptop,jx,iy,kz,idate,iomega)
 
      deallocate(work)
      end subroutine get_cam85
!
!-----------------------------------------------------------------------
! 
      subroutine head_cam85
      implicit none
!
! Local variables
!
      integer :: i , k , kr
      real, dimension(jlat) :: fixlat
!
      data fixlat/ - 88.928 , -87.539 , -86.141 , -84.742 , -83.343 ,   &
         & -81.942 , -80.543 , -79.142 , -77.742 , -76.341 , -74.941 ,  &
         & -73.539 , -72.139 , -70.739 , -69.338 , -67.937 , -66.536 ,  &
         & -65.136 , -63.735 , -62.334 , -60.934 , -59.534 , -58.132 ,  &
         & -56.731 , -55.331 , -53.930 , -52.529 , -51.128 , -49.728 ,  &
         & -48.327 , -46.926 , -45.525 , -44.125 , -42.724 , -41.323 ,  &
         & -39.922 , -38.522 , -37.121 , -35.720 , -34.319 , -32.919 ,  &
         & -31.518 , -30.117 , -28.716 , -27.315 , -25.915 , -24.514 ,  &
         & -23.113 , -21.712 , -20.312 , -18.911 , -17.510 , -16.109 ,  &
         & -14.709 , -13.308 , -11.907 , -10.506 , -9.105 , -7.705 ,    &
         & -6.304 , -4.903 , -3.502 , -2.102 , -0.701 , 0.701 , 2.102 , &
         & 3.502 , 4.903 , 6.304 , 7.705 , 9.105 , 10.506 , 11.907 ,    &
         & 13.308 , 14.709 , 16.109 , 17.510 , 18.911 , 20.312 ,        &
         & 21.712 , 23.113 , 24.514 , 25.915 , 27.315 , 28.716 ,        &
         & 30.117 , 31.518 , 32.919 , 34.319 , 35.720 , 37.121 ,        &
         & 38.522 , 39.922 , 41.323 , 42.724 , 44.125 , 45.525 ,        &
         & 46.926 , 48.327 , 49.728 , 51.128 , 52.529 , 53.930 ,        &
         & 55.331 , 56.731 , 58.132 , 59.533 , 60.934 , 62.334 ,        &
         & 63.735 , 65.136 , 66.536 , 67.937 , 69.338 , 70.739 ,        &
         & 72.139 , 73.539 , 74.941 , 76.341 , 77.742 , 79.142 ,        &
         & 80.543 , 81.942 , 83.343 , 84.742 , 86.141 , 87.539 , 88.928/
 
      do i = 1 , ilon
        glon(i) = float(i-1)*1.40625
      end do
      do i = 1 , jlat
        glat(i) = fixlat(i)
      end do
 
      pplev(1) = 30.
      pplev(2) = 50.
      pplev(3) = 70.
      pplev(4) = 100.
      pplev(5) = 150.
      pplev(6) = 200.
      pplev(7) = 250.
      pplev(8) = 300.
      pplev(9) = 350.
      pplev(10) = 420.
      pplev(11) = 500.
      pplev(12) = 600.
      pplev(13) = 700.
      pplev(14) = 780.
      pplev(15) = 850.
      pplev(16) = 920.
      pplev(17) = 960.
      pplev(18) = 1000.
 
!
      do k = 1 , npl
        sigmar(k) = pplev(k)*0.001
      end do
 
      do k = 1 , npl
        kr = npl - k + 1
        sigma1(k) = sigmar(kr)
      end do
      psref = 1000.
 
!     Ak is acatually Ak*Po
 
      ak(1) = 2.194067
      ak(2) = 4.895209
      ak(3) = 9.882418
      ak(4) = 18.05201
      ak(5) = 29.83724
      ak(6) = 44.62334
      ak(7) = 61.60587
      ak(8) = 78.51243
      ak(9) = 77.31271
      ak(10) = 75.90131
      ak(11) = 74.24086
      ak(12) = 72.28744
      ak(13) = 69.98932
      ak(14) = 67.28574
      ak(15) = 64.10509
      ak(16) = 60.36322
      ak(17) = 55.96111
      ak(18) = 50.78225
      ak(19) = 44.68960
      ak(20) = 37.52190
      ak(21) = 29.08949
      ak(22) = 20.84739
      ak(23) = 13.34443
      ak(24) = 7.084990
      ak(25) = 2.521360
      ak(26) = 0.0
      ak(27) = 0.0
 
      bk(1) = 0.0
      bk(2) = 0.0
      bk(3) = 0.0
      bk(4) = 0.0
      bk(5) = 0.0
      bk(6) = 0.0
      bk(7) = 0.0
      bk(8) = 0.0
      bk(9) = 0.0150530
      bk(10) = 0.0327622
      bk(11) = 0.0535962
      bk(12) = 0.0781062
      bk(13) = 0.1069411
      bk(14) = 0.1408637
      bk(15) = 0.1807720
      bk(16) = 0.2277220
      bk(17) = 0.2829562
      bk(18) = 0.3479364
      bk(19) = 0.4243822
      bk(20) = 0.5143168
      bk(21) = 0.6201202
      bk(22) = 0.7235355
      bk(23) = 0.8176768
      bk(24) = 0.8962153
      bk(25) = 0.9534761
      bk(26) = 0.9851122
      bk(27) = 1.
 
      end subroutine head_cam85
!
!-----------------------------------------------------------------------
! 
      subroutine cam85(idate,idate0,xlon,xlat,glat,jx,iy,i0,i1,j0,j1)
      use netcdf
      implicit none
!
! Dummy arguments
!
      integer :: i0 , i1 , idate , idate0 , iy , j0 , j1 , jx
      real , dimension(jlat) :: glat
      real , dimension(jx,iy) :: xlat , xlon
      intent (in) glat , idate , idate0 , iy , jx , xlat , xlon
      intent (inout) i0 , i1 , j0 , j1
!
! Local variables
!
      integer :: i , ii , ilev , inet , it , ivar , j , jj , jmax ,     &
               & jmin , kkrec , latid , latlen , lonid , lonlen ,       &
               & month , nday , nhour , nlev , nyear , istatus ,        &
               & timid , timlen
      character(25) :: inname
      real :: nlat0 , nlat1 , nlon0 , nlon1
      character(3) , dimension(12) :: nmonth
      character(39) :: pathaddname
      logical :: there
      character(2) , dimension(6) :: varname
      real , allocatable , dimension(:,:,:) :: work
      real , allocatable , dimension(:) :: work1 , work2
      real , dimension(6) :: checklat , checklon , checktim
!
      data nmonth/'JAN' , 'FEB' , 'MAR' , 'APR' , 'MAY' , 'JUN' ,       &
         & 'JUL' , 'AUG' , 'SEP' , 'OCT' , 'NOV' , 'DEC'/
 
      data varname/'T' , 'Z3' , 'Q' , 'U' , 'V' , 'PS'/
 
      nyear = idate/1000000
      month = idate/10000 - nyear*100
      nday = idate/100 - nyear*10000 - month*100
      nhour = idate - nyear*1000000 - month*10000 - nday*100
      do kkrec = 1 , 6
        if ( kkrec<=5 ) nlev = klev
        if ( kkrec==6 ) nlev = 0
        if ( idate==idate0 .or.                                         &
           & ((mod(idate,100000)==10100 .and. mod(idate,1000000)        &
           & /=110100)) .or. mod(idate,10000)==100 ) then
          if ( kkrec==1 ) then
            write (inname,99001) nyear , 'air' , nmonth(month) , nyear
          else if ( kkrec==2 ) then
            write (inname,99001) nyear , 'hgt' , nmonth(month) , nyear
          else if ( kkrec==3 ) then
            write (inname,99002) nyear , 'shum' , nmonth(month) , nyear
          else if ( kkrec==4 ) then
            write (inname,99002) nyear , 'uwnd' , nmonth(month) , nyear
          else if ( kkrec==5 ) then
            write (inname,99002) nyear , 'vwnd' , nmonth(month) , nyear
          else if ( kkrec==6 ) then
            write (inname,99002) nyear , 'pres' , nmonth(month) , nyear
          else
          end if
 
          pathaddname = '../DATA/CAM85/'//inname
          inquire (file=pathaddname,exist=there)
          if ( .not.there ) then
            print * , pathaddname , ' is not available'
            stop
          end if
 
          istatus = nf90_open(pathaddname,nf90_nowrite,inet6(kkrec))
          if ( istatus/=nf90_noerr ) call handle_err(istatus)
          istatus = nf90_inq_varid(inet6(kkrec),varname(kkrec),         &
                 &                 ivar6(kkrec))
          if ( istatus/=nf90_noerr ) call handle_err(istatus)
          write (*,*) inet6(kkrec) , pathaddname , ivar6(kkrec) ,       &
                    & varname(kkrec)
        end if
 
        istatus = nf90_inq_dimid(inet6(kkrec),'lon',lonid)
        if ( istatus/=nf90_noerr ) call handle_err(istatus)
        istatus = nf90_inq_dimid(inet6(kkrec),'lat',latid)
        if ( istatus/=nf90_noerr ) call handle_err(istatus)
        istatus = nf90_inq_dimid(inet6(kkrec),'time',timid)
        if ( istatus/=nf90_noerr ) call handle_err(istatus)
        istatus = nf90_inquire_dimension(inet6(kkrec),lonid,len=lonlen)
        if ( istatus/=nf90_noerr ) call handle_err(istatus)
        istatus = nf90_inquire_dimension(inet6(kkrec),latid,len=latlen)
        if ( istatus/=nf90_noerr ) call handle_err(istatus)
        istatus = nf90_inquire_dimension(inet6(kkrec),timid,len=timlen)
        if ( istatus/=nf90_noerr ) call handle_err(istatus)
 
        checklon(kkrec) = lonlen
        checklat(kkrec) = latlen
        checktim(kkrec) = timlen
        if ( kkrec>1 ) then
          if ( checklat(kkrec)/=checklat(kkrec-1) .or. checklon(kkrec)  &
             & /=checklon(kkrec-1) .or. checktim(kkrec)                 &
             & /=checktim(kkrec-1) ) then
            print * , 'DOMAIN DIMENSIONS DO NOT MATCH'
            print * , 'LAT for' , varname(kkrec+1) , '=' ,              &
                & checklat(kkrec)
            print * , 'LAT for' , varname(kkrec) , '=' , checklat(kkrec)
            print * , 'LON for' , varname(kkrec+1) , '=' ,              &
                & checklon(kkrec)
            print * , 'LON for' , varname(kkrec) , '=' , checklon(kkrec)
            print * , 'TIME for' , varname(kkrec+1) , '=' ,             &
                & checktim(kkrec)
            print * , 'TIME for' , varname(kkrec) , '=' ,               &
                & checktim(kkrec)
            stop 'Check CCSM Data Files'
          end if
        end if
        if ( kkrec==1 ) then
            allocate (work(lonlen,latlen,klev))
            allocate (work1(lonlen))
            allocate (work2(latlen))
        end if
        if ( xlon(1,1)<0.0 ) then
          nlon0 = xlon(1,1) + 360.
        else
          nlon0 = xlon(1,1)
        end if
        if ( xlon(jx,1)<0.0 ) then
          nlon1 = xlon(jx,1) + 360.
        else
          nlon1 = xlon(jx,1)
        end if
        nlat0 = xlat(1,1)
        nlat1 = xlat(1,iy)
        istart(1) = 1
        icount(1) = lonlen
        istatus = nf90_inq_varid(inet6(kkrec),'lon',lonid)
        if ( istatus/=nf90_noerr ) call handle_err(istatus)
        istatus = nf90_inq_varid(inet6(kkrec),'lat',latid)
        if ( istatus/=nf90_noerr ) call handle_err(istatus)
        istatus = nf90_get_var(inet6(kkrec),lonid,work1,istart,icount)
        icount(1) = latlen
        istatus = nf90_get_var(inet6(kkrec),latid,work2,istart,icount)
 
        if ( nlon0<work1(1) .or. nlon1>work1(lonlen) .or. nlat0<work2(1)&
           & .or. nlat1>work2(latlen) ) then
          print * , 'DOMAIN DIMENSIONS DO NOT MATCH'
          print * , 'CCSM Window LON min=' , work1(1) , 'max=' ,        &
              & work1(lonlen)
          print * , 'RCM  Domain LON min=' , nlon0 , 'max=' , nlon1
          print * , 'CCSM Window LAT min=' , work2(1) , 'max=' ,        &
              & work2(latlen)
          print * , 'RCM  Domain LAT min=' , nlat0 , 'max=' , nlat1
          stop 'Correct Domain Parameters in domain.param'
        end if
 
        if ( idate==idate0 ) then
          i0 = work1(1)/1.40625 + 1
          i1 = work1(lonlen)/1.40625 + 1
          do i = 1 , jlat
            jmin = nint(glat(i)-work2(1))
            jmax = nint(glat(i)-work2(latlen))
            if ( jmin==0 ) j0 = i
            if ( jmax==0 ) j1 = i
          end do
        end if
 
        it = (nday-1)*4 + nhour/6 + 1
        if ( month==2 ) it = it + 31*4
        if ( month==3 ) it = it + 59*4
        if ( month==4 ) it = it + 90*4
        if ( month==5 ) it = it + 120*4
        if ( month==6 ) it = it + 151*4
        if ( month==7 ) it = it + 181*4
        if ( month==8 ) it = it + 212*4
        if ( month==9 ) it = it + 243*4
        if ( month==10 ) it = it + 273*4
        if ( month==11 ) it = it + 304*4
        if ( month==12 ) it = it + 334*4
 
        icount(1) = lonlen
        icount(2) = latlen
        icount(3) = klev
        icount(4) = 1
        istart(1) = 1
        istart(2) = 1
        istart(3) = 1
        istart(4) = it
        inet = inet6(kkrec)
        ivar = ivar6(kkrec)
        if ( nlev>0 ) then
          icount(3) = nlev
          istatus = nf90_get_var(inet,ivar,work,istart,icount)
          do ilev = 1 , nlev
            if ( kkrec==1 ) then
              do jj = j0 , j1
                j = jj - j0 + 1
                if ( i0>i1 ) then
                  do ii = i0 , ilon
                    i = ii - i0 + 1
                    tvar(ii,jj,ilev) = work(i,j,ilev)
                  end do
                  do ii = 1 , i1
                    i = ii + (ilon-i0) + 1
                    tvar(ii,jj,ilev) = work(i,j,ilev)
                  end do
                else
                  do ii = i0 , i1
                    i = ii - i0 + 1
                    tvar(ii,jj,ilev) = work(i,j,ilev)
                  end do
                end if
              end do
            else if ( kkrec==2 ) then
              do jj = j0 , j1
                j = jj - j0 + 1
                if ( i0>i1 ) then
                  do ii = i0 , ilon
                    i = ii - i0 + 1
                    hvar(ii,jj,ilev) = work(i,j,ilev)
                  end do
                  do ii = 1 , i1
                    i = ii + (ilon-i0) + 1
                    hvar(ii,jj,ilev) = work(i,j,ilev)
                  end do
                else
                  do ii = i0 , i1
                    i = ii - i0 + 1
                    hvar(ii,jj,ilev) = work(i,j,ilev)
                  end do
                end if
              end do
            else if ( kkrec==3 ) then
              do jj = j0 , j1
                j = jj - j0 + 1
                if ( i0>i1 ) then
                  do ii = i0 , ilon
                    i = ii - i0 + 1
                    qvar(ii,jj,ilev) = work(i,j,ilev)
                  end do
                  do ii = 1 , i1
                    i = ii + (ilon-i0) + 1
                    qvar(ii,jj,ilev) = work(i,j,ilev)
                  end do
                else
                  do ii = i0 , i1
                    i = ii - i0 + 1
                    qvar(ii,jj,ilev) = work(i,j,ilev)
                  end do
                end if
              end do
            else if ( kkrec==4 ) then
              do jj = j0 , j1
                j = jj - j0 + 1
                if ( i0>i1 ) then
                  do ii = i0 , ilon
                    i = ii - i0 + 1
                    uvar(ii,jj,ilev) = work(i,j,ilev)
                  end do
                  do ii = 1 , i1
                    i = ii + (ilon-i0) + 1
                    uvar(ii,jj,ilev) = work(i,j,ilev)
                  end do
                else
                  do ii = i0 , i1
                    i = ii - i0 + 1
                    uvar(ii,jj,ilev) = work(i,j,ilev)
                  end do
                end if
              end do
            else if ( kkrec==5 ) then
              do jj = j0 , j1
                j = jj - j0 + 1
                if ( i0>i1 ) then
                  do ii = i0 , ilon
                    i = ii - i0 + 1
                    vvar(ii,jj,ilev) = work(i,j,ilev)
                  end do
                  do ii = 1 , i1
                    i = ii + (ilon-i0) + 1
                    vvar(ii,jj,ilev) = work(i,j,ilev)
                  end do
                else
                  do ii = i0 , i1
                    i = ii - i0 + 1
                    vvar(ii,jj,ilev) = work(i,j,ilev)
                  end do
                end if
              end do
            else
            end if
          end do
 
        else if ( nlev==0 ) then
          icount(3) = 1
          istart(3) = it
          istatus = nf90_get_var(inet,ivar,work,istart,icount)
          if ( kkrec==6 ) then
            do jj = j0 , j1
              j = jj - j0 + 1
              if ( i0>i1 ) then
                do ii = i0 , ilon
                  i = ii - i0 + 1
                  psvar(ii,jj) = work(i,j,1)*0.01
                end do
                do ii = 1 , i1
                  i = ii + (ilon-i0) + 1
                  psvar(ii,jj) = work(i,j,1)*0.01
                end do
              else
                do ii = i0 , i1
                  i = ii - i0 + 1
                  psvar(ii,jj) = work(i,j,1)*0.01
                end do
              end if
            end do
          end if
        else
        end if
      end do
      do kkrec = 1 , 6
        istatus = nf90_close(inet6(kkrec))
      end do
      deallocate (work)
      deallocate (work1)
      deallocate (work2)
99001 format (i4,'/','ccsm.',a3,a3,'.',i4,'.nc')
99002 format (i4,'/','ccsm.',a4,a3,'.',i4,'.nc')
      end subroutine cam85
!
!-----------------------------------------------------------------------
!
      subroutine camclndr(mdate,nyrp,nmop,wt)
        implicit none
!
! Dummy arguments
!
        integer :: mdate , nmop , nyrp
        real :: wt
        intent (in) mdate
        intent (out) nyrp , wt
        intent (inout) nmop
!
! Local variables
!
        real :: fdemon , fnumer
        integer :: idate , iday , imo , iyr , j , julday , nmo , nyr
        integer , dimension(12) :: jprev , julmid , lenmon , midmon
!
        data lenmon/31 , 28 , 31 , 30 , 31 , 30 , 31 , 31 , 30 , 31 ,   &
                  & 30 , 31/
        data midmon/16 , 14 , 16 , 15 , 16 , 15 , 16 , 16 , 15 , 16 ,   &
                  & 15 , 16/

!       Initialize nmop, nyrp

        nmop = 1
        nyrp = 0
 
        idate = mdate/100
        iyr = idate/10000
        imo = (idate-iyr*10000)/100
        iday = mod(idate,100)
 
        jprev(1) = 0
        do j = 2 , 12
          jprev(j) = jprev(j-1) + lenmon(j-1)
        end do
        do j = 1 , 12
          julmid(j) = jprev(j) + midmon(j)
        end do
        julday = iday + jprev(imo)
 
        do nyr = 1000 , 2100
          do nmo = 1 , 12
            if ( (nyr==iyr) .and. (julmid(nmo)>julday) ) go to 100
            if ( nyr>iyr ) go to 100
            nmop = nmo
            nyrp = nyr
          end do
        end do
 
 100    continue
        fnumer = float(julday-julmid(nmop))
        if ( fnumer<0. ) fnumer = fnumer + 365.
        fdemon = float(julmid(nmo)-julmid(nmop))
        if ( fdemon<=0. ) fdemon = fdemon + 365.
        wt = fnumer/fdemon
 
      end subroutine camclndr
!
!     Error Handler for NETCDF Calls
!
      subroutine handle_err(istatus)
        use netcdf
        implicit none
!
! Dummy arguments
!
        integer :: istatus
        intent (in) :: istatus
!
        print * , nf90_strerror(istatus)
        stop 'Netcdf Error'

      end subroutine handle_err
!
      end module mod_ccsm

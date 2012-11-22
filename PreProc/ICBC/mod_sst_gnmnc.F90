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

module mod_sst_gnmnc

  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_dynparam
  use mod_memutil
  use mod_sst_grid
  use mod_interp
  use mod_message
  use mod_nchelper
  use netcdf

  private

  integer(ik4) :: ilon , jlat 
!
  integer(ik4) :: inet1
  integer(ik4) , dimension(2) :: ivar2
  integer(ik4) :: timlen
  integer(ik4) :: timid
  integer(ik4) :: istatus
  integer(ik4) , dimension(3) :: istart , icount
  integer(ik4) , dimension(1) :: istart_t , icount_t
  real(rk8) , pointer ::  work1(:)
  real(rk8) , pointer , dimension (:, :) :: work2 , work3
  real(rk8) , pointer , dimension(:,:) :: sst
  type(rcm_time_and_date) , save :: fidate1
  character(len=64) :: cunit , ccal
  character(len=256) :: inpfile
  character(len=8), dimension(2) :: varname
!
  data varname/'time', 'TOBESET'/
!
  public :: sst_gnmnc

  contains

  subroutine sst_gnmnc
!
!*******************************************************************************
!
! This is a package of subroutines to read SST data on regular latlon grid in
! NETCDF format and interpolate on RCM grid.
! Need to configure datapath and variable names.
!
! Written By Moetasim Ashfaq Dec-2005 @ PURDUE.EDU for POP dataset prepared
! on the regular 1x1 degree resolution using NCL+NCO programs.
!
! Modified to be generic reader for monthly netCDF by Graziano Giuliani 2011
!
!******************************************************************************

  implicit none
!
  real(rk8) , pointer , dimension(:) :: glat
  real(rk8) , pointer , dimension(:) :: glon
  real(rk8) , pointer , dimension(:,:) :: glat2
  real(rk8) , pointer , dimension(:,:) :: glon2
  type(rcm_time_and_date) :: idate , idatef , idateo
  integer(ik4) :: i , j , k , nsteps , latid , lonid
  integer(ik4) :: year , month , day , hour , y1 , y2
  real(rk8) :: ufac
!
  call split_idate(globidate2, year, month, day, hour)  
!
  ufac = 0.0
  if ( ssttyp == "CAM4N" ) then
    inpfile = trim(inpglob)//'/SST/sst_HadOIBl_bc_0.9x1.25_1870_2008_c091020.nc'
    varname(2) = 'SST_cpl'
    ufac = 273.15
  else if ( ssttyp == "CCSST" ) then
    inpfile = trim(inpglob)//'/SST/ccsm_mn.sst.nc'
    varname(2) = 'SST'
    ufac = 273.15
  else if ( ssttyp == "CA_RF" ) then
    inpfile = trim(inpglob)//'/SST/ts_Amon_CanESM2_historical_r1i1p1_185001-200512.nc'
    varname(2) = 'ts'
  else if ( ssttyp == "CA_26" ) then
    inpfile = trim(inpglob)//'/SST/ts_Amon_CanESM2_rcp26_r1i1p1_200601-210012.nc'
    varname(2) = 'ts'
  else if ( ssttyp == "CA_45" ) then
    inpfile = trim(inpglob)//'/SST/ts_Amon_CanESM2_rcp45_r1i1p1_200601-210012.nc'
    varname(2) = 'ts'
  else if ( ssttyp == "CA_85" ) then
    inpfile = trim(inpglob)//'/SST/ts_Amon_CanESM2_rcp85_r1i1p1_200601-210012.nc'
    varname(2) = 'ts'
  else if ( ssttyp == "HA_RF" ) then
   if ( year < 1959 ) then
     inpfile = trim(inpglob)// &
           '/SST/tos_Omon_HadGEM2-ES_historical_r1i1p1_185912-195911.nc'
   else
     inpfile = trim(inpglob)// &
           '/SST/tos_Omon_HadGEM2-ES_historical_r1i1p1_195912-200512.nc'
   end if
    varname(2) = 'tos'
  else if ( ssttyp == "HA_26" ) then
    inpfile = trim(inpglob)//'/SST/tos_Omon_HadGEM2-ES_rcp26_r1i1p1_200512-209911.nc'
    varname(2) = 'tos'
  else if ( ssttyp == "HA_45" ) then
    inpfile = trim(inpglob)//'/SST/tos_Omon_HadGEM2-ES_rcp45_r1i1p1_200512-209911.nc'
    varname(2) = 'tos'
  else if ( ssttyp == "HA_85" ) then
    inpfile = trim(inpglob)//'/SST/tos_Omon_HadGEM2-ES_rcp85_r1i1p1_200512-209912.nc'
    varname(2) = 'tos'
  else if ( ssttyp == "CS_RF" ) then
   inpfile = trim(inpglob)// &
         '/SST/tos_Omon_CSIRO-Mk3-6-0_historical_r1i1p1_185001-200512.nc'
    varname(2) = 'tos'
  else if ( ssttyp == "CS_26" ) then
    inpfile = trim(inpglob)//'/SST/tos_Omon_CSIRO-Mk3-6-0_rcp26_r1i1p1_200601-210012.nc'
    varname(2) = 'tos'
  else if ( ssttyp == "CS_45" ) then
    inpfile = trim(inpglob)//'/SST/tos_Omon_CSIRO-Mk3-6-0_rcp45_r1i1p1_200601-210012.nc'
    varname(2) = 'tos'
  else if ( ssttyp == "CS_85" ) then
    inpfile = trim(inpglob)//'/SST/tos_Omon_CSIRO-Mk3-6-0_rcp85_r1i1p1_200601-210012.nc'
    varname(2) = 'tos'
  else if ( ssttyp == "EC_RF" ) then
    inpfile = trim(inpglob)//'/SST/EC-EARTH/RF/ich1_sst_1950-2009.nc'
    varname(2) = 'sst'
  else if ( ssttyp == "EC_45" ) then
    inpfile = trim(inpglob)//'/SST/EC-EARTH/RCP45/ic41_sst_2006-2100.nc'
    varname(2) = 'sst'
  else if ( ssttyp == "EC_85" ) then
    inpfile = trim(inpglob)//'/SST/EC-EARTH/RCP85/ic81_sst_2006-2100.nc'
    varname(2) = 'sst'
  else if ( ssttyp == 'IP_RF' ) then
    inpfile = trim(inpglob)//'/SST/tos_Omon_IPSL-CM5A-LR_historical_r1i1p1_185001-200512.nc'
    varname(2) = 'tos'
  else if ( ssttyp == 'IP_45' ) then
    inpfile = trim(inpglob)//'/SST/tos_Omon_IPSL-CM5A-LR_rcp45_r1i1p1_200601-230012.nc'
    varname(2) = 'tos'
  else if ( ssttyp == 'IP_85' ) then
    inpfile = trim(inpglob)//'/SST/tos_Omon_IPSL-CM5A-LR_rcp85_r1i1p1_200601-230012.nc'
    varname(2) = 'tos'
  else if ( ssttyp(1:3) == 'GF_' ) then
    call split_idate(globidate1, year, month, day, hour)  
    y1 = (year-1)/5*5+1
    y2 = y1+4
    if ( year == y1 .and. month == 1 .and. day == 1 .and. hour == 0 ) then
      y1 = y1 - 5
      y2 = y1 + 4
    end if
    if ( ssttyp(4:5) == 'RF' ) then
      write(inpfile,'(a,i4,a,i4,a)') &
         trim(inpglob)//'/SST/ts_Amon_GFDL-ESM2M_historical_r1i1p1_', &
         y1 , '01-', y2, '12.nc'
    else if ( ssttyp(4:5) == '26' ) then
      write(inpfile,'(a,i4,a,i4,a)') &
         trim(inpglob)//'/SST/ts_Amon_GFDL-ESM2M_rcp26_r1i1p1_', &
         y1 , '01-', y2, '12.nc'
    else if ( ssttyp(4:5) == '45' ) then
      write(inpfile,'(a,i4,a,i4,a)') &
         trim(inpglob)//'/SST/ts_Amon_GFDL-ESM2M_rcp45_r1i1p1_', &
         y1 , '01-', y2, '12.nc'
    else if ( ssttyp(4:5) == '85' ) then
      write(inpfile,'(a,i4,a,i4,a)') &
         trim(inpglob)//'/SST/ts_Amon_GFDL-ESM2M_rcp85_r1i1p1_', &
         y1 , '01-', y2, '12.nc'
    else
      call die('gnmnc_sst','Unknown ssttyp: '//ssttyp,1)
    end if
    varname(2) = 'ts'
  else if ( ssttyp(1:3) == 'CN_' ) then
    call split_idate(globidate1, year, month, day, hour)  
    if ( year < 2006 ) then
      y1 = (year)/10*10
      if ( y1 == 2000 ) then
        y2 = 2005
      else
        y2 = y1 + 9
      end if
    else
      y1 = (year-6)/10*10+6
      y2 = y1 + 9
    end if
    if ( year /= 2006 .and. month /= 1 .and. day /= 1 .and. hour /= 0 ) then
      if ( year == y1 .and. month == 1 .and. day == 1 .and. hour == 0 ) then
        y1 = y1 - 10
        y2 = y1 + 9
      end if
    end if
    if ( ssttyp(4:5) == 'RF' ) then
      write(inpfile,'(a,i4,a,i4,a)') &
         trim(inpglob)//'/SST/tos_Omon_CNRM-CM5_historical_r1i1p1_', &
         y1 , '01-', y2, '12.nc'
    else if ( ssttyp(4:5) == '26' ) then
      write(inpfile,'(a,i4,a,i4,a)') &
         trim(inpglob)//'/SST/tos_Omon_CNRM-CM5_rcp26_r1i1p1_', &
         y1 , '01-', y2, '12.nc'
    else if ( ssttyp(4:5) == '45' ) then
      write(inpfile,'(a,i4,a,i4,a)') &
         trim(inpglob)//'/SST/tos_Omon_CNRM-CM5_rcp45_r1i1p1_', &
         y1 , '01-', y2, '12.nc'
    else if ( ssttyp(4:5) == '85' ) then
      write(inpfile,'(a,i4,a,i4,a)') &
         trim(inpglob)//'/SST/tos_Omon_CNRM-CM5_rcp85_r1i1p1_', &
         y1 , '01-', y2, '12.nc'
    else
      call die('gnmnc_sst','Unknown ssttyp: '//ssttyp,1)
    end if
    varname(2) = 'tos'
  else
    call die('gnmnc_sst','Unknown ssttyp: '//ssttyp,1)
  end if

  istatus = nf90_open(inpfile,nf90_nowrite,inet1)
  call checkncerr(istatus,__FILE__,__LINE__,'Error opening '//trim(inpfile))
  write (stdout,*) inet1 , trim(inpfile)

! GET DIMENSION IDs
  if ( ssttyp(1:3) /= 'IP_' .and. ssttyp(1:3) /= 'CN_' ) then
    istatus = nf90_inq_dimid(inet1,'lat',latid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find dim lat')
    istatus = nf90_inq_dimid(inet1,'lon',lonid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find dim lon')
    istatus = nf90_inquire_dimension(inet1,latid,len=jlat)
    call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dim lat')
    istatus = nf90_inquire_dimension(inet1,lonid,len=ilon)
    call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dim lon')
  else
    istatus = nf90_inq_dimid(inet1,'j',latid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find dim j')
    istatus = nf90_inq_dimid(inet1,'i',lonid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find dim i')
    istatus = nf90_inquire_dimension(inet1,latid,len=jlat)
    call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dim j')
    istatus = nf90_inquire_dimension(inet1,lonid,len=ilon)
    call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dim i')
  end if
  istatus = nf90_inq_dimid(inet1,'time',timid)
  call checkncerr(istatus,__FILE__,__LINE__,'Error find dim time')

! GET DIMENSION LENGTHS
  istatus = nf90_inquire_dimension(inet1,timid,len=timlen)
  call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dim time')

  call getmem1d(work1,1,timlen,'mod_gnmnc_sst:work1')

! GET VARIABLE IDs
  istatus = nf90_inq_varid(inet1,'lat',latid)
  call checkncerr(istatus,__FILE__,__LINE__,'Error find var lat')
  istatus = nf90_inq_varid(inet1,'lon',lonid)
  call checkncerr(istatus,__FILE__,__LINE__,'Error find var lon')
  istatus = nf90_inq_varid(inet1,varname(1),ivar2(1))
  call checkncerr(istatus,__FILE__,__LINE__,'Error find var '//varname(1))
  istatus = nf90_inq_varid(inet1,varname(2),ivar2(2))
  call checkncerr(istatus,__FILE__,__LINE__,'Error find var '//varname(2))

  if ( ssttyp(1:3) /= 'CA_' .and. ssttyp(1:3) /= 'CN_' .and. &
       ssttyp(1:3) /= 'CS_' .and. ssttyp(1:3) /= 'GF_' .and. &
       ssttyp(1:3) /= 'IP_' .and. ssttyp(1:3) /= 'EC_' .and. &
       ssttyp(1:3) /= 'HA_' ) then
    call getmem1d(glat,1,jlat,'mod_gnmnc_sst:glat')
    call getmem1d(glon,1,ilon,'mod_gnmnc_sst:glon')
!   GET LATITUDE AND LONGITUDE
    istatus = nf90_get_var(inet1,latid,glat)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var lat')
    istatus = nf90_get_var(inet1,lonid,glon)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var lon')
  else
    if ( ssttyp(1:3) /= 'CN_' .and. ssttyp(1:3) /= 'IP_' ) then
      call getmem1d(glat,1,jlat,'mod_gnmnc_sst:glat')
      call getmem1d(glon,1,ilon,'mod_gnmnc_sst:glon')
      call getmem2d(glat2,1,ilon,1,jlat,'mod_gnmnc_sst:glat2')
      call getmem2d(glon2,1,ilon,1,jlat,'mod_gnmnc_sst:glon2')
!     GET LATITUDE AND LONGITUDE
      istatus = nf90_get_var(inet1,latid,glat)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read var lat')
      istatus = nf90_get_var(inet1,lonid,glon)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read var lon')
      do j = 1 , ilon
        glat2(j,:) = glat(:)
      end do
      do i = 1 , jlat
        glon2(:,i) = glon(:)
      end do
      where (glon2 >= 180.0)
        glon2 = glon2-360.0
      end where
    else
      call getmem2d(glat2,1,ilon,1,jlat,'mod_gnmnc_sst:glat2')
      call getmem2d(glon2,1,ilon,1,jlat,'mod_gnmnc_sst:glon2')
!     GET LATITUDE AND LONGITUDE
      istatus = nf90_get_var(inet1,latid,glat2)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read var lat')
      istatus = nf90_get_var(inet1,lonid,glon2)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read var lon')
      where (glon2 >= 180.0)
        glon2 = glon2-360.0
      end where
    end if
  end if
  call getmem2d(work2,1,ilon,1,jlat,'mod_gnmnc_sst:work2')
  call getmem2d(work3,1,ilon,1,jlat,'mod_gnmnc_sst:work3')
  call getmem2d(sst,1,ilon,1,jlat,'mod_gnmnc_sst:sst')
  
! GET TIME VALUES
  istart_t(1) = 1
  icount_t(1) = timlen
  istatus = nf90_get_var(inet1,ivar2(1),work1,istart_t,icount_t)
  call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(1))
! CHECK FOR THE REQUIRED RECORD IN DATA FILE  
  istatus = nf90_get_att(inet1,ivar2(1),'units',cunit)
  call checkncerr(istatus,__FILE__,__LINE__, &
                  'Error read var '//varname(1)//' units')
  istatus = nf90_get_att(inet1,ivar2(1),'calendar',ccal)
  call checkncerr(istatus,__FILE__,__LINE__, &
                  'Error read var '//varname(1)//' calendar')
  fidate1 = timeval2date(work1(1),cunit,ccal)

  idateo = monfirst(globidate1)
  idatef = monfirst(globidate2)
  if (idatef < globidate2) then
    idatef = nextmon(idatef)
  end if
  nsteps = imondiff(idatef,idateo) + 1
 
  call open_sstfile(idateo)
 
  idate = idateo
  do k = 1 , nsteps

    call gnmnc_sst(idate)
    if ( ssttyp(1:3) == 'CA_' .or. ssttyp(1:3) == 'CN_' .or. &
         ssttyp(1:3) == 'CS_' .or. ssttyp(1:3) == 'GF_' .or. &
         ssttyp(1:3) == 'IP_' .or. ssttyp(1:3) == 'EC_' .or. &
         ssttyp(1:3) == 'HA_' ) then
      call distwgtcr(sstmm,sst,xlon,xlat,glon2,glat2,jx,iy,ilon,jlat)
    else
      call bilinx(sst,sstmm,xlon,xlat,glon,glat,ilon,jlat,jx,iy,1)
    end if

    do i = 2 , iy-1
      do j = 2 , jx-1
        sstmm(j,i) = sstmm(j,i) + ufac
      end do
    end do

    call writerec(idate)
    write (stdout,*) 'WRITEN OUT SST DATA : ' , tochar(idate)

    idate = nextmon(idate)

  end do
 

  end subroutine sst_gnmnc
!
!-----------------------------------------------------------------------
!
!     Subroutine to read required records from SST data file
!
  subroutine gnmnc_sst(idate)

  implicit none

  type(rcm_time_and_date) , intent (in) :: idate
  real(rk8) :: wt1 , wt2
  type(rcm_time_and_date) :: prev , next

  integer(ik4) :: it , i , j
  integer(ik4) :: year , month , day , hour , y1 , y2
  type(rcm_time_interval) :: tdiff1 , tdiff2

  icount(1) = ilon
  icount(2) = jlat
  icount(3) = 1
  istart(1) = 1
  istart(2) = 1
  istart(3) = 1

  it = imondiff(idate,fidate1) + 1
  if ( ssttyp(1:3) == 'GF_' ) then
    if ( it > timlen ) then
      istart(3) = it-1
      istatus = nf90_get_var(inet1,ivar2(2),work2,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(2))
      prev = timeval2date(work1(it-1),cunit,ccal)
      ! Try switching to next file
      istatus = nf90_close(inet1)
      call checkncerr(istatus,__FILE__,__LINE__,'Error Close file')
      call split_idate(idate, year, month, day, hour)  
      y1 = (year-1)/5*5+1
      y2 = y1+4
      if ( ssttyp(4:5) == 'RF' ) then
        write(inpfile,'(a,i4,a,i4,a)') &
           trim(inpglob)//'/SST/ts_Amon_GFDL-ESM2M_historical_r1i1p1_', &
           y1 , '01-', y2, '12.nc'
      else if ( ssttyp(4:5) == '26' ) then
        write(inpfile,'(a,i4,a,i4,a)') &
           trim(inpglob)//'/SST/ts_Amon_GFDL-ESM2M_rcp26_r1i1p1_', &
           y1 , '01-', y2, '12.nc'
      else if ( ssttyp(4:5) == '45' ) then
        write(inpfile,'(a,i4,a,i4,a)') &
           trim(inpglob)//'/SST/ts_Amon_GFDL-ESM2M_rcp45_r1i1p1_', &
           y1 , '01-', y2, '12.nc'
      else if ( ssttyp(4:5) == '85' ) then
        write(inpfile,'(a,i4,a,i4,a)') &
           trim(inpglob)//'/SST/ts_Amon_GFDL-ESM2M_rcp85_r1i1p1_', &
           y1 , '01-', y2, '12.nc'
      end if
      istatus = nf90_open(inpfile,nf90_nowrite,inet1)
      call checkncerr(istatus,__FILE__,__LINE__,'Error opening '//trim(inpfile))
      write (stdout,*) inet1 , trim(inpfile)
      istatus = nf90_inq_varid(inet1,varname(1),ivar2(1))
      call checkncerr(istatus,__FILE__,__LINE__,'Error find var '//varname(1))
      istatus = nf90_inq_varid(inet1,varname(2),ivar2(2))
      call checkncerr(istatus,__FILE__,__LINE__,'Error find var '//varname(2))
      istatus = nf90_inq_dimid(inet1,'time',timid)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find dim time')
      istatus = nf90_inquire_dimension(inet1,timid,len=timlen)
      call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dim time')
      call getmem1d(work1,1,timlen,'mod_gnmnc_sst:work1')
      istart_t(1) = 1
      icount_t(1) = timlen
      istatus = nf90_get_var(inet1,ivar2(1),work1,istart_t,icount_t)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(1))
      istatus = nf90_get_att(inet1,ivar2(1),'units',cunit)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//varname(1)//' units')
      istatus = nf90_get_att(inet1,ivar2(1),'calendar',ccal)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//varname(1)//' calendar')
      fidate1 = timeval2date(work1(1),cunit,ccal)
      it = imondiff(idate,fidate1) + 1
      istart(3) = it
      istatus = nf90_get_var(inet1,ivar2(2),work3,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(2))
      next = timeval2date(work1(it),cunit,ccal)
    else
      istart(3) = it-1
      istatus = nf90_get_var(inet1,ivar2(2),work2,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(2))
      istart(3) = it
      istatus = nf90_get_var(inet1,ivar2(2),work3,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(2))
      prev = timeval2date(work1(it-1),cunit,ccal)
      next = timeval2date(work1(it),cunit,ccal)
    end if
  else if ( ssttyp(1:3) == 'CN_' ) then
    if ( it > timlen ) then
      istart(3) = it-1
      istatus = nf90_get_var(inet1,ivar2(2),work2,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(2))
      prev = timeval2date(work1(it-1),cunit,ccal)
      ! Try switching to next file
      istatus = nf90_close(inet1)
      call checkncerr(istatus,__FILE__,__LINE__,'Error Close file')
      call split_idate(idate, year, month, day, hour)  
      if ( year < 2006 ) then
        y1 = (year)/10*10
        if ( y1 == 2000 ) then
          y2 = 2005
        else
          y2 = y1 + 9
        end if
      else
        y1 = (year-6)/10*10+6
        y2 = y1 + 9
      end if
      if ( year /= 2006 .and. month /= 1 .and. day /= 1 .and. hour /= 0 ) then
        if ( year == y1 .and. month == 1 .and. day == 1 .and. hour == 0 ) then
          y1 = y1 - 10
          y2 = y1 + 9
        end if
      end if
      if ( ssttyp(4:5) == 'RF' ) then
        write(inpfile,'(a,i4,a,i4,a)') &
           trim(inpglob)//'/SST/tos_Omon_CNRM-CM5_historical_r1i1p1_', &
           y1 , '01-', y2, '12.nc'
      else if ( ssttyp(4:5) == '26' ) then
        write(inpfile,'(a,i4,a,i4,a)') &
           trim(inpglob)//'/SST/tos_Omon_CNRM-CM5_rcp26_r1i1p1_', &
           y1 , '01-', y2, '12.nc'
      else if ( ssttyp(4:5) == '45' ) then
        write(inpfile,'(a,i4,a,i4,a)') &
           trim(inpglob)//'/SST/tos_Omon_CNRM-CM5_rcp45_r1i1p1_', &
           y1 , '01-', y2, '12.nc'
      else if ( ssttyp(4:5) == '85' ) then
        write(inpfile,'(a,i4,a,i4,a)') &
           trim(inpglob)//'/SST/tos_Omon_CNRM-CM5_rcp85_r1i1p1_', &
           y1 , '01-', y2, '12.nc'
      end if
      istatus = nf90_open(inpfile,nf90_nowrite,inet1)
      call checkncerr(istatus,__FILE__,__LINE__,'Error opening '//trim(inpfile))
      write (stdout,*) inet1 , trim(inpfile)
      istatus = nf90_inq_varid(inet1,varname(1),ivar2(1))
      call checkncerr(istatus,__FILE__,__LINE__,'Error find var '//varname(1))
      istatus = nf90_inq_varid(inet1,varname(2),ivar2(2))
      call checkncerr(istatus,__FILE__,__LINE__,'Error find var '//varname(2))
      istatus = nf90_inq_dimid(inet1,'time',timid)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find dim time')
      istatus = nf90_inquire_dimension(inet1,timid,len=timlen)
      call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dim time')
      call getmem1d(work1,1,timlen,'mod_gnmnc_sst:work1')
      istart_t(1) = 1
      icount_t(1) = timlen
      istatus = nf90_get_var(inet1,ivar2(1),work1,istart_t,icount_t)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(1))
      istatus = nf90_get_att(inet1,ivar2(1),'units',cunit)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//varname(1)//' units')
      istatus = nf90_get_att(inet1,ivar2(1),'calendar',ccal)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//varname(1)//' calendar')
      fidate1 = timeval2date(work1(1),cunit,ccal)
      it = imondiff(idate,fidate1) + 1
      istart(3) = it
      istatus = nf90_get_var(inet1,ivar2(2),work3,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(2))
      next = timeval2date(work1(it),cunit,ccal)
    else
      istart(3) = it-1
      istatus = nf90_get_var(inet1,ivar2(2),work2,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(2))
      istart(3) = it
      istatus = nf90_get_var(inet1,ivar2(2),work3,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(2))
      prev = timeval2date(work1(it-1),cunit,ccal)
      next = timeval2date(work1(it),cunit,ccal)
    end if
  else
    istart(3) = it-1
    istatus = nf90_get_var(inet1,ivar2(2),work2,istart,icount)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(2))
    istart(3) = it
    istatus = nf90_get_var(inet1,ivar2(2),work3,istart,icount)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(2))
    prev = timeval2date(work1(it-1),cunit,ccal)
    next = timeval2date(work1(it),cunit,ccal)
  end if
  tdiff1 = next-idate
  tdiff2 = next-prev
  wt1 = real(tohours(tdiff1)/tohours(tdiff2))
  wt2 = 1.0 - wt1
  do j = 1 , jlat
    do i = 1 , ilon
      if (work2(i,j) < 0.9E+20 .and. work3(i,j) < 0.9E+20 ) then
        sst(i,j) = work2(i,j)*wt2+work3(i,j)*wt1
      else
        sst(i,j) = -9999.0
      end if
    end do
  end do

  end subroutine gnmnc_sst
!
end module mod_sst_gnmnc

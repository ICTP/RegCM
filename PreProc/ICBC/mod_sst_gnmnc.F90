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

  use netcdf
  use mod_realkinds
  use mod_stdio
  use mod_dynparam
  use mod_memutil
  use mod_sst_grid
  use mod_interp
  use mod_message
  use mod_nchelper

  private

  integer :: ilon , jlat 
!
  integer :: inet1
  integer , dimension(2) :: ivar2
  integer :: timlen
  integer :: timid
  integer :: istatus
  integer , dimension(3) :: istart , icount
  real(dp) , pointer ::  work1(:)
  real(sp) , pointer , dimension (:, :) :: work2 , work3
  real(sp) , pointer , dimension(:,:) :: sst
  type(rcm_time_and_date) , save :: fidate1
  character(64) :: cunit , ccal
  character(256) :: inpfile
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
  real(sp) , pointer , dimension(:) :: glat
  real(sp) , pointer , dimension(:) :: glon
  real(sp) , pointer , dimension(:,:) :: glat2
  real(sp) , pointer , dimension(:,:) :: glon2
  type(rcm_time_and_date) :: idate , idatef , idateo
  integer :: i , j , k , nsteps , latid , lonid
  integer :: year , month , day , hour , y1 , y2
  real(sp) :: ufac
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
    if (year*1000000+month*10000+day*100+hour > 2005110100) then
!     use modified file (r45 first time step is added to hist last time step)
      inpfile = trim(inpglob)// &
           '/SST/ts_Amon_HadGEM2-ES_historical_r1i1p1_193412-200512.nc'
    else
      inpfile = trim(inpglob)// &
           '/SST/ts_Amon_HadGEM2-ES_historical_r1i1p1_193412-200511.nc'
    end if
!    inpfile = trim(inpglob)// &
!          '/SST/tos_Omon_HadGEM2-ES_historical_r1i1p1_195912-200512.nc'
    varname(2) = 'ts'
!    varname(2) = 'tos'
  else if ( ssttyp == "HA_26" ) then
    inpfile = trim(inpglob)//'/SST/ts_Amon_HadGEM2-ES_rcp26_r1i1p1_200512-209911.nc'
    varname(2) = 'ts'
  else if ( ssttyp == "HA_45" ) then
    inpfile = trim(inpglob)//'/SST/ts_Amon_HadGEM2-ES_rcp45_r1i1p1_200511-209911.nc'
    varname(2) = 'ts'
  else if ( ssttyp == "HA_85" ) then
    inpfile = trim(inpglob)//'/SST/ts_Amon_HadGEM2-ES_rcp85_r1i1p1_200512-209911.nc'
    varname(2) = 'ts'
  else if ( ssttyp == "EC_RF" ) then
    inpfile = trim(inpglob)//'/SST/EC-EARTH/RF/ich1_sst_1950-2009.nc'
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
  else
    call die('gnmnc_sst','Unknown ssttyp: '//ssttyp,1)
  end if

  istatus = nf90_open(inpfile,nf90_nowrite,inet1)
  call checkncerr(istatus,__FILE__,__LINE__,'Error opening '//trim(inpfile))
  write (stdout,*) inet1 , trim(inpfile)

! GET DIMENSION IDs
  if ( ssttyp(1:3) /= 'IP_' ) then
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

  if ( ssttyp(1:3) /= 'IP_' ) then
    call getmem1d(glat,1,jlat,'mod_gnmnc_sst:glat')
    call getmem1d(glon,1,ilon,'mod_gnmnc_sst:glon')
!   GET LATITUDE AND LONGITUDE
    istatus = nf90_get_var(inet1,latid,glat)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var lat')
    istatus = nf90_get_var(inet1,lonid,glon)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var lon')
  else  
    call getmem2d(glat2,1,ilon,1,jlat,'mod_gnmnc_sst:glat2')
    call getmem2d(glon2,1,ilon,1,jlat,'mod_gnmnc_sst:glon2')
!   GET LATITUDE AND LONGITUDE
    istatus = nf90_get_var(inet1,latid,glat2)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var lat')
    istatus = nf90_get_var(inet1,lonid,glon2)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var lon')
!    call getmem1d(glat,1,jlat,'mod_gnmnc_sst:glat')
!    call getmem1d(glon,1,ilon,'mod_gnmnc_sst:glon')
!    glat(:) = glat2(1,:)
!    glon(:) = glon2(:,1)
!    call relmem2d(glat2)
!    call relmem2d(glon2)
    where (glon2 >= 180.0)
      glon2 = glon2-360.0
    end where
  end if
  call getmem2d(work2,1,ilon,1,jlat,'mod_gnmnc_sst:work2')
  call getmem2d(work3,1,ilon,1,jlat,'mod_gnmnc_sst:work3')
  call getmem2d(sst,1,ilon,1,jlat,'mod_gnmnc_sst:sst')
  
! GET TIME VALUES
  istart(1) = 1
  icount(1) = timlen
  istatus = nf90_get_var(inet1,ivar2(1),work1,istart,icount)
  call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(1))
! CHECK FOR THE REQUIRED RECORD IN DATA FILE  
  istatus = nf90_get_att(inet1,ivar2(1),'units',cunit)
  call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(1)//' units')
  istatus = nf90_get_att(inet1,ivar2(1),'calendar',ccal)
  call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(1)//' calendar')
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
    if ( ssttyp(1:3) == 'IP_' ) then
      call distwgtcr(sstmm,sst,xlon,xlat,glon2,glat2,jx,iy,ilon,jlat)
    else
      call bilinx(sst,sstmm,xlon,xlat,glon,glat,ilon,jlat,iy,jx,1)
    end if

    do j = 2 , jx-1
      do i = 2 , iy-1
        sstmm(i,j) = sstmm(i,j) + ufac
      end do
    end do

    call writerec(idate,.false.)
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
  real(sp) :: wt1 , wt2
  type(rcm_time_and_date) :: prev , next

  integer :: it , i , j
  integer :: year , month , day , hour , y1 , y2
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
      istart(1) = 1
      icount(1) = timlen
      istatus = nf90_get_var(inet1,ivar2(1),work1,istart,icount)
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
      if (work2(i,j) < 0.9E+20) then
        sst(i,j) = work2(i,j)*wt2+work3(i,j)*wt1
      else
        sst(i,j) = -9999.0
      end if
    end do
  end do

  end subroutine gnmnc_sst
!
end module mod_sst_gnmnc

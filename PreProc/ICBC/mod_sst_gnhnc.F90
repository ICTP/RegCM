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

module mod_sst_gnhnc

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
  real(sp) , pointer , dimension (:, :) :: work2
  real(sp) , pointer , dimension(:,:) :: sst
  type(rcm_time_and_date) , save :: fidate1
  character(64) :: cunit , ccal
  character(256) :: inpfile
  character(len=8), dimension(2) :: varname
!
  data varname/'time', 'TOBESET'/
!
  public :: sst_gnhnc

  contains

  subroutine sst_gnhnc
!
!*******************************************************************************
!
! This is a package of subroutines to read SST data on regular latlon grid in
! NETCDF format and interpolate on RCM grid.
! Need to configure datapath and variable names.
!
!******************************************************************************

  implicit none
!
  real(sp) , pointer , dimension(:) :: glat
  real(sp) , pointer , dimension(:) :: glon
  real(sp) , pointer , dimension(:,:) :: glat2
  real(sp) , pointer , dimension(:,:) :: glon2
  type(rcm_time_and_date) :: idate , idatef , idateo
  type(rcm_time_interval) :: tdif
  integer :: i , j , k , nsteps , latid , lonid
  integer :: year , month , day , hour , y1 , y2 , m1 , m2
!
  call split_idate(globidate1, year, month, day, hour)  
  y1 = year
  m1 = month
  y2 = year
  m2 = month+1
  if ( m2 > 12 ) then
    m2 = 1
    y2 = y2 + 1
  end if
!
  if ( ssttyp == "MP_RF" ) then
    write(inpfile,'(a,i9,a,i9,a)') &
      trim(inpglob)//'/MPI-ESM-MR/SST/'// &
            'tos_6hrLev_MPI-ESM-MR_historical_r1i1p1_', &
      y1*100000+m1*1000+10,'000-',y2*100000+m2*1000+10,'000.nc'
    varname(2) = 'tos'
  else if ( ssttyp == "MP_45" ) then
    write(inpfile,'(a,i9,a,i9,a)') &
      trim(inpglob)//'/MPI-ESM-MR/SST/tos_6hrLev_MPI-ESM-MR_rcp45_r1i1p1_', &
      y1*100000+m1*1000+10,'000-',y2*100000+m2*1000+10,'000.nc'
    varname(2) = 'tos'
  else if ( ssttyp == "MP_85" ) then
    write(inpfile,'(a,i9,a,i9,a)') &
      trim(inpglob)//'/MPI-ESM-MR/SST/tos_6hrLev_MPI-ESM-MR_rcp85_r1i1p1_', &
      y1*100000+m1*1000+10,'000-',y2*100000+m2*1000+10,'000.nc'
    varname(2) = 'tos'
  else
    call die('gnhnc_sst','Unknown ssttyp: '//ssttyp,1)
  end if

  istatus = nf90_open(inpfile,nf90_nowrite,inet1)
  call checkncerr(istatus,__FILE__,__LINE__,'Error opening '//trim(inpfile))
  write (stdout,*) inet1 , trim(inpfile)

! GET DIMENSION IDs
  istatus = nf90_inq_dimid(inet1,'lat',latid)
  call checkncerr(istatus,__FILE__,__LINE__,'Error find dim lat')
  istatus = nf90_inq_dimid(inet1,'lon',lonid)
  call checkncerr(istatus,__FILE__,__LINE__,'Error find dim lon')
  istatus = nf90_inquire_dimension(inet1,latid,len=jlat)
  call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dim lat')
  istatus = nf90_inquire_dimension(inet1,lonid,len=ilon)
  call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dim lon')

! GET DIMENSION LENGTHS
  istatus = nf90_inq_dimid(inet1,'time',timid)
  call checkncerr(istatus,__FILE__,__LINE__,'Error find dim lon')
  istatus = nf90_inquire_dimension(inet1,timid,len=timlen)
  call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dim time')

  call getmem1d(work1,1,timlen,'mod_gnhnc_sst:work1')

! GET VARIABLE IDs
  istatus = nf90_inq_varid(inet1,'lat',latid)
  call checkncerr(istatus,__FILE__,__LINE__,'Error find var lat')
  istatus = nf90_inq_varid(inet1,'lon',lonid)
  call checkncerr(istatus,__FILE__,__LINE__,'Error find var lon')
  istatus = nf90_inq_varid(inet1,varname(1),ivar2(1))
  call checkncerr(istatus,__FILE__,__LINE__,'Error find var '//varname(1))
  istatus = nf90_inq_varid(inet1,varname(2),ivar2(2))
  call checkncerr(istatus,__FILE__,__LINE__,'Error find var '//varname(2))

  call getmem1d(glat,1,jlat,'mod_gnhnc_sst:glat')
  call getmem1d(glon,1,ilon,'mod_gnhnc_sst:glon')
  call getmem2d(glat2,1,ilon,1,jlat,'mod_gnhnc_sst:glat2')
  call getmem2d(glon2,1,ilon,1,jlat,'mod_gnhnc_sst:glon2')
! GET LATITUDE AND LONGITUDE
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
  call getmem2d(work2,1,ilon,1,jlat,'mod_gnhnc_sst:work2')
  call getmem2d(sst,1,ilon,1,jlat,'mod_gnhnc_sst:sst')
  
! GET TIME VALUES
  istart(1) = 1
  icount(1) = timlen
  istatus = nf90_get_var(inet1,ivar2(1),work1,istart,icount)
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
  idatef = globidate2
  tdif = idatef-idateo
  nsteps = tohours(tdif)/6 + 1
  write (stdout,*) 'NSTEPS = ', nsteps
 
  call open_sstfile(idateo)
 
  idate = idateo
  tdif = 6*3600
  do k = 1 , nsteps
    call gnhnc_sst(idate)
    call distwgtcr(sstmm,sst,xlon,xlat,glon2,glat2,jx,iy,ilon,jlat)
    call writerec(idate,.false.)
    write (stdout,*) 'WRITEN OUT SST DATA : ' , tochar(idate)
    idate = idate + tdif
  end do

  end subroutine sst_gnhnc
!
!-----------------------------------------------------------------------
!
!     Subroutine to read required records from SST data file
!
  subroutine gnhnc_sst(idate)
  implicit none

  type(rcm_time_and_date) , intent (in) :: idate
  integer :: it , i , j
  integer :: year , month , day , hour , y1 , y2 , m1 , m2
  type(rcm_time_interval) :: tdif

  istart(3) = 1

  tdif = idate-fidate1
  it = tohours(tdif)/6 + 1

  if ( it > timlen ) then
    ! Try switching to next file
    istatus = nf90_close(inet1)
    call checkncerr(istatus,__FILE__,__LINE__,'Error Close file')
    call split_idate(idate, year, month, day, hour)  
    y1 = year
    m1 = month
    y2 = year
    m2 = month+1
    if ( m2 > 12 ) then
      m2 = 1
      y2 = y2 + 1
    end if
    if ( ssttyp == "MP_RF" ) then
      write(inpfile,'(a,i9,a,i9,a)') &
        trim(inpglob)//'/MPI-ESM-MR/SST/'// &
          'tos_6hrLev_MPI-ESM-MR_historical_r1i1p1_', &
        y1*100000+m1*1000+10,'000-',y2*100000+m2*1000+10,'000.nc'
    else if ( ssttyp == "MP_45" ) then
      write(inpfile,'(a,i9,a,i9,a)') &
        trim(inpglob)//'/MPI-ESM-MR/SST/tos_6hrLev_MPI-ESM-MR_rcp45_r1i1p1_', &
        y1*100000+m1*1000+10,'000-',y2*100000+m2*1000+10,'000.nc'
    else if ( ssttyp == "MP_85" ) then
      write(inpfile,'(a,i9,a,i9,a)') &
        trim(inpglob)//'/MPI-ESM-MR/SST/tos_6hrLev_MPI-ESM-MR_rcp85_r1i1p1_', &
        y1*100000+m1*1000+10,'000-',y2*100000+m2*1000+10,'000.nc'
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
    call getmem1d(work1,1,timlen,'mod_gnhnc_sst:work1')
    istart(1) = 1
    icount(1) = timlen
    istatus = nf90_get_var(inet1,ivar2(1),work1,istart(1:1),icount(1:1))
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(1))
    istatus = nf90_get_att(inet1,ivar2(1),'units',cunit)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var '//varname(1)//' units')
    istatus = nf90_get_att(inet1,ivar2(1),'calendar',ccal)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var '//varname(1)//' calendar')
    fidate1 = timeval2date(work1(1),cunit,ccal)
    tdif = idate-fidate1
    it = tohours(tdif)/6 + 1
  end if
  icount(1) = ilon
  icount(2) = jlat
  icount(3) = 1
  istart(1) = 1
  istart(2) = 1
  istart(3) = it
  istatus = nf90_get_var(inet1,ivar2(2),work2,istart,icount)
  call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(2))
  do j = 1 , jlat
    do i = 1 , ilon
      if (work2(i,j) < 0.9E+20 ) then
        sst(i,j) = work2(i,j)
      else
        sst(i,j) = -9999.0
      end if
    end do
  end do
  end subroutine gnhnc_sst
!
end module mod_sst_gnhnc

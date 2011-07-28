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

  use m_realkinds
  use m_die
  use m_stdio
  use m_mall
  use m_zeit
  use mod_dynparam
  use mod_memutil
  use mod_sst_grid
  use mod_interp
  use netcdf

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
  type(rcm_time_and_date) :: fidate1
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
  type(rcm_time_and_date) :: idate , idatef , idateo
  integer :: i , j , k , ludom , lumax , iv , nsteps , latid , lonid
  integer , dimension(20) :: lund
  real(sp) :: ufac
!
  call zeit_ci('sst_gnmnc')

  ufac = 0.0
  if ( ssttyp == "CAM2N" ) then
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
    inpfile = trim(inpglob)// &
         '/SST/ts_Amon_HadGEM2-ES_historical_r1i1p1_193412-200511.nc'
    varname(2) = 'ts'
  else if ( ssttyp == "HA_26" ) then
    inpfile = trim(inpglob)//'/SST/ts_Amon_HadGEM2-ES_rcp26_r1i1p1_200512-209911.nc'
    varname(2) = 'ts'
  else if ( ssttyp == "HA_45" ) then
    inpfile = trim(inpglob)//'/SST/ts_Amon_HadGEM2-ES_rcp45_r1i1p1_200512-209911.nc'
    varname(2) = 'ts'
  else if ( ssttyp == "HA_85" ) then
    inpfile = trim(inpglob)//'/SST/ts_Amon_HadGEM2-ES_rcp85_r1i1p1_200512-209911.nc'
    varname(2) = 'ts'
  else
    call die('gnmnc_sst','Unknown ssttyp: '//ssttyp,1)
  end if

  istatus = nf90_open(inpfile,nf90_nowrite,inet1)
  if ( istatus /= nf90_noerr ) then
    call die('gnmnc_sst','Error opening '//trim(inpfile),1, &
               nf90_strerror(istatus),istatus)
  end if
  write (stdout,*) inet1 , trim(inpfile)
! GET DIMENSION IDs
  istatus = nf90_inq_dimid(inet1,'lat',latid)
  if ( istatus /= nf90_noerr ) then
    call die('gnmnc_sst','Error '//trim(inpfile)//' dim lat',1,  &
             nf90_strerror(istatus),istatus)
  end if
  istatus = nf90_inq_dimid(inet1,'lon',lonid)
  if ( istatus /= nf90_noerr ) then
    call die('gnmnc_sst','Error '//trim(inpfile)//' dim lon',1,  &
             nf90_strerror(istatus),istatus)
  end if
  istatus = nf90_inq_dimid(inet1,'time',timid)
  if ( istatus /= nf90_noerr ) then
    call die('gnmnc_sst','Error '//trim(inpfile)//' dim time',1, &
             nf90_strerror(istatus),istatus)
  end if

! GET DIMENSION LENGTHS
  istatus = nf90_inquire_dimension(inet1,latid,len=jlat)
  if ( istatus /= nf90_noerr ) then
    call die('gnmnc_sst','Error '//trim(inpfile)//' dim lat',1,  &
             nf90_strerror(istatus),istatus)
  end if
  istatus = nf90_inquire_dimension(inet1,lonid,len=ilon)
  if ( istatus /= nf90_noerr ) then
    call die('gnmnc_sst','Error '//trim(inpfile)//' dim lon',1,  &
             nf90_strerror(istatus),istatus)
  end if
  istatus = nf90_inquire_dimension(inet1,timid,len=timlen)
  if ( istatus /= nf90_noerr ) then
    call die('gnmnc_sst','Error '//trim(inpfile)//' dim time',1, &
             nf90_strerror(istatus),istatus)
  end if

  call getmem1d(work1,1,timlen,'mod_gnmnc_sst:work1')
  call getmem1d(glat,1,jlat,'mod_gnmnc_sst:glat')
  call getmem1d(glon,1,ilon,'mod_gnmnc_sst:glon')
  call getmem2d(work2,1,ilon,1,jlat,'mod_gnmnc_sst:work2')
  call getmem2d(work3,1,ilon,1,jlat,'mod_gnmnc_sst:work3')
  call getmem2d(sst,1,ilon,1,jlat,'mod_gnmnc_sst:sst')
  
! GET VARIABLE IDs
  istatus = nf90_inq_varid(inet1,'lat',latid)
  if ( istatus /= nf90_noerr ) then
    call die('gnmnc_sst','Error '//trim(inpfile)//':'// &
             'lat',1,nf90_strerror(istatus),istatus)
  end if
  istatus = nf90_inq_varid(inet1,'lon',lonid)
  if ( istatus /= nf90_noerr ) then
    call die('gnmnc_sst','Error '//trim(inpfile)//':'// &
             'lon',1,nf90_strerror(istatus),istatus)
  end if
  istatus = nf90_inq_varid(inet1,varname(1),ivar2(1))
  if ( istatus /= nf90_noerr ) then
    call die('gnmnc_sst','Error '//trim(inpfile)//':'// &
             varname(1),1,nf90_strerror(istatus),istatus)
  end if
  istatus = nf90_inq_varid(inet1,varname(2),ivar2(2))
  if ( istatus /= nf90_noerr ) then
    call die('gnmnc_sst','Error '//trim(inpfile)//':'// &
             varname(2),1,nf90_strerror(istatus),istatus)
  end if
! GET LATITUDE AND LONGITUDE
  istatus = nf90_get_var(inet1,latid,glat)
  if ( istatus /= nf90_noerr ) then
    call die('gnmnc_sst','Error '//trim(inpfile)//':'// &
             'lat read',1,nf90_strerror(istatus),istatus)
  end if
  istatus = nf90_get_var(inet1,lonid,glon)
  if ( istatus /= nf90_noerr ) then
    call die('gnmnc_sst','Error '//trim(inpfile)//':'// &
             'lon read',1,nf90_strerror(istatus),istatus)
  end if
! GET TIME VALUES
  istart(1) = 1
  icount(1) = timlen
  istatus = nf90_get_var(inet1,ivar2(1),work1,istart,icount)
  if ( istatus /= nf90_noerr ) then
    call die('gnmnc_sst','Error '//trim(inpfile)//':'// &
             varname(1)//' read',1,nf90_strerror(istatus),istatus)
  end if
! CHECK FOR THE REQUIRED RECORD IN DATA FILE  
  istatus = nf90_get_att(inet1,ivar2(1),'units',cunit)
  if ( istatus /= nf90_noerr ) then
    call die('gnmnc_sst','Error '//trim(inpfile)//':'// &
            varname(1)//':units',1,nf90_strerror(istatus),istatus)
  end if
  istatus = nf90_get_att(inet1,ivar2(1),'calendar',ccal)
  if ( istatus /= nf90_noerr ) then
    call die('gnmnc_sst','Error '//trim(inpfile)//':'// &
            varname(1)//':calendar',1,nf90_strerror(istatus),istatus)
  end if
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
    call bilinx(sst,sstmm,xlon,xlat,glon,glat,ilon,jlat,iy,jx,1)

    do j = 1 , jx
      do i = 1 , iy
        if ( sstmm(i,j) < -5000 .and.      &
            (lu(i,j) > 13.5 .and. lu(i,j) < 15.5) ) then
          do iv = 1 , 20
            lund(iv) = 0
          end do
          lund(nint(lu(i-1,j-1))) = lund(nint(lu(i-1,j-1))) + 2
          lund(nint(lu(i-1,j))) = lund(nint(lu(i-1,j))) + 3
          lund(nint(lu(i-1,j+1))) = lund(nint(lu(i-1,j+1))) + 2
          lund(nint(lu(i,j-1))) = lund(nint(lu(i,j-1))) + 3
          lund(nint(lu(i,j+1))) = lund(nint(lu(i,j+1))) + 3
          lund(nint(lu(i+1,j-1))) = lund(nint(lu(i+1,j-1))) + 2
          lund(nint(lu(i+1,j))) = lund(nint(lu(i+1,j))) + 3
          lund(nint(lu(i+1,j+1))) = lund(nint(lu(i+1,j+1))) + 2
          ludom = 18
          lumax = 0
          do iv = 1 , 20
            if ( iv <= 13 .or. iv >= 16 ) then
              if ( lund(iv) > lumax ) then
                ludom = k
                lumax = lund(iv)
              end if
            end if
          end do
          lu(i,j) = float(ludom)
        end if
        sstmm(i,j) = sstmm(i,j) + ufac
      end do
    end do

    call writerec(idate,.false.)
    write (stdout,*) 'WRITEN OUT SST DATA : ' , idate%tostring()

    idate = nextmon(idate)

  end do
 
  call zeit_co('sst_gnmnc')

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
  type(rcm_time_interval) :: tdiff1 , tdiff2

  call zeit_ci('read_gnmnc')

  it = imondiff(idate,fidate1) + 1
  icount(1) = ilon
  icount(2) = jlat
  icount(3) = 1
  istart(1) = 1
  istart(2) = 1
  istart(3) = 1
  istart(3) = it-1
  istatus = nf90_get_var(inet1,ivar2(2),work2,istart,icount)
  if ( istatus /= nf90_noerr ) then
    call die('gnmnc_sst','Error '//trim(inpfile)//':'// &
              varname(2)//' read',1,nf90_strerror(istatus),istatus)
  end if
  istart(3) = it
  istatus = nf90_get_var(inet1,ivar2(2),work3,istart,icount)
  if ( istatus /= nf90_noerr ) then
    call die('gnmnc_sst','Error '//trim(inpfile)//':'// &
              varname(2)//' read',1,nf90_strerror(istatus),istatus)
  end if

  prev = timeval2date(work1(it-1),cunit,ccal)
  next = timeval2date(work1(it),cunit,ccal)
  tdiff1 = next-idate
  tdiff2 = next-prev
  wt1 = real(tdiff1%hours( )/tdiff2%hours())
  wt2 = 1.0 - wt1
  do j = 1 , jlat
    do i = 1 , ilon
      sst(i,j) = work2(i,j)*wt2+work3(i,j)*wt1
    end do
  end do

  call zeit_co('read_gnmnc')
  end subroutine gnmnc_sst
!
end module mod_sst_gnmnc

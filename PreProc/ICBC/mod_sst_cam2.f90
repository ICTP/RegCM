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

module mod_sst_cam2

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
  type(rcm_time_and_date) :: cssidate1
  character(64) :: cunit , ccal
  character(256) :: inpfile
  character(len=8), dimension(2) :: varname
!
  data varname/'time','SST_cpl'/

  public :: sst_cam2

  contains

  subroutine sst_cam2
!
!*******************************************************************************
!
! This is a package of subroutines to read CCSM SST 1x1 degree data in
! NETCDF format and interpolate SST at RCM grid
! Written By Moetasim Ashfaq Dec-2005 @ PURDUE.EDU
!
!******************************************************************************
!******************************************************************************
!       DATA PREPARATION
!       Dataset required to use this code can be preapred using NCO utilitiles
!       such as NCKS, NCRCAT etc.
!       We need top level of SST for SSTs which can be extarcted as following:
!       ncks -v time,SST -d z_t,0 input.nc output.nc
!       Files can be further concatenated using 'ncrcat'
!       Finally, the POP grid can be converted into lat/lon grid at 1x1 degree
!       resolution  using PopLatLon function in NCL
!******************************************************************************
!     NAMING CONVENTION (Global Data File)
!       cam2_mn.sst.nc  for SST
!     PATH /DATA/SST/
!
!******************************************************************************

  implicit none
!
  real(sp) , pointer , dimension(:) :: glat
  real(sp) , pointer , dimension(:) :: glon
  type(rcm_time_and_date) :: idate , idatef , idateo
  integer :: i , j , k , ludom , lumax , iv , nsteps , latid , lonid
  integer , dimension(20) :: lund
!
  call zeit_ci('sst_cam2')

  inpfile = trim(inpglob)//'/CAM2/sst_HadOIBl_bc_0.9x1.25_1870_2008_c091020.nc'

  istatus = nf90_open(inpfile,nf90_nowrite,inet1)
  if ( istatus /= nf90_noerr ) then
    call die('cam2_sst','Error opening '//trim(inpfile),1, &
               nf90_strerror(istatus),istatus)
  end if
  write (stdout,*) inet1 , trim(inpfile)
! GET DIMENSION IDs
  istatus = nf90_inq_dimid(inet1,'lat',latid)
  if ( istatus /= nf90_noerr ) then
    call die('cam2_sst','Error '//trim(inpfile)//' dim lat',1,  &
             nf90_strerror(istatus),istatus)
  end if
  istatus = nf90_inq_dimid(inet1,'lon',lonid)
  if ( istatus /= nf90_noerr ) then
    call die('cam2_sst','Error '//trim(inpfile)//' dim lon',1,  &
             nf90_strerror(istatus),istatus)
  end if
  istatus = nf90_inq_dimid(inet1,'time',timid)
  if ( istatus /= nf90_noerr ) then
    call die('cam2_sst','Error '//trim(inpfile)//' dim time',1, &
             nf90_strerror(istatus),istatus)
  end if

! GET DIMENSION LENGTHS
  istatus = nf90_inquire_dimension(inet1,latid,len=jlat)
  if ( istatus /= nf90_noerr ) then
    call die('cam2_sst','Error '//trim(inpfile)//' dim lat',1,  &
             nf90_strerror(istatus),istatus)
  end if
  istatus = nf90_inquire_dimension(inet1,lonid,len=ilon)
  if ( istatus /= nf90_noerr ) then
    call die('cam2_sst','Error '//trim(inpfile)//' dim lon',1,  &
             nf90_strerror(istatus),istatus)
  end if
  istatus = nf90_inquire_dimension(inet1,timid,len=timlen)
  if ( istatus /= nf90_noerr ) then
    call die('cam2_sst','Error '//trim(inpfile)//' dim time',1, &
             nf90_strerror(istatus),istatus)
  end if

  call getmem1d(work1,1,timlen,'mod_cam2_sst:work1')
  call getmem1d(glat,1,jlat,'mod_cam2_sst:glat')
  call getmem1d(glon,1,ilon,'mod_cam2_sst:glon')
  call getmem2d(work2,1,ilon,1,jlat,'mod_cam2_sst:work2')
  call getmem2d(work3,1,ilon,1,jlat,'mod_cam2_sst:work3')
  call getmem2d(sst,1,ilon,1,jlat,'mod_cam2_sst:sst')
  
! GET VARIABLE IDs
  istatus = nf90_inq_varid(inet1,'lat',latid)
  if ( istatus /= nf90_noerr ) then
    call die('cam2_sst','Error '//trim(inpfile)//':'// &
             'lat',1,nf90_strerror(istatus),istatus)
  end if
  istatus = nf90_inq_varid(inet1,'lon',lonid)
  if ( istatus /= nf90_noerr ) then
    call die('cam2_sst','Error '//trim(inpfile)//':'// &
             'lon',1,nf90_strerror(istatus),istatus)
  end if
  istatus = nf90_inq_varid(inet1,varname(1),ivar2(1))
  if ( istatus /= nf90_noerr ) then
    call die('cam2_sst','Error '//trim(inpfile)//':'// &
             varname(1),1,nf90_strerror(istatus),istatus)
  end if
  istatus = nf90_inq_varid(inet1,varname(2),ivar2(2))
  if ( istatus /= nf90_noerr ) then
    call die('cam2_sst','Error '//trim(inpfile)//':'// &
             varname(2),1,nf90_strerror(istatus),istatus)
  end if
! GET LATITUDE AND LONGITUDE
  istatus = nf90_get_var(inet1,latid,glat)
  if ( istatus /= nf90_noerr ) then
    call die('cam2_sst','Error '//trim(inpfile)//':'// &
             'lat read',1,nf90_strerror(istatus),istatus)
  end if
  istatus = nf90_get_var(inet1,lonid,glon)
  if ( istatus /= nf90_noerr ) then
    call die('cam2_sst','Error '//trim(inpfile)//':'// &
             'lon read',1,nf90_strerror(istatus),istatus)
  end if
! GET TIME VALUES
  istart(1) = 1
  icount(1) = timlen
  istatus = nf90_get_var(inet1,ivar2(1),work1,istart,icount)
  if ( istatus /= nf90_noerr ) then
    call die('cam2_sst','Error '//trim(inpfile)//':'// &
             varname(1)//' read',1,nf90_strerror(istatus),istatus)
  end if
! CHECK FOR THE REQUIRED RECORD IN DATA FILE  
  istatus = nf90_get_att(inet1,ivar2(1),'units',cunit)
  if ( istatus /= nf90_noerr ) then
    call die('cam2_sst','Error '//trim(inpfile)//':'// &
            varname(1)//':units',1,nf90_strerror(istatus),istatus)
  end if
  istatus = nf90_get_att(inet1,ivar2(1),'calendar',ccal)
  if ( istatus /= nf90_noerr ) then
    call die('cam2_sst','Error '//trim(inpfile)//':'// &
            varname(1)//':calendar',1,nf90_strerror(istatus),istatus)
  end if
  cssidate1 = timeval2date(work1(1),cunit,ccal)

  idateo = monfirst(globidate1)
  idatef = monfirst(globidate2)
  if (idatef < globidate2) then
    idatef = nextmon(idatef)
  end if
  nsteps = imondiff(idatef,idateo) + 1
 
  call open_sstfile(idateo)
 
  idate = idateo
  do k = 1 , nsteps

    call cam2_sst(idate)
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
        if ( sstmm(i,j) > -100. ) then
          sstmm(i,j) = sstmm(i,j) + 273.15
        else
          sstmm(i,j) = -9999.
        end if
      end do
    end do

    call writerec(idate,.false.)
    write (stdout,*) 'WRITEN OUT SST DATA : ' , idate%tostring()

    idate = nextmon(idate)

  end do
 
  call zeit_co('sst_cam2')

  end subroutine sst_cam2
!
!-----------------------------------------------------------------------
!
!     Subroutine to read required records from SST data file
!
  subroutine cam2_sst(idate)

  implicit none

  type(rcm_time_and_date) , intent (in) :: idate
  real(sp) :: wt1 , wt2
  type(rcm_time_and_date) :: prev , next

  integer :: it , i , j
  type(rcm_time_interval) :: tdiff1 , tdiff2

  call zeit_ci('read_cam2')

  it = imondiff(idate,cssidate1) + 1
  icount(1) = ilon
  icount(2) = jlat
  icount(3) = 1
  istart(1) = 1
  istart(2) = 1
  istart(3) = 1
  istart(3) = it-1
  istatus = nf90_get_var(inet1,ivar2(2),work2,istart,icount)
  if ( istatus /= nf90_noerr ) then
    call die('cam2_sst','Error '//trim(inpfile)//':'// &
              varname(2)//' read',1,nf90_strerror(istatus),istatus)
  end if
  istart(3) = it
  istatus = nf90_get_var(inet1,ivar2(2),work3,istart,icount)
  if ( istatus /= nf90_noerr ) then
    call die('cam2_sst','Error '//trim(inpfile)//':'// &
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
      if (work2(i,j) > -100.0 .and. work2(i,j) < 100.0) then
        sst(i,j) = work2(i,j)*wt2+work3(i,j)*wt1
      else
        sst(i,j) = -9999.0
      end if
    end do
  end do

  call zeit_co('read_cam2')
  end subroutine cam2_sst
!
end module mod_sst_cam2

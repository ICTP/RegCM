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

module mod_sst_ersst

  use netcdf
  use mod_realkinds
  use mod_stdio
  use mod_dynparam
  use mod_sst_grid
  use mod_interp
  use mod_message
  use mod_memutil
  use mod_nchelper

  private

  public :: sst_ersst

  contains

  subroutine sst_ersst

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Comments on dataset sources and location:                          !
!                                                                    !
! ERAIN    ERAIN_SST is provided by ERA-Interim project              !
!          6 hourly frequncy, 1.5x1.5 degree resolution              !
!          from 1989010100 to 2009053118.                            !
!          'ERSST' for using the sea surface temperature;            !
!          'ERSKT' for using the skin temperature.                   !
!                                                                    !
!          ML = 1 is   0.0; ML = 2 is   1.5; => ML = 240 is 358.5E   !
!          NL = 1 is  90.0; ML = 2 is  88.5; => ML = 121 is -90.     !
!                                                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  implicit none
!
  integer :: ilon , jlat
  integer , parameter :: idtbc = 6
!
  integer :: i , it , j
  integer :: istatus , inet
  integer :: year , month , day , hour , isyear
  integer :: dimi , vari
  real(sp) , pointer , dimension(:) :: lati
  real(sp) , pointer , dimension(:) :: loni
  integer :: ierrec , nsteps
  type(rcm_time_and_date) :: idate , ierastart
  type(rcm_time_interval) :: tdiff , itbc
  real(sp) , pointer , dimension(:,:) :: sst
  character(len=256) :: inpfile
  logical :: lfirst

  data lfirst /.true./
!
  if ( ssttyp /= 'ERSST' .and. ssttyp /= 'ERSKT' ) then
    write (stderr,*) 'PLEASE SET right SSTTYP in regcm.in'
    write (stderr,*) 'Supported types are ERSST ERSKT'
    call die('sst_ersst')
  end if

  itbc = rcm_time_interval(idtbc,uhrs)
  tdiff = globidate2-globidate1
  nsteps = idnint(tohours(tdiff))/idtbc + 1

  write (stdout,*) 'GLOBIDATE1 : ' , tochar(globidate1)
  write (stdout,*) 'GLOBIDATE2 : ' , tochar(globidate2)
  write (stdout,*) 'NSTEPS     : ' , nsteps

  call open_sstfile(globidate1)
 
  ! SET UP LONGITUDES AND LATITUDES FOR SST DATA

  call split_idate(globidate1,year,month,day,hour)
  if ( year > 1978 .and. year < 1989 ) then
    isyear = 1979
    ierastart = 1979010100
    if ( ssttyp == 'ERSST' ) then
      inpfile=trim(inpglob)//'/SST/sstERAIN.1979-1989.nc'
    else if ( ssttyp == 'ERSKT' ) then
      inpfile=trim(inpglob)//'/SST/tskinERAIN.1979-1989.nc'
    end if
  else if ( year > 1988 .and. year < 2009 ) then
    isyear = 1989
    ierastart = 1989010100
    if ( ssttyp == 'ERSST' ) then
      inpfile=trim(inpglob)//'/SST/sstERAIN.1989-2009.nc'
    else if ( ssttyp == 'ERSKT' ) then
      inpfile=trim(inpglob)//'/SST/tskinERAIN.1989-2009.nc'
    end if
  else if ( year > 2008 ) then
    isyear = 2009
    ierastart = 2009010100
    if ( ssttyp == 'ERSST' ) then
      inpfile=trim(inpglob)//'/SST/sstERAIN.2009-present.nc'
    else if ( ssttyp == 'ERSKT' ) then
      inpfile=trim(inpglob)//'/SST/tskinERAIN.2009-present.nc'
    end if
  else
    call die('sst_ersst','The dataset is prepared only for 1979-present',1)
  end if

  write (stdout,*) trim(inpfile)
  istatus = nf90_open(inpfile,nf90_nowrite,inet)
  call checkncerr(istatus,__FILE__,__LINE__,'Cannot open file '//trim(inpfile))
  istatus = nf90_inq_dimid(inet,'longitude',dimi)
  call checkncerr(istatus,__FILE__,__LINE__,'Error find lon dim')
  istatus = nf90_inquire_dimension(inet,dimi,len=ilon)
  call checkncerr(istatus,__FILE__,__LINE__,'Error inquire lon dim')
  istatus = nf90_inq_dimid(inet,'latitude',dimi)
  call checkncerr(istatus,__FILE__,__LINE__,'Error find lat dim')
  istatus = nf90_inquire_dimension(inet,dimi,len=jlat)
  call checkncerr(istatus,__FILE__,__LINE__,'Error inquire lat dim')

  call getmem1d(loni,1,ilon,'sst_ersst:loni')
  call getmem1d(lati,1,jlat,'sst_ersst:lati')
  call getmem2d(sst,1,ilon,1,jlat,'sst_ersst:sst')

  istatus = nf90_inq_varid(inet,'longitude',vari)
  call checkncerr(istatus,__FILE__,__LINE__,'Error find lon var')
  istatus = nf90_get_var(inet,vari,loni)
  call checkncerr(istatus,__FILE__,__LINE__,'Error read lon var')
  istatus = nf90_inq_varid(inet,'latitude',vari)
  call checkncerr(istatus,__FILE__,__LINE__,'Error find lat var')
  istatus = nf90_get_var(inet,vari,lati)
  call checkncerr(istatus,__FILE__,__LINE__,'Error read lat var')

  idate = globidate1
  do it = 1 , nsteps

    call split_idate(idate,year,month,day,hour)
    if ( year > 1988 .and. isyear == 1979 ) then
      ierastart = 1989010100
      if ( ssttyp == 'ERSST' ) then
        inpfile=trim(inpglob)//'/SST/sstERAIN.1989-2009.nc'
      else if ( ssttyp == 'ERSKT' ) then
        inpfile=trim(inpglob)//'/SST/tskinERAIN.1989-2009.nc'
      end if
      istatus = nf90_close(inet)
      call checkncerr(istatus,__FILE__,__LINE__,'Cannot close 1979 file')
      istatus = nf90_open(inpfile,nf90_nowrite,inet)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Cannot open file '//trim(inpfile))
      write (stdout,*) trim(inpfile)
      lfirst = .true.
    else if ( year > 2008 .and. isyear == 2009 ) then
      ierastart = 2009010100
      if ( ssttyp == 'ERSST' ) then
        inpfile=trim(inpglob)//'/SST/sstERAIN.2009-present.nc'
      else if ( ssttyp == 'ERSKT' ) then
        inpfile=trim(inpglob)//'/SST/tskinERAIN.2009-present.nc'
      end if
      istatus = nf90_close(inet)
      call checkncerr(istatus,__FILE__,__LINE__,'Cannot close 1989 file')
      istatus = nf90_open(inpfile,nf90_nowrite,inet)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Cannot open file '//trim(inpfile))
      write (stdout,*) trim(inpfile)
      lfirst = .true.
    end if

    tdiff = idate-ierastart
    ierrec = idnint(tohours(tdiff))/idtbc+1

    if ( ssttyp == 'ERSST' ) then
      call sst_erain(ierrec,ilon,jlat,sst,inet,lfirst,1)
    else if ( ssttyp == 'ERSKT' ) then
      call sst_erain(ierrec,ilon,jlat,sst,inet,lfirst,2)
    end if

    call bilinx(sst,sstmm,xlon,xlat,loni,lati,ilon,jlat,iy,jx,1)
    write(stdout,*) 'XLON,XLAT,SST = ' , xlon(1,1) , xlat(1,1) , sstmm(1,1)
 
    ! WRITE OUT SST DATA ON REGCM GRID
    call writerec(idate,.false.)
    write(stdout,*) 'WRITING OUT MM4 SST DATA:' , tochar(idate)

    idate = idate + itbc
  end do

  end subroutine sst_ersst
!
!-----------------------------------------------------------------------
!
  subroutine sst_erain(it,ilon,jlat,sst,inet,lfirst,itype)
    use netcdf
    implicit none
!
    integer , intent(in) :: it , ilon , jlat , inet , itype
    logical , intent(inout) :: lfirst
    real(sp) , dimension(ilon,jlat) ,intent(out) :: sst
!
    integer :: i , j , n
    character(4) , dimension(2) :: varname
    integer(2) , dimension(ilon,jlat) :: work
    integer :: istatus
!
    integer , save :: ivar
    real(dp) , save :: xadd , xscale , xmiss
    integer , dimension(3) , save :: icount , istart
!
! This is the latitude, longitude dimension of the grid to be read.
! This corresponds to the lat and lon dimension variables in the
! netCDF file.
!
! The data are packed into short integers (INTEGER*2).  The array
! work will be used to hold the packed integers. The array 'sst'
! will contain the unpacked data.
!
! DATA ARRAY AND WORK ARRAY
!
    data varname/'sst','skt'/
!
    if ( lfirst ) then
      istatus = nf90_inq_varid(inet,varname(itype),ivar)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Cannot find variable '//trim(varname(itype)))
      istatus = nf90_get_att(inet,ivar,'scale_factor',xscale)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Cannot get attribute scale_factor')
      istatus = nf90_get_att(inet,ivar,'add_offset',xadd)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Cannot get attribute add_offset')
      istatus = nf90_get_att(inet,ivar,'_FillValue',xmiss)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Cannot get attribute _FillValue')
      istart(1) = 1
      istart(2) = 1
      icount(1) = ilon
      icount(2) = jlat
      lfirst = .false.
    end if
!
    istart(3) = it
    icount(3) = 1
    istatus = nf90_get_var(inet,ivar,work,istart,icount)
    call checkncerr(istatus,__FILE__,__LINE__, &
                   'Cannot read '//trim(varname(itype)))
!
    do j = 1 , jlat
      do i = 1 , ilon
        if (work(i,j) /= xmiss) then
          sst(i,j) = real(dble(work(i,j))*xscale + xadd)
        else
          sst(i,j) = -9999.0
        end if
      end do
    end do
!
  end subroutine sst_erain
!
end module mod_sst_ersst

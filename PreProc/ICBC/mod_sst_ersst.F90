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

  use m_realkinds
  use m_die
  use m_stdio
  use m_zeit
  use netcdf
  use mod_dynparam
  use mod_sst_grid
  use mod_interp
  use mod_message

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
  integer , parameter :: ilon = 240 , jlat = 121
  integer , parameter :: idtbc = 6
!
  integer :: i , it , j
  real(sp) , dimension(jlat) :: lati
  real(sp) , dimension(ilon) :: loni
  integer :: ierrec , nsteps
  type(rcm_time_and_date) :: idate , ierastart
  type(rcm_time_interval) :: tdiff , itbc
  real(sp) , dimension(ilon,jlat) :: sst
  character(256) :: inpfile
!
  call zeit_ci('sst_ersst')
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
 
!     ******    SET UP LONGITUDES AND LATITUDES FOR SST DATA
  do i = 1 , ilon
    loni(i) = float(i-1)*1.5
  end do
  do j = 1 , jlat
    lati(j) = -90. + 1.5*float(j-1)
  end do
 
  idate = globidate1
  ierastart = 1989010100
  do it = 1 , nsteps

    tdiff = idate-ierastart
    ierrec = idnint(tohours(tdiff))/idtbc+1

    if ( ssttyp == 'ERSST' ) then
      inpfile=trim(inpglob)//'/SST/sstERAIN.1989-2009.nc'
      call sst_erain(ierrec,ilon,jlat,sst,inpfile,1)
    else if ( ssttyp == 'ERSKT' ) then
      inpfile=trim(inpglob)//'/SST/tskinERAIN.1989-2009.nc'
      call sst_erain(ierrec,ilon,jlat,sst,inpfile,2)
    end if

    call bilinx(sst,sstmm,xlon,xlat,loni,lati,ilon,jlat,iy,jx,1)
    write(stdout,*) 'XLON,XLAT,SST = ' , xlon(1,1) , xlat(1,1) , sstmm(1,1)
 
!       ******           WRITE OUT SST DATA ON MM4 GRID
    call writerec(idate,.false.)
    write(stdout,*) 'WRITING OUT MM4 SST DATA:' , tochar(idate)

    idate = idate + itbc

  end do

  call zeit_co('sst_ersst')
 
  end subroutine sst_ersst
!
!-----------------------------------------------------------------------
!
  subroutine sst_erain(it,ilon,jlat,sst,pathaddname,itype)
  use netcdf
  implicit none
!
  integer :: it , ilon , jlat , itype
  character(256) :: pathaddname
  intent (in) it , ilon , jlat , pathaddname , itype
  real(sp) , dimension(ilon,jlat) :: sst
  intent (out) :: sst
!
  integer :: i , j , n
  character(4) , dimension(2) :: varname
  integer(2) , dimension(ilon,jlat) :: work
  integer :: istatus
!
  integer , save :: inet , ivar
  real(dp) , save :: xadd , xscale , xmiss
  integer , dimension(10) , save :: icount , istart
  logical , save :: lfirst
!
!     This is the latitude, longitude dimension of the grid to be read.
!     This corresponds to the lat and lon dimension variables in the
!     netCDF file.
!
!     The data are packed into short integers (INTEGER*2).  The array
!     work will be used to hold the packed integers. The array 'sst'
!     will contain the unpacked data.
!
!     DATA ARRAY AND WORK ARRAY
!
  data varname/'sst','skt'/
  data lfirst /.true./
!
  call zeit_ci('read_sst_era')
  if ( lfirst ) then
    istatus = nf90_open(pathaddname,nf90_nowrite,inet)
    call checkncerr(istatus,__FILE__,__LINE__,'Cannot open file '//trim(pathaddname))
    istatus = nf90_inq_varid(inet,varname(itype),ivar)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Cannot find variable '//trim(varname(itype)))
    istatus = nf90_get_att(inet,ivar,'scale_factor',xscale)
    call checkncerr(istatus,__FILE__,__LINE__,'Cannot get attribute scale_factor')
    istatus = nf90_get_att(inet,ivar,'add_offset',xadd)
    call checkncerr(istatus,__FILE__,__LINE__,'Cannot get attribute add_offset')
    istatus = nf90_get_att(inet,ivar,'_FillValue',xmiss)
    call checkncerr(istatus,__FILE__,__LINE__,'Cannot get attribute _FillValue')
    istart(1) = 1
    istart(2) = 1
    icount(1) = 240
    icount(2) = 121
    do n = 4 , 10
      istart(n) = 0
      icount(n) = 0
    end do
    lfirst = .false.
  end if
!
  istart(3) = it
  icount(3) = 1
  istatus = nf90_get_var(inet,ivar,work,istart,icount)
  call checkncerr(istatus,__FILE__,__LINE__,'Cannot read '//trim(varname(itype)))
!
  do j = 1 , jlat
    do i = 1 , ilon
      if (work(i,j) /= xmiss) then
        sst(i,jlat+1-j) = real(dble(work(i,j))*xscale + xadd)
      else
        sst(i,jlat+1-j) = -9999.0
      end if
    end do
  end do
!
  call zeit_co('read_sst_era')
!
  end subroutine sst_erain
!
end module mod_sst_ersst

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

module mod_gn6hnc

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx!
!
! This is a package of subroutines to read 6 hourly data in
! NETCDF format and to prepare Initial and boundary conditions for RegCM.
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx!

  use mod_realkinds
  use mod_stdio
  use mod_dynparam
  use mod_memutil
  use mod_constants
  use mod_grid
  use mod_write
  use mod_interp
  use mod_vertint
  use mod_hgt
  use mod_humid
  use mod_mksst
  use mod_uvrot
  use mod_vectutil
  use mod_message
  use mod_nchelper

  private

  ! Dimension of input read from input files
  integer :: nlon , nlat , klev

  ! Pressure levels to interpolate to if dataset is on model sigma levels.
  integer , parameter :: nipl = 18
  real(sp) , target , dimension(nipl) :: fplev = &
   (/  30.0,   50.0,   70.0,  100.0,  150.0,  200.0,  250.0, &
      300.0,  350.0,  420.0,  500.0,  600.0,  700.0,  780.0, &
      850.0,  920.0,  960.0, 1000.0 /)

  integer :: npl , nrhlev
  real(sp) , pointer , dimension(:) :: pplev
  real(sp) , pointer , dimension(:) :: sigmar

  ! Whole space
  real(sp) , pointer , dimension(:,:,:) :: b2
  real(sp) , pointer , dimension(:,:,:) :: d2
  real(sp) , pointer , dimension(:,:,:) :: b3
  real(sp) , pointer , dimension(:,:,:) :: d3

  real(sp) , pointer :: u3(:,:,:) , v3(:,:,:)
  real(sp) , pointer :: h3(:,:,:) , q3(:,:,:) , t3(:,:,:)
  real(sp) , pointer :: up(:,:,:) , vp(:,:,:)
  real(sp) , pointer :: hp(:,:,:) , qp(:,:,:) , tp(:,:,:)

  ! Input space
  real(sp) :: p0
  real(sp) , pointer , dimension(:,:) :: psvar , zsvar , pmslvar
  real(sp) , pointer , dimension(:) :: ak , bk
  real(sp) , pointer , dimension(:) :: glat , gltemp
  real(sp) , pointer , dimension(:) :: glon
  real(sp) , pointer , dimension(:,:,:) :: hvar , qvar , tvar , &
                                           uvar , vvar , pp3d , &
                                           vwork
  integer :: timlen , pstimlen
  type(rcm_time_and_date) , pointer , dimension(:) :: itimes
  type(rcm_time_and_date) , pointer , dimension(:) :: ipstimes
  real(dp) , pointer , dimension(:) :: xtimes

  ! Shared by netcdf I/O routines
  integer , dimension(4) :: icount , istart
  ! We will need 6 files (is just one for CAM2)
  integer , parameter :: nvars = 6
  integer , parameter :: nfiles = nvars
  integer , dimension(nvars) :: inet
  integer , dimension(nvars) :: ivar

  public :: get_gn6hnc , headgn6hnc

  character(256) :: pathaddname
  type(rcm_time_and_date) , save :: refdate
  type(rcm_time_and_date) , save :: filedate

  data inet /nvars*-1/

  character(32) :: cambase = 'sococa.ts1.r1.cam2.h1.'
  character(64) :: habase1  = '_6hrLev_HadGEM2-ES_historical'
  character(64) :: habase2  = '_6hrLev_HadGEM2-ES_rcp'
  character(64) :: hapbase1 = '_6hrPlev_HadGEM2-ES_historical'
  character(64) :: hapbase2 = '_6hrPlev_HadGEM2-ES_rcp'
  character(64) :: habase3 = '_r1i1p1_'
  character(64) :: cabase1 = '_6hrLev_CanESM2_historical'
  character(64) :: cabase2 = '_6hrLev_CanESM2_rcp'
  character(64) :: cabase3 = '_r1i1p1_'
  character(64) :: ipbase1 = '_6hrLev_IPSL-CM5A-LR_historical'
  character(64) :: ipbase2 = '_6hrLev_IPSL-CM5A-LR_rcp'
  character(64) :: ipbase3 = '_r1i1p1_'
  character(64) :: gfdlbase1 = '_6hrLev_GFDL-ESM2M_historical'
  character(64) :: gfdlbase2 = '_6hrLev_GFDL-ESM2M_rcp'
  character(64) :: gfdlbase3 = '_r1i1p1_'
  character(64) :: cnrmbase1 = '_6hrLev_CNRM-CM5_historical'
  character(64) :: cnrmbase2 = '_6hrLev_CNRM-CM5_rcp'
  character(64) :: cnrmbase3 = '_r1i1p1_'
  character(64) :: csirbase1 = '_6hrLev_CSIRO-Mk3-6-0_historical'
  character(64) :: csirbase2 = '_6hrLev_CSIRO-Mk3-6-0_rcp'
  character(64) :: csirbase3 = '_r1i1p1_'

  character(3) , target , dimension(nvars) :: cam2vars = &
                         (/'T  ' , 'Z3 ' , 'Q  ' , 'U  ' , 'V  ' , 'PS '/)
  character(3) , target , dimension(nvars) :: ccsmvars = &
                         (/'T  ' , 'Z3 ' , 'Q  ' , 'U  ' , 'V  ' , 'PS '/)
  character(3) , target , dimension(nvars) :: havars = &
                         (/'ta ' , 'XXX' , 'hus' , 'ua ' , 'va ' , 'psl'/)
  character(3) , target , dimension(nvars) :: cavars = &
                         (/'ta ' , 'XXX' , 'hus' , 'ua ' , 'va ' , 'ps '/)
  character(3) , target , dimension(nvars) :: ipvars = &
                         (/'ta ' , 'XXX' , 'hus' , 'ua ' , 'va ' , 'ps '/)
  character(3) , target , dimension(nvars) :: gfdlvars = &
                         (/'ta ' , 'XXX' , 'hus' , 'ua ' , 'va ' , 'ps '/)
  character(3) , target , dimension(nvars) :: cnrmvars = &
                         (/'ta ' , 'XXX' , 'hus' , 'ua ' , 'va ' , 'ps '/)
  character(3) , target , dimension(nvars) :: gfsvars = &
                         (/'ta ' , 'hga' , 'rha' , 'ua ' , 'va ' , 'ps '/)
  character(3) , target , dimension(nvars) :: echvars = &
                         (/'t  ' , 'z  ' , 'q  ' , 'u  ' , 'v  ' , '   '/)
  character(3) , target , dimension(nvars) :: csirvars = &
                         (/'ta ' , 'XXX' , 'hus' , 'ua ' , 'va ' , 'ps '/)

  character(4) , dimension(nvars) :: ccsmfname = &
                         (/'air ', 'hgt ', 'shum', 'uwnd', 'vwnd', 'pres'/)

  character(3) , dimension(12) :: mname = &
                         (/'JAN','FEB','MAR','APR','MAY','JUN', &
                           'JUL','AUG','SEP','OCT','NOV','DEC'/)

  character(3) , dimension(:) , pointer :: varname

  contains
!
  subroutine headgn6hnc
    use netcdf
!
    implicit none
!
    integer :: istatus , ivar1 , inet1 , jdim , j , k
    character(256) :: pathaddname
    real(8) :: dp0
!
!
    if ( dattyp == 'CAM4N' ) then
      pathaddname = trim(inpglob)// &
            '/CAM2/USGS-gtopo30_0.9x1.25_remap_c051027.nc'
    else if ( dattyp == 'CCSMN' ) then
      pathaddname = trim(inpglob)//'/CCSM/ccsm_ht.nc'
    else if ( dattyp(1:3) == 'HA_' ) then
      ! Vertical info are not stored in the fixed orography file.
      pathaddname = trim(inpglob)// &
            '/HadGEM2/RF/ta/ta_6hrLev_HadGEM2-ES_historical_'// &
            'r1i1p1_199012010600-199103010000.nc'
    else if ( dattyp(1:3) == 'CA_' ) then
      ! Vertical info are not stored in the fixed orography file.
      ! Read part of the info from first T file.
      pathaddname = trim(inpglob)// &
            '/CanESM2/RF/ta/ta_6hrLev_CanESM2_historical_'// &
            'r1i1p1_195001010000-195012311800.nc'
    else if ( dattyp(1:3) == 'IP_' ) then
      ! Vertical info are not stored in the fixed orography file.
      ! Read part of the info from first T file.
      pathaddname = trim(inpglob)// &
            '/IPSL-CM5A-LR/HIST/ta/ta_6hrLev_IPSL-CM5A-LR_historical_'// &
            'r1i1p1_195001010300-195912312100.nc'
    else if ( dattyp(1:3) == 'GF_' ) then
      ! Vertical info are not stored in the fixed orography file.
      ! Read part of the info from first T file.
      pathaddname = trim(inpglob)// &
            '/GFDL-ESM2M/RF/ta/ta_6hrLev_GFDL-ESM2M_historical_'// &
            'r1i1p1_1951010100-1955123123.nc'
    else if ( dattyp(1:3) == 'CN_' ) then
      ! Vertical info are not stored in the fixed orography file.
      ! Read part of the info from first T file.
      pathaddname = trim(inpglob)// &
            '/CNRM-CM5/RF/ta/ta_6hrLev_CNRM-CM5_historical_'// &
            'r1i1p1_195001010600-195002010000.nc'
    else if ( dattyp(1:3) == 'CS_' ) then
      ! Vertical info are not stored in the fixed orography file.
      ! Read part of the info from first T file.
      pathaddname = trim(inpglob)// &
            '/CSIRO-MK36/RF/ta/ta_6hrLev_CSIRO-Mk3-6-0_historical_'// &
            'r1i1p1_195001010600-195101010000.nc'
    else if ( dattyp == 'GFS11' ) then
      pathaddname = trim(inpglob)//'/GFS11/fixed/fixed_orography.nc'
    else if ( dattyp(1:3) == 'EC_' ) then
      pathaddname = trim(inpglob)//'/EC-EARTH/fixed/ecearth.nc'
    else
      call die('Unknown dattyp in generic 6h NetCDF driver.')
    end if

    istatus = nf90_open(pathaddname,nf90_nowrite,inet1)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error open '//trim(pathaddname))

    istatus = nf90_inq_dimid(inet1,'lon',jdim)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find lon dim')
    istatus = nf90_inquire_dimension(inet1,jdim,len=nlon)
    call checkncerr(istatus,__FILE__,__LINE__,'Error inquire lon dim')
    istatus = nf90_inq_dimid(inet1,'lat',jdim)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find lat dim')
    istatus = nf90_inquire_dimension(inet1,jdim,len=nlat)
    call checkncerr(istatus,__FILE__,__LINE__,'Error inquire lat dim')
    istatus = nf90_inq_dimid(inet1,'lev',jdim)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find lev dim')
    istatus = nf90_inquire_dimension(inet1,jdim,len=klev)
    call checkncerr(istatus,__FILE__,__LINE__,'Error inquire lev dim')

    ! Input layer and pressure interpolated values

    call getmem1d(glat,1,nlat,'mod_gn6hnc:glat')
    call getmem1d(glon,1,nlon,'mod_gn6hnc:glon')
    call getmem2d(zsvar,1,nlon,1,nlat,'mod_gn6hnc:zsvar')
    call getmem2d(psvar,1,nlon,1,nlat,'mod_gn6hnc:psvar')

    if ( dattyp /= 'GFS11' .and. dattyp(1:3) /= 'EC_' ) then
      call getmem3d(qvar,1,nlon,1,nlat,1,klev,'mod_gn6hnc:qvar')
      call getmem3d(tvar,1,nlon,1,nlat,1,klev,'mod_gn6hnc:tvar')
      call getmem3d(hvar,1,nlon,1,nlat,1,klev,'mod_gn6hnc:hvar')
      call getmem3d(uvar,1,nlon,1,nlat,1,klev,'mod_gn6hnc:uvar')
      call getmem3d(vvar,1,nlon,1,nlat,1,klev,'mod_gn6hnc:vvar')
      call getmem3d(pp3d,1,nlon,1,nlat,1,klev,'mod_gn6hnc:pp3d')
      if ( dattyp(1:3) == 'HA_' ) then
        call getmem3d(vwork,1,nlon,1,nlat-1,1,klev,'mod_gn6hnc:vwork')
      end if
      call getmem1d(ak,1,klev,'mod_gn6hnc:ak')
      call getmem1d(bk,1,klev,'mod_gn6hnc:bk')
    else
      call getmem1d(gltemp,1,nlat,'mod_gn6hnc:gltemp')
      call getmem3d(vwork,1,nlon,1,nlat,1,klev,'mod_gn6hnc:vwork')
      call getmem3d(b2,1,nlon,1,nlat,1,klev*3,'mod_gn6hnc:b2')
      call getmem3d(d2,1,nlon,1,nlat,1,klev*2,'mod_gn6hnc:d2')
      uvar => d2(:,:,1:klev)
      vvar => d2(:,:,klev+1:2*klev)
      tvar => b2(:,:,1:klev)
      hvar => b2(:,:,klev+1:2*klev)
      qvar => b2(:,:,2*klev+1:3*klev)
      if ( dattyp == 'GFS11' ) then
        istatus = nf90_inq_dimid(inet1,'rhlev',jdim)
        call checkncerr(istatus,__FILE__,__LINE__,'Error find rhlev dim')
        istatus = nf90_inquire_dimension(inet1,jdim,len=nrhlev)
        call checkncerr(istatus,__FILE__,__LINE__,'Error inquire rhlev dim')
      else
        call getmem3d(pp3d,1,nlon,1,nlat,1,klev,'mod_gn6hnc:pp3d')
      end if
    end if
 
    istatus = nf90_inq_varid(inet1,'lat',ivar1)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find lat var')
    istatus = nf90_get_var(inet1,ivar1,glat)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read lat var')
    istatus = nf90_inq_varid(inet1,'lon',ivar1)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find lon var')
    istatus = nf90_get_var(inet1,ivar1,glon)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read lon var')

    icount(1) = nlon
    icount(2) = nlat
    icount(3) = 1
    istart(1) = 1
    istart(2) = 1
    istart(3) = 1

    npl = nipl

    if ( dattyp == 'CAM4N' .or. dattyp == 'CCSMN' ) then
      istatus = nf90_inq_varid(inet1,'hyam',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find hyam var')
      istatus = nf90_get_var(inet1,ivar1,ak)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read hyam var')
      istatus = nf90_inq_varid(inet1,'hybm',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find hybm var')
      istatus = nf90_get_var(inet1,ivar1,bk)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read hybm var')
      istatus = nf90_inq_varid(inet1,'P0',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find P0 var')
      istatus = nf90_get_var(inet1,ivar1,dp0)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read P0 var')
      p0 = real(dp0)
      istatus = nf90_inq_varid(inet1,'PHIS',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find PHIS var')
      istatus = nf90_get_var(inet1,ivar1,zsvar,istart(1:3),icount(1:3))
      call checkncerr(istatus,__FILE__,__LINE__,'Error read PHYS var')
      zsvar = zsvar/real(egrav)
      where (zsvar < 0.0) zsvar = 0.0
    else if ( dattyp(1:3) == 'HA_' ) then
      istatus = nf90_inq_varid(inet1,'lev',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find lev var')
      istatus = nf90_get_var(inet1,ivar1,ak)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read lev var')
      istatus = nf90_inq_varid(inet1,'b',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find b var')
      istatus = nf90_get_var(inet1,ivar1,bk)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read b var')
      istatus = nf90_inq_varid(inet1,'orog',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find orog var')
      istatus = nf90_get_var(inet1,ivar1,zsvar,istart(1:3),icount(1:3))
      call checkncerr(istatus,__FILE__,__LINE__,'Error read orog var')
      call getmem2d(pmslvar,1,nlon,1,nlat,'mod_gn6hnc:pmslvar')
    else if ( dattyp(1:3) == 'CA_' ) then
      istatus = nf90_inq_varid(inet1,'ap',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find ap var')
      istatus = nf90_get_var(inet1,ivar1,ak)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read ap var')
      istatus = nf90_inq_varid(inet1,'b',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find b var')
      istatus = nf90_get_var(inet1,ivar1,bk)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read b var')
      ! Close the T file, get just orography from fixed file.
      istatus = nf90_close(inet1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error close file '//trim(pathaddname))
      ! This one contains just orography.
      pathaddname = trim(inpglob)// &
            '/CanESM2/fixed/orog_fx_CanESM2_historical_r0i0p0.nc'
      istatus = nf90_open(pathaddname,nf90_nowrite,inet1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error open '//trim(pathaddname))
      istatus = nf90_inq_varid(inet1,'orog',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find orog var')
      istatus = nf90_get_var(inet1,ivar1,zsvar,istart(1:3),icount(1:3))
      call checkncerr(istatus,__FILE__,__LINE__,'Error read orog var')
    else if ( dattyp(1:3) == 'IP_' ) then
      istatus = nf90_inq_varid(inet1,'ap',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find ap var')
      istatus = nf90_get_var(inet1,ivar1,ak)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read ap var')
      istatus = nf90_inq_varid(inet1,'b',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find b var')
      istatus = nf90_get_var(inet1,ivar1,bk)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read b var')
      ! Close the T file, get just orography from fixed file.
      istatus = nf90_close(inet1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error close file '//trim(pathaddname))
      ! This one contains just orography.
      pathaddname = trim(inpglob)// &
            '/IPSL-CM5A-LR/fixed/orog_fx_IPSL-CM5A-LR_historical_r0i0p0.nc'
      istatus = nf90_open(pathaddname,nf90_nowrite,inet1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error open '//trim(pathaddname))
      istatus = nf90_inq_varid(inet1,'orog',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find orog var')
      istatus = nf90_get_var(inet1,ivar1,zsvar,istart(1:3),icount(1:3))
      call checkncerr(istatus,__FILE__,__LINE__,'Error read orog var')
    else if ( dattyp(1:3) == 'GF_' ) then
      istatus = nf90_inq_varid(inet1,'a',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find a var')
      istatus = nf90_get_var(inet1,ivar1,ak)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read a var')
      istatus = nf90_inq_varid(inet1,'b',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find b var')
      istatus = nf90_get_var(inet1,ivar1,bk)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read b var')
      istatus = nf90_inq_varid(inet1,'p0',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find p0 var')
      istatus = nf90_get_var(inet1,ivar1,dp0)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read p0 var')
      p0 = real(dp0)
      ! Close the T file, get just orography from fixed file.
      istatus = nf90_close(inet1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error close file '//trim(pathaddname))
      ! This one contains just orography.
      pathaddname = trim(inpglob)// &
            '/GFDL-ESM2M/fixed/orog_fx_GFDL-ESM2M_historical_r0i0p0.nc'
      istatus = nf90_open(pathaddname,nf90_nowrite,inet1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error open '//trim(pathaddname))
      istatus = nf90_inq_varid(inet1,'orog',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find orog var')
      istatus = nf90_get_var(inet1,ivar1,zsvar,istart(1:3),icount(1:3))
      call checkncerr(istatus,__FILE__,__LINE__,'Error read orog var')
    else if ( dattyp(1:3) == 'CN_' ) then
      istatus = nf90_inq_varid(inet1,'a',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find a var')
      istatus = nf90_get_var(inet1,ivar1,ak)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read a var')
      istatus = nf90_inq_varid(inet1,'b',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find b var')
      istatus = nf90_get_var(inet1,ivar1,bk)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read b var')
      istatus = nf90_inq_varid(inet1,'p0',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find p0 var')
      istatus = nf90_get_var(inet1,ivar1,dp0)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read p0 var')
      p0 = real(dp0)
      ! Close the T file, get just orography from fixed file.
      istatus = nf90_close(inet1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error close file '//trim(pathaddname))
      ! This one contains just orography.
      pathaddname = trim(inpglob)// &
            '/CNRM-CM5/fixed/orog_fx_CNRM-CM5_historical_r0i0p0.nc'
      istatus = nf90_open(pathaddname,nf90_nowrite,inet1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error open '//trim(pathaddname))
      istatus = nf90_inq_varid(inet1,'orog',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find orog var')
      istatus = nf90_get_var(inet1,ivar1,zsvar,istart(1:3),icount(1:3))
      call checkncerr(istatus,__FILE__,__LINE__,'Error read orog var')
    else if ( dattyp(1:3) == 'CS_' ) then
      istatus = nf90_inq_varid(inet1,'a',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find a var')
      istatus = nf90_get_var(inet1,ivar1,ak)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read a var')
      istatus = nf90_inq_varid(inet1,'b',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find b var')
      istatus = nf90_get_var(inet1,ivar1,bk)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read b var')
      istatus = nf90_inq_varid(inet1,'p0',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find p0 var')
      istatus = nf90_get_var(inet1,ivar1,dp0)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read p0 var')
      p0 = real(dp0)
      ! Close the T file, get just orography from fixed file.
      istatus = nf90_close(inet1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error close file '//trim(pathaddname))
      ! This one contains just orography.
      pathaddname = trim(inpglob)// &
            '/CSIRO-MK36/fixed/orog_fx_CSIRO-Mk3-6-0_historical_r0i0p0.nc'
      istatus = nf90_open(pathaddname,nf90_nowrite,inet1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error open '//trim(pathaddname))
      istatus = nf90_inq_varid(inet1,'orog',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find orog var')
      istatus = nf90_get_var(inet1,ivar1,zsvar,istart(1:3),icount(1:3))
      call checkncerr(istatus,__FILE__,__LINE__,'Error read orog var')
    else if ( dattyp == 'GFS11' ) then
      npl = klev ! Data are on pressure levels
      call getmem1d(pplev,1,klev,'mod_gn6hnc:pplev')
      istatus = nf90_inq_varid(inet1,'lev',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find lev var')
      istatus = nf90_get_var(inet1,ivar1,pplev)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read lev var')
      istatus = nf90_inq_varid(inet1,'orog',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find orog var')
      istatus = nf90_get_var(inet1,ivar1,zsvar,istart(1:3),icount(1:3))
      call checkncerr(istatus,__FILE__,__LINE__,'Error read orog var')
      do j = 1 , nlat
        gltemp(nlat-j+1) = glat(j)
      end do
      glat(:) = gltemp(:)
      call relmem1d(gltemp)
    else if ( dattyp(1:3) == 'EC_' ) then
      npl = klev ! Data are on pressure levels
      call getmem1d(pplev,1,klev,'mod_gn6hnc:pplev')
      istatus = nf90_inq_varid(inet1,'lev',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find lev var')
      istatus = nf90_get_var(inet1,ivar1,pplev)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read lev var')
      istatus = nf90_inq_varid(inet1,'geo',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find geo var')
      istatus = nf90_get_var(inet1,ivar1,zsvar,istart(1:3),icount(1:3))
      call checkncerr(istatus,__FILE__,__LINE__,'Error read geo var')
      ! Transform geopotential to elevation
      zsvar(:,:) = zsvar(:,:)/real(egrav)
      ! We need latitude south/north, in file is north/south
      do j = 1 , nlat
        gltemp(nlat-j+1) = glat(j)
      end do
      glat(:) = gltemp(:)
      call relmem1d(gltemp)
    end if

    istatus = nf90_close(inet1)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error close file '//trim(pathaddname))
!
    call getmem1d(sigmar,1,npl,'mod_gn6hnc:sigmar')
    call getmem3d(b3,1,jx,1,iy,1,npl*3,'mod_gn6hnc:b3')
    call getmem3d(d3,1,jx,1,iy,1,npl*2,'mod_gn6hnc:d3')

    if ( dattyp /= 'GFS11' .and. dattyp(1:3) /= 'EC_' ) then
      call getmem3d(b2,1,nlon,1,nlat,1,npl*3,'mod_gn6hnc:b2')
      call getmem3d(d2,1,nlon,1,nlat,1,npl*2,'mod_gn6hnc:d2')
      up => d2(:,:,1:npl)
      vp => d2(:,:,npl+1:2*npl)
      tp => b2(:,:,1:npl)
      hp => b2(:,:,npl+1:2*npl)
      qp => b2(:,:,2*npl+1:3*npl)
      pplev => fplev
    end if

!   Set up pointers
 
    u3 => d3(:,:,1:npl)
    v3 => d3(:,:,npl+1:2*npl)
    t3 => b3(:,:,1:npl)
    h3 => b3(:,:,npl+1:2*npl)
    q3 => b3(:,:,2*npl+1:3*npl)

    if ( dattyp == 'CAM4N' ) then
      refdate = 1989122700
      call setcal(refdate,noleap)
    end if
    timlen = 1
    call getmem1d(itimes,1,1,'mod_gn6hnc:itimes')
    itimes(1) = 1870010100 ! This set to a "Prehistorical" date
    if ( dattyp(1:3) == 'HA_' ) then
      ! HadGEM datasets has different times for PS and vertical variables.
      pstimlen = 1
      call getmem1d(ipstimes,1,1,'mod_gn6hnc:ipstimes')
      ipstimes(1) = 1870010100 ! This set to a "Prehistorical" date
      call setcal(itimes(1), y360)
      call setcal(ipstimes(1), y360)
    else if ( dattyp(1:3) == 'CS_' ) then
      ! CSIRO datasets has different times for PS and vertical variables.
      pstimlen = 1
      call getmem1d(ipstimes,1,1,'mod_gn6hnc:ipstimes')
      ipstimes(1) = 1870010100 ! This set to a "Prehistorical" date
      call setcal(itimes(1), noleap)
      call setcal(ipstimes(1), noleap)
    else if ( dattyp(1:3) == 'GFS' .or. dattyp(1:3) == 'EC_' .or. &
              dattyp(1:3) == 'CN_' ) then
      call setcal(itimes(1), gregorian)
    else
      call setcal(itimes(1), noleap)
    end if

    if ( dattyp(1:3) == 'EC_' ) then
      do k = 1 , klev
        sigmar(k) = pplev(klev-k+1)*0.00001
      end do
    else
      sigmar(:) = pplev(:)*0.001
    end if

    write (stdout,*) 'Read in Static fields OK'

  end subroutine headgn6hnc
!
!-----------------------------------------------------------------------
! 
  subroutine get_gn6hnc(idate)

    use netcdf

    implicit none
!
    type(rcm_time_and_date) , intent(in) :: idate
!
!
    call readgn6hnc(idate)
    write (stderr,*) 'Read in fields at Date: ', tochar(idate)
 
    ! GFS and EC-EARTH grids are already on pressure levels.
    if ( dattyp /= 'GFS11' .and. dattyp(1:3) /= 'EC_' ) then

      ! All processing assumes dataset in top -> bottom
      ! HadGEM is read bottom -> top
      if ( dattyp(1:3) == 'HA_' ) then
        call top2btm(tvar,nlon,nlat,klev)
        call top2btm(qvar,nlon,nlat,klev)
        call top2btm(uvar,nlon,nlat,klev)
        call top2btm(vvar,nlon,nlat,klev)
        call top2btm(pp3d,nlon,nlat,klev)
        call top2btm(hvar,nlon,nlat,klev)
      end if
 
      ! All processing assumes dataset in top -> bottom
      ! CanESM and IPSL are read bottom -> top
      if ( dattyp(1:3) == 'CA_' .or. dattyp(1:3) == 'IP_' .or. &
           dattyp(1:3) == 'GF_' .or. dattyp(1:3) == 'CN_' .or. &
           dattyp(1:3) == 'CS_' ) then
        call top2btm(tvar,nlon,nlat,klev)
        call top2btm(qvar,nlon,nlat,klev)
        call top2btm(uvar,nlon,nlat,klev)
        call top2btm(vvar,nlon,nlat,klev)
        call top2btm(pp3d,nlon,nlat,klev)
        call htsig(tvar,hvar,pp3d,psvar,zsvar,nlon,nlat,klev)
      end if

      ! Calculate HGT on Pressure levels
      call height(hp,hvar,tvar,psvar,pp3d,zsvar,nlon,nlat,klev,pplev,npl)

      ! Interpolate vertically on Pressure levels
      call intlin(up,uvar,psvar,pp3d,nlon,nlat,klev,pplev,npl)
      call intlin(vp,vvar,psvar,pp3d,nlon,nlat,klev,pplev,npl)
      call intlog(tp,tvar,psvar,pp3d,nlon,nlat,klev,pplev,npl)
      call intlin(qp,qvar,psvar,pp3d,nlon,nlat,klev,pplev,npl)
    end if
 
    ! Horizontal interpolation on RegCM grid
    call bilinx2(b3,b2,xlon,xlat,glon,glat,nlon,nlat,jx,iy,npl*3)
    call bilinx2(d3,d2,dlon,dlat,glon,glat,nlon,nlat,jx,iy,npl*2)
    ! Rotate winds
    call uvrot4(u3,v3,dlon,dlat,clon,clat,xcone,jx,iy,npl,plon,plat,iproj)
 
    ! Go to bottom->top
    if ( dattyp(1:3) /= 'EC_' ) then
      call top2btm(t3,jx,iy,npl)
      call top2btm(q3,jx,iy,npl)
      call top2btm(h3,jx,iy,npl)
      call top2btm(u3,jx,iy,npl)
      call top2btm(v3,jx,iy,npl)
    end if
 
    ! Recalculate pressure on RegCM orography
    call intgtb(pa,za,tlayer,topogm,t3,h3,sigmar,jx,iy,npl)
    call intpsn(ps4,topogm,pa,za,tlayer,ptop,jx,iy)
 
    if ( i_band == 1 ) then
       call p1p2_band(b3pd,ps4,jx,iy)
    else
       call p1p2(b3pd,ps4,jx,iy)
    endif
 
    ! Recalculate temperature on RegCM orography
    call intv3(ts4,t3,ps4,sigmar,ptop,jx,iy,npl)
    ! Replace it with SST on water points
    call readsst(ts4,idate)

    ! Vertically interpolate on RegCM sigma levels
    call intv1(u4,u3,b3pd,sigma2,sigmar,ptop,jx,iy,kz,npl)
    call intv1(v4,v3,b3pd,sigma2,sigmar,ptop,jx,iy,kz,npl)
    call intv2(t4,t3,ps4,sigma2,sigmar,ptop,jx,iy,kz,npl)
    call intv1(q4,q3,ps4,sigma2,sigmar,ptop,jx,iy,kz,npl)

    ! Get back to specific humidity
    call humid2(t4,q4,ps4,ptop,sigma2,jx,iy,kz)
 
    ! Calculate geopotential for RegCM using internal formula
    call hydrost(h4,t4,topogm,ps4,ptop,sigmaf,sigma2,dsigma,jx,iy,kz)
 
  end subroutine get_gn6hnc
!
!-----------------------------------------------------------------------
! 
  subroutine readgn6hnc(idate)
!
    use netcdf
!
    implicit none
!
    type(rcm_time_and_date) , intent(in) :: idate
!
    integer :: istatus
    integer :: i , it , itps , j , k , timid , imon1 , iyear1 , imon2 , iyear2
    character(256) :: inname

    integer :: kkrec
    character(64) :: cunit , ccal
    type(rcm_time_interval) :: tdif
    type(rcm_time_and_date) :: pdate
    integer :: year , month , day , hour , y1 , y2 , m1 , m2
    integer :: fyear , fmonth , fday , fhour
!
!
    call split_idate(idate,year,month,day,hour)
!
!   This is simpler case: just one file for each timestep with
!   all variables already on pressure levels.
!
    if ( dattyp == 'GFS11' ) then
      if ( inet(1) > 0 ) then
        istatus = nf90_close(inet(1))
        call checkncerr(istatus,__FILE__,__LINE__,'Error close file')
      end if
      write (inname,99001) year, pthsep, 'fnl_', &
           year, month, day, '_', hour, '_00.nc'
      pathaddname = trim(inpglob)//'/GFS11/'//inname
      istatus = nf90_open(pathaddname,nf90_nowrite,inet(1))
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error open '//trim(pathaddname))
      write (stdout,*) inet(1), trim(pathaddname)
      varname => gfsvars
      do kkrec = 1 , 6
        istatus = nf90_inq_varid(inet(1),trim(varname(kkrec)),ivar(kkrec))
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error find var '//trim(varname(kkrec)))
      end do
      istatus = nf90_get_var(inet(1),ivar(6),vwork(:,:,1))
      call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(6))
      do j = 1 , nlat
        psvar(:,nlat-j+1) = vwork(:,j,1)*0.01 ! Go to mb
      end do
      istatus = nf90_get_var(inet(1),ivar(1),vwork)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(1))
      do j = 1 , nlat
        tvar(:,nlat-j+1,:) = vwork(:,j,:)
      end do
      istatus = nf90_get_var(inet(1),ivar(2),vwork)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(2))
      do j = 1 , nlat
        hvar(:,nlat-j+1,:) = vwork(:,j,:)
      end do
      vwork = 0.0
      istatus = nf90_get_var(inet(1),ivar(3),vwork(:,:,klev-nrhlev+1:klev))
      call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(3))
      do j = 1 , nlat
        qvar(:,nlat-j+1,:) = vwork(:,j,:)*0.01
      end do
      istatus = nf90_get_var(inet(1),ivar(4),vwork)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(4))
      do j = 1 , nlat
        uvar(:,nlat-j+1,:) = vwork(:,j,:)
      end do
      istatus = nf90_get_var(inet(1),ivar(5),vwork)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(5))
      do j = 1 , nlat
        vvar(:,nlat-j+1,:) = vwork(:,j,:)
      end do

    ! More difficult. Multiple files per variable and per year
    else if ( dattyp(1:3) == 'EC_' ) then
      if ( idate < itimes(1) .or. idate > itimes(timlen) ) then
        do kkrec = 1 , 5
          if ( inet(kkrec) > 0 ) then
            istatus = nf90_close(inet(kkrec))
            call checkncerr(istatus,__FILE__,__LINE__,'Error close file')
          end if
        end do
        varname => echvars
        do kkrec = 1 , 5
          if ( dattyp(4:4) == 'R' ) then
            write (inname,99006) 'RF', pthsep, year, pthsep, 'ich1_', &
                  trim(varname(kkrec))//'_', year, '.nc'
          else
            write (inname,99006) ('RCP'//dattyp(4:5)), pthsep, year,  &
                  pthsep, 'ich1_',                                    &
                  trim(varname(kkrec))//'_', year, '.nc'
          end if
          pathaddname = trim(inpglob)//'/EC-EARTH/'//inname
          istatus = nf90_open(pathaddname,nf90_nowrite,inet(kkrec))
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error open '//trim(pathaddname))
          istatus = nf90_inq_varid(inet(kkrec),trim(varname(kkrec)),ivar(kkrec))
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error find var '//trim(varname(kkrec)))
          write (stdout,*) inet(kkrec), trim(pathaddname)
        end do
        istatus = nf90_inq_dimid(inet(1),'time',timid)
        call checkncerr(istatus,__FILE__,__LINE__,'Error find dim time')
        istatus = nf90_inquire_dimension(inet(1),timid, len=timlen)
        call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dim time')
        istatus = nf90_inq_varid(inet(1),'time',timid)
        call checkncerr(istatus,__FILE__,__LINE__,'Error find var time')
        istatus = nf90_get_att(inet(1),timid,'units',cunit)
        call checkncerr(istatus,__FILE__,__LINE__,'Error read time units')
        istatus = nf90_get_att(inet(1),timid,'calendar',ccal)
        call checkncerr(istatus,__FILE__,__LINE__,'Error read time calendar')
        call getmem1d(itimes,1,timlen,'mod_gn6hnc:itimes')
        call getmem1d(xtimes,1,timlen,'mod_gn6hnc:xtimes')
        if (istatus /= 0) call die('mod_gn6hnc','Allocation error on itimes',1)
        istatus = nf90_get_var(inet(1),timid,xtimes)
        call checkncerr(istatus,__FILE__,__LINE__,'Error read time')
        do it = 1 , timlen
          itimes(it) = timeval2date(xtimes(it),cunit,ccal)
        end do
      end if
      tdif = idate - itimes(1)
      it = idnint(tohours(tdif))/6 + 1
      icount(1) = nlon
      icount(2) = nlat
      icount(3) = klev
      icount(4) = 1
      istart(1) = 1
      istart(2) = 1
      istart(3) = 1
      istart(4) = it
      istatus = nf90_get_var(inet(1),ivar(1),vwork,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(1))
      do j = 1 , nlat
        tvar(:,nlat-j+1,:) = vwork(:,j,:)
      end do
      istatus = nf90_get_var(inet(2),ivar(2),vwork,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(2))
      do j = 1 , nlat
        hvar(:,nlat-j+1,:) = vwork(:,j,:)/real(egrav) ! Get to m
      end do
      where ( hvar < 0.0 )
        hvar = 0.0 ! We do not want negative hgt, don't we?
      end where
      ! This is specific humidity
      istatus = nf90_get_var(inet(3),ivar(3),vwork,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(3))
      do j = 1 , nlat
        qvar(:,nlat-j+1,:) = vwork(:,j,:)
      end do
      istatus = nf90_get_var(inet(4),ivar(4),vwork,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(4))
      do j = 1 , nlat
        uvar(:,nlat-j+1,:) = vwork(:,j,:)
      end do
      istatus = nf90_get_var(inet(5),ivar(5),vwork,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(5))
      do j = 1 , nlat
        vvar(:,nlat-j+1,:) = vwork(:,j,:)
      end do
      do k = 1, klev
        pp3d(:,:,k) = pplev(k)*0.01 ! Get in hPa
      end do
      ! Replace with relative humidity for internal calculation
      call humid1fv(tvar,qvar,pp3d,nlon,nlat,klev)
    else
      ! Even more difficult. Each data type has its own quirks.
      tdif = rcm_time_interval(180,uhrs)
      ! HadGEM dataset has PS files with different times.
      if ( dattyp(1:3) == 'HA_' ) then
        if ( idate < ipstimes(1) .or. idate > ipstimes(pstimlen) ) then
          if ( inet(6) > 0 ) then
            istatus = nf90_close(inet(6))
            call checkncerr(istatus,__FILE__,__LINE__,'Error close file')
          end if
          iyear1 = year
          if ( month == 12 .and. day >= 1 .and. hour > 0 ) then
            iyear1 = iyear1 + 1
          end if
          if ( dattyp(4:4) == 'R' ) then
            write (inname,99005) 'RF', pthsep, trim(havars(6)), pthsep, &
                 trim(havars(6)), trim(hapbase1)//trim(habase3), &
                 iyear1-1, '12010600-', iyear1, '12010000.nc'
          else
            if (year*1000000+month*10000+day*100+hour == 2005120100) then 
              write (inname,99005) ('RCP'//dattyp(4:5)), pthsep, &
                trim(havars(6)), pthsep, trim(havars(6)), &
                trim(hapbase2)//dattyp(4:5)//trim(habase3), &
                iyear1, '12010000-', iyear1+1, '12010000.nc'
            else
              write (inname,99005) ('RCP'//dattyp(4:5)), pthsep, &
                trim(havars(6)), pthsep, trim(havars(6)), &
                trim(hapbase2)//dattyp(4:5)//trim(habase3), &
                iyear1-1, '12010600-', iyear1, '12010000.nc'
            end if
          end if
          pathaddname = trim(inpglob)//'/HadGEM2/'//inname
          istatus = nf90_open(pathaddname,nf90_nowrite,inet(6))
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error open '//trim(pathaddname))
          istatus = nf90_inq_dimid(inet(6),'time',timid)
          call checkncerr(istatus,__FILE__,__LINE__,'Error find dim time')
          istatus = nf90_inquire_dimension(inet(6),timid, len=pstimlen)
          call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dim time')
          istatus = nf90_inq_varid(inet(6),'time',timid)
          call checkncerr(istatus,__FILE__,__LINE__,'Error find var time')
          istatus = nf90_get_att(inet(6),timid,'units',cunit)
          call checkncerr(istatus,__FILE__,__LINE__,'Error read time units')
          istatus = nf90_get_att(inet(6),timid,'calendar',ccal)
          call checkncerr(istatus,__FILE__,__LINE__,'Error read time calendar')
          call getmem1d(ipstimes,1,pstimlen,'mod_gn6hnc:ipstimes')
          call getmem1d(xtimes,1,pstimlen,'mod_gn6hnc:xtimes')
          if (istatus /= 0) call die('mod_gn6hnc','Allocation error on itimes',1)
          istatus = nf90_get_var(inet(6),timid,xtimes)
          call checkncerr(istatus,__FILE__,__LINE__,'Error read time')
          do it = 1 , pstimlen
            ipstimes(it) = timeval2date(xtimes(it),cunit,ccal)
          end do
          istatus = nf90_inq_varid(inet(6), trim(havars(6)), ivar(6))
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error find var '//trim(havars(6)))
          write (stdout,*) inet(6), trim(pathaddname)
        end if
      end if
      ! CSIRO dataset has PS files with different times.
      if ( dattyp(1:3) == 'CS_' ) then
        if ( idate < ipstimes(1) .or. idate > ipstimes(pstimlen) ) then
          if ( inet(6) > 0 ) then
            istatus = nf90_close(inet(6))
            call checkncerr(istatus,__FILE__,__LINE__,'Error close file')
          end if
          iyear1 = year/5*5
          if ( mod(year,5) == 0 .and. month == 1 .and. &
               day == 1 .and. hour == 0 ) then
            iyear1 = iyear1 - 5
          end if
          if ( dattyp(4:4) == 'R' ) then
            write (inname,99005) 'RF', pthsep, trim(csirvars(6)), pthsep, &
                 trim(csirvars(6)), trim(csirbase1)//trim(csirbase3), &
                 iyear1, '01010600-', iyear1+5, '01010000.nc'
          else
            if (year*1000000+month*10000+day*100+hour == 2005120100) then 
              write (inname,99005) ('RCP'//dattyp(4:5)), pthsep, &
                trim(csirvars(6)), pthsep, trim(csirvars(6)), &
                trim(csirbase2)//dattyp(4:5)//trim(csirbase3), &
                iyear1, '01010600-', iyear1+5, '01010000.nc'
            else
              write (inname,99005) ('RCP'//dattyp(4:5)), pthsep, &
                trim(csirvars(6)), pthsep, trim(csirvars(6)), &
                trim(csirbase2)//dattyp(4:5)//trim(csirbase3), &
                iyear1, '01010600-', iyear1+5, '01010000.nc'
            end if
          end if
          pathaddname = trim(inpglob)//'/CSIRO-MK36/'//inname
          istatus = nf90_open(pathaddname,nf90_nowrite,inet(6))
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error open '//trim(pathaddname))
          istatus = nf90_inq_dimid(inet(6),'time',timid)
          call checkncerr(istatus,__FILE__,__LINE__,'Error find dim time')
          istatus = nf90_inquire_dimension(inet(6),timid, len=pstimlen)
          call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dim time')
          istatus = nf90_inq_varid(inet(6),'time',timid)
          call checkncerr(istatus,__FILE__,__LINE__,'Error find var time')
          istatus = nf90_get_att(inet(6),timid,'units',cunit)
          call checkncerr(istatus,__FILE__,__LINE__,'Error read time units')
          istatus = nf90_get_att(inet(6),timid,'calendar',ccal)
          call checkncerr(istatus,__FILE__,__LINE__,'Error read time calendar')
          call getmem1d(ipstimes,1,pstimlen,'mod_gn6hnc:ipstimes')
          call getmem1d(xtimes,1,pstimlen,'mod_gn6hnc:xtimes')
          if (istatus /= 0) call die('mod_gn6hnc','Allocation error on itimes',1)
          istatus = nf90_get_var(inet(6),timid,xtimes)
          call checkncerr(istatus,__FILE__,__LINE__,'Error read time')
          do it = 1 , pstimlen
            ipstimes(it) = timeval2date(xtimes(it),cunit,ccal)
          end do
          istatus = nf90_inq_varid(inet(6), trim(csirvars(6)), ivar(6))
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error find var '//trim(csirvars(6)))
          write (stdout,*) inet(6), trim(pathaddname)
        end if
      end if

      if ( idate < itimes(1) .or. idate > itimes(timlen) ) then
        if (inet(1) > 0) then
          if ( dattyp == 'CAM4N' ) then
            istatus = nf90_close(inet(1))
            call checkncerr(istatus,__FILE__,__LINE__,'Error close file')
            filedate = filedate + tdif
          else
            if ( dattyp(1:3) == 'HA_' .or. dattyp(1:3) == 'CS_' ) then
              do i = 1 , nfiles-1
                if ( havars(i) /= 'XXX' ) then
                  istatus = nf90_close(inet(i))
                  call checkncerr(istatus,__FILE__,__LINE__,'Error close file')
                end if
              end do
            else
              do i = 1 , nfiles
                if ( havars(i) /= 'XXX' ) then
                  istatus = nf90_close(inet(i))
                  call checkncerr(istatus,__FILE__,__LINE__,'Error close file')
                end if
              end do
            end if
          end if
        else
          if ( dattyp == 'CAM4N' ) then
            pdate = refdate
            filedate = refdate
            do while (idate >= pdate)
              filedate = pdate
              pdate = pdate + tdif
            end do
          end if
        end if

        if ( dattyp == 'CAM4N' ) then
          call split_idate(filedate,fyear,fmonth,fday,fhour)
          ! All variables just in a single file. Simpler case.
          write (inname,99002) trim(cambase) , fyear, &
                    fmonth, fday, filedate%second_of_day
          pathaddname = trim(inpglob)//'/CAM2/'//trim(inname)
          istatus = nf90_open(pathaddname,nf90_nowrite,inet(1))
          call checkncerr(istatus,__FILE__,__LINE__, &
                           'Error open '//trim(pathaddname))
          write (stdout,*) inet(1), trim(pathaddname)
          inet(2:nfiles) = inet(1)
          varname => cam2vars
        else if ( dattyp == 'CCSMN' ) then
          ! Dataset prepared in mothly files, one for each variable
          do i = 1 , nfiles
            write (inname,99003) year, trim(ccsmfname(i)), &
                      mname(month), year
            pathaddname = trim(inpglob)//'/CCSM/'//trim(inname)
            istatus = nf90_open(pathaddname,nf90_nowrite,inet(i))
            call checkncerr(istatus,__FILE__,__LINE__, &
                            'Error open '//trim(pathaddname))
            write (stdout,*) inet(i), trim(pathaddname)
          end do
          varname => ccsmvars
        else if ( dattyp(1:3) == 'HA_' ) then
          ! 3 months dataset per file.
          do i = 1 , nfiles-1
            if ( havars(i) /= 'XXX' ) then
              iyear1 = year
              imon1 = month
              imon1 = imon1 / 3 * 3
              if ( day == 1 .and. hour == 0 ) then
                if ( mod(month,3) == 0 ) then
                  imon1 = imon1 - 3
                end if
              end if
              if ( imon1 == 0 ) then
                imon1 = 12
                iyear1 = iyear1 - 1
              end if
              if ( imon1 == 12 ) then
                iyear2 = iyear1 + 1
                imon2  = 3
              else
                iyear2 = iyear1
                imon2 = imon1 + 3
              end if
              if ( dattyp(4:4) == 'R' ) then
                write (inname,99004) 'RF', pthsep, trim(havars(i)), pthsep, &
                  trim(havars(i)), trim(habase1)//trim(habase3), &
                  iyear1, imon1 , '010600-', iyear2, imon2, '010000.nc'
              else
                if (iyear1*1000000+imon1*10000 == 2005090000) then 
                  iyear1 = 2005
                  imon1 = 12
                  iyear2 = 2006
                  imon2 = 3 
                  write (inname,99004) ('RCP'//dattyp(4:5)), pthsep, &
                    trim(havars(i)), pthsep, trim(havars(i)), &
                    trim(habase2)//dattyp(4:5)//trim(habase3), &
                    iyear1, imon1, '010000-', iyear2, imon2, '010000.nc'
                else
                  write (inname,99004) ('RCP'//dattyp(4:5)), pthsep, &
                    trim(havars(i)), pthsep, trim(havars(i)), &
                    trim(habase2)//dattyp(4:5)//trim(habase3), &
                    iyear1, imon1, '010600-', iyear2, imon2, '010000.nc'
                end if
              end if
              pathaddname = trim(inpglob)//'/HadGEM2/'//inname
              istatus = nf90_open(pathaddname,nf90_nowrite,inet(i))
              call checkncerr(istatus,__FILE__,__LINE__, &
                              'Error open '//trim(pathaddname))
              write (stdout,*) inet(i), trim(pathaddname)
            end if
          end do
          varname => havars
        else if ( dattyp(1:3) == 'CA_' ) then
          ! yearly files, one for each variable
          do i = 1 , nfiles
            if ( cavars(i) /= 'XXX' ) then
              if ( dattyp(4:4) == 'R' ) then
                write (inname,99005) 'RF', pthsep, trim(cavars(i)), pthsep, &
                  trim(cavars(i)), trim(cabase1)//trim(cabase3), &
                  year, '01010000-', year, '12311800.nc'
              else
                write (inname,99005) ('RCP'//dattyp(4:5)), pthsep, &
                  trim(cavars(i)), pthsep, trim(cavars(i)), &
                  trim(cabase2)//dattyp(4:5)//trim(cabase3), &
                  year, '01010000-', year, '12311800.nc'
              end if
              pathaddname = trim(inpglob)//'/CanESM2/'//inname
              istatus = nf90_open(pathaddname,nf90_nowrite,inet(i))
              call checkncerr(istatus,__FILE__,__LINE__, &
                              'Error open '//trim(pathaddname))
            end if
            write (stdout,*) inet(i), trim(pathaddname)
          end do
          varname => cavars
        else if ( dattyp(1:3) == 'IP_' ) then
          ! 10 yearly files, one for each variable
          do i = 1 , nfiles-1
            if ( ipvars(i) /= 'XXX' ) then
              if ( dattyp(4:4) == 'R' ) then
                if ( year > 2000 ) then
                  write (inname,99005) 'HIST',pthsep,trim(ipvars(i)), pthsep, &
                    trim(ipvars(i)), trim(ipbase1)//trim(ipbase3), &
                    2000, '01010300-', 2005, '12312100.nc'
                else
                  write (inname,99005) 'HIST',pthsep,trim(ipvars(i)), pthsep, &
                    trim(ipvars(i)), trim(ipbase1)//trim(ipbase3), &
                    (year/10*10), '01010300-', (year/10*10+9), '12312100.nc'
                end if
              else
                write (inname,99005) ('RCP'//dattyp(4:5)), pthsep, &
                  trim(ipvars(i)), pthsep, trim(ipvars(i)), &
                  trim(ipbase2)//dattyp(4:5)//trim(ipbase3), &
                  (year/10*10+6), '01010300-', (year/10*10+15), '12312100.nc'
              end if
              pathaddname = trim(inpglob)//'/IPSL-CM5A-LR/'//inname
              istatus = nf90_open(pathaddname,nf90_nowrite,inet(i))
              call checkncerr(istatus,__FILE__,__LINE__, &
                              'Error open '//trim(pathaddname))
              write (stdout,*) inet(i), trim(pathaddname)
            end if
          end do
          ! PS is in 50 years file...
          if ( dattyp(4:4) == 'R' ) then
            if ( year < 1950 ) then
              write (inname,99005) 'HIST',pthsep,trim(ipvars(6)), pthsep, &
                    trim(ipvars(6)), trim(ipbase1)//trim(ipbase3), &
                    1900, '01010300-', 1949, '12312100.nc'
            else if ( year < 2000 ) then
              write (inname,99005) 'HIST',pthsep,trim(ipvars(6)), pthsep, &
                    trim(ipvars(6)), trim(ipbase1)//trim(ipbase3), &
                    1950, '01010300-', 1999, '12312100.nc'
            else
              write (inname,99005) 'HIST',pthsep,trim(ipvars(6)), pthsep, &
                    trim(ipvars(6)), trim(ipbase1)//trim(ipbase3), &
                    2000, '01010300-', 2005, '12312100.nc'
            end if
          else
            write (inname,99005) ('RCP'//dattyp(4:5)), pthsep, &
               trim(ipvars(6)), pthsep, trim(ipvars(6)), &
               trim(ipbase2)//dattyp(4:5)//trim(ipbase3), &
               (year/50*50+6), '01010300-', (year/50*50+55), '12312100.nc'
          end if
          pathaddname = trim(inpglob)//'/IPSL-CM5A-LR/'//inname
          istatus = nf90_open(pathaddname,nf90_nowrite,inet(6))
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error open '//trim(pathaddname))
          write (stdout,*) inet(6), trim(pathaddname)
          varname => ipvars
        else if ( dattyp(1:3) == 'GF_' ) then
          ! 5 yearly files, one for each variable
          y1 = (year-1)/5*5+1
          y2 = y1+4
          if ( year == y1 .and. month == 1 .and. &
               day == 1 .and. hour == 0 ) then
            y1 = y1-5
            y2 = y1+4
          end if
          do i = 1 , nfiles
            if ( gfdlvars(i) /= 'XXX' ) then
              if ( dattyp(4:4) == 'R' ) then
                write (inname,99005) 'RF',pthsep,trim(gfdlvars(i)), pthsep, &
                  trim(gfdlvars(i)), trim(gfdlbase1)//trim(gfdlbase3), &
                  y1, '010100-', y2, '123123.nc'
              else
              if (year*1000000+month*10000+day*100+hour == 2006010100) then
                ! use modified file (r45 first time step is added to hist last time step)
                  write (inname,99005) ('RCP'//dattyp(4:4)//'.'//dattyp(5:5)), &
                  pthsep, trim(gfdlvars(i)), pthsep, trim(gfdlvars(i)), &
                  trim(gfdlbase2)//dattyp(4:5)//trim(gfdlbase3), &
                  2006, '010100-', 2010, '123123_fixed.nc'
              else
                write (inname,99005) ('RCP'//dattyp(4:4)//'.'//dattyp(5:5)), &
                  pthsep, trim(gfdlvars(i)), pthsep, trim(gfdlvars(i)), &
                  trim(gfdlbase2)//dattyp(4:5)//trim(gfdlbase3), &
                  y1, '010100-', y2, '123123.nc'
              end if
              end if
              pathaddname = trim(inpglob)//'/GFDL-ESM2M/'//inname
              istatus = nf90_open(pathaddname,nf90_nowrite,inet(i))
              call checkncerr(istatus,__FILE__,__LINE__, &
                              'Error open '//trim(pathaddname))
              write (stdout,*) inet(i), trim(pathaddname)
            end if
          end do
          varname => gfdlvars
        else if ( dattyp(1:3) == 'CN_' ) then
          ! monthly files, one for each variable
          if ( day == 1 .and. hour == 0 ) then
            m1 = month-1
            y1 = year
            if ( m1 < 1 ) then
              m1 = 12
              y1 = year-1
            end if
            m2 = m1+1
            y2 = y1
          else
            m1 = month
            y1 = year
            m2 = m1+1
            y2 = y1
          end if
          if ( m2 > 12 ) then
            m2 = 1
            y2 = y2+1
          end if
          do i = 1 , nfiles
            if ( cnrmvars(i) /= 'XXX' ) then
              if ( dattyp(4:4) == 'R' ) then
                write (inname,99007) 'RF',pthsep,trim(cnrmvars(i)), pthsep, &
                  trim(cnrmvars(i)), trim(cnrmbase1)//trim(cnrmbase3), &
                  y1, m1, '010600-', y2, m2, '010000.nc'
              else
                write (inname,99007) ('RCP'//dattyp(4:4)//'.'//dattyp(5:5)), &
                  pthsep, trim(cnrmvars(i)), pthsep, trim(cnrmvars(i)), &
                  trim(cnrmbase2)//dattyp(4:5)//trim(cnrmbase3), &
                  y1, m1, '010600-', y2, m2, '010000.nc'
              end if
              pathaddname = trim(inpglob)//'/CNRM-CM5/'//inname
              istatus = nf90_open(pathaddname,nf90_nowrite,inet(i))
              call checkncerr(istatus,__FILE__,__LINE__, &
                              'Error open '//trim(pathaddname))
              write (stdout,*) inet(i), trim(pathaddname)
            end if
          end do
          varname => cnrmvars
        else if ( dattyp(1:3) == 'CS_' ) then
          ! yearly files, one for each variable
          y1 = year
          if ( month == 1 .and. day == 1 .and. hour == 0 ) then
            y1 = year-1
          end if
          y2 = y1+1
          do i = 1 , nfiles-1
            if ( csirvars(i) /= 'XXX' ) then
              if ( dattyp(4:4) == 'R' ) then
                write (inname,99005) 'RF',pthsep,trim(csirvars(i)), pthsep, &
                  trim(csirvars(i)), trim(csirbase1)//trim(csirbase3), &
                  y1, '01010600-', y2, '01010000.nc'
              else
                write (inname,99005) ('RCP'//dattyp(4:4)//'.'//dattyp(5:5)), &
                  pthsep, trim(csirvars(i)), pthsep, trim(csirvars(i)), &
                  trim(csirbase2)//dattyp(4:5)//trim(csirbase3), &
                  y1, '01010600-', y2, '01010000.nc'
              end if
              pathaddname = trim(inpglob)//'/CSIRO-MK36/'//inname
              istatus = nf90_open(pathaddname,nf90_nowrite,inet(i))
              call checkncerr(istatus,__FILE__,__LINE__, &
                              'Error open '//trim(pathaddname))
              write (stdout,*) inet(i), trim(pathaddname)
            end if
          end do
          varname => csirvars
        end if
      end if

      istatus = nf90_inq_dimid(inet(1),'time',timid)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find dim time')
      istatus = nf90_inquire_dimension(inet(1),timid, len=timlen)
      call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dim time')
      istatus = nf90_inq_varid(inet(1),'time',timid)
      call checkncerr(istatus,__FILE__,__LINE__,'Error find var time')
      istatus = nf90_get_att(inet(1),timid,'units',cunit)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read time units')
      istatus = nf90_get_att(inet(1),timid,'calendar',ccal)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read time calendar')
      call getmem1d(itimes,1,timlen,'mod_gn6hnc:itimes')
      call getmem1d(xtimes,1,timlen,'mod_gn6hnc:xtimes')
      if (istatus /= 0) call die('mod_gn6hnc','Allocation error on itimes',1)
      istatus = nf90_get_var(inet(1),timid,xtimes)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read time')
      do it = 1 , timlen
        itimes(it) = timeval2date(xtimes(it),cunit,ccal)
      end do
      do kkrec = 1 , 6
        if ( varname(kkrec) /= 'XXX' ) then
          istatus = nf90_inq_varid(inet(kkrec), &
                          trim(varname(kkrec)), ivar(kkrec))
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error find var '//trim(varname(kkrec)))
        end if
      end do

      tdif = idate - itimes(1)
      it = idnint(tohours(tdif))/6 + 1

      icount(1) = nlon
      icount(2) = nlat
      icount(3) = 1
      istart(1) = 1
      istart(2) = 1
      istart(3) = it
      if ( dattyp(1:3) == 'HA_' ) then
        tdif = idate - ipstimes(1)
        itps = idnint(tohours(tdif))/6 + 1
        istart(3) = itps
        istatus = nf90_get_var(inet(6),ivar(6),pmslvar,istart(1:3),icount(1:3))
        call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(6))
        pmslvar(:,:) = pmslvar(:,:)*0.01
      else if ( dattyp(1:3) == 'CS_' ) then
        tdif = idate - ipstimes(1)
        itps = idnint(tohours(tdif))/6 + 1
        istart(3) = itps
        istatus = nf90_get_var(inet(6),ivar(6),psvar,istart(1:3),icount(1:3))
        call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(6))
      else
        istatus = nf90_get_var(inet(6),ivar(6),psvar,istart(1:3),icount(1:3))
        call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(6))
      end if
      icount(1) = nlon
      icount(2) = nlat
      icount(3) = klev
      icount(4) = 1
      istart(1) = 1
      istart(2) = 1
      istart(3) = 1
      istart(4) = it
      istatus = nf90_get_var(inet(1),ivar(1),tvar,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(1))
      if ( dattyp == 'CAM4N' .or. dattyp == 'CCSMN' ) then
        ! We have geopotential HGT in m, on hybrid sigma pressure levels
        istatus = nf90_get_var(inet(2),ivar(2),hvar,istart,icount)
        call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(2))
        do k = 1 , klev
          do j = 1 , nlat
            do i = 1 , nlon
              pp3d(i,j,k) = (ak(k)*p0 + bk(k)*psvar(i,j))
            end do
          end do
        end do
        psvar(:,:) = psvar(:,:)*0.01
        pp3d(:,:,:) = pp3d(:,:,:)*0.01
      else if ( dattyp(1:3) == 'HA_' ) then
        ! Data are on sigma Z levels, and we have MSLP instead of PS
        do k = 1 , klev
          hvar(:,:,k) = ak(k) + bk(k)*zsvar(:,:)
        end do
        call mslp2ps(hvar,tvar,pmslvar,zsvar,psvar,nlon,nlat,klev)
        call psig(tvar,hvar,pp3d,psvar,zsvar,nlon,nlat,klev)
      else if ( dattyp(1:3) == 'CA_' .or. dattyp(1:3) == 'IP_' ) then
        ! Data on sigma P levels, the factor for ak is ipotized.
        ! Units in file is Pa, but calculations suggest mPa
        do k = 1, klev
          pp3d(:,:,k) = ak(k)*0.001 + bk(k)*psvar(:,:)
        end do
        psvar(:,:) = psvar(:,:)*0.01
        pp3d(:,:,:) = pp3d(:,:,:)*0.01
      else if ( dattyp(1:3) == 'GF_' .or. dattyp(1:3) == 'CN_' .or. &
                dattyp(1:3) == 'CS_' ) then
        do k = 1, klev
          pp3d(:,:,k) = ak(k)*p0 + bk(k)*psvar(:,:)
        end do
        ! Go to hPa
        psvar(:,:) = psvar(:,:)*0.01
        pp3d(:,:,:) = pp3d(:,:,:)*0.01
      end if
      istatus = nf90_get_var(inet(3),ivar(3),qvar,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(3))

      ! Replace with relative humidity for internal calculation
      call humid1fv(tvar,qvar,pp3d,nlon,nlat,klev)

      istatus = nf90_get_var(inet(4),ivar(4),uvar,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(4))
      if ( dattyp(1:3) == 'HA_' ) then
        ! Data is missing on poles.
        icount(2) = nlat-1
        istatus = nf90_get_var(inet(5),ivar(5),vwork,istart,icount)
        call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(5))
        do k = 1 , klev
          vvar(:,1,k) = vwork(:,1,k)
          do j = 1 , nlat-2
            vvar(:,j,k) = 0.5*(vwork(:,j,k)+vwork(:,j+1,k))
          end do
          vvar(:,nlat,k) = vwork(:,nlat-1,k)
        end do
      else
        istatus = nf90_get_var(inet(5),ivar(5),vvar,istart,icount)
        call checkncerr(istatus,__FILE__,__LINE__,'Error read var '//varname(5))
      end if

    end if ! Other data types not the GFS11 or E_ICH1

99001   format (i0.4,a,a,i0.4,i0.2,i0.2,a,i0.2,a)
99002   format (a,i0.4,'-',i0.2,'-',i0.2,'-',i0.5,'.nc')
99003   format (i0.4,'/','ccsm.',a,a,'.',i0.4,'.nc')
99004   format (a,a,a,a,a,a,i0.4,i0.2,a,i0.4,i0.2,a)
99005   format (a,a,a,a,a,a,i0.4,a,i0.4,a)
99006   format (a,a,i0.4,a,a,a,i0.4,a)
99007   format (a,a,a,a,a,a,i0.4,i0.2,a,i0.4,i0.2,a)

  end subroutine readgn6hnc
!
end module mod_gn6hnc

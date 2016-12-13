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

  use mod_intkinds
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
  use mod_ccsm3_helper
  use mod_hadgem_helper
  use mod_csiro_helper
  use mod_canesm_helper
  use mod_miroc_helper
  use mod_ipsl_helper
  use mod_gfdl_helper
  use mod_cnrm_helper
  use mod_mpiesm_helper

  private

  ! Dimension of input read from input files
  integer(ik4) :: nlon , nlat , nulon , nvlat , klev

  ! Pressure levels to interpolate to if dataset is on model sigma levels.
  integer(ik4) , parameter :: nipl = 38
  real(rkx) , target , dimension(nipl) :: fplev = &
   (/  1.0_rkx,   2.0_rkx,   3.0_rkx,   5.0_rkx,   7.0_rkx,  10.0_rkx, &
      20.0_rkx,  30.0_rkx,  50.0_rkx,  70.0_rkx, 100.0_rkx, 125.0_rkx, &
     150.0_rkx, 175.0_rkx, 200.0_rkx, 225.0_rkx, 250.0_rkx, 300.0_rkx, &
     350.0_rkx, 400.0_rkx, 425.0_rkx, 450.0_rkx, 500.0_rkx, 550.0_rkx, &
     600.0_rkx, 650.0_rkx, 700.0_rkx, 750.0_rkx, 775.0_rkx, 800.0_rkx, &
     825.0_rkx, 850.0_rkx, 875.0_rkx, 900.0_rkx, 925.0_rkx, 950.0_rkx, &
     975.0_rkx, 1000.0_rkx /)

  integer(ik4) :: npl , nrhlev
  real(rkx) , pointer , dimension(:) :: pplev
  real(rkx) , pointer , dimension(:) :: sigmar
  real(rkx) :: pss

  ! Whole space
  real(rkx) , pointer , dimension(:,:,:) :: b2
  real(rkx) , pointer , dimension(:,:,:) :: d2
  real(rkx) , pointer , dimension(:,:,:) :: b3
  real(rkx) , pointer , dimension(:,:,:) :: d3
  real(rkx) , pointer , dimension(:,:,:) :: ha_d2_1
  real(rkx) , pointer , dimension(:,:,:) :: ha_d2_2

  real(rkx) , pointer :: u3(:,:,:) , v3(:,:,:)
  real(rkx) , pointer :: h3(:,:,:) , q3(:,:,:) , t3(:,:,:)
  real(rkx) , pointer :: up(:,:,:) , vp(:,:,:)
  real(rkx) , pointer :: hp(:,:,:) , qp(:,:,:) , tp(:,:,:)

  ! Input space
  real(rkx) :: p0
  real(rkx) , pointer , dimension(:,:) :: psvar , zsvar ! , pmslvar
  real(rkx) , pointer , dimension(:) :: ak , bk
  real(rkx) , pointer , dimension(:) :: glat , gltemp , uglon
  real(rkx) , pointer , dimension(:) :: glon , vglat
  real(rkx) , pointer , dimension(:,:) :: glat2
  real(rkx) , pointer , dimension(:,:) :: glon2
  real(rkx) , pointer , dimension(:,:,:) :: hvar , qvar , tvar , &
                                           uvar , vvar , pp3d , vwork
  integer(ik4) :: timlen , pstimlen
  type(rcm_time_and_date) , pointer , dimension(:) :: itimes
  type(rcm_time_and_date) , pointer , dimension(:) :: ipstimes
  real(rkx) , pointer , dimension(:) :: xtimes

  ! Shared by netcdf I/O routines
  integer(ik4) , dimension(4) :: icount , istart
  ! We will need 6 files (is just one for CAM2)
  integer(ik4) , parameter :: nvars = 6
  integer(ik4) , parameter :: nfiles = nvars
  integer(ik4) , dimension(nvars) :: inet
  integer(ik4) , dimension(nvars) :: ivar

  public :: get_gn6hnc , headgn6hnc

  character(len=256) :: pathaddname
  type(rcm_time_and_date) , save :: refdate
  type(rcm_time_and_date) , save :: filedate

  data inet /nvars*-1/

  character(len=32) :: cambase = 'sococa.ts1.r1.cam2.h1.'

  character(len=3) , target , dimension(nvars) :: cam2vars = &
                         (/'T  ' , 'Z3 ' , 'Q  ' , 'U  ' , 'V  ' , 'PS '/)
  character(len=3) , target , dimension(nvars) :: ccsmvars = &
                         (/'T  ' , 'Z3 ' , 'Q  ' , 'U  ' , 'V  ' , 'PS '/)
  character(len=3) , target , dimension(nvars) :: echvars = &
                         (/'t  ' , 'z  ' , 'q  ' , 'u  ' , 'v  ' , 'XXX'/)

  character(len=3) , target , dimension(nvars) :: gfsvars = &
                         (/'ta ' , 'hga' , 'rha' , 'ua ' , 'va ' , 'ps '/)
  character(len=3) , target , dimension(nvars) :: ec5vars = &
                         (/'ta ' , 'gpa' , 'rha' , 'ua ' , 'va ' , 'XXX'/)
  character(len=5) , target , dimension(nvars) :: jra55vars = &
               (/'var11' , 'var7 ' , 'var52' , 'var33' , 'var34' , 'XXX  '/)

  character(len=4) , dimension(nvars) :: ccsmfname = &
               (/'air ' , 'hgt ' , 'shum' , 'uwnd' , 'vwnd' , 'pres'/)

  character(len=6) , target , dimension(nvars) :: ec5name = &
    (/'STP   ' , 'GPH   ' , 'RELHUM' , 'U     ' , 'V     ' , 'XXX   '/)
  character(len=8) , target , dimension(nvars) :: jra55name = &
    (/'011_tmp ','007_hgt ','052_rh  ','033_ugrd','034_vgrd','XXXXXXXX'/)

  character(len=3) , dimension(12) :: mname = &
                         (/'JAN','FEB','MAR','APR','MAY','JUN', &
                           'JUL','AUG','SEP','OCT','NOV','DEC'/)

  character(len=3) , dimension(:) , pointer :: varname

  contains

  subroutine headgn6hnc
    use netcdf
    implicit none
    integer(ik4) :: istatus , ivar1 , inet1 , inet2 , inet3 , jdim , i , j , k
    character(len=256) :: pathaddname
    real(8) :: dp0

    if ( dattyp == 'CAM4N' ) then
      pathaddname = trim(inpglob)// &
            '/CAM2/USGS-gtopo30_0.9x1.25_remap_c051027.nc'
    else if ( dattyp == 'CCSMN' ) then
      pathaddname = trim(inpglob)//'/CCSM/ccsm_ht.nc'
    else if ( dattyp(1:3) == 'HA_' ) then
      ! Vertical info are not stored in the fixed orography file.
      ! Read part of the info from a T file.
      call find_hadgem_dim(pathaddname)
    else if ( dattyp(1:3) == 'CA_' ) then
      ! Vertical info are not stored in the fixed orography file.
      ! Read part of the info from a T file.
      call find_canesm_dim(pathaddname)
    else if ( dattyp(1:3) == 'IP_' ) then
      ! Vertical info are not stored in the fixed orography file.
      ! Read part of the info from a T file.
      call find_ipsl_dim(pathaddname)
    else if ( dattyp(1:3) == 'GF_' ) then
      ! Vertical info are not stored in the fixed orography file.
      ! Read part of the info from a T file.
      call find_gfdl_dim(pathaddname)
    else if ( dattyp(1:3) == 'CN_' ) then
      ! Vertical info are not stored in the fixed orography file.
      ! Read part of the info from a T file.
      call find_cnrm_dim(pathaddname)
    else if ( dattyp(1:3) == 'CS_' ) then
      ! Vertical info are not stored in the fixed orography file.
      ! Read part of the info from a T file.
      call find_csiro_dim(pathaddname)
    else if ( dattyp(1:3) == 'MI_' ) then
      ! Vertical info are not stored in the fixed orography file.
      ! Read part of the info from a T file.
      call find_miroc_dim(pathaddname)
    else if ( dattyp(1:2) == 'E5' ) then
      ! Vertical info are stored in the fixed vertinfo file.
      pathaddname = trim(inpglob)//'/ECHAM5/fixed/EH5_OM_1_VERTINFO.nc'
    else if ( dattyp(1:3) == 'MP_' ) then
      ! Vertical info are not stored in the fixed orography file.
      ! Read part of the info from a T file.
      call find_mpiesm_dim(pathaddname)
    else if ( dattyp == 'GFS11' ) then
      pathaddname = trim(inpglob)//'/GFS11/fixed/fixed_orography.nc'
    else if ( dattyp(1:3) == 'EC_' ) then
      pathaddname = trim(inpglob)//'/EC-EARTH/fixed/ecearth.nc'
    else if ( dattyp == 'CCSM3' ) then
      call find_ccsm3_topo(pathaddname)
    else if ( dattyp == 'JRA55' ) then
      pathaddname = trim(inpglob)//'/JRA55/fixed/ll125.006_gp.nc'
    else
      call die('Unknown dattyp in generic 6h NetCDF driver.')
    end if

    istatus = nf90_open(pathaddname,nf90_nowrite,inet1)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error open '//trim(pathaddname))
    istatus = nf90_inq_dimid(inet1,'lon',jdim)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find lon dim')
    istatus = nf90_inquire_dimension(inet1,jdim,len=nlon)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire lon dim')
    istatus = nf90_inq_dimid(inet1,'lat',jdim)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find lat dim')
    istatus = nf90_inquire_dimension(inet1,jdim,len=nlat)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire lat dim')
    istatus = nf90_inq_dimid(inet1,'lev',jdim)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find lev dim')
    istatus = nf90_inquire_dimension(inet1,jdim,len=klev)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire lev dim')

    ! Input layer and pressure interpolated values

    call getmem1d(glat,1,nlat,'mod_gn6hnc:glat')
    call getmem1d(glon,1,nlon,'mod_gn6hnc:glon')
    if ( dattyp(1:2) == 'HA' ) then
      call find_hadgem_ufile(pathaddname)
      istatus = nf90_open(pathaddname,nf90_nowrite,inet2)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error open '//trim(pathaddname))
      istatus = nf90_inq_dimid(inet2,'lon',jdim)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find lon dim')
      istatus = nf90_inquire_dimension(inet2,jdim,len=nulon)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error inquire lon dim')
      call find_hadgem_vfile(pathaddname)
      istatus = nf90_open(pathaddname,nf90_nowrite,inet3)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error open '//trim(pathaddname))
      istatus = nf90_inq_dimid(inet3,'lat',jdim)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find lat dim')
      istatus = nf90_inquire_dimension(inet3,jdim,len=nvlat)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error inquire lat dim')
      call getmem1d(vglat,1,nvlat,'mod_gn6hnc:vglat')
      call getmem1d(uglon,1,nulon,'mod_gn6hnc:uglon')
    end if

    call getmem2d(zsvar,1,nlon,1,nlat,'mod_gn6hnc:zsvar')
    call getmem2d(psvar,1,nlon,1,nlat,'mod_gn6hnc:psvar')

    if ( dattyp /= 'GFS11' .and. dattyp(1:3) /= 'EC_' .and. &
         dattyp(1:2) /= 'E5' .and. dattyp /= 'JRA55' ) then
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
      if ( dattyp == 'GFS11' ) then
        call getmem1d(gltemp,1,nlat,'mod_gn6hnc:gltemp')
        call getmem3d(vwork,1,nlon,1,nlat,1,klev,'mod_gn6hnc:vwork')
      end if
      if ( dattyp == 'JRA55' ) then
        call getmem3d(vwork,1,nlon,1,nlat,1,klev,'mod_gn6hnc:vwork')
      end if
      call getmem3d(b2,1,nlon,1,nlat,1,klev*3,'mod_gn6hnc:b2')
      call getmem3d(d2,1,nlon,1,nlat,1,klev*2,'mod_gn6hnc:d2')
      uvar => d2(:,:,1:klev)
      vvar => d2(:,:,klev+1:2*klev)
      tvar => b2(:,:,1:klev)
      hvar => b2(:,:,klev+1:2*klev)
      qvar => b2(:,:,2*klev+1:3*klev)
      if ( dattyp == 'GFS11' ) then
        istatus = nf90_inq_dimid(inet1,'rhlev',jdim)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error find rhlev dim')
        istatus = nf90_inquire_dimension(inet1,jdim,len=nrhlev)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error inquire rhlev dim')
      else
        call getmem3d(pp3d,1,nlon,1,nlat,1,klev,'mod_gn6hnc:pp3d')
      end if
    end if

    istatus = nf90_inq_varid(inet1,'lat',ivar1)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find lat var')
    istatus = nf90_get_var(inet1,ivar1,glat)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read lat var')
    istatus = nf90_inq_varid(inet1,'lon',ivar1)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error find lon var')
    istatus = nf90_get_var(inet1,ivar1,glon)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read lon var')
    if ( dattyp(1:2) == 'HA' ) then
      istatus = nf90_inq_varid(inet2,'lon',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find lon var')
      istatus = nf90_get_var(inet2,ivar1,uglon)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read lon var')
      istatus = nf90_inq_varid(inet3,'lat',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find lat var')
      istatus = nf90_get_var(inet3,ivar1,vglat)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read lat var')
      istatus = nf90_close(inet2)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error close ua file')
      istatus = nf90_close(inet3)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error close va file')
      call getmem3d(ha_d2_1,1,nulon,1,nlat,1,klev,'mod_gn6hnc:ha_d2_1')
      call getmem3d(ha_d2_2,1,nlon,1,nvlat,1,klev,'mod_gn6hnc:ha_d2_2')
    end if

    icount(1) = nlon
    icount(2) = nlat
    icount(3) = 1
    istart(1) = 1
    istart(2) = 1
    istart(3) = 1

    npl = nipl

    if ( dattyp == 'CAM4N' .or. dattyp == 'CCSMN' .or. &
         dattyp == 'CCSM3' ) then
      istatus = nf90_inq_varid(inet1,'hyam',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find hyam var')
      istatus = nf90_get_var(inet1,ivar1,ak)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read hyam var')
      istatus = nf90_inq_varid(inet1,'hybm',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find hybm var')
      istatus = nf90_get_var(inet1,ivar1,bk)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read hybm var')
      istatus = nf90_inq_varid(inet1,'P0',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find P0 var')
      istatus = nf90_get_var(inet1,ivar1,dp0)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read P0 var')
      p0 = real(dp0)
      istatus = nf90_inq_varid(inet1,'PHIS',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find PHIS var')
      istatus = nf90_get_var(inet1,ivar1,zsvar,istart(1:3),icount(1:3))
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read PHYS var')
      zsvar(:,:) = zsvar(:,:)*real(regrav)
      where (zsvar < 0.0) zsvar = 0.0
    else if ( dattyp == 'JRA55' ) then
      npl = klev ! Data are on pressure levels
      call getmem1d(pplev,1,klev,'mod_gn6hnc:pplev')
      istatus = nf90_inq_varid(inet1,'lev',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find lev var')
      istatus = nf90_get_var(inet1,ivar1,pplev)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read lev var')
      istatus = nf90_inq_varid(inet1,'var6',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find var6 var')
      istatus = nf90_get_var(inet1,ivar1,zsvar,istart(1:3),icount(1:3))
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var6 var')
      zsvar(:,:) = zsvar(:,:)*real(regrav)
      where (zsvar < 0.0) zsvar = 0.0
    else if ( dattyp(1:3) == 'HA_' ) then
      istatus = nf90_inq_varid(inet1,'lev',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find lev var')
      istatus = nf90_get_var(inet1,ivar1,ak)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read lev var')
      istatus = nf90_inq_varid(inet1,'b',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find b var')
      istatus = nf90_get_var(inet1,ivar1,bk)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read b var')
      istatus = nf90_inq_varid(inet1,'orog',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find orog var')
      istatus = nf90_get_var(inet1,ivar1,zsvar,istart(1:3),icount(1:3))
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read orog var')
      ! call getmem2d(pmslvar,1,nlon,1,nlat,'mod_gn6hnc:pmslvar')
    else if ( dattyp(1:2) == 'E5' ) then
      npl = klev ! Data are on pressure levels
      call getmem1d(pplev,1,klev,'mod_gn6hnc:pplev')
      istatus = nf90_inq_varid(inet1,'lev',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find lev var')
      istatus = nf90_get_var(inet1,ivar1,pplev)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read lev var')
      ! Close the T file, get just orography from fixed file.
      istatus = nf90_close(inet1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error close file '//trim(pathaddname))
      ! This one contains just orography.
      pathaddname = trim(inpglob)// &
            '/ECHAM5/fixed/EH5_OM_1_GEOSP_19000101.nc'
      istatus = nf90_open(pathaddname,nf90_nowrite,inet1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error open '//trim(pathaddname))
      istatus = nf90_inq_varid(inet1,'geosp',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find geosporog var')
      istatus = nf90_get_var(inet1,ivar1,zsvar,istart(1:3),icount(1:3))
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read orog var')
      zsvar(:,:) = zsvar(:,:)*real(regrav)
      call getmem2d(glat2,1,nlon,1,nlat,'mod_gn6hnc:glat2')
      call getmem2d(glon2,1,nlon,1,nlat,'mod_gn6hnc:glon2')
      do j = 1 , nlon
        glat2(j,:) = glat(:)
      end do
      do i = 1 , nlat
        glon2(:,i) = glon(:)
      end do
    else if ( dattyp(1:3) == 'MP_' ) then
      istatus = nf90_inq_varid(inet1,'hyam',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find hyam var')
      istatus = nf90_get_var(inet1,ivar1,ak)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read hyam var')
      istatus = nf90_inq_varid(inet1,'hybm',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find b var')
      istatus = nf90_get_var(inet1,ivar1,bk)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read b var')
      ! Close the T file, get just orography from fixed file.
      istatus = nf90_close(inet1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error close file '//trim(pathaddname))
      ! This one contains just orography.
      call find_mpiesm_topo(pathaddname)
      istatus = nf90_open(pathaddname,nf90_nowrite,inet1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error open '//trim(pathaddname))
      istatus = nf90_inq_varid(inet1,'geosp',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find geosporog var')
      istatus = nf90_get_var(inet1,ivar1,zsvar,istart(1:3),icount(1:3))
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read orog var')
      zsvar(:,:) = zsvar(:,:)*real(regrav)
      call getmem2d(glat2,1,nlon,1,nlat,'mod_gn6hnc:glat2')
      call getmem2d(glon2,1,nlon,1,nlat,'mod_gn6hnc:glon2')
      do j = 1 , nlon
        glat2(j,:) = glat(:)
      end do
      do i = 1 , nlat
        glon2(:,i) = glon(:)
      end do
    else if ( dattyp(1:3) == 'CA_' ) then
      istatus = nf90_inq_varid(inet1,'ap',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find ap var')
      istatus = nf90_get_var(inet1,ivar1,ak)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read ap var')
      istatus = nf90_inq_varid(inet1,'b',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find b var')
      istatus = nf90_get_var(inet1,ivar1,bk)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read b var')
      ! Close the T file, get just orography from fixed file.
      istatus = nf90_close(inet1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error close file '//trim(pathaddname))
      ! This one contains just orography.
      call find_canesm_topo(pathaddname)
      istatus = nf90_open(pathaddname,nf90_nowrite,inet1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error open '//trim(pathaddname))
      istatus = nf90_inq_varid(inet1,'orog',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find orog var')
      istatus = nf90_get_var(inet1,ivar1,zsvar,istart(1:3),icount(1:3))
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read orog var')
      call getmem2d(glat2,1,nlon,1,nlat,'mod_gn6hnc:glat2')
      call getmem2d(glon2,1,nlon,1,nlat,'mod_gn6hnc:glon2')
      do j = 1 , nlon
        glat2(j,:) = glat(:)
      end do
      do i = 1 , nlat
        glon2(:,i) = glon(:)
      end do
    else if ( dattyp(1:3) == 'IP_' ) then
      istatus = nf90_inq_varid(inet1,'ap',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find ap var')
      istatus = nf90_get_var(inet1,ivar1,ak)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read ap var')
      istatus = nf90_inq_varid(inet1,'b',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find b var')
      istatus = nf90_get_var(inet1,ivar1,bk)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read b var')
      ! Close the T file, get just orography from fixed file.
      istatus = nf90_close(inet1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error close file '//trim(pathaddname))
      ! This one contains just orography.
      call find_ipsl_topo(pathaddname)
      istatus = nf90_open(pathaddname,nf90_nowrite,inet1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error open '//trim(pathaddname))
      istatus = nf90_inq_varid(inet1,'orog',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find orog var')
      istatus = nf90_get_var(inet1,ivar1,zsvar,istart(1:3),icount(1:3))
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read orog var')
      call getmem2d(glat2,1,nlon,1,nlat,'mod_gn6hnc:glat2')
      call getmem2d(glon2,1,nlon,1,nlat,'mod_gn6hnc:glon2')
      do j = 1 , nlon
        glat2(j,:) = glat(:)
      end do
      do i = 1 , nlat
        glon2(:,i) = glon(:)
      end do
    else if ( dattyp(1:3) == 'GF_' ) then
      istatus = nf90_inq_varid(inet1,'a',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find a var')
      istatus = nf90_get_var(inet1,ivar1,ak)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read a var')
      istatus = nf90_inq_varid(inet1,'b',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find b var')
      istatus = nf90_get_var(inet1,ivar1,bk)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read b var')
      istatus = nf90_inq_varid(inet1,'p0',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find p0 var')
      istatus = nf90_get_var(inet1,ivar1,dp0)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read p0 var')
      p0 = real(dp0)
      ! Close the T file, get just orography from fixed file.
      istatus = nf90_close(inet1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error close file '//trim(pathaddname))
      ! This one contains just orography.
      call find_gfdl_topo(pathaddname)
      istatus = nf90_open(pathaddname,nf90_nowrite,inet1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error open '//trim(pathaddname))
      istatus = nf90_inq_varid(inet1,'orog',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find orog var')
      istatus = nf90_get_var(inet1,ivar1,zsvar,istart(1:3),icount(1:3))
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read orog var')
      call getmem2d(glat2,1,nlon,1,nlat,'mod_gn6hnc:glat2')
      call getmem2d(glon2,1,nlon,1,nlat,'mod_gn6hnc:glon2')
      do j = 1 , nlon
        glat2(j,:) = glat(:)
      end do
      do i = 1 , nlat
        glon2(:,i) = glon(:)
      end do
    else if ( dattyp(1:3) == 'MI_' ) then
      istatus = nf90_inq_varid(inet1,'a',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find a var')
      istatus = nf90_get_var(inet1,ivar1,ak)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read a var')
      istatus = nf90_inq_varid(inet1,'b',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find b var')
      istatus = nf90_get_var(inet1,ivar1,bk)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read b var')
      istatus = nf90_inq_varid(inet1,'p0',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find p0 var')
      istatus = nf90_get_var(inet1,ivar1,dp0)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read p0 var')
      p0 = real(dp0)
      ! Close the T file, get just orography from fixed file.
      istatus = nf90_close(inet1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error close file '//trim(pathaddname))
      ! This one contains just orography.
      call find_miroc_topo(pathaddname)
      istatus = nf90_open(pathaddname,nf90_nowrite,inet1)
      call checkncerr(istatus,__FILE__,__LINE__, &
        'Error open '//trim(pathaddname))
      istatus = nf90_inq_varid(inet1,'orog',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find orog var')
      istatus = nf90_get_var(inet1,ivar1,zsvar,istart(1:3),icount(1:3))
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read orog var')
      call getmem2d(glat2,1,nlon,1,nlat,'mod_gn6hnc:glat2')
      call getmem2d(glon2,1,nlon,1,nlat,'mod_gn6hnc:glon2')
      do j = 1 , nlon
        glat2(j,:) = glat(:)
      end do
      do i = 1 , nlat
        glon2(:,i) = glon(:)
      end do
    else if ( dattyp(1:3) == 'CN_' ) then
      istatus = nf90_inq_varid(inet1,'a',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find a var')
      istatus = nf90_get_var(inet1,ivar1,ak)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read a var')
      istatus = nf90_inq_varid(inet1,'b',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find b var')
      istatus = nf90_get_var(inet1,ivar1,bk)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read b var')
      istatus = nf90_inq_varid(inet1,'p0',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find p0 var')
      istatus = nf90_get_var(inet1,ivar1,dp0)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read p0 var')
      p0 = real(dp0)
      ! Close the T file, get just orography from fixed file.
      istatus = nf90_close(inet1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error close file '//trim(pathaddname))
      ! This one contains just orography.
      call find_cnrm_topo(pathaddname)
      istatus = nf90_open(pathaddname,nf90_nowrite,inet1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error open '//trim(pathaddname))
      istatus = nf90_inq_varid(inet1,'orog',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find orog var')
      istatus = nf90_get_var(inet1,ivar1,zsvar,istart(1:3),icount(1:3))
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read orog var')
      call getmem2d(glat2,1,nlon,1,nlat,'mod_gn6hnc:glat2')
      call getmem2d(glon2,1,nlon,1,nlat,'mod_gn6hnc:glon2')
      do j = 1 , nlon
        glat2(j,:) = glat(:)
      end do
      do i = 1 , nlat
        glon2(:,i) = glon(:)
      end do
    else if ( dattyp(1:3) == 'CS_' ) then
      istatus = nf90_inq_varid(inet1,'a',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find a var')
      istatus = nf90_get_var(inet1,ivar1,ak)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read a var')
      istatus = nf90_inq_varid(inet1,'b',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find b var')
      istatus = nf90_get_var(inet1,ivar1,bk)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read b var')
      istatus = nf90_inq_varid(inet1,'p0',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find p0 var')
      istatus = nf90_get_var(inet1,ivar1,dp0)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read p0 var')
      p0 = real(dp0)
      ! Close the T file, get just orography from fixed file.
      istatus = nf90_close(inet1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error close file '//trim(pathaddname))
      ! This one contains just orography.
      call find_csiro_topo(pathaddname)
      istatus = nf90_open(pathaddname,nf90_nowrite,inet1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error open '//trim(pathaddname))
      istatus = nf90_inq_varid(inet1,'orog',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find orog var')
      istatus = nf90_get_var(inet1,ivar1,zsvar,istart(1:3),icount(1:3))
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read orog var')
      call getmem2d(glat2,1,nlon,1,nlat,'mod_gn6hnc:glat2')
      call getmem2d(glon2,1,nlon,1,nlat,'mod_gn6hnc:glon2')
      do j = 1 , nlon
        glat2(j,:) = glat(:)
      end do
      do i = 1 , nlat
        glon2(:,i) = glon(:)
      end do
    else if ( dattyp == 'GFS11' ) then
      npl = klev ! Data are on pressure levels
      call getmem1d(pplev,1,klev,'mod_gn6hnc:pplev')
      istatus = nf90_inq_varid(inet1,'lev',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find lev var')
      istatus = nf90_get_var(inet1,ivar1,pplev)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read lev var')
      istatus = nf90_inq_varid(inet1,'orog',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find orog var')
      istatus = nf90_get_var(inet1,ivar1,zsvar,istart(1:3),icount(1:3))
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read orog var')
      do j = 1 , nlat
        gltemp(nlat-j+1) = glat(j)
      end do
      glat(:) = gltemp(:)
      call relmem1d(gltemp)
    else if ( dattyp(1:3) == 'EC_' ) then
      npl = klev ! Data are on pressure levels
      call getmem1d(pplev,1,klev,'mod_gn6hnc:pplev')
      istatus = nf90_inq_varid(inet1,'lev',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find lev var')
      istatus = nf90_get_var(inet1,ivar1,pplev)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read lev var')
      istatus = nf90_inq_varid(inet1,'geo',ivar1)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find geo var')
      istatus = nf90_get_var(inet1,ivar1,zsvar,istart(1:3),icount(1:3))
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read geo var')
      ! Transform geopotential to elevation
      zsvar(:,:) = zsvar(:,:)/real(egrav)
      call getmem2d(glat2,1,nlon,1,nlat,'mod_gn6hnc:glat2')
      call getmem2d(glon2,1,nlon,1,nlat,'mod_gn6hnc:glon2')
      do j = 1 , nlon
        glat2(j,:) = glat(:)
      end do
      do i = 1 , nlat
        glon2(:,i) = glon(:)
      end do
    end if

    istatus = nf90_close(inet1)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error close file '//trim(pathaddname))

    call getmem1d(sigmar,1,npl,'mod_gn6hnc:sigmar')
    call getmem3d(b3,1,jx,1,iy,1,npl*3,'mod_gn6hnc:b3')
    call getmem3d(d3,1,jx,1,iy,1,npl*2,'mod_gn6hnc:d3')

    if ( dattyp /= 'GFS11' .and. dattyp(1:3) /= 'EC_' .and. &
         dattyp(1:2) /= 'E5' .and. dattyp /= 'JRA55' ) then
      call getmem3d(b2,1,nlon,1,nlat,1,npl*3,'mod_gn6hnc:b2')
      call getmem3d(d2,1,nlon,1,nlat,1,npl*2,'mod_gn6hnc:d2')
      up => d2(:,:,1:npl)
      vp => d2(:,:,npl+1:2*npl)
      tp => b2(:,:,1:npl)
      hp => b2(:,:,npl+1:2*npl)
      qp => b2(:,:,2*npl+1:3*npl)
      pplev => fplev
    end if

    ! Set up pointers

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
      if ( glon(1) > glon(nulon) ) then
        do i = 1 , nlon
          if ( glon(i) > 180.0 ) then
            glon(i) = glon(i) - 360.0
          end if
        end do
      end if
    else if ( dattyp(1:3) == 'CS_' ) then
      ! CSIRO datasets has different times for PS and vertical variables.
      pstimlen = 1
      call getmem1d(ipstimes,1,1,'mod_gn6hnc:ipstimes')
      ipstimes(1) = 1870010100 ! This set to a "Prehistorical" date
      call setcal(itimes(1), noleap)
      call setcal(ipstimes(1), noleap)
    else if ( dattyp(1:3) == 'MI_' ) then
      ! MIROC5 datasets has different times for PS and vertical variables.
      pstimlen = 1
      call getmem1d(ipstimes,1,1,'mod_gn6hnc:ipstimes')
      ipstimes(1) = 1870010100 ! This set to a "Prehistorical" date
      call setcal(itimes(1), noleap)
      call setcal(ipstimes(1), noleap)
    else if ( dattyp(1:3) == 'GFS' .or. dattyp(1:3) == 'EC_' .or. &
              dattyp(1:3) == 'CN_' .or. dattyp(1:3) == 'MP_' .or. &
              dattyp(1:2) == 'E5' .or. dattyp == 'JRA55' ) then
      call setcal(itimes(1), gregorian)
    else
      call setcal(itimes(1), noleap)
    end if

    if ( dattyp(1:3) == 'EC_' .or. dattyp == 'JRA55') then
      do k = 1 , npl
        sigmar(k) = pplev(npl-k+1)/pplev(1)
      end do
      pss = pplev(1)/1000.0_rkx ! Pa -> cb
    else
      do k = 1 , npl
        sigmar(k) = pplev(k)/pplev(npl)
      end do
      pss = pplev(npl)/10.0_rkx ! mb -> cb
    end if

    write (stdout,*) 'Read in Static fields OK'

  end subroutine headgn6hnc
  !
  !-----------------------------------------------------------------------
  !
  subroutine get_gn6hnc(idate)
    use netcdf
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate

    call readgn6hnc(idate)
    write (stdout,*) 'Read in fields at Date: ', tochar(idate)

    ! JRA55, GFS, EC-EARTH and ECHAM5 grids are already on pressure levels.
    if ( dattyp /= 'GFS11' .and. dattyp(1:3) /= 'EC_' .and. &
         dattyp(1:2) /= 'E5' .and. dattyp /= 'JRA55' ) then

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
           dattyp(1:3) == 'CS_' .or. dattyp(1:3) == 'MI_' ) then
        call top2btm(tvar,nlon,nlat,klev)
        call top2btm(qvar,nlon,nlat,klev)
        call top2btm(uvar,nlon,nlat,klev)
        call top2btm(vvar,nlon,nlat,klev)
        call top2btm(pp3d,nlon,nlat,klev)
        call htsig(tvar,hvar,pp3d,psvar,zsvar,nlon,nlat,klev)
      end if

      if ( dattyp(1:3) == 'MP_' ) then
        ! Calculate HGT on model hybrid sigma levels
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

    if ( dattyp == 'JRA55' ) then
      call top2btm(tvar,nlon,nlat,klev)
      call top2btm(qvar,nlon,nlat,klev)
      call top2btm(uvar,nlon,nlat,klev)
      call top2btm(vvar,nlon,nlat,klev)
      call top2btm(pp3d,nlon,nlat,klev)
      call top2btm(hvar,nlon,nlat,klev)
    end if

    ! Horizontal interpolation on RegCM grid
    if ( dattyp(1:3) /= 'MP_' .and. dattyp(1:3) /= 'CA_' .and. &
         dattyp(1:3) /= 'CN_' .and. dattyp(1:3) /= 'CS_' .and. &
         dattyp(1:3) /= 'GF_' .and. dattyp(1:3) /= 'IP_' .and. &
         dattyp(1:3) /= 'EC_' .and. dattyp(1:3) /= 'MI_' .and. &
         dattyp(1:2) /= 'E5' ) then
      call bilinx(b3,b2,xlon,xlat,glon,glat,nlon,nlat,jx,iy,npl*3)
      call bilinx(d3,d2,dlon,dlat,glon,glat,nlon,nlat,jx,iy,npl*2)
    else ! Gaussian grid
      call cressmcr(b3,b2,xlon,xlat,glon2,glat2,jx,iy,nlon,nlat,npl,3)
      call cressmdt(d3,d2,dlon,dlat,glon2,glat2,jx,iy,nlon,nlat,npl,2)
    end if

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
    call intgtb(pa,za,tlayer,topogm,t3,h3,pss,sigmar,jx,iy,npl)
    call intpsn(ps4,topogm,pa,za,tlayer,ptop,jx,iy)
    call crs2dot(pd4,ps4,jx,iy,i_band)

    ! Recalculate temperature on RegCM orography
    call intv3(ts4,t3,ps4,pss,sigmar,ptop,jx,iy,npl)
    ! Replace it with SST on water points
    call readsst(ts4,idate)

    ! Vertically interpolate on RegCM sigma levels
    call intv1(u4,u3,pd4,sigmah,pss,sigmar,ptop,jx,iy,kz,npl)
    call intv1(v4,v3,pd4,sigmah,pss,sigmar,ptop,jx,iy,kz,npl)
    call intv2(t4,t3,ps4,sigmah,pss,sigmar,ptop,jx,iy,kz,npl)
    call intv1(q4,q3,ps4,sigmah,pss,sigmar,ptop,jx,iy,kz,npl)
  end subroutine get_gn6hnc
  !
  !-----------------------------------------------------------------------
  !
  subroutine readgn6hnc(idate)
    use netcdf
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    integer(ik4) :: istatus
    integer(ik4) :: i , it , itps , j , k , timid
    character(len=256) :: inname

    integer(ik4) :: kkrec
    character(len=64) :: cunit , ccal
    type(rcm_time_interval) :: tdif
    type(rcm_time_and_date) :: pdate
    integer(ik4) :: year , month , day , hour , y1 , y2 , m1
    integer(ik4) :: fyear , fmonth , fday , fhour

    call split_idate(idate,year,month,day,hour)
    !
    ! This is simpler case: just one file for each timestep with
    ! all variables already on pressure levels.
    !
    if ( dattyp == 'GFS11' ) then
      if ( inet(1) > 0 ) then
        istatus = nf90_close(inet(1))
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error close file')
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
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//varname(6))
      do j = 1 , nlat
        psvar(:,nlat-j+1) = vwork(:,j,1)*0.01 ! Go to mb
      end do
      istatus = nf90_get_var(inet(1),ivar(1),vwork)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//varname(1))
      do j = 1 , nlat
        tvar(:,nlat-j+1,:) = vwork(:,j,:)
      end do
      istatus = nf90_get_var(inet(1),ivar(2),vwork)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//varname(2))
      do j = 1 , nlat
        hvar(:,nlat-j+1,:) = vwork(:,j,:)
      end do
      vwork = 0.0
      istatus = nf90_get_var(inet(1),ivar(3),vwork(:,:,klev-nrhlev+1:klev))
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//varname(3))
      do j = 1 , nlat
        qvar(:,nlat-j+1,:) = vwork(:,j,:)*0.01
      end do
      istatus = nf90_get_var(inet(1),ivar(4),vwork)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//varname(4))
      do j = 1 , nlat
        uvar(:,nlat-j+1,:) = vwork(:,j,:)
      end do
      istatus = nf90_get_var(inet(1),ivar(5),vwork)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//varname(5))
      do j = 1 , nlat
        vvar(:,nlat-j+1,:) = vwork(:,j,:)
      end do
      call rh2mxr(tvar,qvar,pplev,nlon,nlat,klev)
    else if ( dattyp == 'JRA55' ) then
      varname => ec5vars
      if ( idate < itimes(1) .or. idate > itimes(timlen) ) then
        do kkrec = 1 , 5
          if ( inet(kkrec) > 0 ) then
            istatus = nf90_close(inet(kkrec))
            call checkncerr(istatus,__FILE__,__LINE__, &
                            'Error close file')
          end if
        end do
        ! monthly files, one for each variable
        do i = 1 , nfiles
          y1 = year
          m1 = month
          if ( jra55vars(i) /= 'XXX' ) then
            write (inname,99006) dattyp,pthsep,y1,pthsep, &
               'anl_p125.',trim(jra55name(i)), '.', &
               y1, m1,'0100_',y1,m1,ndaypm(y1,m1,gregorian),'18.nc'
            pathaddname = trim(inpglob)//pthsep//inname
            istatus = nf90_open(pathaddname,nf90_nowrite,inet(i))
            call checkncerr(istatus,__FILE__,__LINE__, &
                            'Error open '//trim(pathaddname))
            istatus = nf90_inq_varid(inet(i),trim(jra55vars(i)),ivar(i))
            call checkncerr(istatus,__FILE__,__LINE__, &
                            'Error find var '//trim(jra55vars(i)))
            write (stdout,*) inet(i), trim(pathaddname)
          end if
        end do
        istatus = nf90_inq_dimid(inet(1),'time',timid)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error find dim time')
        istatus = nf90_inquire_dimension(inet(1),timid, len=timlen)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error inquire dim time')
        istatus = nf90_inq_varid(inet(1),'time',timid)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error find var time')
        istatus = nf90_get_att(inet(1),timid,'units',cunit)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error read time units')
        istatus = nf90_get_att(inet(1),timid,'calendar',ccal)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error read time calendar')
        call getmem1d(itimes,1,timlen,'mod_gn6hnc:itimes')
        call getmem1d(xtimes,1,timlen,'mod_gn6hnc:xtimes')
        istatus = nf90_get_var(inet(1),timid,xtimes)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error read time')
        do it = 1 , timlen
          itimes(it) = timeval2date(xtimes(it),cunit,ccal)
        end do
      end if
      tdif = idate - itimes(1)
      it = nint(tohours(tdif))/6 + 1
      icount(1) = nlon
      icount(2) = nlat
      icount(3) = klev
      icount(4) = 1
      istart(1) = 1
      istart(2) = 1
      istart(3) = 1
      istart(4) = it
      istatus = nf90_get_var(inet(1),ivar(1),vwork,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//jra55vars(1))
      tvar(:,:,:) = vwork(:,:,:)
      istatus = nf90_get_var(inet(2),ivar(2),vwork,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//jra55vars(2))
      hvar(:,:,:) = vwork(:,:,:)
      istatus = nf90_get_var(inet(4),ivar(4),vwork,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//jra55vars(4))
      uvar(:,:,:) = vwork(:,:,:)
      istatus = nf90_get_var(inet(5),ivar(5),vwork,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//jra55vars(5))
      vvar(:,:,:) = vwork(:,:,:)
      do k = 1, klev
        pp3d(:,:,k) = pplev(k)*0.01 ! Get in hPa
      end do
      ! Less levels for relative humidity !
      icount(3) = 27
      vwork = 0.0_rkx
      istatus = nf90_get_var(inet(3),ivar(3),vwork,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//jra55vars(3))
      qvar(:,:,:) = vwork(:,:,:)*0.01_rkx
      call rh2mxr(tvar,qvar,pplev,nlon,nlat,klev)
    else if ( dattyp(1:2) == 'E5' ) then
      varname => ec5vars
      if ( idate < itimes(1) .or. idate > itimes(timlen) ) then
        do kkrec = 1 , 5
          if ( inet(kkrec) > 0 ) then
            istatus = nf90_close(inet(kkrec))
            call checkncerr(istatus,__FILE__,__LINE__, &
                            'Error close file')
          end if
        end do
        ! yearly files, one for each variable
        y1 = year
        y2 = y1+1
        do i = 1 , nfiles
          if ( ec5vars(i) /= 'XXX' ) then
            write (inname,99005) dattyp(4:5),pthsep, &
               'EH5_OM_'//dattyp(4:5)//'_1_', &
               trim(ec5name(i)), '_', y1, '010100-',y1+1,'010100.nc'
            pathaddname = trim(inpglob)//'/ECHAM5/'//inname
            istatus = nf90_open(pathaddname,nf90_nowrite,inet(i))
            call checkncerr(istatus,__FILE__,__LINE__, &
                            'Error open '//trim(pathaddname))
            istatus = nf90_inq_varid(inet(i),trim(varname(i)),ivar(i))
            call checkncerr(istatus,__FILE__,__LINE__, &
                            'Error find var '//trim(varname(i)))
            write (stdout,*) inet(i), trim(pathaddname)
          end if
        end do
        istatus = nf90_inq_dimid(inet(1),'time',timid)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error find dim time')
        istatus = nf90_inquire_dimension(inet(1),timid, len=timlen)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error inquire dim time')
        istatus = nf90_inq_varid(inet(1),'time',timid)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error find var time')
        istatus = nf90_get_att(inet(1),timid,'units',cunit)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error read time units')
        istatus = nf90_get_att(inet(1),timid,'calendar',ccal)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error read time calendar')
        call getmem1d(itimes,1,timlen,'mod_gn6hnc:itimes')
        call getmem1d(xtimes,1,timlen,'mod_gn6hnc:xtimes')
        istatus = nf90_get_var(inet(1),timid,xtimes)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error read time')
        do it = 1 , timlen
          itimes(it) = timeval2date(xtimes(it),cunit,ccal)
        end do
      end if
      tdif = idate - itimes(1)
      it = nint(tohours(tdif))/6 + 1
      icount(1) = nlon
      icount(2) = nlat
      icount(3) = klev
      icount(4) = 1
      istart(1) = 1
      istart(2) = 1
      istart(3) = 1
      istart(4) = it
      istatus = nf90_get_var(inet(1),ivar(1),tvar,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//varname(1))
      istatus = nf90_get_var(inet(2),ivar(2),hvar,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//varname(2))
      istatus = nf90_get_var(inet(3),ivar(3),qvar,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//varname(3))
      istatus = nf90_get_var(inet(4),ivar(4),uvar,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//varname(4))
      istatus = nf90_get_var(inet(5),ivar(5),vvar,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//varname(5))
      do k = 1, klev
        pp3d(:,:,k) = pplev(k)*0.01 ! Get in hPa
      end do
      call rh2mxr(tvar,qvar,pplev,nlon,nlat,klev)
    ! More difficult. Multiple files per variable and per year
    else if ( dattyp(1:3) == 'EC_' ) then
      if ( idate < itimes(1) .or. idate > itimes(timlen) ) then
        do kkrec = 1 , 5
          if ( inet(kkrec) > 0 ) then
            istatus = nf90_close(inet(kkrec))
            call checkncerr(istatus,__FILE__,__LINE__, &
                            'Error close file')
          end if
        end do
        varname => echvars
        do kkrec = 1 , 5
          if ( .not. date_in_scenario(idate,5,.true.) ) then
            write (inname,99004) 'RF', pthsep, year, pthsep, 'ich1_', &
                  trim(varname(kkrec))//'_', year, '.nc'
          else
            write (inname,99004) ('RCP'//dattyp(4:5)), pthsep, year,  &
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
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error find dim time')
        istatus = nf90_inquire_dimension(inet(1),timid, len=timlen)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error inquire dim time')
        istatus = nf90_inq_varid(inet(1),'time',timid)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error find var time')
        istatus = nf90_get_att(inet(1),timid,'units',cunit)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error read time units')
        istatus = nf90_get_att(inet(1),timid,'calendar',ccal)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error read time calendar')
        call getmem1d(itimes,1,timlen,'mod_gn6hnc:itimes')
        call getmem1d(xtimes,1,timlen,'mod_gn6hnc:xtimes')
        istatus = nf90_get_var(inet(1),timid,xtimes)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error read time')
        do it = 1 , timlen
          itimes(it) = timeval2date(xtimes(it),cunit,ccal)
        end do
      end if
      tdif = idate - itimes(1)
      it = nint(tohours(tdif))/6 + 1
      icount(1) = nlon
      icount(2) = nlat
      icount(3) = klev
      icount(4) = 1
      istart(1) = 1
      istart(2) = 1
      istart(3) = 1
      istart(4) = it
      istatus = nf90_get_var(inet(1),ivar(1),tvar,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//varname(1))
      istatus = nf90_get_var(inet(2),ivar(2),hvar,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//varname(2))
      hvar = hvar/egrav ! Get to m
      where ( hvar < 0.0 )
        hvar = 0.0 ! We do not want negative hgt, don't we?
      end where
      ! This is specific humidity
      istatus = nf90_get_var(inet(3),ivar(3),qvar,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//varname(3))
      istatus = nf90_get_var(inet(4),ivar(4),uvar,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//varname(4))
      istatus = nf90_get_var(inet(5),ivar(5),vvar,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//varname(5))
      do k = 1, klev
        pp3d(:,:,k) = pplev(k)*0.01 ! Get in hPa
      end do
    else
      ! Even more difficult. Each data type has its own quirks.
      tdif = rcm_time_interval(180,uhrs)
      ! HadGEM dataset has PS files with different times.
      if ( dattyp(1:3) == 'HA_' ) then
        if ( idate < ipstimes(1) .or. idate > ipstimes(pstimlen) ) then
          if ( inet(6) > 0 ) then
            istatus = nf90_close(inet(6))
            call checkncerr(istatus,__FILE__,__LINE__, &
                            'Error close file')
          end if
          call find_hadgem_file(pathaddname,'ps',idate)
          istatus = nf90_open(pathaddname,nf90_nowrite,inet(6))
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error open '//trim(pathaddname))
          istatus = nf90_inq_dimid(inet(6),'time',timid)
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error find dim time')
          istatus = nf90_inquire_dimension(inet(6),timid, len=pstimlen)
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error inquire dim time')
          istatus = nf90_inq_varid(inet(6),'time',timid)
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error find var time')
          istatus = nf90_get_att(inet(6),timid,'units',cunit)
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error read time units')
          istatus = nf90_get_att(inet(6),timid,'calendar',ccal)
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error read time calendar')
          call getmem1d(ipstimes,1,pstimlen,'mod_gn6hnc:ipstimes')
          call getmem1d(xtimes,1,pstimlen,'mod_gn6hnc:xtimes')
          istatus = nf90_get_var(inet(6),timid,xtimes)
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error read time')
          do it = 1 , pstimlen
            ipstimes(it) = timeval2date(xtimes(it),cunit,ccal)
          end do
          istatus = nf90_inq_varid(inet(6), trim(havars(6)), ivar(6))
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error find var '//trim(havars(6)))
          write (stdout,*) 'Open file ', trim(pathaddname)
        end if
      end if
      ! CSIRO dataset has PS files with different times.
      if ( dattyp(1:3) == 'CS_' ) then
        if ( idate < ipstimes(1) .or. idate > ipstimes(pstimlen) ) then
          if ( inet(6) > 0 ) then
            istatus = nf90_close(inet(6))
            call checkncerr(istatus,__FILE__,__LINE__, &
                            'Error close file')
          end if
          call find_csiro_file(pathaddname,csirvars(6),idate)
          istatus = nf90_open(pathaddname,nf90_nowrite,inet(6))
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error open '//trim(pathaddname))
          istatus = nf90_inq_dimid(inet(6),'time',timid)
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error find dim time')
          istatus = nf90_inquire_dimension(inet(6),timid, len=pstimlen)
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error inquire dim time')
          istatus = nf90_inq_varid(inet(6),'time',timid)
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error find var time')
          istatus = nf90_get_att(inet(6),timid,'units',cunit)
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error read time units')
          istatus = nf90_get_att(inet(6),timid,'calendar',ccal)
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error read time calendar')
          call getmem1d(ipstimes,1,pstimlen,'mod_gn6hnc:ipstimes')
          call getmem1d(xtimes,1,pstimlen,'mod_gn6hnc:xtimes')
          istatus = nf90_get_var(inet(6),timid,xtimes)
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error read time')
          do it = 1 , pstimlen
            ipstimes(it) = timeval2date(xtimes(it),cunit,ccal)
          end do
          istatus = nf90_inq_varid(inet(6), trim(csirvars(6)), ivar(6))
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error find var '//trim(csirvars(6)))
          write (stdout,*) 'Open file ', trim(pathaddname)
        end if
      end if
      ! MIROC dataset has PS files with different times
      if ( dattyp(1:3) == 'MI_' ) then
        if ( idate < ipstimes(1) .or. idate > ipstimes(pstimlen) ) then
          if ( inet(6) > 0 ) then
            istatus = nf90_close(inet(6))
            call checkncerr(istatus,__FILE__,__LINE__, &
                            'Error close file')
          end if
          call find_miroc_file(pathaddname,'ps',idate)
          istatus = nf90_open(pathaddname,nf90_nowrite,inet(6))
          call checkncerr(istatus,__FILE__,__LINE__, &
            'Error open '//trim(pathaddname))
          istatus = nf90_inq_dimid(inet(6),'time',timid)
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error find dim time')
          istatus = nf90_inquire_dimension(inet(6),timid, len=pstimlen)
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error inquire dim time')
          istatus = nf90_inq_varid(inet(6),'time',timid)
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error find var time')
          istatus = nf90_get_att(inet(6),timid,'units',cunit)
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error read time units')
          istatus = nf90_get_att(inet(6),timid,'calendar',ccal)
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error read time calendar')
          call getmem1d(ipstimes,1,pstimlen,'mod_gn6hnc:ipstimes')
          call getmem1d(xtimes,1,pstimlen,'mod_gn6hnc:xtimes')
          istatus = nf90_get_var(inet(6),timid,xtimes)
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error read time')
          do it = 1 , pstimlen
            ipstimes(it) = timeval2date(xtimes(it),cunit,ccal)
          end do
          istatus = nf90_inq_varid(inet(6), trim(mirocvars(6)), ivar(6))
          call checkncerr(istatus,__FILE__,__LINE__, &
             'Error find var '//trim(mirocvars(6)))
          write (stdout,*) inet(6), trim(pathaddname)
        end if
      end if

      if ( idate < itimes(1) .or. idate > itimes(timlen) ) then
        if (inet(1) > 0) then
          if ( dattyp == 'CAM4N' ) then
            istatus = nf90_close(inet(1))
            call checkncerr(istatus,__FILE__,__LINE__, &
                            'Error close file')
            filedate = filedate + tdif
          else if ( dattyp == 'CCSM3' ) then
            istatus = nf90_close(inet(1))
            call checkncerr(istatus,__FILE__,__LINE__, &
                            'Error close file')
          else
            if ( dattyp(1:3) == 'HA_' .or. &
                 dattyp(1:3) == 'CS_' .or. &
                 dattyp(1:3) == 'MI_' .or. &
                 dattyp(1:3) == 'IP_' ) then
              do i = 1 , nfiles-1
                if ( havars(i) /= 'XXX' ) then
                  istatus = nf90_close(inet(i))
                  call checkncerr(istatus,__FILE__,__LINE__, &
                                  'Error close file')
                end if
              end do
            else
              do i = 1 , nfiles
                if ( havars(i) /= 'XXX' ) then
                  istatus = nf90_close(inet(i))
                  call checkncerr(istatus,__FILE__,__LINE__, &
                                  'Error close file')
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
        else if ( dattyp == 'CCSM3' ) then
          call find_ccsm3_file(pathaddname,year,month,day,hour)
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
          do i = 1 , nfiles-1
            if ( havars(i) /= 'XXX' ) then
              call find_hadgem_file(pathaddname,havars(i),idate)
              istatus = nf90_open(pathaddname,nf90_nowrite,inet(i))
              call checkncerr(istatus,__FILE__,__LINE__, &
                              'Error open '//trim(pathaddname))
              write (stdout,*) 'Open file ', trim(pathaddname)
            end if
          end do
          varname => havars
        else if ( dattyp(1:3) == 'CA_' ) then
          ! yearly files, one for each variable
          do i = 1 , nfiles
            if ( cavars(i) /= 'XXX' ) then
              call find_canesm_file(pathaddname,cavars(i),idate)
              istatus = nf90_open(pathaddname,nf90_nowrite,inet(i))
              call checkncerr(istatus,__FILE__,__LINE__, &
                              'Error open '//trim(pathaddname))
            end if
            write (stdout,*) 'Open file ', trim(pathaddname)
          end do
          varname => cavars
        else if ( dattyp(1:3) == 'IP_' ) then
          do i = 1 , nfiles
            if ( ipvars(i) /= 'XXX' ) then
              call find_ipsl_file(pathaddname,ipvars(i),idate)
              istatus = nf90_open(pathaddname,nf90_nowrite,inet(i))
              call checkncerr(istatus,__FILE__,__LINE__, &
                              'Error open '//trim(pathaddname))
              write (stdout,*) 'Open file ', trim(pathaddname)
            end if
          end do
          varname => ipvars
        else if ( dattyp(1:3) == 'GF_' ) then
          ! 5 yearly files, one for each variable
          do i = 1 , nfiles
            if ( gfdlvars(i) /= 'XXX' ) then
              call find_gfdl_file(pathaddname,gfdlvars(i),idate)
              istatus = nf90_open(pathaddname,nf90_nowrite,inet(i))
              call checkncerr(istatus,__FILE__,__LINE__, &
                              'Error open '//trim(pathaddname))
              write (stdout,*) 'Open file ', trim(pathaddname)
            end if
          end do
          varname => gfdlvars
        else if ( dattyp(1:3) == 'CN_' ) then
          do i = 1 , nfiles
            if ( cnrmvars(i) /= 'XXX' ) then
              call find_cnrm_file(pathaddname,cnrmvars(i),idate)
              istatus = nf90_open(pathaddname,nf90_nowrite,inet(i))
              call checkncerr(istatus,__FILE__,__LINE__, &
                              'Error open '//trim(pathaddname))
              write (stdout,*) 'Open file ', trim(pathaddname)
            end if
          end do
          varname => cnrmvars
        else if ( dattyp(1:3) == 'CS_' ) then
          ! yearly files, one for each variable
          do i = 1 , nfiles-1
            if ( csirvars(i) /= 'XXX' ) then
              call find_csiro_file(pathaddname,csirvars(i),idate)
              istatus = nf90_open(pathaddname,nf90_nowrite,inet(i))
              call checkncerr(istatus,__FILE__,__LINE__, &
                              'Error open '//trim(pathaddname))
              write (stdout,*) 'Open file ', trim(pathaddname)
            end if
          end do
          varname => csirvars
        else if ( dattyp(1:3) == 'MI_' ) then
          do i = 1 , nfiles-1
            if ( mirocvars(i) /= 'XXX' ) then
              call find_miroc_file(pathaddname,csirvars(i),idate)
              istatus = nf90_open(pathaddname,nf90_nowrite,inet(i))
              call checkncerr(istatus,__FILE__,__LINE__, &
                'Error open '//trim(pathaddname))
              write (stdout,*) 'Open file ', trim(pathaddname)
            end if
          end do
          varname => mirocvars
        else if ( dattyp(1:3) == 'MP_' ) then
          do i = 1 , nfiles
            if ( mpievars(i) /= 'XXX' ) then
              call find_mpiesm_file(pathaddname,mpievars(i),idate)
              istatus = nf90_open(pathaddname,nf90_nowrite,inet(i))
              call checkncerr(istatus,__FILE__,__LINE__, &
                              'Error open '//trim(pathaddname))
              write (stdout,*) 'Open file ', trim(pathaddname)
            end if
          end do
          varname => mpievars
        end if
      end if

      istatus = nf90_inq_dimid(inet(1),'time',timid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find dim time')
      istatus = nf90_inquire_dimension(inet(1),timid, len=timlen)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error inquire dim time')
      istatus = nf90_inq_varid(inet(1),'time',timid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find var time')
      istatus = nf90_get_att(inet(1),timid,'units',cunit)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read time units')
      istatus = nf90_get_att(inet(1),timid,'calendar',ccal)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read time calendar')
      call getmem1d(itimes,1,timlen,'mod_gn6hnc:itimes')
      call getmem1d(xtimes,1,timlen,'mod_gn6hnc:xtimes')
      istatus = nf90_get_var(inet(1),timid,xtimes)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read time')
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
      it = nint(tohours(tdif))/6 + 1

      icount(1) = nlon
      icount(2) = nlat
      icount(3) = 1
      istart(1) = 1
      istart(2) = 1
      istart(3) = it
      if ( dattyp(1:3) == 'HA_' ) then
        tdif = idate - ipstimes(1)
        itps = nint(tohours(tdif))/6 + 1
        istart(3) = itps
        !istatus = nf90_get_var(inet(6),ivar(6),pmslvar,istart(1:3),icount(1:3))
        !call checkncerr(istatus,__FILE__,__LINE__, &
        !                'Error read var '//varname(6))
        !pmslvar(:,:) = pmslvar(:,:)*0.01_rkx
        istatus = nf90_get_var(inet(6),ivar(6),psvar,istart(1:3),icount(1:3))
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error read var '//varname(6))
        psvar(:,:) = psvar(:,:)*0.01_rkx
      else if ( dattyp(1:3) == 'CS_' .or. &
                dattyp(1:3) == 'MI_' ) then
        tdif = idate - ipstimes(1)
        itps = nint(tohours(tdif))/6 + 1
        istart(3) = itps
        istatus = nf90_get_var(inet(6),ivar(6),psvar,istart(1:3),icount(1:3))
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error read var '//varname(6))
      else
        istatus = nf90_get_var(inet(6),ivar(6),psvar,istart(1:3),icount(1:3))
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error read var '//varname(6))
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
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//varname(1))
      if ( dattyp == 'CAM4N' .or. dattyp == 'CCSMN' .or. &
           dattyp == 'CCSM3' ) then
        ! We have geopotential HGT in m, on hybrid sigma pressure levels
        istatus = nf90_get_var(inet(2),ivar(2),hvar,istart,icount)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error read var '//varname(2))
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
        ! Data are on sigma Z levels
        do k = 1 , klev
          hvar(:,:,k) = ak(k) + bk(k)*zsvar(:,:)
        end do
        ! If we have MSLP instead of PS
        ! call mslp2ps(hvar,tvar,pmslvar,zsvar,psvar,nlon,nlat,klev)
        call psig(tvar,hvar,pp3d,psvar,zsvar,nlon,nlat,klev)
      else if ( dattyp(1:3) == 'CA_' .or. dattyp(1:3) == 'IP_' ) then
        ! Data on sigma P levels, the factor for ak is ipotized.
        ! Units in file is Pa, but calculations suggest mPa
        do k = 1, klev
          pp3d(:,:,k) = ak(k)*0.001 + bk(k)*psvar(:,:)
        end do
        psvar(:,:) = psvar(:,:)*0.01
        pp3d(:,:,:) = pp3d(:,:,:)*0.01
      else if ( dattyp(1:3) == 'MP_' ) then
        ! Data on hybrid sigma P levels
        do k = 1, klev
          pp3d(:,:,k) = ak(k) + bk(k)*psvar(:,:)
        end do
        psvar(:,:) = psvar(:,:)*0.01
        pp3d(:,:,:) = pp3d(:,:,:)*0.01
      else if ( dattyp(1:3) == 'GF_' .or. dattyp(1:3) == 'CN_' .or. &
                dattyp(1:3) == 'CS_' .or. dattyp(1:3) == 'MI_' ) then
        do k = 1, klev
          pp3d(:,:,k) = ak(k)*p0 + bk(k)*psvar(:,:)
        end do
        ! Go to hPa
        psvar(:,:) = psvar(:,:)*0.01
        pp3d(:,:,:) = pp3d(:,:,:)*0.01
      end if
      istatus = nf90_get_var(inet(3),ivar(3),qvar,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//varname(3))

      if ( dattyp(1:3) == 'HA_' .or. dattyp(1:3) == 'CA_' .or. &
           dattyp(1:3) == 'IP_' .or. dattyp(1:3) == 'GF_' .or. &
           dattyp(1:3) == 'CN_' .or. dattyp(1:3) == 'CS_' .or. &
           dattyp(1:3) == 'MI_' .or. dattyp(1:3) == 'MP' ) then
        call sph2mxr(qvar,nlon,nlat,klev)
      end if

      if ( dattyp(1:3) == 'HA_' ) then
        icount(1) = nulon
        icount(2) = nlat
        istatus = nf90_get_var(inet(4),ivar(4),ha_d2_1,istart,icount)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error read var '//varname(4))
        do k = 1 , klev
          do j = 1 , nlat
            uvar(1,j,k) = ha_d2_1(1,j,k)
          end do
        end do
        do k = 1 , klev
          do j = 1 , nlat
            do i = 1 , nulon-1
              uvar(i+1,j,k) = 0.5*(ha_d2_1(i,j,k) + ha_d2_1(i+1,j,k))
            end do
          end do
        end do
        do k = 1 , klev
          do j = 1 , nlat
            uvar(nlon,j,k) = ha_d2_1(nulon,j,k)
          end do
        end do
        icount(1) = nlon
        icount(2) = nvlat
        istatus = nf90_get_var(inet(5),ivar(5),ha_d2_2,istart,icount)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error read var '//varname(5))
        do k = 1 , klev
          do i = 1 , nlon
            vvar(i,1,k) = ha_d2_2(i,1,k)
          end do
        end do
        do k = 1 , klev
          do j = 1 , nvlat-1
            do i = 1 , nlon
              vvar(i,j+1,k) = 0.5*(ha_d2_2(i,j,k) + ha_d2_2(i,j+1,k))
            end do
          end do
        end do
        do k = 1 , klev
          do i = 1 , nlon
            vvar(i,nlat,k) = ha_d2_2(i,nvlat,k)
          end do
        end do
      else
        istatus = nf90_get_var(inet(4),ivar(4),uvar,istart,icount)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error read var '//varname(4))
        istatus = nf90_get_var(inet(5),ivar(5),vvar,istart,icount)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error read var '//varname(5))
      end if

    end if ! Other data types not the GFS11 or E_ICH1

99001   format (i0.4,a,a,i0.4,i0.2,i0.2,a,i0.2,a)
99002   format (a,i0.4,'-',i0.2,'-',i0.2,'-',i0.5,'.nc')
99003   format (i0.4,'/','ccsm.',a,a,'.',i0.4,'.nc')
99004   format (a,a,i0.4,a,a,a,i0.4,a)
99005   format (a,a,a,a,a,i0.4,a,i0.4,a)
99006   format (a,a,i0.4,a,a,a,a,i0.4,i0.2,a,i0.4,i0.2,i0.2,a)

  end subroutine readgn6hnc

end module mod_gn6hnc

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

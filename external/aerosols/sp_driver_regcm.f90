!>
!! @par Copyright
!!
!! @brief Program to drive the MACv2-SP module "mo_simple_plumes"
!!
!! @author Karsten Peters, Stephanie Fiedler and Bjorn Stevens MPI-M, Hamburg (2015-12-05)
!!
!! @change-log:
!!          - 2016-12-05: beta release (KP, SF and BS, MPI-Met)
!!          - 2016-10-15: modified treatment of orography (BS)
!!          - 2016-11-09: changed lambda from  517.46442 to 550 nm 517.46442 nm used to make plots in paper (BS)
!!          - 2016-11-09: changed nlev from 48 to 80 and stretched orography, 48 used to make plots in paper (BS)
!!          - 2016-11-09: capability to output 2D fields (SF)
!!          - 2016-11-09: increasing level thickness with height, constant dz=500m used in plots for paper (KP)
!!
!! $ID: n/a$
!!
!!
!! @par Copyright
!!
!
PROGRAM MACv2SP

  USE mo_simple_plumes, ONLY: sp_aop_profile
  USE netcdf

  IMPLICIT NONE

  INTEGER, PARAMETER :: nmon = 12

  REAL, PARAMETER    :: mons(nmon) = (/1., 2., 3., 4., 5., 6., 7., 8., 9.,10.,11.,12./)
  REAL, PARAMETER    :: lambda     = 550.0

  INTEGER ::           &
       nlat          , & !< number of latitudes, defined by orography file
       nlon          , & !< number of longitudes, defined by orography file
       nlev          , & !< number of vertical levels
       iret          , & !< netCDF reading return variable
       ncid          , & !< netCDF file ID
       latID         , & !< pointer to latitude dimension in netCDF file
       lonID         , & !< pointer to longitude dimension in netCDF file
       VarID         , & !< pointer to generic dimension in netCDF file
       levID         , & !< pointer to level dimension in netCDF file
       monID         , & !< pointer to month dimension in netCDF file
       var_t_ID      , & !< pointer to time variable in netCDF file
       var_lat_ID    , & !< pointer to latitude variable in netCDF file
       var_lon_ID    , & !< pointer to longitude variable in netCDF file
       var_z_ID      , & !< pointer to height variable in netCDF file
       var_aod_ID    , & !< pointer to AOD variable in netCDF file
       var_asy_ID    , & !< pointer to asymmetry parameter variable in netCDF file
       var_ssa_ID    , & !< pointer to single-scattering albedo variable in netCDF file
       var_dNovrN_ID , & !< pointer to change in droplet number variable in netCDF file
       var_aod_ID_2d , & !< pointer to 2D AOD variable in netCDF file
       var_asy_ID_2d , & !< pointer to 2D asymmetry parameter variable in netCDF file
       var_ssa_ID_2d , & !< pointer to 2D single-scattering albedo variable in netCDF file
       iyear         , & !< index for year loop
       k             , & !< index for height loop
       i             , & !< index for latitude loop
       j             , & !< index for longitude loop
       n             , & !< index for bin loop
       m             , & !< index for month loop
       ifile             !< output file index

  REAL ::                    &
       year_fr

  REAL, ALLOCATABLE ::       &
       lat(:,:)            , & !< array of latitudes allocated based on nlat
       lon(:,:)            , & !< array of longitudes allocated based on nlon
       oro(:,:)            , & !< array of orographic heights (lon,lat)
       aod(:,:,:,:)        , & !< array for output aerosol optical depth (lon, lat, level, time)
       ssa(:,:,:,:)        , & !< array for output aerosol single scattering albedo (lon, lat, level, time)
       asy(:,:,:,:)        , & !< array for output aerosol asymmetry parameter (lon, lat, level, time)
       aod_int(:,:,:)      , & !< array for output 2D aerosol optical depth (lon, lat, time)
       ssa_int(:,:,:)      , & !< array for output 2D aerosol single scattering albedo (lon, lat, time)
       asy_int(:,:,:)      , & !< array for output 2D aerosol asymmetry parameter (lon, lat, time)
       dNovrN(:,:,:)       , & !< array for output of change in droplet number (lon, lat, level, time)
       aod_prof(:,:)       , & !< aerosol optical depth profiles by colums (ncolumns,level)
       ssa_prof(:,:)       , & !< aerosol single scattering albedo profiles by colums (ncolumns,level)
       asy_prof(:,:)       , & !< aerosol asymmetry parameter by colums (ncolumns,level)
       dNovrN_prof(:)      , & !< droplet number parameter by colums (ncolumns,level)
       z(:,:)            , & !< heights by colums (ncolumns,level)
       dz(:,:)           , & !< layer thicknesses (ncolumns,level)
       col_lat(:)          , & !< latitudes of columns (ncolumns)
       col_lon(:)          , & !< longitude of columns (ncolumns)
       col_oro(:)              !< longitude of columns (ncolumns)
  !
  ! read in the orography data and use this to define the input (lat/lon) grid
  !
  character(len=256) :: prgname , orography
  call get_command_argument(0,value=prgname)
  call get_command_argument(1,value=orography)

  iret = nf90_open(orography, NF90_NOWRITE, ncid)
  IF (iret /= NF90_NOERR) STOP 'Error in opening orography file'
  iret = nf90_inq_dimid(ncid, "iy", VarID)
  iret = nf90_inquire_dimension(ncid, VarID, len = nlat)
  iret = nf90_inq_dimid(ncid, "jx", VarID)
  iret = nf90_inquire_dimension(ncid, VarID, len = nlon)
  iret = nf90_inq_dimid(ncid, "kz", VarID)
  iret = nf90_inquire_dimension(ncid, VarID, len = nlev)

  ALLOCATE (lat(nlon,nlat))
  ALLOCATE (lon(nlon,nlat))
  ALLOCATE (oro(nlon,nlat))
  ALLOCATE (aod(nlon,nlat,nlev,nmon))
  ALLOCATE (ssa(nlon,nlat,nlev,nmon))
  ALLOCATE (asy(nlon,nlat,nlev,nmon))
  ALLOCATE (dNovrN(nlon,nlat,nmon))
  ALLOCATE (aod_int(nlon,nlat,nmon))
  ALLOCATE (ssa_int(nlon,nlat,nmon))
  ALLOCATE (asy_int(nlon,nlat,nmon))
  ALLOCATE (z(nlon*nlat,nlev))
  ALLOCATE (dz(nlon*nlat,nlev))
  ALLOCATE (col_lat(nlon*nlat))
  ALLOCATE (col_lon(nlon*nlat))
  ALLOCATE (col_oro(nlon*nlat))
  ALLOCATE (dNovrN_prof(nlon*nlat))
  ALLOCATE (aod_prof(nlon*nlat,nlev))
  ALLOCATE (ssa_prof(nlon*nlat,nlev))
  ALLOCATE (asy_prof(nlon*nlat,nlev))

  iret = nf90_inq_varid(ncid, "xlat", VarID)
  iret = nf90_get_var(ncid, VarID, lat  , start=(/1,1/)  ,count=(/nlon,nlat/))
  IF (iret /= NF90_NOERR) STOP 'Error in reading latitudes'
  iret = nf90_inq_varid(ncid, "xlon", VarID)
  iret = nf90_get_var(ncid, VarID, lon  , start=(/1,1/)  ,count=(/nlon,nlat/))
  IF (iret /= NF90_NOERR) STOP 'Error in reading longitudes'
  iret = nf90_inq_varid(ncid, "topo", VarID)
  iret = nf90_get_var(ncid, VarID, oro, start=(/1,1/),count=(/nlon,nlat/))
  IF (iret /= NF90_NOERR) STOP 'Error in reading orographic height'
  iret = nf90_close  (ncid)
  !
  ! define the height array
  !
  z(:,1)  = 0.
  DO k = 2,nlev
      dz(:,k) = 50. * exp(0.03*(k-1))
      z(:,k) = z(:,k-1) + dz(:,k)
  END DO
  n = 1
  do i = 1 , nlat
    do j = 1 , nlon
      col_lat(n) = lat(j,i)
      col_lon(n) = lon(j,i)
      col_oro(n) = oro(j,i)
      n = n + 1
     end do
  end do
  !
  ! loop over annual cycle for given year, or year rangeand calculate aerosol optical properties at each lat and
  ! lon point
  !
  DO iyear = 2005,2005
    DO m = 1,nmon
      year_fr = iyear + (m-0.5) / nmon

      !
      ! here we process ncolumns = nlon in one call to sp_aop_profile, output is written to 2D arrays and then
      ! transfered to 4D output arrays to avoid striding through memory in call to sp_aop_profile
      !
      CALL sp_aop_profile ( &
           'MACv2.0-SP_v1.nc','MACv2.0-SP_v1.nc',nlev,nlon*nlat,lambda, &
           col_oro,col_lon,col_lat,year_fr,z,dz,dNovrN_prof,aod_prof,   &
           ssa_prof,asy_prof )
      DO k=1,nlev
        n = 1
        DO i = 1,nlat
          DO j = 1,nlon
            aod(j,i,k,m) = aod_prof(n,k)
            ssa(j,i,k,m) = ssa_prof(n,k)
            asy(j,i,k,m) = asy_prof(n,k)
            n = n + 1
          END DO
        END DO
      END DO
      n = 1
      DO i = 1,nlat
        DO j = 1,nlon
          dNovrN(j,i,m) = dNovrN_prof(n)
          n = n + 1
        END DO
      END DO
      !
      ! 2D vertically integrated values
      !
      DO i = 1,nlat
        DO j = 1,nlon
          aod_int(j,i,m) = SUM(aod(j,i,:,m))
          ssa_int(j,i,m) = SUM(aod(j,i,:,m)*ssa(j,i,:,m))/aod_int(j,i,m)
          asy_int(j,i,m) = SUM(aod(j,i,:,m)*ssa(j,i,:,m)*asy(j,i,:,m))/(ssa_int(j,i,m)*aod_int(j,i,m))
        END DO
      END DO
    END DO
    !
    ! write output as a netcdf file
    !
    IF (iyear == 2005) THEN
      iret = NF90_NOERR
      m = index(orography,'/',back=.true.)+1
      n = len_trim(orography)
      iret = iret + nf90_create("./MACv2-SP_"//orography(m:n), NF90_CLOBBER, ncid)
      iret = iret + nf90_def_dim(ncid, 'lat'   ,nlat , latID)
      iret = iret + nf90_def_dim(ncid, 'lon'   ,nlon , lonID)
      iret = iret + nf90_def_dim(ncid, 'time'  ,nmon , monID)
      iret = iret + nf90_def_dim(ncid, 'z'     ,nlev , levID)
      IF (iret /= 6*NF90_NOERR) STOP 'Error in Creating File Dimensions'
      !
      iret = NF90_NOERR
      iret = iret + nf90_def_var(ncid, 'time'         , NF90_FLOAT, monID, var_t_ID)
      iret = iret + nf90_def_var(ncid, 'lat'          , NF90_FLOAT, (/lonID,latID/), var_lat_ID)
      iret = iret + nf90_def_var(ncid, 'lon'          , NF90_FLOAT, (/lonID,latID/), var_lon_ID)
      iret = iret + nf90_def_var(ncid, 'z'            , NF90_FLOAT, levID, var_z_ID)
      iret = iret + nf90_def_var(ncid, 'aod'          , NF90_FLOAT, (/lonID,latID,levID,monID/), var_aod_ID)
      iret = iret + nf90_def_var(ncid, 'ssa'          , NF90_FLOAT, (/lonID,latID,levID,monID/), var_ssa_ID)
      iret = iret + nf90_def_var(ncid, 'asy'          , NF90_FLOAT, (/lonID,latID,levID,monID/), var_asy_ID)
      iret = iret + nf90_def_var(ncid, 'dNovrN'       , NF90_FLOAT, (/lonID,latID,monID/), var_dNovrN_ID)
      iret = iret + nf90_def_var(ncid, 'aod_2D'       , NF90_FLOAT, (/lonID,latID,monID/), var_aod_ID_2d)
      iret = iret + nf90_def_var(ncid, 'ssa_2D'       , NF90_FLOAT, (/lonID,latID,monID/), var_ssa_ID_2d)
      iret = iret + nf90_def_var(ncid, 'asy_2D'       , NF90_FLOAT, (/lonID,latID,monID/), var_asy_ID_2d)

      iret = iret + nf90_put_att(ncid, var_t_ID       , "long_name", "month of year")
      iret = iret + nf90_put_att(ncid, var_t_ID       , "units"    , "month")
      iret = iret + nf90_put_att(ncid, var_lat_ID     , "long_name", "latitude")
      iret = iret + nf90_put_att(ncid, var_lat_ID     , "units"    , "degrees_north")
      iret = iret + nf90_put_att(ncid, var_lon_ID     , "long_name", "longitude")
      iret = iret + nf90_put_att(ncid, var_lon_ID     , "units"    , "degrees_east")
      iret = iret + nf90_put_att(ncid, var_z_ID       , "long_name", "height above sea-level")
      iret = iret + nf90_put_att(ncid, var_z_ID       , "units"    , "m")
      iret = iret + nf90_put_att(ncid, var_z_ID       , "coordinates"    , "lon lat")
      iret = iret + nf90_put_att(ncid, var_aod_ID     , "long_name", "aerosol optical depth")
      iret = iret + nf90_put_att(ncid, var_aod_ID     , "units"    ," ")
      iret = iret + nf90_put_att(ncid, var_aod_ID     , "coordinates"    , "lon lat")
      iret = iret + nf90_put_att(ncid, var_ssa_ID     , "long_name", "aerosol single scattering albedo")
      iret = iret + nf90_put_att(ncid, var_ssa_ID     , "units"    ," ")
      iret = iret + nf90_put_att(ncid, var_ssa_ID     , "coordinates"    , "lon lat")
      iret = iret + nf90_put_att(ncid, var_asy_ID     , "long_name", "aerosol asymmetry parameter")
      iret = iret + nf90_put_att(ncid, var_asy_ID     , "units"    ," ")
      iret = iret + nf90_put_att(ncid, var_asy_ID     , "coordinates"    , "lon lat")
      iret = iret + nf90_put_att(ncid, var_dNovrN_ID  , "long_name", "normalized change in drop number")
      iret = iret + nf90_put_att(ncid, var_dNovrN_ID  , "units"    ," ")
      iret = iret + nf90_put_att(ncid, var_dNovrN_ID  , "coordinates"    , "lon lat")
      iret = iret + nf90_put_att(ncid, var_aod_ID_2d  , "long_name", "2D aerosol optical depth")
      iret = iret + nf90_put_att(ncid, var_aod_ID_2d  , "units"    ," ")
      iret = iret + nf90_put_att(ncid, var_aod_ID_2d  , "coordinates"    , "lon lat")
      iret = iret + nf90_put_att(ncid, var_ssa_ID_2d  , "long_name", "2D aerosol single scattering albedo")
      iret = iret + nf90_put_att(ncid, var_ssa_ID_2d  , "units"    ," ")
      iret = iret + nf90_put_att(ncid, var_ssa_ID_2d  , "coordinates"    , "lon lat")
      iret = iret + nf90_put_att(ncid, var_asy_ID_2d  , "long_name", "2D aerosol asymmetry parameter")
      iret = iret + nf90_put_att(ncid, var_asy_ID_2d  , "units"    ," ")
      iret = iret + nf90_put_att(ncid, var_asy_ID_2d  , "coordinates"    , "lon lat")
      iret = iret + nf90_enddef(ncid)
      IF (iret /= 26*NF90_NOERR) STOP 'Error in creating file variables'
      !
      iret = NF90_NOERR
      iret = iret + nf90_put_var(ncid, var_t_ID       , values=mons)
      iret = iret + nf90_put_var(ncid, var_lat_ID     , values=lat)
      iret = iret + nf90_put_var(ncid, var_lon_ID     , values=lon)
      iret = iret + nf90_put_var(ncid, var_z_ID       , values=z(1,:))
      iret = iret + nf90_put_var(ncid, var_aod_ID     , values=aod)
      iret = iret + nf90_put_var(ncid, var_ssa_ID     , values=ssa)
      iret = iret + nf90_put_var(ncid, var_asy_ID     , values=asy)
      iret = iret + nf90_put_var(ncid, var_dNovrN_ID  , values=dNovrN)
      iret = iret + nf90_put_var(ncid, var_aod_ID_2d  , values=aod_int)
      iret = iret + nf90_put_var(ncid, var_ssa_ID_2d  , values=ssa_int)
      iret = iret + nf90_put_var(ncid, var_asy_ID_2d  , values=asy_int)
      iret = iret + nf90_close(ncid)
      IF (iret /= 10*NF90_NOERR) STOP 'error writing data or in closing file'
    END IF

  END DO

END PROGRAM MACv2SP


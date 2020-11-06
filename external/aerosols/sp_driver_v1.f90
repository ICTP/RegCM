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

  INTEGER, PARAMETER :: nlev = 80, nmon = 12

  REAL, PARAMETER    :: mons(nmon) = (/1., 2., 3., 4., 5., 6., 7., 8., 9.,10.,11.,12./)
  REAL, PARAMETER    :: lambda     = 550.0

  INTEGER ::           &
       nlat          , & !< number of latitudes, defined by orography file
       nlon          , & !< number of longitudes, defined by orography file
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
       m             , & !< index for month loop
       ifile             !< output file index

  REAL ::                    &
       lat_wght            , & !< latitudinal weights for global average
       global_aod          , & !< globally averaged aerosol optical depth (diagnosis)
       global_dNovrN       , & !< global averaged change in droplet number (diagnosis)
       year_fr

  REAL, ALLOCATABLE ::       &
       lat(:)              , & !< array of latitudes allocated based on nlat
       lon(:)              , & !< array of longitudes allocated based on nlon
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
       z(:,:)              , & !< heights by colums (ncolumns,level)
       dz(:,:)             , & !< layer thicknesses (ncolumns,level)
       col_lat(:)          , & !< latitudes of columns (ncolumns)
       col_lon(:)              !< longitude of columns (ncolumns)
  !
  ! read in the orography data and use this to define the input (lat/lon) grid
  !
  iret = nf90_open("./orography_T63.nc", NF90_NOWRITE, ncid)
  IF (iret /= NF90_NOERR) STOP 'Error in opening orography file'
  iret = nf90_inq_dimid(ncid, "lat", VarID)
  iret = nf90_inquire_dimension(ncid, VarID, len = nlat)
  iret = nf90_inq_dimid(ncid, "lon", VarID)
  iret = nf90_inquire_dimension(ncid, VarID, len = nlon)

  ALLOCATE (lat(nlat))
  ALLOCATE (lon(nlon))
  ALLOCATE (oro(nlon,nlat))
  ALLOCATE (aod(nlon,nlat,nlev,nmon))
  ALLOCATE (ssa(nlon,nlat,nlev,nmon))
  ALLOCATE (asy(nlon,nlat,nlev,nmon))
  ALLOCATE (dNovrN(nlon,nlat,nmon))
  ALLOCATE (aod_int(nlon,nlat,nmon))
  ALLOCATE (ssa_int(nlon,nlat,nmon))
  ALLOCATE (asy_int(nlon,nlat,nmon))
  ALLOCATE (z(nlon,nlev))
  ALLOCATE (dz(nlon,nlev))
  ALLOCATE (col_lat(nlon))
  ALLOCATE (col_lon(nlon))
  ALLOCATE (aod_prof(nlon,nlev))
  ALLOCATE (ssa_prof(nlon,nlev))
  ALLOCATE (asy_prof(nlon,nlev))

  iret = nf90_inq_varid(ncid, "lat", VarID)
  iret = nf90_get_var(ncid, VarID, lat(:)  , start=(/1/)  ,count=(/nlat/))
  IF (iret /= NF90_NOERR) STOP 'Error in reading latitudes'
  iret = nf90_inq_varid(ncid, "lon", VarID)
  iret = nf90_get_var(ncid, VarID, lon(:)  , start=(/1/)  ,count=(/nlon/))
  IF (iret /= NF90_NOERR) STOP 'Error in reading longitudes'
  iret = nf90_inq_varid(ncid, "asl", VarID)
  iret = nf90_get_var(ncid, VarID, oro(:,:), start=(/1,1/),count=(/nlon,nlat/))
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
  !
  ! loop over annual cycle for given year, or year rangeand calculate aerosol optical properties at each lat and
  ! lon point
  ! 
  DO iyear = 2005,2005
    DO m = 1,nmon
      year_fr = iyear + (m-0.5) / nmon

      DO i = 1,nlat
        !
        ! here we process ncolumns = nlon in one call to sp_aop_profile, output is written to 2D arrays and then
        ! transfered to 4D output arrays to avoid striding through memory in call to sp_aop_profile
        !
        col_lat(:) = lat(i)
        col_lon(:) = lon(:)
        CALL sp_aop_profile                                                               ( &
             nlev         ,nlon       ,lambda     ,oro(:,i)     ,col_lon      ,col_lat    , &
             year_fr      ,z(:,:)     ,dz(:,:)    ,dNovrN(:,i,m),aod_prof     ,ssa_prof   , &
             asy_prof     )
        DO j = 1,nlon
          DO k=1,nlev
            aod(j,i,k,m) = aod_prof(j,k)
            ssa(j,i,k,m) = ssa_prof(j,k)
            asy(j,i,k,m) = asy_prof(j,k)
          END DO
        END DO
        !
        ! 2D vertically integrated values
        !     
        DO j = 1,nlon
          aod_int(j,i,m) = SUM(aod(j,i,:,m))
          ssa_int(j,i,m) = SUM(aod(j,i,:,m)*ssa(j,i,:,m))/aod_int(j,i,m)
          asy_int(j,i,m) = SUM(aod(j,i,:,m)*ssa(j,i,:,m)*asy(j,i,:,m))/(ssa_int(j,i,m)*aod_int(j,i,m))
        END DO
      END DO
    END DO
    !
    ! compute global and monthly average of selected quantities for diagnostic output
    !
    global_aod    = 0.
    global_dNovrN = 0.
    lat_wght      = 0.
    DO i=1,nlat
      lat_wght      = lat_wght      + COS(lat(i)*ASIN(1.)/90.)
      global_aod    = global_aod    + SUM(aod(:,i,:,:))*COS(lat(i)*ASIN(1.)/90.)
      global_dNovrN = global_dNovrN + SUM(dNovrN(:,i,:))*COS(lat(i)*ASIN(1.)/90.)
    END DO

    global_aod    = 0.
    global_dNovrN = 0.
    lat_wght      = 0.
    DO i=1,nlat/2
      lat_wght      = lat_wght      + COS(lat(i)*ASIN(1.)/90.)
      global_aod    = global_aod    + SUM(aod(:,i,:,:))*COS(lat(i)*ASIN(1.)/90.)
      global_dNovrN = global_dNovrN + SUM(dNovrN(:,i,:))*COS(lat(i)*ASIN(1.)/90.)
    END DO
    !
    ! write output as a netcdf file
    !
    IF (iyear == 2005) THEN
      iret = NF90_NOERR
      iret = iret + nf90_create("./MACv2-SP_3Dfields_T63_1km_2005.nc", NF90_CLOBBER, ncid)
      iret = iret + nf90_def_dim(ncid, 'lat'   ,nlat , latID)
      iret = iret + nf90_def_dim(ncid, 'lon'   ,nlon , lonID)
      iret = iret + nf90_def_dim(ncid, 'time'  ,nmon , monID)
      iret = iret + nf90_def_dim(ncid, 'z'     ,nlev , levID)
      IF (iret /= 6*NF90_NOERR) STOP 'Error in Creating File Dimensions'
      !
      iret = NF90_NOERR
      iret = iret + nf90_def_var(ncid, 'time'         , NF90_FLOAT, monID, var_t_ID)
      iret = iret + nf90_def_var(ncid, 'lat'          , NF90_FLOAT, latID, var_lat_ID)
      iret = iret + nf90_def_var(ncid, 'lon'          , NF90_FLOAT, lonID, var_lon_ID)
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
      iret = iret + nf90_put_att(ncid, var_aod_ID     , "long_name", "aerosol optical depth")
      iret = iret + nf90_put_att(ncid, var_aod_ID     , "units"    ," ")
      iret = iret + nf90_put_att(ncid, var_ssa_ID     , "long_name", "aerosol single scattering albedo")
      iret = iret + nf90_put_att(ncid, var_ssa_ID     , "units"    ," ")
      iret = iret + nf90_put_att(ncid, var_asy_ID     , "long_name", "aerosol asymmetry parameter")
      iret = iret + nf90_put_att(ncid, var_asy_ID     , "units"    ," ")
      iret = iret + nf90_put_att(ncid, var_dNovrN_ID  , "long_name", "normalized change in drop number")
      iret = iret + nf90_put_att(ncid, var_dNovrN_ID  , "units"    ," ")
      iret = iret + nf90_put_att(ncid, var_aod_ID_2d  , "long_name", "2D aerosol optical depth")
      iret = iret + nf90_put_att(ncid, var_aod_ID_2d  , "units"    ," ")
      iret = iret + nf90_put_att(ncid, var_ssa_ID_2d  , "long_name", "2D aerosol single scattering albedo")
      iret = iret + nf90_put_att(ncid, var_ssa_ID_2d  , "units"    ," ")
      iret = iret + nf90_put_att(ncid, var_asy_ID_2d  , "long_name", "2D aerosol asymmetry parameter")
      iret = iret + nf90_put_att(ncid, var_asy_ID_2d  , "units"    ," ")
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


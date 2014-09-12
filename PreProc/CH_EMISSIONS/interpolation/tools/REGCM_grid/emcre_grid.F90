! ******************************************************************
! ------------------------------------------------------------------
PROGRAM EMCRE
!EMISSIONS CREATION
! ------------------------------------------------------------------
! Author: Andrea Pozzer, ICTP, Trieste, October 2011
! BASED ON work of
! Patrick Joeckel, MPICH, Mainz, October 2004
! ******************************************************************
!
! CONVERT IPCC EMISSION FILES 
!
! -----------------------------------------------------------------

  USE mo_f2kcli                    ! command line interface
  USE emcre_tools
  USE emcre_netcdf
  USE netcdf

  IMPLICIT NONE

  ! VERSION
  CHARACTER(LEN=*), PARAMETER :: VERSION = '0.4'

  ! FOR COMMAND LINE
  CHARACTER(LEN=256) :: EXE          ! program name
  CHARACTER(LEN=256) :: CMD          ! argument
  INTEGER            :: NARG         ! number of arguments


  CHARACTER(LEN=str_vlong)   :: fname   ! filename


  ! ... DATA AND GRID
  INTEGER                                   :: ntime = 1 ! time steps

  !OUTPUT FILE SPECIFICATIONS
  INTEGER, PARAMETER :: bounds = 4 
  REAL(DP), DIMENSION(:,:),     ALLOCATABLE :: xlon    ! longitudes [deg]
  REAL(DP), DIMENSION(:,:,:),   ALLOCATABLE :: xlonb    ! longitudes [deg]
  REAL(DP), DIMENSION(:,:),     ALLOCATABLE :: xlat    ! latitudes  [deg]
  REAL(DP), DIMENSION(:,:,:),   ALLOCATABLE :: xlatb    ! latitudes  [deg]
  REAL(DP), DIMENSION(:,:),     ALLOCATABLE :: xlonc    ! longitudes [deg]
  REAL(DP), DIMENSION(:,:),     ALLOCATABLE :: xlatc    ! longitudes [deg]
  REAL(DP), DIMENSION(:,:),     ALLOCATABLE :: xlond    ! longitudes [deg]
  REAL(DP), DIMENSION(:,:),     ALLOCATABLE :: xlatd    ! longitudes [deg]
  REAL(DP), DIMENSION(:),       ALLOCATABLE :: xtime   ! time [1]
  REAL(DP), DIMENSION(:,:),     ALLOCATABLE :: var

  INTEGER  :: nlon   ! number of longitude intervals 
  INTEGER  :: nlat   ! number of latitude  intervals 
  REAL(DP),DIMENSION(:,:), ALLOCATABLE :: dhlonx  ! half longitude width  [deg]
  REAL(DP),DIMENSION(:,:), ALLOCATABLE :: dhlatx  ! half longitude width  [deg]
  REAL(DP),DIMENSION(:,:), ALLOCATABLE :: dhlony  ! half longitude width  [deg]
  REAL(DP),DIMENSION(:,:), ALLOCATABLE :: dhlaty  ! half longitude width  [deg]
  REAL(DP),DIMENSION(:,:), ALLOCATABLE :: dhlont  ! half longitude width  [deg]
  REAL(DP),DIMENSION(:,:), ALLOCATABLE :: dhlatt  ! half longitude width  [deg]

  ! VARIABLES
  INTEGER                         :: ismthlev
  INTEGER                         :: status       ! status flag
  INTEGER                         :: jc,jf        ! class/file counter
  INTEGER                         :: ji, jj, jk, jl, jt, jz, ii
  INTEGER                         :: nc           ! act. number of classes
  INTEGER                         :: nf           ! act. number of files
  LOGICAL                         :: file_exists  ! checking existence of file
  ! ... CONVERSION
  REAL(DP) :: conv                      ! conversion

  LOGICAL :: smthbdy , lakedpth, fudge_lnd , fudge_lnd_s , fudge_tex , &
             ltexture , fudge_lak_s , fudge_tex_s , fudge_lak , h2ohgt
  REAL(DP) :: h2opct
  CHARACTER(len=64) :: domname
  CHARACTER(len=256) :: dirter , inpter
  CHARACTER(len=1) :: pthsep='/'
  INTEGER :: ipunit = 101
  namelist /terrainparam/ domname , smthbdy , lakedpth , ltexture ,   &
                fudge_lnd , fudge_lnd_s , fudge_tex , fudge_tex_s ,   &
                fudge_lak, fudge_lak_s , h2opct , h2ohgt , ismthlev , &
                dirter , inpter

  ! (1) READ COMMAND LINE
  NARG = COMMAND_ARGUMENT_COUNT()    ! number of arguments
  CALL GET_COMMAND_ARGUMENT(0,EXE)   ! program name
  !
  IF (NARG > 1) THEN
     WRITE(*,*) 'COMMAND-LINE ERROR: TOO MANY ARGUMENTS !'
     STOP
  END IF
  !
  IF (NARG == 0) THEN
     WRITE(*,*) 'COMMAND-LINE ERROR: MISSING REGCM FILE INPUT'
     STOP
  END IF
  !
  CALL GET_COMMAND_ARGUMENT(1,CMD)
  open(ipunit, file=CMD, status='old', action='read', err=100)

  read(ipunit, terrainparam, err=101)

  ! (2) INIT
  fname=TRIM(DIRTER)//pthsep//TRIM(domname)//'_DOMAIN000.nc'
  write(*,*) "INPUT FROM :", TRIM(fname)

  ! (3) OPEN NETCDF-FILE DIMENSION
  CALL nc_read_dim

  ! (4) ALLOCATE FIELDS
  ALLOCATE(var(nlon,nlat))
  var(:,:)=1.0_dp
  ALLOCATE(xlon(nlon,nlat))
  ALLOCATE(xlonb(4,nlon,nlat))
  ALLOCATE(xlat(nlon,nlat))
  ALLOCATE(xlatb(4,nlon,nlat))
  xlon(:,:)=0.0_dp
  xlat(:,:)=0.0_dp
  xlonb(:,:,:)=0.0_dp
  xlatb(:,:,:)=0.0_dp
 
  ! REGCM specific
  ALLOCATE(xlonc(nlon,nlat)) ! cross
  ALLOCATE(xlond(nlon,nlat)) ! dot
  ALLOCATE(xlatc(nlon,nlat))
  ALLOCATE(xlatd(nlon,nlat))
 
  ! (6) OPEN NETCDF-FILE
  CALL nc_read

  ! (7) prepare output fields
  ! FROM 0 TO 360
  DO ji=1, nlon
    DO jj=1, nlat
     xlon(ji,jj) = xlonc(ji,jj) 
     xlat(ji,jj) = xlatc(ji,jj)
   END DO
  END DO

!  DO ji=1, nlon-1
!    DO jj=1, nlat-1
!     xlonb(1,ji,jj)= xlond(ji,jj) 
!     xlonb(2,ji,jj)= xlond(ji+1,jj) 
!     xlonb(3,ji,jj)= xlond(ji+1,jj+1) 
!     xlonb(4,ji,jj)= xlond(ji,jj+1) 
!     xlatb(1,ji,jj)= xlatd(ji,jj) 
!     xlatb(2,ji,jj)= xlatd(ji+1,jj) 
!     xlatb(3,ji,jj)= xlatd(ji+1,jj+1) 
!     xlatb(4,ji,jj)= xlatd(ji,jj+1) 
!   END DO
!  END DO

  ALLOCATE(dhlonx(nlon,nlat))
  ALLOCATE(dhlatx(nlon,nlat))
  ALLOCATE(dhlony(nlon,nlat))
  ALLOCATE(dhlaty(nlon,nlat))
  ALLOCATE(dhlont(nlon,nlat))
  ALLOCATE(dhlatt(nlon,nlat))
  DO ji=1, nlon-1
  DO jj=1, nlat-1
     dhlonx(ji,jj)=xlond(ji+1,jj)-xlond(ji,jj)
     dhlatx(ji,jj)=xlatd(ji+1,jj)-xlatd(ji,jj)
     dhlony(ji,jj)=xlond(ji,jj+1)-xlond(ji,jj)
     dhlaty(ji,jj)=xlatd(ji,jj+1)-xlatd(ji,jj)
     dhlont(ji,jj)=xlond(ji+1,jj+1)-xlond(ji,jj)
     dhlatt(ji,jj)=xlatd(ji+1,jj+1)-xlatd(ji,jj)
  ENDDO
  ENDDO
  dhlonx(:,nlat)=dhlonx(:,nlat-1)
  dhlatx(:,nlat)=dhlatx(:,nlat-1)
  dhlony(:,nlat)=dhlony(:,nlat-1)
  dhlaty(:,nlat)=dhlaty(:,nlat-1)
  dhlont(:,nlat)=dhlont(:,nlat-1)
  dhlatt(:,nlat)=dhlatt(:,nlat-1)
  dhlonx(nlon,:)=dhlonx(nlon-1,:)
  dhlatx(nlon,:)=dhlatx(nlon-1,:)
  dhlony(nlon,:)=dhlony(nlon-1,:)
  dhlaty(nlon,:)=dhlaty(nlon-1,:)
  dhlont(nlon,:)=dhlont(nlon-1,:)
  dhlatt(nlon,:)=dhlatt(nlon-1,:)

  dhlonx(nlon,nlat)=dhlonx(nlon-1,nlat)
  dhlatx(nlon,nlat)=dhlatx(nlon-1,nlat)
  dhlony(nlon,nlat)=dhlony(nlon,nlat-1)
  dhlaty(nlon,nlat)=dhlaty(nlon,nlat-1)
  dhlont(nlon,nlat)=dhlont(nlon-1,nlat-1)
  dhlatt(nlon,nlat)=dhlatt(nlon-1,nlat-1)

  DO ji=1, nlon
  DO jj=1, nlat
     xlonb(1,ji,jj)= xlond(ji,jj) 
     xlatb(1,ji,jj)= xlatd(ji,jj)
     xlonb(2,ji,jj)= xlond(ji,jj)+dhlonx(ji,jj)
     xlatb(2,ji,jj)= xlatd(ji,jj)+dhlatx(ji,jj)
     xlonb(3,ji,jj)= xlond(ji,jj)+dhlont(ji,jj)
     xlatb(3,ji,jj)= xlatd(ji,jj)+dhlatt(ji,jj)
     xlonb(4,ji,jj)= xlond(ji,jj)+dhlony(ji,jj)
     xlatb(4,ji,jj)= xlatd(ji,jj)+dhlaty(ji,jj)
  ENDDO
  ENDDO

!  xlonb(1,1,1)= xlond(1,1) 
!  xlatb(1,1,1)= xlatd(1,1)
!  xlonb(2,1,1)= xlond(1,2)
!  xlatb(2,1,1)= xlatd(1,2)
!  xlonb(3,1,1)= xlond(2,2)
!  xlatb(3,1,1)= xlatd(2,2)
!  xlonb(4,1,1)= xlond(1,1)
!  xlatb(4,1,1)= xlatd(1,1)
!  xlatb(2,1,1)= xlatd(1,2)




  ! (8) OUTPUT TO NETCDF-FILE
  CALL nc_dump

  ! (9) CLEAN
  IF (ALLOCATED(xlond)) DEALLOCATE(xlond)
  IF (ALLOCATED(xlatd)) DEALLOCATE(xlatd)
  IF (ALLOCATED(xlonc)) DEALLOCATE(xlonc)
  IF (ALLOCATED(xlatc)) DEALLOCATE(xlatc)

  IF (ALLOCATED(xlon)) DEALLOCATE(xlon)
  IF (ALLOCATED(xlat)) DEALLOCATE(xlat)
  IF (ALLOCATED(xlonb)) DEALLOCATE(xlonb)
  IF (ALLOCATED(xlatb)) DEALLOCATE(xlatb)
  IF (ALLOCATED(xtime)) DEALLOCATE(xtime)
  IF (ALLOCATED(var)) DEALLOCATE(var)

  stop

100 write ( 6, * ) 'Cannot read namelist file ', trim(CMD)
    stop
101 write ( 6, * ) 'Cannot read namelist stanza: terrainparam ', trim(CMD)
    stop

! --------------------------------------------------------------------------
! ##########################################################################
! --------------------------------------------------------------------------
CONTAINS

  ! ------------------------------------------------------------------------
  SUBROUTINE nc_read_dim
    INTRINSIC :: TRIM, ADJUSTL, NINT

    
    ! LOCAL
    INTEGER,SAVE                :: ncid   ! netCDF-ID
    INTEGER                     :: dimid_jx, dimid_iy

    CHARACTER(LEN=str_vlong)    :: name_dim   ! line


    INQUIRE(FILE=TRIM(fname), EXIST=file_exists)
    IF (.not.file_exists) THEN
      WRITE(*,*) 'File',TRIM(fname),' does NOT exist! STOP....' 
      STOP
    ENDIF

    ! OPEN FILE
    CALL NFERR( status, &
         nf90_open(TRIM(fname), NF90_NOWRITE, ncid) &
         ,1)
    ! latitude dimension check
    CALL  NFERR( status, &
         nf90_inq_dimid(ncid, 'jx', dimid_jx ) &
         ,2)
    CALL  NFERR( status, &
         nf90_Inquire_Dimension(ncid, dimid_jx, name_dim, nlon ) &
         ,3)
    ! longitude dimension check
    CALL  NFERR( status, &
         nf90_inq_dimid(ncid, 'iy', dimid_iy ) &
         ,4)
    CALL  NFERR( status, &
         nf90_Inquire_Dimension(ncid, dimid_iy, name_dim, nlat ) &
         ,5)
    WRITE (*,*) "REGCM grid dimension:"
    write(*,*) nlon,'x',nlat

    !CLOSE FILE
    CALL NFERR( status, &
         nf90_close(ncid) &
         ,24)

    ! RETURN
    status = 0

  END SUBROUTINE nc_read_dim
  ! ------------------------------------------------------------------------
  SUBROUTINE nc_read
    INTRINSIC :: TRIM, ADJUSTL, NINT

    
    ! LOCAL
    INTEGER,SAVE                :: ncid   ! netCDF-ID
    INTEGER                     :: dimid_jx, dimid_iy
    INTEGER                     :: varid_xlonc, varid_xlatc
    INTEGER                     :: varid_xlond, varid_xlatd

    CHARACTER(LEN=str_vlong)    :: name_dim   ! line


    INQUIRE(FILE=TRIM(fname), EXIST=file_exists)
    IF (.not.file_exists) THEN
      WRITE(*,*) 'File',TRIM(fname),' does NOT exist! STOP....' 
      STOP
    ENDIF

    ! OPEN FILE
    CALL NFERR( status, &
         nf90_open(TRIM(fname), NF90_NOWRITE, ncid) &
         ,21)

    CALL  NFERR( status, &
         nf90_inq_varid(ncid, 'xlat', varid_xlatc ) &
         ,22)

    CALL  NFERR( status, &
         nf90_get_var(ncid, varid_xlatc, xlatc ) &
         ,23)

    CALL  NFERR( status, &
         nf90_inq_varid(ncid, 'xlon', varid_xlonc ) &
         ,24)

    CALL  NFERR( status, &
         nf90_get_var(ncid, varid_xlonc, xlonc ) &
         ,25)

    CALL  NFERR( status, &
         nf90_inq_varid(ncid, 'dlat', varid_xlatd ) &
         ,26)

    CALL  NFERR( status, &
         nf90_get_var(ncid, varid_xlatd, xlatd ) &
         ,27)

    CALL  NFERR( status, &
         nf90_inq_varid(ncid, 'dlon', varid_xlond ) &
         ,28)

    CALL  NFERR( status, &
         nf90_get_var(ncid, varid_xlond, xlond ) &
         ,29)

    !CLOSE FILE
    CALL NFERR( status, &
         nf90_close(ncid) &
         ,30)

    ! RETURN
    status = 0

  END SUBROUTINE nc_read
  ! ------------------------------------------------------------------------
  SUBROUTINE nc_dump

!    USE netcdf

    IMPLICIT NONE

    INTRINSIC :: DATE_AND_TIME, CHAR

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'nc_dump'
    INTEGER :: ncid      ! netCDF-ID
    INTEGER :: dimid_lat, dimid_lon, dimid_bounds,dimid_latb, dimid_lonb, dimid_time
    INTEGER :: varid_lat, varid_lon, varid_bounds,varid_latb, varid_lonb, varid_time
    INTEGER :: varid_var, dimid_rank
    INTEGER :: dimid_grid, varid_rank
    INTEGER, DIMENSION (2) :: rank_var
    INTEGER                  :: jc
    CHARACTER(LEN=str_long)  :: timestr = ''
    CHARACTER(LEN=4)         :: yrstr = '    '
    CHARACTER(LEN=4)         :: yrstr_end = '    '
    CHARACTER(LEN=4)         :: yrstr_start = '    '
    CHARACTER(LEN=2000 + MAX_NCLASS*50) :: nmlstr = ''
    CHARACTER(LEN=4)         :: jcstr = '    '
    CHARACTER(LEN=4)         :: levstr = '    '
    CHARACTER(LEN=str_long)  :: scalestr = ''
    CHARACTER(LEN=str_long)  :: mmassstr = ''
    CHARACTER(LEN=str_long)  :: globssstr = ''
    CHARACTER(LEN=str_vlong) :: heightstr = ''
    CHARACTER(LEN=str_vlong) :: pressstr = ''
    CHARACTER(LEN=str_vlong) :: ipressstr = ''
    !
    CHARACTER(LEN=8)         :: date
    CHARACTER(LEN=10)        :: time
    CHARACTER(LEN=5)         :: zone

    ! INIT
    ! CONVERT TO STRING
    CALL DATE_AND_TIME(date, time, zone)

    ! CREATE NEW FILE
    CALL NFERR(status, &
         nf90_create(trim(dirter)//pthsep//'REGCM_grid.nc', &
                     NF90_CLOBBER, ncid) ,51)

    ! ADD GLOBALE ATTRIBUTES
    ! - VERSION
    CALL NFERR(status, &
         nf90_put_att(ncid, NF90_GLOBAL, 'title',     &
         'CMIP5')  &
         ,52)
    CALL NFERR(status, &
         nf90_put_att(ncid, NF90_GLOBAL, 'gridtype',     &
         'cell')  &
         ,52)
    ! - DATE AND TIME
    CALL NFERR(status, &
         nf90_put_att(ncid, NF90_GLOBAL, 'date', date) &
         ,53)
    CALL NFERR(status, &
         nf90_put_att(ncid, NF90_GLOBAL, 'time', TRIM(time)//TRIM(zone)) &
         ,54)

    ! - NAMELIST

    ! DEFINE DIMENSIONS
    CALL NFERR(status, &
         nf90_def_dim(ncid, 'grid_size', nlon*nlat, dimid_grid) &
         ,57)
    CALL NFERR(status, &
         nf90_def_dim(ncid, 'grid_xsize', nlon, dimid_lon) &
         ,56)
    CALL NFERR(status, &
         nf90_def_dim(ncid, 'grid_ysize', nlat, dimid_lat) &
         ,57)
    CALL NFERR(status, &
         nf90_def_dim(ncid, 'grid_corners', bounds, dimid_bounds) &
         ,58)
    CALL NFERR(status, &
         nf90_def_dim(ncid, 'grid_rank', 2, dimid_rank) &
         ,59)
!    CALL NFERR(status, &
!         nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid_time) &
!         ,59)

    ! DEFINE COORDINATE VARIABLES WITH ATTRIBUTES
    CALL NFERR(status, &
         nf90_def_var(ncid, 'grid_center_lon', NF90_FLOAT, (/ dimid_lon, dimid_lat /), varid_lon) &
         ,60)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_lon, 'long_name', 'longitude') &
         ,61)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_lon, 'units', 'degrees_east') &
         ,62)
!    CALL NFERR(status, &
!         nf90_put_att(ncid, varid_lon, 'bounds', 'xbounds') &
!         ,63)

    CALL NFERR(status, &
         nf90_def_var(ncid, 'grid_center_lat', NF90_FLOAT, (/ dimid_lon, dimid_lat /), varid_lat) &
         ,63)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_lat, 'long_name', 'latitude') &
         ,63)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_lat, 'units', 'degrees_north') &
         ,64)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_lat, 'bounds', 'ybounds') &
         ,65)

    CALL NFERR(status, &
         nf90_def_var(ncid, 'grid_corner_lon', NF90_FLOAT, (/ dimid_bounds, dimid_lon, dimid_lat /), varid_lonb) &
         ,63)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_lonb, 'long_name', 'longitude_bounds') &
         ,63)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_lonb, 'units', 'degrees_east') &
         ,64)

    CALL NFERR(status, &
         nf90_def_var(ncid, 'grid_corner_lat', NF90_FLOAT, (/ dimid_bounds, dimid_lon, dimid_lat /), varid_latb) &
         ,63)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_latb, 'long_name', 'latitude_bounds') &
         ,63)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_latb, 'units', 'degrees_north') &
         ,64)
    
   CALL NFERR(status, &
         nf90_def_var(ncid, 'grid_dims', NF90_FLOAT, (/ dimid_rank /), varid_rank) &
         ,65)

!    CALL NFERR(status, &
!         nf90_def_var(ncid, 'time', NF90_FLOAT, (/ dimid_time /), varid_time) &
!         ,68)
!    CALL NFERR(status, &
!         nf90_put_att(ncid, varid_time, 'long_name', 'time') &
!         ,69)
!    !
!    timestr = 'month since '//yrstr_start//'-01-15 00:00:00'
!    CALL NFERR(status, &
!         nf90_put_att(ncid, varid_time, 'units', timestr) &
!         ,70)
!

! VARIABLE FAKE =1
!    ! - flux
    CALL NFERR(status, &
         nf90_def_var(ncid, 'fixed_value', NF90_FLOAT  &
         , (/ dimid_lon, dimid_lat /), varid_var) &
         ,78)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_var, 'units', ' ') &
         ,80)


    ! SWITCH MODUS
    CALL NFERR(status, &
         nf90_enddef(ncid) &
         ,82)

    ! SAVE COORDINATE VARIBLES
    CALL NFERR(status, &
         nf90_put_var(ncid, varid_lon, xlon) &
         ,83)
    CALL NFERR(status, &
         nf90_put_var(ncid, varid_lat, xlat) &
         ,84)
!    CALL NFERR(status, &
!         nf90_put_var(ncid, varid_time, xtime) &
!         ,86)

    CALL NFERR(status, &
         nf90_put_var(ncid, varid_lonb, xlonb) &
         ,89)
    CALL NFERR(status, &
         nf90_put_var(ncid, varid_latb, xlatb) &
         ,89)

    rank_var(1)=nlon
    rank_var(2)=nlat

    CALL NFERR(status, &
         nf90_put_var(ncid, varid_rank, rank_var) &
         ,89)

    CALL NFERR(status, &
         nf90_put_var(ncid, varid_var, var) &
         ,90)

    ! CLOSE FILE
    CALL NFERR(status, &
         nf90_close(ncid) &
         ,99)

  END SUBROUTINE nc_dump
  ! ------------------------------------------------------------------------

  ! ------------------------------------------------------------------
  SUBROUTINE NFERR(status,command,pos)

    USE netcdf, ONLY: NF90_NOERR, nf90_strerror
    
    IMPLICIT NONE
    
    ! I/O
    INTEGER,          INTENT(OUT) :: status
    INTEGER,          INTENT(IN) :: command
    INTEGER,          INTENT(IN) :: pos
    
    status=command
    IF (status /= NF90_NOERR) THEN
       WRITE(*,*) 'netCDF ERROR at position: ', pos
       WRITE(*,*) 'netCDF ERROR status     : ',status
       WRITE(*,*) 'netCDF ERROR            : ',nf90_strerror(status)
    END IF
  
  END SUBROUTINE NFERR
  ! ------------------------------------------------------------------



! --------------------------------------------------------------------------
! ##########################################################################
END PROGRAM EMCRE
! ##########################################################################
! --------------------------------------------------------------------------

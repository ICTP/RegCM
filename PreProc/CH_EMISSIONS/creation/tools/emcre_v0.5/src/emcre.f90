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

  IMPLICIT NONE

  ! VERSION
  CHARACTER(LEN=*), PARAMETER :: VERSION = '0.4'

  ! FOR COMMAND LINE
  CHARACTER(LEN=256) :: EXE          ! program name
  CHARACTER(LEN=80)  :: CMD          ! argument
  INTEGER            :: NARG         ! number of arguments


  ! NAMELIST
  TYPE SOURCE_IO
     CHARACTER(LEN=str_long) :: name = ''
  END TYPE SOURCE_IO
  TYPE YEAR_IO
     INTEGER, DIMENSION(MAX_SPEC):: year = 0
  END TYPE YEAR_IO
  TYPE FNAME_IO
     CHARACTER(LEN=str_long), DIMENSION(MAX_SPEC) :: fname = ''
  END TYPE FNAME_IO
  TYPE FRAC_IO
     REAL(DP), DIMENSION(MAX_HEIGHTS) :: scale = -999.0_DP
  END TYPE FRAC_IO

  !
  ! NAMELIST VARIABLES
  LOGICAL                               :: VERBOSE   = .TRUE.
  CHARACTER(LEN=str_long)               :: OUTPUT    = ''     ! netCDF-file
  CHARACTER(LEN=str_long)               :: SPECIES   = ''
  CHARACTER(LEN=str_long)               :: TIME_FREQUENCY   = ''
  CHARACTER(LEN=str_long)               :: OUT_UNIT   = ''
  REAL(DP)                              :: MOLARMASS = 0.0_DP ! [g/mol]
  ! global scaling factor
  REAL(DP)                              :: GLOBALSCALE = 0.0_DP
  ! Interpolation options
  CHARACTER(LEN=str_long)               :: INTERPOLATION='linear'
  ! special case for aircraft emissions (3D)
  LOGICAL                               :: L_AIRCRAFT = .FALSE.
  INTEGER                               :: YEAR_START    = 0        
  INTEGER                               :: YEAR_END      = 0        
  REAL(DP), DIMENSION(MAX_HEIGHTS)      :: HEIGHT              ! [m]
  REAL(DP), DIMENSION(MAX_HEIGHTS)      :: PRESS               ! [Pa]
  REAL(DP), DIMENSION(MAX_HEIGHTS)      :: IPRESS              ! [Pa]
  REAL(DP)                              :: delta               ! [m] -> thickness of layers
  CHARACTER(LEN=str_vlong)              :: INPUTPATH = ''
  TYPE(SOURCE_IO), DIMENSION(MAX_NCLASS):: SOURCE 
  TYPE(YEAR_IO),   DIMENSION(MAX_NCLASS):: YEAR 
  TYPE(FNAME_IO),  DIMENSION(MAX_NCLASS):: FILE_NAME
  TYPE(FRAC_IO),  DIMENSION(MAX_NCLASS) :: FRAC              ! emission class

  ! INPUT FROM FILE:
  CHARACTER(LEN=str_long)             :: freq_file = ''
  CHARACTER(LEN=str_long)             :: unit_file = ''
  INTEGER                             :: nlon_file    ! longitudes [deg]
  INTEGER                             :: nlat_file    ! latitudes  [deg]
  INTEGER                             :: ntime_file = 1 ! time steps
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: data_file
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: data_file_air

  ! ... DATA AND GRID
  INTEGER                                   :: nlev  = 0 ! levels
  INTEGER                                   :: nlev_tmp = 0 ! levels test
  INTEGER                                   :: ntime = 1 ! time steps
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE   :: emisflux_class
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: emisflux_air ! special case aircraft
  REAL(DP), DIMENSION(:,:,:,:), ALLOCATABLE :: emisflux, emistmp

  ! Interpolation needed
  INTEGER                             :: month_offset  
  INTEGER                             :: nyears_added  = 0        
  INTEGER, DIMENSION(MAX_SPEC)        :: years_location  = 0        
  REAL(DP), DIMENSION(:,:),   ALLOCATABLE :: slope
  REAL(DP), DIMENSION(:),     ALLOCATABLE :: x_vector
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: y_vector
  REAL(DP), DIMENSION(:,:,:), ALLOCATABLE :: spline_inout
  REAL(DP)                              :: xval
  REAL(DP), DIMENSION(:,:), ALLOCATABLE :: yval

  !OUTPUT FILE SPECIFICATIONS
  REAL(DP), DIMENSION(:),       ALLOCATABLE :: xlon    ! longitudes [deg]
  REAL(DP), DIMENSION(:),       ALLOCATABLE :: xlat    ! latitudes  [deg]
  REAL(DP), DIMENSION(:),       ALLOCATABLE :: xlev    ! levels [1]
  REAL(DP), DIMENSION(:),       ALLOCATABLE :: xilev   ! interface levels [xlev+1]
  REAL(DP), DIMENSION(:),       ALLOCATABLE :: xtime   ! time [1]

  INTEGER, PARAMETER  :: nlon = 720   ! number of longitude intervals (1 deg)
  INTEGER, PARAMETER  :: nlat = 360   ! number of latitude  intervals (1 deg)
  REAL(DP) :: dhlat  ! half latitude width  [deg]
  REAL(DP) :: dhlon  ! half longitude width  [deg]

  ! VARIABLES
  INTEGER                         :: status       ! status flag
  INTEGER                         :: jc,jf        ! class/file counter
  INTEGER                         :: ji, jj, jk, jl, jt, jz, ii
  INTEGER                         :: nc           ! act. number of classes
  INTEGER                         :: nf           ! act. number of files
  LOGICAL                         :: file_exists  ! checking existence of file
  ! ... CONVERSION
  REAL(DP) :: conv                      ! conversion

  ! (1) READ COMMAND LINE
  NARG = COMMAND_ARGUMENT_COUNT()    ! number of arguments
  CALL GET_COMMAND_ARGUMENT(0,EXE)   ! program name
  !
  IF (NARG > 1) THEN
     WRITE(*,*) 'COMMAND-LINE ERROR: TOO MANY ARGUMENTS !'
     CALL USAGE(TRIM(EXE))
     STOP
  END IF
  !
  IF (NARG == 0) THEN
     CALL USAGE(TRIM(EXE))
     STOP
  END IF
  !
  CALL GET_COMMAND_ARGUMENT(1,CMD)

  ! (2) INIT
  nc = 0                 ! act. number of classes
  HEIGHT(:) = -999.0_dp  ! emission height [m]
  PRESS(:)  = -999.0_dp  ! emission height [m]
  IPRESS(:) = -999.0_dp  ! emission height [m]
  ALLOCATE(xlon(nlon))
  dhlat=180._dp/(nlat*2)
  dhlon=360._dp/(nlon*2)
  ! FROM 0 TO 360
  DO ji=1, nlon
     xlon(ji) = 0.0_DP -dhlon + 2*dhlon*REAL(ji, DP)
  END DO
  ALLOCATE(xlat(nlat))
  DO jj=1, nlat
     xlat(jj) = -90.0_DP -dhlat + 2*dhlat*REAL(jj, DP)
  END DO

  ! (3) READ NAMELIST-FILE
  CALL read_nml(status, iou, TRIM(CMD))
  IF (status /= 0) STOP

  !(3a) INIT TIME AXIS
  if (TRIM(TIME_FREQUENCY)=='monthly') then 
     ntime= (year_end-year_start+1)*12
  elseif (TRIM(TIME_FREQUENCY)=='annual') then
     ntime= (year_end-year_start+1)
  else
     WRITE(*,*) 'UNKNOWN OUTPUT FREQUENCY ',TRIM(TIME_FREQUENCY)
     STOP
  endif

  ALLOCATE(xtime(ntime))
  DO jl=1, ntime
     xtime(jl) = REAL(jl-1, DP)
  END DO
  write (*,*) 'NUM. TIME STEPS: ', ntime
    
  ! (4) GET NUMBER OF EMISSION LEVELS
  DO jc = 1, MAX_HEIGHTS
     if (HEIGHT(jc).ne.-999.0_dp) nlev = nlev+1
  ENDDO
  IF (MINVAL(HEIGHT(1:nlev)).eq.-999.0_dp) THEN
    WRITE(*,*) 'PROBLEM IN HEIGHT DEFINITION (missing levels?)'
    STOP
  ENDIF
  IF (L_AIRCRAFT) THEN
  ! PRESSURE and INTEFACE PRESSURE (IPRESS) CALCULATION 
  ! BASED ON THE INTERNATIONAL STANDARD ATMOSPHERE: see home.anadolu.edu.tr/~mcavcar/common/ISAweb.pdf : 
  ! p=p0*(1-0.0065*(h/T0))^(5.2561)  with  p0=1013325 Pa, T0 = 288.15 K, h in meters
    DO jc = 1, nlev
      PRESS(jc)= 1013325.0*(1-0.0065*(HEIGHT(jc)/288.15))**(5.2561) 
    ENDDO
    delta=HEIGHT(1)
    DO jc = 1, nlev
      IPRESS(jc)= 1013325*(1-0.0065*((HEIGHT(jc)-delta)/288.15))**(5.2561) 
    ENDDO
    IPRESS(nlev+1)=1013325*(1-0.0065*((HEIGHT(jc)+delta)/288.15))**(5.2561)
    write (*,*) 'press',press(1:nlev)
    write (*,*) 'ipress',ipress(1:nlev+1)
    
     
  ENDIF
  write (*,*) 'NUM.LEVEL: ', nlev

  ALLOCATE(xlev(nlev))
  DO jk=1, nlev
     xlev(jk) = REAL(jk, DP)
  END DO
  IF (L_AIRCRAFT) THEN
    ALLOCATE(xilev(nlev+1))
    DO jk=1, nlev+1
       xilev(jk) = REAL(jk, DP)
    END DO
  ENDIF

  ! ALLOCATE GLOBAL EMISSIONS FLUXES!
  ALLOCATE(emisflux(nlon, nlat, nlev, ntime))
  emisflux(:,:,:,:) = 0.0

  ! (5) LOOP OVER EMISSION-CLASSES
  class_loop: DO jc = 1, MAX_NCLASS

     ! SKIP IF NAME IS EMPTY
     IF (TRIM(source(jc)%name) == '') CYCLE

     ! ALLOCATE EMISSIONS FLUXES FOR THIS CLASS
     ALLOCATE(emisflux_class(nlon, nlat, ntime))
     emisflux_class(:,:,:) = 0.0
     IF (L_AIRCRAFT) THEN
       ALLOCATE(emisflux_air(nlon, nlat, nlev, ntime))
       emisflux_air(:,:,:,:) = 0.0
     ENDIF

     
     write(*,*) 'SOURCE NAME:', source(jc)%name
 
     nlev_tmp=0
     DO jk=1, MAX_HEIGHTS
       if (FRAC(jc)%scale(jk).ne.-999.0_dp) nlev_tmp = nlev_tmp+1
     ENDDO
     IF (nlev_tmp.ne.nlev) THEN
       WRITE(*,*) 'PROBLEM IN FRACTION DEFINITION (missing levels?)'
       WRITE(*,*) ' LEVELS DIFFERES FROM DEFINED HEIGHTS !!'
       STOP
     ENDIF
     IF (SUM(FRAC(jc)%scale(1:nlev)).ne.1.0_dp) THEN
       WRITE(*,*) 'WARNING LEVEL FRACTION SUM <> 1 !!!'
     ENDIF

     file_loop: DO jf=1,MAX_SPEC
        
        IF (YEAR(jc)%year(jf)==0.or.FILE_NAME(jc)%fname(jf)=="") CYCLE 

        IF (YEAR(jc)%year(jf).lt.YEAR_START.or.YEAR(jc)%year(jf).gt.YEAR_END) THEN
           WRITE(*,*) "YEAR NOT ADDED (OUTSIDE RANGE YEAR_START:YEAR_END)"
           CYCLE
        ENDIF 
        
        WRITE(*,*) '==========================================================='
        write(*,*) 'YEAR : ', YEAR(jc)%year(jf)  
        write(*,*) 'file : ', FILE_NAME(jc)%fname(jf)  
        
        !INQUIRE DIMENSION OF INPUT FILE
        CALL  inquire_ipcc(status, file_exists, TRIM(INPUTPATH)//'/'//TRIM(ADJUSTL(FILE_NAME(jc)%fname(jf))), &
                  TRIM(source(jc)%name), freq_file, unit_file,       &
                  nlon_file, nlat_file, ntime_file)
!        IF (status > 0) STOP   ! ERROR
        
        !ALLOCATE VARIABLE FROM FILENAME FOR READING THE DATA
        IF (L_AIRCRAFT) THEN
          ALLOCATE(data_file_air(nlon_file,nlat_file,nlev,ntime_file))
          data_file_air(:,:,:,:) = 0.0_dp
        ELSE
          ALLOCATE(data_file(nlon_file,nlat_file,ntime_file))
          data_file(:,:,:) = 0.0_dp
        ENDIF

        IF (file_exists) THEN
           IF (L_AIRCRAFT) THEN
              CALL  read_ipcc_air(status, TRIM(INPUTPATH)//'/'//TRIM(ADJUSTL(FILE_NAME(jc)%fname(jf))), &
                      TRIM(source(jc)%name), data_file_air)
              IF(VERBOSE) THEN
                 WRITE(*,*) '-----------------------------------------------------------'
                 WRITE(*,*) ' FILE SUMMARY: '
                 WRITE(*,*) 'UNIT           : ',TRIM(unit_file)
                 WRITE(*,*) 'FREQUENCY      : ',TRIM(freq_file)
                 WRITE(*,*) 'MIN            : ',MINVAL(data_file_air)
                 WRITE(*,*) 'MAX            : ',MAXVAL(data_file_air)
                 WRITE(*,*) 'SUM            : ',SUM(data_file_air)
                 WRITE(*,*) '-----------------------------------------------------------'
              ENDIF
           ELSE
             !READ VARIABLE FROM FILEN 
              CALL  read_ipcc(status, TRIM(INPUTPATH)//'/'//TRIM(ADJUSTL(FILE_NAME(jc)%fname(jf))), &
                      TRIM(source(jc)%name), data_file)
              IF(VERBOSE) THEN
                 WRITE(*,*) '-----------------------------------------------------------'
                 WRITE(*,*) ' FILE SUMMARY: '
                 WRITE(*,*) 'UNIT           : ',TRIM(unit_file)
                 WRITE(*,*) 'FREQUENCY      : ',TRIM(freq_file)
                 WRITE(*,*) 'MIN            : ',MINVAL(data_file)
                 WRITE(*,*) 'MAX            : ',MAXVAL(data_file)
                 WRITE(*,*) 'SUM            : ',SUM(data_file)
                 WRITE(*,*) '-----------------------------------------------------------'
              ENDIF
           ENDIF
        ENDIF

      !ASSIGN READ DATA TO THE RELATIVE YEAR IN THE "LOCAL/CLASS" ARRAY
        
      !OUTPUT FREQUENCY ANNUAL/MONTHLY
      IF (L_AIRCRAFT) THEN
        SELECT CASE (TRIM(TIME_FREQUENCY))
        CASE ('annual')
           jt=(YEAR(jc)%year(jf)-YEAR_START)+1
             !INPUT FREQUENCY ANNUAL/MONTHLY
             SELECT CASE(TRIM(ADJUSTL(freq_file)))
             !NO VERTICAL LEVEL YET... AFTER INTERPOLATION TO ACCELERATE CODE
             CASE('annual')
              emisflux_air(:,:,:,jt) = emisflux_air(:,:,:,jt) +  &
                                          data_file_air(:,:,:,1)
             CASE('monthly')
              ! brutal monthly averages (not considering different n. days/mont
              ! the error is of ~ 0.21 %
              emisflux_air(:,:,:,jt) = emisflux_air(:,:,:,jt) +  &
                                          SUM(data_file_air(:,:,:,1:12),DIM=4)/12.0
             END SELECT
             ! ADDED ONE LOCATION TO LOCALL ARRAY
             nyears_added = nyears_added+1
             years_location(nyears_added)=jt
        CASE('monthly')
           jt=(YEAR(jc)%year(jf)-YEAR_START)*12+1
             !INPUT FREQUENCY ANNUAL/MONTHLY
             SELECT CASE(TRIM(ADJUSTL(freq_file)))
             CASE('annual')
                DO ii=1,12
                  emisflux_air(:,:,:,jt+ii-1) = emisflux_air(:,:,:,jt+ii-1) +  &
                                                data_file_air(:,:,:,1)
                ENDDO
             CASE('monthly')
              emisflux_air(:,:,:,jt:jt+12) = emisflux_air(:,:,:,jt:jt+12) +  &
                                          data_file_air(:,:,:,1:12)
             END SELECT
             ! ADDED ONE LOCATION TO LOCALL ARRAY
             nyears_added = nyears_added+1
             years_location(nyears_added)=jt
        END SELECT
      ELSE
        SELECT CASE (TRIM(TIME_FREQUENCY))
        CASE ('annual')
           jt=(YEAR(jc)%year(jf)-YEAR_START)+1
             !INPUT FREQUENCY ANNUAL/MONTHLY
             SELECT CASE(TRIM(ADJUSTL(freq_file)))
             !NO VERTICAL LEVEL YET... AFTER INTERPOLATION TO ACCELERATE CODE
             CASE('annual')
              emisflux_class(:,:,jt) = emisflux_class(:,:,jt) +  &
                                          data_file(:,:,1)
             CASE('monthly')
              ! brutal monthly averages (not considering different n. days/mont
              ! the error is of ~ 0.21 %
!FAB 
             print*,size(data_file,3)
             emisflux_class(:,:,jt) = emisflux_class(:,:,jt) +  &
                                          SUM(data_file(:,:,1:12),DIM=3)/12.0
             END SELECT
             ! ADDED ONE LOCATION TO LOCALL ARRAY
             nyears_added = nyears_added+1
             years_location(nyears_added)=jt
        CASE('monthly')
           jt=(YEAR(jc)%year(jf)-YEAR_START)*12+1
             !INPUT FREQUENCY ANNUAL/MONTHLY
             SELECT CASE(TRIM(ADJUSTL(freq_file)))
             !NO VERTICAL LEVEL YET... AFTER INTERPOLATION TO ACCELERATE CODE
             CASE('annual')
                DO ii=1,12
                  emisflux_class(:,:,jt+ii-1) = emisflux_class(:,:,jt+ii-1) +  &
                                                data_file(:,:,1)
                ENDDO
             CASE('monthly')
!FAB              emisflux_class(:,:,jt:jt+12) = emisflux_class(:,:,jt:jt+12) +  &
!                                          data_file(:,:,1:12)
             emisflux_class(:,:,jt:jt+11) = emisflux_class(:,:,jt:jt+11) +  &
                                          data_file(:,:,1:12)

             END SELECT
             ! ADDED ONE LOCATION TO LOCALL ARRAY
             nyears_added = nyears_added+1
             years_location(nyears_added)=jt
        END SELECT
      ENDIF

      !DEALLOCATE VARIABLE FOR READING FROM FILE
      IF (L_AIRCRAFT) THEN
        DEALLOCATE(data_file_air)
      ELSE
        DEALLOCATE(data_file)
      ENDIF

    END DO file_loop 
    
   IF (L_AIRCRAFT) THEN
     DO jk=1, nlev
       CALL INTERP(emisflux_air(:,:,jk,:), nyears_added, years_location(:))
     ENDDO
   ELSE
     CALL INTERP(emisflux_class(:,:,:), nyears_added, years_location(:))
   ENDIF
   nyears_added=0
   years_location(:)=0

   do jk=1,nlev 
     IF (L_AIRCRAFT) THEN
       emisflux(:,:,jk,:) = emisflux(:,:,jk,:) + emisflux_air(:,:,jk,:)*FRAC(jc)%scale(jk)
     ELSE
       emisflux(:,:,jk,:) = emisflux(:,:,jk,:) + emisflux_class(:,:,:)*FRAC(jc)%scale(jk)
     ENDIF
   enddo

   !COPY EMISSIONS CLASS TO GLOBAL ARRAY
   ! AND DISTERIBUTE ON VERTICAL LEVELS!

   ! RESET LOCAL ARRAY FOR NEXT CLASS
   DEALLOCATE(emisflux_class)

  END DO class_loop

  ! (6) CONVERT TO mcl/(m^2s) and scale
  ! Kg m^-2 s^-1  -->  mcl m^-22 s^-1
   IF ((OUT_UNIT)=='mcl m-2 s-1') then
      conv   = ( N_A * 1000.0_DP )/(MOLARMASS)
      emisflux(:,:,:,:) = emisflux(:,:,:,:) * GLOBALSCALE * conv 
   ELSEIF ((OUT_UNIT)=='Kg m-2 s-1') then
      conv   = 1.0_dp
      emisflux(:,:,:,:) = emisflux(:,:,:,:) * GLOBALSCALE * conv 
   ELSEIF ((OUT_UNIT)=='mcl m-3 s-1') then
      write (*,*) 'AIRCRAFT EMISSIONS => mlc m-3 s-1!'
      conv   = ( N_A * 1000.0_DP )/(MOLARMASS)
      emisflux(:,:,:,:) = emisflux(:,:,:,:) * GLOBALSCALE * conv
   ELSE
      WRITE(*,*) " UNKNOWN FLUX UNIT"
      WRITE(*,*) " ADD NEW CONVERSION!"
      STOP
   ENDIF
!FAB version AMMA
!grid -180 180 to 0 360 as defined for xlo,xlat 
   If (1==2) then 
   ALLOCATE(emistmp(nlon, nlat, nlev, ntime))
    emistmp(361:720,:,:,:) =  emisflux(1:360,:,:,:)
    emistmp(1:360,:,:,:)= emisflux(361:720,:,:,:)
    emisflux(:,:,:,:) = emistmp(:,:,:,:)
   end if
!
  ! (7) DIAGNOSTIC OUTPUT
  IF (VERBOSE) THEN
    WRITE(*,*) '==========================================================='
    WRITE(*,*) 'VERBOSE        : ', VERBOSE
    WRITE(*,*) 'OUTPUT         : ', TRIM(OUTPUT)
    WRITE(*,*) 'SPECIES        : ', TRIM(SPECIES)
    WRITE(*,*) 'TIME FREQUENCY : ', TRIM(TIME_FREQUENCY)
    WRITE(*,*) 'MOLAR MASS     : ', MOLARMASS
    WRITE(*,*) 'UNIT           : ', OUT_UNIT
    WRITE(*,*) 'GLOBAL SCALING : ', GLOBALSCALE
    WRITE(*,*) 'YEAR_START     : ', YEAR_START
    WRITE(*,*) 'YEAR_END       : ', YEAR_END
    WRITE(*,*) 'HEIGHT         : ', HEIGHT(1:nlev)
    IF (L_AIRCRAFT) THEN
      WRITE(*,*) 'PRESS        : ', PRESS(1:nlev)
      WRITE(*,*) 'IPRESS       : ', IPRESS(1:nlev+1)
    ENDIF
    WRITE(*,*) 'INPUTPATH      : ', TRIM(INPUTPATH)
    WRITE(*,*)
    WRITE(*,*) 'FLUX           : '
    WRITE(*,*) ' LBOUND        : ',LBOUND(emisflux)
    WRITE(*,*) ' UBOUND        : ',UBOUND(emisflux)
    WRITE(*,*) ' MIN           : ',MINVAL(emisflux)
    WRITE(*,*) ' MAX           : ',MAXVAL(emisflux)
    WRITE(*,*) '==========================================================='
  ENDIF

  ! (8) OUTPUT TO NETCDF-FILE
  CALL nc_dump

  ! (9) CLEAN
  IF (ALLOCATED(emisflux)) DEALLOCATE(emisflux)
  IF (ALLOCATED(xlon)) DEALLOCATE(xlon)
  IF (ALLOCATED(xlat)) DEALLOCATE(xlat)
  IF (ALLOCATED(xlev)) DEALLOCATE(xlev)
  IF (ALLOCATED(xtime)) DEALLOCATE(xtime)
  IF (L_AIRCRAFT) THEN
    IF (ALLOCATED(xilev)) DEALLOCATE(xilev)
  ENDIF

! --------------------------------------------------------------------------
! ##########################################################################
! --------------------------------------------------------------------------
CONTAINS

  ! ------------------------------------------------------------------------
  SUBROUTINE USAGE(EXE)
    CHARACTER (LEN=*) :: EXE
    WRITE(*,*) '--------------------------------------------'
    WRITE(*,*) 'EDGAR2NC Version ',VERSION
    WRITE(*,*) 'Author: Patrick Joeckel, MPICH, October 2004'
    WRITE(*,*) '--------------------------------------------'
    WRITE(*,*) 'Usage: '//TRIM(EXE)//' <namelist-file>'
    WRITE(*,*) '--------------------------------------------'
  END SUBROUTINE USAGE
  ! ------------------------------------------------------------------------

  ! ------------------------------------------------------------------------
  SUBROUTINE read_nml(status, iou, fname)

    IMPLICIT NONE

    ! I/O
    INTEGER,          INTENT(OUT) :: status
    INTEGER,          INTENT(IN)  :: iou
    CHARACTER(LEN=*), INTENT(IN)  :: fname

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'read_nml'
    LOGICAL :: lex   ! file exists
    INTEGER :: fstat ! file status

    NAMELIST /CTRL/ VERBOSE, OUTPUT, SPECIES, MOLARMASS, TIME_FREQUENCY, OUT_UNIT  &
                  , GLOBALSCALE, INTERPOLATION, L_AIRCRAFT   & 
                  , YEAR_START, YEAR_END, HEIGHT, INPUTPATH &
                  , SOURCE, FRAC, YEAR, FILE_NAME
  
    !DEFAULT
    TIME_FREQUENCY ='monthly'
    OUT_UNIT ='mcl m-2 s-1'
    INTERPOLATION ='linear'
    L_AIRCRAFT=.FALSE.
    !
 
    status = 1 ! ERROR

    WRITE(*,*) '==========================================================='

    ! CHECK IF FILE EXISTS
    INQUIRE(file=TRIM(fname), exist=lex)
    IF (.NOT.lex) THEN
       WRITE(*,*) substr,': FILE DOES NOT EXIST (',TRIM(fname),')'
       status = 1
       RETURN
    END IF

    ! OPEN FILE
    OPEN(iou,file=TRIM(fname))

    ! READ NEMELIST
    WRITE(*,*) 'READING NAMELIST ''CTRL'''//&
         &' FROM '''//TRIM(fname),''' (unit ',iou,') ...'
    !
    READ(iou, NML=CTRL, IOSTAT=fstat)
    !
    IF (fstat /= 0) THEN
       WRITE(*,*) substr,': READ ERROR IN NAMELIST ''CTRL'' (',TRIM(fname),')'
       status = 3  ! READ ERROR IN NAMELIST
       RETURN
    END IF

    WRITE(*,*) ' OUTPUT     : ', TRIM(OUTPUT)
    WRITE(*,*) ' SPECIES    : ', TRIM(SPECIES)
    WRITE(*,*) ' YEAR_START : ', YEAR_START
    WRITE(*,*) ' YEAR_END   : ', YEAR_END
    WRITE(*,*) ' INPUTPATH  : ', TRIM(INPUTPATH)

    ! CLOSE FILE
    CLOSE(iou)

    WRITE(*,*) '==========================================================='

    status = 0

  END SUBROUTINE read_nml
  ! ------------------------------------------------------------------------

  SUBROUTINE inquire_ipcc(status, file_exists, fname, var,  &
             freq_file, unit_file,         &
             nlon_file, nlat_file, ntime_file   )

    USE netcdf

    IMPLICIT NONE

    INTRINSIC :: TRIM, ADJUSTL, NINT

    ! I/O
    INTEGER, INTENT(OUT)                   :: status
    CHARACTER(LEN=*),          INTENT(IN)  :: fname   ! filename
    CHARACTER(LEN=*),          INTENT(IN)  :: var     ! variable name
    CHARACTER(LEN=str_short),  INTENT(OUT) :: freq_file    ! emissions frequency
    CHARACTER(LEN=str_short),  INTENT(OUT) :: unit_file
    INTEGER,                   INTENT(OUT) :: nlon_file  ! IPCC dimension lenght 
    INTEGER,                   INTENT(OUT) :: nlat_file  ! IPCC dimension length
    INTEGER,                   INTENT(OUT) :: ntime_file ! IPCC dimension length
    LOGICAL,                   INTENT(OUT) :: file_exists! existence of file 
    
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'inquire_ipcc'
    INTEGER,SAVE                :: ncid   ! netCDF-ID
    INTEGER                     :: dimid, varid
 
    INTEGER, DIMENSION(:), ALLOCATABLE  :: date_file
  
    CHARACTER(LEN=str_vlong)    :: name_dim   ! line


    INQUIRE(FILE=TRIM(fname), EXIST=file_exists)
    IF (.not.file_exists) THEN
      WRITE(*,*) 'File',TRIM(fname),' does NOT exist! skipping....' 
      RETURN
    ENDIF

    ! OPEN FILE
    CALL NFERR( status, &
         nf90_open(TRIM(fname), NF90_NOWRITE, ncid) &
         ,1)
    ! latitude dimension check
    CALL  NFERR( status, &
         nf90_inq_dimid(ncid, 'lat', dimid ) &
         ,2)
    CALL  NFERR( status, &
         nf90_Inquire_Dimension(ncid, dimid, name_dim, nlat_file ) &
         ,3)
    IF (nlat_file.ne.nlat) THEN 
      WRITE (*,*) "LAT differs between grid definition and imported grid"
      write(*,*) nlat_file, nlat
      STOP
    ENDIF
    ! longitude dimension check
    CALL  NFERR( status, &
         nf90_inq_dimid(ncid, 'lon', dimid ) &
         ,4)
    CALL  NFERR( status, &
         nf90_Inquire_Dimension(ncid, dimid, name_dim, nlon_file ) &
         ,5)
    IF (nlon_file.ne.nlon) THEN 
      WRITE (*,*) "LON differs between grid definition and imported grid"
      write(*,*) nlon_file, nlon
      STOP
    ENDIF
    ! time dimension check
    CALL  NFERR( status, &
         nf90_inq_dimid(ncid, 'time', dimid ) &
         ,6)
    CALL  NFERR( status, &
         nf90_Inquire_Dimension(ncid, dimid, name_dim, ntime_file ) &
         ,7)
    IF (ntime_file.eq.12) THEN 
      freq_file = 'monthly'
    ELSE IF (ntime_file.eq.1) THEN 
      freq_file = 'annual'
    ELSE
       WRITE(*,*) substr, &
            ': UNKNOWN DATA FREQUENCY',TRIM(ADJUSTL(freq_file))
       status = 5  ! UNKNOWN FREQUENCY
       RETURN
    ENDIF


    !CLOSE FILE
    CALL NFERR( status, &
         nf90_close(ncid) &
         ,14)

    ! RETURN
    status = 0
    
  END SUBROUTINE inquire_ipcc
!  ! ------------------------------------------------------------------------

  ! ------------------------------------------------------------------------
  SUBROUTINE read_ipcc(status, fname, var,  data_file)

    USE netcdf

    IMPLICIT NONE

    INTRINSIC :: TRIM, ADJUSTL, NINT

    ! I/O
    INTEGER, INTENT(OUT)                   :: status
    CHARACTER(LEN=*),          INTENT(IN)  :: fname   ! filename
    CHARACTER(LEN=*),          INTENT(IN)  :: var     ! variable name
    REAL(DP), DIMENSION(:,:,:), INTENT(OUT):: data_file    ! INTENT(OUT)
    
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'read_ipcc'
    INTEGER,SAVE                :: ncid   ! netCDF-ID
    INTEGER                     :: dimid, varid

    CHARACTER(LEN=str_vlong)    :: name_dim   ! line

    ! OPEN FILE
    CALL NFERR( status, &
         nf90_open(TRIM(fname), NF90_NOWRITE, ncid) &
         ,21)

    CALL  NFERR( status, &
         nf90_inq_varid(ncid, var, varid ) &
         ,22)
    IF (status.ne.0) then
        write(*,*) "variable not found in NetCDF file, skipping"
        RETURN
    ENDIF

    CALL  NFERR( status, &
         nf90_get_var(ncid, varid, data_file ) &
         ,23)

    !CLOSE FILE
    CALL NFERR( status, &
         nf90_close(ncid) &
         ,24)

    ! RETURN
    status = 0
    
  END SUBROUTINE read_ipcc
  ! ------------------------------------------------------------------------
  
  ! ------------------------------------------------------------------------
  SUBROUTINE read_ipcc_air(status, fname, var,  data_file)

    USE netcdf

    IMPLICIT NONE

    INTRINSIC :: TRIM, ADJUSTL, NINT

    ! I/O
    INTEGER, INTENT(OUT)                   :: status
    CHARACTER(LEN=*),          INTENT(IN)  :: fname   ! filename
    CHARACTER(LEN=*),          INTENT(IN)  :: var     ! variable name
    REAL(DP), DIMENSION(:,:,:,:), INTENT(OUT):: data_file    ! INTENT(OUT)
    
    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'read_ipcc'
    INTEGER,SAVE                :: ncid   ! netCDF-ID
    INTEGER                     :: dimid, varid

    CHARACTER(LEN=str_vlong)    :: name_dim   ! line

    ! OPEN FILE
    CALL NFERR( status, &
         nf90_open(TRIM(fname), NF90_NOWRITE, ncid) &
         ,21)

    CALL  NFERR( status, &
         nf90_inq_varid(ncid, var, varid ) &
         ,22)
    IF (status.ne.0) then
        write(*,*) "variable not found in NetCDF file, skipping"
        RETURN
    ENDIF

    CALL  NFERR( status, &
         nf90_get_var(ncid, varid, data_file ) &
         ,23)

    !CLOSE FILE
    CALL NFERR( status, &
         nf90_close(ncid) &
         ,24)

    ! RETURN
    status = 0
    
  END SUBROUTINE read_ipcc_air
  ! ------------------------------------------------------------------------
  SUBROUTINE interp(emisflux_class, nyears_added, years_location)
    
   INTEGER, DIMENSION(:), INTENT(IN) :: years_location 
   INTEGER,               INTENT(IN) :: nyears_added
   REAL(DP), DIMENSION(:,:,:), INTENT(INOUT) :: emisflux_class

   SELECT CASE (TRIM(INTERPOLATION))
   CASE ('linear')
       IF (nyears_added.le.1) then
         WRITE (*,*) "NO MULTIPLE YEAR: NO INTERPOLATION NEEDED!"
       ELSE
         SELECT CASE (TRIM(TIME_FREQUENCY))
         CASE ('monthly')
            ALLOCATE (slope(nlon,nlat))
              DO jt=1,ntime
                month_offset=MOD(jt-1,12)
                ! search closest added years
                !DY/DX= (Y(x+1)-Y(x))/Dx
                DO jl=1,nyears_added-1
                  IF (years_location(jl).le.jt.and.years_location(jl+1).gt.jt) THEN
                  slope(:,:)=(emisflux_class(:,:,years_location(jl+1)+month_offset) &
                             -emisflux_class(:,:,years_location(jl)+month_offset)) &
                             / (years_location(jl+1)-years_location(jl))
                  emisflux_class(:,:,jt) = emisflux_class(:,:,years_location(jl)+month_offset)+ &
                        slope(:,:)*(jt-years_location(jl))
                  ENDIF
                ENDDO
              ENDDO
            DEALLOCATE (slope)
         CASE ('annual')
            ALLOCATE (slope(nlon,nlat))
              DO jt=1,ntime
                 ! search closest added years
                 !DY/DX= (Y(x+1)-Y(x))/Dx
                DO jl=1,nyears_added-1
                  IF (years_location(jl).le.jt.and.years_location(jl+1).ge.jt) THEN
                    slope=(emisflux_class(:,:,years_location(jl+1))-emisflux_class(:,:,years_location(jl))) &
                          / (years_location(jl+1)-years_location(jl))
                    emisflux_class(:,:,jt) = emisflux_class(:,:,years_location(jl))+slope(:,:)*(jt-years_location(jl))
                  ENDIF
                ENDDO
              ENDDO
            DEALLOCATE (slope)
         END SELECT
       ENDIF
   CASE ('spline')
       !INTERPOLATE THE EMISSION CLASS
       IF (nyears_added.le.1) then
         WRITE (*,*) "NO MULTIPLE YEAR: NO INTERPOLATION NEEDED!"
       ELSE
         SELECT CASE (TRIM(TIME_FREQUENCY))
         CASE ('monthly')
            ALLOCATE (x_vector(nyears_added))
            ALLOCATE (y_vector(nlon,nlat,nyears_added))
            ALLOCATE (spline_inout(nlon,nlat,nyears_added))
            ALLOCATE (yval(nlon,nlat))
             DO month_offset=0,11
              DO jt=1,nyears_added
                x_vector(jt)=years_location(jt)+month_offset
                y_vector(:,:,jt)=emisflux_class(:,:,years_location(jt)+month_offset)
              ENDDO
              CALL spline_interpolation_base(nyears_added,x_vector(:),y_vector(:,:,:), &
                          spline_inout(:,:,:)) 
              DO jt=1+month_offset,(ntime-12+1+month_offset),12
                xval=jt
                CALL spline_cubic_val(nyears_added,x_vector(:),y_vector(:,:,:),spline_inout(:,:,:),xval,yval(:,:))
                emisflux_class(:,:,jt) = yval(:,:)
              ENDDO
             ENDDO
            DEALLOCATE (yval)
            DEALLOCATE (x_vector)
            DEALLOCATE (y_vector)
            DEALLOCATE (spline_inout)
         CASE ('annual')
            ALLOCATE (x_vector(nyears_added))
            ALLOCATE (y_vector(nlon,nlat,nyears_added))
            ALLOCATE (spline_inout(nlon,nlat,nyears_added))
            ALLOCATE (yval(nlon,nlat))
              DO jt=1,nyears_added
                x_vector(jt)=years_location(jt)
                y_vector(:,:,jt)=emisflux_class(:,:,years_location(jt))
              ENDDO
              CALL spline_interpolation_base(nyears_added,x_vector(:),y_vector(:,:,:), &
                          spline_inout(:,:,:)) 
              DO jt=1,ntime
                xval=jt
                CALL spline_cubic_val(nyears_added,x_vector(:),y_vector(:,:,:),spline_inout(:,:,:),xval,yval(:,:))
                emisflux_class(:,:,jt) = yval(:,:)
              ENDDO
            DEALLOCATE (yval)
            DEALLOCATE (x_vector)
            DEALLOCATE (y_vector)
            DEALLOCATE (spline_inout)
         END SELECT
       ENDIF
   CASE DEFAULT
       WRITE(*,*) 'No valid interpolation method! Check namelist'
       STOP
   END SELECT

  
  END SUBROUTINE

  ! ------------------------------------------------------------------------
  SUBROUTINE nc_dump

    USE netcdf

    IMPLICIT NONE

    INTRINSIC :: DATE_AND_TIME, CHAR

    ! LOCAL
    CHARACTER(LEN=*), PARAMETER :: substr = 'nc_dump'
    INTEGER :: ncid      ! netCDF-ID
    INTEGER :: dimid_lat, dimid_lon, dimid_lev, dimid_ilev, dimid_time
    INTEGER :: varid_lat, varid_lon, varid_lev, varid_ilev, varid_time
    INTEGER :: varid_flux, varid_height, varid_press, varid_ipress
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
    WRITE(yrstr_start,'(i4)') YEAR_START
    WRITE(yrstr_end,'(i4)') YEAR_END
    CALL DATE_AND_TIME(date, time, zone)
    WRITE(mmassstr,*)  MOLARMASS
    WRITE(globssstr,*)  GLOBALSCALE
    IF (L_AIRCRAFT) THEN
      WRITE(pressstr,*) PRESS(1:nlev)
      WRITE(ipressstr,*) IPRESS(1:nlev+1)
    ENDIF
    WRITE(heightstr,*) HEIGHT(1:nlev)

    ! CREATE NEW FILE
    CALL NFERR(status, &
         nf90_create(TRIM(OUTPUT), NF90_CLOBBER, ncid) &
         ,51)

    ! ADD GLOBALE ATTRIBUTES
    ! - VERSION
    CALL NFERR(status, &
         nf90_put_att(ncid, NF90_GLOBAL, 'created_with',     &
         'emcre Version '//TRIM(VERSION)) &
         ,52)
    ! - DATE AND TIME
    CALL NFERR(status, &
         nf90_put_att(ncid, NF90_GLOBAL, 'date', date) &
         ,53)
    CALL NFERR(status, &
         nf90_put_att(ncid, NF90_GLOBAL, 'time', TRIM(time)//TRIM(zone)) &
         ,54)

    ! - NAMELIST
    nmlstr = CHAR(10)//'OUTPUT= '''//TRIM(OUTPUT)//''', '//CHAR(10)
    nmlstr = TRIM(nmlstr)//'SPECIES= '''//TRIM(SPECIES)//''', '//CHAR(10)
    nmlstr = TRIM(nmlstr)//'TIME FREQUENCY= '''//TRIM(TIME_FREQUENCY)//''', '//CHAR(10)
    nmlstr = TRIM(nmlstr)//'OUT UNIT= '''//OUT_UNIT//''', '//CHAR(10)
    nmlstr = TRIM(nmlstr)//'GLOBAL_SCALE= '''//TRIM(globssstr)//''', '//CHAR(10)
    nmlstr = TRIM(nmlstr)//'INTERPOLATION= '''//TRIM(INTERPOLATION)//''', '//CHAR(10)
    nmlstr = TRIM(nmlstr)//'MOLARMASS= '''//TRIM(mmassstr)//''', '//CHAR(10)
    nmlstr = TRIM(nmlstr)//'YEAR_START= '//yrstr_start//', '//CHAR(10)
    nmlstr = TRIM(nmlstr)//'YEAR_END= '//yrstr_end//', '//CHAR(10)
    nmlstr = TRIM(nmlstr)//'HEIGHT= '//TRIM(heightstr)//', '//CHAR(10)
    nmlstr = TRIM(nmlstr)//'INPUTPATH= '''//TRIM(INPUTPATH)//''', '//CHAR(10)
    class_loop: DO jc = 1, MAX_NCLASS
    
      ! SKIP IF NAME IS EMPTY
      IF (TRIM(source(jc)%name) == '') CYCLE

       WRITE(jcstr,'(i4)') jc
       nmlstr=TRIM(nmlstr)//'SOURCE('//jcstr//') : '''//', '//CHAR(10)

       DO jk=1,nlev
         WRITE(levstr,'(i4)') jk
         WRITE(scalestr,'(f6.3)') FRAC(jc)%scale(jk)
         nmlstr = TRIM(nmlstr)//'FRAC('//levstr//')='''//TRIM(scalestr)//''', '//CHAR(10)
       ENDDO

       file_loop: DO jf=1,MAX_SPEC
          ! SKIP IF FILE WAS NOT PRESENT
          IF (YEAR(jc)%year(jf)==0.or.FILE_NAME(jc)%fname(jf)=="") CYCLE 
          WRITE(yrstr,'(i4)') YEAR(jc)%year(jf)
          nmlstr=TRIM(nmlstr)//TRIM(yrstr)//''','//TRIM(FILE_NAME(jc)%fname(jf))&
               &//', '//CHAR(10)
        ENDDO file_loop

    END DO class_loop
    !
    WRITE(*,*) TRIM(nmlstr)
    !
    CALL NFERR(status, &
         nf90_put_att(ncid, NF90_GLOBAL, 'namelist', nmlstr) &
         ,55)

    ! DEFINE DIMENSIONS
    CALL NFERR(status, &
         nf90_def_dim(ncid, 'lon', nlon, dimid_lon) &
         ,56)
    CALL NFERR(status, &
         nf90_def_dim(ncid, 'lat', nlat, dimid_lat) &
         ,57)
    CALL NFERR(status, &
         nf90_def_dim(ncid, 'lev', nlev, dimid_lev) &
         ,58)
    IF (L_AIRCRAFT) THEN
    CALL NFERR(status, &
        nf90_def_dim(ncid, 'ilev', nlev+1, dimid_ilev) &
        , 58)
    ENDIF
    CALL NFERR(status, &
         nf90_def_dim(ncid, 'time', NF90_UNLIMITED, dimid_time) &
         ,59)

    ! DEFINE COORDINATE VARIABLES WITH ATTRIBUTES
    CALL NFERR(status, &
         nf90_def_var(ncid, 'lon', NF90_FLOAT, (/ dimid_lon /), varid_lon) &
         ,60)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_lon, 'long_name', 'longitude') &
         ,61)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_lon, 'units', 'degrees_east') &
         ,62)

    CALL NFERR(status, &
         nf90_def_var(ncid, 'lat', NF90_FLOAT, (/ dimid_lat /), varid_lat) &
         ,63)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_lat, 'long_name', 'latitude') &
         ,63)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_lat, 'units', 'degrees_north') &
         ,64)

    CALL NFERR(status, &
         nf90_def_var(ncid, 'lev', NF90_FLOAT, (/ dimid_lev /), varid_lev) &
         ,65)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_lev, 'long_name', 'level index') &
         ,66)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_lev, 'units', 'level') &
         ,67)
    IF (L_AIRCRAFT) THEN
        ! interface values
       CALL NFERR(status, &
            nf90_def_var(ncid, 'ilev', NF90_FLOAT, (/ dimid_ilev /), varid_ilev) &
            ,15)
       CALL NFERR(status, &
            nf90_put_att(ncid, varid_ilev, 'long_name', 'interface level index') &
            ,16)
       CALL NFERR(status,  &
            nf90_put_att(ncid, varid_ilev, 'units', 'level') &
            ,17)
    ENDIF

    CALL NFERR(status, &
         nf90_def_var(ncid, 'time', NF90_FLOAT, (/ dimid_time /), varid_time) &
         ,68)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_time, 'long_name', 'time') &
         ,69)
    !
    SELECT CASE(TRIM(TIME_FREQUENCY))
    CASE('annual')
       timestr = 'year since '//yrstr_start//'-07-01 00:00:00'
    !CASE('seasonal')
    !   timestr = 'season since '//yrstr//'-01-01 00:00:00'
    CASE('monthly')
       timestr = 'month since '//yrstr_start//'-01-15 00:00:00'
    CASE DEFAULT
       WRITE(*,*) substr,': UNKNOWN DATA FREQUENCY ',TRIM(TIME_FREQUENCY)
       STOP
    END SELECT
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_time, 'units', timestr) &
         ,70)

    ! DEFINE VARIABLES
    IF (L_AIRCRAFT) THEN
       ! - emission pressure mid layer
       CALL NFERR(status, &
            nf90_def_var(ncid, 'press', NF90_FLOAT  &
            , (/ dimid_lev /), varid_press) &
            ,71)
       CALL NFERR(status, &
            nf90_put_att(ncid, varid_press, 'long_name' &
            , 'emission pressure') &
            ,72)
       CALL NFERR(status, &
            nf90_put_att(ncid, varid_press, 'units', 'Pa') &
            ,73)
            ! - emission pressure interface layer
       CALL NFERR(status, &
            nf90_def_var(ncid, 'ipress', NF90_FLOAT  &
            , (/ dimid_ilev /), varid_ipress) &
            ,74)
       CALL NFERR(status, &
            nf90_put_att(ncid, varid_ipress, 'long_name' &
            , 'emission interface pressure') &
            ,75)
       CALL NFERR(status, &
            nf90_put_att(ncid, varid_ipress, 'units', 'Pa') &
            ,76)
    ELSE
        ! - emission height
        CALL NFERR(status, &
             nf90_def_var(ncid, 'height', NF90_FLOAT  &
             , (/ dimid_lev /), varid_height) &
             ,71)
        CALL NFERR(status, &
             nf90_put_att(ncid, varid_height, 'long_name' &
             , 'emission height') &
             ,72)
        CALL NFERR(status, &
             nf90_put_att(ncid, varid_height, 'units', 'm') &
             ,73)
     ENDIF

!    ! - flux
    CALL NFERR(status, &
         nf90_def_var(ncid, TRIM(SPECIES)//'_flux', NF90_FLOAT  &
         , (/ dimid_lon, dimid_lat, dimid_lev, dimid_time /), varid_flux) &
         ,78)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_flux, 'long_name' &
         , 'flux of '//TRIM(SPECIES)) &
         ,79)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_flux, 'units', OUT_UNIT) &
         ,80)
    CALL NFERR(status, &
         nf90_put_att(ncid, varid_flux, 'molar_mass', MOLARMASS) &
         ,81)

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
    CALL NFERR(status, &
         nf90_put_var(ncid, varid_lev, xlev) &
         ,85)
    CALL NFERR(status, &
         nf90_put_var(ncid, varid_time, xtime) &
         ,86)
    IF (L_AIRCRAFT) THEN
         CALL NFERR(status, &
              nf90_put_var(ncid, varid_ilev, xilev ) &
              ,37)
         CALL NFERR(status, &
              nf90_put_var(ncid, varid_press, press(1:nlev)) &
              ,38)
         CALL NFERR(status, &
              nf90_put_var(ncid, varid_ipress, ipress(1:nlev+1)) &
              ,39)
    ELSE
       CALL NFERR(status, &
            nf90_put_var(ncid, varid_height, HEIGHT(1:nlev)) &
            ,87)
    ENDIF

    CALL NFERR(status, &
         nf90_put_var(ncid, varid_flux, emisflux) &
         ,89)

    ! CLOSE FILE
    CALL NFERR(status, &
         nf90_close(ncid) &
         ,90)

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

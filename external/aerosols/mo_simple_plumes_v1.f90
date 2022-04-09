!>
!!
!! @brief Module MO_SIMPLE_PLUMES: provides anthropogenic aerosol optical properties as a function of lat, lon
!!   height, time, and wavelength
!!
!! @remarks
!!
!! @author Bjorn Stevens, Stephanie Fiedler and Karsten Peters MPI-Met, Hamburg (v1 release 2016-11-10)
!!
!! @change-log:
!!          - 2016-12-05: beta release (BS, SF and KP, MPI-Met)
!!          - 2016-09-28: revised representation of Twomey effect (SF, MPI-Met)
!!          - 2015-09-28: bug fixes  (SF, MPI-Met)
!!          - 2016-10-12: revised maximum longitudinal extent of European plume (KP, SF, MPI-Met)
!! $ID: n/a$
!!
!! @par Origin
!!   Based on code originally developed at the MPI-Met by Karsten Peters, Bjorn Stevens, Stephanie Fiedler
!!   and Stefan Kinne with input from Thorsten Mauritsen and Robert Pincus
!!
!! @par Copyright
!!
!
MODULE MO_SIMPLE_PLUMES

  USE netcdf

  IMPLICIT NONE

  integer , parameter :: wp  = selected_real_kind(2*precision(1.0))

  INTEGER, PARAMETER ::                        &
       nplumes   = 9                          ,& !< Number of plumes
       nfeatures = 2                          ,& !< Number of features per plume
       ntimes    = 52                         ,& !< Number of times resolved per year (52 => weekly resolution)
       nyears    = 251                           !< Number of years of available forcing

  LOGICAL, SAVE ::                             &
       sp_initialized = .FALSE.                  !< parameter determining whether input needs to be read

  REAL(kind=wp) ::                             &
       plume_lat      (nplumes)               ,& !< latitude of plume center (AOD maximum)
       plume_lon      (nplumes)               ,& !< longitude of plume center (AOD maximum)
       beta_a         (nplumes)               ,& !< parameter a for beta function vertical profile
       beta_b         (nplumes)               ,& !< parameter b for beta function vertical profile
       aod_spmx       (nplumes)               ,& !< anthropogenic AOD maximum at 550 for plumes
       aod_fmbg       (nplumes)               ,& !< anthropogenic AOD at 550 for fine-mode natural background (idealized to mimic Twomey effect)
       asy550         (nplumes)               ,& !< asymmetry parameter at 550nm for plume
       ssa550         (nplumes)               ,& !< single scattering albedo at 550nm for plume
       angstrom       (nplumes)               ,& !< Angstrom parameter for plume
       sig_lon_E      (nfeatures,nplumes)     ,& !< Eastward extent of plume feature
       sig_lon_W      (nfeatures,nplumes)     ,& !< Westward extent of plume feature
       sig_lat_E      (nfeatures,nplumes)     ,& !< Southward extent of plume feature
       sig_lat_W      (nfeatures,nplumes)     ,& !< Northward extent of plume feature
       theta          (nfeatures,nplumes)     ,& !< Rotation angle of plume feature
       ftr_weight     (nfeatures,nplumes)     ,& !< Feature weights
       time_weight    (nfeatures,nplumes)     ,& !< Time weights
       time_weight_bg (nfeatures,nplumes)     ,& !< as time_weight but for natural background in Twomey effect
       year_weight    (nyears,nplumes)        ,& !< Yearly weight for plume
       ann_cycle      (nfeatures,ntimes,nplumes) !< annual cycle for plume feature

  PUBLIC sp_aop_profile

CONTAINS
  !
  ! ------------------------------------------------------------------------------------------------------------------------
  ! SP_SETUP:  This subroutine should be called at initialization to read the netcdf data that describes the simple plume
  ! climatology.  The information needs to be either read by each processor or distributed to processors.
  !
  SUBROUTINE sp_setup(history,scenario)
    IMPLICIT NONE
    CHARACTER(LEN=*) , INTENT(IN) :: history , scenario
    !
    ! ----------
    !
    INTEGER :: iret, ncid, DimID, VarID, xdmy
    REAL(kind=wp) :: temp(nyears,nplumes)
    !
    ! ----------
    !
    iret = nf90_open(history, NF90_NOWRITE, ncid)
    IF (iret /= NF90_NOERR) THEN
      write (0,*) 'Cannot find file '//trim(history)
      STOP 'NetCDF File not opened'
    END IF
    !
    ! read dimensions and make sure file conforms to expected size
    !
    iret = nf90_inq_dimid(ncid, "plume_number"  , DimId)
    iret = nf90_inquire_dimension(ncid, DimId, len = xdmy)
    IF (xdmy /= nplumes) STOP 'NetCDF improperly dimensioned -- plume_number'

    iret = nf90_inq_dimid(ncid, "plume_feature", DimId)
    iret = nf90_inquire_dimension(ncid, DimId, len = xdmy)
    IF (xdmy /= nfeatures) STOP 'NetCDF improperly dimensioned -- plume_feature'

    iret = nf90_inq_dimid(ncid, "year_fr"   , DimId)
    iret = nf90_inquire_dimension(ncid, DimID, len = xdmy)
    IF (xdmy /= ntimes) STOP 'NetCDF improperly dimensioned -- year_fr'

    iret = nf90_inq_dimid(ncid, "years"   , DimId)
    iret = nf90_inquire_dimension(ncid, DimID, len = xdmy)
    IF (xdmy /= nyears) STOP 'NetCDF improperly dimensioned -- years'
    !
    ! read variables that define the simple plume climatology
    !
    iret = nf90_inq_varid(ncid, "plume_lat", VarId)
    iret = nf90_get_var(ncid, VarID, plume_lat(:), start=(/1/),count=(/nplumes/))
    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading plume_lat'
    iret = nf90_inq_varid(ncid, "plume_lon", VarId)
    iret = nf90_get_var(ncid, VarID, plume_lon(:), start=(/1/),count=(/nplumes/))
    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading plume_lon'
    iret = nf90_inq_varid(ncid, "beta_a"   , VarId)
    iret = nf90_get_var(ncid, VarID, beta_a(:)   , start=(/1/),count=(/nplumes/))
    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading beta_a'
    iret = nf90_inq_varid(ncid, "beta_b"   , VarId)
    iret = nf90_get_var(ncid, VarID, beta_b(:)   , start=(/1/),count=(/nplumes/))
    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading beta_b'
    iret = nf90_inq_varid(ncid, "aod_spmx" , VarId)
    iret = nf90_get_var(ncid, VarID, aod_spmx(:)  , start=(/1/),count=(/nplumes/))
    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading aod_spmx'
    iret = nf90_inq_varid(ncid, "aod_fmbg" , VarId)
    iret = nf90_get_var(ncid, VarID, aod_fmbg(:)  , start=(/1/),count=(/nplumes/))
    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading aod_fmbg'
    iret = nf90_inq_varid(ncid, "ssa550"   , VarId)
    iret = nf90_get_var(ncid, VarID, ssa550(:)  , start=(/1/),count=(/nplumes/))
    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading ssa550'
    iret = nf90_inq_varid(ncid, "asy550"   , VarId)
    iret = nf90_get_var(ncid, VarID, asy550(:)  , start=(/1/),count=(/nplumes/))
    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading asy550'
    iret = nf90_inq_varid(ncid, "angstrom" , VarId)
    iret = nf90_get_var(ncid, VarID, angstrom(:), start=(/1/),count=(/nplumes/))
    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading angstrom'

    iret = nf90_inq_varid(ncid, "sig_lat_W"     , VarId)
    iret = nf90_get_var(ncid, VarID, sig_lat_W(:,:)    , start=(/1,1/),count=(/nfeatures,nplumes/))
    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading sig_lat_W'
    iret = nf90_inq_varid(ncid, "sig_lat_E"     , VarId)
    iret = nf90_get_var(ncid, VarID, sig_lat_E(:,:)    , start=(/1,1/),count=(/nfeatures,nplumes/))
    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading sig_lat_E'
    iret = nf90_inq_varid(ncid, "sig_lon_E"     , VarId)
    iret = nf90_get_var(ncid, VarID, sig_lon_E(:,:)    , start=(/1,1/),count=(/nfeatures,nplumes/))
    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading sig_lon_E'
    iret = nf90_inq_varid(ncid, "sig_lon_W"     , VarId)
    iret = nf90_get_var(ncid, VarID, sig_lon_W(:,:)    , start=(/1,1/),count=(/nfeatures,nplumes/))
    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading sig_lon_W'
    iret = nf90_inq_varid(ncid, "theta"         , VarId)
    iret = nf90_get_var(ncid, VarID, theta(:,:)        , start=(/1,1/),count=(/nfeatures,nplumes/))
    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading theta'
    iret = nf90_inq_varid(ncid, "ftr_weight"    , VarId)
    iret = nf90_get_var(ncid, VarID, ftr_weight(:,:)   , start=(/1,1/),count=(/nfeatures,nplumes/))
    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading plume_lat'
    iret = nf90_inq_varid(ncid, "year_weight"   , VarId)
    iret = nf90_get_var(ncid, VarID, temp, start=(/1,1/),count=(/nyears,nplumes/))
    year_weight(1:165,1:nplumes) = temp(1:165,1:nplumes)
    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading year_weight'
    iret = nf90_inq_varid(ncid, "ann_cycle"     , VarId)
    iret = nf90_get_var(ncid, VarID, ann_cycle(:,:,:)  , start=(/1,1,1/),count=(/nfeatures,ntimes,nplumes/))
    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading ann_cycle'

    iret = nf90_close(ncid)

    iret = nf90_open(scenario, NF90_NOWRITE, ncid)
    IF (iret /= NF90_NOERR) THEN
      write (0,*) 'Cannot find file '//trim(scenario)
      STOP 'NetCDF File not opened'
    END IF
    iret = nf90_inq_varid(ncid, "year_weight"   , VarId)
    iret = nf90_get_var(ncid, VarID, temp, start=(/1,1/),count=(/nyears,nplumes/))
    IF (iret /= NF90_NOERR) STOP 'NetCDF Error reading year_weight'
    year_weight(166:nyears,1:nplumes) = temp(166:nyears,1:nplumes)

    iret = nf90_close(ncid)

    !
    sp_initialized = .TRUE.

    RETURN
  END SUBROUTINE sp_setup
  !
  ! ------------------------------------------------------------------------------------------------------------------------
  ! SET_TIME_WEIGHT:  The simple plume model assumes that meteorology constrains plume shape and that only source strength
  ! influences the amplitude of a plume associated with a given source region.   This routine retrieves the temporal weights
  ! for the plumes.  Each plume feature has its own temporal weights which varies yearly.  The annual cycle is indexed by
  ! week in the year and superimposed on the yearly mean value of the weight.
  !
  SUBROUTINE set_time_weight(year_fr)
    !
    ! ----------
    !
    REAL(kind=wp), INTENT(IN) ::  &
         year_fr           !< Fractional Year (1850.0 - 2100.99)

    INTEGER          ::  &
         iyear          ,& !< Integer year values between 1 and 156 (1850-2100)
         iweek          ,& !< Integer index (between 1 and ntimes); for ntimes=52 this corresponds to weeks (roughly)
         iplume            ! plume number
    !
    ! ----------
    !
    iyear = FLOOR(year_fr) - 1849
    iweek = FLOOR((year_fr - FLOOR(year_fr)) * ntimes) + 1

    IF ((iweek > ntimes) .OR. (iweek < 1) .OR. (iyear > nyears) .OR. (iyear < 1)) STOP 'Time out of bounds in set_time_weight'
    DO iplume=1,nplumes
      time_weight(1,iplume) = year_weight(iyear,iplume) * ann_cycle(1,iweek,iplume)
      time_weight(2,iplume) = year_weight(iyear,iplume) * ann_cycle(2,iweek,iplume)
      time_weight_bg(1,iplume) = ann_cycle(1,iweek,iplume)
      time_weight_bg(2,iplume) = ann_cycle(2,iweek,iplume)
    END DO

    RETURN
  END SUBROUTINE set_time_weight
  !
  ! ------------------------------------------------------------------------------------------------------------------------
  ! SP_AOP_PROFILE:  This subroutine calculates the simple plume aerosol and cloud active optical properties based on the
  ! the simple plume fit to the MPI Aerosol Climatology (Version 2).  It sums over nplumes to provide a profile of aerosol
  ! optical properties on a host models vertical grid.
  !
  SUBROUTINE sp_aop_profile                                                                           ( &
       historic       ,scenario       ,nlevels        ,ncol           , &
       lambda         ,oro            ,lon            ,lat            , &
       year_fr        ,z              ,dz             ,dNovrN         , &
       aod_prof       ,ssa_prof       ,asy_prof       )
    IMPLICIT NONE
    !
    ! ----------
    !
    CHARACTER(LEN=*), INTENT(IN) :: historic, scenario
    INTEGER, INTENT(IN)        :: &
         nlevels,                 & !< number of levels
         ncol                       !< number of columns

    REAL(kind=wp), INTENT(IN)     :: &
         lambda,                  & !< wavelength
         year_fr,                 & !< Fractional Year (1903.0 is the 0Z on the first of January 1903, Gregorian)
         oro(ncol),               & !< orographic height (m)
         lon(ncol),               & !< longitude
         lat(ncol),               & !< latitude
         z (ncol,nlevels),        & !< height above sea-level (m)
         dz(ncol,nlevels)           !< level thickness (difference between half levels) (m)

    REAL(kind=wp), INTENT(OUT) :: &
         dNovrN(ncol)           , & !< anthropogenic increase in cloud drop number concentration (factor)
         aod_prof(ncol,nlevels) , & !< profile of aerosol optical depth
         ssa_prof(ncol,nlevels) , & !< profile of single scattering albedo
         asy_prof(ncol,nlevels)     !< profile of asymmetry parameter

    INTEGER                    :: iplume, icol, k

    REAL(kind=wp)              ::  &
         eta(ncol,nlevels),        & !< normalized height (by 15 km)
         z_beta(ncol,nlevels),     & !< profile for scaling column optical depth
         prof(ncol,nlevels),       & !< scaled profile (by beta function)
         beta_sum(ncol),           & !< vertical sum of beta function
         ssa(ncol),                & !< single scattering albedo
         asy(ncol),                & !< asymmetry parameter
         cw_an(ncol),              & !< column weight for simple plume (anthropogenic) AOD at 550 nm
         cw_bg(ncol),              & !< column weight for fine-mode natural background AOD at 550 nm
         caod_sp(ncol),            & !< column simple plume anthropogenic AOD at 550 nm
         caod_bg(ncol),            & !< column fine-mode natural background AOD at 550 nm
         a_plume1,                 & !< gaussian longitude factor for feature 1
         a_plume2,                 & !< gaussian longitude factor for feature 2
         b_plume1,                 & !< gaussian latitude factor for feature 1
         b_plume2,                 & !< gaussian latitude factor for feature 2
         delta_lat,                & !< latitude offset
         delta_lon,                & !< longitude offset
         delta_lon_t,              & !< threshold for maximum longitudinal plume extent used in transition from 360 to 0 degrees
         lon1,                     & !< rotated longitude for feature 1
         lat1,                     & !< rotated latitude for feature 2
         lon2,                     & !< rotated longitude for feature 1
         lat2,                     & !< rotated latitude for feature 2
         f1,                       & !< contribution from feature 1
         f2,                       & !< contribution from feature 2
         f3,                       & !< contribution from feature 1 in natural background of Twomey effect
         f4,                       & !< contribution from feature 2 in natural background of Twomey effect
         arg,                      & !< exponential "guard"
         aod_550,                  & !< aerosol optical depth at 550nm
         aod_lmd,                  & !< aerosol optical depth at input wavelength
         lfactor                     !< factor to compute wavelength dependence of optical properties
    !
    ! ----------

    !
    ! initialize input data (by calling setup at first instance)
    !
    IF (.NOT.sp_initialized) CALL sp_setup(historic,scenario)
    !
    ! get time weights
    !
    CALL set_time_weight(year_fr)
    !
    ! initialize variables, including output
    !
    DO k=1,nlevels
      DO icol=1,ncol
        aod_prof(icol,k) = 0.0_wp
        ssa_prof(icol,k) = 0.0_wp
        asy_prof(icol,k) = 0.0_wp
!        z_beta(icol,k)   = MERGE(1.0_wp, 0.0_wp, z(icol,k) >= oro(icol))
!FAB in sigma-coordinates the first atm level altitude should always be above
!topography !
        z_beta(icol,k) = 1.!
        eta(icol,k)      = MAX(0.0_wp,MIN(1.0_wp,z(icol,k)/15000.0_wp))
      END DO
    END DO

    DO icol=1,ncol
      dNovrN(icol)   = 1.0_wp
      caod_sp(icol)  = 0.0_wp
      caod_bg(icol)  = 0.02_wp
    END DO
    !
    ! sum contribution from plumes to construct composite profiles of aerosol optical properties
    !
    DO iplume=1,nplumes
      !
      ! calculate vertical distribution function from parameters of beta distribution
      !
      DO icol=1,ncol
        beta_sum(icol) = 0.
      END DO
      DO k=1,nlevels
        DO icol=1,ncol
          prof(icol,k)   = (eta(icol,k)**(beta_a(iplume)-1.) * (1.-eta(icol,k))**(beta_b(iplume)-1.)) * dz(icol,k)
          beta_sum(icol) = beta_sum(icol) + prof(icol,k)
        END DO
      END DO
      DO k=1,nlevels
        DO icol=1,ncol
          prof(icol,k)   = ( prof(icol,k) / beta_sum(icol) ) * z_beta(icol,k)
        END DO
      END DO
      !
      ! calculate plume weights
      !
      DO icol=1,ncol
        !
        ! get plume-center relative spatial parameters for specifying amplitude of plume at given lat and lon
        !
        delta_lat   = lat(icol) - plume_lat(iplume)
        delta_lon   = lon(icol) - plume_lon(iplume)
        delta_lon_t = MERGE (260.0_wp, 180.0_wp, iplume == 1)
        delta_lon   = MERGE ( delta_lon-SIGN(360.0_wp,delta_lon) , delta_lon , ABS(delta_lon) > delta_lon_t)

        a_plume1  = 0.5_wp / (MERGE(sig_lon_E(1,iplume), sig_lon_W(1,iplume), delta_lon > 0)**2)
        b_plume1  = 0.5_wp / (MERGE(sig_lat_E(1,iplume), sig_lat_W(1,iplume), delta_lon > 0)**2)
        a_plume2  = 0.5_wp / (MERGE(sig_lon_E(2,iplume), sig_lon_W(2,iplume), delta_lon > 0)**2)
        b_plume2  = 0.5_wp / (MERGE(sig_lat_E(2,iplume), sig_lat_W(2,iplume), delta_lon > 0)**2)
        !
        ! adjust for a plume specific rotation which helps match plume state to climatology.
        !
        lon1 =   COS(theta(1,iplume))*(delta_lon) + SIN(theta(1,iplume))*(delta_lat)
        lat1 = - SIN(theta(1,iplume))*(delta_lon) + COS(theta(1,iplume))*(delta_lat)
        lon2 =   COS(theta(2,iplume))*(delta_lon) + SIN(theta(2,iplume))*(delta_lat)
        lat2 = - SIN(theta(2,iplume))*(delta_lon) + COS(theta(2,iplume))*(delta_lat)
        !
        ! calculate contribution to plume from its different features, to get a column weight for the anthropogenic
        ! (cw_an) and the fine-mode natural background aerosol (cw_bg)
        !
        arg = (a_plume1 * ((lon1)**2) + (b_plume1 * ((lat1)**2)))
        if ( arg < 25.0_wp ) then
          f1 = time_weight(1,iplume) * ftr_weight(1,iplume) * EXP(-arg)
          f3 = time_weight_bg(1,iplume) * ftr_weight(1,iplume) * EXP(-arg)
        else
          f1 = 0.0_wp
          f3 = 0.0_wp
        end if
        arg = (a_plume2 * ((lon2)**2) + (b_plume2 * ((lat2)**2)))
        if ( arg < 25.0_wp ) then
          f2 = time_weight(2,iplume) * ftr_weight(2,iplume) * EXP(-arg)
          f4 = time_weight_bg(2,iplume) * ftr_weight(2,iplume) * EXP(-arg)
        else
          f2 = 0.0_wp
          f4 = 0.0_wp
        end if
        cw_an(icol) = f1 * aod_spmx(iplume) + f2 * aod_spmx(iplume)
        cw_bg(icol) = f3 * aod_fmbg(iplume) + f4 * aod_fmbg(iplume)
        !
        ! calculate wavelength-dependent scattering properties
        !
        lfactor   = MIN(1.0_wp,700.0_wp/lambda)
        ssa(icol) = (ssa550(iplume) * lfactor**4) / ((ssa550(iplume) * lfactor**4) + ((1-ssa550(iplume)) * lfactor))
        asy(icol) =  asy550(iplume) * SQRT(lfactor)
      END DO
      !
      ! distribute plume optical properties across its vertical profile weighting by optical depth and scaling for
      ! wavelength using the angstrom parameter.
      !
      lfactor = EXP(-angstrom(iplume) * LOG(lambda/550.0_wp))
      !FAB TEST
      DO k=1,nlevels
        DO icol = 1,ncol
          aod_550          = prof(icol,k)     * cw_an(icol)
          aod_lmd          = aod_550          * lfactor
          caod_sp(icol)    = caod_sp(icol)    + aod_550
          caod_bg(icol)    = caod_bg(icol)    + prof(icol,k) * cw_bg(icol)
          asy_prof(icol,k) = asy_prof(icol,k) + aod_lmd * ssa(icol) * asy(icol)
          ssa_prof(icol,k) = ssa_prof(icol,k) + aod_lmd * ssa(icol)
          aod_prof(icol,k) = aod_prof(icol,k) + aod_lmd
        END DO
      END DO
    END DO
    !
    ! complete optical depth weighting
    !
    DO k=1,nlevels
      DO icol = 1,ncol
        asy_prof(icol,k) = MERGE(asy_prof(icol,k)/ssa_prof(icol,k), 0.0_wp, ssa_prof(icol,k) > TINY(1.))
        ssa_prof(icol,k) = MERGE(ssa_prof(icol,k)/aod_prof(icol,k), 1.0_wp, aod_prof(icol,k) > TINY(1.))
      END DO
    END DO
    !
    ! calculate effective radius normalization (divisor) factor
    !
    DO icol=1,ncol
      dNovrN(icol) = LOG((1000.0_wp * (caod_sp(icol) + caod_bg(icol))) + 1.0_wp)/LOG((1000.0_wp * caod_bg(icol)) + 1.0_wp)
    END DO

    RETURN
  END SUBROUTINE sp_aop_profile

END MODULE MO_SIMPLE_PLUMES

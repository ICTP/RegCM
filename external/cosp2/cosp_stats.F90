! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
! Copyright (c) 2015, Regents of the University of Colorado
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without modification, are 
! permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice, this list of 
!    conditions and the following disclaimer.
!
! 2. Redistributions in binary form must reproduce the above copyright notice, this list
!    of conditions and the following disclaimer in the documentation and/or other 
!    materials provided with the distribution.
!
! 3. Neither the name of the copyright holder nor the names of its contributors may be 
!    used to endorse or promote products derived from this software without specific prior
!    written permission.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY 
! EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF 
! MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL 
! THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, 
! SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT 
! OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS 
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT
! LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE
! OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
!
! History:
! Jul 2007 - A. Bodas-Salcedo - Initial version
! Jul 2008 - A. Bodas-Salcedo - Added capability of producing outputs in standard grid
! Oct 2008 - J.-L. Dufresne   - Bug fixed. Assignment of Npoints,Nlevels,Nhydro,Ncolumns 
!                               in COSP_STATS
! Oct 2008 - H. Chepfer       - Added PARASOL reflectance arguments
! Jun 2010 - T. Yokohata, T. Nishimura and K. Ogochi - Added NEC SXs optimisations
! Jan 2013 - G. Cesana        - Added betaperp and temperature arguments 
!                             - Added phase 3D/3Dtemperature/Map output variables in diag_lidar 
! May 2015 - D. Swales        - Modified for cosp2.0 
! Nov 2018 - T. Michibata     - Added CloudSat+MODIS Warmrain Diagnostics
! Jun 2025 - J.K. Shaw.       - Added COSP-RTTOV integration and swathing
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE MOD_COSP_STATS
  USE COSP_KINDS, ONLY: wp
  USE MOD_COSP_CONFIG, &
      ONLY: R_UNDEF,          R_GROUND,         &
            SGCLD_CLR,        SGCLD_ST,         &
            SGCLD_CUM,                          &
            CWP_THRESHOLD,    COT_THRESHOLD,    &
            CFODD_NDBZE,      CFODD_NICOD,      &
            CFODD_BNDRE,      CFODD_BNDZE,      &
            CFODD_NCLASS,                       &
            CFODD_DBZE_MIN,   CFODD_DBZE_MAX,   &
            CFODD_ICOD_MIN,   CFODD_ICOD_MAX,   &
            CFODD_DBZE_WIDTH, CFODD_ICOD_WIDTH, &
            CFODD_HISTDBZE,   CFODD_HISTICOD,   &
            WR_NREGIME,       VGRID_ZU,         &
            VGRID_ZL,         VGRID_Z,          &
            DZ
  USE COSP_PHYS_CONSTANTS,  ONLY: tmelt, radius_earth
  USE COSP_MATH_CONSTANTS, ONLY: pi
  USE MOD_COSP_RTTOV_UTIL, ONLY: rttov_cfg
  IMPLICIT NONE

  ! ######################################################################################
  ! Quickbeam parameters
  integer,parameter :: &
       maxhclass     = 20,  & ! Qucikbeam maximum number of hydrometeor classes.
       nRe_types     = 550, & ! Quickbeam maximum number or Re size bins allowed in N and Z_scaled look up table.
       nd            = 85,  & ! Qucikbeam number of discrete particles used in construction DSDs.
       mt_ntt        = 39,  & ! Quickbeam number of temperatures in mie LUT.
       Re_BIN_LENGTH = 10,  & ! Quickbeam minimum Re interval in scale LUTs  
       Re_MAX_BIN    = 250    ! Quickbeam maximum Re interval in scale LUTs
  real(wp),parameter :: &
       dmin          = 0.1, & ! Quickbeam minimum size of discrete particle
       dmax          = 10000. ! Quickbeam maximum size of discrete particle

  ! ######################################################################################
  ! TYPE cosp_column_inputs
  ! ######################################################################################
  type cosp_column_inputs
     integer :: &
          Npoints,             & ! Number of gridpoints.
          Ncolumns,            & ! Number of columns.
          Nlevels                ! Number of levels.

     integer,allocatable,dimension(:) :: &
          sunlit                 ! Sunlit flag                            (0-1)

     real(wp),allocatable,dimension(:,:) :: &
          at,                  & ! Temperature                            (K)
          pfull,               & ! Pressure                               (Pa)
          phalf,               & ! Pressure at half-levels                (Pa)
          qv,                  & ! Specific humidity                      (kg/kg)
          co2,                 & ! CO2                                    (kg/kg)
          ch4,                 & ! Methane                                (kg/kg)
          n2o,                 & ! N2O                                    (kg/kg)
          co,                  & ! CO                                     (kg/kg)
          so2,                 & ! SO2                                    (kg/kg)
          hgt_matrix,          & ! Height of atmosphere layer             (km)
          hgt_matrix_half        ! Height of bottom interface of atm layer(km)
                                 ! First level contains the bottom of the top layer.
                                 ! Last level contains the bottom of the surface layer.

     real(wp),allocatable,dimension(:) :: &
          land,                & ! Land/Sea mask                          (0 for ocean, 1 for land)
          skt,                 & ! Surface temperature                    (K)
          surfelev               ! Surface Elevation                      (m)
     ! Fields used ONLY by RTTOV
     real(wp),allocatable,dimension(:) :: &
          u_sfc,               & ! Surface u-wind                         (m/s)
          v_sfc,               & ! Surface v-wind                         (m/s)
          t2m,                 & ! 2-meter temperature                    (K)
          q2m,                 & ! 2-meter specific humidity              (kg/kg)
          lat,                 & ! Latitude                               (deg)
          lon,                 & ! Longitude                              (deg)
          sza,                 & ! Solar Zenith Angle in degrees
          psfc,                & ! Surface pressure                       (Pa)
          rttov_sfcmask          ! Mask for RTTOV surface types (0 for ocean, 1 for land, 2 for sea ice)

      real(wp),allocatable,dimension(:,:) :: &
          o3,                  & ! Ozone                                  (kg/kg)
          tca,                 & ! Total layer cloud fraction             (0-1)
          cloudIce,            & ! Cloud ice water mixing ratio           (kg/kg)
          cloudLiq,            & ! Cloud liquid water mixing ratio        (kg/kg)
          DeffLiq,             & ! Cloud liquid effective diameter        (um)
          DeffIce,             & ! Cloud ice effective diameter           (um)
          rttov_date,          & ! Date of the profile as year (e.g. 2013), month (1-12), and day (1-31)
          rttov_time,          & ! Time of profile as hour, minute, second.
          emis_in,             & ! Surface emissivity (point,channel)     (1)
          refl_in,             & ! Surface reflectance (point,channel)    (1)          
          fl_rain,             & ! Precipitation (rain) flux              (kg/m2/s)
          fl_snow                ! Precipitation (snow) flux              (kg/m2/s)
  end type cosp_column_inputs


  ! ######################################################################################
  ! TYPE swath_inputs
  ! ######################################################################################  
  type swath_inputs

     integer ::                            &
          N_inst_swaths = 0
     real(wp),dimension(20) ::      &
          inst_localtimes,                 &
          inst_localtime_widths

  end type swath_inputs

  type radar_cfg
     ! Radar properties
     real(wp) :: freq,k2
     integer  :: nhclass               ! Number of hydrometeor classes in use
     integer  :: use_gas_abs, do_ray
     logical  :: radar_at_layer_one    ! If true radar is assume to be at the edge 
                                       ! of the first layer, if the first layer is the
                                       ! surface than a ground-based radar.   If the
                                       ! first layer is the top-of-atmosphere, then
                                       ! a space borne radar.
     
     ! Variables used to store Z scale factors
     character(len=240)                             :: scale_LUT_file_name
     logical                                        :: load_scale_LUTs, update_scale_LUTs
     logical, dimension(maxhclass,nRe_types)        :: N_scale_flag
     logical, dimension(maxhclass,mt_ntt,nRe_types) :: Z_scale_flag,Z_scale_added_flag
     real(wp),dimension(maxhclass,mt_ntt,nRe_types) :: Ze_scaled,Zr_scaled,kr_scaled
     real(wp),dimension(maxhclass,nd,nRe_types)     :: fc, rho_eff
     real(wp),dimension(Re_MAX_BIN)                 :: base_list,step_list

  end type radar_cfg

  ! ######################################################################################
  ! TYPE cosp_optical_inputs
  ! ######################################################################################
  type cosp_optical_inputs
     integer :: &
          Npoints,             & ! Number of gridpoints.
          Ncolumns,            & ! Number of columns.
          Nlevels,             & ! Number of levels.
          Npart,               & ! Number of cloud meteors for LIDAR simulators.
          Nrefl,               & ! Number of reflectances for PARASOL simulator
          Ninst_rttov            ! Number of RTTOV instruments
     real(wp),pointer :: &
          emis_grey => null()    ! Greybody (spectrally flat) emissivity value for RTTOV          
     real(wp) :: &
          emsfc_lw               ! Surface emissivity @ 11micron
     real(wp),allocatable,dimension(:,:,:) :: &
          frac_out,            & ! Cloud fraction
          tau_067,             & ! Optical depth @ 0.67micron
          emiss_11,            & ! Emissivity @ 11 micron
          fracLiq,             & ! Fraction of optical-depth due to liquid (MODIS)
          asym,                & ! Assymetry parameter @ 3.7micron (MODIS)
          ss_alb,              & ! Single-scattering albedo @ 3.7micron (MODIS)
          betatot_calipso,     & ! Lidar backscatter coefficient (calipso @ 532nm)
          betatot_grLidar532,  & ! Lidar backscatter coefficient (ground-lidar @ 532nm)
          betatot_atlid,       & ! Lidar backscatter coefficient (atlid @ 355nm)
          betatot_ice_calipso, & ! Lidar backscatter coefficient ICE (calipso @ 532nm)
          betatot_liq_calipso, & ! Lidar backscatter coefficient LIQUID (calipso @ 532nm)
          tautot_calipso,      & ! Lidar Optical thickness (calipso @ 532nm)
          tautot_grLidar532,   & ! Lidar Optical thickness (ground-lidar @ 532nm)  
          tautot_atlid,        & ! Lidar Optical thickness (atlid @ 355nm)
          tautot_ice_calipso,  & ! Lidar Ice Optical thickness (calipso @ 532nm)
          tautot_liq_calipso,  & ! Lidar Liquid Optical thickness (calipso @ 532nm)
          z_vol_cloudsat,      & ! Effective reflectivity factor (mm^6/m^3)
          kr_vol_cloudsat,     & ! Attenuation coefficient hydro (dB/km) 
          g_vol_cloudsat         ! Attenuation coefficient gases (dB/km)
     real(wp),allocatable,dimension(:,:) :: &
          beta_mol_calipso,    & ! Lidar molecular backscatter coefficient (calipso @ 532nm)
          beta_mol_grLidar532, & ! Lidar molecular backscatter coefficient (ground-lidar @ 532nm) 
          beta_mol_atlid,      & ! Lidar molecular backscatter coefficient (atlid @ 355nm)
          tau_mol_calipso,     & ! Lidar molecular optical depth (calipso @ 532nm)
          tau_mol_grLidar532,  & ! Lidar molecular optical depth (ground-lidar @ 532nm) 
          tau_mol_atlid,       & ! Lidar molecular optical depth (atlid @ 355nm)
          tautot_S_liq,        & ! Parasol Liquid water optical thickness, from TOA to SFC
          tautot_S_ice,        & ! Parasol Ice water optical thickness, from TOA to SFC
          fracPrecipIce          ! Fraction of precipitation which is frozen (1).
     type(radar_cfg) :: &
          rcfg_cloudsat          ! Radar configuration information (CLOUDSAT)
     type(rttov_cfg),dimension(:),pointer :: &
          cfg_rttov              ! RTTOV configuration information (multiple instruments)
     type(swath_inputs),dimension(6) :: & ! Could be a pointer but fine
          cospswathsIN 
  end type cosp_optical_inputs


CONTAINS
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !---------- SUBROUTINE COSP_CHANGE_VERTICAL_GRID ----------------
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
SUBROUTINE COSP_CHANGE_VERTICAL_GRID(Npoints,Ncolumns,Nlevels,zfull,zhalf,y,Nglevels,newgrid_bot,newgrid_top,r,log_units)
   implicit none
   ! Input arguments
   integer,intent(in) :: Npoints  !# of grid points
   integer,intent(in) :: Nlevels  !# of levels
   integer,intent(in) :: Ncolumns !# of columns
   real(wp),dimension(Npoints,Nlevels),intent(in) :: zfull ! Height at model levels [m] (Bottom of model layer)
   real(wp),dimension(Npoints,Nlevels),intent(in) :: zhalf ! Height at half model levels [m] (Bottom of model layer)
   real(wp),dimension(Npoints,Ncolumns,Nlevels),intent(in) :: y     ! Variable to be changed to a different grid
   integer,intent(in) :: Nglevels  !# levels in the new grid
   real(wp),dimension(Nglevels),intent(in) :: newgrid_bot ! Lower boundary of new levels  [m]
   real(wp),dimension(Nglevels),intent(in) :: newgrid_top ! Upper boundary of new levels  [m]
   logical,optional,intent(in) :: log_units ! log units, need to convert to linear units
   ! Output
   real(wp),dimension(Npoints,Ncolumns,Nglevels),intent(out) :: r ! Variable on new grid

   ! Local variables
   integer :: i,j,k
   logical :: lunits
   integer :: l
   real(wp) :: w ! Weight
   real(wp) :: dbb, dtb, dbt, dtt ! Distances between edges of both grids
   integer :: Nw  ! Number of weights
   real(wp) :: wt  ! Sum of weights
   real(wp),dimension(Nlevels) :: oldgrid_bot,oldgrid_top ! Lower and upper boundaries of model grid
   real(wp) :: yp ! Local copy of y at a particular point.
              ! This allows for change of units.

   lunits=.false.
   if (present(log_units)) lunits=log_units

   r = 0._wp

   do i=1,Npoints
     ! Calculate tops and bottoms of new and old grids
     oldgrid_bot = zhalf(i,:)
     oldgrid_top(1:Nlevels-1) = oldgrid_bot(2:Nlevels)
     oldgrid_top(Nlevels) = zfull(i,Nlevels) +  zfull(i,Nlevels) - zhalf(i,Nlevels) ! Top level symmetric
     l = 0 ! Index of level in the old grid
     ! Loop over levels in the new grid
     do k = 1,Nglevels
       Nw = 0 ! Number of weigths
       wt = 0._wp ! Sum of weights
       ! Loop over levels in the old grid and accumulate total for weighted average
       do
         l = l + 1
         w = 0.0 ! Initialise weight to 0
         if (l > Nlevels) exit
         ! Distances between edges of both grids
         dbb = oldgrid_bot(l) - newgrid_bot(k)
         dtb = oldgrid_top(l) - newgrid_bot(k)
         dbt = oldgrid_bot(l) - newgrid_top(k)
         dtt = oldgrid_top(l) - newgrid_top(k)
         if (dbt >= 0.0) exit ! Do next level in the new grid
         if (dtb > 0.0) then
           if (dbb <= 0.0) then
             if (dtt <= 0) then
               w = dtb
             else
               w = newgrid_top(k) - newgrid_bot(k)
             endif
           else
             if (dtt <= 0) then
               w = oldgrid_top(l) - oldgrid_bot(l)
             else
               w = -dbt
             endif
           endif
           ! If layers overlap (w/=0), then accumulate
           if (w /= 0.0) then
             Nw = Nw + 1
             wt = wt + w
             do j=1,Ncolumns
               if (lunits) then
                 if (y(i,j,l) /= R_UNDEF) then
                   yp = 10._wp**(y(i,j,l)/10._wp)
                 else
                   yp = 0._wp
                 endif
               else
                 yp = y(i,j,l)
               endif
               r(i,j,k) = r(i,j,k) + w*yp
             enddo
           endif
         endif
       enddo
       l = l - 2
       if (l < 1) l = 0
       ! Calculate average in new grid
       if (Nw > 0) then
         do j=1,Ncolumns
           r(i,j,k) = r(i,j,k)/wt
         enddo
       endif
     enddo
   enddo

   ! Set points under surface to R_UNDEF, and change to dBZ if necessary
   do k=1,Nglevels
     do j=1,Ncolumns
       do i=1,Npoints
         if (newgrid_top(k) > zhalf(i,1)) then ! Level above model bottom level
           if (lunits) then
             if (r(i,j,k) <= 0.0) then
               r(i,j,k) = R_UNDEF
             else
               r(i,j,k) = 10._wp*log10(r(i,j,k))
             endif
           endif
         else ! Level below surface
           r(i,j,k) = R_GROUND
         endif
       enddo
     enddo
   enddo

END SUBROUTINE COSP_CHANGE_VERTICAL_GRID

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !------------- SUBROUTINE COSP_LIDAR_ONLY_CLOUD -----------------
  ! (c) 2008, Lawrence Livermore National Security Limited Liability Corporation.
  ! All rights reserved.
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_LIDAR_ONLY_CLOUD(Npoints, Ncolumns, Nlevels, beta_tot, beta_mol,       &
     Ze_tot, lidar_only_freq_cloud, tcc, radar_tcc, radar_tcc2)
    ! Inputs
    integer,intent(in) :: &
         Npoints,       & ! Number of horizontal gridpoints
         Ncolumns,      & ! Number of subcolumns
         Nlevels          ! Number of vertical layers
    real(wp),dimension(Npoints,Nlevels),intent(in) :: &
         beta_mol         ! Molecular backscatter
    real(wp),dimension(Npoints,Ncolumns,Nlevels),intent(in) :: &
         beta_tot,      & ! Total backscattered signal
         Ze_tot           ! Radar reflectivity
    ! Outputs
    real(wp),dimension(Npoints,Nlevels),intent(out) :: &
         lidar_only_freq_cloud
    real(wp),dimension(Npoints),intent(out) ::&
         tcc,       & !
         radar_tcc, & !
         radar_tcc2   !
    
    ! local variables
    real(wp) :: sc_ratio
    real(wp),parameter :: &
         s_cld=5.0, &
         s_att=0.01
    integer :: flag_sat,flag_cld,pr,i,j,flag_radarcld,flag_radarcld_no1km,j_1km
    
    lidar_only_freq_cloud = 0._wp
    tcc = 0._wp
    radar_tcc = 0._wp
    radar_tcc2 = 0._wp
    do pr=1,Npoints
       do i=1,Ncolumns
          flag_sat = 0
          flag_cld = 0
          flag_radarcld = 0 !+JEK
          flag_radarcld_no1km=0 !+JEK
          ! look for j_1km from bottom to top
          j = 1
          do while (Ze_tot(pr,i,j) .eq. R_GROUND)
             j = j+1
          enddo
          j_1km = j+1  !this is the vertical index of 1km above surface  
          
          do j=1,Nlevels
             sc_ratio = beta_tot(pr,i,j)/beta_mol(pr,j)
             if ((sc_ratio .le. s_att) .and. (flag_sat .eq. 0)) flag_sat = j
             if (Ze_tot(pr,i,j) .lt. -30.) then  !radar can't detect cloud
                if ( (sc_ratio .gt. s_cld) .or. (flag_sat .eq. j) ) then  !lidar sense cloud
                   lidar_only_freq_cloud(pr,j)=lidar_only_freq_cloud(pr,j)+1. !top->surf
                   flag_cld=1
                endif
             else  !radar sense cloud (z%Ze_tot(pr,i,j) .ge. -30.)
                flag_cld=1
                flag_radarcld=1
                if (j .gt. j_1km) flag_radarcld_no1km=1              
             endif
          enddo !levels
          if (flag_cld .eq. 1) tcc(pr)=tcc(pr)+1._wp
          if (flag_radarcld .eq. 1) radar_tcc(pr)=radar_tcc(pr)+1.
          if (flag_radarcld_no1km .eq. 1) radar_tcc2(pr)=radar_tcc2(pr)+1.        
       enddo !columns
    enddo !points
    lidar_only_freq_cloud=lidar_only_freq_cloud/Ncolumns
    tcc=tcc/Ncolumns
    radar_tcc=radar_tcc/Ncolumns
    radar_tcc2=radar_tcc2/Ncolumns
    
    ! Unit conversion
    where(lidar_only_freq_cloud /= R_UNDEF) &
            lidar_only_freq_cloud = lidar_only_freq_cloud*100._wp
    where(tcc /= R_UNDEF) tcc = tcc*100._wp
    where(radar_tcc /= R_UNDEF) radar_tcc = radar_tcc*100._wp
    where(radar_tcc2 /= R_UNDEF) radar_tcc2 = radar_tcc2*100._wp
    
  END SUBROUTINE COSP_LIDAR_ONLY_CLOUD
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !-----------------------  SUBROUTINE COSP_DIAG_WARMRAIN ------------------------
  ! (c) 2018, Research Institute for Applied Mechanics (RIAM), Kyushu Univ.
  ! All rights reserved.
  ! * Purpose:    1) Diagnose Contoured Frequency by Optical Depth Diagram (CFODD)
  !                  from CloudSat Radar and MODIS retrievals.
  !               2) Diagnose Warm-Rain Occurrence Frequency (nonprecip/drizzle/rain)
  !                  from CloudSat Radar.
  ! * History:    Jan 2018 (T.Michibata): Test run
  !               May 2018 (T.Michibata): update for COSP2
  !               Sep 2018 (T.Michibata): modified I/O
  !               Nov 2018 (T.Michibata): minor revisions
  !               May 2020 (T.Michibata and X.Jing): bug-fix for frac_out dimsize
  !               Apr 2022 (T.Michibata): bug-fix for non-sunlit columns
  ! * References: Michibata et al. (GMD'19, doi:10.5194/gmd-12-4297-2019)
  !               Michibata et al. (GRL'20, doi:10.1029/2020GL088340)
  !               Suzuki et al. (JAS'10, doi:10.1175/2010JAS3463.1)
  !               Suzuki et al. (JAS'15, doi:10.1175/JAS-D-14-0265.1)
  !               Jing et al. (JCLIM'19, doi:10.1175/jcli-d-18-0789.1)
  ! * Contact:    Takuro Michibata (RIAM, Kyushu University, Japan).
  !               E-mail: michibata@riam.kyushu-u.ac.jp
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_DIAG_WARMRAIN( Npoints, Ncolumns, Nlevels,           & !! in
                                 temp,    zlev,                        & !! in
                                 lwp,     liqcot,   liqreff, liqcfrc,  & !! in
                                 iwp,     icecot,   icereff, icecfrc,  & !! in
                                 fracout, dbze,                        & !! in
                                 cfodd_ntotal,                         & !! inout
                                 wr_occfreq_ntotal                     ) !! inout

    ! Inputs
    integer, intent(in) :: &
         Npoints,          & ! Number of horizontal gridpoints
         Ncolumns,         & ! Number of subcolumns
         Nlevels             ! Number of vertical layers
    real(wp), dimension(Npoints), intent(in) :: &
         lwp,              & ! MODIS LWP [kg m-2]
         liqcot,           & ! MODIS liq. COT
         liqreff,          & ! MODIS liq. Reff [m]
         liqcfrc             ! MODIS liq. cloud fraction
    real(wp), dimension(Npoints), intent(in) :: &
         iwp,              & ! MODIS IWP [kg m-2]
         icecot,           & ! MODIS ice COT
         icereff,          & ! MODIS ice Reff [m]
         icecfrc             ! MODIS ice cloud fraction
    real(wp), dimension(Npoints,Nlevels), intent(in) :: &
         zlev                ! altitude [m] for model level
    real(wp), dimension(Npoints,1,Nlevels),intent(in) :: &
         temp                ! temperature [K]
    real(wp), dimension(Npoints,Ncolumns,Nlevels),intent(in) :: &
         fracout,          & ! SCOPS subcolumn retrieval
         dbze                ! Radar reflectivity [dBZe]

    ! Outputs
    real(wp),dimension(Npoints,CFODD_NDBZE,CFODD_NICOD,CFODD_NCLASS),intent(inout) :: &
         cfodd_ntotal        ! # of CFODD samples for each ICOD/dBZe bin
    real(wp),dimension(Npoints,WR_NREGIME),intent(inout) :: &
         wr_occfreq_ntotal   ! # of SLWC samples

    ! Local variables
    integer  :: i, j, k
    integer  :: kctop, kcbtm
    integer  :: icls
    integer  :: iregime
    real(wp) :: cmxdbz
    real(wp) :: diagcgt   !! diagnosed cloud geometric thickness [m]
    real(wp) :: diagicod  !! diagnosed in-cloud optical depth
    real(wp) :: cbtmh     !! diagnosed in-cloud optical depth
    real(wp), dimension(Npoints,Ncolumns,Nlevels) :: icod  !! in-cloud optical depth (ICOD)
    logical  :: octop, ocbtm, oslwc
    integer, dimension(Npoints,Ncolumns,Nlevels) :: fracout_int  !! fracout (decimal to integer)
    fracout_int(:,:,:) = 0
    where ((fracout(:,:,:) .ge. 0.5_wp) .and. (fracout(:,:,:) .lt. 1.5_wp )) fracout_int(:,:,:) = 1
    where ((fracout(:,:,:) .ge. 1.5_wp) .and. (fracout(:,:,:) .lt. 2.5_wp )) fracout_int(:,:,:) = 2
    !! initialize
    do i = 1, Npoints
       if ( lwp(i) .eq. R_UNDEF ) then  ! for non-sunlit columns
          cfodd_ntotal(i,:,:,:) = R_UNDEF
          wr_occfreq_ntotal(i,:) = R_UNDEF
          icod(i,:,:) = R_UNDEF
       else
          cfodd_ntotal(i,:,:,:)  = 0._wp
          wr_occfreq_ntotal(i,:) = 0._wp
          icod(i,:,:) = 0._wp
       endif
    enddo

    do i = 1, Npoints
       !! check by MODIS retrieval
       if ( ( lwp(i)   .le.  CWP_THRESHOLD  .and. lwp(i) .ne. R_UNDEF  ) .or.  &
          & liqcot(i)  .le.  COT_THRESHOLD  .or.  &
          & liqreff(i) .lt.  CFODD_BNDRE(1) .or.  &
          & liqreff(i) .gt.  CFODD_BNDRE(4) .or.  &
          & iwp(i)     .gt.  CWP_THRESHOLD  .or.  &       !! exclude ice cloud
          & icecot(i)  .gt.  COT_THRESHOLD  .or.  &       !! exclude ice cloud
          & icereff(i) .gt.  CFODD_BNDRE(1)       ) then  !! exclude ice cloud
          cycle
       else
          !! retrieve the CFODD array based on Reff
          icls = 0
          if (    liqreff(i) .ge. CFODD_BNDRE(1) .and. liqreff(i) .lt. CFODD_BNDRE(2) ) then
             icls = 1
          elseif( liqreff(i) .ge. CFODD_BNDRE(2) .and. liqreff(i) .lt. CFODD_BNDRE(3) ) then
             icls = 2
          elseif( liqreff(i) .ge. CFODD_BNDRE(3) .and. liqreff(i) .le. CFODD_BNDRE(4) ) then
             icls = 3
          endif
       endif

       !CDIR NOLOOPCHG
       do j = 1, Ncolumns, 1
          octop = .true.  !! initialize
          ocbtm = .true.  !! initialize
          kcbtm =     0   !! initialize
          kctop =     0   !! initialize
          !CDIR NOLOOPCHG
          do k = Nlevels, 1, -1  !! scan from cloud-bottom to cloud-top
             if ( dbze(i,j,k) .eq. R_GROUND .or. &
                  dbze(i,j,k) .eq. R_UNDEF       ) cycle
             if ( ocbtm                              .and. &
                & fracout_int(i,j,k) .ne. SGCLD_CLR  .and. &
                & dbze(i,j,k)        .ge. CFODD_DBZE_MIN   ) then
                ocbtm = .false.  !! cloud bottom detected
                kcbtm = k
                kctop = k
             endif
             if (       octop                        .and. &  !! scan cloud-top
                & .not. ocbtm                        .and. &  !! cloud-bottom already detected
                & fracout_int(i,j,k) .ne. SGCLD_CLR  .and. &  !! exclude clear sky
                & dbze(i,j,k)        .ge. CFODD_DBZE_MIN   ) then
                kctop = k  !! update
             endif
          enddo  !! k loop
          if( ocbtm ) cycle  !! cloud wasn't detected in this subcolumn
          !! check SLWC?
          if( temp(i,1,kctop) .lt. tmelt ) cycle  !! return to the j (subcolumn) loop
          oslwc = .true.
          cmxdbz = CFODD_DBZE_MIN  !! initialized

          !CDIR NOLOOPCHG
          do k = kcbtm, kctop, -1
             cmxdbz = max( cmxdbz, dbze(i,j,k) )  !! column maximum dBZe update
             if ( fracout_int(i,j,k) .eq. SGCLD_CLR  .or.  &
                & fracout_int(i,j,k) .eq. SGCLD_CUM  .or.  &
                & dbze       (i,j,k) .lt. CFODD_DBZE_MIN   ) then
                oslwc = .false.
             endif
          enddo
          if ( .not. oslwc ) cycle  !! return to the j (subcolumn) loop

          !! warm-rain occurrence frequency
          iregime = 0
          if( cmxdbz .lt. CFODD_BNDZE(1) ) then
             iregime = 1  !! non-precipitating
          elseif( cmxdbz .ge. CFODD_BNDZE(1) .and. cmxdbz .lt. CFODD_BNDZE(2) ) then
             iregime = 2  !! drizzling
          elseif( cmxdbz .ge. CFODD_BNDZE(2) ) then
             iregime = 3  !! raining
          endif
          wr_occfreq_ntotal(i,iregime) = wr_occfreq_ntotal(i,iregime) + 1._wp

          !! retrievals of ICOD and dBZe bin for CFODD plane
          diagcgt = zlev(i,kctop) - zlev(i,kcbtm)
          cbtmh   = zlev(i,kcbtm)
          !CDIR NOLOOPCHG
          do k = kcbtm, kctop, -1
             if( k .eq. kcbtm ) then
                diagicod = liqcot(i)
             else
                diagicod = liqcot(i) * ( 1._wp - ( (zlev(i,k)-cbtmh)/diagcgt)**(5._wp/3._wp) )
             endif
             icod(i,j,k) = min( diagicod, CFODD_ICOD_MAX )
          enddo

       enddo  ! j (Ncolumns)

       !! # of samples for CFODD (joint 2d-histogram dBZe vs ICOD)
       call hist2d( dbze(i,1:Ncolumns,1:Nlevels), icod(i,1:Ncolumns,1:Nlevels), &
                  & Ncolumns*Nlevels,                                           &
                  & CFODD_HISTDBZE, CFODD_NDBZE, CFODD_HISTICOD, CFODD_NICOD,   &
                  & cfodd_ntotal( i, 1:CFODD_NDBZE, 1:CFODD_NICOD, icls )       )

    enddo     ! i (Npoints)

  RETURN
  END SUBROUTINE COSP_DIAG_WARMRAIN

  ! ######################################################################################
  ! FUNCTION hist1D
  ! ######################################################################################
  function hist1d(Npoints,var,nbins,bins)
    ! Inputs
    integer,intent(in) :: &
         Npoints, & ! Number of points in input array
         Nbins      ! Number of bins for sorting
    real(wp),intent(in),dimension(Npoints) :: &
         var        ! Input variable to be sorted
    real(wp),intent(in),dimension(Nbins+1) :: &
         bins       ! Histogram bins [lowest,binTops]  
    ! Outputs
    real(wp),dimension(Nbins) :: &
         hist1d     ! Output histogram      
    ! Local variables
    integer :: ij
    
    do ij=2,Nbins+1  
       hist1D(ij-1) = count(var .ge. bins(ij-1) .and. var .lt. bins(ij))
       if (count(var .eq. R_GROUND) .ge. 1) hist1D(ij-1)=R_UNDEF
    enddo
    
  end function hist1D
  
  ! ######################################################################################
  ! SUBROUTINE hist2D
  ! ######################################################################################
  subroutine hist2D(var1,var2,npts,bin1,nbin1,bin2,nbin2,jointHist)
    implicit none
    
    ! INPUTS
    integer, intent(in) :: &
         npts,  & ! Number of data points to be sorted
         nbin1, & ! Number of bins in histogram direction 1 
         nbin2    ! Number of bins in histogram direction 2
    real(wp),intent(in),dimension(npts) :: &
         var1,  & ! Variable 1 to be sorted into bins
         var2     ! variable 2 to be sorted into bins
    real(wp),intent(in),dimension(nbin1+1) :: &
         bin1     ! Histogram bin 1 boundaries
    real(wp),intent(in),dimension(nbin2+1) :: &
         bin2     ! Histogram bin 2 boundaries
    ! OUTPUTS
    real(wp),intent(out),dimension(nbin1,nbin2) :: &
         jointHist
    
    ! LOCAL VARIABLES
    integer :: ij,ik
    
    do ij=2,nbin1+1
       do ik=2,nbin2+1
          jointHist(ij-1,ik-1)=count(var1 .ge. bin1(ij-1) .and. var1 .lt. bin1(ij) .and. &
               var2 .ge. bin2(ik-1) .and. var2 .lt. bin2(ik))        
       enddo
    enddo
  end subroutine hist2D

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE cosp_cleanUp
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine cosp_cleanUp()  
    deallocate(vgrid_zl,vgrid_zu,vgrid_z,dz)
  end subroutine cosp_cleanUp

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE compute_orbitmasks
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  subroutine compute_orbitmasks(Npoints,Nlocaltimes,localtimes,localtime_widths,      &
                                lat,lon,month,day,hour,minute,swath_mask_out,Nswathed_out)

    ! Inputs
    integer,intent(in) :: &
    Npoints,         &
    Nlocaltimes

    real(wp),dimension(Nlocaltimes),intent(in) :: &
    localtimes,        &
    localtime_widths

    real(wp),dimension(Npoints),intent(in) :: &
    lat,     &
    lon,     &
    month,   &
    day,     &
    hour,    &
    minute

    ! Output
    logical,dimension(Npoints),intent(out) :: &
    swath_mask_out    ! Mask of reals over all gridcells
    integer,intent(out) :: &
    Nswathed_out

    ! Local variables
    integer :: i   ! iterator

    real(wp),dimension(Npoints,Nlocaltimes) :: &
    sat_lon,        & ! Central longitude of the instrument.
    dlon,           & ! distance to satellite longitude in degrees
    dx                ! distance to satellite longitude in km?       

    logical,dimension(Npoints,Nlocaltimes) :: &
    swath_mask_all    ! Mask of logicals over all local times, gridcells  

    integer, dimension(Npoints) :: &
    rttov_DOY         ! Array of day of year values
    real(wp), dimension(Npoints) :: &
    localtime_offsets ! Offset values to avoid striping with hourly RT calls. [hours]
    ! Compute the day of the year and determine the localtime offset
    do i=1,Npoints
      call get_DOY(int(month(i)), int(day(i)), rttov_DOY(i))
    end do
    localtime_offsets = (mod(rttov_DOY(:), 5) - 2) / 5.0  ! Need to cast to real

    ! Iterate over local times
    swath_mask_all(:,:) = .false.
    do i=1,Nlocaltimes
      ! Calculate the central longitude for each gridcell and orbit
      sat_lon(:,i) = 15.0 * (localtimes(i) + localtime_offsets - (hour + minute / 60))
      ! Calculate distance (in degrees) from each grid cell to the satellite central long
      dlon(:,i) = mod((lon - sat_lon(:,i) + 180.0), 360.0) - 180.0             
      ! calculate distance to satellite in km. Remember to convert to radians for cos/sine calls
      dx(:,i)   = dlon(:,i) * (pi/180.0) * COS(lat * pi / 180) * radius_earth
      ! Determine if a gridcell falls in the swath width
      where (abs(dx(:,i))<(localtime_widths(i)*0.5))
        swath_mask_all(:,i) = .true.
      end where        
    end do

    ! Mask is true where values should be calculated
    swath_mask_out = ANY( swath_mask_all(:,:),2) ! Compute mask by collapsing the localtimes dimension
    Nswathed_out   = count(swath_mask_out) ! Number of gridcells that should be calculated.

  end subroutine compute_orbitmasks

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE get_DOY
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
  subroutine get_DOY(month, day, DOY)

   integer,intent(in) :: &
       month,   &
       day
   integer,intent(out) :: &
       DOY

   ! This subroutine does not handle leap years because it is not relevant to the purpose.
   ! Simple look-up table for DOY.
   if (month .eq. 1)  DOY = day
   if (month .eq. 2)  DOY = 31 + day
   if (month .eq. 3)  DOY = 59 + day
   if (month .eq. 4)  DOY = 90 + day
   if (month .eq. 5)  DOY = 120 + day
   if (month .eq. 6)  DOY = 151 + day
   if (month .eq. 7)  DOY = 181 + day
   if (month .eq. 8)  DOY = 212 + day
   if (month .eq. 9)  DOY = 243 + day
   if (month .eq. 10) DOY = 273 + day
   if (month .eq. 11) DOY = 304 + day
   if (month .eq. 12) DOY = 334 + day

 end subroutine get_DOY

END MODULE MOD_COSP_STATS

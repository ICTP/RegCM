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
! May 2015- D. Swales - Original version
! Mar 2018- R. Guzman - Added OPAQ diagnostics and GLID simulator
! Apr 2018- R. Guzman - Added ATLID simulator
! Nov 2018- T. Michibata - Added CloudSat+MODIS Warmrain Diagnostics
! Mar 2024- C. Wall - Added MODIS joint-histogram diagnostics
!
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


MODULE MOD_COSP
  USE COSP_KINDS,                  ONLY: wp
  USE MOD_COSP_CONFIG,             ONLY: R_UNDEF,PARASOL_NREFL,LIDAR_NCAT,LIDAR_NTYPE, SR_BINS,&
                                         N_HYDRO,numMISRHgtBins,      &
                                         cloudsat_DBZE_BINS,LIDAR_NTEMP,calipso_histBsct,&
                                         use_vgrid,Nlvgrid,vgrid_zu,vgrid_zl,vgrid_z,dz, &
                                         WR_NREGIME, CFODD_NCLASS,                       &
                                         CFODD_NDBZE,   CFODD_NICOD,                     &
                                         numMODISTauBins,numMODISPresBins,               &
                                         numMODISReffIceBins,numMODISReffLiqBins,        &
                                         numMODISLWPBins,numMODISIWPBins,                &
                                         numISCCPTauBins,numISCCPPresBins,numMISRTauBins,&
                                         ntau,modis_histTau,tau_binBounds,               &
                                         modis_histTauEdges,tau_binEdges,nCloudsatPrecipClass,&
                                         modis_histTauCenters,tau_binCenters,            &
                                         cloudsat_preclvl,grLidar532_histBsct,atlid_histBsct
  USE MOD_COSP_MODIS_INTERFACE,      ONLY: cosp_modis_init,       modis_IN, &
                                           COSP_ASSIGN_modisIN
  USE MOD_COSP_RTTOV_INTERFACE,      ONLY: cosp_rttov_init,       cosp_rttov_simulate
  USE MOD_COSP_RTTOV_UTIL,           ONLY: rttov_cfg,             rttov_output
  USE MOD_COSP_MISR_INTERFACE,       ONLY: cosp_misr_init,        misr_IN, &
                                           COSP_ASSIGN_misrIN,    COSP_ASSIGN_misrIN_clean
  USE MOD_COSP_ISCCP_INTERFACE,      ONLY: cosp_isccp_init,       isccp_IN, &
                                           COSP_ASSIGN_isccpIN,   COSP_ASSIGN_isccpIN_clean
  USE MOD_COSP_CALIPSO_INTERFACE,    ONLY: cosp_calipso_init,     calipso_IN, &
                                           COSP_ASSIGN_calipsoIN, COSP_ASSIGN_calipsoIN_clean
  USE MOD_COSP_ATLID_INTERFACE,      ONLY: cosp_atlid_init,       atlid_IN, &
                                           COSP_ASSIGN_atlidIN,   COSP_ASSIGN_atlidIN_clean
  USE MOD_COSP_GRLIDAR532_INTERFACE, ONLY: cosp_grLidar532_init,  grLidar532_IN
  USE MOD_COSP_PARASOL_INTERFACE,    ONLY: cosp_parasol_init,     parasol_in, &
                                           COSP_ASSIGN_parasolIN, COSP_ASSIGN_parasolIN_clean
  USE MOD_COSP_CLOUDSAT_INTERFACE,   ONLY: cosp_cloudsat_init,    cloudsat_IN, &
                                           COSP_ASSIGN_cloudsatIN,COSP_ASSIGN_cloudsatIN_clean
  USE quickbeam,                     ONLY: quickbeam_subcolumn,   quickbeam_column
  USE MOD_ICARUS,                    ONLY: icarus_subcolumn,      icarus_column
  USE MOD_MISR_SIMULATOR,            ONLY: misr_subcolumn,        misr_column
  USE MOD_LIDAR_SIMULATOR,           ONLY: lidar_subcolumn,       lidar_column
  USE MOD_MODIS_SIM,                 ONLY: modis_subcolumn,       modis_column
  USE MOD_PARASOL,                   ONLY: parasol_subcolumn,     parasol_column
  USE MOD_COSP_RTTOV,                ONLY: rttov_IN
  USE MOD_COSP_STATS,                ONLY: COSP_LIDAR_ONLY_CLOUD,COSP_CHANGE_VERTICAL_GRID, &
                                           COSP_DIAG_WARMRAIN, COMPUTE_ORBITMASKS,          &
                                           cosp_optical_inputs, cosp_column_inputs,         &
                                           swath_inputs, radar_cfg
  USE COSP_PHYS_CONSTANTS,           ONLY: radius_earth
  USE COSP_MATH_CONSTANTS,           ONLY: pi

  IMPLICIT NONE

  logical :: linitialization ! Initialization flag


  ! ######################################################################################
  ! TYPE cosp_outputs
  ! ######################################################################################
  type cosp_outputs

     ! CALIPSO outputs
     real(wp),dimension(:,:,:),pointer :: &
          calipso_betaperp_tot => null(),  & ! Total backscattered signal
          calipso_beta_tot => null(),      & ! Total backscattered signal
          calipso_lidarcldphase => null(), & ! 3D "lidar" phase cloud fraction
          calipso_lidarcldtype => null(),  & ! 3D "lidar" OPAQ type fraction
          calipso_cldlayerphase => null(), & ! low, mid, high-level lidar phase cloud cover
          calipso_lidarcldtmp => null(),   & ! 3D "lidar" phase cloud temperature
          calipso_cfad_sr => null()          ! CFAD of scattering ratio
     real(wp), dimension(:,:),pointer :: &
          calipso_lidarcld => null(),      & ! 3D "lidar" cloud fraction 
          calipso_cldlayer => null(),      & ! low, mid, high-level, total lidar cloud cover
          calipso_cldtype => null(),       & ! opaque and thin lidar cloud cover + z_opaque altitude
          calipso_cldtypetemp => null(),   & ! opaque and thin cloud temperature  
          calipso_cldtypemeanz => null(),  & ! opaque and thin cloud altitude 
          calipso_cldtypemeanzse => null(),& ! same as just above with respect to SE
          calipso_beta_mol => null()         ! Molecular backscatter
     real(wp), dimension(:),pointer :: &
          calipso_cldthinemis => null(),   & ! thin cloud emissivity 
          calipso_srbval => null()           ! SR bins in cfad_sr

     ! GROUND LIDAR outputs
     real(wp),dimension(:,:,:),pointer :: &  
          grLidar532_beta_tot => null(),   & ! Total GROUND LIDAR backscattered signal 
          grLidar532_cfad_sr => null()       ! CFAD of GROUND LIDAR scattering ratio 
     real(wp), dimension(:,:),pointer :: &  
          grLidar532_lidarcld => null(),   & ! 3D GROUND "lidar" cloud fraction  
          grLidar532_cldlayer => null(),   & ! low, mid, high-level, total GROUND lidar cloud cover
          grLidar532_beta_mol => null()      ! GROUND LIDAR Molecular backscatter 
     real(wp), dimension(:),pointer :: &  
          grLidar532_srbval => null()        ! SR bins in cfad_sr

     ! ATLID outputs
     real(wp),dimension(:,:,:),pointer :: & 
          atlid_beta_tot => null(),   & ! Total ATLID backscattered signal 
          atlid_cfad_sr => null()       ! CFAD of ATLID scattering ratio 
     real(wp), dimension(:,:),pointer :: & 
          atlid_lidarcld => null(),   & ! 3D ATLID cloud fraction
          atlid_cldlayer => null(),   & ! low, mid, high-level, total ATLID cloud cover 
          atlid_beta_mol => null()      ! ATLID Molecular backscatter 
     real(wp), dimension(:),pointer :: & 
          atlid_srbval => null()        ! SR bins in cfad_sr

     ! PARASOL outputs
     real(wp),dimension(:,:,:),pointer :: &
          parasolPix_refl => null()            ! PARASOL reflectances (subcolumn)
     real(wp),dimension(:,:),pointer :: &
          parasolGrid_refl => null()           ! PARASOOL reflectances (column)

     ! CLOUDSAT outputs
     real(wp),dimension(:,:,:),pointer :: &
          cloudsat_Ze_tot => null(),         & ! Effective reflectivity factor (Npoints,Ncolumns,Nlevels)
          cloudsat_cfad_ze => null()           ! Ze CFAD(Npoints,dBZe_bins,Nlevels)
     real(wp), dimension(:,:),pointer :: &
          lidar_only_freq_cloud => null(),   & ! (Npoints,Nlevels)
          cloudsat_precip_cover => null()      ! Radar total cloud amount by CloudSat precip flag (Npoints,dBZe_bins)
     real(wp),dimension(:),pointer :: &
          cloudsat_tcc => null(),             &
          cloudsat_tcc2 => null(),            &          
          radar_lidar_tcc => null(),         & ! Radar&lidar total cloud amount, grid-box scale (Npoints)
          cloudsat_pia => null()               ! Radar path integrated attenuation (Npoints)
          
     ! ISCCP outputs       
     real(wp),dimension(:),pointer :: &
          isccp_totalcldarea => null(), & ! The fraction of model grid box columns with cloud
           		                  ! somewhere in them. (%)
          isccp_meantb => null(),       & ! Mean all-sky 10.5 micron brightness temperature. (K)
          isccp_meantbclr => null(),    & ! Mean clear-sky 10.5 micron brightness temperature. (K)
          isccp_meanptop => null(),     & ! Mean cloud top pressure (mb).
          isccp_meantaucld => null(),   & ! Mean optical thickness. (1)
          isccp_meanalbedocld => null()   ! Mean cloud albedo. (1)
     real(wp),dimension(:,:),pointer ::&
          isccp_boxtau => null(),       & ! Optical thickness in each column. (1)
          isccp_boxptop => null()         ! Cloud top pressure in each column. (mb)
     real(wp),dimension(:,:,:),pointer :: &
          isccp_fq  => null()             ! The fraction of the model grid box covered by each of
                                          ! the 49 ISCCP D level cloud types. (%)

     ! MISR outptus
     real(wp),dimension(:,:,:),pointer ::   & !
          misr_fq => null()          ! Fraction of the model grid box covered by each of the MISR
                           ! cloud types
     real(wp),dimension(:,:),pointer ::   & !
          misr_dist_model_layertops => null() !
     real(wp),dimension(:),pointer ::   & !
          misr_meanztop => null(), & ! Mean MISR cloud top height
          misr_cldarea => null()     ! Mean MISR cloud cover area

     ! MODIS outptus
     real(wp),pointer,dimension(:) ::      & !
          modis_Cloud_Fraction_Total_Mean => null(),       & ! L3 MODIS retrieved cloud fraction (total)
          modis_Cloud_Fraction_Water_Mean => null(),       & ! L3 MODIS retrieved cloud fraction (liq)
          modis_Cloud_Fraction_Ice_Mean => null(),         & ! L3 MODIS retrieved cloud fraction (ice)
          modis_Cloud_Fraction_High_Mean => null(),        & ! L3 MODIS retrieved cloud fraction (high)
          modis_Cloud_Fraction_Mid_Mean => null(),         & ! L3 MODIS retrieved cloud fraction (middle)
          modis_Cloud_Fraction_Low_Mean => null(),         & ! L3 MODIS retrieved cloud fraction (low )
          modis_Optical_Thickness_Total_Mean => null(),    & ! L3 MODIS retrieved optical thickness (tot)
          modis_Optical_Thickness_Water_Mean => null(),    & ! L3 MODIS retrieved optical thickness (liq)
          modis_Optical_Thickness_Ice_Mean => null(),      & ! L3 MODIS retrieved optical thickness (ice)
          modis_Optical_Thickness_Total_LogMean => null(), & ! L3 MODIS retrieved log10 optical thickness
          modis_Optical_Thickness_Water_LogMean => null(), & ! L3 MODIS retrieved log10 optical thickness
          modis_Optical_Thickness_Ice_LogMean => null(),   & ! L3 MODIS retrieved log10 optical thickness
          modis_Cloud_Particle_Size_Water_Mean => null(),  & ! L3 MODIS retrieved particle size (liquid)
          modis_Cloud_Particle_Size_Ice_Mean => null(),    & ! L3 MODIS retrieved particle size (ice)
          modis_Cloud_Top_Pressure_Total_Mean => null(),   & ! L3 MODIS retrieved cloud top pressure
          modis_Liquid_Water_Path_Mean => null(),          & ! L3 MODIS retrieved liquid water path
          modis_Ice_Water_Path_Mean => null()                ! L3 MODIS retrieved ice water path
     real(wp),pointer,dimension(:,:,:) ::  &
          modis_Optical_Thickness_vs_Cloud_Top_Pressure => null(),     & ! Tau/Pressure joint histogram
          modis_Optical_Thickness_vs_Cloud_Top_Pressure_Liq => null(), & ! Tau/Pressure Liq joint histogram
          modis_Optical_Thickness_vs_Cloud_Top_Pressure_Ice => null(), & ! Tau/Pressure Ice joint histogram
          modis_Optical_Thickness_vs_ReffICE => null(),                & ! Tau/ReffICE joint histogram
          modis_Optical_Thickness_vs_ReffLIQ => null(),                & ! Tau/ReffLIQ joint histogram
          modis_LWP_vs_ReffLIQ => null(),                              & ! LWP/ReffLIQ joint histogram
          modis_IWP_vs_ReffICE => null()                                 ! IWP/ReffICE joint histogram

     ! Joint CloudSat+MODIS simulators outputs
     real(wp),dimension(:,:,:,:),pointer :: &
          cfodd_ntotal => null()       ! # of CFODD (Npoints,CFODD_NDBZE,CFODD_NICOD,CFODD_NCLASS)
     real(wp),dimension(:,:),    pointer :: &
          wr_occfreq_ntotal => null()  ! # of nonprecip/drizzle/precip (Npoints,WR_NREGIME)
     integer                    :: &
         Ninst_rttov
     type(rttov_output),dimension(:),allocatable :: &
         rttov_outputs

  end type cosp_outputs

CONTAINS
  ! ######################################################################################
  ! FUNCTION cosp_simulator
  ! ######################################################################################
  function COSP_SIMULATOR(cospIN,cospgridIN,cospOUT,start_idx,stop_idx,debug)
    type(cosp_optical_inputs),intent(in),target :: cospIN     ! Optical inputs to COSP simulator
    type(cosp_column_inputs), intent(in),target :: cospgridIN ! Host model inputs to COSP

    ! Inputs into the simulators
    type(isccp_IN)    :: isccpIN    ! Input to the ISCCP simulator
    type(misr_IN)     :: misrIN     ! Input to the LIDAR simulator
    type(calipso_IN)  :: calipsoIN  ! Input to the LIDAR simulator
    type(grLidar532_IN) :: grLidar532IN ! Input to the GROUND LIDAR simulator
    type(atlid_IN)    :: atlidIN    ! Input to the ATLID simulator 
    type(parasol_IN)  :: parasolIN  ! Input to the PARASOL simulator
    type(cloudsat_IN) :: cloudsatIN ! Input to the CLOUDSAT radar simulator
    type(modis_IN)    :: modisIN    ! Input to the MODIS simulator
    type(rttov_IN)    :: rttovIN    ! Input to the RTTOV simulator

    ! Outputs from the simulators (nested simulator output structure)
    type(cosp_outputs), intent(inout) :: cospOUT
    
    integer,optional  :: start_idx,stop_idx
    logical,optional  :: debug
    
    character(len=256),dimension(100) :: cosp_simulator

    ! Local variables
    integer :: &
         i,icol,ij,ik,nError
    integer :: k
    integer,target :: &
         Npoints
    logical :: &
         Lisccp_subcolumn,     & ! On/Off switch for subcolumn ISCCP simulator
         Lmisr_subcolumn,      & ! On/Off switch for subcolumn MISR simulator
         Lcalipso_subcolumn,   & ! On/Off switch for subcolumn CALIPSO simulator
         LgrLidar532_subcolumn,& ! On/Off switch for subcolumn GROUND LIDAR simulator
         Latlid_subcolumn,     & ! On/Off switch for subcolumn ATLID simulator
         Lparasol_subcolumn,   & ! On/Off switch for subcolumn PARASOL simulator
         Lcloudsat_subcolumn,  & ! On/Off switch for subcolumn CLOUDSAT simulator
         Lmodis_subcolumn,     & ! On/Off switch for subcolumn MODIS simulator
         Lisccp_column,        & ! On/Off switch for column ISCCP simulator
         Lmisr_column,         & ! On/Off switch for column MISR simulator
         Lcalipso_column,      & ! On/Off switch for column CALIPSO simulator
         LgrLidar532_column,   & ! On/Off switch for column GROUND LIDAR simulator
         Latlid_column,        & ! On/Off switch for column ATLID simulator
         Lparasol_column,      & ! On/Off switch for column PARASOL simulator
         Lcloudsat_column,     & ! On/Off switch for column CLOUDSAT simulator
         Lmodis_column,        & ! On/Off switch for column MODIS simulator
         Lrttov_column,        & ! On/Off switch for column RTTOV simulator
         Lradar_lidar_tcc,     & ! On/Off switch from joint Calipso/Cloudsat product
         Lcloudsat_tcc,        & !
         Lcloudsat_tcc2,       & !         
         Llidar_only_freq_cloud, & ! On/Off switch from joint Calipso/Cloudsat product
         Lcloudsat_modis_wr      ! On/Off switch from joint CloudSat/MODIS warm rain product
    logical :: &
         ok_lidar_cfad    = .false.,         &
         ok_lidar_cfad_grLidar532 = .false., & 
         ok_lidar_cfad_atlid = .false.,      &
         verbose = .false.
    integer, dimension(:,:),allocatable  :: &
         modisRetrievedPhase,isccpLEVMATCH
    real(wp), dimension(:),  allocatable  :: &
         modisCfTotal,modisCfLiquid,modisMeanIceWaterPath, isccp_meantbclr,     &
         modisCfIce, modisCfHigh, modisCfMid, modisCfLow,modisMeanTauTotal,     &
         modisMeanTauLiquid, modisMeanTauIce, modisMeanLogTauTotal,             &
         modisMeanLogTauLiquid, modisMeanLogTauIce, modisMeanSizeLiquid,        &
         modisMeanSizeIce, modisMeanCloudTopPressure, modisMeanLiquidWaterPath, &
         radar_lidar_tcc, cloudsat_tcc, cloudsat_tcc2
    REAL(WP), dimension(:,:),allocatable  :: &
         modisRetrievedCloudTopPressure,modisRetrievedTau,modisRetrievedSize,   &
         misr_boxtau,misr_boxztop,misr_dist_model_layertops,isccp_boxtau,       &
         isccp_boxttop,isccp_boxptop,calipso_beta_mol,lidar_only_freq_cloud,    &
         grLidar532_beta_mol,atlid_beta_mol,                                    &
         rttov_bt_total,rttov_bt_clear,                                         & ! RTTOV brightness temps
         rttov_rad_total,rttov_rad_clear,rttov_rad_cloudy,                      & ! RTTOV radiances
         rttov_refl_total,rttov_refl_clear                                        ! RTTOV reflectances
    REAL(WP), dimension(:,:,:),allocatable :: &
         modisJointHistogram,modisJointHistogramIce,modisJointHistogramLiq,     &
         modisJointHistogram_CtpCodLiq,modisJointHistogram_CtpCodIce,           &
         modisJointHistogram_LwpRefLiq,modisJointHistogram_IwpRefIce,           &
         calipso_beta_tot,calipso_betaperp_tot, cloudsatDBZe,parasolPix_refl,   &
         grLidar532_beta_tot,atlid_beta_tot,cloudsatZe_non
    real(wp),dimension(:),allocatable,target :: &
         out1D_1,out1D_2,out1D_3,out1D_4,out1D_5,out1D_6,out1D_7,out1D_8,       &
         out1D_9,out1D_10,out1D_11,out1D_12 
    real(wp),dimension(:,:,:),allocatable :: &
       betamol_in,betamoli,pnormi,ze_toti
    real(wp),dimension(:,:,:),allocatable :: &
         t_in,tempI,frac_outI      ! subscript "I": vertical interpolation (use_vgrid=.true.)
    real(wp), allocatable ::     &
         zlev   (:,:),           & ! altitude (used only when use_vgrid=.true.)
         cfodd_ntotal (:,:,:,:), & ! # of total samples for CFODD (Npoints,CFODD_NDBZE,CFODD_NICOD,CFODD_NCLASS)
         wr_occfreq_ntotal(:,:)    ! # of warm-rain (nonprecip/drizzle/precip) (Npoints,WR_NREGIME)

    ! Fields used in orbit swathing of gridcells.
    logical,dimension(:),allocatable :: & ! Mask of reals over all local times
         CSCAL_SWATH_MASK,    & ! Needed for MODIS CSCAL
         MODIS_SWATH_MASK       ! Needed for MODIS CSCAL + normal MODIS
    integer,dimension(:),allocatable :: & ! Arrays containing the indices of the swath masks
         ISCCP_MASK_INDICES,    &
         MISR_MASK_INDICES,     &
         CSCAL_MASK_INDICES,    &
         PARASOL_MASK_INDICES,  &
         ATLID_MASK_INDICES,    &
         MODIS_CSCAL_MASK_INDICES
    integer        :: j
    ! ISCCP swathing variables  
    real(wp),dimension(:),target,allocatable :: &
         temp_isccp_meanalbedocld,  &
         temp_isccp_meanptop,       &
         temp_isccp_meantaucld,     &
         temp_isccp_totalcldarea,   &
         temp_isccp_meantb
    real(wp),dimension(:,:,:),target,allocatable :: &
         temp_isccp_fq         
    ! MISR swathing variables
    real(wp),dimension(:),target,allocatable :: &
         temp_misr_cldarea,   &
         temp_misr_meanztop
    real(wp),dimension(:,:,:),target,allocatable :: &
         temp_misr_fq
    ! CALIPSO swathing variables
    real(wp),dimension(:),target,allocatable :: &
         temp_calipso_cldthinemis
    real(wp),dimension(:,:),target,allocatable :: &
         temp_calipso_lidarcld,    &
         temp_calipso_cldlayer,    &
         temp_calipso_cldtype,     &
         temp_calipso_cldtypetemp, &
         temp_calipso_cldtypemeanz,&
         temp_calipso_cldtypemeanzse
    real(wp),dimension(:,:,:),target,allocatable :: &
         temp_calipso_cfad_sr,         &
         temp_calipso_lidarcldphase,   &
         temp_calipso_lidarcldtype,    &
         temp_calipso_cldlayerphase,   &
         temp_calipso_lidarcldtmp
    ! ATLID swathing variables
    real(wp),dimension(:,:),target,allocatable :: &
         temp_atlid_lidarcld,          &
         temp_atlid_cldlayer            
    real(wp),dimension(:,:,:),target,allocatable :: &
         temp_atlid_cfad_sr
    ! PARASOL swathing variables
    real(wp),dimension(:,:),target,allocatable :: &
        temp_parasolGrid_refl
    ! CLOUDSAT swathing variables
    real(wp),dimension(:),target,allocatable :: &
        temp_cloudsat_pia
    real(wp),dimension(:,:),target,allocatable :: &
        temp_cloudsat_precip_cover
    real(wp),dimension(:,:,:),target,allocatable :: &
        temp_cloudsat_cfad_ze        
    ! MODIS swathing variables.
    real(wp),dimension(:,:),allocatable :: &
        modis_boxptop,                 &
        modis_boxttop,                 &
        modis_boxtau
    integer,dimension(:,:),allocatable :: &
        modisLEVMATCH
    real(wp),dimension(:),target,allocatable :: &
        modis_meantbclr

    ! Initialize error reporting for output
    cosp_simulator(:)=''
    if (present(debug)) verbose = debug
    
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 1) Determine if using full inputs or subset
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (present(start_idx) .and. present(stop_idx)) then
       ij=start_idx
       ik=stop_idx
    else
       ij=1
       ik=cospIN%Npoints
    endif
    Npoints = ik-ij+1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 2a) Determine which simulators to run and which statistics to compute
    !    - If any of the subcolumn fields are allocated, then run the subcolumn simulators.
    !    - If any of the column fields are allocated, then compute the statistics for that
    !      simulator, but only save the requested fields.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Start with all simulators and joint-diagnostics off
    Lisccp_subcolumn    = .false.
    Lmisr_subcolumn     = .false.
    Lcalipso_subcolumn  = .false.
    LgrLidar532_subcolumn = .false.
    Latlid_subcolumn    = .false. 
    Lparasol_subcolumn  = .false.
    Lcloudsat_subcolumn = .false.
    Lmodis_subcolumn    = .false.
    Lisccp_column       = .false.
    Lmisr_column        = .false.
    Lcalipso_column     = .false.
    LgrLidar532_column = .false.
    Latlid_column       = .false.
    Lparasol_column     = .false.
    Lcloudsat_column    = .false.
    Lmodis_column       = .false.
    Lrttov_column       = .false.
    Lradar_lidar_tcc    = .false.
    Llidar_only_freq_cloud = .false.
    Lcloudsat_tcc       = .false.
    Lcloudsat_tcc2      = .false.
    Lcloudsat_modis_wr  = .false.

    ! CLOUDSAT subcolumn
    if (associated(cospOUT%cloudsat_Ze_tot)) Lcloudsat_subcolumn = .true.

    ! MODIS subcolumn
    if (associated(cospOUT%modis_Cloud_Fraction_Water_Mean)                .or.          &
        associated(cospOUT%modis_Cloud_Fraction_Total_Mean)                .or.          &
        associated(cospOUT%modis_Cloud_Fraction_Ice_Mean)                  .or.          &
        associated(cospOUT%modis_Cloud_Fraction_High_Mean)                 .or.          &
        associated(cospOUT%modis_Cloud_Fraction_Mid_Mean)                  .or.          &
        associated(cospOUT%modis_Cloud_Fraction_Low_Mean)                  .or.          &
        associated(cospOUT%modis_Optical_Thickness_Total_Mean)             .or.          &
        associated(cospOUT%modis_Optical_Thickness_Water_Mean)             .or.          &
        associated(cospOUT%modis_Optical_Thickness_Ice_Mean)               .or.          &
        associated(cospOUT%modis_Optical_Thickness_Total_LogMean)          .or.          &
        associated(cospOUT%modis_Optical_Thickness_Water_LogMean)          .or.          &
        associated(cospOUT%modis_Optical_Thickness_Ice_LogMean)            .or.          &
        associated(cospOUT%modis_Cloud_Particle_Size_Water_Mean)           .or.          &
        associated(cospOUT%modis_Cloud_Particle_Size_Ice_Mean)             .or.          &
        associated(cospOUT%modis_Cloud_Top_Pressure_Total_Mean)            .or.          &
        associated(cospOUT%modis_Liquid_Water_Path_Mean)                   .or.          &
        associated(cospOUT%modis_Ice_Water_Path_Mean)                      .or.          &
        associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure)  .or.          &
        associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Liq) .or.       &
        associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Ice))           &
       Lmodis_subcolumn    = .true.

    ! ISCCP subcolumn
    if (associated(cospOUT%isccp_boxtau)                                   .or.          &
        associated(cospOUT%isccp_boxptop))                                               &
       Lisccp_subcolumn    = .true.

    ! MISR subcolumn
    if (associated(cospOUT%misr_dist_model_layertops))                                   &
       Lmisr_subcolumn     = .true.

    ! CALIPOSO subcolumn
    if (associated(cospOUT%calipso_beta_mol)                               .or.          &
        associated(cospOUT%calipso_betaperp_tot)                           .or.          &
        associated(cospOUT%calipso_beta_tot))                                            &
       Lcalipso_subcolumn  = .true.

    ! GROUND LIDAR subcolumn 
    if (associated(cospOUT%grLidar532_beta_mol)                           .or.          &
        associated(cospOUT%grLidar532_beta_tot))                                        &
       LgrLidar532_subcolumn  = .true. 

    ! ATLID subcolumn
    if (associated(cospOUT%atlid_beta_mol)                                 .or.          & 
        associated(cospOUT%atlid_beta_tot))                                              &
       Latlid_subcolumn  = .true.  

    ! PARASOL subcolumn
    if (associated(cospOUT%parasolPix_refl))                                             &
       Lparasol_subcolumn  = .true.
       
    ! RTTOV Column
    if (allocated(cospOUT%rttov_outputs)) then
       Lrttov_column    = .true.
    endif

    ! ISCCP column
    if (associated(cospOUT%isccp_fq)                                       .or.          &
        associated(cospOUT%isccp_meanalbedocld)                            .or.          &
        associated(cospOUT%isccp_meanptop)                                 .or.          &
        associated(cospOUT%isccp_meantaucld)                               .or.          &
        associated(cospOUT%isccp_totalcldarea)                             .or.          &
        associated(cospOUT%isccp_meantb)) then
       Lisccp_column    = .true.
       Lisccp_subcolumn = .true.
    endif

    ! MISR column
    if (associated(cospOUT%misr_cldarea)                                   .or.          &
        associated(cospOUT%misr_meanztop)                                  .or.          &
        associated(cospOUT%misr_fq)) then
       Lmisr_column    = .true.
       Lmisr_subcolumn = .true.
    endif

    ! CALIPSO column
    if (associated(cospOUT%calipso_cfad_sr)                                .or.          &
        associated(cospOUT%calipso_lidarcld)                               .or.          &
        associated(cospOUT%calipso_lidarcldphase)                          .or.          &
        associated(cospOUT%calipso_lidarcldtype)                           .or.          &
        associated(cospOUT%calipso_cldlayer)                               .or.          &
        associated(cospOUT%calipso_cldtype)                                .or.          & 
        associated(cospOUT%calipso_cldtypetemp)                            .or.          & 
        associated(cospOUT%calipso_cldtypemeanz)                           .or.          & 
        associated(cospOUT%calipso_cldtypemeanzse)                         .or.          & 
        associated(cospOUT%calipso_cldthinemis)                            .or.          &
        associated(cospOUT%calipso_cldlayerphase)                          .or.          &
        associated(cospOUT%calipso_lidarcldtmp)) then
       Lcalipso_column    = .true.
       Lcalipso_subcolumn = .true.
    endif

    ! GROUND LIDAR column 
    if (associated(cospOUT%grLidar532_cfad_sr)                            .or.          & 
        associated(cospOUT%grLidar532_lidarcld)                           .or.          & 
        associated(cospOUT%grLidar532_cldlayer)) then 
       LgrLidar532_column    = .true. 
       LgrLidar532_subcolumn = .true. 
    endif

    ! ATLID column 
    if (associated(cospOUT%atlid_cfad_sr)                                  .or.          & 
        associated(cospOUT%atlid_lidarcld)                                 .or.          & 
        associated(cospOUT%atlid_cldlayer)) then  
       Latlid_column    = .true. 
       Latlid_subcolumn = .true.
    endif

    ! PARASOL column
    if (associated(cospOUT%parasolGrid_refl)) then
       Lparasol_column    = .true.
       Lparasol_subcolumn = .true.
    endif

    ! CLOUDSAT column
    if (associated(cospOUT%cloudsat_cfad_ze)) then
       Lcloudsat_column    = .true.
       Lcloudsat_subcolumn = .true.
    endif

    ! MODIS column
    if (associated(cospOUT%modis_Cloud_Fraction_Total_Mean)                .or.          &
        associated(cospOUT%modis_Cloud_Fraction_Water_Mean)                .or.          &
        associated(cospOUT%modis_Cloud_Fraction_Ice_Mean)                  .or.          &
        associated(cospOUT%modis_Cloud_Fraction_High_Mean)                 .or.          &
        associated(cospOUT%modis_Cloud_Fraction_Mid_Mean)                  .or.          &
        associated(cospOUT%modis_Cloud_Fraction_Low_Mean)                  .or.          &
        associated(cospOUT%modis_Optical_Thickness_Total_Mean)             .or.          &
        associated(cospOUT%modis_Optical_Thickness_Water_Mean)             .or.          &
        associated(cospOUT%modis_Optical_Thickness_Ice_Mean)               .or.          &
        associated(cospOUT%modis_Optical_Thickness_Total_LogMean)          .or.          &
        associated(cospOUT%modis_Optical_Thickness_Water_LogMean)          .or.          &
        associated(cospOUT%modis_Optical_Thickness_Ice_LogMean)            .or.          &
        associated(cospOUT%modis_Cloud_Particle_Size_Water_Mean)           .or.          &
        associated(cospOUT%modis_Cloud_Particle_Size_Ice_Mean)             .or.          &
        associated(cospOUT%modis_Cloud_Top_Pressure_Total_Mean)            .or.          &
        associated(cospOUT%modis_Liquid_Water_Path_Mean)                   .or.          &
        associated(cospOUT%modis_Ice_Water_Path_Mean)                      .or.          &
        associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure)  .or.          &
        associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Liq) .or.       &
        associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Ice)) then
       Lmodis_column    = .true.
       Lmodis_subcolumn = .true.
    endif

    ! Joint simulator products
    if (associated(cospOUT%lidar_only_freq_cloud) .or. associated(cospOUT%radar_lidar_tcc) .or. &
        associated(cospOUT%cloudsat_tcc) .or. associated(cospOUT%cloudsat_tcc2)) then
       Lcalipso_column     = .true.
       Lcalipso_subcolumn  = .true.
       Lcloudsat_column    = .true.
       Lcloudsat_subcolumn = .true.
       Lradar_lidar_tcc    = .true.
       Llidar_only_freq_cloud = .true.
       Lcloudsat_tcc       = .true.
       Lcloudsat_tcc2      = .true.
    endif

    ! CloudSat+MODIS joint simulator product
    if ( associated(cospOUT%cfodd_ntotal) .or. associated(cospOUT%wr_occfreq_ntotal) ) then
       Lmodis_column       = .true.
       Lmodis_subcolumn    = .true.
       Lcloudsat_column    = .true.
       Lcloudsat_subcolumn = .true.
       Lcloudsat_modis_wr  = .true. ! WR: warm rain product
    endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 2b) Error Checking
    !     Enforce bounds on input fields. If input field is out-of-bounds, report error
    !     and turn off simulator
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    call cosp_errorCheck(cospgridIN, cospIN, Lisccp_subcolumn, Lisccp_column,            &
         Lmisr_subcolumn, Lmisr_column, Lmodis_subcolumn, Lmodis_column,                 &
         Lcloudsat_subcolumn, Lcloudsat_column, Lcalipso_subcolumn, Lcalipso_column,     &
         Latlid_subcolumn, Latlid_column, LgrLidar532_subcolumn, LgrLidar532_column,     &
         Lrttov_column, Lparasol_subcolumn, Lparasol_column,                             &
         Lradar_lidar_tcc, Llidar_only_freq_cloud, Lcloudsat_tcc,Lcloudsat_tcc2,         &
         Lcloudsat_modis_wr, cospOUT, cosp_simulator, nError)

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 3) Populate instrument simulator inputs
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! Indexing order for "cospIN % cospswathsIN" is ISCCP, MISR, CLOUDSAT-CALIPSO, ATLID, PARASOL, MODIS
    if (Lisccp_subcolumn .or. Lmodis_subcolumn) then
       call COSP_ASSIGN_isccpIN(cospIN,cospgridIN,Npoints,isccpIN,ISCCP_MASK_INDICES) !COSP_ASSIGN_isccpIN
    endif

    if (Lmisr_subcolumn) then
       call COSP_ASSIGN_misrIN(cospIN,cospgridIN,Npoints,misrIN,MISR_MASK_INDICES)
    endif

    if (Lcalipso_subcolumn) then
       call COSP_ASSIGN_calipsoIN(cospIN,cospgridIN,Npoints,calipsoIN,CSCAL_MASK_INDICES,CSCAL_SWATH_MASK)
    endif

    if (LgrLidar532_subcolumn) then 
       grLidar532IN%Npoints  => Npoints
       grLidar532IN%Ncolumns => cospIN%Ncolumns 
       grLidar532IN%Nlevels  => cospIN%Nlevels
       grLidar532IN%beta_mol => cospIN%beta_mol_grLidar532 
       grLidar532IN%betatot  => cospIN%betatot_grLidar532 
       grLidar532IN%tau_mol  => cospIN%tau_mol_grLidar532 
       grLidar532IN%tautot   => cospIN%tautot_grLidar532  
    endif
    
    if (Latlid_subcolumn) then 
       call COSP_ASSIGN_atlidIN(cospIN,cospgridIN,Npoints,atlidIN,ATLID_MASK_INDICES)
    endif 
    
    if (Lparasol_subcolumn) then
       call COSP_ASSIGN_parasolIN(cospIN,cospgridIN,Npoints,parasolIN,PARASOL_MASK_INDICES)
    endif

    if (Lcloudsat_subcolumn) then
       call COSP_ASSIGN_cloudsatIN(cospIN,cospgridIN,Npoints,cloudsatIN,CSCAL_MASK_INDICES,CSCAL_SWATH_MASK)
    endif

    if (Lmodis_subcolumn) then
       call COSP_ASSIGN_modisIN(cospIN,cospgridIN,Npoints,modisIN,CSCAL_SWATH_MASK,MODIS_CSCAL_MASK_INDICES)
    endif

    if (Lrttov_column) then
       rttovIN%nPoints    => Npoints
       rttovIN%nLevels    => cospIN%nLevels ! This is actually "nlayers" in the RTTOV sense
       rttovIN%nSubCols   => cospIN%nColumns
       rttovIN%co2        => cospgridIN%co2
       rttovIN%ch4        => cospgridIN%ch4
       rttovIN%n2o        => cospgridIN%n2o
       rttovIN%co         => cospgridIN%co
       rttovIN%h_surf     => cospgridIN%surfelev
       rttovIN%u_surf     => cospgridIN%u_sfc
       rttovIN%v_surf     => cospgridIN%v_sfc
       rttovIN%t_skin     => cospgridIN%skt
       rttovIN%p_surf     => cospgridIN%psfc ! Lower boundary of lowest layer may not be the surface.
       if (associated(cospIN%emis_grey)) rttovIN%emis_grey  => cospIN%emis_grey
!       rttovIN%surfem     => cospgridIN%emis_in
!       rttovIN%refl_in     => cospgridIN%refl_in       
       if (allocated(cospgridIN%q2m)) then
          rttovIN%q2m    => cospgridIN%q2m(:)
       else
          rttovIN%q2m    => cospgridIN%qv(:,cospIN%Nlevels)
       end if
       if (allocated(cospgridIN%t2m)) then
          rttovIN%t2m    => cospgridIN%t2m(:)
       else
          rttovIN%t2m    => cospgridIN%at(:,cospIN%Nlevels)
       end if       
       rttovIN%sfcmask    => cospgridIN%rttov_sfcmask
       rttovIN%latitude   => cospgridIN%lat
       rttovIN%longitude  => cospgridIN%lon
       rttovIN%p          => cospgridIN%pfull
       rttovIN%ph         => cospgridIN%phalf
       rttovIN%t          => cospgridIN%at
       rttovIN%q          => cospgridIN%qv
       rttovIN%o3         => cospgridIN%o3
       ! Below only needed for all-sky RTTOV calculation
       rttovIN%rttov_date => cospgridIN%rttov_date
       rttovIN%rttov_time => cospgridIN%rttov_time
       rttovIN%sza        => cospgridIN%sza
       rttovIN%tca        => cospgridIN%tca
       rttovIN%cldLiq     => cospgridIN%cloudLiq
       rttovIN%cldIce     => cospgridIN%cloudIce
       rttovIN%DeffLiq    => cospgridIN%DeffLiq
       rttovIN%DeffIce    => cospgridIN%DeffIce
       rttovIN%fl_rain    => cospgridIN%fl_rain ! Keep in case of RTTOV-SCATT implementation.
       rttovIN%fl_snow    => cospgridIN%fl_snow ! Keep in case of RTTOV-SCATT implementation.
    endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 4) Call subcolumn simulators
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! ISCCP (icarus) subcolumn simulator
    if (Lisccp_subcolumn .or. Lmodis_subcolumn) then
      if (isccpIN%Npoints .gt. 0) then
         allocate(isccpLEVMATCH(isccpIN%Npoints,isccpIN%Ncolumns),  &
                  isccp_boxttop(isccpIN%Npoints,isccpIN%Ncolumns),  &
                  isccp_boxptop(isccpIN%Npoints,isccpIN%Ncolumns),  &
                  isccp_boxtau(isccpIN%Npoints,isccpIN%Ncolumns),   &
                  isccp_meantbclr(isccpIN%Npoints))
            ! Call simulator
         call icarus_subcolumn(isccpIN%npoints,isccpIN%ncolumns,isccpIN%nlevels,           &
                               isccpIN%sunlit,isccpIN%dtau,isccpIN%dem,isccpIN%skt,        &
                               isccpIN%emsfc_lw,isccpIN%qv,isccpIN%at,isccpIN%pfull,       &
                               isccpIN%phalf,isccpIN%frac_out,isccpLEVMATCH,               &
                               isccp_boxtau(:,:),isccp_boxptop(:,:),                       &
                               isccp_boxttop(:,:),isccp_meantbclr(:))
         if (cospIN % cospswathsIN(1) % N_inst_swaths .gt. 0) then
             ! Decode outputs from swaths when reading into cospOUT fields
            cospOUT%isccp_boxtau(ij:ik,:)  = R_UNDEF
            cospOUT%isccp_boxptop(ij:ik,:) = R_UNDEF
            cospOUT%isccp_meantbclr(ij:ik) = R_UNDEF
            if (associated(cospOUT%isccp_boxtau)) cospOUT%isccp_boxtau(ij+int(ISCCP_MASK_INDICES)-1,:) = isccp_boxtau(:,:)
            if (associated(cospOUT%isccp_boxptop)) cospOUT%isccp_boxptop(ij+int(ISCCP_MASK_INDICES)-1,:) = isccp_boxptop(:,:)
            if (associated(cospOUT%isccp_meantbclr)) cospOUT%isccp_meantbclr(ij+int(ISCCP_MASK_INDICES)-1) = isccp_meantbclr(:)
         else
            if (associated(cospOUT%isccp_boxtau))    cospOUT%isccp_boxtau(ij:ik,:)  = isccp_boxtau
            if (associated(cospOUT%isccp_boxptop))   cospOUT%isccp_boxptop(ij:ik,:) = isccp_boxptop
            if (associated(cospOUT%isccp_meantbclr)) cospOUT%isccp_meantbclr(ij:ik) = isccp_meantbclr
         end if
      else
         cospOUT%isccp_boxtau(ij:ik,:) = R_UNDEF
         cospOUT%isccp_boxptop(ij:ik,:) = R_UNDEF
         cospOUT%isccp_meantbclr(ij:ik) = R_UNDEF
      end if
    end if

    ! MISR subcolumn simulator
    if (Lmisr_subcolumn) then
      if (misrIN%Npoints .gt. 0) then
         ! Allocate space for local variables.
         allocate(misr_boxztop(misrIN%Npoints,misrIN%Ncolumns),                                   &
                  misr_boxtau(misrIN%Npoints,misrIN%Ncolumns),                                    &
                  misr_dist_model_layertops(misrIN%Npoints,numMISRHgtBins))     
         ! Call simulator
         call misr_subcolumn(misrIN%Npoints,misrIN%Ncolumns,misrIN%Nlevels,misrIN%dtau,    &
                              misrIN%zfull,misrIN%at,misrIN%sunlit,misr_boxtau,             &
                              misr_dist_model_layertops,misr_boxztop)
         if (cospIN % cospswathsIN(2) % N_inst_swaths .gt. 0) then ! Trigger use of swathed arrays
            ! Decode outputs from swaths when reading into cospOUT fields
            cospOUT%misr_dist_model_layertops(ij:ik,:) = R_UNDEF
            if (associated(cospOUT%misr_dist_model_layertops)) cospOUT%misr_dist_model_layertops(ij+int(MISR_MASK_INDICES)-1,:) = misr_dist_model_layertops(:,:)                             
         else
            if (associated(cospOUT%misr_dist_model_layertops)) cospOUT%misr_dist_model_layertops(ij:ik,:) = misr_dist_model_layertops
         end if
      else
         cospOUT%misr_dist_model_layertops(ij:ik,:) = R_UNDEF
      end if
    end if

    ! Calipso subcolumn simulator
    if (Lcalipso_subcolumn) then
       if (calipsoIN%Npoints .gt. 0) then
          ! Allocate space for local variables
          allocate(calipso_beta_mol(calipsoIN%Npoints,calipsoIN%Nlevels),                   &
                   calipso_beta_tot(calipsoIN%Npoints,calipsoIN%Ncolumns,calipsoIN%Nlevels),&
                   calipso_betaperp_tot(calipsoIN%Npoints,calipsoIN%Ncolumns,calipsoIN%Nlevels))
          ! Call simulator
          call lidar_subcolumn(calipsoIN%npoints, calipsoIN%ncolumns, calipsoIN%nlevels, .false., &
                     calipsoIN%beta_mol, calipsoIN%tau_mol, calipsoIN%betatot, calipsoIN%tautot,  &
                     calipso_beta_mol(:,:), calipso_beta_tot(:,:,:), calipsoIN%betatot_ice,       &
                     calipsoIN%tautot_ice, calipsoIN%betatot_liq, calipsoIN%tautot_liq,           &
                     calipso_betaperp_tot(:,:,:))                      
          if (cospIN % cospswathsIN(3) % N_inst_swaths .gt. 0) then ! Trigger use of swathed arrays
             ! Decode outputs from swaths when reading into cospOUT fields
             cospOUT%calipso_beta_mol(ij:ik,:)       = R_UNDEF
             cospOUT%calipso_beta_tot(ij:ik,:,:)     = R_UNDEF
             cospOUT%calipso_betaperp_tot(ij:ik,:,:) = R_UNDEF
             if (associated(cospOUT%calipso_beta_mol))     cospOUT%calipso_beta_mol(ij+int(CSCAL_MASK_INDICES)-1,calipsoIN%Nlevels:1:-1)     = calipso_beta_mol(:,:)
             if (associated(cospOUT%calipso_beta_tot))     cospOUT%calipso_beta_tot(ij+int(CSCAL_MASK_INDICES)-1,:,calipsoIN%Nlevels:1:-1)   = calipso_beta_tot(:,:,:)
             if (associated(cospOUT%calipso_betaperp_tot)) cospOUT%calipso_betaperp_tot(ij+int(CSCAL_MASK_INDICES)-1,:,:)                    = calipso_betaperp_tot(:,:,:)         
          else ! Proceed normally 
             ! Store output (if requested)
             if (associated(cospOUT%calipso_beta_mol))     cospOUT%calipso_beta_mol(ij:ik,calipsoIN%Nlevels:1:-1) = calipso_beta_mol
             if (associated(cospOUT%calipso_beta_tot))     cospOUT%calipso_beta_tot(ij:ik,:,calipsoIN%Nlevels:1:-1) = calipso_beta_tot
             if (associated(cospOUT%calipso_betaperp_tot)) cospOUT%calipso_betaperp_tot(ij:ik,:,:) = calipso_betaperp_tot
          endif
       else 
          cospOUT%calipso_beta_mol(ij:ik,:)       = R_UNDEF
          cospOUT%calipso_beta_tot(ij:ik,:,:)     = R_UNDEF
          cospOUT%calipso_betaperp_tot(ij:ik,:,:) = R_UNDEF
       endif  
    endif

    ! GROUND LIDAR subcolumn simulator 
    if (LgrLidar532_subcolumn) then  
       ! Allocate space for local variables  
       allocate(grLidar532_beta_mol(grLidar532IN%Npoints,grLidar532IN%Nlevels),       &  
                grLidar532_beta_tot(grLidar532IN%Npoints,grLidar532IN%Ncolumns,grLidar532IN%Nlevels))
       ! Call simulator  
       call lidar_subcolumn(grLidar532IN%npoints, grLidar532IN%ncolumns, grLidar532IN%nlevels,&
            .true., grLidar532IN%beta_mol, grLidar532IN%tau_mol, grLidar532IN%betatot,&
            grLidar532IN%tautot, grLidar532_beta_mol(:,:), grLidar532_beta_tot(:,:,:))
       ! Store output (if requested) 
       if (associated(cospOUT%grLidar532_beta_mol))                                      & 
            cospOUT%grLidar532_beta_mol(ij:ik,grLidar532IN%Nlevels:1:-1) = grLidar532_beta_mol 
       if (associated(cospOUT%grLidar532_beta_tot))                                         & 
            cospOUT%grLidar532_beta_tot(ij:ik,:,grLidar532IN%Nlevels:1:-1) = grLidar532_beta_tot 
    endif 

    ! ATLID subcolumn simulator
    if (Latlid_subcolumn) then
       if (atlidIN%Npoints .gt. 0) then
          ! Allocate space for local variables
          allocate(atlid_beta_mol(atlidIN%Npoints,atlidIN%Nlevels),                      & 
                   atlid_beta_tot(atlidIN%Npoints,atlidIN%Ncolumns,atlidIN%Nlevels)) 
          ! Call simulator 
          call lidar_subcolumn(atlidIN%npoints, atlidIN%ncolumns, atlidIN%nlevels,            &
               .false., atlidIN%beta_mol_atlid, atlidIN%tau_mol_atlid, atlidIN%betatot_atlid, &
               atlidIN%tautot_atlid, atlid_beta_mol(:,:), atlid_beta_tot(:,:,:))
          ! Decode outputs from swaths when reading into cospOUT fields
          if (cospIN % cospswathsIN(4) % N_inst_swaths .gt. 0) then ! Trigger use of swathed arrays
             cospOUT%atlid_beta_mol(ij:ik,:)   = R_UNDEF
             cospOUT%atlid_beta_tot(ij:ik,:,:) = R_UNDEF
             if (associated(cospOUT%atlid_beta_mol))     cospOUT%atlid_beta_mol(ij+int(ATLID_MASK_INDICES)-1,atlidIN%Nlevels:1:-1)     = atlid_beta_mol(:,:)
             if (associated(cospOUT%atlid_beta_tot))     cospOUT%atlid_beta_tot(ij+int(ATLID_MASK_INDICES)-1,:,atlidIN%Nlevels:1:-1)   = atlid_beta_tot(:,:,:)    
          else
             ! Store output (if requested)
             if (associated(cospOUT%atlid_beta_mol))                                        & 
                   cospOUT%atlid_beta_mol(ij:ik,atlidIN%Nlevels:1:-1) = atlid_beta_mol 
             if (associated(cospOUT%atlid_beta_tot))                                        &
                   cospOUT%atlid_beta_tot(ij:ik,:,atlidIN%Nlevels:1:-1) = atlid_beta_tot   
          end if
       else
          cospOUT%atlid_beta_mol(ij:ik,:)   = R_UNDEF
          cospOUT%atlid_beta_tot(ij:ik,:,:) = R_UNDEF
       endif 
    endif

    ! PARASOL subcolumn simulator
    if (Lparasol_subcolumn) then
       if (parasolIN%Npoints .gt. 0) then
          ! Allocate space for local variables
          allocate(parasolPix_refl(parasolIN%Npoints,parasolIN%Ncolumns,PARASOL_NREFL))
          ! Call simulator
          do icol=1,parasolIN%Ncolumns
             call parasol_subcolumn(parasolIN%npoints, PARASOL_NREFL,                       &
                                    parasolIN%tautot_S_liq(1:parasolIN%Npoints,icol),       &
                                    parasolIN%tautot_S_ice(1:parasolIN%Npoints,icol),       &
                                    parasolPix_refl(:,icol,1:PARASOL_NREFL))
             ! Store output (if requested)
             if (cospIN % cospswathsIN(5) % N_inst_swaths .gt. 0) then ! Trigger use of swathed arrays
                cospOUT%parasolPix_refl(ij:ik,icol,1:PARASOL_NREFL) = R_UNDEF
                if (associated(cospOUT%parasolPix_refl))     cospOUT%parasolPix_refl(ij+int(PARASOL_MASK_INDICES)-1,icol,1:PARASOL_NREFL)     = parasolPix_refl(:,icol,1:PARASOL_NREFL)    
             else
                if (associated(cospOUT%parasolPix_refl))     cospOUT%parasolPix_refl(ij:ik,icol,1:PARASOL_NREFL) = parasolPix_refl(:,icol,1:PARASOL_NREFL)
             endif
          enddo
       else
          cospOUT%parasolPix_refl(ij:ik,:,1:PARASOL_NREFL) = R_UNDEF
       endif
    endif

    ! Cloudsat (quickbeam) subcolumn simulator
    if (Lcloudsat_subcolumn) then
       ! Allocate space for local variables
       if (cloudsatIN%Npoints .gt. 0) then
          allocate(cloudsatDBZe(cloudsatIN%Npoints,cloudsatIN%Ncolumns,cloudsatIN%Nlevels), &
                   cloudsatZe_non(cloudsatIN%Npoints,cloudsatIN%Ncolumns,cloudsatIN%Nlevels))
          do icol=1,cloudsatIN%ncolumns
             call quickbeam_subcolumn(cloudsatIN%rcfg,cloudsatIN%Npoints,cloudsatIN%Nlevels,&
                                      cloudsatIN%hgt_matrix/1000._wp,                       &
                                      cloudsatIN%z_vol(:,icol,:),                           &
                                      cloudsatIN%kr_vol(:,icol,:),                          &
                                      cloudsatIN%g_vol(:,1,:),cloudsatDBze(:,icol,:),cloudsatZe_non(:,icol,:))
          enddo
          ! Store output (if requested)
          if (cospIN % cospswathsIN(3) % N_inst_swaths .gt. 0) then ! Trigger use of swathed arrays
             cospOUT%cloudsat_Ze_tot(ij:ik,:,:) = R_UNDEF
             if (associated(cospOUT%cloudsat_Ze_tot))  cospOUT%cloudsat_Ze_tot(ij+int(CSCAL_MASK_INDICES)-1,:,1:cloudsatIN%Nlevels)     = cloudsatDBZe(:,:,1:cloudsatIN%Nlevels)    
          else
             if (associated(cospOUT%cloudsat_Ze_tot))  cospOUT%cloudsat_Ze_tot(ij:ik,:,:) = cloudsatDBZe(:,:,1:cloudsatIN%Nlevels)
          endif
       else
          cospOUT%cloudsat_Ze_tot(ij:ik,:,:) = R_UNDEF
       endif 
    endif

    ! MODIS subcolumn simulator
    if (Lmodis_subcolumn) then
       if (modisiN%nSunlit > 0) then
          ! Allocate space for local variables
          allocate(modisRetrievedTau(modisIN%nSunlit,modisIN%nColumns),                  &
                   modisRetrievedSize(modisIN%nSunlit,modisIN%nColumns),                 &
                   modisRetrievedPhase(modisIN%nSunlit,modisIN%nColumns),                &
                   modisRetrievedCloudTopPressure(modisIN%nSunlit,modisIN%nColumns))
          if ((Lisccp_subcolumn .or. Lmodis_subcolumn) .and. (cospIN % cospswathsIN(1) % N_inst_swaths .gt. 0)) then ! If ISCCP is swathed must run ICARUS separately for MODIS
             allocate(modisLEVMATCH(modisIN%Npoints,modisIN%Ncolumns),     &
                      modis_boxptop(modisIN%Npoints,modisIN%Ncolumns),     &
                      modis_boxttop(modisIN%Npoints,modisIN%Ncolumns),     &
                      modis_boxtau(modisIN%Npoints,modisIN%Ncolumns),      &
                      modis_meantbclr(modisIN%Npoints))
             if (.not. allocated(MODIS_SWATH_MASK)) then ! Allows to run when there is no swathing
                allocate(MODIS_SWATH_MASK(Npoints))
                MODIS_SWATH_MASK(:) = .true.
             end if
             call icarus_subcolumn(modisIN%npoints,modisIN%ncolumns,modisIN%nlevels,                    &
                                   int(MERGE(1,0,MASK=(cospgridIN%sunlit > 0) .and. MODIS_SWATH_MASK)), &
                                   modisIN%tau,cospIN%emiss_11,cospgridIN%skt,                          &
                                   cospIN%emsfc_lw,cospgridIN%qv,cospgridIN%at,cospgridIN%pfull,        &
                                   modisIN%pres,cospIN%frac_out,modisLEVMATCH,                          &    
                                   modis_boxtau(:,:),modis_boxptop(:,:),                                &
                                   modis_boxttop(:,:),modis_meantbclr(:))
             deallocate(modis_boxtau,modis_boxttop,modis_meantbclr,modisLEVMATCH)
          else ! If ISCCP runs fully, just use the cospOUT field indexed appropriately
             allocate(modis_boxptop(modisIN%Npoints,modisIN%Ncolumns))
             modis_boxptop = cospOUT%isccp_boxptop(ij:ik,:)
          end if
          ! Call simulator one column at a time on sunlit columns
          do i = 1, modisIN%nSunlit ! Just run on the sunlit columns, even though the modisIN DDT includes everything orbit swathed
             call modis_subcolumn(modisIN%Ncolumns,modisIN%Nlevels,                      &
                                  modisIN%pres(int(modisIN%sunlit(i)),:),                &
                                  modisIN%tau(int(modisIN%sunlit(i)),:,:),               &
                                  modisIN%liqFrac(int(modisIN%sunlit(i)),:,:),           &
                                  modisIN%g(int(modisIN%sunlit(i)),:,:),                 &
                                  modisIN%w0(int(modisIN%sunlit(i)),:,:),                &
                                  modis_boxptop(int(modisIN%sunlit(i)),:),               &
                                  modisRetrievedPhase(i,:),                              &
                                  modisRetrievedCloudTopPressure(i,:),                   &
                                  modisRetrievedTau(i,:),modisRetrievedSize(i,:))
          end do
          deallocate(modis_boxptop)
       endif
    endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 5) Call column simulators
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! ISCCP
    if (Lisccp_column) then
       ! Check to see which outputs are requested. If not requested, use a local dummy array
       if(.not. associated(cospOUT%isccp_meanalbedocld)) then
          allocate(out1D_1(Npoints))
          cospOUT%isccp_meanalbedocld(ij:ik) => out1D_1
       endif
       if(.not. associated(cospOUT%isccp_meanptop)) then
          allocate(out1D_2(Npoints))
          cospOUT%isccp_meanptop(ij:ik) => out1D_2
       endif
       if(.not. associated(cospOUT%isccp_meantaucld)) then
          allocate(out1D_3(Npoints))
          cospOUT%isccp_meantaucld(ij:ik) => out1D_3
       endif
       if(.not. associated(cospOUT%isccp_totalcldarea)) then
          allocate(out1D_4(Npoints))
          cospOUT%isccp_totalcldarea(ij:ik) => out1D_4
       endif
       if(.not. associated(cospOUT%isccp_meantb)) then
          allocate(out1D_5(Npoints))
          cospOUT%isccp_meantb(ij:ik) => out1D_5
       endif
       if(.not. associated(cospOUT%isccp_fq)) then
          allocate(out1D_6(Npoints*numISCCPTauBins*numISCCPPresBins))
          cospOUT%isccp_fq(ij:ik,1:numISCCPTauBins,1:numISCCPPresBins) => out1D_6
       endif

       ! Call simulator
       if (isccpIN%Npoints .gt. 0) then
          allocate(temp_isccp_fq(isccpIN%Npoints,numISCCPTauBins,numISCCPPresBins),   &
                   temp_isccp_meanalbedocld(isccpIN%Npoints),                         &
                   temp_isccp_meanptop(isccpIN%Npoints),                              &
                   temp_isccp_meantaucld(isccpIN%Npoints),                            &
                   temp_isccp_totalcldarea(isccpIN%Npoints),                          &
                   temp_isccp_meantb(isccpIN%Npoints))
          call icarus_column(isccpIN%npoints, isccpIN%ncolumns,isccp_boxtau(:,:),     &
                             isccp_boxptop(:,:)/100._wp,isccpIN%sunlit,isccp_boxttop(:,:), &
                             temp_isccp_fq,                                           &
                             temp_isccp_meanalbedocld,                                &
                             temp_isccp_meanptop,temp_isccp_meantaucld,               &
                             temp_isccp_totalcldarea,temp_isccp_meantb)
          if (cospIN % cospswathsIN(1) % N_inst_swaths .gt. 0) then ! Trigger use of swathed arrays
             ! Decode back to the cospOUT shapes
             cospOUT%isccp_fq(ij:ik,:,:)        = R_UNDEF
             cospOUT%isccp_meanalbedocld(ij:ik) = R_UNDEF
             cospOUT%isccp_meanptop(ij:ik)      = R_UNDEF
             cospOUT%isccp_meantaucld(ij:ik)    = R_UNDEF
             cospOUT%isccp_totalcldarea(ij:ik)  = R_UNDEF
             cospOUT%isccp_meantb(ij:ik)        = R_UNDEF
             if (associated(cospOUT%isccp_fq))            cospOUT%isccp_fq(ij+int(ISCCP_MASK_INDICES)-1,:,:)        = temp_isccp_fq(:,:,:)
             if (associated(cospOUT%isccp_meanalbedocld)) cospOUT%isccp_meanalbedocld(ij+int(ISCCP_MASK_INDICES)-1) = temp_isccp_meanalbedocld(:)
             if (associated(cospOUT%isccp_meanptop))      cospOUT%isccp_meanptop(ij+int(ISCCP_MASK_INDICES)-1)      = temp_isccp_meanptop(:)
             if (associated(cospOUT%isccp_meantaucld))    cospOUT%isccp_meantaucld(ij+int(ISCCP_MASK_INDICES)-1)    = temp_isccp_meantaucld(:)
             if (associated(cospOUT%isccp_totalcldarea))  cospOUT%isccp_totalcldarea(ij+int(ISCCP_MASK_INDICES)-1)  = temp_isccp_totalcldarea(:)
             if (associated(cospOUT%isccp_meantb))        cospOUT%isccp_meantb(ij+int(ISCCP_MASK_INDICES)-1)        = temp_isccp_meantb(:)
          else
             if (associated(cospOUT%isccp_fq))            cospOUT%isccp_fq(ij:ik,:,:)        = temp_isccp_fq(:,:,:)
             if (associated(cospOUT%isccp_meanalbedocld)) cospOUT%isccp_meanalbedocld(ij:ik) = temp_isccp_meanalbedocld(:)
             if (associated(cospOUT%isccp_meanptop))      cospOUT%isccp_meanptop(ij:ik)      = temp_isccp_meanptop(:)
             if (associated(cospOUT%isccp_meantaucld))    cospOUT%isccp_meantaucld(ij:ik)    = temp_isccp_meantaucld(:)
             if (associated(cospOUT%isccp_totalcldarea))  cospOUT%isccp_totalcldarea(ij:ik)  = temp_isccp_totalcldarea(:)
             if (associated(cospOUT%isccp_meantb))        cospOUT%isccp_meantb(ij:ik)        = temp_isccp_meantb(:)
          end if
          deallocate(temp_isccp_fq,temp_isccp_meanalbedocld,temp_isccp_meanptop,temp_isccp_meantaucld,temp_isccp_totalcldarea,temp_isccp_meantb)
       else
          cospOUT%isccp_fq(ij:ik,:,:)        = R_UNDEF
          cospOUT%isccp_meanalbedocld(ij:ik) = R_UNDEF
          cospOUT%isccp_meanptop(ij:ik)      = R_UNDEF
          cospOUT%isccp_meantaucld(ij:ik)    = R_UNDEF
          cospOUT%isccp_totalcldarea(ij:ik)  = R_UNDEF
          cospOUT%isccp_meantb(ij:ik)        = R_UNDEF
       end if
       cospOUT%isccp_fq(ij:ik,:,:) = cospOUT%isccp_fq(ij:ik,:,7:1:-1)

       ! Check if there is any value slightly greater than 1
       where ((cospOUT%isccp_totalcldarea > 1.0-1.e-5) .and.                             &
              (cospOUT%isccp_totalcldarea < 1.0+1.e-5))
              cospOUT%isccp_totalcldarea = 1.0
       endwhere

       ! Clear up memory (if necessary)
       if (allocated(isccp_boxttop))   deallocate(isccp_boxttop)
       if (allocated(isccp_boxptop))   deallocate(isccp_boxptop)
       if (allocated(isccp_boxtau))    deallocate(isccp_boxtau)
       if (allocated(isccp_meantbclr)) deallocate(isccp_meantbclr)
       if (allocated(isccpLEVMATCH))   deallocate(isccpLEVMATCH)
       if (allocated(ISCCP_MASK_INDICES)) deallocate(ISCCP_MASK_INDICES)
       if (allocated(out1D_1)) then
          deallocate(out1D_1)
          nullify(cospOUT%isccp_meanalbedocld)
       endif
       if (allocated(out1D_2)) then
          deallocate(out1D_2)
          nullify(cospOUT%isccp_meanptop)
       endif
       if (allocated(out1D_3)) then
          deallocate(out1D_3)
          nullify(cospOUT%isccp_meantaucld)
       endif
       if (allocated(out1D_4)) then
          deallocate(out1D_4)
          nullify(cospOUT%isccp_totalcldarea)
       endif
       if (allocated(out1D_5)) then
          deallocate(out1D_5)
          nullify(cospOUT%isccp_meantb)
       endif
       if (allocated(out1D_6)) then
          deallocate(out1D_6)
          nullify(cospOUT%isccp_fq)
       endif
    endif

    ! MISR
    if (Lmisr_column) then
       ! Check to see which outputs are requested. If not requested, use a local dummy array
       if (.not. associated(cospOUT%misr_cldarea)) then
          allocate(out1D_1(Npoints))
          cospOUT%misr_cldarea(ij:ik) => out1D_1
       endif
       if (.not. associated(cospOUT%misr_meanztop)) then
          allocate(out1D_2(Npoints))
          cospOUT%misr_meanztop(ij:ik) => out1D_2
       endif
       if (.not. associated(cospOUT%misr_fq)) then
          allocate(out1D_3(Npoints*numMISRTauBins*numMISRHgtBins))
          cospOUT%misr_fq(ij:ik,1:numMISRTauBins,1:numMISRHgtBins) => out1D_3
       endif

       ! Call simulator
       if (cospIN % cospswathsIN(2) % N_inst_swaths .gt. 0) then ! Trigger use of swathed arrays
          if (misrIN%Npoints .gt. 0) then
             ! Operate at the masked format.
             allocate(temp_misr_cldarea(misrIN%Npoints),                          &
                      temp_misr_meanztop(misrIN%Npoints),                         &
                      temp_misr_fq(misrIN%Npoints,numMISRTauBins,numMISRHgtBins))
             call misr_column(misrIN%Npoints,misrIN%Ncolumns,misr_boxztop,misrIN%sunlit, &
                              misr_boxtau,temp_misr_cldarea(:),                          &
                              temp_misr_meanztop(:),temp_misr_fq(:,:,:))                  
             ! Decode back to the cospOUT shapes
             cospOUT%misr_cldarea(ij:ik) = R_UNDEF
             cospOUT%misr_meanztop(ij:ik) = R_UNDEF
             cospOUT%misr_fq(ij:ik,:,:) = R_UNDEF
             if (associated(cospOUT%misr_cldarea))    cospOUT%misr_cldarea(ij+int(MISR_MASK_INDICES)-1)    = temp_misr_cldarea(:)
             if (associated(cospOUT%misr_meanztop))   cospOUT%misr_meanztop(ij+int(MISR_MASK_INDICES)-1)   = temp_misr_meanztop(:)
             if (associated(cospOUT%misr_fq))         cospOUT%misr_fq(ij+int(MISR_MASK_INDICES)-1,:,:)     = temp_misr_fq(:,:,:)          
             deallocate(temp_misr_cldarea,temp_misr_meanztop,temp_misr_fq)
          else
             cospOUT%misr_cldarea(ij:ik) = R_UNDEF
             cospOUT%misr_meanztop(ij:ik) = R_UNDEF
             cospOUT%misr_fq(ij:ik,:,:) = R_UNDEF           
          endif
       else
          call misr_column(misrIN%Npoints,misrIN%Ncolumns,misr_boxztop,misrIN%sunlit,&
                           misr_boxtau,cospOUT%misr_cldarea(ij:ik),                  &
                           cospOUT%misr_meanztop(ij:ik),cospOUT%misr_fq(ij:ik,:,:))
       endif
       ! Clear up memory
       if (allocated(misr_boxtau))               deallocate(misr_boxtau)
       if (allocated(misr_boxztop))              deallocate(misr_boxztop)
       if (allocated(misr_dist_model_layertops)) deallocate(misr_dist_model_layertops)
       if (allocated(MISR_MASK_INDICES))         deallocate(MISR_MASK_INDICES)
       if (allocated(out1D_1)) then
          deallocate(out1D_1)
          nullify(cospOUT%misr_cldarea)
       endif
       if (allocated(out1D_2)) then
          deallocate(out1D_2)
          nullify(cospOUT%misr_meanztop)
       endif
       if (allocated(out1D_3)) then
          deallocate(out1D_3)
          nullify(cospOUT%misr_fq)
       endif
    endif

    ! CALIPSO LIDAR Simulator
    if (Lcalipso_column) then
       ! Check to see which outputs are requested. If not requested, use a local dummy array
       if (.not. associated(cospOUT%calipso_cfad_sr)) then
          allocate(out1D_1(Npoints*SR_BINS*Nlvgrid))
          cospOUT%calipso_cfad_sr(ij:ik,1:SR_BINS,1:Nlvgrid) => out1D_1
       endif
       if (.not. associated(cospOUT%calipso_lidarcld)) then
          allocate(out1D_2(Npoints*Nlvgrid))
          cospOUT%calipso_lidarcld(ij:ik,1:Nlvgrid) => out1D_2
       endif
       if (.not. associated(cospOUT%calipso_lidarcldphase)) then
          allocate(out1D_3(Npoints*Nlvgrid*6))
          cospOUT%calipso_lidarcldphase(ij:ik,1:Nlvgrid,1:6) => out1D_3
       endif
       if (.not. associated(cospOUT%calipso_cldlayer)) then
          allocate(out1D_4(Npoints*LIDAR_NCAT))
          cospOUT%calipso_cldlayer(ij:ik,1:LIDAR_NCAT) => out1D_4
       endif
       if (.not. associated(cospOUT%calipso_cldlayerphase)) then
          allocate(out1D_5(Npoints*LIDAR_NCAT*6))
          cospOUT%calipso_cldlayerphase(ij:ik,1:LIDAR_NCAT,1:6) => out1D_5
       endif
       if (.not. associated(cospOUT%calipso_lidarcldtmp)) then
          allocate(out1D_6(Npoints*40*5))
          cospOUT%calipso_lidarcldtmp(ij:ik,1:40,1:5) => out1D_6
       endif   
       if (.not. associated(cospOUT%calipso_lidarcldtype)) then
          allocate(out1D_7(Npoints*Nlvgrid*4))
          cospOUT%calipso_lidarcldtype(ij:ik,1:Nlvgrid,1:4) => out1D_7
       endif
       if (.not. associated(cospOUT%calipso_cldtype)) then 
          allocate(out1D_8(Npoints*LIDAR_NTYPE)) 
          cospOUT%calipso_cldtype(ij:ik,1:LIDAR_NTYPE) => out1D_8  
       endif 
       if (.not. associated(cospOUT%calipso_cldtypetemp)) then 
          allocate(out1D_9(Npoints*LIDAR_NTYPE))  
          cospOUT%calipso_cldtypetemp(ij:ik,1:LIDAR_NTYPE) => out1D_9 
       endif 
       if (.not. associated(cospOUT%calipso_cldtypemeanz)) then 
          allocate(out1D_10(Npoints*2))
          cospOUT%calipso_cldtypemeanz(ij:ik,1:2) => out1D_10 
       endif 
       if (.not. associated(cospOUT%calipso_cldtypemeanzse)) then
          allocate(out1D_12(Npoints*3)) 
          cospOUT%calipso_cldtypemeanzse(ij:ik,1:3) => out1D_12
       endif
       if (.not. associated(cospOUT%calipso_cldthinemis)) then
          allocate(out1D_11(Npoints)) 
          cospOUT%calipso_cldthinemis(ij:ik) => out1D_11 
       endif
       
       ! Call simulator
       if (cospIN % cospswathsIN(3) % N_inst_swaths .gt. 0) then ! Trigger use of swathed arrays
          ! Operate at the masked format.
          if (calipsoIN%Npoints .gt. 0) then 
             allocate(temp_calipso_cfad_sr(calipsoIN%Npoints,SR_BINS,Nlvgrid),              &
                      temp_calipso_lidarcld(calipsoIN%Npoints,Nlvgrid),                     &
                      temp_calipso_cldlayer(calipsoIN%Npoints,LIDAR_NCAT),                  &
                      temp_calipso_lidarcldphase(calipsoIN%Npoints,Nlvgrid,6),              &
                      temp_calipso_lidarcldtype(calipsoIN%Npoints,Nlvgrid,LIDAR_NTYPE+1),   &
                      temp_calipso_cldtype(calipsoIN%Npoints,LIDAR_NTYPE),                  &
                      temp_calipso_cldtypetemp(calipsoIN%Npoints,LIDAR_NTYPE),              &
                      temp_calipso_cldtypemeanz(calipsoIN%Npoints,2),                       &
                      temp_calipso_cldtypemeanzse(calipsoIN%Npoints,3),                     &
                      temp_calipso_cldthinemis(calipsoIN%Npoints),                          &
                      temp_calipso_cldlayerphase(calipsoIN%Npoints,LIDAR_NCAT,6),           &
                      temp_calipso_lidarcldtmp(calipsoIN%Npoints,LIDAR_NTEMP,5))
             ok_lidar_cfad=.true.
             call lidar_column(calipsoIN%Npoints, calipsoIN%Ncolumns, calipsoIN%Nlevels,                            &
                  Nlvgrid, SR_BINS, LIDAR_NTYPE, 'calipso',calipso_beta_tot(:,:,:), calipso_beta_mol(:,:),          &
                  cospgridIN%phalf(int(CSCAL_MASK_INDICES),2:calipsoIN%Nlevels+1),                                  &
                  cospgridIN%hgt_matrix(int(CSCAL_MASK_INDICES),:),                                                 &
                  cospgridIN%hgt_matrix_half(int(CSCAL_MASK_INDICES),:), vgrid_z(:), ok_lidar_cfad, LIDAR_NCAT,     &
                  temp_calipso_cfad_sr(:,:,:), temp_calipso_lidarcld(:,:), temp_calipso_cldlayer(:,:),              &
                  cospgridIN%at(int(CSCAL_MASK_INDICES),:), calipso_betaperp_tot(:,:,:),                            &
                  cospgridIN%surfelev(int(CSCAL_MASK_INDICES)),                                                     &
                  temp_calipso_lidarcldphase(:,:,:),                             &
                  temp_calipso_lidarcldtype(:,:,:),  temp_calipso_cldtype(:,:),  &
                  temp_calipso_cldtypetemp(:,:), temp_calipso_cldtypemeanz(:,:), & 
                  temp_calipso_cldtypemeanzse(:,:), temp_calipso_cldthinemis(:), &
                  temp_calipso_cldlayerphase(:,:,:), temp_calipso_lidarcldtmp(:,:,:))
             ! Decode back to the cospOUT shapes
             cospOUT%calipso_cfad_sr(ij:ik,:,:)       = R_UNDEF
             cospOUT%calipso_lidarcld(ij:ik,:)        = R_UNDEF
             cospOUT%calipso_cldlayer(ij:ik,:)        = R_UNDEF
             cospOUT%calipso_lidarcldphase(ij:ik,:,:) = R_UNDEF
             cospOUT%calipso_lidarcldtype(ij:ik,:,:)  = R_UNDEF
             cospOUT%calipso_cldtype(ij:ik,:)         = R_UNDEF
             cospOUT%calipso_cldtypetemp(ij:ik,:)     = R_UNDEF
             cospOUT%calipso_cldtypemeanz(ij:ik,:)    = R_UNDEF
             cospOUT%calipso_cldtypemeanzse(ij:ik,:)  = R_UNDEF
             cospOUT%calipso_cldthinemis(ij:ik)       = R_UNDEF
             cospOUT%calipso_cldlayerphase(ij:ik,:,:) = R_UNDEF
             cospOUT%calipso_lidarcldtmp(ij:ik,:,:)   = R_UNDEF
             ! Note: The active simulators use R_UNDEF to indicate unretrieved and clear-sky values.
             ! Decisions on how to treat these values are delegeate to the host user, and fields
             ! are often set to zero.
             ! Using R_UNDEF in this way is incompatible with the swath approach, where R_UNDEF means
             ! that a field should not be used in statistically averaging, etc.
             ! To avoid this error, when swathing we set R_UNDEF value to zero here.
             where (temp_calipso_cfad_sr(:,:,:) == R_UNDEF) temp_calipso_cfad_sr(:,:,:) = 0._wp
             where (temp_calipso_lidarcldphase(:,:,:) == R_UNDEF) temp_calipso_lidarcldphase(:,:,:) = 0._wp
             where (temp_calipso_lidarcld(:,:) == R_UNDEF) temp_calipso_lidarcld(:,:) = 0._wp
             ! Unpack into the full lat-lon structure
             if (associated(cospOUT%calipso_cfad_sr))          cospOUT%calipso_cfad_sr(ij+int(CSCAL_MASK_INDICES)-1,:,:)         = temp_calipso_cfad_sr(:,:,:)
             if (associated(cospOUT%calipso_lidarcld))         cospOUT%calipso_lidarcld(ij+int(CSCAL_MASK_INDICES)-1,:)          = temp_calipso_lidarcld(:,:)
             if (associated(cospOUT%calipso_cldlayer))         cospOUT%calipso_cldlayer(ij+int(CSCAL_MASK_INDICES)-1,:)          = temp_calipso_cldlayer(:,:)
             if (associated(cospOUT%calipso_lidarcldphase))    cospOUT%calipso_lidarcldphase(ij+int(CSCAL_MASK_INDICES)-1,:,:)   = temp_calipso_lidarcldphase(:,:,:)
             if (associated(cospOUT%calipso_lidarcldtype))     cospOUT%calipso_lidarcldtype(ij+int(CSCAL_MASK_INDICES)-1,:,:)    = temp_calipso_lidarcldtype(:,:,:)
             if (associated(cospOUT%calipso_cldtype))          cospOUT%calipso_cldtype(ij+int(CSCAL_MASK_INDICES)-1,:)           = temp_calipso_cldtype(:,:)
             if (associated(cospOUT%calipso_cldtypetemp))      cospOUT%calipso_cldtypetemp(ij+int(CSCAL_MASK_INDICES)-1,:)       = temp_calipso_cldtypetemp(:,:)
             if (associated(cospOUT%calipso_cldtypemeanz))     cospOUT%calipso_cldtypemeanz(ij+int(CSCAL_MASK_INDICES)-1,:)      = temp_calipso_cldtypemeanz(:,:)
             if (associated(cospOUT%calipso_cldtypemeanzse))   cospOUT%calipso_cldtypemeanzse(ij+int(CSCAL_MASK_INDICES)-1,:)    = temp_calipso_cldtypemeanzse(:,:)
             if (associated(cospOUT%calipso_cldthinemis))      cospOUT%calipso_cldthinemis(ij+int(CSCAL_MASK_INDICES)-1)         = temp_calipso_cldthinemis(:)
             if (associated(cospOUT%calipso_cldlayerphase))    cospOUT%calipso_cldlayerphase(ij+int(CSCAL_MASK_INDICES)-1,:,:)   = temp_calipso_cldlayerphase(:,:,:)
             if (associated(cospOUT%calipso_lidarcldtmp))      cospOUT%calipso_lidarcldtmp(ij+int(CSCAL_MASK_INDICES)-1,:,:)     = temp_calipso_lidarcldtmp(:,:,:)
             deallocate(temp_calipso_cfad_sr,temp_calipso_lidarcld,temp_calipso_cldlayer,temp_calipso_lidarcldphase,   &
                        temp_calipso_lidarcldtype,temp_calipso_cldtype,temp_calipso_cldtypetemp,temp_calipso_cldtypemeanz,                &
                        temp_calipso_cldtypemeanzse,temp_calipso_cldthinemis,temp_calipso_cldlayerphase,temp_calipso_lidarcldtmp)
          else
             cospOUT%calipso_cfad_sr(ij:ik,:,:)       = R_UNDEF
             cospOUT%calipso_lidarcld(ij:ik,:)        = R_UNDEF
             cospOUT%calipso_cldlayer(ij:ik,:)        = R_UNDEF
             cospOUT%calipso_lidarcldphase(ij:ik,:,:) = R_UNDEF
             cospOUT%calipso_lidarcldtype(ij:ik,:,:)  = R_UNDEF
             cospOUT%calipso_cldtype(ij:ik,:)         = R_UNDEF
             cospOUT%calipso_cldtypetemp(ij:ik,:)     = R_UNDEF
             cospOUT%calipso_cldtypemeanz(ij:ik,:)    = R_UNDEF
             cospOUT%calipso_cldtypemeanzse(ij:ik,:)  = R_UNDEF
             cospOUT%calipso_cldthinemis(ij:ik)       = R_UNDEF
             cospOUT%calipso_cldlayerphase(ij:ik,:,:) = R_UNDEF
             cospOUT%calipso_lidarcldtmp(ij:ik,:,:)   = R_UNDEF         
          endif
       else    
          ok_lidar_cfad=.true.
          call lidar_column(calipsoIN%Npoints, calipsoIN%Ncolumns, calipsoIN%Nlevels,       &
               Nlvgrid, SR_BINS, LIDAR_NTYPE, 'calipso',calipso_beta_tot(:,:,:), calipso_beta_mol(:,:),&
               cospgridIN%phalf(:,2:calipsoIN%Nlevels+1),cospgridIN%hgt_matrix,             &
               cospgridIN%hgt_matrix_half, vgrid_z(:), ok_lidar_cfad, LIDAR_NCAT,           &
               cospOUT%calipso_cfad_sr(ij:ik,:,:), cospOUT%calipso_lidarcld(ij:ik,:),       &
               cospOUT%calipso_cldlayer(ij:ik,:),                                           &
               cospgridIN%at(:,:), calipso_betaperp_tot(:,:,:), cospgridIN%surfelev,        & 
               cospOUT%calipso_lidarcldphase(ij:ik,:,:),                       &
               cospOUT%calipso_lidarcldtype(ij:ik,:,:),  cospOUT%calipso_cldtype(ij:ik,:),  &
               cospOUT%calipso_cldtypetemp(ij:ik,:), cospOUT%calipso_cldtypemeanz(ij:ik,:), & 
               cospOUT%calipso_cldtypemeanzse(ij:ik,:), cospOUT%calipso_cldthinemis(ij:ik), &
               cospOUT%calipso_cldlayerphase(ij:ik,:,:), cospOUT%calipso_lidarcldtmp(ij:ik,:,:))
       endif 

       if (associated(cospOUT%calipso_srbval)) cospOUT%calipso_srbval = calipso_histBsct

       ! Free up memory (if necessary)
       if (allocated(out1D_1)) then
          deallocate(out1D_1)
          nullify(cospOUT%calipso_cfad_sr)
       endif
       if (allocated(out1D_2)) then
          deallocate(out1D_2)
          nullify(cospOUT%calipso_lidarcld)
       endif
       if (allocated(out1D_3)) then
          deallocate(out1D_3)
          nullify(cospOUT%calipso_lidarcldphase)
       endif
       if (allocated(out1D_4)) then
          deallocate(out1D_4)
          nullify(cospOUT%calipso_cldlayer)
       endif
       if (allocated(out1D_5)) then
          deallocate(out1D_5)
          nullify(cospOUT%calipso_cldlayerphase)
       endif
       if (allocated(out1D_6)) then
          deallocate(out1D_6)
          nullify(cospOUT%calipso_lidarcldtmp)
       endif
       if (allocated(out1D_7)) then
          deallocate(out1D_7)  
          nullify(cospOUT%calipso_lidarcldtype)
       endif  
       if (allocated(out1D_8)) then
          deallocate(out1D_8)
          nullify(cospOUT%calipso_cldtype)
       endif 
       if (allocated(out1D_9)) then 
          deallocate(out1D_9) 
          nullify(cospOUT%calipso_cldtypetemp)
       endif 
       if (allocated(out1D_10)) then 
          deallocate(out1D_10)
          nullify(cospOUT%calipso_cldtypemeanz)
       endif 
       if (allocated(out1D_12)) then 
          deallocate(out1D_12) 
          nullify(cospOUT%calipso_cldtypemeanzse)
       endif 
       if (allocated(out1D_11)) then 
          deallocate(out1D_11) 
          nullify(cospOUT%calipso_cldthinemis)
       endif

    endif

    ! GROUND LIDAR Simulator
    if (LgrLidar532_column) then
       ! Check to see which outputs are requested. If not requested, use a local dummy array
       if (.not. associated(cospOUT%grLidar532_cfad_sr)) then
          allocate(out1D_1(Npoints*SR_BINS*Nlvgrid)) 
          cospOUT%grLidar532_cfad_sr(ij:ik,1:SR_BINS,1:Nlvgrid) => out1D_1
       endif 
       if (.not. associated(cospOUT%grLidar532_lidarcld)) then
          allocate(out1D_2(Npoints*Nlvgrid)) 
          cospOUT%grLidar532_lidarcld(ij:ik,1:Nlvgrid) => out1D_2
       endif
       if (.not. associated(cospOUT%grLidar532_cldlayer)) then
          allocate(out1D_3(Npoints*LIDAR_NCAT))
          cospOUT%grLidar532_cldlayer(ij:ik,1:LIDAR_NCAT) => out1D_3
       endif
       
       ! Call simulator 
       ok_lidar_cfad_grLidar532=.true.
       call lidar_column(grLidar532IN%Npoints, grLidar532IN%Ncolumns, grLidar532IN%Nlevels,       &
            Nlvgrid, SR_BINS, LIDAR_NTYPE, 'grlidar532',grLidar532_beta_tot(:,:,:), grLidar532_beta_mol(:,:),&
            cospgridIN%phalf(:,2:grLidar532IN%Nlevels+1),cospgridIN%hgt_matrix,               &
            cospgridIN%hgt_matrix_half, vgrid_z(:), ok_lidar_cfad_grLidar532, LIDAR_NCAT,      &
            cospOUT%grLidar532_cfad_sr(ij:ik,:,:), cospOUT%grLidar532_lidarcld(ij:ik,:),       &
            cospOUT%grLidar532_cldlayer(ij:ik,:))

       if (associated(cospOUT%grLidar532_srbval)) cospOUT%grLidar532_srbval = grLidar532_histBsct

       ! Free up memory (if necessary) 
       if (allocated(out1D_1)) then
          deallocate(out1D_1)
          nullify(cospOUT%grLidar532_cfad_sr)
       endif 
       if (allocated(out1D_2)) then
          deallocate(out1D_2)
          nullify(cospOUT%grLidar532_lidarcld)
       endif 
       if (allocated(out1D_3)) then
          deallocate(out1D_3)
          nullify(cospOUT%grLidar532_cldlayer)
       endif

    endif

    ! ATLID Simulator
    if (Latlid_column) then
       ! Check to see which outputs are requested. If not requested, use a local dummy array
       if (.not. associated(cospOUT%atlid_cfad_sr)) then
          allocate(out1D_1(Npoints*SR_BINS*Nlvgrid))
          cospOUT%atlid_cfad_sr(ij:ik,1:SR_BINS,1:Nlvgrid) => out1D_1
       endif                                    
       if (.not. associated(cospOUT%atlid_lidarcld)) then 
          allocate(out1D_2(Npoints*Nlvgrid))                 
          cospOUT%atlid_lidarcld(ij:ik,1:Nlvgrid) => out1D_2 
       endif                                                   
       if (.not. associated(cospOUT%atlid_cldlayer)) then   
          allocate(out1D_3(Npoints*LIDAR_NCAT))                 
          cospOUT%atlid_cldlayer(ij:ik,1:LIDAR_NCAT) => out1D_3 
       endif                                                      
       ! Call simulator                                                             
       ok_lidar_cfad_atlid=.true.
       if (cospIN % cospswathsIN(4) % N_inst_swaths .gt. 0) then ! Trigger use of swathed arrays
          if (atlidIN%Npoints .gt. 0) then
             allocate(temp_atlid_cfad_sr(atlidIN%Npoints,SR_BINS,Nlvgrid), &
                      temp_atlid_lidarcld(atlidIN%Npoints,Nlvgrid),        &
                      temp_atlid_cldlayer(atlidIN%Npoints,LIDAR_NCAT))
             call lidar_column(atlidIN%Npoints, atlidIN%Ncolumns, atlidIN%Nlevels,                     &
                 Nlvgrid, SR_BINS, LIDAR_NTYPE, 'atlid',atlid_beta_tot(:,:,:),                         &
                 atlid_beta_mol(:,:), cospgridIN%phalf(int(ATLID_MASK_INDICES),2:atlidIN%Nlevels+1),   & 
                 cospgridIN%hgt_matrix(int(ATLID_MASK_INDICES),:),                                     &
                 cospgridIN%hgt_matrix_half(int(ATLID_MASK_INDICES),:), vgrid_z(:),                    &
                 ok_lidar_cfad_atlid, LIDAR_NCAT, temp_atlid_cfad_sr(:,:,:),                           &
                 temp_atlid_lidarcld(:,:), temp_atlid_cldlayer(:,:))
             ! Decode back to the cospOUT shapes
             cospOUT%atlid_cfad_sr(ij:ik,:,:) = R_UNDEF
             cospOUT%atlid_lidarcld(ij:ik,:)  = R_UNDEF
             cospOUT%atlid_cldlayer(ij:ik,:)  = R_UNDEF
             if (associated(cospOUT%atlid_cfad_sr))          cospOUT%atlid_cfad_sr(ij+int(ATLID_MASK_INDICES)-1,:,:)         = temp_atlid_cfad_sr(:,:,:)
             if (associated(cospOUT%atlid_lidarcld))         cospOUT%atlid_lidarcld(ij+int(ATLID_MASK_INDICES)-1,:)          = temp_atlid_lidarcld(:,:)
             if (associated(cospOUT%atlid_cldlayer))         cospOUT%atlid_cldlayer(ij+int(ATLID_MASK_INDICES)-1,:)          = temp_atlid_cldlayer(:,:)
             deallocate(temp_atlid_cfad_sr,temp_atlid_lidarcld,temp_atlid_cldlayer)
          else 
             cospOUT%atlid_cfad_sr(ij:ik,:,:) = R_UNDEF
             cospOUT%atlid_lidarcld(ij:ik,:)  = R_UNDEF
             cospOUT%atlid_cldlayer(ij:ik,:)  = R_UNDEF
          endif 
       else
          call lidar_column(atlidIN%Npoints, atlidIN%Ncolumns, atlidIN%Nlevels,       &
              Nlvgrid, SR_BINS, LIDAR_NTYPE, 'atlid',atlid_beta_tot(:,:,:),     &
              atlid_beta_mol(:,:), cospgridIN%phalf(:,2:atlidIN%Nlevels+1),            &
              cospgridIN%hgt_matrix, cospgridIN%hgt_matrix_half, vgrid_z(:),         &
              ok_lidar_cfad_atlid, LIDAR_NCAT, cospOUT%atlid_cfad_sr(ij:ik,:,:),     &
              cospOUT%atlid_lidarcld(ij:ik,:), cospOUT%atlid_cldlayer(ij:ik,:))
       endif
       if (associated(cospOUT%atlid_srbval)) cospOUT%atlid_srbval = atlid_histBsct 

       ! Free up memory (if necessary)
       if (allocated(ATLID_MASK_INDICES))               deallocate(ATLID_MASK_INDICES)
       if (allocated(out1D_1)) then
          deallocate(out1D_1)                 
          nullify(cospOUT%atlid_cfad_sr)
       endif                                  
       if (allocated(out1D_2)) then          
          deallocate(out1D_2)                
          nullify(cospOUT%atlid_lidarcld)
       endif                                  
       if (allocated(out1D_3)) then          
          deallocate(out1D_3)
          nullify(cospOUT%atlid_cldlayer)
       endif          

    endif

    ! PARASOL
    if (Lparasol_column) then
       if (cospIN % cospswathsIN(5) % N_inst_swaths .gt. 0) then ! Trigger use of swathed arrays
          if (parasolIN%Npoints .gt. 0) then 
             allocate(temp_parasolGrid_refl(parasolIN%Npoints,PARASOL_NREFL))
             call parasol_column(parasolIN%Npoints,PARASOL_NREFL,parasolIN%Ncolumns,                &
                                 cospgridIN%land(int(PARASOL_MASK_INDICES)),parasolPix_refl(:,:,:), &
                                 temp_parasolGrid_refl(:,:))
             ! Decode back to the cospOUT shapes
             cospOUT%parasolGrid_refl(ij:ik,:) = R_UNDEF
             ! Note: The active simulators use R_UNDEF to indicate unretrieved and clear-sky values.
             ! Decisions on how to treat these values are delegeate to the host user, and fields
             ! are often set to zero.
             ! Using R_UNDEF in this way is incompatible with the swath approach, where R_UNDEF means
             ! that a field should not be used in statistically averaging, etc.
             ! To avoid this error, when swathing we set R_UNDEF value to zero here.
             where (temp_parasolGrid_refl(:,:) == R_UNDEF) temp_parasolGrid_refl(:,:) = 0._wp
             if (associated(cospOUT%parasolGrid_refl))   cospOUT%parasolGrid_refl(ij+int(PARASOL_MASK_INDICES)-1,:)   = temp_parasolGrid_refl(:,:)
             deallocate(temp_parasolGrid_refl)
          else
             cospOUT%parasolGrid_refl(ij:ik,:) = R_UNDEF
          endif
       else
          call parasol_column(parasolIN%Npoints,PARASOL_NREFL,parasolIN%Ncolumns,           &
                               cospgridIN%land(:),parasolPix_refl(:,:,:),                   &
                               cospOUT%parasolGrid_refl(ij:ik,:))
         !  print*,'cospOUT%parasolGrid_refl:   ',cospOUT%parasolGrid_refl ! Values not zeroed here.
       endif
       if (allocated(parasolPix_refl))      deallocate(parasolPix_refl)
       if (allocated(PARASOL_MASK_INDICES)) deallocate(PARASOL_MASK_INDICES)        
    endif

    ! CLOUDSAT
    if (Lcloudsat_column) then
       ! Check to see which outputs are requested. If not requested, use a local dummy array
       if (.not. associated(cospOUT%cloudsat_cfad_ze)) then
          allocate(out1D_1(Npoints*cloudsat_DBZE_BINS*Nlvgrid))
          cospOUT%cloudsat_cfad_ze(ij:ik,1:cloudsat_DBZE_BINS,1:Nlvgrid) => out1D_1
       endif

       if (.not. associated(cospOUT%cloudsat_pia)) then
          allocate(out1D_2(Npoints))
          cospOUT%cloudsat_pia(ij:ik) => out1D_2
       endif
       if (.not. associated(cospOUT%cloudsat_precip_cover)) then
          allocate(out1D_3(Npoints*nCloudsatPrecipClass))
          cospOUT%cloudsat_precip_cover(ij:ik,1:nCloudsatPrecipClass) => out1D_3
       endif
          
       ! Call simulator
       if (cospIN % cospswathsIN(3) % N_inst_swaths .gt. 0) then ! Trigger use of swathed arrays
          if (cloudsatIN%Npoints .gt. 0) then
             allocate(temp_cloudsat_cfad_ze(cloudsatIN%Npoints,cloudsat_DBZE_BINS,Nlvgrid),    &
                      temp_cloudsat_precip_cover(cloudsatIN%Npoints,cloudsat_DBZE_BINS),       &
                      temp_cloudsat_pia(cloudsatIN%Npoints))
             call quickbeam_column(cloudsatIN%Npoints, cloudsatIN%Ncolumns, cloudsatIN%Nlevels,            &
                  Nlvgrid, cloudsat_DBZE_BINS, 'cloudsat', cloudsatDBZe, cloudsatZe_non,                   &
                  cospgridIN%land(int(CSCAL_MASK_INDICES)), cospgridIN%surfelev(int(CSCAL_MASK_INDICES)),  &
                  cospgridIN%at(int(CSCAL_MASK_INDICES),cospIN%Nlevels), cospIN%fracPrecipIce(int(CSCAL_MASK_INDICES),:),      &
                  cospgridIN%hgt_matrix(int(CSCAL_MASK_INDICES),:), cospgridIN%hgt_matrix_half(int(CSCAL_MASK_INDICES),:),     &
                  temp_cloudsat_cfad_ze(:,:,:), temp_cloudsat_precip_cover(:,:), temp_cloudsat_pia(:))
             ! Decode back to the cospOUT shapes
             cospOUT%cloudsat_cfad_ze(ij:ik,:,:)    = R_UNDEF
             cospOUT%cloudsat_precip_cover(ij:ik,:) = R_UNDEF
             cospOUT%cloudsat_pia(ij:ik)            = R_UNDEF
             ! Note: The active simulators use R_UNDEF to indicate unretrieved and clear-sky values.
             ! Decisions on how to treat these values are delegeate to the host user, and fields
             ! are often set to zero.
             ! Using R_UNDEF in this way is incompatible with the swath approach, where R_UNDEF means
             ! that a field should not be used in statistically averaging, etc.
             ! To avoid this error, when swathing we set R_UNDEF value to zero here.
             where (temp_cloudsat_cfad_ze(:,:,:) == R_UNDEF) temp_cloudsat_cfad_ze(:,:,:) = 0._wp
             if (associated(cospOUT%cloudsat_cfad_ze))        cospOUT%cloudsat_cfad_ze(ij+int(CSCAL_MASK_INDICES)-1,:,:)     = temp_cloudsat_cfad_ze(:,:,:)
             if (associated(cospOUT%cloudsat_precip_cover))   cospOUT%cloudsat_precip_cover(ij+int(CSCAL_MASK_INDICES)-1,:)  = temp_cloudsat_precip_cover(:,:)
             if (associated(cospOUT%cloudsat_pia))            cospOUT%cloudsat_pia(ij+int(CSCAL_MASK_INDICES)-1)             = temp_cloudsat_pia(:)
             deallocate(temp_cloudsat_cfad_ze,temp_cloudsat_precip_cover,temp_cloudsat_pia)
             if (allocated(CSCAL_SWATH_MASK)) deallocate(CSCAL_SWATH_MASK)
          else 
             cospOUT%cloudsat_cfad_ze(ij:ik,:,:)    = R_UNDEF
             cospOUT%cloudsat_precip_cover(ij:ik,:) = R_UNDEF
             cospOUT%cloudsat_pia(ij:ik)            = R_UNDEF
          endif
       else
          call quickbeam_column(cloudsatIN%Npoints, cloudsatIN%Ncolumns, cloudsatIN%Nlevels,&
               Nlvgrid, cloudsat_DBZE_BINS, 'cloudsat', cloudsatDBZe, cloudsatZe_non,       &
               cospgridIN%land(:), cospgridIN%surfelev(:), cospgridIN%at(:,cospIN%Nlevels), &
               cospIN%fracPrecipIce, cospgridIN%hgt_matrix, cospgridIN%hgt_matrix_half,     &
               cospOUT%cloudsat_cfad_ze(ij:ik,:,:), cospOUT%cloudsat_precip_cover(ij:ik,:), &
               cospOUT%cloudsat_pia(ij:ik))
       endif
       ! Free up memory  (if necessary)
       if (allocated(out1D_1)) then
          deallocate(out1D_1)
          nullify(cospOUT%cloudsat_cfad_ze)
       endif
       if (allocated(out1D_2)) then
          deallocate(out1D_2)
          nullify(cospOUT%cloudsat_pia)
       endif
       if (allocated(out1D_3)) then
          deallocate(out1D_3)
          nullify(cospOUT%cloudsat_precip_cover)
       endif
    endif

    ! MODIS
    if (Lmodis_column) then
       if (modisIN%nSunlit > 0) then
          ! Allocate space for local variables
          allocate(modisCftotal(modisIN%nSunlit), modisCfLiquid(modisIN%nSunlit),        &
                   modisCfIce(modisIN%nSunlit),modisCfHigh(modisIN%nSunlit),             &
                   modisCfMid(modisIN%nSunlit),modisCfLow(modisIN%nSunlit),              &
                   modisMeanTauTotal(modisIN%nSunlit),                                   &
                   modisMeanTauLiquid(modisIN%nSunlit),modisMeanTauIce(modisIN%nSunlit), &
                   modisMeanLogTauTotal(modisIN%nSunlit),                                &
                   modisMeanLogTauLiquid(modisIN%nSunlit),                               &
                   modisMeanLogTauIce(modisIN%nSunlit),                                  &
                   modisMeanSizeLiquid(modisIN%nSunlit),                                 &
                   modisMeanSizeIce(modisIN%nSunlit),                                    &
                   modisMeanCloudTopPressure(modisIN%nSunlit),                           &
                   modisMeanLiquidWaterPath(modisIN%nSunlit),                            &
                   modisMeanIceWaterPath(modisIN%nSunlit),                               &
                   modisJointHistogram(modisIN%nSunlit,numMODISTauBins,numMODISPresBins),&
                   modisJointHistogramIce(modisIN%nSunlit,numModisTauBins,numMODISReffIceBins),&
                   modisJointHistogramLiq(modisIN%nSunlit,numModisTauBins,numMODISReffLiqBins),&
                   modisJointHistogram_CtpCodLiq(modisIN%nSunlit,numMODISTauBins,numMODISPresBins),&
                   modisJointHistogram_CtpCodIce(modisIN%nSunlit,numMODISTauBins,numMODISPresBins),&
                   modisJointHistogram_LwpRefLiq(modisIN%nSunlit,numMODISLWPBins,numMODISReffLiqBins), &
                   modisJointHistogram_IwpRefIce(modisIN%nSunlit,numMODISIWPBins,numMODISReffIceBins) &
                   )
          ! Call simulator
          call modis_column(modisIN%nSunlit, modisIN%Ncolumns,modisRetrievedPhase,       &
                             modisRetrievedCloudTopPressure,modisRetrievedTau,           &
                             modisRetrievedSize, modisCfTotal, modisCfLiquid, modisCfIce,&
                             modisCfHigh, modisCfMid, modisCfLow, modisMeanTauTotal,     &
                             modisMeanTauLiquid, modisMeanTauIce, modisMeanLogTauTotal,  &
                             modisMeanLogTauLiquid, modisMeanLogTauIce,                  &
                             modisMeanSizeLiquid, modisMeanSizeIce,                      &
                             modisMeanCloudTopPressure, modisMeanLiquidWaterPath,        &
                             modisMeanIceWaterPath, modisJointHistogram,                 &
                             modisJointHistogramIce,modisJointHistogramLiq,              &
                             modisJointHistogram_CtpCodLiq,                              &
                             modisJointHistogram_CtpCodIce,                              &
                             modisJointHistogram_LwpRefLiq,                              &
                             modisJointHistogram_IwpRefIce                               &
                             )
          ! Store data (if requested)
          if (associated(cospOUT%modis_Cloud_Fraction_Total_Mean)) then
             cospOUT%modis_Cloud_Fraction_Total_Mean(ij+int(modisIN%sunlit(:))-1)   =    &
                  modisCfTotal
          endif
          if (associated(cospOUT%modis_Cloud_Fraction_Water_Mean)) then
             cospOUT%modis_Cloud_Fraction_Water_Mean(ij+int(modisIN%sunlit(:))-1)   =    &
                  modisCfLiquid
          endif
          if (associated(cospOUT%modis_Cloud_Fraction_Ice_Mean)) then
             cospOUT%modis_Cloud_Fraction_Ice_Mean(ij+int(modisIN%sunlit(:))-1)     =    &
                  modisCfIce
          endif
          if (associated(cospOUT%modis_Cloud_Fraction_High_Mean)) then
             cospOUT%modis_Cloud_Fraction_High_Mean(ij+int(modisIN%sunlit(:))-1)    =    &
                  modisCfHigh
          endif
          if (associated(cospOUT%modis_Cloud_Fraction_Mid_Mean)) then
             cospOUT%modis_Cloud_Fraction_Mid_Mean(ij+int(modisIN%sunlit(:))-1)     =    &
                  modisCfMid
          endif
          if (associated(cospOUT%modis_Cloud_Fraction_Low_Mean)) then
             cospOUT%modis_Cloud_Fraction_Low_Mean(ij+int(modisIN%sunlit(:))-1)     =    &
                  modisCfLow
          endif
          if (associated(cospOUT%modis_Optical_Thickness_Total_Mean)) then
             cospOUT%modis_Optical_Thickness_Total_Mean(ij+int(modisIN%sunlit(:))-1) =   &
                  modisMeanTauTotal
          endif
          if (associated(cospOUT%modis_Optical_Thickness_Water_Mean)) then
             cospOUT%modis_Optical_Thickness_Water_Mean(ij+int(modisIN%sunlit(:))-1) =   &
                  modisMeanTauLiquid
          endif
          if (associated(cospOUT%modis_Optical_Thickness_Ice_Mean)) then
             cospOUT%modis_Optical_Thickness_Ice_Mean(ij+int(modisIN%sunlit(:))-1)  =    &
                  modisMeanTauIce
          endif
          if (associated(cospOUT%modis_Optical_Thickness_Total_LogMean)) then
             cospOUT%modis_Optical_Thickness_Total_LogMean(ij+int(modisIN%sunlit(:))-1)= &
                  modisMeanLogTauTotal
          endif
          if (associated(cospOUT%modis_Optical_Thickness_Water_LogMean)) then
             cospOUT%modis_Optical_Thickness_Water_LogMean(ij+int(modisIN%sunlit(:))-1) = &
                  modisMeanLogTauLiquid
          endif
          if (associated(cospOUT%modis_Optical_Thickness_Ice_LogMean)) then
             cospOUT%modis_Optical_Thickness_Ice_LogMean(ij+int(modisIN%sunlit(:))-1) =  &
                  modisMeanLogTauIce
          endif
          if (associated(cospOUT%modis_Cloud_Particle_Size_Water_Mean)) then
             cospOUT%modis_Cloud_Particle_Size_Water_Mean(ij+int(modisIN%sunlit(:))-1) = &
                  modisMeanSizeLiquid
          endif
          if (associated(cospOUT%modis_Cloud_Particle_Size_Ice_Mean)) then
             cospOUT%modis_Cloud_Particle_Size_Ice_Mean(ij+int(modisIN%sunlit(:))-1) =   &
                  modisMeanSizeIce
          endif
          if (associated(cospOUT%modis_Cloud_Top_Pressure_Total_Mean)) then
             cospOUT%modis_Cloud_Top_Pressure_Total_Mean(ij+int(modisIN%sunlit(:))-1) =  &
                  modisMeanCloudTopPressure
          endif
          if (associated(cospOUT%modis_Liquid_Water_Path_Mean)) then
             cospOUT%modis_Liquid_Water_Path_Mean(ij+int(modisIN%sunlit(:))-1)      =    &
                  modisMeanLiquidWaterPath
          endif
          if (associated(cospOUT%modis_Ice_Water_Path_Mean)) then
              cospOUT%modis_Ice_Water_Path_Mean(ij+int(modisIN%sunlit(:))-1)         =   &
                  modisMeanIceWaterPath
          endif
          if (associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure)) then
             cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure(ij+            &
                  int(modisIN%sunlit(:))-1, 1:numModisTauBins, :) = modisJointHistogram(:, :, :)
             ! Reorder pressure bins in joint histogram to go from surface to TOA
             cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure(ij:ik,:,:) = &
                  cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure(ij:ik,:,numMODISPresBins:1:-1)
          endif

          if (associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Liq)) then
             cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Liq(ij+            &
                  int(modisIN%sunlit(:))-1, 1:numModisTauBins, :) = modisJointHistogram_CtpCodLiq(:, :, :)
             ! Reorder pressure bins in joint histogram to go from surface to TOA
             cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Liq(ij:ik,:,:) = &
                  cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Liq(ij:ik,:,numMODISPresBins:1:-1)
          endif
          if (associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Ice)) then
             cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Ice(ij+            &
                  int(modisIN%sunlit(:))-1, 1:numModisTauBins, :) = modisJointHistogram_CtpCodIce(:, :, :)
             ! Reorder pressure bins in joint histogram to go from surface to TOA
             cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Ice(ij:ik,:,:) = &
                  cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Ice(ij:ik,:,numMODISPresBins:1:-1)
          endif
          if (associated(cospOUT%modis_LWP_vs_ReffLIQ)) then
             cospOUT%modis_LWP_vs_ReffLIQ(ij+int(modisIN%sunlit(:))-1, 1:numMODISLWPBins,:) = &
                modisJointHistogram_LwpRefLiq(:,:,:)
          endif
          if (associated(cospOUT%modis_IWP_vs_ReffICE)) then
             cospOUT%modis_IWP_vs_ReffICE(ij+int(modisIN%sunlit(:))-1, 1:numMODISIWPBins,:) = &
                modisJointHistogram_IwpRefIce(:,:,:)
          endif

          if (associated(cospOUT%modis_Optical_Thickness_vs_ReffIce)) then
             cospOUT%modis_Optical_Thickness_vs_ReffIce(ij:ik,1:numMODISTauBins,:) = 0.0
             cospOUT%modis_Optical_Thickness_vs_ReffIce(ij+int(modisIN%sunlit(:))-1, 1:numMODISTauBins,:) = &
                modisJointHistogramIce(:,:,:)
          endif
          if (associated(cospOUT%modis_Optical_Thickness_vs_ReffLiq)) then
             cospOUT%modis_Optical_Thickness_vs_ReffLiq(ij:ik,:,:) = 0.0
             cospOUT%modis_Optical_Thickness_vs_ReffLiq(ij+int(modisIN%sunlit(:))-1, 1:numMODISTauBins,:) = &
                modisJointHistogramLiq(:,:,:)
          endif

          if(modisIN%nSunlit < modisIN%Npoints) then
             ! Where it's night and we haven't done the retrievals the values are undefined
             if (associated(cospOUT%modis_Cloud_Fraction_Total_Mean))                    &
                cospOUT%modis_Cloud_Fraction_Total_Mean(ij+int(modisIN%notSunlit(:))-1) = R_UNDEF
             if (associated(cospOUT%modis_Cloud_Fraction_Water_Mean))                    &
                cospOUT%modis_Cloud_Fraction_Water_Mean(ij+int(modisIN%notSunlit(:))-1) = R_UNDEF
             if (associated(cospOUT%modis_Cloud_Fraction_Ice_Mean))                      &
                cospOUT%modis_Cloud_Fraction_Ice_Mean(ij+int(modisIN%notSunlit(:))-1) = R_UNDEF
             if (associated(cospOUT%modis_Cloud_Fraction_High_Mean))                     &
                cospOUT%modis_Cloud_Fraction_High_Mean(ij+int(modisIN%notSunlit(:))-1) = R_UNDEF
             if (associated(cospOUT%modis_Cloud_Fraction_Mid_Mean))                      &
                cospOUT%modis_Cloud_Fraction_Mid_Mean(ij+int(modisIN%notSunlit(:))-1) = R_UNDEF
             if (associated(cospOUT%modis_Cloud_Fraction_Low_Mean))                      &
                cospOUT%modis_Cloud_Fraction_Low_Mean(ij+int(modisIN%notSunlit(:))-1) = R_UNDEF
             if (associated(cospOUT%modis_Optical_Thickness_Total_Mean))                 &
                cospOUT%modis_Optical_Thickness_Total_Mean(ij+int(modisIN%notSunlit(:))-1) = R_UNDEF
             if (associated(cospOUT%modis_Optical_Thickness_Water_Mean))                 &
                cospOUT%modis_Optical_Thickness_Water_Mean(ij+int(modisIN%notSunlit(:))-1) = R_UNDEF
             if (associated(cospOUT%modis_Optical_Thickness_Ice_Mean))                   &
                cospOUT%modis_Optical_Thickness_Ice_Mean(ij+int(modisIN%notSunlit(:))-1) = R_UNDEF
             if (associated(cospOUT%modis_Optical_Thickness_Total_LogMean))              &
                cospOUT%modis_Optical_Thickness_Total_LogMean(ij+int(modisIN%notSunlit(:))-1) = R_UNDEF
             if (associated(cospOUT%modis_Optical_Thickness_Water_LogMean))              &
                cospOUT%modis_Optical_Thickness_Water_LogMean(ij+int(modisIN%notSunlit(:))-1) = R_UNDEF
             if (associated(cospOUT%modis_Optical_Thickness_Ice_LogMean))                &
                cospOUT%modis_Optical_Thickness_Ice_LogMean(ij+int(modisIN%notSunlit(:))-1) = R_UNDEF
             if (associated(cospOUT%modis_Cloud_Particle_Size_Water_Mean))               &
                cospOUT%modis_Cloud_Particle_Size_Water_Mean(ij+int(modisIN%notSunlit(:))-1) = R_UNDEF
             if (associated(cospOUT%modis_Cloud_Particle_Size_Ice_Mean))                 &
                cospOUT%modis_Cloud_Particle_Size_Ice_Mean(ij+int(modisIN%notSunlit(:))-1) = R_UNDEF
             if (associated(cospOUT%modis_Cloud_Top_Pressure_Total_Mean))                &
                cospOUT%modis_Cloud_Top_Pressure_Total_Mean(ij+int(modisIN%notSunlit(:))-1) = R_UNDEF
             if (associated(cospOUT%modis_Liquid_Water_Path_Mean))                       &
                cospOUT%modis_Liquid_Water_Path_Mean(ij+int(modisIN%notSunlit(:))-1) = R_UNDEF
             if (associated(cospOUT%modis_Ice_Water_Path_Mean))                          &
                cospOUT%modis_Ice_Water_Path_Mean(ij+int(modisIN%notSunlit(:))-1) = R_UNDEF
             if (associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure))      &
                cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure(ij+int(modisIN%notSunlit(:))-1, :, :) = R_UNDEF
             if (associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Liq))      &
                cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Liq(ij+int(modisIN%notSunlit(:))-1, :, :) = R_UNDEF
             if (associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Ice))      &
                cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Ice(ij+int(modisIN%notSunlit(:))-1, :, :) = R_UNDEF
          end if
       else
          ! It's nightime everywhere - everything is undefined
          if (associated(cospOUT%modis_Cloud_Fraction_Total_Mean))                       &
             cospOUT%modis_Cloud_Fraction_Total_Mean(ij:ik) = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_Water_Mean))                       &
             cospOUT%modis_Cloud_Fraction_Water_Mean(ij:ik) = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_Ice_Mean))                         &
             cospOUT%modis_Cloud_Fraction_Ice_Mean(ij:ik) = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_High_Mean))                        &
             cospOUT%modis_Cloud_Fraction_High_Mean(ij:ik) = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_Mid_Mean))                         &
             cospOUT%modis_Cloud_Fraction_Mid_Mean(ij:ik) = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_Low_Mean))                         &
             cospOUT%modis_Cloud_Fraction_Low_Mean(ij:ik) = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Total_Mean))                    &
             cospOUT%modis_Optical_Thickness_Total_Mean(ij:ik) = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Water_Mean))                    &
             cospOUT%modis_Optical_Thickness_Water_Mean(ij:ik) = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Ice_Mean))                      &
             cospOUT%modis_Optical_Thickness_Ice_Mean(ij:ik) = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Total_LogMean))                 &
             cospOUT%modis_Optical_Thickness_Total_LogMean(ij:ik) = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Water_LogMean))                 &
             cospOUT%modis_Optical_Thickness_Water_LogMean(ij:ik) = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Ice_LogMean))                   &
             cospOUT%modis_Optical_Thickness_Ice_LogMean(ij:ik) = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Particle_Size_Water_Mean))                  &
             cospOUT%modis_Cloud_Particle_Size_Water_Mean(ij:ik) = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Particle_Size_Ice_Mean))                    &
              cospOUT%modis_Cloud_Particle_Size_Ice_Mean(ij:ik) = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Top_Pressure_Total_Mean))                   &
             cospOUT%modis_Cloud_Top_Pressure_Total_Mean(ij:ik) = R_UNDEF
          if (associated(cospOUT%modis_Liquid_Water_Path_Mean))                          &
             cospOUT%modis_Liquid_Water_Path_Mean(ij:ik) = R_UNDEF
          if (associated(cospOUT%modis_Ice_Water_Path_Mean))                             &
             cospOUT%modis_Ice_Water_Path_Mean(ij:ik) = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure))         &
             cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure(ij:ik, :, :) = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Liq))      &
             cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Liq(ij+int(modisIN%notSunlit(:))-1, :, :) = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Ice))      &
             cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Ice(ij+int(modisIN%notSunlit(:))-1, :, :) = R_UNDEF
       endif
       ! Free up memory (if necessary)
       if (allocated(modisRetrievedTau))               deallocate(modisRetrievedTau)
       if (allocated(modisRetrievedSize))              deallocate(modisRetrievedSize)
       if (allocated(modisRetrievedPhase))             deallocate(modisRetrievedPhase)
       if (allocated(modisRetrievedCloudTopPressure))  deallocate(modisRetrievedCloudTopPressure)
       if (allocated(modisCftotal))                    deallocate(modisCftotal)
       if (allocated(modisCfLiquid))                   deallocate(modisCfLiquid)
       if (allocated(modisCfIce))                      deallocate(modisCfIce)
       if (allocated(modisCfHigh))                     deallocate(modisCfHigh)
       if (allocated(modisCfMid))                      deallocate(modisCfMid)
       if (allocated(modisCfLow))                      deallocate(modisCfLow)
       if (allocated(modisMeanTauTotal))               deallocate(modisMeanTauTotal)
       if (allocated(modisMeanTauLiquid))              deallocate(modisMeanTauLiquid)
       if (allocated(modisMeanTauIce))                 deallocate(modisMeanTauIce)
       if (allocated(modisMeanLogTauTotal))            deallocate(modisMeanLogTauTotal)
       if (allocated(modisMeanLogTauLiquid))           deallocate(modisMeanLogTauLiquid)
       if (allocated(modisMeanLogTauIce))              deallocate(modisMeanLogTauIce)
       if (allocated(modisMeanSizeLiquid))             deallocate(modisMeanSizeLiquid)
       if (allocated(modisMeanSizeIce))                deallocate(modisMeanSizeIce)
       if (allocated(modisMeanCloudTopPressure))       deallocate(modisMeanCloudTopPressure)
       if (allocated(modisMeanLiquidWaterPath))        deallocate(modisMeanLiquidWaterPath)
       if (allocated(modisMeanIceWaterPath))           deallocate(modisMeanIceWaterPath)
       if (allocated(modisJointHistogram))             deallocate(modisJointHistogram)
       if (allocated(modisJointHistogram_CtpCodLiq))   deallocate(modisJointHistogram_CtpCodLiq)
       if (allocated(modisJointHistogram_CtpCodIce))   deallocate(modisJointHistogram_CtpCodIce)
       if (allocated(modisJointHistogram_LwpRefLiq))   deallocate(modisJointHistogram_LwpRefLiq)
       if (allocated(modisJointHistogram_IwpRefIce))   deallocate(modisJointHistogram_IwpRefIce)
       if (allocated(modisJointHistogramIce))          deallocate(modisJointHistogramIce)
       if (allocated(modisJointHistogramLiq))          deallocate(modisJointHistogramLiq)
       if (allocated(isccp_boxttop))                   deallocate(isccp_boxttop)
       if (allocated(isccp_boxptop))                   deallocate(isccp_boxptop)
       if (allocated(isccp_boxtau))                    deallocate(isccp_boxtau)
       if (allocated(isccp_meantbclr))                 deallocate(isccp_meantbclr)
       if (allocated(isccpLEVMATCH))                   deallocate(isccpLEVMATCH)
       if (allocated(MODIS_SWATH_MASK))                deallocate(MODIS_SWATH_MASK)
    endif

    ! RTTOV multi-instrument
    if (Lrttov_column) then
        do i=1,cospIN%Ninst_rttov
           ! Allocate memory for the outputs - I won't need all of these in every situation.
           ! Only allocate clear-sky memory when PC-RTTOV is run.
           if (cospIN % cfg_rttov(i) % Lrttov_pc) then 
               allocate(rttov_bt_clear(rttovIN%Npoints,cospIN  % cfg_rttov(i) % nchan_out)) ! all-sky brightness temp
               allocate(rttov_rad_clear(rttovIN%Npoints,cospIN % cfg_rttov(i) % nchan_out)) ! all-sky radiance
               rttov_bt_clear(:,:)  = R_UNDEF
               rttov_rad_clear(:,:) = R_UNDEF
               ! Run simulator
               call cosp_rttov_simulate(rttovIN,cospIN%cfg_rttov(i),                   & ! Inputs
                                        cosp_simulator(nError+1),                      & ! Error message holder
                                        bt_clear=rttov_bt_clear,                       & ! Clear-sky BT
                                        rad_clear=rttov_rad_clear)                       ! Clear-sky radiance
           else 
               allocate(rttov_bt_total(rttovIN%Npoints,cospIN   % cfg_rttov(i) % nchan_out))   ! all-sky brightness temp
               allocate(rttov_bt_clear(rttovIN%Npoints,cospIN   % cfg_rttov(i) % nchan_out))   ! clear-sky brightness temp
               allocate(rttov_rad_total(rttovIN%Npoints,cospIN  % cfg_rttov(i) % nchan_out))  ! all-sky brightness temp
               allocate(rttov_rad_clear(rttovIN%Npoints,cospIN  % cfg_rttov(i) % nchan_out))  ! clear-sky brightness temp
               allocate(rttov_rad_cloudy(rttovIN%Npoints,cospIN % cfg_rttov(i) % nchan_out)) ! cloudy-sky brightness temp
               allocate(rttov_refl_total(rttovIN%Npoints,cospIN % cfg_rttov(i) % nchan_out)) ! all-sky Bi-directional reflectance factor
               allocate(rttov_refl_clear(rttovIN%Npoints,cospIN % cfg_rttov(i) % nchan_out)) ! clear-sky Bi-directional reflectance factor
               ! Init to R_UNDEF
               rttov_bt_total(:,:)   = R_UNDEF
               rttov_bt_clear(:,:)   = R_UNDEF               
               rttov_rad_total(:,:)  = R_UNDEF               
               rttov_rad_clear(:,:)  = R_UNDEF               
               rttov_rad_cloudy(:,:) = R_UNDEF               
               rttov_refl_total(:,:) = R_UNDEF               
               rttov_refl_clear(:,:) = R_UNDEF               
               ! Run simulator
               call cosp_rttov_simulate(rttovIN,cospIN%cfg_rttov(i),                   & ! Inputs
                                        cosp_simulator(nError+1),                      & ! Error message holder
                                        bt_total=rttov_bt_total,                       & ! Brightness Temp Outputs
                                        bt_clear=rttov_bt_clear,                       &
                                        rad_total=rttov_rad_total,                     & ! Radiance Outputs
                                        rad_clear=rttov_rad_clear,                     &
                                        rad_cloudy=rttov_rad_cloudy,                   & 
                                        refl_total=rttov_refl_total,                   & ! Reflectance Outputs
                                        refl_clear=rttov_refl_clear,                   &
                                        debug=verbose)
           endif

           ! Write to cospOUT
           if (associated(cospOUT % rttov_outputs(i) % channel_indices))                    &
               cospOUT % rttov_outputs(i) % channel_indices(:) = cospIN % cfg_rttov(i) % iChannel         
           if (cospIN % cfg_rttov(i) % Lrttov_pc) then
               if (associated(cospOUT % rttov_outputs(i) % bt_total_pc))                    &
                  cospOUT % rttov_outputs(i) % bt_total_pc(ij:ik,:) = rttov_bt_clear
               if (associated(cospOUT % rttov_outputs(i) % rad_total_pc))                   &
                  cospOUT % rttov_outputs(i) % rad_total_pc(ij:ik,:) = rttov_rad_clear
           else
               if (associated(cospOUT % rttov_outputs(i) % bt_total))                    &
                  cospOUT % rttov_outputs(i) % bt_total(ij:ik,:) = rttov_bt_total
               if (associated(cospOUT % rttov_outputs(i) % bt_clear))                    &
                  cospOUT % rttov_outputs(i) % bt_clear(ij:ik,:) = rttov_bt_clear
               if (associated(cospOUT % rttov_outputs(i) % rad_total))                   &
                  cospOUT % rttov_outputs(i) % rad_total(ij:ik,:) = rttov_rad_total 
               if (associated(cospOUT % rttov_outputs(i) % rad_clear))                   &
                  cospOUT % rttov_outputs(i) % rad_clear(ij:ik,:) = rttov_rad_clear
               if (associated(cospOUT % rttov_outputs(i) % rad_cloudy))                  &
                  cospOUT % rttov_outputs(i) % rad_cloudy(ij:ik,:) = rttov_rad_cloudy
               if (associated(cospOUT % rttov_outputs(i) % refl_total))                  &
                  cospOUT % rttov_outputs(i) % refl_total(ij:ik,:) = rttov_refl_total
               if (associated(cospOUT % rttov_outputs(i) % refl_clear))                  &
                  cospOUT % rttov_outputs(i) % refl_clear(ij:ik,:) = rttov_refl_clear
           endif 

           ! Free up memory from output (if necessary)
           if (allocated(rttov_bt_total))               deallocate(rttov_bt_total)          
           if (allocated(rttov_bt_clear))               deallocate(rttov_bt_clear)          
           if (allocated(rttov_rad_total))              deallocate(rttov_rad_total)          
           if (allocated(rttov_rad_clear))              deallocate(rttov_rad_clear)          
           if (allocated(rttov_rad_cloudy))             deallocate(rttov_rad_cloudy)          
           if (allocated(rttov_refl_total))             deallocate(rttov_refl_total)          
           if (allocated(rttov_refl_clear))             deallocate(rttov_refl_clear)   

        end do
    
    endif    

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 6) Compute multi-instrument products
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    ! CLOUDSAT/CALIPSO products
    if (Lradar_lidar_tcc .or. Llidar_only_freq_cloud .or. Lcloudsat_tcc .or. Lcloudsat_tcc2) then
       if (calipsoIN%Npoints .gt. 0) then
          if (use_vgrid) then
             allocate(lidar_only_freq_cloud(cloudsatIN%Npoints,Nlvgrid),                    &
                  radar_lidar_tcc(cloudsatIN%Npoints), cloudsat_tcc(cloudsatIN%Npoints),    &
                  cloudsat_tcc2(cloudsatIN%Npoints))
             allocate(betamol_in(cloudsatIN%Npoints,1,cloudsatIN%Nlevels),                  &
                      betamolI(cloudsatIN%Npoints,1,Nlvgrid),                               &
                      pnormI(cloudsatIN%Npoints,cloudsatIN%Ncolumns,Nlvgrid),               &
                      Ze_totI(cloudsatIN%Npoints,cloudsatIN%Ncolumns,Nlvgrid))

             ! Regrid in the vertical (*NOTE* This routine requires SFC-2-TOA ordering, so flip
             ! inputs and outputs to maintain TOA-2-SFC ordering convention in COSP2.)
             ! Use CLOUDSAT masking array here (it is the same as calipso)
             if (cospIN % cospswathsIN(3) % N_inst_swaths .gt. 0) then ! Trigger use of swathed arrays
                betamol_in(:,1,:) = calipso_beta_mol(:,cloudsatIN%Nlevels:1:-1)
                call cosp_change_vertical_grid(cloudsatIN%Npoints,1,cloudsatIN%Nlevels,        &
                     cospgridIN%hgt_matrix(int(CSCAL_MASK_INDICES(:)),cloudsatIN%Nlevels:1:-1),                         &
                     cospgridIN%hgt_matrix_half(int(CSCAL_MASK_INDICES(:)),cloudsatIN%Nlevels:1:-1),betamol_in,         &
                     Nlvgrid,vgrid_zl(Nlvgrid:1:-1),vgrid_zu(Nlvgrid:1:-1),                    &
                     betamolI(:,1,Nlvgrid:1:-1))
                call cosp_change_vertical_grid(cloudsatIN%Npoints,cloudsatIN%Ncolumns,         &
                     cloudsatIN%Nlevels,cospgridIN%hgt_matrix(int(CSCAL_MASK_INDICES(:)),cloudsatIN%Nlevels:1:-1),      &
                     cospgridIN%hgt_matrix_half(int(CSCAL_MASK_INDICES(:)),cloudsatIN%Nlevels:1:-1),                    &
                     calipso_beta_tot(:,:,cloudsatIN%Nlevels:1:-1),Nlvgrid,                    &
                     vgrid_zl(Nlvgrid:1:-1),vgrid_zu(Nlvgrid:1:-1),pnormI(:,:,Nlvgrid:1:-1))
                call cosp_change_vertical_grid(cloudsatIN%Npoints,cloudsatIN%Ncolumns,         &
                     cloudsatIN%Nlevels,cospgridIN%hgt_matrix(int(CSCAL_MASK_INDICES(:)),cloudsatIN%Nlevels:1:-1),      &
                     cospgridIN%hgt_matrix_half(int(CSCAL_MASK_INDICES(:)),cloudsatIN%Nlevels:1:-1),                    &
                     cloudsatDBZe(:,:,cloudsatIN%Nlevels:1:-1),Nlvgrid,vgrid_zl(Nlvgrid:1:-1), &
                     vgrid_zu(Nlvgrid:1:-1),Ze_totI(:,:,Nlvgrid:1:-1),log_units=.true.)
             else
                betamol_in(:,1,:) = calipso_beta_mol(:,cloudsatIN%Nlevels:1:-1)
                call cosp_change_vertical_grid(cloudsatIN%Npoints,1,cloudsatIN%Nlevels,        &
                     cospgridIN%hgt_matrix(:,cloudsatIN%Nlevels:1:-1),                         &
                     cospgridIN%hgt_matrix_half(:,cloudsatIN%Nlevels:1:-1),betamol_in,         &
                     Nlvgrid,vgrid_zl(Nlvgrid:1:-1),vgrid_zu(Nlvgrid:1:-1),                    &
                     betamolI(:,1,Nlvgrid:1:-1))

                call cosp_change_vertical_grid(cloudsatIN%Npoints,cloudsatIN%Ncolumns,         &
                     cloudsatIN%Nlevels,cospgridIN%hgt_matrix(:,cloudsatIN%Nlevels:1:-1),      &
                     cospgridIN%hgt_matrix_half(:,cloudsatIN%Nlevels:1:-1),                    &
                     calipso_beta_tot(:,:,cloudsatIN%Nlevels:1:-1),Nlvgrid,                    &
                     vgrid_zl(Nlvgrid:1:-1),vgrid_zu(Nlvgrid:1:-1),pnormI(:,:,Nlvgrid:1:-1))

                call cosp_change_vertical_grid(cloudsatIN%Npoints,cloudsatIN%Ncolumns,         &
                     cloudsatIN%Nlevels,cospgridIN%hgt_matrix(:,cloudsatIN%Nlevels:1:-1),      &
                     cospgridIN%hgt_matrix_half(:,cloudsatIN%Nlevels:1:-1),                    &
                     cloudsatDBZe(:,:,cloudsatIN%Nlevels:1:-1),Nlvgrid,vgrid_zl(Nlvgrid:1:-1), &
                     vgrid_zu(Nlvgrid:1:-1),Ze_totI(:,:,Nlvgrid:1:-1),log_units=.true.)
             end if
             call cosp_lidar_only_cloud(cloudsatIN%Npoints, cloudsatIN%Ncolumns, Nlvgrid,   &
                  pnormI, betamolI, Ze_totI, lidar_only_freq_cloud, radar_lidar_tcc,          &
                  cloudsat_tcc, cloudsat_tcc2)        
             deallocate(betamol_in,betamolI,pnormI,ze_totI)
          else
             allocate(lidar_only_freq_cloud(cloudsatIN%Npoints,cloudsatIN%Nlevels),         &
                  radar_lidar_tcc(cloudsatIN%Npoints), cloudsat_tcc(cloudsatIN%Npoints),    &
                  cloudsat_tcc2(cloudsatIN%Npoints))
             call cosp_lidar_only_cloud(cloudsatIN%Npoints,cloudsatIN%Ncolumns,             &
                  cospIN%Nlevels,calipso_beta_tot(:,:,cloudsatIN%Nlevels:1:-1),             &
                  calipso_beta_mol(:,cloudsatIN%Nlevels:1:-1),                              &
                  cloudsatDBZe(:,:,cloudsatIN%Nlevels:1:-1),lidar_only_freq_cloud,          &
                  radar_lidar_tcc, cloudsat_tcc, cloudsat_tcc2)
          endif
       endif

       ! Store, when necessary
       if (cospIN % cospswathsIN(3) % N_inst_swaths .gt. 0) then ! Trigger use of swathed arrays
          if (associated(cospOUT%lidar_only_freq_cloud))   then 
             cospOUT%lidar_only_freq_cloud(ij:ik,:)  = R_UNDEF
             cospOUT%lidar_only_freq_cloud(ij+int(CSCAL_MASK_INDICES(:))-1,:)  = lidar_only_freq_cloud(:,:)
          endif
          if (associated(cospOUT%radar_lidar_tcc))         then
             cospOUT%radar_lidar_tcc(ij:ik)  = R_UNDEF
             cospOUT%radar_lidar_tcc(ij+int(CSCAL_MASK_INDICES(:))-1)          = radar_lidar_tcc(:)
          endif
          if (associated(cospOUT%cloudsat_tcc))            then
             cospOUT%cloudsat_tcc(ij:ik)  = R_UNDEF
             cospOUT%cloudsat_tcc(ij+int(CSCAL_MASK_INDICES(:))-1)             = cloudsat_tcc(:)
          endif
          if (associated(cospOUT%cloudsat_tcc2))           then
             cospOUT%cloudsat_tcc2(ij:ik)  = R_UNDEF
             cospOUT%cloudsat_tcc2(ij+int(CSCAL_MASK_INDICES(:))-1)            = cloudsat_tcc2(:)
          endif
       else
          if (associated(cospOUT%lidar_only_freq_cloud)) then
             cospOUT%lidar_only_freq_cloud(ij:ik,:) = lidar_only_freq_cloud
          endif
          if (associated(cospOUT%radar_lidar_tcc)) then
             cospOUT%radar_lidar_tcc(ij:ik) = radar_lidar_tcc
          endif
          if (associated(cospOUT%cloudsat_tcc)) then
             cospOUT%cloudsat_tcc(ij:ik) = cloudsat_tcc
          endif
          if (associated(cospOUT%cloudsat_tcc2)) then
             cospOUT%cloudsat_tcc2(ij:ik) = cloudsat_tcc2
          endif
       endif
    endif

    ! CloudSat/MODIS joint products (CFODDs and Occurrence Frequency of Warm Clouds)
    if (Lcloudsat_modis_wr) then
       if (cloudsatIN%Npoints .gt. 0) then
          allocate( cfodd_ntotal(cloudsatIN%Npoints, CFODD_NDBZE, CFODD_NICOD, CFODD_NCLASS) )
          allocate( wr_occfreq_ntotal(cloudsatIN%Npoints, WR_NREGIME) )

          if ( use_vgrid ) then
             !! interporation for fixed vertical grid:
             allocate( zlev(cloudsatIN%Npoints,Nlvgrid),                         &
                       t_in(cloudsatIN%Npoints,1,cloudsatIN%Nlevels),            &
                       tempI(cloudsatIN%Npoints,1,Nlvgrid),                      &
                       Ze_totI(cloudsatIN%Npoints,cloudsatIN%Ncolumns,Nlvgrid),  &
                       frac_outI(cloudsatIN%Npoints,cloudsatIN%Ncolumns,Nlvgrid) )
             do k = 1, Nlvgrid
                zlev(:,k) = vgrid_zu(k)
             enddo
             ! Use CLOUDSAT masking array here (it is the same as calipso)
             if (cospIN % cospswathsIN(3) % N_inst_swaths .gt. 0) then ! Trigger use of swathed arrays
                t_in(:,1,:) = cospgridIN%at(int(CSCAL_MASK_INDICES(:)),:)
                call cosp_change_vertical_grid (                                    &
                     cloudsatIN%Npoints, 1, cloudsatIN%Nlevels,                     &
                     cospgridIN%hgt_matrix(int(CSCAL_MASK_INDICES(:)),cloudsatIN%Nlevels:1:-1),              &
                     cospgridIN%hgt_matrix_half(int(CSCAL_MASK_INDICES(:)),cloudsatIN%Nlevels:1:-1),         &
                     t_in(:,:,cloudsatIN%Nlevels:1:-1), Nlvgrid,                    &
                     vgrid_zl(Nlvgrid:1:-1), vgrid_zu(Nlvgrid:1:-1),                &
                     tempI(:,:,Nlvgrid:1:-1)                                        )
                call cosp_change_vertical_grid (                                    &
                     cloudsatIN%Npoints, cloudsatIN%Ncolumns, cloudsatIN%Nlevels,   &
                     cospgridIN%hgt_matrix(int(CSCAL_MASK_INDICES(:)),cloudsatIN%Nlevels:1:-1),              &
                     cospgridIN%hgt_matrix_half(int(CSCAL_MASK_INDICES(:)),cloudsatIN%Nlevels:1:-1),         &
                     cloudsatDBZe(:,:,cloudsatIN%Nlevels:1:-1), Nlvgrid,            &
                     vgrid_zl(Nlvgrid:1:-1), vgrid_zu(Nlvgrid:1:-1),                &
                     Ze_totI(:,:,Nlvgrid:1:-1), log_units=.true.                    )
                call cosp_change_vertical_grid (                                    &
                     cloudsatIN%Npoints, cloudsatIN%Ncolumns, cloudsatIN%Nlevels,   &
                     cospgridIN%hgt_matrix(int(CSCAL_MASK_INDICES(:)),cloudsatIN%Nlevels:1:-1),              &
                     cospgridIN%hgt_matrix_half(int(CSCAL_MASK_INDICES(:)),cloudsatIN%Nlevels:1:-1),         &
                     cospIN%frac_out(int(CSCAL_MASK_INDICES(:)),:,cloudsatIN%Nlevels:1:-1), Nlvgrid,         &
                     vgrid_zl(Nlvgrid:1:-1), vgrid_zu(Nlvgrid:1:-1),                &
                     frac_outI(:,:,Nlvgrid:1:-1)                                    )
                call cosp_diag_warmrain(                                            &
                     cloudsatIN%Npoints, cloudsatIN%Ncolumns, Nlvgrid,              & !! in
                     tempI, zlev,                                                   & !! in
                     cospOUT%modis_Liquid_Water_Path_Mean(ij+int(CSCAL_MASK_INDICES(:))-1),                          & !! in
                     cospOUT%modis_Optical_Thickness_Water_Mean(ij+int(CSCAL_MASK_INDICES(:))-1),                    & !! in
                     cospOUT%modis_Cloud_Particle_Size_Water_Mean(ij+int(CSCAL_MASK_INDICES(:))-1),                  & !! in
                     cospOUT%modis_Cloud_Fraction_Water_Mean(ij+int(CSCAL_MASK_INDICES(:))-1),                       & !! in
                     cospOUT%modis_Ice_Water_Path_Mean(ij+int(CSCAL_MASK_INDICES(:))-1),                             & !! in
                     cospOUT%modis_Optical_Thickness_Ice_Mean(ij+int(CSCAL_MASK_INDICES(:))-1),                      & !! in
                     cospOUT%modis_Cloud_Particle_Size_Ice_Mean(ij+int(CSCAL_MASK_INDICES(:))-1),                    & !! in
                     cospOUT%modis_Cloud_Fraction_Ice_Mean(ij+int(CSCAL_MASK_INDICES(:))-1),                         & !! in
                     frac_outI,                                                     & !! in
                     Ze_totI,                                                       & !! in
                     cfodd_ntotal, wr_occfreq_ntotal                                ) !! inout           
             else
                t_in(:,1,:) = cospgridIN%at(:,:)
                call cosp_change_vertical_grid (                                    &
                     cloudsatIN%Npoints, 1, cloudsatIN%Nlevels,                     &
                     cospgridIN%hgt_matrix(:,cloudsatIN%Nlevels:1:-1),              &
                     cospgridIN%hgt_matrix_half(:,cloudsatIN%Nlevels:1:-1),         &
                     t_in(:,:,cloudsatIN%Nlevels:1:-1), Nlvgrid,                    &
                     vgrid_zl(Nlvgrid:1:-1), vgrid_zu(Nlvgrid:1:-1),                &
                     tempI(:,:,Nlvgrid:1:-1)                                        )
                call cosp_change_vertical_grid (                                    &
                     cloudsatIN%Npoints, cloudsatIN%Ncolumns, cloudsatIN%Nlevels,   &
                     cospgridIN%hgt_matrix(:,cloudsatIN%Nlevels:1:-1),              &
                     cospgridIN%hgt_matrix_half(:,cloudsatIN%Nlevels:1:-1),         &
                     cloudsatDBZe(:,:,cloudsatIN%Nlevels:1:-1), Nlvgrid,            &
                     vgrid_zl(Nlvgrid:1:-1), vgrid_zu(Nlvgrid:1:-1),                &
                     Ze_totI(:,:,Nlvgrid:1:-1), log_units=.true.                    )
                call cosp_change_vertical_grid (                                    &
                     cloudsatIN%Npoints, cloudsatIN%Ncolumns, cloudsatIN%Nlevels,   &
                     cospgridIN%hgt_matrix(:,cloudsatIN%Nlevels:1:-1),              &
                     cospgridIN%hgt_matrix_half(:,cloudsatIN%Nlevels:1:-1),         &
                     cospIN%frac_out(:,:,cloudsatIN%Nlevels:1:-1), Nlvgrid,         &
                     vgrid_zl(Nlvgrid:1:-1), vgrid_zu(Nlvgrid:1:-1),                &
                     frac_outI(:,:,Nlvgrid:1:-1)                                    )
                call cosp_diag_warmrain(                                            &
                     cloudsatIN%Npoints, cloudsatIN%Ncolumns, Nlvgrid,              & !! in
                     tempI, zlev,                                                   & !! in
                     cospOUT%modis_Liquid_Water_Path_Mean(ij:ik),                   & !! in
                     cospOUT%modis_Optical_Thickness_Water_Mean(ij:ik),             & !! in
                     cospOUT%modis_Cloud_Particle_Size_Water_Mean(ij:ik),           & !! in
                     cospOUT%modis_Cloud_Fraction_Water_Mean(ij:ik),                & !! in
                     cospOUT%modis_Ice_Water_Path_Mean(ij:ik),                      & !! in
                     cospOUT%modis_Optical_Thickness_Ice_Mean(ij:ik),               & !! in
                     cospOUT%modis_Cloud_Particle_Size_Ice_Mean(ij:ik),             & !! in
                     cospOUT%modis_Cloud_Fraction_Ice_Mean(ij:ik),                  & !! in
                     frac_outI,                                                     & !! in
                     Ze_totI,                                                       & !! in
                     cfodd_ntotal, wr_occfreq_ntotal                                ) !! inout
             endif
             deallocate( zlev, t_in, tempI, frac_outI, Ze_totI )
          else  ! do not use vgrid interporation ---------------------------------------!
             !! original model grid
             if (cospIN % cospswathsIN(3) % N_inst_swaths .gt. 0) then ! Trigger use of swathed arrays
                call cosp_diag_warmrain(                                                             &
                      cloudsatIN%Npoints, cloudsatIN%Ncolumns, cospIN%Nlevels,                       & !! in
                      cospgridIN%at(int(CSCAL_MASK_INDICES(:)),:),                                   & !! in
                      cospgridIN%hgt_matrix(int(CSCAL_MASK_INDICES(:)),:),                           & !! in
                      cospOUT%modis_Liquid_Water_Path_Mean(ij+int(CSCAL_MASK_INDICES(:))-1),         & !! in
                      cospOUT%modis_Optical_Thickness_Water_Mean(ij+int(CSCAL_MASK_INDICES(:))-1),   & !! in
                      cospOUT%modis_Cloud_Particle_Size_Water_Mean(ij+int(CSCAL_MASK_INDICES(:))-1), & !! in
                      cospOUT%modis_Cloud_Fraction_Water_Mean(ij+int(CSCAL_MASK_INDICES(:))-1),      & !! in
                      cospOUT%modis_Ice_Water_Path_Mean(ij+int(CSCAL_MASK_INDICES(:))-1),            & !! in
                      cospOUT%modis_Optical_Thickness_Ice_Mean(ij+int(CSCAL_MASK_INDICES(:))-1),     & !! in
                      cospOUT%modis_Cloud_Particle_Size_Ice_Mean(ij+int(CSCAL_MASK_INDICES(:))-1),   & !! in
                      cospOUT%modis_Cloud_Fraction_Ice_Mean(ij+int(CSCAL_MASK_INDICES(:))-1),        & !! in
                      cospIN%frac_out(int(CSCAL_MASK_INDICES(:)),:,:),                               & !! in
                      cloudsatDBZe,                                                                  & !! in
                      cfodd_ntotal, wr_occfreq_ntotal                                                ) !! inout
             else
                call cosp_diag_warmrain(                                            &
                      cloudsatIN%Npoints, cloudsatIN%Ncolumns, cospIN%Nlevels,       & !! in
                      cospgridIN%at, cospgridIN%hgt_matrix,                          & !! in
                      cospOUT%modis_Liquid_Water_Path_Mean(ij:ik),                   & !! in
                      cospOUT%modis_Optical_Thickness_Water_Mean(ij:ik),             & !! in
                      cospOUT%modis_Cloud_Particle_Size_Water_Mean(ij:ik),           & !! in
                      cospOUT%modis_Cloud_Fraction_Water_Mean(ij:ik),                & !! in
                      cospOUT%modis_Ice_Water_Path_Mean(ij:ik),                      & !! in
                      cospOUT%modis_Optical_Thickness_Ice_Mean(ij:ik),               & !! in
                      cospOUT%modis_Cloud_Particle_Size_Ice_Mean(ij:ik),             & !! in
                      cospOUT%modis_Cloud_Fraction_Ice_Mean(ij:ik),                  & !! in
                      cospIN%frac_out,                                               & !! in
                      cloudsatDBZe,                                                  & !! in
                      cfodd_ntotal, wr_occfreq_ntotal                                ) !! inout
             endif 
          endif  !! use_vgrid or not

          ! Store, when necessary
          if (associated(cospOUT%lidar_only_freq_cloud)) then
             if (cospIN % cospswathsIN(3) % N_inst_swaths .gt. 0) then ! Trigger use of swathed arrays 
                cospOUT%lidar_only_freq_cloud(ij+int(CSCAL_MASK_INDICES(:))-1,:)  = lidar_only_freq_cloud(:,:)
             else
                cospOUT%lidar_only_freq_cloud(ij:ik,:) = lidar_only_freq_cloud
             endif
          endif
          
          if (cospIN % cospswathsIN(3) % N_inst_swaths .gt. 0) then ! Trigger use of swathed arrays for CSCAL
             if (cospIN % cospswathsIN(6) % N_inst_swaths .gt. 0) then ! If MODIS is also swathed then use the joint mask MODIS_CSCAL_MASK_INDICES for setting R_UNDEF
                if ( associated(cospOUT%cfodd_ntotal) ) then
                   cospOUT%cfodd_ntotal(ij+int(CSCAL_MASK_INDICES(:))-1,:,:,:) = cfodd_ntotal
                   cospOUT%cfodd_ntotal(ij+int(MODIS_CSCAL_MASK_INDICES(:))-1,:,:,:) = R_UNDEF
                endif
                if ( associated(cospOUT%wr_occfreq_ntotal) ) then
                   cospOUT%wr_occfreq_ntotal(ij+int(CSCAL_MASK_INDICES(:))-1,:) = wr_occfreq_ntotal
                   cospOUT%wr_occfreq_ntotal(ij+int(MODIS_CSCAL_MASK_INDICES(:))-1,:) = R_UNDEF
                endif         
             else
                if ( associated(cospOUT%cfodd_ntotal) ) then
                   cospOUT%cfodd_ntotal(ij:ik,:,:,:) = R_UNDEF
                   cospOUT%cfodd_ntotal(ij+int(CSCAL_MASK_INDICES(:))-1,:,:,:) = cfodd_ntotal
                endif
                if ( associated(cospOUT%wr_occfreq_ntotal) ) then
                   cospOUT%wr_occfreq_ntotal(ij:ik,:) = R_UNDEF
                   cospOUT%wr_occfreq_ntotal(ij+int(CSCAL_MASK_INDICES(:))-1,:) = wr_occfreq_ntotal
                endif
             endif
          else 
             if ( associated(cospOUT%cfodd_ntotal) ) then
                cospOUT%cfodd_ntotal(ij:ik,:,:,:) = cfodd_ntotal
             endif
             if ( associated(cospOUT%wr_occfreq_ntotal) ) then
                cospOUT%wr_occfreq_ntotal(ij:ik,:) = wr_occfreq_ntotal
             endif
          endif
       else
          if ( associated(cospOUT%cfodd_ntotal) )      cospOUT%cfodd_ntotal(ij:ik,:,:,:)  = R_UNDEF     
          if ( associated(cospOUT%wr_occfreq_ntotal) ) cospOUT%wr_occfreq_ntotal(ij:ik,:) = R_UNDEF
       endif
    endif
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! 7) Cleanup
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (Lisccp_subcolumn .or. Lmodis_subcolumn) then
       nullify(isccpIN%Ncolumns,isccpIN%Nlevels,isccpIN%emsfc_lw,                        &
               isccpIN%skt,isccpIN%qv,isccpIN%at,isccpIN%frac_out,isccpIN%dtau,          &
               isccpIN%dem,isccpIN%phalf,isccpIN%sunlit,isccpIN%pfull)
       call COSP_ASSIGN_isccpIN_CLEAN()
    endif

    if (Lmisr_subcolumn) then
       nullify(misrIN%Ncolumns,misrIN%Nlevels,misrIN%dtau,misrIN%sunlit,  &
               misrIN%zfull,misrIN%at)
       call COSP_ASSIGN_misrIN_CLEAN()
    endif

    if (Lcalipso_subcolumn) then
       nullify(calipsoIN%Ncolumns,calipsoIN%Nlevels,calipsoIN%beta_mol,&
               calipsoIN%betatot,calipsoIN%betatot_liq,calipsoIN%betatot_ice,            &
               calipsoIN%tau_mol,calipsoIN%tautot,calipsoIN%tautot_liq,calipsoIN%tautot_ice)
       if (allocated(CSCAL_MASK_INDICES)) deallocate(CSCAL_MASK_INDICES)
       if (allocated(CSCAL_SWATH_MASK)) deallocate(CSCAL_SWATH_MASK)
       call COSP_ASSIGN_calipsoIN_CLEAN()
    endif

    if (LgrLidar532_subcolumn) then 
       nullify(grLidar532IN%Npoints,grLidar532IN%Ncolumns,grLidar532IN%Nlevels,grLidar532IN%beta_mol, &
               grLidar532IN%betatot,grLidar532IN%tau_mol,grLidar532IN%tautot) 
    endif 

    if (Latlid_subcolumn) then
       nullify(atlidIN%Ncolumns,atlidIN%Nlevels,atlidIN%beta_mol_atlid, &
               atlidIN%betatot_atlid,atlidIN%tau_mol_atlid,atlidIN%tautot_atlid)
       call COSP_ASSIGN_atlidIN_CLEAN()
    endif 

    if (Lparasol_subcolumn) then
       nullify(parasolIN%Nlevels,parasolIN%Ncolumns,parasolIN%Nrefl,   &
               parasolIN%tautot_S_liq,parasolIN%tautot_S_ice)
       call COSP_ASSIGN_parasolIN_CLEAN()
    endif

    if (Lcloudsat_subcolumn) then
       nullify(cloudsatIN%Nlevels,cloudsatIN%Ncolumns,cloudsatIN%rcfg,&
               cloudsatIN%kr_vol,cloudsatIN%g_vol,cloudsatIN%z_vol,cloudsatIN%hgt_matrix)
       call COSP_ASSIGN_cloudsatIN_CLEAN()
       if (allocated(CSCAL_MASK_INDICES)) deallocate(CSCAL_MASK_INDICES)
       if (allocated(CSCAL_SWATH_MASK))   deallocate(CSCAL_SWATH_MASK)
    endif

    if (Lmodis_subcolumn) then
       nullify(modisIN%Npoints,modisIN%Ncolumns,modisIN%Nlevels,modisIN%tau,modisIN%g,   &
               modisIN%liqFrac,modisIN%w0)
       if (allocated(modisIN%sunlit))    deallocate(modisIN%sunlit)
       if (allocated(modisIN%notSunlit)) deallocate(modisIN%notSunlit)
       if (allocated(modisIN%pres))      deallocate(modisIN%pres)
       if (allocated(MODIS_CSCAL_MASK_INDICES)) deallocate(MODIS_CSCAL_MASK_INDICES)
    endif
    
    if (Lrttov_column) then
       nullify(rttovIN%nPoints,rttovIN%nLevels,rttovIN%nSubCols,rttovIN%co2,rttovIN%ch4, &
               rttovIN%n2o,rttovIN%co,rttovIN%h_surf,rttovIN%u_surf,rttovIN%v_surf,      &
               rttovIN%t_skin,rttovIN%p_surf,rttovIN%sfcmask,rttovIN%latitude,           &
               rttovIN%longitude,rttovIN%p,rttovIN%ph,rttovIN%q2m,rttovIN%t2m,rttovIN%t, &
               rttovIN%q,rttovIN%o3,rttovIN%rttov_date,rttovIN%rttov_time,rttovIN%tca,   &
               rttovIN%cldLiq,rttovIN%cldIce,rttovIN%DeffLiq,rttovIN%DeffIce,            &
               rttovIN%fl_rain,rttovIN%fl_snow)
    endif    

    if (allocated(calipso_beta_tot))      deallocate(calipso_beta_tot)
    if (allocated(grLidar532_beta_tot))  deallocate(grLidar532_beta_tot)
    if (allocated(atlid_beta_tot))        deallocate(atlid_beta_tot) 
    if (allocated(calipso_beta_mol))      deallocate(calipso_beta_mol)
    if (allocated(grLidar532_beta_mol))  deallocate(grLidar532_beta_mol)
    if (allocated(atlid_beta_mol))        deallocate(atlid_beta_mol)
    if (allocated(calipso_betaperp_tot))  deallocate(calipso_betaperp_tot)
    if (allocated(cloudsatDBZe))          deallocate(cloudsatDBZe)
    if (allocated(lidar_only_freq_cloud)) deallocate(lidar_only_freq_cloud)
    if (allocated(radar_lidar_tcc))       deallocate(radar_lidar_tcc)
    if (allocated(cloudsat_tcc))          deallocate(cloudsat_tcc)
    if (allocated(cloudsat_tcc2))         deallocate(cloudsat_tcc2)
    if (allocated(cfodd_ntotal))          deallocate(cfodd_ntotal)
    if (allocated(wr_occfreq_ntotal))     deallocate(wr_occfreq_ntotal)

  end function COSP_SIMULATOR
  ! ######################################################################################
  ! SUBROUTINE cosp_init
  ! ######################################################################################
  SUBROUTINE COSP_INIT(Lisccp, Lmodis, Lmisr, Lcloudsat, Lcalipso, LgrLidar532, Latlid, Lparasol, Lrttov,     &
       cloudsat_radar_freq, cloudsat_k2, cloudsat_use_gas_abs, cloudsat_do_ray,          &
       isccp_top_height, isccp_top_height_direction, surface_radar, rcfg, lusevgrid,     &
       luseCSATvgrid, Nvgrid, Nlevels, cloudsat_micro_scheme,                            &
       rttov_Ninstruments, rttov_instrument_namelists,rttov_configs,unitn,debug)

    ! INPUTS
    logical,intent(in)    :: Lisccp,Lmodis,Lmisr,Lcloudsat,Lcalipso,LgrLidar532,Latlid,Lparasol
    logical,intent(inout) :: Lrttov
    integer,intent(in)  :: &
         cloudsat_use_gas_abs,       & !
         cloudsat_do_ray,            & !
         isccp_top_height,           & !
         isccp_top_height_direction, & !
         Nlevels,                    & !
         Nvgrid,                     & ! Number of levels for new L3 grid
         surface_radar,              & !
         rttov_Ninstruments
    real(wp),intent(in) :: &
         cloudsat_radar_freq,        & !
         cloudsat_k2                   !
    logical,intent(in) :: &
         lusevgrid,                  & ! Switch to use different vertical grid
         luseCSATvgrid                 ! Switch to use CLOUDSAT grid spacing for new
                                       ! vertical grid
    character(len=64),intent(in) :: &
       cloudsat_micro_scheme           ! Microphysical scheme used by CLOUDSAT
        
    type(character(len=256)), dimension(rttov_Ninstruments) :: & 
        rttov_instrument_namelists   ! Array of paths to RTTOV instrument namelists
    
    ! OUTPUTS
    type(radar_cfg) :: rcfg
    type(rttov_cfg), dimension(:), allocatable,intent(inout) :: rttov_configs
    
    ! Optional args
    integer,intent(in),Optional :: unitn ! Used for io limits
    logical,intent(in),Optional :: debug
    logical :: verbose = .false.

    ! Local variables
    integer  :: i
    real(wp) :: zstep
    
    ! Init debug print statements:
    if (present(debug)) verbose = debug

    ! Initialize MODIS optical-depth bin boundaries for joint-histogram. (defined in cosp_config.F90)
    if (.not. allocated(modis_histTau)) then
       allocate(modis_histTau(ntau+1),modis_histTauEdges(2,ntau),modis_histTauCenters(ntau))
       numMODISTauBins      = ntau
       modis_histTau        = tau_binBounds
       modis_histTauEdges   = tau_binEdges
       modis_histTauCenters = tau_binCenters
    endif

    ! Set up vertical grid used by CALIPSO and CLOUDSAT L3
    use_vgrid = lusevgrid

    if (use_vgrid) then
      Nlvgrid  = Nvgrid
       allocate(vgrid_zl(Nlvgrid),vgrid_zu(Nlvgrid),vgrid_z(Nlvgrid),dz(Nlvgrid))
       ! CloudSat grid requested
       if (luseCSATvgrid)       zstep = 480._wp
       ! Other grid requested. Constant vertical spacing with top at 20 km
       if (.not. luseCSATvgrid) zstep = 20000._wp/Nvgrid
       do i=1,Nvgrid
          vgrid_zl(Nlvgrid-i+1) = (i-1)*zstep
          vgrid_zu(Nlvgrid-i+1) = i*zstep
       enddo
       vgrid_z = (vgrid_zl+vgrid_zu)/2._wp
       dz = zstep
    else
       Nlvgrid = Nlevels
       allocate(vgrid_zl(Nlvgrid),vgrid_zu(Nlvgrid),vgrid_z(Nlvgrid),dz(Nlvgrid))
       vgrid_zl = 0._wp
       vgrid_zu = 0._wp
       vgrid_z  = 0._wp
       dz       = 0._wp
    endif

    ! Initialize simulators
    if (Lisccp) call cosp_isccp_init(isccp_top_height,isccp_top_height_direction)
    if (Lmodis) call cosp_modis_init()
    if (Lmisr)  call cosp_misr_init()
    
    if (Lrttov) then
        if (present(unitn)) then
            call cosp_rttov_init(Lrttov,Nlevels,rttov_Ninstruments, &
                                 rttov_instrument_namelists,        &
                                 rttov_configs,unitn=unitn,         &
                                 debug=verbose)
        else
            call cosp_rttov_init(Lrttov,Nlevels,rttov_Ninstruments, &
                                 rttov_instrument_namelists,        &
                                 rttov_configs,debug=verbose)
        end if
    endif

    if (Lcloudsat) call cosp_cloudsat_init(cloudsat_radar_freq,cloudsat_k2,              &
         cloudsat_use_gas_abs,cloudsat_do_ray,R_UNDEF,N_HYDRO, surface_radar,            &
         rcfg,cloudsat_micro_scheme)
    if (Lcalipso) call cosp_calipso_init()
    if (LgrLidar532) call cosp_grLidar532_init()
    if (Latlid) call cosp_atlid_init()
    if (Lparasol) call cosp_parasol_init()

    linitialization = .FALSE.
  END SUBROUTINE COSP_INIT


  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE cosp_errorCheck
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine cosp_errorCheck(cospgridIN, cospIN, Lisccp_subcolumn, Lisccp_column,           &
       Lmisr_subcolumn, Lmisr_column, Lmodis_subcolumn, Lmodis_column, Lcloudsat_subcolumn, &
       Lcloudsat_column, Lcalipso_subcolumn, Lcalipso_column, Latlid_subcolumn,             &
       Latlid_column, LgrLidar532_subcolumn, LgrLidar532_column,                            &
       Lrttov_column, Lparasol_subcolumn, Lparasol_column, Lradar_lidar_tcc,                &
       Llidar_only_freq_cloud, Lcloudsat_tcc, Lcloudsat_tcc2, Lcloudsat_modis_wr,           &
       cospOUT, errorMessage, nError)
    
    ! Inputs
    type(cosp_column_inputs),intent(in) :: &
         cospgridIN       ! Model grid inputs to COSP
    type(cosp_optical_inputs),intent(in) :: &
         cospIN           ! Derived (optical) input to COSP
    
    ! Outputs
    logical,intent(inout) :: &
         Lisccp_subcolumn,    & ! ISCCP subcolumn simulator on/off switch
         Lisccp_column,       & ! ISCCP column simulator on/off switch
         Lmisr_subcolumn,     & ! MISR subcolumn simulator on/off switch
         Lmisr_column,        & ! MISR column simulator on/off switch
         Lmodis_subcolumn,    & ! MODIS subcolumn simulator on/off switch
         Lmodis_column,       & ! MODIS column simulator on/off switch
         Lcloudsat_subcolumn, & ! CLOUDSAT subcolumn simulator on/off switch
         Lcloudsat_column,    & ! CLOUDSAT column simulator on/off switch
         Lcalipso_subcolumn,  & ! CALIPSO subcolumn simulator on/off switch
         Lcalipso_column,     & ! CALIPSO column simulator on/off switch
         Latlid_subcolumn,    & ! EarthCare subcolumn simulator on/off switch
         Latlid_column,       & ! EarthCare column simulator on/off switch
         LgrLidar532_subcolumn, & ! Ground Lidar subcolumn simulator on/off switch
         LgrLidar532_column,  & ! Ground Lidar column simulator on/off switch
         Lparasol_subcolumn,  & ! PARASOL subcolumn simulator on/off switch
         Lparasol_column,     & ! PARASOL column simulator on/off switch
         Lrttov_column,       & ! RTTOV column simulator on/off switch
         Lcloudsat_tcc,       & !
         Lcloudsat_tcc2,      & !
         Lradar_lidar_tcc,    & ! On/Off switch for joint Calipso/Cloudsat product
         Llidar_only_freq_cloud, & ! On/Off switch for joint Calipso/Cloudsat product
         Lcloudsat_modis_wr     ! On/Off switch for joint CloudSat/MODIS warm rain product
    type(cosp_outputs),intent(inout) :: &
         cospOUT                ! COSP Outputs
    character(len=256),dimension(100) :: errorMessage
    integer,intent(out) :: nError
    
    ! Local variables
    logical :: alloc_status
    integer :: i
    
    nError = 0
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! PART 0: Ensure that the inputs needed by the requested simulators are allocated.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! ISCCP simulator
    if (Lisccp_subcolumn .or. Lisccp_column) then
       alloc_status = .true.
       if (.not. allocated(cospgridIN%skt)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (ISCCP simulator): cospgridIN%skt has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%qv)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (ISCCP simulator): cospgridIN%qv has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%at)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (ISCCP simulator): cospgridIN%at has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospIN%frac_out)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (ISCCP simulator): cospIN%frac_out has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospIN%tau_067)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (ISCCP simulator): cospIN%tau_067 has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospIN%emiss_11)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (ISCCP simulator): cospIN%emiss_11 has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%phalf)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (ISCCP simulator): cospgridIN%phalf has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%sunlit)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (ISCCP simulator): cospgridIN%sunlit has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%pfull)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (ISCCP simulator): cospgridIN%pfull has not been allocated'
          alloc_status = .false.
       endif
       if (.not. alloc_status) then
          Lisccp_subcolumn = .false.
          Lisccp_column    = .false.
          if (associated(cospOUT%isccp_totalcldarea))  cospOUT%isccp_totalcldarea(:)  = R_UNDEF
          if (associated(cospOUT%isccp_meantb))        cospOUT%isccp_meantb(:)        = R_UNDEF
          if (associated(cospOUT%isccp_meantbclr))     cospOUT%isccp_meantbclr(:)     = R_UNDEF
          if (associated(cospOUT%isccp_meanptop))      cospOUT%isccp_meanptop(:)      = R_UNDEF
          if (associated(cospOUT%isccp_meantaucld))    cospOUT%isccp_meantaucld(:)    = R_UNDEF
          if (associated(cospOUT%isccp_meanalbedocld)) cospOUT%isccp_meanalbedocld(:) = R_UNDEF
          if (associated(cospOUT%isccp_boxtau))        cospOUT%isccp_boxtau(:,:)      = R_UNDEF
          if (associated(cospOUT%isccp_boxptop))       cospOUT%isccp_boxptop(:,:)     = R_UNDEF
          if (associated(cospOUT%isccp_fq))            cospOUT%isccp_fq(:,:,:)        = R_UNDEF
       endif
    endif
    
    ! MISR simulator
    if (Lmisr_subcolumn .or. Lmisr_column) then
       alloc_status = .true.
       if (.not. allocated(cospIN%tau_067)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (MISR simulator): cospIN%tau_067 has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%sunlit)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (MISR simulator): cospgridIN%sunlit has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%hgt_matrix)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (MISR simulator): cospgridIN%hgt_matrix has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%at)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (MISR simulator): cospgridIN%at has not been allocated'
          alloc_status = .false.
       endif
       if (.not. alloc_status) then
          Lmisr_subcolumn = .false.
          Lmisr_column    = .false.
          if (associated(cospOUT%misr_fq))                   cospOUT%misr_fq(:,:,:)                 = R_UNDEF
          if (associated(cospOUT%misr_dist_model_layertops)) cospOUT%misr_dist_model_layertops(:,:) = R_UNDEF
          if (associated(cospOUT%misr_meanztop))             cospOUT%misr_meanztop(:)               = R_UNDEF
          if (associated(cospOUT%misr_cldarea))              cospOUT%misr_cldarea(:)                = R_UNDEF          
       endif
    endif

    ! EarthCare Lidar simulator.
    if (Latlid_subcolumn .or. Latlid_column) then
       alloc_status = .true.
       if (.not. allocated(cospIN%beta_mol_atlid)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (EarthCare Lidar simulator): cospIN%beta_mol_atlid has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospIN%betatot_atlid)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (EarthCare Lidar simulator): cospIN%betatot_atlid has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospIN%tau_mol_atlid)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (EarthCare Lidar simulator): cospIN%tau_mol_atlid has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospIN%tautot_atlid)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (EarthCare Lidar simulator): cospIN%tautot_atlid has not been allocated'
          alloc_status = .false.
       endif
       if (.not. alloc_status) then
          Latlid_subcolumn = .false.
          Latlid_column    = .false.
          if (associated(cospOUT%atlid_cfad_sr))  cospOUT%atlid_cfad_sr(:,:,:)  = R_UNDEF
          if (associated(cospOUT%atlid_lidarcld)) cospOUT%atlid_lidarcld(:,:)   = R_UNDEF
          if (associated(cospOUT%atlid_cldlayer)) cospOUT%atlid_cldlayer(:,:)   = R_UNDEF
          if (associated(cospOUT%atlid_beta_mol)) cospOUT%atlid_beta_mol(:,:)   = R_UNDEF
          if (associated(cospOUT%atlid_beta_tot)) cospOUT%atlid_beta_tot(:,:,:) = R_UNDEF
       endif
       
       ! EarthCare column simulator requires additional inputs not required by the subcolumn simulator.
       alloc_status = .true.
       if (.not. allocated(cospgridIN%hgt_matrix)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (EarthCare Lidar simulator): cospgridIN%hgt_matrix has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%hgt_matrix_half)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (EarthCare Lidar simulator): cospgridIN%hgt_matrix_half has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%phalf)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (EarthCare Lidar simulator): cospgridIN%phalf has not been allocated'
          alloc_status = .false.
       endif
       if (.not. alloc_status) then
          Latlid_column  = .false.
          if (associated(cospOUT%atlid_cfad_sr)) cospOUT%atlid_cfad_sr(:,:,:) = R_UNDEF
       endif
    endif

    ! Ground Lidar simulator.
    if (LgrLidar532_subcolumn .or. LgrLidar532_column) then
       alloc_status = .true.
       if (.not. allocated(cospIN%beta_mol_grLidar532)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (Ground Lidar simulator): cospIN%beta_mol_grLidar532 has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospIN%betatot_grLidar532)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (Ground Lidar simulator): cospIN%betatot_grLidar532 has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospIN%tau_mol_grLidar532)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (Ground Lidar simulator): cospIN%tau_mol_grLidar532 has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospIN%tautot_grLidar532)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (Ground Lidar simulator): cospIN%tautot_grLidar532 has not been allocated'
          alloc_status = .false.
       endif
       if (.not. alloc_status) then
          LgrLidar532_subcolumn = .false.
          LgrLidar532_column    = .false.
          if (associated(cospOUT%grLidar532_cfad_sr))  cospOUT%grLidar532_cfad_sr(:,:,:)  = R_UNDEF
          if (associated(cospOUT%grLidar532_lidarcld)) cospOUT%grLidar532_lidarcld(:,:)   = R_UNDEF
          if (associated(cospOUT%grLidar532_cldlayer)) cospOUT%grLidar532_cldlayer(:,:)   = R_UNDEF
          if (associated(cospOUT%grLidar532_beta_mol)) cospOUT%grLidar532_beta_mol(:,:)   = R_UNDEF
          if (associated(cospOUT%grLidar532_beta_tot)) cospOUT%grLidar532_beta_tot(:,:,:) = R_UNDEF
       endif
       
       ! Ground Lidar column simulator requires additional inputs not required by the subcolumn simulator.
       alloc_status = .true.
       if (.not. allocated(cospgridIN%hgt_matrix)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (Ground Lidar simulator): cospgridIN%hgt_matrix has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%hgt_matrix_half)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (Ground Lidar simulator): cospgridIN%hgt_matrix_half has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%phalf)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (Ground Lidar simulator): cospgridIN%phalf has not been allocated'
          alloc_status = .false.
       endif
       if (.not. alloc_status) then
          LgrLidar532_column  = .false.
          if (associated(cospOUT%grLidar532_cfad_sr)) cospOUT%grLidar532_cfad_sr(:,:,:) = R_UNDEF
       endif
    endif
    
    ! Calipso Lidar simulator
    if (Lcalipso_subcolumn .or. Lcalipso_column) then
       alloc_status = .true.
       if (.not. allocated(cospIN%beta_mol_calipso)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (Calipso Lidar simulator): cospIN%beta_mol_calipso has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospIN%betatot_calipso)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (Calipso Lidar simulator): cospIN%betatot_calipso has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospIN%betatot_liq_calipso)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (Calipso Lidar simulator):'//&
               ' cospIN%betatot_liq_calipso has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospIN%betatot_ice_calipso)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (Calipso Lidar simulator):'//&
               ' cospIN%betatot_ice_calipso has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospIN%tau_mol_calipso)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (Calipso Lidar simulator): cospIN%tau_mol_calipso has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospIN%tautot_calipso)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (Calipso Lidar simulator): cospIN%tautot_calipso has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospIN%tautot_liq_calipso)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (Calipso Lidar simulator):'//&
               ' cospIN%tautot_liq has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospIN%tautot_ice_calipso)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (Calipso Lidar simulator):'//&
               ' cospIN%tautot_ice has not been allocated'
          alloc_status = .false.
       endif
       if (.not. alloc_status) then
          Lcalipso_subcolumn = .false.
          Lcalipso_column    = .false.
          if (associated(cospOUT%calipso_cfad_sr))       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF
          if (associated(cospOUT%calipso_lidarcld))      cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldphase)) cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
          if (associated(cospOUT%calipso_cldlayer))      cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
          if (associated(cospOUT%calipso_cldlayerphase)) cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldtmp))   cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF
          if (associated(cospOUT%calipso_beta_mol))      cospOUT%calipso_beta_mol(:,:)        = R_UNDEF
          if (associated(cospOUT%calipso_beta_tot))      cospOUT%calipso_beta_tot(:,:,:)      = R_UNDEF
          if (associated(cospOUT%calipso_betaperp_tot))  cospOUT%calipso_betaperp_tot(:,:,:)  = R_UNDEF
          ! Also, turn-off joint-products 
          if (Lradar_lidar_tcc) then
             Lradar_lidar_tcc = .false.
             if (associated(cospOUT%radar_lidar_tcc)) cospOUT%radar_lidar_tcc(:) = R_UNDEF
          endif
          if (Lcloudsat_tcc) then
             Lcloudsat_tcc = .false.
             if (associated(cospOUT%cloudsat_tcc)) cospOUT%cloudsat_tcc(:) = R_UNDEF
          endif
          if (Lcloudsat_tcc2) then
             Lcloudsat_tcc2 = .false.
             if (associated(cospOUT%cloudsat_tcc2)) cospOUT%cloudsat_tcc2(:) = R_UNDEF
          endif
          if (Llidar_only_freq_cloud) then
             Llidar_only_freq_cloud = .false.
             if (associated(cospOUT%lidar_only_freq_cloud)) cospOUT%lidar_only_freq_cloud(:,:) = R_UNDEF
          endif
       endif
       
       ! Calipso column simulator requires additional inputs not required by the subcolumn simulator.
       alloc_status = .true.
       if (.not. allocated(cospgridIN%hgt_matrix)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (Calipso Lidar simulator):'//&
               ' cospgridIN%hgt_matrix has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%hgt_matrix_half)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (Calipso Lidar simulator):'//&
               ' cospgridIN%hgt_matrix_half has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%at)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (Calipso Lidar simulator): cospgridIN%at has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%surfelev)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (Calipso Lidar simulator): cospgridIN%surfelev has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%phalf)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (Calipso Lidar simulator): cospgridIN%phalf has not been allocated'
          alloc_status = .false.
       endif
       if (.not. alloc_status) then
          Lcalipso_column  = .false.
          if (associated(cospOUT%calipso_cfad_sr))       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF
          if (associated(cospOUT%calipso_lidarcld))      cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldphase)) cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
          if (associated(cospOUT%calipso_cldlayer))      cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
          if (associated(cospOUT%calipso_cldlayerphase)) cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldtmp))   cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF
          if (Lcloudsat_tcc) then
             Lcloudsat_tcc = .false.
             if (associated(cospOUT%cloudsat_tcc)) cospOUT%cloudsat_tcc(:) = R_UNDEF
          endif
          if (Lcloudsat_tcc2) then
             Lcloudsat_tcc2 = .false.
             if (associated(cospOUT%cloudsat_tcc2)) cospOUT%cloudsat_tcc2(:) = R_UNDEF
          endif
          ! Also, turn-off joint-products 
          if (Lradar_lidar_tcc) then
             Lradar_lidar_tcc = .false.
             if (associated(cospOUT%radar_lidar_tcc)) cospOUT%radar_lidar_tcc(:) = R_UNDEF
          endif
          if (Llidar_only_freq_cloud) then
             Llidar_only_freq_cloud = .false.
             if (associated(cospOUT%lidar_only_freq_cloud)) cospOUT%lidar_only_freq_cloud(:,:) = R_UNDEF
          endif
       endif
    endif
    
    ! PARASOL simulator
    if (Lparasol_subcolumn .or. Lparasol_column) then
       alloc_status = .true.
       if (.not. allocated(cospIN%tautot_S_liq)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (PARASOL simulator): cospIN%tautot_S_liq has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospIN%tautot_S_ice)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (PARASOL simulator): cospIN%tautot_S_ice has not been allocated'
          alloc_status = .false.
       endif
       if (.not. alloc_status) then
          Lparasol_subcolumn  = .false.
          Lparasol_column     = .false.
          if (associated(cospOUT%parasolPix_refl))  cospOUT%parasolPix_refl(:,:,:) = R_UNDEF
          if (associated(cospOUT%parasolGrid_refl)) cospOUT%parasolGrid_refl(:,:)  = R_UNDEF
       endif
       
       ! PARASOL column simulator requires additional inputs not required by the subcolumn simulator.
       alloc_status = .true.
       if (.not. allocated(cospgridIN%land)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (PARASOL simulator): cospgridIN%land has not been allocated'
          alloc_status = .false.
       endif
       if (.not. alloc_status) then
          Lparasol_column  = .false.
          if (associated(cospOUT%parasolGrid_refl)) cospOUT%parasolGrid_refl(:,:)  = R_UNDEF
       endif
    endif
    
    ! Cloudsat radar simulator
    if (Lcloudsat_subcolumn .or. Lcloudsat_column) then
       alloc_status = .true.
       if (.not. allocated(cospIN%z_vol_cloudsat)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (Cloudsat radar simulator):'//&
               ' cospIN%z_vol_cloudsat has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospIN%kr_vol_cloudsat)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (Cloudsat radar simulator):'//&
               ' cospIN%kr_vol_cloudsat has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospIN%g_vol_cloudsat)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (Cloudsat radar simulator):'//&
               ' cospIN%g_vol_cloudsat has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%hgt_matrix)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (Cloudsat radar simulator):'//&
               ' cospgridIN%hgt_matrix has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%surfelev)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (Cloudsat radar simulator):'//&
               ' cospgridIN%surfelev has not been allocated'
          alloc_status = .false.
       endif
       if (.not. alloc_status) then
          Lcloudsat_subcolumn  = .false.
          Lcloudsat_column     = .false.
          if (associated(cospOUT%cloudsat_cfad_ze)) cospOUT%cloudsat_cfad_ze(:,:,:) = R_UNDEF
          if (associated(cospOUT%cloudsat_Ze_tot))  cospOUT%cloudsat_Ze_tot(:,:,:)  = R_UNDEF
          if (Lcloudsat_tcc) then
             Lcloudsat_tcc = .false.
             if (associated(cospOUT%cloudsat_tcc)) cospOUT%cloudsat_tcc(:) = R_UNDEF
          endif
          if (Lcloudsat_tcc2) then
             Lcloudsat_tcc2 = .false.
             if (associated(cospOUT%cloudsat_tcc2)) cospOUT%cloudsat_tcc2(:) = R_UNDEF
          endif
          ! Also, turn-off joint-products 
          if (Lradar_lidar_tcc) then
             Lradar_lidar_tcc = .false.
             if (associated(cospOUT%radar_lidar_tcc)) cospOUT%radar_lidar_tcc(:) = R_UNDEF
          endif
          if (Llidar_only_freq_cloud) then
             Llidar_only_freq_cloud = .false.
             if (associated(cospOUT%lidar_only_freq_cloud)) cospOUT%lidar_only_freq_cloud(:,:) = R_UNDEF
          endif
          if (Lcloudsat_modis_wr) then
             Lcloudsat_modis_wr = .false.
             if (associated(cospOUT%cfodd_ntotal)) cospOUT%cfodd_ntotal(:,:,:,:) = R_UNDEF
             if (associated(cospOUT%wr_occfreq_ntotal)) cospOUT%wr_occfreq_ntotal(:,:) = R_UNDEF
          endif
       endif
       
       ! Cloudsat column simulator requires additional inputs not required by the subcolumn simulator.
       alloc_status = .true.
       if (.not. allocated(cospgridIN%hgt_matrix_half)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (Cloudsat radar simulator):'//&
               ' cospgridIN%hgt_matrix_half has not been allocated'
          alloc_status = .false.
       endif
       if (.not. alloc_status) then
          Lcloudsat_column  = .false.
          if (associated(cospOUT%cloudsat_cfad_ze))      cospOUT%cloudsat_cfad_ze(:,:,:)    = R_UNDEF
          if (Lcloudsat_tcc) then
             Lcloudsat_tcc = .false.
             if (associated(cospOUT%cloudsat_tcc)) cospOUT%cloudsat_tcc(:) = R_UNDEF
          endif
          if (Lcloudsat_tcc2) then
             Lcloudsat_tcc2 = .false.
             if (associated(cospOUT%cloudsat_tcc2)) cospOUT%cloudsat_tcc2(:) = R_UNDEF
          endif
          ! Also, turn-off joint-products 
          if (Lradar_lidar_tcc) then
             Lradar_lidar_tcc = .false.
             if (associated(cospOUT%radar_lidar_tcc)) cospOUT%radar_lidar_tcc(:) = R_UNDEF
          endif
          if (Llidar_only_freq_cloud) then
             Llidar_only_freq_cloud = .false.
             if (associated(cospOUT%lidar_only_freq_cloud)) cospOUT%lidar_only_freq_cloud(:,:) = R_UNDEF
          endif
          if (Lcloudsat_modis_wr) then
             Lcloudsat_modis_wr = .false.
             if (associated(cospOUT%cfodd_ntotal)) cospOUT%cfodd_ntotal(:,:,:,:) = R_UNDEF
             if (associated(cospOUT%wr_occfreq_ntotal)) cospOUT%wr_occfreq_ntotal(:,:) = R_UNDEF
          endif
       endif
    endif
    
    ! MODIS simulator
    if (Lmodis_subcolumn .or. Lmodis_column) then
       alloc_status = .true.
       if (.not. allocated(cospIN%fracLiq)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (MODIS simulator): cospIN%fracLiq has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospIN%tau_067)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (MODIS simulator): cospIN%tau_067 has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospIN%asym)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (MODIS simulator): cospIN%asym has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospIN%ss_alb)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (MODIS simulator): cospIN%ss_alb has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%sunlit)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (MODIS simulator): cospgridIN%sunlit has not been allocated'
          alloc_status = .false.
       endif
       if (.not. alloc_status) then
          Lmodis_subcolumn  = .false.
          Lmodis_column     = .false.
          if (associated(cospOUT%modis_Cloud_Fraction_Total_Mean))                          &
               cospOUT%modis_Cloud_Fraction_Total_Mean(:)                   = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_Water_Mean))                          &
               cospOUT%modis_Cloud_Fraction_Water_Mean(:)                   = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_Ice_Mean))                            &
               cospOUT%modis_Cloud_Fraction_Ice_Mean(:)                     = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_High_Mean))                           &
               cospOUT%modis_Cloud_Fraction_High_Mean(:)                    = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_Mid_Mean))                            &
               cospOUT%modis_Cloud_Fraction_Mid_Mean(:)                     = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_Low_Mean))                            &
               cospOUT%modis_Cloud_Fraction_Low_Mean(:)                     = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Total_Mean))                       &
               cospOUT%modis_Optical_Thickness_Total_Mean(:)                = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Water_Mean))                       &
               cospOUT%modis_Optical_Thickness_Water_Mean(:)                = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Ice_Mean))                         &
               cospOUT%modis_Optical_Thickness_Ice_Mean(:)                  = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Total_LogMean))                    &
               cospOUT%modis_Optical_Thickness_Total_LogMean(:)             = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Water_LogMean))                    &
               cospOUT%modis_Optical_Thickness_Water_LogMean(:)             = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Ice_LogMean))                      &
               cospOUT%modis_Optical_Thickness_Ice_LogMean(:)               = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Particle_Size_Water_Mean))                     &
               cospOUT%modis_Cloud_Particle_Size_Water_Mean(:)              = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Particle_Size_Ice_Mean))                       &
               cospOUT%modis_Cloud_Particle_Size_Ice_Mean(:)                = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Top_Pressure_Total_Mean))                      &
               cospOUT%modis_Cloud_Top_Pressure_Total_Mean(:)               = R_UNDEF
          if (associated(cospOUT%modis_Liquid_Water_Path_Mean))                             &
               cospOUT%modis_Liquid_Water_Path_Mean(:)                      = R_UNDEF
          if (associated(cospOUT%modis_Ice_Water_Path_Mean))                                &
               cospOUT%modis_Ice_Water_Path_Mean(:)                         = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure))            &
               cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure(:,:,:) = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Liq))            &
               cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Liq(:,:,:) = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Ice))            &
               cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Ice(:,:,:) = R_UNDEF
          if (associated(cospOUT%modis_LWP_vs_ReffLIQ))            &
               cospOUT%modis_LWP_vs_ReffLIQ(:,:,:) = R_UNDEF
          if (associated(cospOUT%modis_IWP_vs_ReffICE))            &
               cospOUT%modis_IWP_vs_ReffICE(:,:,:) = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_vs_ReffICE))                       &
               cospOUT%modis_Optical_Thickness_vs_ReffICE(:,:,:)            = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_vs_ReffLIQ))                       &
               cospOUT%modis_Optical_Thickness_vs_ReffLIQ(:,:,:)            = R_UNDEF          
          ! Also, turn-off joint-products 
          if (Lcloudsat_modis_wr) then
             Lcloudsat_modis_wr = .false.
             if (associated(cospOUT%cfodd_ntotal)) cospOUT%cfodd_ntotal(:,:,:,:) = R_UNDEF
             if (associated(cospOUT%wr_occfreq_ntotal)) cospOUT%wr_occfreq_ntotal(:,:) = R_UNDEF
          endif
       endif
    endif
    
    ! RTTOV
    if (Lrttov_column) then
       alloc_status = .true.
!       if (.not. allocated(cospgridIN%emis_in)) then
!          nError=nError+1
!          errorMessage(nError) = 'ERROR: COSP input variable (RTTOV): cospgridIN%emis_in has not been allocated'
!          alloc_status = .false.
!       endif
!       if (.not. allocated(cospgridIN%refl_in)) then
!          nError=nError+1
!          errorMessage(nError) = 'ERROR: COSP input variable (RTTOV): cospgridIN%refl_in has not been allocated'
!          alloc_status = .false.
!       endif
       if (.not. allocated(cospgridIN%hgt_matrix_half)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (RTTOV): cospgridIN%hgt_matrix_half has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%u_sfc)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (RTTOV): cospgridIN%u_sfc has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%v_sfc)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (RTTOV): cospgridIN%v_sfc has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%skt)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (RTTOV): cospgridIN%skt has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%phalf)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (RTTOV): cospgridIN%phalf has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%qv)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (RTTOV): cospgridIN%qv has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%at)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (RTTOV): cospgridIN%at has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%rttov_sfcmask)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (RTTOV simulator): cospgridIN%rttov_sfcmask has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%lat)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (RTTOV): cospgridIN%lat has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%lon)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (RTTOV): cospgridIN%lon has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%pfull)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (RTTOV): cospgridIN%pfull has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%phalf)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (RTTOV): cospgridIN%phalf has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%at)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (RTTOV): cospgridIN%at has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%qv)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (RTTOV): cospgridIN%qv has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%o3)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (RTTOV): cospgridIN%o3 has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%tca)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (RTTOV): cospgridIN%tca has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%sza)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (RTTOV): cospgridIN%sza has not been allocated'
          alloc_status = .false.
       endif        
       if (.not. allocated(cospgridIN%cloudIce)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (RTTOV): cospgridIN%cloudIce has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%cloudLiq)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (RTTOV): cospgridIN%cloudLiq has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%fl_rain)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (RTTOV): cospgridIN%fl_rain has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%fl_snow)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (RTTOV): cospgridIN%fl_snow has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%rttov_date)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (RTTOV): cospgridIN%rttov_date has not been allocated'
          alloc_status = .false.
       endif
       if (.not. allocated(cospgridIN%rttov_time)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable (RTTOV): cospgridIN%rttov_time has not been allocated'
          alloc_status = .false.
       endif       
       if (.not. alloc_status) then
          Lrttov_column     = .false.
         if (allocated(cospOUT%rttov_outputs)) then
             do i=1,cospOUT % Ninst_rttov ! Iterate over each instrument
                 if (associated(cospOUT%rttov_outputs(i)%channel_indices))      cospOUT%rttov_outputs(i)%channel_indices(:)  = 0
                 if (associated(cospOUT%rttov_outputs(i)%bt_total))             cospOUT%rttov_outputs(i)%bt_total(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_clear))             cospOUT%rttov_outputs(i)%bt_clear(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total))            cospOUT%rttov_outputs(i)%rad_total(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_clear))            cospOUT%rttov_outputs(i)%rad_clear(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_cloudy))           cospOUT%rttov_outputs(i)%rad_cloudy(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_total))           cospOUT%rttov_outputs(i)%refl_total(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_clear))           cospOUT%rttov_outputs(i)%refl_clear(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_total_pc))          cospOUT%rttov_outputs(i)%bt_total_pc(:,:)    = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total_pc))         cospOUT%rttov_outputs(i)%rad_total_pc(:,:)   = R_UNDEF                
             end do
         end if          
       endif
    endif

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! PART 1: Check input array values for out-of-bounds values. When an out-of-bound value
    !         is encountered, COSP outputs that are dependent on that input are filled with
    !         an undefined value (set in cosp_config.f90) and if necessary, that simulator
    !         is turned off.
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if (any([Lisccp_subcolumn, Lisccp_column, Lmisr_subcolumn, Lmisr_column, &
             Lmodis_subcolumn, Lmodis_column, Lcloudsat_modis_wr])) then
       if (any(cospgridIN%sunlit .lt. 0)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospgridIN%sunlit contains values out of range (0 or 1)'
          Lisccp_subcolumn = .false.
          Lisccp_column    = .false.
          Lmisr_subcolumn  = .false.
          Lmisr_column     = .false.
          Lmodis_subcolumn = .false.
          Lmodis_column    = .false.
          Lcloudsat_modis_wr = .false.
          if (associated(cospOUT%isccp_totalcldarea))  cospOUT%isccp_totalcldarea(:)  = R_UNDEF
          if (associated(cospOUT%isccp_meantb))        cospOUT%isccp_meantb(:)        = R_UNDEF
          if (associated(cospOUT%isccp_meantbclr))     cospOUT%isccp_meantbclr(:)     = R_UNDEF
          if (associated(cospOUT%isccp_meanptop))      cospOUT%isccp_meanptop(:)      = R_UNDEF
          if (associated(cospOUT%isccp_meantaucld))    cospOUT%isccp_meantaucld(:)    = R_UNDEF
          if (associated(cospOUT%isccp_meanalbedocld)) cospOUT%isccp_meanalbedocld(:) = R_UNDEF
          if (associated(cospOUT%isccp_boxtau))        cospOUT%isccp_boxtau(:,:)      = R_UNDEF
          if (associated(cospOUT%isccp_boxptop))       cospOUT%isccp_boxptop(:,:)     = R_UNDEF
          if (associated(cospOUT%isccp_fq))            cospOUT%isccp_fq(:,:,:)        = R_UNDEF
          if (associated(cospOUT%misr_fq))                   cospOUT%misr_fq(:,:,:)                 = R_UNDEF
          if (associated(cospOUT%misr_dist_model_layertops)) cospOUT%misr_dist_model_layertops(:,:) = R_UNDEF
          if (associated(cospOUT%misr_meanztop))             cospOUT%misr_meanztop(:)               = R_UNDEF
          if (associated(cospOUT%misr_cldarea))              cospOUT%misr_cldarea(:)                = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_Total_Mean))                          &
               cospOUT%modis_Cloud_Fraction_Total_Mean(:)                   = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_Water_Mean))                          &
               cospOUT%modis_Cloud_Fraction_Water_Mean(:)                   = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_Ice_Mean))                            &
               cospOUT%modis_Cloud_Fraction_Ice_Mean(:)                     = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_High_Mean))                           &
               cospOUT%modis_Cloud_Fraction_High_Mean(:)                    = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_Mid_Mean))                            &
               cospOUT%modis_Cloud_Fraction_Mid_Mean(:)                     = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_Low_Mean))                            &
               cospOUT%modis_Cloud_Fraction_Low_Mean(:)                     = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Total_Mean))                       &
               cospOUT%modis_Optical_Thickness_Total_Mean(:)                = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Water_Mean))                       &
               cospOUT%modis_Optical_Thickness_Water_Mean(:)                = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Ice_Mean))                         &
               cospOUT%modis_Optical_Thickness_Ice_Mean(:)                  = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Total_LogMean))                    &
               cospOUT%modis_Optical_Thickness_Total_LogMean(:)             = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Water_LogMean))                    &
               cospOUT%modis_Optical_Thickness_Water_LogMean(:)             = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Ice_LogMean))                      &
               cospOUT%modis_Optical_Thickness_Ice_LogMean(:)               = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Particle_Size_Water_Mean))                     &
               cospOUT%modis_Cloud_Particle_Size_Water_Mean(:)              = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Particle_Size_Ice_Mean))                       &
               cospOUT%modis_Cloud_Particle_Size_Ice_Mean(:)                = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Top_Pressure_Total_Mean))                      &
               cospOUT%modis_Cloud_Top_Pressure_Total_Mean(:)               = R_UNDEF
          if (associated(cospOUT%modis_Liquid_Water_Path_Mean))                             &
               cospOUT%modis_Liquid_Water_Path_Mean(:)                      = R_UNDEF
          if (associated(cospOUT%modis_Ice_Water_Path_Mean))                                &
               cospOUT%modis_Ice_Water_Path_Mean(:)                         = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure))            &
               cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure(:,:,:) = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Liq))            &
               cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Liq(:,:,:) = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Ice))            &
               cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Ice(:,:,:) = R_UNDEF
          if (associated(cospOUT%modis_LWP_vs_ReffLIQ))            &
               cospOUT%modis_LWP_vs_ReffLIQ(:,:,:) = R_UNDEF
          if (associated(cospOUT%modis_IWP_vs_ReffICE))            &
               cospOUT%modis_IWP_vs_ReffICE(:,:,:) = R_UNDEF    
          if (associated(cospOUT%modis_Optical_Thickness_vs_ReffICE))                       &
               cospOUT%modis_Optical_Thickness_vs_ReffICE(:,:,:)            = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_vs_ReffLIQ))                       &
               cospOUT%modis_Optical_Thickness_vs_ReffLIQ(:,:,:)            = R_UNDEF
          if (associated(cospOUT%cfodd_ntotal)) cospOUT%cfodd_ntotal(:,:,:,:) = R_UNDEF
          if (associated(cospOUT%wr_occfreq_ntotal)) cospOUT%wr_occfreq_ntotal(:,:) = R_UNDEF
       endif
    endif

    if (any([Lisccp_subcolumn, Lisccp_column, Lmisr_subcolumn, Lmisr_column, Lrttov_column,&
         Lcalipso_column, Lcloudsat_column, Lradar_lidar_tcc,Llidar_only_freq_cloud, &
         Lcloudsat_tcc, Lcloudsat_tcc2, Lcloudsat_modis_wr])) then
       if (any(cospgridIN%at .lt. 0)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospgridIN%at contains values out of range (at<0), expected units (K)'
          Lisccp_subcolumn = .false.
          Lisccp_column    = .false.
          Lmisr_subcolumn  = .false.
          Lmisr_column     = .false.
          Lrttov_column    = .false.
          Lcalipso_column  = .false.
          Lcloudsat_column = .false.
          Lradar_lidar_tcc = .false.
          Llidar_only_freq_cloud = .false.
          Lcloudsat_tcc    = .false.
          Lcloudsat_tcc2   = .false.
          Lcloudsat_modis_wr = .false.
          if (allocated(cospOUT%rttov_outputs)) then
             do i=1,cospOUT % Ninst_rttov ! Iterate over each instrument
                 if (associated(cospOUT%rttov_outputs(i)%channel_indices))      cospOUT%rttov_outputs(i)%channel_indices(:)  = 0
                 if (associated(cospOUT%rttov_outputs(i)%bt_total))             cospOUT%rttov_outputs(i)%bt_total(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_clear))             cospOUT%rttov_outputs(i)%bt_clear(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total))            cospOUT%rttov_outputs(i)%rad_total(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_clear))            cospOUT%rttov_outputs(i)%rad_clear(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_cloudy))           cospOUT%rttov_outputs(i)%rad_cloudy(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_total))           cospOUT%rttov_outputs(i)%refl_total(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_clear))           cospOUT%rttov_outputs(i)%refl_clear(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_total_pc))          cospOUT%rttov_outputs(i)%bt_total_pc(:,:)    = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total_pc))         cospOUT%rttov_outputs(i)%rad_total_pc(:,:)   = R_UNDEF                
             end do
          end if
          if (associated(cospOUT%isccp_totalcldarea))  cospOUT%isccp_totalcldarea(:)  = R_UNDEF
          if (associated(cospOUT%isccp_meantb))        cospOUT%isccp_meantb(:)        = R_UNDEF
          if (associated(cospOUT%isccp_meantbclr))     cospOUT%isccp_meantbclr(:)     = R_UNDEF
          if (associated(cospOUT%isccp_meanptop))      cospOUT%isccp_meanptop(:)      = R_UNDEF
          if (associated(cospOUT%isccp_meantaucld))    cospOUT%isccp_meantaucld(:)    = R_UNDEF
          if (associated(cospOUT%isccp_meanalbedocld)) cospOUT%isccp_meanalbedocld(:) = R_UNDEF
          if (associated(cospOUT%isccp_boxtau))        cospOUT%isccp_boxtau(:,:)      = R_UNDEF
          if (associated(cospOUT%isccp_boxptop))       cospOUT%isccp_boxptop(:,:)     = R_UNDEF
          if (associated(cospOUT%isccp_fq))            cospOUT%isccp_fq(:,:,:)        = R_UNDEF
          if (associated(cospOUT%misr_fq))                   cospOUT%misr_fq(:,:,:)                 = R_UNDEF
          if (associated(cospOUT%misr_dist_model_layertops)) cospOUT%misr_dist_model_layertops(:,:) = R_UNDEF
          if (associated(cospOUT%misr_meanztop))             cospOUT%misr_meanztop(:)               = R_UNDEF
          if (associated(cospOUT%misr_cldarea))              cospOUT%misr_cldarea(:)                = R_UNDEF
          if (associated(cospOUT%calipso_cfad_sr))       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF
          if (associated(cospOUT%calipso_lidarcld))      cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldphase)) cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
          if (associated(cospOUT%calipso_cldlayer))      cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
          if (associated(cospOUT%calipso_cldlayerphase)) cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldtmp))   cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldtype))  cospOUT%calipso_lidarcldtype(:,:,:)  = R_UNDEF
          if (associated(cospOUT%calipso_cldtype))       cospOUT%calipso_cldtype(:,:)         = R_UNDEF 
          if (associated(cospOUT%calipso_cldtypetemp))   cospOUT%calipso_cldtypetemp(:,:)     = R_UNDEF 
          if (associated(cospOUT%calipso_cldtypemeanz))  cospOUT%calipso_cldtypemeanz(:,:)    = R_UNDEF 
          if (associated(cospOUT%calipso_cldtypemeanzse)) cospOUT%calipso_cldtypemeanzse(:,:) = R_UNDEF 
          if (associated(cospOUT%calipso_cldthinemis))   cospOUT%calipso_cldthinemis(:)       = R_UNDEF
          if (associated(cospOUT%cloudsat_cfad_ze))      cospOUT%cloudsat_cfad_ze(:,:,:)      = R_UNDEF
          if (associated(cospOUT%lidar_only_freq_cloud)) cospOUT%lidar_only_freq_cloud(:,:)   = R_UNDEF
          if (associated(cospOUT%radar_lidar_tcc))       cospOUT%radar_lidar_tcc(:)           = R_UNDEF
          if (associated(cospOUT%cloudsat_tcc)) cospOUT%cloudsat_tcc(:) = R_UNDEF
          if (associated(cospOUT%cloudsat_tcc2)) cospOUT%cloudsat_tcc2(:) = R_UNDEF
          if (associated(cospOUT%cfodd_ntotal)) cospOUT%cfodd_ntotal(:,:,:,:) = R_UNDEF
          if (associated(cospOUT%wr_occfreq_ntotal)) cospOUT%wr_occfreq_ntotal(:,:) = R_UNDEF
       endif
    endif
    if (any([Lisccp_subcolumn, Lisccp_column, Lrttov_column])) then
       if (any(cospgridIN%pfull .lt. 0)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospgridIN%pfull contains values out of range'
          Lisccp_subcolumn = .false.
          Lisccp_column    = .false.
          Lrttov_column    = .false.
          if (allocated(cospOUT%rttov_outputs)) then
             do i=1,cospOUT % Ninst_rttov ! Iterate over each instrument
                 if (associated(cospOUT%rttov_outputs(i)%channel_indices))      cospOUT%rttov_outputs(i)%channel_indices(:)  = 0
                 if (associated(cospOUT%rttov_outputs(i)%bt_total))             cospOUT%rttov_outputs(i)%bt_total(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_clear))             cospOUT%rttov_outputs(i)%bt_clear(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total))            cospOUT%rttov_outputs(i)%rad_total(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_clear))            cospOUT%rttov_outputs(i)%rad_clear(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_cloudy))           cospOUT%rttov_outputs(i)%rad_cloudy(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_total))           cospOUT%rttov_outputs(i)%refl_total(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_clear))           cospOUT%rttov_outputs(i)%refl_clear(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_total_pc))          cospOUT%rttov_outputs(i)%bt_total_pc(:,:)    = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total_pc))         cospOUT%rttov_outputs(i)%rad_total_pc(:,:)   = R_UNDEF                
             end do
          end if
          if (associated(cospOUT%isccp_totalcldarea))  cospOUT%isccp_totalcldarea(:)  = R_UNDEF
          if (associated(cospOUT%isccp_meantb))        cospOUT%isccp_meantb(:)        = R_UNDEF
          if (associated(cospOUT%isccp_meantbclr))     cospOUT%isccp_meantbclr(:)     = R_UNDEF
          if (associated(cospOUT%isccp_meanptop))      cospOUT%isccp_meanptop(:)      = R_UNDEF
          if (associated(cospOUT%isccp_meantaucld))    cospOUT%isccp_meantaucld(:)    = R_UNDEF
          if (associated(cospOUT%isccp_meanalbedocld)) cospOUT%isccp_meanalbedocld(:) = R_UNDEF
          if (associated(cospOUT%isccp_boxtau))        cospOUT%isccp_boxtau(:,:)      = R_UNDEF
          if (associated(cospOUT%isccp_boxptop))       cospOUT%isccp_boxptop(:,:)     = R_UNDEF
          if (associated(cospOUT%isccp_fq))            cospOUT%isccp_fq(:,:,:)        = R_UNDEF
       endif
    endif
    if (any([Lisccp_subcolumn,Lisccp_column,Lmodis_subcolumn,Lmodis_column,Lcalipso_column,Lrttov_column,&
             LgrLidar532_column,Latlid_column])) then
       if (any(cospgridIN%phalf .lt. 0)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospgridIN%phalf contains values out of range'
          Lisccp_subcolumn = .false.
          Lisccp_column    = .false.
          Lmodis_subcolumn = .false.
          Lmodis_column    = .false.
          Lcalipso_column  = .false.
          Lrttov_column    = .false.
          Latlid_column    = .false.
          LgrLidar532_column = .false.
          if (allocated(cospOUT%rttov_outputs)) then
             do i=1,cospOUT % Ninst_rttov ! Iterate over each instrument
                 if (associated(cospOUT%rttov_outputs(i)%channel_indices))      cospOUT%rttov_outputs(i)%channel_indices(:)  = 0
                 if (associated(cospOUT%rttov_outputs(i)%bt_total))             cospOUT%rttov_outputs(i)%bt_total(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_clear))             cospOUT%rttov_outputs(i)%bt_clear(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total))            cospOUT%rttov_outputs(i)%rad_total(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_clear))            cospOUT%rttov_outputs(i)%rad_clear(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_cloudy))           cospOUT%rttov_outputs(i)%rad_cloudy(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_total))           cospOUT%rttov_outputs(i)%refl_total(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_clear))           cospOUT%rttov_outputs(i)%refl_clear(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_total_pc))          cospOUT%rttov_outputs(i)%bt_total_pc(:,:)    = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total_pc))         cospOUT%rttov_outputs(i)%rad_total_pc(:,:)   = R_UNDEF                
             end do
          end if
          if (associated(cospOUT%isccp_totalcldarea))  cospOUT%isccp_totalcldarea(:)  = R_UNDEF
          if (associated(cospOUT%isccp_meantb))        cospOUT%isccp_meantb(:)        = R_UNDEF
          if (associated(cospOUT%isccp_meantbclr))     cospOUT%isccp_meantbclr(:)     = R_UNDEF
          if (associated(cospOUT%isccp_meanptop))      cospOUT%isccp_meanptop(:)      = R_UNDEF
          if (associated(cospOUT%isccp_meantaucld))    cospOUT%isccp_meantaucld(:)    = R_UNDEF
          if (associated(cospOUT%isccp_meanalbedocld)) cospOUT%isccp_meanalbedocld(:) = R_UNDEF
          if (associated(cospOUT%isccp_boxtau))        cospOUT%isccp_boxtau(:,:)      = R_UNDEF
          if (associated(cospOUT%isccp_boxptop))       cospOUT%isccp_boxptop(:,:)     = R_UNDEF
          if (associated(cospOUT%isccp_fq))            cospOUT%isccp_fq(:,:,:)        = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_Total_Mean))                          &
               cospOUT%modis_Cloud_Fraction_Total_Mean(:)                   = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_Water_Mean))                          &
               cospOUT%modis_Cloud_Fraction_Water_Mean(:)                   = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_Ice_Mean))                            &
               cospOUT%modis_Cloud_Fraction_Ice_Mean(:)                     = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_High_Mean))                           &
               cospOUT%modis_Cloud_Fraction_High_Mean(:)                    = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_Mid_Mean))                            &
               cospOUT%modis_Cloud_Fraction_Mid_Mean(:)                     = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_Low_Mean))                            &
               cospOUT%modis_Cloud_Fraction_Low_Mean(:)                     = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Total_Mean))                       &
               cospOUT%modis_Optical_Thickness_Total_Mean(:)                = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Water_Mean))                       &
               cospOUT%modis_Optical_Thickness_Water_Mean(:)                = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Ice_Mean))                         &
               cospOUT%modis_Optical_Thickness_Ice_Mean(:)                  = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Total_LogMean))                    &
               cospOUT%modis_Optical_Thickness_Total_LogMean(:)             = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Water_LogMean))                    &
               cospOUT%modis_Optical_Thickness_Water_LogMean(:)             = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Ice_LogMean))                      &
               cospOUT%modis_Optical_Thickness_Ice_LogMean(:)               = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Particle_Size_Water_Mean))                     &
               cospOUT%modis_Cloud_Particle_Size_Water_Mean(:)              = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Particle_Size_Ice_Mean))                       &
               cospOUT%modis_Cloud_Particle_Size_Ice_Mean(:)                = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Top_Pressure_Total_Mean))                      &
               cospOUT%modis_Cloud_Top_Pressure_Total_Mean(:)               = R_UNDEF
          if (associated(cospOUT%modis_Liquid_Water_Path_Mean))                             &
               cospOUT%modis_Liquid_Water_Path_Mean(:)                      = R_UNDEF
          if (associated(cospOUT%modis_Ice_Water_Path_Mean))                                &
               cospOUT%modis_Ice_Water_Path_Mean(:)                         = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure))            &
               cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure(:,:,:) = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Liq))            &
               cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Liq(:,:,:) = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Ice))            &
               cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Ice(:,:,:) = R_UNDEF
          if (associated(cospOUT%modis_LWP_vs_ReffLIQ))            &
               cospOUT%modis_LWP_vs_ReffLIQ(:,:,:) = R_UNDEF
          if (associated(cospOUT%modis_IWP_vs_ReffICE))            &
               cospOUT%modis_IWP_vs_ReffICE(:,:,:) = R_UNDEF    
          if (associated(cospOUT%modis_Optical_Thickness_vs_ReffICE))                       &
               cospOUT%modis_Optical_Thickness_vs_ReffICE(:,:,:)            = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_vs_ReffLIQ))                       &
               cospOUT%modis_Optical_Thickness_vs_ReffLIQ(:,:,:)            = R_UNDEF
          if (associated(cospOUT%calipso_cfad_sr))       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF
          if (associated(cospOUT%calipso_lidarcld))      cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldphase)) cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
          if (associated(cospOUT%calipso_cldlayer))      cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
          if (associated(cospOUT%calipso_cldlayerphase)) cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldtmp))   cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF
          if (associated(cospOUT%atlid_cfad_sr))         cospOUT%atlid_cfad_sr(:,:,:)         = R_UNDEF
          if (associated(cospOUT%atlid_lidarcld))        cospOUT%atlid_lidarcld(:,:)          = R_UNDEF
          if (associated(cospOUT%atlid_cldlayer))        cospOUT%atlid_cldlayer(:,:)          = R_UNDEF
          if (associated(cospOUT%grLidar532_cfad_sr))   cospOUT%grLidar532_cfad_sr(:,:,:)   = R_UNDEF
          if (associated(cospOUT%grLidar532_lidarcld))  cospOUT%grLidar532_lidarcld(:,:)    = R_UNDEF
          if (associated(cospOUT%grLidar532_cldlayer))  cospOUT%grLidar532_cldlayer(:,:)    = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldtype))  cospOUT%calipso_lidarcldtype(:,:,:)  = R_UNDEF 
          if (associated(cospOUT%calipso_cldtype))       cospOUT%calipso_cldtype(:,:)         = R_UNDEF 
          if (associated(cospOUT%calipso_cldtypetemp))   cospOUT%calipso_cldtypetemp(:,:)     = R_UNDEF
          if (associated(cospOUT%calipso_cldtypemeanz))  cospOUT%calipso_cldtypemeanz(:,:)    = R_UNDEF
          if (associated(cospOUT%calipso_cldtypemeanzse)) cospOUT%calipso_cldtypemeanzse(:,:) = R_UNDEF
          if (associated(cospOUT%calipso_cldthinemis))   cospOUT%calipso_cldthinemis(:)       = R_UNDEF
       endif
    endif
    if (any([Lisccp_subcolumn,Lisccp_column,Lrttov_column])) then
       if (any(cospgridIN%qv .lt. 0)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospgridIN%qv contains values out of range'
          Lisccp_subcolumn = .false.
          Lisccp_column    = .false.
          Lrttov_column    = .false.
          if (allocated(cospOUT%rttov_outputs)) then
             do i=1,cospOUT % Ninst_rttov ! Iterate over each instrument
                 if (associated(cospOUT%rttov_outputs(i)%channel_indices))      cospOUT%rttov_outputs(i)%channel_indices(:)  = 0
                 if (associated(cospOUT%rttov_outputs(i)%bt_total))             cospOUT%rttov_outputs(i)%bt_total(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_clear))             cospOUT%rttov_outputs(i)%bt_clear(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total))            cospOUT%rttov_outputs(i)%rad_total(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_clear))            cospOUT%rttov_outputs(i)%rad_clear(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_cloudy))           cospOUT%rttov_outputs(i)%rad_cloudy(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_total))           cospOUT%rttov_outputs(i)%refl_total(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_clear))           cospOUT%rttov_outputs(i)%refl_clear(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_total_pc))          cospOUT%rttov_outputs(i)%bt_total_pc(:,:)    = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total_pc))         cospOUT%rttov_outputs(i)%rad_total_pc(:,:)   = R_UNDEF                
             end do
          end if         
          if (associated(cospOUT%isccp_totalcldarea))  cospOUT%isccp_totalcldarea(:)  = R_UNDEF
          if (associated(cospOUT%isccp_meantb))        cospOUT%isccp_meantb(:)        = R_UNDEF
          if (associated(cospOUT%isccp_meantbclr))     cospOUT%isccp_meantbclr(:)     = R_UNDEF
          if (associated(cospOUT%isccp_meanptop))      cospOUT%isccp_meanptop(:)      = R_UNDEF
          if (associated(cospOUT%isccp_meantaucld))    cospOUT%isccp_meantaucld(:)    = R_UNDEF
          if (associated(cospOUT%isccp_meanalbedocld)) cospOUT%isccp_meanalbedocld(:) = R_UNDEF
          if (associated(cospOUT%isccp_boxtau))        cospOUT%isccp_boxtau(:,:)      = R_UNDEF
          if (associated(cospOUT%isccp_boxptop))       cospOUT%isccp_boxptop(:,:)     = R_UNDEF
          if (associated(cospOUT%isccp_fq))            cospOUT%isccp_fq(:,:,:)        = R_UNDEF
       endif
    endif
    if (any([Lmisr_subcolumn,Lmisr_column,Lcloudsat_subcolumn,Lcloudsat_column,Lcalipso_column,Lradar_lidar_tcc,&
         Llidar_only_freq_cloud,LgrLidar532_column,Latlid_column,Lcloudsat_tcc, Lcloudsat_tcc2, &
         Lcloudsat_modis_wr])) then
       if (any(cospgridIN%hgt_matrix .lt. -300)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospgridIN%hgt_matrix contains values out of range'
          Lmisr_subcolumn     = .false.
          Lmisr_column        = .false.
          Lcloudsat_subcolumn = .false.
          Lcloudsat_column    = .false.
          Lcalipso_column     = .false.
          Lradar_lidar_tcc    = .false.
          Llidar_only_freq_cloud = .false.
          Lcloudsat_tcc       = .false.
          Lcloudsat_tcc2      = .false.
          Latlid_column       = .false.
          LgrLidar532_column  = .false.
          Lcloudsat_modis_wr  = .false.
          if (associated(cospOUT%misr_fq))                   cospOUT%misr_fq(:,:,:)                 = R_UNDEF
          if (associated(cospOUT%misr_dist_model_layertops)) cospOUT%misr_dist_model_layertops(:,:) = R_UNDEF
          if (associated(cospOUT%misr_meanztop))             cospOUT%misr_meanztop(:)               = R_UNDEF
          if (associated(cospOUT%misr_cldarea))              cospOUT%misr_cldarea(:)                = R_UNDEF
          if (associated(cospOUT%calipso_cfad_sr))           cospOUT%calipso_cfad_sr(:,:,:)         = R_UNDEF
          if (associated(cospOUT%calipso_lidarcld))          cospOUT%calipso_lidarcld(:,:)          = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldphase))     cospOUT%calipso_lidarcldphase(:,:,:)   = R_UNDEF
          if (associated(cospOUT%calipso_cldlayer))          cospOUT%calipso_cldlayer(:,:)          = R_UNDEF
          if (associated(cospOUT%calipso_cldlayerphase))     cospOUT%calipso_cldlayerphase(:,:,:)   = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldtmp))       cospOUT%calipso_lidarcldtmp(:,:,:)     = R_UNDEF
          if (associated(cospOUT%cloudsat_cfad_ze))          cospOUT%cloudsat_cfad_ze(:,:,:)        = R_UNDEF
          if (associated(cospOUT%cloudsat_Ze_tot))           cospOUT%cloudsat_Ze_tot(:,:,:)         = R_UNDEF
          if (associated(cospOUT%lidar_only_freq_cloud))     cospOUT%lidar_only_freq_cloud(:,:)     = R_UNDEF
          if (associated(cospOUT%radar_lidar_tcc))           cospOUT%radar_lidar_tcc(:)             = R_UNDEF
          if (associated(cospOUT%cloudsat_tcc))              cospOUT%cloudsat_tcc(:)                = R_UNDEF
          if (associated(cospOUT%cloudsat_tcc2))             cospOUT%cloudsat_tcc2(:)               = R_UNDEF
          if (associated(cospOUT%atlid_cfad_sr))             cospOUT%atlid_cfad_sr(:,:,:)           = R_UNDEF
          if (associated(cospOUT%atlid_lidarcld))            cospOUT%atlid_lidarcld(:,:)            = R_UNDEF
          if (associated(cospOUT%atlid_cldlayer))            cospOUT%atlid_cldlayer(:,:)            = R_UNDEF
          if (associated(cospOUT%grLidar532_cfad_sr))        cospOUT%grLidar532_cfad_sr(:,:,:)      = R_UNDEF
          if (associated(cospOUT%grLidar532_lidarcld))       cospOUT%grLidar532_lidarcld(:,:)       = R_UNDEF
          if (associated(cospOUT%grLidar532_cldlayer))       cospOUT%grLidar532_cldlayer(:,:)       = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldtype))      cospOUT%calipso_lidarcldtype(:,:,:)    = R_UNDEF
          if (associated(cospOUT%calipso_cldtype))           cospOUT%calipso_cldtype(:,:)           = R_UNDEF 
          if (associated(cospOUT%calipso_cldtypetemp))       cospOUT%calipso_cldtypetemp(:,:)       = R_UNDEF 
          if (associated(cospOUT%calipso_cldtypemeanz))      cospOUT%calipso_cldtypemeanz(:,:)      = R_UNDEF 
          if (associated(cospOUT%calipso_cldtypemeanzse))    cospOUT%calipso_cldtypemeanzse(:,:)    = R_UNDEF
          if (associated(cospOUT%calipso_cldthinemis))       cospOUT%calipso_cldthinemis(:)         = R_UNDEF
          if (associated(cospOUT%cfodd_ntotal))              cospOUT%cfodd_ntotal(:,:,:,:)          = R_UNDEF
          if (associated(cospOUT%wr_occfreq_ntotal))         cospOUT%wr_occfreq_ntotal(:,:)         = R_UNDEF
       endif
    endif
    if (any([Lrttov_column,Lcloudsat_column,Lcalipso_column,Lradar_lidar_tcc,Llidar_only_freq_cloud, &
             LgrLidar532_column, Latlid_column, Lcloudsat_tcc, Lcloudsat_tcc2, Lcloudsat_modis_wr])) then
       if (any(cospgridIN%hgt_matrix_half .lt. -300)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospgridIN%hgt_matrix_half contains values out of range'
          Lrttov_column    = .false.
          Lcloudsat_column = .false.
          Lcalipso_column  = .false.
          Lradar_lidar_tcc = .false.
          Llidar_only_freq_cloud = .false.
          Lcloudsat_tcc    = .false.
          Lcloudsat_tcc2   = .false.
          Latlid_column       = .false.
          LgrLidar532_column = .false.
          Lcloudsat_modis_wr = .false.
          if (allocated(cospOUT%rttov_outputs)) then
             do i=1,cospOUT % Ninst_rttov ! Iterate over each instrument
                 if (associated(cospOUT%rttov_outputs(i)%channel_indices))      cospOUT%rttov_outputs(i)%channel_indices(:)  = 0
                 if (associated(cospOUT%rttov_outputs(i)%bt_total))             cospOUT%rttov_outputs(i)%bt_total(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_clear))             cospOUT%rttov_outputs(i)%bt_clear(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total))            cospOUT%rttov_outputs(i)%rad_total(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_clear))            cospOUT%rttov_outputs(i)%rad_clear(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_cloudy))           cospOUT%rttov_outputs(i)%rad_cloudy(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_total))           cospOUT%rttov_outputs(i)%refl_total(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_clear))           cospOUT%rttov_outputs(i)%refl_clear(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_total_pc))          cospOUT%rttov_outputs(i)%bt_total_pc(:,:)    = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total_pc))         cospOUT%rttov_outputs(i)%rad_total_pc(:,:)   = R_UNDEF                
             end do
          end if
          if (associated(cospOUT%calipso_cfad_sr))        cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF
          if (associated(cospOUT%calipso_lidarcld))       cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldphase))  cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
          if (associated(cospOUT%calipso_cldlayer))       cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
          if (associated(cospOUT%calipso_cldlayerphase))  cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldtmp))    cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF
          if (associated(cospOUT%cloudsat_cfad_ze))       cospOUT%cloudsat_cfad_ze(:,:,:)      = R_UNDEF
          if (associated(cospOUT%lidar_only_freq_cloud))  cospOUT%lidar_only_freq_cloud(:,:)   = R_UNDEF
          if (associated(cospOUT%radar_lidar_tcc))        cospOUT%radar_lidar_tcc(:)           = R_UNDEF
          if (associated(cospOUT%cloudsat_tcc))           cospOUT%cloudsat_tcc(:)              = R_UNDEF
          if (associated(cospOUT%cloudsat_tcc2))          cospOUT%cloudsat_tcc2(:)             = R_UNDEF          
          if (associated(cospOUT%atlid_cfad_sr))          cospOUT%atlid_cfad_sr(:,:,:)         = R_UNDEF
          if (associated(cospOUT%atlid_lidarcld))         cospOUT%atlid_lidarcld(:,:)          = R_UNDEF
          if (associated(cospOUT%atlid_cldlayer))         cospOUT%atlid_cldlayer(:,:)          = R_UNDEF
          if (associated(cospOUT%grLidar532_cfad_sr))     cospOUT%grLidar532_cfad_sr(:,:,:)    = R_UNDEF
          if (associated(cospOUT%grLidar532_lidarcld))    cospOUT%grLidar532_lidarcld(:,:)     = R_UNDEF
          if (associated(cospOUT%grLidar532_cldlayer))    cospOUT%grLidar532_cldlayer(:,:)     = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldtype))   cospOUT%calipso_lidarcldtype(:,:,:)  = R_UNDEF
          if (associated(cospOUT%calipso_cldtype))        cospOUT%calipso_cldtype(:,:)         = R_UNDEF 
          if (associated(cospOUT%calipso_cldtypetemp))    cospOUT%calipso_cldtypetemp(:,:)     = R_UNDEF 
          if (associated(cospOUT%calipso_cldtypemeanz))   cospOUT%calipso_cldtypemeanz(:,:)    = R_UNDEF 
          if (associated(cospOUT%calipso_cldtypemeanzse)) cospOUT%calipso_cldtypemeanzse(:,:)  = R_UNDEF 
          if (associated(cospOUT%calipso_cldthinemis))    cospOUT%calipso_cldthinemis(:)       = R_UNDEF
          if (associated(cospOUT%cfodd_ntotal))           cospOUT%cfodd_ntotal(:,:,:,:)        = R_UNDEF
          if (associated(cospOUT%wr_occfreq_ntotal))      cospOUT%wr_occfreq_ntotal(:,:)       = R_UNDEF
       endif
    endif
    if (any([Lrttov_column,Lcalipso_column,Lparasol_column])) then
       if (any(cospgridIN%land .lt. 0)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospgridIN%land contains values out of range'
          Lrttov_column    = .false.
          Lcalipso_column  = .false.
          Lparasol_column  = .false.
          if (allocated(cospOUT%rttov_outputs)) then
             do i=1,cospOUT % Ninst_rttov ! Iterate over each instrument
                 if (associated(cospOUT%rttov_outputs(i)%channel_indices))      cospOUT%rttov_outputs(i)%channel_indices(:)  = 0
                 if (associated(cospOUT%rttov_outputs(i)%bt_total))             cospOUT%rttov_outputs(i)%bt_total(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_clear))             cospOUT%rttov_outputs(i)%bt_clear(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total))            cospOUT%rttov_outputs(i)%rad_total(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_clear))            cospOUT%rttov_outputs(i)%rad_clear(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_cloudy))           cospOUT%rttov_outputs(i)%rad_cloudy(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_total))           cospOUT%rttov_outputs(i)%refl_total(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_clear))           cospOUT%rttov_outputs(i)%refl_clear(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_total_pc))          cospOUT%rttov_outputs(i)%bt_total_pc(:,:)    = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total_pc))         cospOUT%rttov_outputs(i)%rad_total_pc(:,:)   = R_UNDEF                
             end do
          end if         
          if (associated(cospOUT%calipso_cfad_sr))       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF
          if (associated(cospOUT%calipso_lidarcld))      cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldphase)) cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
          if (associated(cospOUT%calipso_cldlayer))      cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
          if (associated(cospOUT%calipso_cldlayerphase)) cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldtmp))   cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldtype))  cospOUT%calipso_lidarcldtype(:,:,:)  = R_UNDEF
          if (associated(cospOUT%calipso_cldtype))       cospOUT%calipso_cldtype(:,:)         = R_UNDEF 
          if (associated(cospOUT%calipso_cldtypetemp))   cospOUT%calipso_cldtypetemp(:,:)     = R_UNDEF 
          if (associated(cospOUT%calipso_cldtypemeanz))  cospOUT%calipso_cldtypemeanz(:,:)    = R_UNDEF
          if (associated(cospOUT%calipso_cldtypemeanzse)) cospOUT%calipso_cldtypemeanzse(:,:) = R_UNDEF
          if (associated(cospOUT%calipso_cldthinemis))   cospOUT%calipso_cldthinemis(:)       = R_UNDEF
          if (associated(cospOUT%parasolGrid_refl))      cospOUT%parasolGrid_refl(:,:)        = R_UNDEF
       endif
    endif
    if (any([Lisccp_subcolumn,Lisccp_column,Lrttov_column])) then
       if (any(cospgridIN%skt .lt. 0)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospgridIN%skt contains values out of range'
          Lisccp_subcolumn = .false.
          Lisccp_column    = .false.
          Lrttov_column    = .false.
          if (allocated(cospOUT%rttov_outputs)) then
             do i=1,cospOUT % Ninst_rttov ! Iterate over each instrument
                 if (associated(cospOUT%rttov_outputs(i)%channel_indices))      cospOUT%rttov_outputs(i)%channel_indices(:)  = 0
                 if (associated(cospOUT%rttov_outputs(i)%bt_total))             cospOUT%rttov_outputs(i)%bt_total(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_clear))             cospOUT%rttov_outputs(i)%bt_clear(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total))            cospOUT%rttov_outputs(i)%rad_total(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_clear))            cospOUT%rttov_outputs(i)%rad_clear(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_cloudy))           cospOUT%rttov_outputs(i)%rad_cloudy(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_total))           cospOUT%rttov_outputs(i)%refl_total(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_clear))           cospOUT%rttov_outputs(i)%refl_clear(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_total_pc))          cospOUT%rttov_outputs(i)%bt_total_pc(:,:)    = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total_pc))         cospOUT%rttov_outputs(i)%rad_total_pc(:,:)   = R_UNDEF                
             end do
          end if
          if (associated(cospOUT%isccp_totalcldarea))  cospOUT%isccp_totalcldarea(:)  = R_UNDEF
          if (associated(cospOUT%isccp_meantb))        cospOUT%isccp_meantb(:)        = R_UNDEF
          if (associated(cospOUT%isccp_meantbclr))     cospOUT%isccp_meantbclr(:)     = R_UNDEF
          if (associated(cospOUT%isccp_meanptop))      cospOUT%isccp_meanptop(:)      = R_UNDEF
          if (associated(cospOUT%isccp_meantaucld))    cospOUT%isccp_meantaucld(:)    = R_UNDEF
          if (associated(cospOUT%isccp_meanalbedocld)) cospOUT%isccp_meanalbedocld(:) = R_UNDEF
          if (associated(cospOUT%isccp_boxtau))        cospOUT%isccp_boxtau(:,:)      = R_UNDEF
          if (associated(cospOUT%isccp_boxptop))       cospOUT%isccp_boxptop(:,:)     = R_UNDEF
          if (associated(cospOUT%isccp_fq))            cospOUT%isccp_fq(:,:,:)        = R_UNDEF
       endif
    endif
    
    ! RTTOV Inputs
    if (Lrttov_column) then
       if (any(cospgridIN%co2 .lt. 0)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospgridIN%co2 contains values out of range'
          Lrttov_column = .false.
          if (allocated(cospOUT%rttov_outputs)) then
             do i=1,cospOUT % Ninst_rttov ! Iterate over each instrument
                 if (associated(cospOUT%rttov_outputs(i)%channel_indices))      cospOUT%rttov_outputs(i)%channel_indices(:)  = 0
                 if (associated(cospOUT%rttov_outputs(i)%bt_total))             cospOUT%rttov_outputs(i)%bt_total(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_clear))             cospOUT%rttov_outputs(i)%bt_clear(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total))            cospOUT%rttov_outputs(i)%rad_total(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_clear))            cospOUT%rttov_outputs(i)%rad_clear(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_cloudy))           cospOUT%rttov_outputs(i)%rad_cloudy(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_total))           cospOUT%rttov_outputs(i)%refl_total(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_clear))           cospOUT%rttov_outputs(i)%refl_clear(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_total_pc))          cospOUT%rttov_outputs(i)%bt_total_pc(:,:)    = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total_pc))         cospOUT%rttov_outputs(i)%rad_total_pc(:,:)   = R_UNDEF                
             end do
          end if
       endif
       if (any(cospgridIN%ch4 .lt. 0)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospgridIN%ch4 contains values out of range'
          Lrttov_column = .false.
          if (allocated(cospOUT%rttov_outputs)) then
             do i=1,cospOUT % Ninst_rttov ! Iterate over each instrument
                 if (associated(cospOUT%rttov_outputs(i)%channel_indices))      cospOUT%rttov_outputs(i)%channel_indices(:)  = 0
                 if (associated(cospOUT%rttov_outputs(i)%bt_total))             cospOUT%rttov_outputs(i)%bt_total(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_clear))             cospOUT%rttov_outputs(i)%bt_clear(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total))            cospOUT%rttov_outputs(i)%rad_total(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_clear))            cospOUT%rttov_outputs(i)%rad_clear(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_cloudy))           cospOUT%rttov_outputs(i)%rad_cloudy(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_total))           cospOUT%rttov_outputs(i)%refl_total(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_clear))           cospOUT%rttov_outputs(i)%refl_clear(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_total_pc))          cospOUT%rttov_outputs(i)%bt_total_pc(:,:)    = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total_pc))         cospOUT%rttov_outputs(i)%rad_total_pc(:,:)   = R_UNDEF                
             end do
          end if
       endif
       if (any(cospgridIN%n2o .lt. 0)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospgridIN%n2o contains values out of range'
          Lrttov_column = .false.
          if (allocated(cospOUT%rttov_outputs)) then
             do i=1,cospOUT % Ninst_rttov ! Iterate over each instrument
                 if (associated(cospOUT%rttov_outputs(i)%channel_indices))      cospOUT%rttov_outputs(i)%channel_indices(:)  = 0
                 if (associated(cospOUT%rttov_outputs(i)%bt_total))             cospOUT%rttov_outputs(i)%bt_total(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_clear))             cospOUT%rttov_outputs(i)%bt_clear(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total))            cospOUT%rttov_outputs(i)%rad_total(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_clear))            cospOUT%rttov_outputs(i)%rad_clear(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_cloudy))           cospOUT%rttov_outputs(i)%rad_cloudy(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_total))           cospOUT%rttov_outputs(i)%refl_total(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_clear))           cospOUT%rttov_outputs(i)%refl_clear(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_total_pc))          cospOUT%rttov_outputs(i)%bt_total_pc(:,:)    = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total_pc))         cospOUT%rttov_outputs(i)%rad_total_pc(:,:)   = R_UNDEF                
             end do
          end if
       endif
       if (any(cospgridIN%co.lt. 0)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospgridIN%co contains values out of range'
          Lrttov_column = .false.
          if (allocated(cospOUT%rttov_outputs)) then
             do i=1,cospOUT % Ninst_rttov ! Iterate over each instrument
                 if (associated(cospOUT%rttov_outputs(i)%channel_indices))      cospOUT%rttov_outputs(i)%channel_indices(:)  = 0
                 if (associated(cospOUT%rttov_outputs(i)%bt_total))             cospOUT%rttov_outputs(i)%bt_total(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_clear))             cospOUT%rttov_outputs(i)%bt_clear(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total))            cospOUT%rttov_outputs(i)%rad_total(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_clear))            cospOUT%rttov_outputs(i)%rad_clear(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_cloudy))           cospOUT%rttov_outputs(i)%rad_cloudy(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_total))           cospOUT%rttov_outputs(i)%refl_total(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_clear))           cospOUT%rttov_outputs(i)%refl_clear(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_total_pc))          cospOUT%rttov_outputs(i)%bt_total_pc(:,:)    = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total_pc))         cospOUT%rttov_outputs(i)%rad_total_pc(:,:)   = R_UNDEF                
             end do
          end if
       endif
       if (any(cospgridIN%o3 .lt. 0)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospgridIN%o3 contains values out of range'
          Lrttov_column = .false.
          if (allocated(cospOUT%rttov_outputs)) then
             do i=1,cospOUT % Ninst_rttov ! Iterate over each instrument
                 if (associated(cospOUT%rttov_outputs(i)%channel_indices))      cospOUT%rttov_outputs(i)%channel_indices(:)  = 0
                 if (associated(cospOUT%rttov_outputs(i)%bt_total))             cospOUT%rttov_outputs(i)%bt_total(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_clear))             cospOUT%rttov_outputs(i)%bt_clear(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total))            cospOUT%rttov_outputs(i)%rad_total(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_clear))            cospOUT%rttov_outputs(i)%rad_clear(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_cloudy))           cospOUT%rttov_outputs(i)%rad_cloudy(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_total))           cospOUT%rttov_outputs(i)%refl_total(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_clear))           cospOUT%rttov_outputs(i)%refl_clear(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_total_pc))          cospOUT%rttov_outputs(i)%bt_total_pc(:,:)    = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total_pc))         cospOUT%rttov_outputs(i)%rad_total_pc(:,:)   = R_UNDEF                
             end do
          end if
       endif
!       if (any(cospgridIN%emis_in .lt. 0. .OR. cospgridIN%emis_in .gt. 1)) then
!          nError=nError+1
!          errorMessage(nError) = 'ERROR: COSP input variable: cospgridIN%emis_in contains values out of range'
!          Lrttov_column = .false.
!          if (allocated(cospOUT%rttov_outputs)) then
!             do i=1,cospOUT % Ninst_rttov ! Iterate over each instrument
!                 if (associated(cospOUT%rttov_outputs(i)%channel_indices))      cospOUT%rttov_outputs(i)%channel_indices(:)  = 0
!                 if (associated(cospOUT%rttov_outputs(i)%bt_total))             cospOUT%rttov_outputs(i)%bt_total(:,:)       = R_UNDEF
!                 if (associated(cospOUT%rttov_outputs(i)%bt_clear))             cospOUT%rttov_outputs(i)%bt_clear(:,:)       = R_UNDEF
!                 if (associated(cospOUT%rttov_outputs(i)%rad_total))            cospOUT%rttov_outputs(i)%rad_total(:,:)      = R_UNDEF
!                 if (associated(cospOUT%rttov_outputs(i)%rad_clear))            cospOUT%rttov_outputs(i)%rad_clear(:,:)      = R_UNDEF
!                 if (associated(cospOUT%rttov_outputs(i)%rad_cloudy))           cospOUT%rttov_outputs(i)%rad_cloudy(:,:)     = R_UNDEF
!                 if (associated(cospOUT%rttov_outputs(i)%refl_total))           cospOUT%rttov_outputs(i)%refl_total(:,:)     = R_UNDEF
!                 if (associated(cospOUT%rttov_outputs(i)%refl_clear))           cospOUT%rttov_outputs(i)%refl_clear(:,:)     = R_UNDEF
!                 if (associated(cospOUT%rttov_outputs(i)%bt_total_pc))          cospOUT%rttov_outputs(i)%bt_total_pc(:,:)    = R_UNDEF
!                 if (associated(cospOUT%rttov_outputs(i)%rad_total_pc))         cospOUT%rttov_outputs(i)%rad_total_pc(:,:)   = R_UNDEF                
!             end do
!          end if
!       endif
!       if (any(cospgridIN%refl_in .lt. 0. .OR. cospgridIN%refl_in .gt. 1)) then
!          nError=nError+1
!          errorMessage(nError) = 'ERROR: COSP input variable: cospgridIN%refl_in contains values out of range'
!          Lrttov_column = .false.
!          if (allocated(cospOUT%rttov_outputs)) then
!             do i=1,cospOUT % Ninst_rttov ! Iterate over each instrument
!                 if (associated(cospOUT%rttov_outputs(i)%channel_indices))      cospOUT%rttov_outputs(i)%channel_indices(:)  = 0
!                 if (associated(cospOUT%rttov_outputs(i)%bt_total))             cospOUT%rttov_outputs(i)%bt_total(:,:)       = R_UNDEF
!                 if (associated(cospOUT%rttov_outputs(i)%bt_clear))             cospOUT%rttov_outputs(i)%bt_clear(:,:)       = R_UNDEF
!                 if (associated(cospOUT%rttov_outputs(i)%rad_total))            cospOUT%rttov_outputs(i)%rad_total(:,:)      = R_UNDEF
!                 if (associated(cospOUT%rttov_outputs(i)%rad_clear))            cospOUT%rttov_outputs(i)%rad_clear(:,:)      = R_UNDEF
!                 if (associated(cospOUT%rttov_outputs(i)%rad_cloudy))           cospOUT%rttov_outputs(i)%rad_cloudy(:,:)     = R_UNDEF
!                 if (associated(cospOUT%rttov_outputs(i)%refl_total))           cospOUT%rttov_outputs(i)%refl_total(:,:)     = R_UNDEF
!                 if (associated(cospOUT%rttov_outputs(i)%refl_clear))           cospOUT%rttov_outputs(i)%refl_clear(:,:)     = R_UNDEF
!                 if (associated(cospOUT%rttov_outputs(i)%bt_total_pc))          cospOUT%rttov_outputs(i)%bt_total_pc(:,:)    = R_UNDEF
!                 if (associated(cospOUT%rttov_outputs(i)%rad_total_pc))         cospOUT%rttov_outputs(i)%rad_total_pc(:,:)   = R_UNDEF                
!             end do
!          end if
!       endif
       if (any(cospgridIN%rttov_sfcmask .lt. 0 .or. cospgridIN%rttov_sfcmask .gt. 2)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospgridIN%rttov_sfcmask contains values out of range'
          Lrttov_column = .false.
          if (allocated(cospOUT%rttov_outputs)) then
             do i=1,cospOUT % Ninst_rttov ! Iterate over each instrument
                 if (associated(cospOUT%rttov_outputs(i)%channel_indices))      cospOUT%rttov_outputs(i)%channel_indices(:)  = 0
                 if (associated(cospOUT%rttov_outputs(i)%bt_total))             cospOUT%rttov_outputs(i)%bt_total(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_clear))             cospOUT%rttov_outputs(i)%bt_clear(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total))            cospOUT%rttov_outputs(i)%rad_total(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_clear))            cospOUT%rttov_outputs(i)%rad_clear(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_cloudy))           cospOUT%rttov_outputs(i)%rad_cloudy(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_total))           cospOUT%rttov_outputs(i)%refl_total(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_clear))           cospOUT%rttov_outputs(i)%refl_clear(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_total_pc))          cospOUT%rttov_outputs(i)%bt_total_pc(:,:)    = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total_pc))         cospOUT%rttov_outputs(i)%rad_total_pc(:,:)   = R_UNDEF                
             end do
          end if
       endif
       if (any(cospgridIN%u_sfc .lt. -100. .OR. cospgridIN%u_sfc .gt. 100.)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospIN%u_sfc contains values out of range'
          Lrttov_column = .false.
          if (allocated(cospOUT%rttov_outputs)) then
             do i=1,cospOUT % Ninst_rttov ! Iterate over each instrument
                 if (associated(cospOUT%rttov_outputs(i)%channel_indices))      cospOUT%rttov_outputs(i)%channel_indices(:)  = 0
                 if (associated(cospOUT%rttov_outputs(i)%bt_total))             cospOUT%rttov_outputs(i)%bt_total(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_clear))             cospOUT%rttov_outputs(i)%bt_clear(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total))            cospOUT%rttov_outputs(i)%rad_total(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_clear))            cospOUT%rttov_outputs(i)%rad_clear(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_cloudy))           cospOUT%rttov_outputs(i)%rad_cloudy(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_total))           cospOUT%rttov_outputs(i)%refl_total(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_clear))           cospOUT%rttov_outputs(i)%refl_clear(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_total_pc))          cospOUT%rttov_outputs(i)%bt_total_pc(:,:)    = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total_pc))         cospOUT%rttov_outputs(i)%rad_total_pc(:,:)   = R_UNDEF                
             end do
          end if
       endif
       if (any(cospgridIN%v_sfc .lt. -100. .OR. cospgridIN%v_sfc .gt. 100.)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospIN%v_sfc contains values out of range'
          Lrttov_column = .false.
          if (allocated(cospOUT%rttov_outputs)) then
             do i=1,cospOUT % Ninst_rttov ! Iterate over each instrument
                 if (associated(cospOUT%rttov_outputs(i)%channel_indices))      cospOUT%rttov_outputs(i)%channel_indices(:)  = 0
                 if (associated(cospOUT%rttov_outputs(i)%bt_total))             cospOUT%rttov_outputs(i)%bt_total(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_clear))             cospOUT%rttov_outputs(i)%bt_clear(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total))            cospOUT%rttov_outputs(i)%rad_total(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_clear))            cospOUT%rttov_outputs(i)%rad_clear(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_cloudy))           cospOUT%rttov_outputs(i)%rad_cloudy(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_total))           cospOUT%rttov_outputs(i)%refl_total(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_clear))           cospOUT%rttov_outputs(i)%refl_clear(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_total_pc))          cospOUT%rttov_outputs(i)%bt_total_pc(:,:)    = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total_pc))         cospOUT%rttov_outputs(i)%rad_total_pc(:,:)   = R_UNDEF                
             end do
          end if
       endif
       if (any(cospgridIN%lat .lt. -90 .OR. cospgridIN%lat .gt. 90)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospIN%lat contains values out of range'
          Lrttov_column = .false.
          if (allocated(cospOUT%rttov_outputs)) then
             do i=1,cospOUT % Ninst_rttov ! Iterate over each instrument
                 if (associated(cospOUT%rttov_outputs(i)%channel_indices))      cospOUT%rttov_outputs(i)%channel_indices(:)  = 0
                 if (associated(cospOUT%rttov_outputs(i)%bt_total))             cospOUT%rttov_outputs(i)%bt_total(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_clear))             cospOUT%rttov_outputs(i)%bt_clear(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total))            cospOUT%rttov_outputs(i)%rad_total(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_clear))            cospOUT%rttov_outputs(i)%rad_clear(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_cloudy))           cospOUT%rttov_outputs(i)%rad_cloudy(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_total))           cospOUT%rttov_outputs(i)%refl_total(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_clear))           cospOUT%rttov_outputs(i)%refl_clear(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_total_pc))          cospOUT%rttov_outputs(i)%bt_total_pc(:,:)    = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total_pc))         cospOUT%rttov_outputs(i)%rad_total_pc(:,:)   = R_UNDEF                
             end do
          end if
       endif
       if (any(cospgridIN%lon .lt. -180 .OR. cospgridIN%lon .gt. 360)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospIN%lon contains values out of range'
          Lrttov_column = .false.
          if (allocated(cospOUT%rttov_outputs)) then
             do i=1,cospOUT % Ninst_rttov ! Iterate over each instrument
                 if (associated(cospOUT%rttov_outputs(i)%channel_indices))      cospOUT%rttov_outputs(i)%channel_indices(:)  = 0
                 if (associated(cospOUT%rttov_outputs(i)%bt_total))             cospOUT%rttov_outputs(i)%bt_total(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_clear))             cospOUT%rttov_outputs(i)%bt_clear(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total))            cospOUT%rttov_outputs(i)%rad_total(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_clear))            cospOUT%rttov_outputs(i)%rad_clear(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_cloudy))           cospOUT%rttov_outputs(i)%rad_cloudy(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_total))           cospOUT%rttov_outputs(i)%refl_total(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_clear))           cospOUT%rttov_outputs(i)%refl_clear(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_total_pc))          cospOUT%rttov_outputs(i)%bt_total_pc(:,:)    = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total_pc))         cospOUT%rttov_outputs(i)%rad_total_pc(:,:)   = R_UNDEF                
             end do
          end if
       endif       
       if (any(cospgridIN%tca .lt. 0 .OR. cospgridIN%tca .gt. 1)) then ! tca on [0,1]
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospIN%tca contains values out of range'
          Lrttov_column = .false.
          if (allocated(cospOUT%rttov_outputs)) then
             do i=1,cospOUT % Ninst_rttov ! Iterate over each instrument
                 if (associated(cospOUT%rttov_outputs(i)%channel_indices))      cospOUT%rttov_outputs(i)%channel_indices(:)  = 0
                 if (associated(cospOUT%rttov_outputs(i)%bt_total))             cospOUT%rttov_outputs(i)%bt_total(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_clear))             cospOUT%rttov_outputs(i)%bt_clear(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total))            cospOUT%rttov_outputs(i)%rad_total(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_clear))            cospOUT%rttov_outputs(i)%rad_clear(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_cloudy))           cospOUT%rttov_outputs(i)%rad_cloudy(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_total))           cospOUT%rttov_outputs(i)%refl_total(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_clear))           cospOUT%rttov_outputs(i)%refl_clear(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_total_pc))          cospOUT%rttov_outputs(i)%bt_total_pc(:,:)    = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total_pc))         cospOUT%rttov_outputs(i)%rad_total_pc(:,:)   = R_UNDEF                
             end do
          end if
       endif 
       if (any(cospgridIN%rttov_time(:,1) .lt. 0 .OR. cospgridIN%rttov_time(:,1) .gt. 24)) then ! rttov_time(1), hour on [0,24]
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospIN%rttov_time(1) contains values out of range'
          Lrttov_column = .false.
          if (allocated(cospOUT%rttov_outputs)) then
             do i=1,cospOUT % Ninst_rttov ! Iterate over each instrument
                 if (associated(cospOUT%rttov_outputs(i)%channel_indices))      cospOUT%rttov_outputs(i)%channel_indices(:)  = 0
                 if (associated(cospOUT%rttov_outputs(i)%bt_total))             cospOUT%rttov_outputs(i)%bt_total(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_clear))             cospOUT%rttov_outputs(i)%bt_clear(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total))            cospOUT%rttov_outputs(i)%rad_total(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_clear))            cospOUT%rttov_outputs(i)%rad_clear(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_cloudy))           cospOUT%rttov_outputs(i)%rad_cloudy(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_total))           cospOUT%rttov_outputs(i)%refl_total(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_clear))           cospOUT%rttov_outputs(i)%refl_clear(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_total_pc))          cospOUT%rttov_outputs(i)%bt_total_pc(:,:)    = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total_pc))         cospOUT%rttov_outputs(i)%rad_total_pc(:,:)   = R_UNDEF                
             end do
          end if
       endif
       if (any(cospgridIN%rttov_time(:,2) .lt. 0 .OR. cospgridIN%rttov_time(:,2) .gt. 60)) then ! rttov_time(2), minute on [0,60]
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospIN%rttov_time(2) contains values out of range'
          Lrttov_column = .false.
          if (allocated(cospOUT%rttov_outputs)) then
             do i=1,cospOUT % Ninst_rttov ! Iterate over each instrument
                 if (associated(cospOUT%rttov_outputs(i)%channel_indices))      cospOUT%rttov_outputs(i)%channel_indices(:)  = 0
                 if (associated(cospOUT%rttov_outputs(i)%bt_total))             cospOUT%rttov_outputs(i)%bt_total(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_clear))             cospOUT%rttov_outputs(i)%bt_clear(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total))            cospOUT%rttov_outputs(i)%rad_total(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_clear))            cospOUT%rttov_outputs(i)%rad_clear(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_cloudy))           cospOUT%rttov_outputs(i)%rad_cloudy(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_total))           cospOUT%rttov_outputs(i)%refl_total(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_clear))           cospOUT%rttov_outputs(i)%refl_clear(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_total_pc))          cospOUT%rttov_outputs(i)%bt_total_pc(:,:)    = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total_pc))         cospOUT%rttov_outputs(i)%rad_total_pc(:,:)   = R_UNDEF                
             end do
          end if
       endif 
       if (any(cospgridIN%rttov_time(:,3) .lt. 0 .OR. cospgridIN%rttov_time(:,3) .gt. 60)) then ! rttov_time(3), second on [0,60]
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospIN%rttov_time(3) contains values out of range'
          Lrttov_column = .false.
          if (allocated(cospOUT%rttov_outputs)) then
             do i=1,cospOUT % Ninst_rttov ! Iterate over each instrument
                 if (associated(cospOUT%rttov_outputs(i)%channel_indices))      cospOUT%rttov_outputs(i)%channel_indices(:)  = 0
                 if (associated(cospOUT%rttov_outputs(i)%bt_total))             cospOUT%rttov_outputs(i)%bt_total(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_clear))             cospOUT%rttov_outputs(i)%bt_clear(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total))            cospOUT%rttov_outputs(i)%rad_total(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_clear))            cospOUT%rttov_outputs(i)%rad_clear(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_cloudy))           cospOUT%rttov_outputs(i)%rad_cloudy(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_total))           cospOUT%rttov_outputs(i)%refl_total(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_clear))           cospOUT%rttov_outputs(i)%refl_clear(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_total_pc))          cospOUT%rttov_outputs(i)%bt_total_pc(:,:)    = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total_pc))         cospOUT%rttov_outputs(i)%rad_total_pc(:,:)   = R_UNDEF                
             end do
          end if
       endif   
       if (any(cospgridIN%rttov_date(:,1) .lt. 0)) then ! rttov_date(1), year on [0,inf]
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospIN%rttov_date(1) contains values out of range'
          Lrttov_column = .false.
          if (allocated(cospOUT%rttov_outputs)) then
             do i=1,cospOUT % Ninst_rttov ! Iterate over each instrument
                 if (associated(cospOUT%rttov_outputs(i)%channel_indices))      cospOUT%rttov_outputs(i)%channel_indices(:)  = 0
                 if (associated(cospOUT%rttov_outputs(i)%bt_total))             cospOUT%rttov_outputs(i)%bt_total(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_clear))             cospOUT%rttov_outputs(i)%bt_clear(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total))            cospOUT%rttov_outputs(i)%rad_total(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_clear))            cospOUT%rttov_outputs(i)%rad_clear(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_cloudy))           cospOUT%rttov_outputs(i)%rad_cloudy(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_total))           cospOUT%rttov_outputs(i)%refl_total(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_clear))           cospOUT%rttov_outputs(i)%refl_clear(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_total_pc))          cospOUT%rttov_outputs(i)%bt_total_pc(:,:)    = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total_pc))         cospOUT%rttov_outputs(i)%rad_total_pc(:,:)   = R_UNDEF                
             end do
          end if
       endif
       if (any(cospgridIN%rttov_date(:,2) .lt. 0 .OR. cospgridIN%rttov_date(:,2) .gt. 12)) then ! rttov_date(2), month on [0,12]
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospIN%rttov_date(2) contains values out of range'
          Lrttov_column = .false.
          if (allocated(cospOUT%rttov_outputs)) then
             do i=1,cospOUT % Ninst_rttov ! Iterate over each instrument
                 if (associated(cospOUT%rttov_outputs(i)%channel_indices))      cospOUT%rttov_outputs(i)%channel_indices(:)  = 0
                 if (associated(cospOUT%rttov_outputs(i)%bt_total))             cospOUT%rttov_outputs(i)%bt_total(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_clear))             cospOUT%rttov_outputs(i)%bt_clear(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total))            cospOUT%rttov_outputs(i)%rad_total(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_clear))            cospOUT%rttov_outputs(i)%rad_clear(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_cloudy))           cospOUT%rttov_outputs(i)%rad_cloudy(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_total))           cospOUT%rttov_outputs(i)%refl_total(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_clear))           cospOUT%rttov_outputs(i)%refl_clear(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_total_pc))          cospOUT%rttov_outputs(i)%bt_total_pc(:,:)    = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total_pc))         cospOUT%rttov_outputs(i)%rad_total_pc(:,:)   = R_UNDEF                
             end do
          end if
       endif
       if (any(cospgridIN%rttov_date(:,3) .lt. 1 .OR. cospgridIN%rttov_date(:,3) .gt. 31)) then ! rttov_date(3), day on [1,31]
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospIN%rttov_date(3) contains values out of range'
          Lrttov_column = .false.
          if (allocated(cospOUT%rttov_outputs)) then
             do i=1,cospOUT % Ninst_rttov ! Iterate over each instrument
                 if (associated(cospOUT%rttov_outputs(i)%channel_indices))      cospOUT%rttov_outputs(i)%channel_indices(:)  = 0
                 if (associated(cospOUT%rttov_outputs(i)%bt_total))             cospOUT%rttov_outputs(i)%bt_total(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_clear))             cospOUT%rttov_outputs(i)%bt_clear(:,:)       = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total))            cospOUT%rttov_outputs(i)%rad_total(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_clear))            cospOUT%rttov_outputs(i)%rad_clear(:,:)      = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_cloudy))           cospOUT%rttov_outputs(i)%rad_cloudy(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_total))           cospOUT%rttov_outputs(i)%refl_total(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%refl_clear))           cospOUT%rttov_outputs(i)%refl_clear(:,:)     = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%bt_total_pc))          cospOUT%rttov_outputs(i)%bt_total_pc(:,:)    = R_UNDEF
                 if (associated(cospOUT%rttov_outputs(i)%rad_total_pc))         cospOUT%rttov_outputs(i)%rad_total_pc(:,:)   = R_UNDEF                
             end do
          end if
       endif       
    endif

    ! COSP_INPUTS
    if (any([Lisccp_subcolumn,Lisccp_column])) then
       if (cospIN%emsfc_lw .lt. 0. .OR. cospIN%emsfc_lw .gt. 1.) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospIN%emsfc_lw contains values out of range'
          Lisccp_subcolumn = .false.
          Lisccp_column    = .false.
          if (associated(cospOUT%isccp_totalcldarea))  cospOUT%isccp_totalcldarea(:)  = R_UNDEF
          if (associated(cospOUT%isccp_meantb))        cospOUT%isccp_meantb(:)        = R_UNDEF
          if (associated(cospOUT%isccp_meantbclr))     cospOUT%isccp_meantbclr(:)     = R_UNDEF
          if (associated(cospOUT%isccp_meanptop))      cospOUT%isccp_meanptop(:)      = R_UNDEF
          if (associated(cospOUT%isccp_meantaucld))    cospOUT%isccp_meantaucld(:)    = R_UNDEF
          if (associated(cospOUT%isccp_meanalbedocld)) cospOUT%isccp_meanalbedocld(:) = R_UNDEF
          if (associated(cospOUT%isccp_boxtau))        cospOUT%isccp_boxtau(:,:)      = R_UNDEF
          if (associated(cospOUT%isccp_boxptop))       cospOUT%isccp_boxptop(:,:)     = R_UNDEF
          if (associated(cospOUT%isccp_fq))            cospOUT%isccp_fq(:,:,:)        = R_UNDEF
       endif
    endif
    if (any([Lisccp_subcolumn,Lisccp_column,Lmisr_subcolumn,Lmisr_column,Lmodis_subcolumn,Lmodis_column])) then
       if (any(cospIN%tau_067 .lt. 0)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospIN%tau_067 contains values out of range'
          Lisccp_subcolumn = .false.
          Lisccp_column    = .false.
          Lmisr_subcolumn  = .false.
          Lmisr_column     = .false.
          Lmodis_subcolumn = .false.
          Lmodis_column    = .false.
          if (associated(cospOUT%isccp_totalcldarea))  cospOUT%isccp_totalcldarea(:)  = R_UNDEF
          if (associated(cospOUT%isccp_meantb))        cospOUT%isccp_meantb(:)        = R_UNDEF
          if (associated(cospOUT%isccp_meantbclr))     cospOUT%isccp_meantbclr(:)     = R_UNDEF
          if (associated(cospOUT%isccp_meanptop))      cospOUT%isccp_meanptop(:)      = R_UNDEF
          if (associated(cospOUT%isccp_meantaucld))    cospOUT%isccp_meantaucld(:)    = R_UNDEF
          if (associated(cospOUT%isccp_meanalbedocld)) cospOUT%isccp_meanalbedocld(:) = R_UNDEF
          if (associated(cospOUT%isccp_boxtau))        cospOUT%isccp_boxtau(:,:)      = R_UNDEF
          if (associated(cospOUT%isccp_boxptop))       cospOUT%isccp_boxptop(:,:)     = R_UNDEF
          if (associated(cospOUT%isccp_fq))            cospOUT%isccp_fq(:,:,:)        = R_UNDEF
          if (associated(cospOUT%misr_fq))                   cospOUT%misr_fq(:,:,:)                 = R_UNDEF
          if (associated(cospOUT%misr_dist_model_layertops)) cospOUT%misr_dist_model_layertops(:,:) = R_UNDEF
          if (associated(cospOUT%misr_meanztop))             cospOUT%misr_meanztop(:)               = R_UNDEF
          if (associated(cospOUT%misr_cldarea))              cospOUT%misr_cldarea(:)                = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_Total_Mean))                          &
               cospOUT%modis_Cloud_Fraction_Total_Mean(:)                   = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_Water_Mean))                          &
               cospOUT%modis_Cloud_Fraction_Water_Mean(:)                   = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_Ice_Mean))                            &
               cospOUT%modis_Cloud_Fraction_Ice_Mean(:)                     = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_High_Mean))                           &
               cospOUT%modis_Cloud_Fraction_High_Mean(:)                    = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_Mid_Mean))                            &
               cospOUT%modis_Cloud_Fraction_Mid_Mean(:)                     = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_Low_Mean))                            &
               cospOUT%modis_Cloud_Fraction_Low_Mean(:)                     = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Total_Mean))                       &
               cospOUT%modis_Optical_Thickness_Total_Mean(:)                = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Water_Mean))                       &
               cospOUT%modis_Optical_Thickness_Water_Mean(:)                = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Ice_Mean))                         &
               cospOUT%modis_Optical_Thickness_Ice_Mean(:)                  = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Total_LogMean))                    &
               cospOUT%modis_Optical_Thickness_Total_LogMean(:)             = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Water_LogMean))                    &
               cospOUT%modis_Optical_Thickness_Water_LogMean(:)             = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Ice_LogMean))                      &
               cospOUT%modis_Optical_Thickness_Ice_LogMean(:)               = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Particle_Size_Water_Mean))                     &
               cospOUT%modis_Cloud_Particle_Size_Water_Mean(:)              = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Particle_Size_Ice_Mean))                       &
               cospOUT%modis_Cloud_Particle_Size_Ice_Mean(:)                = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Top_Pressure_Total_Mean))                      &
               cospOUT%modis_Cloud_Top_Pressure_Total_Mean(:)               = R_UNDEF
          if (associated(cospOUT%modis_Liquid_Water_Path_Mean))                             &
               cospOUT%modis_Liquid_Water_Path_Mean(:)                      = R_UNDEF
          if (associated(cospOUT%modis_Ice_Water_Path_Mean))                                &
               cospOUT%modis_Ice_Water_Path_Mean(:)                         = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure))            &
               cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure(:,:,:) = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Liq))            &
               cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Liq(:,:,:) = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Ice))            &
               cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Ice(:,:,:) = R_UNDEF
          if (associated(cospOUT%modis_LWP_vs_ReffLIQ))            &
               cospOUT%modis_LWP_vs_ReffLIQ(:,:,:) = R_UNDEF
          if (associated(cospOUT%modis_IWP_vs_ReffICE))            &
               cospOUT%modis_IWP_vs_ReffICE(:,:,:) = R_UNDEF     
          if (associated(cospOUT%modis_Optical_Thickness_vs_ReffICE))                       &
               cospOUT%modis_Optical_Thickness_vs_ReffICE(:,:,:)            = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_vs_ReffLIQ))                       &
               cospOUT%modis_Optical_Thickness_vs_ReffLIQ(:,:,:)            = R_UNDEF
       endif
    endif
    if (any([Lisccp_subcolumn,Lisccp_column])) then
       if (any(cospIN%emiss_11 .lt. 0. .OR. cospIN%emiss_11 .gt. 1)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospIN%emiss_11 contains values out of range'
          Lisccp_subcolumn = .false.
          Lisccp_column    = .false.
          if (associated(cospOUT%isccp_totalcldarea))  cospOUT%isccp_totalcldarea(:)  = R_UNDEF
          if (associated(cospOUT%isccp_meantb))        cospOUT%isccp_meantb(:)        = R_UNDEF
          if (associated(cospOUT%isccp_meantbclr))     cospOUT%isccp_meantbclr(:)     = R_UNDEF
          if (associated(cospOUT%isccp_meanptop))      cospOUT%isccp_meanptop(:)      = R_UNDEF
          if (associated(cospOUT%isccp_meantaucld))    cospOUT%isccp_meantaucld(:)    = R_UNDEF
          if (associated(cospOUT%isccp_meanalbedocld)) cospOUT%isccp_meanalbedocld(:) = R_UNDEF
          if (associated(cospOUT%isccp_boxtau))        cospOUT%isccp_boxtau(:,:)      = R_UNDEF
          if (associated(cospOUT%isccp_boxptop))       cospOUT%isccp_boxptop(:,:)     = R_UNDEF
          if (associated(cospOUT%isccp_fq))            cospOUT%isccp_fq(:,:,:)        = R_UNDEF
       endif
    endif
    if (any([Lmodis_subcolumn,Lmodis_column])) then
       if (any(cospIN%asym .lt. -1. .OR. cospIN%asym .gt. 1)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospIN%asym contains values out of range'
          Lmodis_subcolumn = .false.
          Lmodis_column    = .false.
          if (associated(cospOUT%modis_Cloud_Fraction_Total_Mean))                          &
               cospOUT%modis_Cloud_Fraction_Total_Mean(:)                   = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_Water_Mean))                          &
               cospOUT%modis_Cloud_Fraction_Water_Mean(:)                   = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_Ice_Mean))                            &
               cospOUT%modis_Cloud_Fraction_Ice_Mean(:)                     = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_High_Mean))                           &
               cospOUT%modis_Cloud_Fraction_High_Mean(:)                    = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_Mid_Mean))                            &
               cospOUT%modis_Cloud_Fraction_Mid_Mean(:)                     = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_Low_Mean))                            &
               cospOUT%modis_Cloud_Fraction_Low_Mean(:)                     = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Total_Mean))                       &
               cospOUT%modis_Optical_Thickness_Total_Mean(:)                = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Water_Mean))                       &
               cospOUT%modis_Optical_Thickness_Water_Mean(:)                = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Ice_Mean))                         &
               cospOUT%modis_Optical_Thickness_Ice_Mean(:)                  = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Total_LogMean))                    &
               cospOUT%modis_Optical_Thickness_Total_LogMean(:)             = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Water_LogMean))                    &
               cospOUT%modis_Optical_Thickness_Water_LogMean(:)             = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Ice_LogMean))                      &
               cospOUT%modis_Optical_Thickness_Ice_LogMean(:)               = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Particle_Size_Water_Mean))                     &
               cospOUT%modis_Cloud_Particle_Size_Water_Mean(:)              = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Particle_Size_Ice_Mean))                       &
               cospOUT%modis_Cloud_Particle_Size_Ice_Mean(:)                = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Top_Pressure_Total_Mean))                      &
               cospOUT%modis_Cloud_Top_Pressure_Total_Mean(:)               = R_UNDEF
          if (associated(cospOUT%modis_Liquid_Water_Path_Mean))                             &
               cospOUT%modis_Liquid_Water_Path_Mean(:)                      = R_UNDEF
          if (associated(cospOUT%modis_Ice_Water_Path_Mean))                                &
               cospOUT%modis_Ice_Water_Path_Mean(:)                         = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure))            &
               cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure(:,:,:) = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Liq))            &
               cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Liq(:,:,:) = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Ice))            &
               cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Ice(:,:,:) = R_UNDEF
          if (associated(cospOUT%modis_LWP_vs_ReffLIQ))            &
               cospOUT%modis_LWP_vs_ReffLIQ(:,:,:) = R_UNDEF
          if (associated(cospOUT%modis_IWP_vs_ReffICE))            &
               cospOUT%modis_IWP_vs_ReffICE(:,:,:) = R_UNDEF    
          if (associated(cospOUT%modis_Optical_Thickness_vs_ReffICE))                       &
               cospOUT%modis_Optical_Thickness_vs_ReffICE(:,:,:)            = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_vs_ReffLIQ))                       &
               cospOUT%modis_Optical_Thickness_vs_ReffLIQ(:,:,:)            = R_UNDEF
       endif
       if (any(cospIN%ss_alb .lt. 0 .OR. cospIN%ss_alb .gt. 1)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospIN%ss_alb contains values out of range'
          Lmodis_subcolumn = .false.
          Lmodis_column    = .false.
          if (associated(cospOUT%modis_Cloud_Fraction_Total_Mean))                          &
               cospOUT%modis_Cloud_Fraction_Total_Mean(:)                   = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_Water_Mean))                          &
               cospOUT%modis_Cloud_Fraction_Water_Mean(:)                   = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_Ice_Mean))                            &
               cospOUT%modis_Cloud_Fraction_Ice_Mean(:)                     = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_High_Mean))                           &
               cospOUT%modis_Cloud_Fraction_High_Mean(:)                    = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_Mid_Mean))                            &
               cospOUT%modis_Cloud_Fraction_Mid_Mean(:)                     = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Fraction_Low_Mean))                            &
               cospOUT%modis_Cloud_Fraction_Low_Mean(:)                     = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Total_Mean))                       &
               cospOUT%modis_Optical_Thickness_Total_Mean(:)                = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Water_Mean))                       &
               cospOUT%modis_Optical_Thickness_Water_Mean(:)                = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Ice_Mean))                         &
               cospOUT%modis_Optical_Thickness_Ice_Mean(:)                  = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Total_LogMean))                    &
               cospOUT%modis_Optical_Thickness_Total_LogMean(:)             = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Water_LogMean))                    &
               cospOUT%modis_Optical_Thickness_Water_LogMean(:)             = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_Ice_LogMean))                      &
               cospOUT%modis_Optical_Thickness_Ice_LogMean(:)               = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Particle_Size_Water_Mean))                     &
               cospOUT%modis_Cloud_Particle_Size_Water_Mean(:)              = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Particle_Size_Ice_Mean))                       &
               cospOUT%modis_Cloud_Particle_Size_Ice_Mean(:)                = R_UNDEF
          if (associated(cospOUT%modis_Cloud_Top_Pressure_Total_Mean))                      &
               cospOUT%modis_Cloud_Top_Pressure_Total_Mean(:)               = R_UNDEF
          if (associated(cospOUT%modis_Liquid_Water_Path_Mean))                             &
               cospOUT%modis_Liquid_Water_Path_Mean(:)                      = R_UNDEF
          if (associated(cospOUT%modis_Ice_Water_Path_Mean))                                &
               cospOUT%modis_Ice_Water_Path_Mean(:)                         = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure))            &
               cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure(:,:,:) = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Liq))            &
               cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Liq(:,:,:) = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Ice))            &
               cospOUT%modis_Optical_Thickness_vs_Cloud_Top_Pressure_Ice(:,:,:) = R_UNDEF
          if (associated(cospOUT%modis_LWP_vs_ReffLIQ))            &
               cospOUT%modis_LWP_vs_ReffLIQ(:,:,:) = R_UNDEF
          if (associated(cospOUT%modis_IWP_vs_ReffICE))            &
               cospOUT%modis_IWP_vs_ReffICE(:,:,:) = R_UNDEF     
          if (associated(cospOUT%modis_Optical_Thickness_vs_ReffICE))                       &
               cospOUT%modis_Optical_Thickness_vs_ReffICE(:,:,:)            = R_UNDEF
          if (associated(cospOUT%modis_Optical_Thickness_vs_ReffLIQ))                       &
               cospOUT%modis_Optical_Thickness_vs_ReffLIQ(:,:,:)            = R_UNDEF
       endif
    endif
    if (any([Latlid_subcolumn,Latlid_column])) then
       if (any(cospIN%betatot_atlid .lt. 0)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospIN%betatot_atlid contains values out of range'
          Latlid_subcolumn = .false.
          Latlid_column    = .false.
          if (associated(cospOUT%atlid_cfad_sr))       cospOUT%atlid_cfad_sr(:,:,:)  = R_UNDEF
          if (associated(cospOUT%atlid_lidarcld))      cospOUT%atlid_lidarcld(:,:)   = R_UNDEF
          if (associated(cospOUT%atlid_cldlayer))      cospOUT%atlid_cldlayer(:,:)   = R_UNDEF
          if (associated(cospOUT%atlid_beta_tot))      cospOUT%atlid_beta_tot(:,:,:) = R_UNDEF
          if (associated(cospOUT%atlid_beta_mol))      cospOUT%atlid_beta_mol(:,:)   = R_UNDEF
       endif
       if (any(cospIN%beta_mol_atlid .lt. 0)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospIN%beta_mol_atlid contains values out of range'
          Latlid_subcolumn = .false.
          Latlid_column    = .false.
          if (associated(cospOUT%atlid_cfad_sr))       cospOUT%atlid_cfad_sr(:,:,:)  = R_UNDEF
          if (associated(cospOUT%atlid_lidarcld))      cospOUT%atlid_lidarcld(:,:)   = R_UNDEF
          if (associated(cospOUT%atlid_cldlayer))      cospOUT%atlid_cldlayer(:,:)   = R_UNDEF
          if (associated(cospOUT%atlid_beta_tot))      cospOUT%atlid_beta_tot(:,:,:) = R_UNDEF
          if (associated(cospOUT%atlid_beta_mol))      cospOUT%atlid_beta_mol(:,:)   = R_UNDEF
       endif
       if (any(cospIN%tautot_atlid .lt. 0)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospIN%tautot_atlid contains values out of range'
          Latlid_subcolumn = .false.
          Latlid_column    = .false.
          if (associated(cospOUT%atlid_cfad_sr))       cospOUT%atlid_cfad_sr(:,:,:)  = R_UNDEF
          if (associated(cospOUT%atlid_lidarcld))      cospOUT%atlid_lidarcld(:,:)   = R_UNDEF
          if (associated(cospOUT%atlid_cldlayer))      cospOUT%atlid_cldlayer(:,:)   = R_UNDEF
          if (associated(cospOUT%atlid_beta_tot))      cospOUT%atlid_beta_tot(:,:,:) = R_UNDEF
          if (associated(cospOUT%atlid_beta_mol))      cospOUT%atlid_beta_mol(:,:)   = R_UNDEF
       endif
       if (any(cospIN%tau_mol_atlid .lt. 0)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospIN%tau_mol_atlid contains values out of range'
          Latlid_subcolumn = .false.
          Latlid_column    = .false.
          if (associated(cospOUT%atlid_cfad_sr))       cospOUT%atlid_cfad_sr(:,:,:)  = R_UNDEF
          if (associated(cospOUT%atlid_lidarcld))      cospOUT%atlid_lidarcld(:,:)   = R_UNDEF
          if (associated(cospOUT%atlid_cldlayer))      cospOUT%atlid_cldlayer(:,:)   = R_UNDEF
          if (associated(cospOUT%atlid_beta_tot))      cospOUT%atlid_beta_tot(:,:,:) = R_UNDEF
          if (associated(cospOUT%atlid_beta_mol))      cospOUT%atlid_beta_mol(:,:)   = R_UNDEF
       endif
    endif
    
    if (any([LgrLidar532_subcolumn,LgrLidar532_column])) then
       if (any(cospIN%betatot_grLidar532 .lt. 0)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospIN%betatot_grLidar532 contains values out of range'
          LgrLidar532_subcolumn = .false.
          LgrLidar532_column    = .false.
          if (associated(cospOUT%grLidar532_cfad_sr))       cospOUT%grLidar532_cfad_sr(:,:,:)  = R_UNDEF
          if (associated(cospOUT%grLidar532_lidarcld))      cospOUT%grLidar532_lidarcld(:,:)   = R_UNDEF
          if (associated(cospOUT%grLidar532_cldlayer))      cospOUT%grLidar532_cldlayer(:,:)   = R_UNDEF
          if (associated(cospOUT%grLidar532_beta_tot))      cospOUT%grLidar532_beta_tot(:,:,:) = R_UNDEF
          if (associated(cospOUT%grLidar532_beta_mol))      cospOUT%grLidar532_beta_mol(:,:)   = R_UNDEF
       endif
       if (any(cospIN%beta_mol_grLidar532 .lt. 0)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospIN%beta_mol_grLidar532 contains values out of range'
          LgrLidar532_subcolumn = .false.
          LgrLidar532_column    = .false.
          if (associated(cospOUT%grLidar532_cfad_sr))       cospOUT%grLidar532_cfad_sr(:,:,:)  = R_UNDEF
          if (associated(cospOUT%grLidar532_lidarcld))      cospOUT%grLidar532_lidarcld(:,:)   = R_UNDEF
          if (associated(cospOUT%grLidar532_cldlayer))      cospOUT%grLidar532_cldlayer(:,:)   = R_UNDEF
          if (associated(cospOUT%grLidar532_beta_tot))      cospOUT%grLidar532_beta_tot(:,:,:) = R_UNDEF
          if (associated(cospOUT%grLidar532_beta_mol))      cospOUT%grLidar532_beta_mol(:,:)   = R_UNDEF
       endif
       if (any(cospIN%tautot_grLidar532 .lt. 0)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospIN%tautot_grLidar532 contains values out of range'
          LgrLidar532_subcolumn = .false.
          LgrLidar532_column    = .false.
          if (associated(cospOUT%grLidar532_cfad_sr))       cospOUT%grLidar532_cfad_sr(:,:,:)  = R_UNDEF
          if (associated(cospOUT%grLidar532_lidarcld))      cospOUT%grLidar532_lidarcld(:,:)   = R_UNDEF
          if (associated(cospOUT%grLidar532_cldlayer))      cospOUT%grLidar532_cldlayer(:,:)   = R_UNDEF
          if (associated(cospOUT%grLidar532_beta_tot))      cospOUT%grLidar532_beta_tot(:,:,:) = R_UNDEF
          if (associated(cospOUT%grLidar532_beta_mol))      cospOUT%grLidar532_beta_mol(:,:)   = R_UNDEF
       endif
       if (any(cospIN%tau_mol_grLidar532 .lt. 0)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospIN%tau_mol_grLidar532 contains values out of range'
          LgrLidar532_subcolumn = .false.
          LgrLidar532_column    = .false.
          if (associated(cospOUT%grLidar532_cfad_sr))       cospOUT%grLidar532_cfad_sr(:,:,:)  = R_UNDEF
          if (associated(cospOUT%grLidar532_lidarcld))      cospOUT%grLidar532_lidarcld(:,:)   = R_UNDEF
          if (associated(cospOUT%grLidar532_cldlayer))      cospOUT%grLidar532_cldlayer(:,:)   = R_UNDEF
          if (associated(cospOUT%grLidar532_beta_tot))      cospOUT%grLidar532_beta_tot(:,:,:) = R_UNDEF
          if (associated(cospOUT%grLidar532_beta_mol))      cospOUT%grLidar532_beta_mol(:,:)   = R_UNDEF
       endif
    endif

    if (any([Lcalipso_subcolumn,Lcalipso_column])) then
       if (any(cospIN%betatot_calipso .lt. 0)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospIN%betatot_calipso contains values out of range'
          Lcalipso_subcolumn = .false.
          Lcalipso_column    = .false.
          if (associated(cospOUT%calipso_cfad_sr))       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF
          if (associated(cospOUT%calipso_lidarcld))      cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldphase)) cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
          if (associated(cospOUT%calipso_cldlayer))      cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
          if (associated(cospOUT%calipso_cldlayerphase)) cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldtmp))   cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF
          if (associated(cospOUT%calipso_srbval))        cospOUT%calipso_srbval(:)            = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldtype))  cospOUT%calipso_lidarcldtype(:,:,:)  = R_UNDEF 
          if (associated(cospOUT%calipso_cldtype))       cospOUT%calipso_cldtype(:,:)         = R_UNDEF 
          if (associated(cospOUT%calipso_cldtypetemp))   cospOUT%calipso_cldtypetemp(:,:)     = R_UNDEF 
          if (associated(cospOUT%calipso_cldtypemeanz))  cospOUT%calipso_cldtypemeanz(:,:)    = R_UNDEF 
          if (associated(cospOUT%calipso_cldtypemeanzse)) cospOUT%calipso_cldtypemeanzse(:,:) = R_UNDEF 
          if (associated(cospOUT%calipso_cldthinemis))   cospOUT%calipso_cldthinemis(:)       = R_UNDEF
       endif
       if (any(cospIN%betatot_liq_calipso .lt. 0)) then
          nError=nError+1
          errorMessage(nError) = ('ERROR: COSP input variable: cospIN%betatot_liq_calipso contains values out of range')
          Lcalipso_subcolumn = .false.
          Lcalipso_column    = .false.
          if (associated(cospOUT%calipso_cfad_sr))       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF
          if (associated(cospOUT%calipso_lidarcld))      cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldphase)) cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
          if (associated(cospOUT%calipso_cldlayer))      cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
          if (associated(cospOUT%calipso_cldlayerphase)) cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldtmp))   cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF
          if (associated(cospOUT%calipso_srbval))        cospOUT%calipso_srbval(:)            = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldtype))  cospOUT%calipso_lidarcldtype(:,:,:)  = R_UNDEF
          if (associated(cospOUT%calipso_cldtype))       cospOUT%calipso_cldtype(:,:)         = R_UNDEF
          if (associated(cospOUT%calipso_cldtypetemp))   cospOUT%calipso_cldtypetemp(:,:)     = R_UNDEF 
          if (associated(cospOUT%calipso_cldtypemeanz))  cospOUT%calipso_cldtypemeanz(:,:)    = R_UNDEF 
          if (associated(cospOUT%calipso_cldtypemeanzse)) cospOUT%calipso_cldtypemeanzse(:,:) = R_UNDEF
          if (associated(cospOUT%calipso_cldthinemis))   cospOUT%calipso_cldthinemis(:)       = R_UNDEF
       endif
       if (any(cospIN%betatot_ice_calipso .lt. 0)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospIN%betatot_ice_calipso contains values out of range'
          Lcalipso_subcolumn = .false.
          Lcalipso_column    = .false.
          if (associated(cospOUT%calipso_cfad_sr))       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF
          if (associated(cospOUT%calipso_lidarcld))      cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldphase)) cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
          if (associated(cospOUT%calipso_cldlayer))      cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
          if (associated(cospOUT%calipso_cldlayerphase)) cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldtmp))   cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF
          if (associated(cospOUT%calipso_srbval))        cospOUT%calipso_srbval(:)            = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldtype))  cospOUT%calipso_lidarcldtype(:,:,:)  = R_UNDEF
          if (associated(cospOUT%calipso_cldtype))       cospOUT%calipso_cldtype(:,:)         = R_UNDEF 
          if (associated(cospOUT%calipso_cldtypetemp))   cospOUT%calipso_cldtypetemp(:,:)     = R_UNDEF 
          if (associated(cospOUT%calipso_cldtypemeanz))  cospOUT%calipso_cldtypemeanz(:,:)    = R_UNDEF 
          if (associated(cospOUT%calipso_cldtypemeanzse)) cospOUT%calipso_cldtypemeanzse(:,:) = R_UNDEF
          if (associated(cospOUT%calipso_cldthinemis))   cospOUT%calipso_cldthinemis(:)       = R_UNDEF
       endif
       if (any(cospIN%tautot_calipso .lt. 0)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospIN%tautot_calipso contains values out of range'
          Lcalipso_subcolumn = .false.
          Lcalipso_column    = .false.
          if (associated(cospOUT%calipso_cfad_sr))       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF
          if (associated(cospOUT%calipso_lidarcld))      cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldphase)) cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
          if (associated(cospOUT%calipso_cldlayer))      cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
          if (associated(cospOUT%calipso_cldlayerphase)) cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldtmp))   cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF
          if (associated(cospOUT%calipso_srbval))        cospOUT%calipso_srbval(:)            = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldtype))  cospOUT%calipso_lidarcldtype(:,:,:)  = R_UNDEF
          if (associated(cospOUT%calipso_cldtype))       cospOUT%calipso_cldtype(:,:)         = R_UNDEF 
          if (associated(cospOUT%calipso_cldtypetemp))   cospOUT%calipso_cldtypetemp(:,:)     = R_UNDEF 
          if (associated(cospOUT%calipso_cldtypemeanz))  cospOUT%calipso_cldtypemeanz(:,:)    = R_UNDEF 
          if (associated(cospOUT%calipso_cldtypemeanzse)) cospOUT%calipso_cldtypemeanzse(:,:) = R_UNDEF
          if (associated(cospOUT%calipso_cldthinemis))   cospOUT%calipso_cldthinemis(:)       = R_UNDEF
       endif
       if (any(cospIN%tautot_liq_calipso .lt. 0)) then
          nError=nError+1
          errorMessage(nError) = ('ERROR: COSP input variable: cospIN%tautot_liq_calipso contains values out of range')
          Lcalipso_subcolumn = .false.
          Lcalipso_column    = .false.
          if (associated(cospOUT%calipso_cfad_sr))       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF
          if (associated(cospOUT%calipso_lidarcld))      cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldphase)) cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
          if (associated(cospOUT%calipso_cldlayer))      cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
          if (associated(cospOUT%calipso_cldlayerphase)) cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldtmp))   cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF
          if (associated(cospOUT%calipso_srbval))        cospOUT%calipso_srbval(:)            = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldtype))  cospOUT%calipso_lidarcldtype(:,:,:)  = R_UNDEF
          if (associated(cospOUT%calipso_cldtype))       cospOUT%calipso_cldtype(:,:)         = R_UNDEF
          if (associated(cospOUT%calipso_cldtypetemp))   cospOUT%calipso_cldtypetemp(:,:)     = R_UNDEF
          if (associated(cospOUT%calipso_cldtypemeanz))  cospOUT%calipso_cldtypemeanz(:,:)    = R_UNDEF
          if (associated(cospOUT%calipso_cldtypemeanzse)) cospOUT%calipso_cldtypemeanzse(:,:) = R_UNDEF
          if (associated(cospOUT%calipso_cldthinemis))   cospOUT%calipso_cldthinemis(:)       = R_UNDEF
       endif
       if (any(cospIN%tautot_ice_calipso .lt. 0)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospIN%tautot_ice_calipso contains values out of range'
          Lcalipso_subcolumn = .false.
          Lcalipso_column    = .false.
          if (associated(cospOUT%calipso_cfad_sr))       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF
          if (associated(cospOUT%calipso_lidarcld))      cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldphase)) cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
          if (associated(cospOUT%calipso_cldlayer))      cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
          if (associated(cospOUT%calipso_cldlayerphase)) cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldtmp))   cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF
          if (associated(cospOUT%calipso_srbval))        cospOUT%calipso_srbval(:)            = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldtype))  cospOUT%calipso_lidarcldtype(:,:,:)  = R_UNDEF
          if (associated(cospOUT%calipso_cldtype))       cospOUT%calipso_cldtype(:,:)         = R_UNDEF
          if (associated(cospOUT%calipso_cldtypetemp))   cospOUT%calipso_cldtypetemp(:,:)     = R_UNDEF
          if (associated(cospOUT%calipso_cldtypemeanz))  cospOUT%calipso_cldtypemeanz(:,:)    = R_UNDEF
          if (associated(cospOUT%calipso_cldtypemeanzse)) cospOUT%calipso_cldtypemeanzse(:,:) = R_UNDEF
          if (associated(cospOUT%calipso_cldthinemis))   cospOUT%calipso_cldthinemis(:)       = R_UNDEF
       endif
       if (any(cospIN%tau_mol_calipso .lt. 0)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospIN%tau_mol_calipso contains values out of range'
          Lcalipso_subcolumn = .false.
          Lcalipso_column    = .false.
          if (associated(cospOUT%calipso_cfad_sr))       cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF
          if (associated(cospOUT%calipso_lidarcld))      cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldphase)) cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
          if (associated(cospOUT%calipso_cldlayer))      cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
          if (associated(cospOUT%calipso_cldlayerphase)) cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldtmp))   cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF
          if (associated(cospOUT%calipso_srbval))        cospOUT%calipso_srbval(:)            = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldtype))  cospOUT%calipso_lidarcldtype(:,:,:)  = R_UNDEF 
          if (associated(cospOUT%calipso_cldtype))       cospOUT%calipso_cldtype(:,:)         = R_UNDEF 
          if (associated(cospOUT%calipso_cldtypetemp))   cospOUT%calipso_cldtypetemp(:,:)     = R_UNDEF
          if (associated(cospOUT%calipso_cldtypemeanz))  cospOUT%calipso_cldtypemeanz(:,:)    = R_UNDEF
          if (associated(cospOUT%calipso_cldtypemeanzse)) cospOUT%calipso_cldtypemeanzse(:,:) = R_UNDEF
          if (associated(cospOUT%calipso_cldthinemis))   cospOUT%calipso_cldthinemis(:)       = R_UNDEF
       endif
    endif
    if (any([Lcalipso_subcolumn,Lcalipso_column,Lcloudsat_column,Lradar_lidar_tcc,       &
        Llidar_only_freq_cloud, Lcloudsat_tcc, Lcloudsat_tcc2])) then
       if (any(cospIN%beta_mol_calipso .lt. 0)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospIN%beta_mol_calipso contains values out of range'
          Lcalipso_subcolumn = .false.
          Lcalipso_column    = .false.
          Lcloudsat_column   = .false.
          Lradar_lidar_tcc   = .false.
          Llidar_only_freq_cloud = .false.
          Lcloudsat_tcc      = .false.
          Lcloudsat_tcc2     = .false.
          if (associated(cospOUT%calipso_cfad_sr))        cospOUT%calipso_cfad_sr(:,:,:)       = R_UNDEF
          if (associated(cospOUT%calipso_lidarcld))       cospOUT%calipso_lidarcld(:,:)        = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldphase))  cospOUT%calipso_lidarcldphase(:,:,:) = R_UNDEF
          if (associated(cospOUT%calipso_cldlayer))       cospOUT%calipso_cldlayer(:,:)        = R_UNDEF
          if (associated(cospOUT%calipso_cldlayerphase))  cospOUT%calipso_cldlayerphase(:,:,:) = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldtmp))    cospOUT%calipso_lidarcldtmp(:,:,:)   = R_UNDEF
          if (associated(cospOUT%calipso_srbval))         cospOUT%calipso_srbval(:)            = R_UNDEF
          if (associated(cospOUT%cloudsat_cfad_ze))       cospOUT%cloudsat_cfad_ze(:,:,:)      = R_UNDEF
          if (associated(cospOUT%lidar_only_freq_cloud))  cospOUT%lidar_only_freq_cloud(:,:)   = R_UNDEF
          if (associated(cospOUT%radar_lidar_tcc))        cospOUT%radar_lidar_tcc(:)           = R_UNDEF
          if (associated(cospOUT%cloudsat_tcc))           cospOUT%cloudsat_tcc(:)              = R_UNDEF
          if (associated(cospOUT%cloudsat_tcc2))          cospOUT%cloudsat_tcc2(:)             = R_UNDEF
          if (associated(cospOUT%calipso_lidarcldtype))   cospOUT%calipso_lidarcldtype(:,:,:)  = R_UNDEF
          if (associated(cospOUT%calipso_cldtype))        cospOUT%calipso_cldtype(:,:)         = R_UNDEF
          if (associated(cospOUT%calipso_cldtypetemp))    cospOUT%calipso_cldtypetemp(:,:)     = R_UNDEF
          if (associated(cospOUT%calipso_cldtypemeanz))   cospOUT%calipso_cldtypemeanz(:,:)    = R_UNDEF
          if (associated(cospOUT%calipso_cldtypemeanzse)) cospOUT%calipso_cldtypemeanzse(:,:)  = R_UNDEF
          if (associated(cospOUT%calipso_cldthinemis))    cospOUT%calipso_cldthinemis(:)       = R_UNDEF
       endif
    endif
    if (any([Lparasol_subcolumn,Lparasol_column])) then
       if (any(cospIN%tautot_S_liq .lt. 0)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospIN%tautot_S_liq contains values out of range'
          Lparasol_subcolumn = .false.
          Lparasol_column    = .false.
          if (associated(cospOUT%parasolPix_refl))  cospOUT%parasolPix_refl(:,:,:) = R_UNDEF
          if (associated(cospOUT%parasolGrid_refl)) cospOUT%parasolGrid_refl(:,:)  = R_UNDEF
       endif
       if (any(cospIN%tautot_S_ice .lt. 0)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospIN%tautot_S_ice contains values out of range'
          Lparasol_subcolumn = .false.
          Lparasol_column    = .false.
          if (associated(cospOUT%parasolPix_refl))  cospOUT%parasolPix_refl(:,:,:) = R_UNDEF
          if (associated(cospOUT%parasolGrid_refl)) cospOUT%parasolGrid_refl(:,:)  = R_UNDEF
       endif
    endif
    if (any([Lcloudsat_subcolumn,Lcloudsat_column,Lradar_lidar_tcc,Llidar_only_freq_cloud, &
        Lcloudsat_tcc, Lcloudsat_tcc2, Lcloudsat_modis_wr])) then
       if (any(cospIN%z_vol_cloudsat .lt. 0)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospIN%z_vol_cloudsat contains values out of range'
          Lcloudsat_subcolumn = .false.
          Lcloudsat_column    = .false.
          Lradar_lidar_tcc    = .false.
          Llidar_only_freq_cloud = .false.
          Lcloudsat_tcc       = .false.
          Lcloudsat_tcc2      = .false.
          Lcloudsat_modis_wr  = .false.
          if (associated(cospOUT%cloudsat_cfad_ze))          cospOUT%cloudsat_cfad_ze(:,:,:)        = R_UNDEF
          if (associated(cospOUT%cloudsat_Ze_tot))           cospOUT%cloudsat_Ze_tot(:,:,:)         = R_UNDEF
          if (associated(cospOUT%lidar_only_freq_cloud))     cospOUT%lidar_only_freq_cloud(:,:)     = R_UNDEF
          if (associated(cospOUT%radar_lidar_tcc))           cospOUT%radar_lidar_tcc(:)             = R_UNDEF
          if (associated(cospOUT%cloudsat_tcc)) cospOUT%cloudsat_tcc(:) = R_UNDEF
          if (associated(cospOUT%cloudsat_tcc2)) cospOUT%cloudsat_tcc2(:) = R_UNDEF
          if (associated(cospOUT%cfodd_ntotal)) cospOUT%cfodd_ntotal(:,:,:,:) = R_UNDEF
          if (associated(cospOUT%wr_occfreq_ntotal)) cospOUT%wr_occfreq_ntotal(:,:) = R_UNDEF
       endif
       if (any(cospIN%kr_vol_cloudsat .lt. 0)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospIN%kr_vol_cloudsat contains values out of range'
          Lcloudsat_subcolumn = .false.
          Lcloudsat_column    = .false.
          Lradar_lidar_tcc    = .false.
          Llidar_only_freq_cloud = .false.
          Lcloudsat_tcc       = .false.
          Lcloudsat_tcc2      = .false.
          Lcloudsat_modis_wr  = .false.
          if (associated(cospOUT%cloudsat_cfad_ze))          cospOUT%cloudsat_cfad_ze(:,:,:)        = R_UNDEF
          if (associated(cospOUT%cloudsat_Ze_tot))           cospOUT%cloudsat_Ze_tot(:,:,:)         = R_UNDEF
          if (associated(cospOUT%lidar_only_freq_cloud))     cospOUT%lidar_only_freq_cloud(:,:)     = R_UNDEF
          if (associated(cospOUT%radar_lidar_tcc))           cospOUT%radar_lidar_tcc(:)             = R_UNDEF
          if (associated(cospOUT%cloudsat_tcc)) cospOUT%cloudsat_tcc(:) = R_UNDEF
          if (associated(cospOUT%cloudsat_tcc2)) cospOUT%cloudsat_tcc2(:) = R_UNDEF
          if (associated(cospOUT%cfodd_ntotal)) cospOUT%cfodd_ntotal(:,:,:,:) = R_UNDEF
          if (associated(cospOUT%wr_occfreq_ntotal)) cospOUT%wr_occfreq_ntotal(:,:) = R_UNDEF
       endif
       if (any(cospIN%g_vol_cloudsat .lt. 0)) then
          nError=nError+1
          errorMessage(nError) = 'ERROR: COSP input variable: cospIN%g_vol_cloudsat contains values out of range'
          Lcloudsat_subcolumn = .false.
          Lcloudsat_column    = .false.
          Lradar_lidar_tcc    = .false.
          Llidar_only_freq_cloud = .false.
          Lcloudsat_tcc       = .false.
          Lcloudsat_tcc2      = .false.          
          Lcloudsat_modis_wr  = .false.
          if (associated(cospOUT%cloudsat_cfad_ze))          cospOUT%cloudsat_cfad_ze(:,:,:)        = R_UNDEF
          if (associated(cospOUT%cloudsat_Ze_tot))           cospOUT%cloudsat_Ze_tot(:,:,:)         = R_UNDEF
          if (associated(cospOUT%lidar_only_freq_cloud))     cospOUT%lidar_only_freq_cloud(:,:)     = R_UNDEF
          if (associated(cospOUT%radar_lidar_tcc))           cospOUT%radar_lidar_tcc(:)             = R_UNDEF
          if (associated(cospOUT%cloudsat_tcc)) cospOUT%cloudsat_tcc(:) = R_UNDEF
          if (associated(cospOUT%cloudsat_tcc2)) cospOUT%cloudsat_tcc2(:) = R_UNDEF          
          if (associated(cospOUT%cfodd_ntotal)) cospOUT%cfodd_ntotal(:,:,:,:) = R_UNDEF
          if (associated(cospOUT%wr_occfreq_ntotal)) cospOUT%wr_occfreq_ntotal(:,:) = R_UNDEF
       endif
    endif
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Part 2: Check input fields array size for consistency. This needs to be done for each
    !         simulator
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! ISCCP
    if (Lisccp_subcolumn .or. Lisccp_column) then
       if (size(cospIN%frac_out,1)  .ne. cospIN%Npoints .OR. &
           size(cospIN%tau_067,1)   .ne. cospIN%Npoints .OR. &
           size(cospIN%emiss_11,1)  .ne. cospIN%Npoints .OR. &
           size(cospgridIN%skt)     .ne. cospIN%Npoints .OR. &
           size(cospgridIN%qv,1)    .ne. cospIN%Npoints .OR. &
           size(cospgridIN%at,1)    .ne. cospIN%Npoints .OR. &
           size(cospgridIN%phalf,1) .ne. cospIN%Npoints .OR. &
           size(cospgridIN%sunlit)  .ne. cospIN%Npoints .OR. &
           size(cospgridIN%pfull,1) .ne. cospIN%Npoints) then
          Lisccp_subcolumn = .false.
          Lisccp_column    = .false.
          nError=nError+1
          errorMessage(nError) = 'ERROR(isccp_simulator): The number of points in the input fields are inconsistent'
       endif
       if (size(cospIN%frac_out,2) .ne. cospIN%Ncolumns .OR. &
           size(cospIN%tau_067,2)  .ne. cospIN%Ncolumns .OR. &
           size(cospIN%emiss_11,2) .ne. cospIN%Ncolumns) then
          Lisccp_subcolumn = .false.
          Lisccp_column    = .false.
          nError=nError+1
          errorMessage(nError) = 'ERROR(isccp_simulator): The number of sub-columns in the input fields are inconsistent'
       endif
       if (size(cospIN%frac_out,3)  .ne. cospIN%Nlevels .OR. &
           size(cospIN%tau_067,3)   .ne. cospIN%Nlevels .OR. &
           size(cospIN%emiss_11,3)  .ne. cospIN%Nlevels .OR. &
           size(cospgridIN%qv,2)    .ne. cospIN%Nlevels .OR. &
           size(cospgridIN%at,2)    .ne. cospIN%Nlevels .OR. &
           size(cospgridIN%pfull,2) .ne. cospIN%Nlevels .OR. &
           size(cospgridIN%phalf,2) .ne. cospIN%Nlevels+1) then
          Lisccp_subcolumn = .false.
          Lisccp_column    = .false.
          nError=nError+1
          errorMessage(nError) = 'ERROR(isccp_simulator): The number of levels in the input fields are inconsistent'
       endif
    endif
    
    ! MISR
    if (Lmisr_subcolumn .or. Lmisr_column) then
       if (size(cospIN%tau_067,1)        .ne. cospIN%Npoints .OR. &
           size(cospgridIN%sunlit)       .ne. cospIN%Npoints .OR. &
           size(cospgridIN%hgt_matrix,1) .ne. cospIN%Npoints .OR. &
           size(cospgridIN%at,1)         .ne. cospIN%Npoints) then
          Lmisr_subcolumn = .false.
          Lmisr_column    = .false.
          nError=nError+1
          errorMessage(nError) = 'ERROR(misr_simulator): The number of points in the input fields are inconsistent'
       endif
       if (size(cospIN%tau_067,2) .ne. cospIN%Ncolumns) then
          Lmisr_subcolumn = .false.
          Lmisr_column    = .false.
          nError=nError+1
          errorMessage(nError) = 'ERROR(misr_simulator): The number of sub-columns in the input fields are inconsistent'
       endif
       if (size(cospIN%tau_067,3)        .ne. cospIN%Nlevels .OR. &
           size(cospgridIN%hgt_matrix,2) .ne. cospIN%Nlevels .OR. &
           size(cospgridIN%at,2)         .ne. cospIN%Nlevels) then
          Lmisr_subcolumn = .false.
          Lmisr_column    = .false.
          nError=nError+1
          errorMessage(nError) = 'ERROR(misr_simulator): The number of levels in the input fields are inconsistent'
       endif
    endif
    
    ! MODIS
    if (Lmodis_subcolumn .or. Lmodis_column) then
       if (size(cospIN%fracLiq,1) .ne. cospIN%Npoints .OR. &
           size(cospIN%tau_067,1) .ne. cospIN%Npoints .OR. &
           size(cospIN%asym,1)    .ne. cospIN%Npoints .OR. &
           size(cospIN%ss_alb,1)  .ne. cospIN%Npoints) then
          Lmodis_subcolumn = .false.
          Lmodis_column    = .false.
          nError=nError+1
          errorMessage(nError) = 'ERROR(modis_simulator): The number of points in the input fields are inconsistent'
       endif
       if (size(cospIN%fracLiq,2) .ne. cospIN%Ncolumns .OR. &
           size(cospIN%tau_067,2) .ne. cospIN%Ncolumns .OR. &
           size(cospIN%asym,2)    .ne. cospIN%Ncolumns .OR. &
           size(cospIN%ss_alb,2)  .ne. cospIN%Ncolumns) then
          Lmodis_subcolumn = .false.
          Lmodis_column    = .false.
          nError=nError+1
          errorMessage(nError) = 'ERROR(modis_simulator): The number of sub-columns in the input fields are inconsistent'
       endif
       if (size(cospIN%fracLiq,3) .ne. cospIN%Nlevels .OR. &
           size(cospIN%tau_067,3) .ne. cospIN%Nlevels .OR. &
           size(cospIN%asym,3)    .ne. cospIN%Nlevels .OR. &
           size(cospIN%ss_alb,3)  .ne. cospIN%Nlevels) then
          Lmodis_subcolumn = .false.
          Lmodis_column    = .false.
          nError=nError+1
          errorMessage(nError) = 'ERROR(modis_simulator): The number of levels in the input fields are inconsistent'
       endif
    endif
    
    ! CLOUDSAT
    if (Lcloudsat_subcolumn .or. Lcloudsat_column) then
       if (size(cospIN%z_vol_cloudsat,1)   .ne. cospIN%Npoints .OR. &
           size(cospIN%kr_vol_cloudsat,1)  .ne. cospIN%Npoints .OR. &
           size(cospIN%g_vol_cloudsat,1)   .ne. cospIN%Npoints .OR. &
           size(cospgridIN%hgt_matrix,1)   .ne. cospIN%Npoints) then
          Lcloudsat_subcolumn = .false.
          Lcloudsat_column    = .false.
          nError=nError+1
          errorMessage(nError) = 'ERROR(cloudsat_simulator): The number of points in the input fields are inconsistent'
       endif
       if (size(cospIN%z_vol_cloudsat,2)  .ne. cospIN%Ncolumns .OR. &
           size(cospIN%kr_vol_cloudsat,2) .ne. cospIN%Ncolumns .OR. &
           size(cospIN%g_vol_cloudsat,2)  .ne. cospIN%Ncolumns) then
          Lcloudsat_subcolumn = .false.
          Lcloudsat_column    = .false.
          nError=nError+1
          errorMessage(nError) = 'ERROR(cloudsat_simulator): The number of sub-columns in the input fields are inconsistent'
       endif
       if (size(cospIN%z_vol_cloudsat,3)  .ne. cospIN%Nlevels .OR. &
           size(cospIN%kr_vol_cloudsat,3) .ne. cospIN%Nlevels .OR. &
           size(cospIN%g_vol_cloudsat,3)  .ne. cospIN%Nlevels .OR. &
           size(cospgridIN%hgt_matrix,2)  .ne. cospIN%Nlevels) then
          Lcloudsat_subcolumn = .false.
          Lcloudsat_column    = .false.
          nError=nError+1
          errorMessage(nError) = 'ERROR(cloudsat_simulator): The number of levels in the input fields are inconsistent'
       endif
    endif

    ! GROUND LIDAR @ 532nm
    if (LgrLidar532_subcolumn .or. LgrLidar532_column) then
       if (size(cospIN%beta_mol_grLidar532,1)    .ne. cospIN%Npoints .OR. & 
           size(cospIN%betatot_grLidar532,1)     .ne. cospIN%Npoints .OR. &
           size(cospIN%tau_mol_grLidar532,1)     .ne. cospIN%Npoints .OR. &
           size(cospIN%tautot_grLidar532,1)      .ne. cospIN%Npoints) then
          LgrLidar532_subcolumn = .false.
          LgrLidar532_column    = .false.
          nError=nError+1
          errorMessage(nError) = 'ERROR(grLidar532_simulator): The number of points in the input fields are inconsistent'
       endif
       if (size(cospIN%betatot_grLidar532,2)    .ne. cospIN%Ncolumns .OR. & 
           size(cospIN%tautot_grLidar532,2)     .ne. cospIN%Ncolumns) then 
          LgrLidar532_subcolumn = .false. 
          LgrLidar532_column    = .false. 
          nError=nError+1
          errorMessage(nError) = 'ERROR(grLidar532_simulator): The number of sub-columns in the input fields are inconsistent'
       endif
       if (size(cospIN%beta_mol_grLidar532,2)    .ne. cospIN%Nlevels .OR. &
           size(cospIN%betatot_grLidar532,3)     .ne. cospIN%Nlevels .OR. &
           size(cospIN%tau_mol_grLidar532,2)     .ne. cospIN%Nlevels .OR. &
           size(cospIN%tautot_grLidar532,3)      .ne. cospIN%Nlevels) then
          LgrLidar532_subcolumn = .false. 
          LgrLidar532_column    = .false. 
          nError=nError+1
          errorMessage(nError) = 'ERROR(grLidar532_simulator): The number of levels in the input fields are inconsistent' 
       endif
    endif
    
    ! ATLID
    if (Latlid_subcolumn .or. Latlid_column) then
       if (size(cospIN%beta_mol_atlid,1)    .ne. cospIN%Npoints .OR. &
           size(cospIN%betatot_atlid,1)     .ne. cospIN%Npoints .OR. &
           size(cospIN%tau_mol_atlid,1)     .ne. cospIN%Npoints .OR. & 
           size(cospIN%tautot_atlid,1)      .ne. cospIN%Npoints) then 
          Latlid_subcolumn = .false. 
          Latlid_column    = .false. 
          nError=nError+1             
          errorMessage(nError) = 'ERROR(atlid_simulator): The number of points in the input fields are inconsistent'
       endif
       if (size(cospIN%betatot_atlid,2)    .ne. cospIN%Ncolumns .OR. &
           size(cospIN%tautot_atlid,2)     .ne. cospIN%Ncolumns) then 
          Latlid_subcolumn = .false.
          Latlid_column    = .false.
          nError=nError+1              
          errorMessage(nError) = 'ERROR(atlid_simulator): The number of sub-columns in the input fields are inconsistent'
       endif
       if (size(cospIN%beta_mol_atlid,2)    .ne. cospIN%Nlevels .OR. &
           size(cospIN%betatot_atlid,3)     .ne. cospIN%Nlevels .OR. & 
           size(cospIN%tau_mol_atlid,2)     .ne. cospIN%Nlevels .OR. &
           size(cospIN%tautot_atlid,3)      .ne. cospIN%Nlevels) then 
          Latlid_subcolumn = .false.
          Latlid_column    = .false. 
          nError=nError+1 
          errorMessage(nError) = 'ERROR(atlid_simulator): The number of levels in the input fields are inconsistent'
       endif
    endif

    ! CALIPSO
    if (Lcalipso_subcolumn .or. Lcalipso_column) then
       if (size(cospIN%beta_mol_calipso,1)    .ne. cospIN%Npoints .OR. &
           size(cospIN%betatot_calipso,1)     .ne. cospIN%Npoints .OR. &
           size(cospIN%betatot_liq_calipso,1) .ne. cospIN%Npoints .OR. &
           size(cospIN%betatot_ice_calipso,1) .ne. cospIN%Npoints .OR. &
           size(cospIN%tau_mol_calipso,1)     .ne. cospIN%Npoints .OR. &
           size(cospIN%tautot_calipso,1)      .ne. cospIN%Npoints .OR. &
           size(cospIN%tautot_liq_calipso,1)  .ne. cospIN%Npoints .OR. &
           size(cospIN%tautot_ice_calipso,1)  .ne. cospIN%Npoints) then
          Lcalipso_subcolumn = .false.
          Lcalipso_column    = .false.
          nError=nError+1
          errorMessage(nError) = 'ERROR(calipso_simulator): The number of points in the input fields are inconsistent'
       endif
       if (size(cospIN%betatot_calipso,2)     .ne. cospIN%Ncolumns .OR. &
           size(cospIN%betatot_liq_calipso,2) .ne. cospIN%Ncolumns .OR. &
           size(cospIN%betatot_ice_calipso,2) .ne. cospIN%Ncolumns .OR. &
           size(cospIN%tautot_calipso,2)      .ne. cospIN%Ncolumns .OR. &
           size(cospIN%tautot_liq_calipso,2)  .ne. cospIN%Ncolumns .OR. &
           size(cospIN%tautot_ice_calipso,2)  .ne. cospIN%Ncolumns) then
          Lcalipso_subcolumn = .false.
          Lcalipso_column    = .false.
          nError=nError+1
          errorMessage(nError) = 'ERROR(calipso_simulator): The number of sub-columns in the input fields are inconsistent'
       endif
       if (size(cospIN%beta_mol_calipso,2)    .ne. cospIN%Nlevels .OR. &
           size(cospIN%betatot_calipso,3)     .ne. cospIN%Nlevels .OR. &
           size(cospIN%betatot_liq_calipso,3) .ne. cospIN%Nlevels .OR. &
           size(cospIN%betatot_ice_calipso,3) .ne. cospIN%Nlevels .OR. &
           size(cospIN%tau_mol_calipso,2)     .ne. cospIN%Nlevels .OR. &
           size(cospIN%tautot_calipso,3)      .ne. cospIN%Nlevels .OR. &
           size(cospIN%tautot_liq_calipso,3)  .ne. cospIN%Nlevels .OR. &
           size(cospIN%tautot_ice_calipso,3)  .ne. cospIN%Nlevels) then
          Lcalipso_subcolumn = .false.
          Lcalipso_column    = .false.
          nError=nError+1
          errorMessage(nError) = 'ERROR(calipso_simulator): The number of levels in the input fields are inconsistent'
       endif
    endif
    
    ! PARASOL
    if (Lparasol_subcolumn .or. Lparasol_column) then
       if (size(cospIN%tautot_S_liq,1) .ne. cospIN%Npoints .OR. &
           size(cospIN%tautot_S_ice,1) .ne. cospIN%Npoints) then
          Lparasol_subcolumn = .false.
          Lparasol_column    = .false.
          nError=nError+1
          errorMessage(nError) = 'ERROR(parasol_simulator): The number of points in the input fields are inconsistent'
       endif
       if (size(cospIN%tautot_S_liq,2) .ne. cospIN%Ncolumns .OR. &
           size(cospIN%tautot_S_ice,2) .ne. cospIN%Ncolumns) then
          Lparasol_subcolumn = .false.
          Lparasol_column    = .false.
          nError=nError+1
          errorMessage(nError) = 'ERROR(parasol_simulator): The number of levels in the input fields are inconsistent'
       endif
    endif

    ! RTTOV
    if (Lrttov_column) then
       if (size(cospgridIN%pfull,1)           .ne. cospIN%Npoints .OR. &
           size(cospgridIN%at,1)              .ne. cospIN%Npoints .OR. &
           size(cospgridIN%qv,1)              .ne. cospIN%Npoints .OR. &
           size(cospgridIN%hgt_matrix_half,1) .ne. cospIN%Npoints .OR. &
           size(cospgridIN%u_sfc)             .ne. cospIN%Npoints .OR. &
           size(cospgridIN%v_sfc)             .ne. cospIN%Npoints .OR. &
           size(cospgridIN%skt)               .ne. cospIN%Npoints .OR. &
           size(cospgridIN%phalf,1)           .ne. cospIN%Npoints .OR. &
           size(cospgridIN%cloudIce,1)        .ne. cospIN%Npoints .OR. &
           size(cospgridIN%cloudLiq,1)        .ne. cospIN%Npoints .OR. &
           size(cospgridIN%rttov_sfcmask)     .ne. cospIN%Npoints .OR. &
           size(cospgridIN%lon)               .ne. cospIN%Npoints .OR. &
           size(cospgridIN%lat)               .ne. cospIN%Npoints) then
          Lrttov_column    = .false.
          nError=nError+1
          errorMessage(nError) = 'ERROR(rttov_simulator): The number of points in the input fields are inconsistent'
       endif
       if (size(cospgridIN%pfull,2)           .ne. cospIN%Nlevels   .OR. &
           size(cospgridIN%at,2)              .ne. cospIN%Nlevels   .OR. &
           size(cospgridIN%qv,2)              .ne. cospIN%Nlevels   .OR. &
           size(cospgridIN%hgt_matrix_half,2) .ne. cospIN%Nlevels   .OR. &
           size(cospgridIN%phalf,2)           .ne. cospIN%Nlevels+1) then
          Lrttov_column    = .false.
          nError=nError+1
          errorMessage(nError) = 'ERROR(rttov_simulator): The number of levels in the input fields are inconsistent'
       endif
    endif
  end subroutine cosp_errorCheck

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! END MODULE
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE MOD_COSP

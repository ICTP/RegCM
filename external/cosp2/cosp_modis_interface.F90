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
! History
! May 2015 - D. Swales - Original version
! Jun 2025 - J.K. Shaw - Added COSP-RTTOV integration and swathing
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE MOD_COSP_Modis_INTERFACE
  USE COSP_KINDS,      ONLY: wp
  USE MOD_COSP_CONFIG, ONLY: R_UNDEF
  use mod_modis_sim,   ONLY: num_trial_res,min_OpticalThickness,CO2Slicing_PressureLimit,&
                             CO2Slicing_TauLimit,phase_TauLimit,size_TauLimit,re_fill,   &
                             phaseDiscrimination_Threshold,re_water_min,     &
                             re_water_max,re_ice_min,re_ice_max,               &
                             highCloudPressureLimit,lowCloudPressureLimit,phaseIsNone,   &
                             phaseIsLiquid,phaseIsIce,phaseIsUndetermined,trial_re_w,    &
                             trial_re_i,g_w,g_i,w0_w,w0_i, get_g_nir,get_ssa_nir
  use mod_cosp_stats,  ONLY: compute_orbitmasks,cosp_optical_inputs,cosp_column_inputs
  implicit none

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !  TYPE modis_in
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  type modis_IN
     integer,pointer :: &
          Npoints,        & ! Number of horizontal gridpoints
          Ncolumns,       & ! Number of subcolumns
          Nlevels           ! Number of vertical levels
     integer :: &
          Nsunlit           ! Number of sunlit lit pixels
     real(wp),allocatable,dimension(:) :: &
          sunlit,         & ! Sunlit scenes
          notSunlit         ! Dark scenes
     real(wp),allocatable,dimension(:,:) :: &
          pres              ! Gridmean pressure at layer edges (Pa) 
     real(wp),pointer ::  &
          tau(:,:,:),     & ! Subcolumn optical thickness @ 0.67 microns.
          liqFrac(:,:,:), & ! Liquid water fraction
          g(:,:,:),       & ! Subcolumn assymetry parameter  
          w0(:,:,:)         ! Subcolumn single-scattering albedo
  end type modis_IN
contains
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROTUINE cosp_modis_init
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_MODIS_INIT()
    integer :: i 
    
    ! Retrieval parameters
    min_OpticalThickness          = 0.3_wp     ! Minimum detectable optical thickness
    CO2Slicing_PressureLimit      = 70000._wp  ! Cloud with higher pressures use thermal 
                                               ! methods, units Pa
    CO2Slicing_TauLimit           = 1._wp      ! How deep into the cloud does CO2 slicing 
                                               ! see? 
    phase_TauLimit                = 1._wp      ! How deep into the cloud does the phase 
                                               ! detection see?
    size_TauLimit                 = 2._wp      ! Depth of the re retreivals
    phaseDiscrimination_Threshold = 0.7_wp     ! What fraction of total extincton needs to 
                                               ! be in a single category to make phase 
                                               ! discrim. work? 
    re_fill                       = -999._wp   ! Fill value
    re_water_min                  = 4._wp      ! Minimum effective radius (liquid)
    re_water_max                  = 30._wp     ! Maximum effective radius (liquid)
    re_ice_min                    = 5._wp      ! Minimum effective radius (ice)
    re_ice_max                    = 90._wp     ! Minimum effective radius (ice)
    highCloudPressureLimit        = 44000._wp  ! High cloud pressure limit (Pa)
    lowCloudPressureLimit         = 68000._wp  ! Low cloud pressure limit (Pa)
    phaseIsNone                   = 0          ! No cloud
    phaseIsLiquid                 = 1          ! Liquid cloud
    phaseIsIce                    = 2          ! Ice cloud
    phaseIsUndetermined           = 3          ! Undetermined cloud

    ! Precompute near-IR optical params vs size for retrieval scheme    
    trial_re_w(1:num_trial_res) = re_water_min + (re_water_max - re_water_min) /         &
         (num_trial_res-1) * (/(i, i=0, num_trial_res-1)/)
    trial_re_i(1:num_trial_res) = re_ice_min   + (re_ice_max -   re_ice_min) /           &
         (num_trial_res-1) * (/(i, i=0, num_trial_res-1)/)

    ! Initialize estimates for size retrieval
    g_w(1:num_trial_res)  = get_g_nir(phaseIsLiquid,trial_re_w(1:num_trial_res))
    w0_w(1:num_trial_res) = get_ssa_nir(phaseIsLiquid,trial_re_w(1:num_trial_res))
    g_i(1:num_trial_res)  = get_g_nir(phaseIsIce,trial_re_i(1:num_trial_res))
    w0_i(1:num_trial_res) = get_ssa_nir(phaseIsIce,trial_re_i(1:num_trial_res))

  END SUBROUTINE COSP_MODIS_INIT

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !  							SUBROUTINE COSP_ASSIGN_modisIN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_ASSIGN_modisIN(cospIN,cospgridIN,Npoints,modisIN,CSCAL_SWATH_MASK,MODIS_CSCAL_MASK_INDICES)
     type(cosp_optical_inputs),intent(in),target :: cospIN     ! Optical inputs to COSP simulator
     type(cosp_column_inputs), intent(in),target :: cospgridIN ! Host model inputs to COSP
     integer,intent(in),target :: &
         Npoints
     type(modis_IN),intent(inout) :: &
         modisIN
     logical,dimension(:),allocatable,intent(in) :: & ! Mask of reals over all local times
         CSCAL_SWATH_MASK
     integer,dimension(:),allocatable,intent(out) :: & ! Array containing the indices of the swath masks, already allocated?
         MODIS_CSCAL_MASK_INDICES
     ! Local variables
     logical,dimension(:),allocatable :: & ! Mask of reals over all local times
         MODIS_SWATH_MASK, &
         MODIS_CSCAL_SWATH_MASK
     integer, target :: &
         N_MODIS_SWATHED,  &
         i

     if (cospIN % cospswathsIN(6) % N_inst_swaths .gt. 0) then
         allocate(MODIS_SWATH_MASK(Npoints))
         ! Do swathing to figure out which cells to simulate on
         call compute_orbitmasks(Npoints,                                                &
                                 cospIN % cospswathsIN(6) % N_inst_swaths,               &
                                 cospIN % cospswathsIN(6) % inst_localtimes,             &
                                 cospIN % cospswathsIN(6) % inst_localtime_widths,       &
                                 cospgridIN%lat, cospgridIN%lon,                         &
                                 cospgridIN%rttov_date(:,2), cospgridIN%rttov_date(:,3), & ! Time fields: month, dayofmonth
                                 cospgridIN%rttov_time(:,1), cospgridIN%rttov_time(:,2), & ! Time fields: hour, minute
                                 MODIS_SWATH_MASK,N_MODIS_SWATHED) ! Output: logical mask array

         ! Operate a little differently with MODIS because there is already some swathing at play.
         ! modisIN is over all variables rather than just swathed variables
         modisIN%Ncolumns  => cospIN%Ncolumns
         modisIN%Nlevels   => cospIN%Nlevels
         modisIN%Npoints   => Npoints
         modisIN%liqFrac   => cospIN%fracLiq
         modisIN%tau       => cospIN%tau_067
         modisIN%g         => cospIN%asym
         modisIN%w0        => cospIN%ss_alb

         modisIN%Nsunlit   = count((cospgridIN%sunlit > 0) .and. MODIS_SWATH_MASK) ! Sunlit mask and indices now include swathing as well
         if (modisIN%Nsunlit .gt. 0) then
            allocate(modisIN%sunlit(modisIN%Nsunlit), modisIN%pres(modisIN%Npoints,cospIN%Nlevels+1))
            modisIN%pres      = cospgridIN%phalf
            modisIN%sunlit    = pack((/ (i, i = 1, modisIN%Npoints ) /),mask = ((cospgridIN%sunlit > 0) .and. MODIS_SWATH_MASK)) ! Indices of columns to operate on in modisIN
         endif          
         if (modisIN%Npoints - modisIN%Nsunlit .gt. 0) then ! If more than zero tiles are not sunlit and swathed, create array to mask out these gridcells in cospOUT
            allocate(modisIN%notSunlit(modisIN%Npoints - modisIN%Nsunlit))
            modisIN%notSunlit = pack((/ (i, i = 1, Npoints ) /),mask = (.not. ((cospgridIN%sunlit > 0) .and. MODIS_SWATH_MASK))) ! Create an array with the indices of the non-sunlit tiles
         endif  
         ! Create a CSCAL-MODIS joint mask for the combined product.
         if (allocated(CSCAL_SWATH_MASK)) then
            allocate(MODIS_CSCAL_SWATH_MASK(Npoints))
            MODIS_CSCAL_SWATH_MASK = (.not. (MODIS_SWATH_MASK .and. CSCAL_SWATH_MASK)) ! Gridcells not seen by both MODIS and CSCAL should be set to zero
            MODIS_CSCAL_MASK_INDICES = pack((/ (i, i = 1, Npoints ) /),mask = MODIS_CSCAL_SWATH_MASK)
         end if
      else
         modisIN%Ncolumns  => cospIN%Ncolumns
         modisIN%Nlevels   => cospIN%Nlevels
         modisIN%Npoints   => Npoints
         modisIN%liqFrac   => cospIN%fracLiq
         modisIN%tau       => cospIN%tau_067
         modisIN%g         => cospIN%asym
         modisIN%w0        => cospIN%ss_alb        
         modisIN%Nsunlit   = count(cospgridIN%sunlit > 0)
         if (modisIN%Nsunlit .gt. 0) then
            allocate(modisIN%sunlit(modisIN%Nsunlit),modisIN%pres(modisIN%Npoints,cospIN%Nlevels+1))
            modisIN%pres      = cospgridIN%phalf
            modisIN%sunlit    = pack((/ (i, i = 1, Npoints ) /),mask = cospgridIN%sunlit > 0)
         endif       
         if (count(cospgridIN%sunlit <= 0) .gt. 0) then ! If more than zero tiles are not sunlit a.k.a. if there are dark tiles
            allocate(modisIN%notSunlit(count(cospgridIN%sunlit <= 0)))
            modisIN%notSunlit = pack((/ (i, i = 1, Npoints ) /),mask = .not. cospgridIN%sunlit > 0) ! Create an array with the indices of the non-sunlit tiles
         endif               
     end if

  END SUBROUTINE COSP_ASSIGN_modisIN

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! END MODULE MOD_COSP_Modis_INTERFACE
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE MOD_COSP_Modis_INTERFACE

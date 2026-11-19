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
MODULE MOD_COSP_CLOUDSAT_INTERFACE
  USE COSP_KINDS,      ONLY: wp
  USE quickbeam,       ONLY: quickbeam_init,Re_MAX_BIN,Re_BIN_LENGTH, &                                                                                                                                                       
                             maxhclass, nRe_types, nd, mt_ntt
  use mod_cosp_stats,  ONLY: radar_cfg,compute_orbitmasks,cosp_optical_inputs,cosp_column_inputs
  IMPLICIT NONE

  ! Directory where LUTs will be stored
  character(len=120) :: RADAR_SIM_LUT_DIRECTORY = './'
  logical :: RADAR_SIM_LOAD_scale_LUTs_flag   = .false.
  logical :: RADAR_SIM_UPDATE_scale_LUTs_flag = .false.
  
  ! Module variables
  real(wp),dimension(:,:),target,allocatable :: &
     temp_hgt_matrix
  real(wp),dimension(:,:,:),target,allocatable :: &
     temp_z_vol_cloudsat,   &
     temp_kr_vol_cloudsat,  &
     temp_g_vol_cloudsat

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! TYPE cloudsat_IN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  type cloudsat_IN
     integer         ::            &
          Npoints            ! Number of horizontal grid-points
     integer,pointer ::            &
          Nlevels,         & ! Number of vertical levels
          Ncolumns           ! Number of subcolumns
     real(wp),pointer ::   &
          hgt_matrix(:,:),   & ! Height of hydrometeors (km)
          z_vol(:,:,:),      & ! Effective reflectivity factor (mm^6/m^3)
          kr_vol(:,:,:),     & ! Attenuation coefficient hydro (dB/km)
          g_vol(:,:,:),      & ! Attenuation coefficient gases (dB/km)
          g_to_vol_in(:,:)     ! Gaseous atteunation, radar to vol (dB)
     type(radar_cfg),pointer :: rcfg   ! Radar simulator configuration
  end type cloudsat_IN

CONTAINS

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !  SUBROUTINE cosp_cloudsat_in
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_CLOUDSAT_INIT(radar_freq,k2,use_gas_abs,do_ray,undef,nhydro,   &
                                surface_radar,rcfg,cloudsat_micro_scheme,load_LUT)
    ! INPUTS
    real(wp),intent(in) :: &
         radar_freq,  & ! Radar frequency (GHz)
         k2,          & ! |K|^2, the dielectric constant
         undef          ! Undefined
    integer,intent(in) :: &
         use_gas_abs, & ! 1 = do gaseous abs calcs, 0=no gasesous absorbtion calculated,
                        ! 2 = calculate absorption for first profile on all profiles
         do_ray,      & !
         nhydro,      & !
         surface_radar
    logical,intent(in),optional :: &
         load_LUT
    character(len=64),intent(in) :: &
         cloudsat_micro_scheme

    ! OUTPUTS
    type(radar_cfg) :: &
         rcfg           !

    ! LOCAL VARIABLES
    character(len=240) :: LUT_file_name
    logical       :: local_load_LUT
    integer       :: j

    if (present(load_LUT)) then
       local_load_LUT = load_LUT
    else
       local_load_LUT = RADAR_SIM_LOAD_scale_LUTs_flag
    endif
    
    ! LUT file name
    LUT_file_name = trim(RADAR_SIM_LUT_DIRECTORY) // &
         trim(cloudsat_micro_scheme)

    ! Initialize for NEW radar-configurarion derived type (radar_cfg)
    rcfg%freq                = radar_freq
    rcfg%k2                  = k2
    rcfg%use_gas_abs         = use_gas_abs
    rcfg%do_ray              = do_ray
    rcfg%nhclass             = nhydro
    rcfg%load_scale_LUTs     = local_load_LUT
    rcfg%update_scale_LUTs   = .false.
    rcfg%scale_LUT_file_name = LUT_file_name
    rcfg%N_scale_flag        = .false.
    rcfg%fc                  = undef
    rcfg%rho_eff             = undef
    rcfg%Z_scale_flag        = .false.
    rcfg%Ze_scaled           = 0._wp
    rcfg%Zr_scaled           = 0._wp
    rcfg%kr_scaled           = 0._wp

    ! Set up Re bin "structure" for z_scaling
    rcfg%base_list(1)=0
    do j=1,Re_MAX_BIN
       rcfg%step_list(j)=0.1_wp+0.1_wp*((j-1)**1.5_wp)
       if(rcfg%step_list(j)>Re_BIN_LENGTH) then
          rcfg%step_list(j)=Re_BIN_LENGTH
       endif
       if(j>1) then
          rcfg%base_list(j)=rcfg%base_list(j-1)+floor(Re_BIN_LENGTH/rcfg%step_list(j-1))
       endif
    enddo

    ! Set flag denoting position of radar
    if (surface_radar == 1) then
       rcfg%radar_at_layer_one = .false.
    else
       rcfg%radar_at_layer_one = .true.
    endif

  END SUBROUTINE COSP_CLOUDSAT_INIT

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE COSP_ASSIGN_cloudsatIN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_ASSIGN_cloudsatIN(cospIN, cospgridIN, Npoints, cloudsatIN,                 &
                                CSCAL_MASK_INDICES, CSCAL_SWATH_MASK)
    type(cosp_optical_inputs),intent(in),target :: cospIN     ! Optical inputs to COSP simulator
    type(cosp_column_inputs), intent(in),target :: cospgridIN ! Host model inputs to COSP
    integer,intent(in),target :: &
         Npoints
    type(cloudsat_IN),intent(inout) :: &
         cloudsatIN
    integer,dimension(:),allocatable,intent(out) :: & ! Array containing the indices of the swath masks, already allocated?
         CSCAL_MASK_INDICES
    logical,dimension(:),allocatable,intent(inout) :: & ! Mask of reals over all local times
         CSCAL_SWATH_MASK

    ! Local variables
    integer, target :: &
         N_CLOUDSAT_SWATHED,  &
         i

    if (cospIN % cospswathsIN(3) % N_inst_swaths .gt. 0) then
       if (.not.allocated(CSCAL_SWATH_MASK)) allocate(CSCAL_SWATH_MASK(Npoints))
       ! Do swathing to figure out which cells to simulate on
       call compute_orbitmasks(Npoints,                                                &
                               cospIN % cospswathsIN(3) % N_inst_swaths,               &
                               cospIN % cospswathsIN(3) % inst_localtimes,             &
                               cospIN % cospswathsIN(3) % inst_localtime_widths,       &
                               cospgridIN%lat, cospgridIN%lon,                         &
                               cospgridIN%rttov_date(:,2), cospgridIN%rttov_date(:,3), & ! Time fields: month, dayofmonth
                               cospgridIN%rttov_time(:,1), cospgridIN%rttov_time(:,2), & ! Time fields: hour, minute
                               CSCAL_SWATH_MASK,N_CLOUDSAT_SWATHED) ! Output: logical mask array
       cloudsatIN%Npoints = N_CLOUDSAT_SWATHED
       cloudsatIN%Ncolumns => cospIN%Ncolumns
       if (.not. allocated(CSCAL_MASK_INDICES)) allocate(CSCAL_MASK_INDICES(cloudsatIN%Npoints))
       CSCAL_MASK_INDICES = pack((/ (i, i = 1, Npoints ) /),mask = CSCAL_SWATH_MASK)
       if (cloudsatIN%Npoints .gt. 0) then
          ! Allocate swathed arrays.
          allocate(temp_z_vol_cloudsat(cloudsatIN%Npoints,cospIN%Ncolumns,cospIN%Nlevels),   &
                   temp_kr_vol_cloudsat(cloudsatIN%Npoints,cospIN%Ncolumns,cospIN%Nlevels),  &
                   temp_g_vol_cloudsat(cloudsatIN%Npoints,cospIN%Ncolumns,cospIN%Nlevels),   &
                   temp_hgt_matrix(cloudsatIN%Npoints,cospIN%Nlevels))
          ! Encode step: Read only appropriate values into the new temp arrays. 
          temp_z_vol_cloudsat      = cospIN%z_vol_cloudsat(int(CSCAL_MASK_INDICES),:,:)
          temp_kr_vol_cloudsat     = cospIN%kr_vol_cloudsat(int(CSCAL_MASK_INDICES),:,:)
          temp_g_vol_cloudsat      = cospIN%g_vol_cloudsat(int(CSCAL_MASK_INDICES),:,:)
          temp_hgt_matrix          = cospgridIN%hgt_matrix(int(CSCAL_MASK_INDICES),:)
          ! Reassign swathed values.  
          cloudsatIN%Nlevels    => cospIN%Nlevels
          cloudsatIN%z_vol      => temp_z_vol_cloudsat
          cloudsatIN%kr_vol     => temp_kr_vol_cloudsat
          cloudsatIN%g_vol      => temp_g_vol_cloudsat
          cloudsatIN%rcfg       => cospIN%rcfg_cloudsat
          cloudsatIN%hgt_matrix => temp_hgt_matrix
       end if
    else
       cloudsatIN%Npoints    =  Npoints
       cloudsatIN%Ncolumns   => cospIN%Ncolumns
       cloudsatIN%Nlevels    => cospIN%Nlevels
       cloudsatIN%z_vol      => cospIN%z_vol_cloudsat
       cloudsatIN%kr_vol     => cospIN%kr_vol_cloudsat
       cloudsatIN%g_vol      => cospIN%g_vol_cloudsat
       cloudsatIN%rcfg       => cospIN%rcfg_cloudsat
       cloudsatIN%hgt_matrix => cospgridIN%hgt_matrix    
    end if

  END SUBROUTINE COSP_ASSIGN_cloudsatIN

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE COSP_ASSIGN_cloudsatIN_CLEAN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_ASSIGN_cloudsatIN_CLEAN()
    ! Deallocate temporary arrays
    if (allocated(temp_z_vol_cloudsat))     deallocate(temp_z_vol_cloudsat)
    if (allocated(temp_kr_vol_cloudsat))    deallocate(temp_kr_vol_cloudsat)
    if (allocated(temp_g_vol_cloudsat))     deallocate(temp_g_vol_cloudsat)
    if (allocated(temp_hgt_matrix))         deallocate(temp_hgt_matrix)

  END SUBROUTINE COSP_ASSIGN_cloudsatIN_CLEAN

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! END MODULE
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE MOD_COSP_CLOUDSAT_INTERFACE

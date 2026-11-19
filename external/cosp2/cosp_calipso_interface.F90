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
! Jul 2017 - R. Guzman - Added Ground LIDar variables (GLID)
! Jun 2025 - J.K. Shaw - Added COSP-RTTOV integration and swathing
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE MOD_COSP_CALIPSO_INTERFACE
  USE COSP_KINDS,              ONLY: wp
  USE MOD_LIDAR_SIMULATOR,     ONLY: alpha,beta,gamma
  use mod_cosp_stats,          ONLY: compute_orbitmasks,cosp_optical_inputs,cosp_column_inputs
  IMPLICIT NONE

  ! Module variables
  real(wp),dimension(:,:),target,allocatable :: &
        temp_beta_mol_calipso,    &
        temp_tau_mol_calipso
  real(wp),dimension(:,:,:),target,allocatable :: &
        temp_betatot_calipso,         &
        temp_tautot_calipso,          &
        temp_betatot_liq_calipso,     &
        temp_tautot_liq_calipso,      &
        temp_betatot_ice_calipso,     &
        temp_tautot_ice_calipso

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! TYPE calipso_in
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  type calipso_IN
     integer         ::       &
          Npoints         ! Number of gridpoints.
     integer,pointer ::       &
          Ncolumns,     & ! Number of columns.
          Nlevels         ! Number of levels.
     real(wp),dimension(:,:),pointer :: &
          beta_mol,     & ! Molecular backscatter coefficient
          tau_mol         ! Molecular optical depth
     real(wp),dimension(:,:,:),pointer :: &
          betatot,      & ! 
          tautot,       & ! Optical thickess integrated from top
          betatot_ice,  & ! Backscatter coefficient for ice particles
          betatot_liq,  & ! Backscatter coefficient for liquid particles
          tautot_ice,   & ! Total optical thickness of ice
          tautot_liq      ! Total optical thickness of liq
     real(wp),dimension(:,:,:,:),pointer :: &
          taupart
  end type calipso_IN

CONTAINS
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE cosp_calipso_init
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine cosp_calipso_init() 
    
    ! Polynomial coefficients (Alpha, Beta, Gamma) which allow to compute the 
    ! ATBperpendicular as a function of the ATB for ice or liquid cloud particles 
    ! derived from CALIPSO-GOCCP observations at 120m vertical grid 
    ! (Cesana and Chepfer, JGR, 2013).
    !
    ! Relationship between ATBice and ATBperp,ice for ice particles:
    !                ATBperp,ice = Alpha*ATBice 
    ! Relationship between ATBice and ATBperp,ice for liquid particles:
    !          ATBperp,ice = Beta*ATBice^2 + Gamma*ATBice
    Alpha = 0.2904_wp
    Beta  = 0.4099_wp
    Gamma = 0.009_wp    
    
  end subroutine cosp_calipso_init

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !  							SUBROUTINE COSP_ASSIGN_calipsoIN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_ASSIGN_calipsoIN(cospIN,cospgridIN,Npoints,calipsoIN,CSCAL_MASK_INDICES,CSCAL_SWATH_MASK)
     type(cosp_optical_inputs),intent(in),target :: cospIN     ! Optical inputs to COSP simulator
     type(cosp_column_inputs), intent(in),target :: cospgridIN ! Host model inputs to COSP
     integer,intent(in),target :: &
         Npoints
     type(calipso_IN),intent(inout) :: &
         calipsoIN
     integer,dimension(:),allocatable,intent(out) :: & ! Array containing the indices of the swath masks, already allocated?
         CSCAL_MASK_INDICES
     logical,dimension(:),allocatable,intent(inout) :: & ! Mask of reals over all local times
         CSCAL_SWATH_MASK

     ! Local variables
     integer, target :: &
         N_CALIPSO_SWATHED,  &
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
                                 CSCAL_SWATH_MASK,N_CALIPSO_SWATHED) ! Output: logical mask array
         calipsoIN%Npoints = N_CALIPSO_SWATHED
         calipsoIN%Ncolumns => cospIN%Ncolumns
         allocate(CSCAL_MASK_INDICES(N_CALIPSO_SWATHED))
         CSCAL_MASK_INDICES = pack((/ (i, i = 1, Npoints ) /),mask = CSCAL_SWATH_MASK)
         if (calipsoIN%Npoints .gt. 0) then
             ! Allocate swathed arrays.
             allocate(temp_beta_mol_calipso(calipsoIN%Npoints,cospIN%Nlevels),                        &
                      temp_tau_mol_calipso(calipsoIN%Npoints,cospIN%Nlevels),                         &
                      temp_betatot_calipso(calipsoIN%Npoints,cospIN%Ncolumns,cospIN%Nlevels),         &
                      temp_tautot_calipso(calipsoIN%Npoints,cospIN%Ncolumns,cospIN%Nlevels),          &
                      temp_betatot_liq_calipso(calipsoIN%Npoints,cospIN%Ncolumns,cospIN%Nlevels),     &
                      temp_tautot_liq_calipso(calipsoIN%Npoints,cospIN%Ncolumns,cospIN%Nlevels),      &
                      temp_betatot_ice_calipso(calipsoIN%Npoints,cospIN%Ncolumns,cospIN%Nlevels),     &
                      temp_tautot_ice_calipso(calipsoIN%Npoints,cospIN%Ncolumns,cospIN%Nlevels))
             ! Encode step: Read only appropriate values into the new temp arrays. 
             temp_beta_mol_calipso(:,:)         = cospIN%beta_mol_calipso(int(CSCAL_MASK_INDICES),:)
             temp_tau_mol_calipso(:,:)          = cospIN%tau_mol_calipso(int(CSCAL_MASK_INDICES),:)
             temp_betatot_calipso(:,:,:)        = cospIN%betatot_calipso(int(CSCAL_MASK_INDICES),:,:)
             temp_tautot_calipso(:,:,:)         = cospIN%tautot_calipso(int(CSCAL_MASK_INDICES),:,:)
             temp_betatot_liq_calipso(:,:,:)    = cospIN%betatot_liq_calipso(int(CSCAL_MASK_INDICES),:,:)
             temp_tautot_liq_calipso(:,:,:)     = cospIN%tautot_liq_calipso(int(CSCAL_MASK_INDICES),:,:)
             temp_betatot_ice_calipso(:,:,:)    = cospIN%betatot_ice_calipso(int(CSCAL_MASK_INDICES),:,:)
             temp_tautot_ice_calipso(:,:,:)     = cospIN%tautot_ice_calipso(int(CSCAL_MASK_INDICES),:,:)

             calipsoIN%Nlevels     => cospIN%Nlevels
             calipsoIN%beta_mol    => temp_beta_mol_calipso
             calipsoIN%betatot     => temp_betatot_calipso
             calipsoIN%betatot_liq => temp_betatot_liq_calipso
             calipsoIN%betatot_ice => temp_betatot_ice_calipso
             calipsoIN%tau_mol     => temp_tau_mol_calipso
             calipsoIN%tautot      => temp_tautot_calipso
             calipsoIN%tautot_liq  => temp_tautot_liq_calipso
             calipsoIN%tautot_ice  => temp_tautot_ice_calipso
         end if
     else
          calipsoIN%Npoints     =  Npoints
          calipsoIN%Ncolumns    => cospIN%Ncolumns
          calipsoIN%Nlevels     => cospIN%Nlevels
          calipsoIN%beta_mol    => cospIN%beta_mol_calipso
          calipsoIN%betatot     => cospIN%betatot_calipso
          calipsoIN%betatot_liq => cospIN%betatot_liq_calipso
          calipsoIN%betatot_ice => cospIN%betatot_ice_calipso
          calipsoIN%tau_mol     => cospIN%tau_mol_calipso
          calipsoIN%tautot      => cospIN%tautot_calipso
          calipsoIN%tautot_liq  => cospIN%tautot_liq_calipso
          calipsoIN%tautot_ice  => cospIN%tautot_ice_calipso 
     end if

  END SUBROUTINE COSP_ASSIGN_calipsoIN

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !  							SUBROUTINE COSP_ASSIGN_calipsoIN_CLEAN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_ASSIGN_calipsoIN_CLEAN()
    ! Deallocate temporary arrays
    if (allocated(temp_beta_mol_calipso))       deallocate(temp_beta_mol_calipso)
    if (allocated(temp_tau_mol_calipso))        deallocate(temp_tau_mol_calipso)
    if (allocated(temp_betatot_calipso))        deallocate(temp_betatot_calipso)
    if (allocated(temp_tautot_calipso))         deallocate(temp_tautot_calipso)
    if (allocated(temp_betatot_liq_calipso))    deallocate(temp_betatot_liq_calipso)
    if (allocated(temp_tautot_liq_calipso))     deallocate(temp_tautot_liq_calipso)
    if (allocated(temp_betatot_ice_calipso))    deallocate(temp_betatot_ice_calipso)
    if (allocated(temp_tautot_ice_calipso))     deallocate(temp_tautot_ice_calipso)

  END SUBROUTINE COSP_ASSIGN_calipsoIN_CLEAN

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !	END MODULE
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE MOD_COSP_CALIPSO_INTERFACE

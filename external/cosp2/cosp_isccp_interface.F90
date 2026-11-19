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
MODULE MOD_COSP_ISCCP_INTERFACE
  USE COSP_KINDS,      ONLY: wp
  USE mod_icarus,      ONLY: isccp_top_height,isccp_top_height_direction
  use mod_cosp_stats,  ONLY: compute_orbitmasks,cosp_optical_inputs,cosp_column_inputs

  IMPLICIT NONE
  ! Module variables
  integer,dimension(:),target,allocatable :: &
    temp_isccp_sunlit
  real(wp),dimension(:),target,allocatable :: &
    temp_isccp_skt
  real(wp),dimension(:,:),target,allocatable :: &
    temp_isccp_qv,       &         
    temp_isccp_at,       &
    temp_isccp_phalf,    &
    temp_isccp_pfull   
  real(wp),dimension(:,:,:),target,allocatable :: &
    temp_isccp_frac_out, &
    temp_isccp_tau_067,  &
    temp_isccp_emiss_11

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !                                  TYPE isccp_in
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! Derived input type for ISCCP simulator
  type isccp_IN
     integer          ::       &
          Npoints                ! Number of gridpoints.
     integer,pointer  ::       &
          Ncolumns,            & ! Number of columns.
          Nlevels,             & ! Number of levels.
          top_height,          & !
          top_height_direction   !
     integer,pointer ::        &
          sunlit(:)              ! Sunlit points (npoints)
     real(wp),pointer ::       &
          emsfc_lw
     real(wp),pointer ::       &
          skt(:)                 ! Surface temperature (npoints)
     real(wp),pointer ::       &
          at(:,:),             & ! Temperature (npoint,nlev)
          pfull(:,:),          & ! Pressure (npoints,nlev)
          qv(:,:)                ! Specific humidity (npoints,nlev)
     real(wp),pointer  ::       &          
          phalf(:,:)             ! Pressure at half levels (npoints,nlev+1)
     real(wp),pointer ::       &
          frac_out(:,:,:),     & ! Cloud fraction (npoints,ncolumns,nlevels)
          dtau(:,:,:),         & ! Optical depth (npoints,ncolumns,nlevels)
          dem(:,:,:)             ! Emissivity (npoints,ncolumns,nlevels)
  end type isccp_IN

CONTAINS

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !  							SUBROUTINE cosp_isccp_init
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_ISCCP_INIT(top_height,top_height_direction)
     integer,intent(in) :: &
         top_height, &
         top_height_direction

    ! Cloud-top height determination
    isccp_top_height           = top_height
    isccp_top_height_direction = top_height_direction

  END SUBROUTINE COSP_ISCCP_INIT

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !  							SUBROUTINE COSP_ASSIGN_isccpIN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_ASSIGN_isccpIN(cospIN,cospgridIN,Npoints,isccpIN,ISCCP_MASK_INDICES)
     type(cosp_optical_inputs),intent(in),target :: cospIN     ! Optical inputs to COSP simulator
     type(cosp_column_inputs), intent(in),target :: cospgridIN ! Host model inputs to COSP
     integer,intent(in),target :: &
         Npoints
     type(isccp_IN),intent(inout) :: &
         isccpIN
     integer,dimension(:),allocatable,intent(out) :: & ! Array containing the indices of the swath masks, already allocated?
         ISCCP_MASK_INDICES

     ! Local variables
     logical,dimension(:),allocatable :: & ! Mask of reals over all local times
         ISCCP_SWATH_MASK
     integer, target :: &
         N_ISCCP_SWATHED,  &
         i

     if (cospIN % cospswathsIN(1) % N_inst_swaths .gt. 0) then
         allocate(ISCCP_SWATH_MASK(Npoints))
         ! Do swathing to figure out which cells to simulate on
         call compute_orbitmasks(Npoints,                                                &
                                 cospIN % cospswathsIN(1) % N_inst_swaths,               &
                                 cospIN % cospswathsIN(1) % inst_localtimes,             &
                                 cospIN % cospswathsIN(1) % inst_localtime_widths,       &
                                 cospgridIN%lat, cospgridIN%lon,                         &
                                 cospgridIN%rttov_date(:,2), cospgridIN%rttov_date(:,3), & ! Time fields: month, dayofmonth
                                 cospgridIN%rttov_time(:,1), cospgridIN%rttov_time(:,2), & ! Time fields: hour, minute
                                 ISCCP_SWATH_MASK,N_ISCCP_SWATHED) ! Output: logical mask array
         isccpIN%Npoints = N_ISCCP_SWATHED
         isccpIN%Ncolumns => cospIN%Ncolumns
         allocate(ISCCP_MASK_INDICES(N_ISCCP_SWATHED))
         ISCCP_MASK_INDICES = pack((/ (i, i = 1, Npoints ) /),mask = ISCCP_SWATH_MASK)
         deallocate(ISCCP_SWATH_MASK)
         if (isccpIN%Npoints .gt. 0) then
             ! Allocate swathed arrays.
             allocate(temp_isccp_skt(isccpIN%Npoints),temp_isccp_qv(isccpIN%Npoints,cospIN%Nlevels),temp_isccp_at(isccpIN%Npoints,cospIN%Nlevels),              &
                      temp_isccp_frac_out(isccpIN%Npoints,cospIN%Ncolumns,cospIN%Nlevels),temp_isccp_tau_067(isccpIN%Npoints,cospIN%Ncolumns,cospIN%Nlevels),   &
                      temp_isccp_emiss_11(isccpIN%Npoints,cospIN%Ncolumns,cospIN%Nlevels),temp_isccp_phalf(isccpIN%Npoints,cospIN%Nlevels+1),                   &
                      temp_isccp_pfull(isccpIN%Npoints,cospIN%Nlevels),temp_isccp_sunlit(isccpIN%Npoints))
             ! Encode step: Read only appropriate values into the new temp arrays. 
             temp_isccp_skt(:)          = cospgridIN%skt(int(ISCCP_MASK_INDICES))
             temp_isccp_qv(:,:)         = cospgridIN%qv(int(ISCCP_MASK_INDICES),:)
             temp_isccp_at(:,:)         = cospgridIN%at(int(ISCCP_MASK_INDICES),:)
             temp_isccp_frac_out(:,:,:) = cospIN%frac_out(int(ISCCP_MASK_INDICES),:,:)
             temp_isccp_tau_067(:,:,:)  = cospIN%tau_067(int(ISCCP_MASK_INDICES),:,:)
             temp_isccp_emiss_11(:,:,:) = cospIN%emiss_11(int(ISCCP_MASK_INDICES),:,:)
             temp_isccp_phalf(:,:)      = cospgridIN%phalf(int(ISCCP_MASK_INDICES),:)
             temp_isccp_pfull(:,:)      = cospgridIN%pfull(int(ISCCP_MASK_INDICES),:)
             temp_isccp_sunlit(:)       = cospgridIN%sunlit(int(ISCCP_MASK_INDICES))

             isccpIN%Nlevels  => cospIN%Nlevels
             isccpIN%emsfc_lw => cospIN%emsfc_lw
             isccpIN%skt      => temp_isccp_skt
             isccpIN%qv       => temp_isccp_qv
             isccpIN%at       => temp_isccp_at
             isccpIN%frac_out => temp_isccp_frac_out
             isccpIN%dtau     => temp_isccp_tau_067
             isccpIN%dem      => temp_isccp_emiss_11
             isccpIN%phalf    => temp_isccp_phalf
             isccpIN%pfull    => temp_isccp_pfull
             isccpIN%sunlit   => temp_isccp_sunlit
         end if
     else
         isccpIN%Npoints  = Npoints
         isccpIN%Ncolumns => cospIN%Ncolumns
         isccpIN%Nlevels  => cospIN%Nlevels
         isccpIN%emsfc_lw => cospIN%emsfc_lw
         isccpIN%skt      => cospgridIN%skt
         isccpIN%qv       => cospgridIN%qv
         isccpIN%at       => cospgridIN%at
         isccpIN%frac_out => cospIN%frac_out
         isccpIN%dtau     => cospIN%tau_067
         isccpIN%dem      => cospIN%emiss_11
         isccpIN%phalf    => cospgridIN%phalf
         isccpIN%pfull    => cospgridIN%pfull
         isccpIN%sunlit   => cospgridIN%sunlit
     end if

  END SUBROUTINE COSP_ASSIGN_isccpIN

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !  							SUBROUTINE COSP_ASSIGN_isccpIN_CLEAN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_ASSIGN_isccpIN_CLEAN()
    ! Deallocate temporary arrays
    if (allocated(temp_isccp_sunlit))   deallocate(temp_isccp_sunlit)
    if (allocated(temp_isccp_skt))      deallocate(temp_isccp_skt)
    if (allocated(temp_isccp_qv))       deallocate(temp_isccp_qv)
    if (allocated(temp_isccp_at))       deallocate(temp_isccp_at)
    if (allocated(temp_isccp_phalf))    deallocate(temp_isccp_phalf)
    if (allocated(temp_isccp_pfull))    deallocate(temp_isccp_pfull)
    if (allocated(temp_isccp_frac_out)) deallocate(temp_isccp_frac_out)
    if (allocated(temp_isccp_tau_067))  deallocate(temp_isccp_tau_067)
    if (allocated(temp_isccp_emiss_11)) deallocate(temp_isccp_emiss_11)

  END SUBROUTINE COSP_ASSIGN_isccpIN_CLEAN

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !                                    END MODULE
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE MOD_COSP_ISCCP_INTERFACE

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
MODULE MOD_COSP_MISR_INTERFACE
  USE COSP_KINDS,  ONLY: wp
  use mod_cosp_stats,  ONLY: compute_orbitmasks,cosp_optical_inputs,cosp_column_inputs

  IMPLICIT NONE

  ! Module variables
  integer,dimension(:),target,allocatable :: &
      temp_misr_sunlit
  real(wp),dimension(:,:),target,allocatable :: &
      temp_misr_zfull,   &
      temp_misr_at
  real(wp),dimension(:,:,:),target,allocatable :: &
      temp_misr_dtau

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! 								TYPE misr_in
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  type misr_IN
     integer          ::  &
          Npoints           ! Number of gridpoints.
     integer,pointer  ::  &
          Ncolumns,       & ! Number of columns.
          Nlevels           ! Number of levels.
     integer,pointer ::   &
          sunlit(:)         ! Sunlit points (npoints).
     real(wp),pointer ::  &
          zfull(:,:),     & ! Height of full model levels (i.e. midpoints). (npoints,nlev)
          at(:,:)           ! Temperature. (npoints,nlev)
     real(wp),pointer ::  &                           
          dtau(:,:,:)       ! Optical depth. (npoints,ncolumns,nlev)
          
  end type misr_IN

CONTAINS

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !  							SUBROUTINE cosp_misr_init
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_MISR_INIT()

  END SUBROUTINE COSP_MISR_INIT

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !  							SUBROUTINE COSP_ASSIGN_misrIN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_ASSIGN_misrIN(cospIN,cospgridIN,Npoints,misrIN,MISR_MASK_INDICES)
     type(cosp_optical_inputs),intent(in),target :: cospIN     ! Optical inputs to COSP simulator
     type(cosp_column_inputs), intent(in),target :: cospgridIN ! Host model inputs to COSP
     integer,intent(in),target :: &
         Npoints
     type(misr_IN),intent(inout) :: &
         misrIN
     integer,dimension(:),allocatable,intent(out) :: & ! Array containing the indices of the swath masks, already allocated?
         MISR_MASK_INDICES

     ! Local variables
     logical,dimension(:),allocatable :: & ! Mask of reals over all local times
         MISR_SWATH_MASK
     integer, target :: &
         N_MISR_SWATHED,  &
         i

     if (cospIN % cospswathsIN(2) % N_inst_swaths .gt. 0) then
         allocate(MISR_SWATH_MASK(Npoints))
         ! Do swathing to figure out which cells to simulate on
         call compute_orbitmasks(Npoints,                                                &
                                 cospIN % cospswathsIN(2) % N_inst_swaths,               &
                                 cospIN % cospswathsIN(2) % inst_localtimes,             &
                                 cospIN % cospswathsIN(2) % inst_localtime_widths,       &
                                 cospgridIN%lat, cospgridIN%lon,                         &
                                 cospgridIN%rttov_date(:,2), cospgridIN%rttov_date(:,3), & ! Time fields: month, dayofmonth
                                 cospgridIN%rttov_time(:,1), cospgridIN%rttov_time(:,2), & ! Time fields: hour, minute
                                 MISR_SWATH_MASK,N_MISR_SWATHED) ! Output: logical mask array
         misrIN%Npoints = N_MISR_SWATHED
         misrIN%Ncolumns => cospIN%Ncolumns
         allocate(MISR_MASK_INDICES(N_MISR_SWATHED))
         MISR_MASK_INDICES = pack((/ (i, i = 1, Npoints ) /),mask = MISR_SWATH_MASK)
         deallocate(MISR_SWATH_MASK)
         if (misrIN%Npoints .gt. 0) then
             ! Allocate swathed arrays.
             allocate(temp_misr_dtau(misrIN%Npoints,cospIN%Ncolumns,cospIN%Nlevels),   &
                      temp_misr_sunlit(misrIN%Npoints),                                &
                      temp_misr_zfull(misrIN%Npoints,cospIN%Nlevels),                  &
                      temp_misr_at(misrIN%Npoints,cospIN%Nlevels))
             ! Encode step: Read only appropriate values into the new temp arrays. 
             temp_misr_dtau(:,:,:)     = cospIN%tau_067(int(MISR_MASK_INDICES),:,:)
             temp_misr_at(:,:)         = cospgridIN%at(int(MISR_MASK_INDICES),:)
             temp_misr_zfull(:,:)      = cospgridIN%hgt_matrix(int(MISR_MASK_INDICES),:)
             temp_misr_sunlit(:)       = cospgridIN%sunlit(int(MISR_MASK_INDICES))

             misrIN%Nlevels  => cospIN%Nlevels
             misrIN%dtau     => temp_misr_dtau
             misrIN%sunlit   => temp_misr_sunlit
             misrIN%zfull    => temp_misr_zfull
             misrIN%at       => temp_misr_at
         end if
     else
          misrIN%Npoints  =  Npoints
          misrIN%Ncolumns => cospIN%Ncolumns
          misrIN%Nlevels  => cospIN%Nlevels
          misrIN%dtau     => cospIN%tau_067
          misrIN%sunlit   => cospgridIN%sunlit
          misrIN%zfull    => cospgridIN%hgt_matrix
          misrIN%at       => cospgridIN%at
     end if

  END SUBROUTINE COSP_ASSIGN_misrIN

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !  							SUBROUTINE COSP_ASSIGN_misrIN_CLEAN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_ASSIGN_misrIN_CLEAN()
    ! Deallocate temporary arrays
    if (allocated(temp_misr_sunlit))  deallocate(temp_misr_sunlit)
    if (allocated(temp_misr_zfull))   deallocate(temp_misr_zfull)
    if (allocated(temp_misr_at))      deallocate(temp_misr_at)
    if (allocated(temp_misr_dtau))    deallocate(temp_misr_dtau)

  END SUBROUTINE COSP_ASSIGN_misrIN_CLEAN

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! 							    	END MODULE
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE MOD_COSP_MISR_INTERFACE

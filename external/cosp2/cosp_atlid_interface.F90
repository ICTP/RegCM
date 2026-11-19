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
! Apr 2018 - R. Guzman - Original version
! Jun 2025 - J.K. Shaw - Added COSP-RTTOV integration and swathing
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE MOD_COSP_ATLID_INTERFACE
  USE COSP_KINDS,              ONLY: wp
  use mod_cosp_stats,  ONLY: compute_orbitmasks,cosp_optical_inputs,cosp_column_inputs
  IMPLICIT NONE

  ! Module variables
  real(wp),dimension(:,:),target,allocatable :: &
        temp_beta_mol_atlid,          &
        temp_tau_mol_atlid
  real(wp),dimension(:,:,:),target,allocatable :: &
        temp_betatot_atlid,           &
        temp_tautot_atlid

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! TYPE atlid_in
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
  type atlid_IN
     integer         ::       &
          Npoints         ! Number of gridpoints.
     integer,pointer ::       &
          Ncolumns,     & ! Number of columns.
          Nlevels         ! Number of levels.

     real(wp),dimension(:,:),pointer :: &
          beta_mol_atlid,     & ! Molecular backscatter coefficient
          tau_mol_atlid         ! Molecular optical depth
     real(wp),dimension(:,:,:),pointer :: &
          betatot_atlid,      & ! 
          tautot_atlid          ! Optical thickess integrated from top
     real(wp),dimension(:,:,:,:),pointer :: &
          taupart_atlid
  end type atlid_IN

CONTAINS
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE cosp_atlid_init
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  subroutine cosp_atlid_init() 

  end subroutine cosp_atlid_init

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !  							SUBROUTINE COSP_ASSIGN_atlidIN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_ASSIGN_atlidIN(cospIN,cospgridIN,Npoints,atlidIN,ATLID_MASK_INDICES)
     type(cosp_optical_inputs),intent(in),target :: cospIN     ! Optical inputs to COSP simulator
     type(cosp_column_inputs), intent(in),target :: cospgridIN ! Host model inputs to COSP
     integer,intent(in),target :: &
         Npoints
     type(atlid_IN),intent(inout) :: &
         atlidIN
     integer,dimension(:),allocatable,intent(out) :: & ! Array containing the indices of the swath masks, already allocated?
         ATLID_MASK_INDICES

     ! Local variables
     logical,dimension(:),allocatable :: & ! Mask of reals over all local times
         ATLID_SWATH_MASK
     integer, target :: &
         N_ATLID_SWATHED,  &
         i

     if (cospIN % cospswathsIN(4) % N_inst_swaths .gt. 0) then
         allocate(ATLID_SWATH_MASK(Npoints))
         ! Do swathing to figure out which cells to simulate on
         call compute_orbitmasks(Npoints,                                                &
                                 cospIN % cospswathsIN(4) % N_inst_swaths,               &
                                 cospIN % cospswathsIN(4) % inst_localtimes,             &
                                 cospIN % cospswathsIN(4) % inst_localtime_widths,       &
                                 cospgridIN%lat, cospgridIN%lon,                         &
                                 cospgridIN%rttov_date(:,2), cospgridIN%rttov_date(:,3), & ! Time fields: month, dayofmonth
                                 cospgridIN%rttov_time(:,1), cospgridIN%rttov_time(:,2), & ! Time fields: hour, minute
                                 ATLID_SWATH_MASK,N_ATLID_SWATHED) ! Output: logical mask array
         atlidIN%Npoints = N_ATLID_SWATHED
         atlidIN%Ncolumns => cospIN%Ncolumns
         allocate(ATLID_MASK_INDICES(N_ATLID_SWATHED))
         ATLID_MASK_INDICES = pack((/ (i, i = 1, Npoints ) /),mask = ATLID_SWATH_MASK)
         deallocate(ATLID_SWATH_MASK)
         if (atlidIN%Npoints .gt. 0) then
             ! Allocate swathed arrays.
             allocate(temp_beta_mol_atlid(atlidIN%Npoints,cospIN%Nlevels),                  &
                      temp_betatot_atlid(atlidIN%Npoints,cospIN%Ncolumns,cospIN%Nlevels),   &
                      temp_tau_mol_atlid(atlidIN%Npoints,cospIN%Nlevels),                   &
                      temp_tautot_atlid(atlidIN%Npoints,cospIN%Ncolumns,cospIN%Nlevels))
             ! Encode step: Read only appropriate values into the new temp arrays. 
             temp_beta_mol_atlid(:,:)          = cospIN%beta_mol_atlid(int(ATLID_MASK_INDICES),:)
             temp_tau_mol_atlid(:,:)           = cospIN%tau_mol_atlid(int(ATLID_MASK_INDICES),:)
             temp_betatot_atlid(:,:,:)         = cospIN%betatot_atlid(int(ATLID_MASK_INDICES),:,:)
             temp_tautot_atlid(:,:,:)          = cospIN%tautot_atlid(int(ATLID_MASK_INDICES),:,:)

             atlidIN%Nlevels        => cospIN%Nlevels
             atlidIN%beta_mol_atlid => temp_beta_mol_atlid
             atlidIN%betatot_atlid  => temp_betatot_atlid 
             atlidIN%tau_mol_atlid  => temp_tau_mol_atlid 
             atlidIN%tautot_atlid   => temp_tautot_atlid    
         end if
     else
          atlidIN%Npoints        =  Npoints
          atlidIN%Ncolumns       => cospIN%Ncolumns 
          atlidIN%Nlevels        => cospIN%Nlevels 
          atlidIN%beta_mol_atlid => cospIN%beta_mol_atlid
          atlidIN%betatot_atlid  => cospIN%betatot_atlid 
          atlidIN%tau_mol_atlid  => cospIN%tau_mol_atlid 
          atlidIN%tautot_atlid   => cospIN%tautot_atlid     
     end if

  END SUBROUTINE COSP_ASSIGN_atlidIN

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !  							SUBROUTINE COSP_ASSIGN_atlidIN_CLEAN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_ASSIGN_atlidIN_CLEAN()
    ! Deallocate temporary arrays
    if (allocated(temp_beta_mol_atlid))   deallocate(temp_beta_mol_atlid)
    if (allocated(temp_tau_mol_atlid))    deallocate(temp_tau_mol_atlid)
    if (allocated(temp_betatot_atlid))    deallocate(temp_betatot_atlid)
    if (allocated(temp_tautot_atlid))     deallocate(temp_tautot_atlid)

  END SUBROUTINE COSP_ASSIGN_atlidIN_CLEAN

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !	END MODULE
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE MOD_COSP_ATLID_INTERFACE

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
MODULE MOD_COSP_PARASOL_INTERFACE
  USE COSP_KINDS,  ONLY: WP
  use mod_cosp_stats,  ONLY: compute_orbitmasks,cosp_optical_inputs,cosp_column_inputs
  implicit none

  ! Module variables
  real(wp),dimension(:,:),target,allocatable :: &
      temp_tautot_S_liq,             &
      temp_tautot_S_ice

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !									TYPE cosp_parasol
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  TYPE PARASOL_SGX
     ! Dimensions
     integer :: &
          Npoints,  & ! Number of gridpoints
          Ncolumns, & ! Number of columns
          Nrefl       ! Number of parasol reflectances
     
     ! Arrays with dimensions (Npoints,Ncolumns,Nrefl)
     real(wp),dimension(:,:,:),pointer :: &
          refl        ! parasol reflectances

  END TYPE PARASOL_SGX
  TYPE PARASOL_GBX
     integer :: &
          Npoints,  & ! Number of gridpoints
          Ncolumns, & ! Number of columns
          Nrefl       ! Number of parasol reflectances
     real(wp), dimension(:,:),pointer :: &
          parasolrefl ! Mean parasol reflectance

  END TYPE PARASOL_GBX
  TYPE COSP_PARASOL 
     type(PARASOL_SGX) :: PARASOL_SGX
     type(PARASOL_GBX) :: PARASOL_GBX
  END TYPE COSP_PARASOL
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !										TYPE parasol_in
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  TYPE parasol_IN
     integer         :: &
        Npoints          ! Number of horizontal gridpoints
     integer,pointer :: &
        Nlevels,       & ! Number of vertical levels
        Ncolumns,      & ! Number of columns
        Nrefl            ! Number of angles for which the reflectance is computed
     real(wp),dimension(:,:),pointer ::   &
        tautot_S_liq,  & ! Liquid water optical thickness, from TOA to SFC
        tautot_S_ice     ! Ice water optical thickness, from TOA to SFC
  END TYPE parasol_IN
  
contains

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !                           SUBROUTINE cosp_parasol_init
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_PARASOL_INIT()
    
  end subroutine COSP_PARASOL_INIT

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !  							SUBROUTINE COSP_ASSIGN_parasolIN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_ASSIGN_parasolIN(cospIN,cospgridIN,Npoints,parasolIN,PARASOL_MASK_INDICES)
     type(cosp_optical_inputs),intent(in),target :: cospIN     ! Optical inputs to COSP simulator
     type(cosp_column_inputs), intent(in),target :: cospgridIN ! Host model inputs to COSP
     integer,intent(in),target :: &
         Npoints
     type(parasol_IN),intent(inout) :: &
         parasolIN
     integer,dimension(:),allocatable,intent(out) :: & ! Array containing the indices of the swath masks, already allocated?
         PARASOL_MASK_INDICES

     ! Local variables
     logical,dimension(:),allocatable :: & ! Mask of reals over all local times
         PARASOL_SWATH_MASK
     integer, target :: &
         N_PARASOL_SWATHED,  &
         i

     if (cospIN % cospswathsIN(5) % N_inst_swaths .gt. 0) then
         allocate(PARASOL_SWATH_MASK(Npoints))
         ! Do swathing to figure out which cells to simulate on
         call compute_orbitmasks(Npoints,                                                &
                                 cospIN % cospswathsIN(5) % N_inst_swaths,               &
                                 cospIN % cospswathsIN(5) % inst_localtimes,             &
                                 cospIN % cospswathsIN(5) % inst_localtime_widths,       &
                                 cospgridIN%lat, cospgridIN%lon,                         &
                                 cospgridIN%rttov_date(:,2), cospgridIN%rttov_date(:,3), & ! Time fields: month, dayofmonth
                                 cospgridIN%rttov_time(:,1), cospgridIN%rttov_time(:,2), & ! Time fields: hour, minute
                                 PARASOL_SWATH_MASK,N_PARASOL_SWATHED) ! Output: logical mask array
         parasolIN%Npoints = N_PARASOL_SWATHED
         parasolIN%Ncolumns => cospIN%Ncolumns
         allocate(PARASOL_MASK_INDICES(N_PARASOL_SWATHED))
         PARASOL_MASK_INDICES = pack((/ (i, i = 1, Npoints ) /),mask = PARASOL_SWATH_MASK)
         deallocate(PARASOL_SWATH_MASK)
         if (parasolIN%Npoints .gt. 0) then
             ! Allocate swathed arrays.
             allocate(temp_tautot_S_liq(parasolIN%Npoints,cospIN%Ncolumns),                  &
                      temp_tautot_S_ice(parasolIN%Npoints,cospIN%Ncolumns))
             ! Encode step: Read only appropriate values into the new temp arrays.
             temp_tautot_S_liq(:,:)          = cospIN%tautot_S_liq(int(PARASOL_MASK_INDICES),:)
             temp_tautot_S_ice(:,:)          = cospIN%tautot_S_ice(int(PARASOL_MASK_INDICES),:)

             parasolIN%Nlevels      => cospIN%Nlevels
             parasolIN%Nrefl        => cospIN%Nrefl
             parasolIN%tautot_S_liq => temp_tautot_S_liq
             parasolIN%tautot_S_ice => temp_tautot_S_ice
         endif
     else  
          parasolIN%Npoints      =  Npoints
          parasolIN%Ncolumns => cospIN%Ncolumns
          parasolIN%Nlevels      => cospIN%Nlevels
          parasolIN%Nrefl        => cospIN%Nrefl
          parasolIN%tautot_S_liq => cospIN%tautot_S_liq
          parasolIN%tautot_S_ice => cospIN%tautot_S_ice
     end if

  END SUBROUTINE COSP_ASSIGN_parasolIN

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  !  							SUBROUTINE COSP_ASSIGN_parasolIN_CLEAN
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_ASSIGN_parasolIN_CLEAN()
    ! Deallocate temporary arrays
    if (allocated(temp_tautot_S_liq))   deallocate(temp_tautot_S_liq)
    if (allocated(temp_tautot_S_ice))    deallocate(temp_tautot_S_ice)

  END SUBROUTINE COSP_ASSIGN_parasolIN_CLEAN

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! 								    END MODULE
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end module MOD_COSP_PARASOL_INTERFACE

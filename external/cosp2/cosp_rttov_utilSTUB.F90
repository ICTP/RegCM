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
! Jun 2025 - J.K. Shaw - Initial version
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

MODULE MOD_COSP_RTTOV_UTIL
  USE COSP_KINDS, ONLY: wp
  IMPLICIT NONE


  ! DDT for each instrument being simulated. Values to be assigned during the cosp_rttov_init subroutine
  type rttov_cfg
      logical                      :: &
          Lrttov_bt,           &
          Lrttov_rad,          &
          Lrttov_refl,         &
          Lrttov_cld,          &
          Lrttov_aer,          &
          Lrttov_pc
      character(len=256)           :: &
          rttov_srcDir,        &
          rttov_coefDir,       &
          OD_coef_filepath,    &
          aer_coef_filepath,   &
          cld_coef_filepath,   &
          PC_coef_filepath
      integer                      :: &
          nchanprof,             &
          rttov_direct_nthreads, &
          nchan_out,             &
          nchannels_rec,         &
          rttov_Nlocaltime
      real(wp)                     :: &
          CO2_mr,              &
          CH4_mr,              &
          CO_mr,               &
          N2O_mr,              &
          SO2_mr,              &
          ZenAng
      integer,allocatable          :: &
          iChannel(:),      &  ! Requested channel indices
          iChannel_out(:)      ! Passing out the channel indices (actual output channels)
      real(kind=wp),allocatable    :: &
          emisChannel(:),           &      ! RTTOV channel emissivity
          reflChannel(:),           &      ! RTTOV channel reflectivity
          rttov_localtime(:),       &      ! RTTOV localtime
          rttov_localtime_width(:) 
      logical, allocatable  :: &
          swath_mask(:)          
  end type rttov_cfg

  type rttov_output
      integer             :: &
          nchan_out
      integer,pointer     :: &
          channel_indices(:)
      real(wp),pointer    :: &
          bt_total(:,:),    &
          bt_clear(:,:),    &
          rad_total(:,:),   &
          rad_clear(:,:),   &
          rad_cloudy(:,:),  &
          refl_total(:,:),  &
          refl_clear(:,:),  &
          bt_total_pc(:,:), &
          rad_total_pc(:,:)
  end type rttov_output

  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! END MODULE
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE MOD_COSP_RTTOV_UTIL
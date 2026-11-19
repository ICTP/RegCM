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
! Apr 2015 - D. Swales - Modified for RTTOVv11.3
! Jun 2025 - J.K. Shaw - Added RTTOVv13.2 integration and swathing. rttov_cfg moved to cosp_rttov_util.F90
! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
MODULE MOD_COSP_RTTOV_INTERFACE
  USE COSP_KINDS,          ONLY: wp
  USE MOD_COSP_RTTOV,      ONLY: rttov_IN
  USE MOD_COSP_RTTOV_UTIL, ONLY: rttov_cfg,             rttov_output
  IMPLICIT NONE


CONTAINS

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE cosp_rttov_init
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_RTTOV_INIT(Lrttov,Nlevels,Ninstruments,instrument_namelists,       &
                             rttov_configs,unitn,debug)

      logical,intent(inout) :: &
          Lrttov
      integer,intent(in) :: &
          Nlevels,   &
          Ninstruments
      type(character(len=128)), dimension(Ninstruments)     :: & 
          instrument_namelists   ! Array of paths to RTTOV instrument namelists      
      type(rttov_cfg), dimension(:), intent(out), allocatable :: & ! intent(out)?
          rttov_configs
      integer,intent(in),Optional :: unitn ! Used for io limits
      logical,intent(in),Optional :: debug
          
      Lrttov = .false.
      allocate(rttov_configs(Ninstruments))
      
      print*,'Running COSP_RTTOV_INIT from STUB files.', &
        'To run RTTOV, compile COSP after setting environmental variable "RTTOV"'
       
  END SUBROUTINE COSP_RTTOV_INIT
  
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE DESTROY_RTTOV_CONFIG
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE DESTROY_RTTOV_CONFIG(rttovConfig)
  
      type(rttov_cfg),intent(inout) :: &
          rttovConfig
    
      if (allocated(rttovConfig % iChannel))              deallocate(rttovConfig % iChannel)
      if (allocated(rttovConfig % iChannel_out))          deallocate(rttovConfig % iChannel_out)
      if (allocated(rttovConfig % emisChannel))           deallocate(rttovConfig % emisChannel)
      if (allocated(rttovConfig % reflChannel))           deallocate(rttovConfig % reflChannel)
      if (allocated(rttovConfig % rttov_localtime))       deallocate(rttovConfig % rttov_localtime)
      if (allocated(rttovConfig % rttov_localtime_width)) deallocate(rttovConfig % rttov_localtime_width)
      if (allocated(rttovConfig % swath_mask))            deallocate(rttovConfig % swath_mask)

  END SUBROUTINE DESTROY_RTTOV_CONFIG  

  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! SUBROUTINE cosp_rttov_simulate
  !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  SUBROUTINE COSP_RTTOV_SIMULATE(rttovIN,rttovConfig,error,                        & ! Inputs
                                 bt_total,bt_clear,                                & ! Brightness Temp Outputs
                                 rad_total,rad_clear,rad_cloudy,                   & ! Radiance Outputs
                                 refl_total,refl_clear,                            & ! Reflectance Outputs
                                 debug)      

    type(rttov_in),intent(in) :: &
        rttovIN
    type(rttov_cfg),intent(inout) :: &
        rttovConfig
    character(len=128) :: &
        error     ! Error messages (only populated if error encountered)         
    real(wp),intent(inout),dimension(rttovIN%nPoints,rttovConfig%nchan_out),optional :: &
        bt_total,                          &        ! All-sky
        bt_clear,                          &        ! Clear-sky
        rad_total,                         &        ! All-sky
        rad_clear,                         &        ! Clear-sky
        rad_cloudy,                        &        ! Cloudy-sky
        refl_total,                        &        ! All-sky
        refl_clear                                  ! Clear-sky        
    logical,intent(in),optional :: &
        debug        

    print*,'Running COSP_RTTOV_SIMULATE from STUB files.', &
             'To run RTTOV, compile COSP after setting environmental variable "RTTOV"'

  END SUBROUTINE COSP_RTTOV_SIMULATE
  
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  ! END MODULE
  ! %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
END MODULE MOD_COSP_RTTOV_INTERFACE

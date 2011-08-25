!===============================================================================
! SVN $Id: dshr_const.F90 554 2006-03-28 21:11:44Z kauff $
! SVN $URL: https://svn-ccsm-models.cgd.ucar.edu/csm_share/trunk_tags/share3_070321/dshr/dshr_const.F90 $
!===============================================================================
!BOP ===========================================================================
!
! !MODULE: dshr_const -- defines/provides common constants for data models
!
! !DESCRIPTION:
!     dshr shared constansts
!
! !REVISION HISTORY:
!     2004-Dec-10 - J. Schramm - first version
!
! !INTERFACE: ------------------------------------------------------------------

module dshr_const
 
! !USES:

   use dshr_kind      ! defines F90 kinds
   use shr_const_mod  ! shared physical constants

   implicit none

! !PUBLIC TYPES:

  ! none

! !PUBLIC MEMBER FUNCTIONS:

  ! none

! !PUBLIC DATA MEMBERS:

   public

   !----------------------------------------------------------------------------
   ! physical constants
   !----------------------------------------------------------------------------
   real(R8),parameter :: dshr_const_pi      = SHR_CONST_PI           ! pi
   real(R8),parameter :: dshr_const_rearth  = SHR_CONST_REARTH       ! radius of earth ~ m
   real(R8),parameter :: dshr_const_rearth2 = SHR_CONST_REARTH*SHR_CONST_REARTH ! rad**2
   real(R8),parameter :: dshr_const_g       = SHR_CONST_G            ! gravity
   real(R8),parameter :: dshr_const_deg2rad = dshr_const_pi/180.0_R8 ! deg to rads
   real(R8),parameter :: dshr_const_rad2deg = 180.0_R8/dshr_const_pi ! rad to degs
   real(R8),parameter :: dshr_const_cDay    = SHR_CONST_CDAY         ! sec in siderial day ~ sec 

   real(R8),parameter :: dshr_const_pstd    = SHR_CONST_PSTD    ! standard pressure ~ Pa
   real(R8),parameter :: dshr_const_rdair   = SHR_CONST_RDAIR   ! Dry air gas constant ~ J/K/kg
   real(R8),parameter :: dshr_const_rhosw   = SHR_CONST_RHOSW   ! density of sea water ~ kg/m^3
   real(R8),parameter :: dshr_const_cpdair  = SHR_CONST_CPDAIR  ! spec heat of dry air
   real(R8),parameter :: dshr_const_cpwv    = SHR_CONST_CPWV    ! spec heat of h2o vapor
   real(R8),parameter :: dshr_const_cpsw    = SHR_CONST_CPSW    ! specific heat of sea h2o ~ J/kg/K
   real(R8),parameter :: dshr_const_cpvir   = dshr_const_cpwv/dshr_const_cpdair - 1.0_R8
   real(R8),parameter :: dshr_const_zvir    = SHR_CONST_ZVIR    ! rh2o/rair   - 1.0
   real(R8),parameter :: dshr_const_latvap  = SHR_CONST_LATVAP  ! latent heat of evap
   real(R8),parameter :: dshr_const_latice  = SHR_CONST_LATICE  ! latent heat of fusion
   real(R8),parameter :: dshr_const_stebol  = SHR_CONST_STEBOL  ! Stefan-Boltzmann
   real(R8),parameter :: dshr_const_karman  = SHR_CONST_KARMAN  ! Von Karman constant
   real(R8),parameter :: dshr_const_tkfrz   = SHR_CONST_TKFRZ   ! freezing T of fresh water ~ K 
   real(R8),parameter :: dshr_const_tkfrzsw = SHR_CONST_TKFRZSW ! freezing T of salt water  ~ K

   real(R8),parameter :: dshr_const_ocn_ref_sal = SHR_CONST_OCN_REF_SAL ! ocn ref salt
   real(R8),parameter :: dshr_const_ice_ref_sal = SHR_CONST_ICE_REF_SAL ! ice ref salt

   real(R8),parameter :: dshr_const_spval       = SHR_CONST_SPVAL       ! special value

   !----------------------------------------------------------------------------
   ! commonly used values
   !----------------------------------------------------------------------------
   real(R8),parameter :: dshr_const_c0 = 0.0_R8        ! zero
   real(R8),parameter :: dshr_const_c1 = 1.0_R8        ! one

!EOP

!===============================================================================

end module dshr_const


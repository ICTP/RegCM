#include <misc.h>
#include <preproc.h>

module clm_varsur

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_varsur
!
! !DESCRIPTION:
! Module containing 2-d surface boundary data information
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use mod_bats_param , only : slmo , xmopor , iexsol
!
! !PUBLIC TYPES:
  implicit none
  save
!
! land model grid - moved to domainMod
!
! surface boundary data, these are all "gdc" local
!
  integer , allocatable :: vegxy(:,:) ! vegetation type
  real(r8), allocatable,target :: wtxy(:,:)  ! subgrid weights

  real(r8) ,allocatable :: pctspec(:)         ! percent of spec lunits wrt gcell
  real(r8) ,allocatable :: pctspecB(:)        ! percent of spec lunits wrt gcell

! abt below
  integer,  pointer :: landmask(:,:)     ! land mask: 1 = land. 0 = ocean. 3 = land/ocean
  real(r8), allocatable :: landfrac(:,:)     ! fractional land
  real(r8), pointer :: ht_rcm(:,:)       ! elevation from regcm
  real(r8), pointer :: init_tgb(:,:)     ! ICBC temperature used for soil temp initialization
  real(r8), pointer :: init_snow(:,:)     ! Snow amount used for soil initialization
  logical  :: init_grid                      ! call clm initialization or not true=do init
  real(r8) :: numdays                        ! number of days per year (used in shr_orb_mod)
  integer,  allocatable :: clm_soitex(:,:)   ! Used only for the RegCM Dust Model
  real(r8), allocatable :: glonc(:)          ! center longitude (used only for output)
  real(r8), allocatable :: glatc(:)          ! center latitude  (used only for output)
  real(r8), pointer :: satbrt_clm(:,:)   ! Landuse from BATS read in from DOMAIN.INFO
                                             ! for transfer to CLM
  integer :: r2cimask                        ! landmask evaluation method
                                             ! 1 = Using DOMAIN.INFO landuse type
                                             ! 2 = Using landfraction from RCMnavyoro.nc
  integer :: r2cilawrence_albedo             ! Apply lawrence 2007 albedo modifications
  integer :: r2coutfrq                       ! output frequency for clm in hours
                                             ! variable retrieved from regcm.in (clmfrq)
!  real(r8), allocatable :: c2rcosz1d(:)      ! used for gathering and sending to regcm
  real(r8), allocatable :: c2r_allout(:)     ! used for gathering and sending to regcm
  integer , allocatable :: omap_i(:)         ! used for gathering and sending to regcm
  integer , allocatable :: omap_j(:)         ! used for gathering and sending to regcm
  real(r8), allocatable :: clm2bats_veg(:,:) ! Only used for RegCM Dust model
  real(r8), allocatable :: clm_fracveg(:,:)  ! fraction of grid covered by vegetation

  integer :: cgaschem                        ! Gas-phase chem flag( cgaschem=1 for CBMZ)
  integer :: caerosol                        ! Flag for sending aerosol vars to RegCM
                                             ! particularly clm snow cover and fracveg

! abt above
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! Revised by Ahmed Tawfik
! 2005-11-01 Moved grid to domainMod, T Craig
!
!EOP
!-----------------------------------------------------------------------

end module clm_varsur

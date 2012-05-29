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

  real(r8) ,allocatable :: pctspec(:)        ! percent of spec lunits wrt gcell

! abt below
  integer,  pointer :: landmask(:,:)     ! land mask: 1 = land. 0 = ocean. 3 = land/ocean
  real(r8), allocatable :: landfrac(:,:)     ! fractional land
  real(r8), allocatable :: ht_rcm(:,:)       ! elevation from regcm
  real(r8), allocatable :: init_tgb(:,:)     ! ICBC temperature used for soil temp initialization
  logical  :: init_grid                      ! call clm initialization or not true=do init
  real(r8) :: numdays                        ! number of days per year (used in shr_orb_mod)
  integer,  allocatable :: clm_soitex(:,:)   ! Used only for the RegCM Dust Model
  real(r8), allocatable :: glonc(:)          ! center longitude (used only for output)
  real(r8), allocatable :: glatc(:)          ! center latitude  (used only for output)
  real(r8), allocatable :: satbrt_clm(:,:)   ! Landuse from BATS read in from DOMAIN.INFO
                                             ! for transfer to CLM
  integer :: r2cimask                        ! landmask evaluation method
                                             ! 1 = Using DOMAIN.INFO landuse type
                                             ! 2 = Using landfraction from RCMnavyoro.nc
  integer :: r2coutfrq                       ! output frequency for clm in hours
                                             ! variable retrieved from regcm.in (clmfrq)
!  real(r8), allocatable :: c2rcosz1d(:)      ! used for gathering and sending to regcm
  real(r8), allocatable :: c2r_allout(:)     ! used for gathering and sending to regcm
  integer , allocatable :: omap_i(:)         ! used for gathering and sending to regcm
  integer , allocatable :: omap_j(:)         ! used for gathering and sending to regcm
  real(r8), allocatable :: clm2bats_veg(:,:) ! Only used for RegCM Dust model 
  real(r8), allocatable :: clm_fracveg(:,:)  ! fraction of grid covered by vegetation 


! soil moisture initialization
  real(r8) :: slmo(22)                       ! surface moisture availability in fraction of one
  integer  :: iexsol(22)                     ! soil texture type based on BATS landuse types     
  real(r8) :: xmopor(12)                     ! fraction of soil that is voids                           
  data xmopor/.33,.36,.39,.42,.45,.48,.51,.54,.57,.6,.63,.66/
  data iexsol/6,6,6,6,7,8,6,3,6,6,5,12,6,6,6,6,5,6,6,6,12,8/
  data slmo/0.65,0.45,0.6,0.6,0.65,0.65,0.55,0.1, &
  0.9,0.8,0.2,0.9,0.9,1.0,1.0,0.5,0.5,0.65,0.6,0.6,0.5,0.6/
!
! abt above
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! Revised by Ahmed Tawfik
! 2005-11-01 Moved grid to domainMod, T Craig
!
!EOP
!-----------------------------------------------------------------------

end module clm_varsur

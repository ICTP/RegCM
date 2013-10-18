
module mod_clm_varsur

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_varsur
!
! !DESCRIPTION:
! Module containing 2-d surface boundary data information
!
! !USES:
  use mod_realkinds
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
  real(rk8), allocatable,target :: wtxy(:,:)  ! subgrid weights

  real(rk8),allocatable :: pctspec(:)         ! percent of spec lunits wrt gcell

  real(rk8), allocatable,target :: topoxy(:,:)  ! subgrid glacier_mec sfc elevation
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
! 2005-11-01 Moved grid to domainMod, T Craig
!
!EOP
!-----------------------------------------------------------------------

end module mod_clm_varsur

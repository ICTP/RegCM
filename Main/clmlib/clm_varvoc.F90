#include <misc.h>
#include <preproc.h>

module clm_varvoc

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: clm_varvoc
!
! !DESCRIPTION:
! Module containing information used for VOC flux calculations using
! MEGANv2.03 (non-isoprene) and Guenther 2006 equations (isoprene)
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use clm_varpar  , only : nvoc
!
! !PUBLIC TYPES:
  implicit none
  save
!  ! biogenic emission factor maps
  real(r8), pointer :: ef_iso(:)             ! isoprene emission factor map
  real(r8), pointer :: ef_apin(:)            ! a-pinene emission factor map
  real(r8), pointer :: ef_bpin(:)            ! b-pinene emission factor map
  real(r8), pointer :: ef_mbo(:)             ! methylbutenol emission factor map
  real(r8), pointer :: ef_limo(:)            ! limonene emission factor map
  real(r8), pointer :: ef_sabi(:)            ! sabinene emission factor map
  real(r8), pointer :: ef_myrc(:)            ! myrcene emission factor map
  real(r8), pointer :: ef_acar(:)            ! a-3carene emission factor map
  real(r8), pointer :: ef_ocim(:)            ! ocimene emission factor map
  real(r8), pointer :: ef_omtp(:)            ! other monoterps emission factor map
  real(r8), pointer :: ef_farn(:)            ! farnicene emission factor map
  real(r8), pointer :: ef_bcar(:)            ! b-caryophyllene emission factor map
  real(r8), pointer :: ef_osqt(:)            ! other sesquiterps emission factor map
  real(r8), pointer :: ef_meoh(:)            ! methanol emission factor map
  real(r8), pointer :: ef_acto(:)            ! acetone emission factor map
  real(r8), pointer :: ef_meth(:)            ! methane emission factor map
  real(r8), pointer :: ef_no(:)              ! NO,N2O,NH3 emission factor map
  real(r8), pointer :: ef_acta(:)            ! acetaldehyde-ethanol emission factor map
  real(r8), pointer :: ef_form(:)            ! formic_acid/acetic_acid/formaldehyde emission factor map
  real(r8), pointer :: ef_co(:)              ! carbon monoxide emission factor map
!
!  ! time accumulation indices/variables for VOCemissionMOD (gamma calc)
  integer :: n24                             ! number of timesteps in 24 hours
  integer :: n240                            ! number of timesteps in 240 hours
  integer :: c24                             ! counter for timesteps in 24 hours
  integer :: c240                            ! counter for timesteps in 240 hours
  integer :: n_start                         ! stores the first time step.  Used in VOCmod
!
  integer :: r2cefmap(nvoc)                  ! 1=use emission factor map
                                             ! 2=use pft to deteremine emission factors

  real(r8), pointer :: gamma_t_out(:)        ! Gamma temperature for output only
  real(r8), pointer :: gamma_p_out(:)        ! Gamma PAR for output only
  real(r8), pointer :: gamma_sm_out(:)       ! Gamma soil moisture for output only
  real(r8), pointer :: Eopt_out(:)
  real(r8), pointer :: Topt_out(:)
  real(r8), pointer :: Cpsun_out(:)
  real(r8), pointer :: Cpsh_out(:)
  real(r8), pointer :: alphasu_out(:)
  real(r8), pointer :: alphash_out(:)
  real(r8), pointer :: t24_out(:)
  real(r8), pointer :: t240_out(:)
  real(r8), pointer :: p240sh_out(:)
  real(r8), pointer :: p240su_out(:)
  real(r8), pointer :: p24sh_out(:)
  real(r8), pointer :: p24su_out(:)
  real(r8), pointer :: gamma_age_out(:)


! Species dependent variables used in gamma calculation
  real(r8) :: LDF(nvoc)                      ! surface moisture availability in fraction of one
  real(r8) :: beta_V(nvoc)                     ! soil texture type based on BATS landuse types
  real(r8) :: Anew(nvoc)                     ! relative emission activity for new foilage
  real(r8) :: Aold(nvoc)                     ! relative emission activity for old foilage
  real(r8) :: Agro(nvoc)                     ! relative emission activity for growing foilage
  real(r8) :: Amat(nvoc)                     ! relative emission activity for mature foilage
!  data LDF/0.9999_r8,0.05_r8,0.1_r8,0.5_r8,0.5_r8,0.9999_r8,0.1_r8,0.1_r8,0.8_r8,0.05_r8,0.1_r8, &
!           0.5_r8,0.5_r8,0.5_r8,0.75_r8,0.25_r8,0.75_r8,0.0_r8,0.5_r8,0.5_r8/
!  data beta_V/0.09_r8,0.1_r8,0.1_r8,0.1_r8,0.09_r8,0.09_r8,0.1_r8,0.1_r8,0.1_r8,0.1_r8,0.1_r8, &
!            0.17_r8,0.17_r8,0.17_r8,0.08_r8,0.11_r8,0.05_r8,0.11_r8,0.13_r8,0.09_r8/
!  data Anew/0.05_r8,2._r8,2._r8,2._r8,1._r8,0.05_r8,2._r8,2._r8,2._r8,2._r8,2._r8,0.4_r8,0.4_r8, &
!            0.4_r8,3._r8,0.05_r8,0.05_r8,0.05_r8,0.05_r8,0.05_r8/
!  data Aold/1._r8,1._r8,1._r8,1._r8,1._r8,1._r8,1._r8,1._r8,1._r8,1._r8,1._r8,1._r8,1._r8,1._r8, &
!            1._r8,1._r8,1._r8,1._r8,1._r8,1._r8/
!  data Agro/0.6_r8,1.8_r8,1.8_r8,1.8_r8,1._r8,0.6_r8,1.8_r8,1.8_r8,1.8_r8,1.8_r8,1.8_r8,0.6_r8,0.6_r8, &
!            0.6_r8,2.6_r8,0.6_r8,0.6_r8,0.6_r8,0.6_r8,0.6_r8/
!  data Amat/1.125_r8,0.95_r8,0.95_r8,0.95_r8,1._r8,1.125_r8,0.95_r8,0.95_r8,0.95_r8,0.95_r8,0.95_r8, &
!            1.075_r8,1.075_r8,1.075_r8,0.85_r8,1.125_r8,1.125_r8,1.125_r8,1.125_r8,1.125_r8/
! For above data /isoprene,myrcene,sabinene,limonene,CO,MBO,apin,bpin,ocimene,a-3carene,othermonoterps,
!                 farnicene,b-caryophyllene,other-sesquiterps,methanol,acetone,methane,N2O-NO-NH3,
!                 acetaldehyde-ethanol,formic_acid-formaldehyde-acetic_acid /
! !REVISION HISTORY:
! Created by Ahmed Tawfik
!
!EOP
!-----------------------------------------------------------------------

end module clm_varvoc
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

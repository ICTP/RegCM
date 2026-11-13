!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#ifdef OASIS

module mod_oasis_interface

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_stdio
  use mod_message
  use mod_service
  use mod_mppparam
  use mod_runparams, only : isocean, dtsec, alarm_out_sav, rcmtimer
  use mod_bats_common, only : rdnnsg
  use mod_atm_interface, only : atms, sfs, flwd, flw, fsw, sinc, mddom
  use mod_lm_interface, only : lms
  use mod_date, only : lfdomonth, lmidnight
#ifdef ECLM
  use mod_memutil
  use mod_runparams, only : iqv, ich4, iisop, ichem, dtsrf
  use mod_atm_interface, only : mdsub, solvs, solvsd, solvl, solvld
#endif

  use mod_oasis
  use mod_oasis_params
  use mod_oasis_signature
  use mod_oasis_generic

  implicit none

  private

  public :: comp_name, comp_id  ! -> mod_oasis_params
  public :: oasis_lag           ! -> mod_oasis_params

#ifdef ECLM
  ! work space on the land subgrid (nnsg,jci1:jci2,ici1:ici2)
  real(rkx), dimension(:,:,:), contiguous, pointer :: tmp_lnd => null( )
  ! index of isoprene in the RCM_FLXVOC bundle (0 if absent)
  integer(ik4) :: kisop = 0
#endif

  !--------------------------------------------------------------------
  ! for grid and partition activation
  ! (dot/cross, external/internal)
  logical :: l_make_grdde, l_make_grddi
  logical :: l_make_grdce, l_make_grdci
  type(infogrd), target, allocatable :: grdde, grddi, grdce, grdci

#ifdef ECLM
  ! cross, internal, land grid of the eCLM coupling: the ocean is masked
  logical :: l_make_grdlnd
  type(infogrd), target, allocatable :: grdlnd
#endif

  !--------------------------------------------------------------------
  ! below the namelist parameters (namelist oasisparam)
  public :: l_write_grids, write_restart_option ! -> mod_oasis_params
  public :: oasis_sync_lag                      ! -> mod_oasis_params
#ifdef ECLM
  public :: oasis_nlevsoi, oasis_nddvel         ! -> mod_oasis_params
  public :: oasis_megan_cmp, oasis_max_megan    ! -> mod_oasis_params
#endif
  logical, public :: l_cpl_im_sst,  &
                     l_cpl_im_wz0,  & ! not tested
                     l_cpl_im_wust, & ! not tested
                     l_cpl_ex_u10m, & ! not tested
                     l_cpl_ex_v10m, & ! not tested
                     l_cpl_ex_wspd, & ! not tested
                     l_cpl_ex_wdir, & ! not tested
                     l_cpl_ex_t2m,  & ! not tested
                     l_cpl_ex_q2m,  & ! not tested
                     l_cpl_ex_slp,  &
                     l_cpl_ex_taux, &
                     l_cpl_ex_tauy, &
                     l_cpl_ex_z0,   & ! not tested
                     l_cpl_ex_ustr, & ! not tested
                     l_cpl_ex_evap, &
                     l_cpl_ex_prec, &
                     l_cpl_ex_nuwa, &
                     l_cpl_ex_ulhf, &
                     l_cpl_ex_ushf, &
                     l_cpl_ex_uwlw, &
                     l_cpl_ex_dwlw, &
                     l_cpl_ex_nulw, &
                     l_cpl_ex_uwsw, &
                     l_cpl_ex_dwsw, &
                     l_cpl_ex_ndsw, &
                     l_cpl_ex_rhoa !,
!                    l_cpl_im_sit,  & ! not coded
!                    l_cpl_ex_t10m, & ! not coded
!                    l_cpl_ex_q10m, & ! not coded
!                    l_cpl_ex_pr      ! not coded
#ifdef ECLM
  logical, public ::    l_cpl_ex_tatm,   l_cpl_ex_uatm,    l_cpl_ex_vatm,    &
      l_cpl_ex_qatm,    l_cpl_ex_hgt,    l_cpl_ex_patm,    l_cpl_ex_dwrlwf,  &
      l_cpl_ex_swdir,   l_cpl_ex_swdif,  l_cpl_im_tgbb,    l_cpl_im_t2m,     &
      l_cpl_im_q2m,     l_cpl_im_u10m,   l_cpl_im_sent,    l_cpl_im_evpr,    &
      l_cpl_im_ram1,    l_cpl_im_rah1,   l_cpl_im_tauy,    l_cpl_im_taux,    &
      l_cpl_im_sncv,    l_cpl_im_wt,     l_cpl_im_zo,      l_cpl_im_tlef,    &
      l_cpl_ex_rainf,   l_cpl_ex_snow,   l_cpl_im_albei,   l_cpl_im_albed,   &
      l_cpl_im_flxvoc,  l_cpl_im_flxdst, l_cpl_im_fluxch4, l_cpl_im_ddvel,   &
      l_cpl_im_tlai,    l_cpl_im_roff,   l_cpl_im_srfroff, l_cpl_im_snwmelt, &
      l_cpl_im_tsoi,    l_cpl_im_h2ovol, l_cpl_im_h2oliq,  l_cpl_im_h2oice,  &
      l_cpl_im_h2o10cm, l_cpl_im_tgrd,   l_cpl_im_albei_2, l_cpl_im_albed_2
#endif
  ! OASIS field +++

  !--------------------------------------------------------------------
  ! below the variables for the fields' information
  type(infofld), allocatable :: im_sst,  &
                                im_wz0,  &
                                im_wust, &
                                ex_u10m, &
                                ex_v10m, &
                                ex_wspd, &
                                ex_wdir, &
                                ex_t2m,  &
                                ex_q2m,  &
                                ex_slp,  &
                                ex_taux, &
                                ex_tauy, &
                                ex_z0,   &
                                ex_ustr, &
                                ex_evap, &
                                ex_prec, &
                                ex_nuwa, &
                                ex_ulhf, &
                                ex_ushf, &
                                ex_uwlw, &
                                ex_dwlw, &
                                ex_nulw, &
                                ex_uwsw, &
                                ex_dwsw, &
                                ex_ndsw, &
                                ex_rhoa  !,
!                               im_sit,  & ! not coded
!                               ex_t10m, & ! not coded
!                               ex_q10m    ! not coded
#ifdef ECLM
  type(infofld), allocatable :: ex_tatm, ex_uatm, ex_vatm, ex_qatm,  &
    ex_hgt, ex_patm, ex_dwrlwf, ex_swdir, ex_swdif, im_tgbb, im_t2m, &
    im_q2m, im_u10m, im_sent, im_evpr, im_ram1, im_rah1, im_taux,    &
    im_tauy, im_wt, im_sncv, im_zo, im_tlef, ex_rainf, ex_snow,      &
    im_fluxch4, im_tlai, im_tgrd, im_roff, im_srfroff, im_snwmelt,   &
    im_albei, im_albed, im_albei_2, im_albed_2, im_h2o10cm
  type(info3dfld), allocatable :: im_flxdst, im_flxvoc, im_ddvel,    &
    im_h2oliq, im_h2oice, im_h2ovol, im_tsoi
#endif
  ! OASIS field +++

  ! OASIS needs double precision arrays. These variables are optional.
  ! Required for the fields to import; recommended for the fields to
  ! export which need a bit of work from the model variables.
  real(rkx), dimension(:,:), allocatable :: cpl_sst,  &
                                            cpl_wz0,  &
                                            cpl_wust, &
                                            cpl_wdir
!                                           cpl_sit
#ifdef ECLM
  real(rkx), dimension(:,:), allocatable :: cpl_t2m, cpl_q2m, cpl_u10m,   &
    cpl_sent, cpl_evpr, cpl_ram1, cpl_rah1, cpl_taux, cpl_tauy, cpl_sncv, &
    cpl_wt, cpl_zo, cpl_albei, cpl_albed, cpl_albei_2, cpl_albed_2,       &
    cpl_tlef, cpl_tgbb, cpl_fluxch4, cpl_tlai, cpl_tgrd, cpl_roff,        &
    cpl_srfroff, cpl_snwmelt, cpl_h2o10cm
  real(rkx), dimension(:,:,:), allocatable :: cpl_flxdst, cpl_flxvoc, &
    cpl_h2oliq, cpl_ddvel, cpl_h2oice, cpl_h2ovol, cpl_tsoi
#endif

  ! OASIS field +++

  ! OASIS needs double precision arrays. These variables are optional.
  ! They are used when calling the subroutine oasisxregcm_setup_field_array(),
  ! when an array is to be initialized. If not given, the array will be
  ! initialized with zeros.
  real(rkx), parameter :: init_sst = tzero + 25.0
  ! OASIS field +++

  !--------------------------------------------------------------------
  ! accessibility from oasislib:
  public :: oasisxregcm_init, oasisxregcm_finalize ! -> mod_oasis_generic
  public :: oasisxregcm_header                      ! -> mod_oasis_signature

  !--------------------------------------------------------------------
  ! from this module:
  public :: oasisxregcm_params
  public :: oasisxregcm_release
  public :: oasisxregcm_def
  public :: oasisxregcm_rcv_all
  public :: oasisxregcm_snd_all
  public :: oasisxregcm_sync_wait

  contains

  subroutine oasisxregcm_params
    implicit none
#ifdef ECLM
    integer(ik4) :: nmegan
#endif
    !--------------------------------------------------------------------------
    oasis_lag = 0
    ! state which grids have to be defined,
    ! depending on what field is activated.
    l_make_grdde   = .false.
    ! OASIS field +++
    l_make_grddi   = .false.
    ! OASIS field +++
    l_make_grdce   = l_cpl_ex_slp
    ! OASIS field +++
    l_make_grdci   = l_cpl_im_sst  .or. &
                     l_cpl_im_wz0  .or. &
                     l_cpl_im_wust .or. &
                     l_cpl_ex_u10m .or. &
                     l_cpl_ex_v10m .or. &
                     l_cpl_ex_wspd .or. &
                     l_cpl_ex_wdir .or. &
                     l_cpl_ex_t2m  .or. &
                     l_cpl_ex_q2m  .or. &
                     l_cpl_ex_taux .or. &
                     l_cpl_ex_tauy .or. &
                     l_cpl_ex_z0   .or. &
                     l_cpl_ex_ustr .or. &
                     l_cpl_ex_evap .or. &
                     l_cpl_ex_prec .or. &
                     l_cpl_ex_nuwa .or. &
                     l_cpl_ex_ulhf .or. &
                     l_cpl_ex_ushf .or. &
                     l_cpl_ex_uwlw .or. &
                     l_cpl_ex_dwlw .or. &
                     l_cpl_ex_nulw .or. &
                     l_cpl_ex_uwsw .or. &
                     l_cpl_ex_dwsw .or. &
                     l_cpl_ex_ndsw .or. &
                     l_cpl_ex_rhoa ! .or.
!                    l_cpl_im_sit  .or. &
!                    l_cpl_ex_t10m .or. &
!                    l_cpl_ex_q10m
#ifdef ECLM
    l_make_grdlnd = l_cpl_ex_tatm    .or. l_cpl_ex_uatm    .or. &
                    l_cpl_ex_vatm    .or. l_cpl_ex_qatm    .or. &
                    l_cpl_ex_hgt     .or. l_cpl_ex_patm    .or. &
                    l_cpl_ex_dwrlwf  .or. l_cpl_ex_swdir   .or. &
                    l_cpl_ex_swdif   .or. l_cpl_ex_rainf   .or. &
                    l_cpl_ex_snow    .or. l_cpl_im_tgbb    .or. &
                    l_cpl_im_t2m     .or. l_cpl_im_q2m     .or. &
                    l_cpl_im_u10m    .or. l_cpl_im_sent    .or. &
                    l_cpl_im_evpr    .or. l_cpl_im_ram1    .or. &
                    l_cpl_im_rah1    .or. l_cpl_im_taux    .or. &
                    l_cpl_im_tauy    .or. l_cpl_im_sncv    .or. &
                    l_cpl_im_wt      .or. l_cpl_im_zo      .or. &
                    l_cpl_im_tgrd    .or. l_cpl_im_tlef    .or. &
                    l_cpl_im_tlai    .or. l_cpl_im_roff    .or. &
                    l_cpl_im_srfroff .or. l_cpl_im_snwmelt .or. &
                    l_cpl_im_h2o10cm .or. l_cpl_im_albei   .or. &
                    l_cpl_im_albed   .or. l_cpl_im_albei_2 .or. &
                    l_cpl_im_albed_2 .or. l_cpl_im_fluxch4 .or. &
                    l_cpl_im_flxdst  .or. l_cpl_im_flxvoc  .or. &
                    l_cpl_im_ddvel   .or. l_cpl_im_h2oliq  .or. &
                    l_cpl_im_h2oice  .or. l_cpl_im_tsoi    .or. &
                    l_cpl_im_h2ovol
#endif
    ! OASIS field +++
    !
    ! initialize grids: grid variable, name with no mask, name with mask,
    !                   local j1, j2, i1, i2, global j1, j2, i1, i2, number of
    !                   corners per cell
    if ( l_make_grdde ) then
      call oasisxregcm_setup_grid(grdde, 'rden', 'rdem', &
            jde1, jde2, ide1, ide2, 1, jx, 1, iy, 4)
    end if
    if ( l_make_grddi ) then
      call oasisxregcm_setup_grid(grddi, 'rdin', 'rdim', &
            jdi1, jdi2, idi1, idi2, 2, jx-1, 2, iy-1, 4)
    end if
    if ( l_make_grdce ) then
      call oasisxregcm_setup_grid(grdce, 'rcen', 'rcem', &
            jce1, jce2, ice1, ice2, 1, jx-1, 1, iy-1, 4)
    end if
    if ( l_make_grdci ) then
      call oasisxregcm_setup_grid(grdci, 'rcin', 'rcim', &
            jci1, jci2, ici1, ici2, 2, jx-2, 2, iy-2, 4)
    end if
#ifdef ECLM
    ! Same points as grdci, but the mask selects the land: its own names,
    ! so that ocean and land couplings can coexist.
    if ( l_make_grdlnd ) then
      call oasisxregcm_setup_grid(grdlnd, 'rlin', 'rlim', &
            jci1, jci2, ici1, ici2, 2, jx-2, 2, iy-2, 4)
      call getmem(tmp_lnd,1,nnsg,jci1,jci2,ici1,ici2,'oasis:tmp_lnd')
    end if
#endif
    !
    ! initialize fields: field variable, name, grid, field array (optional),
    ! initialization value (optional; 0 otherwise)
    if ( l_cpl_im_sst ) then
      call oasisxregcm_setup_field(im_sst, 'RCM_SST', grdci, cpl_sst, init_sst)
    end if
    if ( l_cpl_im_wz0 ) then
      call oasisxregcm_setup_field(im_wz0, 'RCM_WZ0', grdci, cpl_wz0)
    end if
    if ( l_cpl_im_wust ) then
      call oasisxregcm_setup_field(im_wust, 'RCM_WUST', grdci, cpl_wust)
    end if
!   if ( l_cpl_im_sit ) then
!     call oasisxregcm_setup_field(im_sit, 'RCM_SIT', grdci, cpl_sit)
!   end if
    if ( l_cpl_ex_u10m ) then
      call oasisxregcm_setup_field(ex_u10m, 'RCM_U10M', grdci)
    end if
    if ( l_cpl_ex_v10m ) then
      call oasisxregcm_setup_field(ex_v10m, 'RCM_V10M', grdci)
    end if
    if ( l_cpl_ex_wspd ) then
      call oasisxregcm_setup_field(ex_wspd, 'RCM_WSPD', grdci)
    end if
    if ( l_cpl_ex_wdir ) then
      call oasisxregcm_setup_field(ex_wdir, 'RCM_WDIR', grdci, cpl_wdir)
    end if
    if ( l_cpl_ex_t2m ) then
      call oasisxregcm_setup_field(ex_t2m, 'RCM_T2M', grdci)
    end if
    if ( l_cpl_ex_q2m ) then
      call oasisxregcm_setup_field(ex_q2m, 'RCM_Q2M', grdci)
    end if
!   if ( l_cpl_ex_t10m ) then
!     call oasisxregcm_setup_field(ex_t10m, 'RCM_T10M', grdci)
!   end if
!   if ( l_cpl_ex_q10m ) then
!     call oasisxregcm_setup_field(ex_q10m, 'RCM_Q10M', grdci)
!   end if
    if ( l_cpl_ex_slp ) then
      call oasisxregcm_setup_field(ex_slp, 'RCM_SLP', grdce)
    end if
    if ( l_cpl_ex_taux ) then
      call oasisxregcm_setup_field(ex_taux, 'RCM_TAUX', grdci)
    end if
    if ( l_cpl_ex_tauy ) then
      call oasisxregcm_setup_field(ex_tauy, 'RCM_TAUY', grdci)
    end if
    if ( l_cpl_ex_z0 ) then
      call oasisxregcm_setup_field(ex_z0,   'RCM_Z0',   grdci)
    end if
    if ( l_cpl_ex_ustr ) then
      call oasisxregcm_setup_field(ex_ustr, 'RCM_USTR', grdci)
    end if
    if ( l_cpl_ex_evap ) then
      call oasisxregcm_setup_field(ex_evap, 'RCM_EVAP', grdci)
    end if
    if ( l_cpl_ex_prec ) then
      call oasisxregcm_setup_field(ex_prec, 'RCM_PREC', grdci)
    end if
    if ( l_cpl_ex_nuwa ) then
      call oasisxregcm_setup_field(ex_nuwa, 'RCM_NUWA', grdci)
    end if
    if ( l_cpl_ex_ulhf ) then
      call oasisxregcm_setup_field(ex_ulhf, 'RCM_ULHF', grdci)
    end if
    if ( l_cpl_ex_ushf ) then
      call oasisxregcm_setup_field(ex_ushf, 'RCM_USHF', grdci)
    end if
    if ( l_cpl_ex_uwlw ) then
      call oasisxregcm_setup_field(ex_uwlw, 'RCM_UWLW', grdci)
    end if
    if ( l_cpl_ex_dwlw ) then
      call oasisxregcm_setup_field(ex_dwlw, 'RCM_DWLW', grdci)
    end if
    if ( l_cpl_ex_nulw ) then
      call oasisxregcm_setup_field(ex_nulw, 'RCM_NULW', grdci)
    end if
    if ( l_cpl_ex_uwsw ) then
      call oasisxregcm_setup_field(ex_uwsw, 'RCM_UWSW', grdci)
    end if
    if ( l_cpl_ex_dwsw ) then
      call oasisxregcm_setup_field(ex_dwsw, 'RCM_DWSW', grdci)
    end if
    if ( l_cpl_ex_ndsw ) then
      call oasisxregcm_setup_field(ex_ndsw, 'RCM_NDSW', grdci)
    end if
    if ( l_cpl_ex_rhoa ) then
      call oasisxregcm_setup_field(ex_rhoa, 'RCM_RHOA', grdci)
    end if
#ifdef ECLM
    if ( l_cpl_ex_tatm ) then
      call oasisxregcm_setup_field(ex_tatm, 'RCM_TATM', grdlnd)
    end if
    if ( l_cpl_ex_uatm ) then
      call oasisxregcm_setup_field(ex_uatm, 'RCM_UATM', grdlnd)
    end if
    if ( l_cpl_ex_vatm ) then
      call oasisxregcm_setup_field(ex_vatm, 'RCM_VATM', grdlnd)
    end if
    if ( l_cpl_ex_qatm ) then
      call oasisxregcm_setup_field(ex_qatm, 'RCM_QATM', grdlnd)
    end if
    if ( l_cpl_ex_hgt ) then
      call oasisxregcm_setup_field(ex_hgt, 'RCM_HGT', grdlnd)
    end if
    if ( l_cpl_ex_patm ) then
      call oasisxregcm_setup_field(ex_patm, 'RCM_PATM', grdlnd)
    end if
    if ( l_cpl_ex_dwrlwf ) then
      call oasisxregcm_setup_field(ex_dwrlwf, 'RCM_DWRLWF', grdlnd)
    end if
    if ( l_cpl_ex_swdir ) then
      call oasisxregcm_setup_field(ex_swdir, 'RCM_SWDIR', grdlnd)
    end if
    if ( l_cpl_ex_swdif ) then
      call oasisxregcm_setup_field(ex_swdif, 'RCM_SWDIF', grdlnd)
    end if
    if ( l_cpl_ex_rainf ) then
      call oasisxregcm_setup_field(ex_rainf, 'RCM_RAINF', grdlnd)
    end if
    if ( l_cpl_ex_snow ) then
      call oasisxregcm_setup_field(ex_snow, 'RCM_SNOW', grdlnd)
    end if
    if ( l_cpl_im_t2m ) then
      call oasisxregcm_setup_field(im_t2m, 'RCM_T2M_In', grdlnd, cpl_t2m)
    end if
    if ( l_cpl_im_q2m ) then
      call oasisxregcm_setup_field(im_q2m, 'RCM_Q2M_In', grdlnd, cpl_q2m)
    end if
    if ( l_cpl_im_u10m ) then
      call oasisxregcm_setup_field(im_u10m, 'RCM_U10M_In', grdlnd, cpl_u10m)
    end if
    if ( l_cpl_im_sent ) then
      call oasisxregcm_setup_field(im_sent, 'RCM_SENT', grdlnd, cpl_sent)
    end if
    if ( l_cpl_im_evpr ) then
      call oasisxregcm_setup_field(im_evpr, 'RCM_EVPR', grdlnd, cpl_evpr)
    end if
    if ( l_cpl_im_ram1 ) then
      call oasisxregcm_setup_field(im_ram1, 'RCM_RAM1', grdlnd, cpl_ram1)
    end if
    if ( l_cpl_im_rah1 ) then
      call oasisxregcm_setup_field(im_rah1, 'RCM_RAH1', grdlnd, cpl_rah1)
    end if
    if ( l_cpl_im_taux ) then
      call oasisxregcm_setup_field(im_taux, 'RCM_TAUX_In', grdlnd, cpl_taux)
    end if
    if ( l_cpl_im_tauy ) then
      call oasisxregcm_setup_field(im_tauy, 'RCM_TAUY_In', grdlnd, cpl_tauy)
    end if
    if ( l_cpl_im_sncv ) then
      call oasisxregcm_setup_field(im_sncv, 'RCM_SNCV', grdlnd, cpl_sncv)
    end if
    if ( l_cpl_im_wt ) then
      call oasisxregcm_setup_field(im_wt, 'RCM_WT', grdlnd, cpl_wt)
    end if
    if ( l_cpl_im_zo ) then
      call oasisxregcm_setup_field(im_zo, 'RCM_ZO', grdlnd, cpl_zo)
    end if
    if ( l_cpl_im_tlef ) then
      call oasisxregcm_setup_field(im_tlef, 'RCM_TLEF', grdlnd, cpl_tlef)
    end if
    if ( l_cpl_im_albei ) then
      call oasisxregcm_setup_field(im_albei, 'RCM_ALBEI', grdlnd, cpl_albei)
    end if
    if ( l_cpl_im_albed ) then
      call oasisxregcm_setup_field(im_albed, 'RCM_ALBED', grdlnd, cpl_albed)
    end if
    if ( l_cpl_im_albei_2 ) then
      call oasisxregcm_setup_field(im_albei_2, 'RCM_ALBEI_2', &
              grdlnd, cpl_albei_2)
    end if
    if ( l_cpl_im_albed_2 ) then
      call oasisxregcm_setup_field(im_albed_2, 'RCM_ALBED_2', &
              grdlnd, cpl_albed_2)
    end if
    if ( l_cpl_im_fluxch4 ) then
      call oasisxregcm_setup_field(im_fluxch4, 'RCM_FLUXCH4', &
              grdlnd, cpl_fluxch4)
    end if
    if ( l_cpl_im_tlai ) then
      call oasisxregcm_setup_field(im_tlai, 'RCM_TLAI', grdlnd, cpl_tlai)
    end if
    if ( l_cpl_im_roff ) then
      call oasisxregcm_setup_field(im_roff, 'RCM_ROFF', grdlnd, cpl_roff)
    end if
    if ( l_cpl_im_srfroff ) then
      call oasisxregcm_setup_field(im_srfroff, 'RCM_SRFROFF', &
              grdlnd, cpl_srfroff)
    end if
    if ( l_cpl_im_snwmelt ) then
      call oasisxregcm_setup_field(im_snwmelt, 'RCM_SNWMELT', &
              grdlnd, cpl_snwmelt)
    end if
    if ( l_cpl_im_tgrd ) then
      call oasisxregcm_setup_field(im_tgrd, 'RCM_TGRD', grdlnd, cpl_tgrd)
    end if
    if ( l_cpl_im_tgbb ) then
      call oasisxregcm_setup_field(im_tgbb, 'RCM_TGBB', grdlnd, cpl_tgbb)
    end if
    if ( l_cpl_im_h2o10cm ) then
      call oasisxregcm_setup_field(im_h2o10cm, 'RCM_H2O10CM', &
              grdlnd, cpl_h2o10cm)
    end if
    ! Bundles: their sizes must match the eCLM configuration
    nmegan = count(oasis_megan_cmp /= ' ')
    kisop = findloc(oasis_megan_cmp,'ISOP',1)
    if ( l_cpl_im_flxvoc ) then
      call oasisxregcm_setup_3dfield(im_flxvoc, 'RCM_FLXVOC', grdlnd, &
              nmegan, cpl_flxvoc)
    end if
    if ( l_cpl_im_flxdst ) then
      call oasisxregcm_setup_3dfield(im_flxdst, 'RCM_FLXDST', grdlnd, &
              oasis_ndust, cpl_flxdst)
    end if
    if ( l_cpl_im_ddvel ) then
      call oasisxregcm_setup_3dfield(im_ddvel, 'RCM_DDVEL', grdlnd, &
              oasis_nddvel, cpl_ddvel)
    end if
    if ( l_cpl_im_tsoi ) then
      call oasisxregcm_setup_3dfield(im_tsoi, 'RCM_TSOI', grdlnd, &
              oasis_nlevsoi, cpl_tsoi)
    end if
    if ( l_cpl_im_h2ovol ) then
      call oasisxregcm_setup_3dfield(im_h2ovol, 'RCM_H2OVOL', grdlnd, &
              oasis_nlevsoi, cpl_h2ovol)
    end if
    if ( l_cpl_im_h2oliq ) then
      call oasisxregcm_setup_3dfield(im_h2oliq, 'RCM_H2OLIQ', grdlnd, &
              oasis_nlevsoi, cpl_h2oliq)
    end if
    if ( l_cpl_im_h2oice ) then
      call oasisxregcm_setup_3dfield(im_h2oice, 'RCM_H2OICE', grdlnd, &
              oasis_nlevsoi, cpl_h2oice)
    end if
#endif
    ! OASIS field +++
  end subroutine oasisxregcm_params

  ! call all definition subroutines for setting up OASIS
  subroutine oasisxregcm_def(localComm)
    implicit none
    integer(kind=ip_intwp_p), intent(in) :: localComm  !< MPI communicator
    !--------------------------------------------------------------------------
    ! partition definition
#ifdef DEBUG
    write(ndebug,*) oasis_prefix, 'definition phase: partitions'
#endif
    if ( l_make_grdde ) call oasisxregcm_def_partition(grdde)
    if ( l_make_grddi ) call oasisxregcm_def_partition(grddi)
    if ( l_make_grdce ) call oasisxregcm_def_partition(grdce)
    if ( l_make_grdci ) call oasisxregcm_def_partition(grdci)
#ifdef ECLM
    if ( l_make_grdlnd ) call oasisxregcm_def_partition(grdlnd)
#endif
    ! grid definition
#ifdef DEBUG
    write(ndebug,*) oasis_prefix, 'definition phase: grids'
#endif
    if ( l_write_grids ) then
      call oasisxregcm_def_grid
    else
#ifdef DEBUG
      write(ndebug,*) oasis_prefix, '               >> skipped!'
#endif
      if ( myid == italk ) write(stdout,*) 'Note: OASIS grid files', &
                                   ' have not been written online.'
    end if
    ! variable definitions
#ifdef DEBUG
    write(ndebug,*) oasis_prefix, 'definition phase: fields'
#endif
    ! field variable, OASIS integer for the direction:
    !                 importing(OASIS_In)/exporting(OASIS_Out)
    if ( l_cpl_im_sst )  call oasisxregcm_def_field(im_sst,  OASIS_In)
    if ( l_cpl_im_wz0 )  call oasisxregcm_def_field(im_wz0,  OASIS_In)
    if ( l_cpl_im_wust ) call oasisxregcm_def_field(im_wust, OASIS_In)
    if ( l_cpl_ex_u10m ) call oasisxregcm_def_field(ex_u10m, OASIS_Out)
    if ( l_cpl_ex_v10m ) call oasisxregcm_def_field(ex_v10m, OASIS_Out)
    if ( l_cpl_ex_wspd ) call oasisxregcm_def_field(ex_wspd, OASIS_Out)
    if ( l_cpl_ex_wdir ) call oasisxregcm_def_field(ex_wdir, OASIS_Out)
    if ( l_cpl_ex_t2m )  call oasisxregcm_def_field(ex_t2m,  OASIS_Out)
    if ( l_cpl_ex_q2m )  call oasisxregcm_def_field(ex_q2m,  OASIS_Out)
    if ( l_cpl_ex_slp )  call oasisxregcm_def_field(ex_slp,  OASIS_Out)
    if ( l_cpl_ex_taux ) call oasisxregcm_def_field(ex_taux, OASIS_Out)
    if ( l_cpl_ex_tauy ) call oasisxregcm_def_field(ex_tauy, OASIS_Out)
    if ( l_cpl_ex_z0 )   call oasisxregcm_def_field(ex_z0,   OASIS_Out)
    if ( l_cpl_ex_ustr ) call oasisxregcm_def_field(ex_ustr, OASIS_Out)
    if ( l_cpl_ex_evap ) call oasisxregcm_def_field(ex_evap, OASIS_Out)
    if ( l_cpl_ex_prec ) call oasisxregcm_def_field(ex_prec, OASIS_Out)
    if ( l_cpl_ex_nuwa ) call oasisxregcm_def_field(ex_nuwa, OASIS_Out)
    if ( l_cpl_ex_ulhf ) call oasisxregcm_def_field(ex_ulhf, OASIS_Out)
    if ( l_cpl_ex_ushf ) call oasisxregcm_def_field(ex_ushf, OASIS_Out)
    if ( l_cpl_ex_uwlw ) call oasisxregcm_def_field(ex_uwlw, OASIS_Out)
    if ( l_cpl_ex_dwlw ) call oasisxregcm_def_field(ex_dwlw, OASIS_Out)
    if ( l_cpl_ex_nulw ) call oasisxregcm_def_field(ex_nulw, OASIS_Out)
    if ( l_cpl_ex_uwsw ) call oasisxregcm_def_field(ex_uwsw, OASIS_Out)
    if ( l_cpl_ex_dwsw ) call oasisxregcm_def_field(ex_dwsw, OASIS_Out)
    if ( l_cpl_ex_ndsw ) call oasisxregcm_def_field(ex_ndsw, OASIS_Out)
    if ( l_cpl_ex_rhoa ) call oasisxregcm_def_field(ex_rhoa, OASIS_Out)
!   if ( l_cpl_im_sit )  call oasisxregcm_def_field(im_sit,  OASIS_In)
!   if ( l_cpl_ex_t10m ) call oasisxregcm_def_field(ex_t10m, OASIS_Out)
!   if ( l_cpl_ex_q10m ) call oasisxregcm_def_field(ex_q10m, OASIS_Out)
#ifdef ECLM
    if ( l_cpl_ex_tatm ) call oasisxregcm_def_field(ex_tatm, OASIS_Out)
    if ( l_cpl_ex_uatm ) call oasisxregcm_def_field(ex_uatm, OASIS_Out)
    if ( l_cpl_ex_vatm ) call oasisxregcm_def_field(ex_vatm, OASIS_Out)
    if ( l_cpl_ex_qatm ) call oasisxregcm_def_field(ex_qatm, OASIS_Out)
    if ( l_cpl_ex_hgt ) call oasisxregcm_def_field(ex_hgt, OASIS_Out)
    if ( l_cpl_ex_patm ) call oasisxregcm_def_field(ex_patm, OASIS_Out)
    if ( l_cpl_ex_dwrlwf ) call oasisxregcm_def_field(ex_dwrlwf, OASIS_Out)
    if ( l_cpl_ex_swdir ) call oasisxregcm_def_field(ex_swdir, OASIS_Out)
    if ( l_cpl_ex_swdif ) call oasisxregcm_def_field(ex_swdif, OASIS_Out)
    if ( l_cpl_ex_rainf ) call oasisxregcm_def_field(ex_rainf, OASIS_Out)
    if ( l_cpl_ex_snow ) call oasisxregcm_def_field(ex_snow, OASIS_Out)
    if ( l_cpl_im_t2m ) call oasisxregcm_def_field(im_t2m, OASIS_In)
    if ( l_cpl_im_q2m ) call oasisxregcm_def_field(im_q2m, OASIS_In)
    if ( l_cpl_im_u10m ) call oasisxregcm_def_field(im_u10m, OASIS_In)
    if ( l_cpl_im_sent ) call oasisxregcm_def_field(im_sent, OASIS_In)
    if ( l_cpl_im_evpr ) call oasisxregcm_def_field(im_evpr, OASIS_In)
    if ( l_cpl_im_ram1 ) call oasisxregcm_def_field(im_ram1, OASIS_In)
    if ( l_cpl_im_rah1 ) call oasisxregcm_def_field(im_rah1, OASIS_In)
    if ( l_cpl_im_taux ) call oasisxregcm_def_field(im_taux, OASIS_In)
    if ( l_cpl_im_tauy ) call oasisxregcm_def_field(im_tauy, OASIS_In)
    if ( l_cpl_im_sncv ) call oasisxregcm_def_field(im_sncv, OASIS_In)
    if ( l_cpl_im_wt ) call oasisxregcm_def_field(im_wt, OASIS_In)
    if ( l_cpl_im_zo ) call oasisxregcm_def_field(im_zo, OASIS_In)
    if ( l_cpl_im_tlef ) call oasisxregcm_def_field(im_tlef, OASIS_In)
    if ( l_cpl_im_albei ) call oasisxregcm_def_field(im_albei, OASIS_In)
    if ( l_cpl_im_albed ) call oasisxregcm_def_field(im_albed, OASIS_In)
    if ( l_cpl_im_albei_2 ) call oasisxregcm_def_field(im_albei_2, OASIS_In)
    if ( l_cpl_im_albed_2 ) call oasisxregcm_def_field(im_albed_2, OASIS_In)
    if ( l_cpl_im_fluxch4 ) call oasisxregcm_def_field(im_fluxch4, OASIS_In)
    if ( l_cpl_im_tlai ) call oasisxregcm_def_field(im_tlai, OASIS_In)
    if ( l_cpl_im_roff ) call oasisxregcm_def_field(im_roff, OASIS_In)
    if ( l_cpl_im_srfroff ) call oasisxregcm_def_field(im_srfroff, OASIS_In)
    if ( l_cpl_im_snwmelt ) call oasisxregcm_def_field(im_snwmelt, OASIS_In)
    if ( l_cpl_im_tgrd ) call oasisxregcm_def_field(im_tgrd, OASIS_In)
    if ( l_cpl_im_tgbb ) call oasisxregcm_def_field(im_tgbb, OASIS_In)
    if ( l_cpl_im_h2o10cm ) call oasisxregcm_def_field(im_h2o10cm, OASIS_In)
    if ( l_cpl_im_flxvoc ) call oasisxregcm_def_3dfield(im_flxvoc, OASIS_In)
    if ( l_cpl_im_flxdst ) call oasisxregcm_def_3dfield(im_flxdst, OASIS_In)
    if ( l_cpl_im_ddvel ) call oasisxregcm_def_3dfield(im_ddvel, OASIS_In)
    if ( l_cpl_im_tsoi ) call oasisxregcm_def_3dfield(im_tsoi, OASIS_In)
    if ( l_cpl_im_h2ovol ) call oasisxregcm_def_3dfield(im_h2ovol, OASIS_In)
    if ( l_cpl_im_h2oliq ) call oasisxregcm_def_3dfield(im_h2oliq, OASIS_In)
    if ( l_cpl_im_h2oice ) call oasisxregcm_def_3dfield(im_h2oice, OASIS_In)
#endif
    ! OASIS field +++
    ! termination of definition phase
#ifdef DEBUG
    write(ndebug,*) oasis_prefix, 'definition phase: end'
#endif
    call oasisxregcm_end_def(localComm)

  contains

  ! define OASIS grids
  subroutine oasisxregcm_def_grid
    implicit none
    ! dot / cross degree coordinates
    real(rkx), pointer, contiguous, dimension(:,:) :: dlon, dlat
    real(rkx), pointer, contiguous, dimension(:,:) :: xlon, xlat
    ! land category (15 for ocean)
    real(rkx), pointer, contiguous, dimension(:,:) :: lndcat
    real(rkx), pointer, contiguous, dimension(:,:) :: oasisgrid_lon
    real(rkx), pointer, contiguous, dimension(:,:) :: oasisgrid_lat
    ! note: 4 corners (always the case?) in the third dimension,
    !       counterclockwisely.
    real(rkx), pointer, contiguous, dimension(:,:,:) :: oasisgrid_clon
    real(rkx), pointer, contiguous, dimension(:,:,:) :: oasisgrid_clat
    ! surface of the grid meshes m2
    real(rkx), pointer, contiguous, dimension(:,:) :: oasisgrid_srf
    ! masks, 0 = valid, 1 = mask (OASIS convention)
    integer(ik4), pointer, contiguous, dimension(:,:) :: oasisgrid_mask
    integer(ik4) :: il_flag ! flag for grid writing by proc 0
    integer(ik4) :: ierror
#ifdef DEBUG
    !--------------------------------------------------------------------------
    write(ndebug,*) oasis_prefix, 'start collecting grid data'
    write(ndebug,*) oasis_prefix, 'grids writing is done by the in/out cpu only'
#endif
    ! collect the domain definition information in RegCM
    allocate(dlon(jx,iy),stat=ierror)
    if ( ierror /= 0 ) then
      write(stderr,*) 'error allocating dlon'
      call fatal(__FILE__,__LINE__,'OASIS GRID')
    end if
    allocate(dlat(jx,iy),stat=ierror)
    if ( ierror /= 0 ) then
      write(stderr,*) 'error allocating dlat'
      call fatal(__FILE__,__LINE__,'OASIS GRID')
    end if
    allocate(xlon(jx,iy),stat=ierror)
    if ( ierror /= 0 ) then
      write(stderr,*) 'error allocating xlon'
      call fatal(__FILE__,__LINE__,'OASIS GRID')
    end if
    allocate(xlat(jx,iy),stat=ierror)
    if ( ierror /= 0 ) then
      write(stderr,*) 'error allocating xlat'
      call fatal(__FILE__,__LINE__,'OASIS GRID')
    end if
    allocate(lndcat(jx,iy),stat=ierror)
    if ( ierror /= 0 ) then
      write(stderr,*) 'error allocating lndcat'
      call fatal(__FILE__,__LINE__,'OASIS GRID')
    end if
    call grid_collect(mddom%dlon,dlon,jde1,jde2,ide1,ide2)
    call grid_collect(mddom%dlat,dlat,jde1,jde2,ide1,ide2)
    call grid_collect(mddom%xlon,xlon,jde1,jde2,ide1,ide2)
    call grid_collect(mddom%xlat,xlat,jde1,jde2,ide1,ide2)
    call grid_collect(mddom%lndcat,lndcat,jde1,jde2,ide1,ide2)
    if ( myid == iocpu ) then
      ! start the grid writing process
#ifdef DEBUG
      write(ndebug,*) oasis_prefix, 'start grids writing'
#endif
      call oasis_start_grids_writing(il_flag)
      !
      if ( l_make_grdde .or. l_make_grddi) then ! dots
        !
        call oasisxregcm_make_oasisgrids_d( &
             oasisgrid_lon,oasisgrid_lat,oasisgrid_clon,oasisgrid_clat, &
             oasisgrid_srf,oasisgrid_mask, &
             dlon,dlat,xlon,xlat,lndcat) ! subroutine writen below
        !
        if ( l_make_grdde ) call oasisxregcm_write_oasisgrids(grdde, &
             oasisgrid_lon,oasisgrid_lat,oasisgrid_clon,oasisgrid_clat, &
             oasisgrid_srf,oasisgrid_mask)
        !
        if ( l_make_grddi ) call oasisxregcm_write_oasisgrids(grddi, &
             oasisgrid_lon,oasisgrid_lat,oasisgrid_clon,oasisgrid_clat, &
             oasisgrid_srf,oasisgrid_mask)
        !
        call oasisxregcm_deallocate_oasisgrids( &
             oasisgrid_lon,oasisgrid_lat,oasisgrid_clon,oasisgrid_clat, &
             oasisgrid_srf,oasisgrid_mask)
        !
      end if
      !
#ifdef ECLM
      if ( l_make_grdce .or. l_make_grdci .or. l_make_grdlnd ) then ! crosses
#else
      if ( l_make_grdce .or. l_make_grdci ) then ! crosses
#endif
        !
        call oasisxregcm_make_oasisgrids_c( &
             oasisgrid_lon,oasisgrid_lat,oasisgrid_clon, &
             oasisgrid_clat,oasisgrid_srf,oasisgrid_mask, &
             dlon,dlat,xlon,xlat,lndcat) ! subroutine writen below
        !
        if ( l_make_grdce ) call oasisxregcm_write_oasisgrids(grdce, &
             oasisgrid_lon,oasisgrid_lat,oasisgrid_clon,oasisgrid_clat, &
             oasisgrid_srf,oasisgrid_mask)
        !
        if ( l_make_grdci ) call oasisxregcm_write_oasisgrids(grdci, &
             oasisgrid_lon,oasisgrid_lat,oasisgrid_clon,oasisgrid_clat, &
             oasisgrid_srf,oasisgrid_mask)
        !
#ifdef ECLM
        ! the land grid masks the ocean: invert the ocean grid mask
        if ( l_make_grdlnd ) then
          oasisgrid_mask = 1 - oasisgrid_mask
          call oasisxregcm_write_oasisgrids(grdlnd, &
               oasisgrid_lon,oasisgrid_lat,oasisgrid_clon,oasisgrid_clat, &
               oasisgrid_srf,oasisgrid_mask)
        end if
#endif
        call oasisxregcm_deallocate_oasisgrids( &
             oasisgrid_lon,oasisgrid_lat,oasisgrid_clon,oasisgrid_clat, &
             oasisgrid_srf,oasisgrid_mask)
        !
      end if
      ! terminate the grid writing process
#ifdef DEBUG
      write(ndebug,*) oasis_prefix, 'terminate grids writing'
#endif
      call oasis_terminate_grids_writing()
    end if
    if ( associated(dlon) ) nullify(dlon)
    if ( associated(dlat) ) nullify(dlat)
    if ( associated(xlon) ) nullify(xlon)
    if ( associated(xlat) ) nullify(xlat)
    if ( associated(lndcat) ) nullify(lndcat)
#ifdef DEBUG
    write(ndebug,"(' ',A,A,I3,A)") oasis_prefix, &
         ': please refer to cpu ', iocpu, &
         ' for grid writing debug statements'
#endif
  end subroutine oasisxregcm_def_grid

  subroutine oasisxregcm_make_oasisgrids_d(lon,lat,clon,clat,srf,mask, &
                                           dlon,dlat,xlon,xlat,lndcat)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:), intent(in) :: dlon, dlat
    real(rkx), pointer, contiguous, dimension(:,:), intent(in) :: xlon, xlat
    real(rkx), pointer, contiguous, dimension(:,:), intent(in) :: lndcat
    real(rkx), pointer, contiguous, dimension(:,:), intent(out) :: lon, lat
    real(rkx), pointer, contiguous, dimension(:,:), intent(out) :: srf
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(out) :: clon, clat
    integer(ik4), pointer, contiguous, dimension(:,:), intent(out) :: mask
    integer(ik4) :: i, j
    !--------------------------------------------------------------------------
    !
    !    ^    corner 2        corner 1
    !    |    x(j-1,i)         x(j,i)
    !  i |
    !    |             center
    !  a |             d(j,i)
    !  x |
    !  i |    corner 3        corner 4
    !  s |   x(j-1,i-1)       x(j,i-1)
    !    |
    !    +-------------------------------->
    !                 j axis
    !
    ! note: concretely xlat and xlon are defined on (1:jx,1:iy) as for
    !       dlat and dlon though the cross grid is read on (1:jx-1,1:iy-1)
    ! only consequence: xlon and xlat exist at (jx,jy) -> no exception needed
    ! however, exceptions are needed for extrem west and south values
    !                              0         0
    call oasisxregcm_allocate_oasisgrids(lon,lat,clon,clat,srf,mask, &
                                         jx,iy,4)
    do i = 1, iy
      do j = 1, jx
        ! center
        lon(j,i) = dlon(j,i)
        lat(j,i) = dlat(j,i)
        ! upper right corner: no problem
        clon(j,i,1) = xlon(j,i)
        clat(j,i,1) = xlat(j,i)
        ! left corners:
        !   upper left corner: use of j-1 -> exception for j=1
        !                      lon -> use of the last lon step to the West
        !                      lat -> use of the same value as for j
        !                             (expected to be the same)
        !   lower left corner: same + use of i-1 -> exception for i=1
        !                      lon -> use of the same value as for i
        !                             (expected to be the same)
        !                      lat -> use of the last lat step to the South
        if ( j == 1 ) then
          clon(j,i,2) = 2 * xlon(j,i) - xlon(j+1,i)
          clat(j,i,2) = xlat(j,i)
          if ( i == 1 ) then
            clon(j,i,3) = 2 * xlon(j,i) - xlon(j+1,i)
            clat(j,i,3) = 2 * xlat(j,i) - xlat(j,i+1)
          else
            clon(j,i,3) = 2 * xlon(j,i-1) - xlon(j+1,i-1)
            clat(j,i,3) = xlat(j,i-1)
          end if
        else
          clon(j,i,2) = xlon(j-1,i)
          clat(j,i,2) = xlat(j-1,i)
          if ( i == 1 ) then
            clon(j,i,3) = xlon(j-1,i)
            clat(j,i,3) = 2 * xlat(j-1,i) - xlat(j-1,i+1)
          else
            clon(j,i,3) = xlon(j-1,i-1)
            clat(j,i,3) = xlat(j-1,i-1)
          end if
        end if
        ! lower right corner: use of i-1 -> exception for i=1
        if ( i == 1 ) then
          clon(j,i,4) = xlon(j,i)
          clat(j,i,4) = 2 * xlat(j,i) - xlat(j,i+1)
        else
          clon(j,i,4) = xlon(j,i-1)
          clat(j,i,4) = xlat(j,i-1)
        end if
        ! surface
        srf(j,i) = srf_sqm(clon(j,i,:),clat(j,i,:))
        ! mask
        if ( isocean(lndcat(j,i)) ) then
          mask(j,i) = 0
        else
          mask(j,i) = 1
        end if
      end do
    end do
  end subroutine oasisxregcm_make_oasisgrids_d

  subroutine oasisxregcm_make_oasisgrids_c(lon,lat,clon,clat,srf,mask, &
                                           dlon,dlat,xlon,xlat,lndcat)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:), intent(in) :: dlon, dlat
    real(rkx), pointer, contiguous, dimension(:,:), intent(in) :: xlon, xlat
    real(rkx), pointer, contiguous, dimension(:,:), intent(in) :: lndcat
    real(rkx), pointer, contiguous, dimension(:,:), intent(out) :: lon, lat
    real(rkx), pointer, contiguous, dimension(:,:), intent(out) :: srf
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(out) :: clon, clat
    integer(ik4), pointer, contiguous, dimension(:,:), intent(out) :: mask
    integer(ik4) :: i, j
    !--------------------------------------------------------------------------
    !
    !    ^    corner 2        corner 1
    !    |    d(j,i+1)       d(j+1,i+1)
    !  i |
    !    |             center
    !  a |             x(j,i)
    !  x |
    !  i |    corner 3        corner 4
    !  s |     d(j,i)         d(j+1,i)
    !    |
    !    +-------------------------------->
    !                 j axis
    !
    ! note: the cross case is easier as the cross grid is contained
    !       in the dot one thus all corners are already defined in dlon
    !
    call oasisxregcm_allocate_oasisgrids(lon,lat,clon,clat,srf,mask, &
                                         jx-1,iy-1,4)
    do i = 1, iy-1
      do j = 1, jx-1
        ! center
        lon(j,i) = xlon(j,i)
        lat(j,i) = xlat(j,i)
        ! upper right corner
        clon(j,i,1) = dlon(j+1,i+1)
        clat(j,i,1) = dlat(j+1,i+1)
        ! upper left corner
        clon(j,i,2) = dlon(j,i+1)
        clat(j,i,2) = dlat(j,i+1)
        ! lower left corner
        clon(j,i,3) = dlon(j,i)
        clat(j,i,3) = dlat(j,i)
        ! lower right corner
        clon(j,i,4) = dlon(j+1,i)
        clat(j,i,4) = dlat(j+1,i)
        ! surface
        srf(j,i) = srf_sqm(clon(j,i,:),clat(j,i,:))
        ! mask
        if ( isocean(lndcat(j,i)) ) then
          mask(j,i) = 0
        else
          mask(j,i) = 1
        end if
      end do
    end do
  end subroutine oasisxregcm_make_oasisgrids_c

  end subroutine oasisxregcm_def

  ! call all subroutines consisting of receiving OASIS fields
  ! and optionally reworking them
  subroutine oasisxregcm_rcv_all(time)
    implicit none
    integer(ik4), intent(in) :: time ! execution time
    logical :: l_act
#ifdef ECLM
    integer(ik4) :: i, j, n
    ! molar mass ratio CH4/C: eCLM gives the methane flux as carbon
    real(rkx), parameter :: ch4_per_c = 16.04_rkx/12.011_rkx
#endif
    !--------------------------------------------------------------------------
#ifdef DEBUG
    if ( time == 0 ) then
      write(ndebug,*) oasis_prefix, 'Initialization:'
    else if ( debug_level > 0 ) then
      write(ndebug,*) oasis_prefix, 'OASIS time: ', time
    end if
#endif
    !
    if ( l_cpl_im_sst ) then ! sea surface temperature [K]
      call oasisxregcm_rcv(cpl_sst,im_sst,time,l_act)
      ! fill the ocean parts: out array (2d or 3d with nsg dimension)
      !                       in array, mask, field grid
      call fill_ocean(sfs%tg,cpl_sst,mddom%lndcat,im_sst%grd)
      call fill_ocean(sfs%tgbb,cpl_sst,mddom%lndcat,im_sst%grd)
    end if
    ! NOT IMPLEMENTED YET, BELOW IS A COPY OF THE PROCESS IN
    ! MOD_LM_INTERFACE.F90; IMPORT_DATA_INTO_SURFACE.
!    if ( l_cpl_im_sit ) then ! sea ice thickness [m]
!      call oasisxregcm_rcv(cpl_sit,im_sit,time,l_act)
!!      if ( l_act ) then
!        grd => im_sit%grd
!        ! dans le mask
!        ! tol = 0.5 E 20
!          if ( impfie%sit(j,i) < tol .and. lm%ldmsk(j,i) /= 1 ) then
!            if ( impfie%sit(j,i) > iceminh ) then
!              flag = .false.
!              if ( lm%ldmsk(j,i) == 0 ) flag = .true.
!              ! set land-sea mask
!              lm%ldmsk(j,i) = 2
!              do n = 1, nnsg
!                lm%ldmsk1(n,j,i) = 2
!                ! set sea ice thikness (in meter)
!                lms%sfice(n,j,i) = impfie%sit(j,i)
!              end do
!              ! write debug info
!              if ( flag ) then
!                write(*,f99002) j, i, 'water', 'ice  ', &
!                   lm%ldmsk(j,i), lms%sfice(1,j,i)
!              end if
!            else
!              if ( ldmskb(j,i) == 0 .and. lm%ldmsk(j,i) == 2 ) then
!                do n = 1, nnsg
!                  ! reduce to one tenth surface ice: it should melt away
!                  lms%sfice(n,j,i) = lms%sfice(n,j,i)*d_r10
!                  ! check that sea ice is melted or not
!                  if ( lms%sfice(n,j,i) <= iceminh ) then
!                    if ( ldmskb(j,i) /= lm%ldmsk(j,i) ) flag = .true.
!                    ! set land-sea mask to its original value
!                    lm%ldmsk(j,i) = ldmskb(j,i)
!                    lm%ldmsk1(n,j,i) = ldmskb(j,i)
!                    ! set land-use type to its original value
!                    ! set sea ice thikness (in meter)
!                    lms%sfice(n,j,i) = d_zero
!                  else
!                    flag = .false.
!                  end if
!                end do
!                ! write debug info
!                if ( flag ) then
!                  write(*,f99003) j, i, 'ice  ', 'water',  &
!                    lm%ldmsk(j,i), lms%sfice(1,j,i)
!                end if
!              end if
!            end if
!          end if
!        ???
!        nullify(grd)
!!      end if
!    end if
    !
    if ( l_cpl_im_wz0 ) then ! surface roughness length [m]
      call oasisxregcm_rcv(cpl_wz0,im_wz0,time,l_act)
      call fill_ocean(sfs%zo,cpl_wz0,mddom%lndcat,im_wz0%grd)
    end if
    !
    if ( l_cpl_im_wust ) then ! surface friction velocity [s-1]
      call oasisxregcm_rcv(cpl_wust,im_wust,time,l_act)
      call fill_ocean(sfs%ustar,cpl_wust,mddom%lndcat,im_wust%grd)
    end if
#ifdef ECLM
    if ( l_cpl_im_t2m ) then ! 2 m air temperature [K]
      call oasisxregcm_rcv(cpl_t2m,im_t2m,time,l_act)
      call fill_land(lms%t2m,cpl_t2m,mddom%lndcat,im_t2m%grd)
    end if
    if ( l_cpl_im_q2m ) then ! 2 m specific humidity [kg/kg]
      call oasisxregcm_rcv(cpl_q2m,im_q2m,time,l_act)
      call fill_land(lms%q2m,cpl_q2m,mddom%lndcat,im_q2m%grd)
    end if
    if ( l_cpl_im_u10m ) then ! 10 m wind speed [m/s]
      ! eCLM gives the speed only: keep the direction of the wind at the
      ! lowest model level.
      call oasisxregcm_rcv(cpl_u10m,im_u10m,time,l_act)
      call fill_land(tmp_lnd,cpl_u10m,mddom%lndcat,im_u10m%grd)
      do concurrent ( n = 1:nnsg, j = jci1:jci2, i = ici1:ici2 )
        if ( mdsub%ldmsk(n,j,i) == 1 ) then
          lms%u10m(n,j,i) = tmp_lnd(n,j,i) * atms%ubx3d(j,i,kz) / &
            max(sqrt(atms%ubx3d(j,i,kz)**2+atms%vbx3d(j,i,kz)**2),0.01_rkx)
          lms%v10m(n,j,i) = tmp_lnd(n,j,i) * atms%vbx3d(j,i,kz) / &
            max(sqrt(atms%ubx3d(j,i,kz)**2+atms%vbx3d(j,i,kz)**2),0.01_rkx)
        end if
      end do
    end if
    if ( l_cpl_im_sent ) then ! sensible heat flux [W/m2]
      call oasisxregcm_rcv(cpl_sent,im_sent,time,l_act)
      call fill_land(lms%sent,cpl_sent,mddom%lndcat,im_sent%grd)
    end if
    if ( l_cpl_im_evpr ) then ! evapotranspiration [kg/m2/s]
      call oasisxregcm_rcv(cpl_evpr,im_evpr,time,l_act)
      call fill_land(lms%evpr,cpl_evpr,mddom%lndcat,im_evpr%grd)
    end if
    if ( l_cpl_im_ram1 ) then
      call oasisxregcm_rcv(cpl_ram1,im_ram1,time,l_act)
      call fill_land(lms%ram1,cpl_ram1,mddom%lndcat,im_ram1%grd)
    end if
    if ( l_cpl_im_rah1 ) then
      call oasisxregcm_rcv(cpl_rah1,im_rah1,time,l_act)
      call fill_land(lms%rah1,cpl_rah1,mddom%lndcat,im_rah1%grd)
    end if
    if ( l_cpl_im_taux .and. l_cpl_im_tauy ) then ! wind stress [N/m2]
      call oasisxregcm_rcv(cpl_taux,im_taux,time,l_act)
      call fill_land(lms%taux,cpl_taux,mddom%lndcat,im_taux%grd)
      call oasisxregcm_rcv(cpl_tauy,im_tauy,time,l_act)
      call fill_land(lms%tauy,cpl_tauy,mddom%lndcat,im_tauy%grd)
      ! drag as used by the atmosphere: tau = drag * V
      do concurrent ( n = 1:nnsg, j = jci1:jci2, i = ici1:ici2 )
        if ( mdsub%ldmsk(n,j,i) == 1 ) then
          lms%drag(n,j,i) = &
            sqrt(lms%taux(n,j,i)**2 + lms%tauy(n,j,i)**2) / &
            max(sqrt(atms%ubx3d(j,i,kz)**2+atms%vbx3d(j,i,kz)**2),0.01_rkx)
        end if
      end do
    end if
    if ( l_cpl_im_sncv ) then ! snow water equivalent [kg/m2]
      call oasisxregcm_rcv(cpl_sncv,im_sncv,time,l_act)
      call fill_land(lms%sncv,cpl_sncv,mddom%lndcat,im_sncv%grd)
    end if
    if ( l_cpl_im_wt ) then
      call oasisxregcm_rcv(cpl_wt,im_wt,time,l_act)
      call fill_land(lms%wt,cpl_wt,mddom%lndcat,im_wt%grd)
    end if
    if ( l_cpl_im_zo ) then ! roughness length [m]
      call oasisxregcm_rcv(cpl_zo,im_zo,time,l_act)
      call fill_land(tmp_lnd,cpl_zo,mddom%lndcat,im_zo%grd)
      do concurrent ( n = 1:nnsg, j = jci1:jci2, i = ici1:ici2 )
        if ( mdsub%ldmsk(n,j,i) == 1 ) lms%zo(n,j,i) = tmp_lnd(n,j,i)
      end do
    end if
    if ( l_cpl_im_tlef ) then ! leaf temperature [K]
      call oasisxregcm_rcv(cpl_tlef,im_tlef,time,l_act)
      call fill_land(lms%tlef,cpl_tlef,mddom%lndcat,im_tlef%grd)
    end if
    if ( l_cpl_im_tlai ) then ! leaf area index
      call oasisxregcm_rcv(cpl_tlai,im_tlai,time,l_act)
      call fill_land(lms%xlai,cpl_tlai,mddom%lndcat,im_tlai%grd)
    end if
    if ( l_cpl_im_roff ) then ! total runoff
      call oasisxregcm_rcv(cpl_roff,im_roff,time,l_act)
      call fill_land(lms%trnof,cpl_roff,mddom%lndcat,im_roff%grd)
    end if
    if ( l_cpl_im_srfroff ) then ! surface runoff
      call oasisxregcm_rcv(cpl_srfroff,im_srfroff,time,l_act)
      call fill_land(lms%srnof,cpl_srfroff,mddom%lndcat,im_srfroff%grd)
    end if
    if ( l_cpl_im_snwmelt ) then ! snow melt rate, stored as amount
      call oasisxregcm_rcv(cpl_snwmelt,im_snwmelt,time,l_act)
      call fill_land(tmp_lnd,cpl_snwmelt,mddom%lndcat,im_snwmelt%grd)
      do concurrent ( n = 1:nnsg, j = jci1:jci2, i = ici1:ici2 )
        if ( mdsub%ldmsk(n,j,i) == 1 ) then
          lms%snwm(n,j,i) = tmp_lnd(n,j,i) * dtsrf
        end if
      end do
    end if
    if ( l_cpl_im_h2o10cm ) then ! soil water in the top 10 cm [kg/m2]
      call oasisxregcm_rcv(cpl_h2o10cm,im_h2o10cm,time,l_act)
      call fill_land(lms%tsw,cpl_h2o10cm,mddom%lndcat,im_h2o10cm%grd)
    end if
    ! Albedos. eCLM gives 1 where the sun is down: use a land value there.
    if ( l_cpl_im_albed ) then
      call oasisxregcm_rcv(cpl_albed,im_albed,time,l_act)
      call fill_albedo(lms%swdiralb,cpl_albed,im_albed%grd)
    end if
    if ( l_cpl_im_albed_2 ) then
      call oasisxregcm_rcv(cpl_albed_2,im_albed_2,time,l_act)
      call fill_albedo(lms%lwdiralb,cpl_albed_2,im_albed_2%grd)
    end if
    if ( l_cpl_im_albei ) then
      call oasisxregcm_rcv(cpl_albei,im_albei,time,l_act)
      call fill_albedo(lms%swdifalb,cpl_albei,im_albei%grd)
    end if
    if ( l_cpl_im_albei_2 ) then
      call oasisxregcm_rcv(cpl_albei_2,im_albei_2,time,l_act)
      call fill_albedo(lms%lwdifalb,cpl_albei_2,im_albei_2%grd)
    end if
    if ( l_cpl_im_albed .or. l_cpl_im_albei .or. &
         l_cpl_im_albed_2 .or. l_cpl_im_albei_2 ) then
      do concurrent ( n = 1:nnsg, j = jci1:jci2, i = ici1:ici2 )
        if ( mdsub%ldmsk(n,j,i) == 1 ) then
          lms%swalb(n,j,i) = max(lms%swdiralb(n,j,i),lms%swdifalb(n,j,i))
          lms%lwalb(n,j,i) = max(lms%lwdiralb(n,j,i),lms%lwdifalb(n,j,i))
        end if
      end do
    end if
    if ( l_cpl_im_tgrd ) then ! ground temperature [K]
      call oasisxregcm_rcv(cpl_tgrd,im_tgrd,time,l_act)
      call fill_land(lms%tgrd,cpl_tgrd,mddom%lndcat,im_tgrd%grd)
    end if
    if ( l_cpl_im_tgbb ) then ! radiative temperature [K]
      call oasisxregcm_rcv(cpl_tgbb,im_tgbb,time,l_act)
      call fill_land(lms%tgbb,cpl_tgbb,mddom%lndcat,im_tgbb%grd)
    end if
    ! Soil profiles, oasis_nlevsoi layers
    if ( l_cpl_im_tsoi ) then ! soil temperature [K]
      call oasisxregcm_rcv_3d(cpl_tsoi,im_tsoi,time,l_act)
      call fill_land(lms%tsoi,cpl_tsoi,mddom%lndcat,im_tsoi%grd)
    end if
    if ( l_cpl_im_h2ovol ) then ! volumetric soil water [m3/m3]
      call oasisxregcm_rcv_3d(cpl_h2ovol,im_h2ovol,time,l_act)
      call fill_land(lms%sw_vol,cpl_h2ovol,mddom%lndcat,im_h2ovol%grd)
    end if
    if ( l_cpl_im_h2oliq ) then ! liquid soil water [kg/m2]
      call oasisxregcm_rcv_3d(cpl_h2oliq,im_h2oliq,time,l_act)
    end if
    if ( l_cpl_im_h2oice ) then ! frozen soil water [kg/m2]
      call oasisxregcm_rcv_3d(cpl_h2oice,im_h2oice,time,l_act)
    end if
    if ( l_cpl_im_h2oliq .and. l_cpl_im_h2oice ) then
      call fill_land(lms%h2osoi,cpl_h2oliq+cpl_h2oice, &
                     mddom%lndcat,im_h2oliq%grd)
    else if ( l_cpl_im_h2oliq ) then
      call fill_land(lms%h2osoi,cpl_h2oliq,mddom%lndcat,im_h2oliq%grd)
    else if ( l_cpl_im_h2oice ) then
      call fill_land(lms%h2osoi,cpl_h2oice,mddom%lndcat,im_h2oice%grd)
    end if
    ! Emissions and deposition velocities for the chemistry
    if ( l_cpl_im_flxdst ) then ! dust emission per size bin [kg/m2/s]
      call oasisxregcm_rcv_3d(cpl_flxdst,im_flxdst,time,l_act)
      call fill_land(lms%dustemiss,cpl_flxdst,mddom%lndcat,im_flxdst%grd)
    end if
    if ( l_cpl_im_flxvoc ) then ! MEGAN emissions [kg/m2/s]
      call oasisxregcm_rcv_3d(cpl_flxvoc,im_flxvoc,time,l_act)
      ! Only isoprene has a tracer in the chemistry mechanisms
      if ( ichem == 1 .and. iisop > 0 .and. kisop > 0 ) then
        call fill_land(tmp_lnd,cpl_flxvoc(:,:,kisop), &
                       mddom%lndcat,im_flxvoc%grd)
        do concurrent ( n = 1:nnsg, j = jci1:jci2, i = ici1:ici2 )
          if ( mdsub%ldmsk(n,j,i) == 1 ) then
            lms%vocemiss(n,j,i,iisop) = tmp_lnd(n,j,i)
          end if
        end do
      end if
    end if
    if ( l_cpl_im_fluxch4 ) then ! methane emission [kgC/m2/s]
      call oasisxregcm_rcv(cpl_fluxch4,im_fluxch4,time,l_act)
      if ( ichem == 1 .and. ich4 > 0 ) then
        call fill_land(tmp_lnd,cpl_fluxch4,mddom%lndcat,im_fluxch4%grd)
        do concurrent ( n = 1:nnsg, j = jci1:jci2, i = ici1:ici2 )
          if ( mdsub%ldmsk(n,j,i) == 1 ) then
            lms%vocemiss(n,j,i,ich4) = tmp_lnd(n,j,i) * ch4_per_c
          end if
        end do
      end if
    end if
    ! The order of the species in RCM_DDVEL is not known here: the field
    ! is received but not passed to the chemistry.
    if ( l_cpl_im_ddvel ) then
      call oasisxregcm_rcv_3d(cpl_ddvel,im_ddvel,time,l_act)
    end if
#endif
    ! OASIS field +++
    !
#ifdef ECLM
  contains

    subroutine fill_albedo(alb,cpl_alb,grd)
      implicit none
      real(rkx), dimension(:,:,:), pointer, contiguous, intent(inout) :: alb
      real(rkx), dimension(:,:), intent(in) :: cpl_alb
      type(infogrd), intent(in) :: grd
      integer(ik4) :: ii, jj, nn
      call fill_land(tmp_lnd,cpl_alb,mddom%lndcat,grd)
      do concurrent ( nn = 1:nnsg, jj = jci1:jci2, ii = ici1:ici2 )
        if ( mdsub%ldmsk(nn,jj,ii) == 1 ) then
          if ( tmp_lnd(nn,jj,ii) > 0.9999_rkx ) then
            alb(nn,jj,ii) = 0.16_rkx
          else
            alb(nn,jj,ii) = tmp_lnd(nn,jj,ii)
          end if
        end if
      end do
    end subroutine fill_albedo
#endif
  end subroutine oasisxregcm_rcv_all

  ! call all subroutines consisting of sending OASIS fields
  ! with optional prior reworking
  subroutine oasisxregcm_snd_all(time)
    implicit none
    integer(ik4), intent(in) :: time ! execution time
    type(infogrd), pointer :: grd
    integer(ik4) :: i, j, ishift, jshift
    logical :: l_write_restart
#ifdef ECLM
    real(rkx), dimension(:,:), allocatable :: tmp1, tmp2
#endif
    !--------------------------------------------------------------------------
    l_write_restart = .false.
    if ( associated(alarm_out_sav) ) then
      if ( alarm_out_sav%act( ) ) then
        l_write_restart = .true.
      else
        if ( lfdomonth(rcmtimer%idate) .and. &
             lmidnight(rcmtimer%idate) ) then
          l_write_restart = .true.
        end if
      end if
    end if
    if ( rcmtimer%reached_endtime ) l_write_restart = .true.
    if ( write_restart_option == 1 .and. time == 0 ) l_write_restart = .true.
    if ( l_cpl_ex_u10m ) then ! eatward near-surface wind (10m-height) [m.s-1]
      grd => ex_u10m%grd
      call oasisxregcm_snd( &
           sfs%u10m(grd%j1:grd%j2, grd%i1:grd%i2), &
           ex_u10m, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_v10m ) then ! northward near-surface wind (10m-height) [m.s-1]
      grd => ex_v10m%grd
      call oasisxregcm_snd( &
           sfs%v10m(grd%j1:grd%j2, grd%i1:grd%i2), &
           ex_v10m, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_wspd ) then ! near-surface wind speed (10m-height) [m.s-1]
      grd => ex_wspd%grd
      call oasisxregcm_snd( &
           sqrt(sfs%u10m(grd%j1:grd%j2, grd%i1:grd%i2)**2   &
              + sfs%v10m(grd%j1:grd%j2, grd%i1:grd%i2)**2), &
           ex_wspd, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    ! near-surface wind direction (10m-height) [degree]
    if ( l_cpl_ex_wdir ) then
      grd => ex_wdir%grd
      do i = grd%i1, grd%i2
        ishift = i - grd%i1 + 1
        do j = grd%j1, grd%j2
          jshift = j - grd%j1 + 1
          cpl_wdir(jshift,ishift) = &
                  atan2( sfs%u10m(j,i), sfs%v10m(j,i) ) * raddeg
          if ( cpl_wdir(jshift,ishift) < d_zero ) &
               cpl_wdir(jshift,ishift) = cpl_wdir(jshift,ishift) + 360_rkx
        end do
      end do
      call oasisxregcm_snd(cpl_wdir, ex_wdir, time, &
              .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_t2m ) then ! near-surface air temperature (2m-height) [K]
      grd => ex_t2m%grd
      call oasisxregcm_snd( &
           sum( lms%t2m(:, grd%j1:grd%j2, grd%i1:grd%i2), 1 ) * rdnnsg, &
           ex_t2m, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
!    if ( l_cpl_ex_t10m ) then ! near-surface air temperature (10m-height) [K]
!      grd => ex_t10m%grd
!      call oasisxregcm_snd( &
!           ! ???
!           ex_t10m, time, .false. .or. l_write_restart)
!      nullify(grd)
!    end if
    !
    ! near-surface air specific humidity (2m-height) [1]
    if ( l_cpl_ex_q2m ) then
      grd => ex_q2m%grd
      call oasisxregcm_snd( &
           sfs%q2m(grd%j1:grd%j2, grd%i1:grd%i2), &
           ex_q2m, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
!    ! near-surface air specific humidity (10m-height) [1]
!    if ( l_cpl_ex_q10m ) then
!      grd => ex_q10m%grd
!      call oasisxregcm_snd( &
!           ! ???
!           ex_q10m, time, .false. .or. l_write_restart)
!      nullify(grd)
!    end if
    !
    if ( l_cpl_ex_slp ) then ! sea level pressure [Pa]
      grd => ex_slp%grd
      call oasisxregcm_snd( &
           atms%ps2d(grd%j1:grd%j2, grd%i1:grd%i2), &
           ex_slp, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_taux) then ! surface eastward wind stress [Pa]
      grd => ex_taux%grd
      call oasisxregcm_snd( &
           sum( lms%taux(:, grd%j1:grd%j2, grd%i1:grd%i2), 1 ) * rdnnsg, &
           ex_taux, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_tauy) then ! surface northward wind stress [Pa]
      grd => ex_tauy%grd
      call oasisxregcm_snd( &
           sum( lms%tauy(:, grd%j1:grd%j2, grd%i1:grd%i2), 1 ) * rdnnsg, &
           ex_tauy, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_z0 ) then ! surface roughness length [m]
      grd => ex_z0%grd
      call oasisxregcm_snd( &
           sfs%zo(grd%j1:grd%j2, grd%i1:grd%i2), &
           ex_z0, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_ustr ) then ! surface friction velocity [s-1]
      grd => ex_ustr%grd
      call oasisxregcm_snd( &
           sfs%ustar(grd%j1:grd%j2, grd%i1:grd%i2), &
           ex_ustr, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_evap ) then ! water evaporation flux [kg.m-2.s-1]
      grd => ex_evap%grd
      call oasisxregcm_snd( &
           sfs%qfx(grd%j1:grd%j2, grd%i1:grd%i2), &
           ex_evap, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_prec ) then ! precipitation flux [kg.m-2.s-1]
      grd => ex_prec%grd
      call oasisxregcm_snd( &
           sum( lms%prcp(:, grd%j1:grd%j2, grd%i1:grd%i2), 1 ) * rdnnsg, &
           ex_prec, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_nuwa ) then ! net upward water flux [kg.m-2.s-1]
      grd => ex_nuwa%grd
      call oasisxregcm_snd( &
           sum( lms%evpr(:, grd%j1:grd%j2, grd%i1:grd%i2) - &
                lms%prcp(:, grd%j1:grd%j2, grd%i1:grd%i2), 1 ) * rdnnsg, &
           ex_nuwa, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_ulhf ) then ! surface upward latent heat flux [W.m-2]
      grd => ex_ulhf%grd
      call oasisxregcm_snd( &
           sum( lms%evpr(:, grd%j1:grd%j2, grd%i1:grd%i2) &
           *wlh(lms%t2m(:, grd%j1:grd%j2, grd%i1:grd%i2)), 1 ) * rdnnsg, &
           ex_ulhf, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_ushf ) then ! surface upward sensible heat flux [W.m-2]
      grd => ex_ushf%grd
      call oasisxregcm_snd( &
           sfs%hfx(grd%j1:grd%j2, grd%i1:grd%i2), &
           ex_ushf, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    ! surface upwelling long-wave radiation flux [W.m-2]
    if ( l_cpl_ex_uwlw ) then
      grd => ex_uwlw%grd
      call oasisxregcm_snd( &
           flwd(grd%j1:grd%j2, grd%i1:grd%i2) &
           + flw(grd%j1:grd%j2, grd%i1:grd%i2), &
           ex_uwlw, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    ! surface downwelling long-wave radiation flux [W.m-2]
    if ( l_cpl_ex_dwlw ) then
      grd => ex_dwlw%grd
      call oasisxregcm_snd( &
           flwd(grd%j1:grd%j2, grd%i1:grd%i2), &
           ex_dwlw, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    ! surface net upward long-wave radiation flux [W.m-2]
    if ( l_cpl_ex_nulw ) then
      grd => ex_nulw%grd
      call oasisxregcm_snd( &
           flw(grd%j1:grd%j2, grd%i1:grd%i2), &
           ex_nulw, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    ! surface upwelling short-wave radiation flux [W.m-2]
    if ( l_cpl_ex_uwsw ) then
      grd => ex_uwsw%grd
      call oasisxregcm_snd( &
           sinc(grd%j1:grd%j2, grd%i1:grd%i2) &
           - fsw(grd%j1:grd%j2, grd%i1:grd%i2), &
           ex_uwsw, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    ! surface downwelling short-wave radiation flux [W.m-2]
    if ( l_cpl_ex_dwsw ) then
      grd => ex_dwsw%grd
      call oasisxregcm_snd( &
           sinc(grd%j1:grd%j2, grd%i1:grd%i2), &
           ex_dwsw, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    ! surface net downward short-wave radiation flux [W.m-2]
    if ( l_cpl_ex_ndsw ) then
      grd => ex_ndsw%grd
      call oasisxregcm_snd( &
           fsw(grd%j1:grd%j2, grd%i1:grd%i2), &
           ex_ndsw, time, .false. .or. l_write_restart)
      nullify(grd)
    end if
    !
    if ( l_cpl_ex_rhoa ) then ! surface air density [kg.m-3]
      grd => ex_rhoa%grd
      call oasisxregcm_snd( &
           sum( lms%rhoa(:, grd%j1:grd%j2, grd%i1:grd%i2), 1 ) * rdnnsg, &
           ex_rhoa, time, .false. .or. l_write_restart)
      nullify(grd)
    end if

#ifdef ECLM
    ! Atmospheric forcing of eCLM: lowest model level and surface fluxes
    if ( l_cpl_ex_tatm ) then ! air temperature [K]
      grd => ex_tatm%grd
      call oasisxregcm_snd(atms%tb3d(grd%j1:grd%j2,grd%i1:grd%i2,kz), &
           ex_tatm, time, l_write_restart)
      nullify(grd)
    end if
    if ( l_cpl_ex_uatm ) then ! grid relative zonal wind [m/s]
      grd => ex_uatm%grd
      call oasisxregcm_snd(atms%ubx3d(grd%j1:grd%j2,grd%i1:grd%i2,kz), &
           ex_uatm, time, l_write_restart)
      nullify(grd)
    end if
    if ( l_cpl_ex_vatm ) then ! grid relative meridional wind [m/s]
      grd => ex_vatm%grd
      call oasisxregcm_snd(atms%vbx3d(grd%j1:grd%j2,grd%i1:grd%i2,kz), &
           ex_vatm, time, l_write_restart)
      nullify(grd)
    end if
    if ( l_cpl_ex_qatm ) then ! specific humidity [kg/kg]
      grd => ex_qatm%grd
      call oasisxregcm_snd( &
           atms%qxb3d(grd%j1:grd%j2,grd%i1:grd%i2,kz,iqv) / &
           (d_one + atms%qxb3d(grd%j1:grd%j2,grd%i1:grd%i2,kz,iqv)), &
           ex_qatm, time, l_write_restart)
      nullify(grd)
    end if
    if ( l_cpl_ex_hgt ) then ! height of the lowest level [m]
      grd => ex_hgt%grd
      call oasisxregcm_snd(atms%za(grd%j1:grd%j2,grd%i1:grd%i2,kz), &
           ex_hgt, time, l_write_restart)
      nullify(grd)
    end if
    if ( l_cpl_ex_patm ) then ! pressure [Pa]
      grd => ex_patm%grd
      call oasisxregcm_snd(atms%pb3d(grd%j1:grd%j2,grd%i1:grd%i2,kz), &
           ex_patm, time, l_write_restart)
      nullify(grd)
    end if
    if ( l_cpl_ex_dwrlwf ) then ! downward longwave [W/m2]
      grd => ex_dwrlwf%grd
      call oasisxregcm_snd(flwd(grd%j1:grd%j2,grd%i1:grd%i2), &
           ex_dwrlwf, time, l_write_restart)
      nullify(grd)
    end if
    if ( l_cpl_ex_swdir ) then ! direct shortwave, visible + nir [W/m2]
      grd => ex_swdir%grd
      call oasisxregcm_snd(solvs(grd%j1:grd%j2,grd%i1:grd%i2) + &
                           solvl(grd%j1:grd%j2,grd%i1:grd%i2), &
           ex_swdir, time, l_write_restart)
      nullify(grd)
    end if
    if ( l_cpl_ex_swdif ) then ! diffuse shortwave, visible + nir [W/m2]
      grd => ex_swdif%grd
      call oasisxregcm_snd(solvsd(grd%j1:grd%j2,grd%i1:grd%i2) + &
                           solvld(grd%j1:grd%j2,grd%i1:grd%i2), &
           ex_swdif, time, l_write_restart)
      nullify(grd)
    end if
    ! Precipitation [kg/m2/s], snow below 1 degree C at the lowest level
    if ( l_cpl_ex_rainf .or. l_cpl_ex_snow ) then
      if ( l_cpl_ex_rainf ) then
        grd => ex_rainf%grd
      else
        grd => ex_snow%grd
      end if
      allocate(tmp1(grd%j1:grd%j2,grd%i1:grd%i2))
      allocate(tmp2(grd%j1:grd%j2,grd%i1:grd%i2))
      do i = grd%i1, grd%i2
        do j = grd%j1, grd%j2
          if ( atms%tb3d(j,i,kz) < tzero + d_one ) then
            tmp1(j,i) = d_zero
            tmp2(j,i) = sum(lms%prcp(:,j,i))*rdnnsg
          else
            tmp1(j,i) = sum(lms%prcp(:,j,i))*rdnnsg
            tmp2(j,i) = d_zero
          end if
        end do
      end do
      if ( l_cpl_ex_rainf ) then
        call oasisxregcm_snd(tmp1, ex_rainf, time, l_write_restart)
      end if
      if ( l_cpl_ex_snow ) then
        call oasisxregcm_snd(tmp2, ex_snow, time, l_write_restart)
      end if
      deallocate(tmp1, tmp2)
      nullify(grd)
    end if
#endif
    ! OASIS field +++
    !
    if ( myid == italk ) then
      if (write_restart_option == 1 .and. time == 0) then
        write(stdout,*) 'Note: OASIS restart files written at the' &
                      //' first time.'
      end if
    end if

    contains

#include <wlh.inc>

  end subroutine oasisxregcm_snd_all

  ! wait during the time indicated through oasis_sync_lag
  ! to make RegCM model time not synchronised with
  ! other coupled components
  ! update oasis_lag such that
  !   OASIS time = RegCM time + oasis_lag
  subroutine oasisxregcm_sync_wait(time)
    implicit none
    integer(ik4), intent(in) :: time ! execution time
    !--------------------------------------------------------------------------
    if ( oasis_lag /= 0 ) then
      write(stderr,*) 'error initiating the lag with OASIS'
      call fatal(__FILE__,__LINE__,'SYNC WAIT')
    end if
    ! case number 1: oasis_sync_lag > 0
    ! RegCM actually starts oasis_sync_lag seconds after OASIS
    ! >> at the beginning of the run, OASIS starts while RegCM runs
    ! dummy loops with oasis_lag increasing like 0 -> oasis_sync_lag
    ! at the end of the run, OASIS and RegCM stop synchronously
    if ( oasis_sync_lag > 0 ) then
      do while ( oasis_lag < oasis_sync_lag )
        call oasisxregcm_rcv_all(time + oasis_lag)
        ! dummy loop
#ifdef DEBUG
        if ( oasis_lag == 0 ) then
          write(ndebug,*) oasis_prefix, 'RegCM is creating the lag ', &
                                        'with the other components'
          write(ndebug,*) oasis_prefix, ' oasis_lag: ', oasis_lag
        end if
#endif
        call oasisxregcm_snd_all(time + oasis_lag)
        oasis_lag = oasis_lag + int(dtsec,ik4)
#ifdef DEBUG
        write(ndebug,*) oasis_prefix, ' oasis_lag: ', oasis_lag
#endif
      end do
    ! case number 2: oasis_sync_lag < 0
    ! RegCM starts synchronously with OASIS
    ! however, some other components start with a lag
    ! at the end of the run, OASIS needs oasis_sync_lag seconds to run
    ! with the other components while RegCM has already finished
    ! >> at the end, OASIS keeps going while RegCM runs dummy loops
    ! with oasis_lag increasing like end_time + dt -> -oasis_sync_lag
    else if ( oasis_sync_lag < 0 ) then
#ifdef DEBUG
      write(ndebug,*) oasis_prefix, 'RegCM is catching up ', &
                                    'with the other components'
      write(ndebug,*) oasis_prefix, ' oasis_lag: ', oasis_lag
#endif
      do while ( oasis_lag < -oasis_sync_lag )
        call oasisxregcm_rcv_all(time + int(dtsec,ik4) + oasis_lag)
        ! dummy loop
        call oasisxregcm_snd_all(time + int(dtsec,ik4) + oasis_lag)
        oasis_lag = oasis_lag + int(dtsec,ik4)
#ifdef DEBUG
        write(ndebug,*) oasis_prefix, ' oasis_lag: ', oasis_lag
#endif
      end do
    end if
  end subroutine oasisxregcm_sync_wait

  ! call all subroutines linked with some deallocation of variables used in
  ! the OASIS coupling
  subroutine oasisxregcm_release
    implicit none
    !--------------------------------------------------------------------------
    call oasisxregcm_deallocate_field(im_sst, cpl_sst)
    call oasisxregcm_deallocate_field(im_wz0, cpl_wz0)
    call oasisxregcm_deallocate_field(im_wust, cpl_wust)
    call oasisxregcm_deallocate_field(ex_u10m)
    call oasisxregcm_deallocate_field(ex_v10m)
    call oasisxregcm_deallocate_field(ex_wspd)
    call oasisxregcm_deallocate_field(ex_wdir, cpl_wdir)
    call oasisxregcm_deallocate_field(ex_t2m)
    call oasisxregcm_deallocate_field(ex_q2m)
    call oasisxregcm_deallocate_field(ex_slp)
    call oasisxregcm_deallocate_field(ex_taux)
    call oasisxregcm_deallocate_field(ex_tauy)
    call oasisxregcm_deallocate_field(ex_z0)
    call oasisxregcm_deallocate_field(ex_ustr)
    call oasisxregcm_deallocate_field(ex_evap)
    call oasisxregcm_deallocate_field(ex_prec)
    call oasisxregcm_deallocate_field(ex_nuwa)
    call oasisxregcm_deallocate_field(ex_ulhf)
    call oasisxregcm_deallocate_field(ex_ushf)
    call oasisxregcm_deallocate_field(ex_uwlw)
    call oasisxregcm_deallocate_field(ex_dwlw)
    call oasisxregcm_deallocate_field(ex_nulw)
    call oasisxregcm_deallocate_field(ex_uwsw)
    call oasisxregcm_deallocate_field(ex_dwsw)
    call oasisxregcm_deallocate_field(ex_ndsw)
    call oasisxregcm_deallocate_field(ex_rhoa)
!   call oasisxregcm_deallocate_field(im_sit, cpl_sit)
!   call oasisxregcm_deallocate_field(ex_t10m)
!   call oasisxregcm_deallocate_field(ex_q10m)
#ifdef ECLM
    call oasisxregcm_deallocate_field(ex_tatm)
    call oasisxregcm_deallocate_field(ex_uatm)
    call oasisxregcm_deallocate_field(ex_vatm)
    call oasisxregcm_deallocate_field(ex_qatm)
    call oasisxregcm_deallocate_field(ex_hgt)
    call oasisxregcm_deallocate_field(ex_patm)
    call oasisxregcm_deallocate_field(ex_dwrlwf)
    call oasisxregcm_deallocate_field(ex_swdir)
    call oasisxregcm_deallocate_field(ex_swdif)
    call oasisxregcm_deallocate_field(ex_rainf)
    call oasisxregcm_deallocate_field(ex_snow)
    call oasisxregcm_deallocate_field(im_t2m, cpl_t2m)
    call oasisxregcm_deallocate_field(im_q2m, cpl_q2m)
    call oasisxregcm_deallocate_field(im_u10m, cpl_u10m)
    call oasisxregcm_deallocate_field(im_sent, cpl_sent)
    call oasisxregcm_deallocate_field(im_evpr, cpl_evpr)
    call oasisxregcm_deallocate_field(im_ram1, cpl_ram1)
    call oasisxregcm_deallocate_field(im_rah1, cpl_rah1)
    call oasisxregcm_deallocate_field(im_taux, cpl_taux)
    call oasisxregcm_deallocate_field(im_tauy, cpl_tauy)
    call oasisxregcm_deallocate_field(im_sncv, cpl_sncv)
    call oasisxregcm_deallocate_field(im_wt, cpl_wt)
    call oasisxregcm_deallocate_field(im_zo, cpl_zo)
    call oasisxregcm_deallocate_field(im_tlef, cpl_tlef)
    call oasisxregcm_deallocate_field(im_albei, cpl_albei)
    call oasisxregcm_deallocate_field(im_albed, cpl_albed)
    call oasisxregcm_deallocate_field(im_albei_2, cpl_albei_2)
    call oasisxregcm_deallocate_field(im_albed_2, cpl_albed_2)
    call oasisxregcm_deallocate_field(im_fluxch4, cpl_fluxch4)
    call oasisxregcm_deallocate_field(im_tlai, cpl_tlai)
    call oasisxregcm_deallocate_field(im_roff, cpl_roff)
    call oasisxregcm_deallocate_field(im_srfroff, cpl_srfroff)
    call oasisxregcm_deallocate_field(im_snwmelt, cpl_snwmelt)
    call oasisxregcm_deallocate_field(im_tgrd, cpl_tgrd)
    call oasisxregcm_deallocate_field(im_tgbb, cpl_tgbb)
    call oasisxregcm_deallocate_field(im_h2o10cm, cpl_h2o10cm)
    call oasisxregcm_deallocate_3dfield(im_flxvoc, cpl_flxvoc)
    call oasisxregcm_deallocate_3dfield(im_flxdst, cpl_flxdst)
    call oasisxregcm_deallocate_3dfield(im_ddvel, cpl_ddvel)
    call oasisxregcm_deallocate_3dfield(im_tsoi, cpl_tsoi)
    call oasisxregcm_deallocate_3dfield(im_h2ovol, cpl_h2ovol)
    call oasisxregcm_deallocate_3dfield(im_h2oliq, cpl_h2oliq)
    call oasisxregcm_deallocate_3dfield(im_h2oice, cpl_h2oice)
#endif
    ! OASIS field +++
    if ( allocated(grdde) ) deallocate(grdde)
    if ( allocated(grddi) ) deallocate(grddi)
    if ( allocated(grdce) ) deallocate(grdce)
    if ( allocated(grdci) ) deallocate(grdci)
#ifdef ECLM
    if ( allocated(grdlnd) ) deallocate(grdlnd)
#endif
  end subroutine oasisxregcm_release

end module mod_oasis_interface
!
#endif
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

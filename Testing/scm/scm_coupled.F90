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

! Single column test of the Holtslag boundary layer (ibltyp = 1), the
! Tiedtke (icup = 5, default) or Kain-Fritsch (icup = 6, icumulus = 6)
! cumulus, the NoTo microphysics (ipptls = 2) and the CCM radiation
! (irrtm = 0) schemes integrated together in time.
!
! Usage: scm_coupled [namelist file]       (default: scm_coupled.in)
!
! The columns of the profile file are laid out as a 1 x ncol grid. Every
! time step the schemes are called through their model entry points,
! tiedtkedrv or kfdrv, subex_cldfrac (icldfrac = 0, or xuran_cldfrac with
! icldfrac = 1, on the cloud water and ice) and nogtom, colmod3
! (every dtrad), holtbl, on the same state, as the MOLOCH core does (idynamic
! = 3, no chemistry), and the state is advanced with the sum of their
! tendencies and of the large scale forcing:
!   - the input tendencies of the profile file (tten, qten: radiative
!     heating and horizontal advection), held fixed; with the model
!     radiation (irad = 1, default) tten is replaced by its heating,
!     or, with ladvt = .true., tten is a large scale advection added to
!     the model radiation (FIRE);
!   - the vertical advection by the omega of the profile file, computed
!     on the current state, which is also the omega seen by the cumulus
!     scheme;
!   - a relaxation of temperature and humidity with time scale taurelax
!     (default one day), above zrelax (linear ramp from zrelax/2; zero,
!     the default, relaxes all levels), to the initial profile or, if
!     tgtfile is given (default scm_amazon_target.txt), to the
!     profiles of that file (same format, for each column four profiles
!     at 0, 6, 12 and 18 UTC, linearly interpolated in time, repeated
!     every day).
! With tauwtg > 0 the weak temperature gradient approximation replaces
! the omega of the profile file: above zwtg and below 100 hPa the
! large scale omega is the one whose vertical advection brings the
! temperature back to the target in tauwtg, omega = (Ttarget - T) /
! (tauwtg S), with S = R T / (cp p) - dT/dp the static stability of the
! target profile (at least 1e-4 K/Pa),
! linear in height from zero at the ground to its value at zwtg; there
! the temperature is not relaxed. The convergence of the large scale
! flow (d omega / dp > 0) brings in air with the target humidity.
! The omega of the profile file is multiplied by omfac (default 1; 0
! for radiative convective equilibrium, with taurelax very large).
! The winds relax to the target winds (or the initial ones), times
! windfac, with time scale tauwind; with tauwind = 0 (default) the
! Coriolis force turns them around the initial winds times windfac, or
! with lgeos = .true. around the target winds times windfac, which are
! then the geostrophic wind (BOMEX, RICO).
! Land columns (and ocean columns with iocnflux = 0): the surface fluxes
! follow a half sine diurnal cycle in local solar time (daylight from 6
! to 18), with the daily means of the profile file: at night the
! sensible heat flux is -0.25 and the latent heat flux 0.1 times the
! daily mean, the daytime half sine keeps the mean; ocean columns with
! iocnflux = 0 keep the fluxes of the profile file constant. With
! ustarfix > 0 the friction velocity of the surface drag of the land
! columns and of the ocean columns with iocnflux = 0 is ustarfix (BOMEX)
! instead of the neutral log law one. Ocean columns with
! iocnflux = 1 (Zeng, with izengrough) or 2 (COARE 3.5): the fluxes,
! the drag and the friction velocity come from the ocean scheme of the
! model on the SST of the profile file, every step, with a net
! shortwave radiation at the surface of daily mean ocnsw in a half sine
! and a net longwave ocnlw (used by the COARE cool skin).
! Kain-Fritsch (icumulus = 6, with &kfparam): the vertical velocity of
! its trigger is the one of the large scale omega (-omega/(rho g)), and
! its tendencies of cloud water, ice, rain and snow go to the
! microphysics with the other cumulus tendencies, as in the model (no
! qdetr); the grid spacing is dxkm. The output columns of the Tiedtke
! mass flux and cloud type are zero (ktype: 1 when it rains).
! Microphysics (lnoto = .true., default): NoTo gets the cumulus
! tendencies, the condensate detrained by the cumulus scheme (qdetr),
! the large scale omega and the input temperature tendency as radiative
! heating; it carries cloud water, ice, rain and snow (not advected by
! the large scale flow) and gives the large scale rain. With lnoto =
! .false. the detrained condensate evaporates in the layer where it is
! detrained and supersaturation is removed by a saturation adjustment
! that rains out. The cloud water and ice of NoTo are also seen by the
! cumulus and boundary layer schemes, and diffused by the latter.
! Radiation (irad = 1): the CCM scheme of the model through colmod3,
! every dtrad (1800 s), with absorptivities and emissivities every
! dtabem (64800 s); sun position of the time of the step for the orbit
! of the year of the case; GHG of the built in table of scenario
! (RCP4.5); built in ozone profile; no aerosol; clouds as cucloud and
! cldfrac of the model (convective cloud at most cu_cfmax(dxkm) overlapped
! with the NoTo cloud); surface albedo of the model ocean (iwhitecap =
! 0) or albland on land, emissivity 0.98 (ocean) or 0.96 (land), ground
! temperature the SST or the one of the boundary layer input. Its net
! surface radiation replaces ocnsw and ocnlw for the ocean scheme, and
! its heating is the radiative heating of NoTo.
! The pressure of the levels is held fixed.
!
! The namelist has a &scmparam group and the model groups &tiedtkeparam,
! &kfparam,
! &cldparam (kfac_deep, kfac_shal, k2_const, iconvcf, rh0oce, rh0land,
! tc0, rhmin, rhmax, larcticcorr), &holtslagparam and &microparam, with
! the RegCM defaults. The output file is plain text, read by
! scm_coupled.py: a header with the settings used, then every nout
! steps per column a STEP line with the scalars and one line per level,
! top to bottom.

program scm_coupled
  use mod_intkinds, only : ik4, ik8
  use mod_realkinds, only : rk8, rkx
  use mod_constants, only : cpd, d_half, d_one, d_two, d_zero, degrad, egrav, &
    eomeg, ep1, ep2, lowcld, mathpi, p00, rcpd, regrav, rgas, rk8, rkx, &
    rovcp, rtice, rtwat, rtwat_rtice_r, rwat, tzero, vonkar, wlhs, wlhv
  use mod_dynparam, only : dayspy, dpd, half_dayspy, ice1, ice1ga, ice1gb, &
    ice2, ice2ga, ice2gb, ici1, ici1ga, ici2, ici2ga, idi1, idi2, idynamic, &
    ifopt, ifrad, jce1, jce1ga, jce1gb, jce2, jce2ga, jce2gb, jci1, jci1ga, &
    jci2, jci2ga, jdi1, jdi2, kz, kzm1, kzm2, kzp1, kzp2, kzp3, kzp4, ntr, &
    sixteenth_dayspy, vernal_equinox
  use mod_runparams, only : auto_rate_kessl, auto_rate_khair, &
    auto_rate_klepi, budget_compute, clfrcv, cu_cfmax, declin, detrpen_lnd, &
    detrpen_ocn, dt, dtabem, dtcum, dtrad, dtsec, dx, dxsq, eccen, eccf, &
    entrdd, entrpen_lnd, entrpen_ocn, entshalp, holtth10iter, iaerosol, &
    iautoconv, ibltyp, ichem, iclimaaer, iclimao3, icldfrac, iconvcf, &
    icumcloud, icup, &
    icup_lnd, icup_ocn, idirect, ifaholt, ifaholtth10, iindirect, iocnflux, &
    ipptls, iqc, iqi, iqr, iqs, iqv, irrtm, istochastic, itweak, iwavcpl, &
    izengrough, k2_const, kf_convrate, kf_dpp, kf_entrate, kf_max_dtcape, &
    kf_max_pef, kf_min_dtcape, kf_min_pef, kf_tkemax, kf_wthreshold, &
    kfac_deep, kfac_shal, lambm0, larcticcorr, lepcld, lmfcapebl, &
    lmfcapeforc, lmfdd, lmfdudv, lmfmid, lmfpen, lmfscv, lmfsmooth, lmftrac, &
    lmfuvdis, lmfwstar, lmoistn2, lqtmix, lsctop, lsrfhack, maxprandtl, &
    mvelp, &
    mvelpp, nqx, nssopt, obliq, obliqr, qxcheckval, qxzeroval, rcldiff, &
    rcmtimer, rcovpmin, rcpec_lnd, rcpec_ocn, rcuc_lnd, rcuc_ocn, rdt, &
    rh0land, rh0oce, rhebc_lnd, rhebc_ocn, rhmax, rhmin, ricr_lnd, ricr_ocn, &
    rkconv, rpecons, rprc_lnd, rprc_ocn, rtaucape, rtaufac, scenario, scon, &
    skconv, solcon, tc0, vfqs, zhnew_fac
  use mod_memutil, only : getmem, memory_init
  use mod_regcm_types, only : ik4, ik8, micro_2_mod, mod_2_cum, mod_2_micro, &
    mod_2_pbl, mod_2_rad, pbl_2_mod, rad_2_mod, rk8, rkx
  use mod_cu_common, only : avg_ww, cu_cldfrc, cu_cldlwc, cu_kbot, cu_ktop, &
    cu_prate, cu_qdetr, cu_qten, cu_raincc, cu_srate, cu_tten, cu_uten, &
    cu_vten, cuscheme, init_mod_cumulus, total_precip_points
  use mod_cu_tiedtke, only : allocate_mod_cu_tiedtke, tiedtkedrv, &
                             pmean, ktype, pmfu
  use mod_cu_kf, only : allocate_mod_cu_kf, kfdrv, kf_lutab
  use mod_pbl_common, only : ricr, kmxpbl
  use mod_pbl_holtbl, only : allocate_mod_pbl_holtbl, holtbl
  use mod_ocn_internal, only : nocnp, iocnbeg, iocnend, dtocn, llake, &
    ldcsst, lseaice, mask, icpl, ht, usw, vsw, tatm, qv, patm, sfps,  &
    sfta, tgb, tgrd, tgbrd, hpbl, rswf, rlwf, rhox, um10, sent, evpr, &
    drag, ustr, zoo, taux, tauy, u10m, v10m, t2m, q2m, rhoa, ram1,    &
    rah1, br, sfice
  use mod_ocn_zeng, only : zengocndrv
  use mod_ocn_coare, only : coare3_drv
  use mod_micro_nogtom, only : allocate_mod_nogtom, init_nogtom, nogtom
  use mod_cloud_subex, only : subex_cldfrac
  use mod_cloud_xuran, only : xuran_cldfrac
  use mod_rad_common, only : o3prof, gasabsnxt, gasabstot, gasemstot, &
    absemin, kth, ktf
  use mod_rad_aerosol, only : allocate_mod_rad_aerosol
  use mod_rad_o3blk, only : allocate_mod_rad_o3blk, o3data
  use mod_rad_outrad, only : allocate_mod_rad_outrad, curflx
  use mod_rad_colmod3, only : allocate_mod_rad_colmod3, colmod3
  use mod_ipcc_scenario, only : set_scenario
  use mod_sunorbit, only : orb_params, orb_decl, orb_cosz
  implicit none

  character(len=256) :: nmlfile, infile, outfile, tgtfile
  real(rkx) :: dxkm, dtstep, ndays, taurelax, zrelax, tauwtg, zwtg, &
    tauwind, windfac, ocnsw, ocnlw, omfac, ustarfix
  logical :: lgeos, ladvt
  integer(ik4) :: nout
  integer(ik4) :: ncol, nk, j, k, n, iu, io, ios, istep, nsteps
  character(len=32), allocatable, dimension(:) :: col_name
  character(len=10), allocatable, dimension(:) :: col_date
  real(rkx), allocatable, dimension(:) :: col_lat, col_lon, col_ps, &
    col_shf, col_lhf, col_zpbl, col_cape, col_crr, col_topo, col_z0,  &
    col_hour, rainc, rainl, rainc_acc, rainl_acc, ustar, col_sst
  integer(ik4), allocatable, dimension(:) :: col_land
  logical, allocatable, dimension(:) :: col_ocn
  real(rkx), allocatable, dimension(:,:) :: t0, q0, u0, v0, omega, &
    ften, fqen, za, zf, pa, pf, t, q, u, v, mfu, ctt, cqt, ptt, pqt, &
    wrelax, omls, ttg, qtg
  real(rkx), allocatable, dimension(:,:,:) :: ttgt, qtgt, utgt, vtgt
  real(rkx), allocatable, dimension(:,:) :: tadv
  real(rkx), allocatable, dimension(:,:) :: mtt, mqt, cfr, hrad, &
    cfrad, lwcrad
  integer(ik4) :: irad, iyear, imonth, nrad, nabs, icumulus, kty
  real(rkx), allocatable, dimension(:) :: ocnsw_now, ocnlw_now
  ! TOA reflected short wave, outgoing long wave, all sky and clear sky
  ! (W/m2) of the last radiation call
  real(rkx), allocatable, dimension(:) :: rsut, rlut, rsutcs, rlutcs
  real(rkx) :: albland
  real(rk8) :: cday0
  type(mod_2_rad) :: m2r
  type(rad_2_mod) :: r2m
  real(rkx), allocatable, dimension(:,:,:) :: qcnd
  real(rkx), pointer, contiguous, dimension(:,:,:) :: qvp, rhcp, totcp
  real(rkx), pointer, contiguous, dimension(:,:) :: rh0p
  logical :: lnoto
  type(mod_2_micro) :: mo2mc
  type(micro_2_mod) :: mc2mo
  type(mod_2_cum) :: m2c
  type(mod_2_pbl) :: m2p
  type(pbl_2_mod) :: p2m
  real(rkx) :: dp, hour, utc, wt, tg, qg, sfac, rho, kbfs, hqfac, &
    zpblj, qs, dqs, cond, ws, zl, dtdp, dqdp, time, ug, vg, stab, &
    dwdp, wspd, wr
  integer(ik4) :: kw
  character(len=512) :: line

  namelist /scmparam/ infile, outfile, dtstep, ndays, nout, dxkm, &
    taurelax, zrelax, tgtfile, tauwtg, zwtg, tauwind, windfac,    &
    omfac, lgeos, ustarfix, ladvt, &
    iocnflux, izengrough, ocnsw, ocnlw, lnoto, irad, dtrad, dtabem, &
    albland, scenario, icumulus, icldfrac
  namelist /tiedtkeparam/ entrdd, entrpen_lnd, entrpen_ocn, &
    detrpen_lnd, detrpen_ocn, entshalp, rcuc_lnd, rcuc_ocn, &
    rcpec_lnd, rcpec_ocn, rhebc_lnd, rhebc_ocn, rprc_ocn, rprc_lnd, &
    lmfpen, lmfmid, lmfdd, lepcld, lmfdudv, &
    lmfscv, lmfuvdis, lmftrac, lmfsmooth, lmfwstar, lmfcapebl, &
    lmfcapeforc, rtaucape, rtaufac
  namelist /kfparam/ kf_min_pef, kf_max_pef, kf_entrate, kf_dpp, &
    kf_min_dtcape, kf_max_dtcape, kf_tkemax, kf_convrate, kf_wthreshold
  namelist /cldparam/ kfac_deep, kfac_shal, k2_const, iconvcf, &
    rh0oce, rh0land, tc0, rhmin, rhmax, larcticcorr
  namelist /microparam/ nssopt, iautoconv, vfqs, auto_rate_khair, &
    auto_rate_kessl, auto_rate_klepi, rkconv, skconv, rcovpmin,   &
    rpecons, rcldiff
  namelist /holtslagparam/ ricr_ocn, ricr_lnd, zhnew_fac, &
    ifaholtth10, ifaholt, holtth10iter, lmoistn2, maxprandtl, lsctop, &
    lqtmix

  ! Defaults: as in Main/mod_params.F90
  infile = 'scm_amazon.txt'
  outfile = 'scm_coupled.out'
  tgtfile = 'scm_amazon_target.txt'
  dtstep = 300.0_rkx
  ndays = 5.0_rkx
  nout = 12
  dxkm = 12.0_rkx
  taurelax = 86400.0_rkx
  zrelax = 0.0_rkx
  tauwtg = 0.0_rkx
  zwtg = 2000.0_rkx
  tauwind = 0.0_rkx
  omfac = 1.0_rkx
  windfac = 1.0_rkx
  lgeos = .false.
  ladvt = .false.
  ustarfix = 0.0_rkx
  iocnflux = 1
  izengrough = 1
  ocnsw = 200.0_rkx
  ocnlw = 50.0_rkx
  entrpen_lnd = 0.5e-3_rkx
  entrpen_ocn = 0.5e-3_rkx
  detrpen_lnd = 0.75e-4_rkx
  detrpen_ocn = 0.75e-4_rkx
  entrdd = 3.0e-4_rkx
  entshalp = 2.0_rkx
  rcuc_lnd = 0.05_rkx
  rcuc_ocn = 0.05_rkx
  rcpec_lnd = 5.55e-5_rkx
  rcpec_ocn = 5.55e-5_rkx
  rhebc_lnd = 0.8_rkx
  rhebc_ocn = 0.8_rkx
  rprc_lnd = 1.4e-3_rkx
  rprc_ocn = 1.4e-3_rkx
  lmfpen = .true.
  lmfmid = .true.
  lmfdd = .true.
  lepcld = .true.
  lmfdudv = .true.
  lmfscv = .true.
  lmfuvdis = .true.
  lmftrac = .true.
  lmfsmooth = .false.
  lmfwstar = .false.
  lmfcapebl = .true.
  lmfcapeforc = .false.
  rtaucape = 43200.0_rkx
  rtaufac = 4.0_rkx
  k2_const = 500.0_rkx
  kfac_shal = 0.07_rkx
  kfac_deep = 0.14_rkx
  iconvcf = 1
  lnoto = .true.
  irad = 1
  dtrad = 1800.0_rkx
  dtabem = 64800.0_rkx
  albland = 0.13_rkx
  scenario = 'RCP4.5'
  rhmax = 1.00_rkx
  rhmin = 0.01_rkx
  rh0land = 0.85_rkx
  rh0oce = 0.85_rkx
  tc0 = 238.0_rkx
  larcticcorr = .true.
  nssopt = 1
  iautoconv = 4
  vfqs = 1.0_rkx
  auto_rate_khair = 0.355_rkx
  auto_rate_kessl = 1.e-3_rkx
  auto_rate_klepi = 0.5e-3_rkx
  rkconv = 0.833e-4_rkx
  skconv = 1.0e-3_rkx
  rcldiff = 1.0e-6_rkx
  rcovpmin = 0.1_rkx
  rpecons = 5.547e-5_rkx
  budget_compute = .false.
  ricr_ocn = 0.275_rkx
  ricr_lnd = 0.250_rkx
  zhnew_fac = 0.25_rkx
  ifaholt = 0
  ifaholtth10 = 2
  holtth10iter = 1
  lmoistn2 = .false.
  maxprandtl = 1.0_rkx
  lsctop = .true.
  lqtmix = .false.
  icumulus = 5
  icldfrac = 0
  kf_wthreshold = 0.02_rkx
  kf_entrate = 0.03_rkx
  kf_convrate = 0.03_rkx
  kf_min_pef = 0.2_rkx
  kf_max_pef = 0.9_rkx
  kf_dpp = 150.0_rkx
  kf_tkemax = 5.0_rkx
  kf_min_dtcape = 1800.0_rkx
  kf_max_dtcape = 3600.0_rkx
  istochastic = 0

  nmlfile = 'scm_coupled.in'
  if ( command_argument_count( ) > 0 ) call get_command_argument(1,nmlfile)
  open(newunit=iu,file=nmlfile,status='old',action='read',iostat=ios)
  if ( ios == 0 ) then
    read(iu,nml=scmparam,iostat=ios)
    if ( ios > 0 ) call scm_error('reading scmparam in '//trim(nmlfile))
    rewind(iu)
    read(iu,nml=tiedtkeparam,iostat=ios)
    if ( ios > 0 ) call scm_error('reading tiedtkeparam in '//trim(nmlfile))
    rewind(iu)
    read(iu,nml=cldparam,iostat=ios)
    if ( ios > 0 ) call scm_error('reading cldparam in '//trim(nmlfile))
    rewind(iu)
    read(iu,nml=holtslagparam,iostat=ios)
    if ( ios > 0 ) call scm_error('reading holtslagparam in '//trim(nmlfile))
    rewind(iu)
    read(iu,nml=microparam,iostat=ios)
    if ( ios > 0 ) call scm_error('reading microparam in '//trim(nmlfile))
    rewind(iu)
    read(iu,nml=kfparam,iostat=ios)
    if ( ios > 0 ) call scm_error('reading kfparam in '//trim(nmlfile))
    close(iu)
  else if ( command_argument_count( ) > 0 ) then
    call scm_error('cannot open '//trim(nmlfile))
  end if
  nout = max(nout,1)
  if ( icumulus /= 5 .and. icumulus /= 6 ) &
    call scm_error('icumulus must be 5 (Tiedtke) or 6 (Kain-Fritsch)')

  ! Model configuration seen by the schemes
  idynamic = 3
  icup_lnd = icumulus
  icup_ocn = icumulus
  icup(:) = icumulus
  ibltyp = 1
  ipptls = 2
  nqx = 5
  ichem = 0
  iaerosol = 0
  iindirect = 0
  icumcloud = 1
  ntr = 0
  dx = dxkm * 1000.0_rkx
  dxsq = dx * dx
  clfrcv = cu_cfmax(dxkm)
  dt = dtstep
  rdt = d_one/dt
  dtcum = dtstep
  nsteps = nint(ndays*86400.0_rkx/dtstep)

  ! Read the column profiles
  open(newunit=iu,file=infile,status='old',action='read',iostat=ios)
  if ( ios /= 0 ) call scm_error('cannot open '//trim(infile))
  do
    read(iu,'(a)') line
    if ( line(1:1) /= '#' ) exit
  end do
  read(line,*) ncol, nk
  kz = nk
  kzp1 = nk + 1
  kzm1 = nk - 1
  kzm2 = nk - 2
  kzp2 = nk + 2
  kzp3 = nk + 3
  kzp4 = nk + 4
  call set_grid(ncol)
  call memory_init
  allocate(col_name(ncol), col_date(ncol), col_lat(ncol), col_lon(ncol), &
    col_ps(ncol), col_shf(ncol), col_lhf(ncol), col_zpbl(ncol),           &
    col_cape(ncol), col_crr(ncol), col_topo(ncol), col_land(ncol),        &
    col_z0(ncol), col_hour(ncol), rainc(ncol), rainl(ncol),               &
    rainc_acc(ncol), rainl_acc(ncol), ustar(ncol), col_sst(ncol),         &
    col_ocn(ncol))
  allocate(t0(ncol,nk), q0(ncol,nk), u0(ncol,nk), v0(ncol,nk),       &
    omega(ncol,nk), ften(ncol,nk), fqen(ncol,nk), za(ncol,nk),       &
    zf(ncol,nk+1), pa(ncol,nk), pf(ncol,nk+1), t(ncol,nk),           &
    q(ncol,nk), u(ncol,nk), v(ncol,nk), mfu(ncol,nk), ctt(ncol,nk),  &
    cqt(ncol,nk), ptt(ncol,nk), pqt(ncol,nk), wrelax(ncol,nk),       &
    omls(ncol,nk), ttg(ncol,nk), qtg(ncol,nk))
  call getmem(cuscheme,jci1,jci2,ici1,ici2,'scm:cuscheme')
  cuscheme(:,:) = icumulus
  call init_mod_cumulus
  if ( icumulus == 5 ) then
    call allocate_mod_cu_tiedtke
  else
    call allocate_mod_cu_kf
    call kf_lutab
    ! arrays of Tiedtke the column code uses: they stay zero
    call getmem(cu_srate,jci1,jci2,ici1,ici2,'scm:srate')
    call getmem(cu_uten,jce1gb,jce2gb,ice1gb,ice2gb,1,kz,'scm:uten')
    call getmem(cu_vten,jce1gb,jce2gb,ice1gb,ice2gb,1,kz,'scm:vten')
    call getmem(cu_qdetr,jdi1,jdi2,idi1,idi2,1,kz,'scm:qdetr')
    call getmem(cu_raincc,jdi1,jdi2,idi1,idi2,1,kz,'scm:raincc')
  end if
  call getmem(ricr,jci1,jci2,ici1,ici2,'scm:ricr')
  call getmem(kmxpbl,jci1,jci2,ici1,ici2,'scm:kmxpbl')
  call allocate_mod_pbl_holtbl
  call allocate_fields
  do j = 1, ncol
    read(iu,*) col_name(j), col_lat(j), col_lon(j), col_land(j), &
      col_date(j), col_ps(j), col_shf(j), col_lhf(j), col_zpbl(j), &
      col_cape(j), col_crr(j), col_topo(j), col_sst(j)
    col_ocn(j) = col_land(j) /= 1 .and. col_sst(j) > d_zero .and. &
      iocnflux > 0
    do k = 1, nk
      read(iu,*) za(j,k), pa(j,k), t(j,k), q(j,k), u(j,k), v(j,k), &
        omega(j,k), ften(j,k), fqen(j,k)
      omega(j,k) = omfac*omega(j,k)
    end do
    do k = 1, nk + 1
      read(iu,*) zf(j,k), pf(j,k)
    end do
    ! specific humidity to mixing ratio, as in the model
    do k = 1, nk
      fqen(j,k) = fqen(j,k)/(d_one-q(j,k))**2
      q(j,k) = q(j,k)/(d_one-q(j,k))
    end do
    read(col_date(j)(9:10),*) k
    ! local solar time of the start
    col_hour(j) = modulo(real(k,rkx) + col_lon(j)/15.0_rkx,24.0_rkx)
    col_z0(j) = merge(0.1_rkx,1.0e-3_rkx,col_land(j) == 1)
    ricr(j,1) = merge(ricr_lnd,ricr_ocn,col_land(j) == 1)
    kmxpbl(j,1) = kz
    do k = kzm1, 2, -1
      if ( za(j,k) >= 5500.0_rkx ) exit
      kmxpbl(j,1) = k
    end do
    do k = 1, nk
      if ( zrelax > d_zero ) then
        wrelax(j,k) = min(max((za(j,k)-d_half*zrelax)/(d_half*zrelax), &
          d_zero),d_one)/taurelax
      else
        wrelax(j,k) = d_one/taurelax
      end if
    end do
    m2p%ldmsk(j,1) = col_land(j)
    m2c%ldmsk(j,1) = col_land(j)
    m2p%coriol(j,1) = d_two*eomeg*sin(col_lat(j)*degrad)
    m2c%psb(j,1) = col_ps(j)
  end do
  close(iu)
  t0 = t
  q0 = q
  u = windfac*u
  v = windfac*v
  u0 = u
  v0 = v
  allocate(ttgt(ncol,nk,0:4), qtgt(ncol,nk,0:4), utgt(ncol,nk,0:4), &
    vtgt(ncol,nk,0:4))
  do n = 0, 4
    ttgt(:,:,n) = t0
    qtgt(:,:,n) = q0
    utgt(:,:,n) = u0
    vtgt(:,:,n) = v0
  end do
  if ( len_trim(tgtfile) > 0 ) call read_targets
  if ( any(col_ocn) ) call init_ocean
  allocate(mtt(ncol,nk), mqt(ncol,nk), cfr(ncol,nk), qcnd(ncol,nk,iqc:iqs))
  allocate(hrad(ncol,nk), cfrad(ncol,nk), lwcrad(ncol,nk))
  allocate(ocnsw_now(ncol), ocnlw_now(ncol))
  allocate(rsut(ncol), rlut(ncol), rsutcs(ncol), rlutcs(ncol))
  rsut = d_zero
  rlut = d_zero
  rsutcs = d_zero
  rlutcs = d_zero
  ocnsw_now = -d_one
  ocnlw_now = ocnlw
  cfr = d_zero
  cfrad = d_zero
  lwcrad = d_zero
  hrad = ften
  ! input temperature tendency kept as advection with the model radiation
  allocate(tadv(ncol,nk))
  tadv = d_zero
  if ( irad == 1 .and. ladvt ) tadv = ften
  allocate(rcmtimer)
  rcmtimer%lcount = 0
  if ( irad == 1 ) call init_rad
  mtt = d_zero
  mqt = d_zero
  qcnd = d_zero
  if ( lnoto ) call init_micro
  if ( icumulus == 5 ) then
    do k = 1, nk
      pmean(k) = sum(pa(:,k))/real(ncol,rkx)
    end do
  end if
  do j = 1, ncol
    m2c%zas(j,1,:) = za(j,:)
    m2c%zfs(j,1,:) = zf(j,:)
    do k = 1, nk
      m2c%dzq(j,1,k) = zf(j,k) - zf(j,k+1)
    end do
    m2c%pas(j,1,:) = pa(j,:)
    m2c%pasf(j,1,:) = pf(j,:)
    m2p%za(j,1,:) = za(j,:)
    m2p%patm(j,1,:) = pa(j,:)
    m2p%patmf(j,1,:) = pf(j,:)
  end do
  p2m%zpbl(:,:) = 500.0_rkx
  rainc_acc = d_zero
  rainl_acc = d_zero

  open(newunit=io,file=outfile,status='replace',action='write')
  if ( icumulus == 5 ) then
    write(io,'(a)') '# RegCM Holtslag + Tiedtke coupled single column test'
  else
    write(io,'(a)') '# RegCM Holtslag + Kain-Fritsch coupled single '// &
      'column test'
  end if
  write(io,'(a,a)') '# profiles: ', trim(infile)
  write(io,'(a,a)') '# targets: ', trim(tgtfile)
  write(io,'(a,f8.2,a,f8.3,a,f10.1,a,f8.1,a,f8.3)') '# dt (s) ', dt, &
    '  ndays ', ndays, '  taurelax (s) ', taurelax, '  zrelax (m) ', &
    zrelax, '  dx (km) ', dxkm
  write(io,'(a,f10.1,a,f8.1,a,f10.1,a,f6.3)') '# tauwtg (s) ', tauwtg, &
    '  zwtg (m) ', zwtg, '  tauwind (s) ', tauwind, '  windfac ', windfac
  write(io,'(a,2i3,a,2f8.1)') '# iocnflux izengrough ', iocnflux, &
    izengrough, '  ocnsw ocnlw (W/m2) ', ocnsw, ocnlw
  if ( lgeos .or. ustarfix > d_zero .or. ladvt ) &
    write(io,'(a,l2,a,f8.3,a,l2)') '# lgeos ', lgeos, &
      '  ustarfix (m/s) ', ustarfix, '  ladvt ', ladvt
  if ( icumulus == 5 ) then
    write(io,nml=tiedtkeparam)
  else
    write(io,nml=kfparam)
  end if
  write(io,nml=cldparam)
  write(io,nml=holtslagparam)
  write(io,'(a,l2)') '# lnoto ', lnoto
  write(io,'(a,i2,a,2f10.1,a,f6.3,a,a)') '# irad ', irad, &
    '  dtrad dtabem (s) ', dtrad, dtabem, '  albland ', albland, &
    '  scenario ', scenario
  if ( lnoto ) write(io,nml=microparam)
  write(io,'(a)') '# STEP name time(h) local_hour shf lhf(W/m2) '// &
    'zpbl(m) kpbl hqfac ktype kctop kcbot rainc rainl(mm/day) '// &
    'rainc_acc rainl_acc(mm) wspd(m/s, lowest level) '// &
    'rsut rlut rsutcs rlutcs(W/m2, TOA, last radiation call)'
  write(io,'(a)') '# z(m) p(hPa) T(K) q(g/kg) mfu(kg/m2/s) '// &
    'dTcu dTpbl(K/day) dqcu dqpbl(g/kg/day) omega(Pa/s) '// &
    'dTmic(K/day) dqmic(g/kg/day) qc+qi(g/kg) cldfrac(radiation) '// &
    'dTrad(K/day) lwc(g/m3, in cloud, seen by radiation) '// &
    'cfls(large scale cloud fraction of the cloud scheme)'
  write(io,'(a,3i6)') 'NCOL ', ncol, nk, nsteps/nout

  do istep = 1, nsteps
    time = real(istep-1,rkx)*dt
    rcmtimer%lcount = int(istep,ik8)
    ! Surface fluxes
    do j = 1, ncol
      hour = modulo(col_hour(j) + time/3600.0_rkx,24.0_rkx)
      sfac = sin(mathpi*(hour-6.0_rkx)/12.0_rkx)
      if ( col_land(j) /= 1 ) then
        m2p%hfx(j,1) = col_shf(j)
        m2p%qfx(j,1) = col_lhf(j)/wlhv
      else if ( sfac > d_zero ) then
        m2p%hfx(j,1) = col_shf(j)*mathpi*1.125_rkx*sfac
        m2p%qfx(j,1) = col_lhf(j)*mathpi*0.95_rkx*sfac/wlhv
      else
        m2p%hfx(j,1) = -0.25_rkx*col_shf(j)
        m2p%qfx(j,1) = 0.1_rkx*col_lhf(j)/wlhv
      end if
    end do
    call update_derived
    if ( any(col_ocn) ) call ocean_fluxes

    ! Targets at this time, large scale omega
    do j = 1, ncol
      utc = modulo(col_hour(j) - col_lon(j)/15.0_rkx + &
        time/3600.0_rkx,24.0_rkx)
      n = min(int(utc/6.0_rkx),3)
      wt = (utc - 6.0_rkx*real(n,rkx))/6.0_rkx
      ttg(j,:) = (d_one-wt)*ttgt(j,:,n) + wt*ttgt(j,:,n+1)
      qtg(j,:) = (d_one-wt)*qtgt(j,:,n) + wt*qtgt(j,:,n+1)
    end do
    if ( tauwtg > d_zero ) then
      call wtg_omega
    else
      omls = omega
    end if

    ! Cumulus, on the state at the start of the step
    do j = 1, ncol
      m2c%tas(j,1,:) = t(j,:)
      m2c%qxas(j,1,:,iqv) = q(j,:)
      m2c%qxas(j,1,:,iqc) = qcnd(j,:,iqc)
      m2c%qxas(j,1,:,iqi) = qcnd(j,:,iqi)
      m2c%qxas(j,1,:,iqr) = qcnd(j,:,iqr)
      m2c%qxas(j,1,:,iqs) = qcnd(j,:,iqs)
      do k = 1, nk
        m2c%rhoas(j,1,k) = pa(j,k)/(rgas*t(j,k)* &
          (d_one+ep1*q(j,k)/(d_one+q(j,k))))
      end do
      m2c%uas(j,1,:) = u(j,:)
      m2c%vas(j,1,:) = v(j,:)
      m2c%hfx(j,1) = m2p%hfx(j,1)
      m2c%qfx(j,1) = m2p%qfx(j,1)
      m2c%zpbl(j,1) = p2m%zpbl(j,1)
      m2c%ustar(j,1) = ustar(j)
    end do
    m2c%tten(:,:,:) = d_zero
    m2c%qxten(:,:,:,:) = d_zero
    if ( icumulus == 5 ) then
      ! Tiedtke: the large scale omega (Pa/s)
      avg_ww(:,1,:) = omls
    else
      ! Kain-Fritsch: the large scale vertical velocity (m/s)
      avg_ww(:,1,:) = -omls/(m2c%rhoas(:,1,:)*egrav)
    end if
    cu_prate(:,:) = d_zero
    cu_srate(:,:) = d_zero
    cu_ktop(:,:) = 0
    cu_kbot(:,:) = 0
    cu_tten(:,:,:) = d_zero
    cu_qten(:,:,:,:) = d_zero
    cu_uten(:,:,:) = d_zero
    cu_vten(:,:,:) = d_zero
    cu_cldfrc(:,:,:) = d_zero
    cu_cldlwc(:,:,:) = d_zero
    cu_qdetr(:,:,:) = d_zero
    cu_raincc(:,:,:) = d_zero
    total_precip_points = 0
    if ( icumulus == 5 ) then
      call tiedtkedrv(m2c)
    else
      call kfdrv(m2c)
    end if

    ! Cloud fraction and microphysics, on the same state, with the
    ! cumulus tendencies as input (detrained condensate in qdetr)
    if ( lnoto ) call micro_step
    call rad_clouds

    ! Radiation, every dtrad, on the same state; the heating rate is
    ! held between calls
    if ( irad == 1 .and. mod(istep-1,nrad) == 0 ) call rad_step

    ! Boundary layer, on the same state
    p2m%tten(:,:,:) = d_zero
    p2m%uten(:,:,:) = d_zero
    p2m%vten(:,:,:) = d_zero
    p2m%qxten(:,:,:,:) = d_zero
    call holtbl(m2p,p2m)

    ! Advance the state
    do j = 1, ncol
      utc = modulo(col_hour(j) - col_lon(j)/15.0_rkx + &
        time/3600.0_rkx,24.0_rkx)
      n = min(int(utc/6.0_rkx),3)
      wt = (utc - 6.0_rkx*real(n,rkx))/6.0_rkx
      do k = 1, nk
        if ( lnoto ) then
          ctt(j,k) = cu_tten(j,1,k)
          cqt(j,k) = cu_qten(j,1,k,iqv)
        else
          ! detrained condensate evaporates in the layer
          ctt(j,k) = cu_tten(j,1,k) - &
            (wlhv*(cu_qten(j,1,k,iqc) + cu_qten(j,1,k,iqr)) + &
             wlhs*(cu_qten(j,1,k,iqi) + cu_qten(j,1,k,iqs)))*rcpd
          cqt(j,k) = cu_qten(j,1,k,iqv) + cu_qten(j,1,k,iqc) + &
            cu_qten(j,1,k,iqi) + cu_qten(j,1,k,iqr) + cu_qten(j,1,k,iqs)
        end if
        ptt(j,k) = p2m%tten(j,1,k)
        pqt(j,k) = p2m%qxten(j,1,k,iqv)
        if ( icumulus == 5 ) then
          mfu(j,k) = pmfu(j,k)
        else
          mfu(j,k) = d_zero
        end if
      end do
      do k = 1, nk
        ! vertical advection by omega, centred differences
        if ( k == 1 ) then
          dtdp = (t(j,2)-t(j,1))/(pa(j,2)-pa(j,1))
          dqdp = (q(j,2)-q(j,1))/(pa(j,2)-pa(j,1))
        else if ( k == nk ) then
          dtdp = (t(j,nk)-t(j,nk-1))/(pa(j,nk)-pa(j,nk-1))
          dqdp = (q(j,nk)-q(j,nk-1))/(pa(j,nk)-pa(j,nk-1))
        else
          dtdp = (t(j,k+1)-t(j,k-1))/(pa(j,k+1)-pa(j,k-1))
          dqdp = (q(j,k+1)-q(j,k-1))/(pa(j,k+1)-pa(j,k-1))
        end if
        tg = ttg(j,k)
        qg = qtg(j,k)
        ! temperature relaxation, not where the WTG omega acts
        wr = wrelax(j,k)
        if ( tauwtg > d_zero .and. za(j,k) >= zwtg .and. &
             pa(j,k) > 1.0e4_rkx ) wr = d_zero
        ! inflow of target air where the large scale flow converges
        dwdp = d_zero
        if ( tauwtg > d_zero ) then
          if ( k == nk ) then
            dwdp = (d_zero-omls(j,k))/(pf(j,nk+1)-pa(j,k))
          else
            dwdp = (omls(j,k+1)-omls(j,k))/(pa(j,k+1)-pa(j,k))
          end if
          dwdp = max(dwdp,d_zero)
        end if
        t(j,k) = t(j,k) + dt*(ctt(j,k) + mtt(j,k) + ptt(j,k) + &
          hrad(j,k) + tadv(j,k) - &
          omls(j,k)*(dtdp - rovcp*t(j,k)/pa(j,k)) - &
          wr*(t(j,k)-tg))
        q(j,k) = max(q(j,k) + dt*(cqt(j,k) + mqt(j,k) + pqt(j,k) + &
          fqen(j,k) - &
          omls(j,k)*dqdp - wrelax(j,k)*(q(j,k)-qg) + &
          dwdp*(qg-q(j,k))),1.0e-8_rkx)
        if ( tauwind > d_zero ) then
          ug = windfac*((d_one-wt)*utgt(j,k,n) + wt*utgt(j,k,n+1))
          vg = windfac*((d_one-wt)*vtgt(j,k,n) + wt*vtgt(j,k,n+1))
          u(j,k) = u(j,k) + dt*(cu_uten(j,1,k) + p2m%uten(j,1,k) - &
            (u(j,k)-ug)/tauwind)
          v(j,k) = v(j,k) + dt*(cu_vten(j,1,k) + p2m%vten(j,1,k) - &
            (v(j,k)-vg)/tauwind)
        else if ( lgeos ) then
          ug = windfac*((d_one-wt)*utgt(j,k,n) + wt*utgt(j,k,n+1))
          vg = windfac*((d_one-wt)*vtgt(j,k,n) + wt*vtgt(j,k,n+1))
          u(j,k) = u(j,k) + dt*(cu_uten(j,1,k) + p2m%uten(j,1,k) + &
            m2p%coriol(j,1)*(v(j,k)-vg))
          v(j,k) = v(j,k) + dt*(cu_vten(j,1,k) + p2m%vten(j,1,k) - &
            m2p%coriol(j,1)*(u(j,k)-ug))
        else
          u(j,k) = u(j,k) + dt*(cu_uten(j,1,k) + p2m%uten(j,1,k) + &
            m2p%coriol(j,1)*(v(j,k)-v0(j,k)))
          v(j,k) = v(j,k) + dt*(cu_vten(j,1,k) + p2m%vten(j,1,k) - &
            m2p%coriol(j,1)*(u(j,k)-u0(j,k)))
        end if
      end do
      if ( lnoto ) then
        ! condensate of the microphysics, no large scale advection
        do n = iqc, iqs
          do k = 1, nk
            ! as status_update of the MOLOCH core, with the vertical
            ! diffusion of cloud water and ice of the boundary layer
            qcnd(j,k,n) = qcnd(j,k,n) + dt*(mc2mo%qxten(j,1,k,n) + &
              p2m%qxten(j,1,k,n))
            if ( qcnd(j,k,n) < qxcheckval(n) ) qcnd(j,k,n) = qxzeroval(n)
          end do
        end do
        rainl(j) = mc2mo%trrate(j,1)
      else
        ! saturation adjustment, the condensate rains out
        rainl(j) = d_zero
        do k = 1, nk
          qs = pfwsat(t(j,k),pa(j,k))
          if ( q(j,k) > qs ) then
            dqs = pfqsdt(t(j,k),pa(j,k))
            cond = (q(j,k)-qs)/(d_one + wlhv*rcpd*dqs)
            q(j,k) = q(j,k) - cond
            t(j,k) = t(j,k) + wlhv*rcpd*cond
            rainl(j) = rainl(j) + cond*(pf(j,k+1)-pf(j,k))*regrav
          end if
        end do
        rainl(j) = rainl(j)/dt
      end if
      rainc(j) = cu_prate(j,1)
      rainc_acc(j) = rainc_acc(j) + rainc(j)*dt
      rainl_acc(j) = rainl_acc(j) + rainl(j)*dt
    end do

    if ( mod(istep,nout) == 0 ) then
      time = real(istep,rkx)*dt
      do j = 1, ncol
        hour = modulo(col_hour(j) + time/3600.0_rkx,24.0_rkx)
        rho = pf(j,nk+1)/(rgas*t(j,nk))
        kbfs = m2p%hfx(j,1)/(rho*cpd) + ep1*t(j,nk)*m2p%qfx(j,1)/rho
        zpblj = p2m%zpbl(j,1)
        if ( kbfs > d_zero ) then
          hqfac = max(0.8_rkx - 8.5_rkx*zhnew_fac*kbfs,0.25_rkx)
        else
          hqfac = 0.8_rkx
        end if
        wspd = sqrt(u(j,nk)**2+v(j,nk)**2)
        if ( icumulus == 5 ) then
          kty = ktype(j)
        else
          kty = merge(1,0,cu_ktop(j,1) > 0)
        end if
        write(io,'(a,1x,a,2f9.3,2f9.2,f9.1,i4,f7.3,3i4,4es13.5,f8.3,4f8.2)') &
          'STEP', trim(col_name(j)), time/3600.0_rkx, hour, &
          m2p%hfx(j,1), m2p%qfx(j,1)*wlhv, zpblj, p2m%kpbl(j,1), hqfac, &
          kty, cu_ktop(j,1), cu_kbot(j,1), rainc(j)*86400.0_rkx, &
          rainl(j)*86400.0_rkx, rainc_acc(j), rainl_acc(j), wspd, &
          rsut(j), rlut(j), rsutcs(j), rlutcs(j)
        do k = 1, nk
          write(io,'(f10.2,f10.3,f9.3,15es13.5)') za(j,k), &
            pa(j,k)*0.01_rkx, t(j,k), q(j,k)/(d_one+q(j,k))*1.0e3_rkx, &
            mfu(j,k), ctt(j,k)*86400.0_rkx, ptt(j,k)*86400.0_rkx, &
            cqt(j,k)*8.64e7_rkx, pqt(j,k)*8.64e7_rkx, omls(j,k), &
            mtt(j,k)*86400.0_rkx, mqt(j,k)*8.64e7_rkx, &
            (qcnd(j,k,iqc)+qcnd(j,k,iqi))*1.0e3_rkx, cfrad(j,k), &
            hrad(j,k)*86400.0_rkx, lwcrad(j,k), cfr(j,k)
        end do
      end do
    end if
  end do
  close(io)
  write(*,'(a,i4,a,i8,a,a)') ' scm_coupled: ', ncol, ' columns, ', &
    nsteps, ' steps, columns written to ', trim(outfile)

  contains

  subroutine scm_error(msg)
    implicit none
    character(len=*), intent(in) :: msg
    write(*,'(a,a)') ' scm_coupled error: ', msg
    stop 1
  end subroutine scm_error

  ! Relaxation targets: per column the profiles at 0, 6, 12, 18 UTC
  subroutine read_targets
    implicit none
    integer(ik4) :: iut, nc, nkt, m, kk
    real(rkx) :: x1, x2
    character(len=32) :: cname
    open(newunit=iut,file=tgtfile,status='old',action='read',iostat=ios)
    if ( ios /= 0 ) call scm_error('cannot open '//trim(tgtfile))
    do
      read(iut,'(a)') line
      if ( line(1:1) /= '#' ) exit
    end do
    read(line,*) nc, nkt
    if ( nc /= 4*ncol .or. nkt /= nk ) &
      call scm_error('target file does not match the profiles')
    do j = 1, ncol
      do m = 0, 3
        read(iut,*) cname
        do kk = 1, nk
          read(iut,*) x1, x2, ttgt(j,kk,m), qtgt(j,kk,m), &
            utgt(j,kk,m), vtgt(j,kk,m)
          qtgt(j,kk,m) = qtgt(j,kk,m)/(d_one-qtgt(j,kk,m))
        end do
        do kk = 1, nk + 1
          read(iut,*) x1, x2
        end do
      end do
      ttgt(j,:,4) = ttgt(j,:,0)
      qtgt(j,:,4) = qtgt(j,:,0)
      utgt(j,:,4) = utgt(j,:,0)
      vtgt(j,:,4) = vtgt(j,:,0)
    end do
    close(iut)
  end subroutine read_targets

  ! Weak temperature gradient omega: above zwtg (ramped to zero between
  ! 150 and 100 hPa) the omega whose vertical advection relaxes the
  ! temperature to the target in tauwtg, with the static stability of
  ! the target; below zwtg linear in height.
  subroutine wtg_omega
    implicit none
    real(rkx) :: ramp
    do j = 1, ncol
      omls(j,:) = d_zero
      kw = 0
      do k = 1, nk
        if ( za(j,k) < zwtg .or. pa(j,k) <= 1.0e4_rkx ) cycle
        if ( k == 1 ) then
          dtdp = (ttg(j,2)-ttg(j,1))/(pa(j,2)-pa(j,1))
        else if ( k == nk ) then
          dtdp = (ttg(j,nk)-ttg(j,nk-1))/(pa(j,nk)-pa(j,nk-1))
        else
          dtdp = (ttg(j,k+1)-ttg(j,k-1))/(pa(j,k+1)-pa(j,k-1))
        end if
        stab = max(rovcp*ttg(j,k)/pa(j,k) - dtdp,1.0e-4_rkx)
        ramp = min((pa(j,k)-1.0e4_rkx)/5.0e3_rkx,d_one)
        omls(j,k) = ramp*(ttg(j,k)-t(j,k))/(tauwtg*stab)
        kw = k
      end do
      if ( kw > 0 ) then
        do k = kw + 1, nk
          omls(j,k) = omls(j,kw)*za(j,k)/za(j,kw)
        end do
      end if
    end do
  end subroutine wtg_omega

  ! Ocean scheme of the model, on the lowest level and the SST
  subroutine ocean_fluxes
    implicit none
    real(rkx) :: shour
    do j = 1, ncol
      ht(j) = za(j,nk)
      usw(j) = sign(d_one,u(j,nk))*max(abs(u(j,nk)),0.001_rkx)
      vsw(j) = sign(d_one,v(j,nk))*max(abs(v(j,nk)),0.001_rkx)
      tatm(j) = t(j,nk)
      qv(j) = q(j,nk)
      patm(j) = pa(j,nk)
      sfps(j) = pf(j,nk+1)
      sfta(j) = tatm(j)*(sfps(j)/patm(j))**rovcp
      rhox(j) = m2p%rhox2d(j,1)
      tgb(j) = merge(col_sst(j),sfta(j),col_ocn(j))
      tgrd(j) = tgb(j)
      tgbrd(j) = tgb(j)
      hpbl(j) = p2m%zpbl(j,1)
      shour = modulo(col_hour(j) + time/3600.0_rkx,24.0_rkx)
      if ( ocnsw_now(j) >= d_zero ) then
        ! net radiation of the last call of the radiation
        rswf(j) = ocnsw_now(j)
        rlwf(j) = ocnlw_now(j)
      else
        rswf(j) = ocnsw*mathpi* &
          max(sin(mathpi*(shour-6.0_rkx)/12.0_rkx),d_zero)
        rlwf(j) = ocnlw
      end if
      um10(j) = sqrt(usw(j)**2+vsw(j)**2)
    end do
    select case ( iocnflux )
      case (1)
        call zengocndrv
      case (2)
        call coare3_drv
      case default
        call scm_error('iocnflux must be 0, 1 or 2')
    end select
    do j = 1, ncol
      if ( .not. col_ocn(j) ) cycle
      m2p%hfx(j,1) = sent(j)
      m2p%qfx(j,1) = evpr(j)
      m2p%uvdrag(j,1) = drag(j)
      m2p%tg(j,1) = tgb(j)
      ustar(j) = ustr(j)
    end do
  end subroutine ocean_fluxes

  subroutine init_ocean
    implicit none
    nocnp = ncol
    iocnbeg = 1
    iocnend = ncol
    iwavcpl = 0
    dtocn = dt
    ldcsst = .false.
    llake = .false.
    lseaice = .false.
    call getmem(mask,1,ncol,'scm:mask')
    call getmem(icpl,1,ncol,'scm:icpl')
    call getmem(ht,1,ncol,'scm:ht')
    call getmem(usw,1,ncol,'scm:usw')
    call getmem(vsw,1,ncol,'scm:vsw')
    call getmem(tatm,1,ncol,'scm:tatm')
    call getmem(qv,1,ncol,'scm:qv')
    call getmem(patm,1,ncol,'scm:patm')
    call getmem(sfps,1,ncol,'scm:sfps')
    call getmem(sfta,1,ncol,'scm:sfta')
    call getmem(tgb,1,ncol,'scm:tgb')
    call getmem(tgrd,1,ncol,'scm:tgrd')
    call getmem(tgbrd,1,ncol,'scm:tgbrd')
    call getmem(hpbl,1,ncol,'scm:hpbl')
    call getmem(rswf,1,ncol,'scm:rswf')
    call getmem(rlwf,1,ncol,'scm:rlwf')
    call getmem(rhox,1,ncol,'scm:rhox')
    call getmem(um10,1,ncol,'scm:um10')
    call getmem(sent,1,ncol,'scm:sent')
    call getmem(evpr,1,ncol,'scm:evpr')
    call getmem(drag,1,ncol,'scm:drag')
    call getmem(ustr,1,ncol,'scm:ustr')
    call getmem(zoo,1,ncol,'scm:zoo')
    call getmem(taux,1,ncol,'scm:taux')
    call getmem(tauy,1,ncol,'scm:tauy')
    call getmem(u10m,1,ncol,'scm:u10m')
    call getmem(v10m,1,ncol,'scm:v10m')
    call getmem(t2m,1,ncol,'scm:t2m')
    call getmem(q2m,1,ncol,'scm:q2m')
    call getmem(rhoa,1,ncol,'scm:rhoa')
    call getmem(ram1,1,ncol,'scm:ram1')
    call getmem(rah1,1,ncol,'scm:rah1')
    call getmem(br,1,ncol,'scm:br')
    call getmem(sfice,1,ncol,'scm:sfice')
    mask(:) = 1
    icpl(:) = 0
    sfice(:) = d_zero
  end subroutine init_ocean

  ! NoTo (ipptls = 2) with the SUBEX cloud fraction (icldfrac = 0), as
  ! the MOLOCH core calls them after the cumulus scheme
  subroutine init_micro
    implicit none
    dtsec = dt
    call allocate_mod_nogtom
    call getmem(mo2mc%xlat,1,ncol,1,1,'scm:mxlat')
    call getmem(mo2mc%psb,1,ncol,1,1,'scm:mpsb')
    call getmem(mo2mc%ldmsk,1,ncol,1,1,'scm:mldmsk')
    call getmem(mo2mc%phs,1,ncol,1,1,1,nk,'scm:mphs')
    call getmem(mo2mc%pfs,1,ncol,1,1,1,nk+1,'scm:mpfs')
    call getmem(mo2mc%delz,1,ncol,1,1,1,nk,'scm:mdelz')
    call getmem(mo2mc%t,1,ncol,1,1,1,nk,'scm:mt')
    call getmem(mo2mc%rho,1,ncol,1,1,1,nk,'scm:mrho')
    call getmem(mo2mc%pverv,1,ncol,1,1,1,nk,'scm:mpverv')
    call getmem(mo2mc%heatrt,1,ncol,1,1,1,nk,'scm:mheatrt')
    call getmem(mo2mc%qdetr,1,ncol,1,1,1,nk,'scm:mqdetr')
    call getmem(mo2mc%qs,1,ncol,1,1,1,nk,'scm:mqs')
    call getmem(mo2mc%rh,1,ncol,1,1,1,nk,'scm:mrh')
    call getmem(mo2mc%cldf,1,ncol,1,1,1,nk,'scm:mcldf')
    call getmem(mo2mc%qxx,1,ncol,1,1,1,nk,1,nqx,'scm:mqxx')
    call getmem(qvp,1,ncol,1,1,1,nk,'scm:qvp')
    call getmem(rhcp,1,ncol,1,1,1,nk,'scm:rhcp')
    call getmem(totcp,1,ncol,1,1,1,nk,'scm:totcp')
    call getmem(rh0p,1,ncol,1,1,'scm:rh0p')
    call getmem(mc2mo%tten,1,ncol,1,1,1,nk,'scm:ctten')
    call getmem(mc2mo%qxten,1,ncol,1,1,1,nk,1,nqx,'scm:cqxten')
    call getmem(mc2mo%rainnc,1,ncol,1,1,'scm:crainnc')
    call getmem(mc2mo%lsmrnc,1,ncol,1,1,'scm:clsmrnc')
    call getmem(mc2mo%snownc,1,ncol,1,1,'scm:csnownc')
    call getmem(mc2mo%trrate,1,ncol,1,1,'scm:ctrrate')
    call getmem(mc2mo%rainls,1,ncol,1,1,1,nk,'scm:crainls')
    call getmem(mc2mo%remrat,1,ncol,1,1,1,nk,'scm:cremrat')
    call getmem(mc2mo%rembc,1,ncol,1,1,1,nk,'scm:crembc')
    do j = 1, ncol
      mo2mc%xlat(j,1) = col_lat(j)
      mo2mc%psb(j,1) = col_ps(j)
      mo2mc%ldmsk(j,1) = col_land(j)
      rh0p(j,1) = merge(rh0land,rh0oce,col_land(j) == 1)
      do k = 1, nk
        mo2mc%phs(j,1,k) = pa(j,k)
        mo2mc%delz(j,1,k) = zf(j,k) - zf(j,k+1)
      end do
      do k = 1, nk + 1
        mo2mc%pfs(j,1,k) = pf(j,k)
      end do
    end do
    call init_nogtom(mo2mc%ldmsk)
  end subroutine init_micro

  subroutine micro_step
    implicit none
    do j = 1, ncol
      do k = 1, nk
        mo2mc%t(j,1,k) = t(j,k)
        mo2mc%rho(j,1,k) = pa(j,k)/(rgas*t(j,k))
        mo2mc%pverv(j,1,k) = omls(j,k)
        mo2mc%heatrt(j,1,k) = hrad(j,k)
        mo2mc%qdetr(j,1,k) = cu_qdetr(j,1,k)*dt
        mo2mc%qxx(j,1,k,iqv) = q(j,k)
        do n = iqc, iqs
          mo2mc%qxx(j,1,k,n) = qcnd(j,k,n)
        end do
        mo2mc%qs(j,1,k) = pfwsat(t(j,k),pa(j,k))
        mo2mc%rh(j,1,k) = min(max(q(j,k)/mo2mc%qs(j,1,k),d_zero),d_one)
        qvp(j,1,k) = q(j,k)
        ! cloud water + ice, as cldfrac in mod_micro_interface
        totcp(j,1,k) = &
          merge(qcnd(j,k,iqc),d_zero,qcnd(j,k,iqc) >= 2.0e-7_rkx) + &
          merge(qcnd(j,k,iqi),d_zero,qcnd(j,k,iqi) >= 2.0e-7_rkx)
        ! Tiedtke: the detrained condensate goes to the microphysics
        ! through qdetr; Kain-Fritsch: as tendencies
        if ( icumulus == 5 ) then
          cu_qten(j,1,k,iqc) = d_zero
          cu_qten(j,1,k,iqi) = d_zero
        end if
        mc2mo%tten(j,1,k) = cu_tten(j,1,k)
        do n = 1, nqx
          mc2mo%qxten(j,1,k,n) = cu_qten(j,1,k,n)
        end do
      end do
      mc2mo%rainnc(j,1) = d_zero
      mc2mo%lsmrnc(j,1) = d_zero
      mc2mo%snownc(j,1) = d_zero
      mc2mo%trrate(j,1) = d_zero
    end do
    select case ( icldfrac )
      case (0)
        call subex_cldfrac(mo2mc%t,mo2mc%phs,qvp,mo2mc%rh,tc0,rh0p,rhcp, &
          mo2mc%cldf)
      case (1)
        call xuran_cldfrac(totcp,mo2mc%qs,mo2mc%rh,mo2mc%cldf)
      case default
        call scm_error('icldfrac must be 0 (SUBEX) or 1 (Xu-Randall)')
    end select
    call nogtom(mo2mc,mc2mo)
    do j = 1, ncol
      do k = 1, nk
        cfr(j,k) = mo2mc%cldf(j,1,k)
        mtt(j,k) = mc2mo%tten(j,1,k) - cu_tten(j,1,k)
        mqt(j,k) = mc2mo%qxten(j,1,k,iqv) - cu_qten(j,1,k,iqv)
        ! the condensate tendencies of the cumulus are in those of the
        ! microphysics (mc2mo%qxten), which advance qcnd
      end do
    end do
  end subroutine micro_step

  ! A 1 x n grid: every index range of the model domain is 1:n x 1:1
  subroutine set_grid(n)
    implicit none
    integer(ik4), intent(in) :: n
    jci1 = 1 ; jci2 = n ; ici1 = 1 ; ici2 = 1
    jce1 = 1 ; jce2 = n ; ice1 = 1 ; ice2 = 1
    jdi1 = 1 ; jdi2 = n ; idi1 = 1 ; idi2 = 1
    jci1ga = 1 ; jci2ga = n ; ici1ga = 1 ; ici2ga = 1
    jce1ga = 1 ; jce2ga = n ; ice1ga = 1 ; ice2ga = 1
    jce1gb = 1 ; jce2gb = n ; ice1gb = 1 ; ice2gb = 1
  end subroutine set_grid

  ! Boundary layer input derived from the state: potential temperature,
  ! surface air density, ground temperature (same potential temperature
  ! as the lowest level, used only with ifaholtth10 /= 2) and the surface
  ! drag from a neutral log law friction velocity (or ustarfix).
  subroutine update_derived
    implicit none
    real(rkx) :: ps, tv, wspd
    do j = 1, ncol
      m2p%tatm(j,1,:) = t(j,:)
      m2p%heatrt(j,1,:) = hrad(j,:)
      if ( lnoto ) m2p%fcc(j,1,:) = mo2mc%cldf(j,1,:)
      m2p%qxatm(j,1,:,iqv) = q(j,:)
      m2p%qxatm(j,1,:,iqc) = qcnd(j,:,iqc)
      m2p%qxatm(j,1,:,iqi) = qcnd(j,:,iqi)
      m2p%uxatm(j,1,:) = u(j,:)
      m2p%vxatm(j,1,:) = v(j,:)
      do k = 1, nk
        m2p%thatm(j,1,k) = t(j,k)*(p00/pa(j,k))**rovcp
      end do
      ps = pf(j,nk+1)
      tv = t(j,nk)*(d_one+ep1*q(j,nk)/(d_one+q(j,nk)))
      m2p%rhox2d(j,1) = ps/(rgas*tv)
      m2p%tg(j,1) = m2p%thatm(j,1,nk)*(ps/p00)**rovcp
      wspd = max(sqrt(u(j,nk)**2+v(j,nk)**2),0.1_rkx)
      zl = max(za(j,nk),d_two*col_z0(j))
      ws = vonkar*wspd/log(zl/col_z0(j))
      if ( ustarfix > d_zero ) ws = ustarfix
      ustar(j) = ws
      m2p%uvdrag(j,1) = m2p%rhox2d(j,1)*ws*ws/wspd
    end do
  end subroutine update_derived

  subroutine allocate_fields
    implicit none
    call getmem(m2c%pas,1,ncol,1,1,1,nk,'scm:pas')
    call getmem(m2c%pasf,1,ncol,1,1,1,nk+1,'scm:pasf')
    call getmem(m2c%zas,1,ncol,1,1,1,nk,'scm:zas')
    call getmem(m2c%zfs,1,ncol,1,1,1,nk+1,'scm:zfs')
    call getmem(m2c%dzq,1,ncol,1,1,1,nk,'scm:dzq')
    call getmem(m2c%rhoas,1,ncol,1,1,1,nk,'scm:rhoas')
    call getmem(m2c%tas,1,ncol,1,1,1,nk,'scm:tas')
    call getmem(m2c%uas,1,ncol,1,1,1,nk,'scm:uas')
    call getmem(m2c%vas,1,ncol,1,1,1,nk,'scm:vas')
    call getmem(m2c%qxas,1,ncol,1,1,1,nk,1,nqx,'scm:qxas')
    call getmem(m2c%tten,1,ncol,1,1,1,nk,'scm:ctten')
    call getmem(m2c%uten,1,ncol,1,1,1,nk,'scm:cuten')
    call getmem(m2c%vten,1,ncol,1,1,1,nk,'scm:cvten')
    call getmem(m2c%qxten,1,ncol,1,1,1,nk,1,nqx,'scm:cqxten')
    call getmem(m2c%qfx,1,ncol,1,1,'scm:cqfx')
    call getmem(m2c%hfx,1,ncol,1,1,'scm:chfx')
    call getmem(m2c%zpbl,1,ncol,1,1,'scm:czpbl')
    call getmem(m2c%ustar,1,ncol,1,1,'scm:ustar')
    call getmem(m2c%psb,1,ncol,1,1,'scm:psb')
    call getmem(m2c%ldmsk,1,ncol,1,1,'scm:cldmsk')
    m2c%qxas(:,:,:,:) = d_zero
    m2c%uten(:,:,:) = d_zero
    m2c%vten(:,:,:) = d_zero
    call getmem(m2p%za,1,ncol,1,1,1,nk,'scm:za')
    call getmem(m2p%heatrt,1,ncol,1,1,1,nk,'scm:pheatrt')
    m2p%heatrt(:,:,:) = d_zero
    ! large scale cloud fraction of the model (cldfrac)
    call getmem(m2p%fcc,1,ncol,1,1,1,nk,'scm:pfcc')
    m2p%fcc(:,:,:) = d_zero
    call getmem(m2p%patm,1,ncol,1,1,1,nk,'scm:patm')
    call getmem(m2p%patmf,1,ncol,1,1,1,nk+1,'scm:patmf')
    call getmem(m2p%tatm,1,ncol,1,1,1,nk,'scm:tatm')
    call getmem(m2p%thatm,1,ncol,1,1,1,nk,'scm:thatm')
    call getmem(m2p%uxatm,1,ncol,1,1,1,nk,'scm:uxatm')
    call getmem(m2p%vxatm,1,ncol,1,1,1,nk,'scm:vxatm')
    call getmem(m2p%qxatm,1,ncol,1,1,1,nk,1,nqx,'scm:qxatm')
    call getmem(m2p%hfx,1,ncol,1,1,'scm:hfx')
    call getmem(m2p%qfx,1,ncol,1,1,'scm:qfx')
    call getmem(m2p%uvdrag,1,ncol,1,1,'scm:uvdrag')
    call getmem(m2p%rhox2d,1,ncol,1,1,'scm:rhox2d')
    call getmem(m2p%tg,1,ncol,1,1,'scm:tg')
    call getmem(m2p%coriol,1,ncol,1,1,'scm:coriol')
    call getmem(m2p%ldmsk,1,ncol,1,1,'scm:ldmsk')
    call getmem(p2m%tten,1,ncol,1,1,1,nk,'scm:tten')
    call getmem(p2m%uten,1,ncol,1,1,1,nk,'scm:uten')
    call getmem(p2m%vten,1,ncol,1,1,1,nk,'scm:vten')
    call getmem(p2m%qxten,1,ncol,1,1,1,nk,1,nqx,'scm:qxten')
    call getmem(p2m%zpbl,1,ncol,1,1,'scm:zpbl')
    call getmem(p2m%kpbl,1,ncol,1,1,'scm:kpbl')
    m2p%qxatm(:,:,:,:) = d_zero
  end subroutine allocate_fields

  ! Radiation of the model (irrtm = 0, the CCM scheme of colmod3) on the
  ! columns: allocations as allocate_radiation of mod_rad_interface, no
  ! aerosol, the built in ozone profile (iclimao3 = 0), the GHG of the
  ! scenario table for the year of the case, the orbit of that year.
  subroutine init_rad
    implicit none
    integer(ik4) :: iday, ihour
    integer(ik4), dimension(12), parameter :: mlen = &
      [31,28,31,30,31,30,31,31,30,31,30,31]
    irrtm = 0
    kth = kz
    ktf = kzp1
    iclimao3 = 0
    iclimaaer = 0
    idirect = 0
    itweak = 0
    lsrfhack = .false.
    ifrad = .true.
    ifopt = .false.
    ! Gregorian calendar constants of the orbit (set by the model from
    ! the namelist calendar)
    dayspy = 365.2422_rkx
    vernal_equinox = real((80.50_rk8/dayspy)*365.0_rk8,rkx)
    dpd = 360.0_rkx/dayspy
    half_dayspy = dayspy/2.0_rkx
    sixteenth_dayspy = dayspy/16.0_rkx
    solcon = 1.36822e+03_rk8
    scon = solcon*1000.0_rk8
    nrad = max(nint(dtrad/dt),1)
    nabs = max(nint(dtabem/dt),1)
    read(col_date(1),'(i4,3i2)') iyear, imonth, iday, ihour
    cday0 = real(sum(mlen(1:imonth-1)) + iday,rk8)
    if ( imonth > 2 .and. mod(iyear,4) == 0 .and. &
         (mod(iyear,100) /= 0 .or. mod(iyear,400) == 0) ) then
      cday0 = cday0 + 1.0_rk8
    end if
    cday0 = cday0 + real(ihour,rk8)/24.0_rk8
    call orb_params(iyear,eccen,obliq,mvelp,obliqr,lambm0,mvelpp)
    call set_scenario(scenario,iyear,imonth)
    call getmem(o3prof,jci1,jci2,ici1,ici2,1,kzp1,'scm:o3prof')
    call allocate_mod_rad_aerosol
    call allocate_mod_rad_o3blk
    call allocate_mod_rad_outrad
    call allocate_mod_rad_colmod3
    call getmem(gasabsnxt,jci1,jci2,ici1,ici2,1,kz,1,4,'scm:gasabsnxt')
    call getmem(gasabstot,jci1,jci2,ici1,ici2,1,kzp1,1,kzp1,'scm:gasabstot')
    call getmem(gasemstot,jci1,jci2,ici1,ici2,1,kzp1,'scm:gasemstot')
    call getmem(absemin,jci1,jci2,ici1,ici2,1,5*kz+8,'scm:absemin')
    call getmem(m2r%tatms,1,ncol,1,1,1,nk,'scm:rt')
    call getmem(m2r%qxatms,1,ncol,1,1,1,nk,1,nqx,'scm:rqx')
    call getmem(m2r%rhatms,1,ncol,1,1,1,nk,'scm:rrh')
    call getmem(m2r%rhoatms,1,ncol,1,1,1,nk,'scm:rrho')
    call getmem(m2r%phatms,1,ncol,1,1,1,nk,'scm:rph')
    call getmem(m2r%pfatms,1,ncol,1,1,1,nk+1,'scm:rpf')
    call getmem(m2r%deltaz,1,ncol,1,1,1,nk,'scm:rdz')
    call getmem(m2r%cldfrc,1,ncol,1,1,1,nk,'scm:rcld')
    call getmem(m2r%cldlwc,1,ncol,1,1,1,nk,'scm:rlwc')
    call getmem(m2r%psatms,1,ncol,1,1,'scm:rps')
    call getmem(m2r%tg,1,ncol,1,1,'scm:rtg')
    call getmem(m2r%xlat,1,ncol,1,1,'scm:rlat')
    call getmem(m2r%xlon,1,ncol,1,1,'scm:rlon')
    call getmem(m2r%ptrop,1,ncol,1,1,'scm:rptrop')
    call getmem(m2r%coszrs,1,ncol,1,1,'scm:rcosz')
    call getmem(m2r%albvs,1,ncol,1,1,'scm:ralbvs')
    call getmem(m2r%albvl,1,ncol,1,1,'scm:ralbvl')
    call getmem(m2r%aldirs,1,ncol,1,1,'scm:raldirs')
    call getmem(m2r%aldifs,1,ncol,1,1,'scm:raldifs')
    call getmem(m2r%aldirl,1,ncol,1,1,'scm:raldirl')
    call getmem(m2r%aldifl,1,ncol,1,1,'scm:raldifl')
    call getmem(m2r%emiss,1,ncol,1,1,'scm:remiss')
    call getmem(m2r%ldmsk,1,ncol,1,1,'scm:rldmsk')
    call getmem(r2m%heatrt,1,ncol,1,1,1,nk,'scm:rheat')
    call getmem(r2m%totcf,1,ncol,1,1,'scm:rtotcf')
    call getmem(r2m%solis,1,ncol,1,1,'scm:rsolis')
    call getmem(r2m%fsw,1,ncol,1,1,'scm:rfsw')
    call getmem(r2m%flw,1,ncol,1,1,'scm:rflw')
    call getmem(r2m%flwd,1,ncol,1,1,'scm:rflwd')
    call getmem(r2m%sabveg,1,ncol,1,1,'scm:rsabveg')
    call getmem(r2m%solvs,1,ncol,1,1,'scm:rsolvs')
    call getmem(r2m%solvsd,1,ncol,1,1,'scm:rsolvsd')
    call getmem(r2m%solvl,1,ncol,1,1,'scm:rsolvl')
    call getmem(r2m%solvld,1,ncol,1,1,'scm:rsolvld')
    call getmem(r2m%sinc,1,ncol,1,1,'scm:rsinc')
    do j = 1, ncol
      m2r%xlat(j,1) = col_lat(j)
      m2r%xlon(j,1) = col_lon(j)
      m2r%ldmsk(j,1) = col_land(j)
      m2r%ptrop(j,1) = 10000.0_rkx
      m2r%psatms(j,1) = pf(j,nk+1)
      do k = 1, nk
        m2r%phatms(j,1,k) = pa(j,k)
        m2r%deltaz(j,1,k) = zf(j,k) - zf(j,k+1)
      end do
      do k = 1, nk + 1
        m2r%pfatms(j,1,k) = pf(j,k)
      end do
    end do
    m2r%qxatms(:,:,:,:) = d_zero
    call o3data(m2r)
  end subroutine init_rad

  ! Cloud fraction and in cloud water for the radiation, as cucloud and
  ! cldfrac of the model (icup > 1, iconvlwp = 0): the convective cloud
  ! of the cumulus scheme (at most clfrcv) overlapped with the large
  ! scale cloud of the microphysics, with the Liang (2005) scaling.
  subroutine rad_clouds
    implicit none
    real(rkx) :: cuf, cwc, fls, totc, w2, cfq
    do j = 1, ncol
      do k = 1, nk
        cuf = min(cu_cldfrc(j,1,k),clfrcv)
        cwc = cu_cldlwc(j,1,k)
        if ( cwc <= d_zero ) cwc = clwfromt(t(j,k))
        fls = d_zero
        totc = d_zero
        if ( lnoto ) then
          fls = min(max(mo2mc%cldf(j,1,k),d_zero),d_one)
          totc = merge(qcnd(j,k,iqc),d_zero,qcnd(j,k,iqc) >= 2.0e-7_rkx) + &
                 merge(qcnd(j,k,iqi),d_zero,qcnd(j,k,iqi) >= 2.0e-7_rkx)
          totc = totc*1000.0_rkx*pa(j,k)/(rgas*t(j,k))
        end if
        if ( cuf > lowcld ) then
          if ( fls > lowcld ) then
            cfrad(j,k) = max(fls,cuf)
            w2 = fls/cfrad(j,k)
            lwcrad(j,k) = w2*totc/fls + (d_one-w2)*cwc
          else
            cfrad(j,k) = cuf
            lwcrad(j,k) = cwc
          end if
        else if ( fls > lowcld ) then
          cfrad(j,k) = fls
          lwcrad(j,k) = totc/fls
        else
          cfrad(j,k) = d_zero
          lwcrad(j,k) = d_zero
        end if
        if ( lwcrad(j,k) > d_zero ) then
          cfq = real(int(cfrad(j,k)*255.0_rkx),rkx)/255.0_rkx
          lwcrad(j,k) = lwcrad(j,k) * &
            (0.97_rkx*exp(-((cfq-0.098_rkx)**2)/0.0365_rkx)+0.255_rkx)
        end if
      end do
    end do
  end subroutine rad_clouds

  ! One call of the radiation on the current state: sun position for
  ! the time of the step, surface albedo (ocean: the model formula with
  ! iwhitecap = 0; land: albland) and emissivity, ground temperature of
  ! the surface (SST over the ocean).
  subroutine rad_step
    implicit none
    real(rk8) :: cday
    real(rkx) :: cz, alb
    logical :: labsem
    cday = cday0 + real(time,rk8)/86400.0_rk8
    call orb_decl(real(int(cday),rk8),eccen,mvelpp,lambm0,obliqr, &
      declin,eccf)
    do j = 1, ncol
      cz = real(orb_cosz(cday,real(col_lat(j),rk8),real(col_lon(j),rk8), &
        declin),rkx)
      cz = min(max(cz,d_zero),d_one)
      m2r%coszrs(j,1) = cz
      if ( col_land(j) == 1 ) then
        m2r%aldirs(j,1) = albland
        m2r%aldifs(j,1) = albland
        m2r%emiss(j,1) = 0.96_rkx
      else
        alb = 0.05_rkx/(cz+0.15_rkx)
        m2r%aldirs(j,1) = alb
        m2r%aldifs(j,1) = 0.08_rkx
        m2r%emiss(j,1) = 0.98_rkx
      end if
      m2r%aldirl(j,1) = m2r%aldirs(j,1)
      m2r%aldifl(j,1) = m2r%aldifs(j,1)
      m2r%albvs(j,1) = m2r%aldirs(j,1)
      m2r%albvl(j,1) = m2r%aldirl(j,1)
      m2r%tg(j,1) = m2p%tg(j,1)
      do k = 1, nk
        m2r%tatms(j,1,k) = t(j,k)
        m2r%qxatms(j,1,k,iqv) = q(j,k)
        m2r%qxatms(j,1,k,iqc) = qcnd(j,k,iqc)
        m2r%qxatms(j,1,k,iqi) = qcnd(j,k,iqi)
        m2r%rhatms(j,1,k) = min(max(q(j,k)/pfwsat(t(j,k),pa(j,k)), &
          d_zero),d_one)
        m2r%rhoatms(j,1,k) = pa(j,k)/(rgas*t(j,k))
        m2r%cldfrc(j,1,k) = cfrad(j,k)
        m2r%cldlwc(j,1,k) = lwcrad(j,k)
      end do
    end do
    labsem = mod(istep-1,nabs) == 0
    call colmod3(iyear,imonth,.false.,labsem,m2r,r2m)
    do j = 1, ncol
      do k = 1, nk
        hrad(j,k) = r2m%heatrt(j,1,k)
      end do
      if ( col_ocn(j) ) then
        ocnsw_now(j) = r2m%fsw(j,1)
        ocnlw_now(j) = r2m%flw(j,1)
      end if
      ! TOA fluxes: incoming minus net (clear sky) short wave, long wave
      rsut(j) = real(curflx(j,8),rkx)
      rlut(j) = real(curflx(j,9),rkx)
      rsutcs(j) = real(curflx(j,7)-curflx(j,3),rkx)
      rlutcs(j) = real(curflx(j,5),rkx)
    end do
  end subroutine rad_step

#include <clwfromt.inc>

#include <saturation.inc>

end program scm_coupled

! The library calls myabort (Main/abort.F90) on fatal errors.
subroutine myabort
  implicit none
  stop 2
end subroutine myabort

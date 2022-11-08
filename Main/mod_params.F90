!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_params

  use mod_runparams
  use mod_mppparam
  use mod_mpmessage
  use mod_domain
  use mod_service
  use mod_cu_interface
  use mod_lm_interface
  use mod_atm_interface
  use mod_che_interface
  use mod_rad_interface
  use mod_pbl_interface
#ifdef CLM45
  use mod_clm_regcm
#endif
  use mod_micro_interface
  use mod_split
  use mod_slice
  use mod_bdycod
  use mod_ncio
  use mod_tendency
  use mod_ncout
  use mod_advection , only : init_advection
  use mod_sladvection , only : init_sladvection
  use mod_diffusion , only : allocate_mod_diffusion
  use mod_savefile
  use mod_slabocean
  use mod_sldepparam
  use mod_sound
  use mod_timer
  use mod_zita
  use mod_moloch
  use mod_timefilter

  implicit none

  private

  real(rkx) , parameter :: mindt = 1.0_rkx

  public :: param

  contains
  !
  ! This subroutine defines the various model parameters.
  !
  subroutine param
    implicit none
    real(rkx) :: afracl , afracs , bb , cc , dlargc , dsmalc , dxtemc , &
               qk , qkp1 , sig700 , ssum , vqmax , wk , wkp1 , xbot ,   &
               xtop , xx , yy , mo_c1 , mo_c2 , dl , minfrq
    real(rkx) , dimension(kzp1) :: fak , fbk
    integer(ik4) :: kbmax
    integer(ik4) :: iretval
    integer(ik4) :: i , j , k , kbase , ktop , ns
    integer(ik8) :: mdate0 , mdate1 , mdate2
    integer(ik4) :: hspan , ipunit
    integer(ik4) :: n , len_path
    character(len=32) :: appdat
    type(rcm_time_interval) :: bdif
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'param'
    integer(ik4) , save :: idindx = 0
#endif
    !
    ! namelist:
    !
    namelist /restartparam/ ifrest , mdate0 , mdate1 , mdate2

    namelist /timeparam/ dtrad , dtsrf , dtcum , dtche , dtabem , dt

    namelist /outparam/ prestr , ifsave , ifatm , ifrad , ifsrf , ifsub , &
      iflak , ifshf , ifsts , ifchem , ifopt , outnwf , savfrq , atmfrq , &
      srffrq , subfrq , lakfrq , radfrq , chemfrq ,optfrq, dirout ,       &
      uvrotate , enable_atm_vars , enable_srf_vars , enable_rad_vars ,    &
      enable_sub_vars , enable_sts_vars , enable_lak_vars ,               &
      enable_opt_vars , enable_che_vars , enable_shf_vars ,               &
      lsync , idiag , icosp , deflate_level , do_parallel_netcdf_in ,     &
      do_parallel_netcdf_out , deflate_level , chechgact

    namelist /physicsparam/ ibltyp , iboudy , isladvec , iqmsl ,         &
      icup_lnd , icup_ocn , ipgf , iemiss , lakemod , ipptls , idiffu ,  &
      iocnflx , iocncpl , iwavcpl , icopcpl , iocnrough , iocnzoq ,      &
      ichem ,  scenario ,  idcsst , ipcpcool , iwhitecap , iseaice ,     &
      idesseas , iconvlwp , icldmstrat , icldfrac , irrtm , iclimao3 ,   &
      iclimaaer , isolconst , icumcloud , islab_ocean , itweak ,         &
      temp_tend_maxval , wind_tend_maxval , ghg_year_const , ifixsolar , &
      fixedsolarval , irceideal , year_offset , radclimpath

    namelist /dynparam/ gnu1 , gnu2 , diffu_hgtf , ckh , adyndif , &
      upstream_mode , uoffc , stability_enhance , t_extrema ,      &
      q_rel_extrema

    namelist /hydroparam/ nsplit , lstand

    namelist /nonhydroparam/ ifupr , nhbet , nhxkd ,       &
      ifrayd , rayndamp , rayalpha0 , rayhd , itopnudge ,  &
      mo_anu2 , mo_nadv , mo_nsound , mo_wmax , mo_nzfilt

    namelist /rrtmparam/ inflgsw , iceflgsw , liqflgsw , inflglw ,    &
      iceflglw , liqflglw , icld , irng , imcica , nradfo

    namelist /cldparam/ ncld , rhmax , rhmin , rh0oce , rh0land , tc0 ,  &
      cllwcv , clfrcvmax , cftotmax , kfac_shal , kfac_deep , k2_const , &
      lsrfhack , larcticcorr , rcrit , coef_ccn , abulk

    namelist /subexparam/ qck1land , qck1oce , gulland , guloce ,  &
      cevaplnd , cevapoce , caccrlnd , caccroce , conf

    namelist /microparam/ stats , budget_compute , nssopt ,  &
      iautoconv , vfqr , vfqi , vfqs , auto_rate_khair ,     &
      auto_rate_kessl , auto_rate_klepi , rkconv , skconv ,  &
      rcovpmin , rpecons , rcldiff

    namelist /grellparam/ igcc , shrmin , shrmax , edtmin , &
      edtmax , edtmino , edtmaxo , edtminx , edtmaxx , pbcmax ,    &
      mincld , htmin , htmax , skbmax , dtauc, shrmin_ocn ,        &
      shrmax_ocn , edtmin_ocn, edtmax_ocn, edtmino_ocn ,           &
      edtmaxo_ocn , edtminx_ocn , edtmaxx_ocn

    namelist /emanparam/ minorig , elcrit_ocn , elcrit_lnd , tlcrit , &
      entp , sigd , sigs , omtrain , omtsnow , coeffr , coeffs , cu , &
      betae , dtmax , alphae , damp , epmax_ocn , epmax_lnd ,         &
      istochastic

    namelist /emanstochastic/ epmax_lnd_min , epmax_lnd_max , &
      elcrit_lnd_min , elcrit_lnd_max ,                       &
      sigs_min , sigs_max , sigd_min , sigd_max

    namelist /tiedtkeparam/ iconv , entrmax , entrdd , entrpen_lnd , &
      entrpen_ocn , entrscv , entrmid , cprcon , detrpen_lnd ,       &
      detrpen_ocn , entshalp , rcuc_lnd , rcuc_ocn , rcpec_lnd ,     &
      rcpec_ocn , rhebc_lnd , rhebc_ocn , rprc_ocn , rprc_lnd ,      &
      revap_lnd , revap_ocn , cmtcape , lmfpen , lmfmid , lmfdd ,    &
      lepcld , lmfdudv , lmfscv , lmfuvdis , lmftrac , lmfsmooth ,   &
      lmfwstar

    namelist /kfparam/ kf_min_pef , kf_max_pef , kf_entrate , kf_dpp , &
      kf_min_dtcape , kf_max_dtcape , kf_tkemax , kf_convrate ,        &
      kf_wthreshold

    namelist /chemparam/ chemsimtype , ichremlsc , ichremcvc , ichdrdepo , &
      ichcumtra , ichsolver , idirect , iindirect , ichdustemd ,           &
      ichdiag , ichsursrc , ichebdy , rdstemfac , ichjphcld , ichbion ,    &
      ismoke , rocemfac, ichlinox , isnowdark, ichdustparam , ichecold ,   &
      carb_aging_control

    namelist /uwparam/ iuwvadv , atwo , rstbl , czero , nuk

    namelist /holtslagparam/ ricr_ocn , ricr_lnd , zhnew_fac , &
      ifaholtth10 , ifaholt , holtth10iter

#ifdef CLM
    namelist /clmparam/ dirclm , imask , clmfrq , ilawrence_albedo
#endif

    namelist /cplparam/ cpldt , zomax , ustarmax

    namelist /slabocparam/ do_qflux_adj , do_restore_sst , &
      sst_restore_timescale , mixed_layer_depth

    namelist /tweakparam/ itweak_temperature , itweak_solar_irradiance , &
            itweak_sst , itweak_greenhouse_gases , temperature_tweak ,   &
            sst_tweak , solar_tweak , gas_tweak_factors

#ifdef DEBUG
    call time_begin(subroutine_name,idindx)
#endif
    !
    ! default values for all the options:
    !     (can be overwritten by namelist input).
    !
    ! restartparam ;
    !
    ifrest = .false.     ! Restart?:  t=true; f=false
    idate0 = 1900010100  ! Start date of the initial simulation
    idate1 = 1900010100  ! Start date of this simulation
    idate2 = 1900010100  ! End Date this simulation
    !
    ! timeparam ;
    !
    dt = 100.0_rkx  ! time step in seconds
    dtrad = 0.0_rkx ! time interval in min solar rad caluclated
    dtsrf = 0.0_rkx ! time interval at which bats is called (secs)
    dtcum = 0.0_rkx ! time interval at which cumulus is called (secs)
    dtabem = 0.0_rkx ! time interval absorption-emission calculated (hours)
    dtche = 900.0_rkx ! time interval at which bats is called (secs)
    !
    ! outparam ;
    !
    prestr = ''
    ifsave = .true.
    ifatm  = .true.
    ifrad  = .true.
    ifsrf  = .true.
    ifsts  = .true.
    ifshf  = .false.
    ifsub  = .false.
    iflak  = .false.
    ifopt  = .false.
    ifchem = .false.
    outnwf  = 0.0_rkx ! Frequency in days to open new files.
    savfrq  = 0.0_rkx ! time interval for disposing sav output (days)
    atmfrq  = 6.0_rkx ! time interval for disposing atm output (hrs)
    radfrq  = 6.0_rkx ! time interval for disposing rad output (hrs)
    srffrq  = 3.0_rkx ! time interval for disposing srf output (hrs)
    lakfrq  = 6.0_rkx ! time interval for disposing lake output (hrs)
    subfrq  = 6.0_rkx ! time interval for disposing lake output (hrs)
    chemfrq = 6.0_rkx ! time interval for disposing chem output (hrs)
    optfrq =  6.0_rkx ! time interval for disposing opt output (hrs)
    enable_atm_vars(:) = .true.
    enable_srf_vars(:) = .true.
    enable_sts_vars(:) = .true.
    enable_sub_vars(:) = .true.
    enable_lak_vars(:) = .true.
    enable_rad_vars(:) = .true.
    enable_opt_vars(:) = .true.
    enable_che_vars(:) = .true.
    enable_shf_vars(:) = .true.
    dirout = './output'
    lsync = .true.
    uvrotate = .false.
    do_parallel_netcdf_in = .false.
    do_parallel_netcdf_out = .false.
    chechgact = .false.
    idiag = 0
    icosp = 0
    !
    ! physicsparam ;
    !
    iboudy = 5
    ibltyp = 1
    isladvec = 0
    iqmsl = 1
    icup_lnd = 4
    icup_ocn = 4
    ipptls = 1
    idiffu = 1
    ipgf = 0
    iemiss = 0
    iocnflx = 2
    iocnrough = 1
    iocnzoq = 1
    iocncpl = 0
    iwavcpl = 0
    icopcpl = 0
    lakemod = 0
    ichem = 0
    scenario = 'RCP4.5'
    ghg_year_const = 1950
    idcsst = 0
    ipcpcool = 0
    iwhitecap = 0
    iseaice = 0
    idesseas = 0
    iconvlwp = 1
    icldfrac = 0
    icldmstrat = 0
    irrtm = 0
    islab_ocean = 0
    iclimao3 = 0
    iclimaaer = 0
    radclimpath = 'OPPMONTH'
    isolconst = 0
    year_offset = 0
    ifixsolar = 0
    fixedsolarval = 343.0_rkx
    irceideal = 0
    icumcloud = 1
    temp_tend_maxval = 5.0_rkx*(dt/secpm)
    wind_tend_maxval = 5.0_rkx*(dt/secpm)
    !
    ! Hydrostatic param ;
    !
    nsplit = 2
    lstand = .true.
    !
    ! Non hydrostatic param ;
    !
    ifupr = 1
    nhbet = 0.4_rkx   ! Arakawa beta (MM5 manual, Sec. 2.5.1)
    nhxkd = 0.1_rkx
    itopnudge = 0
    ifrayd = 0
    rayndamp = 5
    rayalpha0 = 1.0_rkx/86400.0_rkx
    rayhd = 10000.0_rkx
    mo_wmax = 150.0_rkx
    mo_nadv = 3
    mo_nsound = 5
    mo_anu2 = 0.05_rkx
    mo_nzfilt = 0
    !
    ! Rrtm radiation param ;
    !
    inflgsw  = 2
    iceflgsw = 3
    liqflgsw = 1
    inflglw  = 2
    iceflglw = 3
    liqflglw = 1
    icld  = 1
    imcica = 1
    irng = 1
    nradfo = 4
    !
    ! Subexparam ;
    ! From Pal et al, 2000
    !
    qck1land  = 0.0005_rkx ! Autoconversion Rate for Land
    qck1oce   = 0.0005_rkx ! Autoconversion Rate for Ocean
    gulland   = 0.65_rkx   ! Fract of Gultepe eqn (qcth) when prcp occurs (land)
    guloce    = 0.30_rkx   ! Fract of Gultepe eqn (qcth) for ocean
    cevaplnd  = 1.0e-5_rkx ! Raindrop ev rate coef land [[(kg m-2 s-1)-1/2]/s]
    cevapoce  = 1.0e-5_rkx ! Raindrop ev rate coef ocean [[(kg m-2 s-1)-1/2]/s]
    caccrlnd  = 6.0_rkx    ! Raindrop accretion rate land  [m3/kg/s]
    caccroce  = 4.0_rkx    ! Raindrop accretion rate ocean [m3/kg/s]
    conf      = 1.00_rkx   ! Condensation efficiency
    !
    ! Cloud fraction control algorithm
    !
    ncld      = 1        ! # of bottom model levels with no clouds (rad only)
    rhmax     = 1.01_rkx   ! RH at whicn FCC = 1.0
    rhmin     = 0.01_rkx   ! RH min value
    rh0land   = 0.80_rkx   ! Relative humidity threshold for land
    rh0oce    = 0.90_rkx   ! Relative humidity threshold for ocean
    tc0       = 238.0_rkx  ! Below this temp, rh0 begins to approach unity
    cllwcv    = 0.3e-3_rkx ! Cloud liquid water content for convective precip.
    clfrcvmax = 0.75_rkx   ! Max cloud fractional cover for convective precip.
    cftotmax  = 0.75_rkx   ! Max total cover cloud fraction for radiation
    k2_const  = 500.0_rkx  ! K2 CF factor relation with updraft mass flux
    kfac_shal = 0.07_rkx   ! Conv. cf factor in relation with updraft mass flux
    kfac_deep = 0.14_rkx   ! Conv. cf factor in relation with updraft mass flux
    lsrfhack  = .false.    ! Surface radiation hack
    larcticcorr = .true.   ! Vavrus and Waliser Arctic cloud correction
    rcrit     = 13.5_rkx   ! Mean critical radius
    coef_ccn  = 2.5e+20_rkx ! Coefficient determined by assuming a lognormal PMD
    abulk     = 0.9_rkx    ! Bulk activation ratio
    !
    ! microparam ;
    ! From original Nogerotto settings
    !
    stats = .false.
    budget_compute = .false. ! Verify enthalpy and moisture conservation
    nssopt = 1 ! Supersaturation Computation
               ! 0 => No scheme
               ! 1 => Tompkins
               ! 2 => Lohmann and Karcher
               ! 3 => Gierens
    iautoconv = 4 !  Choose the autoconversion paramaterization
                  ! => 1 Klein & Pincus (2000)
                  ! => 2 Khairoutdinov and Kogan (2000)
                  ! => 3 Kessler (1969)
                  ! => 4 Sundqvist
    vfqr = 4.0_rkx
    vfqi = 0.01_rkx
    vfqs = 1.0_rkx
    auto_rate_khair = 0.355_rkx
    auto_rate_kessl = 1.e-3_rkx
    auto_rate_klepi = 0.5e-3_rkx
    rkconv = 1.666e-4_rkx ! 1.0/6000.0
    skconv = 1.0e-3_rkx
    rcldiff = 1.0e-6_rkx
    rcovpmin = 0.1_rkx
    rpecons = 5.547e-5_rkx
    !
    ! grellparam ;
    ! Taken from MM5 Grell implementation
    !
    igcc        = 2        ! Closure scheme
    edtmin      = 0.20_rkx ! Minimum Precipitation Efficiency land
    edtmin_ocn  = 0.20_rkx ! Minimum Precipitation Efficiency ocean
    edtmax      = 0.80_rkx ! Maximum Precipitation Efficiency land
    edtmax_ocn  = 0.80_rkx ! Maximum Precipitation Efficiency ocean
    edtmino     = 0.20_rkx ! Minimum Tendency Efficiency (o var) land
    edtmino_ocn = 0.20_rkx ! Minimum Tendency Efficiency (o var) ocean
    edtmaxo     = 0.80_rkx ! Maximum Tendency Efficiency (o var) land
    edtmaxo_ocn = 0.80_rkx ! Maximum Tendency Efficiency (o var) ocean
    edtminx     = 0.20_rkx ! Minimum Tendency Efficiency (x var) land
    edtminx_ocn = 0.20_rkx ! Minimum Tendency Efficiency (x var) ocean
    edtmaxx     = 0.80_rkx ! Maximum Tendency Efficiency (x var) land
    edtmaxx_ocn = 0.80_rkx ! Maximum Tendency Efficiency (x var) ocean
    shrmin      = 0.30_rkx ! Minimum Shear effect on precip eff. land
    shrmin_ocn  = 0.30_rkx ! Minimum Shear effect on precip eff. ocean
    shrmax      = 0.90_rkx ! Maximum Shear effect on precip eff. land
    shrmax_ocn  = 0.90_rkx ! Maximum Shear effect on precip eff. ocean
    pbcmax =  50.0_rkx     ! Max depth (mb) of stable layer b/twn LCL & LFC
    mincld = 50.0_rkx      ! Min cloud depth (mb).
    htmin = -250.0_rkx     ! Min convective heating
    htmax = 500.0_rkx      ! Max convective heating
    skbmax = 0.4_rkx       ! Max cloud base height in sigma
    dtauc = 60.0_rkx ! Fritsch & Chappell (1980) ABE Removal Timescale (min)
    !
    ! emanparam ;
    ! From Kerry Emanuel convect 4.3c original code
    !
    minorig = 1
    elcrit_ocn = 0.0011_rkx ! Autoconversion threshold water content (gm/gm)
    elcrit_lnd = 0.0011_rkx ! Autoconversion threshold water content (gm/gm)
    tlcrit = -55.0_rkx    ! Below tlcrit auto-conversion threshold is zero
    entp = 0.50_rkx       ! Coefficient of mixing in the entrainment formulation
    sigd = 0.05_rkx       ! Fractional area covered by unsaturated dndraft
    sigs = 0.12_rkx       ! Fraction of precipitation falling outside of cloud
    omtrain = 50.0_rkx    ! Fall speed of rain (P/s)
    omtsnow = 5.5_rkx     ! Fall speed of snow (P/s)
    coeffr = 1.0_rkx      ! Coefficient governing the rate of rain evaporation
    coeffs = 0.8_rkx      ! Coefficient governing the rate of snow evaporation
    cu = 0.7_rkx          ! Coefficient governing convective momentum transport
    betae = 10.0_rkx      ! Controls downdraft velocity scale
    dtmax = 0.90_rkx    ! Max negative parcel temperature perturbation below LFC
    alphae = 0.01_rkx   ! Controls the approach rate to quasi-equilibrium
    damp = 0.1_rkx      ! Controls the approach rate to quasi-equilibrium
    epmax_ocn = 0.999_rkx ! Maximum precipitation efficiency over land
    epmax_lnd = 0.999_rkx ! Maximum precipitation efficiency over ocean
    !
    ! tiedtkeparam ;
    ! Taken from MPI Echam settings
    !
    iconv    = 4            ! Selects the actual scheme
    entrmax  = 1.75e-3_rkx  ! Max entrainment iconv=[1,2,3]
    entrdd   = 3.0e-4_rkx   ! Entrainment rate for cumulus downdrafts
    entrpen_lnd  = 1.75e-3_rkx  ! Entrainment rate for penetrative convection
    entrpen_ocn  = 1.75e-3_rkx  ! Entrainment rate for penetrative convection
    entrscv  = 3.0e-4_rkx   ! Entrainment rate for shallow convn iconv=[1,2,3]
    entrmid  = 1.0e-4_rkx   ! Entrainment rate for midlevel convn iconv=[1,2,3]
    cprcon   = 1.0e-4_rkx   ! Coefficient for determine conversion iconv=[1,2,3]
    detrpen_lnd = 0.75e-4_rkx   ! Detrainment rate for penetrative convection
    detrpen_ocn = 0.75e-4_rkx   ! Detrainment rate for penetrative convection
    entshalp = 2.0_rkx      ! shallow entrainment factor for entrpen
    rcuc_lnd = 0.05_rkx     ! Convective cloud cover for rain evporation
    rcuc_ocn = 0.05_rkx     ! Convective cloud cover for rain evporation
    rcpec_lnd = 5.55e-5_rkx ! Coefficient for rain evaporation below cloud
    rcpec_ocn = 5.55e-5_rkx ! Coefficient for rain evaporation below cloud
    rhebc_lnd = 0.8_rkx     ! Critical relative humidity below
                            ! cloud at which evaporation starts for land
    rhebc_ocn = 0.8_rkx     ! Critical relative humidity below
                            ! cloud at which evaporation starts for ocean
    rprc_lnd = 1.4e-3_rkx   ! coefficient for conversion from cloud water
    rprc_ocn = 1.4e-3_rkx   ! coefficient for conversion from cloud water
    revap_lnd = 1.0e-5_rkx  ! coefficient for evaporation over land
    revap_ocn = 1.0e-5_rkx  ! coefficient for evaporation over ocean
    cmtcape = 3600.0_rkx   ! CAPE adjustment timescale
    lmfpen    = .true.  ! penetrative conv is switched on
    lmfmid    = .true.  ! midlevel conv is switched on
    lmfdd     = .true.  ! cumulus downdraft is switched on
    lepcld    = .true.  ! prognostic cloud scheme is on
    lmfdudv   = .true.  ! cumulus friction is switched on
    lmfscv    = .true.  ! shallow convection is switched on
    lmfuvdis  = .true.  ! use kinetic energy dissipation
    lmftrac   = .true.  ! chemical tracer transport is on
    lmfsmooth = .false. ! smoot of mass fluxes for tracers
    lmfwstar  = .false. ! Grant w* closure for shallow conv
    !
    ! kfparam ;
    ! Taken from WRF KFeta parametrization
    !
    kf_wthreshold = 0.02_rkx ! Kain Fritsch vertical velocity threshold
    kf_entrate = 0.03_rkx    ! Kain Fritsch entrainment rate
    kf_convrate = 0.03_rkx   ! Kain Fritsch condensate to rain conversion rate
    kf_min_pef = 0.2_rkx  ! Minimum precipitation efficiency
    kf_max_pef = 0.9_rkx  ! Maximum precipitation efficiency
    kf_dpp = 150.0_rkx    ! Starting height of downdraft above updraft (mb)
    kf_tkemax = 5.0_rkx   ! Maximum turbolent kinetic energy in sub cloud layer
    kf_min_dtcape = 1800.0_rkx ! Consumption time of CAPE low limit
    kf_max_dtcape = 3600.0_rkx ! Consumption time of CAPE high limit
    !
    ! uwparam ;
    ! Original settings from Travis O'Brian
    !
    iuwvadv = 0
    atwo  = 10.0_rkx
    rstbl = 1.5_rkx
    czero = 5.869_rkx
    nuk   = 5.0_rkx
    !
    ! holtslagparam ;
    ! Settings from C. Torma
    !
    ricr_ocn = 0.25_rkx
    ricr_lnd = 0.25_rkx
    zhnew_fac = 0.25_rkx
    ifaholt = 0
    ifaholtth10 = 2
    holtth10iter = 1
    !
    ! slabocparam ;
    !
    mixed_layer_depth     = 50.0_rkx
    sst_restore_timescale = 5.0_rkx !days
    do_restore_sst = .true.
    do_qflux_adj = .false.
    !
    ! tweakparam ;
    !
    itweak = 0
    itweak_sst = 0
    itweak_temperature = 0
    itweak_solar_irradiance = 0
    itweak_greenhouse_gases = 0
    sst_tweak = 0.0_rkx
    temperature_tweak = 0.0_rkx
    solar_tweak = 0.0_rkx
    gas_tweak_factors(:) = 1.0_rkx
    !
    ! chemparam ; ( 0 = none, 1 = activated)
    !
    ichsolver = 1     ! enable chem solver
    ismoke = 0        ! consider emissions from fires (smoke tracer)
    ichremlsc = 1     ! tracer removal by large scale clouds
    ichremcvc = 1     ! tracer removal by convective clouds
    ichdrdepo = 1     ! tracer dry deposition
    ichcumtra = 1     ! tracer convective transport
    ichdustemd = 1    ! dust emission distribution (1 = alfaro, 2 =kok)
    ichdustparam = 1  ! read dust emission scheme surface parameters
    ichjphcld = 1     ! impact of cloud aod on photolysis coef
    ichecold = 0      ! chemistry cold start (restart without chem data in SAV)
    idirect = 0       ! tracer direct effect
#ifdef CLM45
    isnowdark = 1     ! Snow darkening by CARB/DUST
#else
    isnowdark = 0     ! Snow darkening by CARB/DUST
#endif
    iindirect = 0
    ichdiag = 0       ! chem tend outputs
    ichsursrc = 1
    ichebdy = 1
    ichlinox = 1
    rdstemfac = d_one
    ichbion = 0
    rocemfac = 1.33_rkx
    carb_aging_control = .false.

    ntr = 0
    nbin = 0
    igaschem = 0
    iaerosol = 0
    iisoropia = 0
    ioxclim = 0

#ifdef CLM
    !
    ! clmparam ; (read in case clm surface model compiled in)
    !
    imask = 1
    ilawrence_albedo = 1
    clmfrq = 24.0_rkx
#endif
    !
    ! cplparam ;
    !
    cpldt = 21600.0_rkx  ! coupling time step in seconds (seconds)
    zomax = 0.02_rkx     ! maximum allowed surface roughness from wave comp.
    ustarmax = 0.02_rkx  ! maximum allowed friction velocity from wave comp.

#ifdef CLM
    if ( myid == italk ) then
      if (nsg /= 1 ) then
        write (stderr,*) 'Running SUBGRID with CLM: not implemented'
        write (stderr,*) 'Please set nsg to 1 in regcm.in'
        call fatal(__FILE__,__LINE__, &
                   'CLM & SUBGRID TOGETHER')
      end if
    end if
#endif

    if ( myid == iocpu ) then
      open(newunit=ipunit, file=namelistfile, status='old', &
                   action='read', iostat=iretval)
      if ( iretval /= 0 ) then
        write(stderr,*) 'Error opening input namelist file ',trim(namelistfile)
        call fatal(__FILE__,__LINE__, &
                   'INPUT NAMELIST OPEN ERROR')
#ifdef DEBUG
      else
        write(stdout,*) 'Open ',trim(namelistfile),' OK'
#endif
      end if

      write(stdout,*) 'Reading model namelist file'

      rewind(ipunit)
      read (ipunit, nml=restartparam, iostat=iretval, err=100)
      if ( iretval /= 0 ) then
        write(stderr,*) 'Error reading restartparam namelist'
        call fatal(__FILE__,__LINE__, &
                   'INPUT NAMELIST READ ERROR')
#ifdef DEBUG
      else
        write(stdout,*) 'Read restartparam OK'
#endif
      end if

      idate0 = i8wcal(mdate0,ical)
      idate1 = i8wcal(mdate1,ical)
      idate2 = i8wcal(mdate2,ical)
      bdif = idate2 - idate1
      hspan = nint(tohours(bdif))
      if ( mod(hspan,24) /= 0 ) then
        call fatal(__FILE__,__LINE__,  &
                   'Runtime increments must be modulus 24 hours')
      end if

      rewind(ipunit)
      read (ipunit, nml=timeparam, iostat=iretval, err=101)
      if ( iretval /= 0 ) then
        write(stderr,*) 'Error reading timeparam namelist'
        call fatal(__FILE__,__LINE__, &
                   'INPUT NAMELIST READ ERROR')
#ifdef DEBUG
      else
        write(stdout,*) 'Read timeparam OK'
#endif
      end if

      rewind(ipunit)
      read (ipunit, nml=outparam, iostat=iretval, err=102)
      if ( iretval /= 0 ) then
        write(stderr,*) 'Error reading outparam namelist'
        call fatal(__FILE__,__LINE__, &
                   'INPUT NAMELIST READ ERROR')
#ifdef DEBUG
      else
        write(stdout,*) 'Read outparam OK'
#endif
      end if

      len_path = len(trim(dirout))
      if ( dirout(len_path:len_path) /= '/' ) dirout = trim(dirout)//'/'
      rewind(ipunit)
      read (ipunit, nml=physicsparam, iostat=iretval, err=103)
      if ( iretval /= 0 ) then
        write(stderr,*) 'Error reading physicsparam namelist'
        call fatal(__FILE__,__LINE__, &
                   'INPUT NAMELIST READ ERROR')
#ifdef DEBUG
      else
        write(stdout,*) 'Read physicsparam OK'
#endif
      end if

      if ( idynamic < 3 ) then
        upstream_mode = .true.
        stability_enhance = .true.
        if ( idynamic == 2 ) then
          gnu1 = 0.1000_rkx
          gnu2 = 0.1000_rkx
          diffu_hgtf = 0
        else if ( idynamic == 1 ) then
          gnu1 = 0.0625_rkx
          gnu2 = 0.0625_rkx
          diffu_hgtf = 1
        end if
        ckh = 1.0_rkx
        adyndif = 1.0_rkx
        uoffc = 0.250_rkx
        t_extrema = 5.0_rkx
        q_rel_extrema = 0.20_rkx
        rewind(ipunit)
        read (ipunit, nml=dynparam, iostat=iretval, err=104)
        if ( iretval /= 0 ) then
          write(stdout,*) 'Using default dynamical parameters.'
#ifdef DEBUG
        else
          write(stdout,*) 'Read dynparam OK'
#endif
        end if

        if ( idynamic == 2 ) then
          rewind(ipunit)
          read (ipunit, nml=nonhydroparam, iostat=iretval, err=105)
          if ( iretval /= 0 ) then
            write(stdout,*) 'Using default non-hydrostatc parameters.'
#ifdef DEBUG
          else
            write(stdout,*) 'Read nonhydroparam OK'
#endif
          end if
        else if ( idynamic == 1 ) then
          rewind(ipunit)
          read (ipunit, nml=hydroparam, iostat=iretval, err=106)
          if ( iretval /= 0 ) then
            write(stdout,*) 'Using default hydrostatc parameters.'
#ifdef DEBUG
          else
            write(stdout,*) 'Read hydroparam OK'
#endif
          end if
        end if
      else
        rewind(ipunit)
        read (ipunit, nml=nonhydroparam, iostat=iretval, err=105)
        if ( iretval /= 0 ) then
          write(stdout,*) 'Using default non-hydrostatc parameters.'
#ifdef DEBUG
        else
          write(stdout,*) 'Read nonhydroparam OK'
#endif
        end if
      end if

      ! Hack. permanently disable seasonal albedo.
      idesseas = 0

      icup(1) = icup_lnd
      icup(2) = icup_ocn
      if ( any(icup == 1) ) then
        if ( idynamic == 3 ) then
          write(stderr,*) &
            'ERROR: Kuo scheme cannot be used with MOLOCH dynamical scheme'
          call fatal(__FILE__,__LINE__, &
                     'INPUT NAMELIST ICUP INCONSISTENT')
        end if
        if ( icup_lnd /= icup_ocn ) then
          write(stderr,*) &
            'ERROR: Kuo scheme MUST be used on both Land and Ocean'
          call fatal(__FILE__,__LINE__, &
                     'INPUT NAMELIST ICUP INCONSISTENT')
        end if
      end if
      if ( any(icup == 3) ) then
        if ( icup_lnd /= icup_ocn ) then
          write(stderr,*) 'ERROR: BM scheme MUST be used on both Land and Ocean'
          call fatal(__FILE__,__LINE__, &
                     'INPUT NAMELIST ICUP INCONSISTENT')
        end if
      end if

      rewind(ipunit)
      read (ipunit, nml=cldparam, iostat=iretval, err=107)
      if ( iretval /= 0 ) then
        write(stdout,*) 'Using default cloud parameter.'
#ifdef DEBUG
      else
        write(stdout,*) 'Read cldparam OK'
#endif
      end if
      if ( cftotmax < 0.0 ) then
        cftotmax = 0.1_rkx
      else if ( cftotmax > 1.0_rkx ) then
        cftotmax = 1.00_rkx
      end if

      if ( ipptls == 1 ) then
        rewind(ipunit)
        read (ipunit, nml=subexparam, iostat=iretval, err=108)
        if ( iretval /= 0 ) then
          write(stdout,*) 'Using default subex parameter.'
#ifdef DEBUG
        else
          write(stdout,*) 'Read subexparam OK'
#endif
        end if
      else if ( ipptls == 2 ) then
        rewind(ipunit)
        read (ipunit, nml=microparam, iostat=iretval, err=109)
        if ( iretval /= 0 ) then
          write(stdout,*) 'Using default microphysical parameter.'
#ifdef DEBUG
        else
          write(stdout,*) 'Read microparam OK'
#endif
        end if
        if ( budget_compute ) then
          write(stdout,*) 'Will check the total enthalpy and moisture'
        end if
      end if

      if ( any(icup == 2) ) then
        rewind(ipunit)
        read (ipunit, nml=grellparam, iostat=iretval, err=110)
        if ( iretval /= 0 ) then
          write(stdout,*) 'Using default Grell parameter.'
#ifdef DEBUG
        else
          write(stdout,*) 'Read grellparam OK'
#endif
        end if
      end if
      if ( any(icup == 4) ) then
        rewind(ipunit)
        read (ipunit, nml=emanparam, iostat=iretval, err=111)
        if ( iretval /= 0 ) then
          write(stdout,*) 'Using default MIT parameter.'
#ifdef DEBUG
        else
          write(stdout,*) 'Read emanparam OK'
#endif
        end if
        if ( istochastic == 1 ) then
          rewind(ipunit)
          read (ipunit, nml=emanstochastic, iostat=iretval, err=111)
          if ( iretval /= 0 ) then
            write(stdout,*) 'No emanstochastic namelist.'
            write(stdout,*) 'Setting Emanuel to deterministic.'
            istochastic = 0
#ifdef DEBUG
          else
            write(stdout,*) 'Activate Emanuel Stochastic parameterization.'
#endif
          end if
        end if
      end if
      if ( any(icup == 5) ) then
        rewind(ipunit)
        read (ipunit, nml=tiedtkeparam, iostat=iretval, err=112)
        if ( iretval /= 0 ) then
          write(stdout,*) 'Using default Tiedtke parameter.'
#ifdef DEBUG
        else
          write(stdout,*) 'Read tiedtkeparam OK'
#endif
        end if
      end if
      if ( any(icup == 6) ) then
        rewind(ipunit)
        read (ipunit, nml=kfparam, iostat=iretval, err=113)
        if ( iretval /= 0 ) then
          write(stdout,*) 'Using default Kain Fritsch parameter.'
#ifdef DEBUG
        else
          write(stdout,*) 'Read kfparam OK'
#endif
        end if
        if ( kf_min_dtcape < 600.0_rkx ) then
          write(stdout,*) 'Resetting kf_min_dtcape to 600 s'
          kf_min_dtcape = 600.0_rkx
        end if
        if ( kf_max_dtcape > 7200.0_rkx ) then
          write(stdout,*) 'Resetting kf_max_dtcape to 7200 s'
          kf_max_dtcape = 7200.0_rkx
        end if
        if ( kf_tkemax > 12.0_rkx ) then
          write(stdout,*) 'Resetting kf_tkemax to 12 m2 s-2'
          kf_tkemax = 12.0_rkx
        end if
        if ( kf_tkemax < 3.0_rkx ) then
          write(stdout,*) 'Resetting kf_tkemax to 3 m2 s-2'
          kf_tkemax = 3.0_rkx
        end if
      end if
      if ( iocnflx < 1 .or. iocnflx > 3 ) then
        call fatal(__FILE__,__LINE__, &
                   'UNSUPPORTED OCEAN FLUX SCHEME.')
      end if
      if ( any(icup < -1) .or. any(icup > 6) ) then
        call fatal(__FILE__,__LINE__, &
                   'UNSUPPORTED CUMULUS SCHEME')
      end if
      if ( ibltyp < 0 .or. ibltyp > 4 ) then
        call fatal(__FILE__,__LINE__, &
                   'UNSUPPORTED PBL SCHEME.')
      end if
#ifdef CLM
      if ( ibltyp > 2 ) then
        call fatal(__FILE__,__LINE__, &
                   'UNSUPPORTED PBL SCHEME FOR CLM 3.5')
      end if
#endif
      if ( ibltyp == 1 ) then
        rewind(ipunit)
        read (ipunit, nml=holtslagparam, iostat=iretval, err=114)
        if ( iretval /= 0 ) then
          write(stdout,*) 'Using default Holtslag parameter.'
#ifdef DEBUG
        else
          write(stdout,*) 'Read holtslagparam OK'
#endif
        end if
      end if
      if ( ibltyp == 2 ) then
        rewind(ipunit)
        read (ipunit, nml=uwparam, iostat=iretval, err=115)
        if ( iretval /= 0 ) then
          write(stdout,*) 'Using default UW PBL parameter.'
#ifdef DEBUG
        else
          write(stdout,*) 'Read uwparam OK'
#endif
        end if
      end if
      if ( irrtm == 1 ) then
        rewind(ipunit)
        read (ipunit, nml=rrtmparam, iostat=iretval, err=116)
        if ( iretval /= 0 ) then
          write(stdout,*) 'Using default RRTM parameter.'
#ifdef DEBUG
        else
          write(stdout,*) 'Read rrtmparam OK'
#endif
        end if
      end if

      if ( islab_ocean == 1 ) then
        rewind(ipunit)
        read (ipunit, nml=slabocparam, iostat=iretval, err=117)
        if ( iretval /= 0 ) then
          write(stdout,*) 'Using default SLAB Ocean parameter.'
#ifdef DEBUG
        else
          write(stdout,*) 'Read slabocparam OK'
#endif
        end if
        if ( do_qflux_adj .eqv. do_restore_sst ) then
          write (stderr,*) 'do_qflux_adj   = ' , do_qflux_adj
          write (stderr,*) 'do_restore_sst = ' , do_restore_sst
          write (stderr,*) 'THESE OPTION CANNOT BE EQUAL !!'
          write (stderr,*) 'FIRST DO A RESTORE SST RUN AND THEN AN ADJUST RUN!'
          call fatal(__FILE__,__LINE__, &
                     'SLABOCEAN INPUT INCONSISTENCY')
        end if
        if ( idcsst == 1 ) then
          write(stderr,*) 'The SLAB Ocean model disables the diurnal SST'
          idcsst = 0
        end if
        if ( iseaice == 1 ) then
          write(stderr,*) 'The SLAB Ocean model disables the SeaIce model'
          iseaice = 0
        end if
      end if

      if ( ichem == 1 ) then
        !if ( iclimaaer == 1 ) then
        !  write(stderr,*) 'Cannot define both ichem and iclimaaer'
        !  call fatal(__FILE__,__LINE__, &
        !             'INPUT NAMELIST INCONSISTENT AEROSOLS')
        !end if
        rewind(ipunit)
        read (ipunit, chemparam, iostat=iretval, err=118)
        if ( iretval /= 0 ) then
          write(stderr,*) 'Error reading chemparam namelist'
          call fatal(__FILE__,__LINE__, &
                     'INPUT NAMELIST READ ERROR')
#ifdef DEBUG
        else
          write(stdout,*) 'Read chemparam OK'
#endif
        end if
      end if
#ifndef CLM45
      if ( isnowdark > 0 ) then
        write(stderr,*) 'Snow Darkening effect active only with CLM45'
        write(stderr,*) 'Reset it to zero.'
        isnowdark = 0
      end if
#endif
#ifdef CLM
      rewind(ipunit)
      read (ipunit , clmparam, iostat=iretval, err=119)
      if ( iretval /= 0 ) then
        write(stdout,*) 'Using default CLM parameter.'
#ifdef DEBUG
      else
        write(stdout,*) 'Read clmparam OK'
#endif
      end if
#endif
      if ( iocncpl == 1 .or. iwavcpl == 1 .or. icopcpl == 1 ) then
        rewind(ipunit)
        read (ipunit , cplparam, iostat=iretval, err=120)
        if ( iretval /= 0 ) then
          write(stdout,*) 'Using default Coupling parameter.'
#ifdef DEBUG
        else
          write(stdout,*) 'Read cplparam OK'
#endif
        end if
      end if

      if ( itweak == 1 ) then
        rewind(ipunit)
        read (ipunit , tweakparam, iostat=iretval, err=121)
        if ( iretval /= 0 ) then
          write(stdout,*) 'Tweak parameters absent.'
          write(stdout,*) 'Disable tweaking.'
          itweak = 0
#ifdef DEBUG
        else
          write(stdout,*) 'Read tweakparam OK'
#endif
        end if
        if ( itweak_sst == 0 .and.              &
             itweak_temperature == 0 .and.      &
             itweak_solar_irradiance == 0 .and. &
             itweak_greenhouse_gases == 0 ) then
          write(stdout,*) 'Tweak parameters not enabled.'
          write(stdout,*) 'Disable tweaking.'
          itweak = 0
        end if
      end if

      close(ipunit)

      if ( dt < mindt ) then
        write(stderr,*) 'Minimum dt allowed is ',mindt,' seconds.'
        call fatal(__FILE__,__LINE__, &
                   'DT TOO SMALL')
      end if

      dt = check_against_outparams(dt,mindt)

      minfrq = 86400.0_rkx

      ! Check user input
      if ( srffrq <= 0.0_rkx ) srffrq = 3.0_rkx
      if ( atmfrq <= 0.0_rkx ) atmfrq = 6.0_rkx
      if ( radfrq <= 0.0_rkx ) radfrq = 6.0_rkx
      if ( optfrq <= 0.0_rkx ) optfrq = 6.0_rkx
      if ( lakfrq <= 0.0_rkx ) lakfrq = 6.0_rkx
      if ( subfrq <= 0.0_rkx ) subfrq = 3.0_rkx
      if ( chemfrq <= 0.0_rkx ) chemfrq = 6.0_rkx

      if ( ifsrf ) minfrq = min(minfrq,max(srffrq*3600.0_rkx,dt))
      if ( ifatm ) minfrq = min(minfrq,max(atmfrq*3600.0_rkx,dt))
      if ( ifrad ) minfrq = min(minfrq,max(radfrq*3600.0_rkx,dt))
      if ( ifopt ) minfrq = min(minfrq,max(optfrq*3600.0_rkx,dt))
      if ( ifshf ) minfrq = min(minfrq,3600.0_rkx)
      if ( ichem == 1 ) then
        if ( ifchem ) minfrq = min(minfrq,max(chemfrq*3600.0_rkx,dt))
      end if
      if ( lakemod == 1 ) then
        if ( iflak ) minfrq = min(minfrq,max(lakfrq*3600.0_rkx,dt))
      end if
      if ( nsg > 1 ) then
        if ( ifsub ) minfrq = min(minfrq,max(subfrq*3600.0_rkx,dt))
      end if

      if ( dtsrf <= 0.0_rkx ) dtsrf = 600.0_rkx
      if ( dtcum <= 0.0_rkx ) dtcum = 300.0_rkx
      if ( dtche <= 0.0_rkx ) dtche = 900.0_rkx
      if ( dtrad <= 0.0_rkx ) then
        dtrad = 1800.0_rkx
      else
        dtrad = dtrad * 60.0_rkx
      end if
      if ( dtabem <= 0.0_rkx ) then
        dtabem = 64800.0_rkx
      else
        dtabem = dtabem * 3600.0_rkx
      end if

      if ( dtcum < dt ) dtcum = dt
      if ( dtsrf < dt ) dtsrf = dt
      if ( dtche < dt ) dtche = dt
      if ( dtrad < dt ) dtrad = dt
      if ( dtabem < dt ) dtabem = dt

      dtsrf = int(dtsrf / dt) * dt
      do while ( mod(minfrq,dtsrf) > d_zero )
        dtsrf = dtsrf - dt
      end do

      dtcum = int(dtcum / dt) * dt
      do while ( mod(minfrq,dtcum) > d_zero )
        dtcum = dtcum - dt
      end do
      dtcum = max(int(dtcum / (0.5_rkx*dtsrf)),1) * (0.5_rkx*dtsrf)

      dtrad = int(dtrad / dt) * dt
      do while ( mod(minfrq,dtrad) > d_zero )
        dtrad = dtrad - dt
      end do
      dtrad = max(int(dtrad / (3.0_rkx*dtsrf)),1) * (3.0_rkx*dtsrf)

      dtabem = max(int(dtabem / (36_rkx*dtrad)),1) * (36.0_rkx*dtrad)

      dtche = int(dtche / dt) * dt

      if ( iseaice == 1 ) then
        select case (ssttyp)
          case ('EIN15','EIN75','EIXXX','ERA5 ')
            icetriggert = 271.465_rkx
          case default
            icetriggert = 271.355_rkx
        end select
      end if
    end if
    !
    ! communicate to all processors
    !
    call bcast(lsmoist)
    call bcast(ifrest)
    call bcast(hspan)
    call bcast(idate0)
    call bcast(idate1)
    call bcast(idate2)
    call bcast(globidate1)
    call bcast(globidate2)

    call bcast(dt)
    call bcast(dtrad)
    call bcast(dtsrf)
    call bcast(dtcum)
    call bcast(dtche)
    call bcast(dtabem)

    call bcast(prestr,64)
    call bcast(ifsave)
    call bcast(ifatm)
    call bcast(ifrad)
    call bcast(ifshf)
    call bcast(ifsrf)
    call bcast(ifsub)
    call bcast(iflak)
    call bcast(ifsts)
    call bcast(ifopt)
    call bcast(ifchem)
    call bcast(outnwf)
    call bcast(savfrq)
    call bcast(atmfrq)
    call bcast(radfrq)
    call bcast(srffrq)
    call bcast(lakfrq)
    call bcast(subfrq)
    call bcast(chemfrq)
    call bcast(optfrq)
    call bcast(enable_atm_vars)
    call bcast(enable_rad_vars)
    call bcast(enable_srf_vars)
    call bcast(enable_shf_vars)
    call bcast(enable_sub_vars)
    call bcast(enable_sts_vars)
    call bcast(enable_lak_vars)
    call bcast(enable_opt_vars)
    call bcast(enable_che_vars)
    call bcast(lsync)
    call bcast(uvrotate)
    call bcast(idiag)
    call bcast(icosp)
    call bcast(do_parallel_netcdf_in)
    call bcast(do_parallel_netcdf_out)
    call bcast(chechgact)
#ifdef NETCDF4_HDF5
    call bcast(deflate_level)
#endif

    do_parallel_save = (do_parallel_netcdf_in .and. do_parallel_netcdf_out)

    ! Reset the NEEDED 2D vars.
    enable_atm_vars(1:6) = .true.
    enable_rad_vars(1:6) = .true.
    enable_opt_vars(1:6) = .true.
    enable_che_vars(1:6) = .true.
    ! These do not have p0, no vertical field.
    enable_srf_vars(1:5) = .true.
    enable_shf_vars(1:5) = .true.
    enable_sub_vars(1:5) = .true.
    enable_sts_vars(1:5) = .true.
    enable_lak_vars(1:5) = .true.

    call bcast(gnu1)
    call bcast(gnu2)
    call bcast(diffu_hgtf)
    call bcast(ckh)
    call bcast(adyndif)
    call bcast(upstream_mode)
    call bcast(uoffc)
    call bcast(stability_enhance)
    call bcast(t_extrema)
    call bcast(q_rel_extrema)

    call bcast(iboudy)
    call bcast(isladvec)
    call bcast(iqmsl)
    call bcast(ibltyp)
    call bcast(icup_lnd)
    call bcast(icup_ocn)
    call bcast(icup)
    call bcast(ipptls)
    call bcast(idiffu)
    call bcast(iocnflx)
    call bcast(iocncpl)
    call bcast(iwavcpl)
    call bcast(icopcpl)
    call bcast(iocnrough)
    call bcast(iocnzoq)
    call bcast(ipgf)
    call bcast(iemiss)
    call bcast(lakemod)
    call bcast(ichem)
    call bcast(iclimaaer)
    call bcast(radclimpath,256)

    if ( idynamic == 3 ) then
      if ( isladvec == 1 ) then
        write(stderr,*) 'Moloch core does not work SL advection code'
        isladvec = 0
      end if
      if ( any(icup ==  1) ) then
        write(stderr,*) 'Moloch core does not work with Kuo convection scheme'
        call fatal(__FILE__,__LINE__, &
                   'MOLOCH DOES NOT WORK WITH KUO')
      end if
      ! Moloch paramters here
      call bcast(mo_anu2)
      call bcast(mo_wmax)
      call bcast(mo_nzfilt)
      call bcast(mo_nadv)
      call bcast(mo_nsound)
      call bcast(ifrayd)
      call bcast(rayndamp)
      call bcast(rayalpha0)
      call bcast(rayhd)
    else if ( idynamic == 2 ) then
      call bcast(base_state_pressure)
      call bcast(logp_lrate)
      call bcast(ifupr)
      call bcast(nhbet)
      call bcast(nhxkd)
      call bcast(itopnudge)
      call bcast(ifrayd)
      call bcast(rayndamp)
      call bcast(rayalpha0)
      call bcast(rayhd)
    else
      call bcast(nsplit)
      call bcast(lstand)
    end if

    if ( iboudy == 4 ) then
      nspgd = max(6,nspgd)
      nspgx = max(5,nspgx)
    end if

    ! Check if really do output

#ifdef CLM
    if ( lakemod /= 0 ) then
      if ( myid == italk ) then
        write(stderr,*) 'Disabling BATS lake model, this is a CLM run'
      end if
      lakemod = 0
    end if
#endif
#if defined(CLM45) || defined(CLM)
    if ( iemiss /= 0 ) then
      if ( myid == italk ) then
        write(stderr,*) 'Using CLM Radiant Temperature'
      end if
      iemiss = 0
    end if
#endif
    if ( lakemod /= 1 ) then
      iflak = .false.
    end if
    if ( nsg < 2 ) then
      ifsub = .false.
    end if

    if ( ichem /= 1 ) then
      ifchem = .false.
      if (iclimaaer == 0) ifopt = .false.
    end if

    !
    ! Force the correct scenario from dattyp in CMIP5
    !
    if ( myid == iocpu ) then
      if ( dattyp(4:5) == '26' ) then
        if ( scenario /= 'RCP3PD' .or. &
             scenario /= 'RCP2.6' .or. &
             scenario /= 'RCP26' ) then
          write(stderr,*) 'Forcing scenario from dattyp to RCP2.6'
          scenario = 'RCP2.6'
        end if
      else if ( dattyp(4:5) == '45' ) then
        if ( scenario /= 'RCP4.5' .or. scenario /= 'RCP45' ) then
          write(stderr,*) 'Forcing scenario from dattyp to RCP4.5'
          scenario = 'RCP4.5'
        end if
      else if ( dattyp(4:5) == '60' ) then
        if ( scenario /= 'RCP6' .or.  &
             scenario /= 'RCP60' .or. &
             scenario /= 'RCP6.0' ) then
          write(stderr,*) 'Forcing scenario from dattyp to RCP6.0'
          scenario = 'RCP6.0'
        end if
      else if ( dattyp(4:5) == '85' ) then
        if ( scenario /= 'RCP8.5' .or. scenario /= 'RCP85' ) then
          write(stderr,*) 'Forcing scenario from dattyp to RCP8.5'
          scenario = 'RCP8.5'
        end if
      end if
    end if

    if ( iocncpl == 1 .or. iwavcpl == 1 .or. icopcpl == 1 ) then
      call bcast(cpldt)
      call bcast(zomax)
      call bcast(ustarmax)
    end if

    call bcast(scenario,8)
    call bcast(ghg_year_const)
    call bcast(idcsst)
    call bcast(ipcpcool)
    call bcast(iwhitecap)
    call bcast(iseaice)
    call bcast(icetriggert)
    call bcast(idesseas)
    call bcast(iconvlwp)
    call bcast(icldfrac)
    call bcast(icldmstrat)
    call bcast(irrtm)
    call bcast(iclimao3)
    call bcast(isolconst)
    call bcast(ifixsolar)
    call bcast(fixedsolarval)
    call bcast(irceideal)
    if ( irceideal == 1 ) then
      if ( i_band /= 1 ) then
        call fatal(__FILE__,__LINE__, &
          'irceideal can run only if i_band has been set to 1')
      end if
      if ( i_crm /= 1 ) then
        call fatal(__FILE__,__LINE__, &
          'irceideal can run only if i_crm has been set to 1')
      end if
    end if
    call bcast(year_offset)
    call bcast(icumcloud)
    call bcast(islab_ocean)
    call bcast(itweak)
    if ( idcsst == 1 .and. iocnflx /= 2 ) then
      if ( myid == italk ) then
        write(stderr,*) &
          'Cannot enable diurnal cycle sst without Zheng ocean flux'
        write(stderr,*) 'Disabling idcsst.'
      end if
      idcsst = 0
    else if ( idcsst == 1 .and. iocncpl == 1 ) then
      if ( myid == italk ) then
        write(stderr,*) &
          'Cannot enable diurnal cycle sst with coupled ocean'
        write(stderr,*) 'Disabling idcsst.'
      end if
      idcsst = 0
    end if

#ifdef CLM
    call bcast(dirclm,256)
    call bcast(imask)
    call bcast(clmfrq)
    call bcast(ilawrence_albedo)
#endif

    !
    ! Cloud parameters
    !
    call bcast(ncld)
    call bcast(rhmax)
    call bcast(rhmin)
    call bcast(rh0oce)
    call bcast(rh0land)
    call bcast(tc0)
    call bcast(cftotmax)
    call bcast(clfrcvmax)
    call bcast(cllwcv)
    call bcast(rcrit)
    call bcast(coef_ccn)
    call bcast(abulk)
    call bcast(lsrfhack)
    call bcast(larcticcorr)
    call bcast(k2_const)
    call bcast(kfac_shal)
    call bcast(kfac_deep)

    if ( ipptls == 1 ) then
      call bcast(qck1land)
      call bcast(qck1oce)
      call bcast(gulland)
      call bcast(guloce)
      call bcast(cevaplnd)
      call bcast(cevapoce)
      call bcast(caccrlnd)
      call bcast(caccroce)
      call bcast(conf)
    else if ( ipptls == 2 ) then
      call bcast(stats)
      call bcast(budget_compute)
      call bcast(nssopt)
      call bcast(iautoconv)
      call bcast(vfqr)
      call bcast(vfqi)
      call bcast(vfqs)
      call bcast(auto_rate_khair)
      call bcast(auto_rate_kessl)
      call bcast(auto_rate_klepi)
      call bcast(rcovpmin)
      call bcast(rpecons)
      call bcast(rkconv)
      call bcast(skconv)
      call bcast(rcldiff)
    end if

    if ( ipptls > 1 ) then
      nqx = 5
      iqfrst = iqc
      iqlst  = iqs
    else
      nqx = 2
      iqfrst = iqc
      iqlst  = iqc
    end if

    if ( irrtm == 1 ) then
      call bcast(inflgsw)
      call bcast(iceflgsw)
      call bcast(liqflgsw)
      call bcast(inflglw)
      call bcast(iceflglw)
      call bcast(liqflglw)
      call bcast(icld)
      call bcast(irng)
      call bcast(imcica)
      call bcast(nradfo)
      if ( imcica == 0 .and. inflgsw == 2 ) then
        if ( myid == italk ) then
          write(stderr,*) &
            'Cannot use inflgsw RRTM cloud optical properties options '// &
            'when mcICA is not enabled imcica = 0'
          write(stderr,*) 'Running setting back inflgsw to zero.'
        end if
        inflgsw = 0
      end if
    end if

    if ( any(icup == 2) ) then
      call bcast(igcc)
      call bcast(shrmin)
      call bcast(shrmax)
      call bcast(edtmin)
      call bcast(edtmax)
      call bcast(edtmino)
      call bcast(edtmaxo)
      call bcast(edtminx)
      call bcast(edtmaxx)
      call bcast(shrmin_ocn)
      call bcast(shrmax_ocn)
      call bcast(edtmin_ocn)
      call bcast(edtmax_ocn)
      call bcast(edtmino_ocn)
      call bcast(edtmaxo_ocn)
      call bcast(edtminx_ocn)
      call bcast(edtmaxx_ocn)
      call bcast(pbcmax)
      call bcast(mincld)
      call bcast(htmin)
      call bcast(htmax)
      call bcast(skbmax)
      call bcast(dtauc)
    end if

    if ( any(icup == 4) ) then
      call bcast(minorig)
      call bcast(elcrit_ocn)
      call bcast(elcrit_lnd)
      call bcast(tlcrit)
      call bcast(entp)
      call bcast(sigd)
      call bcast(sigs)
      call bcast(omtrain)
      call bcast(omtsnow)
      call bcast(coeffr)
      call bcast(coeffs)
      call bcast(cu)
      call bcast(betae)
      call bcast(dtmax)
      call bcast(alphae)
      call bcast(damp)
      call bcast(epmax_ocn)
      call bcast(epmax_ocn)
      call bcast(istochastic)
      if ( istochastic == 1 ) then
        call bcast(sigs_min)
        call bcast(sigs_max)
        call bcast(sigd_min)
        call bcast(sigd_max)
        call bcast(epmax_lnd_min)
        call bcast(epmax_lnd_max)
        call bcast(elcrit_lnd_min)
        call bcast(elcrit_lnd_max)
      end if
    end if

    if ( any(icup == 5) ) then
      call bcast(iconv)
      call bcast(entrmax)
      call bcast(entrdd)
      call bcast(entrpen_lnd)
      call bcast(entrpen_ocn)
      call bcast(entrscv)
      call bcast(entrmid)
      call bcast(cprcon)
      call bcast(detrpen_lnd)
      call bcast(detrpen_ocn)
      call bcast(entshalp)
      call bcast(rcuc_lnd)
      call bcast(rcuc_ocn)
      call bcast(rcpec_lnd)
      call bcast(rcpec_ocn)
      call bcast(rhebc_lnd)
      call bcast(rhebc_ocn)
      call bcast(rprc_lnd)
      call bcast(rprc_ocn)
      call bcast(revap_lnd)
      call bcast(revap_ocn)
      call bcast(cmtcape)
      call bcast(lmfpen)
      call bcast(lmfmid)
      call bcast(lmfdd)
      call bcast(lepcld)
      call bcast(lmfdudv)
      call bcast(lmfscv)
      call bcast(lmfuvdis)
      call bcast(lmftrac)
      call bcast(lmfsmooth)
      call bcast(lmfwstar)
    end if

    if ( any(icup == 6) ) then
      call bcast(kf_wthreshold)
      call bcast(kf_entrate)
      call bcast(kf_convrate)
      call bcast(kf_min_pef)
      call bcast(kf_max_pef)
      call bcast(kf_dpp)
      call bcast(kf_min_dtcape)
      call bcast(kf_max_dtcape)
      call bcast(kf_tkemax)
    end if

    if ( ibltyp == 1 ) then
      call bcast(ricr_ocn)
      call bcast(ricr_lnd)
      call bcast(zhnew_fac)
      call bcast(ifaholtth10)
      call bcast(ifaholt)
      call bcast(holtth10iter)
    end if
    if ( ibltyp == 2 ) then
      call bcast(iuwvadv)
      call bcast(atwo)
      call bcast(rstbl)
      call bcast(czero)
      call bcast(nuk)
    end if

    if ( islab_ocean == 1 ) then
      call bcast(do_qflux_adj)
      call bcast(do_restore_sst)
      call bcast(sst_restore_timescale)
      call bcast(mixed_layer_depth)
      ! Save the input restore flux file for the adjust run
      if ( do_restore_sst ) ifslaboc = .true.
    end if

    if ( itweak == 1 ) then
      call bcast(itweak_sst)
      call bcast(itweak_temperature)
      call bcast(itweak_solar_irradiance)
      call bcast(itweak_greenhouse_gases)
      call bcast(sst_tweak)
      call bcast(temperature_tweak)
      call bcast(solar_tweak)
      call bcast(gas_tweak_factors)
    end if

    if ( ichem == 1 ) then
      call bcast(chemsimtype,8)
      call bcast(ichremlsc)
      call bcast(ichremcvc)
      call bcast(ichdrdepo)
      call bcast(ichcumtra)
      call bcast(idirect)
      call bcast(ichecold)
      call bcast(isnowdark)
      call bcast(iindirect)
      call bcast(ichsolver)
      call bcast(ichjphcld)
      call bcast(ichdustemd)
      call bcast(ichdustparam)
      call bcast(rdstemfac)
      call bcast(rocemfac)
      call bcast(ichdiag)
      call bcast(ichsursrc)
      call bcast(ichebdy)
      call bcast(ichlinox)
      call bcast(ichbion)
      call bcast(ismoke)
      call bcast(carb_aging_control)

      ! Set chemistry dimensions and tracer names
      call chem_config

    end if

    rcmtimer => rcm_timer(idate0,idate1,idate2,dt)

    if ( iclimaaer == 1 ) then
      call init_aerclima
    end if
    !
    ! ALLOCATE NEEDED SPACE
    !
    call allocate_mod_runparams

    call allocate_mod_atm_interface

    call allocate_mod_tend

    call allocate_mod_bdycon

    call allocate_pblscheme

    call allocate_micro

    call allocate_mod_split

    call allocate_mod_savefile

    call allocate_surface_model

    call allocate_radiation

    call allocate_mod_che_common
    call allocate_mod_che_mppio
    call allocate_mod_che_dust
    call allocate_mod_che_bdyco
    call allocate_mod_che_bionit

    if ( isladvec == 1 ) then
      call allocate_mod_sldepparam
    end if

    if ( idynamic == 2 ) then
      call allocate_mod_sound
    else if ( idynamic == 3 ) then
      call allocate_moloch
    end if

    call allocate_mod_diffusion

    if ( myid == italk ) then
      if ( ifsrf ) then
        if ( ifsts .and. srffrq > 24.0_rkx ) then
          call fatal(__FILE__,__LINE__, &
                     'NEED SRF FREQUENCY LESS THAN 24H FOR STS OUTPUT')
        end if
      else
        if ( ifsts ) then
          call fatal(__FILE__,__LINE__, &
                     'TO ENABLE STS, ENABLE SRF OUTPUT IS REQUIRED')
        end if
      end if
      if ( ichem == 1 ) then
        if ( chemfrq <= d_zero ) then
          write (stderr,*) 'CHEMFRQ=', chemfrq
          call fatal(__FILE__,__LINE__, &
                     'CHEMFRQ CANNOT BE ZERO')
        end if
      end if
      if ( isladvec == 1 ) then
        if ( jxp < 5 .or. iyp < 5 ) then
          write (stderr,*) 'To use Semi-Lagrangian Advection Scheme reduce'
          write (stderr,*) 'the number of processors !!!!'
          write (stderr,*) 'Minimum number of points is 25 (5x5) per processor'
          call fatal(__FILE__,__LINE__, &
                     'ISLADVEC WITH PPROC < 5x5')
        end if
      end if
    end if

    if ( ifrest ) then
      doing_restart = .true.
    end if
    !
    ! Calculate the time step in minutes.
    !
    dtsec = dt
    if ( idynamic == 3 ) then
      dtbat = dtsrf
    else
      dtbat = dt
    end if
    rdt   = d_one/dt
    dtbdys = real(ibdyfrq,rkx)*secph
    !
    ! Reset the options/calculate variables using namelist info
    !
    bdydate1 = idate1

    alarm_hour => rcm_alarm(rcmtimer,3600.0_rkx)
    alarm_day => rcm_alarm(rcmtimer,86400.0_rkx)
    alarm_in_bdy => rcm_alarm(rcmtimer,dtbdys)

    if ( abs(outnwf) > d_zero ) then
      alarm_out_nwf => rcm_alarm(rcmtimer,secpd*abs(outnwf))
    else
      alarm_out_nwf => null( )
    end if
    if ( abs(savfrq) > d_zero ) then
      alarm_out_sav => rcm_alarm(rcmtimer,secpd*abs(savfrq))
    else
      alarm_out_sav => null( )
    end if
    alarm_out_atm => rcm_alarm(rcmtimer,secph*atmfrq)
    alarm_out_rad => rcm_alarm(rcmtimer,secph*radfrq)
    alarm_out_srf => rcm_alarm(rcmtimer,secph*srffrq)
    alarm_out_shf => alarm_hour
    alarm_out_sts => alarm_day
    if ( lakemod == 1 ) then
      alarm_out_lak => rcm_alarm(rcmtimer,secph*lakfrq)
    end if
    if ( ichem == 1 ) then
      alarm_out_che => rcm_alarm(rcmtimer,secph*chemfrq)
    end if
    alarm_out_opt => rcm_alarm(rcmtimer,secph*optfrq)
    if ( nsg > 1 ) then
      alarm_out_sub => rcm_alarm(rcmtimer,secph*subfrq)
    end if

    syncro_rep => rcm_syncro(rcmtimer,3.0_rkx*3600.0_rkx)
    syncro_srf => rcm_syncro(rcmtimer,dtsrf)
    syncro_cum => rcm_syncro(rcmtimer,dtcum)
    syncro_rad => rcm_syncro(rcmtimer,dtrad)
    syncro_emi => rcm_syncro(rcmtimer,dtabem)
    syncro_che => rcm_syncro(rcmtimer,dtche)
    if ( debug_level > 0 ) then
      syncro_dbg => rcm_syncro(rcmtimer,secph*dbgfrq)
    end if
    if ( irrtm == 1 ) then
      syncro_radfor => rcm_syncro(rcmtimer,dtrad*nradfo)
    end if
    if ( iocncpl == 1 .or. iwavcpl == 1 .or. icopcpl == 1 ) then
      syncro_cpl => rcm_syncro(rcmtimer,cpldt)
    end if

    if ( idynamic == 1 ) then
      do ns = 1 , nsplit
        dtsplit(ns) = dt*(d_half/real(nsplit-ns+1,rkx))
        dtau(ns) = dtsplit(ns)
      end do
    end if
    dtsq = dt*dt
    dtcb = dt*dt*dt
    if ( idynamic < 3 ) then
      dt2 = d_two*dt
    else
      dt2 = dt
    end if

    intbdy = rcm_time_interval(ibdyfrq,uhrs)
    intsom = rcm_time_interval(1,umnt)

    if ( myid == italk ) then
      appdat = tochar(idate0)
      write(stdout,*) 'Initial date of the global simulation: ', appdat
      appdat = tochar(idate1)
      write(stdout,*) 'Initial date of this run             : ', appdat
      appdat = tochar(idate2)
      write(stdout,*) 'Final date of this run               : ', appdat
      write(stdout,*) 'Total simulation lenght              : ', hspan, ' hours'
      write(stdout,'(a,f11.6)') ' Timestep in seconds = ', dtsec
      if ( idynamic == 1 ) then
        write(stdout,'(a,2f11.6)') ' Split explicit dtau = ', dtau
      end if
    end if

    call bcast(dirter,256)
    call bcast(dirglob,256)
    call bcast(dirout,256)
    call bcast(domname,64)
    if ( irceideal == 1 ) then
      mddom%ht = 0.0_rkx
      mddom%lndcat = 15.0_rkx
      mddom%lndtex = 14.0_rkx
      mddom%mask = 0.0_rkx
      mddom%msfx = 1.0_rkx
      dl = raddeg * (ds*d_1000)/earthrad
      do i = ide1 , ide2
        do j = jde1 , jde2
          mddom%xlat(j,i) = clat - dl * (real(iy,rkx)*d_half - i + 0.5_rkx)
          mddom%xlon(j,i) = clon - dl * (real(jx,rkx)*d_half - j + 0.5_rkx)
          mddom%coriol(j,i) = eomeg2*sin(mddom%xlat(j,i)*degrad)
        end do
      end do
      if ( idynamic == 3 ) then
        do i = ide1 , ide2
          do j = jde1 , jde2
            mddom%ulat(j,i) = mddom%xlat(j,i)
            mddom%ulon(j,i) = clon - dl * (real(jx,rkx)*d_half - j + 1.0_rkx)
            mddom%vlat(j,i) = clat - dl * (real(iy,rkx)*d_half - i + 1.0_rkx)
            mddom%vlon(j,i) = mddom%xlon(j,i)
          end do
        end do
        mddom%msfu = 1.0_rkx
        mddom%msfv = 1.0_rkx
      else
        do i = ide1 , ide2
          do j = jde1 , jde2
            mddom%dlat(j,i) = clat - dl * (real(iy,rkx)*d_half - i + 1.0_rkx)
            mddom%dlon(j,i) = clon - dl * (real(jx,rkx)*d_half - j + 1.0_rkx)
          end do
        end do
        mddom%msfd = 1.0_rkx
      end if
      if ( idynamic == 2 ) then
        base_state_ts0 = 288.15_rkx
      end if
    else
      call read_domain_info(mddom%ht,mddom%lndcat,mddom%lndtex,            &
                            mddom%mask,mddom%area,                         &
                            mddom%xlat,mddom%xlon,mddom%dlat,mddom%dlon,   &
                            mddom%ulat,mddom%ulon,mddom%vlat,mddom%vlon,   &
                            mddom%msfx,mddom%msfd,mddom%msfu,mddom%msfv,   &
                            mddom%coriol,mddom%snowam,mddom%smoist,        &
                            mddom%rmoist,mddom%dhlake,base_state_ts0)
    end if
    if ( moloch_do_test_1 ) then
      ifrayd = 0
      mddom%ht = 0.0_rkx
      mddom%lndcat = 15.0_rkx
      mddom%lndtex = 14.0_rkx
      mddom%mask = 0.0_rkx
    end if
    if ( moloch_do_test_2 ) then
      !mddom%ht(jde1:jde2,ide1:ide2) = 100.0_rkx * &
      !              abs(sin(mddom%xlat(jde1:jde2,ide1:ide2)*degrad))
      mddom%ht = 0.0_rkx
      mddom%lndcat = 15.0_rkx
      mddom%lndtex = 14.0_rkx
      mddom%mask = 0.0_rkx
      mddom%msfu = d_one
      mddom%msfv = d_one
      mddom%msfx = d_one
      mddom%coriol = d_one
    end if

    if ( idynamic == 3 ) then
      ptop = 0.1_rkx ! assume 1 mbar (.1 cbar)
    end if
    call bcast(ds)
    call bcast(ptop)
    call bcast(xcone)

    dx = ds * d_1000
    dx2 = d_two*dx
    dx4 = d_four*dx
    dx8 = 8.0_rkx*dx
    dx16 = 16.0_rkx*dx
    dxsq = dx*dx
    rdx = 1.0_rkx/dx
    rdxsq = 1.0_rkx/dxsq

    if ( idynamic == 3 ) then
      mo_c1 = sqrt(d_two)*150.0_rkx*dtsec/dx/real(mo_nadv,rkx)
      mo_c2 = sqrt(d_two)*sqrt(cpd/cvd*rgas*313.16_rkx)* &
              dtsec/dx/real(mo_nadv,rkx)/real(mo_nsound,rkx)
      if ( myid == italk ) then
        write(stdout,'(a, f9.3, a, i2, a)') &
           ' Advection timestep = ', dtsec/real(mo_nadv,rkx), &
           ' (factor ', mo_nadv, ')'
        write(stdout,'(a, f9.4)') &
           ' Max. Courant number for horizontal advection = ', mo_c1
        write(stdout,'(a, f9.3, a, i2, a)') ' Sound waves timestep = ', &
           dtsec/real(mo_nadv,rkx)/real(mo_nsound,rkx), &
           ' (factor ', mo_nsound, ')'
        write(stdout,'(a, f9.4)') &
           ' Courant number of horizontal sound waves = ', mo_c2
      end if
    end if

    !
    ! Calculate boundary areas per processor
    !
    call setup_boundaries(cross,cross,ba_cr)
    if ( idynamic == 3 ) then
      call setup_boundaries(dot,cross,ba_ut)
      call setup_boundaries(cross,dot,ba_vt)
    else
      call setup_boundaries(dot,dot,ba_dt)
    end if

    call allocate_v2dbound(xpsb,cross)
    call allocate_v2dbound(xtsb,cross)
    call allocate_v3dbound(xtb,kz,cross)
    call allocate_v3dbound(xqb,kz,cross)
    call allocate_v3dbound(xub,kz,dot)
    call allocate_v3dbound(xvb,kz,dot)
    if ( idynamic == 2 ) then
      call allocate_v3dbound(xppb,kz,cross)
      call allocate_v3dbound(xwwb,kzp1,cross)
    else if ( idynamic == 3 ) then
      call allocate_v3dbound(xpaib,kz,cross)
      call allocate_v3dbound(xwwb,kzp1,cross)
    end if

    if ( myid == italk ) then
      write(stdout,*) 'Setting IPCC scenario to ', scenario
      if ( scenario == 'CONST' ) then
        write(stdout,*) 'Selected value at year ', ghg_year_const
      end if
    end if

    call set_scenario(scenario,rcmtimer%year,rcmtimer%month)

    if ( myid == italk ) then
      if ( ifrest .and. idate0 == idate1 ) then
        write(stderr,*) 'Error in parameter set.'
        write(stderr,*) 'Cannot set idate0 == idate1 on restart run'
        write(stderr,*) 'Correct idate0.'
        call fatal(__FILE__,__LINE__, &
                   'IDATE0==IDATE1 ON RESTART')
      else if ( .not. ifrest .and. idate0 /= idate1 ) then
        write(stderr,*) 'Error in parameter set.'
        write(stderr,*) 'Cannot set idate0 /= idate1 on non restart run'
        write(stderr,*) 'Correct idate1.'
        call fatal(__FILE__,__LINE__, &
                   'IDATE0/=IDATE1 ON NON RESTART')
      end if
    end if

    if ( myid == italk ) then
      write(stdout,*) 'Create SAV files : ' , ifsave
      write(stdout,*) 'Create ATM files : ' , ifatm
      write(stdout,*) 'Create RAD files : ' , ifrad
      write(stdout,*) 'Create SRF files : ' , ifsrf
      write(stdout,*) 'Create STS files : ' , ifsts
      write(stdout,*) 'Create SHF files : ' , ifshf
      if ( nsg > 1 ) write(stdout,*) 'Create SUB files : ' , ifsub
      if ( lakemod == 1 ) write(stdout,*) 'Create LAK files : ' , iflak
      if ( ichem == 1 ) then
        write(stdout,*) 'Create CHE files : ' , ifchem
        write(stdout,*) 'Create OPT files : ' , ifopt
      end if
      if ( .not. associated(alarm_out_nwf) ) then
        write(stdout,'(a,f6.1)') ' Monthly new files are created'
      else
        write(stdout,'(a,f6.1)') ' Frequency in days for new files: ' , outnwf
      end if
      if ( int(savfrq) == 0 ) then
        write(stdout,'(a,f6.1)') ' Monthly SAV files are written'
      else if ( savfrq > d_zero ) then
        write(stdout,'(a,f6.1)') ' Frequency in days to create SAV : ' , savfrq
      else
        write(stdout,'(a,f6.1)') ' Monthly SAV files are written'
        write(stdout,'(a,f6.1)') ' Frequency in days to create SAV : ' , -savfrq
      end if
      write(stdout,'(a,f6.1)') ' Frequency in hours to create ATM : ' , atmfrq
      write(stdout,'(a,f6.1)') ' Frequency in hours to create RAD : ' , radfrq
      write(stdout,'(a,f6.1)') ' Frequency in hours to create SRF : ' , srffrq
      if ( nsg > 1 ) &
        write(stdout,'(a,f6.1)') ' Frequency in hours to create SUB : ' , subfrq
      if ( lakemod == 1 ) &
        write(stdout,'(a,f6.1)') ' Frequency in hours to create LAK : ' , lakfrq
      if ( ichem == 1 ) then
        write(stdout,'(a,f6.1)') &
          ' Frequency in hours to create CHE : ' , chemfrq
        write(stdout,'(a,f6.1)') &
          ' Frequency in hours to create OPT : ' , optfrq
      end if

      write(stdout,*) 'Physical Parameterizations'
      write(stdout,'(a,i2)') '  Lateral Boundary conditions : ' , iboudy
      write(stdout,'(a,i2)') '  Semi-Lagrangian Advection   : ' , isladvec
      if ( isladvec == 1 ) then
        write(stdout,'(a,i2)') '  QMSL algorithm used         : ' , iqmsl
      end if
      if ( any(icup == -1) ) then
        icup(:) = -1
        write(stdout,'(a)') '  Shallow cumulus scheme '
      else
        write(stdout,'(a,i2)') '  Land cumulus conv. scheme   : ' , icup_lnd
        write(stdout,'(a,i2)') '  Ocean cumulus conv. scheme  : ' , icup_ocn
      end if
      write(stdout,'(a,i2)') '  Moisture schem              : ' , ipptls
      write(stdout,'(a,i2)') '  Ocean Flux scheme           : ' , iocnflx
      if ( iocnflx == 2 ) then
        write(stdout,'(a,i2)') '  Zeng roughness formula      : ' , iocnrough
        write(stdout,'(a,i2)') '  Zeng roughness method       : ' , iocnzoq
      end if
      write(stdout,'(a,i2)') '  Coupling with ocean         : ' , iocncpl
      write(stdout,'(a,i2)') '  Coupling with wave          : ' , iwavcpl
      write(stdout,'(a,i2)') '  Coupling with COP           : ' , icopcpl
      write(stdout,'(a,i2)') '  Pressure gradient force     : ' , ipgf
      write(stdout,'(a,i2)') '  Prescribed LW emissivity    : ' , iemiss
#ifndef CLM
      write(stdout,'(a,i2)') '  Lake model in BATS          : ' , lakemod
      write(stdout,'(a,i2)') '  Simulate diurnal sst cycle  : ' , idcsst
      write(stdout,'(a,i2)') '  Simulate sea ice cover      : ' , iseaice
      write(stdout,'(a,i2)') '  Simulate desert seasons     : ' , idesseas
#endif
      write(stdout,'(a,i2)') '  Enable chem/aerosol model   : ' , ichem
      write(stdout,'(a,i2)') '  Large scale LWP as convect. : ' , iconvlwp
      write(stdout,'(a,i2)') '  Cloud fraction scheme       : ' , icldfrac
      write(stdout,'(a,i2)') '  Marine stratocumulus        : ' , icldmstrat
      write(stdout,'(a,i2)') '  Climate O3 dataset          : ' , iclimao3
      write(stdout,'(a,i2)') '  Climate Aerosol dataset     : ' , iclimaaer
      write(stdout,*) 'Boundary Pameterizations'
      write(stdout,'(a,i2)') '  Num. of bndy points cross  : ', nspgx
      write(stdout,'(a,i2)') '  Num. of bndy points dot    : ', nspgd
      write(stdout,'(a,f9.6)') '  Nudge value high range     : ', high_nudge
      write(stdout,'(a,f9.6)') '  Nudge value medium range   : ', medium_nudge
      write(stdout,'(a,f9.6)') '  Nudge value low range      : ', low_nudge
      write(stdout,'(a,f9.6)') '  Nm paramter                : ', bdy_nm
      write(stdout,'(a,f9.6)') '  Dm paramter                : ', bdy_dm
#ifdef CLM
      write(stdout,*) 'CLM Pameterizations'
      write(stdout,'(a,i2)' ) '  CLM imask                       : ' , imask
      write(stdout,'(a,f9.6)') '  Frequency in hours to write CLM : ', clmfrq
#endif
      write(stdout,*) 'Model Timestep Pameterizations'
      write(stdout,'(a,f12.6)') '  time step for dynamical '// &
            'model in seconds : ' , dt
      write(stdout,'(a,f12.6)') '  time step for surface   '// &
            'model in seconds : ' , dtsrf
      write(stdout,'(a,f12.6)') '  time step for cumulus   '// &
            'model in seconds : ' , dtcum
      write(stdout,'(a,f12.6)') '  time step for radiation '// &
            'model in seconds : ' , dtrad
      write(stdout,'(a,f12.6)') '  time step for emission  '// &
            'model in seconds : ' , dtabem
      if ( ichem == 1 ) then
        write(stdout,'(a,f12.6)') '  time step for chemistry '// &
              'model in seconds : ' , dtche
      end if
    end if

    if ( nsg > 1 ) then
      if ( irceideal == 1 ) then
        call fatal(__FILE__,__LINE__, &
          'The idealized cases can run only with nsg == 1')
      end if
      call read_subdomain_info(mdsub%ht,mdsub%lndcat,mdsub%lndtex,mdsub%mask, &
               mdsub%area,mdsub%xlat,mdsub%xlon,mdsub%dhlake)
      mdsub%ht = mdsub%ht*egrav
    else
      do i = ici1 , ici2
        do j = jci1 , jci2
          mdsub%ht(1,j,i) = mddom%ht(j,i)*egrav
          mdsub%lndcat(1,j,i) = mddom%lndcat(j,i)
          mdsub%lndtex(1,j,i) = mddom%lndtex(j,i)
          mdsub%xlat(1,j,i) = mddom%xlat(j,i)
          mdsub%xlon(1,j,i) = mddom%xlon(j,i)
          mdsub%mask(1,j,i) = mddom%mask(j,i)
          mdsub%area(1,j,i) = mddom%area(j,i)
        end do
      end do
      if ( lakemod == 1 ) then
        do i = ici1 , ici2
          do j = jci1 , jci2
            mdsub%dhlake(1,j,i) = mddom%dhlake(j,i)
          end do
        end do
      end if
    end if
    !
    ! Compute diffusion halo for idynamic /= 3
    !
    if ( idynamic /= 3 ) then
      if ( idiffu == 1 ) then
        idif = 2
      else if ( idiffu == 2 ) then
        idif = 1
      else if ( idiffu == 3 ) then
        idif = 3
      else
        call fatal(__FILE__,__LINE__, &
          'INCORRECT DIFFUSION SCHEME! CHECK IDIFFU IN NAMELIST')
      end if
    end if
    !
    !------invert mapscale factors and convert hgt to geopotential
    !
    do i = ide1 , ide2
      do j = jde1 , jde2
        mddom%ht(j,i)   = mddom%ht(j,i)*egrav
      end do
    end do
    if ( idynamic /= 3 ) then
      do i = ide1 , ide2
        do j = jde1 , jde2
          mddom%msfd(j,i) = d_one/mddom%msfd(j,i)
          mddom%msfx(j,i) = d_one/mddom%msfx(j,i)
        end do
      end do
      do i = idi1 , idi2
        do j = jdi1 , jdi2
          mddom%dmsf(j,i) = d_one/(mddom%msfd(j,i)*mddom%msfd(j,i)*dx16)
          mddom%xmsf(j,i) = d_one/(mddom%msfx(j,i)*mddom%msfx(j,i)*dx4)
        end do
      end do
      call exchange(mddom%msfx,idif,jde1,jde2,ide1,ide2)
      call exchange(mddom%msfd,idif,jde1,jde2,ide1,ide2)
    end if
    !
    !-----compute half sigma levels.
    !
    if ( idynamic == 3 ) then
      call model_zitaf(zita)
      call model_zitah(zitah)
      mo_dzita = zita(kz)
      sigma = sigmazita(zita)
      hsigma = sigmazita(zitah)
      fak = md_ak(zita)
      fbk = md_bk(zita)
      ak = md_ak(zitah)
      bk = md_bk(zitah)
      do k = 1 , kz
        dsigma(k) = (sigma(k+1) - sigma(k))
      end do
    else
      do k = 1 , kz
        hsigma(k) = (sigma(k+1) + sigma(k))*d_half
        dsigma(k) = (sigma(k+1) - sigma(k))
      end do
    end if

    call exchange(mddom%xlat,1,jde1,jde2,ide1,ide2)
    call exchange(mddom%xlon,1,jde1,jde2,ide1,ide2)
    call exchange(mddom%ht,2,jde1,jde2,ide1,ide2)
    !
    !-----compute land/water mask
    !
    do i = ici1 , ici2
      do j = jci1 , jci2
        if ( mddom%mask(j,i) > 0.1_rkx ) then
          mddom%ldmsk(j,i) = 1
        else
          mddom%ldmsk(j,i) = 0
        end if
      end do
    end do
    !
    !-----compute land/water mask on subgrid space
    !
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          if ( mdsub%mask(n,j,i) > 0.1_rkx ) then
            mdsub%ldmsk(n,j,i) = 1
          else
            mdsub%ldmsk(n,j,i) = 0
          end if
        end do
      end do
    end do
    !
    ! Initializations
    !
    if ( idynamic /= 3 ) then
      call init_advection
      if ( isladvec == 1 ) then
        call init_sladvection
      end if
    end if
    call init_slice
    call init_micro
    if ( ichem == 1 ) then
      call init_chem
      do n = 1 , ntr
        call bcast(chtrname(n),6)
      end do
    end if
    call init_surface_model
    call init_cumulus
    call init_radiation
    if ( islab_ocean == 1 ) then
      call allocate_mod_slabocean
      call init_slabocean(sfs,mddom%lndcat,fsw,flw,mddom%xlon,mddom%xlat)
    end if
    !
    ! Setup Boundary condition routines.
    !
    call setup_bdycon
    if ( ichem == 1 ) call setup_che_bdycon

    if ( idynamic == 2 ) then
      call make_reference_atmosphere
      call compute_full_coriolis_coefficients
    else if ( idynamic == 3 ) then
      call compute_moloch_static
    end if

    if ( iboudy < 0 .or. iboudy > 6 ) then
      call fatal(__FILE__,__LINE__, &
                 'UNSUPPORTED BDY SCHEME.')
    end if

    if ( myid == italk ) then
      write(stdout,*) 'Domain grid parameters:'
      write(stdout,'(a,a)') '  Map Projection        : ',iproj
      write(stdout,'(a,i4,a,i4,a,i3)') &
        '  Dot Grid Full Extent  : ',jx,'x',iy,'x',kz
      if ( idynamic /= 3 ) then
        write(stdout,'(a,f11.6,a)') '  Model Top Pressure    : ',ptop,' cb'
      end if
      write(stdout,'(a,f11.6,a)') '  Model Grid Spacing    : ',ds,' km'
      write(stdout,'(a,f11.6,a)') '  Proj Center Latitude  : ',clat,' deg'
      write(stdout,'(a,f11.6,a)') '  Proj Center longitude : ',clon,' deg'
      if ( iproj == 'ROTMER' .or. iproj == 'ROTLLR' ) then
        write(stdout,'(a,f11.6,a)') '  Pole Latitude         : ',plat,' deg'
        write(stdout,'(a,f11.6,a)') '  Pole longitude        : ',plon,' deg'
      else if ( iproj == 'LAMCON' ) then
        write(stdout,'(a,f11.6,a)') '  True Latitude 1       : ',truelatl,' deg'
        write(stdout,'(a,f11.6,a)') '  True Latitude 2       : ',truelath,' deg'
      end if
    end if

    if ( ipptls > 0 ) then
      cevaplnd = max(cevaplnd,d_zero)
      caccrlnd = max(caccrlnd,d_zero)
      cevapoce = max(cevapoce,d_zero)
      caccroce = max(caccroce,d_zero)
    end if

    if ( myid == italk ) then
      if ( idynamic /= 3 ) then
        if ( idiffu == 1 ) then
          write(stdout,*) 'Fourth order / second order boundary diffusion'
        else if ( idiffu == 2 ) then
          write(stdout,*) 'Fourth order diffusion LeVeque Laplacian'
        else if ( idiffu == 3 ) then
          write(stdout,*) 'Sixth order diffusion with flux limiter'
        end if
      end if
      if ( ibltyp == 1 ) then
        write(stdout,*) 'Holtslag PBL Scheme'
        write(stdout,'(a,f11.6)') '  ricr_ocn     = ', ricr_ocn
        write(stdout,'(a,f11.6)') '  ricr_lnd     = ', ricr_lnd
        write(stdout,'(a,f11.6)') '  zhnew_fac    = ', zhnew_fac
        write(stdout,'(a,i2)') '  ifaholt      = ', ifaholt
        write(stdout,'(a,i2)') '  ifaholtth10  = ', ifaholtth10
        write(stdout,'(a,i4)') '  holtth10iter = ', holtth10iter
      else if ( ibltyp == 2 ) then
        write(stdout,*) 'UW TCM Parameters'
        write(stdout,'(a,f11.6)') '  rstbl     = ', rstbl
        write(stdout,'(a,f11.6)') '  atwo      = ', atwo
        write(stdout,'(a,f11.6)') '  czero     = ', czero
        write(stdout,'(a,f11.6)') '  nuk       = ', nuk
        write(stdout,'(a,i3)')    '  iuwvadv   = ', iuwvadv
      else if ( ibltyp == 3 ) then
        write(stdout,*) 'GFS PBL Scheme'
      else if ( ibltyp == 4 ) then
        write(stdout,*) 'MYJ PBL Scheme'
      else
        write(stdout,*) &
          'Model frictionless and insulated for the lower boundary.'
      end if
      write(stdout,*) 'Cloud fraction schemes parameters'
      write(stdout,'(a,i2)' )   '  # of bottom no cloud model levels : ',ncld
      write(stdout,'(a,f11.6)') '  Maximum relative humidity         : ' ,rhmax
      write(stdout,'(a,f11.6)') '  Minimum relative humidity         : ' ,rhmin
      if ( icldfrac == 0 ) then
        write(stdout,'(a,f11.6)') '  RH0 temperature threshold         : ' ,tc0
        write(stdout,'(a,f11.6,a,f11.6)')                      &
            '  Relative humidity thresholds: Land = ',rh0land, &
            ' Ocean = ',rh0oce
      end if
      write(stdout,'(a,f11.6)') &
          '  Maximum total cloud cover for rad : ', cftotmax
      write(stdout,'(a,f11.6)') &
          '  Conv. CF UMF factor deep cumulus  : ', kfac_deep
      write(stdout,'(a,f11.6)') &
          '  Conv. CF UMF factor deep cumulus  : ', kfac_deep
      write(stdout,'(a,f11.6)') &
          '  Conv. CF UMF factor k2            : ', k2_const
      write(stdout,'(a,l11)') &
          '  Surface radiation hack            : ', lsrfhack
      write(stdout,'(a,l11)') &
          '  Vavrus-Waliser Arctic cloud fix   : ', larcticcorr
      if ( ipptls == 1 ) then
        !
        ! specify the constants used in the model.
        !     conf  : condensation efficiency
        !     qcth  : threshold for the onset of autoconversion
        !     qck1  : constant autoconversion rate
        ! all the other constants are used to compute the cloud
        ! microphysical parameterization (ref. orville & kopp, 1977 jas).
        !
        write(stdout,*) 'SUBEX large scale precipitation parameters'
        write(stdout,'(a,f11.6,a,f11.6)')                      &
            '  Auto-conversion rate:         Land = ',qck1land,&
            ' Ocean = ',qck1oce
        write(stdout,'(a,f11.6,a,f11.6)')                      &
            '  Gultepe factors:              Land = ',gulland ,&
            ' Ocean = ',guloce
        if ( cevaplnd <= d_zero ) then
          write (stdout,*) '  Land raindrop evaporation not included'
        else
          write(stdout,'(a,f11.6)') &
            '  Land Raindrop Evaporation Rate    : ',cevaplnd
        end if
        if ( cevapoce <= d_zero ) then
          write (stdout,*) '  Ocean raindrop evaporation not included'
        else
          write(stdout,'(a,f11.6)') &
            '  Ocean Raindrop Evaporation Rate   : ', cevapoce
        end if
        if ( caccrlnd <= d_zero ) then
          write(stdout, *) '  Land raindrop accretion not included'
        else
          write(stdout,'(a,f11.6)') &
            '  Land Raindrop Accretion Rate      : ', caccrlnd
        end if
        if ( caccroce <= d_zero ) then
          write(stdout, *) '  Ocean raindrop accretion not included'
        else
          write(stdout,'(a,f11.6)') &
            '  Ocean Raindrop Accretion Rate     : ', caccroce
        end if
        write(stdout,'(a,f11.6)') &
            '  Condensation efficiency           : ', conf
      end if
    end if

    call allocate_cumulus
    !
    !-----compute the vertical interpolation coefficients for t and qv.
    !
    twt(1,1) = d_zero
    twt(1,2) = d_zero
    qcon(1) = d_zero
    do k = 2 , kz
      twt(k,1) = (sigma(k)-hsigma(k-1))/(hsigma(k)-hsigma(k-1))
      twt(k,2) = d_one - twt(k,1)
      qcon(k) = (sigma(k)-hsigma(k))/(hsigma(k-1)-hsigma(k))
    end do

    if ( any(icup == 1) ) then
      !
      ! specify heating profile (twght) and weighted function
      ! for moisture fluxes due to convection (vqflx)
      ! assume base of cloud varies as  < kbase = 5,kz >
      ! top  of cloud varies as  < ktop  = 1,kbase-3 >
      ! exceptions to this are treated explicitly in subroutine
      ! "cupara".
      !
      do kbase = 5 , kz
        do ktop = 1 , kbase - 3
          do k = 1 , kz
            twght(k,kbase,ktop) = d_zero
            vqflx(k,kbase,ktop) = d_zero
          end do
          !
          ! get twght from 1/2 level sigma values
          !
          bb = log(hsigma(ktop)) + log(hsigma(kbase))
          cc = log(hsigma(ktop))*log(hsigma(kbase))
          ssum = d_zero
          do k = ktop , kbase
            xx = log(hsigma(k))
            twght(k,kbase,ktop) = (xx*xx) - (bb*xx) + cc
            ssum = ssum + twght(k,kbase,ktop)*dsigma(k)
          end do
          do k = ktop , kbase
            twght(k,kbase,ktop) = twght(k,kbase,ktop)/ssum
          end do
          !
          ! get vqflx from  d(w*q) / dsigma on full levels
          ! do computations in p to avoid sigma=0. discontinuity
          !
          xtop = log((d_100-ptop)*sigma(ktop)+ptop)
          xbot = log((d_100-ptop)*sigma(kbase+1)+ptop)
          bb = xtop + xbot
          cc = xtop*xbot
          vqmax = d_zero
          ssum = d_zero
          xx = xtop
          yy = xbot
          wk = (xx*xx) - (bb*xx) + cc
          qk = -((yy*yy)-(bb*yy)+cc)
          do k = ktop , kbase
            xx = log((d_100-ptop)*sigma(k+1)+ptop)
            yy = log((d_100-ptop) * &
                 (sigma(ktop)+sigma(kbase+1)-sigma(k+1))+ptop)
            wkp1 = (xx*xx) - (bb*xx) + cc
            qkp1 = -((yy*yy)-(bb*yy)+cc)
            vqflx(k,kbase,ktop) = -((wkp1*qkp1)-(wk*qk))/dsigma(k)
            ssum = ssum + vqflx(k,kbase,ktop)
            if ( abs(vqflx(k,kbase,ktop)) > vqmax ) &
                 vqmax = abs(vqflx(k,kbase,ktop))
            wk = wkp1
            qk = qkp1
          end do
          do k = ktop , kbase
            vqflx(k,kbase,ktop) = vqflx(k,kbase,ktop)*vqrang/vqmax
          end do
        end do
      end do
      if ( myid == italk ) then
        write(stdout, *) 'Anthes-Kuo Convection Scheme used.'
      end if
    end if
    if ( any(icup == 2) ) then
      kbmax = kz
      do k = 1 , kz - 1
        if ( hsigma(k) <= skbmax ) kbmax = kz - k
      end do
      if ( myid == italk ) then
        write(stdout,*) 'Grell Convection Scheme used.'
        write(stdout,'(a,f11.6)') '  Max Shear       : ' , shrmax
        write(stdout,'(a,f11.6)') '  Min Shear       : ' , shrmin
        write(stdout,'(a,f11.6)') '  Max PPT eff     : ' , edtmax
        write(stdout,'(a,f11.6)') '  Min PPT eff     : ' , edtmin
        write(stdout,'(a,f11.6)') '  Max PPT eff(o)  : ' , edtmaxo
        write(stdout,'(a,f11.6)') '  Min PPT eff(o)  : ' , edtmino
        write(stdout,'(a,f11.6)') '  Max PPT eff(x)  : ' , edtmaxx
        write(stdout,'(a,f11.6)') '  Min PPT eff(x)  : ' , edtminx
        write(stdout,'(a,f11.6)') '  Max PBC         : ' , pbcmax
        write(stdout,'(a,f11.6)') '  Min Cloud Depth : ' , mincld
        write(stdout,'(a,f11.6)') '  Max Cloud Base  : ' , skbmax
        write(stdout,'(a,f11.6)') '  Max Heating     : ' , htmax
        write(stdout,'(a,f11.6)') '  Min Heating     : ' , htmin
        if ( igcc == 1 ) then
          write(stdout,*) ' Arakawa-Schubert (1974) Closure Assumption'
        else if ( igcc == 2 ) then
          write(stdout,*) ' Fritsch-Chappell (1980) Closure Assumption'
          write(stdout,'(a,f11.6)') '  ABE removal timescale : ',dtauc
        else
          write(stderr,*) 'Unknown Closure Assumption for Grell.'
          call fatal(__FILE__,__LINE__, &
                    'PARAMETER ERROR')
        end if
      end if

      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( isocean(mddom%lndcat(j,i)) ) then
            shrmax2d(j,i) = shrmax_ocn
            shrmin2d(j,i) = shrmin_ocn
            edtmax2d(j,i) = edtmax_ocn
            edtmin2d(j,i) = edtmin_ocn
            edtmaxo2d(j,i) = edtmaxo_ocn
            edtmino2d(j,i) = edtmino_ocn
            edtmaxx2d(j,i) = edtmaxx_ocn
            edtminx2d(j,i) = edtminx_ocn
          else
            shrmax2d(j,i) = shrmax
            shrmin2d(j,i) = shrmin
            edtmax2d(j,i) = edtmax
            edtmin2d(j,i) = edtmin
            edtmaxo2d(j,i) = edtmaxo
            edtmino2d(j,i) = edtmino
            edtmaxx2d(j,i) = edtmaxx
            edtminx2d(j,i) = edtminx
          end if
          pbcmax2d(j,i) = pbcmax
          mincld2d(j,i) = mincld
          kbmax2d(j,i) = kbmax
          htmax2d(j,i) = htmax
          htmin2d(j,i) = htmin
          dtauc2d(j,i) = dtauc*secpm
        end do
      end do
    end if
    if ( any(icup == 3) ) then
      if ( myid == italk ) then
        write(stderr,*) &
          'WARNING : The Betts-Miller Convection scheme is not ', &
          'properly implemented'
      end if
    end if
    if ( any(icup == 4) ) then
      if ( myid == italk ) then
        write(stdout,*) 'Emanuel (1991) Convection Scheme V4.3C used.'
        write(stdout,'(a,i3)')    '  Min Convection origin             : ', &
          minorig
        write(stdout,'(a,f11.6)') '  Autoconversion Threshold (ocean)  : ', &
          elcrit_ocn
        write(stdout,'(a,f11.6)') '  Autoconversion Threshold (land)   : ', &
          elcrit_lnd
        write(stdout,'(a,f11.6)') &
          '  Autoconversion Threshold to zero  : ',tlcrit
        write(stdout,'(a,f11.6)') '  Entrainment Coefficient           : ',entp
        write(stdout,'(a,f11.6)') '  Fractional area of uns. downdraft : ',sigd
        write(stdout,'(a,f11.6)') '  Fractional area of uns. downdraft : ',sigd
        write(stdout,'(a,f11.6)') &
          '  Fall speed of rain                : ',omtrain
        write(stdout,'(a,f11.6)') &
          '  Fall speed of snow                : ',omtsnow
        write(stdout,'(a,f11.6)') &
          '  Rain evaporation coefficient      : ',coeffr
        write(stdout,'(a,f11.6)') &
          '  Snow evaporation coefficient      : ',coeffs
        write(stdout,'(a,f11.6)') '  Convective momentum transport coef: ',cu
        write(stdout,'(a,f11.6)') '  Downdraft velocity scale          : ',betae
        write(stdout,'(a,f11.6)') '  Max negative perturbation blw LFC : ',dtmax
        write(stdout,'(a,f11.6)') &
          '  Quasi-equilibrium approach rate 1 : ',alphae
        write(stdout,'(a,f11.6)') '  Quasi-equilibrium approach rate 2 : ',damp
        write(stdout,'(a,f11.6)') '  Precipitation efficienct (land)   : ', &
          epmax_lnd
        write(stdout,'(a,f11.6)') '  Precipitation efficienct (ocean)  : ', &
          epmax_ocn
      end if
    end if
    if ( any(icup == 5)  ) then
      if ( myid == italk ) then
        write(stdout,*) 'Tiedtke (1986) Convection Scheme used.'
        write(stdout,'(a,i2)')    &
          '  Used Scheme                       : ',iconv
        write(stdout,'(a,f11.6)') &
          '  Entrainment rate pen. conv. land  : ',entrpen_lnd
        write(stdout,'(a,f11.6)') &
          '  Entrainment rate pen. conv. ocean : ',entrpen_ocn
        write(stdout,'(a,f11.6)') &
          '  Entrainment rate shallow conv.    : ',entrscv
        write(stdout,'(a,f11.6)') &
          '  Entrainment rate midlev conv.     : ',entrmid
        write(stdout,'(a,f11.6)') &
          '  Entrainment rate cumulus downdraft: ',entrdd
        write(stdout,'(a,f11.6)') &
          '  CLW to rain conversion factor     : ',cprcon
        write(stdout,'(a,f11.6)') &
          '  Max entrainment                   : ',entrmax
        write(stdout,'(a,f12.6)') &
          '  CAPE adjustment timescale         : ',cmtcape
      end if
    end if
    if ( any(icup == 6)  ) then
      if ( myid == italk ) then
        write(stdout,*) 'Kain Fritsch scheme used.'
        write(stdout,'(a,f11.6)') &
          '  Vertical velocity threshold : ', kf_wthreshold
        write(stdout,'(a,f11.6)') &
          '  Entrainment rate            : ', kf_entrate
        write(stdout,'(a,f11.6)') &
          '  Conversion rate             : ', kf_convrate
        write(stdout,'(a,f11.6)') &
          '  Maximum prec. efficiency    : ', kf_max_pef
        write(stdout,'(a,f11.6)') &
          '  Minimum prec. efficiency    : ', kf_min_pef
        write(stdout,'(a,f11.6)') &
          '  Updraft start elevation     : ', kf_dpp
        write(stdout,'(a,f11.6)') &
          '  CAPE removal time min.      : ', kf_min_dtcape
        write(stdout,'(a,f11.6)') &
          '  CAPE removal time max.      : ', kf_max_dtcape
        write(stdout,'(a,f11.6)') &
          '  TKE maximum in sub cloud    : ', kf_tkemax
      end if
    end if

    call init_pblscheme
    !
    ! Convective Cloud Cover
    !
    afracl = 0.25_rkx  ! frac. cover for conv. precip. when dx=dxlarg
    afracs = clfrcvmax !   "     "    "    "      "     "   dx=dxsmal
    dlargc = 100.0_rkx
    dsmalc = 10.0_rkx
    dxtemc = min(max(ds,dsmalc),dlargc)
    clfrcv = afracl + (afracs-afracl)*((dlargc-dxtemc)/(dlargc-dsmalc))**2
    clfrcv = min(clfrcv,afracs)
    clfrcv = max(clfrcv,afracl)
    if ( myid == italk ) then
      write(stdout,*) &
        'Convective Cloud Cover parameters after resolution scaling'
      write(stdout,'(a,f11.6)') '  Maximum Convective Cloud Cover : ',clfrcv
      write(stdout,'(a,f11.6)') '  Convective Cloud Water         : ',cllwcv
    end if

    do i = ici1 , ici2
      do j = jci1 , jci2
        if ( isocean(mddom%lndcat(j,i)) ) then
          rh0(j,i)   = rh0oce
        else
          rh0(j,i)   = rh0land
        end if
      end do
    end do

    if ( ipptls == 1 ) then
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( isocean(mddom%lndcat(j,i)) ) then
            qck1(j,i)  = qck1oce  ! OCEAN
            cgul(j,i)  = guloce
            cevap(j,i) = cevapoce
            caccr(j,i) = caccroce
          else
            qck1(j,i)  = qck1land ! LAND
            cgul(j,i)  = gulland
            cevap(j,i) = cevaplnd
            caccr(j,i) = caccrlnd
          end if
        end do
      end do
    end if
    !
    ! Setup Holtslag PBL Critical Richardson Number
    !
    if ( ibltyp == 1 ) then
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( isocean(mddom%lndcat(j,i)) ) then
            ricr(j,i) = ricr_lnd
          else
            ricr(j,i) = ricr_ocn
          end if
        end do
      end do
    end if
    if ( any(icup == 1) ) then
      sig700 = (70.0_rkx-ptop)/(d_100-ptop)
      do k = 1 , kz
        k700 = k
        if ( sig700 <= sigma(k+1) .and. sig700 > sigma(k) ) exit
      end do
    end if
    if ( myid == italk ) then
      write(stdout,*) &
        'The surface energy budget is used to calculate the ground temperature.'
      write(stdout,'(a,f5.0,a)') &
        ' The radiation is computed every ',dtrad/60.0_rkx,' minutes.'

      if ( iboudy == 0 ) then
        write(stdout,*) 'The lateral boundary conditions are fixed.'
      else if ( iboudy == 1 ) then
        write(stdout,*) 'Relaxation boundary conditions (linear method)'
      else if ( iboudy == 2 ) then
        write(stdout,*) 'Time dependent boundary conditions are used.'
      else if ( iboudy == 3 ) then
        write(stdout,*) 'Inflow/outflow boundary conditions are used.'
      else if ( iboudy == 4 ) then
        write(stdout,*) 'Sponge boundary conditions are used.'
      else if ( iboudy == 5 ) then
        write(stdout,*) 'Relaxation boundary conditions (exponential method)'
      else if ( iboudy == 6 ) then
        write(stdout,*) 'Relaxation boundary conditions (sinusoidal method)'
      end if

      if ( idynamic /= 3 ) then
        write(stdout,'(a,7x,a,11x,a,6x,a,7x,a,7x,a,9x,a)') '# k','sigma','a',&
          'dsigma','twt(1)','twt(2)','qcon'

        do k = 1 , kz
          write(stdout, &
            '(1x,i2,5x,f7.4,5x,f7.4,5x,f7.4,5x,f8.4,5x,f8.4,5x,f8.4)') &
            k , sigma(k) , hsigma(k) , dsigma(k) , twt(k,1) , twt(k,2) , qcon(k)
        end do
        write(stdout,'(1x,i2,5x,f7.4)') kzp1 , sigma(kzp1)
      end if
    end if

    if ( itweak == 1 ) then
      if ( myid == italk ) then
        write(stdout,*) 'TWEAKING OF DATA ENABLED!'
        write(stdout,*) 'THIS RUN IS TO BE CONSIDERED A NON STANDARD SCENARIO!'
        if ( itweak_temperature == 1 ) then
          write(stdout,'(a,f11.6,a)') ' Value added to temperature      : ', &
                  temperature_tweak , ' K'
        end if
        if ( itweak_solar_irradiance == 1 ) then
          write(stdout,'(a,f11.6,a)') ' Value added to solar irradiance : ', &
                  solar_tweak , ' W m-2'
        end if
        if ( itweak_greenhouse_gases == 1 ) then
          write(stdout,'(a,f11.6)') ' CO2 concentration factor        : ', &
                  gas_tweak_factors(1)
          write(stdout,'(a,f11.6)') ' CH4 concentration factor        : ', &
                  gas_tweak_factors(2)
          write(stdout,'(a,f11.6)') ' N2O concentration factor        : ', &
                  gas_tweak_factors(3)
          write(stdout,'(a,f11.6)') ' CFC11 concentration factor      : ', &
                  gas_tweak_factors(4)
          write(stdout,'(a,f11.6)') ' CFC12 concentration factor      : ', &
                  gas_tweak_factors(5)
        end if
      end if
    end if

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif

    return

100 call fatal(__FILE__,__LINE__, 'Error reading RESTARTPARAM')
101 call fatal(__FILE__,__LINE__, 'Error reading TIMEPARAM')
102 call fatal(__FILE__,__LINE__, 'Error reading OUTPARAM')
103 call fatal(__FILE__,__LINE__, 'Error reading PHYSICSPARAM')
104 call fatal(__FILE__,__LINE__, 'Error reading DYNPARAM')
105 call fatal(__FILE__,__LINE__, 'Error reading NONHYDROPARAM')
106 call fatal(__FILE__,__LINE__, 'Error reading HYDROPARAM')
107 call fatal(__FILE__,__LINE__, 'Error reading CLDPARAM')
108 call fatal(__FILE__,__LINE__, 'Error reading SUBEXPARAM')
109 call fatal(__FILE__,__LINE__, 'Error reading MICROPARAM')
110 call fatal(__FILE__,__LINE__, 'Error reading GRELLPARAM')
111 call fatal(__FILE__,__LINE__, 'Error reading EMANPARAM')
112 call fatal(__FILE__,__LINE__, 'Error reading TIEDTKEPARAM')
113 call fatal(__FILE__,__LINE__, 'Error reading KFPARAM')
114 call fatal(__FILE__,__LINE__, 'Error reading HOLTSLAGPARAM')
115 call fatal(__FILE__,__LINE__, 'Error reading UWPARAM')
116 call fatal(__FILE__,__LINE__, 'Error reading RRTMPARAM')
117 call fatal(__FILE__,__LINE__, 'Error reading SLABOCPARAM')
118 call fatal(__FILE__,__LINE__, 'Error reading CHEMPARAM')
#ifdef CLM
119 call fatal(__FILE__,__LINE__, 'Error reading CLMPARAM')
#endif
120 call fatal(__FILE__,__LINE__, 'Error reading CPLPARAM')
121 call fatal(__FILE__,__LINE__, 'Error reading TWEAKPARAM')

    contains

      subroutine make_reference_atmosphere
        use mod_nhinterp
        implicit none
        integer(ik4) :: i , j , k
        real(rkx) :: ztop
        call nhsetup(ptop,base_state_pressure,logp_lrate,base_state_ts0)
        mddom%ht = mddom%ht * regrav
        call nhbase(ice1,ice2,jce1,jce2,kz,hsigma,mddom%ht, &
                    atm0%ps,atm0%pr,atm0%t,atm0%rho,atm0%z)
        call nhbase(ice1,ice2,jce1,jce2,kzp1,sigma,mddom%ht, &
                    atm0%ps,atm0%pf,atm0%tf,atm0%rhof,atm0%zf)
        mddom%ht = mddom%ht * egrav
        call exchange(atm0%ps,1,jce1,jce2,ice1,ice2)
        call exchange(atm0%pr,1,jce1,jce2,ice1,ice2,1,kz)
        call exchange(atm0%t,1,jce1,jce2,ice1,ice2,1,kz)
        call exchange(atm0%z,1,jce1,jce2,ice1,ice2,1,kz)
        call exchange(atm0%zf,1,jce1,jce2,ice1,ice2,1,kzp1)
        call exchange(atm0%rho,1,jce1,jce2,ice1,ice2,1,kz)
        call psc2psd(atm0%ps,atm0%psdot)
        call exchange(atm0%psdot,1,jde1,jde2,ide1,ide2)
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              atm0%dzf(j,i,k) = atm0%zf(j,i,k) - atm0%zf(j,i,k+1)
            end do
          end do
        end do
        do k = 1 , kz
          do i = idi1 , idi2
            do j = jdi1 , jdi2
              atm0%zd(j,i,k) = d_rfour * &
                   (atm0%z(j,i,k) + atm0%z(j-1,i,k) + &
                    atm0%z(j,i-1,k) + atm0%z(j-1,i-1,k))
            end do
          end do
        end do
        do i = ice1 , ice2
          do j = jci1 , jci2
            dpsdxm(j,i) = (atm0%ps(j+1,i) - atm0%ps(j-1,i)) / &
                          (atm0%ps(j,i)*dx8*mddom%msfx(j,i))
          end do
        end do
        if ( ma%has_bdyleft ) then
          do i = ice1 , ice2
            dpsdxm(jce1,i) = (atm0%ps(jci1,i) - atm0%ps(jce1,i)) / &
                        (atm0%ps(jce1,i)*dx8*mddom%msfx(jce1,i))
          end do
        end if
        if ( ma%has_bdyright ) then
          do i = ice1 , ice2
            dpsdxm(jce2,i) = (atm0%ps(jce2,i) - atm0%ps(jci2,i)) / &
                        (atm0%ps(jce2,i)*dx8*mddom%msfx(jce2,i))
          end do
        end if
        do i = ici1 , ici2
          do j = jce1 , jce2
            dpsdym(j,i) = (atm0%ps(j,i+1) - atm0%ps(j,i-1)) / &
                          (atm0%ps(j,i)*dx8*mddom%msfx(j,i))
          end do
        end do
        if ( ma%has_bdybottom ) then
          do j = jce1 , jce2
            dpsdym(j,ice1) = (atm0%ps(j,ici1) - atm0%ps(j,ice1)) / &
                        (atm0%ps(j,ice1)*dx8*mddom%msfx(j,ice1))
          end do
        end if
        if ( ma%has_bdytop ) then
          do j = jce1 , jce2
            dpsdym(j,ice2) = (atm0%ps(j,ice2) - atm0%ps(j,ici2)) / &
                        (atm0%ps(j,ice2)*dx8*mddom%msfx(j,ice2))
          end do
        end if
        do k = 1 , kz
          do i = idi1 , idi2
            do j = jdi1 , jdi2
              atm0%dprddx(j,i,k) = atm0%pr(j,i,k)   - atm0%pr(j-1,i,k) + &
                                   atm0%pr(j,i-1,k) - atm0%pr(j-1,i-1,k)
              atm0%dprddy(j,i,k) = atm0%pr(j,i,k)   - atm0%pr(j,i-1,k) + &
                                   atm0%pr(j-1,i,k) - atm0%pr(j-1,i-1,k)
            end do
          end do
        end do
        if ( myid == italk ) then
          write(stdout,*) 'Reference atmosphere calculated.'
        end if
        ztop = maxval(atm0%zf)
        call maxall(ztop,rayzd)
        if ( myid == italk ) then
          write(stdout,*) 'Model top at ',rayzd,' m'
        end if
      end subroutine make_reference_atmosphere

      subroutine compute_full_coriolis_coefficients
        implicit none
        integer(ik4) :: i , j
        real(rkx) :: rotang , dlat , dlatdy , dlondy
        do i = idi1 , idi2
          do j = jdi1 , jdi2
            dlat = mddom%dlat(j,i)
            dlatdy = d_half * (mddom%xlat(j-1,i)   + mddom%xlat(j,i) - &
                               mddom%xlat(j-1,i-1) - mddom%xlat(j,i-1))
            if ( abs(dlatdy) < 1.0e-8_rkx ) dlatdy = sign(1.0e-8_rkx,dlatdy)
            dlondy = d_half * (mddom%xlon(j-1,i)   + mddom%xlon(j,i) - &
                               mddom%xlon(j-1,i-1) - mddom%xlon(j,i-1))
            if ( dlondy >  180.0_rkx ) dlondy = dlondy - 360.0_rkx
            if ( dlondy < -180.0_rkx ) dlondy = dlondy + 360.0_rkx
            rotang = -atan(dlondy/dlatdy*cos(degrad*dlat))
            if ( dlatdy < 0.0 ) rotang = rotang + mathpi
            mddom%ef(j,i) = eomeg2*cos(degrad*dlat)
            mddom%ddx(j,i) = cos(rotang)
            mddom%ddy(j,i) = sin(rotang)
            mddom%dmdx(j,i) = -d_half * &
              (mddom%msfx(j,i) + mddom%msfx(j,i-1) - &
               mddom%msfx(j-1,i) - mddom%msfx(j-1,i-1)) / &
               (dx*mddom%msfd(j,i)*mddom%msfd(j,i))
            mddom%dmdy(j,i) = -d_half * &
              (mddom%msfx(j,i) + mddom%msfx(j-1,i) - &
               mddom%msfx(j,i-1) - mddom%msfx(j-1,i-1)) / &
               (dx*mddom%msfd(j,i)*mddom%msfd(j,i))
          end do
        end do
        call exchange(mddom%ef,1,jdi1,jdi2,idi1,idi2)
        call exchange(mddom%ddx,1,jdi1,jdi2,idi1,idi2)
        call exchange(mddom%ddy,1,jdi1,jdi2,idi1,idi2)
        do i = ici1 , ici2
          do j = jci1 , jci2
            mddom%ex(j,i) = d_rfour*(mddom%ef(j,i)   + mddom%ef(j,i+1) + &
                                     mddom%ef(j+1,i) + mddom%ef(j+1,i+1))
            mddom%crx(j,i) = d_rfour*(mddom%ddx(j,i)   + mddom%ddx(j,i+1) + &
                                      mddom%ddx(j+1,i) + mddom%ddx(j+1,i+1))
            mddom%cry(j,i) = d_rfour*(mddom%ddy(j,i)   + mddom%ddy(j,i+1) + &
                                      mddom%ddy(j+1,i) + mddom%ddy(j+1,i+1))
          end do
        end do
        if ( myid == italk ) then
          write(stdout,*) 'Full Coriolis coefficients computed.'
        end if
      end subroutine compute_full_coriolis_coefficients

      subroutine compute_moloch_static
        implicit none
        integer(ik4) :: i , j
        call exchange_lr(mddom%msfu,1,jde1,jde2,ide1,ide2)
        call exchange_bt(mddom%msfv,1,jde1,jde2,ide1,ide2)
        do i = ice1 , ice2
          do j = jdi1 , jdi2
            mddom%hx(j,i) = (mddom%ht(j,i) - mddom%ht(j-1,i)) * &
                             mddom%msfu(j,i) * rdx * regrav
          end do
        end do
        if ( iproj == 'ROTLLR' ) then
          do i = idi1 , idi2
            do j = jce1 , jce2
              mddom%hy(j,i) = (mddom%ht(j,i) - mddom%ht(j,i-1)) * &
                               rdx * regrav
            end do
          end do
        else
          do i = idi1 , idi2
            do j = jce1 , jce2
              mddom%hy(j,i) = (mddom%ht(j,i) - mddom%ht(j,i-1)) * &
                               mddom%msfv(j,i) * rdx * regrav
            end do
          end do
        end if
        call exchange_lr(mddom%hx,1,jde1,jde2,ice1,ice2)
        call exchange_bt(mddom%hy,1,jce1,jce2,ide1,ide2)
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              mo_atm%zeta(j,i,k) = md_zeta(zitah(k),mddom%ht(j,i))
              mo_atm%fmz(j,i,k) = md_fmz(zitah(k),mddom%ht(j,i))
            end do
          end do
        end do
        call exchange_lrbt(mo_atm%fmz,1,jce1,jce2,ice1,ice2,1,kz)
        call exchange_lrbt(mo_atm%zeta,2,jce1,jce2,ice1,ice2,1,kz)
        do k = 1 , kzp1
          do i = ice1 , ice2
            do j = jce1 , jce2
              mo_atm%fmzf(j,i,k) = md_fmz(zita(k),mddom%ht(j,i))
              mo_atm%zetaf(j,i,k) = md_zeta(zita(k),mddom%ht(j,i))
            end do
          end do
        end do
        do concurrent ( j = jce1:jce2 , i = ice1:ice2 , k = 1:kz )
          mo_atm%dz(j,i,k) = mo_atm%zetaf(j,i,k) - mo_atm%zetaf(j,i,k+1)
        end do
        rayzd = mo_ztop
        if ( myid == italk ) then
          write(stdout,*) 'Model top at ',mo_ztop,' m'
        end if
      end subroutine compute_moloch_static

      recursive integer function gcd_rec(u,v) result(gcd)
        implicit none
        integer , intent(in) :: u , v
        if ( mod(u,v) /= 0 ) then
          gcd = gcd_rec(v,mod(u,v))
        else
          gcd = v
        end if
      end function gcd_rec

      real(rkx) function check_against_outparams(dt,dec) result(newdt)
        implicit none
        real(rkx) , intent(in) :: dt , dec
        newdt = int(dt/dec)*dec
        if ( ifshf ) then
          do
            if ( gcd_rec(int(newdt), int(secph)) < newdt ) then
              newdt = newdt + dec
              cycle
            end if
            exit
          end do
        else
          do
            if ( gcd_rec(int(newdt), int(secpd)) < newdt ) then
              newdt = newdt + dec
              cycle
            end if
            exit
          end do
        end if
      end function check_against_outparams

  end subroutine param

end module mod_params

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

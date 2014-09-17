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
  use mod_precip
  use mod_cloud_s1
  use mod_split
  use mod_slice
  use mod_bdycod
  use mod_ncio
  use mod_tendency
  use mod_ncout
  use mod_advection , only : init_advection
  use mod_savefile
  use mod_slabocean
  use mod_sldepparam

  implicit none

  private

  public :: param

  contains

  subroutine param

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine defines various model parameters.               c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  implicit none
!
  real(rk8) :: afracl , afracs , bb , cc , chibot , delsig , &
             dlargc , dsmalc , dxtemc , pk , ptmb , pz , qk ,       &
             qkp1 , sig700 , sigtbl , ssum , vqmax , vqrang , wk ,  &
             wkp1 , xbot , xtop , xx , yy
  integer(ik4) :: kbmax
  integer(ik4) :: iretval
  real(rk8) , dimension(nsplit) :: dtsplit
  integer(ik4) :: i , j , k , kbase , ktop , ns , mdate0 , mdate1 , mdate2
  integer(ik4) :: hspan
  integer(ik8) :: ndbgfrq , nsavfrq , natmfrq , nradfrq , nchefrq , nsrffrq
  integer(ik8) :: nlakfrq , nsubfrq , nbdyfrq , nslabfrq
  integer(ik4) :: n , len_path
  character(len=32) :: appdat
  type(rcm_time_interval) :: bdif
#ifdef DEBUG
  character(len=dbgslen) :: subroutine_name = 'param'
  integer(ik4) , save :: idindx = 0
#endif
!
!----------------------------------------------------------------------
!-----vqrang is the range limit on vqflx.
!
  data vqrang /5.0D-4/
!
!----------------------------------------------------------------------
!-----namelist:
!
  namelist /restartparam/ ifrest , mdate0 , mdate1 , mdate2

  namelist /timeparam/ dtrad , dtsrf , dtabem , dt

  namelist /outparam/ ifsave , ifatm , ifrad , ifsrf , ifsub , iflak , &
    ifsts , ifchem , ifopt , savfrq , atmfrq , srffrq , subfrq ,       &
    lakfrq , radfrq , chemfrq , enable_atm_vars ,                      &
    enable_srf_vars , enable_rad_vars , enable_sub_vars ,              &
    enable_sts_vars , enable_lak_vars , enable_opt_vars ,              &
    enable_che_vars , dirout , lsync , do_parallel_netcdf_in ,         &
    do_parallel_netcdf_out , idiag

  namelist /physicsparam/ ibltyp , iboudy , isladvec , icup_lnd ,    &
    icup_ocn , igcc , ipgf , iemiss , lakemod , ipptls , iocnflx ,   &
    iocncpl , iocnrough , ichem , scenario , idcsst , iseaice ,      &
    idesseas , iconvlwp , irrtm , iclimao3 , isolconst , icumcloud , &
    islab_ocean , itweak , temp_tend_maxval , wind_tend_maxval

  namelist /rrtmparam/ inflgsw , iceflgsw , liqflgsw , inflglw ,    &
    iceflglw , liqflglw , icld , irng , idrv

  namelist /subexparam/ ncld , qck1land , qck1oce , gulland , guloce ,     &
    rhmax , rh0oce , rh0land , cevaplnd , cevapoce , caccrlnd , caccroce , &
    tc0 , cllwcv , clfrcvmax , cftotmax , conf

  namelist /microparam/ stats , budget_compute , nssopt , kautoconv ,  &
    ksemi , vqxr , vqxi , vqxs , zauto_rate_khair , zauto_rate_kessl , &
    zauto_rate_klepi , rkconv , rcovpmin , rpecons

  namelist /grellparam/ shrmin , shrmax , edtmin , edtmax ,         &
    edtmino , edtmaxo , edtminx , edtmaxx , pbcmax , mincld ,       &
    htmin , htmax , skbmax , dtauc, shrmin_ocn , shrmax_ocn ,       &
    edtmin_ocn, edtmax_ocn, edtmino_ocn , edtmaxo_ocn ,             &
    edtminx_ocn , edtmaxx_ocn

  namelist /emanparam/ minsig , elcrit_ocn , elcrit_lnd , tlcrit ,  &
    entp , sigd , sigs , omtrain , omtsnow , coeffr , coeffs , cu , &
    betae , dtmax , alphae , damp , epmax_ocn , epmax_lnd

  namelist /tiedtkeparam/ iconv , entrmax , entrdd , entrpen ,   &
    entrscv , entrmid , cprcon , detrpen , entshalp , rcuc_lnd , &
    rcuc_ocn , rcpec_lnd , rcpec_ocn , rhebc_lnd , rhebc_ocn ,   &
    rprc_ocn , rprc_lnd

  namelist /kfparam/ kf_trigger , kf_entrate

  namelist /chemparam/ chemsimtype , ichremlsc , ichremcvc , ichdrdepo , &
         ichcumtra , ichsolver , idirect , iindirect , ichdustemd ,      &
         ichdiag , ichsursrc , ichebdy , rdstemfac, ichjphcld, ichbion

  namelist /uwparam/ iuwvadv , ilenparam , atwo , rstbl

  namelist /holtslagparam/ ricr_ocn , ricr_lnd , zhnew_fac , &
         ifaholtth10 , ifaholtmax , ifaholtmin

#ifdef CLM
  namelist /clmparam/ dirclm , imask , clmfrq , ilawrence_albedo
#endif

  namelist /cplparam/ cpldt

  namelist /slabocparam/ do_qflux_adj , do_restore_sst , &
    sst_restore_timescale , mixed_layer_depth

  namelist /tweakparam/ itweak_temperature , itweak_solar_irradiance , &
          itweak_sst , itweak_greenhouse_gases , temperature_tweak ,   &
          sst_tweak , solar_tweak , gas_tweak_factors
!
#ifdef DEBUG
  call time_begin(subroutine_name,idindx)
#endif
!
!----------------------------------------------------------------------
!-----specify the parameters used in the model:
!
!     rfstrt : whether the model is to be started with refined time
!     steps in the first hour. if ifrest=.true., this parameter
!     is never used.
!     = .true. ; yes
!     = .false. ; no
!
!     ifrest : whether this run is restarting from a saved tape
!     (unit 14).
!     = .true. ; yes
!     = .false. ; no
!
!     ifsave : specify whether a saved tape (unit 24) will be written
!     for restarting.
!     = .true. ; yes
!     = .false. ; no
!
!     savfrq : if ifsave=.true., specify the interval in hours
!     between save operations.
!
!     ibltyp : specify whether bats pbl or holtslag's pbl
!     parameterization is to be used in the model.
!     = 0  ; frictionless
!     = 1  ; Holtslag PBL (Holtslag, 1990)
!     = 2  ; UW PBL (Bretherton and McCaa, 2004)
!     = 99 ; Holtslag PBL, with UW in diag. mode
!
!     dtrad : specify the frequency in
!     minutes, the solar radiation will be computed in
!     subroutine "sfcrd". a typical value would be 30 minutes.
!
!     lakemod: inclusion of the hostetler lake model
!     = 0 ; no
!     = 1 ; yes
!
!     kxout  : the k level of the horizontal slice for printer output.
!
!     jxsex  : the j index of the north-south vertical slice for printer
!     output.
!
!     icup_lnd , icup_ocn   : type of cumulus parameterization over
!                             land and ocean
!     = 1 ; kuo
!     = 2 ; grell
!     = 3 ; betts-miller (1986)
!     = 4 ; emanuel (1991)
!     = 5 ; Tiedtke (1986) - version from ECHAM 5.4
!
!     igcc   : Grell Scheme Convective Closure Assumption
!     = 1 ; Arakawa & Schubert (1974)
!     = 2 ; Fritsch & Chappell (1980)
!
!     ipptls : type of moisture scheme
!     = 1 ; explicit moisture (SUBEX; Pal et al 2000)
!     = 2 ; new microphysics
!
!     iocnflx: type of ocean flux parameterization
!     = 1 ; BATS
!     = 2 ; Zeng et al.
!
!     iocncpl: controls the coupling with driver, RegESM
!     = 0 ; no coupling
!     = 1 ; activate coupling
!
!     iocnrough: Zeng Ocean Model Roughness model
!     = 1 ; (0.0065*ustar*ustar)/egrav (the RegCM V3 one)
!     = 2 ; (0.013*ustar*ustar)/egrav+0.11*visa/ustar (the Zeng one)
!
!     iboudy : specify the laterial boundary conditions.
!     = 0 ; fixed.
!     = 1 ; relaxation, linear technique.
!     = 2 ; time-dependent (from observations or large-scale
!     model).
!     = 3 ; time and inflow/outflow dependent.
!     = 4 ; sponge (perkey & kreitzberg, mwr 1976).
!     = 5 ; relaxation, exponential technique.
!
!     ifatm : whether you want output in a saved tape (unit 20) for
!     analyses in dataflow.
!     = .true. ; yes
!     = .false. ; no
!
!     atmfrq : if ifatm=true., specify the output interval in hours.
!
!     ifrad  : whether you want radiation output saved
!     = .true. ; yes
!     = .false. ; no
!
!     radfrq : if ifrad=1, specify the output interval in hours.
!
!     imask : Type of land surface parameterization
!     1= using DOMAIN.INFO for landmask (same as BATS);
!     2= using mksrf_navyoro file landfraction for landmask and perform
!        a weighted average over ocean/land gridcells; for example:
!        tgb = tgb_ocean*(1 - landfraction) + tgb_land*landfraction
!
!----------------------------------------------------------------------
!-----default values for all the options:
!     (can be overwritten by namelist input).
!------namelist restart param:
!
  ifrest = .false.     ! Restart?:  t=true; f=false
  idate0 = 1900010100  ! Start date of the initial simulation
  idate1 = 1900010100  ! Start date of this simulation
  idate2 = 1900010100  ! End Date this simulation
  ! note: beginning/end forecast time set in restart.mm4

!------namelist timeparam:
!
  dtrad = 30.0D0 ! time interval in min solar rad caluclated
  dtsrf = 600.0D0 ! time interval at which bats is called (secs)
  dtabem = 12.0D0  ! time interval absorption-emission calculated (hours)
  dt = 200.0D0    ! time step in seconds

!-----namelist out      note: * signifies currently not in namelist
!
  rfstrt = .false.      ! *
  ifsave = .false.
  ifatm  = .true.
  ifrad  = .true.
  ifsrf  = .true.
  ifsts  = .true.
  ifsub  = .false.
  iflak  = .false.
  ifopt  = .false.
  ifchem = .false.
  savfrq  = 24.0D0  ! time interval for disposing sav output (hrs)
  atmfrq  = 6.0D0   ! time interval for disposing atm output (hrs)
  radfrq  = 6.0D0   ! time interval for disposing rad output (hrs)
  srffrq  = 3.0D0   ! time interval for disposing srf output (hrs)
  lakfrq  = 6.0D0   ! time interval for disposing lake output (hrs)
  subfrq  = 6.0D0   ! time interval for disposing lake output (hrs)
  chemfrq = 6.0D0   ! time interval for disposing chem output (hrs)
  enable_atm_vars(:) = .true.
  enable_srf_vars(:) = .true.
  enable_sts_vars(:) = .true.
  enable_sub_vars(:) = .true.
  enable_lak_vars(:) = .true.
  enable_rad_vars(:) = .true.
  enable_opt_vars(:) = .true.
  enable_che_vars(:) = .true.
  dirout = './output'
  lsync = .true.
  do_parallel_netcdf_in = .false.
  do_parallel_netcdf_out = .false.

  idiag = 0
!
!----------------------------------------------------------------------
!-----namelist physicsparam:
!
  ibltyp = 1
  iboudy = 5
  isladvec = 0
  icup_lnd = 5
  icup_ocn = 5
  ipptls = 1
  igcc = 1
  ipgf = 1
  iemiss = 0
  iocnflx = 2
  iocncpl = 0
  iocnrough = 1
  lakemod = 0
  ichem = 0
  scenario = 'RF'
  idcsst = 0
  iseaice = 0
  idesseas = 0
  iconvlwp = 1
  irrtm = 0
  islab_ocean = 0
  iclimao3 = 0
  isolconst = 1
  icumcloud = 2
  temp_tend_maxval = 0.5*(dt/secpm)
  wind_tend_maxval = 0.5*(dt/secpm)
!----------------------------------------------------------------------
!-----rrtm radiation namelist param:
!
  inflgsw  = 2
  iceflgsw = 3
  liqflgsw = 1
  inflglw  = 2
  iceflglw = 3
  liqflglw = 1
  idrv = 0
  icld  = 1
  irng = 1
!----------------------------------------------------------------------
!------namelist subexparam:
  ncld = 1             ! # of bottom model levels with no clouds (rad only)
  qck1land = 0.00025D0 ! Autoconversion Rate for Land
  qck1oce = 0.00025D0  ! Autoconversion Rate for Ocean
  gulland = 0.4D0      ! Fract of Gultepe eqn (qcth) when prcp occurs (land)
  guloce = 0.4D0       ! Fract of Gultepe eqn (qcth) for ocean
  rhmax = 1.01D0       ! RH at whicn FCC = 1.0
  rh0oce = 0.90D0      ! Relative humidity threshold for ocean
  rh0land = 0.80D0     ! Relative humidity threshold for land
  tc0 = 238.0D0        ! Below this temp, rh0 begins to approach unity
  cevaplnd = 1.0D-3    ! Raindrop evap rate coef land [[(kg m-2 s-1)-1/2]/s]
  cevapoce = 1.0D-3    ! Raindrop evap rate coef ocean [[(kg m-2 s-1)-1/2]/s]
  caccrlnd = 3.0D0     ! Raindrop accretion rate land [m3/kg/s]
  caccroce = 3.0D0     ! Raindrop accretion rate ocean [m3/kg/s]
  cllwcv = 0.3D-3      ! Cloud liquid water content for convective precip.
  clfrcvmax = 1.00D0   ! Max cloud fractional cover for convective precip.
  cftotmax = 0.75D0    ! Max total cover cloud fraction for radiation
  conf = d_one         !   Condensation threshold

  ! namelist microparam
  stats = .false.
  budget_compute = .false. ! Verify enthalpy and moisture conservation
  nssopt = 1 ! Supersaturation Computation
             ! 0 => No scheme
             ! 1 => Tompkins
             ! 2 => Lohmann and Karcher
             ! 3 => Gierens
  kautoconv = 4 !  Choose the autoconversion paramaterization
                ! => 1 Klein & Pincus (2000)
                ! => 2 Khairoutdinov and Kogan (2000)
                ! => 3 Kessler (1969)
                ! => 4 Sundqvist
  ksemi = d_one ! 0 => fully explicit
                ! 1 => fully implict
                ! 0<ksemi<1 => semi-implicit
  vqxr = 4.0D0
  vqxi = 0.15D0
  vqxs = 1.0D0
  zauto_rate_khair = 0.355D0
  zauto_rate_kessl = 1.D-3
  zauto_rate_klepi = 0.5D-3
  rkconv = 1.666D-4 ! d_one/6000.0D0
  rcovpmin = 0.1D0
  rpecons = 5.547D-5
!
!------namelist grellparam:
  shrmin = 0.25D0       ! Minimum Shear effect on precip eff.
  shrmax = 0.50D0       ! Maximum Shear effect on precip eff.
  edtmin = 0.25D0       ! Minimum Precipitation Efficiency
  edtmax = 0.50D0       ! Maximum Precipitation Efficiency
  edtmino = 0.25D0      ! Minimum Precipitation Efficiency (o var)
  edtmaxo = 0.50D0      ! Maximum Precipitation Efficiency (o var)
  edtminx = 0.25D0      ! Minimum Precipitation Efficiency (x var)
  edtmaxx = 0.50D0      ! Maximum Precipitation Efficiency (x var)
  shrmin_ocn = 0.25D0   ! Minimum Shear effect on precip eff.
  shrmax_ocn = 0.50D0   ! Maximum Shear effect on precip eff.
  edtmin_ocn = 0.25D0   ! Minimum Precipitation Efficiency
  edtmax_ocn = 0.50D0   ! Maximum Precipitation Efficiency
  edtmino_ocn = 0.25D0  ! Minimum Precipitation Efficiency (o var)
  edtmaxo_ocn = 0.50D0  ! Maximum Precipitation Efficiency (o var)
  edtminx_ocn = 0.25D0  ! Minimum Precipitation Efficiency (x var)
  edtmaxx_ocn = 0.50D0  ! Maximum Precipitation Efficiency (x var)
  pbcmax = 150.0D0       ! Max depth (mb) of stable layer b/twn LCL & LFC
  mincld = 150.0D0       ! Min cloud depth (mb).
  htmin = -250.0D0       ! Min convective heating
  htmax = 500.0D0        ! Max convective heating
  skbmax = 0.4D0        ! Max cloud base height in sigma
  dtauc = 60.0D0         ! Fritsch & Chappell (1980)
!
!------namelist emanparam:
  minsig = 0.95D0     ! Lowest sigma level from which convection can originate
  elcrit_ocn = 0.0011D0 ! Autoconversion threshold water content (gm/gm)
  elcrit_lnd = 0.0011D0 ! Autoconversion threshold water content (gm/gm)
  tlcrit = -55.0D0    ! Below tlcrit auto-conversion threshold is zero
  entp = 1.5D0        ! Coefficient of mixing in the entrainment formulation
  sigd = 0.05D0       ! Fractional area covered by unsaturated dndraft
  sigs = 0.12D0       ! Fraction of precipitation falling outside of cloud
  omtrain = 50.0D0    ! Fall speed of rain (P/s)
  omtsnow = 5.5D0     ! Fall speed of snow (P/s)
  coeffr = 1.0D0      ! Coefficient governing the rate of rain evaporation
  coeffs = 0.8D0      ! Coefficient governing the rate of snow evaporation
  cu = 0.7D0          ! Coefficient governing convective momentum transport
  betae = 10.0D0      ! Controls downdraft velocity scale
  dtmax = 0.9D0       ! Max negative parcel temperature perturbation below LFC
  alphae = 0.2D0      ! Controls the approach rate to quasi-equilibrium
  damp = 0.1D0        ! Controls the approach rate to quasi-equilibrium
  epmax_ocn = 0.999D0 ! Maximum precipitation efficiency over land
  epmax_lnd = 0.999D0 ! Maximum precipitation efficiency over ocean
!
!------namelist tiedtkeparam:
  iconv    = 4        ! Selects the actual scheme
  entrmax  = 1.75D-3  ! Max entrainment iconv=[1,2,3]
  entrdd   = 3.0D-4   ! Entrainment rate for cumulus downdrafts
  entrpen  = 1.75D-3  ! Entrainment rate for penetrative convection
  entrscv  = 3.0D-4   ! Entrainment rate for shallow convection iconv=[1,2,3]
  entrmid  = 1.0D-4   ! Entrainment rate for midlevel convection iconv=[1,2,3]
  cprcon   = 1.0D-4   ! Coefficient for determining conversion iconv=[1,2,3]
  detrpen = 0.75D-4   ! Detrainment rate for penetrative convection
  entshalp = 2.0D0    ! shallow entrainment factor for entrpen
  rcuc_lnd = 0.05D0   ! Convective cloud cover for rain evporation
  rcuc_ocn = 0.05D0   ! Convective cloud cover for rain evporation
  rcpec_ocn = 5.55D-5 ! Coefficient for rain evaporation below cloud
  rcpec_lnd = 5.55D-5 ! Coefficient for rain evaporation below cloud
  rhebc_lnd = 0.7D0   ! Critical relative humidity below
                      ! cloud at which evaporation starts for land
  rhebc_ocn = 0.9D0   ! Critical relative humidity below
                      ! cloud at which evaporation starts for ocean
  rprc_lnd = 1.4D-3   ! coefficient for determining conversion from cloud water
  rprc_ocn = 1.4D-3   ! coefficient for determining conversion from cloud water

!c------namelist kfparam ;

  kf_trigger = 3
  kf_entrate = 0.03D0 ! Kain Fritsch entrainment rate

!c------namelist uwparam ;
  iuwvadv = 0
  ilenparam = 0
  atwo = 15.0D0
  rstbl = 1.5D0

!c------namelist holtslagparam ;

  ricr_ocn = 0.25D0
  ricr_lnd = 0.25D0
  zhnew_fac = 0.25D0
  ifaholtth10 = 1
  ifaholtmax = 1
  ifaholtmin = 0

!c-----namelist slabocparam ;

  mixed_layer_depth     = 50.0D0
  sst_restore_timescale = 5.0D0 !days
  do_restore_sst = .true.
  do_qflux_adj = .false.

!c-----namelist tweakparam ;

  itweak = 0
  itweak_sst = 0
  itweak_temperature = 0
  itweak_solar_irradiance = 0
  itweak_greenhouse_gases = 0
  sst_tweak = 0.0D0
  temperature_tweak = 0.0D0
  solar_tweak = 0.0D0
  gas_tweak_factors(:) = 1.0D0

!c------namelist chemparam ; ( 0= none, 1= activated)
  ichsolver = 1     ! enable chem solver
  ichremlsc = 1     ! tracer removal by large scale clouds
  ichremcvc = 1     ! tracer removal by convective clouds
  ichdrdepo = 1     ! tracer dry deposition
  ichcumtra = 1     ! tracer convective transport
  ichdustemd = 1    ! dust emission distribution (1 = alfaro, 2 =kok)
  ichjphcld = 1     ! impact of cloud aod on photolysis coef
  idirect = 1       ! tracer direct effect
  iindirect = 0
  ichdiag = 0       ! chem tend outputs
  ichsursrc = 1
  ichebdy = 1
  rdstemfac = d_one
  ichbion = 0
!c------namelist clmparam ; (read in case clm surface model compiled in)
#ifdef CLM
  imask = 1
  ilawrence_albedo = 1
  clmfrq = 24.0D0
#endif
!c------namelist cplparam ;
  cpldt = 21600.0D0       ! coupling time step in seconds (seconds)
!
!---------------------------------------------------------------------
!
#ifdef CLM
  if ( myid == italk ) then
    if (nsg /= 1 ) then
      write (stderr,*) 'Running SUBGRID with CLM: not implemented'
      write (stderr,*) 'Please set nsg to 1 in regcm.in'
      call fatal(__FILE__,__LINE__,'CLM & SUBGRID TOGETHER')
    end if
  end if
#endif

  if ( myid == iocpu ) then
    open(ipunit, file=namelistfile, status='old', &
                 action='read', iostat=iretval)
    if ( iretval /= 0 ) then
      write(stderr,*) 'Error opening input namelist file ',trim(namelistfile)
      call fatal(__FILE__,__LINE__,'INPUT NAMELIST OPEN ERROR')
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
      call fatal(__FILE__,__LINE__,'INPUT NAMELIST READ ERROR')
#ifdef DEBUG
    else
      write(stdout,*) 'Read restartparam OK'
#endif
    end if

    idate0 = mdate0
    idate1 = mdate1
    idate2 = mdate2
    call setcal(idate0,ical)
    call setcal(idate1,ical)
    call setcal(idate2,ical)
    bdif = idate2 - idate1
    hspan = idnint(tohours(bdif))
    if ( mod(hspan,24) /= 0 ) then
      call fatal(__FILE__,__LINE__,  &
                 'Runtime increments must be modulus 24 hours')
    end if

    rewind(ipunit)
    read (ipunit, nml=timeparam, iostat=iretval, err=101)
    if ( iretval /= 0 ) then
      write(stderr,*) 'Error reading timeparam namelist'
      call fatal(__FILE__,__LINE__,'INPUT NAMELIST READ ERROR')
#ifdef DEBUG
    else
      write(stdout,*) 'Read timeparam OK'
#endif
    end if

    rewind(ipunit)
    read (ipunit, nml=outparam, iostat=iretval, err=102)
    if ( iretval /= 0 ) then
      write(stderr,*) 'Error reading outparam namelist'
      call fatal(__FILE__,__LINE__,'INPUT NAMELIST READ ERROR')
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
      call fatal(__FILE__,__LINE__,'INPUT NAMELIST READ ERROR')
#ifdef DEBUG
    else
      write(stdout,*) 'Read physicsparam OK'
#endif
    end if

    ! Hack. permanently disable seasonal albedo.
    idesseas = 0

    icup(1) = icup_lnd
    icup(2) = icup_ocn
    if ( any(icup == 1) ) then
      if ( icup_lnd /= icup_ocn ) then
        write(stderr,*) 'ERROR: Kuo scheme MUST be used on both Land and Ocean'
        call fatal(__FILE__,__LINE__,'INPUT NAMELIST ICUP INCONSISTENT')
      end if
    end if
    if ( any(icup == 3) ) then
      if ( icup_lnd /= icup_ocn ) then
        write(stderr,*) 'ERROR: BM scheme MUST be used on both Land and Ocean'
        call fatal(__FILE__,__LINE__,'INPUT NAMELIST ICUP INCONSISTENT')
      end if
    end if

    if ( ipptls > 0 ) then
      rewind(ipunit)
      read (ipunit, nml=subexparam, iostat=iretval, err=104)
      if ( iretval /= 0 ) then
        write(stdout,*) 'Using default subex parameter.'
#ifdef DEBUG
      else
        write(stdout,*) 'Read subexparam OK'
#endif
      end if
      if ( ipptls == 2 ) then
        rewind(ipunit)
        read (ipunit, nml=microparam, iostat=iretval, err=105)
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
    end if

    if ( any(icup == 2) ) then
      rewind(ipunit)
      read (ipunit, nml=grellparam, iostat=iretval, err=106)
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
      read (ipunit, nml=emanparam, iostat=iretval, err=107)
      if ( iretval /= 0 ) then
        write(stdout,*) 'Using default MIT parameter.'
#ifdef DEBUG
      else
        write(stdout,*) 'Read emanparam OK'
#endif
      end if
    end if
    if ( any(icup == 5) ) then
      rewind(ipunit)
      read (ipunit, nml=tiedtkeparam, iostat=iretval, err=108)
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
      read (ipunit, nml=kfparam, iostat=iretval, err=109)
      if ( iretval /= 0 ) then
        write(stdout,*) 'Using default Kain Fritsch parameter.'
#ifdef DEBUG
      else
        write(stdout,*) 'Read kfparam OK'
#endif
      end if
    end if
    if ( any(icup < 0) .or. any(icup > 6) ) then
      call fatal(__FILE__,__LINE__,'UNSUPPORTED CUMULUS SCHEME')
    end if
    if ( ibltyp < 0 .or. (ibltyp > 2 .and. ibltyp /= 99) ) then
      call fatal(__FILE__,__LINE__,'UNSUPPORTED PBL SCHEME.')
    end if
    if ( ibltyp == 1 .or. ibltyp == 99 ) then
      rewind(ipunit)
      read (ipunit, nml=holtslagparam, iostat=iretval, err=110)
      if ( iretval /= 0 ) then
        write(stdout,*) 'Using default Holtslag parameter.'
#ifdef DEBUG
      else
        write(stdout,*) 'Read holtslagparam OK'
#endif
      end if
    end if
    if ( ibltyp == 2 .or. ibltyp == 99 ) then
      rewind(ipunit)
      read (ipunit, nml=uwparam, iostat=iretval, err=111)
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
      read (ipunit, nml=rrtmparam, iostat=iretval, err=112)
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
      read (ipunit, nml=slabocparam, iostat=iretval, err=113)
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
        call fatal(__FILE__,__LINE__,'SLABOCEAN INPUT INCONSISTENCY')
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
      rewind(ipunit)
      read (ipunit, chemparam, iostat=iretval, err=114)
      if ( iretval /= 0 ) then
        write(stderr,*) 'Error reading chemparam namelist'
        call fatal(__FILE__,__LINE__,'INPUT NAMELIST READ ERROR')
#ifdef DEBUG
      else
        write(stdout,*) 'Read chemparam OK'
#endif
      end if
      ! force option 2 for drydep in the case of UW PBL
      if ( ibltyp == 2 .or. ibltyp == 99 ) ichdrdepo = 2
    else
      ichem = 0
      ntr = 0
    end if
#ifdef CLM
    rewind(ipunit)
    read (ipunit , clmparam, iostat=iretval, err=115)
    if ( iretval /= 0 ) then
      write(stdout,*) 'Using default CLM parameter.'
#ifdef DEBUG
    else
      write(stdout,*) 'Read clmparam OK'
#endif
    end if
#endif
    if ( iocncpl == 1 ) then
      rewind(ipunit)
      read (ipunit , cplparam, iostat=iretval, err=116)
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
      read (ipunit , tweakparam, iostat=iretval, err=117)
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
  end if
!
!  communicate to all processors
!
  call bcast(ifrest)
  call bcast(hspan)
  call bcast(idate0)
  call bcast(idate1)
  call bcast(idate2)
  call bcast(globidate1)
  call bcast(globidate2)

  call bcast(dtrad)
  call bcast(dtabem)
  call bcast(dtsrf)
  call bcast(dt)

  call bcast(ifsave)
  call bcast(savfrq)
  call bcast(ifatm)
  call bcast(atmfrq)
  call bcast(ifrad)
  call bcast(radfrq)
  call bcast(ifsrf)
  call bcast(ifsub)
  call bcast(iflak)
  call bcast(ifsts)
  call bcast(ifopt)
  call bcast(ifchem)
  call bcast(srffrq)
  call bcast(lakfrq)
  call bcast(subfrq)
  call bcast(chemfrq)
  call bcast(lsync)
  call bcast(idiag)
  call bcast(do_parallel_netcdf_in)
#ifdef NETCDF4_HDF5
  call bcast(do_parallel_netcdf_out)
#else
  do_parallel_netcdf_out = .false.
#endif

  call bcast(iboudy)
  call bcast(isladvec)
  call bcast(ibltyp)
  call bcast(icup_lnd)
  call bcast(icup_ocn)
  call bcast(icup)
  call bcast(igcc)
  call bcast(ipptls)
  call bcast(iocnflx)
  call bcast(iocncpl)
  call bcast(iocnrough)
  call bcast(ipgf)
  call bcast(iemiss)
  call bcast(lakemod)
  call bcast(ichem)
  call bcast(ntr)

  if ( ipptls == 2 ) then
    nqx = 5
    iqfrst = iqc
    iqlst  = iqs
  else
    nqx = 2
    iqfrst = iqc
    iqlst  = iqc
  end if

  ! Check if really do output

#ifdef CLM
  if ( lakemod /= 0 ) then
    if ( myid == italk ) then
      write(stderr,*) 'Disabling BATS lake model, this is a CLM run'
    end if
    lakemod = 0
  end if
  if ( iemiss /= 0 ) then
    if ( myid == italk ) then
      write(stderr,*) 'Disabling Surface Emissivity, this is a CLM run'
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
    ifopt = .false.
  end if
  !
  ! Force the correct scenario from dattyp in CMIP5
  !
  if ( myid == iocpu ) then
    if ( dattyp(4:5) == '26' ) then
      if ( scenario /= 'RCP3PD' .or. scenario /= 'RCP2.6' ) then
        write(stderr,*) 'Forcing scenario from dattyp to RCP2.6'
        scenario = 'RCP3PD'
      end if
    else if ( dattyp(4:5) == '45' ) then
      if ( scenario /= 'RCP4.5') then
        write(stderr,*) 'Forcing scenario from dattyp to RCP4.5'
        scenario = 'RCP4.5'
      end if
    else if ( dattyp(4:5) == '60') then
      if ( scenario /= 'RCP6' .or.  scenario /= 'RCP60' .or.  &
           scenario /= 'RCP60') then
        write(stderr,*) 'Forcing scenario from dattyp to RCP6.0'
        scenario = 'RCP6'
      end if
    else if ( dattyp(4:5) == '85') then
      if ( scenario /= 'RCP8.5') then
        write(stderr,*) 'Forcing scenario from dattyp to RCP8.5'
        scenario = 'RCP8.5'
      end if
    end if
  end if

  call bcast(scenario,8)
  call bcast(idcsst)
  call bcast(iseaice)
  call bcast(idesseas)
  call bcast(iconvlwp)
  call bcast(irrtm)
  call bcast(iclimao3)
  call bcast(isolconst)
  call bcast(icumcloud)
  call bcast(islab_ocean)
  call bcast(itweak)
  if ( idcsst == 1 .and. iocnflx /= 2 ) then
    if ( myid == italk ) then
      write(stderr,*) 'Cannot enable diurnal cycle sst without Zheng ocean flux'
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

  if (iocncpl == 1) then
    call bcast(cpldt)
  end if

  if ( ipptls > 0 ) then
    call bcast(ncld)
    call bcast(qck1land)
    call bcast(qck1oce)
    call bcast(gulland)
    call bcast(guloce)
    call bcast(rhmax)
    call bcast(rh0oce)
    call bcast(rh0land)
    call bcast(tc0)
    call bcast(cevaplnd)
    call bcast(cevapoce)
    call bcast(caccrlnd)
    call bcast(caccroce)
    call bcast(cftotmax)
    call bcast(conf)
    if ( ipptls == 2 ) then
      call bcast(stats)
      call bcast(budget_compute)
      call bcast(nssopt)
      call bcast(kautoconv)
      call bcast(ksemi)
      call bcast(vqxr)
      call bcast(vqxi)
      call bcast(vqxs)
      call bcast(zauto_rate_khair)
      call bcast(zauto_rate_kessl)
      call bcast(zauto_rate_klepi)
      call bcast(rcovpmin)
      call bcast(rpecons)
      call bcast(rkconv)
    end if
  end if

  if ( irrtm == 1 ) then
    call bcast(inflgsw)
    call bcast(iceflgsw)
    call bcast(liqflgsw)
    call bcast(inflglw)
    call bcast(iceflglw)
    call bcast(liqflglw)
    call bcast(idrv)
    call bcast(icld)
    call bcast(irng)
  end if

  call bcast(clfrcvmax)
  call bcast(cllwcv)

  if ( any(icup == 2) ) then
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
    call bcast(minsig)
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
  end if

  if ( any(icup == 5) ) then
    call bcast(iconv)
    call bcast(entrmax)
    call bcast(entrdd)
    call bcast(entrpen)
    call bcast(entrscv)
    call bcast(entrmid)
    call bcast(cprcon)
    call bcast(detrpen)
    call bcast(entshalp)
    call bcast(rcuc_lnd)
    call bcast(rcuc_ocn)
    call bcast(rcpec_lnd)
    call bcast(rcpec_ocn)
    call bcast(rhebc_lnd)
    call bcast(rhebc_ocn)
    call bcast(rprc_lnd)
    call bcast(rprc_ocn)
  end if

  if ( any(icup == 6) ) then
    call bcast(kf_entrate)
  end if

  if ( ibltyp == 1 .or. ibltyp == 99 ) then
    call bcast(ricr_ocn)
    call bcast(ricr_lnd)
    call bcast(zhnew_fac)
    call bcast(ifaholtth10)
    call bcast(ifaholtmax)
    call bcast(ifaholtmin)
  end if
  if ( ibltyp == 2 .or. ibltyp == 99 ) then
    call bcast(iuwvadv)
    call bcast(ilenparam)
    call bcast(atwo)
    call bcast(rstbl)
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
    call bcast(iindirect)
    call bcast(ichsolver)
    call bcast(ichjphcld)
    call bcast(ichdustemd)
    call bcast(rdstemfac)
    call bcast(ichdiag)
    call bcast(ichsursrc)
    call chem_config
    call bcast(ntr)
    call bcast(iaerosol)
    call bcast(ioxclim)
    call bcast(igaschem)
    call bcast(ichebdy)
    call bcast(ichbion)
  end if
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------ALLOCATE NEEDED SPACE---------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!

  call allocate_mod_runparams

  call allocate_mod_atm_interface(ibltyp,isladvec)

  call allocate_mod_tend

  call allocate_mod_bdycon(iboudy)

  call allocate_pblscheme

  if ( ipptls > 0 ) then
    call allocate_mod_precip(ichem)
  end if

  call allocate_mod_split

  call allocate_mod_savefile

  call allocate_surface_model

  if ( ipptls == 2 ) then
    call allocate_mod_cloud_s1
  end if

  call allocate_radiation

  call allocate_mod_che_common(isladvec)
  call allocate_mod_che_mppio
  call allocate_mod_che_dust
  call allocate_mod_che_bdyco
  call allocate_mod_che_bionit
  if ( isladvec == 1 ) then
    call allocate_mod_sldepparam
  end if
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------ALLOCATE NEEDED SPACE---------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
  if ( myid == italk ) then
    if ( mod(idnint(dtrad*60.0D0),idnint(dt)) /= 0 ) then
      write (stderr,*) 'DTRAD=' , dtrad , 'DT=' , dt
      call fatal(__FILE__,__LINE__, &
                 'INCONSISTENT RADIATION TIMESTEPS SPECIFIED')
    end if
    if ( mod(idnint(dtsrf),idnint(dt)) /= 0 ) then
      write (stderr,*) 'DTSRF=' , dtsrf , 'DT=' , dt
      call fatal(__FILE__,__LINE__, &
                 'INCONSISTENT SURFACE TIMESTEPS SPECIFIED')
    end if
    if ( ifsrf ) then
      if ( mod(idnint(srffrq*secph),idnint(dtsrf)) /= 0 ) then
        write (stderr,*) 'BATFRQ=' , srffrq , 'DTSRF=' , dtsrf
        call fatal(__FILE__,__LINE__, &
                   'INCONSISTENT SURFACE OUTPUT FREQUENCY SPECIFIED')
      end if
      if ( ifsts .and. srffrq > 24.0D0 ) then
        call fatal(__FILE__,__LINE__, &
                   'NEED SRF FREQUENCY LESS THAN 24H FOR STS OUTPUT')
      end if
    else
      if ( ifsts ) then
        call fatal(__FILE__,__LINE__, &
                   'TO ENABLE STS, ENABLE SRF OUTPUT IS REQUIRED')
      end if
    end if
    if ( mod(idnint(dtabem*secph),idnint(dt)) /= 0 ) then
      write (stderr,*) 'DTABEM=' , dtabem , 'DT=' , dt
      call fatal(__FILE__,__LINE__, &
                 'INCONSISTENT ABS/EMS TIMESTEPS SPECIFIED')
    end if
    if ( mod(idnint(dtabem*60.0D0),idnint(dtrad)) /= 0 ) then
      write (stderr,*) 'DTABEM=' , dtabem , 'DTRAD=' , dtrad
      call fatal(__FILE__,__LINE__,                                   &
                 'INCONSISTENT LONGWAVE/SHORTWAVE RADIATION'//        &
                 ' TIMESTEPS SPECIFIED')
    end if
    if ( ichem == 1 .and. chemfrq <=  d_zero) then
      write (stderr,*) 'CHEMFRQ=' ,chemfrq
      call fatal(__FILE__,__LINE__,'CHEMFRQ CANNOT BE ZERO')
    end if
    if ( isladvec == 1 ) then
      if ( jxp < 5 .or. iyp < 5 ) then
        write (stderr,*) 'To use Semi-Lagrangian Advection Scheme reduce'
        write (stderr,*) 'the number of processors !!!!'
        write (stderr,*) 'Minimum number of points is 25 (5x5) per processor'
        call fatal(__FILE__,__LINE__,'ISLADVEC WITH PPROC < 5x5')
      end if
    end if
  end if
!
  if ( ifrest ) then
    doing_restart = .true.
  end if
!
!.....calculate the time step in minutes.
!
  dtsec = dt
  dtbat = dt
  rdt   = d_one/dt
  dtbdys = dble(ibdyfrq)*secph
  ntsec = idnint(dt)
!
!-----reset the options/calculate variables using namelist info:
!
  bdydate1 = idate1

  nsavfrq = idnint(secph*savfrq)
  natmfrq = idnint(secph*atmfrq)
  nradfrq = idnint(secph*radfrq)
  ndbgfrq = idnint(secph*dbgfrq)
  nsrffrq = idnint(secph*srffrq)
  nlakfrq = idnint(secph*lakfrq)
  nsubfrq = idnint(secph*subfrq)
  nchefrq = idnint(secph*chemfrq)
  nslabfrq = idnint(dtbdys)
  nbdyfrq = idnint(dtbdys)

  cfdout =  dtsec/(secph*chemfrq)
  afdout =  dtsec/(secph*atmfrq)

  ntsrf = idnint(dtsrf/dtsec)
  rtsrf = d_one/dble(ntsrf)
  ntrad = idnint((dtrad*secpm)/dtsec)
  rtrad = d_one/dble(ntrad)

  ktau = 0

  khour = 3600_8/idnint(dtsec)
  kday  = 86400_8/idnint(dtsec)
  krep  = khour*3
  kbdy  = nbdyfrq/idnint(dtsec)
  katm  = natmfrq/idnint(dtsec)
  ksrf  = nsrffrq/idnint(dtsec)
  klak  = nlakfrq/idnint(dtsec)
  ksub  = nsubfrq/idnint(dtsec)
  ksts  = khour*24
  krad  = nradfrq/idnint(dtsec)
  kche  = nchefrq/idnint(dtsec)
  kdbg  = ndbgfrq/idnint(dtsec)
  ksav  = nsavfrq/idnint(dtsec)

  rnsrf_for_srffrq = d_one/(dble(ksrf)*rtsrf)
  rnsrf_for_lakfrq = d_one/(dble(klak)*rtsrf)
  rnsrf_for_subfrq = d_one/(dble(ksub)*rtsrf)
  rnsrf_for_day = d_one/(dble(kday)*rtsrf)
  rnrad_for_radfrq = d_one/(dble(krad)*rtrad)
  rnrad_for_chem = dble(ntrad)/dble(kche)

  rsrf_in_atm = dble(ntsrf)/dble(katm)
  rsrffrq_sec = d_one/(srffrq*secph)

  mtau = idnint((hspan*secph)/dt)

  do ns = 1 , nsplit
    dtsplit(ns) = dt*(d_half/dble(nsplit-ns+1))
    dtau(ns) = dtsplit(ns)
  end do
  ntabem = idnint(secph*dtabem/dt) !dtabem is time interval abs./emis. calc.
  dt2 = d_two*dt
  dtsq = dt*dt
  dtcb = dt*dt*dt
!
  intmdl = rcm_time_interval(idnint(dt),usec)
  intbdy = rcm_time_interval(ibdyfrq,uhrs)
  intsom = rcm_time_interval(1,umnt)
  deltmx = dt
!
  if ( myid == italk ) then
    appdat = tochar(idate0)
    write(stdout,*) 'Initial date of the global simulation: ',appdat
    appdat = tochar(idate1)
    write(stdout,*) 'Initial date of this run             : ',appdat
    appdat = tochar(idate2)
    write(stdout,*) 'Final date of this run               : ',appdat
    write(stdout,*) 'Total simulation lenght              : ',hspan,' hours'
    write(stdout,'(a,f11.6)') ' Timestep in seconds = ',dtsec
    write(stdout,'(a,2f11.6)') ' Split explicit dtau = ',dtau
  end if

  call bcast(dirter,256)
  call bcast(dirglob,256)
  call bcast(dirout,256)
  call bcast(domname,64)
  call read_domain_info(mddom%ht,mddom%lndcat,mddom%mask, &
                        mddom%xlat,mddom%xlon,mddom%dlat,mddom%dlon, &
                        mddom%msfx,mddom%msfd,mddom%coriol, &
                        mddom%snowam,mddom%dhlake)
  call bcast(ds)
  call bcast(ptop)

  dx = ds * d_1000
!
!-----specify the constants used in the model.
!     conf   : condensation threshold.
!     qcth   : threshold for the onset of autoconversion.
!     qck1oce  : constant autoconversion rate for ocean.
!     qck1land : constant autoconversion rate for land.
!     all the other constants are used to compute the cloud
!     microphysical parameterization (ref. orville & kopp, 1977 jas).
!
  dx2 = d_two*dx
  dx4 = d_four*dx
  dx8 = 8.0D0*dx
  dx16 = 16.0D0*dx
  dxsq = dx*dx
  c200 = vonkar*vonkar*dx/(d_four*(d_100-ptop))
  rdxsq = 1.0D0/dxsq
  xkhz = 1.5D-3*dxsq/dt
  xkhmax = dxsq/(64.0D0*dt)
  akht1 = dxsq/tauht
  akht2 = dxsq/tauht
!
  ! Calculate boundary areas per processor

  call setup_boundaries(cross,ba_cr)
  call setup_boundaries(dot,ba_dt)

  call allocate_v2dbound(xpsb,cross)
  call allocate_v3dbound(xtb,kz,cross)
  call allocate_v3dbound(xqb,kz,cross)
  call allocate_v3dbound(xub,kz,dot)
  call allocate_v3dbound(xvb,kz,dot)
!
  if ( myid == italk ) write(stdout,*) 'Setting IPCC scenario to ', scenario
  call set_scenario(scenario)

!
  if ( myid == italk ) then
    if ( ifrest .and. idate0 == idate1 ) then
      write(stderr,*) 'Error in parameter set.'
      write(stderr,*) 'Cannot set idate0 == idate1 on restart run'
      write(stderr,*) 'Correct idate0.'
      call fatal(__FILE__,__LINE__,'IDATE0==IDATE1 ON RESTART')
    else if ( .not. ifrest .and. idate0 /= idate1 ) then
      write(stderr,*) 'Error in parameter set.'
      write(stderr,*) 'Cannot set idate0 /= idate1 on non restart run'
      write(stderr,*) 'Correct idate1.'
      call fatal(__FILE__,__LINE__,'IDATE0/=IDATE1 ON NON RESTART')
    end if
  end if
!
  idatex = idate1
  call split_idate(idatex,xyear,xmonth,xday,xhour)
  kstsoff = khour*xhour
!
  if ( myid == italk ) then
    write(stdout,*) 'Create SAV files : ' , ifsave
    write(stdout,*) 'Create ATM files : ' , ifatm
    write(stdout,*) 'Create RAD files : ' , ifrad
    write(stdout,*) 'Create SRF files : ' , ifsrf
    write(stdout,*) 'Create STS files : ' , ifsts
    if ( nsg > 1 ) write(stdout,*) 'Create SUB files : ' , ifsub
    if ( lakemod == 1 ) write(stdout,*) 'Create LAK files : ' , iflak
    if ( ichem == 1 ) then
      write(stdout,*) 'Create CHE files : ' , ifchem
      write(stdout,*) 'Create OPT files : ' , ifopt
    end if
    write(stdout,'(a,f6.1)') ' Frequency in hours to create SAV : ' , savfrq
    write(stdout,'(a,f6.1)') ' Frequency in hours to create ATM : ' , atmfrq
    write(stdout,'(a,f6.1)') ' Frequency in hours to create RAD : ' , radfrq
    write(stdout,'(a,f6.1)') ' Frequency in hours to create SRF : ' , srffrq
    if ( nsg > 1 ) &
      write(stdout,'(a,f6.1)') ' Frequency in hours to create SUB : ' , subfrq
    if ( lakemod == 1 ) &
      write(stdout,'(a,f6.1)') ' Frequency in hours to create LAK : ' , lakfrq
    if ( ichem == 1 ) then
      write(stdout,'(a,f6.1)') ' Frequency in hours to create CHE : ' , chemfrq
      write(stdout,'(a,f6.1)') ' Frequency in hours to create OPT : ' , chemfrq
    end if

    write(stdout,*) 'Physical Parameterizations'
    write(stdout,'(a,i2)') '  Lateral Boundary conditions : ' , iboudy
    write(stdout,'(a,i2)') '  Semi-Lagrangian Advection   : ' , isladvec
    write(stdout,'(a,i2)') '  Land cumulus conv. scheme   : ' , icup_lnd
    write(stdout,'(a,i2)') '  Ocean cumulus conv. scheme  : ' , icup_ocn
    if ( any(icup == 2) ) then
      write(stdout,'(a,i2)') '  Grell closure scheme        : ' , igcc
    end if
    write(stdout,'(a,i2)') '  Moisture schem              : ' , ipptls
    write(stdout,'(a,i2)') '  Ocean Flux scheme           : ' , iocnflx
    if ( iocnflx == 2 ) then
      write(stdout,'(a,i2)') '  Zeng roughness formula      : ' , iocnrough
    end if
    write(stdout,'(a,i2)') '  Coupling with ROMS ocean    : ' , iocncpl
    write(stdout,'(a,i2)') '  Pressure gradient force     : ' , ipgf
    write(stdout,'(a,i2)') '  Prescribed LW emissivity    : ' , iemiss
#ifndef CLM
    write(stdout,'(a,i2)') '  Lake model in BATS          : ' , lakemod
    write(stdout,'(a,i2)') '  Simulate diurnal sst cycle  : ' , idcsst
    write(stdout,'(a,i2)') '  Simulate sea ice cover      : ' , iseaice
    write(stdout,'(a,i2)') '  Simulate desert seasons     : ' , idesseas
#endif
    write(stdout,'(a,i2)') '  Enable chem/aerosol model   : ' , ichem
    write(stdout,'(a,i2)') '  Convective LWP scheme       : ' , iconvlwp
    write(stdout,*) 'Boundary Pameterizations'
    write(stdout,'(a,f9.6)') '  Nudge value high range     : ', high_nudge
    write(stdout,'(a,f9.6)') '  Nudge value medium range   : ', medium_nudge
    write(stdout,'(a,f9.6)') '  Nudge value low range      : ', low_nudge
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
    write(stdout,'(a,f12.6)') '  time step for radiation '// &
          'model in minutes : ' , dtrad
    write(stdout,'(a,f12.6)') '  time step for emission  '// &
          'model in hours   : ' , dtabem
  end if

  if ( nsg > 1 ) then
    call read_subdomain_info(mdsub%ht,mdsub%lndcat,mdsub%mask, &
             mdsub%xlat,mdsub%xlon,mdsub%dhlake)
    mdsub%ht = mdsub%ht*egrav
  else
    do i = ici1 , ici2
      do j = jci1 , jci2
        mdsub%ht(1,j,i) = mddom%ht(j,i)*egrav
        mdsub%lndcat(1,j,i) = mddom%lndcat(j,i)
        mdsub%xlat(1,j,i) = mddom%xlat(j,i)
        mdsub%xlon(1,j,i) = mddom%xlon(j,i)
        mdsub%mask(1,j,i) = mddom%mask(j,i)
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
!------invert mapscale factors and convert hgt to geopotential
!
  do i = ide1 , ide2
    do j = jde1 , jde2
      mddom%ht(j,i)   = mddom%ht(j,i)*egrav
      mddom%msfd(j,i) = d_one/mddom%msfd(j,i)
      mddom%msfx(j,i) = d_one/mddom%msfx(j,i)
    end do
  end do

  call exchange(mddom%ht,1,jde1,jde2,ide1,ide2)
  call exchange(mddom%msfx,2,jde1,jde2,ide1,ide2)
  call exchange(mddom%msfd,2,jde1,jde2,ide1,ide2)
!

  call init_advection(mddom,sfs,atm1,qdot,kpbl)
  call init_precip(atms,atm2,aten,sfs,pptnc,cldfra,cldlwc)
  call init_surface_model
  call init_cumulus

  if ( ichem == 1 ) then
#if defined CLM45
    call init_chem(atms,mddom,sfs,xpsb,ba_cr,fcc,cldfra,rembc,remrat,    &
                   coszrs,svegfrac2d,sxlai2d,sfracv2d,sfracb2d,sfracs2d, &
                   solis,sdelt,sdelq,ssw2da,convpr,icumtop,              &
                   icumbot,taucldsp,lms)
#else
#if defined CLM
    call init_chem(atms,mddom,sfs,xpsb,ba_cr,fcc,cldfra,rembc,remrat,    &
                   coszrs,svegfrac2d,sxlai2d,sfracv2d,sfracb2d,sfracs2d, &
                   solis,sdelt,sdelq,ssw2da,convpr,icumtop,              &
                   icumbot,taucldsp,voc_em,voc_em1,voc_em2,dep_vels)
#else
    call init_chem(atms,mddom,sfs,xpsb,ba_cr,fcc,cldfra,rembc,remrat,    &
                   coszrs,svegfrac2d,sxlai2d,sfracv2d,sfracb2d,sfracs2d, &
                   solis,sdelt,sdelq,ssw2da,convpr,icumtop,icumbot,taucldsp)
#endif
#endif
    do n = 1 , ntr
      call bcast(chtrname(n),6)
    end do
  end if
  call init_radiation
  if ( islab_ocean == 1 ) then
    call allocate_mod_slabocean
    call init_slabocean(sfs,mddom%ldmsk,fsw,flw)
  end if
!
  if ( myid == italk ) then
    write(stdout,*) 'Domain grid parameters:'
    write(stdout,'(a,a)') '  Map Projection        : ',iproj
    write(stdout,'(a,i4,a,i4,a,i3)') &
      '  Dot Grid Full Extent  : ',jx,'x',iy,'x',kz
    write(stdout,'(a,f11.6,a)') '  Model Top Pressure    : ',ptop,' cb'
    write(stdout,'(a,f11.6,a)') '  Model Grid Spacing    : ',ds,' km'
    write(stdout,'(a,f11.6,a)') '  Grid Center Latitude  : ',clat,' deg'
    write(stdout,'(a,f11.6,a)') '  Grid Center longitude : ',clon,' deg'
    if ( iproj == 'ROTMER' ) then
      write(stdout,'(a,f11.6,a)') '  Pole Latitude         : ',plat,' deg'
      write(stdout,'(a,f11.6,a)') '  Pole longitude        : ',plon,' deg'
    else if ( iproj == 'LAMCON' ) then
      write(stdout,'(a,f11.6,a)') '  True Latitude 1       : ',truelatl,' deg'
      write(stdout,'(a,f11.6,a)') '  True Latitude 2       : ',truelath,' deg'
    end if
  end if
!
!-----compute land/water mask on subgrid space
!
   do i = ici1 , ici2
     do j = jci1 , jci2
       if ( mddom%mask(j,i) > 0.1D0 ) then
         mddom%ldmsk(j,i) = 1
       else
         mddom%ldmsk(j,i) = 0
       end if
     end do
   end do
   do i = ici1 , ici2
     do j = jci1 , jci2
       do n = 1 , nnsg
         if ( mdsub%mask(n,j,i) > 0.1D0 ) then
           mdsub%ldmsk(n,j,i) = 1
         else
           mdsub%ldmsk(n,j,i) = 0
         end if
       end do
     end do
   end do
!
!-----compute dsigma and half sigma levels.
!
  do k = 1 , kz
    dsigma(k) = (sigma(k+1) - sigma(k))
    hsigma(k) = (sigma(k+1) + sigma(k))*d_half
  end do
!
!----calculate max no of pbl levels: kmxpbl=k at highest allowed pbl level
!-----1. caluclate sigma level at 700mb, assuming 600mb highest
!-----   allowed pbl, and 1013mb as standard surface pressure. (sigtbl)
!-----2. find first model sigma level above sigtbl.
!
  sigtbl = (70.0D0-ptop)/(101.3D0-ptop)
  kmxpbl = kz
  do k = kz , 1 , -1
    delsig = hsigma(k) - sigtbl
    if ( delsig <= d_zero ) then
      kmxpbl = k
      exit
    end if
  end do
  if ( ipptls > 0 ) then
    cevaplnd = max(cevaplnd,d_zero)
    caccrlnd = max(caccrlnd,d_zero)
    cevapoce = max(cevapoce,d_zero)
    caccroce = max(caccroce,d_zero)
  end if

  if ( myid == italk ) then
    if ( ibltyp == 1 ) then
      write(stdout,*) 'PBL limit for Holtstag'
      write(stdout,'(a,i3)') '  Index of highest allowed pbl : ',kmxpbl
    end if
    if ( ipptls > 0 ) then
      write(stdout,*) 'SUBEX large scale precipitation parameters'
      write(stdout,'(a,i2)' )   '  # of bottom no cloud model levels : ',ncld
      write(stdout,'(a,f11.6,a,f11.6)')                      &
          '  Auto-conversion rate:         Land = ',qck1land,' Ocean = ',qck1oce
      write(stdout,'(a,f11.6,a,f11.6)')                      &
          '  Relative humidity thresholds: Land = ',rh0land, ' Ocean = ',rh0oce
      write(stdout,'(a,f11.6,a,f11.6)')                      &
          '  Gultepe factors:              Land = ',gulland ,' Ocean = ',guloce
      write(stdout,'(a,f11.6)') &
          '  Maximum relative humidity         : ' , rhmax
      write(stdout,'(a,f11.6)') &
          '  RH0 temperature threshold         : ' , tc0
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
          '  Maximum total cloud cover for rad : ', cftotmax
      write(stdout,'(a,f11.6)') &
          '  Condensation threshold            : ', conf
    end if
  end if

  call allocate_cumulus

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
        bb = dlog(hsigma(ktop)) + dlog(hsigma(kbase))
        cc = dlog(hsigma(ktop))*dlog(hsigma(kbase))
        ssum = d_zero
        do k = ktop , kbase
          xx = dlog(hsigma(k))
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
        xtop = dlog((d_100-ptop)*sigma(ktop)+ptop)
        xbot = dlog((d_100-ptop)*sigma(kbase+1)+ptop)
        bb = xtop + xbot
        cc = xtop*xbot
        vqmax = d_zero
        ssum = d_zero
        xx = xtop
        yy = xbot
        wk = (xx*xx) - (bb*xx) + cc
        qk = -((yy*yy)-(bb*yy)+cc)
        do k = ktop , kbase
          xx = dlog((d_100-ptop)*sigma(k+1)+ptop)
          yy = dlog((d_100-ptop) * &
               (sigma(ktop)+sigma(kbase+1)-sigma(k+1))+ptop)
          wkp1 = (xx*xx) - (bb*xx) + cc
          qkp1 = -((yy*yy)-(bb*yy)+cc)
          vqflx(k,kbase,ktop) = -((wkp1*qkp1)-(wk*qk))/dsigma(k)
          ssum = ssum + vqflx(k,kbase,ktop)
          if ( dabs(vqflx(k,kbase,ktop)) > vqmax ) &
               vqmax = dabs(vqflx(k,kbase,ktop))
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
      end if
    end if

    do i = ici1 , ici2
      do j = jci1 , jci2
        if ( mddom%lndcat(j,i) > 14.5D0 .and. &
             mddom%lndcat(j,i) < 15.5D0) then
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
      write(stderr,*) 'WARNING : The Betts-Miller Convection scheme is not ', &
                      'properly implemented'
    end if
  end if
  if ( any(icup == 4) ) then
    minorig = kz
    do k = 1 , kz
      if ( hsigma(k) <= minsig ) minorig = kz - k
    end do
    if ( myid == italk ) then
      write(stdout,*) 'Emanuel (1991) Convection Scheme V4.3C used.'
      write(stdout,'(a,f11.6)') '  Min Convection origin             : ',minsig
      write(stdout,'(a,f11.6)') '  Autoconversion Threshold (ocean)  : ', &
        elcrit_ocn
      write(stdout,'(a,f11.6)') '  Autoconversion Threshold (land)   : ', &
        elcrit_lnd
      write(stdout,'(a,f11.6)') '  Autoconversion Threshold to zero  : ',tlcrit
      write(stdout,'(a,f11.6)') '  Entrainment Coefficient           : ',entp
      write(stdout,'(a,f11.6)') '  Fractional area of uns. downdraft : ',sigd
      write(stdout,'(a,f11.6)') '  Fractional area of uns. downdraft : ',sigd
      write(stdout,'(a,f11.6)') '  Fall speed of rain                : ',omtrain
      write(stdout,'(a,f11.6)') '  Fall speed of snow                : ',omtsnow
      write(stdout,'(a,f11.6)') '  Rain evaporation coefficient      : ',coeffr
      write(stdout,'(a,f11.6)') '  Snow evaporation coefficient      : ',coeffs
      write(stdout,'(a,f11.6)') '  Convective momentum transport coef: ',cu
      write(stdout,'(a,f11.6)') '  Downdraft velocity scale          : ',betae
      write(stdout,'(a,f11.6)') '  Max negative perturbation blw LFC : ',dtmax
      write(stdout,'(a,f11.6)') '  Quasi-equilibrium approach rate 1 : ',alphae
      write(stdout,'(a,f11.6)') '  Quasi-equilibrium approach rate 2 : ',damp
      write(stdout,'(a,f11.6)') '  Precipitation efficienct (land)   : ', &
        epmax_lnd
      write(stdout,'(a,f11.6)') '  Precipitation efficienct (ocean)  : ', &
        epmax_ocn
    end if
  end if
  if ( any(icup == 5)  ) then
    if ( myid == italk .and. iconv /= 4 ) then
      write(stdout,*) 'Tiedtke (1986) Convection Scheme ECHAM 5.4 used.'
      write(stdout,'(a,i2)')    '  Used Scheme                       : ',iconv
      write(stdout,'(a,f11.6)') '  Entrainment rate penetrative conv.: ',entrpen
      write(stdout,'(a,f11.6)') '  Entrainment rate shallow conv.    : ',entrscv
      write(stdout,'(a,f11.6)') '  Entrainment rate midlev conv.     : ',entrmid
      write(stdout,'(a,f11.6)') '  Entrainment rate cumulus downdraft: ',entrdd
      write(stdout,'(a,f11.6)') '  CLW to rain conversion factor     : ',cprcon
      write(stdout,'(a,f11.6)') '  Max entrainment                   : ',entrmax
    end if
    if ( myid == italk .and. iconv == 4 ) then
      write(stdout,*) 'Tiedtke (1986) Convection Scheme ECMWF 38R2 used.'
    end if
    cevapu = cevaplnd
  end if
  if ( any(icup == 6)  ) then
    if ( myid == italk ) then
      write(stdout,*) 'Kain Fritsch scheme used.'
      write(stdout,'(a,f11.6)') '  Entrainment rate : ', kf_entrate
    end if
  end if

  if ( ipptls == 2 ) then
    call init_cloud_s1(atms,aten,heatrt,sfs,q_detr,pptnc)
  end if

  call init_pblscheme

  !
  ! Convective Cloud Cover
  !
  afracl = 0.25D0    ! frac. cover for conv. precip. when dx=dxlarg
  afracs = clfrcvmax !   "     "    "    "      "     "   dx=dxsmal
  dlargc = 100.0D0
  dsmalc = 10.0D0
  dxtemc = dmin1(dmax1(ds,dsmalc),dlargc)
  clfrcv = afracl + (afracs-afracl)*((dlargc-dxtemc)/(dlargc-dsmalc))**2
  clfrcv = dmin1(clfrcv,d_one)
  if ( myid == italk ) then
    write(stdout,*) 'Convective Cloud Cover parameters after resolution scaling'
    write(stdout,'(a,f11.6)') '  Maximum Convective Cloud Cover : ',clfrcv
    write(stdout,'(a,f11.6)') '  Convective Cloud Water         : ',cllwcv
  end if
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

  chibot = 450.0D0
  ptmb = d_10*ptop
  pz = hsigma(1)*(d_1000-ptmb) + ptmb
  if ( pz > chibot ) call fatal(__FILE__,__LINE__,                 &
                                 'VERTICAL INTERPOLATION ERROR')
  do k = 1 , kz
    pk = hsigma(k)*(d_1000-ptmb) + ptmb
    if ( pk <= chibot ) kchi = k
  end do

  if ( ipptls > 0 ) then
    do i = ici1 , ici2
      do j = jci1 , jci2
        if ( mddom%lndcat(j,i) > 14.5D0 .and. &
             mddom%lndcat(j,i) < 15.5D0) then
          qck1(j,i)  = qck1oce  ! OCEAN
          cgul(j,i)  = guloce
          rh0(j,i)   = rh0oce
          cevap(j,i) = cevapoce
          caccr(j,i) = caccroce
        else
          qck1(j,i)  = qck1land ! LAND
          cgul(j,i)  = gulland
          rh0(j,i)   = rh0land
          cevap(j,i) = cevaplnd
          caccr(j,i) = caccrlnd
        end if
      end do
    end do
  end if
  !
  ! Setup Holtslag PBL Critical Richardson Number
  !
  if ( ibltyp == 1 .or. ibltyp == 99 ) then
    do i = ici1 , ici2
      do j = jci1 , jci2
        if ( mddom%ldmsk(j,i) == 1 ) then
          ricr(j,i) = ricr_lnd
        else
          ricr(j,i) = ricr_ocn
        end if
      end do
    end do
  end if
  !
  ! Setup Land/Ocean MIT autoconversion
  !
  if ( any(icup == 4) ) then
    do i = ici1 , ici2
      do j = jci1 , jci2
        if ( mddom%ldmsk(j,i) == 1 ) then
          elcrit2d(j,i) = elcrit_lnd
          epmax2d(j,i) = epmax_lnd
        else
          elcrit2d(j,i) = elcrit_ocn
          epmax2d(j,i) = epmax_ocn
        end if
      end do
    end do
  end if
  if ( any(icup == 1) ) then
    sig700 = (70.0D0-ptop)/(d_100-ptop)
    do k = 1 , kz
      k700 = k
      if ( sig700 <= sigma(k+1) .and. sig700 > sigma(k) ) exit
    end do
  end if
  !
  ! Setup Boundary condition routines.
  !
  call setup_bdycon(hsigma)
  if ( ichem == 1 ) call setup_che_bdycon
!
  if ( any(icup == 3) ) call lutbl(ptop)
!
  if ( iboudy < 0 .or. iboudy > 5 ) then
    call fatal(__FILE__,__LINE__,'UNSUPPORTED BDY SCHEME.')
  end if

  if ( myid == italk ) then
    if ( ibltyp == 0 ) &
      write(stdout,*) 'Model frictionless and insulated for the lower boundary.'

    write(stdout,*) &
      'The surface energy budget is used to calculate the ground temperature.'
    write(stdout,'(a,i4,a)') &
      ' The radiation is computed every ',ntrad,' time steps.'
!
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
    end if

    write(stdout,'(a,7x,a,11x,a,6x,a,7x,a,7x,a,9x,a)') '# k','sigma','a',&
      'dsigma','twt(1)','twt(2)','qcon'
!
    do k = 1 , kz
      write(stdout,'(1x,i2,5x,f7.4,5x,f7.4,5x,f7.4,5x,f8.4,5x,f8.4,5x,f8.4)') &
        k , sigma(k) , hsigma(k) , dsigma(k) , twt(k,1) , twt(k,2) , qcon(k)
    end do
    write(stdout,'(1x,i2,5x,f7.4)') kzp1 , sigma(kzp1)
    write(stdout,'(a,e13.6,a)') &
      '  Constant hor. diff. coef. = ',xkhz,' m^2 s-1'
    write(stdout,'(a,e13.6,a)') &
      '  Maximumt hor. diff. coef. = ',xkhmax,' m^2 s-1'
  end if

  if (itweak == 1 ) then
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
104 call fatal(__FILE__,__LINE__, 'Error reading SUBEXPARAM')
105 call fatal(__FILE__,__LINE__, 'Error reading MICROPARAM')
106 call fatal(__FILE__,__LINE__, 'Error reading GRELLPARAM')
107 call fatal(__FILE__,__LINE__, 'Error reading EMANPARAM')
108 call fatal(__FILE__,__LINE__, 'Error reading TIEDTKEPARAM')
109 call fatal(__FILE__,__LINE__, 'Error reading KFPARAM')
110 call fatal(__FILE__,__LINE__, 'Error reading HOLTSLAGPARAM')
111 call fatal(__FILE__,__LINE__, 'Error reading UWPARAM')
112 call fatal(__FILE__,__LINE__, 'Error reading RRTMPARAM')
113 call fatal(__FILE__,__LINE__, 'Error reading SLABOCPARAM')
114 call fatal(__FILE__,__LINE__, 'Error reading CHEMPARAM')
#ifdef CLM
115 call fatal(__FILE__,__LINE__, 'Error reading CLMPARAM')
#endif
116 call fatal(__FILE__,__LINE__, 'Error reading CPLPARAM')
117 call fatal(__FILE__,__LINE__, 'Error reading TWEAKPARAM')

  end subroutine param
!
end module mod_params

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
  use mod_cu_interface
  use mod_lm_interface
  use mod_atm_interface
  use mod_che_interface
  use mod_rad_interface
  use mod_pbl_interface
  use mod_precip
  use mod_cloud_s1
  use mod_split
  use mod_slice
  use mod_bdycod
  use mod_ncio
  use mod_tendency
  use mod_ncio
  use mod_advection , only : init_advection
  use mod_mppio
#ifdef CLM
  use mod_clm
  use clm_varsur , only : landmask
#endif

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
  include 'mpif.h'
!
  real(dp) :: afracl , afracs , bb , cc , chibot , delsig , &
             dlargc , dsmalc , dxtemc , pk , ptmb , pz , qk ,       &
             qkp1 , sig700 , sigtbl , ssum , vqmax , vqrang , wk ,  &
             wkp1 , xbot , xtop , xx , yy
  real(dp) :: gulland , guloce , mincld , qck1land , qck1oce , &
             rh0land , rh0oce , skbmax , clfrcvmax
  real(dp) :: shrmax_ocn , shrmin_ocn , edtmax_ocn , edtmin_ocn , &
             edtmaxo_ocn , edtmino_ocn , edtmaxx_ocn , edtminx_ocn
  real(dp) :: shrmax , shrmin , edtmax , edtmin , edtmaxo , &
             edtmino , edtmaxx , edtminx , dtauc , pbcmax , htmax , htmin
  integer :: kbmax
  real(dp) , dimension(nsplit) :: dtsplit
  integer :: i , j , k , kbase , ktop , ns , mdate0 , mdate1 , mdate2
  integer :: hspan
  logical , dimension(n_atmvar) :: atm_enablevar
  logical , dimension(n_srfvar) :: srf_enablevar
  logical , dimension(n_stsvar) :: sts_enablevar
  logical , dimension(n_subvar) :: sub_enablevar
  logical , dimension(n_radvar) :: rad_enablevar
  logical , dimension(n_lakvar) :: lak_enablevar
  integer(8) :: ndbgfrq , nsavfrq , natmfrq , nradfrq , nchefrq , nsrffrq
  integer(8) :: nbdyfrq
  integer :: n , len_path
  integer :: ierr
  character(len=32) :: appdat
  type(rcm_time_interval) :: bdif
#ifndef CLM
  real(dp) :: clmfrq
#endif
  character (len=64) :: subroutine_name='param'
  integer :: idindx=0
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
 
  namelist /outparam/ ifsave , savfrq , ifatm , atmfrq , ifrad ,   &
    radfrq , ifsrf , ifsub , iflak , ifsts , srffrq , lakfrq ,    &
    ifchem , chemfrq , atm_enablevar , srf_enablevar ,             &
    rad_enablevar , sub_enablevar , sts_enablevar ,                &
    lak_enablevar , dirout

  namelist /physicsparam/ ibltyp , iboudy , icup , igcc , ipgf ,    &
    iemiss , lakemod , ipptls , iocnflx , iocncpl , iocnrough ,     &
    ichem , scenario , idcsst , iseaice , idesseas , iconvlwp ,     &
    irrtm , iclimao3

  namelist /rrtmparam/ inflgsw , iceflgsw , liqflgsw , inflglw ,    &
    iceflglw , liqflglw , icld , irng , idrv

  namelist /subexparam/ ncld , fcmax , qck1land , qck1oce ,         &
    gulland , guloce , rhmax , rh0oce , rh0land , cevap , caccr ,   &
    tc0 , cllwcv , clfrcvmax , cftotmax

  namelist /grellparam/ shrmin , shrmax , edtmin , edtmax ,         &
    edtmino , edtmaxo , edtminx , edtmaxx , pbcmax , mincld ,       &
    htmin , htmax , skbmax , dtauc, shrmin_ocn , shrmax_ocn ,       &
    edtmin_ocn, edtmax_ocn, edtmino_ocn , edtmaxo_ocn ,             &
    edtminx_ocn , edtmaxx_ocn 
 
  namelist /emanparam/ minsig , elcrit , tlcrit , entp , sigd ,     &
    sigs , omtrain , omtsnow , coeffr , coeffs , cu , betae ,       &
    dtmax , alphae , damp
 
  namelist /tiedtkeparam/ iconv , entrpen , entrscv , entrmid ,  &
    entrdd , cmfcmax , cmfcmin , cmfdeps , cmfctop , rhcdd ,     &
    cmtcape , zdlev , cprcon , ctrigger , nmctop , lmfpen ,      &
    lmfscv , lmfmid , lmfdd , lmfdudv

  namelist /chemparam/ chemsimtype , ichremlsc , ichremcvc , ichdrdepo , &
         ichcumtra , ichsolver , idirect , ichdustemd , ichdiag

  namelist /uwparam/ iuwvadv , ilenparam , atwo , rstbl

#ifdef CLM
  namelist /clmparam/ dirclm , imask , clmfrq
#endif

  namelist /cplparam/ cpldt, cpldbglevel
!
  call time_begin(subroutine_name,idindx)
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
!     icup   : type of cumulus parameterization
!     = 1 ; kuo
!     = 2 ; grell
!     = 3 ; betts-miller (1986)
!     = 4 ; emanuel (1991)
!     = 5 ; Tiedtke (1986) - version from ECHAM 5.4 
!     = 98; variable: emanuel over land and grell over ocean
!     = 99; variable: grell over land and emanuel over ocean
!
!     igcc   : Grell Scheme Convective Closure Assumption
!     = 1 ; Arakawa & Schubert (1974)
!     = 2 ; Fritsch & Chappell (1980)
!
!     ipptls : type of moisture scheme
!     = 1 ; explicit moisture (SUBEX; Pal et al 2000)
!
!     iocnflx: type of ocean flux parameterization
!     = 1 ; BATS
!     = 2 ; Zeng et al.
!
!     iocncpl: controls the coupling with ROMS
!     = 0 ; no coupling
!     = 1 ; activate ROMS coupling
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
  idate0 = 1993063000  ! Start date of the initial simulation
  idate1 = 1993063000  ! Start date of this simulation
  idate2 = 1993080100  ! End Date this simulation
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
  savfrq = 6.0D0
  ifatm = .true.
  atmfrq = 6.0D0
  ifrad = .true.
  radfrq = 6.0D0     ! time interval for disposing rad output (hrs)
  ifsrf = .true.
  ifsts = .true.
  ifsub = .false.
  iflak = .false.
  srffrq = 1.0D0    ! time interval for disposing bats output (hrs)
  lakfrq = -1.0D0   ! time interval for disposing lake output (hrs)
  atm_enablevar(:) = .true.
  srf_enablevar(:) = .true.
  sts_enablevar(:) = .true.
  sub_enablevar(:) = .true.
  lak_enablevar(:) = .true.
  rad_enablevar(:) = .true.
  dirout = './output' 
!chem2
  ifchem = .false.
  chemfrq = 6.0D0   ! time interval for disposeing chem output (hrs)
!chem2_
  clmfrq = 12.0D0
!
!----------------------------------------------------------------------
!-----namelist physicsparam:
!
  ibltyp = 1
  iboudy = 1
  icup = 1
  ipptls = 1
  igcc = 1
  ipgf = 1
  iemiss = 1
  iocnflx = 1
  iocncpl = 0
  iocnrough = 1
  lakemod = 0
  ichem = 0
  scenario = 'A1B'
  idcsst = 0
  iseaice = 0
  idesseas = 1
  iconvlwp = 1
  irrtm = 0
  iclimao3 = 0
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
  fcmax = 0.80D0       ! Maximum cloud fraction cover (rad only)
! qck1land = 0.0005D0  ! Autoconversion Rate for Land
! qck1oce  = 0.0005D0  ! Autoconversion Rate for Ocean
  qck1land = 0.00025D0 ! Autoconversion Rate for Land
  qck1oce = 0.00025D0  ! Autoconversion Rate for Ocean
  gulland = 0.4D0      ! Fract of Gultepe eqn (qcth) when prcp occurs (land)
  guloce = 0.4D0       ! Fract of Gultepe eqn (qcth) for ocean
  rhmax = 1.01D0       ! RH at whicn FCC = 1.0
  rh0oce = 0.90D0      ! Relative humidity threshold for ocean
  rh0land = 0.80D0     ! Relative humidity threshold for land
  tc0 = 238.0D0        ! Below this temp, rh0 begins to approach unity
! cevap    = 0.2D-4    ! Raindrop evap rate coef [[(kg m-2 s-1)-1/2]/s]
  cevap = 1.0D-3     ! Raindrop evap rate coef [[(kg m-2 s-1)-1/2]/s]
! caccr    = 6.0D0   ! Raindrop accretion rate [m3/kg/s]
  caccr = 3.0D0      ! Raindrop accretion rate [m3/kg/s]
  cllwcv = 0.5D-4    ! Cloud liquid water content for convective precip.
  clfrcvmax = 0.25D0 ! Max cloud fractional cover for convective precip.
  cftotmax = 0.75D0  ! Max total cover cloud fraction for radiation
 
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
  dtauc = 30.0D0         ! Fritsch & Chappell (1980) 
! 
!------namelist emanparam:
  minsig = 0.95D0   ! Lowest sigma level from which convection can originate
  elcrit = 0.0011D0 ! Autoconversion threshold water content (gm/gm)
  tlcrit = -55.0D0  ! Below tlcrit auto-conversion threshold is zero
  entp = 1.5D0      ! Coefficient of mixing in the entrainment formulation
  sigd = 0.05D0     ! Fractional area covered by unsaturated dndraft
  sigs = 0.12D0     ! Fraction of precipitation falling outside of cloud
  omtrain = 50.0D0  ! Fall speed of rain (P/s)
  omtsnow = 5.5D0   ! Fall speed of snow (P/s)
  coeffr = 1.0D0    ! Coefficient governing the rate of rain evaporation
  coeffs = 0.8D0    ! Coefficient governing the rate of snow evaporation
  cu = 0.7D0        ! Coefficient governing convective momentum transport
  betae = 10.0D0    ! Controls downdraft velocity scale
  dtmax = 0.9D0     ! Max negative parcel temperature perturbation below LFC
  alphae = 0.2D0    ! Controls the approach rate to quasi-equilibrium
  damp = 0.1D0      ! Controls the approach rate to quasi-equilibrium
!
!------namelist tiedtkeparam:
  iconv    = 1  ! Selects the actual scheme
  entrpen  = 1.0D-4
  entrscv  = 3.0D-4
  entrmid  = 1.0D-4
  entrdd   = 2.0D-4
  cmfcmax  = 1.0D0
  cmfcmin  = 1.0D-10
  cmfdeps  = 0.3D0
  rhcdd    = 1.0D0
  cmtcape  = 40.0D0
  zdlev    = 1.5D4
! THOSE ARE FUNCTION OF GRID AND VERTICAL RESOLUTION
  nmctop   = 4
  cmfctop  = 0.35D0
  cprcon   = 1.0D-4
  ctrigger = -1.1D-0
! Control switch flags
  lmfpen   = .true.
  lmfscv   = .true.
  lmfmid   = .true.
  lmfdd    = .true.
  lmfdudv  = .true.
!
!c------namelist uwparam ;
  iuwvadv = 0
  ilenparam = 0
  atwo = 15.0D0
  rstbl = 1.5D0

!c------namelist chemparam ; ( 0= none, 1= activated)
  ichsolver = 1     ! enable chem solver
  ichremlsc = 1     ! tracer removal by large scale clouds
  ichremcvc = 1     ! tracer removal by convective clouds
  ichdrdepo = 1     ! tracer dry deposition
  ichcumtra = 1     ! tracer convective transport
  ichdustemd = 1    ! dust emission distribution (1 = alfaro, 2 =kok)
  idirect = 1       ! tracer direct effect
  ichdiag = 0       ! chem tend outputs 
#ifdef CLM
!c------CLM Specific
  imask = 1
#endif
!------namelist cplparam ;
! cpldbglevel:
! 0 = no debugging
! 1 = only informative print
! 2 = previous + write grid information in VTK format
! 3 = previous + write exchange fields into NetCDF
! 4 = previous + write exchange fileds into ASCII
!
  cpldt = 21600.0D0       ! coupling time step in seconds (seconds)
  cpldbglevel = 1         ! debugging level
!
!---------------------------------------------------------------------
!
#ifdef CLM
  if ( myid == italk ) then
    if (nsg /= 1 ) then
      write (6,*) 'Running SUBGRID with CLM: not implemented'
      write (6,*) 'Please set nsg to 1 in regcm.in'
      call fatal(__FILE__,__LINE__,'CLM & SUBGRID TOGETHER')
    end if
  end if
#endif

  if ( myid == iocpu ) then
!  
!-----read in namelist variables:
    read (ipunit, restartparam)
    print * , 'RESTARTPARAM namelist READ IN'

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

    read (ipunit, timeparam)
    print * , 'TIMEPARAM namelist READ IN'
    read (ipunit, outparam)
    print * , 'OUTPARAM namelist READ IN'
    len_path = len(trim(dirout))
    if ( dirout(len_path:len_path) /= '/' ) dirout = trim(dirout)//'/'
    if ( lakfrq < d_zero ) lakfrq = srffrq
    read (ipunit, physicsparam)
    print * , 'PHYSICSPARAM namelist READ IN'
    if ( ipptls == 1 ) then
      read (ipunit, subexparam)
      print * , 'SUBEXPARAM namelist READ IN'
    end if
    if ( icup == 2 .or. icup == 99 .or. icup == 98 ) then
      read (ipunit, grellparam)
      print * , 'GRELLPARAM namelist READ IN'
    end if
    if ( icup == 4 .or. icup == 99 .or. icup == 98 ) then
      read (ipunit, emanparam)
      print * , 'EMANPARAM namelist READ IN'
    end if
    if ( icup == 5 ) then
      read (ipunit, tiedtkeparam)
      print * , 'TIEDTKEPARAM namelist READ IN'
    end if
    if ( icup < 0 .or. (icup > 5 .and. icup < 98) .or. icup > 99 ) then
      call fatal(__FILE__,__LINE__,'UNSUPPORTED CUMULUS SCHEME')
    end if
    if ( ibltyp == 2 .or. ibltyp == 99 ) then
      read (ipunit, uwparam)
      print * , 'UWPARAM namelist READ IN'
    end if
    if ( irrtm == 1 ) then
      read (ipunit, rrtmparam)
      print * , 'RRTMPARAM namelist READ IN'
    end if
    if ( ichem == 1 ) then
      read (ipunit, chemparam)
      print * , 'CHEMPARAM namelist READ IN'
    else
      ichem = 0
      ntr = 0
      ifchem = .false.
    end if
#ifdef CLM
    read (ipunit , clmparam)
    print * , 'CLMPARAM namelist READ IN'
#endif
    if (iocncpl == 1) then
      read (ipunit , cplparam)
      print *, 'CPLPARAM namelist READ IN'
    end if
  end if 

  call mpi_barrier(mycomm,ierr) 
!
!  communicate to all processors 
!
  call mpi_bcast(ifrest,1,mpi_logical,iocpu,mycomm,ierr)
  call mpi_bcast(hspan,1,mpi_integer,iocpu,mycomm,ierr)
  call date_bcast(idate0,iocpu,mycomm,ierr)
  call date_bcast(idate1,iocpu,mycomm,ierr)
  call date_bcast(idate2,iocpu,mycomm,ierr)
 
  call mpi_bcast(dtrad,1,mpi_real8,iocpu,mycomm,ierr)
  call mpi_bcast(dtabem,1,mpi_real8,iocpu,mycomm,ierr)
  call mpi_bcast(dtsrf,1,mpi_real8,iocpu,mycomm,ierr)
  call mpi_bcast(dt,1,mpi_real8,iocpu,mycomm,ierr)
 
  call mpi_bcast(ifsave,1,mpi_logical,iocpu,mycomm,ierr)
  call mpi_bcast(savfrq,1,mpi_real8,iocpu,mycomm,ierr)
  call mpi_bcast(ifatm,1,mpi_logical,iocpu,mycomm,ierr)
  call mpi_bcast(atmfrq,1,mpi_real8,iocpu,mycomm,ierr)
  call mpi_bcast(ifrad,1,mpi_logical,iocpu,mycomm,ierr)
  call mpi_bcast(radfrq,1,mpi_real8,iocpu,mycomm,ierr)
  call mpi_bcast(ifsrf,1,mpi_logical,iocpu,mycomm,ierr)
  call mpi_bcast(ifsub,1,mpi_logical,iocpu,mycomm,ierr)
  call mpi_bcast(iflak,1,mpi_logical,iocpu,mycomm,ierr)
  call mpi_bcast(ifsts,1,mpi_logical,iocpu,mycomm,ierr)
  call mpi_bcast(srffrq,1,mpi_real8,iocpu,mycomm,ierr)
  call mpi_bcast(lakfrq,1,mpi_real8,iocpu,mycomm,ierr)
  call mpi_bcast(ifchem,1,mpi_logical,iocpu,mycomm,ierr)
  call mpi_bcast(chemfrq,1,mpi_real8,iocpu,mycomm,ierr)
 
  call mpi_bcast(iboudy,1,mpi_integer,iocpu,mycomm,ierr)
  call mpi_bcast(ibltyp,1,mpi_integer,iocpu,mycomm,ierr)
  call mpi_bcast(icup,1,mpi_integer,iocpu,mycomm,ierr)
  call mpi_bcast(igcc,1,mpi_integer,iocpu,mycomm,ierr)
  call mpi_bcast(ipptls,1,mpi_integer,iocpu,mycomm,ierr)
  call mpi_bcast(iocnflx,1,mpi_integer,iocpu,mycomm,ierr)
  call mpi_bcast(iocncpl,1,mpi_integer,iocpu,mycomm,ierr)
  call mpi_bcast(iocnrough,1,mpi_integer,iocpu,mycomm,ierr)
  call mpi_bcast(ipgf,1,mpi_integer,iocpu,mycomm,ierr)
  call mpi_bcast(iemiss,1,mpi_integer,iocpu,mycomm,ierr)
  call mpi_bcast(lakemod,1,mpi_integer,iocpu,mycomm,ierr)
  call mpi_bcast(ichem,1,mpi_integer,iocpu,mycomm,ierr)
  call mpi_bcast(ntr,1,mpi_integer,iocpu,mycomm,ierr)

  ! Force the correct scenario from dattyp in CMIP5
  if ( myid == iocpu ) then
    if ( dattyp(4:5) == '26' ) then
      if ( scenario /= 'RCP3PD' .or. scenario /= 'RCP2.6' ) then
        print *, 'Forcing scenario from dattyp to RCP2.6'
        scenario = 'RCP3PD'
      end if
    else if ( dattyp(4:5) == '45' ) then
      if ( scenario /= 'RCP4.5') then
        print *, 'Forcing scenario from dattyp to RCP4.5'
        scenario = 'RCP4.5'
      end if
    else if ( dattyp(4:5) == '60') then
      if ( scenario /= 'RCP6' .or.  scenario /= 'RCP60' .or.  &
           scenario /= 'RCP60') then
        print *, 'Forcing scenario from dattyp to RCP6.0'
        scenario = 'RCP6'
      end if
    else if ( dattyp(4:5) == '85') then
      if ( scenario /= 'RCP8.5') then
        print *, 'Forcing scenario from dattyp to RCP8.5'
        scenario = 'RCP8.5'
      end if
    end if
  end if
  call mpi_bcast(scenario,8,mpi_character,iocpu,mycomm,ierr)
  call mpi_bcast(idcsst,1,mpi_integer,iocpu,mycomm,ierr)
  call mpi_bcast(iseaice,1,mpi_integer,iocpu,mycomm,ierr)
  call mpi_bcast(idesseas,1,mpi_integer,iocpu,mycomm,ierr)
  call mpi_bcast(iconvlwp,1,mpi_integer,iocpu,mycomm,ierr)
  call mpi_bcast(irrtm,1,mpi_integer,iocpu,mycomm,ierr)
  call mpi_bcast(iclimao3,1,mpi_integer,iocpu,mycomm,ierr)

#ifdef CLM
  call mpi_bcast(dirclm,256,mpi_character,iocpu,mycomm,ierr)
  call mpi_bcast(imask,1,mpi_integer,iocpu,mycomm,ierr)
  call mpi_bcast(clmfrq,1,mpi_real8,iocpu,mycomm,ierr)
#endif

  if (iocncpl == 1) then
    call mpi_bcast(cpldt,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(cpldbglevel,1,mpi_integer,iocpu,mycomm,ierr)
  end if

  if ( ipptls == 1 ) then
    call mpi_bcast(ncld,1,mpi_integer,iocpu,mycomm,ierr)
    call mpi_bcast(fcmax,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(qck1land,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(qck1oce,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(gulland,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(guloce,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(rhmax,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(rh0oce,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(rh0land,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(tc0,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(cevap,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(caccr,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(cftotmax,1,mpi_real8,iocpu,mycomm,ierr)
  end if

  if ( irrtm == 1 ) then
    call mpi_bcast(inflgsw,1,mpi_integer,iocpu,mycomm,ierr)
    call mpi_bcast(iceflgsw,1,mpi_integer,iocpu,mycomm,ierr)
    call mpi_bcast(liqflgsw,1,mpi_integer,iocpu,mycomm,ierr)
    call mpi_bcast(inflglw,1,mpi_integer,iocpu,mycomm,ierr)
    call mpi_bcast(iceflglw,1,mpi_integer,iocpu,mycomm,ierr)
    call mpi_bcast(liqflglw,1,mpi_integer,iocpu,mycomm,ierr)
    call mpi_bcast(idrv,1,mpi_integer,iocpu,mycomm,ierr)
    call mpi_bcast(icld,1,mpi_integer,iocpu,mycomm,ierr)
    call mpi_bcast(irng,1,mpi_integer,iocpu,mycomm,ierr)
  end if

  call mpi_bcast(clfrcvmax,1,mpi_real8,iocpu,mycomm,ierr)
  call mpi_bcast(cllwcv,1,mpi_real8,iocpu,mycomm,ierr)
 
  if ( icup == 2 .or. icup == 99 .or. icup == 98 ) then
    call mpi_bcast(shrmin,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(shrmax,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(edtmin,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(edtmax,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(edtmino,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(edtmaxo,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(edtminx,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(edtmaxx,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(shrmin_ocn,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(shrmax_ocn,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(edtmin_ocn,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(edtmax_ocn,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(edtmino_ocn,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(edtmaxo_ocn,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(edtminx_ocn,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(edtmaxx_ocn,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(pbcmax,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(mincld,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(htmin,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(htmax,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(skbmax,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(dtauc,1,mpi_real8,iocpu,mycomm,ierr)
  end if
 
  if ( icup == 4 .or. icup == 99 .or. icup == 98 ) then
    call mpi_bcast(minsig,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(elcrit,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(tlcrit,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(entp,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(sigd,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(sigs,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(omtrain,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(omtsnow,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(coeffr,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(coeffs,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(cu,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(betae,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(dtmax,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(alphae,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(damp,1,mpi_real8,iocpu,mycomm,ierr)
  end if
 
  if ( icup.eq.5 ) then
    call init_convect_tables
    call mpi_bcast(iconv,1,mpi_integer,iocpu,mycomm,ierr)
    call mpi_bcast(entrpen,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(entrscv,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(entrmid,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(entrdd,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(cmfctop,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(cmfcmax,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(cmfcmin,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(cmfdeps,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(rhcdd,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(cmtcape,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(zdlev,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(cprcon,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(ctrigger,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(nmctop,1,mpi_integer,iocpu,mycomm,ierr)
    call mpi_bcast(lmfpen,1,mpi_logical,iocpu,mycomm,ierr)
    call mpi_bcast(lmfscv,1,mpi_logical,iocpu,mycomm,ierr)
    call mpi_bcast(lmfmid,1,mpi_logical,iocpu,mycomm,ierr)
    call mpi_bcast(lmfdd,1,mpi_logical,iocpu,mycomm,ierr)
    call mpi_bcast(lmfdudv,1,mpi_logical,iocpu,mycomm,ierr)
  end if

  if ( ibltyp == 2 .or. ibltyp == 99 ) then
    call mpi_bcast(iuwvadv,1,mpi_integer,iocpu,mycomm,ierr)
    call mpi_bcast(ilenparam,1,mpi_integer,iocpu,mycomm,ierr)
    call mpi_bcast(atwo,1,mpi_real8,iocpu,mycomm,ierr)
    call mpi_bcast(rstbl,1,mpi_real8,iocpu,mycomm,ierr)
  end if

  if ( ichem == 1 ) then
    call mpi_bcast(chemsimtype,1,mpi_integer,iocpu,mycomm,ierr)
    call mpi_bcast(ichremlsc,1,mpi_integer,iocpu,mycomm,ierr)
    call mpi_bcast(ichremcvc,1,mpi_integer,iocpu,mycomm,ierr)
    call mpi_bcast(ichdrdepo,1,mpi_integer,iocpu,mycomm,ierr)
    call mpi_bcast(ichcumtra,1,mpi_integer,iocpu,mycomm,ierr)
    call mpi_bcast(idirect,1,mpi_integer,iocpu,mycomm,ierr)
    call mpi_bcast(ichsolver,1,mpi_integer,iocpu,mycomm,ierr)
    call mpi_bcast(ichdustemd,1,mpi_integer,iocpu,mycomm,ierr)
    call mpi_bcast(ichdiag,1,mpi_integer,iocpu,mycomm,ierr)
    call chem_config
    call mpi_bcast(ntr,1,mpi_integer,iocpu,mycomm,ierr)
    call mpi_bcast(iaerosol,1,mpi_integer,iocpu,mycomm,ierr)  
    call mpi_bcast(ioxclim,1,mpi_integer,iocpu,mycomm,ierr)
    call mpi_bcast(igaschem,1,mpi_integer,iocpu,mycomm,ierr)  
  end if
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------ALLOCATE NEEDED SPACE---------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!

  call allocate_mod_runparams

  call allocate_mod_atm_interface(ibltyp)

  call allocate_mod_tend

  call allocate_mod_bdycon(iboudy)

  call allocate_mod_pbl_common(ichem)

  call allocate_mod_cu_common(ichem)

  call allocate_mod_precip(ichem)

  call allocate_mod_split

  call allocate_mod_mppio

  call allocate_mod_bats_common(ichem,idcsst,lakemod)
#ifndef CLM
  call allocate_mod_bats_mppio(lakemod)
#else
  call allocate_mod_clm(ntr,igaschem)
#endif

  call allocate_mod_cloud_s1
  call allocate_mod_rad_common(ichem)
  call allocate_mod_rad_aerosol(ichem)
  call allocate_mod_rad_o3blk
  call allocate_mod_rad_outrad
  if ( irrtm == 1 ) then
    call allocate_mod_rad_rrtmg
  else
    call allocate_mod_rad_radiation 
    call allocate_mod_rad_colmod3(ichem)
  end if

  call allocate_mod_che_common(ichem)
  call allocate_mod_che_mppio
  call allocate_mod_che_dust(ichem)
  call allocate_mod_che_bdyco
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------ALLOCATE NEEDED SPACE---------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
! 
  if ( myid == italk ) then
    write (aline,*) 'starting first checks' 
    call say
    if ( mod(idnint(dtrad*60.0D0),idnint(dt)) /= 0 ) then
      write (aline,*) 'DTRAD=' , dtrad , 'DT=' , dt
      call say
      call fatal(__FILE__,__LINE__, &
                 'INCONSISTENT RADIATION TIMESTEPS SPECIFIED')
    end if
    if ( mod(idnint(dtsrf),idnint(dt)) /= 0 ) then
      write (aline,*) 'DTSRF=' , dtsrf , 'DT=' , dt
      call say
      call fatal(__FILE__,__LINE__, &
                 'INCONSISTENT SURFACE TIMESTEPS SPECIFIED')
    end if
    if ( ifsrf ) then
      if ( mod(idnint(srffrq*secph),idnint(dtsrf)) /= 0 ) then
        write (aline,*) 'BATFRQ=' , srffrq , 'DTSRF=' , dtsrf
        call say
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
    if ( lakemod == 1 ) then
      if ( lakfrq < srffrq .or. &
           mod(idnint(lakfrq),idnint(srffrq)) /= 0 ) then
        write (aline,*) 'BATFRQ=' , srffrq , ' LAKFRQ=' , lakfrq
        call say
        write (aline,*) 'Lake frequency needs to be an integer ',&
                        ' multiple of srffrq.'
        call say
        if ( myid == italk ) then
          call fatal(__FILE__,__LINE__, &
                    'INCONSISTENT LAKE/SURFACE TIMESTEPS SPECIFIED')
        end if
      end if
    end if
    if ( mod(idnint(dtabem*secph),idnint(dt)) /= 0 ) then
      write (aline,*) 'DTABEM=' , dtabem , 'DT=' , dt
      call say
      call fatal(__FILE__,__LINE__, &
                 'INCONSISTENT ABS/EMS TIMESTEPS SPECIFIED')
    end if
    if ( mod(idnint(dtabem*60.0D0),idnint(dtrad)) /= 0 ) then
      write (aline,*) 'DTABEM=' , dtabem , 'DTRAD=' , dtrad
      call fatal(__FILE__,__LINE__,                                   &
                 'INCONSISTENT LONGWAVE/SHORTWAVE RADIATION'//        &
                 ' TIMESTEPS SPECIFIED')
    end if
    if ( ichem == 1 .and. chemfrq <=  d_zero) then
      write (aline,*) 'CHEMFRQ=' ,chemfrq
      call say
      call fatal(__FILE__,__LINE__,'CHEMFRQ CANNOT BE ZERO')
    end if
  end if
!
  if ( ifrest ) then
    doing_restart = .true.
  end if

  if ( myid == iocpu ) then
    do i = 1 , n_atmvar
      atm_variables(i)%enabled = atm_enablevar(i)
    end do
    do i = 1 , n_srfvar
      srf_variables(i)%enabled = srf_enablevar(i)
    end do
    do i = 1 , n_stsvar
      sts_variables(i)%enabled = sts_enablevar(i)
    end do
    do i = 1 , n_lakvar
      lak_variables(i)%enabled = lak_enablevar(i)
    end do
    do i = 1 , n_subvar
      sub_variables(i)%enabled = sub_enablevar(i)
    end do
    do i = 1 , n_radvar
      rad_variables(i)%enabled = rad_enablevar(i)
    end do
    if ( lakemod /= 1 .and. iseaice /= 1 ) then
      srf_variables(ivarname_lookup('SRF','seaice'))%enabled = .false.
    end if
  end if
!
!.....calculate the time step in minutes.
!
  dtsec = dt
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
  nchefrq = idnint(secph*chemfrq)
  nbdyfrq = idnint(dtbdys)
  klak = idnint(lakfrq/srffrq)
  chfrovrradfr = chemfrq/radfrq

  ntsrf = idnint(dtsrf/dtsec)
  ntrad = idnint(dtrad/(dtsec/secpm))

  ktau = 0

  khour = 3600_8/idnint(dtsec)
  krep  = khour*3
  kbdy  = nbdyfrq/idnint(dtsec)
  katm  = natmfrq/idnint(dtsec)
  ksrf  = nsrffrq/idnint(dtsec)
  ksts  = khour*24
  krad  = nradfrq/idnint(dtsec)
  kche  = nchefrq/idnint(dtsec)
  kdbg  = ndbgfrq/idnint(dtsec)
  ksav  = nsavfrq/idnint(dtsec)

  fdaysrf = real(secpd/dtsrf)

  mtau = idnint((hspan*secph)/dt)

  nbdytime = 0
  xbctime = d_zero

  do ns = 1 , nsplit
    dtsplit(ns) = dt*(d_half/dble(nsplit-ns+1))
    dtau(ns) = dtsplit(ns)
  end do
  write (aline, *) 'dtau = ' , dtau
  call say
  ntabem = idnint(secph*dtabem/dt) !dtabem is time interval abs./emis. calc.
  dt2 = d_two*dt
!
  intmdl = rcm_time_interval(idnint(dt),usec)
  intbdy = rcm_time_interval(ibdyfrq,uhrs)
  deltmx = dt
!.....compute the time steps for radiation computation.
!sb   lake model mods
!.....compute the time steps for lake model call.
  dtlake = dtsrf
!sb   end lake model mods
!

  if ( myid == iocpu ) then
    call read_domain_info
    sigma(:) = mddom_io%sigma(:)
  end if

  call mpi_bcast(ds,1,mpi_real8,iocpu,mycomm,ierr)
  call mpi_bcast(ptop,1,mpi_real8,iocpu,mycomm,ierr)
  call mpi_bcast(sigma,kzp1,mpi_real8,iocpu,mycomm,ierr)

  dx = ds * d_1000
 
!rst-fix
  appdat = tochar(idate0)
  write (aline, *) 'initial date of the global '// &
                   'simulation: idate  = ' , appdat
  call say
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
  conf = d_one

  ! Calculate boundary areas per processor

  call setup_boundaries(cross,ba_cr)
  call setup_boundaries(dot,ba_dt)

  call set_scenario(scenario)

  call init_advection(mddom,sfs,atm1,qdot,kpbl)
  call init_precip(atms,atm2,aten,sfs,pptnc,cldfra,cldlwc)
  call init_cloud_s1(atms)
#ifdef CLM
  allocate(landmask(jx,iy))
  call init_clm(dtsec,ksrf,ichem,iemiss,mddom,mddom_io, &
                atms,sfs,ts0_io,zpbl,landmask)
#else
  call init_bats(dtsec,ksrf,ichem,iemiss,mddom,atms,sfs,zpbl)
#endif
  call init_cuscheme(ichem,dtsec,ntsrf,mddom,atm1,aten,atms,chiten,  &
                     sfs,qdot,pptc,ldmsk,sigma,hsigma,dsigma,qcon,cldfra,cldlwc)
  if ( ichem == 1 ) then
#ifdef CLM
    call init_chem(ifrest,idirect,dtsec,dx,chemfrq,dtrad,dsigma,atms,    &
                   mddom,sfs,ba_cr,fcc,cldfra,rembc,remrat,hsigma,anudg, &
                   twt,ptop,coszrs,iveg,svegfrac2d,sfracv2d,sfracb2d,    &
                   sfracs2d,solis,sdeltk2d,sdelqk2d,ssw2da,convpr,       &
                   icumtop,icumbot,taucldsp,voc_em,dep_vels)
#else
    call init_chem(ifrest,idirect,dtsec,dx,chemfrq,dtrad,dsigma,atms,    &
                   mddom,sfs,ba_cr,fcc,cldfra,rembc,remrat,hsigma,anudg, &
                   twt,ptop,coszrs,iveg,svegfrac2d,sfracv2d,sfracb2d,    &
                   sfracs2d,solis,sdeltk2d,sdelqk2d,ssw2da,convpr,       &
                   icumtop,icumbot,taucldsp)
#endif
  end if
  call init_rad(ichem,ptop,hsigma,sigma,twt,atms,sfs,mddom,sabveg,solis, &
                coszrs,aldirs,aldifs,aldirl,aldifl,albvs,albvl,aemiss,   &
                sinc,solvs,solvd,fsw,flw,flwd,ldmsk1,chia,chtrname)
#ifdef CLM
  call init_rad_clm(sols2d,soll2d,solsd2d,solld2d)
#endif
!
  if ( myid == italk ) then
    if ( ifrest .and. idate0 == idate1 ) then
      write (6,*) 'Error in parameter set.'
      write (6,*) 'Cannot set idate0 == idate1 on restart run'
      write (6,*) 'Correct idate0.'
      call fatal(__FILE__,__LINE__,'IDATE0==IDATE1 ON RESTART')
    end if
  end if
!
  appdat = tochar(idate1)
  write (aline,*) 'initial date of this '//'simulation: ' , appdat
  call say
  appdat = tochar(idate2)
  write (aline,*) '  final date of this '//'simulation: ' , appdat
  call say
  write (aline,'(a,i10,a)')  &
       'total simulation lenght ' , hspan, ' hours'
  call say
  write (aline,'(a,f9.4)')  &
       'dtsec (timestep in seconds)' , dtsec
  call say
  idatex = idate1
  call split_idate(idatex,xyear,xmonth,xday,xhour)
  kstsoff = khour*xhour
!
  write (aline, *) 'input/output parameters '
  call say
  write (aline,*) 'if true(T) create SAV files for '// &
                  'restart: ifsave = ' , ifsave  
  call say 
  write (aline,*) 'Frequency in hours to create SAV: savfrq = ' , savfrq
  call say 
  write (aline,*) 'if true (T) Output ATM files:  ifatm = ' , ifatm 
  call say 
  write (aline,*) 'Frequency in hours to write  ATM: atmfrq = ' , atmfrq 
  call say  
  write (aline,*) 'Frequency in hours to write  RAD: radfrq = ' , radfrq 
  call say
  write (aline,*) 'Frequency in hours to write  SRF: srffrq = ' , srffrq  
  call say
  if ( lakemod == 1 ) then
    write (aline,*) 'Frequency in hours to write  LAK: lakfrq = ' , lakfrq
  end if
  write (aline,*) 'if true (T) output CHEM files:  ifchem = ' , ifchem 
  call say 
  write (aline,*) 'Frequency in hours to write CHEM: chemfrq =' , chemfrq
  call say  
  write (aline,*) 'Frequency in hours to write CLM: clmfrq = ', clmfrq
  call say
  write (aline,*) ' '
  call say
  write (aline,*) 'physical parameterizations '
  call say
  write (aline,'(a,i2)') ' Lateral Boundary conditions '// &
                          'scheme: iboudy = ' , iboudy
  call say  
  write (aline,'(a,i2)') ' Cumulus convection scheme: icup = ' , icup
  call say
  write  (aline,'(a,i2)') ' Grell Scheme Cumulus closure '//  &
                          'scheme: igcc =' , igcc 
  call say
  write  (aline,'(a,i2)') ' Moisture scheme: ipptls = ' , ipptls 
  call say
  write  (aline,'(a,i2)') ' Ocean Flux scheme: iocnflx = ' , iocnflx
  call say
  if ( iocnflx == 2 ) then
    write  (aline,'(a,i2)') ' Zeng roughness formula: iocnrough = ', iocnrough
    call say
  end if
  write  (aline,'(a,i2)') ' Coupling with ROMS ocean model: iocncpl = ', iocncpl
  call say
  write  (aline,'(a,i2)') ' Pressure gradient force scheme: ipgf = ' , ipgf 
  call say
  write  (aline,'(a,i2)') ' Prescribed a surface LW emissivity: '// &
                          'iemiss = ' , iemiss 
  call say 
  write  (aline,'(a,i2)') ' Use lake model lakemod = ' , lakemod 
  call say
  write  (aline,'(a,i2)') ' Use Chemical/aerosol model '// &
                          '(0=no,1=yes):  ichem =' , ichem 
  call say
  write  (aline,'(a,i2)') ' Use diurnal sst cycle effect '// &
                          '(0=no,1=yes):  idcsst =' , idcsst 
  call say 
  write  (aline,'(a,i2)') ' Use sea ice effect '// &
                          '(0=no,1=yes):  iseaice =' , iseaice 
  call say
  write  (aline,'(a,i2)') ' Use desert seasonal effect '// &
                          '(0=no,1=yes):  idesseas =' , idesseas 

  call say
  write  (aline,'(a,i2)') ' Use convective LWP parameterization '// &
                     'for large-scale '// &
                     'clouds (0=no,1=yes):  iconvlwp =' , iconvlwp 
  call say
  write  (aline,'(a,f9.6)') ' Nudge value high range   =', high_nudge 
  call say
  write  (aline,'(a,f9.6)') ' Nudge value medium range =', medium_nudge 
  call say
  write  (aline,'(a,f9.6)') ' Nudge value low range    =', low_nudge 
  call say
#ifdef CLM 
   write  (aline,'(a,i2)' ) '  imask=' , imask 
  call say
#endif
  write (aline, *) ' '
  call say
  write (aline, *) 'model parameters '
  call say
  write (aline,'(a,f12.6)')  ' time step for radiation '// &
          'model in minutes:  dtrad = ' , dtrad 
  call say 
  write (aline,'(a,f12.6)')  ' time step for land surface '// &
          'model in seconds :  dtsrf  = ' , dtsrf 
  call say
  write (aline,'(a,f12.6)')  ' time step for LW absorption/'// &
          'emissivity in hours:  dtabem  = ' , dtabem 
  call say
  write (aline,'(a,f12.6)')  ' time step for atmosphere model '// &
          ' in seconds:  dt     = ' , dt
  call say
  write (aline, *) ' '
  call say
  write (aline,'(a,i2)' ) '# of bottom model levels '// &
                          'with no clouds:  ncld = ' , ncld
  call say
  write (aline, *) ' '
  call say

  if (ichem==1) then
    do n = 1 , ntr
      call mpi_bcast(chtrname(n),6,mpi_character,iocpu,mycomm,ierr)
    end do
  end if

  write (aline, *) 'Reading in DOMAIN data'
  call say

  if ( myid == iocpu ) then
    if ( nsg > 1 ) then
      call read_subdomain(ht1_io,lndcat1_io,xlat1_io,xlon1_io)
      ht1_io(:,:,:) = ht1_io(:,:,:)*egrav
#ifndef CLM
      if ( lakemod == 1 ) call read_subdomain_lake(dhlake1_io)
#endif
    else
#ifndef CLM
      if ( lakemod == 1 ) call read_domain_lake(dhlake1_io)
#endif
      do i = idot1 , idot2
        do j = jdot1 , jdot2
          ht1_io(1,j,i) = mddom_io%ht(j,i)*egrav
          lndcat1_io(1,j,i) = mddom_io%lndcat(j,i)
          xlat1_io(1,j,i) = mddom_io%xlat(j,i)
          xlon1_io(1,j,i) = mddom_io%xlon(j,i)
        end do
      end do
    end if
    call close_domain
!
!------invert mapscale factors and convert hgt to geopotential
!
    do i = idot1 , idot2
      do j = jdot1 , jdot2
        mddom_io%ht(j,i)   = mddom_io%ht(j,i)*egrav
        mddom_io%msfd(j,i) = d_one/mddom_io%msfd(j,i)
        mddom_io%msfx(j,i) = d_one/mddom_io%msfx(j,i)
      end do
    end do

    print * , ' '
    print * ,                                                     &
          '***************************************************'
    print * ,                                                     &
          '***************************************************'
    print * ,                                                     &
          '**** RegCM IS BEING RUN ON THE FOLLOWING GRID: ****'
    print * , '****     Map Projection: ' , iproj ,               &
          '                ****'
    print * , '****     IY=' , iy , ' JX=' , jx , ' KX=' , kz ,   &
          '             ****'
    print * , '****     PTOP=' , ptop , ' DX=' , ds ,             &
          '       ****'
    print * , '****     CLAT= ' , clat , ' CLON=' , clon ,        &
          '    ****'
    if ( iproj == 'ROTMER' ) print * , '****     PLAT= ' , plat , &
                                  ' PLON=' , plon , '    ****'
    print * ,                                                     &
          '***************************************************'
    print * , ' '

  end if  ! end if (myid == iocpu)

  call grid_distribute(mddom_io%ht,mddom%ht,jde1,jde2,ide1,ide2)
  call grid_distribute(mddom_io%lndcat,mddom%lndcat,jde1,jde2,ide1,ide2)
  call grid_distribute(mddom_io%xlat,mddom%xlat,jde1,jde2,ide1,ide2)
  call grid_distribute(mddom_io%xlon,mddom%xlon,jde1,jde2,ide1,ide2)
  call grid_distribute(mddom_io%dlat,mddom%dlat,jde1,jde2,ide1,ide2)
  call grid_distribute(mddom_io%dlon,mddom%dlon,jde1,jde2,ide1,ide2)
  call grid_distribute(mddom_io%msfx,mddom%msfx,jde1,jde2,ide1,ide2)
  call grid_distribute(mddom_io%msfd,mddom%msfd,jde1,jde2,ide1,ide2)
  call grid_distribute(mddom_io%coriol,mddom%coriol,jde1,jde2,ide1,ide2)

  call subgrid_distribute(ht1_io,ht1,jde1,jde2,ide1,ide2)
  call subgrid_distribute(lndcat1_io,lndcat1,jde1,jde2,ide1,ide2)
  call subgrid_distribute(xlat1_io,xlat1,jde1,jde2,ide1,ide2)
  call subgrid_distribute(xlon1_io,xlon1,jde1,jde2,ide1,ide2)

#ifndef CLM
  if ( lakemod == 1 ) then
    call subgrid_distribute(dhlake1_io,dhlake1,jci1,jci2,ici1,ici2)
  endif
#endif
!
  call exchange(mddom%ht,1,jde1,jde2,ide1,ide2)
  call exchange(mddom%msfx,2,jde1,jde2,ide1,ide2)
  call exchange(mddom%msfd,2,jde1,jde2,ide1,ide2)
!
!-----compute land/water mask on subgrid space
!
   do i = ici1 , ici2
     do j = jci1 , jci2
       if ( mddom%lndcat(j,i) > 13.5D0 .and. &
            mddom%lndcat(j,i) < 15.5D0 ) then
         ldmsk(j,i) = 0
         do n = 1, nnsg
           ldmsk1(n,j,i) = 0
         end do
       else
         ldmsk(j,i) = 1
         do n = 1, nnsg
           ldmsk1(n,j,i) = 1
         end do
       end if
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
  write (aline, '(a,i3)') 'Index of highest allowed pbl:'// &
                          '  kmxpbl = ' , kmxpbl
  call say
  write (aline, *) ' '
  call say
!
  if ( ipptls == 1 ) then
    write (aline, '(a,f11.6,a,f11.6)')                        &
             'Auto-conversion rate:  Land = ' , qck1land ,    &
             '                      Ocean = ' , qck1oce
    call say
    write (aline, *) 'Relative humidity thresholds:  Land = ' ,   &
             rh0land , '                              Ocean = ' , &
             rh0oce
    call say
    write (aline, *) 'Gultepe factors:  Land=' , gulland ,   &
            '                 Ocean=' , guloce
    call say
    write (aline, *) 'Maximum cloud cover for radiation: ' , fcmax
    call say
    write (aline, *) 'Maximum relative humidity: ' , rhmax
    call say
    write (aline, *) 'rh0 temperature threshold: ' , tc0
    call say
    if ( cevap <= d_zero ) then
      write (aline, *) 'Raindrop evaporation not included'
      call say
      cevap = d_zero
    end if
    cevapu = cevap
    if ( caccr <= d_zero ) then
      write (aline, *) 'Raindrop accretion not included'
      call say
    end if
    write (aline, *) 'Raindrop Evaporation Rate' , cevap
    call say
    write (aline, *) 'Raindrop Accretion Rate' , caccr
    call say
    write (aline, *) 'Maximum total cloud cover for radiation: ' ,  &
         cftotmax
    call say
  end if
 
  write (aline, *) ' '
  call say
 
  if (icup == 99) then
    write (aline,*) 'Variable cumulus scheme: will use Grell '// &
         'over land and Emanuel over ocean.'
    call say
    do i = ici1 , ici2
      do j = jci1 , jci2
        if ( mddom%lndcat(j,i) > 14.5D0 .and. &
             mddom%lndcat(j,i) < 15.5D0 ) then
          cucontrol(j,i) = 4
        else
          cucontrol(j,i) = 2
        end if
      end do
    end do
  end if
  if (icup == 98) then
    write (aline,*) 'Variable cumulus scheme: will use Emanuel '// &
         'over land and Grell over ocean.'
    call say
    do i = ici1 , ici2
      do j = jci1 , jci2
        if ( mddom%lndcat(j,i) > 14.5D0 .and. &
             mddom%lndcat(j,i) < 15.5D0 ) then
          cucontrol(j,i) = 2
        else
          cucontrol(j,i) = 4
        end if
      end do
    end do
  end if

  if ( icup == 1 ) then
    call allocate_mod_cu_kuo
    write (aline, *) '*********************************'
    call say
    write (aline, *) '***** Anthes-Kuo Convection *****'
    call say
    write (aline, *) '*********************************'
    call say
!
!   specify heating profile (twght) and weighted function
!   for moisture fluxes due to convection (vqflx)
!   assume base of cloud varies as  < kbase = 5,kz >
!   top  of cloud varies as  < ktop  = 1,kbase-3 >
!   exceptions to this are treated explicitly in subroutine
!   "cupara".
!
    do kbase = 5 , kz
      do ktop = 1 , kbase - 3
        do k = 1 , kz
          twght(k,kbase,ktop) = d_zero
          vqflx(k,kbase,ktop) = d_zero
        end do
!
!       get twght from 1/2 level sigma values
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
!       get vqflx from  d(w*q) / dsigma on full levels
!       do computations in p to avoid sigma=0. discontinuity
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
  end if
  if ( icup == 2 .or. icup == 99 .or. icup == 98 ) then
    call allocate_mod_cu_grell
    kbmax = kz
    do k = 1 , kz - 1
      if ( hsigma(k) <= skbmax ) kbmax = kz - k
    end do
    write (aline, *) '*********************************' 
    call say
    write (aline, *) '*****    Grell Convection   *****'
    call say
    write (aline, *) '*********************************'
    call say
    write (aline, *) '   Max Shear:' , shrmax
    call say
    write (aline, *) '   Min Shear:' , shrmin
    call say
    write (aline, *) '   Max PPT eff:' , edtmax
    call say
    write (aline, *) '   Min PPT eff:' , edtmin
    call say
    write (aline, *) '   Max PPT eff(o):' , edtmaxo
    call say
    write (aline, *) '   Min PPT eff(o):' , edtmino
    call say
    write (aline, *) '   Max PPT eff(x):' , edtmaxx
    call say
    write (aline, *) '   Min PPT eff(x):' , edtminx
    call say
    write (aline, *) '   Max PBC (pbcmax):' , pbcmax
    call say
    write (aline, *) '   Min Cloud Depth:' , mincld
    call say
    write (aline, *) '   Max Cloud Base:' , kbmax , skbmax
    call say
    write (aline, *) '   Max Heating:' , htmax
    call say
    write (aline, *) '   Min Heating:' , htmin
    call say
    if ( igcc == 1 ) then
      write (aline, *)                                              &
            '   Arakawa-Schubert (1974) Closure Assumption'
      call say
    else if ( igcc == 2 ) then
      write (aline, *)                                              &
            '   Fritsch-Chappell (1980) Closure Assumption'
      call say
      write (aline, *) '     ABE removal timescale: dtauc=' , dtauc
      call say
    end if

    write (aline, *) '*********************************'
    call say
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
  if ( icup == 3 ) then
    call allocate_mod_cu_bm
    write (aline,*) ' The Betts-Miller Convection scheme is not' ,  &
                    ' properly implemented'
    call say
    !call fatal(__FILE__,__LINE__,'BETTS-MILLER NOT WORKING')
  end if
  if ( icup == 4 .or. icup == 99 .or. icup == 98 ) then
    call allocate_mod_cu_em
    minorig = kz
    do k = 1 , kz
      if ( hsigma(k) <= minsig ) minorig = kz - k
    end do
    write (aline, *) ' '
    call say
    write (aline, *) 'EMANUEL (1991) CONVECTION V4.3C (20 May, 2002)'
    call say
    write (aline, *) '  MIN CONVECTION ORIGIN (minsig/orig): ', minsig, minorig
    call say
    write (aline, *)'  AUTOCONVERSION THERSHOLD (elcrit): ' , elcrit
    call say
    write (aline, *)'  AUTOCONVERSION THRESHOLD TO ZERO (tlcrit): ', tlcrit
    call say
    write (aline, *)'  ENTRAINMENT COEFFICIENT (entp): ' , entp
    call say
    write (aline, *) '  FRACTIONAL AREA OF UNSATURATED DNDRAFT (sigd): ', sigd
    call say
    write (aline, *)'  PRECIP FRACTION OUTSIDE OF CLOUD (sigs): ', sigs
    call say
    write (aline, *)'  FALL SPEED OF RAIN (omtrain): ' , omtrain
    call say
    write (aline, *)'  FALL SPEED OF SNOW (omtsnow): ' , omtsnow
    call say
    write (aline, *)'  RAIN EVAPORATION COEFFICIENT (coeffr): ' , coeffr
    call say
    write (aline, *)'  SNOW EVAPORATION COEFFICIENT (coeffs): ' , coeffs
    call say
    write (aline, *) '  CONVECTIVE MOMENTUM TRANSPORT COEFFICIENT (cu): ', cu
    call say
    write (aline, *)'  DOWNDRAFT VELOCITY SCALE (betae): ' , betae
    call say
    write (aline, *) '  MAX NEGATIVE PERTURBATION BELOW LFC (dtmax): ' , dtmax
    call say
    write (aline, *)'  QUASI-EQUILIBRIUM APPROACH RATE (alphae): ' , alphae
    call say
    write (aline, *)'  QUASI-EQUILIBRIUM APPROACH RATE (damp): ' , damp
    call say
    write (aline, *) ' '
    call say
    write (aline, *) ' '
    call say
  end if
  if ( icup == 5 ) then
    call allocate_mod_cu_tiedtke
    write (aline, *) ' '
    call say
    write (aline, *) 'TIEDTKE (1986) CONVECTION SCHEME FROM ECHAM 5.4'
    call say
    write (aline, *) '  USED SCHEME: ', iconv
    call say
    write (aline, *) '  ENTRAINMENT RATE FOR PENETRATIVE CONVECTION: ', entrpen
    call say
    write (aline, *) '  ENTRAINMENT RATE FOR SHALLOW CONVECTION: ', entrscv
    call say
    write (aline, *) '  ENTRAINMENT RATE FOR MIDLEVEL CONVECTION: ', entrmid
    call say
    write (aline, *) '  ENTRAINMENT RATE FOR CUMULUS DOWNDRAFT: ', entrdd
    call say
    write (aline, *) '  RELAT. CLOUD MASSFLUX AT LEVEL ABV NONBUOY: ', cmfctop
    call say
    write (aline, *) '  MAXIMUM MASSFLUX VALUE ALLOWED: ', cmfcmax
    call say
    write (aline, *) '  MINIMUM MASSFLUX VALUE ALLOWED: ', cmfcmin
    call say
    write (aline, *) '  FRACTIONAL MASSFLUX FOR DOWNDRAFTS AT LFS: ', cmfdeps
    call say
    write (aline, *) '  RELATIVE SATURATION IN DOWNDRAFTS: ', rhcdd
    call say
    write (aline, *) '  CAPE ADJUSTMENT TIMESCALE PARAMETER: ', cmtcape
    call say
    write (aline, *) '  RESTRICT RAINFALL LEVEL: ', zdlev
    call say
    write (aline, *) '  COEFF. FOR CONV. FROM CLOUD WATER TO RAIN: ', cprcon
    call say
    write (aline, *) '  MAX. LEVEL FOR CLOUD BASE OF MID LEVEL CONV.: ', nmctop
    call say
    write (aline, *) '  PENETRATIVE CONVECTION IS SWITCHED ON: ', lmfpen
    call say
    write (aline, *) '  SHALLOW CONVECTION IS SWITCHED ON: ', lmfscv
    call say
    write (aline, *) '  MIDLEVEL CONVECTION IS SWITCHED ON: ', lmfmid
    call say
    write (aline, *) '  CUMULUS DOWNDRAFT IS SWITCHED ON: ', lmfdd
    call say
    write (aline, *) '  CUMULUS FRICTON IS SWITCHED ON: ', lmfdudv
    call say
  end if

  if ( ibltyp == 1 .or. ibltyp == 99 ) then
    call allocate_mod_pbl_holtbl(ichem,ichdrdepo)
  end if
  if ( ibltyp == 2 .or. ibltyp == 99 ) then
    call init_mod_pbl_uwtcm
  end if

  call init_pbl(atm2,atms,aten,holtten,uwten,adf,heatrt,chiten,remdrd,   &
                cchifxuw,psdot,sfs,mddom,ldmsk,hsigma,sigma,dsigma,ptop, &
                chtrdpv,chtrname,ichem,ichdrdepo,dt)
 
!     Convective Cloud Cover
  afracl = 0.3D0 ! frac. cover for conv. precip. when dx=dxlarg
  afracs = 1.0D0 !   "     "    "    "      "     "   dx=dxsmal
  dlargc = 200.0D0
  dsmalc = 10.0D0
  dxtemc = dmin1(dmax1(dx,dsmalc),dlargc)
  clfrcv = afracl + (afracs-afracl)                                 &
           *((dlargc-dxtemc)/(dlargc-dsmalc))**d_two
  clfrcv = dmin1(clfrcv,clfrcvmax)
  write (aline, *) ' '
  call say
  write (aline, *) 'CONVECTIVE CLOUD FRACTION/WATER'
  call say
  write (aline, *) '   Maximum Convective Cloud Cover'
  call say
  write (aline, *) '     before resolution scaling: ' , clfrcvmax
  call say
  write (aline, *) '   Maximum Convective Cloud Cover'
  call say
  write (aline, *) '     after resolution scaling: ' , clfrcv
  call say
  write (aline, *) '   Convective Cloud Water: ' , cllwcv
  call say
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
 
!
!-----compute the k level under which the maximum equivalent potential
!     temperature will be regarded as the origin of air parcel that
!     produces cloud (used in the cumulus parameterization scheme).
!
  sig700 = (70.0D0-ptop)/(d_100-ptop)
  do k = 1 , kz
    k700 = k
    if ( sig700 <= sigma(k+1) .and. sig700 > sigma(k) ) exit
  end do
!
!-----specify the coefficients for sponge boundary conditions.
!
!
! Move back form init : the land categories are read just here,
! they are no more in the restart file.
!
  if ( ipptls == 1 ) then
    do i = ici1 , ici2
      do j = jci1 , jci2
        if ( mddom%lndcat(j,i) > 14.5D0 .and. &
             mddom%lndcat(j,i) < 15.5D0) then
          qck1(j,i) = qck1oce  ! OCEAN
          cgul(j,i) = guloce   ! OCEAN
          rh0(j,i)  = rh0oce    ! OCEAN
        else
          qck1(j,i) = qck1land ! LAND
          cgul(j,i) = gulland  ! LAND
          rh0(j,i)  = rh0land   ! LAND
        end if
      end do
    end do
  end if
!
! Setup Boundary condition routines.
!
  call setup_bdycon(hsigma)
!
  if ( icup == 3 ) call lutbl(ptop)
!
!-----print out the parameters specified in the model.
!
  if ( myid == italk ) then
    if ( ibltyp < 0 .or. ibltyp > 2 ) then
      call fatal(__FILE__,__LINE__,'UNSUPPORTED PBL SCHEME.')
    end if
    if ( ibltyp == 0 ) print 99002
    print 99003 , ntrad
!
    if ( iboudy < 0 .or. iboudy > 5 ) then
      call fatal(__FILE__,__LINE__,'UNSUPPORTED BDY SCHEME.')
    end if
    if ( iboudy == 0 ) print 99004
    if ( iboudy == 1 ) print 99005 , fnudge , gnudge
    if ( iboudy == 2 ) print 99007
    if ( iboudy == 3 ) print 99008
    if ( iboudy == 4 ) print 99009
    if ( iboudy == 5 ) print 99006 , fnudge , gnudge
 
    print 99010
!
    do k = 1 , kz
      print 99011 , k , sigma(k) , hsigma(k) , dsigma(k) , twt(k,1) ,    &
            twt(k,2) , qcon(k)
    end do
    print 99012 , kzp1 , sigma(kzp1)
    print 99014 , dt
    print 99015 , dx
    print 99016 , jx , iy
    print 99017 , kz
    print 99018 , xkhz
    print 99019 , xkhmax
  end if

  call allocate_v2dbound(xpsb,cross)
  call allocate_v3dbound(xtb,kz,cross)
  call allocate_v3dbound(xqb,kz,cross)
  call allocate_v3dbound(xub,kz,dot)
  call allocate_v3dbound(xvb,kz,dot)

  call mpi_barrier(mycomm,ierr) 
!
  call time_end(subroutine_name,idindx)
!
99002 format (/'   frictionless and insulated for the lower boundary.')
99003 format (                                                            &
      '     the surface energy budget is used to calculate the ground '// &
      'temperature. the radiation is computed every ',i4,' time steps.')
!1100 format('     heat and moisture fluxes from the ground are turned o
!     1ff.')
99004 format (/'   the lateral boundary conditions are fixed.')
99005 format (/'   relaxation boundary conditions (linear method)',     &
          ' are used.  fnudge = ',e15.5,'  gnudge = ',e15.5)
99006 format (/'   relaxation boudnary conditions (exponential method)',&
          ' are used. fnudge = ',e15.5,'  gnudge = ',e15.5)
99007 format (/'   time dependent boundary conditions are used.')
99008 format (/'   inflow/outflow boundary conditions are used.')
99009 format (/'   sponge boundary conditions are used.')
99010 format ('0 k',4x,'sigma(k)',3x,'  a(k)',5x,'dsigma(k)',4x,        &
          'twt(k,1)',5x,'twt(k,2)',5x,'qcon(k)'/)
99011 format (1x,i2,5x,f7.4,5x,f7.4,5x,f7.4,5x,f8.4,5x,f8.4,5x,f8.4)
99012 format (1x,i2,5x,f7.4)
99014 format (' time step = ',f7.2,' seconds')
99015 format (' dx = ',f7.0,' meters')
99016 format (' grid points (x,y) = (',i4,',',i4,')')
99017 format (' Number of levels = ',i2)
99018 format (' constant hor. diff. coef. = ',e12.5,' m*m/s')
99019 format (' maximum  hor. diff. coef. = ',e12.5,' m*m/s')
!
  end subroutine param
!
end module mod_params

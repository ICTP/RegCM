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

      module mod_param

      use mod_constants
      use mod_dynparam
      use mod_runparams
      use mod_pmoist
      use mod_bats
      use mod_main
      use mod_trachem
      use mod_date
      use mod_message
      use mod_cu_bm
      use mod_cu_em
      use mod_rad
      use mod_split
      use mod_slice
      use mod_pbldim
      use mod_outrad
      use mod_holtbl
      use mod_che_semdde
      use mod_aerosol
      use mod_radiation
      use mod_dust
      use mod_bdycod
      use mod_mainchem
      use mod_cvaria
      use mod_leaftemp
      use mod_o3blk
      use mod_ncio
      use mod_scenarios
      use mod_diagnosis

#ifdef MPP1
      use mod_mppio
#ifdef CLM
      use mod_clm
#endif
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
#ifdef MPP1
#ifndef IBM
      use mpi
#else 
      include 'mpif.h'
#endif 
#endif
      implicit none
!
! Local variables
!
      real(8) :: afracl , afracs , bb , cc , chibot , delsig , &
               & dlargc , dsmalc , dxtemc , pk , ptmb , pz , qk ,       &
               & qkp1 , sig700 , sigtbl , ssum , vqmax , vqrang , wk ,  &
               & wkp1 , xbot , xtop , xx , yy
      real(8) , dimension(nsplit) :: dtsplit
      integer :: i , j , k , kbase , ktop , ns
      character(5) , dimension(maxntr) :: inpchtrname
      real(8) , dimension(maxntr) :: inpchtrsol
      real(8) , dimension(maxntr,2) :: inpchtrdpv
      real(8) , dimension(maxnbin,2) :: inpdustbsiz
      integer :: len_path
      logical :: lband , lmpi
#ifdef MPP1
      integer :: n , ierr
#ifndef CLM
      integer :: imask
      real(8) :: clmfrq
#endif
#else
      integer :: imask
      real(8) :: clmfrq
#endif
!
!
!----------------------------------------------------------------------
!-----vqrang is the range limit on vqflx.
!
      data vqrang /5.0E-4/
!
!----------------------------------------------------------------------
!-----namelist:
!
      namelist /restartparam/ ifrest , idate0 , idate1 , idate2
 
      namelist /timeparam/ radfrq , abatm , abemh , dt
 
!chem2
      namelist /outparam/ ifsave , savfrq , iftape , tapfrq , ifrad ,   &
      & radisp , ifbat , ifsub , batfrq , ifchem , chemfrq , dirout
!chem2
      namelist /physicsparam/ ibltyp , iboudy , icup , igcc , ipgf ,    &
      & iemiss , lakemod , ipptls , iocnflx , ichem, high_nudge,        &
      & medium_nudge, low_nudge , scenario , idcsst , iseaice ,         &
      & idesseas
!chem2_
      namelist /subexparam/ ncld , fcmax , qck1land , qck1oce ,         &
      & gulland , guloce , rhmax , rh0oce , rh0land , cevap , caccr ,   &
      & tc0 , cllwcv , clfrcvmax
 
      namelist /grellparam/ shrmin , shrmax , edtmin , edtmax ,         &
      & edtmino , edtmaxo , edtminx , edtmaxx , pbcmax , mincld ,       &
      & htmin , htmax , skbmax , dtauc
 
      namelist /emanparam/ minsig , elcrit , tlcrit , entp , sigd ,     &
      & sigs , omtrain , omtsnow , coeffr , coeffs , cu , betae ,       &
      & dtmax , alphae , damp
 
!chem2
      namelist /chemparam/ ichremlsc , ichremcvc , ichdrdepo ,          &
      & ichcumtra , idirect , mixtype , inpchtrname , inpchtrsol ,      &
      & inpchtrdpv , inpdustbsiz
!chem2_

#ifdef CLM
      namelist /clmparam/ dirclm , imask , clmfrq
#endif
!
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
!     ktaur  : if ifrest=.false., ktaur=0.
!     if ifrest=.true., set ktaur=ktau, which is the time-step
!     the model output is saved.
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
!     = 0 ; frictionless
!     = 1 ; holtslag pbl (holtslag, 1990)
!
!     radfrq : specify the frequency in
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
!     iboudy : specify the laterial boundary conditions.
!     = 0 ; fixed.
!     = 1 ; relaxation, linear technique.
!     = 2 ; time-dependent (from observations or large-scale
!     model).
!     = 3 ; time and inflow/outflow dependent.
!     = 4 ; sponge (perkey & kreitzberg, mwr 1976).
!     = 5 ; relaxation, exponential technique.
!
!     iftape : whether you want output in a saved tape (unit 20) for
!     analyses in dataflow.
!     = .true. ; yes
!     = .false. ; no
!
!     tapfrq : if iftape=true., specify the output interval in hours.
!
!     ifrad  : whether you want radiation output saved
!     = .true. ; yes
!     = .false. ; no
!
!     radisp : if ifrad=1, specify the output interval in hours.
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
      radfrq = 30.    ! time interval in min solar rad caluclated
      abatm = 600.    ! time interval at which bats is called (secs)
      abemh = 12.     ! time interval absorption-emission calculated (hours)
      dt = 200.       ! time step in seconds
 
!-----namelist outparam:       note: * signifies currently not in namelist
!
      rfstrt = .false.      ! *
      ifsave = .false.
      savfrq = 6.
      iftape = .true.
      tapfrq = 6.
      ifrad = .true.
      radisp = 6.0       ! time interval for disposing rad output (hrs)
      ifbat = .true.
      batfrq = 1.0      ! time interval for disposing bats output (hrs)
      ifsub = .true.
      dirout = './output' 
!chem2
      ifchem = .false.
      chemfrq = 6.0   ! time interval for disposeing chem output (hrs)
!chem2_
      clmfrq = 12.0
 
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
      lakemod = 1
      ichem = 0
      imask = 1
      scenario = 'A1B'
      idcsst = 0
      iseaice = 0
      idesseas = 0
      high_nudge = 3 
      medium_nudge=2
      low_nudge=1   
!
!----------------------------------------------------------------------
!------namelist subexparam:
      ncld = 1           ! # of bottom model levels with no clouds (rad only)
      fcmax = 0.80       ! Maximum cloud fraction cover (rad only)
!     qck1land = 0.0005  ! Autoconversion Rate for Land
!     qck1oce  = 0.0005  ! Autoconversion Rate for Ocean
      qck1land = 0.00025 ! Autoconversion Rate for Land
      qck1oce = 0.00025  ! Autoconversion Rate for Ocean
      gulland = 0.4      ! Fract of Gultepe eqn (qcth) when precip occurs (land)
      guloce = 0.4       ! Fract of Gultepe eqn (qcth) for ocean
      rhmax = 1.01       ! RH at whicn FCC = 1.0
      rh0oce = 0.90      ! Relative humidity threshold for ocean
      rh0land = 0.80     ! Relative humidity threshold for land
      tc0 = 238.0        ! Below this temperature, rh0 begins to approach unity
!     cevap    = 0.2e-4  ! Raindrop evap rate coef [[(kg m-2
!     s-1)-1/2]/s]
      cevap = 1.0E-3     ! Raindrop evap rate coef [[(kg m-2 s-1)-1/2]/s]
!     caccr    = 6.0     ! Raindrop accretion rate [m3/kg/s]
      caccr = 3.0        ! Raindrop accretion rate [m3/kg/s]
      cllwcv = 0.3E-3    ! Cloud liquid water content for convective precip.
      clfrcvmax = 0.25   ! Max cloud fractional cover for convective precip.
 
!------namelist grellparam:
      shrmin = 0.25     ! Minimum Shear effect on precip eff.
      shrmax = 0.50     ! Maximum Shear effect on precip eff.
      edtmin = 0.25     ! Minimum Precipitation Efficiency
      edtmax = 1.00     ! Maximum Precipitation Efficiency
      edtmino = 0.0     ! Minimum Precipitation Efficiency (o var)
      edtmaxo = 1.00    ! Maximum Precipitation Efficiency (o var)
      edtminx = 0.25    ! Minimum Precipitation Efficiency (x var)
      edtmaxx = 1.00    ! Maximum Precipitation Efficiency (x var)
      pbcmax = 150.     ! Max depth (mb) of stable layer b/twn LCL & LFC
      mincld = 150.     ! Min cloud depth (mb).
      htmin = -250.     ! Min convective heating
      htmax = 500.      ! Max convective heating
      skbmax = 0.4      ! Max cloud base height in sigma
      dtauc = 30.       ! Fritsch & Chappell (1980) ABE Removal Timescale (min)
 
!------namelist emanparam:
      minsig = 0.95     ! Lowest sigma level from which convection can originate
      elcrit = 0.0011   ! AUTOCONVERSION THRESHOLD WATER CONTENT (gm/gm)
      tlcrit = -55.0    ! BELOW TLCRIT AUTO-CONVERSION THRESHOLD IS ZERO
      entp = 1.5        ! COEFFICIENT OF MIXING IN THE ENTRAINMENT FORMULATION
      sigd = 0.05       ! FRACTIONAL AREA COVERED BY UNSATURATED DNDRAFT
      sigs = 0.12       ! FRACTION OF PRECIPITATION FALLING OUTSIDE OF CLOUD
      omtrain = 50.0    ! FALL SPEED OF RAIN (P/s)
      omtsnow = 5.5     ! FALL SPEED OF SNOW (P/s)
      coeffr = 1.0      ! COEFFICIENT GOVERNING THE RATE OF RAIN EVAPORATION
      coeffs = 0.8      ! COEFFICIENT GOVERNING THE RATE OF SNOW EVAPORATION
      cu = 0.7          ! COEFFICIENT GOVERNING CONVECTIVE MOMENTUM TRANSPORT
      betae = 10.0      ! CONTROLS DOWNDRAFT VELOCITY SCALE
      dtmax = 0.9       ! MAX NEGATIVE PARCEL TEMPERATURE PERTURBATION BELOW LFC
      alphae = 0.2      ! CONTROLS THE APPROACH RATE TO QUASI-EQUILIBRIUM
      damp = 0.1        ! CONTROLS THE APPROACH RATE TO QUASI-EQUILIBRIUM
!
!c------namelist chemparam ; ( 0= none, 1= activated)
      ichremlsc = 1     ! tracer removal by large scale clouds
      ichremcvc = 1     ! tracer removal by convective clouds
      ichdrdepo = 1     ! tracer dry deposition
      ichcumtra = 1     ! tracer convective transport
      idirect = 1       ! tracer direct effect
      mixtype = 1
#ifdef CLM
!c------CLM Specific
      imask = 2
#endif

#ifdef MPP1
      lmpi = .true.
#else
      lmpi = .false.
#endif
#ifdef BAND
      lband = .true.
#else
      lband = .false.
#endif

!---------------------------------------------------------------------
!
#ifdef MPP1
      if ( myid.eq.0 ) then
#endif 
  
!  
!-----read in namelist variables:
      read (ipunit, restartparam)
      print * , 'param: RESTARTPARAM namelist READ IN'
      read (ipunit, timeparam)
      print * , 'param: TIMEPARAM namelist READ IN'
      read (ipunit, outparam)
      print * , 'param: OUTPARAM namelist READ IN'
      len_path = len(trim(dirout))
      if ( dirout(len_path:len_path).ne.'/' ) dirout = trim(dirout)//'/'
      read (ipunit, physicsparam)
      print * , 'param: PHYSICSPARAM namelist READ IN'
      if ( ipptls.eq.1 ) then
        read (ipunit, subexparam)
        print * , 'param: SUBEXPARAM namelist READ IN'
      end if
      if ( icup.eq.2 ) then
        read (ipunit, grellparam)
        print * , 'param: GRELLPARAM namelist READ IN'
      else if ( icup.eq.4 ) then
        read (ipunit, emanparam)
        print * , 'param: EMANPARAM namelist READ IN'
      else
      end if
      if ( ichem.eq.1 ) then
        read (ipunit, chemparam)
        print * , 'param: CHEMPARAM namelist READ IN'
      else
        ntr = 0
      end if
#ifdef CLM
      read (ipunit , clmparam)
      print * , 'param: CLMPARAM namelist READ IN'
#endif
#ifdef MPP1
      end if 

      call mpi_barrier(mpi_comm_world,ierr) 
!
!  communicate to all processors 
!
      call mpi_bcast(ifrest,1,mpi_logical,0,mpi_comm_world,ierr)
      call mpi_bcast(idate0,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(idate1,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(idate2,1,mpi_integer,0,mpi_comm_world,ierr)
 
      call mpi_bcast(radfrq,1,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_bcast(abemh,1,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_bcast(abatm,1,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_bcast(dt,1,mpi_real8,0,mpi_comm_world,ierr)
 
      call mpi_bcast(ifsave,1,mpi_logical,0,mpi_comm_world,ierr)
      call mpi_bcast(savfrq,1,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_bcast(iftape,1,mpi_logical,0,mpi_comm_world,ierr)
      call mpi_bcast(tapfrq,1,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_bcast(ifrad,1,mpi_logical,0,mpi_comm_world,ierr)
      call mpi_bcast(radisp,1,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_bcast(ifbat,1,mpi_logical,0,mpi_comm_world,ierr)
      call mpi_bcast(ifsub,1,mpi_logical,0,mpi_comm_world,ierr)
      call mpi_bcast(batfrq,1,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_bcast(ifchem,1,mpi_logical,0,mpi_comm_world,ierr)
      call mpi_bcast(chemfrq,1,mpi_real8,0,mpi_comm_world,ierr)
 
      call mpi_bcast(iboudy,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(ibltyp,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(icup,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(igcc,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(ipptls,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(iocnflx,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(ipgf,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(iemiss,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(lakemod,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(ichem,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(scenario,3,mpi_character,0,mpi_comm_world,ierr)
      call mpi_bcast(idcsst,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(iseaice,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(idesseas,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(high_nudge,1,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_bcast(medium_nudge,1,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_bcast(low_nudge,1,mpi_real8,0,mpi_comm_world,ierr)
#ifdef CLM
      call mpi_bcast(dirclm,256,mpi_character,0,mpi_comm_world,ierr)
      call mpi_bcast(imask,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(clmfrq,1,mpi_real8,0,mpi_comm_world,ierr)
#endif

      if ( ipptls.eq.1 ) then
        call mpi_bcast(ncld,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_bcast(fcmax,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(qck1land,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(qck1oce,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(gulland,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(guloce,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(rhmax,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(rh0oce,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(rh0land,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(tc0,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(cevap,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(caccr,1,mpi_real8,0,mpi_comm_world,ierr)
      end if
 
      if ( icup.eq.2 ) then
        call mpi_bcast(shrmin,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(shrmax,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(edtmin,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(edtmax,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(edtmino,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(edtmaxo,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(edtminx,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(edtmaxx,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(pbcmax,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(mincld,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(htmin,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(htmax,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(skbmax,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(dtauc,1,mpi_real8,0,mpi_comm_world,ierr)
 
      else if ( icup.eq.4 ) then
        call mpi_bcast(minsig,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(elcrit,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(tlcrit,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(entp,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(sigd,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(sigs,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(omtrain,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(omtsnow,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(coeffr,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(coeffs,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(cu,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(betae,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(dtmax,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(alphae,1,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(damp,1,mpi_real8,0,mpi_comm_world,ierr)
      else
      end if
 
      if ( ichem.eq.1 ) then
        call mpi_bcast(ichremlsc,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_bcast(ichremcvc,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_bcast(ichdrdepo,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_bcast(ichcumtra,1,mpi_integer,0,mpi_comm_world,ierr)
        call mpi_bcast(idirect,1,mpi_integer,0,mpi_comm_world,ierr)
      end if
#endif
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------ALLOCATE NEEDED SPACE---------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      call allocate_mod_che_semdde
      call allocate_mod_aerosol
      call allocate_mod_bats(lmpi,lband)
      call allocate_mod_bdycon
      call allocate_mod_holtbl
      call allocate_mod_cvaria(lmpi)
      call allocate_mod_dust
      call allocate_mod_leaftemp
      call allocate_mod_main(lmpi)
      call allocate_mod_mainchem(lmpi)
      call allocate_mod_outrad
      call allocate_mod_o3blk
      call allocate_mod_pbldim(lmpi)
      call allocate_mod_pmoist(lmpi)
      call allocate_mod_radiation 
      call allocate_mod_rad(lmpi,lband)
      call allocate_mod_slice
      call allocate_mod_split
      call allocate_mod_trachem
      call allocate_mod_runparams
#ifdef MPP1
      call allocate_mod_mppio(lband)
#ifdef CLM
      call allocate_mod_clm(lmpi,lband)
#endif
#endif
#ifndef BAND
      if (debug_level > 2) call allocate_mod_diagnosis
#endif
!
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!-----------------------ALLOCATE NEEDED SPACE---------------------------
!-----------------------------------------------------------------------
!-----------------------------------------------------------------------
!
      write (aline,*) 'param: starting first checks' 
      call say
      if ( mod(anint(radfrq*60.),anint(dt)).ne.0 ) then
        write (aline,*) 'RADFRQ=' , radfrq , 'DT=' , dt
        call say
        call fatal(__FILE__,__LINE__,                                   &
                  &'INCONSISTENT RADIATION TIMESTEPS SPECIFIED')
      end if
      if ( mod(anint(abatm),anint(dt)).ne.0 ) then
        write (aline,*) 'ABATM=' , abatm , 'DT=' , dt
        call say
        call fatal(__FILE__,__LINE__,                                   &
                  &'INCONSISTENT SURFACE TIMESTEPS SPECIFIED')
      end if
      if ( mod(anint(batfrq*3600.),anint(abatm)).ne.0 ) then
        write (aline,*) 'BATFRQ=' , batfrq , 'ABATM=' , abatm
        call say
        call fatal(__FILE__,__LINE__,                                   &
                  &'INCONSISTENT SURFACE/RADIATION TIMESTEPS SPECIFIED')
      end if
      if ( mod(anint(abemh*3600.),anint(dt)).ne.0 ) then
        write (aline,*) 'ABEMH=' , abemh , 'DT=' , dt
        call say
        call fatal(__FILE__,__LINE__,                                   &
                  &'INCONSISTENT ABS/EMS TIMESTEPS SPECIFIED')
      end if
      if ( mod(anint(abemh*60.),anint(radfrq)).ne.0 ) then
        write (aline,*) 'ABEMH=' , abemh , 'RADFRQ=' , radfrq
        call fatal(__FILE__,__LINE__,                                   &
                  &'INCONSISTENT LONGWAVE/SHORTWAVE RADIATION'//        &
                  &' TIMESTEPS SPECIFIED')
      end if

      if ( ichem.eq.0 ) ifchem = .false.
      if ( ichem.eq.1 .and. chemfrq.eq. 0.) then
        write (aline,*) 'CHEMFRQ=' ,chemfrq
        call say
        call fatal(__FILE__,__LINE__,'CHEMFRQ CANNOT BE ZERO')
      end if
!
!-----reset the options/calculate variables using namelist info:
!
      ndate0 = idate1
      nsavfrq = nint(3600.*savfrq)
      ntapfrq = nint(3600.*tapfrq)
      ndbgfrq = nint(3600.*dbgfrq)
      ktau = 0
      ktaur = 0
      xtime = 0.
      ntime = 0
      dtsplit(2) = dt/2.
      dtsplit(1) = dt/4.
      do ns = 1 , nsplit
        dtau(ns) = dtsplit(ns)
      end do
      write (aline, *) 'param: dtau = ' , dtau
      call say
      dt0 = dt      !store original dt
      nradisp = nint(radisp*3600)
                                !convert radisp to time steps
      ifrabe = nint(3600.*abemh/dt)
                                   !abemh is time interval abs./emis. calc.
      kbats = nint(3600.*batfrq)
      nbatst = nint(abatm/dt)
      dt2 = 2.*dt
!chem2
      kchem = nint(3600.*chemfrq)  ! convert chemfrq to time steps
!chem2_
!.....calculate the time step in minutes.
      dtmin = dt/60.
      deltmx = dt
!.....compute the time steps for radiation computation.
      ntrad = nint(radfrq/dtmin)
!sb   lake model mods
!.....compute the time steps for lake model call.
      dtlake = 3600.
      klake = nint(dtlake/dt)
!sb   end lake model mods
!
      call set_scenario(scenario)
!
      nstrt0 = 0
      nstart = idatediff(idate1,idate0)/ibdyfrq
      nnnend = idatediff(idate2,idate0)/ibdyfrq
      nnnchk = nstart
! 
!     ktaur = nint((NSTART-NSTRT0)*ibdyfrq*60/dtmin)
      write (aline,*) 'param: initial date of this '// &
                      'simulation: IDATE1',idate1
      call say
      write (aline,*) 'param:   final date of this '// &
                      'simulation: IDATE2' , idate2
      call say
      write (aline,'(a,i10,a)')  &
           'param: total simulation lenght ' ,  &
            idatediff(idate2,idate1) , ' hours'
      call say
      write (aline,'(a,f9.4)')  &
           'param: dtmin (timestep in minutes)' , dtmin
      call say
      ldatez = idate1
      nnnnnn = nstart
      lyear = ldatez/1000000
      lmonth = (ldatez-lyear*1000000)/10000
      lday = (ldatez-lyear*1000000-lmonth*10000)/100
      lhour = mod(ldatez,100)
      idatex = ldatez
      jyear = lyear
      jyearr = jyear
!
!-----specify the julian date and gmt of the initial data.
!     dectim : is the time in minutes after which the solar declination
!     angle must be recalculated.
!
#ifdef MPP1
      if ( myid.eq.0 ) then
#endif              
        call open_domain(r8pt,dx,sigma)
#ifdef MPP1        
      end if
      call mpi_bcast(clat,1,mpi_real,0,mpi_comm_world,ierr)
      call mpi_bcast(plon,1,mpi_real,0,mpi_comm_world,ierr)
      call mpi_bcast(r8pt,1,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_bcast(dx,1,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_bcast(sigma,kzp1,mpi_real8,0,mpi_comm_world,ierr)
#endif 
 
!rst-fix
      mdate0 = idate0
      write (aline, *) 'param: initial date of the global '// &
                       'simulation: mdate  = ' , mdate0
      call say
      call split_idate(mdate0, myear, mmonth, mday, mhour)
      gmt = mhour
      jyear0 = myear
!
!.....find the julian day of the year and calulate dectim
!
      julday = idayofyear(mdate0)
      dectim = (1440.-gmt*60.)
 
!-----specify the constants used in the model.
!     conf   : condensation threshold.
!     qcth   : threshold for the onset of autoconversion.
!     qck1oce  : constant autoconversion rate for ocean.
!     qck1land : constant autoconversion rate for land.
!     all the other constants are used to compute the cloud
!     microphysical parameterization (ref. orville & kopp, 1977 jas).
!
      dx2 = 2.*dx
      dx4 = 4.*dx
      dx8 = 8.*dx
      dx16 = 16.*dx
      dxsq = dx*dx
      c200 = vonkar*vonkar*dx/(4.*(100.-r8pt))
      c203 = 1./dxsq
      xkhz = 1.5E-3*dxsq/dt
      xkhmax = dxsq/(64.*dt)
      akht1 = dxsq/tauht
      akht2 = dxsq/tauht
!
      conf = 1.
 
      write (aline, *) 'param: input/output parameters '
      call say
      write (aline,*) 'if true(T) create SAV files for '// &
                      'restart: ifsave = ' , ifsave  
      call say 
      write (aline,*) 'Frequency in hours to create SAV: savfrq = ' , &
                      savfrq
      call say 
      write (aline,*) 'if true (T) Output ATM files:  iftape = ' , &
                       iftape 
      call say 
      write (aline,*) 'Frequency in hours to write  ATM: tapfrq = ' , &
                      tapfrq 
      call say  
      write (aline,*) 'Frequency in hours to write  RAD: radisp = ' , &
                      radisp 
      call say
      write (aline,*) 'Frequency in hours to write  SRF: batfrq = ' , &
                      batfrq  
      call say
      write (aline,*) 'if true (T) output CHEM files:  ifchem = ' , &
                      ifchem 
      call say 
      write (aline,*) 'Frequency in hours to write CHEM: chemfrq =' , &
                      chemfrq
      call say  
      write (aline,*) 'Frequency in hours to write CLM: clmfrq = ', &
                      clmfrq
      call say
      write (aline,*) ' '
      call say
      write (aline,*) 'param: physical parameterizations '
      call say
      write (aline,'(a,i2)') ' Lateral Boundary conditions '// &
                              'scheme: iboudy = ' , iboudy
      call say  
      write (aline,'(a,i2)') ' Cumulus convection scheme: icup = ' , &
                             icup
      call say
      write  (aline,'(a,i2)') ' Grell Scheme Cumulus closure '// &
                              'scheme: igcc =' , igcc 
      call say
      write  (aline,'(a,i2)') ' Moisture scheme: ipptls = ' , ipptls 
      call say
      write  (aline,'(a,i2)') ' Ocean Flux scheme: iocnflx = ' , iocnflx
      call say
      write  (aline,'(a,i2)') ' Pressure gradient force scheme: '// &
                              'ipgf = ' , ipgf 
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
      write  (aline,'(a,f9.6)') ' Nudge value high range   =', &
                                high_nudge 
      call say
      write  (aline,'(a,f9.6)') ' Nudge value medium range =', &
                                medium_nudge 
      call say
      write  (aline,'(a,f9.6)') ' Nudge value low range    =', &
                                low_nudge 
      call say
#ifdef CLM 
       write  (aline,'(a,i2)' ) '  imask=' , imask 
      call say
#endif
      write (aline, *) ' '
      call say
      write (aline, *) 'param: model parameters '
      call say
      write (aline,'(a,f12.6)')  ' time step for radiation '// &
              'model in minutes:  radfrq = ' , radfrq 
      call say 
      write (aline,'(a,f12.6)')  ' time step for land surface '// &
              'model in seconds:  abatm  = ' , abatm 
      call say
      write (aline,'(a,f12.6)')  ' time step for LW absorption/'// &
              'emissivity in hours:  abemh  = ' , abemh 
      call say
      write (aline,'(a,f12.6)')  ' time step for atmosphere model '// &
              ' in seconds:  dt     = ' , dt
      call say
      write (aline, *) ' '
      call say
      write (aline,'(a,i2)' ) 'param: # of bottom model levels '// &
                              'with no clouds:  ncld = ' , ncld
      call say
      write (aline, *) ' '
      call say

      if ( ichem.eq.1 ) then
        chtrname = inpchtrname(1:ntr)
        chtrdpv = inpchtrdpv(1:ntr,:)
        dustbsiz = inpdustbsiz(1:nbin,:)
        chtrsol = inpchtrsol(1:ntr)
#ifdef MPP1
        do n = 1 , ntr
          call mpi_bcast(chtrname(n),5,mpi_character,0,mpi_comm_world,  &
                       & ierr)
        end do
        call mpi_bcast(chtrsol,ntr,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(chtrdpv,ntr*2,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(dustbsiz,nbin*2,mpi_real8,0,mpi_comm_world,ierr)
#endif
      end if

#ifdef MPP1
      if ( .not.ifrest ) then

        write (aline, *) 'param: Reading in DOMAIN data'
        call say

        if ( myid.eq.0 ) then
          call read_domain(mddom_io%ht,mddom_io%htsd,mddom_io%satbrt, &
                           mddom_io%xlat,mddom_io%xlong,mddom_io%msfx,&
                           mddom_io%msfd,mddom_io%f,snowc_io)
          if ( nsg.gt.1 ) then
            call read_subdomain(ht1_io,satbrt1_io,xlat1_io,xlon1_io)
          else
            do j = 1 , jx
              do i = 1 , iy
                ht1_io(1,i,j) = mddom_io%ht(i,j)*gti
                satbrt1_io(1,i,j) = mddom_io%satbrt(i,j)
              end do
            end do
          end if
          call close_domain
 
          do j = 1 , jx
            do i = 1 , iy
              inisrf_0(i,1,j) = mddom_io%ht(i,j)
              inisrf_0(i,2,j) = mddom_io%htsd(i,j)
              inisrf_0(i,3,j) = mddom_io%satbrt(i,j)
              inisrf_0(i,4,j) = mddom_io%xlat(i,j)
              inisrf_0(i,5,j) = mddom_io%xlong(i,j)
              inisrf_0(i,6,j) = mddom_io%msfx(i,j)
              inisrf_0(i,7,j) = mddom_io%msfd(i,j)
              inisrf_0(i,8,j) = mddom_io%f(i,j)
            end do
            do n = 1 , nnsg
              do i = 1 , iy
                inisrf_0(i,8+n,j) = ht1_io(n,i,j)
                inisrf_0(i,8+nnsg+n,j) = satbrt1_io(n,i,j)
                inisrf_0(i,8+nnsg*2+n,j) = snowc_io(n,i,j)
              end do
            end do
          end do
          do j = 1 , jx
            do i = 1 , iy
              mddom_io%ht(i,j) = mddom_io%ht(i,j)*gti
            end do
          end do
        end if  ! end if (myid.eq.0)
 
        call mpi_barrier(mpi_comm_world,ierr)
        call mpi_scatter(inisrf_0,iy*(nnsg*3+8)*jxp,mpi_real8,   &
                       & inisrf0, iy*(nnsg*3+8)*jxp,mpi_real8,   &
                       & 0,mpi_comm_world,ierr)
        call mpi_barrier(mpi_comm_world,ierr)
        do j = 1 , jxp
          do i = 1 , iy
            mddom%ht(i,j) = inisrf0(i,1,j)
            mddom%htsd(i,j) = inisrf0(i,2,j)
            mddom%satbrt(i,j) = inisrf0(i,3,j)
            mddom%xlat(i,j) = inisrf0(i,4,j)
            mddom%xlong(i,j) = inisrf0(i,5,j)
            mddom%msfx(i,j) = inisrf0(i,6,j)
            mddom%msfd(i,j) = inisrf0(i,7,j)
            mddom%f(i,j) = inisrf0(i,8,j)
          end do
          do n = 1 , nnsg
            do i = 1 , iy
              ht1(n,i,j) = inisrf0(i,8+n,j)
              satbrt1(n,i,j) = inisrf0(i,8+nnsg+n,j)
              snowc(n,i,j) = inisrf0(i,8+nnsg*2+n,j)
            end do
          end do
        end do
!
!------invert mapscale factors:
!
        do j = 1 , jendl
          do i = 1 , iy
            mddom%msfd(i,j) = 1./mddom%msfd(i,j)
          end do
        end do
        do j = 1 , jendl
          do i = 1 , iy
            mddom%msfx(i,j) = 1./mddom%msfx(i,j)
            mddom%ht(i,j) = mddom%ht(i,j)*gti
          end do
        end do
        if ( myid.eq.0 ) then
          print * , ' '
          print * ,                                                     &
               &'***************************************************'
          print * ,                                                     &
               &'***************************************************'
          print * ,                                                     &
               &'**** RegCM IS BEING RUN ON THE FOLLOWING GRID: ****'
          print * , '****     Map Projection: ' , iproj ,               &
               &'                ****'
          print * , '****     IY=' , iy , ' JX=' , jx , ' KX=' , kz ,   &
               &'             ****'
          print * , '****     PTOP=' , r8pt , ' DX=' , ds ,             &
               &'       ****'
          print * , '****     CLAT= ' , clat , ' CLON=' , clon ,        &
               &'    ****'
          if ( iproj.eq.'ROTMER' ) print * , '****     PLAT= ' , plat , &
                                       &' PLON=' , plon , '    ****'
          print * ,                                                     &
               &'***************************************************'
          print * , ' '
 
        end if
        call mpi_sendrecv(mddom%ht(1,jxp),iy,mpi_real8,ieast,1,        &
                        & mddom%ht(1,0),iy,mpi_real8,iwest,1,          &
                        & mpi_comm_world,mpi_status_ignore,ierr)
        call mpi_sendrecv(mddom%ht(1,1),iy,mpi_real8,iwest,2,          &
                        & mddom%ht(1,jxp+1),iy,mpi_real8,ieast,2,      &
                        & mpi_comm_world,mpi_status_ignore,ierr)
        call mpi_sendrecv(mddom%msfx(1,jxp-1),iy*2,mpi_real8,ieast,    &
                        & 1,mddom%msfx(1,-1),iy*2,mpi_real8,iwest,     &
                        & 1,mpi_comm_world,mpi_status_ignore,ierr)
        call mpi_sendrecv(mddom%msfx(1,1),iy*2,mpi_real8,iwest,2,      &
                        & mddom%msfx(1,jxp+1),iy*2,mpi_real8,ieast,    &
                        & 2,mpi_comm_world,mpi_status_ignore,ierr)
        call mpi_sendrecv(mddom%msfd(1,jxp-1),iy*2,mpi_real8,ieast,    &
                        & 1,mddom%msfd(1,-1),iy*2,mpi_real8,iwest,     &
                        & 1,mpi_comm_world,mpi_status_ignore,ierr)
        call mpi_sendrecv(mddom%msfd(1,1),iy*2,mpi_real8,iwest,2,      &
                        & mddom%msfd(1,jxp+1),iy*2,mpi_real8,ieast,    &
                        & 2,mpi_comm_world,mpi_status_ignore,ierr)
      end if
#else
      if ( .not.ifrest ) then
        write (aline, *) 'Reading in DOMAIN data'
        call say
        call read_domain(mddom%ht,mddom%htsd,mddom%satbrt, &
                         mddom%xlat,mddom%xlong,mddom%msfx,&
                         mddom%msfd,mddom%f,snowc)
        if ( nsg.gt.1 ) then
          call read_subdomain(ht1,satbrt1,xlat1,xlon1)
        else
          do j = 1 , jx
            do i = 1 , iy
              ht1(1,i,j) = mddom%ht(i,j)*gti
              satbrt1(1,i,j) = mddom%satbrt(i,j)
            end do
          end do
        end if
        call close_domain

!------invert mapscale factors:

        do j = 1 , jx
          do i = 1 , iy
            mddom%msfd(i,j) = 1./mddom%msfd(i,j)
          end do
        end do
        do j = 1 , jx
          do i = 1 , iy
            mddom%msfx(i,j) = 1./mddom%msfx(i,j)
            mddom%ht(i,j) = mddom%ht(i,j)*gti
          end do
        end do
        print * , ' '
        print * , '***************************************************'
        print * , '***************************************************'
        print * , '**** RegCM IS BEING RUN ON THE FOLLOWING GRID: ****'
        print * , '****     Map Projection: ' , iproj ,                 &
             &'                ****'
        print * , '****     IY=' , iy , ' JX=' , jx , ' KX=' , kz ,     &
             &'             ****'
        print * , '****     PTOP=' , r8pt , ' DX=' , ds ,               &
             &'       ****'
        print * , '****     CLAT= ' , clat , ' CLON=' , clon ,          &
             &'    ****'
        if ( iproj.eq.'ROTMER' ) print * , '****     PLAT= ' , plat ,   &
                                     &' PLON=' , plon , '    ****'
        print * , '***************************************************'
        print * , ' '

      end if
#endif

!
!-----compute dsigma and half sigma levels.
!
      do k = 1 , kz
        dsigma(k) = sigma(k+1) - sigma(k)
        a(k) = 0.5*(sigma(k+1)+sigma(k))
      end do
 
      do k = 1 , kz
        if ( a(k).lt.0.4 ) then
          anudg(k) = high_nudge
        else if ( a(k).lt.0.8 ) then
          anudg(k) = medium_nudge
        else
          anudg(k) = low_nudge
        end if
      end do
!
!-----specify heating profile (twght) and weighted function
!     for moisture fluxes due to convection (vqflx)
!     assume base of cloud varies as  < kbase = 5,kz >
!     top  of cloud varies as  < ktop  = 1,kbase-3 >
!     exceptions to this are treated explicitly in subroutine
!     "cupara".
!
      do kbase = 5 , kz
        do ktop = 1 , kbase - 3
          do k = 1 , kz
            twght(k,kbase,ktop) = 0.
            vqflx(k,kbase,ktop) = 0.
          end do
!
!......get twght from 1/2 level sigma values
!
          bb = dlog(a(ktop)) + dlog(a(kbase))
          cc = dlog(a(ktop))*dlog(a(kbase))
          ssum = 0.
          do k = ktop , kbase
            xx = dlog(a(k))
            twght(k,kbase,ktop) = (xx*xx) - (bb*xx) + cc
            ssum = ssum + twght(k,kbase,ktop)*dsigma(k)
          end do
          do k = ktop , kbase
            twght(k,kbase,ktop) = twght(k,kbase,ktop)/ssum
          end do
!
!......get vqflx from  d(w*q) / dsigma on full levels
!         do computations in p to avoid sigma=0. discontinuity
!
          xtop = dlog((100.-r8pt)*sigma(ktop)+r8pt)
          xbot = dlog((100.-r8pt)*sigma(kbase+1)+r8pt)
          bb = xtop + xbot
          cc = xtop*xbot
          vqmax = 0.
          ssum = 0.
          xx = xtop
          yy = xbot
          wk = (xx*xx) - (bb*xx) + cc
          qk = -((yy*yy)-(bb*yy)+cc)
          do k = ktop , kbase
            xx = dlog((100.-r8pt)*sigma(k+1)+r8pt)
            yy = dlog((100.-r8pt)                                       &
               & *(sigma(ktop)+sigma(kbase+1)-sigma(k+1))+r8pt)
            wkp1 = (xx*xx) - (bb*xx) + cc
            qkp1 = -((yy*yy)-(bb*yy)+cc)
            vqflx(k,kbase,ktop) = -((wkp1*qkp1)-(wk*qk))/dsigma(k)
            ssum = ssum + vqflx(k,kbase,ktop)
            if ( dabs(vqflx(k,kbase,ktop)).gt.vqmax )                   &
               & vqmax = dabs(vqflx(k,kbase,ktop))
            wk = wkp1
            qk = qkp1
          end do
          do k = ktop , kbase
            vqflx(k,kbase,ktop) = vqflx(k,kbase,ktop)*vqrang/vqmax
          end do
!
        end do
      end do
!
!----calculate max no of pbl levels: kt=k at highest allowed pbl level
!-----1. caluclate sigma level at 700mb, assuming 600mb highest
!-----   allowed pbl, and 1013mb as standard surface pressure. (sigtbl)
!-----2. find first model sigma level above sigtbl.
!
      sigtbl = (70.-r8pt)/(101.3-r8pt)
      kt = 1
      do k = kz , 1 , -1
        delsig = a(k) - sigtbl
        if ( delsig.le.0. ) then
          kt = k
          exit
        end if
      end do
      write (aline, '(a,i3)') 'param: Index of highest allowed pbl:'// &
                              '  kt = ' , kt
      call say
      write (aline, *) ' '
      call say
!
      if ( ipptls.eq.1 ) then
        write (aline, '(a,f11.6,a,f11.6)')                              &
               & 'AUTO-CONVERSION RATE:  LAND=' , qck1land ,    &
               & '                      OCEAN=' , qck1oce
        call say
        write (aline, *) 'RELATIVE HUMIDITY THRESHOLDS:  LAND=' ,       &
               & rh0land , '                              OCEAN=' ,     &
               & rh0oce
        call say
        write (aline, *) 'GULTEPE FACTORS:  LAND=' , gulland ,          &
               &'                 OCEAN=' , guloce
        call say
        write (aline, *) 'MAXIMUM CLOUD COVER FOR RADIATION: ' , fcmax
        call say
        write (aline, *) 'MAXIMUM RELATIVE HUMIDITY: ' , rhmax
        call say
        write (aline, *) 'rh0 temperature threshold: ' , tc0
        call say
        if ( cevap.le.0.0 ) then
          write (aline, *) 'RAINDROP EVAPORATION NOT INCLUDED'
          call say
        end if
        if ( caccr.le.0.0 ) then
          write (aline, *) 'RAINDROP ACCRETION NOT INCLUDED'
          call say
        end if
        write (aline, *) 'Raindrop Evaporation Rate' , cevap
        call say
        write (aline, *) 'Raindrop Accretion Rate' , caccr
        call say
      end if
 
      write (aline, *) ' '
      call say
 
      if ( icup.eq.1 ) then
        write (aline, *) '*********************************'
        call say
        write (aline, *) '***** Anthes-Kuo Convection *****'
        call say
        write (aline, *) '*********************************'
        call say
      else if ( icup.eq.2 ) then
        kbmax = kz
        do k = 1 , kz - 1
          if ( a(k).le.skbmax ) kbmax = kz - k
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
        if ( igcc.eq.1 ) then
          write (aline, *)                                              &
               &'   Arakawa-Schubert (1974) Closure Assumption'
          call say
        else if ( igcc.eq.2 ) then
          write (aline, *)                                              &
               &'   Fritsch-Chappell (1980) Closure Assumption'
          call say
          write (aline, *) '     ABE removal timescale: dtauc=' , dtauc
          call say
        else
        end if
        write (aline, *) '*********************************'
        call say
#ifdef MPP1
        do j = 1 , jendx
#else
#ifdef BAND
        do j = 1 , jx
#else
        do j = 1 , jxm1
#endif
#endif
          do i = 1 , iym1
            shrmax2d(i,j) = shrmax
            shrmin2d(i,j) = shrmin
            edtmax2d(i,j) = edtmax
            edtmin2d(i,j) = edtmin
            edtmaxo2d(i,j) = edtmaxo
            edtmino2d(i,j) = edtmino
            edtmaxx2d(i,j) = edtmaxx
            edtminx2d(i,j) = edtminx
            pbcmax2d(i,j) = pbcmax
            mincld2d(i,j) = mincld
            kbmax2d(i,j) = kbmax
            htmax2d(i,j) = htmax
            htmin2d(i,j) = htmin
            dtauc2d(i,j) = dtauc*60.
          end do
        end do
      else if ( icup.eq.3 ) then
        write (aline,*) ' The Betts-Miller Convection scheme is not' ,  &
                       &' properly implemented'
        call say
        call fatal(__FILE__,__LINE__,'BETTS-MILLER NOT WORKING')
        call allocate_mod_cu_bm(lmpi)
      else if ( icup.eq.4 ) then
        cllwcv = 0.5E-4    ! Cloud liquid water content for convective precip.
        clfrcvmax = 0.25   ! Max cloud fractional cover for convective precip.
        minorig = kz
        do k = 1 , kz
          if ( a(k).le.minsig ) minorig = kz - k
        end do
        write (aline, *) ' '
        call say
        write (aline, *)                                                &
              &'EMANUEL (1991) CONVECTION V4.3C (20 May, 2002)'
        call say
        write (aline, *)                                                &
              &'  MIN CONVECTION ORIGIN (minsig/orig): ' , minsig ,     &
              & minorig
        call say
        write (aline, *)'  AUTOCONVERSION THERSHOLD (elcrit): ' , elcrit
        call say
        write (aline, *)'  AUTOCONVERSION THRESHOLD TO ZERO (tlcrit): ',&
              & tlcrit
        call say
        write (aline, *)'  ENTRAINMENT COEFFICIENT (entp): ' , entp
        call say
        write (aline, *)                                                &
              &'  FRACTIONAL AREA OF UNSATURATED DNDRAFT (sigd): '      &
              & , sigd
        call say
        write (aline, *)'  PRECIP FRACTION OUTSIDE OF CLOUD (sigs): ' , &
              & sigs
        call say
        write (aline, *)'  FALL SPEED OF RAIN (omtrain): ' , omtrain
        call say
        write (aline, *)'  FALL SPEED OF SNOW (omtsnow): ' , omtsnow
        call say
        write (aline, *)'  RAIN EVAPORATION COEFFICIENT (coeffr): ' ,   &
              & coeffr
        call say
        write (aline, *)'  SNOW EVAPORATION COEFFICIENT (coeffs): ' ,   &
              & coeffs
        call say
        write (aline, *)                                                &
              &'  CONVECTIVE MOMENTUM TRANSPORT COEFFICIENT (cu): '     &
              & , cu
        call say
        write (aline, *)'  DOWNDRAFT VELOCITY SCALE (betae): ' , betae
        call say
        write (aline, *)                                                &
              &'  MAX NEGATIVE PERTURBATION BELOW LFC (dtmax): ' ,      &
              & dtmax
        call say
        write (aline, *)'  QUASI-EQUILIBRIUM APPROACH RATE (alphae): ' ,&
              & alphae
        call say
        write (aline, *)'  QUASI-EQUILIBRIUM APPROACH RATE (damp): ' ,  &
              & damp
        call say
        write (aline, *) ' '
        call say
        write (aline, *) ' '
        call say
      else
      end if
 
!     Convective Cloud Cover
      afracl = 0.3 ! frac. cover for conv. precip. when dx=dxlarg
      afracs = 1.0 !   "     "    "    "      "     "   dx=dxsmal
      dlargc = 200.0
      dsmalc = 10.0
      dxtemc = dmin1(dmax1(dx,dsmalc),dlargc)
      clfrcv = afracl + (afracs-afracl)                                 &
             & *((dlargc-dxtemc)/(dlargc-dsmalc))**2
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
      twt(1,1) = 0.
      twt(1,2) = 0.
      qcon(1) = 0.
      do k = 2 , kz
        twt(k,1) = (sigma(k)-a(k-1))/(a(k)-a(k-1))
        twt(k,2) = 1. - twt(k,1)
        qcon(k) = (sigma(k)-a(k))/(a(k-1)-a(k))
      end do
 
      chibot = 450.
      ptmb = 10.*r8pt
      pz = a(1)*(1000.-ptmb) + ptmb
      if ( pz.gt.chibot ) call fatal(__FILE__,__LINE__,                 &
                                    &'VERTICAL INTERPOLATION ERROR')
      do k = 1 , kz
        pk = a(k)*(1000.-ptmb) + ptmb
        if ( pk.le.chibot ) kchi = k
      end do
 
!
!-----compute the k level under which the maximum equivalent potential
!     temperature will be regarded as the origin of air parcel that
!     produces cloud (used in the cumulus parameterization scheme).
!
      sig700 = (70.-r8pt)/(100.-r8pt)
      do k = 1 , kz
        k700 = k
        if ( sig700.le.sigma(k+1) .and. sig700.gt.sigma(k) ) exit
      end do
!
!-----specify the coefficients for sponge boundary conditions.
!
      ispgd = nspgd - 1
      ispgx = nspgx - 1
!.....for dot point variables:
      if ( iboudy.eq.4 ) then
        wgtd(1) = 0.
        wgtd(2) = 0.2
        wgtd(3) = 0.55
        wgtd(4) = 0.8
        wgtd(5) = 0.95
        do k = 4 , nspgx
          wgtd(k) = 1.
        end do
!.....for cross point variables:
        wgtx(1) = 0.
        wgtx(2) = 0.4
        wgtx(3) = 0.7
        wgtx(4) = 0.9
        do k = 5 , nspgx
          wgtx(k) = 1.
        end do
      end if
!
!-----specify the coefficients for nudging boundary conditions:
!
!.....for large domain:
      if ( iboudy.eq.1 .or. iboudy.eq.5 ) then
        fnudge = 0.1/dt2
        gnudge = (dxsq/dt)/50.
      end if
      if ( icup.eq.3 ) call lutbl(r8pt)
!
!-----print out the parameters specified in the model.
!
#ifdef MPP1
      if ( myid.eq.0 ) then
        if ( ibltyp.eq.0 ) print 99002
        print 99003 , julday , gmt , ntrad
!
        if ( iboudy.eq.0 ) print 99004
        if ( iboudy.eq.1 ) print 99005 , fnudge , gnudge
        if ( iboudy.eq.2 ) print 99007
        if ( iboudy.eq.3 ) print 99008
        if ( iboudy.eq.4 ) print 99009
        if ( iboudy.eq.5 ) print 99006 , fnudge , gnudge
 
        print 99010
!
        do k = 1 , kz
          print 99011 , k , sigma(k) , a(k) , dsigma(k) , twt(k,1) ,    &
              & twt(k,2) , qcon(k)
        end do
        print 99012 , kzp1 , sigma(kzp1)
        print 99014 , dt
        print 99015 , dx
        print 99016 , jx , iy
        print 99017 , kz
        print 99018 , xkhz
        print 99019 , xkhmax
      end if
#else
      if ( ibltyp.eq.0 ) print 99002
      print 99003 , julday , gmt , ntrad
!
      if ( iboudy.eq.0 ) print 99004
      if ( iboudy.eq.1 ) print 99005 , fnudge , gnudge
      if ( iboudy.eq.2 ) print 99007
      if ( iboudy.eq.3 ) print 99008
      if ( iboudy.eq.4 ) print 99009
      if ( iboudy.eq.5 ) print 99006 , fnudge , gnudge

      print 99010
!
      do k = 1 , kz
        print 99011 , k , sigma(k) , a(k) , dsigma(k) , twt(k,1) ,      &
            & twt(k,2) , qcon(k)
      end do
      print 99012 , kzp1 , sigma(kzp1)
      print 99014 , dt
      print 99015 , dx
      print 99016 , jx , iy
      print 99017 , kz
      print 99018 , xkhz
      print 99019 , xkhmax
#endif
99002 format (/'   frictionless and insulated for the lower boundary.')
99003 format (                                                          &
     &'     the surface energy budget is used to calculate the ground te&
     &mperature.   julday = ',i3,'   gmt = ',f4.1/10x,                  &
     &'the radiation is computed every ',i4,' time steps.')
!1100 format('     heat and moisture fluxes from the ground are turned o
!     1ff.')
99004 format (/'   the lateral boundary conditions are fixed.')
99005 format (/'   relaxation boundary conditions (linear method)',     &
             &' are used.  fnudge = ',e15.5,'  gnudge = ',e15.5)
99006 format (/'   relaxation boudnary conditions (exponential method)',&
             &' are used. fnudge = ',e15.5,'  gnudge = ',e15.5)
99007 format (/'   time dependent boundary conditions are used.')
99008 format (/'   inflow/outflow boundary conditions are used.')
99009 format (/'   sponge boundary conditions are used.')
99010 format ('0 k',4x,'sigma(k)',3x,'  a(k)',5x,'dsigma(k)',4x,        &
             &'twt(k,1)',5x,'twt(k,2)',5x,'qcon(k)'/)
99011 format (1x,i2,5x,f6.4,5x,f6.4,5x,f6.4,5x,f8.4,5x,f8.4,5x,f8.4)
99012 format (1x,i2,5x,f6.4)
99014 format (' time step = ',f7.2,' seconds')
99015 format (' dx = ',f7.0,' meters')
99016 format (' grid points (x,y) = (',i4,',',i4,')')
99017 format (' number of levels = ',i2)
99018 format (' constant hor. diff. coef. = ',e12.5,' m*m/s')
99019 format (' maximum  hor. diff. coef. = ',e12.5,' m*m/s')
!
      end subroutine param
!
      end module mod_param

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
      subroutine param

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine defines various model parameters.               c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      use mod_regcm_param
      use mod_param1
      use mod_param2
      use mod_param3 , only : wgtx , sigma , dsigma , a , anudg , twt , &
                   & qcon , wgtd , akht1 , akht2 , kt , kxout , ncld ,  &
                   & ptop , ptop4 , kchi , k700 , jxsex , ispgx , ispgd
      use mod_iunits
      use mod_pmoist
      use mod_bats
      use mod_main
      use mod_trachem
      use mod_convect
      use mod_date , only : mdate0 , ntimax , julday , gmt , deltmx ,   &
                   & dectim , nnnnnn , nstart , idate0 , idate1 ,       &
                   & idate2 , nstrt0 , nyear , nnnend , nmonth ,        &
                   & ndate0 , nnnchk , lyear , lmonth , lday , lhour ,  &
                   & ldatez , idatex , jyear , jyear0 , jyearr , ntime ,&
                   & ktau , ktaur , xtime
      use mod_message
      use mod_grads
      use mod_constants , only : mathpi , gti , rgti , rgas , vonkar ,  &
                               & cpd , tauht
#ifdef MPP1
      use mod_mppio
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
      real(8) :: afracl , afracs , bb , cc , chibot , daymax , delsig , &
               & dlargc , dsmalc , dxtemc , pk , ptmb , pz , qk ,       &
               & qkp1 , sig700 , sigtbl , ssum , vqmax , vqrang , wk ,  &
               & wkp1 , xbot , xtop , xx , yy
      real(8) , dimension(nsplit) :: dtsplit
      character(7) :: finm
      real(4) :: grdfac
      integer :: i , ibigend , ierr1 , igrads , ii , j , jj , k ,       &
               & kbase , ktop , kzz , m , mdate1 , mday , mmon , my1 ,  &
               & my2 , my3 , myear , n , ns , nxx , nyy
      integer , dimension(12) :: mmd
      real(4) , dimension(kxp1) :: sp1d
      real(4) , dimension(ix,jx) :: sp2d
      real(4) , dimension(ix*nsg,jx*nsg) :: sp2d1
#ifdef MPP1
      integer :: ierr
#endif
!
!
!----------------------------------------------------------------------
!-----vqrang is the range limit on vqflx.
!.....qdcrit is the precipitation threshold for moisture convergence.
!
      data vqrang/5.0E-4/
      data (mmd(i),i=1,12)/31 , 28 , 31 , 30 , 31 , 30 , 31 , 31 , 30 , &
          & 31 , 30 , 31/
 
! afracl   - frac. cover for conv. precip. when dx=dxlarg
! afracs   -   "     "    "    "      "     "   dx=dxsmal
!
!----------------------------------------------------------------------
!-----namelist:
!
      namelist /restartparam/ ifrest , idate0 , idate1 , idate2 , nslice
 
      namelist /timeparam/ radfrq , abatm , abemh , dt , ibdyfrq
 
!chem2
#ifdef CLM
      namelist /outparam/ ifsave , savfrq , iftape , tapfrq , ifprt ,   &
      & prtfrq , kxout , jxsex , ifrad , radisp , ifbat , ifsub ,       &
      & batfrq , iotyp , ibintyp , ifchem , chemfrq , clmfrq
#else
      namelist /outparam/ ifsave , savfrq , iftape , tapfrq , ifprt ,   &
      & prtfrq , kxout , jxsex , ifrad , radisp , ifbat , ifsub ,       &
      & batfrq , iotyp , ibintyp , ifchem , chemfrq
#endif
!
!chem2
#ifdef CLM
      namelist /physicsparam/ ibltyp , iboudy , icup , igcc , ipgf ,    &
      & iemiss , lakemod , ipptls , iocnflx , ichem , imask
#else
      namelist /physicsparam/ ibltyp , iboudy , icup , igcc , ipgf ,    &
      & iemiss , lakemod , ipptls , iocnflx , ichem
#endif
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
      & ichcumtra , idirect , mixtype , chtrname , chtrsol , chtrdpv ,  &
      & dustbsiz
!chem2_
!
!----------------------------------------------------------------------
!-----specify unit numbers for input/output.
!     units 10,101,14  are input.
!     iutin : input initial data for large domain.
!     iutbc : input slab temperature for large domain ,
!     boundary values and tendencies for large doamin.
!     iutrs : input saved file for large domain for restart from
!     previous forecast.
!     iutdat: output for dataflow, if iftape=.true.
!     iutsav: output saved file for restart, if ifsave=.true. or
!
      qdcrit = 3.0E-7                      ! real*4 for compact input
      iutin = 10
      iutin1 = 11
      iutbc = 101
!_sgi iutbc = 71
      iutrs = 14
!**   initialize output file unit number
      iutdat = 50
!**   initialize save file unit number
      iutsav = 52
!**   initialize bats file unit number
      iutbat = 54
      iutsub = 55
!**   initialize radiation file unit number
      iutrad = 56
!**   initalize lake file unit number
      iutlak = 58
!**   initialize chemistry file number
!chem2
!     emission file
      iutopt = 12
      iutchsrc = 13
!     out file
      iutchem = 60
!chem2_
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
!     ifprt  : whether you want printer output or not.
!     = 0 ; no
!     = 1 ; yes
!
!     prtfrq : if ifprt=1, specify the output interval in mimutes.
!
!     iotyp  : Type of output files,
!     1=direct access (GrADS); 2=sequential w/ time listing
!
!     ibintyp : Type of binary output for GrADS ctl files
!     1=big endian; 2=little endian
!
!     maschk : specify the frequency in time steps, the mass-
!     conservation information will be printed out.
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
      nslice = 4      ! # of days for next model run (used with restart)
      ! note: beginning/end forecast time set in restart.mm4
 
!------namelist timeparam:
!
      radfrq = 30.    ! time interval in min solar rad caluclated
      abatm = 600.    ! time interval at which bats is called (secs)
      abemh = 12.     ! time interval absorption-emission calculated (hours)
      dt = 200.       ! time step in seconds
      ibdyfrq = 6     ! boundary condition interval (hours)
 
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
      ifprt = .true.
      prtfrq = 12.
      kxout = kx
      jxsex = 25
      iotyp = 1
      ibintyp = 1
      maschk = 10           ! * defined below
 
!chem2
      ifchem = .false.
      chemfrq = 6.0   ! time interval for disposeing chem output (hrs)
!chem2_
 
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
      clmfrq = 0
      imask = 2
#endif

!---------------------------------------------------------------------
!-----read in namelist variables:
!
#ifdef MPP1
      if ( myid.eq.0 ) then
        open (2,file='regcm.in',status='old')
        read (2,restartparam)
        print * , 'RESTARTPARAM READ IN'
        read (2,timeparam)
        print * , 'TIMEPARAM READ IN'
        read (2,outparam)
        print * , 'OUTPARAM READ IN'
        read (2,physicsparam)
        print * , 'PHYSICSPARAM READ IN'
        if ( ipptls.eq.1 ) then
          read (2,subexparam)
          print * , 'SUBEXPARAM READ IN'
        end if
        if ( icup.eq.2 ) then
          read (2,grellparam)
          print * , 'GRELLPARAM READ IN'
        else if ( icup.eq.4 ) then
          read (2,emanparam)
          print * , 'EMANPARAM READ IN'
        else
        end if
        if ( ichem.eq.1 ) then
          read (2,chemparam)
          print * , 'CHEMPARAM READ IN'
        end if
        close (2)
      end if
      call mpi_bcast(ifrest,1,mpi_logical,0,mpi_comm_world,ierr)
      call mpi_bcast(idate0,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(idate1,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(idate2,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(nslice,1,mpi_integer,0,mpi_comm_world,ierr)
 
      call mpi_bcast(radfrq,1,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_bcast(abemh,1,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_bcast(abatm,1,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_bcast(dt,1,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_bcast(ibdyfrq,1,mpi_integer,0,mpi_comm_world,ierr)
 
      call mpi_bcast(ifsave,1,mpi_logical,0,mpi_comm_world,ierr)
      call mpi_bcast(savfrq,1,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_bcast(iftape,1,mpi_logical,0,mpi_comm_world,ierr)
      call mpi_bcast(tapfrq,1,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_bcast(ifrad,1,mpi_logical,0,mpi_comm_world,ierr)
      call mpi_bcast(radisp,1,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_bcast(ifbat,1,mpi_logical,0,mpi_comm_world,ierr)
      call mpi_bcast(ifsub,1,mpi_logical,0,mpi_comm_world,ierr)
      call mpi_bcast(batfrq,1,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_bcast(ifprt,1,mpi_logical,0,mpi_comm_world,ierr)
      call mpi_bcast(prtfrq,1,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_bcast(kxout,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(jxsex,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(iotyp,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(ibintyp,1,mpi_integer,0,mpi_comm_world,ierr)
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
 
#ifdef CLM
      call mpi_bcast(clmfrq,1,mpi_integer,0,mpi_comm_world,ierr)
      call mpi_bcast(imask,1,mpi_integer,0,mpi_comm_world,ierr)
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
        do n = 1 , ntr
          call mpi_bcast(chtrname(n),5,mpi_character,0,mpi_comm_world,  &
                       & ierr)
        end do
        call mpi_bcast(chtrsol,ntr,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(chtrdpv,ntr*2,mpi_real8,0,mpi_comm_world,ierr)
        call mpi_bcast(dustbsiz,nbin*2,mpi_real8,0,mpi_comm_world,ierr)
      end if
#else
      read (*,restartparam)
      print * , 'RESTARTPARAM READ IN'
      read (*,timeparam)
      print * , 'TIMEPARAM READ IN'
      read (*,outparam)
      print * , 'OUTPARAM READ IN'
      read (*,physicsparam)
      print * , 'PHYSICSPARAM READ IN'
      if ( ipptls.eq.1 ) then
        read (*,subexparam)
        print * , 'SUBEXPARAM READ IN'
      end if
      if ( icup.eq.2 ) then
        read (*,grellparam)
        print * , 'GRELLPARAM READ IN'
      else if ( icup.eq.4 ) then
        read (*,emanparam)
        print * , 'EMANPARAM READ IN'
      else
      end if
      if ( ichem.eq.1 ) then
        read (*,chemparam)
        print * , 'CHEMPARAM READ IN'
      end if
#endif
 
      if ( ichem.eq.0 ) ifchem = .false.
!as
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
 
!-----reset the options/calculate variables using namelist info:
!
      ndate0 = idate1
      nsavfrq = nint(3600.*savfrq)
      ntapfrq = nint(3600.*tapfrq)
      nprtfrq = nint(3600.*prtfrq)
      ktau = 0
      ktaur = 0
      xtime = 0.
      ntime = 0
      dtsplit(2) = dt/2.
      dtsplit(1) = dt/4.
      do ns = 1 , nsplit
        dtau(ns) = dtsplit(ns)
      end do
      write (aline, *) ' dtau = ' , dtau
      call say
      dt0 = dt      !store original dt
      maschk = nint(prtfrq*3600./dt)
                                   !convert prtfrq to #of time steps
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
      call initdate
      call finddate(nstrt0,idate0)
      call finddate(nstart,idate1)
      call finddate(nnnend,idate2)
      nnnchk = nstart
 
!     ktaur = nint((NSTART-NSTRT0)*ibdyfrq*60/dtmin)
      ntimax = (nnnend-nstrt0)*ibdyfrq*60
      write (aline, *) 'IDATE1, IDATE2, dtmin, ktaur = ' ,              &
                           & idate1 , idate2 , dtmin , ktaur
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
        print * , 'READING HEADER FILE'
        write (finm,99001) iutin
        if ( nsg.gt.1 ) open (iutin1,file='fort.11',form='unformatted', &
                            & status='old',access='direct',             &
                            & recl=ix*jx*nnsg*ibyte)
        open (iutin,file=finm,form='unformatted',status='old',          &
             &access='direct',recl=ix*jx*ibyte)
        read (iutin,rec=1,iostat=ierr1) nyy , nxx , kzz , dxsp , clat , &
                                      & clon , plat , plon , grdfac ,   &
                                      & proj , sp1d , ptsp , igrads ,   &
                                      & ibigend , truelatl , truelath
        print * , 'DIMS' , nyy , nxx , kzz
        print * , 'DOMAIN' , dxsp , clat , clon , plat , plon , grdfac
        print * , 'PROJ' , proj
        print * , 'SIGMA' , sp1d
        print * , 'PTOP' , ptsp
        print * , 'OUTPUT' , igrads , ibigend
        ptop = ptsp
        dx = dxsp
        if ( nyy.ne.ix .or. nxx.ne.jx .or. kzz.ne.kx ) then
          write (aline,*) '  SET IN regcm.param:  IX=' , ix , ' JX=' ,  &
                        & jx , ' KX=' , kx
          call say
          write (aline,*) '  SET IN TERRAIN: NYY=' , nyy , ' NXX=' ,    &
                        & nxx , ' NZZ=' , kzz
          call say
          write (aline,*) '  Also check ibyte in regcm.param: ibyte = ' &
                        & , ibyte
          call fatal(__FILE__,__LINE__,                                 &
                    &'IMPROPER DIMENSION SPECIFICATION')
        end if
        do k = 1 , kxp1
          sigma(k) = dble(sp1d(k))
        end do
      end if
      call mpi_bcast(clat,1,mpi_real,0,mpi_comm_world,ierr)
      call mpi_bcast(plon,1,mpi_real,0,mpi_comm_world,ierr)
      call mpi_bcast(ptop,1,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_bcast(dx,1,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_bcast(sigma,kxp1,mpi_real8,0,mpi_comm_world,ierr)
#else
      print * , 'READING HEADER FILE'
      write (finm,99001) iutin
      if ( nsg.gt.1 ) open (iutin1,file='fort.11',form='unformatted',   &
                          & status='old',access='direct',               &
                          & recl=ix*jx*nnsg*ibyte)
      open (iutin,file=finm,form='unformatted',status='old',            &
           &access='direct',recl=ix*jx*ibyte)
      read (iutin,rec=1,iostat=ierr1) nyy , nxx , kzz , dxsp , clat ,   &
                                    & clon , plat , plon , grdfac ,     &
                                    & proj , sp1d , ptsp , igrads ,     &
                                    & ibigend , truelatl , truelath
      print * , 'DIMS' , nyy , nxx , kzz
      print * , 'DOMAIN' , dxsp , clat , clon , plat , plon , grdfac
      print * , 'PROJ' , proj
      print * , 'SIGMA' , sp1d
      print * , 'PTOP' , ptsp
      print * , 'OUTPUT' , igrads , ibigend
      ptop = ptsp
      dx = dxsp
      if ( nyy.ne.ix .or. nxx.ne.jx .or. kzz.ne.kx ) then
        write (aline,*) '  SET IN regcm.param:  IX=' , ix , ' JX=' ,    &
                      & jx , ' KX=' , kx
        call say
        write (aline,*) '  SET IN TERRAIN: NYY=' , nyy , ' NXX=' , nxx ,&
                       &' NZZ=' , kzz
        call say
        write (aline,*) '  Also check ibyte in regcm.param: ibyte = ' , &
                      & ibyte
        call fatal(__FILE__,__LINE__,'IMPROPER DIMENSION SPECIFICATION')
      end if
      do k = 1 , kxp1
        sigma(k) = dble(sp1d(k))
      end do
#endif
 
!rst-fix
      mdate0 = idate0
      write (aline, *) '***** mdate = ' , mdate0
      call say
      myear = mdate0/1000000
      nyear = idate1/1000000
      nmonth = (idate1-nyear*1000000)/10000
      mdate1 = mdate0 - myear*1000000
      mmon = mdate1/10000
      mday = (mdate1-mmon*10000)/100
      gmt = mdate1 - mmon*10000 - mday*100
      jyear0 = myear
!rst-fix_
!
!.....find the year is a leap year or not (assume the initial data
!     is within this century.):
      if ( myear.lt.100 ) myear = 1900 + myear
      my1 = mod(myear,4)
      my2 = mod(myear,100)
      my3 = mod(myear,400)
      if ( my1.eq.0 .and. my2.ne.0 .or. my3.eq.0 ) mmd(2) = 29
      julday = mday
      do m = 1 , mmon - 1
        julday = julday + mmd(m)
      end do
      dectim = (1440.-gmt*60.)
 
!-----specify the constants used in the model.
!     conf   : condensation threshold.
!     qcth   : threshold for the onset of autoconversion.
!     qck1oce  : constant autoconversion rate for ocean.
!     qck1land : constant autoconversion rate for land.
!     all the other constants are used to compute the cloud
!     microphysical parameterization (ref. orville & kopp, 1977 jas).
!
      ptop4 = 4.*ptop
      dx2 = 2.*dx
      dx4 = 4.*dx
      dx8 = 8.*dx
      dx16 = 16.*dx
      dxsq = dx*dx
      c200 = vonkar*vonkar*dx/(4.*(100.-ptop))
      c201 = (100.-ptop)/dxsq
      c203 = 1./dxsq
      xkhz = 1.5E-3*dxsq/dt
      xkhmax = dxsq/(64.*dt)
      akht1 = dxsq/tauht
      akht2 = dxsq/tauht
!
      conf = 1.
 
      write (aline, *) ' input/output parameters '
      call say
      write (aline, *) ' ifsave = ' , ifsave , ' savfrq = ' , savfrq ,  &
             &' iftape = ' , iftape , ' tapfrq = ' , tapfrq ,           &
             &' ifprt  = ' , ifprt , ' prtfrq = ' , prtfrq ,            &
             &' kxout  = ' , kxout , ' jxsex  = ' , jxsex ,             &
             &' radisp = ' , radisp , ' batfrq = ' , batfrq ,           &
             &' nslice = ' , nslice , ' ifchem = ' , ifchem ,           &
             &' chemfrq =' , chemfrq
      call say
      write (aline, *) ' '
      call say
      write (aline, *) ' physical parameterizations '
      call say
      write (aline, *) ' iboudy = ' , iboudy , ' icup = ' , icup ,      &
            & ' igcc =' , igcc , ' ipptls = ' , ipptls , ' iocnflx = ' ,&
            & iocnflx , ' ipgf = ' , ipgf , 'iemiss = ' , iemiss ,      &
            &' lakemod = ' , lakemod , ' ichem =' , ichem , ' imask = ',&
            & imask
      call say
      write (aline, *) ' '
      call say
      write (aline, *) ' model parameters '
      call say
      write (aline, *) ' radfrq = ' , radfrq , ' abatm = ' , abatm ,    &
            &' abemh = ' , abemh , ' dt = ' , dt
      call say
      write (aline, *) ' '
      call say
      write (aline, *) ' ncld = ' , ncld
      call say
      write (aline, *) ' '
      call say
!
#ifdef MPP1
      if ( .not.ifrest ) then
        if ( myid.eq.0 ) then
          print * , 'HT'
          read (iutin,rec=2,iostat=ierr1) ((sp2d(i,j),j=1,jx),i=1,ix)
          do j = 1 , jx
            do i = 1 , ix
              ht_io(i,j) = dble(sp2d(i,j))
            end do
          end do
          if ( nsg.gt.1 ) then
            read (iutin1,rec=2,iostat=ierr1)                            &
                & ((sp2d1(i,j),j=1,jx*nsg),i=1,ix*nsg)
            do j = 1 , jx*nsg
              do i = 1 , ix*nsg
                jj = mod(j,nsg)
                if ( jj.eq.0 ) jj = nsg
                ii = mod(i,nsg)
                if ( ii.eq.0 ) ii = nsg
                k = (jj-1)*nsg + ii
                jj = (j+nsg-1)/nsg
                ii = (i+nsg-1)/nsg
                ht1_io(k,ii,jj) = sp2d1(i,j)*gti
              end do
            end do
          else
            do j = 1 , jx
              do i = 1 , ix
                ht1_io(1,i,j) = sp2d(i,j)*gti
              end do
            end do
          end if
 
          print * , 'HTSD'
          read (iutin,rec=3,iostat=ierr1) ((sp2d(i,j),j=1,jx),i=1,ix)
          do j = 1 , jx
            do i = 1 , ix
              htsd_io(i,j) = dble(sp2d(i,j))
            end do
          end do
          print * , 'SATBRT'
          read (iutin,rec=4,iostat=ierr1) ((sp2d(i,j),j=1,jx),i=1,ix)
          do j = 1 , jx
            do i = 1 , ix
              satbrt_io(i,j) = dble(sp2d(i,j))
            end do
          end do
          if ( nsg.gt.1 ) then
            read (iutin1,rec=4,iostat=ierr1)                            &
                & ((sp2d1(i,j),j=1,jx*nsg),i=1,ix*nsg)
            do j = 1 , jx*nsg
              do i = 1 , ix*nsg
                jj = mod(j,nsg)
                if ( jj.eq.0 ) jj = nsg
                ii = mod(i,nsg)
                if ( ii.eq.0 ) ii = nsg
                k = (jj-1)*nsg + ii
                jj = (j+nsg-1)/nsg
                ii = (i+nsg-1)/nsg
                satbrt1_io(k,ii,jj) = sp2d1(i,j)
              end do
            end do
          else
            do j = 1 , jx
              do i = 1 , ix
                satbrt1_io(1,i,j) = satbrt_io(i,j)
              end do
            end do
          end if
          print * , 'XLAT'
          read (iutin,rec=5,iostat=ierr1) ((sp2d(i,j),j=1,jx),i=1,ix)
          do j = 1 , jx
            do i = 1 , ix
              xlat_io(i,j) = dble(sp2d(i,j))
            end do
          end do
          print * , 'XLONG'
          read (iutin,rec=6,iostat=ierr1) ((sp2d(i,j),j=1,jx),i=1,ix)
          do j = 1 , jx
            do i = 1 , ix
              xlong_io(i,j) = dble(sp2d(i,j))
            end do
          end do
          print * , 'MSFX'
          read (iutin,rec=9,iostat=ierr1) ((sp2d(i,j),j=1,jx),i=1,ix)
          do j = 1 , jx
            do i = 1 , ix
              msfx_io(i,j) = dble(sp2d(i,j))
            end do
          end do
          print * , 'MSFD'
          read (iutin,rec=10,iostat=ierr1) ((sp2d(i,j),j=1,jx),i=1,ix)
          do j = 1 , jx
            do i = 1 , ix
              msfd_io(i,j) = dble(sp2d(i,j))
            end do
          end do
          print * , 'F'
          read (iutin,rec=11,iostat=ierr1) ((sp2d(i,j),j=1,jx),i=1,ix)
          do j = 1 , jx
            do i = 1 , ix
              f_io(i,j) = dble(sp2d(i,j))
            end do
          end do
          print * , 'SNOWC'
          read (iutin,rec=12,iostat=ierr1) ((sp2d(i,j),j=1,jx),i=1,ix)
          do j = 1 , jx
            do i = 1 , ix
              do n = 1 , nnsg
                snowc_io(n,i,j) = dble(sp2d(i,j))
              end do
            end do
          end do
          if ( ierr1.ne.0 ) then
            write (aline,*) '  Check ibyte in parameter.inc: ibyte = ' ,&
                          & ibyte
            call fatal(__FILE__,__LINE__,'REACHED EOF')
          end if
          do j = 1 , jx
            do i = 1 , ix
              inisrf_0(i,1,j) = ht_io(i,j)
              inisrf_0(i,2,j) = htsd_io(i,j)
              inisrf_0(i,3,j) = satbrt_io(i,j)
              inisrf_0(i,4,j) = xlat_io(i,j)
              inisrf_0(i,5,j) = xlong_io(i,j)
              inisrf_0(i,6,j) = msfx_io(i,j)
              inisrf_0(i,7,j) = msfd_io(i,j)
              inisrf_0(i,8,j) = f_io(i,j)
            end do
            do n = 1 , nnsg
              do i = 1 , ix
                inisrf_0(i,8+n,j) = ht1_io(n,i,j)
                inisrf_0(i,8+nnsg+n,j) = satbrt1_io(n,i,j)
                inisrf_0(i,8+nnsg*2+n,j) = snowc_io(n,i,j)
              end do
            end do
          end do
          do j = 1 , jx
            do i = 1 , ix
              ht_io(i,j) = ht_io(i,j)*gti
            end do
          end do
        end if                 ! end if (myid.eq.0)
 
        call mpi_scatter(inisrf_0(1,1,1),ix*(nnsg*3+8)*jxp,mpi_real8,   &
                       & inisrf0(1,1,1), ix*(nnsg*3+8)*jxp,mpi_real8,   &
                       & 0,mpi_comm_world,ierr)
        do j = 1 , jxp
          do i = 1 , ix
            ht(i,j) = inisrf0(i,1,j)
            htsd(i,j) = inisrf0(i,2,j)
            satbrt(i,j) = inisrf0(i,3,j)
            xlat(i,j) = inisrf0(i,4,j)
            xlong(i,j) = inisrf0(i,5,j)
            msfx(i,j) = inisrf0(i,6,j)
            msfd(i,j) = inisrf0(i,7,j)
            f(i,j) = inisrf0(i,8,j)
          end do
          do n = 1 , nnsg
            do i = 1 , ix
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
          do i = 1 , ix
            msfd(i,j) = 1./msfd(i,j)
          end do
        end do
        do j = 1 , jendl
          do i = 1 , ix
            msfx(i,j) = 1./msfx(i,j)
            ht(i,j) = ht(i,j)*gti
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
          print * , '****     Map Projection: ' , proj ,                &
               &'                ****'
          print * , '****     IX=' , ix , ' JX=' , jx , ' KX=' , kx ,   &
               &'             ****'
          print * , '****     PTOP=' , ptsp , ' DX=' , dxsp ,           &
               &'       ****'
          print * , '****     CLAT= ' , clat , ' CLON=' , clon ,        &
               &'    ****'
          if ( proj.eq.'MERCAT' ) print * , '****     PLAT= ' , plat ,  &
                                       &' PLON=' , plon , '    ****'
          print * ,                                                     &
               &'***************************************************'
          print * , ' '
 
        end if
        call mpi_sendrecv(ht(1,jxp),ix,mpi_real8,ieast,1,               &
                        & ht(1,0),ix,mpi_real8,iwest,1,                 &
                        & mpi_comm_world,mpi_status_ignore,ierr)
        call mpi_sendrecv(ht(1,1),ix,mpi_real8,iwest,2,                 &
                        & ht(1,jxp+1),ix,mpi_real8,ieast,2,             &
                        & mpi_comm_world,mpi_status_ignore,ierr)
        call mpi_sendrecv(msfx(1,jxp-1),ix*2,mpi_real8,ieast,           &
                        & 1,msfx(1,-1),ix*2,mpi_real8,iwest,            &
                        & 1,mpi_comm_world,mpi_status_ignore,ierr)
        call mpi_sendrecv(msfx(1,1),ix*2,mpi_real8,iwest,2,             &
                        & msfx(1,jxp+1),ix*2,mpi_real8,ieast,           &
                        & 2,mpi_comm_world,mpi_status_ignore,ierr)
        call mpi_sendrecv(msfd(1,jxp-1),ix*2,mpi_real8,ieast,           &
                        & 1,msfd(1,-1),ix*2,mpi_real8,iwest,            &
                        & 1,mpi_comm_world,mpi_status_ignore,ierr)
        call mpi_sendrecv(msfd(1,1),ix*2,mpi_real8,iwest,2,             &
                        & msfd(1,jxp+1),ix*2,mpi_real8,ieast,           &
                        & 2,mpi_comm_world,mpi_status_ignore,ierr)
      end if
#else
      if ( .not.ifrest ) then
        print * , 'HT'
        read (iutin,rec=2,iostat=ierr1) ((sp2d(i,j),j=1,jx),i=1,ix)
        do j = 1 , jx
          do i = 1 , ix
            ht(i,j) = dble(sp2d(i,j))
          end do
        end do
        if ( nsg.gt.1 ) then
          read (iutin1,rec=2,iostat=ierr1)                              &
              & ((sp2d1(i,j),j=1,jx*nsg),i=1,ix*nsg)
          do j = 1 , jx*nsg
            do i = 1 , ix*nsg
              jj = mod(j,nsg)
              if ( jj.eq.0 ) jj = nsg
              ii = mod(i,nsg)
              if ( ii.eq.0 ) ii = nsg
              k = (jj-1)*nsg + ii
              jj = (j+nsg-1)/nsg
              ii = (i+nsg-1)/nsg
              ht1(k,ii,jj) = sp2d1(i,j)*gti
            end do
          end do
        else
          do j = 1 , jx
            do i = 1 , ix
              ht1(1,i,j) = sp2d(i,j)*gti
            end do
          end do
        end if

        print * , 'HTSD'
        read (iutin,rec=3,iostat=ierr1) ((sp2d(i,j),j=1,jx),i=1,ix)
        do j = 1 , jx
          do i = 1 , ix
            htsd(i,j) = dble(sp2d(i,j))
          end do
        end do
        print * , 'SATBRT'
        read (iutin,rec=4,iostat=ierr1) ((sp2d(i,j),j=1,jx),i=1,ix)
        do j = 1 , jx
          do i = 1 , ix
            satbrt(i,j) = dble(sp2d(i,j))
          end do
        end do
        if ( nsg.gt.1 ) then
          read (iutin1,rec=4,iostat=ierr1)                              &
              & ((sp2d1(i,j),j=1,jx*nsg),i=1,ix*nsg)
          do j = 1 , jx*nsg
            do i = 1 , ix*nsg
              jj = mod(j,nsg)
              if ( jj.eq.0 ) jj = nsg
              ii = mod(i,nsg)
              if ( ii.eq.0 ) ii = nsg
              k = (jj-1)*nsg + ii
              jj = (j+nsg-1)/nsg
              ii = (i+nsg-1)/nsg
              satbrt1(k,ii,jj) = sp2d1(i,j)
            end do
          end do
        else
          do j = 1 , jx
            do i = 1 , ix
              satbrt1(1,i,j) = satbrt(i,j)
            end do
          end do
        end if
        print * , 'XLAT'
        read (iutin,rec=5,iostat=ierr1) ((sp2d(i,j),j=1,jx),i=1,ix)
        do j = 1 , jx
          do i = 1 , ix
            xlat(i,j) = dble(sp2d(i,j))
          end do
        end do
        print * , 'XLONG'
        read (iutin,rec=6,iostat=ierr1) ((sp2d(i,j),j=1,jx),i=1,ix)
        do j = 1 , jx
          do i = 1 , ix
            xlong(i,j) = dble(sp2d(i,j))
          end do
        end do
        print * , 'MSFX'
        read (iutin,rec=9,iostat=ierr1) ((sp2d(i,j),j=1,jx),i=1,ix)
        do j = 1 , jx
          do i = 1 , ix
            msfx(i,j) = dble(sp2d(i,j))
          end do
        end do
        print * , 'MSFD'
        read (iutin,rec=10,iostat=ierr1) ((sp2d(i,j),j=1,jx),i=1,ix)
        do j = 1 , jx
          do i = 1 , ix
            msfd(i,j) = dble(sp2d(i,j))
          end do
        end do
        print * , 'F'
        read (iutin,rec=11,iostat=ierr1) ((sp2d(i,j),j=1,jx),i=1,ix)
        do j = 1 , jx
          do i = 1 , ix
            f(i,j) = dble(sp2d(i,j))
          end do
        end do
        print * , 'SNOWC'
        read (iutin,rec=12,iostat=ierr1) ((sp2d(i,j),j=1,jx),i=1,ix)
        do j = 1 , jx
          do i = 1 , ix
            do n = 1 , nnsg
              snowc(n,i,j) = dble(sp2d(i,j))
            end do
          end do
        end do
        if ( ierr1.ne.0 ) then
          write (aline,*) '  Check ibyte in parameter.inc: ibyte = ' ,  &
                        & ibyte
          call fatal(__FILE__,__LINE__,'REACHED EOF')
        end if
!------invert mapscale factors:
        do j = 1 , jx
          do i = 1 , ix
            msfd(i,j) = 1./msfd(i,j)
          end do
        end do
        do j = 1 , jx
          do i = 1 , ix
            msfx(i,j) = 1./msfx(i,j)
            ht(i,j) = ht(i,j)*gti
          end do
        end do
        print * , ' '
        print * , '***************************************************'
        print * , '***************************************************'
        print * , '**** RegCM IS BEING RUN ON THE FOLLOWING GRID: ****'
        print * , '****     Map Projection: ' , proj ,                  &
             &'                ****'
        print * , '****     IX=' , ix , ' JX=' , jx , ' KX=' , kx ,     &
             &'             ****'
        print * , '****     PTOP=' , ptsp , ' DX=' , dxsp ,             &
             &'       ****'
        print * , '****     CLAT= ' , clat , ' CLON=' , clon ,          &
             &'    ****'
        if ( proj.eq.'MERCAT' ) print * , '****     PLAT= ' , plat ,    &
                                     &' PLON=' , plon , '    ****'
        print * , '***************************************************'
        print * , ' '

      end if
#endif

!
!-----compute dsigma and half sigma levels.
!
      do k = 1 , kx
        dsigma(k) = sigma(k+1) - sigma(k)
        a(k) = 0.5*(sigma(k+1)+sigma(k))
      end do
 
      do k = 1 , kx
        if ( a(k).lt.0.4 ) then
          anudg(k) = 3.
        else if ( a(k).lt.0.8 ) then
          anudg(k) = 2.
        else
          anudg(k) = 1.
        end if
      end do
!
!-----specify heating profile (twght) and weighted function
!     for moisture fluxes due to convection (vqflx)
!     assume base of cloud varies as  < kbase = 5,kx >
!     top  of cloud varies as  < ktop  = 1,kbase-3 >
!     exceptions to this are treated explicitly in subroutine
!     "cupara".
!
      do kbase = 5 , kx
        do ktop = 1 , kbase - 3
          do k = 1 , kx
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
          xtop = dlog((100.-ptop)*sigma(ktop)+ptop)
          xbot = dlog((100.-ptop)*sigma(kbase+1)+ptop)
          bb = xtop + xbot
          cc = xtop*xbot
          vqmax = 0.
          ssum = 0.
          xx = xtop
          yy = xbot
          wk = (xx*xx) - (bb*xx) + cc
          qk = -((yy*yy)-(bb*yy)+cc)
          do k = ktop , kbase
            xx = dlog((100.-ptop)*sigma(k+1)+ptop)
            yy = dlog((100.-ptop)                                       &
               & *(sigma(ktop)+sigma(kbase+1)-sigma(k+1))+ptop)
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
      sigtbl = (70.-ptop)/(101.3-ptop)
      kt = 1
      do k = kx , 1 , -1
        delsig = a(k) - sigtbl
        if ( delsig.le.0. ) then
          kt = k
          exit
        end if
      end do
      write (aline, *) ' Index of highest allowed pbl:  kt = ' , kt
      call say
      write (aline, *) ' '
      call say
!
      if ( ipptls.eq.1 ) then
        write (aline, *) 'AUTO-CONVERSION RATE:  LAND=' , qck1land ,    &
               &'                      OCEAN=' , qck1oce
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
        kbmax = kx
        do k = 1 , kx - 1
          if ( a(k).le.skbmax ) kbmax = kx - k
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
          do i = 1 , ixm1
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
#else
        do j = 1 , jxm1
          do i = 1 , ixm1
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
#endif
      else if ( icup.eq.3 ) then
        write (aline,*) ' The Betts-Miller Convection scheme is not' ,  &
                       &' properly implemented'
        call say
        call fatal(__FILE__,__LINE__,'BETTS-MILLER NOT WORKING')
      else if ( icup.eq.4 ) then
        cllwcv = 0.5E-4    ! Cloud liquid water content for convective precip.
        clfrcvmax = 0.25   ! Max cloud fractional cover for convective precip.
        minorig = kx
        do k = 1 , kx
          if ( a(k).le.minsig ) minorig = kx - k
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
      do k = 2 , kx
        twt(k,1) = (sigma(k)-a(k-1))/(a(k)-a(k-1))
        twt(k,2) = 1. - twt(k,1)
        qcon(k) = (sigma(k)-a(k))/(a(k-1)-a(k))
      end do
 
      chibot = 450.
      ptmb = 10.*ptop
      pz = a(1)*(1000.-ptmb) + ptmb
      if ( pz.gt.chibot ) call fatal(__FILE__,__LINE__,                 &
                                    &'VERTICAL INTERPOLATION ERROR')
      do k = 1 , kx
        pk = a(k)*(1000.-ptmb) + ptmb
        if ( pk.le.chibot ) kchi = k
      end do
 
!
!-----compute the k level under which the maximum equivalent potential
!     temperature will be regarded as the origin of air parcel that
!     produces cloud (used in the cumulus parameterization scheme).
!
      sig700 = (70.-ptop)/(100.-ptop)
      do k = 1 , kx
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
      if ( icup.eq.3 ) call lutbl(ptop)
!
!-----print out the parameters specified in the model.
!
!     print outparam
!     print physicsparam
!     print timeparam
!.....for large domain:
!
      daymax = ntimax/1440.
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
        do k = 1 , kx
          print 99011 , k , sigma(k) , a(k) , dsigma(k) , twt(k,1) ,    &
              & twt(k,2) , qcon(k)
        end do
        print 99012 , kxp1 , sigma(kxp1)
        print 99013 , daymax
        print 99014 , dt
        print 99015 , dx
        print 99016 , jx , ix
        print 99017 , kx
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
      do k = 1 , kx
        print 99011 , k , sigma(k) , a(k) , dsigma(k) , twt(k,1) ,      &
            & twt(k,2) , qcon(k)
      end do
      print 99012 , kxp1 , sigma(kxp1)
      print 99013 , daymax
      print 99014 , dt
      print 99015 , dx
      print 99016 , jx , ix
      print 99017 , kx
      print 99018 , xkhz
      print 99019 , xkhmax
#endif
99001 format ('fort.',i2)
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
99013 format (//' maximum time = ',f8.3,' days.')
99014 format (' time step = ',f7.2,' seconds')
99015 format (' dx = ',f7.0,' meters')
99016 format (' grid points (x,y) = (',i3,',',i3,')')
99017 format (' number of levels = ',i2)
99018 format (' constant hor. diff. coef. = ',e12.5,' m*m/s')
99019 format (' maximum  hor. diff. coef. = ',e12.5,' m*m/s')
!
      end subroutine param

C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C    This file is part of RegCM model.
C
C    RegCM model is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    RegCM model is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
C
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine param
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine defines various model parameters.               c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#     include "include/param1.cb"
#     include "include/param2.cb"
#     include "include/param3.cb"
#     include "include/iunits.cb"
#     include "include/pmoist.cb"
#     include "include/bats.cb2"
#     include "include/main.cb"
#     include "include/bdycod.cb"
#     include "include/trachem.cb"
#     include "include/cupeman.cb"
#     include "include/date.cb"
#     include "include/tiedparam0.cb"
#     include "include/tied_constants.cb"
#     include "include/culookup.cb"
Cchem2
Cchem2_
#ifdef MPP1
#     include "include/io.cb"
      include 'mpiregcm.h'
      integer status(MPI_STATUS_SIZE)
      integer ierr
      real(kind=8)  inisrf0(ix,NNSG*3+8,jxp),inisrf_0(ix,NNSG*3+8,mjx)
      COMMON /TMPSRF/ inisrf0,inisrf_0
#endif
c
c----------------------------------------------------------------------
c-----vqrang is the range limit on vqflx.
c.....qdcrit is the precipitation threshold for moisture convergence.
c
      real(kind=8)  dtsplit(nsplit)
c      character*1 xchar
      real(kind=8)  vqrang
      data vqrang /5.0e-4/
      integer mmd(12),i,j,n,ii,jj
      data (mmd(i),i=1,12) /31,28,31,30,31,30,31,31,30,31,30,31/
      integer igrads,ibigend

      real(kind=8) afracl ! frac. cover for conv. precip. when dx=dxlarg
      real(kind=8) afracs !   "     "    "    "      "     "   dx=dxsmal
      real(kind=8) dlargc
      real(kind=8) dsmalc
      real(kind=8) dxtemc
c
c----------------------------------------------------------------------
c-----namelist:
c
      namelist/restartparam/ifrest,idate0,idate1,idate2,nslice

      namelist/timeparam/radfrq,abatm,abemh,dt,ibdyfrq

      namelist/outparam/ ifsave,savfrq,iftape,tapfrq,ifprt,prtfrq,kxout
     &        ,jxsex,ifrad,radisp,ifbat,ifsub,batfrq,iotyp,ibintyp
Cchem2
     &        ,ifchem, chemfrq
c
      namelist/physicsparam/ibltyp,iboudy,icup,igcc,ipgf,iemiss,lakemod
     &        ,ipptls,iocnflx,
Cchem2
     &         ichem    
Cchem2_
      namelist/subexparam/ncld,fcmax,qck1land,qck1oce,gulland,guloce
     &        ,rhmax,rh0oce,rh0land,cevap,caccr,tc0,cllwcv,clfrcvmax

      namelist/grellparam/ shrmin,shrmax,edtmin,edtmax,edtmino
     &        ,edtmaxo,edtminx,edtmaxx,pbcmax,mincld,htmin,htmax
     &        ,skbmax,dtauc

      namelist/emanparam/ minsig, elcrit, tlcrit, entp, sigd, sigs
     &                  , omtrain, omtsnow, coeffr, coeffs, cu
     &                  , betae, dtmax, alphae, damp

Cchem2
      namelist/chemparam/ ichremlsc,ichremcvc,ichdrdepo,ichcumtra,
     &         idirect,
     &         mixtype,chtrname,chtrsol,chtrdpv,dustbsiz
Cchem2_
c
      integer ns,mdate1,myear,mmon,mday,my1,my2,my3
      integer m,k,kbase,ktop
      real(kind=8)  eqoft,bb,cc,ssum,xx,xtop,xbot,vqmax
      real(kind=8)  yy,wk,qk,wkp1,qkp1,sigtbl,delsig,pi,chibot
      real(kind=8)  ptmb,pz,pk,sig700,daymax

      real(kind=4)  sp1d(kxp1)             ! real*4 for compact input
#ifdef MPP1
      real(kind=4)  sp2d(ix,mjx),sp2d1(ix*NSG,mjx*NSG)
#else
      real(kind=4)  sp2d(ix,jx),sp2d1(ix*NSG,jx*NSG)
#endif
      COMMON /TRAPUPA/ sp2d,sp2d1
      real(kind=4)  clat,clon,dxsp,ptsp,plat,plon
      character proj*6
      COMMON /GRADS/clat,clon,dxsp,ptsp,plat,plon,proj
      real(kind=4)  TRUELATL,TRUELATH
      COMMON /SAVEPAR/ TRUELATL,TRUELATH
      integer nyy,nxx,kzz
      
      integer ierr1
      real(kind=4)  grdfac
      character*7 finm
c----------------------------------------------------------------------
c-----specify unit numbers for input/output.
c     units 10,101,14  are input.
c     iutin : input initial data for large domain.
c     iutbc : input slab temperature for large domain ,
c             boundary values and tendencies for large doamin.
c     iutrs : input saved file for large domain for restart from
c             previous forecast.
c     iutdat: output for dataflow, if iftape=.true.
c     iutsav: output saved file for restart, if ifsave=.true. or
c
      qdcrit =3.0e-7
      iutin = 10
      iutin1= 11
      iutbc = 101
C_sgi iutbc = 71
      iutrs = 14
c** initialize output file unit number
      iutdat = 50
c** initialize save file unit number
      iutsav = 52
c** initialize bats file unit number
      iutbat = 54
      iutsub = 55
c** initialize radiation file unit number
      iutrad = 56
c** initalize lake file unit number
      iutlak = 58
c** initialize chemistry file number    
Cchem2 
c      emission file       
       iutopt = 12
       iutchsrc = 13
c      out file
       iutchem = 60 
Cchem2_
c
c----------------------------------------------------------------------
c-----specify the parameters used in the model:
c
c     rfstrt : whether the model is to be started with refined time
c              steps in the first hour. if ifrest=.true., this parameter
c              is never used.
c              = .true. ; yes
c              = .false. ; no
c
c     ifrest : whether this run is restarting from a saved tape
c              (unit 14).
c              = .true. ; yes
c              = .false. ; no
c
c     ktaur  : if ifrest=.false., ktaur=0.
c              if ifrest=.true., set ktaur=ktau, which is the time-step
c                                the model output is saved.
c
c     ifsave : specify whether a saved tape (unit 24) will be written
c              for restarting.
c              = .true. ; yes
c              = .false. ; no
c
c     savfrq : if ifsave=.true., specify the interval in hours
c              between save operations.
c
c     ibltyp : specify whether bats pbl or holtslag's pbl
c              parameterization is to be used in the model.
c              = 0 ; frictionless
c              = 1 ; holtslag pbl (holtslag, 1990)
c
c     radfrq : specify the frequency in
c              minutes, the solar radiation will be computed in
c              subroutine "sfcrd". a typical value would be 30 minutes.
c
c     lakemod: inclusion of the hostetler lake model
c              = 0 ; no
c              = 1 ; yes
c
c     kxout  : the k level of the horizontal slice for printer output.
c
c     jxsex  : the j index of the north-south vertical slice for printer
c              output.
c
c     icup   : type of cumulus parameterization
c              = 1 ; kuo
c              = 2 ; grell
c              = 3 ; betts-miller (1986)
c              = 4 ; emanuel (1991)
c
c       igcc   : Grell Scheme Convective Closure Assumption
c                = 1 ; Arakawa & Schubert (1974)
c                = 2 ; Fritsch & Chappell (1980)
c
c     ipptls : type of moisture scheme
c              = 1 ; explicit moisture (SUBEX; Pal et al 2000)
c
c     iocnflx: type of ocean flux parameterization
c              = 1 ; BATS
c              = 2 ; Zeng et al.
c
c     iboudy : specify the laterial boundary conditions.
c              = 0 ; fixed.
c              = 1 ; relaxation, linear technique.
c              = 2 ; time-dependent (from observations or large-scale
c                                    model).
c              = 3 ; time and inflow/outflow dependent.
c              = 4 ; sponge (perkey & kreitzberg, mwr 1976).
c              = 5 ; relaxation, exponential technique.
c
c     iftape : whether you want output in a saved tape (unit 20) for
c              analyses in dataflow.
c              = .true. ; yes
c              = .false. ; no
c
c         tapfrq : if iftape=true., specify the output interval in hours.
c
c     ifrad  : whether you want radiation output saved
c              = .true. ; yes
c              = .false. ; no
c
c         radisp : if ifrad=1, specify the output interval in hours.
c
c     ifprt  : whether you want printer output or not.
c              = 0 ; no
c              = 1 ; yes
c
c     prtfrq : if ifprt=1, specify the output interval in mimutes.
c
c     iotyp  : Type of output files,
c                 1=direct access (GrADS); 2=sequential w/ time listing
c
c     ibintyp : Type of binary output for GrADS ctl files
c                 1=big endian; 2=little endian
c
c     maschk : specify the frequency in time steps, the mass-
c              conservation information will be printed out.
c
c----------------------------------------------------------------------
c-----default values for all the options:
c     (can be overwritten by namelist input).
c------namelist restart param:
c
      ifrest = .false.     ! Restart?:  t=true; f=false
      idate0 = 1993063000  ! Start date of the initial simulation
      idate1 = 1993063000  ! Start date of this simulation
      idate2 = 1993080100  ! End Date this simulation
      nslice = 4      ! # of days for next model run (used with restart)
      ! note: beginning/end forecast time set in restart.mm4

c------namelist timeparam:
c
      radfrq = 30.    ! time interval in min solar rad caluclated
      abatm = 600.    ! time interval at which bats is called (secs)
      abemh = 12.     ! time interval absorption-emission calculated (hours)
      dt =  200.      ! time step in seconds
      ibdyfrq = 6     ! boundary condition interval (hours)

c-----namelist outparam:       note: * signifies currently not in namelist
c
      rfstrt = .false.      ! *
      ifsave = .false.
         savfrq = 6.
      iftape = .true.
         tapfrq = 6.
      ifrad = .true.
         radisp = 6.0    ! time interval for disposing rad output (hrs)
      ifbat = .true.
        batfrq = 1.0    ! time interval for disposing bats output (hrs)
      ifsub = .true.
      ifprt  = .true.
         prtfrq =  12.
         kxout   = kx
         jxsex   = 25
      iotyp = 1
        ibintyp = 1
      maschk = 10           ! * defined below

Cchem2
      ifchem=.false.
      chemfrq = 6.0   ! time interval for disposeing chem output (hrs)
Cchem2_

c
c----------------------------------------------------------------------
c-----namelist physicsparam:
c
      ibltyp = 1
      iboudy = 1
      icup   = 1
      ipptls = 1
      igcc   = 1
      ipgf   = 1
      iemiss = 1
      iocnflx = 1
      lakemod = 1
      ichem = 0 
c
c----------------------------------------------------------------------
c------namelist subexparam:
      ncld     = 1       ! # of bottom model levels with no clouds (rad only)
      fcmax    = 0.80    ! Maximum cloud fraction cover (rad only)
!     qck1land = 0.0005  ! Autoconversion Rate for Land
!     qck1oce  = 0.0005  ! Autoconversion Rate for Ocean
      qck1land = 0.00025 ! Autoconversion Rate for Land
      qck1oce  = 0.00025 ! Autoconversion Rate for Ocean
      gulland  = 0.4     ! Fract of Gultepe eqn (qcth) when precip occurs (land)
      guloce   = 0.4     ! Fract of Gultepe eqn (qcth) for ocean
      rhmax    = 1.01    ! RH at whicn FCC = 1.0
      rh0oce   = 0.90    ! Relative humidity threshold for ocean
      rh0land  = 0.80    ! Relative humidity threshold for land
      tc0      = 238.0   ! Below this temperature, rh0 begins to approach unity
!     cevap    = 0.2e-4  ! Raindrop evap rate coef [[(kg m-2 s-1)-1/2]/s]
      cevap    = 1.0e-3  ! Raindrop evap rate coef [[(kg m-2 s-1)-1/2]/s]
!     caccr    = 6.0     ! Raindrop accretion rate [m3/kg/s]
      caccr    = 3.0     ! Raindrop accretion rate [m3/kg/s]
      cllwcv   = 0.3e-3  ! Cloud liquid water content for convective precip.
      clfrcvmax= 0.25    ! Max cloud fractional cover for convective precip.

c------namelist grellparam:
      shrmin  = 0.25    ! Minimum Shear effect on precip eff.
      shrmax  = 0.50    ! Maximum Shear effect on precip eff.
      edtmin  = 0.25    ! Minimum Precipitation Efficiency
      edtmax  = 1.00    ! Maximum Precipitation Efficiency
      edtmino = 0.0     ! Minimum Precipitation Efficiency (o var)
      edtmaxo = 1.00    ! Maximum Precipitation Efficiency (o var)
      edtminx = 0.25    ! Minimum Precipitation Efficiency (x var)
      edtmaxx = 1.00    ! Maximum Precipitation Efficiency (x var)
      pbcmax  = 150.    ! Max depth (mb) of stable layer b/twn LCL & LFC
      mincld  = 150.    ! Min cloud depth (mb).
      htmin   = -250.   ! Min convective heating
      htmax   = 500.    ! Max convective heating
      skbmax  = 0.4     ! Max cloud base height in sigma
      dtauc   = 30.     ! Fritsch & Chappell (1980) ABE Removal Timescale (min)

c------namelist emanparam:
      minsig  = 0.95    ! Lowest sigma level from which convection can originate
      elcrit  = 0.0011  ! AUTOCONVERSION THRESHOLD WATER CONTENT (gm/gm)
      tlcrit  = -55.0   ! BELOW TLCRIT AUTO-CONVERSION THRESHOLD IS ZERO
      entp    = 1.5     ! COEFFICIENT OF MIXING IN THE ENTRAINMENT FORMULATION
      sigd    = 0.05    ! FRACTIONAL AREA COVERED BY UNSATURATED DNDRAFT
      sigs    = 0.12    ! FRACTION OF PRECIPITATION FALLING OUTSIDE OF CLOUD
      omtrain = 50.0    ! FALL SPEED OF RAIN (P/s)
      omtsnow = 5.5     ! FALL SPEED OF SNOW (P/s)
      coeffr  = 1.0     ! COEFFICIENT GOVERNING THE RATE OF RAIN EVAPORATION
      coeffs  = 0.8     ! COEFFICIENT GOVERNING THE RATE OF SNOW EVAPORATION
      cu      = 0.7     ! COEFFICIENT GOVERNING CONVECTIVE MOMENTUM TRANSPORT
      betae   = 10.0    ! CONTROLS DOWNDRAFT VELOCITY SCALE
      dtmax   = 0.9     ! MAX NEGATIVE PARCEL TEMPERATURE PERTURBATION BELOW LFC
      alphae  = 0.2     ! CONTROLS THE APPROACH RATE TO QUASI-EQUILIBRIUM
      damp    = 0.1     ! CONTROLS THE APPROACH RATE TO QUASI-EQUILIBRIUM
c
cc------namelist chemparam ; ( 0= none, 1= activated)
      ichremlsc = 1     ! tracer removal by large scale clouds 
      ichremcvc = 1     ! tracer removal by convective clouds
      ichdrdepo = 1     ! tracer dry deposition
      ichcumtra = 1     ! tracer convective transport
      idirect   = 1     ! tracer direct effect 
      mixtype   = 1
c---------------------------------------------------------------------
c-----read in namelist variables:
c
#ifdef MPP1
      if(myid .eq. 0) then
        open (2,file='regcm.in',status='old')
        read(2,restartparam)
        print*,'RESTARTPARAM READ IN'
        read(2,timeparam)
        print*,'TIMEPARAM READ IN'
        read(2,outparam)
        print*,'OUTPARAM READ IN'
        read(2,physicsparam)
        print*,'PHYSICSPARAM READ IN'
        if(ipptls.eq.1) then
          read(2,subexparam)
          print*,'SUBEXPARAM READ IN'
        endif
        if(icup.eq.2) then
          read(2,grellparam)
          print*,'GRELLPARAM READ IN'
        else if(icup.eq.4) then
          read(2,emanparam)
          print*,'EMANPARAM READ IN'
        endif
        if(ichem.eq.1) then
          read(2,chemparam)
          print*,'CHEMPARAM READ IN'
        endif
        close(2)
      endif
#else
      read (*,restartparam)
      print*,'RESTARTPARAM READ IN'
      read (*,timeparam)
      print*,'TIMEPARAM READ IN'
      read (*,outparam)
      print*,'OUTPARAM READ IN'
      read (*,physicsparam)
      print*,'PHYSICSPARAM READ IN'
      if(ipptls.eq.1) then
        read (*,subexparam)
        print*,'SUBEXPARAM READ IN'
      endif
      if(icup.eq.2) then
        read (*,grellparam)
        print*,'GRELLPARAM READ IN'
      else if(icup.eq.4) then
        read (*,emanparam)
        print*,'EMANPARAM READ IN'
      endif
      if(ichem.eq.1) then
        read (*,chemparam)
        print*,'CHEMPARAM READ IN'
      endif
#endif
#ifdef MPP1
      call MPI_BCAST(ifrest,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(idate0,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(idate1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(idate2,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(nslice,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(radfrq,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(abemh,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(abatm,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dt,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ibdyfrq,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(ifsave,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(savfrq,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(iftape,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(tapfrq,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ifrad,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(radisp,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ifbat,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ifsub,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(batfrq,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ifprt,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(prtfrq,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(kxout,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(jxsex,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(iotyp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ibintyp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ifchem,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(chemfrq,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)

      call MPI_BCAST(iboudy,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ibltyp,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(icup,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(igcc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ipptls,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(iocnflx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ipgf,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(iemiss,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(lakemod,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ichem,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)

      if(ipptls.eq.1) then
      call MPI_BCAST(ncld,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(fcmax,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(qck1land,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(qck1oce,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(gulland,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(guloce,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(rhmax,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(rh0oce,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(rh0land,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(tc0,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(cevap,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(caccr,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      endif

      if(icup.eq.2) then
      call MPI_BCAST(shrmin,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(shrmax,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(edtmin,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(edtmax,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(edtmino,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(edtmaxo,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(edtminx,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(edtmaxx,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(pbcmax,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(mincld,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(htmin,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(htmax,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(skbmax,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dtauc,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)

      elseif(icup.eq.4) then
      call MPI_BCAST(minsig,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(elcrit,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(tlcrit,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(entp,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(sigd,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(sigs,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(omtrain,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(omtsnow,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(coeffr,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(coeffs,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(cu,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(betae,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dtmax,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(alphae,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(damp,1,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      endif

      if(ichem.eq.1) then
      call MPI_BCAST(ichremlsc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ichremcvc,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ichdrdepo,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ichcumtra,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(idirect,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      do n=1,ntr
      call MPI_BCAST(chtrname(n),5,MPI_CHARACTER,0,MPI_COMM_WORLD,ierr)
      enddo
      call MPI_BCAST(chtrsol,ntr,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(chtrdpv,ntr*2,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dustbsiz,nbin*2,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      endif
      
#endif
      if(ichem.eq.0) ifchem=.false.
cas
      if (mod(anint(radfrq*60.),anint(dt)).ne.0) then
#ifdef MPP1
         if(myid .eq. 0) then
           print*,'INCONSISTENT TIMESTEPS SPECIFIED'
           print*,'RADFRQ=',radfrq,'DT=',dt
           print*,'stop Radiation Timestep (param.f)'
         endif
         call MPI_ABORT(MPI_COMM_WORLD,ierr)
#else
         print*,'INCONSISTENT TIMESTEPS SPECIFIED'
         print*,'RADFRQ=',radfrq,'DT=',dt
         stop 'Radiation Timestep (param.f)'
#endif
      end if
      if (mod(anint(abatm),anint(dt)).ne.0) then
#ifdef MPP1
         if(myid .eq. 0) then
           print*,'INCONSISTENT TIMESTEPS SPECIFIED'
           print*,'ABATM=',abatm,'DT=',dt
           print*,'stop Surface Timestep (param.f)'
         endif
         call MPI_ABORT(MPI_COMM_WORLD,ierr)
#else
         print*,'INCONSISTENT TIMESTEPS SPECIFIED'
         print*,'ABATM=',abatm,'DT=',dt
         stop 'Surface Timestep (param.f)'
#endif
      end if
      if (mod(anint(batfrq*3600.),anint(abatm)).ne.0) then
#ifdef MPP1
         if(myid .eq. 0) then
           print*,'INCONSISTENT TIMESTEPS SPECIFIED'
           print*,'BATFRQ=',batfrq,'ABATM=',abatm
           print*,'stop Surface/Radiation Timestep (param.f)'
         endif
         call MPI_ABORT(MPI_COMM_WORLD,ierr)
#else
         print*,'INCONSISTENT TIMESTEPS SPECIFIED'
         print*,'BATFRQ=',batfrq,'ABATM=',abatm
         stop 'Surface/Radiation Timestep (param.f)'
#endif
      end if
      if (mod(anint(abemh*3600.),anint(dt)).ne.0) then
#ifdef MPP1
         if(myid .eq. 0) then
           print*,'INCONSISTENT TIMESTEPS SPECIFIED'
           print*,'ABEMH=',abemh,'DT=',dt
           print*,'stop Abs-Ems Timestep (param.f)'
         endif
         call MPI_ABORT(MPI_COMM_WORLD,ierr)
#else
         print*,'INCONSISTENT TIMESTEPS SPECIFIED'
         print*,'ABEMH=',abemh,'DT=',dt
         stop 'Abs-Ems Timestep (param.f)'
#endif
      end if
      if (mod(anint(abemh*60.),anint(radfrq)).ne.0) then
#ifdef MPP1
         if(myid .eq. 0) then
           print*,'INCONSISTENT TIMESTEPS SPECIFIED'
           print*,'ABEMH=',abemh,'RADFRQ=',radfrq
           print*,'stop Longwave/shortwave Radiation Timestep (param.f)'
         endif
         call MPI_ABORT(MPI_COMM_WORLD,ierr)
#else
         print*,'INCONSISTENT TIMESTEPS SPECIFIED'
         print*,'ABEMH=',abemh,'RADFRQ=',radfrq
         stop 'Longwave/shortwave Radiation Timestep (param.f)'
#endif
      end if

c-----reset the options/calculate variables using namelist info:
c
      ndate0  = idate1
      nsavfrq = nint(3600.*savfrq)
      ntapfrq = nint(3600.*tapfrq)
      nprtfrq = nint(3600.*prtfrq)
      ktau = 0
      ktaur= 0
      xtime = 0.
      ntime = 0
      dtsplit(2) = dt/2.
      dtsplit(1) = dt/4.
      do ns=1,nsplit
         dtau(ns)=dtsplit(ns)
      enddo
#ifdef MPP1
      if (myid.eq.0) print *,' dtau = ', dtau
#else
      print *,' dtau = ', dtau
#endif
      dt0 = dt      !store original dt
      maschk=nint(prtfrq*3600./dt) !convert prtfrq to #of time steps
      nradisp=nint(radisp*3600) !convert radisp to time steps
      ifrabe=nint(3600.*abemh/dt)  !abemh is time interval abs./emis. calc.
      kbats = nint(3600.*batfrq)
      nbatst=nint(abatm/dt)
      dt2 = 2.*dt
Cchem2
      kchem = nint(3600.*chemfrq)  ! convert chemfrq to time steps
Cchem2_
c.....calculate the time step in minutes.
      dtmin = dt/60.
      deltmx = dt
c.....compute the time steps for radiation computation.
      ntrad = nint(radfrq/dtmin)
csb lake model mods
c.....compute the time steps for lake model call.
      dtlake = 3600.
      klake = nint(dtlake/dt)
csb end lake model mods
c
      CALL INITDATE
      CALL FINDDATE(NSTRT0,IDATE0)
      CALL FINDDATE(NSTART,IDATE1)
      CALL FINDDATE(NNNEND,IDATE2)
      NNNCHK=NSTART

c     ktaur = nint((NSTART-NSTRT0)*ibdyfrq*60/dtmin)
      ntimax=      (NNNEND-NSTRT0)*ibdyfrq*60
#ifdef MPP1
      if (myid .eq. 0)
     &print *, 'IDATE1, IDATE2, dtmin, ktaur = '
     &       ,  IDATE1, IDATE2, dtmin, ktaur
#else
      print *, 'IDATE1, IDATE2, dtmin, ktaur = '
     a       ,  IDATE1, IDATE2, dtmin, ktaur
#endif
      ldatez= IDATE1
      NNNNNN= NSTART
      lyear = ldatez/1000000
      lmonth=(ldatez-lyear*1000000)/10000
      lday  =(ldatez-lyear*1000000-lmonth*10000)/100
      lhour = mod(ldatez,100)
      idatex= ldatez
      jyear = lyear
      jyearr= jyear
c
c-----specify the julian date and gmt of the initial data.
c     dectim : is the time in minutes after which the solar declination
c              angle must be recalculated.
c
#ifdef MPP1
      if (myid .eq. 0) then
#endif
      print*,'READING HEADER FILE'
      write(finm,101) iutin
101   format('fort.',I2)
      if(NSG.gt.1) THEN
#ifdef MPP1
        open(iutin1,file='fort.11',form='unformatted',status='old'
     &             ,access='direct',recl=ix*mjx*NNSG*ibyte)
#else
        open(iutin1,file='fort.11',form='unformatted',status='old'
     &             ,access='direct',recl=ix*jx*NNSG*ibyte)
#endif
      endif
#ifdef MPP1
      open(iutin,file=finm,form='unformatted',status='old'
     &     ,access='direct',recl=ix*mjx*ibyte)
#else
      open(iutin,file=finm,form='unformatted',status='old'
     &     ,access='direct',recl=ix*jx*ibyte)
#endif
      read(iutin,rec=1,iostat=ierr1)nyy,nxx,kzz,dxsp,clat,clon,plat,plon
     &                             ,grdfac,proj,sp1d,ptsp,igrads,ibigend
     &                             ,TRUELATL,TRUELATH
      print*,'DIMS',nyy,nxx,kzz
      print*,'DOMAIN',dxsp,clat,clon,plat,plon,grdfac
      print*,'PROJ',proj
      print*,'SIGMA',sp1d
      print*,'PTOP',ptsp
      print*,'OUTPUT',igrads,ibigend
      ptop = ptsp
      dx = dxsp
#ifdef MPP1
      if (nyy.ne.ix .or. nxx.ne.mjx .or. kzz.ne.kx) then
        print*,'IMPROPER DIMENSION SPECIFICATION (SUBROUTINE PARAM)'
        print*,'  SET IN regcm.param:  IX=',ix,' JX=',mjx, ' KX=',kx
        print*,'  SET IN TERRAIN: NYY=',nyy,' NXX=',nxx,' NZZ=',kzz
        print*,'  Also check ibyte in regcm.param: ibyte = ',ibyte
        print*,'stop dimensions (param.f)'
        call MPI_ABORT(MPI_COMM_WORLD,ierr)
#else
      if (nyy.ne.ix .or. nxx.ne.jx .or. kzz.ne.kx) then
        print*,'IMPROPER DIMENSION SPECIFICATION (SUBROUTINE PARAM)'
        print*,'  SET IN regcm.param:  IX=',ix,' JX=',jx, ' KX=',kx
        print*,'  SET IN TERRAIN: NYY=',nyy,' NXX=',nxx,' NZZ=',kzz
        print*,'  Also check ibyte in regcm.param: ibyte = ',ibyte
        stop 'dimensions (param.f)'
#endif
      end if
      do k=1,kxp1
        sigma(k) = dble(sp1d(k))
      enddo
#ifdef MPP1
      endif
      call MPI_BCAST(clat,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(plon,1,MPI_REAL,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ptop,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(dx,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(sigma,kxp1,MPI_DOUBLE_PRECISION,0
     &                         ,MPI_COMM_WORLD,ierr)
#endif

Crst-fix
      mdate0 = IDATE0
#ifdef MPP1
      if(myid.eq.0) then
#endif
      print *,'***** mdate = ',mdate0
#ifdef MPP1
      endif
#endif
      myear = mdate0/1000000
      nyear = IDATE1/1000000
      nmonth= (IDATE1-nyear*1000000)/10000
      mdate1 = mdate0 - myear*1000000
      mmon = mdate1/10000
      mday = (mdate1 - mmon*10000)/100
      gmt = mdate1-mmon*10000-mday*100
      jyear0 = myear
Crst-fix_
c
c.....find the year is a leap year or not (assume the initial data
c     is within this century.):
      if(myear.lt.100) myear = 1900+myear
      my1 = mod(myear,4)
      my2 = mod(myear,100)
      my3 = mod(myear,400)
      if (my1 .eq. 0 .and. my2 .ne. 0 .or. my3 .eq. 0) mmd(2) = 29
      julday = mday
      do 85 m=1,mmon-1
       julday = julday+mmd(m)
85    continue
      dectim = (1440.-gmt*60.)

c-----specify the constants used in the model.
c        conf   : condensation threshold.
c        qcth   : threshold for the onset of autoconversion.
c        qck1oce  : constant autoconversion rate for ocean.
c        qck1land : constant autoconversion rate for land.
c        xlv    : latent heat of condensation at 0 c in j/kg.
c        ep1    : constant used in computing virture temperature.
c        ep2    : constant used in computing saturation mixing ratio.
c        svp1, svp2, svp3 : constants used in computing saturation
c                           vapor pressure.
c        avt, bvt : constants used in computing terminal velocity of
c                   raindrops.
c        all the other constants are used to compute the cloud
c        microphysical parameterization (ref. orville & kopp, 1977 jas).
c
      degrad = 0.0174533
      eomeg = 7.292e-5
      dpd = 360./365.25   !0.986301  (degrees/day)
      eqoft = 0.
      solcon = 1.3956e3
      stbolt = 5.76383e-8
      r = 287.
      cp = 1004.
      rovcp = r/cp
      g = 9.8
      rovg = r/g
      karman = 0.4
      rhos = 1.16
      ptop4 = 4.*ptop
      cd = 0.002
      cdsea = 0.0015
      ch = 0.002
      chsea = 0.0015
      alpha = 0.             ! =.2495 in brown-campana; =0. in split explicit
      beta = 1.-2.*alpha
      gnu = 0.10             !  default = 0.1
      omu = 1.-2.*gnu
      gnuhf = 0.5*gnu
      omuhf = 1.-2.*gnuhf
      dx2 = 2.*dx
      dx4 = 4.*dx
      dx8 = 8.*dx
      dx16 = 16.*dx
      dxsq = dx*dx
      c200 = karman*karman*dx/(4.*(100.-ptop))
      c201 = (100.-ptop)/dxsq
      c203 = 1./dxsq
      xkhz = 1.5e-3*dxsq/dt
      xkhmax = dxsq/(64.*dt)
      rv = 461.5
      xlv = 2.5e6
      ep1 = 0.608
      ep2 = 0.622
      xlvocp = xlv/cp
cex5 ---------------------------------
cex5      svp1 = 0.611
cex5      svp2 = 19.84659
cex5      svp3 = 5418.12
      svp1 = 0.6112
      svp2 = 17.67
      svp3 = 29.65
cex5 ---------------------------------

      thrlh1=0.0002
      thrlh2=0.0003
      tauht=7200.
      akht1=dxsq/tauht
      akht2=dxsq/tauht
c
      conf = 1.

#ifdef MPP1
      if (myid .eq. 0) then
#endif
      print *,' input/output parameters '
      print *,' ifsave = ',ifsave,' savfrq = ',savfrq
     &       ,' iftape = ',iftape,' tapfrq = ',tapfrq
     &       ,' ifprt  = ',ifprt,' prtfrq = ',prtfrq
     &       ,' kxout  = ',kxout,' jxsex  = ',jxsex
     &       ,' radisp = ',radisp,
     5        ' batfrq = ',batfrq,' nslice = ',nslice
     &       ,' ifchem = ',ifchem,' chemfrq =',chemfrq
      print *, ' '
      print *, ' physical parameterizations '
      print *,' iboudy = ',iboudy,' icup = ',icup,' igcc =',igcc,
     1        ' ipptls = ',ipptls,' iocnflx = ',iocnflx,
     1        ' ipgf = ',ipgf, 'iemiss = ',iemiss,
     1        ' lakemod = ',lakemod, ' ichem =',ichem
      print *, ' '
      print *,' model parameters '
      print *, ' radfrq = ',radfrq,' abatm = ',abatm,
     &         ' abemh = ',abemh,  ' dt = ',dt

      print *, ' '
      print *, ' ncld = ',ncld

      print *, ' '
#ifdef MPP1
      endif
#endif
c
      if (.not.ifrest) then
#ifdef MPP1
        if (myid .eq. 0) then
#endif
        print*,'HT'
#ifdef MPP1
        read(iutin,rec=2,iostat=ierr1) ((sp2d(i,j),j=1,mjx),i=1,ix)
        do j=1,mjx
        do i=1,ix
          ht_io(i,j)=dble(sp2d(i,j))
        enddo
        enddo
#else
        read(iutin,rec=2,iostat=ierr1) ((sp2d(i,j),j=1,jx),i=1,ix)
        do j=1,jx
        do i=1,ix
          ht(i,j)=dble(sp2d(i,j))
        enddo
        enddo
#endif
        if(NSG.gt.1) THEN
#ifdef MPP1
          read(iutin1,rec=2,iostat=ierr1)((sp2d1(i,j),j=1,mjx*NSG)
     &                                               ,i=1,ix*NSG)
          do j=1,mjx*NSG
          do i=1,ix*NSG
            jj=mod(j,NSG)
            if(jj.eq.0) jj=NSG
            ii=mod(i,NSG)
            if(ii.eq.0) ii=NSG
            k=(jj-1)*NSG+ii
            jj=(j+NSG-1)/NSG
            ii=(i+NSG-1)/NSG
            ht1_io(k,ii,jj)=sp2d1(i,j)*g
          enddo
          enddo
#else
          read(iutin1,rec=2,iostat=ierr1)((sp2d1(i,j),j=1,jx*NSG)
     &                                               ,i=1,ix*NSG)
          do j=1,jx*NSG
          do i=1,ix*NSG
            jj=mod(j,NSG)
            if(jj.eq.0) jj=NSG
            ii=mod(i,NSG)
            if(ii.eq.0) ii=NSG
            k=(jj-1)*NSG+ii
            jj=(j+NSG-1)/NSG
            ii=(i+NSG-1)/NSG
            ht1(k,ii,jj)=sp2d1(i,j)*g
          enddo
          enddo
#endif
        else
#ifdef MPP1
          do j=1,mjx
          do i=1,ix
            ht1_io(1,i,j)=sp2d(i,j)*g
          enddo
          enddo
#else
          do j=1,jx
          do i=1,ix
            ht1(1,i,j)=sp2d(i,j)*g
          enddo
          enddo
#endif
        endif

        print*,'HTSD'
#ifdef MPP1
        read(iutin,rec=3,iostat=ierr1) ((sp2d(i,j),j=1,mjx),i=1,ix)
        do j=1,mjx
        do i=1,ix
          htsd_io(i,j)=dble(sp2d(i,j))
        enddo
        enddo
#else
        read(iutin,rec=3,iostat=ierr1) ((sp2d(i,j),j=1,jx),i=1,ix)
        do j=1,jx
        do i=1,ix
          htsd(i,j)=dble(sp2d(i,j))
        enddo
        enddo
#endif
        print*,'SATBRT'
#ifdef MPP1
        read(iutin,rec=4,iostat=ierr1) ((sp2d(i,j),j=1,mjx),i=1,ix)
        do j=1,mjx
        do i=1,ix
          satbrt_io(i,j)=dble(sp2d(i,j))
        enddo
        enddo
#else
        read(iutin,rec=4,iostat=ierr1) ((sp2d(i,j),j=1,jx),i=1,ix)
        do j=1,jx
        do i=1,ix
          satbrt(i,j)=dble(sp2d(i,j))
        enddo
        enddo
#endif
        IF(NSG.gt.1) then
#ifdef MPP1
           read(iutin1,rec=4,iostat=ierr1)((sp2d1(i,j),j=1,mjx*NSG)
     &                                                ,i=1,ix*NSG)
           do j=1,mjx*NSG
           do i=1,ix*NSG
              jj=mod(j,NSG)
              if(jj.eq.0) jj=NSG
              ii=mod(i,NSG)
              if(ii.eq.0) ii=NSG
              k=(jj-1)*NSG+ii
              jj=(j+NSG-1)/NSG
              ii=(i+NSG-1)/NSG
              satbrt1_io(k,ii,jj)=sp2d1(i,j)
           enddo
           enddo
#else
           read(iutin1,rec=4,iostat=ierr1)((sp2d1(i,j),j=1,jx*NSG)
     &                                                ,i=1,ix*NSG)
           do j=1,jx*NSG
           do i=1,ix*NSG
              jj=mod(j,NSG)
              if(jj.eq.0) jj=NSG
              ii=mod(i,NSG)
              if(ii.eq.0) ii=NSG
              k=(jj-1)*NSG+ii
              jj=(j+NSG-1)/NSG
              ii=(i+NSG-1)/NSG
              satbrt1(k,ii,jj)=sp2d1(i,j)
           enddo
           enddo
#endif
        else
#ifdef MPP1
           do j=1,mjx
           do i=1,ix
              satbrt1_io(1,i,j)=satbrt_io(i,j)
           enddo
           enddo
#else
           do j=1,jx
           do i=1,ix
              satbrt1(1,i,j)=satbrt(i,j)
           enddo
           enddo
#endif
        endif
        print*,'XLAT'
#ifdef MPP1
        read(iutin,rec=5,iostat=ierr1) ((sp2d(i,j),j=1,mjx),i=1,ix)
        do j=1,mjx
        do i=1,ix
          xlat_io(i,j)=dble(sp2d(i,j))
        enddo
        enddo
#else
        read(iutin,rec=5,iostat=ierr1) ((sp2d(i,j),j=1,jx),i=1,ix)
        do j=1,jx
        do i=1,ix
          xlat(i,j)=dble(sp2d(i,j))
        enddo
        enddo
#endif
        print*,'XLONG'
#ifdef MPP1
        read(iutin,rec=6,iostat=ierr1) ((sp2d(i,j),j=1,mjx),i=1,ix)
        do j=1,mjx
        do i=1,ix
          xlong_io(i,j)=dble(sp2d(i,j))
        enddo
        enddo
#else
        read(iutin,rec=6,iostat=ierr1) ((sp2d(i,j),j=1,jx),i=1,ix)
        do j=1,jx
        do i=1,ix
          xlong(i,j)=dble(sp2d(i,j))
        enddo
        enddo
#endif
        print*,'MSFX'
#ifdef MPP1
        read(iutin,rec=9,iostat=ierr1) ((sp2d(i,j),j=1,mjx),i=1,ix)
        do j=1,mjx
        do i=1,ix
          msfx_io(i,j)=dble(sp2d(i,j))
        enddo
        enddo
#else
        read(iutin,rec=9,iostat=ierr1) ((sp2d(i,j),j=1,jx),i=1,ix)
        do j=1,jx
        do i=1,ix
          msfx(i,j)=dble(sp2d(i,j))
        enddo
        enddo
#endif
        print*,'MSFD'
#ifdef MPP1
        read(iutin,rec=10,iostat=ierr1) ((sp2d(i,j),j=1,mjx),i=1,ix)
        do j=1,mjx
        do i=1,ix
          msfd_io(i,j)=dble(sp2d(i,j))
        enddo
        enddo
#else
        read(iutin,rec=10,iostat=ierr1) ((sp2d(i,j),j=1,jx),i=1,ix)
        do j=1,jx
        do i=1,ix
          msfd(i,j)=dble(sp2d(i,j))
        enddo
        enddo
#endif
        print*,'F'
#ifdef MPP1
        read(iutin,rec=11,iostat=ierr1) ((sp2d(i,j),j=1,mjx),i=1,ix)
        do j=1,mjx
        do i=1,ix
          f_io(i,j)=dble(sp2d(i,j))
        enddo
        enddo
#else
        read(iutin,rec=11,iostat=ierr1) ((sp2d(i,j),j=1,jx),i=1,ix)
        do j=1,jx
        do i=1,ix
          f(i,j)=dble(sp2d(i,j))
        enddo
        enddo
#endif
        print*,'SNOWC'
#ifdef MPP1
        read(iutin,rec=12,iostat=ierr1) ((sp2d(i,j),j=1,mjx),i=1,ix)
        do j=1,mjx
        do i=1,ix
        do n=1,NNSG
          snowc_io(n,i,j) = dble(sp2d(i,j))
        enddo
        enddo
        enddo
#else
        read(iutin,rec=12,iostat=ierr1) ((sp2d(i,j),j=1,jx),i=1,ix)
        do j=1,jx
        do i=1,ix
        do n=1,NNSG
          snowc(n,i,j) = dble(sp2d(i,j))
        enddo
        enddo
        enddo
#endif
        if (ierr1.ne.0) then
#ifdef MPP1
          print*,'END OF FILE REACHED (SUBROUTINE PARAM)'
          print*,'  Check ibyte in regcm.param: ibyte = ',ibyte
          print*,'stop EOF (param.f)'
          call MPI_ABORT(MPI_COMM_WORLD,ierr)
#else
          print*,'END OF FILE REACHED (SUBROUTINE PARAM)'
          print*,'  Check ibyte in regcm.param: ibyte = ',ibyte
          stop 'EOF (param.f)'
#endif
        endif
#ifdef MPP1
        do j=1,mjx
          do i=1,ix
            inisrf_0(i,1,j) = ht_io(i,j)
            inisrf_0(i,2,j) = htsd_io(i,j)
            inisrf_0(i,3,j) = satbrt_io(i,j)
            inisrf_0(i,4,j) = xlat_io(i,j)
            inisrf_0(i,5,j) = xlong_io(i,j)
            inisrf_0(i,6,j) = msfx_io(i,j)
            inisrf_0(i,7,j) = msfd_io(i,j)
            inisrf_0(i,8,j) = f_io(i,j)
          enddo
          do n=1,NNSG
            do i=1,ix
              inisrf_0(i,8+n       ,j) = ht1_io(n,i,j)
              inisrf_0(i,8+NNSG+n  ,j) = satbrt1_io(n,i,j)
              inisrf_0(i,8+NNSG*2+n,j) = snowc_io(n,i,j)
            enddo
          enddo
        enddo
        do j=1,mjx
        do i=1,ix
          ht_io(i,j)=ht_io(i,j)*g
        enddo
        enddo
        endif                 ! endif (myid.eq.0) 
            
        call MPI_SCATTER(inisrf_0(1,1,1),ix*(NNSG*3+8)*jxp,
     &       MPI_DOUBLE_PRECISION,inisrf0(1,1,1),ix*(NNSG*3+8)*jxp,
     &       MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        do j=1,jxp
          do i=1,ix
            ht(i,j)     = inisrf0(i,1,j)
            htsd(i,j)   = inisrf0(i,2,j)
            satbrt(i,j) = inisrf0(i,3,j)
            xlat(i,j)   = inisrf0(i,4,j)
            xlong(i,j)  = inisrf0(i,5,j)
            msfx(i,j)   = inisrf0(i,6,j)
            msfd(i,j)   = inisrf0(i,7,j)
            f(i,j)      = inisrf0(i,8,j)
          enddo
          do n=1,NNSG
            do i=1,ix
              ht1(n,i,j)     = inisrf0(i,8+n       ,j)
              satbrt1(n,i,j) = inisrf0(i,8+NNSG+n  ,j)
              snowc(n,i,j)   = inisrf0(i,8+NNSG*2+n,j)
            enddo
          enddo
        enddo
#endif
c------invert mapscale factors:
#ifdef MPP1
        do j=1,jendl
#else
        do j=1,jx
#endif
        do i=1,ix
          msfd(i,j)=1./msfd(i,j)
        end do
        end do
#ifdef MPP1
        do j=1,jendl
#else
        do j=1,jx
#endif
        do i=1,ix
          msfx(i,j)=1./msfx(i,j)
          ht(i,j)=ht(i,j)*g
        end do
        end do
#ifdef MPP1
        if(myid .eq. 0) then
#endif
        print*,' '
        print*,'***************************************************'
        print*,'***************************************************'
        print*,'**** RegCM IS BEING RUN ON THE FOLLOWING GRID: ****'
        print*,'****     Map Projection: ',proj,'                ****'
#ifdef MPP1
        print*,'****     IX=',ix,' JX=',mjx,' KX=',kx
     &       ,'             ****'
#else
        print*,'****     IX=',ix,' JX=',jx,' KX=',kx
     &       ,'             ****'
#endif
        print*,'****     PTOP=',ptsp,' DX=',dxsp,'       ****'
        print*,'****     CLAT= ',clat,' CLON=',clon,'    ****'
        if (proj.eq.'MERCAT')
     &  print*,'****     PLAT= ',plat,' PLON=',plon,'    ****'
        print*,'***************************************************'
        print*,' '

#ifdef MPP1
        endif
        call MPI_SENDRECV(ht(1,jxp),ix,MPI_DOUBLE_PRECISION,ieast,1,
     &                    ht(1,0),ix,MPI_DOUBLE_PRECISION,iwest,1,
     &                    MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(ht(1,1),ix,MPI_DOUBLE_PRECISION,iwest,2,
     &                    ht(1,jxp+1),ix,MPI_DOUBLE_PRECISION,ieast,2,
     &                    MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(msfx(1,jxp-1),ix*2,MPI_DOUBLE_PRECISION,
     &            ieast,1,msfx(1,-1),ix*2,MPI_DOUBLE_PRECISION,
     &            iwest,1,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(msfx(1,1),ix*2,MPI_DOUBLE_PRECISION,
     &            iwest,2,msfx(1,jxp+1),ix*2,MPI_DOUBLE_PRECISION,
     &            ieast,2,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(msfd(1,jxp-1),ix*2,MPI_DOUBLE_PRECISION,
     &            ieast,1,msfd(1,-1),ix*2,MPI_DOUBLE_PRECISION,
     &            iwest,1,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(msfd(1,1),ix*2,MPI_DOUBLE_PRECISION,
     &            iwest,2,msfd(1,jxp+1),ix*2,MPI_DOUBLE_PRECISION,
     &            ieast,2,MPI_COMM_WORLD,status,ierr)
#endif
      end if
c
c-----compute dsigma and half sigma levels.
c
      do 100 k=1,kx
       dsigma(k)=sigma(k+1)-sigma(k)
       a(k)=0.5*(sigma(k+1)+sigma(k))
100    continue

      do k=1,kx
        if (a(k).lt.0.4) then
          anudg(k) = 3.
        elseif (a(k).lt.0.8) then
          anudg(k) = 2.
        else
          anudg(k) = 1.
        end if
      end do
c
c-----specify heating profile (twght) and weighted function
c     for moisture fluxes due to convection (vqflx)
c     assume base of cloud varies as  < kbase = 5,kx >
c     top  of cloud varies as  < ktop  = 1,kbase-3 >
c     exceptions to this are treated explicitly in subroutine
c     "cupara".
c
      do 108 kbase = 5,kx
      do 108 ktop  = 1,kbase-3
       do 102 k=1,kx
        twght(k,kbase,ktop) = 0.
        vqflx(k,kbase,ktop) = 0.
102    continue
c
c......get twght from 1/2 level sigma values
c
       bb = dlog(a(ktop)) + dlog(a(kbase))
       cc = dlog(a(ktop)) * dlog(a(kbase))
       ssum = 0.
       do 103 k=ktop,kbase
        xx = dlog(a(k))
        twght(k,kbase,ktop) = (xx * xx) - (bb * xx) + cc
        ssum = ssum + twght(k,kbase,ktop)*dsigma(k)
103    continue
       do 104 k=ktop,kbase
        twght(k,kbase,ktop) = twght(k,kbase,ktop)/ssum
104    continue
c
c......get vqflx from  d(w*q) / dsigma on full levels
c      do computations in p to avoid sigma=0. discontinuity
c
       xtop = dlog((100.-ptop)*sigma(ktop) + ptop)
       xbot = dlog((100.-ptop)*sigma(kbase+1) + ptop)
       bb = xtop + xbot
       cc = xtop * xbot
       vqmax = 0.
       ssum  = 0.
       xx = xtop
       yy = xbot
       wk =   (xx * xx) - (bb * xx) + cc
       qk = -((yy * yy) - (bb * yy) + cc)
       do 105 k=ktop,kbase
        xx = dlog((100.-ptop)*sigma(k+1) + ptop)
        yy = dlog((100.-ptop)*(sigma(ktop)+sigma(kbase+1)-sigma(k+1))
     &       + ptop)
        wkp1 =   (xx * xx) - (bb * xx) + cc
        qkp1 = -((yy * yy) - (bb * yy) + cc)
        vqflx(k,kbase,ktop) = - ( (wkp1*qkp1) - (wk*qk) ) / dsigma(k)
        ssum = ssum + vqflx(k,kbase,ktop)
        if (dabs(vqflx(k,kbase,ktop)) .gt. vqmax)
     &           vqmax = dabs(vqflx(k,kbase,ktop))
        wk = wkp1
        qk = qkp1
105    continue
       do 106 k=ktop,kbase
        vqflx(k,kbase,ktop) = vqflx(k,kbase,ktop) * vqrang / vqmax
106    continue
c
108   continue
c
c----calculate max no of pbl levels: kt=k at highest allowed pbl level
c-----1. caluclate sigma level at 700mb, assuming 600mb highest
c-----   allowed pbl, and 1013mb as standard surface pressure. (sigtbl)
c-----2. find first model sigma level above sigtbl.
c
      sigtbl = (70. - ptop)/(101.3 -ptop)
      kt = 1
      do k = kx,1,-1
       delsig = a(k)-sigtbl
       if(delsig.le.0.) then
         kt = k
         goto 111
       endif
      enddo
111   continue
#ifdef MPP1
      if(myid.eq.0) then
#endif
      print *, ' Index of highest allowed pbl:  kt = ', kt
      print *,' '
#ifdef MPP1
      endif
#endif
c
      if(ipptls.eq.1) then
#ifdef MPP1
        if(myid .eq. 0) then
#endif
        print*,'AUTO-CONVERSION RATE:  LAND=',qck1land
     &        ,'                      OCEAN=',qck1oce
        print*,'RELATIVE HUMIDITY THRESHOLDS:  LAND=',rh0land
     &        ,'                              OCEAN=',rh0oce
        print*,'GULTEPE FACTORS:  LAND=',gulland
     &        ,'                 OCEAN=',guloce
        print*,'MAXIMUM CLOUD COVER FOR RADIATION: ',fcmax
        print*,'MAXIMUM RELATIVE HUMIDITY: ',rhmax
        print*,'rh0 temperature threshold: ',tc0
        if (cevap.le.0.0) print*,'RAINDROP EVAPORATION NOT INCLUDED'
        if (caccr.le.0.0) print*,'RAINDROP ACCRETION NOT INCLUDED'
        print*,'Raindrop Evaporation Rate',cevap
        print*,'Raindrop Accretion Rate',caccr
#ifdef MPP1
        endif 
#endif
      endif

#ifdef MPP1
      if (myid .eq. 0) print*,' '
#else
      print*,' '
#endif

      if (icup.eq.1) then
#ifdef MPP1
        if (myid .eq. 0) then
#endif
        print*,'*********************************'
        print*,'***** Anthes-Kuo Convection *****'
        print*,'*********************************'
#ifdef MPP1
        endif
#endif
      elseif (icup.eq.2) then
        kbmax = kx
        do k=1,kx-1
          if (a(k).le.skbmax) kbmax=kx-k
        end do
#ifdef MPP1
        if (myid .eq. 0) then
#endif
        print*,'*******************************************************'
        print*,' Grell Convection:'
        print*,'   Max Shear:',shrmax
        print*,'   Min Shear:',shrmin
        print*,'   Max PPT eff:',edtmax
        print*,'   Min PPT eff:',edtmin
        print*,'   Max PPT eff(o):',edtmaxo
        print*,'   Min PPT eff(o):',edtmino
        print*,'   Max PPT eff(x):',edtmaxx
        print*,'   Min PPT eff(x):',edtminx
        print*,'   Max PBC (pbcmax):',pbcmax
        print*,'   Min Cloud Depth:',mincld
        print*,'   Max Cloud Base:',kbmax,skbmax
        print*,'   Max Heating:',htmax
        print*,'   Min Heating:',htmin
        if (igcc.eq.1) then
          print*,'   Arakawa-Schubert (1974) Closure Assumption   *****'
        elseif (igcc.eq.2) then
          print*,'   Fritsch-Chappell (1980) Closure Assumption   *****'
          print*,'     ABE removal timescale: dtauc=',dtauc
        end if
        print*,'*******************************************************'
#ifdef MPP1
        endif
        do j=1,jendx
#else
        do j=1,jlx
#endif
        do i=1,ilx
          shrmax2d(i,j)  = shrmax
          shrmin2d(i,j)  = shrmin
          edtmax2d(i,j)  = edtmax
          edtmin2d(i,j)  = edtmin
          edtmaxo2d(i,j) = edtmaxo
          edtmino2d(i,j) = edtmino
          edtmaxx2d(i,j) = edtmaxx
          edtminx2d(i,j) = edtminx
          pbcmax2d(i,j)  = pbcmax
          mincld2d(i,j)  = mincld
          kbmax2d(i,j)   = kbmax
          htmax2d(i,j)   = htmax
          htmin2d(i,j)   = htmin
          dtauc2d(i,j)   = dtauc*60.
        end do
        end do
      elseif (icup.eq.3) then
#ifdef MPP1
        if (myid .eq. 0) then
#endif
        print*,' '
        print*,' The Betts-Miller Convection scheme is not properly'
        print*,' implemented.  If want to work with it, comment out'
        print*,' the following stop statement in param.f'
#ifdef MPP1
        print*,' STOP (SUBROUTINE PARAM):  BETTS-MILLER NOT WORKING'
        endif
        call MPI_ABORT(MPI_COMM_WORLD,ierr)
#else
        stop 'STOP (SUBROUTINE PARAM):  BETTS-MILLER NOT WORKING'
#endif
      elseif (icup.eq.4) then
        cllwcv    = 0.5e-4 ! Cloud liquid water content for convective precip.
        clfrcvmax = 0.25   ! Max cloud fractional cover for convective precip.
        minorig = kx
        do k=1,kx
          if (a(k).le.minsig) minorig=kx-k
        end do
#ifdef MPP1
        if (myid .eq. 0) then
#endif
        print*,' '
        print*,'EMANUEL (1991) CONVECTION V4.3C (20 May, 2002)'
        print*,'  MIN CONVECTION ORIGIN (minsig/orig): ',minsig,minorig
        print*,'  AUTOCONVERSION THERSHOLD (elcrit): ',elcrit
        print*,'  AUTOCONVERSION THRESHOLD TO ZERO (tlcrit): ',tlcrit
        print*,'  ENTRAINMENT COEFFICIENT (entp): ',entp
        print*,'  FRACTIONAL AREA OF UNSATURATED DNDRAFT (sigd): ',sigd
        print*,'  PRECIP FRACTION OUTSIDE OF CLOUD (sigs): ',sigs
        print*,'  FALL SPEED OF RAIN (omtrain): ',omtrain
        print*,'  FALL SPEED OF SNOW (omtsnow): ',omtsnow
        print*,'  RAIN EVAPORATION COEFFICIENT (coeffr): ',coeffr
        print*,'  SNOW EVAPORATION COEFFICIENT (coeffs): ',coeffs
        print*,'  CONVECTIVE MOMENTUM TRANSPORT COEFFICIENT (cu): ',cu
        print*,'  DOWNDRAFT VELOCITY SCALE (betae): ',betae
        print*,'  MAX NEGATIVE PERTURBATION BELOW LFC (dtmax): ',dtmax
        print*,'  QUASI-EQUILIBRIUM APPROACH RATE (alphae): ',alphae
        print*,'  QUASI-EQUILIBRIUM APPROACH RATE (damp): ',damp
        print*,' '
        print*,' '
#ifdef MPP1
        endif
#endif
      end if

C Convective Cloud Cover
      afracl = 0.3 ! frac. cover for conv. precip. when dx=dxlarg
      afracs = 1.0 !   "     "    "    "      "     "   dx=dxsmal
      dlargc = 200.0
      dsmalc = 10.0
      dxtemc=dmin1(dmax1(dx,dsmalc),dlargc)
      clfrcv=afracl+(afracs-afracl)*((dlargc-dxtemc)/
     1       (dlargc-dsmalc))**2
      clfrcv = dmin1(clfrcv,clfrcvmax)
#ifdef MPP1
        if (myid .eq. 0) then
#endif
      print*,' '
      print*,'CONVECTIVE CLOUD FRACTION/WATER'
      print*,'   Maximum Convective Cloud Cover'
      print*,'     before resolution scaling: ',clfrcvmax
      print*,'   Maximum Convective Cloud Cover'
      print*,'     after resolution scaling: ',clfrcv
      print*,'   Convective Cloud Water: ',cllwcv
#ifdef MPP1
        endif
#endif

      pi = 4.*datan(1.d0)
      avt = 841.99667
      bvt = 0.8
      trel= 3000.
      g4pb = 17.837825
      g3pb = g4pb/(3.+bvt)
      g5pb = 1.8273
      n0r = 8.e6
      ppi = 1./(pi*n0r)
      vtc = avt*g4pb/6.
      prac = pi*n0r*avt*g3pb*0.25
      prec1 = 2.*pi*n0r*0.78
      prec2 = 2.*pi*n0r*0.32*avt**0.5*g5pb
c
c-----compute the vertical interpolation coefficients for t and qv.
c
      twt(1,1) = 0.
      twt(1,2) = 0.
      qcon(1) = 0.
      do 150 k=2,kx
       twt(k,1) = (sigma(k)-a(k-1))/(a(k)-a(k-1))
       twt(k,2) = 1.-twt(k,1)
       qcon(k) = (sigma(k)-a(k))/(a(k-1)-a(k))
150   continue

      chibot = 450.
      ptmb = 10.*ptop
      pz = a(1)*(1000.-ptmb)+ptmb
      if (pz .gt. chibot) stop 160
      do 160 k=1,kx
      pk = a(k)*(1000.-ptmb)+ptmb
      if (pk .gt. chibot) go to 160
      kchi = k
160   continue

c
c-----compute the k level under which the maximum equivalent potential
c     temperature will be regarded as the origin of air parcel that
c     produces cloud (used in the cumulus parameterization scheme).
c
      sig700 = (70. - ptop) / (100. - ptop)
      do 192 k=1,kx
       k700 = k
       if (sig700 .le. sigma(k+1) .and. sig700 .gt. sigma(k))
     1       go to 194
192   continue
194   continue
c
c-----specify the coefficients for sponge boundary conditions.
c
      ispgd = nspgd-1
      ispgx = nspgx-1
c.....for dot point variables:
      if (iboudy .eq. 4 ) then
       wgtd(1) = 0.
       wgtd(2) = 0.2
       wgtd(3) = 0.55
       wgtd(4) = 0.8
       wgtd(5) = 0.95
       do k=4,nspgx
         wgtd(k)=1.
       enddo
c.....for cross point variables:
       wgtx(1) = 0.
       wgtx(2) = 0.4
       wgtx(3) = 0.7
       wgtx(4) = 0.9
       do k=5,nspgx
         wgtx(k)=1.
       enddo
      endif
c
c-----specify the coefficients for nudging boundary conditions:
c
c.....for large domain:
      if (iboudy .eq. 1. or. iboudy.eq.5) then
       fnudge = 0.1/dt2
       gnudge = (dxsq/dt)/50.
      endif
      if (icup.eq.3) CALL LUTBL(PTOP)
c
c-----print out the parameters specified in the model.
c
c     print outparam
c     print physicsparam
c     print timeparam
c.....for large domain:
c
      daymax = ntimax / 1440.
#ifdef MPP1
      if (myid .eq. 0) then
#endif
      if (ibltyp .eq. 0) print 1040
      print 1070,julday,gmt,ntrad
c
      if (iboudy .eq. 0) print 1180
      if (iboudy .eq. 1) print 1182,fnudge,gnudge
      if (iboudy .eq. 2) print 1184
      if (iboudy .eq. 3) print 1190
      if (iboudy .eq. 4) print 1200
      if (iboudy .eq. 5) print 1183,fnudge,gnudge

      print 1210
c
      do 250 k=1,kx
250   print 1220,k,sigma(k),a(k),dsigma(k),twt(k,1),twt(k,2),qcon(k)
      print 1230,kxp1,sigma(kxp1)
      print 1240, daymax
      print 1250,dt
      print 1260,dx
#ifdef MPP1
      print 1270,mjx,ix
#else
      print 1270,jx,ix
#endif
      print 1280,kx
      print 1290,xkhz
      print 1300,xkhmax
#ifdef MPP1
      endif
#endif
1040  format(/'   frictionless and insulated for the lower boundary.')
1070  format('     the surface energy budget is used to calculate the gr
     1ound temperature.   julday = ',i3,'   gmt = ',f4.1/10x,
     2'the radiation is computed every ',i4,' time steps.')
C1100  format('     heat and moisture fluxes from the ground are turned o
C     1ff.')
1180  format(/'   the lateral boundary conditions are fixed.')
1182  format(/'   relaxation boundary conditions (linear method)',
     1       ' are used.  fnudge = ', e15.5,'  gnudge = ',e15.5)
1183  format(/'   relaxation boudnary conditions (exponential method)',
     1        ' are used. fnudge = ', e15.5,'  gnudge = ',e15.5)
1184  format(/'   time dependent boundary conditions are used.')
1190  format(/'   inflow/outflow boundary conditions are used.')
1200  format(/'   sponge boundary conditions are used.')
1210  format('0 k',4x,'sigma(k)',3x,'  a(k)',5x,'dsigma(k)',4x,'twt(k,1)
     1',5x,'twt(k,2)',5x,'qcon(k)'/)
1220  format(1x,i2,5x,f6.4,5x,f6.4,5x,f6.4,5x,f8.4,5x,f8.4,5x,f8.4)
1230  format(1x,i2,5x,f6.4)
1240  format(//' maximum time = ',f8.3,' days.')
1250  format(' time step = ',f7.2,' seconds')
1260  format(' dx = ',f7.0,' meters')
1270  format(' grid points (x,y) = (',i3,',',i3,')')
1280  format(' number of levels = ',i2)
1290  format(' constant hor. diff. coef. = ',e12.5,' m*m/s')
1300  format(' maximum  hor. diff. coef. = ',e12.5,' m*m/s')
c
      return
      end

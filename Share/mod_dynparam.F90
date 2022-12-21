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
!
module mod_dynparam

  use mod_stdio
  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_date
#ifndef PNETCDF
  use netcdf
#endif

  implicit none

  private
  !
  ! PARAMETER definitions
  !
  !################### GRID DIMENSION ####################################
  !

  ! Point in Y (latitude) direction

  integer(ik4) , public :: iy

  ! Point in X (longitude) direction

  integer(ik4) , public :: jx

  ! Point in vertical

  integer(ik4) , public :: kz

  ! If not 14 , 18 or 23 (precalculated), hint for custom calculation

  real(rkx) , public :: dsmax , dsmin

  ! Sub grid decomposition

  integer(ik4) , public :: nsg

  ! Dynamical core

  integer(ik4) , public :: idynamic

  ! Projection
  !
  ! One in : 'LAMCON', Lambert conformal
  !          'POLSTR', Polar stereographic
  !          'NORMER', Normal  Mercator (ROTMER w/ plat = clat
  !          'ROTMER', Rotated Mercator
  !
  character(len=6) , public :: iproj

  ! Control flag for tropical band option.

  integer(ik4) , public :: i_band
  integer(ik4) , public :: i_crm

  ! Control flag for creating bathymetry for lake model
  !    (Hostetler, etal. 1991, 1993a,b, 1995)

  logical , public :: lakedpth = .false.

  ! Control flag for crating initial soil moisture dataset

  logical , public :: lsmoist = .false.

  ! Grid point horizontal resolution in km

  real(rkx) , public :: ds

  ! Pressure of model top in cbar

  real(rkx) , public :: ptop

  ! Central latitude  of model projection in degrees, north hem. is positive
  ! Projection center location in I

  real(rkx) , public :: clat
  real(rkx) , public :: cntri

  ! Central longitude of model projection in degrees, west is negative
  ! Projection center location in J

  real(rkx) , public :: clon
  real(rkx) , public :: cntrj

  ! Pole latitude (only for rotated Mercator Proj, else set = clat)

  real(rkx) , public :: plat

  ! Pole longitude (only for rotated Mercator Proj, else set = clon)

  real(rkx) , public :: plon

  ! Lambert / Polar Cone factor

  real(rkx) , public :: xcone

  ! Lambert true latitude (low latitude side)

  real(rkx) , public :: truelatl = 0.0_rkx

  ! Lambert true latitude (high latitude side)

  real(rkx) , public :: truelath = 0.0_rkx

  ! Smoothness level

  integer(ik4) , public :: ismthlev

  !###################### DEBUG I/O control flag #########################

  ! Set amount of printout (still unused, sorry)

  integer(ik4) , public :: debug_level = 0
  integer(ik4) , public :: dbgfrq = 24

  !###################### I/O control flag ###############################

  ! Buffer Zone Depth
  ! nspgx-1,nspgd-1 represent the number of cross/dot point slices
  ! on the boundary sponge or relaxation boundary conditions.

  integer(ik4) , public :: nspgx = 12
  integer(ik4) , public :: nspgd = 12

  real(rkx) , public :: bdy_nm = -1.0_rkx ! 0.0033_rkx
  real(rkx) , public :: bdy_dm = -1.0_rkx ! 0.0001_rkx

  ! Nudge control coefficients

  real(rkx) , public :: high_nudge   = 3.0_rkx
  real(rkx) , public :: medium_nudge = 2.0_rkx
  real(rkx) , public :: low_nudge    = 1.0_rkx

  ! Type of global analysis datasets used in Pre processing
  ! One in: ECMWF,ERA40,ERAIN,EIN75,EIN15,EIM25,ERAHI,NNRP1,NNRP2,
  !         NRP2W,GFS11,FVGCM,FNEST,EH5OM

  character(len=5) , public :: dattyp
  character(len=256) , public :: cmip6_inp = &
              'https://esgf3.dkrz.de/thredds/dodsC'
  character(len=16) , public :: cmip6_model = 'MPI-ESM1-2-HR'
  character(len=12) , public :: cmip6_variant = 'r1i1p1f1'
  character(len=6) , public :: cmip6_ssp = 'ssp585'
  character(len=12) , public :: cmip6_grid = 'gn'

  !Type of Global chemistry boundary conditions
  !      MZ6HR is for MOZART 6 hourly boundary conditions
  !      MZCLM is for MOZART climatology

  character(len=5) , public :: chemtyp

  ! Type of Sea Surface Temperature used
  ! One in: GISST,OISST,OI2ST,OI_WK,OI2WK,FV_RF,FV_A2,FV_B2,EH5RF,
  !         EH5A2,EH5B1,EHA1B,ERSST,ERSKT

  character(len=5) , public :: ssttyp

  ! Land Surface Legend number

  integer(ik4) , public :: nveg

  ! Tracer parameters: number of tracers

  integer(ik4) , public :: ntr = 0  ! Total number of chemical tracers

  ! Base state atmosphere for non-hydrostatic MM5

  real(rkx) , public :: base_state_pressure ! Base state reference pressure
  real(rkx) , public :: logp_lrate  ! Logp lapse rate d(T)/d(ln P) [K/ln(Pa)]

  ! Moloch dynamical vertical profile

  real(rkx) , public :: mo_ztop = 30000.0_rkx
  real(rkx) , public :: mo_h = 8000.0_rkx
  real(rkx) , public :: mo_a0 = 0.0_rkx

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  ! End of configureation. Below this point things are
  !    calculated from above or should be considered as fixed
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer(ik4) , public :: iym1
  integer(ik4) , public :: iym2
  integer(ik4) , public :: iym3
  integer(ik4) , public :: jxm1
  integer(ik4) , public :: jxm2
  integer(ik4) , public :: jxm3
  integer(ik4) , public :: kzm1
  integer(ik4) , public :: kzm2
  integer(ik4) , public :: kzp1
  integer(ik4) , public :: kzp2
  integer(ik4) , public :: kzp3
  integer(ik4) , public :: kzp4
  integer(ik4) , public :: iysg
  integer(ik4) , public :: jxsg
  integer(ik4) , public :: iym1sg
  integer(ik4) , public :: jxm1sg
  integer(ik4) , public :: iym2sg
  integer(ik4) , public :: jxm2sg
  integer(ik4) , public :: iym3sg
  integer(ik4) , public :: jxm3sg
  integer(ik4) , public :: nnsg

  integer(ik4) , public :: njcross , njdot , njout , njoutsg
  integer(ik4) , public :: nicross , nidot , niout , nioutsg

  integer(ik4) , public :: jcross1 , icross1
  integer(ik4) , public :: jcross2 , icross2
  integer(ik4) , public :: jdot1 , idot1
  integer(ik4) , public :: jdot2 , idot2
  integer(ik4) , public :: jout1 , iout1
  integer(ik4) , public :: jout2 , iout2
  integer(ik4) , public :: joutsg1 , ioutsg1
  integer(ik4) , public :: joutsg2 , ioutsg2

  ! D stands for DOT
  ! External i (included bdy) (latitude)
  integer(ik4) , public :: ide1 , ide2
  ! External j (included bdy) (longitude)
  integer(ik4) , public :: jde1 , jde2
  ! External i SUB (included bdy) (latitude)
  integer(ik4) , public :: ide1sg , ide2sg
  ! External j SUB (included bdy) (longitude)
  integer(ik4) , public :: jde1sg , jde2sg
  ! Internal (excluded first and last line) i
  integer(ik4) , public :: idi1 , idi2
  ! Internal (excluded first and last column) j
  integer(ik4) , public :: jdi1 , jdi2
  ! Internal (excluded 2 lines and cols) i
  integer(ik4) , public :: idii1 , idii2
  ! Internal (excluded 2 lines and cols) j
  integer(ik4) , public :: jdii1 , jdii2

  ! C stands for CROSS
  ! External (included bdy) i (latitude)
  integer(ik4) , public :: ice1 , ice2
  ! External (included bdy) j (longitude)
  integer(ik4) , public :: jce1 , jce2
  ! Internal (excluded first and last line) i
  integer(ik4) , public :: ici1 , ici2
  ! Internal (excluded first and last column) j
  integer(ik4) , public :: jci1 , jci2
  ! Internal (excluded 2 lines and cols) i
  integer(ik4) , public :: icii1 , icii2
  ! Internal (excluded 2 lines and cols) j
  integer(ik4) , public :: jcii1 , jcii2

  ! Ghost points A

  integer(ik4) , public :: ici1ga , ici2ga , jci1ga , jci2ga
  integer(ik4) , public :: ice1ga , ice2ga , jce1ga , jce2ga
  integer(ik4) , public :: ide1ga , ide2ga , jde1ga , jde2ga
  integer(ik4) , public :: idi1ga , idi2ga , jdi1ga , jdi2ga

  ! Ghost points B

  integer(ik4) , public :: ici1gb , ici2gb , jci1gb , jci2gb
  integer(ik4) , public :: ice1gb , ice2gb , jce1gb , jce2gb
  integer(ik4) , public :: ide1gb , ide2gb , jde1gb , jde2gb
  integer(ik4) , public :: idi1gb , idi2gb , jdi1gb , jdi2gb

  ! Ghost points C

  integer(ik4) , public :: ici1gc , ici2gc , jci1gc , jci2gc
  integer(ik4) , public :: ice1gc , ice2gc , jce1gc , jce2gc
  integer(ik4) , public :: ide1gc , ide2gc , jde1gc , jde2gc
  integer(ik4) , public :: idi1gc , idi2gc , jdi1gc , jdi2gc

  ! SL stencil points

  integer(ik4) , public :: ici1sl , ici2sl , jci1sl , jci2sl
  integer(ik4) , public :: ice1sl , ice2sl , jce1sl , jce2sl
  integer(ik4) , public :: ide1sl , ide2sl , jde1sl , jde2sl
  integer(ik4) , public :: idi1sl , idi2sl , jdi1sl , jdi2sl

  ! J index Dot points Full Domain  = jde1 : begin , jde2 : end
  ! I index Cross points Internal Domain = ici1 : begin , ici2 : end

  !####################### MPI parameters ################################

  integer(ik4) , public :: mycomm
  integer(ik4) , public :: nproc , nprocshm
  integer(ik4) , public :: myid , myidshm
  integer(ik4) , public :: njxcpus , niycpus
  integer(ik4) , public :: iyp , jxp
  integer(ik4) , public :: iypsg , jxpsg

  !####################### MPI parameters ################################

  ! Surface minimum H2O percent to be considered water

  real(rkx) , public :: h2opct

  ! Allow water pixels to have an elevation

  logical , public :: h2ohgt

  ! Do a first resampling befor final interpolation

  logical , public :: lresamp

  ! Interpolation radius in ds unit for topography

  real(rkx) , public :: roidem

  ! Smoothing Control flag
  !     true  -> Perform extra smoothing in boundaries

  logical , public :: smthbdy

  ! Fudging for landuse and texture for grid and subgrid

  logical , public :: fudge_lnd
  logical , public :: fudge_lnd_s
  logical , public :: fudge_tex
  logical , public :: fudge_tex_s
  logical , public :: fudge_lak
  logical , public :: fudge_lak_s

  ! Terrain output files

  character(len=64) , public :: domname
  character(len=64) , public :: prestr

  ! Global Begin and End date for Input Pre processing

  type(rcm_time_and_date) , save , public :: globidate1 ! BEGIN
  type(rcm_time_and_date) , save , public :: globidate2 ! END

  ! Days per year and degrees per day

  character(len=12) , public :: calendar
  integer(ik4) , public :: ical
  real(rkx) , public :: dayspy
  real(rkx) , public :: vernal_equinox
  real(rkx) , public :: half_dayspy
  real(rkx) , public :: sixteenth_dayspy
  real(rkx) , public :: dpd

  ! Fixed dimensions

  integer(ik4) , public , parameter :: mpy = 12         ! Months per Year

  ! Number of Soil texture categories, leave it to 17

  integer(ik4) , public , parameter :: ntex = 17
  ! Should be ntex-5. Soil classes.
  integer(ik4) , public , parameter :: nats = 12

  ! Maximum number of depths in lake model

  integer(ik4) , public , parameter :: ndpmax = 200

  ! Number of bins in solar spectra

  integer(ik4) , public , parameter :: nspi = 19

  ! Number of PGW pressure levels

  integer(ik4) , public , parameter :: npgwlev = 17

#ifdef CLM45
  ! Soil layer thickness discretization (m)
  real(rkx) , public , parameter :: scalez = 0.025_rkx
  integer(ik4) , public , parameter :: num_soil_layers = 10
#else
  integer(ik4) , public , parameter :: num_soil_layers = 3
#endif

  ! Shall we use this to port?

  character(len=1) , public , parameter :: pthsep = '/'

  ! Paths

  character(len=256) , public :: dirter , inpter
  character(len=256) , public :: dirglob , inpglob
  character(len=256) , public :: dirout
  character(len=256) , public :: moist_filename
  character(len=8)   , public :: tersrc , smsrc
#ifdef NETCDF4_HDF5
  integer(ik4) , public :: iomode = ior(nf90_clobber, nf90_netcdf4)
#else
#ifndef PNETCDF
#ifdef NETCDF_CDF5
  integer(ik4) , public :: iomode = ior(nf90_clobber, nf90_cdf5)
#else
  integer(ik4) , public :: iomode = nf90_clobber
#endif
#else
  integer(ik4) , public :: iomode = 0
#endif
#endif

  integer(ik4) , public :: deflate_level = 1

  ! Model output control parameters

  logical , public :: ifsave
  logical , public :: ifatm
  logical , public :: ifshf
  logical , public :: ifrad
  logical , public :: ifsrf
  logical , public :: ifsub
  logical , public :: ifsts
  logical , public :: iflak
  logical , public :: ifopt
  logical , public :: ifchem

  real(rkx) , public :: outnwf
  real(rkx) , public :: savfrq
  real(rkx) , public :: atmfrq
  real(rkx) , public :: radfrq
  real(rkx) , public :: lakfrq
  real(rkx) , public :: subfrq
  real(rkx) , public :: srffrq
  real(rkx) , public :: chemfrq
  real(rkx) , public :: optfrq
  integer(ik4) , public :: ibdyfrq

  logical , public :: ensemble_run

  logical , public :: lperturb_topo
  logical , public :: lperturb_ts
  logical , public :: lperturb_ps
  logical , public :: lperturb_t
  logical , public :: lperturb_q
  logical , public :: lperturb_u
  logical , public :: lperturb_v

  real(rkx) , public :: perturb_frac_topo
  real(rkx) , public :: perturb_frac_ts
  real(rkx) , public :: perturb_frac_ps
  real(rkx) , public :: perturb_frac_t
  real(rkx) , public :: perturb_frac_q
  real(rkx) , public :: perturb_frac_u
  real(rkx) , public :: perturb_frac_v

#ifdef CLM45
  logical , public :: enable_megan_emission = .false.
  logical , public :: enable_urban_landunit = .true.
  logical , public :: enable_more_crop_pft = .false.
  logical , public :: enable_dv_baresoil = .false.
  logical , public :: enable_cru_precip = .false.
#endif

  public :: initparam , init_fnestparam , init_globwindow

  contains

  subroutine initparam(filename, ierr)
    implicit none
    character (len=*) , intent(in) :: filename
    integer(ik4) , intent(out) :: ierr
    integer(ik8) :: gdate1 , gdate2
    integer(ik4) :: iresult
    integer(ik4) :: ipunit

    namelist /dimparam/ iy , jx , kz , dsmax , dsmin , nsg , njxcpus , niycpus
    namelist /coreparam/ idynamic
    namelist /molochparam/ mo_a0 , mo_ztop , mo_h
    namelist /geoparam/ iproj , ds , ptop , clat , clon , plat ,    &
      plon , cntri , cntrj , truelatl , truelath , i_band , i_crm
    namelist /terrainparam/ domname , lresamp , smthbdy , lakedpth,   &
      lsmoist , fudge_lnd , fudge_lnd_s , fudge_tex , fudge_tex_s ,   &
      fudge_lak , fudge_lak_s , h2opct , h2ohgt , ismthlev , dirter , &
      inpter , moist_filename , tersrc , smsrc , roidem
    namelist /debugparam/ debug_level , dbgfrq
    namelist /boundaryparam/ nspgx , nspgd , high_nudge , &
      medium_nudge , low_nudge , bdy_nm , bdy_dm
    namelist /globdatparam/ dattyp , chemtyp, ssttyp , gdate1 , gdate2 , &
      dirglob , inpglob , calendar , ibdyfrq , ensemble_run
    namelist /cmip6param/ cmip6_inp , cmip6_model , cmip6_ssp , &
      cmip6_variant , cmip6_grid
    namelist /perturbparam/ lperturb_ts , perturb_frac_ts ,         &
      lperturb_topo , perturb_frac_topo ,         &
      lperturb_ps , perturb_frac_ps , lperturb_t , perturb_frac_t , &
      lperturb_q , perturb_frac_q , lperturb_u , perturb_frac_u ,   &
      lperturb_v , perturb_frac_v
#ifdef CLM45
    namelist /clm_regcm/ enable_megan_emission , enable_urban_landunit, &
      enable_more_crop_pft , enable_dv_baresoil , enable_cru_precip
#endif
    namelist /referenceatm/ base_state_pressure , logp_lrate

    open(newunit=ipunit, file=filename, status='old', &
         action='read', iostat=iresult)
    if ( iresult /= 0 ) then
      write (stderr,*) 'Error opening input namelist file ',trim(filename)
      ierr = 1
      return
    end if

    dsmax = 0.05_rkx
    dsmin = 0.01_rkx
    njxcpus = -1
    niycpus = -1

    rewind(ipunit)
    read(ipunit, nml=dimparam, iostat=iresult)
    if ( iresult /= 0 ) then
      write (stderr,*) 'Error reading dimparam namelist in ',trim(filename)
      ierr = 2
      return
    end if
    if ( nsg < 1 ) then
      nsg = 1
    end if

    idynamic = 1
    rewind(ipunit)
    read(ipunit, nml=coreparam, iostat=iresult)

    if ( idynamic < 1 .or. idynamic > 3 ) then
      write (stderr,*) 'Error reading coreparam namelist in ',trim(filename)
      write (stderr,*) 'Unknown dynamical core: ', idynamic
      ierr = 3
      return
    end if

    if ( idynamic == 2 ) then
      base_state_pressure = 101325.0_rkx ! Base state reference pressure
      logp_lrate = 47.70_rkx  ! Logp lapse rate d(T)/d(ln P) [K/ln(Pa)]
      rewind(ipunit)
      read(ipunit, nml=referenceatm, iostat=iresult)
    end if

    ! Moloch NH
    if ( idynamic == 3 ) then
      rewind(ipunit)
      read(ipunit, nml=molochparam, iostat=iresult)
    end if

    i_band = 0
    i_crm = 0
    cntri = -1.0_rkx
    cntrj = -1.0_rkx
    rewind(ipunit)
    read(ipunit, nml=geoparam, iostat=iresult)
    if ( iresult /= 0 ) then
      write (stderr,*) 'Error reading geoparam namelist in ',trim(filename)
      ierr = 4
      return
    end if
    if ( ds < 0.0_rkx ) then
      ds = -erkm*ds*degrad
    end if
    if ( iproj == 'LAMCON' ) then
      if ( abs(truelatl) < epsilon(1.0_rkx) .and. &
           abs(truelath) < epsilon(1.0_rkx) ) then
        write(stderr,*) 'SETTING TRUE LATITUDE TO DOMAIN CENTER!!!'
        write(stderr,*) 'If this is not what you want, set the values in ', &
          trim(filename)
        truelatl = clat
        truelath = clat
      end if
    end if
    if ( cntri < 0.0_rkx ) then
      cntri = real(iy,rkx)/d_two
    end if
    if ( cntrj < 0.0_rkx ) then
      cntrj = real(jx,rkx)/d_two
    end if

!   Ensure that band mode is active if CRM mode is active
    if ( i_crm == 1 ) then
      i_band = 1
    end if

    ! No need of ptop here.
    if ( idynamic == 3 ) then
      ptop = 0.0_rkx
    end if

!   Setup all convenience dimensions

    iym1 = iy - 1
    iym2 = iy - 2
    iym3 = iy - 3
    jxm1 = jx - 1
    jxm2 = jx - 2
    jxm3 = jx - 3
    kzm1 = kz - 1
    kzm2 = kz - 2
    kzp1 = kz + 1
    kzp2 = kz + 2
    kzp3 = kz + 3
    kzp4 = kz + 4
    iysg = iy * nsg
    jxsg = jx * nsg
    iym1sg = iym1 * nsg
    jxm1sg = jxm1 * nsg
    iym2sg = iym2 * nsg
    jxm2sg = jxm2 * nsg
    iym3sg = iym3 * nsg
    jxm3sg = jxm3 * nsg
    nnsg = nsg*nsg
    jdot1 = 1
    jdot2 = jx
    jcross1 = 1
    if ( i_band == 1 ) then
      jcross2 = jx
      jout1 = 1
      jout2 = jx
      joutsg1 = 1
      joutsg2 = jx*nsg
    else
      jcross2 = jxm1
      jout1 = 2
      jout2 = jxm2
      joutsg1 = nsg+1
      joutsg2 = jxm2*nsg
    end if
    idot1 = 1
    idot2 = iy
    icross1 = 1
    if ( i_crm == 1 ) then
      icross2 = iy
      iout1 = 1
      iout2 = iy
      ioutsg1 = 1
      ioutsg2 = iy*nsg
    else
      icross2 = iym1
      iout1 = 2
      iout2 = iym2
      ioutsg1 = nsg+1
      ioutsg2 = iym2*nsg
    end if
    njcross = jcross2-jcross1+1
    nicross = icross2-icross1+1
    njdot = jdot2-jdot1+1
    nidot = idot2-idot1+1
    njout = jout2-jout1+1
    niout = iout2-iout1+1
    njoutsg = joutsg2-joutsg1+1
    nioutsg = ioutsg2-ioutsg1+1

    nveg = 22

    if ( i_crm == 1 ) then
      iproj = 'NORMER'
      clat  =   0.0_rkx
      clon  = 180.0_rkx
    else
      if ( i_band.eq.1 ) then
        ds = real((twopi*erkm)/real(jx,rk8),rkx)
        iproj = 'NORMER'
        clat  =   0.0_rkx
        clon  = 180.0_rkx
      end if
    end if

    ! Defaults to have SAME behaviour of V3 if not specified
    inpter  = '../DATA'
    inpglob = '../DATA'
    dirter  = '../../Input'
    dirglob = '../../Input'
    moist_filename = 'moist.nc'
    tersrc = 'GMTED'
    smsrc = 'ESACCI'

    lresamp = .false.
    smthbdy = .false.
    h2ohgt = .true.
    h2opct = 50.0_rkx
    roidem = 1.5_rkx
    ismthlev = 1
    rewind(ipunit)
    read(ipunit, nml=terrainparam, iostat=iresult)
    if ( iresult /= 0 ) then
      write (stderr,*) 'Error reading terrainparam namelist in ',trim(filename)
      ierr = 5
      return
    end if

    ! Set convenient defaults for debug I/O parameters
    rewind(ipunit)
    read(ipunit, nml=debugparam, iostat=iresult)
    if ( iresult /= 0 ) then
      ! We have defaults
      continue
    end if

    ! We assume 12 points is "OK" for a resolution of 50 km.
    ! Let us assume a "default" for the selected ds, not getting less
    ! than the OK number of points. If the domain is REALLY small,
    ! use 1/4 of the overall points, or AT LEAST 3 points...
    nspgx = max(min(max(int(real(nspgx*50,rkx)/ds),nspgx),min(jx,iy)/4),3)
    nspgd = max(min(max(int(real(nspgd*50,rkx)/ds),nspgd),min(jx,iy)/4),3)
    ! Anyway the user specify this...
    rewind(ipunit)
    read(ipunit, nml=boundaryparam, iostat=iresult)
    if ( iresult /= 0 ) then
      ! We have defaults
      continue
    end if
    ! Just double check ;)
    nspgx = max(nspgx,3)
    nspgd = max(nspgd,3)

    ibdyfrq = 6 ! Convenient default
    calendar = 'gregorian'
    ensemble_run = .false.
    chemtyp = 'MZCLM'
    rewind(ipunit)
    read(ipunit, nml=globdatparam, iostat=iresult)
    if ( iresult /= 0 ) then
      write (stderr,*) 'Error reading globdatparam namelist in ',trim(filename)
      ierr = 6
      return
    end if
    if ( dattyp == 'CMIP6' ) then
      rewind(ipunit)
      read(ipunit, nml=cmip6param, iostat=iresult)
      if ( iresult /= 0 ) then
        write (stderr,*) 'Error reading cmip6param namelist in ',trim(filename)
        ierr = 7
        return
      end if
    end if
    if (calendar == 'gregorian') then
      dayspy = 365.2422_rkx
      ical = gregorian
    else if (calendar == 'noleap' .or. calendar == '365_day') then
      dayspy = 365.0_rkx
      ical = noleap
    else if (calendar == '360_day') then
      dayspy = 360.0_rkx
      ical = y360
    else
      write(stderr,*) 'No calendar specified. Assuming gregorian'
      dayspy = 365.2422_rkx
      ical = gregorian
    end if
    vernal_equinox = (80.50_rk8/dayspy)*365.0_rk8
    dpd = 360.0_rkx/dayspy
    half_dayspy = dayspy/2.0_rkx
    sixteenth_dayspy = dayspy/16.0_rkx
    globidate1 = i8wcal(gdate1,ical)
    globidate2 = i8wcal(gdate2,ical)
    if ( ensemble_run ) then
      lperturb_topo = .false.
      lperturb_ts = .false.
      lperturb_ps = .false.
      lperturb_t = .false.
      lperturb_q = .false.
      lperturb_u = .false.
      lperturb_v = .false.
      perturb_frac_topo = d_r1000
      perturb_frac_ts = d_r1000
      perturb_frac_ps = d_r1000
      perturb_frac_t = d_r1000
      perturb_frac_q = d_r1000
      perturb_frac_u = d_r1000
      perturb_frac_v = d_r1000
      rewind(ipunit)
      read(ipunit, nml=perturbparam, iostat=iresult)
      if ( iresult /= 0 ) then
        write (stderr,*) 'Error reading perturbparam namelist in ' , &
          trim(filename)
        ierr = 8
        return
      end if
    end if

#ifdef CLM45
    rewind(ipunit)
    read(ipunit, nml=clm_regcm, iostat=iresult)
    if ( iresult /= 0 ) then
      ! No error here. we have defualts. ;)
      continue
    end if
#endif

    close(ipunit)
    ierr = 0
    return
  end subroutine initparam

  subroutine init_fnestparam(filename,coarse_outdir,coarse_domname)
    implicit none
    character(len=*) , intent(in) :: filename
    character(len=256) , intent(out) :: coarse_outdir , coarse_domname
    integer(ik4) :: iresult
    integer(ik4) :: ipunit
    namelist /fnestparam/ coarse_outdir , coarse_domname

    coarse_outdir = '        '
    coarse_domname = '        '
    open(newunit=ipunit, file=filename, status='old', &
                 action='read', iostat=iresult)
    if ( iresult /= 0 ) then
      write (stderr,*) 'Error opening input namelist file ',trim(filename)
      return
    end if
    read(ipunit, nml=fnestparam, iostat=iresult)
    if ( iresult /= 0 ) then
      write (stdout,*) 'fnestparam not present: Using defaults'
    end if
    close(ipunit)
  end subroutine init_fnestparam

  subroutine init_globwindow(filename,lat0,lon0,lat1,lon1)
    implicit none
    character(len=*) , intent(in) :: filename
    real(rkx) , intent(out) :: lat0 , lat1 , lon0 , lon1
    integer(ik4) :: iresult
    integer(ik4) :: ipunit
    namelist /globwindow/ lat0 , lat1 , lon0 , lon1

    lat0 = 0.0_rkx
    lon0 = 0.0_rkx
    lat1 = 0.0_rkx
    lon1 = 0.0_rkx

    open(newunit=ipunit, file=filename, status='old', &
                 action='read', iostat=iresult)
    if ( iresult /= 0 ) then
      write (stderr,*) 'Error opening input namelist file ',trim(filename)
      return
    end if
    read(ipunit, nml=globwindow, iostat=iresult)
    if ( iresult /= 0 ) then
      write (stdout,*) 'Globwindow namelist not present: Assuming Global'
    end if
    close(ipunit)
  end subroutine init_globwindow

end module mod_dynparam
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

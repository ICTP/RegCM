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

  use mod_constants
  use mod_date

  public
!
! PARAMETER definitions
!
  integer , parameter :: ipunit = 255
!
!################### GRID DIMENSION ####################################
!

! Point in Y (latitude) direction

  integer :: iy

! Point in X (longitude) direction

  integer :: jx

! Point in vertical

  integer :: kz

! If not 14 , 18 or 23 (precalculated), hint for custom calculation

  real(8) :: dsmax , dsmin

! Sub grid decomposition

  integer :: nsg

! Projection
!
! One in : 'LAMCON', Lambert conformal
!          'POLSTR', Polar stereographic
!          'NORMER', Normal  Mercator (ROTMER w/ plat = clat
!          'ROTMER', Rotated Mercator
!
  character(6) :: iproj
 
! Control flag for tropical band option.
 
  integer :: i_band

! Control flag for lake model (Hostetler, etal. 1991, 1993a,b, 1995)
 
  logical :: lakedpth

! Grid point horizontal resolution in km

  real(8) :: ds

! Pressure of model top in cbar

  real(8) :: ptop

! Central latitude  of model domain in degrees, north hem. is positive

  real(8) :: clat

! Central longitude of model domain in degrees, west is negative

  real(8) :: clon

! Pole latitude (only for rotated Mercator Proj, else set = clat)

  real(8) :: plat

! Pole longitude (only for rotated Mercator Proj, else set = clon)

  real(8) :: plon

! Lambert true latitude (low latitude side)

  real(8) :: truelatl

! Lambert true latitude (high latitude side)

  real(8) :: truelath

!###################### I/O control flag ###############################

! Number of bytes in reclen. Usually 4

  integer :: ibyte

! Set amount of printout (still unused, sorry)

  integer :: debug_level
  integer :: dbgfrq

!###################### I/O control flag ###############################

! Buffer Zone Depth
! nspgx-1,nspgd-1 represent the number of cross/dot point slices
! on the boundary sponge or relaxation boundary conditions.
!
  integer :: nspgx
  integer :: nspgd

! Nudge control coefficients
  real(8) :: high_nudge
  real(8) :: medium_nudge
  real(8) :: low_nudge

! Number od split exp modes

  integer :: nsplit

! Type of global analysis datasets used in Pre processing
!
! One in: ECMWF,ERA40,ERAIN,EIN75,EIN15,EIM25,ERAHI,NNRP1,NNRP2,
!         NRP2W,GFS11,FVGCM,FNEST,EH5OM
!
  character(5) :: dattyp

! Type of Sea Surface Temperature used
!
! One in: GISST,OISST,OI2ST,OI_WK,OI2WK,FV_RF,FV_A2,FV_B2,EH5RF,
!         EH5A2,EH5B1,EHA1B,ERSST,ERSKT
!
  character(5) :: ssttyp

! Land Surface Legend number

  integer :: nveg

! Aerosol dataset used
!
! One in : AER00D0 -> Neither aerosol, nor dust used
!          AER01D0 -> Biomass, SO2 + BC + OC, no dust
!          AER10D0 -> Anthropogenic, SO2 + BC + OC, no dust
!          AER11D0 -> Anthropogenic+Biomass, SO2 + BC + OC, no dust
!          AER00D1 -> No aerosol, with dust
!          AER01D1 -> Biomass, SO2 + BC + OC, with dust
!          AER10D1 -> Anthropogenic, SO2 + BC + OC, with dust
!          AER11D1 -> Anthropogenic+Biomass, SO2 + BC + OC, with dust

  character(7) :: aertyp

! Tracer parameters: number of tracers and bins number for dust and sea salt

  integer :: ntr   ! Total number of chemical tracers

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! End of configureation. Below this point things are
!    calculated from above or should be considered as fixed
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  integer :: iym1
  integer :: iym2
  integer :: iym3
  integer :: jxm1
  integer :: jxm2
  integer :: jxm3
  integer :: kzm1
  integer :: kzm2
  integer :: kzp1
  integer :: kzp2
  integer :: kzp3
  integer :: kzp4
  integer :: iysg
  integer :: jxsg
  integer :: iym1sg
  integer :: jxm1sg
  integer :: iym2sg
  integer :: jxm2sg
  integer :: iym3sg
  integer :: jxm3sg
  integer :: nnsg
  integer :: nspgv
  integer :: nspgp
!
  integer :: njcross , njdot , njout
  integer :: nicross , nidot , niout
!
  integer :: jcross1 , icross1
  integer :: jcross2 , icross2
  integer :: jdot1 , idot1
  integer :: jdot2 , idot2
  integer :: jout1 , iout1
  integer :: jout2 , iout2
!
  ! D stands for DOT
  integer :: ide1 , ide2 ! External i (included bdy) (latitude)
  integer :: jde1 , jde2 ! External j (included bdy) (longitude)
  integer :: idi1 , idi2 ! Internal (excluded first and last line) i
  integer :: jdi1 , jdi2 ! Internal (excluded first and last column) j
  integer :: idii1 , idii2 ! Internal (excluded 2 lines and cols) i
  integer :: jdii1 , jdii2 ! Internal (excluded 2 lines and cols) j

  ! C stands for CROSS
  integer :: ice1 , ice2 ! External (included bdy) i (latitude)
  integer :: jce1 , jce2 ! External (included bdy) j (longitude)
  integer :: ici1 , ici2 ! Internal (excluded first and last line) i
  integer :: jci1 , jci2 ! Internal (excluded first and last column) j
  integer :: icii1 , icii2 ! Internal (excluded 2 lines and cols) i
  integer :: jcii1 , jcii2 ! Internal (excluded 2 lines and cols) j

  ! J index Dot points Full Domain  = jde1 : begin , jde2 : end
  ! I index Cross points Internal Domain = ici1 : begin , ici2 : end

!####################### MPI parameters ################################

  integer :: mycomm
  integer :: nproc
  integer :: myid
  integer :: jxp
  integer :: jxpsg

!####################### MPI parameters ################################

! Surface minimum H2O percent to be considered water

  real(8) :: h2opct

! Resolution of the global terrain and landuse data be used
!
!     Use 60, for  1  degree resolution
!         30, for 30 minutes resolution
!         10, for 10 minutes resolution
!          5, for  5 minutes resolution
!          3, for  3 minutes resolution
!          2, for  2 minutes resolution

  integer :: ntypec

! Same for subgrid (Used only if nsg > 1)

  integer :: ntypec_s

! Smoothing Control flag
!
!     true  -> Perform extra smoothing in boundaries

  logical :: smthbdy

! Fudging for landuse and texture for grid and subgrid

  logical :: fudge_lnd
  logical :: fudge_lnd_s
  logical :: fudge_tex
  logical :: fudge_tex_s
  logical :: fudge_lak
  logical :: fudge_lak_s

! Terrain output files

  character(64) :: domname

! Global Begin and End date for Input Pre processing

  type(rcm_time_and_date) , save :: globidate1 ! BEGIN
  type(rcm_time_and_date) , save :: globidate2 ! END

! Days per year and degrees per day

  character(12) :: calendar
  integer :: ical
  real(8) :: dayspy
  real(8) :: dpd

! Fixed dimensions

  integer , parameter :: numsts = 10
  integer , parameter :: numbat = 24 + numsts
  integer , parameter :: numsub = 16

  integer , parameter :: mpy = 12         ! Months per Year

! Number of Soil texture categories, leave it to 17

  integer , parameter :: ntex = 17 
  integer , parameter :: nats = 12 ! Should be ntex-5. Soil classes.

! Maximum number of depths in lake model

  integer , parameter :: ndpmax = 400 ! This means 400 m max depth

! Number of bins in solar spectra

  integer , parameter :: nspi = 19

! Shall we use this to port?

  character(1), parameter :: pthsep = '/'

! Paths

  character(256) :: dirter , inpter
  character(256) :: dirglob , inpglob
  character(256) :: dirout
  character(256) :: dirclm

! Model output control parameters

  logical :: ifsave
  real(8) :: savfrq

  logical :: ifatm
  real(8) :: atmfrq

  logical :: ifrad
  real(8) :: radfrq

  logical :: ifsrf
  logical :: ifsub
  logical :: ifsts
  logical :: iflak
  real(8) :: lakfrq
  real(8) :: srffrq

  logical :: ifchem
  real(8) :: chemfrq

  integer :: ibdyfrq

  contains

  subroutine initparam(filename, ierr)
    implicit none
    character (len=*) , intent(in) :: filename
    integer , intent(out) :: ierr
    integer :: gdate1 , gdate2

    namelist /geoparam/ iproj , ds , ptop , clat , clon , plat ,    &
                   plon , truelatl, truelath , i_band
    namelist /terrainparam/ domname , ntypec , ntypec_s ,           &
                  smthbdy , lakedpth, fudge_lnd , fudge_lnd_s ,     &
                  fudge_tex , fudge_tex_s , fudge_lak, fudge_lak_s ,&
                  h2opct , dirter , inpter
    namelist /dimparam/ iy , jx , kz , dsmax , dsmin , nsg
    namelist /ioparam/ ibyte
    namelist /debugparam/ debug_level , dbgfrq
    namelist /boundaryparam/ nspgx , nspgd , high_nudge , &
                   medium_nudge , low_nudge
    namelist /modesparam/ nsplit
    namelist /globdatparam/ dattyp , ssttyp , gdate1 , gdate2 , &
                   dirglob , inpglob , calendar , ibdyfrq
    namelist /aerosolparam/ aertyp , ntr
! , nbin, sbin

    open(ipunit, file=filename, status='old', &
                 action='read', err=100)
!
    dsmax = d_zero
    dsmin = d_zero

    read(ipunit, dimparam, err=101)

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
    iym1sg = (iy-1) * nsg
    jxm1sg = (jx-1) * nsg
    iym2sg = (iy-2) * nsg
    jxm2sg = (jx-2) * nsg
    iym3sg = (iy-3) * nsg
    jxm3sg = (jx-3) * nsg

    jdot1 = 1
    jdot2 = jx
    jcross1 = 1
#ifdef BAND
    jcross2 = jx
    jout1 = 1
    jout2 = jx
#else
    jout1 = 2
    jcross2 = jxm1
    jout2 = jxm2
#endif
    idot1 = 1
    idot2 = iy
    icross1 = 1
    icross2 = iym1
    iout1 = 2
    iout2 = iym2
    njcross = jcross2-jcross1+1
    nicross = icross2-icross1+1
    njdot = jdot2-jdot1+1
    nidot = idot2-idot1+1
    njout = jout2-jout1+1
    niout = iout2-iout1+1

    nnsg = nsg*nsg
    nveg = 22

    i_band = 0

    read(ipunit, geoparam, err=102)

    if ( i_band.eq.1 ) then
      ds = (2.0D0*mathpi*erkm)/dble(jx)
      iproj = 'NORMER'
      clat  =   0.0D0
      clon  = 180.0D0
    end if

    ! Defaults to have SAME behaviour of V3 if not specified
    inpter  = '../DATA'
    inpglob = '../DATA'
    dirter  = '../../Input'
    dirglob = '../../Input'

    read(ipunit, terrainparam, err=103)

    ! Set convenient defaults for I/O parameters
    ibyte = 4
    read(ipunit, ioparam, err=104)
    dbgfrq = 3600
    read(ipunit, debugparam, err=105)

    high_nudge = 3.0D0
    medium_nudge = 2.0D0
    low_nudge = 1.0D0
    read(ipunit, boundaryparam, err=106)

    nspgv = (nspgd+nspgx)*8 + 8
    nspgp = nspgx*4

    read(ipunit, modesparam, err=107)

    ibdyfrq = 6 ! Convenient default
    calendar = 'gregorian'
    read(ipunit, globdatparam, err=109)
    if (calendar == 'gregorian') then
      dayspy = 365.2422D+00
      ical = gregorian
    else if (calendar == 'noleap' .or. calendar == '365_day') then
      dayspy = 365.0D+00
      ical = noleap
    else if (calendar == '360_day') then
      dayspy = 360.0D+00
      ical = y360
    else
      dayspy = 365.2422D+00
      ical = gregorian
    end if
    dpd = 360.0D0/dayspy
    globidate1 = gdate1
    globidate2 = gdate2
    call setcal(globidate1,ical)
    call setcal(globidate2,ical)

    ntr  = 0
    read(ipunit, aerosolparam, err=111)

    ierr = 0
    return

  100   write ( 6, * ) 'Cannot read namelist file ', trim(filename)
    ierr = 1 
    return 
  101   write ( 6, * ) 'Cannot read namelist stanza: dimparam       ',  &
        & trim(filename)
    ierr = 1
    return
  102   write ( 6, * ) 'Cannot read namelist stanza: geoparam       ',  &
        & trim(filename)
    ierr = 1
    return
  103   write ( 6, * ) 'Cannot read namelist stanza: terrainparam   ',  &
        & trim(filename)
    ierr = 1
    return
  104   write ( 6, * ) 'Cannot read namelist stanza: ioparam        ',  &
        & trim(filename)
    ierr = 1
    return
  105   write ( 6, * ) 'Cannot read namelist stanza: debugparam     ',  &
        & trim(filename)
    ierr = 1
    return
  106   write ( 6, * ) 'Cannot read namelist stanza: boundaryparam  ',  &
        & trim(filename)
    ierr = 1
    return
  107   write ( 6, * ) 'Cannot read namelist stanza: modesparam     ',  &
        & trim(filename)
    ierr = 1
    return
  109   write ( 6, * ) 'Cannot read namelist stanza: globdatparam   ',  &
        & trim(filename)
    ierr = 1
    return
  111   write ( 6, * ) 'Cannot read namelist stanza: aereosolparam  ',  &
        & trim(filename)
    ierr = 1

  end subroutine initparam

  subroutine init_globwindow(lat0,lon0,lat1,lon1)
    implicit none
    real(8) , intent(out) :: lat0 , lat1 , lon0 , lon1
    namelist /globwindow/ lat0 , lat1 , lon0 , lon1

    lat0 = 0.0D0
    lon0 = 0.0D0
    lat1 = 0.0D0
    lon1 = 0.0D0

    rewind(ipunit)
    read(ipunit, globwindow,err=101)
    return
  101   print *, 'Globwindow not present: Assuming Global data input'
    return
  end subroutine init_globwindow

  subroutine init_outparam
    implicit none

    namelist /outparam/ ifsave , savfrq , ifatm , atmfrq ,       &
      ifrad , radfrq , ifsrf , ifsub , iflak , ifsts , srffrq , &
      lakfrq , ifchem , chemfrq

    read(ipunit, outparam, err=100)
    return

  100   write ( 6, * ) 'Cannot read namelist stanza: outparam'

  end subroutine init_outparam

end module mod_dynparam

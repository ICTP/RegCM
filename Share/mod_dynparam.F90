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

  use mod_realkinds
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

  real(dp) :: dsmax , dsmin

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

! Control flag for creating bathymetry for lake model
!    (Hostetler, etal. 1991, 1993a,b, 1995)
 
  logical :: lakedpth = .false.

! Control flag for crating teture dataset for aerosol dust
!

  logical :: ltexture = .false.

! Grid point horizontal resolution in km

  real(dp) :: ds

! Pressure of model top in cbar

  real(dp) :: ptop

! Central latitude  of model domain in degrees, north hem. is positive

  real(dp) :: clat

! Central longitude of model domain in degrees, west is negative

  real(dp) :: clon

! Pole latitude (only for rotated Mercator Proj, else set = clat)

  real(dp) :: plat

! Pole longitude (only for rotated Mercator Proj, else set = clon)

  real(dp) :: plon

! Lambert / Polar Cone factor

  real(dp) :: xcone

! Lambert true latitude (low latitude side)

  real(dp) :: truelatl

! Lambert true latitude (high latitude side)

  real(dp) :: truelath

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
  real(dp) :: high_nudge
  real(dp) :: medium_nudge
  real(dp) :: low_nudge

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
!
! Tracer parameters: number of tracers

  integer :: ntr = 0  ! Total number of chemical tracers

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

  ! Global reference in global grid jx*iy of dot points
  ! The CROSS grid is contained within
  integer :: global_istart
  integer :: global_iend
  integer :: global_jstart
  integer :: global_jend

!####################### MPI parameters ################################

  integer :: mycomm
  integer :: nproc
  integer :: myid
  integer :: iyp , jxp
  integer :: iypsg , jxpsg

!####################### MPI parameters ################################

! Surface minimum H2O percent to be considered water

  real(dp) :: h2opct

! Allow water pixels to have an elevation

  logical :: h2ohgt 
!
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
  real(dp) :: dayspy
  real(dp) :: half_dayspy
  real(dp) :: sixteenth_dayspy
  real(dp) :: dpd

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
  real(dp) :: savfrq

  logical :: ifatm
  real(dp) :: atmfrq

  logical :: ifrad
  real(dp) :: radfrq

  logical :: ifsrf
  logical :: ifsub
  logical :: ifsts
  logical :: iflak
  real(dp) :: lakfrq
  real(dp) :: srffrq

  logical :: ifchem
  real(dp) :: chemfrq

  integer :: ibdyfrq

  contains

  subroutine initparam(filename, ierr)
    implicit none
    character (len=*) , intent(in) :: filename
    integer , intent(out) :: ierr
    integer :: gdate1 , gdate2

    namelist /geoparam/ iproj , ds , ptop , clat , clon , plat ,    &
                   plon , truelatl, truelath , i_band
    namelist /terrainparam/ domname , smthbdy , ltexture , lakedpth,  &
                  fudge_lnd , fudge_lnd_s , fudge_tex , fudge_tex_s , &
                  fudge_lak,  fudge_lak_s , h2opct , h2ohgt , dirter , inpter
    namelist /dimparam/ iy , jx , kz , dsmax , dsmin , nsg
    namelist /ioparam/ ibyte
    namelist /debugparam/ debug_level , dbgfrq
    namelist /boundaryparam/ nspgx , nspgd , high_nudge , &
                   medium_nudge , low_nudge
    namelist /globdatparam/ dattyp , ssttyp , gdate1 , gdate2 , &
                   dirglob , inpglob , calendar , ibdyfrq
! , nbin, sbin

    open(ipunit, file=filename, status='old', &
                 action='read', err=100)
!
    dsmax = d_zero
    dsmin = d_zero

    read(ipunit, dimparam, err=101)

    i_band = 0

    read(ipunit, geoparam, err=102)

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
    else
      jcross2 = jxm1
      jout1 = 2
      jout2 = jxm2
    end if
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

    nveg = 22

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

    h2ohgt = .false.
    h2opct = 50.0D0
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

    ! Removed modesparam. It must not be changed, so avoid cluttering.
    nsplit = 2

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
    half_dayspy = dayspy/2.0D0
    sixteenth_dayspy = dayspy/16.0D0
    globidate1 = gdate1
    globidate2 = gdate2
    call setcal(globidate1,ical)
    call setcal(globidate2,ical)

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
  109   write ( 6, * ) 'Cannot read namelist stanza: globdatparam   ',  &
        & trim(filename)
    ierr = 1
    return

  end subroutine initparam

  subroutine init_globwindow(lat0,lon0,lat1,lon1)
    implicit none
    real(dp) , intent(out) :: lat0 , lat1 , lon0 , lon1
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

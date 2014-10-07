#include <misc.h>
#include <preproc.h>

module controlMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: controlMod
!
! !DESCRIPTION:
! Module which initializes run control variables. The following possible
! namelist variables are set default values and possibly read in on startup
!
! === define run =======================
!
!    o caseid     = 256 character case name
!    o ctitle     = 256 character case title
!    o nsrest     = integer flag. 0: initial run. 1: restart: 3: branch
!
! === model time =======================
!
!    o dtime      = integer model time step (s)
!    o calendar   = Calendar to use in date calculations.
!                  'no_leap' (default) or 'gregorian'
!    o start_ymd  = Starting date for run encoded in yearmmdd format.
!                   Default value is read from initial conditions file.
!    o start_tod  = Starting time of day for run in seconds since 0Z.
!                   Default value is read from initial conditions file.
!    o stop_ymd   = Stopping date for run encoded in yearmmdd format.
!                   No default.
!    o stop_tod   = Stopping time of day for run in seconds since 0Z.
!                   Default: 0.
!    o nelapse    = nnn, Specify the ending time for the run as an interval
!                   starting at the current time in either timesteps
!                   (if positive) or days (if negative).
!                   Either nestep or (stop_ymd,stop_tod) take precedence.
!    o nestep     = nnnn, Specify the ending time for the run as an interval
!                   starting at (start_ymd,start_tod) in either timesteps
!                   (if positive) or days (if negative).
!                   (stop_ymd,stop_tod) takes precedence if set.
!    o ref_ymd    = Reference date for time coordinate encoded in yearmmdd format.
!                   Default value is start_ymd.
!    o ref_tod    = Reference time of day for time coordinate in seconds since 0Z.
!                   Default value is start_tod.
!
! === input data ===
!
!    o finidat         = 256 character initial conditions file name
!    o fsurdat         = 256 character surface data file name
!    o flndtopo        = 256 character land topography file name
!    o fatmgrid        = 256 character atmosphere grid data file name
!    o fatmlndfrc      = 256 character landfrac (on atm grid) file name
!    o fatmtopo        = 256 character atmosphere topography file name
!    o fndepdat        = 254 character nitrogen deposition data file name (netCDF)
!    o fpftcon         = 256 character data file with PFT physiological constants
!    o frivinp_rtm     = 256 character input data file for rtm
!    o nrevsn          = 256 character restart file name for use with branch run
! !!!!!!!! abt rcm below
!    o mksrf_fsoicol   = 256 character soil color file name
!    o mksrf_furban    = 256 character urban fraction file name
!    o mksrf_flai      = 256 character leaf area index file name
!    o mksrf_fsoitex   = 256 character soil texture data file name
!    o mksrf_fglacier  = 256 character glacier fraction file name
!    o mksrf_fvegtyp   = 256 character vegetation type file name
!    o mksrf_flanwat   = 256 character lake and water file name
!    o mksrf_fisop     = 256 character Isoprene file name
!    o mksrf_fmbo      = 256 character Methylbutenol file name
!    o mksrf_fbpin     = 256 character B-pinene file name
!    o mksrf_fapin     = 256 character A-pinene file name
!    o mksrf_fmyrc     = 256 character Myrcene file name
!    o mksrf_flimo     = 256 character Limonene file name
!    o mksrf_fsabi     = 256 character Sabinene file name
!    o mksrf_fco       = 256 character Carbonmonoxide file name
!    o mksrf_fomtp     = 256 character Other Monoterpenes file name
!    o mksrf_facar     = 256 character A-3carnene file name
!    o mksrf_fbcar     = 256 character B-Caryophyllene file name
!    o mksrf_ffarn     = 256 character Farnicene file name
!    o mksrf_fmeoh     = 256 character Methanol file name
!    o mksrf_fmeth     = 256 character Methane file name
!    o mksrf_facto     = 256 character Acetone file name
!    o mksrf_facta     = 256 character Acetaldehyde/Ethanol file name
!    o mksrf_fform     = 256 character Formic Acid Acetic Acid and Formaldehyde file name
!    o mksrf_focim     = 256 character Ocimene file name
!    o mksrf_fno       = 256 character NO2,N2O,NH3 file name
!    o mksrf_fosqt     = 256 character Other Sesquiterpenes file name
!    o mksrf_fmax      = 256 character saturation fraction file name
!    o mksrf_offline_fnavyoro = 256 character land fraction and oro file name
! !!!!!!!! abt rcm above
!
! === offline forcing data ===
!
!    o offline_atmdir  = 256 character directory for input atm data files (can be Mass Store)
!
! === history and restart files ===
!
!    o hist_ndens    = integer, can have value of 1 (nc_double) or 2 (nf_float)
!    o hist_dov2xy   = true if want grid-average history field (false = vector)
!    o hist_nhtfrq   = integer history interval (+ = iterations,  - = hours, 0=monthly ave)
!    o hist_mfilt    = integer number of time samples per history file
!    o hist_fincl1   = 10 character name of fields for first  auxillary history file
!    o hist_fincl2   = 10 character name of fields for second auxillary history file
!    o hist_fincl3   = 10 character name of fields for first  auxillary history file
!    o hist_fincl4   = 10 character name of fields for second auxillary history file
!    o hist_fincl5   = 10 character name of fields for first  auxillary history file
!    o hist_fincl6   = 10 character name of fields for second auxillary history file
!    o hist_fexcl1   = 8  character name of fields for first  auxillary history file
!    o hist_fexcl2   = 8  character name of fields for second auxillary history file
!    o hist_fexcl3   = 8  character name of fields for first  auxillary history file
!    o hist_fexcl4   = 8  character name of fields for second auxillary history file
!    o hist_fexcl5   = 8  character name of fields for first  auxillary history file
!    o hist_fexcl6   = 8  character name of fields for second auxillary history file
!    o hist_crtinic  = 8  character frequency to generate initial dataset
!                         ['6-HOURLY','DAILY','MONTHLY','YEARLY','NONE']
!    o rest_flag     = logical, turns off restart file writing [.TRUE., .FALSE.]
!    o rpntpath      = 256 character full UNIX pathname of the local restart pointer file.
!                      This file must exist when the model is restarted.
!                      This file is overwritten every time new restart data files are output.
!
! === long term archiving =====
!
!    o archive_dir = 256 character long term archive directory (can be MSS directory)
!    o mss_irt     = integer mass store retention period (days)
!    o mss_wpass   = 8 character mass store write password for output data sets
!
! === model physics ===
!
!    o irad         = integer solar radiation frequency (+ = iteration. - = hour)
!    o wrtdia       = true if want output written
!    o csm_doflxave = true => flux averaging is to be performed (only used for csm mode)
!
! === rtm control variables ===
!
!    o rtm_nsteps  = if > 1, average rtm over rtm_nsteps time steps
!    o nsegspc     = number of segments per clump for decomposition
!
! When coupled to CAM: base calendar info, nstep, nestep, nsrest, and time
! step are input to the land model from CAM. The values in the clm_inparm namelist
! are not used.
!
! !USES:
  use shr_kind_mod , only : r8 => shr_kind_r8, SHR_KIND_CL
  use shr_sys_mod  , only : shr_sys_getenv
  use clm_varpar   , only : maxpatch_pft, numpft
  use clm_varctl
  use spmdMod
  use decompMod    , only : clump_pproc
  use histFileMod  , only : max_tapes, max_namlen, &
                            hist_empty_htapes, hist_dov2xy, &
                            hist_avgflag_pertape, hist_type1d_pertape, &
                            hist_nhtfrq, hist_ndens, hist_mfilt, &
                            hist_fincl1, hist_fincl2, hist_fincl3, &
                            hist_fincl4, hist_fincl5, hist_fincl6, &
                            hist_fexcl1, hist_fexcl2, hist_fexcl3, &
                            hist_fexcl4, hist_fexcl5, hist_fexcl6
  use restFileMod  , only : rest_flag
  use shr_const_mod, only : SHR_CONST_CDAY
  use abortutils   , only : endrun
  use clm_varsur   , only : r2coutfrq
  use mod_clm
  use mod_runparams , only : dirclm
!
! !PUBLIC TYPES:
  implicit none
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: control_setNL ! Set namelist filename
  public :: control_init  ! initial run control information
  public :: control_print ! print run control information
!
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! PRIVATE TYPES:
! Namelist variables only used locally
  character(len=256) :: rpntpath                         ! full UNIX pathname of restart pointer file
  character(len=  7) :: runtyp(4)                        ! run type
  character(len=SHR_KIND_CL) :: NLFilename = 'lnd.stdin' ! Namelist filename
#if (defined _OPENMP)
   integer, external :: omp_get_max_threads  ! max number of threads that can execute
                                             ! concurrently in a single parallel region
#endif
!-----------------------------------------------------------------------

contains

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: control_setNL
!
! !INTERFACE:
  subroutine control_setNL( NLfile )

    implicit none
!
! !DESCRIPTION:
! Set the namelist filename to use
!
!
! !ARGUMENTS:
  character(len=*), intent(IN) :: NLFile ! Namelist filename
!
! !REVISION HISTORY:
! Created by Erik Kluzek
!
!EOP
!
! !LOCAL VARIABLES:
    character(len=32) :: subname = 'control_setNL'  ! subroutine name
    logical :: lexist                               ! File exists

    ! Error checking...
    if ( len_trim(NLFile) == 0 )then
       call endrun( subname//' error: nlfilename entered is not set' )
    end if
    inquire (file = trim(NLFile), exist = lexist)
    if ( .not. lexist )then
       call endrun( subname//' error: NLfilename entered does NOT exist:'//trim(NLFile) )
    end if
    if ( len_trim(NLFile) > len(NLFilename) )then
       call endrun( subname//' error: entered NLFile is too long' )
    end if
    ! Set the filename
    NLFilename = NLFile
  end subroutine control_setNL

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: control_init
!
! !INTERFACE:
  subroutine control_init( CCSMInit )
!
! !DESCRIPTION:
! Initialize CLM run control information
!
! !USES:
    use clm_time_manager
#if (defined CASA)
    use CASAMod          , only : lnpp, lalloc, q10, spunup, fcpool
#endif
    use shr_inputinfo_mod, only : shr_inputInfo_initType,       &
                                  shr_inputInfo_initGetData,    &
                                  shr_inputInfo_initIsBranch,   &
                                  shr_inputInfo_initIsContinue, &
                                  shr_inputInfo_initIsStartup
    use fileutils        , only : getavu, relavu
    use shr_string_mod   , only : shr_string_getParentDir

    implicit none
!
! !ARGUMENTS:
    type(shr_InputInfo_initType), intent(in), optional :: CCSMInit   ! Input CCSM info

    include 'netcdf.inc'
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    character(len=256) :: homedir   ! full UNIX filepath name of home directory
    character(len=256) :: logid     ! logid part of file path name
    character(len=256) :: cap       ! upper case logid
    character(len=  1) :: ctmp      ! character temporary
    character(len=256) :: drvarchdir! driver archive directory
    integer :: i,j,n                ! loop indices
    integer :: iundef               ! integer undefined value
    real(r8):: rundef               ! real undefined value
    integer :: ierr                 ! error code
    integer :: unitn                ! unit for namelist file
    character(len=32) :: subname = 'control_init'  ! subroutine name
!------------------------------------------------------------------------

    ! ----------------------------------------------------------------------
    ! Namelist Variables
    ! ----------------------------------------------------------------------

    ! clm time manager info

#if (defined OFFLINE) || (defined COUP_CSM)
    namelist /clm_inparm/  &
         ctitle, caseid, nsrest,  &
         calendar, nelapse, nestep, start_ymd, start_tod,  &
         stop_ymd, stop_tod, ref_ymd, ref_tod

    ! Archive options
    namelist /clm_inparm/ archive_dir, mss_wpass, mss_irt, &
         brnch_retain_casename
#endif

    ! clm input datasets

    namelist / clm_inparm/ &
	 dtime

    namelist /clm_inparm/  &
         finidat, fsurdat, fatmgrid, fatmlndfrc, fatmtopo, flndtopo, &
         fpftcon, frivinp_rtm,  &
         fpftdyn, fndepdat, fndepdyn, nrevsn, offline_atmdir

!!!!!!!! abt rcm below
    namelist /clm_inparm/  &
         mksrf_fvegtyp, mksrf_fsoicol, mksrf_fsoitex, mksrf_flanwat, &
         mksrf_fglacier, mksrf_flai,  &
         mksrf_furban, mksrf_fisop, mksrf_offline_fnavyoro, mksrf_fmax, &
         mksrf_fbpin,mksrf_fapin,mksrf_fmbo, &
         mksrf_flimo,mksrf_fsabi,mksrf_fmyrc, &
         mksrf_fco,mksrf_focim,mksrf_facar, &
         mksrf_fbcar,mksrf_ffarn,mksrf_fomtp, &
         mksrf_fosqt,mksrf_fmeoh,mksrf_facto, &
         mksrf_facta,mksrf_fno,mksrf_fform, &
         mksrf_fmeth

!!!!!!! abt rcm above


    ! clm history, restart, archive options

    namelist /clm_inparm/  &
         hist_empty_htapes, hist_dov2xy, &
         hist_avgflag_pertape, hist_type1d_pertape, &
         hist_nhtfrq,  hist_ndens, hist_mfilt, &
         hist_fincl1,  hist_fincl2, hist_fincl3, &
         hist_fincl4,  hist_fincl5, hist_fincl6, &
         hist_fexcl1,  hist_fexcl2, hist_fexcl3, &
         hist_fexcl4,  hist_fexcl5, hist_fexcl6, &
         hist_crtinic, rpntpath, rest_flag

    ! clm bgc info

#if (defined CASA)
    namelist /clm_inparm/  &
         lnpp, lalloc, q10, spunup, fcpool
#endif
    namelist /clm_inparm / &
         co2_type

    ! clm other options

    namelist /clm_inparm/  &
         clump_pproc, irad, wrtdia, csm_doflxave, rtm_nsteps, pertlim, &
         create_crop_landunit, nsegspc

    ! ----------------------------------------------------------------------
    ! Default values
    ! ----------------------------------------------------------------------

    if (masterproc) then
       write(6,*) 'Attempting to initialize run control settings .....'
    endif

    runtyp(0 + 1) = 'initial'
    runtyp(1 + 1) = 'restart'
    runtyp(3 + 1) = 'branch '

    iundef = -9999999
    rundef = -9999999._r8

    ! control variables

    caseid  = 'clmoutput'
    ctitle  = 'clmoutput'
    nsrest  = iundef

    ! initial data

    fsurdat     = ' '
    fatmgrid    = ' '
    fatmlndfrc  = ' '
    fatmtopo    = ' '
    flndtopo    = ' '
    fndepdat    = ' '
    fndepdyn    = ' '
    finidat     = ' '
    fpftcon     = trim(dirclm)//pthsep//'pft-physiology.c070207'
    frivinp_rtm = ' '
    fpftdyn     = ' '
    nrevsn      = ' '

!!!!!!! abt rcm below
    mksrf_fvegtyp   = trim(dirclm)//pthsep//trim(domname)//'_RCMpft.nc'
    mksrf_fsoicol   = trim(dirclm)//pthsep//trim(domname)//'_RCMsoicol_clm2.nc'
    mksrf_fsoitex   = trim(dirclm)//pthsep//trim(domname)//'_RCMsoitex.10level.nc'
    mksrf_flanwat   = trim(dirclm)//pthsep//trim(domname)//'_RCMlanwat.nc'
    mksrf_fglacier  = trim(dirclm)//pthsep//trim(domname)//'_RCMglacier.nc'
    mksrf_flai      = trim(dirclm)//pthsep//trim(domname)//'_RCMlai.nc'
    mksrf_furban    = trim(dirclm)//pthsep//trim(domname)//'_RCMurban.nc'
    mksrf_fisop     = trim(dirclm)//pthsep//trim(domname)//'_RCMiso.nc'
    mksrf_fmbo      = trim(dirclm)//pthsep//trim(domname)//'_RCMmbo.nc'
    mksrf_fbpin     = trim(dirclm)//pthsep//trim(domname)//'_RCMpinb.nc'
    mksrf_fapin     = trim(dirclm)//pthsep//trim(domname)//'_RCMpina.nc'
    mksrf_flimo     = trim(dirclm)//pthsep//trim(domname)//'_RCMlimo.nc'
    mksrf_fsabi     = trim(dirclm)//pthsep//trim(domname)//'_RCMsabi.nc'
    mksrf_fmyrc     = trim(dirclm)//pthsep//trim(domname)//'_RCMmyrc.nc'
    mksrf_fco       = trim(dirclm)//pthsep//trim(domname)//'_RCMco.nc'
    mksrf_focim     = trim(dirclm)//pthsep//trim(domname)//'_RCMocim.nc'
    mksrf_facar     = trim(dirclm)//pthsep//trim(domname)//'_RCMacar.nc'
    mksrf_fbcar     = trim(dirclm)//pthsep//trim(domname)//'_RCMbcar.nc'
    mksrf_ffarn     = trim(dirclm)//pthsep//trim(domname)//'_RCMfarn.nc'
    mksrf_fomtp     = trim(dirclm)//pthsep//trim(domname)//'_RCMomtp.nc'
    mksrf_fosqt     = trim(dirclm)//pthsep//trim(domname)//'_RCMosqt.nc'
    mksrf_fmeoh     = trim(dirclm)//pthsep//trim(domname)//'_RCMmeoh.nc'
    mksrf_facto     = trim(dirclm)//pthsep//trim(domname)//'_RCMacto.nc'
    mksrf_facta     = trim(dirclm)//pthsep//trim(domname)//'_RCMacta.nc'
    mksrf_fno       = trim(dirclm)//pthsep//trim(domname)//'_RCMno.nc'
    mksrf_fform     = trim(dirclm)//pthsep//trim(domname)//'_RCMform.nc'
    mksrf_fmeth     = trim(dirclm)//pthsep//trim(domname)//'_RCMmeth.nc'
    mksrf_fmax      = trim(dirclm)//pthsep//trim(domname)//'_RCMfmax.nc'
    mksrf_offline_fnavyoro = trim(dirclm)//pthsep//trim(domname)//'_RCMnavyoro_20min.nc'
!!!!!!! abt rcm above


    ! offline mode

    offline_atmdir   = ' '

    ! landunit generation

    create_crop_landunit = .false.
    if (maxpatch_pft == numpft+1) then
       allocate_all_vegpfts = .true.
    else
       allocate_all_vegpfts = .false.
    end if

    ! bgc

    co2_type = 'constant'

    ! long term archive settings

    archive_dir = ' '
    mss_irt = 0
    mss_wpass = ' '

    ! history file variables

    hist_crtinic = 'NONE'
    rpntpath = 'not_specified'

    ! other namelist variables

    irad = -1
    wrtdia = .true.
    csm_doflxave = .true.
    pertlim = 0._r8
    single_column=.false.
    scmlat=-999.
    scmlon=-999.
    nsegspc = 20

#if (defined RTM)
    ! If rtm_nsteps is not set in the namelist then
    ! will be given default value below

    rtm_nsteps = -999
#endif

#if (defined CASA)
    lnpp = 2
    lalloc = 1
    q10 = 2.0_r8          ! set Q10 to 2.0  03/11/19
    spunup = 0
    fcpool = ' '
#endif

    ! Set clumps per procoessor

    clump_pproc = 1
#if (defined _OPENMP)
    clump_pproc = omp_get_max_threads()
#else
#if (defined UNICOSMP)
#if (defined SSP)
    clump_pproc = 1
#else
    clump_pproc = 1 ! 4 when using CSDs in driver.F90; 1 otherwise
#endif
#endif
#endif

    if (masterproc) then

       ! ----------------------------------------------------------------------
       ! Read namelist from standard input.
       ! ----------------------------------------------------------------------

!abt commented out reading the clm namelist
!       if ( len_trim(NLFilename) == 0  )then
!          call endrun( subname//' error: nlfilename not set' )
!       end if
!       unitn = getavu()
!       write(6,*) 'Read in clm_inparm namelist from: ', trim(NLFilename)
!       open( unitn, file=trim(NLFilename), status='old' )
!       ierr = 1
!       do while ( ierr /= 0 )
!          read(unitn, clm_inparm, iostat=ierr)
!          if (ierr < 0) then
!             call endrun( subname//' encountered end-of-file on namelist read' )
!          endif
!       end do
!       call relavu( unitn )
!abt commented out reading the namelist

! abt ------- force RegCM namelist variables ----------
       nsrest    = r2cnsrest
       dtime     = r2cdtime
       start_ymd = r2cstart_ymd
       start_tod = r2cstart_tod
       stop_ymd  = r2cstop_ymd
       stop_tod  = r2cstop_tod
       calendar  = r2cclndr
       wrtdia    = r2cwrtdia
       irad      = r2cirad
       mss_irt   = r2cmss_irt
       if(r2coutfrq .gt. 0) then
         hist_nhtfrq = -r2coutfrq
       else
         hist_nhtfrq = r2coutfrq
       end if
       hist_mfilt = 1

      print *, 'nsrest = ',nsrest,' dtime = ',dtime
      print *,'start_ymd = ',start_ymd,' start_tod = ',start_tod
      print *,'stop_ymd = ',stop_ymd,' stop_tod = ',stop_tod
      print *,'calendar = ',calendar,'wrtdia = ',wrtdia
      print *,'irad = ',irad
      print *, 'hist_nhtfrq = ', hist_nhtfrq

       ! ----------------------------------------------------------------------
       ! Consistency checks on input namelist.
       ! ----------------------------------------------------------------------

       ! Consistency settings for co2 type

       if (co2_type /= 'constant' .and. co2_type /= 'prognostic' .and. co2_type /= 'diagnostic') then
          write(6,*)'co2_type = ',co2_type,' is not supported'
          write(6,*)'choices are constant, prognostic or diagnostic'
          call endrun()
       end if
#if (defined OFFLINE)
       if (co2_type /= 'constant') then
          write(6,*)'co2_type = ',co2_type,' is not supported in offline mode'
          write(6,*)'choice is only constant'
          call endrun()
       end if
#endif

       ! Consistency settings for dynamic land use, etc.

       if (fpftdyn /= ' ' .and. create_crop_landunit) then
          write(6,*)'dynamic landuse is currently not supported with create_crop_landunit option'
          call endrun()
       end if
       if (fpftdyn /= ' ') then
#if (defined DGVM)
          write(6,*)'dynamic landuse is currently not supported with DGVM option'
          call endrun()
#elif (defined CASA)
          write(6,*)'dynamic landuse is currently not supported with CASA option'
          call endrun()
#endif
       end if

#if (defined SEQ_MCT) || (defined SEQ_ESMF)
       ! Override select set of namelist values with sequential driver input

       if ( .not. present(CCSMInit) )then
          call endrun( subname//' error CCSMInit not present but is '// &
                       'required when linking with CAM' )
       end if
       call shr_inputInfo_initGetData( CCSMInit, case_name=caseid,        &
                                       case_desc=ctitle, mss_irt=mss_irt, &
                                       mss_wpass=mss_wpass,               &
                                       brnch_retain_casename=brnch_retain_casename,  &
                                       archive_dir=drvarchdir,&
                                       single_column=single_column ,      &
				       scmlat=scmlat,scmlon=scmlon)
       archive_dir = shr_string_getParentDir( drvarchdir )//'/lnd/'
       if (      shr_inputInfo_initIsStartup(  CCSMInit ) )then
          nsrest = 0
       else if ( shr_inputInfo_initIsContinue( CCSMInit ) )then
          nsrest = 1
       else if ( shr_inputInfo_initIsBranch(   CCSMInit ) )then
          nsrest = 3
       end if
#if (defined RTM) || (defined DGVM)
       if (is_perpetual()) then
          write(6,*)'RTM or DGVM cannot be defined in perpetual mode'
          call endrun()
       end if
#endif
       if (is_perpetual()) then
          if (finidat == ' ') then
             write(6,*)'must specify initial dataset for perpetual mode'
             call endrun()
          end if
       end if
#endif

       ! Check that if archive directory not input in namelist, set default from caseid

       if (archive_dir == ' ') then
          logid  = ' '
          call shr_sys_getenv('LOGNAME', logid, ierr)
          if (ierr /= 0) then
             write (6,*) 'error: logname not defined'
             call endrun
          end if
          cap = ' '
          do i = 1, len_trim(logid)
             cap(i:i) = logid(i:i)
             ctmp = cap(i:i)
             if (ichar(logid(i:i))>=97 .and. ichar(logid(i:i))<=122) then
                cap(i:i) = char(ichar(ctmp) - 32)
             endif
          end do
          archive_dir = 'mss:/' // trim(cap) // '/csm/' // trim(caseid) // '/lnd'
       end if

#if (defined RTM)
       ! If rtm_nsteps was not entered in the namelist, give it the following default value

       if (rtm_nsteps == -999) then
          rtm_nsteps = (3600*3)/dtime ! 3 hours
       endif
#endif

       ! Check that hist_type_1d is not set for primary tape

       if (hist_type1d_pertape(1) /= ' ') then
          write(6,*)'CONTROL_INIT error: hist_type1d_pertape can only be set for tapes 2-6'
          call endrun()
       end if

       ! Check on run type

       if (nsrest == iundef) then
          write(6,*) 'error: must set nsrest'
          call endrun()
       end if
       if (nsrest == 3 .and. nrevsn == ' ') then
          write(6,*) 'error: need to set restart data file name'
          call endrun()
       end if

       ! Check on offline mode

#if (defined OFFLINE)
!       if (offline_atmdir == ' ') then
!          write(6,*)'error: atmos  input data file must be specified'; call endrun()
!       end if
#endif

#if (defined OFFLINE)


!!!!!!!! abt rcm below
       if (mksrf_fvegtyp == ' ') then
          write(6,*)'error: vegtype  input data file must be specified'; call endrun()
       end if
       if (mksrf_fsoitex == ' ') then
          write(6,*)'error: soil texture  input data file must be specified'; call endrun()
       end if
       if (mksrf_fsoicol == ' ') then
          write(6,*)'error: soil color  input data file must be specified'; call endrun()
       end if
       if (mksrf_flanwat == ' ') then
          write(6,*)'error: lake/water  input data file must be specified'; call endrun()
       end if
       if (mksrf_fglacier == ' ') then
          write(6,*)'error: glacier  input data file must be specified'; call endrun()
       end if
       if (mksrf_furban == ' ') then
          write(6,*)'error: urban  input data file must be specified'; call endrun()
       end if
       if (mksrf_flai == ' ') then
          write(6,*)'error: leaf area index  input data file must be specified'; call endrun()
       end if
       if (mksrf_fmax == ' ') then
          write(6,*)'error: fmax  input data file must be specified'; call endrun()
       end if
       if (mksrf_offline_fnavyoro == ' ') then
          write(6,*)'error: navy file  input data file must be specified'; call endrun()
       end if
#endif
!!!!!!!!! abt rcm above


       ! Check on ccsm mode

#if (defined COUP_CSM)
       if (csm_doflxave .and. irad ==1 ) then
          write(6,*)'error: irad must be greater that one if', &
            ' flux averaging option is enabled'
          call endrun
       end if
#endif

       ! Check on nitrogen deposition dataset

       if (fndepdat /= ' ' .and. fndepdyn /= ' ') then
          write(6,*)'namelist error: only one of fndepdat or fndepdyn can be defined'
          call endrun()
       end if

       ! Model physics

       if (irad < 0) irad = nint(-irad*3600._r8/dtime)

       ! History and restart files

       mss_irt = min(mss_irt,1825)

       do i = 1, max_tapes
          if (hist_nhtfrq(i) == 0) then
             hist_mfilt(i) = 1
          else if (hist_nhtfrq(i) < 0) then
             hist_nhtfrq(i) = nint(-hist_nhtfrq(i)*SHR_CONST_CDAY/(24._r8*dtime))
          endif
       end do

       if (rpntpath == 'not_specified') then
!abt commented          call shr_sys_getenv('HOME', homedir, ierr)
!abt commented          rpntpath = trim(homedir)//'/lnd.'//trim(caseid)//'.rpointer'
!abt use Input folder as default location for clm restart pointer file
         rpntpath = './lnd.'//trim(caseid)//'.rpointer'
       endif

       if (nsrest == 0) nrevsn = ' '
       if (nsrest == 1) nrevsn = 'set by restart pointer file file'

#if (defined DGVM)
       hist_crtinic = 'YEARLY'
#else
       if (trim(hist_crtinic) /= 'MONTHLY'  .and. trim(hist_crtinic) /= 'YEARLY' .and. &
            trim(hist_crtinic) /= '6-HOURLY' .and. trim(hist_crtinic) /= 'DAILY'  ) then
          hist_crtinic = 'NONE'
       endif
#endif

       ! Split the full pathname of the restart pointer file into a
       ! directory name and a file name

       rpntdir = ' '
       rpntfil = ' '
       do n = len_trim(rpntpath),1,-1
          if (rpntpath(n:n) ==  '/') then
             rpntdir = rpntpath(1:n-1)
             rpntfil = rpntpath(n+1:len_trim(rpntpath))
             go to 100
          endif
       enddo
       rpntdir = '.'        ! no "/" found, set path = "."
       rpntfil = rpntpath   ! no "/" found, use whole input string.
100    continue

    endif   ! end of if-masterproc if-block

    ! ----------------------------------------------------------------------
    ! Broadcast all control information if appropriate
    ! ----------------------------------------------------------------------

    call control_spmd()

    if (masterproc) then
       write(6,*) 'Successfully initialized run control settings'
       write(6,*)
    endif

  end subroutine control_init


!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: control_spmd
!
! !INTERFACE:
  subroutine control_spmd()
!
! !DESCRIPTION:
! Distribute namelist data all processors. The cpp SPMD definition
! provides for the funnelling of all program i/o through the master
! processor. Processor 0 either reads restart/history data from the
! disk and distributes it to all processors, or collects data from
! all processors and writes it to disk.
!
! !USES:
!
    use clm_time_manager
#if (defined CASA)
    use CASAMod, only : lnpp, lalloc, q10, spunup
#endif
    use spmdMod, only : mpicom
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer ier       !error code
!-----------------------------------------------------------------------

    ! run control variables

    call mpi_bcast (caseid, len(caseid), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (ctitle, len(ctitle), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (nsrest,           1, MPI_INTEGER  , 0, mpicom, ier)

    call mpi_bcast (dtime    , 1, MPI_INTEGER  , 0, mpicom, ier)

#if (defined OFFLINE) || (defined COUP_CSM)
    call mpi_bcast (nestep   , 1, MPI_INTEGER  , 0, mpicom, ier)
    call mpi_bcast (nelapse  , 1, MPI_INTEGER  , 0, mpicom, ier)
    call mpi_bcast (start_ymd, 1, MPI_INTEGER  , 0, mpicom, ier)
    call mpi_bcast (start_tod, 1, MPI_INTEGER  , 0, mpicom, ier)
    call mpi_bcast (stop_ymd , 1, MPI_INTEGER  , 0, mpicom, ier)
    call mpi_bcast (stop_tod , 1, MPI_INTEGER  , 0, mpicom, ier)
    call mpi_bcast (ref_ymd  , 1, MPI_INTEGER  , 0, mpicom, ier)
    call mpi_bcast (ref_tod  , 1, MPI_INTEGER  , 0, mpicom, ier)
    call mpi_bcast (calendar ,len(calendar), MPI_CHARACTER, 0, mpicom, ier)
#endif

    ! initial file variables

    call mpi_bcast (nrevsn  , len(nrevsn)  , MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (finidat , len(finidat) , MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fsurdat , len(fsurdat) , MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fatmgrid, len(fatmgrid), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fatmlndfrc,len(fatmlndfrc),MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fatmtopo, len(fatmtopo) ,MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (flndtopo, len(flndtopo) ,MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fndepdat, len(fndepdat), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fndepdyn, len(fndepdyn), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (fpftdyn , len(fpftdyn) , MPI_CHARACTER, 0, mpicom, ier)

!!!!!!!! abt rcm below
    call mpi_bcast (mksrf_offline_fnavyoro , len(mksrf_offline_fnavyoro) , &
                    MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (mksrf_fsoitex, len(mksrf_fsoitex), MPI_CHARACTER, &
                    0, mpicom, ier)
    call mpi_bcast (mksrf_fsoicol, len(mksrf_fsoicol), MPI_CHARACTER, &
                    0, mpicom, ier)
    call mpi_bcast (mksrf_flanwat, len(mksrf_flanwat), MPI_CHARACTER, &
                    0, mpicom, ier)
    call mpi_bcast (mksrf_fglacier,len(mksrf_fglacier), MPI_CHARACTER, &
                    0, mpicom, ier)
    call mpi_bcast (mksrf_furban, len(mksrf_furban), MPI_CHARACTER, 0, &
                    mpicom, ier)
    call mpi_bcast (mksrf_flai,   len(mksrf_flai)  , MPI_CHARACTER, 0, &
                    mpicom, ier)
#if (defined VOC)
    call mpi_bcast (mksrf_fisop,   len(mksrf_fisop)  , MPI_CHARACTER, 0, &
                    mpicom, ier)
    call mpi_bcast (mksrf_fbpin,   len(mksrf_fbpin)  , MPI_CHARACTER, 0, &
                    mpicom, ier)
    call mpi_bcast (mksrf_fapin,   len(mksrf_fapin)  , MPI_CHARACTER, 0, &
                    mpicom, ier)
    call mpi_bcast (mksrf_fmbo ,   len(mksrf_fmbo)   , MPI_CHARACTER, 0, &
                    mpicom, ier)
    call mpi_bcast (mksrf_fmyrc,   len(mksrf_fmyrc)  , MPI_CHARACTER, 0, &
                    mpicom, ier)
    call mpi_bcast (mksrf_fsabi,   len(mksrf_fsabi)  , MPI_CHARACTER, 0, &
                    mpicom, ier)
    call mpi_bcast (mksrf_flimo,   len(mksrf_flimo)  , MPI_CHARACTER, 0, &
                    mpicom, ier)
    call mpi_bcast (mksrf_fco ,   len(mksrf_fco)   , MPI_CHARACTER, 0, &
                    mpicom, ier)
    call mpi_bcast (mksrf_focim,   len(mksrf_focim)  , MPI_CHARACTER, 0, &
                    mpicom, ier)
    call mpi_bcast (mksrf_facar,   len(mksrf_facar)  , MPI_CHARACTER, 0, &
                    mpicom, ier)
    call mpi_bcast (mksrf_fomtp,   len(mksrf_fomtp)  , MPI_CHARACTER, 0, &
                    mpicom, ier)
    call mpi_bcast (mksrf_ffarn ,   len(mksrf_ffarn)   , MPI_CHARACTER, 0, &
                    mpicom, ier)
    call mpi_bcast (mksrf_fbcar,   len(mksrf_fbcar)  , MPI_CHARACTER, 0, &
                    mpicom, ier)
    call mpi_bcast (mksrf_fosqt,   len(mksrf_fosqt)  , MPI_CHARACTER, 0, &
                    mpicom, ier)
    call mpi_bcast (mksrf_fmeoh,   len(mksrf_fmeoh)  , MPI_CHARACTER, 0, &
                    mpicom, ier)
    call mpi_bcast (mksrf_facto ,   len(mksrf_facto)   , MPI_CHARACTER, 0, &
                    mpicom, ier)
    call mpi_bcast (mksrf_fmeth,   len(mksrf_fmeth)  , MPI_CHARACTER, 0, &
                    mpicom, ier)
    call mpi_bcast (mksrf_fno,   len(mksrf_fno)  , MPI_CHARACTER, 0, &
                    mpicom, ier)
    call mpi_bcast (mksrf_facta,   len(mksrf_facta)  , MPI_CHARACTER, 0, &
                    mpicom, ier)
    call mpi_bcast (mksrf_fform ,   len(mksrf_fform)   , MPI_CHARACTER, 0, &
                    mpicom, ier)
#endif
    call mpi_bcast (mksrf_fvegtyp, len(mksrf_fvegtyp) , MPI_CHARACTER, &
                    0, mpicom, ier)
    call mpi_bcast (mksrf_fmax,   len(mksrf_fmax)  , MPI_CHARACTER, 0, &
                    mpicom, ier)
!!!!!!! abt rcm above



#if (defined RTM)
    call mpi_bcast (frivinp_rtm, len(frivinp_rtm), MPI_CHARACTER, 0, mpicom, ier)
#endif

    ! Landunit generation

    call mpi_bcast(create_crop_landunit, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast(allocate_all_vegpfts, 1, MPI_LOGICAL, 0, mpicom, ier)

    ! BGC

    call mpi_bcast (co2_type, len(co2_type), MPI_CHARACTER, 0, mpicom, ier)

    ! physics variables

    call mpi_bcast (irad        , 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (csm_doflxave, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (rtm_nsteps  , 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (nsegspc     , 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (wrtdia      , 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (single_column,1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (scmlat,       1, MPI_REAL8,   0, mpicom, ier)
    call mpi_bcast (scmlon,       1, MPI_REAL8,   0, mpicom, ier)

    ! history file variables

    call mpi_bcast (hist_empty_htapes, 1, MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (hist_dov2xy, size(hist_dov2xy), MPI_LOGICAL, 0, mpicom, ier)
    call mpi_bcast (hist_nhtfrq, size(hist_nhtfrq), MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (hist_mfilt, size(hist_mfilt), MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (hist_ndens, size(hist_ndens), MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (hist_crtinic, len(hist_crtinic), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_avgflag_pertape, size(hist_avgflag_pertape), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_type1d_pertape, max_namlen*size(hist_type1d_pertape), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl1, max_namlen*size(hist_fexcl1), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl2, max_namlen*size(hist_fexcl2), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl3, max_namlen*size(hist_fexcl3), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl4, max_namlen*size(hist_fexcl4), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl5, max_namlen*size(hist_fexcl5), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fexcl6, max_namlen*size(hist_fexcl6), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl1, (max_namlen+2)*size(hist_fincl1), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl2, (max_namlen+2)*size(hist_fincl2), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl3, (max_namlen+2)*size(hist_fincl3), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl4, (max_namlen+2)*size(hist_fincl4), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl5, (max_namlen+2)*size(hist_fincl5), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (hist_fincl6, (max_namlen+2)*size(hist_fincl6), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (rest_flag, 1, MPI_LOGICAL, 0, mpicom, ier)

    ! restart file variables

    call mpi_bcast (rpntpath, len(rpntpath), MPI_CHARACTER, 0, mpicom, ier)

    ! clump decomposition variables

    call mpi_bcast (clump_pproc, 1, MPI_INTEGER, 0, mpicom, ier)

    ! long term archiving variables

    call mpi_bcast (mss_irt, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (mss_wpass, len(mss_wpass), MPI_CHARACTER, 0, mpicom, ier)
    call mpi_bcast (archive_dir, len(archive_dir), MPI_CHARACTER, 0, mpicom, ier)

    ! error growth perturbation limit
    call mpi_bcast (pertlim, 1, MPI_REAL8, 0, mpicom, ier)

#if (defined CASA)
    call mpi_bcast (lnpp  , 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (lalloc, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (spunup, 1, MPI_INTEGER, 0, mpicom, ier)
    call mpi_bcast (q10   , 1, MPI_REAL8  , 0, mpicom, ier)
#endif

  end subroutine control_spmd

!------------------------------------------------------------------------
!BOP
!
! !IROUTINE: control_print
!
! !INTERFACE:
  subroutine control_print ()
!
! !DESCRIPTION:
! Write out run control variables
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!
! !LOCAL VARIABLES:
    integer i  !loop index
!------------------------------------------------------------------------

    write(6,*) 'define run:'
    write(6,*) '   run type              = ',runtyp(nsrest+1)
    write(6,*) '   case title            = ',trim(ctitle)
    write(6,*) 'input data files:'
    write(6,*) '   PFT physiology = ',trim(fpftcon)
    if (fsurdat == ' ') then
       write(6,*) '   fsurdat, surface dataset not set'
    else
       write(6,*) '   surface data   = ',trim(fsurdat)
    end if
    if (flndtopo == ' ') then
       write(6,*) '   flndtopo not set'
    else
       write(6,*) '   land topographic data = ',trim(flndtopo)
    end if
    if (fatmgrid == ' ') then
       write(6,*) '   fatmgrid not set, using fsurdat'
       fatmgrid = fsurdat
       write(6,*) '   atm grid data  = ',trim(fatmgrid)
    else
       write(6,*) '   atm grid data  = ',trim(fatmgrid)
    end if
    if (fatmlndfrc == ' ') then
       write(6,*) '   fatmlndfrc not set, using fatmgrid'
       fatmlndfrc = fatmgrid
       write(6,*) '   land frac data = ',trim(fatmlndfrc)
    else
       write(6,*) '   land frac data = ',trim(fatmlndfrc)
    end if
    if (fatmtopo == ' ') then
       write(6,*) '   fatmtopo not set'
    else
       write(6,*) '   atm topographic data = ',trim(fatmtopo)
    end if
    if (fndepdat == ' ') then
        write(6,*) '   NOT using input data for nitrogen deposition'
    else
        write(6,*) '   nitrogen deposition data = ',trim(fndepdat)
    endif
    if (fndepdyn == ' ') then
        write(6,*) '   NOT using dynamic input data for nitrogen deposition'
    else
        write(6,*) '   dynamic nitrogen deposition data = ',trim(fndepdyn)
    endif

!!!!!!!! abt rcm below
    if (mksrf_fvegtyp == ' ') then
       write(6,*) '   mksrf_fvegtyp not set'
    else
       write(6,*) '   vegetation type data = ',trim(mksrf_fvegtyp)
    end if

    if (mksrf_fmax == ' ') then
       write(6,*) '   mksrf_fmax not set'
    else
       write(6,*) '   fmax data = ',trim(mksrf_fmax)
    end if

    if (mksrf_fsoitex == ' ') then
       write(6,*) '   mksrf_fsoitex not set'
    else
       write(6,*) '   soil texture data = ',trim(mksrf_fsoitex)
    end if

    if (mksrf_fsoicol == ' ') then
       write(6,*) '   mksrf_fsoicol not set'
    else
       write(6,*) '   soil color data = ',trim(mksrf_fsoicol)
    end if

    if (mksrf_flanwat == ' ') then
       write(6,*) '   mksrf_flanwat not set'
    else
       write(6,*) '   lake/water data = ',trim(mksrf_flanwat)
    end if

    if (mksrf_fglacier == ' ') then
       write(6,*) '   mksrf_fglacier not set'
    else
       write(6,*) '   glacier data = ',trim(mksrf_fglacier)
    end if

    if (mksrf_furban == ' ') then
       write(6,*) '   mksrf_furban not set'
    else
       write(6,*) '   urban data = ',trim(mksrf_furban)
    end if

    if (mksrf_flai == ' ') then
       write(6,*) '   mksrf_flai not set'
    else
       write(6,*) '   LAI data = ',trim(mksrf_flai)
    end if
#if (defined VOC)
    if (mksrf_fisop == ' ') then
       write(6,*) '   mksrf_fisop not set'
    else
       write(6,*) '   Isoprene data = ',trim(mksrf_fisop)
    end if
    if (mksrf_fmbo == ' ') then
       write(6,*) '   mksrf_fmbo not set'
    else
       write(6,*) '   Methylbutenol data = ',trim(mksrf_fmbo)
    end if
    if (mksrf_fbpin == ' ') then
       write(6,*) '   mksrf_fbpin not set'
    else
       write(6,*) '   B-pinene data = ',trim(mksrf_fbpin)
    end if
    if (mksrf_fapin == ' ') then
       write(6,*) '   mksrf_fapin not set'
    else
       write(6,*) '   A-pinene data = ',trim(mksrf_fapin)
    end if
    if (mksrf_flimo == ' ') then
       write(6,*) '   mksrf_flimo not set'
    else
       write(6,*) '   Limonene data = ',trim(mksrf_flimo)
    end if
#endif
    if (mksrf_offline_fnavyoro == ' ') then
       write(6,*) '   mksrf_offline_fnavyoro not set'
    else
       write(6,*) '   navyoro data = ',trim(mksrf_offline_fnavyoro)
    end if
!!!!!!!!!! abt rcm above




    if (nsrest == 0 .and. finidat == ' ') write(6,*) '   initial data created by model'
    if (nsrest == 0 .and. finidat /= ' ') write(6,*) '   initial data   = ',trim(finidat)
    if (nsrest /= 0) write(6,*) '   restart data   = ',trim(nrevsn)
#if (defined OFFLINE)
    if (offline_atmdir /= ' ') then
       write(6,*) '   atmospheric forcing data    = ',trim(offline_atmdir)
    end if
#elif (defined SEQ_MCT) || (defined SEQ_ESMF)
    write(6,*) '   atmospheric forcing data is from sequential ccsm model'
#elif (defined COUP_CSM)
    write(6,*) '   atmospheric forcint data is from ccsm flux coupler'
#endif
#if (defined RTM)
    if (frivinp_rtm /= ' ') write(6,*) '   RTM river data       = ',trim(frivinp_rtm)
#endif
    if (mss_irt /= 0) then
       write(6,*) 'Mass store control values'
       write(6,*)'   mass store path                    = ',trim(archive_dir)
       write(6,*)'   mass store retention (days)        = ',mss_irt
       write(6,*)'   mass store write password          = ',mss_wpass
    endif
    write(6,*) 'Restart parameters:'
    write(6,*)'   restart pointer file directory     = ',trim(rpntdir)
    write(6,*)'   restart pointer file name          = ',trim(rpntfil)
    if (hist_crtinic == 'MONTHLY') then
       write(6,*)'initial datasets will be written monthly'
    else if (hist_crtinic == 'YEARLY') then
       write(6,*)'initial datasets will be written yearly'
    else if (hist_crtinic == 'DAILY') then
       write(6,*)'initial datasets will be written daily'
    else if (hist_crtinic == '6-HOURLY') then
       write(6,*)'initial datasets will be written 6-hourly'
    else
       write(6,*)'initial datasets will not be produced'
    endif
    write(6,*) 'model physics parameters:'
#if (defined PERGRO)
    write(6,*) '   flag for random perturbation test is set'
#else
    write(6,*) '   flag for random perturbation test is not set'
#endif
    write(6,*) '   solar radiation frequency (iterations) = ',irad
#if (defined COUP_CSM)
    write(6,*) 'communication with the flux coupler'
    if (csm_doflxave) then
       write(6,*)'    data will be sent to the flux coupler ', &
            'only when an albedo calculation is performed '
       write(6,*)'     fluxes will be averaged on steps where ', &
            'communication with the flux coupler does not occur'
    else
       write(6,*)'    data will be sent and received to/from ', &
            'the flux coupler at every time step except for nstep=1'
    endif
#endif
#if (defined RTM)
    if (rtm_nsteps > 1) then
       write(6,*)'river runoff calculation performed only every ',rtm_nsteps,' nsteps'
    else
       write(6,*)'river runoff calculation performed every time step'
    endif
#endif
    if (nsrest == 1) then
       write(6,*) 'restart warning:'
       write(6,*) '   Namelist not checked for agreement with initial run.'
       write(6,*) '   Namelist should not differ except for ending time step and run type'
    end if
    if (nsrest == 3) then
       write(6,*) 'branch warning:'
       write(6,*) '   Namelist not checked for agreement with initial run.'
       write(6,*) '   Surface data set and reference date should not differ from initial run'
    end if
#if (defined COUP_CSM)
    write(6,*) '   last time step determined by flux coupler'
#endif
#if (defined PERGRO)
    write(6,*) '   perturbation limit = ',pertlim
#endif
    write(6,*) '   maxpatch_pft         = ',maxpatch_pft
    write(6,*) '   allocate_all_vegpfts = ',allocate_all_vegpfts
    write(6,*) '   nsegspc              = ',nsegspc

  end subroutine control_print

end module controlMod
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

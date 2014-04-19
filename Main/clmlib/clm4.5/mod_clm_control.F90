module mod_clm_control
  !
  ! Module which initializes run control variables. The following possible
  ! namelist variables are set default values and possibly read in on startup
  !
  ! Note: For definitions of namelist variables see
  !       ../../bld/namelist_files/namelist_definition.xml
  !       Display the file in a browser to see it neatly formatted in html.
  !
  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_mpmessage
  use mod_dynparam
  use mod_mppparam
  use mod_runparams
  use mod_clm_varctl , only : clmvarctl_init , set_clmvarctl , &
          nsrStartup , nsrContinue
  use mod_clm_varpar , only : maxpatch_pft , more_vertlayers
  use mod_clm_varctl , only : hostname , model_version=>version , &
          outnc_large_files , finidat , fsurdat , fatmlndfrc ,    &
          fpftdyn , fpftcon , nrevsn ,  create_crop_landunit ,    &
          allocate_all_vegpfts , co2_type , wrtdia , co2_ppmv ,   &
          pertlim , username , fsnowaging , fsnowoptics ,         &
          subgridflag , use_c13 , use_c14 , irrigate ,            &
          spinup_state , override_bgc_restart_mismatch_dump ,     &
          source , enable_megan_emission
  use mod_clm_varpar, only : numrad
  use mod_clm_varctl , only : NLFileName_in , ctitle , caseid , nsrest
  use mod_clm_varcon , only : secspday
  use mod_clm_canopyfluxes , only : perchroot , perchroot_alt
#if (defined LCH4) && (defined CN)
  use mod_clm_varctl , only : anoxia
#ifndef CENTURY_DECOMP
  use mod_clm_cndecompcascadebgc , only : anoxia_wtsat
#else
  use mod_clm_cndecompcascadecentury , only : anoxia_wtsat
#endif
#endif
  ! Lakes
  ! lake_use_old_fcrit_minz0 , lakepuddling , lake_puddle_thick
  ! and lake_no_ed are currently hardwired.
  use mod_clm_slakecon , only : deepmixing_depthcrit , deepmixing_mixfact , &
          lake_melt_icealb
!
  use mod_clm_surfacealbedo , only : albice
#ifdef CN
  use mod_clm_cnallocation , only : suplnitro , suplnNon
  use mod_clm_cnndynamics , only : nfix_timeconst
#endif
  use mod_clm_histfile , only : max_tapes, max_namlen, &
                 hist_empty_htapes, hist_dov2xy, &
                 hist_avgflag_pertape, hist_type1d_pertape, &
                 hist_nhtfrq, hist_ndens, &
                 hist_fincl1, hist_fincl2, hist_fincl3, &
                 hist_fincl4, hist_fincl5, hist_fincl6, &
                 hist_fexcl1, hist_fexcl2, hist_fexcl3, &
                 hist_fexcl4, hist_fexcl5, hist_fexcl6
#ifdef LCH4
  use mod_clm_histfile , only : hist_wrtch4diag
#endif
  use mod_clm_urban , only : urban_hac, urban_traffic

#if (defined CN) && (defined VERTSOILC)
  use mod_clm_cnsoillittverttransp , only: som_diffus , som_adv_flux , &
       cryoturb_diffusion_k , max_altdepth_cryoturbation , max_depth_cryoturb

  use mod_clm_cnverticalprofile , only: exponential_rooting_profile, &
       rootprof_exp , surfprof_exp , pftspecific_rootingprofile

#ifndef CENTURY_DECOMP
  use mod_clm_cndecompcascadebgc , only: decomp_depth_efolding , froz_q10
#else
  use mod_clm_cndecompcascadecentury , only: decomp_depth_efolding , froz_q10
#endif

#endif

#if (defined CN) && (defined NITRIF_DENITRIF)
  use mod_clm_cnnitrifdenitrif , only: no_frozen_nitrif_denitrif
#endif

#if (defined CN)
  !!! C14
  use mod_clm_cnc14decay , only: use_c14_bombspike , atm_c14_filename
#endif

  use mod_clm_surfacealbedo , only : albice
  use mod_clm_hydrology1 , only : Hydrology1_readnl
  use mod_clm_soilhydrology , only : SoilHydrology_readnl
  use mod_clm_megan, only : shr_megan_readnl , shr_megan_mechcomps_n
  implicit none

  private

  save

  public :: control_init  ! initial run control information
  public :: control_print ! print run control information

  character(len=  7) :: runtyp(4)                        ! run type

  contains
  !
  ! Initialize CLM run control information
  !
  subroutine control_init( )
    implicit none
    character(len=32)  :: starttype ! infodata start type
    integer(ik4) :: i , j , n            ! loop indices
    integer(ik4) :: ierr                 ! error code
    integer(ik4) :: ihost
    integer(ik4) :: unitn                ! unit for namelist file
    ! If want to override the startup type sent from driver
    character(len=32) :: subname = 'control_init'  ! subroutine name
    character(len=32) :: hostname = '?'
    character(len=32) :: user = '?'
    integer :: hostnm

    ! ----------------------------------------------------------------------
    ! Namelist Variables
    ! ----------------------------------------------------------------------

    ! CLM namelist settings

    namelist /clm_inparm / &
         fatmlndfrc, finidat, nrevsn

    ! Input datasets

    namelist /clm_inparm/  &
         fpftcon , fpftdyn , fsnowoptics , fsnowaging

    ! History, restart options

    namelist /clm_inparm/  &
         hist_empty_htapes, hist_dov2xy, &
         hist_avgflag_pertape, hist_type1d_pertape, &
         hist_nhtfrq,  hist_ndens, &
         hist_fincl1,  hist_fincl2, hist_fincl3, &
         hist_fincl4,  hist_fincl5, hist_fincl6, &
         hist_fexcl1,  hist_fexcl2, hist_fexcl3, &
         hist_fexcl4,  hist_fexcl5, hist_fexcl6, &
         outnc_large_files
#ifdef LCH4
    namelist /clm_inparm/ hist_wrtch4diag
#endif

    ! BGC info

#if (defined CN)
    namelist /clm_inparm/  &
         suplnitro
    namelist /clm_inparm/ &
         nfix_timeconst
    namelist /clm_inparm/ &
         spinup_state, override_bgc_restart_mismatch_dump
#endif

    namelist /clm_inparm / &
         co2_type

    namelist /clm_inparm / perchroot, perchroot_alt
#ifdef LCH4
    namelist /clm_inparm / anoxia, anoxia_wtsat
#endif
    namelist /clm_inparm / deepmixing_depthcrit, &
            deepmixing_mixfact,  lake_melt_icealb
    ! lake_melt_icealb is of dimension numrad

    ! Other options

    namelist /clm_inparm/  &
         wrtdia, pertlim, &
         create_crop_landunit, co2_ppmv, &
         albice, more_vertlayers, subgridflag, irrigate

    ! Urban options

    namelist /clm_inparm/  &
         urban_hac, urban_traffic

#if (defined CN) && (defined VERTSOILC)
    ! vertical soil mixing variables
    namelist /clm_inparm/  &
         som_diffus, som_adv_flux, cryoturb_diffusion_k, &
         max_altdepth_cryoturbation, max_depth_cryoturb

    ! depth inhibition of decomposition paramters
    namelist /clm_inparm/  &
         decomp_depth_efolding, froz_q10

    ! C and N input vertical profiles
    namelist /clm_inparm/  &
          exponential_rooting_profile, rootprof_exp, &
          surfprof_exp, pftspecific_rootingprofile
#endif

#if (defined CN) && (defined NITRIF_DENITRIF)
    namelist /clm_inparm / no_frozen_nitrif_denitrif
#endif

    namelist /clm_inparm / use_c13, use_c14

    namelist /clm_inparm / enable_megan_emission

#if (defined CN)
    !!! C14
    namelist /clm_inparm/  &
         use_c14_bombspike, atm_c14_filename
#endif

    ! ----------------------------------------------------------------------
    ! Default values
    ! ----------------------------------------------------------------------

    if (myid == italk) then
      write(stdout,*) 'Attempting to initialize run control settings .....'
    end if

    runtyp(:)               = 'missing'
    runtyp(nsrStartup  + 1) = 'initial'
    runtyp(nsrContinue + 1) = 'restart'

    if ( myid == iocpu ) then

      ! ----------------------------------------------------------------------
      ! Read namelist from standard input.
      ! ----------------------------------------------------------------------

      if ( len_trim(namelistfile) == 0  )then
        call fatal(__FILE__,__LINE__,subname//' error: nlfilename not set' )
      end if
      unitn = file_getunit( )
      write(stdout,*) 'Read in clm_inparm namelist from: ', trim(namelistfile)
      open(unitn,file=namelistfile,status='old',action='read',iostat=ierr)
      if (ierr == 0) then
        read(unitn, nml=clm_inparm, iostat=ierr ,err=100)
        if (ierr /= 0) then
          call fatal(__FILE__,__LINE__, &
              subname // ':: ERROR reading clm_inparm namelist')
        end if
      else
        write(stderr,*) 'Cannot open input namelist file ',trim(namelistfile)
        call fatal(__FILE__,__LINE__,'clm now stopping')
      end if
      call file_freeunit( unitn )

      ! ----------------------------------------------------------------------
      ! Consistency checks on input namelist.
      ! ----------------------------------------------------------------------

      if ( urban_traffic ) then
        write(stderr,*)'Urban traffic fluxes are not implemented currently'
        call fatal(__FILE__,__LINE__,'clm now stopping')
      end if

      ! History and restart files

      do i = 1 , max_tapes
        if (hist_nhtfrq(i) < 0) then
          hist_nhtfrq(i) = nint(-hist_nhtfrq(i)*secspday/(24.D0*dtsrf))
        end if
      end do

    end if   ! end if-block

#ifdef IBM
    hostname='ibm platform '
    user= 'Unknown'
#else
    ihost = hostnm(hostname)
    call getlog(user)
#endif
    if ( ifrest ) then
      call set_clmvarctl(domname, 'RegCM driven CLM4.5', nsrContinue, &
                         SVN_REV, hostname, user)
    else
      call set_clmvarctl(domname, 'RegCM driven CLM4.5', nsrStartup, &
                         SVN_REV, hostname, user)
    end if

    call clmvarctl_init

    ! ----------------------------------------------------------------------
    ! Read in other namelists for other modules
    ! ----------------------------------------------------------------------

    call Hydrology1_readnl(    namelistfile )
    call SoilHydrology_readnl( namelistfile )

    if ( shr_megan_mechcomps_n > 0 ) call shr_megan_readnl(namelistfile)
    ! ----------------------------------------------------------------------
    ! Broadcast all control information if appropriate
    ! ----------------------------------------------------------------------

    call control_spmd()

    ! Set input file path in RegCM world

    fpftcon = trim(inpglob)//pthsep//'CLM45'//pthsep// &
            'pftdata'//pthsep//fpftcon
    fsnowoptics = trim(inpglob)//pthsep//'CLM45'//pthsep// &
            'snicardata'//pthsep//fsnowoptics
    fsnowaging = trim(inpglob)//pthsep//'CLM45'//pthsep// &
            'snicardata'//pthsep//fsnowaging
    fsurdat = trim(dirglob)//pthsep//trim(domname)//'_CLM45_surface.nc'
    write(fatmlndfrc,'(a,i0.3,a)') &
      trim(dirglob)//pthsep//trim(domname)//'_DOMAIN',0,'.nc'

    if (myid == italk) then
      write(stdout,*) 'Successfully initialized run control settings'
      write(stdout,*)
    end if
    return
100 write (stderr,*) 'CLM : Error reading input namelist file !'
    call fatal(__FILE__,__LINE__,'clm now stopping')
  end subroutine control_init
  !
  ! Distribute namelist data all processors. All program i/o is
  ! funnelled through the master processor. Processor 0 either
  ! reads restart/history data from the disk and distributes
  ! it to all processors, or collects data from
  ! all processors and writes it to disk.
  !
  subroutine control_spmd()
    implicit none

    ! run control variables
    call bcast(caseid,len(caseid))
    call bcast(ctitle,len(ctitle))
    call bcast(model_version,len(model_version))
    call bcast(hostname,len(hostname))
    call bcast(username,len(username))
    call bcast(nsrest)

    ! initial file variables

    call bcast(nrevsn,len(nrevsn))
    call bcast(finidat,len(finidat))
    call bcast(fsurdat,len(fsurdat))
    call bcast(fatmlndfrc,len(fatmlndfrc))
    call bcast(fpftcon,len(fpftcon))
    call bcast(fpftdyn,len(fpftdyn))
    call bcast(fsnowoptics,len(fsnowoptics))
    call bcast(fsnowaging,len(fsnowaging))

    ! Irrigation

    call bcast(irrigate)

    ! Landunit generation

    call bcast(create_crop_landunit)
    call bcast(allocate_all_vegpfts)

    ! BGC

    call bcast(co2_type,len(co2_type))
#ifdef CN
    call bcast(suplnitro,len(suplnitro))
    call bcast(nfix_timeconst)
    call bcast(spinup_state)
    call bcast(spinup_state)
    call bcast(override_bgc_restart_mismatch_dump)
#endif

    ! isotopes
    call bcast(use_c13)
    call bcast(use_c14)

    call bcast(enable_megan_emission)

#if (defined CN) && (defined VERTSOILC)
    ! vertical soil mixing variables
    call bcast(som_diffus)
    call bcast(som_adv_flux)
    call bcast(cryoturb_diffusion_k)
    call bcast(max_altdepth_cryoturbation)
    call bcast(max_depth_cryoturb)

    ! depth inhibition of decomposition paramters
    call bcast(decomp_depth_efolding)
    call bcast(froz_q10)

    ! C and N input vertical profiles
    call bcast(exponential_rooting_profile)
    call bcast(rootprof_exp)
    call bcast(surfprof_exp)
    call bcast(pftspecific_rootingprofile)
#endif

#if (defined CN) && (defined NITRIF_DENITRIF)
    call bcast(no_frozen_nitrif_denitrif)
#endif

#if (defined CN)
    !!! C14
    call bcast(use_c14_bombspike)
    call bcast(atm_c14_filename,len(atm_c14_filename))
#endif

    call bcast(perchroot)
    call bcast(perchroot_alt)
#ifdef LCH4
    call bcast(anoxia)
    call bcast(anoxia_wtsat)
#endif
! Lakes
    call bcast(deepmixing_depthcrit)
    call bcast(deepmixing_mixfact)
    call bcast(lake_melt_icealb)

    ! physics variables

    call bcast(urban_hac,len(urban_hac))
    call bcast(urban_traffic)
    call bcast(subgridflag)
    call bcast(wrtdia)
    call bcast(co2_ppmv)
    call bcast(albice)
    call bcast(more_vertlayers)

    ! history file variables

    call bcast(outnc_large_files)
    call bcast(hist_empty_htapes)
    call bcast(hist_dov2xy)
    call bcast(hist_nhtfrq)
    call bcast(hist_ndens)
    call bcast(hist_avgflag_pertape,max_namlen)
    call bcast(hist_type1d_pertape,max_namlen)
#ifdef LCH4
    call bcast(hist_wrtch4diag)
#endif
    call bcast(hist_fexcl1,max_namlen)
    call bcast(hist_fexcl2,max_namlen)
    call bcast(hist_fexcl3,max_namlen)
    call bcast(hist_fexcl4,max_namlen)
    call bcast(hist_fexcl5,max_namlen)
    call bcast(hist_fexcl6,max_namlen)
    call bcast(hist_fincl1,max_namlen+2)
    call bcast(hist_fincl2,max_namlen+2)
    call bcast(hist_fincl3,max_namlen+2)
    call bcast(hist_fincl4,max_namlen+2)
    call bcast(hist_fincl5,max_namlen+2)
    call bcast(hist_fincl6,max_namlen+2)

    ! error growth perturbation limit
    call bcast(pertlim)

  end subroutine control_spmd
  !
  ! Write out the clm namelist run control variables
  !
  subroutine control_print ()
    implicit none
    integer(ik4) :: i  !loop index
    character(len=32) :: subname = 'control_print'  ! subroutine name

    if ( myid /= italk ) return

    write(stdout,*) 'define run:'
    write(stdout,*) '   source                = ',trim(source)
    write(stdout,*) '   model_version         = ',trim(model_version)
    write(stdout,*) '   run type              = ',runtyp(nsrest+1)
    write(stdout,*) '   case title            = ',trim(ctitle)
    write(stdout,*) '   username              = ',trim(username)
    write(stdout,*) '   hostname              = ',trim(hostname)
    write(stdout,*) 'input data files:'
    write(stdout,*) '   PFT physiology = ',trim(fpftcon)
    if (fsurdat == ' ') then
      write(stdout,*) '   fsurdat, surface dataset not set'
    else
      write(stdout,*) '   surface data   = ',trim(fsurdat)
    end if
    if (fatmlndfrc == ' ') then
      write(stdout,*) '   fatmlndfrc not set, setting frac/mask to 1'
    else
      write(stdout,*) '   land frac data = ',trim(fatmlndfrc)
    end if
#ifdef CN
    if (suplnitro /= suplnNon)then
       write(stdout,*) &
               '   Supplemental Nitrogen mode is set to run over PFTs: ', &
              trim(suplnitro)
    end if

    if (nfix_timeconst /= 0.D0) then
      write(stdout,*) &
        '   nfix_timeconst, timescale for smoothing npp in N fixation term: ', &
        nfix_timeconst
    else
      write(stdout,*) &
        '   nfix_timeconst == zero, use standard N fixation scheme. '
    end if

    write(stdout,*) &
        '   spinup_state, (0 = normal mode; 1 = AD spinup)         : ', &
        spinup_state
    if ( spinup_state .eq. 0 ) then
      write(stdout,*) '   model is currently NOT in AD spinup mode.'
    else if ( spinup_state .eq. 1 ) then
      write(stdout,*) '   model is currently in AD spinup mode.'
    else
      call fatal(__FILE__,__LINE__, &
      subname//' error: spinup_state can only have integer value of 0 or 1' )
    end if

    write(stdout,*) &
      '   override_bgc_restart_mismatch_dump                     : ', &
      override_bgc_restart_mismatch_dump
#endif

#if (defined CN) && (defined VERTSOILC)
    write(stdout, *) &
      '   som_adv_flux, the advection term in soil mixing (m/s) : ', &
      som_adv_flux
    write(stdout, *) &
      '   som_diffus, the diffusion term in soil mixing (m/s^2) : ', &
      som_diffus
    write(stdout, *) &
      '   cryoturb_diffusion_k  (m/s^2)                         : ', &
      cryoturb_diffusion_k
    write(stdout, *) &
      '   max_altdepth_cryoturbation (m)                        : ', &
      max_altdepth_cryoturbation
    write(stdout, *) &
      '   max_depth_cryoturb (m)                                : ', &
      max_depth_cryoturb
    write(stdout, *) &
      '   decomp_depth_efolding                                 : ', &
      decomp_depth_efolding
    write(stdout, *) &
      '   froz_q10                                              : ', &
      froz_q10
    write(stdout, *) &
      '   exponential_rooting_profile                           : ', &
      exponential_rooting_profile
    write(stdout, *) &
      '   rootprof_exp                                          : ', &
      rootprof_exp
    write(stdout, *) &
      '   surfprof_exp                                          : ', &
      surfprof_exp
    write(stdout, *) &
      '   pftspecific_rootingprofile                            : ', &
      pftspecific_rootingprofile
#endif

#if (defined CN) && (defined NITRIF_DENITRIF)
    write(stdout, *) &
      '   no_frozen_nitrif_denitrif                             : ', &
      no_frozen_nitrif_denitrif
#endif

    write(stdout, *) &
      '   enable_megan_emission                                 : ', &
      enable_megan_emission

#if (defined CN)
    write(stdout, *) &
      '  use_c13                                                : ', &
      use_c13
    write(stdout, *) &
      '  use_c14                                                : ', &
      use_c14
    !!! C14
    write(stdout, *) &
      '  use_c14_bombspike                                      : ', &
      use_c14_bombspike
    write(stdout, *) &
      '  atm_c14_filename                                       : ', &
      atm_c14_filename
#endif

    if (fsnowoptics == ' ') then
      write(stdout,*) '   snow optical properties file NOT set'
    else
      write(stdout,*) '   snow optical properties file = ',trim(fsnowoptics)
    end if
    if (fsnowaging == ' ') then
      write(stdout,*) '   snow aging parameters file NOT set'
    else
      write(stdout,*) '   snow aging parameters file = ',trim(fsnowaging)
    end if

    if (nsrest == nsrStartup ) &
      write(stdout,*) '   initial data is from RegCM atm model'
    if (nsrest /= nsrStartup) &
      write(stdout,*) '   restart data   = ',trim(nrevsn)
    write(stdout,*) '   atmospheric forcing data is from RegCM atm model'
    if ( outnc_large_files ) then
       write(stdout,*)'Large file support for output files is ON'
    end if
    write(stdout,*) 'model physics parameters:'
#if (defined PERGRO)
    write(stdout,*) '   flag for random perturbation test is set'
#else
    write(stdout,*) '   flag for random perturbation test is not set'
#endif
    write(stdout,*) '   CO2 volume mixing ratio   (umol/mol)   = ', co2_ppmv
    write(stdout,*) '   land-ice albedos      (unitless 0-1)   = ', albice
    write(stdout,*) &
         '   urban air conditioning/heating and wasteheat   = ', urban_hac
    write(stdout,*) '   urban traffic flux   = ', urban_traffic
    write(stdout,*) '   more vertical layers = ', more_vertlayers
    if (nsrest == nsrContinue) then
      write(stdout,*) 'restart warning:'
      write(stdout,*) '   Namelist not checked for agreement with initial run.'
      write(stdout,*) '   Namelist should not differ except for &
                      &ending time step and run type'
    end if
    if ( pertlim /= 0.0D0 ) &
      write(stdout,*) '   perturbation limit   = ',pertlim
    write(stdout,*) '   maxpatch_pft         = ',maxpatch_pft
    write(stdout,*) '   allocate_all_vegpfts = ',allocate_all_vegpfts
! New fields
    write(stdout,*) ' perchroot (plant water stress based on unfrozen &
                   &layers only) = ',perchroot
    write(stdout,*) ' perchroot (plant water stress based on &
                   &time-integrated active layer only) = ',perchroot
#ifdef LCH4
    write(stdout,*) &
      ' anoxia (applied to soil decomposition)             = ',anoxia
    write(stdout,*) &
      ' anoxia_wtsat (weight anoxia by inundated fraction) = ',anoxia_wtsat
#endif
! Lakes
    write(stdout,*)
    write(stdout,*) 'Lake Model Namelists:'
    write(stdout,*) &
      'Increased mixing relative to Hostetler wind-driven eddy expression ',&
      'will be used for deep lakes exceeding depth ', deepmixing_depthcrit,&
      ' by a factor of ', deepmixing_mixfact, '.'
    write(stdout,*) &
      'Albedo over melting lakes will approach values (visible, NIR):', &
      lake_melt_icealb, &
      'as compared with 0.60, 0.40 for cold frozen lakes with no snow.'
  end subroutine control_print

end module mod_clm_control

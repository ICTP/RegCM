module mod_clm_cnrest

#if (defined CN)
  !
  ! Read/Write to/from CN info to CLM restart file.
  !
  use mod_realkinds
  use mod_dynparam
  use mod_mppparam
  use mod_mpmessage
  use mod_runparams , only : ktau , ntsrf
  use mod_clm_nchelper
  use mod_clm_type
  use mod_clm_decomp
  use mod_clm_surfrd , only : crop_prog
  use mod_clm_atmlnd, only : clm_a2l
  use mod_clm_varpar, only : numrad , ndecomp_pools , nlevdecomp
  use mod_clm_varpar , only : nlevgrnd
  use mod_clm_varctl , only : override_bgc_restart_mismatch_dump
  use mod_clm_varctl , only : use_c13 , use_c14 , spinup_state
  use mod_clm_varcon , only : c13ratio , c14ratio , spval
  use mod_constants , only : pdbratio

  implicit none

  private

  save

  public :: CNrest

  contains
  !
  ! Read/write CN restart data
  !
  subroutine CNRest ( ncid, flag )
    implicit none
    type(clm_filetype)  :: ncid   ! netcdf id
    character(len=*) , intent(in) :: flag   !'read' or 'write'
    ! typical del13C for C3 photosynthesis (permil, relative to PDB)
    real(rk8) :: c3_del13c
    ! typical del13C for C4 photosynthesis (permil, relative to PDB)
    real(rk8) :: c4_del13c
    ! isotope ratio (13c/12c) for C3 photosynthesis
    real(rk8) :: c3_r1
    ! isotope ratio (13c/12c) for C4 photosynthesis
    real(rk8) :: c4_r1
    ! isotope ratio (13c/[12c+13c]) for C3 photosynthesis
    real(rk8) :: c3_r2
    ! isotope ratio (13c/[12c+13c]) for C4 photosynthesis
    real(rk8) :: c4_r2
!   real(rk8) , pointer :: rc13_annsum_npp(:)
!   real(rk8) , pointer :: rc13_cannsum_npp(:)
    type(pft_cstate_type) , pointer :: pcisos
    type(pft_cstate_type) , pointer :: pcbulks
    integer(ik4) :: c , j , k , i ! indices
    integer(ik4) :: begp , endp ! per-proc beginning and ending pft indices
    integer(ik4) :: begc , endc ! per-proc beginning and ending column indices
    integer(ik4) :: begl , endl ! per-proc beginning and ending landunit indices
    integer(ik4) :: begg , endg   ! per-proc gridcell ending gridcell indices
    real(rk8):: m            ! multiplier for the exit_spinup code
    logical :: lstop         ! determine if variable is on initial file
    character(len=128) :: varname         ! temporary
    type(gridcell_type), pointer :: gptr  ! pointer to gridcell derived subtype
    type(landunit_type), pointer :: lptr  ! pointer to landunit derived subtype
    type(column_type)  , pointer :: cptr  ! pointer to column derived subtype
    type(pft_type)     , pointer :: pptr  ! pointer to pft derived subtype
#if (defined CNDV)
    integer(ik4) , pointer :: iptemp(:) ! pointer to memory to be allocated
    integer(ik4) :: p , ier! indices
#endif
    !temporary arrays for slicing larger arrays
    real(rk8) , pointer :: ptr2d(:,:)
    ! spinup state as read from restart file, for determining whether
    ! to enter or exit spinup mode.
    integer(ik4) :: restart_file_spinup_state
    logical :: exit_spinup = .false.
    logical :: enter_spinup = .false.
    ! flags for comparing the model and restart decomposition cascades
    integer(ik4) :: decomp_cascade_state, restart_file_decomp_cascade_state

    if ( use_c13 ) then
      pcisos => clm3%g%l%c%p%pc13s
      pcbulks => clm3%g%l%c%p%pcs
    end if
    ! Set pointers into derived type

    gptr => clm3%g
    lptr => clm3%g%l
    cptr => clm3%g%l%c
    pptr => clm3%g%l%c%p

    ! Determine necessary subgrid bounds

    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

    if ( use_c13 ) then
      c3_del13c = -28.D0
      c4_del13c = -13.D0
      c3_r1 = pdbratio + ((c3_del13c*pdbratio)/1000.D0)
      c3_r2 = c3_r1/(1.D0 + c3_r1)
      c4_r1 = pdbratio + ((c4_del13c*pdbratio)/1000.D0)
      c4_r2 = c4_r1/(1.D0 + c4_r1)
    end if

    !--------------------------------
    ! pft ecophysiological variables
    !--------------------------------

    ! dormant_flag
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'dormant_flag',(/'pft'/), &
            long_name='dormancy flag',units='unitless' )
    else if ( flag == 'read' ) then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'dormant_flag') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'dormant_flag',pptr%pepv%dormant_flag,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'dormant_flag',pptr%pepv%dormant_flag,gcomm_pft)
    end if

    ! days_active
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'days_active',(/'pft'/), &
            long_name='number of days since last dormancy',units='days' )
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'days_active') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'days_active',pptr%pepv%days_active,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'days_active',pptr%pepv%days_active,gcomm_pft)
    end if

    ! onset_flag
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'onset_flag',(/'pft'/), &
            long_name='flag if critical growing degree-day sum is exceeded', &
              units='unitless' )
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'onset_flag') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'onset_flag',pptr%pepv%onset_flag,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'onset_flag',pptr%pepv%onset_flag,gcomm_pft)
    end if

    ! onset_counter
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'onset_counter',(/'pft'/), &
            long_name='onset days counter',units='sec' )
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'onset_counter') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'onset_counter',pptr%pepv%onset_counter,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'onset_counter',pptr%pepv%onset_counter,gcomm_pft)
    end if

    ! onset_gddflag
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'onset_gddflag',(/'pft'/), &
            long_name='onset flag for growing degree day sum',units='' )
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'onset_gddflag') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'onset_gddflag',pptr%pepv%onset_gddflag,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'onset_gddflag',pptr%pepv%onset_gddflag,gcomm_pft)
    end if

    ! onset_fdd
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'onset_fdd',(/'pft'/), &
            long_name='onset freezing degree days counter',units='days' )
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'onset_fdd') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'onset_fdd',pptr%pepv%onset_fdd,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'onset_fdd',pptr%pepv%onset_fdd,gcomm_pft)
    end if

    ! onset_gdd
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'onset_gdd',(/'pft'/), &
            long_name='onset growing degree days',units='days' )
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'onset_gdd') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
       call clm_readvar(ncid,'onset_gdd',pptr%pepv%onset_gdd,gcomm_pft)
     end if
   else if (flag == 'write') then
     call clm_writevar(ncid,'onset_gdd',pptr%pepv%onset_gdd,gcomm_pft)
   end if

    ! onset_swi
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'onset_swi',(/'pft'/), &
            long_name='onset soil water index',units='days' )
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'onset_swi') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'onset_swi',pptr%pepv%onset_swi,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'onset_swi',pptr%pepv%onset_swi,gcomm_pft)
    end if

    ! offset_flag
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'offset_flag',(/'pft'/), &
            long_name='offset flag',units='unitless' )
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'offset_flag') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'offset_flag',pptr%pepv%offset_flag,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'offset_flag',pptr%pepv%offset_flag,gcomm_pft)
    end if

    ! offset_counter
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'offset_counter',(/'pft'/), &
            long_name='offset days counter',units='sec' )
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'offset_counter') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'offset_counter', &
                pptr%pepv%offset_counter,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'offset_counter', &
              pptr%pepv%offset_counter,gcomm_pft)
    end if

    ! offset_fdd
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'offset_fdd',(/'pft'/), &
            long_name='offset freezing degree days counter',units='days' )
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'offset_fdd') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'offset_fdd',pptr%pepv%offset_fdd,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'offset_fdd',pptr%pepv%offset_fdd,gcomm_pft)
    end if

    ! offset_swi
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'offset_swi',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'offset_swi') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'offset_swi',pptr%pepv%offset_swi,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'offset_swi',pptr%pepv%offset_swi,gcomm_pft)
    end if

#if (defined CROP)
    ! fert_counter
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'fert_counter',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'fert_counter') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'fert_counter',pptr%pepv%fert_counter,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'fert_counter',pptr%pepv%fert_counter,gcomm_pft)
    end if

    ! fert
    if ( crop_prog ) then
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'fert',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( ktau /= 0 .and. .not. clm_check_var(ncid,'fert') ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        else
          call clm_readvar(ncid,'fert',pptr%pnf%fert,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'fert',pptr%pnf%fert,gcomm_pft)
      end if
    end if
#endif

    ! lgsf
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'lgsf',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'lgsf') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'lgsf',pptr%pepv%lgsf,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'lgsf',pptr%pepv%lgsf,gcomm_pft)
    end if

    ! bglfr
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'bglfr',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'bglfr') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'bglfr',pptr%pepv%bglfr,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'bglfr',pptr%pepv%bglfr,gcomm_pft)
    end if

    ! bgtr
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'bgtr',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'bgtr') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'bgtr',pptr%pepv%bgtr,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'bgtr',pptr%pepv%bgtr,gcomm_pft)
    end if

    ! dayl
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'dayl',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'dayl') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'dayl',pptr%pepv%dayl,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'dayl',pptr%pepv%dayl,gcomm_pft)
    end if

    ! prev_dayl
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'prev_dayl',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'prev_dayl') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'prev_dayl',pptr%pepv%prev_dayl,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'prev_dayl',pptr%pepv%prev_dayl,gcomm_pft)
    end if

    ! annavg_t2m
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'annavg_t2m',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'annavg_t2m') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'annavg_t2m',pptr%pepv%annavg_t2m,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'annavg_t2m',pptr%pepv%annavg_t2m,gcomm_pft)
    end if

    ! tempavg_t2m
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'tempavg_t2m',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'tempavg_t2m') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'tempavg_t2m',pptr%pepv%tempavg_t2m,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'tempavg_t2m',pptr%pepv%tempavg_t2m,gcomm_pft)
    end if

    ! gpp
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'gpp_pepv',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'gpp_pepv') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'gpp_pepv',pptr%pepv%gpp,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'gpp_pepv',pptr%pepv%gpp,gcomm_pft)
    end if

    ! availc
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'availc',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'availc') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'availc',pptr%pepv%availc,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'availc',pptr%pepv%availc,gcomm_pft)
    end if

    ! xsmrpool_recover
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'xsmrpool_recover',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'xsmrpool_recover') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'xsmrpool_recover',pptr%pepv%xsmrpool_recover, &
                gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'xsmrpool_recover',pptr%pepv%xsmrpool_recover, &
              gcomm_pft)
    end if

    if ( use_c13 ) then
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'xsmrpool_c13ratio',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( ktau /= 0 .and. &
             .not. clm_check_var(ncid,'xsmrpool_c13ratio') ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        else
          call clm_readvar(ncid,'xsmrpool_c13ratio', &
                  pptr%pepv%xsmrpool_c13ratio, gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'xsmrpool_c13ratio', &
                pptr%pepv%xsmrpool_c13ratio, gcomm_pft)
      end if
    end if

    ! alloc_pnow
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'alloc_pnow',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'alloc_pnow') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'alloc_pnow',pptr%pepv%alloc_pnow,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'alloc_pnow',pptr%pepv%alloc_pnow,gcomm_pft)
    end if

    ! c_allometry
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'c_allometry',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'c_allometry') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'c_allometry',pptr%pepv%c_allometry,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'c_allometry',pptr%pepv%c_allometry,gcomm_pft)
    end if

    ! n_allometry
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'n_allometry',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'n_allometry') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'n_allometry',pptr%pepv%n_allometry,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'n_allometry',pptr%pepv%n_allometry,gcomm_pft)
    end if

    ! plant_ndemand
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'plant_ndemand',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'plant_ndemand') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'plant_ndemand',pptr%pepv%plant_ndemand, &
                gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'plant_ndemand',pptr%pepv%plant_ndemand, &
              gcomm_pft)
    end if

    ! tempsum_potential_gpp
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'tempsum_potential_gpp',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. &
           .not. clm_check_var(ncid,'tempsum_potential_gpp') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'tempsum_potential_gpp', &
                pptr%pepv%tempsum_potential_gpp,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'tempsum_potential_gpp', &
              pptr%pepv%tempsum_potential_gpp,gcomm_pft)
    end if

    !annsum_potential_gpp
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'annsum_potential_gpp',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. &
           .not. clm_check_var(ncid,'annsum_potential_gpp') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'annsum_potential_gpp', &
                pptr%pepv%annsum_potential_gpp,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'annsum_potential_gpp', &
              pptr%pepv%annsum_potential_gpp,gcomm_pft)
    end if

    ! tempmax_retransn
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'tempmax_retransn',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. &
           .not. clm_check_var(ncid,'tempmax_retransn') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'tempmax_retransn', &
                pptr%pepv%tempmax_retransn,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'tempmax_retransn', &
              pptr%pepv%tempmax_retransn,gcomm_pft)
    end if

    ! annmax_retransn
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'annmax_retransn',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. &
           .not. clm_check_var(ncid,'annmax_retransn') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'annmax_retransn', &
                pptr%pepv%annmax_retransn,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'annmax_retransn', &
              pptr%pepv%annmax_retransn,gcomm_pft)
    end if

    ! avail_retransn
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'avail_retransn',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. &
           .not. clm_check_var(ncid,'avail_retransn') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'avail_retransn', &
                pptr%pepv%avail_retransn,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'avail_retransn', &
              pptr%pepv%avail_retransn,gcomm_pft)
    end if

    ! plant_nalloc
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'plant_nalloc',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and.  .not. clm_check_var(ncid,'plant_nalloc') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'plant_nalloc', &
                pptr%pepv%plant_nalloc,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'plant_nalloc', &
              pptr%pepv%plant_nalloc,gcomm_pft)
    end if

    ! plant_calloc
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'plant_calloc',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and.  .not. clm_check_var(ncid,'plant_calloc') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'plant_calloc', &
                pptr%pepv%plant_calloc,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'plant_calloc', &
              pptr%pepv%plant_calloc,gcomm_pft)
    end if

    ! excess_cflux
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'excess_cflux',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and.  .not. clm_check_var(ncid,'excess_cflux') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'excess_cflux', &
                pptr%pepv%excess_cflux,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'excess_cflux', &
              pptr%pepv%excess_cflux,gcomm_pft)
    end if

    ! downreg
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'downreg',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
            if ( ktau /= 0 .and.  .not. clm_check_var(ncid,'downreg') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'downreg',pptr%pepv%downreg,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'downreg',pptr%pepv%downreg,gcomm_pft)
    end if

    ! prev_leafc_to_litter
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'prev_leafc_to_litter',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. &
           .not. clm_check_var(ncid,'prev_leafc_to_litter') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'prev_leafc_to_litter', &
                pptr%pepv%prev_leafc_to_litter,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'prev_leafc_to_litter', &
              pptr%pepv%prev_leafc_to_litter,gcomm_pft)
    end if

    ! prev_frootc_to_litter
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'prev_frootc_to_litter',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. &
           .not. clm_check_var(ncid,'prev_frootc_to_litter') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'prev_frootc_to_litter', &
                pptr%pepv%prev_frootc_to_litter,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'prev_frootc_to_litter', &
              pptr%pepv%prev_frootc_to_litter,gcomm_pft)
    end if

    ! tempsum_npp
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'tempsum_npp',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'tempsum_npp') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'tempsum_npp', &
                pptr%pepv%tempsum_npp,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'tempsum_npp', &
              pptr%pepv%tempsum_npp,gcomm_pft)
    end if

    ! annsum_npp
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'annsum_npp',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'annsum_npp') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'annsum_npp', &
                pptr%pepv%annsum_npp,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'annsum_npp', &
              pptr%pepv%annsum_npp,gcomm_pft)
    end if

    if ( use_c13 ) then
      ! rc13_canair
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'rc13_canair',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( ktau /= 0 .and. .not. clm_check_var(ncid,'rc13_canair') ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        else
          call clm_readvar(ncid,'rc13_canair', &
                  pptr%pepv%rc13_canair,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'rc13_canair', &
                pptr%pepv%rc13_canair,gcomm_pft)
      end if

      ! rc13_psnsun
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'rc13_psnsun',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( ktau /= 0 .and. .not. clm_check_var(ncid,'rc13_psnsun') ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        else
          call clm_readvar(ncid,'rc13_psnsun', &
                  pptr%pepv%rc13_psnsun,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'rc13_psnsun', &
                pptr%pepv%rc13_psnsun,gcomm_pft)
      end if

      ! rc13_psnsha
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'rc13_psnsha',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( ktau /= 0 .and. .not. clm_check_var(ncid,'rc13_psnsha') ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        else
          call clm_readvar(ncid,'rc13_psnsha', &
                  pptr%pepv%rc13_psnsha,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'rc13_psnsha', &
                pptr%pepv%rc13_psnsha,gcomm_pft)
      end if
    end if

#if (defined CROP)
    ! grain_flag
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'grain_flag',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'grain_flag') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'grain_flag', &
                pptr%pepv%grain_flag,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'grain_flag', &
              pptr%pepv%grain_flag,gcomm_pft)
    end if
#endif

    !--------------------------------
    ! pft carbon state variables
    !--------------------------------

    ! leafc
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'leafc',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'leafc') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'leafc',pptr%pcs%leafc,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'leafc',pptr%pcs%leafc,gcomm_pft)
    end if

    ! leafc_storage
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'leafc_storage',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'leafc_storage') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'leafc_storage', &
                pptr%pcs%leafc_storage,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'leafc_storage', &
              pptr%pcs%leafc_storage,gcomm_pft)
    end if

    ! leafc_xfer
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'leafc_xfer',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'leafc_xfer') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'leafc_xfer',pptr%pcs%leafc_xfer,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'leafc_xfer',pptr%pcs%leafc_xfer,gcomm_pft)
    end if

    ! frootc
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'frootc',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'frootc') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'frootc',pptr%pcs%frootc,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'frootc',pptr%pcs%frootc,gcomm_pft)
    end if

    ! frootc_storage
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'frootc_storage',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'frootc_storage') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'frootc_storage', &
                pptr%pcs%frootc_storage,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'frootc_storage', &
              pptr%pcs%frootc_storage,gcomm_pft)
    end if

    !frootc_xfer
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'frootc_xfer',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'frootc_xfer') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'frootc_xfer',pptr%pcs%frootc_xfer,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'frootc_xfer',pptr%pcs%frootc_xfer,gcomm_pft)
    end if

    ! livestemc
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'livestemc',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'livestemc') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'livestemc',pptr%pcs%livestemc,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'livestemc',pptr%pcs%livestemc,gcomm_pft)
    end if

    ! livestemc_storage
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'livestemc_storage',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'livestemc_storage') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'livestemc_storage', &
                pptr%pcs%livestemc_storage,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'livestemc_storage', &
              pptr%pcs%livestemc_storage,gcomm_pft)
    end if

    ! livestemc_xfer
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'livestemc_xfer',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'livestemc_xfer') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'livestemc_xfer', &
                pptr%pcs%livestemc_xfer,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'livestemc_xfer', &
              pptr%pcs%livestemc_xfer,gcomm_pft)
    end if

    ! deadstemc
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'deadstemc',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'deadstemc') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'deadstemc',pptr%pcs%deadstemc,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'deadstemc',pptr%pcs%deadstemc,gcomm_pft)
    end if

    ! deadstemc_storage
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'deadstemc_storage',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'deadstemc_storage') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'deadstemc_storage', &
                pptr%pcs%deadstemc_storage,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'deadstemc_storage', &
              pptr%pcs%deadstemc_storage,gcomm_pft)
    end if

    ! deadstemc_xfer
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'deadstemc_xfer',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'deadstemc_xfer') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'deadstemc_xfer', &
                pptr%pcs%deadstemc_xfer,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'deadstemc_xfer', &
              pptr%pcs%deadstemc_xfer,gcomm_pft)
    end if

    ! livecrootc
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'livecrootc',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'livecrootc') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'livecrootc',pptr%pcs%livecrootc,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'livecrootc',pptr%pcs%livecrootc,gcomm_pft)
    end if

    ! livecrootc_storage
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'livecrootc_storage',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. &
           .not. clm_check_var(ncid,'livecrootc_storage') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'livecrootc_storage', &
                pptr%pcs%livecrootc_storage,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'livecrootc_storage', &
              pptr%pcs%livecrootc_storage,gcomm_pft)
    end if

    ! livecrootc_xfer
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'livecrootc_xfer',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'livecrootc_xfer') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'livecrootc_xfer', &
                pptr%pcs%livecrootc_xfer,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'livecrootc_xfer', &
              pptr%pcs%livecrootc_xfer,gcomm_pft)
    end if

    ! deadcrootc
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'deadcrootc',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'deadcrootc') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'deadcrootc',pptr%pcs%deadcrootc,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'deadcrootc',pptr%pcs%deadcrootc,gcomm_pft)
    end if

    ! deadcrootc_storage
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'deadcrootc_storage',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. &
           .not. clm_check_var(ncid,'deadcrootc_storage') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'deadcrootc_storage', &
                pptr%pcs%deadcrootc_storage,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'deadcrootc_storage', &
              pptr%pcs%deadcrootc_storage,gcomm_pft)
    end if

    ! deadcrootc_xfer
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'deadcrootc_xfer',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'deadcrootc_xfer') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'deadcrootc_xfer', &
                pptr%pcs%deadcrootc_xfer,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'deadcrootc_xfer', &
              pptr%pcs%deadcrootc_xfer,gcomm_pft)
    end if

    ! gresp_storage
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'gresp_storage',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'gresp_storage') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'gresp_storage', &
                pptr%pcs%gresp_storage,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'gresp_storage', &
              pptr%pcs%gresp_storage,gcomm_pft)
    end if

    ! gresp_xfer
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'gresp_xfer',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'gresp_xfer') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'gresp_xfer',pptr%pcs%gresp_xfer,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'gresp_xfer',pptr%pcs%gresp_xfer,gcomm_pft)
    end if

    ! cpool
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'cpool',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'cpool') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'cpool',pptr%pcs%cpool,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'cpool',pptr%pcs%cpool,gcomm_pft)
    end if

    ! xsmrpool
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'xsmrpool',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'xsmrpool') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'xsmrpool',pptr%pcs%xsmrpool,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'xsmrpool',pptr%pcs%xsmrpool,gcomm_pft)
    end if

    ! pft_ctrunc
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'pft_ctrunc',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'pft_ctrunc') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'pft_ctrunc',pptr%pcs%pft_ctrunc,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'pft_ctrunc',pptr%pcs%pft_ctrunc,gcomm_pft)
    end if

    ! totvegc
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'totvegc',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'totvegc') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'totvegc',pptr%pcs%totvegc,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'totvegc',pptr%pcs%totvegc,gcomm_pft)
    end if

    if ( use_c13 ) then
      !--------------------------------
      ! C13 pft carbon state variables
      !--------------------------------

      ! leafc
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'leafc_13',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'leafc_13') ) then
          if ( ktau == 0 ) then
            write(stdout,*) 'initializing C13 leafc with atmospheric c13 value'
            do i = begp , endp
              if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1.D0) then
                pcisos%leafc(i) = pcbulks%leafc(i) * c3_r2
              else
                pcisos%leafc(i) = pcbulks%leafc(i) * c4_r2
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'leafc_13',pptr%pc13s%leafc,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'leafc_13',pptr%pc13s%leafc,gcomm_pft)
      end if

      ! leafc_storage
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'leafc_storage_13',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'leafc_storage_13') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
              'initializing C13 leafc_storage_13 with atmospheric c13 value'
            do i = begp , endp
              if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1.D0) then
                pcisos%leafc_storage(i) = pcbulks%leafc_storage(i) * c3_r2
              else
                pcisos%leafc_storage(i) = pcbulks%leafc_storage(i) * c4_r2
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'leafc_storage_13', &
                  pptr%pc13s%leafc_storage,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'leafc_storage_13', &
                pptr%pc13s%leafc_storage,gcomm_pft)
      end if

      ! leafc_xfer
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'leafc_xfer_13',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'leafc_xfer_13') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
              'initializing C13 leafc_xfer_13 with atmospheric c13 value'
            do i = begp , endp
              if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1.D0) then
                pcisos%leafc_xfer(i) = pcbulks%leafc_xfer(i) * c3_r2
              else
                pcisos%leafc_xfer(i) = pcbulks%leafc_xfer(i) * c4_r2
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'leafc_xfer_13', &
                  pptr%pc13s%leafc_xfer,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'leafc_xfer_13', &
                pptr%pc13s%leafc_xfer,gcomm_pft)
      end if

      ! frootc
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'frootc_13',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'frootc_13') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
              'initializing C13 frootc_13 with atmospheric c13 value'
            do i = begp , endp
              if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1.D0) then
                pcisos%frootc(i) = pcbulks%frootc(i) * c3_r2
              else
                pcisos%frootc(i) = pcbulks%frootc(i) * c4_r2
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'frootc_13',pptr%pc13s%frootc,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'frootc_13',pptr%pc13s%frootc,gcomm_pft)
      end if

      ! frootc_storage
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'frootc_storage_13',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'frootc_storage_13') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
              'initializing C13 frootc_storage_13 with atmospheric c13 value'
            do i = begp , endp
              if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1.D0) then
                pcisos%frootc_storage(i) = pcbulks%frootc_storage(i) * c3_r2
              else
                pcisos%frootc_storage(i) = pcbulks%frootc_storage(i) * c4_r2
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'frootc_storage_13', &
                  pptr%pc13s%frootc_storage,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'frootc_storage_13', &
                pptr%pc13s%frootc_storage,gcomm_pft)
      end if

      !frootc_xfer
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'frootc_xfer_13',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'frootc_xfer_13') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
              'initializing C13 frootc_xfer_13 with atmospheric c13 value'
            do i = begp , endp
              if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1.D0) then
                pcisos%frootc_xfer(i) = pcbulks%frootc_xfer(i) * c3_r2
              else
                pcisos%frootc_xfer(i) = pcbulks%frootc_xfer(i) * c4_r2
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'frootc_xfer_13', &
                  pptr%pc13s%frootc_xfer,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'frootc_xfer_13', &
                pptr%pc13s%frootc_xfer,gcomm_pft)
      end if

      ! livestemc
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'livestemc_13',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'livestemc_13') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
              'initializing C13 livestemc_13 with atmospheric c13 value'
            do i = begp , endp
              if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1.D0) then
                pcisos%livestemc(i) = pcbulks%livestemc(i) * c3_r2
              else
                pcisos%livestemc(i) = pcbulks%livestemc(i) * c4_r2
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'livestemc_13', &
                  pptr%pc13s%livestemc,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'livestemc_13', &
                pptr%pc13s%livestemc,gcomm_pft)
      end if

      ! livestemc_storage
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'livestemc_storage_13',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'livestemc_storage_13') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
              'initializing C13 livestemc_storage_13 with atmospheric c13 value'
            do i = begp , endp
              if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1.D0) then
                pcisos%livestemc_storage(i) = &
                        pcbulks%livestemc_storage(i) * c3_r2
              else
                pcisos%livestemc_storage(i) = &
                        pcbulks%livestemc_storage(i) * c4_r2
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'livestemc_storage_13', &
                  pptr%pc13s%livestemc_storage,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'livestemc_storage_13', &
                pptr%pc13s%livestemc_storage,gcomm_pft)
      end if

      ! livestemc_xfer
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'livestemc_xfer_13',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'livestemc_xfer_13') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
              'initializing C13 livestemc_xfer_13 with atmospheric c13 value'
            do i = begp , endp
              if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1.D0) then
                pcisos%livestemc_xfer(i) = pcbulks%livestemc_xfer(i) * c3_r2
              else
                pcisos%livestemc_xfer(i) = pcbulks%livestemc_xfer(i) * c4_r2
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'livestemc_xfer_13', &
                  pptr%pc13s%livestemc_xfer,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'livestemc_xfer_13', &
                pptr%pc13s%livestemc_xfer,gcomm_pft)
      end if

      ! deadstemc
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'deadstemc_13',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'deadstemc_13') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
              'initializing C13 deadstemc_13 with atmospheric c13 value'
            do i = begp , endp
              if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1.D0) then
                pcisos%deadstemc(i) = pcbulks%deadstemc(i) * c3_r2
              else
                pcisos%deadstemc(i) = pcbulks%deadstemc(i) * c4_r2
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'deadstemc_13', &
                  pptr%pc13s%deadstemc,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'deadstemc_13', &
                pptr%pc13s%deadstemc,gcomm_pft)
      end if

      ! deadstemc_storage
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'deadstemc_storage_13',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'deadstemc_storage_13') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
              'initializing C13 deadstemc_storage_13 with atmospheric c13 value'
            do i = begp , endp
              if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1.D0) then
                pcisos%deadstemc_storage(i) = &
                        pcbulks%deadstemc_storage(i) * c3_r2
              else
                pcisos%deadstemc_storage(i) = &
                        pcbulks%deadstemc_storage(i) * c4_r2
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'deadstemc_storage_13', &
                  pptr%pc13s%deadstemc_storage,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'deadstemc_storage_13', &
                pptr%pc13s%deadstemc_storage,gcomm_pft)
      end if

      ! deadstemc_xfer
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'deadstemc_xfer_13',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'deadstemc_xfer_13') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
              'initializing C13 deadstemc_xfer_13 with atmospheric c13 value'
            do i = begp , endp
              if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1.D0) then
                pcisos%deadstemc_xfer(i) = pcbulks%deadstemc_xfer(i) * c3_r2
              else
                pcisos%deadstemc_xfer(i) = pcbulks%deadstemc_xfer(i) * c4_r2
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'deadstemc_xfer_13', &
                  pptr%pc13s%deadstemc_xfer,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'deadstemc_xfer_13', &
                pptr%pc13s%deadstemc_xfer,gcomm_pft)
      end if

      ! livecrootc
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'livecrootc_13',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'livecrootc_13') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
              'initializing C13 livecrootc_13 with atmospheric c13 value'
            do i = begp , endp
              if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1.D0) then
                pcisos%livecrootc(i) = pcbulks%livecrootc(i) * c3_r2
              else
                pcisos%livecrootc(i) = pcbulks%livecrootc(i) * c4_r2
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'livecrootc_13', &
                  pptr%pc13s%livecrootc,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'livecrootc_13', &
                pptr%pc13s%livecrootc,gcomm_pft)
      end if

      ! livecrootc_storage
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'livecrootc_storage_13',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'livecrootc_storage_13') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
             'initializing C13 livecrootc_storage_13 with atmospheric c13 value'
            do i = begp , endp
              if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1.D0) then
                pcisos%livecrootc_storage(i) = &
                        pcbulks%livecrootc_storage(i) * c3_r2
              else
                pcisos%livecrootc_storage(i) = &
                        pcbulks%livecrootc_storage(i) * c4_r2
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'livecrootc_storage_13', &
                  pptr%pc13s%livecrootc_storage,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'livecrootc_storage_13', &
                pptr%pc13s%livecrootc_storage,gcomm_pft)
      end if

      ! livecrootc_xfer
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'livecrootc_xfer_13',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'livecrootc_xfer_13') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
             'initializing C13 livecrootc_xfer_13 with atmospheric c13 value'
            do i = begp , endp
              if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1.D0) then
                pcisos%livecrootc_xfer(i) = pcbulks%livecrootc_xfer(i) * c3_r2
              else
                pcisos%livecrootc_xfer(i) = pcbulks%livecrootc_xfer(i) * c4_r2
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'livecrootc_xfer_13', &
                  pptr%pc13s%livecrootc_xfer,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'livecrootc_xfer_13', &
                pptr%pc13s%livecrootc_xfer,gcomm_pft)
      end if

      ! deadcrootc
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'deadcrootc_13',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'deadcrootc_13') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
             'initializing C13 deadcrootc_13 with atmospheric c13 value'
            do i = begp , endp
              if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1.D0) then
                pcisos%deadcrootc(i) = pcbulks%deadcrootc(i) * c3_r2
              else
                pcisos%deadcrootc(i) = pcbulks%deadcrootc(i) * c4_r2
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'deadcrootc_13', &
                  pptr%pc13s%deadcrootc,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'deadcrootc_13', &
                pptr%pc13s%deadcrootc,gcomm_pft)
      end if

      ! deadcrootc_storage
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'deadcrootc_storage_13',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'deadcrootc_storage_13') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
             'initializing C13 deadcrootc_storage_13 with atmospheric c13 value'
            do i = begp , endp
              if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1.D0) then
                pcisos%deadcrootc_storage(i) = &
                        pcbulks%deadcrootc_storage(i) * c3_r2
              else
                pcisos%deadcrootc_storage(i) = &
                        pcbulks%deadcrootc_storage(i) * c4_r2
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'deadcrootc_storage_13', &
                  pptr%pc13s%deadcrootc_storage,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'deadcrootc_storage_13', &
                pptr%pc13s%deadcrootc_storage,gcomm_pft)
      end if

      ! deadcrootc_xfer
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'deadcrootc_xfer_13',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'deadcrootc_xfer_13') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
             'initializing C13 deadcrootc_xfer_13 with atmospheric c13 value'
            do i = begp , endp
              if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1.D0) then
                pcisos%deadcrootc_xfer(i) = pcbulks%deadcrootc_xfer(i) * c3_r2
              else
                pcisos%deadcrootc_xfer(i) = pcbulks%deadcrootc_xfer(i) * c4_r2
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'deadcrootc_xfer_13', &
                  pptr%pc13s%deadcrootc_xfer,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'deadcrootc_xfer_13', &
                pptr%pc13s%deadcrootc_xfer,gcomm_pft)
      end if

      ! gresp_storage
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'gresp_storage_13',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'gresp_storage_13') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
             'initializing C13 gresp_storage_13 with atmospheric c13 value'
            do i = begp , endp
              if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1.D0) then
                pcisos%gresp_storage(i) = pcbulks%gresp_storage(i) * c3_r2
              else
                pcisos%gresp_storage(i) = pcbulks%gresp_storage(i) * c4_r2
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'gresp_storage_13', &
                  pptr%pc13s%gresp_storage,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'gresp_storage_13', &
                pptr%pc13s%gresp_storage,gcomm_pft)
      end if

      ! gresp_xfer
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'gresp_xfer_13',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'gresp_xfer_13') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
             'initializing C13 gresp_xfer_13 with atmospheric c13 value'
            do i = begp , endp
              if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1.D0) then
                pcisos%gresp_xfer(i) = pcbulks%gresp_xfer(i) * c3_r2
              else
                pcisos%gresp_xfer(i) = pcbulks%gresp_xfer(i) * c4_r2
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'gresp_xfer_13', &
                  pptr%pc13s%gresp_xfer,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'gresp_xfer_13', &
                pptr%pc13s%gresp_xfer,gcomm_pft)
      end if

      ! cpool
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'cpool_13',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'cpool_13') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
             'initializing C13 cpool_13 with atmospheric c13 value'
            do i = begp , endp
              if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1.D0) then
                pcisos%cpool(i) = pcbulks%cpool(i) * c3_r2
              else
                pcisos%cpool(i) = pcbulks%cpool(i) * c4_r2
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'cpool_13',pptr%pc13s%cpool,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'cpool_13',pptr%pc13s%cpool,gcomm_pft)
      end if

      ! xsmrpool
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'xsmrpool_13',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'xsmrpool_13') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
             'initializing C13 xsmrpool_13 with atmospheric c13 value'
            do i = begp , endp
              if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1.D0) then
                pcisos%xsmrpool(i) = pcbulks%xsmrpool(i) * c3_r2
              else
                pcisos%xsmrpool(i) = pcbulks%xsmrpool(i) * c4_r2
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'xsmrpool_13',pptr%pc13s%xsmrpool,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'xsmrpool_13',pptr%pc13s%xsmrpool,gcomm_pft)
      end if

      ! pft_ctrunc
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'pft_ctrunc_13',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'pft_ctrunc_13') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
             'initializing C13 pft_ctrunc_13 with atmospheric c13 value'
            do i = begp , endp
              if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1.D0) then
                pcisos%pft_ctrunc(i) = pcbulks%pft_ctrunc(i) * c3_r2
              else
                pcisos%pft_ctrunc(i) = pcbulks%pft_ctrunc(i) * c4_r2
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'pft_ctrunc_13',pptr%pc13s%pft_ctrunc,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'pft_ctrunc_13',pptr%pc13s%pft_ctrunc,gcomm_pft)
      end if

      ! totvegc
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'totvegc_13',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'totvegc_13') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
             'initializing C13 totvegc_13 with atmospheric c13 value'
            do i = begp , endp
              if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1.D0) then
                pcisos%totvegc(i) = pcbulks%totvegc(i) * c3_r2
              else
                pcisos%totvegc(i) = pcbulks%totvegc(i) * c4_r2
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'totvegc_13',pptr%pc13s%totvegc,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'totvegc_13',pptr%pc13s%totvegc,gcomm_pft)
      end if
    end if

    if ( use_c14 ) then
      !--------------------------------
      ! C14 pft carbon state variables
      !--------------------------------

      ! leafc
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'leafc_14',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'leafc_14') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
             'initializing C14 leafc_14 with atmospheric c14 value'
            do i = begp , endp
              if (pptr%pcs%leafc(i) /= spval .and. &
                  pptr%pcs%leafc(i) /= nan ) then
                pptr%pc14s%leafc(i) = pptr%pcs%leafc(i) * c14ratio
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'leafc_14',pptr%pc14s%leafc,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'leafc_14',pptr%pc14s%leafc,gcomm_pft)
      end if

      ! leafc_storage
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'leafc_storage_14',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'leafc_storage_14') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
             'initializing C14 leafc_storage_14 with atmospheric c14 value'
            do i = begp , endp
              if (pptr%pcs%leafc_storage(i) /= spval .and. &
                  pptr%pcs%leafc_storage(i) /= nan ) then
                pptr%pc14s%leafc_storage(i) = &
                        pptr%pcs%leafc_storage(i) * c14ratio
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'leafc_storage_14', &
                  pptr%pc14s%leafc_storage,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'leafc_storage_14', &
                pptr%pc14s%leafc_storage,gcomm_pft)
      end if

      ! leafc_xfer
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'leafc_xfer_14',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'leafc_xfer_14') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
             'initializing C14 leafc_xfer_14 with atmospheric c14 value'
            do i = begp , endp
              if (pptr%pcs%leafc_xfer(i) /= spval .and. &
                  pptr%pcs%leafc_xfer(i) /= nan ) then
                pptr%pc14s%leafc_xfer(i) = pptr%pcs%leafc_xfer(i) * c14ratio
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'leafc_xfer_14', &
                  pptr%pc14s%leafc_xfer,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'leafc_xfer_14', &
                pptr%pc14s%leafc_xfer,gcomm_pft)
      end if

      ! frootc
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'frootc_14',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'frootc_14') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
             'initializing C14 frootc_14 with atmospheric c14 value'
            do i = begp , endp
              if (pptr%pcs%frootc(i) /= spval .and. &
                  pptr%pcs%frootc(i) /= nan ) then
                pptr%pc14s%frootc(i) = pptr%pcs%frootc(i) * c14ratio
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'frootc_14', &
                  pptr%pc14s%frootc,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'frootc_14', &
                pptr%pc14s%frootc,gcomm_pft)
      end if

      ! frootc_storage
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'frootc_storage_14',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'frootc_storage_14') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
             'initializing C14 frootc_storage_14 with atmospheric c14 value'
            do i = begp , endp
              if (pptr%pcs%frootc_storage(i) /= spval .and. &
                  pptr%pcs%frootc_storage(i) /= nan ) then
                pptr%pc14s%frootc_storage(i) = &
                        pptr%pcs%frootc_storage(i) * c14ratio
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'frootc_storage_14', &
                  pptr%pc14s%frootc_storage,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'frootc_storage_14', &
                pptr%pc14s%frootc_storage,gcomm_pft)
      end if

      !frootc_xfer
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'frootc_xfer_14',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'frootc_xfer_14') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
             'initializing C14 frootc_xfer_14 with atmospheric c14 value'
            do i = begp , endp
              if (pptr%pcs%frootc_xfer(i) /= spval .and. &
                  pptr%pcs%frootc_xfer(i) /= nan ) then
                pptr%pc14s%frootc_xfer(i) = pptr%pcs%frootc_xfer(i) * c14ratio
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'frootc_xfer_14', &
                  pptr%pc14s%frootc_xfer,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'frootc_xfer_14', &
                pptr%pc14s%frootc_xfer,gcomm_pft)
      end if

      ! livestemc
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'livestemc_14',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'livestemc_14') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
             'initializing C14 livestemc_14 with atmospheric c14 value'
            do i = begp , endp
              if (pptr%pcs%livestemc(i) /= spval .and. &
                  pptr%pcs%livestemc(i) /= nan ) then
                pptr%pc14s%livestemc(i) = pptr%pcs%livestemc(i) * c14ratio
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'livestemc_14', &
                  pptr%pc14s%livestemc,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'livestemc_14', &
                pptr%pc14s%livestemc,gcomm_pft)
      end if

      ! livestemc_storage
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'livestemc_storage_14',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'livestemc_storage_14') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
             'initializing C14 livestemc_storage_14 with atmospheric c14 value'
            do i = begp , endp
              if (pptr%pcs%livestemc_storage(i) /= spval .and. &
                  pptr%pcs%livestemc_storage(i) /= nan ) then
                pptr%pc14s%livestemc_storage(i) = &
                        pptr%pcs%livestemc_storage(i) * c14ratio
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'livestemc_storage_14', &
                  pptr%pc14s%livestemc_storage,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'livestemc_storage_14', &
                pptr%pc14s%livestemc_storage,gcomm_pft)
      end if

      ! livestemc_xfer
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'livestemc_xfer_14',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'livestemc_xfer_14') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
             'initializing C14 livestemc_xfer_14 with atmospheric c14 value'
            do i = begp , endp
              if (pptr%pcs%livestemc_xfer(i) /= spval .and. &
                  pptr%pcs%livestemc_xfer(i) /= nan ) then
                pptr%pc14s%livestemc_xfer(i) = &
                        pptr%pcs%livestemc_xfer(i) * c14ratio
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'livestemc_xfer_14', &
                  pptr%pc14s%livestemc_xfer,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'livestemc_xfer_14', &
                pptr%pc14s%livestemc_xfer,gcomm_pft)
      end if

      ! deadstemc
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'deadstemc_14',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'deadstemc_14') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
             'initializing C14 deadstemc_14 with atmospheric c14 value'
            do i = begp , endp
              if (pptr%pcs%deadstemc(i) /= spval .and. &
                  pptr%pcs%deadstemc(i) /= nan ) then
                pptr%pc14s%deadstemc(i) = pptr%pcs%deadstemc(i) * c14ratio
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'deadstemc_14', &
                  pptr%pc14s%deadstemc,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'deadstemc_14', &
                pptr%pc14s%deadstemc,gcomm_pft)
      end if

      ! deadstemc_storage
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'deadstemc_storage_14',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'deadstemc_storage_14') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
             'initializing C14 deadstemc_storage_14 with atmospheric c14 value'
            do i = begp , endp
              if (pptr%pcs%deadstemc_storage(i) /= spval .and. &
                  pptr%pcs%deadstemc_storage(i) /= nan ) then
                pptr%pc14s%deadstemc_storage(i) = &
                        pptr%pcs%deadstemc_storage(i) * c14ratio
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'deadstemc_storage_14', &
                  pptr%pc14s%deadstemc_storage,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'deadstemc_storage_14', &
                pptr%pc14s%deadstemc_storage,gcomm_pft)
      end if

      ! deadstemc_xfer
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'deadstemc_xfer_14',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'deadstemc_xfer_14') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
             'initializing C14 deadstemc_xfer_14 with atmospheric c14 value'
            do i = begp , endp
              if (pptr%pcs%deadstemc_xfer(i) /= spval .and. &
                  pptr%pcs%deadstemc_xfer(i) /= nan ) then
                pptr%pc14s%deadstemc_xfer(i) = &
                        pptr%pcs%deadstemc_xfer(i) * c14ratio
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'deadstemc_xfer_14', &
                  pptr%pc14s%deadstemc_xfer,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'deadstemc_xfer_14', &
                pptr%pc14s%deadstemc_xfer,gcomm_pft)
      end if

      ! livecrootc
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'livecrootc_14',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'livecrootc_14') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
             'initializing C14 livecrootc_14 with atmospheric c14 value'
            do i = begp , endp
              if (pptr%pcs%livecrootc(i) /= spval .and. &
                  pptr%pcs%livecrootc(i) /= nan ) then
                pptr%pc14s%livecrootc(i) = pptr%pcs%livecrootc(i) * c14ratio
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'livecrootc_14', &
                  pptr%pc14s%livecrootc,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'livecrootc_14', &
                pptr%pc14s%livecrootc,gcomm_pft)
      end if

      ! livecrootc_storage
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'livecrootc_storage_14',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'livecrootc_storage_14') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
             'initializing C14 livecrootc_storage_14 with atmospheric c14 value'
            do i = begp , endp
              if (pptr%pcs%livecrootc_storage(i) /= spval .and. &
                  pptr%pcs%livecrootc_storage(i) /= nan ) then
                pptr%pc14s%livecrootc_storage(i) = &
                        pptr%pcs%livecrootc_storage(i) * c14ratio
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'livecrootc_storage_14', &
                  pptr%pc14s%livecrootc_storage,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'livecrootc_storage_14', &
                pptr%pc14s%livecrootc_storage,gcomm_pft)
      end if

      ! livecrootc_xfer
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'livecrootc_xfer_14',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'livecrootc_xfer_14') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
             'initializing C14 livecrootc_xfer_14 with atmospheric c14 value'
            do i = begp , endp
              if (pptr%pcs%livecrootc_xfer(i) /= spval .and. &
                  pptr%pcs%livecrootc_xfer(i) /= nan ) then
                pptr%pc14s%livecrootc_xfer(i) = &
                        pptr%pcs%livecrootc_xfer(i) * c14ratio
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'livecrootc_xfer_14', &
                  pptr%pc14s%livecrootc_xfer,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'livecrootc_xfer_14', &
                pptr%pc14s%livecrootc_xfer,gcomm_pft)
      end if

      ! deadcrootc
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'deadcrootc_14',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'deadcrootc_14') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
             'initializing C14 deadcrootc_14 with atmospheric c14 value'
            do i = begp , endp
              if (pptr%pcs%deadcrootc(i) /= spval .and. &
                  pptr%pcs%deadcrootc(i) /= nan ) then
                pptr%pc14s%deadcrootc(i) = pptr%pcs%deadcrootc(i) * c14ratio
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'deadcrootc_14', &
                  pptr%pc14s%deadcrootc,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'deadcrootc_14', &
                pptr%pc14s%deadcrootc,gcomm_pft)
      end if

      ! deadcrootc_storage
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'deadcrootc_storage_14',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'deadcrootc_storage_14') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
             'initializing C14 deadcrootc_storage_14 with atmospheric c14 value'
            do i = begp , endp
              if (pptr%pcs%deadcrootc_storage(i) /= spval .and. &
                  pptr%pcs%deadcrootc_storage(i) /= nan ) then
                pptr%pc14s%deadcrootc_storage(i) = &
                        pptr%pcs%deadcrootc_storage(i) * c14ratio
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'deadcrootc_storage_14', &
                  pptr%pc14s%deadcrootc_storage,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'deadcrootc_storage_14', &
                pptr%pc14s%deadcrootc_storage,gcomm_pft)
      end if

      ! deadcrootc_xfer
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'deadcrootc_xfer_14',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'deadcrootc_xfer_14') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
             'initializing C14 deadcrootc_xfer_14 with atmospheric c14 value'
            do i = begp,endp
              if (pptr%pcs%deadcrootc_xfer(i) /= spval .and. &
                  pptr%pcs%deadcrootc_xfer(i) /= nan ) then
                pptr%pc14s%deadcrootc_xfer(i) = &
                        pptr%pcs%deadcrootc_xfer(i) * c14ratio
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'deadcrootc_xfer_14', &
                  pptr%pc14s%deadcrootc_xfer,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'deadcrootc_xfer_14', &
                pptr%pc14s%deadcrootc_xfer,gcomm_pft)
      end if

      ! gresp_storage
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'gresp_storage_14',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'gresp_storage_14') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
             'initializing C14 gresp_storage_14 with atmospheric c14 value'
            do i = begp , endp
              if (pptr%pcs%gresp_storage(i) /= spval .and. &
                  pptr%pcs%gresp_storage(i) /= nan ) then
                pptr%pc14s%gresp_storage(i) = &
                        pptr%pcs%gresp_storage(i) * c14ratio
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'gresp_storage_14', &
                  pptr%pc14s%gresp_storage,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'gresp_storage_14', &
                pptr%pc14s%gresp_storage,gcomm_pft)
      end if

      ! gresp_xfer
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'gresp_xfer_14',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'gresp_xfer_14') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
             'initializing C14 gresp_xfer_14 with atmospheric c14 value'
            do i = begp,endp
              if (pptr%pcs%gresp_xfer(i) /= spval .and. &
                  pptr%pcs%gresp_xfer(i) /= nan ) then
                pptr%pc14s%gresp_xfer(i) = pptr%pcs%gresp_xfer(i) * c14ratio
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'gresp_xfer_14', &
                  pptr%pc14s%gresp_xfer,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'gresp_xfer_14', &
                pptr%pc14s%gresp_xfer,gcomm_pft)
      end if

      ! cpool
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'cpool_14',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'cpool_14') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
             'initializing C14 cpool_14 with atmospheric c14 value'
            do i = begp , endp
              if (pptr%pcs%cpool(i) /= spval .and. &
                  pptr%pcs%cpool(i) /= nan ) then
                pptr%pc14s%cpool(i) = pptr%pcs%cpool(i) * c14ratio
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'cpool_14',pptr%pc14s%cpool,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'cpool_14',pptr%pc14s%cpool,gcomm_pft)
      end if

      ! xsmrpool
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'xsmrpool_14',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'xsmrpool_14') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
             'initializing C14 xsmrpool_14 with atmospheric c14 value'
            do i = begp , endp
              if (pptr%pcs%xsmrpool(i) /= spval .and. &
                  pptr%pcs%xsmrpool(i) /= nan ) then
                pptr%pc14s%xsmrpool(i) = pptr%pcs%xsmrpool(i) * c14ratio
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'xsmrpool_14',pptr%pc14s%xsmrpool,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'xsmrpool_14',pptr%pc14s%xsmrpool,gcomm_pft)
      end if

      ! pft_ctrunc
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'pft_ctrunc_14',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'pft_ctrunc_14') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
             'initializing C14 pft_ctrunc_14 with atmospheric c14 value'
            do i = begp,endp
              if (pptr%pcs%pft_ctrunc(i) /= spval .and. &
                  pptr%pcs%pft_ctrunc(i) /= nan ) then
                pptr%pc14s%pft_ctrunc(i) = pptr%pcs%pft_ctrunc(i) * c14ratio
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'pft_ctrunc_14',pptr%pc14s%pft_ctrunc,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'pft_ctrunc_14',pptr%pc14s%pft_ctrunc,gcomm_pft)
      end if

      ! totvegc
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'totvegc_14',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'totvegc_14') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
             'initializing C14 totvegc_14 with atmospheric c14 value'
            do i = begp,endp
              if (pptr%pcs%totvegc(i) /= spval .and. &
                  pptr%pcs%totvegc(i) /= nan ) then
                pptr%pc14s%totvegc(i) = pptr%pcs%totvegc(i) * c14ratio
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'totvegc_14',pptr%pc14s%totvegc,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'totvegc_14',pptr%pc14s%totvegc,gcomm_pft)
      end if

      ! rc14_atm
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'rc14_atm',(/'pft'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( .not. clm_check_var(ncid,'rc14_atm') ) then
          if ( ktau == 0 ) then
            write(stdout,*) &
             'initializing C14 rc14_atm with atmospheric c14 value'
            do i = begp,endp
              pptr%pepv%rc14_atm(i) = c14ratio
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'rc14_atm',pptr%pepv%rc14_atm,gcomm_pft)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'rc14_atm',pptr%pepv%rc14_atm,gcomm_pft)
      end if
    end if ! use C14

    !--------------------------------
    ! pft nitrogen state variables
    !--------------------------------

    ! leafn
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'leafn',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. &
           .not. clm_check_var(ncid,'leafn') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'leafn',pptr%pns%leafn,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'leafn',pptr%pns%leafn,gcomm_pft)
    end if

    ! leafn_storage
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'leafn_storage',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. &
           .not. clm_check_var(ncid,'leafn_storage') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'leafn_storage',pptr%pns%leafn_storage,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'leafn_storage',pptr%pns%leafn_storage,gcomm_pft)
    end if

    ! leafn_xfer
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'leafn_xfer',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. &
           .not. clm_check_var(ncid,'leafn_xfer') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'leafn_xfer',pptr%pns%leafn_xfer,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'leafn_xfer',pptr%pns%leafn_xfer,gcomm_pft)
    end if

    ! frootn
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'frootn',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. &
           .not. clm_check_var(ncid,'frootn') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'frootn',pptr%pns%frootn,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'frootn',pptr%pns%frootn,gcomm_pft)
    end if

    ! frootn_storage
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'frootn_storage',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. &
           .not. clm_check_var(ncid,'frootn_storage') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'frootn_storage', &
                pptr%pns%frootn_storage,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'frootn_storage', &
              pptr%pns%frootn_storage,gcomm_pft)
    end if

    ! frootn_xfer
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'frootn_xfer',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. &
           .not. clm_check_var(ncid,'frootn_xfer') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'frootn_xfer',pptr%pns%frootn_xfer,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'frootn_xfer',pptr%pns%frootn_xfer,gcomm_pft)
    end if

    ! livestemn
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'livestemn',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. &
           .not. clm_check_var(ncid,'livestemn') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'livestemn',pptr%pns%livestemn,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'livestemn',pptr%pns%livestemn,gcomm_pft)
    end if

    ! livestemn_storage
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'livestemn_storage',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. &
           .not. clm_check_var(ncid,'livestemn_storage') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'livestemn_storage', &
                pptr%pns%livestemn_storage,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'livestemn_storage', &
              pptr%pns%livestemn_storage,gcomm_pft)
    end if

    ! livestemn_xfer
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'livestemn_xfer',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. &
           .not. clm_check_var(ncid,'livestemn_xfer') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'livestemn_xfer', &
                pptr%pns%livestemn_xfer,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'livestemn_xfer', &
              pptr%pns%livestemn_xfer,gcomm_pft)
    end if

    ! deadstemn
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'deadstemn',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. &
           .not. clm_check_var(ncid,'deadstemn') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'deadstemn',pptr%pns%deadstemn,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'deadstemn',pptr%pns%deadstemn,gcomm_pft)
    end if

    !deadstemn_storage
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'deadstemn_storage',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. &
           .not. clm_check_var(ncid,'deadstemn_storage') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'deadstemn_storage', &
                pptr%pns%deadstemn_storage,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'deadstemn_storage', &
              pptr%pns%deadstemn_storage,gcomm_pft)
    end if

    !deadstemn_xfer
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'deadstemn_xfer',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. &
           .not. clm_check_var(ncid,'deadstemn_xfer') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'deadstemn_xfer', &
                pptr%pns%deadstemn_xfer,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'deadstemn_xfer', &
              pptr%pns%deadstemn_xfer,gcomm_pft)
    end if

    ! livecrootn
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'livecrootn',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. &
           .not. clm_check_var(ncid,'livecrootn') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'livecrootn',pptr%pns%livecrootn,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'livecrootn',pptr%pns%livecrootn,gcomm_pft)
    end if

    ! livecrootn_storage
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'livecrootn_storage',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. &
           .not. clm_check_var(ncid,'livecrootn_storage') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'livecrootn_storage', &
                pptr%pns%livecrootn_storage,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'livecrootn_storage', &
              pptr%pns%livecrootn_storage,gcomm_pft)
    end if

    ! livecrootn_xfer
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'livecrootn_xfer',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. &
           .not. clm_check_var(ncid,'livecrootn_xfer') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'livecrootn_xfer', &
                pptr%pns%livecrootn_xfer,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'livecrootn_xfer', &
              pptr%pns%livecrootn_xfer,gcomm_pft)
    end if

    ! deadcrootn
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'deadcrootn',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. &
           .not. clm_check_var(ncid,'deadcrootn') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'deadcrootn',pptr%pns%deadcrootn,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'deadcrootn',pptr%pns%deadcrootn,gcomm_pft)
    end if

    ! deadcrootn_storage
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'deadcrootn_storage',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. &
           .not. clm_check_var(ncid,'deadcrootn_storage') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'deadcrootn_storage', &
                pptr%pns%deadcrootn_storage,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'deadcrootn_storage', &
              pptr%pns%deadcrootn_storage,gcomm_pft)
    end if

    ! deadcrootn_xfer
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'deadcrootn_xfer',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. &
           .not. clm_check_var(ncid,'deadcrootn_xfer') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'deadcrootn_xfer', &
                pptr%pns%deadcrootn_xfer,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'deadcrootn_xfer', &
              pptr%pns%deadcrootn_xfer,gcomm_pft)
    end if

    !retransn
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'retransn',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. &
           .not. clm_check_var(ncid,'retransn') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'retransn',pptr%pns%retransn,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'retransn',pptr%pns%retransn,gcomm_pft)
    end if

    ! npool
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'npool',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. &
           .not. clm_check_var(ncid,'npool') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'npool',pptr%pns%npool,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'npool',pptr%pns%npool,gcomm_pft)
    end if

    ! pft_ntrunc
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'pft_ntrunc',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. &
           .not. clm_check_var(ncid,'pft_ntrunc') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'pft_ntrunc',pptr%pns%pft_ntrunc,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'pft_ntrunc',pptr%pns%pft_ntrunc,gcomm_pft)
    end if

    !--------------------------------
    ! column physical state variables
    !--------------------------------

    ! decl
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'decl',(/'column'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. &
           .not. clm_check_var(ncid,'decl') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'decl',cptr%cps%decl,gcomm_column)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'decl',cptr%cps%decl,gcomm_column)
    end if

    ! fpi
    call cnrest_addfld_decomp(ncid=ncid, varname='fpi',     &
           longname='fraction of potential immobilization', &
           units='unitless', flag=flag, data_rl=cptr%cps%fpi_vr)

#ifdef VERTSOILC
    ! som_adv_coef
    call cnrest_addfld_decomp(ncid=ncid, varname='som_adv_coef', &
           longname='SOM advective flux', units='m/s',           &
           flag=flag, data_rl=cptr%cps%som_adv_coef)

    ! som_diffus_coef
    call cnrest_addfld_decomp(ncid=ncid, varname='som_diffus_coef', &
           longname='SOM diffusivity due to bio/cryo-turbation',   &
           units='m^2/s',  flag=flag, data_rl=cptr%cps%som_diffus_coef)
#endif

!#ifdef NITRIF_DENITRIF
!    ! tmean_monthly_max
!    call cnrest_addfld_decomp(ncid=ncid, varname='tmean_monthly_max', &
!              longname='', units='', flag=flag, &
!              data_rl=cptr%cps%tmean_monthly_max)
!
!    ! tmean_monthly
!    call cnrest_addfld_decomp(ncid=ncid, varname='tmean_monthly', &
!              longname='', units='', flag=flag, &
!              data_rl=cptr%cps%tmean_monthly)
!#end if

    ! fpg
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'fpg',(/'column'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'fpg') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'fpg',cptr%cps%fpg,gcomm_column)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'fpg',cptr%cps%fpg,gcomm_column)
    end if

    ! annsum_counter
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'annsum_counter',(/'column'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'annsum_counter') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'annsum_counter', &
                cptr%cps%annsum_counter,gcomm_column)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'annsum_counter', &
              cptr%cps%annsum_counter,gcomm_column)
    end if

    ! cannsum_npp
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'cannsum_npp',(/'column'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'cannsum_npp') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'cannsum_npp', &
                cptr%cps%cannsum_npp,gcomm_column)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'cannsum_npp', &
              cptr%cps%cannsum_npp,gcomm_column)
    end if

    ! col_lag_npp
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'col_lag_npp',(/'column'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'col_lag_npp') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'col_lag_npp', &
                cptr%cps%col_lag_npp,gcomm_column)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'col_lag_npp', &
              cptr%cps%col_lag_npp,gcomm_column)
    end if

    ! cannavg_t2m
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'cannavg_t2m',(/'column'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'cannavg_t2m') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'cannavg_t2m', &
                cptr%cps%cannavg_t2m,gcomm_column)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'cannavg_t2m', &
              cptr%cps%cannavg_t2m,gcomm_column)
    end if

    ! for fire model changed by F. Li and S. Levis
    !  burndate
    if (flag == 'define') then
      call clm_addvar(clmvar_integer,ncid,'burndate',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'burndate') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'burndate',pptr%pps%burndate,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'burndate',pptr%pps%burndate,gcomm_pft)
    end if

    !lfc
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'lfc',(/'column'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'lfc') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'lfc',cptr%cps%lfc,gcomm_column)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'lfc',cptr%cps%lfc,gcomm_column)
    end if

    !wf
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'wf',(/'column'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'wf') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'wf',cptr%cps%wf,gcomm_column)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'wf',cptr%cps%wf,gcomm_column)
    end if

    !btran2
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'btran2',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'btran2') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'btran2',pptr%pps%btran2,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'btran2',pptr%pps%btran2,gcomm_pft)
    end if

    !farea_burned
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'farea_burned',(/'column'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'farea_burned') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'farea_burned',cptr%cps%farea_burned,gcomm_column)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'farea_burned',cptr%cps%farea_burned,gcomm_column)
    end if

    !baf_crop
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'baf_crop',(/'column'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'baf_crop') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'baf_crop',cptr%cps%baf_crop,gcomm_column)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'baf_crop',cptr%cps%baf_crop,gcomm_column)
    end if

    !baf_peatf
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'baf_peatf',(/'column'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'baf_peatf') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'baf_peatf',cptr%cps%baf_peatf,gcomm_column)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'baf_peatf',cptr%cps%baf_peatf,gcomm_column)
    end if

    !fbac
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'fbac',(/'column'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'fbac') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'fbac',cptr%cps%fbac,gcomm_column)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'fbac',cptr%cps%fbac,gcomm_column)
    end if

    !fbac1
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'fbac1',(/'column'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'fbac1') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'fbac1',cptr%cps%fbac1,gcomm_column)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'fbac1',cptr%cps%fbac1,gcomm_column)
    end if

    !--------------------------------
    ! column carbon state variables
    !--------------------------------

    do k = 1, ndecomp_pools
      ptr2d => cptr%ccs%decomp_cpools_vr(:,:,k)
      varname = trim(decomp_cascade_con%decomp_pool_name_restart(k))//'c'
      call cnrest_addfld_decomp(ncid=ncid, varname=varname, longname='', &
              units='', flag=flag, data_rl=ptr2d)
    end do

    ! col_ctrunc
    varname = 'col_ctrunc'
    call cnrest_addfld_decomp(ncid=ncid, varname=varname, longname='', &
            units='', flag=flag, data_rl=cptr%ccs%col_ctrunc_vr)

    ! ! nfixation_prof
    ! call cnrest_addfld_decomp(ncid=ncid, varname='nfixation_prof', &
    !        longname='', units='', flag=flag, &
    !        data_rl=cptr%cps%nfixation_prof)

    ! ! ndep_prof
    ! call cnrest_addfld_decomp(ncid=ncid, varname='ndep_prof', &
    !        longname='', units='', flag=flag, &
    !        data_rl=cptr%cps%ndep_prof)

    ! altmax
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'altmax',(/'column'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'altmax') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'altmax',cptr%cps%altmax,gcomm_column)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'altmax',cptr%cps%altmax,gcomm_column)
    end if

    ! altmax_lastyear
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'altmax_lastyear',(/'column'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'altmax_lastyear') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'altmax_lastyear', &
                cptr%cps%altmax_lastyear,gcomm_column)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'altmax_lastyear', &
              cptr%cps%altmax_lastyear,gcomm_column)
    end if

    ! altmax_indx
    if (flag == 'define') then
      call clm_addvar(clmvar_integer,ncid,'altmax_indx',(/'column'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'altmax_indx') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'altmax_indx', &
                cptr%cps%altmax_indx,gcomm_column)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'altmax_indx', &
              cptr%cps%altmax_indx,gcomm_column)
    end if

    ! altmax_lastyear_indx
    if (flag == 'define') then
      call clm_addvar(clmvar_integer,ncid,'altmax_lastyear_indx',(/'column'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. &
           clm_check_var(ncid,'altmax_lastyear_indx') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'altmax_lastyear_indx', &
                cptr%cps%altmax_lastyear_indx,gcomm_column)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'altmax_lastyear_indx', &
              cptr%cps%altmax_lastyear_indx,gcomm_column)
    end if

    ! seedc
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'seedc',(/'column'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'seedc') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'seedc',cptr%ccs%seedc,gcomm_column)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'seedc',cptr%ccs%seedc,gcomm_column)
    end if

    ! totlitc
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'totlitc',(/'column'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'totlitc') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'totlitc',cptr%ccs%totlitc,gcomm_column)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'totlitc',cptr%ccs%totlitc,gcomm_column)
    end if

    ! totcolc
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'totcolc',(/'column'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'totcolc') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'totcolc',cptr%ccs%totcolc,gcomm_column)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'totcolc',cptr%ccs%totcolc,gcomm_column)
    end if
 
    ! prod10c
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'prod10c',(/'column'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'prod10c') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'prod10c',cptr%ccs%prod10c,gcomm_column)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'prod10c',cptr%ccs%prod10c,gcomm_column)
    end if

    ! prod100c
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'prod100c',(/'column'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'prod100c') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'prod100c',cptr%ccs%prod100c,gcomm_column)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'prod100c',cptr%ccs%prod100c,gcomm_column)
    end if

    ! totsomc
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'totsomc',(/'column'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'totsomc') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'totsomc',cptr%ccs%totsomc,gcomm_column)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'totsomc',cptr%ccs%totsomc,gcomm_column)
    end if

    if ( use_c13 ) then
      !--------------------------------
      ! C13 column carbon state variables
      !--------------------------------

      do k = 1 , ndecomp_pools
        ptr2d => cptr%cc13s%decomp_cpools_vr(:,:,k)
        lstop = .false.
        varname=trim(decomp_cascade_con%decomp_pool_name_restart(k))//'c_13'
        call cnrest_addfld_decomp(ncid=ncid, varname=varname, longname='', &
                units='', flag=flag, data_rl=ptr2d, lstop=lstop)
        if ( flag == 'read' .and. lstop ) then
          write(stdout,*) 'init decomp_cpools_vr with c13 value for: '//varname
          do i = begc , endc
            do j = 1 , nlevdecomp
              if (cptr%ccs%decomp_cpools_vr(i,j,k) /= spval .and. &
                  cptr%ccs%decomp_cpools_vr(i,j,k) /= nan ) then
                cptr%cc13s%decomp_cpools_vr(i,j,k) = &
                        cptr%ccs%decomp_cpools_vr(i,j,k) * c3_r2
              end if
            end do
          end do
        end if
      end do

      ! seedc
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'seedc_13',(/'column'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( ktau == 0 ) then
          if ( .not. clm_check_var(ncid,'seedc_13') ) then
            do i = begc,endc
              if (cptr%ccs%seedc(i) /= spval .and. &
                  cptr%ccs%seedc(i) /= nan ) then
                cptr%cc13s%seedc(i) = cptr%ccs%seedc(i) * c3_r2
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'seedc_13',cptr%cc13s%seedc,gcomm_column)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'seedc_13',cptr%cc13s%seedc,gcomm_column)
      end if

      ! col_ctrunc_13
      call cnrest_addfld_decomp(ncid=ncid, varname='col_ctrunc_13_vr', &
              longname='', units='', flag=flag, &
              data_rl=cptr%cc13s%col_ctrunc_vr)

      ! totlitc
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'totlitc_13',(/'column'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( ktau == 0 ) then
          if ( .not. clm_check_var(ncid,'totlitc_13') ) then
            do i = begc , endc
              if (cptr%ccs%totlitc(i) /= spval .and. &
                  cptr%ccs%totlitc(i) /= nan ) then
                cptr%cc13s%totlitc(i) = cptr%ccs%totlitc(i) * c3_r2
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'totlitc_13',cptr%cc13s%totlitc,gcomm_column)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'totlitc_13',cptr%cc13s%totlitc,gcomm_column)
      end if

      ! totcolc
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'totcolc_13',(/'column'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( ktau == 0 ) then
          if ( .not. clm_check_var(ncid,'totcolc_13') ) then
            do i = begc , endc
              if (cptr%ccs%totcolc(i) /= spval .and. &
                  cptr%ccs%totcolc(i) /= nan ) then
                cptr%cc13s%totcolc(i) = cptr%ccs%totcolc(i) * c3_r2
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'totcolc_13',cptr%cc13s%totcolc,gcomm_column)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'totcolc_13',cptr%cc13s%totcolc,gcomm_column)
      end if

      ! prod10c
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'prod10c_13',(/'column'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( ktau == 0 ) then
          if ( .not. clm_check_var(ncid,'prod10c_13') ) then
            do i = begc , endc
              if (cptr%ccs%prod10c(i) /= spval .and. &
                  cptr%ccs%prod10c(i) /= nan ) then
                cptr%cc13s%prod10c(i) = cptr%ccs%prod10c(i) * c3_r2
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'prod10c_13',cptr%cc13s%prod10c,gcomm_column)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'prod10c_13',cptr%cc13s%prod10c,gcomm_column)
      end if

      ! prod100c
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'prod100c_13',(/'column'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( ktau == 0 ) then
          if ( .not. clm_check_var(ncid,'prod100c_13') ) then
            do i = begc , endc
              if (cptr%ccs%prod100c(i) /= spval .and. &
                  cptr%ccs%prod100c(i) /= nan ) then
                cptr%cc13s%prod100c(i) = cptr%ccs%prod100c(i) * c3_r2
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'prod100c_13',cptr%cc13s%prod100c,gcomm_column)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'prod100c_13',cptr%cc13s%prod100c,gcomm_column)
      end if

    end if

    if ( use_c14 ) then
      !--------------------------------
      ! C14 column carbon state variables
      !--------------------------------

      do k = 1 , ndecomp_pools
        ptr2d => cptr%cc14s%decomp_cpools_vr(:,:,k)
        lstop = .false.
        varname=trim(decomp_cascade_con%decomp_pool_name_restart(k))//'c_14'
        call cnrest_addfld_decomp(ncid=ncid, varname=varname, &
                longname='', units='', flag=flag, data_rl=ptr2d, lstop=lstop)
        if ( flag == 'read' .and. lstop ) then
          write(stdout,*) 'init decomp_cpools_vr with C14 value for: '//varname
          do i = begc , endc
            do j = 1 , nlevdecomp
              if (cptr%ccs%decomp_cpools_vr(i,j,k) /= spval .and. &
                  cptr%ccs%decomp_cpools_vr(i,j,k) /= nan ) then
                cptr%cc14s%decomp_cpools_vr(i,j,k) = &
                        cptr%ccs%decomp_cpools_vr(i,j,k) * c14ratio
              end if
            end do
          end do
        end if
      end do

      ! seedc
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'seedc_14',(/'column'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( ktau == 0 ) then
          if ( .not. clm_check_var(ncid,'seedc_14') ) then
            do i = begc , endc
              if (cptr%ccs%seedc(i) /= spval .and. &
                  cptr%ccs%seedc(i) /= nan ) then
                cptr%cc14s%seedc(i) = cptr%ccs%seedc(i) * c14ratio
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'seedc_14',cptr%cc14s%seedc,gcomm_column)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'seedc_14',cptr%cc14s%seedc,gcomm_column)
      end if

      ! col_ctrunc_c14
      call cnrest_addfld_decomp(ncid=ncid, varname='col_ctrunc_14_vr', &
              longname='', units='', flag=flag, &
              data_rl=cptr%cc14s%col_ctrunc_vr)

      ! totlitc
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'totlitc_14',(/'column'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( ktau == 0 ) then
          if ( .not. clm_check_var(ncid,'totlitc_14') ) then
            do i = begc , endc
              if (cptr%ccs%totlitc(i) /= spval .and. &
                  cptr%ccs%totlitc(i) /= nan ) then
                cptr%cc14s%totlitc(i) = cptr%ccs%totlitc(i) * c14ratio
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'totlitc_14',cptr%cc14s%totlitc,gcomm_column)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'totlitc_14',cptr%cc14s%totlitc,gcomm_column)
      end if

      ! totcolc
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'totcolc_14',(/'column'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( ktau == 0 ) then
          if ( .not. clm_check_var(ncid,'totcolc_14') ) then
            do i = begc , endc
              if (cptr%ccs%totcolc(i) /= spval .and. &
                  cptr%ccs%totcolc(i) /= nan ) then
                cptr%cc14s%totcolc(i) = cptr%ccs%totcolc(i) * c14ratio
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'totcolc_14',cptr%cc14s%totcolc,gcomm_column)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'totcolc_14',cptr%cc14s%totcolc,gcomm_column)
      end if

      ! prod10c
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'prod10c_14',(/'column'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( ktau == 0 ) then
          if ( .not. clm_check_var(ncid,'prod10c_14') ) then
            do i = begc , endc
              if (cptr%ccs%prod10c(i) /= spval .and. &
                  cptr%ccs%prod10c(i) /= nan ) then
                cptr%cc14s%prod10c(i) = cptr%ccs%prod10c(i) * c14ratio
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'prod10c_14',cptr%cc14s%prod10c,gcomm_column)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'prod10c_14',cptr%cc14s%prod10c,gcomm_column)
      end if

      ! prod100c
      if (flag == 'define') then
        call clm_addvar(clmvar_double,ncid,'prod100c_14',(/'column'/), &
              long_name='',units='')
      else if (flag == 'read') then
        if ( ktau == 0 ) then
          if ( .not. clm_check_var(ncid,'prod100c_14') ) then
            do i = begc , endc
              if (cptr%ccs%prod100c(i) /= spval .and. &
                  cptr%ccs%prod100c(i) /= nan ) then
                cptr%cc14s%prod100c(i) = cptr%ccs%prod100c(i) * c14ratio
              end if
            end do
          else
            call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
        else
          call clm_readvar(ncid,'prod100c_14',cptr%cc14s%prod100c,gcomm_column)
        end if
      else if (flag == 'write') then
        call clm_writevar(ncid,'prod100c_14',cptr%cc14s%prod100c,gcomm_column)
      end if

    end if

    !--------------------------------
    ! column nitrogen state variables
    !--------------------------------

    ! sminn
    varname='sminn'
    call cnrest_addfld_decomp(ncid=ncid, varname=varname, longname='', &
            units='', flag=flag, data_rl=cptr%cns%sminn_vr)

    ! decomposing N pools
    do k = 1 , ndecomp_pools
      ptr2d => cptr%cns%decomp_npools_vr(:,:,k)
      varname=trim(decomp_cascade_con%decomp_pool_name_restart(k))//'n'
      call cnrest_addfld_decomp(ncid=ncid, varname=varname, longname='', &
              units='', flag=flag, data_rl=ptr2d)
    end do

    ! col_ntrunc
    call cnrest_addfld_decomp(ncid=ncid, varname='col_ntrunc', longname='', &
            units='', flag=flag, data_rl=cptr%cns%col_ntrunc_vr)

#ifdef NITRIF_DENITRIF
    ! f_nit_vr
    call cnrest_addfld_decomp(ncid=ncid, varname='f_nit_vr', &
            longname='soil nitrification flux', units='gN/m3/s', &
            flag=flag, data_rl=cptr%cnf%f_nit_vr)

    ! pot_f_nit_vr
    call cnrest_addfld_decomp(ncid=ncid, varname='pot_f_nit_vr', &
            longname='potential soil nitrification flux', units='gN/m3/s', &
            flag=flag, data_rl=cptr%cnf%pot_f_nit_vr)

    ! smin_no3_vr
    varname = 'smin_no3'
    call cnrest_addfld_decomp(ncid=ncid, varname=varname, longname='', &
            units='', flag=flag, data_rl=cptr%cns%smin_no3_vr)

    ! smin_nh4
    varname = 'smin_nh4'
    call cnrest_addfld_decomp(ncid=ncid, varname=varname, longname='', &
            units='', flag=flag, data_rl=cptr%cns%smin_nh4_vr)
#endif

    ! totcoln
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'totcoln',(/'column'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'totcoln') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'totcoln',cptr%cns%totcoln,gcomm_column)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'totcoln',cptr%cns%totcoln,gcomm_column)
    end if

    ! seedn
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'seedn',(/'column'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'seedn') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'seedn',cptr%cns%seedn,gcomm_column)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'seedn',cptr%cns%seedn,gcomm_column)
    end if

    ! prod10n
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'prod10n',(/'column'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'prod10n') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'prod10n',cptr%cns%prod10n,gcomm_column)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'prod10n',cptr%cns%prod10n,gcomm_column)
    end if

    ! prod100n
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'prod100n',(/'column'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'prod100n') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'prod100n',cptr%cns%prod100n,gcomm_column)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'prod100n',cptr%cns%prod100n,gcomm_column)
    end if

    ! decomp_cascade_state
    ! the purpose of this is to check to make sure the bgc used matches
    ! what the restart file was generated with.
    ! add info about the SOM decomposition cascade
#ifdef CENTURY_DECOMP
    decomp_cascade_state = 1
#else
    decomp_cascade_state = 0
#endif
    ! add info about the nitrification / denitrification state
#ifdef NITRIF_DENITRIF
    decomp_cascade_state = decomp_cascade_state + 10
#endif

    if (flag == 'define') then
      call clm_addvar(clmvar_integer,ncid,'decomp_cascade_state', &
            long_name='BGC of the model that wrote this restart file:&
            &1s column: 0 = CLM-CN cascade, 1 = Century cascade;&
            &10s column: 0 = CLM-CN denitrification, 10 = Century&
            &denitrification',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. &
           .not. clm_check_var(ncid,'decomp_cascade_state') ) then
        !!! assume, for sake of backwards compatibility, that if
        !!! decomp_cascade_state is not in the restart file, then
        !!! the current model state is the same as the prior model state
        restart_file_decomp_cascade_state = decomp_cascade_state
        write(stdout,*) 'Assuming decomp_cascade_state set : ', &
                decomp_cascade_state
      else
        call clm_readvar(ncid,'decomp_cascade_state', &
                restart_file_decomp_cascade_state)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'decomp_cascade_state', decomp_cascade_state)
    end if

    if ( flag == 'read' .and. &
         decomp_cascade_state /= restart_file_decomp_cascade_state ) then
      if ( myid == italk ) then
        write(stderr,*) 'CNRest: ERROR--the decomposition cascade differs&
                & between the current model state and the model that wrote&
                & the restart file. '
        write(stderr,*) 'This means that the model will be horribly out &
                &of equilibrium until after a lengthy spinup. '
        write(stderr,*) 'Stopping here since this is probably an error in&
               & configuring the run. If you really wish to proceed, '
        write(stderr,*) 'then override by setting &
               &override_bgc_restart_mismatch_dump to .true. in the namelist'
        if ( .not. override_bgc_restart_mismatch_dump ) then
          call fatal(__FILE__,__LINE__, &
                ' CNRest: Stopping. Decomposition cascade mismatch error.')
        end if
      end if
    end if

    ! spinup_state
    if (flag == 'define') then
      call clm_addvar(clmvar_integer,ncid,'spinup_state', &
            long_name='Spinup state of the model that wrote this &
            &restart file: 0 = normal model mode, 1 = AD spinup',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'spinup_state') ) then
        !!! assume, for sake of backwards compatibility, that if
        !!! decomp_cascade_state is not in the restart file, then
        !!! the current model state is the same as the prior model state
        restart_file_spinup_state = spinup_state
        write(stdout,*) 'Assuming spinup_state set : ', spinup_state
      else
        call clm_readvar(ncid,'spinup_state',restart_file_spinup_state)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'spinup_state',spinup_state)
    end if

    ! now compare the model and restart file spinup states, and either
    ! take the model into spinup mode or out of it if they are not identical
    ! taking model out of spinup mode requires multiplying each decomposing
    ! pool by the associated AD factor.
    ! putting model into spinup mode requires dividing each decomposing
    ! pool by the associated AD factor.
    if ( flag == 'read' .and. spinup_state /= restart_file_spinup_state ) then
      if ( spinup_state == 0 .and. restart_file_spinup_state == 1 ) then
        if ( myid == italk ) then
          write(stdout,*) ' CNRest: taking SOM pools out of AD spinup mode'
        end if
        exit_spinup = .true.
      else if ( spinup_state == 1 .and. restart_file_spinup_state == 0 ) then
        if ( myid == italk ) then
          write(stdout,*) ' CNRest: taking SOM pools into AD spinup mode'
        end if
        enter_spinup = .true.
      else
        call fatal(__FILE__,__LINE__,&
            ' CNRest: error in entering/exiting spinup.  &
            &spinup_state != restart_file_spinup_state, but &
            &do not know what to do')
      end if
      if ( ktau >= ntsrf ) then
        call fatal(__FILE__,__LINE__,&
            ' CNRest: error in entering/exiting spinup. this should &
            &occur only when nstep = 1 ')
      end if
      do k = 1 , ndecomp_pools
        if ( exit_spinup ) then
          m = decomp_cascade_con%spinup_factor(k)
        else if ( enter_spinup ) then
          m = d_one / decomp_cascade_con%spinup_factor(k)
        end if
        do c = begc , endc
          do j = 1 , nlevdecomp
            clm3%g%l%c%ccs%decomp_cpools_vr(c,j,k) = &
                    clm3%g%l%c%ccs%decomp_cpools_vr(c,j,k) * m
            if ( use_c13 ) then
              clm3%g%l%c%cc13s%decomp_cpools_vr(c,j,k) = &
                      clm3%g%l%c%cc13s%decomp_cpools_vr(c,j,k) * m
            end if
            if ( use_c14 ) then
              clm3%g%l%c%cc14s%decomp_cpools_vr(c,j,k) = &
                      clm3%g%l%c%cc14s%decomp_cpools_vr(c,j,k) * m
            end if
            clm3%g%l%c%cns%decomp_npools_vr(c,j,k) = &
                    clm3%g%l%c%cns%decomp_npools_vr(c,j,k) * m
          end do
        end do
      end do
    end if

    if ( ktau == 0 ) then
      do i = begp , endp
        if (pftcon%c3psn(clm3%g%l%c%p%itype(i)) == 1.D0) then
          pcisos%grainc(i) = pcbulks%grainc(i) * c3_r2
          pcisos%grainc_storage(i) = pcbulks%grainc_storage(i) * c3_r2
          pcisos%grainc_xfer(i) = pcbulks%grainc_xfer(i) * c3_r2
          pcisos%dispvegc(i) = pcbulks%dispvegc(i) * c3_r2
          pcisos%storvegc(i) = pcbulks%storvegc(i) * c3_r2
          pcisos%totvegc(i) = pcbulks%totvegc(i) * c3_r2
          pcisos%totpftc(i) = pcbulks%totpftc(i) * c3_r2
          pcisos%woodc(i) = pcbulks%woodc(i) * c3_r2
        else
          pcisos%grainc(i) = pcbulks%grainc(i) * c4_r2
          pcisos%grainc_storage(i) = pcbulks%grainc_storage(i) * c4_r2
          pcisos%grainc_xfer(i) = pcbulks%grainc_xfer(i) * c4_r2
          pcisos%dispvegc(i) = pcbulks%dispvegc(i) * c4_r2
          pcisos%storvegc(i) = pcbulks%storvegc(i) * c4_r2
          pcisos%totvegc(i) = pcbulks%totvegc(i) * c4_r2
          pcisos%totpftc(i) = pcbulks%totpftc(i) * c4_r2
          pcisos%woodc(i) = pcbulks%woodc(i) * c4_r2
        end if
      end do
    end if

#if (defined CNDV)
    ! pft type dgvm physical state - crownarea
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'CROWNAREA',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'CROWNAREA') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'CROWNAREA',pptr%pdgvs%crownarea,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'CROWNAREA',pptr%pdgvs%crownarea,gcomm_pft)
    end if

    ! tempsum_litfall
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'tempsum_litfall',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'tempsum_litfall') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'tempsum_litfall', &
                pptr%pepv%tempsum_litfall,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'tempsum_litfall', &
              pptr%pepv%tempsum_litfall,gcomm_pft)
    end if

    ! annsum_litfall
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'annsum_litfall',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'annsum_litfall') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'annsum_litfall', &
                pptr%pepv%annsum_litfall,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'annsum_litfall', &
              pptr%pepv%annsum_litfall,gcomm_pft)
    end if

    ! nind
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'nind',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'nind') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'nind',pptr%pdgvs%nind,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'nind',pptr%pdgvs%nind,gcomm_pft)
    end if

    ! fpcgrid
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'fpcgrid',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'fpcgrid') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'fpcgrid',pptr%pdgvs%fpcgrid,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'fpcgrid',pptr%pdgvs%fpcgrid,gcomm_pft)
    end if

    ! fpcgridold
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'fpcgridold',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'fpcgridold') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'fpcgridold',pptr%pdgvs%fpcgridold,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'fpcgridold',pptr%pdgvs%fpcgridold,gcomm_pft)
    end if

    ! gridcell type dgvm physical state - tmomin20
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'TMOMIN20',(/'gridcell'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'TMOMIN20') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'TMOMIN20',gptr%gdgvs%tmomin20,gcomm_gridcell)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'TMOMIN20',gptr%gdgvs%tmomin20,gcomm_gridcell)
    end if

    ! gridcell type dgvm physical state - agdd20
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'AGDD20',(/'gridcell'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'AGDD20') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'AGDD20',gptr%gdgvs%agdd20,gcomm_gridcell)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'AGDD20',gptr%gdgvs%agdd20,gcomm_gridcell)
    end if

    ! pft type dgvm physical state - t_mo_min
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'T_MO_MIN',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'T_MO_MIN') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'T_MO_MIN',pptr%pdgvs%t_mo_min,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'T_MO_MIN',pptr%pdgvs%t_mo_min,gcomm_pft)
    end if

    ! present
    if (flag == 'define') then
      call clm_addvar(clmvar_integer,ncid,'present',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'present') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        allocate (iptemp(begp:endp), stat=ier)
        call clm_readvar(ncid,'present',iptemp,gcomm_pft)
        do p = begp , endp
          pptr%pdgvs%present(p) = .false.
          if (iptemp(p) == 1) pptr%pdgvs%present(p) = .true.
        end do
        deallocate(iptemp)
      end if
    else if (flag == 'write') then
      allocate (iptemp(begp:endp), stat=ier)
      do p = begp , endp
        iptemp(p) = 0
        if (pptr%pdgvs%present(p)) iptemp(p) = 1
      end do
      call clm_writevar(ncid,'present',iptemp,gcomm_pft)
      deallocate(iptemp)
    end if

    ! leafcmax
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'leafcmax',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'leafcmax') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'leafcmax',pptr%pcs%leafcmax,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'leafcmax',pptr%pcs%leafcmax,gcomm_pft)
    end if

    ! heatstress
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'heatstress',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'heatstress') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'heatstress',pptr%pdgvs%heatstress,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'heatstress',pptr%pdgvs%heatstress,gcomm_pft)
    end if

    ! greffic
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'greffic',(/'pft'/), &
            long_name='',units='')
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'greffic') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'greffic',pptr%pdgvs%greffic,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'greffic',pptr%pdgvs%greffic,gcomm_pft)
    end if

#endif

  end subroutine CNRest
  !
  ! Read/write CN restart data, for vertical decomp grid that can be set
  ! to have length = 1 or nlevgrnd
  !
  subroutine cnrest_addfld_decomp(ncid,varname,longname,units,flag,data_rl, &
                                  lstop)
    implicit none
    type(clm_filetype)  :: ncid   ! netcdf id
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: longname
    character(len=*), intent(in) :: units
    character(len=*), intent(in) :: flag   !'read' or 'write'
    real(rk8), optional, pointer :: data_rl(:,:)
    logical , intent(inout) , optional :: lstop
    ! true => variable is on initial dataset (read only)
    real(rk8), pointer :: ptr1d(:)
    character(len=256) :: name_vr

    name_vr = trim(varname)//'_vr'
#ifdef VERTSOILC
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,name_vr,         &
              (/'column ','levgrnd'/), long_name=longname, &
              units=units,fill_value=1, switchdim=.true.)
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,name_vr) ) then
        if ( present(lstop) ) then
          if ( .not. lstop ) then
            call fatal(__FILE__,__LINE__,'clm now stopping')
          else
            lstop = .true.
          end if
        else
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,name_vr,data_rl,gcomm_column, switchdim=.true.)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,name_vr,data_rl,gcomm_column, switchdim=.true.)
    end if
#else
    !! nlevdecomp = 1; so treat as 1D variable
    ptr1d => data_rl(:,1)
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,varname, &
              (/'column'/), long_name=longname,   &
              units=units,fill_value=1)
    else if (flag == 'read') then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,varname) ) then
        if ( present(lstop) ) then
          if ( .not. lstop ) then
            call fatal(__FILE__,__LINE__,'clm now stopping')
          else
            lstop = .true.
          end if
        else
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,varname,ptr1d,gcomm_column)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,varname,ptr1d,gcomm_column)
    end if
#endif
  end subroutine cnrest_addfld_decomp

#endif

end module mod_clm_cnrest
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

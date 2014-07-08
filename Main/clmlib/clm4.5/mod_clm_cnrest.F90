module mod_clm_cnrest

#if (defined CN)
  !
  ! Read/Write to/from CN info to CLM restart file.
  !
  use mod_realkinds
  use mod_dynparam
  use mod_mppparam
  use mod_mpmessage
  use mod_clm_nchelper
  use mod_clm_type
  use mod_clm_decomp
  use mod_clm_atmlnd, only : clm_a2l
  use mod_clm_varpar, only : numrad , ndecomp_pools , nlevdecomp
  use mod_clm_decomp , only : get_proc_bounds
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
    integer :: c , p , j , k , i , l ! indices
    integer :: begp , endp   ! per-proc beginning and ending pft indices
    integer :: begc , endc   ! per-proc beginning and ending column indices
    integer :: begl , endl   ! per-proc beginning and ending landunit indices
    integer :: begg , endg   ! per-proc gridcell ending gridcell indices
    real(rk8):: m            ! multiplier for the exit_spinup code
    logical :: readvar       ! determine if variable is on initial file
    character(len=128) :: varname         ! temporary
    type(gridcell_type), pointer :: gptr  ! pointer to gridcell derived subtype
    type(landunit_type), pointer :: lptr  ! pointer to landunit derived subtype
    type(column_type)  , pointer :: cptr  ! pointer to column derived subtype
    type(pft_type)     , pointer :: pptr  ! pointer to pft derived subtype
    integer , pointer :: iptemp(:) ! pointer to memory to be allocated
    integer :: ier                 ! error status
    !temporary arrays for slicing larger arrays
    real(rk8) , pointer :: ptr1d(:), ptr2d(:,:)
    ! spinup state as read from restart file, for determining whether
    ! to enter or exit spinup mode.
    integer  :: restart_file_spinup_state
    logical :: exit_spinup = .false.
    logical :: enter_spinup = .false.
    ! flags for comparing the model and restart decomposition cascades
    integer :: decomp_cascade_state, restart_file_decomp_cascade_state

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
              if (pptr%pcs%leafc(i) .ne. spval .and. &
                  pptr%pcs%leafc(i) .ne. nan ) then
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
              if (pptr%pcs%leafc_storage(i) .ne. spval .and. &
                  pptr%pcs%leafc_storage(i) .ne. nan ) then
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
              if (pptr%pcs%leafc_xfer(i) .ne. spval .and. &
                  pptr%pcs%leafc_xfer(i) .ne. nan ) then
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
              if (pptr%pcs%frootc(i) .ne. spval .and. &
                  pptr%pcs%frootc(i) .ne. nan ) then
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
              if (pptr%pcs%frootc_storage(i) .ne. spval .and. &
                  pptr%pcs%frootc_storage(i) .ne. nan ) then
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
              if (pptr%pcs%frootc_xfer(i) .ne. spval .and. &
                  pptr%pcs%frootc_xfer(i) .ne. nan ) then
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
              if (pptr%pcs%livestemc(i) .ne. spval .and. &
                  pptr%pcs%livestemc(i) .ne. nan ) then
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
              if (pptr%pcs%livestemc_storage(i) .ne. spval .and. &
                  pptr%pcs%livestemc_storage(i) .ne. nan ) then
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
              if (pptr%pcs%livestemc_xfer(i) .ne. spval .and. &
                  pptr%pcs%livestemc_xfer(i) .ne. nan ) then
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
              if (pptr%pcs%deadstemc(i) .ne. spval .and. &
                  pptr%pcs%deadstemc(i) .ne. nan ) then
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
              if (pptr%pcs%deadstemc_storage(i) .ne. spval .and. &
                  pptr%pcs%deadstemc_storage(i) .ne. nan ) then
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
              if (pptr%pcs%deadstemc_xfer(i) .ne. spval .and. &
                  pptr%pcs%deadstemc_xfer(i) .ne. nan ) then
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
              if (pptr%pcs%livecrootc(i) .ne. spval .and. &
                  pptr%pcs%livecrootc(i) .ne. nan ) then
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
              if (pptr%pcs%livecrootc_storage(i) .ne. spval .and. &
                  pptr%pcs%livecrootc_storage(i) .ne. nan ) then
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
              if (pptr%pcs%livecrootc_xfer(i) .ne. spval .and. &
                  pptr%pcs%livecrootc_xfer(i) .ne. nan ) then
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
           call ncd_defvar(ncid=ncid, varname='deadcrootc_14', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='deadcrootc_14', data=pptr%pc14s%deadcrootc, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
           if (flag=='read' .and. .not. readvar) then
              write(stdout,*) 'initializing pptr%pc14s%deadcrootc with atmospheric c14 value'
              do i = begp,endp
                 if (pptr%pcs%deadcrootc(i) .ne. spval .and. pptr%pcs%deadcrootc(i) .ne. nan ) then
                    pptr%pc14s%deadcrootc(i) = pptr%pcs%deadcrootc(i) * c14ratio
                 end if
              end do
              if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
           end if
        end if

        ! deadcrootc_storage
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='deadcrootc_storage_14', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='deadcrootc_storage_14', data=pptr%pc14s%deadcrootc_storage, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
           if (flag=='read' .and. .not. readvar) then
              write(stdout,*) 'initializing pptr%pc14s%deadcrootc_storage with atmospheric c14 value'
              do i = begp,endp
                 if (pptr%pcs%deadcrootc_storage(i) .ne. spval .and. pptr%pcs%deadcrootc_storage(i) .ne. nan ) then
                    pptr%pc14s%deadcrootc_storage(i) = pptr%pcs%deadcrootc_storage(i) * c14ratio
                 end if
              end do
              if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
           end if
        end if

        ! deadcrootc_xfer
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='deadcrootc_xfer_14', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='deadcrootc_xfer_14', data=pptr%pc14s%deadcrootc_xfer, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
           if (flag=='read' .and. .not. readvar) then
              write(stdout,*) 'initializing pptr%pc14s%deadcrootc_xfer with atmospheric c14 value'
              do i = begp,endp
                 if (pptr%pcs%deadcrootc_xfer(i) .ne. spval .and. pptr%pcs%deadcrootc_xfer(i) .ne. nan ) then
                    pptr%pc14s%deadcrootc_xfer(i) = pptr%pcs%deadcrootc_xfer(i) * c14ratio
                 end if
              end do
              if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
           end if
        end if

        ! gresp_storage
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='gresp_storage_14', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='gresp_storage_14', data=pptr%pc14s%gresp_storage, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
           if (flag=='read' .and. .not. readvar) then
              write(stdout,*) 'initializing pptr%pc14s%gresp_storage with atmospheric c14 value'
              do i = begp,endp
                 if (pptr%pcs%gresp_storage(i) .ne. spval .and. pptr%pcs%gresp_storage(i) .ne. nan ) then
                    pptr%pc14s%gresp_storage(i) = pptr%pcs%gresp_storage(i) * c14ratio
                 end if
              end do
              if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
           end if
        end if

        ! gresp_xfer
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='gresp_xfer_14', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='gresp_xfer_14', data=pptr%pc14s%gresp_xfer, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
           if (flag=='read' .and. .not. readvar) then
              write(stdout,*) 'initializing pptr%pc14s%gresp_xfer with atmospheric c14 value'
              do i = begp,endp
                 if (pptr%pcs%gresp_xfer(i) .ne. spval .and. pptr%pcs%gresp_xfer(i) .ne. nan ) then
                    pptr%pc14s%gresp_xfer(i) = pptr%pcs%gresp_xfer(i) * c14ratio
                 end if
              end do
              if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
           end if
        end if

        ! cpool
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='cpool_14', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='cpool_14', data=pptr%pc14s%cpool, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
           if (flag=='read' .and. .not. readvar) then
              write(stdout,*) 'initializing pptr%pc14s%cpool with atmospheric c14 value'
              do i = begp,endp
                 if (pptr%pcs%cpool(i) .ne. spval .and. pptr%pcs%cpool(i) .ne. nan ) then
                    pptr%pc14s%cpool(i) = pptr%pcs%cpool(i) * c14ratio
                 end if
              end do
              if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
           end if
        end if

        ! xsmrpool
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='xsmrpool_14', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='xsmrpool_14', data=pptr%pc14s%xsmrpool, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
           if (flag=='read' .and. .not. readvar) then
              write(stdout,*) 'initializing pptr%pc14s%xsmrpool with atmospheric c14 value'
              do i = begp,endp
                 if (pptr%pcs%xsmrpool(i) .ne. spval .and. pptr%pcs%xsmrpool(i) .ne. nan ) then
                    pptr%pc14s%xsmrpool(i) = pptr%pcs%xsmrpool(i) * c14ratio
                 end if
              end do
              if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
           end if
        end if

        ! pft_ctrunc
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='pft_ctrunc_14', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='pft_ctrunc_14', data=pptr%pc14s%pft_ctrunc, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
           if (flag=='read' .and. .not. readvar) then
              write(stdout,*) 'initializing pptr%pc14s%pft_ctrunc with atmospheric c14 value'
              do i = begp,endp
                 if (pptr%pcs%pft_ctrunc(i) .ne. spval .and. pptr%pcs%pft_ctrunc(i) .ne. nan ) then
                    pptr%pc14s%pft_ctrunc(i) = pptr%pcs%pft_ctrunc(i) * c14ratio
                 end if
              end do
              if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
           end if
        end if

        ! totvegc
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='totvegc_14', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='totvegc_14', data=pptr%pc14s%totvegc, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
           if (flag=='read' .and. .not. readvar) then
              write(stdout,*) 'initializing pptr%pc14s%totvegc with atmospheric c14 value'
              do i = begp,endp
                 if (pptr%pcs%totvegc(i) .ne. spval .and. pptr%pcs%totvegc(i) .ne. nan ) then
                    pptr%pc14s%totvegc(i) = pptr%pcs%totvegc(i) * c14ratio
                 end if
              end do
              if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
           end if
        end if

        ! rc14_atm
        if (flag == 'define') then
           call ncd_defvar(ncid=ncid, varname='rc14_atm', xtype=ncd_double,  &
                dim1name='pft',long_name='',units='')
        else if (flag == 'read' .or. flag == 'write') then
           call ncd_io(varname='rc14_atm', data=pptr%pepv%rc14_atm, &
                dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
           if (flag=='read' .and. .not. readvar) then
              write(stdout,*) 'initializing cptr%cc14s%rc14_atm with atmospheric c14 value'
              do i = begp,endp
                 pptr%pepv%rc14_atm(i) = c14ratio
              end do
              if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
           end if
        end if

     end if

    !--------------------------------
    ! pft nitrogen state variables
    !--------------------------------

    ! leafn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafn', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='leafn', data=pptr%pns%leafn, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! leafn_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafn_storage', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='leafn_storage', data=pptr%pns%leafn_storage, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! leafn_xfer
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafn_xfer', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='leafn_xfer', data=pptr%pns%leafn_xfer, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! frootn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frootn', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='frootn', data=pptr%pns%frootn, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! frootn_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frootn_storage', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='frootn_storage', data=pptr%pns%frootn_storage, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! frootn_xfer
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='frootn_xfer', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='frootn_xfer', data=pptr%pns%frootn_xfer, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! livestemn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livestemn', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='livestemn', data=pptr%pns%livestemn, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! livestemn_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livestemn_storage', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='livestemn_storage', data=pptr%pns%livestemn_storage, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! livestemn_xfer
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livestemn_xfer', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='livestemn_xfer', data=pptr%pns%livestemn_xfer, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! deadstemn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadstemn', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='deadstemn', data=pptr%pns%deadstemn, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    !deadstemn_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadstemn_storage', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='deadstemn_storage', data=pptr%pns%deadstemn_storage, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    !deadstemn_xfer
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadstemn_xfer', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='deadstemn_xfer', data=pptr%pns%deadstemn_xfer, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! livecrootn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livecrootn', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='livecrootn', data=pptr%pns%livecrootn, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! livecrootn_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livecrootn_storage', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='livecrootn_storage', data=pptr%pns%livecrootn_storage, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    !livecrootn_xfer
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='livecrootn_xfer', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='livecrootn_xfer', data=pptr%pns%livecrootn_xfer, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! deadcrootn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadcrootn', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='deadcrootn', data=pptr%pns%deadcrootn, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! deadcrootn_storage
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadcrootn_storage', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='deadcrootn_storage', data=pptr%pns%deadcrootn_storage, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! deadcrootn_xfer
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='deadcrootn_xfer', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='deadcrootn_xfer', data=pptr%pns%deadcrootn_xfer, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    !retransn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='retransn', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='retransn', data=pptr%pns%retransn, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! npool
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='npool', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='npool', data=pptr%pns%npool, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! pft_ntrunc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='pft_ntrunc', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='pft_ntrunc', data=pptr%pns%pft_ntrunc, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    !--------------------------------
    ! column physical state variables
    !--------------------------------

    ! decl
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='decl', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='decl', data=cptr%cps%decl, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! fpi
    call cnrest_addfld_decomp(ncid=ncid, varname='fpi',                        &
                              longname='fraction of potential immobilization', &
                              units='unitless', flag=flag,                     &
                              data_rl=cptr%cps%fpi_vr, readvar=readvar)

#ifdef VERTSOILC
    ! som_adv_coef
    call cnrest_addfld_decomp(ncid=ncid, varname='som_adv_coef',   &
                      longname='SOM advective flux', units='m/s',  &
                      flag=flag, fill_value=1,                 &
                      data_rl=cptr%cps%som_adv_coef, readvar=readvar)

    ! som_diffus_coef
    call cnrest_addfld_decomp(ncid=ncid, varname='som_diffus_coef',          &
                      longname='SOM diffusivity due to bio/cryo-turbation',  &
                      units='m^2/s',  flag=flag, fill_value=1,           &
                      data_rl=cptr%cps%som_diffus_coef, readvar=readvar)
#endif

! #ifdef NITRIF_DENITRIF
!     ! tmean_monthly_max
!     call cnrest_addfld_decomp(ncid=ncid, varname='tmean_monthly_max', longname='', units='', flag=flag, data_rl=cptr%cps%tmean_monthly_max, readvar=readvar)

!     ! tmean_monthly
!     call cnrest_addfld_decomp(ncid=ncid, varname='tmean_monthly', longname='', units='', flag=flag, data_rl=cptr%cps%tmean_monthly, readvar=readvar)
! #end if
    ! fpg
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='fpg', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='fpg', data=cptr%cps%fpg, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! annsum_counter
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='annsum_counter', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='annsum_counter', data=cptr%cps%annsum_counter, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! cannsum_npp
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cannsum_npp', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='cannsum_npp', data=cptr%cps%cannsum_npp, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! col_lag_npp
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='col_lag_npp', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='col_lag_npp', data=cptr%cps%col_lag_npp, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! cannavg_t2m
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='cannavg_t2m', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='cannavg_t2m', data=cptr%cps%cannavg_t2m, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

  ! for fire model changed by F. Li and S. Levis
  !  burndate
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='burndate', xtype=ncd_int,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='burndate', data=pptr%pps%burndate, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

   !lfc
     if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='lfc', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='lfc', data=cptr%cps%lfc, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    !wf
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='wf', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='wf', data=cptr%cps%wf, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if


    !btran2
      if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='btran2', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='btran2', data=pptr%pps%btran2, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    !--------------------------------
    ! column carbon state variables
    !--------------------------------

    do k = 1, ndecomp_pools
       ptr2d => cptr%ccs%decomp_cpools_vr(:,:,k)
       varname=trim(decomp_cascade_con%decomp_pool_name_restart(k))//'c'
       call cnrest_addfld_decomp(ncid=ncid, varname=varname, longname='', units='', flag=flag, data_rl=ptr2d, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          call fatal(__FILE__,__LINE__, &
            'ERROR:: '//trim(varname)//' is required on an initialization dataset' )
       end if
    end do

    ! col_ctrunc
    varname = 'col_ctrunc'
    call cnrest_addfld_decomp(ncid=ncid, varname=varname, longname='', units='', flag=flag, data_rl=cptr%ccs%col_ctrunc_vr, readvar=readvar)
    if (flag=='read' .and. .not. readvar) then
       call fatal(__FILE__,__LINE__, &
         'ERROR:: '//trim(varname)//' is required on an initialization dataset' )
    end if

    ! ! nfixation_prof
    ! call cnrest_addfld_decomp(ncid=ncid, varname='nfixation_prof', longname='', units='', flag=flag, data_rl=cptr%cps%nfixation_prof, readvar=readvar)

    ! ! ndep_prof
    ! call cnrest_addfld_decomp(ncid=ncid, varname='ndep_prof', longname='', units='', flag=flag, data_rl=cptr%cps%ndep_prof, readvar=readvar)

    ! altmax
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='altmax', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='altmax', data=cptr%cps%altmax, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! altmax_lastyear
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='altmax_lastyear', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='altmax_lastyear', data=cptr%cps%altmax_lastyear, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! altmax_indx
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='altmax_indx', xtype=ncd_int,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='altmax_indx', data=cptr%cps%altmax_indx, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! altmax_lastyear_indx
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='altmax_lastyear_indx', xtype=ncd_int,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='altmax_lastyear_indx', data=cptr%cps%altmax_lastyear_indx, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! seedc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='seedc', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='seedc', data=cptr%ccs%seedc, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! totlitc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='totlitc', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='totlitc', data=cptr%ccs%totlitc, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! totcolc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='totcolc', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='totcolc', data=cptr%ccs%totcolc, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! prod10c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='prod10c', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='prod10c', data=cptr%ccs%prod10c, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! prod100c
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='prod100c', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='prod100c', data=cptr%ccs%prod100c, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! totsomc
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='totsomc', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='totsomc', data=cptr%ccs%totsomc, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if


    if ( use_c13 ) then
       !--------------------------------
       ! C13 column carbon state variables
       !--------------------------------

       do k = 1, ndecomp_pools
          ptr2d => cptr%cc13s%decomp_cpools_vr(:,:,k)
          varname=trim(decomp_cascade_con%decomp_pool_name_restart(k))//'c_13'
          call cnrest_addfld_decomp(ncid=ncid, varname=varname, longname='', units='', flag=flag, data_rl=ptr2d, readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             write(stdout,*) 'initializing cptr%cc13s%decomp_cpools_vr with atmospheric c13 value for: '//varname
             do i = begc,endc
                do j = 1, nlevdecomp
                   if (cptr%ccs%decomp_cpools_vr(i,j,k) .ne. spval .and. cptr%ccs%decomp_cpools_vr(i,j,k) .ne. nan ) then
                      cptr%cc13s%decomp_cpools_vr(i,j,k) = cptr%ccs%decomp_cpools_vr(i,j,k) * c3_r2
                   end if
                end do
             end do
          end if
       end do

       ! seedc
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='seedc_13', xtype=ncd_double,  &
               dim1name='column',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='seedc_13', data=cptr%cc13s%seedc, &
               dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             do i = begc,endc
                if (cptr%ccs%seedc(i) .ne. spval .and. cptr%ccs%seedc(i) .ne. nan ) then
                   cptr%cc13s%seedc(i) = cptr%ccs%seedc(i) * c3_r2
                end if
             end do
             if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
       end if

       ! col_ctrunc_13
       call cnrest_addfld_decomp(ncid=ncid, varname='col_ctrunc_13_vr', longname='', units='', flag=flag, data_rl=cptr%cc13s%col_ctrunc_vr, readvar=readvar)

       ! ! col_ctrunc
       ! if (flag == 'define') then
       !    call ncd_defvar(ncid=ncid, varname='col_ctrunc_13', xtype=ncd_double,  &
       !         dim1name='column',long_name='',units='')
       ! else if (flag == 'read' .or. flag == 'write') then
       !    call ncd_io(varname='col_ctrunc_13', data=cptr%cc13s%col_ctrunc, &
       !         dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
       !    if (flag=='read' .and. .not. readvar) then
       !       if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now
       !       stopping')
       !    end if
       ! end if

       ! totlitc
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='totlitc_13', xtype=ncd_double,  &
               dim1name='column',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='totlitc_13', data=cptr%cc13s%totlitc, &
               dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             do i = begc,endc
                if (cptr%ccs%totlitc(i) .ne. spval .and. cptr%ccs%totlitc(i) .ne. nan ) then
                   cptr%cc13s%totlitc(i) = cptr%ccs%totlitc(i) * c3_r2
                end if
             end do
             if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
       end if

       ! totcolc
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='totcolc_13', xtype=ncd_double,  &
               dim1name='column',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='totcolc_13', data=cptr%cc13s%totcolc, &
               dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             do i = begc,endc
                if (cptr%ccs%totcolc(i) .ne. spval .and. cptr%ccs%totcolc(i) .ne. nan ) then
                   cptr%cc13s%totcolc(i) = cptr%ccs%totcolc(i) * c3_r2
                end if
             end do
             if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
       end if

       ! prod10c
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='prod10c_13', xtype=ncd_double,  &
               dim1name='column',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='prod10c_13', data=cptr%cc13s%prod10c, &
               dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             do i = begc,endc
                if (cptr%ccs%prod10c(i) .ne. spval .and. cptr%ccs%prod10c(i) .ne. nan ) then
                   cptr%cc13s%prod10c(i) = cptr%ccs%prod10c(i) * c3_r2
                end if
             end do
             if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
       end if

       ! prod100c
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='prod100c_13', xtype=ncd_double,  &
               dim1name='column',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='prod100c_13', data=cptr%cc13s%prod100c, &
               dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             do i = begc,endc
                if (cptr%ccs%prod100c(i) .ne. spval .and. cptr%ccs%prod100c(i) .ne. nan ) then
                   cptr%cc13s%prod100c(i) = cptr%ccs%prod100c(i) * c3_r2
                end if
             end do
             if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
       end if
    end if

    if ( use_c14 ) then
       !--------------------------------
       ! C14 column carbon state variables
       !--------------------------------

       do k = 1, ndecomp_pools
          ptr2d => cptr%cc14s%decomp_cpools_vr(:,:,k)
          varname=trim(decomp_cascade_con%decomp_pool_name_restart(k))//'c_14'
          call cnrest_addfld_decomp(ncid=ncid, varname=varname, longname='', units='', flag=flag, data_rl=ptr2d, readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             write(stdout,*) 'initializing cptr%cc14s%decomp_cpools_vr with atmospheric c14 value for: '//varname
             do i = begc,endc
                do j = 1, nlevdecomp
                   if (cptr%ccs%decomp_cpools_vr(i,j,k) .ne. spval .and. cptr%ccs%decomp_cpools_vr(i,j,k) .ne. nan ) then
                      cptr%cc14s%decomp_cpools_vr(i,j,k) = cptr%ccs%decomp_cpools_vr(i,j,k) * c14ratio
                   end if
                end do
             end do
          end if
       end do

       ! seedc
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='seedc_14', xtype=ncd_double,  &
               dim1name='column',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='seedc_14', data=cptr%cc14s%seedc, &
               dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             write(stdout,*) 'initializing cptr%cc14s%seedc with atmospheric c14 value'
             do i = begc,endc
                if (cptr%ccs%seedc(i) .ne. spval .and. cptr%ccs%seedc(i) .ne. nan ) then
                   cptr%cc14s%seedc(i) = cptr%ccs%seedc(i) * c14ratio
                end if
             end do
             if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
       end if

       ! col_ctrunc_c14
       call cnrest_addfld_decomp(ncid=ncid, varname='col_ctrunc_14_vr', longname='', units='', flag=flag, data_rl=cptr%cc14s%col_ctrunc_vr, readvar=readvar)

       ! ! col_ctrunc
       ! if (flag == 'define') then
       !    call ncd_defvar(ncid=ncid, varname='col_ctrunc_14', xtype=ncd_double,  &
       !         dim1name='column',long_name='',units='')
       ! else if (flag == 'read' .or. flag == 'write') then
       !    call ncd_io(varname='col_ctrunc_14', data=cptr%cc14s%col_ctrunc, &
       !         dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
       !    if (flag=='read' .and. .not. readvar) then
       !       write(stdout,*) 'initializing cptr%cc14s%col_ctrunc with atmospheric c14 value'
       !       do i = begc,endc
       !          if (cptr%ccs%col_ctrunc(i) .ne. spval .and. cptr%ccs%col_ctrunc(i) .ne. nan ) then
       !             cptr%cc14s%col_ctrunc(i) = cptr%ccs%col_ctrunc(i) * c14ratio
       !          end if
       !       end do
       !       if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now
       !       stopping')
       !    end if
       ! end if

       ! totlitc
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='totlitc_14', xtype=ncd_double,  &
               dim1name='column',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='totlitc_14', data=cptr%cc14s%totlitc, &
               dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             write(stdout,*) 'initializing cptr%cc14s%totlitc with atmospheric c14 value'
             do i = begc,endc
                if (cptr%ccs%totlitc(i) .ne. spval .and. cptr%ccs%totlitc(i) .ne. nan ) then
                   cptr%cc14s%totlitc(i) = cptr%ccs%totlitc(i) * c14ratio
                end if
             end do
             if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
       end if

       ! totcolc
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='totcolc_14', xtype=ncd_double,  &
               dim1name='column',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='totcolc_14', data=cptr%cc14s%totcolc, &
               dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             write(stdout,*) 'initializing cptr%cc14s%totcolc with atmospheric c14 value'
             do i = begc,endc
                if (cptr%ccs%totcolc(i) .ne. spval .and. cptr%ccs%totcolc(i) .ne. nan ) then
                   cptr%cc14s%totcolc(i) = cptr%ccs%totcolc(i) * c14ratio
                end if
             end do
             if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
       end if

       ! prod10c
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='prod10c_14', xtype=ncd_double,  &
               dim1name='column',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='prod10c_14', data=cptr%cc14s%prod10c, &
               dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             write(stdout,*) 'initializing cptr%cc14s%prod10c with atmospheric c14 value'
             do i = begc,endc
                if (cptr%ccs%prod10c(i) .ne. spval .and. cptr%ccs%prod10c(i) .ne. nan ) then
                   cptr%cc14s%prod10c(i) = cptr%ccs%prod10c(i) * c14ratio
                end if
             end do
             if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
       end if

       ! prod100c
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname='prod100c_14', xtype=ncd_double,  &
               dim1name='column',long_name='',units='')
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname='prod100c_14', data=cptr%cc14s%prod100c, &
               dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             write(stdout,*) 'initializing cptr%cc14s%prod100c with atmospheric c14 value'
             do i = begc,endc
                if (cptr%ccs%prod100c(i) .ne. spval .and. cptr%ccs%prod100c(i) .ne. nan ) then
                   cptr%cc14s%prod100c(i) = cptr%ccs%prod100c(i) * c14ratio
                end if
             end do
             if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
       end if
    end if

    !--------------------------------
    ! column nitrogen state variables
    !--------------------------------

    ! sminn
    varname='sminn'
    call cnrest_addfld_decomp(ncid=ncid, varname=varname, longname='', units='', flag=flag, data_rl=cptr%cns%sminn_vr, readvar=readvar)
    if (flag=='read' .and. .not. readvar) then
       call fatal(__FILE__,__LINE__, &
         'ERROR:: '//trim(varname)//' is required on an initialization dataset' )
    end if

    ! decomposing N pools
    do k = 1, ndecomp_pools
       ptr2d => cptr%cns%decomp_npools_vr(:,:,k)
       varname=trim(decomp_cascade_con%decomp_pool_name_restart(k))//'n'
       call cnrest_addfld_decomp(ncid=ncid, varname=varname, longname='', units='', flag=flag, data_rl=ptr2d, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          call fatal(__FILE__,__LINE__,&
            'ERROR:: '//trim(varname)//' is required on an initialization dataset' )
       end if
    end do


    ! col_ntrunc
    call cnrest_addfld_decomp(ncid=ncid, varname='col_ntrunc', longname='', units='', flag=flag, data_rl=cptr%cns%col_ntrunc_vr, readvar=readvar)

#ifdef NITRIF_DENITRIF
    ! f_nit_vr
    call cnrest_addfld_decomp(ncid=ncid, varname='f_nit_vr', &
                              longname='soil nitrification flux', &
                              units='gN/m3/s', &
                              flag=flag, data_rl=cptr%cnf%f_nit_vr, readvar=readvar)

    ! pot_f_nit_vr
    call cnrest_addfld_decomp(ncid=ncid, varname='pot_f_nit_vr', &
                              longname='potential soil nitrification flux', &
                              units='gN/m3/s', &
                              flag=flag, data_rl=cptr%cnf%pot_f_nit_vr, readvar=readvar)

    ! smin_no3_vr
    varname = 'smin_no3'
    call cnrest_addfld_decomp(ncid=ncid, varname=varname, longname='', units='', flag=flag, data_rl=cptr%cns%smin_no3_vr, readvar=readvar)
    if (flag=='read' .and. .not. readvar) then
       call fatal(__FILE__,__LINE__,&
         'ERROR:: '//trim(varname)//' is required on an initialization dataset' )
    end if

    ! smin_nh4
    varname = 'smin_nh4'
    call cnrest_addfld_decomp(ncid=ncid, varname=varname, longname='', units='', flag=flag, data_rl=cptr%cns%smin_nh4_vr, readvar=readvar)
    if (flag=='read' .and. .not. readvar) then
       call fatal(__FILE__,__LINE__, &
         'ERROR:: '//trim(varname)//' is required on an initialization dataset' )
    end if

#endif

    ! totcoln
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='totcoln', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='totcoln', data=cptr%cns%totcoln, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! seedn
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='seedn', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='seedn', data=cptr%cns%seedn, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! prod10n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='prod10n', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='prod10n', data=cptr%cns%prod10n, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! prod100n
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='prod100n', xtype=ncd_double,  &
            dim1name='column',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='prod100n', data=cptr%cns%prod100n, &
            dim1name=namec, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
      if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! decomp_cascade_state
    ! the purpose of this is to check to make sure the bgc used matches what the restart file was generated with.
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
       call ncd_defvar(ncid=ncid, varname='decomp_cascade_state', xtype=ncd_int,  &
            long_name='BGC of the model that wrote this restart file: 1s column: 0 = CLM-CN cascade, 1 = Century cascade;' &
            // ' 10s column: 0 = CLM-CN denitrification, 10 = Century denitrification',units='')
    else if (flag == 'read') then
       call ncd_io(varname='decomp_cascade_state', data=restart_file_decomp_cascade_state, &
            ncid=ncid, flag=flag, readvar=readvar)
       if ( .not. readvar) then
      !!! assume, for sake of backwards compatibility, that if decomp_cascade_state is not in the restart file, then the current model state is the same as the prior model state
          restart_file_decomp_cascade_state = decomp_cascade_state
          if ( myid == italk ) then
            write(stderr,*) ' CNRest: WARNING!  Restart file does not contain info on decomp_cascade_state used to generate the restart file.  '
            write(stderr,*) '   Assuming the same as current setting: ', decomp_cascade_state
          end if
       end if
    else if (flag == 'write') then
       call ncd_io(varname='decomp_cascade_state', data=decomp_cascade_state, &
            ncid=ncid, flag=flag, readvar=readvar)
    end if
    if ( flag == 'read' .and. decomp_cascade_state .ne. restart_file_decomp_cascade_state ) then
       if ( myid == italk ) then
           write(stderr,*) 'CNRest: ERROR--the decomposition cascade differs between the current model state and the model that wrote the restart file. '
           write(stderr,*) 'This means that the model will be horribly out of equilibrium until after a lengthy spinup. '
           write(stderr,*) 'Stopping here since this is probably an error in configuring the run. If you really wish to proceed, '
           write(stderr,*) 'then override by setting override_bgc_restart_mismatch_dump to .true. in the namelist'
           if ( .not. override_bgc_restart_mismatch_dump ) then
              call fatal(__FILE__,__LINE__, &
                ' CNRest: Stopping. Decomposition cascade mismatch error.')
           end if
        end if
    end if

    ! spinup_state
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='spinup_state', xtype=ncd_int,  &
            long_name='Spinup state of the model that wrote this restart file: 0 = normal model mode, 1 = AD spinup',units='')
    else if (flag == 'read') then
       call ncd_io(varname='spinup_state', data=restart_file_spinup_state, &
            ncid=ncid, flag=flag, readvar=readvar)
       if ( .not. readvar) then
      !!! assume, for sake of backwards compatibility, that if spinup_state is not in the restart file, then the current model state is the same as the prior model state
          restart_file_spinup_state = spinup_state
          if ( myid == italk ) then
            write(stderr,*) ' CNRest: WARNING!  Restart file does not contain info on spinup state used to generate the restart file. '
            write(stderr,*) '   Assuming the same as current setting: ', spinup_state
          end if
       end if
    else if (flag == 'write') then
       call ncd_io(varname='spinup_state', data=spinup_state, &
            ncid=ncid, flag=flag, readvar=readvar)
    end if

    ! now compare the model and restart file spinup states, and either take the model into spinup mode or out of it if they are not identical
    ! taking model out of spinup mode requires multiplying each decomposing pool by the associated AD factor.
    ! putting model into spinup mode requires dividing each decomposing pool by the associated AD factor.
    if (flag == 'read' .and. spinup_state .ne. restart_file_spinup_state ) then
       if (spinup_state .eq. 0 .and. restart_file_spinup_state .eq. 1 ) then
          if ( myid == italk ) then
            write(stderr,*) ' CNRest: taking SOM pools out of AD spinup mode'
          end if
          exit_spinup = .true.
       else if (spinup_state .eq. 1 .and. restart_file_spinup_state .eq. 0 ) then
          if ( myid == italk ) then
            write(stderr,*) ' CNRest: taking SOM pools into AD spinup mode'
          end if
          enter_spinup = .true.
       else
          call fatal(__FILE__,__LINE__,&
            ' CNRest: error in entering/exiting spinup.  spinup_state != restart_file_spinup_state, but do not know what to do')
       end if
       if (ktau >= ntsrf) then
          call fatal(__FILE__,__LINE__,&
            ' CNRest: error in entering/exiting spinup. this should occur only when nstep = 1 ')
       end if
       do k = 1, ndecomp_pools
          if ( exit_spinup ) then
             m = decomp_cascade_con%spinup_factor(k)
          else if ( enter_spinup ) then
             m = 1. / decomp_cascade_con%spinup_factor(k)
          end if
          do c = begc, endc
             do j = 1, nlevdecomp
                clm3%g%l%c%ccs%decomp_cpools_vr(c,j,k) = clm3%g%l%c%ccs%decomp_cpools_vr(c,j,k) * m

                if ( use_c13 ) then
                   clm3%g%l%c%cc13s%decomp_cpools_vr(c,j,k) = clm3%g%l%c%cc13s%decomp_cpools_vr(c,j,k) * m
                end if

                if ( use_c14 ) then
                   clm3%g%l%c%cc14s%decomp_cpools_vr(c,j,k) = clm3%g%l%c%cc14s%decomp_cpools_vr(c,j,k) * m
                end if

                clm3%g%l%c%cns%decomp_npools_vr(c,j,k) = clm3%g%l%c%cns%decomp_npools_vr(c,j,k) * m
             end do
          end do
       end do
    end if

    if ( ktau == 0 ) then
       do i = begp, endp
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
       call ncd_defvar(ncid=ncid, varname='CROWNAREA', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='CROWNAREA', data=pptr%pdgvs%crownarea, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! tempsum_litfall
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='tempsum_litfall', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='tempsum_litfall', data=pptr%pepv%tempsum_litfall, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! annsum_litfall
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='annsum_litfall', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='annsum_litfall', data=pptr%pepv%annsum_litfall, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! nind
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='nind', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='nind', data=pptr%pdgvs%nind, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! fpcgrid
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='fpcgrid', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='fpcgrid', data=pptr%pdgvs%fpcgrid, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! fpcgridold
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='fpcgridold', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='fpcgridold', data=pptr%pdgvs%fpcgridold, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! gridcell type dgvm physical state - tmomin20
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='TMOMIN20', xtype=ncd_double,  &
            dim1name='gridcell',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='TMOMIN20', data=gptr%gdgvs%tmomin20, &
            dim1name=nameg, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! gridcell type dgvm physical state - agdd20
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='AGDD20', xtype=ncd_double,  &
            dim1name='gridcell',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='AGDD20', data=gptr%gdgvs%agdd20, &
            dim1name=nameg, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! pft type dgvm physical state - t_mo_min
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='T_MO_MIN', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='T_MO_MIN', data=pptr%pdgvs%t_mo_min, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! present
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='present', xtype=ncd_int,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       allocate (iptemp(begp:endp), stat=ier)
       if (ier /= 0) then
          call fatal(__FILE__,__LINE__,'CNrest: allocation error ')
       end if
       if (flag == 'write') then
          do p = begp,endp
             iptemp(p) = 0
             if (pptr%pdgvs%present(p)) iptemp(p) = 1
          end do
       end if
       call ncd_io(varname='present', data=iptemp, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read') then
          if (.not. readvar) then
             if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
          else
             do p = begp,endp
                pptr%pdgvs%present(p) = .false.
                if (iptemp(p) == 1) pptr%pdgvs%present(p) = .true.
             end do
          end if
       end if
       deallocate (iptemp)
    end if

    ! leafcmax
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='leafcmax', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='leafcmax', data=pptr%pcs%leafcmax, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! heatstress
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='heatstress', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='heatstress', data=pptr%pdgvs%heatstress, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if

    ! greffic
    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='greffic', xtype=ncd_double,  &
            dim1name='pft',long_name='',units='')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_io(varname='greffic', data=pptr%pdgvs%greffic, &
            dim1name=namep, ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
       end if
    end if
#endif

  end subroutine CNRest
  !
  ! Read/write CN restart data, for vertical decomp grid that can be set
  ! to have length = 1 or nlevgrnd
  !
  subroutine cnrest_addfld_decomp(ncid,varname,longname,units,flag,data_rl, &
                                  fill_value,readvar)
    use mod_clm_varpar , only : nlevgrnd
    use mod_clm_type , only : namec
    implicit none
    type(clm_filetype)  :: ncid   ! netcdf id
    character(len=*), intent(in) :: varname
    character(len=*), intent(in) :: longname
    character(len=*), intent(in) :: units
    character(len=*), intent(in) :: flag   !'read' or 'write'
    real(rk8), optional :: fill_value
    real(rk8), optional, pointer :: data_rl(:,:)
    ! true => variable is on initial dataset (read only)
    logical, optional, intent(out):: readvar
    real(rk8), pointer :: ptr1d(:)

#ifdef VERTSOILC
    if (flag == 'define') then
      call ncd_defvar(ncid=ncid, varname=trim(varname)//'_vr', &
              xtype=ncd_double,  &
         dim1name='column',dim2name='levgrnd', &
         long_name=longname,units=units)
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname=trim(varname)//'_vr', data=data_rl, &
               dim1name=namec,ncid=ncid, flag=flag, readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
       end if
#else
       !! nlevdecomp = 1; so treat as 1D variable
       ptr1d => data_rl(:,1)
       if (flag == 'define') then
          call ncd_defvar(ncid=ncid, varname=trim(varname), xtype=ncd_double,  &
               dim1name='column',long_name=longname,units=units)
       else if (flag == 'read' .or. flag == 'write') then
          call ncd_io(varname=trim(varname), data=ptr1d, &
               dim1name=namec,ncid=ncid, flag=flag, readvar=readvar)
          if (flag=='read' .and. .not. readvar) then
             if (ktau /= 0) call fatal(__FILE__,__LINE__,'clm now stopping')
          end if
       end if
#endif
  end subroutine cnrest_addfld_decomp

#endif

end module mod_clm_cnrest

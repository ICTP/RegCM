module mod_clm_croprest

#if (defined CN)
!----------------------------------------------------------------------- 
!BOP
!
! !MODULE: CropRestMod
! 
! !DESCRIPTION: 
! Read/Write to/from Crop info to CLM restart file. 
!
! !USES:
  use mod_realkinds
  use mod_mppparam
  use mod_dynparam
  use mod_mpmessage
  use mod_clm_nchelper
!
! !PUBLIC TYPES:
  implicit none
  private
  save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: CropRest        ! Restart prognostic crop model
  public :: CropRestYear    ! Get the number of years crop has spunup
  public :: CropRestIncYear ! Increment the crop spinup years
!
! !REVISION HISTORY:
! Module created by slevis following CNRestMod by Peter Thornton
!

! !PRIVATE DATA MEMBERS:
   integer, parameter :: unset = -999  ! Flag that restart year is not set
   integer :: restyear = unset         ! Restart year from the initial conditions file

!EOP
!----------------------------------------------------------------------- 

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CropRest
!
! !INTERFACE:
  subroutine CropRest ( ncid, flag )
!
! !DESCRIPTION: 
! Read/write Crop restart data
!
! !USES:
    use mod_clm_type
    use mod_clm_atmlnd      , only : clm_a2l
    use mod_clm_varpar      , only : numrad
    use mod_clm_decomp      , only : get_proc_bounds
    use mod_clm_time_manager, only : is_restart
!
! !ARGUMENTS:
    implicit none
    type(clm_filetype) :: ncid             ! netcdf id
    character(len=*), intent(in) :: flag   !'read' or 'write'
!
! !CALLED FROM:
! subroutine restart in module restFileMod
!
! !REVISION HISTORY:
! Author: slevis
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: c,p,j                      ! indices 
    integer :: begp, endp                 ! per-proc beginning and ending pft indices
    integer :: begc, endc                 ! per-proc beginning and ending column indices 
    integer :: begl, endl                 ! per-proc beginning and ending landunit indices
    integer :: begg, endg                 ! per-proc gridcell ending gridcell indices
    real(rk8):: m                          ! multiplier for the exit_spinup code
    character(len=128) :: varname         ! temporary
    type(gridcell_type), pointer :: gptr  ! pointer to gridcell derived subtype
    type(landunit_type), pointer :: lptr  ! pointer to landunit derived subtype
    type(column_type)  , pointer :: cptr  ! pointer to column derived subtype
    type(pft_type)     , pointer :: pptr  ! pointer to pft derived subtype
    integer , pointer :: iptemp(:)        ! pointer to memory to be allocated
    integer :: ier                        ! error status
!-----------------------------------------------------------------------

    ! Prognostic crop restart year
    if (flag == 'define') then
      call clm_addvar(clmvar_integer,ncid,'restyear', &
          long_name='Number of years prognostic crop ran', units="years", &
          missing_value=1,fill_value=1)
    else if (flag == 'read' ) then
      if ( is_restart() .and. .not. clm_check_var(ncid,'restyear') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'restyear',restyear)
        call checkDates( )
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'restyear',restyear)
    end if

    ! Set pointers into derived type

    gptr => clm3%g
    lptr => clm3%g%l
    cptr => clm3%g%l%c
    pptr => clm3%g%l%c%p

    ! Determine necessary subgrid bounds

    call get_proc_bounds(begg, endg, begl, endl, begc, endc, begp, endp)

    !--------------------------------
    ! pft physical state variables 
    !--------------------------------

    ! peaklai
    if (flag == 'define') then
       call clm_addvar(clmvar_integer,ncid,'peaklai',cdims=(/'pft'/), &
            long_name='Flag if at max allowed LAI or not', &
            flag_values=(/0,1/), valid_range=(/0,1/),     &
            flag_meanings=(/'NOT-at-peak', 'AT_peak-LAI' /) )
    else if (flag == 'read' ) then
      if ( is_restart() .and. .not. clm_check_var(ncid,'peaklai') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'peaklai',pptr%pps%peaklai)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'peaklai',pptr%pps%peaklai)
    end if

    ! idop
    if (flag == 'define') then
       call clm_addvar(clmvar_integer,ncid,'idop',cdims=(/'pft'/), &
            long_name='Date of planting',units='jday', &
            valid_range=(/1,366/) )
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'idop') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'idop',pptr%pps%idop)
      end if
    else if ( flag == 'write') then
       call clm_writevar(ncid,'idop',pptr%pps%idop)
    end if

    ! aleaf
    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'aleaf',cdims=(/'pft'/), &
            long_name='leaf allocation coefficient',units='')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'aleaf') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'aleaf',pptr%pps%aleaf)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'aleaf',pptr%pps%aleaf)
    end if

    ! aleafi
    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'aleafi',cdims=(/'pft'/), &
            long_name='Saved leaf allocation coefficient from phase 2', &
            units='')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'aleafi') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'aleafi',pptr%pps%aleafi)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'aleafi',pptr%pps%aleafi)
    end if

    ! astem
    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'astem',cdims=(/'pft'/), &
            long_name='stem allocation coefficient',units='')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'astem') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'astem',pptr%pps%astem)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'astem',pptr%pps%astem)
    end if

    ! astemi
    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'astemi',cdims=(/'pft'/), &
            long_name='Saved stem allocation coefficient from phase 2',&
            units='')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'astemi') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'astemi',pptr%pps%astemi)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'astemi',pptr%pps%astemi)
    end if

    ! htmx 
    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'htmx',cdims=(/'pft'/), &
            long_name='max height attained by a crop during year',&
            units='m')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'htmx') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'htmx',pptr%pps%htmx)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'htmx',pptr%pps%htmx)
    end if

    ! hdidx
    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'hdidx',cdims=(/'pft'/), &
            long_name='cold hardening index',units='')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'hdidx') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'hdidx',pptr%pps%hdidx)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'hdidx',pptr%pps%hdidx)
    end if

    ! vf
    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'vf',cdims=(/'pft'/), &
            long_name='vernalization factor',units='')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'vf') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'vf',pptr%pps%vf)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'vf',pptr%pps%vf)
    end if

    ! cumvd
    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'cumvd',cdims=(/'pft'/), &
            long_name='cumulative vernalization d',units='')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'cumvd') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'cumvd',pptr%pps%cumvd)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'cumvd',pptr%pps%cumvd)
    end if

    ! croplive
    if (flag == 'define') then
       call clm_addvar(clmvar_logical,ncid,'croplive',cdims=(/'pft'/), &
            long_name='Flag that crop is alive, but not harvested')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'croplive') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'croplive',pptr%pps%croplive)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'croplive',pptr%pps%croplive)
    end if

    ! cropplant
    if (flag == 'define') then
       call clm_addvar(clmvar_logical,ncid,'cropplant',cdims=(/'pft'/), &
            long_name='Flag that crop is planted, but not harvested' )
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'cropplant') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'cropplant',pptr%pps%cropplant)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'cropplant',pptr%pps%cropplant)
    end if

    ! harvdate
    if (flag == 'define') then
       call clm_addvar(clmvar_integer,ncid,'harvdate',cdims=(/'pft'/), &
            long_name='harvest date',units='jday',valid_range=(/1,366/) )
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'harvdate') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'harvdate',pptr%pps%harvdate)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'harvdate',pptr%pps%harvdate)
    end if

    ! gdd1020
    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'gdd1020',cdims=(/'pft'/), &
            long_name='20 year average of growing degree-days base 10C from planting', &
            units='ddays')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'gdd1020') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'gdd1020',pptr%pps%gdd1020)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'gdd1020',pptr%pps%gdd1020)
    end if

    ! gdd820
    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'gdd820',cdims=(/'pft'/), &
            long_name='20 year average of growing degree-days base 8C from planting', &
            units='ddays')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'gdd820') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'gdd820',pptr%pps%gdd820)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'gdd820',pptr%pps%gdd820)
    end if

    ! gdd020
    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'gdd020',cdims=(/'pft'/), &
            long_name='20 year average of growing degree-days base 0C from planting', &
            units='ddays')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'gdd020') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'gdd020',pptr%pps%gdd020)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'gdd020',pptr%pps%gdd020)
    end if

    ! gddmaturity
    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'gddmaturity',cdims=(/'pft'/), &
            long_name='Growing degree days needed to harvest',units='ddays')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'gddmaturity') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'gddmaturity',pptr%pps%gddmaturity)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'gddmaturity',pptr%pps%gddmaturity)
    end if

    ! huileaf
    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'huileaf',cdims=(/'pft'/), &
            long_name='heat unit index needed from planting to leaf emergence',units='')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'huileaf') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'huileaf',pptr%pps%huileaf)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'huileaf',pptr%pps%huileaf)
    end if

    ! huigrain
    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'huigrain',cdims=(/'pft'/), &
            long_name='heat unit index needed to reach vegetative maturity', &
            units='')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'huigrain') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'huigrain',pptr%pps%huigrain)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'huigrain',pptr%pps%huigrain)
    end if

    ! grainc
    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'grainc',cdims=(/'pft'/), &
            long_name='grain C',units='gC/m2')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'grainc') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'grainc',pptr%pcs%grainc)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'grainc',pptr%pcs%grainc)
    end if

    ! grainc_storage
    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'grainc_storage',cdims=(/'pft'/), &
            long_name='grain C storage',units='gC/m2')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'grainc_storage') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'grainc_storage',pptr%pcs%grainc_storage)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'grainc_storage',pptr%pcs%grainc_storage)
    end if

    ! grainc_xfer
    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'grainc_xfer',cdims=(/'pft'/), &
            long_name='grain C transfer',units='gC/m2')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'grainc_xfer') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'grainc_xfer',pptr%pcs%grainc_xfer)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'grainc_xfer',pptr%pcs%grainc_xfer)
    end if

    ! grainn
    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'grainn',cdims=(/'pft'/), &
            long_name='grain N',units='gN/m2')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'grainn') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'grainn',pptr%pns%grainn)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'grainn',pptr%pns%grainn)
    end if

    ! grainn_storage
    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'grainn_storage',cdims=(/'pft'/), &
            long_name='grain N storage',units='gN/m2')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'grainn_storage') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'grainn_storage',pptr%pns%grainn_storage)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'grainn_storage',pptr%pns%grainn_storage)
    end if

    ! grainn_xfer
    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'grainn_xfer',cdims=(/'pft'/), &
            long_name='grain N transfer',units='gN/m2')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'grainn_xfer') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'grainn_xfer',pptr%pns%grainn_xfer)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'grainn_xfer',pptr%pns%grainn_xfer)
    end if

    ! grainc_xfer_to_grainc
    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'grainc_xfer_to_grainc',cdims=(/'pft'/), &
            long_name='grain C growth from storage',units='gC/m2/s')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'grainc_xfer_to_grainc') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'grainc_xfer_to_grainc',pptr%pcf%grainc_xfer_to_grainc)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'grainc_xfer_to_grainc',pptr%pcf%grainc_xfer_to_grainc)
    end if

    ! livestemc_to_litter
    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'livestemc_to_litter',cdims=(/'pft'/), &
            long_name='live stem C litterfall',units='gC/m2/s')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'livestemc_to_litter') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'livestemc_to_litter',pptr%pcf%livestemc_to_litter)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'livestemc_to_litter',pptr%pcf%livestemc_to_litter)
    end if

    ! grainc_to_food
    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'grainc_to_food',cdims=(/'pft'/), &
            long_name='grain C to food',units='gC/m2/s')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'grainc_to_food') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'grainc_to_food',pptr%pcf%grainc_to_food)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'grainc_to_food',pptr%pcf%grainc_to_food)
    end if

    ! grainn_xfer_to_grainn
    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'grainn_xfer_to_grainn',cdims=(/'pft'/), &
            long_name='grain N growth from storage',units='gN/m2/s')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'grainn_xfer_to_grainn') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'grainn_xfer_to_grainn',pptr%pnf%grainn_xfer_to_grainn)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'grainn_xfer_to_grainn',pptr%pnf%grainn_xfer_to_grainn)
    end if

    ! livestemn_to_litter
    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'livestemn_to_litter',cdims=(/'pft'/), &
            long_name='livestem N to litter',units='gN/m2/s')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'livestemn_to_litter') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'livestemn_to_litter',pptr%pnf%livestemn_to_litter)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'livestemn_to_litter',pptr%pnf%livestemn_to_litter)
    end if

    ! grainn_to_food
    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'grainn_to_food',cdims=(/'pft'/), &
            long_name='grain N to food',units='gN/m2/s')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'grainn_to_food') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'grainn_to_food',pptr%pnf%grainn_to_food)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'grainn_to_food',pptr%pnf%grainn_to_food)
    end if

    ! cpool_to_grainc
    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'cpool_to_grainc',cdims=(/'pft'/), &
            long_name='allocation to grain C',units='gC/m2/s')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'cpool_to_grainc') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'cpool_to_grainc',pptr%pcf%cpool_to_grainc)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'cpool_to_grainc',pptr%pcf%cpool_to_grainc)
    end if

    ! cpool_to_grainc_storage
    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'cpool_to_grainc_storage',cdims=(/'pft'/), &
            long_name='allocation to grain C storage',units='gC/m2/s')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'cpool_to_grainc_storage') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'cpool_to_grainc_storage',pptr%pcf%cpool_to_grainc_storage)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'cpool_to_grainc_storage',pptr%pcf%cpool_to_grainc_storage)
    end if

    ! npool_to_grainn
    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'npool_to_grainn',cdims=(/'pft'/), &
            long_name='allocation to grain N',units='gN/m2/s')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'npool_to_grainn') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'npool_to_grainn',pptr%pnf%npool_to_grainn)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'npool_to_grainn',pptr%pnf%npool_to_grainn)
    end if

    ! npool_to_grainn_storage
    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'npool_to_grainn_storage',cdims=(/'pft'/), &
            long_name='allocation to grain N storage',units='gN/m2/s')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'npool_to_grainn_storage') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'npool_to_grainn_storage',pptr%pnf%npool_to_grainn_storage)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'npool_to_grainn_storage',pptr%pnf%npool_to_grainn_storage)
    end if

    ! cpool_grain_gr
    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'cpool_grain_gr',cdims=(/'pft'/), &
            long_name='grain growth respiration',units='gC/m2/s')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'cpool_grain_gr') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'cpool_grain_gr',pptr%pcf%cpool_grain_gr)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'cpool_grain_gr',pptr%pcf%cpool_grain_gr)
    end if

    ! cpool_grain_storage_gr
    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'cpool_grain_storage_gr',cdims=(/'pft'/), &
            long_name='grain growth respiration to storage',units='gC/m2/s')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'cpool_grain_storage_gr') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'cpool_grain_storage_gr',pptr%pcf%cpool_grain_storage_gr)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'cpool_grain_storage_gr',pptr%pcf%cpool_grain_storage_gr)
    end if

    ! transfer_grain_gr
    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'transfer_grain_gr',cdims=(/'pft'/), &
            long_name='grain growth respiration from storage',units='gC/m2/s')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'transfer_grain_gr') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'transfer_grain_gr',pptr%pcf%transfer_grain_gr)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'transfer_grain_gr',pptr%pcf%transfer_grain_gr)
    end if

    ! grainc_storage_to_xfer
    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'grainc_storage_to_xfer',cdims=(/'pft'/), &
            long_name='grain C shift storage to transfer',units='gC/m2/s')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'grainc_storage_to_xfer') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'grainc_storage_to_xfer',pptr%pcf%grainc_storage_to_xfer)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'grainc_storage_to_xfer',pptr%pcf%grainc_storage_to_xfer)
    end if

    ! grainn_storage_to_xfer
    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'grainn_storage_to_xfer',cdims=(/'pft'/), &
            long_name='grain N shift storage to transfer',units='gN/m2/s')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'grainn_storage_to_xfer') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'grainn_storage_to_xfer',pptr%pnf%grainn_storage_to_xfer)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'grainn_storage_to_xfer',pptr%pnf%grainn_storage_to_xfer)
    end if

  end subroutine CropRest

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CropRestYear
!
! !INTERFACE:
  integer function CropRestYear ( )
!
! !DESCRIPTION: 
! Return the restart year for prognostic crop
!
! !USES:
!
! !ARGUMENTS:
    implicit none
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!EOP
!
! !LOCAL VARIABLES:
     CropRestYear = restyear
     if ( CropRestYear == unset )then
        CropRestYear = 0
     end if
  end function CropRestYear

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: CropRestIncYear
!
! !INTERFACE:
  subroutine CropRestIncYear ( nyrs )
!
! !DESCRIPTION: 
! Increment the crop restart year
!
! !USES:
!
! !ARGUMENTS:
    implicit none
    integer, intent(out) :: nyrs ! Number of years crop has run
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!EOP
!
! !LOCAL VARIABLES:
      if ( restyear == unset ) restyear = 0
      restyear = restyear + 1
      nyrs     = restyear
  end subroutine CropRestIncYear

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: checkDates
!
! !INTERFACE:
  subroutine checkDates( )
!
! !DESCRIPTION: 
! Make sure the dates are compatible. The date given to startup the model
! and the date on the restart file must be the same although years can be
! different. The dates need to be checked when the restart file is being
! read in for a startup or branch case (they are NOT allowed to be different
! for a restart case).
!
! For the prognostic crop model the date of planting is tracked and growing
! degree days is tracked (with a 20 year mean) -- so shifting the start dates
! messes up these bits of saved information.
!
! !USES:
!
! !ARGUMENTS:
    use mod_clm_time_manager, only : get_driver_start_ymd, get_start_date
    use mod_clm_varctl      , only : nsrest, nsrBranch, nsrStartup
!
! !REVISION HISTORY:
! Author: Erik Kluzek
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: stymd       ! Start date YYYYMMDD from driver
    integer :: styr        ! Start year from driver
    integer :: stmon_day   ! Start date MMDD from driver
    integer :: rsmon_day   ! Restart date MMDD from restart file
    integer :: rsyr        ! Restart year from restart file
    integer :: rsmon       ! Restart month from restart file
    integer :: rsday       ! Restart day from restart file
    integer :: tod         ! Restart time of day from restart file
    character(len=*), parameter :: formDate = '(A,i4.4,"/",i2.2,"/",i2.2)' ! log output format
    character(len=32) :: subname = 'CropRest::checkDates'
    !
    ! If branch or startup make sure the startdate is compatible with the date
    ! on the restart file.
    !
    if ( nsrest == nsrBranch .or. nsrest == nsrStartup )then
       stymd       = get_driver_start_ymd()
       styr        = stymd / 10000
       stmon_day   = stymd - styr*10000
       call get_start_date( rsyr, rsmon, rsday, tod )
       rsmon_day = rsmon*100 + rsday
       if ( myid == italk ) &
       write(stdout,formDate) 'Date on the restart file is: ', rsyr, rsmon, rsday
       if ( stmon_day /= rsmon_day )then
          write(stdout,formDate) 'Start date is: ', styr, stmon_day/100, &
                                 (stmon_day - stmon_day/100)
          call fatal(__FILE__,__LINE__,trim(subname)// &
          ' ERROR: For prognostic crop to work correctly, the start date (month and day)'// &
          ' and the date on the restart file needs to match (years can be different)' )
       end if
    end if

  end subroutine checkDates

#endif

end module mod_clm_croprest

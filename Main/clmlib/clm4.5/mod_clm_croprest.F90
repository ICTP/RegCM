module mod_clm_croprest

#if (defined CN)
  !
  ! Read/Write to/from Crop info to CLM restart file.
  !
  use mod_intkinds
  use mod_realkinds
  use mod_mppparam
  use mod_dynparam
  use mod_mpmessage
  use mod_clm_nchelper

  implicit none

  private

  save

  public :: CropRest        ! Restart prognostic crop model
  public :: CropRestYear    ! Get the number of years crop has spunup
  public :: CropRestIncYear ! Increment the crop spinup years

  ! Flag that restart year is not set
  integer(ik4), parameter :: unset = -999
  ! Restart year from the initial conditions file
  integer(ik4) :: restyear = unset

  contains
  !
  ! Read/write Crop restart data
  !
  subroutine CropRest ( ncid, flag )
    use mod_clm_type
    use mod_clm_atmlnd , only : clm_a2l
    use mod_clm_varpar , only : numrad
    use mod_clm_decomp , only : get_proc_bounds , gcomm_pft
    use mod_clm_time_manager , only : is_restart

    implicit none
    type(clm_filetype) :: ncid             ! netcdf id
    character(len=*), intent(in) :: flag   !'read' or 'write'

    integer(ik4) :: begp, endp  ! per-proc beginning and ending pft indices
    integer(ik4) :: begc, endc  ! per-proc beginning and ending column indices
    integer(ik4) :: begl, endl  ! per-proc beginning and ending landunit indices
    integer(ik4) :: begg, endg  ! per-proc gridcell ending gridcell indices
    type(gridcell_type), pointer :: gptr  ! pointer to gridcell derived subtype
    type(landunit_type), pointer :: lptr  ! pointer to landunit derived subtype
    type(column_type)  , pointer :: cptr  ! pointer to column derived subtype
    type(pft_type)     , pointer :: pptr  ! pointer to pft derived subtype

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
        call clm_readvar(ncid,'peaklai',pptr%pps%peaklai,gcomm_pft)
      end if
    else if (flag == 'write') then
       call clm_writevar(ncid,'peaklai',pptr%pps%peaklai,gcomm_pft)
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
        call clm_readvar(ncid,'idop',pptr%pps%idop,gcomm_pft)
      end if
    else if ( flag == 'write') then
      call clm_writevar(ncid,'idop',pptr%pps%idop,gcomm_pft)
    end if

    ! aleaf
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'aleaf',cdims=(/'pft'/), &
            long_name='leaf allocation coefficient',units='')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'aleaf') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'aleaf',pptr%pps%aleaf,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'aleaf',pptr%pps%aleaf,gcomm_pft)
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
        call clm_readvar(ncid,'aleafi',pptr%pps%aleafi,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'aleafi',pptr%pps%aleafi,gcomm_pft)
    end if

    ! astem
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'astem',cdims=(/'pft'/), &
            long_name='stem allocation coefficient',units='')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'astem') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'astem',pptr%pps%astem,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'astem',pptr%pps%astem,gcomm_pft)
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
        call clm_readvar(ncid,'astemi',pptr%pps%astemi,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'astemi',pptr%pps%astemi,gcomm_pft)
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
        call clm_readvar(ncid,'htmx',pptr%pps%htmx,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'htmx',pptr%pps%htmx,gcomm_pft)
    end if

    ! hdidx
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'hdidx',cdims=(/'pft'/), &
            long_name='cold hardening index',units='')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'hdidx') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'hdidx',pptr%pps%hdidx,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'hdidx',pptr%pps%hdidx,gcomm_pft)
    end if

    ! vf
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'vf',cdims=(/'pft'/), &
            long_name='vernalization factor',units='')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'vf') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'vf',pptr%pps%vf,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'vf',pptr%pps%vf,gcomm_pft)
    end if

    ! cumvd
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'cumvd',cdims=(/'pft'/), &
            long_name='cumulative vernalization d',units='')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'cumvd') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'cumvd',pptr%pps%cumvd,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'cumvd',pptr%pps%cumvd,gcomm_pft)
    end if

    ! croplive
    if (flag == 'define') then
      call clm_addvar(clmvar_logical,ncid,'croplive',cdims=(/'pft'/), &
            long_name='Flag that crop is alive, but not harvested')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'croplive') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'croplive',pptr%pps%croplive,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'croplive',pptr%pps%croplive,gcomm_pft)
    end if

    ! cropplant
    if (flag == 'define') then
      call clm_addvar(clmvar_logical,ncid,'cropplant',cdims=(/'pft'/), &
            long_name='Flag that crop is planted, but not harvested' )
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'cropplant') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'cropplant',pptr%pps%cropplant,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'cropplant',pptr%pps%cropplant,gcomm_pft)
    end if

    ! harvdate
    if (flag == 'define') then
      call clm_addvar(clmvar_integer,ncid,'harvdate',cdims=(/'pft'/), &
            long_name='harvest date',units='jday',valid_range=(/1,366/) )
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'harvdate') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'harvdate',pptr%pps%harvdate,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'harvdate',pptr%pps%harvdate,gcomm_pft)
    end if

    ! gdd1020
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'gdd1020',cdims=(/'pft'/), &
      long_name = '20 year average of growing degree-days base &
                  &10C from planting', &
      units='ddays')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'gdd1020') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'gdd1020',pptr%pps%gdd1020,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'gdd1020',pptr%pps%gdd1020,gcomm_pft)
    end if

    ! gdd820
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'gdd820',cdims=(/'pft'/), &
       long_name = '20 year average of growing degree-days base 8C &
                   &from planting', &
       units='ddays')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'gdd820') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'gdd820',pptr%pps%gdd820,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'gdd820',pptr%pps%gdd820,gcomm_pft)
    end if

    ! gdd020
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'gdd020',cdims=(/'pft'/), &
       long_name = '20 year average of growing degree-days base 0C &
                   &from planting', &
       units='ddays')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'gdd020') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'gdd020',pptr%pps%gdd020,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'gdd020',pptr%pps%gdd020,gcomm_pft)
    end if

    ! gddmaturity
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'gddmaturity',cdims=(/'pft'/), &
            long_name='Growing degree days needed to harvest',units='ddays')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'gddmaturity') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'gddmaturity',pptr%pps%gddmaturity,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'gddmaturity',pptr%pps%gddmaturity,gcomm_pft)
    end if

    ! huileaf
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'huileaf',cdims=(/'pft'/), &
        long_name='heat unit index needed from planting to leaf emergence',&
        units='')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'huileaf') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'huileaf',pptr%pps%huileaf,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'huileaf',pptr%pps%huileaf,gcomm_pft)
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
        call clm_readvar(ncid,'huigrain',pptr%pps%huigrain,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'huigrain',pptr%pps%huigrain,gcomm_pft)
    end if

    ! grainc
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'grainc',cdims=(/'pft'/), &
            long_name='grain C',units='gC/m2')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'grainc') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'grainc',pptr%pcs%grainc,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'grainc',pptr%pcs%grainc,gcomm_pft)
    end if

    ! grainc_storage
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'grainc_storage',cdims=(/'pft'/), &
            long_name='grain C storage',units='gC/m2')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'grainc_storage') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'grainc_storage', &
                pptr%pcs%grainc_storage,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'grainc_storage', &
              pptr%pcs%grainc_storage,gcomm_pft)
    end if

    ! grainc_xfer
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'grainc_xfer',cdims=(/'pft'/), &
            long_name='grain C transfer',units='gC/m2')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'grainc_xfer') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'grainc_xfer',pptr%pcs%grainc_xfer,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'grainc_xfer',pptr%pcs%grainc_xfer,gcomm_pft)
    end if

    ! grainn
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'grainn',cdims=(/'pft'/), &
            long_name='grain N',units='gN/m2')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'grainn') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'grainn',pptr%pns%grainn,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'grainn',pptr%pns%grainn,gcomm_pft)
    end if

    ! grainn_storage
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'grainn_storage',cdims=(/'pft'/), &
            long_name='grain N storage',units='gN/m2')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'grainn_storage') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'grainn_storage', &
                pptr%pns%grainn_storage,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'grainn_storage', &
              pptr%pns%grainn_storage,gcomm_pft)
    end if

    ! grainn_xfer
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'grainn_xfer',cdims=(/'pft'/), &
            long_name='grain N transfer',units='gN/m2')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'grainn_xfer') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'grainn_xfer',pptr%pns%grainn_xfer,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'grainn_xfer',pptr%pns%grainn_xfer,gcomm_pft)
    end if

    ! grainc_xfer_to_grainc
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'grainc_xfer_to_grainc', &
         cdims=(/'pft'/),long_name='grain C growth from storage', &
         units='gC/m2/s')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. &
           clm_check_var(ncid,'grainc_xfer_to_grainc') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'grainc_xfer_to_grainc', &
                pptr%pcf%grainc_xfer_to_grainc,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'grainc_xfer_to_grainc', &
               pptr%pcf%grainc_xfer_to_grainc,gcomm_pft)
    end if

    ! livestemc_to_litter
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'livestemc_to_litter', &
            cdims=(/'pft'/), long_name='live stem C litterfall',units='gC/m2/s')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. &
           clm_check_var(ncid,'livestemc_to_litter') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'livestemc_to_litter', &
                pptr%pcf%livestemc_to_litter,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'livestemc_to_litter', &
              pptr%pcf%livestemc_to_litter,gcomm_pft)
    end if

    ! grainc_to_food
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'grainc_to_food',cdims=(/'pft'/), &
            long_name='grain C to food',units='gC/m2/s')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'grainc_to_food') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'grainc_to_food', &
                pptr%pcf%grainc_to_food,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'grainc_to_food', &
              pptr%pcf%grainc_to_food,gcomm_pft)
    end if

    ! grainn_xfer_to_grainn
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'grainn_xfer_to_grainn', &
         cdims=(/'pft'/),long_name='grain N growth from storage', &
         units='gN/m2/s')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. &
           clm_check_var(ncid,'grainn_xfer_to_grainn') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'grainn_xfer_to_grainn', &
                pptr%pnf%grainn_xfer_to_grainn,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'grainn_xfer_to_grainn', &
               pptr%pnf%grainn_xfer_to_grainn,gcomm_pft)
    end if

    ! livestemn_to_litter
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'livestemn_to_litter', &
            cdims=(/'pft'/), long_name='livestem N to litter',units='gN/m2/s')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. &
           clm_check_var(ncid,'livestemn_to_litter') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'livestemn_to_litter', &
                pptr%pnf%livestemn_to_litter,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'livestemn_to_litter', &
              pptr%pnf%livestemn_to_litter,gcomm_pft)
    end if

    ! grainn_to_food
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'grainn_to_food',cdims=(/'pft'/), &
            long_name='grain N to food',units='gN/m2/s')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'grainn_to_food') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'grainn_to_food', &
                pptr%pnf%grainn_to_food,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'grainn_to_food', &
              pptr%pnf%grainn_to_food,gcomm_pft)
    end if

    ! cpool_to_grainc
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'cpool_to_grainc',cdims=(/'pft'/), &
            long_name='allocation to grain C',units='gC/m2/s')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'cpool_to_grainc') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'cpool_to_grainc', &
                pptr%pcf%cpool_to_grainc,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'cpool_to_grainc', &
              pptr%pcf%cpool_to_grainc,gcomm_pft)
    end if

    ! cpool_to_grainc_storage
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'cpool_to_grainc_storage', &
        cdims=(/'pft'/),long_name='allocation to grain C storage',  &
        units='gC/m2/s')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. &
           clm_check_var(ncid,'cpool_to_grainc_storage') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'cpool_to_grainc_storage', &
                pptr%pcf%cpool_to_grainc_storage,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'cpool_to_grainc_storage', &
               pptr%pcf%cpool_to_grainc_storage,gcomm_pft)
    end if

    ! npool_to_grainn
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'npool_to_grainn',cdims=(/'pft'/), &
            long_name='allocation to grain N',units='gN/m2/s')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'npool_to_grainn') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'npool_to_grainn', &
                pptr%pnf%npool_to_grainn,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'npool_to_grainn', &
              pptr%pnf%npool_to_grainn,gcomm_pft)
    end if

    ! npool_to_grainn_storage
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'npool_to_grainn_storage', &
        cdims=(/'pft'/), &
        long_name='allocation to grain N storage',units='gN/m2/s')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. &
           clm_check_var(ncid,'npool_to_grainn_storage') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'npool_to_grainn_storage', &
                pptr%pnf%npool_to_grainn_storage,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'npool_to_grainn_storage', &
               pptr%pnf%npool_to_grainn_storage,gcomm_pft)
    end if

    ! cpool_grain_gr
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'cpool_grain_gr',cdims=(/'pft'/), &
            long_name='grain growth respiration',units='gC/m2/s')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. clm_check_var(ncid,'cpool_grain_gr') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'cpool_grain_gr', &
                pptr%pcf%cpool_grain_gr,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'cpool_grain_gr', &
              pptr%pcf%cpool_grain_gr,gcomm_pft)
    end if

    ! cpool_grain_storage_gr
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'cpool_grain_storage_gr', &
            cdims=(/'pft'/), &
            long_name='grain growth respiration to storage',units='gC/m2/s')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. &
           clm_check_var(ncid,'cpool_grain_storage_gr') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'cpool_grain_storage_gr', &
                pptr%pcf%cpool_grain_storage_gr,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'cpool_grain_storage_gr', &
               pptr%pcf%cpool_grain_storage_gr,gcomm_pft)
    end if

    ! transfer_grain_gr
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'transfer_grain_gr',cdims=(/'pft'/), &
            long_name='grain growth respiration from storage',units='gC/m2/s')
    else if (flag == 'read') then
      if ( is_restart() .and. &
              .not. clm_check_var(ncid,'transfer_grain_gr') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'transfer_grain_gr', &
                pptr%pcf%transfer_grain_gr,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'transfer_grain_gr', &
              pptr%pcf%transfer_grain_gr,gcomm_pft)
    end if

    ! grainc_storage_to_xfer
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'grainc_storage_to_xfer', &
            cdims=(/'pft'/), &
            long_name='grain C shift storage to transfer',units='gC/m2/s')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. &
           clm_check_var(ncid,'grainc_storage_to_xfer') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'grainc_storage_to_xfer', &
                pptr%pcf%grainc_storage_to_xfer,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'grainc_storage_to_xfer', &
               pptr%pcf%grainc_storage_to_xfer,gcomm_pft)
    end if

    ! grainn_storage_to_xfer
    if (flag == 'define') then
      call clm_addvar(clmvar_double,ncid,'grainn_storage_to_xfer', &
            cdims=(/'pft'/), &
            long_name='grain N shift storage to transfer',units='gN/m2/s')
    else if (flag == 'read') then
      if ( is_restart() .and. .not. &
           clm_check_var(ncid,'grainn_storage_to_xfer') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'grainn_storage_to_xfer', &
                pptr%pnf%grainn_storage_to_xfer,gcomm_pft)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'grainn_storage_to_xfer', &
               pptr%pnf%grainn_storage_to_xfer,gcomm_pft)
    end if

  end subroutine CropRest
  !
  ! Return the restart year for prognostic crop
  !
  integer(ik4) function CropRestYear ( )
    implicit none
    CropRestYear = restyear
    if ( CropRestYear == unset )then
      CropRestYear = 0
    end if
  end function CropRestYear
  !
  ! Increment the crop restart year
  !
  subroutine CropRestIncYear ( nyrs )
    implicit none
    integer(ik4), intent(out) :: nyrs ! Number of years crop has run
    if ( restyear == unset ) restyear = 0
    restyear = restyear + 1
    nyrs     = restyear
  end subroutine CropRestIncYear
  !
  ! Make sure the dates are compatible. The date given to startup the model
  ! and the date on the restart file must be the same although years can be
  ! different. The dates need to be checked when the restart file is being
  ! read in for a startup (they are NOT allowed to be different
  ! for a restart case).
  !
  ! For the prognostic crop model the date of planting is tracked and growing
  ! degree days is tracked (with a 20 year mean) -- so shifting the start dates
  ! messes up these bits of saved information.
  !
  subroutine checkDates( )
    use mod_clm_time_manager, only : get_driver_start_ymd, get_start_date
    use mod_clm_varctl , only : nsrest, nsrStartup
    integer(ik4) :: stymd       ! Start date YYYYMMDD from driver
    integer(ik4) :: styr        ! Start year from driver
    integer(ik4) :: stmon_day   ! Start date MMDD from driver
    integer(ik4) :: rsmon_day   ! Restart date MMDD from restart file
    integer(ik4) :: rsyr        ! Restart year from restart file
    integer(ik4) :: rsmon       ! Restart month from restart file
    integer(ik4) :: rsday       ! Restart day from restart file
    integer(ik4) :: tod         ! Restart time of day from restart file
    ! log output format
    character(len=*), parameter :: formDate = '(A,i4.4,"/",i2.2,"/",i2.2)'
    character(len=32) :: subname = 'CropRest::checkDates'
    !
    ! If startup make sure the startdate is compatible with the date
    ! on the restart file.
    !
    if ( nsrest == nsrStartup ) then
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
         ' ERROR: For prognostic crop to work correctly, the start &
         &date (month and day) and the date on the restart file &
         &needs to match (years can be different)' )
      end if
    end if
  end subroutine checkDates

#endif

end module mod_clm_croprest

module mod_clm_ch4rest
#ifdef LCH4
  !
  ! Reads from or writes restart data
  !
  use mod_intkinds
  use mod_realkinds
  use mod_mpmessage
  use mod_stdio
  use mod_runparams
  use mod_clm_type
  use mod_clm_decomp
  use mod_clm_nchelper
  use mod_clm_varctl , only : nsrest

  implicit none

  private

  save

  public :: ch4Rest

  contains
  !
  ! Read/Write biogeophysics information to/from restart file.
  !
  subroutine ch4Rest( ncid, flag )
    implicit none
    type(clm_filetype), intent(inout) :: ncid ! netcdf id
    character(len=*), intent(in) :: flag     ! 'read' or 'write'
    integer(ik4) :: c , l , g , j ! indices
    integer(ik4) :: begp , endp ! per-proc beginning and ending pft indices
    integer(ik4) :: begc , endc ! per-proc beginning and ending column indices
    integer(ik4) :: begl , endl ! per-proc beginning and ending landunit indices
    integer(ik4) :: begg , endg ! per-proc gridcell ending gridcell indices
    logical :: readvar  ! determine if variable is on initial file
    character(len=128) :: varname         ! temporary
    type(gridcell_type) , pointer :: gptr ! pointer to gridcell derived subtype
    type(landunit_type) , pointer :: lptr ! pointer to landunit derived subtype
    type(column_type) , pointer :: cptr   ! pointer to column derived subtype
    type(pft_type) , pointer :: pptr      ! pointer to pft derived subtype

    ! Set pointers into derived type

    gptr => clm3%g
    lptr => clm3%g%l
    cptr => clm3%g%l%c
    pptr => clm3%g%l%c%p

    ! column ch4 state variable - conc_ch4_sat

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CONC_CH4_SAT', &
              cdims=(/'column ','levgrnd'/), &
              long_name='methane soil concentration', units='mol/m^3')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CONC_CH4_SAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'CONC_CH4_SAT', &
                cptr%cch4%conc_ch4_sat,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CONC_CH4_SAT', &
              cptr%cch4%conc_ch4_sat,gcomm_column)
    end if

    ! column ch4 state variable - conc_ch4_unsat

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CONC_CH4_UNSAT', &
              cdims=(/'column ','levgrnd'/), &
              long_name='methane soil concentration', units='mol/m^3')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CONC_CH4_UNSAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'CONC_CH4_UNSAT', &
                cptr%cch4%conc_ch4_unsat,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CONC_CH4_UNSAT', &
              cptr%cch4%conc_ch4_unsat,gcomm_column)
    end if

    ! column ch4 state variable - conc_o2_sat

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CONC_O2_SAT', &
              cdims=(/'column ','levgrnd'/), &
              long_name='oxygen soil concentration', units='mol/m^3')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CONC_O2_SAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'CONC_O2_SAT', &
                cptr%cch4%conc_o2_sat,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CONC_O2_SAT', &
              cptr%cch4%conc_o2_sat,gcomm_column)
    end if

    ! column ch4 state variable - conc_o2_unsat

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'CONC_O2_UNSAT', &
              cdims=(/'column ','levgrnd'/), &
              long_name='oxygen soil concentration', units='mol/m^3')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'CONC_O2_UNSAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'CONC_O2_UNSAT', &
                cptr%cch4%conc_o2_unsat,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'CONC_O2_UNSAT', &
              cptr%cch4%conc_o2_unsat,gcomm_column)
    end if

    ! column ch4 flux variable - o2stress_sat (used in CNDecompCascade)

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'O2STRESS_SAT', &
              cdims=(/'column ','levgrnd'/), &
              long_name='oxygen stress fraction')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'O2STRESS_SAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'O2STRESS_SAT', &
                cptr%cch4%o2stress_sat,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'O2STRESS_SAT', &
              cptr%cch4%o2stress_sat,gcomm_column)
    end if

    ! column ch4 flux variable - o2stress_unsat (used in CNDecompCascade)

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'O2STRESS_UNSAT', &
              cdims=(/'column ','levgrnd'/), &
              long_name='oxygen stress fraction')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'O2STRESS_UNSAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'O2STRESS_UNSAT', &
                cptr%cch4%o2stress_unsat,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'O2STRESS_UNSAT', &
              cptr%cch4%o2stress_unsat,gcomm_column)
    end if

    ! column ch4 state variable - layer_sat_lag

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'LAYER_SAT_LAG', &
              cdims=(/'column ','levgrnd'/), &
              long_name='lagged saturation status of layer in unsat. zone')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'LAYER_SAT_LAG') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'LAYER_SAT_LAG', &
                cptr%cch4%layer_sat_lag,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'LAYER_SAT_LAG', &
              cptr%cch4%layer_sat_lag,gcomm_column)
    end if

    ! column ch4 state variable - qflx_surf_lag

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'QFLX_SURF_LAG', &
              cdims=(/'column '/), &
              long_name='time-lagged surface runoff', units='mm/s')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'QFLX_SURF_LAG') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'QFLX_SURF_LAG', &
                cptr%cch4%qflx_surf_lag,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'QFLX_SURF_LAG', &
              cptr%cch4%qflx_surf_lag,gcomm_column)
    end if

    ! column ch4 state variable - finundated_lag

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'FINUNDATED_LAG', &
              cdims=(/'column '/), &
              long_name='time-lagged inundated fraction')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'FINUNDATED_LAG') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'FINUNDATED_LAG', &
                cptr%cch4%finundated_lag,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'FINUNDATED_LAG', &
              cptr%cch4%finundated_lag,gcomm_column)
    end if

    ! column ch4 state variable - fsat_bef
    ! fsat_bef = finundated except inside methane code
    ! Only necessary for bit-for-bit restarts

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'FINUNDATED', &
              cdims=(/'column '/),long_name='inundated fraction')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'FINUNDATED') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'FINUNDATED',cptr%cch4%fsat_bef,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'FINUNDATED',cptr%cch4%fsat_bef,gcomm_column)
    end if

#ifdef CN

    ! column ch4 state variable - annavg_somhr
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'annavg_somhr', &
              cdims=(/'column '/), &
              long_name='Annual Average SOMHR',units='gC/m^2/s')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'annavg_somhr') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'annavg_somhr', &
                cptr%cch4%annavg_somhr,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'annavg_somhr', &
              cptr%cch4%annavg_somhr,gcomm_column)
    end if

    ! column ch4 state variable - annavg_finrw
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'annavg_finrw', &
              cdims=(/'column '/), &
              long_name='Annual Average Respiration-Weighted FINUNDATED')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'annavg_finrw') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'annavg_finrw', &
                cptr%cch4%annavg_finrw,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'annavg_finrw', &
              cptr%cch4%annavg_finrw,gcomm_column)
    end if

    ! column ch4 state variable - annsum_counter
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'annsum_counter_ch4', &
              cdims=(/'column '/), &
              long_name='CH4 Ann. Sum Time Counter')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'annsum_counter_ch4') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'annsum_counter_ch4', &
                cptr%cch4%annsum_counter,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'annsum_counter_ch4', &
              cptr%cch4%annsum_counter,gcomm_column)
    end if

    ! column ch4 state variable - tempavg_somhr
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'tempavg_somhr', &
              cdims=(/'column '/), &
              long_name='Temp. Average SOMHR',units='gC/m^2/s')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'tempavg_somhr') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'tempavg_somhr', &
                cptr%cch4%tempavg_somhr,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'tempavg_somhr', &
              cptr%cch4%tempavg_somhr,gcomm_column)
    end if

    ! column ch4 state variable - tempavg_finrwi
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'tempavg_finrw', &
              cdims=(/'column '/), &
              long_name='Temp. Average Respiration-Weighted FINUNDATED')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'tempavg_finrw') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'tempavg_finrw', &
                cptr%cch4%tempavg_finrw,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'tempavg_finrw', &
              cptr%cch4%tempavg_finrw,gcomm_column)
    end if

    ! pft ch4 state variable - tempavg_agnpp
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'tempavg_agnpp', &
              cdims=(/'pft '/), &
              long_name='Temp. Average AGNPP',units='gC/m^2/s')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'tempavg_agnpp') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'tempavg_agnpp',pptr%pcf%tempavg_agnpp,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'tempavg_agnpp',pptr%pcf%tempavg_agnpp,gcomm_pft)
    end if

    ! pft ch4 state variable - tempavg_bgnpp
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'tempavg_bgnpp', &
              cdims=(/'pft '/), &
              long_name='Temp. Average BGNPP',units='gC/m^2/s')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'tempavg_bgnpp') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'tempavg_bgnpp',pptr%pcf%tempavg_bgnpp,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'tempavg_bgnpp',pptr%pcf%tempavg_bgnpp,gcomm_pft)
    end if

    ! pft ch4 state variable - annavg_agnpp
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'annavg_agnpp', &
              cdims=(/'pft '/), &
              long_name='Ann. Average AGNPP',units='gC/m^2/s')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'annavg_agnpp') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'annavg_agnpp',pptr%pcf%annavg_agnpp,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'annavg_agnpp',pptr%pcf%annavg_agnpp,gcomm_pft)
    end if

    ! pft ch4 state variable - annavg_bgnpp
    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'annavg_bgnpp', &
              cdims=(/'pft '/), &
              long_name='Ann. Average BGNPP',units='gC/m^2/s')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'annavg_bgnpp') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'annavg_bgnpp',pptr%pcf%annavg_bgnpp,gcomm_pft)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'annavg_bgnpp',pptr%pcf%annavg_bgnpp,gcomm_pft)
    end if

   ! column ch4 flux variable - o2_decomp_depth_sat (used in CNNitrifDenitrif)

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'O2_DECOMP_DEPTH_SAT', &
              cdims=(/'column ','levgrnd'/), &
              long_name='O2 consumption during decomposition',units='mol/m3/s')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'O2_DECOMP_DEPTH_SAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'O2_DECOMP_DEPTH_SAT', &
                cptr%cch4%o2_decomp_depth_sat,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'O2_DECOMP_DEPTH_SAT', &
              cptr%cch4%o2_decomp_depth_sat,gcomm_column)
    end if

    ! column ch4 flux variable - o2_decomp_depth_unsat
    !   (used in CNNitrifDenitrif)

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'O2_DECOMP_DEPTH_UNSAT', &
              cdims=(/'column ','levgrnd'/), &
              long_name='O2 consumption during decomposition',units='mol/m3/s')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'O2_DECOMP_DEPTH_UNSAT') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'O2_DECOMP_DEPTH_UNSAT', &
                cptr%cch4%o2_decomp_depth_unsat,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'O2_DECOMP_DEPTH_UNSAT', &
              cptr%cch4%o2_decomp_depth_unsat,gcomm_column)
    end if

#endif

    ! column ch4 state variable - lake_soilc

    if ( flag == 'define' ) then
      call clm_addvar(clmvar_double,ncid,'LAKE_SOILC', &
              cdims=(/'column ','levgrnd'/), &
              long_name='lake soil carbon concentration', units='g/m^3')
    else if ( flag == 'read' ) then
      if ( .not. clm_check_var(ncid,'LAKE_SOILC') ) then
        if ( ktau > 0 ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        end if
      else
        call clm_readvar(ncid,'LAKE_SOILC', &
                cptr%cch4%lake_soilc,gcomm_column)
      end if
    else if ( flag == 'write' ) then
      call clm_writevar(ncid,'LAKE_SOILC', &
                cptr%cch4%lake_soilc,gcomm_column)
    end if

  end subroutine ch4Rest

#endif

end module mod_clm_ch4rest
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

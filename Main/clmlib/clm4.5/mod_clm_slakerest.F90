module mod_clm_slakerest
  !
  ! Reads from or writes restart data
  !
  use mod_intkinds
  use mod_realkinds
  use mod_mpmessage
  use mod_clm_nchelper
  use mod_clm_type
  use mod_clm_decomp
  use mod_runparams

  implicit none

  private

  save

  public :: SLakeRest

  contains
  !
  ! Read/Write biogeophysics information to/from restart file.
  !
  subroutine SLakeRest( ncid, flag )
    implicit none
    type(clm_filetype) , intent(inout) :: ncid ! netcdf id
    character(len=*) , intent(in) :: flag      ! 'read' or 'write'
    type(gridcell_type) , pointer :: gptr ! pointer to gridcell derived subtype
    type(landunit_type) , pointer :: lptr ! pointer to landunit derived subtype
    type(column_type) , pointer :: cptr   ! pointer to column derived subtype
    type(pft_type) , pointer :: pptr      ! pointer to pft derived subtype

    ! Set pointers into derived type

    gptr => clm3%g
    lptr => clm3%g%l
    cptr => clm3%g%l%c
    pptr => clm3%g%l%c%p

    ! Note t_lake is already in BiogeophysRest.

    ! column water state variable - lake_icefrac

    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'LAKE_ICEFRAC', &
            cdims=(/'column','levlak'/), &
            long_name='lake layer ice fraction',units='kg/kg', switchdim=.true.)
    else if (flag == 'read' ) then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'LAKE_ICEFRAC') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'LAKE_ICEFRAC', &
                cptr%cws%lake_icefrac,gcomm_column, switchdim=.true.)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'LAKE_ICEFRAC', &
              cptr%cws%lake_icefrac,gcomm_column, switchdim=.true.)
    end if

    ! column physical state variable - savedtke1

    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'SAVEDTKE1', &
            cdims=(/'column'/), &
            long_name='top lake layer eddy conductivity', units='W/(m K)')
    else if (flag == 'read' ) then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'SAVEDTKE1') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'SAVEDTKE1',cptr%cps%savedtke1,gcomm_column)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'SAVEDTKE1',cptr%cps%savedtke1,gcomm_column)
    end if

    ! column physical state variable - ust_lake

    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'USTLAKE', &
            cdims=(/'column'/), &
            long_name='friction velocity for lakes', units='m/s')
    else if (flag == 'read' ) then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'USTLAKE') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'USTLAKE',cptr%cps%ust_lake,gcomm_column)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'USTLAKE',cptr%cps%ust_lake,gcomm_column)
    end if

    ! column physical state variable - z0mg

    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'Z0MG', cdims=(/'column'/), &
            long_name='ground momentum roughness length', units='m')
    else if (flag == 'read' ) then
      if ( ktau /= 0 .and. .not. clm_check_var(ncid,'Z0MG') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'Z0MG',cptr%cps%z0mg,gcomm_column)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'Z0MG',cptr%cps%z0mg,gcomm_column)
    end if
  end subroutine SLakeRest

end module mod_clm_slakerest
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

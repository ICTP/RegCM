module mod_clm_slakerest

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: SLakeRestMod
!
! !DESCRIPTION:
! Reads from or writes restart data
!
! !USES:
  use mod_realkinds
  use mod_mpmessage
  use mod_clm_nchelper
!
! !PUBLIC TYPES:
  implicit none
  private
! save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: SLakeRest
!
! !REVISION HISTORY:
! 2009, June: Created by Zack Subin
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: SLakeRest
!
! !INTERFACE:
  subroutine SLakeRest( ncid, flag )
!
! !DESCRIPTION:
! Read/Write biogeophysics information to/from restart file.
!
! !USES:
    use mod_clm_type
    use mod_clm_decomp     , only : get_proc_bounds
    use mod_clm_time_manager , only : is_restart
!
! !ARGUMENTS:
    implicit none
    type(clm_filetype), intent(inout) :: ncid ! netcdf id
    character(len=*), intent(in) :: flag     ! 'read' or 'write'
!
! !CALLED FROM:
!
! !REVISION HISTORY:
! Author: Mariana Vertenstein
! 12/11/2003, Peter Thornton: Added cps%coszen, pps%gdir, and pps%omega
!   for new sunlit/shaded canopy algorithm (in SUNSHA ifdef block)
! 4/25/2005, Peter Thornton: Removed the SUNSHA ifdefs, since this is now the
!   default code behavior.
! 6/12/2005, Moved to netcdf format and renamed file
! 6/2009, Zack Subin: Adapted for S Lake physics.
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: c,l,g,j      ! indices
    integer :: begp, endp   ! per-proc beginning and ending pft indices
    integer :: begc, endc   ! per-proc beginning and ending column indices
    integer :: begl, endl   ! per-proc beginning and ending landunit indices
    integer :: begg, endg   ! per-proc gridcell ending gridcell indices
    logical :: readvar      ! determine if variable is on initial file
    character(len=128) :: varname         ! temporary
    type(gridcell_type), pointer :: gptr  ! pointer to gridcell derived subtype
    type(landunit_type), pointer :: lptr  ! pointer to landunit derived subtype
    type(column_type)  , pointer :: cptr  ! pointer to column derived subtype
    type(pft_type)     , pointer :: pptr  ! pointer to pft derived subtype
!-----------------------------------------------------------------------

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
            long_name='lake layer ice fraction',units='kg/kg')
    else if (flag == 'read' ) then
      if ( is_restart() .and. .not. clm_check_var(ncid,'LAKE_ICEFRAC') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'LAKE_ICEFRAC',cptr%cws%lake_icefrac)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'LAKE_ICEFRAC',cptr%cws%lake_icefrac)
    end if

    ! column physical state variable - savedtke1

    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'SAVEDTKE1', &
            cdims=(/'column'/), &
            long_name='top lake layer eddy conductivity', units='W/(m K)')
    else if (flag == 'read' ) then
      if ( is_restart() .and. .not. clm_check_var(ncid,'SAVEDTKE1') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'SAVEDTKE1',cptr%cps%savedtke1)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'SAVEDTKE1',cptr%cps%savedtke1)
    end if

    ! column physical state variable - ust_lake

    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'USTLAKE', &
            cdims=(/'column'/), &
            long_name='friction velocity for lakes', units='m/s')
    else if (flag == 'read' ) then
      if ( is_restart() .and. .not. clm_check_var(ncid,'USTLAKE') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'USTLAKE',cptr%cps%ust_lake)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'USTLAKE',cptr%cps%ust_lake)
    end if

    ! column physical state variable - z0mg

    if (flag == 'define') then
       call clm_addvar(clmvar_double,ncid,'Z0MG', cdims=(/'column'/), &
            long_name='ground momentum roughness length', units='m')
    else if (flag == 'read' ) then
      if ( is_restart() .and. .not. clm_check_var(ncid,'Z0MG') ) then
        call fatal(__FILE__,__LINE__,'clm now stopping')
      else
        call clm_readvar(ncid,'Z0MG',cptr%cps%z0mg)
      end if
    else if (flag == 'write') then
      call clm_writevar(ncid,'Z0MG',cptr%cps%z0mg)
    end if

  end subroutine SLakeRest

end module mod_clm_slakerest

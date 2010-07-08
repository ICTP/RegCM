#include <misc.h>
#include <preproc.h>

module VOCRestMod

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: VOCRestMod
!
! !DESCRIPTION:
! Reads from or volatile organic compound emissions restart/initial data
!
! !USES:
  use shr_kind_mod, only : r8 => shr_kind_r8
  use abortutils,   only : endrun
!
! !PUBLIC TYPES:
  implicit none
! save
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: VOCRest
!
! !REVISION HISTORY:
! 2005-06-12: Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: VOCRest
!
! !INTERFACE:
  subroutine VOCRest( ncid, flag )
!
! !DESCRIPTION:
! Read/Write volatile organic compound information to/from restart file.
!
! !USES:
    use clmtype
    use ncdio
    use decompMod     , only : get_proc_bounds
    use clm_varcon    , only : denice, denh2o
    use clm_varctl    , only : allocate_all_vegpfts, nsrest
!abt added below
    use clm_varvoc    , only : c24,c240,n24,n240
!abt added above
!
! !ARGUMENTS:
    implicit none
    integer, intent(in) :: ncid           ! netcdf id
    character(len=*), intent(in) :: flag  ! 'read' or 'write'
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
! 12/15/2008, modified biogeophysics restart to voc restart
!
!EOP
!
! !LOCAL VARIABLES:
    integer :: p,c,l,g,j    ! indices
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
    

    ! accumulation variable - c24/c240

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='c24', xtype=nf_int, &
               long_name='current timestep count until 24hr')
       call ncd_defvar(ncid=ncid, varname='c240', xtype=nf_int,  &
               long_name='current timestep count until 240hr')
    else if (flag == 'write' .or. flag == 'read') then
       call ncd_ioglobal(varname='c24', data=c24, ncid=ncid, flag=flag, readvar=readvar)
       call ncd_ioglobal(varname='c240' , data=c240 , ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft lai for current and previous month variable - monlai

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='monlai', xtype=nf_double,  &
            dim1name='pft', dim2name='numrad', &
            long_name='leaf area index for current and previous month', units=' ')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='monlai', data=pptr%pva%monlai, &
            dim1name='pft', dim2name='numrad', &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft previous temperature/ppfd variable - t_sum24

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='TEMP_24', xtype=nf_double,  &
            dim1name='pft', dim2name='nday', &
            long_name='past 24 hour temp', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='TEMP_24', data=pptr%pva%t_sum24, &
!            dim1name='nday', dim2name='pft', &
            dim1name='pft', dim2name='nday', &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft previous temperature/ppfd variable - t_sum240

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='TEMP_240', xtype=nf_double,  &
!            dim1name='nten', dim2name='pft', &
            dim1name='pft', dim2name='nten', &
            long_name='past 240 hour temp', units='K')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='TEMP_240', data=pptr%pva%t_sum240, &
!            dim1name='nten', dim2name='pft', &
            dim1name='pft', dim2name='nten', &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft previous temperature/ppfd variable - p_sum24su

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='PPFD_24su', xtype=nf_double,  &
!            dim1name='nday', dim2name='pft', &
            dim1name='pft', dim2name='nday', &
            long_name='past 24 hour ppfd for sunlit portion', units='umol/m2/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='PPFD_24su', data=pptr%pva%p_sum24su, &
!            dim1name='nday', dim2name='pft', &
            dim1name='pft', dim2name='nday', &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft previous temperature/ppfd variable - p_sum240su

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='PPFD_240su', xtype=nf_double,  &
!            dim1name='nten', dim2name='pft', &
            dim1name='pft', dim2name='nten', &
            long_name='past 240 hour ppfd for sunlit', units='umol/m2/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='PPFD_240su', data=pptr%pva%p_sum240su, &
!            dim1name='nten', dim2name='pft', &
            dim1name='pft', dim2name='nten', &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft previous temperature/ppfd variable - p_sum24sh

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='PPFD_24sh', xtype=nf_double,  &
            dim1name='pft', dim2name='nday', &
            long_name='past 24 hour ppfd for shaded', units='umol/m2/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='PPFD_24sh', data=pptr%pva%p_sum24sh, &
            dim1name='pft', dim2name='nday', &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if

    ! pft previous temperature/ppfd variable - p_sum240sh

    if (flag == 'define') then
       call ncd_defvar(ncid=ncid, varname='PPFD_240sh', xtype=nf_double,  &
            dim1name='pft', dim2name='nten', &
            long_name='past 240 hour ppfd for shaded', units='umol/m2/s')
    else if (flag == 'read' .or. flag == 'write') then
       call ncd_iolocal(varname='PPFD_240sh', data=pptr%pva%p_sum240sh, &
            dim1name='pft', dim2name='nten', &
            ncid=ncid, flag=flag, readvar=readvar)
       if (flag=='read' .and. .not. readvar) then
          if (is_restart()) call endrun()
       end if
    end if



  end subroutine VOCRest

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: is_restart
!
! !INTERFACE:
  logical function is_restart( )
!
! !DESCRIPTION:
! Determine if restart run
!
! !USES:
    use clm_varctl, only : nsrest
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! subroutine initialize in this module
!
! !REVISION HISTORY:
! Created by Mariana Vertenstein
!
!EOP
!-----------------------------------------------------------------------

    if (nsrest == 1) then
       is_restart = .true.
    else
       is_restart = .false.
    end if

  end function is_restart

end module VOCRestMod

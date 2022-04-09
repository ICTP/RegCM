module mod_clm_accumul
  !
  ! This module contains generic subroutines that can be used to
  ! define, accumulate and extract  user-specified fields over
  ! user-defined intervals. Each interval  and accumulation type is
  ! unique to each field processed.
  ! Subroutine [init_accumulator] defines the values of the accumulated
  ! field data structure. Subroutine [update_accum_field] does
  ! the actual accumulation for a given field.
  ! Four types of accumulations are possible:
  ! - Average over time interval. Time average fields are only
  !   valid at the end of the averaging interval.
  ! - Running mean over time interval. Running means are valid once the
  !   length of the simulation exceeds the
  ! - Running accumulation over time interval. Accumulated fields are
  !   continuously accumulated. The trigger value "-99999." resets
  !   the accumulation to zero.
  !
  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_mpmessage
  use mod_dynparam
  use mod_runparams
  use mod_mppparam
  use mod_clm_nchelper
  use mod_clm_decomp
  use mod_clm_varcon , only : spval , secspday

  implicit none

  private

  save

  ! Write/read restart of accumulation fields
  public :: accumulRest
  ! Initialize an accumulator field
  public :: init_accum_field
  ! Print info about accumulator fields
  public :: print_accum_fields
  ! Extracts the current value of an accumulator field
  public :: extract_accum_field
  ! Update the current value of an accumulator field
  public :: update_accum_field

  interface extract_accum_field
    ! Extract current val of single-level accumulator field
    module procedure extract_accum_field_sl
    ! Extract current val of multi-level accumulator field
    module procedure extract_accum_field_ml
  end interface

  ! Updates the current value of an accumulator field
  interface update_accum_field
    ! Update single-level accumulator field
    module procedure update_accum_field_sl
    ! Update multi-level accumulator field
    module procedure update_accum_field_ml
  end interface

  type accum_field
    character(len=  8) :: fname     !field name
    character(len=128) :: desc     !field description
    character(len=  8) :: units    !field units
    !accumulation type: ["timeavg","runmean","runaccum"]
    character(len=  8) :: acctype
    !subgrid type: ["gridcell","landunit","column" or "pft"]
    character(len=  8) :: type1d
    !type2d ('','levsoi','numrad',..etc. )
    character(len=  8) :: type2d
    integer(ik4) :: beg1d  !subgrid type beginning index
    integer(ik4) :: end1d  !subgrid type ending index
    integer(ik4) :: num1d  !total subgrid points
    integer(ik4) :: numlev !number of vertical levels in field
    real(rk8):: initval    !initial value of accumulated field
    real(rk8), pointer :: val(:,:)  !accumulated field
    integer(ik8) :: period  !field accumulation period (in model time steps)
    type(subgrid_type) , pointer :: gcomm
  end type accum_field

  !maximum number of accumulated fields
  integer(ik4) , parameter :: max_accum = 100
  type (accum_field) :: accum(max_accum)   !array accumulated fields
  integer(ik4) :: naccflds = 0             !accumulator field counter

  contains
  !
  ! Initialize accumulation fields. This subroutine sets:
  ! o name  of accumulated field
  ! o units of accumulated field
  ! o accumulation type of accumulated field
  ! o description of accumulated fields: accdes
  ! o accumulation period for accumulated field (in iterations)
  ! o initial value of accumulated field
  !
  subroutine init_accum_field (fname, units, desc, &
       accum_type, accum_period, numlev, subgrid_type, init_value, type2d)
    implicit none
    character(len=*), intent(in) :: fname             !field name
    character(len=*), intent(in) :: units            !field units
    character(len=*), intent(in) :: desc             !field description
    !field type: tavg, runm, runa, ins
    character(len=*), intent(in) :: accum_type
    !field accumulation period
    integer(ik4) , intent(in) :: accum_period
    !["gridcell","landunit","column" or "pft"]
    character(len=*), intent(in) :: subgrid_type
    !number of vertical levels
    integer(ik4) , intent(in) :: numlev
    !field initial or reset value
    real(rk8), intent(in) :: init_value
    !level type (optional) - needed if numlev > 1
    character(len=*), intent(in), optional :: type2d
    integer(ik4) :: nf             ! field index
    integer(ik4) :: beg1d , end1d  ! beggining and end subgrid indices
    integer(ik4) :: num1d          ! total number subgrid indices
    integer(ik4) :: begp , endp ! per-proc beginning and ending pft indices
    integer(ik4) :: begc , endc ! per-proc beginning and ending column indices
    integer(ik4) :: begl , endl ! per-proc beginning and ending landunit indices
    integer(ik4) :: begg , endg ! per-proc gridcell ending gridcell indices
    integer(ik4) :: numg    ! total number of gridcells across all processors
    integer(ik4) :: numl    ! total number of landunits across all processors
    integer(ik4) :: numc    ! total number of columns across all processors
    integer(ik4) :: nump    ! total number of pfts across all processors

    ! Determine necessary indices

    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)
    call get_proc_global(numg,numl,numc,nump)

    ! update field index
    ! Consistency check that number of accumulated does not exceed maximum.

    naccflds = naccflds + 1
    if ( naccflds > max_accum ) then
      write(stderr,*) 'ACCUMULINIT error: user-defined accumulation fields ', &
           'equal to ',naccflds,' exceeds max_accum'
      call fatal(__FILE__,__LINE__,'clm now stopping')
    end if
    nf = naccflds

    ! Note accumulation period must be converted from days
    ! to number of iterations

    accum(nf)%fname    = trim(fname)
    accum(nf)%units   = trim(units)
    accum(nf)%desc    = trim(desc)
    accum(nf)%acctype = trim(accum_type)
    accum(nf)%initval = init_value
    accum(nf)%period  = accum_period
    if (accum(nf)%period < 0) then
      accum(nf)%period = -accum(nf)%period * nint(secspday/dtsrf)
    end if
    select case (trim(subgrid_type))
      case ('gridcell')
        beg1d = begg
        end1d = endg
        num1d = numg
        accum(nf)%gcomm => gcomm_gridcell
      case ('landunit')
        beg1d = begl
        end1d = endl
        num1d = numl
        accum(nf)%gcomm => gcomm_landunit
      case ('column')
        beg1d = begc
        end1d = endc
        num1d = numc
        accum(nf)%gcomm => gcomm_column
      case ('pft')
        beg1d = begp
        end1d = endp
        num1d = nump
        accum(nf)%gcomm => gcomm_pft
      case default
        write(stderr,*)'ACCUMULINIT: unknown subgrid type ',subgrid_type
        call fatal(__FILE__,__LINE__,'clm now stopping')
    end select

    accum(nf)%type1d = trim(subgrid_type)
    accum(nf)%beg1d = beg1d
    accum(nf)%end1d = end1d
    accum(nf)%num1d = num1d
    accum(nf)%numlev = numlev

    if (present(type2d)) then
      accum(nf)%type2d = type2d
    else
      accum(nf)%type2d = ' '
    end if

    ! Allocate and initialize accumulation field

    allocate(accum(nf)%val(beg1d:end1d,numlev))
    accum(nf)%val(beg1d:end1d,numlev) = init_value
  end subroutine init_accum_field
  !
  ! Diagnostic printout of accumulated fields
  !
  subroutine print_accum_fields()
    implicit none
    integer(ik4) :: i , nf   !indices
    ! Do not bloat output !
    if ( debug_level > 3 ) then
      if ( myid == italk ) then
        write(stdout,*)
        write(stdout,*) 'Initializing variables for time accumulation .....'
        write(stdout,'(72a1)') ("-",i=1,60)
        write(stdout,*) 'Accumulated fields'
        write(stdout,1002)
        write(stdout,'(72a1)') ("_",i=1,71)
        do nf = 1 , naccflds
          if (accum(nf)%period /= bigint) then
            write(stdout,1003) nf,accum(nf)%fname,accum(nf)%units,&
                  accum(nf)%acctype, accum(nf)%period, accum(nf)%initval, &
                  accum(nf)%desc
          else
            write(stdout,1004) nf,accum(nf)%fname,accum(nf)%units,&
                  accum(nf)%acctype, accum(nf)%initval, accum(nf)%desc
          end if
        end do
        write(stdout,'(72a1)') ("_",i=1,71)
        write(stdout,*)
        write(stdout,'(72a1)') ("-",i=1,60)
     end if
     write(stdout,*) 'Successfully initialized variables for accumulation'
    end if
1002 format(' No',' Name    ',' Units   ',' Type    ', &
            'Period',' Inival',' Description')
1003 format((1x,i2),(1x,a8),(1x,a8),(1x,a8), (1x,i5),(1x,f4.0),(1x,a40))
1004 format((1x,i2),(1x,a8),(1x,a8),(1x,a8),'  N.A.',(1x,f4.0),(1x,a40))
  end subroutine print_accum_fields
  !
  ! Extract single-level accumulated field.
  ! This routine extracts the field values from the multi-level
  ! accumulation field. It extracts the current value except if
  ! the field type is a time average. In this case, an absurd value
  ! is assigned to  indicate the time average is not yet valid.
  !
  subroutine extract_accum_field_sl (fname, field, nstep)
    implicit none
    character(len=*) , intent(in) :: fname     !field name
    !field values for current time step
    real(rk8) , pointer , dimension(:) :: field
    integer(ik8) , intent(in) :: nstep         !timestep index
    integer(ik4) :: i , k , nf        !indices
    integer(ik4) :: ibeg , iend         !subgrid beginning,ending indices
!------------------------------------------------------------------------

    ! find field index. return if "name" is not on list

    nf = 0
    do i = 1 , naccflds
       if (fname == accum(i)%fname) nf = i
    end do
    if ( nf == 0 ) then
      write(stderr,*) &
              'EXTRACT_ACCUM_FIELD_SL error: field name ',fname,' not found'
      call fatal(__FILE__,__LINE__,'clm now stopping')
    end if

    ! error check

    ibeg = accum(nf)%beg1d
    iend = accum(nf)%end1d
    if ( size(field,dim=1) /= iend-ibeg+1 ) then
      write(stderr,*)'ERROR in extract_accum_field for field ',accum(nf)%fname
      write(stderr,*)'size of first dimension of field is ',&
            size(field,dim=1),' and should be ',iend-ibeg+1
      call fatal(__FILE__,__LINE__,'clm now stopping')
    end if
    ! extract field
    if (accum(nf)%acctype == 'timeavg' .and. &
        mod(nstep,accum(nf)%period) /= 0) then
      do k = ibeg , iend
        field(k) = spval  !assign absurd value when avg not ready
      end do
    else
      do k = ibeg , iend
        field(k) = accum(nf)%val(k,1)
      end do
    end if
  end subroutine extract_accum_field_sl
  !
  ! Extract mutli-level accumulated field.
  ! This routine extracts the field values from the multi-level
  ! accumulation field. It extracts the current value except if
  ! the field type is a time average. In this case, an absurd value
  ! is assigned to  indicate the time average is not yet valid.
  !
  subroutine extract_accum_field_ml (fname, field, nstep)
    implicit none
    character(len=*) , intent(in) :: fname       !field name
    !field values for current time step
    real(rk8) , pointer , dimension(:,:) :: field
    integer(ik8) , intent(in) :: nstep           !timestep index
    integer(ik4) :: i , j , k , nf     !indices
    integer(ik4) :: ibeg , iend        !subgrid beginning,ending indices
    integer(ik4) :: numlev             !number of vertical levels

    ! find field index. return if "name" is not on list

    nf = 0
    do i = 1, naccflds
      if (fname == accum(i)%fname) nf = i
    end do
    if ( nf == 0 ) then
      write(stderr,*) &
              'EXTRACT_ACCUM_FIELD_ML error: field name ',fname,' not found'
      call fatal(__FILE__,__LINE__,'clm now stopping')
    end if

    ! error check

    numlev = accum(nf)%numlev
    ibeg = accum(nf)%beg1d
    iend = accum(nf)%end1d
    if ( size(field,dim=1) /= iend-ibeg+1 ) then
      write(stderr,*)'ERROR in extract_accum_field for field ',accum(nf)%fname
      write(stderr,*)'size of first dimension of field is ',&
            size(field,dim=1),' and should be ',iend-ibeg+1
      call fatal(__FILE__,__LINE__,'clm now stopping')
    else if ( size(field,dim=2) /= numlev ) then
      write(stderr,*)'ERROR in extract_accum_field for field ',accum(nf)%fname
      write(stderr,*)'size of second dimension of field iis ',&
            size(field,dim=2),' and should be ',numlev
      call fatal(__FILE__,__LINE__,'clm now stopping')
    end if

    !extract field

    if ( accum(nf)%acctype == 'timeavg' .and. &
         mod(nstep,accum(nf)%period) /= 0 ) then
      do j = 1 , numlev
        do k = ibeg , iend
          field(k,j) = spval  !assign absurd value when avg not ready
        end do
      end do
    else
      do j = 1 , numlev
        do k = ibeg , iend
          field(k,j) = accum(nf)%val(k,j)
        end do
      end do
    end if
  end subroutine extract_accum_field_ml
  !
  ! Accumulate single level field over specified time interval.
  ! The appropriate field is accumulated in the array [accval].
  !
  subroutine update_accum_field_sl (fname, field, nstep)
    implicit none
    character(len=*) , intent(in) :: fname     !field name
    !field values for current time step
    real(rk8) , pointer , dimension(:) :: field
    integer(ik8) , intent(in) :: nstep   !time step index
    integer(ik4) :: i , k , nf           !indices
    integer(ik4) :: accper               !temporary accumulation period
    integer(ik4) :: ibeg , iend          !subgrid beginning,ending indices

    ! find field index. return if "name" is not on list

    nf = 0
    do i = 1 , naccflds
      if ( fname == accum(i)%fname ) nf = i
    end do
    if ( nf == 0 ) then
      write(stderr,*) &
              'UPDATE_ACCUM_FIELD_SL error: field name ',fname,' not found'
      call fatal(__FILE__,__LINE__,'clm now stopping')
    end if

    ! error check

    ibeg = accum(nf)%beg1d
    iend = accum(nf)%end1d
    if ( size(field,dim=1) /= iend-ibeg+1 ) then
      write(stderr,*)'ERROR in UPDATE_ACCUM_FIELD_SL for field ',accum(nf)%fname
      write(stderr,*)'size of first dimension of field is ',size(field,dim=1),&
            ' and should be ',iend-ibeg+1
      call fatal(__FILE__,__LINE__,'clm now stopping')
    end if

    ! accumulate field

    if (accum(nf)%acctype /= 'timeavg' .AND. &
        accum(nf)%acctype /= 'runmean' .AND. &
        accum(nf)%acctype /= 'runaccum') then
      write(stderr,*) &
              'UPDATE_ACCUM_FIELD_SL error: incorrect accumulation type'
      write(stderr,*) ' was specified for field ',fname
      write(stderr,*)' accumulation type specified is ',accum(nf)%acctype
      write(stderr,*) &
              ' only [timeavg, runmean, runaccum] are currently acceptable'
      call fatal(__FILE__,__LINE__,'clm now stopping')
    end if

    ! reset accumulated field value if necessary and  update
    ! accumulation field
    ! running mean never reset

    if ( accum(nf)%acctype == 'timeavg' ) then

      !time average field reset every accumulation period
      !normalize at end of accumulation period

      if ( (mod(nstep,accum(nf)%period) == 1 .or. &
            accum(nf)%period == 1) .and. (nstep /= 0) ) then
        accum(nf)%val(ibeg:iend,1) = 0._rk8
      end if
      accum(nf)%val(ibeg:iend,1) = accum(nf)%val(ibeg:iend,1) + field(ibeg:iend)
      if ( mod(nstep,accum(nf)%period) == 0 ) then
        accum(nf)%val(ibeg:iend,1) = &
                accum(nf)%val(ibeg:iend,1) / accum(nf)%period
      end if

    else if (accum(nf)%acctype == 'runmean') then

      !running mean - reset accumulation period until greater than nstep

      accper = int(min (nstep,accum(nf)%period), ik4)
      accum(nf)%val(ibeg:iend,1) = &
           ((accper-1)*accum(nf)%val(ibeg:iend,1) + field(ibeg:iend)) / accper

    else if (accum(nf)%acctype == 'runaccum') then

      !running accumulation field reset at trigger -99999

      do k = ibeg , iend
        if ( nint(field(k)) == -99999 ) then
          accum(nf)%val(k,1) = 0._rk8
        end if
      end do
      accum(nf)%val(ibeg:iend,1) = &
        min(max(accum(nf)%val(ibeg:iend,1) + field(ibeg:iend), 0._rk8), 99999._rk8)
    end if
  end subroutine update_accum_field_sl
  !
  ! Accumulate multi level field over specified time interval.
  !
  subroutine update_accum_field_ml (fname, field, nstep)
    implicit none
    character(len=*) , intent(in) :: fname        !field name
    !field values for current time step
    real(rk8) , pointer , dimension(:,:) :: field
    integer(ik8) , intent(in) :: nstep    !time step index
    integer(ik4) :: i , j , k , nf        !indices
    integer(ik4) :: accper                !temporary accumulation period
    integer(ik4) :: ibeg , iend           !subgrid beginning,ending indices
    integer(ik4) :: numlev                !number of vertical levels

    ! find field index. return if "name" is not on list

    nf = 0
    do i = 1 , naccflds
      if ( fname == accum(i)%fname ) nf = i
    end do
    if ( nf == 0 ) then
      write(stderr,*) &
              'UPDATE_ACCUM_FIELD_ML error: field name ',fname,' not found'
      call fatal(__FILE__,__LINE__,'clm now stopping')
    end if

    ! error check

    numlev = accum(nf)%numlev
    ibeg = accum(nf)%beg1d
    iend = accum(nf)%end1d
    if ( size(field,dim=1) /= iend-ibeg+1 ) then
      write(stderr,*) &
            'ERROR in UPDATE_ACCUM_FIELD_ML for field ',accum(nf)%fname
      write(stderr,*) &
            'size of first dimension of field is ',size(field,dim=1),&
            ' and should be ',iend-ibeg+1
      call fatal(__FILE__,__LINE__,'clm now stopping')
    else if ( size(field,dim=2) /= numlev ) then
      write(stderr,*) &
            'ERROR in UPDATE_ACCUM_FIELD_ML for field ',accum(nf)%fname
      write(stderr,*) &
            'size of second dimension of field is ',size(field,dim=2),&
            ' and should be ',numlev
      call fatal(__FILE__,__LINE__,'clm now stopping')
    end if

    ! accumulate field

    if (accum(nf)%acctype /= 'timeavg' .AND. &
        accum(nf)%acctype /= 'runmean' .AND. &
        accum(nf)%acctype /= 'runaccum') then
      write(stderr,*) &
              'UPDATE_ACCUM_FIELD_ML error: incorrect accumulation type'
      write(stderr,*) ' was specified for field ',fname
      write(stderr,*) &
              ' accumulation type specified is ',accum(nf)%acctype
      write(stderr,*) &
              ' only [timeavg, runmean, runaccum] are currently acceptable'
      call fatal(__FILE__,__LINE__,'clm now stopping')
    end if

    ! accumulate field

    ! reset accumulated field value if necessary and  update
    ! accumulation field
    ! running mean never reset

    if ( accum(nf)%acctype == 'timeavg' ) then

      !time average field reset every accumulation period
      !normalize at end of accumulation period

      if (( mod(nstep,accum(nf)%period) == 1 .or. &
            accum(nf)%period == 1) .and. (nstep /= 0) ) then
        accum(nf)%val(ibeg:iend,1:numlev) = 0._rk8
      end if
      accum(nf)%val(ibeg:iend,1:numlev) = &
              accum(nf)%val(ibeg:iend,1:numlev) + field(ibeg:iend,1:numlev)
      if ( mod(nstep,accum(nf)%period) == 0 ) then
        accum(nf)%val(ibeg:iend,1:numlev) = &
                accum(nf)%val(ibeg:iend,1:numlev) / accum(nf)%period
      end if

    else if ( accum(nf)%acctype == 'runmean' ) then

      !running mean - reset accumulation period until greater than nstep

      accper = int(min (nstep,accum(nf)%period) , ik4)
      accum(nf)%val(ibeg:iend,1:numlev) = &
            ((accper-1)*accum(nf)%val(ibeg:iend,1:numlev) + &
              field(ibeg:iend,1:numlev)) / accper

    else if ( accum(nf)%acctype == 'runaccum' ) then

      !running accumulation field reset at trigger -99999

      do j = 1 , numlev
        do k = ibeg , iend
          if (nint(field(k,j)) == -99999) then
            accum(nf)%val(k,j) = 0._rk8
          end if
        end do
      end do
      accum(nf)%val(ibeg:iend,1:numlev) = &
            min(max(accum(nf)%val(ibeg:iend,1:numlev) + &
              field(ibeg:iend,1:numlev), 0._rk8), 99999._rk8)
    end if
  end subroutine update_accum_field_ml
  !
  ! Read/write accumulation restart data
  !
  subroutine accumulRest( ncid, flag )
    implicit none
    type(clm_filetype), intent(inout) :: ncid   !netcdf unit
    character(len=*) , intent(in) :: flag   !'define','read', or 'write'
    integer(ik4) :: nf , iper       ! indices
    integer(ik4) :: beg1d , end1d       ! buffer bounds
    real(rk8), pointer :: rbuf1d(:) ! temporary 1d buffer
    character(len=128) :: varname  ! temporary
    logical :: lstart

    lstart = rcmtimer%integrating( )

    do nf = 1 , naccflds
      ! Note = below need to allocate rbuf for single level variables, since
      ! accum(nf)%val is always 2d
      varname = trim(accum(nf)%fname) // '_VALUE'
      if ( flag == 'define' ) then
        if ( accum(nf)%numlev == 1 ) then
          call clm_addvar(clmvar_double,ncid,varname, &
                        [accum(nf)%type1d],accum(nf)%desc, &
                        accum(nf)%units)
        else
          call clm_addvar(clmvar_double,ncid,varname, &
                  [accum(nf)%type1d,accum(nf)%type2d], &
                   accum(nf)%desc,accum(nf)%units)
        end if
      else if ( flag == 'read' ) then
        if ( lstart .and. .not. clm_check_var(ncid,varname) ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        else
          beg1d = accum(nf)%beg1d
          end1d = accum(nf)%end1d
          allocate(rbuf1d(beg1d:end1d))
          call clm_readvar(ncid,varname,rbuf1d,accum(nf)%gcomm)
          accum(nf)%val(beg1d:end1d,1) = rbuf1d(beg1d:end1d)
          deallocate(rbuf1d)
        end if
      else if ( flag == 'write' ) then
        beg1d = accum(nf)%beg1d
        end1d = accum(nf)%end1d
        allocate(rbuf1d(beg1d:end1d))
        rbuf1d(beg1d:end1d) = accum(nf)%val(beg1d:end1d,1)
        call clm_writevar(ncid,varname,rbuf1d,accum(nf)%gcomm)
        deallocate(rbuf1d)
      end if

      varname = trim(accum(nf)%fname) // '_PERIOD'
      if ( flag == 'define' ) then
        call clm_addvar(clmvar_integer,ncid,varname, &
                        long_name='',units='time steps')
      else if ( flag == 'read' ) then
        if ( lstart .and. .not. clm_check_var(ncid,varname) ) then
          call fatal(__FILE__,__LINE__,'clm now stopping')
        else
          call clm_readvar(ncid,varname,iper)
          accum(nf)%period = iper
        end if
      else if ( flag == 'write' ) then
        iper = int(accum(nf)%period, ik4)
        call clm_writevar(ncid,varname,iper)
      end if
    end do
  end subroutine accumulRest

end module mod_clm_accumul
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

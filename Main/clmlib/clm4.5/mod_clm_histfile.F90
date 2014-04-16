module mod_clm_histfile
  !
  ! Module containing methods to for CLM history file handling.
  !
  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_mpmessage
  use mod_dynparam
  use mod_runparams
  use mod_mppparam
  use mod_date
  use mod_clm_type
  use mod_clm_nchelper
  use mod_clm_decomp
  use mod_clm_varcon , only : spval , ispval
  use mod_clm_varcon , only : dzsoi_decomp
  use mod_clm_subgridave , only : p2g , c2g , l2g
  use mod_clm_time_manager , only : get_prev_time
  use mod_clm_varcon , only : zsoi , zlak
  use mod_clm_varcon , only : secspday
  use mod_clm_varpar , only : nlevgrnd , nlevlak , nlevurb , numrad , &
         nlevdecomp_full
  use mod_clm_varctl , only : caseid , ctitle , fsurdat , finidat , fpftcon , &
         version , hostname , username , conventions , source , inst_suffix , &
         nsrest , nsrStartup
  use mod_clm_domain , only : ldomain
  use mod_clm_time_manager , only : get_prev_date , getdatetime

  implicit none

  private

  save

  ! max number of history tapes
  integer(ik4) , public , parameter :: max_tapes = 6
  ! max number of history fields
  integer(ik4) , public , parameter :: max_flds = 2500
  ! maximum number of characters for field name
  integer(ik4) , public , parameter :: max_namlen = 32
  !
  ! Counters
  !
  ! index of max history file requested
  integer(ik4) , public :: ntapes = 0
  !
  ! Namelist
  !
  integer(ik4) :: ni ! implicit index below
  ! namelist: flag indicates no default history fields
  logical , public :: hist_empty_htapes  = .false.
  ! namelist: output density of netcdf history files
  integer(ik4) , public :: hist_ndens(max_tapes) = 2
  ! namelist: number of time samples per tape
  integer(ik4) , public :: hist_mfilt(max_tapes) = 30
  ! namelist: true=> do grid averaging
  logical , public :: hist_dov2xy(max_tapes) = &
          (/.true.,(.true.,ni=2,max_tapes)/)
  ! namelist: history write freq(0=monthly)
  integer(ik4) , public :: hist_nhtfrq(max_tapes) = &
          (/0, (-24, ni=2,max_tapes)/)
  ! namelist: per tape averaging flag
  character(len=1) , public :: hist_avgflag_pertape(max_tapes) = &
          (/(' ',ni=1,max_tapes)/)
  ! namelist: per tape type1d
  character(len=max_namlen) , public :: hist_type1d_pertape(max_tapes) = &
          (/(' ',ni=1,max_tapes)/)
#ifdef LCH4
  ! namelist: write CH4 extra diagnostic output
  logical , public :: hist_wrtch4diag = .false.
#endif

  ! namelist-equivalence list of fields to add
  character(len=max_namlen+2) , public :: fincl(max_flds,max_tapes)

  ! namelist: list of fields to add
  character(len=max_namlen+2) , public :: hist_fincl1(max_flds) = ' '
  ! namelist: list of fields to add
  character(len=max_namlen+2) , public :: hist_fincl2(max_flds) = ' '
  ! namelist: list of fields to add
  character(len=max_namlen+2) , public :: hist_fincl3(max_flds) = ' '
  ! namelist: list of fields to add
  character(len=max_namlen+2) , public :: hist_fincl4(max_flds) = ' '
  ! namelist: list of fields to add
  character(len=max_namlen+2) , public :: hist_fincl5(max_flds) = ' '
  ! namelist: list of fields to add
  character(len=max_namlen+2) , public :: hist_fincl6(max_flds) = ' '

  ! namelist-equivalence list of fields to remove
  character(len=max_namlen+2) , public :: fexcl(max_flds,max_tapes)

  ! namelist: list of fields to remove
  character(len=max_namlen+2) , public :: hist_fexcl1(max_flds) = ' '
  ! namelist: list of fields to remove
  character(len=max_namlen+2) , public :: hist_fexcl2(max_flds) = ' '
  ! namelist: list of fields to remove
  character(len=max_namlen+2) , public :: hist_fexcl3(max_flds) = ' '
  ! namelist: list of fields to remove
  character(len=max_namlen+2) , public :: hist_fexcl4(max_flds) = ' '
  ! namelist: list of fields to remove
  character(len=max_namlen+2) , public :: hist_fexcl5(max_flds) = ' '
  ! namelist: list of fields to remove
  character(len=max_namlen+2) , public :: hist_fexcl6(max_flds) = ' '
  !
  ! Restart
  !
  logical , private :: if_disphist(max_tapes)   ! true => save history file

  ! Add a 1d single-level field to the master field list
  public :: hist_addfld1d
  ! Add a 2d multi-level field to the master field list
  public :: hist_addfld2d
  ! Add a 2d subscript dimension
  public :: hist_add_subscript
  ! Print summary of master field list
  public :: hist_printflds
  ! Initialize history file handler for initial or continue run
  public :: hist_htapes_build
  ! Updates history buffer for all fields and tapes
  public :: hist_update_hbuf
  ! Write history tape(s)
  public :: hist_htapes_wrapup
  ! Read/write history file restart data
  public :: hist_restart_ncd
  ! Define the contents of each history file based on namelist
  public :: htapes_fieldlist

  ! Add a field to a history file default "on" list
  private :: masterlist_make_active
  ! Add a field to the master field list
  private :: masterlist_addfld
  ! Override default history tape contents for specific tape
  private :: masterlist_change_timeavg
  ! Add a field to the active list for a history tape
  private :: htape_addfld
  ! Define contents of history file t
  private :: htape_create
  ! Write time constant values to history tape
  private :: htape_timeconst
  ! Write time constant 3D values to primary history tape
  private :: htape_timeconst3D
  ! Normalize history file fields by number of accumulations
  private :: hfields_normalize
  ! Zero out accumulation and hsitory buffers for a tape
  private :: hfields_zero
  ! Write a variable to a history tape
  private :: hfields_write
  ! Define/output 1d subgrid info if appropriate
  private :: hfields_1dinfo
  ! Updates history buffer for specific field and tape
  private :: hist_update_hbuf_field_1d
  ! Updates history buffer for specific field and tape
  private :: hist_update_hbuf_field_2d
  ! Find index of field in exclude list
  private :: list_index
  ! Determine history dataset filenames
  private :: set_hist_filename
  ! Retrieve name portion of input "inname"
  private :: getname
  ! Retrieve flag
  private :: getflag
  ! Track data pointer indices
  private :: pointer_index
  ! The max number of fields on any tape
  private :: max_nFields

  ! max chars for char variables
  integer(ik4) , parameter :: max_chars = 256
  !
  ! Subscript dimensions
  !
  ! max number of subscripts
  integer(ik4) , parameter :: max_subs = 100
  ! actual number of subscripts
  integer(ik4) :: num_subs = 0
  ! name of subscript
  character(len=32) :: subs_name(max_subs)
  ! dimension of subscript
  integer(ik4) :: subs_dim(max_subs)
  !
  ! Derived types
  !
  type field_info
    ! field name
    character(len=max_namlen) :: name
    ! long name
    character(len=max_chars) :: long_name
    ! units
    character(len=max_chars) :: units
    ! clm pointer first dimension type from clmtype (nameg, etc)
    character(len=8) :: type1d
    ! hbuf first dimension type from clmtype (nameg, etc)
    character(len=8) :: type1d_out
    ! hbuf second dimension type ["levgrnd","levlak",
    !                             "numrad","subname(n)"]
    character(len=8) :: type2d
    ! on-node 1d clm pointer start index
    integer(ik4) :: beg1d
    ! on-node 1d clm pointer end index
    integer(ik4) :: end1d
    ! size of clm pointer first dimension (all nodes)
    integer(ik4) :: num1d
    ! on-node 1d hbuf pointer start index
    integer(ik4) :: beg1d_out
    ! on-node 1d hbuf pointer end index
    integer(ik4) :: end1d_out
    ! size of hbuf first dimension (all nodes)
    integer(ik4) :: num1d_out
    ! size of hbuf second dimension (e.g. number of vertical levels)
    integer(ik4) :: num2d
    ! history pointer index
    integer(ik4) :: hpindex
    ! scale factor when averaging pft to column
    character(len=8) :: p2c_scale_type
    ! scale factor when averaging column to landunit
    character(len=8) :: c2l_scale_type
    ! scale factor when averaging landunit to gridcell
    character(len=8) :: l2g_scale_type
    type(subgrid_type) , pointer :: gcomm
  end type field_info

  type master_entry
    ! field information
    type (field_info) :: field
    ! active/inactive flag
    logical :: actflag(max_tapes)
    ! time averaging flag ("X","A","M" or "I",)
    character(len=1) :: avgflag(max_tapes)
  end type master_entry

  type history_entry
    ! field information
    type (field_info) :: field
    ! time averaging flag
    character(len=1) :: avgflag
    ! history buffer (dimensions: dim1d x num2d)
    real(rk8) , pointer :: hbuf(:,:)
    ! accumulation counter (dimensions: dim1d x num2d)
    integer(ik4) , pointer :: nacs(:,:)
  end type history_entry

  type history_tape
    ! number of active fields on tape
    integer(ik4) :: nflds
    ! current number of time samples on tape
    integer(ik4) :: ntimes
    ! maximum number of time samples per tape
    integer(ik4) :: mfilt
    ! number of time samples per tape
    integer(ik4) :: nhtfrq
    ! netcdf output precision
    integer(ik4) :: ncprec
    ! true => do xy average for all fields
    logical :: dov2xy
    ! true => current time step is end of history interval
    logical :: is_endhist
    ! time at beginning of history averaging interval
    real(rk8) :: begtime
    ! array of active history tape entries
    type(history_entry) :: hlist(max_flds)
  end type history_tape

  integer(ik4) , dimension(max_tapes) :: tapes_ntimes
  integer(ik4) , dimension(max_tapes) :: tapes_mfilt

  type clmpoint_rs  ! Pointer to real scalar data (1D)
    real(rk8) , pointer :: ptr(:)
  end type clmpoint_rs

  type clmpoint_ra  ! Pointer to real array data (2D)
    real(rk8) , pointer :: ptr(:,:)
  end type clmpoint_ra
  !
  ! Pointers into clmtype arrays
  !
  ! Maximum number of fields to track
  integer(ik4) , parameter :: max_mapflds = 2500
  type(clmpoint_rs) :: clmptr_rs(max_mapflds) ! Real scalar data (1D)
  type(clmpoint_ra) :: clmptr_ra(max_mapflds) ! Real array data (2D)
  !
  ! Master list: an array of master_entry entities
  !
  type(master_entry) :: masterlist(max_flds)  ! master field list
  !
  ! History tape: an array of history_tape entities (only active fields)
  !
  type(history_tape) :: tape(max_tapes)       ! array history tapes
  !
  ! Namelist input
  !
  ! Counters
  !
  ! number of fields in master field list
  integer(ik4) :: nfmaster = 0
  !
  ! Other variables
  !
  ! local history file names
  character(len=max_chars) :: locfnh(max_tapes)
  ! local history restart file names
  character(len=max_chars) :: locfnhr(max_tapes)
  ! flag indicates history contents have been defined
  logical :: htapes_defined = .false.
  !
  ! NetCDF  Id's
  !
  ! file ids
  type(clm_filetype) :: nfid(max_tapes)
  ! file ids for history restart files
  type(clm_filetype) :: ncid_hist(max_tapes)
  ! string dimension id
  integer(ik4) :: strlen_dimid

  contains
  !
  ! Print summary of master field list.
  !
  subroutine hist_printflds()
    implicit none
    integer(ik4) :: nf
    character(len=*) ,parameter :: subname = 'CLM_hist_printflds'

    ! This is just bloating the screen
    if ( myid == italk ) then
      write(stdout,*) trim(subname),' : number of master fields = ',nfmaster
    end if

    if ( debug_level > 3 ) then
      if ( myid == italk ) then
        write(stdout,*)' ******* MASTER FIELD LIST *******'
        do nf = 1 , nfmaster
          write(stdout,9000) nf, masterlist(nf)%field%name, &
                  masterlist(nf)%field%units
9000      format (i5,1x,a32,1x,a16)
        end do
      end if
    end if
  end subroutine hist_printflds
  !
  ! Add a field to the master field list. Put input arguments of
  ! field name, units, number of levels, averaging flag, and long name
  ! into a type entry in the global master field list (masterlist).
  !
  subroutine masterlist_addfld (fname, type1d, type1d_out, &
        type2d, num2d, units, avgflag, long_name, hpindex, &
        p2c_scale_type, c2l_scale_type, l2g_scale_type)
    implicit none
    character(len=*) , intent(in) :: fname        ! field name
    character(len=*) , intent(in) :: type1d       ! 1d data type
    character(len=*) , intent(in) :: type1d_out   ! 1d output type
    character(len=*) , intent(in) :: type2d       ! 2d output type
    ! size of second dimension (e.g. number of vertical levels)
    integer(ik4) , intent(in) :: num2d
    character(len=*) , intent(in) :: units        ! units of field
    character(len=1) , intent(in) :: avgflag      ! time averaging flag
    character(len=*) , intent(in) :: long_name    ! long name of field
    ! clmtype index for history buffer output
    integer(ik4) , intent(in) :: hpindex
    ! scale type for subgrid averaging of pfts to column
    character(len=*) , intent(in) :: p2c_scale_type
    ! scale type for subgrid averaging of columns to landunits
    character(len=*) , intent(in) :: c2l_scale_type
    ! scale type for subgrid averaging of landunits to gridcells
    character(len=*) , intent(in) :: l2g_scale_type
    integer(ik4) :: n            ! loop index
    integer(ik4) :: f            ! masterlist index
    integer(ik4) :: begp , endp  ! per-proc beginning and ending pft indices
    integer(ik4) :: begc , endc  ! per-proc beginning and ending column indices
    integer(ik4) :: begl , endl  ! per-proc beginning and ending ldunit indices
    integer(ik4) :: begg , endg  ! per-proc gridcell ending gridcell indices
    integer(ik4) :: numa         ! total number of atm cells across all proc
    integer(ik4) :: numg         ! total number of gridcells across all proc
    integer(ik4) :: numl         ! total number of landunits across all proc
    integer(ik4) :: numc         ! total number of columns across all proc
    integer(ik4) :: nump         ! total number of pfts across all proc
    character(len=*) , parameter :: subname = 'masterlist_addfld'

    ! Determine bounds

    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)
    call get_proc_global(numg,numl,numc,nump)
    ! Ensure that new field is not all blanks

    if ( fname == ' ' ) then
      write(stderr,*) trim(subname),' ERROR: blank field name not allowed'
      call fatal(__FILE__,__LINE__,'clm now stopping.')
    end if

    ! Ensure that new field doesn't already exist

    do n = 1 , nfmaster
      if ( masterlist(n)%field%name == fname ) then
        write(stderr,*) trim(subname),' ERROR:', fname, ' already on list'
        call fatal(__FILE__,__LINE__,'clm now stopping.')
      end if
    end do

    ! Increase number of fields on master field list

    nfmaster = nfmaster + 1
    f = nfmaster

    ! Check number of fields in master list against maximum number
    ! for master list

    if ( nfmaster > max_flds ) then
      write(stderr,*) trim(subname), &
            ' ERROR: too many fields for primary history file ', &
            '-- max_flds,nfmaster=', max_flds, nfmaster
      call fatal(__FILE__,__LINE__,'clm now stopping.')
    end if

    ! Add field to master list

    masterlist(f)%field%name           = fname
    masterlist(f)%field%long_name      = long_name
    masterlist(f)%field%units          = units
    masterlist(f)%field%type1d         = type1d
    masterlist(f)%field%type1d_out     = type1d_out
    masterlist(f)%field%type2d         = type2d
    masterlist(f)%field%num2d          = num2d
    masterlist(f)%field%hpindex        = hpindex
    masterlist(f)%field%p2c_scale_type = p2c_scale_type
    masterlist(f)%field%c2l_scale_type = c2l_scale_type
    masterlist(f)%field%l2g_scale_type = l2g_scale_type

    select case (type1d)
      case (grlnd)
        masterlist(f)%field%beg1d = begg
        masterlist(f)%field%end1d = endg
        masterlist(f)%field%num1d = numg
        masterlist(f)%field%gcomm => gcomm_gridcell
      case (nameg)
        masterlist(f)%field%beg1d = begg
        masterlist(f)%field%end1d = endg
        masterlist(f)%field%num1d = numg
        masterlist(f)%field%gcomm => gcomm_gridcell
      case (namel)
        masterlist(f)%field%beg1d = begl
        masterlist(f)%field%end1d = endl
        masterlist(f)%field%num1d = numl
        masterlist(f)%field%gcomm => gcomm_landunit
      case (namec)
        masterlist(f)%field%beg1d = begc
        masterlist(f)%field%end1d = endc
        masterlist(f)%field%num1d = numc
        masterlist(f)%field%gcomm => gcomm_column
      case (namep)
        masterlist(f)%field%beg1d = begp
        masterlist(f)%field%end1d = endp
        masterlist(f)%field%num1d = nump
        masterlist(f)%field%gcomm => gcomm_pft
      case default
        write(stderr,*) trim(subname),' ERROR: unknown 1d output type= ',type1d
        call fatal(__FILE__,__LINE__,'clm now stopping.')
    end select

    ! The following two fields are used only in master field list,
    ! NOT in the runtime active field list
    ! ALL FIELDS IN THE MASTER LIST ARE INITIALIZED WITH THE ACTIVE
    ! FLAG SET TO FALSE
    masterlist(f)%avgflag(:) = avgflag
    masterlist(f)%actflag(:) = .false.
  end subroutine masterlist_addfld
  !
  ! Initialize history file for initial or continuation run.  For example,
  ! on an initial run, this routine initializes ``ntapes'' history files.
  ! On a restart run, this routine only initializes history files declared
  ! beyond what existed on the previous run.  Files which already existed on
  ! the previous run have already been initialized (i.e. named and opened)
  ! in routine restart\_history.  Loop over tapes and fields per tape setting
  ! appropriate variables and calling appropriate routines
  !
  subroutine hist_htapes_build ()
    implicit none
    integer(ik4) :: i            ! index
    integer(ik4) :: ier          ! error code
    integer(ik4) :: t , f        ! tape, field indices
    integer(ik4) :: day , sec    ! day and seconds from base date
    character(len=*) , parameter :: subname = 'hist_htapes_build'

    if ( myid == italk ) then
      write(stdout,*)  trim(subname),' Initializing clm2 history files'
      write(stdout,'(72a1)') ("-",i=1,60)
    end if

    ! Define field list information for all history files.
    ! Update ntapes to reflect number of active history files

    call htapes_fieldlist()

    ! Determine if gridcell (xy) averaging is done for all fields on tape

    do t = 1 , ntapes
      tape(t)%dov2xy = hist_dov2xy(t)
      if ( myid == italk ) then
        write(stdout,*) trim(subname),' hist tape = ',t,&
              ' written with dov2xy= ',tape(t)%dov2xy
      end if
    end do

    ! Set number of time samples in each history file and
    ! Note - the following entries will be overwritten by history restart
    ! Note - with netcdf, only 1 (ncd_double) and 2 (ncd_float) are allowed

    do t = 1 , ntapes
      tape(t)%ntimes = 0
      tape(t)%dov2xy = hist_dov2xy(t)
      tape(t)%nhtfrq = hist_nhtfrq(t)
      tape(t)%mfilt = hist_mfilt(t)
      if ( hist_ndens(t) == 1 ) then
        tape(t)%ncprec = clmvar_double
      else
        tape(t)%ncprec = clmvar_real
      end if
    end do

    ! Set time of beginning of current averaging interval
    ! First etermine elapased time since reference date

    call get_prev_time(day, sec)
    do t = 1 , ntapes
      tape(t)%begtime = day + sec/secspday
    end do

    if (myid == italk) then
      write(stdout,*) trim(subname), &
              ' Successfully initialized clm2 history files'
      write(stdout,'(72a1)') ("-",i=1,60)
    end if
  end subroutine hist_htapes_build
  !
  ! Add a field to the default ``on'' list for a given history file.
  ! Also change the default time averaging flag if requested.
  !
  subroutine masterlist_make_active (name, tape_index, avgflag)
    implicit none
    character(len=*) , intent(in) :: name          ! field name
    integer(ik4) , intent(in) :: tape_index        ! history tape index
    character(len=1) , intent(in) , optional :: avgflag  ! time averaging flag
    integer(ik4) :: f     ! field index
    logical :: found      ! flag indicates field found in masterlist
    character(len=*) , parameter :: subname = 'masterlist_make_active'

    ! Check validity of input arguments

    if (tape_index > max_tapes) then
      write(stderr,*) trim(subname), &
              ' ERROR: tape index=', tape_index, ' is too big'
      call fatal(__FILE__,__LINE__,'clm now stopping.')
    end if

    if ( present(avgflag) ) then
      if ( avgflag /= ' ' .and. &
           avgflag /= 'A' .and. avgflag /= 'I' .and. &
           avgflag /= 'X' .and. avgflag /= 'M') then
        write(stderr,*) trim(subname), &
                ' ERROR: unknown averaging flag=', avgflag
        call fatal(__FILE__,__LINE__,'clm now stopping.')
      end if
    end if

    ! Look through master list for input field name.
    ! When found, set active flag for that tape to true.
    ! Also reset averaging flag if told to use other than default.

    found = .false.
    do f = 1 , nfmaster
      if ( trim(name) == trim(masterlist(f)%field%name) ) then
        masterlist(f)%actflag(tape_index) = .true.
        if ( present(avgflag) ) then
          if ( avgflag /= ' ' ) masterlist(f)%avgflag(tape_index) = avgflag
        end if
        found = .true.
        exit
      end if
    end do
    if ( .not. found ) then
      write(stderr,*) trim(subname),' ERROR: field=', name, ' not found'
      call fatal(__FILE__,__LINE__,'clm now stopping.')
    end if
  end subroutine masterlist_make_active
  !
  ! Override default history tape contents for a specific tape.
  ! Copy the flag into the master field list.
  !
  subroutine masterlist_change_timeavg (t)
    implicit none
    integer(ik4) , intent(in) :: t  ! history tape index
    integer(ik4) :: f               ! field index
    character(len=1) :: avgflag     ! lcl equiv of hist_avgflag_pertape(t)
    character(len=*) ,parameter :: subname = 'masterlist_change_timeavg'

    avgflag = hist_avgflag_pertape(t)

    do f = 1 , nfmaster
      select case (avgflag)
        case ('A')
          masterlist(f)%avgflag(t) = avgflag
        case ('I')
          masterlist(f)%avgflag(t) = avgflag
        case ('X')
          masterlist(f)%avgflag(t) = avgflag
        case ('M')
          masterlist(f)%avgflag(t) = avgflag
        case default
          write(stderr,*) trim(subname),' ERROR: unknown avgflag=',avgflag
          call fatal(__FILE__,__LINE__,'clm now stopping.')
      end select
    end do
  end subroutine masterlist_change_timeavg
  !
  ! Define the contents of each history file based on namelist
  ! input for initial, and restart data if a restart run.
  ! Use arrays fincl and fexcl to modify default history tape contents.
  ! Then sort the result alphanumerically.
  !
  subroutine htapes_fieldlist()
    implicit none
    integer(ik4) :: t , f   ! tape, field indices
    integer(ik4) :: ff      ! index into include, exclude and fprec list
    ! field name portion of fincl (i.e. no avgflag separator)
    character(len=max_namlen) :: name
    character(len=max_namlen) :: mastername ! name from masterlist field
    character(len=1) :: avgflag    ! averaging flag
    character(len=1) :: prec_acc   ! history buffer precision flag
    character(len=1) :: prec_wrt   ! history buffer write precision flag
    type (history_entry) :: tmp    ! temporary used for swapping
    character(len=*) , parameter :: subname = 'htapes_fieldlist'

    ! Override averaging flag for all fields on a particular tape
    ! if namelist input so specifies

    do t = 1 , max_tapes
      if ( hist_avgflag_pertape(t) /= ' ' ) then
        call masterlist_change_timeavg (t)
      end if
    end do

    fincl(:,1) = hist_fincl1(:)
    fincl(:,2) = hist_fincl2(:)
    fincl(:,3) = hist_fincl3(:)
    fincl(:,4) = hist_fincl4(:)
    fincl(:,5) = hist_fincl5(:)
    fincl(:,6) = hist_fincl6(:)

    fexcl(:,1) = hist_fexcl1(:)
    fexcl(:,2) = hist_fexcl2(:)
    fexcl(:,3) = hist_fexcl3(:)
    fexcl(:,4) = hist_fexcl4(:)
    fexcl(:,5) = hist_fexcl5(:)
    fexcl(:,6) = hist_fexcl6(:)

    ! First ensure contents of fincl and fexcl are valid names

    do t = 1 , max_tapes
      f = 1
      do while (f < max_flds .and. fincl(f,t) /= ' ')
        name = getname (fincl(f,t))
        do ff = 1 , nfmaster
          mastername = masterlist(ff)%field%name
          if ( name == mastername ) exit
        end do
        if ( name /= mastername ) then
          write(stderr,*) trim(subname), &
                 ' ERROR: ', trim(name), ' in fincl(', f, ') ',&
                 'for history tape ',t,' not found'
          call fatal(__FILE__,__LINE__,'clm now stopping.')
        end if
        f = f + 1
      end do

      f = 1
      do while (f < max_flds .and. fexcl(f,t) /= ' ')
        do ff = 1 , nfmaster
          mastername = masterlist(ff)%field%name
          if ( fexcl(f,t) == mastername ) exit
        end do
        if ( fexcl(f,t) /= mastername ) then
          write(stderr,*) trim(subname), &
                 ' ERROR: ', fexcl(f,t), ' in fexcl(', f, ') ', &
                 'for history tape ',t,' not found'
          call fatal(__FILE__,__LINE__,'clm now stopping.')
        end if
        f = f + 1
      end do
    end do

    tape(:)%nflds = 0
    do t = 1 , max_tapes
      ! Loop through the masterlist set of field names and determine if
      ! any of those are in the FINCL or FEXCL arrays
      ! The call to list_index determines the index in the FINCL or FEXCL
      ! arrays that the masterlist field corresponds to
      ! Add the field to the tape if specified via namelist
      ! (FINCL[1-max_tapes]), or if it is on by default and was not
      ! excluded via namelist (FEXCL[1-max_tapes]).

      do f = 1 , nfmaster
        mastername = masterlist(f)%field%name
        call list_index (fincl(1,t), mastername, ff)
        if ( ff > 0 ) then
          ! if field is in include list, ff > 0 and htape_addfld
          ! will not be called for field
          avgflag = getflag (fincl(ff,t))
          call htape_addfld (t, f, avgflag)
        else if ( .not. hist_empty_htapes ) then
          ! find index of field in exclude list
          call list_index (fexcl(1,t), mastername, ff)
          ! if field is in exclude list, ff > 0 and htape_addfld
          ! will not be called for field
          ! if field is not in exclude list, ff =0 and htape_addfld
          ! will be called for field (note that htape_addfld will be
          ! called below only if field is not in exclude list OR in
          ! include list
          if ( ff == 0 .and. masterlist(f)%actflag(t) ) then
            call htape_addfld (t, f, ' ')
          end if
        end if
      end do

      ! Specification of tape contents now complete.
      ! Sort each list of active entries

      do f = tape(t)%nflds-1 , 1 , -1
        do ff = 1 , f
          if ( tape(t)%hlist(ff)%field%name > &
               tape(t)%hlist(ff+1)%field%name ) then
            tmp = tape(t)%hlist(ff)
            tape(t)%hlist(ff  ) = tape(t)%hlist(ff+1)
            tape(t)%hlist(ff+1) = tmp
          else if (tape(t)%hlist(ff)%field%name == &
                   tape(t)%hlist(ff+1)%field%name ) then
            write(stderr,*) trim(subname),' ERROR: Duplicate field ', &
                   tape(t)%hlist(ff)%field%name, &
                   't,ff,name=',t,ff,tape(t)%hlist(ff+1)%field%name
            call fatal(__FILE__,__LINE__,'clm now stopping.')
          end if
        end do
      end do

      if ( debug_level > 3 ) then
        if ( myid == italk ) then
          if ( tape(t)%nflds > 0 ) then
            write(stdout,*) trim(subname), &
                    ' : Included fields tape ',t,'=',tape(t)%nflds
          end if
          do f = 1 , tape(t)%nflds
            write(stdout,*) f,' ',tape(t)%hlist(f)%field%name, &
                  tape(t)%hlist(f)%field%num2d,' ',tape(t)%hlist(f)%avgflag
          end do
        end if
      end if
    end do

    ! Determine total number of active history tapes

    ntapes = 0
    do t = max_tapes , 1 , -1
      if ( tape(t)%nflds > 0 ) then
        ntapes = t
        exit
      end if
    end do

    ! Ensure there are no "holes" in tape specification, i.e. empty tapes.
    ! Enabling holes should not be difficult if necessary.

    do t = 1 , ntapes
      if (tape(t)%nflds  ==  0) then
        write(stderr,*) trim(subname),' ERROR: Tape ',t,' is empty'
        call fatal(__FILE__,__LINE__,'clm now stopping.')
      end if
    end do

    ! Check that the number of history files declared does not exceed
    ! the maximum allowed.

    if (ntapes > max_tapes) then
      write(stderr,*) trim(subname), &
               ' ERROR: Too many history files declared, max_tapes=',max_tapes
      call fatal(__FILE__,__LINE__,'clm now stopping.')
    end if

    ! Change 1d output per tape output flag if requested - only for history
    ! tapes where 2d xy averaging is not enabled

    do t = 1 , ntapes
      if ( hist_type1d_pertape(t) /= ' ' .and. (.not. hist_dov2xy(t)) ) then
        select case (trim(hist_type1d_pertape(t)))
          case ('PFTS','COLS', 'LAND', 'GRID')
            if ( myid == italk ) &
            write(stdout,*) 'history tape ',t, &
              ' will have 1d output type of ',hist_type1d_pertape(t)
          case default
            write(stderr,*) trim(subname), &
                     ' ERROR: unknown namelist type1d per tape=', &
                      hist_type1d_pertape(t)
            call fatal(__FILE__,__LINE__,'clm now stopping.')
        end select
      end if
    end do

    if ( myid == italk ) then
      write(stdout,*) 'There will be a total of ',ntapes,' history tapes'
      do t = 1 , ntapes
        write(stdout,*)
        if ( hist_nhtfrq(t) == 0 ) then
          write(stdout,*) 'History tape ',t,' write frequency is MONTHLY'
        else
          write(stdout,*)'History tape ',t,' write frequency = ',hist_nhtfrq(t)
        end if
        if ( hist_dov2xy(t) ) then
          write(stdout,*)'All fields on history tape ',t,' are grid averaged'
        else
          write(stdout,*)'All fields on history tape ',t, &
                  ' are not grid averaged'
        end if
        write(stdout,*)'Number of time samples on history tape ',t, &
                ' is ',hist_mfilt(t)
        write(stdout,*)'Output precision on history tape ',t,'=',hist_ndens(t)
        write(stdout,*)
      end do
    end if

    ! Set flag indicating h-tape contents are now defined
    ! (needed by masterlist_addfld)
    htapes_defined = .true.
  end subroutine htapes_fieldlist
  !
  ! Add a field to the active list for a history tape. Copy the data from
  ! the master field list to the active list for the tape.
  !
  subroutine htape_addfld (t, f, avgflag)
    implicit none
    integer(ik4) , intent(in) :: t  ! history tape index
    integer(ik4) , intent(in) :: f  ! field index from master field list
    character(len=1) , intent(in) :: avgflag  ! time averaging flag
    integer(ik4) :: n                    ! field index on defined tape
    character(len=8) :: type1d      ! clm pointer 1d type
    character(len=8) :: type1d_out  ! history buffer 1d type
    integer(ik4) :: begp , endp ! per-proc beginning and ending pft indices
    integer(ik4) :: begc , endc ! per-proc beginning and ending column indices
    integer(ik4) :: begl , endl ! per-proc beginning and ending landunit indices
    integer(ik4) :: begg , endg ! per-proc gridcell ending gridcell indices
    integer(ik4) :: numa        ! total number of atm cells across all proc
    integer(ik4) :: numg        ! total number of gridcells across all proc
    integer(ik4) :: numl        ! total number of landunits across all proc
    integer(ik4) :: numc        ! total number of columns across all proc
    integer(ik4) :: nump        ! total number of pfts across all proc
    ! size of second dimension (e.g. .number of vertical levels)
    integer(ik4) :: num2d
    ! history output per-proc 1d beginning and ending indices
    integer(ik4) :: beg1d_out , end1d_out
    integer(ik4) :: num1d_out  ! history output 1d size
    character(len=*) , parameter :: subname = 'htape_addfld'

    ! Ensure that it is not to late to add a field to the history tape

    if ( htapes_defined ) then
      write(stderr,*) trim(subname),' ERROR: attempt to add field ', &
           masterlist(f)%field%name, ' after history files are set'
      call fatal(__FILE__,__LINE__,'clm now stopping.')
    end if

    tape(t)%nflds = tape(t)%nflds + 1
    n = tape(t)%nflds

    ! Copy field information

    tape(t)%hlist(n)%field = masterlist(f)%field

    ! Determine bounds

    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)
    call get_proc_global(numg,numl,numc,nump)

    ! Modify type1d_out if necessary

    if ( hist_dov2xy(t) ) then

      ! If xy output averaging is requested, set output 1d type to grlnd
      ! ***NOTE- the following logic is what permits non lat/lon grids to
      ! be written to clm history file

      type1d = tape(t)%hlist(n)%field%type1d

      if ( type1d == nameg .or. &
           type1d == namel .or. &
           type1d == namec .or. &
           type1d == namep) then
        tape(t)%hlist(n)%field%type1d_out = grlnd(1:8)
        tape(t)%hlist(n)%field%gcomm => gcomm_gridcell
      end if
      if ( type1d == grlnd ) then
        tape(t)%hlist(n)%field%type1d_out = grlnd(1:8)
        tape(t)%hlist(n)%field%gcomm => gcomm_gridcell
      end if

    else if ( hist_type1d_pertape(t) /= ' ' ) then

      ! Set output 1d type  based on namelist setting of  hist_type1d_pertape
      ! Only applies to tapes when xy output is not required

      type1d = tape(t)%hlist(n)%field%type1d

      select case (trim(hist_type1d_pertape(t)))
        case('GRID')
          tape(t)%hlist(n)%field%type1d_out = nameg(1:8)
          tape(t)%hlist(n)%field%gcomm => gcomm_gridcell
        case('LAND')
          tape(t)%hlist(n)%field%type1d_out = namel(1:8)
          tape(t)%hlist(n)%field%gcomm => gcomm_landunit
        case('COLS')
          tape(t)%hlist(n)%field%type1d_out = namec(1:8)
          tape(t)%hlist(n)%field%gcomm => gcomm_column
        case ('PFTS')
          tape(t)%hlist(n)%field%type1d_out = namep(1:8)
          tape(t)%hlist(n)%field%gcomm => gcomm_pft
        case default
          write(stderr,*) trim(subname), &
                  ' ERROR: unknown input hist_type1d_pertape= ', &
                  hist_type1d_pertape(t)
          call fatal(__FILE__,__LINE__,'clm now stopping.')
      end select
    end if

    ! Determine output 1d dimensions

    type1d_out = tape(t)%hlist(n)%field%type1d_out
    if ( type1d_out == grlnd ) then
      beg1d_out = begg
      end1d_out = endg
      num1d_out = numg
    else if ( type1d_out == nameg ) then
      beg1d_out = begg
      end1d_out = endg
      num1d_out = numg
    else if ( type1d_out == namel ) then
      beg1d_out = begl
      end1d_out = endl
      num1d_out = numl
    else if ( type1d_out == namec ) then
      beg1d_out = begc
      end1d_out = endc
      num1d_out = numc
    else if ( type1d_out == namep ) then
      beg1d_out = begp
      end1d_out = endp
      num1d_out = nump
    else
      write(stderr,*) trim(subname), &
              ' ERROR: incorrect value of type1d_out= ',type1d_out
      call fatal(__FILE__,__LINE__,'clm now stopping.')
    end if

    tape(t)%hlist(n)%field%beg1d_out = beg1d_out
    tape(t)%hlist(n)%field%end1d_out = end1d_out
    tape(t)%hlist(n)%field%num1d_out = num1d_out

    ! Alloccate and initialize history buffer and related info

    num2d = tape(t)%hlist(n)%field%num2d
    allocate (tape(t)%hlist(n)%hbuf(beg1d_out:end1d_out,num2d))
    allocate (tape(t)%hlist(n)%nacs(beg1d_out:end1d_out,num2d))
    tape(t)%hlist(n)%hbuf(:,:) = 0.D0
    tape(t)%hlist(n)%nacs(:,:) = 0

    ! Set time averaging flag based on masterlist setting or
    ! override the default averaging flag with namelist setting

    select case (avgflag)
      case (' ')
        tape(t)%hlist(n)%avgflag = masterlist(f)%avgflag(t)
      case ('A','I','X','M')
        tape(t)%hlist(n)%avgflag = avgflag
      case default
        write(stderr,*) trim(subname),' ERROR: unknown avgflag=', avgflag
        call fatal(__FILE__,__LINE__,'clm now stopping.')
    end select
  end subroutine htape_addfld
  !
  ! Accumulate (or take min, max, etc. as appropriate) input field
  ! into its history buffer for appropriate tapes.
  !
  subroutine hist_update_hbuf()
    implicit none
    integer(ik4) :: t           ! tape index
    integer(ik4) :: f           ! field index
    integer(ik4) :: begp , endp ! per-proc beginning and ending pft indices
    integer(ik4) :: begc , endc ! per-proc beginning and ending column indices
    integer(ik4) :: begl , endl ! per-proc beginning and ending landunit indices
    integer(ik4) :: begg , endg ! per-proc gridcell ending gridcell indices
    ! size of second dimension (e.g. number of vertical levels)
    integer(ik4) :: num2d
    character(len=*) , parameter :: subname = 'hist_update_hbuf'

    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

    do t = 1 , ntapes
      do f = 1 , tape(t)%nflds
        num2d = tape(t)%hlist(f)%field%num2d
        if ( num2d == 1 ) then
          call hist_update_hbuf_field_1d(t,f,begp,endp,begc,endc, &
                                         begl,endl,begg,endg)
        else
          call hist_update_hbuf_field_2d(t,f,begp,endp,begc,endc, &
                                         begl,endl,begg,endg,num2d)
        end if
      end do
    end do
  end subroutine hist_update_hbuf
  !
  ! Accumulate (or take min, max, etc. as appropriate) input field
  ! into its history buffer for appropriate tapes.
  !
  subroutine hist_update_hbuf_field_1d (t,f,begp,endp,begc,endc, &
                                            begl,endl,begg,endg)
    implicit none
    integer(ik4) , intent(in) :: t   ! tape index
    integer(ik4) , intent(in) :: f   ! field index
    ! per-proc beginning and ending pft indices
    integer(ik4) , intent(in) :: begp , endp
    ! per-proc beginning and ending column indices
    integer(ik4) , intent(in) :: begc , endc
    ! per-proc beginning and ending landunit indices
    integer(ik4) , intent(in) :: begl , endl
    ! per-proc gridcell ending gridcell indices
    integer(ik4) , intent(in) :: begg , endg
    integer(ik4) :: hpindex   ! history pointer index
    integer(ik4) :: k         ! gridcell, landunit, column or pft index
    integer(ik4) :: beg1d , end1d    ! beginning and ending indices
    ! true => check 'active' flag of each point
    ! (this refers to a point being active, NOT a history field being active)
    logical :: check_active
    logical :: valid       ! true => history operation is valid
    logical :: map2gcell   ! true => map clm pointer field to gridcell
    ! 1d clm pointerr type   ["gridcell","landunit","column","pft"]
    character(len=8) :: type1d
    ! 1d history buffer type ["gridcell","landunit","column","pft"]
    character(len=8) :: type1d_out
    ! time averaging flag
    character(len=1) :: avgflag
    ! scale type for subgrid averaging of pfts to column
    character(len=8) :: p2c_scale_type
    ! scale type for subgrid averaging of columns to landunits
    character(len=8) :: c2l_scale_type
    ! scale type for subgrid averaging of landunits to gridcells
    character(len=8) :: l2g_scale_type
    real(rk8) , pointer :: hbuf(:,:)      ! history buffer
    integer(ik4) , pointer :: nacs(:,:)   ! accumulation counter
    real(rk8) , pointer :: field(:)       ! clm 1d pointer field
    ! flag saying whether each point is active
    ! (used for type1d = landunit/column/pft)
    ! (this refers to a point being active, NOT a history field being active)
    logical , pointer :: active(:)
    ! gricell level field (used if mapping to gridcell is done)
    real(rk8) :: field_gcell(begg:endg)
    integer(ik4) :: j
    character(len=*) , parameter :: subname = 'hist_update_hbuf_field_1d'
    ! offset for mapping sliced subarray pointers when outputting variables
    ! in PFT/col vector form
    integer(ik4) :: k_offset

    avgflag        =  tape(t)%hlist(f)%avgflag
    nacs           => tape(t)%hlist(f)%nacs
    hbuf           => tape(t)%hlist(f)%hbuf
    beg1d          =  tape(t)%hlist(f)%field%beg1d
    end1d          =  tape(t)%hlist(f)%field%end1d
    type1d         =  tape(t)%hlist(f)%field%type1d
    type1d_out     =  tape(t)%hlist(f)%field%type1d_out
    p2c_scale_type =  tape(t)%hlist(f)%field%p2c_scale_type
    c2l_scale_type =  tape(t)%hlist(f)%field%c2l_scale_type
    l2g_scale_type =  tape(t)%hlist(f)%field%l2g_scale_type
    hpindex        =  tape(t)%hlist(f)%field%hpindex
    field          => clmptr_rs(hpindex)%ptr

    ! set variables to check weights when allocate all pfts

    map2gcell = .false.
    if ( type1d_out == nameg .or. type1d_out == grlnd ) then
      if (type1d == namep ) then
        call p2g(begp,endp,begc,endc,begl,endl,begg,endg,field,field_gcell, &
               p2c_scale_type,c2l_scale_type,l2g_scale_type)
        map2gcell = .true.
      else if ( type1d == namec ) then
        call c2g(begc,endc,begl,endl,begg,endg,field,field_gcell, &
               c2l_scale_type,l2g_scale_type)
        map2gcell = .true.
      else if ( type1d == namel ) then
        call l2g(begl,endl,begg,endg,field,field_gcell, &
               l2g_scale_type)
        map2gcell = .true.
      end if
    end if

    if ( map2gcell ) then  ! Map to gridcell
      ! note that in this case beg1d = begg and end1d=endg
      select case (avgflag)
        case ('I') ! Instantaneous
          do k = begg , endg
            if ( field_gcell(k) /= spval ) then
              hbuf(k,1) = field_gcell(k)
            else
              hbuf(k,1) = spval
            end if
            nacs(k,1) = 1
          end do
        case ('A') ! Time average
          do k = begg , endg
            if ( field_gcell(k) /= spval ) then
              if ( nacs(k,1) == 0 ) hbuf(k,1) = 0.D0
              hbuf(k,1) = hbuf(k,1) + field_gcell(k)
              nacs(k,1) = nacs(k,1) + 1
            else
              if ( nacs(k,1) == 0 ) hbuf(k,1) = spval
            end if
          end do
        case ('X') ! Maximum over time
          do k = begg , endg
            if ( field_gcell(k) /= spval ) then
              if ( nacs(k,1) == 0 ) hbuf(k,1) = -1.D50
              hbuf(k,1) = max( hbuf(k,1), field_gcell(k) )
            else
              hbuf(k,1) = spval
            end if
            nacs(k,1) = 1
          end do
        case ('M') ! Minimum over time
          do k = begg , endg
            if ( field_gcell(k) /= spval ) then
              if ( nacs(k,1) == 0 ) hbuf(k,1) = +1.D50
              hbuf(k,1) = min( hbuf(k,1), field_gcell(k) )
            else
              hbuf(k,1) = spval
            end if
            nacs(k,1) = 1
          end do
        case default
          write(stderr,*) trim(subname), &
                  ' ERROR: invalid time averaging flag ', avgflag
          call fatal(__FILE__,__LINE__,'clm now stopping.')
      end select
    else  ! Do not map to gridcell

      ! For data defined on the pft, col or landunit, we need to
      ! check if a point is active to determine whether that point
      ! should be assigned spval
      if ( type1d == namep ) then
        check_active = .true.
        active => clm3%g%l%c%p%active
      else if ( type1d == namec ) then
        check_active = .true.
        active => clm3%g%l%c%active
      else if ( type1d == namel ) then
        check_active = .true.
        active => clm3%g%l%active
      else
        check_active = .false.
      end if

      select case (avgflag)
        case ('I') ! Instantaneous
          do k = beg1d , end1d
            valid = .true.
            if ( check_active ) then
              if ( .not. active(k) ) valid = .false.
            end if
            if ( valid ) then
              if ( field(k) /= spval ) then
                hbuf(k,1) = field(k)
              else
                hbuf(k,1) = spval
              end if
            else
              hbuf(k,1) = spval
            end if
            nacs(k,1) = 1
          end do
        case ('A') ! Time average
          ! create mappings for array slice pointers (which go from
          ! 1 to size(field) rather than beg1d to end1d)
          if ( end1d .eq. ubound(field,1) ) then
            k_offset = 0
          else
            k_offset = 1 - beg1d
          end if
          do k = beg1d , end1d
            valid = .true.
            if ( check_active ) then
              if ( .not. active(k) ) valid = .false.
            end if
            if ( valid ) then
              if ( field(k+k_offset) /= spval ) then   ! add k_offset
                if ( nacs(k,1) == 0 ) hbuf(k,1) = 0.D0
                hbuf(k,1) = hbuf(k,1) + field(k+k_offset)   ! add k_offset
                nacs(k,1) = nacs(k,1) + 1
              else
                if ( nacs(k,1) == 0 ) hbuf(k,1) = spval
              end if
            else
              if ( nacs(k,1) == 0 ) hbuf(k,1) = spval
            end if
          end do
        case ('X') ! Maximum over time
          do k = beg1d , end1d
            valid = .true.
            if ( check_active ) then
              if ( .not. active(k) ) valid = .false.
            end if
            if ( valid ) then
              if ( field(k) /= spval ) then
                if ( nacs(k,1) == 0 ) hbuf(k,1) = -1.D50
                hbuf(k,1) = max( hbuf(k,1), field(k) )
              else
                if ( nacs(k,1) == 0 ) hbuf(k,1) = spval
              end if
            else
              if ( nacs(k,1) == 0 ) hbuf(k,1) = spval
            end if
            nacs(k,1) = 1
          end do
        case ('M') ! Minimum over time
          do k = beg1d , end1d
            valid = .true.
            if ( check_active ) then
              if ( .not. active(k) ) valid = .false.
            end if
            if ( valid ) then
              if ( field(k) /= spval ) then
                if ( nacs(k,1) == 0 ) hbuf(k,1) = +1.D50
                hbuf(k,1) = min( hbuf(k,1), field(k) )
              else
                if ( nacs(k,1) == 0 ) hbuf(k,1) = spval
              end if
            else
              if ( nacs(k,1) == 0 ) hbuf(k,1) = spval
            end if
            nacs(k,1) = 1
          end do
        case default
          write(stderr,*) trim(subname), &
                  ' ERROR: invalid time averaging flag ', avgflag
          call fatal(__FILE__,__LINE__,'clm now stopping.')
      end select
    end if
  end subroutine hist_update_hbuf_field_1d
  !
  ! Accumulate (or take min, max, etc. as appropriate) input field
  ! into its history buffer for appropriate tapes.
  !
  subroutine hist_update_hbuf_field_2d(t,f,begp,endp,begc,endc, &
                                           begl,endl,begg,endg,num2d)
    implicit none
    integer(ik4) , intent(in) :: t   ! tape index
    integer(ik4) , intent(in) :: f   ! field index
    ! per-proc beginning and ending pft indices
    integer(ik4) , intent(in) :: begp , endp
    ! per-proc beginning and ending column indices
    integer(ik4) , intent(in) :: begc , endc
    ! per-proc beginning and ending landunit indices
    integer(ik4) , intent(in) :: begl , endl
    ! per-proc gridcell ending gridcell indices
    integer(ik4) , intent(in) :: begg , endg
    ! size of second dimension
    integer(ik4) , intent(in) :: num2d
    integer(ik4) :: hpindex   ! history pointer index
    integer(ik4) :: k         ! gridcell, landunit, column or pft index
    integer(ik4) :: j         ! level index
    integer(ik4) :: beg1d , end1d   ! beginning and ending indices
    ! true => check 'active' flag of each point
    ! (this refers to a point being active, NOT a history field being active)
    logical :: check_active
    logical :: valid        ! true => history operation is valid
    logical :: map2gcell    ! true => map clm pointer field to gridcell
    ! 1d clm pointerr type   ["gridcell","landunit","column","pft"]
    character(len=8) :: type1d
    ! 1d history buffer type ["gridcell","landunit","column","pft"]
    character(len=8) :: type1d_out
    ! time averaging flag
    character(len=1) :: avgflag
    ! scale type for subgrid averaging of pfts to column
    character(len=8) :: p2c_scale_type
    ! scale type for subgrid averaging of columns to landunits
    character(len=8) :: c2l_scale_type
    ! scale type for subgrid averaging of landunits to gridcells
    character(len=8) :: l2g_scale_type
    real(rk8) , pointer :: hbuf(:,:)      ! history buffer
    integer(ik4) , pointer :: nacs(:,:)   ! accumulation counter
    real(rk8) , pointer :: field(:,:)     ! clm 2d pointer field
    ! flag saying whether each point is active
    ! (used for type1d = landunit/column/pft)
    ! (this refers to a point being active, NOT a history field being active)
    logical , pointer :: active(:)
    ! gricell level field (used if mapping to gridcell is done)
    real(rk8) :: field_gcell(begg:endg,num2d)
    character(len=*) , parameter :: subname = 'hist_update_hbuf_field_2d'

    avgflag        =  tape(t)%hlist(f)%avgflag
    nacs           => tape(t)%hlist(f)%nacs
    hbuf           => tape(t)%hlist(f)%hbuf
    beg1d          =  tape(t)%hlist(f)%field%beg1d
    end1d          =  tape(t)%hlist(f)%field%end1d
    type1d         =  tape(t)%hlist(f)%field%type1d
    type1d_out     =  tape(t)%hlist(f)%field%type1d_out
    p2c_scale_type =  tape(t)%hlist(f)%field%p2c_scale_type
    c2l_scale_type =  tape(t)%hlist(f)%field%c2l_scale_type
    l2g_scale_type =  tape(t)%hlist(f)%field%l2g_scale_type
    hpindex        =  tape(t)%hlist(f)%field%hpindex
    field          => clmptr_ra(hpindex)%ptr(:,1:num2d)

    ! set variables to check weights when allocate all pfts

    map2gcell = .false.
    if ( type1d_out == nameg .or. type1d_out == grlnd ) then
      if ( type1d == namep ) then
        call p2g(begp,endp,begc,endc,begl,endl,begg,endg,num2d,   &
                 field,field_gcell,p2c_scale_type,c2l_scale_type, &
                 l2g_scale_type)
        map2gcell = .true.
      else if ( type1d == namec ) then
        call c2g(begc,endc,begl,endl,begg,endg,num2d,field,field_gcell, &
                 c2l_scale_type,l2g_scale_type)
        map2gcell = .true.
      else if ( type1d == namel ) then
        call l2g(begl,endl,begg,endg,num2d,field,field_gcell, &
                 l2g_scale_type)
        map2gcell = .true.
      end if
    end if

    if ( map2gcell ) then  ! Map to gridcell

      ! note that in this case beg1d = begg and end1d=endg
      select case (avgflag)
        case ('I') ! Instantaneous
          do j = 1 , num2d
            do k = begg , endg
              if ( field_gcell(k,j) /= spval ) then
                hbuf(k,j) = field_gcell(k,j)
              else
                hbuf(k,j) = spval
              end if
              nacs(k,j) = 1
            end do
          end do
        case ('A') ! Time average
          do j = 1 , num2d
            do k = begg,endg
              if ( field_gcell(k,j) /= spval ) then
                if ( nacs(k,j) == 0 ) hbuf(k,j) = 0.D0
                hbuf(k,j) = hbuf(k,j) + field_gcell(k,j)
                nacs(k,j) = nacs(k,j) + 1
              else
                if ( nacs(k,j) == 0 ) hbuf(k,j) = spval
              end if
            end do
          end do
        case ('X') ! Maximum over time
          do j = 1 , num2d
            do k = begg , endg
              if ( field_gcell(k,j) /= spval ) then
                if ( nacs(k,j) == 0 ) hbuf(k,j) = -1.D50
                hbuf(k,j) = max( hbuf(k,j), field_gcell(k,j) )
              else
                hbuf(k,j) = spval
              end if
              nacs(k,j) = 1
            end do
          end do
        case ('M') ! Minimum over time
          do j = 1 , num2d
            do k = begg , endg
              if ( field_gcell(k,j) /= spval ) then
                if ( nacs(k,j) == 0 ) hbuf(k,j) = +1.D50
                hbuf(k,j) = min( hbuf(k,j), field_gcell(k,j) )
              else
                hbuf(k,j) = spval
              end if
              nacs(k,j) = 1
            end do
          end do
        case default
          write(stderr,*) trim(subname), &
                  ' ERROR: invalid time averaging flag ', avgflag
          call fatal(__FILE__,__LINE__,'clm now stopping.')
      end select
    else  ! Do not map to gridcell

      ! For data defined on the pft, col or landunit, we need to
      ! check if a point is active to determine whether that point
      ! should be assigned spval
      if ( type1d == namep ) then
        check_active = .true.
        active => clm3%g%l%c%p%active
      else if ( type1d == namec ) then
        check_active = .true.
        active => clm3%g%l%c%active
      else if ( type1d == namel ) then
        check_active = .true.
        active => clm3%g%l%active
      else
        check_active = .false.
      end if

      ! Note that since field points to an array section the
      ! bounds are field(1:end1d-beg1d+1, num2d) - therefore
      ! need to do the shifting below

      select case (avgflag)
        case ('I') ! Instantaneous
          do j = 1 , num2d
            do k = beg1d , end1d
              valid = .true.
              if ( check_active ) then
                if ( .not. active(k) ) valid = .false.
              end if
              if ( valid ) then
                if ( field(k-beg1d+1,j) /= spval ) then
                  hbuf(k,j) = field(k-beg1d+1,j)
                else
                  hbuf(k,j) = spval
                end if
              else
                hbuf(k,j) = spval
              end if
              nacs(k,j) = 1
            end do
          end do
        case ('A') ! Time average
          do j = 1 , num2d
            do k = beg1d , end1d
              valid = .true.
              if ( check_active ) then
                if ( .not. active(k) ) valid = .false.
              end if
              if ( valid ) then
                if ( field(k-beg1d+1,j) /= spval ) then
                  if ( nacs(k,j) == 0 ) hbuf(k,j) = 0.D0
                  hbuf(k,j) = hbuf(k,j) + field(k-beg1d+1,j)
                  nacs(k,j) = nacs(k,j) + 1
                else
                  if ( nacs(k,j) == 0 ) hbuf(k,j) = spval
                end if
              else
                if ( nacs(k,j) == 0 ) hbuf(k,j) = spval
              end if
            end do
          end do
        case ('X') ! Maximum over time
          do j = 1 , num2d
            do k = beg1d , end1d
              valid = .true.
              if ( check_active ) then
                if ( .not. active(k) ) valid = .false.
              end if
              if ( valid ) then
                if ( field(k-beg1d+1,j) /= spval ) then
                  if ( nacs(k,j) == 0 ) hbuf(k,j) = -1.D50
                  hbuf(k,j) = max( hbuf(k,j), field(k-beg1d+1,j) )
                else
                  if ( nacs(k,j) == 0 ) hbuf(k,j) = spval
                end if
              else
                if ( nacs(k,j) == 0 ) hbuf(k,j) = spval
              end if
              nacs(k,j) = 1
            end do
          end do
        case ('M') ! Minimum over time
          do j = 1 , num2d
            do k = beg1d , end1d
              valid = .true.
              if ( check_active ) then
                if ( .not. active(k) ) valid = .false.
              end if
              if ( valid ) then
                if ( field(k-beg1d+1,j) /= spval ) then
                  if ( nacs(k,j) == 0 ) hbuf(k,j) = +1.D50
                  hbuf(k,j) = min( hbuf(k,j), field(k-beg1d+1,j))
                else
                  if ( nacs(k,j) == 0 ) hbuf(k,j) = spval
                end if
              else
                if ( nacs(k,j) == 0 ) hbuf(k,j) = spval
              end if
              nacs(k,j) = 1
            end do
          end do
        case default
          write(stderr,*) trim(subname), &
                  ' ERROR: invalid time averaging flag ', avgflag
          call fatal(__FILE__,__LINE__,'clm now stopping.')
      end select
    end if
  end subroutine hist_update_hbuf_field_2d
  !
  ! Normalize fields on a history file by the number of accumulations.
  ! Loop over fields on the tape.  Need averaging flag and number of
  ! accumulations to perform normalization.
  !
  subroutine hfields_normalize (t)
    implicit none
    integer(ik4) , intent(in) :: t ! tape index
    integer(ik4) :: f              ! field index
    integer(ik4) :: k              ! 1d index
    integer(ik4) :: j              ! 2d index
    logical :: aflag               ! averaging flag
    ! hbuf 1d beginning and ending indices
    integer(ik4) :: beg1d_out , end1d_out
    ! hbuf size of second dimension (e.g. number of vertical levels)
    integer(ik4) :: num2d
    character(len=1) :: avgflag         ! averaging flag
    real(rk8) , pointer :: hbuf(:,:)    ! history buffer
    integer(ik4) , pointer :: nacs(:,:) ! accumulation counter
    character(len=*) , parameter :: subname = 'hfields_normalize'

    ! Normalize by number of accumulations for time averaged case

    do f = 1 , tape(t)%nflds
      avgflag   =  tape(t)%hlist(f)%avgflag
      beg1d_out =  tape(t)%hlist(f)%field%beg1d_out
      end1d_out =  tape(t)%hlist(f)%field%end1d_out
      num2d     =  tape(t)%hlist(f)%field%num2d
      nacs      => tape(t)%hlist(f)%nacs
      hbuf      => tape(t)%hlist(f)%hbuf

      if (avgflag == 'A') then
        aflag = .true.
      else
        aflag = .false.
      end if

      do j = 1 , num2d
        do k = beg1d_out , end1d_out
          if ( aflag .and. nacs(k,j) /= 0 ) then
            hbuf(k,j) = hbuf(k,j) / float(nacs(k,j))
          end if
        end do
      end do
    end do
  end subroutine hfields_normalize
  !
  ! Zero out accumulation and history buffers for a given history tape.
  ! Loop through fields on the tape.
  !
  subroutine hfields_zero (t)
    implicit none
    integer(ik4) , intent(in) :: t    ! tape index
    integer(ik4) :: f                 ! field index
    character(len=*) , parameter :: subname = 'hfields_zero'
    do f = 1,tape(t)%nflds
      tape(t)%hlist(f)%hbuf(:,:) = 0.D0
      tape(t)%hlist(f)%nacs(:,:) = 0
    end do
  end subroutine hfields_zero
  !
  ! Define contents of history file t. Issue the required netcdf
  ! wrapper calls to define the history file contents.
  !
  subroutine htape_create (t, histrest)
    implicit none
    integer(ik4) , intent(in) :: t    ! tape index
    ! if creating the history restart file
    logical , intent(in) , optional :: histrest
    integer(ik4) :: f                   ! field index
    integer(ik4) :: p , c , l , n       ! indices
    integer(ik4) :: ier                 ! error code
    ! size of second dimension (e.g. number of vertical levels)
    integer(ik4) :: num2d
    integer(ik4) :: dimid      ! dimension id temporary
    integer(ik4) :: dim1id(1)  ! netCDF dimension id
    integer(ik4) :: dim2id(2)  ! netCDF dimension id
    integer(ik4) :: ndims      ! dimension counter
    integer(ik4) :: omode      ! returned mode from netCDF call
    integer(ik4) :: ncprec     ! output netCDF write precision
    integer(ik4) :: ret        ! netCDF error status
    integer(ik4) :: nump       ! total number of pfts across all processors
    integer(ik4) :: numc       ! total number of columns across all processors
    integer(ik4) :: numl       ! total number of landunits across all processors
    integer(ik4) :: numg       ! total number of gridcells across all processors
    integer(ik4) :: numa       ! total number of atm cells across all processors
    logical :: lhistrest       ! local history restart flag
    type(clm_filetype) :: lnfid    ! local file id
    character(len=  8) :: curdate  ! current date
    character(len=  8) :: curtime  ! current time
    character(len=256) :: name     ! name of attribute
    character(len=256) :: units    ! units of attribute
    character(len=256) :: str      ! global attribute string
    character(len=  1) :: avgflag  ! time averaging flag
    character(len=*) , parameter :: subname = 'htape_create'

    if ( present(histrest) )then
      lhistrest = histrest
    else
      lhistrest = .false.
    end if

    ! Determine necessary indices

    call get_proc_global(numg,numl,numc,nump)

    ! define output write precsion for tape

    ncprec = tape(t)%ncprec

    ! Create new netCDF file. It will be in define mode

    if ( .not. lhistrest ) then
      if (myid == italk) then
        write(stdout,*) trim(subname), &
                ' : Opening netcdf htape ', trim(locfnh(t))
      end if
      call clm_createfile(trim(locfnh(t)),lnfid)
      call clm_addatt(lnfid,'title','CLM History file information')
      call clm_addatt(lnfid,'comment', &
          "NOTE: None of the variables are weighted by land fraction!" )
    else
      if ( myid == italk ) then
        write(stdout,*) trim(subname), &
                ' : Opening netcdf rhtape ',trim(locfnhr(t))
      end if
      call clm_createfile(trim(locfnhr(t)),lnfid)
      call clm_addatt(lnfid,'title', &
          'CLM Restart History information, required to continue a simulation' )
      call clm_addatt(lnfid,'comment', &
          'This entire file NOT needed for startup simulations')
    end if

    ! Create global attributes. Attributes are used to store information
    ! about the data set. Global attributes are information about the
    ! data set as a whole, as opposed to a single variable

    call clm_addatt(lnfid,'Conventions', trim(conventions))
    call getdatetime(curdate,curtime)
    str = 'created on ' // curdate // ' ' // curtime
    call clm_addatt(lnfid, 'history' , trim(str))
    call clm_addatt(lnfid, 'source'  , trim(source))
    call clm_addatt(lnfid, 'hostname', trim(hostname))
    call clm_addatt(lnfid, 'username', trim(username))
    call clm_addatt(lnfid, 'version' , trim(version))

    str = SVN_REV
    call clm_addatt(lnfid, 'revision_id', trim(str))
    call clm_addatt(lnfid, 'case_title', trim(ctitle))
    call clm_addatt(lnfid, 'case_id', trim(caseid))
    call clm_addatt(lnfid, 'Surface_dataset', trim(fsurdat))
    if (finidat == ' ') then
      str = 'arbitrary initialization'
    else
      str = finidat
    end if
    call clm_addatt(lnfid, 'Initial_conditions_dataset', trim(str))
    call clm_addatt(lnfid, 'PFT_physiological_constants_dataset', trim(fpftcon))

    ! Define dimensions.
    ! Time is an unlimited dimension. Character string is treated as
    ! an array of characters.

    ! Global uncompressed dimensions (including non-land points)
    call clm_adddim(lnfid, trim(grlnd), ldomain%ns)

    ! Global compressed dimensions (not including non-land points)
    call clm_adddim(lnfid, trim(nameg), numg)
    call clm_adddim(lnfid, trim(namel), numl)
    call clm_adddim(lnfid, trim(namec), numc)
    call clm_adddim(lnfid, trim(namep), nump)

    ! "level" dimensions
    call clm_adddim(lnfid, 'levgrnd', nlevgrnd)
    if ( nlevurb > 0 ) then
       call clm_adddim(lnfid, 'levurb', nlevurb)
    end if
    call clm_adddim(lnfid, 'levlak', nlevlak)
    call clm_adddim(lnfid, 'numrad', numrad)

    do n = 1 , num_subs
      call clm_adddim(lnfid, subs_name(n), subs_dim(n))
    end do
    call clm_adddim(lnfid, 'string_length', 8)
    call clm_adddim(lnfid, 'levdcmp', nlevdecomp_full)
    if ( .not. lhistrest ) then
      call clm_adddim(lnfid, 'hist_interval', 2)
      call clm_adddim(lnfid, 'time', clmvar_unlim)
      nfid(t) = lnfid
      if ( myid == italk ) then
        write(stdout,*) trim(subname), &
           ' : Successfully defined netcdf history file ',t
      end if
    else
      ncid_hist(t) = lnfid
      if ( myid == italk ) then
        write(stdout,*) trim(subname), &
           ' : Successfully defined netcdf restart history file ',t
      end if
    end if
  end subroutine htape_create
  !
  ! Write time constant 3D variables to history tapes.
  ! Only write out when this subroutine is called (normally only for
  ! primary history files at very first time-step, nstep=0).
  ! Issue the required netcdf wrapper calls to define the history file
  ! contents.
  !
  subroutine htape_timeconst3D(t, mode)
    implicit none
    integer(ik4) , intent(in) :: t  ! tape index
    character(len=*) , intent(in) :: mode  ! 'define' or 'write'
    integer(ik4) :: c , l , lev , ifld ! indices
    integer(ik4) :: ier                ! error status
    integer(ik4) :: begp , endp ! per-proc beginning and ending pft indices
    integer(ik4) :: begc , endc ! per-proc beginning and ending column indices
    integer(ik4) :: begl , endl ! per-proc beginning and ending landunit indices
    integer(ik4) :: begg , endg ! per-proc gridcell ending gridcell indices
    character(len=max_chars) :: long_name ! variable long name
    character(len=max_namlen) :: varname  ! variable name
    character(len=max_namlen) :: units    ! variable units
    ! scale type for subgrid averaging of landunits to grid cells
    character(len=8) :: l2g_scale_type
    real(rk8) , pointer :: histi(:,:)     ! temporary
    real(rk8) , pointer :: histo(:,:)     ! temporary
    type(landunit_type) , pointer :: lptr ! pointer to landunit derived subtype
    type(column_type) , pointer :: cptr   ! pointer to column derived subtype
    integer(ik4) , parameter :: nflds = 6 ! Number of 3D time-constant fields
    character(len=*) , parameter :: subname = 'htape_timeconst3D'
    character(len=*) , parameter :: varnames(nflds) = (/ &
                                                        'ZSOI  ', &
                                                        'DZSOI ', &
                                                        'WATSAT', &
                                                        'SUCSAT', &
                                                        'BSW   ', &
                                                        'HKSAT '  &
                                                    /)
    real(rk8) , pointer :: histil(:,:)      ! temporary
    real(rk8) , pointer :: histol(:,:)
    integer(ik4) , parameter :: nfldsl = 2
    character(len=*) , parameter :: varnamesl(nfldsl) = (/ &
                                                          'ZLAKE ', &
                                                          'DZLAKE' &
                                                        /)
    !
    ! Non-time varying 3D fields
    ! Only write out when this subroutine is called
    ! Normally only called once for primary tapes
    !
    if ( mode == 'define' ) then
      do ifld = 1 , nflds
        ! Field indices MUST match varnames array order above!
        if ( ifld == 1 ) then
          long_name='soil depth'
          units = 'm'
        else if (ifld == 2) then
          long_name='soil thickness'
          units = 'm'
        else if (ifld == 3) then
          long_name='saturated soil water content (porosity)'
          units = 'mm3/mm3'
        else if (ifld == 4) then
          long_name='saturated soil matric potential'
          units = 'mm'
        else if (ifld == 5) then
          long_name='slope of soil water retention curve'
          units = 'unitless'
        else if (ifld == 6) then
          long_name='saturated hydraulic conductivity'
          units = 'unitless'
        else
          call fatal(__FILE__,__LINE__, &
               trim(subname)//' ERROR: bad 3D time-constant field index' )
        end if
        if ( tape(t)%dov2xy ) then
          call clm_addvar(clmvar_double,ncid=nfid(t), &
                   varname=trim(varnames(ifld)),      &
                   cdims=(/grlnd,'levgrnd'/),         &
                   long_name=long_name, units=units,  &
                   missing_value=1, fill_value=1)
        else
          call clm_addvar(clmvar_double,ncid=nfid(t), &
                  varname=trim(varnames(ifld)),       &
                  cdims=(/namec,'levgrnd'/),          &
                  long_name=long_name, units=units,   &
                  missing_value=1, fill_value=1)
        end if
      end do
    else if (mode == 'write') then

      ! Set pointers into derived type and get necessary bounds

      lptr => clm3%g%l
      cptr => clm3%g%l%c

      call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

      allocate(histi(begc:endc,nlevgrnd), stat=ier)
      if ( ier /= 0 ) then
        call fatal(__FILE__,__LINE__, &
                trim(subname)//' ERROR: allocation error for histi')
      end if

      ! Write time constant fields

      if ( tape(t)%dov2xy ) then
        allocate(histo(begg:endg,nlevgrnd), stat=ier)
        if ( ier /= 0 ) then
          call fatal(__FILE__,__LINE__, &
             trim(subname)//' ERROR: allocation error for histo')
        end if
      end if

      do ifld = 1 , nflds
        ! WJS (10-25-11): Note about l2g_scale_type in the following:
        ! ZSOI & DZSOI are currently constant in space, except for
        ! urban points, so their scale type doesn't matter at the
        ! moment as long as it excludes urban points. I am using
        ! 'nonurb' so that the values are output everywhere where
        ! the fields are constant (i.e., everywhere except urban points).
        ! For the other fields, I am using 'veg' to be consistent with
        ! the l2g_scale_type that is now used for many of the 3-d
        ! time-variant fields; in theory, though, one might want versions of
        ! these variables output for different landunits.

        ! Field indices MUST match varnames array order above!
        if ( ifld == 1 ) then  ! ZSOI
          l2g_scale_type = 'nonurb'
        else if ( ifld == 2 ) then  ! DZSOI
          l2g_scale_type = 'nonurb'
        else if ( ifld == 3 ) then  ! WATSAT
          l2g_scale_type = 'veg'
        else if ( ifld == 4 ) then  ! SUCSAT
          l2g_scale_type = 'veg'
        else if ( ifld == 5 ) then  ! BSW
          l2g_scale_type = 'veg'
        else if ( ifld == 6 ) then  ! HKSAT
          l2g_scale_type = 'veg'
        end if

        histi(:,:) = spval
        do lev = 1 , nlevgrnd
          do c = begc , endc
            l = cptr%landunit(c)
            ! Field indices MUST match varnames array order above!
            if ( ifld == 1 ) histi(c,lev) = cptr%cps%z(c,lev)
            if ( ifld == 2 ) histi(c,lev) = cptr%cps%dz(c,lev)
            if ( ifld == 3 ) histi(c,lev) = cptr%cps%watsat(c,lev)
            if ( ifld == 4 ) histi(c,lev) = cptr%cps%sucsat(c,lev)
            if ( ifld == 5 ) histi(c,lev) = cptr%cps%bsw(c,lev)
            if ( ifld == 6 ) histi(c,lev) = cptr%cps%hksat(c,lev)
          end do
        end do
        if ( tape(t)%dov2xy ) then
          histo(:,:) = spval
          call c2g(begc,endc,begl,endl,begg,endg,nlevgrnd,histi,histo, &
                   c2l_scale_type='unity',l2g_scale_type=l2g_scale_type)
          call clm_writevar(nfid(t),trim(varnames(ifld)),histo,gcomm_gridcell)
        else
          call clm_writevar(nfid(t),trim(varnames(ifld)),histi,gcomm_gridcell)
        end if
      end do

      if ( tape(t)%dov2xy ) deallocate(histo)
      deallocate(histi)

    end if  ! (define/write mode

    if ( mode == 'define' ) then
      do ifld = 1 , nfldsl
        ! Field indices MUST match varnamesl array order above!
        if ( ifld == 1 ) then
          long_name='lake layer node depth'
          units = 'm'
        else if ( ifld == 2 ) then
          long_name='lake layer thickness'
          units = 'm'
        else
          call fatal(__FILE__,__LINE__, &
              trim(subname)//' ERROR: bad 3D time-constant field index' )
        end if
        if ( tape(t)%dov2xy ) then
          call clm_addvar(clmvar_double,ncid=nfid(t), &
                   varname=trim(varnamesl(ifld)),     &
                   cdims=(/grlnd,'levlak'/),          &
                   long_name=long_name, units=units,  &
                   missing_value=1, fill_value=1)
        else
          call clm_addvar(clmvar_double,ncid=nfid(t), &
                   varname=trim(varnamesl(ifld)),     &
                   cdims=(/namec,'levlak'/),          &
                   long_name=long_name, units=units,  &
                   missing_value=1, fill_value=1)
        end if
      end do
    else if ( mode == 'write' ) then

      ! Set pointers into derived type and get necessary bounds

      lptr => clm3%g%l
      cptr => clm3%g%l%c

      call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

      allocate(histil(begc:endc,nlevlak), stat=ier)
      if ( ier /= 0 ) then
        call fatal(__FILE__,__LINE__, &
           trim(subname)//' ERROR: allocation error for histil')
      end if

      ! Write time constant fields

      if (tape(t)%dov2xy) then
        allocate(histol(begg:endg,nlevlak), stat=ier)
        if ( ier /= 0 ) then
          call fatal(__FILE__,__LINE__, &
              trim(subname)//' ERROR: allocation error for histol')
        end if
      end if

      do ifld = 1 , nfldsl
        histil(:,:) = spval
        do lev = 1 , nlevlak
          do c = begc , endc
            l = cptr%landunit(c)
            if ( lptr%lakpoi(l) ) then
              ! Field indices MUST match varnamesl array order above!
              if ( ifld == 1 ) histil(c,lev) = cptr%cps%z_lake(c,lev)
              if ( ifld == 2 ) histil(c,lev) = cptr%cps%dz_lake(c,lev)
            end if
          end do
        end do
        if ( tape(t)%dov2xy ) then
          histol(:,:) = spval
          call c2g(begc,endc,begl,endl,begg,endg,nlevlak,histil,histol, &
                   c2l_scale_type='unity', l2g_scale_type='lake')
          call clm_writevar(nfid(t),trim(varnamesl(ifld)),histol, &
             gcomm_gridcell)
        else
          call clm_writevar(nfid(t),trim(varnamesl(ifld)),histil,gcomm_gridcell)
        end if
      end do

      if ( tape(t)%dov2xy ) deallocate(histol)
      deallocate(histil)
    end if  ! (define/write mode
  end subroutine htape_timeconst3D
  !
  ! Write time constant values to primary history tape.
  ! Issue the required netcdf wrapper calls to define the history file
  ! contents.
  !
  subroutine htape_timeconst(t, mode)
    implicit none
    integer(ik4) , intent(in) :: t          ! tape index
    character(len=*) , intent(in) :: mode   ! 'define' or 'write'
    integer(ik4) :: vid , n , i , j , m     ! indices
    integer(ik4) :: mcsec                   ! seconds of current date
    integer(ik4) :: mdcur                   ! current day
    integer(ik4) :: mscur                   ! seconds of current day
    integer(ik4) :: mcdate                  ! current date
    integer(ik4) :: yr , mon , day , nbsec  ! year,month,day,seconds
    integer(ik4) :: hours , minutes , secs  ! hours,minutes,seconds of hh:mm:ss
    character(len= 10) :: basedate          ! base date (yyyymmdd)
    character(len=  8) :: basesec           ! base seconds
    character(len=  8) :: cdate             ! system date
    character(len=  8) :: ctime             ! system time
    real(rk8) :: time                       ! current time
    real(rk8) :: timedata(2)                ! time interval boundaries
    integer(ik4) :: dim1id(1)               ! netCDF dimension id
    integer(ik4) :: dim2id(2)               ! netCDF dimension id
    integer(ik4) :: varid                   ! netCDF variable id
    integer(ik4) :: begp , endp ! per-proc beginning and ending pft indices
    integer(ik4) :: begc , endc ! per-proc beginning and ending column indices
    integer(ik4) :: begl , endl ! per-proc beginning and ending landunit indices
    integer(ik4) :: begg , endg ! per-proc gridcell ending gridcell indices
    character(len=max_chars) :: long_name ! variable long name
    character(len=max_namlen) :: varname  ! variable name
    character(len=max_namlen) :: units    ! variable units
    character(len=max_namlen) :: cal      ! calendar from the time-manager
    character(len=max_namlen) :: caldesc  ! calendar description to put on file
    character(len=256) :: str             ! global attribute string
    real(rk8) , pointer :: histo(:,:)     ! temporary
    type(landunit_type) , pointer :: lptr ! pointer to landunit derived subtype
    type(column_type) , pointer :: cptr   ! pointer to column derived subtype
    integer(ik4) :: status
    real(rk8) :: zsoi_1d(1)

    character(len=*) , parameter :: subname = 'htape_timeconst'

    ! Time constant grid variables only on first time-sample of file

    if ( tape(t)%ntimes == 1 ) then
      if (mode == 'define') then
        call clm_addvar(clmvar_double,ncid=nfid(t), &
                     varname='levgrnd',cdims=(/'levgrnd'/) , &
                     long_name='coordinate soil levels', units='m')
        call clm_addvar(clmvar_double,ncid=nfid(t), &
                     varname='levlak',cdims=(/'levlak'/) , &
                     long_name='coordinate lake levels', units='m')
        call clm_addvar(clmvar_double,ncid=nfid(t), &
                     varname='levdcmp',cdims=(/'levlak'/) , &
                     long_name='coordinate soil levels', units='m')
      else if ( mode == 'write' ) then
        call clm_writevar(nfid(t),'levgrnd',zsoi)
        call clm_writevar(nfid(t),'levlak',zlak)
#ifdef VERTSOILC
        call clm_writevar(nfid(t),'levdcmp',zsoi)
#else
        zsoi_1d(1) = 1.D0
        call clm_writevar(nfid(t),'levdcmp',zsoi_1d)
#endif
      end if
    end if

    ! Time definition variables

    ! For define mode -- only do this for first time-sample
    if (mode == 'define' .and. tape(t)%ntimes == 1) then
      call ref_date(yr, mon, day, nbsec)
      hours   = nbsec / 3600
      minutes = (nbsec - hours*3600) / 60
      secs    = (nbsec - hours*3600 - minutes*60)
      write(basedate,80) yr,mon,day
80    format(i4.4,'-',i2.2,'-',i2.2)
      write(basesec ,90) hours, minutes, secs
90    format(i2.2,':',i2.2,':',i2.2)
      str = 'hours since ' // basedate // " " // basesec
      call clm_addvar(clmvar_double,ncid=nfid(t), &
                      varname='time',cdims=(/'time'/), &
                      long_name='time', units=str)
      call clm_addatt(nfid(t),'calendar',trim(calendar_str(idatex)),cvar='time')
      call clm_addatt(nfid(t),'bounds','time_bounds',cvar='time')
      call clm_addvar(clmvar_double,ncid=nfid(t), &
                      varname='time_bounds',      &
                      cdims=(/'hist_interval','time         '/), &
                      long_name='history time interval endpoints')
    else if (mode == 'write') then
      call curr_time(idatex,mdcur,mscur)
      time = mdcur*24.0D0 + mscur/secph
      call clm_writevar(nfid(t),'time',time,tape(t)%ntimes)
      timedata(1) = tape(t)%begtime
      timedata(2) = time
      call clm_writevar(nfid(t),'time_bounds',timedata,tape(t)%ntimes)
    end if

    ! Grid definition variables
    ! For define mode -- only do this for first time-sample

    if ( mode == 'define' .and. tape(t)%ntimes == 1 ) then
      call clm_addvar(clmvar_double,ncid=nfid(t), &
                    varname='lon',cdims=(/grlnd/), &
                    long_name='coordinate longitude', units='degrees_east')
      call clm_addvar(clmvar_double,ncid=nfid(t), &
                    varname='lat',cdims=(/grlnd/), &
                    long_name='coordinate latitude', units='degrees_north')
      call clm_addvar(clmvar_double,ncid=nfid(t), &
                    varname='area',cdims=(/grlnd/), &
                    long_name='grid cell areas', units='km^2')
      call clm_addvar(clmvar_double,ncid=nfid(t), &
                    varname='topo',cdims=(/grlnd/), &
                    long_name='grid cell topography', units='m')
      call clm_addvar(clmvar_double,ncid=nfid(t), &
                    varname='landfrac',cdims=(/grlnd/), &
                    long_name='land fraction', units='1')
      call clm_addvar(clmvar_double,ncid=nfid(t), &
                    varname='landmask',cdims=(/grlnd/), &
                    long_name='land/ocean mask', units='1')
      call clm_addvar(clmvar_double,ncid=nfid(t), &
                    varname='pftmask',cdims=(/grlnd/), &
                    long_name='pft real/fake mask', units='1')

    else if ( mode == 'write' ) then

      ! Most of this is constant and only needs to be done on tape(t)%ntimes=1
      ! But, some may change for dynamic PFT mode for example
      ! Set pointers into derived type and get necessary bounds

      lptr => clm3%g%l
      cptr => clm3%g%l%c

      call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

      call clm_writevar(nfid(t),'lon',ldomain%lonc,gcomm_gridcell)
      call clm_writevar(nfid(t),'lat',ldomain%latc,gcomm_gridcell)
      call clm_writevar(nfid(t),'area',ldomain%area,gcomm_gridcell)
      call clm_writevar(nfid(t),'landfrac',ldomain%frac,gcomm_gridcell)
      call clm_writevar(nfid(t),'landmask',ldomain%mask,gcomm_gridcell)
      call clm_writevar(nfid(t),'pftmask',ldomain%pftm,gcomm_gridcell)

    end if  ! (define/write mode

  end subroutine htape_timeconst
  !
  ! Write history tape.  Issue the call to write the variable.
  !
  subroutine hfields_write(t, mode)
    implicit none
    integer(ik4) , intent(in) :: t        ! tape index
    character(len=*) , intent(in) :: mode ! 'define' or 'write'
    integer(ik4) :: f         ! field index
    integer(ik4) :: k         ! 1d index
    integer(ik4) :: c , l , p ! indices
    integer(ik4) :: beg1d_out ! on-node 1d hbuf pointer start index
    integer(ik4) :: end1d_out ! on-node 1d hbuf pointer end index
    integer(ik4) :: num1d_out ! size of hbuf first dimension (overall all nodes)
    integer(ik4) :: num2d     ! hbuf second dimension size
    integer(ik4) :: nt        ! time index
    integer(ik4) :: ier       ! error status
    character(len=1) :: avgflag          ! time averaging flag
    character(len=max_chars) :: long_name! long name
    character(len=max_chars) :: units    ! units
    character(len=max_namlen):: varname  ! variable name
    character(len=32) :: avgstr          ! time averaging type
    character(len=8) :: type1d_out       ! history output 1d type
    character(len=8) :: type2d           ! history output 2d type
    character(len=32) :: dim1name        ! temporary
    character(len=32) :: dim2name        ! temporary
    real(rk8) , pointer :: histo(:,:)    ! temporary
    real(rk8) , pointer :: hist1do(:)    ! temporary
    character(len=*) , parameter :: subname = 'hfields_write'

    ! Write/define 1d topological info

    if ( .not. tape(t)%dov2xy ) then
      if ( mode == 'define' ) then
        call hfields_1dinfo(t, mode='define')
      else if ( mode == 'write' ) then
        call hfields_1dinfo(t, mode='write')
      end if
    end if

    ! Define time-dependent variables create variables and
    ! attributes for field list

    do f = 1 , tape(t)%nflds

      ! Set history field variables

      varname    = tape(t)%hlist(f)%field%name
      long_name  = tape(t)%hlist(f)%field%long_name
      units      = tape(t)%hlist(f)%field%units
      avgflag    = tape(t)%hlist(f)%avgflag
      type1d_out = tape(t)%hlist(f)%field%type1d_out
      beg1d_out  = tape(t)%hlist(f)%field%beg1d_out
      end1d_out  = tape(t)%hlist(f)%field%end1d_out
      num1d_out  = tape(t)%hlist(f)%field%num1d_out
      type2d     = tape(t)%hlist(f)%field%type2d
      num2d      = tape(t)%hlist(f)%field%num2d
      nt         = tape(t)%ntimes

      if ( mode == 'define' ) then

        select case (avgflag)
          case ('A')
            avgstr = 'mean'
          case ('I')
            avgstr = 'instantaneous'
          case ('X')
            avgstr = 'maximum'
          case ('M')
            avgstr = 'minimum'
          case default
            write(stderr,*) trim(subname), &
                  ' ERROR: unknown time averaging flag (avgflag)=',avgflag
            call fatal(__FILE__,__LINE__,'clm now stopping.')
        end select

        if ( type1d_out == grlnd ) then
          dim1name = trim(grlnd)
          dim2name = 'undefined'
        else
          dim1name = type1d_out
          dim2name = 'undefined'
        end if

        if ( dim2name == 'undefined' ) then
          if ( num2d == 1 ) then
            call clm_addvar(tape(t)%ncprec,nfid(t),varname, &
              cdims=(/dim1name,'time'/), long_name=long_name, &
              units=units, cell_method=avgstr, missing_value=1, fill_value=1)
          else
            call clm_addvar(tape(t)%ncprec,nfid(t),varname, &
              cdims=(/dim1name,type2d,'time'/), long_name=long_name, &
              units=units, cell_method=avgstr, missing_value=1, fill_value=1)
          end if
        else
          if ( num2d == 1 ) then
            call clm_addvar(tape(t)%ncprec,nfid(t),varname, &
              cdims=(/dim1name,dim2name,'time'/), long_name=long_name, &
              units=units, cell_method=avgstr, missing_value=1, fill_value=1)
          else
            call clm_addvar(tape(t)%ncprec,nfid(t),varname, &
              cdims=(/dim1name,dim2name,type2d,'time'/), long_name=long_name, &
              units=units, cell_method=avgstr, missing_value=1, fill_value=1)
          end if
        end if

      else if ( mode == 'write' ) then

        ! Determine output buffer

        histo => tape(t)%hlist(f)%hbuf

        ! Allocate dynamic memory

        if ( num2d == 1 ) then
          allocate(hist1do(beg1d_out:end1d_out), stat=ier)
          if ( ier /= 0 ) then
            write(stderr,*) trim(subname),' ERROR: allocation'
            call fatal(__FILE__,__LINE__,'clm now stopping.')
          end if
          hist1do(beg1d_out:end1d_out) = histo(beg1d_out:end1d_out,1)
        end if

        ! Write history output.  Always output land and ocean runoff on xy grid.

        if ( num2d == 1 ) then
          call clm_writevar(nfid(t),varname,hist1do,procinfo,nt)
        else
          call clm_writevar(nfid(t),varname,histo,procinfo,nt)
        end if

        ! Deallocate dynamic memory

        if ( num2d == 1 ) then
          deallocate(hist1do)
        end if
      end if
    end do
  end subroutine hfields_write
  !
  ! Write/define 1d info for history tape.
  !
  subroutine hfields_1dinfo(t, mode)
    implicit none
    integer(ik4) , intent(in) :: t         ! tape index
    character(len=*) , intent(in) :: mode  ! 'define' or 'write'
    integer(ik4) :: f                ! field index
    integer(ik4) :: k                ! 1d index
    integer(ik4) :: g , c , l , p    ! indices
    integer(ik4) :: begp , endp ! per-proc beginning and ending pft indices
    integer(ik4) :: begc , endc ! per-proc beginning and ending column indices
    integer(ik4) :: begl , endl ! per-proc beginning and ending landunit indices
    integer(ik4) :: begg , endg ! per-proc gridcell ending gridcell indices
    integer(ik4) :: ier         ! errir status
    real(rk8) , pointer :: rgarr(:)        ! temporary
    real(rk8) , pointer :: rcarr(:)        ! temporary
    real(rk8) , pointer :: rlarr(:)        ! temporary
    real(rk8) , pointer :: rparr(:)        ! temporary
    integer(ik4) , pointer :: igarr(:)     ! temporary
    integer(ik4) , pointer :: icarr(:)     ! temporary
    integer(ik4) , pointer :: ilarr(:)     ! temporary
    integer(ik4) , pointer :: iparr(:)     ! temporary
    type(gridcell_type), pointer :: gptr ! pointer to gridcell derived subtype
    type(landunit_type), pointer :: lptr ! pointer to landunit derived subtype
    type(column_type)  , pointer :: cptr ! pointer to column derived subtype
    type(pft_type)     , pointer :: pptr ! pointer to pft derived subtype
    character(len=*) , parameter :: subname = 'hfields_1dinfo'

    if ( mode == 'define' ) then

      ! Define gridcell info

      call clm_addvar(clmvar_double,nfid(t),'grid1d_lon',cdims=(/nameg/), &
              long_name='gridcell longitude', units='degrees_east')
      call clm_addvar(clmvar_double,nfid(t),'grid1d_lat',cdims=(/nameg/), &
              long_name='gridcell latitude', units='degrees_north')

      ! Define landunit info
      call clm_addvar(clmvar_double,nfid(t),'land1d_lon',cdims=(/namel/), &
              long_name='ladunit longitude', units='degrees_east')
      call clm_addvar(clmvar_double,nfid(t),'land1d_lat',cdims=(/namel/), &
              long_name='ladunit latitude', units='degrees_north')
      call clm_addvar(clmvar_double,nfid(t),'land1d_wtgcell',cdims=(/namel/), &
              long_name='landunit weight relative to corresponding gridcell')
      call clm_addvar(clmvar_integer,nfid(t),'land1d_ityplunit', &
              cdims=(/namel/), &
              long_name='landunit type (vegetated,urban,lake,wetland,glacier)')
      call clm_addvar(clmvar_logical,nfid(t),'land1d_active', &
              cdims=(/namel/), &
              long_name='true => do computations on this landunit')

      ! Define column info
      call clm_addvar(clmvar_double,nfid(t),'cols1d_lon',cdims=(/namec/), &
              long_name='column longitude', units='degrees_east')
      call clm_addvar(clmvar_double,nfid(t),'cols1d_lat',cdims=(/namec/), &
              long_name='column latitude', units='degrees_north')
      call clm_addvar(clmvar_double,nfid(t),'cols1d_wtgcell',cdims=(/namec/), &
              long_name='column weight relative to corresponding gridcell')
      call clm_addvar(clmvar_double,nfid(t),'cols1d_wtlunit',cdims=(/namec/), &
              long_name='column weight relative to corresponding landunit')
      call clm_addvar(clmvar_integer,nfid(t),'cols1d_itype_lunit', &
        cdims=(/namec/), &
        long_name='column landunit type (vegetated,urban,lake,wetland,glacier)')
      call clm_addvar(clmvar_logical,nfid(t),'cols1d_active', &
        cdims=(/namec/), long_name='true => do computations on this column')

      ! Define pft info
      call clm_addvar(clmvar_double,nfid(t),'pfts1d_lon',cdims=(/namep/), &
              long_name='column longitude', units='degrees_east')
      call clm_addvar(clmvar_double,nfid(t),'pfts1d_lat',cdims=(/namep/), &
              long_name='column latitude', units='degrees_north')
      call clm_addvar(clmvar_double,nfid(t),'pfts1d_wtgcell',cdims=(/namep/), &
              long_name='pft weight relative to corresponding gridcell')
      call clm_addvar(clmvar_double,nfid(t),'pfts1d_wtlunit',cdims=(/namep/), &
              long_name='pft weight relative to corresponding landunit')
      call clm_addvar(clmvar_double,nfid(t),'pfts1d_wtcol',cdims=(/namep/), &
              long_name='pft weight relative to corresponding column')
      call clm_addvar(clmvar_integer,nfid(t),'pfts1d_itype_veg', &
        cdims=(/namep/),long_name='pft vegetation type')
      call clm_addvar(clmvar_integer,nfid(t),'pfts1d_itype_lunit', &
        cdims=(/namep/), &
        long_name='pft landunit type (vegetated,urban,lake,wetland,glacier)')
      call clm_addvar(clmvar_logical,nfid(t),'pfts1d_active', &
        cdims=(/namep/), &
        long_name='true => do computations on this pft')

    else if (mode == 'write') then

      ! Set pointers into derived type

      gptr => clm3%g
      lptr => clm3%g%l
      cptr => clm3%g%l%c
      pptr => clm3%g%l%c%p

      ! Determine bounds

      call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

      allocate(rgarr(begg:endg),rlarr(begl:endl), &
               rcarr(begc:endc),rparr(begp:endp),stat=ier)
      if ( ier /= 0 ) then
        call fatal(__FILE__,__LINE__,'hfields_1dinfo allocation error of rarrs')
      end if

      allocate(igarr(begg:endg),ilarr(begl:endl), &
               icarr(begc:endc),iparr(begp:endp),stat=ier)
      if ( ier /= 0 ) then
        call fatal(__FILE__,__LINE__,'hfields_1dinfo allocation error of iarrs')
      end if

      ! Write gridcell info

      call clm_writevar(nfid(t),'grid1d_lon',gptr%londeg,gcomm_gridcell)
      call clm_writevar(nfid(t),'grid1d_lat',gptr%latdeg,gcomm_gridcell)

      ! Write landunit info

      do l = begl , endl
        rlarr(l) = gptr%londeg(lptr%gridcell(l))
      end do
      call clm_writevar(nfid(t),'land1d_lon',rlarr,gcomm_landunit)
      do l = begl , endl
        rlarr(l) = gptr%latdeg(lptr%gridcell(l))
      end do
      call clm_writevar(nfid(t),'land1d_lat',rlarr,gcomm_landunit)
      call clm_writevar(nfid(t),'land1d_wtgcell',lptr%wtgcell,gcomm_landunit)
      call clm_writevar(nfid(t),'land1d_ityplunit',lptr%itype,gcomm_landunit)
      call clm_writevar(nfid(t),'land1d_active',lptr%active,gcomm_landunit)

      ! Write column info

      do c = begc , endc
        rcarr(c) = gptr%londeg(cptr%gridcell(c))
      end do
      call clm_writevar(nfid(t),'cols1d_lon',rcarr,gcomm_column)
      do c = begc , endc
        rcarr(c) = gptr%latdeg(cptr%gridcell(c))
      end do
      call clm_writevar(nfid(t),'cols1d_lat',rcarr,gcomm_column)
      call clm_writevar(nfid(t),'cols1d_wtgcell',cptr%wtgcell,gcomm_column)
      call clm_writevar(nfid(t),'cols1d_wtlunit',cptr%wtlunit,gcomm_column)
      do c = begc , endc
        icarr(c) = lptr%itype(cptr%landunit(c))
      end do
      call clm_writevar(nfid(t),'cols1d_itype_lunit',icarr,gcomm_column)
      call clm_writevar(nfid(t),'cols1d_active',cptr%active,gcomm_column)

      ! Write pft info

      do p = begp , endp
        rparr(p) = gptr%londeg(pptr%gridcell(p))
      end do
      call clm_writevar(nfid(t),'pfts1d_lon',rparr,gcomm_pft)
      do p = begp , endp
        rparr(p) = gptr%latdeg(pptr%gridcell(p))
      end do
      call clm_writevar(nfid(t),'pfts1d_lat',rparr,gcomm_pft)
      call clm_writevar(nfid(t),'pfts1d_wtgcell',pptr%wtgcell,gcomm_pft)
      call clm_writevar(nfid(t),'pfts1d_wtlunit',pptr%wtlunit,gcomm_pft)
      call clm_writevar(nfid(t),'pfts1d_wtcol',pptr%wtcol,gcomm_pft)
      call clm_writevar(nfid(t),'pfts1d_itype_veg',pptr%itype,gcomm_pft)
      do p = begp , endp
        iparr(p) = lptr%itype(pptr%landunit(p))
      end do
      call clm_writevar(nfid(t),'pfts1d_itype_lunit',iparr,gcomm_pft)
      call clm_writevar(nfid(t),'pfts1d_active',pptr%active,gcomm_pft)

      deallocate(rgarr,rlarr,rcarr,rparr)
      deallocate(igarr,ilarr,icarr,iparr)

    end if
  end subroutine hfields_1dinfo
  !
  ! Write history tape(s)
  ! Determine if next time step is beginning of history interval and if so:
  !   increment the current time sample counter, open a new history file
  !   and if needed (i.e., when ntim = 1), write history data to current
  !   history file, reset field accumulation counters to zero.
  ! If primary history file is full or at the last time step of the simulation,
  !   write restart dataset and close all history fiels.
  ! If history file is full or at the last time step of the simulation:
  !   close history file
  !   and reset time sample counter to zero if file is full.
  ! Daily-averaged data for the first day in September are written on
  !   date = 00/09/02 with mscur = 0.
  ! Daily-averaged data for the first day in month mm are written on
  !   date = yyyy/mm/02 with mscur = 0.
  ! Daily-averaged data for the 30th day (last day in September) are written
  !   on date = 0000/10/01 mscur = 0.
  ! Daily-averaged data for the last day in month mm are written on
  !   date = yyyy/mm+1/01 with mscur = 0.
  !
  subroutine hist_htapes_wrapup( rstwr, nlend )
    implicit none
    logical , intent(in) :: rstwr    ! true => write restart file this step
    logical , intent(in) :: nlend    ! true => end of run on this step
    integer(ik4) :: t                ! tape index
    integer(ik4) :: f                ! field index
    integer(ik4) :: ier              ! error code
    integer(ik4) :: day              ! current day (1 -> 31)
    integer(ik4) :: mon              ! current month (1 -> 12)
    integer(ik4) :: yr               ! current year (0 -> ...)
    integer(ik4) :: mdcur            ! current day
    integer(ik4) :: mscur            ! seconds of current day
    integer(ik4) :: mcsec            ! current time of day [seconds]
    integer(ik4) :: daym1            ! nstep-1 day (1 -> 31)
    integer(ik4) :: monm1            ! nstep-1 month (1 -> 12)
    integer(ik4) :: yrm1             ! nstep-1 year (0 -> ...)
    integer(ik4) :: mcsecm1          ! nstep-1 time of day [seconds]
    integer(ik8) :: temp
    real(rk8):: time                 ! current time
    character(len=256) :: str        ! global attribute string
    logical :: if_stop               ! true => last time step of run
    ! true => write out 3D time-constant data
    logical , save :: do_3Dtconst = .true.
    character(len=*),parameter :: subname = 'hist_htapes_wrapup'

    ! Set calendar for current time step

    call curr_date(idatex, yr, mon, day, mcsec)
    call curr_time(idatex, mdcur, mscur)
    time = mdcur + mscur/secspday

    ! Set calendar for current for previous time step

    call get_prev_date(yrm1, monm1, daym1, mcsecm1)

    ! Loop over active history tapes, create new history files if necessary
    ! and write data to history files if end of history interval.
    do t = 1 , ntapes

      ! Skip nstep=0 if monthly average

      if ( ktau == 0 .and. tape(t)%nhtfrq == 0 ) cycle

      ! Determine if end of history interval
      tape(t)%is_endhist = .false.
      if ( tape(t)%nhtfrq == 0 ) then   !monthly average
        if ( mon /= monm1 ) tape(t)%is_endhist = .true.
      else
        temp = tape(t)%nhtfrq
        if ( mod(ktau,temp) == 0 ) tape(t)%is_endhist = .true.
      end if

      ! If end of history interval

      if ( tape(t)%is_endhist ) then

        ! Normalize history buffer if time averaged

        call hfields_normalize(t)

        ! Increment current time sample counter.

        tape(t)%ntimes = tape(t)%ntimes + 1

        ! Create history file if appropriate and build time comment

        ! If first time sample, generate unique history file name, open file,
        ! define dims, vars, etc.

        if ( tape(t)%ntimes == 1 ) then
          locfnh(t) = set_hist_filename(hist_freq=tape(t)%nhtfrq, &
                                        hist_mfilt=tape(t)%mfilt, &
                                        hist_file=t)
          if ( myid == italk ) then
            write(stdout,*) trim(subname), &
              ' : Creating history file ', trim(locfnh(t)), &
              ' at nstep = ',ktau
            write(stdout,*)'calling htape_create for file t = ',t
          end if
          call htape_create (t)

          ! Define time-constant field variables
          call htape_timeconst(t, mode='define')

          ! Define 3D time-constant field variables only to first primary tape

          if ( do_3Dtconst .and. t == 1 ) then
            call htape_timeconst3D(t, mode='define')
          end if

          ! Define model field variables

          call hfields_write(t, mode='define')

          ! Exit define model
          call clm_enddef(nfid(t))

        end if

        ! Write time constant history variables
        call htape_timeconst(t, mode='write')

        ! Write 3D time constant history variables only to first primary tape
        if ( do_3Dtconst .and. t == 1 .and. tape(t)%ntimes == 1 ) then
          call htape_timeconst3D(t, mode='write')
          do_3Dtconst = .false.
        end if

        if ( myid == italk ) then
          write(stdout,*)
          write(stdout,*) trim(subname), &
            ' : Writing current time sample to local history file ', &
            trim(locfnh(t)),' at nstep = ',ktau, &
            ' for history time interval beginning at ', tape(t)%begtime, &
            ' and ending at ',time
          write(stdout,*)
        end if

        ! Update beginning time of next interval

        tape(t)%begtime = time

        ! Write history time samples

        call hfields_write(t, mode='write')

        ! Zero necessary history buffers

        call hfields_zero(t)

      end if
    end do  ! end loop over history tapes

    ! Determine if file needs to be closed

    tapes_ntimes = tape(:)%ntimes
    tapes_mfilt = tape(:)%mfilt
    call hist_do_disp(ntapes, tapes_ntimes , tapes_mfilt, &
                      if_stop, if_disphist, rstwr, nlend)

    ! Close open history file
    ! Auxilary files may have been closed and saved off without being full,
    ! must reopen the files

    do t = 1 , ntapes
      if ( if_disphist(t) ) then
        if ( tape(t)%ntimes /= 0 ) then
          if ( myid == italk ) then
            write(stdout,*)
            write(stdout,*)  trim(subname),' : Closing local history file ',&
               trim(locfnh(t)),' at nstep = ', ktau
            write(stdout,*)
          end if
          call clm_closefile(nfid(t))
          if ( .not. if_stop .and. (tape(t)%ntimes/=tape(t)%mfilt) ) then
            call clm_openfile(trim(locfnh(t)), nfid(t), clm_readwrite)
          end if
        else
          if ( myid == italk ) then
            write(stdout,*) trim(subname), &
              ' : history tape ',t,': no open file to close'
          end if
        end if
      end if
    end do

    ! Reset number of time samples to zero if file is full

    do t = 1 , ntapes
      if ( if_disphist(t) .and. tape(t)%ntimes==tape(t)%mfilt ) then
        tape(t)%ntimes = 0
      end if
    end do
  end subroutine hist_htapes_wrapup
  !
  ! Read/write history file restart data.
  ! If the current history file(s) are not full, file(s) are opened
  ! so that subsequent time samples are added until the file is full.
  !
  subroutine hist_restart_ncd (ncid, flag, rdate)
    implicit none
    type(clm_filetype), intent(inout) :: ncid  ! netcdf file
    character(len=*) , intent(in) :: flag      !'read' or 'write'
    ! restart file time stamp for name
    character(len=*) , intent(in) , optional :: rdate
    integer(ik4) :: max_nflds ! Max number of fields
    ! 1d size, beginning and ending indices
    integer(ik4) :: num1d , beg1d , end1d
    ! 1d size, beginning and ending indices
    integer(ik4) :: num1d_out , beg1d_out , end1d_out
    ! 2d size (e.g. number of vertical levels)
    integer(ik4) :: num2d
    integer(ik4) :: begp , endp ! per-proc beginning and ending pft indices
    integer(ik4) :: begc , endc ! per-proc beginning and ending column indices
    integer(ik4) :: begl , endl ! per-proc beginning and ending landunit indices
    integer(ik4) :: begg , endg ! per-proc gridcell ending gridcell indices
    integer(ik4) :: numa      ! total number of atm cells across all processors
    integer(ik4) :: numg      ! total number of gridcells across all processors
    integer(ik4) :: numl      ! total number of landunits across all processors
    integer(ik4) :: numc      ! total number of columns across all processors
    integer(ik4) :: nump      ! total number of pfts across all processors
    character(len=max_namlen) :: name       ! variable name
    character(len=max_namlen) :: name_acc   ! accumulator variable name
    character(len=max_namlen) :: long_name  ! long name of variable
    character(len=max_chars) :: long_name_acc   ! long name for accumulator
    character(len=max_chars) :: units           ! units of variable
    character(len=max_chars) :: units_acc       ! accumulator units
    character(len=max_chars) :: fname           ! full name of history file
    ! local history restart file names
    character(len=max_chars) :: locrest(max_tapes)

    character(len=max_namlen) , allocatable :: tname(:)
    character(len=max_chars) , allocatable :: tunits(:) , tlongname(:)
    character(len=8) , allocatable :: tmpstr(:,:)
    character(len=1) , allocatable :: tavgflag(:)
    integer(ik4) :: start(2)

    character(len=1) :: hnum            ! history file index
    character(len=8) :: type1d          ! clm pointer 1d type
    character(len=8) :: type1d_out      ! history buffer 1d type
    character(len=8) :: type2d          ! history buffer 2d type
    character(len=32) :: dim1name       ! temporary
    character(len=32) :: dim2name       ! temporary
    character(len=32) :: dim3name       ! temporary
    character(len=32) :: dim4name       ! temporary
    integer(ik4) :: status ! error status
    integer(ik4) :: dimid  ! dimension ID
    integer(ik4) :: k      ! 1d index
    ! number of history tapes on the restart file
    integer(ik4) :: ntapes_onfile
    ! number of history fields on the restart file
    integer(ik4) :: nflds_onfile
    integer(ik4) :: t  ! tape index
    integer(ik4) :: f  ! field index
    integer(ik4) :: varid  ! variable id
    integer(ik4) , allocatable :: itemp2d(:,:)     ! 2D temporary
    real(rk8) , pointer :: hbuf(:,:)               ! history buffer
    real(rk8) , pointer :: hbuf1d(:)               ! 1d history buffer
    integer(ik4) , pointer :: nacs(:,:)            ! accumulation counter
    integer(ik4) , pointer :: nacs1d(:)            ! 1d accumulation counter
    character(len=*) , parameter :: subname = 'hist_restart_ncd'
    type(subgrid_type) , pointer :: gcomm

    if ( flag == 'read' ) then
      ! If startup run just return
      if ( nsrest == nsrStartup ) then
        return
      end if
    end if

    ! Read history file data only for restart run
    !
    ! First when writing out and in define mode, create files and define
    ! all variables
    !
    !================================================
    if ( flag == 'define' ) then
    !================================================

      if (.not. present(rdate)) then
        call fatal(__FILE__,__LINE__, &
                  'variable rdate must be present for writing restart files')
      end if

      !
      ! On master restart file add ntapes/max_chars dimension
      ! and then add the history and history restart filenames
      !
      call clm_adddim(ncid, 'ntapes', ntapes)
      call clm_adddim(ncid, 'max_chars', max_chars)

      call clm_addvar(clmvar_text,ncid,'locfnh', &
              cdims=(/'max_chars','ntapes   '/), &
              long_name='History filename',      &
              comment='This variable NOT needed for startup&
                     & simulations')
      call clm_addvar(clmvar_text,ncid,'locfnhr', &
              cdims=(/'max_chars','ntapes   '/), &
              long_name='Restart history filename',      &
              comment='This variable NOT needed for startup&
                     & simulations')

      ! max_nflds is the maximum number of fields on any tape
      ! max_flds is the maximum number possible number of fields

      max_nflds = max_nFields()

      call get_proc_global(numg,numl,numc,nump)

      ! Loop over tapes - write out namelist information to each
      ! restart-history tape only read/write accumulators and
      ! counters if needed

      do t = 1 , ntapes

        !
        ! Create the restart history filename and open it
        !
        write(hnum,'(i1.1)') t-1
        locfnhr(t) = trim(dirout)//trim(caseid)//".clm."//trim(inst_suffix) &
                        // ".rh" // hnum //"."// trim(rdate) //".nc"

        call htape_create( t, histrest=.true. )

        !
        ! Add read/write accumultators and counters if needed
        !
        if ( .not. tape(t)%is_endhist ) then
          do f = 1 , tape(t)%nflds
            name          =  tape(t)%hlist(f)%field%name
            long_name     =  tape(t)%hlist(f)%field%long_name(1:max_namlen)
            units         =  tape(t)%hlist(f)%field%units
            name_acc      =  trim(name) // "_acc"
            units_acc     =  "unitless positive integer(ik4)"
            long_name_acc =  trim(long_name) // " accumulator number of samples"
            type1d_out    =  tape(t)%hlist(f)%field%type1d_out
            type2d        =  tape(t)%hlist(f)%field%type2d
            num2d         =  tape(t)%hlist(f)%field%num2d
            nacs          => tape(t)%hlist(f)%nacs
            hbuf          => tape(t)%hlist(f)%hbuf

            if (type1d_out == grlnd) then
              dim1name = trim(grlnd)
              dim2name = 'undefined'
            else
              dim1name = type1d_out
              dim2name = 'undefined'
            end if

            if ( dim2name == 'undefined' ) then
              if ( num2d == 1 ) then
                call clm_addvar(clmvar_double,ncid_hist(t),trim(name), &
                  cdims=(/dim1name/),long_name=trim(long_name), &
                  units=trim(units))
                call clm_addvar(clmvar_integer,ncid_hist(t),trim(name_acc), &
                  cdims=(/dim1name/),long_name=trim(long_name_acc), &
                  units=trim(units_acc))
              else
                dim2name = type2d
                call clm_addvar(clmvar_double,ncid_hist(t),trim(name), &
                  cdims=(/dim1name,dim2name/),long_name=trim(long_name), &
                  units=trim(units))
                call clm_addvar(clmvar_integer,ncid_hist(t),trim(name_acc), &
                  cdims=(/dim1name,dim2name/),long_name=trim(long_name_acc), &
                  units=trim(units_acc))
              end if
            else
              if ( num2d == 1 ) then
                call clm_addvar(clmvar_double,ncid_hist(t),trim(name), &
                  cdims=(/dim1name,dim2name/),long_name=trim(long_name), &
                  units=trim(units))
                call clm_addvar(clmvar_integer,ncid_hist(t),trim(name_acc), &
                  cdims=(/dim1name,dim2name/),long_name=trim(long_name_acc), &
                  units=trim(units_acc))
              else
                dim3name = type2d
                call clm_addvar(clmvar_double,ncid_hist(t),trim(name), &
                  cdims=(/dim1name,dim2name,dim2name/), &
                  long_name=trim(long_name), units=trim(units))
                call clm_addvar(clmvar_integer,ncid_hist(t),trim(name_acc), &
                  cdims=(/dim1name,dim2name,dim2name/), &
                  long_name=trim(long_name_acc), units=trim(units_acc))
              end if
            end if
          end do
        end if

        !
        ! Add namelist information to each restart history tape
        !
        call clm_adddim(ncid_hist(t),'fname_lenp2',max_namlen+2)
        call clm_adddim(ncid_hist(t),'fname_len',max_namlen)
        call clm_adddim(ncid_hist(t),'len1',1)
        call clm_adddim(ncid_hist(t),'scalar',1)
        call clm_adddim(ncid_hist(t),'max_chars',max_chars)
        call clm_adddim(ncid_hist(t),'max_nflds',max_nflds)
        call clm_adddim(ncid_hist(t),'max_flds',max_flds)

        call clm_addvar(clmvar_integer,ncid_hist(t),'nhtfrq', &
                  long_name="Frequency of history writes",    &
                  comment="Namelist item", &
                  units="absolute value of negative is in hours, 0=monthly,"&
                       &" positive is time-steps")
        call clm_addvar(clmvar_integer,ncid_hist(t),'mfilt', &
                  long_name="Number of history time samples on a file", &
                  comment="Namelist item")
        call clm_addvar(clmvar_integer,ncid_hist(t),'ncprec', &
                  long_name="Flag for data precision", &
                  flag_values=(/1,2/),valid_range=(/1,2/), &
                  flag_meanings=(/"single-precision", "double-precision"/))
        call clm_addvar(clmvar_logical,ncid_hist(t),'dov2xy', &
                  long_name="Output on 2D grid format (TRUE) or vector"&
                           &" format (FALSE)", comment="Namelist item")
        call clm_addvar(clmvar_text,ncid_hist(t),'fincl', &
                  cdims=(/'fname_lenp2','max_flds   '/), &
                  long_name="Fieldnames to include",comment="Namelist item")
        call clm_addvar(clmvar_text,ncid_hist(t),'fexcl', &
                  cdims=(/'fname_lenp2','max_flds   '/), &
                  long_name="Fieldnames to exclude",comment="Namelist item")
        call clm_addvar(clmvar_integer,ncid_hist(t),'nflds', &
                  long_name="Number of fields on file")
        call clm_addvar(clmvar_integer,ncid_hist(t),'ntimes', &
                  long_name="Number of time steps on file")
        call clm_addvar(clmvar_logical,ncid_hist(t),'is_endhist', &
                long_name="End of history file")
        call clm_addvar(clmvar_double,ncid_hist(t),'begtime', &
                long_name="Beginning time")
        call clm_addvar(clmvar_integer,ncid_hist(t),'num2d', &
                cdims=(/'max_nflds'/), long_name="Size of second dimension")
        call clm_addvar(clmvar_integer,ncid_hist(t),'hpindex', &
                cdims=(/'max_nflds'/), long_name="History pointer index")
        call clm_addvar(clmvar_text,ncid_hist(t),'avgflag', &
                cdims=(/'len1     ','max_nflds'/), &
                long_name="Averaging flag", &
                comment="A=Average, X=Maximum, M=Minimum, I=Instantaneous")
        call clm_addvar(clmvar_text,ncid_hist(t),'name', &
                cdims=(/'fname_len','max_nflds'/), long_name="Fieldnames")
        call clm_addvar(clmvar_text,ncid_hist(t),'long_name', &
                cdims=(/'max_chars','max_nflds'/), &
                long_name="Long descriptive names for fields")
        call clm_addvar(clmvar_text,ncid_hist(t),'units', &
                cdims=(/'max_chars','max_nflds'/), &
                long_name="Units for each history field output")
        call clm_addvar(clmvar_text,ncid_hist(t),'type1d', &
                cdims=(/'string_length','max_nflds    '/), &
                long_name="1st dimension type")
        call clm_addvar(clmvar_text,ncid_hist(t),'type1d_out', &
                cdims=(/'string_length','max_nflds    '/), &
                long_name="1st output dimension type")
        call clm_addvar(clmvar_text,ncid_hist(t),'type2d', &
                cdims=(/'string_length','max_nflds    '/), &
                long_name="2nd dimension type")
        call clm_addvar(clmvar_text,ncid_hist(t),'p2c_scale_type', &
                cdims=(/'string_length','max_nflds    '/), &
                long_name="PFT to column scale type")
        call clm_addvar(clmvar_text,ncid_hist(t),'c2l_scale_type', &
                cdims=(/'string_length','max_nflds    '/), &
                long_name="column to landunit scale type")
        call clm_addvar(clmvar_text,ncid_hist(t),'l2g_scale_type', &
                cdims=(/'string_length','max_nflds    '/), &
                long_name="landunit to gridpoint scale type")

        call clm_enddef(ncid_hist(t))

      end do   ! end of ntapes loop

      return

    !================================================
    else if ( flag == 'write' ) then
    !================================================

      !
      ! First write out namelist information to each restart history file
      !
      ! Add history filenames to master restart file
      do t = 1 , ntapes
        call clm_writevar(ncid,'locfnh',locfnh(1:ntapes))
        call clm_writevar(ncid,'locfnhr',locfnhr(1:ntapes))
      end do

      fincl(:,1) = hist_fincl1(:)
      fincl(:,2) = hist_fincl2(:)
      fincl(:,3) = hist_fincl3(:)
      fincl(:,4) = hist_fincl4(:)
      fincl(:,5) = hist_fincl5(:)
      fincl(:,6) = hist_fincl6(:)

      fexcl(:,1) = hist_fexcl1(:)
      fexcl(:,2) = hist_fexcl2(:)
      fexcl(:,3) = hist_fexcl3(:)
      fexcl(:,4) = hist_fexcl4(:)
      fexcl(:,5) = hist_fexcl5(:)
      fexcl(:,6) = hist_fexcl6(:)

      max_nflds = max_nFields()

      start(1)=1

      allocate(itemp2d(max_nflds,ntapes))

      !
      ! Add history namelist data to each history restart tape
      !
      do t = 1 , ntapes
        call clm_writevar(ncid_hist(t),'fincl',fincl(:,t))
        call clm_writevar(ncid_hist(t),'fexcl',fexcl(:,t))
        call clm_writevar(ncid_hist(t),'is_endhist',tape(t)%is_endhist)
        call clm_writevar(ncid_hist(t),'dov2xy',tape(t)%dov2xy)
        itemp2d(:,:) = 0
        do f = 1 , tape(t)%nflds
          itemp2d(f,t) = tape(t)%hlist(f)%field%num2d
        end do
        call clm_writevar(ncid_hist(t),'num2d',itemp2d(:,t))
        itemp2d(:,:) = 0
        do f = 1 , tape(t)%nflds
          itemp2d(f,t) = tape(t)%hlist(f)%field%hpindex
        end do
        call clm_writevar(ncid_hist(t),'hpindex',itemp2d(:,t))
        call clm_writevar(ncid_hist(t),'nflds',tape(t)%nflds)
        call clm_writevar(ncid_hist(t),'ntimes',tape(t)%ntimes)
        call clm_writevar(ncid_hist(t),'nhtfrq',tape(t)%nhtfrq)
        call clm_writevar(ncid_hist(t),'mfilt',tape(t)%mfilt)
        call clm_writevar(ncid_hist(t),'ncprec',tape(t)%ncprec)
        call clm_writevar(ncid_hist(t),'begtime',tape(t)%begtime)
        allocate(tmpstr(tape(t)%nflds,6 ),tname(tape(t)%nflds), &
                 tavgflag(tape(t)%nflds),tunits(tape(t)%nflds), &
                 tlongname(tape(t)%nflds))
        do f = 1 , tape(t)%nflds
          tname(f)  = tape(t)%hlist(f)%field%name
          tunits(f) = tape(t)%hlist(f)%field%units
          tlongname(f) = tape(t)%hlist(f)%field%long_name
          tmpstr(f,1) = tape(t)%hlist(f)%field%type1d
          tmpstr(f,2) = tape(t)%hlist(f)%field%type1d_out
          tmpstr(f,3) = tape(t)%hlist(f)%field%type2d
          tavgflag(f) = tape(t)%hlist(f)%avgflag
          tmpstr(f,4) = tape(t)%hlist(f)%field%p2c_scale_type
          tmpstr(f,5) = tape(t)%hlist(f)%field%c2l_scale_type
          tmpstr(f,6) = tape(t)%hlist(f)%field%l2g_scale_type
        end do
        call clm_writevar(ncid_hist(t),'name',tname)
        call clm_writevar(ncid_hist(t),'long_name',tlongname)
        call clm_writevar(ncid_hist(t),'units',tunits)
        call clm_writevar(ncid_hist(t),'type1d',tmpstr(:,1))
        call clm_writevar(ncid_hist(t),'type1d_out',tmpstr(:,2))
        call clm_writevar(ncid_hist(t),'type2d',tmpstr(:,3))
        call clm_writevar(ncid_hist(t),'avgflag',tavgflag)
        call clm_writevar(ncid_hist(t),'p2c_scale_type',tmpstr(:,4))
        call clm_writevar(ncid_hist(t),'c2l_scale_type',tmpstr(:,5))
        call clm_writevar(ncid_hist(t),'l2g_scale_type',tmpstr(:,6))
        deallocate(tname,tlongname,tunits,tmpstr,tavgflag)
      end do
      deallocate(itemp2d)

    !================================================
    else if ( flag == 'read' ) then
    !================================================
      !
      ! Read in namelist information
      !
      call clm_inqdim(ncid,'ntapes',ntapes_onfile)
      if ( ktau > 0 .and. ntapes_onfile /= ntapes ) then
        write(stderr,*) &
            'ntapes = ', ntapes, ' ntapes_onfile = ', ntapes_onfile
        call fatal(__FILE__,__LINE__, &
          trim(subname)//' ERROR: number of ntapes different '// &
          'than on restart file!, you can NOT change history options on '// &
          'restart!' )
      end if
      if ( ktau > 0 .and. ntapes > 0 ) then
        call clm_readvar(ncid,'locfnh',locfnh(1:ntapes))
        call clm_readvar(ncid,'locfnhr',locrest(1:ntapes))
        do t = 1,ntapes
          call strip_null(locrest(t))
          call strip_null(locfnh(t))
        end do
      end if

      ! Determine necessary indices - the following is needed if model
      ! decomposition is different on restart

      call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)
      call get_proc_global(numg,numl,numc,nump)

      if ( ktau > 0 ) then
        do t = 1 , ntapes
          call clm_openfile(locrest(t),ncid_hist(t))
          if ( t == 1 ) then
            call clm_inqdim(ncid_hist(t),'max_nflds',max_nflds)
            allocate(itemp2d(max_nflds,ntapes))
          end if
          call clm_readvar(ncid_hist(t),'fincl',fincl(:,t))
          call clm_readvar(ncid_hist(t),'fexcl',fexcl(:,t))
          call clm_readvar(ncid_hist(t),'nflds',nflds_onfile)
          if ( nflds_onfile /= tape(t)%nflds ) then
            write(stderr,*) &
               'nflds = ', tape(t)%nflds, ' nflds_onfile = ', nflds_onfile
            call fatal(__FILE__,__LINE__, &
                  trim(subname)//' ERROR: number of fields '// &
                  'different than on restart file!, you can NOT change '// &
                  'history options on restart!' )
          end if
          call clm_readvar(ncid_hist(t),'ntimes',tape(t)%ntimes)
          call clm_readvar(ncid_hist(t),'nhtfrq',tape(t)%nhtfrq)
          call clm_readvar(ncid_hist(t),'mfilt',tape(t)%mfilt)
          call clm_readvar(ncid_hist(t),'ncprec',tape(t)%ncprec)
          call clm_readvar(ncid_hist(t),'begtime',tape(t)%begtime)
          call clm_readvar(ncid_hist(t),'is_endhist',tape(t)%is_endhist)
          call clm_readvar(ncid_hist(t),'dov2xy',tape(t)%dov2xy)
          call clm_readvar(ncid_hist(t),'num2d',itemp2d(:,t))
          do f = 1 , tape(t)%nflds
            tape(t)%hlist(f)%field%num2d = itemp2d(f,t)
          end do
          call clm_readvar(ncid_hist(t),'hpindex',itemp2d(:,t))
          do f = 1 , tape(t)%nflds
            tape(t)%hlist(f)%field%hpindex = itemp2d(f,t)
          end do
          do f = 1 , tape(t)%nflds
            call clm_readvar(ncid_hist(t),'name',tape(t)%hlist(f)%field%name,f)
            call clm_readvar(ncid_hist(t),'long_name', &
                    tape(t)%hlist(f)%field%long_name,f)
            call clm_readvar(ncid_hist(t),'units', &
                    tape(t)%hlist(f)%field%units,f)
            call clm_readvar(ncid_hist(t),'type1d', &
                    tape(t)%hlist(f)%field%type1d,f)
            call clm_readvar(ncid_hist(t),'type1d_out', &
                    tape(t)%hlist(f)%field%type1d_out,f)
            call clm_readvar(ncid_hist(t),'type2d', &
                    tape(t)%hlist(f)%field%type2d,f)
            call clm_readvar(ncid_hist(t),'avgflag', &
                    tape(t)%hlist(f)%avgflag,f)
            call clm_readvar(ncid_hist(t),'p2c_scale_type', &
                    tape(t)%hlist(f)%field%p2c_scale_type,f)
            call clm_readvar(ncid_hist(t),'c2l_scale_type', &
                    tape(t)%hlist(f)%field%c2l_scale_type,f)
            call clm_readvar(ncid_hist(t),'l2g_scale_type', &
                    tape(t)%hlist(f)%field%l2g_scale_type,f)
            call strip_null(tape(t)%hlist(f)%field%name)
            call strip_null(tape(t)%hlist(f)%field%long_name)
            call strip_null(tape(t)%hlist(f)%field%units)
            call strip_null(tape(t)%hlist(f)%field%type1d)
            call strip_null(tape(t)%hlist(f)%field%type1d_out)
            call strip_null(tape(t)%hlist(f)%field%type2d)
            call strip_null(tape(t)%hlist(f)%field%p2c_scale_type)
            call strip_null(tape(t)%hlist(f)%field%c2l_scale_type)
            call strip_null(tape(t)%hlist(f)%field%l2g_scale_type)
            call strip_null(tape(t)%hlist(f)%avgflag)
            type1d_out = trim(tape(t)%hlist(f)%field%type1d_out)
            select case (trim(type1d_out))
              case (grlnd)
                num1d_out = numg
                beg1d_out = begg
                end1d_out = endg
              case (nameg)
                num1d_out = numg
                beg1d_out = begg
                end1d_out = endg
              case (namel)
                num1d_out = numl
                beg1d_out = begl
                end1d_out = endl
              case (namec)
                num1d_out = numc
                beg1d_out = begc
                end1d_out = endc
              case (namep)
                num1d_out = nump
                beg1d_out = begp
                end1d_out = endp
              case default
                write(stderr,*) trim(subname), &
                       ' ERROR: read unknown 1d output type=',trim(type1d_out)
                call fatal(__FILE__,__LINE__,'clm now stopping.')
            end select

            tape(t)%hlist(f)%field%num1d_out = num1d_out
            tape(t)%hlist(f)%field%beg1d_out = beg1d_out
            tape(t)%hlist(f)%field%end1d_out = end1d_out

            num2d  = tape(t)%hlist(f)%field%num2d
            allocate (tape(t)%hlist(f)%hbuf(beg1d_out:end1d_out,num2d), &
                      tape(t)%hlist(f)%nacs(beg1d_out:end1d_out,num2d), &
                      stat=status)
            if (status /= 0) then
              write(stderr,*) trim(subname), &
                    ' ERROR: allocation error for hbuf,nacs at t,f=',t,f
              call fatal(__FILE__,__LINE__,'clm now stopping.')
            end if
            tape(t)%hlist(f)%hbuf(:,:) = 0.D0
            tape(t)%hlist(f)%nacs(:,:) = 0

            type1d = tape(t)%hlist(f)%field%type1d
            select case (type1d)
              case (grlnd)
                num1d = numg
                beg1d = begg
                end1d = endg
              case (nameg)
                num1d = numg
                beg1d = begg
                end1d = endg
              case (namel)
                num1d = numl
                beg1d = begl
                end1d = endl
              case (namec)
                num1d = numc
                beg1d = begc
                end1d = endc
              case (namep)
                num1d = nump
                beg1d = begp
                end1d = endp
              case default
                 write(stderr,*) trim(subname), &
                           ' ERROR: read unknown 1d type=',type1d
                 call fatal(__FILE__,__LINE__,'clm now stopping.')
            end select

            tape(t)%hlist(f)%field%num1d = num1d
            tape(t)%hlist(f)%field%beg1d = beg1d
            tape(t)%hlist(f)%field%end1d = end1d

          end do   ! end of flds loop

          ! If history file is not full, open it

          if (tape(t)%ntimes /= 0) then
            call clm_openfile(trim(locfnh(t)),nfid(t),clm_readwrite)
          end if

        end do  ! end of tapes loop

        hist_fincl1(:) = fincl(:,1)
        hist_fincl2(:) = fincl(:,2)
        hist_fincl3(:) = fincl(:,3)
        hist_fincl4(:) = fincl(:,4)
        hist_fincl5(:) = fincl(:,5)
        hist_fincl6(:) = fincl(:,6)

        hist_fexcl1(:) = fexcl(:,1)
        hist_fexcl2(:) = fexcl(:,2)
        hist_fexcl3(:) = fexcl(:,3)
        hist_fexcl4(:) = fexcl(:,4)
        hist_fexcl5(:) = fexcl(:,5)
        hist_fexcl6(:) = fexcl(:,6)

      end if

      if ( allocated(itemp2d) ) deallocate(itemp2d)

    end if

    !======================================================================
    ! Read/write history file restart data.
    ! If the current history file(s) are not full, file(s) are opened
    ! so that subsequent time samples are added until the file is full.
    !======================================================================

    if ( flag == 'write' ) then

      do t = 1 , ntapes
        if ( .not. tape(t)%is_endhist ) then
          do f = 1 , tape(t)%nflds
            name       =  tape(t)%hlist(f)%field%name
            name_acc   =  trim(name) // "_acc"
            type1d_out =  tape(t)%hlist(f)%field%type1d_out
            gcomm      => tape(t)%hlist(f)%field%gcomm
            type2d     =  tape(t)%hlist(f)%field%type2d
            num2d      =  tape(t)%hlist(f)%field%num2d
            beg1d_out  =  tape(t)%hlist(f)%field%beg1d_out
            end1d_out  =  tape(t)%hlist(f)%field%end1d_out
            nacs       => tape(t)%hlist(f)%nacs
            hbuf       => tape(t)%hlist(f)%hbuf

            if ( num2d == 1 ) then
              allocate(hbuf1d(beg1d_out:end1d_out), &
                       nacs1d(beg1d_out:end1d_out), stat=status)
              if (status /= 0) then
                write(stderr,*) trim(subname),' ERROR: allocation'
                call fatal(__FILE__,__LINE__,'clm now stopping.')
              end if
              hbuf1d(beg1d_out:end1d_out) = hbuf(beg1d_out:end1d_out,1)
              nacs1d(beg1d_out:end1d_out) = nacs(beg1d_out:end1d_out,1)
              call clm_writevar(ncid_hist(t),trim(name),hbuf1d,gcomm)
              call clm_writevar(ncid_hist(t),trim(name_acc),nacs1d,gcomm)
              deallocate(hbuf1d)
              deallocate(nacs1d)
            else
              call clm_writevar(ncid_hist(t),trim(name),hbuf,gcomm)
              call clm_writevar(ncid_hist(t),trim(name_acc),nacs,gcomm)
            end if
          end do
        end if  ! end of is_endhist block
        call clm_closefile(ncid_hist(t))
      end do   ! end of ntapes loop

    else if (flag == 'read') then

      ! Read history restart information if history files are not full

      do t = 1 , ntapes
        if ( .not. tape(t)%is_endhist ) then
          do f = 1 , tape(t)%nflds
            name       =  tape(t)%hlist(f)%field%name
            name_acc   =  trim(name) // "_acc"
            type1d_out =  tape(t)%hlist(f)%field%type1d_out
            gcomm      => tape(t)%hlist(f)%field%gcomm
            type2d     =  tape(t)%hlist(f)%field%type2d
            num2d      =  tape(t)%hlist(f)%field%num2d
            beg1d_out  =  tape(t)%hlist(f)%field%beg1d_out
            end1d_out  =  tape(t)%hlist(f)%field%end1d_out
            nacs       => tape(t)%hlist(f)%nacs
            hbuf       => tape(t)%hlist(f)%hbuf

            if ( num2d == 1 ) then
              allocate(hbuf1d(beg1d_out:end1d_out), &
                       nacs1d(beg1d_out:end1d_out), stat=status)
              if (status /= 0) then
                write(stderr,*) trim(subname),' ERROR: allocation'
                call fatal(__FILE__,__LINE__,'clm now stopping.')
              end if
              call clm_readvar(ncid_hist(t),trim(name),hbuf1d,gcomm)
              call clm_readvar(ncid_hist(t),trim(name_acc),nacs1d,gcomm)
              hbuf(beg1d_out:end1d_out,1) = hbuf1d(beg1d_out:end1d_out)
              nacs(beg1d_out:end1d_out,1) = nacs1d(beg1d_out:end1d_out)
              deallocate(hbuf1d)
              deallocate(nacs1d)
            else
              call clm_readvar(ncid_hist(t),trim(name),hbuf,gcomm)
              call clm_readvar(ncid_hist(t),trim(name_acc),nacs,gcomm)
            end if
          end do
        end if
        call clm_closefile(ncid_hist(t))
      end do
    end if
  end subroutine hist_restart_ncd
  !
  ! Initialize a single level history field. The pointer, ptrhist,
  ! is a pointer to the clmtype array that the history buffer will use.
  ! The value of type1d passed to masterlist\_add\_fld determines which of the
  ! 1d type of the output and the beginning and ending indices the history
  ! buffer field). Default history contents for given field on all tapes
  ! are set by calling [masterlist\_make\_active] for the appropriate tape.
  ! After the masterlist is built, routine [htapes\_build] is called for an
  ! initial run to initialize the actual history tapes.
  !
  subroutine hist_addfld1d (fname, units, avgflag, long_name, type1d_out, &
                        ptr_gcell, ptr_lunit, ptr_col, ptr_pft, ptr_lnd, &
                        ptr_atm, p2c_scale_type, c2l_scale_type, &
                        l2g_scale_type, set_lake, set_nolake, set_urb, &
                        set_nourb, set_spec, default)
    implicit none
    character(len=*) , intent(in) :: fname       ! field name
    character(len=*) , intent(in) :: units       ! units of field
    character(len=1) , intent(in) :: avgflag     ! time averaging flag
    character(len=*) , intent(in)  :: long_name  ! long name of field
    ! output type (from clmtype)
    character(len=*) , optional , intent(in) :: type1d_out
    real(rk8) , optional , pointer :: ptr_gcell(:)   ! pointer to gridcell array
    real(rk8) , optional , pointer :: ptr_lunit(:)   ! pointer to landunit array
    real(rk8) , optional , pointer :: ptr_col(:)     ! pointer to column array
    real(rk8) , optional , pointer :: ptr_pft(:)     ! pointer to pft array
    real(rk8) , optional , pointer :: ptr_lnd(:)     ! pointer to lnd array
    real(rk8) , optional , pointer :: ptr_atm(:)     ! pointer to atm array
    real(rk8) , optional , intent(in) :: set_lake    ! value to set lakes to
    real(rk8) , optional , intent(in) :: set_nolake  ! value to set non-lakes to
    real(rk8) , optional , intent(in) :: set_urb     ! value to set urban to
    real(rk8) , optional , intent(in) :: set_nourb   ! value to set non-urban to
    ! value to set non-glacier_mec to
    real(rk8) , optional , intent(in) :: set_spec     ! value to set special to
    ! scale type for subgrid averaging of pfts to column
    character(len=*) , optional , intent(in) :: p2c_scale_type
    ! scale type for subgrid averaging of columns to landunits
    character(len=*) , optional , intent(in) :: c2l_scale_type
    ! scale type for subgrid averaging of landunits to gridcells
    character(len=*) , optional , intent(in) :: l2g_scale_type
    ! if set to 'inactive, field will not appear on primary tape
    character(len=*) , optional , intent(in) :: default
    integer(ik4) :: p , c , l , g  ! indices
    integer(ik4) :: begp , endp ! per-proc beginning and ending pft indices
    integer(ik4) :: begc , endc ! per-proc beginning and ending column indices
    integer(ik4) :: begl , endl ! per-proc beginning and ending landunit indices
    integer(ik4) :: begg , endg ! per-proc gridcell ending gridcell indices
    integer(ik4) :: hpindex     ! history buffer pointer index
    character(len=8) :: l_type1d       ! 1d data type
    character(len=8) :: l_type1d_out   ! 1d output type
    ! scale type for subgrid averaging of pfts to column
    character(len=8) :: scale_type_p2c
    ! scale type for subgrid averaging of columns to landunits
    character(len=8) :: scale_type_c2l
    ! scale type for subgrid averaging of landunits to gridcells
    character(len=8) :: scale_type_l2g
    character(len=*),parameter :: subname = 'hist_addfld1d'

    ! Determine processor bounds

    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

    ! History buffer pointer

    hpindex = pointer_index()

    if ( present(ptr_lnd) ) then
      l_type1d = grlnd(1:8)
      l_type1d_out = grlnd(1:8)
      clmptr_rs(hpindex)%ptr => ptr_lnd
    else if ( present(ptr_gcell) ) then
      l_type1d = nameg(1:8)
      l_type1d_out = nameg(1:8)
      clmptr_rs(hpindex)%ptr => ptr_gcell
    else if ( present(ptr_lunit) ) then
      l_type1d = namel(1:8)
      l_type1d_out = namel(1:8)
      clmptr_rs(hpindex)%ptr => ptr_lunit
      if ( present(set_lake) ) then
        do l = begl , endl
          if ( clm3%g%l%lakpoi(l) ) ptr_lunit(l) = set_lake
        end do
      end if
      if ( present(set_nolake) ) then
        do l = begl , endl
          if ( .not. (clm3%g%l%lakpoi(l)) ) ptr_lunit(l) = set_nolake
        end do
      end if
      if ( present(set_urb) ) then
        do l = begl , endl
          if ( clm3%g%l%urbpoi(l) ) ptr_lunit(l) = set_urb
        end do
      end if
      if ( present(set_nourb) ) then
        do l = begl , endl
          if ( .not. (clm3%g%l%urbpoi(l)) ) ptr_lunit(l) = set_nourb
        end do
      end if
      if ( present(set_spec) ) then
        do l = begl , endl
          if ( clm3%g%l%ifspecial(l) ) ptr_lunit(l) = set_spec
        end do
      end if
    else if ( present(ptr_col) ) then
      l_type1d = namec(1:8)
      l_type1d_out = namec(1:8)
      clmptr_rs(hpindex)%ptr => ptr_col
      if ( present(set_lake) ) then
        do c = begc , endc
          l = clm3%g%l%c%landunit(c)
          if ( clm3%g%l%lakpoi(l) ) ptr_col(c) = set_lake
        end do
      end if
      if ( present(set_nolake) ) then
        do c = begc , endc
          l = clm3%g%l%c%landunit(c)
          if ( .not. (clm3%g%l%lakpoi(l)) ) ptr_col(c) = set_nolake
        end do
      end if
      if ( present(set_urb) ) then
        do c = begc , endc
          l = clm3%g%l%c%landunit(c)
          if ( clm3%g%l%urbpoi(l) ) ptr_col(c) = set_urb
        end do
      end if
      if ( present(set_nourb) ) then
        do c = begc , endc
          l = clm3%g%l%c%landunit(c)
          if ( .not. (clm3%g%l%urbpoi(l)) ) ptr_col(c) = set_nourb
        end do
      end if
      if ( present(set_spec) ) then
        do c = begc , endc
          l = clm3%g%l%c%landunit(c)
          if ( clm3%g%l%ifspecial(l) ) ptr_col(c) = set_spec
        end do
      end if
    else if ( present(ptr_pft) ) then
      l_type1d = namep(1:8)
      l_type1d_out = namep(1:8)
      clmptr_rs(hpindex)%ptr => ptr_pft
      if ( present(set_lake) ) then
        do p = begp , endp
          l = clm3%g%l%c%p%landunit(p)
          if ( clm3%g%l%lakpoi(l) ) ptr_pft(p) = set_lake
        end do
      end if
      if ( present(set_nolake) ) then
        do p = begp , endp
          l = clm3%g%l%c%p%landunit(p)
          if ( .not. (clm3%g%l%lakpoi(l)) ) ptr_pft(p) = set_nolake
        end do
      end if
      if ( present(set_urb) ) then
        do p = begp , endp
          l = clm3%g%l%c%p%landunit(p)
          if ( clm3%g%l%urbpoi(l) ) ptr_pft(p) = set_urb
        end do
      end if
      if ( present(set_nourb) ) then
        do p = begp , endp
          l = clm3%g%l%c%p%landunit(p)
          if ( .not. (clm3%g%l%urbpoi(l)) ) ptr_pft(p) = set_nourb
        end do
      end if
      if (present(set_spec)) then
        do p = begp , endp
          l = clm3%g%l%c%p%landunit(p)
          if ( clm3%g%l%ifspecial(l) ) ptr_pft(p) = set_spec
        end do
      end if
    else
      write(stderr,*) trim(subname), &
          ' ERROR: must specify a valid pointer index, choices are ', &
          '[ptr_atm, ptr_lnd, ptr_gcell, ptr_lunit, ptr_col, ptr_pft] '
      call fatal(__FILE__,__LINE__,'clm now stopping.')
    end if

    ! Set scaling factor

    scale_type_p2c = 'unity'
    scale_type_c2l = 'unity'
    scale_type_l2g = 'unity'

    if ( present(p2c_scale_type) ) scale_type_p2c = p2c_scale_type
    if ( present(c2l_scale_type) ) scale_type_c2l = c2l_scale_type
    if ( present(l2g_scale_type) ) scale_type_l2g = l2g_scale_type
    if ( present(type1d_out) ) l_type1d_out = type1d_out

    ! Add field to masterlist

    call masterlist_addfld(fname=trim(fname), &
      type1d=l_type1d, type1d_out=l_type1d_out, type2d='unset', num2d=1,  &
      units=units, avgflag=avgflag, long_name=long_name, hpindex=hpindex, &
      p2c_scale_type=scale_type_p2c, c2l_scale_type=scale_type_c2l,       &
      l2g_scale_type=scale_type_l2g)

    if ( present(default) ) then
      if ( trim(default) == 'inactive' ) return
    else
      call masterlist_make_active(name=trim(fname),tape_index=1)
    end if
  end subroutine hist_addfld1d
  !
  ! Initialize a single level history field. The pointer, ptrhist,
  ! is a pointer to the clmtype array that the history buffer will use.
  ! The value of type1d passed to masterlist\_add\_fld determines which of the
  ! 1d type of the output and the beginning and ending indices the history
  ! buffer field). Default history contents for given field on all tapes
  ! are set by calling [masterlist\_make\_active] for the appropriatae tape.
  ! After the masterlist is built, routine [htapes\_build] is called for an
  ! initial run to initialize the actual history tapes.
  !
  subroutine hist_addfld2d(fname,type2d,units,avgflag,long_name,type1d_out, &
                        ptr_gcell,ptr_lunit,ptr_col,ptr_pft,ptr_lnd,ptr_atm, &
                        p2c_scale_type,c2l_scale_type,l2g_scale_type, &
                        set_lake,set_nolake,set_urb,set_nourb,set_spec,default)
    implicit none
    character(len=*) , intent(in) :: fname   ! field name
    character(len=*) , intent(in) :: type2d  ! 2d output type
    character(len=*) , intent(in) :: units   ! units of field
    ! time averaging flag
    character(len=1) , intent(in) :: avgflag
    ! long name of field
    character(len=*) , intent(in) :: long_name
    ! output type (from clmtype)
    character(len=*) , optional , intent(in) :: type1d_out
    real(rk8) , optional , pointer :: ptr_atm(:,:)   ! pointer to atm array
    real(rk8) , optional , pointer :: ptr_lnd(:,:)   ! pointer to lnd array
    real(rk8) , optional , pointer :: ptr_gcell(:,:) ! pointer to gridcell array
    real(rk8) , optional , pointer :: ptr_lunit(:,:) ! pointer to landunit array
    real(rk8) , optional , pointer :: ptr_col(:,:)   ! pointer to column array
    real(rk8) , optional , pointer :: ptr_pft(:,:)   ! pointer to pft array
    real(rk8) , optional , intent(in) :: set_lake    ! value to set lakes to
    real(rk8) , optional , intent(in) :: set_nolake  ! value to set non-lakes to
    real(rk8) , optional , intent(in) :: set_urb     ! value to set urban to
    real(rk8) , optional , intent(in) :: set_nourb   ! value to set non-urban to
    real(rk8) , optional , intent(in) :: set_spec    ! value to set special to
    ! scale type for subgrid averaging of pfts to column
    character(len=*) , optional , intent(in) :: p2c_scale_type
    ! scale type for subgrid averaging of columns to landunits
    character(len=*) , optional , intent(in) :: c2l_scale_type
    ! scale type for subgrid averaging of landunits to gridcells
    character(len=*) , optional , intent(in) :: l2g_scale_type
    ! if set to 'inactive, field will not appear on primary tape
    character(len=*) , optional , intent(in) :: default
    integer(ik4) :: p , c , l , g ! indices
    ! size of second dimension (e.g. number of vertical levels)
    integer(ik4) :: num2d
    integer(ik4) :: begp , endp ! per-proc beginning and ending pft indices
    integer(ik4) :: begc , endc ! per-proc beginning and ending column indices
    integer(ik4) :: begl , endl ! per-proc beginning and ending landunit indices
    integer(ik4) :: begg , endg ! per-proc gridcell ending gridcell indices
    integer(ik4) :: hpindex     ! history buffer index
    character(len=8) :: l_type1d      ! 1d data type
    character(len=8) :: l_type1d_out  ! 1d output type
    ! scale type for subgrid averaging of pfts to column
    character(len=8) :: scale_type_p2c
    ! scale type for subgrid averaging of columns to landunits
    character(len=8) :: scale_type_c2l
    ! scale type for subgrid averaging of landunits to gridcells
    character(len=8) :: scale_type_l2g
    character(len=*),parameter :: subname = 'hist_addfld2d'

    ! Determine second dimension size

    select case (type2d)
      case ('levgrnd')
        num2d = nlevgrnd
      case ('levlak')
        num2d = nlevlak
      case ('numrad')
        num2d = numrad
      case ('levdcmp')
        num2d = nlevdecomp_full
      case default
        write(stderr,*) trim(subname),' ERROR: unsupported 2d type ',type2d, &
          ' currently supported types for multi level fields are ', &
          '[levgrnd,levlak,numrad,levdcmp]'
        call fatal(__FILE__,__LINE__,'clm now stopping.')
    end select

    ! Determine processor bounds

    call get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)

    ! History buffer pointer

    hpindex = pointer_index()

    if ( present(ptr_lnd) ) then
      l_type1d = grlnd(1:8)
      l_type1d_out = grlnd(1:8)
      clmptr_ra(hpindex)%ptr => ptr_lnd
    else if ( present(ptr_gcell) ) then
      l_type1d = nameg(1:8)
      l_type1d_out = nameg(1:8)
      clmptr_ra(hpindex)%ptr => ptr_gcell
    else if ( present(ptr_lunit) ) then
      l_type1d = namel(1:8)
      l_type1d_out = namel(1:8)
      clmptr_ra(hpindex)%ptr => ptr_lunit
      if ( present(set_lake) ) then
        do l = begl , endl
          if ( clm3%g%l%lakpoi(l) ) ptr_lunit(l,:) = set_lake
        end do
      end if
      if ( present(set_nolake) ) then
        do l = begl , endl
          if ( .not. (clm3%g%l%lakpoi(l)) ) ptr_lunit(l,:) = set_nolake
        end do
      end if
      if ( present(set_urb) ) then
        do l = begl , endl
          if ( clm3%g%l%urbpoi(l) ) ptr_lunit(l,:) = set_urb
        end do
      end if
      if ( present(set_nourb) ) then
        do l = begl , endl
          if ( .not. (clm3%g%l%urbpoi(l)) ) ptr_lunit(l,:) = set_nourb
        end do
      end if
      if ( present(set_spec) ) then
        do l = begl , endl
          if ( clm3%g%l%ifspecial(l) ) ptr_lunit(l,:) = set_spec
        end do
      end if
    else if (present(ptr_col)) then
      l_type1d = namec(1:8)
      l_type1d_out = namec(1:8)
      clmptr_ra(hpindex)%ptr => ptr_col
      if ( present(set_lake) ) then
        do c = begc , endc
          l = clm3%g%l%c%landunit(c)
          if ( clm3%g%l%lakpoi(l) ) ptr_col(c,:) = set_lake
        end do
      end if
      if ( present(set_nolake) ) then
        do c = begc , endc
          l = clm3%g%l%c%landunit(c)
          if ( .not. (clm3%g%l%lakpoi(l)) ) ptr_col(c,:) = set_nolake
        end do
      end if
      if ( present(set_urb) ) then
        do c = begc , endc
          l = clm3%g%l%c%landunit(c)
          if ( clm3%g%l%urbpoi(l) ) ptr_col(c,:) = set_urb
        end do
      end if
      if ( present(set_nourb) ) then
        do c = begc , endc
          l = clm3%g%l%c%landunit(c)
          if ( .not. (clm3%g%l%urbpoi(l)) ) ptr_col(c,:) = set_nourb
        end do
      end if
      if ( present(set_spec) ) then
        do c = begc , endc
          l = clm3%g%l%c%landunit(c)
          if ( clm3%g%l%ifspecial(l) ) ptr_col(c,:) = set_spec
        end do
      end if
    else if ( present(ptr_pft) ) then
      l_type1d = namep(1:8)
      l_type1d_out = namep(1:8)
      clmptr_ra(hpindex)%ptr => ptr_pft
      if ( present(set_lake) ) then
        do p = begp , endp
          l = clm3%g%l%c%p%landunit(p)
          if ( clm3%g%l%lakpoi(l) ) ptr_pft(p,:) = set_lake
        end do
      end if
      if ( present(set_nolake) ) then
        do p = begp , endp
          l = clm3%g%l%c%p%landunit(p)
          if ( .not. (clm3%g%l%lakpoi(l)) ) ptr_pft(p,:) = set_nolake
        end do
      end if
      if ( present(set_urb) ) then
        do p = begp , endp
          l = clm3%g%l%c%p%landunit(p)
          if ( clm3%g%l%urbpoi(l) ) ptr_pft(p,:) = set_urb
        end do
      end if
      if ( present(set_nourb) ) then
        do p = begp , endp
          l = clm3%g%l%c%p%landunit(p)
          if ( .not. (clm3%g%l%urbpoi(l)) ) ptr_pft(p,:) = set_nourb
        end do
      end if
      if ( present(set_spec) ) then
        do p = begp , endp
          l = clm3%g%l%c%p%landunit(p)
          if ( clm3%g%l%ifspecial(l) ) ptr_pft(p,:) = set_spec
        end do
      end if
    else
      write(stderr,*) trim(subname), &
         ' ERROR: must specify a valid pointer index,', &
         ' choices are ptr_atm, ptr_lnd, ptr_gcell, ptr_lunit, ptr_col, ptr_pft'
      call fatal(__FILE__,__LINE__,'clm now stopping.')
    end if

    ! Set scaling factor

    scale_type_p2c = 'unity'
    scale_type_c2l = 'unity'
    scale_type_l2g = 'unity'

    if ( present(p2c_scale_type) ) scale_type_p2c = p2c_scale_type
    if ( present(c2l_scale_type) ) scale_type_c2l = c2l_scale_type
    if ( present(l2g_scale_type) ) scale_type_l2g = l2g_scale_type
    if ( present(type1d_out) ) l_type1d_out = type1d_out

    ! Add field to masterlist

    call masterlist_addfld(fname=trim(fname), &
         type1d=l_type1d,type1d_out=l_type1d_out,type2d=type2d, num2d=num2d, &
         units=units, avgflag=avgflag, long_name=long_name, hpindex=hpindex, &
         p2c_scale_type=scale_type_p2c, c2l_scale_type=scale_type_c2l,       &
         l2g_scale_type=scale_type_l2g)

    if (present(default)) then
      if ( trim(default) == 'inactive' ) return
    else
      call masterlist_make_active(name=trim(fname), tape_index=1)
    end if
  end subroutine hist_addfld2d
  !
  ! Get the maximum number of fields on all tapes.
  !
  integer(ik4) function max_nFields()
    implicit none
    integer(ik4) :: t  ! index
    max_nFields = 0
    do t = 1 , ntapes
      max_nFields = max(max_nFields,tape(t)%nflds)
    end do
  end function max_nFields
  !
  ! Retrieve name portion of inname. If an averaging flag separater character
  ! is present (:) in inname, lop it off.
  !
   character(len=max_namlen) function getname (inname)
     implicit none
     character(len=*) , intent(in) :: inname
     integer(ik4) :: length
     integer(ik4) :: i
     character(len=*) , parameter :: subname = 'getname'

     length = len(inname)

     if ( length < max_namlen .or. length > max_namlen+2 ) then
       write(stderr,*) trim(subname),' ERROR: bad length=',length
       call fatal(__FILE__,__LINE__,'clm now stopping.')
     end if

     getname = ' '
     do i = 1 , max_namlen
       if (inname(i:i) == ':') exit
       getname(i:i) = inname(i:i)
     end do
   end function getname
   !
   ! Retrieve flag portion of inname. If an averaging flag separater character
   ! is present (:) in inname, return the character after it as the flag
   !
   character(len=1) function getflag (inname)
     implicit none
     character(len=*) :: inname ! character string
     integer(ik4) :: length     ! length of inname
     integer(ik4) :: i          ! loop index
     character(len=*) , parameter :: subname = 'getflag'
     length = len(inname)
     if ( length < max_namlen .or. length > max_namlen+2 ) then
       write(stderr,*) trim(subname),' ERROR: bad length=',length
       call fatal(__FILE__,__LINE__,'clm now stopping.')
     end if
     getflag = ' '
     do i = 1 , length
       if ( inname(i:i) == ':' ) then
         getflag = inname(i+1:i+1)
         exit
       end if
     end do
  end function getflag

  subroutine list_index (list, sname, iindex)
    implicit none
    ! input list of names, possibly ":" delimited
    character(len=*) , intent(in) :: list(max_flds)
    ! name to be searched for
    character(len=max_namlen) , intent(in) :: sname
    ! index of "name" in "list"
    integer(ik4) , intent(out) :: iindex
    ! input name with ":" stripped off.
    character(len=max_namlen) :: listname
    integer(ik4) :: f  ! field index
    character(len=*) , parameter :: subname = 'list_index'
    ! Only list items
    iindex = 0
    do f = 1 , max_flds
      listname = getname (list(f))
      if (listname == ' ') exit
      if (listname == sname) then
        iindex = f
        exit
      end if
    end do
  end subroutine list_index
  !
  ! Determine history dataset filenames.
  !
  character(len=256) function set_hist_filename(hist_freq,hist_mfilt,hist_file)
    implicit none
    integer(ik4) , intent(in) :: hist_freq  !history file frequency
    integer(ik4) , intent(in) :: hist_mfilt !history file number of time-samples
    integer(ik4) , intent(in) :: hist_file  !history file index
    character(len=256) :: cdate       !date char string
    character(len=  1) :: hist_index  !p,1 or 2 (currently)
    integer(ik4) :: day !day (1 -> 31)
    integer(ik4) :: mon !month (1 -> 12)
    integer(ik4) :: yr  !year (0 -> ...)
    integer(ik4) :: sec !seconds into current day
    character(len=*) , parameter :: subname = 'set_hist_filename'

    if (hist_freq == 0 .and. hist_mfilt == 1) then   !monthly
      call get_prev_date (yr, mon, day, sec)
      write(cdate,'(i4.4,"-",i2.2)') yr,mon
    else                        !other
      call curr_date (idatex, yr, mon, day, sec)
      write(cdate,'(i4.4,"-",i2.2,"-",i2.2,"-",i5.5)') yr,mon,day,sec
    end if
    write(hist_index,'(i1.1)') hist_file - 1
    set_hist_filename = trim(dirout)//trim(caseid)//".clm."// &
                        trim(inst_suffix)//".h"//hist_index// &
                        "."//trim(cdate)//".nc"
  end function set_hist_filename
  !
  ! Set the current pointer index and increment the value of the index.
  !
  integer(ik4) function pointer_index ()
    implicit none
    integer(ik4) , save :: lastindex = 1
    character(len=*) , parameter :: subname = 'pointer_index'
    pointer_index = lastindex
    lastindex = lastindex + 1
    if (lastindex > max_mapflds) then
      write(stderr,*) trim(subname),' ERROR: ',&
           ' lastindex = ',lastindex,' greater than max_mapflds= ',max_mapflds
      call fatal(__FILE__,__LINE__,'clm now stopping.')
    end if
  end function pointer_index
  !
  ! Add a history variable to the output history tape.
  !
  subroutine hist_add_subscript(name, dim)
    implicit none
    character(len=*) , intent(in) :: name ! name of subscript
    integer(ik4) , intent(in) :: dim      ! dimension of subscript
    character(len=*) , parameter :: subname = 'hist_add_subscript'
    num_subs = num_subs + 1
    if ( num_subs > max_subs ) then
      write(stderr,*) trim(subname),' ERROR: ',&
            ' num_subs = ',num_subs,' greater than max_subs= ',max_subs
      call fatal(__FILE__,__LINE__,'clm now stopping.')
    end if
    subs_name(num_subs) = name
    subs_dim(num_subs) =  dim
  end subroutine hist_add_subscript

  subroutine strip_null(str)
    implicit none
    character(len=*) , intent(inout) :: str
    integer(ik4) :: i
    do i = 1 , len(str)
      if ( ichar(str(i:i)) == 0 ) str(i:i) = ' '
    end do
  end subroutine strip_null
  !
  ! Determine logic for closeing and/or disposing history file
  ! Sets values for if_disphist, if_stop (arguments)
  ! Remove history files unless this is end of run or
  ! history file is not full.
  !
  subroutine hist_do_disp(ntapes,hist_ntimes,hist_mfilt,if_stop, &
                  if_disphist,rstwr,nlend)
    implicit none
    !actual number of history tapes
    integer(ik4) , intent(in) :: ntapes
    !current numbers of time samples on history tape
    integer(ik4) , intent(in) :: hist_ntimes(ntapes)
    !maximum number of time samples per tape
    integer(ik4) , intent(in) :: hist_mfilt(ntapes)
    !true => last time step of run
    logical , intent(out) :: if_stop
    !true => save and dispose history file
    logical , intent(out) :: if_disphist(ntapes)
    logical , intent(in) :: rstwr
    logical , intent(in) :: nlend
    integer(ik4) :: t         ! history tape index
    logical :: rest_now       ! temporary
    logical :: stop_now       ! temporary
    rest_now = .false.
    stop_now = .false.
    if ( nlend ) stop_now = .true.
    if ( rstwr ) rest_now = .true.
    if_stop = stop_now
    if ( stop_now ) then
      ! End of run -  dispose all history files
      if_disphist(1:ntapes) = .true.
    else if ( rest_now ) then
      ! Restart - dispose all history files
      do t = 1 , ntapes
        if_disphist(t) = .true.
      end do
    else
      ! Dispose
      if_disphist(1:ntapes) = .false.
      do t = 1 , ntapes
        if ( hist_ntimes(t) ==  hist_mfilt(t) ) then
          if_disphist(t) = .true.
        end if
      end do
    end if
  end subroutine hist_do_disp

end module mod_clm_histfile

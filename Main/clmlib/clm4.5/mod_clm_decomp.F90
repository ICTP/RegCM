module mod_clm_decomp
  !
  ! Module provides a descomposition into a data structure which can
  ! be mapped back to atmosphere physics chunks.
  ! This has been simplified, as in regional model the communication
  ! ATM <-> LND is managed in the RegCM code.
  !
  use mod_stdio
  use mod_realkinds
  use mod_intkinds
  use mod_mpmessage
  use mod_regcm_types
  use mod_clm_type , only : grlnd , nameg , namel , namec , namep
  use mod_clm_domain , only : ldomain

  implicit none

  private

  ! this processor beg and end gridcell , landunit , column,pft
  public :: get_proc_bounds

  ! total number of gridcells , landunits , columns and pfts for any processor
  public :: get_proc_total

  ! total gridcells , landunits , columns , pfts across all processors
  public :: get_proc_global

  ! get global size associated with clmlevel
  public :: get_clmlevel_gsize

  ! total number of gridcells on all procs
  integer(ik4) , public :: numg
  ! total number of landunits on all procs
  integer(ik4) , public :: numl
  ! total number of columns on all procs
  integer(ik4) , public :: numc
  ! total number of pfts on all procs
  integer(ik4) , public :: nump

  !---global information on each pe
  type processor_type
    type (masked_comm) , pointer :: cl
    logical , dimension(:,:) , pointer :: gcmask
    integer(ik4) :: icomm
    integer(ik4) :: ncells           ! number of gridcells in proc
    integer(ik4) :: nlunits          ! number of landunits in proc
    integer(ik4) :: ncols            ! number of columns in proc
    integer(ik4) :: npfts            ! number of pfts in proc
    integer(ik4) :: begg , endg      ! beginning and ending gridcell index
    integer(ik4) :: begl , endl      ! beginning and ending landunit index
    integer(ik4) :: begc , endc      ! beginning and ending column index
    integer(ik4) :: begp , endp      ! beginning and ending pft index
    integer(ik4) , pointer , dimension(:) :: gc
    integer(ik4) , pointer , dimension(:) :: gd
    integer(ik4) , pointer , dimension(:) :: lc
    integer(ik4) , pointer , dimension(:) :: ld
    integer(ik4) , pointer , dimension(:) :: cc
    integer(ik4) , pointer , dimension(:) :: cd
    integer(ik4) , pointer , dimension(:) :: pc
    integer(ik4) , pointer , dimension(:) :: pd
  end type processor_type

  type subgrid_type
    integer(ik4) :: icomm
    integer(ik4) :: ns , is , ie
    integer(ik4) , pointer , dimension(:) :: ic
    integer(ik4) , pointer , dimension(:) :: id
  end type subgrid_type

  public processor_type , subgrid_type

  type(processor_type) , public :: procinfo

  type(subgrid_type) , public , target :: gcomm_gridcell
  type(subgrid_type) , public , target :: gcomm_landunit
  type(subgrid_type) , public , target :: gcomm_column
  type(subgrid_type) , public , target :: gcomm_pft

  contains
    !
    ! Retrieve gridcell, landunit, column, and pft bounds for process.
    !
    subroutine get_proc_bounds(begg,endg,begl,endl,begc,endc,begp,endp)
      implicit none
      ! proc beg and end pft indices
      integer(ik4) , optional , intent(out) :: begp , endp
      ! proc beg and end column indices
      integer(ik4) , optional , intent(out) :: begc , endc
      ! proc beg and end landunit indices
      integer(ik4) , optional , intent(out) :: begl , endl
      ! proc beg and end gridcell indices
      integer(ik4) , optional , intent(out) :: begg , endg
      if ( present(begp) ) begp = procinfo%begp
      if ( present(endp) ) endp = procinfo%endp
      if ( present(begc) ) begc = procinfo%begc
      if ( present(endc) ) endc = procinfo%endc
      if ( present(begl) ) begl = procinfo%begl
      if ( present(endl) ) endl = procinfo%endl
      if ( present(begg) ) begg = procinfo%begg
      if ( present(endg) ) endg = procinfo%endg
    end subroutine get_proc_bounds
    !
    ! Count up gridcells, landunits, columns, and pfts on process.
    !
    subroutine get_proc_total(pid, ncells, nlunits, ncols, npfts)
      implicit none
      integer(ik4) , intent(in)  :: pid     ! proc id
      ! total number of gridcells on the processor
      integer(ik4) , intent(out) :: ncells
      ! total number of landunits on the processor
      integer(ik4) , intent(out) :: nlunits
      ! total number of columns on the processor
      integer(ik4) , intent(out) :: ncols
      ! total number of pfts on the processor
      integer(ik4) , intent(out) :: npfts
      ncells   = procinfo%ncells
      nlunits  = procinfo%nlunits
      ncols    = procinfo%ncols
      npfts    = procinfo%npfts
    end subroutine get_proc_total
    !
    ! Return number of gridcells, landunits, columns, and pfts across all
    ! processes.
    !
    subroutine get_proc_global(ng,nl,nc,np)
      implicit none
      ! total number of gridcells across all processors
      integer(ik4) , intent(out) :: ng
      ! total number of landunits across all processors
      integer(ik4) , intent(out) :: nl
      ! total number of columns across all processors
      integer(ik4) , intent(out) :: nc
      ! total number of pfts across all processors
      integer(ik4) , intent(out) :: np
      np = nump
      nc = numc
      nl = numl
      ng = numg
    end subroutine get_proc_global
    !
    ! Determine 1d size from clmlevel
    !
    integer(ik4) function get_clmlevel_gsize(clmlevel)
      implicit none
      character(len=*) , intent(in) :: clmlevel ! type of clm 1d array
      select case (clmlevel)
        case(grlnd)
          get_clmlevel_gsize = ldomain%ns
        case(nameg)
          get_clmlevel_gsize = numg
        case(namel)
          get_clmlevel_gsize = numl
        case(namec)
          get_clmlevel_gsize = numc
        case(namep)
          get_clmlevel_gsize = nump
        case default
          write(stderr,*) &
            'get_clmlevel_gsize does not match clmlevel type: ', trim(clmlevel)
          call fatal(__FILE__,__LINE__,'clm will now stop')
      end select
    end function get_clmlevel_gsize

end module mod_clm_decomp

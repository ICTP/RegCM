module mod_clm_decomp
  !
  ! Module provides a descomposition into a clumped data structure which can
  ! be mapped back to atmosphere physics chunks.
  !
  use mod_stdio
  use mod_realkinds
  use mod_intkinds
  use mod_clm_type , only : grlnd , nameg , namel , namec , namep
  use mod_clm_domain , only : ldomain

  implicit none

  private

  integer(ik4) , public :: clump_pproc ! number of clumps per MPI process

  ! clump beg and end gridcell,landunit,column,pft
  public get_clump_bounds
  ! number of clumps for this processor
  public get_proc_clumps
  ! this processor beg and end gridcell,landunit,column,pft
  public get_proc_bounds
  ! total number of gridcells, landunits, columns and pfts for any processor
  public get_proc_total
  ! total gridcells, landunits, columns, pfts across all processors
  public get_proc_global
  ! get global size associated with clmlevel
  public get_clmlevel_gsize

  ! total number of clumps across all processors
  integer(ik4) , public :: nclumps
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
    integer(ik4)  :: nclumps          ! number of clumps for processor_type iam
    integer(ik4) , pointer , dimension(:) :: cid  ! clump indices
    integer(ik4)  :: ncells           ! number of gridcells in proc
    integer(ik4)  :: nlunits          ! number of landunits in proc
    integer(ik4)  :: ncols            ! number of columns in proc
    integer(ik4)  :: npfts            ! number of pfts in proc
    integer(ik4)  :: begg, endg       ! beginning and ending gridcell index
    integer(ik4)  :: begl, endl       ! beginning and ending landunit index
    integer(ik4)  :: begc, endc       ! beginning and ending column index
    integer(ik4)  :: begp, endp       ! beginning and ending pft index
  end type processor_type

  public processor_type
  type(processor_type) , public :: procinfo

  !---global information on each pe
  type clump_type
    integer(ik4)  :: owner            ! process id owning clump
    integer(ik4)  :: ncells           ! number of gridcells in clump
    integer(ik4)  :: nlunits          ! number of landunits in clump
    integer(ik4)  :: ncols            ! number of columns in clump
    integer(ik4)  :: npfts            ! number of pfts in clump
    integer(ik4)  :: begg , endg      ! beginning and ending gridcell index
    integer(ik4)  :: begl , endl      ! beginning and ending landunit index
    integer(ik4)  :: begc , endc      ! beginning and ending column index
    integer(ik4)  :: begp , endp      ! beginning and ending pft index
  end type clump_type

  public clump_type
  type(clump_type) , public , allocatable , dimension(:) :: clumps

  !---global information on each pe
  !--- i,j = 2d global
  !--- glo = 1d global sn ordered
  !--- gsn = 1d global sn ordered compressed
  !--- gdc = 1d global dc ordered compressed
  type decomp_type
     integer(ik4) , pointer , dimension(:) :: glo2gdc    ! 1d glo to 1d gdc
     integer(ik4) , pointer , dimension(:) :: gdc2glo    ! 1d gdc to 1d glo
  end type decomp_type

  public decomp_type
  type(decomp_type) , public , target :: ldecomp

!  type(mct_gsMap)  ,public,target :: gsMap_lnd_gdc2glo
!  type(mct_gsMap)  ,public,target :: gsMap_gce_gdc2glo
!  type(mct_gsMap)  ,public,target :: gsMap_lun_gdc2glo
!  type(mct_gsMap)  ,public,target :: gsMap_col_gdc2glo
!  type(mct_gsMap)  ,public,target :: gsMap_pft_gdc2glo

  contains
    !
    ! Determine clump beginning and ending pft, column, landunit and
    ! gridcell indices.
    !
    subroutine get_clump_bounds(n,begg,endg,begl,endl,begc,endc,begp,endp)
      implicit none
      integer(ik4) , intent(in)  :: n           ! proc clump index
      ! clump beg and end pft indices
      integer(ik4) , intent(out) :: begp , endp
      ! clump beg and end column indices
      integer(ik4) , intent(out) :: begc , endc
      ! clump beg and end landunit indices
      integer(ik4) , intent(out) :: begl , endl
      ! clump beg and end gridcell indices
      integer(ik4) , intent(out) :: begg , endg
      integer(ik4)  :: cid ! clump id
      cid  = procinfo%cid(n)
      begp = clumps(cid)%begp
      endp = clumps(cid)%endp
      begc = clumps(cid)%begc
      endc = clumps(cid)%endc
      begl = clumps(cid)%begl
      endl = clumps(cid)%endl
      begg = clumps(cid)%begg
      endg = clumps(cid)%endg
    end subroutine get_clump_bounds
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

      integer(ik4)  :: cid       ! clump index

      npfts   = 0
      nlunits = 0
      ncols   = 0
      ncells  = 0
      do cid = 1 , nclumps
        if ( clumps(cid)%owner == pid ) then
          ncells  = ncells  + clumps(cid)%ncells
          nlunits = nlunits + clumps(cid)%nlunits
          ncols   = ncols   + clumps(cid)%ncols
          npfts   = npfts   + clumps(cid)%npfts
        end if
      end do
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
    ! Return the number of clumps.
    !
    integer(ik4) function get_proc_clumps()
      implicit none
      get_proc_clumps = procinfo%nclumps
    end function get_proc_clumps
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

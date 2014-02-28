module mod_clm_decompinit
  !
  ! Module provides a descomposition into a data structure which can
  ! be mapped back to atmosphere physics chunks.
  !
  use mod_realkinds
  use mod_intkinds
  use mod_memutil
  use mod_mpmessage
  use mod_stdio
  use mod_dynparam
  use mod_mppparam
  use mod_regcm_types
  use mod_service
  use mod_clm_decomp
  use mod_clm_subgrid , only : subgrid_get_gcellinfo
  use mpi

  implicit none

  private

  ! initializes atm grid decomposition into processors
  public :: decompInit_lnd
  ! initializes g , l , c , p decomp info
  public :: decompInit_glcp

  contains
  !
  ! This subroutine initializes the land surface decomposition into a
  ! processor_type data structure.
  !
  subroutine decompInit_lnd(cl)
    implicit none
    type (masked_comm) , intent(in) , target :: cl
    integer(ik4) , pointer , dimension(:) :: gcount
    integer(ik4) :: nout , mynumg , np
    logical , dimension(:,:,:) , pointer  :: lsubgrid

    allocate(procinfo%gc(nproc))
    allocate(procinfo%gd(nproc))
    allocate(procinfo%lc(nproc))
    allocate(procinfo%ld(nproc))
    allocate(procinfo%cc(nproc))
    allocate(procinfo%cd(nproc))
    allocate(procinfo%pc(nproc))
    allocate(procinfo%pd(nproc))

    procinfo%cl => cl
    procinfo%ncells = cl%linear_npoint_sg(myid+1)

    procinfo%icomm = cl%linear_communicator
    gcomm_gridcell%icomm = cl%linear_communicator
    gcomm_landunit%icomm = cl%linear_communicator
    gcomm_column%icomm = cl%linear_communicator
    gcomm_pft%icomm = cl%linear_communicator

    allocate(gcount(procinfo%ncells))
    gcount(:) = 1
    mynumg  = sum(gcount)
    call sumall(mynumg,numg)
    call allgather_i(procinfo%gc,mynumg)

    if ( myid > 0 ) then
      procinfo%begg = sum(procinfo%gc(1:myid))+1
      procinfo%endg = procinfo%begg + mynumg - 1
    else
      procinfo%begg = 1
      procinfo%endg = mynumg
    end if

    procinfo%gd(:) = 0
    do np = 2 , nproc
      procinfo%gd(np) = sum(procinfo%gc(1:np-1))
    end do

    deallocate(gcount)

    gcomm_gridcell%ns = numg
    gcomm_gridcell%is = procinfo%begg
    gcomm_gridcell%ie = procinfo%endg
    gcomm_gridcell%ic => procinfo%gc
    gcomm_gridcell%id => procinfo%gd

    ! Set default for landunits , column and pfts

    procinfo%nlunits = 0
    procinfo%ncols   = 0
    procinfo%npfts   = 0

    procinfo%begl    = 1
    procinfo%begc    = 1
    procinfo%begp    = 1

    procinfo%endl    = 0
    procinfo%endc    = 0
    procinfo%endp    = 0
  end subroutine decompInit_lnd
  !
  ! This subroutine initializes the land surface decomposition into a
  ! data structure. 
  !
  subroutine decompInit_glcp
    implicit none
    integer(ik4) :: begg , endg  ! beg , end gridcells
    integer(ik4) :: ln           ! lnd num gridcells
    integer(ik4) :: mynumg , mynumc , mynump , mynuml
    integer(ik4) :: np , ilunits , icols , ipfts
    integer(ik4) , pointer , dimension(:) :: lcount
    integer(ik4) , pointer , dimension(:) :: ccount
    integer(ik4) , pointer , dimension(:) :: pcount

    allocate(lcount(procinfo%begg:procinfo%endg))
    allocate(ccount(procinfo%begg:procinfo%endg))
    allocate(pcount(procinfo%begg:procinfo%endg))

    lcount = 0
    ccount = 0
    pcount = 0

    begg = procinfo%begg
    endg = procinfo%endg
    write(0,*) 'MYID = ',myid, 'begg = ', begg
    write(0,*) 'MYID = ',myid, 'endg = ', endg
    do ln = begg , endg
      call subgrid_get_gcellinfo(ln,nlunits=ilunits,ncols=icols,npfts=ipfts)
      lcount(ln) = ilunits
      ccount(ln) = icols
      pcount(ln) = ipfts
      write(0,*) 'MYID = ',myid, ln , ilunits , icols , ipfts
    end do

    mynuml  = sum(lcount)
    mynumc  = sum(ccount)
    mynump  = sum(pcount)

    procinfo%nlunits = mynuml
    procinfo%ncols   = mynumc
    procinfo%npfts   = mynump

    call sumall(mynuml,numl)
    call sumall(mynumc,numc)
    call sumall(mynump,nump)

    call allgather_i(procinfo%lc,mynuml)
    call allgather_i(procinfo%cc,mynumc)
    call allgather_i(procinfo%pc,mynump)

    write(0,*) 'MYID = ',myid, 'mynuml = ', mynuml
    write(0,*) 'MYID = ',myid, 'mynumc = ', mynumc
    write(0,*) 'MYID = ',myid, 'mynump = ', mynump

    call fatal(__FILE__,__LINE__,'clm now stopping')

    if ( myid > 0 ) then
      procinfo%begl = sum(procinfo%lc(1:myid))+1
      procinfo%endl = procinfo%begl + mynuml - 1
      procinfo%begc = sum(procinfo%cc(1:myid))+1
      procinfo%endc = procinfo%begc + mynumc - 1
      procinfo%begp = sum(procinfo%pc(1:myid))+1
      procinfo%endp = procinfo%begp + mynump - 1
    else
      procinfo%begl = 1
      procinfo%endl = mynuml
      procinfo%begc = 1
      procinfo%endc = mynumc
      procinfo%begp = 1
      procinfo%endp = mynump
    end if

    procinfo%ld(:) = 0
    procinfo%cd(:) = 0
    procinfo%pd(:) = 0
    do np = 2 , nproc
      procinfo%ld(np) = sum(procinfo%lc(1:np-1))
      procinfo%cd(np) = sum(procinfo%cc(1:np-1))
      procinfo%pd(np) = sum(procinfo%pc(1:np-1))
    end do

    deallocate(lcount)
    deallocate(ccount)
    deallocate(pcount)

    gcomm_landunit%ns = numl
    gcomm_landunit%is = procinfo%begl
    gcomm_landunit%ie = procinfo%endl
    gcomm_landunit%ic => procinfo%lc
    gcomm_landunit%id => procinfo%ld

    gcomm_column%ns = numc
    gcomm_column%is = procinfo%begc
    gcomm_column%ie = procinfo%endc
    gcomm_column%ic => procinfo%cc
    gcomm_column%id => procinfo%cd

    gcomm_pft%ns = nump
    gcomm_pft%is = procinfo%begp
    gcomm_pft%ie = procinfo%endp
    gcomm_pft%ic => procinfo%pc
    gcomm_pft%id => procinfo%pd

    ! Diagnostic output

    if ( myid == italk ) then
      write(stdout,*)' Surface Grid Characteristics'
      write(stdout,*)'   total number of gridcells = ',numg
      write(stdout,*)'   total number of landunits = ',numl
      write(stdout,*)'   total number of columns   = ',numc
      write(stdout,*)'   total number of pfts      = ',nump
    end if

#ifdef DEBUG
    write(ndebug+myid,*)'proc= ',myid,&
               ' beg gridcell= ',procinfo%begg, &
               ' end gridcell= ',procinfo%endg,                   &
               ' total gridcells per proc= ',procinfo%ncells
    write(ndebug+myid,*)'proc= ',myid,&
               ' beg landunit= ',procinfo%begl, &
               ' end landunit= ',procinfo%endl,                   &
               ' total landunits per proc= ',procinfo%nlunits
    write(ndebug+myid,*)'proc= ',myid,&
               ' beg column  = ',procinfo%begc, &
               ' end column  = ',procinfo%endc,                   &
               ' total columns per proc  = ',procinfo%ncols
    write(ndebug+myid,*)'proc= ',myid,&
               ' beg pft     = ',procinfo%begp, &
               ' end pft     = ',procinfo%endp,                   &
               ' total pfts per proc     = ',procinfo%npfts
#endif

  end subroutine decompInit_glcp

end module mod_clm_decompinit

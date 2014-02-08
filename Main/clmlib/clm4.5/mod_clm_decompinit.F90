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

    allocate(procinfo%gc(nproc))
    allocate(procinfo%gd(nproc))
    allocate(procinfo%lc(nproc))
    allocate(procinfo%ld(nproc))
    allocate(procinfo%cc(nproc))
    allocate(procinfo%cd(nproc))
    allocate(procinfo%pc(nproc))
    allocate(procinfo%pd(nproc))

    call getmem2d(procinfo%gcmask,jout1,jout2,iout1,iout2,'clm decomp')
    call grid_collect(cl%gmask,procinfo%gcmask,jci1,jci2,ici1,ici2)

    procinfo%cl => cl
    procinfo%ncells = cl%linear_npoint_sg(myid+1)

    procinfo%icomm = cl%linear_communicator
    gcomm_gridcell%icomm = cl%linear_communicator
    gcomm_landunit%icomm = cl%linear_communicator
    gcomm_column%icomm = cl%linear_communicator
    gcomm_pft%icomm = cl%linear_communicator

    procinfo%nlunits = 0
    procinfo%ncols   = 0
    procinfo%npfts   = 0

    procinfo%begg    = 1
    procinfo%begl    = 1
    procinfo%begc    = 1
    procinfo%begp    = 1

    procinfo%endg    = 0
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
    integer(ik4) , pointer , dimension(:) :: gcount
    integer(ik4) , pointer , dimension(:) :: lcount
    integer(ik4) , pointer , dimension(:) :: ccount
    integer(ik4) , pointer , dimension(:) :: pcount

    allocate(gcount(procinfo%ncells))
    allocate(lcount(procinfo%ncells))
    allocate(ccount(procinfo%ncells))
    allocate(pcount(procinfo%ncells))

    gcount = 0
    lcount = 0
    ccount = 0
    pcount = 0
    begg = 1
    endg = procinfo%ncells

    do ln = begg , endg
      call subgrid_get_gcellinfo(ln,nlunits=ilunits,ncols=icols,npfts=ipfts)
      gcount(ln) = 1
      lcount(ln) = ilunits
      ccount(ln) = icols
      pcount(ln) = ipfts
    end do

    mynumg  = sum(gcount)
    mynuml  = sum(lcount)
    mynumc  = sum(ccount)
    mynump  = sum(pcount)

    procinfo%nlunits = mynuml
    procinfo%ncols   = mynumc
    procinfo%npfts   = mynump

    call sumall(mynumg,numg)
    call sumall(mynuml,numl)
    call sumall(mynumc,numc)
    call sumall(mynump,nump)

    call allgather_i(procinfo%gc,mynumg)
    call allgather_i(procinfo%lc,mynuml)
    call allgather_i(procinfo%cc,mynumc)
    call allgather_i(procinfo%pc,mynump)

    if ( myid > 1 ) then
      procinfo%begg = sum(procinfo%gc(1:myid))
      procinfo%endg = procinfo%begg + mynumg
      procinfo%begl = sum(procinfo%lc(1:myid))
      procinfo%endl = procinfo%begl + mynuml
      procinfo%begc = sum(procinfo%cc(1:myid))
      procinfo%endc = procinfo%begc + mynumc
      procinfo%begp = sum(procinfo%pc(1:myid))
      procinfo%endp = procinfo%begp + mynump
    else
      procinfo%begg = 1
      procinfo%endg = mynumg
      procinfo%begl = 1
      procinfo%endl = mynuml
      procinfo%begc = 1
      procinfo%endc = mynumc
      procinfo%begp = 1
      procinfo%endp = mynump
    end if

    procinfo%gd(:) = 0
    procinfo%ld(:) = 0
    procinfo%cd(:) = 0
    procinfo%pd(:) = 0
    do np = 2 , nproc
      procinfo%gd(np) = sum(procinfo%gc(1:np-1))
      procinfo%ld(np) = sum(procinfo%lc(1:np-1))
      procinfo%cd(np) = sum(procinfo%cc(1:np-1))
      procinfo%pd(np) = sum(procinfo%pc(1:np-1))
    end do

    deallocate(gcount)
    deallocate(lcount)
    deallocate(ccount)
    deallocate(pcount)

    gcomm_gridcell%ns = numg
    gcomm_landunit%ns = numl
    gcomm_column%ns = numc
    gcomm_pft%ns = nump
    gcomm_gridcell%is = procinfo%begg
    gcomm_landunit%is = procinfo%begl
    gcomm_column%is = procinfo%begc
    gcomm_pft%is = procinfo%begp
    gcomm_gridcell%ie = procinfo%endg
    gcomm_landunit%ie = procinfo%endl
    gcomm_column%ie = procinfo%endc
    gcomm_pft%ie = procinfo%endp
    gcomm_gridcell%ic => procinfo%gc
    gcomm_landunit%ic => procinfo%lc
    gcomm_column%ic => procinfo%cc
    gcomm_pft%ic => procinfo%pc
    gcomm_gridcell%id => procinfo%gd
    gcomm_landunit%id => procinfo%ld
    gcomm_column%id => procinfo%cd
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

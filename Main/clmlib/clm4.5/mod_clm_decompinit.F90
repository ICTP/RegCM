module decompInitMod
  !
  ! Module provides a descomposition into a data structure which can
  ! be mapped back to atmosphere physics chunks.
  !
  use mod_realkinds
  use mod_intkinds
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

  type (masked_comm) , pointer :: clm_cl

  contains
  !
  ! This subroutine initializes the land surface decomposition into a
  ! processor_type data structure.
  !
  subroutine decompInit_lnd(cl)
    implicit none
    type (masked_comm) , intent(in) , target :: cl

    clm_cl => cl
    procinfo%ncells  = cl%linear_npoint_sg(myid+1)

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
    integer(ik4) :: anumg        ! lnd num gridcells
    integer(ik4) :: ln           ! lnd num gridcells
    integer(ik4) :: mynumg , mynumc , mynump , mynuml
    integer(ik4) :: ilunits , icols , ipfts
    integer(ik4) , pointer , dimension(:) :: gcount
    integer(ik4) , pointer , dimension(:) :: lcount
    integer(ik4) , pointer , dimension(:) :: ccount
    integer(ik4) , pointer , dimension(:) :: pcount
    integer(ik4) , pointer , dimension(:,:) :: xstart , xend

    allocate(gcount(procinfo%ncells))
    allocate(lcount(procinfo%ncells))
    allocate(ccount(procinfo%ncells))
    allocate(pcount(procinfo%ncells))
    allocate(xstart(nproc,4))

    gcount = 0
    lcount = 0
    ccount = 0
    pcount = 0
    begg = 1
    endg = procinfo%ncells

    do anumg = begg , endg
      ln  = anumg
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

    call gather_i(xstart(:,1),mynumg)
    call gather_i(xstart(:,2),mynuml)
    call gather_i(xstart(:,3),mynumc)
    call gather_i(xstart(:,4),mynump)

    if ( myid > 1 ) then
      procinfo%begg = sum(xstart(1:myid,1))
      procinfo%endg = procinfo%begg + mynumg
      procinfo%begl = sum(xstart(1:myid,2))
      procinfo%endl = procinfo%begl + mynuml
      procinfo%begc = sum(xstart(1:myid,3))
      procinfo%endc = procinfo%begc + mynumc
      procinfo%begp = sum(xstart(1:myid,4))
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

    deallocate(gcount)
    deallocate(lcount)
    deallocate(ccount)
    deallocate(pcount)
    deallocate(xstart)

    ! Diagnostic output

    if ( myid == italk ) then
      write(stdout,*)' Surface Grid Characteristics'
      write(stdout,*)'   total number of gridcells = ',numg
      write(stdout,*)'   total number of landunits = ',numl
      write(stdout,*)'   total number of columns   = ',numc
      write(stdout,*)'   total number of pfts      = ',nump
    end if

#ifdef DEBUG
    write(ndebug+myid,*)'proc= ',pid,&
               ' beg gridcell= ',procinfo%begg, &
               ' end gridcell= ',procinfo%endg,                   &
               ' total gridcells per proc= ',procinfo%ncells
    write(ndebug+myid,*)'proc= ',pid,&
               ' beg landunit= ',procinfo%begl, &
               ' end landunit= ',procinfo%endl,                   &
               ' total landunits per proc= ',procinfo%nlunits
    write(ndebug+myid,*)'proc= ',pid,&
               ' beg column  = ',procinfo%begc, &
               ' end column  = ',procinfo%endc,                   &
               ' total columns per proc  = ',procinfo%ncols
    write(ndebug+myid,*)'proc= ',pid,&
               ' beg pft     = ',procinfo%begp, &
               ' end pft     = ',procinfo%endp,                   &
               ' total pfts per proc     = ',procinfo%npfts
#endif

  end subroutine decompInit_glcp

end module decompInitMod

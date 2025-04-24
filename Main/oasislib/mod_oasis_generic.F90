!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

#ifdef OASIS

module mod_oasis_generic

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_stdio
  use mod_message
  use mod_memutil
  use mod_service
  use mod_mppparam
  use mod_runparams

  use mod_oasis
  use mod_oasis_params
  use mod_oasis_signature

  implicit none

  private

  integer(ik4) :: ierror

  interface fill_ocean
    module procedure fill_ocean_2d, fill_ocean_3d
  end interface fill_ocean

  public :: oasisxregcm_init, oasisxregcm_finalize
  public :: oasisxregcm_setup_grid
  public :: oasisxregcm_setup_field, oasisxregcm_deallocate_field
  public :: oasisxregcm_def_partition
  public :: oasisxregcm_allocate_oasisgrids, srf_sqm
  public :: oasisxregcm_write_oasisgrids, oasisxregcm_deallocate_oasisgrids
  public :: oasisxregcm_def_field, oasisxregcm_end_def
  public :: oasisxregcm_rcv, fill_ocean
  public :: oasisxregcm_snd

  contains

  ! call initialization OASIS subroutines
  subroutine oasisxregcm_init(localComm)
    implicit none
    integer(ik4), intent(out) :: localComm
    character(len=*), parameter :: sub_name = 'oasisxregcm_init'
    !--------------------------------------------------------------------------
    ! initialize component
    call oasis_init_comp(comp_id,comp_name,ierror)
    if ( ierror /= 0 ) then
      write(stderr,*) 'oasis_init_comp abort compid ', comp_id
      call oasis_abort(comp_id,sub_name,'Problem in oasis_init_comp call')
    end if
    ! get the mpi communicator
    call oasis_get_localcomm(localComm,ierror)
    if ( ierror /= 0 ) then
      write(stderr,*) 'oasis_get_localcomm abort compid ', comp_id
      call oasis_abort(comp_id,sub_name,'Problem in oasis_get_localcomm call')
    end if
  end subroutine oasisxregcm_init

  ! terminate OASIS
  subroutine oasisxregcm_finalize
    implicit none
    character(len=*), parameter :: sub_name = 'oasisxregcm_finalize'
    !--------------------------------------------------------------------------
    call oasis_terminate(ierror)
    if ( ierror /= 0 ) then
      write(stderr,*) 'oasis_terminate abort compid ', comp_id
      call oasis_abort(comp_id,sub_name,'Problem in oasis_terminate call')
    end if
    !call oasisxregcm_endscreen
  end subroutine oasisxregcm_finalize

  ! initialize a type(infogrd) variable
  subroutine oasisxregcm_setup_grid(grd,naNM,naWM,j1,j2,i1,i2, &
                                                  ja,jb,ia,ib,nc)
    implicit none
    character(len=4), intent(in) :: naNM, naWM
    integer, intent(in) :: j1, j2, i1, i2, &
                            ja, jb, ia, ib, nc
    type(infogrd), allocatable, intent(out) :: grd
    !--------------------------------------------------------------------------
    allocate(grd)
    ! names
    grd%naNM = naNM
    grd%naWM = naWM
    ! local indexes
    grd%j1 = j1
    grd%j2 = j2
    grd%jll = j2 - j1 + 1
    grd%i1 = i1
    grd%i2 = i2
    grd%ill = i2 - i1 + 1
    ! global indexes
    grd%ja = ja
    grd%jb = jb
    grd%jgl = jb - ja + 1
    grd%ia = ia
    grd%ib = ib
    grd%igl = ib - ia + 1
    ! number of corner per cell
    grd%nc = nc
  end subroutine oasisxregcm_setup_grid

  ! initialize a type(infofld) variable
  subroutine oasisxregcm_setup_field(fld, na, grd, array, init_val)
    implicit none
    character(len=*), intent(in) :: na
    type(infogrd), target, intent(in) :: grd
    real(rkx), intent(in), optional :: init_val
    type(infofld), allocatable, intent(out) :: fld
    real(rkx), dimension(:,:), allocatable, intent(out), optional :: array
    !--------------------------------------------------------------------------
    allocate(fld)
    allocate(character(len=len(na)) :: fld%na)
    fld%na = na
    fld%grd => grd
    if ( present(array) ) then
      allocate(array( grd%jll, grd%ill ))
      if ( present(init_val) ) then
        array(:,:) = init_val
      else
        array(:,:) = 0.0
      end if
    end if
  end subroutine oasisxregcm_setup_field

  ! deallocate a type(infofld) variable
  subroutine oasisxregcm_deallocate_field(fld, array)
    implicit none
    type(infofld), allocatable, intent(inout) :: fld
    real(rkx), dimension(:,:), allocatable, intent(inout), optional :: array
    !--------------------------------------------------------------------------
    !
    if ( allocated(fld) ) then
      nullify(fld%grd)
      deallocate(fld)
    end if
    if ( present(array) ) then
      if ( allocated(array) ) deallocate(array)
    end if
  end subroutine oasisxregcm_deallocate_field

  ! define an OASIS partition
  subroutine oasisxregcm_def_partition(grd)
    implicit none
    type(infogrd), intent(inout) :: grd
    integer, dimension(:), allocatable :: il_paral ! OASIS partition instructions
    character(len=*), parameter :: sub_name = 'oasisxregcm_def_partition'
    !--------------------------------------------------------------------------
#ifdef DEBUG
    write(ndebug,*) oasis_prefix, '>> ', grd%naNM, '/', grd%naWM, ' partitions'
#endif
    call oasisxregcm_box_partition(il_paral, grd)
    call oasis_def_partition(grd%id, il_paral, ierror)
    if ( ierror /= 0 ) then
      write(stderr,*) 'oasis_def_partition (', grd%naNM, '/', grd%naWM, ') abort compid ', comp_id
      call oasis_abort(comp_id, sub_name, &
        'Problem in oasis_def_partition call for '//grd%naNM//'/'//grd%naWM)
    end if
#ifdef DEBUG
    write(ndebug,"(' ',A,A,I2)") oasis_prefix, '-- id: ', grd%id
#endif

  contains

  ! fill the OASIS partition instruction array in the case
  ! of a 'box' partition:
  ! paral(1) = 2 indicates a box partition
  ! paral(2) is the global offset (here lower left)
  ! paral(3) is the local extent in x
  ! paral(4) is the local extent in y
  ! paral(5) is the global extent in x
  subroutine oasisxregcm_box_partition(paral, grd)
    implicit none
    type(infogrd), intent(in) :: grd
    integer, dimension(:), allocatable, intent(out) :: paral
    !--------------------------------------------------------------------------
    allocate(paral(5))
    paral(1) = 2
    paral(2) = (grd%i1 - grd%ia) * grd%jgl + (grd%j1 - grd%ja)
    paral(3) = grd%jll
    paral(4) = grd%ill
    paral(5) = grd%jgl
#ifdef DEBUG
    write(ndebug,*) oasis_prefix, 'Box definition:'
    write(ndebug,*) oasis_prefix, 'gl off |   lo x |   lo y |   gl x'
    write(ndebug,"(' ',A,3(I6,' | '),I6)") oasis_prefix, paral(2), paral(3), paral(4), paral(5)
#endif
  end subroutine oasisxregcm_box_partition

  end subroutine oasisxregcm_def_partition

  ! allocate oasisgrids temporary pointers
  subroutine oasisxregcm_allocate_oasisgrids(lon,lat,clon,clat,srf,mask,jsize,isize,csize)
    implicit none
    integer(ik4), intent(in) :: jsize, isize, csize
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: lon, lat, srf
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: clon, clat
    integer(ik4), pointer, contiguous, dimension(:,:), intent(inout) :: mask
    !--------------------------------------------------------------------------
    allocate(lon(jsize,isize),stat=ierror)
    if ( ierror /= 0 ) then
      write(stderr,*) 'error allocating lon'
      call fatal(__FILE__,__LINE__,'OASIS GRID')
    end if
    allocate(lat(jsize,isize),stat=ierror)
    if ( ierror /= 0 ) then
      write(stderr,*) 'error allocating lat'
      call fatal(__FILE__,__LINE__,'OASIS GRID')
    end if
    allocate(clon(jsize,isize,csize),stat=ierror)
    if ( ierror /= 0 ) then
      write(stderr,*) 'error allocating clon'
      call fatal(__FILE__,__LINE__,'OASIS GRID')
    end if
    allocate(clat(jsize,isize,csize),stat=ierror)
    if ( ierror /= 0 ) then
      write(stderr,*) 'error allocating clat'
      call fatal(__FILE__,__LINE__,'OASIS GRID')
    end if
    allocate(srf(jsize,isize),stat=ierror)
    if ( ierror /= 0 ) then
      write(stderr,*) 'error allocating srf'
      call fatal(__FILE__,__LINE__,'OASIS GRID')
    end if
    allocate(mask(jsize,isize),stat=ierror)
    if ( ierror /= 0 ) then
      write(stderr,*) 'error allocating mask'
      call fatal(__FILE__,__LINE__,'OASIS GRID')
    end if
  end subroutine oasisxregcm_allocate_oasisgrids

  ! return the surface of the mesh in square meters
  ! (only for number of corner == 4)
  real(rkx) function srf_sqm(clon,clat)
    implicit none
    real(rkx), intent(in), dimension(4) :: clon, &
                                             clat ! degree coordinates of the corners
    real(rkx) :: a, b, h, lat_a, lat_b
    real(rk8), parameter :: factor = degrad * earthrad ! m.deg^(-1)
    !------------------------------------------------------------------------
    ! trapezoid area
    ! base 1: north side
    lat_a = ( clat(1) + clat(2) ) / 2
    a = (clon(1) - clon(2)) * factor * cos( degrad*lat_a )
    ! base 2: south side
    lat_b = ( clat(3) + clat(4) ) / 2
    b = (clon(4) - clon(3)) * factor * cos( degrad*lat_b )
    ! height
    h = abs(lat_a - lat_b) * factor
    ! area
    srf_sqm = h * (a+b) / 2
  end function srf_sqm

  ! give to OASIS the information about a specific grid for writing
  subroutine oasisxregcm_write_oasisgrids(grd,lon,lat,clon,clat,srf,mask)
    implicit none
    type(infogrd), intent(in) :: grd
    real(rkx), pointer, contiguous, dimension(:,:), intent(in) :: lon, lat, srf
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(in) :: clon, clat
    integer(ik4), pointer, contiguous, dimension(:,:), intent(in) :: mask
    integer(ik4), allocatable, dimension(:,:) :: mask0
    !--------------------------------------------------------------------------
    ! No Mask
#ifdef DEBUG
    write(ndebug,*) oasis_prefix, '>>       name: ', grd%naNM
    write(ndebug,"(' ',A,A,I4,A,I4)") oasis_prefix, '-- dimensions: ', grd%jgl, 'x', grd%igl
#endif
    call oasis_write_grid(grd%naNM,grd%jgl,grd%igl, &
         lon(grd%ja:grd%jb, grd%ia:grd%ib), lat(grd%ja:grd%jb, grd%ia:grd%ib))
    !
    call oasis_write_corner(grd%naNM,grd%jgl,grd%igl,grd%nc, &
         clon(grd%ja:grd%jb, grd%ia:grd%ib, :), clat(grd%ja:grd%jb, grd%ia:grd%ib, :))
    !
    call oasis_write_area(grd%naNM,grd%jgl,grd%igl, &
         srf(grd%ja:grd%jb, grd%ia:grd%ib))
    !
    allocate(mask0( lbound(mask,1):ubound(mask,1), lbound(mask,2):ubound(mask,2) ))
    mask0(:,:) = 0
    call oasis_write_mask(grd%naNM,grd%jgl,grd%igl, &
         mask0(grd%ja:grd%jb, grd%ia:grd%ib))
    deallocate(mask0)
    !
    ! With Mask
#ifdef DEBUG
    write(ndebug,*) oasis_prefix, '>>       name: ', grd%naWM
    write(ndebug,"(' ',A,A,I4,A,I4)") oasis_prefix, '-- dimensions: ', grd%jgl, 'x', grd%igl
#endif
    call oasis_write_grid(grd%naWM,grd%jgl,grd%igl, &
         lon(grd%ja:grd%jb, grd%ia:grd%ib), lat(grd%ja:grd%jb, grd%ia:grd%ib))
    !
    call oasis_write_corner(grd%naWM,grd%jgl,grd%igl,grd%nc, &
         clon(grd%ja:grd%jb, grd%ia:grd%ib, :), clat(grd%ja:grd%jb, grd%ia:grd%ib, :))
    !
    call oasis_write_area(grd%naWM,grd%jgl,grd%igl, &
         srf(grd%ja:grd%jb, grd%ia:grd%ib))
    !
    call oasis_write_mask(grd%naWM,grd%jgl,grd%igl, &
         mask(grd%ja:grd%jb, grd%ia:grd%ib))
    !
  end subroutine oasisxregcm_write_oasisgrids

  ! deallocate oasisgrids temporary pointers
  subroutine oasisxregcm_deallocate_oasisgrids(lon,lat,clon,clat,srf,mask)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: lon, lat, srf
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: clon, clat
    integer(ik4), pointer, contiguous, dimension(:,:), intent(inout) :: mask
    !--------------------------------------------------------------------------
    if ( associated(lon) ) nullify(lon)
    if ( associated(lat) ) nullify(lat)
    if ( associated(clon) ) nullify(clon)
    if ( associated(clat) ) nullify(clat)
    if ( associated(srf) ) nullify(srf)
    if ( associated(mask) ) nullify(mask)
  end subroutine oasisxregcm_deallocate_oasisgrids

  ! define OASIS variables to be used in the coupling
  subroutine oasisxregcm_def_field(fld, kinout)
    implicit none
    integer, intent(in) :: kinout ! OASIS_Out or OASIS_In
    type(infofld), intent(inout) :: fld ! field information
    integer, dimension(2) :: var_nodims, & ! not useful but still
                            var_actual_shape ! arguments of oasis_def_var
    integer :: var_type ! type of coupling field array and
    character(len=*), parameter :: sub_name = 'oasisxregcm_def_field'
#ifdef DEBUG
    character(len=6) :: kinout_char
#endif
    data var_nodims / 1, 1 /
    data var_actual_shape / 1, 1 /
    data var_type / OASIS_Real / ! always the case?
    !--------------------------------------------------------------------------
#ifdef DEBUG
    if ( kinout == OASIS_Out ) then
      kinout_char = 'export'
    else
      kinout_char = 'import'
    end if
    write(ndebug,*) oasis_prefix, '>>      name: ', fld%na
    write(ndebug,*) oasis_prefix, '--      grid: ', fld%grd%naNM, '/', fld%grd%naWM
    write(ndebug,*) oasis_prefix, '-- direction: ', kinout_char
    !write(ndebug,*) oasis_prefix, '--      type: ', 'real'
#endif
    call oasis_def_var(fld%id, fld%na, fld%grd%id, &
                       var_nodims, kinout, var_actual_shape,     &
                       var_type, ierror)
    if ( ierror /= 0 ) then
      write(stderr,*) 'oasis_def_var (', fld%na, ') abort compid ', comp_id
      call oasis_abort(comp_id, sub_name, &
        'Problem in oasis_def_var call for field '//fld%na)
    end if
#ifdef DEBUG
    write(ndebug,"(' ',A,A,I2)") oasis_prefix, '--        id: ', fld%id
#endif
  end subroutine oasisxregcm_def_field

  ! terminate the definition phase (posterior to partition and variable definitions)
  subroutine oasisxregcm_end_def
    implicit none
    character(len=*), parameter :: sub_name = 'oasisxregcm_end_def'
    !--------------------------------------------------------------------------
    call oasis_enddef(ierror)
    if ( ierror /= 0 ) then
      write(stderr,*) 'oasis_enddef abort compid ', comp_id
      call oasis_abort(comp_id,sub_name,'Problem in oasis_enddef call')
    end if
  end subroutine oasisxregcm_end_def

  ! receive a single field with id fld_id through oasis_get
  subroutine oasisxregcm_rcv(array,fld,time,l_act)
    implicit none
    type(infofld), intent(in) :: fld
    integer(ik4), intent(in) :: time ! execution time
    real(rkx), dimension(:,:), intent(out) :: array
    logical, intent(out) :: l_act
    character(len=*), parameter :: sub_name  = 'oasisxregcm_rcv'
    !--------------------------------------------------------------------------
    call oasis_get(fld%id,time,array,ierror)
    if ( ierror .ne. OASIS_Ok .and. ierror .lt. OASIS_Recvd ) then
      write(stderr,*) 'oasis_get (', fld%na, ') abort compid ', comp_id
      call oasis_abort(comp_id,sub_name,'Problem in oasis_get call for '//fld%na)
    end if
    l_act = ( ( ierror == OASIS_Recvd ) .or. &
              ( ierror >= OASIS_FromRest .and. ierror <= OASIS_FromRestOut ) .or. &
              ( ierror == OASIS_WaitGroup ) )
#ifdef DEBUG
    if ( l_act ) then
      if ( time == 0 .or. debug_level > 0 ) then
        write(ndebug,"(' ',A,A11,A,ES12.5E2,A,A8,A)") oasis_prefix, &
        getput_status(ierror), ': ', minval(array), ' < ', fld%na, ' < '
        write(ndebug,"(' ',A,A11,A,ES12.5E2)") oasis_prefix, &
        '',                    '  ', maxval(array)
      end if
    end if
#endif
  end subroutine oasisxregcm_rcv

  ! fill the ocean parts of array_out with array_in
  subroutine fill_ocean_2d(array_out,array_in,lndcat,grd)
    implicit none
    real(rkx), dimension(:,:), intent(in) :: array_in
    real(rkx), dimension(:,:), pointer, contiguous, intent(in) :: lndcat
    type(infogrd), intent(in) :: grd
    real(rkx), dimension(:,:), pointer, contiguous, intent(inout) :: array_out
    integer(ik4) :: i, j, ishift, jshift
    !--------------------------------------------------------------------------
    ! It seems that whether it's on crosses or dots,
    ! mddom%lndcat is used.
    do i = grd%i1, grd%i2
      ishift = i - grd%i1 + 1
      do j = grd%j1, grd%j2
        jshift = j - grd%j1 + 1
        if ( isocean(lndcat(j,i)) ) array_out(j,i) = array_in(jshift,ishift)
      end do
    end do
  end subroutine fill_ocean_2d

  ! fill the ocean parts of array_out with array_in (with a subgrid dimension)
  subroutine fill_ocean_3d(array_out,array_in,lndcat,grd)
    implicit none
    real(rkx), dimension(:,:), intent(in) :: array_in
    real(rkx), dimension(:,:), pointer, contiguous, intent(in) :: lndcat
    type(infogrd), intent(in) :: grd
    real(rkx), dimension(:,:,:), pointer, contiguous, intent(inout) :: array_out
    integer(ik4) :: i, j, ishift, jshift
    !--------------------------------------------------------------------------
    ! It seems that whether it's on crosses or dots,
    ! mddom%lndcat is used.
    do i = grd%i1, grd%i2
      ishift = i - grd%i1 + 1
      do j = grd%j1, grd%j2
        jshift = j - grd%j1 + 1
        if ( isocean(lndcat(j,i)) ) array_out(:,j,i) = array_in(jshift,ishift)
      end do
    end do
  end subroutine fill_ocean_3d

  ! send a single field with id fld_id through oasis_put
  subroutine oasisxregcm_snd(array,fld,time,write_out)
    implicit none
    real(rkx), dimension(:,:), intent(in) :: array
    type(infofld), intent(in) :: fld ! field information
    integer(ik4), intent(in) :: time ! execution time
    logical, intent(in) :: write_out
    character(len=*), parameter :: sub_name = 'oasisxregcm_snd'
    logical :: l_act
    !--------------------------------------------------------------------------
    call oasis_put(fld%id,time,array,ierror,write_restart=write_out)
    if ( ierror .ne. OASIS_Ok .and. ierror .lt. OASIS_Sent ) then
      write(stderr,*) 'oasis_put (', fld%na, ') abort compid ', comp_id
      call oasis_abort(comp_id,sub_name,'Problem in oasis_put call for '//fld%na)
    end if
    l_act = ( ( ierror == OASIS_Sent ) .or. &
              ( ierror >= OASIS_ToRest .and. ierror <= OASIS_ToRestOut ) .or. &
              ( ierror == OASIS_WaitGroup ) )
#ifdef DEBUG
    if ( l_act ) then
      if ( debug_level > 0 ) then
        write(ndebug,"(' ',A,A11,A,ES12.5E2,A,A8,A)") oasis_prefix, &
        getput_status(ierror), ': ', minval(array), ' < ', fld%na, ' < '
        write(ndebug,"(' ',A,A11,A,ES12.5E2)") oasis_prefix, &
        '',                    '  ', maxval(array)
      end if
    end if
#endif
  end subroutine oasisxregcm_snd

end module mod_oasis_generic
!
#endif
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

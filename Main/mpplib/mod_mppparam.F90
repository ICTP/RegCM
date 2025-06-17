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

module mod_mppparam

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_constants
  use mod_runparams, only : namelistfile, prgname
  use mod_mpmessage
  use mod_memutil
  use mod_date
  use mod_stdio
  use netcdf
  use mod_regcm_types
#ifndef MPI_SERIAL
  use mpi
#endif

  implicit none

  private

  public :: get_cartcomm

  integer(ik4), public :: global_cross_istart, global_cross_iend
  integer(ik4), public :: global_cross_jstart, global_cross_jend
  integer(ik4), public :: global_dot_istart, global_dot_iend
  integer(ik4), public :: global_dot_jstart, global_dot_jend

  logical, parameter :: lreorder = .false.

  type(masked_comm), public :: lndcomm
  type(masked_comm), public :: ocncomm

  integer(ik4), public, parameter :: iocpu = 0 ! The id of the cpu doing I/O
  integer(ik4), public, parameter :: italk = 0 ! Who is doing the print ?

#ifdef MPI_SERIAL
  include 'mpif.h'
  integer(ik4) mpi_status_ignore(mpi_status_size)
  integer(ik4) mpi_statuses_ignore(mpi_status_size,4)
  integer(ik4), parameter :: mpi_proc_null = -2
  integer(ik4), parameter :: mpi_comm_type_shared = 0
  interface mpi_sendrecv
    module procedure mpi_sendrecvr8
    module procedure mpi_sendrecvr4
  end interface mpi_sendrecv
#endif

  public :: set_nproc, broadcast_params

  integer(ik4) :: cartesian_communicator
  !integer(ik4) :: node_local_communicator
  integer(ik4) :: ccid, ccio

  integer(ik4), public :: ncout_mpi_info = mpi_info_null

  type grid_nc_var2d
    character(len=64) :: varname
    integer(ik4) :: irec = -1
    integer(ik4) :: ncid = -1
    integer(ik4) :: varid = -1
    integer(ik4) :: nx = 0
    integer(ik4) :: ny = 0
    integer(ik4) :: mynx1 = 0
    integer(ik4) :: mynx2 = 0
    integer(ik4) :: myny1 = 0
    integer(ik4) :: myny2 = 0
    real(rk8), pointer, contiguous, dimension(:,:) :: val => null()
    real(rk8), pointer, contiguous, dimension(:,:) :: iobuf => null()
  end type grid_nc_var2d

  type grid_nc_var3d
    character(len=64) :: varname
    integer(ik4) :: irec = -1
    integer(ik4) :: ncid = -1
    integer(ik4) :: varid = -1
    integer(ik4) :: nx = 0
    integer(ik4) :: ny = 0
    integer(ik4) :: mynx1 = 0
    integer(ik4) :: mynx2 = 0
    integer(ik4) :: myny1 = 0
    integer(ik4) :: myny2 = 0
    integer(ik4) :: nz = 0
    real(rk8), pointer, contiguous, dimension(:,:,:) :: val => null()
    real(rk8), pointer, contiguous, dimension(:,:,:) :: iobuf => null()
  end type grid_nc_var3d

  type grid_nc_var4d
    character(len=64) :: varname
    integer(ik4) :: irec = -1
    integer(ik4) :: ncid = -1
    integer(ik4) :: varid = -1
    integer(ik4) :: nx = 0
    integer(ik4) :: ny = 0
    integer(ik4) :: mynx1 = 0
    integer(ik4) :: mynx2 = 0
    integer(ik4) :: myny1 = 0
    integer(ik4) :: myny2 = 0
    integer(ik4) :: nz = 0
    integer(ik4) :: nl = 0
    real(rk8), pointer, contiguous, dimension(:,:,:,:) :: val => null()
    real(rk8), pointer, contiguous, dimension(:,:,:,:) :: iobuf => null()
  end type grid_nc_var4d

  public :: grid_nc_var2d, grid_nc_var3d, grid_nc_var4d

  interface exchange_array
    module procedure exchange_array_r8, &
                     exchange_array_r4
  end interface exchange_array

  interface cyclic_exchange_array
    module procedure cyclic_exchange_array_r8, &
                     cyclic_exchange_array_r4
  end interface cyclic_exchange_array

  interface grid_nc_create
    module procedure grid_nc_create_var2d, &
                     grid_nc_create_var3d, &
                     grid_nc_create_var4d
  end interface grid_nc_create

  interface grid_nc_write
    module procedure grid_nc_write_var2d, &
                     grid_nc_write_var3d, &
                     grid_nc_write_var4d
  end interface grid_nc_write

  interface grid_nc_destroy
    module procedure grid_nc_destroy_var2d, &
                     grid_nc_destroy_var3d, &
                     grid_nc_destroy_var4d
  end interface grid_nc_destroy

  public :: grid_nc_create, grid_nc_write, grid_nc_destroy

  interface grid_distribute
    module procedure real8_2d_distribute,   &
                     real8_3d_distribute,   &
                     real8_4d_distribute,   &
                     real4_2d_distribute,   &
                     real4_3d_distribute,   &
                     real4_4d_distribute,   &
                     integer_2d_distribute, &
                     integer_3d_distribute, &
                     integer_4d_distribute, &
                     logical_2d_distribute, &
                     logical_3d_distribute, &
                     logical_4d_distribute
  end interface grid_distribute

  interface subgrid_distribute
    module procedure real8_2d_sub_distribute,   &
                     real8_3d_sub_distribute,   &
                     real4_2d_sub_distribute,   &
                     real4_3d_sub_distribute,   &
                     integer_2d_sub_distribute, &
                     integer_3d_sub_distribute, &
                     logical_2d_sub_distribute
  end interface subgrid_distribute

  interface grid_collect
    module procedure real8_2d_collect,    &
                     real8_2d_3d_collect, &
                     real8_3d_collect,    &
                     real8_3d_2d_collect, &
                     real8_4d_collect,    &
                     real8_4d_2d_collect, &
                     real4_2d_collect,    &
                     real4_2d_3d_collect, &
                     real4_3d_collect,    &
                     real4_3d_2d_collect, &
                     real4_4d_collect,    &
                     real4_4d_2d_collect, &
                     integer_2d_collect,  &
                     integer_3d_collect,  &
                     integer_4d_collect,  &
                     logical_2d_collect
  end interface grid_collect

  interface subgrid_collect
    module procedure real8_2d_sub_collect,   &
                     real8_3d_sub_collect,   &
                     real4_2d_sub_collect,   &
                     real4_3d_sub_collect,   &
                     integer_2d_sub_collect, &
                     integer_3d_sub_collect, &
                     logical_2d_sub_collect
  end interface subgrid_collect

  interface exchange
    module procedure real8_2d_exchange, &
                     real8_3d_exchange, &
                     real8_4d_exchange, &
                     real4_2d_exchange, &
                     real4_3d_exchange, &
                     real4_4d_exchange
  end interface exchange

  interface exchange_lrbt
    module procedure real8_2d_exchange_left_right_bottom_top, &
                     real8_3d_exchange_left_right_bottom_top, &
                     real8_4d_exchange_left_right_bottom_top, &
                     real4_2d_exchange_left_right_bottom_top, &
                     real4_3d_exchange_left_right_bottom_top, &
                     real4_4d_exchange_left_right_bottom_top
  end interface exchange_lrbt

  interface exchange_lr
    module procedure real8_2d_exchange_left_right, &
                     real8_3d_exchange_left_right, &
                     real8_4d_exchange_left_right, &
                     real4_2d_exchange_left_right, &
                     real4_3d_exchange_left_right, &
                     real4_4d_exchange_left_right
  end interface exchange_lr

  interface exchange_bt
    module procedure real8_2d_exchange_bottom_top, &
                     real8_3d_exchange_bottom_top, &
                     real8_4d_exchange_bottom_top, &
                     real4_2d_exchange_bottom_top, &
                     real4_3d_exchange_bottom_top, &
                     real4_4d_exchange_bottom_top
  end interface exchange_bt

  interface exchange_lb
    module procedure real8_2d_exchange_left_bottom, &
                     real8_3d_exchange_left_bottom, &
                     real8_4d_exchange_left_bottom, &
                     real4_2d_exchange_left_bottom, &
                     real4_3d_exchange_left_bottom, &
                     real4_4d_exchange_left_bottom
  end interface exchange_lb

  interface exchange_rt
    module procedure real8_2d_exchange_right_top, &
                     real8_3d_exchange_right_top, &
                     real8_4d_exchange_right_top, &
                     real4_2d_exchange_right_top, &
                     real4_3d_exchange_right_top, &
                     real4_4d_exchange_right_top
  end interface exchange_rt

  interface exchange_bdy_lr
    module procedure real8_bdy_exchange_left_right, &
                     real4_bdy_exchange_left_right
  end interface exchange_bdy_lr

  interface exchange_bdy_bt
    module procedure real8_bdy_exchange_bottom_top, &
                     real4_bdy_exchange_bottom_top
  end interface exchange_bdy_bt

  interface grid_fill
    module procedure real8_2d_grid_fill_extend1, &
                     real8_2d_grid_fill_extend2, &
                     real4_2d_grid_fill_extend1, &
                     real4_2d_grid_fill_extend2
  end interface grid_fill

  interface bcast
    module procedure bcast_logical,           &
                     bcast_int4,              &
                     bcast_int8,              &
                     bcast_real4,             &
                     bcast_real8,             &
#ifdef QUAD_PRECISION
                     bcast_real16,            &
#endif
                     bcast_arr_logical,       &
                     bcast_arr_character,     &
                     bcast_arr_text_list,     &
                     bcast_arr_int4,          &
                     bcast_arr_int8,          &
                     bcast_arr_real4,         &
                     bcast_arr_real8,         &
                     bcast_matr_real8,        &
                     bcast_matr_real4,        &
                     bcast_rcm_time_and_date, &
                     bcast_arr_rcm_time_and_date
  end interface bcast

#ifdef QUAD_PRECISION
  interface sumall
    module procedure sumall_real16, &
                     sumall_real8, &
                     sumall_real4, &
                     sumall_int4,  &
                     sumall_int4_array
  end interface sumall
#else
  interface sumall
    module procedure sumall_real8, &
                     sumall_real4, &
                     sumall_int4,  &
                     sumall_int4_array
  end interface sumall
#endif

#ifndef USE_MPI3
  interface send_array
    module procedure send_array_logical, &
                     send_array_int4,    &
                     send_array_real4,   &
                     send_array_real8
  end interface send_array

  interface recv_array
    module procedure recv_array_logical, &
                     recv_array_int4,    &
                     recv_array_real4,   &
                     recv_array_real8
  end interface recv_array
#endif

  interface c2l_ss
    module procedure cartesian_to_linear_integer_subgrid_subgrid,  &
                     cartesian_to_linear_logical_subgrid_subgrid,  &
                     cartesian_to_linear_real8_subgrid_subgrid,    &
                     cartesian_to_linear_real8_subgrid_subgrid_4d, &
                     cartesian_to_linear_real4_subgrid_subgrid,    &
                     cartesian_to_linear_real4_subgrid_subgrid_4d
  end interface c2l_ss

  interface c2l_gs
    module procedure cartesian_to_linear_integer_grid_subgrid, &
                     cartesian_to_linear_logical_grid_subgrid, &
                     cartesian_to_linear_real8_grid_subgrid,   &
                     cartesian_to_linear_real4_grid_subgrid
  end interface c2l_gs

  interface l2c_ss
    module procedure linear_to_cartesian_integer_subgrid_subgrid,  &
                     linear_to_cartesian_logical_subgrid_subgrid,  &
                     linear_to_cartesian_real8_subgrid_subgrid,    &
                     linear_to_cartesian_real8_subgrid_subgrid_4d, &
                     linear_to_cartesian_real4_subgrid_subgrid,    &
                     linear_to_cartesian_real4_subgrid_subgrid_4d
  end interface l2c_ss

  interface glb_c2l_ss
    module procedure global_to_linear_integer_subgrid_subgrid,       &
                     global_to_linear_logical_subgrid_subgrid,       &
                     global_to_linear_real8_subgrid_subgrid,         &
                     global_to_linear_real8_subgrid_subgrid_4d,      &
                     global_to_linear_real4_subgrid_subgrid,         &
                     global_to_linear_real4_subgrid_subgrid_4d,      &
                     global_to_linear_real4_real8_subgrid_subgrid,   &
                     global_to_linear_real4_real8_subgrid_subgrid_4d
  end interface glb_c2l_ss

  interface glb_c2l_gs
    module procedure global_to_linear_integer_grid_subgrid,    &
                     global_to_linear_logical_grid_subgrid,    &
                     global_to_linear_real8_grid_subgrid,      &
                     global_to_linear_real4_grid_subgrid,      &
                     global_to_linear_real4_real8_grid_subgrid
  end interface glb_c2l_gs

  interface glb_l2c_ss
    module procedure linear_to_global_integer_subgrid_subgrid,      &
                     linear_to_global_logical_subgrid_subgrid,      &
                     linear_to_global_real8_subgrid_subgrid,        &
                     linear_to_global_real8_subgrid_subgrid_4d,     &
                     linear_to_global_real4_subgrid_subgrid,        &
                     linear_to_global_real4_subgrid_subgrid_4d,     &
                     linear_to_global_real8_real4_subgrid_subgrid,  &
                     linear_to_global_real8_real4_subgrid_subgrid_4d
  end interface glb_l2c_ss

  interface reorder_glb_subgrid
     module procedure reorder_logical_global_subgrid_2d
  end interface reorder_glb_subgrid

  interface reorder_subgrid
    module procedure reorder_subgrid_2d_real8,   &
                     reorder_subgrid_2d3d_real8, &
                     reorder_subgrid_3d_real8,   &
                     reorder_subgrid_4d_real8,   &
                     reorder_subgrid_2d_real4,   &
                     reorder_subgrid_2d3d_real4, &
                     reorder_subgrid_3d_real4,   &
                     reorder_subgrid_4d_real4,   &
                     reorder_subgrid_2d_logical
  end interface reorder_subgrid

  interface reorder_add_subgrid
    module procedure reorder_add_subgrid_2d_real8,  &
                     reorder_add_subgrid_2d3d_real8, &
                     reorder_add_subgrid_3d_real8,   &
                     reorder_add_subgrid_2d_real4,   &
                     reorder_add_subgrid_2d3d_real4, &
                     reorder_add_subgrid_3d_real4
  end interface reorder_add_subgrid

  interface input_reorder
    module procedure input_reorder_real8, &
                     input_reorder_real4
  end interface input_reorder

  interface mypack
    module procedure mypack_logical_grid,        &
                     mypack_logical_subgrid,     &
                     mypack_integer_grid,        &
                     mypack_integer_subgrid,     &
                     mypack_real8_grid,          &
                     mypack_real8_subgrid,       &
                     mypack_real8_subgrid_4d,    &
                     mypack_real8_subgrid_slice, &
                     mypack_real4_grid,          &
                     mypack_real4_subgrid,       &
                     mypack_real4_subgrid_4d,    &
                     mypack_real4_subgrid_slice
  end interface mypack

  interface myunpack
    module procedure myunpack_logical_grid,              &
                     myunpack_logical_subgrid,           &
                     myunpack_integer_grid,              &
                     myunpack_integer_subgrid,           &
                     myunpack_real8_grid,                &
                     myunpack_real8_subgrid,             &
                     myunpack_real8_subgrid_4d,          &
                     myunpack_real8_subgrid_slice,       &
                     myunpack_real4_grid,                &
                     myunpack_real4_subgrid,             &
                     myunpack_real4_subgrid_4d,          &
                     myunpack_real4_subgrid_slice,       &
                     myunpack_real8_real4_grid,          &
                     myunpack_real8_real4_subgrid,       &
                     myunpack_real8_real4_subgrid_4d,    &
                     myunpack_real8_real4_subgrid_slice
  end interface myunpack

  interface mypack_global
    module procedure mypack_global_logical_grid,            &
                     mypack_global_logical_subgrid,         &
                     mypack_global_integer_grid,            &
                     mypack_global_integer_subgrid,         &
                     mypack_global_real8_grid,              &
                     mypack_global_real8_subgrid,           &
                     mypack_global_real8_subgrid_4d,        &
                     mypack_global_real8_subgrid_slice,     &
                     mypack_global_real4_grid,              &
                     mypack_global_real4_subgrid,           &
                     mypack_global_real4_subgrid_4d,        &
                     mypack_global_real4_subgrid_slice,     &
                     mypack_global_real4_real8_grid,        &
                     mypack_global_real4_real8_subgrid,     &
                     mypack_global_real4_real8_subgrid_4d,  &
                     mypack_global_real4_real8_subgrid_slice
  end interface mypack_global

  interface myunpack_global
    module procedure myunpack_global_logical_grid,              &
                     myunpack_global_logical_subgrid,           &
                     myunpack_global_integer_grid,              &
                     myunpack_global_integer_subgrid,           &
                     myunpack_global_real8_grid,                &
                     myunpack_global_real8_subgrid,             &
                     myunpack_global_real8_subgrid_4d,          &
                     myunpack_global_real8_subgrid_slice,       &
                     myunpack_global_real4_grid,                &
                     myunpack_global_real4_subgrid,             &
                     myunpack_global_real4_subgrid_4d,          &
                     myunpack_global_real4_subgrid_slice,       &
                     myunpack_global_real8_real4_grid,          &
                     myunpack_global_real8_real4_subgrid,       &
                     myunpack_global_real8_real4_subgrid_4d,    &
                     myunpack_global_real8_real4_subgrid_slice
  end interface myunpack_global

  interface cl_setup
    module procedure cl_setup_real8
    module procedure cl_setup_real4
  end interface cl_setup

  interface maxall
    module procedure maxall_real8
    module procedure maxall_real4
    module procedure maxall_integer4
  end interface maxall

  interface minall
    module procedure minall_real8
    module procedure minall_real4
    module procedure minall_integer4
  end interface minall

  interface meanall
    module procedure meanall_1D_real8
    module procedure meanall_real8
    module procedure meanall_real4
  end interface meanall

  interface cross2dot
    module procedure cross2dot2d
    module procedure cross2dot3d
  end interface cross2dot

  type(model_area), public :: ma

  real(rk8), pointer, contiguous, dimension(:) :: r8vector1
  real(rk8), pointer, contiguous, dimension(:) :: r8vector2
  real(rk4), pointer, contiguous, dimension(:) :: r4vector1
  real(rk4), pointer, contiguous, dimension(:) :: r4vector2
  integer(ik4), pointer, contiguous, dimension(:) :: i4vector1
  integer(ik4), pointer, contiguous, dimension(:) :: i4vector2
  logical, pointer, contiguous, dimension(:) :: lvector1
  logical, pointer, contiguous, dimension(:) :: lvector2
  integer(ik4), dimension(4) :: window
  integer(ik4), pointer, contiguous, dimension(:) :: windows
  integer(ik4), dimension(1), target :: ifake
  integer(ik4), pointer, contiguous, dimension(:) :: wincount
  integer(ik4), pointer, contiguous, dimension(:) :: windispl
  integer(ik4) :: mpierr
  real(rk8), pointer, contiguous, dimension(:,:,:) :: r8subgrid
  real(rk4), pointer, contiguous, dimension(:,:,:) :: r4subgrid
  integer(ik4), pointer, contiguous, dimension(:,:,:) :: i4subgrid
  logical, pointer, contiguous, dimension(:,:,:) :: lsubgrid
  real(rk8), pointer, contiguous, dimension(:,:,:) :: global_r8subgrid
  real(rk4), pointer, contiguous, dimension(:,:,:) :: global_r4subgrid
  integer(ik4), pointer, contiguous, dimension(:,:,:) :: global_i4subgrid
  logical, pointer, contiguous, dimension(:,:,:) :: global_lsubgrid
  real(rk8), pointer, contiguous, dimension(:,:) :: global_r8grid
  real(rk4), pointer, contiguous, dimension(:,:) :: global_r4grid
  integer(ik4), pointer, contiguous, dimension(:,:) :: global_i4grid
  logical, pointer, contiguous, dimension(:,:) :: global_lgrid

  integer(ik4), parameter :: tag_bt = 1 ! FROM bottom TO top
  integer(ik4), parameter :: tag_tb = 2 ! FROM top TO bottom
  integer(ik4), parameter :: tag_lr = 3 ! FROM left TO right
  integer(ik4), parameter :: tag_rl = 4 ! FROM right TO left
  integer(ik4), parameter :: tag_brtl = 5 ! FROM bottomrigth TO topleft
  integer(ik4), parameter :: tag_tlbr = 6 ! FROM topleft TO bottomright
  integer(ik4), parameter :: tag_bltr = 7 ! FROM bottomleft TO topright
  integer(ik4), parameter :: tag_trbl = 8 ! FROM topright TO bottomleft

  public :: exchange, exchange_lb, exchange_rt
  public :: exchange_lr, exchange_bt, exchange_lrbt
  public :: exchange_bdy_lr, exchange_bdy_bt
  public :: grid_distribute, grid_collect, grid_fill
  public :: subgrid_distribute, subgrid_collect
  public :: uvcross2dot, uvdot2cross, cross2dot, psc2psd
  public :: tenxtouvten, uvtentotenx
  public :: bcast, sumall, maxall, minall, meanall
  public :: gather_r, gather_i
  public :: allgather_r, allgather_i
  public :: reorder_subgrid, reorder_glb_subgrid, reorder_add_subgrid
  public :: input_reorder
  public :: trueforall
  public :: allsync
  public :: cl_setup, cl_dispose
  public :: c2l_gs, c2l_ss, l2c_ss
  public :: glb_c2l_gs, glb_c2l_ss, glb_l2c_ss

  logical, save, public :: on_device = .false.

  contains

#ifdef MPI_SERIAL

  subroutine mpi_comm_split_type(comm, split_type, key, info, newcomm, ierror)
    implicit none
    integer(ik4) :: comm, split_type, key, info, newcomm, ierror
  end subroutine mpi_comm_split_type

  subroutine mpi_cart_shift(comm, direction, shift, source, dest, ierror)
    implicit none
    integer(ik4) :: comm, direction, shift, source, dest, ierror
  end subroutine mpi_cart_shift

  subroutine mpi_sendrecvr4(sendbuf, sendcount, sendtype, dest, sendtag, &
                            recvbuf, recvcount, recvtype, source, recvtag, &
                            comm, status, ierror)
    implicit none
    real(rk4), dimension(:) :: sendbuf, recvbuf
    integer(ik4) :: sendcount, sendtype, dest, sendtag
    integer(ik4) :: recvcount, recvtype, source, recvtag, comm
    integer(ik4) :: status(mpi_status_size)
    integer(ik4) :: ierror
    recvbuf(1:recvcount) = sendbuf(1:sendcount)
  end subroutine mpi_sendrecvr4

  subroutine mpi_sendrecvr8(sendbuf, sendcount, sendtype, dest, sendtag, &
                            recvbuf, recvcount, recvtype, source, recvtag, &
                            comm, status, ierror)
    implicit none
    real(rk8), dimension(:) :: sendbuf, recvbuf
    integer(ik4) :: sendcount, sendtype, dest, sendtag
    integer(ik4) :: recvcount, recvtype, source, recvtag, comm
    integer(ik4) :: status(mpi_status_size)
    integer(ik4) :: ierror
    recvbuf(1:recvcount) = sendbuf(1:sendcount)
  end subroutine mpi_sendrecvr8

  subroutine mpi_cart_create(comm_old,ndims,dims,periods,reorder, &
                             comm_cart,ierror)
    implicit none
    integer(ik4) :: comm_old, ndims, comm_cart, ierror
    integer(ik4), dimension(:) :: dims
    logical :: reorder
    logical, dimension(:) :: periods
  end subroutine mpi_cart_create

  subroutine mpi_cart_coords(comm,rank,maxdims,coords,ierror)
    implicit none
    integer(ik4) :: comm, rank, maxdims, ierror
    integer(ik4), dimension(:) :: coords
  end subroutine mpi_cart_coords

  subroutine mpi_cart_rank(comm,coords,rank,ierror)
    implicit none
    integer(ik4) :: comm, rank, ierror
    integer(ik4), dimension(:) :: coords
  end subroutine mpi_cart_rank
#endif

  subroutine bcast_logical(lval)
    implicit none
    logical, intent(inout) :: lval
    call mpi_bcast(lval,1,mpi_logical,iocpu,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
#endif
  end subroutine bcast_logical

  subroutine bcast_int4(ival)
    implicit none
    integer(ik4), intent(inout) :: ival
    call mpi_bcast(ival,1,mpi_integer4,iocpu,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
#endif
  end subroutine bcast_int4

  subroutine bcast_int8(ival)
    implicit none
    integer(rk8), intent(inout) :: ival
    call mpi_bcast(ival,1,mpi_integer8,iocpu,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
#endif
  end subroutine bcast_int8

  subroutine bcast_real4(rval)
    implicit none
    real(rk4), intent(inout) :: rval
    call mpi_bcast(rval,1,mpi_real4,iocpu,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
#endif
  end subroutine bcast_real4

  subroutine bcast_real8(rval)
    implicit none
    real(rk8), intent(inout) :: rval
    call mpi_bcast(rval,1,mpi_real8,iocpu,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
#endif
  end subroutine bcast_real8

#ifdef QUAD_PRECISION
  subroutine bcast_real16(rval)
    implicit none
    real(rk16), intent(inout) :: rval
    call mpi_bcast(rval,1,mpi_real16,iocpu,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
#endif
  end subroutine bcast_real16
#endif

  subroutine bcast_arr_logical(lval)
    implicit none
    logical, dimension(:), intent(inout) :: lval
    call mpi_bcast(lval,size(lval),mpi_logical,iocpu,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
#endif
  end subroutine bcast_arr_logical

  subroutine bcast_arr_character(cval,is)
    implicit none
    character(len=*), intent(inout) :: cval
    integer(ik4), intent(in) :: is
    call mpi_bcast(cval,is,mpi_character,iocpu,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
#endif
  end subroutine bcast_arr_character

  subroutine bcast_arr_text_list(cval,is)
    implicit none
    character(len=*), intent(inout), dimension(:) :: cval
    integer(ik4), intent(in) :: is
    call mpi_bcast(cval,is*size(cval),mpi_character,iocpu,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
#endif
  end subroutine bcast_arr_text_list

  subroutine bcast_arr_int4(ival)
    implicit none
    integer(ik4), dimension(:), intent(inout) :: ival
    call mpi_bcast(ival,size(ival),mpi_integer4,iocpu,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
#endif
  end subroutine bcast_arr_int4

  subroutine bcast_arr_int8(ival)
    implicit none
    integer(rk8), dimension(:), intent(inout) :: ival
    call mpi_bcast(ival,size(ival),mpi_integer8,iocpu,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
#endif
  end subroutine bcast_arr_int8

  subroutine bcast_arr_real4(rval)
    implicit none
    real(rk4), dimension(:), intent(inout) :: rval
    call mpi_bcast(rval,size(rval),mpi_real4,iocpu,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
#endif
  end subroutine bcast_arr_real4

  subroutine bcast_arr_real8(rval)
    implicit none
    real(rk8), dimension(:), intent(inout) :: rval
    call mpi_bcast(rval,size(rval),mpi_real8,iocpu,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
#endif
  end subroutine bcast_arr_real8

  subroutine bcast_matr_real8(rval)
    implicit none
    real(rk8), dimension(:,:), intent(inout) :: rval
    call mpi_bcast(rval,size(rval,1)*size(rval,2), &
                   mpi_real8,iocpu,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
#endif
  end subroutine bcast_matr_real8

  subroutine bcast_matr_real4(rval)
    implicit none
    real(rk4), dimension(:,:), intent(inout) :: rval
    call mpi_bcast(rval,size(rval,1)*size(rval,2), &
                   mpi_real4,iocpu,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
#endif
  end subroutine bcast_matr_real4

  subroutine bcast_arr_rcm_time_and_date(x)
    implicit none
    type (rcm_time_and_date), dimension(:), intent(inout) :: x
    integer(ik4) :: n
    do n = 1, size(x)
      call bcast(x(n)%calendar)
      call bcast(x(n)%days_from_reference)
      call bcast(x(n)%second_of_day)
    end do
  end subroutine bcast_arr_rcm_time_and_date

  subroutine bcast_rcm_time_and_date(x)
    implicit none
    type (rcm_time_and_date), intent(inout) :: x
    call bcast(x%calendar)
    call bcast(x%days_from_reference)
    call bcast(x%second_of_day)
  end subroutine bcast_rcm_time_and_date

  subroutine trueforall(rlval,rtval)
    implicit none
    logical, intent(in) :: rlval
    logical, intent(out) :: rtval
    call mpi_allreduce(rlval,rtval,1,mpi_logical,mpi_lor,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_allreduce error.')
    end if
#endif
  end subroutine trueforall

#ifdef QUAD_PRECISION
  subroutine sumall_real16(rlval,rtval)
    implicit none
    real(rk16), intent(in) :: rlval
    real(rk16), intent(out) :: rtval
    call mpi_allreduce(rlval,rtval,1,mpi_real16,mpi_sum,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_allreduce error.')
    end if
#endif
  end subroutine sumall_real16
#endif

  subroutine sumall_real8(rlval,rtval)
    implicit none
    real(rk8), intent(in) :: rlval
    real(rk8), intent(out) :: rtval
    call mpi_allreduce(rlval,rtval,1,mpi_real8,mpi_sum,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_allreduce error.')
    end if
#endif
  end subroutine sumall_real8

  subroutine sumall_real4(rlval,rtval)
    implicit none
    real(rk4), intent(in) :: rlval
    real(rk4), intent(out) :: rtval
    call mpi_allreduce(rlval,rtval,1,mpi_real4,mpi_sum,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_allreduce error.')
    end if
#endif
  end subroutine sumall_real4

  subroutine maxall_real8(rlval,rtval)
    implicit none
    real(rk8), intent(in) :: rlval
    real(rk8), intent(out) :: rtval
    call mpi_allreduce(rlval,rtval,1,mpi_real8,mpi_max,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_allreduce error.')
    end if
#endif
  end subroutine maxall_real8

  subroutine maxall_real4(rlval,rtval)
    implicit none
    real(rk4), intent(in) :: rlval
    real(rk4), intent(out) :: rtval
    call mpi_allreduce(rlval,rtval,1,mpi_real4,mpi_max,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_allreduce error.')
    end if
#endif
  end subroutine maxall_real4

  subroutine maxall_integer4(rlval,rtval)
    implicit none
    integer(ik4), intent(in) :: rlval
    integer(ik4), intent(out) :: rtval
    call mpi_allreduce(rlval,rtval,1,mpi_integer4,mpi_max,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_allreduce error.')
    end if
#endif
  end subroutine maxall_integer4

  subroutine meanall_1D_real8(rlval,rtval,elem)
    implicit none
    integer(ik4), intent(in) :: elem
    real(rk8), dimension(:), intent(in) :: rlval
    real(rk8), dimension(:), intent(out) :: rtval
    integer(ik4) :: i
    !$acc data copyin(rlval(1:elem)) copyout(rtval(1:elem))
    !$acc host_data use_device(rlval,rtval)
    call mpi_allreduce(rlval,rtval,elem,mpi_real8,mpi_sum,mycomm,mpierr)
    !$acc end host_data
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_allreduce error.')
    end if
#endif
    do concurrent( i = 1:elem )
      rtval(i) = rtval(i)/real(nproc,rk8)
    end do
    !$acc end data
  end subroutine meanall_1D_real8

  subroutine meanall_real8(rlval,rtval)
    implicit none
    real(rk8), intent(in) :: rlval
    real(rk8), intent(out) :: rtval
    call mpi_allreduce(rlval,rtval,1,mpi_real8,mpi_sum,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_allreduce error.')
    end if
#endif
    rtval = rtval/real(nproc,rk8)
  end subroutine meanall_real8

  subroutine meanall_real4(rlval,rtval)
    implicit none
    real(rk4), intent(in) :: rlval
    real(rk4), intent(out) :: rtval
    call mpi_allreduce(rlval,rtval,1,mpi_real4,mpi_sum,mycomm,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_allreduce error.')
    end if
    rtval = rtval/real(nproc,rk4)
  end subroutine meanall_real4

  subroutine minall_real8(rlval,rtval)
    implicit none
    real(rk8), intent(in) :: rlval
    real(rk8), intent(out) :: rtval
    call mpi_allreduce(rlval,rtval,1,mpi_real8,mpi_min,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_allreduce error.')
    end if
#endif
  end subroutine minall_real8

  subroutine minall_real4(rlval,rtval)
    implicit none
    real(rk4), intent(in) :: rlval
    real(rk4), intent(out) :: rtval
    call mpi_allreduce(rlval,rtval,1,mpi_real4,mpi_min,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_allreduce error.')
    end if
#endif
  end subroutine minall_real4

  subroutine minall_integer4(rlval,rtval)
    implicit none
    integer(ik4), intent(in) :: rlval
    integer(ik4), intent(out) :: rtval
    call mpi_allreduce(rlval,rtval,1,mpi_integer4,mpi_min,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_allreduce error.')
    end if
#endif
  end subroutine minall_integer4

  subroutine sumall_int4(ilval,itval)
    implicit none
    integer(ik4), intent(in) :: ilval
    integer(ik4), intent(out) :: itval
    call mpi_allreduce(ilval,itval,1,mpi_integer4,mpi_sum,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_allreduce error.')
    end if
#endif
  end subroutine sumall_int4

  subroutine sumall_int4_array(ilval,itval)
    implicit none
    integer(ik4), dimension(:), intent(in) :: ilval
    integer(ik4), dimension(:), intent(out) :: itval
    call mpi_allreduce(ilval,itval,size(itval),mpi_integer4, &
                       mpi_sum,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_allreduce error.')
    end if
#endif
  end subroutine sumall_int4_array

#ifndef USE_MPI3
  subroutine send_array_logical(lval,isize,icpu,itag,req)
    implicit none
    integer(ik4), intent(in) :: isize, icpu, itag
    logical, dimension(isize), intent(in) :: lval
    integer(ik4), intent(out) :: req
    call mpi_isend(lval,isize,mpi_logical,icpu,itag, &
                  cartesian_communicator,req,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      write(stderr, *) 'At line :', itag
      call fatal(__FILE__,__LINE__,'mpi_isend error.')
    end if
#endif
  end subroutine send_array_logical

  subroutine send_array_int4(ival,isize,icpu,itag,req)
    implicit none
    integer(ik4), intent(in) :: isize, icpu, itag
    integer(ik4), dimension(isize), intent(in) :: ival
    integer(ik4), intent(out) :: req
    call mpi_isend(ival,isize,mpi_integer4,icpu,itag, &
                  cartesian_communicator,req,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      write(stderr, *) 'At line :', itag
      call fatal(__FILE__,__LINE__,'mpi_isend error.')
    end if
#endif
  end subroutine send_array_int4

  subroutine send_array_real4(rval,isize,icpu,itag,req)
    implicit none
    integer(ik4), intent(in) :: isize, icpu, itag
    real(rk4), dimension(isize), intent(in) :: rval
    integer(ik4), intent(out) :: req
    call mpi_isend(rval,isize,mpi_real4,icpu,itag, &
                  cartesian_communicator,req,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      write(stderr, *) 'At line :', itag
      call fatal(__FILE__,__LINE__,'mpi_isend error.')
    end if
#endif
  end subroutine send_array_real4

  subroutine send_array_real8(rval,isize,icpu,itag,req)
    implicit none
    integer(ik4), intent(in) :: isize, icpu, itag
    real(rk8), dimension(isize), intent(in) :: rval
    integer(ik4), intent(out) :: req
    call mpi_isend(rval,isize,mpi_real8,icpu,itag, &
                  cartesian_communicator,req,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      write(stderr, *) 'At line :', itag
      call fatal(__FILE__,__LINE__,'mpi_isend error.')
    end if
#endif
  end subroutine send_array_real8

  subroutine recv_array_logical(lval,isize,icpu,itag,req)
    implicit none
    integer(ik4), intent(in) :: isize, icpu, itag
    logical, dimension(isize), intent(out) :: lval
    integer(ik4), intent(out) :: req
    call mpi_irecv(lval,isize,mpi_logical,icpu,itag, &
                   cartesian_communicator,req,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      write(stderr, *) 'At line :', itag
      call fatal(__FILE__,__LINE__,'mpi_irecv error.')
    end if
#endif
  end subroutine recv_array_logical

  subroutine recv_array_int4(ival,isize,icpu,itag,req)
    implicit none
    integer(ik4), intent(in) :: isize, icpu, itag
    integer(ik4), dimension(isize), intent(out) :: ival
    integer(ik4), intent(out) :: req
    call mpi_irecv(ival,isize,mpi_integer4,icpu,itag, &
                   cartesian_communicator,req,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      write(stderr, *) 'At line :', itag
      call fatal(__FILE__,__LINE__,'mpi_irecv error.')
    end if
#endif
  end subroutine recv_array_int4

  subroutine recv_array_real4(rval,isize,icpu,itag,req)
    implicit none
    integer(ik4), intent(in) :: isize, icpu, itag
    real(rk4), dimension(isize), intent(out) :: rval
    integer(ik4), intent(out) :: req
    call mpi_irecv(rval,isize,mpi_real4,icpu,itag, &
                   cartesian_communicator,req,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      write(stderr, *) 'At line :', itag
      call fatal(__FILE__,__LINE__,'mpi_irecv error.')
    end if
#endif
  end subroutine recv_array_real4

  subroutine recv_array_real8(rval,isize,icpu,itag,req)
    implicit none
    integer(ik4), intent(in) :: isize, icpu, itag
    real(rk8), dimension(isize), intent(out) :: rval
    integer(ik4), intent(out) :: req
    call mpi_irecv(rval,isize,mpi_real8,icpu,itag, &
                   cartesian_communicator,req,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      write(stderr, *) 'At line :', itag
      call fatal(__FILE__,__LINE__,'mpi_irecv error.')
    end if
#endif
  end subroutine recv_array_real8
#endif

  subroutine exchange_array_r8(rv1,rv2,isize,icpu,tag1,tag2,srq,rrq)
    implicit none
    integer(ik4), intent(in) :: isize, icpu, tag1, tag2
    real(rk8), dimension(isize), intent(in) :: rv1
    real(rk8), dimension(isize), intent(inout) :: rv2
    integer(ik4), intent(out) :: srq, rrq
    call mpi_irecv(rv2,isize,mpi_real8,icpu,tag1, &
                   cartesian_communicator,rrq,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_irecv error.')
    end if
#endif
    call mpi_isend(rv1,isize,mpi_real8,icpu,tag2, &
                  cartesian_communicator,srq,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_isend error.')
    end if
#endif
  end subroutine exchange_array_r8

  subroutine exchange_array_r4(rv1,rv2,isize,icpu,tag1,tag2,srq,rrq)
    implicit none
    integer(ik4), intent(in) :: isize, icpu, tag1, tag2
    real(rk4), dimension(isize), intent(in) :: rv1
    real(rk4), dimension(isize), intent(inout) :: rv2
    integer(ik4), intent(out) :: srq, rrq
    call mpi_irecv(rv2,isize,mpi_real4,icpu,tag1, &
                   cartesian_communicator,rrq,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_irecv error.')
    end if
#endif
    call mpi_isend(rv1,isize,mpi_real4,icpu,tag2, &
                  cartesian_communicator,srq,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_isend error.')
    end if
#endif
  end subroutine exchange_array_r4

  subroutine cyclic_exchange_array_r8(rv1,rv2,isize,icpu1,icpu2,itag)
    implicit none
    real(rk8), pointer, contiguous, dimension(:), intent(in) :: rv1
    real(rk8), pointer, contiguous, dimension(:), intent(inout) :: rv2
    integer(ik4), intent(in) :: isize, icpu1, icpu2, itag
    call mpi_sendrecv(rv1,isize,mpi_real8,icpu1,itag, &
                      rv2,isize,mpi_real8,icpu2,itag, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      write(stderr, *) 'At line :', itag
      call fatal(__FILE__,__LINE__,'mpi_sendrecv error.')
    end if
#endif
  end subroutine cyclic_exchange_array_r8

  subroutine cyclic_exchange_array_r4(rv1,rv2,isize,icpu1,icpu2,itag)
    implicit none
    real(rk4), pointer, contiguous, dimension(:), intent(in) :: rv1
    real(rk4), pointer, contiguous, dimension(:), intent(inout) :: rv2
    integer(ik4), intent(in) :: isize, icpu1, icpu2, itag
    call mpi_sendrecv(rv1,isize,mpi_real4,icpu1,itag, &
                      rv2,isize,mpi_real4,icpu2,itag, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      write(stderr, *) 'At line :', itag
      call fatal(__FILE__,__LINE__,'mpi_sendrecv error.')
    end if
#endif
  end subroutine cyclic_exchange_array_r4

  subroutine set_nproc
    implicit none
    integer(ik4), dimension(2) :: cpus_per_dim
    logical, dimension(2) :: dim_period
    integer(ik4), dimension(2) :: isearch
    integer(ik4) :: imaxcpus, imax1, imax2, imiss
    integer(ik4) :: maximum_buffer_size
    data dim_period /.false.,.false./

    ma%bandflag    = (i_band == 1 .or. i_crm  == 1)
    ma%crmflag     = (i_crm  == 1)
    ma%top         = mpi_proc_null
    ma%bottom      = mpi_proc_null
    ma%left        = mpi_proc_null
    ma%right       = mpi_proc_null
    ma%bottomleft  = mpi_proc_null
    ma%bottomright = mpi_proc_null
    ma%topleft     = mpi_proc_null
    ma%topright    = mpi_proc_null

    ma%has_bdyright       = .false.
    ma%has_bdyleft        = .false.
    ma%has_bdytop         = .false.
    ma%has_bdybottom      = .false.
    ma%has_bdytopleft     = .false.
    ma%has_bdytopright    = .false.
    ma%has_bdybottomleft  = .false.
    ma%has_bdybottomright = .false.

    if ( nproc == 1 ) then
      cpus_per_dim(1) = 1
      cpus_per_dim(2) = 1
      jxp =  jx
      iyp =  iy

      global_dot_jstart = 1
      global_dot_istart = 1
      global_dot_jend = jx
      global_dot_iend = iy

      global_cross_jstart = 1
      global_cross_istart = 1
      if ( ma%bandflag ) then
        global_cross_jend = jx
      else
        global_cross_jend = jx-1
      end if
      if ( ma%crmflag ) then
        global_cross_iend = iy
      else
        global_cross_iend = iy-1
      end if

      if ( ma%crmflag ) then
        dim_period(1) = .true.
        dim_period(2) = .true.
        ma%left  = myid
        ma%right = myid
        ma%top  = myid
        ma%bottom = myid
        ma%bottomleft  = myid
        ma%bottomright = myid
        ma%topleft     = myid
        ma%topright    = myid
        ma%has_bdy = .false.
      else
        dim_period(1) = .true.
        ma%has_bdytop    = .true.
        ma%has_bdybottom = .true.
        if ( ma%bandflag ) then
          ma%left  = myid
          ma%right = myid
        else
          ma%has_bdyright       = .true.
          ma%has_bdyleft        = .true.
          ma%has_bdytopleft     = .true.
          ma%has_bdytopright    = .true.
          ma%has_bdybottomleft  = .true.
          ma%has_bdybottomright = .true.
        end if
        ma%has_bdy = .true.
      end if

      call mpi_cart_create(mycomm,2,cpus_per_dim,dim_period,lreorder, &
                           cartesian_communicator,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_cart_create error.')
      end if
#endif
      call mpi_comm_rank(cartesian_communicator,ccid,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_comm_rank error.')
      end if
#endif
      call mpi_cart_coords(cartesian_communicator,ccid,2,ma%location,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_cart_coords error.')
      end if
#endif
      ccio = iocpu

    else

      if ( ma%bandflag ) dim_period(1) = .true.
      if ( ma%crmflag  ) then
        dim_period(1) = .true.
        dim_period(2) = .true.
      end if
      if ( njxcpus > 0 .or. niycpus > 0 ) then
        ! Force just the number of CPUs in J direction
        if ( njxcpus > 0 .and. niycpus <= 0 ) then
          cpus_per_dim(1) = njxcpus
          cpus_per_dim(2) = nproc/njxcpus
        else if ( njxcpus <= 0 .and. niycpus > 0 ) then
          cpus_per_dim(2) = niycpus
          cpus_per_dim(1) = nproc/niycpus
        else
          cpus_per_dim(1) = njxcpus
          cpus_per_dim(2) = niycpus
        end if
        if ( cpus_per_dim(1) * cpus_per_dim(2) /= nproc ) then
          write (stderr,*) 'Requested ', cpus_per_dim(1), 'x', &
                           cpus_per_dim(2), ' CPUS'
          write (stderr,*) 'Available from MPI commandline ', nproc
          call fatal(__FILE__,__LINE__,'CPU/RCPU mismatch')
        end if
      else
        if ( nproc < 4 ) then
          cpus_per_dim(1) = nproc
          cpus_per_dim(2) = 1
        else if ( nproc >= 4 ) then
          cpus_per_dim(1) = (nint(sqrt(dble(nproc)))/2)*2
          if ( iy > int(1.5*dble(jx)) ) then
            cpus_per_dim(1) = cpus_per_dim(1) - 1
            do while ( mod(nproc,cpus_per_dim(1)) /= 0 )
              cpus_per_dim(1) = cpus_per_dim(1) - 1
            end do
          else if ( jx > int(1.5*dble(iy)) ) then
            cpus_per_dim(1) = cpus_per_dim(1) + 1
            do while ( mod(nproc,cpus_per_dim(1)) /= 0 )
              cpus_per_dim(1) = cpus_per_dim(1) + 1
            end do
          else
            do while ( mod(nproc,cpus_per_dim(1)) /= 0 )
              cpus_per_dim(1) = cpus_per_dim(1) + 1
            end do
          end if
          cpus_per_dim(2) = nproc/cpus_per_dim(1)
          imaxcpus = cpus_per_dim(1)*cpus_per_dim(2)
          if ( mod(nproc,imaxcpus) /= 0 ) then
            write(stderr,*) 'Work does not evenly divide.'
            write(stderr,*) 'I have calculated : '
            write(stderr,*) 'CPUS DIM1 = ', cpus_per_dim(1)
            write(stderr,*) 'CPUS DIM2 = ', cpus_per_dim(2)
            imax1 = ((jx/3)/2)*2
            imax2 = ((iy/3)/2)*2
            write(stderr,*) 'Suggested maximum number of CPUS jx: ', imax1
            write(stderr,*) 'Suggested maximum number of CPUS iy: ', imax2
            write(stderr,*) 'Closest number : ', imaxcpus
            call fatal(__FILE__,__LINE__,'CPU/WORK mismatch')
          end if
        end if
      end if

      call mpi_cart_create(mycomm,2,cpus_per_dim,dim_period,lreorder, &
                           cartesian_communicator,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_cart_create error.')
      end if
#endif
      call mpi_comm_rank(cartesian_communicator,ccid,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_comm_rank error.')
      end if
#endif
      call mpi_cart_coords(cartesian_communicator,ccid,2,ma%location,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_cart_coords error.')
      end if
#endif
      !call mpi_comm_split_type(cartesian_communicator,mpi_comm_type_shared, &
      !                         0, mpi_info_null,node_local_communicator, &
      !                         mpierr)
#ifdef DEBUG
      !if ( mpierr /= mpi_success ) then
      !  call fatal(__FILE__,__LINE__,'mpi_comm_split_type error.')
      !end if
#endif
      !call mpi_comm_rank(node_local_communicator,myidshm,mpierr)
#ifdef DEBUG
      !if ( mpierr /= mpi_success ) then
      !  call fatal(__FILE__,__LINE__,'mpi_comm_rank error.')
      !end if
#endif
      !call mpi_comm_size(node_local_communicator, nprocshm, mpierr)
#ifdef DEBUG
      !if ( mpierr /= mpi_success ) then
      !  call fatal(__FILE__,__LINE__,'mpi_comm_size error.')
      !end if
#endif

      if ( myid == iocpu ) ccio = ccid

      call bcast(ccio)

      call mpi_cart_shift(cartesian_communicator, 0, 1, &
        ma%left, ma%right, mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_cart_shift error.')
      end if
#endif
      call mpi_cart_shift(cartesian_communicator, 1, 1, &
        ma%bottom, ma%top, mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_cart_shift error.')
      end if
#endif

      ! Set coordinates in the grid for the other processors
      if ( ma%top /= mpi_proc_null .and. ma%right /= mpi_proc_null ) then
        isearch(1) = ma%location(1)+1
        isearch(2) = ma%location(2)+1
        call mpi_cart_rank(cartesian_communicator,isearch,ma%topright,mpierr)
#ifdef DEBUG
        if ( mpierr /= mpi_success ) then
          call fatal(__FILE__,__LINE__,'mpi_cart_rank error.')
        end if
#endif
      end if
      if ( ma%top /= mpi_proc_null .and. ma%left /= mpi_proc_null ) then
        isearch(1) = ma%location(1)-1
        isearch(2) = ma%location(2)+1
        call mpi_cart_rank(cartesian_communicator,isearch,ma%topleft,mpierr)
#ifdef DEBUG
        if ( mpierr /= mpi_success ) then
          call fatal(__FILE__,__LINE__,'mpi_cart_rank error.')
        end if
#endif
      end if
      if ( ma%bottom /= mpi_proc_null .and. ma%right /= mpi_proc_null ) then
        isearch(1) = ma%location(1)+1
        isearch(2) = ma%location(2)-1
        call mpi_cart_rank(cartesian_communicator,isearch,ma%bottomright,mpierr)
#ifdef DEBUG
        if ( mpierr /= mpi_success ) then
          call fatal(__FILE__,__LINE__,'mpi_cart_rank error.')
        end if
#endif
      end if
      if ( ma%bottom /= mpi_proc_null .and. ma%left /= mpi_proc_null ) then
        isearch(1) = ma%location(1)-1
        isearch(2) = ma%location(2)-1
        call mpi_cart_rank(cartesian_communicator,isearch,ma%bottomleft,mpierr)
#ifdef DEBUG
        if ( mpierr /= mpi_success ) then
          call fatal(__FILE__,__LINE__,'mpi_cart_rank error.')
        end if
#endif
      end if

      ma%has_bdytop    = (ma%top    == mpi_proc_null)
      ma%has_bdybottom = (ma%bottom == mpi_proc_null)
      ma%has_bdyright  = (ma%right  == mpi_proc_null)
      ma%has_bdyleft   = (ma%left   == mpi_proc_null)
      ma%has_bdy       = ( ma%has_bdytop .or. ma%has_bdybottom .or. &
                           ma%has_bdyright .or. ma%has_bdyleft )
      ma%has_bdytopleft     = (ma%has_bdytop .and. ma%has_bdyleft)
      ma%has_bdytopright    = (ma%has_bdytop .and. ma%has_bdyright)
      ma%has_bdybottomleft  = (ma%has_bdybottom .and. ma%has_bdyleft)
      ma%has_bdybottomright = (ma%has_bdybottom .and. ma%has_bdyright)

      jxp =  jx/cpus_per_dim(1)
      iyp =  iy/cpus_per_dim(2)

      global_dot_jstart = ma%location(1)*jxp+1
      global_dot_istart = ma%location(2)*iyp+1

      if ( jxp * cpus_per_dim(1) < jx ) then
        imiss = jx - jxp * cpus_per_dim(1)
        if ( ma%location(1) < imiss ) then
          global_dot_jstart = global_dot_jstart + ma%location(1)
          jxp = jxp + 1
        else
          global_dot_jstart = global_dot_jstart + imiss
        end if
      end if
      if ( iyp * cpus_per_dim(2) < iy ) then
        imiss = iy - iyp * cpus_per_dim(2)
        if ( ma%location(2) < imiss ) then
          global_dot_istart = global_dot_istart + ma%location(2)
          iyp = iyp + 1
        else
          global_dot_istart = global_dot_istart + imiss
        end if
      end if

      global_dot_jend = global_dot_jstart+jxp-1
      global_dot_iend = global_dot_istart+iyp-1
      if ( global_dot_iend > iy .or. global_dot_jend > jx ) then
        write(stderr,*) 'Cannot evenly divide!!!!'
        write(stderr,*) 'Processor ',myid,' has I : ', global_dot_istart, &
                                                       global_dot_iend
        write(stderr,*) 'Processor ',myid,' has J : ', global_dot_jstart, &
                                                       global_dot_jend
        call fatal(__FILE__,__LINE__,'DECOMPOSITION ERROR')
      end if

      global_cross_istart = global_dot_istart
      global_cross_iend = global_dot_iend
      if ( .not. ma%crmflag ) then
        if ( global_dot_iend == iy ) then
          ! South-North direction: Cross grid is one internal to the dot one.
          global_cross_iend = global_cross_iend - 1
        end if
      end if

      global_cross_jstart = global_dot_jstart
      global_cross_jend = global_dot_jend
      if ( .not. ma%bandflag ) then
        if ( global_dot_jend == jx ) then
          ! West-East direction: Cross grid is one internal to the dot one.
          global_cross_jend = global_cross_jend - 1
        end if
      end if
    end if
    !
    ! Check the results to be fit (minum for the advection is to have 3 points
    !
    if ( jxp < 3 .or. iyp < 3 ) then
      write(stderr,*) 'CPUS DIM1 = ', cpus_per_dim(1)
      write(stderr,*) 'CPUS DIM2 = ', cpus_per_dim(2)
      write(stderr,*) 'Cannot have one processor with less than 3x3 points.'
      write(stderr,*) 'Processor ',myid,' has ',jxp*iyp,' (',jxp,'x',iyp,')'
      call fatal(__FILE__,__LINE__,'Too much processors')
    end if

    jxpsg  = jxp * nsg
    iypsg  = iyp * nsg

    if ( myid == italk ) then
      write(stdout,*) 'CPUS DIM1 = ', cpus_per_dim(1)
      write(stdout,*) 'CPUS DIM2 = ', cpus_per_dim(2)
      write(stdout,*)
    end if
    if ( nproc > 1 ) then
      if ( myid == ccio ) then
        call getmem1d(windows,1,nproc*4,'set_nproc:windows')
      else
        windows => ifake
      end if
      call getmem1d(wincount,1,nproc*4,'set_nproc:wincount')
      call getmem1d(windispl,1,nproc*4,'set_nproc:windispl')
      ! Allocate to something should fit all
      maximum_buffer_size = jxsg*iysg
      maximum_buffer_size = max(maximum_buffer_size,jxpsg*iypsg*kzp1)
      call getmem1d(r8vector1,1,maximum_buffer_size,'set_nproc:r8vector1')
      call getmem1d(r8vector2,1,maximum_buffer_size,'set_nproc:r8vector2')
      call getmem1d(r4vector1,1,maximum_buffer_size,'set_nproc:r4vector1')
      call getmem1d(r4vector2,1,maximum_buffer_size,'set_nproc:r4vector2')
      call getmem1d(i4vector1,1,maximum_buffer_size,'set_nproc:i4vector1')
      call getmem1d(i4vector2,1,maximum_buffer_size,'set_nproc:i4vector2')
      call getmem1d(lvector1,1,maximum_buffer_size,'set_nproc:lvector1')
      call getmem1d(lvector2,1,maximum_buffer_size,'set_nproc:lvector2')
    end if
  end subroutine set_nproc

  subroutine broadcast_params
    implicit none

    call bcast(namelistfile,256)
    call bcast(prgname,256)

    call bcast(iy)
    call bcast(jx)
    call bcast(kz)
    call bcast(nsg)
    call bcast(njxcpus)
    call bcast(niycpus)
    call bcast(nveg)

    call bcast(idynamic)

    call bcast(iproj,6)
    call bcast(ds)
    call bcast(ptop)
    call bcast(clat)
    call bcast(clon)
    call bcast(cntri)
    call bcast(cntrj)
    call bcast(plat)
    call bcast(plon)
    call bcast(truelatl)
    call bcast(truelath)
    call bcast(i_band)
    call bcast(i_crm)

    call bcast(domname,64)

    call bcast(debug_level)
    call bcast(dbgfrq)
    call bcast(inpglob,256)

    call bcast(nspgx)
    call bcast(nspgd)
    call bcast(high_nudge)
    call bcast(medium_nudge)
    call bcast(low_nudge)
    call bcast(bdy_nm)
    call bcast(bdy_dm)

    call bcast(calendar,12)
    call bcast(ical)
    call bcast(dayspy)
    call bcast(vernal_equinox)

    call bcast(ibdyfrq)

    if ( idynamic == 3 ) then
      call bcast(mo_a0)
      call bcast(mo_h)
      call bcast(mo_ztop)
    end if

    ! Setup all convenience dimensions
    ! The IOCPU has performed this previously.

    if ( myid /= iocpu ) then
      iym1 = iy - 1
      iym2 = iy - 2
      iym3 = iy - 3
      jxm1 = jx - 1
      jxm2 = jx - 2
      jxm3 = jx - 3
      kzm1 = kz - 1
      kzm2 = kz - 2
      kzp1 = kz + 1
      kzp2 = kz + 2
      kzp3 = kz + 3
      kzp4 = kz + 4
      iysg = iy * nsg
      jxsg = jx * nsg
      iym1sg = iym1 * nsg
      jxm1sg = jxm1 * nsg
      iym2sg = iym2 * nsg
      jxm2sg = jxm2 * nsg
      iym3sg = iym3 * nsg
      jxm3sg = jxm3 * nsg
      nnsg = nsg*nsg
      jdot1 = 1
      jdot2 = jx
      jcross1 = 1
      if ( i_band == 1 ) then
        jcross2 = jx
        jout1 = 1
        jout2 = jx
        joutsg1 = 1
        joutsg2 = jx*nsg
      else
        jcross2 = jxm1
        jout1 = 2
        jout2 = jxm2
        joutsg1 = nsg+1
        joutsg2 = jxm2*nsg
      end if
      idot1 = 1
      idot2 = iy
      icross1 = 1
      if ( i_crm == 1) then
        icross2 = iy
        iout1 = 1
        iout2 = iy
        ioutsg1 = 1
        ioutsg2 = iy*nsg
      else
        icross2 = iym1
        iout1 = 2
        iout2 = iym2
        ioutsg1 = nsg+1
        ioutsg2 = iym2*nsg
      end if
      njcross = jcross2-jcross1+1
      nicross = icross2-icross1+1
      njdot = jdot2-jdot1+1
      nidot = idot2-idot1+1
      njout = jout2-jout1+1
      niout = iout2-iout1+1
      njoutsg = joutsg2-joutsg1+1
      nioutsg = ioutsg2-ioutsg1+1
      dpd = 360.0_rkx/dayspy
      half_dayspy = dayspy/2.0_rkx
      sixteenth_dayspy = dayspy/16.0_rkx
    end if
  end subroutine broadcast_params

  integer(ik4) function glosplitw(j1,j2,i1,i2,ls) result(tsize)
    implicit none
    integer(ik4), intent(in) :: j1, j2, i1, i2
    logical, intent(in), optional :: ls
    integer(ik4) :: isize, jsize, lsize, icpu, isub
    tsize = 0
    if ( nproc == 1 ) return
    isub = 1
    if ( present(ls) ) then
      if ( ls ) isub = nnsg
    end if
    isize = i2-i1+1
    jsize = j2-j1+1
    tsize = isize*jsize*isub
    window(1) = i1
    window(2) = window(1)+isize-1
    window(3) = j1
    window(4) = window(3)+jsize-1
    call mpi_gather(window,4,mpi_integer4, &
                    windows,4,mpi_integer4,ccio,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gather error.')
    end if
#endif
    if ( ccid == ccio ) then
      do icpu = 0, nproc-1
        window = windows(icpu*4+1:icpu*4+4)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        lsize = isize*jsize*isub
        wincount(icpu+1) = lsize
        windispl(icpu+1) = sum(wincount(1:icpu))
      end do
    end if
  end function glosplitw

  subroutine real8_2d_do_distribute(mg,ml,j1,j2,i1,i2,tsize)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:), intent(in) :: mg  ! model global
    real(rk8), pointer, contiguous, dimension(:,:), intent(inout) :: ml ! model local
    integer(ik4), intent(in) :: j1, j2, i1, i2, tsize
    integer(ik4) :: ib, i, j, icpu
    if ( nproc == 1 ) then
      do concurrent ( j = j1:j2, i = i1:i2 )
        ml(j,i) = mg(j,i)
      end do
      return
    end if
    if ( ccid == ccio ) then
      ib = 1
      do icpu = 0, nproc-1
        window = windows(icpu*4+1:icpu*4+4)
        do i = window(1), window(2)
          do j = window(3), window(4)
            r8vector1(ib) = mg(j,i)
            ib = ib + 1
          end do
        end do
      end do
    end if
    call mpi_scatterv(r8vector1,wincount,windispl,mpi_real8, &
                      r8vector2,tsize,mpi_real8,ccio,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_scatterv error.')
    end if
#endif
    !ml(j1:j2,i1:i2) = reshape(r8vector2,shape(ml(j1:j2,i1:i2)))
    ib = 1
    do i = i1, i2
      do j = j1, j2
        ml(j,i) = r8vector2(ib)
        ib = ib + 1
      end do
    end do
  end subroutine real8_2d_do_distribute

  subroutine real4_2d_do_distribute(mg,ml,j1,j2,i1,i2,tsize)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:), intent(in) :: mg  ! model global
    real(rk4), pointer, contiguous, dimension(:,:), intent(inout) :: ml ! model local
    integer(ik4), intent(in) :: j1, j2, i1, i2, tsize
    integer(ik4) :: ib, i, j, icpu
    if ( nproc == 1 ) then
      do concurrent ( j = j1:j2, i = i1:i2 )
        ml(j,i) = mg(j,i)
      end do
      return
    end if
    if ( ccid == ccio ) then
      ib = 1
      do icpu = 0, nproc-1
        window = windows(icpu*4+1:icpu*4+4)
        do i = window(1), window(2)
          do j = window(3), window(4)
            r4vector1(ib) = mg(j,i)
            ib = ib + 1
          end do
        end do
      end do
    end if
    call mpi_scatterv(r4vector1,wincount,windispl,mpi_real4, &
                      r4vector2,tsize,mpi_real4,ccio,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_scatterv error.')
    end if
#endif
    ib = 1
    do i = i1, i2
      do j = j1, j2
        ml(j,i) = r4vector2(ib)
        ib = ib + 1
      end do
    end do
  end subroutine real4_2d_do_distribute

  subroutine integer4_2d_do_distribute(mg,ml,j1,j2,i1,i2,tsize)
    implicit none
    integer(ik4), pointer, contiguous, dimension(:,:), intent(in) :: mg  ! model global
    integer(ik4), pointer, contiguous, dimension(:,:), intent(inout) :: ml ! model local
    integer(ik4), intent(in) :: j1, j2, i1, i2, tsize
    integer(ik4) :: ib, i, j, icpu
    if ( nproc == 1 ) then
      do concurrent ( j = j1:j2, i = i1:i2 )
        ml(j,i) = mg(j,i)
      end do
      return
    end if
    if ( ccid == ccio ) then
      ib = 1
      do icpu = 0, nproc-1
        window = windows(icpu*4+1:icpu*4+4)
        do i = window(1), window(2)
          do j = window(3), window(4)
            i4vector1(ib) = mg(j,i)
            ib = ib + 1
          end do
        end do
      end do
    end if
    call mpi_scatterv(i4vector1,wincount,windispl,mpi_integer4, &
                      i4vector2,tsize,mpi_integer4,ccio,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_scatterv error.')
    end if
#endif
    ib = 1
    do i = i1, i2
      do j = j1, j2
        ml(j,i) = i4vector2(ib)
        ib = ib + 1
      end do
    end do
  end subroutine integer4_2d_do_distribute

  subroutine logical_2d_do_distribute(mg,ml,j1,j2,i1,i2,tsize)
    implicit none
    logical, pointer, contiguous, dimension(:,:), intent(in) :: mg  ! model global
    logical, pointer, contiguous, dimension(:,:), intent(inout) :: ml ! model local
    integer(ik4), intent(in) :: j1, j2, i1, i2, tsize
    integer(ik4) :: ib, i, j, icpu
    if ( nproc == 1 ) then
      do concurrent ( j = j1:j2, i = i1:i2 )
        ml(j,i) = mg(j,i)
      end do
      return
    end if
    if ( ccid == ccio ) then
      ib = 1
      do icpu = 0, nproc-1
        window = windows(icpu*4+1:icpu*4+4)
        do i = window(1), window(2)
          do j = window(3), window(4)
            lvector1(ib) = mg(j,i)
            ib = ib + 1
          end do
        end do
      end do
    end if
    call mpi_scatterv(lvector1,wincount,windispl,mpi_logical, &
                      lvector2,tsize,mpi_logical,ccio,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_scatterv error.')
    end if
#endif
    ib = 1
    do i = i1, i2
      do j = j1, j2
        ml(j,i) = lvector2(ib)
        ib = ib + 1
      end do
    end do
  end subroutine logical_2d_do_distribute

  subroutine real8_2d_distribute(mg,ml,j1,j2,i1,i2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:), intent(in) :: mg  ! model global
    real(rk8), pointer, contiguous, dimension(:,:), intent(inout) :: ml ! model local
    integer(ik4), intent(in) :: j1, j2, i1, i2
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2)
    call real8_2d_do_distribute(mg,ml,j1,j2,i1,i2,tsize)
  end subroutine real8_2d_distribute

  subroutine real8_3d_distribute(mg,ml,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(in) :: mg  ! model global
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(inout) :: ml ! model local
    integer(ik4), intent(in) :: j1, j2, i1, i2, k1, k2
    real(rk8), pointer, contiguous, dimension(:,:) :: mg2 => null()
    real(rk8), pointer, contiguous, dimension(:,:) :: ml2 => null()
    integer(ik4) :: tsize, k
    tsize = glosplitw(j1,j2,i1,i2)
    do k = k1, k2
      if ( ccid == ccio ) call assignpnt(mg,mg2,k)
      call assignpnt(ml,ml2,k)
      call real8_2d_do_distribute(mg2,ml2,j1,j2,i1,i2,tsize)
    end do
  end subroutine real8_3d_distribute

  subroutine real8_4d_distribute(mg,ml,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:,:), intent(in) :: mg  ! model global
    real(rk8), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml ! model local
    integer(ik4), intent(in) :: j1, j2, i1, i2, k1, k2, n1, n2
    real(rk8), pointer, contiguous, dimension(:,:) :: mg2 => null()
    real(rk8), pointer, contiguous, dimension(:,:) :: ml2 => null()
    integer(ik4) :: tsize, k, n
    tsize = glosplitw(j1,j2,i1,i2)
    do n = n1, n2
      do k = k1, k2
        if ( ccid == ccio ) call assignpnt(mg,mg2,k,n)
        call assignpnt(ml,ml2,k,n)
        call real8_2d_do_distribute(mg2,ml2,j1,j2,i1,i2,tsize)
      end do
    end do
  end subroutine real8_4d_distribute

  subroutine real4_2d_distribute(mg,ml,j1,j2,i1,i2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:), intent(in) :: mg  ! model global
    real(rk4), pointer, contiguous, dimension(:,:), intent(inout) :: ml ! model local
    integer(ik4), intent(in) :: j1, j2, i1, i2
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2)
    call real4_2d_do_distribute(mg,ml,j1,j2,i1,i2,tsize)
  end subroutine real4_2d_distribute

  subroutine real4_3d_distribute(mg,ml,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(in) :: mg  ! model global
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(inout) :: ml ! model local
    integer(ik4), intent(in) :: j1, j2, i1, i2, k1, k2
    real(rk4), pointer, contiguous, dimension(:,:) :: mg2 => null( )
    real(rk4), pointer, contiguous, dimension(:,:) :: ml2 => null( )
    integer(ik4) :: tsize, k
    tsize = glosplitw(j1,j2,i1,i2)
    do k = k1, k2
      if ( ccid == ccio ) call assignpnt(mg,mg2,k)
      call assignpnt(ml,ml2,k)
      call real4_2d_do_distribute(mg2,ml2,j1,j2,i1,i2,tsize)
    end do
  end subroutine real4_3d_distribute

  subroutine real4_4d_distribute(mg,ml,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(in) :: mg  ! model global
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml ! model local
    integer(ik4), intent(in) :: j1, j2, i1, i2, k1, k2, n1, n2
    real(rk4), pointer, contiguous, dimension(:,:) :: mg2 => null( )
    real(rk4), pointer, contiguous, dimension(:,:) :: ml2 => null( )
    integer(ik4) :: tsize, k, n
    tsize = glosplitw(j1,j2,i1,i2)
    do n = n1, n2
      do k = k1, k2
        if ( ccid == ccio ) call assignpnt(mg,mg2,k,n)
        call assignpnt(ml,ml2,k,n)
        call real4_2d_do_distribute(mg2,ml2,j1,j2,i1,i2,tsize)
      end do
    end do
  end subroutine real4_4d_distribute

  subroutine integer_2d_distribute(mg,ml,j1,j2,i1,i2)
    implicit none
    integer(ik4), pointer, contiguous, dimension(:,:), intent(in) :: mg  ! model global
    integer(ik4), pointer, contiguous, dimension(:,:), intent(inout) :: ml ! model local
    integer(ik4), intent(in) :: j1, j2, i1, i2
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2)
    call integer4_2d_do_distribute(mg,ml,j1,j2,i1,i2,tsize)
  end subroutine integer_2d_distribute

  subroutine integer_3d_distribute(mg,ml,j1,j2,i1,i2,k1,k2)
    implicit none
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(in) :: mg  ! model global
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(inout) :: ml !model local
    integer(ik4), intent(in) :: j1, j2, i1, i2, k1, k2
    integer(ik4), pointer, contiguous, dimension(:,:) :: mg2 => null()
    integer(ik4), pointer, contiguous, dimension(:,:) :: ml2 => null()
    integer(ik4) :: tsize, k
    tsize = glosplitw(j1,j2,i1,i2)
    do k = k1, k2
      if ( ccid == ccio ) call assignpnt(mg,mg2,k)
      call assignpnt(ml,ml2,k)
      call integer4_2d_do_distribute(mg2,ml2,j1,j2,i1,i2,tsize)
    end do
  end subroutine integer_3d_distribute

  subroutine integer_4d_distribute(mg,ml,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    integer(ik4), pointer, contiguous, dimension(:,:,:,:), intent(in) :: mg  ! model glob
    integer(ik4), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml !model loc
    integer(ik4), intent(in) :: j1, j2, i1, i2, k1, k2, n1, n2
    integer(ik4), pointer, contiguous, dimension(:,:) :: mg2 => null( )
    integer(ik4), pointer, contiguous, dimension(:,:) :: ml2 => null( )
    integer(ik4) :: tsize, k, n
    tsize = glosplitw(j1,j2,i1,i2)
    do n = n1, n2
      do k = k1, k2
        if ( ccid == ccio ) call assignpnt(mg,mg2,k,n)
        call assignpnt(ml,ml2,k,n)
        call integer4_2d_do_distribute(mg2,ml2,j1,j2,i1,i2,tsize)
      end do
    end do
  end subroutine integer_4d_distribute

  subroutine logical_2d_distribute(mg,ml,j1,j2,i1,i2)
    implicit none
    logical, pointer, contiguous, dimension(:,:), intent(in) :: mg  ! model global
    logical, pointer, contiguous, dimension(:,:), intent(inout) :: ml ! model local
    integer(ik4), intent(in) :: j1, j2, i1, i2
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2)
    call logical_2d_do_distribute(mg,ml,j1,j2,i1,i2,tsize)
  end subroutine logical_2d_distribute

  subroutine logical_3d_distribute(mg,ml,j1,j2,i1,i2,k1,k2)
    implicit none
    logical, pointer, contiguous, dimension(:,:,:), intent(in) :: mg  ! model global
    logical, pointer, contiguous, dimension(:,:,:), intent(inout) :: ml !model local
    integer(ik4), intent(in) :: j1, j2, i1, i2, k1, k2
    logical, pointer, contiguous, dimension(:,:) :: mg2 => null()
    logical, pointer, contiguous, dimension(:,:) :: ml2 => null()
    integer(ik4) :: tsize, k
    tsize = glosplitw(j1,j2,i1,i2)
    do k = k1, k2
      if ( ccid == ccio ) call assignpnt(mg,mg2,k)
      call assignpnt(ml,ml2,k)
      call logical_2d_do_distribute(mg2,ml2,j1,j2,i1,i2,tsize)
    end do
  end subroutine logical_3d_distribute

  subroutine logical_4d_distribute(mg,ml,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    logical, pointer, contiguous, dimension(:,:,:,:), intent(in) :: mg  ! model glob
    logical, pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml !model loc
    integer(ik4), intent(in) :: j1, j2, i1, i2, k1, k2, n1, n2
    logical, pointer, contiguous, dimension(:,:) :: mg2 => null( )
    logical, pointer, contiguous, dimension(:,:) :: ml2 => null( )
    integer(ik4) :: tsize, k, n
    tsize = glosplitw(j1,j2,i1,i2)
    do n = n1, n2
      do k = k1, k2
        if ( ccid == ccio ) call assignpnt(mg,mg2,k,n)
        call assignpnt(ml,ml2,k,n)
        call logical_2d_do_distribute(mg2,ml2,j1,j2,i1,i2,tsize)
      end do
    end do
  end subroutine logical_4d_distribute

  subroutine real8_2d_do_sub_distribute(mg,ml,j1,j2,i1,i2,tsize,mask)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(in) :: mg  ! model global
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(inout) :: ml ! model local
    logical, pointer, contiguous, dimension(:,:,:), intent(in), optional :: mask
    integer(ik4), intent(in) :: j1, j2, i1, i2, tsize
    integer(ik4) :: ib, i, j, n, icpu
    if ( nproc == 1 ) then
      if ( present(mask) ) then
        do concurrent ( n = 1:nnsg, j = j1:j2, i = i1:i2 )
          if ( mask(n,j,i) ) ml(n,j,i) = mg(n,j,i)
        end do
      else
        do concurrent ( n = 1:nnsg, j = j1:j2, i = i1:i2 )
          ml(n,j,i) = mg(n,j,i)
        end do
      end if
      return
    end if
    if ( ccid == ccio ) then
      ib = 1
      do icpu = 0, nproc-1
        window = windows(icpu*4+1:icpu*4+4)
        do i = window(1), window(2)
          do j = window(3), window(4)
            do n = 1, nnsg
              r8vector1(ib) = mg(n,j,i)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
    call mpi_scatterv(r8vector1,wincount,windispl,mpi_real8, &
                      r8vector2,tsize,mpi_real8,ccio,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_scatterv error.')
    end if
#endif
    ib = 1
    if ( present(mask) ) then
      do i = i1, i2
        do j = j1, j2
          do n = 1, nnsg
            if ( mask(n,j,i) ) ml(n,j,i) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    else
      do i = i1, i2
        do j = j1, j2
          do n = 1, nnsg
            ml(n,j,i) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine real8_2d_do_sub_distribute

  subroutine real4_2d_do_sub_distribute(mg,ml,j1,j2,i1,i2,tsize,mask)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(in) :: mg  ! model global
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(inout) :: ml ! model local
    logical, pointer, contiguous, dimension(:,:,:), intent(in), optional :: mask
    integer(ik4), intent(in) :: j1, j2, i1, i2, tsize
    integer(ik4) :: ib, i, j, n, icpu
    if ( nproc == 1 ) then
      if ( present(mask) ) then
        do concurrent ( n = 1:nnsg, j = j1:j2, i = i1:i2 )
          if ( mask(n,j,i) ) ml(n,j,i) = mg(n,j,i)
        end do
      else
        do concurrent ( n = 1:nnsg, j = j1:j2, i = i1:i2 )
          ml(n,j,i) = mg(n,j,i)
        end do
      end if
      return
    end if
    if ( ccid == ccio ) then
      ib = 1
      do icpu = 0, nproc-1
        window = windows(icpu*4+1:icpu*4+4)
        do i = window(1), window(2)
          do j = window(3), window(4)
            do n = 1, nnsg
              r4vector1(ib) = mg(n,j,i)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
    call mpi_scatterv(r4vector1,wincount,windispl,mpi_real4, &
                      r4vector2,tsize,mpi_real4,ccio,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_scatterv error.')
    end if
#endif
    ib = 1
    if ( present(mask) ) then
      do i = i1, i2
        do j = j1, j2
          do n = 1, nnsg
            if ( mask(n,j,i) ) ml(n,j,i) = r4vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    else
      do i = i1, i2
        do j = j1, j2
          do n = 1, nnsg
            ml(n,j,i) = r4vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine real4_2d_do_sub_distribute

  subroutine integer4_2d_do_sub_distribute(mg,ml,j1,j2,i1,i2,tsize,mask)
    implicit none
    integer(rk4), pointer, contiguous, dimension(:,:,:), intent(in) :: mg  ! model glb
    integer(rk4), pointer, contiguous, dimension(:,:,:), intent(inout) :: ml ! model loc
    logical, pointer, contiguous, dimension(:,:,:), intent(in), optional :: mask
    integer(ik4), intent(in) :: j1, j2, i1, i2, tsize
    integer(ik4) :: ib, i, j, n, icpu
    if ( nproc == 1 ) then
      if ( present(mask) ) then
        do concurrent ( n = 1:nnsg, j = j1:j2, i = i1:i2 )
          if ( mask(n,j,i) ) ml(n,j,i) = mg(n,j,i)
        end do
      else
        do concurrent ( n = 1:nnsg, j = j1:j2, i = i1:i2 )
          ml(n,j,i) = mg(n,j,i)
        end do
      end if
      return
    end if
    if ( ccid == ccio ) then
      ib = 1
      do icpu = 0, nproc-1
        window = windows(icpu*4+1:icpu*4+4)
        do i = window(1), window(2)
          do j = window(3), window(4)
            do n = 1, nnsg
              i4vector1(ib) = mg(n,j,i)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
    call mpi_scatterv(i4vector1,wincount,windispl,mpi_integer4, &
                      i4vector2,tsize,mpi_integer4,ccio,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_scatterv error.')
    end if
#endif
    ib = 1
    if ( present(mask) ) then
      do i = i1, i2
        do j = j1, j2
          do n = 1, nnsg
            if ( mask(n,j,i) ) ml(n,j,i) = i4vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    else
      do i = i1, i2
        do j = j1, j2
          do n = 1, nnsg
            ml(n,j,i) = i4vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine integer4_2d_do_sub_distribute

  subroutine logical_2d_do_sub_distribute(mg,ml,j1,j2,i1,i2,tsize,mask)
    implicit none
    logical, pointer, contiguous, dimension(:,:,:), intent(in) :: mg  ! model glb
    logical, pointer, contiguous, dimension(:,:,:), intent(inout) :: ml ! model loc
    logical, pointer, contiguous, dimension(:,:,:), intent(in), optional :: mask
    integer(ik4), intent(in) :: j1, j2, i1, i2, tsize
    integer(ik4) :: ib, i, j, n, icpu
    if ( nproc == 1 ) then
      if ( present(mask) ) then
        do concurrent ( n = 1:nnsg, j = j1:j2, i = i1:i2 )
          if ( mask(n,j,i) ) ml(n,j,i) = mg(n,j,i)
        end do
      else
        do concurrent ( n = 1:nnsg, j = j1:j2, i = i1:i2 )
          ml(n,j,i) = mg(n,j,i)
        end do
      end if
      return
    end if
    if ( ccid == ccio ) then
      ib = 1
      do icpu = 0, nproc-1
        window = windows(icpu*4+1:icpu*4+4)
        do i = window(1), window(2)
          do j = window(3), window(4)
            do n = 1, nnsg
              lvector1(ib) = mg(n,j,i)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
    call mpi_scatterv(lvector1,wincount,windispl,mpi_logical, &
                      lvector2,tsize,mpi_logical,ccio,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_scatterv error.')
    end if
#endif
    ib = 1
    if ( present(mask) ) then
      do i = i1, i2
        do j = j1, j2
          do n = 1, nnsg
            if ( mask(n,j,i) ) ml(n,j,i) = lvector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    else
      do i = i1, i2
        do j = j1, j2
          do n = 1, nnsg
            ml(n,j,i) = lvector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine logical_2d_do_sub_distribute

  subroutine real8_2d_sub_distribute(mg,ml,j1,j2,i1,i2,mask)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(in) :: mg  ! model global
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(inout) :: ml ! model local
    logical, pointer, contiguous, dimension(:,:,:), intent(in), optional :: mask
    integer(ik4), intent(in) :: j1, j2, i1, i2
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2,.true.)
    call real8_2d_do_sub_distribute(mg,ml,j1,j2,i1,i2,tsize,mask)
  end subroutine real8_2d_sub_distribute

  subroutine real4_2d_sub_distribute(mg,ml,j1,j2,i1,i2,mask)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(in) :: mg  ! model global
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(inout) :: ml ! model local
    logical, pointer, contiguous, dimension(:,:,:), intent(in), optional :: mask
    integer(ik4), intent(in) :: j1, j2, i1, i2
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2,.true.)
    call real4_2d_do_sub_distribute(mg,ml,j1,j2,i1,i2,tsize,mask)
  end subroutine real4_2d_sub_distribute

  subroutine real8_3d_sub_distribute(mg,ml,j1,j2,i1,i2,k1,k2,mask)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:,:), intent(in) :: mg  ! model global
    real(rk8), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml ! model local
    logical, pointer, contiguous, dimension(:,:,:), intent(in), optional :: mask
    integer(ik4), intent(in) :: j1, j2, i1, i2, k1, k2
    real(rk8), pointer, contiguous, dimension(:,:,:) :: mg2 => null()
    real(rk8), pointer, contiguous, dimension(:,:,:) :: ml2 => null()
    integer(ik4) :: tsize, k
    tsize = glosplitw(j1,j2,i1,i2,.true.)
    do k = k1, k2
      if ( ccid == ccio ) call assignpnt(mg,mg2,k)
      call assignpnt(ml,ml2,k)
      call real8_2d_do_sub_distribute(mg2,ml2,j1,j2,i1,i2,tsize,mask)
    end do
  end subroutine real8_3d_sub_distribute

  subroutine real4_3d_sub_distribute(mg,ml,j1,j2,i1,i2,k1,k2,mask)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(in) :: mg  ! model global
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml ! model local
    logical, pointer, contiguous, dimension(:,:,:), intent(in), optional :: mask
    integer(ik4), intent(in) :: j1, j2, i1, i2, k1, k2
    real(rk4), pointer, contiguous, dimension(:,:,:) :: mg2 => null()
    real(rk4), pointer, contiguous, dimension(:,:,:) :: ml2 => null()
    integer(ik4) :: tsize, k
    tsize = glosplitw(j1,j2,i1,i2,.true.)
    do k = k1, k2
      if ( ccid == ccio ) call assignpnt(mg,mg2,k)
      call assignpnt(ml,ml2,k)
      call real4_2d_do_sub_distribute(mg2,ml2,j1,j2,i1,i2,tsize,mask)
    end do
  end subroutine real4_3d_sub_distribute

  subroutine logical_2d_sub_distribute(mg,ml,j1,j2,i1,i2,mask)
    implicit none
    logical, pointer, contiguous, dimension(:,:,:), intent(in) :: mg  ! model global
    logical, pointer, contiguous, dimension(:,:,:), intent(inout) :: ml ! model local
    logical, pointer, contiguous, dimension(:,:,:), intent(in), optional :: mask
    integer(ik4), intent(in) :: j1, j2, i1, i2
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2,.true.)
    call logical_2d_do_sub_distribute(mg,ml,j1,j2,i1,i2,tsize,mask)
  end subroutine logical_2d_sub_distribute

  subroutine integer_2d_sub_distribute(mg,ml,j1,j2,i1,i2,mask)
    implicit none
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(in) :: mg  ! model global
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(inout) :: ml ! model locl
    logical, pointer, contiguous, dimension(:,:,:), intent(in), optional :: mask
    integer(ik4), intent(in) :: j1, j2, i1, i2
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2,.true.)
    call integer4_2d_do_sub_distribute(mg,ml,j1,j2,i1,i2,tsize,mask)
  end subroutine integer_2d_sub_distribute

  subroutine integer_3d_sub_distribute(mg,ml,j1,j2,i1,i2,k1,k2,mask)
    implicit none
    integer(ik4), pointer, contiguous, dimension(:,:,:,:), intent(in) :: mg  ! model glob
    integer(ik4), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml ! modl loc
    logical, pointer, contiguous, dimension(:,:,:), intent(in), optional :: mask
    integer(ik4), intent(in) :: j1, j2, i1, i2, k1, k2
    integer(ik4), pointer, contiguous, dimension(:,:,:) :: mg2 => null()
    integer(ik4), pointer, contiguous, dimension(:,:,:) :: ml2 => null()
    integer(ik4) :: tsize, k
    tsize = glosplitw(j1,j2,i1,i2,.true.)
    do k = k1, k2
      if ( ccid == ccio ) call assignpnt(mg,mg2,k)
      call assignpnt(ml,ml2,k)
      call integer4_2d_do_sub_distribute(mg2,ml2,j1,j2,i1,i2,tsize,mask)
    end do
  end subroutine integer_3d_sub_distribute

  subroutine real8_2d_do_collect(ml,mg,j1,j2,i1,i2,tsize)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:), intent(in) :: ml  ! model local
    real(rk8), pointer, contiguous, dimension(:,:), intent(inout) :: mg ! model global
    integer(ik4), intent(in) :: j1, j2, i1, i2, tsize
    integer(ik4) :: ib, i, j, icpu
    if ( nproc == 1 ) then
      do concurrent ( j = j1:j2, i = i1:i2 )
        mg(j,i) = ml(j,i)
      end do
      return
    end if
    ib = 1
    do i = i1, i2
      do j = j1, j2
        r8vector2(ib) = ml(j,i)
        ib = ib + 1
      end do
    end do
    call mpi_gatherv(r8vector2,tsize,mpi_real8, &
                     r8vector1,wincount,windispl,mpi_real8, &
                     ccio,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
#endif
    if ( ccid == ccio ) then
      ib = 1
      do icpu = 0, nproc-1
        window = windows(icpu*4+1:icpu*4+4)
        do i = window(1), window(2)
          do j = window(3), window(4)
            mg(j,i) = r8vector1(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine real8_2d_do_collect

  subroutine real4_2d_do_collect(ml,mg,j1,j2,i1,i2,tsize)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:), intent(in) :: ml  ! model local
    real(rk4), pointer, contiguous, dimension(:,:), intent(inout) :: mg ! model global
    integer(ik4), intent(in) :: j1, j2, i1, i2, tsize
    integer(ik4) :: ib, i, j, icpu
    if ( nproc == 1 ) then
      do concurrent ( j = j1:j2, i = i1:i2 )
        mg(j,i) = ml(j,i)
      end do
      return
    end if
    ib = 1
    do i = i1, i2
      do j = j1, j2
        r4vector2(ib) = ml(j,i)
        ib = ib + 1
      end do
    end do
    call mpi_gatherv(r4vector2,tsize,mpi_real4, &
                     r4vector1,wincount,windispl,mpi_real4, &
                     ccio,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
#endif
    if ( ccid == ccio ) then
      ib = 1
      do icpu = 0, nproc-1
        window = windows(icpu*4+1:icpu*4+4)
        do i = window(1), window(2)
          do j = window(3), window(4)
            mg(j,i) = r4vector1(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine real4_2d_do_collect

  subroutine integer4_2d_do_collect(ml,mg,j1,j2,i1,i2,tsize)
    implicit none
    integer(ik4), pointer, contiguous, dimension(:,:), intent(in) :: ml  ! model local
    integer(ik4), pointer, contiguous, dimension(:,:), intent(inout) :: mg ! model global
    integer(ik4), intent(in) :: j1, j2, i1, i2, tsize
    integer(ik4) :: ib, i, j, icpu
    if ( nproc == 1 ) then
      do concurrent ( j = j1:j2, i = i1:i2 )
        mg(j,i) = ml(j,i)
      end do
      return
    end if
    ib = 1
    do i = i1, i2
      do j = j1, j2
        i4vector2(ib) = ml(j,i)
        ib = ib + 1
      end do
    end do
    call mpi_gatherv(i4vector2,tsize,mpi_integer4, &
                     i4vector1,wincount,windispl,mpi_integer4, &
                     ccio,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
#endif
    if ( ccid == ccio ) then
      ib = 1
      do icpu = 0, nproc-1
        window = windows(icpu*4+1:icpu*4+4)
        do i = window(1), window(2)
          do j = window(3), window(4)
            mg(j,i) = i4vector1(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine integer4_2d_do_collect

  subroutine logical_2d_do_collect(ml,mg,j1,j2,i1,i2,tsize)
    implicit none
    logical, pointer, contiguous, dimension(:,:), intent(in) :: ml  ! model local
    logical, pointer, contiguous, dimension(:,:), intent(inout) :: mg ! model global
    integer(ik4), intent(in) :: j1, j2, i1, i2, tsize
    integer(ik4) :: ib, i, j, icpu
    if ( nproc == 1 ) then
      do concurrent ( j = j1:j2, i = i1:i2 )
        mg(j,i) = ml(j,i)
      end do
      return
    end if
    ib = 1
    do i = i1, i2
      do j = j1, j2
        lvector2(ib) = ml(j,i)
        ib = ib + 1
      end do
    end do
    call mpi_gatherv(lvector2,tsize,mpi_logical, &
                     lvector1,wincount,windispl,mpi_logical, &
                     ccio,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
#endif
    if ( ccid == ccio ) then
      ib = 1
      do icpu = 0, nproc-1
        window = windows(icpu*4+1:icpu*4+4)
        do i = window(1), window(2)
          do j = window(3), window(4)
            mg(j,i) = lvector1(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine logical_2d_do_collect

  subroutine real8_2d_collect(ml,mg,j1,j2,i1,i2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:), intent(in) :: ml  ! model local
    real(rk8), pointer, contiguous, dimension(:,:), intent(inout) :: mg ! model global
    integer(ik4), intent(in) :: j1, j2, i1, i2
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2)
    call real8_2d_do_collect(ml,mg,j1,j2,i1,i2,tsize)
  end subroutine real8_2d_collect

  subroutine real8_2d_3d_collect(ml,mg,j1,j2,i1,i2,k)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:), intent(in) :: ml    ! model local
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(inout) :: mg ! model global
    integer(ik4), intent(in) :: j1, j2, i1, i2
    integer(ik4), intent(in), optional :: k
    real(rk8), pointer, contiguous, dimension(:,:) :: mg2 => null( )
    integer(ik4) :: kk
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2)
    kk = 1
    if ( present(k) ) kk = k
    if ( ccid == ccio )  call assignpnt(mg,mg2,kk)
    call real8_2d_do_collect(ml,mg2,j1,j2,i1,i2,tsize)
  end subroutine real8_2d_3d_collect

  subroutine real8_3d_collect(ml,mg,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(in) :: ml  ! model local
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(inout) :: mg ! model global
    integer(ik4), intent(in) :: j1, j2, i1, i2, k1, k2
    real(rk8), pointer, contiguous, dimension(:,:) :: ml2 => null( )
    real(rk8), pointer, contiguous, dimension(:,:) :: mg2 => null( )
    integer(ik4) :: tsize, k
    tsize = glosplitw(j1,j2,i1,i2)
    do k = k1, k2
      if ( ccid == ccio )  call assignpnt(mg,mg2,k)
      call assignpnt(ml,ml2,k)
      call real8_2d_do_collect(ml2,mg2,j1,j2,i1,i2,tsize)
    end do
  end subroutine real8_3d_collect

  subroutine real8_3d_2d_collect(ml,mg,j1,j2,i1,i2,k)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(in) :: ml  ! model local
    real(rk8), pointer, contiguous, dimension(:,:), intent(inout) :: mg   ! model global
    integer(ik4), intent(in) :: j1, j2, i1, i2, k
    real(rk8), pointer, contiguous, dimension(:,:) :: ml2 => null( )
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2)
    call assignpnt(ml,ml2,k)
    call real8_2d_do_collect(ml2,mg,j1,j2,i1,i2,tsize)
  end subroutine real8_3d_2d_collect

  subroutine real8_4d_collect(ml,mg,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:,:), intent(in) :: ml  ! model local
    real(rk8), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: mg ! model glob
    integer(ik4), intent(in) :: j1, j2, i1, i2, k1, k2, n1, n2
    real(rk8), pointer, contiguous, dimension(:,:) :: ml2 => null( )
    real(rk8), pointer, contiguous, dimension(:,:) :: mg2 => null( )
    integer(ik4) :: tsize, k, n
    tsize = glosplitw(j1,j2,i1,i2)
    do n = n1, n2
      do k = k1, k2
        if ( ccid == ccio ) call assignpnt(mg,mg2,k,n)
        call assignpnt(ml,ml2,k,n)
        call real8_2d_do_collect(ml2,mg2,j1,j2,i1,i2,tsize)
      end do
    end do
  end subroutine real8_4d_collect

  subroutine real8_4d_2d_collect(ml,mg,j1,j2,i1,i2,k,n)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:,:), intent(in) :: ml ! model local
    real(rk8), pointer, contiguous, dimension(:,:), intent(inout) :: mg ! model global
    integer(ik4), intent(in) :: j1, j2, i1, i2, k, n
    real(rk8), pointer, contiguous, dimension(:,:) :: ml2 => null( )
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2)
    call assignpnt(ml,ml2,k,n)
    call real8_2d_do_collect(ml2,mg,j1,j2,i1,i2,tsize)
  end subroutine real8_4d_2d_collect

  subroutine real4_2d_collect(ml,mg,j1,j2,i1,i2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:), intent(in) :: ml  ! model local
    real(rk4), pointer, contiguous, dimension(:,:), intent(inout) :: mg ! model global
    integer(ik4), intent(in) :: j1, j2, i1, i2
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2)
    call real4_2d_do_collect(ml,mg,j1,j2,i1,i2,tsize)
  end subroutine real4_2d_collect

  subroutine real4_2d_3d_collect(ml,mg,j1,j2,i1,i2,k)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:), intent(in) :: ml    ! model local
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(inout) :: mg ! model global
    integer(ik4), intent(in) :: j1, j2, i1, i2
    integer(ik4), intent(in), optional :: k
    real(rk4), pointer, contiguous, dimension(:,:) :: mg2 => null( )
    integer(ik4) :: kk
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2)
    kk = 1
    if ( present(k) ) kk = k
    if ( ccid == ccio )  call assignpnt(mg,mg2,kk)
    call real4_2d_do_collect(ml,mg2,j1,j2,i1,i2,tsize)
  end subroutine real4_2d_3d_collect

  subroutine real4_3d_collect(ml,mg,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(in) :: ml  ! model local
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(inout) :: mg ! model global
    integer(ik4), intent(in) :: j1, j2, i1, i2, k1, k2
    real(rk4), pointer, contiguous, dimension(:,:) :: ml2 => null( )
    real(rk4), pointer, contiguous, dimension(:,:) :: mg2 => null( )
    integer(ik4) :: tsize, k
    tsize = glosplitw(j1,j2,i1,i2)
    do k = k1, k2
      if ( ccid == ccio )  call assignpnt(mg,mg2,k)
      call assignpnt(ml,ml2,k)
      call real4_2d_do_collect(ml2,mg2,j1,j2,i1,i2,tsize)
    end do
  end subroutine real4_3d_collect

  subroutine real4_3d_2d_collect(ml,mg,j1,j2,i1,i2,k)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(in) :: ml  ! model local
    real(rk4), pointer, contiguous, dimension(:,:), intent(inout) :: mg   ! model global
    integer(ik4), intent(in) :: j1, j2, i1, i2, k
    real(rk4), pointer, contiguous, dimension(:,:) :: ml2 => null( )
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2)
    call assignpnt(ml,ml2,k)
    call real4_2d_do_collect(ml2,mg,j1,j2,i1,i2,tsize)
  end subroutine real4_3d_2d_collect

  subroutine real4_4d_collect(ml,mg,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(in) :: ml  ! model local
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: mg ! model glob
    integer(ik4), intent(in) :: j1, j2, i1, i2, k1, k2, n1, n2
    real(rk4), pointer, contiguous, dimension(:,:) :: ml2 => null( )
    real(rk4), pointer, contiguous, dimension(:,:) :: mg2 => null( )
    integer(ik4) :: tsize, k, n
    tsize = glosplitw(j1,j2,i1,i2)
    do n = n1, n2
      do k = k1, k2
        if ( ccid == ccio ) call assignpnt(mg,mg2,k,n)
        call assignpnt(ml,ml2,k,n)
        call real4_2d_do_collect(ml2,mg2,j1,j2,i1,i2,tsize)
      end do
    end do
  end subroutine real4_4d_collect

  subroutine real4_4d_2d_collect(ml,mg,j1,j2,i1,i2,k,n)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(in) :: ml ! model local
    real(rk4), pointer, contiguous, dimension(:,:), intent(inout) :: mg ! model global
    integer(ik4), intent(in) :: j1, j2, i1, i2, k, n
    real(rk4), pointer, contiguous, dimension(:,:) :: ml2 => null( )
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2)
    call assignpnt(ml,ml2,k,n)
    call real4_2d_do_collect(ml2,mg,j1,j2,i1,i2,tsize)
  end subroutine real4_4d_2d_collect

  subroutine logical_2d_collect(ml,mg,j1,j2,i1,i2)
    implicit none
    logical, pointer, contiguous, dimension(:,:), intent(in) :: ml  ! model local
    logical, pointer, contiguous, dimension(:,:), intent(inout) :: mg ! model global
    integer(ik4), intent(in) :: j1, j2, i1, i2
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2)
    call logical_2d_do_collect(ml,mg,j1,j2,i1,i2,tsize)
  end subroutine logical_2d_collect

  subroutine integer_2d_collect(ml,mg,j1,j2,i1,i2)
    implicit none
    integer(ik4), pointer, contiguous, dimension(:,:), intent(in) :: ml  ! model local
    integer(ik4), pointer, contiguous, dimension(:,:), intent(inout) :: mg ! model global
    integer(ik4), intent(in) :: j1, j2, i1, i2
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2)
    call integer4_2d_do_collect(ml,mg,j1,j2,i1,i2,tsize)
  end subroutine integer_2d_collect

  subroutine integer_3d_collect(ml,mg,j1,j2,i1,i2,k1,k2)
    implicit none
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(in) :: ml  ! model local
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(inout) :: mg ! model glbl
    integer(ik4), intent(in) :: j1, j2, i1, i2, k1, k2
    integer(ik4), pointer, contiguous, dimension(:,:) :: ml2 => null()
    integer(ik4), pointer, contiguous, dimension(:,:) :: mg2 => null()
    integer(ik4) :: tsize, k
    tsize = glosplitw(j1,j2,i1,i2)
    do k = k1, k2
      if ( ccid == ccio )  call assignpnt(mg,mg2,k)
      call assignpnt(ml,ml2,k)
      call integer4_2d_do_collect(ml2,mg2,j1,j2,i1,i2,tsize)
    end do
  end subroutine integer_3d_collect

  subroutine integer_4d_collect(ml,mg,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    integer(ik4), pointer, contiguous, dimension(:,:,:,:), intent(in) :: ml  ! mdl local
    integer(ik4), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: mg ! mdl glob
    integer(ik4), intent(in) :: j1, j2, i1, i2, k1, k2, n1, n2
    integer(ik4), pointer, contiguous, dimension(:,:) :: ml2 => null( )
    integer(ik4), pointer, contiguous, dimension(:,:) :: mg2 => null( )
    integer(ik4) :: tsize, k, n
    tsize = glosplitw(j1,j2,i1,i2)
    do n = n1, n2
      do k = k1, k2
        if ( ccid == ccio ) call assignpnt(mg,mg2,k,n)
        call assignpnt(ml,ml2,k,n)
        call integer4_2d_do_collect(ml2,mg2,j1,j2,i1,i2,tsize)
      end do
    end do
  end subroutine integer_4d_collect

  subroutine real8_2d_do_sub_collect(ml,mg,j1,j2,i1,i2,tsize)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(in) :: ml  ! model local
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(inout) :: mg ! model global
    integer(ik4), intent(in) :: j1, j2, i1, i2, tsize
    integer(ik4) :: ib, i, j, n, icpu
    if ( nproc == 1 ) then
      do concurrent ( n = 1:nnsg, j = j1:j2, i = i1:i2 )
        mg(n,j,i) = ml(n,j,i)
      end do
      return
    end if
    ib = 1
    do i = i1, i2
      do j = j1, j2
        do n = 1, nnsg
          r8vector2(ib) = ml(n,j,i)
          ib = ib + 1
        end do
      end do
    end do
    call mpi_gatherv(r8vector2,tsize,mpi_real8, &
                     r8vector1,wincount,windispl,mpi_real8, &
                     ccio,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
#endif
    if ( ccid == ccio ) then
      ib = 1
      do icpu = 0, nproc-1
        window = windows(icpu*4+1:icpu*4+4)
        do i = window(1), window(2)
          do j = window(3), window(4)
            do n = 1, nnsg
              mg(n,j,i) = r8vector1(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
  end subroutine real8_2d_do_sub_collect

  subroutine real4_2d_do_sub_collect(ml,mg,j1,j2,i1,i2,tsize)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(in) :: ml  ! model local
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(inout) :: mg ! model global
    integer(ik4), intent(in) :: j1, j2, i1, i2, tsize
    integer(ik4) :: ib, i, j, n, icpu
    if ( nproc == 1 ) then
      do concurrent ( n = 1:nnsg, j = j1:j2, i = i1:i2 )
        mg(n,j,i) = ml(n,j,i)
      end do
      return
    end if
    ib = 1
    do i = i1, i2
      do j = j1, j2
        do n = 1, nnsg
          r4vector2(ib) = ml(n,j,i)
          ib = ib + 1
        end do
      end do
    end do
    call mpi_gatherv(r4vector2,tsize,mpi_real4, &
                     r4vector1,wincount,windispl,mpi_real4, &
                     ccio,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
#endif
    if ( ccid == ccio ) then
      ib = 1
      do icpu = 0, nproc-1
        window = windows(icpu*4+1:icpu*4+4)
        do i = window(1), window(2)
          do j = window(3), window(4)
            do n = 1, nnsg
              mg(n,j,i) = r4vector1(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
  end subroutine real4_2d_do_sub_collect

  subroutine integer4_2d_do_sub_collect(ml,mg,j1,j2,i1,i2,tsize)
    implicit none
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(in) :: ml  ! model loc
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(inout) :: mg ! model glb
    integer(ik4), intent(in) :: j1, j2, i1, i2, tsize
    integer(ik4) :: ib, i, j, n, icpu
    if ( nproc == 1 ) then
      do concurrent ( n = 1:nnsg, j = j1:j2, i = i1:i2 )
        mg(n,j,i) = ml(n,j,i)
      end do
      return
    end if
    ib = 1
    do i = i1, i2
      do j = j1, j2
        do n = 1, nnsg
          i4vector2(ib) = ml(n,j,i)
          ib = ib + 1
        end do
      end do
    end do
    call mpi_gatherv(i4vector2,tsize,mpi_integer4, &
                     i4vector1,wincount,windispl,mpi_integer4, &
                     ccio,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
#endif
    if ( ccid == ccio ) then
      ib = 1
      do icpu = 0, nproc-1
        window = windows(icpu*4+1:icpu*4+4)
        do i = window(1), window(2)
          do j = window(3), window(4)
            do n = 1, nnsg
              mg(n,j,i) = i4vector1(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
  end subroutine integer4_2d_do_sub_collect

  subroutine logical_2d_do_sub_collect(ml,mg,j1,j2,i1,i2,tsize)
    implicit none
    logical, pointer, contiguous, dimension(:,:,:), intent(in) :: ml  ! model loc
    logical, pointer, contiguous, dimension(:,:,:), intent(inout) :: mg ! model glb
    integer(ik4), intent(in) :: j1, j2, i1, i2, tsize
    integer(ik4) :: ib, i, j, n, icpu
    if ( nproc == 1 ) then
      do concurrent ( n = 1:nnsg, j = j1:j2, i = i1:i2 )
        mg(n,j,i) = ml(n,j,i)
      end do
      return
    end if
    ib = 1
    do i = i1, i2
      do j = j1, j2
        do n = 1, nnsg
          lvector2(ib) = ml(n,j,i)
          ib = ib + 1
        end do
      end do
    end do
    call mpi_gatherv(lvector2,tsize,mpi_logical, &
                     lvector1,wincount,windispl,mpi_logical, &
                     ccio,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
#endif
    if ( ccid == ccio ) then
      ib = 1
      do icpu = 0, nproc-1
        window = windows(icpu*4+1:icpu*4+4)
        do i = window(1), window(2)
          do j = window(3), window(4)
            do n = 1, nnsg
              mg(n,j,i) = lvector1(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
  end subroutine logical_2d_do_sub_collect

  subroutine real8_2d_sub_collect(ml,mg,j1,j2,i1,i2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(in) :: ml  ! model local
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(inout) :: mg ! model global
    integer(ik4), intent(in) :: j1, j2, i1, i2
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2,.true.)
    call real8_2d_do_sub_collect(ml,mg,j1,j2,i1,i2,tsize)
  end subroutine real8_2d_sub_collect

  subroutine real8_3d_sub_collect(ml,mg,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:,:), intent(in) :: ml  ! model local
    real(rk8), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: mg ! model globl
    integer(ik4), intent(in) :: j1, j2, i1, i2, k1, k2
    real(rk8), pointer, contiguous, dimension(:,:,:) :: ml2 => null( )
    real(rk8), pointer, contiguous, dimension(:,:,:) :: mg2 => null( )
    integer(ik4) :: tsize, k
    tsize = glosplitw(j1,j2,i1,i2,.true.)
    do k = k1, k2
      if ( ccid == ccio ) call assignpnt(mg,mg2,k)
      call assignpnt(ml,ml2,k)
      call real8_2d_do_sub_collect(ml2,mg2,j1,j2,i1,i2,tsize)
    end do
  end subroutine real8_3d_sub_collect

  subroutine real4_2d_sub_collect(ml,mg,j1,j2,i1,i2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(in) :: ml  ! model local
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(inout) :: mg ! model global
    integer(ik4), intent(in) :: j1, j2, i1, i2
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2,.true.)
    call real4_2d_do_sub_collect(ml,mg,j1,j2,i1,i2,tsize)
  end subroutine real4_2d_sub_collect

  subroutine real4_3d_sub_collect(ml,mg,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(in) :: ml  ! model local
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: mg ! model globl
    integer(ik4), intent(in) :: j1, j2, i1, i2, k1, k2
    real(rk4), pointer, contiguous, dimension(:,:,:) :: ml2 => null( )
    real(rk4), pointer, contiguous, dimension(:,:,:) :: mg2 => null( )
    integer(ik4) :: tsize, k
    tsize = glosplitw(j1,j2,i1,i2,.true.)
    do k = k1, k2
      if ( ccid == ccio ) call assignpnt(mg,mg2,k)
      call assignpnt(ml,ml2,k)
      call real4_2d_do_sub_collect(ml2,mg2,j1,j2,i1,i2,tsize)
    end do
  end subroutine real4_3d_sub_collect

  subroutine integer_2d_sub_collect(ml,mg,j1,j2,i1,i2)
    implicit none
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(in) :: ml  ! model local
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(inout) :: mg ! model glob
    integer(ik4), intent(in) :: j1, j2, i1, i2
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2,.true.)
    call integer4_2d_do_sub_collect(ml,mg,j1,j2,i1,i2,tsize)
  end subroutine integer_2d_sub_collect

  subroutine integer_3d_sub_collect(ml,mg,j1,j2,i1,i2,k1,k2)
    implicit none
    integer(ik4), pointer, contiguous, dimension(:,:,:,:), intent(in) :: ml  ! mdl local
    integer(ik4), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: mg ! mdl glob
    integer(ik4), intent(in) :: j1, j2, i1, i2, k1, k2
    integer(ik4), pointer, contiguous, dimension(:,:,:) :: ml2 => null()
    integer(ik4), pointer, contiguous, dimension(:,:,:) :: mg2 => null()
    integer(ik4) :: tsize, k
    tsize = glosplitw(j1,j2,i1,i2,.true.)
    do k = k1, k2
      if ( ccid == ccio ) call assignpnt(mg,mg2,k)
      call assignpnt(ml,ml2,k)
      call integer4_2d_do_sub_collect(ml2,mg2,j1,j2,i1,i2,tsize)
    end do
  end subroutine integer_3d_sub_collect

  subroutine logical_2d_sub_collect(ml,mg,j1,j2,i1,i2)
    implicit none
    logical, pointer, contiguous, dimension(:,:,:), intent(in) :: ml  ! model local
    logical, pointer, contiguous, dimension(:,:,:), intent(inout) :: mg ! model glob
    integer(ik4), intent(in) :: j1, j2, i1, i2
    integer(ik4) :: tsize
    tsize = glosplitw(j1,j2,i1,i2,.true.)
    call logical_2d_do_sub_collect(ml,mg,j1,j2,i1,i2,tsize)
  end subroutine logical_2d_sub_collect

  subroutine logical_3d_sub_collect(ml,mg,j1,j2,i1,i2,k1,k2)
    implicit none
    logical, pointer, contiguous, dimension(:,:,:,:), intent(in) :: ml  ! mdl local
    logical, pointer, contiguous, dimension(:,:,:,:), intent(inout) :: mg ! mdl glob
    integer(ik4), intent(in) :: j1, j2, i1, i2, k1, k2
    logical, pointer, contiguous, dimension(:,:,:) :: ml2 => null()
    logical, pointer, contiguous, dimension(:,:,:) :: mg2 => null()
    integer(ik4) :: tsize, k
    tsize = glosplitw(j1,j2,i1,i2,.true.)
    do k = k1, k2
      if ( ccid == ccio ) call assignpnt(mg,mg2,k)
      call assignpnt(ml,ml2,k)
      call logical_2d_do_sub_collect(ml2,mg2,j1,j2,i1,i2,tsize)
    end do
  end subroutine logical_3d_sub_collect

#ifdef USE_MPI3
  subroutine real8_2d_exchange(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2
    integer(ik4) :: nx, ny
    integer(ik4) :: lb, rb
    integer(ik4) :: ndx, ndy, tx, ty, sizex, sizey

    nx = j2-j1+1
    ny = i2-i1+1
    lb = nex
    rb = nex
    if ( ma%left  == mpi_proc_null ) lb = 0
    if ( ma%right == mpi_proc_null ) rb = 0
    tx = ny
    ty = nx+lb+rb
    sizex = nex*tx
    sizey = nex*ty
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    real(rk8), dimension(ndx) :: sdatax
    real(rk8), dimension(ndy) :: sdatay
    real(rk8), dimension(ndx), volatile :: rdatax
    real(rk8), dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1, ib2, iex

    ib2 = 0
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2)
    end do
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2)
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdatax, counts, displs, mpi_real8, &
        rdatax, counts, displs, mpi_real8, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%left /= mpi_proc_null ) then
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        ml(j1-iex,i1:i2) = rdatax(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( ma%right /= mpi_proc_null ) then
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        ml(j2+iex,i1:i2) = rdatax(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizex
    end if

    ib2 = 0
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdatay(ib1:ib2) = ml(j1-lb:j2+rb,i1+(iex-1))
    end do
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdatay(ib1:ib2) = ml(j1-lb:j2+rb,i2-(iex-1))
    end do
    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdatay, counts, displs, mpi_real8, &
        rdatay, counts, displs, mpi_real8, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%bottom /= mpi_proc_null ) then
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        ml(j1-lb:j2+rb,i1-iex) = rdatay(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizey
    end if
    if ( ma%top /= mpi_proc_null ) then
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        ml(j1-lb:j2+rb,i2+iex) = rdatay(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange
  end subroutine real8_2d_exchange

  subroutine real8_3d_exchange(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2
    integer(ik4) :: nx, ny, nk
    integer(ik4) :: lb, rb
    integer(ik4) :: ndx, ndy, tx, ty, sizex, sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    lb = nex
    rb = nex
    if ( ma%left  == mpi_proc_null ) lb = 0
    if ( ma%right == mpi_proc_null ) rb = 0
    tx = ny
    ty = nx+lb+rb
    sizex = nex*tx*nk
    sizey = nex*ty*nk
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    real(rk8), dimension(ndx) :: sdatax
    real(rk8), dimension(ndy) :: sdatay
    real(rk8), dimension(ndx), volatile :: rdatax
    real(rk8), dimension(ndy), volatile :: rdatay
    integer(ik4) :: id, ib1, ib2, iex, k

    !$acc data copy(ml) create(sdatax,sdatay,rdatax,rdatay)

    do concurrent ( iex = 1:nex, k = k1:k2 )
      id = (k-k1)*nex + iex
      ib1 = 1 + (id-1)*(tx)
      ib2 = ib1 + tx - 1
      sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2,k)
    end do
    do concurrent ( iex = 1:nex, k = k1:k2 )
      id = (k-k1)*nex + iex
      ib1 = 1 + (id-1)*(tx) + sizex
      ib2 = ib1 + tx - 1
      sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2,k)
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    !$acc host_data use_device(sdatax,rdatax)
    call mpi_neighbor_alltoallv(sdatax, counts, displs, mpi_real8, &
        rdatax, counts, displs, mpi_real8, cartesian_communicator, mpierr)
    !$acc end host_data
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    if ( ma%left /= mpi_proc_null ) then
      do concurrent ( iex = 1:nex, k = k1:k2 )
        id = (k-k1)*nex + iex
        ib1 = 1 + (id-1)*(tx)
        ib2 = ib1 + tx - 1
        ml(j1-iex,i1:i2,k) = rdatax(ib1:ib2)
      end do
    end if
    if ( ma%right /= mpi_proc_null ) then
      do concurrent ( iex = 1:nex, k = k1:k2 )
        id = (k-k1)*nex + iex
        ib1 = 1 + (id-1)*(tx) + sizex
        ib2 = ib1 + tx - 1
        ml(j2+iex,i1:i2,k) = rdatax(ib1:ib2)
      end do
    end if

    do concurrent ( iex = 1:nex, k = k1:k2 )
      id = (k-k1)*nex + iex
      ib1 = 1 + (id-1)*(ty)
      ib2 = ib1 + ty - 1
      sdatay(ib1:ib2) = ml(j1-lb:j2+rb,i1+(iex-1),k)
    end do
    do concurrent ( iex = 1:nex, k = k1:k2 )
      id = (k-k1)*nex + iex
      ib1 = 1 + (id-1)*(ty) + sizey
      ib2 = ib1 + ty - 1
      sdatay(ib1:ib2) = ml(j1-lb:j2+rb,i2-(iex-1),k)
    end do
    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    !$acc host_data use_device(sdatay,rdatay)
    call mpi_neighbor_alltoallv(sdatay, counts, displs, mpi_real8, &
        rdatay, counts, displs, mpi_real8, cartesian_communicator, mpierr)
    !$acc end host_data
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    if ( ma%bottom /= mpi_proc_null ) then
      do concurrent ( iex = 1:nex, k = k1:k2 )
        id = (k-k1)*nex + iex
        ib1 = 1 + (id-1)*(ty)
        ib2 = ib1 + ty - 1
        ml(j1-lb:j2+rb,i1-iex,k) = rdatay(ib1:ib2)
      end do
    end if
    if ( ma%top /= mpi_proc_null ) then
      do concurrent ( iex = 1:nex, k = k1:k2 )
        id = (k-k1)*nex + iex
        ib1 = 1 + (id-1)*(ty) + sizey
        ib2 = ib1 + ty - 1
        ml(j1-lb:j2+rb,i2+iex,k) = rdatay(ib1:ib2)
      end do
    end if
    !$acc end data
    end block exchange

  end subroutine real8_3d_exchange

  subroutine real8_4d_exchange(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2, n1, n2
    integer(ik4) :: nx, ny, nk, nn
    integer(ik4) :: lb, rb
    integer(ik4) :: ndx, ndy, tx, ty, sizex, sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    nn = n2-n1+1
    lb = nex
    rb = nex
    if ( ma%left == mpi_proc_null ) lb = 0
    if ( ma%right == mpi_proc_null ) rb = 0
    tx = ny
    ty = nx+lb+rb
    sizex = nex*tx*nk*nn
    sizey = nex*ty*nk*nn
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    real(rk8), dimension(ndx) :: sdatax
    real(rk8), dimension(ndy) :: sdatay
    real(rk8), dimension(ndx), volatile :: rdatax
    real(rk8), dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1, ib2, iex, k, n

    ib2 = 0
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2,k,n)
        end do
      end do
    end do
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2,k,n)
        end do
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdatax, counts, displs, mpi_real8, &
        rdatax, counts, displs, mpi_real8, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%left /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + tx - 1
            ml(j1-iex,i1:i2,k,n) = rdatax(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( ma%right /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + tx - 1
            ml(j2+iex,i1:i2,k,n) = rdatax(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if

    ib2 = 0
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdatay(ib1:ib2) = ml(j1-lb:j2+rb,i1+(iex-1),k,n)
        end do
      end do
    end do
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdatay(ib1:ib2) = ml(j1-lb:j2+rb,i2-(iex-1),k,n)
        end do
      end do
    end do
    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdatay, counts, displs, mpi_real8, &
        rdatay, counts, displs, mpi_real8, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%bottom /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + ty - 1
            ml(j1-lb:j2+rb,i1-iex,k,n) = rdatay(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if
    if ( ma%top /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + ty - 1
            ml(j1-lb:j2+rb,i2+iex,k,n) = rdatay(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real8_4d_exchange

  subroutine real4_2d_exchange(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2
    integer(ik4) :: nx, ny
    integer(ik4) :: lb, rb
    integer(ik4) :: ndx, ndy, tx, ty, sizex, sizey

    nx = j2-j1+1
    ny = i2-i1+1
    lb = nex
    rb = nex
    if ( ma%left  == mpi_proc_null ) lb = 0
    if ( ma%right == mpi_proc_null ) rb = 0
    tx = ny
    ty = nx+lb+rb
    sizex = nex*tx
    sizey = nex*ty
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    real(rk4), dimension(ndx) :: sdatax
    real(rk4), dimension(ndy) :: sdatay
    real(rk4), dimension(ndx), volatile :: rdatax
    real(rk4), dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1, ib2, iex

    ib2 = 0
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2)
    end do
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2)
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdatax, counts, displs, mpi_real4, &
        rdatax, counts, displs, mpi_real4, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%left /= mpi_proc_null ) then
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        ml(j1-iex,i1:i2) = rdatax(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( ma%right /= mpi_proc_null ) then
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        ml(j2+iex,i1:i2) = rdatax(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizex
    end if

    ib2 = 0
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdatay(ib1:ib2) = ml(j1-lb:j2+rb,i1+(iex-1))
    end do
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdatay(ib1:ib2) = ml(j1-lb:j2+rb,i2-(iex-1))
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdatay, counts, displs, mpi_real4, &
        rdatay, counts, displs, mpi_real4, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%bottom /= mpi_proc_null ) then
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        ml(j1-lb:j2+rb,i1-iex) = rdatay(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizey
    end if
    if ( ma%top /= mpi_proc_null ) then
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        ml(j1-lb:j2+rb,i2+iex) = rdatay(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange
  end subroutine real4_2d_exchange

  subroutine real4_3d_exchange(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2
    integer(ik4) :: nx, ny, nk
    integer(ik4) :: lb, rb
    integer(ik4) :: ndx, ndy, tx, ty, sizex, sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    lb = nex
    rb = nex
    if ( ma%left  == mpi_proc_null ) lb = 0
    if ( ma%right == mpi_proc_null ) rb = 0
    tx = ny
    ty = nx+lb+rb
    sizex = nex*tx*nk
    sizey = nex*ty*nk
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    real(rk4), dimension(ndx) :: sdatax
    real(rk4), dimension(ndy) :: sdatay
    real(rk4), dimension(ndx), volatile :: rdatax
    real(rk4), dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1, ib2, iex, k

    ib2 = 0
    do k = k1, k2
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2,k)
      end do
    end do
    do k = k1, k2
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2,k)
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdatax, counts, displs, mpi_real4, &
        rdatax, counts, displs, mpi_real4, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_ineighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%left /= mpi_proc_null ) then
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          ml(j1-iex,i1:i2,k) = rdatax(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( ma%right /= mpi_proc_null ) then
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          ml(j2+iex,i1:i2,k) = rdatax(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if

    ib2 = 0
    do k = k1, k2
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        sdatay(ib1:ib2) = ml(j1-lb:j2+rb,i1+(iex-1),k)
      end do
    end do
    do k = k1, k2
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        sdatay(ib1:ib2) = ml(j1-lb:j2+rb,i2-(iex-1),k)
      end do
    end do
    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdatay, counts, displs, mpi_real8, &
        rdatay, counts, displs, mpi_real8, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%bottom /= mpi_proc_null ) then
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          ml(j1-lb:j2+rb,i1-iex,k) = rdatay(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if
    if ( ma%top /= mpi_proc_null ) then
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          ml(j1-lb:j2+rb,i2+iex,k) = rdatay(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real4_3d_exchange

  subroutine real4_4d_exchange(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2, n1, n2
    integer(ik4) :: nx, ny, nk, nn
    integer(ik4) :: lb, rb
    integer(ik4) :: ndx, ndy, tx, ty, sizex, sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    nn = n2-n1+1
    lb = nex
    rb = nex
    if ( ma%left == mpi_proc_null ) lb = 0
    if ( ma%right == mpi_proc_null ) rb = 0
    tx = ny
    ty = nx+lb+rb
    sizex = nex*tx*nk*nn
    sizey = nex*ty*nk*nn
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    real(rk4), dimension(ndx) :: sdatax
    real(rk4), dimension(ndy) :: sdatay
    real(rk4), dimension(ndx), volatile :: rdatax
    real(rk4), dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1, ib2, iex, k, n

    ib2 = 0
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2,k,n)
        end do
      end do
    end do
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2,k,n)
        end do
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdatax, counts, displs, mpi_real4, &
        rdatax, counts, displs, mpi_real4, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%left /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + tx - 1
            ml(j1-iex,i1:i2,k,n) = rdatax(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( ma%right /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + tx - 1
            ml(j2+iex,i1:i2,k,n) = rdatax(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if

    ib2 = 0
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdatay(ib1:ib2) = ml(j1-lb:j2+rb,i1+(iex-1),k,n)
        end do
      end do
    end do
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdatay(ib1:ib2) = ml(j1-lb:j2+rb,i2-(iex-1),k,n)
        end do
      end do
    end do
    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdatay, counts, displs, mpi_real4, &
        rdatay, counts, displs, mpi_real4, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%bottom /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + ty - 1
            ml(j1-lb:j2+rb,i1-iex,k,n) = rdatay(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if
    if ( ma%top /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + ty - 1
            ml(j1-lb:j2+rb,i2+iex,k,n) = rdatay(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real4_4d_exchange

  subroutine real8_2d_exchange_left_right_bottom_top(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2
    integer(ik4) :: ndx, ndy, nx, ny, tx, ty, sizex, sizey

    nx = j2-j1+1
    ny = i2-i1+1
    tx = ny
    ty = nx
    sizex = nex*tx
    sizey = nex*ty
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    real(rk8), dimension(ndx+ndy) :: sdata
    real(rk8), dimension(ndx+ndy), volatile :: rdata
    integer(ik4) :: ib1, ib2, iex

    ib2 = 0
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2)
    end do
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2)
    end do
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1))
    end do
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1))
    end do

    counts = [ sizex, sizex, sizey, sizey ]
    displs = [ 0, sizex, 2*sizex, 2*sizex+sizey ]
    call mpi_neighbor_alltoallv(sdata, counts, displs, mpi_real8, &
        rdata, counts, displs, mpi_real8, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%left /= mpi_proc_null ) then
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        ml(j1-iex,i1:i2) = rdata(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( ma%right /= mpi_proc_null ) then
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        ml(j2+iex,i1:i2) = rdata(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( ma%bottom /= mpi_proc_null ) then
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        ml(j1:j2,i1-iex) = rdata(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizey
    end if
    if ( ma%top /= mpi_proc_null ) then
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        ml(j1:j2,i2+iex) = rdata(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real8_2d_exchange_left_right_bottom_top

  subroutine real8_3d_exchange_left_right_bottom_top(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2
    integer(ik4) :: ndx, ndy, nx, ny, nk, tx, ty, sizex, sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    tx = ny
    ty = nx
    sizex = nex*tx*nk
    sizey = nex*ty*nk
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    real(rk8), dimension(ndx+ndy) :: sdata
    real(rk8), dimension(ndx+ndy), volatile :: rdata
    integer(ik4) :: id, ib1, ib2, iex, k

    !$acc data copy(ml) create(sdata,rdata)

    do concurrent ( iex = 1:nex, k = k1:k2 )
      id = (k-k1)*nex + iex
      ib1 = 1 + (id-1)*(tx)
      ib2 = ib1 + tx - 1
      sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2,k)
    end do
    do concurrent ( iex = 1:nex, k = k1:k2 )
      id = (k-k1)*nex + iex
      ib1 = 1 + (id-1)*(tx) + sizex
      ib2 = ib1 + tx - 1
      sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2,k)
    end do
    do concurrent ( iex = 1:nex, k = k1:k2 )
      id = (k-k1)*nex + iex
      ib1 = 1 + (id-1)*(ty) + sizex + sizex
      ib2 = ib1 + ty - 1
      sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1),k)
    end do
    do concurrent ( iex = 1:nex, k = k1:k2 )
      id = (k-k1)*nex + iex
      ib1 = 1 + (id-1)*(ty) + sizex + sizex + sizey
      ib2 = ib1 + ty - 1
      sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1),k)
    end do

    counts = [ sizex, sizex, sizey, sizey ]
    displs = [ 0, sizex, 2*sizex, 2*sizex+sizey ]
    !$acc host_data use_device(sdata,rdata)
    call mpi_neighbor_alltoallv(sdata, counts, displs, mpi_real8, &
        rdata, counts, displs, mpi_real8, cartesian_communicator, mpierr)
    !$acc end host_data
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    if ( ma%left /= mpi_proc_null ) then
      do concurrent ( iex = 1:nex, k = k1:k2 )
        id = (k-k1)*nex + iex
        ib1 = 1 + (id-1)*(tx)
        ib2 = ib1 + tx - 1
        ml(j1-iex,i1:i2,k) = rdata(ib1:ib2)
      end do
    end if
    if ( ma%right /= mpi_proc_null ) then
      do concurrent ( iex = 1:nex, k = k1:k2 )
        id = (k-k1)*nex + iex
        ib1 = 1 + (id-1)*(tx) + sizex
        ib2 = ib1 + tx - 1
        ml(j2+iex,i1:i2,k) = rdata(ib1:ib2)
      end do
    end if
    if ( ma%bottom /= mpi_proc_null ) then
      do concurrent ( iex = 1:nex, k = k1:k2 )
        id = (k-k1)*nex + iex
        ib1 = 1 + (id-1)*(ty) + sizex + sizex
        ib2 = ib1 + ty - 1
        ml(j1:j2,i1-iex,k) = rdata(ib1:ib2)
      end do
    end if
    if ( ma%top /= mpi_proc_null ) then
      do concurrent ( iex = 1:nex, k = k1:k2 )
        id = (k-k1)*nex + iex
        ib1 = 1 + (id-1)*(ty) + sizex + sizex + sizey
        ib2 = ib1 + ty - 1
        ml(j1:j2,i2+iex,k) = rdata(ib1:ib2)
      end do
    end if
    !$acc end data
    end block exchange

  end subroutine real8_3d_exchange_left_right_bottom_top

  subroutine real8_4d_exchange_left_right_bottom_top(ml,nex,j1,j2, &
                                                     i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2, n1, n2
    integer(ik4) :: ndx, ndy, nx, ny, nk, nn, tx, ty, sizex, sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    nn = n2-n1+1
    tx = ny
    ty = nx
    sizex = nex*tx*nk*nn
    sizey = nex*ty*nk*nn
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    real(rk8), dimension(ndx+ndy) :: sdata
    real(rk8), dimension(ndx+ndy), volatile :: rdata
    integer(ik4) :: ib1, ib2, iex, k, n

    ib2 = 0
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2,k,n)
        end do
      end do
    end do
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2,k,n)
        end do
      end do
    end do
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1),k,n)
        end do
      end do
    end do
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1),k,n)
        end do
      end do
    end do

    counts = [ sizex, sizex, sizey, sizey ]
    displs = [ 0, sizex, 2*sizex, 2*sizex+sizey ]
    call mpi_neighbor_alltoallv(sdata, counts, displs, mpi_real8, &
        rdata, counts, displs, mpi_real8, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%left /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + tx - 1
            ml(j1-iex,i1:i2,k,n) = rdata(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( ma%right /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + tx - 1
            ml(j2+iex,i1:i2,k,n) = rdata(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( ma%bottom /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + ty - 1
            ml(j1:j2,i1-iex,k,n) = rdata(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if
    if ( ma%top /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + ty - 1
            ml(j1:j2,i2+iex,k,n) = rdata(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real8_4d_exchange_left_right_bottom_top

  subroutine real4_2d_exchange_left_right_bottom_top(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2
    integer(ik4) :: ndx, ndy, nx, ny, tx, ty, sizex, sizey

    nx = j2-j1+1
    ny = i2-i1+1
    tx = ny
    ty = nx
    sizex = nex*tx
    sizey = nex*ty
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    real(rk4), dimension(ndx+ndy) :: sdata
    real(rk4), dimension(ndx+ndy), volatile :: rdata
    integer(ik4) :: ib1, ib2, iex

    ib2 = 0
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2)
    end do
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2)
    end do
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1))
    end do
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1))
    end do

    counts = [ sizex, sizex, sizey, sizey ]
    displs = [ 0, sizex, 2*sizex, 2*sizex+sizey ]
    call mpi_neighbor_alltoallv(sdata, counts, displs, mpi_real4, &
        rdata, counts, displs, mpi_real4, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%left /= mpi_proc_null ) then
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        ml(j1-iex,i1:i2) = rdata(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( ma%right /= mpi_proc_null ) then
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        ml(j2+iex,i1:i2) = rdata(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( ma%bottom /= mpi_proc_null ) then
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        ml(j1:j2,i1-iex) = rdata(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizey
    end if
    if ( ma%top /= mpi_proc_null ) then
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        ml(j1:j2,i2+iex) = rdata(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real4_2d_exchange_left_right_bottom_top

  subroutine real4_3d_exchange_left_right_bottom_top(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2
    integer(ik4) :: ndx, ndy, nx, ny, nk, tx, ty, sizex, sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    tx = ny
    ty = nx
    sizex = nex*tx*nk
    sizey = nex*ty*nk
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    real(rk4), dimension(ndx+ndy) :: sdata
    real(rk4), dimension(ndx+ndy), volatile :: rdata
    integer(ik4) :: ib1, ib2, iex, k

    ib2 = 0
    do k = k1, k2
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2,k)
      end do
    end do
    do k = k1, k2
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2,k)
      end do
    end do
    do k = k1, k2
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1),k)
      end do
    end do
    do k = k1, k2
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1),k)
      end do
    end do

    counts = [ sizex, sizex, sizey, sizey ]
    displs = [ 0, sizex, 2*sizex, 2*sizex+sizey ]
    call mpi_neighbor_alltoallv(sdata, counts, displs, mpi_real4, &
        rdata, counts, displs, mpi_real4, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%left /= mpi_proc_null ) then
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          ml(j1-iex,i1:i2,k) = rdata(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( ma%right /= mpi_proc_null ) then
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          ml(j2+iex,i1:i2,k) = rdata(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( ma%bottom /= mpi_proc_null ) then
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          ml(j1:j2,i1-iex,k) = rdata(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if
    if ( ma%top /= mpi_proc_null ) then
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          ml(j1:j2,i2+iex,k) = rdata(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real4_3d_exchange_left_right_bottom_top

  subroutine real4_4d_exchange_left_right_bottom_top(ml,nex,j1,j2, &
                                                     i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2, n1, n2
    integer(ik4) :: ndx, ndy, nx, ny, nk, nn, tx, ty, sizex, sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    nn = n2-n1+1
    tx = ny
    ty = nx
    sizex = nex*tx*nk*nn
    sizey = nex*ty*nk*nn
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    real(rk4), dimension(ndx+ndy) :: sdata
    real(rk4), dimension(ndx+ndy), volatile :: rdata
    integer(ik4) :: ib1, ib2, iex, k, n

    ib2 = 0
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2,k,n)
        end do
      end do
    end do
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2,k,n)
        end do
      end do
    end do
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1),k,n)
        end do
      end do
    end do
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1),k,n)
        end do
      end do
    end do

    counts = [ sizex, sizex, sizey, sizey ]
    displs = [ 0, sizex, 2*sizex, 2*sizex+sizey ]
    call mpi_neighbor_alltoallv(sdata, counts, displs, mpi_real4, &
        rdata, counts, displs, mpi_real4, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%left /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + tx - 1
            ml(j1-iex,i1:i2,k,n) = rdata(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( ma%right /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + tx - 1
            ml(j2+iex,i1:i2,k,n) = rdata(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( ma%bottom /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + ty - 1
            ml(j1:j2,i1-iex,k,n) = rdata(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if
    if ( ma%top /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + ty - 1
            ml(j1:j2,i2+iex,k,n) = rdata(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real4_4d_exchange_left_right_bottom_top

  subroutine real8_2d_exchange_left_right(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2
    integer(ik4) :: ndx, ny, tx, sizex

    ny = i2-i1+1
    tx = ny
    sizex = nex*tx
    ndx = 2*sizex

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    real(rk8), dimension(ndx) :: sdata
    real(rk8), dimension(ndx), volatile :: rdata
    integer(ik4) :: ib1, ib2, iex

    ib2 = 0
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2)
    end do
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2)
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdata, counts, displs, mpi_real8, &
        rdata, counts, displs, mpi_real8, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%left /= mpi_proc_null ) then
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        ml(j1-iex,i1:i2) = rdata(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( ma%right /= mpi_proc_null ) then
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        ml(j2+iex,i1:i2) = rdata(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizex
    end if

    end block exchange

  end subroutine real8_2d_exchange_left_right

  subroutine real8_3d_exchange_left_right(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2
    integer(ik4) :: ndx, ny, nk, tx, sizex

    ny = i2-i1+1
    nk = k2-k1+1
    tx = ny
    sizex = nex*tx*nk
    ndx = 2*sizex

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    real(rk8), dimension(ndx) :: sdata
    real(rk8), dimension(ndx), volatile :: rdata
    integer(ik4) :: id, ib1, ib2, iex, k

    !$acc data copy(ml) create(sdata,rdata)

    do concurrent ( iex = 1:nex, k = k1:k2 )
      id = (k-k1)*nex + iex
      ib1 = 1 + (id-1)*(tx)
      ib2 = ib1 + tx - 1
      sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2,k)
    end do
    do concurrent ( iex = 1:nex, k = k1:k2 )
      id = (k-k1)*nex + iex
      ib1 = 1 + (id-1)*(tx) + sizex
      ib2 = ib1 + tx - 1
      sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2,k)
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    !$acc host_data use_device(sdata,rdata)
    call mpi_neighbor_alltoallv(sdata, counts, displs, mpi_real8, &
        rdata, counts, displs, mpi_real8, cartesian_communicator, mpierr)
    !$acc end host_data
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    if ( ma%left /= mpi_proc_null ) then
      do concurrent ( iex = 1:nex, k = k1:k2 )
        id = (k-k1)*nex + iex
        ib1 = 1 + (id-1)*(tx)
        ib2 = ib1 + tx - 1
        ml(j1-iex,i1:i2,k) = rdata(ib1:ib2)
      end do
    end if
    if ( ma%right /= mpi_proc_null ) then
      do concurrent ( iex = 1:nex, k = k1:k2 )
        id = (k-k1)*nex + iex
        ib1 = 1 + (id-1)*(tx) + sizex
        ib2 = ib1 + tx - 1
        ml(j2+iex,i1:i2,k) = rdata(ib1:ib2)
      end do
    end if
    !$acc end data
    end block exchange

  end subroutine real8_3d_exchange_left_right

  subroutine real8_4d_exchange_left_right(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2, n1, n2
    integer(ik4) :: ndx, ny, nk, nn, tx, sizex

    ny = i2-i1+1
    nk = k2-k1+1
    nn = n2-n1+1
    tx = ny
    sizex = nex*tx*nk*nn
    ndx = 2*sizex

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    real(rk8), dimension(ndx) :: sdata
    real(rk8), dimension(ndx), volatile :: rdata
    integer(ik4) :: ib1, ib2, iex, k, n

    ib2 = 0
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2,k,n)
        end do
      end do
    end do
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2,k,n)
        end do
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdata, counts, displs, mpi_real8, &
        rdata, counts, displs, mpi_real8, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%left /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + tx - 1
            ml(j1-iex,i1:i2,k,n) = rdata(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( ma%right /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + tx - 1
            ml(j2+iex,i1:i2,k,n) = rdata(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if

    end block exchange

  end subroutine real8_4d_exchange_left_right

  subroutine real4_2d_exchange_left_right(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2
    integer(ik4) :: ndx, ny, tx, sizex

    ny = i2-i1+1
    tx = ny
    sizex = nex*tx
    ndx = 2*sizex

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    real(rk4), dimension(ndx) :: sdata
    real(rk4), dimension(ndx), volatile :: rdata
    integer(ik4) :: ib1, ib2, iex

    ib2 = 0
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2)
    end do
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2)
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdata, counts, displs, mpi_real4, &
        rdata, counts, displs, mpi_real4, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%left /= mpi_proc_null ) then
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        ml(j1-iex,i1:i2) = rdata(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( ma%right /= mpi_proc_null ) then
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        ml(j2+iex,i1:i2) = rdata(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizex
    end if

    end block exchange

  end subroutine real4_2d_exchange_left_right

  subroutine real4_3d_exchange_left_right(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2
    integer(ik4) :: ndx, ny, nk, tx, sizex

    ny = i2-i1+1
    nk = k2-k1+1
    tx = ny
    sizex = nex*tx*nk
    ndx = 2*sizex

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    real(rk4), dimension(ndx) :: sdata
    real(rk4), dimension(ndx), volatile :: rdata
    integer(ik4) :: ib1, ib2, iex, k

    ib2 = 0
    do k = k1, k2
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2,k)
      end do
    end do
    do k = k1, k2
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2,k)
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdata, counts, displs, mpi_real4, &
        rdata, counts, displs, mpi_real4, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%left /= mpi_proc_null ) then
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          ml(j1-iex,i1:i2,k) = rdata(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( ma%right /= mpi_proc_null ) then
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          ml(j2+iex,i1:i2,k) = rdata(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if

    end block exchange

  end subroutine real4_3d_exchange_left_right

  subroutine real4_4d_exchange_left_right(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2, n1, n2
    integer(ik4) :: ndx, ny, nk, nn, tx, sizex

    ny = i2-i1+1
    nk = k2-k1+1
    nn = n2-n1+1
    tx = ny
    sizex = nex*tx*nk*nn
    ndx = 2*sizex

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    real(rk4), dimension(ndx) :: sdata
    real(rk4), dimension(ndx), volatile :: rdata
    integer(ik4) :: ib1, ib2, iex, k, n

    ib2 = 0
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2,k,n)
        end do
      end do
    end do
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2,k,n)
        end do
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdata, counts, displs, mpi_real4, &
        rdata, counts, displs, mpi_real4, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%left /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + tx - 1
            ml(j1-iex,i1:i2,k,n) = rdata(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if
    if ( ma%right /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + tx - 1
            ml(j2+iex,i1:i2,k,n) = rdata(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if

    end block exchange

  end subroutine real4_4d_exchange_left_right

  subroutine real8_2d_exchange_bottom_top(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2
    integer(ik4) :: ndy, nx, ty, sizey

    nx = j2-j1+1
    ty = nx
    sizey = nex*ty
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    real(rk8), dimension(ndy) :: sdata
    real(rk8), dimension(ndy), volatile :: rdata
    integer(ik4) :: ib1, ib2, iex

    ib2 = 0
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1))
    end do
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1))
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdata, counts, displs, mpi_real8, &
        rdata, counts, displs, mpi_real8, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%bottom /= mpi_proc_null ) then
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        ml(j1:j2,i1-iex) = rdata(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizey
    end if
    if ( ma%top /= mpi_proc_null ) then
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        ml(j1:j2,i2+iex) = rdata(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real8_2d_exchange_bottom_top

  subroutine real8_3d_exchange_bottom_top(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2
    integer(ik4) :: ndy, nx, nk, ty, sizey

    nx = j2-j1+1
    nk = k2-k1+1
    ty = nx
    sizey = nex*ty*nk
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    real(rk8), dimension(ndy) :: sdata
    real(rk8), dimension(ndy), volatile :: rdata
    integer(ik4) :: id, ib1, ib2, iex, k

    !$acc data copy(ml) create(sdata,rdata)

    do concurrent ( iex = 1:nex, k = k1:k2 )
      id = (k-k1)*nex + iex
      ib1 = 1 + (id-1)*(ty)
      ib2 = ib1 + ty - 1
      sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1),k)
    end do
    do concurrent ( iex = 1:nex, k = k1:k2 )
      id = (k-k1)*nex + iex
      ib1 = 1 + (id-1)*(ty) + sizey
      ib2 = ib1 + ty - 1
      sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1),k)
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    !$acc host_data use_device(sdata,rdata)
    call mpi_neighbor_alltoallv(sdata, counts, displs, mpi_real8, &
        rdata, counts, displs, mpi_real8, cartesian_communicator, mpierr)
    !$acc end host_data
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    if ( ma%bottom /= mpi_proc_null ) then
      do concurrent ( iex = 1:nex, k = k1:k2 )
        id = (k-k1)*nex + iex
        ib1 = 1 + (id-1)*(ty)
        ib2 = ib1 + ty - 1
        ml(j1:j2,i1-iex,k) = rdata(ib1:ib2)
      end do
    end if
    if ( ma%top /= mpi_proc_null ) then
      do concurrent ( iex = 1:nex, k = k1:k2 )
        id = (k-k1)*nex + iex
        ib1 = 1 + (id-1)*(ty) + sizey
        ib2 = ib1 + ty - 1
        ml(j1:j2,i2+iex,k) = rdata(ib1:ib2)
      end do
    end if
    !$acc end data
    end block exchange

  end subroutine real8_3d_exchange_bottom_top

  subroutine real8_4d_exchange_bottom_top(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2, n1, n2
    integer(ik4) :: ndy, nx, nk, nn, ty, sizey

    nx = j2-j1+1
    nk = k2-k1+1
    nn = n2-n1+1
    ty = nx
    sizey = nex*ty*nk*nn
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    real(rk8), dimension(ndy) :: sdata
    real(rk8), dimension(ndy), volatile :: rdata
    integer(ik4) :: ib1, ib2, iex, k, n

    ib2 = 0
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1),k,n)
        end do
      end do
    end do
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1),k,n)
        end do
      end do
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdata, counts, displs, mpi_real8, &
        rdata, counts, displs, mpi_real8, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%bottom /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + ty - 1
            ml(j1:j2,i1-iex,k,n) = rdata(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if
    if ( ma%top /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + ty - 1
            ml(j1:j2,i2+iex,k,n) = rdata(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real8_4d_exchange_bottom_top

  subroutine real4_2d_exchange_bottom_top(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2
    integer(ik4) :: ndy, nx, ty, sizey

    nx = j2-j1+1
    ty = nx
    sizey = nex*ty
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    real(rk4), dimension(ndy) :: sdata
    real(rk4), dimension(ndy), volatile :: rdata
    integer(ik4) :: ib1, ib2, iex

    ib2 = 0
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1))
    end do
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1))
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdata, counts, displs, mpi_real4, &
        rdata, counts, displs, mpi_real4, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%bottom /= mpi_proc_null ) then
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        ml(j1:j2,i1-iex) = rdata(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizey
    end if
    if ( ma%top /= mpi_proc_null ) then
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        ml(j1:j2,i2+iex) = rdata(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real4_2d_exchange_bottom_top

  subroutine real4_3d_exchange_bottom_top(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2
    integer(ik4) :: ndy, nx, nk, ty, sizey

    nx = j2-j1+1
    nk = k2-k1+1
    ty = nx
    sizey = nex*ty*nk
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    real(rk4), dimension(ndy) :: sdata
    real(rk4), dimension(ndy), volatile :: rdata
    integer(ik4) :: id, ib1, ib2, iex, k

    !$acc data copy(ml) create(sdata,rdata)

    do concurrent ( iex = 1:nex, k = k1:k2 )
      id = (k-k1)*nex + iex
      ib1 = 1 + (id-1)*(ty)
      ib2 = ib1 + ty - 1
      sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1),k)
    end do
    do concurrent ( iex = 1:nex, k = k1:k2 )
      id = (k-k1)*nex + iex
      ib1 = 1 + (id-1)*(ty) + sizey
      ib2 = ib1 + ty - 1
      sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1),k)
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    !$acc host_data use_device(sdata,rdata)
    call mpi_neighbor_alltoallv(sdata, counts, displs, mpi_real4, &
        rdata, counts, displs, mpi_real4, cartesian_communicator, mpierr)
    !$acc end host_data
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    if ( ma%bottom /= mpi_proc_null ) then
      do concurrent ( iex = 1:nex, k = k1:k2 )
        id = (k-k1)*nex + iex
        ib1 = 1 + (id-1)*(ty)
        ib2 = ib1 + ty - 1
        ml(j1:j2,i1-iex,k) = rdata(ib1:ib2)
      end do
    end if
    if ( ma%top /= mpi_proc_null ) then
      do concurrent ( iex = 1:nex, k = k1:k2 )
        id = (k-k1)*nex + iex
        ib1 = 1 + (id-1)*(ty) + sizey
        ib2 = ib1 + ty - 1
        ml(j1:j2,i2+iex,k) = rdata(ib1:ib2)
      end do
    end if
    !$acc end data
    end block exchange

  end subroutine real4_3d_exchange_bottom_top

  subroutine real4_4d_exchange_bottom_top(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2, n1, n2
    integer(ik4) :: ndy, nx, nk, nn, ty, sizey

    nx = j2-j1+1
    nk = k2-k1+1
    nn = n2-n1+1
    ty = nx
    sizey = nex*ty*nk*nn
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    real(rk4), dimension(ndy) :: sdata
    real(rk4), dimension(ndy), volatile :: rdata
    integer(ik4) :: ib1, ib2, iex, k, n

    ib2 = 0
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1),k,n)
        end do
      end do
    end do
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1),k,n)
        end do
      end do
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdata, counts, displs, mpi_real4, &
        rdata, counts, displs, mpi_real4, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%bottom /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + ty - 1
            ml(j1:j2,i1-iex,k,n) = rdata(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if
    if ( ma%top /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + ty - 1
            ml(j1:j2,i2+iex,k,n) = rdata(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real4_4d_exchange_bottom_top

  subroutine real8_2d_exchange_left_bottom(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2
    integer(ik4) :: nx, ny
    integer(ik4) :: lb
    integer(ik4) :: ndx, ndy, tx, ty, sizex, sizey

    nx = j2-j1+1
    ny = i2-i1+1
    lb = nex
    if ( ma%left == mpi_proc_null ) lb = 0
    tx = ny
    ty = nx+lb
    sizex = nex*tx
    sizey = nex*ty
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    real(rk8), dimension(ndx) :: sdatax
    real(rk8), dimension(ndy) :: sdatay
    real(rk8), dimension(ndx), volatile :: rdatax
    real(rk8), dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1, ib2, iex

    !$acc data copy(ml) create(sdatax,sdatay,rdatax,rdatay)

    do concurrent ( iex = 1:nex )
      ib1 = 1 + (iex-1)*(tx)
      ib2 = ib1 + tx - 1
      sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2)
    end do
    do concurrent ( iex = 1:nex )
      ib1 = 1 + (iex-1)*(tx) + sizex
      ib2 = ib1 + tx - 1
      sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2)
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    !$acc host_data use_device(sdatax,rdatax)
    call mpi_neighbor_alltoallv(sdatax, counts, displs, mpi_real8, &
        rdatax, counts, displs, mpi_real8, cartesian_communicator, mpierr)
    !$acc end host_data
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    if ( ma%left /= mpi_proc_null ) then
      do concurrent ( iex = 1:nex )
        ib1 = 1 + (iex-1)*(tx)
        ib2 = ib1 + tx - 1
        ml(j1-iex,i1:i2) = rdatax(ib1:ib2)
      end do
    end if

    do concurrent ( iex = 1:nex )
      ib1 = 1 + (iex-1)*(ty)
      ib2 = ib1 + ty - 1
      sdatay(ib1:ib2) = ml(j1-lb:j2,i1+(iex-1))
    end do
    do concurrent ( iex = 1:nex )
      ib1 = 1 + (iex-1)*(ty) + sizey
      ib2 = ib1 + ty - 1
      sdatay(ib1:ib2) = ml(j1-lb:j2,i2-(iex-1))
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    !$acc host_data use_device(sdatay,rdatay)
    call mpi_neighbor_alltoallv(sdatay, counts, displs, mpi_real8, &
        rdatay, counts, displs, mpi_real8, cartesian_communicator, mpierr)
    !$acc end host_data
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    if ( ma%bottom /= mpi_proc_null ) then
      do concurrent ( iex = 1:nex )
        ib1 = 1 + (iex-1)*(ty)
        ib2 = ib1 + ty - 1
        ml(j1-lb:j2,i1-iex) = rdatay(ib1:ib2)
      end do
    end if
    !$acc end data
    end block exchange

  end subroutine real8_2d_exchange_left_bottom

  subroutine real8_3d_exchange_left_bottom(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2
    integer(ik4) :: nx, ny, nk
    integer(ik4) :: lb
    integer(ik4) :: ndx, ndy, tx, ty, sizex, sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    lb = nex
    if ( ma%left == mpi_proc_null ) lb = 0
    tx = ny
    ty = nx+lb
    sizex = nex*tx*nk
    sizey = nex*ty*nk
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    real(rk8), dimension(ndx) :: sdatax
    real(rk8), dimension(ndy) :: sdatay
    real(rk8), dimension(ndx), volatile :: rdatax
    real(rk8), dimension(ndy), volatile :: rdatay
    integer(ik4) :: id, ib1, ib2, iex, k

    !$acc data copy(ml) create(sdatax,sdatay,rdatax,rdatay)

    do concurrent ( iex = 1:nex, k = k1:k2 )
      id = (k-k1)*nex + iex
      ib1 = 1 + (id-1)*(tx)
      ib2 = ib1 + tx - 1
      sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2,k)
    end do
    do concurrent ( iex = 1:nex, k = k1:k2 )
      id = (k-k1)*nex + iex
      ib1 = 1 + (id-1)*(tx) + sizex
      ib2 = ib1 + tx - 1
      sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2,k)
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    !$acc host_data use_device(sdatax,rdatax)
    call mpi_neighbor_alltoallv(sdatax, counts, displs, mpi_real8, &
        rdatax, counts, displs, mpi_real8, cartesian_communicator, mpierr)
    !$acc end host_data
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    if ( ma%left /= mpi_proc_null ) then
      do concurrent ( iex = 1:nex, k = k1:k2 )
        id = (k-k1)*nex + iex
        ib1 = 1 + (id-1)*(tx)
        ib2 = ib1 + tx - 1
        ml(j1-iex,i1:i2,k) = rdatax(ib1:ib2)
      end do
    end if

    do concurrent ( iex = 1:nex, k = k1:k2 )
      id = (k-k1)*nex + iex
      ib1 = 1 + (id-1)*(ty)
      ib2 = ib1 + ty - 1
      sdatay(ib1:ib2) = ml(j1-lb:j2,i1+(iex-1),k)
    end do
    do concurrent ( iex = 1:nex, k = k1:k2 )
      id = (k-k1)*nex + iex
      ib1 = 1 + (id-1)*(ty) + sizey
      ib2 = ib1 + ty - 1
      sdatay(ib1:ib2) = ml(j1-lb:j2,i2-(iex-1),k)
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    !$acc host_data use_device(sdatay,rdatay)
    call mpi_neighbor_alltoallv(sdatay, counts, displs, mpi_real8, &
        rdatay, counts, displs, mpi_real8, cartesian_communicator, mpierr)
    !$acc end host_data
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    if ( ma%bottom /= mpi_proc_null ) then
      do concurrent ( iex = 1:nex, k = k1:k2 )
        id = (k-k1)*nex + iex
        ib1 = 1 + (id-1)*(ty)
        ib2 = ib1 + ty - 1
        ml(j1-lb:j2,i1-iex,k) = rdatay(ib1:ib2)
      end do
    end if
    !$acc end data
    end block exchange

  end subroutine real8_3d_exchange_left_bottom

  subroutine real8_4d_exchange_left_bottom(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2, n1, n2
    integer(ik4) :: nx, ny, nk, nn
    integer(ik4) :: lb
    integer(ik4) :: ndx, ndy, tx, ty, sizex, sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    nn = n2-n1+1
    lb = nex
    if ( ma%left == mpi_proc_null ) lb = 0
    tx = ny
    ty = nx+lb
    sizex = nex*tx*nk
    sizey = nex*ty*nk
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    real(rk8), dimension(ndx) :: sdatax
    real(rk8), dimension(ndy) :: sdatay
    real(rk8), dimension(ndx), volatile :: rdatax
    real(rk8), dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1, ib2, iex, k, n

    ib2 = 0
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2,k,n)
        end do
      end do
    end do
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2,k,n)
        end do
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdatax, counts, displs, mpi_real8, &
        rdatax, counts, displs, mpi_real8, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%left /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + tx - 1
            ml(j1-iex,i1:i2,k,n) = rdatax(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if

    ib2 = 0
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdatay(ib1:ib2) = ml(j1-lb:j2,i1+(iex-1),k,n)
        end do
      end do
    end do
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdatay(ib1:ib2) = ml(j1-lb:j2,i2-(iex-1),k,n)
        end do
      end do
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdatay, counts, displs, mpi_real8, &
        rdatay, counts, displs, mpi_real8, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%bottom /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + ty - 1
            ml(j1-lb:j2,i1-iex,k,n) = rdatay(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real8_4d_exchange_left_bottom

  subroutine real4_2d_exchange_left_bottom(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2
    integer(ik4) :: nx, ny
    integer(ik4) :: lb
    integer(ik4) :: ndx, ndy, tx, ty, sizex, sizey

    nx = j2-j1+1
    ny = i2-i1+1
    lb = nex
    if ( ma%left == mpi_proc_null ) lb = 0
    tx = ny
    ty = nx+lb
    sizex = nex*tx
    sizey = nex*ty
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    real(rk4), dimension(ndx) :: sdatax
    real(rk4), dimension(ndy) :: sdatay
    real(rk4), dimension(ndx), volatile :: rdatax
    real(rk4), dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1, ib2, iex

    ib2 = 0
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2)
    end do
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2)
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdatax, counts, displs, mpi_real4, &
        rdatax, counts, displs, mpi_real4, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%left /= mpi_proc_null ) then
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        ml(j1-iex,i1:i2) = rdatax(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizex
    end if

    ib2 = 0
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdatay(ib1:ib2) = ml(j1-lb:j2,i1+(iex-1))
    end do
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdatay(ib1:ib2) = ml(j1-lb:j2,i2-(iex-1))
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdatay, counts, displs, mpi_real4, &
        rdatay, counts, displs, mpi_real4, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%bottom /= mpi_proc_null ) then
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        ml(j1-lb:j2,i1-iex) = rdatay(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real4_2d_exchange_left_bottom

  subroutine real4_3d_exchange_left_bottom(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2
    integer(ik4) :: nx, ny, nk
    integer(ik4) :: lb
    integer(ik4) :: ndx, ndy, tx, ty, sizex, sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    lb = nex
    if ( ma%left == mpi_proc_null ) lb = 0
    tx = ny
    ty = nx+lb
    sizex = nex*tx*nk
    sizey = nex*ty*nk
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    real(rk4), dimension(ndx) :: sdatax
    real(rk4), dimension(ndy) :: sdatay
    real(rk4), dimension(ndx), volatile :: rdatax
    real(rk4), dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1, ib2, iex, k

    ib2 = 0
    do k = k1, k2
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2,k)
      end do
    end do
    do k = k1, k2
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2,k)
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdatax, counts, displs, mpi_real4, &
        rdatax, counts, displs, mpi_real4, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%left /= mpi_proc_null ) then
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          ml(j1-iex,i1:i2,k) = rdatax(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if

    ib2 = 0
    do k = k1, k2
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        sdatay(ib1:ib2) = ml(j1-lb:j2,i1+(iex-1),k)
      end do
    end do
    do k = k1, k2
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        sdatay(ib1:ib2) = ml(j1-lb:j2,i2-(iex-1),k)
      end do
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdatay, counts, displs, mpi_real4, &
        rdatay, counts, displs, mpi_real4, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%bottom /= mpi_proc_null ) then
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          ml(j1-lb:j2,i1-iex,k) = rdatay(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real4_3d_exchange_left_bottom

  subroutine real4_4d_exchange_left_bottom(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2, n1, n2
    integer(ik4) :: nx, ny, nk, nn
    integer(ik4) :: lb
    integer(ik4) :: ndx, ndy, tx, ty, sizex, sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    nn = n2-n1+1
    lb = nex
    if ( ma%left == mpi_proc_null ) lb = 0
    tx = ny
    ty = nx+lb
    sizex = nex*tx*nk
    sizey = nex*ty*nk
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    real(rk4), dimension(ndx) :: sdatax
    real(rk4), dimension(ndy) :: sdatay
    real(rk4), dimension(ndx), volatile :: rdatax
    real(rk4), dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1, ib2, iex, k, n

    ib2 = 0
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2,k,n)
        end do
      end do
    end do
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2,k,n)
        end do
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdatax, counts, displs, mpi_real4, &
        rdatax, counts, displs, mpi_real4, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%left /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + tx - 1
            ml(j1-iex,i1:i2,k,n) = rdatax(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if

    ib2 = 0
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdatay(ib1:ib2) = ml(j1-lb:j2,i1+(iex-1),k,n)
        end do
      end do
    end do
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdatay(ib1:ib2) = ml(j1-lb:j2,i2-(iex-1),k,n)
        end do
      end do
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdatay, counts, displs, mpi_real4, &
        rdatay, counts, displs, mpi_real4, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = 0
    if ( ma%bottom /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + ty - 1
            ml(j1-lb:j2,i1-iex,k,n) = rdatay(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real4_4d_exchange_left_bottom

  subroutine real8_2d_exchange_right_top(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2
    integer(ik4) :: nx, ny
    integer(ik4) :: rb
    integer(ik4) :: ndx, ndy, tx, ty, sizex, sizey

    nx = j2-j1+1
    ny = i2-i1+1
    rb = nex
    if ( ma%right == mpi_proc_null ) rb = 0
    tx = ny
    ty = nx+rb
    sizex = nex*tx
    sizey = nex*ty
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    real(rk8), dimension(ndx) :: sdatax
    real(rk8), dimension(ndy) :: sdatay
    real(rk8), dimension(ndx), volatile :: rdatax
    real(rk8), dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1, ib2, iex

    ib2 = 0
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2)
    end do
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2)
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdatax, counts, displs, mpi_real8, &
        rdatax, counts, displs, mpi_real8, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = sizex
    if ( ma%right /= mpi_proc_null ) then
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        ml(j2+iex,i1:i2) = rdatax(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizex
    end if

    ib2 = 0
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdatay(ib1:ib2) = ml(j1:j2+rb,i1+(iex-1))
    end do
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdatay(ib1:ib2) = ml(j1:j2+rb,i2-(iex-1))
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdatay, counts, displs, mpi_real8, &
        rdatay, counts, displs, mpi_real8, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = sizey
    if ( ma%top /= mpi_proc_null ) then
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        ml(j1:j2+rb,i2+iex) = rdatay(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real8_2d_exchange_right_top

  subroutine real8_3d_exchange_right_top(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2
    integer(ik4) :: nx, ny, nk
    integer(ik4) :: rb
    integer(ik4) :: ndx, ndy, tx, ty, sizex, sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    rb = nex
    if ( ma%right == mpi_proc_null ) rb = 0
    tx = ny
    ty = nx+rb
    sizex = nex*tx*nk
    sizey = nex*ty*nk
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    real(rk8), dimension(ndx) :: sdatax
    real(rk8), dimension(ndy) :: sdatay
    real(rk8), dimension(ndx), volatile :: rdatax
    real(rk8), dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1, ib2, iex, k

    ib2 = 0
    do k = k1, k2
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2,k)
      end do
    end do
    do k = k1, k2
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2,k)
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdatax, counts, displs, mpi_real8, &
        rdatax, counts, displs, mpi_real8, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = sizex
    if ( ma%right /= mpi_proc_null ) then
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          ml(j2+iex,i1:i2,k) = rdatax(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if

    ib2 = 0
    do k = k1, k2
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        sdatay(ib1:ib2) = ml(j1:j2+rb,i1+(iex-1),k)
      end do
    end do
    do k = k1, k2
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        sdatay(ib1:ib2) = ml(j1:j2+rb,i2-(iex-1),k)
      end do
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdatay, counts, displs, mpi_real8, &
        rdatay, counts, displs, mpi_real8, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = sizey
    if ( ma%top /= mpi_proc_null ) then
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          ml(j1:j2+rb,i2+iex,k) = rdatay(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real8_3d_exchange_right_top

  subroutine real8_4d_exchange_right_top(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2, n1, n2
    integer(ik4) :: nx, ny, nk, nn
    integer(ik4) :: rb
    integer(ik4) :: ndx, ndy, tx, ty, sizex, sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    nn = n2-n1+1
    rb = nex
    if ( ma%right == mpi_proc_null ) rb = 0
    tx = ny
    ty = nx+rb
    sizex = nex*tx*nk*nn
    sizey = nex*ty*nk*nn
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    real(rk8), dimension(ndx) :: sdatax
    real(rk8), dimension(ndy) :: sdatay
    real(rk8), dimension(ndx), volatile :: rdatax
    real(rk8), dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1, ib2, iex, k, n

    ib2 = 0
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2,k,n)
        end do
      end do
    end do
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2,k,n)
        end do
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdatax, counts, displs, mpi_real8, &
        rdatax, counts, displs, mpi_real8, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = sizex
    if ( ma%right /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + tx - 1
            ml(j2+iex,i1:i2,k,n) = rdatax(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if

    ib2 = 0
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdatay(ib1:ib2) = ml(j1:j2+rb,i1+(iex-1),k,n)
        end do
      end do
    end do
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdatay(ib1:ib2) = ml(j1:j2+rb,i2-(iex-1),k,n)
        end do
      end do
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdatay, counts, displs, mpi_real8, &
        rdatay, counts, displs, mpi_real8, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = sizey
    if ( ma%top /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + ty - 1
            ml(j1:j2+rb,i2+iex,k,n) = rdatay(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real8_4d_exchange_right_top

  subroutine real4_2d_exchange_right_top(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2
    integer(ik4) :: nx, ny
    integer(ik4) :: rb
    integer(ik4) :: ndx, ndy, tx, ty, sizex, sizey

    nx = j2-j1+1
    ny = i2-i1+1
    rb = nex
    if ( ma%right == mpi_proc_null ) rb = 0
    tx = ny
    ty = nx+rb
    sizex = nex*tx
    sizey = nex*ty
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    real(rk4), dimension(ndx) :: sdatax
    real(rk4), dimension(ndy) :: sdatay
    real(rk4), dimension(ndx), volatile :: rdatax
    real(rk4), dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1, ib2, iex

    ib2 = 0
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2)
    end do
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2)
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdatax, counts, displs, mpi_real4, &
        rdatax, counts, displs, mpi_real4, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = sizex
    if ( ma%right /= mpi_proc_null ) then
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        ml(j2+iex,i1:i2) = rdatax(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizex
    end if

    ib2 = 0
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdatay(ib1:ib2) = ml(j1:j2+rb,i1+(iex-1))
    end do
    do iex = 1, nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      sdatay(ib1:ib2) = ml(j1:j2+rb,i2-(iex-1))
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdatay, counts, displs, mpi_real4, &
        rdatay, counts, displs, mpi_real4, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = sizey
    if ( ma%top /= mpi_proc_null ) then
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        ml(j1:j2+rb,i2+iex) = rdatay(ib1:ib2)
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real4_2d_exchange_right_top

  subroutine real4_3d_exchange_right_top(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2
    integer(ik4) :: nx, ny, nk
    integer(ik4) :: rb
    integer(ik4) :: ndx, ndy, tx, ty, sizex, sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    rb = nex
    if ( ma%right == mpi_proc_null ) rb = 0
    tx = ny
    ty = nx+rb
    sizex = nex*tx*nk
    sizey = nex*ty*nk
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    real(rk4), dimension(ndx) :: sdatax
    real(rk4), dimension(ndy) :: sdatay
    real(rk4), dimension(ndx), volatile :: rdatax
    real(rk4), dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1, ib2, iex, k

    ib2 = 0
    do k = k1, k2
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2,k)
      end do
    end do
    do k = k1, k2
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2,k)
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdatax, counts, displs, mpi_real4, &
        rdatax, counts, displs, mpi_real4, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = sizex
    if ( ma%right /= mpi_proc_null ) then
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          ml(j2+iex,i1:i2,k) = rdatax(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if

    ib2 = 0
    do k = k1, k2
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        sdatay(ib1:ib2) = ml(j1:j2+rb,i1+(iex-1),k)
      end do
    end do
    do k = k1, k2
      do iex = 1, nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        sdatay(ib1:ib2) = ml(j1:j2+rb,i2-(iex-1),k)
      end do
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdatay, counts, displs, mpi_real4, &
        rdatay, counts, displs, mpi_real4, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = sizey
    if ( ma%top /= mpi_proc_null ) then
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          ml(j1:j2+rb,i2+iex,k) = rdatay(ib1:ib2)
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real4_3d_exchange_right_top

  subroutine real4_4d_exchange_right_top(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2, n1, n2
    integer(ik4) :: nx, ny, nk, nn
    integer(ik4) :: rb
    integer(ik4) :: ndx, ndy, tx, ty, sizex, sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    nn = n2-n1+1
    rb = nex
    if ( ma%right == mpi_proc_null ) rb = 0
    tx = ny
    ty = nx+rb
    sizex = nex*tx*nk*nn
    sizey = nex*ty*nk*nn
    ndx = 2*sizex
    ndy = 2*sizey

    exchange : block

    integer(ik4), dimension(4) :: counts, displs
    real(rk4), dimension(ndx) :: sdatax
    real(rk4), dimension(ndy) :: sdatay
    real(rk4), dimension(ndx), volatile :: rdatax
    real(rk4), dimension(ndy), volatile :: rdatay
    integer(ik4) :: ib1, ib2, iex, k, n

    ib2 = 0
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2,k,n)
        end do
      end do
    end do
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          sdatax(ib1:ib2) = ml(j2-(iex-1),i1:i2,k,n)
        end do
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_neighbor_alltoallv(sdatax, counts, displs, mpi_real4, &
        rdatax, counts, displs, mpi_real4, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = sizex
    if ( ma%right /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + tx - 1
            ml(j2+iex,i1:i2,k,n) = rdatax(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizex
    end if

    ib2 = 0
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdatay(ib1:ib2) = ml(j1:j2+rb,i1+(iex-1),k,n)
        end do
      end do
    end do
    do n = n1, n2
      do k = k1, k2
        do iex = 1, nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          sdatay(ib1:ib2) = ml(j1:j2+rb,i2-(iex-1),k,n)
        end do
      end do
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_neighbor_alltoallv(sdatay, counts, displs, mpi_real4, &
        rdatay, counts, displs, mpi_real4, cartesian_communicator, mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_neighbor_alltoallv error.')
    end if
#endif

    ib2 = sizey
    if ( ma%top /= mpi_proc_null ) then
      do n = n1, n2
        do k = k1, k2
          do iex = 1, nex
            ib1 = ib2 + 1
            ib2 = ib1 + ty - 1
            ml(j1:j2+rb,i2+iex,k,n) = rdatay(ib1:ib2)
          end do
        end do
      end do
    else
      ib2 = ib2 + sizey
    end if

    end block exchange

  end subroutine real4_4d_exchange_right_top

#else

  subroutine real8_2d_exchange(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2
    integer(ik4) :: isize, jsize, ssize, j, i, ib, irc, ipos
    integer(ik4), dimension(16) :: req

    isize = i2-i1+1
    jsize = j2-j1+1

    irc = 0
    ipos = 1
    req = mpi_request_null

    ! Total max possible communication size
    ! set the size of the dummy exchange vector
    ! nex is the width of the exchange stencil
    ! isize is the height of a block in the N-S direction
    ! jsize is the height of a block in the W-E direction
    ssize = nex * (4*nex+2*isize+2*jsize)
    if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'real8_2d_exchange')
    end if
    if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'real8_2d_exchange')
    end if

    if ( ma%bandflag ) then
      ssize = nex*isize
      !********************************************************
      ! (1) Send right boundary to left side of right neighbor
      !********************************************************
      ! loop over the given latitudes and number of exchange blocks
      ! and unravel the right boundary into a vector
      ib = 1
      do i = i1, i2
        do j = 1, nex
          r8vector1(ib) = ml(j2-j+1,i)
          ib = ib + 1
        end do
      end do
      ! do the send-right exchange
      call cyclic_exchange_array(r8vector1,r8vector2, &
                                 ssize,ma%right,ma%left,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the left boundary of this block
      ib = 1
      do i = i1, i2
        do j = 1, nex
          ml(j1-j,i) = r8vector2(ib)
          ib = ib + 1
        end do
      end do

      !********************************************************
      ! (2) Send left boundary to right side of left neighbor
      !********************************************************
      ! loop over the given latitudes and number of exchange blocks
      ! and unravel the left boundary into a vector
      ib = 1
      do i = i1, i2
        do j = 1, nex
          r8vector1(ib) = ml(j1+j-1,i)
          ib = ib + 1
        end do
      end do
      ! do the send-left exchange
      call cyclic_exchange_array(r8vector1,r8vector2, &
                                 ssize,ma%left,ma%right,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the right boundary of this block
      ib = 1
      do i = i1, i2
        do j = 1, nex
          ml(j2+j,i) = r8vector2(ib)
          ib = ib + 1
        end do
      end do

      if ( ma%crmflag ) then
        !*********************************************************
        ! (3) Send bottom boundary to top side of bottom neighbor
        !*********************************************************
        ssize = nex*jsize
        ! loop over the given longitudes and number of exchange blocks
        ! and unravel the bottom boundary into a vector
        ib = 1
        do i = 1, nex
          do j = j1, j2
            r8vector1(ib) = ml(j,i1+i-1)
            ib = ib + 1
          end do
        end do

        ! do the send-up exchange
        call cyclic_exchange_array(r8vector1,r8vector2, &
                                   ssize,ma%bottom,ma%top,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the top boundary of this block
        ib = 1
        do i = 1, nex
          do j = j1, j2
            ml(j,i2+i) = r8vector2(ib)
            ib = ib + 1
          end do
        end do

        !******************************************************
        ! (4) Send top boundary to bottom side of top neighbor
        !******************************************************
        ! loop over the given longitudes and number of exchange blocks
        ! and unravel the top boundary into a vector
        ib = 1
        do i = 1, nex
          do j = j1, j2
            r8vector1(ib) = ml(j,i2-i+1)
            ib = ib + 1
          end do
        end do

        ! do the send-down exchange
        call cyclic_exchange_array(r8vector1,r8vector2, &
                                   ssize,ma%top,ma%bottom,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the bottom boundary of this block
        ib = 1
        do i = 1, nex
          do j = j1, j2
            ml(j,i1-i) = r8vector2(ib)
            ib = ib + 1
          end do
        end do

        !*******************************************************************
        ! (5) Send bottom-right boundary to top-left side of bottom-right
        !     neighbor
        !*******************************************************************
        ! set the size of the exchange vector to the
        ! a square the length/width of the exchange stencil
        ssize = nex*nex
        ! loop over the exchange block and unravel it into a vector
        ib = 1
        do i = 1, nex
          do j = 1, nex
            r8vector1(ib) = ml(j2-j+1,i1+i-1)
            ib = ib + 1
          end do
        end do

        ! do the send down-right exchange
        call cyclic_exchange_array(r8vector1,r8vector2, &
                                   ssize,ma%bottomright,ma%topleft,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the top-left boundary of this block
        ib = 1
        do i = 1, nex
          do j = 1, nex
            ml(j1-j,i2+i) = r8vector2(ib)
            ib = ib + 1
          end do
        end do

        !**********************************************************************
        ! (6) Send top-left boundary to bottom-right side of top-left neighbor
        !**********************************************************************
        ! loop over the exchange block and unravel it into a vector
        ib = 1
        do i = 1, nex
          do j = 1, nex
            r8vector1(ib) = ml(j1+j-1,i2-i+1)
            ib = ib + 1
          end do
        end do

        ! do the send up-left exchange
        call cyclic_exchange_array(r8vector1,r8vector2, &
                                   ssize,ma%topleft,ma%bottomright,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the bottom-right boundary of this block
        ib = 1
        do i = 1, nex
          do j = 1, nex
            ml(j2+j,i1-i) = r8vector2(ib)
            ib = ib + 1
          end do
        end do

        !******************************************************************
        ! (7) Send bottom-left boundary to top-right side of bottom-left
        !     neighbor
        !******************************************************************
        ! loop over the exchange block and unravel it into a vector
        ib = 1
        do i = 1, nex
          do j = 1, nex
            r8vector1(ib) = ml(j1+j-1,i1+i-1)
            ib = ib + 1
          end do
        end do

        ! do the send down-left exchange
        call cyclic_exchange_array(r8vector1,r8vector2, &
                                   ssize,ma%bottomleft,ma%topright,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the top-right boundary of this block
        ib = 1
        do i = 1, nex
          do j = 1, nex
            ml(j2+j,i2+i) = r8vector2(ib)
            ib = ib + 1
          end do
        end do

        !*****************************************************************
        ! (8) Send top-right boundary to bottom-left side of top-right
        !     neighbor
        !******************************************************************
        ! loop over the exchange block and unravel it into a vector
        ib = 1
        do i = 1, nex
          do j = 1, nex
            r8vector1(ib) = ml(j2-j+1,i2-i+1)
            ib = ib + 1
          end do
        end do

        ! do the send up-right exchange
        call cyclic_exchange_array(r8vector1,r8vector2, &
                                   ssize,ma%topright,ma%bottomleft,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the bottom-left boundary of this block
        ib = 1
        do i = 1, nex
          do j = 1, nex
            ml(j1-j,i1-i) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end if

    else

      ! RIGHT AND LEFT EXCHANGE

      if ( ma%right /= mpi_proc_null) then
        ssize = nex*isize
        ib = ipos
        do i = i1, i2
          do j = 1, nex
            r8vector1(ib) = ml(j2-j+1,i)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%right,tag_rl,tag_lr, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if
      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*isize
        ib = ipos
        do i = i1, i2
          do j = 1, nex
            r8vector1(ib) = ml(j1+j-1,i)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%left,tag_lr,tag_rl, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if

    end if

    if ( .not. ma%crmflag ) then

      ! BOTTOM AND TOP EXCHANGE

      if ( ma%top /= mpi_proc_null) then
        ssize = nex*jsize
        ib = ipos
        do i = 1, nex
          do j = j1, j2
            r8vector1(ib) = ml(j,i2-i+1)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%top,tag_tb,tag_bt, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if
      if ( ma%bottom /= mpi_proc_null) then
        ssize = nex*jsize
        ib = ipos
        do i = 1, nex
          do j = j1, j2
            r8vector1(ib) = ml(j,i1+i-1)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%bottom,tag_bt,tag_tb, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if

      ! TOP LEFT WITH BOTTOM RIGHT EXCHANGE

      if ( ma%topleft /= mpi_proc_null ) then
        ssize = nex*nex
        ib = ipos
        do i = 1, nex
          do j = 1, nex
            r8vector1(ib) = ml(j1+j-1,i2-i+1)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%topleft,tag_tlbr,tag_brtl, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if
      if ( ma%bottomright /= mpi_proc_null ) then
        ssize = nex*nex
        ib = ipos
        do i = 1, nex
          do j = 1, nex
            r8vector1(ib) = ml(j2-j+1,i1+i-1)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%bottomright,tag_brtl,tag_tlbr, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if

      ! BOTTOM LEFT WITH TOP RIGHT EXCHANGE

      if ( ma%topright /= mpi_proc_null ) then
        ssize = nex*nex
        ib = ipos
        do i = 1, nex
          do j = 1, nex
            r8vector1(ib) = ml(j2-j+1,i2-i+1)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%topright,tag_trbl,tag_bltr, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if
      if ( ma%bottomleft /= mpi_proc_null ) then
        ssize = nex*nex
        ib = ipos
        do i = 1, nex
          do j = 1, nex
            r8vector1(ib) = ml(j1+j-1,i1+i-1)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%bottomleft,tag_bltr,tag_trbl, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if
    end if

    ! Finalize non cyclic comms

    if ( irc /= 0 ) then
      ipos = 1
      call mpi_waitall(16,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      if ( .not. ma%bandflag ) then
        if ( ma%right /= mpi_proc_null) then
          ib = ipos
          ssize = nex*isize
          do i = i1, i2
            do j = 1, nex
              ml(j2+j,i) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%left /= mpi_proc_null ) then
          ib = ipos
          ssize = nex*isize
          do i = i1, i2
            do j = 1, nex
              ml(j1-j,i) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
      if ( .not. ma%crmflag ) then
        if ( ma%top /= mpi_proc_null) then
          ib = ipos
          ssize = nex*jsize
          do i = 1, nex
            do j = j1, j2
              ml(j,i2+i) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%bottom /= mpi_proc_null) then
          ib = ipos
          ssize = nex*jsize
          do i = 1, nex
            do j = j1, j2
              ml(j,i1-i) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%topleft /= mpi_proc_null ) then
          ib = ipos
          ssize = nex*nex
          do i = 1, nex
            do j = 1, nex
              ml(j1-j,i2+i) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%bottomright /= mpi_proc_null ) then
          ib = ipos
          ssize = nex*nex
          do i = 1, nex
            do j = 1, nex
              ml(j2+j,i1-i) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%topright /= mpi_proc_null ) then
          ib = ipos
          ssize =nex*nex
          do i = 1, nex
            do j = 1, nex
              ml(j2+j,i2+i) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%bottomleft /= mpi_proc_null ) then
          ib = ipos
          ssize = nex*nex
          do i = 1, nex
            do j = 1, nex
              ml(j1-j,i1-i) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
    end if
  end subroutine real8_2d_exchange

  subroutine real8_3d_exchange(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2
    integer(ik4) :: isize, jsize, ksize, ssize, j, i, k, ib, irc, ipos
    integer(ik4), dimension(16) :: req

    isize = i2-i1+1
    jsize = j2-j1+1
    ksize = k2-k1+1

    irc = 0
    ipos = 1
    req = mpi_request_null

    ! Total max possible communication size
    ! set the size of the dummy exchange vector
    ! nex is the width of the exchange stencil
    ! isize is the height of a block in the N-S direction
    ! jsize is the height of a block in the W-E direction
    ! ksize is the height of a block in the T-B direction
    ssize = nex * ksize * (4*nex+2*isize+2*jsize)
    if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'real8_3d_exchange')
    end if
    if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'real8_3d_exchange')
    end if

    if ( ma%bandflag ) then
      ! set the size of the dummy exchange vector
      ! nex is the width of the exchange stencil
      ! isize is the height of a block in the N-S direction
      ssize = nex*isize*ksize

      !********************************************************
      ! (1) Send right boundary to left side of right neighbor
      !********************************************************
      ! loop over the given latitudes and number of exchange blocks
      ! and unravel the right boundary into a vector
      do k = k1, k2
        do i = i1, i2
          do j = 1, nex
            ib = j + (i - i1) * nex + (k - k1) * (i2 - i1 + 1) * nex
            r8vector1(ib) = ml(j2-j+1,i,k)
          end do
        end do
      end do
      ! do the send-right exchange
      call cyclic_exchange_array(r8vector1,r8vector2, &
                                 ssize,ma%right,ma%left,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the left boundary of this block
      do k = k1, k2
        do i = i1, i2
          do j = 1, nex
            ib = j + (i - i1) * nex + (k - k1) * (i2 - i1 + 1) * nex
            ml(j1-j,i,k) = r8vector2(ib)
          end do
        end do
      end do

      !********************************************************
      ! (2) Send left boundary to right side of left neighbor
      !********************************************************
      ! loop over the given latitudes and number of exchange blocks
      ! and unravel the left boundary into a vector
      do k = k1, k2
        do i = i1, i2
          do j = 1, nex
            ib = j + (i - i1) * nex + (k - k1) * (i2 - i1 + 1) * nex
            r8vector1(ib) = ml(j1+j-1,i,k)
          end do
        end do
      end do
      ! do the send-left exchange
      call cyclic_exchange_array(r8vector1,r8vector2, &
                                 ssize,ma%left,ma%right,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the right boundary of this block
      do k = k1, k2
        do i = i1, i2
          do j = 1, nex
            ib = j + (i - i1) * nex + (k - k1) * (i2 - i1 + 1) * nex
            ml(j2+j,i,k) = r8vector2(ib)
          end do
        end do
      end do

      if ( ma%crmflag ) then
        !*********************************************************
        ! (3) Send bottom boundary to top side of bottom neighbor
        !*********************************************************
        ! set the size of the dummy exchange vector
        ! nex is the width of the exchange stencil
        ! jsize is the width of a block in the E-W direction
        ssize = nex*jsize*ksize
        ! loop over the given longitudes and number of exchange blocks
        ! and unravel the bottom boundary into a vector
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              ib = (j - j1 + 1) + (i - 1) * (j2 - j1 + 1) + &
                (k - k1) * (j2 - j1 + 1) * nex
              r8vector1(ib) = ml(j,i1+i-1,k)
            end do
          end do
        end do

        ! do the send-up exchange
        call cyclic_exchange_array(r8vector1,r8vector2, &
                                   ssize,ma%bottom,ma%top,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the top boundary of this block
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              ib = (j - j1 + 1) + (i - 1) * (j2 - j1 + 1) + &
                (k - k1) * (j2 - j1 + 1) * nex
              ml(j,i2+i,k) = r8vector2(ib)
            end do
          end do
        end do

        !******************************************************
        ! (4) Send top boundary to bottom side of top neighbor
        !******************************************************
        ! loop over the given longitudes and number of exchange blocks
        ! and unravel the top boundary into a vector
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              ib = (j - j1 + 1) + (i - 1) * (j2 - j1 + 1) + &
                (k - k1) * (j2 - j1 + 1) * nex
              r8vector1(ib) = ml(j,i2-i+1,k)
            end do
          end do
        end do

        ! do the send-down exchange
        call cyclic_exchange_array(r8vector1,r8vector2, &
                                   ssize,ma%top,ma%bottom,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the bottom boundary of this block
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              ib = (j - j1 + 1) + (i - 1) * (j2 - j1 + 1) + &
                (k - k1) * (j2 - j1 + 1) * nex
              ml(j,i1-i,k) = r8vector2(ib)
            end do
          end do
        end do

        !*******************************************************************
        ! (5) Send bottom-right boundary to top-left side of bottom-right
        !     neighbor
        !*******************************************************************
        ! set the size of the exchange vector to the
        ! a square the length/width of the exchange stencil
        ssize = nex*nex*ksize
        ! loop over the exchange block and unravel it into a vector
        do k = k1, k2
          do i = 1, nex
            do j = 1, nex
              ib = j + (i - 1) * nex + (k - k1) * nex * nex
              r8vector1(ib) = ml(j2-j+1,i1+i-1,k)
            end do
          end do
        end do

        ! do the send down-right exchange
        call cyclic_exchange_array(r8vector1,r8vector2, &
                                   ssize,ma%bottomright,ma%topleft,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the top-left boundary of this block
        do k = k1, k2
          do i = 1, nex
            do j = 1, nex
              ib = j + (i - 1) * nex + (k - k1) * nex * nex
              ml(j1-j,i2+i,k) = r8vector2(ib)
            end do
          end do
        end do

        !**********************************************************************
        ! (6) Send top-left boundary to bottom-right side of top-left neighbor
        !**********************************************************************
        ! loop over the exchange block and unravel it into a vector
        do k = k1, k2
          do i = 1, nex
            do j = 1, nex
              ib = j + (i - 1) * nex + (k - k1) * nex * nex
              r8vector1(ib) = ml(j1+j-1,i2-i+1,k)
            end do
          end do
        end do

        ! do the send up-left exchange
        call cyclic_exchange_array(r8vector1,r8vector2, &
                                   ssize,ma%topleft,ma%bottomright,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the bottom-right boundary of this block
        do k = k1, k2
          do i = 1, nex
            do j = 1, nex
              ib = j + (i - 1) * nex + (k - k1) * nex * nex
              ml(j2+j,i1-i,k) = r8vector2(ib)
            end do
          end do
        end do

        !******************************************************************
        ! (7) Send bottom-left boundary to top-right side of bottom-left
        !     neighbor
        !******************************************************************
        ! loop over the exchange block and unravel it into a vector
        do k = k1, k2
          do i = 1, nex
            do j = 1, nex
              ib = j + (i - 1) * nex + (k - k1) * nex * nex
              r8vector1(ib) = ml(j1+j-1,i1+i-1,k)
            end do
          end do
        end do

        ! do the send down-left exchange
        call cyclic_exchange_array(r8vector1,r8vector2, &
                                   ssize,ma%bottomleft,ma%topright,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the top-right boundary of this block
        do k = k1, k2
          do i = 1, nex
            do j = 1, nex
              ib = j + (i - 1) * nex + (k - k1) * nex * nex
              ml(j2+j,i2+i,k) = r8vector2(ib)
            end do
          end do
        end do

        !*****************************************************************
        ! (8) Send top-right boundary to bottom-left side of top-right
        !     neighbor
        !******************************************************************
        ! loop over the exchange block and unravel it into a vector
        do k = k1, k2
          do i = 1, nex
            do j = 1, nex
              ib = j + (i - 1) * nex + (k - k1) * nex * nex
              r8vector1(ib) = ml(j2-j+1,i2-i+1,k)
            end do
          end do
        end do

        ! do the send up-right exchange
        call cyclic_exchange_array(r8vector1,r8vector2, &
                                   ssize,ma%topright,ma%bottomleft,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the bottom-left boundary of this block
        do k = k1, k2
          do i = 1, nex
            do j = 1, nex
              ib = j + (i - 1) * nex + (k - k1) * nex * nex
              ml(j1-j,i1-i,k) = r8vector2(ib)
            end do
          end do
        end do
      end if

    else

      ! RIGHT AND LEFT EXCHANGE

      if ( ma%right /= mpi_proc_null) then
        ssize = nex*isize*ksize
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              ib = j + (i - i1) * nex + &
                (k - k1) * (i2 - i1 + 1) * nex + ipos - 1
              r8vector1(ib) = ml(j2-j+1,i,k)
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%right,tag_rl,tag_lr, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if
      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*isize*ksize
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              ib = j + (i - i1) * nex + &
                (k - k1) * (i2 - i1 + 1) * nex + ipos - 1
              r8vector1(ib) = ml(j1+j-1,i,k)
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%left,tag_lr,tag_rl, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if

    end if

    if ( .not. ma%crmflag ) then

      ! BOTTOM AND TOP EXCHANGE

      if ( ma%top /= mpi_proc_null) then
        ssize = nex*jsize*ksize
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              ib = (j - j1 + 1) + (i - 1) * (j2 - j1 + 1) + &
                (k - k1) * (j2 - j1 + 1) * nex + ipos - 1
              r8vector1(ib) = ml(j,i2-i+1,k)
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%top,tag_tb,tag_bt, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if
      if ( ma%bottom /= mpi_proc_null) then
        ssize = nex*jsize*ksize
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              ib = (j - j1 + 1) + (i - 1) * (j2 - j1 + 1) + &
                (k - k1) * (j2 - j1 + 1) * nex + ipos - 1
              r8vector1(ib) = ml(j,i1+i-1,k)
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%bottom,tag_bt,tag_tb, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if

      ! TOP LEFT WITH BOTTOM RIGHT EXCHANGE

      if ( ma%topleft /= mpi_proc_null ) then
        ssize = nex*nex*ksize
        do k = k1, k2
          do i = 1, nex
            do j = 1, nex
              ib = j + (i - 1) * nex + (k - k1) * nex * nex + ipos - 1
              r8vector1(ib) = ml(j1+j-1,i2-i+1,k)
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%topleft,tag_tlbr,tag_brtl, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if
      if ( ma%bottomright /= mpi_proc_null ) then
        ssize = nex*nex*ksize
        do k = k1, k2
          do i = 1, nex
            do j = 1, nex
              ib = j + (i - 1) * nex + (k - k1) * nex * nex + ipos - 1
              r8vector1(ib) = ml(j2-j+1,i1+i-1,k)
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%bottomright,tag_brtl,tag_tlbr, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if

      ! BOTTOM LEFT WITH TOP RIGHT EXCHANGE

      if ( ma%topright /= mpi_proc_null ) then
        ssize = nex*nex*ksize
        do k = k1, k2
          do i = 1, nex
            do j = 1, nex
              ib = j + (i - 1) * nex + (k - k1) * nex * nex + ipos - 1
              r8vector1(ib) = ml(j2-j+1,i2-i+1,k)
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%topright,tag_trbl,tag_bltr, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if
      if ( ma%bottomleft /= mpi_proc_null ) then
        ssize = nex*nex*ksize
        do k = k1, k2
          do i = 1, nex
            do j = 1, nex
              ib = j + (i - 1) * nex + (k - k1) * nex * nex + ipos - 1
              r8vector1(ib) = ml(j1+j-1,i1+i-1,k)
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%bottomleft,tag_bltr,tag_trbl, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if
    end if

    ! Finalize non cyclic comms

    if ( irc /= 0 ) then
      ipos = 1
      call mpi_waitall(16,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      if ( .not. ma%bandflag ) then
        if ( ma%right /= mpi_proc_null) then
          ssize = nex*isize*ksize
          do k = k1, k2
            do i = i1, i2
              do j = 1, nex
                ib = j + (i - i1) * nex + &
                  (k - k1) * (i2 - i1 + 1) * nex + ipos - 1
                ml(j2+j,i,k) = r8vector2(ib)
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%left /= mpi_proc_null ) then
          ssize = nex*isize*ksize
          do k = k1, k2
            do i = i1, i2
              do j = 1, nex
                ib = j + (i - i1) * nex + &
                  (k - k1) * (i2 - i1 + 1) * nex + ipos - 1
                ml(j1-j,i,k) = r8vector2(ib)
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
      if ( .not. ma%crmflag ) then
        if ( ma%top /= mpi_proc_null) then
          ssize = nex*jsize*ksize
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                ib = (j - j1 + 1) + (i - 1) * (j2 - j1 + 1) + &
                  (k - k1) * (j2 - j1 + 1) * nex + ipos - 1
                ml(j,i2+i,k) = r8vector2(ib)
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%bottom /= mpi_proc_null) then
          ssize = nex*jsize*ksize
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                ib = (j - j1 + 1) + (i - 1) * (j2 - j1 + 1) + &
                  (k - k1) * (j2 - j1 + 1) * nex + ipos - 1
                ml(j,i1-i,k) = r8vector2(ib)
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%topleft /= mpi_proc_null ) then
          ssize = nex*nex*ksize
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                ib = j + (i - 1) * nex + (k - k1) * nex * nex + ipos - 1
                ml(j1-j,i2+i,k) = r8vector2(ib)
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%bottomright /= mpi_proc_null ) then
          ssize = nex*nex*ksize
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                ib = j + (i - 1) * nex + (k - k1) * nex * nex + ipos - 1
                ml(j2+j,i1-i,k) = r8vector2(ib)
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%topright /= mpi_proc_null ) then
          ssize =nex*nex*ksize
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                ib = j + (i - 1) * nex + (k - k1) * nex * nex + ipos - 1
                ml(j2+j,i2+i,k) = r8vector2(ib)
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%bottomleft /= mpi_proc_null ) then
          ssize = nex*nex*ksize
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                ib = j + (i - 1) * nex + (k - k1) * nex * nex + ipos - 1
                ml(j1-j,i1-i,k) = r8vector2(ib)
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
    end if
  end subroutine real8_3d_exchange

  subroutine real8_4d_exchange(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2, n1, n2
    integer(ik4) :: isize, jsize, ksize, nsize, ssize
    integer(ik4) :: j, i, k, n, ib, irc, ipos
    integer(ik4), dimension(16) :: req

    isize = i2-i1+1
    jsize = j2-j1+1
    ksize = k2-k1+1
    nsize = n2-n1+1

    irc = 0
    ipos = 1
    req = mpi_request_null

    ! Total max possible communication size
    ! set the size of the dummy exchange vector
    ! nex is the width of the exchange stencil
    ! isize is the height of a block in the N-S direction
    ! jsize is the height of a block in the W-E direction
    ! ksize is the height of a block in the T-B direction
    ! nsize is the height of a block in the tracer direction
    ssize = nex * ksize * nsize * (4*nex+2*isize+2*jsize)
    if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'real8_4d_exchange')
    end if
    if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'real8_4d_exchange')
    end if

    if ( ma%bandflag ) then
      ! set the size of the dummy exchange vector
      ! nex is the width of the exchange stencil
      ! isize is the height of a block in the N-S direction
      ssize = nex*isize*ksize*nsize

      !********************************************************
      ! (1) Send right boundary to left side of right neighbor
      !********************************************************
      ! loop over the given latitudes and number of exchange blocks
      ! and unravel the right boundary into a vector
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              r8vector1(ib) = ml(j2-j+1,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      ! do the send-right exchange
      call cyclic_exchange_array(r8vector1,r8vector2, &
                                 ssize,ma%right,ma%left,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the left boundary of this block
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              ml(j1-j,i,k,n) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do

      !********************************************************
      ! (2) Send left boundary to right side of left neighbor
      !********************************************************
      ! loop over the given latitudes and number of exchange blocks
      ! and unravel the left boundary into a vector
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              r8vector1(ib) = ml(j1+j-1,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      ! do the send-left exchange
      call cyclic_exchange_array(r8vector1,r8vector2, &
                                 ssize,ma%left,ma%right,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the right boundary of this block
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              ml(j2+j,i,k,n) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do

      if ( ma%crmflag ) then
        !*********************************************************
        ! (3) Send bottom boundary to top side of bottom neighbor
        !*********************************************************
        ! set the size of the dummy exchange vector
        ! nex is the width of the exchange stencil
        ! jsize is the width of a block in the E-W direction
        ssize = nex*jsize*ksize*nsize
        ! loop over the given longitudes and number of exchange blocks
        ! and unravel the bottom boundary into a vector
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                r8vector1(ib) = ml(j,i1+i-1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do

        ! do the send-up exchange
        call cyclic_exchange_array(r8vector1,r8vector2, &
                                   ssize,ma%bottom,ma%top,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the top boundary of this block
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                ml(j,i2+i,k,n) = r8vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do

        !******************************************************
        ! (4) Send top boundary to bottom side of top neighbor
        !******************************************************
        ! loop over the given longitudes and number of exchange blocks
        ! and unravel the top boundary into a vector
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                r8vector1(ib) = ml(j,i2-i+1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do

        ! do the send-down exchange
        call cyclic_exchange_array(r8vector1,r8vector2, &
                                   ssize,ma%top,ma%bottom,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the bottom boundary of this block
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                ml(j,i1-i,k,n) = r8vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do

        !*******************************************************************
        ! (5) Send bottom-right boundary to top-left side of bottom-right
        !     neighbor
        !*******************************************************************
        ! set the size of the exchange vector to the
        ! a square the length/width of the exchange stencil
        ssize = nex*nex*ksize*nsize
        ! loop over the exchange block and unravel it into a vector
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                r8vector1(ib) = ml(j2-j+1,i1+i-1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do

        ! do the send down-right exchange
        call cyclic_exchange_array(r8vector1,r8vector2, &
                                   ssize,ma%bottomright,ma%topleft,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the top-left boundary of this block
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                ml(j1-j,i2+i,k,n) = r8vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do

        !**********************************************************************
        ! (6) Send top-left boundary to bottom-right side of top-left neighbor
        !**********************************************************************
        ! loop over the exchange block and unravel it into a vector
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                r8vector1(ib) = ml(j1+j-1,i2-i+1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do

        ! do the send up-left exchange
        call cyclic_exchange_array(r8vector1,r8vector2, &
                                   ssize,ma%topleft,ma%bottomright,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the bottom-right boundary of this block
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                ml(j2+j,i1-i,k,n) = r8vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do

        !******************************************************************
        ! (7) Send bottom-left boundary to top-right side of bottom-left
        !     neighbor
        !******************************************************************
        ! loop over the exchange block and unravel it into a vector
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                r8vector1(ib) = ml(j1+j-1,i1+i-1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do

        ! do the send down-left exchange
        call cyclic_exchange_array(r8vector1,r8vector2, &
                                   ssize,ma%bottomleft,ma%topright,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the top-right boundary of this block
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                ml(j2+j,i2+i,k,n) = r8vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do

        !*****************************************************************
        ! (8) Send top-right boundary to bottom-left side of top-right
        !     neighbor
        !******************************************************************
        ! loop over the exchange block and unravel it into a vector
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                r8vector1(ib) = ml(j2-j+1,i2-i+1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do

        ! do the send up-right exchange
        call cyclic_exchange_array(r8vector1,r8vector2, &
                                   ssize,ma%topright,ma%bottomleft,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the bottom-left boundary of this block
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                ml(j1-j,i1-i,k,n) = r8vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do
      end if

    else

      ! RIGHT AND LEFT EXCHANGE

      if ( ma%right /= mpi_proc_null) then
        ssize = nex*isize*ksize*nsize
        ib = ipos
        do n = n1, n2
          do k = k1, k2
            do i = i1, i2
              do j = 1, nex
                r8vector1(ib) = ml(j2-j+1,i,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%right,tag_rl,tag_lr, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if
      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*isize*ksize*nsize
        ib = ipos
        do n = n1, n2
          do k = k1, k2
            do i = i1, i2
              do j = 1, nex
                r8vector1(ib) = ml(j1+j-1,i,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%left,tag_lr,tag_rl, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if

    end if

    if ( .not. ma%crmflag ) then

      ! BOTTOM AND TOP EXCHANGE

      if ( ma%top /= mpi_proc_null) then
        ssize = nex*jsize*ksize*nsize
        ib = ipos
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                r8vector1(ib) = ml(j,i2-i+1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%top,tag_tb,tag_bt, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if
      if ( ma%bottom /= mpi_proc_null) then
        ssize = nex*jsize*ksize*nsize
        ib = ipos
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                r8vector1(ib) = ml(j,i1+i-1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%bottom,tag_bt,tag_tb, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if

      ! TOP LEFT WITH BOTTOM RIGHT EXCHANGE

      if ( ma%topleft /= mpi_proc_null ) then
        ssize = nex*nex*ksize*nsize
        ib = ipos
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                r8vector1(ib) = ml(j1+j-1,i2-i+1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%topleft,tag_tlbr,tag_brtl, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if
      if ( ma%bottomright /= mpi_proc_null ) then
        ssize = nex*nex*ksize*nsize
        ib = ipos
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                r8vector1(ib) = ml(j2-j+1,i1+i-1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%bottomright,tag_brtl,tag_tlbr, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if

      ! BOTTOM LEFT WITH TOP RIGHT EXCHANGE

      if ( ma%topright /= mpi_proc_null ) then
        ssize = nex*nex*ksize*nsize
        ib = ipos
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                r8vector1(ib) = ml(j2-j+1,i2-i+1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%topright,tag_trbl,tag_bltr, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if
      if ( ma%bottomleft /= mpi_proc_null ) then
        ssize = nex*nex*ksize*nsize
        ib = ipos
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                r8vector1(ib) = ml(j1+j-1,i1+i-1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%bottomleft,tag_bltr,tag_trbl, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if
    end if

    ! Finalize non cyclic comms

    if ( irc /= 0 ) then
      ipos = 1
      call mpi_waitall(16,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      if ( .not. ma%bandflag ) then
        if ( ma%right /= mpi_proc_null) then
          ib = ipos
          ssize = nex*isize*ksize*nsize
          do n = n1, n2
            do k = k1, k2
              do i = i1, i2
                do j = 1, nex
                  ml(j2+j,i,k,n) = r8vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%left /= mpi_proc_null ) then
          ib = ipos
          ssize = nex*isize*ksize*nsize
          do n = n1, n2
            do k = k1, k2
              do i = i1, i2
                do j = 1, nex
                  ml(j1-j,i,k,n) = r8vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
      if ( .not. ma%crmflag ) then
        if ( ma%top /= mpi_proc_null) then
          ib = ipos
          ssize = nex*jsize*ksize*nsize
          do n = n1, n2
            do k = k1, k2
              do i = 1, nex
                do j = j1, j2
                  ml(j,i2+i,k,n) = r8vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%bottom /= mpi_proc_null) then
          ib = ipos
          ssize = nex*jsize*ksize*nsize
          do n = n1, n2
            do k = k1, k2
              do i = 1, nex
                do j = j1, j2
                  ml(j,i1-i,k,n) = r8vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%topleft /= mpi_proc_null ) then
          ib = ipos
          ssize = nex*nex*ksize*nsize
          do n = n1, n2
            do k = k1, k2
             do i = 1, nex
                do j = 1, nex
                  ml(j1-j,i2+i,k,n) = r8vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%bottomright /= mpi_proc_null ) then
          ib = ipos
          ssize = nex*nex*ksize*nsize
          do n = n1, n2
            do k = k1, k2
              do i = 1, nex
                do j = 1, nex
                  ml(j2+j,i1-i,k,n) = r8vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%topright /= mpi_proc_null ) then
          ib = ipos
          ssize =nex*nex*ksize*nsize
          do n = n1, n2
            do k = k1, k2
              do i = 1, nex
                do j = 1, nex
                  ml(j2+j,i2+i,k,n) = r8vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%bottomleft /= mpi_proc_null ) then
          ib = ipos
          ssize = nex*nex*ksize*nsize
          do n = n1, n2
            do k = k1, k2
              do i = 1, nex
                do j = 1, nex
                  ml(j1-j,i1-i,k,n) = r8vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
    end if
  end subroutine real8_4d_exchange

  subroutine real4_2d_exchange(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2
    integer(ik4) :: isize, jsize, ssize, j, i, ib, irc, ipos
    integer(ik4), dimension(16) :: req

    isize = i2-i1+1
    jsize = j2-j1+1

    irc = 0
    ipos = 1
    req = mpi_request_null

    ! Total max possible communication size
    ! set the size of the dummy exchange vector
    ! nex is the width of the exchange stencil
    ! isize is the height of a block in the N-S direction
    ! jsize is the height of a block in the W-E direction
    ssize = nex * (4*nex+2*isize+2*jsize)
    if ( size(r4vector1) < ssize ) then
      call getmem1d(r4vector1,1,ssize,'real4_2d_exchange')
    end if
    if ( size(r4vector2) < ssize ) then
      call getmem1d(r4vector2,1,ssize,'real4_2d_exchange')
    end if

    if ( ma%bandflag ) then
      ! set the size of the dummy exchange vector
      ! nex is the width of the exchange stencil
      ! isize is the height of a block in the N-S direction
      ssize = nex*isize

      !********************************************************
      ! (1) Send right boundary to left side of right neighbor
      !********************************************************
      ! loop over the given latitudes and number of exchange blocks
      ! and unravel the right boundary into a vector
      ib = 1
      do i = i1, i2
        do j = 1, nex
          r4vector1(ib) = ml(j2-j+1,i)
          ib = ib + 1
        end do
      end do
      ! do the send-right exchange
      call cyclic_exchange_array(r4vector1,r4vector2, &
                                 ssize,ma%right,ma%left,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the left boundary of this block
      ib = 1
      do i = i1, i2
        do j = 1, nex
          ml(j1-j,i) = r4vector2(ib)
          ib = ib + 1
        end do
      end do

      !********************************************************
      ! (2) Send left boundary to right side of left neighbor
      !********************************************************
      ! loop over the given latitudes and number of exchange blocks
      ! and unravel the left boundary into a vector
      ib = 1
      do i = i1, i2
        do j = 1, nex
          r4vector1(ib) = ml(j1+j-1,i)
          ib = ib + 1
        end do
      end do
      ! do the send-left exchange
      call cyclic_exchange_array(r4vector1,r4vector2, &
                                 ssize,ma%left,ma%right,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the right boundary of this block
      ib = 1
      do i = i1, i2
        do j = 1, nex
          ml(j2+j,i) = r4vector2(ib)
          ib = ib + 1
        end do
      end do

      if ( ma%crmflag ) then
        !*********************************************************
        ! (3) Send bottom boundary to top side of bottom neighbor
        !*********************************************************
        ! set the size of the dummy exchange vector
        ! nex is the width of the exchange stencil
        ! jsize is the width of a block in the E-W direction
        ssize = nex*jsize
        ! loop over the given longitudes and number of exchange blocks
        ! and unravel the bottom boundary into a vector
        ib = 1
        do i = 1, nex
          do j = j1, j2
            r4vector1(ib) = ml(j,i1+i-1)
            ib = ib + 1
          end do
        end do

        ! do the send-up exchange
        call cyclic_exchange_array(r4vector1,r4vector2, &
                                   ssize,ma%bottom,ma%top,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the top boundary of this block
        ib = 1
        do i = 1, nex
          do j = j1, j2
            ml(j,i2+i) = r4vector2(ib)
            ib = ib + 1
          end do
        end do

        !******************************************************
        ! (4) Send top boundary to bottom side of top neighbor
        !******************************************************
        ! loop over the given longitudes and number of exchange blocks
        ! and unravel the top boundary into a vector
        ib = 1
        do i = 1, nex
          do j = j1, j2
            r4vector1(ib) = ml(j,i2-i+1)
            ib = ib + 1
          end do
        end do

        ! do the send-down exchange
        call cyclic_exchange_array(r4vector1,r4vector2, &
                                   ssize,ma%top,ma%bottom,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the bottom boundary of this block
        ib = 1
        do i = 1, nex
          do j = j1, j2
            ml(j,i1-i) = r4vector2(ib)
            ib = ib + 1
          end do
        end do

        !*******************************************************************
        ! (5) Send bottom-right boundary to top-left side of bottom-right
        !     neighbor
        !*******************************************************************
        ! set the size of the exchange vector to the
        ! a square the length/width of the exchange stencil
        ssize = nex*nex
        ! loop over the exchange block and unravel it into a vector
        ib = 1
        do i = 1, nex
          do j = 1, nex
            r4vector1(ib) = ml(j2-j+1,i1+i-1)
            ib = ib + 1
          end do
        end do

        ! do the send down-right exchange
        call cyclic_exchange_array(r4vector1,r4vector2, &
                                   ssize,ma%bottomright,ma%topleft,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the top-left boundary of this block
        ib = 1
        do i = 1, nex
          do j = 1, nex
            ml(j1-j,i2+i) = r4vector2(ib)
            ib = ib + 1
          end do
        end do

        !**********************************************************************
        ! (6) Send top-left boundary to bottom-right side of top-left neighbor
        !**********************************************************************
        ! loop over the exchange block and unravel it into a vector
        ib = 1
        do i = 1, nex
          do j = 1, nex
            r4vector1(ib) = ml(j1+j-1,i2-i+1)
            ib = ib + 1
          end do
        end do

        ! do the send up-left exchange
        call cyclic_exchange_array(r4vector1,r4vector2, &
                                   ssize,ma%topleft,ma%bottomright,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the bottom-right boundary of this block
        ib = 1
        do i = 1, nex
          do j = 1, nex
            ml(j2+j,i1-i) = r4vector2(ib)
            ib = ib + 1
          end do
        end do

        !******************************************************************
        ! (7) Send bottom-left boundary to top-right side of bottom-left
        !     neighbor
        !******************************************************************
        ! loop over the exchange block and unravel it into a vector
        ib = 1
        do i = 1, nex
          do j = 1, nex
            r4vector1(ib) = ml(j1+j-1,i1+i-1)
            ib = ib + 1
          end do
        end do

        ! do the send down-left exchange
        call cyclic_exchange_array(r4vector1,r4vector2, &
                                   ssize,ma%bottomleft,ma%topright,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the top-right boundary of this block
        ib = 1
        do i = 1, nex
          do j = 1, nex
            ml(j2+j,i2+i) = r4vector2(ib)
            ib = ib + 1
          end do
        end do

        !*****************************************************************
        ! (8) Send top-right boundary to bottom-left side of top-right
        !     neighbor
        !******************************************************************
        ! loop over the exchange block and unravel it into a vector
        ib = 1
        do i = 1, nex
          do j = 1, nex
            r4vector1(ib) = ml(j2-j+1,i2-i+1)
            ib = ib + 1
          end do
        end do

        ! do the send up-right exchange
        call cyclic_exchange_array(r4vector1,r4vector2, &
                                   ssize,ma%topright,ma%bottomleft,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the bottom-left boundary of this block
        ib = 1
        do i = 1, nex
          do j = 1, nex
            ml(j1-j,i1-i) = r4vector2(ib)
            ib = ib + 1
          end do
        end do
      end if

    else

      ! RIGHT AND LEFT EXCHANGE

      if ( ma%right /= mpi_proc_null) then
        ssize = nex*isize
        ib = ipos
        do i = i1, i2
          do j = 1, nex
            r4vector1(ib) = ml(j2-j+1,i)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%right,tag_rl,tag_lr, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if
      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*isize
        ib = ipos
        do i = i1, i2
          do j = 1, nex
            r4vector1(ib) = ml(j1+j-1,i)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%left,tag_lr,tag_rl, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if

    end if

    if ( .not. ma%crmflag ) then

      ! BOTTOM AND TOP EXCHANGE

      if ( ma%top /= mpi_proc_null) then
        ssize = nex*jsize
        ib = ipos
        do i = 1, nex
          do j = j1, j2
            r4vector1(ib) = ml(j,i2-i+1)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%top,tag_tb,tag_bt, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if
      if ( ma%bottom /= mpi_proc_null) then
        ssize = nex*jsize
        ib = ipos
        do i = 1, nex
          do j = j1, j2
            r4vector1(ib) = ml(j,i1+i-1)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%bottom,tag_bt,tag_tb, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if

      ! TOP LEFT WITH BOTTOM RIGHT EXCHANGE

      if ( ma%topleft /= mpi_proc_null ) then
        ssize = nex*nex
        ib = ipos
        do i = 1, nex
          do j = 1, nex
            r4vector1(ib) = ml(j1+j-1,i2-i+1)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%topleft,tag_tlbr,tag_brtl, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if
      if ( ma%bottomright /= mpi_proc_null ) then
        ssize = nex*nex
        ib = ipos
        do i = 1, nex
          do j = 1, nex
            r4vector1(ib) = ml(j2-j+1,i1+i-1)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%bottomright,tag_brtl,tag_tlbr, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if

      ! BOTTOM LEFT WITH TOP RIGHT EXCHANGE

      if ( ma%topright /= mpi_proc_null ) then
        ssize = nex*nex
        ib = ipos
        do i = 1, nex
          do j = 1, nex
            r4vector1(ib) = ml(j2-j+1,i2-i+1)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%topright,tag_trbl,tag_bltr, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if
      if ( ma%bottomleft /= mpi_proc_null ) then
        ssize = nex*nex
        ib = ipos
        do i = 1, nex
          do j = 1, nex
            r4vector1(ib) = ml(j1+j-1,i1+i-1)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%bottomleft,tag_bltr,tag_trbl, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if
    end if

    ! Finalize non cyclic comms

    if ( irc /= 0 ) then
      ipos = 1
      call mpi_waitall(16,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      if ( .not. ma%bandflag ) then
        if ( ma%right /= mpi_proc_null) then
          ib = ipos
          ssize = nex*isize
          do i = i1, i2
            do j = 1, nex
              ml(j2+j,i) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%left /= mpi_proc_null ) then
          ib = ipos
          ssize = nex*isize
          do i = i1, i2
            do j = 1, nex
              ml(j1-j,i) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
      if ( .not. ma%crmflag ) then
        if ( ma%top /= mpi_proc_null) then
          ib = ipos
          ssize = nex*jsize
          do i = 1, nex
            do j = j1, j2
              ml(j,i2+i) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%bottom /= mpi_proc_null) then
          ib = ipos
          ssize = nex*jsize
          do i = 1, nex
            do j = j1, j2
              ml(j,i1-i) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%topleft /= mpi_proc_null ) then
          ib = ipos
          ssize = nex*nex
          do i = 1, nex
            do j = 1, nex
              ml(j1-j,i2+i) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%bottomright /= mpi_proc_null ) then
          ib = ipos
          ssize = nex*nex
          do i = 1, nex
            do j = 1, nex
              ml(j2+j,i1-i) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%topright /= mpi_proc_null ) then
          ib = ipos
          ssize =nex*nex
          do i = 1, nex
            do j = 1, nex
              ml(j2+j,i2+i) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%bottomleft /= mpi_proc_null ) then
          ib = ipos
          ssize = nex*nex
          do i = 1, nex
            do j = 1, nex
              ml(j1-j,i1-i) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
    end if
  end subroutine real4_2d_exchange

  subroutine real4_3d_exchange(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2
    integer(ik4) :: isize, jsize, ksize, ssize, j, i, k, ib, irc, ipos
    integer(ik4), dimension(16) :: req

    isize = i2-i1+1
    jsize = j2-j1+1
    ksize = k2-k1+1

    irc = 0
    ipos = 1
    req = mpi_request_null

    ! Total max possible communication size
    ! set the size of the dummy exchange vector
    ! nex is the width of the exchange stencil
    ! isize is the height of a block in the N-S direction
    ! jsize is the height of a block in the W-E direction
    ! ksize is the height of a block in the T-B direction
    ssize = nex * ksize * (4*nex+2*isize+2*jsize)
    if ( size(r4vector1) < ssize ) then
      call getmem1d(r4vector1,1,ssize,'real4_3d_exchange')
    end if
    if ( size(r4vector2) < ssize ) then
      call getmem1d(r4vector2,1,ssize,'real4_3d_exchange')
    end if

    if ( ma%bandflag ) then
      ! set the size of the dummy exchange vector
      ! nex is the width of the exchange stencil
      ! isize is the height of a block in the N-S direction
      ssize = nex*isize*ksize

      !********************************************************
      ! (1) Send right boundary to left side of right neighbor
      !********************************************************
      ! loop over the given latitudes and number of exchange blocks
      ! and unravel the right boundary into a vector
      ib = 1
      do k = k1, k2
        do i = i1, i2
          do j = 1, nex
            r4vector1(ib) = ml(j2-j+1,i,k)
            ib = ib + 1
          end do
        end do
      end do
      ! do the send-right exchange
      call cyclic_exchange_array(r4vector1,r4vector2, &
                                 ssize,ma%right,ma%left,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the left boundary of this block
      ib = 1
      do k = k1, k2
        do i = i1, i2
          do j = 1, nex
            ml(j1-j,i,k) = r4vector2(ib)
            ib = ib + 1
          end do
        end do
      end do

      !********************************************************
      ! (2) Send left boundary to right side of left neighbor
      !********************************************************
      ! loop over the given latitudes and number of exchange blocks
      ! and unravel the left boundary into a vector
      ib = 1
      do k = k1, k2
        do i = i1, i2
          do j = 1, nex
            r4vector1(ib) = ml(j1+j-1,i,k)
            ib = ib + 1
          end do
        end do
      end do
      ! do the send-left exchange
      call cyclic_exchange_array(r4vector1,r4vector2, &
                                 ssize,ma%left,ma%right,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the right boundary of this block
      ib = 1
      do k = k1, k2
        do i = i1, i2
          do j = 1, nex
            ml(j2+j,i,k) = r4vector2(ib)
            ib = ib + 1
          end do
        end do
      end do

      if ( ma%crmflag ) then
        !*********************************************************
        ! (3) Send bottom boundary to top side of bottom neighbor
        !*********************************************************
        ! set the size of the dummy exchange vector
        ! nex is the width of the exchange stencil
        ! jsize is the width of a block in the E-W direction
        ssize = nex*jsize*ksize
        ! loop over the given longitudes and number of exchange blocks
        ! and unravel the bottom boundary into a vector
        ib = 1
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              r4vector1(ib) = ml(j,i1+i-1,k)
              ib = ib + 1
            end do
          end do
        end do

        ! do the send-up exchange
        call cyclic_exchange_array(r4vector1,r4vector2, &
                                   ssize,ma%bottom,ma%top,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the top boundary of this block
        ib = 1
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              ml(j,i2+i,k) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
        end do

        !******************************************************
        ! (4) Send top boundary to bottom side of top neighbor
        !******************************************************
        ! loop over the given longitudes and number of exchange blocks
        ! and unravel the top boundary into a vector
        ib = 1
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              r4vector1(ib) = ml(j,i2-i+1,k)
              ib = ib + 1
            end do
          end do
        end do

        ! do the send-down exchange
        call cyclic_exchange_array(r4vector1,r4vector2, &
                                   ssize,ma%top,ma%bottom,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the bottom boundary of this block
        ib = 1
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              ml(j,i1-i,k) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
        end do

        !*******************************************************************
        ! (5) Send bottom-right boundary to top-left side of bottom-right
        !     neighbor
        !*******************************************************************
        ! set the size of the exchange vector to the
        ! a square the length/width of the exchange stencil
        ssize = nex*nex*ksize
        ! loop over the exchange block and unravel it into a vector
        ib = 1
        do k = k1, k2
          do i = 1, nex
            do j = 1, nex
              r4vector1(ib) = ml(j2-j+1,i1+i-1,k)
              ib = ib + 1
            end do
          end do
        end do

        ! do the send down-right exchange
        call cyclic_exchange_array(r4vector1,r4vector2, &
                                   ssize,ma%bottomright,ma%topleft,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the top-left boundary of this block
        ib = 1
        do k = k1, k2
          do i = 1, nex
            do j = 1, nex
              ml(j1-j,i2+i,k) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
        end do

        !**********************************************************************
        ! (6) Send top-left boundary to bottom-right side of top-left neighbor
        !**********************************************************************
        ! loop over the exchange block and unravel it into a vector
        ib = 1
        do k = k1, k2
          do i = 1, nex
            do j = 1, nex
              r4vector1(ib) = ml(j1+j-1,i2-i+1,k)
              ib = ib + 1
            end do
          end do
        end do

        ! do the send up-left exchange
        call cyclic_exchange_array(r4vector1,r4vector2, &
                                   ssize,ma%topleft,ma%bottomright,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the bottom-right boundary of this block
        ib = 1
        do k = k1, k2
          do i = 1, nex
            do j = 1, nex
              ml(j2+j,i1-i,k) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
        end do

        !******************************************************************
        ! (7) Send bottom-left boundary to top-right side of bottom-left
        !     neighbor
        !******************************************************************
        ! loop over the exchange block and unravel it into a vector
        ib = 1
        do k = k1, k2
          do i = 1, nex
            do j = 1, nex
              r4vector1(ib) = ml(j1+j-1,i1+i-1,k)
              ib = ib + 1
            end do
          end do
        end do

        ! do the send down-left exchange
        call cyclic_exchange_array(r4vector1,r4vector2, &
                                   ssize,ma%bottomleft,ma%topright,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the top-right boundary of this block
        ib = 1
        do k = k1, k2
          do i = 1, nex
            do j = 1, nex
              ml(j2+j,i2+i,k) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
        end do

        !*****************************************************************
        ! (8) Send top-right boundary to bottom-left side of top-right
        !     neighbor
        !******************************************************************
        ! loop over the exchange block and unravel it into a vector
        ib = 1
        do k = k1, k2
          do i = 1, nex
            do j = 1, nex
              r4vector1(ib) = ml(j2-j+1,i2-i+1,k)
              ib = ib + 1
            end do
          end do
        end do

        ! do the send up-right exchange
        call cyclic_exchange_array(r4vector1,r4vector2, &
                                   ssize,ma%topright,ma%bottomleft,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the bottom-left boundary of this block
        ib = 1
        do k = k1, k2
          do i = 1, nex
            do j = 1, nex
              ml(j1-j,i1-i,k) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end if

    else

      ! RIGHT AND LEFT EXCHANGE

      if ( ma%right /= mpi_proc_null) then
        ssize = nex*isize*ksize
        ib = ipos
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              r4vector1(ib) = ml(j2-j+1,i,k)
              ib = ib + 1
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%right,tag_rl,tag_lr, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if
      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*isize*ksize
        ib = ipos
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              r4vector1(ib) = ml(j1+j-1,i,k)
              ib = ib + 1
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%left,tag_lr,tag_rl, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if

    end if

    if ( .not. ma%crmflag ) then

      ! BOTTOM AND TOP EXCHANGE

      if ( ma%top /= mpi_proc_null) then
        ssize = nex*jsize*ksize
        ib = ipos
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              r4vector1(ib) = ml(j,i2-i+1,k)
              ib = ib + 1
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%top,tag_tb,tag_bt, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if
      if ( ma%bottom /= mpi_proc_null) then
        ssize = nex*jsize*ksize
        ib = ipos
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              r4vector1(ib) = ml(j,i1+i-1,k)
              ib = ib + 1
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%bottom,tag_bt,tag_tb, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if

      ! TOP LEFT WITH BOTTOM RIGHT EXCHANGE

      if ( ma%topleft /= mpi_proc_null ) then
        ssize = nex*nex*ksize
        ib = ipos
        do k = k1, k2
          do i = 1, nex
            do j = 1, nex
              r4vector1(ib) = ml(j1+j-1,i2-i+1,k)
              ib = ib + 1
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%topleft,tag_tlbr,tag_brtl, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if
      if ( ma%bottomright /= mpi_proc_null ) then
        ssize = nex*nex*ksize
        ib = ipos
        do k = k1, k2
          do i = 1, nex
            do j = 1, nex
              r4vector1(ib) = ml(j2-j+1,i1+i-1,k)
              ib = ib + 1
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%bottomright,tag_brtl,tag_tlbr, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if

      ! BOTTOM LEFT WITH TOP RIGHT EXCHANGE

      if ( ma%topright /= mpi_proc_null ) then
        ssize = nex*nex*ksize
        ib = ipos
        do k = k1, k2
          do i = 1, nex
            do j = 1, nex
              r4vector1(ib) = ml(j2-j+1,i2-i+1,k)
              ib = ib + 1
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%topright,tag_trbl,tag_bltr, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if
      if ( ma%bottomleft /= mpi_proc_null ) then
        ssize = nex*nex*ksize
        ib = ipos
        do k = k1, k2
          do i = 1, nex
            do j = 1, nex
              r4vector1(ib) = ml(j1+j-1,i1+i-1,k)
              ib = ib + 1
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%bottomleft,tag_bltr,tag_trbl, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if
    end if

    ! Finalize non cyclic comms

    if ( irc /= 0 ) then
      ipos = 1
      call mpi_waitall(16,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      if ( .not. ma%bandflag ) then
        if ( ma%right /= mpi_proc_null) then
          ib = ipos
          ssize = nex*isize*ksize
          do k = k1, k2
            do i = i1, i2
              do j = 1, nex
                ml(j2+j,i,k) = r4vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%left /= mpi_proc_null ) then
          ib = ipos
          ssize = nex*isize*ksize
          do k = k1, k2
            do i = i1, i2
              do j = 1, nex
                ml(j1-j,i,k) = r4vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
      if ( .not. ma%crmflag ) then
        if ( ma%top /= mpi_proc_null) then
          ib = ipos
          ssize = nex*jsize*ksize
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                ml(j,i2+i,k) = r4vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%bottom /= mpi_proc_null) then
          ib = ipos
          ssize = nex*jsize*ksize
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                ml(j,i1-i,k) = r4vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%topleft /= mpi_proc_null ) then
          ib = ipos
          ssize = nex*nex*ksize
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                ml(j1-j,i2+i,k) = r4vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%bottomright /= mpi_proc_null ) then
          ib = ipos
          ssize = nex*nex*ksize
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                ml(j2+j,i1-i,k) = r4vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%topright /= mpi_proc_null ) then
          ib = ipos
          ssize =nex*nex*ksize
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                ml(j2+j,i2+i,k) = r4vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%bottomleft /= mpi_proc_null ) then
          ib = ipos
          ssize = nex*nex*ksize
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                ml(j1-j,i1-i,k) = r4vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
    end if
  end subroutine real4_3d_exchange

  subroutine real4_4d_exchange(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2, n1, n2
    integer(ik4) :: isize, jsize, ksize, nsize, ssize
    integer(ik4) :: j, i, k, n, ib, irc, ipos
    integer(ik4), dimension(16) :: req

    isize = i2-i1+1
    jsize = j2-j1+1
    ksize = k2-k1+1
    nsize = n2-n1+1

    irc = 0
    ipos = 1
    req = mpi_request_null

    ! Total max possible communication size
    ! set the size of the dummy exchange vector
    ! nex is the width of the exchange stencil
    ! isize is the height of a block in the N-S direction
    ! jsize is the height of a block in the W-E direction
    ! ksize is the height of a block in the T-B direction
    ! nsize is the height of a block in the tracer direction
    ssize = nex * ksize * nsize * (4*nex+2*isize+2*jsize)
    if ( size(r4vector1) < ssize ) then
      call getmem1d(r4vector1,1,ssize,'real4_4d_exchange')
    end if
    if ( size(r4vector2) < ssize ) then
      call getmem1d(r4vector2,1,ssize,'real4_4d_exchange')
    end if

    if ( ma%bandflag ) then
      ! set the size of the dummy exchange vector
      ! nex is the width of the exchange stencil
      ! isize is the height of a block in the N-S direction
      ssize = nex*isize*ksize*nsize

      !********************************************************
      ! (1) Send right boundary to left side of right neighbor
      !********************************************************
      ! loop over the given latitudes and number of exchange blocks
      ! and unravel the right boundary into a vector
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              r4vector1(ib) = ml(j2-j+1,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      ! do the send-right exchange
      call cyclic_exchange_array(r4vector1,r4vector2, &
                                 ssize,ma%right,ma%left,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the left boundary of this block
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              ml(j1-j,i,k,n) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do

      !********************************************************
      ! (2) Send left boundary to right side of left neighbor
      !********************************************************
      ! loop over the given latitudes and number of exchange blocks
      ! and unravel the left boundary into a vector
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              r4vector1(ib) = ml(j1+j-1,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      ! do the send-left exchange
      call cyclic_exchange_array(r4vector1,r4vector2, &
                                 ssize,ma%left,ma%right,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the right boundary of this block
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              ml(j2+j,i,k,n) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do

      if ( ma%crmflag ) then
        !*********************************************************
        ! (3) Send bottom boundary to top side of bottom neighbor
        !*********************************************************
        ! set the size of the dummy exchange vector
        ! nex is the width of the exchange stencil
        ! jsize is the width of a block in the E-W direction
        ssize = nex*jsize*ksize*nsize
        ! loop over the given longitudes and number of exchange blocks
        ! and unravel the bottom boundary into a vector
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                r4vector1(ib) = ml(j,i1+i-1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do

        ! do the send-up exchange
        call cyclic_exchange_array(r4vector1,r4vector2, &
                                   ssize,ma%bottom,ma%top,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the top boundary of this block
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                ml(j,i2+i,k,n) = r4vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do

        !******************************************************
        ! (4) Send top boundary to bottom side of top neighbor
        !******************************************************
        ! loop over the given longitudes and number of exchange blocks
        ! and unravel the top boundary into a vector
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                r4vector1(ib) = ml(j,i2-i+1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do

        ! do the send-down exchange
        call cyclic_exchange_array(r4vector1,r4vector2, &
                                   ssize,ma%top,ma%bottom,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the bottom boundary of this block
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                ml(j,i1-i,k,n) = r4vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do

        !*******************************************************************
        ! (5) Send bottom-right boundary to top-left side of bottom-right
        !     neighbor
        !*******************************************************************
        ! set the size of the exchange vector to the
        ! a square the length/width of the exchange stencil
        ssize = nex*nex*ksize*nsize
        ! loop over the exchange block and unravel it into a vector
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                r4vector1(ib) = ml(j2-j+1,i1+i-1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do

        ! do the send down-right exchange
        call cyclic_exchange_array(r4vector1,r4vector2, &
                                   ssize,ma%bottomright,ma%topleft,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the top-left boundary of this block
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                ml(j1-j,i2+i,k,n) = r4vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do

        !**********************************************************************
        ! (6) Send top-left boundary to bottom-right side of top-left neighbor
        !**********************************************************************
        ! loop over the exchange block and unravel it into a vector
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                r4vector1(ib) = ml(j1+j-1,i2-i+1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do

        ! do the send up-left exchange
        call cyclic_exchange_array(r4vector1,r4vector2, &
                                   ssize,ma%topleft,ma%bottomright,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the bottom-right boundary of this block
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                ml(j2+j,i1-i,k,n) = r4vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do

        !******************************************************************
        ! (7) Send bottom-left boundary to top-right side of bottom-left
        !     neighbor
        !******************************************************************
        ! loop over the exchange block and unravel it into a vector
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                r4vector1(ib) = ml(j1+j-1,i1+i-1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do

        ! do the send down-left exchange
        call cyclic_exchange_array(r4vector1,r4vector2, &
                                   ssize,ma%bottomleft,ma%topright,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the top-right boundary of this block
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                ml(j2+j,i2+i,k,n) = r4vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do

        !*****************************************************************
        ! (8) Send top-right boundary to bottom-left side of top-right
        !     neighbor
        !******************************************************************
        ! loop over the exchange block and unravel it into a vector
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                r4vector1(ib) = ml(j2-j+1,i2-i+1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do

        ! do the send up-right exchange
        call cyclic_exchange_array(r4vector1,r4vector2, &
                                   ssize,ma%topright,ma%bottomleft,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the bottom-left boundary of this block
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                ml(j1-j,i1-i,k,n) = r4vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do
      end if

    else

      ! RIGHT AND LEFT EXCHANGE

      if ( ma%right /= mpi_proc_null) then
        ssize = nex*isize*ksize*nsize
        ib = ipos
        do n = n1, n2
          do k = k1, k2
            do i = i1, i2
              do j = 1, nex
                r4vector1(ib) = ml(j2-j+1,i,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%right,tag_rl,tag_lr, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if
      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*isize*ksize*nsize
        ib = ipos
        do n = n1, n2
          do k = k1, k2
            do i = i1, i2
              do j = 1, nex
                r4vector1(ib) = ml(j1+j-1,i,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%left,tag_lr,tag_rl, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if

    end if

    if ( .not. ma%crmflag ) then

      ! BOTTOM AND TOP EXCHANGE

      if ( ma%top /= mpi_proc_null) then
        ssize = nex*jsize*ksize*nsize
        ib = ipos
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                r4vector1(ib) = ml(j,i2-i+1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%top,tag_tb,tag_bt, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if
      if ( ma%bottom /= mpi_proc_null) then
        ssize = nex*jsize*ksize*nsize
        ib = ipos
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                r4vector1(ib) = ml(j,i1+i-1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%bottom,tag_bt,tag_tb, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if

      ! TOP LEFT WITH BOTTOM RIGHT EXCHANGE

      if ( ma%topleft /= mpi_proc_null ) then
        ssize = nex*nex*ksize*nsize
        ib = ipos
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                r4vector1(ib) = ml(j1+j-1,i2-i+1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%topleft,tag_tlbr,tag_brtl, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if
      if ( ma%bottomright /= mpi_proc_null ) then
        ssize = nex*nex*ksize*nsize
        ib = ipos
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                r4vector1(ib) = ml(j2-j+1,i1+i-1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%bottomright,tag_brtl,tag_tlbr, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if

      ! BOTTOM LEFT WITH TOP RIGHT EXCHANGE

      if ( ma%topright /= mpi_proc_null ) then
        ssize = nex*nex*ksize*nsize
        ib = ipos
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                r4vector1(ib) = ml(j2-j+1,i2-i+1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%topright,tag_trbl,tag_bltr, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if
      if ( ma%bottomleft /= mpi_proc_null ) then
        ssize = nex*nex*ksize*nsize
        ib = ipos
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                r4vector1(ib) = ml(j1+j-1,i1+i-1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%bottomleft,tag_bltr,tag_trbl, &
                            req(irc),req(irc+8))
        ipos = ipos + ssize
      end if
    end if

    ! Finalize non cyclic comms

    if ( irc /= 0 ) then
      ipos = 1
      call mpi_waitall(16,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      if ( .not. ma%bandflag ) then
        if ( ma%right /= mpi_proc_null) then
          ib = ipos
          ssize = nex*isize*ksize*nsize
          do n = n1, n2
            do k = k1, k2
              do i = i1, i2
                do j = 1, nex
                  ml(j2+j,i,k,n) = r4vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%left /= mpi_proc_null ) then
          ib = ipos
          ssize = nex*isize*ksize*nsize
          do n = n1, n2
            do k = k1, k2
              do i = i1, i2
                do j = 1, nex
                  ml(j1-j,i,k,n) = r4vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
      if ( .not. ma%crmflag ) then
        if ( ma%top /= mpi_proc_null) then
          ib = ipos
          ssize = nex*jsize*ksize*nsize
          do n = n1, n2
            do k = k1, k2
              do i = 1, nex
                do j = j1, j2
                  ml(j,i2+i,k,n) = r4vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%bottom /= mpi_proc_null) then
          ib = ipos
          ssize = nex*jsize*ksize*nsize
          do n = n1, n2
            do k = k1, k2
              do i = 1, nex
                do j = j1, j2
                  ml(j,i1-i,k,n) = r4vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%topleft /= mpi_proc_null ) then
          ib = ipos
          ssize = nex*nex*ksize*nsize
          do n = n1, n2
            do k = k1, k2
             do i = 1, nex
                do j = 1, nex
                  ml(j1-j,i2+i,k,n) = r4vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%bottomright /= mpi_proc_null ) then
          ib = ipos
          ssize = nex*nex*ksize*nsize
          do n = n1, n2
            do k = k1, k2
              do i = 1, nex
                do j = 1, nex
                  ml(j2+j,i1-i,k,n) = r4vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%topright /= mpi_proc_null ) then
          ib = ipos
          ssize =nex*nex*ksize*nsize
          do n = n1, n2
            do k = k1, k2
              do i = 1, nex
                do j = 1, nex
                  ml(j2+j,i2+i,k,n) = r4vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%bottomleft /= mpi_proc_null ) then
          ib = ipos
          ssize = nex*nex*ksize*nsize
          do n = n1, n2
            do k = k1, k2
              do i = 1, nex
                do j = 1, nex
                  ml(j1-j,i1-i,k,n) = r4vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
    end if
  end subroutine real4_4d_exchange

  subroutine real8_2d_exchange_left_right_bottom_top(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2
    integer(ik4) :: isize, jsize, ssize, j, i, ib, irc, ipos
    integer(ik4), dimension(8) :: req

    isize = i2-i1+1
    jsize = j2-j1+1

    irc = 0
    ipos = 1
    req = mpi_request_null

    ! Total max possible communication size
    ! set the size of the dummy exchange vector
    ! nex is the width of the exchange stencil
    ! isize is the height of a block in the N-S direction
    ! jsize is the height of a block in the W-E direction
    ssize = nex * (2*isize+2*jsize)
    if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'real8_2d_exchange_left_right_bottom_top')
    end if
    if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'real8_2d_exchange_left_right_bottom_top')
    end if

    if ( ma%bandflag ) then
      ssize = nex*isize
      !********************************************************
      ! (1) Send right boundary to left side of right neighbor
      !********************************************************
      ! loop over the given latitudes and number of exchange blocks
      ! and unravel the right boundary into a vector
      do i = i1, i2
        do j = 1, nex
          ib = j + (i - i1) * nex
          r8vector1(ib) = ml(j2-j+1,i)
        end do
      end do
      ! do the send-right exchange
      call cyclic_exchange_array(r8vector1,r8vector2, &
                                 ssize,ma%right,ma%left,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the left boundary of this block
      do i = i1, i2
        do j = 1, nex
          ib = j + (i - i1) * nex
          ml(j1-j,i) = r8vector2(ib)
        end do
      end do

      !********************************************************
      ! (2) Send left boundary to right side of left neighbor
      !********************************************************
      ! loop over the given latitudes and number of exchange blocks
      ! and unravel the left boundary into a vector
      do i = i1, i2
        do j = 1, nex
          ib = j + (i - i1) * nex
          r8vector1(ib) = ml(j1+j-1,i)
        end do
      end do
      ! do the send-left exchange
      call cyclic_exchange_array(r8vector1,r8vector2, &
                                 ssize,ma%left,ma%right,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the right boundary of this block
      do i = i1, i2
        do j = 1, nex
          ib = j + (i - i1) * nex
          ml(j2+j,i) = r8vector2(ib)
        end do
      end do

      if ( ma%crmflag ) then
        !*********************************************************
        ! (3) Send bottom boundary to top side of bottom neighbor
        !*********************************************************
        ssize = nex*jsize
        ! loop over the given longitudes and number of exchange blocks
        ! and unravel the bottom boundary into a vector
        do i = 1, nex
          do j = j1, j2
            ib = (j - j1 + 1) + (i - 1) * (j2 - j1 + 1)
            r8vector1(ib) = ml(j,i1+i-1)
          end do
        end do

        ! do the send-up exchange
        call cyclic_exchange_array(r8vector1,r8vector2, &
                                   ssize,ma%bottom,ma%top,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the top boundary of this block
        do i = 1, nex
          do j = j1, j2
            ib = (j - j1 + 1) + (i - 1) * (j2 - j1 + 1)
            ml(j,i2+i) = r8vector2(ib)
          end do
        end do

        !******************************************************
        ! (4) Send top boundary to bottom side of top neighbor
        !******************************************************
        ! loop over the given longitudes and number of exchange blocks
        ! and unravel the top boundary into a vector
        do i = 1, nex
          do j = j1, j2
            ib = (j - j1 + 1) + (i - 1) * (j2 - j1 + 1)
            r8vector1(ib) = ml(j,i2-i+1)
          end do
        end do

        ! do the send-down exchange
        call cyclic_exchange_array(r8vector1,r8vector2, &
                                   ssize,ma%top,ma%bottom,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the bottom boundary of this block
        do i = 1, nex
          do j = j1, j2
            ib = (j - j1 + 1) + (i - 1) * (j2 - j1 + 1)
            ml(j,i1-i) = r8vector2(ib)
          end do
        end do

      end if

    else

      ! RIGHT AND LEFT EXCHANGE

      if ( ma%right /= mpi_proc_null) then
        ssize = nex*isize
        do i = i1, i2
          do j = 1, nex
            ib = j + (i - i1) * nex + ipos - 1
            r8vector1(ib) = ml(j2-j+1,i)
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%right,tag_rl,tag_lr, &
                            req(irc),req(irc+4))
        ipos = ipos + ssize
      end if
      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*isize
        do i = i1, i2
          do j = 1, nex
            ib = j + (i - i1) * nex + ipos - 1
            r8vector1(ib) = ml(j1+j-1,i)
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%left,tag_lr,tag_rl, &
                            req(irc),req(irc+4))
        ipos = ipos + ssize
      end if

    end if

    if ( .not. ma%crmflag ) then

      ! BOTTOM AND TOP EXCHANGE

      if ( ma%top /= mpi_proc_null) then
        ssize = nex*jsize
        do i = 1, nex
          do j = j1, j2
            ib = (j - j1 + 1) + (i - 1) * (j2 - j1 + 1) + ipos - 1
            r8vector1(ib) = ml(j,i2-i+1)
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%top,tag_tb,tag_bt, &
                            req(irc),req(irc+4))
        ipos = ipos + ssize
      end if
      if ( ma%bottom /= mpi_proc_null) then
        ssize = nex*jsize
        do i = 1, nex
          do j = j1, j2
            ib = (j - j1 + 1) + (i - 1) * (j2 - j1 + 1) + ipos - 1
            r8vector1(ib) = ml(j,i1+i-1)
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%bottom,tag_bt,tag_tb, &
                            req(irc),req(irc+4))
        ipos = ipos + ssize
      end if

    end if

    ! Finalize non cyclic comms

    if ( irc /= 0 ) then
      ipos = 1
      call mpi_waitall(8,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      if ( .not. ma%bandflag ) then
        if ( ma%right /= mpi_proc_null) then
          ssize = nex*isize
          do i = i1, i2
            do j = 1, nex
              ib = j + (i - i1) * nex + ipos - 1
              ml(j2+j,i) = r8vector2(ib)
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%left /= mpi_proc_null ) then
          ssize = nex*isize
          do i = i1, i2
            do j = 1, nex
              ib = j + (i - i1) * nex + ipos - 1
              ml(j1-j,i) = r8vector2(ib)
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
      if ( .not. ma%crmflag ) then
        if ( ma%top /= mpi_proc_null) then
          ssize = nex*jsize
          do i = 1, nex
            do j = j1, j2
              ib = (j - j1 + 1) + (i - 1) * (j2 - j1 + 1) + ipos - 1
              ml(j,i2+i) = r8vector2(ib)
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%bottom /= mpi_proc_null) then
          ssize = nex*jsize
          do i = 1, nex
            do j = j1, j2
              ib = (j - j1 + 1) + (i - 1) * (j2 - j1 + 1) + ipos - 1
              ml(j,i1-i) = r8vector2(ib)
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
    end if
  end subroutine real8_2d_exchange_left_right_bottom_top

  subroutine real8_3d_exchange_left_right_bottom_top(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2
    integer(ik4) :: isize, jsize, ksize, ssize, j, i, k, ib, irc, ipos
    integer(ik4), dimension(8) :: req

    isize = i2-i1+1
    jsize = j2-j1+1
    ksize = k2-k1+1

    irc = 0
    ipos = 1
    req = mpi_request_null

    ! Total max possible communication size
    ! set the size of the dummy exchange vector
    ! nex is the width of the exchange stencil
    ! isize is the height of a block in the N-S direction
    ! jsize is the height of a block in the W-E direction
    ! ksize is the height of a block in the T-B direction
    ssize = nex * ksize * (2*isize+2*jsize)
    if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'real8_3d_exchange_left_right_bottom_top')
    end if
    if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'real8_3d_exchange_left_right_bottom_top')
    end if

    if ( ma%bandflag ) then
      ! set the size of the dummy exchange vector
      ! nex is the width of the exchange stencil
      ! isize is the height of a block in the N-S direction
      ssize = nex*isize*ksize

      !********************************************************
      ! (1) Send right boundary to left side of right neighbor
      !********************************************************
      ! loop over the given latitudes and number of exchange blocks
      ! and unravel the right boundary into a vector
      do k = k1, k2
        do i = i1, i2
          do j = 1, nex
            ib = j + (i - i1) * nex + (k - k1) * (i2 - i1 + 1) * nex
            r8vector1(ib) = ml(j2-j+1,i,k)
          end do
        end do
      end do
      ! do the send-right exchange
      call cyclic_exchange_array(r8vector1,r8vector2, &
                                 ssize,ma%right,ma%left,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the left boundary of this block
      do k = k1, k2
        do i = i1, i2
          do j = 1, nex
            ib = j + (i - i1) * nex + (k - k1) * (i2 - i1 + 1) * nex
            ml(j1-j,i,k) = r8vector2(ib)
          end do
        end do
      end do

      !********************************************************
      ! (2) Send left boundary to right side of left neighbor
      !********************************************************
      ! loop over the given latitudes and number of exchange blocks
      ! and unravel the left boundary into a vector
      do k = k1, k2
        do i = i1, i2
          do j = 1, nex
            ib = j + (i - i1) * nex + (k - k1) * (i2 - i1 + 1) * nex
            r8vector1(ib) = ml(j1+j-1,i,k)
          end do
        end do
      end do
      ! do the send-left exchange
      call cyclic_exchange_array(r8vector1,r8vector2, &
                                 ssize,ma%left,ma%right,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the right boundary of this block
      do k = k1, k2
        do i = i1, i2
          do j = 1, nex
            ib = j + (i - i1) * nex + (k - k1) * (i2 - i1 + 1) * nex
            ml(j2+j,i,k) = r8vector2(ib)
          end do
        end do
      end do

      if ( ma%crmflag ) then
        !*********************************************************
        ! (3) Send bottom boundary to top side of bottom neighbor
        !*********************************************************
        ! set the size of the dummy exchange vector
        ! nex is the width of the exchange stencil
        ! jsize is the width of a block in the E-W direction
        ssize = nex*jsize*ksize
        ! loop over the given longitudes and number of exchange blocks
        ! and unravel the bottom boundary into a vector
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              ib = (j - j1 + 1) + (i - 1) * (j2 - j1 + 1) + &
                (k - k1) * (j2 - j1 + 1) * nex
              r8vector1(ib) = ml(j,i1+i-1,k)
            end do
          end do
        end do

        ! do the send-up exchange
        call cyclic_exchange_array(r8vector1,r8vector2, &
                                   ssize,ma%bottom,ma%top,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the top boundary of this block
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              ib = (j - j1 + 1) + (i - 1) * (j2 - j1 + 1) + &
                (k - k1) * (j2 - j1 + 1) * nex
              ml(j,i2+i,k) = r8vector2(ib)
            end do
          end do
        end do

        !******************************************************
        ! (4) Send top boundary to bottom side of top neighbor
        !******************************************************
        ! loop over the given longitudes and number of exchange blocks
        ! and unravel the top boundary into a vector
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              ib = (j - j1 + 1) + (i - 1) * (j2 - j1 + 1) + &
                (k - k1) * (j2 - j1 + 1) * nex
              r8vector1(ib) = ml(j,i2-i+1,k)
            end do
          end do
        end do

        ! do the send-down exchange
        call cyclic_exchange_array(r8vector1,r8vector2, &
                                   ssize,ma%top,ma%bottom,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the bottom boundary of this block
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              ib = (j - j1 + 1) + (i - 1) * (j2 - j1 + 1) + &
                (k - k1) * (j2 - j1 + 1) * nex
              ml(j,i1-i,k) = r8vector2(ib)
            end do
          end do
        end do
      end if

    else

      ! RIGHT AND LEFT EXCHANGE

      if ( ma%right /= mpi_proc_null) then
        ssize = nex*isize*ksize
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              ib = j + (i - i1) * nex + &
                (k - k1) * (i2 - i1 + 1) * nex + ipos - 1
              r8vector1(ib) = ml(j2-j+1,i,k)
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%right,tag_rl,tag_lr, &
                            req(irc),req(irc+4))
        ipos = ipos + ssize
      end if
      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*isize*ksize
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              ib = j + (i - i1) * nex + &
                (k - k1) * (i2 - i1 + 1) * nex + ipos - 1
              r8vector1(ib) = ml(j1+j-1,i,k)
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%left,tag_lr,tag_rl, &
                            req(irc),req(irc+4))
        ipos = ipos + ssize
      end if

    end if

    if ( .not. ma%crmflag ) then

      ! BOTTOM AND TOP EXCHANGE

      if ( ma%top /= mpi_proc_null) then
        ssize = nex*jsize*ksize
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              ib = (j - j1 + 1) + (i - 1) * (j2 - j1 + 1) + &
                (k - k1) * (j2 - j1 + 1) * nex + ipos - 1
              r8vector1(ib) = ml(j,i2-i+1,k)
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%top,tag_tb,tag_bt, &
                            req(irc),req(irc+4))
        ipos = ipos + ssize
      end if
      if ( ma%bottom /= mpi_proc_null) then
        ssize = nex*jsize*ksize
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              ib = (j - j1 + 1) + (i - 1) * (j2 - j1 + 1) + &
                (k - k1) * (j2 - j1 + 1) * nex + ipos - 1
              r8vector1(ib) = ml(j,i1+i-1,k)
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%bottom,tag_bt,tag_tb, &
                            req(irc),req(irc+4))
        ipos = ipos + ssize
      end if

    end if

    ! Finalize non cyclic comms

    if ( irc /= 0 ) then
      ipos = 1
      call mpi_waitall(8,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      if ( .not. ma%bandflag ) then
        if ( ma%right /= mpi_proc_null) then
          ssize = nex*isize*ksize
          do k = k1, k2
            do i = i1, i2
              do j = 1, nex
                ib = j + (i - i1) * nex + &
                  (k - k1) * (i2 - i1 + 1) * nex + ipos - 1
                ml(j2+j,i,k) = r8vector2(ib)
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%left /= mpi_proc_null ) then
          ssize = nex*isize*ksize
          do k = k1, k2
            do i = i1, i2
              do j = 1, nex
                ib = j + (i - i1) * nex + &
                  (k - k1) * (i2 - i1 + 1) * nex + ipos - 1
                ml(j1-j,i,k) = r8vector2(ib)
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
      if ( .not. ma%crmflag ) then
        if ( ma%top /= mpi_proc_null) then
          ssize = nex*jsize*ksize
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                ib = (j - j1 + 1) + (i - 1) * (j2 - j1 + 1) + &
                  (k - k1) * (j2 - j1 + 1) * nex + ipos - 1
                ml(j,i2+i,k) = r8vector2(ib)
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%bottom /= mpi_proc_null) then
          ssize = nex*jsize*ksize
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                ib = (j - j1 + 1) + (i - 1) * (j2 - j1 + 1) + &
                  (k - k1) * (j2 - j1 + 1) * nex + ipos - 1
                ml(j,i1-i,k) = r8vector2(ib)
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
    end if
  end subroutine real8_3d_exchange_left_right_bottom_top

  subroutine real8_4d_exchange_left_right_bottom_top(ml,nex,j1,j2, &
                                                     i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2, n1, n2
    integer(ik4) :: isize, jsize, ksize, nsize, ssize
    integer(ik4) :: j, i, k, n, ib, irc, ipos
    integer(ik4), dimension(8) :: req

    isize = i2-i1+1
    jsize = j2-j1+1
    ksize = k2-k1+1
    nsize = n2-n1+1

    irc = 0
    ipos = 1
    req = mpi_request_null

    ! Total max possible communication size
    ! set the size of the dummy exchange vector
    ! nex is the width of the exchange stencil
    ! isize is the height of a block in the N-S direction
    ! jsize is the height of a block in the W-E direction
    ! ksize is the height of a block in the T-B direction
    ! nsize is the height of a block in the tracer direction
    ssize = nex * ksize * nsize * (2*isize+2*jsize)
    if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'real8_4d_exchange_left_right_bottom_top')
    end if
    if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'real8_4d_exchange_left_right_bottom_top')
    end if

    if ( ma%bandflag ) then
      ! set the size of the dummy exchange vector
      ! nex is the width of the exchange stencil
      ! isize is the height of a block in the N-S direction
      ssize = nex*isize*ksize*nsize

      !********************************************************
      ! (1) Send right boundary to left side of right neighbor
      !********************************************************
      ! loop over the given latitudes and number of exchange blocks
      ! and unravel the right boundary into a vector
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              r8vector1(ib) = ml(j2-j+1,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      ! do the send-right exchange
      call cyclic_exchange_array(r8vector1,r8vector2, &
                                 ssize,ma%right,ma%left,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the left boundary of this block
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              ml(j1-j,i,k,n) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do

      !********************************************************
      ! (2) Send left boundary to right side of left neighbor
      !********************************************************
      ! loop over the given latitudes and number of exchange blocks
      ! and unravel the left boundary into a vector
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              r8vector1(ib) = ml(j1+j-1,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      ! do the send-left exchange
      call cyclic_exchange_array(r8vector1,r8vector2, &
                                 ssize,ma%left,ma%right,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the right boundary of this block
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              ml(j2+j,i,k,n) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do

      if ( ma%crmflag ) then
        !*********************************************************
        ! (3) Send bottom boundary to top side of bottom neighbor
        !*********************************************************
        ! set the size of the dummy exchange vector
        ! nex is the width of the exchange stencil
        ! jsize is the width of a block in the E-W direction
        ssize = nex*jsize*ksize*nsize
        ! loop over the given longitudes and number of exchange blocks
        ! and unravel the bottom boundary into a vector
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                r8vector1(ib) = ml(j,i1+i-1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do

        ! do the send-up exchange
        call cyclic_exchange_array(r8vector1,r8vector2, &
                                   ssize,ma%bottom,ma%top,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the top boundary of this block
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                ml(j,i2+i,k,n) = r8vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do

        !******************************************************
        ! (4) Send top boundary to bottom side of top neighbor
        !******************************************************
        ! loop over the given longitudes and number of exchange blocks
        ! and unravel the top boundary into a vector
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                r8vector1(ib) = ml(j,i2-i+1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do

        ! do the send-down exchange
        call cyclic_exchange_array(r8vector1,r8vector2, &
                                   ssize,ma%top,ma%bottom,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the bottom boundary of this block
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                ml(j,i1-i,k,n) = r8vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do
      end if

    else

      ! RIGHT AND LEFT EXCHANGE

      if ( ma%right /= mpi_proc_null) then
        ssize = nex*isize*ksize*nsize
        ib = ipos
        do n = n1, n2
          do k = k1, k2
            do i = i1, i2
              do j = 1, nex
                r8vector1(ib) = ml(j2-j+1,i,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%right,tag_rl,tag_lr, &
                            req(irc),req(irc+4))
        ipos = ipos + ssize
      end if
      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*isize*ksize*nsize
        ib = ipos
        do n = n1, n2
          do k = k1, k2
            do i = i1, i2
              do j = 1, nex
                r8vector1(ib) = ml(j1+j-1,i,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%left,tag_lr,tag_rl, &
                            req(irc),req(irc+4))
        ipos = ipos + ssize
      end if

    end if

    if ( .not. ma%crmflag ) then

      ! BOTTOM AND TOP EXCHANGE

      if ( ma%top /= mpi_proc_null) then
        ssize = nex*jsize*ksize*nsize
        ib = ipos
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                r8vector1(ib) = ml(j,i2-i+1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%top,tag_tb,tag_bt, &
                            req(irc),req(irc+4))
        ipos = ipos + ssize
      end if
      if ( ma%bottom /= mpi_proc_null) then
        ssize = nex*jsize*ksize*nsize
        ib = ipos
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                r8vector1(ib) = ml(j,i1+i-1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%bottom,tag_bt,tag_tb, &
                            req(irc),req(irc+4))
        ipos = ipos + ssize
      end if

    end if

    ! Finalize non cyclic comms

    if ( irc /= 0 ) then
      ipos = 1
      call mpi_waitall(8,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      if ( .not. ma%bandflag ) then
        if ( ma%right /= mpi_proc_null) then
          ib = ipos
          ssize = nex*isize*ksize*nsize
          do n = n1, n2
            do k = k1, k2
              do i = i1, i2
                do j = 1, nex
                  ml(j2+j,i,k,n) = r8vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%left /= mpi_proc_null ) then
          ib = ipos
          ssize = nex*isize*ksize*nsize
          do n = n1, n2
            do k = k1, k2
              do i = i1, i2
                do j = 1, nex
                  ml(j1-j,i,k,n) = r8vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
      if ( .not. ma%crmflag ) then
        if ( ma%top /= mpi_proc_null) then
          ib = ipos
          ssize = nex*jsize*ksize*nsize
          do n = n1, n2
            do k = k1, k2
              do i = 1, nex
                do j = j1, j2
                  ml(j,i2+i,k,n) = r8vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%bottom /= mpi_proc_null) then
          ib = ipos
          ssize = nex*jsize*ksize*nsize
          do n = n1, n2
            do k = k1, k2
              do i = 1, nex
                do j = j1, j2
                  ml(j,i1-i,k,n) = r8vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
    end if
  end subroutine real8_4d_exchange_left_right_bottom_top

  subroutine real4_2d_exchange_left_right_bottom_top(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2
    integer(ik4) :: isize, jsize, ssize, j, i, ib, irc, ipos
    integer(ik4), dimension(8) :: req

    isize = i2-i1+1
    jsize = j2-j1+1

    irc = 0
    ipos = 1
    req = mpi_request_null

    ! Total max possible communication size
    ! set the size of the dummy exchange vector
    ! nex is the width of the exchange stencil
    ! isize is the height of a block in the N-S direction
    ! jsize is the height of a block in the W-E direction
    ssize = nex * (2*isize+2*jsize)
    if ( size(r4vector1) < ssize ) then
      call getmem1d(r4vector1,1,ssize,'real4_2d_exchange_left_right_bottom_top')
    end if
    if ( size(r4vector2) < ssize ) then
      call getmem1d(r4vector2,1,ssize,'real4_2d_exchange_left_right_bottom_top')
    end if

    if ( ma%bandflag ) then
      ! set the size of the dummy exchange vector
      ! nex is the width of the exchange stencil
      ! isize is the height of a block in the N-S direction
      ssize = nex*isize

      !********************************************************
      ! (1) Send right boundary to left side of right neighbor
      !********************************************************
      ! loop over the given latitudes and number of exchange blocks
      ! and unravel the right boundary into a vector
      ib = 1
      do i = i1, i2
        do j = 1, nex
          r4vector1(ib) = ml(j2-j+1,i)
          ib = ib + 1
        end do
      end do
      ! do the send-right exchange
      call cyclic_exchange_array(r4vector1,r4vector2, &
                                 ssize,ma%right,ma%left,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the left boundary of this block
      ib = 1
      do i = i1, i2
        do j = 1, nex
          ml(j1-j,i) = r4vector2(ib)
          ib = ib + 1
        end do
      end do

      !********************************************************
      ! (2) Send left boundary to right side of left neighbor
      !********************************************************
      ! loop over the given latitudes and number of exchange blocks
      ! and unravel the left boundary into a vector
      ib = 1
      do i = i1, i2
        do j = 1, nex
          r4vector1(ib) = ml(j1+j-1,i)
          ib = ib + 1
        end do
      end do
      ! do the send-left exchange
      call cyclic_exchange_array(r4vector1,r4vector2, &
                                 ssize,ma%left,ma%right,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the right boundary of this block
      ib = 1
      do i = i1, i2
        do j = 1, nex
          ml(j2+j,i) = r4vector2(ib)
          ib = ib + 1
        end do
      end do

      if ( ma%crmflag ) then
        !*********************************************************
        ! (3) Send bottom boundary to top side of bottom neighbor
        !*********************************************************
        ! set the size of the dummy exchange vector
        ! nex is the width of the exchange stencil
        ! jsize is the width of a block in the E-W direction
        ssize = nex*jsize
        ! loop over the given longitudes and number of exchange blocks
        ! and unravel the bottom boundary into a vector
        ib = 1
        do i = 1, nex
          do j = j1, j2
            r4vector1(ib) = ml(j,i1+i-1)
            ib = ib + 1
          end do
        end do

        ! do the send-up exchange
        call cyclic_exchange_array(r4vector1,r4vector2, &
                                   ssize,ma%bottom,ma%top,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the top boundary of this block
        ib = 1
        do i = 1, nex
          do j = j1, j2
            ml(j,i2+i) = r4vector2(ib)
            ib = ib + 1
          end do
        end do

        !******************************************************
        ! (4) Send top boundary to bottom side of top neighbor
        !******************************************************
        ! loop over the given longitudes and number of exchange blocks
        ! and unravel the top boundary into a vector
        ib = 1
        do i = 1, nex
          do j = j1, j2
            r4vector1(ib) = ml(j,i2-i+1)
            ib = ib + 1
          end do
        end do

        ! do the send-down exchange
        call cyclic_exchange_array(r4vector1,r4vector2, &
                                   ssize,ma%top,ma%bottom,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the bottom boundary of this block
        ib = 1
        do i = 1, nex
          do j = j1, j2
            ml(j,i1-i) = r4vector2(ib)
            ib = ib + 1
          end do
        end do

      end if

    else

      ! RIGHT AND LEFT EXCHANGE

      if ( ma%right /= mpi_proc_null) then
        ssize = nex*isize
        ib = ipos
        do i = i1, i2
          do j = 1, nex
            r4vector1(ib) = ml(j2-j+1,i)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%right,tag_rl,tag_lr, &
                            req(irc),req(irc+4))
        ipos = ipos + ssize
      end if
      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*isize
        ib = ipos
        do i = i1, i2
          do j = 1, nex
            r4vector1(ib) = ml(j1+j-1,i)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%left,tag_lr,tag_rl, &
                            req(irc),req(irc+4))
        ipos = ipos + ssize
      end if

    end if

    if ( .not. ma%crmflag ) then

      ! BOTTOM AND TOP EXCHANGE

      if ( ma%top /= mpi_proc_null) then
        ssize = nex*jsize
        ib = ipos
        do i = 1, nex
          do j = j1, j2
            r4vector1(ib) = ml(j,i2-i+1)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%top,tag_tb,tag_bt, &
                            req(irc),req(irc+4))
        ipos = ipos + ssize
      end if
      if ( ma%bottom /= mpi_proc_null) then
        ssize = nex*jsize
        ib = ipos
        do i = 1, nex
          do j = j1, j2
            r4vector1(ib) = ml(j,i1+i-1)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%bottom,tag_bt,tag_tb, &
                            req(irc),req(irc+4))
        ipos = ipos + ssize
      end if

    end if

    ! Finalize non cyclic comms

    if ( irc /= 0 ) then
      ipos = 1
      call mpi_waitall(8,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      if ( .not. ma%bandflag ) then
        if ( ma%right /= mpi_proc_null) then
          ib = ipos
          ssize = nex*isize
          do i = i1, i2
            do j = 1, nex
              ml(j2+j,i) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%left /= mpi_proc_null ) then
          ib = ipos
          ssize = nex*isize
          do i = i1, i2
            do j = 1, nex
              ml(j1-j,i) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
      if ( .not. ma%crmflag ) then
        if ( ma%top /= mpi_proc_null) then
          ib = ipos
          ssize = nex*jsize
          do i = 1, nex
            do j = j1, j2
              ml(j,i2+i) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%bottom /= mpi_proc_null) then
          ib = ipos
          ssize = nex*jsize
          do i = 1, nex
            do j = j1, j2
              ml(j,i1-i) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
    end if
  end subroutine real4_2d_exchange_left_right_bottom_top

  subroutine real4_3d_exchange_left_right_bottom_top(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2
    integer(ik4) :: isize, jsize, ksize, ssize, j, i, k, ib, irc, ipos
    integer(ik4), dimension(8) :: req

    isize = i2-i1+1
    jsize = j2-j1+1
    ksize = k2-k1+1

    irc = 0
    ipos = 1
    req = mpi_request_null

    ! Total max possible communication size
    ! set the size of the dummy exchange vector
    ! nex is the width of the exchange stencil
    ! isize is the height of a block in the N-S direction
    ! jsize is the height of a block in the W-E direction
    ! ksize is the height of a block in the T-B direction
    ssize = nex * ksize * (2*isize+2*jsize)
    if ( size(r4vector1) < ssize ) then
      call getmem1d(r4vector1,1,ssize,'real4_3d_exchange_left_right_bottom_top')
    end if
    if ( size(r4vector2) < ssize ) then
      call getmem1d(r4vector2,1,ssize,'real4_3d_exchange_left_right_bottom_top')
    end if

    if ( ma%bandflag ) then
      ! set the size of the dummy exchange vector
      ! nex is the width of the exchange stencil
      ! isize is the height of a block in the N-S direction
      ssize = nex*isize*ksize

      !********************************************************
      ! (1) Send right boundary to left side of right neighbor
      !********************************************************
      ! loop over the given latitudes and number of exchange blocks
      ! and unravel the right boundary into a vector
      ib = 1
      do k = k1, k2
        do i = i1, i2
          do j = 1, nex
            r4vector1(ib) = ml(j2-j+1,i,k)
            ib = ib + 1
          end do
        end do
      end do
      ! do the send-right exchange
      call cyclic_exchange_array(r4vector1,r4vector2, &
                                 ssize,ma%right,ma%left,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the left boundary of this block
      ib = 1
      do k = k1, k2
        do i = i1, i2
          do j = 1, nex
            ml(j1-j,i,k) = r4vector2(ib)
            ib = ib + 1
          end do
        end do
      end do

      !********************************************************
      ! (2) Send left boundary to right side of left neighbor
      !********************************************************
      ! loop over the given latitudes and number of exchange blocks
      ! and unravel the left boundary into a vector
      ib = 1
      do k = k1, k2
        do i = i1, i2
          do j = 1, nex
            r4vector1(ib) = ml(j1+j-1,i,k)
            ib = ib + 1
          end do
        end do
      end do
      ! do the send-left exchange
      call cyclic_exchange_array(r4vector1,r4vector2, &
                                 ssize,ma%left,ma%right,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the right boundary of this block
      ib = 1
      do k = k1, k2
        do i = i1, i2
          do j = 1, nex
            ml(j2+j,i,k) = r4vector2(ib)
            ib = ib + 1
          end do
        end do
      end do

      if ( ma%crmflag ) then
        !*********************************************************
        ! (3) Send bottom boundary to top side of bottom neighbor
        !*********************************************************
        ! set the size of the dummy exchange vector
        ! nex is the width of the exchange stencil
        ! jsize is the width of a block in the E-W direction
        ssize = nex*jsize*ksize
        ! loop over the given longitudes and number of exchange blocks
        ! and unravel the bottom boundary into a vector
        ib = 1
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              r4vector1(ib) = ml(j,i1+i-1,k)
              ib = ib + 1
            end do
          end do
        end do

        ! do the send-up exchange
        call cyclic_exchange_array(r4vector1,r4vector2, &
                                   ssize,ma%bottom,ma%top,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the top boundary of this block
        ib = 1
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              ml(j,i2+i,k) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
        end do

        !******************************************************
        ! (4) Send top boundary to bottom side of top neighbor
        !******************************************************
        ! loop over the given longitudes and number of exchange blocks
        ! and unravel the top boundary into a vector
        ib = 1
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              r4vector1(ib) = ml(j,i2-i+1,k)
              ib = ib + 1
            end do
          end do
        end do

        ! do the send-down exchange
        call cyclic_exchange_array(r4vector1,r4vector2, &
                                   ssize,ma%top,ma%bottom,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the bottom boundary of this block
        ib = 1
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              ml(j,i1-i,k) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end if

    else

      ! RIGHT AND LEFT EXCHANGE

      if ( ma%right /= mpi_proc_null) then
        ssize = nex*isize*ksize
        ib = ipos
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              r4vector1(ib) = ml(j2-j+1,i,k)
              ib = ib + 1
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%right,tag_rl,tag_lr, &
                            req(irc),req(irc+4))
        ipos = ipos + ssize
      end if
      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*isize*ksize
        ib = ipos
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              r4vector1(ib) = ml(j1+j-1,i,k)
              ib = ib + 1
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%left,tag_lr,tag_rl, &
                            req(irc),req(irc+4))
        ipos = ipos + ssize
      end if

    end if

    if ( .not. ma%crmflag ) then

      ! BOTTOM AND TOP EXCHANGE

      if ( ma%top /= mpi_proc_null) then
        ssize = nex*jsize*ksize
        ib = ipos
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              r4vector1(ib) = ml(j,i2-i+1,k)
              ib = ib + 1
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%top,tag_tb,tag_bt, &
                            req(irc),req(irc+4))
        ipos = ipos + ssize
      end if
      if ( ma%bottom /= mpi_proc_null) then
        ssize = nex*jsize*ksize
        ib = ipos
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              r4vector1(ib) = ml(j,i1+i-1,k)
              ib = ib + 1
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%bottom,tag_bt,tag_tb, &
                            req(irc),req(irc+4))
        ipos = ipos + ssize
      end if

    end if

    ! Finalize non cyclic comms

    if ( irc /= 0 ) then
      ipos = 1
      call mpi_waitall(8,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      if ( .not. ma%bandflag ) then
        if ( ma%right /= mpi_proc_null) then
          ib = ipos
          ssize = nex*isize*ksize
          do k = k1, k2
            do i = i1, i2
              do j = 1, nex
                ml(j2+j,i,k) = r4vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%left /= mpi_proc_null ) then
          ib = ipos
          ssize = nex*isize*ksize
          do k = k1, k2
            do i = i1, i2
              do j = 1, nex
                ml(j1-j,i,k) = r4vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
      if ( .not. ma%crmflag ) then
        if ( ma%top /= mpi_proc_null) then
          ib = ipos
          ssize = nex*jsize*ksize
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                ml(j,i2+i,k) = r4vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%bottom /= mpi_proc_null) then
          ib = ipos
          ssize = nex*jsize*ksize
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                ml(j,i1-i,k) = r4vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
    end if
  end subroutine real4_3d_exchange_left_right_bottom_top

  subroutine real4_4d_exchange_left_right_bottom_top(ml,nex,j1,j2, &
                                                     i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2, n1, n2
    integer(ik4) :: isize, jsize, ksize, nsize, ssize
    integer(ik4) :: j, i, k, n, ib, irc, ipos
    integer(ik4), dimension(8) :: req

    isize = i2-i1+1
    jsize = j2-j1+1
    ksize = k2-k1+1
    nsize = n2-n1+1

    irc = 0
    ipos = 1
    req = mpi_request_null

    ! Total max possible communication size
    ! set the size of the dummy exchange vector
    ! nex is the width of the exchange stencil
    ! isize is the height of a block in the N-S direction
    ! jsize is the height of a block in the W-E direction
    ! ksize is the height of a block in the T-B direction
    ! nsize is the height of a block in the tracer direction
    ssize = nex * ksize * nsize * (2*isize+2*jsize)
    if ( size(r4vector1) < ssize ) then
      call getmem1d(r4vector1,1,ssize,'real4_4d_exchange_left_right_bottom_top')
    end if
    if ( size(r4vector2) < ssize ) then
      call getmem1d(r4vector2,1,ssize,'real4_4d_exchange_left_right_bottom_top')
    end if

    if ( ma%bandflag ) then
      ! set the size of the dummy exchange vector
      ! nex is the width of the exchange stencil
      ! isize is the height of a block in the N-S direction
      ssize = nex*isize*ksize*nsize

      !********************************************************
      ! (1) Send right boundary to left side of right neighbor
      !********************************************************
      ! loop over the given latitudes and number of exchange blocks
      ! and unravel the right boundary into a vector
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              r4vector1(ib) = ml(j2-j+1,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      ! do the send-right exchange
      call cyclic_exchange_array(r4vector1,r4vector2, &
                                 ssize,ma%right,ma%left,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the left boundary of this block
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              ml(j1-j,i,k,n) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do

      !********************************************************
      ! (2) Send left boundary to right side of left neighbor
      !********************************************************
      ! loop over the given latitudes and number of exchange blocks
      ! and unravel the left boundary into a vector
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              r4vector1(ib) = ml(j1+j-1,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      ! do the send-left exchange
      call cyclic_exchange_array(r4vector1,r4vector2, &
                                 ssize,ma%left,ma%right,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the right boundary of this block
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              ml(j2+j,i,k,n) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do

      if ( ma%crmflag ) then
        !*********************************************************
        ! (3) Send bottom boundary to top side of bottom neighbor
        !*********************************************************
        ! set the size of the dummy exchange vector
        ! nex is the width of the exchange stencil
        ! jsize is the width of a block in the E-W direction
        ssize = nex*jsize*ksize*nsize
        ! loop over the given longitudes and number of exchange blocks
        ! and unravel the bottom boundary into a vector
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                r4vector1(ib) = ml(j,i1+i-1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do

        ! do the send-up exchange
        call cyclic_exchange_array(r4vector1,r4vector2, &
                                   ssize,ma%bottom,ma%top,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the top boundary of this block
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                ml(j,i2+i,k,n) = r4vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do

        !******************************************************
        ! (4) Send top boundary to bottom side of top neighbor
        !******************************************************
        ! loop over the given longitudes and number of exchange blocks
        ! and unravel the top boundary into a vector
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                r4vector1(ib) = ml(j,i2-i+1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do

        ! do the send-down exchange
        call cyclic_exchange_array(r4vector1,r4vector2, &
                                   ssize,ma%top,ma%bottom,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the bottom boundary of this block
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                ml(j,i1-i,k,n) = r4vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do
      end if

    else

      ! RIGHT AND LEFT EXCHANGE

      if ( ma%right /= mpi_proc_null) then
        ssize = nex*isize*ksize*nsize
        ib = ipos
        do n = n1, n2
          do k = k1, k2
            do i = i1, i2
              do j = 1, nex
                r4vector1(ib) = ml(j2-j+1,i,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%right,tag_rl,tag_lr, &
                            req(irc),req(irc+4))
        ipos = ipos + ssize
      end if
      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*isize*ksize*nsize
        ib = ipos
        do n = n1, n2
          do k = k1, k2
            do i = i1, i2
              do j = 1, nex
                r4vector1(ib) = ml(j1+j-1,i,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%left,tag_lr,tag_rl, &
                            req(irc),req(irc+4))
        ipos = ipos + ssize
      end if

    end if

    if ( .not. ma%crmflag ) then

      ! BOTTOM AND TOP EXCHANGE

      if ( ma%top /= mpi_proc_null) then
        ssize = nex*jsize*ksize*nsize
        ib = ipos
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                r4vector1(ib) = ml(j,i2-i+1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%top,tag_tb,tag_bt, &
                            req(irc),req(irc+4))
        ipos = ipos + ssize
      end if
      if ( ma%bottom /= mpi_proc_null) then
        ssize = nex*jsize*ksize*nsize
        ib = ipos
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                r4vector1(ib) = ml(j,i1+i-1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%bottom,tag_bt,tag_tb, &
                            req(irc),req(irc+4))
        ipos = ipos + ssize
      end if

    end if

    ! Finalize non cyclic comms

    if ( irc /= 0 ) then
      ipos = 1
      call mpi_waitall(8,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      if ( .not. ma%bandflag ) then
        if ( ma%right /= mpi_proc_null) then
          ib = ipos
          ssize = nex*isize*ksize*nsize
          do n = n1, n2
            do k = k1, k2
              do i = i1, i2
                do j = 1, nex
                  ml(j2+j,i,k,n) = r4vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%left /= mpi_proc_null ) then
          ib = ipos
          ssize = nex*isize*ksize*nsize
          do n = n1, n2
            do k = k1, k2
              do i = i1, i2
                do j = 1, nex
                  ml(j1-j,i,k,n) = r4vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
      if ( .not. ma%crmflag ) then
        if ( ma%top /= mpi_proc_null) then
          ib = ipos
          ssize = nex*jsize*ksize*nsize
          do n = n1, n2
            do k = k1, k2
              do i = 1, nex
                do j = j1, j2
                  ml(j,i2+i,k,n) = r4vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%bottom /= mpi_proc_null) then
          ib = ipos
          ssize = nex*jsize*ksize*nsize
          do n = n1, n2
            do k = k1, k2
              do i = 1, nex
                do j = j1, j2
                  ml(j,i1-i,k,n) = r4vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
    end if
  end subroutine real4_4d_exchange_left_right_bottom_top

  subroutine real8_2d_exchange_left_right(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2
    integer(ik4) :: isize, ssize, j, i, ib, irc, ipos
    integer(ik4), dimension(4) :: req

    isize = i2-i1+1

    irc = 0
    ipos = 1
    req = mpi_request_null

    ! Total max possible communication size
    ! set the size of the dummy exchange vector
    ! nex is the width of the exchange stencil
    ! isize is the height of a block in the N-S direction
    ssize = nex * 2*isize
    if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'real8_2d_exchange_left_right')
    end if
    if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'real8_2d_exchange_left_right')
    end if

    if ( ma%bandflag ) then
      ssize = nex*isize
      !********************************************************
      ! (1) Send right boundary to left side of right neighbor
      !********************************************************
      ! loop over the given latitudes and number of exchange blocks
      ! and unravel the right boundary into a vector
      ib = 1
      do i = i1, i2
        do j = 1, nex
          r8vector1(ib) = ml(j2-j+1,i)
          ib = ib + 1
        end do
      end do
      ! do the send-right exchange
      call cyclic_exchange_array(r8vector1,r8vector2, &
                                 ssize,ma%right,ma%left,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the left boundary of this block
      ib = 1
      do i = i1, i2
        do j = 1, nex
          ml(j1-j,i) = r8vector2(ib)
          ib = ib + 1
        end do
      end do

      !********************************************************
      ! (2) Send left boundary to right side of left neighbor
      !********************************************************
      ! loop over the given latitudes and number of exchange blocks
      ! and unravel the left boundary into a vector
      ib = 1
      do i = i1, i2
        do j = 1, nex
          r8vector1(ib) = ml(j1+j-1,i)
          ib = ib + 1
        end do
      end do
      ! do the send-left exchange
      call cyclic_exchange_array(r8vector1,r8vector2, &
                                 ssize,ma%left,ma%right,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the right boundary of this block
      ib = 1
      do i = i1, i2
        do j = 1, nex
          ml(j2+j,i) = r8vector2(ib)
          ib = ib + 1
        end do
      end do

    else

      ! RIGHT AND LEFT EXCHANGE

      if ( ma%right /= mpi_proc_null) then
        ssize = nex*isize
        ib = ipos
        do i = i1, i2
          do j = 1, nex
            r8vector1(ib) = ml(j2-j+1,i)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%right,tag_rl,tag_lr, &
                            req(irc),req(irc+2))
        ipos = ipos + ssize
      end if
      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*isize
        ib = ipos
        do i = i1, i2
          do j = 1, nex
            r8vector1(ib) = ml(j1+j-1,i)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%left,tag_lr,tag_rl, &
                            req(irc),req(irc+2))
        ipos = ipos + ssize
      end if

    end if

    ! Finalize non cyclic comms

    if ( irc /= 0 ) then
      ipos = 1
      call mpi_waitall(4,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      if ( .not. ma%bandflag ) then
        if ( ma%right /= mpi_proc_null) then
          ib = ipos
          ssize = nex*isize
          do i = i1, i2
            do j = 1, nex
              ml(j2+j,i) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%left /= mpi_proc_null ) then
          ib = ipos
          ssize = nex*isize
          do i = i1, i2
            do j = 1, nex
              ml(j1-j,i) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
    end if
  end subroutine real8_2d_exchange_left_right

  subroutine real8_3d_exchange_left_right(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2
    integer(ik4) :: isize, ksize, ssize, j, i, k, ib, irc, ipos
    integer(ik4), dimension(4) :: req

    isize = i2-i1+1
    ksize = k2-k1+1

    irc = 0
    ipos = 1
    req = mpi_request_null

    ! Total max possible communication size
    ! set the size of the dummy exchange vector
    ! nex is the width of the exchange stencil
    ! isize is the height of a block in the N-S direction
    ! ksize is the height of a block in the T-B direction
    ssize = nex * ksize * 2*isize
    if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'real8_3d_exchange_left_right')
    end if
    if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'real8_3d_exchange_left_right')
    end if

    if ( ma%bandflag ) then
      ! set the size of the dummy exchange vector
      ! nex is the width of the exchange stencil
      ! isize is the height of a block in the N-S direction
      ssize = nex*isize*ksize

      !********************************************************
      ! (1) Send right boundary to left side of right neighbor
      !********************************************************
      ! loop over the given latitudes and number of exchange blocks
      ! and unravel the right boundary into a vector
      do k = k1, k2
        do i = i1, i2
          do j = 1, nex
            ib = j + (i - i1) * nex + (k - k1) * (i2 - i1 + 1) * nex
            r8vector1(ib) = ml(j2-j+1,i,k)
          end do
        end do
      end do
      ! do the send-right exchange
      call cyclic_exchange_array(r8vector1,r8vector2, &
                                 ssize,ma%right,ma%left,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the left boundary of this block
      do k = k1, k2
        do i = i1, i2
          do j = 1, nex
            ib = j + (i - i1) * nex + (k - k1) * (i2 - i1 + 1) * nex
            ml(j1-j,i,k) = r8vector2(ib)
          end do
        end do
      end do

      !********************************************************
      ! (2) Send left boundary to right side of left neighbor
      !********************************************************
      ! loop over the given latitudes and number of exchange blocks
      ! and unravel the left boundary into a vector
      do k = k1, k2
        do i = i1, i2
          do j = 1, nex
            ib = j + (i - i1) * nex + (k - k1) * (i2 - i1 + 1) * nex
            r8vector1(ib) = ml(j1+j-1,i,k)
          end do
        end do
      end do
      ! do the send-left exchange
      call cyclic_exchange_array(r8vector1,r8vector2, &
                                 ssize,ma%left,ma%right,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the right boundary of this block
      do k = k1, k2
        do i = i1, i2
          do j = 1, nex
            ib = j + (i - i1) * nex + (k - k1) * (i2 - i1 + 1) * nex
            ml(j2+j,i,k) = r8vector2(ib)
          end do
        end do
      end do

    else

      ! RIGHT AND LEFT EXCHANGE

      if ( ma%right /= mpi_proc_null) then
        ssize = nex*isize*ksize
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              ib = j + (i - i1) * nex + &
                (k - k1) * (i2 - i1 + 1) * nex + ipos - 1
              r8vector1(ib) = ml(j2-j+1,i,k)
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%right,tag_rl,tag_lr, &
                            req(irc),req(irc+2))
        ipos = ipos + ssize
      end if
      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*isize*ksize
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              ib = j + (i - i1) * nex + &
                (k - k1) * (i2 - i1 + 1) * nex + ipos - 1
              r8vector1(ib) = ml(j1+j-1,i,k)
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%left,tag_lr,tag_rl, &
                            req(irc),req(irc+2))
        ipos = ipos + ssize
      end if

    end if

    ! Finalize non cyclic comms

    if ( irc /= 0 ) then
      ipos = 1
      call mpi_waitall(4,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      if ( .not. ma%bandflag ) then
        if ( ma%right /= mpi_proc_null) then
          ssize = nex*isize*ksize
          do k = k1, k2
            do i = i1, i2
              do j = 1, nex
                ib = j + (i - i1) * nex + &
                  (k - k1) * (i2 - i1 + 1) * nex + ipos - 1
                ml(j2+j,i,k) = r8vector2(ib)
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%left /= mpi_proc_null ) then
          ssize = nex*isize*ksize
          do k = k1, k2
            do i = i1, i2
              do j = 1, nex
                ib = j + (i - i1) * nex + &
                  (k - k1) * (i2 - i1 + 1) * nex + ipos - 1
                ml(j1-j,i,k) = r8vector2(ib)
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
    end if
  end subroutine real8_3d_exchange_left_right

  subroutine real8_4d_exchange_left_right(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2, n1, n2
    integer(ik4) :: isize, ksize, nsize, ssize
    integer(ik4) :: j, i, k, n, ib, irc, ipos
    integer(ik4), dimension(4) :: req

    isize = i2-i1+1
    ksize = k2-k1+1
    nsize = n2-n1+1

    irc = 0
    ipos = 1
    req = mpi_request_null

    ! Total max possible communication size
    ! set the size of the dummy exchange vector
    ! nex is the width of the exchange stencil
    ! isize is the height of a block in the N-S direction
    ! ksize is the height of a block in the T-B direction
    ! nsize is the height of a block in the tracer direction
    ssize = nex * ksize * nsize * 2*isize
    if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'real8_4d_exchange_left_right')
    end if
    if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'real8_4d_exchange_left_right')
    end if

    if ( ma%bandflag ) then
      ! set the size of the dummy exchange vector
      ! nex is the width of the exchange stencil
      ! isize is the height of a block in the N-S direction
      ssize = nex*isize*ksize*nsize

      !********************************************************
      ! (1) Send right boundary to left side of right neighbor
      !********************************************************
      ! loop over the given latitudes and number of exchange blocks
      ! and unravel the right boundary into a vector
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              r8vector1(ib) = ml(j2-j+1,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      ! do the send-right exchange
      call cyclic_exchange_array(r8vector1,r8vector2, &
                                 ssize,ma%right,ma%left,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the left boundary of this block
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              ml(j1-j,i,k,n) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do

      !********************************************************
      ! (2) Send left boundary to right side of left neighbor
      !********************************************************
      ! loop over the given latitudes and number of exchange blocks
      ! and unravel the left boundary into a vector
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              r8vector1(ib) = ml(j1+j-1,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      ! do the send-left exchange
      call cyclic_exchange_array(r8vector1,r8vector2, &
                                 ssize,ma%left,ma%right,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the right boundary of this block
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              ml(j2+j,i,k,n) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do

    else

      ! RIGHT AND LEFT EXCHANGE

      if ( ma%right /= mpi_proc_null) then
        ssize = nex*isize*ksize*nsize
        ib = ipos
        do n = n1, n2
          do k = k1, k2
            do i = i1, i2
              do j = 1, nex
                r8vector1(ib) = ml(j2-j+1,i,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%right,tag_rl,tag_lr, &
                            req(irc),req(irc+2))
        ipos = ipos + ssize
      end if
      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*isize*ksize*nsize
        ib = ipos
        do n = n1, n2
          do k = k1, k2
            do i = i1, i2
              do j = 1, nex
                r8vector1(ib) = ml(j1+j-1,i,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%left,tag_lr,tag_rl, &
                            req(irc),req(irc+2))
        ipos = ipos + ssize
      end if

    end if

    ! Finalize non cyclic comms

    if ( irc /= 0 ) then
      ipos = 1
      call mpi_waitall(4,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      if ( .not. ma%bandflag ) then
        if ( ma%right /= mpi_proc_null) then
          ib = ipos
          ssize = nex*isize*ksize*nsize
          do n = n1, n2
            do k = k1, k2
              do i = i1, i2
                do j = 1, nex
                  ml(j2+j,i,k,n) = r8vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%left /= mpi_proc_null ) then
          ib = ipos
          ssize = nex*isize*ksize*nsize
          do n = n1, n2
            do k = k1, k2
              do i = i1, i2
                do j = 1, nex
                  ml(j1-j,i,k,n) = r8vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
    end if
  end subroutine real8_4d_exchange_left_right

  subroutine real4_2d_exchange_left_right(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2
    integer(ik4) :: isize, ssize, j, i, ib, irc, ipos
    integer(ik4), dimension(4) :: req

    isize = i2-i1+1

    irc = 0
    ipos = 1
    req = mpi_request_null

    ! Total max possible communication size
    ! set the size of the dummy exchange vector
    ! nex is the width of the exchange stencil
    ! isize is the height of a block in the N-S direction
    ssize = nex * 2*isize
    if ( size(r4vector1) < ssize ) then
      call getmem1d(r4vector1,1,ssize,'real4_2d_exchange_left_right')
    end if
    if ( size(r4vector2) < ssize ) then
      call getmem1d(r4vector2,1,ssize,'real4_2d_exchange_left_right')
    end if

    if ( ma%bandflag ) then
      ssize = nex*isize
      !********************************************************
      ! (1) Send right boundary to left side of right neighbor
      !********************************************************
      ! loop over the given latitudes and number of exchange blocks
      ! and unravel the right boundary into a vector
      ib = 1
      do i = i1, i2
        do j = 1, nex
          r4vector1(ib) = ml(j2-j+1,i)
          ib = ib + 1
        end do
      end do
      ! do the send-right exchange
      call cyclic_exchange_array(r4vector1,r4vector2, &
                                 ssize,ma%right,ma%left,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the left boundary of this block
      ib = 1
      do i = i1, i2
        do j = 1, nex
          ml(j1-j,i) = r4vector2(ib)
          ib = ib + 1
        end do
      end do

      !********************************************************
      ! (2) Send left boundary to right side of left neighbor
      !********************************************************
      ! loop over the given latitudes and number of exchange blocks
      ! and unravel the left boundary into a vector
      ib = 1
      do i = i1, i2
        do j = 1, nex
          r4vector1(ib) = ml(j1+j-1,i)
          ib = ib + 1
        end do
      end do
      ! do the send-left exchange
      call cyclic_exchange_array(r4vector1,r4vector2, &
                                 ssize,ma%left,ma%right,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the right boundary of this block
      ib = 1
      do i = i1, i2
        do j = 1, nex
          ml(j2+j,i) = r4vector2(ib)
          ib = ib + 1
        end do
      end do

    else

      ! RIGHT AND LEFT EXCHANGE

      if ( ma%right /= mpi_proc_null) then
        ssize = nex*isize
        ib = ipos
        do i = i1, i2
          do j = 1, nex
            r4vector1(ib) = ml(j2-j+1,i)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%right,tag_rl,tag_lr, &
                            req(irc),req(irc+2))
        ipos = ipos + ssize
      end if
      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*isize
        ib = ipos
        do i = i1, i2
          do j = 1, nex
            r4vector1(ib) = ml(j1+j-1,i)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%left,tag_lr,tag_rl, &
                            req(irc),req(irc+2))
        ipos = ipos + ssize
      end if

    end if

    ! Finalize non cyclic comms

    if ( irc /= 0 ) then
      ipos = 1
      call mpi_waitall(4,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      if ( .not. ma%bandflag ) then
        if ( ma%right /= mpi_proc_null) then
          ib = ipos
          ssize = nex*isize
          do i = i1, i2
            do j = 1, nex
              ml(j2+j,i) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%left /= mpi_proc_null ) then
          ib = ipos
          ssize = nex*isize
          do i = i1, i2
            do j = 1, nex
              ml(j1-j,i) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
    end if
  end subroutine real4_2d_exchange_left_right

  subroutine real4_3d_exchange_left_right(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2
    integer(ik4) :: isize, ksize, ssize, j, i, k, ib, irc, ipos
    integer(ik4), dimension(4) :: req

    isize = i2-i1+1
    ksize = k2-k1+1

    irc = 0
    ipos = 1
    req = mpi_request_null

    ! Total max possible communication size
    ! set the size of the dummy exchange vector
    ! nex is the width of the exchange stencil
    ! isize is the height of a block in the N-S direction
    ! ksize is the height of a block in the T-B direction
    ssize = nex * ksize * 2*isize
    if ( size(r4vector1) < ssize ) then
      call getmem1d(r4vector1,1,ssize,'real4_3d_exchange_left_right')
    end if
    if ( size(r4vector2) < ssize ) then
      call getmem1d(r4vector2,1,ssize,'real4_3d_exchange_left_right')
    end if

    if ( ma%bandflag ) then
      ! set the size of the dummy exchange vector
      ! nex is the width of the exchange stencil
      ! isize is the height of a block in the N-S direction
      ssize = nex*isize*ksize

      !********************************************************
      ! (1) Send right boundary to left side of right neighbor
      !********************************************************
      ! loop over the given latitudes and number of exchange blocks
      ! and unravel the right boundary into a vector
      ib = 1
      do k = k1, k2
        do i = i1, i2
          do j = 1, nex
            r4vector1(ib) = ml(j2-j+1,i,k)
            ib = ib + 1
          end do
        end do
      end do
      ! do the send-right exchange
      call cyclic_exchange_array(r4vector1,r4vector2, &
                                 ssize,ma%right,ma%left,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the left boundary of this block
      ib = 1
      do k = k1, k2
        do i = i1, i2
          do j = 1, nex
            ml(j1-j,i,k) = r4vector2(ib)
            ib = ib + 1
          end do
        end do
      end do

      !********************************************************
      ! (2) Send left boundary to right side of left neighbor
      !********************************************************
      ! loop over the given latitudes and number of exchange blocks
      ! and unravel the left boundary into a vector
      ib = 1
      do k = k1, k2
        do i = i1, i2
          do j = 1, nex
            r4vector1(ib) = ml(j1+j-1,i,k)
            ib = ib + 1
          end do
        end do
      end do
      ! do the send-left exchange
      call cyclic_exchange_array(r4vector1,r4vector2, &
                                 ssize,ma%left,ma%right,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the right boundary of this block
      ib = 1
      do k = k1, k2
        do i = i1, i2
          do j = 1, nex
            ml(j2+j,i,k) = r4vector2(ib)
            ib = ib + 1
          end do
        end do
      end do

    else

      ! RIGHT AND LEFT EXCHANGE

      if ( ma%right /= mpi_proc_null) then
        ssize = nex*isize*ksize
        ib = ipos
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              r4vector1(ib) = ml(j2-j+1,i,k)
              ib = ib + 1
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%right,tag_rl,tag_lr, &
                            req(irc),req(irc+2))
        ipos = ipos + ssize
      end if
      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*isize*ksize
        ib = ipos
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              r4vector1(ib) = ml(j1+j-1,i,k)
              ib = ib + 1
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%left,tag_lr,tag_rl, &
                            req(irc),req(irc+2))
        ipos = ipos + ssize
      end if

    end if

    ! Finalize non cyclic comms

    if ( irc /= 0 ) then
      ipos = 1
      call mpi_waitall(4,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      if ( .not. ma%bandflag ) then
        if ( ma%right /= mpi_proc_null) then
          ib = ipos
          ssize = nex*isize*ksize
          do k = k1, k2
            do i = i1, i2
              do j = 1, nex
                ml(j2+j,i,k) = r4vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%left /= mpi_proc_null ) then
          ib = ipos
          ssize = nex*isize*ksize
          do k = k1, k2
            do i = i1, i2
              do j = 1, nex
                ml(j1-j,i,k) = r4vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
    end if
  end subroutine real4_3d_exchange_left_right

  subroutine real4_4d_exchange_left_right(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2, n1, n2
    integer(ik4) :: isize, ksize, nsize, ssize
    integer(ik4) :: j, i, k, n, ib, irc, ipos
    integer(ik4), dimension(4) :: req

    isize = i2-i1+1
    ksize = k2-k1+1
    nsize = n2-n1+1

    irc = 0
    ipos = 1
    req = mpi_request_null

    ! Total max possible communication size
    ! set the size of the dummy exchange vector
    ! nex is the width of the exchange stencil
    ! isize is the height of a block in the N-S direction
    ! ksize is the height of a block in the T-B direction
    ! nsize is the height of a block in the tracer direction
    ssize = nex * ksize * nsize * 2*isize
    if ( size(r4vector1) < ssize ) then
      call getmem1d(r4vector1,1,ssize,'real4_4d_exchange_left_right')
    end if
    if ( size(r4vector2) < ssize ) then
      call getmem1d(r4vector2,1,ssize,'real4_4d_exchange_left_right')
    end if

    if ( ma%bandflag ) then
      ! set the size of the dummy exchange vector
      ! nex is the width of the exchange stencil
      ! isize is the height of a block in the N-S direction
      ssize = nex*isize*ksize*nsize

      !********************************************************
      ! (1) Send right boundary to left side of right neighbor
      !********************************************************
      ! loop over the given latitudes and number of exchange blocks
      ! and unravel the right boundary into a vector
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              r4vector1(ib) = ml(j2-j+1,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      ! do the send-right exchange
      call cyclic_exchange_array(r4vector1,r4vector2, &
                                 ssize,ma%right,ma%left,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the left boundary of this block
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              ml(j1-j,i,k,n) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do

      !********************************************************
      ! (2) Send left boundary to right side of left neighbor
      !********************************************************
      ! loop over the given latitudes and number of exchange blocks
      ! and unravel the left boundary into a vector
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              r4vector1(ib) = ml(j1+j-1,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      ! do the send-left exchange
      call cyclic_exchange_array(r4vector1,r4vector2, &
                                 ssize,ma%left,ma%right,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the right boundary of this block
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              ml(j2+j,i,k,n) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do

    else

      ! RIGHT AND LEFT EXCHANGE

      if ( ma%right /= mpi_proc_null) then
        ssize = nex*isize*ksize*nsize
        ib = ipos
        do n = n1, n2
          do k = k1, k2
            do i = i1, i2
              do j = 1, nex
                r4vector1(ib) = ml(j2-j+1,i,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%right,tag_rl,tag_lr, &
                            req(irc),req(irc+2))
        ipos = ipos + ssize
      end if
      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*isize*ksize*nsize
        ib = ipos
        do n = n1, n2
          do k = k1, k2
            do i = i1, i2
              do j = 1, nex
                r4vector1(ib) = ml(j1+j-1,i,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%left,tag_lr,tag_rl, &
                            req(irc),req(irc+2))
        ipos = ipos + ssize
      end if

    end if

    ! Finalize non cyclic comms

    if ( irc /= 0 ) then
      ipos = 1
      call mpi_waitall(4,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      if ( .not. ma%bandflag ) then
        if ( ma%right /= mpi_proc_null) then
          ib = ipos
          ssize = nex*isize*ksize*nsize
          do n = n1, n2
            do k = k1, k2
              do i = i1, i2
                do j = 1, nex
                  ml(j2+j,i,k,n) = r4vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%left /= mpi_proc_null ) then
          ib = ipos
          ssize = nex*isize*ksize*nsize
          do n = n1, n2
            do k = k1, k2
              do i = i1, i2
                do j = 1, nex
                  ml(j1-j,i,k,n) = r4vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
    end if
  end subroutine real4_4d_exchange_left_right

  subroutine real8_2d_exchange_bottom_top(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2
    integer(ik4) :: jsize, ssize, j, i, ib, irc, ipos
    integer(ik4), dimension(4) :: req

    jsize = j2-j1+1

    irc = 0
    ipos = 1
    req = mpi_request_null

    ! Total max possible communication size
    ! set the size of the dummy exchange vector
    ! nex is the width of the exchange stencil
    ! jsize is the height of a block in the W-E direction
    ssize = nex * 2*jsize
    if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'real8_2d_exchange_bottom_top')
    end if
    if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'real8_2d_exchange_bottom_top')
    end if

    if ( ma%crmflag ) then
      !*********************************************************
      ! (1) Send bottom boundary to top side of bottom neighbor
      !*********************************************************
      ssize = nex*jsize
      ! loop over the given longitudes and number of exchange blocks
      ! and unravel the bottom boundary into a vector
      ib = 1
      do i = 1, nex
        do j = j1, j2
          r8vector1(ib) = ml(j,i1+i-1)
          ib = ib + 1
        end do
      end do

      ! do the send-up exchange
      call cyclic_exchange_array(r8vector1,r8vector2, &
                                 ssize,ma%bottom,ma%top,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the top boundary of this block
      ib = 1
      do i = 1, nex
        do j = j1, j2
          ml(j,i2+i) = r8vector2(ib)
          ib = ib + 1
        end do
      end do

      !******************************************************
      ! (2) Send top boundary to bottom side of top neighbor
      !******************************************************
      ! loop over the given longitudes and number of exchange blocks
      ! and unravel the top boundary into a vector
      ib = 1
      do i = 1, nex
        do j = j1, j2
          r8vector1(ib) = ml(j,i2-i+1)
          ib = ib + 1
        end do
      end do

      ! do the send-down exchange
      call cyclic_exchange_array(r8vector1,r8vector2, &
                                 ssize,ma%top,ma%bottom,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the bottom boundary of this block
      ib = 1
      do i = 1, nex
        do j = j1, j2
          ml(j,i1-i) = r8vector2(ib)
          ib = ib + 1
        end do
      end do

    else

      ! BOTTOM AND TOP EXCHANGE

      if ( ma%top /= mpi_proc_null) then
        ssize = nex*jsize
        ib = ipos
        do i = 1, nex
          do j = j1, j2
            r8vector1(ib) = ml(j,i2-i+1)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%top,tag_tb,tag_bt, &
                            req(irc),req(irc+2))
        ipos = ipos + ssize
      end if
      if ( ma%bottom /= mpi_proc_null) then
        ssize = nex*jsize
        ib = ipos
        do i = 1, nex
          do j = j1, j2
            r8vector1(ib) = ml(j,i1+i-1)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%bottom,tag_bt,tag_tb, &
                            req(irc),req(irc+2))
        ipos = ipos + ssize
      end if

    end if

    ! Finalize non cyclic comms

    if ( irc /= 0 ) then
      ipos = 1
      call mpi_waitall(4,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      if ( .not. ma%crmflag ) then
        if ( ma%top /= mpi_proc_null) then
          ib = ipos
          ssize = nex*jsize
          do i = 1, nex
            do j = j1, j2
              ml(j,i2+i) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%bottom /= mpi_proc_null) then
          ib = ipos
          ssize = nex*jsize
          do i = 1, nex
            do j = j1, j2
              ml(j,i1-i) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
    end if
  end subroutine real8_2d_exchange_bottom_top

  subroutine real8_3d_exchange_bottom_top(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2
    integer(ik4) :: jsize, ksize, ssize, j, i, k, ib, irc, ipos
    integer(ik4), dimension(4) :: req

    jsize = j2-j1+1
    ksize = k2-k1+1

    irc = 0
    ipos = 1
    req = mpi_request_null

    ! Total max possible communication size
    ! set the size of the dummy exchange vector
    ! nex is the width of the exchange stencil
    ! jsize is the height of a block in the W-E direction
    ! ksize is the height of a block in the T-B direction
    ssize = nex * ksize * 2*jsize
    if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'real8_3d_exchange_bottom_top')
    end if
    if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'real8_3d_exchange_bottom_top')
    end if

    if ( ma%crmflag ) then
      !*********************************************************
      ! (1) Send bottom boundary to top side of bottom neighbor
      !*********************************************************
      ! set the size of the dummy exchange vector
      ! nex is the width of the exchange stencil
      ! jsize is the width of a block in the E-W direction
      ssize = nex*jsize*ksize
      ! loop over the given longitudes and number of exchange blocks
      ! and unravel the bottom boundary into a vector
      do k = k1, k2
        do i = 1, nex
          do j = j1, j2
            ib = (j - j1 + 1) + (i - 1) * (j2 - j1 + 1) + &
              (k - k1) * (j2 - j1 + 1) * nex
            r8vector1(ib) = ml(j,i1+i-1,k)
          end do
        end do
      end do

      ! do the send-up exchange
      call cyclic_exchange_array(r8vector1,r8vector2, &
                                 ssize,ma%bottom,ma%top,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the top boundary of this block
      do k = k1, k2
        do i = 1, nex
          do j = j1, j2
            ib = (j - j1 + 1) + (i - 1) * (j2 - j1 + 1) + &
              (k - k1) * (j2 - j1 + 1) * nex
            ml(j,i2+i,k) = r8vector2(ib)
          end do
        end do
      end do

      !******************************************************
      ! (2) Send top boundary to bottom side of top neighbor
      !******************************************************
      ! loop over the given longitudes and number of exchange blocks
      ! and unravel the top boundary into a vector
      do k = k1, k2
        do i = 1, nex
          do j = j1, j2
            ib = (j - j1 + 1) + (i - 1) * (j2 - j1 + 1) + &
              (k - k1) * (j2 - j1 + 1) * nex
            r8vector1(ib) = ml(j,i2-i+1,k)
          end do
        end do
      end do

      ! do the send-down exchange
      call cyclic_exchange_array(r8vector1,r8vector2, &
                                 ssize,ma%top,ma%bottom,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the bottom boundary of this block
      do k = k1, k2
        do i = 1, nex
          do j = j1, j2
            ib = (j - j1 + 1) + (i - 1) * (j2 - j1 + 1) + &
              (k - k1) * (j2 - j1 + 1) * nex
            ml(j,i1-i,k) = r8vector2(ib)
          end do
        end do
      end do

    else

      ! BOTTOM AND TOP EXCHANGE

      if ( ma%top /= mpi_proc_null) then
        ssize = nex*jsize*ksize
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              ib = (j - j1 + 1) + (i - 1) * (j2 - j1 + 1) + &
                (k - k1) * (j2 - j1 + 1) * nex + ipos - 1
              r8vector1(ib) = ml(j,i2-i+1,k)
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%top,tag_tb,tag_bt, &
                            req(irc),req(irc+2))
        ipos = ipos + ssize
      end if
      if ( ma%bottom /= mpi_proc_null) then
        ssize = nex*jsize*ksize
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              ib = (j - j1 + 1) + (i - 1) * (j2 - j1 + 1) + &
                (k - k1) * (j2 - j1 + 1) * nex + ipos - 1
              r8vector1(ib) = ml(j,i1+i-1,k)
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%bottom,tag_bt,tag_tb, &
                            req(irc),req(irc+2))
        ipos = ipos + ssize
      end if

    end if

    ! Finalize non cyclic comms

    if ( irc /= 0 ) then
      ipos = 1
      call mpi_waitall(4,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      if ( .not. ma%crmflag ) then
        if ( ma%top /= mpi_proc_null) then
          ssize = nex*jsize*ksize
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                ib = (j - j1 + 1) + (i - 1) * (j2 - j1 + 1) + &
                  (k - k1) * (j2 - j1 + 1) * nex + ipos - 1
                ml(j,i2+i,k) = r8vector2(ib)
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%bottom /= mpi_proc_null) then
          ssize = nex*jsize*ksize
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                ib = (j - j1 + 1) + (i - 1) * (j2 - j1 + 1) + &
                  (k - k1) * (j2 - j1 + 1) * nex + ipos - 1
                ml(j,i1-i,k) = r8vector2(ib)
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
    end if
  end subroutine real8_3d_exchange_bottom_top

  subroutine real8_4d_exchange_bottom_top(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2, n1, n2
    integer(ik4) :: jsize, ksize, nsize, ssize
    integer(ik4) :: j, i, k, n, ib, irc, ipos
    integer(ik4), dimension(4) :: req

    jsize = j2-j1+1
    ksize = k2-k1+1
    nsize = n2-n1+1

    irc = 0
    ipos = 1
    req = mpi_request_null

    ! Total max possible communication size
    ! set the size of the dummy exchange vector
    ! nex is the width of the exchange stencil
    ! jsize is the height of a block in the W-E direction
    ! ksize is the height of a block in the T-B direction
    ! nsize is the height of a block in the tracer direction
    ssize = nex * ksize * nsize * 2*jsize
    if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'real8_4d_exchange_bottom_top')
    end if
    if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'real8_4d_exchange_bottom_top')
    end if

    if ( ma%crmflag ) then
      !*********************************************************
      ! (1) Send bottom boundary to top side of bottom neighbor
      !*********************************************************
      ! set the size of the dummy exchange vector
      ! nex is the width of the exchange stencil
      ! jsize is the width of a block in the E-W direction
      ssize = nex*jsize*ksize*nsize
      ! loop over the given longitudes and number of exchange blocks
      ! and unravel the bottom boundary into a vector
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              r8vector1(ib) = ml(j,i1+i-1,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do

      ! do the send-up exchange
      call cyclic_exchange_array(r8vector1,r8vector2, &
                                 ssize,ma%bottom,ma%top,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the top boundary of this block
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              ml(j,i2+i,k,n) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do

      !******************************************************
      ! (2) Send top boundary to bottom side of top neighbor
      !******************************************************
      ! loop over the given longitudes and number of exchange blocks
      ! and unravel the top boundary into a vector
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              r8vector1(ib) = ml(j,i2-i+1,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do

      ! do the send-down exchange
      call cyclic_exchange_array(r8vector1,r8vector2, &
                                 ssize,ma%top,ma%bottom,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the bottom boundary of this block
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              ml(j,i1-i,k,n) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do

    else

      ! BOTTOM AND TOP EXCHANGE

      if ( ma%top /= mpi_proc_null) then
        ssize = nex*jsize*ksize*nsize
        ib = ipos
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                r8vector1(ib) = ml(j,i2-i+1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%top,tag_tb,tag_bt, &
                            req(irc),req(irc+2))
        ipos = ipos + ssize
      end if
      if ( ma%bottom /= mpi_proc_null) then
        ssize = nex*jsize*ksize*nsize
        ib = ipos
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                r8vector1(ib) = ml(j,i1+i-1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ssize), &
                            r8vector2(ipos:ipos+ssize),ssize, &
                            ma%bottom,tag_bt,tag_tb, &
                            req(irc),req(irc+2))
        ipos = ipos + ssize
      end if

    end if

    ! Finalize non cyclic comms

    if ( irc /= 0 ) then
      ipos = 1
      call mpi_waitall(4,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      if ( .not. ma%crmflag ) then
        if ( ma%top /= mpi_proc_null) then
          ib = ipos
          ssize = nex*jsize*ksize*nsize
          do n = n1, n2
            do k = k1, k2
              do i = 1, nex
                do j = j1, j2
                  ml(j,i2+i,k,n) = r8vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%bottom /= mpi_proc_null) then
          ib = ipos
          ssize = nex*jsize*ksize*nsize
          do n = n1, n2
            do k = k1, k2
              do i = 1, nex
                do j = j1, j2
                  ml(j,i1-i,k,n) = r8vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
    end if
  end subroutine real8_4d_exchange_bottom_top

  subroutine real4_2d_exchange_bottom_top(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2
    integer(ik4) :: jsize, ssize, j, i, ib, irc, ipos
    integer(ik4), dimension(4) :: req

    jsize = j2-j1+1

    irc = 0
    ipos = 1
    req = mpi_request_null

    ! Total max possible communication size
    ! set the size of the dummy exchange vector
    ! nex is the width of the exchange stencil
    ! jsize is the height of a block in the W-E direction
    ssize = nex * 2*jsize
    if ( size(r4vector1) < ssize ) then
      call getmem1d(r4vector1,1,ssize,'real4_2d_exchange_bottom_top')
    end if
    if ( size(r4vector2) < ssize ) then
      call getmem1d(r4vector2,1,ssize,'real4_2d_exchange_bottom_top')
    end if

    if ( ma%crmflag ) then
      !*********************************************************
      ! (1) Send bottom boundary to top side of bottom neighbor
      !*********************************************************
      ssize = nex*jsize
      ! loop over the given longitudes and number of exchange blocks
      ! and unravel the bottom boundary into a vector
      ib = 1
      do i = 1, nex
        do j = j1, j2
          r4vector1(ib) = ml(j,i1+i-1)
          ib = ib + 1
        end do
      end do

      ! do the send-up exchange
      call cyclic_exchange_array(r4vector1,r4vector2, &
                                 ssize,ma%bottom,ma%top,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the top boundary of this block
      ib = 1
      do i = 1, nex
        do j = j1, j2
          ml(j,i2+i) = r4vector2(ib)
          ib = ib + 1
        end do
      end do

      !******************************************************
      ! (2) Send top boundary to bottom side of top neighbor
      !******************************************************
      ! loop over the given longitudes and number of exchange blocks
      ! and unravel the top boundary into a vector
      ib = 1
      do i = 1, nex
        do j = j1, j2
          r4vector1(ib) = ml(j,i2-i+1)
          ib = ib + 1
        end do
      end do

      ! do the send-down exchange
      call cyclic_exchange_array(r4vector1,r4vector2, &
                                 ssize,ma%top,ma%bottom,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the bottom boundary of this block
      ib = 1
      do i = 1, nex
        do j = j1, j2
          ml(j,i1-i) = r4vector2(ib)
          ib = ib + 1
        end do
      end do

    else

      ! BOTTOM AND TOP EXCHANGE

      if ( ma%top /= mpi_proc_null) then
        ssize = nex*jsize
        ib = ipos
        do i = 1, nex
          do j = j1, j2
            r4vector1(ib) = ml(j,i2-i+1)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%top,tag_tb,tag_bt, &
                            req(irc),req(irc+2))
        ipos = ipos + ssize
      end if
      if ( ma%bottom /= mpi_proc_null) then
        ssize = nex*jsize
        ib = ipos
        do i = 1, nex
          do j = j1, j2
            r4vector1(ib) = ml(j,i1+i-1)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%bottom,tag_bt,tag_tb, &
                            req(irc),req(irc+2))
        ipos = ipos + ssize
      end if

    end if

    ! Finalize non cyclic comms

    if ( irc /= 0 ) then
      ipos = 1
      call mpi_waitall(4,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      if ( .not. ma%crmflag ) then
        if ( ma%top /= mpi_proc_null) then
          ib = ipos
          ssize = nex*jsize
          do i = 1, nex
            do j = j1, j2
              ml(j,i2+i) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%bottom /= mpi_proc_null) then
          ib = ipos
          ssize = nex*jsize
          do i = 1, nex
            do j = j1, j2
              ml(j,i1-i) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
    end if
  end subroutine real4_2d_exchange_bottom_top

  subroutine real4_3d_exchange_bottom_top(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2
    integer(ik4) :: jsize, ksize, ssize, j, i, k, ib, irc, ipos
    integer(ik4), dimension(4) :: req

    jsize = j2-j1+1
    ksize = k2-k1+1

    irc = 0
    ipos = 1
    req = mpi_request_null

    ! Total max possible communication size
    ! set the size of the dummy exchange vector
    ! nex is the width of the exchange stencil
    ! jsize is the height of a block in the W-E direction
    ! ksize is the height of a block in the T-B direction
    ssize = nex * ksize * 2*jsize
    if ( size(r4vector1) < ssize ) then
      call getmem1d(r4vector1,1,ssize,'real4_3d_exchange_bottom_top')
    end if
    if ( size(r4vector2) < ssize ) then
      call getmem1d(r4vector2,1,ssize,'real4_3d_exchange_bottom_top')
    end if

    if ( ma%crmflag ) then
      !*********************************************************
      ! (1) Send bottom boundary to top side of bottom neighbor
      !*********************************************************
      ! set the size of the dummy exchange vector
      ! nex is the width of the exchange stencil
      ! jsize is the width of a block in the E-W direction
      ssize = nex*jsize*ksize
      ! loop over the given longitudes and number of exchange blocks
      ! and unravel the bottom boundary into a vector
      ib = 1
      do k = k1, k2
        do i = 1, nex
          do j = j1, j2
            r4vector1(ib) = ml(j,i1+i-1,k)
            ib = ib + 1
          end do
        end do
      end do

      ! do the send-up exchange
      call cyclic_exchange_array(r4vector1,r4vector2, &
                                 ssize,ma%bottom,ma%top,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the top boundary of this block
      ib = 1
      do k = k1, k2
        do i = 1, nex
          do j = j1, j2
            ml(j,i2+i,k) = r4vector2(ib)
            ib = ib + 1
          end do
        end do
      end do

      !******************************************************
      ! (2) Send top boundary to bottom side of top neighbor
      !******************************************************
      ! loop over the given longitudes and number of exchange blocks
      ! and unravel the top boundary into a vector
      ib = 1
      do k = k1, k2
        do i = 1, nex
          do j = j1, j2
            r4vector1(ib) = ml(j,i2-i+1,k)
            ib = ib + 1
          end do
        end do
      end do

      ! do the send-down exchange
      call cyclic_exchange_array(r4vector1,r4vector2, &
                                 ssize,ma%top,ma%bottom,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the bottom boundary of this block
      ib = 1
      do k = k1, k2
        do i = 1, nex
          do j = j1, j2
            ml(j,i1-i,k) = r4vector2(ib)
            ib = ib + 1
          end do
        end do
      end do

    else

      ! BOTTOM AND TOP EXCHANGE

      if ( ma%top /= mpi_proc_null) then
        ssize = nex*jsize*ksize
        ib = ipos
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              r4vector1(ib) = ml(j,i2-i+1,k)
              ib = ib + 1
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%top,tag_tb,tag_bt, &
                            req(irc),req(irc+2))
        ipos = ipos + ssize
      end if
      if ( ma%bottom /= mpi_proc_null) then
        ssize = nex*jsize*ksize
        ib = ipos
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              r4vector1(ib) = ml(j,i1+i-1,k)
              ib = ib + 1
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%bottom,tag_bt,tag_tb, &
                            req(irc),req(irc+2))
        ipos = ipos + ssize
      end if

    end if

    ! Finalize non cyclic comms

    if ( irc /= 0 ) then
      ipos = 1
      call mpi_waitall(4,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      if ( .not. ma%crmflag ) then
        if ( ma%top /= mpi_proc_null) then
          ib = ipos
          ssize = nex*jsize*ksize
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                ml(j,i2+i,k) = r4vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%bottom /= mpi_proc_null) then
          ib = ipos
          ssize = nex*jsize*ksize
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                ml(j,i1-i,k) = r4vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
    end if
  end subroutine real4_3d_exchange_bottom_top

  subroutine real4_4d_exchange_bottom_top(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2, n1, n2
    integer(ik4) :: jsize, ksize, nsize, ssize
    integer(ik4) :: j, i, k, n, ib, irc, ipos
    integer(ik4), dimension(4) :: req

    jsize = j2-j1+1
    ksize = k2-k1+1
    nsize = n2-n1+1

    irc = 0
    ipos = 1
    req = mpi_request_null

    ! Total max possible communication size
    ! set the size of the dummy exchange vector
    ! nex is the width of the exchange stencil
    ! jsize is the height of a block in the W-E direction
    ! ksize is the height of a block in the T-B direction
    ! nsize is the height of a block in the tracer direction
    ssize = nex * ksize * nsize * 2*jsize
    if ( size(r4vector1) < ssize ) then
      call getmem1d(r4vector1,1,ssize,'real4_4d_exchange_bottom_top')
    end if
    if ( size(r4vector2) < ssize ) then
      call getmem1d(r4vector2,1,ssize,'real4_4d_exchange_bottom_top')
    end if

    if ( ma%crmflag ) then
      !*********************************************************
      ! (1) Send bottom boundary to top side of bottom neighbor
      !*********************************************************
      ! set the size of the dummy exchange vector
      ! nex is the width of the exchange stencil
      ! jsize is the width of a block in the E-W direction
      ssize = nex*jsize*ksize*nsize
      ! loop over the given longitudes and number of exchange blocks
      ! and unravel the bottom boundary into a vector
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              r4vector1(ib) = ml(j,i1+i-1,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do

      ! do the send-up exchange
      call cyclic_exchange_array(r4vector1,r4vector2, &
                                 ssize,ma%bottom,ma%top,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the top boundary of this block
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              ml(j,i2+i,k,n) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do

      !******************************************************
      ! (2) Send top boundary to bottom side of top neighbor
      !******************************************************
      ! loop over the given longitudes and number of exchange blocks
      ! and unravel the top boundary into a vector
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              r4vector1(ib) = ml(j,i2-i+1,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do

      ! do the send-down exchange
      call cyclic_exchange_array(r4vector1,r4vector2, &
                                 ssize,ma%top,ma%bottom,__LINE__)
      ! loop over the receiving dummy vector and ravel it into
      ! the bottom boundary of this block
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              ml(j,i1-i,k,n) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do

    else

      ! BOTTOM AND TOP EXCHANGE

      if ( ma%top /= mpi_proc_null) then
        ssize = nex*jsize*ksize*nsize
        ib = ipos
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                r4vector1(ib) = ml(j,i2-i+1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%top,tag_tb,tag_bt, &
                            req(irc),req(irc+2))
        ipos = ipos + ssize
      end if
      if ( ma%bottom /= mpi_proc_null) then
        ssize = nex*jsize*ksize*nsize
        ib = ipos
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                r4vector1(ib) = ml(j,i1+i-1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ssize), &
                            r4vector2(ipos:ipos+ssize),ssize, &
                            ma%bottom,tag_bt,tag_tb, &
                            req(irc),req(irc+2))
        ipos = ipos + ssize
      end if

    end if

    ! Finalize non cyclic comms

    if ( irc /= 0 ) then
      ipos = 1
      call mpi_waitall(4,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      if ( .not. ma%crmflag ) then
        if ( ma%top /= mpi_proc_null) then
          ib = ipos
          ssize = nex*jsize*ksize*nsize
          do n = n1, n2
            do k = k1, k2
              do i = 1, nex
                do j = j1, j2
                  ml(j,i2+i,k,n) = r4vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%bottom /= mpi_proc_null) then
          ib = ipos
          ssize = nex*jsize*ksize*nsize
          do n = n1, n2
            do k = k1, k2
              do i = 1, nex
                do j = j1, j2
                  ml(j,i1-i,k,n) = r4vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
    end if
  end subroutine real4_4d_exchange_bottom_top

  subroutine real8_2d_exchange_left_bottom(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2
    integer(ik4) :: isize, jsize, ssize, j, i, ib
    integer(ik4) :: spos, rpos, ipos, irc
    integer(ik4), dimension(6) :: req
    isize = i2-i1+1
    jsize = j2-j1+1
    irc = 0
    spos = 1
    rpos = 1
    req = mpi_request_null

    ssize = nex*(isize+jsize+nex)
    if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'real8_2d_exchange_left_bottom')
    end if
    if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'real8_2d_exchange_left_bottom')
    end if

    if ( ma%bandflag ) then
      ssize = nex*isize
      ib = 1
      do i = i1, i2
        do j = 1, nex
          r8vector1(ib) = ml(j2-j+1,i)
          ib = ib + 1
        end do
      end do
      call cyclic_exchange_array(r8vector1,r8vector2, &
                                 ssize,ma%right,ma%left,__LINE__)
      ib = 1
      do i = i1, i2
        do j = 1, nex
          ml(j1-j,i) = r8vector2(ib)
          ib = ib + 1
        end do
      end do

      if ( ma%crmflag ) then
        ! set the size of the dummy exchange vector
        ! nex is the width of the exchange stencil
        ! jsize is the width of a block in the E-W direction
        ssize = nex*jsize

        ! (1) Send top boundary to bottom side of top neighbor
        ! loop over the given longitudes and number of exchange blocks
        ! and unravel the top boundary into a vector
        ib = 1

        do i = 1, nex
          do j = j1, j2
            r8vector1(ib) = ml(j,i2-i+1)
            ib = ib + 1
          end do
        end do

        ! do the send-down exchange
        call cyclic_exchange_array(r8vector1,r8vector2, &
                                   ssize,ma%top,ma%bottom,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the bottom boundary of this block
        ib = 1

        do i = 1, nex
          do j = j1, j2
            ml(j,i1-i) = r8vector2(ib)
            ib = ib + 1
          end do
        end do

        ! (2) Send top-right boundary to bottom-left side of top-right
        ! neighbor: loop over the exchange block and unravel it into
        ! a vector
        ssize = nex*nex
        ib = 1
        do i = 1, nex
          do j = 1, nex
            r8vector1(ib) = ml(j2-j+1,i2-i+1)
            ib = ib + 1
          end do
        end do

        ! do the send up-right exchange
        call cyclic_exchange_array(r8vector1,r8vector2, &
                                   ssize,ma%topright,ma%bottomleft,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the bottom-left boundary of this block
        ib = 1

        do i = 1, nex
          do j = 1, nex
            ml(j1-j,i1-i) = r8vector2(ib)
            ib = ib + 1
          end do
        end do

      end if

    else

      if ( ma%right /= mpi_proc_null) then
        ssize = nex*isize
        ib = spos
        do i = i1, i2
          do j = 1, nex
            r8vector1(ib) = ml(j2-j+1,i)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call send_array(r8vector1(spos:spos+ssize), &
                        ssize,ma%right,tag_lr,req(irc))
        spos = spos + ssize
      end if
      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*isize
        irc = irc + 1
        call recv_array(r8vector2(rpos:rpos+ssize), &
                        ssize,ma%left,tag_lr,req(irc))
        rpos = rpos + ssize
      end if

    end if

    if ( .not. ma%crmflag ) then
      if ( ma%top /= mpi_proc_null) then
        ssize = nex*jsize
        ib = spos
        do i = 1, nex
          do j = j1, j2
            r8vector1(ib) = ml(j,i2-i+1)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call send_array(r8vector1(spos:spos+ssize), &
                        ssize,ma%top,tag_bt,req(irc))
        spos = spos + ssize
      end if
      if ( ma%bottom /= mpi_proc_null) then
        ssize = nex*jsize
        irc = irc + 1
        call recv_array(r8vector2(rpos:rpos+ssize), &
                        ssize,ma%bottom,tag_bt,req(irc))
        rpos = rpos + ssize
      end if
      if ( ma%topright /= mpi_proc_null ) then
        ssize = nex*nex
        ib = spos
        do i = 1, nex
          do j = 1, nex
            r8vector1(ib) = ml(j2-j+1,i2-i+1)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call send_array(r8vector1(spos:spos+ssize), &
                        ssize,ma%topright,tag_bltr,req(irc))
        spos = spos + ssize
      end if
      if ( ma%bottomleft /= mpi_proc_null ) then
        ssize = nex*nex
        irc = irc + 1
        call recv_array(r8vector2(rpos:rpos+ssize), &
                        ssize,ma%bottomleft,tag_bltr,req(irc))
        rpos = rpos + ssize
      end if
    end if
    if ( irc /= 0 ) then
      ipos = 1
      call mpi_waitall(6,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      if ( .not. ma%bandflag ) then
        if ( ma%left /= mpi_proc_null ) then
          ib = ipos
          ssize = nex * isize
          do i = i1, i2
            do j = 1, nex
              ml(j1-j,i) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
      if ( .not. ma%crmflag ) then
        if ( ma%bottom /= mpi_proc_null) then
          ib = ipos
          ssize = nex * jsize
          do i = 1, nex
            do j = j1, j2
              ml(j,i1-i) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%bottomleft /= mpi_proc_null ) then
          ib = ipos
          ssize = nex * nex
          do i = 1, nex
            do j = 1, nex
              ml(j1-j,i1-i) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
    end if
  end subroutine real8_2d_exchange_left_bottom

  subroutine real8_3d_exchange_left_bottom(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2
    integer(ik4) :: isize, jsize, ksize, ssize, j, i, k, ib
    integer(ik4) :: spos, rpos, ipos, irc
    integer(ik4), dimension(6) :: req
    isize = i2-i1+1
    jsize = j2-j1+1
    ksize = k2-k1+1
    irc = 0
    spos = 1
    rpos = 1
    req = mpi_request_null

    ssize = nex*ksize*(isize+jsize+nex)
    if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'real8_3d_exchange_left_bottom')
    end if
    if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'real8_3d_exchange_left_bottom')
    end if

    if ( ma%bandflag ) then
      ssize = nex*isize*ksize
      ib = 1
      do k = k1, k2
        do i = i1, i2
          do j = 1, nex
            r8vector1(ib) = ml(j2-j+1,i,k)
            ib = ib + 1
          end do
        end do
      end do
      call cyclic_exchange_array(r8vector1,r8vector2, &
                                 ssize,ma%right,ma%left,__LINE__)
      ib = 1
      do k = k1, k2
        do i = i1, i2
          do j = 1, nex
            ml(j1-j,i,k) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end do

      if ( ma%crmflag ) then
        ! set the size of the dummy exchange vector
        ! nex is the width of the exchange stencil
        ! jsize is the width of a block in the E-W direction
        ssize = nex*jsize*ksize
        ! (1) Send top boundary to bottom side of top neighbor
        ! loop over the given longitudes and number of exchange blocks
        ! and unravel the top boundary into a vector
        ib = 1

        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              r8vector1(ib) = ml(j,i2-i+1,k)
              ib = ib + 1
            end do
          end do
        end do

        ! do the send-down exchange
        call cyclic_exchange_array(r8vector1,r8vector2, &
                                   ssize,ma%top,ma%bottom,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the bottom boundary of this block
        ib = 1

        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              ml(j,i1-i,k) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do

        ! (2) Send top-right boundary to bottom-left side of
        ! top-right neighbor
        ! loop over the exchange block and unravel it into a vector
        ssize = nex*nex*ksize
        ib = 1

        do k = k1, k2
          do i = 1, nex
            do j = 1, nex
              r8vector1(ib) = ml(j2-j+1,i2-i+1,k)
              ib = ib + 1
            end do
          end do
        end do

        ! do the send up-right exchange
        call cyclic_exchange_array(r8vector1,r8vector2, &
                                   ssize,ma%topright,ma%bottomleft,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the bottom-left boundary of this block
        ib = 1

        do k = k1, k2
          do i = 1, nex
            do j = 1, nex
              ml(j1-j,i1-i,k) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do

      end if

    else

      if ( ma%right /= mpi_proc_null) then
        ssize = nex*ksize*isize
        ib = spos
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              r8vector1(ib) = ml(j2-j+1,i,k)
              ib = ib + 1
            end do
          end do
        end do
        irc = irc + 1
        call send_array(r8vector1(spos:spos+ssize), &
                        ssize,ma%right,tag_lr,req(irc))
        spos = spos + ssize
      end if
      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*ksize*isize
        irc = irc + 1
        call recv_array(r8vector2(rpos:rpos+ssize), &
                        ssize,ma%left,tag_lr,req(irc))
        rpos = rpos + ssize
      end if

    end if

    if ( .not. ma%crmflag ) then
      if ( ma%top /= mpi_proc_null) then
        ssize = nex*ksize*jsize
        ib = spos
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              r8vector1(ib) = ml(j,i2-i+1,k)
              ib = ib + 1
            end do
          end do
        end do
        irc = irc + 1
        call send_array(r8vector1(spos:spos+ssize), &
                        ssize,ma%top,tag_bt,req(irc))
        spos = spos + ssize
      end if
      if ( ma%bottom /= mpi_proc_null) then
        ssize = nex*ksize*jsize
        irc = irc + 1
        call recv_array(r8vector2(rpos:rpos+ssize), &
                        ssize,ma%bottom,tag_bt,req(irc))
        rpos = rpos + ssize
      end if
      if ( ma%topright /= mpi_proc_null ) then
        ssize = nex*ksize*nex
        ib = spos
        do k = k1, k2
          do i = 1, nex
            do j = 1, nex
              r8vector1(ib) = ml(j2-j+1,i2-i+1,k)
              ib = ib + 1
            end do
          end do
        end do
        irc = irc + 1
        call send_array(r8vector1(spos:spos+ssize), &
                        ssize,ma%topright,tag_bltr,req(irc))
        spos = spos + ssize
      end if
      if ( ma%bottomleft /= mpi_proc_null ) then
        ssize = nex*ksize*nex
        irc = irc + 1
        call recv_array(r8vector2(rpos:rpos+ssize), &
                        ssize,ma%bottomleft,tag_bltr,req(irc))
        rpos = rpos + ssize
      end if
    end if
    if ( irc /= 0 ) then
      ipos = 1
      call mpi_waitall(6,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      if ( .not. ma%bandflag ) then
        if ( ma%left /= mpi_proc_null ) then
          ib = ipos
          ssize = nex * isize * ksize
          do k = k1, k2
            do i = i1, i2
              do j = 1, nex
                ml(j1-j,i,k) = r8vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
      if ( .not. ma%crmflag ) then
        if ( ma%bottom /= mpi_proc_null) then
          ib = ipos
          ssize = nex * jsize * ksize
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                ml(j,i1-i,k) = r8vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%bottomleft /= mpi_proc_null ) then
          ib = ipos
          ssize = nex * nex * ksize
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                ml(j1-j,i1-i,k) = r8vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
    end if
  end subroutine real8_3d_exchange_left_bottom

  subroutine real8_4d_exchange_left_bottom(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2, n1, n2
    integer(ik4) :: isize, jsize, ksize, nsize, ssize, j, i, k, n, ib
    integer(ik4) :: spos, rpos, ipos, irc
    integer(ik4), dimension(6) :: req
    isize = i2-i1+1
    jsize = j2-j1+1
    ksize = k2-k1+1
    nsize = n2-n1+1
    irc = 0
    spos = 1
    rpos = 1
    req = mpi_request_null

    ssize = nex*ksize*nsize*(isize+jsize+nex)
    if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'real8_4d_exchange_left_bottom')
    end if
    if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'real8_4d_exchange_left_bottom')
    end if

    if ( ma%bandflag ) then
      ssize = nex*isize*ksize*nsize
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              r8vector1(ib) = ml(j2-j+1,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      call cyclic_exchange_array(r8vector1,r8vector2, &
                                 ssize,ma%right,ma%left,__LINE__)
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              ml(j1-j,i,k,n) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do

      if ( ma%crmflag ) then
        ! set the size of the dummy exchange vector
        ! nex is the width of the exchange stencil
        ! jsize is the width of a block in the E-W direction
        ssize = nex*jsize*ksize*nsize
        ! (1) Send top boundary to bottom side of top neighbor
        ! loop over the given longitudes and number of exchange blocks
        ! and unravel the top boundary into a vector
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                r8vector1(ib) = ml(j,i2-i+1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        ! do the send-down exchange
        call cyclic_exchange_array(r8vector1,r8vector2, &
                                   ssize,ma%top,ma%bottom,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the bottom boundary of this block
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                ml(j,i1-i,k,n) = r8vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do

        ! (2) Send top-right boundary to bottom-left side of top-right neighbor
        ! loop over the exchange block and unravel it into a vector
        ssize = nex*nex*ksize*nsize
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                r8vector1(ib) = ml(j2-j+1,i2-i+1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        ! do the send up-right exchange
        call cyclic_exchange_array(r8vector1,r8vector2, &
                                   ssize,ma%topright,ma%bottomleft,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the bottom-left boundary of this block
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                ml(j1-j,i1-i,k,n) = r8vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do

      end if

    else

      if ( ma%right /= mpi_proc_null) then
        ssize = nex*ksize*nsize*isize
        ib = spos
        do n = n1, n2
          do k = k1, k2
            do i = i1, i2
              do j = 1, nex
                r8vector1(ib) = ml(j2-j+1,i,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call send_array(r8vector1(spos:spos+ssize), &
                        ssize,ma%right,tag_lr,req(irc))
        spos = spos + ssize
      end if
      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*ksize*nsize*isize
        irc = irc + 1
        call recv_array(r8vector2(rpos:rpos+ssize), &
                        ssize,ma%left,tag_lr,req(irc))
        rpos = rpos + ssize
      end if
    end if
    if ( .not. ma%crmflag ) then
      if ( ma%top /= mpi_proc_null) then
        ssize = nex*ksize*nsize*jsize
        ib = spos
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                r8vector1(ib) = ml(j,i2-i+1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call send_array(r8vector1(spos:spos+ssize), &
                        ssize,ma%top,tag_bt,req(irc))
        spos = spos + ssize
      end if
      if ( ma%bottom /= mpi_proc_null) then
        ssize = nex*ksize*nsize*jsize
        irc = irc + 1
        call recv_array(r8vector2(rpos:rpos+ssize), &
                        ssize,ma%bottom,tag_bt,req(irc))
        rpos = rpos + ssize
      end if
      if ( ma%topright /= mpi_proc_null ) then
        ssize = nex*ksize*nsize*nex
        ib = spos
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                r8vector1(ib) = ml(j2-j+1,i2-i+1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call send_array(r8vector1(spos:spos+ssize), &
                        ssize,ma%topright,tag_bltr,req(irc))
        spos = spos + ssize
      end if
      if ( ma%bottomleft /= mpi_proc_null ) then
        ssize = nex*ksize*nsize*nex
        irc = irc + 1
        call recv_array(r8vector2(rpos:rpos+ssize), &
                        ssize,ma%bottomleft,tag_bltr,req(irc))
        rpos = rpos + ssize
      end if
    end if
    if ( irc /= 0 ) then
      ipos = 1
      call mpi_waitall(6,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      if ( .not. ma%bandflag ) then
        if ( ma%left /= mpi_proc_null ) then
          ib = ipos
          ssize = nex * isize * ksize * nsize
          do n = n1, n2
            do k = k1, k2
              do i = i1, i2
                do j = 1, nex
                  ml(j1-j,i,k,n) = r8vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
      if ( .not. ma%crmflag ) then
        if ( ma%bottom /= mpi_proc_null) then
          ib = ipos
          ssize = nex * jsize * ksize * nsize
          do n = n1, n2
            do k = k1, k2
              do i = 1, nex
                do j = j1, j2
                  ml(j,i1-i,k,n) = r8vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%bottomleft /= mpi_proc_null ) then
          ib = ipos
          ssize = nex * nex * ksize * nsize
          do n = n1, n2
            do k = k1, k2
              do i = 1, nex
                do j = 1, nex
                  ml(j1-j,i1-i,k,n) = r8vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
    end if
  end subroutine real8_4d_exchange_left_bottom

  subroutine real4_2d_exchange_left_bottom(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2
    integer(ik4) :: isize, jsize, ssize, j, i, ib
    integer(ik4) :: spos, rpos, ipos, irc
    integer(ik4), dimension(6) :: req
    isize = i2-i1+1
    jsize = j2-j1+1
    irc = 0
    spos = 1
    rpos = 1
    req = mpi_request_null

    ssize = nex*(isize+jsize+nex)
    if ( size(r4vector1) < ssize ) then
      call getmem1d(r4vector1,1,ssize,'real4_2d_exchange_left_bottom')
    end if
    if ( size(r4vector2) < ssize ) then
      call getmem1d(r4vector2,1,ssize,'real4_2d_exchange_left_bottom')
    end if

    if ( ma%bandflag ) then
      ssize = nex*isize
      ib = 1
      do i = i1, i2
        do j = 1, nex
          r4vector1(ib) = ml(j2-j+1,i)
          ib = ib + 1
        end do
      end do
      call cyclic_exchange_array(r4vector1,r4vector2, &
                                 ssize,ma%right,ma%left,__LINE__)
      ib = 1
      do i = i1, i2
        do j = 1, nex
          ml(j1-j,i) = r4vector2(ib)
          ib = ib + 1
        end do
      end do

      if ( ma%crmflag ) then
        ! set the size of the dummy exchange vector
        ! nex is the width of the exchange stencil
        ! jsize is the width of a block in the E-W direction
        ssize = nex*jsize
        ! (1) Send top boundary to bottom side of top neighbor
        ! loop over the given longitudes and number of exchange blocks
        ! and unravel the top boundary into a vector
        ib = 1

        do i = 1, nex
          do j = j1, j2
            r4vector1(ib) = ml(j,i2-i+1)
            ib = ib + 1
          end do
        end do

        ! do the send-down exchange
        call cyclic_exchange_array(r4vector1,r4vector2, &
                                   ssize,ma%top,ma%bottom,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the bottom boundary of this block
        ib = 1

        do i = 1, nex
          do j = j1, j2
            ml(j,i1-i) = r4vector2(ib)
            ib = ib + 1
          end do
        end do

        ! (2) Send top-right boundary to bottom-left side of top-right
        ! neighbor: loop over the exchange block and unravel it into
        ! a vector
        ssize = nex*nex
        ib = 1
        do i = 1, nex
          do j = 1, nex
            r4vector1(ib) = ml(j2-j+1,i2-i+1)
            ib = ib + 1
          end do
        end do

        ! do the send up-right exchange
        call cyclic_exchange_array(r4vector1,r4vector2, &
                                   ssize,ma%topright,ma%bottomleft,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the bottom-left boundary of this block
        ib = 1

        do i = 1, nex
          do j = 1, nex
            ml(j1-j,i1-i) = r4vector2(ib)
            ib = ib + 1
          end do
        end do

      end if

    else

      if ( ma%right /= mpi_proc_null) then
        ssize = nex*isize
        ib = spos
        do i = i1, i2
          do j = 1, nex
            r4vector1(ib) = ml(j2-j+1,i)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call send_array(r4vector1(spos:spos+ssize), &
                        ssize,ma%right,tag_lr,req(irc))
        spos = spos + ssize
      end if
      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*isize
        irc = irc + 1
        call recv_array(r4vector2(rpos:rpos+ssize), &
                        ssize,ma%left,tag_lr,req(irc))
        rpos = rpos + ssize
      end if

    end if

    if ( .not. ma%crmflag ) then
      if ( ma%top /= mpi_proc_null) then
        ssize = nex*jsize
        ib = spos
        do i = 1, nex
          do j = j1, j2
            r4vector1(ib) = ml(j,i2-i+1)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call send_array(r4vector1(spos:spos+ssize), &
                        ssize,ma%top,tag_bt,req(irc))
        spos = spos + ssize
      end if
      if ( ma%bottom /= mpi_proc_null) then
        ssize = nex*jsize
        irc = irc + 1
        call recv_array(r4vector2(rpos:rpos+ssize), &
                        ssize,ma%bottom,tag_bt,req(irc))
        rpos = rpos + ssize
      end if
      if ( ma%topright /= mpi_proc_null ) then
        ssize = nex*nex
        ib = spos
        do i = 1, nex
          do j = 1, nex
            r4vector1(ib) = ml(j2-j+1,i2-i+1)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call send_array(r4vector1(spos:spos+ssize), &
                        ssize,ma%topright,tag_bltr,req(irc))
        spos = spos + ssize
      end if
      if ( ma%bottomleft /= mpi_proc_null ) then
        ssize = nex*nex
        irc = irc + 1
        call recv_array(r4vector2(rpos:rpos+ssize), &
                        ssize,ma%bottomleft,tag_bltr,req(irc))
        rpos = rpos + ssize
      end if
    end if
    if ( irc /= 0 ) then
      ipos = 1
      call mpi_waitall(6,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      if ( .not. ma%bandflag ) then
        if ( ma%left /= mpi_proc_null ) then
          ib = ipos
          ssize = nex * isize
          do i = i1, i2
            do j = 1, nex
              ml(j1-j,i) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
      if ( .not. ma%crmflag ) then
        if ( ma%bottom /= mpi_proc_null) then
          ib = ipos
          ssize = nex * jsize
          do i = 1, nex
            do j = j1, j2
              ml(j,i1-i) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%bottomleft /= mpi_proc_null ) then
          ib = ipos
          ssize = nex * nex
          do i = 1, nex
            do j = 1, nex
              ml(j1-j,i1-i) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
    end if
  end subroutine real4_2d_exchange_left_bottom

  subroutine real4_3d_exchange_left_bottom(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2
    integer(ik4) :: isize, jsize, ksize, ssize, j, i, k, ib
    integer(ik4) :: spos, rpos, ipos, irc
    integer(ik4), dimension(6) :: req
    isize = i2-i1+1
    jsize = j2-j1+1
    ksize = k2-k1+1
    irc = 0
    spos = 1
    rpos = 1
    req = mpi_request_null

    ssize = nex*ksize*(isize+jsize+nex)
    if ( size(r4vector1) < ssize ) then
      call getmem1d(r4vector1,1,ssize,'real4_3d_exchange_left_bottom')
    end if
    if ( size(r4vector2) < ssize ) then
      call getmem1d(r4vector2,1,ssize,'real4_3d_exchange_left_bottom')
    end if

    if ( ma%bandflag ) then
      ssize = nex*isize*ksize
      ib = 1
      do k = k1, k2
        do i = i1, i2
          do j = 1, nex
            r4vector1(ib) = ml(j2-j+1,i,k)
            ib = ib + 1
          end do
        end do
      end do
      call cyclic_exchange_array(r4vector1,r4vector2, &
                                 ssize,ma%right,ma%left,__LINE__)
      ib = 1
      do k = k1, k2
        do i = i1, i2
          do j = 1, nex
            ml(j1-j,i,k) = r4vector2(ib)
            ib = ib + 1
          end do
        end do
      end do

      if ( ma%crmflag ) then
        ! set the size of the dummy exchange vector
        ! nex is the width of the exchange stencil
        ! jsize is the width of a block in the E-W direction
        ssize = nex*jsize*ksize

        ! (1) Send top boundary to bottom side of top neighbor
        ! loop over the given longitudes and number of exchange blocks
        ! and unravel the top boundary into a vector
        ib = 1

        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              r4vector1(ib) = ml(j,i2-i+1,k)
              ib = ib + 1
            end do
          end do
        end do

        ! do the send-down exchange
        call cyclic_exchange_array(r4vector1,r4vector2, &
                                   ssize,ma%top,ma%bottom,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the bottom boundary of this block
        ib = 1

        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              ml(j,i1-i,k) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
        end do

        ! (2) Send top-right boundary to bottom-left side of
        ! top-right neighbor
        ! loop over the exchange block and unravel it into a vector
        ssize = nex*nex*ksize
        ib = 1

        do k = k1, k2
          do i = 1, nex
            do j = 1, nex
              r4vector1(ib) = ml(j2-j+1,i2-i+1,k)
              ib = ib + 1
            end do
          end do
        end do

        ! do the send up-right exchange
        call cyclic_exchange_array(r4vector1,r4vector2, &
                                   ssize,ma%topright,ma%bottomleft,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the bottom-left boundary of this block
        ib = 1

        do k = k1, k2
          do i = 1, nex
            do j = 1, nex
              ml(j1-j,i1-i,k) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
        end do

      end if

    else

      if ( ma%right /= mpi_proc_null) then
        ssize = nex*ksize*isize
        ib = spos
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              r4vector1(ib) = ml(j2-j+1,i,k)
              ib = ib + 1
            end do
          end do
        end do
        irc = irc + 1
        call send_array(r4vector1(spos:spos+ssize), &
                        ssize,ma%right,tag_lr,req(irc))
        spos = spos + ssize
      end if
      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*ksize*isize
        irc = irc + 1
        call recv_array(r4vector2(rpos:rpos+ssize), &
                        ssize,ma%left,tag_lr,req(irc))
        rpos = rpos + ssize
      end if
    end if
    if ( .not. ma%crmflag ) then
      if ( ma%top /= mpi_proc_null) then
        ssize = nex*ksize*jsize
        ib = spos
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              r4vector1(ib) = ml(j,i2-i+1,k)
              ib = ib + 1
            end do
          end do
        end do
        irc = irc + 1
        call send_array(r4vector1(spos:spos+ssize), &
                        ssize,ma%top,tag_bt,req(irc))
        spos = spos + ssize
      end if
      if ( ma%bottom /= mpi_proc_null) then
        ssize = nex*ksize*jsize
        irc = irc + 1
        call recv_array(r4vector2(rpos:rpos+ssize), &
                        ssize,ma%bottom,tag_bt,req(irc))
        rpos = rpos + ssize
      end if
      if ( ma%topright /= mpi_proc_null ) then
        ssize = nex*ksize*nex
        ib = spos
        do k = k1, k2
          do i = 1, nex
            do j = 1, nex
              r4vector1(ib) = ml(j2-j+1,i2-i+1,k)
              ib = ib + 1
            end do
          end do
        end do
        irc = irc + 1
        call send_array(r4vector1(spos:spos+ssize), &
                        ssize,ma%topright,tag_bltr,req(irc))
        spos = spos + ssize
      end if
      if ( ma%bottomleft /= mpi_proc_null ) then
        ssize = nex*ksize*nex
        irc = irc + 1
        call recv_array(r4vector2(rpos:rpos+ssize), &
                        ssize,ma%bottomleft,tag_bltr,req(irc))
        rpos = rpos + ssize
      end if
    end if
    if ( irc /= 0 ) then
      ipos = 1
      call mpi_waitall(6,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      if ( .not. ma%bandflag ) then
        if ( ma%left /= mpi_proc_null ) then
          ib = ipos
          ssize = nex * isize * ksize
          do k = k1, k2
            do i = i1, i2
              do j = 1, nex
                ml(j1-j,i,k) = r4vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
      if ( .not. ma%crmflag ) then
        if ( ma%bottom /= mpi_proc_null) then
          ib = ipos
          ssize = nex * jsize * ksize
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                ml(j,i1-i,k) = r4vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%bottomleft /= mpi_proc_null ) then
          ib = ipos
          ssize = nex * nex * ksize
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                ml(j1-j,i1-i,k) = r4vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
    end if
  end subroutine real4_3d_exchange_left_bottom

  subroutine real4_4d_exchange_left_bottom(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2, n1, n2
    integer(ik4) :: isize, jsize, ksize, nsize, ssize, j, i, k, n, ib
    integer(ik4) :: spos, rpos, ipos, irc
    integer(ik4), dimension(6) :: req
    isize = i2-i1+1
    jsize = j2-j1+1
    ksize = k2-k1+1
    nsize = n2-n1+1
    irc = 0
    spos = 1
    rpos = 1
    req = mpi_request_null

    ssize = nex*ksize*nsize*(isize+jsize+nex)
    if ( size(r4vector1) < ssize ) then
      call getmem1d(r4vector1,1,ssize,'real4_4d_exchange_left_bottom')
    end if
    if ( size(r4vector2) < ssize ) then
      call getmem1d(r4vector2,1,ssize,'real4_4d_exchange_left_bottom')
    end if

    if ( ma%bandflag ) then
      ssize = nex*isize*ksize*nsize
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              r4vector1(ib) = ml(j2-j+1,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      call cyclic_exchange_array(r4vector1,r4vector2, &
                                 ssize,ma%right,ma%left,__LINE__)
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              ml(j1-j,i,k,n) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do

      if ( ma%crmflag ) then
        ! set the size of the dummy exchange vector
        ! nex is the width of the exchange stencil
        ! jsize is the width of a block in the E-W direction
        ssize = nex*jsize*ksize*nsize

        ! (1) Send top boundary to bottom side of top neighbor
        ! loop over the given longitudes and number of exchange blocks
        ! and unravel the top boundary into a vector
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                r4vector1(ib) = ml(j,i2-i+1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        ! do the send-down exchange
        call cyclic_exchange_array(r4vector1,r4vector2, &
                                   ssize,ma%top,ma%bottom,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the bottom boundary of this block
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                ml(j,i1-i,k,n) = r4vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do

        ! (2) Send top-right boundary to bottom-left side of top-right neighbor
        ! loop over the exchange block and unravel it into a vector
        ssize = nex*nex*ksize*nsize
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                r4vector1(ib) = ml(j2-j+1,i2-i+1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        ! do the send up-right exchange
        call cyclic_exchange_array(r4vector1,r4vector2, &
                                   ssize,ma%topright,ma%bottomleft,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the bottom-left boundary of this block
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                ml(j1-j,i1-i,k,n) = r4vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do

      end if

    else

      if ( ma%right /= mpi_proc_null) then
        ssize = nex*ksize*nsize*isize
        ib = spos
        do n = n1, n2
          do k = k1, k2
            do i = i1, i2
              do j = 1, nex
                r4vector1(ib) = ml(j2-j+1,i,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call send_array(r4vector1(spos:spos+ssize), &
                        ssize,ma%right,tag_lr,req(irc))
        spos = spos + ssize
      end if
      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*ksize*nsize*isize
        irc = irc + 1
        call recv_array(r4vector2(rpos:rpos+ssize), &
                        ssize,ma%left,tag_lr,req(irc))
        rpos = rpos + ssize
      end if
    end if
    if ( .not. ma%crmflag ) then
      if ( ma%top /= mpi_proc_null) then
        ssize = nex*ksize*nsize*jsize
        ib = spos
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                r4vector1(ib) = ml(j,i2-i+1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call send_array(r4vector1(spos:spos+ssize), &
                        ssize,ma%top,tag_bt,req(irc))
        spos = spos + ssize
      end if
      if ( ma%bottom /= mpi_proc_null) then
        ssize = nex*ksize*nsize*jsize
        irc = irc + 1
        call recv_array(r4vector2(rpos:rpos+ssize), &
                        ssize,ma%bottom,tag_bt,req(irc))
        rpos = rpos + ssize
      end if
      if ( ma%topright /= mpi_proc_null ) then
        ssize = nex*ksize*nsize*nex
        ib = spos
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                r4vector1(ib) = ml(j2-j+1,i2-i+1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call send_array(r4vector1(spos:spos+ssize), &
                        ssize,ma%topright,tag_bltr,req(irc))
        spos = spos + ssize
      end if
      if ( ma%bottomleft /= mpi_proc_null ) then
        ssize = nex*ksize*nsize*nex
        irc = irc + 1
        call recv_array(r4vector2(rpos:rpos+ssize), &
                        ssize,ma%bottomleft,tag_bltr,req(irc))
        rpos = rpos + ssize
      end if
    end if
    if ( irc /= 0 ) then
      ipos = 1
      call mpi_waitall(6,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      if ( .not. ma%bandflag ) then
        if ( ma%left /= mpi_proc_null ) then
          ib = ipos
          ssize = nex * isize * ksize * nsize
          do n = n1, n2
            do k = k1, k2
              do i = i1, i2
                do j = 1, nex
                  ml(j1-j,i,k,n) = r4vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
      if ( .not. ma%crmflag ) then
        if ( ma%bottom /= mpi_proc_null) then
          ib = ipos
          ssize = nex * jsize * ksize * nsize
          do n = n1, n2
            do k = k1, k2
              do i = 1, nex
                do j = j1, j2
                  ml(j,i1-i,k,n) = r4vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%bottomleft /= mpi_proc_null ) then
          ib = ipos
          ssize = nex * nex * ksize * nsize
          do n = n1, n2
            do k = k1, k2
              do i = 1, nex
                do j = 1, nex
                  ml(j1-j,i1-i,k,n) = r4vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
    end if
  end subroutine real4_4d_exchange_left_bottom

  subroutine real8_2d_exchange_right_top(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2
    integer(ik4) :: isize, jsize, ssize, j, i, ib
    integer(ik4) :: spos, rpos, ipos, irc
    integer(ik4), dimension(6) :: req
    isize = i2-i1+1
    jsize = j2-j1+1
    irc = 0
    spos = 1
    rpos = 1
    req = mpi_request_null

    ssize = nex*(isize+jsize+nex)
    if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'real8_2d_exchange_right_top')
    end if
    if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'real8_2d_exchange_right_top')
    end if

    if ( ma%bandflag ) then
      ssize = nex*isize
      ib = 1
      do i = i1, i2
        do j = 1, nex
          r8vector1(ib) = ml(j1+j-1,i)
          ib = ib + 1
        end do
      end do
      call cyclic_exchange_array(r8vector1,r8vector2, &
                                 ssize,ma%left,ma%right,__LINE__)
      ib = 1
      do i = i1, i2
        do j = 1, nex
          ml(j2+j,i) = r8vector2(ib)
          ib = ib + 1
        end do
      end do

      if ( ma%crmflag ) then
        ! (1) Send bottom boundary to top side of bottom neighbor
        ! nex is the width of the exchange stencil
        ! jsize is the width of a block in the E-W direction
        ssize = nex*jsize
        ! loop over the given longitudes and number of exchange blocks
        ! and unravel the bottom boundary into a vector
        ib = 1
        do i = 1, nex
          do j = j1, j2
            r8vector1(ib) = ml(j,i1+i-1)
            ib = ib + 1
          end do
        end do

        ! do the send-up exchange
        call cyclic_exchange_array(r8vector1,r8vector2, &
                                   ssize,ma%bottom,ma%top,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the top boundary of this block
        ib = 1

        do i = 1, nex
          do j = j1, j2
            ml(j,i2+i) = r8vector2(ib)
            ib = ib + 1
          end do
        end do

        ! (2) Send bottom-left boundary to top-right side of
        ! bottom-left neighbor
        ! loop over the exchange block and unravel it into a vector
        ssize = nex*nex
        ib = 1

        do i = 1, nex
          do j = 1, nex
            r8vector1(ib) = ml(j1+j-1,i1+i-1)
            ib = ib + 1
          end do
        end do

        ! do the send down-left exchange
        call cyclic_exchange_array(r8vector1,r8vector2, &
                                    ssize,ma%bottomleft,ma%topright,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the top-right boundary of this block
        ib = 1

        do i = 1, nex
          do j = 1, nex
            ml(j2+j,i2+i) = r8vector2(ib)
            ib = ib + 1
          end do
        end do

      end if

    else

      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*isize
        ib = spos
        do i = i1, i2
          do j = 1, nex
            r8vector1(ib) = ml(j1+j-1,i)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call send_array(r8vector1(spos:spos+ssize), &
                        ssize,ma%left,tag_rl,req(irc))
        spos = spos + ssize
      end if
      if ( ma%right /= mpi_proc_null) then
        ssize = nex*isize
        irc = irc + 1
        call recv_array(r8vector2(rpos:rpos+ssize), &
                        ssize,ma%right,tag_rl,req(irc))
        rpos = rpos + ssize
      end if
    end if

    if ( .not. ma%crmflag ) then
      if ( ma%bottom /= mpi_proc_null) then
        ssize = nex*jsize
        ib = spos
        do i = 1, nex
          do j = j1, j2
            r8vector1(ib) = ml(j,i1+i-1)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call send_array(r8vector1(spos:spos+ssize), &
                        ssize,ma%bottom,tag_tb,req(irc))
        spos = spos + ssize
      end if
      if ( ma%top /= mpi_proc_null) then
        ssize = nex*jsize
        irc = irc + 1
        call recv_array(r8vector2(rpos:rpos+ssize), &
                        ssize,ma%top,tag_tb,req(irc))
        rpos = rpos + ssize
      end if
      if ( ma%bottomleft /= mpi_proc_null ) then
        ssize = nex*nex
        ib = spos
        do i = 1, nex
          do j = 1, nex
            r8vector1(ib) = ml(j1+j-1,i1+i-1)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call send_array(r8vector1(spos:spos+ssize), &
                        ssize,ma%bottomleft,tag_trbl,req(irc))
        spos = spos + ssize
      end if
      if ( ma%topright /= mpi_proc_null ) then
        ssize = nex*nex
        irc = irc + 1
        call recv_array(r8vector2(rpos:rpos+ssize), &
                        ssize,ma%topright,tag_trbl,req(irc))
        rpos = rpos + ssize
      end if
    end if
    if ( irc /= 0 ) then
      ipos = 1
      call mpi_waitall(6,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      if ( .not. ma%bandflag ) then
        if ( ma%right /= mpi_proc_null) then
          ssize = nex*isize
          ib = ipos
          do i = i1, i2
            do j = 1, nex
              ml(j2+j,i) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
      if ( .not. ma%crmflag ) then
        if ( ma%top /= mpi_proc_null) then
          ssize = nex*jsize
          ib = ipos
          do i = 1, nex
            do j = j1, j2
              ml(j,i2+i) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%topright /= mpi_proc_null ) then
          ssize = nex*nex
          ib = ipos
          do i = 1, nex
            do j = 1, nex
              ml(j2+j,i2+i) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
    end if
  end subroutine real8_2d_exchange_right_top

  subroutine real8_3d_exchange_right_top(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2
    integer(ik4) :: isize, jsize, ksize, ssize, j, i, k, ib
    integer(ik4) :: spos, rpos, ipos, irc
    integer(ik4), dimension(6) :: req
    isize = i2-i1+1
    jsize = j2-j1+1
    ksize = k2-k1+1
    irc = 0
    spos = 1
    rpos = 1
    req = mpi_request_null

    ssize = nex*ksize*(isize+jsize+nex)
    if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'real8_3d_exchange_right_top')
    end if
    if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'real8_3d_exchange_right_top')
    end if

    if ( ma%bandflag ) then
      ssize = nex*isize*ksize
      ib = 1
      do k = k1, k2
        do i = i1, i2
          do j = 1, nex
            r8vector1(ib) = ml(j1+j-1,i,k)
            ib = ib + 1
          end do
        end do
      end do
      call cyclic_exchange_array(r8vector1,r8vector2, &
                                 ssize,ma%left,ma%right,__LINE__)
      ib = 1
      do k = k1, k2
        do i = i1, i2
          do j = 1, nex
            ml(j2+j,i,k) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end do

      if ( ma%crmflag ) then
        ! (1) Send bottom boundary to top side of bottom neighbor
        ! nex is the width of the exchange stencil
        ! jsize is the width of a block in the E-W direction
        ssize = nex*jsize*ksize
        ! loop over the given longitudes and number of exchange blocks
        ! and unravel the bottom boundary into a vector
        ib = 1
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              r8vector1(ib) = ml(j,i1+i-1,k)
              ib = ib + 1
            end do
          end do
        end do

        ! do the send-up exchange
        call cyclic_exchange_array(r8vector1,r8vector2, &
                                   ssize,ma%bottom,ma%top,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the top boundary of this block
        ib = 1
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              ml(j,i2+i,k) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do

        ! (2) Send bottom-left boundary to top-right side of
        ! bottom-left neighbor
        ! loop over the exchange block and unravel it into a vector
        ssize = nex*nex*ksize
        ib = 1
        do k = k1, k2
          do i = 1, nex
            do j = 1, nex
              r8vector1(ib) = ml(j1+j-1,i1+i-1,k)
              ib = ib + 1
            end do
          end do
        end do

        ! do the send down-left exchange
        call cyclic_exchange_array(r8vector1,r8vector2, &
                                    ssize,ma%bottomleft,ma%topright,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the top-right boundary of this block
        ib = 1
        do k = k1, k2
          do i = 1, nex
            do j = 1, nex
              ml(j2+j,i2+i,k) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do

      end if

    else

      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*isize*ksize
        ib = spos
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              r8vector1(ib) = ml(j1+j-1,i,k)
              ib = ib + 1
            end do
          end do
        end do
        irc = irc + 1
        call send_array(r8vector1(spos:spos+ssize), &
                        ssize,ma%left,tag_rl,req(irc))
        spos = spos + ssize
      end if
      if ( ma%right /= mpi_proc_null) then
        ssize = nex*isize*ksize
        irc = irc + 1
        call recv_array(r8vector2(rpos:rpos+ssize), &
                        ssize,ma%right,tag_rl,req(irc))
        rpos = rpos + ssize
      end if
    end if

    if ( .not. ma%crmflag ) then
      if ( ma%bottom /= mpi_proc_null) then
        ssize = nex*jsize*ksize
        ib = spos
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              r8vector1(ib) = ml(j,i1+i-1,k)
              ib = ib + 1
            end do
          end do
        end do
        irc = irc + 1
        call send_array(r8vector1(spos:spos+ssize), &
                        ssize,ma%bottom,tag_tb,req(irc))
        spos = spos + ssize
      end if
      if ( ma%top /= mpi_proc_null) then
        ssize = nex*jsize*ksize
        irc = irc + 1
        call recv_array(r8vector2(rpos:rpos+ssize), &
                        ssize,ma%top,tag_tb,req(irc))
        rpos = rpos + ssize
      end if
      if ( ma%bottomleft /= mpi_proc_null ) then
        ssize = nex*nex*ksize
        ib = spos
        do k = k1, k2
          do i = 1, nex
            do j = 1, nex
              r8vector1(ib) = ml(j1+j-1,i1+i-1,k)
              ib = ib + 1
            end do
          end do
        end do
        irc = irc + 1
        call send_array(r8vector1(spos:spos+ssize), &
                        ssize,ma%bottomleft,tag_trbl,req(irc))
        spos = spos + ssize
      end if
      if ( ma%topright /= mpi_proc_null ) then
        ssize = nex*nex*ksize
        irc = irc + 1
        call recv_array(r8vector2(rpos:rpos+ssize), &
                        ssize,ma%topright,tag_trbl,req(irc))
        rpos = rpos + ssize
      end if
    end if
    if ( irc /= 0 ) then
      ipos = 1
      call mpi_waitall(6,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      if ( .not. ma%bandflag ) then
        if ( ma%right /= mpi_proc_null) then
          ssize = nex*isize*ksize
          ib = ipos
          do k = k1, k2
            do i = i1, i2
              do j = 1, nex
                ml(j2+j,i,k) = r8vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
      if ( .not. ma%crmflag ) then
        if ( ma%top /= mpi_proc_null) then
          ssize = nex*jsize*ksize
          ib = ipos
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                ml(j,i2+i,k) = r8vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%topright /= mpi_proc_null ) then
          ssize = nex*nex*ksize
          ib = ipos
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                ml(j2+j,i2+i,k) = r8vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
    end if
  end subroutine real8_3d_exchange_right_top

  subroutine real8_4d_exchange_right_top(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2, n1, n2
    integer(ik4) :: isize, jsize, ksize, nsize, ssize, j, i, k, n, ib
    integer(ik4) :: spos, rpos, ipos, irc
    integer(ik4), dimension(6) :: req
    isize = i2-i1+1
    jsize = j2-j1+1
    ksize = k2-k1+1
    nsize = n2-n1+1
    irc = 0
    spos = 1
    rpos = 1
    req = mpi_request_null

    ssize = nex*ksize*nsize*(isize+jsize+nex)
    if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'real8_4d_exchange_right_top')
    end if
    if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'real8_4d_exchange_right_top')
    end if

    if ( ma%bandflag ) then
      ssize = nex*isize*ksize*nsize
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              r8vector1(ib) = ml(j1+j-1,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      call cyclic_exchange_array(r8vector1,r8vector2, &
                                 ssize,ma%left,ma%right,__LINE__)
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              ml(j2+j,i,k,n) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do

      if ( ma%crmflag ) then
        ! (1) Send bottom boundary to top side of bottom neighbor
        ! nex is the width of the exchange stencil
        ! jsize is the width of a block in the E-W direction
        ssize = nex*jsize*ksize*nsize
        ! loop over the given longitudes and number of exchange blocks
        ! and unravel the bottom boundary into a vector
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                r8vector1(ib) = ml(j,i1+i-1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        ! do the send-up exchange
        call cyclic_exchange_array(r8vector1,r8vector2, &
                                   ssize,ma%bottom,ma%top,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the top boundary of this block
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                ml(j,i2+i,k,n) = r8vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do

        ! (2) Send bottom-left boundary to top-right side
        !     of bottom-left neighbor
        ! loop over the exchange block and unravel it into a vector
        ssize = nex*nex*ksize*nsize
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                r8vector1(ib) = ml(j1+j-1,i1+i-1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        ! do the send down-left exchange
        call cyclic_exchange_array(r8vector1,r8vector2, &
                                    ssize,ma%bottomleft,ma%topright,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the top-right boundary of this block
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                ml(j2+j,i2+i,k,n) = r8vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do

      end if

    else

      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*isize*ksize*nsize
        ib = spos
        do n = n1, n2
          do k = k1, k2
            do i = i1, i2
              do j = 1, nex
                r8vector1(ib) = ml(j1+j-1,i,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call send_array(r8vector1(spos:spos+ssize), &
                        ssize,ma%left,tag_rl,req(irc))
        spos = spos + ssize
      end if
      if ( ma%right /= mpi_proc_null) then
        ssize = nex*isize*ksize*nsize
        irc = irc + 1
        call recv_array(r8vector2(rpos:rpos+ssize), &
                        ssize,ma%right,tag_rl,req(irc))
        rpos = rpos + ssize
      end if
    end if

    if ( .not. ma%crmflag ) then
      if ( ma%bottom /= mpi_proc_null) then
        ssize = nex*jsize*ksize*nsize
        ib = spos
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                r8vector1(ib) = ml(j,i1+i-1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call send_array(r8vector1(spos:spos+ssize), &
                        ssize,ma%bottom,tag_tb,req(irc))
        spos = spos + ssize
      end if
      if ( ma%top /= mpi_proc_null) then
        ssize = nex*jsize*ksize*nsize
        irc = irc + 1
        call recv_array(r8vector2(rpos:rpos+ssize), &
                        ssize,ma%top,tag_tb,req(irc))
        rpos = rpos + ssize
      end if
      if ( ma%bottomleft /= mpi_proc_null ) then
        ssize = nex*nex*ksize*nsize
        ib = spos
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                r8vector1(ib) = ml(j1+j-1,i1+i-1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call send_array(r8vector1(spos:spos+ssize), &
                        ssize,ma%bottomleft,tag_trbl,req(irc))
        spos = spos + ssize
      end if
      if ( ma%topright /= mpi_proc_null ) then
        ssize = nex*nex*ksize*nsize
        irc = irc + 1
        call recv_array(r8vector2(rpos:rpos+ssize), &
                        ssize,ma%topright,tag_trbl,req(irc))
        rpos = rpos + ssize
      end if
    end if
    if ( irc /= 0 ) then
      ipos = 1
      call mpi_waitall(6,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      if ( .not. ma%bandflag ) then
        if ( ma%right /= mpi_proc_null) then
          ssize = nex*isize*ksize*nsize
          ib = ipos
          do n = n1, n2
            do k = k1, k2
              do i = i1, i2
                do j = 1, nex
                  ml(j2+j,i,k,n) = r8vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
      if ( .not. ma%crmflag ) then
        if ( ma%top /= mpi_proc_null) then
          ssize = nex*jsize*ksize*nsize
          ib = ipos
          do n = n1, n2
            do k = k1, k2
              do i = 1, nex
                do j = j1, j2
                  ml(j,i2+i,k,n) = r8vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%topright /= mpi_proc_null ) then
          ssize = nex*nex*ksize*nsize
          ib = ipos
          do n = n1, n2
            do k = k1, k2
              do i = 1, nex
                do j = 1, nex
                  ml(j2+j,i2+i,k,n) = r8vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
    end if
  end subroutine real8_4d_exchange_right_top

  subroutine real4_2d_exchange_right_top(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2
    integer(ik4) :: isize, jsize, ssize, j, i, ib
    integer(ik4) :: spos, rpos, ipos, irc
    integer(ik4), dimension(6) :: req
    isize = i2-i1+1
    jsize = j2-j1+1
    irc = 0
    spos = 1
    rpos = 1
    req = mpi_request_null

    ssize = nex*(isize+jsize+nex)
    if ( size(r4vector1) < ssize ) then
      call getmem1d(r4vector1,1,ssize,'real4_2d_exchange_right_top')
    end if
    if ( size(r4vector2) < ssize ) then
      call getmem1d(r4vector2,1,ssize,'real4_2d_exchange_right_top')
    end if

    if ( ma%bandflag ) then
      ssize = nex*isize
      ib = 1
      do i = i1, i2
        do j = 1, nex
          r4vector1(ib) = ml(j1+j-1,i)
          ib = ib + 1
        end do
      end do
      call cyclic_exchange_array(r4vector1,r4vector2, &
                                 ssize,ma%left,ma%right,__LINE__)
      ib = 1
      do i = i1, i2
        do j = 1, nex
          ml(j2+j,i) = r4vector2(ib)
          ib = ib + 1
        end do
      end do

      if ( ma%crmflag ) then
        ! (1) Send bottom boundary to top side of bottom neighbor
        ! nex is the width of the exchange stencil
        ! jsize is the width of a block in the E-W direction
        ssize = nex*jsize
        ! loop over the given longitudes and number of exchange blocks
        ! and unravel the bottom boundary into a vector
        ib = 1
        do i = 1, nex
          do j = j1, j2
            r4vector1(ib) = ml(j,i1+i-1)
            ib = ib + 1
          end do
        end do

        ! do the send-up exchange
        call cyclic_exchange_array(r4vector1,r4vector2, &
                                   ssize,ma%bottom,ma%top,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the top boundary of this block
        ib = 1

        do i = 1, nex
          do j = j1, j2
            ml(j,i2+i) = r4vector2(ib)
            ib = ib + 1
          end do
        end do

        ! (2) Send bottom-left boundary to top-right side of
        ! bottom-left neighbor
        ! loop over the exchange block and unravel it into a vector
        ssize = nex*nex
        ib = 1

        do i = 1, nex
          do j = 1, nex
            r4vector1(ib) = ml(j1+j-1,i1+i-1)
            ib = ib + 1
          end do
        end do

        ! do the send down-left exchange
        call cyclic_exchange_array(r4vector1,r4vector2, &
                                    ssize,ma%bottomleft,ma%topright,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the top-right boundary of this block
        ib = 1
        do i = 1, nex
          do j = 1, nex
            ml(j2+j,i2+i) = r4vector2(ib)
            ib = ib + 1
          end do
        end do

      end if

    else

      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*isize
        ib = spos
        do i = i1, i2
          do j = 1, nex
            r4vector1(ib) = ml(j1+j-1,i)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call send_array(r4vector1(spos:spos+ssize), &
                        ssize,ma%left,tag_rl,req(irc))
        spos = spos + ssize
      end if
      if ( ma%right /= mpi_proc_null) then
        ssize = nex*isize
        irc = irc + 1
        call recv_array(r4vector2(rpos:rpos+ssize), &
                        ssize,ma%right,tag_rl,req(irc))
        rpos = rpos + ssize
      end if
    end if

    if ( .not. ma%crmflag ) then
      if ( ma%bottom /= mpi_proc_null) then
        ssize = nex*jsize
        ib = spos
        do i = 1, nex
          do j = j1, j2
            r4vector1(ib) = ml(j,i1+i-1)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call send_array(r4vector1(spos:spos+ssize), &
                        ssize,ma%bottom,tag_tb,req(irc))
        spos = spos + ssize
      end if
      if ( ma%top /= mpi_proc_null) then
        ssize = nex*jsize
        irc = irc + 1
        call recv_array(r4vector2(rpos:rpos+ssize), &
                        ssize,ma%top,tag_tb,req(irc))
        rpos = rpos + ssize
      end if
      if ( ma%bottomleft /= mpi_proc_null ) then
        ssize = nex*nex
        ib = spos
        do i = 1, nex
          do j = 1, nex
            r4vector1(ib) = ml(j1+j-1,i1+i-1)
            ib = ib + 1
          end do
        end do
        irc = irc + 1
        call send_array(r4vector1(spos:spos+ssize), &
                        ssize,ma%bottomleft,tag_trbl,req(irc))
        spos = spos + ssize
      end if
      if ( ma%topright /= mpi_proc_null ) then
        ssize = nex*nex
        irc = irc + 1
        call recv_array(r4vector2(rpos:rpos+ssize), &
                        ssize,ma%topright,tag_trbl,req(irc))
        rpos = rpos + ssize
      end if
    end if
    if ( irc /= 0 ) then
      ipos = 1
      call mpi_waitall(6,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      if ( .not. ma%bandflag ) then
        if ( ma%right /= mpi_proc_null) then
          ssize = nex*isize
          ib = ipos
          do i = i1, i2
            do j = 1, nex
              ml(j2+j,i) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
      if ( .not. ma%crmflag ) then
        if ( ma%top /= mpi_proc_null) then
          ssize = nex*jsize
          ib = ipos
          do i = 1, nex
            do j = j1, j2
              ml(j,i2+i) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%topright /= mpi_proc_null ) then
          ssize = nex*nex
          ib = ipos
          do i = 1, nex
            do j = 1, nex
              ml(j2+j,i2+i) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
    end if
  end subroutine real4_2d_exchange_right_top

  subroutine real4_3d_exchange_right_top(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2
    integer(ik4) :: isize, jsize, ksize, ssize, j, i, k, ib
    integer(ik4) :: spos, rpos, ipos, irc
    integer(ik4), dimension(6) :: req
    isize = i2-i1+1
    jsize = j2-j1+1
    ksize = k2-k1+1
    irc = 0
    spos = 1
    rpos = 1
    req = mpi_request_null

    ssize = nex*ksize*(isize+jsize+nex)
    if ( size(r4vector1) < ssize ) then
      call getmem1d(r4vector1,1,ssize,'real4_3d_exchange_right_top')
    end if
    if ( size(r4vector2) < ssize ) then
      call getmem1d(r4vector2,1,ssize,'real4_3d_exchange_right_top')
    end if

    if ( ma%bandflag ) then
      ssize = nex*isize*ksize
      ib = 1
      do k = k1, k2
        do i = i1, i2
          do j = 1, nex
            r4vector1(ib) = ml(j1+j-1,i,k)
            ib = ib + 1
          end do
        end do
      end do
      call cyclic_exchange_array(r4vector1,r4vector2, &
                                 ssize,ma%left,ma%right,__LINE__)
      ib = 1
      do k = k1, k2
        do i = i1, i2
          do j = 1, nex
            ml(j2+j,i,k) = r4vector2(ib)
            ib = ib + 1
          end do
        end do
      end do

      if ( ma%crmflag ) then
        ! (1) Send bottom boundary to top side of bottom neighbor
        ! nex is the width of the exchange stencil
        ! jsize is the width of a block in the E-W direction
        ssize = nex*jsize*ksize
        ! loop over the given longitudes and number of exchange blocks
        ! and unravel the bottom boundary into a vector
        ib = 1
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              r4vector1(ib) = ml(j,i1+i-1,k)
              ib = ib + 1
            end do
          end do
        end do

        ! do the send-up exchange
        call cyclic_exchange_array(r4vector1,r4vector2, &
                                   ssize,ma%bottom,ma%top,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the top boundary of this block
        ib = 1
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              ml(j,i2+i,k) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
        end do

        ! (2) Send bottom-left boundary to top-right side of
        ! bottom-left neighbor
        ! loop over the exchange block and unravel it into a vector
        ssize = nex*nex*ksize
        ib = 1
        do k = k1, k2
          do i = 1, nex
            do j = 1, nex
              r4vector1(ib) = ml(j1+j-1,i1+i-1,k)
              ib = ib + 1
            end do
          end do
        end do

        ! do the send down-left exchange
        call cyclic_exchange_array(r4vector1,r4vector2, &
                                    ssize,ma%bottomleft,ma%topright,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the top-right boundary of this block
        ib = 1
        do k = k1, k2
          do i = 1, nex
            do j = 1, nex
              ml(j2+j,i2+i,k) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
        end do

      end if

    else

      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*isize*ksize
        ib = spos
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              r4vector1(ib) = ml(j1+j-1,i,k)
              ib = ib + 1
            end do
          end do
        end do
        irc = irc + 1
        call send_array(r4vector1(spos:spos+ssize), &
                        ssize,ma%left,tag_rl,req(irc))
        spos = spos + ssize
      end if
      if ( ma%right /= mpi_proc_null) then
        ssize = nex*isize*ksize
        irc = irc + 1
        call recv_array(r4vector2(rpos:rpos+ssize), &
                        ssize,ma%right,tag_rl,req(irc))
        rpos = rpos + ssize
      end if
    end if

    if ( .not. ma%crmflag ) then
      if ( ma%bottom /= mpi_proc_null) then
        ssize = nex*jsize*ksize
        ib = spos
        do k = k1, k2
          do i = 1, nex
            do j = j1, j2
              r4vector1(ib) = ml(j,i1+i-1,k)
              ib = ib + 1
            end do
          end do
        end do
        irc = irc + 1
        call send_array(r4vector1(spos:spos+ssize), &
                        ssize,ma%bottom,tag_tb,req(irc))
        spos = spos + ssize
      end if
      if ( ma%top /= mpi_proc_null) then
        ssize = nex*jsize*ksize
        irc = irc + 1
        call recv_array(r4vector2(rpos:rpos+ssize), &
                        ssize,ma%top,tag_tb,req(irc))
        rpos = rpos + ssize
      end if
      if ( ma%bottomleft /= mpi_proc_null ) then
        ssize = nex*nex*ksize
        ib = spos
        do k = k1, k2
          do i = 1, nex
            do j = 1, nex
              r4vector1(ib) = ml(j1+j-1,i1+i-1,k)
              ib = ib + 1
            end do
          end do
        end do
        irc = irc + 1
        call send_array(r4vector1(spos:spos+ssize), &
                        ssize,ma%bottomleft,tag_trbl,req(irc))
        spos = spos + ssize
      end if
      if ( ma%topright /= mpi_proc_null ) then
        ssize = nex*nex*ksize
        irc = irc + 1
        call recv_array(r4vector2(rpos:rpos+ssize), &
                        ssize,ma%topright,tag_trbl,req(irc))
        rpos = rpos + ssize
      end if
    end if
    if ( irc /= 0 ) then
      ipos = 1
      call mpi_waitall(6,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      if ( .not. ma%bandflag ) then
        if ( ma%right /= mpi_proc_null) then
          ssize = nex*isize*ksize
          ib = ipos
          do k = k1, k2
            do i = i1, i2
              do j = 1, nex
                ml(j2+j,i,k) = r4vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
      if ( .not. ma%crmflag ) then
        if ( ma%top /= mpi_proc_null) then
          ssize = nex*jsize*ksize
          ib = ipos
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                ml(j,i2+i,k) = r4vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%topright /= mpi_proc_null ) then
          ssize = nex*nex*ksize
          ib = ipos
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                ml(j2+j,i2+i,k) = r4vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
    end if
  end subroutine real4_3d_exchange_right_top

  subroutine real4_4d_exchange_right_top(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: nex, j1, j2 , i1, i2, k1, k2, n1, n2
    integer(ik4) :: isize, jsize, ksize, nsize, ssize, j, i, k, n, ib
    integer(ik4) :: spos, rpos, ipos, irc
    integer(ik4), dimension(6) :: req
    isize = i2-i1+1
    jsize = j2-j1+1
    ksize = k2-k1+1
    nsize = n2-n1+1
    irc = 0
    spos = 1
    rpos = 1
    req = mpi_request_null

    ssize = nex*ksize*nsize*(isize+jsize+nex)
    if ( size(r4vector1) < ssize ) then
      call getmem1d(r4vector1,1,ssize,'real4_4d_exchange_right_top')
    end if
    if ( size(r4vector2) < ssize ) then
      call getmem1d(r4vector2,1,ssize,'real4_4d_exchange_right_top')
    end if

    if ( ma%bandflag ) then
      ssize = nex*isize*ksize*nsize
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              r4vector1(ib) = ml(j1+j-1,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      call cyclic_exchange_array(r4vector1,r4vector2, &
                                 ssize,ma%left,ma%right,__LINE__)
      ib = 1
      do n = n1, n2
        do k = k1, k2
          do i = i1, i2
            do j = 1, nex
              ml(j2+j,i,k,n) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do

      if ( ma%crmflag ) then
        ! (1) Send bottom boundary to top side of bottom neighbor
        ! nex is the width of the exchange stencil
        ! jsize is the width of a block in the E-W direction
        ssize = nex*jsize*ksize*nsize
        ! loop over the given longitudes and number of exchange blocks
        ! and unravel the bottom boundary into a vector
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                r4vector1(ib) = ml(j,i1+i-1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        ! do the send-up exchange
        call cyclic_exchange_array(r4vector1,r4vector2, &
                                   ssize,ma%bottom,ma%top,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the top boundary of this block
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                ml(j,i2+i,k,n) = r4vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do

        ! (2) Send bottom-left boundary to top-right side
        !     of bottom-left neighbor
        ! loop over the exchange block and unravel it into a vector
        ssize = nex*nex*ksize*nsize
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                r4vector1(ib) = ml(j1+j-1,i1+i-1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        ! do the send down-left exchange
        call cyclic_exchange_array(r4vector1,r4vector2, &
                                    ssize,ma%bottomleft,ma%topright,__LINE__)
        ! loop over the receiving dummy vector and ravel it into
        ! the top-right boundary of this block
        ib = 1
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                ml(j2+j,i2+i,k,n) = r4vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do

      end if

    else

      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*isize*ksize*nsize
        ib = spos
        do n = n1, n2
          do k = k1, k2
            do i = i1, i2
              do j = 1, nex
                r4vector1(ib) = ml(j1+j-1,i,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call send_array(r4vector1(spos:spos+ssize), &
                        ssize,ma%left,tag_rl,req(irc))
        spos = spos + ssize
      end if
      if ( ma%right /= mpi_proc_null) then
        ssize = nex*isize*ksize*nsize
        irc = irc + 1
        call recv_array(r4vector2(rpos:rpos+ssize), &
                        ssize,ma%right,tag_rl,req(irc))
        rpos = rpos + ssize
      end if
    end if

    if ( .not. ma%crmflag ) then
      if ( ma%bottom /= mpi_proc_null) then
        ssize = nex*jsize*ksize*nsize
        ib = spos
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = j1, j2
                r4vector1(ib) = ml(j,i1+i-1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call send_array(r4vector1(spos:spos+ssize), &
                        ssize,ma%bottom,tag_tb,req(irc))
        spos = spos + ssize
      end if
      if ( ma%top /= mpi_proc_null) then
        ssize = nex*jsize*ksize*nsize
        irc = irc + 1
        call recv_array(r4vector2(rpos:rpos+ssize), &
                        ssize,ma%top,tag_tb,req(irc))
        rpos = rpos + ssize
      end if
      if ( ma%bottomleft /= mpi_proc_null ) then
        ssize = nex*nex*ksize*nsize
        ib = spos
        do n = n1, n2
          do k = k1, k2
            do i = 1, nex
              do j = 1, nex
                r4vector1(ib) = ml(j1+j-1,i1+i-1,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        irc = irc + 1
        call send_array(r4vector1(spos:spos+ssize), &
                        ssize,ma%bottomleft,tag_trbl,req(irc))
        spos = spos + ssize
      end if
      if ( ma%topright /= mpi_proc_null ) then
        ssize = nex*nex*ksize*nsize
        irc = irc + 1
        call recv_array(r4vector2(rpos:rpos+ssize), &
                        ssize,ma%topright,tag_trbl,req(irc))
        rpos = rpos + ssize
      end if
    end if
    if ( irc /= 0 ) then
      ipos = 1
      call mpi_waitall(6,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      if ( .not. ma%bandflag ) then
        if ( ma%right /= mpi_proc_null) then
          ssize = nex*isize*ksize*nsize
          ib = ipos
          do n = n1, n2
            do k = k1, k2
              do i = i1, i2
                do j = 1, nex
                  ml(j2+j,i,k,n) = r4vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%top /= mpi_proc_null) then
          ssize = nex*jsize*ksize*nsize
          ib = ipos
          do n = n1, n2
            do k = k1, k2
              do i = 1, nex
                do j = j1, j2
                  ml(j,i2+i,k,n) = r4vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
        if ( ma%topright /= mpi_proc_null ) then
          ssize = nex*nex*ksize*nsize
          ib = ipos
          do n = n1, n2
            do k = k1, k2
              do i = 1, nex
                do j = 1, nex
                  ml(j2+j,i2+i,k,n) = r4vector2(ib)
                  ib = ib + 1
                end do
              end do
            end do
          end do
          ipos = ipos + ssize
        end if
      end if
    end if
  end subroutine real4_4d_exchange_right_top

#endif

  subroutine real8_bdy_exchange_left_right(ml,k1,k2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: k1, k2
    integer(ik4) :: ssize, ksize, k, ib, irc, ipos
    integer(ik4), dimension(4) :: req
    req =mpi_request_null
    ksize = k2-k1+1
    irc = 0
    ipos = 1
    ssize = 2*ksize
    if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'real8_bdy_exchange_left_right')
    end if
    if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'real8_bdy_exchange_left_right')
    end if
    if ( ma%bandflag ) then
      ib = 1
      do k = k1, k2
        r8vector1(ib) = ml(jde2,k)
        ib = ib + 1
      end do
      call cyclic_exchange_array(r8vector1,r8vector2, &
                                 ksize,ma%left,ma%right,__LINE__)
      ib = 1
      do k = k1, k2
        ml(jde1-1,k) = r8vector2(ib)
        ib = ib + 1
      end do
      ib = 1
      do k = k1, k2
        r8vector1(ib) = ml(jde1,k)
        ib = ib + 1
      end do
      call cyclic_exchange_array(r8vector1,r8vector2, &
                                 ksize,ma%right,ma%left,__LINE__)
      ib = 1
      do k = k1, k2
        ml(jde2+1,k) = r8vector2(ib)
        ib = ib + 1
      end do
    else
      if ( ma%right /= mpi_proc_null) then
        ib = ipos
        do k = k1, k2
          r8vector1(ib) = ml(jde2,k)
          ib = ib + 1
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ksize), &
                            r8vector2(ipos:ipos+ksize),ksize, &
                            ma%right,tag_rl,tag_lr, &
                            req(irc),req(irc+2))
        ipos = ipos + ksize
      end if
      if ( ma%left /= mpi_proc_null ) then
        ib = ipos
        do k = k1, k2
          r8vector1(ib) = ml(jde1,k)
          ib = ib + 1
        end do
        irc = irc + 1
        call exchange_array(r8vector1(ipos:ipos+ksize), &
                            r8vector2(ipos:ipos+ksize),ksize, &
                            ma%left,tag_lr,tag_rl, &
                            req(irc),req(irc+2))
        ipos = ipos + ksize
      end if
      if ( irc /= 0 ) then
        call mpi_waitall(4,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
        if ( mpierr /= mpi_success ) then
          call fatal(__FILE__,__LINE__,'mpi_waitall error.')
        end if
#endif
        ipos = 1
        if ( ma%right /= mpi_proc_null) then
          ib = ipos
          do k = k1, k2
            ml(jde2+1,k) = r8vector2(ib)
            ib = ib + 1
          end do
          ipos = ipos + ksize
        end if
        if ( ma%left /= mpi_proc_null ) then
          ib = ipos
          do k = k1, k2
            ml(jde1-1,k) = r8vector2(ib)
            ib = ib + 1
          end do
          ipos = ipos + ksize
        end if
      end if
    end if
  end subroutine real8_bdy_exchange_left_right

  subroutine real4_bdy_exchange_left_right(ml,k1,k2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: k1, k2
    integer(ik4) :: ssize, ksize, k, ib, irc, ipos
    integer(ik4), dimension(4) :: req
    req = mpi_request_null
    ksize = k2-k1+1
    irc = 0
    ipos = 1
    ssize = 2*ksize
    if ( size(r4vector1) < ssize ) then
      call getmem1d(r4vector1,1,ssize,'real4_bdy_exchange_left_right')
    end if
    if ( size(r4vector2) < ssize ) then
      call getmem1d(r4vector2,1,ssize,'real4_bdy_exchange_left_right')
    end if
    if ( ma%bandflag ) then
      ib = 1
      do k = k1, k2
        r4vector1(ib) = ml(jde2,k)
        ib = ib + 1
      end do
      call cyclic_exchange_array(r4vector1,r4vector2, &
                                 ksize,ma%left,ma%right,__LINE__)
      ib = 1
      do k = k1, k2
        ml(jde1-1,k) = r4vector2(ib)
        ib = ib + 1
      end do
      ib = 1
      do k = k1, k2
        r4vector1(ib) = ml(jde1,k)
        ib = ib + 1
      end do
      call cyclic_exchange_array(r4vector1,r4vector2, &
                                 ksize,ma%right,ma%left,__LINE__)
      ib = 1
      do k = k1, k2
        ml(jde2+1,k) = r4vector2(ib)
        ib = ib + 1
      end do
    else
      if ( ma%right /= mpi_proc_null) then
        ib = ipos
        do k = k1, k2
          r4vector1(ib) = ml(jde2,k)
          ib = ib + 1
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ksize), &
                            r4vector2(ipos:ipos+ksize),ksize, &
                            ma%right,tag_rl,tag_lr, &
                            req(irc),req(irc+2))
        ipos = ipos + ksize
      end if
      if ( ma%left /= mpi_proc_null ) then
        ib = ipos
        do k = k1, k2
          r4vector1(ib) = ml(jde1,k)
          ib = ib + 1
        end do
        irc = irc + 1
        call exchange_array(r4vector1(ipos:ipos+ksize), &
                            r4vector2(ipos:ipos+ksize),ksize, &
                            ma%left,tag_lr,tag_rl, &
                            req(irc),req(irc+2))
        ipos = ipos + ksize
      end if
      if ( irc /= 0 ) then
        call mpi_waitall(4,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
        if ( mpierr /= mpi_success ) then
          call fatal(__FILE__,__LINE__,'mpi_waitall error.')
        end if
#endif
        ipos = 1
        if ( ma%right /= mpi_proc_null) then
          ib = ipos
          do k = k1, k2
            ml(jde2+1,k) = r4vector2(ib)
            ib = ib + 1
          end do
          ipos = ipos + ksize
        end if
        if ( ma%left /= mpi_proc_null ) then
          ib = ipos
          do k = k1, k2
            ml(jde1-1,k) = r4vector2(ib)
            ib = ib + 1
          end do
          ipos = ipos + ksize
        end if
      end if
    end if
  end subroutine real4_bdy_exchange_left_right

  subroutine real8_bdy_exchange_bottom_top(ml,k1,k2)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: k1, k2
    integer(ik4) :: ssize, ksize, k, ib, irc, ipos
    integer(ik4), dimension(4) :: req
    req = mpi_request_null
    ksize = k2-k1+1
    ssize = 2*ksize
    if ( size(r8vector1) < ssize ) then
      call getmem1d(r8vector1,1,ssize,'real8_bdy_exchange_bottom_top')
    end if
    if ( size(r8vector2) < ssize ) then
      call getmem1d(r8vector2,1,ssize,'real8_bdy_exchange_bottom_top')
    end if
    irc = 0
    ipos = 1
    if ( ma%top /= mpi_proc_null) then
      ib = ipos
      do k = k1, k2
        r8vector1(ib) = ml(ide2,k)
        ib = ib + 1
      end do
      irc = irc + 1
      call exchange_array(r8vector1(ipos:ipos+ksize), &
                          r8vector2(ipos:ipos+ksize),ksize, &
                          ma%top,tag_tb,tag_bt, &
                          req(irc),req(irc+2))
      ipos = ipos + ksize
    end if
    if ( ma%bottom /= mpi_proc_null ) then
      ib = ipos
      do k = k1, k2
        r8vector1(ib) = ml(ide1,k)
        ib = ib + 1
      end do
      irc = irc + 1
      call exchange_array(r8vector1(ipos:ipos+ksize), &
                          r8vector2(ipos:ipos+ksize),ksize, &
                          ma%bottom,tag_bt,tag_tb, &
                          req(irc),req(irc+2))
      ipos = ipos + ksize
    end if
    if ( irc /= 0 ) then
      call mpi_waitall(4,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      ipos = 1
      if ( ma%top /= mpi_proc_null) then
        ib = ipos
        do k = k1, k2
          ml(ide2+1,k) = r8vector2(ib)
          ib = ib + 1
        end do
        ipos = ipos + ksize
      end if
      if ( ma%bottom /= mpi_proc_null ) then
        ib = ipos
        do k = k1, k2
          ml(ide1-1,k) = r8vector2(ib)
          ib = ib + 1
        end do
        ipos = ipos + ksize
      end if
    end if
  end subroutine real8_bdy_exchange_bottom_top

  subroutine real4_bdy_exchange_bottom_top(ml,k1,k2)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:), intent(inout) :: ml
    integer(ik4), intent(in) :: k1, k2
    integer(ik4) :: ssize, ksize, k, ib, irc, ipos
    integer(ik4), dimension(4) :: req
    req = mpi_request_null
    ksize = k2-k1+1
    ssize = 2*ksize
    if ( size(r4vector1) < ssize ) then
      call getmem1d(r4vector1,1,ssize,'real4_bdy_exchange_bottom_top')
    end if
    if ( size(r4vector2) < ssize ) then
      call getmem1d(r4vector2,1,ssize,'real4_bdy_exchange_bottom_top')
    end if
    irc = 0
    ipos = 1
    if ( ma%top /= mpi_proc_null) then
      ib = ipos
      do k = k1, k2
        r4vector1(ib) = ml(ide2,k)
        ib = ib + 1
      end do
      irc = irc + 1
      call exchange_array(r4vector1(ipos:ipos+ksize), &
                          r4vector2(ipos:ipos+ksize),ksize, &
                          ma%top,tag_tb,tag_bt, &
                          req(irc),req(irc+2))
      ipos = ipos + ksize
    end if
    if ( ma%bottom /= mpi_proc_null ) then
      ib = ipos
      do k = k1, k2
        r4vector1(ib) = ml(ide1,k)
        ib = ib + 1
      end do
      irc = irc + 1
      call exchange_array(r4vector1(ipos:ipos+ksize), &
                          r4vector2(ipos:ipos+ksize),ksize, &
                          ma%bottom,tag_bt,tag_tb, &
                          req(irc),req(irc+2))
      ipos = ipos + ksize
    end if
    if ( irc /= 0 ) then
      call mpi_waitall(4,req,mpi_statuses_ignore,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_waitall error.')
      end if
#endif
      ipos = 1
      if ( ma%top /= mpi_proc_null) then
        ib = ipos
        do k = k1, k2
          ml(ide2+1,k) = r4vector2(ib)
          ib = ib + 1
        end do
        ipos = ipos + ksize
      end if
      if ( ma%bottom /= mpi_proc_null ) then
        ib = ipos
        do k = k1, k2
          ml(ide1-1,k) = r4vector2(ib)
          ib = ib + 1
        end do
        ipos = ipos + ksize
      end if
    end if
  end subroutine real4_bdy_exchange_bottom_top

  subroutine real8_2d_grid_fill_extend1(a,b)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:), intent(in) :: a
    real(rk8), pointer, contiguous, dimension(:,:), intent(inout) :: b
    call grid_collect(a,b,max(jce1,lbound(a,1)),min(jce2,ubound(a,1)), &
                          max(ice1,lbound(a,2)),min(ice2,ubound(a,2)))
    ! Extend on a 'fake' filled cross grid (+1)
    b(1,:) = b(2,:)
    b(jx,:) = b(jx-2,:)
    b(jx-1,:) = b(jx-2,:)
    b(:,1) = b(:,2)
    b(:,iy) = b(:,iy-2)
    b(:,iy-1) = b(:,iy-2)
    b(1,1) = b(2,2)
    b(jx,1) = b(jx-1,2)
    b(jx-1,1) = b(jx-1,2)
    b(1,iy) = b(2,iy-2)
    b(1,iy-1) = b(2,iy-2)
    b(jx,iy) = b(jx-2,iy-2)
    b(jx-1,iy-1) = b(jx-2,iy-2)
    call mpi_bcast(b,iy*jx,mpi_real8,iocpu,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
#endif
  end subroutine real8_2d_grid_fill_extend1

  subroutine real4_2d_grid_fill_extend1(a,b)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:), intent(in) :: a
    real(rk4), pointer, contiguous, dimension(:,:), intent(inout) :: b
    call grid_collect(a,b,max(jce1,lbound(a,1)),min(jce2,ubound(a,1)), &
                          max(ice1,lbound(a,2)),min(ice2,ubound(a,2)))
    ! Extend on a 'fake' filled cross grid (+1)
    b(1,:) = b(2,:)
    b(jx,:) = b(jx-2,:)
    b(jx-1,:) = b(jx-2,:)
    b(:,1) = b(:,2)
    b(:,iy) = b(:,iy-2)
    b(:,iy-1) = b(:,iy-2)
    b(1,1) = b(2,2)
    b(jx,1) = b(jx-1,2)
    b(jx-1,1) = b(jx-1,2)
    b(1,iy) = b(2,iy-2)
    b(1,iy-1) = b(2,iy-2)
    b(jx,iy) = b(jx-2,iy-2)
    b(jx-1,iy-1) = b(jx-2,iy-2)
    call mpi_bcast(b,iy*jx,mpi_real4,iocpu,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
#endif
  end subroutine real4_2d_grid_fill_extend1

  subroutine real8_2d_grid_fill_extend2(a,b,i1,i2,j1,j2)
    implicit none
    integer(ik4), intent(in) :: i1, i2, j1, j2
    real(rk8), pointer, contiguous, dimension(:,:), intent(in) :: a
    real(rk8), pointer, contiguous, dimension(:,:), intent(inout) :: b
    call grid_collect(a,b,i1,i2,j1,j2)
    call mpi_bcast(b,product(shape(b)),mpi_real8,iocpu,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
#endif
  end subroutine real8_2d_grid_fill_extend2

  subroutine real4_2d_grid_fill_extend2(a,b,i1,i2,j1,j2)
    implicit none
    integer(ik4), intent(in) :: i1, i2, j1, j2
    real(rk4), pointer, contiguous, dimension(:,:), intent(in) :: a
    real(rk4), pointer, contiguous, dimension(:,:), intent(inout) :: b
    call grid_collect(a,b,i1,i2,j1,j2)
    call mpi_bcast(b,product(shape(b)),mpi_real4,iocpu,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
#endif
  end subroutine real4_2d_grid_fill_extend2

  subroutine uvtentotenx(u,v,ux,vx)
    implicit none
    real(rkx), intent(inout), dimension(:,:,:), pointer, contiguous :: u, v
    real(rkx), intent(inout), dimension(:,:,:), pointer, contiguous :: ux, vx
    integer(ik4) :: i, j, k

    call exchange_lr(u,2,jdi1,jdi2,ici1,ici2,1,kz)
    call exchange_bt(v,2,jci1,jci2,idi1,idi2,1,kz)
    ! Back to wind points
    do concurrent ( j = jcii1:jcii2, i = ici1:ici2, k = 1:kz )
      ux(j,i,k) = 0.5625_rkx * (u(j+1,i,k)+u(j,i,k)) - &
                  0.0625_rkx * (u(j+2,i,k)+u(j-1,i,k))
    end do
    if ( ma%has_bdyleft ) then
      do concurrent ( i = ici1:ici2, k = 1:kz )
        ux(jci1,i,k) = 0.5_rkx * (u(jdi1,i,k)+u(jdii1,i,k))
      end do
    end if
    if ( ma%has_bdyright ) then
      do concurrent ( i = ici1:ici2, k = 1:kz )
        ux(jci2,i,k) = 0.5_rkx*(u(jdi2,i,k) + u(jdii2,i,k))
      end do
    end if
    do concurrent ( j = jci1:jci2, i = icii1:icii2, k = 1:kz )
      vx(j,i,k) = 0.5625_rkx * (v(j,i+1,k)+v(j,i,k)) - &
                  0.0625_rkx * (v(j,i+2,k)+v(j,i-1,k))
    end do
    if ( ma%has_bdybottom ) then
      do concurrent( j = jci1:jci2, k = 1:kz )
        vx(j,ici1,k) = 0.5_rkx * (v(j,idi1,k)+v(j,idii1,k))
      end do
    end if
    if ( ma%has_bdytop ) then
      do concurrent( j = jci1:jci2, k = 1:kz )
        vx(j,ici2,k) = 0.5_rkx*(v(j,idi2,k) + v(j,idii2,k))
      end do
    end if
  end subroutine uvtentotenx

  subroutine tenxtouvten(ux,vx,u,v)
    implicit none
    real(rkx), intent(inout), dimension(:,:,:), pointer, contiguous :: ux, vx
    real(rkx), intent(inout), dimension(:,:,:), pointer, contiguous :: u, v
    integer(ik4) :: i, j, k

    call exchange_lr(ux,2,jce1,jce2,ice1,ice2,1,kz)
    call exchange_bt(vx,2,jce1,jce2,ice1,ice2,1,kz)

    ! Back to wind points: U (fourth order)
    do concurrent ( j = jdii1:jdii2, i = ici1:ici2, k = 1:kz )
      u(j,i,k) = 0.5625_rkx * (ux(j,i,k)  +ux(j-1,i,k)) - &
                 0.0625_rkx * (ux(j+1,i,k)+ux(j-2,i,k))
    end do
    ! Back to wind points: V (fourth order)
    do concurrent ( j = jci1:jci2, i = idii1:idii2, k = 1:kz )
      v(j,i,k) = 0.5625_rkx * (vx(j,i,k)  +vx(j,i-1,k)) - &
                 0.0625_rkx * (vx(j,i+1,k)+vx(j,i-2,k))
    end do
  end subroutine tenxtouvten
  !
  ! Takes u and v tendencies on the cross grid (as t, qv, qc, etc.)
  ! and interpolates the u and v to the dot grid.
  ! This routine sheilds the user of the function from the need to worry
  ! about the details of the domain decomposition.
  !
  ! Written by Travis A. O'Brien 01/04/11.
  !
  subroutine uvcross2dot(ux,vx,ud,vd)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: ux, vx
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: ud, vd
    integer(ik4) :: i, j, k

    call exchange(ux,1,jce1,jce2,ice1,ice2,1,kz)
    call exchange(vx,1,jce1,jce2,ice1,ice2,1,kz)

    !
    !  o     o     o     o     o     o     o
    !
    !     x     x     x     x     x     x
    !
    !  o     o     o     o     o     o     o
    !         (i-1,j-1)     (i,j-1)
    !     x     x     x-----x     x     x
    !                 |(i,j)|
    !  o     o     o  |  o  |  o     o     o
    !                 |     |
    !     x     x     x-----x     x     x
    !           (i-1,j)     (i,j)
    !  o     o     o     o     o     o     o
    !
    !     x     x     x     x     x     x
    !
    !  o     o     o     o     o     o     o

    ! Perform the bilinear interpolation necessary
    ! to put the u and v variables on the dot grid.

    do concurrent ( j = jdi1:jdi2, i = idi1:idi2, k = 1:kz )
      ud(j,i,k) =  d_rfour*(ux(j,i,k) + ux(j-1,i,k) +   &
                            ux(j,i-1,k) + ux(j-1,i-1,k))
      vd(j,i,k) =  d_rfour*(vx(j,i,k) + vx(j-1,i,k) +   &
                            vx(j,i-1,k) + vx(j-1,i-1,k))
    end do
  end subroutine uvcross2dot

  subroutine uvdot2cross(ud,vd,ux,vx)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: ud, vd
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: ux, vx
    integer(ik4) :: i, j, k

    call exchange(ud,1,jdi1,jdi2,idi1,idi2,1,kz)
    call exchange(vd,1,jdi1,jdi2,idi1,idi2,1,kz)

    !
    !     o     o     o     o     o     o
    !
    !        x     x     x     x     x
    !           (i+1,j)   (i+1,j+1)
    !     o     o     o-----o     o     o
    !                 |(i,j)|
    !        x     x  |  x  |  x     x
    !                 |     |
    !     o     o     o-----o     o     o
    !             (i,j)     (i,j+1)
    !        x     x     x     x     x
    !
    !     o     o     o     o     o     o
    !

    ! Perform the bilinear interpolation necessary
    ! to put the u and v variables on the cross grid.

    do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
      ux(j,i,k) =  d_rfour*(ud(j,i,  k) + ud(j+1,i,  k) + &
                            ud(j,i+1,k) + ud(j+1,i+1,k))
      vx(j,i,k) =  d_rfour*(vd(j,i  ,k) + vd(j+1,i  ,k) + &
                            vd(j,i+1,k) + vd(j+1,i+1,k))
    end do
  end subroutine uvdot2cross

  subroutine cross2dot2d(x,d)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: x
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: d
    integer(ik4) :: i, j

    call exchange(x,1,jci1,jci2,ici1,ici2)

    do concurrent ( j = jdii1:jdii2, i = idii1:idii2 )
      d(j,i) =  d_rfour*(x(j,i)   + x(j-1,i)   + &
                           x(j,i-1) + x(j-1,i-1))
    end do
    if ( ma%has_bdyleft ) then
      do i = idii1, idii2
        d(jdi1,i) = d_half*(x(jci1,i) + x(jci1,i-1))
      end do
    end if
    if ( ma%has_bdyright ) then
      do i = idii1, idii2
        d(jdi2,i) = d_half*(x(jci2,i) + x(jci2,i-1))
      end do
    end if
    if ( ma%has_bdytop ) then
      do j = jdii1, jdii2
        d(j,idi2) = d_half*(x(j,ici2) + x(j-1,ici2))
      end do
    end if
    if ( ma%has_bdybottom ) then
      do j = jdii1, jdii2
        d(j,idi1) = d_half*(x(j,ici1) + x(j-1,ici1))
      end do
    end if
    if ( ma%has_bdytopleft ) then
      d(jdi1,idi2) = x(jci1,ici2)
    end if
    if ( ma%has_bdybottomleft ) then
      d(jdi1,idi1) = x(jci1,ici1)
    end if
    if ( ma%has_bdytopright ) then
      d(jdi2,idi2) = x(jci2,ici2)
    end if
    if ( ma%has_bdybottomright ) then
      d(jdi2,idi1) = x(jci2,ici1)
    end if
  end subroutine cross2dot2d

  subroutine cross2dot3d(x,d)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: x
    real(rkx), pointer, contiguous, dimension(:,:,:), intent(inout) :: d
    integer(ik4) :: i, j, k

    call exchange(x,1,jci1,jci2,ici1,ici2,1,kz)

    do concurrent ( j = jdii1:jdii2, i = idii1:idii2, k = 1:kz )
      d(j,i,k) =  d_rfour*(x(j,i,k)   + x(j-1,i,k)   + &
                           x(j,i-1,k) + x(j-1,i-1,k))
    end do
    if ( ma%has_bdyleft ) then
      do k = 1, kz
        do i = idii1, idii2
          d(jdi1,i,k) = d_half*(x(jci1,i,k) + x(jci1,i-1,k))
        end do
      end do
    end if
    if ( ma%has_bdyright ) then
      do k = 1, kz
        do i = idii1, idii2
          d(jdi2,i,k) = d_half*(x(jci2,i,k) + x(jci2,i-1,k))
        end do
      end do
    end if
    if ( ma%has_bdytop ) then
      do k = 1, kz
        do j = jdii1, jdii2
          d(j,idi2,k) = d_half*(x(j,ici2,k) + x(j-1,ici2,k))
        end do
      end do
    end if
    if ( ma%has_bdybottom ) then
      do k = 1, kz
        do j = jdii1, jdii2
          d(j,idi1,k) = d_half*(x(j,ici1,k) + x(j-1,ici1,k))
        end do
      end do
    end if
    if ( ma%has_bdytopleft ) then
      do k = 1, kz
        d(jdi1,idi2,k) = x(jci1,ici2,k)
      end do
    end if
    if ( ma%has_bdybottomleft ) then
      do k = 1, kz
        d(jdi1,idi1,k) = x(jci1,ici1,k)
      end do
    end if
    if ( ma%has_bdytopright ) then
      do k = 1, kz
        d(jdi2,idi2,k) = x(jci2,ici2,k)
      end do
    end if
    if ( ma%has_bdybottomright ) then
      do k = 1, kz
        d(jdi2,idi1,k) = x(jci2,ici1,k)
      end do
    end if
  end subroutine cross2dot3d

  subroutine psc2psd(pc,pd)
    implicit none
    real(rkx), pointer, contiguous, dimension(:,:), intent(in)  :: pc
    real(rkx), pointer, contiguous, dimension(:,:), intent(inout) :: pd
    integer(ik4) :: i, j
    !
    ! Internal points
    !
    do concurrent ( j = jdi1:jdi2, i = idi1:idi2 )
      pd(j,i) = (pc(j,i)+pc(j,i-1)+pc(j-1,i)+pc(j-1,i-1))*d_rfour
    end do
    !
    ! Boundaries
    !
    if ( ma%has_bdytop ) then
      do j = jdi1, jdi2
        pd(j,ide2) = (pc(j,ice2)+pc(j-1,ice2))*d_half
      end do
    end if
    if ( ma%has_bdybottom ) then
      do j = jdi1, jdi2
        pd(j,ide1)  = (pc(j,ice1)+pc(j-1,ice1))*d_half
      end do
    end if
    if ( ma%has_bdyleft ) then
      do i = idi1, idi2
        pd(jde1,i) = (pc(jce1,i)+pc(jce1,i-1))*d_half
      end do
    end if
    if ( ma%has_bdyright ) then
      do i = idi1, idi2
        pd(jde2,i) = (pc(jce2,i)+pc(jce2,i-1))*d_half
      end do
    end if
    !
    ! Corner points
    !
    if ( ma%has_bdybottomleft ) then
      pd(jde1,ide1) = pc(jce1,ice1)
    end if
    if ( ma%has_bdytopleft ) then
      pd(jde1,ide2) = pc(jce1,ice2)
    end if
    if ( ma%has_bdybottomright ) then
      pd(jde2,ide1) = pc(jce2,ice1)
    end if
    if ( ma%has_bdytopright ) then
      pd(jde2,ide2) = pc(jce2,ice2)
    end if
  end subroutine psc2psd

  subroutine grid_nc_create_var2d(varname,ldot,val,xvar)
    implicit none
    character(len=*), intent(in) :: varname
    logical, intent(in) :: ldot
    real(rk8), pointer, contiguous, dimension(:,:), intent(in) :: val
    type (grid_nc_var2d), intent(inout) :: xvar
    integer(ik4) :: istat
    integer(ik4), dimension(3) :: idims

    xvar%varname = varname
    call assignpnt(val,xvar%val)

    if ( ldot ) then
      xvar%nx = njdot
      xvar%ny = nidot
      xvar%mynx1 = jde1
      if ( lbound(val,1) > jde1 ) xvar%mynx1 = lbound(val,1)
      xvar%mynx2 = jde2
      if ( ubound(val,1) < jde2 ) xvar%mynx2 = ubound(val,1)
      xvar%myny1 = ide1
      if ( lbound(val,2) > ide1 ) xvar%myny1 = lbound(val,2)
      xvar%myny2 = ide2
      if ( ubound(val,2) < ide2 ) xvar%myny2 = ubound(val,2)
    else
      xvar%nx = njcross
      xvar%ny = nicross
      xvar%mynx1 = jce1
      if ( lbound(val,1) > jce1 ) xvar%mynx1 = lbound(val,1)
      xvar%mynx2 = jce2
      if ( ubound(val,1) < jce2 ) xvar%mynx2 = ubound(val,1)
      xvar%myny1 = ice1
      if ( lbound(val,2) > ice1 ) xvar%myny1 = lbound(val,2)
      xvar%myny2 = ice2
      if ( ubound(val,2) < ice2 ) xvar%myny2 = ubound(val,2)
    end if
    xvar%irec = 1
    if ( myid /= iocpu ) return
    istat = nf90_create(trim(xvar%varname)//'.nc',nf90_clobber,xvar%ncid)
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_dim(xvar%ncid, 'JX', xvar%nx, idims(1))
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_dim(xvar%ncid, 'IY', xvar%ny, idims(2))
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_dim(xvar%ncid, 'KTAU', nf90_unlimited, idims(3))
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_var(xvar%ncid,xvar%varname,nf90_double,idims,xvar%varid)
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_enddef(xvar%ncid)
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    call getmem2d(xvar%iobuf,1,xvar%nx,1,xvar%ny,'var2d:iobuf')
  end subroutine grid_nc_create_var2d

  subroutine grid_nc_write_var2d(xvar)
    implicit none
    type (grid_nc_var2d), intent(inout) :: xvar
    integer(ik4) :: istat
    integer(ik4), dimension(3) :: istart, icount
    if ( .not. associated(xvar%val) .or. xvar%irec < 1 ) then
      return
    end if
    call grid_collect(xvar%val,xvar%iobuf, &
                      xvar%mynx1,xvar%mynx2,xvar%myny1,xvar%myny2)
    if ( myid == iocpu ) then
      istart(3) = xvar%irec
      istart(2) = 1
      istart(1) = 1
      icount(3) = 1
      icount(2) = xvar%ny
      icount(1) = xvar%nx
      istat = nf90_put_var(xvar%ncid,xvar%varid,xvar%iobuf,istart,icount)
      if ( istat /= nf90_noerr ) then
        write(stderr, *) nf90_strerror(istat)
        return
      end if
      istat = nf90_sync(xvar%ncid)
      if ( istat /= nf90_noerr ) then
        write(stderr, *) nf90_strerror(istat)
        return
      end if
    end if
    xvar%irec = xvar%irec + 1
  end subroutine grid_nc_write_var2d

  subroutine grid_nc_destroy_var2d(xvar)
    implicit none
    type (grid_nc_var2d), intent(inout) :: xvar
    integer(ik4) :: istat
    if ( myid == iocpu ) then
      istat = nf90_close(xvar%ncid)
      if ( istat /= nf90_noerr ) then
        write(stderr, *) nf90_strerror(istat)
        return
      end if
    end if
    xvar%ncid = -1
    xvar%irec = -1
    nullify(xvar%val)
    call relmem2d(xvar%iobuf)
  end subroutine grid_nc_destroy_var2d

  subroutine grid_nc_create_var3d(varname,ldot,val,xvar)
    implicit none
    character(len=*), intent(in) :: varname
    logical, intent(in) :: ldot
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(in) :: val
    type (grid_nc_var3d), intent(inout) :: xvar
    integer(ik4) :: istat
    integer(ik4), dimension(4) :: idims

    xvar%varname = varname
    call assignpnt(val,xvar%val)

    if ( ldot ) then
      xvar%nx = njdot
      xvar%ny = nidot
      xvar%mynx1 = jde1
      if ( lbound(val,1) > jde1 ) xvar%mynx1 = lbound(val,1)
      xvar%mynx2 = jde2
      if ( ubound(val,1) < jde2 ) xvar%mynx2 = ubound(val,1)
      xvar%myny1 = ide1
      if ( lbound(val,2) > ide1 ) xvar%myny1 = lbound(val,2)
      xvar%myny2 = ide2
      if ( ubound(val,2) < ide2 ) xvar%myny2 = ubound(val,2)
    else
      xvar%nx = njcross
      xvar%ny = nicross
      xvar%mynx1 = jce1
      if ( lbound(val,1) > jce1 ) xvar%mynx1 = lbound(val,1)
      xvar%mynx2 = jce2
      if ( ubound(val,1) < jce2 ) xvar%mynx2 = ubound(val,1)
      xvar%myny1 = ice1
      if ( lbound(val,2) > ice1 ) xvar%myny1 = lbound(val,2)
      xvar%myny2 = ice2
      if ( ubound(val,2) < ice2 ) xvar%myny2 = ubound(val,2)
    end if
    xvar%nz = size(xvar%val,3)
    xvar%irec = 1
    if ( myid /= iocpu ) return
    istat = nf90_create(trim(xvar%varname)//'.nc',nf90_clobber,xvar%ncid)
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_dim(xvar%ncid, 'JX', xvar%nx, idims(1))
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_dim(xvar%ncid, 'IY', xvar%ny, idims(2))
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_dim(xvar%ncid, 'KZ', xvar%nz, idims(3))
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_dim(xvar%ncid, 'KTAU', nf90_unlimited, idims(4))
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_var(xvar%ncid,xvar%varname,nf90_double,idims,xvar%varid)
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_enddef(xvar%ncid)
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    call getmem3d(xvar%iobuf,1,xvar%nx,1,xvar%ny,1,xvar%nz,'var3d:iobuf')
  end subroutine grid_nc_create_var3d

  subroutine grid_nc_write_var3d(xvar)
    implicit none
    type (grid_nc_var3d), intent(inout) :: xvar
    integer(ik4) :: istat
    integer(ik4), dimension(4) :: istart, icount
    if ( .not. associated(xvar%val) .or. xvar%irec < 1 ) then
      return
    end if
    call grid_collect(xvar%val,xvar%iobuf, &
                      xvar%mynx1,xvar%mynx2,xvar%myny1,xvar%myny2,1,xvar%nz)
    if ( myid == iocpu ) then
      istart(4) = xvar%irec
      istart(3) = 1
      istart(2) = 1
      istart(1) = 1
      icount(4) = 1
      icount(3) = xvar%nz
      icount(2) = xvar%ny
      icount(1) = xvar%nx
      istat = nf90_put_var(xvar%ncid,xvar%varid,xvar%iobuf,istart,icount)
      if ( istat /= nf90_noerr ) then
        write(stderr, *) nf90_strerror(istat)
        return
      end if
      istat = nf90_sync(xvar%ncid)
      if ( istat /= nf90_noerr ) then
        write(stderr, *) nf90_strerror(istat)
        return
      end if
    end if
    xvar%irec = xvar%irec + 1
  end subroutine grid_nc_write_var3d

  subroutine grid_nc_destroy_var3d(xvar)
    implicit none
    type (grid_nc_var3d), intent(inout) :: xvar
    integer(ik4) :: istat
    if ( myid == iocpu ) then
      istat = nf90_close(xvar%ncid)
      if ( istat /= nf90_noerr ) then
        write(stderr, *) nf90_strerror(istat)
        return
      end if
    end if
    xvar%ncid = -1
    xvar%irec = -1
    nullify(xvar%val)
    call relmem3d(xvar%iobuf)
  end subroutine grid_nc_destroy_var3d

  subroutine grid_nc_create_var4d(varname,ldot,val,xvar)
    implicit none
    character(len=*), intent(in) :: varname
    logical, intent(in) :: ldot
    real(rk8), pointer, contiguous, dimension(:,:,:,:), intent(in) :: val
    type (grid_nc_var4d), intent(inout) :: xvar
    integer(ik4) :: istat
    integer(ik4), dimension(5) :: idims

    xvar%varname = varname
    call assignpnt(val,xvar%val)

    if ( ldot ) then
      xvar%nx = njdot
      xvar%ny = nidot
      xvar%mynx1 = jde1
      if ( lbound(val,1) > jde1 ) xvar%mynx1 = lbound(val,1)
      xvar%mynx2 = jde2
      if ( ubound(val,1) < jde2 ) xvar%mynx2 = ubound(val,1)
      xvar%myny1 = ide1
      if ( lbound(val,2) > ide1 ) xvar%myny1 = lbound(val,2)
      xvar%myny2 = ide2
      if ( ubound(val,2) < ide2 ) xvar%myny2 = ubound(val,2)
    else
      xvar%nx = njcross
      xvar%ny = nicross
      xvar%mynx1 = jce1
      if ( lbound(val,1) > jce1 ) xvar%mynx1 = lbound(val,1)
      xvar%mynx2 = jce2
      if ( ubound(val,1) < jce2 ) xvar%mynx2 = ubound(val,1)
      xvar%myny1 = ice1
      if ( lbound(val,2) > ice1 ) xvar%myny1 = lbound(val,2)
      xvar%myny2 = ice2
      if ( ubound(val,2) < ice2 ) xvar%myny2 = ubound(val,2)
    end if
    xvar%nz = size(xvar%val,3)
    xvar%nl = size(xvar%val,4)
    xvar%irec = 1
    if ( myid /= iocpu ) return
    istat = nf90_create(trim(xvar%varname)//'.nc',nf90_clobber,xvar%ncid)
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_dim(xvar%ncid, 'JX', xvar%nx, idims(1))
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_dim(xvar%ncid, 'IY', xvar%ny, idims(2))
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_dim(xvar%ncid, 'KZ', xvar%nz, idims(3))
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_dim(xvar%ncid, 'LL', xvar%nl, idims(4))
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_dim(xvar%ncid, 'KTAU', nf90_unlimited, idims(5))
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_def_var(xvar%ncid,xvar%varname,nf90_double,idims,xvar%varid)
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    istat = nf90_enddef(xvar%ncid)
    if ( istat /= nf90_noerr ) then
      write(stderr, *) nf90_strerror(istat)
      return
    end if
    call getmem4d(xvar%iobuf,1,xvar%nx,1,xvar%ny,1,xvar%nz, &
                  1,xvar%nl,'var3d:iobuf')
  end subroutine grid_nc_create_var4d

  subroutine grid_nc_write_var4d(xvar)
    implicit none
    type (grid_nc_var4d), intent(inout) :: xvar
    integer(ik4) :: istat
    integer(ik4), dimension(5) :: istart, icount
    if ( .not. associated(xvar%val) .or. xvar%irec < 1 ) then
      return
    end if
    call grid_collect(xvar%val,xvar%iobuf, &
                      xvar%mynx1,xvar%mynx2,xvar%myny1,xvar%myny2, &
                      1,xvar%nz,1,xvar%nl)
    if ( myid == iocpu ) then
      istart(5) = xvar%irec
      istart(4) = 1
      istart(3) = 1
      istart(2) = 1
      istart(1) = 1
      icount(5) = 1
      icount(4) = xvar%nl
      icount(3) = xvar%nz
      icount(2) = xvar%ny
      icount(1) = xvar%nx
      istat = nf90_put_var(xvar%ncid,xvar%varid,xvar%iobuf,istart,icount)
      if ( istat /= nf90_noerr ) then
        write(stderr, *) nf90_strerror(istat)
        return
      end if
      istat = nf90_sync(xvar%ncid)
      if ( istat /= nf90_noerr ) then
        write(stderr, *) nf90_strerror(istat)
        return
      end if
    end if
    xvar%irec = xvar%irec + 1
  end subroutine grid_nc_write_var4d

  subroutine grid_nc_destroy_var4d(xvar)
    implicit none
    type (grid_nc_var4d), intent(inout) :: xvar
    integer(ik4) :: istat
    if ( myid == iocpu ) then
      istat = nf90_close(xvar%ncid)
      if ( istat /= nf90_noerr ) then
        write(stderr, *) nf90_strerror(istat)
        return
      end if
    end if
    xvar%ncid = -1
    xvar%irec = -1
    nullify(xvar%val)
    call relmem4d(xvar%iobuf)
  end subroutine grid_nc_destroy_var4d

  subroutine gather_r(f_collect,f_sub)
    implicit none
    real(rk8), dimension(:), intent(out) :: f_collect
    real(rk8), intent(in) :: f_sub
    real(rk8), dimension(1) :: tmp
    tmp(1) = f_sub
    call mpi_gather(tmp,      1,mpi_real8, &
                    f_collect,1,mpi_real8,iocpu,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) THEN
      call fatal(__FILE__,__LINE__,'error in mpi_gather!!')
    end if
#endif
  end subroutine gather_r

  subroutine gather_i(i_collect,i_sub)
    implicit none
    integer(ik4), dimension(:), intent(out) :: i_collect
    integer(ik4), intent(in) :: i_sub
    integer(ik4), dimension(1) :: tmp
    tmp(1) = i_sub
    call mpi_gather(tmp,      1,mpi_integer4, &
                    i_collect,1,mpi_integer4,iocpu,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) THEN
      call fatal(__FILE__,__LINE__,'error in mpi_gather!!')
    end if
#endif
  end subroutine gather_i

  subroutine allgather_r(f_collect,f_sub)
    implicit none
    real(rk8), dimension(:), intent(out) :: f_collect
    real(rk8), intent(in) :: f_sub
    real(rk8), dimension(1) :: tmp
    tmp(1) = f_sub
    call mpi_allgather(tmp,      1,mpi_real8, &
                       f_collect,1,mpi_real8,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) THEN
      call fatal(__FILE__,__LINE__,'error in mpi_allgather!!')
    end if
#endif
  end subroutine allgather_r

  subroutine allgather_i(i_collect,i_sub)
    implicit none
    integer(ik4), dimension(:), intent(out) :: i_collect
    integer(ik4), intent(in) :: i_sub
    integer(ik4), dimension(1) :: tmp
    tmp(1) = i_sub
    call mpi_allgather(tmp,      1,mpi_integer4, &
                       i_collect,1,mpi_integer4,mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) THEN
      call fatal(__FILE__,__LINE__,'error in mpi_allgather!!')
    end if
#endif
  end subroutine allgather_i

  subroutine reorder_add_subgrid_2d_real8(var3,var2,mask)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(in) :: var3
    real(rk8), pointer, contiguous, dimension(:,:), intent(inout) :: var2
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(in), optional :: mask
    integer(ik4) :: i, j, ii, jj, n1, n2
    if ( present(mask) ) then
      do i = ici1, ici2
        do j = jci1, jci2
          do n2 = 1, nsg
            ii = (i-1) * nsg + n2
            do n1 = 1, nsg
              jj = (j-1) * nsg + n1
              if ( mask((n2-1)*nsg+n1,j,i) > 0 ) then
                var2(jj,ii) = var2(jj,ii) + var3((n2-1)*nsg+n1,j,i)
              else
                var2(jj,ii) = dmissval
              end if
            end do
          end do
        end do
      end do
    else
      do i = ici1, ici2
        do j = jci1, jci2
          do n2 = 1, nsg
            ii = (i-1) * nsg + n2
            do n1 = 1, nsg
              jj = (j-1) * nsg + n1
              var2(jj,ii) = var2(jj,ii) + var3((n2-1)*nsg+n1,j,i)
            end do
          end do
        end do
      end do
    end if
  end subroutine reorder_add_subgrid_2d_real8

  subroutine reorder_add_subgrid_2d_real4(var3,var2,mask)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(in) :: var3
    real(rk4), pointer, contiguous, dimension(:,:), intent(inout) :: var2
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(in), optional :: mask
    integer(ik4) :: i, j, ii, jj, n1, n2
    if ( present(mask) ) then
      do i = ici1, ici2
        do j = jci1, jci2
          do n2 = 1, nsg
            ii = (i-1) * nsg + n2
            do n1 = 1, nsg
              jj = (j-1) * nsg + n1
              if ( mask((n2-1)*nsg+n1,j,i) > 0 ) then
                var2(jj,ii) = var2(jj,ii) + var3((n2-1)*nsg+n1,j,i)
              else
                var2(jj,ii) = smissval
              end if
            end do
          end do
        end do
      end do
    else
      do i = ici1, ici2
        do j = jci1, jci2
          do n2 = 1, nsg
            ii = (i-1) * nsg + n2
            do n1 = 1, nsg
              jj = (j-1) * nsg + n1
              var2(jj,ii) = var2(jj,ii) + var3((n2-1)*nsg+n1,j,i)
            end do
          end do
        end do
      end do
    end if
  end subroutine reorder_add_subgrid_2d_real4

  subroutine reorder_subgrid_2d_real8(var3,var2,mask)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(in) :: var3
    real(rk8), pointer, contiguous, dimension(:,:), intent(inout) :: var2
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(in), optional :: mask
    integer(ik4) :: i, j, ii, jj, n1, n2
    if ( present(mask) ) then
      do i = ici1, ici2
        do j = jci1, jci2
          do n2 = 1, nsg
            ii = (i-1) * nsg + n2
            do n1 = 1, nsg
              jj = (j-1) * nsg + n1
              if ( mask((n2-1)*nsg+n1,j,i) > 0 ) then
                var2(jj,ii) = var3((n2-1)*nsg+n1,j,i)
              else
                var2(jj,ii) = dmissval
              end if
            end do
          end do
        end do
      end do
    else
      do i = ici1, ici2
        do j = jci1, jci2
          do n2 = 1, nsg
            ii = (i-1) * nsg + n2
            do n1 = 1, nsg
              jj = (j-1) * nsg + n1
              var2(jj,ii) = var3((n2-1)*nsg+n1,j,i)
            end do
          end do
        end do
      end do
    end if
  end subroutine reorder_subgrid_2d_real8

  subroutine reorder_subgrid_2d_real4(var3,var2,mask)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(in) :: var3
    real(rk4), pointer, contiguous, dimension(:,:), intent(inout) :: var2
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(in), optional :: mask
    integer(ik4) :: i, j, ii, jj, n1, n2
    if ( present(mask) ) then
      do i = ici1, ici2
        do j = jci1, jci2
          do n2 = 1, nsg
            ii = (i-1) * nsg + n2
            do n1 = 1, nsg
              jj = (j-1) * nsg + n1
              if ( mask((n2-1)*nsg+n1,j,i) > 0 ) then
                var2(jj,ii) = var3((n2-1)*nsg+n1,j,i)
              else
                var2(jj,ii) = smissval
              end if
            end do
          end do
        end do
      end do
    else
      do i = ici1, ici2
        do j = jci1, jci2
          do n2 = 1, nsg
            ii = (i-1) * nsg + n2
            do n1 = 1, nsg
              jj = (j-1) * nsg + n1
              var2(jj,ii) = var3((n2-1)*nsg+n1,j,i)
            end do
          end do
        end do
      end do
    end if
  end subroutine reorder_subgrid_2d_real4

  subroutine reorder_logical_global_subgrid_2d(var3,var2)
    implicit none
    logical, pointer, contiguous, dimension(:,:,:), intent(in) :: var3
    logical, pointer, contiguous, dimension(:,:), intent(inout) :: var2
    integer(ik4) :: i, j, ii, jj, n1, n2
    do i = iout1, iout2
      do j = jout1, jout2
        do n2 = 1, nsg
          ii = (i-1) * nsg + n2
          do n1 = 1, nsg
            jj = (j-1) * nsg + n1
            var2(jj,ii) = var3((n2-1)*nsg+n1,j,i)
          end do
        end do
      end do
    end do
  end subroutine reorder_logical_global_subgrid_2d

  subroutine reorder_subgrid_2d_logical(var3,var2)
    implicit none
    logical, pointer, contiguous, dimension(:,:,:), intent(in) :: var3
    logical, pointer, contiguous, dimension(:,:), intent(inout) :: var2
    integer(ik4) :: i, j, ii, jj, n1, n2
    do i = ici1, ici2
      do j = jci1, jci2
        do n2 = 1, nsg
          ii = (i-1) * nsg + n2
          do n1 = 1, nsg
            jj = (j-1) * nsg + n1
            var2(jj,ii) = var3((n2-1)*nsg+n1,j,i)
          end do
        end do
      end do
    end do
  end subroutine reorder_subgrid_2d_logical

  subroutine reorder_add_subgrid_2d3d_real8(var3,var2_3,l,mask)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(in) :: var3
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(inout) :: var2_3
    integer(ik4), optional, intent(in) :: l
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(in), optional :: mask
    integer(ik4) :: i, j, ii, jj, n1, n2, ll
    ll = 1
    if ( present(l) ) ll = l
    if ( present(mask) ) then
      do i = ici1, ici2
        do j = jci1, jci2
          do n2 = 1, nsg
            ii = (i-1) * nsg + n2
            do n1 = 1, nsg
              jj = (j-1) * nsg + n1
              if ( mask((n2-1)*nsg+n1,j,i) > 0 ) then
                var2_3(jj,ii,ll) = var2_3(jj,ii,ll) + var3((n2-1)*nsg+n1,j,i)
              else
                var2_3(jj,ii,ll) = dmissval
              end if
            end do
          end do
        end do
      end do
    else
      do i = ici1, ici2
        do j = jci1, jci2
          do n2 = 1, nsg
            ii = (i-1) * nsg + n2
            do n1 = 1, nsg
              jj = (j-1) * nsg + n1
              var2_3(jj,ii,ll) = var2_3(jj,ii,ll) + var3((n2-1)*nsg+n1,j,i)
            end do
          end do
        end do
      end do
    end if
  end subroutine reorder_add_subgrid_2d3d_real8

  subroutine reorder_add_subgrid_2d3d_real4(var3,var2_3,l,mask)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(in) :: var3
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(inout) :: var2_3
    integer(ik4), optional, intent(in) :: l
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(in), optional :: mask
    integer(ik4) :: i, j, ii, jj, n1, n2, ll
    ll = 1
    if ( present(l) ) ll = l
    if ( present(mask) ) then
      do i = ici1, ici2
        do j = jci1, jci2
          do n2 = 1, nsg
            ii = (i-1) * nsg + n2
            do n1 = 1, nsg
              jj = (j-1) * nsg + n1
              if ( mask((n2-1)*nsg+n1,j,i) > 0 ) then
                var2_3(jj,ii,ll) = var2_3(jj,ii,ll) + var3((n2-1)*nsg+n1,j,i)
              else
                var2_3(jj,ii,ll) = smissval
              end if
            end do
          end do
        end do
      end do
    else
      do i = ici1, ici2
        do j = jci1, jci2
          do n2 = 1, nsg
            ii = (i-1) * nsg + n2
            do n1 = 1, nsg
              jj = (j-1) * nsg + n1
              var2_3(jj,ii,ll) = var2_3(jj,ii,ll) + var3((n2-1)*nsg+n1,j,i)
            end do
          end do
        end do
      end do
    end if
  end subroutine reorder_add_subgrid_2d3d_real4

  subroutine reorder_subgrid_2d3d_real8(var3,var2_3,l,mask)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(in) :: var3
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(inout) :: var2_3
    integer(ik4), optional, intent(in) :: l
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(in), optional :: mask
    integer(ik4) :: i, j, ii, jj, n1, n2, ll
    ll = 1
    if ( present(l) ) ll = l
    if ( present(mask) ) then
      do i = ici1, ici2
        do j = jci1, jci2
          do n2 = 1, nsg
            ii = (i-1) * nsg + n2
            do n1 = 1, nsg
              jj = (j-1) * nsg + n1
              if ( mask((n2-1)*nsg+n1,j,i) > 0 ) then
                var2_3(jj,ii,ll) = var3((n2-1)*nsg+n1,j,i)
              else
                var2_3(jj,ii,ll) = dmissval
              end if
            end do
          end do
        end do
      end do
    else
      do i = ici1, ici2
        do j = jci1, jci2
          do n2 = 1, nsg
            ii = (i-1) * nsg + n2
            do n1 = 1, nsg
              jj = (j-1) * nsg + n1
              var2_3(jj,ii,ll) = var3((n2-1)*nsg+n1,j,i)
            end do
          end do
        end do
      end do
    end if
  end subroutine reorder_subgrid_2d3d_real8

  subroutine reorder_subgrid_2d3d_real4(var3,var2_3,l,mask)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(in) :: var3
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(inout) :: var2_3
    integer(ik4), optional, intent(in) :: l
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(in), optional :: mask
    integer(ik4) :: i, j, ii, jj, n1, n2, ll
    ll = 1
    if ( present(l) ) ll = l
    if ( present(mask) ) then
      do i = ici1, ici2
        do j = jci1, jci2
          do n2 = 1, nsg
            ii = (i-1) * nsg + n2
            do n1 = 1, nsg
              jj = (j-1) * nsg + n1
              if ( mask((n2-1)*nsg+n1,j,i) > 0 ) then
                var2_3(jj,ii,ll) = var3((n2-1)*nsg+n1,j,i)
              else
                var2_3(jj,ii,ll) = smissval
              end if
            end do
          end do
        end do
      end do
    else
      do i = ici1, ici2
        do j = jci1, jci2
          do n2 = 1, nsg
            ii = (i-1) * nsg + n2
            do n1 = 1, nsg
              jj = (j-1) * nsg + n1
              var2_3(jj,ii,ll) = var3((n2-1)*nsg+n1,j,i)
            end do
          end do
        end do
      end do
    end if
  end subroutine reorder_subgrid_2d3d_real4

  subroutine reorder_add_subgrid_3d_real8(var4,var2,l,mask)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:,:), intent(in) :: var4
    real(rk8), pointer, contiguous, dimension(:,:), intent(inout) :: var2
    integer(ik4), intent(in) :: l
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(in), optional :: mask
    integer(ik4) :: i, j, ii, jj, n1, n2
    if ( present(mask) ) then
      do i = ici1, ici2
        do j = jci1, jci2
          do n2 = 1, nsg
            ii = (i-1) * nsg + n2
            do n1 = 1, nsg
              jj = (j-1) * nsg + n1
              if ( mask((n2-1)*nsg+n1,j,i) > 0 ) then
                var2(jj,ii) = var2(jj,ii) + var4((n2-1)*nsg+n1,j,i,l)
              else
                var2(jj,ii) = dmissval
              end if
            end do
          end do
        end do
      end do
    else
      do i = ici1, ici2
        do j = jci1, jci2
          do n2 = 1, nsg
            ii = (i-1) * nsg + n2
            do n1 = 1, nsg
              jj = (j-1) * nsg + n1
              var2(jj,ii) = var2(jj,ii) + var4((n2-1)*nsg+n1,j,i,l)
            end do
          end do
        end do
      end do
    end if
  end subroutine reorder_add_subgrid_3d_real8

  subroutine reorder_add_subgrid_3d_real4(var4,var2,l,mask)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(in) :: var4
    real(rk4), pointer, contiguous, dimension(:,:), intent(inout) :: var2
    integer(ik4), intent(in) :: l
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(in), optional :: mask
    integer(ik4) :: i, j, ii, jj, n1, n2
    if ( present(mask) ) then
      do i = ici1, ici2
        do j = jci1, jci2
          do n2 = 1, nsg
            ii = (i-1) * nsg + n2
            do n1 = 1, nsg
              jj = (j-1) * nsg + n1
              if ( mask((n2-1)*nsg+n1,j,i) > 0 ) then
                var2(jj,ii) = var2(jj,ii) + var4((n2-1)*nsg+n1,j,i,l)
              else
                var2(jj,ii) = smissval
              end if
            end do
          end do
        end do
      end do
    else
      do i = ici1, ici2
        do j = jci1, jci2
          do n2 = 1, nsg
            ii = (i-1) * nsg + n2
            do n1 = 1, nsg
              jj = (j-1) * nsg + n1
              var2(jj,ii) = var2(jj,ii) + var4((n2-1)*nsg+n1,j,i,l)
            end do
          end do
        end do
      end do
    end if
  end subroutine reorder_add_subgrid_3d_real4

  subroutine reorder_subgrid_3d_real8(var4,var2,l,mask)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:,:), intent(in) :: var4
    real(rk8), pointer, contiguous, dimension(:,:), intent(inout) :: var2
    integer(ik4), intent(in) :: l
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(in), optional :: mask
    integer(ik4) :: i, j, ii, jj, n1, n2
    if ( present(mask) ) then
      do i = ici1, ici2
        do j = jci1, jci2
          do n2 = 1, nsg
            ii = (i-1) * nsg + n2
            do n1 = 1, nsg
              jj = (j-1) * nsg + n1
              if ( mask((n2-1)*nsg+n1,j,i) > 0 ) then
                var2(jj,ii) = var4((n2-1)*nsg+n1,j,i,l)
              else
                var2(jj,ii) = dmissval
              end if
            end do
          end do
        end do
      end do
    else
      do i = ici1, ici2
        do j = jci1, jci2
          do n2 = 1, nsg
            ii = (i-1) * nsg + n2
            do n1 = 1, nsg
              jj = (j-1) * nsg + n1
              var2(jj,ii) = var4((n2-1)*nsg+n1,j,i,l)
            end do
          end do
        end do
      end do
    end if
  end subroutine reorder_subgrid_3d_real8

  subroutine reorder_subgrid_3d_real4(var4,var2,l,mask)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(in) :: var4
    real(rk4), pointer, contiguous, dimension(:,:), intent(inout) :: var2
    integer(ik4), intent(in) :: l
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(in), optional :: mask
    integer(ik4) :: i, j, ii, jj, n1, n2
    if ( present(mask) ) then
      do i = ici1, ici2
        do j = jci1, jci2
          do n2 = 1, nsg
            ii = (i-1) * nsg + n2
            do n1 = 1, nsg
              jj = (j-1) * nsg + n1
              if ( mask((n2-1)*nsg+n1,j,i) > 0 ) then
                var2(jj,ii) = var4((n2-1)*nsg+n1,j,i,l)
              else
                var2(jj,ii) = smissval
              end if
            end do
          end do
        end do
      end do
    else
      do i = ici1, ici2
        do j = jci1, jci2
          do n2 = 1, nsg
            ii = (i-1) * nsg + n2
            do n1 = 1, nsg
              jj = (j-1) * nsg + n1
              var2(jj,ii) = var4((n2-1)*nsg+n1,j,i,l)
            end do
          end do
        end do
      end do
    end if
  end subroutine reorder_subgrid_3d_real4

  subroutine reorder_subgrid_4d_real8(var4,var3,mask)
    implicit none
    real(rk8), pointer, contiguous, dimension(:,:,:,:), intent(in) :: var4
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(inout) :: var3
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(in), optional :: mask
    integer(ik4) :: i, j, l, ii, jj, n1, n2
    if ( present(mask) ) then
      do l = 1, size(var4,4)
        do i = ici1, ici2
          do j = jci1, jci2
            do n2 = 1, nsg
              ii = (i-1) * nsg + n2
              do n1 = 1, nsg
                jj = (j-1) * nsg + n1
                if ( mask((n2-1)*nsg+n1,j,i) > 0 ) then
                  var3(jj,ii,l) = var4((n2-1)*nsg+n1,j,i,l)
                else
                  var3(jj,ii,l) = dmissval
                end if
              end do
            end do
          end do
        end do
      end do
    else
      do l = 1, size(var4,4)
        do i = ici1, ici2
          do j = jci1, jci2
            do n2 = 1, nsg
              ii = (i-1) * nsg + n2
              do n1 = 1, nsg
                jj = (j-1) * nsg + n1
                var3(jj,ii,l) = var4((n2-1)*nsg+n1,j,i,l)
              end do
            end do
          end do
        end do
      end do
    end if
  end subroutine reorder_subgrid_4d_real8

  subroutine reorder_subgrid_4d_real4(var4,var3,mask)
    implicit none
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(in) :: var4
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(inout) :: var3
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(in), optional :: mask
    integer(ik4) :: i, j, l, ii, jj, n1, n2
    if ( present(mask) ) then
      do l = 1, size(var4,4)
        do i = ici1, ici2
          do j = jci1, jci2
            do n2 = 1, nsg
              ii = (i-1) * nsg + n2
              do n1 = 1, nsg
                jj = (j-1) * nsg + n1
                if ( mask((n2-1)*nsg+n1,j,i) > 0 ) then
                  var3(jj,ii,l) = var4((n2-1)*nsg+n1,j,i,l)
                else
                  var3(jj,ii,l) = smissval
                end if
              end do
            end do
          end do
        end do
      end do
    else
      do l = 1, size(var4,4)
        do i = ici1, ici2
          do j = jci1, jci2
            do n2 = 1, nsg
              ii = (i-1) * nsg + n2
              do n1 = 1, nsg
                jj = (j-1) * nsg + n1
                var3(jj,ii,l) = var4((n2-1)*nsg+n1,j,i,l)
              end do
            end do
          end do
        end do
      end do
    end if
  end subroutine reorder_subgrid_4d_real4

  subroutine input_reorder_real8(m1,m2,j1,j2,i1,i2)
    implicit none
    real(rk8), dimension(:,:), intent(in) :: m1
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(inout) :: m2
    integer(ik4), intent(in) :: j1, j2, i1, i2
    integer(ik4) :: i, j, ii, jj, n1, n2
    do i = i1, i2
      do j = j1, j2
        do n2 = 1, nsg
          ii = (i-1) * nsg + n2
          do n1 = 1, nsg
            jj = (j-1) * nsg + n1
            m2((n2-1)*nsg+n1,j,i) = m1(jj,ii)
          end do
        end do
      end do
    end do
  end subroutine input_reorder_real8

  subroutine input_reorder_real4(m1,m2,j1,j2,i1,i2)
    implicit none
    real(rk4), dimension(:,:), intent(in) :: m1
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(inout) :: m2
    integer(ik4), intent(in) :: j1, j2, i1, i2
    integer(ik4) :: i, j, ii, jj, n1, n2
    do i = i1, i2
      do j = j1, j2
        do n2 = 1, nsg
          ii = (i-1) * nsg + n2
          do n1 = 1, nsg
            jj = (j-1) * nsg + n1
            m2((n2-1)*nsg+n1,j,i) = m1(jj,ii)
          end do
        end do
      end do
    end do
  end subroutine input_reorder_real4

  subroutine allsync
    implicit none
    call mpi_barrier(mycomm,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_barrier error.')
    end if
#endif
  end subroutine allsync

  subroutine clset(ncart_tot_g,ncart_tot_sg,cl)
    implicit none
    type(masked_comm), intent(inout) :: cl
    integer(ik4), intent(in) :: ncart_tot_g, ncart_tot_sg
    integer(ik4) :: linp, nrem, np, ntotg
    integer(ik4), dimension(1) :: tmp
    tmp(1) = ncart_tot_g
    call mpi_allgather(tmp,1,mpi_integer4,                   &
                       cl%cartesian_npoint_g,1,mpi_integer4, &
                       cartesian_communicator,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_allgather error.')
    end if
#endif
    cl%cartesian_displ_g(:) = 0
    do np = 2, nproc
      cl%cartesian_displ_g(np) = cl%cartesian_displ_g(np-1) + &
                                 cl%cartesian_npoint_g(np-1)
    end do
    cl%linear_npoint_g(:) = 0
    ntotg = sum(cl%cartesian_npoint_g)
    if ( nproc == 1 ) then
      cl%linear_npoint_g(1) = ntotg
    else if ( ntotg < nproc ) then
      cl%linear_npoint_g(2) = ntotg
    else
      linp = ntotg / nproc
      cl%linear_npoint_g(:) = linp
      nrem = ntotg - linp*nproc
      if ( nrem > 0 ) then
        np = 2
        do while (nrem > 0)
          cl%linear_npoint_g(np) = cl%linear_npoint_g(np) + 1
          nrem = nrem - 1
          np = np + 1
        end do
      end if
    end if
    cl%linear_displ_g(:) = 0
    do np = 2, nproc
      cl%linear_displ_g(np) = cl%linear_displ_g(np-1) + &
                              cl%linear_npoint_g(np-1)
    end do
    if ( nsg > 1 ) then
      tmp(1) = ncart_tot_sg
      call mpi_allgather(tmp,1,mpi_integer4,                    &
                         cl%cartesian_npoint_sg,1,mpi_integer4, &
                         cartesian_communicator,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_allgather error.')
      end if
#endif
      cl%cartesian_displ_sg(:) = 0
      do np = 2, nproc
        cl%cartesian_displ_sg(np) = cl%cartesian_displ_sg(np-1) + &
                                    cl%cartesian_npoint_sg(np-1)
      end do
      cl%linear_npoint_sg(:) = 0
      ntotg = sum(cl%cartesian_npoint_sg)
      if ( ntotg < nproc ) then
        cl%linear_npoint_sg(2) = ntotg
      else
        linp = ntotg / nproc
        cl%linear_npoint_sg(:) = linp
        nrem = ntotg - linp*nproc
        if ( nrem > 0 ) then
          np = 2
          do while (nrem > 0)
            cl%linear_npoint_sg(np) = cl%linear_npoint_sg(np) + 1
            nrem = nrem - 1
            np = np + 1
          end do
        end if
      end if
      cl%linear_displ_sg(:) = 0
      do np = 2, nproc
        cl%linear_displ_sg(np) = cl%linear_displ_sg(np-1) + &
                                 cl%linear_npoint_sg(np-1)
      end do
    else
      cl%cartesian_npoint_sg = cl%cartesian_npoint_g
      cl%cartesian_displ_sg = cl%cartesian_displ_g
      cl%linear_npoint_sg = cl%linear_npoint_g
      cl%linear_displ_sg = cl%linear_displ_g
    end if
    call getmem3d(r8subgrid,1,nnsg,jci1,jci2,ici1,ici2,'mpp:r8subgrid')
    call getmem3d(r4subgrid,1,nnsg,jci1,jci2,ici1,ici2,'mpp:r4subgrid')
    call getmem3d(i4subgrid,1,nnsg,jci1,jci2,ici1,ici2,'mpp:i4subgrid')
    call getmem3d(lsubgrid,1,nnsg,jci1,jci2,ici1,ici2,'mpp:lsubgrid')
    call getmem3d(global_r8subgrid,1,nnsg,jout1,jout2, &
            iout1,iout2,'mpp:global_r8subgrid')
    call getmem3d(global_r4subgrid,1,nnsg,jout1,jout2, &
            iout1,iout2,'mpp:global_r4subgrid')
    call getmem3d(global_i4subgrid,1,nnsg,jout1,jout2, &
            iout1,iout2,'mpp:global_i4subgrid')
    call getmem3d(global_lsubgrid,1,nnsg,jout1,jout2, &
            iout1,iout2,'mpp:global_lsubgrid')
    call getmem2d(global_r8grid,jout1,jout2,iout1,iout2,'mpp:global_r8grid')
    call getmem2d(global_r4grid,jout1,jout2,iout1,iout2,'mpp:global_r4grid')
    call getmem2d(global_i4grid,jout1,jout2,iout1,iout2,'mpp:global_i4grid')
    call getmem2d(global_lgrid,jout1,jout2,iout1,iout2,'mpp:global_lgrid')
  end subroutine clset

  subroutine cl_setup_real8(cl,gmask,sgmask,lrev)
    implicit none
    type(masked_comm), intent(inout) :: cl
    real(rk8), pointer, contiguous, dimension(:,:), intent(in) :: gmask
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(in) :: sgmask
    logical, optional, intent(in) :: lrev
    integer(ik4) :: ncart_tot_g, ncart_tot_sg
    if ( .not. associated(cl%linear_npoint_g) ) then
      call mpi_comm_dup(mycomm,cl%linear_communicator,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_comm_dup error.')
      end if
#endif
      call getmem1d(cl%linear_npoint_g,1,nproc,'cl:linear_npoint_g')
      call getmem1d(cl%linear_displ_g,1,nproc,'cl:linear_displ_g')
      call getmem1d(cl%cartesian_npoint_g,1,nproc,'cl:cartesian_npoint_g')
      call getmem1d(cl%cartesian_displ_g,1,nproc,'cl:cartesian_displ_g')
      call getmem1d(cl%linear_npoint_sg,1,nproc,'cl:linear_npoint_sg')
      call getmem1d(cl%linear_displ_sg,1,nproc,'cl:linear_displ_sg')
      call getmem1d(cl%cartesian_npoint_sg,1,nproc,'cl:cartesian_npoint_sg')
      call getmem1d(cl%cartesian_displ_sg,1,nproc,'cl:cartesian_displ_sg')
      call getmem2d(cl%gmask,jci1,jci2,ici1,ici2,'cl:gmask')
      call getmem3d(cl%sgmask,1,nnsg,jci1,jci2,ici1,ici2,'cl:sgmask')
      call getmem2d(cl%global_gmask,jout1,jout2,iout1,iout2,'cl:global_gmask')
      call getmem3d(cl%global_sgmask,1,nnsg,jout1,jout2, &
                                            iout1,iout2,'cl:global_sgmask')
      call getmem2d(cl%global_out_sgmask,joutsg1,joutsg2, &
                                         ioutsg1,ioutsg2,'cl:global_out_sgmask')
    end if
    cl%gmask = gmask(jci1:jci2,ici1:ici2) > 0.0D0
    cl%sgmask = sgmask(1:nnsg,jci1:jci2,ici1:ici2) > 0.0D0
    if ( present(lrev) ) then
      if ( lrev ) then
        cl%gmask = .not. cl%gmask
        cl%sgmask = .not. cl%sgmask
      end if
    end if
    call grid_collect(cl%gmask,cl%global_gmask,jci1,jci2,ici1,ici2)
    call subgrid_collect(cl%sgmask,cl%global_sgmask,jci1,jci2,ici1,ici2)
    call reorder_glb_subgrid(cl%global_sgmask,cl%global_out_sgmask)
    ncart_tot_g = count(cl%gmask)
    ncart_tot_sg = count(cl%sgmask)
    call clset(ncart_tot_g,ncart_tot_sg,cl)
  end subroutine cl_setup_real8

  subroutine cl_setup_real4(cl,gmask,sgmask,lrev)
    implicit none
    type(masked_comm), intent(inout) :: cl
    real(rk4), pointer, contiguous, dimension(:,:), intent(in) :: gmask
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(in) :: sgmask
    logical, optional, intent(in) :: lrev
    integer(ik4) :: ncart_tot_g, ncart_tot_sg
    if ( .not. associated(cl%linear_npoint_g) ) then
      call mpi_comm_dup(mycomm,cl%linear_communicator,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_comm_dup error.')
      end if
#endif
      call getmem1d(cl%linear_npoint_g,1,nproc,'cl:linear_npoint_g')
      call getmem1d(cl%linear_displ_g,1,nproc,'cl:linear_displ_g')
      call getmem1d(cl%cartesian_npoint_g,1,nproc,'cl:cartesian_npoint_g')
      call getmem1d(cl%cartesian_displ_g,1,nproc,'cl:cartesian_displ_g')
      call getmem1d(cl%linear_npoint_sg,1,nproc,'cl:linear_npoint_sg')
      call getmem1d(cl%linear_displ_sg,1,nproc,'cl:linear_displ_sg')
      call getmem1d(cl%cartesian_npoint_sg,1,nproc,'cl:cartesian_npoint_sg')
      call getmem1d(cl%cartesian_displ_sg,1,nproc,'cl:cartesian_displ_sg')
      call getmem2d(cl%gmask,jci1,jci2,ici1,ici2,'cl:gmask')
      call getmem3d(cl%sgmask,1,nnsg,jci1,jci2,ici1,ici2,'cl:sgmask')
      call getmem2d(cl%global_gmask,jout1,jout2,iout1,iout2,'cl:global_gmask')
      call getmem3d(cl%global_sgmask,1,nnsg,jout1,jout2, &
                                            iout1,iout2,'cl:global_sgmask')
      call getmem2d(cl%global_out_sgmask,joutsg1,joutsg2, &
                                         ioutsg1,ioutsg2,'cl:global_out_sgmask')
    end if
    cl%gmask = gmask(jci1:jci2,ici1:ici2) > 0.0D0
    cl%sgmask = sgmask(1:nnsg,jci1:jci2,ici1:ici2) > 0.0D0
    if ( present(lrev) ) then
      if ( lrev ) then
        cl%gmask = .not. cl%gmask
        cl%sgmask = .not. cl%sgmask
      end if
    end if
    call grid_collect(cl%gmask,cl%global_gmask,jci1,jci2,ici1,ici2)
    call subgrid_collect(cl%sgmask,cl%global_sgmask,jci1,jci2,ici1,ici2)
    call reorder_glb_subgrid(cl%global_sgmask,cl%global_out_sgmask)
    ncart_tot_g = count(cl%gmask)
    ncart_tot_sg = count(cl%sgmask)
    call clset(ncart_tot_g,ncart_tot_sg,cl)
  end subroutine cl_setup_real4

  subroutine cartesian_to_linear_logical_subgrid_subgrid(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    logical, pointer, contiguous, dimension(:,:,:), intent(in) :: matrix
    logical, pointer, contiguous, dimension(:), intent(inout) :: vector
    integer(ik4) :: nval, npt
    nval = cl%cartesian_npoint_sg(ccid+1)
    npt  = cl%linear_npoint_sg(myid+1)
    if ( nproc == 1 ) then
      call mypack(cl,matrix,vector)
      return
    end if
    if ( nval > 0 ) then
      call mypack(cl,matrix,lvector1)
    end if
    call mpi_gatherv(lvector1,nval,mpi_logical,                             &
                     lvector2,cl%cartesian_npoint_sg,cl%cartesian_displ_sg, &
                     mpi_logical,ccio,cartesian_communicator,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
#endif
    call mpi_scatterv(lvector2,cl%linear_npoint_sg,cl%linear_displ_sg, &
                      mpi_logical,vector,npt,mpi_logical,              &
                      iocpu,cl%linear_communicator,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_scatterv error.')
    end if
#endif
  end subroutine cartesian_to_linear_logical_subgrid_subgrid

  subroutine linear_to_cartesian_logical_subgrid_subgrid(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    logical, pointer, contiguous, dimension(:), intent(in) :: vector
    logical, pointer, contiguous, dimension(:,:,:), intent(inout) :: matrix
    integer(ik4) :: nval, npt
    nval = cl%cartesian_npoint_sg(ccid+1)
    npt  = cl%linear_npoint_sg(myid+1)
    if ( nproc == 1 ) then
      call myunpack(cl,vector,matrix)
      return
    end if
    call mpi_gatherv(vector,npt,mpi_logical,                          &
                     lvector2,cl%linear_npoint_sg,cl%linear_displ_sg, &
                     mpi_logical,iocpu,cl%linear_communicator,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
#endif
    call mpi_scatterv(lvector2,cl%cartesian_npoint_sg,   &
                      cl%cartesian_displ_sg,mpi_logical, &
                      lvector1,nval,mpi_logical,         &
                      ccio,cartesian_communicator,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_scatterv error.')
    end if
#endif
    if ( nval > 0 ) then
      call myunpack(cl,lvector1,matrix)
    end if
  end subroutine linear_to_cartesian_logical_subgrid_subgrid

  subroutine cartesian_to_linear_integer_subgrid_subgrid(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(in) :: matrix
    integer(ik4), pointer, contiguous, dimension(:), intent(inout) :: vector
    integer(ik4) :: nval, npt
    nval = cl%cartesian_npoint_sg(ccid+1)
    npt  = cl%linear_npoint_sg(myid+1)
    if ( nproc == 1 ) then
      call mypack(cl,matrix,vector)
      return
    end if
    if ( nval > 0 ) then
      call mypack(cl,matrix,i4vector1)
    end if
    call mpi_gatherv(i4vector1,nval,mpi_integer4,                            &
                     i4vector2,cl%cartesian_npoint_sg,cl%cartesian_displ_sg, &
                     mpi_integer4,ccio,cartesian_communicator,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
#endif
    call mpi_scatterv(i4vector2,cl%linear_npoint_sg,cl%linear_displ_sg, &
                      mpi_integer4,vector,npt,mpi_integer4,             &
                      iocpu,cl%linear_communicator,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_scatterv error.')
    end if
#endif
  end subroutine cartesian_to_linear_integer_subgrid_subgrid

  subroutine linear_to_cartesian_integer_subgrid_subgrid(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    integer(ik4), pointer, contiguous, dimension(:), intent(in) :: vector
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(inout) :: matrix
    integer(ik4) :: nval, npt
    nval = cl%cartesian_npoint_sg(ccid+1)
    npt  = cl%linear_npoint_sg(myid+1)
    if ( nproc == 1 ) then
      call myunpack(cl,vector,matrix)
      return
    end if
    call mpi_gatherv(vector,npt,mpi_integer4,                          &
                     i4vector2,cl%linear_npoint_sg,cl%linear_displ_sg, &
                     mpi_integer4,iocpu,cl%linear_communicator,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
#endif
    call mpi_scatterv(i4vector2,cl%cartesian_npoint_sg,   &
                      cl%cartesian_displ_sg,mpi_integer4, &
                      i4vector1,nval,mpi_integer4,        &
                      ccio,cartesian_communicator,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_scatterv error.')
    end if
#endif
    if ( nval > 0 ) then
      call myunpack(cl,i4vector1,matrix)
    end if
  end subroutine linear_to_cartesian_integer_subgrid_subgrid

  subroutine cartesian_to_linear_real8_subgrid_subgrid_4d(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk8), pointer, contiguous, dimension(:,:,:,:), intent(in) :: matrix
    real(rk8), pointer, contiguous, dimension(:,:), intent(inout) :: vector
    integer(ik4) :: nval, npt, nlev, k
    nval = cl%cartesian_npoint_sg(ccid+1)
    npt  = cl%linear_npoint_sg(myid+1)
    nlev = size(matrix,4)
    if ( nproc == 1 ) then
      call mypack(cl,matrix,vector,nlev)
      return
    end if
    do k = 1, nlev
      if ( nval > 0 ) then
        call mypack(cl,matrix,r8vector1,k)
      end if
      call mpi_gatherv(r8vector1,nval,mpi_real8,                               &
                       r8vector2,cl%cartesian_npoint_sg,cl%cartesian_displ_sg, &
                       mpi_real8,ccio,cartesian_communicator,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
      end if
#endif
      call mpi_scatterv(r8vector2,cl%linear_npoint_sg,cl%linear_displ_sg, &
                        mpi_real8,vector(:,k),npt,mpi_real8,              &
                        iocpu,cl%linear_communicator,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_scatterv error.')
      end if
#endif
    end do
  end subroutine cartesian_to_linear_real8_subgrid_subgrid_4d

  subroutine cartesian_to_linear_real4_subgrid_subgrid_4d(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(in) :: matrix
    real(rk4), pointer, contiguous, dimension(:,:), intent(inout) :: vector
    integer(ik4) :: nval, npt, nlev, k
    nval = cl%cartesian_npoint_sg(ccid+1)
    npt  = cl%linear_npoint_sg(myid+1)
    nlev = size(matrix,4)
    if ( nproc == 1 ) then
      call mypack(cl,matrix,vector,nlev)
      return
    end if
    do k = 1, nlev
      if ( nval > 0 ) then
        call mypack(cl,matrix,r4vector1,k)
      end if
      call mpi_gatherv(r4vector1,nval,mpi_real4,                               &
                       r4vector2,cl%cartesian_npoint_sg,cl%cartesian_displ_sg, &
                       mpi_real4,ccio,cartesian_communicator,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
      end if
#endif
      call mpi_scatterv(r4vector2,cl%linear_npoint_sg,cl%linear_displ_sg, &
                        mpi_real4,vector(:,k),npt,mpi_real4,              &
                        iocpu,cl%linear_communicator,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_scatterv error.')
      end if
#endif
    end do
  end subroutine cartesian_to_linear_real4_subgrid_subgrid_4d

  subroutine linear_to_cartesian_real8_subgrid_subgrid_4d(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk8), pointer, contiguous, dimension(:,:), intent(in) :: vector
    real(rk8), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: matrix
    integer(ik4) :: nval, npt, nlev, k
    nval = cl%cartesian_npoint_sg(ccid+1)
    npt  = cl%linear_npoint_sg(myid+1)
    nlev = size(vector,2)
    if ( nproc == 1 ) then
      call myunpack(cl,vector,matrix,nlev)
      return
    end if
    do k = 1, nlev
      call mpi_gatherv(vector(:,k),npt,mpi_real8,                        &
                       r8vector2,cl%linear_npoint_sg,cl%linear_displ_sg, &
                       mpi_real8,iocpu,cl%linear_communicator,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
      end if
#endif
      call mpi_scatterv(r8vector2,cl%cartesian_npoint_sg, &
                        cl%cartesian_displ_sg,mpi_real8,  &
                        r8vector1,nval,mpi_real8,         &
                        ccio,cartesian_communicator,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_scatterv error.')
      end if
#endif
      if ( nval > 0 ) then
        call myunpack(cl,r8vector1,matrix,k)
      end if
    end do
  end subroutine linear_to_cartesian_real8_subgrid_subgrid_4d

  subroutine linear_to_cartesian_real4_subgrid_subgrid_4d(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk4), pointer, contiguous, dimension(:,:), intent(in) :: vector
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: matrix
    integer(ik4) :: nval, npt, nlev, k
    nval = cl%cartesian_npoint_sg(ccid+1)
    npt  = cl%linear_npoint_sg(myid+1)
    nlev = size(vector,2)
    if ( nproc == 1 ) then
      call myunpack(cl,vector,matrix,nlev)
      return
    end if
    do k = 1, nlev
      call mpi_gatherv(vector(:,k),npt,mpi_real4,                        &
                       r4vector2,cl%linear_npoint_sg,cl%linear_displ_sg, &
                       mpi_real4,iocpu,cl%linear_communicator,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
      end if
#endif
      call mpi_scatterv(r4vector2,cl%cartesian_npoint_sg, &
                        cl%cartesian_displ_sg,mpi_real4,  &
                        r4vector1,nval,mpi_real4,         &
                        ccio,cartesian_communicator,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_scatterv error.')
      end if
#endif
      if ( nval > 0 ) then
        call myunpack(cl,r4vector1,matrix,k)
      end if
    end do
  end subroutine linear_to_cartesian_real4_subgrid_subgrid_4d

  subroutine cartesian_to_linear_real8_subgrid_subgrid(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(in) :: matrix
    real(rk8), pointer, contiguous, dimension(:), intent(inout) :: vector
    integer(ik4) :: nval, npt
    nval = cl%cartesian_npoint_sg(ccid+1)
    npt  = cl%linear_npoint_sg(myid+1)
    if ( nproc == 1 ) then
      call mypack(cl,matrix,vector)
      return
    end if
    if ( nval > 0 ) then
      call mypack(cl,matrix,r8vector1)
    end if
    call mpi_gatherv(r8vector1,nval,mpi_real8,                               &
                     r8vector2,cl%cartesian_npoint_sg,cl%cartesian_displ_sg, &
                     mpi_real8,ccio,cartesian_communicator,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
#endif
    call mpi_scatterv(r8vector2,cl%linear_npoint_sg,cl%linear_displ_sg, &
                      mpi_real8,vector,npt,mpi_real8,                   &
                      iocpu,cl%linear_communicator,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_scatterv error.')
    end if
#endif
  end subroutine cartesian_to_linear_real8_subgrid_subgrid

  subroutine cartesian_to_linear_real4_subgrid_subgrid(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(in) :: matrix
    real(rk4), pointer, contiguous, dimension(:), intent(inout) :: vector
    integer(ik4) :: nval, npt
    nval = cl%cartesian_npoint_sg(ccid+1)
    npt  = cl%linear_npoint_sg(myid+1)
    if ( nproc == 1 ) then
      call mypack(cl,matrix,vector)
      return
    end if
    if ( nval > 0 ) then
      call mypack(cl,matrix,r4vector1)
    end if
    call mpi_gatherv(r4vector1,nval,mpi_real4,                               &
                     r4vector2,cl%cartesian_npoint_sg,cl%cartesian_displ_sg, &
                     mpi_real4,ccio,cartesian_communicator,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
#endif
    call mpi_scatterv(r4vector2,cl%linear_npoint_sg,cl%linear_displ_sg, &
                      mpi_real4,vector,npt,mpi_real4,                   &
                      iocpu,cl%linear_communicator,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_scatterv error.')
    end if
#endif
  end subroutine cartesian_to_linear_real4_subgrid_subgrid

  subroutine linear_to_cartesian_real8_subgrid_subgrid(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk8), pointer, contiguous, dimension(:), intent(in) :: vector
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(inout) :: matrix
    integer(ik4) :: nval, npt
    nval = cl%cartesian_npoint_sg(ccid+1)
    npt  = cl%linear_npoint_sg(myid+1)
    if ( nproc == 1 ) then
      call myunpack(cl,vector,matrix)
      return
    end if
    call mpi_gatherv(vector,npt,mpi_real8,                             &
                     r8vector2,cl%linear_npoint_sg,cl%linear_displ_sg, &
                     mpi_real8,iocpu,cl%linear_communicator,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
#endif
    call mpi_scatterv(r8vector2,cl%cartesian_npoint_sg, &
                      cl%cartesian_displ_sg,mpi_real8,  &
                      r8vector1,nval,mpi_real8,         &
                      ccio,cartesian_communicator,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_scatterv error.')
    end if
#endif
    if ( nval > 0 ) then
      call myunpack(cl,r8vector1,matrix)
    end if
  end subroutine linear_to_cartesian_real8_subgrid_subgrid

  subroutine linear_to_cartesian_real4_subgrid_subgrid(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk4), pointer, contiguous, dimension(:), intent(in) :: vector
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(inout) :: matrix
    integer(ik4) :: nval, npt
    nval = cl%cartesian_npoint_sg(ccid+1)
    npt  = cl%linear_npoint_sg(myid+1)
    if ( nproc == 1 ) then
      call myunpack(cl,vector,matrix)
      return
    end if
    call mpi_gatherv(vector,npt,mpi_real4,                             &
                     r4vector2,cl%linear_npoint_sg,cl%linear_displ_sg, &
                     mpi_real4,iocpu,cl%linear_communicator,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
#endif
    call mpi_scatterv(r4vector2,cl%cartesian_npoint_sg, &
                      cl%cartesian_displ_sg,mpi_real4,  &
                      r4vector1,nval,mpi_real4,         &
                      ccio,cartesian_communicator,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_scatterv error.')
    end if
#endif
    if ( nval > 0 ) then
      call myunpack(cl,r4vector1,matrix)
    end if
  end subroutine linear_to_cartesian_real4_subgrid_subgrid

  subroutine cartesian_to_linear_logical_grid_subgrid(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    logical, pointer, contiguous, dimension(:,:), intent(in) :: matrix
    logical, pointer, contiguous, dimension(:), intent(inout) :: vector
    integer(ik4) :: i, j, n
    do concurrent ( n = 1:nnsg, j = jci1:jci2, i = ici1:ici2 )
      lsubgrid(n,j,i) = matrix(j,i)
    end do
    call cartesian_to_linear_logical_subgrid_subgrid(cl,lsubgrid,vector)
  end subroutine cartesian_to_linear_logical_grid_subgrid

  subroutine cartesian_to_linear_integer_grid_subgrid(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    integer(ik4), pointer, contiguous, dimension(:,:), intent(in) :: matrix
    integer(ik4), pointer, contiguous, dimension(:), intent(inout) :: vector
    integer(ik4) :: i, j, n
    do concurrent ( n = 1:nnsg, j = jci1:jci2, i = ici1:ici2 )
      i4subgrid(n,j,i) = matrix(j,i)
    end do
    call cartesian_to_linear_integer_subgrid_subgrid(cl,i4subgrid,vector)
  end subroutine cartesian_to_linear_integer_grid_subgrid

  subroutine cartesian_to_linear_real8_grid_subgrid(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk8), pointer, contiguous, dimension(:,:), intent(in) :: matrix
    real(rk8), pointer, contiguous, dimension(:), intent(inout) :: vector
    integer(ik4) :: i, j, n
    do concurrent ( n = 1:nnsg, j = jci1:jci2, i = ici1:ici2 )
      r8subgrid(n,j,i) = matrix(j,i)
    end do
    call cartesian_to_linear_real8_subgrid_subgrid(cl,r8subgrid,vector)
  end subroutine cartesian_to_linear_real8_grid_subgrid

  subroutine cartesian_to_linear_real4_grid_subgrid(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk4), pointer, contiguous, dimension(:,:), intent(in) :: matrix
    real(rk4), pointer, contiguous, dimension(:), intent(inout) :: vector
    integer(ik4) :: i, j, n
    do concurrent ( n = 1:nnsg, j = jci1:jci2, i = ici1:ici2 )
      r4subgrid(n,j,i) = matrix(j,i)
    end do
    call cartesian_to_linear_real4_subgrid_subgrid(cl,r4subgrid,vector)
  end subroutine cartesian_to_linear_real4_grid_subgrid

  subroutine global_to_linear_logical_subgrid_subgrid(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    logical, pointer, contiguous, dimension(:,:,:), intent(in) :: matrix
    logical, pointer, contiguous, dimension(:), intent(inout) :: vector
    integer(ik4) :: npt

    if ( nproc == 1 ) then
      call mypack_global(cl,matrix,vector)
      return
    end if
    call subgrid_collect(matrix,global_lsubgrid,jci1,jci2,ici1,ici2)
    call mypack_global(cl,global_lsubgrid,lvector1)
    npt  = cl%linear_npoint_sg(myid+1)
    call mpi_scatterv(lvector1,cl%linear_npoint_sg,cl%linear_displ_sg, &
                      mpi_logical,vector,npt,mpi_logical,              &
                      iocpu,cl%linear_communicator,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_scatterv error.')
    end if
#endif
  end subroutine global_to_linear_logical_subgrid_subgrid

  subroutine linear_to_global_logical_subgrid_subgrid(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    logical, pointer, contiguous, dimension(:), intent(in) :: vector
    logical, pointer, contiguous, dimension(:,:,:), intent(inout) :: matrix
    if ( nproc == 1 ) then
      call myunpack_global(cl,vector,matrix)
      return
    end if
    call mpi_gatherv(vector,cl%linear_npoint_sg(myid+1),mpi_logical,  &
                     lvector1,cl%linear_npoint_sg,cl%linear_displ_sg, &
                     mpi_logical,iocpu,cl%linear_communicator,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
#endif
    call myunpack_global(cl,lvector1,global_lsubgrid)
    call subgrid_distribute(global_lsubgrid,matrix, &
                            jci1,jci2,ici1,ici2,cl%sgmask)
  end subroutine linear_to_global_logical_subgrid_subgrid

  subroutine global_to_linear_integer_subgrid_subgrid(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(in) :: matrix
    integer(ik4), pointer, contiguous, dimension(:), intent(inout) :: vector
    integer(ik4) :: npt
    if ( nproc == 1 ) then
      call mypack_global(cl,matrix,vector)
      return
    end if
    call subgrid_collect(matrix,global_i4subgrid,jci1,jci2,ici1,ici2)
    call mypack_global(cl,global_i4subgrid,i4vector1)
    npt  = cl%linear_npoint_sg(myid+1)
    call mpi_scatterv(i4vector1,cl%linear_npoint_sg,cl%linear_displ_sg, &
                      mpi_integer4,vector,npt,mpi_integer4,             &
                      iocpu,cl%linear_communicator,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_scatterv error.')
    end if
#endif
  end subroutine global_to_linear_integer_subgrid_subgrid

  subroutine linear_to_global_integer_subgrid_subgrid(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    integer(ik4), pointer, contiguous, dimension(:), intent(in) :: vector
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(inout) :: matrix
    if ( nproc == 1 ) then
      call myunpack_global(cl,vector,matrix)
      return
    end if
    call mpi_gatherv(vector,cl%linear_npoint_sg(myid+1),mpi_integer4,  &
                     i4vector1,cl%linear_npoint_sg,cl%linear_displ_sg, &
                     mpi_integer4,iocpu,cl%linear_communicator,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
#endif
    call myunpack_global(cl,i4vector1,global_i4subgrid)
    call subgrid_distribute(global_i4subgrid,matrix, &
            jci1,jci2,ici1,ici2,cl%sgmask)
  end subroutine linear_to_global_integer_subgrid_subgrid

  subroutine global_to_linear_real8_subgrid_subgrid(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(in) :: matrix
    real(rk8), pointer, contiguous, dimension(:), intent(inout) :: vector
    integer(ik4) :: npt
    if ( nproc == 1 ) then
      call mypack_global(cl,matrix,vector)
      return
    end if
    call subgrid_collect(matrix,global_r8subgrid,jci1,jci2,ici1,ici2)
    call mypack_global(cl,global_r8subgrid,r8vector1)
    npt  = cl%linear_npoint_sg(myid+1)
    call mpi_scatterv(r8vector1,cl%linear_npoint_sg,cl%linear_displ_sg, &
                      mpi_real8,vector,npt,mpi_real8,                   &
                      iocpu,cl%linear_communicator,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_scatterv error.')
    end if
#endif
  end subroutine global_to_linear_real8_subgrid_subgrid

  subroutine global_to_linear_real4_subgrid_subgrid(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(in) :: matrix
    real(rk4), pointer, contiguous, dimension(:), intent(inout) :: vector
    integer(ik4) :: npt
    if ( nproc == 1 ) then
      call mypack_global(cl,matrix,vector)
      return
    end if
    call subgrid_collect(matrix,global_r4subgrid,jci1,jci2,ici1,ici2)
    call mypack_global(cl,global_r4subgrid,r4vector1)
    npt  = cl%linear_npoint_sg(myid+1)
    call mpi_scatterv(r4vector1,cl%linear_npoint_sg,cl%linear_displ_sg, &
                      mpi_real4,vector,npt,mpi_real4,                   &
                      iocpu,cl%linear_communicator,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_scatterv error.')
    end if
#endif
  end subroutine global_to_linear_real4_subgrid_subgrid

  subroutine global_to_linear_real4_real8_subgrid_subgrid(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(in) :: matrix
    real(rk8), pointer, contiguous, dimension(:), intent(inout) :: vector
    integer(ik4) :: npt
    if ( nproc == 1 ) then
      call mypack_global(cl,matrix,vector)
      return
    end if
    call subgrid_collect(matrix,global_r4subgrid,jci1,jci2,ici1,ici2)
    call mypack_global(cl,global_r4subgrid,r8vector1)
    npt  = cl%linear_npoint_sg(myid+1)
    call mpi_scatterv(r8vector1,cl%linear_npoint_sg,cl%linear_displ_sg, &
                      mpi_real8,vector,npt,mpi_real8,                   &
                      iocpu,cl%linear_communicator,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_scatterv error.')
    end if
#endif
  end subroutine global_to_linear_real4_real8_subgrid_subgrid

  subroutine linear_to_global_real8_subgrid_subgrid(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk8), pointer, contiguous, dimension(:), intent(in) :: vector
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(inout) :: matrix
    if ( nproc == 1 ) then
      call myunpack_global(cl,vector,matrix)
      return
    end if
    call mpi_gatherv(vector,cl%linear_npoint_sg(myid+1),mpi_real8,  &
                     r8vector1,cl%linear_npoint_sg,cl%linear_displ_sg, &
                     mpi_real8,iocpu,cl%linear_communicator,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
#endif
    call myunpack_global(cl,r8vector1,global_r8subgrid)
    call subgrid_distribute(global_r8subgrid,matrix, &
            jci1,jci2,ici1,ici2,cl%sgmask)
  end subroutine linear_to_global_real8_subgrid_subgrid

  subroutine linear_to_global_real4_subgrid_subgrid(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk4), pointer, contiguous, dimension(:), intent(in) :: vector
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(inout) :: matrix
    if ( nproc == 1 ) then
      call myunpack_global(cl,vector,matrix)
      return
    end if
    call mpi_gatherv(vector,cl%linear_npoint_sg(myid+1),mpi_real4,  &
                     r4vector1,cl%linear_npoint_sg,cl%linear_displ_sg, &
                     mpi_real4,iocpu,cl%linear_communicator,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
#endif
    call myunpack_global(cl,r4vector1,global_r4subgrid)
    call subgrid_distribute(global_r4subgrid,matrix, &
            jci1,jci2,ici1,ici2,cl%sgmask)
  end subroutine linear_to_global_real4_subgrid_subgrid

  subroutine linear_to_global_real8_real4_subgrid_subgrid(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk8), pointer, contiguous, dimension(:), intent(in) :: vector
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(inout) :: matrix
    if ( nproc == 1 ) then
      call myunpack_global(cl,vector,matrix)
      return
    end if
    call mpi_gatherv(vector,cl%linear_npoint_sg(myid+1),mpi_real8,     &
                     r8vector1,cl%linear_npoint_sg,cl%linear_displ_sg, &
                     mpi_real8,iocpu,cl%linear_communicator,mpierr)
#ifdef DEBUG
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
    end if
#endif
    call myunpack_global(cl,r8vector1,global_r4subgrid)
    call subgrid_distribute(global_r4subgrid,matrix, &
            jci1,jci2,ici1,ici2,cl%sgmask)
  end subroutine linear_to_global_real8_real4_subgrid_subgrid

  subroutine global_to_linear_real8_subgrid_subgrid_4d(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk8), pointer, contiguous, dimension(:,:,:,:), intent(in) :: matrix
    real(rk8), pointer, contiguous, dimension(:,:), intent(inout) :: vector
    integer(ik4) :: npt, k, nlev
    nlev = size(matrix,4)
    if ( nproc == 1 ) then
      call mypack_global(cl,matrix,vector,nlev)
      return
    end if
    npt  = cl%linear_npoint_sg(myid+1)
    do k = 1, nlev
      r8subgrid = matrix(:,jci1:jci2,ici1:ici2,k)
      call subgrid_collect(r8subgrid,global_r8subgrid,jci1,jci2,ici1,ici2)
      call mypack_global(cl,global_r8subgrid,r8vector1)
      call mpi_scatterv(r8vector1,cl%linear_npoint_sg,cl%linear_displ_sg, &
                        mpi_real8,vector(:,k),npt,mpi_real8,                   &
                        iocpu,cl%linear_communicator,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_scatterv error.')
      end if
#endif
    end do
  end subroutine global_to_linear_real8_subgrid_subgrid_4d

  subroutine global_to_linear_real4_subgrid_subgrid_4d(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(in) :: matrix
    real(rk4), pointer, contiguous, dimension(:,:), intent(inout) :: vector
    integer(ik4) :: npt, k, nlev
    nlev = size(matrix,4)
    if ( nproc == 1 ) then
      call mypack_global(cl,matrix,vector,nlev)
      return
    end if
    npt  = cl%linear_npoint_sg(myid+1)
    do k = 1, nlev
      r4subgrid = matrix(:,jci1:jci2,ici1:ici2,k)
      call subgrid_collect(r4subgrid,global_r4subgrid,jci1,jci2,ici1,ici2)
      call mypack_global(cl,global_r4subgrid,r4vector1)
      call mpi_scatterv(r4vector1,cl%linear_npoint_sg,cl%linear_displ_sg, &
                        mpi_real4,vector(:,k),npt,mpi_real4,              &
                        iocpu,cl%linear_communicator,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_scatterv error.')
      end if
#endif
    end do
  end subroutine global_to_linear_real4_subgrid_subgrid_4d

  subroutine global_to_linear_real4_real8_subgrid_subgrid_4d(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(in) :: matrix
    real(rk8), pointer, contiguous, dimension(:,:), intent(inout) :: vector
    integer(ik4) :: npt, k, nlev
    nlev = size(matrix,4)
    if ( nproc == 1 ) then
      call mypack_global(cl,matrix,vector,nlev)
      return
    end if
    npt  = cl%linear_npoint_sg(myid+1)
    do k = 1, nlev
      r4subgrid = matrix(:,jci1:jci2,ici1:ici2,k)
      call subgrid_collect(r4subgrid,global_r4subgrid,jci1,jci2,ici1,ici2)
      call mypack_global(cl,global_r4subgrid,r8vector1)
      call mpi_scatterv(r8vector1,cl%linear_npoint_sg,cl%linear_displ_sg, &
                        mpi_real8,vector(:,k),npt,mpi_real8,              &
                        iocpu,cl%linear_communicator,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_scatterv error.')
      end if
#endif
    end do
  end subroutine global_to_linear_real4_real8_subgrid_subgrid_4d

  subroutine linear_to_global_real8_subgrid_subgrid_4d(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk8), pointer, contiguous, dimension(:,:), intent(in) :: vector
    real(rk8), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: matrix
    integer(ik4) :: npt, nlev, k
    nlev = size(vector,2)
    if ( nproc == 1 ) then
      call myunpack(cl,vector,matrix,nlev)
      return
    end if
    npt  = cl%linear_npoint_sg(myid+1)
    do k = 1, nlev
      call mpi_gatherv(vector(:,k),npt,mpi_real8,                        &
                       r8vector1,cl%linear_npoint_sg,cl%linear_displ_sg, &
                       mpi_real8,iocpu,cl%linear_communicator,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
      end if
#endif
      call myunpack_global(cl,r8vector1,global_r8subgrid)
      call subgrid_distribute(global_r8subgrid,r8subgrid, &
              jci1,jci2,ici1,ici2,cl%sgmask)
      where ( cl%sgmask )
        matrix(:,jci1:jci2,ici1:ici2,k) = r8subgrid
      end where
    end do
  end subroutine linear_to_global_real8_subgrid_subgrid_4d

  subroutine linear_to_global_real4_subgrid_subgrid_4d(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk4), pointer, contiguous, dimension(:,:), intent(in) :: vector
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: matrix
    integer(ik4) :: npt, nlev, k
    nlev = size(vector,2)
    if ( nproc == 1 ) then
      call myunpack(cl,vector,matrix,nlev)
      return
    end if
    npt  = cl%linear_npoint_sg(myid+1)
    do k = 1, nlev
      call mpi_gatherv(vector(:,k),npt,mpi_real4,                        &
                       r4vector1,cl%linear_npoint_sg,cl%linear_displ_sg, &
                       mpi_real4,iocpu,cl%linear_communicator,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
      end if
#endif
      call myunpack_global(cl,r4vector1,global_r4subgrid)
      call subgrid_distribute(global_r4subgrid,r4subgrid, &
              jci1,jci2,ici1,ici2,cl%sgmask)
      where ( cl%sgmask )
        matrix(:,jci1:jci2,ici1:ici2,k) = r4subgrid
      end where
    end do
  end subroutine linear_to_global_real4_subgrid_subgrid_4d

  subroutine linear_to_global_real8_real4_subgrid_subgrid_4d(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk8), pointer, contiguous, dimension(:,:), intent(in) :: vector
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: matrix
    integer(ik4) :: npt, nlev, k
    nlev = size(vector,2)
    if ( nproc == 1 ) then
      call myunpack(cl,vector,matrix,nlev)
      return
    end if
    npt  = cl%linear_npoint_sg(myid+1)
    do k = 1, nlev
      call mpi_gatherv(vector(:,k),npt,mpi_real8,                        &
                       r8vector1,cl%linear_npoint_sg,cl%linear_displ_sg, &
                       mpi_real8,iocpu,cl%linear_communicator,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_gatherv error.')
      end if
#endif
      call myunpack_global(cl,r8vector1,global_r4subgrid)
      call subgrid_distribute(global_r4subgrid,r4subgrid, &
              jci1,jci2,ici1,ici2,cl%sgmask)
      where ( cl%sgmask )
        matrix(:,jci1:jci2,ici1:ici2,k) = r4subgrid
      end where
    end do
  end subroutine linear_to_global_real8_real4_subgrid_subgrid_4d

  subroutine global_to_linear_logical_grid_subgrid(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    logical, pointer, contiguous, dimension(:,:), intent(in) :: matrix
    logical, pointer, contiguous, dimension(:), intent(inout) :: vector
    integer(ik4) :: i, j, n
    do concurrent ( n = 1:nnsg, j = jci1:jci2, i = ici1:ici2 )
      lsubgrid(n,j,i) = matrix(j,i)
    end do
    call global_to_linear_logical_subgrid_subgrid(cl,lsubgrid,vector)
  end subroutine global_to_linear_logical_grid_subgrid

  subroutine global_to_linear_integer_grid_subgrid(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    integer(ik4), pointer, contiguous, dimension(:,:), intent(in) :: matrix
    integer(ik4), pointer, contiguous, dimension(:), intent(inout) :: vector
    integer(ik4) :: i, j, n
    do concurrent ( n = 1:nnsg, j = jci1:jci2, i = ici1:ici2 )
      i4subgrid(n,j,i) = matrix(j,i)
    end do
    call global_to_linear_integer_subgrid_subgrid(cl,i4subgrid,vector)
  end subroutine global_to_linear_integer_grid_subgrid

  subroutine global_to_linear_real8_grid_subgrid(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk8), pointer, contiguous, dimension(:,:), intent(in) :: matrix
    real(rk8), pointer, contiguous, dimension(:), intent(inout) :: vector
    integer(ik4) :: i, j, n
    do concurrent ( n = 1:nnsg, j = jci1:jci2, i = ici1:ici2 )
      r8subgrid(n,j,i) = matrix(j,i)
    end do
    call global_to_linear_real8_subgrid_subgrid(cl,r8subgrid,vector)
  end subroutine global_to_linear_real8_grid_subgrid

  subroutine global_to_linear_real4_grid_subgrid(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk4), pointer, contiguous, dimension(:,:), intent(in) :: matrix
    real(rk4), pointer, contiguous, dimension(:), intent(inout) :: vector
    integer(ik4) :: i, j, n
    do concurrent ( n = 1:nnsg, j = jci1:jci2, i = ici1:ici2 )
      r4subgrid(n,j,i) = matrix(j,i)
    end do
    call global_to_linear_real4_subgrid_subgrid(cl,r4subgrid,vector)
  end subroutine global_to_linear_real4_grid_subgrid

  subroutine global_to_linear_real4_real8_grid_subgrid(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk4), pointer, contiguous, dimension(:,:), intent(in) :: matrix
    real(rk8), pointer, contiguous, dimension(:), intent(inout) :: vector
    integer(ik4) :: i, j, n
    do concurrent ( n = 1:nnsg, j = jci1:jci2, i = ici1:ici2 )
      r4subgrid(n,j,i) = matrix(j,i)
    end do
    call global_to_linear_real4_real8_subgrid_subgrid(cl,r4subgrid,vector)
  end subroutine global_to_linear_real4_real8_grid_subgrid

  subroutine cl_dispose(cl)
    implicit none
    type(masked_comm), intent(inout) :: cl
    if ( associated(cl%linear_npoint_g) ) then
      call relmem1d(cl%linear_npoint_g)
      call relmem1d(cl%linear_displ_g)
      call relmem1d(cl%cartesian_npoint_g)
      call relmem1d(cl%cartesian_displ_g)
      call relmem1d(cl%linear_npoint_sg)
      call relmem1d(cl%linear_displ_sg)
      call relmem1d(cl%cartesian_npoint_sg)
      call relmem1d(cl%cartesian_displ_sg)
      call relmem2d(cl%gmask)
      call relmem3d(cl%sgmask)
      call relmem2d(cl%global_gmask)
      call relmem3d(cl%global_sgmask)
      call mpi_comm_free(cl%linear_communicator,mpierr)
#ifdef DEBUG
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_comm_free error.')
      end if
#endif
    end if
  end subroutine cl_dispose

  subroutine mypack_logical_grid(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    logical, pointer, contiguous, dimension(:), intent(inout) :: vector
    logical, pointer, contiguous, dimension(:,:), intent(in) :: matrix
    integer(ik4) :: i, j, iv
    iv = 1
    do i = ici1, ici2
      do j = jci1, jci2
        if ( cl%gmask(j,i) ) then
          vector(iv) = matrix(j,i)
          iv = iv + 1
        end if
      end do
    end do
  end subroutine mypack_logical_grid

  subroutine mypack_logical_subgrid(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    logical, pointer, contiguous, dimension(:), intent(inout) :: vector
    logical, pointer, contiguous, dimension(:,:,:), intent(in) :: matrix
    integer(ik4) :: i, j, n, iv
    iv = 1
    do i = ici1, ici2
      do j = jci1, jci2
        do n = 1, nnsg
          if ( cl%sgmask(n,j,i) ) then
            vector(iv) = matrix(n,j,i)
            iv = iv + 1
          end if
        end do
      end do
    end do
  end subroutine mypack_logical_subgrid

  subroutine mypack_integer_grid(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    integer(ik4), pointer, contiguous, dimension(:), intent(inout) :: vector
    integer(ik4), pointer, contiguous, dimension(:,:), intent(in) :: matrix
    integer(ik4) :: i, j, iv
    iv = 1
    do i = ici1, ici2
      do j = jci1, jci2
        if ( cl%gmask(j,i) ) then
          vector(iv) = matrix(j,i)
          iv = iv + 1
        end if
      end do
    end do
  end subroutine mypack_integer_grid

  subroutine mypack_integer_subgrid(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    integer(ik4), pointer, contiguous, dimension(:), intent(inout) :: vector
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(in) :: matrix
    integer(ik4) :: i, j, n, iv
    iv = 1
    do i = ici1, ici2
      do j = jci1, jci2
        do n = 1, nnsg
          if ( cl%sgmask(n,j,i) ) then
            vector(iv) = matrix(n,j,i)
            iv = iv + 1
          end if
        end do
      end do
    end do
  end subroutine mypack_integer_subgrid

  subroutine mypack_real8_grid(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk8), pointer, contiguous, dimension(:), intent(inout) :: vector
    real(rk8), pointer, contiguous, dimension(:,:), intent(in) :: matrix
    integer(ik4) :: i, j, iv
    iv = 1
    do i = ici1, ici2
      do j = jci1, jci2
        if ( cl%gmask(j,i) ) then
          vector(iv) = matrix(j,i)
          iv = iv + 1
        end if
      end do
    end do
  end subroutine mypack_real8_grid

  subroutine mypack_real4_grid(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk4), pointer, contiguous, dimension(:), intent(inout) :: vector
    real(rk4), pointer, contiguous, dimension(:,:), intent(in) :: matrix
    integer(ik4) :: i, j, iv
    iv = 1
    do i = ici1, ici2
      do j = jci1, jci2
        if ( cl%gmask(j,i) ) then
          vector(iv) = matrix(j,i)
          iv = iv + 1
        end if
      end do
    end do
  end subroutine mypack_real4_grid

  subroutine mypack_real8_subgrid(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk8), pointer, contiguous, dimension(:), intent(inout) :: vector
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(in) :: matrix
    integer(ik4) :: i, j, n, iv
    iv = 1
    do i = ici1, ici2
      do j = jci1, jci2
        do n = 1, nnsg
          if ( cl%sgmask(n,j,i) ) then
            vector(iv) = matrix(n,j,i)
            iv = iv + 1
          end if
        end do
      end do
    end do
  end subroutine mypack_real8_subgrid

  subroutine mypack_real4_subgrid(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk4), pointer, contiguous, dimension(:), intent(inout) :: vector
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(in) :: matrix
    integer(ik4) :: i, j, n, iv
    iv = 1
    do i = ici1, ici2
      do j = jci1, jci2
        do n = 1, nnsg
          if ( cl%sgmask(n,j,i) ) then
            vector(iv) = matrix(n,j,i)
            iv = iv + 1
          end if
        end do
      end do
    end do
  end subroutine mypack_real4_subgrid

  subroutine mypack_real8_subgrid_4d(cl,matrix,vector,klev)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk8), pointer, contiguous, dimension(:,:), intent(inout) :: vector
    real(rk8), pointer, contiguous, dimension(:,:,:,:), intent(in) :: matrix
    integer(ik4), intent(in) :: klev
    integer(ik4) :: i, j, k, n, iv
    do k = 1, klev
      iv = 1
      do i = ici1, ici2
        do j = jci1, jci2
          do n = 1, nnsg
            if ( cl%sgmask(n,j,i) ) then
              vector(iv,k) = matrix(n,j,i,k)
              iv = iv + 1
            end if
          end do
        end do
      end do
    end do
  end subroutine mypack_real8_subgrid_4d

  subroutine mypack_real4_subgrid_4d(cl,matrix,vector,klev)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk4), pointer, contiguous, dimension(:,:), intent(inout) :: vector
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(in) :: matrix
    integer(ik4), intent(in) :: klev
    integer(ik4) :: i, j, k, n, iv
    do k = 1, klev
      iv = 1
      do i = ici1, ici2
        do j = jci1, jci2
          do n = 1, nnsg
            if ( cl%sgmask(n,j,i) ) then
              vector(iv,k) = matrix(n,j,i,k)
              iv = iv + 1
            end if
          end do
        end do
      end do
    end do
  end subroutine mypack_real4_subgrid_4d

  subroutine mypack_real8_subgrid_slice(cl,matrix,vector,k)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk8), pointer, contiguous, dimension(:), intent(inout) :: vector
    real(rk8), pointer, contiguous, dimension(:,:,:,:), intent(in) :: matrix
    integer(ik4), intent(in) :: k
    integer(ik4) :: i, j, n, iv
    iv = 1
    do i = ici1, ici2
      do j = jci1, jci2
        do n = 1, nnsg
          if ( cl%sgmask(n,j,i) ) then
            vector(iv) = matrix(n,j,i,k)
            iv = iv + 1
          end if
        end do
      end do
    end do
  end subroutine mypack_real8_subgrid_slice

  subroutine mypack_real4_subgrid_slice(cl,matrix,vector,k)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk4), pointer, contiguous, dimension(:), intent(inout) :: vector
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(in) :: matrix
    integer(ik4), intent(in) :: k
    integer(ik4) :: i, j, n, iv
    iv = 1
    do i = ici1, ici2
      do j = jci1, jci2
        do n = 1, nnsg
          if ( cl%sgmask(n,j,i) ) then
            vector(iv) = matrix(n,j,i,k)
            iv = iv + 1
          end if
        end do
      end do
    end do
  end subroutine mypack_real4_subgrid_slice

  subroutine myunpack_logical_grid(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    logical, pointer, contiguous, dimension(:), intent(in) :: vector
    logical, pointer, contiguous, dimension(:,:), intent(inout) :: matrix
    integer(ik4) :: i, j, iv
    iv = 1
    do i = ici1, ici2
      do j = jci1, jci2
        if ( cl%gmask(j,i) ) then
          matrix(j,i) = vector(iv)
          iv = iv + 1
        end if
      end do
    end do
  end subroutine myunpack_logical_grid

  subroutine myunpack_logical_subgrid(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    logical, pointer, contiguous, dimension(:), intent(in) :: vector
    logical, pointer, contiguous, dimension(:,:,:), intent(inout) :: matrix
    integer(ik4) :: i, j, n, iv
    iv = 1
    do i = ici1, ici2
      do j = jci1, jci2
        do n = 1, nnsg
          if ( cl%sgmask(n,j,i) ) then
            matrix(n,j,i) = vector(iv)
            iv = iv + 1
          end if
        end do
      end do
    end do
  end subroutine myunpack_logical_subgrid

  subroutine myunpack_integer_grid(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    integer(ik4), pointer, contiguous, dimension(:), intent(in) :: vector
    integer(ik4), pointer, contiguous, dimension(:,:), intent(inout) :: matrix
    integer(ik4) :: i, j, iv
    iv = 1
    do i = ici1, ici2
      do j = jci1, jci2
        if ( cl%gmask(j,i) ) then
          matrix(j,i) = vector(iv)
          iv = iv + 1
        end if
      end do
    end do
  end subroutine myunpack_integer_grid

  subroutine myunpack_integer_subgrid(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    integer(ik4), pointer, contiguous, dimension(:), intent(in) :: vector
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(inout) :: matrix
    integer(ik4) :: i, j, n, iv
    iv = 1
    do i = ici1, ici2
      do j = jci1, jci2
        do n = 1, nnsg
          if ( cl%sgmask(n,j,i) ) then
            matrix(n,j,i) = vector(iv)
            iv = iv + 1
          end if
        end do
      end do
    end do
  end subroutine myunpack_integer_subgrid

  subroutine myunpack_real8_grid(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk8), pointer, contiguous, dimension(:), intent(in) :: vector
    real(rk8), pointer, contiguous, dimension(:,:), intent(inout) :: matrix
    integer(ik4) :: i, j, iv
    iv = 1
    do i = ici1, ici2
      do j = jci1, jci2
        if ( cl%gmask(j,i) ) then
          matrix(j,i) = vector(iv)
          iv = iv + 1
        end if
      end do
    end do
  end subroutine myunpack_real8_grid

  subroutine myunpack_real4_grid(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk4), pointer, contiguous, dimension(:), intent(in) :: vector
    real(rk4), pointer, contiguous, dimension(:,:), intent(inout) :: matrix
    integer(ik4) :: i, j, iv
    iv = 1
    do i = ici1, ici2
      do j = jci1, jci2
        if ( cl%gmask(j,i) ) then
          matrix(j,i) = vector(iv)
          iv = iv + 1
        end if
      end do
    end do
  end subroutine myunpack_real4_grid

  subroutine myunpack_real8_real4_grid(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk8), pointer, contiguous, dimension(:), intent(in) :: vector
    real(rk4), pointer, contiguous, dimension(:,:), intent(inout) :: matrix
    integer(ik4) :: i, j, iv
    iv = 1
    do i = ici1, ici2
      do j = jci1, jci2
        if ( cl%gmask(j,i) ) then
          matrix(j,i) = real(vector(iv),rk4)
          iv = iv + 1
        end if
      end do
    end do
  end subroutine myunpack_real8_real4_grid

  subroutine myunpack_real8_subgrid(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk8), pointer, contiguous, dimension(:), intent(in) :: vector
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(inout) :: matrix
    integer(ik4) :: i, j, n, iv
    iv = 1
    do i = ici1, ici2
      do j = jci1, jci2
        do n = 1, nnsg
          if ( cl%sgmask(n,j,i) ) then
            matrix(n,j,i) = vector(iv)
            iv = iv + 1
          end if
        end do
      end do
    end do
  end subroutine myunpack_real8_subgrid

  subroutine myunpack_real4_subgrid(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk4), pointer, contiguous, dimension(:), intent(in) :: vector
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(inout) :: matrix
    integer(ik4) :: i, j, n, iv
    iv = 1
    do i = ici1, ici2
      do j = jci1, jci2
        do n = 1, nnsg
          if ( cl%sgmask(n,j,i) ) then
            matrix(n,j,i) = vector(iv)
            iv = iv + 1
          end if
        end do
      end do
    end do
  end subroutine myunpack_real4_subgrid

  subroutine myunpack_real8_real4_subgrid(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk8), pointer, contiguous, dimension(:), intent(in) :: vector
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(inout) :: matrix
    integer(ik4) :: i, j, n, iv
    iv = 1
    do i = ici1, ici2
      do j = jci1, jci2
        do n = 1, nnsg
          if ( cl%sgmask(n,j,i) ) then
            matrix(n,j,i) = real(vector(iv),rk4)
            iv = iv + 1
          end if
        end do
      end do
    end do
  end subroutine myunpack_real8_real4_subgrid

  subroutine myunpack_real8_subgrid_4d(cl,vector,matrix,klev)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk8), pointer, contiguous, dimension(:,:), intent(in) :: vector
    real(rk8), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: matrix
    integer(ik4), intent(in) :: klev
    integer(ik4) :: i, j, k, n, iv
    do k = 1, klev
      iv = 1
      do i = ici1, ici2
        do j = jci1, jci2
          do n = 1, nnsg
            if ( cl%sgmask(n,j,i) ) then
              matrix(n,j,i,k) = vector(iv,k)
              iv = iv + 1
            end if
          end do
        end do
      end do
    end do
  end subroutine myunpack_real8_subgrid_4d

  subroutine myunpack_real4_subgrid_4d(cl,vector,matrix,klev)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk4), pointer, contiguous, dimension(:,:), intent(in) :: vector
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: matrix
    integer(ik4), intent(in) :: klev
    integer(ik4) :: i, j, k, n, iv
    do k = 1, klev
      iv = 1
      do i = ici1, ici2
        do j = jci1, jci2
          do n = 1, nnsg
            if ( cl%sgmask(n,j,i) ) then
              matrix(n,j,i,k) = vector(iv,k)
              iv = iv + 1
            end if
          end do
        end do
      end do
    end do
  end subroutine myunpack_real4_subgrid_4d

  subroutine myunpack_real8_real4_subgrid_4d(cl,vector,matrix,klev)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk8), pointer, contiguous, dimension(:,:), intent(in) :: vector
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: matrix
    integer(ik4), intent(in) :: klev
    integer(ik4) :: i, j, k, n, iv
    do k = 1, klev
      iv = 1
      do i = ici1, ici2
        do j = jci1, jci2
          do n = 1, nnsg
            if ( cl%sgmask(n,j,i) ) then
              matrix(n,j,i,k) = real(vector(iv,k),rk4)
              iv = iv + 1
            end if
          end do
        end do
      end do
    end do
  end subroutine myunpack_real8_real4_subgrid_4d

  subroutine myunpack_real8_subgrid_slice(cl,vector,matrix,k)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk8), pointer, contiguous, dimension(:), intent(in) :: vector
    real(rk8), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: matrix
    integer(ik4), intent(in) :: k
    integer(ik4) :: i, j, n, iv
    iv = 1
    do i = ici1, ici2
      do j = jci1, jci2
        do n = 1, nnsg
          if ( cl%sgmask(n,j,i) ) then
            matrix(n,j,i,k) = vector(iv)
            iv = iv + 1
          end if
        end do
      end do
    end do
  end subroutine myunpack_real8_subgrid_slice

  subroutine myunpack_real4_subgrid_slice(cl,vector,matrix,k)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk4), pointer, contiguous, dimension(:), intent(in) :: vector
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: matrix
    integer(ik4), intent(in) :: k
    integer(ik4) :: i, j, n, iv
    iv = 1
    do i = ici1, ici2
      do j = jci1, jci2
        do n = 1, nnsg
          if ( cl%sgmask(n,j,i) ) then
            matrix(n,j,i,k) = vector(iv)
            iv = iv + 1
          end if
        end do
      end do
    end do
  end subroutine myunpack_real4_subgrid_slice

  subroutine myunpack_real8_real4_subgrid_slice(cl,vector,matrix,k)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk8), pointer, contiguous, dimension(:), intent(in) :: vector
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: matrix
    integer(ik4), intent(in) :: k
    integer(ik4) :: i, j, n, iv
    iv = 1
    do i = ici1, ici2
      do j = jci1, jci2
        do n = 1, nnsg
          if ( cl%sgmask(n,j,i) ) then
            matrix(n,j,i,k) = real(vector(iv),rk4)
            iv = iv + 1
          end if
        end do
      end do
    end do
  end subroutine myunpack_real8_real4_subgrid_slice

  subroutine mypack_global_logical_grid(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    logical, pointer, contiguous, dimension(:), intent(inout) :: vector
    logical, pointer, contiguous, dimension(:,:), intent(in) :: matrix
    integer(ik4) :: i, j, iv
    iv = 1
    do i = iout1, iout2
      do j = jout1, jout2
        if ( cl%global_gmask(j,i) ) then
          vector(iv) = matrix(j,i)
          iv = iv + 1
        end if
      end do
    end do
  end subroutine mypack_global_logical_grid

  subroutine mypack_global_logical_subgrid(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    logical, pointer, contiguous, dimension(:), intent(inout) :: vector
    logical, pointer, contiguous, dimension(:,:,:), intent(in) :: matrix
    integer(ik4) :: i, j, n, iv
    iv = 1
    do i = iout1, iout2
      do j = jout1, jout2
        do n = 1, nnsg
          if ( cl%global_sgmask(n,j,i) ) then
            vector(iv) = matrix(n,j,i)
            iv = iv + 1
          end if
        end do
      end do
    end do
  end subroutine mypack_global_logical_subgrid

  subroutine mypack_global_integer_grid(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    integer(ik4), pointer, contiguous, dimension(:), intent(inout) :: vector
    integer(ik4), pointer, contiguous, dimension(:,:), intent(in) :: matrix
    integer(ik4) :: i, j, iv
    iv = 1
    do i = iout1, iout2
      do j = jout1, jout2
        if ( cl%global_gmask(j,i) ) then
          vector(iv) = matrix(j,i)
          iv = iv + 1
        end if
      end do
    end do
  end subroutine mypack_global_integer_grid

  subroutine mypack_global_integer_subgrid(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    integer(ik4), pointer, contiguous, dimension(:), intent(inout) :: vector
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(in) :: matrix
    integer(ik4) :: i, j, n, iv
    iv = 1
    do i = iout1, iout2
      do j = jout1, jout2
        do n = 1, nnsg
          if ( cl%global_sgmask(n,j,i) ) then
            vector(iv) = matrix(n,j,i)
            iv = iv + 1
          end if
        end do
      end do
    end do
  end subroutine mypack_global_integer_subgrid

  subroutine mypack_global_real8_grid(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk8), pointer, contiguous, dimension(:), intent(inout) :: vector
    real(rk8), pointer, contiguous, dimension(:,:), intent(in) :: matrix
    integer(ik4) :: i, j, iv
    iv = 1
    do i = iout1, iout2
      do j = jout1, jout2
        if ( cl%global_gmask(j,i) ) then
          vector(iv) = matrix(j,i)
          iv = iv + 1
        end if
      end do
    end do
  end subroutine mypack_global_real8_grid

  subroutine mypack_global_real4_grid(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk4), pointer, contiguous, dimension(:,:), intent(in) :: matrix
    real(rk4), pointer, contiguous, dimension(:), intent(inout) :: vector
    integer(ik4) :: i, j, iv
    iv = 1
    do i = iout1, iout2
      do j = jout1, jout2
        if ( cl%global_gmask(j,i) ) then
          vector(iv) = matrix(j,i)
          iv = iv + 1
        end if
      end do
    end do
  end subroutine mypack_global_real4_grid

  subroutine mypack_global_real4_real8_grid(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk4), pointer, contiguous, dimension(:,:), intent(in) :: matrix
    real(rk8), pointer, contiguous, dimension(:), intent(inout) :: vector
    integer(ik4) :: i, j, iv
    iv = 1
    do i = iout1, iout2
      do j = jout1, jout2
        if ( cl%global_gmask(j,i) ) then
          vector(iv) = real(matrix(j,i),rk8)
          iv = iv + 1
        end if
      end do
    end do
  end subroutine mypack_global_real4_real8_grid

  subroutine mypack_global_real8_subgrid(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk8), pointer, contiguous, dimension(:), intent(inout) :: vector
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(in) :: matrix
    integer(ik4) :: i, j, n, iv
    iv = 1
    do i = iout1, iout2
      do j = jout1, jout2
        do n = 1, nnsg
          if ( cl%global_sgmask(n,j,i) ) then
            vector(iv) = matrix(n,j,i)
            iv = iv + 1
          end if
        end do
      end do
    end do
  end subroutine mypack_global_real8_subgrid

  subroutine mypack_global_real4_subgrid(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk4), pointer, contiguous, dimension(:), intent(inout) :: vector
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(in) :: matrix
    integer(ik4) :: i, j, n, iv
    iv = 1
    do i = iout1, iout2
      do j = jout1, jout2
        do n = 1, nnsg
          if ( cl%global_sgmask(n,j,i) ) then
            vector(iv) = matrix(n,j,i)
            iv = iv + 1
          end if
        end do
      end do
    end do
  end subroutine mypack_global_real4_subgrid

  subroutine mypack_global_real4_real8_subgrid(cl,matrix,vector)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(in) :: matrix
    real(rk8), pointer, contiguous, dimension(:), intent(inout) :: vector
    integer(ik4) :: i, j, n, iv
    iv = 1
    do i = iout1, iout2
      do j = jout1, jout2
        do n = 1, nnsg
          if ( cl%global_sgmask(n,j,i) ) then
            vector(iv) = real(matrix(n,j,i),rk8)
            iv = iv + 1
          end if
        end do
      end do
    end do
  end subroutine mypack_global_real4_real8_subgrid

  subroutine mypack_global_real8_subgrid_4d(cl,matrix,vector,klev)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk8), pointer, contiguous, dimension(:,:), intent(inout) :: vector
    real(rk8), pointer, contiguous, dimension(:,:,:,:), intent(in) :: matrix
    integer(ik4), intent(in) :: klev
    integer(ik4) :: i, j, k, n, iv
    do k = 1, klev
      iv = 1
      do i = iout1, iout2
        do j = jout1, jout2
          do n = 1, nnsg
            if ( cl%global_sgmask(n,j,i) ) then
              vector(iv,k) = matrix(n,j,i,k)
              iv = iv + 1
            end if
          end do
        end do
      end do
    end do
  end subroutine mypack_global_real8_subgrid_4d

  subroutine mypack_global_real4_subgrid_4d(cl,matrix,vector,klev)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk4), pointer, contiguous, dimension(:,:), intent(inout) :: vector
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(in) :: matrix
    integer(ik4), intent(in) :: klev
    integer(ik4) :: i, j, k, n, iv
    do k = 1, klev
      iv = 1
      do i = iout1, iout2
        do j = jout1, jout2
          do n = 1, nnsg
            if ( cl%global_sgmask(n,j,i) ) then
              vector(iv,k) = matrix(n,j,i,k)
              iv = iv + 1
            end if
          end do
        end do
      end do
    end do
  end subroutine mypack_global_real4_subgrid_4d

  subroutine mypack_global_real4_real8_subgrid_4d(cl,matrix,vector,klev)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(in) :: matrix
    real(rk8), pointer, contiguous, dimension(:,:), intent(inout) :: vector
    integer(ik4), intent(in) :: klev
    integer(ik4) :: i, j, k, n, iv
    do k = 1, klev
      iv = 1
      do i = iout1, iout2
        do j = jout1, jout2
          do n = 1, nnsg
            if ( cl%global_sgmask(n,j,i) ) then
              vector(iv,k) = real(matrix(n,j,i,k),rk8)
              iv = iv + 1
            end if
          end do
        end do
      end do
    end do
  end subroutine mypack_global_real4_real8_subgrid_4d

  subroutine mypack_global_real8_subgrid_slice(cl,matrix,vector,k)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk8), pointer, contiguous, dimension(:), intent(inout) :: vector
    real(rk8), pointer, contiguous, dimension(:,:,:,:), intent(in) :: matrix
    integer(ik4), intent(in) :: k
    integer(ik4) :: i, j, n, iv
    iv = 1
    do i = iout1, iout2
      do j = jout1, jout2
        do n = 1, nnsg
          if ( cl%global_sgmask(n,j,i) ) then
            vector(iv) = matrix(n,j,i,k)
            iv = iv + 1
          end if
        end do
      end do
    end do
  end subroutine mypack_global_real8_subgrid_slice

  subroutine mypack_global_real4_subgrid_slice(cl,matrix,vector,k)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk4), pointer, contiguous, dimension(:), intent(inout) :: vector
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(in) :: matrix
    integer(ik4), intent(in) :: k
    integer(ik4) :: i, j, n, iv
    iv = 1
    do i = iout1, iout2
      do j = jout1, jout2
        do n = 1, nnsg
          if ( cl%global_sgmask(n,j,i) ) then
            vector(iv) = matrix(n,j,i,k)
            iv = iv + 1
          end if
        end do
      end do
    end do
  end subroutine mypack_global_real4_subgrid_slice

  subroutine mypack_global_real4_real8_subgrid_slice(cl,matrix,vector,k)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(in) :: matrix
    real(rk8), pointer, contiguous, dimension(:), intent(inout) :: vector
    integer(ik4), intent(in) :: k
    integer(ik4) :: i, j, n, iv
    iv = 1
    do i = iout1, iout2
      do j = jout1, jout2
        do n = 1, nnsg
          if ( cl%global_sgmask(n,j,i) ) then
            vector(iv) = real(matrix(n,j,i,k),rk8)
            iv = iv + 1
          end if
        end do
      end do
    end do
  end subroutine mypack_global_real4_real8_subgrid_slice

  subroutine myunpack_global_logical_grid(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    logical, pointer, contiguous, dimension(:), intent(in) :: vector
    logical, pointer, contiguous, dimension(:,:), intent(inout) :: matrix
    integer(ik4) :: i, j, iv
    iv = 1
    do i = iout1, iout2
      do j = jout1, jout2
        if ( cl%global_gmask(j,i) ) then
          matrix(j,i) = vector(iv)
          iv = iv + 1
        end if
      end do
    end do
  end subroutine myunpack_global_logical_grid

  subroutine myunpack_global_logical_subgrid(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    logical, pointer, contiguous, dimension(:), intent(in) :: vector
    logical, pointer, contiguous, dimension(:,:,:), intent(inout) :: matrix
    integer(ik4) :: i, j, n, iv
    iv = 1
    do i = iout1, iout2
      do j = jout1, jout2
        do n = 1, nnsg
          if ( cl%global_sgmask(n,j,i) ) then
            matrix(n,j,i) = vector(iv)
            iv = iv + 1
          end if
        end do
      end do
    end do
  end subroutine myunpack_global_logical_subgrid

  subroutine myunpack_global_integer_grid(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    integer(ik4), pointer, contiguous, dimension(:), intent(in) :: vector
    integer(ik4), pointer, contiguous, dimension(:,:), intent(inout) :: matrix
    integer(ik4) :: i, j, iv
    iv = 1
    do i = iout1, iout2
      do j = jout1, jout2
        if ( cl%global_gmask(j,i) ) then
          matrix(j,i) = vector(iv)
          iv = iv + 1
        end if
      end do
    end do
  end subroutine myunpack_global_integer_grid

  subroutine myunpack_global_integer_subgrid(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    integer(ik4), pointer, contiguous, dimension(:), intent(in) :: vector
    integer(ik4), pointer, contiguous, dimension(:,:,:), intent(inout) :: matrix
    integer(ik4) :: i, j, n, iv
    iv = 1
    do i = iout1, iout2
      do j = jout1, jout2
        do n = 1, nnsg
          if ( cl%global_sgmask(n,j,i) ) then
            matrix(n,j,i) = vector(iv)
            iv = iv + 1
          end if
        end do
      end do
    end do
  end subroutine myunpack_global_integer_subgrid

  subroutine myunpack_global_real8_grid(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk8), pointer, contiguous, dimension(:), intent(in) :: vector
    real(rk8), pointer, contiguous, dimension(:,:), intent(inout) :: matrix
    integer(ik4) :: i, j, iv
    iv = 1
    do i = iout1, iout2
      do j = jout1, jout2
        if ( cl%global_gmask(j,i) ) then
          matrix(j,i) = vector(iv)
          iv = iv + 1
        end if
      end do
    end do
  end subroutine myunpack_global_real8_grid

  subroutine myunpack_global_real4_grid(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk4), pointer, contiguous, dimension(:), intent(in) :: vector
    real(rk4), pointer, contiguous, dimension(:,:), intent(inout) :: matrix
    integer(ik4) :: i, j, iv
    iv = 1
    do i = iout1, iout2
      do j = jout1, jout2
        if ( cl%global_gmask(j,i) ) then
          matrix(j,i) = vector(iv)
          iv = iv + 1
        end if
      end do
    end do
  end subroutine myunpack_global_real4_grid

  subroutine myunpack_global_real8_real4_grid(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk8), pointer, contiguous, dimension(:), intent(in) :: vector
    real(rk4), pointer, contiguous, dimension(:,:), intent(inout) :: matrix
    integer(ik4) :: i, j, iv
    iv = 1
    do i = iout1, iout2
      do j = jout1, jout2
        if ( cl%global_gmask(j,i) ) then
          matrix(j,i) = real(vector(iv),rk4)
          iv = iv + 1
        end if
      end do
    end do
  end subroutine myunpack_global_real8_real4_grid

  subroutine myunpack_global_real8_subgrid(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk8), pointer, contiguous, dimension(:), intent(in) :: vector
    real(rk8), pointer, contiguous, dimension(:,:,:), intent(inout) :: matrix
    integer(ik4) :: i, j, n, iv
    iv = 1
    do i = iout1, iout2
      do j = jout1, jout2
        do n = 1, nnsg
          if ( cl%global_sgmask(n,j,i) ) then
            matrix(n,j,i) = vector(iv)
            iv = iv + 1
          end if
        end do
      end do
    end do
  end subroutine myunpack_global_real8_subgrid

  subroutine myunpack_global_real4_subgrid(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk4), pointer, contiguous, dimension(:), intent(in) :: vector
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(inout) :: matrix
    integer(ik4) :: i, j, n, iv
    iv = 1
    do i = iout1, iout2
      do j = jout1, jout2
        do n = 1, nnsg
          if ( cl%global_sgmask(n,j,i) ) then
            matrix(n,j,i) = vector(iv)
            iv = iv + 1
          end if
        end do
      end do
    end do
  end subroutine myunpack_global_real4_subgrid

  subroutine myunpack_global_real8_real4_subgrid(cl,vector,matrix)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk8), pointer, contiguous, dimension(:), intent(in) :: vector
    real(rk4), pointer, contiguous, dimension(:,:,:), intent(inout) :: matrix
    integer(ik4) :: i, j, n, iv
    iv = 1
    do i = iout1, iout2
      do j = jout1, jout2
        do n = 1, nnsg
          if ( cl%global_sgmask(n,j,i) ) then
            matrix(n,j,i) = real(vector(iv),rk4)
            iv = iv + 1
          end if
        end do
      end do
    end do
  end subroutine myunpack_global_real8_real4_subgrid

  subroutine myunpack_global_real8_subgrid_4d(cl,vector,matrix,klev)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk8), pointer, contiguous, dimension(:,:), intent(in) :: vector
    real(rk8), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: matrix
    integer(ik4), intent(in) :: klev
    integer(ik4) :: i, j, k, n, iv
    do k = 1, klev
      iv = 1
      do i = iout1, iout2
        do j = jout1, jout2
          do n = 1, nnsg
            if ( cl%global_sgmask(n,j,i) ) then
              matrix(n,j,i,k) = vector(iv,k)
              iv = iv + 1
            end if
          end do
        end do
      end do
    end do
  end subroutine myunpack_global_real8_subgrid_4d

  subroutine myunpack_global_real4_subgrid_4d(cl,vector,matrix,klev)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk4), pointer, contiguous, dimension(:,:), intent(in) :: vector
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: matrix
    integer(ik4), intent(in) :: klev
    integer(ik4) :: i, j, k, n, iv
    do k = 1, klev
      iv = 1
      do i = iout1, iout2
        do j = jout1, jout2
          do n = 1, nnsg
            if ( cl%global_sgmask(n,j,i) ) then
              matrix(n,j,i,k) = vector(iv,k)
              iv = iv + 1
            end if
          end do
        end do
      end do
    end do
  end subroutine myunpack_global_real4_subgrid_4d

  subroutine myunpack_global_real8_real4_subgrid_4d(cl,vector,matrix,klev)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk8), pointer, contiguous, dimension(:,:), intent(in) :: vector
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: matrix
    integer(ik4), intent(in) :: klev
    integer(ik4) :: i, j, k, n, iv
    do k = 1, klev
      iv = 1
      do i = iout1, iout2
        do j = jout1, jout2
          do n = 1, nnsg
            if ( cl%global_sgmask(n,j,i) ) then
              matrix(n,j,i,k) = real(vector(iv,k),rk4)
              iv = iv + 1
            end if
          end do
        end do
      end do
    end do
  end subroutine myunpack_global_real8_real4_subgrid_4d

  subroutine myunpack_global_real8_subgrid_slice(cl,vector,matrix,k)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk8), pointer, contiguous, dimension(:), intent(in) :: vector
    real(rk8), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: matrix
    integer(ik4), intent(in) :: k
    integer(ik4) :: i, j, n, iv
    iv = 1
    do i = iout1, iout2
      do j = jout1, jout2
        do n = 1, nnsg
          if ( cl%global_sgmask(n,j,i) ) then
            matrix(n,j,i,k) = vector(iv)
            iv = iv + 1
          end if
        end do
      end do
    end do
  end subroutine myunpack_global_real8_subgrid_slice

  subroutine myunpack_global_real4_subgrid_slice(cl,vector,matrix,k)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk4), pointer, contiguous, dimension(:), intent(in) :: vector
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: matrix
    integer(ik4), intent(in) :: k
    integer(ik4) :: i, j, n, iv
    iv = 1
    do i = iout1, iout2
      do j = jout1, jout2
        do n = 1, nnsg
          if ( cl%global_sgmask(n,j,i) ) then
            matrix(n,j,i,k) = vector(iv)
            iv = iv + 1
          end if
        end do
      end do
    end do
  end subroutine myunpack_global_real4_subgrid_slice

  subroutine myunpack_global_real8_real4_subgrid_slice(cl,vector,matrix,k)
    implicit none
    type(masked_comm), intent(in) :: cl
    real(rk8), pointer, contiguous, dimension(:), intent(in) :: vector
    real(rk4), pointer, contiguous, dimension(:,:,:,:), intent(inout) :: matrix
    integer(ik4), intent(in) :: k
    integer(ik4) :: i, j, n, iv
    iv = 1
    do i = iout1, iout2
      do j = jout1, jout2
        do n = 1, nnsg
          if ( cl%global_sgmask(n,j,i) ) then
            matrix(n,j,i,k) = real(vector(iv),rk4)
            iv = iv + 1
          end if
        end do
      end do
    end do
  end subroutine myunpack_global_real8_real4_subgrid_slice

  integer(ik4) function get_cartcomm( ) result(cc)
    implicit none
    cc = cartesian_communicator
  end function get_cartcomm

end module mod_mppparam
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

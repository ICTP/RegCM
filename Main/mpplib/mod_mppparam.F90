!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
module mod_mppparam

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_mpmessage
  use mod_memutil
  use mod_date
  use mod_stdio
  use netcdf

  private

  include 'mpif.h'
#ifdef MPI_SERIAL
  integer(ik4) mpi_status_ignore(mpi_status_size)
  integer(ik4) , parameter :: mpi_proc_null = -2
#endif

  public :: set_nproc , broadcast_params

  integer(ik4) :: cartesian_communicator

  integer(ik4) , public :: ncout_mpi_info = mpi_info_null

  type model_area
    logical :: bandflag
    logical :: has_bdy
    logical :: has_bdyleft , has_bdyright , has_bdytop , has_bdybottom
    logical :: has_bdytopleft , has_bdytopright
    logical :: has_bdybottomleft , has_bdybottomright
    integer(ik4) , dimension(2) :: location
    integer(ik4) :: left , right , top , bottom
    integer(ik4) :: topleft , topright , bottomleft , bottomright
    integer(ik4) :: ibt1 , ibt2 , ibb1 , ibb2
    integer(ik4) :: jbl1 , jbl2 , jbr1 , jbr2
  end type model_area

  type grid_nc_var2d
    character(len=64) :: varname
    integer(ik4) :: irec = -1
    integer(ik4) :: ncid = -1
    integer(ik4) :: varid = 1
    integer(ik4) :: nx = 0
    integer(ik4) :: ny = 0
    integer(ik4) :: mynx1 = 0
    integer(ik4) :: mynx2 = 0
    integer(ik4) :: myny1 = 0
    integer(ik4) :: myny2 = 0
    real(rk8) , pointer , dimension(:,:) :: val => null()
    real(rk8) , pointer , dimension(:,:) :: iobuf => null()
  end type grid_nc_var2d

  type grid_nc_var3d
    character(len=64) :: varname
    integer(ik4) :: irec = -1
    integer(ik4) :: ncid = -1
    integer(ik4) :: varid = 1
    integer(ik4) :: nx = 0
    integer(ik4) :: ny = 0
    integer(ik4) :: mynx1 = 0
    integer(ik4) :: mynx2 = 0
    integer(ik4) :: myny1 = 0
    integer(ik4) :: myny2 = 0
    integer(ik4) :: nz = 0
    real(rk8) , pointer , dimension(:,:,:) :: val => null()
    real(rk8) , pointer , dimension(:,:,:) :: iobuf => null()
  end type grid_nc_var3d

  type grid_nc_var4d
    character(len=64) :: varname
    integer(ik4) :: irec = -1
    integer(ik4) :: ncid = -1
    integer(ik4) :: varid = 1
    integer(ik4) :: nx = 0
    integer(ik4) :: ny = 0
    integer(ik4) :: mynx1 = 0
    integer(ik4) :: mynx2 = 0
    integer(ik4) :: myny1 = 0
    integer(ik4) :: myny2 = 0
    integer(ik4) :: nz = 0
    integer(ik4) :: nl = 0
    real(rk8) , pointer , dimension(:,:,:,:) :: val => null()
    real(rk8) , pointer , dimension(:,:,:,:) :: iobuf => null()
  end type grid_nc_var4d

  public :: grid_nc_var2d , grid_nc_var3d , grid_nc_var4d

  interface grid_nc_create
    module procedure grid_nc_create_var2d , &
                     grid_nc_create_var3d , &
                     grid_nc_create_var4d
  end interface grid_nc_create

  interface grid_nc_write
    module procedure grid_nc_write_var2d , &
                     grid_nc_write_var3d , &
                     grid_nc_write_var4d
  end interface grid_nc_write

  interface grid_nc_destroy
    module procedure grid_nc_destroy_var2d , &
                     grid_nc_destroy_var3d , &
                     grid_nc_destroy_var4d
  end interface grid_nc_destroy

  public :: grid_nc_create , grid_nc_write , grid_nc_destroy

  interface grid_distribute
    module procedure real8_2d_distribute ,   &
                     real8_3d_distribute ,   &
                     real8_4d_distribute ,   &
                     real4_2d_distribute ,   &
                     real4_3d_distribute ,   &
                     real4_4d_distribute ,   &
                     integer_2d_distribute , &
                     integer_3d_distribute , &
                     integer_4d_distribute
  end interface grid_distribute

  interface subgrid_distribute
    module procedure real8_2d_sub_distribute ,   &
                     real8_3d_sub_distribute ,   &
                     real4_2d_sub_distribute ,   &
                     real4_3d_sub_distribute ,   &
                     integer_2d_sub_distribute , &
                     integer_3d_sub_distribute
  end interface subgrid_distribute

  interface grid_collect
    module procedure real8_2d_collect ,    &
                     real8_2d_3d_collect , &
                     real8_3d_collect ,    &
                     real8_3d_2d_collect , &
                     real8_4d_collect ,    &
                     real8_4d_2d_collect , &
                     real4_2d_collect ,    &
                     real4_3d_collect ,    &
                     real4_4d_collect ,    &
                     integer_2d_collect ,  &
                     integer_3d_collect ,  &
                     integer_4d_collect
  end interface grid_collect

  interface subgrid_collect
    module procedure real8_2d_sub_collect ,   &
                     real8_3d_sub_collect ,   &
                     real4_2d_sub_collect ,   &
                     real4_3d_sub_collect ,   &
                     integer_2d_sub_collect , &
                     integer_3d_sub_collect
  end interface subgrid_collect

  interface exchange
    module procedure real8_2d_exchange , &
                     real8_3d_exchange , &
                     real8_4d_exchange
  end interface exchange

  interface exchange_lb
    module procedure real8_2d_exchange_left_bottom , &
                     real8_3d_exchange_left_bottom , &
                     real8_4d_exchange_left_bottom
  end interface exchange_lb

  interface exchange_rt
    module procedure real8_2d_exchange_right_top , &
                     real8_3d_exchange_right_top , &
                     real8_4d_exchange_right_top
  end interface exchange_rt

  interface exchange_bdy_lr
    module procedure real8_bdy_exchange_left_right
  end interface exchange_bdy_lr

  interface exchange_bdy_tb
    module procedure real8_bdy_exchange_top_bottom
  end interface exchange_bdy_tb

  interface grid_fill
    module procedure real8_2d_grid_fill
  end interface grid_fill

  interface bcast
    module procedure bcast_logical,       &
                     bcast_int4,          &
                     bcast_int8,          &
                     bcast_real4,         &
                     bcast_real8,         &
                     bcast_arr_logical,   &
                     bcast_arr_character, &
                     bcast_arr_int4,      &
                     bcast_arr_int8,      &
                     bcast_arr_real4,     &
                     bcast_arr_real8,     &
                     bcast_rcm_time_and_date
  end interface bcast

  interface sumall
    module procedure sumall_real8 , &
                     sumall_int4
  end interface sumall

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

  interface reorder_subgrid
    module procedure reorder_subgrid_2d
    module procedure reorder_subgrid_2d3d
    module procedure reorder_subgrid_3d
  end interface reorder_subgrid

  interface reorder_add_subgrid
    module procedure reorder_add_subgrid_2d
    module procedure reorder_add_subgrid_2d3d
    module procedure reorder_add_subgrid_3d
  end interface reorder_add_subgrid

  public :: model_area
  type(model_area) , public :: ma
!
  real(rk8) , pointer , dimension(:) :: r8vector1
  real(rk8) , pointer , dimension(:) :: r8vector2
  real(rk4) , pointer , dimension(:) :: r4vector1
  real(rk4) , pointer , dimension(:) :: r4vector2
  integer(ik4) , pointer , dimension(:) :: i4vector1
  integer(ik4) , pointer , dimension(:) :: i4vector2
  integer(ik4) , dimension(4) :: window
  integer(ik4) :: mpierr
!
  integer(ik4) , public , parameter :: iocpu = 0 ! The id of the cpu doing I/O
  integer(ik4) , public , parameter :: italk = 0 ! Who is doing the print ?
!
  integer(ik4) , parameter :: tag_bt = 1     ! FROM bottom TO top
  integer(ik4) , parameter :: tag_tb = 2     ! FROM top TO bottom
  integer(ik4) , parameter :: tag_lr = 3     ! FROM left TO right
  integer(ik4) , parameter :: tag_rl = 4     ! FROM right TO left
  integer(ik4) , parameter :: tag_brtl = 5   ! FROM bottomrigth TO topleft
  integer(ik4) , parameter :: tag_tlbr = 6   ! FROM topleft TO bottomright
  integer(ik4) , parameter :: tag_bltr = 7   ! FROM bottomleft TO topright
  integer(ik4) , parameter :: tag_trbl = 8   ! FROM topright TO bottomleft
  integer(ik4) , parameter :: tag_w = 100    ! The global indexes from the cpu
  integer(ik4) , parameter :: tag_base = 200 ! The data to/from the cpu to iocpu
!
  public :: exchange , exchange_lb , exchange_rt
  public :: exchange_bdy_lr , exchange_bdy_tb
  public :: grid_distribute , grid_collect , grid_fill
  public :: subgrid_distribute , subgrid_collect
  public :: uvcross2dot , psc2psd
  public :: bcast , sumall
  public :: gather_r , gather_i
  public :: reorder_subgrid , reorder_add_subgrid
  public :: input_reorder
!
  contains
!
#ifdef MPI_SERIAL
  subroutine mpi_sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, &
                          recvbuf, recvcount, recvtype, source, recvtag, &
                          comm, status, ierror)
    implicit none
    real(rk8) , dimension(:) :: sendbuf , recvbuf
    integer(ik4) :: sendcount , sendtype , dest , sendtag
    integer(ik4) :: recvcount , recvtype , source , recvtag , comm
    integer(ik4) :: status(mpi_status_size)
    integer(ik4) :: ierror
  end subroutine mpi_sendrecv

  subroutine mpi_cart_create(comm_old,ndims,dims,periods,reorder, &
                             comm_cart,ierror)
    implicit none
    integer(ik4) :: comm_old , ndims , comm_cart , ierror
    integer(ik4) , dimension(:) :: dims
    logical :: reorder
    logical , dimension(:) :: periods
  end subroutine mpi_cart_create
  
  subroutine mpi_cart_coords(comm,rank,maxdims,coords,ierror)
    implicit none
    integer(ik4) :: comm , rank , maxdims , ierror
    integer(ik4) , dimension(:) :: coords
  end subroutine mpi_cart_coords

  subroutine mpi_cart_rank(comm,coords,rank,ierror)
    implicit none
    integer(ik4) :: comm , rank , ierror
    integer(ik4) , dimension(:) :: coords
  end subroutine mpi_cart_rank
#endif

  subroutine bcast_logical(lval)
    implicit none
    logical , intent(inout) :: lval
    call mpi_bcast(lval,1,mpi_logical,iocpu,mycomm,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
  end subroutine bcast_logical

  subroutine bcast_int4(ival)
    implicit none
    integer(ik4) , intent(inout) :: ival
    call mpi_bcast(ival,1,mpi_integer4,iocpu,mycomm,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
  end subroutine bcast_int4

  subroutine bcast_int8(ival)
    implicit none
    integer(ik8) , intent(inout) :: ival
    call mpi_bcast(ival,1,mpi_integer8,iocpu,mycomm,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
  end subroutine bcast_int8

  subroutine bcast_real4(rval)
    implicit none
    real(rk4) , intent(inout) :: rval
    call mpi_bcast(rval,1,mpi_real4,iocpu,mycomm,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
  end subroutine bcast_real4

  subroutine bcast_real8(rval)
    implicit none
    real(rk8) , intent(inout) :: rval
    call mpi_bcast(rval,1,mpi_real8,iocpu,mycomm,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
  end subroutine bcast_real8

  subroutine bcast_arr_logical(lval)
    implicit none
    logical , dimension(:) , intent(inout) :: lval
    call mpi_bcast(lval,1,mpi_logical,iocpu,mycomm,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
  end subroutine bcast_arr_logical

  subroutine bcast_arr_character(cval,is)
    implicit none
    character(len=*) , intent(inout) :: cval
    integer(ik4) , intent(in) :: is
    call mpi_bcast(cval,is,mpi_character,iocpu,mycomm,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
  end subroutine bcast_arr_character

  subroutine bcast_arr_int4(ival)
    implicit none
    integer(ik4) , dimension(:) , intent(inout) :: ival
    call mpi_bcast(ival,size(ival),mpi_integer4,iocpu,mycomm,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
  end subroutine bcast_arr_int4

  subroutine bcast_arr_int8(ival)
    implicit none
    integer(ik8) , dimension(:) , intent(inout) :: ival
    call mpi_bcast(ival,size(ival),mpi_integer8,iocpu,mycomm,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
  end subroutine bcast_arr_int8

  subroutine bcast_arr_real4(rval)
    implicit none
    real(rk4) , dimension(:) , intent(inout) :: rval
    call mpi_bcast(rval,size(rval),mpi_real4,iocpu,mycomm,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
  end subroutine bcast_arr_real4

  subroutine bcast_arr_real8(rval)
    implicit none
    real(rk8) , dimension(:) , intent(inout) :: rval
    call mpi_bcast(rval,size(rval),mpi_real8,iocpu,mycomm,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
  end subroutine bcast_arr_real8

  subroutine bcast_rcm_time_and_date(x)
    implicit none
    type (rcm_time_and_date) , intent(inout) :: x
    call bcast(x%calendar)
    call bcast(x%days_from_reference)
    call bcast(x%second_of_day)
  end subroutine bcast_rcm_time_and_date

  subroutine sumall_real8(rlval,rtval)
    implicit none
    real(rk8) , intent(in) :: rlval
    real(rk8) , intent(out) :: rtval
    call mpi_reduce(rlval,rtval,1,mpi_real8,mpi_sum,iocpu,mycomm,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
  end subroutine sumall_real8

  subroutine sumall_int4(ilval,itval)
    implicit none
    integer(ik4) , intent(in) :: ilval
    integer(ik4) , intent(out) :: itval
    call mpi_reduce(ilval,itval,1,mpi_integer4,mpi_sum,iocpu,mycomm,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
  end subroutine sumall_int4

  subroutine send_array_logical(lval,isize,icpu,itag)
    implicit none
    logical , dimension(:) , intent(in) :: lval
    integer(ik4) , intent(in) :: isize , icpu , itag
    call mpi_send(lval,isize,mpi_logical,icpu,itag, &
                  cartesian_communicator,mpi_status_ignore,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_send error.')
    end if
  end subroutine send_array_logical

  subroutine send_array_int4(ival,isize,icpu,itag)
    implicit none
    integer(ik4) , dimension(:) , intent(in) :: ival
    integer(ik4) , intent(in) :: isize , icpu , itag
    call mpi_send(ival,isize,mpi_integer4,icpu,itag, &
                  cartesian_communicator,mpi_status_ignore,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_send error.')
    end if
  end subroutine send_array_int4

  subroutine send_array_real4(rval,isize,icpu,itag)
    implicit none
    real(rk4) , dimension(:) , intent(in) :: rval
    integer(ik4) , intent(in) :: isize , icpu , itag
    call mpi_send(rval,isize,mpi_real4,icpu,itag, &
                  cartesian_communicator,mpi_status_ignore,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_send error.')
    end if
  end subroutine send_array_real4

  subroutine send_array_real8(rval,isize,icpu,itag)
    implicit none
    real(rk8) , dimension(:) , intent(in) :: rval
    integer(ik4), intent(in) :: isize , icpu , itag
    call mpi_send(rval,isize,mpi_real8,icpu,itag, &
                  cartesian_communicator,mpi_status_ignore,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_send error.')
    end if
  end subroutine send_array_real8

  subroutine recv_array_logical(lval,isize,icpu,itag)
    implicit none
    logical , dimension(:) , intent(out) :: lval
    integer(ik4) , intent(in) :: isize , icpu , itag
    call mpi_recv(lval,isize,mpi_logical,icpu,itag, &
                  cartesian_communicator,mpi_status_ignore,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_recv error.')
    end if
  end subroutine recv_array_logical

  subroutine recv_array_int4(ival,isize,icpu,itag)
    implicit none
    integer(ik4) , dimension(:) , intent(out) :: ival
    integer(ik4) , intent(in) :: isize , icpu , itag
    call mpi_recv(ival,isize,mpi_integer4,icpu,itag, &
                  cartesian_communicator,mpi_status_ignore,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_recv error.')
    end if
  end subroutine recv_array_int4

  subroutine recv_array_real4(rval,isize,icpu,itag)
    implicit none
    real(rk4) , dimension(:) , intent(out) :: rval
    integer(ik4) , intent(in) :: isize , icpu , itag
    call mpi_recv(rval,isize,mpi_real4,icpu,itag, &
                  cartesian_communicator,mpi_status_ignore,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_recv error.')
    end if
  end subroutine recv_array_real4

  subroutine recv_array_real8(rval,isize,icpu,itag)
    implicit none
    real(rk8) , dimension(:) , intent(out) :: rval
    integer(ik4), intent(in) :: isize , icpu , itag
    call mpi_recv(rval,isize,mpi_real8,icpu,itag, &
                  cartesian_communicator,mpi_status_ignore,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_recv error.')
    end if
  end subroutine recv_array_real8

  subroutine exchange_array(isize,icpu,tag1,tag2)
    implicit none
    integer(ik4) , intent(in) :: isize , icpu , tag1 , tag2
    integer(ik4) :: ireq
    call mpi_irecv(r8vector2,isize,mpi_real8,icpu,tag1, &
                   cartesian_communicator,ireq,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_irecv error.')
    end if
    call mpi_send(r8vector1,isize,mpi_real8,icpu,tag2, &
                  cartesian_communicator,mpi_status_ignore,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_send error.')
    end if
    call mpi_wait(ireq,mpi_status_ignore,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_wait error.')
    end if
  end subroutine exchange_array

  subroutine cyclic_exchange_array(isize,icpu1,icpu2,itag)
    implicit none
    integer(ik4) , intent(in) :: isize , icpu1 , icpu2 , itag
    call mpi_sendrecv(r8vector1,isize,mpi_real8,icpu1,itag, &
                      r8vector2,isize,mpi_real8,icpu2,itag, &
                      cartesian_communicator,mpi_status_ignore,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_sendrecv error.')
    end if
  end subroutine cyclic_exchange_array

  subroutine set_nproc
    implicit none
    integer(ik4) , dimension(2) :: cpus_per_dim
    logical , dimension(2) :: dim_period
    integer(ik4) , dimension(2) :: isearch
    integer(ik4) :: imaxcpus , imax1 , imax2
    real(rk8) :: dimfac
    data dim_period /.false.,.false./

    ma%bandflag    = (i_band == 1)
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
      cartesian_communicator = mycomm
      ma%location(1) = 0
      ma%location(2) = 0
      global_dot_jstart = 1
      global_dot_istart = 1
      global_dot_jend = jx
      global_dot_iend = iy
      global_cross_jstart = 1
      global_cross_istart = 1
      global_cross_jend = jx-1
      global_cross_iend = iy-1
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
    else
      if ( ma%bandflag ) dim_period(1) = .true.
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
          if ( mod(nproc,2) /= 0 ) then
            write(stderr,*) 'Work does not evenly divide.'
            write(stderr,*) 'Required number of CPUS must be even.'
            call fatal(__FILE__,__LINE__,'CPU/WORK mismatch')
          end if
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
            write(stderr,*) 'Suggested ratio per dimension : ', dimfac
            write(stderr,*) 'Closest number : ' , imaxcpus
            call fatal(__FILE__,__LINE__,'CPU/WORK mismatch')
          end if
        end if
      end if

      call mpi_cart_create(mycomm,2,cpus_per_dim,dim_period,.true., &
                           cartesian_communicator,mpierr)
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_cart_create error.')
      end if
      call mpi_comm_rank(cartesian_communicator,myid,mpierr)
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_comm_rank error.')
      end if
      call mpi_cart_coords(cartesian_communicator,myid,2,ma%location,mpierr)
      if ( mpierr /= mpi_success ) then
        call fatal(__FILE__,__LINE__,'mpi_cart_coords error.')
      end if

      ! Set coordinates in the grid for the other processors
      isearch(1) = ma%location(1)
      isearch(2) = ma%location(2)+1
      if ( isearch(2) < cpus_per_dim(2) ) then
        call mpi_cart_rank(cartesian_communicator,isearch,ma%top,mpierr)
        if ( mpierr /= mpi_success ) then
          call fatal(__FILE__,__LINE__,'mpi_cart_rank error.')
        end if
      end if
      isearch(1) = ma%location(1)
      isearch(2) = ma%location(2)-1
      if ( isearch(2) >= 0 ) then
        call mpi_cart_rank(cartesian_communicator,isearch,ma%bottom,mpierr)
        if ( mpierr /= mpi_success ) then
          call fatal(__FILE__,__LINE__,'mpi_cart_rank error.')
        end if
      end if
      isearch(1) = ma%location(1)-1
      isearch(2) = ma%location(2)
      if ( ma%bandflag .or. ( isearch(1) >= 0 ) ) then
        call mpi_cart_rank(cartesian_communicator,isearch,ma%left,mpierr)
        if ( mpierr /= mpi_success ) then
          call fatal(__FILE__,__LINE__,'mpi_cart_rank error.')
        end if
      end if
      isearch(1) = ma%location(1)+1
      isearch(2) = ma%location(2)
      if ( ma%bandflag .or. ( isearch(1) < cpus_per_dim(1) ) ) then
        call mpi_cart_rank(cartesian_communicator,isearch,ma%right,mpierr)
        if ( mpierr /= mpi_success ) then
          call fatal(__FILE__,__LINE__,'mpi_cart_rank error.')
        end if
      end if
      isearch(1) = ma%location(1)+1
      isearch(2) = ma%location(2)+1
      if ( ( ma%bandflag .or. ( isearch(1) < cpus_per_dim(1) ) ) .and. &
             isearch(2) < cpus_per_dim(2) ) then
        call mpi_cart_rank(cartesian_communicator,isearch,ma%topright,mpierr)
        if ( mpierr /= mpi_success ) then
          call fatal(__FILE__,__LINE__,'mpi_cart_rank error.')
        end if
      end if
      isearch(1) = ma%location(1)-1
      isearch(2) = ma%location(2)+1
      if ( ( ma%bandflag .or. ( isearch(1) >= 0 ) ) .and. &
             isearch(2) < cpus_per_dim(2) ) then
        call mpi_cart_rank(cartesian_communicator,isearch,ma%topleft,mpierr)
        if ( mpierr /= mpi_success ) then
          call fatal(__FILE__,__LINE__,'mpi_cart_rank error.')
        end if
      end if
      isearch(1) = ma%location(1)+1
      isearch(2) = ma%location(2)-1
      if ( ( ma%bandflag .or. ( isearch(1) < cpus_per_dim(1) ) ) .and. &
             isearch(2) >= 0 ) then
        call mpi_cart_rank(cartesian_communicator,isearch,ma%bottomright,mpierr)
        if ( mpierr /= mpi_success ) then
          call fatal(__FILE__,__LINE__,'mpi_cart_rank error.')
        end if
      end if
      isearch(1) = ma%location(1)-1
      isearch(2) = ma%location(2)-1
      if ( ( ma%bandflag .or. ( isearch(1) >= 0 ) ) .and. isearch(2) >= 0 ) then
        call mpi_cart_rank(cartesian_communicator,isearch,ma%bottomleft,mpierr)
        if ( mpierr /= mpi_success ) then
          call fatal(__FILE__,__LINE__,'mpi_cart_rank error.')
        end if
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
      !
      ! Topmost and rightmost processors are doing what's left
      !
      if ( ma%location(2)+1 == cpus_per_dim(2) ) then
        if ( mod(iy,iyp) /= 0 ) then
          global_dot_iend = global_dot_istart+iyp-1
          if ( global_dot_iend < iy ) then
            iyp = iy-global_dot_istart+1
          else
            iyp = mod(iy,iyp)
          end if
        end if
      end if
      if ( ma%location(1)+1 == cpus_per_dim(1) ) then
        if ( mod(jx,jxp) /= 0 ) then
          if ( global_dot_jend < jx ) then
            jxp = jx-global_dot_jstart+1
          else
            jxp = mod(jx,jxp)
          end if
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
      global_cross_jstart = global_dot_jstart
      global_cross_istart = global_dot_istart
      global_cross_jend   = global_dot_jend
      global_cross_iend   = global_dot_iend
      if ( global_dot_jend == jx ) global_cross_jend = global_cross_jend - 1
      if ( global_dot_iend == iy ) global_cross_iend = global_cross_iend - 1
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
    ! Allocate to something should fit all
    if ( nproc > 1 ) then
      call getmem1d(r8vector1,1,nsg*jx*iy*kz,'set_nproc:r8vector1')
      call getmem1d(r8vector2,1,nsg*jx*iy*kz,'set_nproc:r8vector2')
      call getmem1d(r4vector1,1,nsg*jx*iy*kz,'set_nproc:r4vector1')
      call getmem1d(r4vector2,1,nsg*jx*iy*kz,'set_nproc:r4vector2')
      call getmem1d(i4vector1,1,nsg*jx*iy*kz,'set_nproc:i4vector1')
      call getmem1d(i4vector2,1,nsg*jx*iy*kz,'set_nproc:i4vector2')
    end if
  end subroutine set_nproc

  subroutine broadcast_params
    implicit none

    call bcast(iy)
    call bcast(jx)
    call bcast(kz)
    call bcast(nsg)
    call bcast(njxcpus)
    call bcast(niycpus)
    call bcast(nveg)
    call bcast(iproj,6)
    call bcast(ds)
    call bcast(ptop)
    call bcast(clat)
    call bcast(clon)
    call bcast(plat)
    call bcast(plon)
    call bcast(truelatl)
    call bcast(truelath)
    call bcast(i_band)

    call bcast(domname,64)
    call bcast(ibyte)

    call bcast(debug_level)
    call bcast(dbgfrq)

    call bcast(nspgx)
    call bcast(nspgd)
    call bcast(high_nudge)
    call bcast(medium_nudge)
    call bcast(low_nudge)

    call bcast(calendar,12)
    call bcast(ical)
    call bcast(dayspy)
    call bcast(dpd)

    call bcast(nsplit)

    call bcast(ibdyfrq)

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
      else
        jcross2 = jxm1
        jout1 = 2
        jout2 = jxm2
      end if
      idot1 = 1
      idot2 = iy
      icross1 = 1
      icross2 = iym1
      iout1 = 2
      iout2 = iym2
      njcross = jcross2-jcross1+1
      nicross = icross2-icross1+1
      njdot = jdot2-jdot1+1
      nidot = idot2-idot1+1
      njout = jout2-jout1+1
      niout = iout2-iout1+1
    end if
  end subroutine broadcast_params
!
  subroutine real8_2d_distribute(mg,ml,j1,j2,i1,i2)
    implicit none
    real(rk8) , pointer , dimension(:,:) , intent(in) :: mg  ! model global
    real(rk8) , pointer , dimension(:,:) , intent(out) :: ml ! model local
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2
    integer(ik4) :: ib , i , j , isize , jsize , lsize , icpu
    if ( myid == iocpu ) then
      ! Copy in memory my piece.
      do i = i1 , i2
        do j = j1 , j2
          ml(j,i) = mg(global_dot_jstart+j-1,global_dot_istart+i-1)
        end do
      end do
      ! Send to other nodes the piece they request.
      do icpu = 1 , nproc-1
        call recv_array(window,4,icpu,tag_w)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        lsize = isize*jsize
        if ( size(r8vector1) < lsize ) then
          call getmem1d(r8vector1,1,lsize,'real8_2d_distribute')
        end if
        ib = 1
        do i = window(1) , window(2)
          do j = window(3) , window(4)
            r8vector1(ib) = mg(j,i)
            ib = ib + 1
          end do
        end do
        call send_array(r8vector1,lsize,icpu,tag_base)
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      lsize = isize*jsize
      if ( size(r8vector2) < lsize ) then
        call getmem1d(r8vector2,1,lsize,'real8_2d_distribute')
      end if
      window(1) = global_dot_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_dot_jstart+j1-1
      window(4) = window(3)+jsize-1
      call send_array(window,4,iocpu,tag_w)
      call recv_array(r8vector2,lsize,iocpu,tag_base)
      ib = 1
      do i = i1 , i2
        do j = j1 , j2
          ml(j,i) = r8vector2(ib)
          ib = ib + 1
        end do
      end do
    end if
  end subroutine real8_2d_distribute
!
  subroutine real8_3d_distribute(mg,ml,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:) , intent(in) :: mg  ! model global
    real(rk8) , pointer , dimension(:,:,:) , intent(out) :: ml ! model local
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer(ik4) :: ib , i , j , k , isize , jsize , ksize , lsize , icpu
    if ( myid == iocpu ) then
      ! Copy in memory my piece.
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            ml(j,i,k) = mg(global_dot_jstart+j-1,global_dot_istart+i-1,k)
          end do
        end do
      end do
      ! Send to other nodes the piece they request.
      do icpu = 1 , nproc-1
        call recv_array(window,4,icpu,tag_w)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        ksize = k2-k1+1
        lsize = isize*jsize*ksize
        if ( size(r8vector1) < lsize ) then
          call getmem1d(r8vector1,1,lsize,'real8_3d_distribute')
        end if
        ib = 1
        do k = k1 , k2
          do i = window(1) , window(2)
            do j = window(3) , window(4)
              r8vector1(ib) = mg(j,i,k)
              ib = ib + 1
            end do
          end do
        end do
        call send_array(r8vector1,lsize,icpu,tag_base)
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      ksize = k2-k1+1
      lsize = isize*jsize*ksize
      if ( size(r8vector2) < lsize ) then
        call getmem1d(r8vector2,1,lsize,'real8_3d_distribute')
      end if
      window(1) = global_dot_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_dot_jstart+j1-1
      window(4) = window(3)+jsize-1
      call send_array(window,4,iocpu,tag_w)
      call recv_array(r8vector2,lsize,iocpu,tag_base)
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            ml(j,i,k) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine real8_3d_distribute
!
  subroutine real8_4d_distribute(mg,ml,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:,:) , intent(in) :: mg  ! model global
    real(rk8) , pointer , dimension(:,:,:,:) , intent(out) :: ml ! model local
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2 , n1 , n2
    integer(ik4) :: ib , i , j , k , n , isize , jsize , &
                    ksize , nsize , lsize , icpu
    if ( myid == iocpu ) then
      ! Copy in memory my piece.
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = j1 , j2
              ml(j,i,k,n) = mg(global_dot_jstart+j-1,global_dot_istart+i-1,k,n)
            end do
          end do
        end do
      end do
      ! Send to other nodes the piece they request.
      do icpu = 1 , nproc-1
        call recv_array(window,4,icpu,tag_w)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        ksize = k2-k1+1
        nsize = n2-n1+1
        lsize = isize*jsize*ksize*nsize
        if ( size(r8vector1) < lsize ) then
          call getmem1d(r8vector1,1,lsize,'real8_4d_distribute')
        end if
        ib = 1
        do n = n1 , n2
          do k = k1 , k2
            do i = window(1) , window(2)
              do j = window(3) , window(4)
                r8vector1(ib) = mg(j,i,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        call send_array(r8vector1,lsize,icpu,tag_base)
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      ksize = k2-k1+1
      nsize = n2-n1+1
      lsize = isize*jsize*ksize*nsize
      if ( size(r8vector2) < lsize ) then
        call getmem1d(r8vector2,1,lsize,'real8_4d_distribute')
      end if
      window(1) = global_dot_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_dot_jstart+j1-1
      window(4) = window(3)+jsize-1
      call send_array(window,4,iocpu,tag_w)
      call recv_array(r8vector2,lsize,iocpu,tag_base)
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = j1 , j2
              ml(j,i,k,n) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
  end subroutine real8_4d_distribute
!
  subroutine real4_2d_distribute(mg,ml,j1,j2,i1,i2)
    implicit none
    real(rk4) , pointer , dimension(:,:) , intent(in) :: mg  ! model global
    real(rk4) , pointer , dimension(:,:) , intent(out) :: ml ! model local
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2
    integer(ik4) :: ib , i , j , isize , jsize , lsize , icpu
    if ( myid == iocpu ) then
      ! Copy in memory my piece.
      do i = i1 , i2
        do j = j1 , j2
          ml(j,i) = mg(global_dot_jstart+j-1,global_dot_istart+i-1)
        end do
      end do
      ! Send to other nodes the piece they request.
      do icpu = 1 , nproc-1
        call recv_array(window,4,icpu,tag_w)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        lsize = isize*jsize
        if ( size(r4vector1) < lsize ) then
          call getmem1d(r4vector1,1,lsize,'real4_2d_distribute')
        end if
        ib = 1
        do i = window(1) , window(2)
          do j = window(3) , window(4)
            r4vector1(ib) = mg(j,i)
            ib = ib + 1
          end do
        end do
        call send_array(r4vector1,lsize,icpu,tag_base)
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      lsize = isize*jsize
      if ( size(r4vector2) < lsize ) then
        call getmem1d(r4vector2,1,lsize,'real4_2d_distribute')
      end if
      window(1) = global_dot_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_dot_jstart+j1-1
      window(4) = window(3)+jsize-1
      call send_array(window,4,iocpu,tag_w)
      call recv_array(r4vector2,lsize,iocpu,tag_base)
      ib = 1
      do i = i1 , i2
        do j = j1 , j2
          ml(j,i) = r4vector2(ib)
          ib = ib + 1
        end do
      end do
    end if
  end subroutine real4_2d_distribute
!
  subroutine real4_3d_distribute(mg,ml,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk4) , pointer , dimension(:,:,:) , intent(in) :: mg  ! model global
    real(rk4) , pointer , dimension(:,:,:) , intent(out) :: ml ! model local
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer(ik4) :: ib , i , j , k , isize , jsize , ksize , lsize , icpu
    if ( myid == iocpu ) then
      ! Copy in memory my piece.
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            ml(j,i,k) = mg(global_dot_jstart+j-1,global_dot_istart+i-1,k)
          end do
        end do
      end do
      ! Send to other nodes the piece they request.
      do icpu = 1 , nproc-1
        call recv_array(window,4,icpu,tag_w)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        ksize = k2-k1+1
        lsize = isize*jsize*ksize
        if ( size(r4vector1) < lsize ) then
          call getmem1d(r4vector1,1,lsize,'real4_3d_distribute')
        end if
        ib = 1
        do k = k1 , k2
          do i = window(1) , window(2)
            do j = window(3) , window(4)
              r4vector1(ib) = mg(j,i,k)
              ib = ib + 1
            end do
          end do
        end do
        call send_array(r4vector1,lsize,icpu,tag_base)
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      ksize = k2-k1+1
      lsize = isize*jsize*ksize
      if ( size(r4vector2) < lsize ) then
        call getmem1d(r4vector2,1,lsize,'real4_3d_distribute')
      end if
      window(1) = global_dot_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_dot_jstart+j1-1
      window(4) = window(3)+jsize-1
      call send_array(window,4,iocpu,tag_w)
      call recv_array(r4vector2,lsize,iocpu,tag_base)
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            ml(j,i,k) = r4vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine real4_3d_distribute
!
  subroutine real4_4d_distribute(mg,ml,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk4) , pointer , dimension(:,:,:,:) , intent(in) :: mg  ! model global
    real(rk4) , pointer , dimension(:,:,:,:) , intent(out) :: ml ! model local
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2 , n1 , n2
    integer(ik4) :: ib , i , j , k , n , isize , &
                    jsize , ksize , nsize , lsize , icpu
    if ( myid == iocpu ) then
      ! Copy in memory my piece.
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = j1 , j2
              ml(j,i,k,n) = mg(global_dot_jstart+j-1,global_dot_istart+i-1,k,n)
            end do
          end do
        end do
      end do
      ! Send to other nodes the piece they request.
      do icpu = 1 , nproc-1
        call recv_array(window,4,icpu,tag_w)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        ksize = k2-k1+1
        nsize = n2-n1+1
        lsize = isize*jsize*ksize*nsize
        if ( size(r4vector1) < lsize ) then
          call getmem1d(r4vector1,1,lsize,'real4_4d_distribute')
        end if
        ib = 1
        do n = n1 , n2
          do k = k1 , k2
            do i = window(1) , window(2)
              do j = window(3) , window(4)
                r4vector1(ib) = mg(j,i,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        call send_array(r4vector1,lsize,icpu,tag_base)
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      ksize = k2-k1+1
      nsize = n2-n1+1
      lsize = isize*jsize*ksize*nsize
      if ( size(r4vector2) < lsize ) then
        call getmem1d(r4vector2,1,lsize,'real4_4d_distribute')
      end if
      window(1) = global_dot_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_dot_jstart+j1-1
      window(4) = window(3)+jsize-1
      call send_array(window,4,iocpu,tag_w)
      call recv_array(r4vector2,lsize,iocpu,tag_base)
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = j1 , j2
              ml(j,i,k,n) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
  end subroutine real4_4d_distribute
!
  subroutine integer_2d_distribute(mg,ml,j1,j2,i1,i2)
    implicit none
    integer(ik4) , pointer , dimension(:,:) , intent(in) :: mg  ! model global
    integer(ik4) , pointer , dimension(:,:) , intent(out) :: ml ! model local
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2
    integer(ik4) :: ib , i , j , isize , jsize , lsize , icpu
    if ( myid == iocpu ) then
      ! Copy in memory my piece.
      do i = i1 , i2
        do j = j1 , j2
          ml(j,i) = mg(global_dot_jstart+j-1,global_dot_istart+i-1)
        end do
      end do
      ! Send to other nodes the piece they request.
      do icpu = 1 , nproc-1
        call recv_array(window,4,icpu,tag_w)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        lsize = isize*jsize
        if ( size(i4vector1) < lsize ) then
          call getmem1d(i4vector1,1,lsize,'integer_2d_distribute')
        end if
        ib = 1
        do i = window(1) , window(2)
          do j = window(3) , window(4)
            i4vector1(ib) = mg(j,i)
            ib = ib + 1
          end do
        end do
        call send_array(i4vector1,lsize,icpu,tag_base)
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      lsize = isize*jsize
      if ( size(i4vector2) < lsize ) then
        call getmem1d(i4vector2,1,lsize,'integer_2d_distribute')
      end if
      window(1) = global_dot_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_dot_jstart+j1-1
      window(4) = window(3)+jsize-1
      call send_array(window,4,iocpu,tag_w)
      call recv_array(i4vector2,lsize,iocpu,tag_base)
      ib = 1
      do i = i1 , i2
        do j = j1 , j2
          ml(j,i) = i4vector2(ib)
          ib = ib + 1
        end do
      end do
    end if
  end subroutine integer_2d_distribute
!
  subroutine integer_3d_distribute(mg,ml,j1,j2,i1,i2,k1,k2)
    implicit none
    integer(ik4) , pointer , dimension(:,:,:) , intent(in) :: mg  ! model global
    integer(ik4) , pointer , dimension(:,:,:) , intent(out) :: ml ! model local
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer(ik4) :: ib , i , j , k , isize , jsize , ksize , lsize , icpu
    if ( myid == iocpu ) then
      ! Copy in memory my piece.
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            ml(j,i,k) = mg(global_dot_jstart+j-1,global_dot_istart+i-1,k)
          end do
        end do
      end do
      ! Send to other nodes the piece they request.
      do icpu = 1 , nproc-1
        call recv_array(window,4,icpu,tag_w)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        ksize = k2-k1+1
        lsize = isize*jsize*ksize
        if ( size(i4vector1) < lsize ) then
          call getmem1d(i4vector1,1,lsize,'integer_3d_distribute')
        end if
        ib = 1
        do k = k1 , k2
          do i = window(1) , window(2)
            do j = window(3) , window(4)
              i4vector1(ib) = mg(j,i,k)
              ib = ib + 1
            end do
          end do
        end do
        call send_array(i4vector1,lsize,icpu,tag_base)
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      ksize = k2-k1+1
      lsize = isize*jsize*ksize
      if ( size(i4vector2) < lsize ) then
        call getmem1d(i4vector2,1,lsize,'integer_3d_distribute')
      end if
      window(1) = global_dot_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_dot_jstart+j1-1
      window(4) = window(3)+jsize-1
      call send_array(window,4,iocpu,tag_w)
      call recv_array(i4vector2,lsize,iocpu,tag_base)
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            ml(j,i,k) = i4vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine integer_3d_distribute
!
  subroutine integer_4d_distribute(mg,ml,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    integer(ik4) , pointer , dimension(:,:,:,:) , intent(in) :: mg  ! model glob
    integer(ik4) , pointer , dimension(:,:,:,:) , intent(out) :: ml ! model loc
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2 , n1 , n2
    integer(ik4) :: ib , i , j , k , n , isize , jsize , &
                    ksize , nsize , lsize , icpu
    if ( myid == iocpu ) then
      ! Copy in memory my piece.
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = j1 , j2
              ml(j,i,k,n) = mg(global_dot_jstart+j-1,global_dot_istart+i-1,k,n)
            end do
          end do
        end do
      end do
      ! Send to other nodes the piece they request.
      do icpu = 1 , nproc-1
        call recv_array(window,4,icpu,tag_w)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        ksize = k2-k1+1
        nsize = n2-n1+1
        lsize = isize*jsize*ksize*nsize
        if ( size(i4vector1) < lsize ) then
          call getmem1d(i4vector1,1,lsize,'integer_4d_distribute')
        end if
        ib = 1
        do n = n1 , n2
          do k = k1 , k2
            do i = window(1) , window(2)
              do j = window(3) , window(4)
                i4vector1(ib) = mg(j,i,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        call send_array(i4vector1,lsize,icpu,tag_base)
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      ksize = k2-k1+1
      nsize = n2-n1+1
      lsize = isize*jsize*ksize*nsize
      if ( size(i4vector2) < lsize ) then
        call getmem1d(i4vector2,1,lsize,'integer_4d_distribute')
      end if
      window(1) = global_dot_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_dot_jstart+j1-1
      window(4) = window(3)+jsize-1
      call send_array(window,4,iocpu,tag_w)
      call recv_array(i4vector2,lsize,iocpu,tag_base)
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = j1 , j2
              ml(j,i,k,n) = i4vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
  end subroutine integer_4d_distribute
!
  subroutine real8_2d_sub_distribute(mg,ml,j1,j2,i1,i2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:) , intent(in) :: mg  ! model global
    real(rk8) , pointer , dimension(:,:,:) , intent(out) :: ml ! model local
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2
    integer(ik4) :: ib , i , j , n , isize , jsize , lsize , icpu
    if ( myid == iocpu ) then
      ! Copy in memory my piece.
      do i = i1 , i2
        do j = j1 , j2
          do n = 1 , nnsg
            ml(n,j,i) = mg(n,global_dot_jstart+j-1,global_dot_istart+i-1)
          end do
        end do
      end do
      ! Send to other nodes the piece they request.
      do icpu = 1 , nproc-1
        call recv_array(window,4,icpu,tag_w)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        lsize = isize*jsize*nnsg
        if ( size(r8vector1) < lsize ) then
          call getmem1d(r8vector1,1,lsize,'real8_2d_sub_distribute')
        end if
        ib = 1
        do i = window(1) , window(2)
          do j = window(3) , window(4)
            do n = 1 , nnsg
              r8vector1(ib) = mg(n,j,i)
              ib = ib + 1
            end do
          end do
        end do
        call send_array(r8vector1,lsize,icpu,tag_base)
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      lsize = isize*jsize*nnsg
      if ( size(r8vector2) < lsize ) then
        call getmem1d(r8vector2,1,lsize,'real8_2d_sub_distribute')
      end if
      window(1) = global_dot_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_dot_jstart+j1-1
      window(4) = window(3)+jsize-1
      call send_array(window,4,iocpu,tag_w)
      call recv_array(r8vector2,lsize,iocpu,tag_base)
      ib = 1
      do i = i1 , i2
        do j = j1 , j2
          do n = 1 , nnsg
            ml(n,j,i) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine real8_2d_sub_distribute
!
  subroutine real8_3d_sub_distribute(mg,ml,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:,:) , intent(in) :: mg  ! model global
    real(rk8) , pointer , dimension(:,:,:,:) , intent(out) :: ml ! model local
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer(ik4) :: ib , i , j , k , n , isize , jsize , ksize , lsize , icpu
    if ( myid == iocpu ) then
      ! Copy in memory my piece.
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            do n = 1 , nnsg
              ml(n,j,i,k) = mg(n,global_dot_jstart+j-1,global_dot_istart+i-1,k)
            end do
          end do
        end do
      end do
      ! Send to other nodes the piece they request.
      do icpu = 1 , nproc-1
        call recv_array(window,4,icpu,tag_w)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        ksize = k2-k1+1
        lsize = isize*jsize*ksize*nnsg
        if ( size(r8vector1) < lsize ) then
          call getmem1d(r8vector1,1,lsize,'real8_3d_sub_distribute')
        end if
        ib = 1
        do k = k1 , k2
          do i = window(1) , window(2)
            do j = window(3) , window(4)
              do n = 1 , nnsg
                r8vector1(ib) = mg(n,j,i,k)
                ib = ib + 1
              end do
            end do
          end do
        end do
        call send_array(r8vector1,lsize,icpu,tag_base)
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      ksize = k2-k1+1
      lsize = isize*jsize*ksize*nnsg
      if ( size(r8vector2) < lsize ) then
        call getmem1d(r8vector2,1,lsize,'real8_3d_sub_distribute')
      end if
      window(1) = global_dot_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_dot_jstart+j1-1
      window(4) = window(3)+jsize-1
      call send_array(window,4,iocpu,tag_w)
      call recv_array(r8vector2,lsize,iocpu,tag_base)
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            do n = 1 , nnsg
              ml(n,j,i,k) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
  end subroutine real8_3d_sub_distribute
!
  subroutine real4_2d_sub_distribute(mg,ml,j1,j2,i1,i2)
    implicit none
    real(rk4) , pointer , dimension(:,:,:) , intent(in) :: mg  ! model global
    real(rk4) , pointer , dimension(:,:,:) , intent(out) :: ml ! model local
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2
    integer(ik4) :: ib , i , j , n , isize , jsize , lsize , icpu
    if ( myid == iocpu ) then
      ! Copy in memory my piece.
      do i = i1 , i2
        do j = j1 , j2
          do n = 1 , nnsg
            ml(n,j,i) = mg(n,global_dot_jstart+j-1,global_dot_istart+i-1)
          end do
        end do
      end do
      ! Send to other nodes the piece they request.
      do icpu = 1 , nproc-1
        call recv_array(window,4,icpu,tag_w)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        lsize = isize*jsize*nnsg
        if ( size(r4vector1) < lsize ) then
          call getmem1d(r4vector1,1,lsize,'real4_2d_sub_distribute')
        end if
        ib = 1
        do i = window(1) , window(2)
          do j = window(3) , window(4)
            do n = 1 , nnsg
              r4vector1(ib) = mg(n,j,i)
              ib = ib + 1
            end do
          end do
        end do
        call send_array(r4vector1,lsize,icpu,tag_base)
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      lsize = isize*jsize*nnsg
      if ( size(r4vector2) < lsize ) then
        call getmem1d(r4vector2,1,lsize,'real4_2d_sub_distribute')
      end if
      window(1) = global_dot_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_dot_jstart+j1-1
      window(4) = window(3)+jsize-1
      call send_array(window,4,iocpu,tag_w)
      call recv_array(r4vector2,lsize,iocpu,tag_base)
      ib = 1
      do i = i1 , i2
        do j = j1 , j2
          do n = 1 , nnsg
            ml(n,j,i) = r4vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine real4_2d_sub_distribute
!
  subroutine real4_3d_sub_distribute(mg,ml,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk4) , pointer , dimension(:,:,:,:) , intent(in) :: mg  ! model global
    real(rk4) , pointer , dimension(:,:,:,:) , intent(out) :: ml ! model local
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer(ik4) :: ib , i , j , k , n , isize , jsize , ksize , lsize , icpu
    if ( myid == iocpu ) then
      ! Copy in memory my piece.
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            do n = 1 , nnsg
              ml(n,j,i,k) = mg(n,global_dot_jstart+j-1,global_dot_istart+i-1,k)
            end do
          end do
        end do
      end do
      ! Send to other nodes the piece they request.
      do icpu = 1 , nproc-1
        call recv_array(window,4,icpu,tag_w)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        ksize = k2-k1+1
        lsize = isize*jsize*ksize*nnsg
        if ( size(r4vector1) < lsize ) then
          call getmem1d(r4vector1,1,lsize,'real4_3d_sub_distribute')
        end if
        ib = 1
        do k = k1 , k2
          do i = window(1) , window(2)
            do j = window(3) , window(4)
              do n = 1 , nnsg
                r4vector1(ib) = mg(n,j,i,k)
                ib = ib + 1
              end do
            end do
          end do
        end do
        call send_array(r4vector1,lsize,icpu,tag_base)
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      ksize = k2-k1+1
      lsize = isize*jsize*ksize*nnsg
      if ( size(r4vector2) < lsize ) then
        call getmem1d(r4vector2,1,lsize,'real4_3d_sub_distribute')
      end if
      window(1) = global_dot_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_dot_jstart+j1-1
      window(4) = window(3)+jsize-1
      call send_array(window,4,iocpu,tag_w)
      call recv_array(r4vector2,lsize,iocpu,tag_base)
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            do n = 1 , nnsg
              ml(n,j,i,k) = r4vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
  end subroutine real4_3d_sub_distribute
!
  subroutine integer_2d_sub_distribute(mg,ml,j1,j2,i1,i2)
    implicit none
    integer(ik4) , pointer , dimension(:,:,:) , intent(in) :: mg  ! model global
    integer(ik4) , pointer , dimension(:,:,:) , intent(out) :: ml ! model local
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2
    integer(ik4) :: ib , i , j , n , isize , jsize , lsize , icpu
    if ( myid == iocpu ) then
      ! Copy in memory my piece.
      do i = i1 , i2
        do j = j1 , j2
          do n = 1 , nnsg
            ml(n,j,i) = mg(n,global_dot_jstart+j-1,global_dot_istart+i-1)
          end do
        end do
      end do
      ! Send to other nodes the piece they request.
      do icpu = 1 , nproc-1
        call recv_array(window,4,icpu,tag_w)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        lsize = isize*jsize*nnsg
        if ( size(i4vector1) < lsize ) then
          call getmem1d(i4vector1,1,lsize,'integer_2d_sub_distribute')
        end if
        ib = 1
        do i = window(1) , window(2)
          do j = window(3) , window(4)
            do n = 1 , nnsg
              i4vector1(ib) = mg(n,j,i)
              ib = ib + 1
            end do
          end do
        end do
        call send_array(i4vector1,lsize,icpu,tag_base)
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      lsize = isize*jsize*nnsg
      if ( size(i4vector2) < lsize ) then
        call getmem1d(i4vector2,1,lsize,'integer_2d_sub_distribute')
      end if
      window(1) = global_dot_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_dot_jstart+j1-1
      window(4) = window(3)+jsize-1
      call send_array(window,4,iocpu,tag_w)
      call recv_array(i4vector2,lsize,iocpu,tag_base)
      ib = 1
      do i = i1 , i2
        do j = j1 , j2
          do n = 1 , nnsg
            ml(n,j,i) = i4vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine integer_2d_sub_distribute
!
  subroutine integer_3d_sub_distribute(mg,ml,j1,j2,i1,i2,k1,k2)
    implicit none
    integer(ik4) , pointer , dimension(:,:,:,:) , intent(in) :: mg  ! model glob
    integer(ik4) , pointer , dimension(:,:,:,:) , intent(out) :: ml ! model loc
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer(ik4) :: ib , i , j , k , n , isize , jsize , ksize , lsize , icpu
    if ( myid == iocpu ) then
      ! Copy in memory my piece.
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            do n = 1 , nnsg
              ml(n,j,i,k) = mg(n,global_dot_jstart+j-1,global_dot_istart+i-1,k)
            end do
          end do
        end do
      end do
      ! Send to other nodes the piece they request.
      do icpu = 1 , nproc-1
        call recv_array(window,4,icpu,tag_w)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        ksize = k2-k1+1
        lsize = isize*jsize*ksize*nnsg
        if ( size(i4vector1) < lsize ) then
          call getmem1d(i4vector1,1,lsize,'integer_3d_sub_distribute')
        end if
        ib = 1
        do k = k1 , k2
          do i = window(1) , window(2)
            do j = window(3) , window(4)
              do n = 1 , nnsg
                i4vector1(ib) = mg(n,j,i,k)
                ib = ib + 1
              end do
            end do
          end do
        end do
        call send_array(i4vector1,lsize,icpu,tag_base)
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      ksize = k2-k1+1
      lsize = isize*jsize*ksize*nnsg
      if ( size(i4vector2) < lsize ) then
        call getmem1d(i4vector2,1,lsize,'integer_3d_sub_distribute')
      end if
      window(1) = global_dot_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_dot_jstart+j1-1
      window(4) = window(3)+jsize-1
      call send_array(window,4,iocpu,tag_w)
      call recv_array(i4vector2,lsize,iocpu,tag_base)
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            do n = 1 , nnsg
              ml(n,j,i,k) = i4vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
  end subroutine integer_3d_sub_distribute
!
  subroutine real8_2d_3d_collect(ml,mg,j1,j2,i1,i2,k,lsub)
    implicit none
    real(rk8) , pointer , dimension(:,:) , intent(in) :: ml    ! model local
    real(rk8) , pointer , dimension(:,:,:) , intent(out) :: mg ! model global
    logical , intent(in) , optional :: lsub
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 
    integer(ik4) , intent(in) , optional :: k
    integer(ik4) :: ib , i , j , kk , isize , jsize , lsize , icpu
    integer(ik4) :: jsub1 , isub1
    if ( .not. present(lsub) .or. .not. lsub ) then
      jsub1 = global_dot_jstart
      isub1 = global_dot_istart
    else
      jsub1 = (global_dot_jstart-1)*nsg+1
      isub1 = (global_dot_istart-1)*nsg+1
    end if
    if ( myid == iocpu ) then
      kk = 1
      if ( present(k) ) kk = k
      ! Copy in memory my piece.
      do i = i1 , i2
        do j = j1 , j2
          mg(jsub1+j-1,isub1+i-1,kk) = ml(j,i)
        end do
      end do
      ! Receive from other nodes the piece they have
      do icpu = 1 , nproc-1
        call recv_array(window,4,icpu,tag_w)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        lsize = isize*jsize
        if ( size(r8vector1) < lsize ) then
          call getmem1d(r8vector1,1,lsize,'real8_2d_3d_collect')
        end if
        call recv_array(r8vector1,lsize,icpu,tag_base)
        ib = 1
        do i = window(1) , window(2)
          do j = window(3) , window(4)
            mg(j,i,kk) = r8vector1(ib)
            ib = ib + 1
          end do
        end do
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      lsize = isize*jsize
      if ( size(r8vector2) < lsize ) then
        call getmem1d(r8vector2,1,lsize,'real8_2d_3d_collect')
      end if
      window(1) = isub1+i1-1
      window(2) = window(1)+isize-1
      window(3) = jsub1+j1-1
      window(4) = window(3)+jsize-1
      call send_array(window,4,iocpu,tag_w)
      ib = 1
      do i = i1 , i2
        do j = j1 , j2
          r8vector2(ib) = ml(j,i)
          ib = ib + 1
         end do
      end do
      call send_array(r8vector2,lsize,iocpu,tag_base)
    end if
  end subroutine real8_2d_3d_collect
!
  subroutine real8_2d_collect(ml,mg,j1,j2,i1,i2,lsub)
    implicit none
    real(rk8) , pointer , dimension(:,:) , intent(in) :: ml  ! model local
    real(rk8) , pointer , dimension(:,:) , intent(out) :: mg ! model global
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2
    integer(ik4) :: ib , i , j , isize , jsize , lsize , icpu
    logical , intent(in) , optional :: lsub
    integer(ik4) :: jsub1 , isub1
    if ( .not. present(lsub) .or. .not. lsub ) then
      jsub1 = global_dot_jstart
      isub1 = global_dot_istart
    else
      jsub1 = (global_dot_jstart-1)*nsg+1
      isub1 = (global_dot_istart-1)*nsg+1
    end if
    if ( myid == iocpu ) then
      ! Copy in memory my piece.
      do i = i1 , i2
        do j = j1 , j2
          mg(jsub1+j-1,isub1+i-1) = ml(j,i)
        end do
      end do
      ! Receive from other nodes the piece they have
      do icpu = 1 , nproc-1
        call recv_array(window,4,icpu,tag_w)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        lsize = isize*jsize
        if ( size(r8vector1) < lsize ) then
          call getmem1d(r8vector1,1,lsize,'real8_2d_collect')
        end if
        call recv_array(r8vector1,lsize,icpu,tag_base)
        ib = 1
        do i = window(1) , window(2)
          do j = window(3) , window(4)
            mg(j,i) = r8vector1(ib)
            ib = ib + 1
          end do
        end do
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      lsize = isize*jsize
      if ( size(r8vector2) < lsize ) then
        call getmem1d(r8vector2,1,lsize,'real8_2d_collect')
      end if
      window(1) = isub1+i1-1
      window(2) = window(1)+isize-1
      window(3) = jsub1+j1-1
      window(4) = window(3)+jsize-1
      call send_array(window,4,iocpu,tag_w)
      ib = 1
      do i = i1 , i2
        do j = j1 , j2
          r8vector2(ib) = ml(j,i)
          ib = ib + 1
        end do
      end do
      call send_array(r8vector2,lsize,iocpu,tag_base)
    end if
  end subroutine real8_2d_collect
!
  subroutine real8_3d_collect(ml,mg,j1,j2,i1,i2,k1,k2,lsub)
    implicit none
    real(rk8) , pointer , dimension(:,:,:) , intent(in) :: ml  ! model local
    real(rk8) , pointer , dimension(:,:,:) , intent(out) :: mg ! model global
    logical , intent(in) , optional :: lsub
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer(ik4) :: ib , i , j , k , isize , jsize , ksize , lsize , icpu
    integer(ik4) :: isub1 , jsub1
    if ( .not. present(lsub) .or. .not. lsub ) then
      jsub1 = global_dot_jstart
      isub1 = global_dot_istart
    else
      jsub1 = (global_dot_jstart-1)*nsg+1
      isub1 = (global_dot_istart-1)*nsg+1
    end if
    if ( myid == iocpu ) then
      ! Copy in memory my piece.
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            mg(jsub1+j-1,isub1+i-1,k) = ml(j,i,k)
          end do
        end do
      end do
      ! Receive from other nodes the piece they have
      do icpu = 1 , nproc-1
        call recv_array(window,4,icpu,tag_w)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        ksize = k2-k1+1
        lsize = isize*jsize*ksize
        if ( size(r8vector1) < lsize ) then
          call getmem1d(r8vector1,1,lsize,'real8_3d_collect')
        end if
        call recv_array(r8vector1,lsize,icpu,tag_base)
        ib = 1
        do k = k1 , k2
          do i = window(1) , window(2)
            do j = window(3) , window(4)
              mg(j,i,k) = r8vector1(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      ksize = k2-k1+1
      lsize = isize*jsize*ksize
      if ( size(r8vector2) < lsize ) then
        call getmem1d(r8vector2,1,lsize,'real8_3d_collect')
      end if
      window(1) = isub1+i1-1
      window(2) = window(1)+isize-1
      window(3) = jsub1+j1-1
      window(4) = window(3)+jsize-1
      call send_array(window,4,iocpu,tag_w)
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            r8vector2(ib) = ml(j,i,k)
            ib = ib + 1
          end do
        end do
      end do
      call send_array(r8vector2,lsize,iocpu,tag_base)
    end if
  end subroutine real8_3d_collect
!
  subroutine real8_3d_2d_collect(ml,mg,j1,j2,i1,i2,k)
    implicit none
    real(rk8) , pointer , dimension(:,:,:) , intent(in) :: ml  ! model local
    real(rk8) , pointer , dimension(:,:) , intent(out) :: mg   ! model global
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k
    integer(ik4) :: ib , i , j , isize , jsize , lsize , icpu
    if ( myid == iocpu ) then
      ! Copy in memory my piece.
      do i = i1 , i2
        do j = j1 , j2
          mg(global_dot_jstart+j-1,global_dot_istart+i-1) = ml(j,i,k)
        end do
      end do
      ! Receive from other nodes the piece they have
      do icpu = 1 , nproc-1
        call recv_array(window,4,icpu,tag_w)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        lsize = isize*jsize
        if ( size(r8vector1) < lsize ) then
          call getmem1d(r8vector1,1,lsize,'real8_3d_collect')
        end if
        call recv_array(r8vector1,lsize,icpu,tag_base)
        ib = 1
        do i = window(1) , window(2)
          do j = window(3) , window(4)
            mg(j,i) = r8vector1(ib)
            ib = ib + 1
          end do
        end do
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      lsize = isize*jsize
      if ( size(r8vector2) < lsize ) then
        call getmem1d(r8vector2,1,lsize,'real8_3d_collect')
      end if
      window(1) = global_dot_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_dot_jstart+j1-1
      window(4) = window(3)+jsize-1
      call send_array(window,4,iocpu,tag_w)
      ib = 1
      do i = i1 , i2
        do j = j1 , j2
          r8vector2(ib) = ml(j,i,k)
          ib = ib + 1
        end do
      end do
      call send_array(r8vector2,lsize,iocpu,tag_base)
    end if
  end subroutine real8_3d_2d_collect
!
  subroutine real8_4d_collect(ml,mg,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:,:) , intent(in) :: ml  ! model local
    real(rk8) , pointer , dimension(:,:,:,:) , intent(out) :: mg ! model global
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2 , n1 , n2
    integer(ik4) :: ib , i , j , k , n , isize , jsize , &
                    ksize , nsize , lsize , icpu
    if ( myid == iocpu ) then
      ! Copy in memory my piece.
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = j1 , j2
              mg(global_dot_jstart+j-1,global_dot_istart+i-1,k,n) = ml(j,i,k,n)
            end do
          end do
        end do
      end do
      ! Receive from other nodes the piece they have
      do icpu = 1 , nproc-1
        call recv_array(window,4,icpu,tag_w)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        ksize = k2-k1+1
        nsize = n2-n1+1
        lsize = isize*jsize*ksize*nsize
        if ( size(r8vector1) < lsize ) then
          call getmem1d(r8vector1,1,lsize,'real8_4d_collect')
        end if
        call recv_array(r8vector1,lsize,icpu,tag_base)
        ib = 1
        do n = n1 , n2
          do k = k1 , k2
            do i = window(1) , window(2)
              do j = window(3) , window(4)
                mg(j,i,k,n) = r8vector1(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      ksize = k2-k1+1
      nsize = n2-n1+1
      lsize = isize*jsize*ksize*nsize
      if ( size(r8vector2) < lsize ) then
        call getmem1d(r8vector2,1,lsize,'real8_4d_collect')
      end if
      window(1) = global_dot_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_dot_jstart+j1-1
      window(4) = window(3)+jsize-1
      call send_array(window,4,iocpu,tag_w)
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = j1 , j2
              r8vector2(ib) = ml(j,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      call send_array(r8vector2,lsize,iocpu,tag_base)
    end if
  end subroutine real8_4d_collect
!
  subroutine real8_4d_2d_collect(ml,mg,j1,j2,i1,i2,k,n)
    implicit none
    real(rk8) , pointer , dimension(:,:,:,:) , intent(in) :: ml  ! model local
    real(rk8) , pointer , dimension(:,:) , intent(out) :: mg     ! model global
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k , n
    integer(ik4) :: ib , i , j , isize , jsize , lsize , icpu
    if ( myid == iocpu ) then
      ! Copy in memory my piece.
      do i = i1 , i2
        do j = j1 , j2
          mg(global_dot_jstart+j-1,global_dot_istart+i-1) = ml(j,i,k,n)
        end do
      end do
      ! Receive from other nodes the piece they have
      do icpu = 1 , nproc-1
        call recv_array(window,4,icpu,tag_w)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        lsize = isize*jsize
        if ( size(r8vector1) < lsize ) then
          call getmem1d(r8vector1,1,lsize,'real8_4d_collect')
        end if
        call recv_array(r8vector1,lsize,icpu,tag_base)
        ib = 1
        do i = window(1) , window(2)
          do j = window(3) , window(4)
            mg(j,i) = r8vector1(ib)
            ib = ib + 1
          end do
        end do
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      lsize = isize*jsize
      if ( size(r8vector2) < lsize ) then
        call getmem1d(r8vector2,1,lsize,'real8_4d_collect')
      end if
      window(1) = global_dot_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_dot_jstart+j1-1
      window(4) = window(3)+jsize-1
      call send_array(window,4,iocpu,tag_w)
      ib = 1
      do i = i1 , i2
        do j = j1 , j2
          r8vector2(ib) = ml(j,i,k,n)
          ib = ib + 1
        end do
      end do
      call send_array(r8vector2,lsize,iocpu,tag_base)
    end if
  end subroutine real8_4d_2d_collect
!
  subroutine real4_2d_collect(ml,mg,j1,j2,i1,i2)
    implicit none
    real(rk4) , pointer , dimension(:,:) , intent(in) :: ml  ! model local
    real(rk4) , pointer , dimension(:,:) , intent(out) :: mg ! model global
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2
    integer(ik4) :: ib , i , j , isize , jsize , lsize , icpu
    if ( myid == iocpu ) then
      ! Copy in memory my piece.
      do i = i1 , i2
        do j = j1 , j2
          mg(global_dot_jstart+j-1,global_dot_istart+i-1) = ml(j,i)
        end do
      end do
      ! Receive from other nodes the piece they have
      do icpu = 1 , nproc-1
        call recv_array(window,4,icpu,tag_w)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        lsize = isize*jsize
        if ( size(r4vector1) < lsize ) then
          call getmem1d(r4vector1,1,lsize,'real4_2d_collect')
        end if
        call recv_array(r4vector1,lsize,icpu,tag_base)
        ib = 1
        do i = window(1) , window(2)
          do j = window(3) , window(4)
            mg(j,i) = r4vector1(ib)
            ib = ib + 1
          end do
        end do
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      lsize = isize*jsize
      if ( size(r4vector2) < lsize ) then
        call getmem1d(r4vector2,1,lsize,'real4_2d_collect')
      end if
      window(1) = global_dot_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_dot_jstart+j1-1
      window(4) = window(3)+jsize-1
      call send_array(window,4,iocpu,tag_w)
      ib = 1
      do i = i1 , i2
        do j = j1 , j2
          r4vector2(ib) = ml(j,i)
          ib = ib + 1
        end do
      end do
      call send_array(r4vector2,lsize,iocpu,tag_base)
    end if
  end subroutine real4_2d_collect
!
  subroutine real4_3d_collect(ml,mg,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk4) , pointer , dimension(:,:,:) , intent(in) :: ml  ! model local
    real(rk4) , pointer , dimension(:,:,:) , intent(out) :: mg ! model global
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer(ik4) :: ib , i , j , k , isize , jsize , ksize , lsize , icpu
    if ( myid == iocpu ) then
      ! Copy in memory my piece.
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            mg(global_dot_jstart+j-1,global_dot_istart+i-1,k) = ml(j,i,k)
          end do
        end do
      end do
      ! Receive from other nodes the piece they have
      do icpu = 1 , nproc-1
        call recv_array(window,4,icpu,tag_w)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        ksize = k2-k1+1
        lsize = isize*jsize*ksize
        if ( size(r4vector1) < lsize ) then
          call getmem1d(r4vector1,1,lsize,'real4_3d_collect')
        end if
        call recv_array(r4vector1,lsize,icpu,tag_base)
        ib = 1
        do k = k1 , k2
          do i = window(1) , window(2)
            do j = window(3) , window(4)
              mg(j,i,k) = r4vector1(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      ksize = k2-k1+1
      lsize = isize*jsize*ksize
      if ( size(r4vector2) < lsize ) then
        call getmem1d(r4vector2,1,lsize,'real4_3d_collect')
      end if
      window(1) = global_dot_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_dot_jstart+j1-1
      window(4) = window(3)+jsize-1
      call send_array(window,4,iocpu,tag_w)
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            r4vector2(ib) = ml(j,i,k)
            ib = ib + 1
          end do
        end do
      end do
      call send_array(r4vector2,lsize,iocpu,tag_base)
    end if
  end subroutine real4_3d_collect
!
  subroutine real4_4d_collect(ml,mg,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk4) , pointer , dimension(:,:,:,:) , intent(in) :: ml  ! model local
    real(rk4) , pointer , dimension(:,:,:,:) , intent(out) :: mg ! model global
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2 , n1 , n2
    integer(ik4) :: ib , i , j , k , n , isize , jsize , &
                    ksize , nsize , lsize , icpu
    if ( myid == iocpu ) then
      ! Copy in memory my piece.
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = j1 , j2
              mg(global_dot_jstart+j-1,global_dot_istart+i-1,k,n) = ml(j,i,k,n)
            end do
          end do
        end do
      end do
      ! Receive from other nodes the piece they have
      do icpu = 1 , nproc-1
        call recv_array(window,4,icpu,tag_w)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        ksize = k2-k1+1
        nsize = n2-n1+1
        lsize = isize*jsize*ksize*nsize
        if ( size(r4vector1) < lsize ) then
          call getmem1d(r4vector1,1,lsize,'real4_4d_collect')
        end if
        call recv_array(r4vector1,lsize,icpu,tag_base)
        ib = 1
        do n = n1 , n2
          do k = k1 , k2
            do i = window(1) , window(2)
              do j = window(3) , window(4)
                mg(j,i,k,n) = r4vector1(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      ksize = k2-k1+1
      nsize = n2-n1+1
      lsize = isize*jsize*ksize*nsize
      if ( size(r4vector2) < lsize ) then
        call getmem1d(r4vector2,1,lsize,'real4_4d_collect')
      end if
      window(1) = global_dot_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_dot_jstart+j1-1
      window(4) = window(3)+jsize-1
      call send_array(window,4,iocpu,tag_w)
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = j1 , j2
              r4vector2(ib) = ml(j,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      call send_array(r4vector2,lsize,iocpu,tag_base)
    end if
  end subroutine real4_4d_collect
!
  subroutine integer_2d_collect(ml,mg,j1,j2,i1,i2)
    implicit none
    integer(ik4) , pointer , dimension(:,:) , intent(in) :: ml  ! model local
    integer(ik4) , pointer , dimension(:,:) , intent(out) :: mg ! model global
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2
    integer(ik4) :: ib , i , j , isize , jsize , lsize , icpu
    if ( myid == iocpu ) then
      ! Copy in memory my piece.
      do i = i1 , i2
        do j = j1 , j2
          mg(global_dot_jstart+j-1,global_dot_istart+i-1) = ml(j,i)
        end do
      end do
      ! Receive from other nodes the piece they have
      do icpu = 1 , nproc-1
        call recv_array(window,4,icpu,tag_w)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        lsize = isize*jsize
        if ( size(i4vector1) < lsize ) then
          call getmem1d(i4vector1,1,lsize,'integer_2d_collect')
        end if
        call recv_array(i4vector1,lsize,icpu,tag_base)
        ib = 1
        do i = window(1) , window(2)
          do j = window(3) , window(4)
            mg(j,i) = i4vector1(ib)
            ib = ib + 1
          end do
        end do
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      lsize = isize*jsize
      if ( size(i4vector2) < lsize ) then
        call getmem1d(i4vector2,1,lsize,'integer_2d_collect')
      end if
      window(1) = global_dot_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_dot_jstart+j1-1
      window(4) = window(3)+jsize-1
      call send_array(window,4,iocpu,tag_w)
      ib = 1
      do i = i1 , i2
        do j = j1 , j2
          i4vector2(ib) = ml(j,i)
          ib = ib + 1
        end do
      end do
      call send_array(i4vector2,lsize,iocpu,tag_base)
    end if
  end subroutine integer_2d_collect
!
  subroutine integer_3d_collect(ml,mg,j1,j2,i1,i2,k1,k2)
    implicit none
    integer(ik4) , pointer , dimension(:,:,:) , intent(in) :: ml  ! model local
    integer(ik4) , pointer , dimension(:,:,:) , intent(out) :: mg ! model global
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer(ik4) :: ib , i , j , k , isize , jsize , ksize , lsize , icpu
    if ( myid == iocpu ) then
      ! Copy in memory my piece.
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            mg(global_dot_jstart+j-1,global_dot_istart+i-1,k) = ml(j,i,k)
          end do
        end do
      end do
      ! Receive from other nodes the piece they have
      do icpu = 1 , nproc-1
        call recv_array(window,4,icpu,tag_w)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        ksize = k2-k1+1
        lsize = isize*jsize*ksize
        if ( size(i4vector1) < lsize ) then
          call getmem1d(i4vector1,1,lsize,'integer_3d_collect')
        end if
        call recv_array(i4vector1,lsize,icpu,tag_base)
        ib = 1
        do k = k1 , k2
          do i = window(1) , window(2)
            do j = window(3) , window(4)
              mg(j,i,k) = i4vector1(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      ksize = k2-k1+1
      lsize = isize*jsize*ksize
      if ( size(i4vector2) < lsize ) then
        call getmem1d(i4vector2,1,lsize,'integer_3d_collect')
      end if
      window(1) = global_dot_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_dot_jstart+j1-1
      window(4) = window(3)+jsize-1
      call send_array(window,4,iocpu,tag_w)
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            i4vector2(ib) = ml(j,i,k)
            ib = ib + 1
          end do
        end do
      end do
      call send_array(i4vector2,lsize,iocpu,tag_base)
    end if
  end subroutine integer_3d_collect
!
  subroutine integer_4d_collect(ml,mg,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    integer(ik4) , pointer , dimension(:,:,:,:) , intent(in) :: ml  ! model loc
    integer(ik4) , pointer , dimension(:,:,:,:) , intent(out) :: mg ! model glob
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2 , n1 , n2
    integer(ik4) :: ib , i , j , k , n , isize , jsize , &
                    ksize , nsize , lsize , icpu
    if ( myid == iocpu ) then
      ! Copy in memory my piece.
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = j1 , j2
              mg(global_dot_jstart+j-1,global_dot_istart+i-1,k,n) = ml(j,i,k,n)
            end do
          end do
        end do
      end do
      ! Receive from other nodes the piece they have
      do icpu = 1 , nproc-1
        call recv_array(window,4,icpu,tag_w)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        ksize = k2-k1+1
        nsize = n2-n1+1
        lsize = isize*jsize*ksize*nsize
        if ( size(i4vector1) < lsize ) then
          call getmem1d(i4vector1,1,lsize,'integer_4d_collect')
        end if
        call recv_array(i4vector1,lsize,icpu,tag_base)
        ib = 1
        do n = n1 , n2
          do k = k1 , k2
            do i = window(1) , window(2)
              do j = window(3) , window(4)
                mg(j,i,k,n) = i4vector1(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      ksize = k2-k1+1
      nsize = n2-n1+1
      lsize = isize*jsize*ksize*nsize
      if ( size(i4vector2) < lsize ) then
        call getmem1d(i4vector2,1,lsize,'integer_4d_collect')
      end if
      window(1) = global_dot_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_dot_jstart+j1-1
      window(4) = window(3)+jsize-1
      call send_array(window,4,iocpu,tag_w)
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = j1 , j2
              i4vector2(ib) = ml(j,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      call send_array(i4vector2,lsize,iocpu,tag_base)
    end if
  end subroutine integer_4d_collect
!
  subroutine real8_2d_sub_collect(ml,mg,j1,j2,i1,i2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:) , intent(in) :: ml  ! model local
    real(rk8) , pointer , dimension(:,:,:) , intent(out) :: mg ! model global
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2
    integer(ik4) :: ib , i , j , n , isize , jsize , lsize , icpu
    if ( myid == iocpu ) then
      ! Copy in memory my piece.
      do i = i1 , i2
        do j = j1 , j2
          do n = 1 , nnsg
            mg(n,global_dot_jstart+j-1,global_dot_istart+i-1) = ml(n,j,i)
          end do
        end do
      end do
      ! Receive from other nodes the piece they have
      do icpu = 1 , nproc-1
        call recv_array(window,4,icpu,tag_w)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        lsize = isize*jsize*nnsg
        if ( size(r8vector1) < lsize ) then
          call getmem1d(r8vector1,1,lsize,'real8_2d_sub_collect')
        end if
        call recv_array(r8vector1,lsize,icpu,tag_base)
        ib = 1
        do i = window(1) , window(2)
          do j = window(3) , window(4)
            do n = 1 , nnsg
              mg(n,j,i) = r8vector1(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      lsize = isize*jsize*nnsg
      if ( size(r8vector2) < lsize ) then
        call getmem1d(r8vector2,1,lsize,'real8_2d_sub_collect')
      end if
      window(1) = global_dot_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_dot_jstart+j1-1
      window(4) = window(3)+jsize-1
      call send_array(window,4,iocpu,tag_w)
      ib = 1
      do i = i1 , i2
        do j = j1 , j2
          do n = 1 , nnsg
            r8vector2(ib) = ml(n,j,i)
            ib = ib + 1
          end do
        end do
      end do
      call send_array(r8vector2,lsize,iocpu,tag_base)
    end if
  end subroutine real8_2d_sub_collect
!
  subroutine real8_3d_sub_collect(ml,mg,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:,:) , intent(in) :: ml  ! model local
    real(rk8) , pointer , dimension(:,:,:,:) , intent(out) :: mg ! model global
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer(ik4) :: ib , i , j , k , n , isize , jsize , ksize , lsize , icpu
    if ( myid == iocpu ) then
      ! Copy in memory my piece.
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            do n = 1 , nnsg
              mg(n,global_dot_jstart+j-1,global_dot_istart+i-1,k) = ml(n,j,i,k)
            end do
          end do
        end do
      end do
      ! Receive from other nodes the piece they have
      do icpu = 1 , nproc-1
        call recv_array(window,4,icpu,tag_w)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        ksize = k2-k1+1
        lsize = isize*jsize*ksize*nnsg
        if ( size(r8vector1) < lsize ) then
          call getmem1d(r8vector1,1,lsize,'real8_3d_sub_collect')
        end if
        call recv_array(r8vector1,lsize,icpu,tag_base)
        ib = 1
        do k = k1 , k2
          do i = window(1) , window(2)
            do j = window(3) , window(4)
              do n = 1 , nnsg
                mg(n,j,i,k) = r8vector1(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      ksize = k2-k1+1
      lsize = isize*jsize*ksize*nnsg
      if ( size(r8vector2) < lsize ) then
        call getmem1d(r8vector2,1,lsize,'real8_3d_sub_collect')
      end if
      window(1) = global_dot_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_dot_jstart+j1-1
      window(4) = window(3)+jsize-1
      call send_array(window,4,iocpu,tag_w)
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            do n = 1 , nnsg
              r8vector2(ib) = ml(n,j,i,k)
              ib = ib + 1
            end do
          end do
        end do
      end do
      call send_array(r8vector2,lsize,iocpu,tag_base)
    end if
  end subroutine real8_3d_sub_collect
!
  subroutine real4_2d_sub_collect(ml,mg,j1,j2,i1,i2)
    implicit none
    real(rk4) , pointer , dimension(:,:,:) , intent(in) :: ml  ! model local
    real(rk4) , pointer , dimension(:,:,:) , intent(out) :: mg ! model global
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2
    integer(ik4) :: ib , i , j , n , isize , jsize , lsize , icpu
    if ( myid == iocpu ) then
      ! Copy in memory my piece.
      do i = i1 , i2
        do j = j1 , j2
          do n = 1 , nnsg
            mg(n,global_dot_jstart+j-1,global_dot_istart+i-1) = ml(n,j,i)
          end do
        end do
      end do
      ! Receive from other nodes the piece they have
      do icpu = 1 , nproc-1
        call recv_array(window,4,icpu,tag_w)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        lsize = isize*jsize*nnsg
        if ( size(r4vector1) < lsize ) then
          call getmem1d(r4vector1,1,lsize,'real4_2d_sub_collect')
        end if
        call recv_array(r4vector1,lsize,icpu,tag_base)
        ib = 1
        do i = window(1) , window(2)
          do j = window(3) , window(4)
            do n = 1 , nnsg
              mg(n,j,i) = r4vector1(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      lsize = isize*jsize*nnsg
      if ( size(r4vector2) < lsize ) then
        call getmem1d(r4vector2,1,lsize,'real4_2d_sub_collect')
      end if
      window(1) = global_dot_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_dot_jstart+j1-1
      window(4) = window(3)+jsize-1
      call send_array(window,4,iocpu,tag_w)
      ib = 1
      do i = i1 , i2
        do j = j1 , j2
          do n = 1 , nnsg
            r4vector2(ib) = ml(n,j,i)
            ib = ib + 1
          end do
        end do
      end do
      call send_array(r4vector2,lsize,iocpu,tag_base)
    end if
  end subroutine real4_2d_sub_collect
!
  subroutine real4_3d_sub_collect(ml,mg,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk4) , pointer , dimension(:,:,:,:) , intent(in) :: ml  ! model local
    real(rk4) , pointer , dimension(:,:,:,:) , intent(out) :: mg ! model global
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer(ik4) :: ib , i , j , k , n , isize , jsize , ksize , lsize , icpu
    if ( myid == iocpu ) then
      ! Copy in memory my piece.
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            do n = 1 , nnsg
              mg(n,global_dot_jstart+j-1,global_dot_istart+i-1,k) = ml(n,j,i,k)
            end do
          end do
        end do
      end do
      ! Receive from other nodes the piece they have
      do icpu = 1 , nproc-1
        call recv_array(window,4,icpu,tag_w)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        ksize = k2-k1+1
        lsize = isize*jsize*ksize*nnsg
        if ( size(r4vector1) < lsize ) then
          call getmem1d(r4vector1,1,lsize,'real4_3d_sub_collect')
        end if
        call recv_array(r4vector1,lsize,icpu,tag_base)
        ib = 1
        do k = k1 , k2
          do i = window(1) , window(2)
            do j = window(3) , window(4)
              do n = 1 , nnsg
                mg(n,j,i,k) = r4vector1(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      ksize = k2-k1+1
      lsize = isize*jsize*ksize*nnsg
      if ( size(r4vector2) < lsize ) then
        call getmem1d(r4vector2,1,lsize,'real4_3d_sub_collect')
      end if
      window(1) = global_dot_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_dot_jstart+j1-1
      window(4) = window(3)+jsize-1
      call send_array(window,4,iocpu,tag_w)
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            do n = 1 , nnsg
              r4vector2(ib) = ml(n,j,i,k)
              ib = ib + 1
            end do
          end do
        end do
      end do
      call send_array(r4vector2,lsize,iocpu,tag_base)
    end if
  end subroutine real4_3d_sub_collect
!
  subroutine integer_2d_sub_collect(ml,mg,j1,j2,i1,i2)
    implicit none
    integer(ik4) , pointer , dimension(:,:,:) , intent(in) :: ml  ! model local
    integer(ik4) , pointer , dimension(:,:,:) , intent(out) :: mg ! model global
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2
    integer(ik4) :: ib , i , j , n , isize , jsize , lsize , icpu
    if ( myid == iocpu ) then
      ! Copy in memory my piece.
      do i = i1 , i2
        do j = j1 , j2
          do n = 1 , nnsg
            mg(n,global_dot_jstart+j-1,global_dot_istart+i-1) = ml(n,j,i)
          end do
        end do
      end do
      ! Receive from other nodes the piece they have
      do icpu = 1 , nproc-1
        call recv_array(window,4,icpu,tag_w)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        lsize = isize*jsize*nnsg
        if ( size(i4vector1) < lsize ) then
          call getmem1d(i4vector1,1,lsize,'integer_2d_sub_collect')
        end if
        call recv_array(i4vector1,lsize,icpu,tag_base)
        ib = 1
        do i = window(1) , window(2)
          do j = window(3) , window(4)
            do n = 1 , nnsg
              mg(n,j,i) = i4vector1(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      lsize = isize*jsize*nnsg
      if ( size(i4vector2) < lsize ) then
        call getmem1d(i4vector2,1,lsize,'integer_2d_sub_collect')
      end if
      window(1) = global_dot_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_dot_jstart+j1-1
      window(4) = window(3)+jsize-1
      call send_array(window,4,iocpu,tag_w)
      ib = 1
      do i = i1 , i2
        do j = j1 , j2
          do n = 1 , nnsg
            i4vector2(ib) = ml(n,j,i)
            ib = ib + 1
          end do
        end do
      end do
      call send_array(i4vector2,lsize,iocpu,tag_base)
    end if
  end subroutine integer_2d_sub_collect
!
  subroutine integer_3d_sub_collect(ml,mg,j1,j2,i1,i2,k1,k2)
    implicit none
    integer(ik4) , pointer , dimension(:,:,:,:) , intent(in) :: ml  ! model loc
    integer(ik4) , pointer , dimension(:,:,:,:) , intent(out) :: mg ! model glob
    integer(ik4) , intent(in) :: j1 , j2 , i1 , i2 , k1 , k2
    integer(ik4) :: ib , i , j , k , n , isize , jsize , ksize , lsize , icpu
    if ( myid == iocpu ) then
      ! Copy in memory my piece.
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            do n = 1 , nnsg
              mg(n,global_dot_jstart+j-1,global_dot_istart+i-1,k) = ml(n,j,i,k)
            end do
          end do
        end do
      end do
      ! Receive from other nodes the piece they have
      do icpu = 1 , nproc-1
        call recv_array(window,4,icpu,tag_w)
        isize = window(2)-window(1)+1
        jsize = window(4)-window(3)+1
        ksize = k2-k1+1
        lsize = isize*jsize*ksize*nnsg
        if ( size(i4vector1) < lsize ) then
          call getmem1d(i4vector1,1,lsize,'integer_3d_sub_collect')
        end if
        call recv_array(i4vector1,lsize,icpu,tag_base)
        ib = 1
        do k = k1 , k2
          do i = window(1) , window(2)
            do j = window(3) , window(4)
              do n = 1 , nnsg
                mg(n,j,i,k) = i4vector1(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do
      end do
    else
      isize = i2-i1+1
      jsize = j2-j1+1
      ksize = k2-k1+1
      lsize = isize*jsize*ksize*nnsg
      if ( size(i4vector2) < lsize ) then
        call getmem1d(i4vector2,1,lsize,'integer_3d_sub_collect')
      end if
      window(1) = global_dot_istart+i1-1
      window(2) = window(1)+isize-1
      window(3) = global_dot_jstart+j1-1
      window(4) = window(3)+jsize-1
      call send_array(window,4,iocpu,tag_w)
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = j1 , j2
            do n = 1 , nnsg
              i4vector2(ib) = ml(n,j,i,k)
              ib = ib + 1
            end do
          end do
        end do
      end do
      call send_array(i4vector2,lsize,iocpu,tag_base)
    end if
  end subroutine integer_3d_sub_collect
!
  subroutine real8_2d_exchange(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk8) , pointer , dimension(:,:) , intent(out) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2
    integer(ik4) :: isize , jsize , ssize , j , i , ib
    isize = i2-i1+1
    jsize = j2-j1+1
    if ( ma%bandflag ) then
      ssize = nex*isize
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange')
      end if
      if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange')
      end if
      ib = 1
      do i = i1 , i2
        do j = 1 , nex
          r8vector1(ib) = ml(j2-j+1,i)
          ib = ib + 1
        end do
      end do
      call cyclic_exchange_array(ssize,ma%right,ma%left,tag_lr)
      ib = 1
      do i = i1 , i2
        do j = 1 , nex
          ml(j1-j,i) = r8vector2(ib)
          ib = ib + 1
        end do
      end do
      ib = 1
      do i = i1 , i2
        do j = 1 , nex
          r8vector1(ib) = ml(j,i)
          ib = ib + 1
        end do
      end do
      call cyclic_exchange_array(ssize,ma%left,ma%right,tag_rl)
      ib = 1
      do i = i1 , i2
        do j = 1 , nex
          ml(j2+j,i) = r8vector2(ib)
          ib = ib + 1
        end do
      end do
    else
      if ( ma%right /= mpi_proc_null) then
        ssize = nex*isize
        if ( size(r8vector1) < ssize ) then
          call getmem1d(r8vector1,1,ssize,'real8_2d_exchange')
        end if
        if ( size(r8vector2) < ssize ) then
          call getmem1d(r8vector2,1,ssize,'real8_2d_exchange')
        end if
        ib = 1
        do i = i1 , i2
          do j = 1 , nex
            r8vector1(ib) = ml(j2-j+1,i)
            ib = ib + 1
          end do
        end do
        call exchange_array(ssize,ma%right,tag_rl,tag_lr)
        ib = 1
        do i = i1 , i2
          do j = 1 , nex
            ml(j2+j,i) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end if
      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*isize
        if ( size(r8vector2) < ssize ) then
          call getmem1d(r8vector2,1,ssize,'real8_2d_exchange')
        end if
        if ( size(r8vector1) < ssize ) then
          call getmem1d(r8vector1,1,ssize,'real8_2d_exchange')
        end if
        ib = 1
        do i = i1 , i2
          do j = 1 , nex
            r8vector1(ib) = ml(j,i)
            ib = ib + 1
          end do
        end do
        call exchange_array(ssize,ma%left,tag_lr,tag_rl)
        ib = 1
        do i = i1 , i2
          do j = 1 , nex
            ml(j1-j,i) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end if
    end if
    if ( ma%top /= mpi_proc_null) then
      ssize = nex*jsize
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange')
      end if
      if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange')
      end if
      ib = 1
      do i = 1 , nex
        do j = j1 , j2
          r8vector1(ib) = ml(j,i2-i+1)
          ib = ib + 1
        end do
      end do
      call exchange_array(ssize,ma%top,tag_tb,tag_bt)
      ib = 1
      do i = 1 , nex
        do j = j1 , j2
          ml(j,i2+i) = r8vector2(ib)
          ib = ib + 1
        end do
      end do
    end if
    if ( ma%bottom /= mpi_proc_null) then
      ssize = nex*jsize
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange')
      end if
      if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange')
      end if
      ib = 1
      do i = 1 , nex
        do j = j1 , j2
          r8vector1(ib) = ml(j,i)
          ib = ib + 1
        end do
      end do
      call exchange_array(ssize,ma%bottom,tag_bt,tag_tb)
      ib = 1
      do i = 1 , nex
        do j = j1 , j2
          ml(j,i1-i) = r8vector2(ib)
          ib = ib + 1
        end do
      end do
    end if
    if ( ma%topleft /= mpi_proc_null ) then
      ssize = nex*nex
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange')
      end if
      if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange')
      end if
      ib = 1
      do i = 1 , nex
        do j = 1 , nex
          r8vector1(ib) = ml(j,i2-i+1)
          ib = ib + 1
        end do
      end do
      call exchange_array(ssize,ma%topleft,tag_tlbr,tag_brtl)
      ib = 1
      do i = 1 , nex
        do j = 1 , nex
          ml(j1-j,i2+i) = r8vector2(ib)
          ib = ib + 1
        end do
      end do
    end if
    if ( ma%bottomright /= mpi_proc_null ) then
      ssize = nex*nex
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange')
      end if
      if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange')
      end if
      ib = 1
      do i = 1 , nex
        do j = 1 , nex
          r8vector1(ib) = ml(j2-j+1,i)
          ib = ib + 1
        end do
      end do
      call exchange_array(ssize,ma%bottomright,tag_brtl,tag_tlbr)
      ib = 1
      do i = 1 , nex
        do j = 1 , nex
          ml(j2+j,i1-i) = r8vector2(ib)
          ib = ib + 1
        end do
      end do
    end if
    if ( ma%topright /= mpi_proc_null ) then
      ssize = nex*nex
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange')
      end if
      if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange')
      end if
      ib = 1
      do i = 1 , nex
        do j = 1 , nex
          r8vector1(ib) = ml(j2-j+1,i2-i+1)
          ib = ib + 1
        end do
      end do
      call exchange_array(ssize,ma%topright,tag_trbl,tag_bltr)
      ib = 1
      do i = 1 , nex
        do j = 1 , nex
          ml(j2+j,i2+i) = r8vector2(ib)
          ib = ib + 1
        end do
      end do
    end if
    if ( ma%bottomleft /= mpi_proc_null ) then
      ssize = nex*nex
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange')
      end if
      if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange')
      end if
      ib = 1
      do i = 1 , nex
        do j = 1 , nex
          r8vector1(ib) = ml(j,i)
          ib = ib + 1
        end do
      end do
      call exchange_array(ssize,ma%bottomleft,tag_bltr,tag_trbl)
      ib = 1
      do i = 1 , nex
        do j = 1 , nex
          ml(j1-j,i1-i) = r8vector2(ib)
          ib = ib + 1
        end do
      end do
    end if
  end subroutine real8_2d_exchange
!
  subroutine real8_3d_exchange(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:) , intent(out) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2
    integer(ik4) :: isize , jsize , ksize , ssize , j , i , k , ib
    isize = i2-i1+1
    jsize = j2-j1+1
    ksize = k2-k1+1
    if ( ma%bandflag ) then
      ssize = nex*isize*ksize
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange')
      end if
      if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange')
      end if
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = 1 , nex
            r8vector1(ib) = ml(j2-j+1,i,k)
            ib = ib + 1
          end do
        end do
      end do
      call cyclic_exchange_array(ssize,ma%right,ma%left,tag_lr)
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = 1 , nex
            ml(j1-j,i,k) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = 1 , nex
            r8vector1(ib) = ml(j,i,k)
            ib = ib + 1
          end do
        end do
      end do
      call cyclic_exchange_array(ssize,ma%left,ma%right,tag_rl)
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = 1 , nex
            ml(j2+j,i,k) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    else
      if ( ma%right /= mpi_proc_null) then
        ssize = nex*isize*ksize
        if ( size(r8vector1) < ssize ) then
          call getmem1d(r8vector1,1,ssize,'real8_3d_exchange')
        end if
        if ( size(r8vector2) < ssize ) then
          call getmem1d(r8vector2,1,ssize,'real8_3d_exchange')
        end if
        ib = 1
        do k = k1 , k2
          do i = i1 , i2
            do j = 1 , nex
              r8vector1(ib) = ml(j2-j+1,i,k)
              ib = ib + 1
            end do
          end do
        end do
        call exchange_array(ssize,ma%right,tag_rl,tag_lr)
        ib = 1
        do k = k1 , k2
          do i = i1 , i2
            do j = 1 , nex
              ml(j2+j,i,k) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end if
      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*isize*ksize
        if ( size(r8vector2) < ssize ) then
          call getmem1d(r8vector2,1,ssize,'real8_3d_exchange')
        end if
        if ( size(r8vector1) < ssize ) then
          call getmem1d(r8vector1,1,ssize,'real8_3d_exchange')
        end if
        ib = 1
        do k = k1 , k2
          do i = i1 , i2
            do j = 1 , nex
              r8vector1(ib) = ml(j,i,k)
              ib = ib + 1
            end do
          end do
        end do
        call exchange_array(ssize,ma%left,tag_lr,tag_rl)
        ib = 1
        do k = k1 , k2
          do i = i1 , i2
            do j = 1 , nex
              ml(j1-j,i,k) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end if
    end if
    if ( ma%top /= mpi_proc_null) then
      ssize = nex*jsize*ksize
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange')
      end if
      if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange')
      end if
      ib = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = j1 , j2
            r8vector1(ib) = ml(j,i2-i+1,k)
            ib = ib + 1
          end do
        end do
      end do
      call exchange_array(ssize,ma%top,tag_tb,tag_bt)
      ib = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = j1 , j2
            ml(j,i2+i,k) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
    if ( ma%bottom /= mpi_proc_null) then
      ssize = nex*jsize*ksize
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange')
      end if
      if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange')
      end if
      ib = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = j1 , j2
            r8vector1(ib) = ml(j,i,k)
            ib = ib + 1
          end do
        end do
      end do
      call exchange_array(ssize,ma%bottom,tag_bt,tag_tb)
      ib = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = j1 , j2
            ml(j,i1-i,k) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
    if ( ma%topleft /= mpi_proc_null ) then
      ssize = nex*nex*ksize
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange')
      end if
      if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange')
      end if
      ib = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = 1 , nex
            r8vector1(ib) = ml(j,i2-i+1,k)
            ib = ib + 1
          end do
        end do
      end do
      call exchange_array(ssize,ma%topleft,tag_tlbr,tag_brtl)
      ib = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = 1 , nex
            ml(j1-j,i2+i,k) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
    if ( ma%bottomright /= mpi_proc_null ) then
      ssize = nex*nex*ksize
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange')
      end if
      if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange')
      end if
      ib = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = 1 , nex
            r8vector1(ib) = ml(j2-j+1,i,k)
            ib = ib + 1
          end do
        end do
      end do
      call exchange_array(ssize,ma%bottomright,tag_brtl,tag_tlbr)
      ib = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = 1 , nex
            ml(j2+j,i1-i,k) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
    if ( ma%topright /= mpi_proc_null ) then
      ssize = nex*nex*ksize
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange')
      end if
      if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange')
      end if
      ib = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = 1 , nex
            r8vector1(ib) = ml(j2-j+1,i2-i+1,k)
            ib = ib + 1
          end do
        end do
      end do
      call exchange_array(ssize,ma%topright,tag_trbl,tag_bltr)
      ib = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = 1 , nex
            ml(j2+j,i2+i,k) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
    if ( ma%bottomleft /= mpi_proc_null ) then
      ssize = nex*nex*ksize
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange')
      end if
      if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange')
      end if
      ib = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = 1 , nex
            r8vector1(ib) = ml(j,i,k)
            ib = ib + 1
          end do
        end do
      end do
      call exchange_array(ssize,ma%bottomleft,tag_bltr,tag_trbl)
      ib = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = 1 , nex
            ml(j1-j,i1-i,k) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine real8_3d_exchange
!
  subroutine real8_4d_exchange(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:,:) , intent(out) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2 , n1 , n2
    integer(ik4) :: isize , jsize , ksize , nsize , ssize
    integer(ik4) :: j , i , k , n , ib
    isize = i2-i1+1
    jsize = j2-j1+1
    ksize = k2-k1+1
    nsize = n2-n1+1
    if ( ma%bandflag ) then
      ssize = nex*isize*ksize*nsize
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange')
      end if
      if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange')
      end if
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = 1 , nex
              r8vector1(ib) = ml(j2-j+1,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      call cyclic_exchange_array(ssize,ma%right,ma%left,tag_lr)
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = 1 , nex
              ml(j1-j,i,k,n) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = 1 , nex
              r8vector1(ib) = ml(j,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      call cyclic_exchange_array(ssize,ma%left,ma%right,tag_rl)
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = 1 , nex
              ml(j2+j,i,k,n) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    else
      if ( ma%right /= mpi_proc_null) then
        ssize = nex*isize*ksize*nsize
        if ( size(r8vector1) < ssize ) then
          call getmem1d(r8vector1,1,ssize,'real8_4d_exchange')
        end if
        if ( size(r8vector2) < ssize ) then
          call getmem1d(r8vector2,1,ssize,'real8_4d_exchange')
        end if
        ib = 1
        do n = n1 , n2
          do k = k1 , k2
            do i = i1 , i2
              do j = 1 , nex
                r8vector1(ib) = ml(j2-j+1,i,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        call exchange_array(ssize,ma%right,tag_rl,tag_lr)
        ib = 1
        do n = n1 , n2
          do k = k1 , k2
            do i = i1 , i2
              do j = 1 , nex
                ml(j2+j,i,k,n) = r8vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do
      end if
      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*isize*ksize*nsize
        if ( size(r8vector2) < ssize ) then
          call getmem1d(r8vector2,1,ssize,'real8_4d_exchange')
        end if
        if ( size(r8vector1) < ssize ) then
          call getmem1d(r8vector1,1,ssize,'real8_4d_exchange')
        end if
        ib = 1
        do n = n1 , n2
          do k = k1 , k2
            do i = i1 , i2
              do j = 1 , nex
                r8vector1(ib) = ml(j,i,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        call exchange_array(ssize,ma%left,tag_lr,tag_rl)
        ib = 1
        do n = n1 , n2
          do k = k1 , k2
            do i = i1 , i2
              do j = 1 , nex
                ml(j1-j,i,k,n) = r8vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do
      end if
    end if
    if ( ma%top /= mpi_proc_null) then
      ssize = nex*jsize*ksize*nsize
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange')
      end if
      if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange')
      end if
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = j1 , j2
              r8vector1(ib) = ml(j,i2-i+1,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      call exchange_array(ssize,ma%top,tag_tb,tag_bt)
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = j1 , j2
              ml(j,i2+i,k,n) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
    if ( ma%bottom /= mpi_proc_null) then
      ssize = nex*jsize*ksize*nsize
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange')
      end if
      if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange')
      end if
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = j1 , j2
              r8vector1(ib) = ml(j,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      call exchange_array(ssize,ma%bottom,tag_bt,tag_tb)
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = j1 , j2
              ml(j,i1-i,k,n) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
    if ( ma%topleft /= mpi_proc_null ) then
      ssize = nex*nex*ksize*nsize
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange')
      end if
      if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange')
      end if
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = 1 , nex
              r8vector1(ib) = ml(j,i2-i+1,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      call exchange_array(ssize,ma%topleft,tag_tlbr,tag_brtl)
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = 1 , nex
              ml(j1-j,i2+i,k,n) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
    if ( ma%bottomright /= mpi_proc_null ) then
      ssize = nex*nex*ksize*nsize
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange')
      end if
      if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange')
      end if
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = 1 , nex
              r8vector1(ib) = ml(j2-j+1,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      call exchange_array(ssize,ma%bottomright,tag_brtl,tag_tlbr)
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = 1 , nex
              ml(j2+j,i1-i,k,n) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
    if ( ma%topright /= mpi_proc_null ) then
      ssize = nex*nex*ksize*nsize
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange')
      end if
      if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange')
      end if
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = 1 , nex
              r8vector1(ib) = ml(j2-j+1,i2-i+1,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      call exchange_array(ssize,ma%topright,tag_trbl,tag_bltr)
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = 1 , nex
              ml(j2+j,i2+i,k,n) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
    if ( ma%bottomleft /= mpi_proc_null ) then
      ssize = nex*nex*ksize*nsize
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange')
      end if
      if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange')
      end if
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = 1 , nex
              r8vector1(ib) = ml(j,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      call exchange_array(ssize,ma%bottomleft,tag_bltr,tag_trbl)
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = 1 , nex
              ml(j1-j,i1-i,k,n) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
  end subroutine real8_4d_exchange
!
  subroutine real8_2d_exchange_left_bottom(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk8) , pointer , dimension(:,:) , intent(out) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2
    integer(ik4) :: isize , jsize , ssize , j , i , ib
    isize = i2-i1+1
    jsize = j2-j1+1
    if ( ma%bandflag ) then
      ssize = nex*isize
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange_left_bottom')
      end if
      if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange_left_bottom')
      end if
      ib = 1
      do i = i1 , i2
        do j = 1 , nex
          r8vector1(ib) = ml(j2-j+1,i)
          ib = ib + 1
        end do
      end do
      call cyclic_exchange_array(ssize,ma%right,ma%left,tag_lr)
      ib = 1
      do i = i1 , i2
        do j = 1 , nex
          ml(j1-j,i) = r8vector2(ib)
          ib = ib + 1
        end do
      end do
    else
      if ( ma%right /= mpi_proc_null) then
        ssize = nex*isize
        if ( size(r8vector1) < ssize ) then
          call getmem1d(r8vector1,1,ssize,'real8_2d_exchange_left_bottom')
        end if
        ib = 1
        do i = i1 , i2
          do j = 1 , nex
            r8vector1(ib) = ml(j2-j+1,i)
            ib = ib + 1
          end do
        end do
        call send_array(r8vector1,ssize,ma%right,tag_lr)
      end if
      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*isize
        if ( size(r8vector2) < ssize ) then
          call getmem1d(r8vector2,1,ssize,'real8_2d_exchange_left_bottom')
        end if
        call recv_array(r8vector2,ssize,ma%left,tag_lr)
        ib = 1
        do i = i1 , i2
          do j = 1 , nex
            ml(j1-j,i) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end if
    end if
    if ( ma%top /= mpi_proc_null) then
      ssize = nex*jsize
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange_left_bottom')
      end if
      ib = 1
      do i = 1 , nex
        do j = j1 , j2
          r8vector1(ib) = ml(j,i2-i+1)
          ib = ib + 1
        end do
      end do
      call send_array(r8vector1,ssize,ma%top,tag_bt)
    end if
    if ( ma%bottom /= mpi_proc_null) then
      ssize = nex*jsize
      if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange_left_bottom')
      end if
      call recv_array(r8vector2,ssize,ma%bottom,tag_bt)
      ib = 1
      do i = 1 , nex
        do j = j1 , j2
          ml(j,i1-i) = r8vector2(ib)
          ib = ib + 1
        end do
      end do
    end if
    if ( ma%topright /= mpi_proc_null ) then
      ssize = nex*nex
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange_left_bottom')
      end if
      ib = 1
      do i = 1 , nex
        do j = 1 , nex
          r8vector1(ib) = ml(j2-j+1,i2-i+1)
          ib = ib + 1
        end do
      end do
      call send_array(r8vector1,ssize,ma%topright,tag_bltr)
    end if
    if ( ma%bottomleft /= mpi_proc_null ) then
      ssize = nex*nex
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange_left_bottom')
      end if
      call recv_array(r8vector2,ssize,ma%bottomleft,tag_bltr)
      ib = 1
      do i = 1 , nex
        do j = 1 , nex
          ml(j1-j,i1-i) = r8vector2(ib)
          ib = ib + 1
        end do
      end do
    end if
  end subroutine real8_2d_exchange_left_bottom
!
  subroutine real8_3d_exchange_left_bottom(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:) , intent(out) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2
    integer(ik4) :: isize , jsize , ksize , ssize , j , i , k , ib
    isize = i2-i1+1
    jsize = j2-j1+1
    ksize = k2-k1+1
    if ( ma%bandflag ) then
      ssize = nex*isize*ksize
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange_left_bottom')
      end if
      if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange_left_bottom')
      end if
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = 1 , nex
            r8vector1(ib) = ml(j2-j+1,i,k)
            ib = ib + 1
          end do
        end do
      end do
      call cyclic_exchange_array(ssize,ma%right,ma%left,tag_lr)
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = 1 , nex
            ml(j1-j,i,k) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    else
      if ( ma%right /= mpi_proc_null) then
        ssize = nex*isize*ksize
        if ( size(r8vector1) < ssize ) then
          call getmem1d(r8vector1,1,ssize,'real8_3d_exchange_left_bottom')
        end if
        ib = 1
        do k = k1 , k2
          do i = i1 , i2
            do j = 1 , nex
              r8vector1(ib) = ml(j2-j+1,i,k)
              ib = ib + 1
            end do
          end do
        end do
        call send_array(r8vector1,ssize,ma%right,tag_lr)
      end if
      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*isize*ksize
        if ( size(r8vector2) < ssize ) then
          call getmem1d(r8vector2,1,ssize,'real8_3d_exchange_left_bottom')
        end if
        call recv_array(r8vector2,ssize,ma%left,tag_lr)
        ib = 1
        do k = k1 , k2
          do i = i1 , i2
            do j = 1 , nex
              ml(j1-j,i,k) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end if
    end if
    if ( ma%top /= mpi_proc_null) then
      ssize = nex*jsize*ksize
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange_left_bottom')
      end if
      ib = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = j1 , j2
            r8vector1(ib) = ml(j,i2-i+1,k)
            ib = ib + 1
          end do
        end do
      end do
      call send_array(r8vector1,ssize,ma%top,tag_bt)
    end if
    if ( ma%bottom /= mpi_proc_null) then
      ssize = nex*jsize*ksize
      if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange_left_bottom')
      end if
      call recv_array(r8vector2,ssize,ma%bottom,tag_bt)
      ib = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = j1 , j2
            ml(j,i1-i,k) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
    if ( ma%topright /= mpi_proc_null ) then
      ssize = nex*nex*ksize
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange_left_bottom')
      end if
      ib = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = 1 , nex
            r8vector1(ib) = ml(j2-j+1,i2-i+1,k)
            ib = ib + 1
          end do
        end do
      end do
      call send_array(r8vector1,ssize,ma%topright,tag_bltr)
    end if
    if ( ma%bottomleft /= mpi_proc_null ) then
      ssize = nex*nex*ksize
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange_left_bottom')
      end if
      call recv_array(r8vector2,ssize,ma%bottomleft,tag_bltr)
      ib = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = 1 , nex
            ml(j1-j,i1-i,k) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine real8_3d_exchange_left_bottom
!
  subroutine real8_4d_exchange_left_bottom(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:,:) , intent(out) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2 , n1 , n2
    integer(ik4) :: isize , jsize , ksize , nsize , ssize , j , i , k , n , ib
    isize = i2-i1+1
    jsize = j2-j1+1
    ksize = k2-k1+1
    nsize = n2-n1+1
    if ( ma%bandflag ) then
      ssize = nex*isize*ksize*nsize
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange_left_bottom')
      end if
      if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange_left_bottom')
      end if
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = 1 , nex
              r8vector1(ib) = ml(j2-j+1,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      call cyclic_exchange_array(ssize,ma%right,ma%left,tag_lr)
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = 1 , nex
              ml(j1-j,i,k,n) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    else
      if ( ma%right /= mpi_proc_null) then
        ssize = nex*isize*ksize*nsize
        if ( size(r8vector1) < ssize ) then
          call getmem1d(r8vector1,1,ssize,'real8_4d_exchange_left_bottom')
        end if
        ib = 1
        do n = n1 , n2
          do k = k1 , k2
            do i = i1 , i2
              do j = 1 , nex
                r8vector1(ib) = ml(j2-j+1,i,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        call send_array(r8vector1,ssize,ma%right,tag_lr)
      end if
      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*isize*ksize*nsize
        if ( size(r8vector2) < ssize ) then
          call getmem1d(r8vector2,1,ssize,'real8_4d_exchange_left_bottom')
        end if
        call recv_array(r8vector2,ssize,ma%left,tag_lr)
        ib = 1
        do n = n1 , n2
          do k = k1 , k2
            do i = i1 , i2
              do j = 1 , nex
                ml(j1-j,i,k,n) = r8vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do
      end if
    end if
    if ( ma%top /= mpi_proc_null) then
      ssize = nex*jsize*ksize*nsize
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange_left_bottom')
      end if
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = j1 , j2
              r8vector1(ib) = ml(j,i2-i+1,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      call send_array(r8vector1,ssize,ma%top,tag_bt)
    end if
    if ( ma%bottom /= mpi_proc_null) then
      ssize = nex*jsize*ksize*nsize
      if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange_left_bottom')
      end if
      call recv_array(r8vector2,ssize,ma%bottom,tag_bt)
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = j1 , j2
              ml(j,i1-i,k,n) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
    if ( ma%topright /= mpi_proc_null ) then
      ssize = nex*nex*ksize*nsize
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange_left_bottom')
      end if
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = 1 , nex
              r8vector1(ib) = ml(j2-j+1,i2-i+1,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      call send_array(r8vector1,ssize,ma%topright,tag_bltr)
    end if
    if ( ma%bottomleft /= mpi_proc_null ) then
      ssize = nex*nex*ksize*nsize
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange_left_bottom')
      end if
      call recv_array(r8vector2,ssize,ma%bottomleft,tag_bltr)
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = 1 , nex
              ml(j1-j,i1-i,k,n) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
  end subroutine real8_4d_exchange_left_bottom
!
  subroutine real8_2d_exchange_right_top(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk8) , pointer , dimension(:,:) , intent(out) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2
    integer(ik4) :: isize , jsize , ssize , j , i , ib
    isize = i2-i1+1
    jsize = j2-j1+1
    if ( ma%bandflag ) then
      ssize = nex*isize
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange_right_top')
      end if
      if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange_right_top')
      end if
      ib = 1
      do i = i1 , i2
        do j = 1 , nex
          r8vector1(ib) = ml(j,i)
          ib = ib + 1
        end do
      end do
      call cyclic_exchange_array(ssize,ma%left,ma%right,tag_rl)
      ib = 1
      do i = i1 , i2
        do j = 1 , nex
          ml(j2+j,i) = r8vector2(ib)
          ib = ib + 1
        end do
      end do
    else
      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*isize
        if ( size(r8vector1) < ssize ) then
          call getmem1d(r8vector1,1,ssize,'real8_2d_exchange_right_top')
        end if
        ib = 1
        do i = i1 , i2
          do j = 1 , nex
            r8vector1(ib) = ml(j,i)
            ib = ib + 1
          end do
        end do
        call send_array(r8vector1,ssize,ma%left,tag_rl)
      end if
      if ( ma%right /= mpi_proc_null) then
        ssize = nex*isize
        if ( size(r8vector2) < ssize ) then
          call getmem1d(r8vector2,1,ssize,'real8_2d_exchange_right_top')
        end if
        call recv_array(r8vector2,ssize,ma%right,tag_rl)
        ib = 1
        do i = i1 , i2
          do j = 1 , nex
            ml(j2+j,i) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end if
    end if
    if ( ma%bottom /= mpi_proc_null) then
      ssize = nex*jsize
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange_right_top')
      end if
      ib = 1
      do i = 1 , nex
        do j = j1 , j2
          r8vector1(ib) = ml(j,i)
          ib = ib + 1
        end do
      end do
      call send_array(r8vector1,ssize,ma%bottom,tag_tb)
    end if
    if ( ma%top /= mpi_proc_null) then
      ssize = nex*jsize
      if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange_right_top')
      end if
      call recv_array(r8vector2,ssize,ma%top,tag_tb)
      ib = 1
      do i = 1 , nex
        do j = j1 , j2
          ml(j,i2+i) = r8vector2(ib)
          ib = ib + 1
        end do
      end do
    end if
    if ( ma%bottomleft /= mpi_proc_null ) then
      ssize = nex*nex
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_2d_exchange_right_top')
      end if
      ib = 1
      do i = 1 , nex
        do j = 1 , nex
          r8vector1(ib) = ml(j,i)
          ib = ib + 1
        end do
      end do
      call send_array(r8vector1,ssize,ma%bottomleft,tag_trbl)
    end if
    if ( ma%topright /= mpi_proc_null ) then
      ssize = nex*nex
      if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_2d_exchange_right_top')
      end if
      call recv_array(r8vector2,ssize,ma%topright,tag_trbl)
      ib = 1
      do i = 1 , nex
        do j = 1 , nex
          ml(j2+j,i2+i) = r8vector2(ib)
          ib = ib + 1
        end do
      end do
    end if
  end subroutine real8_2d_exchange_right_top
!
  subroutine real8_3d_exchange_right_top(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:) , intent(out) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2
    integer(ik4) :: isize , jsize , ksize , ssize , j , i , k , ib
    isize = i2-i1+1
    jsize = j2-j1+1
    ksize = k2-k1+1
    if ( ma%bandflag ) then
      ssize = nex*isize*ksize
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange_right_top')
      end if
      if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange_right_top')
      end if
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = 1 , nex
            r8vector1(ib) = ml(j,i,k)
            ib = ib + 1
          end do
        end do
      end do
      call cyclic_exchange_array(ssize,ma%left,ma%right,tag_rl)
      ib = 1
      do k = k1 , k2
        do i = i1 , i2
          do j = 1 , nex
            ml(j2+j,i,k) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    else
      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*isize*ksize
        if ( size(r8vector1) < ssize ) then
          call getmem1d(r8vector1,1,ssize,'real8_3d_exchange_right_top')
        end if
        ib = 1
        do k = k1 , k2
          do i = i1 , i2
            do j = 1 , nex
              r8vector1(ib) = ml(j,i,k)
              ib = ib + 1
            end do
          end do
        end do
        call send_array(r8vector1,ssize,ma%left,tag_rl)
      end if
      if ( ma%right /= mpi_proc_null) then
        ssize = nex*isize*ksize
        if ( size(r8vector2) < ssize ) then
          call getmem1d(r8vector2,1,ssize,'real8_3d_exchange_right_top')
        end if
        call recv_array(r8vector2,ssize,ma%right,tag_rl)
        ib = 1
        do k = k1 , k2
          do i = i1 , i2
            do j = 1 , nex
              ml(j2+j,i,k) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end if
    end if
    if ( ma%bottom /= mpi_proc_null) then
      ssize = nex*jsize*ksize
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange_right_top')
      end if
      ib = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = j1 , j2
            r8vector1(ib) = ml(j,i,k)
            ib = ib + 1
          end do
        end do
      end do
      call send_array(r8vector1,ssize,ma%bottom,tag_tb)
    end if
    if ( ma%top /= mpi_proc_null) then
      ssize = nex*jsize*ksize
      if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange_right_top')
      end if
      call recv_array(r8vector2,ssize,ma%top,tag_tb)
      ib = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = j1 , j2
            ml(j,i2+i,k) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
    if ( ma%bottomleft /= mpi_proc_null ) then
      ssize = nex*nex*ksize
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_3d_exchange_right_top')
      end if
      ib = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = 1 , nex
            r8vector1(ib) = ml(j,i,k)
            ib = ib + 1
          end do
        end do
      end do
      call send_array(r8vector1,ssize,ma%bottomleft,tag_trbl)
    end if
    if ( ma%topright /= mpi_proc_null ) then
      ssize = nex*nex*ksize
      if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_3d_exchange_right_top')
      end if
      call recv_array(r8vector2,ssize,ma%topright,tag_trbl)
      ib = 1
      do k = k1 , k2
        do i = 1 , nex
          do j = 1 , nex
            ml(j2+j,i2+i,k) = r8vector2(ib)
            ib = ib + 1
          end do
        end do
      end do
    end if
  end subroutine real8_3d_exchange_right_top
!
  subroutine real8_4d_exchange_right_top(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:,:) , intent(out) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2 , n1 , n2
    integer(ik4) :: isize , jsize , ksize , nsize , ssize , j , i , k , n , ib
    isize = i2-i1+1
    jsize = j2-j1+1
    ksize = k2-k1+1
    nsize = n2-n1+1
    if ( ma%bandflag ) then
      ssize = nex*isize*ksize*nsize
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange_right_top')
      end if
      if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange_right_top')
      end if
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = 1 , nex
              r8vector1(ib) = ml(j,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      call cyclic_exchange_array(ssize,ma%left,ma%right,tag_rl)
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = i1 , i2
            do j = 1 , nex
              ml(j2+j,i,k,n) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    else
      if ( ma%left /= mpi_proc_null ) then
        ssize = nex*isize*ksize*nsize
        if ( size(r8vector1) < ssize ) then
          call getmem1d(r8vector1,1,ssize,'real8_4d_exchange_right_top')
        end if
        ib = 1
        do n = n1 , n2
          do k = k1 , k2
            do i = i1 , i2
              do j = 1 , nex
                r8vector1(ib) = ml(j,i,k,n)
                ib = ib + 1
              end do
            end do
          end do
        end do
        call send_array(r8vector1,ssize,ma%left,tag_rl)
      end if
      if ( ma%right /= mpi_proc_null) then
        ssize = nex*isize*ksize*nsize
        if ( size(r8vector2) < ssize ) then
          call getmem1d(r8vector2,1,ssize,'real8_4d_exchange_right_top')
        end if
        call recv_array(r8vector2,ssize,ma%right,tag_rl)
        ib = 1
        do n = n1 , n2
          do k = k1 , k2
            do i = i1 , i2
              do j = 1 , nex
                ml(j2+j,i,k,n) = r8vector2(ib)
                ib = ib + 1
              end do
            end do
          end do
        end do
      end if
    end if
    if ( ma%bottom /= mpi_proc_null) then
      ssize = nex*jsize*ksize*nsize
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange_right_top')
      end if
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = j1 , j2
              r8vector1(ib) = ml(j,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      call send_array(r8vector1,ssize,ma%bottom,tag_tb)
    end if
    if ( ma%top /= mpi_proc_null) then
      ssize = nex*jsize*ksize*nsize
      if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange_right_top')
      end if
      call recv_array(r8vector2,ssize,ma%top,tag_tb)
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = j1 , j2
              ml(j,i2+i,k,n) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
    if ( ma%bottomleft /= mpi_proc_null ) then
      ssize = nex*nex*ksize*nsize
      if ( size(r8vector1) < ssize ) then
        call getmem1d(r8vector1,1,ssize,'real8_4d_exchange_right_top')
      end if
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = 1 , nex
              r8vector1(ib) = ml(j,i,k,n)
              ib = ib + 1
            end do
          end do
        end do
      end do
      call send_array(r8vector1,ssize,ma%bottomleft,tag_trbl)
    end if
    if ( ma%topright /= mpi_proc_null ) then
      ssize = nex*nex*ksize*nsize
      if ( size(r8vector2) < ssize ) then
        call getmem1d(r8vector2,1,ssize,'real8_4d_exchange_right_top')
      end if
      call recv_array(r8vector2,ssize,ma%topright,tag_trbl)
      ib = 1
      do n = n1 , n2
        do k = k1 , k2
          do i = 1 , nex
            do j = 1 , nex
              ml(j2+j,i2+i,k,n) = r8vector2(ib)
              ib = ib + 1
            end do
          end do
        end do
      end do
    end if
  end subroutine real8_4d_exchange_right_top
!
  subroutine real8_bdy_exchange_left_right(ml,k1,k2)
    implicit none
    real(rk8) , pointer , dimension(:,:) , intent(out) :: ml
    integer(ik4) , intent(in) :: k1 , k2
    integer(ik4) :: ksize , k , ib
    ksize = k2-k1+1
    if ( ma%bandflag ) then
      if ( size(r8vector1) < ksize ) then
        call getmem1d(r8vector1,1,ksize,'real8_bdy_exchange_left_right')
      end if
      if ( size(r8vector2) < ksize ) then
        call getmem1d(r8vector2,1,ksize,'real8_bdy_exchange_left_right')
      end if
      ib = 1
      do k = k1 , k2
        r8vector1(ib) = ml(jxp,k)
        ib = ib + 1
      end do
      call cyclic_exchange_array(ksize,ma%left,ma%right,tag_lr)
      ib = 1
      do k = k1 , k2
        ml(0,k) = r8vector2(ib)
        ib = ib + 1
      end do
      ib = 1
      do k = k1 , k2
        r8vector1(ib) = ml(1,k)
        ib = ib + 1
      end do
      call cyclic_exchange_array(ksize,ma%right,ma%left,tag_rl)
      ib = 1
      do k = k1 , k2
        ml(jxp+1,k) = r8vector2(ib)
        ib = ib + 1
      end do
    else
      if ( ma%right /= mpi_proc_null) then
        if ( size(r8vector1) < ksize ) then
          call getmem1d(r8vector1,1,ksize,'real8_bdy_exchange_left_right')
        end if
        if ( size(r8vector2) < ksize ) then
          call getmem1d(r8vector2,1,ksize,'real8_bdy_exchange_left_right')
        end if
        ib = 1
        do k = k1 , k2
          r8vector1(ib) = ml(jxp,k)
          ib = ib + 1
        end do
        call exchange_array(ksize,ma%right,tag_rl,tag_lr)
        ib = 1
        do k = k1 , k2
          ml(jxp+1,k) = r8vector2(ib)
          ib = ib + 1
        end do
      end if
      if ( ma%left /= mpi_proc_null ) then
        if ( size(r8vector2) < ksize ) then
          call getmem1d(r8vector2,1,ksize,'real8_bdy_exchange_left_right')
        end if
        if ( size(r8vector1) < ksize ) then
          call getmem1d(r8vector1,1,ksize,'real8_bdy_exchange_left_right')
        end if
        ib = 1
        do k = k1 , k2
          r8vector1(ib) = ml(1,k)
          ib = ib + 1
        end do
        call exchange_array(ksize,ma%left,tag_lr,tag_rl)
        ib = 1
        do k = k1 , k2
          ml(0,k) = r8vector2(ib)
          ib = ib + 1
        end do
      end if
    end if
  end subroutine real8_bdy_exchange_left_right
!
  subroutine real8_bdy_exchange_top_bottom(ml,k1,k2)
    implicit none
    real(rk8) , pointer , dimension(:,:) , intent(out) :: ml
    integer(ik4) , intent(in) :: k1 , k2
    integer(ik4) :: ksize , k , ib
    ksize = k2-k1+1
    if ( ma%top /= mpi_proc_null) then
      if ( size(r8vector1) < ksize ) then
        call getmem1d(r8vector1,1,ksize,'real8_bdy_exchange_top_bottom')
      end if
      if ( size(r8vector2) < ksize ) then
        call getmem1d(r8vector2,1,ksize,'real8_bdy_exchange_top_bottom')
      end if
      ib = 1
      do k = k1 , k2
        r8vector1(ib) = ml(iyp,k)
        ib = ib + 1
      end do
      call exchange_array(ksize,ma%top,tag_tb,tag_bt)
      ib = 1
      do k = k1 , k2
        ml(iyp+1,k) = r8vector2(ib)
        ib = ib + 1
      end do
    end if
    if ( ma%bottom /= mpi_proc_null ) then
      if ( size(r8vector2) < ksize ) then
        call getmem1d(r8vector2,1,ksize,'real8_bdy_exchange_top_bottom')
      end if
      if ( size(r8vector1) < ksize ) then
        call getmem1d(r8vector1,1,ksize,'real8_bdy_exchange_top_bottom')
      end if
      ib = 1
      do k = k1 , k2
        r8vector1(ib) = ml(1,k)
        ib = ib + 1
      end do
      call exchange_array(ksize,ma%bottom,tag_bt,tag_tb)
      ib = 1
      do k = k1 , k2
        ml(0,k) = r8vector2(ib)
        ib = ib + 1
      end do
    end if
  end subroutine real8_bdy_exchange_top_bottom
!
  subroutine real8_2d_grid_fill(a,b)
    implicit none
    real(rk8) , pointer , dimension(:,:) , intent(in) :: a
    real(rk8) , pointer , dimension(:,:) , intent(out) :: b
    call grid_collect(a,b,1,jxp,1,iyp)
    call mpi_bcast(b,iy*jx,mpi_real8,iocpu,cartesian_communicator,mpierr)
    if ( mpierr /= mpi_success ) then
      call fatal(__FILE__,__LINE__,'mpi_bcast error.')
    end if
  end subroutine real8_2d_grid_fill
!
! Takes u and v on the cross grid (the same grid as t, qv, qc, etc.)
! and interpolates the u and v to the dot grid.
! This routine sheilds the user of the function from the need to worry
! about the details of the domain decomposition.
!
! Written by Travis A. O'Brien 01/04/11.
!
  subroutine uvcross2dot(ux,vx,ud,vd)
    implicit none
    real(rk8) , pointer , dimension(:,:,:) , intent(out) :: ux , vx
    real(rk8) , pointer , dimension(:,:,:) , intent(out) :: ud , vd
    integer(ik4) :: i , j

    ! TODO:  It might make sense to encapsulate the following code
    ! in to a standard routine, since this boundary sending code is
    ! ubiquitous throughout the RegCM code and it is domain
    ! decomposition-dependent.

    ! Send the right-edge of the u/v tendencies to the left
    ! edge of the next process's u/v tendencies (so that
    ! invar%u(i,k,0) holds invar%u(i,k,jxp) of the parallel
    ! chunk next door)

    call exchange_lb(ux,1,jde1,jde2,ide1,ide2,1,kz)
    call exchange_lb(vx,1,jde1,jde2,ide1,ide2,1,kz)

    !
    !     x     x     x     x     x     x
    !
    !        o     o     o     o     o
    !         (i-1,j-1)     (i,j-1)
    !     x     x     x-----x     x     x
    !                 |(i,j)|
    !        o     o  |  o  |  o     o
    !                 |     |
    !     x     x     x-----x     x     x
    !           (i-1,j)     (i,j)
    !
    !        o     o     o     o     o
    !
    !     x     x     x     x     x     x
    !

    ! Perform the bilinear interpolation necessary
    ! to put the u and v variables on the dot grid.

    do i = idi1 , idi2
      do j = jdi1 , jdi2
        ud(j,i,:) =  ud(j,i,:) +             &
          d_rfour*(ux(j,i,:) + ux(j-1,i,:) +   &
                   ux(j,i-1,:) + ux(j-1,i-1,:))
        vd(j,i,:) =  vd(j,i,:) +             &
          d_rfour*(vx(j,i,:) + vx(j-1,i,:) +   &
                   vx(j,i-1,:) + vx(j-1,i-1,:))
      end do
    end do
  end subroutine uvcross2dot
!
  subroutine psc2psd(pc,pd)
    implicit none
    real(rk8) , pointer , dimension(:,:) , intent(in)  :: pc
    real(rk8) , pointer , dimension(:,:) , intent(out) :: pd
    integer(ik4) :: i , j
    !
    ! Internal points
    !
    do i = idi1 , idi2
      do j = jdi1 , jdi2
        pd(j,i) = (pc(j,i)+pc(j,i-1)+pc(j-1,i)+pc(j-1,i-1))*d_rfour
      end do
    end do
    !
    ! Boundaries
    !
    if ( ma%has_bdytop ) then
      do j = jdi1 , jdi2
        pd(j,ide2) = (pc(j,ice2)+pc(j-1,ice2))*d_half
      end do
    end if
    if ( ma%has_bdybottom ) then
      do j = jdi1 , jdi2
        pd(j,ide1)  = (pc(j,ice1)+pc(j-1,ice1))*d_half
      end do
    end if
    if ( ma%has_bdyleft ) then
      do i = idi1 , idi2
        pd(jde1,i) = (pc(jce1,i)+pc(jce1,i-1))*d_half
      end do
    end if
    if ( ma%has_bdyright ) then
      do i = idi1 , idi2
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
    character(len=*) , intent(in) :: varname
    logical , intent(in) :: ldot
    real(rk8) , pointer , dimension(:,:) , intent(in) :: val
    type (grid_nc_var2d) , intent(inout) :: xvar
    integer(ik4) :: istat
    integer(ik4) , dimension(3) :: idims

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
    type (grid_nc_var2d) , intent(inout) :: xvar
    integer(ik4) :: istat
    integer(ik4) , dimension(3) :: istart , icount
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
    type (grid_nc_var2d) , intent(inout) :: xvar
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
    character(len=*) , intent(in) :: varname
    logical , intent(in) :: ldot
    real(rk8) , pointer , dimension(:,:,:) , intent(in) :: val
    type (grid_nc_var3d) , intent(inout) :: xvar
    integer(ik4) :: istat
    integer(ik4) , dimension(4) :: idims

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
    type (grid_nc_var3d) , intent(inout) :: xvar
    integer(ik4) :: istat
    integer(ik4) , dimension(4) :: istart , icount
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
    type (grid_nc_var3d) , intent(inout) :: xvar
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
    character(len=*) , intent(in) :: varname
    logical , intent(in) :: ldot
    real(rk8) , pointer , dimension(:,:,:,:) , intent(in) :: val
    type (grid_nc_var4d) , intent(inout) :: xvar
    integer(ik4) :: istat
    integer(ik4) , dimension(5) :: idims

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
    type (grid_nc_var4d) , intent(inout) :: xvar
    integer(ik4) :: istat
    integer(ik4) , dimension(5) :: istart , icount
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
    type (grid_nc_var4d) , intent(inout) :: xvar
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
    real(rk8) , dimension(:) :: f_collect 
    real(rk8) :: f_sub
    call mpi_allgather(f_sub,    1,mpi_real8, &
                       f_collect,1,mpi_real8,mycomm,mpierr)
    if ( mpierr /= mpi_success ) THEN
      call fatal(__FILE__,__LINE__,'error in mpi_allgather!!')
    end if
  end subroutine gather_r

  subroutine gather_i(i_collect,i_sub)
    implicit none 
    integer(ik4) , dimension(:) :: i_collect 
    integer(ik4) :: i_sub
    call mpi_allgather(i_sub,    1,mpi_integer, &
                       i_collect,1,mpi_integer,mycomm,mpierr)
    if ( mpierr /= mpi_success ) THEN
      call fatal(__FILE__,__LINE__,'error in mpi_allgather!!')
    end if
  end subroutine gather_i

  subroutine reorder_add_subgrid_2d(var3,var2,mask)
    implicit none
    real(rk8) , pointer , dimension(:,:,:) , intent(in) :: var3
    real(rk8) , pointer , dimension(:,:) , intent(out) :: var2
    integer(ik4) , pointer , dimension(:,:,:) , intent(in) , optional :: mask
    integer(ik4) :: i , j , ii , jj , n1 , n2
    if ( present(mask) ) then
      do i = ici1 , ici2
        do j = jci1 , jci2
          do n2 = 1 , nsg
            ii = (i-1) * nsg + n2
            do n1 = 1 , nsg
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
      do i = ici1 , ici2
        do j = jci1 , jci2
          do n2 = 1 , nsg
            ii = (i-1) * nsg + n2
            do n1 = 1 , nsg
              jj = (j-1) * nsg + n1
              var2(jj,ii) = var2(jj,ii) + var3((n2-1)*nsg+n1,j,i)
            end do
          end do
        end do
      end do
    end if
  end subroutine reorder_add_subgrid_2d

  subroutine reorder_subgrid_2d(var3,var2,mask)
    implicit none
    real(rk8) , pointer , dimension(:,:,:) , intent(in) :: var3
    real(rk8) , pointer , dimension(:,:) , intent(out) :: var2
    integer(ik4) , pointer , dimension(:,:,:) , intent(in) , optional :: mask
    integer(ik4) :: i , j , ii , jj , n1 , n2
    if ( present(mask) ) then
      do i = ici1 , ici2
        do j = jci1 , jci2
          do n2 = 1 , nsg
            ii = (i-1) * nsg + n2
            do n1 = 1 , nsg
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
      do i = ici1 , ici2
        do j = jci1 , jci2
          do n2 = 1 , nsg
            ii = (i-1) * nsg + n2
            do n1 = 1 , nsg
              jj = (j-1) * nsg + n1
              var2(jj,ii) = var3((n2-1)*nsg+n1,j,i)
            end do
          end do
        end do
      end do
    end if
  end subroutine reorder_subgrid_2d

  subroutine reorder_add_subgrid_2d3d(var3,var2_3,l,mask)
    implicit none
    real(rk8) , pointer , dimension(:,:,:) , intent(in) :: var3
    real(rk8) , pointer , dimension(:,:,:) , intent(out) :: var2_3
    integer(ik4) , optional , intent(in) :: l
    integer(ik4) , pointer , dimension(:,:,:) , intent(in) , optional :: mask
    integer(ik4) :: i , j , ii , jj , n1 , n2 , ll
    ll = 1
    if ( present(l) ) ll = l
    if ( present(mask) ) then
      do i = ici1 , ici2
        do j = jci1 , jci2
          do n2 = 1 , nsg
            ii = (i-1) * nsg + n2
            do n1 = 1 , nsg
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
      do i = ici1 , ici2
        do j = jci1 , jci2
          do n2 = 1 , nsg
            ii = (i-1) * nsg + n2
            do n1 = 1 , nsg
              jj = (j-1) * nsg + n1
              var2_3(jj,ii,ll) = var2_3(jj,ii,ll) + var3((n2-1)*nsg+n1,j,i)
            end do
          end do
        end do
      end do
    end if
  end subroutine reorder_add_subgrid_2d3d

  subroutine reorder_subgrid_2d3d(var3,var2_3,l,mask)
    implicit none
    real(rk8) , pointer , dimension(:,:,:) , intent(in) :: var3
    real(rk8) , pointer , dimension(:,:,:) , intent(out) :: var2_3
    integer(ik4) , optional , intent(in) :: l
    integer(ik4) , pointer , dimension(:,:,:) , intent(in) , optional :: mask
    integer(ik4) :: i , j , ii , jj , n1 , n2 , ll
    ll = 1
    if ( present(l) ) ll = l
    if ( present(mask) ) then
      do i = ici1 , ici2
        do j = jci1 , jci2
          do n2 = 1 , nsg
            ii = (i-1) * nsg + n2
            do n1 = 1 , nsg
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
      do i = ici1 , ici2
        do j = jci1 , jci2
          do n2 = 1 , nsg
            ii = (i-1) * nsg + n2
            do n1 = 1 , nsg
              jj = (j-1) * nsg + n1
              var2_3(jj,ii,ll) = var3((n2-1)*nsg+n1,j,i)
            end do
          end do
        end do
      end do
    end if
  end subroutine reorder_subgrid_2d3d

  subroutine reorder_add_subgrid_3d(var4,var2,l,mask)
    implicit none
    real(rk8) , pointer , dimension(:,:,:,:) , intent(in) :: var4
    real(rk8) , pointer , dimension(:,:) , intent(out) :: var2
    integer(ik4) , intent(in) :: l
    integer(ik4) , pointer , dimension(:,:,:) , intent(in) , optional :: mask
    integer(ik4) :: i , j , ii , jj , n1 , n2
    if ( present(mask) ) then
      do i = ici1 , ici2
        do j = jci1 , jci2
          do n2 = 1 , nsg
            ii = (i-1) * nsg + n2
            do n1 = 1 , nsg
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
      do i = ici1 , ici2
        do j = jci1 , jci2
          do n2 = 1 , nsg
            ii = (i-1) * nsg + n2
            do n1 = 1 , nsg
              jj = (j-1) * nsg + n1
              var2(jj,ii) = var2(jj,ii) + var4((n2-1)*nsg+n1,j,i,l)
            end do
          end do
        end do
      end do
    end if
  end subroutine reorder_add_subgrid_3d

  subroutine reorder_subgrid_3d(var4,var2,l,mask)
    implicit none
    real(rk8) , pointer , dimension(:,:,:,:) , intent(in) :: var4
    real(rk8) , pointer , dimension(:,:) , intent(out) :: var2
    integer(ik4) , intent(in) :: l
    integer(ik4) , pointer , dimension(:,:,:) , intent(in) , optional :: mask
    integer(ik4) :: i , j , ii , jj , n1 , n2
    if ( present(mask) ) then
      do i = ici1 , ici2
        do j = jci1 , jci2
          do n2 = 1 , nsg
            ii = (i-1) * nsg + n2
            do n1 = 1 , nsg
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
      do i = ici1 , ici2
        do j = jci1 , jci2
          do n2 = 1 , nsg
            ii = (i-1) * nsg + n2
            do n1 = 1 , nsg
              jj = (j-1) * nsg + n1
              var2(jj,ii) = var4((n2-1)*nsg+n1,j,i,l)
            end do
          end do
        end do
      end do
    end if
  end subroutine reorder_subgrid_3d

  subroutine input_reorder(m1,m2)
    implicit none
    real(rk4) , dimension(:,:) , intent(in) :: m1
    real(rk8) , pointer , dimension(:,:,:) , intent(out) :: m2
    integer(ik4) :: i , j , ii , jj , n1 , n2
    do i = idot1 , idot2
      do j = jdot1 , jdot2
        do n2 = 1 , nsg
          ii = (i-1) * nsg + n2
          do n1 = 1 , nsg
            jj = (j-1) * nsg + n1
            m2((n2-1)*nsg+n1,j,i) = m1(jj,ii)
          end do
        end do
      end do
    end do
  end subroutine input_reorder

end module mod_mppparam

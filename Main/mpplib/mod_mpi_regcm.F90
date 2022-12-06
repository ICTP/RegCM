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

module mod_mpi_regcm
  use, intrinsic :: iso_fortran_env
  use mpi_f08
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_constants
  use mod_runparams , only : namelistfile , prgname
  use mod_mpmessage
  use mod_memutil
  use mod_date
  use mod_stdio
  use netcdf
  use mod_regcm_types

  implicit none

  private

  integer(ik4) , public , parameter :: iocpu = 0 ! The id of the cpu doing I/O
  integer(ik4) , public , parameter :: italk = 0 ! Who is doing the print ?

  public :: rcmpi
  public :: rcm_mpi_init
  public :: rcm_bcast
  public :: rcm_true
  public :: rcm_sum , rcm_max , rcm_min , rcm_mean
  public :: rcm_exch , rcm_exch_lrbt , rcm_exch_lr , rcm_exch_bt
  public :: rcm_exch_lb , rcm_exch_rt , rcm_exch_bdy_lr , rcm_exch_bdy_bt

  type regcm_mpi
    type(mpi_comm) :: globalcom
    type(mpi_comm) :: nodecom
    type(mpi_comm) :: cartcom
    integer :: myrank = 0
    integer :: cartrank = 0
    integer :: mysize = 0
    integer :: north = mpi_proc_null
    integer :: south = mpi_proc_null
    integer :: east  = mpi_proc_null
    integer :: west  = mpi_proc_null
    integer , dimension(2) :: pdims = [0,0]
    integer , dimension(2) :: pplace = [0,0]
    logical , dimension(2) :: period = [.false.,.false.]
  end type regcm_mpi

  type(regcm_mpi) :: rcmpi

  type regcm_exchange_handler
    type(mpi_request) :: req
  end type regcm_exchange_handler

  interface rcm_bcast
    module procedure bcast_logical,           &
                     bcast_int4,              &
                     bcast_int8,              &
                     bcast_real4,             &
                     bcast_real8,             &
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
  end interface rcm_bcast

  interface rcm_sum
    module procedure sumall_real4, &
                     sumall_real8 , &
#ifdef QUAD_PRECISION
                     sumall_real16, &
#endif
                     sumall_int4 ,  &
                     sumall_int8 ,  &
                     sumall_int4_array
  end interface rcm_sum

  interface rcm_max
    module procedure maxall_real8
    module procedure maxall_real4
    module procedure maxall_integer8
    module procedure maxall_integer4
  end interface rcm_max

  interface rcm_min
    module procedure minall_real8
    module procedure minall_real4
    module procedure minall_integer8
    module procedure minall_integer4
  end interface rcm_min

  interface rcm_mean
    module procedure meanall_real8
    module procedure meanall_real4
  end interface rcm_mean

  interface rcm_exch
    module procedure real8_2d_exchange, &
                     real8_3d_exchange, &
                     real8_4d_exchange, &
                     real4_2d_exchange, &
                     real4_3d_exchange, &
                     real4_4d_exchange
  end interface rcm_exch

  interface rcm_exch_lrbt
    module procedure real8_2d_exchange_lrbt, &
                     real8_3d_exchange_lrbt, &
                     real8_4d_exchange_lrbt, &
                     real4_2d_exchange_lrbt, &
                     real4_3d_exchange_lrbt, &
                     real4_4d_exchange_lrbt
  end interface rcm_exch_lrbt

  interface rcm_exch_lr
    module procedure real8_2d_exchange_lr, &
                     real8_3d_exchange_lr, &
                     real8_4d_exchange_lr, &
                     real4_2d_exchange_lr, &
                     real4_3d_exchange_lr, &
                     real4_4d_exchange_lr
  end interface rcm_exch_lr

  interface rcm_exch_bt
    module procedure real8_2d_exchange_bt, &
                     real8_3d_exchange_bt, &
                     real8_4d_exchange_bt, &
                     real4_2d_exchange_bt, &
                     real4_3d_exchange_bt, &
                     real4_4d_exchange_bt
  end interface rcm_exch_bt

  interface rcm_exch_lb
    module procedure real8_2d_exchange_lb, &
                     real8_3d_exchange_lb, &
                     real8_4d_exchange_lb, &
                     real4_2d_exchange_lb, &
                     real4_3d_exchange_lb, &
                     real4_4d_exchange_lb
  end interface rcm_exch_lb

  interface rcm_exch_rt
    module procedure real8_2d_exchange_rt, &
                     real8_3d_exchange_rt, &
                     real8_4d_exchange_rt, &
                     real4_2d_exchange_rt, &
                     real4_3d_exchange_rt, &
                     real4_4d_exchange_rt
  end interface rcm_exch_rt

  interface rcm_exch_bdy_lr
    module procedure real8_exchange_bdy_lr, &
                     real4_exchange_bdy_lr
  end interface rcm_exch_bdy_lr

  interface rcm_exch_bdy_bt
    module procedure real8_exchange_bdy_bt, &
                     real4_exchange_bdy_bt
  end interface rcm_exch_bdy_bt

  contains

  subroutine rcm_mpi_init(startcomm)
    implicit none
    type(mpi_comm), intent(in) :: startcomm

    if ( startcomm == MPI_COMM_WORLD ) then
      call mpi_comm_dup(startcomm,rcmpi%globalcom)
    else
      rcmpi%globalcom = startcomm
    end if
    if ( i_band == 1 ) rcmpi%period(1) = .true.
    if ( i_crm == 1 ) rcmpi%period(:) = .true.
    call mpi_comm_size(rcmpi%globalcom, rcmpi%mysize)
    call mpi_comm_rank(rcmpi%globalcom, rcmpi%myrank)
    call mpi_comm_split_type(rcmpi%globalcom, mpi_comm_type_shared, &
                             0, mpi_info_null, rcmpi%nodecom)
    call mpi_dims_create(rcmpi%mysize,2,rcmpi%pdims)
    call mpi_cart_create(rcmpi%globalcom, 2, rcmpi%pdims, rcmpi%period, &
                         .false., rcmpi%cartcom)
    call mpi_comm_rank(rcmpi%cartcom,rcmpi%cartrank)
    call mpi_cart_coords(rcmpi%cartcom,rcmpi%cartrank,2,rcmpi%pplace)
    call mpi_cart_shift(rcmpi%cartcom, 0, 1, rcmpi%west, rcmpi%east)
    call mpi_cart_shift(rcmpi%cartcom, 1, 1, rcmpi%south, rcmpi%north)

  end subroutine rcm_mpi_init

  subroutine bcast_logical(x)
    implicit none
    logical , intent(inout) :: x
    call mpi_bcast(x,1,mpi_logical,iocpu,rcmpi%globalcom)
  end subroutine bcast_logical
  subroutine bcast_int4(x)
    implicit none
    integer(ik4) , intent(inout) :: x
    call mpi_bcast(x,1,mpi_integer4,iocpu,rcmpi%globalcom)
  end subroutine bcast_int4
  subroutine bcast_int8(x)
    implicit none
    integer(rk8) , intent(inout) :: x
    call mpi_bcast(x,1,mpi_integer8,iocpu,rcmpi%globalcom)
  end subroutine bcast_int8
  subroutine bcast_real4(x)
    implicit none
    real(rk4) , intent(inout) :: x
    call mpi_bcast(x,1,mpi_real4,iocpu,rcmpi%globalcom)
  end subroutine bcast_real4
  subroutine bcast_real8(x)
    implicit none
    real(rk8) , intent(inout) :: x
    call mpi_bcast(x,1,mpi_real8,iocpu,rcmpi%globalcom)
  end subroutine bcast_real8
  subroutine bcast_arr_logical(x)
    implicit none
    logical , dimension(:) , intent(inout) :: x
    call mpi_bcast(x,size(x),mpi_logical,iocpu,rcmpi%globalcom)
  end subroutine bcast_arr_logical
  subroutine bcast_arr_character(cval,is)
    implicit none
    character(len=*) , intent(inout) :: cval
    integer(ik4) , intent(in) :: is
    call mpi_bcast(cval,is,mpi_character,iocpu,rcmpi%globalcom)
  end subroutine bcast_arr_character
  subroutine bcast_arr_text_list(cval,is)
    implicit none
    character(len=*) , intent(inout) , dimension(:) :: cval
    integer(ik4) , intent(in) :: is
    call mpi_bcast(cval,is*size(cval),mpi_character,iocpu,rcmpi%globalcom)
  end subroutine bcast_arr_text_list
  subroutine bcast_arr_int4(x)
    implicit none
    integer(ik4) , dimension(:) , intent(inout) :: x
    call mpi_bcast(x,size(x),mpi_integer4,iocpu,rcmpi%globalcom)
  end subroutine bcast_arr_int4
  subroutine bcast_arr_int8(x)
    implicit none
    integer(rk8) , dimension(:) , intent(inout) :: x
    call mpi_bcast(x,size(x),mpi_integer8,iocpu,rcmpi%globalcom)
  end subroutine bcast_arr_int8
  subroutine bcast_arr_real4(x)
    implicit none
    real(rk4) , dimension(:) , intent(inout) :: x
    call mpi_bcast(x,size(x),mpi_real4,iocpu,rcmpi%globalcom)
  end subroutine bcast_arr_real4
  subroutine bcast_arr_real8(x)
    implicit none
    real(rk8) , dimension(:) , intent(inout) :: x
    call mpi_bcast(x,size(x),mpi_real8,iocpu,rcmpi%globalcom)
  end subroutine bcast_arr_real8
  subroutine bcast_matr_real8(x)
    implicit none
    real(rk8) , dimension(:,:) , intent(inout) :: x
    call mpi_bcast(x,size(x,1)*size(x,2), &
                   mpi_real8,iocpu,rcmpi%globalcom)
  end subroutine bcast_matr_real8
  subroutine bcast_matr_real4(x)
    implicit none
    real(rk4) , dimension(:,:) , intent(inout) :: x
    call mpi_bcast(x,size(x,1)*size(x,2), &
                   mpi_real4,iocpu,rcmpi%globalcom)
  end subroutine bcast_matr_real4
  subroutine bcast_rcm_time_and_date(x)
    implicit none
    type (rcm_time_and_date) , intent(inout) :: x
    call rcm_bcast(x%calendar)
    call rcm_bcast(x%days_from_reference)
    call rcm_bcast(x%second_of_day)
  end subroutine bcast_rcm_time_and_date
  subroutine bcast_arr_rcm_time_and_date(x)
    implicit none
    type (rcm_time_and_date) , dimension(:) , intent(inout) :: x
    integer(ik4) :: n
    do n = 1 , size(x)
      call rcm_bcast(x(n))
    end do
  end subroutine bcast_arr_rcm_time_and_date

  subroutine rcm_true(l,g)
    implicit none
    logical , intent(in) :: l
    logical , intent(out) :: g
    call mpi_allreduce(l,g,1,mpi_logical,mpi_lor,rcmpi%globalcom)
  end subroutine rcm_true

#ifdef QUAD_PRECISION
  subroutine sumall_real16(l,g)
    implicit none
    real(rk16) , intent(in) :: l
    real(rk16) , intent(out) :: g
    call mpi_allreduce(l,g,1,mpi_real16,mpi_sum,rcmpi%globalcom)
  end subroutine sumall_real16
#endif
  subroutine sumall_int8(l,g)
    implicit none
    integer(ik8) , intent(in) :: l
    integer(ik8) , intent(out) :: g
    call mpi_allreduce(l,g,1,mpi_integer8,mpi_sum,rcmpi%globalcom)
  end subroutine sumall_int8
  subroutine sumall_int4(l,g)
    implicit none
    integer(ik4) , intent(in) :: l
    integer(ik4) , intent(out) :: g
    call mpi_allreduce(l,g,1,mpi_integer4,mpi_sum,rcmpi%globalcom)
  end subroutine sumall_int4
  subroutine sumall_real8(l,g)
    implicit none
    real(rk8) , intent(in) :: l
    real(rk8) , intent(out) :: g
    call mpi_allreduce(l,g,1,mpi_real8,mpi_sum,rcmpi%globalcom)
  end subroutine sumall_real8
  subroutine sumall_real4(l,g)
    implicit none
    real(rk4) , intent(in) :: l
    real(rk4) , intent(out) :: g
    call mpi_allreduce(l,g,1,mpi_real4,mpi_sum,rcmpi%globalcom)
  end subroutine sumall_real4
  subroutine sumall_int4_array(l,g)
    implicit none
    integer(ik4) , dimension(:) , intent(in) :: l
    integer(ik4) , dimension(:) , intent(out) :: g
    call mpi_allreduce(l,g,size(g),mpi_integer4, &
                       mpi_sum,rcmpi%globalcom)
  end subroutine sumall_int4_array

  subroutine maxall_real8(l,g)
    implicit none
    real(rk8) , intent(in) :: l
    real(rk8) , intent(out) :: g
    call mpi_allreduce(l,g,1,mpi_real8,mpi_max,rcmpi%globalcom)
  end subroutine maxall_real8
  subroutine maxall_real4(l,g)
    implicit none
    real(rk4) , intent(in) :: l
    real(rk4) , intent(out) :: g
    call mpi_allreduce(l,g,1,mpi_real4,mpi_max,rcmpi%globalcom)
  end subroutine maxall_real4
  subroutine maxall_integer4(l,g)
    implicit none
    integer(ik4) , intent(in) :: l
    integer(ik4) , intent(out) :: g
    call mpi_allreduce(l,g,1,mpi_integer4,mpi_max,rcmpi%globalcom)
  end subroutine maxall_integer4
  subroutine maxall_integer8(l,g)
    implicit none
    integer(ik8) , intent(in) :: l
    integer(ik8) , intent(out) :: g
    call mpi_allreduce(l,g,1,mpi_integer8,mpi_max,rcmpi%globalcom)
  end subroutine maxall_integer8

  subroutine minall_real8(l,g)
    implicit none
    real(rk8) , intent(in) :: l
    real(rk8) , intent(out) :: g
    call mpi_allreduce(l,g,1,mpi_real8,mpi_min,rcmpi%globalcom)
  end subroutine minall_real8
  subroutine minall_real4(l,g)
    implicit none
    real(rk4) , intent(in) :: l
    real(rk4) , intent(out) :: g
    call mpi_allreduce(l,g,1,mpi_real4,mpi_min,rcmpi%globalcom)
  end subroutine minall_real4
  subroutine minall_integer4(l,g)
    implicit none
    integer(ik4) , intent(in) :: l
    integer(ik4) , intent(out) :: g
    call mpi_allreduce(l,g,1,mpi_integer4,mpi_min,rcmpi%globalcom)
  end subroutine minall_integer4
  subroutine minall_integer8(l,g)
    implicit none
    integer(ik8) , intent(in) :: l
    integer(ik8) , intent(out) :: g
    call mpi_allreduce(l,g,1,mpi_integer8,mpi_min,rcmpi%globalcom)
  end subroutine minall_integer8

  subroutine meanall_real8(l,g)
    implicit none
    real(rk8) , intent(in) :: l
    real(rk8) , intent(out) :: g
    call mpi_allreduce(l,g,1,mpi_real8,mpi_sum,rcmpi%globalcom)
    g = g/real(nproc,rk8)
  end subroutine meanall_real8
  subroutine meanall_real4(l,g)
    implicit none
    real(rk4) , intent(in) :: l
    real(rk4) , intent(out) :: g
    call mpi_allreduce(l,g,1,mpi_real4,mpi_sum,rcmpi%globalcom)
    g = g/real(nproc,rk4)
  end subroutine meanall_real4

  subroutine real8_2d_exchange(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk8) , pointer , dimension(:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2
    integer :: nx , ny
    integer :: lb , rb , tb , bb
    integer :: ndx , ndy , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    lb = nex
    rb = nex
    tb = 0
    bb = 0
    if ( rcmpi%west == mpi_proc_null ) lb = 1
    if ( rcmpi%east == mpi_proc_null ) rb = 1
    if ( rcmpi%north == mpi_proc_null ) tb = 1
    if ( rcmpi%south == mpi_proc_null ) bb = 1
    tx = ny+tb+bb
    ty = nx+lb+rb
    sizex = nex*tx
    sizey = nex*ty
    ndx = 2*sizex
    ndy = 2*sizey

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk8), dimension(ndx), asynchronous :: sdatax
    real(rk8), dimension(ndx), asynchronous :: rdatax
    real(rk8), dimension(ndy), asynchronous :: sdatay
    real(rk8), dimension(ndy), asynchronous :: rdatay
    type(mpi_request) :: mrequestx, mrequesty
    type(mpi_status) :: mstatus
    integer :: ib1 , ib2 , iex

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      if ( rcmpi%west /= mpi_proc_null ) &
          sdatax(ib1:ib2) = ml(j1+(iex-1),i1-bb:i2+tb)
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      if ( rcmpi%east /= mpi_proc_null ) &
          sdatax(ib1:ib2) = ml(j2-(iex-1),i1-bb:i2+tb)
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_ineighbor_alltoallv(sdatax, counts, displs, mpi_real8, &
        rdatax, counts, displs, mpi_real8, rcmpi%cartcom, mrequestx)

    call mpi_wait(mrequestx,mstatus)

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
       if ( rcmpi%west /= mpi_proc_null ) &
           ml(j1-iex,i1-bb:i2+tb) = rdatax(ib1:ib2)
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      if ( rcmpi%east /= mpi_proc_null ) &
          ml(j2+iex,i1-bb:i2+tb) = rdatax(ib1:ib2)
    end do

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      if ( rcmpi%south /= mpi_proc_null ) &
          sdatay(ib1:ib2) = ml(:,i1+(iex-1))
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      if ( rcmpi%north /= mpi_proc_null ) &
          sdatay(ib1:ib2) = ml(:,i2-(iex-1))
    end do
    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_ineighbor_alltoallv(sdatay, counts, displs, mpi_real8, &
        rdatay, counts, displs, mpi_real8, rcmpi%cartcom, mrequesty)

    call mpi_wait(mrequesty,mstatus)

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      if ( rcmpi%south /= mpi_proc_null ) &
          ml(:,i1-iex) = rdatay(ib1:ib2)
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      if ( rcmpi%north /= mpi_proc_null ) &
         ml(:,i2+iex) = rdatay(ib1:ib2)
    end do

    end block transmit

  end subroutine real8_2d_exchange

  subroutine real8_3d_exchange(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2
    integer :: nx , ny , nk
    integer :: lb , rb , tb , bb
    integer :: ndx , ndy , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    lb = nex
    rb = nex
    tb = 0
    bb = 0
    if ( rcmpi%west == mpi_proc_null ) lb = 1
    if ( rcmpi%east == mpi_proc_null ) rb = 1
    if ( rcmpi%north == mpi_proc_null ) tb = 1
    if ( rcmpi%south == mpi_proc_null ) bb = 1
    tx = ny+bb+tb
    ty = nx+lb+rb
    sizex = nex*tx*nk
    sizey = nex*ty*nk
    ndx = 2*sizex
    ndy = 2*sizey

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk8), dimension(ndx), asynchronous :: sdatax
    real(rk8), dimension(ndx), asynchronous :: rdatax
    real(rk8), dimension(ndy), asynchronous :: sdatay
    real(rk8), dimension(ndy), asynchronous :: rdatay
    integer :: ib1 , ib2 , iex , k
    type(mpi_request) :: mrequestx, mrequesty
    type(mpi_status) :: mstatus

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        if ( rcmpi%west /= mpi_proc_null ) &
            sdatax(ib1:ib2) = ml(j1+(iex-1),i1-bb:i2+tb,k)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        if ( rcmpi%east /= mpi_proc_null ) &
            sdatax(ib1:ib2) = ml(j2-(iex-1),i1-bb:i2+tb,k)
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_ineighbor_alltoallv(sdatax, counts, displs, mpi_real8, &
        rdatax, counts, displs, mpi_real8, rcmpi%cartcom, mrequestx)

    call mpi_wait(mrequestx,mstatus)

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
         if ( rcmpi%west /= mpi_proc_null ) &
             ml(j1-iex,i1-bb:i2+tb,k) = rdatax(ib1:ib2)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        if ( rcmpi%east /= mpi_proc_null ) &
            ml(j2+iex,i1-bb:i2+tb,k) = rdatax(ib1:ib2)
      end do
    end do

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        if ( rcmpi%south /= mpi_proc_null ) &
            sdatay(ib1:ib2) = ml(:,i1+(iex-1),k)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        if ( rcmpi%north /= mpi_proc_null ) &
            sdatay(ib1:ib2) = ml(:,i2-(iex-1),k)
      end do
    end do
    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_ineighbor_alltoallv(sdatay, counts, displs, mpi_real8, &
        rdatay, counts, displs, mpi_real8, rcmpi%cartcom, mrequesty)

    call mpi_wait(mrequesty,mstatus)

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        if ( rcmpi%south /= mpi_proc_null ) &
            ml(:,i1-iex,k) = rdatay(ib1:ib2)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        if ( rcmpi%north /= mpi_proc_null ) &
           ml(:,i2+iex,k) = rdatay(ib1:ib2)
      end do
    end do

    end block transmit

  end subroutine real8_3d_exchange

  subroutine real8_4d_exchange(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2 , n1 , n2
    integer :: nx , ny , nk , nn
    integer :: lb , rb , tb , bb
    integer :: ndx , ndy , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    nn = n2-n1+1
    lb = nex
    rb = nex
    tb = 0
    bb = 0
    if ( rcmpi%west == mpi_proc_null ) lb = 1
    if ( rcmpi%east == mpi_proc_null ) rb = 1
    if ( rcmpi%north == mpi_proc_null ) tb = 1
    if ( rcmpi%south == mpi_proc_null ) bb = 1
    tx = ny+bb+tb
    ty = nx+lb+rb
    sizex = nex*tx*nk*nn
    sizey = nex*ty*nk*nn
    ndx = 2*sizex
    ndy = 2*sizey

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk8), dimension(ndx), asynchronous :: sdatax
    real(rk8), dimension(ndx), asynchronous :: rdatax
    real(rk8), dimension(ndy), asynchronous :: sdatay
    real(rk8), dimension(ndy), asynchronous :: rdatay
    integer :: ib1 , ib2 , iex , k , n
    type(mpi_request) :: mrequestx, mrequesty
    type(mpi_status) :: mstatus

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          if ( rcmpi%west /= mpi_proc_null ) &
              sdatax(ib1:ib2) = ml(j1+(iex-1),i1-bb:i2+tb,k,n)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          if ( rcmpi%east /= mpi_proc_null ) &
              sdatax(ib1:ib2) = ml(j2-(iex-1),i1-bb:i2+tb,k,n)
        end do
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_ineighbor_alltoallv(sdatax, counts, displs, mpi_real8, &
        rdatax, counts, displs, mpi_real8, rcmpi%cartcom, mrequestx)

    call mpi_wait(mrequestx,mstatus)

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
           if ( rcmpi%west /= mpi_proc_null ) &
               ml(j1-iex,i1-bb:i2+tb,k,n) = rdatax(ib1:ib2)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          if ( rcmpi%east /= mpi_proc_null ) &
              ml(j2+iex,i1-bb:i2+tb,k,n) = rdatax(ib1:ib2)
        end do
      end do
    end do

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          if ( rcmpi%south /= mpi_proc_null ) &
              sdatay(ib1:ib2) = ml(:,i1+(iex-1),k,n)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          if ( rcmpi%north /= mpi_proc_null ) &
              sdatay(ib1:ib2) = ml(:,i2-(iex-1),k,n)
        end do
      end do
    end do
    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_ineighbor_alltoallv(sdatay, counts, displs, mpi_real8, &
        rdatay, counts, displs, mpi_real8, rcmpi%cartcom, mrequesty)

    call mpi_wait(mrequesty,mstatus)

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          if ( rcmpi%south /= mpi_proc_null ) &
              ml(:,i1-iex,k,n) = rdatay(ib1:ib2)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          if ( rcmpi%north /= mpi_proc_null ) &
             ml(:,i2+iex,k,n) = rdatay(ib1:ib2)
        end do
      end do
    end do

    end block transmit

  end subroutine real8_4d_exchange

  subroutine real4_2d_exchange(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk4) , pointer , dimension(:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2
    integer :: nx , ny
    integer :: lb , rb , tb , bb
    integer :: ndx , ndy , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    lb = nex
    rb = nex
    tb = 0
    bb = 0
    if ( rcmpi%west == mpi_proc_null ) lb = 1
    if ( rcmpi%east == mpi_proc_null ) rb = 1
    if ( rcmpi%north == mpi_proc_null ) tb = 1
    if ( rcmpi%south == mpi_proc_null ) bb = 1
    tx = ny+tb+bb
    ty = nx+lb+rb
    sizex = nex*tx
    sizey = nex*ty
    ndx = 2*sizex
    ndy = 2*sizey

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk4), dimension(ndx), asynchronous :: sdatax
    real(rk4), dimension(ndx), asynchronous :: rdatax
    real(rk4), dimension(ndy), asynchronous :: sdatay
    real(rk4), dimension(ndy), asynchronous :: rdatay
    type(mpi_request) :: mrequestx, mrequesty
    type(mpi_status) :: mstatus
    integer :: ib1 , ib2 , iex

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      if ( rcmpi%west /= mpi_proc_null ) &
          sdatax(ib1:ib2) = ml(j1+(iex-1),i1-bb:i2+tb)
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      if ( rcmpi%east /= mpi_proc_null ) &
          sdatax(ib1:ib2) = ml(j2-(iex-1),i1-bb:i2+tb)
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_ineighbor_alltoallv(sdatax, counts, displs, mpi_real4, &
        rdatax, counts, displs, mpi_real4, rcmpi%cartcom, mrequestx)

    call mpi_wait(mrequestx,mstatus)

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
       if ( rcmpi%west /= mpi_proc_null ) &
           ml(j1-iex,i1-bb:i2+tb) = rdatax(ib1:ib2)
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      if ( rcmpi%east /= mpi_proc_null ) &
          ml(j2+iex,i1-bb:i2+tb) = rdatax(ib1:ib2)
    end do

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      if ( rcmpi%south /= mpi_proc_null ) &
          sdatay(ib1:ib2) = ml(:,i1+(iex-1))
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      if ( rcmpi%north /= mpi_proc_null ) &
          sdatay(ib1:ib2) = ml(:,i2-(iex-1))
    end do
    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_ineighbor_alltoallv(sdatay, counts, displs, mpi_real4, &
        rdatay, counts, displs, mpi_real4, rcmpi%cartcom, mrequesty)

    call mpi_wait(mrequesty,mstatus)

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      if ( rcmpi%south /= mpi_proc_null ) &
          ml(:,i1-iex) = rdatay(ib1:ib2)
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      if ( rcmpi%north /= mpi_proc_null ) &
         ml(:,i2+iex) = rdatay(ib1:ib2)
    end do

    end block transmit

  end subroutine real4_2d_exchange

  subroutine real4_3d_exchange(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk4) , pointer , dimension(:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2
    integer :: nx , ny , nk
    integer :: lb , rb , tb , bb
    integer :: ndx , ndy , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    lb = nex
    rb = nex
    tb = 0
    bb = 0
    if ( rcmpi%west == mpi_proc_null ) lb = 1
    if ( rcmpi%east == mpi_proc_null ) rb = 1
    if ( rcmpi%north == mpi_proc_null ) tb = 1
    if ( rcmpi%south == mpi_proc_null ) bb = 1
    tx = ny+bb+tb
    ty = nx+lb+rb
    sizex = nex*tx*nk
    sizey = nex*ty*nk
    ndx = 2*sizex
    ndy = 2*sizey

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk4), dimension(ndx), asynchronous :: sdatax
    real(rk4), dimension(ndx), asynchronous :: rdatax
    real(rk4), dimension(ndy), asynchronous :: sdatay
    real(rk4), dimension(ndy), asynchronous :: rdatay
    integer :: ib1 , ib2 , iex , k
    type(mpi_request) :: mrequestx, mrequesty
    type(mpi_status) :: mstatus

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        if ( rcmpi%west /= mpi_proc_null ) &
            sdatax(ib1:ib2) = ml(j1+(iex-1),i1-bb:i2+tb,k)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        if ( rcmpi%east /= mpi_proc_null ) &
            sdatax(ib1:ib2) = ml(j2-(iex-1),i1-bb:i2+tb,k)
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_ineighbor_alltoallv(sdatax, counts, displs, mpi_real4, &
        rdatax, counts, displs, mpi_real4, rcmpi%cartcom, mrequestx)

    call mpi_wait(mrequestx,mstatus)

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
         if ( rcmpi%west /= mpi_proc_null ) &
             ml(j1-iex,i1-bb:i2+tb,k) = rdatax(ib1:ib2)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        if ( rcmpi%east /= mpi_proc_null ) &
            ml(j2+iex,i1-bb:i2+tb,k) = rdatax(ib1:ib2)
      end do
    end do

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        if ( rcmpi%south /= mpi_proc_null ) &
            sdatay(ib1:ib2) = ml(:,i1+(iex-1),k)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        if ( rcmpi%north /= mpi_proc_null ) &
            sdatay(ib1:ib2) = ml(:,i2-(iex-1),k)
      end do
    end do
    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_ineighbor_alltoallv(sdatay, counts, displs, mpi_real4, &
        rdatay, counts, displs, mpi_real4, rcmpi%cartcom, mrequesty)

    call mpi_wait(mrequesty,mstatus)

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        if ( rcmpi%south /= mpi_proc_null ) &
            ml(:,i1-iex,k) = rdatay(ib1:ib2)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        if ( rcmpi%north /= mpi_proc_null ) &
           ml(:,i2+iex,k) = rdatay(ib1:ib2)
      end do
    end do

    end block transmit

  end subroutine real4_3d_exchange

  subroutine real4_4d_exchange(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk4) , pointer , dimension(:,:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2 , n1 , n2
    integer :: nx , ny , nk , nn
    integer :: lb , rb , tb , bb
    integer :: ndx , ndy , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    nn = n2-n1+1
    lb = nex
    rb = nex
    tb = 0
    bb = 0
    if ( rcmpi%west == mpi_proc_null ) lb = 1
    if ( rcmpi%east == mpi_proc_null ) rb = 1
    if ( rcmpi%north == mpi_proc_null ) tb = 1
    if ( rcmpi%south == mpi_proc_null ) bb = 1
    tx = ny+bb+tb
    ty = nx+lb+rb
    sizex = nex*tx*nk*nn
    sizey = nex*ty*nk*nn
    ndx = 2*sizex
    ndy = 2*sizey

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk4), dimension(ndx), asynchronous :: sdatax
    real(rk4), dimension(ndx), asynchronous :: rdatax
    real(rk4), dimension(ndy), asynchronous :: sdatay
    real(rk4), dimension(ndy), asynchronous :: rdatay
    integer :: ib1 , ib2 , iex , k , n
    type(mpi_request) :: mrequestx, mrequesty
    type(mpi_status) :: mstatus

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          if ( rcmpi%west /= mpi_proc_null ) &
              sdatax(ib1:ib2) = ml(j1+(iex-1),i1-bb:i2+tb,k,n)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          if ( rcmpi%east /= mpi_proc_null ) &
              sdatax(ib1:ib2) = ml(j2-(iex-1),i1-bb:i2+tb,k,n)
        end do
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_ineighbor_alltoallv(sdatax, counts, displs, mpi_real4, &
        rdatax, counts, displs, mpi_real4, rcmpi%cartcom, mrequestx)

    call mpi_wait(mrequestx,mstatus)

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
           if ( rcmpi%west /= mpi_proc_null ) &
               ml(j1-iex,i1-bb:i2+tb,k,n) = rdatax(ib1:ib2)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          if ( rcmpi%east /= mpi_proc_null ) &
              ml(j2+iex,i1-bb:i2+tb,k,n) = rdatax(ib1:ib2)
        end do
      end do
    end do

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          if ( rcmpi%south /= mpi_proc_null ) &
              sdatay(ib1:ib2) = ml(:,i1+(iex-1),k,n)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          if ( rcmpi%north /= mpi_proc_null ) &
              sdatay(ib1:ib2) = ml(:,i2-(iex-1),k,n)
        end do
      end do
    end do
    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_ineighbor_alltoallv(sdatay, counts, displs, mpi_real4, &
        rdatay, counts, displs, mpi_real4, rcmpi%cartcom, mrequesty)

    call mpi_wait(mrequesty,mstatus)

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          if ( rcmpi%south /= mpi_proc_null ) &
              ml(:,i1-iex,k,n) = rdatay(ib1:ib2)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          if ( rcmpi%north /= mpi_proc_null ) &
             ml(:,i2+iex,k,n) = rdatay(ib1:ib2)
        end do
      end do
    end do

    end block transmit

  end subroutine real4_4d_exchange

  subroutine real8_2d_exchange_lrbt(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk8) , pointer , dimension(:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2
    integer :: ndx , ndy , nx , ny , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    tx = ny
    ty = nx
    sizex = nex*tx
    sizey = nex*ty
    ndx = 2*sizex
    ndy = 2*sizey

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk8), dimension(ndx+ndy), asynchronous :: sdata
    real(rk8), dimension(ndx+ndy), asynchronous :: rdata
    type(mpi_request) :: mrequest
    type(mpi_status) :: mstatus
    integer :: ib1 , ib2 , iex

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      if ( rcmpi%west /= mpi_proc_null ) &
          sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2)
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      if ( rcmpi%east /= mpi_proc_null ) &
          sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2)
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      if ( rcmpi%south /= mpi_proc_null ) &
          sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1))
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      if ( rcmpi%north /= mpi_proc_null ) &
          sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1))
    end do

    counts = [ sizex, sizex, sizey, sizey ]
    displs = [ 0, sizex, 2*sizex, 2*sizex+sizey ]
    call mpi_ineighbor_alltoallv(sdata, counts, displs, mpi_real8, &
        rdata, counts, displs, mpi_real8, rcmpi%cartcom, mrequest)

    call mpi_wait(mrequest,mstatus)

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
       if ( rcmpi%west /= mpi_proc_null ) &
           ml(j1-iex,i1:i2) = rdata(ib1:ib2)
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      if ( rcmpi%east /= mpi_proc_null ) &
          ml(j2+iex,i1:i2) = rdata(ib1:ib2)
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      if ( rcmpi%south /= mpi_proc_null ) &
          ml(j1:j2,i1-iex) = rdata(ib1:ib2)
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      if ( rcmpi%north /= mpi_proc_null ) &
         ml(j1:j2,i2+iex) = rdata(ib1:ib2)
    end do

    end block transmit

  end subroutine real8_2d_exchange_lrbt

  subroutine real8_3d_exchange_lrbt(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2
    integer :: ndx , ndy , nx , ny , nk , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    tx = ny
    ty = nx
    sizex = nex*tx*nk
    sizey = nex*ty*nk
    ndx = 2*sizex
    ndy = 2*sizey

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk8), dimension(ndx+ndy), asynchronous :: sdata
    real(rk8), dimension(ndx+ndy), asynchronous :: rdata
    type(mpi_request) :: mrequest
    type(mpi_status) :: mstatus
    integer :: ib1 , ib2 , iex , k

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        if ( rcmpi%west /= mpi_proc_null ) &
            sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2,k)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        if ( rcmpi%east /= mpi_proc_null ) &
            sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2,k)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        if ( rcmpi%south /= mpi_proc_null ) &
            sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1),k)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        if ( rcmpi%north /= mpi_proc_null ) &
            sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1),k)
      end do
    end do

    counts = [ sizex, sizex, sizey, sizey ]
    displs = [ 0, sizex, 2*sizex, 2*sizex+sizey ]
    call mpi_ineighbor_alltoallv(sdata, counts, displs, mpi_real8, &
        rdata, counts, displs, mpi_real8, rcmpi%cartcom, mrequest)

    call mpi_wait(mrequest,mstatus)

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
         if ( rcmpi%west /= mpi_proc_null ) &
             ml(j1-iex,i1:i2,k) = rdata(ib1:ib2)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        if ( rcmpi%east /= mpi_proc_null ) &
            ml(j2+iex,i1:i2,k) = rdata(ib1:ib2)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        if ( rcmpi%south /= mpi_proc_null ) &
            ml(j1:j2,i1-iex,k) = rdata(ib1:ib2)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        if ( rcmpi%north /= mpi_proc_null ) &
           ml(j1:j2,i2+iex,k) = rdata(ib1:ib2)
      end do
    end do

    end block transmit

  end subroutine real8_3d_exchange_lrbt

  subroutine real8_4d_exchange_lrbt(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2 , n1 , n2
    integer :: ndx , ndy , nx , ny , nk , nn , tx , ty , sizex , sizey

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

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk8), dimension(ndx+ndy), asynchronous :: sdata
    real(rk8), dimension(ndx+ndy), asynchronous :: rdata
    type(mpi_request) :: mrequest
    type(mpi_status) :: mstatus
    integer :: ib1 , ib2 , iex , k , n

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          if ( rcmpi%west /= mpi_proc_null ) &
              sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2,k,n)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          if ( rcmpi%east /= mpi_proc_null ) &
              sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2,k,n)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          if ( rcmpi%south /= mpi_proc_null ) &
              sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1),k,n)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          if ( rcmpi%north /= mpi_proc_null ) &
              sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1),k,n)
        end do
      end do
    end do

    counts = [ sizex, sizex, sizey, sizey ]
    displs = [ 0, sizex, 2*sizex, 2*sizex+sizey ]
    call mpi_ineighbor_alltoallv(sdata, counts, displs, mpi_real8, &
        rdata, counts, displs, mpi_real8, rcmpi%cartcom, mrequest)

    call mpi_wait(mrequest,mstatus)

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
           if ( rcmpi%west /= mpi_proc_null ) &
               ml(j1-iex,i1:i2,k,n) = rdata(ib1:ib2)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          if ( rcmpi%east /= mpi_proc_null ) &
              ml(j2+iex,i1:i2,k,n) = rdata(ib1:ib2)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          if ( rcmpi%south /= mpi_proc_null ) &
              ml(j1:j2,i1-iex,k,n) = rdata(ib1:ib2)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          if ( rcmpi%north /= mpi_proc_null ) &
             ml(j1:j2,i2+iex,k,n) = rdata(ib1:ib2)
        end do
      end do
    end do

    end block transmit

  end subroutine real8_4d_exchange_lrbt

  subroutine real4_2d_exchange_lrbt(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk4) , pointer , dimension(:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2
    integer :: ndx , ndy , nx , ny , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    tx = ny
    ty = nx
    sizex = nex*tx
    sizey = nex*ty
    ndx = 2*sizex
    ndy = 2*sizey

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk4), dimension(ndx+ndy), asynchronous :: sdata
    real(rk4), dimension(ndx+ndy), asynchronous :: rdata
    type(mpi_request) :: mrequest
    type(mpi_status) :: mstatus
    integer :: ib1 , ib2 , iex

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      if ( rcmpi%west /= mpi_proc_null ) &
          sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2)
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      if ( rcmpi%east /= mpi_proc_null ) &
          sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2)
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      if ( rcmpi%south /= mpi_proc_null ) &
          sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1))
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      if ( rcmpi%north /= mpi_proc_null ) &
          sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1))
    end do

    counts = [ sizex, sizex, sizey, sizey ]
    displs = [ 0, sizex, 2*sizex, 2*sizex+sizey ]
    call mpi_ineighbor_alltoallv(sdata, counts, displs, mpi_real4, &
        rdata, counts, displs, mpi_real4, rcmpi%cartcom, mrequest)

    call mpi_wait(mrequest,mstatus)

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
       if ( rcmpi%west /= mpi_proc_null ) &
           ml(j1-iex,i1:i2) = rdata(ib1:ib2)
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      if ( rcmpi%east /= mpi_proc_null ) &
          ml(j2+iex,i1:i2) = rdata(ib1:ib2)
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      if ( rcmpi%south /= mpi_proc_null ) &
          ml(j1:j2,i1-iex) = rdata(ib1:ib2)
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      if ( rcmpi%north /= mpi_proc_null ) &
         ml(j1:j2,i2+iex) = rdata(ib1:ib2)
    end do

    end block transmit

  end subroutine real4_2d_exchange_lrbt

  subroutine real4_3d_exchange_lrbt(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk4) , pointer , dimension(:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2
    integer :: ndx , ndy , nx , ny , nk , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    tx = ny
    ty = nx
    sizex = nex*tx*nk
    sizey = nex*ty*nk
    ndx = 2*sizex
    ndy = 2*sizey

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk4), dimension(ndx+ndy), asynchronous :: sdata
    real(rk4), dimension(ndx+ndy), asynchronous :: rdata
    type(mpi_request) :: mrequest
    type(mpi_status) :: mstatus
    integer :: ib1 , ib2 , iex , k

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        if ( rcmpi%west /= mpi_proc_null ) &
            sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2,k)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        if ( rcmpi%east /= mpi_proc_null ) &
            sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2,k)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        if ( rcmpi%south /= mpi_proc_null ) &
            sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1),k)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        if ( rcmpi%north /= mpi_proc_null ) &
            sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1),k)
      end do
    end do

    counts = [ sizex, sizex, sizey, sizey ]
    displs = [ 0, sizex, 2*sizex, 2*sizex+sizey ]
    call mpi_ineighbor_alltoallv(sdata, counts, displs, mpi_real4, &
        rdata, counts, displs, mpi_real4, rcmpi%cartcom, mrequest)

    call mpi_wait(mrequest,mstatus)

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
         if ( rcmpi%west /= mpi_proc_null ) &
             ml(j1-iex,i1:i2,k) = rdata(ib1:ib2)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        if ( rcmpi%east /= mpi_proc_null ) &
            ml(j2+iex,i1:i2,k) = rdata(ib1:ib2)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        if ( rcmpi%south /= mpi_proc_null ) &
            ml(j1:j2,i1-iex,k) = rdata(ib1:ib2)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        if ( rcmpi%north /= mpi_proc_null ) &
           ml(j1:j2,i2+iex,k) = rdata(ib1:ib2)
      end do
    end do

    end block transmit

  end subroutine real4_3d_exchange_lrbt

  subroutine real4_4d_exchange_lrbt(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk4) , pointer , dimension(:,:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2 , n1 , n2
    integer :: ndx , ndy , nx , ny , nk , nn , tx , ty , sizex , sizey

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

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk4), dimension(ndx+ndy), asynchronous :: sdata
    real(rk4), dimension(ndx+ndy), asynchronous :: rdata
    type(mpi_request) :: mrequest
    type(mpi_status) :: mstatus
    integer :: ib1 , ib2 , iex , k , n

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          if ( rcmpi%west /= mpi_proc_null ) &
              sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2,k,n)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          if ( rcmpi%east /= mpi_proc_null ) &
              sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2,k,n)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          if ( rcmpi%south /= mpi_proc_null ) &
              sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1),k,n)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          if ( rcmpi%north /= mpi_proc_null ) &
              sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1),k,n)
        end do
      end do
    end do

    counts = [ sizex, sizex, sizey, sizey ]
    displs = [ 0, sizex, 2*sizex, 2*sizex+sizey ]
    call mpi_ineighbor_alltoallv(sdata, counts, displs, mpi_real4, &
        rdata, counts, displs, mpi_real4, rcmpi%cartcom, mrequest)

    call mpi_wait(mrequest,mstatus)

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
           if ( rcmpi%west /= mpi_proc_null ) &
               ml(j1-iex,i1:i2,k,n) = rdata(ib1:ib2)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          if ( rcmpi%east /= mpi_proc_null ) &
              ml(j2+iex,i1:i2,k,n) = rdata(ib1:ib2)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          if ( rcmpi%south /= mpi_proc_null ) &
              ml(j1:j2,i1-iex,k,n) = rdata(ib1:ib2)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          if ( rcmpi%north /= mpi_proc_null ) &
             ml(j1:j2,i2+iex,k,n) = rdata(ib1:ib2)
        end do
      end do
    end do

    end block transmit

  end subroutine real4_4d_exchange_lrbt

  subroutine real8_2d_exchange_lr(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk8) , pointer , dimension(:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2
    integer :: ndx , ny , tx , sizex

    ny = i2-i1+1
    tx = ny
    sizex = nex*tx
    ndx = 2*sizex

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk8), dimension(ndx), asynchronous :: sdata
    real(rk8), dimension(ndx), asynchronous :: rdata
    type(mpi_request) :: mrequest
    type(mpi_status) :: mstatus
    integer :: ib1 , ib2 , iex

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      if ( rcmpi%west /= mpi_proc_null ) &
          sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2)
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      if ( rcmpi%east /= mpi_proc_null ) &
          sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2)
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_ineighbor_alltoallv(sdata, counts, displs, mpi_real8, &
        rdata, counts, displs, mpi_real8, rcmpi%cartcom, mrequest)

    call mpi_wait(mrequest,mstatus)

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
       if ( rcmpi%west /= mpi_proc_null ) &
           ml(j1-iex,i1:i2) = rdata(ib1:ib2)
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      if ( rcmpi%east /= mpi_proc_null ) &
          ml(j2+iex,i1:i2) = rdata(ib1:ib2)
    end do

    end block transmit

  end subroutine real8_2d_exchange_lr

  subroutine real8_3d_exchange_lr(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2
    integer :: ndx , ny , nk , tx , sizex

    ny = i2-i1+1
    nk = k2-k1+1
    tx = ny
    sizex = nex*tx*nk
    ndx = 2*sizex

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk8), dimension(ndx), asynchronous :: sdata
    real(rk8), dimension(ndx), asynchronous :: rdata
    type(mpi_request) :: mrequest
    type(mpi_status) :: mstatus
    integer :: ib1 , ib2 , iex , k

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        if ( rcmpi%west /= mpi_proc_null ) &
            sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2,k)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        if ( rcmpi%east /= mpi_proc_null ) &
            sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2,k)
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_ineighbor_alltoallv(sdata, counts, displs, mpi_real8, &
        rdata, counts, displs, mpi_real8, rcmpi%cartcom, mrequest)

    call mpi_wait(mrequest,mstatus)

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
         if ( rcmpi%west /= mpi_proc_null ) &
             ml(j1-iex,i1:i2,k) = rdata(ib1:ib2)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        if ( rcmpi%east /= mpi_proc_null ) &
            ml(j2+iex,i1:i2,k) = rdata(ib1:ib2)
      end do
    end do

    end block transmit

  end subroutine real8_3d_exchange_lr

  subroutine real8_4d_exchange_lr(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2 , n1 , n2
    integer :: ndx , ny , nk , nn , tx , sizex

    ny = i2-i1+1
    nk = k2-k1+1
    nn = n2-n1+1
    tx = ny
    sizex = nex*tx*nk*nn
    ndx = 2*sizex

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk8), dimension(ndx), asynchronous :: sdata
    real(rk8), dimension(ndx), asynchronous :: rdata
    type(mpi_request) :: mrequest
    type(mpi_status) :: mstatus
    integer :: ib1 , ib2 , iex , k , n

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          if ( rcmpi%west /= mpi_proc_null ) &
              sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2,k,n)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          if ( rcmpi%east /= mpi_proc_null ) &
              sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2,k,n)
        end do
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_ineighbor_alltoallv(sdata, counts, displs, mpi_real8, &
        rdata, counts, displs, mpi_real8, rcmpi%cartcom, mrequest)

    call mpi_wait(mrequest,mstatus)

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
           if ( rcmpi%west /= mpi_proc_null ) &
               ml(j1-iex,i1:i2,k,n) = rdata(ib1:ib2)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          if ( rcmpi%east /= mpi_proc_null ) &
              ml(j2+iex,i1:i2,k,n) = rdata(ib1:ib2)
        end do
      end do
    end do

    end block transmit

  end subroutine real8_4d_exchange_lr

  subroutine real4_2d_exchange_lr(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk4) , pointer , dimension(:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2
    integer :: ndx , ny , tx , sizex

    ny = i2-i1+1
    tx = ny
    sizex = nex*tx
    ndx = 2*sizex

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk4), dimension(ndx), asynchronous :: sdata
    real(rk4), dimension(ndx), asynchronous :: rdata
    type(mpi_request) :: mrequest
    type(mpi_status) :: mstatus
    integer :: ib1 , ib2 , iex

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      if ( rcmpi%west /= mpi_proc_null ) &
          sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2)
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      if ( rcmpi%east /= mpi_proc_null ) &
          sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2)
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_ineighbor_alltoallv(sdata, counts, displs, mpi_real4, &
        rdata, counts, displs, mpi_real4, rcmpi%cartcom, mrequest)

    call mpi_wait(mrequest,mstatus)

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
       if ( rcmpi%west /= mpi_proc_null ) &
           ml(j1-iex,i1:i2) = rdata(ib1:ib2)
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      if ( rcmpi%east /= mpi_proc_null ) &
          ml(j2+iex,i1:i2) = rdata(ib1:ib2)
    end do

    end block transmit

  end subroutine real4_2d_exchange_lr

  subroutine real4_3d_exchange_lr(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk4) , pointer , dimension(:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2
    integer :: ndx , ny , nk , tx , sizex

    ny = i2-i1+1
    nk = k2-k1+1
    tx = ny
    sizex = nex*tx*nk
    ndx = 2*sizex

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk4), dimension(ndx), asynchronous :: sdata
    real(rk4), dimension(ndx), asynchronous :: rdata
    type(mpi_request) :: mrequest
    type(mpi_status) :: mstatus
    integer :: ib1 , ib2 , iex , k

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        if ( rcmpi%west /= mpi_proc_null ) &
            sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2,k)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        if ( rcmpi%east /= mpi_proc_null ) &
            sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2,k)
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_ineighbor_alltoallv(sdata, counts, displs, mpi_real4, &
        rdata, counts, displs, mpi_real4, rcmpi%cartcom, mrequest)

    call mpi_wait(mrequest,mstatus)

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
         if ( rcmpi%west /= mpi_proc_null ) &
             ml(j1-iex,i1:i2,k) = rdata(ib1:ib2)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        if ( rcmpi%east /= mpi_proc_null ) &
            ml(j2+iex,i1:i2,k) = rdata(ib1:ib2)
      end do
    end do

    end block transmit

  end subroutine real4_3d_exchange_lr

  subroutine real4_4d_exchange_lr(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk4) , pointer , dimension(:,:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2 , n1 , n2
    integer :: ndx , ny , nk , nn , tx , sizex

    ny = i2-i1+1
    nk = k2-k1+1
    nn = n2-n1+1
    tx = ny
    sizex = nex*tx*nk*nn
    ndx = 2*sizex

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk4), dimension(ndx), asynchronous :: sdata
    real(rk4), dimension(ndx), asynchronous :: rdata
    type(mpi_request) :: mrequest
    type(mpi_status) :: mstatus
    integer :: ib1 , ib2 , iex , k , n

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          if ( rcmpi%west /= mpi_proc_null ) &
              sdata(ib1:ib2) = ml(j1+(iex-1),i1:i2,k,n)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          if ( rcmpi%east /= mpi_proc_null ) &
              sdata(ib1:ib2) = ml(j2-(iex-1),i1:i2,k,n)
        end do
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_ineighbor_alltoallv(sdata, counts, displs, mpi_real4, &
        rdata, counts, displs, mpi_real4, rcmpi%cartcom, mrequest)

    call mpi_wait(mrequest,mstatus)

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
           if ( rcmpi%west /= mpi_proc_null ) &
               ml(j1-iex,i1:i2,k,n) = rdata(ib1:ib2)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          if ( rcmpi%east /= mpi_proc_null ) &
              ml(j2+iex,i1:i2,k,n) = rdata(ib1:ib2)
        end do
      end do
    end do

    end block transmit

  end subroutine real4_4d_exchange_lr

  subroutine real8_2d_exchange_bt(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk8) , pointer , dimension(:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2
    integer :: ndy , nx , ty , sizey

    nx = j2-j1+1
    ty = nx
    sizey = nex*ty
    ndy = 2*sizey

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk8), dimension(ndy), asynchronous :: sdata
    real(rk8), dimension(ndy), asynchronous :: rdata
    type(mpi_request) :: mrequest
    type(mpi_status) :: mstatus
    integer :: ib1 , ib2 , iex

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      if ( rcmpi%south /= mpi_proc_null ) &
          sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1))
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      if ( rcmpi%north /= mpi_proc_null ) &
          sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1))
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_ineighbor_alltoallv(sdata, counts, displs, mpi_real8, &
        rdata, counts, displs, mpi_real8, rcmpi%cartcom, mrequest)

    call mpi_wait(mrequest,mstatus)

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      if ( rcmpi%south /= mpi_proc_null ) &
          ml(j1:j2,i1-iex) = rdata(ib1:ib2)
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      if ( rcmpi%north /= mpi_proc_null ) &
         ml(j1:j2,i2+iex) = rdata(ib1:ib2)
    end do

    end block transmit

  end subroutine real8_2d_exchange_bt

  subroutine real8_3d_exchange_bt(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2
    integer :: ndy , nx , nk , ty , sizey

    nx = j2-j1+1
    nk = k2-k1+1
    ty = nx
    sizey = nex*ty*nk
    ndy = 2*sizey

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk8), dimension(ndy), asynchronous :: sdata
    real(rk8), dimension(ndy), asynchronous :: rdata
    type(mpi_request) :: mrequest
    type(mpi_status) :: mstatus
    integer :: ib1 , ib2 , iex , k

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        if ( rcmpi%south /= mpi_proc_null ) &
            sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1),k)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        if ( rcmpi%north /= mpi_proc_null ) &
            sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1),k)
      end do
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_ineighbor_alltoallv(sdata, counts, displs, mpi_real8, &
        rdata, counts, displs, mpi_real8, rcmpi%cartcom, mrequest)

    call mpi_wait(mrequest,mstatus)

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        if ( rcmpi%south /= mpi_proc_null ) &
            ml(j1:j2,i1-iex,k) = rdata(ib1:ib2)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        if ( rcmpi%north /= mpi_proc_null ) &
           ml(j1:j2,i2+iex,k) = rdata(ib1:ib2)
      end do
    end do

    end block transmit

  end subroutine real8_3d_exchange_bt

  subroutine real8_4d_exchange_bt(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2 , n1 , n2
    integer :: ndy , nx , nk , nn , ty , sizey

    nx = j2-j1+1
    nk = k2-k1+1
    nn = n2-n1+1
    ty = nx
    sizey = nex*ty*nk*nn
    ndy = 2*sizey

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk8), dimension(ndy), asynchronous :: sdata
    real(rk8), dimension(ndy), asynchronous :: rdata
    type(mpi_request) :: mrequest
    type(mpi_status) :: mstatus
    integer :: ib1 , ib2 , iex , k , n

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          if ( rcmpi%south /= mpi_proc_null ) &
              sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1),k,n)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          if ( rcmpi%north /= mpi_proc_null ) &
              sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1),k,n)
        end do
      end do
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_ineighbor_alltoallv(sdata, counts, displs, mpi_real8, &
        rdata, counts, displs, mpi_real8, rcmpi%cartcom, mrequest)

    call mpi_wait(mrequest,mstatus)

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          if ( rcmpi%south /= mpi_proc_null ) &
              ml(j1:j2,i1-iex,k,n) = rdata(ib1:ib2)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          if ( rcmpi%north /= mpi_proc_null ) &
             ml(j1:j2,i2+iex,k,n) = rdata(ib1:ib2)
        end do
      end do
    end do

    end block transmit

  end subroutine real8_4d_exchange_bt

  subroutine real4_2d_exchange_bt(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk4) , pointer , dimension(:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2
    integer :: ndy , nx , ty , sizey

    nx = j2-j1+1
    ty = nx
    sizey = nex*ty
    ndy = 2*sizey

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk4), dimension(ndy), asynchronous :: sdata
    real(rk4), dimension(ndy), asynchronous :: rdata
    type(mpi_request) :: mrequest
    type(mpi_status) :: mstatus
    integer :: ib1 , ib2 , iex

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      if ( rcmpi%south /= mpi_proc_null ) &
          sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1))
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      if ( rcmpi%north /= mpi_proc_null ) &
          sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1))
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_ineighbor_alltoallv(sdata, counts, displs, mpi_real4, &
        rdata, counts, displs, mpi_real4, rcmpi%cartcom, mrequest)

    call mpi_wait(mrequest,mstatus)

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      if ( rcmpi%south /= mpi_proc_null ) &
          ml(j1:j2,i1-iex) = rdata(ib1:ib2)
    end do
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      if ( rcmpi%north /= mpi_proc_null ) &
         ml(j1:j2,i2+iex) = rdata(ib1:ib2)
    end do

    end block transmit

  end subroutine real4_2d_exchange_bt

  subroutine real4_3d_exchange_bt(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk4) , pointer , dimension(:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2
    integer :: ndy , nx , nk , ty , sizey

    nx = j2-j1+1
    nk = k2-k1+1
    ty = nx
    sizey = nex*ty*nk
    ndy = 2*sizey

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk4), dimension(ndy), asynchronous :: sdata
    real(rk4), dimension(ndy), asynchronous :: rdata
    type(mpi_request) :: mrequest
    type(mpi_status) :: mstatus
    integer :: ib1 , ib2 , iex , k

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        if ( rcmpi%south /= mpi_proc_null ) &
            sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1),k)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        if ( rcmpi%north /= mpi_proc_null ) &
            sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1),k)
      end do
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_ineighbor_alltoallv(sdata, counts, displs, mpi_real4, &
        rdata, counts, displs, mpi_real4, rcmpi%cartcom, mrequest)

    call mpi_wait(mrequest,mstatus)

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        if ( rcmpi%south /= mpi_proc_null ) &
            ml(j1:j2,i1-iex,k) = rdata(ib1:ib2)
      end do
    end do
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        if ( rcmpi%north /= mpi_proc_null ) &
           ml(j1:j2,i2+iex,k) = rdata(ib1:ib2)
      end do
    end do

    end block transmit

  end subroutine real4_3d_exchange_bt

  subroutine real4_4d_exchange_bt(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk4) , pointer , dimension(:,:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2 , n1 , n2
    integer :: ndy , nx , nk , nn , ty , sizey

    nx = j2-j1+1
    nk = k2-k1+1
    nn = n2-n1+1
    ty = nx
    sizey = nex*ty*nk*nn
    ndy = 2*sizey

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk4), dimension(ndy), asynchronous :: sdata
    real(rk4), dimension(ndy), asynchronous :: rdata
    type(mpi_request) :: mrequest
    type(mpi_status) :: mstatus
    integer :: ib1 , ib2 , iex , k , n

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          if ( rcmpi%south /= mpi_proc_null ) &
              sdata(ib1:ib2) = ml(j1:j2,i1+(iex-1),k,n)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          if ( rcmpi%north /= mpi_proc_null ) &
              sdata(ib1:ib2) = ml(j1:j2,i2-(iex-1),k,n)
        end do
      end do
    end do

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_ineighbor_alltoallv(sdata, counts, displs, mpi_real4, &
        rdata, counts, displs, mpi_real4, rcmpi%cartcom, mrequest)

    call mpi_wait(mrequest,mstatus)

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          if ( rcmpi%south /= mpi_proc_null ) &
              ml(j1:j2,i1-iex,k,n) = rdata(ib1:ib2)
        end do
      end do
    end do
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          if ( rcmpi%north /= mpi_proc_null ) &
             ml(j1:j2,i2+iex,k,n) = rdata(ib1:ib2)
        end do
      end do
    end do

    end block transmit

  end subroutine real4_4d_exchange_bt

  subroutine real8_2d_exchange_lb(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk8) , pointer , dimension(:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2
    integer :: nx , ny
    integer :: lb , bb
    integer :: ndx , ndy , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    lb = nex
    bb = 0
    if ( rcmpi%west == mpi_proc_null ) lb = 1
    if ( rcmpi%south == mpi_proc_null ) bb = 1
    tx = ny+bb
    ty = nx+lb
    sizex = nex*tx
    sizey = nex*ty
    ndx = 2*sizex
    ndy = 2*sizey

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk8), dimension(ndx), asynchronous :: sdatax
    real(rk8), dimension(ndx), asynchronous :: rdatax
    real(rk8), dimension(ndy), asynchronous :: sdatay
    real(rk8), dimension(ndy), asynchronous :: rdatay
    type(mpi_request) :: mrequestx, mrequesty
    type(mpi_status) :: mstatus
    integer :: ib1 , ib2 , iex

    ib2 = sizex
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      if ( rcmpi%east /= mpi_proc_null ) &
          sdatax(ib1:ib2) = ml(j2-(iex-1),i1-bb:i2)
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_ineighbor_alltoallv(sdatax, counts, displs, mpi_real8, &
        rdatax, counts, displs, mpi_real8, rcmpi%cartcom, mrequestx)

    call mpi_wait(mrequestx,mstatus)

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
       if ( rcmpi%west /= mpi_proc_null ) &
           ml(j1-iex,i1-bb:i2) = rdatax(ib1:ib2)
    end do

    ib2 = sizey
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      if ( rcmpi%north /= mpi_proc_null ) &
          sdatay(ib1:ib2) = ml(:j2,i2-(iex-1))
    end do
    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_ineighbor_alltoallv(sdatay, counts, displs, mpi_real8, &
        rdatay, counts, displs, mpi_real8, rcmpi%cartcom, mrequesty)

    call mpi_wait(mrequesty,mstatus)

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      if ( rcmpi%south /= mpi_proc_null ) &
          ml(:j2,i1-iex) = rdatay(ib1:ib2)
    end do

    end block transmit

  end subroutine real8_2d_exchange_lb

  subroutine real8_3d_exchange_lb(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2
    integer :: nx , ny , nk
    integer :: lb , bb
    integer :: ndx , ndy , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    lb = nex
    bb = 0
    if ( rcmpi%west == mpi_proc_null ) lb = 1
    if ( rcmpi%south == mpi_proc_null ) bb = 1
    tx = ny+bb
    ty = nx+lb
    sizex = nex*tx*nk
    sizey = nex*ty*nk
    ndx = 2*sizex
    ndy = 2*sizey

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk8), dimension(ndx), asynchronous :: sdatax
    real(rk8), dimension(ndx), asynchronous :: rdatax
    real(rk8), dimension(ndy), asynchronous :: sdatay
    real(rk8), dimension(ndy), asynchronous :: rdatay
    type(mpi_request) :: mrequestx, mrequesty
    type(mpi_status) :: mstatus
    integer :: ib1 , ib2 , iex , k

    ib2 = sizex
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        if ( rcmpi%east /= mpi_proc_null ) &
            sdatax(ib1:ib2) = ml(j2-(iex-1),i1-bb:i2,k)
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_ineighbor_alltoallv(sdatax, counts, displs, mpi_real8, &
        rdatax, counts, displs, mpi_real8, rcmpi%cartcom, mrequestx)

    call mpi_wait(mrequestx,mstatus)

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
         if ( rcmpi%west /= mpi_proc_null ) &
             ml(j1-iex,i1-bb:i2,k) = rdatax(ib1:ib2)
      end do
    end do

    ib2 = sizey
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        if ( rcmpi%north /= mpi_proc_null ) &
            sdatay(ib1:ib2) = ml(:j2,i2-(iex-1),k)
      end do
    end do
    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_ineighbor_alltoallv(sdatay, counts, displs, mpi_real8, &
        rdatay, counts, displs, mpi_real8, rcmpi%cartcom, mrequesty)

    call mpi_wait(mrequesty,mstatus)

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        if ( rcmpi%south /= mpi_proc_null ) &
            ml(:j2,i1-iex,k) = rdatay(ib1:ib2)
      end do
    end do

    end block transmit

  end subroutine real8_3d_exchange_lb

  subroutine real8_4d_exchange_lb(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2 , n1 , n2
    integer :: nx , ny , nk , nn
    integer :: lb , bb
    integer :: ndx , ndy , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    nn = n2-n1+1
    lb = nex
    bb = 0
    if ( rcmpi%west == mpi_proc_null ) lb = 1
    if ( rcmpi%south == mpi_proc_null ) bb = 1
    tx = ny+bb
    ty = nx+lb
    sizex = nex*tx*nk
    sizey = nex*ty*nk
    ndx = 2*sizex
    ndy = 2*sizey

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk8), dimension(ndx), asynchronous :: sdatax
    real(rk8), dimension(ndx), asynchronous :: rdatax
    real(rk8), dimension(ndy), asynchronous :: sdatay
    real(rk8), dimension(ndy), asynchronous :: rdatay
    type(mpi_request) :: mrequestx, mrequesty
    type(mpi_status) :: mstatus
    integer :: ib1 , ib2 , iex , k , n

    ib2 = sizex
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          if ( rcmpi%east /= mpi_proc_null ) &
              sdatax(ib1:ib2) = ml(j2-(iex-1),i1-bb:i2,k,n)
        end do
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_ineighbor_alltoallv(sdatax, counts, displs, mpi_real8, &
        rdatax, counts, displs, mpi_real8, rcmpi%cartcom, mrequestx)

    call mpi_wait(mrequestx,mstatus)

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          if ( rcmpi%west /= mpi_proc_null ) &
              ml(j1-iex,i1-bb:i2,k,n) = rdatax(ib1:ib2)
        end do
      end do
    end do

    ib2 = sizey
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          if ( rcmpi%north /= mpi_proc_null ) &
              sdatay(ib1:ib2) = ml(:j2,i2-(iex-1),k,n)
        end do
      end do
    end do
    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_ineighbor_alltoallv(sdatay, counts, displs, mpi_real8, &
        rdatay, counts, displs, mpi_real8, rcmpi%cartcom, mrequesty)

    call mpi_wait(mrequesty,mstatus)

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          if ( rcmpi%south /= mpi_proc_null ) &
              ml(:j2,i1-iex,k,n) = rdatay(ib1:ib2)
        end do
      end do
    end do

    end block transmit

  end subroutine real8_4d_exchange_lb

  subroutine real4_2d_exchange_lb(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk4) , pointer , dimension(:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2
    integer :: nx , ny
    integer :: lb , bb
    integer :: ndx , ndy , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    lb = nex
    bb = 0
    if ( rcmpi%west == mpi_proc_null ) lb = 1
    if ( rcmpi%south == mpi_proc_null ) bb = 1
    tx = ny+bb
    ty = nx+lb
    sizex = nex*tx
    sizey = nex*ty
    ndx = 2*sizex
    ndy = 2*sizey

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk4), dimension(ndx), asynchronous :: sdatax
    real(rk4), dimension(ndx), asynchronous :: rdatax
    real(rk4), dimension(ndy), asynchronous :: sdatay
    real(rk4), dimension(ndy), asynchronous :: rdatay
    type(mpi_request) :: mrequestx, mrequesty
    type(mpi_status) :: mstatus
    integer :: ib1 , ib2 , iex

    ib2 = sizex
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      if ( rcmpi%east /= mpi_proc_null ) &
          sdatax(ib1:ib2) = ml(j2-(iex-1),i1-bb:i2)
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_ineighbor_alltoallv(sdatax, counts, displs, mpi_real4, &
        rdatax, counts, displs, mpi_real4, rcmpi%cartcom, mrequestx)

    call mpi_wait(mrequestx,mstatus)

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
       if ( rcmpi%west /= mpi_proc_null ) &
           ml(j1-iex,i1-bb:i2) = rdatax(ib1:ib2)
    end do

    ib2 = sizey
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      if ( rcmpi%north /= mpi_proc_null ) &
          sdatay(ib1:ib2) = ml(:j2,i2-(iex-1))
    end do
    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_ineighbor_alltoallv(sdatay, counts, displs, mpi_real4, &
        rdatay, counts, displs, mpi_real4, rcmpi%cartcom, mrequesty)

    call mpi_wait(mrequesty,mstatus)

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      if ( rcmpi%south /= mpi_proc_null ) &
          ml(:j2,i1-iex) = rdatay(ib1:ib2)
    end do

    end block transmit

  end subroutine real4_2d_exchange_lb

  subroutine real4_3d_exchange_lb(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk4) , pointer , dimension(:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2
    integer :: nx , ny , nk
    integer :: lb , bb
    integer :: ndx , ndy , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    lb = nex
    bb = 0
    if ( rcmpi%west == mpi_proc_null ) lb = 1
    if ( rcmpi%south == mpi_proc_null ) bb = 1
    tx = ny+bb
    ty = nx+lb
    sizex = nex*tx*nk
    sizey = nex*ty*nk
    ndx = 2*sizex
    ndy = 2*sizey

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk4), dimension(ndx), asynchronous :: sdatax
    real(rk4), dimension(ndx), asynchronous :: rdatax
    real(rk4), dimension(ndy), asynchronous :: sdatay
    real(rk4), dimension(ndy), asynchronous :: rdatay
    type(mpi_request) :: mrequestx, mrequesty
    type(mpi_status) :: mstatus
    integer :: ib1 , ib2 , iex , k

    ib2 = sizex
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        if ( rcmpi%east /= mpi_proc_null ) &
            sdatax(ib1:ib2) = ml(j2-(iex-1),i1-bb:i2,k)
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_ineighbor_alltoallv(sdatax, counts, displs, mpi_real4, &
        rdatax, counts, displs, mpi_real4, rcmpi%cartcom, mrequestx)

    call mpi_wait(mrequestx,mstatus)

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
         if ( rcmpi%west /= mpi_proc_null ) &
             ml(j1-iex,i1-bb:i2,k) = rdatax(ib1:ib2)
      end do
    end do

    ib2 = sizey
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        if ( rcmpi%north /= mpi_proc_null ) &
            sdatay(ib1:ib2) = ml(:j2,i2-(iex-1),k)
      end do
    end do
    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_ineighbor_alltoallv(sdatay, counts, displs, mpi_real4, &
        rdatay, counts, displs, mpi_real4, rcmpi%cartcom, mrequesty)

    call mpi_wait(mrequesty,mstatus)

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        if ( rcmpi%south /= mpi_proc_null ) &
            ml(:j2,i1-iex,k) = rdatay(ib1:ib2)
      end do
    end do

    end block transmit

  end subroutine real4_3d_exchange_lb

  subroutine real4_4d_exchange_lb(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk4) , pointer , dimension(:,:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2 , n1 , n2
    integer :: nx , ny , nk , nn
    integer :: lb , bb
    integer :: ndx , ndy , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    nn = n2-n1+1
    lb = nex
    bb = 0
    if ( rcmpi%west == mpi_proc_null ) lb = 1
    if ( rcmpi%south == mpi_proc_null ) bb = 1
    tx = ny+bb
    ty = nx+lb
    sizex = nex*tx*nk
    sizey = nex*ty*nk
    ndx = 2*sizex
    ndy = 2*sizey

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk4), dimension(ndx), asynchronous :: sdatax
    real(rk4), dimension(ndx), asynchronous :: rdatax
    real(rk4), dimension(ndy), asynchronous :: sdatay
    real(rk4), dimension(ndy), asynchronous :: rdatay
    type(mpi_request) :: mrequestx, mrequesty
    type(mpi_status) :: mstatus
    integer :: ib1 , ib2 , iex , k , n

    ib2 = sizex
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          if ( rcmpi%east /= mpi_proc_null ) &
              sdatax(ib1:ib2) = ml(j2-(iex-1),i1-bb:i2,k,n)
        end do
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_ineighbor_alltoallv(sdatax, counts, displs, mpi_real4, &
        rdatax, counts, displs, mpi_real4, rcmpi%cartcom, mrequestx)

    call mpi_wait(mrequestx,mstatus)

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          if ( rcmpi%west /= mpi_proc_null ) &
              ml(j1-iex,i1-bb:i2,k,n) = rdatax(ib1:ib2)
        end do
      end do
    end do

    ib2 = sizey
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          if ( rcmpi%north /= mpi_proc_null ) &
              sdatay(ib1:ib2) = ml(:j2,i2-(iex-1),k,n)
        end do
      end do
    end do
    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_ineighbor_alltoallv(sdatay, counts, displs, mpi_real4, &
        rdatay, counts, displs, mpi_real4, rcmpi%cartcom, mrequesty)

    call mpi_wait(mrequesty,mstatus)

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          if ( rcmpi%south /= mpi_proc_null ) &
              ml(:j2,i1-iex,k,n) = rdatay(ib1:ib2)
        end do
      end do
    end do

    end block transmit

  end subroutine real4_4d_exchange_lb

  subroutine real8_2d_exchange_rt(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk8) , pointer , dimension(:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2
    integer :: nx , ny
    integer :: rb , tb
    integer :: ndx , ndy , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    rb = nex
    tb = 0
    if ( rcmpi%east == mpi_proc_null ) rb = 1
    if ( rcmpi%north == mpi_proc_null ) tb = 1
    tx = ny+tb
    ty = nx+rb
    sizex = nex*tx
    sizey = nex*ty
    ndx = 2*sizex
    ndy = 2*sizey

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk8), dimension(ndx), asynchronous :: sdatax
    real(rk8), dimension(ndx), asynchronous :: rdatax
    real(rk8), dimension(ndy), asynchronous :: sdatay
    real(rk8), dimension(ndy), asynchronous :: rdatay
    type(mpi_request) :: mrequestx, mrequesty
    type(mpi_status) :: mstatus
    integer :: ib1 , ib2 , iex

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      if ( rcmpi%west /= mpi_proc_null ) &
          sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2+tb)
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_ineighbor_alltoallv(sdatax, counts, displs, mpi_real8, &
        rdatax, counts, displs, mpi_real8, rcmpi%cartcom, mrequestx)

    call mpi_wait(mrequestx,mstatus)

    ib2 = sizex
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      if ( rcmpi%east /= mpi_proc_null ) &
          ml(j2+iex,i1:i2+tb) = rdatax(ib1:ib2)
    end do

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      if ( rcmpi%south /= mpi_proc_null ) &
          sdatay(ib1:ib2) = ml(j1:,i1+(iex-1))
    end do
    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_ineighbor_alltoallv(sdatay, counts, displs, mpi_real8, &
        rdatay, counts, displs, mpi_real8, rcmpi%cartcom, mrequesty)

    call mpi_wait(mrequesty,mstatus)

    ib2 = sizey
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      if ( rcmpi%north /= mpi_proc_null ) &
         ml(j1:,i2+iex) = rdatay(ib1:ib2)
    end do

    end block transmit

  end subroutine real8_2d_exchange_rt

  subroutine real8_3d_exchange_rt(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2
    integer :: nx , ny , nk
    integer :: rb , tb
    integer :: ndx , ndy , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    rb = nex
    tb = 0
    if ( rcmpi%east == mpi_proc_null ) rb = 1
    if ( rcmpi%north == mpi_proc_null ) tb = 1
    tx = ny+tb
    ty = nx+rb
    sizex = nex*tx*nk
    sizey = nex*ty*nk
    ndx = 2*sizex
    ndy = 2*sizey

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk8), dimension(ndx), asynchronous :: sdatax
    real(rk8), dimension(ndx), asynchronous :: rdatax
    real(rk8), dimension(ndy), asynchronous :: sdatay
    real(rk8), dimension(ndy), asynchronous :: rdatay
    type(mpi_request) :: mrequestx, mrequesty
    type(mpi_status) :: mstatus
    integer :: ib1 , ib2 , iex , k

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        if ( rcmpi%west /= mpi_proc_null ) &
            sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2+tb,k)
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_ineighbor_alltoallv(sdatax, counts, displs, mpi_real8, &
        rdatax, counts, displs, mpi_real8, rcmpi%cartcom, mrequestx)

    call mpi_wait(mrequestx,mstatus)

    ib2 = sizex
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        if ( rcmpi%east /= mpi_proc_null ) &
            ml(j2+iex,i1:i2+tb,k) = rdatax(ib1:ib2)
      end do
    end do

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        if ( rcmpi%south /= mpi_proc_null ) &
            sdatay(ib1:ib2) = ml(j1:,i1+(iex-1),k)
      end do
    end do
    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_ineighbor_alltoallv(sdatay, counts, displs, mpi_real8, &
        rdatay, counts, displs, mpi_real8, rcmpi%cartcom, mrequesty)

    call mpi_wait(mrequesty,mstatus)

    ib2 = sizey
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        if ( rcmpi%north /= mpi_proc_null ) &
           ml(j1:,i2+iex,k) = rdatay(ib1:ib2)
      end do
    end do

    end block transmit

  end subroutine real8_3d_exchange_rt

  subroutine real8_4d_exchange_rt(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk8) , pointer , dimension(:,:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2 , n1 , n2
    integer :: nx , ny , nk , nn
    integer :: rb , tb
    integer :: ndx , ndy , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    nn = n2-n1+1
    rb = nex
    tb = 0
    if ( rcmpi%east == mpi_proc_null ) rb = 1
    if ( rcmpi%north == mpi_proc_null ) tb = 1
    tx = ny+tb
    ty = nx+rb
    sizex = nex*tx*nk*nn
    sizey = nex*ty*nk*nn
    ndx = 2*sizex
    ndy = 2*sizey

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk8), dimension(ndx), asynchronous :: sdatax
    real(rk8), dimension(ndx), asynchronous :: rdatax
    real(rk8), dimension(ndy), asynchronous :: sdatay
    real(rk8), dimension(ndy), asynchronous :: rdatay
    type(mpi_request) :: mrequestx, mrequesty
    type(mpi_status) :: mstatus
    integer :: ib1 , ib2 , iex , k , n

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          if ( rcmpi%west /= mpi_proc_null ) &
              sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2+tb,k,n)
        end do
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_ineighbor_alltoallv(sdatax, counts, displs, mpi_real8, &
        rdatax, counts, displs, mpi_real8, rcmpi%cartcom, mrequestx)

    call mpi_wait(mrequestx,mstatus)

    ib2 = sizex
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          if ( rcmpi%east /= mpi_proc_null ) &
              ml(j2+iex,i1:i2+tb,k,n) = rdatax(ib1:ib2)
        end do
      end do
    end do

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          if ( rcmpi%south /= mpi_proc_null ) &
              sdatay(ib1:ib2) = ml(j1:,i1+(iex-1),k,n)
        end do
      end do
    end do
    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_ineighbor_alltoallv(sdatay, counts, displs, mpi_real8, &
        rdatay, counts, displs, mpi_real8, rcmpi%cartcom, mrequesty)

    call mpi_wait(mrequesty,mstatus)

    ib2 = sizey
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          if ( rcmpi%north /= mpi_proc_null ) &
             ml(j1:,i2+iex,k,n) = rdatay(ib1:ib2)
        end do
      end do
    end do

    end block transmit

  end subroutine real8_4d_exchange_rt

  subroutine real4_2d_exchange_rt(ml,nex,j1,j2,i1,i2)
    implicit none
    real(rk4) , pointer , dimension(:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2
    integer :: nx , ny
    integer :: rb , tb
    integer :: ndx , ndy , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    rb = nex
    tb = 0
    if ( rcmpi%east == mpi_proc_null ) rb = 1
    if ( rcmpi%north == mpi_proc_null ) tb = 1
    tx = ny+tb
    ty = nx+rb
    sizex = nex*tx
    sizey = nex*ty
    ndx = 2*sizex
    ndy = 2*sizey

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk4), dimension(ndx), asynchronous :: sdatax
    real(rk4), dimension(ndx), asynchronous :: rdatax
    real(rk4), dimension(ndy), asynchronous :: sdatay
    real(rk4), dimension(ndy), asynchronous :: rdatay
    type(mpi_request) :: mrequestx, mrequesty
    type(mpi_status) :: mstatus
    integer :: ib1 , ib2 , iex

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      if ( rcmpi%west /= mpi_proc_null ) &
          sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2+tb)
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_ineighbor_alltoallv(sdatax, counts, displs, mpi_real4, &
        rdatax, counts, displs, mpi_real4, rcmpi%cartcom, mrequestx)

    call mpi_wait(mrequestx,mstatus)

    ib2 = sizex
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + tx - 1
      if ( rcmpi%east /= mpi_proc_null ) &
          ml(j2+iex,i1:i2+tb) = rdatax(ib1:ib2)
    end do

    ib2 = 0
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      if ( rcmpi%south /= mpi_proc_null ) &
          sdatay(ib1:ib2) = ml(j1:,i1+(iex-1))
    end do
    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_ineighbor_alltoallv(sdatay, counts, displs, mpi_real4, &
        rdatay, counts, displs, mpi_real4, rcmpi%cartcom, mrequesty)

    call mpi_wait(mrequesty,mstatus)

    ib2 = sizey
    do iex = 1 , nex
      ib1 = ib2 + 1
      ib2 = ib1 + ty - 1
      if ( rcmpi%north /= mpi_proc_null ) &
         ml(j1:,i2+iex) = rdatay(ib1:ib2)
    end do

    end block transmit

  end subroutine real4_2d_exchange_rt

  subroutine real4_3d_exchange_rt(ml,nex,j1,j2,i1,i2,k1,k2)
    implicit none
    real(rk4) , pointer , dimension(:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2
    integer :: nx , ny , nk
    integer :: rb , tb
    integer :: ndx , ndy , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    rb = nex
    tb = 0
    if ( rcmpi%east == mpi_proc_null ) rb = 1
    if ( rcmpi%north == mpi_proc_null ) tb = 1
    tx = ny+tb
    ty = nx+rb
    sizex = nex*tx*nk
    sizey = nex*ty*nk
    ndx = 2*sizex
    ndy = 2*sizey

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk4), dimension(ndx), asynchronous :: sdatax
    real(rk4), dimension(ndx), asynchronous :: rdatax
    real(rk4), dimension(ndy), asynchronous :: sdatay
    real(rk4), dimension(ndy), asynchronous :: rdatay
    type(mpi_request) :: mrequestx, mrequesty
    type(mpi_status) :: mstatus
    integer :: ib1 , ib2 , iex , k

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        if ( rcmpi%west /= mpi_proc_null ) &
            sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2+tb,k)
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_ineighbor_alltoallv(sdatax, counts, displs, mpi_real4, &
        rdatax, counts, displs, mpi_real4, rcmpi%cartcom, mrequestx)

    call mpi_wait(mrequestx,mstatus)

    ib2 = sizex
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + tx - 1
        if ( rcmpi%east /= mpi_proc_null ) &
            ml(j2+iex,i1:i2+tb,k) = rdatax(ib1:ib2)
      end do
    end do

    ib2 = 0
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        if ( rcmpi%south /= mpi_proc_null ) &
            sdatay(ib1:ib2) = ml(j1:,i1+(iex-1),k)
      end do
    end do
    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_ineighbor_alltoallv(sdatay, counts, displs, mpi_real4, &
        rdatay, counts, displs, mpi_real4, rcmpi%cartcom, mrequesty)

    call mpi_wait(mrequesty,mstatus)

    ib2 = sizey
    do k = k1 , k2
      do iex = 1 , nex
        ib1 = ib2 + 1
        ib2 = ib1 + ty - 1
        if ( rcmpi%north /= mpi_proc_null ) &
           ml(j1:,i2+iex,k) = rdatay(ib1:ib2)
      end do
    end do

    end block transmit

  end subroutine real4_3d_exchange_rt

  subroutine real4_4d_exchange_rt(ml,nex,j1,j2,i1,i2,k1,k2,n1,n2)
    implicit none
    real(rk4) , pointer , dimension(:,:,:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: nex , j1 , j2  , i1 , i2 , k1 , k2 , n1 , n2
    integer :: nx , ny , nk , nn
    integer :: rb , tb
    integer :: ndx , ndy , tx , ty , sizex , sizey

    nx = j2-j1+1
    ny = i2-i1+1
    nk = k2-k1+1
    nn = n2-n1+1
    rb = nex
    tb = 0
    if ( rcmpi%east == mpi_proc_null ) rb = 1
    if ( rcmpi%north == mpi_proc_null ) tb = 1
    tx = ny+tb
    ty = nx+rb
    sizex = nex*tx*nk*nn
    sizey = nex*ty*nk*nn
    ndx = 2*sizex
    ndy = 2*sizey

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk4), dimension(ndx), asynchronous :: sdatax
    real(rk4), dimension(ndx), asynchronous :: rdatax
    real(rk4), dimension(ndy), asynchronous :: sdatay
    real(rk4), dimension(ndy), asynchronous :: rdatay
    type(mpi_request) :: mrequestx, mrequesty
    type(mpi_status) :: mstatus
    integer :: ib1 , ib2 , iex , k , n

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          if ( rcmpi%west /= mpi_proc_null ) &
              sdatax(ib1:ib2) = ml(j1+(iex-1),i1:i2+tb,k,n)
        end do
      end do
    end do

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_ineighbor_alltoallv(sdatax, counts, displs, mpi_real4, &
        rdatax, counts, displs, mpi_real4, rcmpi%cartcom, mrequestx)

    call mpi_wait(mrequestx,mstatus)

    ib2 = sizex
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + tx - 1
          if ( rcmpi%east /= mpi_proc_null ) &
              ml(j2+iex,i1:i2+tb,k,n) = rdatax(ib1:ib2)
        end do
      end do
    end do

    ib2 = 0
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          if ( rcmpi%south /= mpi_proc_null ) &
              sdatay(ib1:ib2) = ml(j1:,i1+(iex-1),k,n)
        end do
      end do
    end do
    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_ineighbor_alltoallv(sdatay, counts, displs, mpi_real4, &
        rdatay, counts, displs, mpi_real4, rcmpi%cartcom, mrequesty)

    call mpi_wait(mrequesty,mstatus)

    ib2 = sizey
    do n = n1 , n2
      do k = k1 , k2
        do iex = 1 , nex
          ib1 = ib2 + 1
          ib2 = ib1 + ty - 1
          if ( rcmpi%north /= mpi_proc_null ) &
             ml(j1:,i2+iex,k,n) = rdatay(ib1:ib2)
        end do
      end do
    end do

    end block transmit

  end subroutine real4_4d_exchange_rt

  subroutine real8_exchange_bdy_lr(ml,j1,j2,k1,k2)
    implicit none
    real(rk8) , pointer , dimension(:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: j1 , j2 , k1 , k2
    integer :: ndx , nk , tx , sizex

    nk = k2-k1+1
    tx = nk
    sizex = tx
    ndx = 2*sizex

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk8), dimension(ndx), asynchronous :: sdata
    real(rk8), dimension(ndx), asynchronous :: rdata
    type(mpi_request) :: mrequest
    type(mpi_status) :: mstatus
    integer :: ib1 , ib2

    ib1 = 1
    ib2 = ib1 + tx - 1
    if ( rcmpi%west /= mpi_proc_null ) sdata(ib1:ib2) = ml(j1,:)
    ib1 = ib2 + 1
    ib2 = ib1 + tx - 1
    if ( rcmpi%east /= mpi_proc_null ) sdata(ib1:ib2) = ml(j2,:)

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_ineighbor_alltoallv(sdata, counts, displs, mpi_real8, &
        rdata, counts, displs, mpi_real8, rcmpi%cartcom, mrequest)

    call mpi_wait(mrequest,mstatus)

    ib1 = 1
    ib2 = ib1 + tx - 1
    if ( rcmpi%west /= mpi_proc_null ) ml(j1-1,:) = rdata(ib1:ib2)
    ib1 = ib2 + 1
    ib2 = ib1 + tx - 1
    if ( rcmpi%east /= mpi_proc_null ) ml(j2+1,:) = rdata(ib1:ib2)

    end block transmit

  end subroutine real8_exchange_bdy_lr

  subroutine real4_exchange_bdy_lr(ml,j1,j2,k1,k2)
    implicit none
    real(rk4) , pointer , dimension(:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: j1 , j2 , k1 , k2
    integer :: ndx , nk , tx , sizex

    nk = k2-k1+1
    tx = nk
    sizex = tx
    ndx = 2*sizex

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk4), dimension(ndx), asynchronous :: sdata
    real(rk4), dimension(ndx), asynchronous :: rdata
    type(mpi_request) :: mrequest
    type(mpi_status) :: mstatus
    integer :: ib1 , ib2

    ib1 = 1
    ib2 = ib1 + tx - 1
    if ( rcmpi%west /= mpi_proc_null ) sdata(ib1:ib2) = ml(j1,:)
    ib1 = ib2 + 1
    ib2 = ib1 + tx - 1
    if ( rcmpi%east /= mpi_proc_null ) sdata(ib1:ib2) = ml(j2,:)

    counts = [ sizex, sizex, 0, 0 ]
    displs = [ 0, sizex, 2*sizex, 2*sizex ]
    call mpi_ineighbor_alltoallv(sdata, counts, displs, mpi_real4, &
        rdata, counts, displs, mpi_real4, rcmpi%cartcom, mrequest)

    call mpi_wait(mrequest,mstatus)

    ib1 = 1
    ib2 = ib1 + tx - 1
    if ( rcmpi%west /= mpi_proc_null ) ml(j1-1,:) = rdata(ib1:ib2)
    ib1 = ib2 + 1
    ib2 = ib1 + tx - 1
    if ( rcmpi%east /= mpi_proc_null ) ml(j2+1,:) = rdata(ib1:ib2)

    end block transmit

  end subroutine real4_exchange_bdy_lr

  subroutine real8_exchange_bdy_bt(ml,i1,i2,k1,k2)
    implicit none
    real(rk8) , pointer , dimension(:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: i1 , i2 , k1 , k2
    integer :: ndy , nk , ty , sizey

    nk = k2-k1+1
    ty = nk
    sizey = ty
    ndy = 2*sizey

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk8), dimension(ndy), asynchronous :: sdata
    real(rk8), dimension(ndy), asynchronous :: rdata
    type(mpi_request) :: mrequest
    type(mpi_status) :: mstatus
    integer :: ib1 , ib2

    ib2 = 0
    ib1 = ib2 + 1
    ib2 = ib1 + ty - 1
    if ( rcmpi%south /= mpi_proc_null ) sdata(ib1:ib2) = ml(i1,:)
    ib1 = ib2 + 1
    ib2 = ib1 + ty - 1
    if ( rcmpi%north /= mpi_proc_null ) sdata(ib1:ib2) = ml(i2,:)

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_ineighbor_alltoallv(sdata, counts, displs, mpi_real8, &
        rdata, counts, displs, mpi_real8, rcmpi%cartcom, mrequest)

    call mpi_wait(mrequest,mstatus)

    ib2 = 0
    ib1 = ib2 + 1
    ib2 = ib1 + ty - 1
    if ( rcmpi%south /= mpi_proc_null ) ml(i1-1,:) = rdata(ib1:ib2)
    ib1 = ib2 + 1
    ib2 = ib1 + ty - 1
    if ( rcmpi%north /= mpi_proc_null ) ml(i2+1,:) = rdata(ib1:ib2)

    end block transmit

  end subroutine real8_exchange_bdy_bt

  subroutine real4_exchange_bdy_bt(ml,i1,i2,k1,k2)
    implicit none
    real(rk4) , pointer , dimension(:,:) , intent(inout) :: ml
    integer(ik4) , intent(in) :: i1 , i2 , k1 , k2
    integer :: ndy , nk , ty , sizey

    nk = k2-k1+1
    ty = nk
    sizey = ty
    ndy = 2*sizey

    transmit : block

    integer, dimension(4), asynchronous :: counts, displs
    real(rk4), dimension(ndy), asynchronous :: sdata
    real(rk4), dimension(ndy), asynchronous :: rdata
    type(mpi_request) :: mrequest
    type(mpi_status) :: mstatus
    integer :: ib1 , ib2

    ib2 = 0
    ib1 = ib2 + 1
    ib2 = ib1 + ty - 1
    if ( rcmpi%south /= mpi_proc_null ) sdata(ib1:ib2) = ml(i1,:)
    ib1 = ib2 + 1
    ib2 = ib1 + ty - 1
    if ( rcmpi%north /= mpi_proc_null ) sdata(ib1:ib2) = ml(i2,:)

    counts = [ 0, 0, sizey, sizey ]
    displs = [ 0, 0, 0, sizey ]
    call mpi_ineighbor_alltoallv(sdata, counts, displs, mpi_real4, &
        rdata, counts, displs, mpi_real4, rcmpi%cartcom, mrequest)

    call mpi_wait(mrequest,mstatus)

    ib2 = 0
    ib1 = ib2 + 1
    ib2 = ib1 + ty - 1
    if ( rcmpi%south /= mpi_proc_null ) ml(i1-1,:) = rdata(ib1:ib2)
    ib1 = ib2 + 1
    ib2 = ib1 + ty - 1
    if ( rcmpi%north /= mpi_proc_null ) ml(i2+1,:) = rdata(ib1:ib2)

    end block transmit

  end subroutine real4_exchange_bdy_bt

end module mod_mpi_regcm

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

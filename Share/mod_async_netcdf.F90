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

module mod_async_netcdf

  use, intrinsic :: iso_c_binding
  use mod_intkinds
  use mod_realkinds
  use netcdf
#ifdef OPENACC
  use cudafor
#endif

  implicit none

  private

  public :: async_netcdf_put_var
  public :: async_netcdf_put_var_device
  public :: async_netcdf_put_var_rkx_as_r4_device
  public :: async_netcdf_get_var_rkx
  public :: async_netcdf_configured_cap_bytes
  public :: async_netcdf_initialize
  public :: async_netcdf_is_enabled
  public :: async_netcdf_sync
  public :: async_netcdf_close
  public :: async_netcdf_wait_all
  public :: async_netcdf_shutdown
  public :: async_netcdf_lock
  public :: async_netcdf_unlock
  public :: async_netcdf_defer_wake_begin
  public :: async_netcdf_defer_wake_end
  public :: async_netcdf_buffer
  public :: async_netcdf_counter
  public :: async_netcdf_try_acquire_buffer
  public :: async_netcdf_release_buffer
  public :: async_netcdf_counter_reset
  public :: async_netcdf_counter_poll
  public :: async_netcdf_counter_wait
  public :: async_netcdf_counter_fail

  integer(ik4), parameter :: async_ok = nf90_noerr
  integer(ik4), parameter :: async_err = -9999

  integer(ik4), parameter :: async_write_r4 = 1
  integer(ik4), parameter :: async_write_r8 = 2
  integer(ik4), parameter :: async_write_i4 = 3
  integer(ik4), parameter :: async_sync = 4
  integer(ik4), parameter :: async_read_r4 = 5
  integer(ik4), parameter :: async_read_r8 = 6
  integer(ik4), parameter :: async_close = 7

  integer(c_int64_t), parameter :: bytes_r4 = 4_c_int64_t
  integer(c_int64_t), parameter :: bytes_r8 = 8_c_int64_t
  integer(c_int64_t), parameter :: bytes_i4 = 4_c_int64_t

  integer(c_int64_t), parameter :: default_cap_bytes = 0_c_int64_t
  integer(c_int64_t), parameter :: pool_alignment = 128_c_int64_t

  type async_netcdf_counter
    integer(ik4) :: pending = 0
    integer(ik4) :: status = async_ok
  end type async_netcdf_counter

  type async_netcdf_buffer
    real(rkx), pointer, contiguous, dimension(:) :: rkx => null()
    integer(c_int64_t) :: bytes = 0_c_int64_t
#ifdef OPENACC
    type(c_ptr) :: cptr = c_null_ptr
    integer(c_int64_t) :: pool_offset = 0_c_int64_t
    logical :: pooled_buffer = .false.
#endif
    logical :: active = .false.
  end type async_netcdf_buffer

  type async_item
    integer(ik4) :: op = 0
    integer(ik4) :: ncid = -1
    integer(ik4) :: varid = -1
    integer(ik4) :: nd = 0
    integer(ik4) :: nvals = 0
    integer(ik4), dimension(5) :: start = 1
    integer(ik4), dimension(5) :: count = 1
    logical :: has_start = .false.
    logical :: has_count = .false.
    logical :: scalar = .false.
    logical :: external_buffer = .false.
    logical :: speculative = .false.
    integer(c_int64_t) :: bytes = 0_c_int64_t
    type(async_netcdf_counter), pointer :: counter => null()
#ifdef OPENACC
    type(c_ptr) :: cptr = c_null_ptr
    integer(c_int64_t) :: buffer_bytes = 0_c_int64_t
    integer(c_int64_t) :: pool_offset = 0_c_int64_t
    logical :: pooled_buffer = .false.
#endif
    real(rk4), pointer, contiguous, dimension(:) :: r4 => null()
    real(rk8), pointer, contiguous, dimension(:) :: r8 => null()
    integer(ik4), pointer, contiguous, dimension(:) :: i4 => null()
    type(async_item), pointer :: next => null()
  end type async_item

  type(async_item), pointer, save :: queue_head => null()
  type(async_item), pointer, save :: queue_tail => null()
  type(async_item), pointer, save :: read_queue_head => null()
  type(async_item), pointer, save :: read_queue_tail => null()
#ifdef OPENACC
  type pool_free_block
    integer(c_int64_t) :: offset = 0_c_int64_t
    integer(c_int64_t) :: bytes = 0_c_int64_t
    type(pool_free_block), pointer :: next => null()
  end type pool_free_block

  type(c_ptr), save :: pool_cptr = c_null_ptr
  integer(c_int64_t), save :: pool_capacity = 0_c_int64_t
  integer(c_int64_t), save :: pool_used = 0_c_int64_t
  type(pool_free_block), pointer, save :: pool_free_head => null()
#endif

  logical, save :: initialized = .false.
  logical, save :: enabled = .true.
  logical, save :: stop_requested = .false.
  logical, save :: worker_running = .false.
  integer(ik4), save :: first_error = async_ok
  integer(ik4), save :: active_items = 0
  integer(ik4), save :: defer_work_signals = 0
  integer(c_int64_t), save :: max_bytes = default_cap_bytes
  integer(c_int64_t), save :: env_max_bytes = default_cap_bytes
  integer(c_int64_t), save :: reserved_bytes = 0_c_int64_t
  logical, save :: memory_cap_read = .false.
  logical, save :: deferred_work_pending = .false.

  interface async_netcdf_put_var
    module procedure async_put_var_r4_1d
    module procedure async_put_var_r8_1d
    module procedure async_put_var_i4_1d
    module procedure async_put_var_r4_scalar
    module procedure async_put_var_r8_scalar
    module procedure async_put_var_i4_scalar
  end interface async_netcdf_put_var

  interface async_netcdf_put_var_device
    module procedure async_put_var_r4_device_1d
    module procedure async_put_var_r8_device_1d
    module procedure async_put_var_i4_device_1d
  end interface async_netcdf_put_var_device

#ifdef ASYNC_NETCDF
  interface
    integer(c_int) function regcm_async_thread_start() bind(C, &
        name='regcm_async_thread_start')
      import :: c_int
    end function regcm_async_thread_start

    integer(c_int) function regcm_async_thread_join() bind(C, &
        name='regcm_async_thread_join')
      import :: c_int
    end function regcm_async_thread_join

    subroutine regcm_async_queue_lock() bind(C, name='regcm_async_queue_lock')
    end subroutine regcm_async_queue_lock

    subroutine regcm_async_queue_unlock() bind(C, &
        name='regcm_async_queue_unlock')
    end subroutine regcm_async_queue_unlock

    subroutine regcm_async_wait_work() bind(C, name='regcm_async_wait_work')
    end subroutine regcm_async_wait_work

    subroutine regcm_async_signal_work() bind(C, name='regcm_async_signal_work')
    end subroutine regcm_async_signal_work

    subroutine regcm_async_broadcast_work() bind(C, &
        name='regcm_async_broadcast_work')
    end subroutine regcm_async_broadcast_work

    subroutine regcm_async_wait_space() bind(C, name='regcm_async_wait_space')
    end subroutine regcm_async_wait_space

    subroutine regcm_async_broadcast_space() bind(C, &
        name='regcm_async_broadcast_space')
    end subroutine regcm_async_broadcast_space

    subroutine regcm_async_wait_idle() bind(C, name='regcm_async_wait_idle')
    end subroutine regcm_async_wait_idle

    subroutine regcm_async_broadcast_idle() bind(C, &
        name='regcm_async_broadcast_idle')
    end subroutine regcm_async_broadcast_idle

    subroutine regcm_async_netcdf_lock() bind(C, name='regcm_async_netcdf_lock')
    end subroutine regcm_async_netcdf_lock

    subroutine regcm_async_netcdf_unlock() bind(C, &
        name='regcm_async_netcdf_unlock')
    end subroutine regcm_async_netcdf_unlock

    subroutine regcm_async_worker_netcdf_lock() bind(C, &
        name='regcm_async_worker_netcdf_lock')
    end subroutine regcm_async_worker_netcdf_lock

    subroutine regcm_async_worker_netcdf_unlock() bind(C, &
        name='regcm_async_worker_netcdf_unlock')
    end subroutine regcm_async_worker_netcdf_unlock

#ifdef OPENACC
    type(c_ptr) function regcm_async_ptr_offset(base,offset) bind(C, &
        name='regcm_async_ptr_offset')
      import :: c_ptr, c_int64_t
      type(c_ptr), value :: base
      integer(c_int64_t), value :: offset
    end function regcm_async_ptr_offset
#endif
  end interface
#endif

  contains

#ifndef ASYNC_NETCDF
    integer(c_int) function regcm_async_thread_start()
      implicit none
      regcm_async_thread_start = 1_c_int
    end function regcm_async_thread_start

    integer(c_int) function regcm_async_thread_join()
      implicit none
      regcm_async_thread_join = 0_c_int
    end function regcm_async_thread_join

    subroutine regcm_async_queue_lock()
      implicit none
    end subroutine regcm_async_queue_lock

    subroutine regcm_async_queue_unlock()
      implicit none
    end subroutine regcm_async_queue_unlock

    subroutine regcm_async_wait_work()
      implicit none
    end subroutine regcm_async_wait_work

    subroutine regcm_async_signal_work()
      implicit none
    end subroutine regcm_async_signal_work

    subroutine regcm_async_broadcast_work()
      implicit none
    end subroutine regcm_async_broadcast_work

    subroutine regcm_async_wait_space()
      implicit none
    end subroutine regcm_async_wait_space

    subroutine regcm_async_broadcast_space()
      implicit none
    end subroutine regcm_async_broadcast_space

    subroutine regcm_async_wait_idle()
      implicit none
    end subroutine regcm_async_wait_idle

    subroutine regcm_async_broadcast_idle()
      implicit none
    end subroutine regcm_async_broadcast_idle

    subroutine regcm_async_netcdf_lock()
      implicit none
    end subroutine regcm_async_netcdf_lock

    subroutine regcm_async_netcdf_unlock()
      implicit none
    end subroutine regcm_async_netcdf_unlock

    subroutine regcm_async_worker_netcdf_lock()
      implicit none
    end subroutine regcm_async_worker_netcdf_lock

    subroutine regcm_async_worker_netcdf_unlock()
      implicit none
    end subroutine regcm_async_worker_netcdf_unlock

#ifdef OPENACC
    type(c_ptr) function regcm_async_ptr_offset(base,offset)
      implicit none
      type(c_ptr), intent(in) :: base
      integer(c_int64_t), intent(in) :: offset

      regcm_async_ptr_offset = base
    end function regcm_async_ptr_offset
#endif
#endif

    subroutine async_netcdf_worker() bind(C, name='regcm_async_netcdf_worker')
      implicit none
      type(async_item), pointer :: item
      integer(ik4) :: ncstat_local

      do
        call regcm_async_queue_lock()
        do while ( ( .not. associated(read_queue_head) .and. &
                   ( .not. associated(queue_head) .or. &
                     deferred_work_pending ) ) .and. .not. stop_requested )
          if ( .not. associated(read_queue_head) .and. &
               .not. associated(queue_head) ) call regcm_async_broadcast_idle()
          call regcm_async_wait_work()
        end do
        if ( stop_requested .and. .not. associated(read_queue_head) .and. &
             .not. associated(queue_head) ) then
          call regcm_async_broadcast_idle()
          call regcm_async_queue_unlock()
          exit
        end if
        if ( associated(read_queue_head) ) then
          item => read_queue_head
          read_queue_head => item%next
          if ( .not. associated(read_queue_head) ) read_queue_tail => null()
        else
          item => queue_head
          queue_head => item%next
          if ( .not. associated(queue_head) ) queue_tail => null()
        end if
        item%next => null()
        active_items = active_items + 1
        call regcm_async_queue_unlock()

        call regcm_async_worker_netcdf_lock()
        ncstat_local = execute_item(item)
        call regcm_async_worker_netcdf_unlock()

        call regcm_async_queue_lock()
        if ( associated(item%counter) ) then
          call complete_counter_locked(item%counter,ncstat_local)
        end if
        if ( .not. item%speculative .and. first_error == async_ok .and. &
             ncstat_local /= async_ok ) then
          first_error = ncstat_local
        end if
        active_items = active_items - 1
        reserved_bytes = reserved_bytes - item%bytes
        call regcm_async_broadcast_space()
        if ( .not. associated(read_queue_head) .and. &
             .not. associated(queue_head) .and. active_items == 0 ) then
          call regcm_async_broadcast_idle()
        end if
        call regcm_async_queue_unlock()

        call free_item(item)
      end do
    end subroutine async_netcdf_worker

    integer(ik4) function async_put_var_r4_1d(ncid,varid,values,start,count)
      implicit none
      integer(ik4), intent(in) :: ncid, varid
      real(rk4), contiguous, dimension(:), intent(in) :: values
      integer(ik4), dimension(:), intent(in), optional :: start, count
      async_put_var_r4_1d = enqueue_r4(ncid,varid,values,.false.,start,count)
    end function async_put_var_r4_1d

    integer(ik4) function async_put_var_r8_1d(ncid,varid,values,start,count)
      implicit none
      integer(ik4), intent(in) :: ncid, varid
      real(rk8), contiguous, dimension(:), intent(in) :: values
      integer(ik4), dimension(:), intent(in), optional :: start, count
      async_put_var_r8_1d = enqueue_r8(ncid,varid,values,.false.,start,count)
    end function async_put_var_r8_1d

    integer(ik4) function async_put_var_i4_1d(ncid,varid,values,start,count)
      implicit none
      integer(ik4), intent(in) :: ncid, varid
      integer(ik4), contiguous, dimension(:), intent(in) :: values
      integer(ik4), dimension(:), intent(in), optional :: start, count
      async_put_var_i4_1d = enqueue_i4(ncid,varid,values,.false.,start,count)
    end function async_put_var_i4_1d

    integer(ik4) function async_put_var_r4_scalar(ncid,varid,value,start,count)
      implicit none
      integer(ik4), intent(in) :: ncid, varid
      real(rk4), intent(in) :: value
      integer(ik4), dimension(:), intent(in), optional :: start, count
      real(rk4), dimension(1) :: values
      values(1) = value
      async_put_var_r4_scalar = enqueue_r4(ncid,varid,values,.true.,start,count)
    end function async_put_var_r4_scalar

    integer(ik4) function async_put_var_r8_scalar(ncid,varid,value,start,count)
      implicit none
      integer(ik4), intent(in) :: ncid, varid
      real(rk8), intent(in) :: value
      integer(ik4), dimension(:), intent(in), optional :: start, count
      real(rk8), dimension(1) :: values
      values(1) = value
      async_put_var_r8_scalar = enqueue_r8(ncid,varid,values,.true.,start,count)
    end function async_put_var_r8_scalar

    integer(ik4) function async_put_var_i4_scalar(ncid,varid,value,start,count)
      implicit none
      integer(ik4), intent(in) :: ncid, varid
      integer(ik4), intent(in) :: value
      integer(ik4), dimension(:), intent(in), optional :: start, count
      integer(ik4), dimension(1) :: values
      values(1) = value
      async_put_var_i4_scalar = enqueue_i4(ncid,varid,values,.true.,start,count)
    end function async_put_var_i4_scalar

    integer(ik4) function async_put_var_r4_device_1d(ncid,varid,values,start,count)
      implicit none
      integer(ik4), intent(in) :: ncid, varid
      real(rk4), contiguous, dimension(:), intent(in) :: values
      integer(ik4), dimension(:), intent(in), optional :: start, count
      async_put_var_r4_device_1d = &
        enqueue_r4_device(ncid,varid,values,start,count)
    end function async_put_var_r4_device_1d

    integer(ik4) function async_put_var_r8_device_1d(ncid,varid,values,start,count)
      implicit none
      integer(ik4), intent(in) :: ncid, varid
      real(rk8), contiguous, dimension(:), intent(in) :: values
      integer(ik4), dimension(:), intent(in), optional :: start, count
      async_put_var_r8_device_1d = &
        enqueue_r8_device(ncid,varid,values,start,count)
    end function async_put_var_r8_device_1d

    integer(ik4) function async_put_var_i4_device_1d(ncid,varid,values,start,count)
      implicit none
      integer(ik4), intent(in) :: ncid, varid
      integer(ik4), contiguous, dimension(:), intent(in) :: values
      integer(ik4), dimension(:), intent(in), optional :: start, count
      async_put_var_i4_device_1d = &
        enqueue_i4_device(ncid,varid,values,start,count)
    end function async_put_var_i4_device_1d

    integer(ik4) function async_netcdf_put_var_rkx_as_r4_device(ncid,varid, &
        values,start,count)
      implicit none
      integer(ik4), intent(in) :: ncid, varid
      real(rkx), contiguous, dimension(:), intent(in), target :: values
      integer(ik4), dimension(:), intent(in), optional :: start, count
      async_netcdf_put_var_rkx_as_r4_device = &
        enqueue_rkx_as_r4_device(ncid,varid,values,start,count)
    end function async_netcdf_put_var_rkx_as_r4_device

    integer(c_int64_t) function async_netcdf_configured_cap_bytes()
      implicit none

      if ( .not. memory_cap_read ) call read_memory_cap()
      async_netcdf_configured_cap_bytes = env_max_bytes
    end function async_netcdf_configured_cap_bytes

    logical function async_netcdf_initialize(limit_bytes)
      implicit none
      integer(c_int64_t), intent(in), optional :: limit_bytes

      if ( present(limit_bytes) ) then
        call async_init(limit_bytes)
      else
        call async_init()
      end if
      async_netcdf_initialize = enabled .and. first_error == async_ok
    end function async_netcdf_initialize

    logical function async_netcdf_is_enabled()
      implicit none

      call async_init()
      async_netcdf_is_enabled = enabled .and. first_error == async_ok
    end function async_netcdf_is_enabled

    integer(ik4) function async_netcdf_get_var_rkx(ncid,varid,values,start, &
        count,counter)
      implicit none
      integer(ik4), intent(in) :: ncid, varid
      real(rkx), contiguous, dimension(:), intent(inout), target :: values
      integer(ik4), dimension(:), intent(in) :: start, count
      type(async_netcdf_counter), target, intent(inout) :: counter

#ifdef SINGLE_PRECISION_REAL
      async_netcdf_get_var_rkx = enqueue_read_r4(ncid,varid,values,start, &
        count,counter)
#else
      async_netcdf_get_var_rkx = enqueue_read_r8(ncid,varid,values,start, &
        count,counter)
#endif
    end function async_netcdf_get_var_rkx

    subroutine async_netcdf_try_acquire_buffer(buffer,bytes,stat)
      implicit none
      type(async_netcdf_buffer), intent(inout) :: buffer
      integer(c_int64_t), intent(in) :: bytes
      integer(ik4), intent(out) :: stat
      integer(c_int64_t) :: aligned_bytes, nvals
#ifdef OPENACC
      logical :: have_buffer
#else
      integer(ik4) :: istat
#endif

      call async_init()
      stat = async_err
      if ( buffer%active ) return
      if ( first_error /= async_ok .or. .not. enabled ) return
      if ( bytes <= 0_c_int64_t .or. max_bytes <= 0_c_int64_t ) return

      aligned_bytes = align_pool_bytes(bytes)
      nvals = aligned_bytes / bytes_rkx()

      call regcm_async_queue_lock()
      if ( reserved_bytes + aligned_bytes <= max_bytes ) then
#ifdef OPENACC
        if ( c_associated(pool_cptr) .and. aligned_bytes <= pool_capacity ) then
          call take_pool_slice_locked(aligned_bytes,buffer%pool_offset, &
            have_buffer)
          if ( have_buffer ) then
            reserved_bytes = reserved_bytes + aligned_bytes
            buffer%bytes = aligned_bytes
            buffer%cptr = regcm_async_ptr_offset(pool_cptr,buffer%pool_offset)
            buffer%pooled_buffer = .true.
            buffer%active = .true.
            stat = async_ok
          end if
        end if
#else
        reserved_bytes = reserved_bytes + aligned_bytes
        buffer%bytes = aligned_bytes
        buffer%active = .true.
        stat = async_ok
#endif
      end if
      call regcm_async_queue_unlock()
      if ( stat /= async_ok ) return

#ifdef OPENACC
      call c_f_pointer(buffer%cptr,buffer%rkx,[int(nvals,ik4)])
#else
      allocate(buffer%rkx(int(nvals,ik4)),stat=istat)
      if ( istat /= 0 ) then
        call async_netcdf_release_buffer(buffer)
        stat = async_err
      end if
#endif
    end subroutine async_netcdf_try_acquire_buffer

    subroutine async_netcdf_release_buffer(buffer)
      implicit none
      type(async_netcdf_buffer), intent(inout) :: buffer
#ifdef OPENACC
      integer(c_int) :: istat
#endif

      if ( .not. buffer%active ) return
#ifdef OPENACC
      if ( buffer%pooled_buffer .and. c_associated(buffer%cptr) ) then
        call release_item_buffer(buffer%cptr,buffer%pool_offset,buffer%bytes)
      else if ( c_associated(buffer%cptr) ) then
        istat = cudaFreeHost(buffer%cptr)
        buffer%cptr = c_null_ptr
      end if
      if ( associated(buffer%rkx) ) buffer%rkx => null()
      buffer%pool_offset = 0_c_int64_t
      buffer%pooled_buffer = .false.
#else
      if ( associated(buffer%rkx) ) deallocate(buffer%rkx)
#endif
      call release_reservation(buffer%bytes)
      buffer%bytes = 0_c_int64_t
      buffer%active = .false.
    end subroutine async_netcdf_release_buffer

    subroutine async_netcdf_counter_reset(counter)
      implicit none
      type(async_netcdf_counter), intent(inout) :: counter

      call regcm_async_queue_lock()
      counter%pending = 0
      counter%status = async_ok
      call regcm_async_queue_unlock()
    end subroutine async_netcdf_counter_reset

    subroutine async_netcdf_counter_poll(counter,pending,status)
      implicit none
      type(async_netcdf_counter), intent(inout) :: counter
      integer(ik4), intent(out) :: pending
      integer(ik4), intent(out) :: status

      call regcm_async_queue_lock()
      pending = counter%pending
      status = counter%status
      call regcm_async_queue_unlock()
    end subroutine async_netcdf_counter_poll

    subroutine async_netcdf_counter_wait(counter,status)
      implicit none
      type(async_netcdf_counter), intent(inout) :: counter
      integer(ik4), intent(out) :: status

      call regcm_async_queue_lock()
      do while ( counter%pending > 0 )
        call regcm_async_wait_space()
      end do
      status = counter%status
      call regcm_async_queue_unlock()
    end subroutine async_netcdf_counter_wait

    subroutine async_netcdf_counter_fail(counter,status)
      implicit none
      type(async_netcdf_counter), intent(inout) :: counter
      integer(ik4), intent(in) :: status

      call regcm_async_queue_lock()
      if ( counter%status == async_ok ) counter%status = status
      call regcm_async_broadcast_space()
      call regcm_async_queue_unlock()
    end subroutine async_netcdf_counter_fail

    subroutine async_netcdf_lock()
      implicit none
      call regcm_async_netcdf_lock()
    end subroutine async_netcdf_lock

    subroutine async_netcdf_unlock()
      implicit none
      call regcm_async_netcdf_unlock()
    end subroutine async_netcdf_unlock

    subroutine async_netcdf_defer_wake_begin()
      implicit none

      call async_init()
      if ( .not. enabled ) return

      call regcm_async_queue_lock()
      defer_work_signals = defer_work_signals + 1
      deferred_work_pending = .true.
      call regcm_async_queue_unlock()
    end subroutine async_netcdf_defer_wake_begin

    subroutine async_netcdf_defer_wake_end()
      implicit none

      if ( .not. initialized .or. .not. enabled ) return

      call regcm_async_queue_lock()
      if ( defer_work_signals > 0 ) then
        defer_work_signals = defer_work_signals - 1
      end if
      if ( defer_work_signals == 0 ) then
        call release_deferred_work_locked()
      end if
      call regcm_async_queue_unlock()
    end subroutine async_netcdf_defer_wake_end

    integer(ik4) function async_netcdf_sync(ncid)
      implicit none
      integer(ik4), intent(in) :: ncid
      type(async_item), pointer :: item
      integer(ik4) :: stat

      call async_init()
      if ( first_error /= async_ok ) then
        async_netcdf_sync = first_error
        return
      end if
      if ( .not. enabled ) then
        async_netcdf_sync = nf90_sync(ncid)
        return
      end if

      allocate(item, stat=stat)
      if ( stat /= 0 ) then
        async_netcdf_sync = async_err
        return
      end if
      item%op = async_sync
      item%ncid = ncid

      call regcm_async_queue_lock()
      call append_item(item)
      call signal_work_locked()
      call regcm_async_queue_unlock()
      async_netcdf_sync = async_ok
    end function async_netcdf_sync

    integer(ik4) function async_netcdf_close(ncid)
      implicit none
      integer(ik4), intent(in) :: ncid
      type(async_item), pointer :: item
      integer(ik4) :: stat

      call async_init()
      if ( first_error /= async_ok ) then
        async_netcdf_close = first_error
        return
      end if
      if ( .not. enabled ) then
        async_netcdf_close = nf90_close(ncid)
        return
      end if

      allocate(item, stat=stat)
      if ( stat /= 0 ) then
        async_netcdf_close = async_err
        return
      end if
      item%op = async_close
      item%ncid = ncid

      call regcm_async_queue_lock()
      call append_item(item)
      call signal_work_locked()
      call regcm_async_queue_unlock()
      async_netcdf_close = async_ok
    end function async_netcdf_close

    integer(ik4) function async_netcdf_wait_all()
      implicit none

      if ( .not. initialized ) then
        async_netcdf_wait_all = first_error
        return
      end if

      call regcm_async_queue_lock()
      call release_deferred_work_locked()
      do while ( associated(read_queue_head) .or. associated(queue_head) .or. &
                 active_items > 0 )
        call regcm_async_wait_idle()
      end do
      async_netcdf_wait_all = first_error
      call regcm_async_queue_unlock()
    end function async_netcdf_wait_all

    integer(ik4) function async_netcdf_shutdown()
      implicit none
      integer(c_int) :: rc

      if ( .not. initialized ) then
        async_netcdf_shutdown = first_error
        return
      end if

      call regcm_async_queue_lock()
      call release_deferred_work_locked()
      do while ( associated(read_queue_head) .or. associated(queue_head) .or. &
                 active_items > 0 )
        call regcm_async_wait_idle()
      end do
      stop_requested = .true.
      call regcm_async_broadcast_work()
      call regcm_async_queue_unlock()

      if ( worker_running ) then
        rc = regcm_async_thread_join()
        if ( rc /= 0 .and. first_error == async_ok ) first_error = async_err
        worker_running = .false.
      end if
#ifdef OPENACC
      call free_pool_buffers()
#endif
      initialized = .false.
      stop_requested = .false.
      defer_work_signals = 0
      deferred_work_pending = .false.
      async_netcdf_shutdown = first_error
    end function async_netcdf_shutdown

    integer(ik4) function enqueue_r4(ncid,varid,values,scalar,start,count)
      implicit none
      integer(ik4), intent(in) :: ncid, varid
      real(rk4), contiguous, dimension(:), intent(in) :: values
      logical, intent(in) :: scalar
      integer(ik4), dimension(:), intent(in), optional :: start, count
      type(async_item), pointer :: item
      integer(ik4) :: nvals, stat
      integer(c_int64_t) :: bytes

      call async_init()
      nvals = size(values)
      bytes = int(nvals,c_int64_t)*bytes_r4
      if ( first_error /= async_ok ) then
        enqueue_r4 = first_error
      else if ( .not. enabled ) then
        enqueue_r4 = put_var_r4_sync(ncid,varid,values,scalar,start,count)
      else
        call reserve_memory(bytes,stat)
        if ( stat /= async_ok ) then
          enqueue_r4 = stat
          return
        end if
        allocate(item, stat=stat)
        if ( stat /= 0 ) then
          call release_reservation(bytes)
          enqueue_r4 = async_err
          return
        end if
        item%op = async_write_r4
        item%ncid = ncid
        item%varid = varid
        item%nvals = nvals
        item%bytes = bytes
        item%scalar = scalar
        call copy_bounds(item,start,count)
        call allocate_item_r4(item,nvals,bytes,stat)
        if ( stat /= async_ok ) then
          call release_reservation(bytes)
          call free_item(item)
          enqueue_r4 = stat
          return
        end if
        item%r4(1:nvals) = values(1:nvals)
        call enqueue_item(item)
        enqueue_r4 = async_ok
      end if
    end function enqueue_r4

    integer(ik4) function enqueue_r8(ncid,varid,values,scalar,start,count)
      implicit none
      integer(ik4), intent(in) :: ncid, varid
      real(rk8), contiguous, dimension(:), intent(in) :: values
      logical, intent(in) :: scalar
      integer(ik4), dimension(:), intent(in), optional :: start, count
      type(async_item), pointer :: item
      integer(ik4) :: nvals, stat
      integer(c_int64_t) :: bytes

      call async_init()
      nvals = size(values)
      bytes = int(nvals,c_int64_t)*bytes_r8
      if ( first_error /= async_ok ) then
        enqueue_r8 = first_error
      else if ( .not. enabled ) then
        enqueue_r8 = put_var_r8_sync(ncid,varid,values,scalar,start,count)
      else
        call reserve_memory(bytes,stat)
        if ( stat /= async_ok ) then
          enqueue_r8 = stat
          return
        end if
        allocate(item, stat=stat)
        if ( stat /= 0 ) then
          call release_reservation(bytes)
          enqueue_r8 = async_err
          return
        end if
        item%op = async_write_r8
        item%ncid = ncid
        item%varid = varid
        item%nvals = nvals
        item%bytes = bytes
        item%scalar = scalar
        call copy_bounds(item,start,count)
        call allocate_item_r8(item,nvals,bytes,stat)
        if ( stat /= async_ok ) then
          call release_reservation(bytes)
          call free_item(item)
          enqueue_r8 = stat
          return
        end if
        item%r8(1:nvals) = values(1:nvals)
        call enqueue_item(item)
        enqueue_r8 = async_ok
      end if
    end function enqueue_r8

    integer(ik4) function enqueue_i4(ncid,varid,values,scalar,start,count)
      implicit none
      integer(ik4), intent(in) :: ncid, varid
      integer(ik4), contiguous, dimension(:), intent(in) :: values
      logical, intent(in) :: scalar
      integer(ik4), dimension(:), intent(in), optional :: start, count
      type(async_item), pointer :: item
      integer(ik4) :: nvals, stat
      integer(c_int64_t) :: bytes

      call async_init()
      nvals = size(values)
      bytes = int(nvals,c_int64_t)*bytes_i4
      if ( first_error /= async_ok ) then
        enqueue_i4 = first_error
      else if ( .not. enabled ) then
        enqueue_i4 = put_var_i4_sync(ncid,varid,values,scalar,start,count)
      else
        call reserve_memory(bytes,stat)
        if ( stat /= async_ok ) then
          enqueue_i4 = stat
          return
        end if
        allocate(item, stat=stat)
        if ( stat /= 0 ) then
          call release_reservation(bytes)
          enqueue_i4 = async_err
          return
        end if
        item%op = async_write_i4
        item%ncid = ncid
        item%varid = varid
        item%nvals = nvals
        item%bytes = bytes
        item%scalar = scalar
        call copy_bounds(item,start,count)
        call allocate_item_i4(item,nvals,bytes,stat)
        if ( stat /= async_ok ) then
          call release_reservation(bytes)
          call free_item(item)
          enqueue_i4 = stat
          return
        end if
        item%i4(1:nvals) = values(1:nvals)
        call enqueue_item(item)
        enqueue_i4 = async_ok
      end if
    end function enqueue_i4

    integer(ik4) function enqueue_r4_device(ncid,varid,values,start,count)
      implicit none
      integer(ik4), intent(in) :: ncid, varid
      real(rk4), contiguous, dimension(:), intent(in) :: values
      integer(ik4), dimension(:), intent(in), optional :: start, count
      type(async_item), pointer :: item
      integer(ik4) :: nvals, stat
      integer(c_int64_t) :: bytes

      call async_init()
      nvals = size(values)
      bytes = int(nvals,c_int64_t)*bytes_r4
      if ( first_error /= async_ok ) then
        enqueue_r4_device = first_error
      else if ( .not. enabled ) then
        enqueue_r4_device = enqueue_r4_device_sync(ncid,varid,values,start,count)
      else
        call reserve_memory(bytes,stat)
        if ( stat /= async_ok ) then
          enqueue_r4_device = stat
          return
        end if
        allocate(item, stat=stat)
        if ( stat /= 0 ) then
          call release_reservation(bytes)
          enqueue_r4_device = async_err
          return
        end if
        item%op = async_write_r4
        item%ncid = ncid
        item%varid = varid
        item%nvals = nvals
        item%bytes = bytes
        call copy_bounds(item,start,count)
        call allocate_item_r4(item,nvals,bytes,stat)
        if ( stat /= async_ok ) then
          call release_reservation(bytes)
          call free_item(item)
          enqueue_r4_device = stat
          return
        end if
        call copy_device_to_item_r4(item,values,nvals,stat)
        if ( stat /= async_ok ) then
          call release_reservation(bytes)
          call free_item(item)
          enqueue_r4_device = stat
          return
        end if
        call enqueue_item(item)
        enqueue_r4_device = async_ok
      end if
    end function enqueue_r4_device

    integer(ik4) function enqueue_r8_device(ncid,varid,values,start,count)
      implicit none
      integer(ik4), intent(in) :: ncid, varid
      real(rk8), contiguous, dimension(:), intent(in) :: values
      integer(ik4), dimension(:), intent(in), optional :: start, count
      type(async_item), pointer :: item
      integer(ik4) :: nvals, stat
      integer(c_int64_t) :: bytes

      call async_init()
      nvals = size(values)
      bytes = int(nvals,c_int64_t)*bytes_r8
      if ( first_error /= async_ok ) then
        enqueue_r8_device = first_error
      else if ( .not. enabled ) then
        enqueue_r8_device = enqueue_r8_device_sync(ncid,varid,values,start,count)
      else
        call reserve_memory(bytes,stat)
        if ( stat /= async_ok ) then
          enqueue_r8_device = stat
          return
        end if
        allocate(item, stat=stat)
        if ( stat /= 0 ) then
          call release_reservation(bytes)
          enqueue_r8_device = async_err
          return
        end if
        item%op = async_write_r8
        item%ncid = ncid
        item%varid = varid
        item%nvals = nvals
        item%bytes = bytes
        call copy_bounds(item,start,count)
        call allocate_item_r8(item,nvals,bytes,stat)
        if ( stat /= async_ok ) then
          call release_reservation(bytes)
          call free_item(item)
          enqueue_r8_device = stat
          return
        end if
        call copy_device_to_item_r8(item,values,nvals,stat)
        if ( stat /= async_ok ) then
          call release_reservation(bytes)
          call free_item(item)
          enqueue_r8_device = stat
          return
        end if
        call enqueue_item(item)
        enqueue_r8_device = async_ok
      end if
    end function enqueue_r8_device

    integer(ik4) function enqueue_i4_device(ncid,varid,values,start,count)
      implicit none
      integer(ik4), intent(in) :: ncid, varid
      integer(ik4), contiguous, dimension(:), intent(in) :: values
      integer(ik4), dimension(:), intent(in), optional :: start, count
      type(async_item), pointer :: item
      integer(ik4) :: nvals, stat
      integer(c_int64_t) :: bytes

      call async_init()
      nvals = size(values)
      bytes = int(nvals,c_int64_t)*bytes_i4
      if ( first_error /= async_ok ) then
        enqueue_i4_device = first_error
      else if ( .not. enabled ) then
        enqueue_i4_device = enqueue_i4_device_sync(ncid,varid,values,start,count)
      else
        call reserve_memory(bytes,stat)
        if ( stat /= async_ok ) then
          enqueue_i4_device = stat
          return
        end if
        allocate(item, stat=stat)
        if ( stat /= 0 ) then
          call release_reservation(bytes)
          enqueue_i4_device = async_err
          return
        end if
        item%op = async_write_i4
        item%ncid = ncid
        item%varid = varid
        item%nvals = nvals
        item%bytes = bytes
        call copy_bounds(item,start,count)
        call allocate_item_i4(item,nvals,bytes,stat)
        if ( stat /= async_ok ) then
          call release_reservation(bytes)
          call free_item(item)
          enqueue_i4_device = stat
          return
        end if
        call copy_device_to_item_i4(item,values,nvals,stat)
        if ( stat /= async_ok ) then
          call release_reservation(bytes)
          call free_item(item)
          enqueue_i4_device = stat
          return
        end if
        call enqueue_item(item)
        enqueue_i4_device = async_ok
      end if
    end function enqueue_i4_device

    integer(ik4) function enqueue_rkx_as_r4_device(ncid,varid,values,start, &
        count)
      implicit none
      integer(ik4), intent(in) :: ncid, varid
      real(rkx), contiguous, dimension(:), intent(in), target :: values
      integer(ik4), dimension(:), intent(in), optional :: start, count
      type(async_item), pointer :: item
      integer(ik4) :: nvals, stat
      integer(c_int64_t) :: bytes

      call async_init()
      nvals = size(values)
      bytes = int(nvals,c_int64_t)*bytes_r4
      if ( first_error /= async_ok ) then
        enqueue_rkx_as_r4_device = first_error
      else if ( .not. enabled ) then
        enqueue_rkx_as_r4_device = &
          enqueue_rkx_as_r4_device_sync(ncid,varid,values,start,count)
      else
        call reserve_memory(bytes,stat)
        if ( stat /= async_ok ) then
          enqueue_rkx_as_r4_device = stat
          return
        end if
        allocate(item, stat=stat)
        if ( stat /= 0 ) then
          call release_reservation(bytes)
          enqueue_rkx_as_r4_device = async_err
          return
        end if
        item%op = async_write_r4
        item%ncid = ncid
        item%varid = varid
        item%nvals = nvals
        item%bytes = bytes
        call copy_bounds(item,start,count)
        call allocate_item_r4(item,nvals,bytes,stat)
        if ( stat /= async_ok ) then
          call release_reservation(bytes)
          call free_item(item)
          enqueue_rkx_as_r4_device = stat
          return
        end if
        call copy_device_to_item_rkx_as_r4(item,values,nvals,stat)
        if ( stat /= async_ok ) then
          call release_reservation(bytes)
          call free_item(item)
          enqueue_rkx_as_r4_device = stat
          return
        end if
        call enqueue_item(item)
        enqueue_rkx_as_r4_device = async_ok
      end if
    end function enqueue_rkx_as_r4_device

    integer(ik4) function enqueue_read_r4(ncid,varid,values,start,count, &
        counter)
      implicit none
      integer(ik4), intent(in) :: ncid, varid
      real(rk4), contiguous, dimension(:), intent(inout), target :: values
      integer(ik4), dimension(:), intent(in) :: start, count
      type(async_netcdf_counter), target, intent(inout) :: counter
      type(async_item), pointer :: item
      integer(ik4) :: nvals, stat

      call async_init()
      nvals = size(values)
      if ( first_error /= async_ok ) then
        enqueue_read_r4 = first_error
      else if ( .not. enabled ) then
        enqueue_read_r4 = async_err
      else
        allocate(item, stat=stat)
        if ( stat /= 0 ) then
          enqueue_read_r4 = async_err
          return
        end if
        item%op = async_read_r4
        item%ncid = ncid
        item%varid = varid
        item%nvals = nvals
        item%external_buffer = .true.
        item%speculative = .true.
        item%counter => counter
        item%r4 => values
        call copy_bounds(item,start,count)
        call enqueue_read_item(item)
        enqueue_read_r4 = async_ok
      end if
    end function enqueue_read_r4

    integer(ik4) function enqueue_read_r8(ncid,varid,values,start,count, &
        counter)
      implicit none
      integer(ik4), intent(in) :: ncid, varid
      real(rk8), contiguous, dimension(:), intent(inout), target :: values
      integer(ik4), dimension(:), intent(in) :: start, count
      type(async_netcdf_counter), target, intent(inout) :: counter
      type(async_item), pointer :: item
      integer(ik4) :: nvals, stat

      call async_init()
      nvals = size(values)
      if ( first_error /= async_ok ) then
        enqueue_read_r8 = first_error
      else if ( .not. enabled ) then
        enqueue_read_r8 = async_err
      else
        allocate(item, stat=stat)
        if ( stat /= 0 ) then
          enqueue_read_r8 = async_err
          return
        end if
        item%op = async_read_r8
        item%ncid = ncid
        item%varid = varid
        item%nvals = nvals
        item%external_buffer = .true.
        item%speculative = .true.
        item%counter => counter
        item%r8 => values
        call copy_bounds(item,start,count)
        call enqueue_read_item(item)
        enqueue_read_r8 = async_ok
      end if
    end function enqueue_read_r8

    subroutine async_init(limit_bytes)
      implicit none
      integer(c_int64_t), intent(in), optional :: limit_bytes
      integer(c_int) :: rc
      integer(ik4) :: stat

      if ( initialized ) return
      call read_memory_cap()
      if ( enabled .and. present(limit_bytes) ) then
        if ( limit_bytes <= 0_c_int64_t ) then
          enabled = .false.
          max_bytes = 0_c_int64_t
        else
          max_bytes = min(max_bytes,limit_bytes)
        end if
      end if
      stop_requested = .false.
      defer_work_signals = 0
      deferred_work_pending = .false.
#ifdef OPENACC
      call reset_pool_state()
      if ( enabled ) then
        call initialize_pool(stat)
        if ( stat /= async_ok ) then
          enabled = .false.
          first_error = async_err
        end if
      end if
#endif
      if ( enabled ) then
        rc = regcm_async_thread_start()
        if ( rc /= 0 ) then
          enabled = .false.
          first_error = async_err
        else
          worker_running = .true.
        end if
      end if
      initialized = .true.
    end subroutine async_init

    subroutine read_memory_cap()
      implicit none
      character(len=64) :: env_value
      integer(ik4) :: env_len, env_stat, read_stat
      real(rk8) :: cap_gb

#ifndef ASYNC_NETCDF
      env_max_bytes = 0_c_int64_t
      max_bytes = 0_c_int64_t
      enabled = .false.
      memory_cap_read = .true.
      return
#endif

      if ( .not. memory_cap_read ) then
        cap_gb = 0.0_rk8
        call get_environment_variable('RCM_ASYNC_OUTPUT_GB',env_value, &
          length=env_len,status=env_stat)
        if ( env_stat == 0 .and. env_len > 0 ) then
          read(env_value(1:env_len),*,iostat=read_stat) cap_gb
          if ( read_stat /= 0 ) cap_gb = 0.0_rk8
        end if
        if ( cap_gb <= 0.0_rk8 ) then
          env_max_bytes = 0_c_int64_t
        else
          env_max_bytes = int(cap_gb*1024.0_rk8*1024.0_rk8* &
            1024.0_rk8,c_int64_t)
        end if
        memory_cap_read = .true.
      end if
      max_bytes = env_max_bytes
      if ( max_bytes <= 0_c_int64_t ) then
        enabled = .false.
      else
        enabled = .true.
      end if
    end subroutine read_memory_cap

    subroutine reserve_memory(bytes,stat)
      implicit none
      integer(c_int64_t), intent(in) :: bytes
      integer(ik4), intent(out) :: stat

      stat = async_ok
      call regcm_async_queue_lock()
      do while ( first_error == async_ok .and. max_bytes > 0_c_int64_t .and. &
                 reserved_bytes > 0_c_int64_t .and. &
                 reserved_bytes+bytes > max_bytes )
        call release_deferred_work_locked()
        call regcm_async_wait_space()
      end do
      if ( first_error /= async_ok ) then
        stat = first_error
      else
        reserved_bytes = reserved_bytes + bytes
      end if
      call regcm_async_queue_unlock()
    end subroutine reserve_memory

    subroutine release_reservation(bytes)
      implicit none
      integer(c_int64_t), intent(in) :: bytes

      call regcm_async_queue_lock()
      reserved_bytes = reserved_bytes - bytes
      call regcm_async_broadcast_space()
      call regcm_async_queue_unlock()
    end subroutine release_reservation

    subroutine enqueue_item(item)
      implicit none
      type(async_item), pointer, intent(inout) :: item

      call regcm_async_queue_lock()
      call append_item(item)
      call signal_work_locked()
      call regcm_async_queue_unlock()
    end subroutine enqueue_item

    subroutine enqueue_read_item(item)
      implicit none
      type(async_item), pointer, intent(inout) :: item

      call regcm_async_queue_lock()
      if ( associated(item%counter) ) then
        item%counter%pending = item%counter%pending + 1
      end if
      call append_read_item(item)
      call regcm_async_signal_work()
      call regcm_async_queue_unlock()
    end subroutine enqueue_read_item

    subroutine signal_work_locked()
      implicit none

      if ( defer_work_signals > 0 ) then
        deferred_work_pending = .true.
      else
        call regcm_async_signal_work()
      end if
    end subroutine signal_work_locked

    subroutine release_deferred_work_locked()
      implicit none

      if ( associated(read_queue_head) ) then
        deferred_work_pending = .false.
        call regcm_async_signal_work()
      else if ( deferred_work_pending .and. associated(queue_head) ) then
        deferred_work_pending = .false.
        call regcm_async_signal_work()
      else if ( .not. associated(queue_head) ) then
        deferred_work_pending = .false.
      end if
    end subroutine release_deferred_work_locked

    subroutine append_item(item)
      implicit none
      type(async_item), pointer, intent(inout) :: item

      item%next => null()
      if ( associated(queue_tail) ) then
        queue_tail%next => item
      else
        queue_head => item
      end if
      queue_tail => item
    end subroutine append_item

    subroutine append_read_item(item)
      implicit none
      type(async_item), pointer, intent(inout) :: item

      item%next => null()
      if ( associated(read_queue_tail) ) then
        read_queue_tail%next => item
      else
        read_queue_head => item
      end if
      read_queue_tail => item
    end subroutine append_read_item

    subroutine complete_counter_locked(counter,status)
      implicit none
      type(async_netcdf_counter), pointer, intent(inout) :: counter
      integer(ik4), intent(in) :: status

      if ( status /= async_ok .and. counter%status == async_ok ) then
        counter%status = status
      end if
      if ( counter%pending > 0 ) counter%pending = counter%pending - 1
      call regcm_async_broadcast_space()
    end subroutine complete_counter_locked

    subroutine copy_bounds(item,start,count)
      implicit none
      type(async_item), intent(inout) :: item
      integer(ik4), dimension(:), intent(in), optional :: start, count
      integer(ik4) :: nd

      if ( present(start) ) then
        nd = min(size(start),size(item%start))
        item%start(1:nd) = start(1:nd)
        item%nd = nd
        item%has_start = .true.
      end if
      if ( present(count) ) then
        nd = min(size(count),size(item%count))
        item%count(1:nd) = count(1:nd)
        if ( item%nd == 0 ) item%nd = nd
        item%has_count = .true.
      end if
    end subroutine copy_bounds

    integer(c_int64_t) function align_pool_bytes(bytes)
      implicit none
      integer(c_int64_t), intent(in) :: bytes

      align_pool_bytes = ((bytes+pool_alignment-1_c_int64_t) / &
        pool_alignment) * pool_alignment
    end function align_pool_bytes

    integer(c_int64_t) function bytes_rkx()
      implicit none

#ifdef SINGLE_PRECISION_REAL
      bytes_rkx = bytes_r4
#else
      bytes_rkx = bytes_r8
#endif
    end function bytes_rkx

#ifdef OPENACC
    subroutine reset_pool_state()
      implicit none

      call free_pool_free_list()
      pool_cptr = c_null_ptr
      pool_capacity = 0_c_int64_t
      pool_used = 0_c_int64_t
    end subroutine reset_pool_state

    subroutine initialize_pool(stat)
      implicit none
      integer(ik4), intent(out) :: stat
      integer :: alloc_stat
      integer(c_int) :: istat

      stat = async_ok
      if ( max_bytes <= 0_c_int64_t ) return
      istat = cudaMallocHost(pool_cptr,int(max_bytes,c_size_t))
      if ( istat == cudaSuccess ) then
        pool_capacity = max_bytes
        pool_used = 0_c_int64_t
        allocate(pool_free_head,stat=alloc_stat)
        if ( alloc_stat == 0 ) then
          pool_free_head%offset = 0_c_int64_t
          pool_free_head%bytes = pool_capacity
          pool_free_head%next => null()
        else
          istat = cudaFreeHost(pool_cptr)
          pool_cptr = c_null_ptr
          pool_capacity = 0_c_int64_t
          pool_used = 0_c_int64_t
          stat = async_err
        end if
      else
        stat = async_err
      end if
    end subroutine initialize_pool

    subroutine allocate_item_buffer(item,bytes,stat)
      implicit none
      type(async_item), intent(inout) :: item
      integer(c_int64_t), intent(in) :: bytes
      integer(ik4), intent(out) :: stat
      integer(c_int) :: istat
      integer(c_int64_t) :: aligned_bytes

      aligned_bytes = align_pool_bytes(bytes)
      if ( enabled .and. c_associated(pool_cptr) .and. &
           aligned_bytes <= pool_capacity ) then
        call acquire_pool_buffer(item,aligned_bytes,stat)
      else
        istat = cudaMallocHost(item%cptr,int(bytes,c_size_t))
        if ( istat == cudaSuccess ) then
          item%buffer_bytes = bytes
          item%pooled_buffer = .false.
          stat = async_ok
        else
          stat = async_err
        end if
      end if
    end subroutine allocate_item_buffer

    subroutine acquire_pool_buffer(item,bytes,stat)
      implicit none
      type(async_item), intent(inout) :: item
      integer(c_int64_t), intent(in) :: bytes
      integer(ik4), intent(out) :: stat
      logical :: have_buffer

      stat = async_ok
      have_buffer = .false.

      call regcm_async_queue_lock()
      do while ( first_error == async_ok )
        call take_pool_slice_locked(bytes,item%pool_offset,have_buffer)
        if ( have_buffer ) exit
        call release_deferred_work_locked()
        call regcm_async_wait_space()
      end do
      if ( first_error /= async_ok ) then
        stat = first_error
      else
        item%buffer_bytes = bytes
        item%cptr = regcm_async_ptr_offset(pool_cptr,item%pool_offset)
        item%pooled_buffer = .true.
      end if
      call regcm_async_queue_unlock()
    end subroutine acquire_pool_buffer

    subroutine take_pool_slice_locked(bytes,offset,found)
      implicit none
      integer(c_int64_t), intent(in) :: bytes
      integer(c_int64_t), intent(out) :: offset
      logical, intent(out) :: found
      type(pool_free_block), pointer :: block, prev

      offset = 0_c_int64_t
      found = .false.

      prev => null()
      block => pool_free_head
      do while ( associated(block) )
        if ( block%bytes >= bytes ) then
          offset = block%offset
          block%offset = block%offset + bytes
          block%bytes = block%bytes - bytes
          if ( block%bytes == 0_c_int64_t ) then
            if ( associated(prev) ) then
              prev%next => block%next
            else
              pool_free_head => block%next
            end if
            deallocate(block)
          end if
          pool_used = pool_used + bytes
          found = .true.
          return
        end if
        prev => block
        block => block%next
      end do
    end subroutine take_pool_slice_locked

    subroutine release_item_buffer(cptr,offset,buffer_bytes)
      implicit none
      type(c_ptr), intent(inout) :: cptr
      integer(c_int64_t), intent(in) :: offset, buffer_bytes

      if ( .not. c_associated(cptr) ) return
      call regcm_async_queue_lock()
      call release_pool_slice_locked(offset,buffer_bytes)
      call regcm_async_broadcast_space()
      call regcm_async_queue_unlock()
      cptr = c_null_ptr
    end subroutine release_item_buffer

    subroutine release_pool_slice_locked(offset,bytes)
      implicit none
      integer(c_int64_t), intent(in) :: offset, bytes
      type(pool_free_block), pointer :: block, prev, next_block
      integer :: alloc_stat

      if ( pool_used >= bytes ) then
        pool_used = pool_used - bytes
      else
        pool_used = 0_c_int64_t
      end if

      allocate(block,stat=alloc_stat)
      if ( alloc_stat /= 0 ) then
        if ( first_error == async_ok ) first_error = async_err
        return
      end if
      block%offset = offset
      block%bytes = bytes
      block%next => null()

      prev => null()
      next_block => pool_free_head
      do while ( associated(next_block) .and. next_block%offset < offset )
        prev => next_block
        next_block => next_block%next
      end do
      block%next => next_block
      if ( associated(prev) ) then
        prev%next => block
      else
        pool_free_head => block
      end if

      if ( associated(block%next) .and. &
           block%offset + block%bytes == block%next%offset ) then
        next_block => block%next
        block%bytes = block%bytes + next_block%bytes
        block%next => next_block%next
        deallocate(next_block)
      end if
      if ( associated(prev) .and. &
           prev%offset + prev%bytes == block%offset ) then
        prev%bytes = prev%bytes + block%bytes
        prev%next => block%next
        deallocate(block)
      end if
    end subroutine release_pool_slice_locked

    subroutine free_pool_buffers()
      implicit none
      integer(c_int) :: istat

      call regcm_async_queue_lock()
      pool_capacity = 0_c_int64_t
      pool_used = 0_c_int64_t
      call free_pool_free_list()
      call regcm_async_broadcast_space()
      call regcm_async_queue_unlock()

      if ( c_associated(pool_cptr) ) istat = cudaFreeHost(pool_cptr)
      pool_cptr = c_null_ptr
    end subroutine free_pool_buffers

    subroutine free_pool_free_list()
      implicit none
      type(pool_free_block), pointer :: block

      do while ( associated(pool_free_head) )
        block => pool_free_head
        pool_free_head => pool_free_head%next
        deallocate(block)
      end do
    end subroutine free_pool_free_list
#endif

    subroutine allocate_item_r4(item,nvals,bytes,stat)
      implicit none
      type(async_item), intent(inout) :: item
      integer(ik4), intent(in) :: nvals
      integer(c_int64_t), intent(in) :: bytes
      integer(ik4), intent(out) :: stat
#ifndef OPENACC
      integer(ik4) :: istat
#endif

#ifdef OPENACC
      call allocate_item_buffer(item,bytes,stat)
      if ( stat == async_ok ) then
        call c_f_pointer(item%cptr,item%r4,[nvals])
      end if
#else
      allocate(item%r4(nvals),stat=istat)
      if ( istat == 0 ) then
        stat = async_ok
      else
        stat = async_err
      end if
#endif
    end subroutine allocate_item_r4

    subroutine allocate_item_r8(item,nvals,bytes,stat)
      implicit none
      type(async_item), intent(inout) :: item
      integer(ik4), intent(in) :: nvals
      integer(c_int64_t), intent(in) :: bytes
      integer(ik4), intent(out) :: stat
#ifndef OPENACC
      integer(ik4) :: istat
#endif

#ifdef OPENACC
      call allocate_item_buffer(item,bytes,stat)
      if ( stat == async_ok ) then
        call c_f_pointer(item%cptr,item%r8,[nvals])
      end if
#else
      allocate(item%r8(nvals),stat=istat)
      if ( istat == 0 ) then
        stat = async_ok
      else
        stat = async_err
      end if
#endif
    end subroutine allocate_item_r8

    subroutine allocate_item_i4(item,nvals,bytes,stat)
      implicit none
      type(async_item), intent(inout) :: item
      integer(ik4), intent(in) :: nvals
      integer(c_int64_t), intent(in) :: bytes
      integer(ik4), intent(out) :: stat
#ifndef OPENACC
      integer(ik4) :: istat
#endif

#ifdef OPENACC
      call allocate_item_buffer(item,bytes,stat)
      if ( stat == async_ok ) then
        call c_f_pointer(item%cptr,item%i4,[nvals])
      end if
#else
      allocate(item%i4(nvals),stat=istat)
      if ( istat == 0 ) then
        stat = async_ok
      else
        stat = async_err
      end if
#endif
    end subroutine allocate_item_i4

    subroutine copy_device_to_item_r4(item,values,nvals,stat)
      implicit none
      type(async_item), intent(inout) :: item
      real(rk4), contiguous, dimension(:), intent(in) :: values
      integer(ik4), intent(in) :: nvals
      integer(ik4), intent(out) :: stat
#ifdef OPENACC
      integer(c_int) :: istat

      istat = cudaMemcpy(item%r4,values,nvals,cudaMemcpyDeviceToHost)
      if ( istat == cudaSuccess ) then
        stat = async_ok
      else
        stat = async_err
      end if
#else
      item%r4(1:nvals) = values(1:nvals)
      stat = async_ok
#endif
    end subroutine copy_device_to_item_r4

    subroutine copy_device_to_item_r8(item,values,nvals,stat)
      implicit none
      type(async_item), intent(inout) :: item
      real(rk8), contiguous, dimension(:), intent(in) :: values
      integer(ik4), intent(in) :: nvals
      integer(ik4), intent(out) :: stat
#ifdef OPENACC
      integer(c_int) :: istat

      istat = cudaMemcpy(item%r8,values,nvals,cudaMemcpyDeviceToHost)
      if ( istat == cudaSuccess ) then
        stat = async_ok
      else
        stat = async_err
      end if
#else
      item%r8(1:nvals) = values(1:nvals)
      stat = async_ok
#endif
    end subroutine copy_device_to_item_r8

    subroutine copy_device_to_item_i4(item,values,nvals,stat)
      implicit none
      type(async_item), intent(inout) :: item
      integer(ik4), contiguous, dimension(:), intent(in) :: values
      integer(ik4), intent(in) :: nvals
      integer(ik4), intent(out) :: stat
#ifdef OPENACC
      integer(c_int) :: istat

      istat = cudaMemcpy(item%i4,values,nvals,cudaMemcpyDeviceToHost)
      if ( istat == cudaSuccess ) then
        stat = async_ok
      else
        stat = async_err
      end if
#else
      item%i4(1:nvals) = values(1:nvals)
      stat = async_ok
#endif
    end subroutine copy_device_to_item_i4

    subroutine copy_device_to_item_rkx_as_r4(item,values,nvals,stat)
      implicit none
      type(async_item), intent(inout) :: item
      real(rkx), contiguous, dimension(:), intent(in), target :: values
      integer(ik4), intent(in) :: nvals
      integer(ik4), intent(out) :: stat
#ifdef OPENACC
      real(rk4), pointer, contiguous, dimension(:) :: r4_1d => null()
      real(rkx), pointer, contiguous, dimension(:) :: rkx_1d => null()
      integer(ik4) :: i

      r4_1d(1:nvals) => item%r4(1:nvals)
      rkx_1d(1:nvals) => values(1:nvals)
      !$acc parallel loop deviceptr(r4_1d)
      do i = 1, nvals
        r4_1d(i) = real(rkx_1d(i),rk4)
      end do
      !$acc end parallel loop
      stat = async_ok
#else
      item%r4(1:nvals) = real(values(1:nvals),rk4)
      stat = async_ok
#endif
    end subroutine copy_device_to_item_rkx_as_r4

    integer(ik4) function enqueue_r4_device_sync(ncid,varid,values,start,count)
      implicit none
      integer(ik4), intent(in) :: ncid, varid
      real(rk4), contiguous, dimension(:), intent(in) :: values
      integer(ik4), dimension(:), intent(in), optional :: start, count
      type(async_item), pointer :: item
      integer(ik4) :: nvals, stat
      integer(c_int64_t) :: bytes

      nvals = size(values)
      bytes = int(nvals,c_int64_t)*bytes_r4
      allocate(item, stat=stat)
      if ( stat /= 0 ) then
        enqueue_r4_device_sync = async_err
        return
      end if
      item%op = async_write_r4
      item%ncid = ncid
      item%varid = varid
      item%nvals = nvals
      item%bytes = bytes
      call copy_bounds(item,start,count)
      call allocate_item_r4(item,nvals,bytes,stat)
      if ( stat == async_ok ) call copy_device_to_item_r4(item,values,nvals,stat)
      if ( stat == async_ok ) then
        enqueue_r4_device_sync = put_var_r4(item)
      else
        enqueue_r4_device_sync = stat
      end if
      call free_item(item)
    end function enqueue_r4_device_sync

    integer(ik4) function enqueue_r8_device_sync(ncid,varid,values,start,count)
      implicit none
      integer(ik4), intent(in) :: ncid, varid
      real(rk8), contiguous, dimension(:), intent(in) :: values
      integer(ik4), dimension(:), intent(in), optional :: start, count
      type(async_item), pointer :: item
      integer(ik4) :: nvals, stat
      integer(c_int64_t) :: bytes

      nvals = size(values)
      bytes = int(nvals,c_int64_t)*bytes_r8
      allocate(item, stat=stat)
      if ( stat /= 0 ) then
        enqueue_r8_device_sync = async_err
        return
      end if
      item%op = async_write_r8
      item%ncid = ncid
      item%varid = varid
      item%nvals = nvals
      item%bytes = bytes
      call copy_bounds(item,start,count)
      call allocate_item_r8(item,nvals,bytes,stat)
      if ( stat == async_ok ) call copy_device_to_item_r8(item,values,nvals,stat)
      if ( stat == async_ok ) then
        enqueue_r8_device_sync = put_var_r8(item)
      else
        enqueue_r8_device_sync = stat
      end if
      call free_item(item)
    end function enqueue_r8_device_sync

    integer(ik4) function enqueue_i4_device_sync(ncid,varid,values,start,count)
      implicit none
      integer(ik4), intent(in) :: ncid, varid
      integer(ik4), contiguous, dimension(:), intent(in) :: values
      integer(ik4), dimension(:), intent(in), optional :: start, count
      type(async_item), pointer :: item
      integer(ik4) :: nvals, stat
      integer(c_int64_t) :: bytes

      nvals = size(values)
      bytes = int(nvals,c_int64_t)*bytes_i4
      allocate(item, stat=stat)
      if ( stat /= 0 ) then
        enqueue_i4_device_sync = async_err
        return
      end if
      item%op = async_write_i4
      item%ncid = ncid
      item%varid = varid
      item%nvals = nvals
      item%bytes = bytes
      call copy_bounds(item,start,count)
      call allocate_item_i4(item,nvals,bytes,stat)
      if ( stat == async_ok ) call copy_device_to_item_i4(item,values,nvals,stat)
      if ( stat == async_ok ) then
        enqueue_i4_device_sync = put_var_i4(item)
      else
        enqueue_i4_device_sync = stat
      end if
      call free_item(item)
    end function enqueue_i4_device_sync

    integer(ik4) function enqueue_rkx_as_r4_device_sync(ncid,varid,values, &
        start,count)
      implicit none
      integer(ik4), intent(in) :: ncid, varid
      real(rkx), contiguous, dimension(:), intent(in), target :: values
      integer(ik4), dimension(:), intent(in), optional :: start, count
      type(async_item), pointer :: item
      integer(ik4) :: nvals, stat
      integer(c_int64_t) :: bytes

      nvals = size(values)
      bytes = int(nvals,c_int64_t)*bytes_r4
      allocate(item, stat=stat)
      if ( stat /= 0 ) then
        enqueue_rkx_as_r4_device_sync = async_err
        return
      end if
      item%op = async_write_r4
      item%ncid = ncid
      item%varid = varid
      item%nvals = nvals
      item%bytes = bytes
      call copy_bounds(item,start,count)
      call allocate_item_r4(item,nvals,bytes,stat)
      if ( stat == async_ok ) call copy_device_to_item_rkx_as_r4(item,values, &
        nvals,stat)
      if ( stat == async_ok ) then
        enqueue_rkx_as_r4_device_sync = put_var_r4(item)
      else
        enqueue_rkx_as_r4_device_sync = stat
      end if
      call free_item(item)
    end function enqueue_rkx_as_r4_device_sync

    integer(ik4) function execute_item(item)
      implicit none
      type(async_item), pointer, intent(in) :: item

      select case ( item%op )
        case ( async_write_r4 )
          execute_item = put_var_r4(item)
        case ( async_write_r8 )
          execute_item = put_var_r8(item)
        case ( async_write_i4 )
          execute_item = put_var_i4(item)
        case ( async_sync )
          execute_item = nf90_sync(item%ncid)
        case ( async_close )
          execute_item = nf90_close(item%ncid)
        case ( async_read_r4 )
          execute_item = get_var_r4(item)
        case ( async_read_r8 )
          execute_item = get_var_r8(item)
        case default
          execute_item = async_err
      end select
    end function execute_item

    integer(ik4) function put_var_r4(item)
      implicit none
      type(async_item), pointer, intent(in) :: item

      if ( item%scalar ) then
        if ( item%has_start .and. item%has_count ) then
          put_var_r4 = nf90_put_var(item%ncid,item%varid,item%r4(1:1), &
            start=item%start(1:item%nd),count=item%count(1:item%nd))
        else if ( item%has_start ) then
          put_var_r4 = nf90_put_var(item%ncid,item%varid,item%r4(1:1), &
            start=item%start(1:item%nd))
        else
          put_var_r4 = nf90_put_var(item%ncid,item%varid,item%r4(1))
        end if
      else
        if ( item%has_start .and. item%has_count ) then
          put_var_r4 = nf90_put_var(item%ncid,item%varid,item%r4, &
            start=item%start(1:item%nd),count=item%count(1:item%nd))
        else if ( item%has_start ) then
          put_var_r4 = nf90_put_var(item%ncid,item%varid,item%r4, &
            start=item%start(1:item%nd))
        else
          put_var_r4 = nf90_put_var(item%ncid,item%varid,item%r4)
        end if
      end if
    end function put_var_r4

    integer(ik4) function put_var_r8(item)
      implicit none
      type(async_item), pointer, intent(in) :: item

      if ( item%scalar ) then
        if ( item%has_start .and. item%has_count ) then
          put_var_r8 = nf90_put_var(item%ncid,item%varid,item%r8(1:1), &
            start=item%start(1:item%nd),count=item%count(1:item%nd))
        else if ( item%has_start ) then
          put_var_r8 = nf90_put_var(item%ncid,item%varid,item%r8(1:1), &
            start=item%start(1:item%nd))
        else
          put_var_r8 = nf90_put_var(item%ncid,item%varid,item%r8(1))
        end if
      else
        if ( item%has_start .and. item%has_count ) then
          put_var_r8 = nf90_put_var(item%ncid,item%varid,item%r8, &
            start=item%start(1:item%nd),count=item%count(1:item%nd))
        else if ( item%has_start ) then
          put_var_r8 = nf90_put_var(item%ncid,item%varid,item%r8, &
            start=item%start(1:item%nd))
        else
          put_var_r8 = nf90_put_var(item%ncid,item%varid,item%r8)
        end if
      end if
    end function put_var_r8

    integer(ik4) function put_var_i4(item)
      implicit none
      type(async_item), pointer, intent(in) :: item

      if ( item%scalar ) then
        if ( item%has_start .and. item%has_count ) then
          put_var_i4 = nf90_put_var(item%ncid,item%varid,item%i4(1:1), &
            start=item%start(1:item%nd),count=item%count(1:item%nd))
        else if ( item%has_start ) then
          put_var_i4 = nf90_put_var(item%ncid,item%varid,item%i4(1:1), &
            start=item%start(1:item%nd))
        else
          put_var_i4 = nf90_put_var(item%ncid,item%varid,item%i4(1))
        end if
      else
        if ( item%has_start .and. item%has_count ) then
          put_var_i4 = nf90_put_var(item%ncid,item%varid,item%i4, &
            start=item%start(1:item%nd),count=item%count(1:item%nd))
        else if ( item%has_start ) then
          put_var_i4 = nf90_put_var(item%ncid,item%varid,item%i4, &
            start=item%start(1:item%nd))
        else
          put_var_i4 = nf90_put_var(item%ncid,item%varid,item%i4)
        end if
      end if
    end function put_var_i4

    integer(ik4) function get_var_r4(item)
      implicit none
      type(async_item), pointer, intent(in) :: item

      if ( item%has_start .and. item%has_count ) then
        get_var_r4 = nf90_get_var(item%ncid,item%varid,item%r4, &
          start=item%start(1:item%nd),count=item%count(1:item%nd))
      else if ( item%has_start ) then
        get_var_r4 = nf90_get_var(item%ncid,item%varid,item%r4, &
          start=item%start(1:item%nd))
      else
        get_var_r4 = nf90_get_var(item%ncid,item%varid,item%r4)
      end if
    end function get_var_r4

    integer(ik4) function get_var_r8(item)
      implicit none
      type(async_item), pointer, intent(in) :: item

      if ( item%has_start .and. item%has_count ) then
        get_var_r8 = nf90_get_var(item%ncid,item%varid,item%r8, &
          start=item%start(1:item%nd),count=item%count(1:item%nd))
      else if ( item%has_start ) then
        get_var_r8 = nf90_get_var(item%ncid,item%varid,item%r8, &
          start=item%start(1:item%nd))
      else
        get_var_r8 = nf90_get_var(item%ncid,item%varid,item%r8)
      end if
    end function get_var_r8

    integer(ik4) function put_var_r4_sync(ncid,varid,values,scalar,start,count)
      implicit none
      integer(ik4), intent(in) :: ncid, varid
      real(rk4), contiguous, dimension(:), intent(in) :: values
      logical, intent(in) :: scalar
      integer(ik4), dimension(:), intent(in), optional :: start, count

      put_var_r4_sync = put_var_r4_direct(ncid,varid,values,scalar,start,count)
    end function put_var_r4_sync

    integer(ik4) function put_var_r8_sync(ncid,varid,values,scalar,start,count)
      implicit none
      integer(ik4), intent(in) :: ncid, varid
      real(rk8), contiguous, dimension(:), intent(in) :: values
      logical, intent(in) :: scalar
      integer(ik4), dimension(:), intent(in), optional :: start, count

      put_var_r8_sync = put_var_r8_direct(ncid,varid,values,scalar,start,count)
    end function put_var_r8_sync

    integer(ik4) function put_var_i4_sync(ncid,varid,values,scalar,start,count)
      implicit none
      integer(ik4), intent(in) :: ncid, varid
      integer(ik4), contiguous, dimension(:), intent(in) :: values
      logical, intent(in) :: scalar
      integer(ik4), dimension(:), intent(in), optional :: start, count

      put_var_i4_sync = put_var_i4_direct(ncid,varid,values,scalar,start,count)
    end function put_var_i4_sync

    integer(ik4) function put_var_r4_direct(ncid,varid,values,scalar,start,count)
      implicit none
      integer(ik4), intent(in) :: ncid, varid
      real(rk4), contiguous, dimension(:), intent(in) :: values
      logical, intent(in) :: scalar
      integer(ik4), dimension(:), intent(in), optional :: start, count

      if ( scalar ) then
        if ( present(start) .and. present(count) ) then
          put_var_r4_direct = nf90_put_var(ncid,varid,values(1:1),start,count)
        else if ( present(start) ) then
          put_var_r4_direct = nf90_put_var(ncid,varid,values(1:1),start)
        else
          put_var_r4_direct = nf90_put_var(ncid,varid,values(1))
        end if
      else
        if ( present(start) .and. present(count) ) then
          put_var_r4_direct = nf90_put_var(ncid,varid,values,start,count)
        else if ( present(start) ) then
          put_var_r4_direct = nf90_put_var(ncid,varid,values,start)
        else
          put_var_r4_direct = nf90_put_var(ncid,varid,values)
        end if
      end if
    end function put_var_r4_direct

    integer(ik4) function put_var_r8_direct(ncid,varid,values,scalar,start,count)
      implicit none
      integer(ik4), intent(in) :: ncid, varid
      real(rk8), contiguous, dimension(:), intent(in) :: values
      logical, intent(in) :: scalar
      integer(ik4), dimension(:), intent(in), optional :: start, count

      if ( scalar ) then
        if ( present(start) .and. present(count) ) then
          put_var_r8_direct = nf90_put_var(ncid,varid,values(1:1),start,count)
        else if ( present(start) ) then
          put_var_r8_direct = nf90_put_var(ncid,varid,values(1:1),start)
        else
          put_var_r8_direct = nf90_put_var(ncid,varid,values(1))
        end if
      else
        if ( present(start) .and. present(count) ) then
          put_var_r8_direct = nf90_put_var(ncid,varid,values,start,count)
        else if ( present(start) ) then
          put_var_r8_direct = nf90_put_var(ncid,varid,values,start)
        else
          put_var_r8_direct = nf90_put_var(ncid,varid,values)
        end if
      end if
    end function put_var_r8_direct

    integer(ik4) function put_var_i4_direct(ncid,varid,values,scalar,start,count)
      implicit none
      integer(ik4), intent(in) :: ncid, varid
      integer(ik4), contiguous, dimension(:), intent(in) :: values
      logical, intent(in) :: scalar
      integer(ik4), dimension(:), intent(in), optional :: start, count

      if ( scalar ) then
        if ( present(start) .and. present(count) ) then
          put_var_i4_direct = nf90_put_var(ncid,varid,values(1:1),start,count)
        else if ( present(start) ) then
          put_var_i4_direct = nf90_put_var(ncid,varid,values(1:1),start)
        else
          put_var_i4_direct = nf90_put_var(ncid,varid,values(1))
        end if
      else
        if ( present(start) .and. present(count) ) then
          put_var_i4_direct = nf90_put_var(ncid,varid,values,start,count)
        else if ( present(start) ) then
          put_var_i4_direct = nf90_put_var(ncid,varid,values,start)
        else
          put_var_i4_direct = nf90_put_var(ncid,varid,values)
        end if
      end if
    end function put_var_i4_direct

    subroutine free_item(item)
      implicit none
      type(async_item), pointer, intent(inout) :: item
#ifdef OPENACC
      integer(c_int) :: istat
#endif

      if ( .not. associated(item) ) return
#ifdef OPENACC
      if ( c_associated(item%cptr) ) then
        if ( item%external_buffer ) then
          item%cptr = c_null_ptr
        else if ( item%pooled_buffer ) then
          call release_item_buffer(item%cptr,item%pool_offset, &
            item%buffer_bytes)
        else
          istat = cudaFreeHost(item%cptr)
          item%cptr = c_null_ptr
        end if
      end if
      if ( associated(item%r4) ) item%r4 => null()
      if ( associated(item%r8) ) item%r8 => null()
      if ( associated(item%i4) ) item%i4 => null()
#else
      if ( associated(item%r4) ) then
        if ( item%external_buffer ) then
          item%r4 => null()
        else
          deallocate(item%r4)
        end if
      end if
      if ( associated(item%r8) ) then
        if ( item%external_buffer ) then
          item%r8 => null()
        else
          deallocate(item%r8)
        end if
      end if
      if ( associated(item%i4) ) then
        if ( item%external_buffer ) then
          item%i4 => null()
        else
          deallocate(item%i4)
        end if
      end if
#endif
      deallocate(item)
    end subroutine free_item

end module mod_async_netcdf

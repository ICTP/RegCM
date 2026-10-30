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

! Unit test of Share/mod_async_netcdf.F90 through its public interface.
! Without ASYNC_NETCDF, or with RCM_ASYNC_OUTPUT_GB unset, every call is
! executed at once; with both, the writes go through the worker thread.
! Either way, records written with start and count, scalars, a sync and a
! close must give the same file, and the record reads must return it.
! Errors are checked last, as the first error stops the queue.

program test_async_netcdf
  use, intrinsic :: iso_c_binding
  use netcdf
  use mod_intkinds
  use mod_realkinds
  use mod_async_netcdf
  use share_test_utils
  implicit none
  character(len=*), parameter :: fname = 'test_async_netcdf.nc'
  integer(ik4), parameter :: nx = 50, nrec = 200
  integer(ik4) :: ncid, idx, idt, ivr4, ivr8, ivi4, ivsc, istat, n, nbad
  real(rk4), dimension(nx) :: r4
  real(rk8), dimension(nx) :: r8
  integer(ik4), dimension(nx) :: i4
  real(rk8) :: sc
  real(rkx), dimension(nx), target :: rr
  type(async_netcdf_counter) :: counter
  type(async_netcdf_buffer) :: buffer
  logical :: lasync

  lasync = async_netcdf_is_enabled()
  write(*,'(a,l1)') 'async worker enabled: ', lasync

  call async_netcdf_lock()
  istat = nf90_create(fname,ior(nf90_clobber,nf90_netcdf4),ncid)
  istat = nf90_def_dim(ncid,'x',nx,idx)
  istat = nf90_def_dim(ncid,'time',nf90_unlimited,idt)
  istat = nf90_def_var(ncid,'r4',nf90_float,[idx,idt],ivr4)
  istat = nf90_def_var(ncid,'r8',nf90_double,[idx,idt],ivr8)
  istat = nf90_def_var(ncid,'i4',nf90_int,[idx,idt],ivi4)
  istat = nf90_def_var(ncid,'sc',nf90_double,[idt],ivsc)
  istat = nf90_enddef(ncid)
  call async_netcdf_unlock()
  call check(istat == nf90_noerr,'file defined under the netCDF lock')

  nbad = 0
  call async_netcdf_defer_wake_begin()
  do n = 1, nrec
    call fill(n)
    if ( async_netcdf_put_var(ncid,ivr4,r4,[1,n],[nx,1]) /= 0 ) &
      nbad = nbad + 1
    if ( async_netcdf_put_var(ncid,ivr8,r8,[1,n],[nx,1]) /= 0 ) &
      nbad = nbad + 1
    if ( async_netcdf_put_var(ncid,ivi4,i4,[1,n],[nx,1]) /= 0 ) &
      nbad = nbad + 1
    if ( async_netcdf_put_var(ncid,ivsc,sc,[n]) /= 0 ) nbad = nbad + 1
    if ( n == nrec/2 ) then
      call async_netcdf_defer_wake_end()
      if ( async_netcdf_sync(ncid) /= 0 ) nbad = nbad + 1
    end if
  end do
  call check(nbad == 0,'writes accepted')
  call check(async_netcdf_wait_all() == nf90_noerr,'all writes completed')
  call check(async_netcdf_close(ncid) == nf90_noerr,'close')
  call check(async_netcdf_wait_all() == nf90_noerr,'close completed')

  ! The file is complete.
  istat = nf90_open(fname,nf90_nowrite,ncid)
  nbad = 0
  do n = 1, nrec
    call fill(n)
    if ( .not. same(n) ) nbad = nbad + 1
  end do
  call check(nbad == 0,'every record written in place')

  ! Reads of records into an external buffer.
  if ( lasync ) then
    call async_netcdf_counter_reset(counter)
    istat = async_netcdf_get_var_rkx(ncid,ivr8,rr,[1,7],[nx,1],counter)
    call async_netcdf_counter_wait(counter,istat)
    call fill(7)
    call check(istat == nf90_noerr .and. all(rr == real(r8,rkx)), &
      'record read by the worker')
    call async_netcdf_try_acquire_buffer(buffer,100_c_int64_t,istat)
    call check(istat == nf90_noerr .and. buffer%active .and. &
      size(buffer%rkx) >= 100/int(storage_size(1.0_rkx)/8), &
      'acquire a buffer')
    call async_netcdf_release_buffer(buffer)
    call check(.not. buffer%active,'release a buffer')
  else
    call async_netcdf_counter_reset(counter)
    istat = async_netcdf_get_var_rkx(ncid,ivr8,rr,[1,7],[nx,1],counter)
    call check(istat /= nf90_noerr,'no queued reads without the worker')
    call async_netcdf_try_acquire_buffer(buffer,100_c_int64_t,istat)
    call check(istat /= nf90_noerr .and. .not. buffer%active, &
      'no buffers without the worker')
  end if
  istat = nf90_close(ncid)

  ! A write to a closed file fails, at once or at the next wait.
  istat = async_netcdf_put_var(ncid,ivr4,r4,[1,1],[nx,1])
  if ( istat == nf90_noerr ) istat = async_netcdf_wait_all()
  call check(istat /= nf90_noerr,'error reported')
  istat = async_netcdf_shutdown()

  call test_end('test_async_netcdf')

  contains

  subroutine fill(n)
    implicit none
    integer(ik4), intent(in) :: n
    integer(ik4) :: i
    do i = 1, nx
      r4(i) = real(i+1000*n,rk4)
      r8(i) = real(i,rk8)/3.0_rk8 + real(n,rk8)
      i4(i) = i - n
    end do
    sc = real(n,rk8)*0.5_rk8
  end subroutine fill

  logical function same(n)
    implicit none
    integer(ik4), intent(in) :: n
    real(rk4), dimension(nx) :: x4
    real(rk8), dimension(nx) :: x8
    integer(ik4), dimension(nx) :: xi
    real(rk8), dimension(1) :: xs
    istat = nf90_get_var(ncid,ivr4,x4,[1,n],[nx,1])
    istat = nf90_get_var(ncid,ivr8,x8,[1,n],[nx,1])
    istat = nf90_get_var(ncid,ivi4,xi,[1,n],[nx,1])
    istat = nf90_get_var(ncid,ivsc,xs,[n],[1])
    same = all(x4 == r4) .and. all(x8 == r8) .and. all(xi == i4) .and. &
           xs(1) == sc
  end function same

end program test_async_netcdf

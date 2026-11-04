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

! Counter based random numbers.
!
! The value is a hash of a stream seed, a counter (the time step) and the
! global grid indices of the point. It does not depend on the domain
! decomposition, on the order in which points are visited, or on the
! history of the run, so a restarted run draws the same numbers as a
! continuous one. The functions are pure and keep no state, so they can
! be called inside do concurrent loops on CPU and GPU.
!
! The seed and the counter are hashed once into a key with hash_key,
! outside the loop on the points; hash_uniform then mixes in the three
! indices of each point.
!
! The hash is the lowbias32 integer mixer (C. Wellons, "Prospecting for
! hash functions", 2018) applied in chain to each input. Unsigned 32 bit
! arithmetic is emulated in 64 bit integers, with the products split in
! 16 bit halves so that no intermediate value overflows.

module mod_hashrng

  use mod_intkinds
  use mod_realkinds

  implicit none

  private

  integer(ik8), parameter :: mask32 = 4294967295_ik8
  integer(ik8), parameter :: mask16 = 65535_ik8
  integer(ik8), parameter :: hmul1 = 2146121005_ik8 ! 0x7feb352d
  integer(ik8), parameter :: hmul2 = 2221713035_ik8 ! 0x846ca68b
  real(rk8), parameter :: two_m24 = 1.0_rk8/16777216.0_rk8

  public :: hash_key, hash_uniform, hash_uint32

  contains

  ! x * c modulo 2**32, with x and c in [0,2**32)
  pure integer(ik8) function mul32(x,c)
!$acc routine seq
    implicit none
    integer(ik8), intent(in) :: x, c
    mul32 = iand(x*iand(c,mask16) + &
                 ishft(iand(x*ishft(c,-16),mask16),16),mask32)
  end function mul32

  pure integer(ik8) function lowbias32(x)
!$acc routine seq
    implicit none
    integer(ik8), intent(in) :: x
    integer(ik8) :: h
    h = iand(x,mask32)
    h = ieor(h,ishft(h,-16))
    h = mul32(h,hmul1)
    h = ieor(h,ishft(h,-15))
    h = mul32(h,hmul2)
    lowbias32 = ieor(h,ishft(h,-16))
  end function lowbias32

  ! The key of a stream (seed) at a counter value (time step), in
  ! [0,2**32).
  pure integer(ik8) function hash_key(seed,n)
!$acc routine seq
    implicit none
    integer(ik4), intent(in) :: seed
    integer(ik8), intent(in) :: n
    integer(ik8) :: h
    h = lowbias32(int(seed,ik8))
    h = lowbias32(h + iand(n,mask32))
    hash_key = lowbias32(h + iand(ishft(n,-32),mask32))
  end function hash_key

  ! A 32 bit hash in [0,2**32) of a key and three indices.
  pure integer(ik8) function hash_uint32(key,k,i,j)
!$acc routine seq
    implicit none
    integer(ik8), intent(in) :: key
    integer(ik4), intent(in) :: k, i, j
    integer(ik8) :: h
    h = lowbias32(key + iand(int(k,ik8),mask32))
    h = lowbias32(h + iand(int(i,ik8),mask32))
    hash_uint32 = lowbias32(h + iand(int(j,ik8),mask32))
  end function hash_uint32

  ! A uniform random number in [0,1) with 24 random bits, so that it is
  ! exact (and below one) also in single precision.
  pure real(rkx) function hash_uniform(key,k,i,j)
!$acc routine seq
    implicit none
    integer(ik8), intent(in) :: key
    integer(ik4), intent(in) :: k, i, j
    hash_uniform = real(real(ishft(hash_uint32(key,k,i,j),-8),rk8) * &
                        two_m24,rkx)
  end function hash_uniform

end module mod_hashrng

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

module mod_grid

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_memutil
  use mod_message
  use mod_stdio

  implicit none

  private

  real(rkx), public :: clatx, clonx

  real(rkx), public, pointer, contiguous, dimension(:,:) :: xlat, xlon, xmask, topo
  real(rkx), public, pointer, contiguous, dimension(:) :: sigx, zita
  real(rk4), public, pointer, contiguous, dimension(:,:) :: rxlat, rxlon
  real(rk4), public, pointer, contiguous, dimension(:) :: rsigx, ax, bx
  logical, public, pointer, contiguous, dimension(:,:,:) :: sgmask
  real(rk4), pointer, contiguous, dimension(:,:,:) :: xtrans_r4
  integer(ik4), pointer, contiguous, dimension(:,:,:) :: xtrans_i4
  real(rk8), pointer, contiguous, dimension(:,:,:) :: xtrans_r8

  interface mypack
    module procedure pack_integer
    module procedure pack_real4
    module procedure pack_real8
  end interface

  interface g2s
    module procedure g2s_i
    module procedure g2s_r4
    module procedure g2s_r8
  end interface g2s

  integer(ik4) :: js, je, is, ie
  public :: init_domain, mypack, setup_pack

  real(rkx), public :: h_missing_value = -9999.0_rkx

  contains

  subroutine init_domain
    implicit none
    call getmem(xlat,1,jxsg,1,iysg,'mod_read_domain:xlat')
    call getmem(xlon,1,jxsg,1,iysg,'mod_read_domain:xlon')
    call getmem(rxlat,1,jxsg,1,iysg,'mod_read_domain:rxlat')
    call getmem(rxlon,1,jxsg,1,iysg,'mod_read_domain:rxlon')
    call getmem(xmask,1,jxsg,1,iysg,'mod_read_domain:xmask')
    call getmem(topo,1,jxsg,1,iysg,'mod_read_domain:topo')
    call getmem(xtrans_i4,1,nnsg,1,jx,1,iy,'mod_read_domain:xtrans_i4')
    call getmem(xtrans_r4,1,nnsg,1,jx,1,iy,'mod_read_domain:xtrans_r4')
    call getmem(xtrans_r8,1,nnsg,1,jx,1,iy,'mod_read_domain:xtrans_r8')
    call getmem(sgmask,1,nnsg,1,jx,1,iy,'mod_read_domain:sgmask')
    call getmem(sigx,1,kzp1,'mod_read_domain:sigx')
    call getmem(zita,1,kzp1,'mod_read_domain:zita')
    call getmem(rsigx,1,kzp1,'mod_read_domain:rsigx')
    call getmem(ax,1,kzp1,'mod_read_domain:ax')
    call getmem(bx,1,kzp1,'mod_read_domain:bx')
  end subroutine init_domain

  subroutine g2s_i(m2,m3)
    implicit none
    integer(ik4), dimension(:,:), intent(in) :: m2
    integer(ik4), dimension(:,:,:), intent(inout) :: m3
    integer(ik4) :: n1, n2, j, i, ii, jj
    do i = 1, iy
      do j = 1, jx
        do n2 = 1, nsg
          ii = (i-1) * nsg + n2
          do n1 = 1, nsg
            jj = (j-1) * nsg + n1
            m3((n2-1)*nsg+n1,j,i) = m2(jj,ii)
          end do
        end do
      end do
    end do
  end subroutine g2s_i

  subroutine g2s_r4(m2,m3)
    implicit none
    real(rk4), dimension(:,:), intent(in) :: m2
    real(rk4), dimension(:,:,:), intent(inout) :: m3
    integer(ik4) :: n1, n2, j, i, ii, jj
    do i = 1, iy
      do j = 1, jx
        do n2 = 1, nsg
          ii = (i-1) * nsg + n2
          do n1 = 1, nsg
            jj = (j-1) * nsg + n1
            m3((n2-1)*nsg+n1,j,i) = m2(jj,ii)
          end do
        end do
      end do
    end do
  end subroutine g2s_r4

  subroutine g2s_r8(m2,m3)
    implicit none
    real(rk8), dimension(:,:), intent(in) :: m2
    real(rk8), dimension(:,:,:), intent(inout) :: m3
    integer(ik4) :: n1, n2, j, i, ii, jj
    do i = 1, iy
      do j = 1, jx
        do n2 = 1, nsg
          ii = (i-1) * nsg + n2
          do n1 = 1, nsg
            jj = (j-1) * nsg + n1
            m3((n2-1)*nsg+n1,j,i) = m2(jj,ii)
          end do
        end do
      end do
    end do
  end subroutine g2s_r8

  subroutine setup_pack(j1,j2,i1,i2)
    implicit none
    integer(ik4), intent(in) :: j1, j2, i1, i2
    integer(ik4) :: n1, n2, j, i, ii, jj
    js = j1
    je = j2
    is = i1
    ie = i2
    sgmask = .false.
    do i = is, ie
      do j = js, je
        do n2 = 1, nsg
          ii = (i-1) * nsg + n2
          do n1 = 1, nsg
            jj = (j-1) * nsg + n1
            sgmask((n2-1)*nsg+n1,j,i) = (xmask(jj,ii) > 0.5_rkx)
          end do
        end do
      end do
    end do
  end subroutine setup_pack

  subroutine pack_integer(matrix,vector)
    implicit none
    integer(ik4), dimension(:,:), intent(in) :: matrix
    integer(ik4), dimension(:), intent(inout) :: vector
    integer(ik4) :: i, j, n, ip
    call g2s(matrix,xtrans_i4)
    ip = 1
    do i = is, ie
      do j = js, je
        do n = 1, nnsg
          if ( sgmask(n,j,i) ) then
            vector(ip) = xtrans_i4(n,j,i)
            ip = ip + 1
          end if
        end do
      end do
    end do
  end subroutine pack_integer

  subroutine pack_real4(matrix,vector)
    implicit none
    real(rk4), dimension(:,:), intent(in) :: matrix
    real(rk4), dimension(:), intent(inout) :: vector
    integer(ik4) :: i, j, n, ip
    call g2s(matrix,xtrans_r4)
    ip = 1
    do i = is, ie
      do j = js, je
        do n = 1, nnsg
          if ( sgmask(n,j,i) ) then
            vector(ip) = xtrans_r4(n,j,i)
            ip = ip + 1
          end if
        end do
      end do
    end do
  end subroutine pack_real4

  subroutine pack_real8(matrix,vector)
    implicit none
    real(rk8), dimension(:,:), intent(in) :: matrix
    real(rk8), dimension(:), intent(inout) :: vector
    integer(ik4) :: i, j, n, ip
    call g2s(matrix,xtrans_r8)
    ip = 1
    do i = is, ie
      do j = js, je
        do n = 1, nnsg
          if ( sgmask(n,j,i) ) then
            vector(ip) = xtrans_r8(n,j,i)
            ip = ip + 1
          end if
        end do
      end do
    end do
  end subroutine pack_real8

end module mod_grid
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

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

module mod_grid

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_memutil
  use mod_message
  use mod_stdio

  implicit none

  private

  real(rk8) , public :: clatx , clonx

  real(rk8) , public , pointer , dimension(:,:) :: xlat , xlon , xmask , topo
  real(rk8) , public , pointer , dimension(:) :: sigx
  real(rk4) , public , pointer , dimension(:,:) :: rxlat , rxlon
  real(rk4) , public , pointer , dimension(:) :: rsigx
  logical , public , pointer , dimension(:,:,:) :: sgmask
  real(rk4) , public , pointer , dimension(:,:,:) :: xtrans_r4
  integer(ik4) , public , pointer , dimension(:,:,:) :: xtrans_i4
  real(rk8) , public , pointer , dimension(:,:,:) :: xtrans_r8

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

  integer(ik4) :: jgstart , jgstop
  public :: init_domain , mypack , setup_pack

  contains

  subroutine init_domain
    implicit none
    call getmem2d(xlat,1,jxsg,1,iysg,'mod_read_domain:xlat')
    call getmem2d(xlon,1,jxsg,1,iysg,'mod_read_domain:xlon')
    call getmem2d(rxlat,1,jxsg,1,iysg,'mod_read_domain:rxlat')
    call getmem2d(rxlon,1,jxsg,1,iysg,'mod_read_domain:rxlon')
    call getmem2d(xmask,1,jxsg,1,iysg,'mod_read_domain:xmask')
    call getmem2d(topo,1,jxsg,1,iysg,'mod_read_domain:topo')
    call getmem3d(xtrans_i4,1,nnsg,1,jx,1,iy,'mod_read_domain:xtrans_i4')
    call getmem3d(xtrans_r4,1,nnsg,1,jx,1,iy,'mod_read_domain:xtrans_r4')
    call getmem3d(xtrans_r8,1,nnsg,1,jx,1,iy,'mod_read_domain:xtrans_r8')
    call getmem3d(sgmask,1,nnsg,1,jx,1,iy,'mod_read_domain:sgmask')
    call getmem1d(sigx,1,kzp1,'mod_read_domain:sigx')
    call getmem1d(rsigx,1,kzp1,'mod_read_domain:rsigx')
  end subroutine init_domain

  subroutine g2s_i(m2,m3)
    implicit none
    integer(ik4) , dimension(:,:) :: m2
    integer(ik4) , dimension(:,:,:) :: m3
    integer(ik4) :: n1 , n2 , j , i , ii , jj
    do i = 1 , iy
      do j = 1 , jx
        do n2 = 1 , nsg
          ii = (i-1) * nsg + n2
          do n1 = 1 , nsg
            jj = (j-1) * nsg + n1
            m3((n2-1)*nsg+n1,j,i) = m2(jj,ii)
          end do
        end do
      end do
    end do
  end subroutine g2s_i

  subroutine g2s_r4(m2,m3)
    implicit none
    real(rk4) , dimension(:,:) :: m2
    real(rk4) , dimension(:,:,:) :: m3
    integer(ik4) :: n1 , n2 , j , i , ii , jj
    do i = 1 , iy
      do j = 1 , jx
        do n2 = 1 , nsg
          ii = (i-1) * nsg + n2
          do n1 = 1 , nsg
            jj = (j-1) * nsg + n1
            m3((n2-1)*nsg+n1,j,i) = m2(jj,ii)
          end do
        end do
      end do
    end do
  end subroutine g2s_r4

  subroutine g2s_r8(m2,m3)
    implicit none
    real(rk8) , dimension(:,:) :: m2
    real(rk8) , dimension(:,:,:) :: m3
    integer(ik4) :: n1 , n2 , j , i , ii , jj
    do i = 1 , iy
      do j = 1 , jx
        do n2 = 1 , nsg
          ii = (i-1) * nsg + n2
          do n1 = 1 , nsg
            jj = (j-1) * nsg + n1
            m3((n2-1)*nsg+n1,j,i) = m2(jj,ii)
          end do
        end do
      end do
    end do
  end subroutine g2s_r8

  subroutine setup_pack(j1,j2)
    implicit none
    integer(ik4) , intent(in) :: j1 , j2
    integer(ik4) :: n1 , n2 , j , i , ii , jj
    jgstart = j1
    jgstop = j2
    do i = 1 , iy
      do j = 1 , jx
        do n2 = 1 , nsg
          ii = (i-1) * nsg + n2
          do n1 = 1 , nsg
            jj = (j-1) * nsg + n1
            sgmask((n2-1)*nsg+n1,j,i) = (xmask(jj,ii) > 0.5D0)
          end do
        end do
      end do
    end do
  end subroutine setup_pack

  subroutine pack_integer(matrix,vector)
    implicit none
    integer(ik4) , dimension(:,:) :: matrix
    integer(ik4) , dimension(:) :: vector
    integer(ik4) :: i , j , n , ip
    call g2s(matrix,xtrans_i4)
    ip = 1
    do i = 2 , iy-2
      do j = jgstart , jgstop
        do n = 1 , nnsg
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
    real(rk4) , dimension(:,:) :: matrix
    real(rk4) , dimension(:) :: vector
    integer(ik4) :: i , j , n , ip
    call g2s(matrix,xtrans_r4)
    ip = 1
    do i = 2 , iy-2
      do j = jgstart , jgstop
        do n = 1 , nnsg
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
    real(rk8) , dimension(:,:) :: matrix
    real(rk8) , dimension(:) :: vector
    integer(ik4) :: i , j , n , ip
    call g2s(matrix,xtrans_r8)
    ip = 1
    do i = 2 , iy-2
      do j = jgstart , jgstop
        do n = 1 , nnsg
          if ( sgmask(n,j,i) ) then
            vector(ip) = xtrans_r8(n,j,i)
            ip = ip + 1
          end if
        end do
      end do
    end do
  end subroutine pack_real8

end module mod_grid

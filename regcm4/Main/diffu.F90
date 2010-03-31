!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
      subroutine diffu_u(ften,xkc,c203,j,ind)
!                                                                     c
!     this subroutine computes the diffusion term for decoupled       c
!     variable on constant sigma surface.                             c
!                                                                     c
!     ften    : tendency for variable ubd3d                           c
!                                                                     c
!     ubd3d   : coupled variable at time t-1                          c
!                                                                     c
!     xkc     : horizontal diffusion coefficient                      c
!                                                                     c
!     msfd    : map scale factor at the points where ubd3d is defined c
!                                                                     c
!     scr     : dummy array used for working space                    c
!                                                                     c
!     c203    : 1./(dx*dx), defined in 'param'                        c
!                                                                     c
!     j       : j'th slice of variable ubd3d                          c
!                                                                     c
!     iend    : = iym2 for cross-point variables                      c
!               = iym1  for dot-point   variables                      c
!                                                                     c
!     jend    : = jxm2 for cross-point variables                      c
!               = jlx  for dot-point   variables                      c
!                                                                     c
!     ind = 1 : ubd3d is already multiplied by map scale factor (msfd)c
!         = 0 : ubd3d is "not"   multiplied by msfd                   c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use mod_regcm_param
      use mod_slice
      use mod_main
      implicit none
!
! Dummy arguments
!
      real(8) :: c203
      integer :: ind , j
      real(8) , dimension(iy,kz) :: ften , xkc
      intent (in) c203 , ind , j , xkc
      intent (inout) ften
!
! Local variables
!
      integer :: i , k
!
!---------------------------------------------------------------------
!
#ifdef MPP1
      if ( (myid.eq.0 .and. j.eq.2) .or.                                &
         & (myid.eq.nproc-1 .and. j.eq.jendx) ) then
#else
      if (j.eq.2 .or. j.eq.jxm1) then
#endif
!
!......second-order scheme for east or west boundary:
        do k = 1 , kz
          do i = 2 , iym1
            if ( ind.eq.0 ) then
              ften(i,k) = ften(i,k) + xkc(i,k)                          &
                        & *c203*(ubd3d(i,k,j+1)+ubd3d(i,k,j-1)          &
                        & +ubd3d(i+1,k,j)+ubd3d(i-1,k,j)-4.*ubd3d(i,k,j)&
                        & )*pdotb(i,j)
            else
              ften(i,k) = ften(i,k) + xkc(i,k)                          &
                        & *c203*(ubd3d(i,k,j+1)/msfd(i,j+1)             &
                        & +ubd3d(i,k,j-1)/msfd(i,j-1)+ubd3d(i+1,k,j)    &
                        & /msfd(i+1,j)+ubd3d(i-1,k,j)/msfd(i-1,j)       &
                        & -4.*ubd3d(i,k,j)/msfd(i,j))*pdotb(i,j)
            end if
          end do
        end do
!
      else
!
!.....fourth-order scheme for interior:
        do k = 1 , kz
          do i = 3 , iym1 - 1
            if ( ind.eq.0 ) then
              ften(i,k) = ften(i,k) - xkc(i,k)                          &
                        & *c203*(ubd3d(i,k,j+2)+ubd3d(i,k,j-2)          &
                        & +ubd3d(i+2,k,j)+ubd3d(i-2,k,j)                &
                        & -4.*(ubd3d(i,k,j+1)+ubd3d(i,k,j-1)            &
                        & +ubd3d(i+1,k,j)+ubd3d(i-1,k,j))               &
                        & +12.*ubd3d(i,k,j))*pdotb(i,j)
            else
              ften(i,k) = ften(i,k) - xkc(i,k)                          &
                        & *c203*(ubd3d(i,k,j+2)/msfd(i,j+2)             &
                        & +ubd3d(i,k,j-2)/msfd(i,j-2)+ubd3d(i+2,k,j)    &
                        & /msfd(i+2,j)+ubd3d(i-2,k,j)/msfd(i-2,j)       &
                        & -4.*(ubd3d(i,k,j+1)/msfd(i,j+1)+ubd3d(i,k,j-1)&
                        & /msfd(i,j-1)+ubd3d(i+1,k,j)/msfd(i+1,j)       &
                        & +ubd3d(i-1,k,j)/msfd(i-1,j))+12.*ubd3d(i,k,j) &
                        & /msfd(i,j))*pdotb(i,j)
            end if
          end do
        end do
!......second-order scheme for north and south boundaries:
        do i = 2 , iym1 , iym1 - 2
          do k = 1 , kz
            if ( ind.eq.0 ) then
              ften(i,k) = ften(i,k) + xkc(i,k)                          &
                        & *c203*(ubd3d(i,k,j+1)+ubd3d(i,k,j-1)          &
                        & +ubd3d(i+1,k,j)+ubd3d(i-1,k,j)-4.*ubd3d(i,k,j)&
                        & )*pdotb(i,j)
            else
              ften(i,k) = ften(i,k) + xkc(i,k)                          &
                        & *c203*(ubd3d(i,k,j+1)/msfd(i,j+1)             &
                        & +ubd3d(i,k,j-1)/msfd(i,j-1)+ubd3d(i+1,k,j)    &
                        & /msfd(i+1,j)+ubd3d(i-1,k,j)/msfd(i-1,j)       &
                        & -4.*ubd3d(i,k,j)/msfd(i,j))*pdotb(i,j)
            end if
          end do
        end do
!
      end if
!
      end subroutine diffu_u
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine diffu_v(ften,xkc,c203,j,ind)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine computes the diffusion term for decoupled       c
!     variable on constant sigma surface.                             c
!                                                                     c
!     ften    : tendency for variable vbd3d                           c
!                                                                     c
!     vbd3d   : coupled variable at time t-1                          c
!                                                                     c
!     xkc     : horizontal diffusion coefficient                      c
!                                                                     c
!     msfd    : map scale factor at the points where vbd3d is defined c
!                                                                     c
!     scr     : dummy array used for working space                    c
!                                                                     c
!     c203    : 1./(dx*dx), defined in 'param'                        c
!                                                                     c
!     j       : j'th slice of variable vbd3d                          c
!                                                                     c
!     iend    : = iym2 for cross-point variables                      c
!               = iym1  for dot-point   variables                      c
!                                                                     c
!     jend    : = jxm2 for cross-point variables                      c
!               = jlx  for dot-point   variables                      c
!                                                                     c
!     ind = 1 : vbd3d is already multiplied by map scale factor (msfd)c
!         = 0 : vbd3d is "not"   multiplied by msfd                   c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use mod_regcm_param
      use mod_slice
      use mod_main
      implicit none
!
! Dummy arguments
!
      real(8) :: c203
      integer :: ind , j
      real(8) , dimension(iy,kz) :: ften , xkc
      intent (in) c203 , ind , j , xkc
      intent (inout) ften
!
! Local variables
!
      integer :: i , k
!
!---------------------------------------------------------------------
!
#ifdef MPP1
      if ( (myid.eq.0 .and. j.eq.2) .or.                                &
         & (myid.eq.nproc-1 .and. j.eq.jendx) ) then
#else
      if (j.eq.2 .or. j.eq.jxm1) then
#endif
!
!......second-order scheme for east or west boundary:
        do k = 1 , kz
          do i = 2 , iym1
            if ( ind.eq.0 ) then
              ften(i,k) = ften(i,k) + xkc(i,k)                          &
                        & *c203*(vbd3d(i,k,j+1)+vbd3d(i,k,j-1)          &
                        & +vbd3d(i+1,k,j)+vbd3d(i-1,k,j)-4.*vbd3d(i,k,j)&
                        & )*pdotb(i,j)
            else
              ften(i,k) = ften(i,k) + xkc(i,k)                          &
                        & *c203*(vbd3d(i,k,j+1)/msfd(i,j+1)             &
                        & +vbd3d(i,k,j-1)/msfd(i,j-1)+vbd3d(i+1,k,j)    &
                        & /msfd(i+1,j)+vbd3d(i-1,k,j)/msfd(i-1,j)       &
                        & -4.*vbd3d(i,k,j)/msfd(i,j))*pdotb(i,j)
            end if
          end do
        end do
!
      else
!
!.....fourth-order scheme for interior:
        do k = 1 , kz
          do i = 3 , iym1 - 1
            if ( ind.eq.0 ) then
              ften(i,k) = ften(i,k) - xkc(i,k)                          &
                        & *c203*(vbd3d(i,k,j+2)+vbd3d(i,k,j-2)          &
                        & +vbd3d(i+2,k,j)+vbd3d(i-2,k,j)                &
                        & -4.*(vbd3d(i,k,j+1)+vbd3d(i,k,j-1)            &
                        & +vbd3d(i+1,k,j)+vbd3d(i-1,k,j))               &
                        & +12.*vbd3d(i,k,j))*pdotb(i,j)
            else
              ften(i,k) = ften(i,k) - xkc(i,k)                          &
                        & *c203*(vbd3d(i,k,j+2)/msfd(i,j+2)             &
                        & +vbd3d(i,k,j-2)/msfd(i,j-2)+vbd3d(i+2,k,j)    &
                        & /msfd(i+2,j)+vbd3d(i-2,k,j)/msfd(i-2,j)       &
                        & -4.*(vbd3d(i,k,j+1)/msfd(i,j+1)+vbd3d(i,k,j-1)&
                        & /msfd(i,j-1)+vbd3d(i+1,k,j)/msfd(i+1,j)       &
                        & +vbd3d(i-1,k,j)/msfd(i-1,j))+12.*vbd3d(i,k,j) &
                        & /msfd(i,j))*pdotb(i,j)
            end if
          end do
        end do
!......second-order scheme for north and south boundaries:
        do i = 2 , iym1 , iym1 - 2
          do k = 1 , kz
            if ( ind.eq.0 ) then
              ften(i,k) = ften(i,k) + xkc(i,k)                          &
                        & *c203*(vbd3d(i,k,j+1)+vbd3d(i,k,j-1)          &
                        & +vbd3d(i+1,k,j)+vbd3d(i-1,k,j)-4.*vbd3d(i,k,j)&
                        & )*pdotb(i,j)
            else
              ften(i,k) = ften(i,k) + xkc(i,k)                          &
                        & *c203*(vbd3d(i,k,j+1)/msfd(i,j+1)             &
                        & +vbd3d(i,k,j-1)/msfd(i,j-1)+vbd3d(i+1,k,j)    &
                        & /msfd(i+1,j)+vbd3d(i-1,k,j)/msfd(i-1,j)       &
                        & -4.*vbd3d(i,k,j)/msfd(i,j))*pdotb(i,j)
            end if
          end do
        end do
!
      end if
!
      end subroutine diffu_v

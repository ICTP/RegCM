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
 
      module mod_diffusion
!
! Diffusion calculations
!
      use mod_dynparam
      use mod_runparams
      use mod_slice
      use mod_main

      private

      public :: diffu_u , diffu_v , diffut_t , diffutqv , diffutqc , &
                diffutch
!
      contains
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     These subroutines computes the diffusion term for decoupled     c
!     variable on constant sigma surface.                             c
!                                                                     c
!     ften    : tendency for variable                                 c
!                                                                     c
!     xkc     : horizontal diffusion coefficient                      c
!                                                                     c
!     j       : j'th slice of variable ubd3d                          c
!                                                                     c
!     ind = 1 : var is already multiplied by map scale factor         c
!         = 0 : var is "not"   multiplied by map scale factor         c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine diffu_u(ften,xkc,j,ind)
!
      implicit none
!
! Dummy arguments
!
      integer :: ind , j
      real(8) , dimension(iy,kz) :: ften , xkc
      intent (in) ind , j , xkc
      intent (inout) ften
!
! Local variables
!
      integer :: i , k
      integer :: jm1 , jm2 , jp1, jp2
!
      jm1 = j - 1
      jm2 = j - 2
      jp1 = j + 1
      jp2 = j + 2
#ifdef BAND
!---------------------------------------------------------------------
#if defined(BAND) && (!defined(MPP1))
      if(jm1.lt.1) jm1 = jm1 + jx
      if(jm2.lt.1) jm2 = jm2 + jx
      if(jp1.gt.jx) jp1 = jp1 -jx
      if(jp2.gt.jx) jp2 = jp2 -jx
#endif
!
!.....fourth-order scheme for interior:
      do k = 1 , kz
        do i = 3 , iym1 - 1
          if ( ind.eq.0 ) then
            ften(i,k) = ften(i,k) - xkc(i,k)                          &
                      & *c203*(ubd3d(i,k,jp2)+ubd3d(i,k,jm2)          &
                      & +ubd3d(i+2,k,j)+ubd3d(i-2,k,j)                &
                      & -4.*(ubd3d(i,k,jp1)+ubd3d(i,k,jm1)            &
                      & +ubd3d(i+1,k,j)+ubd3d(i-1,k,j))               &
                      & +12.*ubd3d(i,k,j))*pdotb(i,j)
          else
            ften(i,k) = ften(i,k) - xkc(i,k)                          &
                      & *c203*(ubd3d(i,k,jp2)/msfd(i,jp2)             &
                      & +ubd3d(i,k,jm2)/msfd(i,jm2)+ubd3d(i+2,k,j)    &
                      & /msfd(i+2,j)+ubd3d(i-2,k,j)/msfd(i-2,j)       &
                      & -4.*(ubd3d(i,k,jp1)/msfd(i,jp1)+ubd3d(i,k,jm1)&
                      & /msfd(i,jm1)+ubd3d(i+1,k,j)/msfd(i+1,j)       &
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
                      & *c203*(ubd3d(i,k,jp1)+ubd3d(i,k,jm1)          &
                      & +ubd3d(i+1,k,j)+ubd3d(i-1,k,j)-4.*ubd3d(i,k,j)&
                      & )*pdotb(i,j)
          else
            ften(i,k) = ften(i,k) + xkc(i,k)                          &
                      & *c203*(ubd3d(i,k,jp1)/msfd(i,jp1)             &
                      & +ubd3d(i,k,jm1)/msfd(i,jm1)+ubd3d(i+1,k,j)    &
                      & /msfd(i+1,j)+ubd3d(i-1,k,j)/msfd(i-1,j)       &
                      & -4.*ubd3d(i,k,j)/msfd(i,j))*pdotb(i,j)
          end if
        end do
      end do
!
#else
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
                        & *c203*(ubd3d(i,k,jp1)+ubd3d(i,k,jm1)          &
                        & +ubd3d(i+1,k,j)+ubd3d(i-1,k,j)-4.*ubd3d(i,k,j)&
                        & )*pdotb(i,j)
            else
              ften(i,k) = ften(i,k) + xkc(i,k)                          &
                        & *c203*(ubd3d(i,k,jp1)/msfd(i,jp1)             &
                        & +ubd3d(i,k,jm1)/msfd(i,jm1)+ubd3d(i+1,k,j)    &
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
                        & *c203*(ubd3d(i,k,jp2)+ubd3d(i,k,jm2)          &
                        & +ubd3d(i+2,k,j)+ubd3d(i-2,k,j)                &
                        & -4.*(ubd3d(i,k,jp1)+ubd3d(i,k,jm1)            &
                        & +ubd3d(i+1,k,j)+ubd3d(i-1,k,j))               &
                        & +12.*ubd3d(i,k,j))*pdotb(i,j)
            else
              ften(i,k) = ften(i,k) - xkc(i,k)                          &
                        & *c203*(ubd3d(i,k,jp2)/msfd(i,jp2)             &
                        & +ubd3d(i,k,jm2)/msfd(i,jm2)+ubd3d(i+2,k,j)    &
                        & /msfd(i+2,j)+ubd3d(i-2,k,j)/msfd(i-2,j)       &
                        & -4.*(ubd3d(i,k,jp1)/msfd(i,jp1)+ubd3d(i,k,jm1)&
                        & /msfd(i,jm1)+ubd3d(i+1,k,j)/msfd(i+1,j)       &
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
                        & *c203*(ubd3d(i,k,jp1)+ubd3d(i,k,jm1)          &
                        & +ubd3d(i+1,k,j)+ubd3d(i-1,k,j)-4.*ubd3d(i,k,j)&
                        & )*pdotb(i,j)
            else
              ften(i,k) = ften(i,k) + xkc(i,k)                          &
                        & *c203*(ubd3d(i,k,jp1)/msfd(i,jp1)             &
                        & +ubd3d(i,k,jm1)/msfd(i,jm1)+ubd3d(i+1,k,j)    &
                        & /msfd(i+1,j)+ubd3d(i-1,k,j)/msfd(i-1,j)       &
                        & -4.*ubd3d(i,k,j)/msfd(i,j))*pdotb(i,j)
            end if
          end do
        end do
!
      end if
!
#endif
      end subroutine diffu_u
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine diffu_v(ften,xkc,j,ind)
!
      implicit none
!
      integer :: ind , j
      real(8) , dimension(iy,kz) :: ften , xkc
      intent (in) ind , j , xkc
      intent (inout) ften
!
      integer :: i , k
      integer :: jm1 , jm2 , jp1, jp2
!
      jm1 = j - 1
      jm2 = j - 2
      jp1 = j + 1
      jp2 = j + 2
#ifdef BAND
!---------------------------------------------------------------------
#if defined(BAND) && (!defined(MPP1))
      if(jm1.lt.1) jm1 = jm1 + jx
      if(jm2.lt.1) jm2 = jm2 + jx
      if(jp1.gt.jx) jp1 = jp1 -jx
      if(jp2.gt.jx) jp2 = jp2 -jx
#endif
!
!.....fourth-order scheme for interior:
      do k = 1 , kz
        do i = 3 , iym1 - 1
          if ( ind.eq.0 ) then
            ften(i,k) = ften(i,k) - xkc(i,k)                          &
                      & *c203*(vbd3d(i,k,jp2)+vbd3d(i,k,jm2)          &
                      & +vbd3d(i+2,k,j)+vbd3d(i-2,k,j)                &
                      & -4.*(vbd3d(i,k,jp1)+vbd3d(i,k,jm1)            &
                      & +vbd3d(i+1,k,j)+vbd3d(i-1,k,j))               &
                      & +12.*vbd3d(i,k,j))*pdotb(i,j)
          else
            ften(i,k) = ften(i,k) - xkc(i,k)                          &
                      & *c203*(vbd3d(i,k,jp2)/msfd(i,jp2)             &
                      & +vbd3d(i,k,jm2)/msfd(i,jm2)+vbd3d(i+2,k,j)    &
                      & /msfd(i+2,j)+vbd3d(i-2,k,j)/msfd(i-2,j)       &
                      & -4.*(vbd3d(i,k,jp1)/msfd(i,jp1)+vbd3d(i,k,jm1)&
                      & /msfd(i,jm1)+vbd3d(i+1,k,j)/msfd(i+1,j)       &
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
                      & *c203*(vbd3d(i,k,jp1)+vbd3d(i,k,jm1)          &
                      & +vbd3d(i+1,k,j)+vbd3d(i-1,k,j)-4.*vbd3d(i,k,j)&
                      & )*pdotb(i,j)
          else
            ften(i,k) = ften(i,k) + xkc(i,k)                          &
                      & *c203*(vbd3d(i,k,jp1)/msfd(i,jp1)             &
                      & +vbd3d(i,k,jm1)/msfd(i,jm1)+vbd3d(i+1,k,j)    &
                      & /msfd(i+1,j)+vbd3d(i-1,k,j)/msfd(i-1,j)       &
                      & -4.*vbd3d(i,k,j)/msfd(i,j))*pdotb(i,j)
          end if
        end do
      end do
#else
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
                        & *c203*(vbd3d(i,k,jp1)+vbd3d(i,k,jm1)          &
                        & +vbd3d(i+1,k,j)+vbd3d(i-1,k,j)-4.*vbd3d(i,k,j)&
                        & )*pdotb(i,j)
            else
              ften(i,k) = ften(i,k) + xkc(i,k)                          &
                        & *c203*(vbd3d(i,k,jp1)/msfd(i,jp1)             &
                        & +vbd3d(i,k,jm1)/msfd(i,jm1)+vbd3d(i+1,k,j)    &
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
                        & *c203*(vbd3d(i,k,jp2)+vbd3d(i,k,jm2)          &
                        & +vbd3d(i+2,k,j)+vbd3d(i-2,k,j)                &
                        & -4.*(vbd3d(i,k,jp1)+vbd3d(i,k,jm1)            &
                        & +vbd3d(i+1,k,j)+vbd3d(i-1,k,j))               &
                        & +12.*vbd3d(i,k,j))*pdotb(i,j)
            else
              ften(i,k) = ften(i,k) - xkc(i,k)                          &
                        & *c203*(vbd3d(i,k,jp2)/msfd(i,jp2)             &
                        & +vbd3d(i,k,jm2)/msfd(i,jm2)+vbd3d(i+2,k,j)    &
                        & /msfd(i+2,j)+vbd3d(i-2,k,j)/msfd(i-2,j)       &
                        & -4.*(vbd3d(i,k,jp1)/msfd(i,jp1)+vbd3d(i,k,jm1)&
                        & /msfd(i,jm1)+vbd3d(i+1,k,j)/msfd(i+1,j)       &
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
                        & *c203*(vbd3d(i,k,jp1)+vbd3d(i,k,jm1)          &
                        & +vbd3d(i+1,k,j)+vbd3d(i-1,k,j)-4.*vbd3d(i,k,j)&
                        & )*pdotb(i,j)
            else
              ften(i,k) = ften(i,k) + xkc(i,k)                          &
                        & *c203*(vbd3d(i,k,jp1)/msfd(i,jp1)             &
                        & +vbd3d(i,k,jm1)/msfd(i,jm1)+vbd3d(i+1,k,j)    &
                        & /msfd(i+1,j)+vbd3d(i-1,k,j)/msfd(i-1,j)       &
                        & -4.*vbd3d(i,k,j)/msfd(i,j))*pdotb(i,j)
            end if
          end do
        end do
!
      end if
!
#endif
      end subroutine diffu_v
!
      subroutine diffut_t(ften,xkc,j)
!
      implicit none
!
      integer :: j
      real(8) , dimension(iy,kz) :: ften , xkc
      intent (in) j , xkc
      intent (inout) ften
!
      integer :: i , k
      integer :: jm1 , jm2 , jp1, jp2
!
!-----compute the diffusion term for t:
!
      jm1 = j - 1
      jm2 = j - 2
      jp1 = j + 1
      jp2 = j + 2
#ifdef BAND
!---------------------------------------------------------------------
#if defined(BAND) && (!defined(MPP1))
      if(jm1.lt.1) jm1 = jm1 + jx
      if(jm2.lt.1) jm2 = jm2 + jx
      if(jp1.gt.jx) jp1 = jp1 -jx
      if(jp2.gt.jx) jp2 = jp2 -jx
#endif
!
!......fourth-order scheme for interior:
      do k = 1 , kz
        do i = 3 , iym3
          ften(i,k) = ften(i,k) - xkc(i,k)                            &
                    & *c203*(tb3d(i,k,jp2)+tb3d(i,k,jm2)+tb3d(i+2,k,j)&
                    & +tb3d(i-2,k,j)                                  &
                    & -4.*(tb3d(i,k,jp1)+tb3d(i,k,jm1)+tb3d(i+1,k,j)  &
                    & +tb3d(i-1,k,j))+12.*tb3d(i,k,j))*psb(i,j)
        end do
      end do
!......second-order scheme for north and south boundaries:
      i = 2
      do k = 1 , kz
        ften(i,k) = ften(i,k) + xkc(i,k)                              &
                  & *c203*(tb3d(i,k,jp1)+tb3d(i,k,jm1)+tb3d(i+1,k,j)  &
                  & +tb3d(i-1,k,j)-4.*tb3d(i,k,j))*psb(i,j)
      end do
      i = iym2
      do k = 1 , kz
        ften(i,k) = ften(i,k) + xkc(i,k)                              &
                  & *c203*(tb3d(i,k,jp1)+tb3d(i,k,jm1)+tb3d(i+1,k,j)  &
                  & +tb3d(i-1,k,j)-4.*tb3d(i,k,j))*psb(i,j)
      end do
!
#else
!----------------------------------------------------------------------
!
#ifdef MPP1
      if ( (myid.eq.0 .and. j.eq.2) .or.                                &
         & (myid.eq.nproc-1 .and. j.eq.jendm) ) then
#else
      if ( j.eq.2 .or. j.eq.jxm2 ) then
#endif
!
!......second-order scheme for east or west boundary:
        do k = 1 , kz
          do i = 2 , iym2
            ften(i,k) = ften(i,k) + xkc(i,k)                            &
                      & *c203*(tb3d(i,k,jp1)+tb3d(i,k,jm1)+tb3d(i+1,k,j)&
                      & +tb3d(i-1,k,j)-4.*tb3d(i,k,j))*psb(i,j)
          end do
        end do
!
      else
!
!......fourth-order scheme for interior:
        do k = 1 , kz
          do i = 3 , iym3
            ften(i,k) = ften(i,k) - xkc(i,k)                            &
                      & *c203*(tb3d(i,k,jp2)+tb3d(i,k,jm2)+tb3d(i+2,k,j)&
                      & +tb3d(i-2,k,j)                                  &
                      & -4.*(tb3d(i,k,jp1)+tb3d(i,k,jm1)+tb3d(i+1,k,j)  &
                      & +tb3d(i-1,k,j))+12.*tb3d(i,k,j))*psb(i,j)
          end do
        end do
!......second-order scheme for north and south boundaries:
        i = 2
        do k = 1 , kz
          ften(i,k) = ften(i,k) + xkc(i,k)                              &
                    & *c203*(tb3d(i,k,jp1)+tb3d(i,k,jm1)+tb3d(i+1,k,j)  &
                    & +tb3d(i-1,k,j)-4.*tb3d(i,k,j))*psb(i,j)
        end do
        i = iym2
        do k = 1 , kz
          ften(i,k) = ften(i,k) + xkc(i,k)                              &
                    & *c203*(tb3d(i,k,jp1)+tb3d(i,k,jm1)+tb3d(i+1,k,j)  &
                    & +tb3d(i-1,k,j)-4.*tb3d(i,k,j))*psb(i,j)
        end do
!
      end if
!
#endif
      end subroutine diffut_t
!
      subroutine diffutqv(ften,xkc,j)
!
      implicit none
!
      integer :: j
      real(8) , dimension(iy,kz) :: ften , xkc
      intent (in) j , xkc
      intent (inout) ften
!
      integer :: i , k
      integer :: jm1 , jm2 , jp1, jp2
!
!-----compute the diffusion term for qv:
!
      jm1 = j - 1
      jm2 = j - 2
      jp1 = j + 1
      jp2 = j + 2
#ifdef BAND
!---------------------------------------------------------------------
#if defined(BAND) && (!defined(MPP1))
      if(jm1.lt.1) jm1 = jm1 + jx
      if(jm2.lt.1) jm2 = jm2 + jx
      if(jp1.gt.jx) jp1 = jp1 -jx
      if(jp2.gt.jx) jp2 = jp2 -jx
#endif
!
!......fourth-order scheme for interior:
      do k = 1 , kz
        do i = 3 , iym3
          ften(i,k) = ften(i,k) - xkc(i,k)                            &
                    & *c203*(qvb3d(i,k,jp2)+qvb3d(i,k,jm2)            &
                    & +qvb3d(i+2,k,j)+qvb3d(i-2,k,j)                  &
                    & -4.*(qvb3d(i,k,jp1)+qvb3d(i,k,jm1)              &
                    & +qvb3d(i+1,k,j)+qvb3d(i-1,k,j))+12.*qvb3d(i,k,j)&
                    & )*psb(i,j)
        end do
      end do
!......second-order scheme for north and south boundaries:
      i = 2
      do k = 1 , kz
        ften(i,k) = ften(i,k) + xkc(i,k)                              &
                  & *c203*(qvb3d(i,k,jp1)+qvb3d(i,k,jm1)              &
                  & +qvb3d(i+1,k,j)+qvb3d(i-1,k,j)-4.*qvb3d(i,k,j))   &
                  & *psb(i,j)
      end do
      i = iym2
      do k = 1 , kz
        ften(i,k) = ften(i,k) + xkc(i,k)                              &
                  & *c203*(qvb3d(i,k,jp1)+qvb3d(i,k,jm1)              &
                  & +qvb3d(i+1,k,j)+qvb3d(i-1,k,j)-4.*qvb3d(i,k,j))   &
                  & *psb(i,j)
      end do
#else
!----------------------------------------------------------------------
!
#ifdef MPP1
      if ( (myid.eq.0 .and. j.eq.2) .or.                                &
         & (myid.eq.nproc-1 .and. j.eq.jendm) ) then
#else
      if ( j.eq.2 .or. j.eq.jxm2 ) then
#endif
!
!......second-order scheme for east or west boundary:
        do k = 1 , kz
          do i = 2 , iym2
            ften(i,k) = ften(i,k) + xkc(i,k)                            &
                      & *c203*(qvb3d(i,k,jp1)+qvb3d(i,k,jm1)            &
                      & +qvb3d(i+1,k,j)+qvb3d(i-1,k,j)-4.*qvb3d(i,k,j)) &
                      & *psb(i,j)
          end do
        end do
!
      else
!
!......fourth-order scheme for interior:
        do k = 1 , kz
          do i = 3 , iym3
            ften(i,k) = ften(i,k) - xkc(i,k)                            &
                      & *c203*(qvb3d(i,k,jp2)+qvb3d(i,k,jm2)            &
                      & +qvb3d(i+2,k,j)+qvb3d(i-2,k,j)                  &
                      & -4.*(qvb3d(i,k,jp1)+qvb3d(i,k,jm1)              &
                      & +qvb3d(i+1,k,j)+qvb3d(i-1,k,j))+12.*qvb3d(i,k,j)&
                      & )*psb(i,j)
          end do
        end do
!......second-order scheme for north and south boundaries:
        i = 2
        do k = 1 , kz
          ften(i,k) = ften(i,k) + xkc(i,k)                              &
                    & *c203*(qvb3d(i,k,jp1)+qvb3d(i,k,jm1)              &
                    & +qvb3d(i+1,k,j)+qvb3d(i-1,k,j)-4.*qvb3d(i,k,j))   &
                    & *psb(i,j)
        end do
        i = iym2
        do k = 1 , kz
          ften(i,k) = ften(i,k) + xkc(i,k)                              &
                    & *c203*(qvb3d(i,k,jp1)+qvb3d(i,k,jm1)              &
                    & +qvb3d(i+1,k,j)+qvb3d(i-1,k,j)-4.*qvb3d(i,k,j))   &
                    & *psb(i,j)
        end do
!
      end if
!
#endif
      end subroutine diffutqv
!
      subroutine diffutqc(ften,xkc,j)
!
      implicit none
!
      integer :: j
      real(8) , dimension(iy,kz) :: ften , xkc
      intent (in) j , xkc
      intent (inout) ften
!
      integer :: i , k
      integer :: jm1 , jm2 , jp1, jp2
!
!-----compute the diffusion term for qc:
!
      jm1 = j - 1
      jm2 = j - 2
      jp1 = j + 1
      jp2 = j + 2
#ifdef BAND
!---------------------------------------------------------------------
#if defined(BAND) && (!defined(MPP1))
      if(jm1.lt.1) jm1 = jm1 + jx
      if(jm2.lt.1) jm2 = jm2 + jx
      if(jp1.gt.jx) jp1 = jp1 -jx
      if(jp2.gt.jx) jp2 = jp2 -jx
#endif
!
!......fourth-order scheme for interior:
      do k = 1 , kz
        do i = 3 , iym3
          ften(i,k) = ften(i,k) - xkc(i,k)                            &
                    & *c203*(qcb3d(i,k,jp2)+qcb3d(i,k,jm2)            &
                    & +qcb3d(i+2,k,j)+qcb3d(i-2,k,j)                  &
                    & -4.*(qcb3d(i,k,jp1)+qcb3d(i,k,jm1)              &
                    & +qcb3d(i+1,k,j)+qcb3d(i-1,k,j))+12.*qcb3d(i,k,j)&
                    & )*psb(i,j)
        end do
      end do
!......second-order scheme for north and south boundaries:
      i = 2
      do k = 1 , kz
        ften(i,k) = ften(i,k) + xkc(i,k)                              &
                  & *c203*(qcb3d(i,k,jp1)+qcb3d(i,k,jm1)              &
                  & +qcb3d(i+1,k,j)+qcb3d(i-1,k,j)-4.*qcb3d(i,k,j))   &
                  & *psb(i,j)
      end do
      i = iym2
      do k = 1 , kz
        ften(i,k) = ften(i,k) + xkc(i,k)                              &
                  & *c203*(qcb3d(i,k,jp1)+qcb3d(i,k,jm1)              &
                  & +qcb3d(i+1,k,j)+qcb3d(i-1,k,j)-4.*qcb3d(i,k,j))   &
                  & *psb(i,j)
      end do
#else
!----------------------------------------------------------------------
!
#ifdef MPP1
      if ( (myid.eq.0 .and. j.eq.2) .or.                                &
         & (myid.eq.nproc-1 .and. j.eq.jendm) ) then
#else
      if ( j.eq.2 .or. j.eq.jxm2 ) then
#endif
!
!......second-order scheme for east or west boundary:
        do k = 1 , kz
          do i = 2 , iym2
            ften(i,k) = ften(i,k) + xkc(i,k)                            &
                      & *c203*(qcb3d(i,k,jp1)+qcb3d(i,k,jm1)            &
                      & +qcb3d(i+1,k,j)+qcb3d(i-1,k,j)-4.*qcb3d(i,k,j)) &
                      & *psb(i,j)
          end do
        end do
!
      else
!
!......fourth-order scheme for interior:
        do k = 1 , kz
          do i = 3 , iym3
            ften(i,k) = ften(i,k) - xkc(i,k)                            &
                      & *c203*(qcb3d(i,k,jp2)+qcb3d(i,k,jm2)            &
                      & +qcb3d(i+2,k,j)+qcb3d(i-2,k,j)                  &
                      & -4.*(qcb3d(i,k,jp1)+qcb3d(i,k,jm1)              &
                      & +qcb3d(i+1,k,j)+qcb3d(i-1,k,j))+12.*qcb3d(i,k,j)&
                      & )*psb(i,j)
          end do
        end do
!......second-order scheme for north and south boundaries:
        i = 2
        do k = 1 , kz
          ften(i,k) = ften(i,k) + xkc(i,k)                              &
                    & *c203*(qcb3d(i,k,jp1)+qcb3d(i,k,jm1)              &
                    & +qcb3d(i+1,k,j)+qcb3d(i-1,k,j)-4.*qcb3d(i,k,j))   &
                    & *psb(i,j)
        end do
        i = iym2
        do k = 1 , kz
          ften(i,k) = ften(i,k) + xkc(i,k)                              &
                    & *c203*(qcb3d(i,k,jp1)+qcb3d(i,k,jm1)              &
                    & +qcb3d(i+1,k,j)+qcb3d(i-1,k,j)-4.*qcb3d(i,k,j))   &
                    & *psb(i,j)
        end do
!
      end if
!
#endif
      end subroutine diffutqc
!
      subroutine diffutch(ften,xkc,n,j)
!
      implicit none
!
      integer :: j , n
      real(8) , dimension(iy,kz) :: ften , xkc
      intent (in) j , n , xkc
      intent (inout) ften
!
      integer :: i , k
      integer :: jm1 , jm2 , jp1, jp2
!
!-----compute the diffusion term for chi:
!
      jm1 = j - 1
      jm2 = j - 2
      jp1 = j + 1
      jp2 = j + 2
#ifdef BAND
!---------------------------------------------------------------------
#if defined(BAND) && (!defined(MPP1))
      if(jm1.lt.1) jm1 = jm1 + jx
      if(jm2.lt.1) jm2 = jm2 + jx
      if(jp1.gt.jx) jp1 = jp1 -jx
      if(jp2.gt.jx) jp2 = jp2 -jx
#endif
!
!......fourth-order scheme for interior:
      do k = 1 , kz
        do i = 3 , iym3
          ften(i,k) = ften(i,k) - xkc(i,k)                            &
                    & *c203*(chib3d(i,k,jp2,n)+chib3d(i,k,jm2,n)      &
                    & +chib3d(i+2,k,j,n)+chib3d(i-2,k,j,n)            &
                    & -4.*(chib3d(i,k,jp1,n)+chib3d(i,k,jm1,n)        &
                    & +chib3d(i+1,k,j,n)+chib3d(i-1,k,j,n))           &
                    & +12.*chib3d(i,k,j,n))*psb(i,j)
        end do
      end do
!......second-order scheme for north and south boundaries:
      i = 2
      do k = 1 , kz
        ften(i,k) = ften(i,k) + xkc(i,k)                              &
                  & *c203*(chib3d(i,k,jp1,n)+chib3d(i,k,jm1,n)        &
                  & +chib3d(i+1,k,j,n)+chib3d(i-1,k,j,n)              &
                  & -4.*chib3d(i,k,j,n))*psb(i,j)
      end do
      i = iym2
      do k = 1 , kz
        ften(i,k) = ften(i,k) + xkc(i,k)                              &
                  & *c203*(chib3d(i,k,jp1,n)+chib3d(i,k,jm1,n)        &
                  & +chib3d(i+1,k,j,n)+chib3d(i-1,k,j,n)              &
                  & -4.*chib3d(i,k,j,n))*psb(i,j)
      end do
#else
!----------------------------------------------------------------------
!
#ifdef MPP1
      if ( (myid.eq.0 .and. j.eq.2) .or.                                &
         & (myid.eq.nproc-1 .and. j.eq.jendm) ) then
#else
      if ( j.eq.2 .or. j.eq.jxm2 ) then
#endif
!
!......second-order scheme for east or west boundary:
        do k = 1 , kz
          do i = 2 , iym2
            ften(i,k) = ften(i,k) + xkc(i,k)                            &
                      & *c203*(chib3d(i,k,jp1,n)+chib3d(i,k,jm1,n)      &
                      & +chib3d(i+1,k,j,n)+chib3d(i-1,k,j,n)            &
                      & -4.*chib3d(i,k,j,n))*psb(i,j)
          end do
        end do
!
      else
!
!......fourth-order scheme for interior:
        do k = 1 , kz
          do i = 3 , iym3
            ften(i,k) = ften(i,k) - xkc(i,k)                            &
                      & *c203*(chib3d(i,k,jp2,n)+chib3d(i,k,jm2,n)      &
                      & +chib3d(i+2,k,j,n)+chib3d(i-2,k,j,n)            &
                      & -4.*(chib3d(i,k,jp1,n)+chib3d(i,k,jm1,n)        &
                      & +chib3d(i+1,k,j,n)+chib3d(i-1,k,j,n))           &
                      & +12.*chib3d(i,k,j,n))*psb(i,j)
          end do
        end do
!......second-order scheme for north and south boundaries:
        i = 2
        do k = 1 , kz
          ften(i,k) = ften(i,k) + xkc(i,k)                              &
                    & *c203*(chib3d(i,k,jp1,n)+chib3d(i,k,jm1,n)        &
                    & +chib3d(i+1,k,j,n)+chib3d(i-1,k,j,n)              &
                    & -4.*chib3d(i,k,j,n))*psb(i,j)
        end do
        i = iym2
        do k = 1 , kz
          ften(i,k) = ften(i,k) + xkc(i,k)                              &
                    & *c203*(chib3d(i,k,jp1,n)+chib3d(i,k,jm1,n)        &
                    & +chib3d(i+1,k,j,n)+chib3d(i-1,k,j,n)              &
                    & -4.*chib3d(i,k,j,n))*psb(i,j)
        end do
!
      end if
!
#endif
      end subroutine diffutch
!
      end module mod_diffusion

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
      use mod_constants
      use mod_dynparam
      use mod_runparams
      use mod_service 
      private

      public :: diffu_d , diffu_x
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
!     j       : j'th slice of variable bd3d on dot points             c
!                                                                     c
!     ind = 1 : var is already multiplied by map scale factor         c
!         = 0 : var is "not"   multiplied by map scale factor         c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
      subroutine diffu_d(ften,bd3d,press,mapf,xkc,j,ind)
!
      implicit none
!
! Dummy arguments
!
      integer :: ind , j
      real(8) , dimension(iy,kz) :: ften , xkc
#ifdef MPP1
      real(8) , dimension(iy,kz,-1:jxp+2) , intent(in) :: bd3d
      real(8) , dimension(iy,-1:jxp+2) , intent(in) :: press
      real(8) , dimension(iy,-1:jxp+2) , intent(in) :: mapf
#else
      real(8) , dimension(iy,kz,jx) , intent(in) :: bd3d
      real(8) , dimension(iy,jx) , intent(in) :: press
      real(8) , dimension(iy,jx) , intent(in) :: mapf
#endif
      intent (in) ind , j , xkc
      intent (inout) ften
!
! Local variables
!
      integer :: i , k
      integer :: jm1 , jm2 , jp1, jp2
!
      character (len=50) :: subroutine_name='diffu_d'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
      jm1 = j - 1
      jm2 = j - 2
      jp1 = j + 1
      jp2 = j + 2
#ifdef BAND
!---------------------------------------------------------------------
#if defined(BAND) && (!defined(MPP1))
      if(jm1 < 1) jm1 = jm1 + jx
      if(jm2 < 1) jm2 = jm2 + jx
      if(jp1 > jx) jp1 = jp1 -jx
      if(jp2 > jx) jp2 = jp2 -jx
#endif
!
!.....fourth-order scheme for interior:
      do k = 1 , kz
        do i = 3 , iym1 - 1
          if ( ind == 0 ) then
            ften(i,k) = ften(i,k) - xkc(i,k)                  &
                      & *c203*(bd3d(i,k,jp2)+bd3d(i,k,jm2)    &
                      & +bd3d(i+2,k,j)+bd3d(i-2,k,j)          &
                      & -d_four*(bd3d(i,k,jp1)+bd3d(i,k,jm1)   &
                      & +bd3d(i+1,k,j)+bd3d(i-1,k,j))         &
                      & +d_twelve*bd3d(i,k,j))*press(i,j)
          else
            ften(i,k) = ften(i,k) - xkc(i,k)                  &
                      & *c203*(bd3d(i,k,jp2)/mapf(i,jp2)      &
                      & +bd3d(i,k,jm2)/mapf(i,jm2)+           &
                      & bd3d(i+2,k,j)/mapf(i+2,j)+            &
                      & bd3d(i-2,k,j)/mapf(i-2,j)             &
                      & -d_four*(bd3d(i,k,jp1)/mapf(i,jp1)+    &
                      & bd3d(i,k,jm1)/mapf(i,jm1)+            &
                      & bd3d(i+1,k,j)/mapf(i+1,j)             &
                      & +bd3d(i-1,k,j)/mapf(i-1,j))+          &
                      & d_twelve*bd3d(i,k,j)/mapf(i,j))*        &
                      & press(i,j)
          end if
        end do
      end do
!......second-order scheme for north and south boundaries:
      do i = 2 , iym1 , iym1 - 2
        do k = 1 , kz
          if ( ind == 0 ) then
            ften(i,k) = ften(i,k) + xkc(i,k)                  &
                      & *c203*(bd3d(i,k,jp1)+bd3d(i,k,jm1)    &
                      & +bd3d(i+1,k,j)+bd3d(i-1,k,j)          &
                      & -d_four*bd3d(i,k,j))*press(i,j)
          else
            ften(i,k) = ften(i,k) + xkc(i,k)                  &
                      & *c203*(bd3d(i,k,jp1)/mapf(i,jp1)      &
                      & +bd3d(i,k,jm1)/mapf(i,jm1)+           &
                      &  bd3d(i+1,k,j)/mapf(i+1,j)+           &
                      &  bd3d(i-1,k,j)/mapf(i-1,j)            &
                      & -d_four*bd3d(i,k,j)/mapf(i,j))*        &
                      & press(i,j)
          end if
        end do
      end do
!
#else
!---------------------------------------------------------------------
!
#ifdef MPP1
      if ( (myid == 0 .and. j == 2) .or.                        &
         & (myid == nproc-1 .and. j == jendx) ) then
#else
      if (j == 2 .or. j == jxm1) then
#endif
!
!......second-order scheme for east or west boundary:
        do k = 1 , kz
          do i = 2 , iym1
            if ( ind == 0 ) then
              ften(i,k) = ften(i,k) + xkc(i,k)                  &
                        & *c203*(bd3d(i,k,jp1)+bd3d(i,k,jm1)    &
                        & +bd3d(i+1,k,j)+bd3d(i-1,k,j)          &
                        & -d_four*bd3d(i,k,j))*press(i,j)
            else
              ften(i,k) = ften(i,k) + xkc(i,k)                  &
                        & *c203*(bd3d(i,k,jp1)/mapf(i,jp1)      &
                        & +bd3d(i,k,jm1)/mapf(i,jm1)+           &
                        &  bd3d(i+1,k,j)/mapf(i+1,j)+           &
                        &  bd3d(i-1,k,j)/mapf(i-1,j)            &
                        & -d_four*bd3d(i,k,j)/mapf(i,j))*        &
                        & press(i,j)
            end if
          end do
        end do
!
      else
!
!.....fourth-order scheme for interior:
        do k = 1 , kz
          do i = 3 , iym1 - 1
            if ( ind == 0 ) then
              ften(i,k) = ften(i,k) - xkc(i,k)                  &
                        & *c203*(bd3d(i,k,jp2)+bd3d(i,k,jm2)    &
                        & +bd3d(i+2,k,j)+bd3d(i-2,k,j)          &
                        & -d_four*(bd3d(i,k,jp1)+bd3d(i,k,jm1)   &
                        & +bd3d(i+1,k,j)+bd3d(i-1,k,j))         &
                        & +d_twelve*bd3d(i,k,j))*press(i,j)
            else
              ften(i,k) = ften(i,k) - xkc(i,k)                  &
                        & *c203*(bd3d(i,k,jp2)/mapf(i,jp2)      &
                        & +bd3d(i,k,jm2)/mapf(i,jm2)+           &
                        &  bd3d(i+2,k,j)/mapf(i+2,j)+           &
                        &  bd3d(i-2,k,j)/mapf(i-2,j)            &
                        & -d_four*(bd3d(i,k,jp1)/mapf(i,jp1)+    &
                        &      bd3d(i,k,jm1)/mapf(i,jm1)+       &
                        &      bd3d(i+1,k,j)/mapf(i+1,j)+       &
                        &      bd3d(i-1,k,j)/mapf(i-1,j))+      &
                        & d_twelve*bd3d(i,k,j)/mapf(i,j))*        &
                        & press(i,j)
            end if
          end do
        end do
!......second-order scheme for north and south boundaries:
        do i = 2 , iym1 , iym1 - 2
          do k = 1 , kz
            if ( ind == 0 ) then
              ften(i,k) = ften(i,k) + xkc(i,k)                  &
                        & *c203*(bd3d(i,k,jp1)+bd3d(i,k,jm1)    &
                        & +bd3d(i+1,k,j)+bd3d(i-1,k,j)          &
                        & -d_four*bd3d(i,k,j))*press(i,j)
            else
              ften(i,k) = ften(i,k) + xkc(i,k)                  &
                        & *c203*(bd3d(i,k,jp1)/mapf(i,jp1)      &
                        & +bd3d(i,k,jm1)/mapf(i,jm1)+           &
                        &  bd3d(i+1,k,j)/mapf(i+1,j)+           &
                        &  bd3d(i-1,k,j)/mapf(i-1,j)            &
                        & -d_four*bd3d(i,k,j)/mapf(i,j))*        &
                        & press(i,j)
            end if
          end do
        end do
!
      end if
!
#endif
      call time_end(subroutine_name,idindx) 
      end subroutine diffu_d
!
      subroutine diffu_x(ften,bc3d,press,xkc,j)
!
      implicit none
!
      integer :: j
      real(8) , dimension(iy,kz) :: ften , xkc
#ifdef MPP1
      real(8) , dimension(iy,kz,-1:jxp+2) , intent(in) :: bc3d
      real(8) , dimension(iy,-1:jxp+2) , intent(in) :: press
#else
      real(8) , dimension(iy,kz,jx) , intent(in) :: bc3d
      real(8) , dimension(iy,jx) , intent(in) :: press
#endif
      intent (in) j , xkc
      intent (inout) ften
!
      integer :: i , k
      integer :: jm1 , jm2 , jp1, jp2
!
      character (len=50) :: subroutine_name='diffu_x'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)
!-----compute the diffusion term for t:
!
      jm1 = j - 1
      jm2 = j - 2
      jp1 = j + 1
      jp2 = j + 2
#ifdef BAND
!---------------------------------------------------------------------
#if defined(BAND) && (!defined(MPP1))
      if(jm1 < 1) jm1 = jm1 + jx
      if(jm2 < 1) jm2 = jm2 + jx
      if(jp1 > jx) jp1 = jp1 -jx
      if(jp2 > jx) jp2 = jp2 -jx
#endif
!
!......fourth-order scheme for interior:
      do k = 1 , kz
        do i = 3 , iym3
          ften(i,k) = ften(i,k) - xkc(i,k) *                &
                    & c203*(bc3d(i,k,jp2)+bc3d(i,k,jm2)+    &
                    &       bc3d(i+2,k,j)+bc3d(i-2,k,j)     &
                    &  -d_four*(bc3d(i,k,jp1)+bc3d(i,k,jm1)+ &
                    &       bc3d(i+1,k,j)+bc3d(i-1,k,j))+   &
                    &   d_twelve*bc3d(i,k,j))*press(i,j)
        end do
      end do
!......second-order scheme for north and south boundaries:
      i = 2
      do k = 1 , kz
        ften(i,k) = ften(i,k) + xkc(i,k) *             &
                  & c203*(bc3d(i,k,jp1)+bc3d(i,k,jm1)+ &
                  &       bc3d(i+1,k,j)+bc3d(i-1,k,j)  &
                  &   -d_four*bc3d(i,k,j))*press(i,j)
      end do
      i = iym2
      do k = 1 , kz
        ften(i,k) = ften(i,k) + xkc(i,k) *             &
                  & c203*(bc3d(i,k,jp1)+bc3d(i,k,jm1)+ &
                  &       bc3d(i+1,k,j)+bc3d(i-1,k,j)  &
                  &   -d_four*bc3d(i,k,j))*press(i,j)
      end do
!
#else
!----------------------------------------------------------------------
!
#ifdef MPP1
      if ( (myid == 0 .and. j == 2) .or.                                &
         & (myid == nproc-1 .and. j == jendm) ) then
#else
      if ( j == 2 .or. j == jxm2 ) then
#endif
!
!......second-order scheme for east or west boundary:
        do k = 1 , kz
          do i = 2 , iym2
            ften(i,k) = ften(i,k) + xkc(i,k) *             &
                      & c203*(bc3d(i,k,jp1)+bc3d(i,k,jm1)+ &
                      &       bc3d(i+1,k,j)+bc3d(i-1,k,j)  &
                      &   -d_four*bc3d(i,k,j))*press(i,j)
          end do
        end do
!
      else
!
!......fourth-order scheme for interior:
        do k = 1 , kz
          do i = 3 , iym3
            ften(i,k) = ften(i,k) - xkc(i,k) *                &
                      & c203*(bc3d(i,k,jp2)+bc3d(i,k,jm2)+    &
                      &       bc3d(i+2,k,j)+bc3d(i-2,k,j)     &
                      &  -d_four*(bc3d(i,k,jp1)+bc3d(i,k,jm1)+ &
                      &       bc3d(i+1,k,j)+bc3d(i-1,k,j))    &
                      &  +d_twelve*bc3d(i,k,j))*press(i,j)
          end do
        end do
!......second-order scheme for north and south boundaries:
        i = 2
        do k = 1 , kz
          ften(i,k) = ften(i,k) + xkc(i,k) *             &
                    & c203*(bc3d(i,k,jp1)+bc3d(i,k,jm1)+ &
                    &       bc3d(i+1,k,j)+bc3d(i-1,k,j)  &
                    &   -d_four*bc3d(i,k,j))*press(i,j)
        end do
        i = iym2
        do k = 1 , kz
          ften(i,k) = ften(i,k) + xkc(i,k) *             &
                    & c203*(bc3d(i,k,jp1)+bc3d(i,k,jm1)+ &
                    &       bc3d(i+1,k,j)+bc3d(i-1,k,j)  &
                    &   -d_four*bc3d(i,k,j))*press(i,j)
        end do
!
      end if
!
#endif
      call time_end(subroutine_name,idindx)
      end subroutine diffu_x
!
      end module mod_diffusion

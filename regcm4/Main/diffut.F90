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
 
      subroutine diffut_t(ften,xkc,c203,j)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine computes the diffusion term for cross-point     c
!     variable 'f' on a constant sigma surface.                       c
!                                                                     c
!     ---fourth-order diffusion scheme is applied for the interior    c
!        and second-order scheme is applied for the boundary.         c
!                                                                     c
!     ften   : tendency for variable f                                c
!                                                                     c
!     tb3d   : variable f at t-1 time step                            c
!                                                                     c
!     psb    : p* at t-1 time step                                    c
!                                                                     c
!     xkc    : horizontal diffusion coefficient                       c
!                                                                     c
!     c203   : = c203 defined in 'param'                              c
!     c203   : = (100.-ptop)/(dx*dx), defined in 'param'              c
!                                                                     c
!     j      : j'th slice of variable tb3d                            c
!                                                                     c
!     iend   : = iym2 for cross-point variables                       c
!              = iym1  for dot-point   variables                       c
!                                                                     c
!     jend   : = jxm2 for cross-point variables                       c
!              = jlx  for dot-point   variables                       c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      use mod_dynparam
      use mod_slice
      use mod_main
      implicit none
!
! Dummy arguments
!
      real(8) :: c203
      integer :: j
      real(8) , dimension(iy,kz) :: ften , xkc
      intent (in) c203 , j , xkc
      intent (inout) ften
!
! Local variables
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
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine diffutqv(ften,xkc,c203,j)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine computes the diffusion term for cross-point     c
!     variable 'f' on a constant sigma surface.                       c
!                                                                     c
!     ---fourth-order diffusion scheme is applied for the interior    c
!        and second-order scheme is applied for the boundary.         c
!                                                                     c
!     ften   : tendency for variable f                                c
!                                                                     c
!     qvb3d  : variable f at t-1 time step                            c
!                                                                     c
!     psb    : p* at t-1 time step                                    c
!                                                                     c
!     xkc    : horizontal diffusion coefficient                       c
!                                                                     c
!     c203   : = c203 defined in 'param'                              c
!     c203   : = (100.-ptop)/(dx*dx), defined in 'param'              c
!                                                                     c
!     j      : j'th slice of variable qvb3d                           c
!                                                                     c
!     iend   : = iym2 for cross-point variables                       c
!              = iym1  for dot-point   variables                       c
!                                                                     c
!     jend   : = jxm2 for cross-point variables                       c
!              = jlx  for dot-point   variables                       c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use mod_dynparam
      use mod_slice
      use mod_main
      implicit none
!
! Dummy arguments
!
      real(8) :: c203
      integer :: j
      real(8) , dimension(iy,kz) :: ften , xkc
      intent (in) c203 , j , xkc
      intent (inout) ften
!
! Local variables
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
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine diffutqc(ften,xkc,c203,j)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine computes the diffusion term for cross-point     c
!     variable 'f' on a constant sigma surface.                       c
!                                                                     c
!     ---fourth-order diffusion scheme is applied for the interior    c
!        and second-order scheme is applied for the boundary.         c
!                                                                     c
!     ften   : tendency for variable f                                c
!                                                                     c
!     qcb3d  : variable f at t-1 time step                            c
!                                                                     c
!     psb    : p* at t-1 time step                                    c
!                                                                     c
!     xkc    : horizontal diffusion coefficient                       c
!                                                                     c
!     c203   : = c203 defined in 'param'                              c
!     c203   : = (100.-ptop)/(dx*dx), defined in 'param'              c
!                                                                     c
!     j      : j'th slice of variable qcb3d                           c
!                                                                     c
!     iend   : = iym2 for cross-point variables                       c
!              = iym1  for dot-point   variables                       c
!                                                                     c
!     jend   : = jxm2 for cross-point variables                       c
!              = jlx  for dot-point   variables                       c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use mod_dynparam
      use mod_slice
      use mod_main
      implicit none
!
! Dummy arguments
!
      real(8) :: c203
      integer :: j
      real(8) , dimension(iy,kz) :: ften , xkc
      intent (in) c203 , j , xkc
      intent (inout) ften
!
! Local variables
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
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine diffutch(ften,xkc,c203,n,j)
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine computes the diffusion term for cross-point     c
!     variable 'f' on a constant sigma surface.                       c
!                                                                     c
!     ---fourth-order diffusion scheme is applied for the interior    c
!        and second-order scheme is applied for the boundary.         c
!                                                                     c
!     ften   : tendency for variable f                                c
!                                                                     c
!     chib3d : variable f at t-1 time step                            c
!                                                                     c
!     psb    : p* at t-1 time step                                    c
!                                                                     c
!     xkc    : horizontal diffusion coefficient                       c
!                                                                     c
!     c203   : = c203 defined in 'param'                              c
!     c203   : = (100.-ptop)/(dx*dx), defined in 'param'              c
!                                                                     c
!     j      : j'th slice of variable chib3d                          c
!                                                                     c
!     iend   : = iym2 for cross-point variables                       c
!              = iym1  for dot-point   variables                       c
!                                                                     c
!     jend   : = jxm2 for cross-point variables                       c
!              = jlx  for dot-point   variables                       c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      use mod_dynparam
      use mod_slice
      use mod_main
      implicit none
!
! Dummy arguments
!
      real(8) :: c203
      integer :: j , n
      real(8) , dimension(iy,kz) :: ften , xkc
      intent (in) c203 , j , n , xkc
      intent (inout) ften
!
! Local variables
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

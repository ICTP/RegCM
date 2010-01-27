C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C    This file is part of RegCM model.
C
C    RegCM model is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    RegCM model is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
C
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      subroutine diffu_U(ften,xkc,c203,j,ind)
c                                                                     c
c     this subroutine computes the diffusion term for decoupled       c
c     variable on constant sigma surface.                             c
c                                                                     c
c     ften    : tendency for variable ubd3d                           c
c                                                                     c
c     ubd3d   : coupled variable at time t-1                          c
c                                                                     c
c     xkc     : horizontal diffusion coefficient                      c
c                                                                     c
c     msfd    : map scale factor at the points where ubd3d is defined c
c                                                                     c
c     scr     : dummy array used for working space                    c
c                                                                     c
c     c203    : 1./(dx*dx), defined in 'param'                        c
c                                                                     c
c     j       : j'th slice of variable ubd3d                          c
c                                                                     c
c     iend    : = ilxm for cross-point variables                      c
c               = ilx  for dot-point   variables                      c
c                                                                     c
c     jend    : = jlxm for cross-point variables                      c
c               = jlx  for dot-point   variables                      c
c                                                                     c
c     ind = 1 : ubd3d is already multiplied by map scale factor (msfd)c
c         = 0 : ubd3d is "not"   multiplied by msfd                   c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer j,ind
#     include "parameter.inc"
#ifdef MPP1
#     include "mpiregcm.inc"
#endif
#     include "include/slice.cb"
#     include "include/main.cb"

      real(kind=8)  ften(ix,kx),xkc(ix,kx)
      real(kind=8)  c203
      integer i,k
c
c---------------------------------------------------------------------
c
#ifdef MPP1
      if ((myid.eq.0 .and. j.eq.2).or.
     &    (myid.eq.nproc-1 .and. j.eq.jendx)) then
#else
      if (j .eq. 2 .or. j .eq. jlx) then
#endif
c
c......second-order scheme for east or west boundary:
       do 20 k=1,kx
       do 20 i=2,ilx
        if(ind.eq.0) then
          ften(i,k)=ften(i,k)+xkc(i,k)*c203*
     &     (ubd3d(i,k,j+1)+ubd3d(i,k,j-1)+ubd3d(i+1,k,j)+ubd3d(i-1,k,j)
     &  -4.*ubd3d(i,k,j))*pdotb(i,j)
        else
          ften(i,k)=ften(i,k)+xkc(i,k)*c203*
     &        (ubd3d(i,k,j+1)/msfd(i,j+1)+ubd3d(i,k,j-1)/msfd(i,j-1)
     &        +ubd3d(i+1,k,j)/msfd(i+1,j)+ubd3d(i-1,k,j)/msfd(i-1,j)
     &     -4.*ubd3d(i,k,j)/msfd(i,j))*pdotb(i,j)
        end if
20     continue
c
      else
c
c.....fourth-order scheme for interior:
       do 50 k=1,kx
       do 50 i=3,ilx-1
        if(ind.eq.0) then
          ften(i,k)=ften(i,k)-xkc(i,k)*c203*
     &     (ubd3d(i,k,j+2)+ubd3d(i,k,j-2)+ubd3d(i+2,k,j)+ubd3d(i-2,k,j)
     & -4.*(ubd3d(i,k,j+1)+ubd3d(i,k,j-1)+ubd3d(i+1,k,j)+ubd3d(i-1,k,j))
     & +12.*ubd3d(i,k,j))*pdotb(i,j)
        else
          ften(i,k)=ften(i,k)-xkc(i,k)*c203*
     &        (ubd3d(i,k,j+2)/msfd(i,j+2)+ubd3d(i,k,j-2)/msfd(i,j-2)
     &        +ubd3d(i+2,k,j)/msfd(i+2,j)+ubd3d(i-2,k,j)/msfd(i-2,j)
     &   - 4.*(ubd3d(i,k,j+1)/msfd(i,j+1)+ubd3d(i,k,j-1)/msfd(i,j-1)
     &        +ubd3d(i+1,k,j)/msfd(i+1,j)+ubd3d(i-1,k,j)/msfd(i-1,j))
     &   + 12.*ubd3d(i,k,j)/msfd(i,j))*pdotb(i,j)
        end if
50     continue
c......second-order scheme for north and south boundaries:
       do 60 i=2,ilx,ilx-2
       do 60 k=1,kx
        if(ind.eq.0) then
          ften(i,k)=ften(i,k)+xkc(i,k)*c203*
     &     (ubd3d(i,k,j+1)+ubd3d(i,k,j-1)+ubd3d(i+1,k,j)+ubd3d(i-1,k,j)
     &  -4.*ubd3d(i,k,j))*pdotb(i,j)
        else
          ften(i,k)=ften(i,k)+xkc(i,k)*c203*
     &        (ubd3d(i,k,j+1)/msfd(i,j+1)+ubd3d(i,k,j-1)/msfd(i,j-1)
     &        +ubd3d(i+1,k,j)/msfd(i+1,j)+ubd3d(i-1,k,j)/msfd(i-1,j)
     &     -4.*ubd3d(i,k,j)/msfd(i,j))*pdotb(i,j)
        end if
60     continue
c
      end if
c
      return
      end subroutine diffu_U
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine diffu_V(ften,xkc,c203,j,ind)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine computes the diffusion term for decoupled       c
c     variable on constant sigma surface.                             c
c                                                                     c
c     ften    : tendency for variable vbd3d                           c
c                                                                     c
c     vbd3d   : coupled variable at time t-1                          c
c                                                                     c
c     xkc     : horizontal diffusion coefficient                      c
c                                                                     c
c     msfd    : map scale factor at the points where vbd3d is defined c
c                                                                     c
c     scr     : dummy array used for working space                    c
c                                                                     c
c     c203    : 1./(dx*dx), defined in 'param'                        c
c                                                                     c
c     j       : j'th slice of variable vbd3d                          c
c                                                                     c
c     iend    : = ilxm for cross-point variables                      c
c               = ilx  for dot-point   variables                      c
c                                                                     c
c     jend    : = jlxm for cross-point variables                      c
c               = jlx  for dot-point   variables                      c
c                                                                     c
c     ind = 1 : vbd3d is already multiplied by map scale factor (msfd)c
c         = 0 : vbd3d is "not"   multiplied by msfd                   c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer j,ind
#     include "parameter.inc"
#     include "include/slice.cb"
#     include "include/main.cb"
#ifdef MPP1
#     include "mpiregcm.inc"
#endif
      real(kind=8)  ften(ix,kx),xkc(ix,kx)
      real(kind=8)  c203
      integer i,k
c
c---------------------------------------------------------------------
c
#ifdef MPP1
      if ((myid.eq.0 .and. j.eq.2).or.
     &    (myid.eq.nproc-1 .and. j.eq.jendx)) then
#else
      if (j .eq. 2 .or. j .eq. jlx) then
#endif
c
c......second-order scheme for east or west boundary:
       do 20 k=1,kx
       do 20 i=2,ilx
        if(ind.eq.0) then
          ften(i,k)=ften(i,k)+xkc(i,k)*c203*
     &     (vbd3d(i,k,j+1)+vbd3d(i,k,j-1)+vbd3d(i+1,k,j)+vbd3d(i-1,k,j)
     &  -4.*vbd3d(i,k,j))*pdotb(i,j)
        else
          ften(i,k)=ften(i,k)+xkc(i,k)*c203*
     &        (vbd3d(i,k,j+1)/msfd(i,j+1)+vbd3d(i,k,j-1)/msfd(i,j-1)
     &        +vbd3d(i+1,k,j)/msfd(i+1,j)+vbd3d(i-1,k,j)/msfd(i-1,j)
     &     -4.*vbd3d(i,k,j)/msfd(i,j))*pdotb(i,j)
        end if
20     continue
c
      else
c
c.....fourth-order scheme for interior:
       do 50 k=1,kx
       do 50 i=3,ilx-1
        if(ind.eq.0) then
          ften(i,k)=ften(i,k)-xkc(i,k)*c203*
     &     (vbd3d(i,k,j+2)+vbd3d(i,k,j-2)+vbd3d(i+2,k,j)+vbd3d(i-2,k,j)
     & -4.*(vbd3d(i,k,j+1)+vbd3d(i,k,j-1)+vbd3d(i+1,k,j)+vbd3d(i-1,k,j))
     & +12.*vbd3d(i,k,j))*pdotb(i,j)
        else
          ften(i,k)=ften(i,k)-xkc(i,k)*c203*
     &        (vbd3d(i,k,j+2)/msfd(i,j+2)+vbd3d(i,k,j-2)/msfd(i,j-2)
     &        +vbd3d(i+2,k,j)/msfd(i+2,j)+vbd3d(i-2,k,j)/msfd(i-2,j)
     &   - 4.*(vbd3d(i,k,j+1)/msfd(i,j+1)+vbd3d(i,k,j-1)/msfd(i,j-1)
     &        +vbd3d(i+1,k,j)/msfd(i+1,j)+vbd3d(i-1,k,j)/msfd(i-1,j))
     &   + 12.*vbd3d(i,k,j)/msfd(i,j))*pdotb(i,j)
        end if
50     continue
c......second-order scheme for north and south boundaries:
       do 60 i=2,ilx,ilx-2
       do 60 k=1,kx
        if(ind.eq.0) then
          ften(i,k)=ften(i,k)+xkc(i,k)*c203*
     &     (vbd3d(i,k,j+1)+vbd3d(i,k,j-1)+vbd3d(i+1,k,j)+vbd3d(i-1,k,j)
     &  -4.*vbd3d(i,k,j))*pdotb(i,j)
        else
          ften(i,k)=ften(i,k)+xkc(i,k)*c203*
     &        (vbd3d(i,k,j+1)/msfd(i,j+1)+vbd3d(i,k,j-1)/msfd(i,j-1)
     &        +vbd3d(i+1,k,j)/msfd(i+1,j)+vbd3d(i-1,k,j)/msfd(i-1,j)
     &     -4.*vbd3d(i,k,j)/msfd(i,j))*pdotb(i,j)
        end if
60     continue
c
      end if
c
      return
      end subroutine diffu_V

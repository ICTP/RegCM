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

      subroutine diffut_T(ften,xkc,c203,j)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine computes the diffusion term for cross-point     c
c     variable 'f' on a constant sigma surface.                       c
c                                                                     c
c     ---fourth-order diffusion scheme is applied for the interior    c
c        and second-order scheme is applied for the boundary.         c
c                                                                     c
c     ften   : tendency for variable f                                c
c                                                                     c
c     tb3d   : variable f at t-1 time step                            c
c                                                                     c
c     psb    : p* at t-1 time step                                    c
c                                                                     c
c     xkc    : horizontal diffusion coefficient                       c
c                                                                     c
c     c203   : = c203 defined in 'param'                              c
c     c203   : = (100.-ptop)/(dx*dx), defined in 'param'              c
c                                                                     c
c     j      : j'th slice of variable tb3d                            c
c                                                                     c
c     iend   : = ilxm for cross-point vaeiables                       c
c              = ilx  for dot-point   variables                       c
c                                                                     c
c     jend   : = jlxm for cross-point variables                       c
c              = jlx  for dot-point   variables                       c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c
      integer j
#     include "parameter.inc"
#ifdef MPP1
#     include "mpiregcm.inc"
#endif
#     include "include/slice.cb"
#     include "include/main.cb"
c
      real(kind=8)  ften(ix,kx),xkc(ix,kx)
      real(kind=8)  c203
      integer i,k
c
c----------------------------------------------------------------------
c-----compute the diffusion term for t:
c
#ifdef MPP1
      if ((myid.eq.0 .and. j.eq.2).or.
     &    (myid.eq.nproc-1 .and. j.eq.jendm)) then
#else
      if (j .eq. 2 .or. j .eq. jlxm) then
#endif
c
c......second-order scheme for east or west boundary:
       do 25 k=1,kx
       do 25 i=2,ilxm
        ften(i,k)=ften(i,k)+xkc(i,k)*c203*(tb3d(i,k,j+1)+tb3d(i,k,j-1)+
     &          tb3d(i+1,k,j)+tb3d(i-1,k,j) - 4.*tb3d(i,k,j))*psb(i,j)
25     continue
c
      else
c
c......fourth-order scheme for interior:
       do 50 k=1,kx
       do 50 i=3,ilxm-1
        ften(i,k)=ften(i,k)-xkc(i,k)*c203*(tb3d(i,k,j+2)+tb3d(i,k,j-2)+
     &          tb3d(i+2,k,j)+tb3d(i-2,k,j)
     &                               - 4.*(tb3d(i,k,j+1)+tb3d(i,k,j-1)+
     &          tb3d(i+1,k,j)+tb3d(i-1,k,j)) + 12.*tb3d(i,k,j))*psb(i,j)
50     continue
c......second-order scheme for north and south boundaries:
       i = 2
       do 60 k=1,kx
        ften(i,k)=ften(i,k)+xkc(i,k)*c203*(tb3d(i,k,j+1)+tb3d(i,k,j-1)+
     &          tb3d(i+1,k,j)+tb3d(i-1,k,j) - 4.*tb3d(i,k,j))*psb(i,j)
60     continue
       i = ilxm
       do 70 k=1,kx
        ften(i,k)=ften(i,k)+xkc(i,k)*c203*(tb3d(i,k,j+1)+tb3d(i,k,j-1)+
     &          tb3d(i+1,k,j)+tb3d(i-1,k,j) - 4.*tb3d(i,k,j))*psb(i,j)
   70  continue
c
      end if
c
      return
      end subroutine diffut_T
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine diffutQV(ften,xkc,c203,j)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine computes the diffusion term for cross-point     c
c     variable 'f' on a constant sigma surface.                       c
c                                                                     c
c     ---fourth-order diffusion scheme is applied for the interior    c
c        and second-order scheme is applied for the boundary.         c
c                                                                     c
c     ften   : tendency for variable f                                c
c                                                                     c
c     qvb3d  : variable f at t-1 time step                            c
c                                                                     c
c     psb    : p* at t-1 time step                                    c
c                                                                     c
c     xkc    : horizontal diffusion coefficient                       c
c                                                                     c
c     c203   : = c203 defined in 'param'                              c
c     c203   : = (100.-ptop)/(dx*dx), defined in 'param'              c
c                                                                     c
c     j      : j'th slice of variable qvb3d                           c
c                                                                     c
c     iend   : = ilxm for cross-point vaeiables                       c
c              = ilx  for dot-point   variables                       c
c                                                                     c
c     jend   : = jlxm for cross-point variables                       c
c              = jlx  for dot-point   variables                       c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c
      integer j
#     include "parameter.inc"
#     include "include/slice.cb"
#     include "include/main.cb"
#ifdef MPP1
#     include "mpiregcm.inc"
#endif
c
      real(kind=8)  ften(ix,kx),xkc(ix,kx)
      real(kind=8)  c203
      integer i,k
c
c----------------------------------------------------------------------
c-----compute the diffusion term for t:
c
#ifdef MPP1
      if ((myid.eq.0 .and. j.eq.2).or.
     &    (myid.eq.nproc-1 .and. j.eq.jendm)) then
#else
      if (j .eq. 2 .or. j .eq. jlxm) then
#endif
c
c......second-order scheme for east or west boundary:
       do 25 k=1,kx
       do 25 i=2,ilxm
        ften(i,k)=ften(i,k)+xkc(i,k)*c203*
     &      (qvb3d(i,k,j+1)+qvb3d(i,k,j-1)+
     &       qvb3d(i+1,k,j)+qvb3d(i-1,k,j) - 4.*qvb3d(i,k,j))*psb(i,j)
25     continue
c
      else
c
c......fourth-order scheme for interior:
       do 50 k=1,kx
       do 50 i=3,ilxm-1
        ften(i,k)=ften(i,k)-xkc(i,k)*c203*
     &      (qvb3d(i,k,j+2)+qvb3d(i,k,j-2)+
     &       qvb3d(i+2,k,j)+qvb3d(i-2,k,j)
     & - 4.*(qvb3d(i,k,j+1)+qvb3d(i,k,j-1)+
     &       qvb3d(i+1,k,j)+qvb3d(i-1,k,j)) + 12.*qvb3d(i,k,j))*psb(i,j)
50     continue
c......second-order scheme for north and south boundaries:
       i = 2
       do 60 k=1,kx
        ften(i,k)=ften(i,k)+xkc(i,k)*c203*
     &      (qvb3d(i,k,j+1)+qvb3d(i,k,j-1)+
     &       qvb3d(i+1,k,j)+qvb3d(i-1,k,j) - 4.*qvb3d(i,k,j))*psb(i,j)
60     continue
       i = ilxm
       do 70 k=1,kx
        ften(i,k)=ften(i,k)+xkc(i,k)*c203*
     &      (qvb3d(i,k,j+1)+qvb3d(i,k,j-1)+
     &       qvb3d(i+1,k,j)+qvb3d(i-1,k,j) - 4.*qvb3d(i,k,j))*psb(i,j)
   70  continue
c
      end if
c
      return
      end subroutine diffutQV
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine diffutQC(ften,xkc,c203,j)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine computes the diffusion term for cross-point     c
c     variable 'f' on a constant sigma surface.                       c
c                                                                     c
c     ---fourth-order diffusion scheme is applied for the interior    c
c        and second-order scheme is applied for the boundary.         c
c                                                                     c
c     ften   : tendency for variable f                                c
c                                                                     c
c     qcb3d  : variable f at t-1 time step                            c
c                                                                     c
c     psb    : p* at t-1 time step                                    c
c                                                                     c
c     xkc    : horizontal diffusion coefficient                       c
c                                                                     c
c     c203   : = c203 defined in 'param'                              c
c     c203   : = (100.-ptop)/(dx*dx), defined in 'param'              c
c                                                                     c
c     j      : j'th slice of variable qcb3d                           c
c                                                                     c
c     iend   : = ilxm for cross-point vaeiables                       c
c              = ilx  for dot-point   variables                       c
c                                                                     c
c     jend   : = jlxm for cross-point variables                       c
c              = jlx  for dot-point   variables                       c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c
      integer j
#     include "parameter.inc"
#     include "include/slice.cb"
#     include "include/main.cb"
#ifdef MPP1
#     include "mpiregcm.inc"
#endif
c
      real(kind=8)  ften(ix,kx),xkc(ix,kx)
      real(kind=8)  c203
      integer i,k
c
c----------------------------------------------------------------------
c-----compute the diffusion term for t:
c
#ifdef MPP1
      if ((myid.eq.0 .and. j.eq.2).or.
     &    (myid.eq.nproc-1 .and. j.eq.jendm)) then
#else
      if (j .eq. 2 .or. j .eq. jlxm) then
#endif
c
c......second-order scheme for east or west boundary:
       do 25 k=1,kx
       do 25 i=2,ilxm
        ften(i,k)=ften(i,k)+xkc(i,k)*c203*
     &      (qcb3d(i,k,j+1)+qcb3d(i,k,j-1)+
     &       qcb3d(i+1,k,j)+qcb3d(i-1,k,j) - 4.*qcb3d(i,k,j))*psb(i,j)
25     continue
c
      else
c
c......fourth-order scheme for interior:
       do 50 k=1,kx
       do 50 i=3,ilxm-1
        ften(i,k)=ften(i,k)-xkc(i,k)*c203*
     &      (qcb3d(i,k,j+2)+qcb3d(i,k,j-2)+
     &       qcb3d(i+2,k,j)+qcb3d(i-2,k,j)
     & - 4.*(qcb3d(i,k,j+1)+qcb3d(i,k,j-1)+
     &       qcb3d(i+1,k,j)+qcb3d(i-1,k,j)) + 12.*qcb3d(i,k,j))*psb(i,j)
50     continue
c......second-order scheme for north and south boundaries:
       i = 2
       do 60 k=1,kx
        ften(i,k)=ften(i,k)+xkc(i,k)*c203*
     &      (qcb3d(i,k,j+1)+qcb3d(i,k,j-1)+
     &       qcb3d(i+1,k,j)+qcb3d(i-1,k,j) - 4.*qcb3d(i,k,j))*psb(i,j)
60     continue
       i = ilxm
       do 70 k=1,kx
        ften(i,k)=ften(i,k)+xkc(i,k)*c203*
     &      (qcb3d(i,k,j+1)+qcb3d(i,k,j-1)+
     &       qcb3d(i+1,k,j)+qcb3d(i-1,k,j) - 4.*qcb3d(i,k,j))*psb(i,j)
   70  continue
c
      end if
c
      return
      end subroutine diffutQC
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine diffutCH(ften,xkc,c203,n,j)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine computes the diffusion term for cross-point     c
c     variable 'f' on a constant sigma surface.                       c
c                                                                     c
c     ---fourth-order diffusion scheme is applied for the interior    c
c        and second-order scheme is applied for the boundary.         c
c                                                                     c
c     ften   : tendency for variable f                                c
c                                                                     c
c     chib3d : variable f at t-1 time step                            c
c                                                                     c
c     psb    : p* at t-1 time step                                    c
c                                                                     c
c     xkc    : horizontal diffusion coefficient                       c
c                                                                     c
c     c203   : = c203 defined in 'param'                              c
c     c203   : = (100.-ptop)/(dx*dx), defined in 'param'              c
c                                                                     c
c     j      : j'th slice of variable chib3d                          c
c                                                                     c
c     iend   : = ilxm for cross-point vaeiables                       c
c              = ilx  for dot-point   variables                       c
c                                                                     c
c     jend   : = jlxm for cross-point variables                       c
c              = jlx  for dot-point   variables                       c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c
      integer j
#     include "parameter.inc"
#     include "include/slice.cb"
#     include "include/main.cb"
#ifdef MPP1
#     include "mpiregcm.inc"
#endif
c
      real(kind=8)  ften(ix,kx),xkc(ix,kx)
      real(kind=8)  c203
      integer n,i,k
c
c----------------------------------------------------------------------
c-----compute the diffusion term for t:
c
#ifdef MPP1
      if ((myid.eq.0 .and. j.eq.2).or.
     &    (myid.eq.nproc-1 .and. j.eq.jendm)) then
#else
      if (j .eq. 2 .or. j .eq. jlxm) then
#endif
c
c......second-order scheme for east or west boundary:
       do 25 k=1,kx
       do 25 i=2,ilxm
        ften(i,k)=ften(i,k)+xkc(i,k)*c203*
     &      (chib3d(i,k,j+1,n)+chib3d(i,k,j-1,n)+
     &       chib3d(i+1,k,j,n)+chib3d(i-1,k,j,n)
     &  - 4.*chib3d(i,k,j,n))*psb(i,j)
25     continue
c
      else
c
c......fourth-order scheme for interior:
       do 50 k=1,kx
       do 50 i=3,ilxm-1
        ften(i,k)=ften(i,k)-xkc(i,k)*c203*
     &      (chib3d(i,k,j+2,n)+chib3d(i,k,j-2,n)+
     &       chib3d(i+2,k,j,n)+chib3d(i-2,k,j,n)
     & - 4.*(chib3d(i,k,j+1,n)+chib3d(i,k,j-1,n)+
     &       chib3d(i+1,k,j,n)+chib3d(i-1,k,j,n))
     & + 12.*chib3d(i,k,j,n))*psb(i,j)
50     continue
c......second-order scheme for north and south boundaries:
       i = 2
       do 60 k=1,kx
        ften(i,k)=ften(i,k)+xkc(i,k)*c203*
     &      (chib3d(i,k,j+1,n)+chib3d(i,k,j-1,n)+
     &       chib3d(i+1,k,j,n)+chib3d(i-1,k,j,n)
     &  - 4.*chib3d(i,k,j,n))*psb(i,j)
60     continue
       i = ilxm
       do 70 k=1,kx
        ften(i,k)=ften(i,k)+xkc(i,k)*c203*
     &      (chib3d(i,k,j+1,n)+chib3d(i,k,j-1,n)+
     &       chib3d(i+1,k,j,n)+chib3d(i-1,k,j,n)
     &  - 4.*chib3d(i,k,j,n))*psb(i,j)
   70  continue
c
      end if
c
      return
      end subroutine diffutCH

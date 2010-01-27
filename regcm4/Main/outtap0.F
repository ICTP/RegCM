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

      subroutine outtap0
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine writes the model output to tape or disk for use c
c     in dataflow analyses.                                           c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
#     include "parameter.inc"
#     include "include/param1.cb"
#     include "include/param2.cb"
#     include "include/param3.cb"
#     include "include/iunits.cb"
#ifdef MPP1
#     include "include/io.cb"
      real(kind=4)  fout(mjx-2,ix-2)
#else
#     include "include/main.cb"
      real(kind=4)  fout(jx-2,ix-2)
#endif
#     include "include/bats.cb2"
#     include "include/cvaria.cb"
#     include "include/date.cb"
      real(kind=4)  sp1d(kxp1),dto,dtb,dtr,dtc
      real(kind=4)  clat,clon,dxsp,ptsp,plat,plon
      character proj*6
      COMMON /GRADS/clat,clon,dxsp,ptsp,plat,plon,proj
c
      integer i,j,k
c
c---------------------------------------------------------------------
c-----output large-domain variables:
c
#ifdef MPP1
      open(20,file='output/OUT_HEAD',form='unformatted'
     &    ,recl=(ix-2)*(mjx-2)*ibyte,access='direct')
#else
      open(20,file='output/OUT_HEAD',form='unformatted'
     &    ,recl=(ix-2)*(jx-2)*ibyte,access='direct')
#endif
      do k=1,kxp1
        sp1d(k) = sigma(kxp1-k+1)
      end do
      dto = tapfrq
      dtb = batfrq
      dtr = radisp
      dtc = chemfrq
#ifdef MPP1
      write(20,rec=1) mdate0,ibltyp,icup,ipptls,iboudy
     &    ,ix,mjx,kx,sp1d,dxsp,ptsp,clat,clon,plat,plon,proj
     &    ,dto,dtb,dtr,dtc,iotyp
#else
      write(20,rec=1) mdate0,ibltyp,icup,ipptls,iboudy
     &    ,ix,jx,kx,sp1d,dxsp,ptsp,clat,clon,plat,plon,proj
     &    ,dto,dtb,dtr,dtc,iotyp
#endif
      do i=1,ix-2
#ifdef MPP1
        do j=1,mjx-2
          fout(j,i) =  ht_io(i+1,j+1)/g
#else
        do j=1,jx-2
          fout(j,i) =  ht(i+1,j+1)/g
#endif
        end do
      end do
      write(20,rec=2) fout
      do i=1,ix-2
#ifdef MPP1
        do j=1,mjx-2
          fout(j,i) =  htsd_io(i+1,j+1)
#else
        do j=1,jx-2
          fout(j,i) =  htsd(i+1,j+1)
#endif
        end do
      end do
      write(20,rec=3) fout
      do i=1,ix-2
#ifdef MPP1
        do j=1,mjx-2
          fout(j,i) =  veg2d_io(i+1,j+1)
#else
        do j=1,jx-2
          fout(j,i) =  veg2d(i+1,j+1)
#endif
        end do
      end do
      write(20,rec=4) fout
      do i=1,ix-2
#ifdef MPP1
        do j=1,mjx-2
          fout(j,i) =  satbrt_io(i+1,j+1)
#else
        do j=1,jx-2
          fout(j,i) =  satbrt(i+1,j+1)
#endif
        end do
      end do
      write(20,rec=5) fout
      do i=1,ix-2
#ifdef MPP1
        do j=1,mjx-2
          fout(j,i) =  xlat_io(i+1,j+1)
#else
        do j=1,jx-2
          fout(j,i) =  xlat(i+1,j+1)
#endif
        end do
      end do
      write(20,rec=6) fout
      do i=1,ix-2
#ifdef MPP1
        do j=1,mjx-2
          fout(j,i) =  xlong_io(i+1,j+1)
#else
        do j=1,jx-2
          fout(j,i) =  xlong(i+1,j+1)
#endif
        end do
      end do
      write(20,rec=7) fout
      do i=1,ix-2
#ifdef MPP1
        do j=1,mjx-2
          fout(j,i) =  1./msfx_io(i+1,j+1)
#else
        do j=1,jx-2
          fout(j,i) =  1./msfx(i+1,j+1)
#endif
        end do
      end do
      write(20,rec=8) fout
      do i=1,ix-2
#ifdef MPP1
        do j=1,mjx-2
          fout(j,i) =  1./msfd_io(i+1,j+1)
#else
        do j=1,jx-2
          fout(j,i) =  1./msfd(i+1,j+1)
#endif
        end do
      end do
      write(20,rec=9) fout
      do i=1,ix-2
#ifdef MPP1
        do j=1,mjx-2
          fout(j,i) =  f_io(i+1,j+1)
#else
        do j=1,jx-2
          fout(j,i) =  f(i+1,j+1)
#endif
        end do
      end do
      write(20,rec=10) fout
      do i=1,ix-2
#ifdef MPP1
        do j=1,mjx-2
          if(satbrt_io(i+1,j+1).gt.13.5.and.
     &       satbrt_io(i+1,j+1).lt.15.5) then
#else
        do j=1,jx-2
          if(satbrt(i+1,j+1).gt.13.5.and.satbrt(i+1,j+1).lt.15.5) then
#endif
            fout(j,i) = 0.
          else
            fout(j,i) = 2.
          end if
        end do
      end do
      write(20,rec=11) fout

      close(20)
c
      return
      end subroutine outtap0

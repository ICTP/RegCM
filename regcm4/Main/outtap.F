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

      subroutine outtap
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
      real(kind=8)  out1(mjx,ix)
      real(kind=4)  fout(mjx-2,ix-2)
#else
#     include "include/main.cb"
#     include "include/cvaria.cb"
#     include "include/bats2.cb"
      real(kind=8)  out1(jx,ix)
      real(kind=4)  fout(jx-2,ix-2)
#endif
#     include "include/date.cb"

      integer nrcout
      COMMON /GENREC/ nrcout
      real(kind=8)  mmpd
c
      integer i,j,k,n,nn
      real(kind=8)  xh
c
c---------------------------------------------------------------------
c-----output RegCM-domain variables:
c
      xh=xtime/60.
      if (iotyp.eq.2) write (iutdat) idatex
c
c ******  write one time of data to mm4 output file.

      do k=kx,1,-1
         do i=1,ix
#ifdef MPP1
           do j=1,mjx
             out1(j,i) = ua_io(i,k,j)
#else
           do j=1,jx
             out1(j,i) = ua(i,k,j)
#endif
           end do
         end do
         do i=1,ix-2
#ifdef MPP1
           do j=1,mjx-2
             fout(j,i) = 0.25*(out1(j+1,i+1)+out1(j+2,i+1)
     &                        +out1(j+1,i+2)+out1(j+2,i+2))
     &                        /psa_io(i+1,j+1)
#else
           do j=1,jx-2
             fout(j,i) = 0.25*(out1(j+1,i+1)+out1(j+2,i+1)
     &                        +out1(j+1,i+2)+out1(j+2,i+2))
     &                        /psa(i+1,j+1)
#endif
           end do
         end do 
         if (iotyp.eq.1) then
           nrcout=nrcout+1
           write(iutdat,rec=nrcout) fout
         else if (iotyp.eq.2) then
           write(iutdat) fout
         end if
      end do
      do k=kx,1,-1
         do i=1,ix
#ifdef MPP1
           do j=1,mjx
             out1(j,i) = va_io(i,k,j)
#else
           do j=1,jx
             out1(j,i) = va(i,k,j)
#endif
           end do
         end do
         do i=1,ix-2
#ifdef MPP1
           do j=1,mjx-2
             fout(j,i) = 0.25*(out1(j+1,i+1)+out1(j+2,i+1)
     &                        +out1(j+1,i+2)+out1(j+2,i+2))
     &                       /psa_io(i+1,j+1)
#else
           do j=1,jx-2
             fout(j,i) = 0.25*(out1(j+1,i+1)+out1(j+2,i+1)
     &                        +out1(j+1,i+2)+out1(j+2,i+2))
     &                       /psa(i+1,j+1)
#endif
           end do
         end do
         if (iotyp.eq.1) then
           nrcout=nrcout+1
           write(iutdat,rec=nrcout) fout
         else if (iotyp.eq.2) then
           write(iutdat) fout
         end if
      end do
      do k=kx,1,-1
        do i=1,ix-2 
#ifdef MPP1
          do j=1,mjx-2
            fout(j,i) = omega_io(i+1,k,j+1)
#else
          do j=1,jx-2
            fout(j,i) = omega(i+1,k,j+1)
#endif
          end do
        end do
        if (iotyp.eq.1) then
          nrcout=nrcout+1
          write(iutdat,rec=nrcout) fout
        else if (iotyp.eq.2) then
          write(iutdat) fout
        end if
      end do
      do k=kx,1,-1
        do i=1,ix-2 
#ifdef MPP1
          do j=1,mjx-2
            fout(j,i) = ta_io(i+1,k,j+1)/psa_io(i+1,j+1)
#else
          do j=1,jx-2
            fout(j,i) = ta(i+1,k,j+1)/psa(i+1,j+1)
#endif
          end do
        end do
        if (iotyp.eq.1) then
          nrcout=nrcout+1
          write(iutdat,rec=nrcout) fout
        else if (iotyp.eq.2) then
          write(iutdat) fout
        end if
      end do
      do k=kx,1,-1
        do i=1,ix-2 
#ifdef MPP1
          do j=1,mjx-2
            fout(j,i) = qva_io(i+1,k,j+1)/psa_io(i+1,j+1)
#else
          do j=1,jx-2
            fout(j,i) = qva(i+1,k,j+1)/psa(i+1,j+1)
#endif
          end do
        end do
        if (iotyp.eq.1) then
          nrcout=nrcout+1
          write(iutdat,rec=nrcout) fout
        else if (iotyp.eq.2) then
          write(iutdat) fout
        end if
      end do
      do k=kx,1,-1
        do i=1,ix-2 
#ifdef MPP1
          do j=1,mjx-2
            fout(j,i) = qca_io(i+1,k,j+1)/psa_io(i+1,j+1)
#else
          do j=1,jx-2
            fout(j,i) = qca(i+1,k,j+1)/psa(i+1,j+1)
#endif
          end do
        end do
        if (iotyp.eq.1) then
          nrcout=nrcout+1
          write(iutdat,rec=nrcout) fout
        else if (iotyp.eq.2) then
          write(iutdat) fout
        end if
      end do
      do i=1,ix-2 
#ifdef MPP1
        do j=1,mjx-2
          fout(j,i) = (psa_io(i+1,j+1)+ptop)*10.
#else
        do j=1,jx-2
          fout(j,i) = (psa(i+1,j+1)+ptop)*10.
#endif
        end do
      end do
      if (iotyp.eq.1) then
        nrcout=nrcout+1
        write(iutdat,rec=nrcout) fout
      else if (iotyp.eq.2) then
        write(iutdat) fout
      end if
      mmpd = 24./tapfrq
      do i=1,ix-2
#ifdef MPP1
        do j=1,mjx-2
          fout(j,i) = (rainc_io(i+1,j+1)+rainnc_io(i+1,j+1))*mmpd
#else
        do j=1,jx-2
          fout(j,i) = (rainc(i+1,j+1)+rainnc(i+1,j+1))*mmpd
#endif
        end do 
      end do 
      if (iotyp.eq.1) then
        nrcout=nrcout+1
        write(iutdat,rec=nrcout) fout
      else if (iotyp.eq.2) then
        write(iutdat) fout
      end if

c ****** write out the following surface fields:
c                  1.  temp of lower soil layer (17)
c                  2.  total soil water in mm h2o (13)
c                  3.  accum infiltration (30)
      do i=1,ix-2
#ifdef MPP1
        do j=1,mjx-2
          fout(j,i) = tgb2d_io(1,i+1,j+1)
          do n=2,NNSG
            fout(j,i) = fout(j,i) + tgb2d_io(n,i+1,j+1)
          end do
          fout(j,i) = fout(j,i)/float(NNSG)
#else
        do j=1,jx-2
          fout(j,i) = tgb2d(1,i+1,j+1)
          do n=2,NNSG
            fout(j,i) = fout(j,i) + tgb2d(n,i+1,j+1)
          end do
          fout(j,i) = fout(j,i)/float(NNSG)
#endif
        end do 
      end do 
      if (iotyp.eq.1) then
        nrcout=nrcout+1
        write(iutdat,rec=nrcout) fout
      else if (iotyp.eq.2) then
        write(iutdat) fout
      end if

      do i=1,ix-2
#ifdef MPP1
        do j=1,mjx-2
          fout(j,i) = 0.0
          nn=0
          do n=1,NNSG
            if (ocld2d_io(n,i+1,j+1).ge.0.5) then
              fout(j,i) = fout(j,i)+swt2d_io(n,i+1,j+1)
              nn=nn+1
            end if
          end do
#else
        do j=1,jx-2
          fout(j,i) = 0.0
          nn=0
          do n=1,NNSG
            if (ocld2d(n,i+1,j+1).ge.0.5) then
              fout(j,i) = fout(j,i)+swt2d(n,i+1,j+1)
              nn=nn+1
            end if
          end do
#endif
          if(nn.ge.max0(NNSG/2,1)) then
            fout(j,i) = fout(j,i)/float(nn)
          else
            fout(j,i) = -1.e34
          end if
        end do 
      end do 
      if (iotyp.eq.1) then
        nrcout=nrcout+1
        write(iutdat,rec=nrcout) fout
      else if (iotyp.eq.2) then
        write(iutdat) fout
      end if

      do i=1,ix-2
#ifdef MPP1
        do j=1,mjx-2
          fout(j,i) = 0.0
          nn=0
          do n=1,NNSG
            if (ocld2d_io(n,i+1,j+1).ge.0.5) then
              fout(j,i) = fout(j,i)+rno2d_io(n,i+1,j+1)
              nn=nn+1
            end if
          end do
#else
        do j=1,jx-2
          fout(j,i) = 0.0
          nn=0
          do n=1,NNSG
            if (ocld2d(n,i+1,j+1).ge.0.5) then
              fout(j,i) = fout(j,i)+rno2d(n,i+1,j+1)
              nn=nn+1
            end if
          end do
#endif
          if(nn.ge.max0(NNSG/2,1)) then
            fout(j,i) = fout(j,i)/float(nn)
          else
            fout(j,i) = -1.e34
          end if
        end do 
      end do 
      if (iotyp.eq.1) then
        nrcout=nrcout+1
        write(iutdat,rec=nrcout) fout
      else if (iotyp.eq.2) then
        write(iutdat) fout
      end if

c     changes for accum infiltration

      do i=1,ix-1
#ifdef MPP1
        do j=1,mjx-1
          do n=1,NNSG
            rno2d_io(n,i,j) = 0.
          end do
          rainc_io(i,j) = 0.
          rainnc_io(i,j) = 0.
#else
        do j=1,jx-1
          do n=1,NNSG
            rno2d(n,i,j) = 0.
          end do
          rainc(i,j) = 0.
          rainnc(i,j) = 0.
#endif
        end do
      end do

      print *,'OUT-history written date = ',ldatez+xtime/1440.

      return
      end subroutine outtap

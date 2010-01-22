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


ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine conadv
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine computes the amounts of dry air and water       c
c     substance advected through the lateral boundaries.              c
c                                                                     c
c     ---the unit from advection is converted to "kg".                c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
#ifdef MPP1
      include 'regcm.param2'
      include 'mpiregcm.h'
#else
      include 'regcm.param'
#endif
      include 'parame'
#     include "include/param1.cb"
#     include "include/diagnosis.cb"
#     include "include/param2.cb"
#     include "include/param3.cb"
#     include "include/main.cb"

#ifdef MPP1
#     include "include/io.cb"
      integer ierr
      real(kind=8)  vaix(kx,jxp),vaix_g(kx,mjx)
      real(kind=8)  va01(kx,jxp),va01_g(kx,mjx)
      real(kind=8)  psailx(jxp),psailx_g(mjx),psa01(jxp),psa01_g(mjx)
      real(kind=8)  qvailx(kx,jxp),qvailx_g(kx,mjx)
      real(kind=8)  qva01(kx,jxp),qva01_g(kx,mjx)
      real(kind=8)  qcailx(kx,jxp),qcailx_g(kx,mjx)
      real(kind=8)  qca01(kx,jxp),qca01_g(kx,mjx)
#endif
c
      real(kind=8)  worka(ilx,kx),workb(ilx,kx)
      integer i,j,k
c
c----------------------------------------------------------------------
c-----advection of dry air through the lateral boundaries:
c
c.....advection through east-west boundaries:
c
#ifdef MPP1
      if (myid .eq. nproc-1) then
      do k=1,kx
      do i=1,ilx
         worka(i,k)=(ua(i+1,k,jendl)+ua(i,k,jendl))
     &             /(msfx(i,jendx)*msfx(i,jendx))
      enddo
      enddo
      endif
      if (myid .eq.0) then
      do k=1,kx
      do i=1,ilx
         workb(i,k)=(ua(i+1,k,1)+ua(i,k,1))/(msfx(i,1)*msfx(i,1))
      enddo
      enddo
      endif
      call MPI_BCAST(worka,ilx*kx,MPI_DOUBLE_PRECISION,nproc-1
     &              ,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(workb,ilx*kx,MPI_DOUBLE_PRECISION,0
     &              ,MPI_COMM_WORLD,ierr)
#else
      do k=1,kx
      do i=1,ilx
         worka(i,k)=(ua(i+1,k,jx)+ua(i,k,jx))/(msfx(i,jlx)*msfx(i,jlx))
         workb(i,k)=(ua(i+1,k,1)+ua(i,k,1))/(msfx(i,1)*msfx(i,1))
      enddo
      enddo
#endif
      do 10 k=1,kx
      do 10 i=1,ilx
       tdadv=tdadv-dtmin*3.e4*dsigma(k)*dx*(worka(i,k)-workb(i,k))/g
10    continue
c
c.....advection through north-south boundaries:
c
#ifdef MPP1
      do j=1,jendl
        do k=1,kx
          vaix(k,j)=va(ix,k,j)
          va01(k,j)=va(1 ,k,j)
        enddo
      enddo
      call MPI_GATHER(vaix(1,1),kx*jxp,MPI_DOUBLE_PRECISION,
     &              vaix_g(1,1),kx*jxp,MPI_DOUBLE_PRECISION,
     &                               0,MPI_COMM_WORLD,ierr)
      call MPI_GATHER(va01(1,1),kx*jxp,MPI_DOUBLE_PRECISION,
     &              va01_g(1,1),kx*jxp,MPI_DOUBLE_PRECISION,
     &                               0,MPI_COMM_WORLD,ierr)
      if(myid.eq.0) then
        do k=1,kx
        do j=1,mjx-1
         tdadv=tdadv-dtmin*3.e4*dsigma(k)*dx*
     &  ((vaix_g(k,j+1)+vaix_g(k,j))/(msfx_io(ilx,j)*msfx_io(ilx,j))
     &  -(va01_g(k,j+1)+va01_g(k,j))/(msfx_io(1,j)*msfx_io(1,j)))/g
        enddo
        enddo
      endif
      call MPI_BCAST(tdadv,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#else
      do k=1,kx
      do j=1,jlx
       tdadv=tdadv-dtmin*3.e4*dsigma(k)*dx*((va(ix,k,j+1)+va(ix,k,j))/
     1            (msfx(ilx,j)*msfx(ilx,j))-(va(1,k,j+1)+va(1,k,j))/
     2            (msfx(1,j)*msfx(1,j)))/g
      enddo
      enddo
#endif
c
c----------------------------------------------------------------------
c-----advection of water vapor through the lateral boundaries:
c
c.....advection through east-west boundaries:
c
#ifdef MPP1
      if (myid .eq. nproc-1) then
      do k=1,kx
      do i=1,ilx
         worka(i,k)=(ua(i+1,k,jendl)+ua(i,k,jendl))
     &             *(qva(i,k,jendx)/psa(i,jendx))
     &             /(msfx(i,jendx)*msfx(i,jendx))
      enddo
      enddo
      endif
      if (myid .eq.0) then
      do k=1,kx
      do i=1,ilx
         workb(i,k)=(ua(i+1,k,1)+ua(i,k,1))*(qva(i,k,1)/psa(i,1))
     &             /(msfx(i,1)*msfx(i,1))
      enddo
      enddo
      endif
      call MPI_BCAST(worka,ilx*kx,MPI_DOUBLE_PRECISION,nproc-1
     &              ,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(workb,ilx*kx,MPI_DOUBLE_PRECISION,0
     &              ,MPI_COMM_WORLD,ierr)
#else
      do k=1,kx
      do i=1,ilx
         worka(i,k)=(ua(i+1,k,jx)+ua(i,k,jx))*(qva(i,k,jlx)/psa(i,jlx))
     &             /(msfx(i,jlx)*msfx(i,jlx))
         workb(i,k)=(ua(i+1,k,1)+ua(i,k,1))*(qva(i,k,1)/psa(i,1))
     &             /(msfx(i,1)*msfx(i,1))
      enddo
      enddo
#endif
      do 30 k=1,kx
      do 30 i=1,ilx
       tqadv=tqadv-dtmin*3.e4*dsigma(k)*dx*(worka(i,k)-workb(i,k))/g
30    continue
c
c....advection through north-south boundaries:
c
#ifdef MPP1
      do j=1,jendl
        do k=1,kx
          qvailx(k,j)=qva(ilx,k,j)
          qva01(k,j) =qva(1  ,k,j)
        enddo
        psailx(j)=psa(ilx,j)
        psa01(j)=psa(1,j)
      enddo
      call MPI_GATHER(qvailx(1,1),kx*jxp,MPI_DOUBLE_PRECISION,
     &              qvailx_g(1,1),kx*jxp,MPI_DOUBLE_PRECISION,
     &                                0,MPI_COMM_WORLD,ierr)
      call MPI_GATHER(qva01(1,1),kx*jxp,MPI_DOUBLE_PRECISION,
     &              qva01_g(1,1),kx*jxp,MPI_DOUBLE_PRECISION,
     &                                0,MPI_COMM_WORLD,ierr)
      call MPI_GATHER(psailx(1),jxp,MPI_DOUBLE_PRECISION,
     &              psailx_g(1),jxp,MPI_DOUBLE_PRECISION,
     &                              0,MPI_COMM_WORLD,ierr)
      call MPI_GATHER(psa01(1),jxp,MPI_DOUBLE_PRECISION,
     &              psa01_g(1),jxp,MPI_DOUBLE_PRECISION,
     &                             0,MPI_COMM_WORLD,ierr)
      if(myid.eq.0) then
        do k=1,kx
        do j=1,mjx-1
       tqadv=tqadv-dtmin*3.e4*dsigma(k)*dx*
     &    (   (vaix_g(k,j+1)+vaix_g(k,j))*
     &      (qvailx_g(k,j)/psailx_g(j))/(msfx_io(ilx,j)*msfx_io(ilx,j))
     &    -   (va01_g(k,j+1)+va01_g(k,j))*
     &      (qva01_g(k,j)/psa01_g(j))/(msfx_io(1,j)*msfx_io(1,j))  )/g
        enddo
        enddo
      endif
      call MPI_BCAST(tqadv,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#else
      do k=1,kx
      do j=1,jlx
       tqadv=tqadv-dtmin*3.e4*dsigma(k)*dx*((va(ix,k,j+1)+va(ix,k,j))*
     1            (qva(ilx,k,j)/psa(ilx,j))/(msfx(ilx,j)*msfx(ilx,j))-
     2            (va(1,k,j+1)+va(1,k,j))*(qva(1,k,j)/psa(1,j))/
     3            (msfx(1,j)*msfx(1,j)))/g
      enddo
      enddo
#endif
c
c-----advection of cloud water and rainwater through lateral boundaries:
c
c.....advection through east-west boundaries:
c
#ifdef MPP1
      if (myid .eq. nproc-1) then
      do k=1,kx
      do i=1,ilx
         worka(i,k)=(ua(i+1,k,jendl)+ua(i,k,jendl))
     &             *(qca(i,k,jendx)/psa(i,jendx))
     &             /(msfx(i,jendx)*msfx(i,jendx))
      enddo
      enddo
      endif
      if (myid .eq.0) then
      do k=1,kx
      do i=1,ilx
         workb(i,k)=(ua(i+1,k,1)+ua(i,k,1))*(qca(i,k,1)/psa(i,1))
     &             /(msfx(i,1)*msfx(i,1))
      enddo
      enddo
      endif
      call MPI_BCAST(worka,ilx*kx,MPI_DOUBLE_PRECISION,nproc-1
     &              ,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(workb,ilx*kx,MPI_DOUBLE_PRECISION,0
     &              ,MPI_COMM_WORLD,ierr)
#else
      do k=1,kx
      do i=1,ilx
         worka(i,k)=(ua(i+1,k,jx)+ua(i,k,jx))*(qca(i,k,jlx)/psa(i,jlx))
     &             /(msfx(i,jlx)*msfx(i,jlx))
         workb(i,k)=(ua(i+1,k,1)+ua(i,k,1))*(qca(i,k,1)/psa(i,1))
     &             /(msfx(i,1)*msfx(i,1))
      enddo
      enddo
#endif
      do 50 k=1,kx
      do 50 i=1,ilx
       tqadv=tqadv-dtmin*3.e4*dsigma(k)*dx*(worka(i,k)-workb(i,k))/g
50    continue
c
c....advection through north-south boundaries:
c
#ifdef MPP1
      do j=1,jendl
        do k=1,kx
          qcailx(k,j)=qca(ilx,k,j)
          qca01(k,j) =qca(1  ,k,j)
        enddo
      enddo
      call MPI_GATHER(qcailx(1,1),kx*jxp,MPI_DOUBLE_PRECISION,
     &              qcailx_g(1,1),kx*jxp,MPI_DOUBLE_PRECISION,
     &                                0,MPI_COMM_WORLD,ierr)
      call MPI_GATHER(qca01(1,1),kx*jxp,MPI_DOUBLE_PRECISION,
     &              qca01_g(1,1),kx*jxp,MPI_DOUBLE_PRECISION,
     &                                0,MPI_COMM_WORLD,ierr)
      if(myid.eq.0) then
        do k=1,kx
        do j=1,mjx-1
       tqadv=tqadv-dtmin*3.e4*dsigma(k)*dx*
     &    (   (vaix_g(k,j+1)+vaix_g(k,j))*
     &      (qcailx_g(k,j)/psailx_g(j))/(msfx_io(ilx,j)*msfx_io(ilx,j))
     &    -   (va01_g(k,j+1)+va01_g(k,j))*
     &      (qca01_g(k,j)/psa01_g(j))/(msfx_io(1,j)*msfx_io(1,j))  )/g
        enddo
        enddo
      endif
      call MPI_BCAST(tqadv,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
#else
      do k=1,kx
      do j=1,jlx
       tqadv=tqadv-dtmin*3.e4*dsigma(k)*dx*((va(ix,k,j+1)+va(ix,k,j))*
     1            (qca(ilx,k,j)/psa(ilx,j))/(msfx(ilx,j)*msfx(ilx,j))-
     2            (va(1,k,j+1)+va(1,k,j))*(qca(1,k,j)/psa(1,j))/
     3            (msfx(1,j)*msfx(1,j)))/g
      enddo
      enddo
#endif
c
      return
      end

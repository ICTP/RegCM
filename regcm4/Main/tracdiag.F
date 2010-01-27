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

      subroutine tracdiag(xkc)

      implicit none
#     include "parameter.inc"
#     include "include/param1.cb"
#     include "include/diagnosis.cb"
#     include "include/param3.cb"
#     include "include/main.cb"
#     include "include/mainchem.cb"

#ifdef MPP1
#     include "include/io.cb"
#     include "mpiregcm.inc"
      integer ierr
      real(kind=8)  xkc(ix,kx,jxp)
      real(kind=8)  vailx(kx,jxp),vailx_g(kx,mjx),va02(kx,jxp)
      real(kind=8)  va02_g(kx,mjx)
      real(kind=8)  chiailx(kx,ntr,jxp),chiailx_g(kx,ntr,mjx)
      real(kind=8)  chiailx1(kx,ntr,jxp),chiailx1_g(kx,ntr,mjx)
      real(kind=8)  psailx(jxp),psailx_g(mjx),psailx1(jxp)
      real(kind=8)  psailx1_g(mjx)
      real(kind=8)  chia01(kx,ntr,jxp),chia01_g(kx,ntr,mjx)
      real(kind=8)  chia02(kx,ntr,jxp),chia02_g(kx,ntr,mjx)
      real(kind=8)  psa01(jxp),psa01_g(mjx),psa02(jxp),psa02_g(mjx)
      real(kind=8)  xkcilx1(kx,jxp),xkcilx1_g(kx,mjx)
      real(kind=8)  xkc02(kx,jxp),xkc02_g(kx,mjx)
#else
      real(kind=8)  xkc(ix,kx,jx)
#endif
      real(kind=8)  worka(ilx,kx,ntr),workb(ilx,kx,ntr)
      real(kind=8)  fact1,fact2,uavg2,fx2,uavg1,fx1,vavg2,vavg1
!      real(kind=8)  chixp1,chixm1,chiyp1,chiym1,chi00,chidx,chidy
!      real(kind=8)  chixp2,chixm2,chiyp2,chiym2
      real(kind=8)  chid1,chid2

      integer i,j,k,n
cccccccccccccccccccccccccccccccccccccccccccccccc
        
C-------------------------
C 1  ADVECTION budgets 
C-----------------------

c-----advection of tracer through lateral boundaries
     
c
c.....advection through east-west boundaries:
c
c the 'relaxed' upstream scheme
      fact1=0.6
      fact2=1-fact1

c   inflow/outflow
      do n=1,ntr
        do k=1,kx
        do i=1+1,ilx-1
#ifdef MPP1
          if (myid .eq. nproc-1) then
          uavg2=0.5*(ua(i+1,k,jendx)+ua(i,k,jendx))
          if(uavg2.lt.0.) then
           worka(i,k,n)=-uavg2*(fact1*chia(i,k,jendx,n)/psa(i,jendx)/
     &          (msfx(i,jendx)*msfx(i,jendx))+fact2*
     &          chia(i,k,jendm,n)/psa(i,jendm)/
     &          (msfx(i,jendm)*msfx(i,jendm)))
          else
           worka(i,k,n)=-uavg2*(fact1*chia(i,k,jendm,n)/psa(i,jendm)/
     &          (msfx(i,jendm)*msfx(i,jendm))+fact2*
     &          chia(i,k,jendx,n)/psa(i,jendx)/
     &          (msfx(i,jendx)*msfx(i,jendx)))
          end if
          end if
#else
          uavg2=0.5*(ua(i+1,k,jx-1)+ua(i,k,jx-1))
          if(uavg2.lt.0.) then
           worka(i,k,n)=-uavg2*(fact1*chia(i,k,jx-1,n)/psa(i,jx-1)/
     1          (msfx(i,jx-1)*msfx(i,jx-1))+fact2*
     2          chia(i,k,jlx-1,n)/psa(i,jlx-1)/
     3          (msfx(i,jlx-1)*msfx(i,jlx-1)))
          else
           worka(i,k,n)=-uavg2*(fact1*chia(i,k,jlx-1,n)/psa(i,jlx-1)/
     1          (msfx(i,jlx-1)*msfx(i,jlx-1))+fact2*
     2          chia(i,k,jx-1,n)/psa(i,jx-1)/
     3          (msfx(i,jx-1)*msfx(i,jx-1)))
          end if
#endif

#ifdef MPP1
          if (myid .eq. 0) then
#endif
          uavg1=0.5*(ua(i+1,k,1+1)+ua(i,k,1+1))
          if(uavg1.gt.0.) then
           workb(i,k,n)=-uavg1*(fact1*chia(i,k,1,n)/psa(i,1)/
     1          (msfx(i,1)*msfx(i,1))+fact2*
     2          chia(i,k,1+1,n)/psa(i,1+1)/
     3          (msfx(i,1+1)*msfx(i,1+1)))
          else
           workb(i,k,n)=-uavg1*(fact1*chia(i,k,1+1,n)/psa(i,1+1)/
     1          (msfx(i,1+1)*msfx(i,1+1))+fact2*
     2          chia(i,k,1,n)/psa(i,1)/
     3          (msfx(i,1)*msfx(i,1)))
          end if
#ifdef MPP1
          end if
#endif

        end do
        end do
      end do
#ifdef MPP1
      call MPI_BCAST(worka,ilx*kx*ntr,MPI_DOUBLE_PRECISION,nproc-1
     &              ,MPI_COMM_WORLD,ierr)
      do j=1,jendl
        do k=1,kx
          vailx(k,j) = va(ilx,k,j)
          va02(k,j)  = va(2  ,k,j)
          xkcilx1(k,j) = xkc(ix-2,k,j)
          xkc02(k,j)   = xkc(2   ,k,j)
          do n=1,ntr
            chiailx(k,n,j) = chia(ilx ,k,j,n)
            chiailx1(k,n,j)= chia(ix-2,k,j,n)
            chia01(k,n,j)  = chia(1   ,k,j,n)
            chia02(k,n,j)  = chia(2   ,k,j,n)
          end do
        end do
        psailx(j)  = psa(ilx ,j)
        psailx1(j) = psa(ix-2,j)
        psa01(j)  = psa(1   ,j)
        psa02(j)  = psa(2   ,j)
      end do
      call MPI_GATHER(vailx(1,1),kx*jxp,MPI_DOUBLE_PRECISION,
     &              vailx_g(1,1),kx*jxp,MPI_DOUBLE_PRECISION,
     &                                0,MPI_COMM_WORLD,ierr)
      call MPI_GATHER(va02(1,1),kx*jxp,MPI_DOUBLE_PRECISION,
     &              va02_g(1,1),kx*jxp,MPI_DOUBLE_PRECISION,
     &                               0,MPI_COMM_WORLD,ierr)
      call MPI_GATHER(xkcilx1(1,1),kx*jxp,MPI_DOUBLE_PRECISION,
     &              xkcilx1_g(1,1),kx*jxp,MPI_DOUBLE_PRECISION,
     &                                  0,MPI_COMM_WORLD,ierr)
      call MPI_GATHER(xkc02(1,1),kx*jxp,MPI_DOUBLE_PRECISION,
     &              xkc02_g(1,1),kx*jxp,MPI_DOUBLE_PRECISION,
     &                                0,MPI_COMM_WORLD,ierr)
      call MPI_GATHER(chiailx(1,1,1),kx*ntr*jxp,MPI_DOUBLE_PRECISION,
     &              chiailx_g(1,1,1),kx*ntr*jxp,MPI_DOUBLE_PRECISION,
     &                                        0,MPI_COMM_WORLD,ierr)
      call MPI_GATHER(chiailx1(1,1,1),kx*ntr*jxp,MPI_DOUBLE_PRECISION,
     &              chiailx1_g(1,1,1),kx*ntr*jxp,MPI_DOUBLE_PRECISION,
     &                                         0,MPI_COMM_WORLD,ierr)
      call MPI_GATHER(chia01(1,1,1),kx*ntr*jxp,MPI_DOUBLE_PRECISION,
     &              chia01_g(1,1,1),kx*ntr*jxp,MPI_DOUBLE_PRECISION,
     &                                       0,MPI_COMM_WORLD,ierr)
      call MPI_GATHER(chia02(1,1,1),kx*ntr*jxp,MPI_DOUBLE_PRECISION,
     &              chia02_g(1,1,1),kx*ntr*jxp,MPI_DOUBLE_PRECISION,
     &                                       0,MPI_COMM_WORLD,ierr)
      call MPI_GATHER(psailx(1),jxp,MPI_DOUBLE_PRECISION,
     &              psailx_g(1),jxp,MPI_DOUBLE_PRECISION,
     &                           0,MPI_COMM_WORLD,ierr)
      call MPI_GATHER(psailx1(1),jxp,MPI_DOUBLE_PRECISION,
     &              psailx1_g(1),jxp,MPI_DOUBLE_PRECISION,
     &                            0,MPI_COMM_WORLD,ierr)
      call MPI_GATHER(psa01(1),jxp,MPI_DOUBLE_PRECISION,
     &              psa01_g(1),jxp,MPI_DOUBLE_PRECISION,
     &                           0,MPI_COMM_WORLD,ierr)
      call MPI_GATHER(psa02(1),jxp,MPI_DOUBLE_PRECISION,
     &              psa02_g(1),jxp,MPI_DOUBLE_PRECISION,
     &                           0,MPI_COMM_WORLD,ierr)
      if(myid.eq.0) then
#endif
      do 44 n=1,ntr
        do 42 k=1,kx
        do 42 i=1+1,ilx-1
         tchiad(n)=tchiad(n)+dtmin*6.e4*dsigma(k)*dx*
     1                  (worka(i,k,n)-workb(i,k,n))/g
42      continue
c.....
c.....advection through north-south boundaries:
c
#ifdef MPP1
        do k=1,kx
        do j=1+1,mjx-2
chy inflow/outflow
          vavg2=0.5*(vailx_g(k,j+1)+vailx_g(k,j))
          if(vavg2.lt.0.) then
           fx2=-vavg2*(fact1*chiailx_g(k,n,j)/psailx_g(j)/
     1          (msfx_io(ilx,j)*msfx_io(ilx,j))+fact2*
     2          chiailx1_g(k,n,j)/psailx1_g(j)/
     3          (msfx_io(ix-2,j)*msfx_io(ix-2,j)))
          else
           fx2=-vavg2*(fact1*chiailx1_g(k,n,j)/psailx1_g(j)/
     1          (msfx_io(ix-2,j)*msfx_io(ix-2,j))+fact2*
     2          chiailx_g(k,n,j)/psailx_g(j)/
     3          (msfx_io(ilx,j)*msfx_io(ilx,j)))
          end if

          vavg1=0.5*(va02_g(k,j+1)+va02_g(k,j))
          if(vavg1.gt.0.) then
           fx1=-vavg1*(fact1*chia01_g(k,n,j)/psa01_g(j)/
     1          (msfx_io(1,j)*msfx_io(1,j))+fact2*
     2          chia02_g(k,n,j)/psa02_g(j)/
     3          (msfx_io(2,j)*msfx_io(2,j)))
          else
           fx1=-vavg1*(fact1*chia02_g(k,n,j)/psa02_g(j)/
     1          (msfx_io(2,j)*msfx_io(2,j))+fact2*
     2          chia01_g(k,n,j)/psa01_g(j)/
     3          (msfx_io(1,j)*msfx_io(1,j)))
          end if

          tchiad(n)=tchiad(n)+dtmin*6.e4*dsigma(k)*dx*
     1          (fx2-fx1)/g

        end do
        end do
#else
        do k=1,kx
        do j=1+1,jlx-1
chy inflow/outflow
          vavg2=0.5*(va(ix-1,k,j+1)+va(ix-1,k,j))
          if(vavg2.lt.0.) then
           fx2=-vavg2*(fact1*chia(ix-1,k,j,n)/psa(ix-1,j)/
     1          (msfx(ix-1,j)*msfx(ix-1,j))+fact2*
     2          chia(ilx-1,k,j,n)/psa(ilx-1,j)/
     3          (msfx(ilx-1,j)*msfx(ilx-1,j)))
          else
           fx2=-vavg2*(fact1*chia(ilx-1,k,j,n)/psa(ilx-1,j)/
     1          (msfx(ilx-1,j)*msfx(ilx-1,j))+fact2*
     2          chia(ix-1,k,j,n)/psa(ix-1,j)/
     3          (msfx(ix-1,j)*msfx(ix-1,j)))
          end if

          vavg1=0.5*(va(1+1,k,j+1)+va(1+1,k,j))
          if(vavg1.gt.0.) then
           fx1=-vavg1*(fact1*chia(1,k,j,n)/psa(1,j)/
     1          (msfx(1,j)*msfx(1,j))+fact2*
     2          chia(1+1,k,j,n)/psa(1+1,j)/
     3          (msfx(1+1,j)*msfx(1+1,j)))
          else
           fx1=-vavg1*(fact1*chia(1+1,k,j,n)/psa(1+1,j)/
     1          (msfx(1+1,j)*msfx(1+1,j))+fact2*
     2          chia(1,k,j,n)/psa(1,j)/
     3          (msfx(1,j)*msfx(1,j)))
          end if

          tchiad(n)=tchiad(n)+dtmin*6.e4*dsigma(k)*dx*
     1          (fx2-fx1)/g

        end do
        end do
#endif
44    continue
#ifdef MPP1
      end if
      call MPI_BCAST(tchiad,ntr,MPI_DOUBLE_PRECISION,
     &                        0,MPI_COMM_WORLD,ierr)
#endif

CC      do 48 n=1,ntr
CCchy net advection at x direction and y direction at every inner grid box
CC
CC        do 45 k=1,kx
CC#ifdef MPP1
CC        do 45 j=jbegin,jendm
CC#else
CC        do 45 j=2,jlx-1
CC#endif
CC        do 45 i=2,ilx-1
CC
CC          uavg2=0.5*(ua(i+1,k,j+1)+ua(i,k,j+1))
CC          uavg1=0.5*(ua(i+1,k,j)+ua(i,k,j))
CC          if(uavg2.lt.0.) then
CC           fx2=-uavg2*(fact1*chia(i,k,j+1,n)/psa(i,j+1)/
CC     1          (msfx(i,j+1)*msfx(i,j+1))+fact2*
CC     2          chia(i,k,j,n)/psa(i,j)/
CC     3          (msfx(i,j)*msfx(i,j)))
CC          else
CC           fx2=-uavg2*(fact1*chia(i,k,j,n)/psa(i,j)/
CC     1          (msfx(i,j)*msfx(i,j))+fact2*
CC     2          chia(i,k,j+1,n)/psa(i,j+1)/
CC     3          (msfx(i,j+1)*msfx(i,j+1)))
CC          end if
CC
CC          if(uavg1.gt.0.) then
CC           fx1=-uavg1*(fact1*chia(i,k,j-1,n)/psa(i,j-1)/
CC     1          (msfx(i,j-1)*msfx(i,j-1))+fact2*
CC     2          chia(i,k,j,n)/psa(i,j)/
CC     3          (msfx(i,j)*msfx(i,j)))
CC          else
CC           fx1=-uavg1*(fact1*chia(i,k,j,n)/psa(i,j)/
CC     1          (msfx(i,j)*msfx(i,j))+fact2*
CC     2          chia(i,k,j-1,n)/psa(i,j-1)/
CC     3          (msfx(i,j-1)*msfx(i,j-1)))
CC          end if
CC
CC          tchifx(i,j,k,n)=tchifx(i,j,k,n)+dtmin*6.e4*dsigma(k)*dx*
CC     1                   (fx2-fx1)/g
CC
CCccc advection at y direction
CC
CC          vavg2=0.5*(va(i+1,k,j+1)+va(i+1,k,j))
CC          vavg1=0.5*(va(i,k,j+1)+va(i,k,j))
CC          if(vavg2.lt.0.) then
CC           fx2=-vavg2*(fact1*chia(i+1,k,j,n)/psa(i+1,j)/
CC     1          (msfx(i+1,j)*msfx(i+1,j))+fact2*
CC     2          chia(i,k,j,n)/psa(i,j)/
CC     3          (msfx(i,j)*msfx(i,j)))
CC          else
CC           fx2=-vavg2*(fact1*chia(i,k,j,n)/psa(i,j)/
CC     1          (msfx(i,j)*msfx(i,j))+fact2*
CC     2          chia(i+1,k,j,n)/psa(i+1,j)/
CC     3          (msfx(i+1,j)*msfx(i+1,j)))
CC          end if
CC
CC          if(vavg1.gt.0.) then
CC           fx1=-vavg1*(fact1*chia(i-1,k,j,n)/psa(i-1,j)/
CC     1          (msfx(i-1,j)*msfx(i-1,j))+fact2*
CC     2          chia(i,k,j,n)/psa(i,j)/
CC     3          (msfx(i,j)*msfx(i,j)))
CC          else
CC           fx1=-vavg1*(fact1*chia(i,k,j,n)/psa(i,j)/
CC     1          (msfx(i,j)*msfx(i,j))+fact2*
CC     2          chia(i-1,k,j,n)/psa(i-1,j)/
CC     3          (msfx(i-1,j)*msfx(i-1,j)))
CC          end if
CC
CC          tchify(i,j,k,n)=tchify(i,j,k,n)+dtmin*6.e4*dsigma(k)*dx*
CC     1                   (fx2-fx1)/g
CC
CC45      continue
CC
CC48    continue
CC
CCC
CCC-----------------------------
CCC 2  DIFFUSION Budgets
CCC-----------------------------
CCC
CC      do 1979 n=1,ntr
CC      do 1979 k=1,kx
CC#ifdef MPP1
CC      do 1979 j=jbegin,jendm
CC#else
CC      do 1979 j=2,jlxm
CC#endif
CC      do 1979 i=2,ilxm
CCc
CCc west/east and south/north boundary
CCc
CC        chixp1=chia(i,k,j+1,n)/psa(i,j+1)
CC        chixm1=chia(i,k,j-1,n)/psa(i,j-1)
CC        chiyp1=chia(i+1,k,j,n)/psa(i+1,j)
CC        chiym1=chia(i-1,k,j,n)/psa(i-1,j)
CC        chi00 =chia(i,k,j,n)/psa(i,j)
CC
CC#ifdef MPP1
CC        if((myid .eq.0 .and.j.eq.2) .or.
CC     &     (myid .eq. nproc-1 .and. j.eq.jendm)
CC     &    .or. i.eq.2 .or. i.eq.ilxm) then
CC#else
CC        if(j.eq.2 .or. j.eq.jlxm .or. i.eq.2 .or. i.eq.ilxm) then
CC#endif
CC
CC         chidx=xkc(i,k,j)*(chixp1+chixm1-2.*chi00)*psa(i,j)
CC         chidy=xkc(i,k,j)*(chiyp1+chiym1-2.*chi00)*psa(i,j)
CC
CC        else
CCc
CCc interior points --- fourth-order form
CCc
CC        chixp2=chia(i,k,j+2,n)/psa(i,j+2)
CC        chixm2=chia(i,k,j-2,n)/psa(i,j-2)
CC        chiyp2=chia(i+2,k,j,n)/psa(i+2,j)
CC        chiym2=chia(i-2,k,j,n)/psa(i-2,j)
CC
CC         chidx=-xkc(i,k,j)*
CC     1    (chixp2+chixm2- 4.*(chixp1+chixm1) + 6.*chi00)*psa(i,j)
CC
CC         chidy=-xkc(i,k,j)*
CC     1    (chixp2+chixm2- 4.*(chixp1+chixm1) + 6.*chi00)*psa(i,j)
CC
CC         chidy=-xkc(i,k,j)*
CC     1    (chiyp2+chiym2- 4.*(chiyp1+chiym1) + 6.*chi00)*psa(i,j)
CC
CC        end if
CC
CC        tchidx(i,j,k,n)=tchidx(i,j,k,n)+dtmin*6.e4*dsigma(k)
CC     1                  *chidx/g
CC        tchidy(i,j,k,n)=tchidy(i,j,k,n)+dtmin*6.e4*dsigma(k)
CC     1                  *chidy/g
CC
CC1979  continue
c
c..... diffusion through east-west boundaries:
c 
c
      do n=1,ntr
        do k=1,kx
        do i=1+1,ilx-1
#ifdef MPP1
          if (myid .eq. nproc-1) then
          worka(i,k,n)=xkc(i,k,jendm)*psa(i,jendm)*
     &          (chia(i,k,jendm,n)/psa(i,jendm)-
     &           chia(i,k,jendx,n)/psa(i,jendx))
          end if
#else
          worka(i,k,n)=xkc(i,k,jlx-1)*psa(i,jlx-1)*
     1          (chia(i,k,jlx-1,n)/psa(i,jlx-1)-
     2           chia(i,k,jlx,n)/psa(i,jlx))
#endif
#ifdef MPP1
          if (myid .eq. 0) then
#endif
          workb(i,k,n)=xkc(i,k,2)*psa(i,2)*
     1          (chia(i,k,2,n)/psa(i,2)-chia(i,k,1,n)/psa(i,1))
#ifdef MPP1
          end if
#endif
        end do
        end do
      end do
#ifdef MPP1
      call MPI_BCAST(worka,ilx*kx*ntr,MPI_DOUBLE_PRECISION,nproc-1
     &              ,MPI_COMM_WORLD,ierr)
      if(myid.eq.0) then
#endif
      do 1300 n=1,ntr
        do k=1,kx
        do i=1+1,ilx-1
          tchitb(n)=tchitb(n)-dtmin*6.e4*dsigma(k)
     &                             *(workb(i,k,n)+worka(i,k,n))/g
        end do
        end do
     
c..... diffusion through north-south boundaries:

#ifdef MPP1
        do k=1,kx
        do j=1+1,mjx-2
          chid1=xkcilx1_g(k,j)*psailx1_g(j)*
     &          (chiailx1_g(k,n,j)/psailx1_g(j)-
     &           chiailx_g(k,n,j)/psailx_g(j))
          chid2=xkc02_g(k,j)*psa02_g(j)*
     &          (chia02_g(k,n,j)/psa02_g(j)-
     &           chia01_g(k,n,j)/psa01_g(j))
          tchitb(n)=tchitb(n)-dtmin*6.e4*dsigma(k)*(chid2+chid1)/g
        end do
        end do
#else
        do k=1,kx
        do j=1+1,jlx-1
          chid1=xkc(ilx-1,k,j)*psa(ilx-1,j)*
     1          (chia(ilx-1,k,j,n)/psa(ilx-1,j)-
     2           chia(ilx,k,j,n)/psa(ilx,j))
          chid2=xkc(2,k,j)*psa(2,j)*
     1          (chia(2,k,j,n)/psa(2,j)-chia(1,k,j,n)/psa(1,j))
          tchitb(n)=tchitb(n)-dtmin*6.e4*dsigma(k)*(chid2+chid1)/g
        end do
        end do
#endif
1300  continue
#ifdef MPP1
      end if
      call MPI_BCAST(tchitb,ntr,MPI_DOUBLE_PRECISION,
     &                        0,MPI_COMM_WORLD,ierr)
#endif
      return
      end subroutine tracdiag

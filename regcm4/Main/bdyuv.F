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

      subroutine bdyuv(ib,dtb)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine sets the boundary values of u and v according   c
c     to the boundary conditions specified.                           c
c                                                                     c
c     ua, va, and psa : variables needed                              c
c                                                                     c
c     ib = 0 : fixed                                                  c
c        = 1 : relaxation, linear technique                           c
c        = 2 : time dependent                                         c
c        = 3 : time dependent and inflow/outflow dependent            c
c        = 4 : sponge                                                 c
c        = 5 : relaxation, exponential technique                      c
c                                                                     c
c     dtb    : elapsed time from the initial boundary values.         c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
#     include "parameter.inc"

#     include "include/bdycod.cb"
#     include "include/cvaria.cb"
#     include "include/main.cb"

#ifdef MPP1
      include 'mpif.h'
      integer status(MPI_STATUS_SIZE)
      integer ierr
      real(kind=8)  var1snd(kx,8),var1rcv(kx,8)
      COMMON /TMP_var1/ var1snd,var1rcv
#endif

      integer ib
      real(kind=8)  dtb
      integer i,j,k
c
c----------------------------------------------------------------------
c-----compute the p* at dot points:
c
C=======================================================================
#ifdef MPP1
      call MPI_SENDRECV(psa(1,jxp),ix,MPI_DOUBLE_PRECISION,ieast,1,
     &                  psa(1,0),ix,MPI_DOUBLE_PRECISION,iwest,1,
     &                  MPI_COMM_WORLD,status,ierr)
#endif
C=======================================================================
c
c-----interior points:
c
#ifdef MPP1
      do j=jbegin,jendx
#else
      do j=2,jlx
#endif
        do i=2,ilx
          pdota(i,j)=0.25*(psa(i,j)+psa(i-1,j)+psa(i,j-1)+psa(i-1,j-1))
        end do
      end do
c
c-----east and west boundaries:
c
      do i=2,ilx
#ifdef MPP1
      if (myid .eq. 0) then
#endif
        pdota(i,1)=0.5*(psa(i,1)+psa(i-1,1))
#ifdef MPP1
      end if
      if (myid .eq. nproc-1) then
        pdota(i,jendl)=0.5*(psa(i,jendx)+psa(i-1,jendx))
      end if
#else
        pdota(i,jx)=0.5*(psa(i,jlx)+psa(i-1,jlx))
#endif
      end do
c
c-----north and south boundaries:
c
#ifdef MPP1
      do j=jbegin,jendx
#else
      do j=2,jlx
#endif
        pdota(1,j)=0.5*(psa(1,j)+psa(1,j-1))
        pdota(ix,j)=0.5*(psa(ilx,j)+psa(ilx,j-1))
      end do
c
c-----corner points:
c
#ifdef MPP1
      if (myid .eq. 0) then
#endif
      pdota(1,1)=psa(1,1)
      pdota(ix,1)=psa(ilx,1)
#ifdef MPP1
      end if
      if (myid .eq. nproc-1) then
        pdota(1,jendl)=psa(1,jendx)
        pdota(ix,jendl)=psa(ilx,jendx)
      end if
#else
      pdota(1,jx)=psa(1,jlx)
      pdota(ix,jx)=psa(ilx,jlx)
#endif
C=======================================================================
c
c-----interior silces:
c
      do 30 k=1,kx
c
c.....for j = 2 and j = jlx :
c
       do 10 i=2,ilx
#ifdef MPP1
        if (myid .eq. 0) then
#endif
        uj2(i,k)=ua(i,k,2)/pdota(i,2)
        vj2(i,k)=va(i,k,2)/pdota(i,2)
#ifdef MPP1
        end if
        if (myid .eq. nproc-1 ) then
          ujlx(i,k)=ua(i,k,jendx)/pdota(i,jendx)
          vjlx(i,k)=va(i,k,jendx)/pdota(i,jendx)
        end if
#else
        ujlx(i,k)=ua(i,k,jlx)/pdota(i,jlx)
        vjlx(i,k)=va(i,k,jlx)/pdota(i,jlx)
#endif
10     continue
c
c.....for i = 2 and i = ilx :
c
#ifdef MPP1
       do 20 j=jbegin,jendx
#else
       do 20 j=2,jlx
#endif
        ui2(k,j)=ua(2,k,j)/pdota(2,j)
        vi2(k,j)=va(2,k,j)/pdota(2,j)
        uilx(k,j)=ua(ilx,k,j)/pdota(ilx,j)
        vilx(k,j)=va(ilx,k,j)/pdota(ilx,j)
20     continue
c
30    continue
c
c----------------------------------------------------------------------
c-----boundary silces:
c
      if (ib .eq. 0) then
c
c-----fixed boundary conditions:
c
       do 60 k=1,kx
c
c.....west (j = 1) and east (j = jx) boundaries:
c
        do 40 i=1,ix
#ifdef MPP1
         if (myid .eq. 0) then
#endif
         uj1(i,k)=uwb(i,k,1)/pdota(i,1)
         vj1(i,k)=vwb(i,k,1)/pdota(i,1)
#ifdef MPP1
         end if
         if (myid .eq. nproc-1) then
           ujl(i,k)=ueb(i,k,1)/pdota(i,jendl)
           vjl(i,k)=veb(i,k,1)/pdota(i,jendl)
         end if
#else
         ujl(i,k)=ueb(i,k,1)/pdota(i,jx)
         vjl(i,k)=veb(i,k,1)/pdota(i,jx)
#endif
40      continue
c
c.....south (i = 1) and north (i = ix) boundaries:
c
#ifdef MPP1
        do 60 j=1,jendl
#else
        do 60 j=1,jx
#endif
         ui1(k,j)=usb(1,k,j)/pdota(1,j)
         vi1(k,j)=vsb(1,k,j)/pdota(1,j)
         uil(k,j)=unb(1,k,j)/pdota(ix,j)
         vil(k,j)=vnb(1,k,j)/pdota(ix,j)
c
60     continue
       go to 110
c
      end if
c
c-----time-dependent boundary conditions:
c
      do 100 k=1,kx
c
c.....west (j = 1) and east (j = jx) boundaries:
c
       do 80 i=1,ix
#ifdef MPP1
        if (myid .eq. 0) then
#endif
        uj1(i,k)=(uwb(i,k,1)+dtb*uwbt(i,k,1))/pdota(i,1)
        vj1(i,k)=(vwb(i,k,1)+dtb*vwbt(i,k,1))/pdota(i,1)
#ifdef MPP1
        end if
        if (myid .eq. nproc-1) then
          ujl(i,k)=(ueb(i,k,1)+dtb*uebt(i,k,1))/pdota(i,jendl)
          vjl(i,k)=(veb(i,k,1)+dtb*vebt(i,k,1))/pdota(i,jendl)
        end if
#else
        ujl(i,k)=(ueb(i,k,1)+dtb*uebt(i,k,1))/pdota(i,jx)
        vjl(i,k)=(veb(i,k,1)+dtb*vebt(i,k,1))/pdota(i,jx)
#endif
80     continue
c
c.....south (i = 1) and north (i = ix) boundaries:
c
#ifdef MPP1
       do  90 j=1,jendl
#else
       do  90 j=1,jx
#endif
        ui1(k,j)=(usb(1,k,j)+dtb*usbt(1,k,j))/pdota(1,j)
        vi1(k,j)=(vsb(1,k,j)+dtb*vsbt(1,k,j))/pdota(1,j)
        uil(k,j)=(unb(1,k,j)+dtb*unbt(1,k,j))/pdota(ix,j)
        vil(k,j)=(vnb(1,k,j)+dtb*vnbt(1,k,j))/pdota(ix,j)
90     continue
c
100   continue
c
c-----fill up the interior silces:
c
110   continue
      do 120 k=1,kx
#ifdef MPP1
       if (myid .eq. 0) then
#endif
       uj2(1,k) = ui1(k,2)
       uj2(ix,k) = uil(k,2)
       ui2(k,1) = uj1(2,k)
       uilx(k,1) = uj1(ilx,k)
       vj2(1,k) = vi1(k,2)
       vj2(ix,k) = vil(k,2)
       vi2(k,1) = vj1(2,k)
       vilx(k,1) = vj1(ilx,k)
#ifdef MPP1
       end if
       if (myid .eq. nproc-1) then
       ujlx(1,k) = ui1(k,jendx)
       ujlx(ix,k) = uil(k,jendx)
       ui2(k,jendl) = ujl(2,k)
       uilx(k,jendl) = ujl(ilx,k)
       vjlx(1,k) = vi1(k,jendx)
       vjlx(ix,k) = vil(k,jendx)
       vi2(k,jendl) = vjl(2,k)
       vilx(k,jendl) = vjl(ilx,k)
       end if
#else
       ujlx(1,k) = ui1(k,jlx)
       ujlx(ix,k) = uil(k,jlx)
       ui2(k,jx) = ujl(2,k)
       uilx(k,jx) = ujl(ilx,k)
       vjlx(1,k) = vi1(k,jlx)
       vjlx(ix,k) = vil(k,jlx)
       vi2(k,jx) = vjl(2,k)
       vilx(k,jx) = vjl(ilx,k)
#endif
120   continue
#ifdef MPP1
      if(myid.ne.nproc-1) then
        do k=1,kx
          var1snd(k,1) = ui1(k,jxp)
          var1snd(k,2) = vi1(k,jxp)
          var1snd(k,3) = ui2(k,jxp)
          var1snd(k,4) = vi2(k,jxp)
          var1snd(k,5) =uilx(k,jxp)
          var1snd(k,6) =vilx(k,jxp)
          var1snd(k,7) = uil(k,jxp)
          var1snd(k,8) = vil(k,jxp)
        end do
      end if
      call MPI_SENDRECV(var1snd(1,1),kx*8,MPI_DOUBLE_PRECISION,
     &          ieast,1,var1rcv(1,1),kx*8,MPI_DOUBLE_PRECISION,
     &          iwest,1,MPI_COMM_WORLD,status,ierr)
      if(myid.ne.0) then
        do k=1,kx
           ui1(k,0) = var1rcv(k,1)
           vi1(k,0) = var1rcv(k,2)
           ui2(k,0) = var1rcv(k,3)
           vi2(k,0) = var1rcv(k,4)
          uilx(k,0) = var1rcv(k,5)
          vilx(k,0) = var1rcv(k,6)
           uil(k,0) = var1rcv(k,7)
           vil(k,0) = var1rcv(k,8)
        end do
      end if
C
      if(myid.ne.0) then
        do k=1,kx
          var1snd(k,1) = ui1(k,1)
          var1snd(k,2) = vi1(k,1)
          var1snd(k,3) = ui2(k,1)
          var1snd(k,4) = vi2(k,1)
          var1snd(k,5) =uilx(k,1)
          var1snd(k,6) =vilx(k,1)
          var1snd(k,7) = uil(k,1)
          var1snd(k,8) = vil(k,1)
        end do
      end if
      call MPI_SENDRECV(var1snd(1,1),kx*8,MPI_DOUBLE_PRECISION,
     &          iwest,2,var1rcv(1,1),kx*8,MPI_DOUBLE_PRECISION,
     &          ieast,2,MPI_COMM_WORLD,status,ierr)
      if(myid.ne.nproc-1) then
        do k=1,kx
           ui1(k,jxp+1) = var1rcv(k,1)
           vi1(k,jxp+1) = var1rcv(k,2)
           ui2(k,jxp+1) = var1rcv(k,3)
           vi2(k,jxp+1) = var1rcv(k,4)
          uilx(k,jxp+1) = var1rcv(k,5)
          vilx(k,jxp+1) = var1rcv(k,6)
           uil(k,jxp+1) = var1rcv(k,7)
           vil(k,jxp+1) = var1rcv(k,8)
        end do
      end if
#endif
c
      return
      end subroutine bdyuv

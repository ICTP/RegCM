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

      subroutine bdyin
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine reads in the boundary conditions.               c
c                                                                     c
c        iunit : is the unit number from which the data are read in.  c
c                                                                     c
c        xtime : is the time in minutes into the forecast.            c
c                                                                     c
c        bdytim : is the time in minutes after which the boundary     c
c                 conditions are needed.                              c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
#     include "parameter.inc"

#     include "include/iunits.cb"
#     include "include/param1.cb"
#     include "include/param2.cb"
#     include "include/param3.cb"
#     include "include/main.cb"
#     include "include/bdycod.cb"
#     include "include/bats.cb2"
#     include "include/trachem.cb"
#     include "include/date.cb"

#ifdef MPP1
#     include "include/io.cb"
#     include "mpiregcm.inc"
      integer status(MPI_STATUS_SIZE)
      integer ierr
      integer nxwb,nxeb,ndwb,ndeb,nkk
      real(kind=8)  tdum(ix,jxp)
      real(kind=8)  psdot(ix,jxp)
      real(kind=8)   sav0(ix,kx*4+2,jxp),sav_0(ix,kx*4+2,mjx)
      real(kind=8)   sav0s(ix,kx,jxp),sav_0s(ix,kx,mjx)
      COMMON /TMPSAV0/ sav0,sav_0,sav0s,sav_0s
      real(kind=4)  io2d(ix,mjx)
#else
      real(kind=8)  tdum(ix,jx)
      real(kind=8)  psdot(ix,jx)
      real(kind=4)  io2d(ix,jx)
#endif
C
      INTEGER         MDATEZ
      COMMON /DATENUM/MDATEZ(289276)
C
      integer nn,nnb,k,ierr1
      real(kind=8)  dtbdys

      real(kind=8)  tinit,tfinal

      integer i,j
      character*8 finm
C_sgi character*7 finm
c
c-----tinit  is the initial time in minutes of the boundary conditions.
c     tfinal is the final   time in minutes of the boundary conditions.
c     the boundary tendencies stored in tape are calculated as:
c         a(at tfinal) - a (at tinit) / (tfinal-tinit)*60.
c     where a is p*, p*u, p*v, p*t, and p*qv.
c
c----------------------------------------------------------------------
      if ( dabs(xtime).gt.0.0001) return
c
10    continue
      dtbdys = IBDYFRQ*60.*60.
#ifdef MPP1
      if (myid .eq. 0) then
#endif
 222  continue
      if (EHSO4) then
         do k=1,kx
#ifdef MPP1
           do j=1,jendl
#else
           do j=1,jx
#endif
             do i=1,ix
                so4(i,k,j) = SO0(i,k,j)
             end do
           end do
         end do
      end if
      mmrec=mmrec+1
      read(iutbc,rec=mmrec,iostat=ierr1) ndate1
      if (ierr1.ne.0) then
        close (iutbc)
        iutbc = iutbc + 1
        write(finm,101) iutbc
 101    format('fort.',I3)
C_sgi 101    format('fort.',I2)
#ifdef MPP1
        open(iutbc,file=finm,form='unformatted',status='old'
     &            ,access='direct',recl=ix*mjx*ibyte)
#else
        open(iutbc,file=finm,form='unformatted',status='old'
     &            ,access='direct',recl=ix*jx*ibyte)
#endif
        mmrec = 0
        print*,'CHANGING BDY UNIT NUMBER:  iutbc=',iutbc
        if (iutbc.gt.999) then
           call fatal(__FILE__,__LINE__,'BDY UNIT MAX EXCEEDED')
        end if
        go to 222
      end if
      if (ndate1.lt.MDATEZ(NNNCHK+1)) then
         if (ndate1.lt.MDATEZ(NNNCHK+1)) then
            print*,'Searching for proper date: ',ndate1,MDATEZ(NNNCHK+1)
            print*,'read in datasets at :',ndate0
            if (EHSO4) then
               if (LSMTYP.ne.'USGS') then
                  mmrec=mmrec+kx*5+2
               else
                  mmrec=mmrec+kx*5+2+13
               end if
            else
               if (LSMTYP.ne.'USGS') then
                  mmrec=mmrec+kx*4+2
               else
                  mmrec=mmrec+kx*4+2+13
               end if
            end if
            go to 222
         end if
      else if (ndate1.gt.MDATEZ(NNNCHK+1)) then
         print*,'DATE IN BC FILE EXCEEDED DATE IN RegCM'
         print*,ndate1,MDATEZ(NNNCHK+1),NNNCHK+1
         call fatal(__FILE__,__LINE__,'ICBC date')
      end if
#ifdef MPP1
c     print*,'UB1'
      do k=kx,1,-1
        mmrec=mmrec+1
        read(iutbc,rec=mmrec) ((io2d(i,j),j=1,mjx),i=1,ix)
        do j=1,mjx
        do i=1,ix
          ub1_io(i,k,j)=dble(io2d(i,j))
        end do
        end do
      end do
c     print*,'VB1'
      do k=kx,1,-1
        mmrec=mmrec+1
        read(iutbc,rec=mmrec) ((io2d(i,j),j=1,mjx),i=1,ix)
        do j=1,mjx
        do i=1,ix
          vb1_io(i,k,j)=dble(io2d(i,j))
        end do
        end do
      end do
c     print*,'TB1'
      do k=kx,1,-1
        mmrec=mmrec+1
        read(iutbc,rec=mmrec) ((io2d(i,j),j=1,mjx),i=1,ix)
        do j=1,mjx
        do i=1,ix
          tb1_io(i,k,j)=dble(io2d(i,j))
        end do
        end do
      end do
c     print*,'QB1'
      do k=kx,1,-1
        mmrec=mmrec+1
        read(iutbc,rec=mmrec) ((io2d(i,j),j=1,mjx),i=1,ix)
        do j=1,mjx
        do i=1,ix
          qb1_io(i,k,j)=dble(io2d(i,j))
        end do
        end do
      end do
c     print*,'PS1'
      mmrec=mmrec+1
      read(iutbc,rec=mmrec) ((io2d(i,j),j=1,mjx),i=1,ix)
      do j=1,mjx
      do i=1,ix
        ps1_io(i,j)=dble(io2d(i,j))
      end do
      end do
c     print*,'TS1'
      mmrec=mmrec+1
      read(iutbc,rec=mmrec) ((io2d(i,j),j=1,mjx),i=1,ix)
      do j=1,mjx
      do i=1,ix
        ts1_io(i,j)=dble(io2d(i,j))
      end do
      end do
      if (EHSO4) then
c        print*,'SO1'
         do k=kx,1,-1
           mmrec=mmrec+1
           read(iutbc,rec=mmrec) ((io2d(i,j),j=1,mjx),i=1,ix)
           do j=1,mjx
           do i=1,ix
             so1_io(i,k,j)=dble(io2d(i,j))
           end do
           end do
         end do
      end if
#else
c     print*,'UB1'
      do k=kx,1,-1
        mmrec=mmrec+1
        read(iutbc,rec=mmrec) ((io2d(i,j),j=1,jx),i=1,ix)
        do j=1,jx
        do i=1,ix
          ub1(i,k,j)=dble(io2d(i,j))
        end do
        end do
      end do
c     print*,'VB1'
      do k=kx,1,-1
        mmrec=mmrec+1
        read(iutbc,rec=mmrec) ((io2d(i,j),j=1,jx),i=1,ix)
        do j=1,jx
        do i=1,ix
          vb1(i,k,j)=dble(io2d(i,j))
        end do
        end do
      end do
c     print*,'TB1'
      do k=kx,1,-1
        mmrec=mmrec+1
        read(iutbc,rec=mmrec) ((io2d(i,j),j=1,jx),i=1,ix)
        do j=1,jx
        do i=1,ix
          tb1(i,k,j)=dble(io2d(i,j))
        end do
        end do
      end do
c     print*,'QB1'
      do k=kx,1,-1
        mmrec=mmrec+1
        read(iutbc,rec=mmrec) ((io2d(i,j),j=1,jx),i=1,ix)
        do j=1,jx
        do i=1,ix
          qb1(i,k,j)=dble(io2d(i,j))
        end do
        end do
      end do
c     print*,'PS1'
      mmrec=mmrec+1
      read(iutbc,rec=mmrec) ((io2d(i,j),j=1,jx),i=1,ix)
      do j=1,jx
      do i=1,ix
        ps1(i,j)=dble(io2d(i,j))
      end do
      end do
c     print*,'TS1'
      mmrec=mmrec+1
      read(iutbc,rec=mmrec) ((io2d(i,j),j=1,jx),i=1,ix)
      do j=1,jx
      do i=1,ix
        ts1(i,j)=dble(io2d(i,j))
      end do
      end do
      if (EHSO4) then
c     print*,'SO1'
         do k=kx,1,-1
           mmrec=mmrec+1
           read(iutbc,rec=mmrec) ((io2d(i,j),j=1,jx),i=1,ix)
           do j=1,jx
           do i=1,ix
             so1(i,k,j)=dble(io2d(i,j))
           end do
           end do
         end do
      end if
#endif
      if (LSMTYP.eq.'USGS') mmrec=mmrec+13
#ifdef MPP1
        do j=1,mjx
          do k=1,kx
            do i=1,ix
              sav_0(i,      k,j) = UB1_io(i,k,j)
              sav_0(i,kx   +k,j) = VB1_io(i,k,j)
              sav_0(i,kx*2 +k,j) = QB1_io(i,k,j)
              sav_0(i,kx*3 +k,j) = TB1_io(i,k,j)
            end do
          end do
          do i=1,ix
            sav_0(i,kx*4+1 ,j) = PS1_io(i,j)
            sav_0(i,kx*4+2 ,j) = TS1_io(i,j)
          end do
        end do
        if (EHSO4) then
           do j=1,mjx
             do k=1,kx
               do i=1,ix
                 sav_0s(i,      k,j) = SO1_io(i,k,j)
               end do
             end do
           end do
        end if
      end if
      call MPI_BCAST(ndate1,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_SCATTER(
     &         sav_0(1,1,1),ix*(kx*4+2)*jxp,MPI_DOUBLE_PRECISION,
     &          sav0(1,1,1),ix*(kx*4+2)*jxp,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
      do j=1,jendl
        do k=1,kx
          do i=1,ix
            UB1(i,k,j) = sav0(i,      k,j)
            VB1(i,k,j) = sav0(i,kx   +k,j)
            QB1(i,k,j) = sav0(i,kx*2 +k,j)
            TB1(i,k,j) = sav0(i,kx*3 +k,j)
          end do
        end do
        do i=1,ix
          PS1(i,j) = sav0(i,kx*4+1 ,j)
          TS1(i,j) = sav0(i,kx*4+2 ,j)
        end do
      end do
      if (EHSO4) then
         call MPI_SCATTER(
     &         sav_0s(1,1,1),ix*kx*jxp,MPI_DOUBLE_PRECISION,
     &          sav0s(1,1,1),ix*kx*jxp,MPI_DOUBLE_PRECISION,
     &                               0,MPI_COMM_WORLD,ierr)
         do j=1,jendl
           do k=1,kx
             do i=1,ix
               SO1(i,k,j) = sav0s(i,      k,j)
             end do
           end do
         end do
      end if
#endif
c Convert surface pressure to pstar
#ifdef MPP1
      do j=1,jendl
#else
      do j=1,jx
#endif
      do i=1,ix
        PS1(i,j) = PS1(i,j) - ptop
      end do
      end do
C=======================================================================
c
c  this routine determines p(.) from p(x) by a 4-point interpolation.
c  on the x-grid, a p(x) point outside the grid domain is assumed to
c  satisfy p(0,j)=p(1,j); p(ix,j)=p(ix-1,j); and similarly for the i's.
c
#ifdef MPP1
      call MPI_SENDRECV(PS1(1,jxp),ix,MPI_DOUBLE_PRECISION,ieast,1,
     &                  PS1(1,0),ix,MPI_DOUBLE_PRECISION,iwest,1,
     &                  MPI_COMM_WORLD,status,ierr)
#endif
#ifdef MPP1
      do j=jbegin,jendx
#else
      do j=2,jlx
#endif
        do i=2,ilx
          psdot(i,j)=0.25*(PS1(i,j)+PS1(i-1,j)+PS1(i,j-1)+PS1(i-1,j-1))
        end do
      end do
c
      do i=2,ilx
#ifdef MPP1
        if (myid .eq. 0) then
#endif
        psdot(i,1)=0.5*(PS1(i,1)+PS1(i-1,1))
#ifdef MPP1
        end if
        if (myid .eq. nproc-1) then
          psdot(i,jendl)=0.5*(PS1(i,jendx)+PS1(i-1,jendx))
        end if
#else
        psdot(i,jx)=0.5*(PS1(i,jlx)+PS1(i-1,jlx))
#endif
      end do
c
#ifdef MPP1
      do j=jbegin,jendx
#else
      do j=2,jlx
#endif
        psdot(1,j)=0.5*(PS1(1,j)+PS1(1,j-1))
        psdot(ix,j)=0.5*(PS1(ilx,j)+PS1(ilx,j-1))
      end do
c
#ifdef MPP1
      if (myid .eq. 0) then
#endif
      psdot(1,1)=PS1(1,1)
      psdot(ix,1)=PS1(ilx,1)
#ifdef MPP1
      end if
      if (myid .eq. nproc-1) then
        psdot(1,jendl)=PS1(1,jendx)
        psdot(ix,jendl)=PS1(ilx,jendx)
      end if
#else
      psdot(1,jx)=PS1(1,jlx)
      psdot(ix,jx)=PS1(ilx,jlx)
#endif
c
C=======================================================================
C Couple pressure u,v,t,q
      do k=1,kx
#ifdef MPP1
      do j=1,jendl
#else
      do j=1,jx
#endif
      do i=1,ix
        UB1(i,k,j) = UB1(i,k,j)*PSDOT(i,j)
        VB1(i,k,j) = VB1(i,k,j)*PSDOT(i,j)
        TB1(i,k,j) = TB1(i,k,j)*PS1(i,j)
        QB1(i,k,j) = QB1(i,k,j)*PS1(i,j)
      end do
      end do
      end do

      mdate =ndate0
      tinit = (NNNCHK-NSTRT0)*6*60
      tfinal= (NNNCHK-NSTRT0+1)*6*60
      NNNCHK=NNNCHK+1

c     print*,'read in datasets at :',ndate1
c
c-----compute boundary conditions for p*:
c

#ifdef MPP1
      if (nspgx.le.jxp) then
         nxwb=nspgx
      else
         nkk=nspgx/jxp
         if (nspgx.eq.nkk*jxp) then
            nxwb=jxp
         else
            nxwb=nspgx-nkk*jxp
         end if
      end if
      if (nxwb+myid*jxp.gt.nspgx) then
         nxwb=0
      else if (nxwb+myid*jxp.lt.nspgx) then
         nxwb=jxp
      end if
c
      if (nspgx.le.jxp-1) then
         nxeb=nspgx
      else
         nkk=(nspgx-jxp+1)/jxp
         if ((nspgx-jxp+1).eq.nkk*jxp) then
            nxeb=jxp
         else
            nxeb=(nspgx-jxp+1)-nkk*jxp
         end if
      end if
      if (mjx-1-(myid*jxp+jxp-nxeb).gt.nspgx) then
         nxeb=0
      else if (mjx-1-(myid*jxp+jxp-nxeb).lt.nspgx) then
         nxeb=MIN(jendx,jxp)
      end if
      do nn=1,nxwb
#else
      do nn=1,nspgx
#endif
         do i=1,ilx 
            pwb(i,nn)=PS0(i,nn)
            pwbt(i,nn)=(PS1(i,nn )-PS0(i,nn ))/dtbdys
         end do
      end do
#ifdef MPP1
      do nn=1,nxeb
         nnb=MIN(jendx,jxp)-nn+1
#else
      do nn=1,nspgx
         nnb=jlx-nn+1
#endif
         do i=1,ilx 
            peb(i,nn)=PS0(i,nnb)
            pebt(i,nn)=(PS1(i,nnb)-PS0(i,nnb))/dtbdys
         end do
      end do
      do nn=1,nspgx
         nnb=ilx-nn+1
#ifdef MPP1
         do j=1,jendx
#else
         do j=1,jlx
#endif
            pnb(nn,j)=PS0(nnb,j)
            pss(nn,j)=PS0(nn ,j)
            pnbt(nn,j)=(PS1(nnb,j)-PS0(nnb,j))/dtbdys
            psbt(nn,j)=(PS1(nn ,j)-PS0(nn ,j))/dtbdys
         end do
      end do
c
c-----compute boundary conditions for p*u and p*v:
c
#ifdef MPP1
      if (nspgd.le.jxp) then
         ndwb=nspgd
      else
         nkk=nspgd/jxp
         if (nspgd.eq.nkk*jxp) then
            ndwb=jxp
         else
            ndwb=nspgd-nkk*jxp
         end if
      end if
      if (ndwb+myid*jxp.gt.nspgd) then
         ndwb=0
      else if (ndwb+myid*jxp.lt.nspgd) then
         ndwb=jxp
      end if
c
      if (nspgd.le.jendl) then
         ndeb=nspgd
      else
         nkk=nspgd/jxp
         if (nspgd.eq.nkk*jxp) then
            ndeb=jxp
         else
            ndeb=nspgd-nkk*jxp
         end if
      end if
      if (mjx-(myid*jxp+jxp-ndeb).gt.nspgd) then
         ndeb=0
      else if (mjx-(myid*jxp+jxp-ndeb).lt.nspgd) then
         ndeb=jxp
      end if
      do nn=1,ndwb
#else
      do nn=1,nspgd
#endif
         do k=1,kx
         do i=1,ix
            uwb(i,k,nn)=UB0(i,k,nn)
            vwb(i,k,nn)=VB0(i,k,nn)
            uwbt(i,k,nn)=(UB1(i,k,nn)-UB0(i,k,nn))/dtbdys
            vwbt(i,k,nn)=(VB1(i,k,nn)-VB0(i,k,nn))/dtbdys
         end do
         end do
      end do
#ifdef MPP1
      do nn=1,ndeb
        nnb=MIN(jendl,jxp)-nn+1
#else
      do nn=1,nspgd
         nnb=jx-nn+1
#endif
         do k=1,kx
         do i=1,ix
            ueb(i,k,nn)=UB0(i,k,nnb)
            veb(i,k,nn)=VB0(i,k,nnb)
            uebt(i,k,nn)=(UB1(i,k,nnb)-UB0(i,k,nnb))/dtbdys
            vebt(i,k,nn)=(VB1(i,k,nnb)-VB0(i,k,nnb))/dtbdys
         end do
         end do
      end do
      do nn=1,nspgd
         nnb=ix-nn+1
         do k=1,kx
#ifdef MPP1
         do j=1,jendl
#else
         do j=1,jx
#endif
            unb(nn,k,j)=UB0(nnb,k,j)
            usb(nn,k,j)=UB0(nn ,k,j)
            vnb(nn,k,j)=VB0(nnb,k,j)
            vsb(nn,k,j)=VB0(nn ,k,j)
            unbt(nn,k,j)=(UB1(nnb,k,j)-UB0(nnb,k,j))/dtbdys
            usbt(nn,k,j)=(UB1(nn ,k,j)-UB0(nn ,k,j))/dtbdys
            vnbt(nn,k,j)=(VB1(nnb,k,j)-VB0(nnb,k,j))/dtbdys
            vsbt(nn,k,j)=(VB1(nn ,k,j)-VB0(nn ,k,j))/dtbdys
         end do
         end do
      end do
c
c-----compute boundary conditions for p*t and p*qv:
c
#ifdef MPP1
      do nn=1,nxwb
#else
      do nn=1,nspgx
#endif
         do k=1,kx
         do i=1,ilx
            twb(i,k,nn)=TB0(i,k,nn)
            qwb(i,k,nn)=QB0(i,k,nn)
            twbt(i,k,nn)=(TB1(i,k,nn)-TB0(i,k,nn))/dtbdys
            qwbt(i,k,nn)=(QB1(i,k,nn)-QB0(i,k,nn))/dtbdys
         end do
         end do
      end do
#ifdef MPP1
      do nn=1,nxeb
        nnb=MIN(jendx,jxp)-nn+1
#else
      do nn=1,nspgx
         nnb=jlx-nn+1
#endif
         do k=1,kx
         do i=1,ilx
            teb(i,k,nn)=TB0(i,k,nnb)
            qeb(i,k,nn)=QB0(i,k,nnb)
            tebt(i,k,nn)=(TB1(i,k,nnb)-TB0(i,k,nnb))/dtbdys
            qebt(i,k,nn)=(QB1(i,k,nnb)-QB0(i,k,nnb))/dtbdys
         end do
         end do
      end do
      do nn=1,nspgx
         nnb=ilx-nn+1
         do k=1,kx
#ifdef MPP1
         do j=1,jendx
#else
         do j=1,jlx
#endif
            tnb(nn,k,j)=TB0(nnb,k,j)
            tsb(nn,k,j)=TB0(nn ,k,j)
            qnb(nn,k,j)=QB0(nnb,k,j)
            qsb(nn,k,j)=QB0(nn ,k,j)
            tnbt(nn,k,j)=(TB1(nnb,k,j)-TB0(nnb,k,j))/dtbdys
            tsbt(nn,k,j)=(TB1(nn ,k,j)-TB0(nn ,k,j))/dtbdys
            qnbt(nn,k,j)=(QB1(nnb,k,j)-QB0(nnb,k,j))/dtbdys
            qsbt(nn,k,j)=(QB1(nn ,k,j)-QB0(nn ,k,j))/dtbdys
         end do
         end do
      end do
#ifdef MPP1
      if (myid .eq. 0) then
#endif
      print*,'BCs are ready from ',ndate0,'  to ',ndate1
#ifdef MPP1
      end if
#endif
      idatex = ndate0
      ndate0=ndate1
#ifdef MPP1
      do j=1,jendx
#else
      do j=1,jlx
#endif
      do i=1,ilx
         tdum(i,j) = TS1(i,j)
      end do
      end do
      do k=1,kx
#ifdef MPP1
      do j=1,jendl
#else
      do j=1,jx
#endif
      do i=1,ix
         UB0(i,k,j) = UB1(i,k,j)
         VB0(i,k,j) = VB1(i,k,j)
         QB0(i,k,j) = QB1(i,k,j)
         TB0(i,k,j) = TB1(i,k,j)
      end do
      end do
      end do
#ifdef MPP1
      do j=1,jendl
#else
      do j=1,jx
#endif
      do i=1,ix
         PS0(i,j) = PS1(i,j)
         TS0(i,j) = TS1(i,j)
      end do
      end do
      if (EHSO4) then
        do k=1,kx
#ifdef MPP1
          do j=1,jendl
#else
          do j=1,jx
#endif
            do i=1,ix
               SO0(i,k,j) = SO1(i,k,j)
            end do
          end do
        end do
      end if
CbxqOCT2001_

      NNBASE= NNNNNN 
      nyear = mdate/1000000
      nmonth= (mdate-nyear*1000000)/10000

c-----------------------------------------------------------------------
      if ( ldatez .ge. ndate1) go to 10

#ifdef MPP1
      do 50 j = 1,jendx
#else
      do 50 j = 1,jx-1
#endif
      do 50 i = 1,ix-1
       if (veg2d(i,j).le.0.00001)then
        tga(i,j) = tdum(i,j)
        tgb(i,j) = tdum(i,j)
Celguindi if (tdum(i,j).le.271.35) then
c       if (tdum(i,j).le.271.38) then
c         print *,'Setting ocld2d to ice at i=',i,' j=',j,' t=',tdum(i,j)
c         do n=1,NNSG
c           ocld2d(n,i,j)=2.
c           sice2d(n,i,j)=1000.
c         end do
c       else
c         do n=1,NNSG
c           ocld2d(n,i,j)=0.
c           sice2d(n,i,j)=0.
c         end do
c       end if
       end if
50    continue

      return
      end subroutine bdyin

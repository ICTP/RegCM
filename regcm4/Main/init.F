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

      subroutine init
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine reads in the initial and boundary conditions.   c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
#     include "parameter.inc"
#     include "include/parrad.cb"
#     include "include/param1.cb"
#     include "include/param2.cb"
#     include "include/param3.cb"
#     include "include/iunits.cb"
#     include "include/bats.cb"
#     include "include/pmoist.cb"
#     include "include/main.cb"
#     include "include/mainchem.cb"
#     include "include/bdycod.cb"
#     include "include/cvaria.cb"
#     include "include/bats2.cb"
#     include "include/rad.cb"
#     include "include/trachem.cb"
#     include "include/message.cb"
#ifdef MPP1
      include 'mpif.h'
#     include "include/io.cb"
#endif
#ifdef DIAG
#     include "include/diagnosis.cb"
      real(kind=8)  tvmass,tcmass,tttmp
#endif
#     include "include/date.cb"
Cchem2
Cchem2_
#ifdef MPP1
      integer status(MPI_STATUS_SIZE)
      integer ierr
      real(kind=8)  psdot(ix,jxp)
      real(kind=8)   sav0(ix,kx*4+2,jxp),sav_0(ix,kx*4+2,mjx)
      real(kind=8)   sav0a(ix,kx+NNSG+5,jxp),sav_0a(ix,kx+NNSG+5,mjx)
      real(kind=8)   sav0c(ix,kx*2,jxp),sav_0c(ix,kx*2,mjx)
      real(kind=8)   sav0b(ix,kx+1,jxp),sav_0b(ix,kx+1,mjx)
      real(kind=8)   sav0s(ix,kx,jxp),sav_0s(ix,kx,mjx)
      COMMON /TMPSAV0/ sav0,sav_0,sav0a,sav_0a,sav0c,sav_0c
     &               ,sav0b,sav_0b,sav0s,sav_0s
      real(kind=8)   sav1(ilx,kx*4+(kx+1)*(kx+2),jxp)
      real(kind=8)  sav_1(ilx,kx*4+(kx+1)*(kx+2),mjx)
      COMMON /TMPSAV1/ sav1,sav_1
      real(kind=8)   sav2(ilx,NNSG*4+4,jxp),sav_2(ilx,NNSG*4+4,mjx)
      real(kind=8)   sav2a(ilx,NNSG*5+1,jxp),sav_2a(ilx,NNSG*5+1,mjx)
      COMMON /TMPSAV2/ sav2,sav_2,sav2a,sav_2a
      real(kind=8)   sav4(ix,ntr*(kx*4+1),jxp),sav4a(ilx,7,jxp)
      real(kind=8)  sav_4(ix,ntr*(kx*4+1),mjx),sav_4a(ilx,7,mjx)
      COMMON /TMPSAV4/ sav4,sav_4,sav4a,sav_4a
      real(kind=8)  inisrf0(ix,NNSG*2+7,jxp),inisrf_0(ix,NNSG*2+7,mjx)
      COMMON /TMPSRF1/ inisrf0,inisrf_0
      integer l,allrec
      real(kind=4)  io2d(ix,mjx)
#else
      real(kind=4)  io2d(ix,jx)
      real(kind=8)  psdot(ix,jx)
#endif
      INTEGER         MDATEZ
      COMMON /DATENUM/MDATEZ(289276)
C
      real(kind=8)  absnxt,              ! Nearest layer absorptivities
     &        abstot,              ! Non-adjacent layer absorptivites
     &        emstot               ! Total emissivity
#ifdef MPP1
      common /radbuf/ absnxt(plond,plev,4,jxp)
     &      , abstot(plond,plevp,plevp,jxp),emstot(plond,plevp,jxp)
      real(kind=8)  absnxt_io,
     &        abstot_io,
     &        emstot_io
      common /radbuf_io/absnxt_io(plond,plev,4,mjx-1)
     &                , abstot_io(plond,plevp,plevp,mjx-1)
     &                , emstot_io(plond,plevp,mjx-1)
#else
      common /radbuf/ absnxt(plond,plev,4,jlx)
     &      , abstot(plond,plevp,plevp,jlx),emstot(plond,plevp,jlx)
#endif
c
      integer lcount,iin,iout,numpts
      common /ldata/ lcount, iin, iout, numpts
      integer ilake, jlake, depth, freeze
      real(kind=8)  hi, hii, hs, eta, tlake(400)
      integer ibdyhr0, ibdyhr1, ibdydiff

c----------------------------------------------------------------------
c-----dimension the arrays for parameterizing the sfc. variables.
c     change the variable surface parameters
c
      integer i,j,k,ist,itr,ibin
      integer nyyy,nxxx,kzzz
c
      real(kind=8)  hg1,hg2,hg3,hg4,hgmax
      integer n,im1h,ip1h,jm1h,jp1h
      real(kind=8)  tlp,ts00
      data    tlp,ts00/50.d0,288.d0/
      character*8 finm
C_sgi character*7 finm
c----------------------------------------------------------------------
#ifdef MPP1
      do j=1,jxp
      do i=1,ix
        peb(i,j)=0.0
        pwb(i,j)=0.0
        pebt(i,j)=0.0
        pwbt(i,j)=0.0
      end do
      end do
      do j=1,jxp
      do k=1,kx
      do i=1,ix
        teb(i,k,j)=0.0
        twb(i,k,j)=0.0
        tebt(i,k,j)=0.0
        twbt(i,k,j)=0.0
        qeb(i,k,j)=0.0
        qwb(i,k,j)=0.0
        qebt(i,k,j)=0.0
        qwbt(i,k,j)=0.0
      end do
      end do
      end do
      do j=1,jxp
      do k=1,kx
      do i=1,ix
        ueb(i,k,j)=0.0
        uwb(i,k,j)=0.0
        uebt(i,k,j)=0.0
        uwbt(i,k,j)=0.0
        veb(i,k,j)=0.0
        vwb(i,k,j)=0.0
        vebt(i,k,j)=0.0
        vwbt(i,k,j)=0.0
      end do
      end do
      end do
      do j=1,jxp
      do i=1,nspgx
        pnb(i,j)=0.0
        pss(i,j)=0.0
        pnbt(i,j)=0.0
        psbt(i,j)=0.0
      end do
      end do
      do j=1,jxp
      do k=1,kx
      do i=1,nspgx
        tnb(i,k,j)=0.0
        tsb(i,k,j)=0.0
        tnbt(i,k,j)=0.0
        tsbt(i,k,j)=0.0
        qnb(i,k,j)=0.0
        qsb(i,k,j)=0.0
        qnbt(i,k,j)=0.0
        qsbt(i,k,j)=0.0
      end do
      end do
      end do
      do j=1,jxp
      do k=1,kx
      do i=1,nspgd
        unb(i,k,j)=0.0
        usb(i,k,j)=0.0
        unbt(i,k,j)=0.0
        usbt(i,k,j)=0.0
        vnb(i,k,j)=0.0
        vsb(i,k,j)=0.0
        vnbt(i,k,j)=0.0
        vsbt(i,k,j)=0.0
      end do
      end do
      end do
      do i=1,ix-2
      do j=1,jxp
        tgmx_o(j,i) = -1.e30
        t2mx_o(j,i) = -1.e30
        tgmn_o(j,i) =  1.e30
        t2mn_o(j,i) =  1.e30
        w10x_o(j,i) = -1.e30
        psmn_o(j,i) =  1.e30
      end do
      end do
#else
      do i=1,ix-2
      do j=1,jx-2
        tgmx_o(j,i) = -1.e30
        t2mx_o(j,i) = -1.e30
        tgmn_o(j,i) =  1.e30
        t2mn_o(j,i) =  1.e30
        w10x_o(j,i) = -1.e30
        psmn_o(j,i) =  1.e30
      end do
      end do
#endif
#ifdef MPP1
      if (myid .eq. 0) then
#endif
      write(finm,101) iutbc
 101  format('fort.',I3)
C_sgi 101  format('fort.',I2)
#ifdef MPP1
      open(iutbc,file=finm,form='unformatted',status='old'
     &          ,access='direct',recl=ix*mjx*ibyte)
#else
      open(iutbc,file=finm,form='unformatted',status='old'
     &          ,access='direct',recl=ix*jx*ibyte)
#endif
      mmrec = 0
#ifdef MPP1
      end if
#endif
c
      if (.not.ifrest) then
c-----for initial run--not using restart
c
c------set rainwater and cloud water equal to zero initially.
c
       do 10 k=1,kx
#ifdef MPP1
       do 10 j=1,jendx
#else
       do 10 j=1,jlx
#endif
       do 10 i=1,ilx
        qca(i,k,j) = 0.
        qcb(i,k,j) = 0.
10     continue
c
#ifdef DIAG
       tdini = 0.
       tdadv = 0.
       tqini = 0.
       tqadv = 0.
       tqeva = 0.
       tqrai = 0.
#endif
c
Cchem2
       if(ichem.eq.1) then
c----- total tracer concs (initial, emission, advected)
#ifdef DIAG
        do 20 itr = 1, ntr
         ttrace(itr,1) = 0.
         ttrace(itr,2) = 0.
         tchie(itr) = 0.
         tchiad(itr) = 0.
         tchitb(itr) = 0.
20      continue
#endif
cqhy tchie, tchitb(replace tchidp:deposition)
c   initialize removal terms
c
         do 615 itr=1,ntr
         do 615 k=1,kx
#ifdef MPP1
         do 615 j=1,jendl
#else
         do 615 j=1,jx
#endif
         do 615 i=1,ix
            remlsc(i,k,j,itr) = 0.
            remcvc(i,k,j,itr) = 0.
cqhy   rxsg et al.
            rxsg(i,k,j,itr) = 0.
            rxsaq1(i,k,j,itr) = 0.
            rxsaq2(i,k,j,itr) = 0.
CC          tchifx(i,j,k,itr) = 0.
CC          tchify(i,j,k,itr) = 0.
CC          tchidx(i,j,k,itr) = 0.
CC          tchidy(i,j,k,itr) = 0.
 615     continue

         do 616 itr=1,ntr
#ifdef MPP1
         do 616 j=1,jendl
#else
         do 616 j=1,jx
#endif
         do 616 i=1,ix
            remdrd(i,j,itr) = 0.
            wdlsc(i,j,ntr) = 0.
 616     continue
       end if
Cchem2_
c------set the variables related to blackadar pbl equal to 0 initially.
c
       if (ibltyp .ne. 0) then
#ifdef MPP1
        do 30 j=1,jendx
#else
        do 30 j=1,jlx
#endif
        do 30 i=1,ilx
         hfx(i,j) = 0.
         qfx(i,j) = 0.
30      continue
       end if
c
       if(icup.eq.1) then
       do 44 k=1,kx
#ifdef MPP1
       do 44 j=1,jendl
#else
       do 44 j=1,jx
#endif
       do 44 i=1,ix
        rsheat(i,k,j)=0.
        rswat(i,k,j)=0.
44     continue
       end if
c
c------read in the initial conditions for large domain:
c      the initial conditions are the output from PREPROC/ICBC.
c
#ifdef MPP1
        if (myid .eq. 0) then
#endif
 111    continue
        mmrec=mmrec+1
        read(iutbc,rec=mmrec) ndate0,nxxx,nyyy,kzzz
#ifdef MPP1
        if (nyyy.ne.ix .or. nxxx.ne.mjx .or. kzzz.ne.kx) then
          write (aline,*) 
     &      'SET IN regcm.param: IX=',ix,' JX=',jxp, ' KX=',kx
#else
        if (nyyy.ne.ix .or. nxxx.ne.jx .or. kzzz.ne.kx) then
          write (aline,*) 
     &      'SET IN regcm.param: IX=',ix,' JX=',jx, ' KX=',kx
#endif
          call say
          write (aline, *)
     &       'SET IN ICBC: NY=',nyyy,' NX=',nxxx,' NZ=',kzzz
          call fatal(__FILE__,__LINE__,
     &               'IMPROPER DIMENSION SPECIFICATION')
        end if
        print*,'READING INITAL CONDITIONS',ndate0
        if (ndate0.lt.MDATEZ(NNNCHK)) then
          print*,ndate0,MDATEZ(NNNCHK),NNNCHK
          print*,'read in datasets at :',ndate0
          if(EHSO4) then
             if(LSMTYP.ne.'USGS') then
                mmrec=mmrec+kx*5+2
             else
                mmrec=mmrec+kx*5+2+13
             end if
          else
             if(LSMTYP.ne.'USGS') then
                mmrec=mmrec+kx*4+2
             else
                mmrec=mmrec+kx*4+2+13
             end if
          end if
          print*,'Searching for proper date: ',ndate1,MDATEZ(NNNCHK+1)
          print*,ndate0,MDATEZ(NNNCHK)
          go to 111
        else if (ndate0.gt.MDATEZ(NNNCHK)) then
          write (aline, *) ndate0,MDATEZ(NNNCHK)
          call say
          call fatal(__FILE__,__LINE__,
     &        'DATE IN ICBC FILE EXCEEDED DATE IN RegCM')
        end if
C
        PRINT*,'U'
        do k=kx,1,-1
          mmrec=mmrec+1
#ifdef MPP1
          read(iutbc,rec=mmrec) ((io2d(i,j),j=1,mjx),i=1,ix)
          do j=1,mjx
          do i=1,ix
            ub0_io(i,k,j)=dble(io2d(i,j))
          end do
          end do
#else
          read(iutbc,rec=mmrec) ((io2d(i,j),j=1,jx),i=1,ix)
          do j=1,jx
          do i=1,ix
            ub0(i,k,j)=dble(io2d(i,j))
          end do
          end do
#endif
        end do
        PRINT*,'V'
        do k=kx,1,-1
          mmrec=mmrec+1
#ifdef MPP1
          read(iutbc,rec=mmrec) ((io2d(i,j),j=1,mjx),i=1,ix)
          do j=1,mjx
          do i=1,ix
            vb0_io(i,k,j)=dble(io2d(i,j))
          end do
          end do
#else
          read(iutbc,rec=mmrec) ((io2d(i,j),j=1,jx),i=1,ix)
          do j=1,jx
          do i=1,ix
            vb0(i,k,j)=dble(io2d(i,j))
          end do
          end do
#endif
        end do
        PRINT*,'TA'
        do k=kx,1,-1
          mmrec=mmrec+1
#ifdef MPP1
          read(iutbc,rec=mmrec) ((io2d(i,j),j=1,mjx),i=1,ix)
          do j=1,mjx
          do i=1,ix
            tb0_io(i,k,j)=dble(io2d(i,j))
          end do
          end do
#else
          read(iutbc,rec=mmrec) ((io2d(i,j),j=1,jx),i=1,ix)
          do j=1,jx
          do i=1,ix
            tb0(i,k,j)=dble(io2d(i,j))
          end do
          end do
#endif
        end do
        PRINT*,'QV'
        do k=kx,1,-1
          mmrec=mmrec+1
#ifdef MPP1
          read(iutbc,rec=mmrec) ((io2d(i,j),j=1,mjx),i=1,ix)
          do j=1,mjx
          do i=1,ix
            qb0_io(i,k,j)=dble(io2d(i,j))
          end do
          end do
#else
          read(iutbc,rec=mmrec) ((io2d(i,j),j=1,jx),i=1,ix)
          do j=1,jx
          do i=1,ix
            qb0(i,k,j)=dble(io2d(i,j))
          end do
          end do
#endif
        end do
        PRINT*,'PS'
        mmrec=mmrec+1
#ifdef MPP1
        read(iutbc,rec=mmrec) ((io2d(i,j),j=1,mjx),i=1,ix)
        do j=1,mjx
        do i=1,ix
          ps0_io(i,j)=dble(io2d(i,j))
        end do
        end do
#else
        read(iutbc,rec=mmrec) ((io2d(i,j),j=1,jx),i=1,ix)
        do j=1,jx
        do i=1,ix
          ps0(i,j)=dble(io2d(i,j))
        end do
        end do
#endif
        PRINT*,'TS'
        mmrec=mmrec+1
#ifdef MPP1
        read(iutbc,rec=mmrec) ((io2d(i,j),j=1,mjx),i=1,ix)
        do j=1,mjx
        do i=1,ix
          ts0_io(i,j)=dble(io2d(i,j))
        end do
        end do
#else
        read(iutbc,rec=mmrec) ((io2d(i,j),j=1,jx),i=1,ix)
        do j=1,jx
        do i=1,ix
          ts0(i,j)=dble(io2d(i,j))
        end do
        end do
#endif
        if(EHSO4) then
          PRINT*,'SO0'
          do k=kx,1,-1
            mmrec=mmrec+1
#ifdef MPP1
            read(iutbc,rec=mmrec) ((io2d(i,j),j=1,mjx),i=1,ix)
            do j=1,mjx
            do i=1,ix
              so0_io(i,k,j)=dble(io2d(i,j))
            end do
            end do
#else
            read(iutbc,rec=mmrec) ((io2d(i,j),j=1,jx),i=1,ix)
            do j=1,jx
            do i=1,ix
              so0(i,k,j)=dble(io2d(i,j))
            end do
            end do
#endif
          end do
        end if
        if(LSMTYP.eq.'USGS') mmrec=mmrec+13
#ifdef MPP1
          do j=1,mjx
            do k=1,kx
              do i=1,ix
                sav_0(i,      k,j) = UB0_io(i,k,j)
                sav_0(i,kx   +k,j) = VB0_io(i,k,j)
                sav_0(i,kx*2 +k,j) = QB0_io(i,k,j)
                sav_0(i,kx*3 +k,j) = TB0_io(i,k,j)
              end do
            end do
            do i=1,ix
              sav_0(i,kx*4+1 ,j) = PS0_io(i,j)
              sav_0(i,kx*4+2 ,j) = TS0_io(i,j)
            end do
          end do
          if(EHSO4) then
            do j=1,mjx
              do k=1,kx
                do i=1,ix
                  sav_0s(i,k,j) = SO0_io(i,k,j)
                end do
              end do
            end do
          end if
        end if            ! end if myid=0
        call MPI_BCAST(ndate0,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(nxxx,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(nyyy,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(kzzz,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
        call MPI_SCATTER(
     &           sav_0(1,1,1),ix*(kx*4+2)*jxp,MPI_DOUBLE_PRECISION,
     &            sav0(1,1,1),ix*(kx*4+2)*jxp,MPI_DOUBLE_PRECISION,
     &                                      0,MPI_COMM_WORLD,ierr)
        if(EHSO4)
     &  call MPI_SCATTER(
     &           sav_0s(1,1,1),ix*kx*jxp,MPI_DOUBLE_PRECISION,
     &            sav0s(1,1,1),ix*kx*jxp,MPI_DOUBLE_PRECISION,
     &                                 0,MPI_COMM_WORLD,ierr)
        do j=1,jendl
          do k=1,kx
            do i=1,ix
              UB0(i,k,j) = sav0(i,      k,j)
              VB0(i,k,j) = sav0(i,kx   +k,j)
              QB0(i,k,j) = sav0(i,kx*2 +k,j)
              TB0(i,k,j) = sav0(i,kx*3 +k,j)
            end do
          end do
          do i=1,ix
            PS0(i,j) = sav0(i,kx*4+1 ,j)
            TS0(i,j) = sav0(i,kx*4+2 ,j)
          end do
          if(EHSO4) then
            do k=1,kx
              do i=1,ix
                SO0(i,k,j) = sav0s(i,k,j)
              end do
            end do
          end if
        end do
#endif
c
c Convert surface pressure to pstar
#ifdef MPP1
        do j=1,jendl
#else
        do j=1,jx
#endif
          do i=1,ix
            PS0(i,j) = PS0(i,j) - ptop
          end do
        end do
C=======================================================================
c
c  this routine determines p(.) from p(x) by a 4-point interpolation.
c  on the x-grid, a p(x) point outside the grid domain is assumed to
c  satisfy p(0,j)=p(1,j); p(ix,j)=p(ix-1,j); and similarly for the i's.
c
#ifdef MPP1
      call MPI_SENDRECV(PS0(1,jxp),ix,MPI_DOUBLE_PRECISION,ieast,1,
     &                  PS0(1,0),ix,MPI_DOUBLE_PRECISION,iwest,1,
     &                  MPI_COMM_WORLD,status,ierr)
#endif
#ifdef MPP1
      do j=jbegin,jendx
#else
      do j=2,jlx
#endif
        do i=2,ilx
          psdot(i,j)=0.25*(PS0(i,j)+PS0(i-1,j)+PS0(i,j-1)+PS0(i-1,j-1))
        end do
      end do
c
      do i=2,ilx
#ifdef MPP1
        if (myid .eq. 0) then
#endif
        psdot(i,1)=0.5*(PS0(i,1)+PS0(i-1,1))
#ifdef MPP1
        end if
        if (myid .eq. nproc-1) then
          psdot(i,jendl)=0.5*(PS0(i,jendx)+PS0(i-1,jendx))
        end if
#else
        psdot(i,jx)=0.5*(PS0(i,jlx)+PS0(i-1,jlx))
#endif
      end do
c
#ifdef MPP1
      do j=jbegin,jendx
#else
      do j=2,jlx
#endif
        psdot(1,j)=0.5*(PS0(1,j)+PS0(1,j-1))
        psdot(ix,j)=0.5*(PS0(ilx,j)+PS0(ilx,j-1))
      end do
c
#ifdef MPP1
      if (myid .eq. 0) then
#endif
      psdot(1,1)=PS0(1,1)
      psdot(ix,1)=PS0(ilx,1)
#ifdef MPP1
      end if
      if (myid .eq. nproc-1) then
        psdot(1,jendl)=PS0(1,jendx)
        psdot(ix,jendl)=PS0(ilx,jendx)
      end if
#else
      psdot(1,jx)=PS0(1,jlx)
      psdot(ix,jx)=PS0(ilx,jlx)
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
          UB0(i,k,j) = UB0(i,k,j)*PSDOT(i,j)
          VB0(i,k,j) = VB0(i,k,j)*PSDOT(i,j)
          TB0(i,k,j) = TB0(i,k,j)*PS0(i,j)
          QB0(i,k,j) = QB0(i,k,j)*PS0(i,j)
        end do
        end do
        end do
#ifdef MPP1
        if (myid .eq. 0) then
#endif
        mmrec = mmrec + 1
        read(iutbc,rec=mmrec) ndate1
        mmrec = mmrec - 1
        ibdyhr0 = ndate0 - (ndate0/100)*100
        ibdyhr1 = ndate1 - (ndate1/100)*100
        if (ibdyhr1.eq.0) ibdyhr1 = 24
        ibdydiff = ibdyhr1-ibdyhr0
        if (ibdydiff.ne.ibdyfrq) then
          write(aline, *) '  ndate0=',ndate0,'ndate1=',ndate1
          call say
          write(aline, *) '  ibdyfrq=',ibdyfrq,'ibdydiff=',ibdydiff
          call say
          write(aline, *) '  ibdyhr0=',ibdyhr0,'ibdyhr1=',ibdyhr1
          call say
          call fatal(__FILE__,__LINE__,
     &        'BOUNDARY CONDITION FREQUENCY INCOMPATIBILITY')
        end if
#ifdef MPP1
        end if
#endif
        mdate=ndate0
c Initialize variables and convert to double precision
        do k=1,kx
#ifdef MPP1
        do j=1,jendl
#else
        do j=1,jx
#endif
        do i=1,ix
           ua(i,k,j)  = ub0(i,k,j)
           ub(i,k,j)  = ub0(i,k,j)
           va(i,k,j)  = vb0(i,k,j)
           vb(i,k,j)  = vb0(i,k,j)
           qva(i,k,j) = qb0(i,k,j)
           qvb(i,k,j) = qb0(i,k,j)
           ta(i,k,j)  = tb0(i,k,j)
           tb(i,k,j)  = tb0(i,k,j)
        end do
        end do
        end do
#ifdef MPP1
        do j=1,jendl
#else
        do j=1,jx
#endif
        do i=1,ix
           psa(i,j) = ps0(i,j)
           psb(i,j) = ps0(i,j)
           tga(i,j) = ts0(i,j)
           tgb(i,j) = ts0(i,j)
        end do
        end do
        do k=1,kx
#ifdef MPP1
        do j=1,jendl
#else
        do j=1,jx
#endif
        do i=1,ix
           TBASE(i,k,j)=TS00+TLP*DLOG((psa(i,j)*a(k)+ptop)/100.)
        end do
        end do
        end do
       if(EHSO4) then
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
c
#ifdef MPP1
       do 129 j = 1,jendx
#else
       do 129 j = 1,jlx
#endif
       do 129 i = 1,ilx
        tga(i,j)=ta(i,kx,j)/psa(i,j)
        tgb(i,j)=tb(i,kx,j)/psb(i,j)
        tgbb(i,j)=tb(i,kx,j)/psb(i,j)
        zpbl(i,j)=500.   ! For Zeng Ocean Flux Scheme
129    continue
#ifdef MPP1
       do j=1,jendx
#else
       do j=1,jlx
#endif
       do i=1,ilx
       do k=1,NNSG
        snowc(k,i,j)=0.
       end do
       end do
       end do
       if(ichem.eq.1) then
#ifdef MPP1
         do j=1,jendx
#else
         do j=1,jlx
#endif
         do i=1,ilx
           ssw2da(i,j) = 0.0
           sdeltk2d(i,j) = 0.0
           sdelqk2d(i,j) = 0.0
           sfracv2d(i,j) = 0.5
           sfracb2d(i,j) = 0.5
           sfracs2d(i,j) = 0.0
           svegfrac2d(i,j) = 0.0
         end do
         end do
       end if
#ifdef DIAG
C=======================================================================
c
c-----dry air (unit = kg):
c
      tdini = 0.
#ifdef MPP1
      call MPI_GATHER(psa(1,1),ix*jxp,MPI_DOUBLE_PRECISION,
     &             psa_io(1,1),ix*jxp,MPI_DOUBLE_PRECISION,
     &                              0,MPI_COMM_WORLD,ierr)
      if(myid.eq.0) then
        do k=1,kx
          tttmp= 0.
          do j=1,mjx-1
            do i=1,ix-1
              tttmp=tttmp+psa_io(i,j)
            end do
          end do
          tdini=tdini+tttmp*dsigma(k)
        end do
        tdini=tdini*dx*dx*1000./g
      end if
      call MPI_BCAST(tdini,1,MPI_DOUBLE_PRECISION,
     &                     0,MPI_COMM_WORLD,ierr)
#else
      do k=1,kx
        tttmp= 0.
        do j=1,jx-1
          do i=1,ix-1
            tttmp=tttmp+psa(i,j)
          end do
        end do
        tdini=tdini+tttmp*dsigma(k)
      end do
      tdini=tdini*dx*dx*1000./g
#endif
c
c-----water substance (unit = kg):
c
      tvmass = 0.
#ifdef MPP1
      call MPI_GATHER(qva(1,1,1),ix*kx*jxp,MPI_DOUBLE_PRECISION,
     &             qva_io(1,1,1),ix*kx*jxp,MPI_DOUBLE_PRECISION,
     &                                   0,MPI_COMM_WORLD,ierr)
      if(myid.eq.0) then
        do k=1,kx
          tttmp= 0.
          do j=1,mjx-1
            do i=1,ix-1
              tttmp=tttmp+qva_io(i,k,j)
            end do
          end do
          tvmass=tvmass+tttmp*dsigma(k)
        end do
        tvmass=tvmass*dx*dx*1000./g
      end if
      call MPI_BCAST(tvmass,1,MPI_DOUBLE_PRECISION,
     &                      0,MPI_COMM_WORLD,ierr)
#else
      do k=1,kx
        tttmp= 0.
        do j=1,jx-1
          do i=1,ix-1
            tttmp=tttmp+qva(i,k,j)
          end do
        end do
        tvmass=tvmass+tttmp*dsigma(k)
      end do
      tvmass=tvmass*dx*dx*1000./g
#endif
c
      tcmass = 0.
#ifdef MPP1
      call MPI_GATHER(qca(1,1,1),ix*kx*jxp,MPI_DOUBLE_PRECISION,
     &             qca_io(1,1,1),ix*kx*jxp,MPI_DOUBLE_PRECISION,
     &                                   0,MPI_COMM_WORLD,ierr)
      if(myid.eq.0) then
        do k=1,kx
          tttmp= 0.
          do j=1,mjx-1
            do i=1,ix-1
              tttmp=tttmp+qca_io(i,k,j)
            end do
          end do
          tcmass=tcmass+tttmp*dsigma(k)
        end do
        tcmass=tcmass*dx*dx*1000./g
      end if
      call MPI_BCAST(tcmass,1,MPI_DOUBLE_PRECISION,
     &                      0,MPI_COMM_WORLD,ierr)
#else
      do k=1,kx
        tttmp= 0.
        do j=1,jx-1
          do i=1,ix-1
            tttmp=tttmp+qca(i,k,j)
          end do
        end do
        tcmass=tcmass+tttmp*dsigma(k)
      end do
      tcmass=tcmass*dx*dx*1000./g
#endif
      tqini=tvmass+tcmass
C=======================================================================
#ifdef MPP1
       if (myid .eq. 0) then
#endif
       print 150,tdini,tqini
  150  format(' *** initial total air = ',e12.5,' kg, total water = ',
     1       e12.5,' kg in large domain.')
#ifdef MPP1
       end if
#endif
#endif
c
Cchem2
      if(ichem.eq.1) then
c-----set tracer concs to 1 (kg/kg) initially. Must convert this to p* mass
c-----mixing ratio to compute tendencies:
CUS test zero concs init input for advection
cqhy  initial chia is 10ppt
chy set the initial tracer concentration 10ppt (1.e-11), 9/4/98

       do 135 itr=1,ntr
       do 135 k=1,kx
#ifdef MPP1
       do 135 j=1,jendx
#else
       do 135 j=1,jlx
#endif
       do 135 i=1,ilx
         chia(i,k,j,itr)=psa(i,j)*0.0d0
         chib(i,k,j,itr)=psb(i,j)*0.0d0
c        chia(i,k,j,itr)=psa(i,j)*1.e-11
c        chib(i,k,j,itr)=psb(i,j)*1.e-11
135    continue

      end if
Cchem2_
c
c------set rainc and rainnc equal to 0. initially
c
#ifdef MPP1
       do 280 j=1,jendx
#else
       do 280 j=1,jlx
#endif
       do 280 i=1,ilx
        rainc(i,j) = 0.
        rainnc(i,j) = 0.
280    continue

       if (icup.eq.4) then
#ifdef MPP1
         do j=1,jendx
#else
         do j=1,jlx
#endif
         do i=1,ilx
           cbmf2d(i,j) = 0.
         end do
         end do
       end if
c
      else
c-----when ifrest=.true., read in the data saved from previous run
c     for large domain from unit 14.
c
#ifdef MPP1
       if (myid .eq. 0) then
#endif
580    read(iutrs) mdate0
       jyear0 = mdate0/1000000
       read(iutrs) ktau,xtime,ldatez,lyear,lmonth,lday,lhour,ntime
       jyear = lyear
       jyearr= jyear
       ktaur = ktau
#ifdef MPP1
       if(EHSO4) then
          read(iutrs) UB0_IO,VB0_IO,QB0_IO,TB0_IO,PS0_IO,TS0_IO,SO0_IO
       else
          read(iutrs) UB0_IO,VB0_IO,QB0_IO,TB0_IO,PS0_IO,TS0_IO
       end if
       read(iutrs) ua_io
       read(iutrs) ub_io
       read(iutrs) va_io
       read(iutrs) vb_io
       read(iutrs) ta_io
       read(iutrs) tb_io
       read(iutrs) qva_io
       read(iutrs) qvb_io
       read(iutrs) qca_io
       read(iutrs) qcb_io
       read(iutrs) psa_io,psb_io,satbrt_io,satbrt1_io,f_io
       read(iutrs) ht_io,ht1_io,msfx_io,msfd_io,xlat_io,xlong_io
       read(iutrs) tga_io,tgb_io,rainc_io,rainnc_io
       if(icup.eq.1) then
         read(iutrs) rsheat_io,  rswat_io
       else if(icup.eq.3) then
         read(iutrs) TBASE_io,CLDEFI_io
       else if(icup.eq.4) then
         read(iutrs) cbmf2d_io
       end if
       read(iutrs) hfx_io,qfx_io,snowc_io,uvdrag_io
#ifdef DIAG
       read(iutrs) tdini,tdadv,tqini,tqadv,tqeva,tqrai
#endif
       read(iutrs) absnxt_io, abstot_io, emstot_io
       if(ipptls.eq.1) read(iutrs) fcc_io
       read(iutrs)  sol2d_io, solvd2d_io, solvs2d_io
     &           ,  flw2d_io,  flwd2d_io, fsw2d_io
     &           , sabv2d_io,  sinc2d_io
       read(iutrs)  taf2d_io, tlef2d_io, tgbb_io
     b           ,  ssw2d_io,  srw2d_io
     c           ,  tg2d_io,   tgb2d_io,  swt2d_io
     d           ,  scv2d_io, gwet2d_io,  veg2d_io,veg2d1_io
     e           ,  sag2d_io, sice2d_io,  dew2d_io
     &           , ircp2d_io
     g           , text2d_io,  col2d_io, ocld2d_io
     h           , heatrt_io, o3prof_io
       read(iutrs) pptnc_io, pptc_io, prca2d_io, prnca2d_io
       if (iocnflx.eq.2) read(iutrs)  zpbl_io
Cchem2---
       if(ichem.eq.1) then
        read(iutrs) chia_io
        read(iutrs) chib_io
c cumul removal terms (3d, 2d)
        read(iutrs) remlsc_io
        read(iutrs) remcvc_io
        read(iutrs) remdrd_io
        read(iutrs) ssw2da_io
        read(iutrs) sdeltk2d_io
        read(iutrs) sdelqk2d_io
        read(iutrs) sfracv2d_io
        read(iutrs) sfracb2d_io
        read(iutrs) sfracs2d_io
        read(iutrs) svegfrac2d_io
c cumul ad, dif, emis terms ( scalar)
#ifdef DIAG
        read(iutrs) tchiad
        read(iutrs) tchitb
        read(iutrs) tchie        
#endif
       end if
#else
       if(EHSO4) then
          read(iutrs) UB0,VB0,QB0,TB0,PS0,TS0,SO0
       else
          read(iutrs) UB0,VB0,QB0,TB0,PS0,TS0
       end if
       read(iutrs) ua
       read(iutrs) ub
       read(iutrs) va
       read(iutrs) vb
       read(iutrs) ta
       read(iutrs) tb
       read(iutrs) qva
       read(iutrs) qvb
       read(iutrs) qca
       read(iutrs) qcb
       read(iutrs) psa,psb,satbrt,satbrt1,f
       read(iutrs) ht,ht1,msfx,msfd,xlat,xlong
       read(iutrs) tga,tgb,rainc,rainnc
       if(icup.eq.1) then
         read(iutrs) rsheat,  rswat
       else if(icup.eq.3) then
         read(iutrs) TBASE,CLDEFI
       else if(icup.eq.4) then
         read(iutrs) cbmf2d
       end if
       read(iutrs) hfx,qfx,snowc,uvdrag
#ifdef DIAG
       read(iutrs) tdini,tdadv,tqini,tqadv,tqeva,tqrai
#endif
       read(iutrs) absnxt, abstot, emstot
       if(ipptls.eq.1) read(iutrs) fcc
       read(iutrs)  sol2d, solvd2d, solvs2d
     &           ,  flw2d,  flwd2d, fsw2d
     &           , sabv2d,  sinc2d
       read(iutrs)  taf2d, tlef2d, tgbb
     b           ,  ssw2d,  srw2d
     c           ,  tg2d,   tgb2d,  swt2d
     d           ,  scv2d, gwet2d,  veg2d,veg2d1
     e           ,  sag2d, sice2d,  dew2d,ircp2d
     g           , text2d,  col2d, ocld2d
     h           , heatrt, o3prof
       read(iutrs) pptnc, pptc, prca2d, prnca2d
       if (iocnflx.eq.2) read(iutrs)  zpbl

Cchem2---
       if(ichem.eq.1) then
        read(iutrs) chia
        read(iutrs) chib
c cumul removal terms (3d, 2d)
        read(iutrs) remlsc
        read(iutrs) remcvc
        read(iutrs) remdrd
        read(iutrs) ssw2da
        read(iutrs) sdeltk2d
        read(iutrs) sdelqk2d
        read(iutrs) sfracv2d
        read(iutrs) sfracb2d
        read(iutrs) sfracs2d
        read(iutrs) svegfrac2d
c cumul ad, dif, emis terms ( scalar)
#ifdef DIAG
        read(iutrs) tchiad
        read(iutrs) tchitb
        read(iutrs) tchie        
#endif
       end if
Cchem2_
#endif
       
c------lake model
       if(lakemod.eq.1)then
        lcount = 0
        iin = 41
        iout = 42
        rewind(iin)
        read (iutrs) numpts
        print *, 'reading lake model restart file. numpts = ', numpts
        print *, 'jyear, ktau, xtime = ', jyear, ktau, xtime
        do n=1,numpts
         read(iutrs) ilake, jlake, depth, freeze, hi ,hii, hs
     &            , eta, (tlake(j),j=1,depth)
         print *, 'reading restart file at i, j = ', ilake, jlake
         write (iin) ilake, jlake, depth, freeze, hi, hii, hs
     &            , eta, (tlake(j),j=1,depth)
        end do
        rewind(iin)
       end if
c
       print *,'ozone profiles restart'
       do 34 k=1,kx+1
#ifdef MPP1
        write(6,1212)o3prof_io(3,3,k)
#else
        write(6,1212)o3prof(3,3,k)
#endif
34     continue
1212   format (1x,7e12.4)
       print 590, xtime,ktau,jyear,iutrs
590    format(' ***** restart file for large domain at time = ',f8.0,
     &       ' minutes, ktau = ',i7,' in year = ',i4,
     &       '  read in from unit ',i2,'.'/)
c
       if (ldatez.ne.idate1) then
          write(*,*) 'INIT: ldatez, idate1=', ldatez, idate1
          go to 580
       end if
#ifdef MPP1
        do j=1,mjx
          do i=1,ix
            inisrf_0(i,1,j) = ht_io(i,j)
            inisrf_0(i,2,j) = satbrt_io(i,j)
            inisrf_0(i,3,j) = xlat_io(i,j)
            inisrf_0(i,4,j) = xlong_io(i,j)
            inisrf_0(i,5,j) = msfx_io(i,j)
            inisrf_0(i,6,j) = msfd_io(i,j)
            inisrf_0(i,7,j) = f_io(i,j)
          end do
          do n=1,NNSG
            do i=1,ix
              inisrf_0(i,7+n     ,j) = ht1_io(n,i,j)
              inisrf_0(i,7+NNSG+n  ,j) = satbrt1_io(n,i,j)
            end do
          end do
        end do
            
          do j=1,mjx
            do k=1,kx
              do i=1,ix
                sav_0(i,      k,j) = UB0_io(i,k,j)
                sav_0(i,kx   +k,j) = VB0_io(i,k,j)
                sav_0(i,kx*2 +k,j) = QB0_io(i,k,j)
                sav_0(i,kx*3 +k,j) = TB0_io(i,k,j)
              end do
            end do
            do i=1,ix
              sav_0(i,kx*4+1 ,j) = PS0_io(i,j)
              sav_0(i,kx*4+2 ,j) = TS0_io(i,j)
            end do
            if(EHSO4) then
              do  k=1,kx
                do i=1,ix
                  sav_0s(i,k,j) = SO0_io(i,k,j)
                end do
              end do
            end if
          end do
       end if
       call MPI_SCATTER(inisrf_0(1,1,1),ix*(NNSG*2+7)*jxp,
     &       MPI_DOUBLE_PRECISION,inisrf0(1,1,1),ix*(NNSG*2+7)*jxp,
     &       MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
       do j=1,jxp
         do i=1,ix
           ht(i,j)     = inisrf0(i,1,j)
           satbrt(i,j) = inisrf0(i,2,j)
           xlat(i,j)   = inisrf0(i,3,j)
           xlong(i,j)  = inisrf0(i,4,j)
           msfx(i,j)   = inisrf0(i,5,j)
           msfd(i,j)   = inisrf0(i,6,j)
           f(i,j)      = inisrf0(i,7,j)
         end do
         do n=1,NNSG
           do i=1,ix
             ht1(n,i,j)     = inisrf0(i,7+n     ,j)
             satbrt1(n,i,j) = inisrf0(i,7+NNSG+n  ,j)
           end do
         end do
       end do
       call MPI_SCATTER(
     &          sav_0(1,1,1),ix*(kx*4+2)*jxp,MPI_DOUBLE_PRECISION,
     &           sav0(1,1,1),ix*(kx*4+2)*jxp,MPI_DOUBLE_PRECISION,
     &                                     0,MPI_COMM_WORLD,ierr)
       if(EHSO4)
     & call MPI_SCATTER(
     &          sav_0s(1,1,1),ix*kx*jxp,MPI_DOUBLE_PRECISION,
     &           sav0s(1,1,1),ix*kx*jxp,MPI_DOUBLE_PRECISION,
     &                                0,MPI_COMM_WORLD,ierr)
       do j=1,jendl
         do k=1,kx
           do i=1,ix
             UB0(i,k,j) = sav0(i,     k,j)
             VB0(i,k,j) = sav0(i,kx  +k,j)
             QB0(i,k,j) = sav0(i,kx*2+k,j)
             TB0(i,k,j) = sav0(i,kx*3+k,j)
           end do
         end do
         do i=1,ix
           PS0(i,j) = sav0(i,kx*4+1 ,j)
           TS0(i,j) = sav0(i,kx*4+2 ,j)
         end do
         if(EHSO4) then
           do k=1,kx
             do i=1,ix
               SO0(i,k,j) = sav0s(i,k,j)
             end do
           end do
         end if
       end do
       if (myid .eq. 0) then
          do j=1,mjx
            do k=1,kx
              do i=1,ix
                sav_0(i,      k,j) = ua_io(i,k,j)
                sav_0(i,kx   +k,j) = ub_io(i,k,j)
                sav_0(i,kx*2 +k,j) = va_io(i,k,j)
                sav_0(i,kx*3 +k,j) = vb_io(i,k,j)
              end do
            end do
            do i=1,ix
              sav_0(i,kx*4+1 ,j) = psa_io(i,j)
              sav_0(i,kx*4+2 ,j) = psb_io(i,j)
            end do
          end do
       end if
       call MPI_SCATTER(
     &          sav_0(1,1,1),ix*(kx*4+2)*jxp,MPI_DOUBLE_PRECISION,
     &           sav0(1,1,1),ix*(kx*4+2)*jxp,MPI_DOUBLE_PRECISION,
     &                                     0,MPI_COMM_WORLD,ierr)
       do j=1,jendl
         do k=1,kx
           do i=1,ix
             ua(i,k,j) = sav0(i,     k,j)
             ub(i,k,j) = sav0(i,kx  +k,j)
             va(i,k,j) = sav0(i,kx*2+k,j)
             vb(i,k,j) = sav0(i,kx*3+k,j)
           end do
         end do
         do i=1,ix
           psa(i,j) = sav0(i,kx*4+1 ,j)
           psb(i,j) = sav0(i,kx*4+2 ,j)
         end do
       end do
       if (myid .eq. 0) then
          do j=1,mjx
            do k=1,kx
              do i=1,ix
                sav_0(i,      k,j) = ta_io(i,k,j)
                sav_0(i,kx   +k,j) = tb_io(i,k,j)
                sav_0(i,kx*2 +k,j) =qva_io(i,k,j)
                sav_0(i,kx*3 +k,j) =qvb_io(i,k,j)
              end do
            end do
            do i=1,ix
              sav_0(i,kx*4+1 ,j) = tga_io(i,j)
              sav_0(i,kx*4+2 ,j) = tgb_io(i,j)
            end do
          end do
       end if
       call MPI_SCATTER(
     &          sav_0(1,1,1),ix*(kx*4+2)*jxp,MPI_DOUBLE_PRECISION,
     &           sav0(1,1,1),ix*(kx*4+2)*jxp,MPI_DOUBLE_PRECISION,
     &                                     0,MPI_COMM_WORLD,ierr)
       do j=1,jendl
         do k=1,kx
           do i=1,ix
             ta(i,k,j) = sav0(i,     k,j)
             tb(i,k,j) = sav0(i,kx  +k,j)
            qva(i,k,j) = sav0(i,kx*2+k,j)
            qvb(i,k,j) = sav0(i,kx*3+k,j)
           end do
         end do
         do i=1,ix
           tga(i,j) = sav0(i,kx*4+1 ,j)
           tgb(i,j) = sav0(i,kx*4+2 ,j)
         end do
       end do
       if (myid .eq. 0) then
          do j=1,mjx
            do k=1,kx
              do i=1,ix
                sav_0(i,      k,j) =qca_io(i,k,j)
                sav_0(i,kx   +k,j) =qcb_io(i,k,j)
                sav_0(i,kx*2 +k,j) =fcc_io(i,k,j)
              end do
            end do
            do i=1,ix
              sav_0(i,kx*4+1 ,j) = rainc_io(i,j)
              sav_0(i,kx*4+2 ,j) = rainnc_io(i,j)
            end do
          end do
          do j=1,mjx-1
            do k=1,kx
              do i=1,ilx
                sav_0(i,kx*3+k,j) =heatrt_io(i,k,j)
              end do
            end do
          end do
       end if
       call MPI_SCATTER(
     &          sav_0(1,1,1),ix*(kx*4+2)*jxp,MPI_DOUBLE_PRECISION,
     &           sav0(1,1,1),ix*(kx*4+2)*jxp,MPI_DOUBLE_PRECISION,
     &                                     0,MPI_COMM_WORLD,ierr)
       do j=1,jendl
         do k=1,kx
           do i=1,ix
            qca(i,k,j) = sav0(i,     k,j)
            qcb(i,k,j) = sav0(i,kx  +k,j)
            fcc(i,k,j) = sav0(i,kx*2+k,j)
           end do
         end do
         do i=1,ix
           rainc(i,j) = sav0(i,kx*4+1 ,j)
           rainnc(i,j) = sav0(i,kx*4+2 ,j)
         end do
       end do
       do j=1,jendx
         do k=1,kx
           do i=1,ilx
             heatrt(i,k,j) = sav0(i,kx*3+k,j)
           end do
         end do
       end do
       if (myid .eq. 0) then
          do j=1,mjx
            do i=1,ix
              sav_0a(i,1,j) = hfx_io(i,j)
              sav_0a(i,2,j) = qfx_io(i,j)
              sav_0a(i,3,j) = uvdrag_io(i,j)
              sav_0a(i,4,j) = tgbb_io(i,j)
            end do
            do n=1,NNSG
              do i=1,ix
                sav_0a(i,4+n,j) = snowc_io(n,i,j)
              end do
            end do
          end do
          do j=1,mjx-1
            do k=1,kx+1
              do i=1,ilx
                sav_0a(i,NNSG+4 +k,j) = o3prof_io(i,k,j)
              end do
            end do
          end do
       end if
       allrec = kx+5+NNSG
       call MPI_SCATTER(
     &          sav_0a(1,1,1),ix*allrec*jxp,MPI_DOUBLE_PRECISION,
     &           sav0a(1,1,1),ix*allrec*jxp,MPI_DOUBLE_PRECISION,
     &                                   0,MPI_COMM_WORLD,ierr)
       do j=1,jendl
         do i=1,ix
           hfx(i,j)    = sav0a(i,1 ,j)
           qfx(i,j)    = sav0a(i,2 ,j)
           uvdrag(i,j) = sav0a(i,3 ,j)
           tgbb(i,j)   = sav0a(i,4 ,j)
         end do
         do n=1,NNSG
           do i=1,ix
             snowc(n,i,j) = sav0a(i,4+n,j)
           end do
         end do
       end do
       do j=1,jendx
         do k=1,kx+1
           do i=1,ilx
             o3prof(i,k,j) = sav0a(i,NNSG+4 +k,j)
           end do
         end do
       end do
       if(iocnflx.eq.2) then
         call MPI_SCATTER(zpbl_io(1,1),ix*jxp,MPI_DOUBLE_PRECISION,
     &                       zpbl(1,1),ix*jxp,MPI_DOUBLE_PRECISION,
     &                                      0,MPI_COMM_WORLD,ierr)
       end if
       if(icup.eq.1) then
         if (myid .eq. 0) then
           do j=1,mjx
             do k=1,kx
               do i=1,ix
                 sav_0c(i,   k,j) = rsheat_io(i,k,j)
                 sav_0c(i,kx+k,j) =  rswat_io(i,k,j)
               end do
             end do
           end do
         end if
         call MPI_SCATTER(
     &            sav_0c(1,1,1),ix*kx*2*jxp,MPI_DOUBLE_PRECISION,
     &             sav0c(1,1,1),ix*kx*2*jxp,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
         do j=1,jendl
           do k=1,kx
             do i=1,ix
               rsheat(i,k,j) = sav0c(i,   k,j)
                rswat(i,k,j) = sav0c(i,kx+k,j)
             end do
           end do
         end do
       else if(icup.eq.3) then
         if (myid .eq. 0) then
           do j=1,mjx
             do k=1,kx
               do i=1,ix
                 sav_0b(i,   k,j) = TBASE_IO(i,k,j)
               end do
             end do
             do i=1,ix
               sav_0b(i,kx+1,j) = CLDEFI_IO(i,j)
             end do
           end do
         end if
         call MPI_SCATTER(
     &            sav_0b(1,1,1),ix*(kx+1)*jxp,MPI_DOUBLE_PRECISION,
     &             sav0b(1,1,1),ix*(kx+1)*jxp,MPI_DOUBLE_PRECISION,
     &                                      0,MPI_COMM_WORLD,ierr)
         do j=1,jendl
           do k=1,kx
             do i=1,ix
               TBASE(i,k,j) = sav0b(i,   k,j)
             end do
           end do
           do i=1,ix
             CLDEFI(i,j) = sav0b(i,kx+1,j)
           end do
         end do
       else if(icup.eq.4) then
         call MPI_SCATTER(cbmf2d_io(1,1),ix*jxp,MPI_DOUBLE_PRECISION,
     &                       cbmf2d(1,1),ix*jxp,MPI_DOUBLE_PRECISION,
     &                                        0,MPI_COMM_WORLD,ierr)
       end if
       if (myid .eq. 0) then
          do j=1,mjx-1
            do l=1,4
              do k=1,kx
                do i=1,ilx
                  sav_1(i,(l-1)*kx+k,j)=absnxt_io(i,k,l,j)
                end do
              end do
            end do
          end do
          allrec = kx*4
          do j=1,mjx-1
            do l=1,kx+1
              do k=1,kx+1
                do i=1,ilx
                  sav_1(i,allrec+(l-1)*(kx+1)+k,j)=abstot_io(i,k,l,j)
                end do
              end do
            end do
          end do
          allrec = allrec+(kx+1)*(kx+1)
          do j=1,mjx-1
            do k=1,kx+1
              do i=1,ilx
                sav_1(i,allrec+k,j)=emstot_io(i,k,j)
              end do
            end do
          end do
          allrec = allrec+kx+1
       end if
       allrec = kx*4+(kx+1)*(kx+2)
       call MPI_SCATTER(
     &          sav_1(1,1,1),ilx*allrec*jxp,MPI_DOUBLE_PRECISION,
     &           sav1(1,1,1),ilx*allrec*jxp,MPI_DOUBLE_PRECISION,
     &                                   0,MPI_COMM_WORLD,ierr)
       do j=1,jendx
         do l=1,4
           do k=1,kx
             do i=1,ilx
               absnxt(i,k,l,j)=sav1(i,(l-1)*kx+k,j)
             end do
           end do
         end do
       end do
       allrec = kx*4
       do j=1,jendx
         do l=1,kx+1
           do k=1,kx+1
             do i=1,ilx
               abstot(i,k,l,j)=sav1(i,allrec+(l-1)*(kx+1)+k,j)
             end do
           end do
         end do
       end do
       allrec = allrec+(kx+1)*(kx+1)
       do j=1,jendx
         do k=1,kx+1
           do i=1,ilx
             emstot(i,k,j)=sav1(i,allrec+k,j)
           end do
         end do
       end do
       if (myid .eq. 0) then
         do j=1,mjx-1
           do n=1,NNSG
             do i=1,ilx
               sav_2(i,       n,j) =  taf2d_io(n,i,j)
               sav_2(i,NNSG  +n,j) = tlef2d_io(n,i,j)
               sav_2(i,NNSG*2+n,j) =  ssw2d_io(n,i,j)
               sav_2(i,NNSG*3+n,j) =  srw2d_io(n,i,j)
             end do
           end do
           do i=1,ilx
             sav_2(i,NNSG*4+1 ,j) =   sol2d_io(i,j)
             sav_2(i,NNSG*4+2 ,j) = solvd2d_io(i,j)
             sav_2(i,NNSG*4+3 ,j) = solvs2d_io(i,j)
             sav_2(i,NNSG*4+4 ,j) =   flw2d_io(i,j)
           end do
         end do
       end if
       allrec = NNSG*4+4
       call MPI_SCATTER(
     &          sav_2(1,1,1),ilx*allrec*jxp,MPI_DOUBLE_PRECISION,
     &           sav2(1,1,1),ilx*allrec*jxp,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
       do j=1,jendx
         do n=1,NNSG
           do i=1,ilx
              taf2d(n,i,j) = sav2(i,       n,j)
             tlef2d(n,i,j) = sav2(i,NNSG  +n,j)
              ssw2d(n,i,j) = sav2(i,NNSG*2+n,j)
              srw2d(n,i,j) = sav2(i,NNSG*3+n,j)
           end do
         end do
         do i=1,ilx
             sol2d(i,j) = sav2(i,NNSG*4+1 ,j)
           solvd2d(i,j) = sav2(i,NNSG*4+2 ,j)
           solvs2d(i,j) = sav2(i,NNSG*4+3 ,j)
             flw2d(i,j) = sav2(i,NNSG*4+4 ,j)
         end do
       end do
       if (myid .eq. 0) then
         do j=1,mjx-1
           do n=1,NNSG
             do i=1,ilx
               sav_2(i,       n,j) =  tgb2d_io(n,i,j)
               sav_2(i,NNSG  +n,j) =  swt2d_io(n,i,j)
               sav_2(i,NNSG*2+n,j) =  scv2d_io(n,i,j)
               sav_2(i,NNSG*3+n,j) = gwet2d_io(n,i,j)
             end do
           end do
           do i=1,ilx
             sav_2(i,NNSG*4+1 ,j) =  flwd2d_io(i,j)
             sav_2(i,NNSG*4+2 ,j) =   fsw2d_io(i,j)
             sav_2(i,NNSG*4+3 ,j) =  sabv2d_io(i,j)
             sav_2(i,NNSG*4+4 ,j) =  sinc2d_io(i,j)
           end do
         end do
       end if
       allrec = NNSG*4+4
       call MPI_SCATTER(
     &          sav_2(1,1,1),ilx*allrec*jxp,MPI_DOUBLE_PRECISION,
     &           sav2(1,1,1),ilx*allrec*jxp,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
       do j=1,jendx
         do n=1,NNSG
           do i=1,ilx
              tgb2d(n,i,j) = sav2(i,       n,j)
              swt2d(n,i,j) = sav2(i,NNSG  +n,j)
              scv2d(n,i,j) = sav2(i,NNSG*2+n,j)
             gwet2d(n,i,j) = sav2(i,NNSG*3+n,j)
           end do
         end do
         do i=1,ilx
            flwd2d(i,j) = sav2(i,NNSG*4+1 ,j)
             fsw2d(i,j) = sav2(i,NNSG*4+2 ,j)
            sabv2d(i,j) = sav2(i,NNSG*4+3 ,j)
            sinc2d(i,j) = sav2(i,NNSG*4+4 ,j)
         end do
       end do
       if (myid .eq. 0) then
         do j=1,mjx-1
           do n=1,NNSG
             do i=1,ilx
               sav_2(i,       n,j) = veg2d1_io(n,i,j)
               sav_2(i,NNSG  +n,j) =  sag2d_io(n,i,j)
               sav_2(i,NNSG*2+n,j) = sice2d_io(n,i,j)
               sav_2(i,NNSG*3+n,j) =  dew2d_io(n,i,j)
             end do
           end do
           do i=1,ilx
             sav_2(i,NNSG*4+1 ,j) =   pptnc_io(i,j)
             sav_2(i,NNSG*4+2 ,j) =    pptc_io(i,j)
             sav_2(i,NNSG*4+3 ,j) =  prca2d_io(i,j)
             sav_2(i,NNSG*4+4 ,j) = prnca2d_io(i,j)
           end do
         end do
       end if
       allrec = NNSG*4+4
       call MPI_SCATTER(
     &          sav_2(1,1,1),ilx*allrec*jxp,MPI_DOUBLE_PRECISION,
     &           sav2(1,1,1),ilx*allrec*jxp,MPI_DOUBLE_PRECISION,
     &                                    0,MPI_COMM_WORLD,ierr)
       do j=1,jendx
         do n=1,NNSG
           do i=1,ilx
             veg2d1(n,i,j) = sav2(i,       n,j)
              sag2d(n,i,j) = sav2(i,NNSG  +n,j)
             sice2d(n,i,j) = sav2(i,NNSG*2+n,j)
              dew2d(n,i,j) = sav2(i,NNSG*3+n,j)
           end do
         end do
         do i=1,ilx
             pptnc(i,j) = sav2(i,NNSG*4+1 ,j)
              pptc(i,j) = sav2(i,NNSG*4+2 ,j)
            prca2d(i,j) = sav2(i,NNSG*4+3 ,j)
           prnca2d(i,j) = sav2(i,NNSG*4+4 ,j)
         end do
       end do
       if (myid .eq. 0) then
         do j=1,mjx-1
           do n=1,NNSG
             do i=1,ilx
               sav_2a(i,       n,j) = ircp2d_io(n,i,j)
               sav_2a(i,NNSG  +n,j) = text2d_io(n,i,j)
               sav_2a(i,NNSG*2+n,j) =  col2d_io(n,i,j)
               sav_2a(i,NNSG*3+n,j) = ocld2d_io(n,i,j)
               sav_2a(i,NNSG*4+n,j) =   tg2d_io(n,i,j)
             end do
           end do
           do i=1,ilx
             sav_2a(i,NNSG*5+1 ,j) =   veg2d_io(i,j)
           end do
         end do
       end if
       allrec = NNSG*5+1
       call MPI_SCATTER(
     &          sav_2a(1,1,1),ilx*allrec*jxp,MPI_DOUBLE_PRECISION,
     &           sav2a(1,1,1),ilx*allrec*jxp,MPI_DOUBLE_PRECISION,
     &                                     0,MPI_COMM_WORLD,ierr)
       do j=1,jendx
         do n=1,NNSG
           do i=1,ilx
             ircp2d(n,i,j) = sav2a(i,       n,j)
             text2d(n,i,j) = sav2a(i,NNSG  +n,j)
              col2d(n,i,j) = sav2a(i,NNSG*2+n,j)
             ocld2d(n,i,j) = sav2a(i,NNSG*3+n,j)
               tg2d(n,i,j) = sav2a(i,NNSG*4+n,j)
           end do
         end do
         do i=1,ilx
             veg2d(i,j) = sav2a(i,NNSG*5+1 ,j)
         end do
       end do
       if(ichem.eq.1) then
         if (myid .eq. 0) then
           do j=1,mjx
             do n=1,ntr
               do k=1,kx
                 do i=1,ix
                   sav_4(i,         (n-1)*kx+k,j)=  chia_io(i,k,j,n)
                   sav_4(i,ntr*kx+  (n-1)*kx+k,j)=  chib_io(i,k,j,n)
                   sav_4(i,ntr*kx*2+(n-1)*kx+k,j)=remlsc_io(i,k,j,n)
                   sav_4(i,ntr*kx*3+(n-1)*kx+k,j)=remcvc_io(i,k,j,n)
                 end do
               end do
             end do
           end do
           allrec = 4*ntr*kx
           do j=1,mjx
             do n=1,ntr
               do i=1,ix
                 sav_4(i,allrec+n,j)=remdrd_io(i,j,n)
               end do
             end do
           end do
           allrec = allrec+ntr
         end if
         allrec = ntr*(kx*4+1)
         call MPI_SCATTER(
     &            sav_4(1,1,1),ix*allrec*jxp,MPI_DOUBLE_PRECISION,
     &             sav4(1,1,1),ix*allrec*jxp,MPI_DOUBLE_PRECISION,
     &                                     0,MPI_COMM_WORLD,ierr)
         do j=1,jendl
           do n=1,ntr
             do k=1,kx
               do i=1,ix
                   chia(i,k,j,n)=sav4(i,         (n-1)*kx+k,j)
                   chib(i,k,j,n)=sav4(i,ntr*kx+  (n-1)*kx+k,j)
                 remlsc(i,k,j,n)=sav4(i,ntr*kx*2+(n-1)*kx+k,j)
                 remcvc(i,k,j,n)=sav4(i,ntr*kx*3+(n-1)*kx+k,j)
               end do
             end do
           end do
         end do
         allrec = 4*ntr*kx
         do j=1,jendl
           do n=1,ntr
             do i=1,ix
               remdrd(i,j,n)=sav4(i,allrec+n,j)
             end do
           end do
         end do
         if(myid.eq.0) then
           do j=1,mjx-1
             do i=1,ilx
               sav_4a(i,1,j) = ssw2da_io(i,j)
               sav_4a(i,2,j) = sdeltk2d_io(i,j)
               sav_4a(i,3,j) = sdelqk2d_io(i,j)
               sav_4a(i,4,j) = sfracv2d_io(i,j)
               sav_4a(i,5,j) = sfracb2d_io(i,j)
               sav_4a(i,6,j) = sfracs2d_io(i,j)
               sav_4a(i,7,j) = svegfrac2d_io(i,j)
             end do
           end do
         end if
         call MPI_SCATTER(
     &          sav_4a(1,1,1),ilx*7*jxp,MPI_DOUBLE_PRECISION,
     &           sav4a(1,1,1),ilx*7*jxp,MPI_DOUBLE_PRECISION,
     &                                0,MPI_COMM_WORLD,ierr)
         do j=1,jendx
           do i=1,ilx
             ssw2da(i,j)     = sav4a(i,1,j)
             sdeltk2d(i,j)   = sav4a(i,2,j)
             sdelqk2d(i,j)   = sav4a(i,3,j)
             sfracv2d(i,j)   = sav4a(i,4,j)
             sfracb2d(i,j)   = sav4a(i,5,j)
             sfracs2d(i,j)   = sav4a(i,6,j)
             svegfrac2d(i,j) = sav4a(i,7,j)
           end do
         end do
       end if
#endif
       dt = dt2
#ifdef MPP1
       call MPI_BCAST(mdate0,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(jyear0,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(ktau,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(jyear,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(xtime,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(ldatez,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(lyear,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(lmonth,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(lday,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(lhour,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(ntime,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(jyearr,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(ktaur,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
#ifdef DIAG
      call MPI_BCAST(tdini,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(tdadv,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(tqini,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(tqadv,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(tqeva,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(tqrai,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
       if(ichem.eq.1) then
         call MPI_BCAST(tchiad,ntr,MPI_DOUBLE_PRECISION,
     &                  0,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(tchitb,ntr,MPI_DOUBLE_PRECISION,
     &                  0,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(tchie,ntr,MPI_DOUBLE_PRECISION,
     &                  0,MPI_COMM_WORLD,ierr)
       end if
#endif
c------lake model
       if(lakemod.eq.1)then
         call MPI_BCAST(ilake,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(jlake,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(depth,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(hs,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
        call MPI_BCAST(eta,1,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(tlake,depth,MPI_DOUBLE_PRECISION
     &                 ,0,MPI_COMM_WORLD,ierr)
       end if
        call MPI_SENDRECV(ht(1,jxp),ix,MPI_DOUBLE_PRECISION,ieast,1,
     &                    ht(1,0),ix,MPI_DOUBLE_PRECISION,iwest,1,
     &                    MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(ht(1,1),ix,MPI_DOUBLE_PRECISION,iwest,2,
     &                    ht(1,jxp+1),ix,MPI_DOUBLE_PRECISION,ieast,2,
     &                    MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(msfx(1,jxp-1),ix*2,MPI_DOUBLE_PRECISION,
     &            ieast,1,msfx(1,-1),ix*2,MPI_DOUBLE_PRECISION,
     &            iwest,1,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(msfx(1,1),ix*2,MPI_DOUBLE_PRECISION,
     &            iwest,2,msfx(1,jxp+1),ix*2,MPI_DOUBLE_PRECISION,
     &            ieast,2,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(msfd(1,jxp-1),ix*2,MPI_DOUBLE_PRECISION,
     &            ieast,1,msfd(1,-1),ix*2,MPI_DOUBLE_PRECISION,
     &            iwest,1,MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(msfd(1,1),ix*2,MPI_DOUBLE_PRECISION,
     &            iwest,2,msfd(1,jxp+1),ix*2,MPI_DOUBLE_PRECISION,
     &            ieast,2,MPI_COMM_WORLD,status,ierr)
#endif
c
c-----end of initial/restart if test
      end if
! Move from param.F to fix the reatart problem found by
!      Zhang DongFeng
      if(ipptls.eq.1) then
#ifdef MPP1
        do j=1,jendx
#else
        do j=1,jlx
#endif
        do i=1,ilx
          if (satbrt(i,j).gt.13.9 .and. satbrt(i,j).lt.15.1) then
            qck1(i,j) = qck1oce    ! OCEAN
            cgul(i,j) = guloce     ! OCEAN
            rh0(i,j)  = rh0oce     ! OCEAN
          else
            qck1(i,j) = qck1land   ! LAND
            cgul(i,j) = gulland    ! LAND
            rh0(i,j)  = rh0land    ! LAND
          end if
        end do
        end do
      end if
Cchem2
      if(ichem.eq.1) then
        iso2 = 0
        iso4 = 0
        ibchl = 0
        ibchb = 0
        iochl = 0
        iochb = 0
        ibin = 0
        do itr =1,ntr
          if  (chtrname(itr).eq. 'SO2') iso2 = itr
          if  (chtrname(itr).eq. 'SO4') iso4 = itr
          if  (chtrname(itr).eq. 'BC_HL') ibchl = itr
          if  (chtrname(itr).eq. 'BC_HB') ibchb = itr
          if  (chtrname(itr).eq. 'OC_HL') iochl = itr
          if  (chtrname(itr).eq. 'OC_HB') iochb = itr
          if  (chtrname(itr).eq. 'DUST') then 
            ibin = ibin + 1
            idust(ibin) = itr 
          end if      
        end do
      end if
Cchem2_
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c ****** initialize and define constants for vector bats

      call bconst
      if (jyear.eq.jyear0.and.ktau.eq.0) call initb
      if(iemiss.eq.1) then
#ifdef MPP1
         do j=1,jendx
#else
         do j=1,jx-1
#endif
         do i=1,ix-1
         do n=1,NNSG
            ist = nint(veg2d1(n,i,j))
            if(ist.eq.0) then
               emiss2d(n,i,j) = 0.955d0
            else if(ist.eq.8) then
               emiss2d(n,i,j) = 0.76d0
            else if(ist.eq.11) then
               emiss2d(n,i,j) = 0.85d0
            else if(ist.eq.12) then
               emiss2d(n,i,j) = 0.97d0
            else
               emiss2d(n,i,j) = 0.99d0-(albvgs(ist)+albvgl(ist))*0.1d0
            end if
!           emiss2d(n,i,j) = 1.0d0
         end do
         end do
         end do
      end if
c ::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
c
c-----read in the boundary conditions for large domain:
c
c-----compute the solar declination angle:
c
      call solar1(xtime)
      call inirad
c
c-----calculating topographical correction to diffusion coefficient
#ifdef MPP1
      do 634 j=1,jendl
#else
      do 634 j=1,jx
#endif
      do 634 i=1,ix
       hgfact(i,j)=1.
634   continue
#ifdef MPP1
      do 635 j=jbegin,jendm
       if (myid .eq. 0) then
         jm1h=max0(j-1,2)
       else
         jm1h=j-1
       end if
       if (myid .eq. nproc-1) then
         jp1h=min0(j+1,jxp-2)
       else
         jp1h=j+1
       end if
#else
      do 635 j=2,jx-2
       jm1h=max0(j-1,2)
       jp1h=min0(j+1,jx-2)
#endif
      do 635 i=2,ix-2
       im1h=max0(i-1,2)
       ip1h=min0(i+1,ix-2)
       hg1=dabs((ht(i,j)-ht(im1h,j))/dx)
       hg2=dabs((ht(i,j)-ht(ip1h,j))/dx)
       hg3=dabs((ht(i,j)-ht(i,jm1h))/dx)
       hg4=dabs((ht(i,j)-ht(i,jp1h))/dx)
       hgmax=dmax1(hg1,hg2,hg3,hg4)/9.8
       hgfact(i,j)=1./(1.+(hgmax/0.001)**2.)
635   continue
c
c-----set up output time:
c
      icnt = 11        ! set counter for safety-file moves on cycad
      dectim= anint (xtime + dectim)
#ifdef MPP1
      if (myid .eq. 0) then
#endif
      write(*,*) 'dectim = ',dectim
#ifdef MPP1
      end if
#endif
c
      return
      end subroutine init

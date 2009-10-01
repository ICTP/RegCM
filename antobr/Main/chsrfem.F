      subroutine chsrfem
      implicit none
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
      include 'Commons/param1.cb'
      include 'Commons/param2.cb'
      include 'Commons/param3.cb'
#ifdef MPP1
      include 'Common2/main.cb'
      include 'Common2/mainchem.cb'
      include 'Common2/trachem.cb'
      include 'Common2/dust.cb'
      include 'Common2/io.cb'
#else
      include 'Commons/main.cb'
      include 'Commons/mainchem.cb'
      include 'Commons/trachem.cb'
      include 'Commons/dust.cb'
#endif
      include 'Commons/iunits.cb'
      include 'Commons/aeroppt.cb'

#ifdef MPP1
      include 'mpiregcm.h'
      real*8  src0(ix,12,ntr,jxp),src_0(ix,12,ntr,mjx)
      COMMON /tmpSRC/ src0,src_0
      real*4  toto(ix,mjx)
#else
      real*4  toto(ix,jx)
#endif
      integer i,j,k,l,m,n,itr
      character*5 AERCTL
      logical there
      logical rd_tex

      rd_tex = .false.
#ifdef MPP1
      if (myid .eq. 0) then
#endif
      do itr=1, ntr
        AERCTL = chtrname(itr)
        if(AERCTL(1:4).eq.'DUST') then
          rd_tex = .true.
          goto 1000
        endif
      enddo
 1000 continue
#ifdef MPP1
      endif
#endif
#ifdef MPP1
      call MPI_BCAST(rd_tex,1,MPI_LOGICAL,0,MPI_COMM_WORLD,ierr)
      if (myid .eq. 0) then
      do itr=1, ntr
         AERCTL = chtrname(itr)
         print *,itr, AERCTL
         if(AERCTL(1:4).ne.'DUST'.and.AERTYP(4:5).ne.'00')then
            open(unit=iutchsrc,file='AERO.dat',status='old'
     &          ,form='unformatted',access='direct',recl=ix*mjx*IBYTE)
            if(AERCTL(1:3).eq.'SO2') then
               if(AERTYP(4:4).eq.'1') then
                  read(iutchsrc,rec=1) ((toto(i,j),j=1,mjx),i=1,ix)
                  do m=1,12
                    do j=1,mjx
                    do i=1,ix
                      chemsrc_io(i,j,m,itr) = toto(i,j)
                    enddo
                    enddo
                  enddo
               else
                 do m=1,12
                   do j=1,mjx
                   do i=1,ix
                      chemsrc_io(i,j,m,itr) = 0.0d0
                    enddo
                    enddo
                  enddo
               endif
               if(AERTYP(5:5).eq.'1') then
                  do m=1,12
                    read(iutchsrc,rec=3+m) ((toto(i,j),j=1,mjx),i=1,ix)
                    do j=1,mjx
                    do i=1,ix
                      chemsrc_io(i,j,m,itr) = chemsrc_io(i,j,m,itr)
     &                                       +toto(i,j)
                    enddo
                    enddo
                  enddo
               endif
            else if(AERCTL(1:2).eq.'BC') then
               if(AERTYP(4:4).eq.'1') then
                  read(iutchsrc,rec=2) ((toto(i,j),j=1,mjx),i=1,ix)
                  do m=1,12
                    do j=1,mjx
                    do i=1,ix
                      chemsrc_io(i,j,m,itr) = toto(i,j)
                    enddo
                    enddo
                  enddo
               else
                 do m=1,12
                   do j=1,mjx
                   do i=1,ix
                      chemsrc_io(i,j,m,itr) = 0.0d0
                    enddo
                    enddo
                  enddo
               endif
               if(AERTYP(5:5).eq.'1') then
                  do m=1,12
                    read(iutchsrc,rec=15+m) ((toto(i,j),j=1,mjx),i=1,ix)
                    do j=1,mjx
                    do i=1,ix
                      chemsrc_io(i,j,m,itr) = chemsrc_io(i,j,m,itr)
     &                                       +toto(i,j)
                    enddo
                    enddo
                  enddo
               endif
            else if(AERCTL(1:2).eq.'OC') then
               if(AERTYP(4:4).eq.'1') then
                  read(iutchsrc,rec=3) ((toto(i,j),j=1,mjx),i=1,ix)
                  do m=1,12
                    do j=1,mjx
                    do i=1,ix
                      chemsrc_io(i,j,m,itr) = toto(i,j)
                    enddo
                    enddo
                  enddo
               else
                 do m=1,12
                   do j=1,mjx
                   do i=1,ix
                      chemsrc_io(i,j,m,itr) = 0.0d0
                    enddo
                    enddo
                  enddo
               endif
               if(AERTYP(5:5).eq.'1') then
                  do m=1,12
                    read(iutchsrc,rec=27+m) ((toto(i,j),j=1,mjx),i=1,ix)
                    do j=1,mjx
                    do i=1,ix
                      chemsrc_io(i,j,m,itr) = chemsrc_io(i,j,m,itr)
     &                                       +toto(i,j)
                    enddo
                    enddo
                  enddo
               endif
            endif
         endif
      enddo
      do j=1,mjx
        do itr=1,ntr
          do m=1,12
            do i=1,ix
              src_0(i,m,itr,j)=chemsrc_io(i,j,m,itr)
            enddo
          enddo
        enddo
      enddo
      endif
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
#else
      do itr=1, ntr
         AERCTL = chtrname(itr)
         print *,itr, AERCTL
         if(AERCTL(1:4).ne.'DUST'.and.AERTYP(4:5).ne.'00')then
            open(unit=iutchsrc,file='AERO.dat',status='old'
     &          ,form='unformatted',access='direct',recl=ix*jx*IBYTE)
            if(AERCTL(1:3).eq.'SO2') then
               if(AERTYP(4:4).eq.'1') then
                  read(iutchsrc,rec=1) ((toto(i,j),j=1,jx),i=1,ix)
                  do m=1,12
                    do j=1,jx
                    do i=1,ix
                      chemsrc(i,j,m,itr) = toto(i,j)
                    enddo
                    enddo
                  enddo
               else
                 do m=1,12
                   do j=1,jx
                   do i=1,ix
                      chemsrc(i,j,m,itr) = 0.0d0
                    enddo
                    enddo
                  enddo
               endif
               if(AERTYP(5:5).eq.'1') then
                  do m=1,12
                    read(iutchsrc,rec=3+m) ((toto(i,j),j=1,jx),i=1,ix)
                    do j=1,jx
                    do i=1,ix
                      chemsrc(i,j,m,itr) = chemsrc(i,j,m,itr)+toto(i,j)
                    enddo
                    enddo
                  enddo
               endif
            else if(AERCTL(1:2).eq.'BC') then
               if(AERTYP(4:4).eq.'1') then
                  read(iutchsrc,rec=2) ((toto(i,j),j=1,jx),i=1,ix)
                  do m=1,12
                    do j=1,jx
                    do i=1,ix
                      chemsrc(i,j,m,itr) = toto(i,j)
                    enddo
                    enddo
                  enddo
               else
                 do m=1,12
                   do j=1,jx
                   do i=1,ix
                      chemsrc(i,j,m,itr) = 0.0d0
                    enddo
                    enddo
                  enddo
               endif
               if(AERTYP(5:5).eq.'1') then
                  do m=1,12
                    read(iutchsrc,rec=15+m) ((toto(i,j),j=1,jx),i=1,ix)
                    do j=1,jx
                    do i=1,ix
                      chemsrc(i,j,m,itr) = chemsrc(i,j,m,itr)+toto(i,j)
                    enddo
                    enddo
                  enddo
               endif
            else if(AERCTL(1:2).eq.'OC') then
               if(AERTYP(4:4).eq.'1') then
                  read(iutchsrc,rec=3) ((toto(i,j),j=1,jx),i=1,ix)
                  do m=1,12
                    do j=1,jx
                    do i=1,ix
                      chemsrc(i,j,m,itr) = toto(i,j)
                    enddo
                    enddo
                  enddo
               else
                 do m=1,12
                   do j=1,jx
                   do i=1,ix
                      chemsrc(i,j,m,itr) = 0.0d0
                    enddo
                    enddo
                  enddo
               endif
               if(AERTYP(5:5).eq.'1') then
                  do m=1,12
                    read(iutchsrc,rec=27+m) ((toto(i,j),j=1,jx),i=1,ix)
                    do j=1,jx
                    do i=1,ix
                      chemsrc(i,j,m,itr) = chemsrc(i,j,m,itr)+toto(i,j)
                    enddo
                    enddo
                  enddo
               endif
            endif
         endif
      enddo
#endif
#ifdef MPP1
      call MPI_SCATTER(src_0(1,1,1,1),ix*12*ntr*jxp,MPI_DOUBLE_PRECISION
     &                 ,src0(1,1,1,1),ix*12*ntr*jxp,MPI_DOUBLE_PRECISION
     &                                           ,0,MPI_COMM_WORLD,ierr)
      do j=1,jendl
        do itr=1,ntr
          do m=1,12
            do i=1,ix
              chemsrc(i,j,m,itr)=src0(i,m,itr,j)
            enddo
          enddo
        enddo
      enddo
#endif

c modification dust : read the soil texture type
#ifdef MPP1
      if (myid .eq. 0) then
      if(rd_tex) then
         if(LSMTYP.eq.'BATS') then
            read(iutin,rec=14) ((toto(i,j),j=1,mjx),i=1,ix)
         else if(LSMTYP.eq.'USGS') then
            read(iutin,rec=43) ((toto(i,j),j=1,mjx),i=1,ix)
         endif
         close(iutin)
         do j=1,mjx
         do i=1,ix
           dustsotex_io(i,j) = dble(toto(i,j))
         enddo
         enddo
      endif
      endif
c     call MPI_BARRIER(MPI_COMM_WORLD,ierr)
      call MPI_SCATTER(dustsotex_io(1,1),ix*jxp,MPI_DOUBLE_PRECISION
     &                   ,dustsotex(1,1),ix*jxp,MPI_DOUBLE_PRECISION
     &                                       ,0,MPI_COMM_WORLD,ierr)
#else
      if(rd_tex) then
         if(LSMTYP.eq.'BATS') then
            read(iutin,rec=14) ((toto(i,j),j=1,jx),i=1,ix)
         else if(LSMTYP.eq.'USGS') then
            read(iutin,rec=43) ((toto(i,j),j=1,jx),i=1,ix)
         endif
         close(iutin)
         do j=1,jx
         do i=1,ix
           dustsotex(i,j) = dble(toto(i,j))
         enddo
         enddo
      endif
#endif
      if(rd_tex) then
         call inidust
      endif

c    sulfates sources
      do m=1,12
#ifdef MPP1
         do j=1,jendl
#else
         do j=1,jx
#endif
         do i=1,ix
            if(iso4.gt.0) chemsrc(i,j,m,iso4)=0.02*chemsrc(i,j,m,iso2)
            if(iso2.gt.0) chemsrc(i,j,m,iso2)=0.98*chemsrc(i,j,m,iso2)

c    partition hydrophilic hydrophonic ( cooke et al.1999)
c    BC
            if(ibchb.gt.0 .and. ibchl.gt.0) then
               chemsrc(i,j,m,ibchl) = 0.2 * chemsrc(i,j,m,ibchb)
               chemsrc(i,j,m,ibchb) = 0.8 * chemsrc(i,j,m,ibchb)
            endif
c    OC
            if(iochb.gt.0 .and. iochl.gt.0) then
               chemsrc(i,j,m,iochl) = 0.5 * chemsrc(i,j,m,iochb)
               chemsrc(i,j,m,iochb) = 0.5 * chemsrc(i,j,m,iochb)
            endif
         enddo
         enddo
      enddo

c    OPtical properties / internal mixing
      if(mixtype.eq.2) then
#ifdef MPP1
      if (myid .eq. 0) then
#endif
         inquire(file='optdat.bin',exist=there)
         if(.not.there) then
            write(*,*) 'For mixtype=2, optdat.bin is required'
            write(*,*) 'ln -s ../Main/Commons/optdat.bin optdat.bin'
#ifdef MPP1
            write(*,*) 'stop in chsrfem.f'
            call MPI_ABORT(MPI_COMM_WORLD,ierr)
#else
            stop
#endif
         endif
         open(iutopt,file='optdat.bin',form='unformatted'
     &              ,recl=4*19*11*11*11*11*IBYTE,access='direct')
         read(iutopt,rec=1) ((((((dextmix(i,j,k,l,m,n)
     &              ,i=1,4),j=1,19),k=1,11),l=1,11),m=1,11),n=1,11)
         read(iutopt,rec=2) ((((((dssamix(i,j,k,l,m,n)
     &              ,i=1,4),j=1,19),k=1,11),l=1,11),m=1,11),n=1,11)
         read(iutopt,rec=3) ((((((dgmix(i,j,k,l,m,n)
     &              ,i=1,4),j=1,19),k=1,11),l=1,11),m=1,11),n=1,11)
         close(iutopt)
#ifdef MPP1
      endif
         call MPI_BCAST(dextmix,4*19*11*11*11*11,MPI_REAL,0
     &       ,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(dssamix,4*19*11*11*11*11,MPI_REAL,0
     &       ,MPI_COMM_WORLD,ierr)
         call MPI_BCAST(dgmix,4*19*11*11*11*11,MPI_REAL,0
     &       ,MPI_COMM_WORLD,ierr)
#endif

c       Check !
        
         do k=1,11
         do l=1,11 
         do m=1,11
         do n=1,11
  
         if (k+l+m+n.eq.14) then
         do i=1,4
         do j=1,19

            if ((dextmix(i,j,k,l,m,n) .lt. 0.) .or.  
     &          (dextmix(i,j,k,l,m,n) .gt. 20.)) then
               print*,'problem in dextmix ', dextmix(i,j,k,l,m,n)
               stop
            endif
    
            if ((dssamix(i,j,k,l,m,n) .lt. 0.) .or.  
     &          (dssamix(i,j,k,l,m,n) .gt. 1.)) then
               print*,'problem in dextmix ', dssamix(i,j,k,l,m,n)
               stop
            endif
     
            if ((dgmix(i,j,k,l,m,n) .lt. 0.) .or.  
     &          (dgmix(i,j,k,l,m,n) .gt. 1.)) then
               print*,'problem in dextmix ', dgmix(i,j,k,l,m,n)
               stop
            endif
     
         end do
         end do
         endif

         end do
         end do
         end do
         end do

#ifdef MPP1
         if(myid .eq. 0) write(*,*) '! OPDATA CHECKED !'
#else
         write(*,*) '! OPDATA CHECKED !'
#endif
      end if
      return
      end

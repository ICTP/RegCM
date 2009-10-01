cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine tracbud(mc)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
c
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
      include 'Commons/param1.cb'
#ifdef DIAG
      include 'Commons/diagnosis.cb'
#endif
      include 'Commons/param2.cb'
      include 'Commons/param3.cb'
#ifdef MPP1
      include 'Common2/main.cb'
      include 'Common2/mainchem.cb'
      include 'Common2/trachem.cb'
      include 'Common2/io.cb'
#else
      include 'Commons/main.cb'
      include 'Commons/mainchem.cb'
      include 'Commons/trachem.cb'
#endif
      include 'Commons/date.cb'
#ifdef MPP1
      include 'mpiregcm.h'
      integer l
#endif
c
      integer mc
      real*8  tttmp
      integer i,j,k,itr
c
c---- add tracers,  tremlsc, tremcvc, tremdrd are total amount mass
c and chemical conversion term
c   removed by large scale, convective precipation and dry deposition
c   respectively

      do 5 itr = 1,ntr
C
#ifdef DIAG
       ttrace(itr,1) = 0.
       tremlsc(itr,1)= 0.
       tremcvc(itr,1)= 0.
       trxsg(itr,1) = 0.
       trxsaq1(itr,1) = 0.
       trxsaq2(itr,1) = 0.
       tremdrd(itr,1)= 0.       
#endif
C
#ifdef MPP1
      do 5 j=1,jendx
#else
      do 5 j=1,jx-1
#endif
      do 5 i=1,ix-1
       dtrace(i,j,itr)=0.0
       wdlsc(i,j,itr)=0.0
       wdcvc(i,j,itr)=0.0
       wxsg(i,j,itr)=0.0
       wxaq(i,j,itr)=0.0
       ddsfc(i,j,itr)=0.0
5     continue
      
           
c-----tracers (unit = kg):
      do 24 itr=1,ntr
#ifdef MPP1
      do 24 j=1,jendx
#else
      do 24 j=1,jx-1
#endif
      do 24 i=1,ix-1
      do  k=1,kx
       dtrace(i,j,itr)=dtrace(i,j,itr)+chia(i,k,j,itr)*dsigma(k)
       wdlsc(i,j,itr)=wdlsc(i,j,itr)+remlsc(i,k,j,itr)*dsigma(k)
       wdcvc(i,j,itr)=wdcvc(i,j,itr)+remcvc(i,k,j,itr)*dsigma(k)
       wxsg(i,j,itr) = wxsg(i,j,itr)+ rxsg(i,k,j,itr)*dsigma(k)
c     sum ls and conv contribution
       wxaq(i,j,itr) = wxaq(i,j,itr)
     &               + (rxsaq1(i,k,j,itr)+rxsaq2(i,k,j,itr))*dsigma(k)
      enddo
      ddsfc(i,j,itr)=ddsfc(i,j,itr)+remdrd(i,j,itr)*dsigma(kx)
c Source cumulated diag(care the unit are alredy .m-2)
      cemtrac(i,j,itr)=cemtr(i,j,itr)
24    continue

#ifdef DIAG
#ifdef MPP1
      do itr=1,ntr
      call MPI_GATHER(chia(1,1,1,itr),ix*kx*jxp,MPI_DOUBLE_PRECISION,
     &             chia_io(1,1,1,itr),ix*kx*jxp,MPI_DOUBLE_PRECISION,
     &                                        0,MPI_COMM_WORLD,ierr)
      call MPI_GATHER(remlsc(1,1,1,itr),ix*kx*jxp,MPI_DOUBLE_PRECISION,
     &             remlsc_io(1,1,1,itr),ix*kx*jxp,MPI_DOUBLE_PRECISION,
     &                                          0,MPI_COMM_WORLD,ierr)
      call MPI_GATHER(remcvc(1,1,1,itr),ix*kx*jxp,MPI_DOUBLE_PRECISION,
     &             remcvc_io(1,1,1,itr),ix*kx*jxp,MPI_DOUBLE_PRECISION,
     &                                          0,MPI_COMM_WORLD,ierr)
      call MPI_GATHER(rxsg(1,1,1,itr),ix*kx*jxp,MPI_DOUBLE_PRECISION,
     &             rxsg_io(1,1,1,itr),ix*kx*jxp,MPI_DOUBLE_PRECISION,
     &                                        0,MPI_COMM_WORLD,ierr)
      call MPI_GATHER(rxsaq1(1,1,1,itr),ix*kx*jxp,MPI_DOUBLE_PRECISION,
     &             rxsaq1_io(1,1,1,itr),ix*kx*jxp,MPI_DOUBLE_PRECISION,
     &                                          0,MPI_COMM_WORLD,ierr)
      call MPI_GATHER(rxsaq2(1,1,1,itr),ix*kx*jxp,MPI_DOUBLE_PRECISION,
     &             rxsaq2_io(1,1,1,itr),ix*kx*jxp,MPI_DOUBLE_PRECISION,
     &                                          0,MPI_COMM_WORLD,ierr)
      call MPI_GATHER(remdrd(1,1,itr),ix*jxp,MPI_DOUBLE_PRECISION,
     &             remdrd_io(1,1,itr),ix*jxp,MPI_DOUBLE_PRECISION,
     &                                       0,MPI_COMM_WORLD,ierr)
      do l=1,12
      call MPI_GATHER(chemsrc(1,1,l,itr),ix*jxp,MPI_DOUBLE_PRECISION,
     &             chemsrc_io(1,1,l,itr),ix*jxp,MPI_DOUBLE_PRECISION,
     &                                       0,MPI_COMM_WORLD,ierr)
      enddo
      enddo
      if(myid.eq.0) then
      do itr=1,ntr
        do k=1,kx
          tttmp= 0.
          do j=1+1,mjx-1-1
          do i=1+1,ix-1-1
            tttmp=tttmp+chia_io(i,k,j,itr)
          enddo
          enddo
          ttrace(itr,1)=ttrace(itr,1)+tttmp*dsigma(k)
          tttmp= 0.
          do j=1+1,mjx-1-1
          do i=1+1,ix-1-1
            tttmp=tttmp+remlsc_io(i,k,j,itr)
          enddo
          enddo
          tremlsc(itr,1)=tremlsc(itr,1)+tttmp*dsigma(k)
          tttmp= 0.
          do j=1+1,mjx-1-1
          do i=1+1,ix-1-1
            tttmp=tttmp+remcvc_io(i,k,j,itr)
          enddo
          enddo
          tremcvc(itr,1)=tremcvc(itr,1)+tttmp*dsigma(k)
          tttmp= 0.
          do j=1+1,mjx-1-1
          do i=1+1,ix-1-1
            tttmp=tttmp+rxsg_io(i,k,j,itr)
          enddo
          enddo
          trxsg(itr,1)=trxsg(itr,1)+tttmp*dsigma(k)
          tttmp= 0.
          do j=1+1,mjx-1-1
          do i=1+1,ix-1-1
            tttmp=tttmp+rxsaq1_io(i,k,j,itr)
          enddo
          enddo
          trxsaq1(itr,1)=trxsaq1(itr,1)+tttmp*dsigma(k)
          tttmp= 0.
          do j=1+1,mjx-1-1
          do i=1+1,ix-1-1
            tttmp=tttmp+rxsaq2_io(i,k,j,itr)
          enddo
          enddo
          trxsaq2(itr,1)=trxsaq2(itr,1)+tttmp*dsigma(k)
        enddo
      enddo

      do itr=1,ntr
        ttrace(itr,1)=ttrace(itr,1)*dx*dx
        tremlsc(itr,1)=tremlsc(itr,1)*dx*dx
        tremcvc(itr,1)=tremcvc(itr,1)*dx*dx
        trxsg(itr,1)=trxsg(itr,1)*dx*dx
        trxsaq1(itr,1)=trxsaq1(itr,1)*dx*dx
        trxsaq2(itr,1)=trxsaq2(itr,1)*dx*dx
      enddo

      do itr=1,ntr
        tttmp= 0.
        do j=1+1,mjx-1-1
        do i=1+1,ix-1-1
          tttmp=tttmp+remdrd_io(i,j,itr)
        end do
        end do
        tremdrd(itr,1)=tremdrd(itr,1)+tttmp*dx*dx*dsigma(kx)
 
c       emissions
        tttmp= 0.
        do j=1+1,mjx-1-1
        do i=1+1,ix-1-1
          tttmp=tttmp+chemsrc_io(i,j,lmonth,itr)*dtmin*60.*dx*dx
        end do
        end do
        tchie(itr)= tchie(itr) + tttmp
      enddo
      endif
      call MPI_BCAST(ttrace(1,1),ntr,MPI_DOUBLE_PRECISION,
     &                             0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(tremlsc(1,1),ntr,MPI_DOUBLE_PRECISION,
     &                              0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(tremcvc(1,1),ntr,MPI_DOUBLE_PRECISION,
     &                              0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(trxsg(1,1),ntr,MPI_DOUBLE_PRECISION,
     &                            0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(trxsaq1(1,1),ntr,MPI_DOUBLE_PRECISION,
     &                              0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(trxsaq2(1,1),ntr,MPI_DOUBLE_PRECISION,
     &                              0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(tremdrd(1,1),ntr,MPI_DOUBLE_PRECISION,
     &                              0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(tchie(1),ntr,MPI_DOUBLE_PRECISION,
     &                          0,MPI_COMM_WORLD,ierr)
#else
      do itr=1,ntr
        do k=1,kx
          tttmp= 0.
          do j=1+1,jx-1-1
          do i=1+1,ix-1-1
            tttmp=tttmp+chia(i,k,j,itr)
          enddo
          enddo
          ttrace(itr,1)=ttrace(itr,1)+tttmp*dsigma(k)
          tttmp= 0.
          do j=1+1,jx-1-1
          do i=1+1,ix-1-1
            tttmp=tttmp+remlsc(i,k,j,itr)
          enddo
          enddo
          tremlsc(itr,1)=tremlsc(itr,1)+tttmp*dsigma(k)
          tttmp= 0.
          do j=1+1,jx-1-1
          do i=1+1,ix-1-1
            tttmp=tttmp+remcvc(i,k,j,itr)
          enddo
          enddo
          tremcvc(itr,1)=tremcvc(itr,1)+tttmp*dsigma(k)
          tttmp= 0.
          do j=1+1,jx-1-1
          do i=1+1,ix-1-1
            tttmp=tttmp+rxsg(i,k,j,itr)
          enddo
          enddo
          trxsg(itr,1)=trxsg(itr,1)+tttmp*dsigma(k)
          tttmp= 0.
          do j=1+1,jx-1-1
          do i=1+1,ix-1-1
            tttmp=tttmp+rxsaq1(i,k,j,itr)
          enddo
          enddo
          trxsaq1(itr,1)=trxsaq1(itr,1)+tttmp*dsigma(k)
          tttmp= 0.
          do j=1+1,jx-1-1
          do i=1+1,ix-1-1
            tttmp=tttmp+rxsaq2(i,k,j,itr)
          enddo
          enddo
          trxsaq2(itr,1)=trxsaq2(itr,1)+tttmp*dsigma(k)
        enddo
      enddo

      do itr=1,ntr
        ttrace(itr,1)=ttrace(itr,1)*dx*dx
        tremlsc(itr,1)=tremlsc(itr,1)*dx*dx
        tremcvc(itr,1)=tremcvc(itr,1)*dx*dx
        trxsg(itr,1)=trxsg(itr,1)*dx*dx
        trxsaq1(itr,1)=trxsaq1(itr,1)*dx*dx
        trxsaq2(itr,1)=trxsaq2(itr,1)*dx*dx
      enddo

      do itr=1,ntr
        tttmp= 0.
        do j=1+1,jx-1-1
        do i=1+1,ix-1-1
          tttmp=tttmp+remdrd(i,j,itr)
        end do
        end do
        tremdrd(itr,1)=tremdrd(itr,1)+tttmp*dx*dx*dsigma(kx)
 
c       emissions
        tttmp= 0.
        do j=1+1,jx-1-1
        do i=1+1,ix-1-1
          tttmp=tttmp+chemsrc(i,j,lmonth,itr)*dtmin*60.*dx*dx
        end do
        end do
        tchie(itr)= tchie(itr) + tttmp
      enddo
#endif
#endif

      do 40 itr=1,ntr
#ifdef DIAG
       ttrace(itr,1)=ttrace(itr,1)*1000./g
       tremlsc(itr,1)=tremlsc(itr,1)*1000./g
       tremcvc(itr,1)=tremcvc(itr,1)*1000./g
       tremdrd(itr,1)=tremdrd(itr,1)*1000./g
       trxsg(itr,1)=trxsg(itr,1)*1000./g
       trxsaq1(itr,1)=trxsaq1(itr,1)*1000./g
       trxsaq2(itr,1)=trxsaq2(itr,1)*1000./g
#endif
c     
#ifdef MPP1
      do 40 j=1,jendx
#else
      do 40 j=1,jx-1
#endif
      do 40 i=1,ix-1
       dtrace(i,j,itr)=1.e6*dtrace(i,j,itr)*1000./g     ! unit: mg/m2
       wdlsc(i,j,itr)=1.e6*wdlsc(i,j,itr)*1000./g
       wdcvc(i,j,itr)=1.e6*wdcvc(i,j,itr)*1000./g
       ddsfc(i,j,itr)=1.e6*ddsfc(i,j,itr)*1000./g
       wxsg(i,j,itr)= 1.e6*wxsg(i,j,itr)*1000./g
       wxaq(i,j,itr)= 1.e6*wxaq(i,j,itr)*1000./g
c emtrac isbuilt from chsurfem so just need the 1e6*dt/2 factor to
c to pass im mg/m2
       cemtrac(i,j,itr) = 1.e6 *cemtrac(i,j,itr)
40    continue

c-----print out the information:
#ifdef DIAG
#ifdef MPP1
      if(myid .eq. 0) then
#endif

        if (ifprt .and. mod(ntime,nprtfrq).eq.0) then    

c---- tracers

      print*,'************************************************'
      print*,' Budgets for tracers (intergrated quantitites)'
      print*,' day = ',ldatez,' *****'
      print*,'************************************************'
     
      do itr=1,ntr

c       errort= ttrace(itr,1)-
c     &     (tchie(itr) - (tchiad(itr)-tchitb(itr)+ 
c     &     tremdrd(itr,1)+ tremlsc(itr,1)+tremcvc(itr,1)))

       print *,'****************************' 
       print 11, itr,tchie(itr)
       print 12, itr, tchiad(itr)
       print 13, itr, tchitb(itr)
       print 14, itr, tremdrd(itr,1)
       print 15, itr, tremlsc(itr,1)
       print 16, itr, tremcvc(itr,1)
       print 17, itr, ttrace(itr,1)
c       print 18, itr, errort 

      end do

11     format('total emission tracer',I1,':',e12.5,' kg')
12     format('total advected tracer',I1,':',e12.5,' kg')
13     format('total diffused tracer',I1,':',e12.5,' kg')
14     format('total dry deposition tracer',I1,':',e12.5,' kg')      
15     format('total large scale wet dep.tracer',I1,':',e12.5,' kg')
16     format('total convective wet dep.tracer',I1,':',e12.5,' kg')  
17     format('total mass (at ktau), tracer',I1,':',e12.5,' kg')
18     format('error trac',I1,':',e12.5,' % ')
       
       end if
#ifdef MPP1
       endif
#endif
#endif

      return
      end 

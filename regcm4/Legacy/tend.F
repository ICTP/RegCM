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

      subroutine tend(iexec)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     this subroutine computes the tendencies of the prognostic       c
c     variables p*, u, v, and t.                                      c
c                                                                     c
c     p*u, p*v, p*t ,p*qv, and p*qc stored in main common block.      c
c                                                                     c
c     all the two-dimension arrays stored in main common block.       c
c                                                                     c
c     east/west boundary conditions stored in common block /bdycod/ . c
c                                                                     c
c     north/south boundary conditions stored in common block          c
c              /bdycod/.                                              c
c                                                                     c
c     all the integers stored in common block /param1/.               c
c                                                                     c
c     all the constants stored in common block /param1/.              c
c                                                                     c
c     iexec  : = 1 ; represents this subroutine is called for the     c
c                    first time in this forecast run.                 c
c              > 1 ; represents subsequent calls to this subroutine.  c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
      integer iexec
#     include "parameter.inc"
#     include "include/param1.cb"
#     include "include/param2.cb"
#     include "include/param3.cb"
#     include "include/main.cb"
#     include "include/mainchem.cb"
#     include "include/bdycod.cb"
#     include "include/cvaria.cb"
#     include "include/pmoist.cb"
#     include "include/rad.cb"
#     include "include/bats2.cb"
#     include "include/trachem.cb"
#     include "include/date.cb"
#     include "include/message.cb"

#ifdef MPP1
#     include "include/slice.cb"
      include 'mpif.h'
      integer ierr
      real(kind=8)  aermm
      common /aerm/ aermm(ix-1,kx,jxp)
#else
      real(kind=8)  aermm
      common /aerm/ aermm(ix-1,kx,jx-1)
#endif
      INTEGER         MDATEZ
      COMMON /DATENUM/MDATEZ(289276)
c
#ifdef MPP1
      integer status(MPI_STATUS_SIZE)
      real(kind=8)  ttld(ix,kx,jxp)
      real(kind=8)  ps4(ix,4,jxp),ps_4(ix,4,mjx)
      COMMON /TMP_DIAG/ps4,ps_4
      real(kind=8)  bdyewsnd(ix,kx*16+4),bdyewrcv(ix,kx*16+4)
      real(kind=8)  bdynssnd(nspgx,kx*16+4),bdynsrcv(nspgx,kx*16+4)
      COMMON /TMP_BDY/bdyewsnd,bdyewrcv,bdynssnd,bdynsrcv
      real(kind=8)  var1snd(ix,kx*(ntr+11)+1),var1rcv(ix,kx*(ntr+11)+1)
      COMMON /TMP_VAR3/var1snd,var1rcv
      real(kind=8)  var2snd(ix,kx*(ntr+5)*2),var2rcv(ix,kx*(ntr+5)*2)
      COMMON /TMP_VAR2/var2snd,var2rcv
      integer numrec,icons_mpi
#else
      real(kind=8)  ttld(ix,kx,jx)
#endif
      real(kind=8)  divl(ix,kx)
      integer i,j,k,n,itr,iptn,lev,icons
      real(kind=8)  t00pg,p00pg,alam,pgfaa1,ptntot,pt2tot,tv,tvavg
      real(kind=8)  xtm1,dudx,dvdx,dudy,dvdy,cell,rovcpm,tvc,tva,tvb
      real(kind=8)  psasum,sigpsa,tv1,tv2,tv3,tv4,rtbar
      real(kind=8)  qvbs,qvas,qcbs,qcas,ptnbar,pt2bar,xday
      real(kind=8)  psabar,xmsf
      real(kind=8)  dto2
Cchem2
      real(kind=8)  chibs, chias
Cchem2_
c
c----for pressure gradient force calculations
c
      t00pg=287.
      p00pg=101.325
      alam=6.5e-3
      pgfaa1=alam*r/g
c
c----------------------------------------------------------------------
c-----fill up the boundary slices:
c
c     if (iexec .eq. 1) then
      if (.not. ifrest .AND. iexec .eq. 1) then
       call bdyval(xtime,iexec)
       iexec = 2
      else
       iexec = 2
      end if
c
c----------------------------------------------------------------------
c*****for large domain, subroutine tend just needed to go through once.
c
c-----multiply ua and va by inverse of mapscale factor at dot point:
c
#ifdef MPP1
      do 30 j=1,jendl
#else
      do 30 j=1,jx
#endif
      do 30 k=1,kx
      do 30 i=1,ix
       ua(i,k,j)=ua(i,k,j)*msfd(i,j)
       va(i,k,j)=va(i,k,j)*msfd(i,j)
30    continue
#ifdef MPP1
      call MPI_SENDRECV(psa(1,jxp),ix,MPI_DOUBLE_PRECISION,ieast,1,
     &                  psa(1,0),ix,MPI_DOUBLE_PRECISION,iwest,1,
     &                  MPI_COMM_WORLD,status,ierr)
      call MPI_SENDRECV(psa(1,1),ix,MPI_DOUBLE_PRECISION,iwest,2,
     &                  psa(1,jxp+1),ix,MPI_DOUBLE_PRECISION,ieast,2,
     &                  MPI_COMM_WORLD,status,ierr)
#endif
c
c-----decouple u, v, t, qv, and qc
c
#ifdef MPP1
      do j=1,jendl
        if(myid .eq. 0 .and. j.eq.1) then
#else
      do j=1,jx
        if(j.eq.1) then
#endif
c-----------lateral slices:
c-----------west boundary:
          do k=1,kx
          do i=1,ix
            u(i,k,j)=uj1(i,k)
            v(i,k,j)=vj1(i,k)
          end do
          end do
          if(iboudy.eq.3.or.iboudy.eq.4) then
c..............inflow/outflow dependence:
            do k=1,kx
            do i=1,ix
              if(u(i,k,j).lt.0.) then
                v(i,k,j)=vj2(i,k)
                u(i,k,j)=uj2(i,k)
              end if
            end do
            end do
          end if
#ifdef MPP1
        else if(myid .eq. 0 .and. j.eq.2) then
#else
        else if(j.eq.2) then
#endif
          do k=1,kx
          do i=1,ix
            u(i,k,j)=uj2(i,k)
            v(i,k,j)=vj2(i,k)
          end do
          end do
          if(iboudy.eq.3.or.iboudy.eq.4) then
c..............inflow/outflow dependence:
            do k=1,kx
c.................south boundary:
              if(v(1,k,j).lt.0.) then
                v(1,k,j)=v(2,k,j)
                u(1,k,j)=u(2,k,j)
              end if
c.................north boundary:
              if(v(ix,k,j).ge.0.)then
                v(ix,k,j)=v(ix-1,k,j)
                u(ix,k,j)=u(ix-1,k,j)
              end if
            end do
          end if
#ifdef MPP1
        else if(myid .eq. nproc-1 .and. j.eq.jendl-1) then
#else
        else if(j.eq.jx-1) then
#endif
          do k=1,kx
          do i=1,ix
            u(i,k,j) = ujlx(i,k)
            v(i,k,j) = vjlx(i,k)
          end do
          end do
          if(iboudy.eq.3.or.iboudy.eq.4) then
c..............inflow/outflow dependence:
            do k=1,kx
c.................south boundary:
              if(v(1,k,j).lt.0.) then
                v(1,k,j)=v(2,k,j)
                u(1,k,j)=u(2,k,j)
              end if
c.................north boundary:
              if(v(ix,k,j).ge.0.)then
                v(ix,k,j)=v(ix-1,k,j)
                u(ix,k,j)=u(ix-1,k,j)
              end if
            end do
          end if
#ifdef MPP1
        else if(myid .eq. nproc-1 .and. j.eq.jendl) then
#else
        else if(j.eq.jx) then
#endif
c-----------east boundary:
c
c...........no inflow/outflow dependence:
c
          do k=1,kx
          do i=1,ix
            u(i,k,j)=ujl(i,k)
            v(i,k,j)=vjl(i,k)
          end do
          end do
          if(iboudy.eq.3.or.iboudy.eq.4) then
c..............inflow/outflow dependence:
            do k=1,kx
            do i=1,ix
              if(u(i,k,j).ge.0.) then
                v(i,k,j)=vjlx(i,k)
                u(i,k,j)=ujlx(i,k)
              end if
            end do
            end do
          end if
        else
c
c-----interior slice:
c-----interior points:
c
          do k=1,kx
          do i=3,ix-2
            psabar=0.25*(psa(i,j)+psa(i,j-1)+psa(i-1,j)+psa(i-1,j-1))
            xmsf=msfd(i,j)
            u(i,k,j)=ua(i,k,j)/(psabar*xmsf)
            v(i,k,j)=va(i,k,j)/(psabar*xmsf)
          end do
          end do
c
c-----------north/south boundary points:
c...........no inflow/outflow dependence:
c
          do k=1,kx
c..............for i=2 and i=ix-1:
            u(2,k,j) = ui2(k,j)
            u(ix-1,k,j) = uilx(k,j)
            v(2,k,j) = vi2(k,j)
            v(ix-1,k,j) = vilx(k,j)
c..............for i=1 and i=ix:
            u(1,k,j)=ui1(k,j)
            u(ix,k,j)=uil(k,j)
            v(1,k,j)=vi1(k,j)
            v(ix,k,j)=vil(k,j)
          end do
          if(iboudy.eq.3.or.iboudy.eq.4) then
c..............inflow/outflow dependence:
            do k=1,kx
c.................south boundary:
              if(v(1,k,j).lt.0.)then
                v(1,k,j)=v(2,k,j)
                u(1,k,j)=u(2,k,j)
              end if
c.................north boundary:
              if(v(ix,k,j).ge.0.)then
                v(ix,k,j)=v(ix-1,k,j)
                u(ix,k,j)=u(ix-1,k,j)
              end if
            end do
          end if
        end if
      end do
c
#ifdef MPP1
      do j=1,jendx
#else
      do j=1,jlx
#endif
      do k=1,kx
      do i=1,ilx
        t(i,k,j)=ta(i,k,j)/psa(i,j)
        qv(i,k,j)=qva(i,k,j)/psa(i,j)
        qc(i,k,j)=qca(i,k,j)/psa(i,j)
      end do
      end do
      end do
Cchem2
      if(ichem.eq.1) then
c
c-----call special tracer decoupling routine for multiple (ntr) species
c
       do n=1,ntr
#ifdef MPP1
       do j=1,jendx
#else
       do j=1,jlx
#endif
       do k=1,kx
       do i=1,ilx
         chi(i,k,j,n)=chia(i,k,j,n)/psa(i,j)
       end do
       end do
       end do
       end do
      end if
Cchem2_
c
C=======================================================================
#ifdef MPP1
      if(myid.ne.nproc-1) then
        do i=1,ix
          var1snd(i,1) = psb(i,jxp)
        end do
        do k=1,kx
        do i=1,ix
          var1snd(i,1      +k) = tb(i,k,jxp)
          var1snd(i,1+kx   +k) =qvb(i,k,jxp)
          var1snd(i,1+kx*2 +k) = ub(i,k,jxp)
          var1snd(i,1+kx*3 +k) = vb(i,k,jxp)
          var1snd(i,1+kx*4 +k) =  u(i,k,jxp)
          var1snd(i,1+kx*5 +k) =  v(i,k,jxp)
          var1snd(i,1+kx*6 +k) =  t(i,k,jxp)
          var1snd(i,1+kx*7 +k) = qv(i,k,jxp)
          var1snd(i,1+kx*8 +k) = qc(i,k,jxp)
          var1snd(i,1+kx*9 +k) = ua(i,k,jxp)
          var1snd(i,1+kx*10+k) = va(i,k,jxp)
        end do
        end do
        if(ichem.eq.1) then
          do n=1,ntr
            do k=1,kx
            do i=1,ix
              var1snd(i,kx*11+1+(n-1)*kx+k) = chi(i,k,jxp,n)
            end do
            end do
          end do
        end if
      end if
      numrec = kx*11+1
      if(ichem.eq.1) numrec = kx*(ntr+11)+1
      call MPI_SENDRECV(var1snd(1,1),ix*numrec,MPI_DOUBLE_PRECISION,
     &          ieast,1,var1rcv(1,1),ix*numrec,MPI_DOUBLE_PRECISION,
     &          iwest,1,MPI_COMM_WORLD,status,ierr)
      if(myid.ne.0) then
        do i=1,ix
          psb(i,0) = var1rcv(i,1)
        end do
        do k=1,kx
        do i=1,ix
           tb(i,k,0) = var1rcv(i,1      +k)
          qvb(i,k,0) = var1rcv(i,1+kx   +k)
           ub(i,k,0) = var1rcv(i,1+kx*2 +k)
           vb(i,k,0) = var1rcv(i,1+kx*3 +k)
            u(i,k,0) = var1rcv(i,1+kx*4 +k)
            v(i,k,0) = var1rcv(i,1+kx*5 +k)
            t(i,k,0) = var1rcv(i,1+kx*6 +k)
           qv(i,k,0) = var1rcv(i,1+kx*7 +k)
           qc(i,k,0) = var1rcv(i,1+kx*8 +k)
           ua(i,k,0) = var1rcv(i,1+kx*9 +k)
           va(i,k,0) = var1rcv(i,1+kx*10+k)
        end do
        end do
        if(ichem.eq.1) then
          do n=1,ntr
            do k=1,kx
            do i=1,ix
              chi(i,k,0,n) = var1rcv(i,kx*11+1+(n-1)*kx+k)
            end do
            end do
          end do
        end if
      end if
C
      if(myid.ne.0) then
        do i=1,ix
          var1snd(i,1) = psb(i,1)
        end do
        do k=1,kx
        do i=1,ix
          var1snd(i,1      +k) = tb(i,k,1)
          var1snd(i,1+kx   +k) =qvb(i,k,1)
          var1snd(i,1+kx*2 +k) = ub(i,k,1)
          var1snd(i,1+kx*3 +k) = vb(i,k,1)
          var1snd(i,1+kx*4 +k) =  u(i,k,1)
          var1snd(i,1+kx*5 +k) =  v(i,k,1)
          var1snd(i,1+kx*6 +k) =  t(i,k,1)
          var1snd(i,1+kx*7 +k) = qv(i,k,1)
          var1snd(i,1+kx*8 +k) = qc(i,k,1)
          var1snd(i,1+kx*9 +k) = ua(i,k,1)
          var1snd(i,1+kx*10+k) = va(i,k,1)
        end do
        end do
        if(ichem.eq.1) then
          do n=1,ntr
            do k=1,kx
            do i=1,ix
              var1snd(i,kx*11+1+(n-1)*kx+k) = chi(i,k,1,n)
            end do
            end do
          end do
        end if
      end if
      numrec = kx*11+1
      if(ichem.eq.1) numrec = kx*(ntr+11)+1
      call MPI_SENDRECV(var1snd(1,1),ix*numrec,MPI_DOUBLE_PRECISION,
     &          iwest,2,var1rcv(1,1),ix*numrec,MPI_DOUBLE_PRECISION,
     &          ieast,2,MPI_COMM_WORLD,status,ierr)
      if(myid.ne.nproc-1) then
        do i=1,ix
          psb(i,jxp+1) = var1rcv(i,1)
        end do
        do k=1,kx
        do i=1,ix
           tb(i,k,jxp+1) = var1rcv(i,1      +k)
          qvb(i,k,jxp+1) = var1rcv(i,1+kx   +k)
           ub(i,k,jxp+1) = var1rcv(i,1+kx*2 +k)
           vb(i,k,jxp+1) = var1rcv(i,1+kx*3 +k)
            u(i,k,jxp+1) = var1rcv(i,1+kx*4 +k)
            v(i,k,jxp+1) = var1rcv(i,1+kx*5 +k)
            t(i,k,jxp+1) = var1rcv(i,1+kx*6 +k)
           qv(i,k,jxp+1) = var1rcv(i,1+kx*7 +k)
           qc(i,k,jxp+1) = var1rcv(i,1+kx*8 +k)
           ua(i,k,jxp+1) = var1rcv(i,1+kx*9 +k)
           va(i,k,jxp+1) = var1rcv(i,1+kx*10+k)
        end do
        end do
        if(ichem.eq.1) then
          do n=1,ntr
            do k=1,kx
            do i=1,ix
              chi(i,k,jxp+1,n) = var1rcv(i,kx*11+1+(n-1)*kx+k)
            end do
            end do
          end do
        end if
      end if
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
          pdotb(i,j)=0.25*(psb(i,j)+psb(i-1,j)+psb(i,j-1)+psb(i-1,j-1))
        end do
      end do
c
c-----east and west boundaries:
c
      do i=2,ilx
#ifdef MPP1
      if (myid .eq. 0) then
#endif
        pdotb(i,1)=0.5*(psb(i,1)+psb(i-1,1))
#ifdef MPP1
      end if
      if (myid .eq. nproc-1) then
        pdotb(i,jendl)=0.5*(psb(i,jendx)+psb(i-1,jendx))
      end if
#else
        pdotb(i,jx)=0.5*(psb(i,jlx)+psb(i-1,jlx))
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
        pdotb(1,j)=0.5*(psb(1,j)+psb(1,j-1))
        pdotb(ix,j)=0.5*(psb(ilx,j)+psb(ilx,j-1))
      end do
c
c-----corner points:
c
#ifdef MPP1
      if (myid .eq. 0) then
#endif
      pdotb(1,1)=psb(1,1)
      pdotb(ix,1)=psb(ilx,1)
#ifdef MPP1
      end if
      if (myid .eq. nproc-1) then
        pdotb(1,jendl)=psb(1,jendx)
        pdotb(ix,jendl)=psb(ilx,jendx)
      end if
#else
      pdotb(1,jx)=psb(1,jlx)
      pdotb(ix,jx)=psb(ilx,jlx)
#endif
C=======================================================================
      call slice3d
#ifdef MPP1
      if(myid.ne.nproc-1) then
        do k=1,kx
        do i=1,ix
          var2snd(i,     +k) = ubd3d(i,k,jxp-1)
          var2snd(i,kx   +k) = ubd3d(i,k,jxp)
          var2snd(i,kx*2 +k) = vbd3d(i,k,jxp-1)
          var2snd(i,kx*3 +k) = vbd3d(i,k,jxp)
          var2snd(i,kx*4 +k) =  tb3d(i,k,jxp-1)
          var2snd(i,kx*5 +k) =  tb3d(i,k,jxp)
          var2snd(i,kx*6 +k) = qvb3d(i,k,jxp-1)
          var2snd(i,kx*7 +k) = qvb3d(i,k,jxp)
          var2snd(i,kx*8 +k) = qcb3d(i,k,jxp-1)
          var2snd(i,kx*9 +k) = qcb3d(i,k,jxp)
        end do
        end do
        if(ichem.eq.1) then
          do n=1,ntr
            do k=1,kx
            do i=1,ix
              var2snd(i,kx*10+(n-1)*2*kx   +k) = chib3d(i,k,jxp-1,n)
              var2snd(i,kx*10+(n-1)*2*kx+kx+k) = chib3d(i,k,jxp  ,n)
            end do
            end do
          end do
        end if
      end if
      numrec = kx*5*2
      if(ichem.eq.1) numrec = kx*(ntr+5)*2
      call MPI_SENDRECV(var2snd(1,1),ix*numrec,MPI_DOUBLE_PRECISION,
     &          ieast,1,var2rcv(1,1),ix*numrec,MPI_DOUBLE_PRECISION,
     &          iwest,1,MPI_COMM_WORLD,status,ierr)
      if(myid.ne.0) then
        do k=1,kx
        do i=1,ix
          ubd3d(i,k,-1)= var2rcv(i,    +k)
          ubd3d(i,k,0) = var2rcv(i,kx  +k)
          vbd3d(i,k,-1)= var2rcv(i,kx*2+k)
          vbd3d(i,k,0) = var2rcv(i,kx*3+k)
           tb3d(i,k,-1)= var2rcv(i,kx*4+k)
           tb3d(i,k,0) = var2rcv(i,kx*5+k)
          qvb3d(i,k,-1)= var2rcv(i,kx*6+k)
          qvb3d(i,k,0) = var2rcv(i,kx*7+k)
          qcb3d(i,k,-1)= var2rcv(i,kx*8+k)
          qcb3d(i,k,0) = var2rcv(i,kx*9+k)
        end do
        end do
        if(ichem.eq.1) then
          do n=1,ntr
            do k=1,kx
            do i=1,ix
              chib3d(i,k,-1,n) = var2rcv(i,kx*10+(n-1)*2*kx   +k)
              chib3d(i,k,0 ,n) = var2rcv(i,kx*10+(n-1)*2*kx+kx+k)
            end do
            end do
          end do
        end if
      end if
C
      if(myid.ne.0) then
        do k=1,kx
        do i=1,ix
          var2snd(i,    +k) = ubd3d(i,k,1)
          var2snd(i,kx  +k) = ubd3d(i,k,2)
          var2snd(i,kx*2+k) = vbd3d(i,k,1)
          var2snd(i,kx*3+k) = vbd3d(i,k,2)
          var2snd(i,kx*4+k) =  tb3d(i,k,1)
          var2snd(i,kx*5+k) =  tb3d(i,k,2)
          var2snd(i,kx*6+k) = qvb3d(i,k,1)
          var2snd(i,kx*7+k) = qvb3d(i,k,2)
          var2snd(i,kx*8+k) = qcb3d(i,k,1)
          var2snd(i,kx*9+k) = qcb3d(i,k,2)
        end do
        end do
        if(ichem.eq.1) then
          do n=1,ntr
            do k=1,kx
            do i=1,ix
              var2snd(i,kx*10+(n-1)*2*kx   +k) = chib3d(i,k,1,n)
              var2snd(i,kx*10+(n-1)*2*kx+kx+k) = chib3d(i,k,2,n)
            end do
            end do
          end do
        end if
      end if
      numrec = kx*5*2
      if(ichem.eq.1) numrec = kx*(ntr+5)*2
      call MPI_SENDRECV(var2snd(1,1),ix*numrec,MPI_DOUBLE_PRECISION,
     &          iwest,2,var2rcv(1,1),ix*numrec,MPI_DOUBLE_PRECISION,
     &          ieast,2,MPI_COMM_WORLD,status,ierr)
      if(myid.ne.nproc-1) then
        do k=1,kx
        do i=1,ix
          ubd3d(i,k,jxp+1)= var2rcv(i,    +k)
          ubd3d(i,k,jxp+2)= var2rcv(i,kx  +k)
          vbd3d(i,k,jxp+1)= var2rcv(i,kx*2+k)
          vbd3d(i,k,jxp+2)= var2rcv(i,kx*3+k)
           tb3d(i,k,jxp+1)= var2rcv(i,kx*4+k)
           tb3d(i,k,jxp+2)= var2rcv(i,kx*5+k)
          qvb3d(i,k,jxp+1)= var2rcv(i,kx*6+k)
          qvb3d(i,k,jxp+2)= var2rcv(i,kx*7+k)
          qcb3d(i,k,jxp+1)= var2rcv(i,kx*8+k)
          qcb3d(i,k,jxp+2)= var2rcv(i,kx*9+k)
        end do
        end do
        if(ichem.eq.1) then
          do n=1,ntr
            do k=1,kx
            do i=1,ix
              chib3d(i,k,jxp+1,n) = var2rcv(i,kx*10+(n-1)*2*kx   +k)
              chib3d(i,k,jxp+2,n) = var2rcv(i,kx*10+(n-1)*2*kx+kx+k)
            end do
            end do
          end do
        end if
      end if
#endif
c
c**********************************************************************
c***** "j" loop begins here:
c
#ifdef MPP1
      do 1070 j=1,jendx
#else
      do 1070 j=1,jlx
#endif
c
       icon(j)=0
#ifdef MPP1
       if( (myid.eq.0       .and. j.eq.1) .or.
     &     (myid.eq.nproc-1 .and. j.eq.jendx) ) then
#else
       if(j .eq. 1 .or. j.eq.jlx) then
#endif
        do k=1,kxp1
        do i=1,ilx
         qdot(i,k,j) = 0.
        end do
        end do
       else
c
c----------------------------------------------------------------------
c**p**compute the pressure tendency:
c
        do 250 i=2,ilxm
         pten(i,j) = 0.
250     continue
        do 260 k=1,kx
        do 260 i=2,ilxm
         divl(i,k)=(ua(i+1,k,j+1)+ua(i,k,j+1)-ua(i+1,k,j)-ua(i,k,j))
     1            +(va(i+1,k,j+1)+va(i+1,k,j)-va(i,k,j+1)-va(i,k,j))
         pten(i,j)=pten(i,j)
     &            -divl(i,k)*dsigma(k)/(dx2*msfx(i,j)*msfx(i,j))
260     continue
c
c..p..compute vertical sigma-velocity (qdot):
c
        do 270 k=1,kxp1
        do 270 i=1,ilx
         qdot(i,k,j) = 0.
270     continue
        do 280 k=2,kx
        do 280 i=2,ilxm
         qdot(i,k,j)=qdot(i,k-1,j)-(pten(i,j)+divl(i,k-1)/
     1            (dx2*msfx(i,j)*msfx(i,j)))*dsigma(k-1)/psa(i,j)
280     continue
       end if
 1070 continue
#ifdef MPP1
      call MPI_SENDRECV(qdot(1,1,jxp),ix*kxp1,MPI_DOUBLE_PRECISION,
     &          ieast,1,qdot(1,1,0)  ,ix*kxp1,MPI_DOUBLE_PRECISION,
     &          iwest,1,MPI_COMM_WORLD,status,ierr)
      call MPI_SENDRECV(qdot(1,1,1),    ix*kxp1,MPI_DOUBLE_PRECISION,
     &          iwest,2,qdot(1,1,jxp+1),ix*kxp1,MPI_DOUBLE_PRECISION,
     &          ieast,2,MPI_COMM_WORLD,status,ierr)
#endif
c
c..p..compute omega
c
#ifdef MPP1
      do 1071 j=1,jendx
       if( (myid.eq.0       .and. j.eq.1) .or.
     &     (myid.eq.nproc-1 .and. j.eq.jendx) ) then
#else
      do 1071 j=1,jlx
       if (j .eq. 1 .or. j.eq.jlx) then
#endif
        do k=1,kx
        do i=1,ilx
         omega(i,k,j) = 0.
        end do
        end do
       else
        do 310 k=1,kx
        do 310 i=2,ilxm
         omega(i,k,j)=0.5*psa(i,j)*(qdot(i,k+1,j)+qdot(i,k,j))+a(k)*
     1          (pten(i,j)+((u(i,k,j)+u(i+1,k,j)+u(i+1,k,j+1)+
     2          u(i,k,j+1))*(psa(i,j+1)-psa(i,j-1))+(v(i,k,j)+
     3          v(i+1,k,j)+v(i+1,k,j+1)+v(i,k,j+1))*(psa(i+1,j)-
     4          psa(i-1,j)))/(dx8*msfx(i,j)))
310     continue
       end if
 1071 continue
#ifdef MPP1
      if(nspgx.ge.jxp) then
      do i=1,ix
        bdyewsnd(i,1) = peb(i,1)
        bdyewsnd(i,2) = pwb(i,jxp)
        bdyewsnd(i,3) = pebt(i,1)
        bdyewsnd(i,4) = pwbt(i,jxp)
      end do
      do k=1,kx
      do i=1,ix
        bdyewsnd(i,4      +k) = teb(i,k,1)
        bdyewsnd(i,4+kx   +k) = twb(i,k,jxp)
        bdyewsnd(i,4+kx*2 +k) = tebt(i,k,1)
        bdyewsnd(i,4+kx*3 +k) = twbt(i,k,jxp)
        bdyewsnd(i,4+kx*4 +k) = qeb(i,k,1)
        bdyewsnd(i,4+kx*5 +k) = qwb(i,k,jxp)
        bdyewsnd(i,4+kx*6 +k) = qebt(i,k,1)
        bdyewsnd(i,4+kx*7 +k) = qwbt(i,k,jxp)
        bdyewsnd(i,4+kx*8 +k) = ueb(i,k,1)
        bdyewsnd(i,4+kx*9 +k) = uwb(i,k,jxp)
        bdyewsnd(i,4+kx*10+k) = uebt(i,k,1)
        bdyewsnd(i,4+kx*11+k) = uwbt(i,k,jxp)
        bdyewsnd(i,4+kx*12+k) = veb(i,k,1)
        bdyewsnd(i,4+kx*13+k) = vwb(i,k,jxp)
        bdyewsnd(i,4+kx*14+k) = vebt(i,k,1)
        bdyewsnd(i,4+kx*15+k) = vwbt(i,k,jxp)
      end do
      end do
      call MPI_SENDRECV(bdyewsnd(1,1),ix*(kx*16+4),MPI_DOUBLE_PRECISION,
     &          ieast,1,bdyewrcv(1,1),ix*(kx*16+4),MPI_DOUBLE_PRECISION,
     &          iwest,1,MPI_COMM_WORLD,status,ierr)
      do i=1,ix
        if(myid.eq.nproc-1) then
          peb(i,jendl) = bdyewrcv(i,1)
          pebt(i,jendl)= bdyewrcv(i,3)
        else
          peb(i,jxp+1) = bdyewrcv(i,1)
          pebt(i,jxp+1)= bdyewrcv(i,3)
        end if
        pwb(i,0) = bdyewrcv(i,2)
        pwbt(i,0)= bdyewrcv(i,4)
      end do
      do k=1,kx
      do i=1,ix
        if(myid.eq.nproc-1) then
          teb(i,k,jendl) = bdyewrcv(i,4      +k)
          tebt(i,k,jendl)= bdyewrcv(i,4+kx*2 +k)
          qeb(i,k,jendl) = bdyewrcv(i,4+kx*4 +k)
          qebt(i,k,jendl)= bdyewrcv(i,4+kx*6 +k)
        else
          teb(i,k,jxp+1) = bdyewrcv(i,4      +k)
          tebt(i,k,jxp+1)= bdyewrcv(i,4+kx*2 +k)
          qeb(i,k,jxp+1) = bdyewrcv(i,4+kx*4 +k)
          qebt(i,k,jxp+1)= bdyewrcv(i,4+kx*6 +k)
        end if
        ueb(i,k,jxp+1) = bdyewrcv(i,4+kx*8 +k)
        uebt(i,k,jxp+1)= bdyewrcv(i,4+kx*10+k)
        veb(i,k,jxp+1) = bdyewrcv(i,4+kx*12+k)
        vebt(i,k,jxp+1)= bdyewrcv(i,4+kx*14+k)
        twb(i,k,0) = bdyewrcv(i,4+kx   +k)
        twbt(i,k,0)= bdyewrcv(i,4+kx*3 +k)
        qwb(i,k,0) = bdyewrcv(i,4+kx*5 +k)
        qwbt(i,k,0)= bdyewrcv(i,4+kx*7 +k)
        uwb(i,k,0) = bdyewrcv(i,4+kx*9 +k)
        uwbt(i,k,0)= bdyewrcv(i,4+kx*11+k)
        vwb(i,k,0) = bdyewrcv(i,4+kx*13+k)
        vwbt(i,k,0)= bdyewrcv(i,4+kx*15+k)
      end do
      end do
      do i=1,ix
        if(myid.eq.nproc-1) then
          bdyewsnd(i,1) = peb(i,jendx)
          bdyewsnd(i,3) = pebt(i,jendx)
        else
          bdyewsnd(i,1) = peb(i,jxp)
          bdyewsnd(i,3) = pebt(i,jxp)
        end if
        bdyewsnd(i,2) = pwb(i,1)
        bdyewsnd(i,4) = pwbt(i,1)
      end do
      do k=1,kx
      do i=1,ix
        if(myid.eq.nproc-1) then
          bdyewsnd(i,4     +k) = teb(i,k,jendx)
          bdyewsnd(i,4+kx*2+k) = tebt(i,k,jendx)
          bdyewsnd(i,4+kx*4+k) = qeb(i,k,jendx)
          bdyewsnd(i,4+kx*6+k) = qebt(i,k,jendx)
        else
          bdyewsnd(i,4     +k) = teb(i,k,jxp)
          bdyewsnd(i,4+kx*2+k) = tebt(i,k,jxp)
          bdyewsnd(i,4+kx*4+k) = qeb(i,k,jxp)
          bdyewsnd(i,4+kx*6+k) = qebt(i,k,jxp)
        end if
        bdyewsnd(i,4+kx*8 +k) = ueb(i,k,jxp)
        bdyewsnd(i,4+kx*10+k) = uebt(i,k,jxp)
        bdyewsnd(i,4+kx*12+k) = veb(i,k,jxp)
        bdyewsnd(i,4+kx*14+k) = vebt(i,k,jxp)
        bdyewsnd(i,4+kx   +k) = twb(i,k,1)
        bdyewsnd(i,4+kx*3 +k) = twbt(i,k,1)
        bdyewsnd(i,4+kx*5 +k) = qwb(i,k,1)
        bdyewsnd(i,4+kx*7 +k) = qwbt(i,k,1)
        bdyewsnd(i,4+kx*9 +k) = uwb(i,k,1)
        bdyewsnd(i,4+kx*11+k) = uwbt(i,k,1)
        bdyewsnd(i,4+kx*13+k) = vwb(i,k,1)
        bdyewsnd(i,4+kx*15+k) = vwbt(i,k,1)
      end do 
      end do 
      call MPI_SENDRECV(bdyewsnd(1,1),ix*(kx*16+4),MPI_DOUBLE_PRECISION,
     &          iwest,2,bdyewrcv(1,1),ix*(kx*16+4),MPI_DOUBLE_PRECISION,
     &          ieast,2,MPI_COMM_WORLD,status,ierr)
      do i=1,ix
        peb(i,0) = bdyewrcv(i,1)
        pebt(i,0) = bdyewrcv(i,3)
        pwb(i,jxp+1) = bdyewrcv(i,2)
        pwbt(i,jxp+1) = bdyewrcv(i,4)
      end do
      do k=1,kx
      do i=1,ix
        teb(i,k,0)      = bdyewrcv(i,4      +k)
        twb(i,k,jxp+1)  = bdyewrcv(i,4+kx   +k)
        tebt(i,k,0)     = bdyewrcv(i,4+kx*2 +k)
        twbt(i,k,jxp+1) = bdyewrcv(i,4+kx*3 +k)
        qeb(i,k,0)      = bdyewrcv(i,4+kx*4 +k)
        qwb(i,k,jxp+1)  = bdyewrcv(i,4+kx*5 +k)
        qebt(i,k,0)     = bdyewrcv(i,4+kx*6 +k)
        qwbt(i,k,jxp+1) = bdyewrcv(i,4+kx*7 +k)
        ueb(i,k,0)      = bdyewrcv(i,4+kx*8 +k)
        uwb(i,k,jxp+1)  = bdyewrcv(i,4+kx*9 +k)
        uebt(i,k,0)     = bdyewrcv(i,4+kx*10+k)
        uwbt(i,k,jxp+1) = bdyewrcv(i,4+kx*11+k)
        veb(i,k,0)      = bdyewrcv(i,4+kx*12+k)
        vwb(i,k,jxp+1)  = bdyewrcv(i,4+kx*13+k)
        vebt(i,k,0)     = bdyewrcv(i,4+kx*14+k)
        vwbt(i,k,jxp+1) = bdyewrcv(i,4+kx*15+k)
      end do 
      end do 
      end if
C
      if(myid.ne.nproc-1) then
      do i=1,nspgx
        bdynssnd(i,1) = pnb(i,jxp)
        bdynssnd(i,2) = pnbt(i,jxp)
        bdynssnd(i,3) = pss(i,jxp)
        bdynssnd(i,4) = psbt(i,jxp)
      end do
      do k=1,kx
      do i=1,nspgx
        bdynssnd(i,4      +k) = tnb(i,k,jxp)
        bdynssnd(i,4+kx   +k) = tnbt(i,k,jxp)
        bdynssnd(i,4+kx*2 +k) = tsb(i,k,jxp)
        bdynssnd(i,4+kx*3 +k) = tsbt(i,k,jxp)
        bdynssnd(i,4+kx*4 +k) = qnb(i,k,jxp)
        bdynssnd(i,4+kx*5 +k) = qnbt(i,k,jxp)
        bdynssnd(i,4+kx*6 +k) = qsb(i,k,jxp)
        bdynssnd(i,4+kx*7 +k) = qsbt(i,k,jxp)
        bdynssnd(i,4+kx*8 +k) = unb(i,k,jxp)
        bdynssnd(i,4+kx*9 +k) = unbt(i,k,jxp)
        bdynssnd(i,4+kx*10+k) = usb(i,k,jxp)
        bdynssnd(i,4+kx*11+k) = usbt(i,k,jxp)
        bdynssnd(i,4+kx*12+k) = vnb(i,k,jxp)
        bdynssnd(i,4+kx*13+k) = vnbt(i,k,jxp)
        bdynssnd(i,4+kx*14+k) = vsb(i,k,jxp)
        bdynssnd(i,4+kx*15+k) = vsbt(i,k,jxp)
      end do
      end do
      end if
      call MPI_SENDRECV(
     &          bdynssnd(1,1),nspgx*(kx*16+4),MPI_DOUBLE_PRECISION,
     &  ieast,1,bdynsrcv(1,1),nspgx*(kx*16+4),MPI_DOUBLE_PRECISION,
     &  iwest,1,MPI_COMM_WORLD,status,ierr)
      if(myid.ne.0) then
      do i=1,nspgx
        pnb(i,0) = bdynsrcv(i,1)
        pnbt(i,0)= bdynsrcv(i,2)
        pss(i,0) = bdynsrcv(i,3)
        psbt(i,0)= bdynsrcv(i,4)
      end do
      do k=1,kx
      do i=1,nspgx
        tnb(i,k,0) = bdynsrcv(i,4      +k)
        tnbt(i,k,0)= bdynsrcv(i,4+kx   +k)
        tsb(i,k,0) = bdynsrcv(i,4+kx*2 +k)
        tsbt(i,k,0)= bdynsrcv(i,4+kx*3 +k)
        qnb(i,k,0) = bdynsrcv(i,4+kx*4 +k)
        qnbt(i,k,0)= bdynsrcv(i,4+kx*5 +k)
        qsb(i,k,0) = bdynsrcv(i,4+kx*6 +k)
        qsbt(i,k,0)= bdynsrcv(i,4+kx*7 +k)
        unb(i,k,0) = bdynsrcv(i,4+kx*8 +k)
        unbt(i,k,0)= bdynsrcv(i,4+kx*9 +k)
        usb(i,k,0) = bdynsrcv(i,4+kx*10+k)
        usbt(i,k,0)= bdynsrcv(i,4+kx*11+k)
        vnb(i,k,0) = bdynsrcv(i,4+kx*12+k)
        vnbt(i,k,0)= bdynsrcv(i,4+kx*13+k)
        vsb(i,k,0) = bdynsrcv(i,4+kx*14+k)
        vsbt(i,k,0)= bdynsrcv(i,4+kx*15+k)
      end do
      end do
      end if
C
      if(myid.ne.0) then
      do i=1,nspgx
        bdynssnd(i,1) = pnb(i,1)
        bdynssnd(i,2) = pnbt(i,1)
        bdynssnd(i,3) = pss(i,1)
        bdynssnd(i,4) = psbt(i,1)
      end do
      do k=1,kx
      do i=1,nspgx
        bdynssnd(i,4      +k) = tnb(i,k,1)
        bdynssnd(i,4+kx   +k) = tnbt(i,k,1)
        bdynssnd(i,4+kx*2 +k) = tsb(i,k,1)
        bdynssnd(i,4+kx*3 +k) = tsbt(i,k,1)
        bdynssnd(i,4+kx*4 +k) = qnb(i,k,1)
        bdynssnd(i,4+kx*5 +k) = qnbt(i,k,1)
        bdynssnd(i,4+kx*6 +k) = qsb(i,k,1)
        bdynssnd(i,4+kx*7 +k) = qsbt(i,k,1)
        bdynssnd(i,4+kx*8 +k) = unb(i,k,1)
        bdynssnd(i,4+kx*9 +k) = unbt(i,k,1)
        bdynssnd(i,4+kx*10+k) = usb(i,k,1)
        bdynssnd(i,4+kx*11+k) = usbt(i,k,1)
        bdynssnd(i,4+kx*12+k) = vnb(i,k,1)
        bdynssnd(i,4+kx*13+k) = vnbt(i,k,1)
        bdynssnd(i,4+kx*14+k) = vsb(i,k,1)
        bdynssnd(i,4+kx*15+k) = vsbt(i,k,1)
      end do
      end do
      end if
      call MPI_SENDRECV(
     &          bdynssnd(1,1),nspgx*(kx*16+4),MPI_DOUBLE_PRECISION,
     &  iwest,2,bdynsrcv(1,1),nspgx*(kx*16+4),MPI_DOUBLE_PRECISION,
     &  ieast,2,MPI_COMM_WORLD,status,ierr)
      if(myid.ne.nproc-1) then
      do i=1,nspgx
        pnb(i,jxp+1) = bdynsrcv(i,1)
        pnbt(i,jxp+1)= bdynsrcv(i,2)
        pss(i,jxp+1) = bdynsrcv(i,3)
        psbt(i,jxp+1)= bdynsrcv(i,4)
      end do
      do k=1,kx
      do i=1,nspgx
        tnb(i,k,jxp+1) = bdynsrcv(i,4      +k)
        tnbt(i,k,jxp+1)= bdynsrcv(i,4+kx   +k)
        tsb(i,k,jxp+1) = bdynsrcv(i,4+kx*2 +k)
        tsbt(i,k,jxp+1)= bdynsrcv(i,4+kx*3 +k)
        qnb(i,k,jxp+1) = bdynsrcv(i,4+kx*4 +k)
        qnbt(i,k,jxp+1)= bdynsrcv(i,4+kx*5 +k)
        qsb(i,k,jxp+1) = bdynsrcv(i,4+kx*6 +k)
        qsbt(i,k,jxp+1)= bdynsrcv(i,4+kx*7 +k)
        unb(i,k,jxp+1) = bdynsrcv(i,4+kx*8 +k)
        unbt(i,k,jxp+1)= bdynsrcv(i,4+kx*9 +k)
        usb(i,k,jxp+1) = bdynsrcv(i,4+kx*10+k)
        usbt(i,k,jxp+1)= bdynsrcv(i,4+kx*11+k)
        vnb(i,k,jxp+1) = bdynsrcv(i,4+kx*12+k)
        vnbt(i,k,jxp+1)= bdynsrcv(i,4+kx*13+k)
        vsb(i,k,jxp+1) = bdynsrcv(i,4+kx*14+k)
        vsbt(i,k,jxp+1)= bdynsrcv(i,4+kx*15+k)
      end do
      end do
      end if
      do 1072 j=jbegin,jendx
       if (.not.(myid .eq. nproc-1 .and.j .eq. jendx) ) then
#else
      do 1072 j=2,jlx
       if (j .ne. jlx) then
#endif
cEES omega change: I broke up the write onto two lines
c and commented out the if statment
c
        if (iboudy .eq. 4) then
c..p..apply sponge boundary conditions to pten:
         call sponge_P(ispgx,wgtx,pten(1,j),j)
c....apply the nudging boundary conditions:
        else if (iboudy .eq. 1. or. iboudy. eq. 5) then
         xtm1=xtime-dtmin
         if(dabs(xtime).lt.0.00001.and.ldatez.gt.IDATE0)
     &   xtm1=-dtmin
         call nudge_P(ispgx,fnudge,gnudge,xtm1,pten(1,j)
     &              ,c203,j,iboudy)
        end if
       end if      !end if(j.ne.jlx) test
 1072 continue
#ifdef MPP1
      do 1073 j=1,jendx
       if( myid.eq.0 .and. j.eq.1 ) then
#else
      do 1073 j=1,jlx
       if (j .eq. 1) then
#endif
        do i=1,ilx
         psc(i,j)=psb(i,j) + dt*pwbt(i,j)
         psd(i,j)=psa(i,j)
        end do
#ifdef MPP1
       else if (myid.eq.nproc-1 .and. j.eq.jendx ) then
#else
       else if (j.eq.jlx) then
#endif
        do i=1,ilx
         psd(i,j)=psa(i,j)
        end do
       else
c
c..p..forecast pressure:
c
        do 360 i=2,ilxm
         psc(i,j)=psb(i,j)+pten(i,j)*dt
360     continue
c
c..p..weighted p* (psd)
c
        do 370 i=2,ilxm
         psd(i,j)=psa(i,j)
370     continue
c
        psc(1,j)=psb(1,j)+dt*psbt(1,j)
        psc(ilx,j)=psb(ilx,j)+dt*pnbt(1,j)
        psd(1,j)=psa(1,j)
        psd(ilx,j)=psa(ilx,j)
       end if
 1073 continue
#ifdef MPP1
      call MPI_SENDRECV(psd(1,jxp),ix,MPI_DOUBLE_PRECISION,ieast,1,
     &                  psd(1,0),ix,MPI_DOUBLE_PRECISION,iwest,1,
     &                  MPI_COMM_WORLD,status,ierr)
#endif
c
c-----compute bleck (1977) noise parameters:
c
#ifdef MPP1
      do j=1,jendl
        do i=1,ix
          ps4(i,1,j) = pten(i,j)
          ps4(i,2,j) = psc(i,j)
          ps4(i,3,j) = psb(i,j)
          ps4(i,4,j) = psa(i,j)
        end do
      end do
      call MPI_GATHER(ps4(1,1,1),ix*4*jxp,MPI_DOUBLE_PRECISION,
     &               ps_4(1,1,1),ix*4*jxp,MPI_DOUBLE_PRECISION,
     &                                  0,MPI_COMM_WORLD,ierr)
#endif
      iptn = 0
      ptntot = 0.
      pt2tot = 0.
#ifdef MPP1
      if(myid.eq.0) then
      do j=2,mjx-2
        if (.not.(jyear.eq.jyear0.and.ktau.eq.0)) then
          do i=2,ilxm
            iptn=iptn+1
            ptntot=ptntot+dabs(ps_4(i,1,j))
            pt2tot=pt2tot+dabs((ps_4(i,2,j)+ps_4(i,3,j)-2.*ps_4(i,4,j))
     1                        /(0.25*dt*dt))
          end do
        end if
      end do
      end if
      call MPI_BCAST(iptn,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      call MPI_BCAST(ptntot,1,MPI_DOUBLE_PRECISION,0,
     &                        MPI_COMM_WORLD,ierr)
      call MPI_BCAST(pt2tot,1,MPI_DOUBLE_PRECISION,0,
     &                        MPI_COMM_WORLD,ierr)
#else
      do 1074 j=2,jlx
       if (j .ne. jlx) then
        if (.not.(jyear.eq.jyear0.and.ktau.eq.0)) then
         do 400 i=2,ilxm
          iptn=iptn+1
          ptntot=ptntot+dabs(pten(i,j))
          pt2tot=pt2tot+dabs((psc(i,j)+psb(i,j)-2.*psa(i,j))
     1           /(0.25*dt*dt))
400      continue
        end if
c
       end if      !end if(j.ne.jlx) test
 1074 continue
#endif
#ifdef MPP1
      do 1075 j=jbegin,jendx
#else
      do 1075 j=2,jlx
#endif
c
c------compute the horizontal diffusion coefficient and stored in xkc:
c      the values are calculated at cross points, but they also used
c      for dot-point variables.

       do 420 k=1,kx
       do 420 i=2,ilx
        dudx=ub(i,k,j+1)+ub(i+1,k,j+1)-ub(i,k,j)-ub(i+1,k,j)
        dvdx=vb(i,k,j+1)+vb(i+1,k,j+1)-vb(i,k,j)-vb(i+1,k,j)
        dudy=ub(i+1,k,j)+ub(i+1,k,j+1)-ub(i,k,j)-ub(i,k,j+1)
        dvdy=vb(i+1,k,j)+vb(i+1,k,j+1)-vb(i,k,j)-vb(i,k,j+1)
cfil    cell=(xkhz*hgfact(i,j)/5.+c200*dsqrt((dudx-dvdy)*(dudx-dvdy)
        cell=(xkhz*hgfact(i,j)+c200*dsqrt((dudx-dvdy)*(dudx-dvdy)
     1     +(dvdx+dudy)*(dvdx+dudy)))
        xkc(i,k,j)=dmin1(cell,xkhmax)
     
420    continue
 1075 continue
C
#ifdef MPP1
      do 1079 j=jbegin,jendx
       if (.not.(myid .eq. nproc-1 .and.j .eq. jendx) ) then
#else
      do 1079 j=2,jlx
       if (j .ne. jlx) then
#endif
c
c---------------------------------------------------------------------
c**t**compute the temperature tendency:
c
        do 430 k=1,kx
        do 430 i=2,ilxm
          tten(i,k,j) = 0.
          qvten(i,k,j)= 0.
          qcten(i,k,j)= 0.
430     continue

c
c..t..compute the horizontal advection term:
c
c       call hadv_T(tten(1,1,j),ua,va,t,msfx,dx4,j,1)
        call hadv_T(tten(1,1,j),dx4,j,1)
c
c..t..compute the vertical advection term:
c
        call vadv(tten(1,1,j),ta(1,1,j),j,1)
c
c..t..compute the adiabatic term:
c
        do 440 k=1,kx
        do 440 i=2,ilxm
         rovcpm=r/(cp*(1.+0.8*(qv(i,k,j))))
         tv=t(i,k,j)*(1.+ep1*(qv(i,k,j)))
         tten(i,k,j)=tten(i,k,j)+(omega(i,k,j)*rovcpm*tv)
     &            /(ptop/psa(i,j)+a(k))
440     continue
c
c..t..compute the diffusion term for t and store in difft:
c
        do 450 k=1,kx
        do 450 i = 1,ilx
         difft(i,k,j) = 0.
         diffq(i,k,j)=0.
450     continue
c
       call diffut_T(difft(1,1,j),xkc(1,1,j),c203,j)
c
c**q**compute the moisture tendencies:
c
c....icup   = 1 : kuo-anthes cumulus parameterizaion scheme
c....icup   = 2 : grell cumulus paramterization scheme
c....icup   = 3 : betts-miller (1986)
c....icup   = 4 : emanuel (1991)
c
       if (icup.ne.1)then
         call hadvQV(qvten(1,1,j),dx4,j,1)
         call vadv(qvten(1,1,j),qva(1,1,j),j,2)
       end if

       if (icup.eq.1)then
         call cupara(j)
       else if(icup.eq.2)then
         call cuparan(tten(1,1,j),qvten(1,1,j),j)
       else if(icup.eq.3) then
         write(aline,*) 'ICTP RegCM team thinks the Betts-Miller code',
     &              ' is not ready for Regional Climate Run yet.'
         call say
         call fatal(__FILE__,__LINE__,
     &        'BETTS MILLER CUMULUS OPTION NOT ALLOWED')
         call bmpara(tten(1,1,j),qvten(1,1,j),j)
       else if(icup.eq.4) then
         call cupemandrv(j)
       end if

       if(ipptls.eq.1) then
         call hadvQC(qcten(1,1,j),dx4,j,1)
c        call hadvQC(qcten(1,1,j),dx,j,2)
Cfix     call vadv(qcten(1,1,j),qca(1,1,j),j,3)
         call vadv(qcten(1,1,j),qca(1,1,j),j,5)
         call pcp(j)
         call cldfrac(j)

c  need also to set diffq to 0 here before calling diffut
         do 808 k=1,kx
         do 808 i=1,ilx
           diffq(i,k,j) = 0.
 808     continue

c-----compute the diffusion terms:
c     the diffusion term for qv is stored in diffq. before completing
c     qvten computation, do not use diffq for other purpose.

         call diffutQV(diffq(1,1,j),xkc(1,1,j),c203,j)
         call diffutQC(qcten(1,1,j),xkc(1,1,j),c203,j)
       end if
c
Cchem2 compute the tracers tendencies
       if(ichem .eq. 1) then 
          call zenitm(coszrs,ix,j)
          call tractend2(j)
       end if 
Cchem2_
c
      end if             !end if(j.ne.jlx) test
c----------------------------------------------------------------------
c*****compute the pbl fluxes:
c     the diffusion and pbl tendencies of t and qv are stored in
c     difft and diffq.
c
       do 490 k=1,kx
       do 490 i=2,ilx
        uten(i,k,j) = 0.
        vten(i,k,j) = 0.
490    continue


c ****** calculate solar zenith angle
      if((jyear.eq.jyear0.and.ktau.eq.0) .or. mod(ktau+1,nbatst).eq.0
     &              .or. mod(ktau+1,ntrad).eq.0) then
        call zenitm(coszrs,ix,j)
        call slice(j)
      end if

c ****** calculate albedo
      if( (jyear.eq.jyear0.and.ktau.eq.0) .or.
     &            mod(ktau+1,ntrad).eq.0 ) then
        call albedov(j,iemiss)
      end if

c ****** call ccm3 radiative transfer package
      if( (jyear.eq.jyear0.and.ktau.eq.0) .or. 
     &            mod(ktau+1,ntrad).eq.0 ) then
        call colmod3(j)
      end if

c ****** call vector bats for surface physics calculations
      if( (jyear.eq.jyear0.and.ktau.eq.0) .or. 
     &           mod(ktau+1,nbatst).eq.0 ) then
        dtbat = dt/2.*nbatst
        if (jyear.eq.jyear0.and.ktau.eq.0) dtbat = dt
        call vecbats(j)
        if (iocnflx.eq.2) then        ! Zeng ocean flux model
          call zengocndrv(j)
        end if
c ****** accumulate quantities for energy and moisture budgets
        call interf(2,j)
      end if

 1079 continue
      if(icup.eq.1) then
        dto2=dt/2
        call htdiff(dto2,dxsq,akht1)
      end if
c     call medium resolution pbl
      if (ibltyp.eq.1) then
        call holtbl
      end if
#ifdef MPP1
      do 1081 j=jbegin,jendx
#else
      do 1081 j=2,jlx
#endif

c   add ccm radiative transfer package-calculated heating rates to
c   temperature tendency
       do 487 k=1,kx
       do 487 i=2,ilxm
        tten(i,k,j)=tten(i,k,j)+psb(i,j)*heatrt(i,k,j) !heating rate in deg/sec
487    continue
c
#ifdef MPP1
       if(.not.(myid.eq.nproc-1.and.j.eq.jendx) ) then
#else
       if (j .ne. jlx) then
#endif
c
c..tq.add horizontal diffusion and pbl tendencies for t and qv to tten
c     and qvten for calculating condensational term in subroutine
c     "condtq".
c
        do 500 k=1,kx
        do 500 i=2,ilxm
         tten(i,k,j)=tten(i,k,j)+difft(i,k,j)
500     continue
c
        do 510 k=1,kx
        do 510 i=2,ilxm
         qvten(i,k,j)=qvten(i,k,j)+diffq(i,k,j)
510     continue
c
c..tq.compute the condensation and precipitation terms for explicit
c     moisture scheme:
c
        call condtq(j)
c
c..tq.subtract horizontal diffusion and pbl tendencies from tten and
c     qvten for appling the sponge boundary conditions on t and qv:
c
        if (iboudy .eq. 4) then
         do 530 k=1,kx
         do 530 i=2,ilxm
          tten(i,k,j)=tten(i,k,j)-difft(i,k,j)
530      continue
         call sponge_T(ispgx,wgtx,tten(1,1,j),j)
         do 540 k=1,kx
         do 540 i=2,ilxm
          tten(i,k,j)=tten(i,k,j)+difft(i,k,j)
540      continue
         do 550 k=1,kx
         do 550 i=2,ilxm
          qvten(i,k,j)=qvten(i,k,j)-diffq(i,k,j)
550      continue
         call spongeQV(ispgx,wgtx,qvten(1,1,j),j)
         do 560 k=1,kx
         do 560 i=2,ilxm
          qvten(i,k,j)=qvten(i,k,j)+diffq(i,k,j)
560      continue
        end if
c
c..tq.apply the nudging boundary conditions:
c
        if (iboudy .eq. 1. or. iboudy. eq. 5) then
         xtm1=xtime-dtmin
         if(dabs(xtime).lt.0.00001.and.ldatez.gt.IDATE0)
     &   xtm1=-dtmin
         call nudge_T(ispgx,fnudge,gnudge,xtm1,tten(1,1,j)
     &             ,c203,j,iboudy)
         call nudgeQV(ispgx,fnudge,gnudge,xtm1,qvten(1,1,j)
     &             ,c203,j,iboudy)
        end if
c
c..tq.forecast t, qv, and qc at tau+1:
c
        do 580 k=1,kx
        do 580 i=2,ilxm
         qvc(i,k,j)=qvb(i,k,j)+dt*qvten(i,k,j)
580     continue
c
        do 610 k=1,kx
        do 610 i=2,ilxm
         qcc(i,k,j)=qcb(i,k,j)+dt*qcten(i,k,j)
610     continue
c
        do 630 k=1,kx
        do 630 i=2,ilxm
         tc(i,k,j)=tb(i,k,j)+dt*tten(i,k,j)
630     continue
c
Cchem2 
c  forecast tracer chi at at tau+1:
        if(ichem.eq.1) then
c   
        do 600 itr=1,ntr
        do 600 k=1,kx
        do 600 i=2,ilxm
         chic(i,k,j,itr)=chib(i,k,j,itr)+dt*chiten(i,k,j,itr)
600     continue
        end if
Cchem2_         
       end if        !end if(j.ne.jlx),else test
 1081 continue
#ifdef MPP1
      do 1082 j=1,jendx
       if(myid.eq.0 .and. j.eq.1) then
#else
      do 1082 j=1,jlx
       if(j.eq.1) then
#endif
        if(ipgf.eq.1)then
         do 70 k=1,kx
         do 70 i=1,ilx
          td(i,k,j) = ta(i,k,j)*(1.+ep1*(qv(i,k,j)))
          ttld(i,k,j)=td(i,k,j)-psa(i,j)*t00pg
     &                          *((a(k)*psa(i,j)+ptop)/p00pg)**pgfaa1
70       continue
        else if(ipgf.eq.0)then
         do 71 k=1,kx
         do 71 i=1,ilx
          td(i,k,j) = ta(i,k,j)*(1.+ep1*(qv(i,k,j)))
71       continue
        end if
c
#ifdef MPP1
       else if (myid.eq.nproc-1 .and. j.eq. jendx) then
#else
       else if (j .eq. jlx) then
#endif
c
c-----set td and psd at j=jlx equal to ta and psa:
c
        if(ipgf.eq.1)then
         do 720 k=1,kx
         do 720 i=1,ilx
          td(i,k,j)=ta(i,k,j)*(1.+ep1*(qv(i,k,j)))
          ttld(i,k,j)=td(i,k,j)-psa(i,j)*t00pg*((a(k)*psa(i,j)+
     1       ptop)/p00pg)**pgfaa1
720      continue
        else if(ipgf.eq.0)then
         do 721 k=1,kx
         do 721 i=1,ilx
          td(i,k,j)=ta(i,k,j)*(1.+ep1*(qv(i,k,j)))
721      continue
        end if
c
       else
c
c..t..compute weighted p*t (td) for use in ssi:
c
        if(ipgf.eq.1)then
c
         do 660 k=1,kx
         do 660 i=2,ilxm
          tvc=tc(i,k,j)*(1.+ep1*(qvc(i,k,j))/psc(i,j))
          tva=ta(i,k,j)*(1.+ep1*(qv(i,k,j)))
          tvb=tb(i,k,j)*(1.+ep1*(qvb(i,k,j))/psb(i,j))
          td(i,k,j)=alpha*(tvc+tvb)+beta*tva
          ttld(i,k,j)=td(i,k,j)-psd(i,j)*t00pg*((a(k)*psd(i,j)+
     1       ptop)/p00pg)**pgfaa1
660      continue
         do 670 k=1,kx
          td(1,k,j)=ta(1,k,j)*(1.+ep1*(qv(1,k,j)))
          ttld(1,k,j)=td(1,k,j)-psa(1,j)*t00pg*((a(k)*psa(1,j)+
     1       ptop)/p00pg)**pgfaa1
          td(ilx,k,j)=ta(ilx,k,j)*(1.+ep1*(qv(ilx,k,j)))
          ttld(ilx,k,j)=td(ilx,k,j)-psa(ilx,j)*t00pg*((a(k)*
     1      psa(ilx,j)+ptop)/p00pg)**pgfaa1
670      continue
c
        else if(ipgf.eq.0)then
c
         do 661 k=1,kx
         do 661 i=2,ilxm
          tvc=tc(i,k,j)*(1.+ep1*(qvc(i,k,j))/psc(i,j))
          tva=ta(i,k,j)*(1.+ep1*(qv(i,k,j)))
          tvb=tb(i,k,j)*(1.+ep1*(qvb(i,k,j))/psb(i,j))
          td(i,k,j)=alpha*(tvc+tvb)+beta*tva
661      continue
         do 671 k=1,kx
          td(1,k,j)=ta(1,k,j)*(1.+ep1*(qv(1,k,j)))
          td(ilx,k,j)=ta(ilx,k,j)*(1.+ep1*(qv(ilx,k,j)))
671      continue
c
        end if  ! ipgf if block

       end if        !end if(j.ne.jlx),else test
 1082 continue
c----------------------------------------------------------------------
c**uv*compute the u and v tendencies:
#ifdef MPP1
      do 1083 j=jbegin,jendx
#else
      do 1083 j=2,jlx
#endif
c
c..uv.compute the diffusion terms:
c     put diffusion and pbl tendencies of u and v in difuu and difuv.
c
       do 830 k=1,kx
       do 830 i=2,ilx
        difuu(i,k,j) = uten(i,k,j)
        difuv(i,k,j) = vten(i,k,j)
830    continue
c
       call diffu_U(difuu(1,1,j),xkc(1,1,j),c203,j,1)
       call diffu_V(difuv(1,1,j),xkc(1,1,j),c203,j,1)
c
c..uv.compute the horizontal advection terms for u and v:
c
       do 840 k=1,kx
       do 840 i=2,ilx
        uten(i,k,j) = 0.
        vten(i,k,j) = 0.
840    continue
c
       call hadv_U(uten(1,1,j),dx16,j,3)
       call hadv_V(vten(1,1,j),dx16,j,3)
c
c..uv.compute coriolis terms:
c
       do 850 k=1,kx
       do 850 i=2,ilx
        uten(i,k,j)=uten(i,k,j)+f(i,j)*va(i,k,j)/msfd(i,j)
        vten(i,k,j)=vten(i,k,j)-f(i,j)*ua(i,k,j)/msfd(i,j)
850    continue
 1083  continue
#ifdef MPP1
      do 1086 j=jbegin,jendx
#else
      do 1086 j=2,jlx
#endif
c
c..uv.compute pressure gradient terms:
c
       if(ipgf.eq.1)then
        do 860 k=1,kx
        do 860 i=2,ilx
         psasum=psd(i,j)+psd(i-1,j)+psd(i,j-1)+psd(i-1,j-1)
         sigpsa=psasum
         tv1=t(i-1,k,j-1)*(1.+ep1*(qv(i-1,k,j-1)))
         tv2=t(i,k,j-1)*(1.+ep1*(qv(i,k,j-1)))
         tv3=t(i-1,k,j)*(1.+ep1*(qv(i-1,k,j)))
         tv4=t(i,k,j)*(1.+ep1*(qv(i,k,j)))
         rtbar=tv1+tv2+tv3+tv4-4.*t00pg*((a(k)*psasum/4.+ptop)
     1    /p00pg)**pgfaa1
         rtbar=r*rtbar*sigpsa/16.
         uten(i,k,j)=uten(i,k,j)-rtbar*(dlog(0.5*(psd(i,j)+psd(i-1,j))
     1    *a(k)+ptop)- dlog(0.5*(psd(i,j-1)+psd(i-1,j-1))*a(k)+ptop))
     1    /(dx*msfd(i,j))
         vten(i,k,j)=vten(i,k,j)-rtbar*(dlog(0.5*(psd(i,j)+psd(i,j-1))
     1    *a(k)+ptop)- dlog(0.5*(psd(i-1,j-1)+psd(i-1,j))*a(k)+ptop))
     1    /(dx*msfd(i,j))
860     continue
       else if(ipgf.eq.0)then
        do 861 k=1,kx
        do 861 i=2,ilx
         psasum=psd(i,j)+psd(i-1,j)+psd(i,j-1)+psd(i-1,j-1)
         sigpsa=psasum
         tv1=t(i-1,k,j-1)*(1.+ep1*(qv(i-1,k,j-1)))
         tv2=t(i,k,j-1)*(1.+ep1*(qv(i,k,j-1)))
         tv3=t(i-1,k,j)*(1.+ep1*(qv(i-1,k,j)))
         tv4=t(i,k,j)*(1.+ep1*(qv(i,k,j)))
         rtbar=r*(tv1+tv2+tv3+tv4)*sigpsa/16.
         uten(i,k,j)=uten(i,k,j)-rtbar*(dlog(0.5*(psd(i,j)+psd(i-1,j))
     1    *a(k)+ptop)- dlog(0.5*(psd(i,j-1)+psd(i-1,j-1))*a(k)+ptop))
     1    /(dx*msfd(i,j))
         vten(i,k,j)=vten(i,k,j)-rtbar*(dlog(0.5*(psd(i,j)+psd(i,j-1))
     1    *a(k)+ptop)- dlog(0.5*(psd(i-1,j-1)+psd(i-1,j))*a(k)+ptop))
     1    /(dx*msfd(i,j))
861     continue
       end if
1086  continue
#ifdef MPP1
      do 1087 j=1,jendx
#else
      do 1087 j=1,jlx
#endif
c
c..uv.compute geopotential height at half-k levels, cross points:
c
       if(ipgf.eq.1)then

        do 870 i=1,ilx
         tv = (ttld(i,kx,j)/psd(i,j))
     a      / ( 1. + qc(i,kx,j)/(1.+qv(i,kx,j)) )
         phi(i,kx,j) = ht(i,j) + r*t00pg/pgfaa1
     a                 *((psd(i,j)+ptop)/p00pg)**pgfaa1
         phi(i,kx,j) = phi(i,kx,j)
     a         - r*tv*dlog((a(kx)+ptop/psd(i,j))/(1.+ptop/psd(i,j)))
  870   continue

        do 890 k=1,kxm
         lev=kx-k
        do 890 i=1,ilx
         tvavg = ( (ttld(i,lev  ,j)*dsigma(lev  )
     a             +ttld(i,lev+1,j)*dsigma(lev+1) )
     b         / (psd(i,j)*(dsigma(lev)+dsigma(lev+1))))
     c    / ( 1. + qc(i,lev,j)/(1.+qv(i,lev,j)) )
         phi(i,lev,j) = phi(i,lev+1,j) - r*tvavg
     a        *dlog((a(lev)+ptop/psd(i,j))/(a(lev+1)+ptop/psd(i,j)))
  890   continue

       else if(ipgf.eq.0)then

        do 871 i=1,ilx
         tv = (td(i,kx,j)/psd(i,j))
     a      / ( 1. + qc(i,kx,j)/(1.+qv(i,kx,j)) )
         phi(i,kx,j) = ht(i,j)
     a         - r*tv*dlog((a(kx)+ptop/psd(i,j))/(1.+ptop/psd(i,j)))
  871   continue

        do 891 k=1,kxm
         lev=kx-k
        do 891 i=1,ilx
         tvavg = ( (td(i,lev  ,j)*dsigma(lev  )
     a             +td(i,lev+1,j)*dsigma(lev+1) )
     b         / (psd(i,j)*(dsigma(lev)+dsigma(lev+1))))
     c    / ( 1. + qc(i,lev,j)/(1.+qv(i,lev,j)) )
         phi(i,lev,j) = phi(i,lev+1,j) - r*tvavg
     a        *dlog((a(lev)+ptop/psd(i,j))/(a(lev+1)+ptop/psd(i,j)))
  891   continue

       end if  ! ipgf if block
1087  continue
#ifdef MPP1
       call MPI_SENDRECV(phi(1,1,jxp),ix*kx,MPI_DOUBLE_PRECISION,
     &           ieast,1,phi(1,1,0)  ,ix*kx,MPI_DOUBLE_PRECISION,
     &           iwest,1,MPI_COMM_WORLD,status,ierr)
#endif
#ifdef MPP1
      do 1088 j=jbegin,jendx
#else
      do 1088 j=2,jlx
#endif
c
c..uv.compute the geopotential gradient terms:
c
       do 950 k=1,kx
       do 950 i=2,ilx
        uten(i,k,j)=uten(i,k,j)
     &             -(psd(i-1,j-1)+psd(i,j-1)+psd(i-1,j)+psd(i,j)
     1         )*(phi(i,k,j)+phi(i-1,k,j)-phi(i,k,j-1)-phi(i-1,k,j-1))
     2          /(dx8*msfd(i,j))
        vten(i,k,j)=vten(i,k,j)
     &             -(psd(i-1,j-1)+psd(i,j-1)+psd(i-1,j)+psd(i,j)
     1         )*(phi(i,k,j)+phi(i,k,j-1)-phi(i-1,k,j)-phi(i-1,k,j-1))
     2          /(dx8*msfd(i,j))
950    continue
1088  continue
#ifdef MPP1
      do 1089 j=jbegin,jendx
#else
      do 1089 j=2,jlx
#endif
c
c..uv.compute teh vertical advection terms:
c
       call vadv(uten(1,1,j),ua(1,1,j),j,4)
       call vadv(vten(1,1,j),va(1,1,j),j,4)
c
c..uv.apply the sponge boundary condition on u and v:
c
       if (iboudy .eq. 4) then
        call sponge_U(ispgd,wgtd,uten(1,1,j),j)
        call sponge_V(ispgd,wgtd,vten(1,1,j),j)
       end if
c
c..uv.apply the nudging boundary conditions:
c
       if (iboudy .eq. 1. or. iboudy .eq. 5) then
        call nudge_U(ispgd,fnudge,gnudge,xtm1,uten(1,1,j)
     &              ,c203,j,iboudy)
        call nudge_V(ispgd,fnudge,gnudge,xtm1,vten(1,1,j)
     &            ,c203,j,iboudy)
       end if
c
c..uv.add the diffusion and pbl tendencies to uten and vten:
c
       do 970 k=1,kx
       do 970 i=2,ilx
        uten(i,k,j)=uten(i,k,j)+difuu(i,k,j)
        vten(i,k,j)=vten(i,k,j)+difuv(i,k,j)
970    continue
c
c..uv.forecast p*u and p*v at tau+1:
c
       do 980 k=1,kx
       do 980 i=2,ilx
        uc(i,k,j)=ub(i,k,j)+dt*uten(i,k,j)
        vc(i,k,j)=vb(i,k,j)+dt*vten(i,k,j)
980    continue
c
c*****end of j loop.
c**********************************************************************
1089  continue
c
c---------------------------------------------------------------------
c-----store the xxa variables in xxb and xxc in xxa:
c     perform time smoothing operations.
c
#ifdef MPP1
      do j=jbegin,jendx
#else
      do j=2,jlx
#endif
        do 1020 k=1,kx
        do 1020 i=2,ilx
         ub(i,k,j)=omuhf*ua(i,k,j)/msfd(i,j)+gnuhf*(ub(i,k,j)+uc(i,k,j))
         vb(i,k,j)=omuhf*va(i,k,j)/msfd(i,j)+gnuhf*(vb(i,k,j)+vc(i,k,j))
         ua(i,k,j)=uc(i,k,j)
         va(i,k,j)=vc(i,k,j)
1020    continue
      end do
#ifdef MPP1
      do j=jbegin,jendm
#else
      do j=2,jlxm
#endif
        do 1060 k=1,kx
         do 1030 i=2,ilxm
          tb(i,k,j)=omuhf*ta(i,k,j) +gnuhf*(tb(i,k,j)+tc(i,k,j))
          ta(i,k,j)=tc(i,k,j)
1030     continue
         do 1040 i=2,ilxm
          qvbs=omuhf*qva(i,k,j)+gnuhf*(qvb(i,k,j)+qvc(i,k,j))
          qvas=qvc(i,k,j)
          qvb(i,k,j)=dmax1(qvbs,1.d-99)
          qva(i,k,j)=dmax1(qvas,1.d-99)
1040     continue
         do 1050 i=2,ilxm
          qcbs=omu*qca(i,k,j) +gnu*(qcb(i,k,j)+qcc(i,k,j))
          qcb(i,k,j)=dmax1(qcbs,0.d0)
1050     continue
         do 1055 i=2,ilxm
          qcas=qcc(i,k,j)
          qca(i,k,j)=dmax1(qcas,0.d0)
1055     continue
Cchem2
         if(ichem.eq.1) then
         do 1045 itr=1,ntr
         do 1045 i=2,ilxm
          chibs=omu*chia(i,k,j,itr)
     &         +gnu*(chib(i,k,j,itr)+chic(i,k,j,itr))
          chib(i,k,j,itr)=dmax1(chibs,0.d0)
          chias=chic(i,k,j,itr)
          chia(i,k,j,itr)=dmax1(chias,0.d0)
1045     continue
         end if         
Cchem2_
1060    continue
        do i=2,ilxm
         psb(i,j)=omuhf*psa(i,j)+gnuhf*(psb(i,j)+psc(i,j))
         psa(i,j)=psc(i,j)
        end do
       end do
       if(EHSO4) then
          do k=1,kx
#ifdef MPP1
            do j=1,jendx
#else
            do j=1,jlx
#endif
              do i=1,ilx
                 aermm(i,k,j) = so4(i,k,j)
              end do
            end do
          end do
       end if
c
c----------------------------------------------------------------------
c-----increment elapsed forecast time:
c
      ktau = ktau+1
      xtime = xtime+dtmin
      ntime = ntime+nint(dtmin*60.)
      if(dabs(xtime-IBDYFRQ*60.).lt.0.00001) then
         lhour = lhour+IBDYFRQ
         if(lhour.eq.24) then
            call FINDDATE(NNNNNN,ldatez)
            ldatez = MDATEZ(NNNNNN+1)
            lyear = ldatez/1000000
            lmonth=(ldatez-lyear*1000000)/10000
            lday  =(ldatez-lyear*1000000-lmonth*10000)/100
            lhour = mod(ldatez,100)
         else
            ldatez = ldatez + IBDYFRQ
         end if
         NNNNNN = NNNNNN + 1
         xtime  = 0.0
         if(mod(ldatez,1000000).eq.10100.and.xtime.lt.0.0001) then
           jyear=ldatez/1000000
           ktau = 0
           ntime= 0
         end if
      end if
      if (.not.(jyear.eq.jyear0.and.ktau.eq.0)) dt=dt2
c
c-----compute the amounts advected through the lateral boundaries:
c     *** note *** we must calculate the amounts advected through
c                  the lateral boundaries before updating the values
c                  at boundary slices.
c
#ifdef DIAG
      call conadv
c
cChem2 
      if (ichem .eq. 1) then
       call tracdiag(xkc)
      end if   
cChem2_
#endif

c-----fill up the boundary values for xxb and xxa variables:
c
      call bdyval(xtime,iexec)
c
c-----compute the nonconvective precipitation:
c
c???      call nconvp(psa,psb,ta,tb,qva,qvb,qca,qcb)
c
Cchem2_
      if(ichem.eq.1 .and. ichcumtra .eq. 1) then
c       do cumulus transport of tracers      
        call cumtran
      end if
 
Cchem2_

c-----trace the mass conservation of dry air and water substance:
c
#ifdef DIAG
      call conmas
#endif
c 
c
c---- budgets for tracers   
      if(ichem.eq.1) then
        call tracbud      
      end if  
c
c-----print out noise parameter:
c
      if (.not.(jyear.eq.jyear0.and.ktau.le.1)) then
       ptnbar=ptntot/dble(iptn)
       pt2bar=pt2tot/dble(iptn)
       icons=0
#ifdef MPP1
       icons_mpi=0
       do 1195 j=jbegin,jendm
#else
       do 1195 j=2,jlxm
#endif
        icons=icons+icon(j)
1195   continue
#ifdef MPP1
       icons_mpi=0
       call MPI_ALLREDUCE(icons,icons_mpi,1,MPI_INTEGER,
     &                    MPI_SUM,MPI_COMM_WORLD,ierr)
#endif
       xday = ((NNNNNN-NSTRT0)*IBDYFRQ*60.+xtime-dtmin) / 1440.
#ifdef MPP1
       if (myid .eq. 0) then
       if(mod(ktau,50).eq.0)print 1170,xday,ktau,ptnbar,pt2bar
     &                                ,icons_mpi
#else
       if(mod(ktau,50).eq.0)print 1170,xday,ktau,ptnbar,pt2bar,icons
#endif
1170   format(5x,'at day = ',f9.4,', ktau = ',i10,
     a       ' :  1st, 2nd time deriv of ps = ',
     a       2e12.5,',  no. of points w/convection = ',i7)
#ifdef MPP1
       end if
#endif
      end if
c
c----------------------------------------------------------------------
c
c-----recalculate solar declination angle if forecast time larger than
c     24 hours:
c
      if(dabs(xtime).lt.0.00001.and.ldatez.ne.IDATE1) then
       call solar1(xtime)
       dectim = anint(1440. + dectim)
#ifdef MPP1
       if (myid .eq. 0) write(*,*) ' dectim = ',dectim
#else
       write(*,*) ' dectim = ',dectim
#endif
      end if
c
      return
      end subroutine tend

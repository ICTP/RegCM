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

      subroutine spinit( ptop, sigma, kv1)
c
c** compute vertical modes.
c
      implicit none
      real(kind=8)  ptop
      integer kv1
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#     include "include/param1.cb"
#     include "include/param2.cb"
#     include "include/iunits.cb"
#     include "include/bdycod.cb"
#     include "include/main.cb"
#     include "include/split.cb"
#ifdef MPP1
#     include "include/io.cb"
      include 'mpiregcm.h'
      integer status(MPI_STATUS_SIZE)
      integer ierr
      real(kind=8)  uuu,vvv
      common /bxq_tmq/ uuu(ix,kx,jxp+1),vvv(ix,kx,jxp+1)
      real(kind=8)   sav0d(ix,nsplit*2,jxp),sav_0d(ix,nsplit*2,mjx)
      COMMON /TMPSAV0d/ sav0d,sav_0d
      real(kind=8)   sav6(kx,8,jxp),sav_6(kx,8,mjx)
      COMMON /TMPSAV6/ sav6,sav_6
      real(kind=8)  var1snd(kx,8),var1rcv(kx,8)
      COMMON /TMP_var1/ var1snd,var1rcv
#else
      real(kind=8)  uuu,vvv
      common /bxq_tmq/ uuu(ix,kx,jx),vvv(ix,kx,jx)
#endif
c
      real(kind=8)  sigma(kv1)
      logical lstand
c
c  lstand = .true. if standard atmosphere t to be used (ignore input
c            tbarh and ps in that case).  otherwise, ps and tbarh must
c            be defined on input.  note that in either case, pt must
c            also be defined on input (common block named cvert).
c
      integer i,j,k,l,n,ns,ijlx
      real(kind=8)  fac,pdlog,eps1,eps
c
c ******dtau = time steps(in sec)for modes in split explicit is
c ******specified in namelist as array dtsplit
c
c** zero new arrays
      do ns=1,nsplit
#ifdef MPP1
      do j=1,jxp
#else
      do j=1,jx
#endif
      do i=1,ix
        dstor(i,j,ns) = 0.
        hstor(i,j,ns) = 0.
      enddo
      enddo
      enddo
c
c** compute m.
      do 10 ns=1,nsplit
       m(ns) = nint( dt / dtau(ns) )
       if(.not.(jyear.eq.jyear0.and.ktau.eq.0)) 
     & m(ns) = nint(.5*dt/dtau(ns) )
10    continue
#ifdef MPP1
      if (myid .eq. 0) print *, 'dt, dtau = ', dt, dtau
#else
      print *, 'dt, dtau = ', dt, dtau
#endif
c
c** compute pt, ps and tbarh for use in vmodes.
      ps = 0.
      do k=1,kx
        tbarh(k) = 0.
      enddo
      pt = ptop
#ifdef MPP1
      ijlx = ilx * jendx
      do 20 j=1,jendx
#else
      ijlx = ilx * jlx
      do 20 j=1,jlx
#endif
      do 20 i=1,ilx
       ps = ps + psa(i,j) / ijlx
20    continue
      do 25 k=1,kx
#ifdef MPP1
      do 25 j=1,jendx
#else
      do 25 j=1,jlx
#endif
      do 25 i=1,ilx
       tbarh(k) = tbarh(k) + ta(i,k,j) / ( psa(i,j) * ijlx )
25    continue
c
c** compute vertical modes.
      lstand = .true.
      if(.not.(jyear.eq.jyear0.and.ktau.eq.0)) lstand = .true.
      call vmodes( lstand, sigma, kv1)
c
c** subract a4 from a for use in computing am.
      do l=1,kx
      do k=1,kx
        a(k,l) = a(k,l) - a4(k,l)
      enddo
      enddo
c
c** compute am and an.
      do 50 n=1,nsplit
       an(n) = 0.
       do 30 l=1,kx
        an(n) = an(n) + dsigma(l) * zmatx(l,n)
30     continue
       do 40 k=1,kx
        am(k,n) = 0.
        tau(n,k) = 0.
40     continue
       do 45 l=1,kx
       do 45 k=1,kx
        am(k,n) = am(k,n) + a(k,l) * zmatx(l,n)
        tau(n,k) = tau(n,k) + r * zmatxr(n,l) * hydros(l,k)
45     continue
c
       do 47 k=1,kxp1
        varpa1(n,k) = 0.
47     continue
       do 48 l=1,kx
       do 48 k=1,kxp1
        varpa1(n,k) = varpa1(n,k) + r * zmatxr(n,l) * hydroc(l,k)
48     continue
50    continue
c
c** multiply am, an and zmatx by factor.
      do 100 l=1,nsplit
       fac = 2. * dt / (2.*dble(m(l))+1.)
       if(.not.(jyear.eq.jyear0.and.ktau.eq.0))
     & fac = dt / (2.*dble(m(l))+1.)
#ifdef MPP1
       if (myid .eq. 0) print *, 'm, fac = ', m(l), fac
#else
       print *, 'm, fac = ', m(l), fac
#endif
       an(l) = an(l) * fac
      do 100 k=1,kx
       zmatx(k,l) = zmatx(k,l) * fac
       am(k,l)    = am(k,l)    * fac
100   continue
c
      IF ( ifrest ) THEN
#ifdef MPP1
       if (myid .eq. 0) then
         read(iutrs) dstor_io
         read(iutrs) hstor_io
         read(iutrs) uj1,uj2,ujlx,ujl
         read(iutrs) ui1_io,ui2_io,uilx_io,uil_io
         read(iutrs) vj1,vj2,vjlx,vjl
         read(iutrs) vi1_io,vi2_io,vilx_io,vil_io
         do j=1,mjx
           do n=1,nsplit
             do i=1,ix
               sav_0d(i,n       ,j) =  dstor_io(i,j,n)
               sav_0d(i,n+nsplit,j) =  hstor_io(i,j,n)
             enddo
           enddo
         enddo
          do j=1,mjx
            do k=1,kx
              sav_6(k,1,j) =  ui1_io(k,j)
              sav_6(k,2,j) =  ui2_io(k,j)
              sav_6(k,3,j) = uilx_io(k,j)
              sav_6(k,4,j) =  uil_io(k,j)
              sav_6(k,5,j) =  vi1_io(k,j)
              sav_6(k,6,j) =  vi2_io(k,j)
              sav_6(k,7,j) = vilx_io(k,j)
              sav_6(k,8,j) =  vil_io(k,j)
            enddo
          enddo
       endif
       call MPI_SCATTER(sav_0d(1,1,1),ix*nsplit*2*jxp,
     &          MPI_DOUBLE_PRECISION,sav0d(1,1,1),ix*nsplit*2*jxp,
     &          MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
       do j=1,jendl
         do n=1,nsplit
           do i=1,ix
             dstor(i,j,n) =  sav0d(i,n       ,j)
             hstor(i,j,n) =  sav0d(i,n+nsplit,j)
           enddo
         enddo
       enddo
       call MPI_SCATTER(sav_6(1,1,1),kx*8*jxp,MPI_DOUBLE_PRECISION,
     &                   sav6(1,1,1),kx*8*jxp,MPI_DOUBLE_PRECISION,
     &                                       0,MPI_COMM_WORLD,ierr)
       do j=1,jendl
         do k=1,kx
           ui1(k,j) = sav6(k,1,j)
           ui2(k,j) = sav6(k,2,j)
           uilx(k,j)= sav6(k,3,j)
           uil(k,j) = sav6(k,4,j)
           vi1(k,j) = sav6(k,5,j)
           vi2(k,j) = sav6(k,6,j)
           vilx(k,j)= sav6(k,7,j)
           vil(k,j) = sav6(k,8,j)
         enddo
       enddo
       call MPI_BCAST(uj1,ix*kx,MPI_DOUBLE_PRECISION,
     &                        0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(uj2,ix*kx,MPI_DOUBLE_PRECISION,
     &                        0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(vj1,ix*kx,MPI_DOUBLE_PRECISION,
     &                        0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(vj2,ix*kx,MPI_DOUBLE_PRECISION,
     &                        0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(ujlx,ix*kx,MPI_DOUBLE_PRECISION,
     &                         0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(ujl,ix*kx,MPI_DOUBLE_PRECISION,
     &                        0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(vjlx,ix*kx,MPI_DOUBLE_PRECISION,
     &                         0,MPI_COMM_WORLD,ierr)
       call MPI_BCAST(vjl,ix*kx,MPI_DOUBLE_PRECISION,
     &                        0,MPI_COMM_WORLD,ierr)
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
        enddo
      endif
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
        enddo
      endif
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
        enddo
      endif
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
        enddo
      endif
#else
       read(iutrs) dstor
       read(iutrs) hstor
       read(iutrs) uj1,uj2,ujlx,ujl
       read(iutrs) ui1,ui2,uilx,uil
       read(iutrs) vj1,vj2,vjlx,vjl
       read(iutrs) vi1,vi2,vilx,vil
#endif
      ELSE
c
C=======================================================================
c******* divergence manipulations (0)
c
c** compute divergence z from u and v
c    ( u must be pstar * u ; similarly for v )
c    ( note: map scale factors have been inverted in model (init) )
c
      do k=1,kx
#ifdef MPP1
        do j=1,jendl
#else
        do j=1,jx
#endif
          do i=1,ix
            uuu(i,k,j) = ub(i,k,j) * msfd(i,j)
            vvv(i,k,j) = vb(i,k,j) * msfd(i,j)
          enddo
        enddo
      enddo
#ifdef MPP1
      call MPI_SENDRECV(uuu(1,1,1),    ix*kx,MPI_DOUBLE_PRECISION,
     &          iwest,2,uuu(1,1,jxp+1),ix*kx,MPI_DOUBLE_PRECISION,
     &          ieast,2,MPI_COMM_WORLD,status,ierr)
      call MPI_SENDRECV(vvv(1,1,1),    ix*kx,MPI_DOUBLE_PRECISION,
     &          iwest,2,vvv(1,1,jxp+1),ix*kx,MPI_DOUBLE_PRECISION,
     &          ieast,2,MPI_COMM_WORLD,status,ierr)
#endif
c
      do l=1,nsplit
#ifdef MPP1
        do j=1,jendl
#else
        do j=1,jx
#endif
          do i=1,ix
            dstor(i,j,l) = 0.
          enddo
        enddo
      enddo
      do l=1,nsplit
        do k=1,kx
#ifdef MPP1
          do j=1,jendx
#else
          do j=1,jlx
#endif
            do i=1,ilx
              fac = dx2 * msfx(i,j) * msfx(i,j)
              dstor(i,j,l) = dstor(i,j,l) + zmatxr(l,k) *
     a    ( -uuu(i+1,k,j)+uuu(i+1,k,j+1)-uuu(i,k,j)+uuu(i,k,j+1)
     a      +vvv(i+1,k,j)+vvv(i+1,k,j+1)-vvv(i,k,j)-vvv(i,k,j+1) )
     a     / fac
            enddo
          enddo
        enddo
      enddo
c
C=======================================================================
c
c******* geopotential manipulations
      do l=1,nsplit
        pdlog=varpa1(l,kxp1)*dlog( sigmah(kxp1) * pd + pt )
        eps1 =varpa1(l,kxp1)*sigmah(kxp1)/( sigmah(kxp1) * pd + pt )
#ifdef MPP1
        do j=1,jendx
#else
        do j=1,jlx
#endif
          do i=1,ilx
            eps = eps1 * (psb(i,j)-pd)
            hstor(i,j,l) = pdlog + eps
          enddo
        enddo
        do k=1,kx
          pdlog=varpa1(l,k) * dlog( sigmah(k) * pd + pt )
          eps1 =varpa1(l,k) * sigmah(k) / ( sigmah(k) * pd + pt )
#ifdef MPP1
          do j=1,jendx
#else
          do j=1,jlx
#endif
            do i=1,ilx
              eps = eps1 * (psb(i,j)-pd)
              hstor(i,j,l)=hstor(i,j,l)
     &                    +pdlog+tau(l,k)*tb(i,k,j)/psb(i,j) + eps
            enddo
          enddo
        enddo
      enddo
      ENDIF
c
      return
      end

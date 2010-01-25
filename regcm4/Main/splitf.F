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

      subroutine splitf( gnuhf )
c
c** compute deld, delh
c** integrate in time and add correction terms appropriately
c
      implicit none
      real(kind=8)  gnuhf
#     include "parameter.inc"
#     include "include/param1.cb"
#     include "include/main.cb"
#     include "include/cvaria.cb"
#     include "include/pbldim.cb"
#     include "include/split.cb"
#ifdef MPP1
#     include "mpiregcm.inc"
      integer status(MPI_STATUS_SIZE)
      integer ierr
      real(kind=8)  uu,vv
      common /bxq_tmp/ uu(ix,jxp+1),vv(ix,jxp+1)
      real(kind=8)  uuu,vvv
      common /bxq_tmq/ uuu(ix,kx,jxp+1),vvv(ix,kx,jxp+1)
      real(kind=8)  deld(ix,jxp,nsplit,3), delh(ix,0:jxp,nsplit,3)
      real(kind=8)  work(ix,jxp,3), psdot(ix,jxp)
      real(kind=8)  wksend(ix*nsplit),wkrecv(ix*nsplit)
      real(kind=8)  ddsum(ix,jxp,nsplit), dhsum(ix,0:jxp,nsplit)
      common /bxq_aaa/ddsum,dhsum,deld,delh,psdot,work
      integer ii
#else
      real(kind=8)  uu,vv
      common /bxq_tmp/ uu(ix,jx),vv(ix,jx)
      real(kind=8)  uuu,vvv
      common /bxq_tmq/ uuu(ix,kx,jx),vvv(ix,kx,jx)
      real(kind=8)  deld(ix,jx,nsplit,3), delh(ix,jx,nsplit,3)
      real(kind=8)  work(ix,jx,3), psdot(ix,jx)
      real(kind=8)  ddsum(ix,jx,nsplit), dhsum(ix,jx,nsplit)
      common /bxq_aaa/ddsum,dhsum,deld,delh,psdot,work
#endif
      integer i,j,k,n,l
      real(kind=8)  fac,pdlog,eps1,eps,gnuan,gnuam,gnuzm,x,y
c
      do l=1,3
      do n=1,nsplit
#ifdef MPP1
      do j=1,jendl
#else
      do j=1,jx
#endif
      do i=1,ix
        deld(i,j,n,l) = 0.
        delh(i,j,n,l) = 0.
      enddo
      enddo
      enddo
      enddo
c
c** compute pressure on dot grid
C=======================================================================
c
c  this routine determines p(.) from p(x) by a 4-point interpolation.
c  on the x-grid, a p(x) point outside the grid domain is assumed to
c  satisfy p(0,j)=p(1,j); p(ix,j)=p(ix-1,j); and similarly for the i's.
c
#ifdef MPP1
      call MPI_SENDRECV(psa(1,jxp),ix,MPI_DOUBLE_PRECISION,ieast,1,
     &                  psa(1,0),ix,MPI_DOUBLE_PRECISION,iwest,1,
     &                  MPI_COMM_WORLD,status,ierr)
#endif
#ifdef MPP1
      do j=jbegin,jendx
#else
      do j=2,jlx
#endif
        do i=2,ilx
          psdot(i,j)=0.25*(psa(i,j)+psa(i-1,j)+psa(i,j-1)+psa(i-1,j-1))
        enddo
      enddo
c
      do i=2,ilx
#ifdef MPP1
        if (myid .eq. 0) then
#endif
        psdot(i,1)=0.5*(psa(i,1)+psa(i-1,1))
#ifdef MPP1
        end if
        if (myid .eq. nproc-1) then
          psdot(i,jendl)=0.5*(psa(i,jendx)+psa(i-1,jendx))
        end if
#else
        psdot(i,jx)=0.5*(psa(i,jlx)+psa(i-1,jlx))
#endif
      enddo
c
#ifdef MPP1
      do j=jbegin,jendx
#else
      do j=2,jlx
#endif
        psdot(1,j)=0.5*(psa(1,j)+psa(1,j-1))
        psdot(ix,j)=0.5*(psa(ilx,j)+psa(ilx,j-1))
      enddo
c
#ifdef MPP1
      if (myid .eq. 0) then
#endif
      psdot(1,1)=psa(1,1)
      psdot(ix,1)=psa(ilx,1)
#ifdef MPP1
      end if
      if (myid .eq. nproc-1) then
        psdot(1,jendl)=psa(1,jendx)
        psdot(ix,jendl)=psa(ilx,jendx)
      end if
#else
      psdot(1,jx)=psa(1,jlx)
      psdot(ix,jx)=psa(ilx,jlx)
#endif
c
C=======================================================================
c
c** get deld(0), delh(0) from storage
      do n=1,nsplit
#ifdef MPP1
      do j=1,jendl
#else
      do j=1,jx
#endif
      do i=1,ix
        deld(i,j,n,1) = dstor(i,j,n)
        delh(i,j,n,1) = hstor(i,j,n)
      enddo
      enddo
      enddo
c
C=======================================================================
c******* divergence manipulations (f)
      do k=1,kx
#ifdef MPP1
        do j=1,jendl
#else
        do j=1,jx
#endif
          do i=1,ix
            uuu(i,k,j) = ua(i,k,j) * msfd(i,j)
            vvv(i,k,j) = va(i,k,j) * msfd(i,j)
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
      do l=1,nsplit
#ifdef MPP1
        do j=1,jendl
#else
        do j=1,jx
#endif
          do i=1,ix
            deld(i,j,l,3) = 0.
          enddo
        enddo
        do k=1,kx
#ifdef MPP1
          do j=1,jendx
#else
          do j=1,jlx
#endif
            do i=1,ilx
              fac = dx2 * msfx(i,j) * msfx(i,j)
              deld(i,j,l,3) = deld(i,j,l,3) + zmatxr(l,k) *
     &      ( -uuu(i+1,k,j)+uuu(i+1,k,j+1)-uuu(i,k,j)+uuu(i,k,j+1)
     &        +vvv(i+1,k,j)+vvv(i+1,k,j+1)-vvv(i,k,j)-vvv(i,k,j+1) )
     &     / fac
            enddo
          enddo
        enddo
      enddo
c
C=======================================================================
      
      do n=1,nsplit
#ifdef MPP1
      do j=1,jendl
#else
      do j=1,jx
#endif
      do i=1,ix
         deld(i,j,n,3) = deld(i,j,n,3) - deld(i,j,n,1)
      enddo
      enddo
      enddo
c
C=======================================================================
c******* divergence manipulations (0)
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
      do l=1,nsplit
#ifdef MPP1
        do j=1,jendl
#else
        do j=1,jx
#endif
          do i=1,ix
            deld(i,j,l,2) = 0.
          enddo
        enddo
        do k=1,kx
#ifdef MPP1
          do j=1,jendx
#else
          do j=1,jlx
#endif
            do i=1,ilx
              fac = dx2 * msfx(i,j) * msfx(i,j)
              deld(i,j,l,2) = deld(i,j,l,2) + zmatxr(l,k) *
     &      ( -uuu(i+1,k,j)+uuu(i+1,k,j+1)-uuu(i,k,j)+uuu(i,k,j+1)
     &        +vvv(i+1,k,j)+vvv(i+1,k,j+1)-vvv(i,k,j)-vvv(i,k,j+1) )
     &     / fac
            enddo
          enddo
        enddo
      enddo
c
C=======================================================================
      do n=1,nsplit
#ifdef MPP1
      do j=1,jendl
#else
      do j=1,jx
#endif
      do i=1,ix
         deld(i,j,n,1) = deld(i,j,n,1) - deld(i,j,n,2)
      enddo
      enddo
      enddo
c
C=======================================================================
c******* geopotential manipulations (f)
      do l=1,nsplit
        pdlog = varpa1(l,kxp1) * dlog( sigmah(kxp1) * pd + pt )
        eps1  = varpa1(l,kxp1) * sigmah(kxp1)/(sigmah(kxp1)*pd + pt )
#ifdef MPP1
        do j=1,jendx
#else
        do j=1,jlx
#endif
          do i=1,ilx
            eps = eps1 * (psa(i,j)-pd)
            delh(i,j,l,3) = pdlog + eps
          enddo
        enddo
        do k=1,kx
        pdlog = varpa1(l,k) * dlog( sigmah(k) * pd + pt )
        eps1  = varpa1(l,k) * sigmah(k) / ( sigmah(k) * pd + pt )
#ifdef MPP1
        do j=1,jendx
#else
        do j=1,jlx
#endif
          do i=1,ilx
            eps = eps1 * (psa(i,j)-pd)
            delh(i,j,l,3)=delh(i,j,l,3)
     &                   +pdlog+tau(l,k)*ta(i,k,j)/psa(i,j) + eps
          enddo
        enddo
        enddo
      enddo
C=======================================================================

      do n=1,nsplit
#ifdef MPP1
      do j=1,jendl
#else
      do j=1,jx
#endif
      do i=1,ix
         delh(i,j,n,3) = delh(i,j,n,3) - delh(i,j,n,1)
      enddo
      enddo
      enddo
c
C=======================================================================
c******* geopotential manipulations (0)
      do l=1,nsplit
        pdlog = varpa1(l,kxp1) * dlog( sigmah(kxp1) * pd + pt )
        eps1  = varpa1(l,kxp1) * sigmah(kxp1)/(sigmah(kxp1)*pd + pt )
#ifdef MPP1
        do j=1,jendx
#else
        do j=1,jlx
#endif
          do i=1,ilx
            eps = eps1 * (psb(i,j)-pd)
            delh(i,j,l,2) = pdlog + eps
          enddo
        enddo
        do k=1,kx
        pdlog = varpa1(l,k) * dlog( sigmah(k) * pd + pt )
        eps1  = varpa1(l,k) * sigmah(k) / ( sigmah(k) * pd + pt )
#ifdef MPP1
        do j=1,jendx
#else
        do j=1,jlx
#endif
          do i=1,ilx
            eps = eps1 * (psb(i,j)-pd)
            delh(i,j,l,2)=delh(i,j,l,2)
     &                   +pdlog+tau(l,k)*tb(i,k,j)/psb(i,j) + eps
          enddo
        enddo
        enddo
      enddo
C=======================================================================
      do n=1,nsplit
#ifdef MPP1
      do j=1,jendl
#else
      do j=1,jx
#endif
      do i=1,ix
         delh(i,j,n,1) = delh(i,j,n,1) - delh(i,j,n,2)
      enddo
      enddo
      enddo
c
c** put deld(0), delh(0) into storage
      do n=1,nsplit
#ifdef MPP1
      do j=1,jendl
#else
      do j=1,jx
#endif
      do i=1,ix
        dstor(i,j,n) = deld(i,j,n,2)
        hstor(i,j,n) = delh(i,j,n,2)
      enddo
      enddo
      enddo
c
c******* split explicit time integration
      call spstep( hbar, dx2, dtau, m)
c
c******* add corrections to t and p;  u and v
C=======================================================================
      do l=1,nsplit
        gnuan = gnuhf *an(l)
#ifdef MPP1
        do j=jbegin,jendm
#else
        do j=2,jlxm
#endif
          do i=2,ilxm
            psa(i,j) = psa(i,j) - an(l) * ddsum(i,j,l)
            psb(i,j) = psb(i,j) - gnuan * ddsum(i,j,l)
          enddo
        enddo
      enddo
      do l=1,nsplit
        do k=1,kx
          gnuam = gnuhf * am(k,l)
#ifdef MPP1
          do j=jbegin,jendm
#else
          do j=2,jlxm
#endif
            do i=2,ilxm
              ta(i,k,j) = ta(i,k,j) + am(k,l) * ddsum(i,j,l)
              tb(i,k,j) = tb(i,k,j) + gnuam   * ddsum(i,j,l)
            enddo
          enddo
        enddo
      enddo
C=======================================================================
#ifdef MPP1
      ii=0
      do l=1,nsplit
      do i=1,ix
        ii=ii+1
        wksend(ii)=dhsum(i,jxp,l)
      enddo
      enddo
      call MPI_SENDRECV(wksend(1),ix*nsplit,MPI_DOUBLE_PRECISION,
     &          ieast,1,wkrecv(1),ix*nsplit,MPI_DOUBLE_PRECISION,
     &          iwest,1,MPI_COMM_WORLD,status,ierr)
      ii=0
      do l=1,nsplit
      do i=1,ix
        ii=ii+1
        dhsum(i,0,l)=wkrecv(ii)
      enddo
      enddo
#endif
      do l=1,nsplit
        do k=1,kx
          gnuzm = gnuhf * zmatx(k,l)
#ifdef MPP1
          do j=jbegin,jendx
#else
          do j=2,jlx
#endif
            do i=2,ilx
              fac = psdot(i,j) / ( dx2 * msfd(i,j) )
              x = fac * ( dhsum(i,j  ,l) + dhsum(i-1,j  ,l)
     a                  - dhsum(i,j-1,l) - dhsum(i-1,j-1,l) )
              y = fac * ( dhsum(i,j  ,l) - dhsum(i-1,j  ,l)
     a                  + dhsum(i,j-1,l) - dhsum(i-1,j-1,l) )
c
              ua(i,k,j) = ua(i,k,j) - zmatx(k,l) * x
              va(i,k,j) = va(i,k,j) - zmatx(k,l) * y
              ub(i,k,j) = ub(i,k,j) - gnuzm * x
              vb(i,k,j) = vb(i,k,j) - gnuzm * y
            enddo
          enddo
        enddo
      enddo
c
C=======================================================================
c
      return
      end

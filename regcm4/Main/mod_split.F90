!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      module mod_split

      use mod_constants
      use mod_dynparam
      use mod_runparams
      use mod_main

      private

      public :: allocate_mod_split , splitf
      public :: sigmah , sdsigma
      public :: thetaf , thetah
      public :: tbarh , tbarf
      public :: hweigh , tweigh , hbar
      public :: hydroc , hydros , hydror
      public :: alpha1 , alpha2
      public :: varpa1 , varpa2
      public :: zmatx , zmatxr
      public :: x1 , iw2 , w1 , w2 , w3
      public :: e1 , e2 , e3
      public :: s1 , s2 , s3
      public :: a0 , a1 , a2 , a3 , a4
      public :: d1 , d2 , g1 , g2
      public :: ps , pd , tau
      public :: cpfac
      public :: dstor , hstor
      public :: m
      public :: uuu , vvv
      public :: am , an

      integer , allocatable , dimension(:) :: m
      real(8) , allocatable , dimension(:,:) :: a0 , a1 , a2 , a3 , a4 ,&
               & d1 , d2 , e1 , e2 , e3 , g1 , g2 , g3 , s1 , s2 , w1 , &
               & w2 , w3 , x1
      integer , allocatable , dimension(:) :: iw2
      real(8) , allocatable , dimension(:) :: tbarf , thetaf
      real(8) , allocatable , dimension(:) :: thetah , tweigh
!
      real(8) :: alpha1 , alpha2 , pd , ps
      real(8) , allocatable , dimension(:) :: cpfac , sdsigma , hbar ,  &
               & hweigh , tbarh
      real(8) , allocatable , dimension(:,:) :: hydroc , varpa1
      real(8) , allocatable , dimension(:,:) :: hydror , hydros , tau , &
               & zmatx , zmatxr
      real(8) , allocatable , dimension(:) :: sigmah
      real(8) , allocatable , dimension(:,:) :: varpa2
!
      real(8) , allocatable , dimension(:,:) :: am
      real(8) , allocatable , dimension(:) :: an

      real(8) , allocatable, dimension(:,:,:) :: dstor , hstor

      real(8) ,allocatable, dimension(:,:,:) :: ddsum
      real(8) ,allocatable, dimension(:,:,:,:) :: deld
      real(8) ,allocatable, dimension(:,:,:,:) :: delh
      real(8) ,allocatable, dimension(:,:,:) :: dhsum
      real(8) ,allocatable, dimension(:,:) :: psdot
      real(8) ,allocatable, dimension(:,:,:) :: work
      real(8) ,allocatable, dimension(:,:) :: uu , vv
      real(8) ,allocatable, dimension(:,:,:) :: uuu , vvv

      contains 

      subroutine allocate_mod_split
        implicit none
#ifdef MPP1
        allocate(dstor(iy,0:jxp+1,nsplit))
        allocate(hstor(iy,0:jxp+1,nsplit))
#else
        allocate(dstor(iy,jx,nsplit))
        allocate(hstor(iy,jx,nsplit))
#endif 
        allocate(m(nsplit))
        allocate(a0(kz,kz))
        allocate(a1(kz,kz))
        allocate(a2(kz,kz))
        allocate(a3(kz,kz))
        allocate(a4(kz,kz))
        allocate(d1(kz,kz))
        allocate(d2(kz,kz))
        allocate(e1(kz,kz))
        allocate(e2(kz,kz))
        allocate(e3(kz,kz))
        allocate(g1(kz,kz))
        allocate(g2(kz,kz))
        allocate(g3(kz,kz))
        allocate(s1(kz,kz))
        allocate(s2(kz,kz))
        allocate(w1(kz,kz))
        allocate(w2(kz,kz))
        allocate(x1(kz,kz))
        allocate(iw2(kz))
        allocate(thetah(kz))
        allocate(tweigh(kz))
        allocate(tbarf(kzp1))
        allocate(thetaf(kzp1))
        allocate(w3(kzp1,kz))
        allocate(cpfac(kz))
        allocate(sdsigma(kz))
        allocate(hbar(kz))
        allocate(hweigh(kz))
        allocate(tbarh(kz))
        allocate(hydroc(kz,kzp1))
        allocate(varpa1(kz,kzp1))
        allocate(hydror(kz,kz))
        allocate(hydros(kz,kz))
        allocate(tau(kz,kz))
        allocate(zmatx(kz,kz))
        allocate(zmatxr(kz,kz))
        allocate(sigmah(kzp1))
        allocate(varpa2(kzp1,kzp1))
        allocate(am(kz,nsplit))
        allocate(an(nsplit))
!
#ifdef MPP1
        allocate(ddsum(iy,jxp,nsplit))
        allocate(deld(iy,jxp,nsplit,3))
        allocate(delh(iy,0:jxp,nsplit,3))
        allocate(dhsum(iy,0:jxp,nsplit))
        allocate(psdot(iy,jxp))
        allocate(work(iy,jxp,3))
        allocate(uu(iy,jxp+1))
        allocate(vv(iy,jxp+1))
        allocate(uuu(iy,kz,jxp+1))
        allocate(vvv(iy,kz,jxp+1))
#else
        allocate(ddsum(iy,jx,nsplit))
        allocate(deld(iy,jx,nsplit,3))
        allocate(delh(iy,jx,nsplit,3))
        allocate(dhsum(iy,jx,nsplit))
        allocate(psdot(iy,jx))
        allocate(work(iy,jx,3))
        allocate(uu(iy,jx))
        allocate(vv(iy,jx))
        allocate(uuu(iy,kz,jx))
        allocate(vvv(iy,kz,jx))
#endif 
        end subroutine allocate_mod_split
!
      subroutine splitf
!
!** compute deld, delh
!** integrate in time and add correction terms appropriately
!
#ifdef MPP1
#ifndef IBM
      use mpi
#else
      include 'mpif.h'
#endif
#endif
      implicit none
!
! Local variables
!
      real(8) :: eps , eps1 , fac , gnuam , gnuan , gnuzm , pdlog , x , &
               & y
      integer :: i , j , k , l , n
      integer :: jm1, jp1
#ifdef MPP1
      integer :: ierr , ii
      real(8) , dimension(iy*nsplit) :: wkrecv , wksend
#endif
!
      do l = 1 , 3
        do n = 1 , nsplit
#ifdef MPP1
          do j = 1 , jendl
            do i = 1 , iy
              deld(i,j,n,l) = 0.
              delh(i,j,n,l) = 0.
            end do
          end do
#else
          do j = 1 , jx
            do i = 1 , iy
              deld(i,j,n,l) = 0.
              delh(i,j,n,l) = 0.
            end do
          end do
#endif
        end do
      end do
!
!**   compute pressure on dot grid
!=======================================================================
!
!     this routine determines p(.) from p(x) by a 4-point interpolation.
!     on the x-grid, a p(x) point outside the grid domain is assumed to
!     satisfy p(0,j)=p(1,j); p(iy,j)=p(iym1,j); and similarly for the
!     i's.
#ifdef MPP1
      call mpi_sendrecv(psa(1,jxp),iy,mpi_real8,ieast,1,                &
                      & psa(1,0),iy,mpi_real8,iwest,1,                  &
                      & mpi_comm_world,mpi_status_ignore,ierr)
#endif
#ifdef MPP1
      do j = jbegin , jendx
#else
#ifdef BAND
      do j = 1 , jx
#else
      do j = 2 , jxm1
#endif
#endif
        jm1 = j-1
#if defined(BAND) && (!defined(MPP1))
        if(jm1.eq.0) jm1=jx
#endif
        do i = 2 , iym1
          psdot(i,j)=0.25*(psa(i,j)+psa(i-1,j)+psa(i,jm1)+psa(i-1,jm1))
        end do
      end do
!
#ifndef BAND
      do i = 2 , iym1
#ifdef MPP1
        if ( myid.eq.0 ) psdot(i,1) = 0.5*(psa(i,1)+psa(i-1,1))
        if ( myid.eq.nproc-1 ) psdot(i,jendl)                           &
           & = 0.5*(psa(i,jendx)+psa(i-1,jendx))
#else
        psdot(i,1) = 0.5*(psa(i,1)+psa(i-1,1))
        psdot(i,jx) = 0.5*(psa(i,jxm1)+psa(i-1,jxm1))
#endif
      end do
#endif
!
#ifdef MPP1
      do j = jbegin , jendx
#else
#ifdef BAND
      do j = 1 , jx
#else
      do j = 2 , jxm1
#endif
#endif
        jm1 = j-1
#if defined(BAND) && (!defined(MPP1))
        if(jm1.eq.0) jm1=jx
#endif
        psdot(1,j) = 0.5*(psa(1,j)+psa(1,jm1))
        psdot(iy,j) = 0.5*(psa(iym1,j)+psa(iym1,jm1))
      end do
!
#ifndef BAND
#ifdef MPP1
      if ( myid.eq.0 ) then
        psdot(1,1) = psa(1,1)
        psdot(iy,1) = psa(iym1,1)
      end if
      if ( myid.eq.nproc-1 ) then
        psdot(1,jendl) = psa(1,jendx)
        psdot(iy,jendl) = psa(iym1,jendx)
      end if
#else
      psdot(1,1) = psa(1,1)
      psdot(iy,1) = psa(iym1,1)
      psdot(1,jx) = psa(1,jxm1)
      psdot(iy,jx) = psa(iym1,jxm1)
#endif
#endif
!
!=======================================================================
!
!**   get deld(0), delh(0) from storage
      do n = 1 , nsplit
#ifdef MPP1
        do j = 1 , jendl
#else
        do j = 1 , jx
#endif
          do i = 1 , iy
            deld(i,j,n,1) = dstor(i,j,n)
            delh(i,j,n,1) = hstor(i,j,n)
          end do
        end do
      end do
!
!=======================================================================
!******* divergence manipulations (f)
      do k = 1 , kz
#ifdef MPP1
        do j = 1 , jendl
#else
        do j = 1 , jx
#endif
          do i = 1 , iy
            uuu(i,k,j) = ua(i,k,j)*msfd(i,j)
            vvv(i,k,j) = va(i,k,j)*msfd(i,j)
          end do
        end do
      end do
#ifdef MPP1
      call mpi_sendrecv(uuu(1,1,1),iy*kz,mpi_real8,iwest,2,             &
                      & uuu(1,1,jxp+1),iy*kz,mpi_real8,ieast,           &
                      & 2,mpi_comm_world,mpi_status_ignore,ierr)
      call mpi_sendrecv(vvv(1,1,1),iy*kz,mpi_real8,iwest,2,             &
                      & vvv(1,1,jxp+1),iy*kz,mpi_real8,ieast,           &
                      & 2,mpi_comm_world,mpi_status_ignore,ierr)
#endif
      do l = 1 , nsplit
#ifdef MPP1
        do j = 1 , jendl
#else
        do j = 1 , jx
#endif
          do i = 1 , iy
            deld(i,j,l,3) = 0.
          end do
        end do

        do k = 1 , kz
#ifdef MPP1
          do j = 1 , jendx
#else
#ifdef BAND
          do j = 1 , jx
#else
          do j = 1 , jxm1
#endif
#endif
            jp1 = j+1
#if defined(BAND) && (!defined(MPP1))
            if(jp1.eq.jx+1) jp1 = 1
#endif
            do i = 1 , iym1
              fac = dx2*msfx(i,j)*msfx(i,j)
              deld(i,j,l,3) = deld(i,j,l,3) + zmatxr(l,k)      &
                   & *(-uuu(i+1,k,j)+uuu(i+1,k,jp1)-uuu(i,k,j) &
                   & +uuu(i,k,jp1)+vvv(i+1,k,j)+vvv(i+1,k,jp1) &
                   & -vvv(i,k,j)-vvv(i,k,jp1))/fac
            end do
          end do
        end do
      end do
!
!=======================================================================
 
      do n = 1 , nsplit
#ifdef MPP1
        do j = 1 , jendl
#else
        do j = 1 , jx
#endif
          do i = 1 , iy
            deld(i,j,n,3) = deld(i,j,n,3) - deld(i,j,n,1)
          end do
        end do
      end do
!
!=======================================================================
!******* divergence manipulations (0)
      do k = 1 , kz
#ifdef MPP1
        do j = 1 , jendl
#else
        do j = 1 , jx
#endif
          do i = 1 , iy
            uuu(i,k,j) = ub(i,k,j)*msfd(i,j)
            vvv(i,k,j) = vb(i,k,j)*msfd(i,j)
          end do
        end do
      end do
#ifdef MPP1
      call mpi_sendrecv(uuu(1,1,1),iy*kz,mpi_real8,iwest,2,             &
                      & uuu(1,1,jxp+1),iy*kz,mpi_real8,ieast,           &
                      & 2,mpi_comm_world,mpi_status_ignore,ierr)
      call mpi_sendrecv(vvv(1,1,1),iy*kz,mpi_real8,iwest,2,             &
                      & vvv(1,1,jxp+1),iy*kz,mpi_real8,ieast,           &
                      & 2,mpi_comm_world,mpi_status_ignore,ierr)
#endif
      do l = 1 , nsplit
#ifdef MPP1
        do j = 1 , jendl
#else
        do j = 1 , jx
#endif
          do i = 1 , iy
            deld(i,j,l,2) = 0.
          end do
        end do
        do k = 1 , kz
#ifdef MPP1
          do j = 1 , jendx
#else
#ifdef BAND
          do j = 1 , jx
#else
          do j = 1 , jxm1
#endif
#endif
            jp1 = j+1
#if defined(BAND) && (!defined(MPP1))
            if(jp1.eq.jx+1) jp1 = 1
#endif
            do i = 1 , iym1
              fac = dx2*msfx(i,j)*msfx(i,j)
              deld(i,j,l,2) = deld(i,j,l,2) + zmatxr(l,k)     &
                  & *(-uuu(i+1,k,j)+uuu(i+1,k,jp1)-uuu(i,k,j) &
                  & +uuu(i,k,jp1)+vvv(i+1,k,j)+vvv(i+1,k,jp1) &
                  & -vvv(i,k,j)-vvv(i,k,jp1))/fac
            end do
          end do
        end do
      end do
!
!=======================================================================
      do n = 1 , nsplit
#ifdef MPP1
        do j = 1 , jendl
#else
        do j = 1 , jx
#endif
          do i = 1 , iy
            deld(i,j,n,1) = deld(i,j,n,1) - deld(i,j,n,2)
          end do
        end do
      end do
!
!=======================================================================
!******* geopotential manipulations (f)
      do l = 1 , nsplit
        pdlog = varpa1(l,kzp1)*dlog(sigmah(kzp1)*pd+r8pt)
        eps1 = varpa1(l,kzp1)*sigmah(kzp1)/(sigmah(kzp1)*pd+r8pt)
#ifdef MPP1
        do j = 1 , jendx
#else
#ifdef BAND
        do j = 1 , jx
#else
        do j = 1 , jxm1
#endif
#endif
          do i = 1 , iym1
            eps = eps1*(psa(i,j)-pd)
            delh(i,j,l,3) = pdlog + eps
          end do
        end do
        do k = 1 , kz
          pdlog = varpa1(l,k)*dlog(sigmah(k)*pd+r8pt)
          eps1 = varpa1(l,k)*sigmah(k)/(sigmah(k)*pd+r8pt)
#ifdef MPP1
          do j = 1 , jendx
#else
#ifdef BAND
          do j = 1 , jx
#else
          do j = 1 , jxm1
#endif
#endif
            do i = 1 , iym1
              eps = eps1*(psa(i,j)-pd)
              delh(i,j,l,3) = delh(i,j,l,3) + pdlog + tau(l,k)*ta(i,k,j)&
                            & /psa(i,j) + eps
            end do
          end do
        end do
      end do
!=======================================================================
 
      do n = 1 , nsplit
#ifdef MPP1
        do j = 1 , jendl
#else
        do j = 1 , jx
#endif
          do i = 1 , iy
            delh(i,j,n,3) = delh(i,j,n,3) - delh(i,j,n,1)
          end do
        end do
      end do
!
!=======================================================================
!******* geopotential manipulations (0)
      do l = 1 , nsplit
        pdlog = varpa1(l,kzp1)*dlog(sigmah(kzp1)*pd+r8pt)
        eps1 = varpa1(l,kzp1)*sigmah(kzp1)/(sigmah(kzp1)*pd+r8pt)
#ifdef MPP1
        do j = 1 , jendx
#else
#ifdef BAND
        do j = 1 , jx
#else
        do j = 1 , jxm1
#endif
#endif
          do i = 1 , iym1
            eps = eps1*(psb(i,j)-pd)
            delh(i,j,l,2) = pdlog + eps
          end do
        end do
        do k = 1 , kz
          pdlog = varpa1(l,k)*dlog(sigmah(k)*pd+r8pt)
          eps1 = varpa1(l,k)*sigmah(k)/(sigmah(k)*pd+r8pt)
#ifdef MPP1
          do j = 1 , jendx
#else
#ifdef BAND
          do j = 1 , jx
#else
          do j = 1 , jxm1
#endif
#endif
            do i = 1 , iym1
              eps = eps1*(psb(i,j)-pd)
              delh(i,j,l,2) = delh(i,j,l,2) + pdlog + tau(l,k)*tb(i,k,j)&
                            & /psb(i,j) + eps
            end do
          end do
        end do
      end do
!=======================================================================
      do n = 1 , nsplit
#ifdef MPP1
        do j = 1 , jendl
#else
        do j = 1 , jx
#endif
          do i = 1 , iy
            delh(i,j,n,1) = delh(i,j,n,1) - delh(i,j,n,2)
          end do
        end do
      end do
!
!**   put deld(0), delh(0) into storage
      do n = 1 , nsplit
#ifdef MPP1
        do j = 1 , jendl
#else
        do j = 1 , jx
#endif
          do i = 1 , iy
            dstor(i,j,n) = deld(i,j,n,2)
            hstor(i,j,n) = delh(i,j,n,2)
          end do
        end do
      end do
!
!******* split explicit time integration
      call spstep(hbar,dx2,dtau,m)
!
!******* add corrections to t and p;  u and v
!=======================================================================
      do l = 1 , nsplit
        gnuan = gnuhf*an(l)
#ifdef MPP1
        do j = jbegin , jendm
#else
#ifdef BAND
        do j = 1 , jx
#else
        do j = 2 , jxm2
#endif
#endif
          do i = 2 , iym2
            psa(i,j) = psa(i,j) - an(l)*ddsum(i,j,l)
            psb(i,j) = psb(i,j) - gnuan*ddsum(i,j,l)
          end do
        end do
      end do
      do l = 1 , nsplit
        do k = 1 , kz
          gnuam = gnuhf*am(k,l)
#ifdef MPP1
          do j = jbegin , jendm
#else
#ifdef BAND
          do j = 1 , jx
#else
          do j = 2 , jxm2
#endif
#endif
            do i = 2 , iym2
              ta(i,k,j) = ta(i,k,j) + am(k,l)*ddsum(i,j,l)
              tb(i,k,j) = tb(i,k,j) + gnuam*ddsum(i,j,l)
            end do
          end do
        end do
      end do
!=======================================================================
#ifdef MPP1
      ii = 0
      do l = 1 , nsplit
        do i = 1 , iy
          ii = ii + 1
          wksend(ii) = dhsum(i,jxp,l)
        end do
      end do
      call mpi_sendrecv(wksend(1),iy*nsplit,mpi_real8,ieast,            &
                      & 1,wkrecv(1),iy*nsplit,mpi_real8,                &
                      & iwest,1,mpi_comm_world,mpi_status_ignore,ierr)
      ii = 0
      do l = 1 , nsplit
        do i = 1 , iy
          ii = ii + 1
          dhsum(i,0,l) = wkrecv(ii)
        end do
      end do
#endif
      do l = 1 , nsplit
        do k = 1 , kz
          gnuzm = gnuhf*zmatx(k,l)
#ifdef MPP1
          do j = jbegin , jendx
#else
#ifdef BAND
          do j = 1 , jx
#else
          do j = 2 , jxm1
#endif
#endif
            jm1 = j-1
#if defined(BAND) && (!defined(MPP1))
            if(jm1.eq.0) jm1 = jx
#endif
            do i = 2 , iym1
              fac = psdot(i,j)/(dx2*msfd(i,j))
              x = fac*(dhsum(i,j,l)+dhsum(i-1,j,l)-dhsum(i,jm1,l) &
                & -dhsum(i-1,jm1,l))
              y = fac*(dhsum(i,j,l)-dhsum(i-1,j,l)+dhsum(i,jm1,l) &
                & -dhsum(i-1,jm1,l))
!
              ua(i,k,j) = ua(i,k,j) - zmatx(k,l)*x
              va(i,k,j) = va(i,k,j) - zmatx(k,l)*y
              ub(i,k,j) = ub(i,k,j) - gnuzm*x
              vb(i,k,j) = vb(i,k,j) - gnuzm*y
            end do
          end do
        end do
      end do
!
!=======================================================================
!
      end subroutine splitf
!
      subroutine spstep(hhbar,dx2,dtau,im)
!
#ifdef MPP1
#ifndef IBM
      use mpi
#else
      include 'mpif.h'
#endif
#endif
      implicit none
!
! Dummy arguments
!
      real(8) :: dx2
      real(8) , dimension(nsplit) :: dtau , hhbar
      integer , dimension(nsplit) :: im
      intent (in) dtau , dx2 , hhbar , im
!
! Local variables
!
      real(8) :: dtau2 , fac
      integer :: i , j , m2 , n , n0 , n1 , n2 , ns , nw
      integer :: jm1, jp1
#ifdef MPP1
      integer :: ierr
      real(8) , dimension(iy*2) :: wkrecv , wksend
#endif
!
      do n = 1 , nsplit
#ifdef MPP1
        do j = 1 , jendl
          do i = 1 , iy
            ddsum(i,j,n) = 0.
            dhsum(i,j,n) = 0.
          end do
        end do
#else
        do j = 1 , jx
          do i = 1 , iy
            ddsum(i,j,n) = 0.
            dhsum(i,j,n) = 0.
          end do
        end do
#endif
      end do
!
      do ns = 1 , nsplit
!
        n0 = 1
        n1 = 2
        n2 = n0
        m2 = im(ns)*2
        dtau2 = dtau(ns)*2.
!
!**     below follows madala(1987)
!c      do 101 j=1,jlx
!c      do 101 i=1,ilx
!       deld, delh: 1,ilx on cross grid
!c      deld(i,j,ns,n1) = deld(i,j,ns,n0)
!c      delh(i,j,ns,n1) = delh(i,j,ns,n0)
!c101   continue
!
#ifdef MPP1
        do j = 1 , jendx
#else
#ifdef BAND
        do j = 1 , jx
#else
        do j = 1 , jxm1
#endif
#endif
          do i = 1 , iym1
!           deld, delh: 1,ilx on cross grid
            ddsum(i,j,ns) = deld(i,j,ns,n0)
            dhsum(i,j,ns) = delh(i,j,ns,n0)
          end do
        end do
!
!**     first step, use forward scheme
!=======================================================================
!
!**     compute gradient of delh;  output = (work1,work2)
!
#ifdef MPP1
        call mpi_sendrecv(delh(1,jxp,ns,n0),iy,mpi_real8,               &
                        & ieast,1,delh(1,0,ns,n0),iy,                   &
                        & mpi_real8,iwest,1,mpi_comm_world,             &
                        & mpi_status_ignore,ierr)
#endif
#ifdef MPP1
        do j = jbegin , jendx
#else
#ifdef BAND
        do j = 1 , jx
#else
        do j = 2 , jxm1
#endif
#endif
          jm1 = j-1
#if defined(BAND) && (!defined(MPP1))
          if(jm1.eq.0) jm1=jx
#endif
          do i = 2 , iym1
            fac = dx2*msfx(i,j)
            work(i,j,1) = (delh(i,j,ns,n0)+delh(i-1,j,ns,n0)  &
                        & -delh(i,jm1,ns,n0)-delh(i-1,jm1,ns,n0))/fac
            work(i,j,2) = (delh(i,j,ns,n0)+delh(i,jm1,ns,n0) &
                        & -delh(i-1,j,ns,n0)-delh(i-1,jm1,ns,n0))/fac
          end do
        end do
!
!=======================================================================
        do nw = 1 , 2
#ifdef MPP1
          do j = jbegin , jendx
#else
#ifdef BAND
          do j = 1 , jx
#else
          do j = 2 , jxm1
#endif
#endif
            do i = 2 , iym1
!             work: 2,ilx on dot grid
              work(i,j,nw) = work(i,j,nw)*psdot(i,j)
            end do
          end do
        end do
!=======================================================================
!
!**     compute divergence z from u and v
!       ( u must be pstar * u ; similarly for v )
!       ( note: map scale factors have been inverted in model (init) )
!
#ifdef MPP1
        do j = jbegin , jendx
#else
#ifdef BAND
        do j = 1 , jx
#else
        do j = 2 , jxm1
#endif
#endif
          do i = 2 , iym1
            uu(i,j) = work(i,j,1)*msfd(i,j)
            vv(i,j) = work(i,j,2)*msfd(i,j)
          end do
        end do
!
#ifdef MPP1
        do i = 1 , iy
          wksend(i) = uu(i,1)
          wksend(i+iy) = vv(i,1)
        end do
        call mpi_sendrecv(wksend(1),2*iy,mpi_real8,iwest,2,             &
                        & wkrecv(1),2*iy,mpi_real8,ieast,2,             &
                        & mpi_comm_world,mpi_status_ignore,ierr)
        do i = 1 , iy
          uu(i,jxp+1) = wkrecv(i)
          vv(i,jxp+1) = wkrecv(i+iy)
        end do
#endif
!
#ifdef MPP1
        do j = jbegin , jendm
#else
#ifdef BAND
        do j = 1 , jx
#else
        do j = 2 , jxm2
#endif
#endif
          jp1 = j+1
#if defined(BAND) && (!defined(MPP1))
          if(jp1.eq.jx+1) jp1=1
#endif
          do i = 2 , iym2
            fac = dx2*msfx(i,j)*msfx(i,j)
            work(i,j,3) = (-uu(i+1,j)+uu(i+1,jp1)-uu(i,j)+uu(i,jp1) &
                         & +vv(i+1,j)+vv(i+1,jp1)-vv(i,j)-vv(i,jp1))/fac
          end do
        end do
!
!=======================================================================
!
#ifdef MPP1
        do j = jbegin , jendm
#else
#ifdef BAND
        do j = 1 , jx
#else
        do j = 2 , jxm2
#endif
#endif
          do i = 2 , iym2
!           work3: 2,iym2 on cross grid
            deld(i,j,ns,n1) = deld(i,j,ns,n0) - dtau(ns)*work(i,j,3)    &
                            & + deld(i,j,ns,3)/m2
            delh(i,j,ns,n1) = delh(i,j,ns,n0) - dtau(ns)*hhbar(ns)      &
                            & *deld(i,j,ns,n0)/psa(i,j) + delh(i,j,ns,3)&
                            & /m2
          end do
        end do
 
!**     not in madala(1987)
        fac = (im(ns)-1.)/im(ns)
#ifndef BAND
        do i = 2 , iym2
#ifdef MPP1
          if ( myid.eq.0 ) delh(i,1,ns,n1) = delh(i,1,ns,n0)*fac
          if ( myid.eq.nproc-1 ) delh(i,jendx,ns,n1)                    &
             & = delh(i,jendx,ns,n0)*fac
#else
          delh(i,1,ns,n1) = delh(i,1,ns,n0)*fac
          delh(i,jxm1,ns,n1) = delh(i,jxm1,ns,n0)*fac
#endif
        end do
#endif
#ifdef MPP1
        do j = 1 , jendx
#else
#ifdef BAND
        do j = 1 , jx
#else
        do j = 1 , jxm1
#endif
#endif
          delh(1,j,ns,n1) = delh(1,j,ns,n0)*fac
          delh(iym1,j,ns,n1) = delh(iym1,j,ns,n0)*fac
        end do
!
#ifdef MPP1
        do j = 1 , jendx
#else
#ifdef BAND
        do j = 1 , jx
#else
        do j = 1 , jxm1
#endif
#endif
          do i = 1 , iym1
            ddsum(i,j,ns) = ddsum(i,j,ns) + deld(i,j,ns,n1)
            dhsum(i,j,ns) = dhsum(i,j,ns) + delh(i,j,ns,n1)
          end do
        end do
!
!**     subsequent steps, use leapfrog scheme
        do n = 2 , m2
!=======================================================================
!
!**       compute gradient of delh;  output = (work1,work2)
!
#ifdef MPP1
          call mpi_sendrecv(delh(1,jxp,ns,n1),iy,mpi_real8,             &
                          & ieast,1,delh(1,0,ns,n1),iy,                 &
                          & mpi_real8,iwest,1,mpi_comm_world,           &
                          & mpi_status_ignore,ierr)
#endif
#ifdef MPP1
          do j = jbegin , jendx
#else
#ifdef BAND
          do j = 1 , jx
#else
          do j = 2 , jxm1
#endif
#endif
            jm1 = j-1
#if defined(BAND) && (!defined(MPP1))
            if(jm1.eq.0) jm1=jx
#endif
            do i = 2 , iym1
              fac = dx2*msfx(i,j)
              work(i,j,1) = (delh(i,j,ns,n1)+delh(i-1,j,ns,n1)  &
                            -delh(i,jm1,ns,n1)-delh(i-1,jm1,ns,n1))/fac
              work(i,j,2) = (delh(i,j,ns,n1)+delh(i,jm1,ns,n1)  &
                            -delh(i-1,j,ns,n1)-delh(i-1,jm1,ns,n1))/fac
            end do
          end do
!=======================================================================
!
          do nw = 1 , 2
#ifdef MPP1
            do j = jbegin , jendx
#else
#ifdef BAND
            do j = 1 , jx
#else
            do j = 2 , jxm1
#endif
#endif
              do i = 2 , iym1
                work(i,j,nw) = work(i,j,nw)*psdot(i,j)
              end do
            end do
          end do
!=======================================================================
!
!**       compute divergence z from u and v
!         ( u must be pstar * u ; similarly for v )
!         ( note: map scale factors have been inverted in model (init) )
!
#ifdef MPP1
          do j = jbegin , jendx
#else
#ifdef BAND
          do j = 1 , jx
#else
          do j = 2 , jxm1
#endif
#endif
            do i = 2 , iym1
              uu(i,j) = work(i,j,1)*msfd(i,j)
              vv(i,j) = work(i,j,2)*msfd(i,j)
            end do
          end do
!
#ifdef MPP1
          do i = 1 , iy
            wksend(i) = uu(i,1)
            wksend(i+iy) = vv(i,1)
          end do
          call mpi_sendrecv(wksend(1),2*iy,mpi_real8,iwest,2,           &
                          & wkrecv(1),2*iy,mpi_real8,ieast,2,           &
                          & mpi_comm_world,mpi_status_ignore,ierr)
          do i = 1 , iy
            uu(i,jxp+1) = wkrecv(i)
            vv(i,jxp+1) = wkrecv(i+iy)
          end do
#endif
!
#ifdef MPP1
          do j = jbegin , jendm
#else
#ifdef BAND
          do j = 1 , jx
#else
          do j = 2 , jxm2
#endif
#endif
            jp1 = j+1
#if defined(BAND) && (!defined(MPP1))
            if(jp1.eq.jx+1) jp1=1
#endif
            do i = 2 , iym2
              fac = dx2*msfx(i,j)*msfx(i,j)
              work(i,j,3) = (-uu(i+1,j)+uu(i+1,jp1)-uu(i,j)+uu(i,jp1) &
                          & +vv(i+1,j)+vv(i+1,jp1)-vv(i,j)-vv(i,jp1)) &
                          & /fac
            end do
          end do
!
!=======================================================================
!
#ifdef MPP1
          do j = jbegin , jendm
#else
#ifdef BAND
          do j = 1 , jx
#else
          do j = 2 , jxm2
#endif
#endif
            do i = 2 , iym2
              deld(i,j,ns,n2) = deld(i,j,ns,n0) - dtau2*work(i,j,3)     &
                              & + deld(i,j,ns,3)/m(ns)
              delh(i,j,ns,n2) = delh(i,j,ns,n0) - dtau2*hhbar(ns)       &
                              & *deld(i,j,ns,n1)/psa(i,j)               &
                              & + delh(i,j,ns,3)/m(ns)
            end do
          end do
!
!**       not in madala(1987)
#ifndef BAND
          do i = 2 , iym2
#ifdef MPP1
            if ( myid.eq.0 ) delh(i,1,ns,n2) = 2.*delh(i,1,ns,n1)       &
               & - delh(i,1,ns,n0)
            if ( myid.eq.nproc-1 ) delh(i,jendx,ns,n2)                  &
               & = 2.*delh(i,jendx,ns,n1) - delh(i,jendx,ns,n0)
#else
            delh(i,1,ns,n2) = 2.*delh(i,1,ns,n1) - delh(i,1,ns,n0)
            delh(i,jxm1,ns,n2) = 2.*delh(i,jxm1,ns,n1) - delh(i,jxm1,ns,n0)
#endif
          end do
#endif
#ifdef MPP1
          do j = 1 , jendx
#else
#ifdef BAND
          do j = 1 , jx
#else
          do j = 1 , jxm1
#endif
#endif
            delh(1,j,ns,n2) = 2.*delh(1,j,ns,n1) - delh(1,j,ns,n0)
            delh(iym1,j,ns,n2) = 2.*delh(iym1,j,ns,n1) - delh(iym1,j,ns,n0)
          end do
!
#ifdef MPP1
          do j = 1 , jendx
#else
#ifdef BAND
          do j = 1 , jx
#else
          do j = 1 , jxm1
#endif
#endif
            do i = 1 , iym1
              ddsum(i,j,ns) = ddsum(i,j,ns) + deld(i,j,ns,n2)
              dhsum(i,j,ns) = dhsum(i,j,ns) + delh(i,j,ns,n2)
            end do
          end do
!
          n0 = n1
          n1 = n2
          n2 = n0
        end do
!
      end do
!
      end subroutine spstep
      end module mod_split

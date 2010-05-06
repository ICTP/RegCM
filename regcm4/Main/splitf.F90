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
 
      subroutine splitf
!
!** compute deld, delh
!** integrate in time and add correction terms appropriately
!
      use mod_dynparam
      use mod_param1
      use mod_main
      use mod_split
      use mod_bxq
      use mod_constants , only : gnuhf
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
      do j = jbegin , jendx
        do i = 2 , iym1
          psdot(i,j) = 0.25*(psa(i,j)+psa(i-1,j)+psa(i,j-1)+psa(i-1,j-1)&
                     & )
        end do
      end do
#else
      do j = 2 , jxm1
        do i = 2 , iym1
          psdot(i,j) = 0.25*(psa(i,j)+psa(i-1,j)+psa(i,j-1)+psa(i-1,j-1)&
                     & )
        end do
      end do
#endif
!
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
!
#ifdef MPP1
      do j = jbegin , jendx
        psdot(1,j) = 0.5*(psa(1,j)+psa(1,j-1))
        psdot(iy,j) = 0.5*(psa(iym1,j)+psa(iym1,j-1))
      end do
#else
      do j = 2 , jxm1
        psdot(1,j) = 0.5*(psa(1,j)+psa(1,j-1))
        psdot(iy,j) = 0.5*(psa(iym1,j)+psa(iym1,j-1))
      end do
#endif
!
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
!
!=======================================================================
!
!**   get deld(0), delh(0) from storage
      do n = 1 , nsplit
#ifdef MPP1
        do j = 1 , jendl
          do i = 1 , iy
            deld(i,j,n,1) = dstor(i,j,n)
            delh(i,j,n,1) = hstor(i,j,n)
          end do
        end do
#else
        do j = 1 , jx
          do i = 1 , iy
            deld(i,j,n,1) = dstor(i,j,n)
            delh(i,j,n,1) = hstor(i,j,n)
          end do
        end do
#endif
      end do
!
!=======================================================================
!******* divergence manipulations (f)
      do k = 1 , kz
#ifdef MPP1
        do j = 1 , jendl
          do i = 1 , iy
            uuu(i,k,j) = ua(i,k,j)*msfd(i,j)
            vvv(i,k,j) = va(i,k,j)*msfd(i,j)
          end do
        end do
#else
        do j = 1 , jx
          do i = 1 , iy
            uuu(i,k,j) = ua(i,k,j)*msfd(i,j)
            vvv(i,k,j) = va(i,k,j)*msfd(i,j)
          end do
        end do
#endif
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
          do i = 1 , iy
            deld(i,j,l,3) = 0.
          end do
        end do
#else
        do j = 1 , jx
          do i = 1 , iy
            deld(i,j,l,3) = 0.
          end do
        end do
#endif
        do k = 1 , kz
#ifdef MPP1
          do j = 1 , jendx
            do i = 1 , iym1
              fac = dx2*msfx(i,j)*msfx(i,j)
              deld(i,j,l,3) = deld(i,j,l,3) + zmatxr(l,k)               &
                            & *(-uuu(i+1,k,j)+uuu(i+1,k,j+1)-uuu(i,k,j) &
                            & +uuu(i,k,j+1)+vvv(i+1,k,j)+vvv(i+1,k,j+1) &
                            & -vvv(i,k,j)-vvv(i,k,j+1))/fac
            end do
          end do
#else
          do j = 1 , jxm1
            do i = 1 , iym1
              fac = dx2*msfx(i,j)*msfx(i,j)
              deld(i,j,l,3) = deld(i,j,l,3) + zmatxr(l,k)               &
                            & *(-uuu(i+1,k,j)+uuu(i+1,k,j+1)-uuu(i,k,j) &
                            & +uuu(i,k,j+1)+vvv(i+1,k,j)+vvv(i+1,k,j+1) &
                            & -vvv(i,k,j)-vvv(i,k,j+1))/fac
            end do
          end do
#endif
        end do
      end do
!
!=======================================================================
 
      do n = 1 , nsplit
#ifdef MPP1
        do j = 1 , jendl
          do i = 1 , iy
            deld(i,j,n,3) = deld(i,j,n,3) - deld(i,j,n,1)
          end do
        end do
#else
        do j = 1 , jx
          do i = 1 , iy
            deld(i,j,n,3) = deld(i,j,n,3) - deld(i,j,n,1)
          end do
        end do
#endif
      end do
!
!=======================================================================
!******* divergence manipulations (0)
      do k = 1 , kz
#ifdef MPP1
        do j = 1 , jendl
          do i = 1 , iy
            uuu(i,k,j) = ub(i,k,j)*msfd(i,j)
            vvv(i,k,j) = vb(i,k,j)*msfd(i,j)
          end do
        end do
#else
        do j = 1 , jx
          do i = 1 , iy
            uuu(i,k,j) = ub(i,k,j)*msfd(i,j)
            vvv(i,k,j) = vb(i,k,j)*msfd(i,j)
          end do
        end do
#endif
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
          do i = 1 , iy
            deld(i,j,l,2) = 0.
          end do
        end do
#else
        do j = 1 , jx
          do i = 1 , iy
            deld(i,j,l,2) = 0.
          end do
        end do
#endif
        do k = 1 , kz
#ifdef MPP1
          do j = 1 , jendx
            do i = 1 , iym1
              fac = dx2*msfx(i,j)*msfx(i,j)
              deld(i,j,l,2) = deld(i,j,l,2) + zmatxr(l,k)               &
                            & *(-uuu(i+1,k,j)+uuu(i+1,k,j+1)-uuu(i,k,j) &
                            & +uuu(i,k,j+1)+vvv(i+1,k,j)+vvv(i+1,k,j+1) &
                            & -vvv(i,k,j)-vvv(i,k,j+1))/fac
            end do
          end do
#else
          do j = 1 , jxm1
            do i = 1 , iym1
              fac = dx2*msfx(i,j)*msfx(i,j)
              deld(i,j,l,2) = deld(i,j,l,2) + zmatxr(l,k)               &
                            & *(-uuu(i+1,k,j)+uuu(i+1,k,j+1)-uuu(i,k,j) &
                            & +uuu(i,k,j+1)+vvv(i+1,k,j)+vvv(i+1,k,j+1) &
                            & -vvv(i,k,j)-vvv(i,k,j+1))/fac
            end do
          end do
#endif
        end do
      end do
!
!=======================================================================
      do n = 1 , nsplit
#ifdef MPP1
        do j = 1 , jendl
          do i = 1 , iy
            deld(i,j,n,1) = deld(i,j,n,1) - deld(i,j,n,2)
          end do
        end do
#else
        do j = 1 , jx
          do i = 1 , iy
            deld(i,j,n,1) = deld(i,j,n,1) - deld(i,j,n,2)
          end do
        end do
#endif
      end do
!
!=======================================================================
!******* geopotential manipulations (f)
      do l = 1 , nsplit
        pdlog = varpa1(l,kzp1)*dlog(sigmah(kzp1)*pd+pt)
        eps1 = varpa1(l,kzp1)*sigmah(kzp1)/(sigmah(kzp1)*pd+pt)
#ifdef MPP1
        do j = 1 , jendx
          do i = 1 , iym1
            eps = eps1*(psa(i,j)-pd)
            delh(i,j,l,3) = pdlog + eps
          end do
        end do
#else
        do j = 1 , jxm1
          do i = 1 , iym1
            eps = eps1*(psa(i,j)-pd)
            delh(i,j,l,3) = pdlog + eps
          end do
        end do
#endif
        do k = 1 , kz
          pdlog = varpa1(l,k)*dlog(sigmah(k)*pd+pt)
          eps1 = varpa1(l,k)*sigmah(k)/(sigmah(k)*pd+pt)
#ifdef MPP1
          do j = 1 , jendx
            do i = 1 , iym1
              eps = eps1*(psa(i,j)-pd)
              delh(i,j,l,3) = delh(i,j,l,3) + pdlog + tau(l,k)*ta(i,k,j)&
                            & /psa(i,j) + eps
            end do
          end do
#else
          do j = 1 , jxm1
            do i = 1 , iym1
              eps = eps1*(psa(i,j)-pd)
              delh(i,j,l,3) = delh(i,j,l,3) + pdlog + tau(l,k)*ta(i,k,j)&
                            & /psa(i,j) + eps
            end do
          end do
#endif
        end do
      end do
!=======================================================================
 
      do n = 1 , nsplit
#ifdef MPP1
        do j = 1 , jendl
          do i = 1 , iy
            delh(i,j,n,3) = delh(i,j,n,3) - delh(i,j,n,1)
          end do
        end do
#else
        do j = 1 , jx
          do i = 1 , iy
            delh(i,j,n,3) = delh(i,j,n,3) - delh(i,j,n,1)
          end do
        end do
#endif
      end do
!
!=======================================================================
!******* geopotential manipulations (0)
      do l = 1 , nsplit
        pdlog = varpa1(l,kzp1)*dlog(sigmah(kzp1)*pd+pt)
        eps1 = varpa1(l,kzp1)*sigmah(kzp1)/(sigmah(kzp1)*pd+pt)
#ifdef MPP1
        do j = 1 , jendx
          do i = 1 , iym1
            eps = eps1*(psb(i,j)-pd)
            delh(i,j,l,2) = pdlog + eps
          end do
        end do
#else
        do j = 1 , jxm1
          do i = 1 , iym1
            eps = eps1*(psb(i,j)-pd)
            delh(i,j,l,2) = pdlog + eps
          end do
        end do
#endif
        do k = 1 , kz
          pdlog = varpa1(l,k)*dlog(sigmah(k)*pd+pt)
          eps1 = varpa1(l,k)*sigmah(k)/(sigmah(k)*pd+pt)
#ifdef MPP1
          do j = 1 , jendx
            do i = 1 , iym1
              eps = eps1*(psb(i,j)-pd)
              delh(i,j,l,2) = delh(i,j,l,2) + pdlog + tau(l,k)*tb(i,k,j)&
                            & /psb(i,j) + eps
            end do
          end do
#else
          do j = 1 , jxm1
            do i = 1 , iym1
              eps = eps1*(psb(i,j)-pd)
              delh(i,j,l,2) = delh(i,j,l,2) + pdlog + tau(l,k)*tb(i,k,j)&
                            & /psb(i,j) + eps
            end do
          end do
#endif
        end do
      end do
!=======================================================================
      do n = 1 , nsplit
#ifdef MPP1
        do j = 1 , jendl
          do i = 1 , iy
            delh(i,j,n,1) = delh(i,j,n,1) - delh(i,j,n,2)
          end do
        end do
#else
        do j = 1 , jx
          do i = 1 , iy
            delh(i,j,n,1) = delh(i,j,n,1) - delh(i,j,n,2)
          end do
        end do
#endif
      end do
!
!**   put deld(0), delh(0) into storage
      do n = 1 , nsplit
#ifdef MPP1
        do j = 1 , jendl
          do i = 1 , iy
            dstor(i,j,n) = deld(i,j,n,2)
            hstor(i,j,n) = delh(i,j,n,2)
          end do
        end do
#else
        do j = 1 , jx
          do i = 1 , iy
            dstor(i,j,n) = deld(i,j,n,2)
            hstor(i,j,n) = delh(i,j,n,2)
          end do
        end do
#endif
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
          do i = 2 , iym2
            psa(i,j) = psa(i,j) - an(l)*ddsum(i,j,l)
            psb(i,j) = psb(i,j) - gnuan*ddsum(i,j,l)
          end do
        end do
#else
        do j = 2 , jxm2
          do i = 2 , iym2
            psa(i,j) = psa(i,j) - an(l)*ddsum(i,j,l)
            psb(i,j) = psb(i,j) - gnuan*ddsum(i,j,l)
          end do
        end do
#endif
      end do
      do l = 1 , nsplit
        do k = 1 , kz
          gnuam = gnuhf*am(k,l)
#ifdef MPP1
          do j = jbegin , jendm
            do i = 2 , iym2
              ta(i,k,j) = ta(i,k,j) + am(k,l)*ddsum(i,j,l)
              tb(i,k,j) = tb(i,k,j) + gnuam*ddsum(i,j,l)
            end do
          end do
#else
          do j = 2 , jxm2
            do i = 2 , iym2
              ta(i,k,j) = ta(i,k,j) + am(k,l)*ddsum(i,j,l)
              tb(i,k,j) = tb(i,k,j) + gnuam*ddsum(i,j,l)
            end do
          end do
#endif
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
            do i = 2 , iym1
              fac = psdot(i,j)/(dx2*msfd(i,j))
              x = fac*(dhsum(i,j,l)+dhsum(i-1,j,l)-dhsum(i,j-1,l)       &
                & -dhsum(i-1,j-1,l))
              y = fac*(dhsum(i,j,l)-dhsum(i-1,j,l)+dhsum(i,j-1,l)       &
                & -dhsum(i-1,j-1,l))
!
              ua(i,k,j) = ua(i,k,j) - zmatx(k,l)*x
              va(i,k,j) = va(i,k,j) - zmatx(k,l)*y
              ub(i,k,j) = ub(i,k,j) - gnuzm*x
              vb(i,k,j) = vb(i,k,j) - gnuzm*y
            end do
          end do
#else
          do j = 2 , jxm1
            do i = 2 , iym1
              fac = psdot(i,j)/(dx2*msfd(i,j))
              x = fac*(dhsum(i,j,l)+dhsum(i-1,j,l)-dhsum(i,j-1,l)       &
                & -dhsum(i-1,j-1,l))
              y = fac*(dhsum(i,j,l)-dhsum(i-1,j,l)+dhsum(i,j-1,l)       &
                & -dhsum(i-1,j-1,l))
!
              ua(i,k,j) = ua(i,k,j) - zmatx(k,l)*x
              va(i,k,j) = va(i,k,j) - zmatx(k,l)*y
              ub(i,k,j) = ub(i,k,j) - gnuzm*x
              vb(i,k,j) = vb(i,k,j) - gnuzm*y
            end do
          end do
#endif
        end do
      end do
!
!=======================================================================
!
      end subroutine splitf

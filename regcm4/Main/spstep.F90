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
 
      subroutine spstep(hbar,dx2,dtau,m)
!
      use mod_dynparam
      use mod_main
      use mod_bxq
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
      real(8) , dimension(nsplit) :: dtau , hbar
      integer , dimension(nsplit) :: m
      intent (in) dtau , dx2 , hbar , m
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
        m2 = m(ns)*2
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
            delh(i,j,ns,n1) = delh(i,j,ns,n0) - dtau(ns)*hbar(ns)       &
                            & *deld(i,j,ns,n0)/psa(i,j) + delh(i,j,ns,3)&
                            & /m2
          end do
        end do
 
!**     not in madala(1987)
        fac = (m(ns)-1.)/m(ns)
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
              delh(i,j,ns,n2) = delh(i,j,ns,n0) - dtau2*hbar(ns)        &
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

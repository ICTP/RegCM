!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
      subroutine tracdiag(xkc)
 
      use mod_regcm_param
      use mod_param1
      use mod_param3 , only : dsigma
      use mod_diagnosis
      use mod_main
      use mod_mainchem
      use mod_constants , only : rgti
#ifdef MPP1
      use mod_mppio
      use mpi
#endif
      implicit none
!
! Dummy arguments
!
#ifdef MPP1
      real(8) , dimension(ix,kx,jxp) :: xkc
#else
      real(8) , dimension(ix,kx,jx) :: xkc
#endif
      intent (in) xkc
!
! Local variables
!
#ifdef MPP1
      integer :: ierr
      real(8) , dimension(kx,ntr,jxp) :: chia01 , chia02 , chiaill ,    &
           & chiaill1
      real(8) , dimension(kx,ntr,jx) :: chia01_g , chia02_g ,           &
           & chiaill1_g , chiaill_g
      real(8) , dimension(jxp) :: psa01 , psa02 , psaill , psaill1
      real(8) , dimension(jx) :: psa01_g , psa02_g , psaill1_g ,        &
                                & psaill_g
      real(8) , dimension(kx,jxp) :: va02 , vaill , xkc02 , xkcill1
      real(8) , dimension(kx,jx) :: va02_g , vaill_g , xkc02_g ,        &
                                   & xkcill1_g
      real(8) :: chid1 , chid2
#endif
      real(8) :: fact1 , fact2 , fx1 , fx2 , uavg1 , uavg2 , vavg1 ,    &
           &  vavg2
      integer :: i , j , k , n
      real(8) , dimension(ixm1,kx,ntr) :: worka , workb
!
!     real(kind=8)  chixp1,chixm1,chiyp1,chiym1,chi00,chidx,chidy
!     real(kind=8)  chixp2,chixm2,chiyp2,chiym2
 
!ccccccccccccccccccccccccccccccccccccccccccccccc
 
!-------------------------
!     1  ADVECTION budgets
!-----------------------
 
!-----advection of tracer through lateral boundaries
 
!
!.....advection through east-west boundaries:
!
!     the 'relaxed' upstream scheme
      fact1 = 0.6
      fact2 = 1 - fact1
 
!     inflow/outflow
#ifdef MPP1
      do n = 1 , ntr
        do k = 1 , kx
          do i = 2 , ixm2
            if ( myid.eq.nproc-1 ) then
              uavg2 = 0.5*(ua(i+1,k,jendx)+ua(i,k,jendx))
              if ( uavg2.lt.0. ) then
                worka(i,k,n) = -uavg2*(fact1*chia(i,k,jendx,n)/psa(i,   &
                             & jendx)/(msfx(i,jendx)*msfx(i,jendx))     &
                             & +fact2*chia(i,k,jendm,n)/psa(i,jendm)    &
                             & /(msfx(i,jendm)*msfx(i,jendm)))
              else
                worka(i,k,n) = -uavg2*(fact1*chia(i,k,jendm,n)/psa(i,   &
                             & jendm)/(msfx(i,jendm)*msfx(i,jendm))     &
                             & +fact2*chia(i,k,jendx,n)/psa(i,jendx)    &
                             & /(msfx(i,jendx)*msfx(i,jendx)))
              end if
            end if
            if ( myid.eq.0 ) then
              uavg1 = 0.5*(ua(i+1,k,1+1)+ua(i,k,1+1))
              if ( uavg1.gt.0. ) then
                workb(i,k,n) = -uavg1*(fact1*chia(i,k,1,n)/psa(i,1)/(   &
                             & msfx(i,1)*msfx(i,1))                     &
                             & +fact2*chia(i,k,1+1,n)/psa(i,1+1)        &
                             & /(msfx(i,1+1)*msfx(i,1+1)))
              else
                workb(i,k,n) = -uavg1*(fact1*chia(i,k,1+1,n)/psa(i,1+1) &
                             & /(msfx(i,1+1)*msfx(i,1+1))               &
                             & +fact2*chia(i,k,1,n)/psa(i,1)            &
                             & /(msfx(i,1)*msfx(i,1)))
              end if
            end if
          end do
        end do
      end do
      call mpi_bcast(worka,ixm1*kx*ntr,mpi_real8,nproc-1,               &
                   & mpi_comm_world,ierr)
#else
      do n = 1 , ntr
        do k = 1 , kx
          do i = 2 , ixm2
            uavg2 = 0.5*(ua(i+1,k,jxm1)+ua(i,k,jxm1))
            if ( uavg2.lt.0. ) then
              worka(i,k,n) = -uavg2*(fact1*chia(i,k,jxm1,n)/psa(i,jxm1) &
                           & /(msfx(i,jxm1)*msfx(i,jxm1))               &
                           & +fact2*chia(i,k,jxm2,n)/psa(i,jxm2)        &
                           & /(msfx(i,jxm2)*msfx(i,jxm2)))
            else
              worka(i,k,n) = -uavg2*(fact1*chia(i,k,jxm2,n)/psa(i,jxm2  &
                           & )/(msfx(i,jxm2)*msfx(i,jxm2))              &
                           & +fact2*chia(i,k,jxm1,n)/psa(i,jxm1)        &
                           & /(msfx(i,jxm1)*msfx(i,jxm1)))
            end if
 
            uavg1 = 0.5*(ua(i+1,k,1+1)+ua(i,k,1+1))
            if ( uavg1.gt.0. ) then
              workb(i,k,n) = -uavg1*(fact1*chia(i,k,1,n)/psa(i,1)/(msfx(&
                           & i,1)*msfx(i,1))+fact2*chia(i,k,1+1,n)      &
                           & /psa(i,1+1)/(msfx(i,1+1)*msfx(i,1+1)))
            else
              workb(i,k,n) = -uavg1*(fact1*chia(i,k,1+1,n)/psa(i,1+1)   &
                           & /(msfx(i,1+1)*msfx(i,1+1))                 &
                           & +fact2*chia(i,k,1,n)/psa(i,1)              &
                           & /(msfx(i,1)*msfx(i,1)))
            end if
          end do
        end do
      end do
#endif 

#ifdef MPP1
      do j = 1 , jendl
        do k = 1 , kx
          vaill(k,j) = va(ixm1,k,j)
          va02(k,j) = va(2,k,j)
          xkcill1(k,j) = xkc(ixm2,k,j)
          xkc02(k,j) = xkc(2,k,j)
          do n = 1 , ntr
            chiaill(k,n,j) = chia(ixm1,k,j,n)
            chiaill1(k,n,j) = chia(ixm2,k,j,n)
            chia01(k,n,j) = chia(1,k,j,n)
            chia02(k,n,j) = chia(2,k,j,n)
          end do
        end do
        psaill(j) = psa(ixm1,j)
        psaill1(j) = psa(ixm2,j)
        psa01(j) = psa(1,j)
        psa02(j) = psa(2,j)
      end do
      call mpi_gather(vaill(1,1),kx*jxp,mpi_real8,                      &
                    & vaill_g(1,1),kx*jxp,mpi_real8,0,                  &
                    & mpi_comm_world,ierr)
      call mpi_gather(va02(1,1),kx*jxp,mpi_real8,va02_g(1,1),           &
                    & kx*jxp,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_gather(xkcill1(1,1),kx*jxp,mpi_real8,                    &
                    & xkcill1_g(1,1),kx*jxp,mpi_real8,0,                &
                    & mpi_comm_world,ierr)
      call mpi_gather(xkc02(1,1),kx*jxp,mpi_real8,                      &
                    & xkc02_g(1,1),kx*jxp,mpi_real8,0,                  &
                    & mpi_comm_world,ierr)
      call mpi_gather(chiaill(1,1,1),kx*ntr*jxp,mpi_real8,              &
                    & chiaill_g(1,1,1),kx*ntr*jxp,mpi_real8,            &
                    & 0,mpi_comm_world,ierr)
      call mpi_gather(chiaill1(1,1,1),kx*ntr*jxp,mpi_real8,             &
                    & chiaill1_g(1,1,1),kx*ntr*jxp,mpi_real8,           &
                    & 0,mpi_comm_world,ierr)
      call mpi_gather(chia01(1,1,1),kx*ntr*jxp,mpi_real8,               &
                    & chia01_g(1,1,1),kx*ntr*jxp,mpi_real8,0,           &
                    & mpi_comm_world,ierr)
      call mpi_gather(chia02(1,1,1),kx*ntr*jxp,mpi_real8,               &
                    & chia02_g(1,1,1),kx*ntr*jxp,mpi_real8,0,           &
                    & mpi_comm_world,ierr)
      call mpi_gather(psaill(1),jxp,mpi_real8,psaill_g(1),              &
                    & jxp,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_gather(psaill1(1),jxp,mpi_real8,psaill1_g(1),            &
                    & jxp,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_gather(psa01(1),jxp,mpi_real8,psa01_g(1),jxp,            &
                    & mpi_real8,0,mpi_comm_world,ierr)
      call mpi_gather(psa02(1),jxp,mpi_real8,psa02_g(1),jxp,            &
                    & mpi_real8,0,mpi_comm_world,ierr)
      if ( myid.eq.0 ) then
        do n = 1 , ntr
          do k = 1 , kx
            do i = 2 , ixm2
              tchiad(n) = tchiad(n) + dtmin*6.E4*dsigma(k)              &
                        & *dx*(worka(i,k,n)-workb(i,k,n))*rgti
            end do
          end do
!.....
!.....advection through north-south boundaries:
!
          do k = 1 , kx
            do j = 2 , jxm2
!hy           inflow/outflow
              vavg2 = 0.5*(vaill_g(k,j+1)+vaill_g(k,j))
              if ( vavg2.lt.0. ) then
                fx2 = -vavg2*(fact1*chiaill_g(k,n,j)/psaill_g(j)        &
                    & /(msfx_io(ixm1,j)*msfx_io(ixm1,j))                &
                    & +fact2*chiaill1_g(k,n,j)/psaill1_g(j)             &
                    & /(msfx_io(ixm2,j)*msfx_io(ixm2,j)))
              else
                fx2 = -vavg2*(fact1*chiaill1_g(k,n,j)/psaill1_g(j)      &
                    & /(msfx_io(ixm2,j)*msfx_io(ixm2,j))                &
                    & +fact2*chiaill_g(k,n,j)/psaill_g(j)               &
                    & /(msfx_io(ixm1,j)*msfx_io(ixm1,j)))
              end if
 
              vavg1 = 0.5*(va02_g(k,j+1)+va02_g(k,j))
              if ( vavg1.gt.0. ) then
                fx1 = -vavg1*(fact1*chia01_g(k,n,j)/psa01_g(j)          &
                    & /(msfx_io(1,j)*msfx_io(1,j))+fact2*chia02_g(k,n,j)&
                    & /psa02_g(j)/(msfx_io(2,j)*msfx_io(2,j)))
              else
                fx1 = -vavg1*(fact1*chia02_g(k,n,j)/psa02_g(j)          &
                    & /(msfx_io(2,j)*msfx_io(2,j))+fact2*chia01_g(k,n,j)&
                    & /psa01_g(j)/(msfx_io(1,j)*msfx_io(1,j)))
              end if
              tchiad(n) = tchiad(n) + dtmin*6.E4*dsigma(k)*dx*(fx2-fx1) &
                        & *rgti
            end do
          end do
        end do
      endif
      call mpi_bcast(tchiad,ntr,mpi_real8,0,mpi_comm_world,ierr)
#else
      do n = 1 , ntr
        do k = 1 , kx
          do i = 2 , ixm2
            tchiad(n) = tchiad(n) + dtmin*6.E4*dsigma(k)                &
                      & *dx*(worka(i,k,n)-workb(i,k,n))/g
          end do
        end do
      end do
!.....
!.....advection through north-south boundaries:
!
      do n = 1 , ntr
        do k = 1 , kx
          do j = 2 , jxm2
!hy         inflow/outflow
            vavg2 = 0.5*(va(ixm1,k,j+1)+va(ixm1,k,j))
            if ( vavg2.lt.0. ) then
              fx2 = -vavg2*(fact1*chia(ixm1,k,j,n)/psa(ixm1,j)          &
                  & /(msfx(ixm1,j)*msfx(ixm1,j))+fact2*chia(ixm2,k,j,n) &
                  & /psa(ixm2,j)/(msfx(ixm2,j)*msfx(ixm2,j)))
            else
              fx2 = -vavg2*(fact1*chia(ixm2,k,j,n)/psa(ixm2,j)          &
                  & /(msfx(ixm2,j)*msfx(ixm2,j))                        &
                  & +fact2*chia(ixm1,k,j,n)/psa(ixm1,j)                 &
                  & /(msfx(ixm1,j)*msfx(ixm1,j)))
            end if
 
            vavg1 = 0.5*(va(1+1,k,j+1)+va(1+1,k,j))
            if ( vavg1.gt.0. ) then
              fx1 = -vavg1*(fact1*chia(1,k,j,n)/psa(1,j)                &
                  & /(msfx(1,j)*msfx(1,j))+fact2*chia(1+1,k,j,n)        &
                  & /psa(1+1,j)/(msfx(1+1,j)*msfx(1+1,j)))
            else
              fx1 = -vavg1*(fact1*chia(1+1,k,j,n)/psa(1+1,j)            &
                  & /(msfx(1+1,j)*msfx(1+1,j))+fact2*chia(1,k,j,n)      &
                  & /psa(1,j)/(msfx(1,j)*msfx(1,j)))
            end if
 
            tchiad(n) = tchiad(n) + dtmin*6.E4*dsigma(k)*dx*(fx2-fx1)/g
 
          end do
        end do
      end do
#endif
 
!
!..... diffusion through east-west boundaries:
!
#ifdef MPP1
      do n = 1 , ntr
        do k = 1 , kx
          do i = 2 , ixm2
            if ( myid.eq.nproc-1 ) worka(i,k,n) = xkc(i,k,jendm)        &
               & *psa(i,jendm)                                          &
               & *(chia(i,k,jendm,n)/psa(i,jendm)-chia(i,k,jendx,n)     &
               & /psa(i,jendx))
            if ( myid.eq.0 ) workb(i,k,n) = xkc(i,k,2)*psa(i,2)         &
               & *(chia(i,k,2,n)/psa(i,2)-chia(i,k,1,n)/psa(i,1))
          end do
        end do
      end do
      call mpi_bcast(worka,ixm1*kx*ntr,mpi_real8,nproc-1,               &
                   & mpi_comm_world,ierr)
#else
      do n = 1 , ntr
        do k = 1 , kx
          do i = 2 , ixm2
            worka(i,k,n) = xkc(i,k,jxm2)*psa(i,jxm2)                    &
                         & *(chia(i,k,jxm2,n)/psa(i,jxm2)               &
                         & -chia(i,k,jxm1,n)/psa(i,jxm1))
            workb(i,k,n) = xkc(i,k,2)*psa(i,2)                          &
                         & *(chia(i,k,2,n)/psa(i,2)-chia(i,k,1,n)       &
                         & /psa(i,1))
          end do
        end do
      end do
#endif

#ifdef MPP1
      if ( myid.eq.0 ) then
        do n = 1 , ntr
          do k = 1 , kx
            do i = 2 , ixm2
              tchitb(n) = tchitb(n) - dtmin*6.E4*dsigma(k)              &
                        & *(workb(i,k,n)+worka(i,k,n))*rgti
            end do
          end do
 
!.....  diffusion through north-south boundaries:
 
          do k = 1 , kx
            do j = 2 , jxm2
              chid1 = xkcill1_g(k,j)*psaill1_g(j)                       &
                    & *(chiaill1_g(k,n,j)/psaill1_g(j)-chiaill_g(k,n,j) &
                    & /psaill_g(j))
              chid2 = xkc02_g(k,j)*psa02_g(j)                           &
                    & *(chia02_g(k,n,j)/psa02_g(j)-chia01_g(k,n,j)      &
                    & /psa01_g(j))
              tchitb(n) = tchitb(n) - dtmin*6.E4*dsigma(k)*(chid2+chid1)&
                        & *rgti
            end do
          end do
        end do
      end if
      call mpi_bcast(tchitb,ntr,mpi_real8,0,mpi_comm_world,ierr)
#else
      do n = 1 , ntr
        do k = 1 , kx
          do i = 2 , ixm2
            tchitb(n) = tchitb(n) - dtmin*6.E4*dsigma(k)                &
                      & *(workb(i,k,n)+worka(i,k,n))*rgti
          end do
        end do
      end do
#endif

      end subroutine tracdiag

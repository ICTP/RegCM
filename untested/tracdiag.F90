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
      use mod_param3
      use mod_diagnosis
      use mod_main
      use mod_mainchem
#ifdef MPP1
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
      real(8) , dimension(kx,ntr,jxp) :: chia01 , chia02 , chiailx ,    &
           & chiailx1
      real(8) , dimension(kx,ntr,mjx) :: chia01_g , chia02_g ,          &
           & chiailx1_g , chiailx_g
      real(8) , dimension(jxp) :: psa01 , psa02 , psailx , psailx1
      real(8) , dimension(mjx) :: psa01_g , psa02_g , psailx1_g ,       &
                                & psailx_g
      real(8) , dimension(kx,jxp) :: va02 , vailx , xkc02 , xkcilx1
      real(8) , dimension(kx,mjx) :: va02_g , vailx_g , xkc02_g ,       &
                                   & xkcilx1_g
      real(8) :: chid1 , chid2
#endif
      real(8) :: fact1 , fact2 , fx1 , fx2 , uavg1 , uavg2 , vavg1 ,    &
           &  vavg2
      integer :: i , j , k , n
      real(8) , dimension(ilx,kx,ntr) :: worka , workb
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
          do i = 1 + 1 , ilx - 1
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
      call mpi_bcast(worka,ilx*kx*ntr,mpi_double_precision,nproc-1,     &
                   & mpi_comm_world,ierr)
#else
      do n = 1 , ntr
        do k = 1 , kx
          do i = 1 + 1 , ilx - 1
            uavg2 = 0.5*(ua(i+1,k,jx-1)+ua(i,k,jx-1))
            if ( uavg2.lt.0. ) then
              worka(i,k,n) = -uavg2*(fact1*chia(i,k,jx-1,n)/psa(i,jx-1) &
                           & /(msfx(i,jx-1)*msfx(i,jx-1))               &
                           & +fact2*chia(i,k,jlx-1,n)/psa(i,jlx-1)      &
                           & /(msfx(i,jlx-1)*msfx(i,jlx-1)))
            else
              worka(i,k,n) = -uavg2*(fact1*chia(i,k,jlx-1,n)/psa(i,jlx-1&
                           & )/(msfx(i,jlx-1)*msfx(i,jlx-1))            &
                           & +fact2*chia(i,k,jx-1,n)/psa(i,jx-1)        &
                           & /(msfx(i,jx-1)*msfx(i,jx-1)))
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
          vailx(k,j) = va(ilx,k,j)
          va02(k,j) = va(2,k,j)
          xkcilx1(k,j) = xkc(ix-2,k,j)
          xkc02(k,j) = xkc(2,k,j)
          do n = 1 , ntr
            chiailx(k,n,j) = chia(ilx,k,j,n)
            chiailx1(k,n,j) = chia(ix-2,k,j,n)
            chia01(k,n,j) = chia(1,k,j,n)
            chia02(k,n,j) = chia(2,k,j,n)
          end do
        end do
        psailx(j) = psa(ilx,j)
        psailx1(j) = psa(ix-2,j)
        psa01(j) = psa(1,j)
        psa02(j) = psa(2,j)
      end do
      call mpi_gather(vailx(1,1),kx*jxp,mpi_double_precision,           &
                    & vailx_g(1,1),kx*jxp,mpi_double_precision,0,       &
                    & mpi_comm_world,ierr)
      call mpi_gather(va02(1,1),kx*jxp,mpi_double_precision,va02_g(1,1),&
                    & kx*jxp,mpi_double_precision,0,mpi_comm_world,ierr)
      call mpi_gather(xkcilx1(1,1),kx*jxp,mpi_double_precision,         &
                    & xkcilx1_g(1,1),kx*jxp,mpi_double_precision,0,     &
                    & mpi_comm_world,ierr)
      call mpi_gather(xkc02(1,1),kx*jxp,mpi_double_precision,           &
                    & xkc02_g(1,1),kx*jxp,mpi_double_precision,0,       &
                    & mpi_comm_world,ierr)
      call mpi_gather(chiailx(1,1,1),kx*ntr*jxp,mpi_double_precision,   &
                    & chiailx_g(1,1,1),kx*ntr*jxp,mpi_double_precision, &
                    & 0,mpi_comm_world,ierr)
      call mpi_gather(chiailx1(1,1,1),kx*ntr*jxp,mpi_double_precision,  &
                    & chiailx1_g(1,1,1),kx*ntr*jxp,mpi_double_precision,&
                    & 0,mpi_comm_world,ierr)
      call mpi_gather(chia01(1,1,1),kx*ntr*jxp,mpi_double_precision,    &
                    & chia01_g(1,1,1),kx*ntr*jxp,mpi_double_precision,0,&
                    & mpi_comm_world,ierr)
      call mpi_gather(chia02(1,1,1),kx*ntr*jxp,mpi_double_precision,    &
                    & chia02_g(1,1,1),kx*ntr*jxp,mpi_double_precision,0,&
                    & mpi_comm_world,ierr)
      call mpi_gather(psailx(1),jxp,mpi_double_precision,psailx_g(1),   &
                    & jxp,mpi_double_precision,0,mpi_comm_world,ierr)
      call mpi_gather(psailx1(1),jxp,mpi_double_precision,psailx1_g(1), &
                    & jxp,mpi_double_precision,0,mpi_comm_world,ierr)
      call mpi_gather(psa01(1),jxp,mpi_double_precision,psa01_g(1),jxp, &
                    & mpi_double_precision,0,mpi_comm_world,ierr)
      call mpi_gather(psa02(1),jxp,mpi_double_precision,psa02_g(1),jxp, &
                    & mpi_double_precision,0,mpi_comm_world,ierr)
      if ( myid.eq.0 ) then
        do n = 1 , ntr
          do k = 1 , kx
            do i = 1 + 1 , ilx - 1
              tchiad(n) = tchiad(n) + dtmin*6.E4*dsigma(k)              &
                        & *dx*(worka(i,k,n)-workb(i,k,n))/g
            end do
          end do
!.....
!.....advection through north-south boundaries:
!
          do k = 1 , kx
            do j = 1 + 1 , mjx - 2
!hy           inflow/outflow
              vavg2 = 0.5*(vailx_g(k,j+1)+vailx_g(k,j))
              if ( vavg2.lt.0. ) then
                fx2 = -vavg2*(fact1*chiailx_g(k,n,j)/psailx_g(j)        &
                    & /(msfx_io(ilx,j)*msfx_io(ilx,j))                  &
                    & +fact2*chiailx1_g(k,n,j)/psailx1_g(j)             &
                    & /(msfx_io(ix-2,j)*msfx_io(ix-2,j)))
              else
                fx2 = -vavg2*(fact1*chiailx1_g(k,n,j)/psailx1_g(j)      &
                    & /(msfx_io(ix-2,j)*msfx_io(ix-2,j))                &
                    & +fact2*chiailx_g(k,n,j)/psailx_g(j)               &
                    & /(msfx_io(ilx,j)*msfx_io(ilx,j)))
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
                        & /g
            end do
          end do
        end do
      endif
      call mpi_bcast(tchiad,ntr,mpi_double_precision,0,mpi_comm_world,  &
                   & ierr)
#else
      do n = 1 , ntr
        do k = 1 , kx
          do i = 1 + 1 , ilx - 1
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
          do j = 1 + 1 , jlx - 1
!hy         inflow/outflow
            vavg2 = 0.5*(va(ix-1,k,j+1)+va(ix-1,k,j))
            if ( vavg2.lt.0. ) then
              fx2 = -vavg2*(fact1*chia(ix-1,k,j,n)/psa(ix-1,j)          &
                  & /(msfx(ix-1,j)*msfx(ix-1,j))+fact2*chia(ilx-1,k,j,n)&
                  & /psa(ilx-1,j)/(msfx(ilx-1,j)*msfx(ilx-1,j)))
            else
              fx2 = -vavg2*(fact1*chia(ilx-1,k,j,n)/psa(ilx-1,j)        &
                  & /(msfx(ilx-1,j)*msfx(ilx-1,j))                      &
                  & +fact2*chia(ix-1,k,j,n)/psa(ix-1,j)                 &
                  & /(msfx(ix-1,j)*msfx(ix-1,j)))
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
          do i = 1 + 1 , ilx - 1
            if ( myid.eq.nproc-1 ) worka(i,k,n) = xkc(i,k,jendm)        &
               & *psa(i,jendm)                                          &
               & *(chia(i,k,jendm,n)/psa(i,jendm)-chia(i,k,jendx,n)     &
               & /psa(i,jendx))
            if ( myid.eq.0 ) workb(i,k,n) = xkc(i,k,2)*psa(i,2)         &
               & *(chia(i,k,2,n)/psa(i,2)-chia(i,k,1,n)/psa(i,1))
          end do
        end do
      end do
      call mpi_bcast(worka,ilx*kx*ntr,mpi_double_precision,nproc-1,     &
                   & mpi_comm_world,ierr)
#else
      do n = 1 , ntr
        do k = 1 , kx
          do i = 1 + 1 , ilx - 1
            worka(i,k,n) = xkc(i,k,jlx-1)*psa(i,jlx-1)                  &
                         & *(chia(i,k,jlx-1,n)/psa(i,jlx-1)             &
                         & -chia(i,k,jlx,n)/psa(i,jlx))
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
            do i = 1 + 1 , ilx - 1
              tchitb(n) = tchitb(n) - dtmin*6.E4*dsigma(k)              &
                        & *(workb(i,k,n)+worka(i,k,n))/g
            end do
          end do
 
!.....  diffusion through north-south boundaries:
 
          do k = 1 , kx
            do j = 1 + 1 , mjx - 2
              chid1 = xkcilx1_g(k,j)*psailx1_g(j)                       &
                    & *(chiailx1_g(k,n,j)/psailx1_g(j)-chiailx_g(k,n,j) &
                    & /psailx_g(j))
              chid2 = xkc02_g(k,j)*psa02_g(j)                           &
                    & *(chia02_g(k,n,j)/psa02_g(j)-chia01_g(k,n,j)      &
                    & /psa01_g(j))
              tchitb(n) = tchitb(n) - dtmin*6.E4*dsigma(k)*(chid2+chid1)&
                        & /g
            end do
          end do
        end do
      end if
      call mpi_bcast(tchitb,ntr,mpi_double_precision,0,mpi_comm_world,  &
                   & ierr)
#else
      do n = 1 , ntr
        do k = 1 , kx
          do i = 1 + 1 , ilx - 1
            tchitb(n) = tchitb(n) - dtmin*6.E4*dsigma(k)                &
                      & *(workb(i,k,n)+worka(i,k,n))/g
          end do
        end do
      end do
#endif

      end subroutine tracdiag

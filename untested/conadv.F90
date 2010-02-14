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
 
      subroutine conadv

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     this subroutine computes the amounts of dry air and water       c
!     substance advected through the lateral boundaries.              c
!                                                                     c
!     ---the unit from advection is converted to "kg".                c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      use mod_regcm_param
      use mod_param1
      use mod_param3
      use mod_diagnosis
      use mod_main
      use mod_constants , only : rgti
#ifdef MPP1
      use mod_mppio
      use mpi
#endif
      implicit none
!
! Local variables
!
#ifdef MPP1
      integer :: ierr
      real(8) , dimension(jxp) :: psa01 , psailx
      real(8) , dimension(jx) :: psa01_g , psailx_g
      real(8) , dimension(kx,jxp) :: qca01 , qcailx , qva01 , qvailx ,  &
                                   & va01 , vaix
      real(8) , dimension(kx,jx) :: qca01_g , qcailx_g , qva01_g ,      &
                                   & qvailx_g , va01_g , vaix_g
#endif
      real(8) , dimension(ixm1,kx) :: worka , workb
      integer :: i ,  j , k
!
!----------------------------------------------------------------------
!-----advection of dry air through the lateral boundaries:
!
!.....advection through east-west boundaries:
!
#ifdef MPP1
      if ( myid.eq.nproc-1 ) then
        do k = 1 , kx
          do i = 1 , ixm1
            worka(i,k) = (ua(i+1,k,jendl)+ua(i,k,jendl))                &
                       & /(msfx(i,jendx)*msfx(i,jendx))
          end do
        end do
      end if
      if ( myid.eq.0 ) then
        do k = 1 , kx
          do i = 1 , ixm1
            workb(i,k) = (ua(i+1,k,1)+ua(i,k,1))/(msfx(i,1)*msfx(i,1))
          end do
        end do
      end if
      call mpi_bcast(worka,ixm1*kx,mpi_real8,nproc-1,                   &
                   & mpi_comm_world,ierr)
      call mpi_bcast(workb,ixm1*kx,mpi_real8,0,                         &
                   & mpi_comm_world,ierr)
#else
      do k = 1 , kx
        do i = 1 , ixm1
          worka(i,k) = (ua(i+1,k,jx)+ua(i,k,jx))                        &
                     & /(msfx(i,jxm1)*msfx(i,jxm1))
          workb(i,k) = (ua(i+1,k,1)+ua(i,k,1))/(msfx(i,1)*msfx(i,1))
        end do
      end do
#endif
      do k = 1 , kx
        do i = 1 , ixm1
          tdadv = tdadv - dtmin*3.E4*dsigma(k)                          &
                & *dx*(worka(i,k)-workb(i,k))*rgti
        end do
      end do
!
!.....advection through north-south boundaries:
!
#ifdef MPP1
      do j = 1 , jendl
        do k = 1 , kx
          vaix(k,j) = va(ix,k,j)
          va01(k,j) = va(1,k,j)
        end do
      end do
      call mpi_gather(vaix(1,1),kx*jxp,mpi_real8,vaix_g(1,1),           &
                    & kx*jxp,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_gather(va01(1,1),kx*jxp,mpi_real8,va01_g(1,1),           &
                    & kx*jxp,mpi_real8,0,mpi_comm_world,ierr)
      if ( myid.eq.0 ) then
        do k = 1 , kx
          do j = 1 , jxm1
            tdadv = tdadv - dtmin*3.E4*dsigma(k)                        &
                  & *dx*((vaix_g(k,j+1)+vaix_g(k,j))                    &
                  & /(msfx_io(ixm1,j)*msfx_io(ixm1,j))                  &
                  & -(va01_g(k,j+1)+va01_g(k,j))                        &
                  & /(msfx_io(1,j)*msfx_io(1,j)))*rgti
          end do
        end do
      end if
      call mpi_bcast(tdadv,1,mpi_real8,0,mpi_comm_world,ierr)
#else
      do k = 1 , kx
        do j = 1 , jxm1
          tdadv = tdadv - dtmin*3.E4*dsigma(k)                          &
                & *dx*((va(ix,k,j+1)+va(ix,k,j))                        &
                & /(msfx(ixm1,j)*msfx(ixm1,j))-(va(1,k,j+1)+va(1,k,j))  &
                & /(msfx(1,j)*msfx(1,j)))*rgti
        end do
      end do
#endif
!
!----------------------------------------------------------------------
!-----advection of water vapor through the lateral boundaries:
!
!.....advection through east-west boundaries:
!
#ifdef MPP1
      if ( myid.eq.nproc-1 ) then
        do k = 1 , kx
          do i = 1 , ixm1
            worka(i,k) = (ua(i+1,k,jendl)+ua(i,k,jendl))                &
                       & *(qva(i,k,jendx)/psa(i,jendx))                 &
                       & /(msfx(i,jendx)*msfx(i,jendx))
          end do
        end do
      end if
      if ( myid.eq.0 ) then
        do k = 1 , kx
          do i = 1 , ixm1
            workb(i,k) = (ua(i+1,k,1)+ua(i,k,1))*(qva(i,k,1)/psa(i,1))  &
                       & /(msfx(i,1)*msfx(i,1))
          end do
        end do
      end if
      call mpi_bcast(worka,ixm1*kx,mpi_real8,nproc-1,                   &
                   & mpi_comm_world,ierr)
      call mpi_bcast(workb,ixm1*kx,mpi_real8,0,                         &
                   & mpi_comm_world,ierr)
#else
      do k = 1 , kx
        do i = 1 , ixm1
          worka(i,k) = (ua(i+1,k,jx)+ua(i,k,jx))                        &
                     & *(qva(i,k,jxm1)/psa(i,jxm1))                     &
                     & /(msfx(i,jxm1)*msfx(i,jxm1))
          workb(i,k) = (ua(i+1,k,1)+ua(i,k,1))*(qva(i,k,1)/psa(i,1))    &
                     & /(msfx(i,1)*msfx(i,1))
        end do
      end do
#endif
      do k = 1 , kx
        do i = 1 , ixm1
          tqadv = tqadv - dtmin*3.E4*dsigma(k)                          &
                & *dx*(worka(i,k)-workb(i,k))*rgti
        end do
      end do
!
!....advection through north-south boundaries:
!
#ifdef MPP1
      do j = 1 , jendl
        do k = 1 , kx
          qvailx(k,j) = qva(ixm1,k,j)
          qva01(k,j) = qva(1,k,j)
        end do
        psailx(j) = psa(ixm1,j)
        psa01(j) = psa(1,j)
      end do
      call mpi_gather(qvailx(1,1),kx*jxp,mpi_real8,                     &
                    & qvailx_g(1,1),kx*jxp,mpi_real8,0,                 &
                    & mpi_comm_world,ierr)
      call mpi_gather(qva01(1,1),kx*jxp,mpi_real8,                      &
                    & qva01_g(1,1),kx*jxp,mpi_real8,0,                  &
                    & mpi_comm_world,ierr)
      call mpi_gather(psailx(1),jxp,mpi_real8,psailx_g(1),              &
                    & jxp,mpi_real8,0,mpi_comm_world,ierr)
      call mpi_gather(psa01(1),jxp,mpi_real8,psa01_g(1),jxp,            &
                    & mpi_real8,0,mpi_comm_world,ierr)
      if ( myid.eq.0 ) then
        do k = 1 , kx
          do j = 1 , jxm1
            tqadv = tqadv - dtmin*3.E4*dsigma(k)                        &
                  & *dx*((vaix_g(k,j+1)+vaix_g(k,j))                    &
                  & *(qvailx_g(k,j)/psailx_g(j))                        &
                  & /(msfx_io(ixm1,j)*msfx_io(ixm1,j))                  &
                  & -(va01_g(k,j+1)+va01_g(k,j))                        &
                  & *(qva01_g(k,j)/psa01_g(j))                          &
                  & /(msfx_io(1,j)*msfx_io(1,j)))*rgti
          end do
        end do
      end if
      call mpi_bcast(tqadv,1,mpi_real8,0,mpi_comm_world,ierr)
#else
      do k = 1 , kx
        do j = 1 , jxm1
          tqadv = tqadv - dtmin*3.E4*dsigma(k)                          &
                & *dx*((va(ix,k,j+1)+va(ix,k,j))                        &
                & *(qva(ixm1,k,j)/psa(ixm1,j))/(msfx(ixm1,j)            &
                & *msfx(ixm1,j))-(va(1,k,j+1)+va(1,k,j))                &
                & *(qva(1,k,j)/psa(1,j))/(msfx(1,j)*msfx(1,j)))*rgti
        end do
      end do
#endif
!
!-----advection of cloud water and rainwater through lateral boundaries:
!
!.....advection through east-west boundaries:
!
#ifdef MPP1
      if ( myid.eq.nproc-1 ) then
        do k = 1 , kx
          do i = 1 , ixm1
            worka(i,k) = (ua(i+1,k,jendl)+ua(i,k,jendl))                &
                       & *(qca(i,k,jendx)/psa(i,jendx))                 &
                       & /(msfx(i,jendx)*msfx(i,jendx))
          end do
        end do
      end if
      if ( myid.eq.0 ) then
        do k = 1 , kx
          do i = 1 , ixm1
            workb(i,k) = (ua(i+1,k,1)+ua(i,k,1))*(qca(i,k,1)/psa(i,1))  &
                       & /(msfx(i,1)*msfx(i,1))
          end do
        end do
      end if
      call mpi_bcast(worka,ixm1*kx,mpi_real8,nproc-1,                   &
                   & mpi_comm_world,ierr)
      call mpi_bcast(workb,ixm1*kx,mpi_real8,0,                         &
                   & mpi_comm_world,ierr)
#else
      do k = 1 , kx
        do i = 1 , ixm1
          worka(i,k) = (ua(i+1,k,jx)+ua(i,k,jx))                        &
                     & *(qca(i,k,jxm1)/psa(i,jxm1))                     &
                     & /(msfx(i,jxm1)*msfx(i,jxm1))
          workb(i,k) = (ua(i+1,k,1)+ua(i,k,1))*(qca(i,k,1)/psa(i,1))    &
                     & /(msfx(i,1)*msfx(i,1))
        end do
      end do
#endif
      do k = 1 , kx
        do i = 1 , ixm1
          tqadv = tqadv - dtmin*3.E4*dsigma(k)                          &
                & *dx*(worka(i,k)-workb(i,k))*rgti
        end do
      end do
!
!....advection through north-south boundaries:
!
#ifdef MPP1
      do j = 1 , jendl
        do k = 1 , kx
          qcailx(k,j) = qca(ixm1,k,j)
          qca01(k,j) = qca(1,k,j)
        end do
      end do
      call mpi_gather(qcailx(1,1),kx*jxp,mpi_real8,                     &
                    & qcailx_g(1,1),kx*jxp,mpi_real8,0,                 &
                    & mpi_comm_world,ierr)
      call mpi_gather(qca01(1,1),kx*jxp,mpi_real8,                      &
                    & qca01_g(1,1),kx*jxp,mpi_real8,0,                  &
                    & mpi_comm_world,ierr)
      if ( myid.eq.0 ) then
        do k = 1 , kx
          do j = 1 , jxm1
            tqadv = tqadv - dtmin*3.E4*dsigma(k)                        &
                  & *dx*((vaix_g(k,j+1)+vaix_g(k,j))                    &
                  & *(qcailx_g(k,j)/psailx_g(j))                        &
                  & /(msfx_io(ixm1,j)*msfx_io(ixm1,j))                  &
                  & -(va01_g(k,j+1)+va01_g(k,j))                        &
                  & *(qca01_g(k,j)/psa01_g(j))                          &
                  & /(msfx_io(1,j)*msfx_io(1,j)))*rgti
          end do
        end do
      end if
      call mpi_bcast(tqadv,1,mpi_real8,0,mpi_comm_world,ierr)
#else
      do k = 1 , kx
        do j = 1 , jxm1
          tqadv = tqadv - dtmin*3.E4*dsigma(k)                          &
                & *dx*((va(ix,k,j+1)+va(ix,k,j))                        &
                & *(qca(ixm1,k,j)/psa(ixm1,j))/(msfx(ixm1,j)            &
                & *msfx(ixm1,j))-(va(1,k,j+1)+va(1,k,j))                &
                & *(qca(1,k,j)/psa(1,j))/(msfx(1,j)*msfx(1,j)))/g
        end do
      end do
#endif
!
      end subroutine conadv

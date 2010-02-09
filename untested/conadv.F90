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

      use regcm_param
      use param1
      use param3
      use diagnosis
      use main
#ifdef MPP1
      use mpi
#endif
      implicit none
!
! Local variables
!
#ifdef MPP1
      integer :: ierr
      real(8) , dimension(jxp) :: psa01 , psailx
      real(8) , dimension(mjx) :: psa01_g , psailx_g
      real(8) , dimension(kx,jxp) :: qca01 , qcailx , qva01 , qvailx ,  &
                                   & va01 , vaix
      real(8) , dimension(kx,mjx) :: qca01_g , qcailx_g , qva01_g ,     &
                                   & qvailx_g , va01_g , vaix_g
      real(8) , dimension(ilx,kx) :: worka , workb
#else
      real(8) , dimension(ilx,kx) :: worka , workb
#endif
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
          do i = 1 , ilx
            worka(i,k) = (ua(i+1,k,jendl)+ua(i,k,jendl))                &
                       & /(msfx(i,jendx)*msfx(i,jendx))
          end do
        end do
      end if
      if ( myid.eq.0 ) then
        do k = 1 , kx
          do i = 1 , ilx
            workb(i,k) = (ua(i+1,k,1)+ua(i,k,1))/(msfx(i,1)*msfx(i,1))
          end do
        end do
      end if
      call mpi_bcast(worka,ilx*kx,mpi_double_precision,nproc-1,         &
                   & mpi_comm_world,ierr)
      call mpi_bcast(workb,ilx*kx,mpi_double_precision,0,mpi_comm_world,&
                   & ierr)
#else
      do k = 1 , kx
        do i = 1 , ilx
          worka(i,k) = (ua(i+1,k,jx)+ua(i,k,jx))                        &
                     & /(msfx(i,jlx)*msfx(i,jlx))
          workb(i,k) = (ua(i+1,k,1)+ua(i,k,1))/(msfx(i,1)*msfx(i,1))
        end do
      end do
#endif
      do k = 1 , kx
        do i = 1 , ilx
          tdadv = tdadv - dtmin*3.E4*dsigma(k)                          &
                & *dx*(worka(i,k)-workb(i,k))/g
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
      call mpi_gather(vaix(1,1),kx*jxp,mpi_double_precision,vaix_g(1,1),&
                    & kx*jxp,mpi_double_precision,0,mpi_comm_world,ierr)
      call mpi_gather(va01(1,1),kx*jxp,mpi_double_precision,va01_g(1,1),&
                    & kx*jxp,mpi_double_precision,0,mpi_comm_world,ierr)
      if ( myid.eq.0 ) then
        do k = 1 , kx
          do j = 1 , mjx - 1
            tdadv = tdadv - dtmin*3.E4*dsigma(k)                        &
                  & *dx*((vaix_g(k,j+1)+vaix_g(k,j))                    &
                  & /(msfx_io(ilx,j)*msfx_io(ilx,j))                    &
                  & -(va01_g(k,j+1)+va01_g(k,j))                        &
                  & /(msfx_io(1,j)*msfx_io(1,j)))/g
          end do
        end do
      end if
      call mpi_bcast(tdadv,1,mpi_double_precision,0,mpi_comm_world,ierr)
#else
      do k = 1 , kx
        do j = 1 , jlx
          tdadv = tdadv - dtmin*3.E4*dsigma(k)                          &
                & *dx*((va(ix,k,j+1)+va(ix,k,j))                        &
                & /(msfx(ilx,j)*msfx(ilx,j))-(va(1,k,j+1)+va(1,k,j))    &
                & /(msfx(1,j)*msfx(1,j)))/g
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
          do i = 1 , ilx
            worka(i,k) = (ua(i+1,k,jendl)+ua(i,k,jendl))                &
                       & *(qva(i,k,jendx)/psa(i,jendx))                 &
                       & /(msfx(i,jendx)*msfx(i,jendx))
          end do
        end do
      end if
      if ( myid.eq.0 ) then
        do k = 1 , kx
          do i = 1 , ilx
            workb(i,k) = (ua(i+1,k,1)+ua(i,k,1))*(qva(i,k,1)/psa(i,1))  &
                       & /(msfx(i,1)*msfx(i,1))
          end do
        end do
      end if
      call mpi_bcast(worka,ilx*kx,mpi_double_precision,nproc-1,         &
                   & mpi_comm_world,ierr)
      call mpi_bcast(workb,ilx*kx,mpi_double_precision,0,mpi_comm_world,&
                   & ierr)
#else
      do k = 1 , kx
        do i = 1 , ilx
          worka(i,k) = (ua(i+1,k,jx)+ua(i,k,jx))                        &
                     & *(qva(i,k,jlx)/psa(i,jlx))                       &
                     & /(msfx(i,jlx)*msfx(i,jlx))
          workb(i,k) = (ua(i+1,k,1)+ua(i,k,1))*(qva(i,k,1)/psa(i,1))    &
                     & /(msfx(i,1)*msfx(i,1))
        end do
      end do
#endif
      do k = 1 , kx
        do i = 1 , ilx
          tqadv = tqadv - dtmin*3.E4*dsigma(k)                          &
                & *dx*(worka(i,k)-workb(i,k))/g
        end do
      end do
!
!....advection through north-south boundaries:
!
#ifdef MPP1
      do j = 1 , jendl
        do k = 1 , kx
          qvailx(k,j) = qva(ilx,k,j)
          qva01(k,j) = qva(1,k,j)
        end do
        psailx(j) = psa(ilx,j)
        psa01(j) = psa(1,j)
      end do
      call mpi_gather(qvailx(1,1),kx*jxp,mpi_double_precision,          &
                    & qvailx_g(1,1),kx*jxp,mpi_double_precision,0,      &
                    & mpi_comm_world,ierr)
      call mpi_gather(qva01(1,1),kx*jxp,mpi_double_precision,           &
                    & qva01_g(1,1),kx*jxp,mpi_double_precision,0,       &
                    & mpi_comm_world,ierr)
      call mpi_gather(psailx(1),jxp,mpi_double_precision,psailx_g(1),   &
                    & jxp,mpi_double_precision,0,mpi_comm_world,ierr)
      call mpi_gather(psa01(1),jxp,mpi_double_precision,psa01_g(1),jxp, &
                    & mpi_double_precision,0,mpi_comm_world,ierr)
      if ( myid.eq.0 ) then
        do k = 1 , kx
          do j = 1 , mjx - 1
            tqadv = tqadv - dtmin*3.E4*dsigma(k)                        &
                  & *dx*((vaix_g(k,j+1)+vaix_g(k,j))                    &
                  & *(qvailx_g(k,j)/psailx_g(j))                        &
                  & /(msfx_io(ilx,j)*msfx_io(ilx,j))                    &
                  & -(va01_g(k,j+1)+va01_g(k,j))                        &
                  & *(qva01_g(k,j)/psa01_g(j))                          &
                  & /(msfx_io(1,j)*msfx_io(1,j)))/g
          end do
        end do
      end if
      call mpi_bcast(tqadv,1,mpi_double_precision,0,mpi_comm_world,ierr)
#else
      do k = 1 , kx
        do j = 1 , jlx
          tqadv = tqadv - dtmin*3.E4*dsigma(k)                          &
                & *dx*((va(ix,k,j+1)+va(ix,k,j))                        &
                & *(qva(ilx,k,j)/psa(ilx,j))/(msfx(ilx,j)*msfx(ilx,j))  &
                & -(va(1,k,j+1)+va(1,k,j))*(qva(1,k,j)/psa(1,j))        &
                & /(msfx(1,j)*msfx(1,j)))/g
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
          do i = 1 , ilx
            worka(i,k) = (ua(i+1,k,jendl)+ua(i,k,jendl))                &
                       & *(qca(i,k,jendx)/psa(i,jendx))                 &
                       & /(msfx(i,jendx)*msfx(i,jendx))
          end do
        end do
      end if
      if ( myid.eq.0 ) then
        do k = 1 , kx
          do i = 1 , ilx
            workb(i,k) = (ua(i+1,k,1)+ua(i,k,1))*(qca(i,k,1)/psa(i,1))  &
                       & /(msfx(i,1)*msfx(i,1))
          end do
        end do
      end if
      call mpi_bcast(worka,ilx*kx,mpi_double_precision,nproc-1,         &
                   & mpi_comm_world,ierr)
      call mpi_bcast(workb,ilx*kx,mpi_double_precision,0,mpi_comm_world,&
                   & ierr)
#else
      do k = 1 , kx
        do i = 1 , ilx
          worka(i,k) = (ua(i+1,k,jx)+ua(i,k,jx))                        &
                     & *(qca(i,k,jlx)/psa(i,jlx))                       &
                     & /(msfx(i,jlx)*msfx(i,jlx))
          workb(i,k) = (ua(i+1,k,1)+ua(i,k,1))*(qca(i,k,1)/psa(i,1))    &
                     & /(msfx(i,1)*msfx(i,1))
        end do
      end do
#endif
      do k = 1 , kx
        do i = 1 , ilx
          tqadv = tqadv - dtmin*3.E4*dsigma(k)                          &
                & *dx*(worka(i,k)-workb(i,k))/g
        end do
      end do
!
!....advection through north-south boundaries:
!
#ifdef MPP1
      do j = 1 , jendl
        do k = 1 , kx
          qcailx(k,j) = qca(ilx,k,j)
          qca01(k,j) = qca(1,k,j)
        end do
      end do
      call mpi_gather(qcailx(1,1),kx*jxp,mpi_double_precision,          &
                    & qcailx_g(1,1),kx*jxp,mpi_double_precision,0,      &
                    & mpi_comm_world,ierr)
      call mpi_gather(qca01(1,1),kx*jxp,mpi_double_precision,           &
                    & qca01_g(1,1),kx*jxp,mpi_double_precision,0,       &
                    & mpi_comm_world,ierr)
      if ( myid.eq.0 ) then
        do k = 1 , kx
          do j = 1 , mjx - 1
            tqadv = tqadv - dtmin*3.E4*dsigma(k)                        &
                  & *dx*((vaix_g(k,j+1)+vaix_g(k,j))                    &
                  & *(qcailx_g(k,j)/psailx_g(j))                        &
                  & /(msfx_io(ilx,j)*msfx_io(ilx,j))                    &
                  & -(va01_g(k,j+1)+va01_g(k,j))                        &
                  & *(qca01_g(k,j)/psa01_g(j))                          &
                  & /(msfx_io(1,j)*msfx_io(1,j)))/g
          end do
        end do
      end if
      call mpi_bcast(tqadv,1,mpi_double_precision,0,mpi_comm_world,ierr)
#else
      do k = 1 , kx
        do j = 1 , jlx
          tqadv = tqadv - dtmin*3.E4*dsigma(k)                          &
                & *dx*((va(ix,k,j+1)+va(ix,k,j))                        &
                & *(qca(ilx,k,j)/psa(ilx,j))/(msfx(ilx,j)*msfx(ilx,j))  &
                & -(va(1,k,j+1)+va(1,k,j))*(qca(1,k,j)/psa(1,j))        &
                & /(msfx(1,j)*msfx(1,j)))/g
        end do
      end do
#endif
!
      end subroutine conadv

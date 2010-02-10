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
 
      subroutine htdiff(dto2,dxsq,akht1)

      use mod_regcm_param
      use mod_pmoist
#ifdef MPP1
      use mpi
#endif
      implicit none
!
! Dummy arguments
!
      real(8) :: akht1 , dto2 , dxsq
      intent (in) akht1 , dto2 , dxsq
!
! Local variables
!
      integer :: i , im1 , ip1 , j , jm1 , jp1 , k
#ifdef MPP1
      integer , dimension(mpi_status_size) :: status
      integer :: ierr
      real(8) , dimension(ix,0:jxp+1) :: wr
#else
      real(8) , dimension(ix,jx) :: wr
#endif
!
#ifdef MPP1
      do k = 1 , kx
        do j = 1 , jendl
          do i = 1 , ix
            wr(i,j) = rsheat(i,k,j)
          end do
        end do
        call mpi_sendrecv(wr(1,jxp),ix,mpi_double_precision,ieast,1,    &
                        & wr(1,0),ix,mpi_double_precision,iwest,1,      &
                        & mpi_comm_world,status,ierr)
        call mpi_sendrecv(wr(1,1),ix,mpi_double_precision,iwest,2,      &
                        & wr(1,jxp+1),ix,mpi_double_precision,ieast,2,  &
                        & mpi_comm_world,status,ierr)
        do j = jbegin , jendm
          if ( myid.eq.0 ) then
            jm1 = max0(j-1,2)
          else
            jm1 = j - 1
          end if
          if ( myid.eq.nproc-1 ) then
            jp1 = min0(j+1,jxp-2)
          else
            jp1 = j + 1
          end if
          do i = 2 , ixm2
            im1 = max0(i-1,2)
            ip1 = min0(i+1,ixm2)
            rsheat(i,k,j) = rsheat(i,k,j)                               &
                          & + akht1*dto2/dxsq*(wr(im1,j)+wr(ip1,j)      &
                          & +wr(i,jm1)+wr(i,jp1)-4.*wr(i,j))
          end do
        end do
      end do
#else
      do k = 1 , kx
        do j = 1 , jx
          do i = 1 , ix
            wr(i,j) = rsheat(i,k,j)
          end do
        end do
        do j = 2 , jxm2
          jm1 = max0(j-1,2)
          jp1 = min0(j+1,jxm2)
          do i = 2 , ixm2
            im1 = max0(i-1,2)
            ip1 = min0(i+1,ixm2)
            rsheat(i,k,j) = rsheat(i,k,j)                               &
                          & + akht1*dto2/dxsq*(wr(im1,j)+wr(ip1,j)      &
                          & +wr(i,jm1)+wr(i,jp1)-4.*wr(i,j))
          end do
        end do
      end do
#endif
!
#ifdef MPP1
      do k = 1 , kx
        do j = 1 , jendl
          do i = 1 , ix
            wr(i,j) = rswat(i,k,j)
          end do
        end do
        call mpi_sendrecv(wr(1,jxp),ix,mpi_double_precision,ieast,1,    &
                        & wr(1,0),ix,mpi_double_precision,iwest,1,      &
                        & mpi_comm_world,status,ierr)
        call mpi_sendrecv(wr(1,1),ix,mpi_double_precision,iwest,2,      &
                        & wr(1,jxp+1),ix,mpi_double_precision,ieast,2,  &
                        & mpi_comm_world,status,ierr)
        do j = jbegin , jendm
          if ( myid.eq.0 ) then
            jm1 = max0(j-1,2)
          else
            jm1 = j - 1
          end if
          if ( myid.eq.nproc-1 ) then
            jp1 = min0(j+1,jxp-2)
          else
            jp1 = j + 1
          end if
          do i = 2 , ixm2
            im1 = max0(i-1,2)
            ip1 = min0(i+1,ixm2)
            rswat(i,k,j) = rswat(i,k,j)                                 &
                         & + akht1*dto2/dxsq*(wr(im1,j)+wr(ip1,j)       &
                         & +wr(i,jm1)+wr(i,jp1)-4.*wr(i,j))
          end do
        end do
      end do
#else
      do k = 1 , kx
        do j = 1 , jx
          do i = 1 , ix
            wr(i,j) = rswat(i,k,j)
          end do
        end do
        do j = 2 , jxm2
          jm1 = max0(j-1,2)
          jp1 = min0(j+1,jxm2)
          do i = 2 , ixm2
            im1 = max0(i-1,2)
            ip1 = min0(i+1,ixm2)
            rswat(i,k,j) = rswat(i,k,j)                                 &
                         & + akht1*dto2/dxsq*(wr(im1,j)+wr(ip1,j)       &
                         & +wr(i,jm1)+wr(i,jp1)-4.*wr(i,j))
          end do
        end do
      end do
#endif
      end subroutine htdiff

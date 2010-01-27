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

      subroutine htdiff(dto2,dxsq,akht1)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
#     include "parameter.inc"
#     include "include/pmoist.cb"
#ifdef MPP1
#     include "mpiregcm.inc"
      integer status(MPI_STATUS_SIZE)
      integer ierr
#endif
      real(kind=8)  dto2,dxsq,akht1
#ifdef MPP1
      real(kind=8)  wr(ix,0:jxp+1)
#else
      real(kind=8)  wr(ix,jx)
#endif
      integer i,j,k,im1,ip1,jm1,jp1
c
      do k=1,kx
#ifdef MPP1
        do j=1,jendl
#else
        do j=1,jx
#endif
          do i=1,ix
            wr(i,j)=rsheat(i,k,j)
          end do
        end do
#ifdef MPP1
        call MPI_SENDRECV(wr(1,jxp),ix,MPI_DOUBLE_PRECISION,ieast,1,
     &                    wr(1,0),ix,MPI_DOUBLE_PRECISION,iwest,1,
     &                    MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(wr(1,1),ix,MPI_DOUBLE_PRECISION,iwest,2,
     &                    wr(1,jxp+1),ix,MPI_DOUBLE_PRECISION,ieast,2,
     &                    MPI_COMM_WORLD,status,ierr)
#endif
#ifdef MPP1
        do j=jbegin,jendm
          if (myid .eq. 0) then
            jm1=max0(j-1,2)
          else
            jm1=j-1
          end if
          if (myid .eq. nproc-1) then
            jp1=min0(j+1,jxp-2)
          else
            jp1=j+1
          end if
#else
        do j=2,jx-2
          jm1=max0(j-1,2)
          jp1=min0(j+1,jx-2)
#endif
          do i=2,ix-2
            im1=max0(i-1,2)
            ip1=min0(i+1,ix-2)
            rsheat(i,k,j)=rsheat(i,k,j)+akht1*dto2/dxsq*
     &         (wr(im1,j)+wr(ip1,j)+wr(i,jm1)+wr(i,jp1)-4.*wr(i,j))
          end do
        end do
      end do
c
      do k=1,kx
#ifdef MPP1
        do j=1,jendl
#else
        do j=1,jx
#endif
          do i=1,ix
            wr(i,j)=rswat(i,k,j)
          end do
        end do
#ifdef MPP1
        call MPI_SENDRECV(wr(1,jxp),ix,MPI_DOUBLE_PRECISION,ieast,1,
     &                    wr(1,0),ix,MPI_DOUBLE_PRECISION,iwest,1,
     &                    MPI_COMM_WORLD,status,ierr)
        call MPI_SENDRECV(wr(1,1),ix,MPI_DOUBLE_PRECISION,iwest,2,
     &                    wr(1,jxp+1),ix,MPI_DOUBLE_PRECISION,ieast,2,
     &                    MPI_COMM_WORLD,status,ierr)
#endif
#ifdef MPP1
        do j=jbegin,jendm
          if (myid .eq. 0) then
            jm1=max0(j-1,2)
          else
            jm1=j-1
          end if
          if (myid .eq. nproc-1) then
            jp1=min0(j+1,jxp-2)
          else
            jp1=j+1
          end if
#else
        do j=2,jx-2
          jm1=max0(j-1,2)
          jp1=min0(j+1,jx-2)
#endif
          do i=2,ix-2
            im1=max0(i-1,2)
            ip1=min0(i+1,ix-2)
            rswat(i,k,j)=rswat(i,k,j)+akht1*dto2/dxsq*
     &         (wr(im1,j)+wr(ip1,j)+wr(i,jm1)+wr(i,jp1)-4.*wr(i,j))
          end do
        end do
      end do
      return
      end subroutine htdiff

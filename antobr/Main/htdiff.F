ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine htdiff(dto2,dxsq,akht1)
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      implicit none
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#ifdef MPP1
      include 'mpiregcm.h'
      include 'Common2/pmoist.cb'
#else
      include 'Commons/pmoist.cb'
#endif
      real*8  dto2,dxsq,akht1
#ifdef MPP1
      real*8  wr(ix,0:jxp+1)
#else
      real*8  wr(ix,jx)
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
          enddo
        enddo
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
          endif
          if (myid .eq. nproc-1) then
            jp1=min0(j+1,jxp-2)
          else
            jp1=j+1
          endif
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
          enddo
        enddo
      enddo
c
      do k=1,kx
#ifdef MPP1
        do j=1,jendl
#else
        do j=1,jx
#endif
          do i=1,ix
            wr(i,j)=rswat(i,k,j)
          enddo
        enddo
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
          endif
          if (myid .eq. nproc-1) then
            jp1=min0(j+1,jxp-2)
          else
            jp1=j+1
          endif
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
          enddo
        enddo
      enddo
      return
      end

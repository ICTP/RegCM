      subroutine invmtrx(a, na, v, nv, n, d, ip, ier, work)
      implicit none
#ifdef MPP1
      include 'mpiregcm.h'
#endif
      integer na, nv, n, ier, info, ip(n)
      real*8  a(n,n), v(n,n), work(n), d(2)
      integer i,j,job
c
c 08/23/91 version 1.0
c 12/10/92 updated to correct bugs
c note: different from cray routine invmtx
c uses subroutines sgefa/sgedi from library linpack
c see dick valent, (SCD, consulting) if problems
c
      if (n .ne. na .or. n .ne. nv) then
        print*, ' stop in valent invmtx: equate n, na, nv'
#ifdef MPP1
        print*, ' stop in subroutine invmtrx'
        call MPI_ABORT(MPI_COMM_WORLD,ierr)
#else
        stop
#endif
      endif
c
      do 10 j=1,n
      do 10 i=1,n
        v(i,j) = a(i,j)
   10 continue
      call sgefa(v,n,n,ip,info)
      if (info .ne. 0) then
        print*, ' stop in valent invmtx: sgefa info=',info
#ifdef MPP1
        print*, ' stop in subroutine invmtrx'
        call MPI_ABORT(MPI_COMM_WORLD,ierr)
#else
        stop
#endif
      endif
      job = 11
      call sgedi(v,n,n,ip,d,work,job)
      ier = info
      return
      end

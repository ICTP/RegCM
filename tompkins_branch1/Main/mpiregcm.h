      include 'mpif.h'
      integer myid,ierr,iwest,ieast
      integer status(MPI_STATUS_SIZE)
      integer jbegin,jendl,jendx,jendm
      common/mpi/ myid,iwest,ieast
      common/jbegin/jbegin,jendl,jendx,jendm

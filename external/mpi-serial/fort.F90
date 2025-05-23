

	subroutine mpi_init(ierror)

        implicit none
        include 'mpif.h'

	integer fint(2)
	logical flog(2)
	real freal(2)
	double precision fdub(2)
	complex fcomp(2)
	integer status(MPI_STATUS_SIZE)

        integer ierror


        !!
        !! Pass values from mpif.h to the C side
        !! to check for consistency mpi.h and hardware sizes.
        !!

        call mpi_init_fort( MPI_COMM_WORLD, &
	                    MPI_ANY_SOURCE, MPI_ANY_TAG, &
                            MPI_COMM_NULL, MPI_REQUEST_NULL, &
                            MPI_GROUP_NULL, MPI_GROUP_EMPTY, &
                            MPI_UNDEFINED, &
                            MPI_MAX_ERROR_STRING, &
                            MPI_MAX_PROCESSOR_NAME, &
                            MPI_STATUS_SIZE, &
                            MPI_SOURCE, MPI_TAG, MPI_ERROR, &
                            status, status(MPI_SOURCE), &
                            status(MPI_TAG), status(MPI_ERROR), &
                            MPI_INTEGER, fint(1), fint(2), &
                            MPI_LOGICAL, flog(1), flog(2), &
                            MPI_REAL, freal(1), freal(2), &
                            MPI_DOUBLE_PRECISION, fdub(1), fdub(2), &
                            MPI_COMPLEX, fcomp(1), fcomp(2), &
                            IERROR )


        return
	end


        subroutine mpi_init_thread(required,provided,ierror)
          implicit none
          integer, intent(in) :: required
          integer, intent(out) :: provided, ierror

          call mpi_init(ierror)
          provided = 0
          return
        end



! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

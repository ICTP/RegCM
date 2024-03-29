
	program test
	implicit none
	use mpi

	integer ier

	integer sreq(10), sreq2(10), rreq(10), rreq2(10)
	integer sbuf(10), sbuf2(10), rbuf(10), rbuf2(10)
	integer tag
	integer status(MPI_STATUS_SIZE,10)
	integer i
	integer comm2;
	logical flag;
	character pname(MPI_MAX_PROCESSOR_NAME)
	integer pnamesize

        integer temp,position


        print *, 'Time=',mpi_wtime()

	call mpi_initialized(flag,ier)
	print *, 'MPI is initialized=',flag

	call mpi_init(ier)

	call mpi_get_processor_name(pname,pnamesize,ier)
	print *, 'proc name: "',pname(1:pnamesize),'"  size:',pnamesize


	call mpi_comm_dup(MPI_COMM_WORLD,comm2,ier)

	call mpi_initialized(flag,ier)
	print *, 'MPI is initialized=',flag




	do i=1,5
	  tag= 100+i
	  print *,  'Post receive tag ',tag

	  call mpi_irecv( rbuf(i),1,MPI_INTEGER,0,tag, &
                          MPI_COMM_WORLD,rreq(i),ier)

	end do
	do i=1,5
!	  tag=1100+i
!	  print *,  'Post receive tag ',tag

	  call mpi_irecv( rbuf2(i),1,MPI_INTEGER, &
                          MPI_ANY_SOURCE, MPI_ANY_TAG, &
                          comm2,rreq2(i),ier)

	end do


	do i=1,5
	  sbuf(i)=10*i
	  tag=100+i
	  print *, 'Send ',sbuf(i),' tag ',tag

	  call mpi_isend( sbuf(i),1,MPI_INTEGER,0,tag, &
	                  MPI_COMM_WORLD,sreq(i),ier)
	end do


	do i=1,5
	  sbuf2(i)=1000+10*i
	  tag=1100+i
	  print *, 'Send ',sbuf2(i),' tag ',tag

	  call mpi_isend( sbuf2(i),1,MPI_INTEGER,0,tag, &
	                  comm2,sreq2(i),ier)
	end do


        print *, 'Time=',mpi_wtime()
	call mpi_waitall(5,sreq,status,ier)
	print *,'sends on MPI_COMM_WORLD done'

	call mpi_waitall(5,rreq,status,ier)
	print *,'recvs on MPI_COMM_WORLD done'

	do i=1,5
          print *, 'Status source=',status(MPI_SOURCE,i), &
                   '  tag=',status(MPI_TAG,i)
	end do

	call mpi_waitall(5,sreq2,status,ier)
	print *,'sends on comm2 done'

	call mpi_waitall(5,rreq2,status,ier)
	print *,'recvs on comm2 done'

	do i=1,5
          print *, 'Status source=',status(MPI_SOURCE,i), &
                   '  tag=',status(MPI_TAG,i)
	end do


! pack/unpack

	position=0
	do i=1,5
          temp=100+i
	  call mpi_pack(temp,1,MPI_INTEGER,sbuf,20,position,MPI_COMM_WORLD,ier)
 	end do

        call mpi_isend(sbuf,position,MPI_PACKED,0,0,MPI_COMM_WORLD,sreq(1),ier)
	call mpi_irecv(rbuf,position,MPI_PACKED,0,0,MPI_COMM_WORLD,rreq(1),ier)
        call mpi_waitall(1,rreq,status,ier)

        print *,"Pack/send/unpack:"

        position=0
	do i=1,5
	  call mpi_unpack( rbuf,20,position,temp,1,MPI_INTEGER, &
                           MPI_COMM_WORLD,ier)
          print *,temp
        end do

!


	call mpi_finalize(ier)

	do i=1,5
          print *, 'Time=',mpi_wtime()
          call sleep(1)
	end do

 	end

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
module mpi

  include 'mpif.h'

  integer mpi_status_ignore(mpi_status_size)
  integer , parameter :: mpi_proc_null = -2

end module mpi

subroutine mpi_sendrecv(sendbuf, sendcount, sendtype, dest, sendtag, &
             recvbuf, recvcount, recvtype, source, recvtag, &
             comm, status, ierror)
  implicit none

  include 'mpif.h'

  real(8) , dimension(:) :: sendbuf , recvbuf
  integer :: sendcount , sendtype , dest , sendtag
  integer :: recvcount , recvtype , source , recvtag , comm
  integer :: status(mpi_status_size)
  integer :: ierror
end subroutine mpi_sendrecv

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
!
module mod_mppparam

  use mpi
  use mod_realkinds
  use mod_dynparam
  use mod_message
  use mod_date

  private

  public :: set_nproc , broadcast_params , date_bcast

  interface master_to_nodes
    module procedure master_to_nodes_r1d
    module procedure master_to_nodes_r2d
    module procedure master_to_nodes_r3d
    module procedure master_to_nodes_r4d
    module procedure master_to_nodes_d1d
    module procedure master_to_nodes_d2d
    module procedure master_to_nodes_d3d
    module procedure master_to_nodes_d4d
    module procedure master_to_nodes_i1d
    module procedure master_to_nodes_i2d
    module procedure master_to_nodes_i3d
    module procedure master_to_nodes_i4d
    module procedure master_to_nodes_l1d
    module procedure master_to_nodes_l2d
    module procedure master_to_nodes_l3d
    module procedure master_to_nodes_l4d
  end interface master_to_nodes

  contains
!
  subroutine set_nproc(ncpu)
    implicit none
    integer , intent(in) :: ncpu
    nproc = ncpu 
    jxp   =  jx/nproc
    jxpsg  = jxp * nsg
  end subroutine set_nproc

  subroutine broadcast_params
    implicit none
    integer :: ierr

    call mpi_barrier(mycomm,ierr)

    call mpi_bcast(iy,1,mpi_integer,0,mycomm,ierr)
    call mpi_bcast(jx,1,mpi_integer,0,mycomm,ierr)
    call mpi_bcast(kz,1,mpi_integer,0,mycomm,ierr)
    call mpi_bcast(nsg,1,mpi_integer,0,mycomm,ierr)
    call mpi_bcast(nveg,1,mpi_integer,0,mycomm,ierr)

    call mpi_bcast(iproj,6,mpi_character,0,mycomm,ierr)
    call mpi_bcast(ds,1,mpi_real8,0,mycomm,ierr)
    call mpi_bcast(ptop,1,mpi_real8,0,mycomm,ierr)
    call mpi_bcast(clat,1,mpi_real8,0,mycomm,ierr)
    call mpi_bcast(clon,1,mpi_real8,0,mycomm,ierr)
    call mpi_bcast(plat,1,mpi_real8,0,mycomm,ierr)
    call mpi_bcast(plon,1,mpi_real8,0,mycomm,ierr)
    call mpi_bcast(truelatl,1,mpi_real8,0,mycomm,ierr)
    call mpi_bcast(truelath,1,mpi_real8,0,mycomm,ierr)
    call mpi_bcast(i_band,1,mpi_integer,0,mycomm,ierr)

    call mpi_bcast(domname,64,mpi_character,0,mycomm,ierr)

    call mpi_bcast(ibyte,1,mpi_integer,0,mycomm,ierr)

    call mpi_bcast(debug_level,1,mpi_integer,0,mycomm,ierr)
    call mpi_bcast(dbgfrq,1,mpi_integer,0,mycomm,ierr)

    call mpi_bcast(nspgx,1,mpi_integer,0,mycomm,ierr)
    call mpi_bcast(nspgd,1,mpi_integer,0,mycomm,ierr)
    call mpi_bcast(high_nudge,1,mpi_real8,0,mycomm,ierr)
    call mpi_bcast(medium_nudge,1,mpi_real8,0,mycomm,ierr)
    call mpi_bcast(low_nudge,1,mpi_real8,0,mycomm,ierr)

    call mpi_bcast(calendar,12,mpi_character,0,mycomm,ierr)
    call mpi_bcast(ical,1,mpi_integer,0,mycomm,ierr)
    call mpi_bcast(dayspy,1,mpi_real8,0,mycomm,ierr)
    call mpi_bcast(dpd,1,mpi_real8,0,mycomm,ierr)

    call mpi_bcast(nsplit,1,mpi_integer,0,mycomm,ierr)

    call mpi_bcast(aertyp,7,mpi_character,0,mycomm,ierr)
    call mpi_bcast(ntr,1,mpi_integer,0,mycomm,ierr)
!    call mpi_bcast(nbin,1,mpi_integer,0,mycomm,ierr)
!    call mpi_bcast(sbin,1,mpi_integer,0,mycomm,ierr)

    call mpi_bcast(ibdyfrq,1,mpi_integer,0,mycomm,ierr)

    ! Setup all convenience dimensions

    if ( myid/= 0) then
      iym1 = iy - 1
      iym2 = iy - 2
      iym3 = iy - 3
      jxp1 = jx + 1
      jxm1 = jx - 1
      jxm2 = jx - 2
      kzm1 = kz - 1
      kzm2 = kz - 2
      kzp1 = kz + 1
      kzp2 = kz + 2
      kzp3 = kz + 3
      kzp4 = kz + 4
      iysg = iy * nsg
      jxsg = jx * nsg
      iym1sg = (iy-1) * nsg
      jxm1sg = (jx-1) * nsg
      iym2sg = (iy-2) * nsg
      jxm2sg = (jx-2) * nsg
      nnsg = nsg*nsg
    end if

    call mpi_barrier(mycomm,ierr)

  end subroutine broadcast_params

  subroutine date_bcast(x,from,comm,ierr)
    type (rcm_time_and_date) , intent(inout) :: x
    integer , intent(in) :: from , comm
    integer , intent(out) :: ierr
    integer :: lerr
    ierr = 0
    call mpi_bcast(x%calendar,1,mpi_integer,from,comm,lerr)
    ierr = ierr+lerr
    call mpi_bcast(x%days_from_reference,1,mpi_integer,from,comm,lerr)
    ierr = ierr+lerr
    call mpi_bcast(x%second_of_day,1,mpi_integer,from,comm,lerr)
    ierr = ierr+lerr
  end subroutine date_bcast

  subroutine master_to_nodes_r1d(mdata,ndata)
    real(sp) , pointer , intent(in) , dimension(:) :: mdata
    real(sp) , pointer , intent(out) , dimension(:) :: ndata
    integer :: ierr , dj , djp
    dj = size(mdata,1)
    djp = size(ndata,1)
    if ( dj /= nproc*djp ) then
      call fatal(__FILE__,__LINE__,'Unmatch dim error')
    end if
    call mpi_scatter(mdata,djp,mpi_real4,ndata,djp,mpi_real4,0,mycomm,ierr)
    if ( ierr /= 0 ) then
      call fatal(__FILE__,__LINE__,'MPI_scatter error')
    end if
  end subroutine master_to_nodes_r1d

  subroutine master_to_nodes_r2d(mdata,ndata)
    real(sp) , pointer , intent(in) , dimension(:,:) :: mdata
    real(sp) , pointer , intent(out) , dimension(:,:) :: ndata
    integer :: ierr , d1 , d2 , dj , djp
    d1 = size(mdata,1)
    d2 = size(ndata,1)
    dj = size(mdata,2)
    djp = size(ndata,2)
    if ( d1 /= d2 .or. dj /= nproc*djp ) then
      call fatal(__FILE__,__LINE__,'Unmatch dim error')
    end if
    call mpi_scatter(mdata,d1*djp,mpi_real4,ndata,d2*djp,mpi_real4,0, &
                     mycomm,ierr)
    if ( ierr /= 0 ) then
      call fatal(__FILE__,__LINE__,'MPI_scatter error')
    end if
  end subroutine master_to_nodes_r2d

  subroutine master_to_nodes_r3d(mdata,ndata)
    real(sp) , pointer , intent(in) , dimension(:,:,:) :: mdata
    real(sp) , pointer , intent(out) , dimension(:,:,:) :: ndata
    integer :: ierr , d1 , d2 , d3 , d4 , dj , djp
    d1 = size(mdata,1)
    d2 = size(ndata,1)
    d3 = size(mdata,2)
    d4 = size(ndata,2)
    dj = size(mdata,3)
    djp = size(ndata,3)
    if ( d1 /= d2 .or. d3 /= d4 .or. dj /= nproc*djp ) then
      call fatal(__FILE__,__LINE__,'Unmatch dim error')
    end if
    call mpi_scatter(mdata,d1*d3*djp,mpi_real4,ndata,d2*d4*djp,mpi_real4,0, &
                     mycomm,ierr)
    if ( ierr /= 0 ) then
      call fatal(__FILE__,__LINE__,'MPI_scatter error')
    end if
  end subroutine master_to_nodes_r3d

  subroutine master_to_nodes_r4d(mdata,ndata)
    real(sp) , pointer , intent(in) , dimension(:,:,:,:) :: mdata
    real(sp) , pointer , intent(out) , dimension(:,:,:,:) :: ndata
    integer :: ierr , d1 , d2 , d3 , d4 , d5 , d6 , dj , djp
    d1 = size(mdata,1)
    d2 = size(ndata,1)
    d3 = size(mdata,2)
    d4 = size(ndata,2)
    d5 = size(mdata,3)
    d6 = size(ndata,3)
    dj = size(mdata,4)
    djp = size(ndata,4)
    if ( d1 /= d2 .or. d3 /= d4 .or. d5 /= d6 .or. dj /= nproc*djp ) then
      call fatal(__FILE__,__LINE__,'Unmatch dim error')
    end if
    call mpi_scatter(mdata,d1*d3*d5*djp,mpi_real4, &
                     ndata,d2*d4*d6*djp,mpi_real4, &
                     0,mycomm,ierr)
    if ( ierr /= 0 ) then
      call fatal(__FILE__,__LINE__,'MPI_scatter error')
    end if
  end subroutine master_to_nodes_r4d

  subroutine master_to_nodes_d1d(mdata,ndata)
    real(dp) , pointer , intent(in) , dimension(:) :: mdata
    real(dp) , pointer , intent(out) , dimension(:) :: ndata
    integer :: ierr , dj , djp
    dj = size(mdata,1)
    djp = size(ndata,1)
    if ( dj /= nproc*djp ) then
      call fatal(__FILE__,__LINE__,'Unmatch dim error')
    end if
    call mpi_scatter(mdata,djp,mpi_real8,ndata,djp,mpi_real8,0,mycomm,ierr)
    if ( ierr /= 0 ) then
      call fatal(__FILE__,__LINE__,'MPI_scatter error')
    end if
  end subroutine master_to_nodes_d1d

  subroutine master_to_nodes_d2d(mdata,ndata)
    real(dp) , pointer , intent(in) , dimension(:,:) :: mdata
    real(dp) , pointer , intent(out) , dimension(:,:) :: ndata
    integer :: ierr , d1 , d2 , dj , djp
    d1 = size(mdata,1)
    d2 = size(ndata,1)
    dj = size(mdata,2)
    djp = size(ndata,2)
    if ( d1 /= d2 .or. dj /= nproc*djp ) then
      call fatal(__FILE__,__LINE__,'Unmatch dim error')
    end if
    call mpi_scatter(mdata,d1*djp,mpi_real8,ndata,d2*djp,mpi_real8,0, &
                     mycomm,ierr)
    if ( ierr /= 0 ) then
      call fatal(__FILE__,__LINE__,'MPI_scatter error')
    end if
  end subroutine master_to_nodes_d2d

  subroutine master_to_nodes_d3d(mdata,ndata)
    real(dp) , pointer , intent(in) , dimension(:,:,:) :: mdata
    real(dp) , pointer , intent(out) , dimension(:,:,:) :: ndata
    integer :: ierr , d1 , d2 , d3 , d4 , dj , djp
    d1 = size(mdata,1)
    d2 = size(ndata,1)
    d3 = size(mdata,2)
    d4 = size(ndata,2)
    dj = size(mdata,3)
    djp = size(ndata,3)
    if ( d1 /= d2 .or. d3 /= d4 .or. dj /= nproc*djp ) then
      call fatal(__FILE__,__LINE__,'Unmatch dim error')
    end if
    call mpi_scatter(mdata,d1*d3*djp,mpi_real8,ndata,d2*d4*djp,mpi_real8,0, &
                     mycomm,ierr)
    if ( ierr /= 0 ) then
      call fatal(__FILE__,__LINE__,'MPI_scatter error')
    end if
  end subroutine master_to_nodes_d3d

  subroutine master_to_nodes_d4d(mdata,ndata)
    real(dp) , pointer , intent(in) , dimension(:,:,:,:) :: mdata
    real(dp) , pointer , intent(out) , dimension(:,:,:,:) :: ndata
    integer :: ierr , d1 , d2 , d3 , d4 , d5 , d6 , dj , djp
    d1 = size(mdata,1)
    d2 = size(ndata,1)
    d3 = size(mdata,2)
    d4 = size(ndata,2)
    d5 = size(mdata,3)
    d6 = size(ndata,3)
    dj = size(mdata,4)
    djp = size(ndata,4)
    if ( d1 /= d2 .or. d3 /= d4 .or. d5 /= d6 .or. dj /= nproc*djp ) then
      call fatal(__FILE__,__LINE__,'Unmatch dim error')
    end if
    call mpi_scatter(mdata,d1*d3*d5*djp,mpi_real8, &
                     ndata,d2*d4*d6*djp,mpi_real8, &
                     0,mycomm,ierr)
    if ( ierr /= 0 ) then
      call fatal(__FILE__,__LINE__,'MPI_scatter error')
    end if
  end subroutine master_to_nodes_d4d

  subroutine master_to_nodes_i1d(mdata,ndata)
    integer , pointer , intent(in) , dimension(:) :: mdata
    integer , pointer , intent(out) , dimension(:) :: ndata
    integer :: ierr , dj , djp
    dj = size(mdata,1)
    djp = size(ndata,1)
    if ( dj /= nproc*djp ) then
      call fatal(__FILE__,__LINE__,'Unmatch dim error')
    end if
    call mpi_scatter(mdata,djp,mpi_integer,ndata,djp,mpi_integer,0,mycomm,ierr)
    if ( ierr /= 0 ) then
      call fatal(__FILE__,__LINE__,'MPI_scatter error')
    end if
  end subroutine master_to_nodes_i1d

  subroutine master_to_nodes_i2d(mdata,ndata)
    integer , pointer , intent(in) , dimension(:,:) :: mdata
    integer , pointer , intent(out) , dimension(:,:) :: ndata
    integer :: ierr , d1 , d2 , dj , djp
    d1 = size(mdata,1)
    d2 = size(ndata,1)
    dj = size(mdata,2)
    djp = size(ndata,2)
    if ( d1 /= d2 .or. dj /= nproc*djp ) then
      call fatal(__FILE__,__LINE__,'Unmatch dim error')
    end if
    call mpi_scatter(mdata,d1*djp,mpi_integer, &
                     ndata,d2*djp,mpi_integer, &
                     0,mycomm,ierr)
    if ( ierr /= 0 ) then
      call fatal(__FILE__,__LINE__,'MPI_scatter error')
    end if
  end subroutine master_to_nodes_i2d

  subroutine master_to_nodes_i3d(mdata,ndata)
    integer , pointer , intent(in) , dimension(:,:,:) :: mdata
    integer , pointer , intent(out) , dimension(:,:,:) :: ndata
    integer :: ierr , d1 , d2 , d3 , d4 , dj , djp
    d1 = size(mdata,1)
    d2 = size(ndata,1)
    d3 = size(mdata,2)
    d4 = size(ndata,2)
    dj = size(mdata,3)
    djp = size(ndata,3)
    if ( d1 /= d2 .or. d3 /= d4 .or. dj /= nproc*djp ) then
      call fatal(__FILE__,__LINE__,'Unmatch dim error')
    end if
    call mpi_scatter(mdata,d1*d3*djp,mpi_integer, &
                     ndata,d2*d4*djp,mpi_integer, &
                     0,mycomm,ierr)
    if ( ierr /= 0 ) then
      call fatal(__FILE__,__LINE__,'MPI_scatter error')
    end if
  end subroutine master_to_nodes_i3d

  subroutine master_to_nodes_i4d(mdata,ndata)
    integer , pointer , intent(in) , dimension(:,:,:,:) :: mdata
    integer , pointer , intent(out) , dimension(:,:,:,:) :: ndata
    integer :: ierr , d1 , d2 , d3 , d4 , d5 , d6 , dj , djp
    d1 = size(mdata,1)
    d2 = size(ndata,1)
    d3 = size(mdata,2)
    d4 = size(ndata,2)
    d5 = size(mdata,3)
    d6 = size(ndata,3)
    dj = size(mdata,4)
    djp = size(ndata,4)
    if ( d1 /= d2 .or. d3 /= d4 .or. d5 /= d6 .or. dj /= nproc*djp ) then
      call fatal(__FILE__,__LINE__,'Unmatch dim error')
    end if
    call mpi_scatter(mdata,d1*d3*d5*djp,mpi_integer, &
                     ndata,d2*d4*d6*djp,mpi_integer, &
                     0,mycomm,ierr)
    if ( ierr /= 0 ) then
      call fatal(__FILE__,__LINE__,'MPI_scatter error')
    end if
  end subroutine master_to_nodes_i4d

  subroutine master_to_nodes_l1d(mdata,ndata)
    logical , pointer , intent(in) , dimension(:) :: mdata
    logical , pointer , intent(out) , dimension(:) :: ndata
    integer :: ierr , dj , djp
    dj = size(mdata,1)
    djp = size(ndata,1)
    if ( dj /= nproc*djp ) then
      call fatal(__FILE__,__LINE__,'Unmatch dim error')
    end if
    call mpi_scatter(mdata,djp,mpi_logical,ndata,djp,mpi_logical,0,mycomm,ierr)
    if ( ierr /= 0 ) then
      call fatal(__FILE__,__LINE__,'MPI_scatter error')
    end if
  end subroutine master_to_nodes_l1d

  subroutine master_to_nodes_l2d(mdata,ndata)
    logical , pointer , intent(in) , dimension(:,:) :: mdata
    logical , pointer , intent(out) , dimension(:,:) :: ndata
    integer :: ierr , d1 , d2 , dj , djp
    d1 = size(mdata,1)
    d2 = size(ndata,1)
    dj = size(mdata,2)
    djp = size(ndata,2)
    if ( d1 /= d2 .or. dj /= nproc*djp ) then
      call fatal(__FILE__,__LINE__,'Unmatch dim error')
    end if
    call mpi_scatter(mdata,d1*djp,mpi_logical, &
                     ndata,d2*djp,mpi_logical, &
                     0,mycomm,ierr)
    if ( ierr /= 0 ) then
      call fatal(__FILE__,__LINE__,'MPI_scatter error')
    end if
  end subroutine master_to_nodes_l2d

  subroutine master_to_nodes_l3d(mdata,ndata)
    logical , pointer , intent(in) , dimension(:,:,:) :: mdata
    logical , pointer , intent(out) , dimension(:,:,:) :: ndata
    integer :: ierr , d1 , d2 , d3 , d4 , dj , djp
    d1 = size(mdata,1)
    d2 = size(ndata,1)
    d3 = size(mdata,2)
    d4 = size(ndata,2)
    dj = size(mdata,3)
    djp = size(ndata,3)
    if ( d1 /= d2 .or. d3 /= d4 .or. dj /= nproc*djp ) then
      call fatal(__FILE__,__LINE__,'Unmatch dim error')
    end if
    call mpi_scatter(mdata,d1*d3*djp,mpi_logical, &
                     ndata,d2*d4*djp,mpi_logical, &
                     0,mycomm,ierr)
    if ( ierr /= 0 ) then
      call fatal(__FILE__,__LINE__,'MPI_scatter error')
    end if
  end subroutine master_to_nodes_l3d

  subroutine master_to_nodes_l4d(mdata,ndata)
    logical , pointer , intent(in) , dimension(:,:,:,:) :: mdata
    logical , pointer , intent(out) , dimension(:,:,:,:) :: ndata
    integer :: ierr , d1 , d2 , d3 , d4 , d5 , d6 , dj , djp
    d1 = size(mdata,1)
    d2 = size(ndata,1)
    d3 = size(mdata,2)
    d4 = size(ndata,2)
    d5 = size(mdata,3)
    d6 = size(ndata,3)
    dj = size(mdata,4)
    djp = size(ndata,4)
    if ( d1 /= d2 .or. d3 /= d4 .or. d5 /= d6 .or. dj /= nproc*djp ) then
      call fatal(__FILE__,__LINE__,'Unmatch dim error')
    end if
    call mpi_scatter(mdata,d1*d3*d5*djp,mpi_logical, &
                     ndata,d2*d4*d6*djp,mpi_logical, &
                     0,mycomm,ierr)
    if ( ierr /= 0 ) then
      call fatal(__FILE__,__LINE__,'MPI_scatter error')
    end if
  end subroutine master_to_nodes_l4d

end module mod_mppparam

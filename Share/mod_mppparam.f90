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
  use mod_dynparam
  use mod_date

  private

  public :: set_nproc , broadcast_params , date_bcast

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

    call mpi_bcast(ehso4,1,mpi_logical,0,mycomm,ierr)

    call mpi_bcast(aertyp,7,mpi_character,0,mycomm,ierr)
    call mpi_bcast(ntr,1,mpi_integer,0,mycomm,ierr)
    call mpi_bcast(nbin,1,mpi_integer,0,mycomm,ierr)
    call mpi_bcast(sbin,1,mpi_integer,0,mycomm,ierr)

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

end module mod_mppparam

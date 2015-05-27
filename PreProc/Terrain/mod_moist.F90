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

module mod_moist
!
  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_interp
  use mod_stdio
  use netcdf
!
  private
!
  public :: read_moist

  contains

  subroutine read_moist(rmoist,smoist,jx,iy)
    implicit none
    real(rk8) , pointer , dimension(:,:) , intent(inout) :: rmoist
    real(rk8) , pointer , dimension(:,:) , intent(inout) :: smoist
    integer(ik4) , intent(in) :: jx , iy
    integer(ik4) :: ncid , istat
    integer(ik4) :: dimid , njx , niy , ntime
    integer(ik4) :: varid
    integer(ik4) :: istart(4) , icount(4)
    real(rk8) , dimension(:,:,:) , allocatable :: moist_in

    istat = nf90_open('moist.nc',nf90_nowrite,ncid)
    if ( istat /= nf90_noerr ) then
      write(stdout,*) 'Skip moisture creation, no moist.nc file to read'
      return
    end if

    istat = nf90_inq_dimid(ncid,'jx',dimid)
    if ( istat /= nf90_noerr ) then
      ! Check if coming from cdo
      istat = nf90_inq_dimid(ncid,'x',dimid)
      if ( istat /= nf90_noerr ) then
        write(stderr,*) 'No jx dimension in file moist.nc'
        return
      end if
    end if
    istat = nf90_inquire_dimension(ncid,dimid,len=njx)
    if ( istat /= nf90_noerr ) then
      write(stderr,*) 'Error reading dimension x in file moist.nc'
      return
    end if
    if ( njx /= jx-3 ) then
      write(stderr,*) 'File moist.nc is not coming from a previous run'
      write(stderr,*) 'JX : ', jx , njx+3
      return
    end if
    istat = nf90_inq_dimid(ncid,'iy',dimid)
    if ( istat /= nf90_noerr ) then
      ! Check if coming from cdo
      istat = nf90_inq_dimid(ncid,'y',dimid)
      if ( istat /= nf90_noerr ) then
        write(stderr,*) 'No iy dimension in file moist.nc'
        return
      end if
    end if
    istat = nf90_inquire_dimension(ncid,dimid,len=niy)
    if ( istat /= nf90_noerr ) then
      write(stderr,*) 'Error reading dimension y in file moist.nc'
      return
    end if
    if ( niy /= iy-3 ) then
      write(stderr,*) 'File moist.nc is not coming from a previous run'
      write(stderr,*) 'IY : ', iy , niy+3
      return
    end if

    allocate(moist_in(njx,niy,2))
    istat = nf90_inq_varid(ncid,'mrso',varid)
    if ( istat /= nf90_noerr ) then
      write(stderr,*) 'Error finding variable mrso in file moist.nc'
      return
    end if

    istat = nf90_inq_dimid(ncid,'time',dimid)
    if ( istat /= nf90_noerr ) then
      istat = nf90_get_var(ncid,varid,moist_in)
      if ( istat /= nf90_noerr ) then
        write(stderr,*) 'Error reading variable mrso in file moist.nc'
        return
      end if
    else
      ! Read last timestep
      istat = nf90_inquire_dimension(ncid,dimid,len=ntime)
      if ( istat /= nf90_noerr ) then
        write(stderr,*) 'Error reading dimension time in file moist.nc'
        return
      end if
      istart(1) = 1
      icount(1) = njx
      istart(2) = 1
      icount(2) = niy
      istart(3) = 1
      icount(3) = 2
      istart(4) = ntime
      icount(4) = 1
      istat = nf90_get_var(ncid,varid,moist_in,istart,icount)
      if ( istat /= nf90_noerr ) then
        write(stderr,*) 'Error reading variable moist in file moist.nc'
        return
      end if
    end if

    smoist(2:jx-2,2:iy-2) = moist_in(:,:,1)
    rmoist(2:jx-2,2:iy-2) = moist_in(:,:,2)

    write(stdout,*) 'Read initial soil moisture from previous run SRF file'
    deallocate(moist_in)

  end subroutine read_moist
!
!
end module mod_moist
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

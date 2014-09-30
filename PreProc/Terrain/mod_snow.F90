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

module mod_snow
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
  public :: read_snow

  contains

  subroutine read_snow(snow,jx,iy)
    implicit none
    real(rk8) , pointer , dimension(:,:) , intent(inout) :: snow
    integer(ik4) , intent(in) :: jx , iy
    integer(ik4) :: ncid , istat
    integer(ik4) :: dimid , njx , niy , ntime
    integer(ik4) :: varid
    integer(ik4) :: istart(3) , icount(3)
    real(rk8) , dimension(:,:) , allocatable :: snow_in

    istat = nf90_open('snow.nc',nf90_nowrite,ncid)
    if ( istat /= nf90_noerr ) then
      write(stderr,*) 'No snow.nc file to read'
      return
    end if

    istat = nf90_inq_dimid(ncid,'jx',dimid)
    if ( istat /= nf90_noerr ) then
      ! Check if coming from cdo
      istat = nf90_inq_dimid(ncid,'x',dimid)
      if ( istat /= nf90_noerr ) then
        write(stderr,*) 'No jx dimension in file snow.nc'
        return
      end if
    end if
    istat = nf90_inquire_dimension(ncid,dimid,len=njx)
    if ( istat /= nf90_noerr ) then
      write(stderr,*) 'Error reading dimension x in file snow.nc'
      return
    end if
    if ( njx /= jx-3 ) then
      write(stderr,*) 'File snow.nc is not coming from a previous run'
      write(stderr,*) 'JX : ', jx , njx+3
      return
    end if
    istat = nf90_inq_dimid(ncid,'iy',dimid)
    if ( istat /= nf90_noerr ) then
      ! Check if coming from cdo
      istat = nf90_inq_dimid(ncid,'y',dimid)
      if ( istat /= nf90_noerr ) then
        write(stderr,*) 'No iy dimension in file snow.nc'
        return
      end if
    end if
    istat = nf90_inquire_dimension(ncid,dimid,len=niy)
    if ( istat /= nf90_noerr ) then
      write(stderr,*) 'Error reading dimension y in file snow.nc'
      return
    end if
    if ( niy /= iy-3 ) then
      write(stderr,*) 'File snow.nc is not coming from a previous run'
      write(stderr,*) 'IY : ', iy , niy+3
      return
    end if

    allocate(snow_in(njx,niy))
    istat = nf90_inq_varid(ncid,'snv',varid)
    if ( istat /= nf90_noerr ) then
      write(stderr,*) 'Error finding variable snv in file snow.nc'
      return
    end if

    istat = nf90_inq_dimid(ncid,'time',dimid)
    if ( istat /= nf90_noerr ) then
      istat = nf90_get_var(ncid,varid,snow_in)
      if ( istat /= nf90_noerr ) then
        write(stderr,*) 'Error reading variable snv in file snow.nc'
        return
      end if
    else
      ! Read last timestep
      istat = nf90_inquire_dimension(ncid,dimid,len=ntime)
      if ( istat /= nf90_noerr ) then
        write(stderr,*) 'Error reading dimension time in file snow.nc'
        return
      end if
      istart(1) = 1
      icount(1) = njx
      istart(2) = 1
      icount(2) = niy
      istart(3) = ntime
      icount(3) = 1
      istat = nf90_get_var(ncid,varid,snow_in,istart,icount)
      if ( istat /= nf90_noerr ) then
        write(stderr,*) 'Error reading variable snv in file snow.nc'
        return
      end if
    end if

    snow(2:jx-2,2:iy-2) = snow_in

    write(stdout,*) 'Read initial snow from previous run SRF file'
    deallocate(snow_in)

  end subroutine read_snow
!
!
end module mod_snow

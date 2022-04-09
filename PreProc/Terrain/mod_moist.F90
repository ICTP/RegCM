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

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_stdio
  use netcdf

  private

  public :: read_moist

  contains

  subroutine read_moist(fname,rmoist,snow,jx,iy,ns,lrmoist)
    implicit none
    character(len=256) , intent(in) :: fname
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: rmoist
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: snow
    integer(ik4) , intent(in) :: jx , iy , ns
    logical , intent(out) :: lrmoist
    integer(ik4) :: ncid , istat
    integer(ik4) :: dimid , njx , niy , nlev , ntime
    integer(ik4) :: idmoist , idsnow
    integer(ik4) :: istart(4) , icount(4)
    real(rkx) , dimension(:,:,:) , allocatable :: moist_in
    real(rkx) , dimension(:,:) , allocatable :: snow_in
    logical :: lsnw , lmrso

    lrmoist = .false.
    lsnw = .false.
    lmrso = .false.

    istat = nf90_open(fname,nf90_nowrite,ncid)
    if ( istat /= nf90_noerr ) then
      write(stdout,*) 'Skip moisture creation, no moisture file to read'
      return
    end if

    istat = nf90_inq_dimid(ncid,'jx',dimid)
    if ( istat /= nf90_noerr ) then
      ! Check if coming from cdo
      istat = nf90_inq_dimid(ncid,'x',dimid)
      if ( istat /= nf90_noerr ) then
        write(stderr,*) 'No jx dimension in moisture file'
        return
      end if
    end if
    istat = nf90_inquire_dimension(ncid,dimid,len=njx)
    if ( istat /= nf90_noerr ) then
      write(stderr,*) 'Error reading dimension x in moisture file'
      return
    end if
    if ( njx /= jx-3 ) then
      write(stderr,*) 'Moisture file is not coming from a previous run'
      write(stderr,*) 'JX : ', jx , njx+3
      return
    end if
    istat = nf90_inq_dimid(ncid,'iy',dimid)
    if ( istat /= nf90_noerr ) then
      ! Check if coming from cdo
      istat = nf90_inq_dimid(ncid,'y',dimid)
      if ( istat /= nf90_noerr ) then
        write(stderr,*) 'No iy dimension in moisture file'
      end if
    end if
    istat = nf90_inquire_dimension(ncid,dimid,len=niy)
    if ( istat /= nf90_noerr ) then
      write(stderr,*) 'Error reading dimension y in moisture file'
      return
    end if
    if ( niy /= iy-3 ) then
      write(stderr,*) 'Moisture file is not coming from a previous run'
      write(stderr,*) 'IY : ', iy , niy+3
      return
    end if

    istat = nf90_inq_dimid(ncid,'soil_layer',dimid)
    if ( istat /= nf90_noerr ) then
      write(stderr,*) 'No soil_layer dimension in moisture file'
      return
    end if
    istat = nf90_inquire_dimension(ncid,dimid,len=nlev)
    if ( istat /= nf90_noerr ) then
      write(stderr,*) 'Error reading dimension soil_layer in moisture file'
      return
    end if

    if ( ns /= nlev ) then
      write(stderr,*) 'Number of soil layers in moisture file inconsistent!'
      write(stderr,*) nlev, ' /= ', ns
      return
    end if

    allocate(moist_in(njx,niy,nlev))
    allocate(snow_in(njx,niy))

    istat = nf90_inq_varid(ncid,'mrso',idmoist)
    if ( istat /= nf90_noerr ) then
      istat = nf90_inq_varid(ncid,'mrsos',idmoist)
      if ( istat /= nf90_noerr ) then
        write(stderr,*) 'Error finding variable mrso in moisture file'
        lmrso = .false.
      else
        lmrso = .true.
      end if
    else
      lmrso = .true.
    end if
    istat = nf90_inq_varid(ncid,'snw',idsnow)
    if ( istat /= nf90_noerr ) then
      write(stderr,*) 'Error finding variable snw in moisture file'
      lsnw = .false.
    else
      lsnw = .true.
    end if

    if ( lmrso .or. lsnw ) then
      istat = nf90_inq_dimid(ncid,'time',dimid)
      if ( istat /= nf90_noerr ) then
        if ( lmrso ) then
          istat = nf90_get_var(ncid,idmoist,moist_in)
          if ( istat /= nf90_noerr ) then
            write(stderr,*) 'Error reading variable mrso in moisture file'
          else
            rmoist(2:jx-2,2:iy-2,:) = moist_in(:,:,:)
          end if
        end if
        if ( lsnw ) then
          istat = nf90_get_var(ncid,idsnow,snow_in)
          if ( istat /= nf90_noerr ) then
             write(stderr,*) 'Error reading variable snv in moisture file'
          else
            snow(2:jx-2,2:iy-2) = snow_in(:,:)
          end if
        end if
      else
        ! Read last timestep
        istat = nf90_inquire_dimension(ncid,dimid,len=ntime)
        if ( istat /= nf90_noerr ) then
         write(stderr,*) 'Error reading dimension time in moisture file'
          return
        end if
        istart(1) = 1
        icount(1) = njx
        istart(2) = 1
        icount(2) = niy
        istart(3) = 1
        icount(3) = nlev
        istart(4) = ntime
        icount(4) = 1
        if ( lmrso ) then
          istat = nf90_get_var(ncid,idmoist,moist_in,istart,icount)
          if ( istat /= nf90_noerr ) then
            write(stderr,*) 'Error reading variable mrso in moisture file'
          else
            rmoist(2:jx-2,2:iy-2,:) = moist_in(:,:,:)
          end if
        end if
        istart(1) = 1
        icount(1) = njx
        istart(2) = 1
        icount(2) = niy
        istart(3) = ntime
        icount(3) = 1
        if ( lsnw ) then
          istat = nf90_get_var(ncid,idsnow,snow_in,istart(1:3),icount(1:3))
          if ( istat /= nf90_noerr ) then
            write(stderr,*) 'Error reading variable snv in moisture file'
          else
            snow(2:jx-2,2:iy-2) = snow_in(:,:)
          end if
        end if
      end if
    end if

    deallocate(moist_in)
    deallocate(snow_in)

    lrmoist = .true.

  end subroutine read_moist

end module mod_moist

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

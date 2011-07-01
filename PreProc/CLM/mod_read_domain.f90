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

module mod_read_domain

  use netcdf
  use m_realkinds
  use mod_dynparam
  use mod_memutil
  use m_stdio
  use m_die

  private

  real(sp) , public :: clatx , clonx
  real(sp) :: dsx , grdfacx , platx , plonx , ptopx
  integer :: iyy , jxx , kzz

  character(6) , public :: iprojx

  real(sp) , public , pointer , dimension(:,:) :: xlat, xlon
  real(sp) , public , pointer , dimension(:) :: xlat1d
  real(sp) , public , pointer , dimension(:) :: xlon1d
  real(sp) , public , pointer , dimension(:) :: sigx

  public :: init_domain , read_domain

  contains

  subroutine init_domain
    implicit none
    call getmem2d(xlat,1,iy,1,jx,'mod_read_domain:xlat')
    call getmem2d(xlon,1,iy,1,jx,'mod_read_domain:xlon')
    call getmem1d(xlat1d,1,iy,'mod_read_domain:xlat1d')
    call getmem1d(xlon1d,1,jx,'mod_read_domain:xlon1d')
    call getmem1d(sigx,1,kzp1,'mod_read_domain:sigx')
  end subroutine init_domain

  subroutine read_domain(terfile)
    implicit none
    character(256) , intent(in) :: terfile
    integer :: istatus
    integer :: incin
    integer :: idimid
    integer :: ivarid
    real(sp) , pointer , dimension(:,:) :: xlat_dum, xlon_dum

    !Open the netcdf file
    call handle_nc_err( nf90_open(terfile, nf90_nowrite, incin),   &
      'Opening', trim(terfile))

    !Read the dimensions from the netcdf file
    call handle_nc_err( nf90_inq_dimid(incin,'iy',idimid),'Finding','iy')
    call handle_nc_err( nf90_inquire_dimension(incin, idimid, len=iyy),'Reading','iy')

    call handle_nc_err( nf90_inq_dimid(incin,'jx',idimid),'Finding','jx')
    call handle_nc_err( nf90_inquire_dimension(incin,idimid,len=jxx),'Reading','jx')

    call handle_nc_err( nf90_inq_dimid(incin,'kz',idimid),'Finding','kz')
    call handle_nc_err( nf90_inquire_dimension(incin,idimid,len=kzz),'Reading','kz')

    !Check for consistency with regcm.in
    if ( iyy/=iy .or. jxx/=jx .or. kzz/=kz+1 ) then
      write(stderr,*) 'DOMAIN.INFO is inconsistent with regcm.in'
      write(stderr,*) '  namelist  : ' , iy , jx , kz
      write(stderr,*) '  DOMAIN    : ' , iyy , jxx , kzz
      call die('read_domain','Dimension mismatch',1)
    end if

    call getmem2d(xlat_dum,1,jx,1,iy,'mod_read_domain:xlat_dum')
    call getmem2d(xlon_dum,1,jx,1,iy,'mod_read_domain:xlon_dum')

    !Read ds
    call handle_nc_err(  &
       nf90_get_att(incin, NF90_GLOBAL,'grid_size_in_meters', dsx),   &
      'Reading','ds')
    !Convert from m to km
    dsx = dsx/1000

    !Read clatx
    call handle_nc_err(  &
      nf90_get_att(incin, NF90_GLOBAL,'latitude_of_projection_origin' &
                   ,clatx),&
      'Reading','latitude_of_projection_origin')

    !Read clonx
    call handle_nc_err(  &
      nf90_get_att(incin, NF90_GLOBAL,'longitude_of_projection_origin', &
                   clonx),&
      'Reading','longitude_of_projection_origin')

    !Read iproj
    call handle_nc_err(  &
       nf90_get_att(incin, NF90_GLOBAL,'projection', iprojx), &
      'Reading','projection')

    !Only if using the Rotated Mercator projection, read the poles
    if(iprojx.eq.'ROTMER')then
      !Read plat
      call handle_nc_err(  &
         nf90_get_att(incin, NF90_GLOBAL,'grid_north_pole_latitude', &
                      platx), &
        'Reading','grid_north_pole_latitude')
      !Read plon
      call handle_nc_err(  &
         nf90_get_att(incin, NF90_GLOBAL,'grid_north_pole_longitude', &
                      plonx), &
        'Reading','grid_north_pole_longitude')
    else
      platx = clatx
      plonx = clonx
    endif

    !Read grdfacx
    call handle_nc_err(  &
       nf90_get_att(incin, NF90_GLOBAL,'grid_factor',grdfacx), &
      'Reading','projection')

    !Read sigx
    call handle_nc_err(  &
      nf90_inq_varid(incin,'sigma',ivarid), &
      'Finding','sigma')
    call handle_nc_err(  &
      nf90_get_var(incin,ivarid,sigx), &
      'Reading','sigma')

    !Read ptopx
    call handle_nc_err(  &
      nf90_inq_varid(incin,'ptop',ivarid), &
      'Finding','ptop')
    call handle_nc_err(  &
      nf90_get_var(incin,ivarid,ptopx), &
      'Reading','ptop')

    !Read ptopx
    call handle_nc_err(  &
      nf90_inq_varid(incin,'ptop',ivarid), &
      'Finding','ptop')
    call handle_nc_err(  &
      nf90_get_var(incin,ivarid,ptopx), &
      'Reading','ptop')

    !Read xlat
    call handle_nc_err(  &
      nf90_inq_varid(incin,'xlat',ivarid), &
      'Finding','xlat')
    call handle_nc_err(  &
      nf90_get_var(incin,ivarid,xlat_dum), &
      'Reading','xlat')
      
    !Read xlon
    call handle_nc_err(  &
      nf90_inq_varid(incin,'xlon',ivarid), &
      'Finding','xlon')
    call handle_nc_err(  &
      nf90_get_var(incin,ivarid,xlon_dum), &
      'Reading','xlon')

    !Set xlat and xlon, swapping the i/j indicies of what was read in
    xlat = transpose(xlat_dum)
    xlon = transpose(xlon_dum)

    !Close the netcdf flie
    istatus = nf90_close(incin)
    if (istatus /= nf90_noerr) then
      write(stderr,*) 'Error closing Domain file ', trim(terfile)
      write(stderr,*) nf90_strerror(istatus)
      call die('read_domain')
    end if

  end subroutine read_domain

  !A routine to simply handle NetCDF errors
  subroutine handle_nc_err(incerr,sAction,sVarname)
    implicit none
    integer,intent(in) :: incerr
    character*(*),intent(in) :: sAction,sVarname

     if (incerr /= nf90_noerr) then
      write(stderr,*) 'Error associated with ', trim(sAction), &
                      ' ',trim(sVarname)
      write(stderr,*) nf90_strerror(incerr)
      call die('read_domain')
    end if

  end subroutine handle_nc_err

end module mod_read_domain

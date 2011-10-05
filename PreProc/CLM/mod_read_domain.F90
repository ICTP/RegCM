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
  use mod_realkinds
  use mod_dynparam
  use mod_memutil
  use mod_message
  use mod_stdio

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
    integer :: incin
    integer :: idimid
    integer :: ivarid
    integer :: istatus
    real(sp) , pointer , dimension(:,:) :: xlat_dum, xlon_dum

    !Open the netcdf file
    istatus = nf90_open(terfile, nf90_nowrite, incin)
    call checkncerr(istatus,__FILE__,__LINE__,'Opening '//trim(terfile))

    !Read the dimensions from the netcdf file
    istatus = nf90_inq_dimid(incin,'iy',idimid)
    call checkncerr(istatus,__FILE__,__LINE__,'Finding dim iy')
    istatus = nf90_inquire_dimension(incin, idimid, len=iyy)
    call checkncerr(istatus,__FILE__,__LINE__,'Inquire dim iy')

    istatus = nf90_inq_dimid(incin,'jx',idimid)
    call checkncerr(istatus,__FILE__,__LINE__,'Finding dim jx')
    istatus = nf90_inquire_dimension(incin,idimid,len=jxx)
    call checkncerr(istatus,__FILE__,__LINE__,'Inquire dim jx')

    istatus = nf90_inq_dimid(incin,'kz',idimid)
    call checkncerr(istatus,__FILE__,__LINE__,'Finding dim kz')
    istatus = nf90_inquire_dimension(incin,idimid,len=kzz)
    call checkncerr(istatus,__FILE__,__LINE__,'Inquire dim kz')

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
    istatus = nf90_get_att(incin,NF90_GLOBAL,'grid_size_in_meters', dsx)
    call checkncerr(istatus,__FILE__,__LINE__,'Read att grid_size_in_meters')
    !Convert from m to km
    dsx = dsx/1000

    !Read clatx
    istatus = nf90_get_att(incin,NF90_GLOBAL,'latitude_of_projection_origin',clatx)
    call checkncerr(istatus,__FILE__,__LINE__,'Read att latitude_of_projection_origin')

    !Read clonx
    istatus = nf90_get_att(incin,NF90_GLOBAL,'longitude_of_projection_origin',clonx)
    call checkncerr(istatus,__FILE__,__LINE__,'Read att longitude_of_projection_origin')

    !Read iproj
    istatus = nf90_get_att(incin,NF90_GLOBAL,'projection',iprojx)
    call checkncerr(istatus,__FILE__,__LINE__,'Read att projection')

    !Only if using the Rotated Mercator projection, read the poles
    if(iprojx.eq.'ROTMER')then
      !Read plat
      istatus = nf90_get_att(incin,NF90_GLOBAL,'grid_north_pole_latitude',platx)
      call checkncerr(istatus,__FILE__,__LINE__,'Read att grid_north_pole_latitude')
      !Read plon
      istatus = nf90_get_att(incin,NF90_GLOBAL,'grid_north_pole_longitude',plonx)
      call checkncerr(istatus,__FILE__,__LINE__,'Read att grid_north_pole_longitude')
    else
      platx = clatx
      plonx = clonx
    endif

    !Read grdfacx
    istatus = nf90_get_att(incin,NF90_GLOBAL,'grid_factor',grdfacx)
    call checkncerr(istatus,__FILE__,__LINE__,'Read att grid_factor')

    !Read sigx
    istatus = nf90_inq_varid(incin,'sigma',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Find var sigma')
    istatus = nf90_get_var(incin,ivarid,sigx)
    call checkncerr(istatus,__FILE__,__LINE__,'Read var sigma')

    !Read ptopx
    istatus = nf90_inq_varid(incin,'ptop',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Find var ptop')
    istatus = nf90_get_var(incin,ivarid,ptopx)
    call checkncerr(istatus,__FILE__,__LINE__,'Read var ptop')

    !Read xlat
    istatus = nf90_inq_varid(incin,'xlat',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Find var xlat')
    istatus = nf90_get_var(incin,ivarid,xlat_dum)
    call checkncerr(istatus,__FILE__,__LINE__,'Read var xlat')
      
    !Read xlon
    istatus = nf90_inq_varid(incin,'xlon',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Find var xlon')
    istatus = nf90_get_var(incin,ivarid,xlon_dum)
    call checkncerr(istatus,__FILE__,__LINE__,'Read var xlon')

    !Set xlat and xlon, swapping the i/j indicies of what was read in
    xlat = transpose(xlat_dum)
    xlon = transpose(xlon_dum)

    !Close the netcdf flie
    istatus = nf90_close(incin)
    call checkncerr(istatus,__FILE__,__LINE__,'Error close '//trim(terfile))

  end subroutine read_domain

end module mod_read_domain

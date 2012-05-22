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

module mod_domain

  use netcdf
  use mod_realkinds
  use mod_dynparam
  use mod_memutil
  use mod_nchelper
  use mod_stdio
  use mod_message

  private

  type domain_io
    real(dp) , pointer , dimension(:) :: sigma
    real(dp) , pointer , dimension(:,:) :: xlat
    real(dp) , pointer , dimension(:,:) :: xlon
    real(dp) , pointer , dimension(:,:) :: dlat
    real(dp) , pointer , dimension(:,:) :: dlon
    real(dp) , pointer , dimension(:,:) :: ht
    real(dp) , pointer , dimension(:,:) :: mask
    real(dp) , pointer , dimension(:,:) :: lndcat
    real(dp) , pointer , dimension(:,:) :: msfx
    real(dp) , pointer , dimension(:,:) :: msfd
    real(dp) , pointer , dimension(:,:) :: coriol
  end type domain_io

  type (domain_io) :: mddom_io

  public :: domain_io , mddom_io , read_domain

  contains

  integer function read_domain(ncid)
    implicit none
    integer , intent(in) :: ncid

    read_domain = check_domain(ncid)
    if ( read_domain /= 0 ) then
      return
    end if
    call allocate_domain( )
    call read_var1d_static(ncid,'sigma',mddom_io%sigma)
    call read_var2d_static(ncid,'xlat',mddom_io%xlat)
    call read_var2d_static(ncid,'xlon',mddom_io%xlon)
    call read_var2d_static(ncid,'dlat',mddom_io%dlat)
    call read_var2d_static(ncid,'dlon',mddom_io%dlon)
    call read_var2d_static(ncid,'topo',mddom_io%ht)
    call read_var2d_static(ncid,'mask',mddom_io%mask)
    call read_var2d_static(ncid,'landuse',mddom_io%lndcat)
    call read_var2d_static(ncid,'xmap',mddom_io%msfx)
    call read_var2d_static(ncid,'dmap',mddom_io%msfd)
    call read_var2d_static(ncid,'coriol',mddom_io%coriol)
  end function read_domain

  subroutine allocate_domain
    implicit none
    call getmem1d(mddom_io%sigma,1,kz+1,'domain:sigma')
    call getmem2d(mddom_io%xlat,1,jx,1,iy,'domain:xlat')
    call getmem2d(mddom_io%xlon,1,jx,1,iy,'domain:xlon')
    call getmem2d(mddom_io%dlat,1,jx,1,iy,'domain:dlat')
    call getmem2d(mddom_io%dlon,1,jx,1,iy,'domain:dlon')
    call getmem2d(mddom_io%ht,1,jx,1,iy,'domain:ht')
    call getmem2d(mddom_io%mask,1,jx,1,iy,'domain:mask')
    call getmem2d(mddom_io%lndcat,1,jx,1,iy,'domain:lndcat')
    call getmem2d(mddom_io%msfx,1,jx,1,iy,'domain:msfx')
    call getmem2d(mddom_io%msfd,1,jx,1,iy,'domain:msfd')
    call getmem2d(mddom_io%coriol,1,jx,1,iy,'domain:coriol')
  end subroutine allocate_domain

  integer function check_domain(ncid)
    implicit none
    integer , intent(in) :: ncid
    integer :: istatus
    integer :: idimid , ivarid
    integer :: iyy , jxx , kzz
    character(6) :: proj
    real(sp) :: dsx , iclat , iclon , ptsp
    check_domain = -1
    istatus = nf90_inq_dimid(ncid, 'jx', idimid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error search dimension JX')
    istatus = nf90_inquire_dimension(ncid, idimid, len=jxx)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read dimension JX')
    if ( jx /= jxx ) then
      write(stderr,*) 'Mismatch: JX in DOMAIN file /= JX in namelist'
      write(stderr,*) 'DOMAIN FILE : ', jxx
      write(stderr,*) 'NAMELIST    : ', jx
      return
    end if
    istatus = nf90_inq_dimid(ncid, 'iy', idimid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error search dimension IY')
    istatus = nf90_inquire_dimension(ncid, idimid, len=iyy)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read dimension IY')
    if ( iy /= iyy ) then
      write(stderr,*) 'Mismatch: IY in DOMAIN file /= IY in namelist'
      write(stderr,*) 'DOMAIN FILE : ', iyy
      write(stderr,*) 'NAMELIST    : ', iy
      return
    end if
    istatus = nf90_inq_dimid(ncid, 'kz', idimid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error search dimension KZ')
    istatus = nf90_inquire_dimension(ncid, idimid, len=kzz)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read dimension KZ')
    if ( kz+1 /= kzz ) then
      write(stderr,*) 'Mismatch: KZ in DOMAIN file /= KZ+1 in namelist'
      write(stderr,*) 'DOMAIN FILE : ', kzz
      write(stderr,*) 'NAMELIST    : ', kz
      return
    end if
    istatus = nf90_inq_varid(ncid, 'ptop', ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error search variable PTOP')
    istatus = nf90_get_var(ncid, ivarid, ptsp)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read variable ptop')
    if ( dabs(dble(ptsp*d_r10)-dble(ptop)) > 0.001D+00 ) then
      write(stderr,*) 'Mismatch: PTOP in DOMAIN file /= PTOP in namelist'
      write(stderr,*) 'DOMAIN FILE : ', ptop
      write(stderr,*) 'NAMELIST    : ', ptsp
      return
    end if
    istatus = nf90_get_att(ncid, nf90_global, 'projection', proj)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read attribute projection')
    if (proj /= iproj) then
      write(stderr,*) 'Mismatch: IPROJ in DOMAIN file /= IPROJ in namelist'
      write(stderr,*) 'DOMAIN FILE : ', proj
      write(stderr,*) 'NAMELIST    : ', iproj
      return
    end if
    istatus = nf90_get_att(ncid, nf90_global,'grid_size_in_meters', dsx)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read attribute grid_size_in_meters')
    if (dabs(dble(dsx*d_r1000)-dble(ds)) > 0.001D+00 ) then 
      write(stderr,*) 'Mismatch: DS in DOMAIN file /= DS in namelist'
      write(stderr,*) 'DOMAIN FILE : ', dsx/1000.0
      write(stderr,*) 'NAMELIST    : ', ds
      return
    end if
    istatus = nf90_get_att(ncid, nf90_global, &
                           'latitude_of_projection_origin', iclat)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read attribute latitude_of_projection_origin')
    if (dabs(dble(iclat)-dble(clat)) > 0.001D+00 ) then
      write(stderr,*) 'Mismatch: CLAT in DOMAIN file /= CLAT in namelist'
      write(stderr,*) 'DOMAIN FILE : ', iclat
      write(stderr,*) 'NAMELIST    : ', clat
      return
    end if
    istatus = nf90_get_att(ncid, nf90_global, &
                           'longitude_of_projection_origin', iclon)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read attribute longitude_of_projection_origin')
    if (dabs(dble(iclon)-dble(clon)) > 0.001D+00 ) then
      write(stderr,*) 'Mismatch: CLON in DOMAIN file /= CLON in namelist'
      write(stderr,*) 'DOMAIN FILE : ', iclon
      write(stderr,*) 'NAMELIST    : ', clon
      return
    end if
    check_domain = 0
    return
  end function check_domain

end module mod_domain

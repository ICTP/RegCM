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
    real(rk8) , pointer , dimension(:) :: sigma
    real(rk8) , pointer , dimension(:,:) :: xlat
    real(rk8) , pointer , dimension(:,:) :: xlon
    real(rk8) , pointer , dimension(:,:) :: dlat
    real(rk8) , pointer , dimension(:,:) :: dlon
    real(rk8) , pointer , dimension(:,:) :: ht
    real(rk8) , pointer , dimension(:,:) :: mask
    real(rk8) , pointer , dimension(:,:) :: lndcat
    real(rk8) , pointer , dimension(:,:) :: msfx
    real(rk8) , pointer , dimension(:,:) :: msfd
    real(rk8) , pointer , dimension(:,:) :: coriol
  end type domain_io

  type (domain_io) :: mddom_io

  interface read_domain
    module procedure read_domain_type
    module procedure read_domain_array
    module procedure read_domain_array_single
  end interface read_domain

  public :: domain_io , mddom_io , read_domain , check_domain

  contains

  subroutine read_domain_type(ncid)
    implicit none
    integer(ik4) , intent(in) :: ncid
    call check_domain(ncid)
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
  end subroutine read_domain_type

  subroutine read_domain_array(ncid,sigma,xlat,xlon,dlat,dlon,ht,mask, &
                               lndcat,msfx,msfd,coriol)
    implicit none
    integer(ik4) , intent(in) :: ncid
    real(rk8) , pointer , dimension(:) , intent(out) :: sigma
    real(rk8) , pointer , dimension(:,:) , intent(out) , optional :: xlat
    real(rk8) , pointer , dimension(:,:) , intent(out) , optional :: xlon
    real(rk8) , pointer , dimension(:,:) , intent(out) , optional :: dlat
    real(rk8) , pointer , dimension(:,:) , intent(out) , optional :: dlon
    real(rk8) , pointer , dimension(:,:) , intent(out) , optional :: ht
    real(rk8) , pointer , dimension(:,:) , intent(out) , optional :: mask
    real(rk8) , pointer , dimension(:,:) , intent(out) , optional :: lndcat
    real(rk8) , pointer , dimension(:,:) , intent(out) , optional :: msfx
    real(rk8) , pointer , dimension(:,:) , intent(out) , optional :: msfd
    real(rk8) , pointer , dimension(:,:) , intent(out) , optional :: coriol
    call check_domain(ncid)
    call read_var1d_static(ncid,'sigma',sigma)
    if ( present(xlat) ) call read_var2d_static(ncid,'xlat',xlat)
    if ( present(xlon) ) call read_var2d_static(ncid,'xlon',xlon)
    if ( present(dlat) ) call read_var2d_static(ncid,'dlat',dlat)
    if ( present(dlon) ) call read_var2d_static(ncid,'dlon',dlon)
    if ( present(ht) ) call read_var2d_static(ncid,'topo',ht)
    if ( present(mask) ) call read_var2d_static(ncid,'mask',mask)
    if ( present(lndcat) ) call read_var2d_static(ncid,'landuse',lndcat)
    if ( present(msfx) ) call read_var2d_static(ncid,'xmap',msfx)
    if ( present(msfd) ) call read_var2d_static(ncid,'dmap',msfd)
    if ( present(coriol) ) call read_var2d_static(ncid,'coriol',coriol)
  end subroutine read_domain_array

  subroutine read_domain_array_single(ncid,sigma,xlat,xlon,dlat,dlon,ht,mask, &
                                      lndcat,msfx,msfd,coriol,ltrans)
    implicit none
    integer(ik4) , intent(in) :: ncid
    real(rk4) , pointer , dimension(:) , intent(out) :: sigma
    real(rk4) , pointer , dimension(:,:) , intent(out) , optional :: xlat
    real(rk4) , pointer , dimension(:,:) , intent(out) , optional :: xlon
    real(rk4) , pointer , dimension(:,:) , intent(out) , optional :: dlat
    real(rk4) , pointer , dimension(:,:) , intent(out) , optional :: dlon
    real(rk4) , pointer , dimension(:,:) , intent(out) , optional :: ht
    real(rk4) , pointer , dimension(:,:) , intent(out) , optional :: mask
    real(rk4) , pointer , dimension(:,:) , intent(out) , optional :: lndcat
    real(rk4) , pointer , dimension(:,:) , intent(out) , optional :: msfx
    real(rk4) , pointer , dimension(:,:) , intent(out) , optional :: msfd
    real(rk4) , pointer , dimension(:,:) , intent(out) , optional :: coriol
    logical , intent(in) , optional :: ltrans
    real(rk4) , pointer , dimension(:,:) :: dum
    logical :: dotrans
    call check_domain(ncid)
    call read_var1d_static(ncid,'sigma',sigma)

    dotrans = .false.
    if ( present(ltrans) ) then
      dotrans = ltrans
    end if
    if ( dotrans ) then
      call getmem2d(dum,1,jx,1,iy,'read_domain_arrays_single:dum')
      if ( present(xlat) ) then
        call read_var2d_static(ncid,'xlat',dum)
        xlat = transpose(dum)
      end if
      if ( present(xlon) ) then
        call read_var2d_static(ncid,'xlon',dum)
        xlon = transpose(dum)
      end if
      if ( present(dlat) ) then
        call read_var2d_static(ncid,'dlat',dum)
        dlat = transpose(dum)
      end if
      if ( present(dlon) ) then
        call read_var2d_static(ncid,'dlon',dum)
        dlon = transpose(dum)
      end if
      if ( present(ht) ) then
        call read_var2d_static(ncid,'topo',dum)
        ht = transpose(dum)
      end if
      if ( present(mask) ) then
        call read_var2d_static(ncid,'mask',dum)
        mask = transpose(dum)
      end if
      if ( present(lndcat) ) then
        call read_var2d_static(ncid,'landuse',dum)
        lndcat = transpose(dum)
      end if
      if ( present(msfx) ) then
        call read_var2d_static(ncid,'xmap',dum)
        msfx = transpose(dum)
      end if
      if ( present(msfd) ) then
        call read_var2d_static(ncid,'dmap',dum)
        msfd = transpose(dum)
      end if
      if ( present(coriol) ) then
        call read_var2d_static(ncid,'coriol',dum)
        coriol = transpose(dum)
      end if
      call relmem2d(dum)
    else
      if ( present(xlat) ) call read_var2d_static(ncid,'xlat',xlat)
      if ( present(xlon) ) call read_var2d_static(ncid,'xlon',xlon)
      if ( present(dlat) ) call read_var2d_static(ncid,'dlat',dlat)
      if ( present(dlon) ) call read_var2d_static(ncid,'dlon',dlon)
      if ( present(ht) ) call read_var2d_static(ncid,'topo',ht)
      if ( present(mask) ) call read_var2d_static(ncid,'mask',mask)
      if ( present(lndcat) ) call read_var2d_static(ncid,'landuse',lndcat)
      if ( present(msfx) ) call read_var2d_static(ncid,'xmap',msfx)
      if ( present(msfd) ) call read_var2d_static(ncid,'dmap',msfd)
      if ( present(coriol) ) call read_var2d_static(ncid,'coriol',coriol)
    end if
  end subroutine read_domain_array_single

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

  subroutine check_domain(ncid,lmod)
    implicit none
    integer(ik4) , intent(in) :: ncid
    logical , optional :: lmod
    integer(ik4) :: istatus
    integer(ik4) :: idimid , ivarid
    integer(ik4) :: iyy , jxx , kzz , kcheck
    character(6) :: proj
    logical :: lh
    real(rk4) :: dsx , iclat , iclon , ptsp

    lh = .false.
    if ( present(lmod) ) lh = lmod
    istatus = nf90_inq_dimid(ncid, 'jx', idimid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error search dimension JX')
    istatus = nf90_inquire_dimension(ncid, idimid, len=jxx)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read dimension JX')
    if ( jx /= jxx ) then
      write(stderr,*) 'DOMAIN FILE : ', jxx
      write(stderr,*) 'NAMELIST    : ', jx
      call die('Mismatch: JX in DOMAIN file /= JX in namelist')
    end if
    istatus = nf90_inq_dimid(ncid, 'iy', idimid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error search dimension IY')
    istatus = nf90_inquire_dimension(ncid, idimid, len=iyy)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read dimension IY')
    if ( iy /= iyy ) then
      write(stderr,*) 'DOMAIN FILE : ', iyy
      write(stderr,*) 'NAMELIST    : ', iy
      call die('Mismatch: IY in DOMAIN file /= IY in namelist')
    end if
    istatus = nf90_inq_dimid(ncid, 'kz', idimid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error search dimension KZ')
    istatus = nf90_inquire_dimension(ncid, idimid, len=kzz)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read dimension KZ')
    kcheck = kzp1
    if ( lh ) kcheck = kz
    if ( kcheck /= kzz ) then
      write(stderr,*) 'DOMAIN FILE : ', kzz
      write(stderr,*) 'NAMELIST    : ', kz
      if ( lh ) then
        call die('Mismatch: KZ in DOMAIN file /= KZ in namelist')
      else
        call die('Mismatch: KZ in DOMAIN file /= KZ+1 in namelist')
      end if 
    end if
    istatus = nf90_inq_varid(ncid, 'ptop', ivarid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error search variable PTOP')
    istatus = nf90_get_var(ncid, ivarid, ptsp)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read variable ptop')
    if ( dabs(dble(ptsp*d_r10)-dble(ptop)) > 0.001D+00 ) then
      write(stderr,*) 'DOMAIN FILE : ', ptop
      write(stderr,*) 'NAMELIST    : ', ptsp
      call die('Mismatch: PTOP in DOMAIN file /= PTOP in namelist')
    end if
    istatus = nf90_get_att(ncid, nf90_global, 'projection', proj)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read attribute projection')
    if (proj /= iproj) then
      write(stderr,*) 'Mismatch: IPROJ in DOMAIN file /= IPROJ in namelist'
      write(stderr,*) 'DOMAIN FILE : ', proj
      write(stderr,*) 'NAMELIST    : ', iproj
      call die('Mismatch: IPROJ in DOMAIN file /= IPROJ in namelist')
    end if
    istatus = nf90_get_att(ncid, nf90_global,'grid_size_in_meters', dsx)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read attribute grid_size_in_meters')
    if (dabs(dble(dsx*d_r1000)-dble(ds)) > 0.001D+00 ) then 
      write(stderr,*) 'DOMAIN FILE : ', dsx/1000.0
      write(stderr,*) 'NAMELIST    : ', ds
      call die('Mismatch: DS in DOMAIN file /= DS in namelist')
    end if
    istatus = nf90_get_att(ncid, nf90_global, &
                           'latitude_of_projection_origin', iclat)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read attribute latitude_of_projection_origin')
    if (dabs(dble(iclat)-dble(clat)) > 0.001D+00 ) then
      write(stderr,*) 'DOMAIN FILE : ', iclat
      write(stderr,*) 'NAMELIST    : ', clat
      call die('Mismatch: CLAT in DOMAIN file /= CLAT in namelist')
    end if
    istatus = nf90_get_att(ncid, nf90_global, &
                           'longitude_of_projection_origin', iclon)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read attribute longitude_of_projection_origin')
    if (dabs(dble(iclon)-dble(clon)) > 0.001D+00 ) then
      write(stderr,*) 'DOMAIN FILE : ', iclon
      write(stderr,*) 'NAMELIST    : ', clon
      call die('Mismatch: CLON in DOMAIN file /= CLON in namelist')
    end if
    istatus = nf90_get_att(ncid,nf90_global,'grid_factor',xcone)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read attribute grid_factor')
  end subroutine check_domain

end module mod_domain

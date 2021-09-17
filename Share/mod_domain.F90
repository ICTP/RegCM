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
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_memutil
  use mod_constants
  use mod_nchelper
  use mod_stdio
  use mod_message

  implicit none

  private

  type domain_io
    real(rkx) , pointer , dimension(:) :: sigma
    real(rkx) , pointer , dimension(:,:) :: xlat
    real(rkx) , pointer , dimension(:,:) :: xlon
    real(rkx) , pointer , dimension(:,:) :: dlat
    real(rkx) , pointer , dimension(:,:) :: dlon
    real(rkx) , pointer , dimension(:,:) :: ulat
    real(rkx) , pointer , dimension(:,:) :: ulon
    real(rkx) , pointer , dimension(:,:) :: vlat
    real(rkx) , pointer , dimension(:,:) :: vlon
    real(rkx) , pointer , dimension(:,:) :: ht
    real(rkx) , pointer , dimension(:,:) :: mask
    real(rkx) , pointer , dimension(:,:) :: lndcat
    real(rkx) , pointer , dimension(:,:) :: msfx
    real(rkx) , pointer , dimension(:,:) :: msfd
    real(rkx) , pointer , dimension(:,:) :: coriol
    real(rkx) , pointer , dimension(:,:) :: snowam
    real(rkx) , pointer , dimension(:,:) :: hlake
  end type domain_io

  type (domain_io) :: mddom_io

  interface read_domain
    module procedure read_domain_type
    module procedure read_domain_array_double
    module procedure read_domain_array_single
  end interface read_domain

  public :: domain_io , mddom_io , read_domain , check_domain
  public :: read_reference_state , read_reference_surface_temp

  contains

  subroutine read_domain_type(ncid)
    implicit none
    integer(ik4) , intent(in) :: ncid
    logical :: has_snow , has_dhlake , has_kz
    has_snow = .true.
    has_dhlake = .true.
    has_kz = .true.
    call check_domain(ncid)
    call allocate_domain( )
    call read_var1d_static(ncid,'kz',mddom_io%sigma,has_kz)
    if ( .not. has_kz ) then
      call read_var1d_static(ncid,'sigma',mddom_io%sigma)
    end if
    call read_var2d_static(ncid,'xlat',mddom_io%xlat)
    call read_var2d_static(ncid,'xlon',mddom_io%xlon)
    if ( idynamic == 3 ) then
      call read_var2d_static(ncid,'ulat',mddom_io%ulat)
      call read_var2d_static(ncid,'ulon',mddom_io%ulon)
      call read_var2d_static(ncid,'vlat',mddom_io%vlat)
      call read_var2d_static(ncid,'vlon',mddom_io%vlon)
    else
      call read_var2d_static(ncid,'dlat',mddom_io%dlat)
      call read_var2d_static(ncid,'dlon',mddom_io%dlon)
      call read_var2d_static(ncid,'xmap',mddom_io%msfx)
      call read_var2d_static(ncid,'dmap',mddom_io%msfd)
    end if
    call read_var2d_static(ncid,'topo',mddom_io%ht)
    call read_var2d_static(ncid,'mask',mddom_io%mask)
    call read_var2d_static(ncid,'landuse',mddom_io%lndcat)
    call read_var2d_static(ncid,'coriol',mddom_io%coriol)
    call read_var2d_static(ncid,'snowam',mddom_io%snowam,has_snow)
    call read_var2d_static(ncid,'dhlake',mddom_io%hlake,has_dhlake)
  end subroutine read_domain_type

  subroutine read_domain_array_double(ncid,sigma,xlat,xlon,dlat,dlon, &
                                      ulat,ulon,vlat,vlon,ht,mask,    &
                                      lndcat,msfx,msfd,coriol,snowam, &
                                      hlake,lsubgrid)
    implicit none
    integer(ik4) , intent(in) :: ncid
    real(rk8) , pointer , dimension(:) , intent(inout) :: sigma
    real(rk8) , pointer , dimension(:,:) , intent(inout) , optional :: xlat
    real(rk8) , pointer , dimension(:,:) , intent(inout) , optional :: xlon
    real(rk8) , pointer , dimension(:,:) , intent(inout) , optional :: dlat
    real(rk8) , pointer , dimension(:,:) , intent(inout) , optional :: dlon
    real(rk8) , pointer , dimension(:,:) , intent(inout) , optional :: ulat
    real(rk8) , pointer , dimension(:,:) , intent(inout) , optional :: ulon
    real(rk8) , pointer , dimension(:,:) , intent(inout) , optional :: vlat
    real(rk8) , pointer , dimension(:,:) , intent(inout) , optional :: vlon
    real(rk8) , pointer , dimension(:,:) , intent(inout) , optional :: ht
    real(rk8) , pointer , dimension(:,:) , intent(inout) , optional :: mask
    real(rk8) , pointer , dimension(:,:) , intent(inout) , optional :: lndcat
    real(rk8) , pointer , dimension(:,:) , intent(inout) , optional :: msfx
    real(rk8) , pointer , dimension(:,:) , intent(inout) , optional :: msfd
    real(rk8) , pointer , dimension(:,:) , intent(inout) , optional :: coriol
    real(rk8) , pointer , dimension(:,:) , intent(inout) , optional :: snowam
    real(rk8) , pointer , dimension(:,:) , intent(inout) , optional :: hlake
    logical , intent(in) , optional :: lsubgrid
    logical :: has_snow , has_dhlake , has_kz
    has_snow = .true.
    has_dhlake = .true.
    has_kz = .true.
    if ( present(lsubgrid) ) then
      call check_domain(ncid,lsubgrid=lsubgrid)
    else
      call check_domain(ncid)
    end if
    call read_var1d_static(ncid,'kz',sigma,has_kz)
    if ( .not. has_kz ) then
      call read_var1d_static(ncid,'sigma',sigma)
    end if
    if ( present(xlat) ) call read_var2d_static(ncid,'xlat',xlat)
    if ( present(xlon) ) call read_var2d_static(ncid,'xlon',xlon)
    if ( idynamic == 3 ) then
      if ( present(ulat) ) call read_var2d_static(ncid,'ulat',ulat)
      if ( present(ulon) ) call read_var2d_static(ncid,'ulon',ulon)
      if ( present(vlat) ) call read_var2d_static(ncid,'vlat',vlat)
      if ( present(vlon) ) call read_var2d_static(ncid,'vlon',vlon)
    else
      if ( present(dlat) ) call read_var2d_static(ncid,'dlat',dlat)
      if ( present(dlon) ) call read_var2d_static(ncid,'dlon',dlon)
      if ( present(msfx) ) call read_var2d_static(ncid,'xmap',msfx)
      if ( present(msfd) ) call read_var2d_static(ncid,'dmap',msfd)
    end if
    if ( present(ht) ) call read_var2d_static(ncid,'topo',ht)
    if ( present(mask) ) call read_var2d_static(ncid,'mask',mask)
    if ( present(lndcat) ) call read_var2d_static(ncid,'landuse',lndcat)
    if ( present(coriol) ) call read_var2d_static(ncid,'coriol',coriol)
    if ( present(snowam) ) call read_var2d_static(ncid,'snowam',snowam,has_snow)
    if ( present(hlake) ) call read_var2d_static(ncid,'dhlake',hlake,has_dhlake)
  end subroutine read_domain_array_double

  subroutine read_domain_array_single(ncid,sigma,xlat,xlon,dlat,dlon, &
                                      ulat,ulon,vlat,vlon,ht,mask,    &
                                      lndcat,msfx,msfd,coriol,snowam, &
                                      hlake,lsubgrid)
    implicit none
    integer(ik4) , intent(in) :: ncid
    real(rk4) , pointer , dimension(:) , intent(inout) :: sigma
    real(rk4) , pointer , dimension(:,:) , intent(inout) , optional :: xlat
    real(rk4) , pointer , dimension(:,:) , intent(inout) , optional :: xlon
    real(rk4) , pointer , dimension(:,:) , intent(inout) , optional :: dlat
    real(rk4) , pointer , dimension(:,:) , intent(inout) , optional :: dlon
    real(rk4) , pointer , dimension(:,:) , intent(inout) , optional :: ulat
    real(rk4) , pointer , dimension(:,:) , intent(inout) , optional :: ulon
    real(rk4) , pointer , dimension(:,:) , intent(inout) , optional :: vlat
    real(rk4) , pointer , dimension(:,:) , intent(inout) , optional :: vlon
    real(rk4) , pointer , dimension(:,:) , intent(inout) , optional :: ht
    real(rk4) , pointer , dimension(:,:) , intent(inout) , optional :: mask
    real(rk4) , pointer , dimension(:,:) , intent(inout) , optional :: lndcat
    real(rk4) , pointer , dimension(:,:) , intent(inout) , optional :: msfx
    real(rk4) , pointer , dimension(:,:) , intent(inout) , optional :: msfd
    real(rk4) , pointer , dimension(:,:) , intent(inout) , optional :: coriol
    real(rk4) , pointer , dimension(:,:) , intent(inout) , optional :: snowam
    real(rk4) , pointer , dimension(:,:) , intent(inout) , optional :: hlake
    logical , intent(in) , optional :: lsubgrid
    logical :: has_snow , has_dhlake , has_kz
    has_snow = .true.
    has_dhlake = .true.
    has_kz = .true.
    if ( present(lsubgrid) ) then
      call check_domain(ncid,lsubgrid=lsubgrid)
    else
      call check_domain(ncid)
    end if
    call read_var1d_static(ncid,'kz',sigma,has_kz)
    if ( .not. has_kz ) then
      call read_var1d_static(ncid,'sigma',sigma)
    end if
    if ( present(xlat) ) call read_var2d_static(ncid,'xlat',xlat)
    if ( present(xlon) ) call read_var2d_static(ncid,'xlon',xlon)
    if ( idynamic == 3 ) then
      if ( present(ulat) ) call read_var2d_static(ncid,'ulat',ulat)
      if ( present(ulon) ) call read_var2d_static(ncid,'ulon',ulon)
      if ( present(vlat) ) call read_var2d_static(ncid,'vlat',vlat)
      if ( present(vlon) ) call read_var2d_static(ncid,'vlon',vlon)
    else
      if ( present(dlat) ) call read_var2d_static(ncid,'dlat',dlat)
      if ( present(dlon) ) call read_var2d_static(ncid,'dlon',dlon)
      if ( present(msfx) ) call read_var2d_static(ncid,'xmap',msfx)
      if ( present(msfd) ) call read_var2d_static(ncid,'dmap',msfd)
    end if
    if ( present(ht) ) call read_var2d_static(ncid,'topo',ht)
    if ( present(mask) ) call read_var2d_static(ncid,'mask',mask)
    if ( present(lndcat) ) call read_var2d_static(ncid,'landuse',lndcat)
    if ( present(coriol) ) call read_var2d_static(ncid,'coriol',coriol)
    if ( present(snowam) ) call read_var2d_static(ncid,'snowam',snowam,has_snow)
    if ( present(hlake) ) call read_var2d_static(ncid,'dhlake',hlake,has_dhlake)
  end subroutine read_domain_array_single

  subroutine read_reference_state(ncid,ps0,pr0,t0,rho0,ts0)
    implicit none
    integer , intent(in) :: ncid
    real(rkx) , pointer , dimension(:,:) , intent(inout) , optional :: ps0
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) , optional :: pr0
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) , optional :: t0
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) , optional :: rho0
    real(rkx) , intent(out) :: ts0
    call read_var2d_static(ncid,'ps0',ps0)
    call read_var3d_static(ncid,'pr0',pr0)
    call read_var3d_static(ncid,'t0',t0)
    call read_var3d_static(ncid,'rho0',rho0)
    call get_attribute(ncid,'base_state_surface_temperature',ts0)
  end subroutine read_reference_state

  subroutine read_reference_surface_temp(ncid,ts0)
    implicit none
    integer , intent(in) :: ncid
    real(rkx) , intent(out) :: ts0
    call get_attribute(ncid,'base_state_surface_temperature',ts0)
  end subroutine read_reference_surface_temp

  subroutine allocate_domain
    implicit none
    call getmem1d(mddom_io%sigma,1,kz+1,'domain:sigma')
    call getmem2d(mddom_io%xlat,1,jx,1,iy,'domain:xlat')
    call getmem2d(mddom_io%xlon,1,jx,1,iy,'domain:xlon')
    if ( idynamic == 3 ) then
      call getmem2d(mddom_io%ulat,1,jx,1,iy,'domain:ulat')
      call getmem2d(mddom_io%ulon,1,jx,1,iy,'domain:ulon')
      call getmem2d(mddom_io%vlat,1,jx,1,iy,'domain:vlat')
      call getmem2d(mddom_io%vlon,1,jx,1,iy,'domain:vlon')
    else
      call getmem2d(mddom_io%dlat,1,jx,1,iy,'domain:dlat')
      call getmem2d(mddom_io%dlon,1,jx,1,iy,'domain:dlon')
      call getmem2d(mddom_io%msfx,1,jx,1,iy,'domain:msfx')
      call getmem2d(mddom_io%msfd,1,jx,1,iy,'domain:msfd')
    end if
    call getmem2d(mddom_io%ht,1,jx,1,iy,'domain:ht')
    call getmem2d(mddom_io%mask,1,jx,1,iy,'domain:mask')
    call getmem2d(mddom_io%lndcat,1,jx,1,iy,'domain:lndcat')
    call getmem2d(mddom_io%coriol,1,jx,1,iy,'domain:coriol')
    call getmem2d(mddom_io%snowam,1,jx,1,iy,'domain:snowam')
    call getmem2d(mddom_io%hlake,1,jx,1,iy,'domain:hlake')
  end subroutine allocate_domain

  subroutine check_domain(ncid,lmod,linternal,lsubgrid)
    implicit none
    integer(ik4) , intent(in) :: ncid
    logical , optional :: lmod , linternal , lsubgrid
    integer(ik4) :: istatus
    integer(ik4) :: idimid , ivarid
    integer(ik4) :: iyy , jxx , kzz , kcheck , jcheck , icheck
    character(len=6) :: proj
    logical :: lh , lb , ls
    real(rkx) :: dsx , iclat , iclon , ptsp
    real(rkx) , dimension(2) :: icntr

    lh = .false.
    lb = .false.
    ls = .false.
    if ( present(lmod) ) lh = lmod
    if ( present(linternal) ) lb = linternal
    if ( present(lsubgrid) ) ls = lsubgrid
    istatus = nf90_inq_dimid(ncid, 'jx', idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error search dimension JX')
    istatus = nf90_inquire_dimension(ncid, idimid, len=jxx)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read dimension JX')
    jcheck = jx
    if ( lb .and. i_band /= 1 ) jcheck = jcheck - 3
    if ( ls ) jcheck = jcheck * nsg
    if ( jcheck /= jxx ) then
      write(stderr,*) 'DOMAIN FILE : ', jxx
      write(stderr,*) 'NAMELIST    : ', jx
      if ( lb ) then
        call die('Mismatch: JX+3 in DOMAIN file /= JX in namelist')
      else
        call die('Mismatch: JX in DOMAIN file /= JX in namelist')
      end if
    end if
    istatus = nf90_inq_dimid(ncid, 'iy', idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error search dimension IY')
    istatus = nf90_inquire_dimension(ncid, idimid, len=iyy)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read dimension IY')
    icheck = iy
    if ( lb .and. i_crm /= 1 ) icheck = icheck - 3
    if ( ls ) icheck = icheck * nsg
    if ( icheck /= iyy ) then
      write(stderr,*) 'DOMAIN FILE : ', iyy
      write(stderr,*) 'NAMELIST    : ', iy
      if ( lb ) then
        call die('Mismatch: IY+3 in DOMAIN file /= IY in namelist')
      else
        call die('Mismatch: IY in DOMAIN file /= IY in namelist')
      end if
    end if
    istatus = nf90_inq_dimid(ncid, 'kz', idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error search dimension KZ')
    istatus = nf90_inquire_dimension(ncid, idimid, len=kzz)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read dimension KZ')
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
    if ( idynamic < 3 ) then
      istatus = nf90_inq_varid(ncid, 'ptop', ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error search variable PTOP')
      istatus = nf90_get_var(ncid, ivarid, ptsp)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read variable ptop')
      if ( abs(real(ptsp*d_r10,rkx)-real(ptop,rkx)) > 0.001_rkx ) then
        write(stderr,*) 'DOMAIN FILE : ', ptsp
        write(stderr,*) 'NAMELIST    : ', ptop
        call die('Mismatch: PTOP in DOMAIN file /= PTOP in namelist')
      end if
    else
      ptop = d_zero
    end if
    istatus = nf90_get_att(ncid, nf90_global, 'projection', proj)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read attribute projection')
    if (proj /= iproj) then
      write(stderr,*) 'Mismatch: IPROJ in DOMAIN file /= IPROJ in namelist'
      write(stderr,*) 'DOMAIN FILE : ', proj
      write(stderr,*) 'NAMELIST    : ', iproj
      call die('Mismatch: IPROJ in DOMAIN file /= IPROJ in namelist')
    end if
    istatus = nf90_get_att(ncid, nf90_global,'grid_size_in_meters', dsx)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read attribute grid_size_in_meters')
    if ( ls ) dsx = dsx * real(nsg,rkx)
    if (abs(real(dsx*d_r1000,rkx)-real(ds,rkx)) > 0.001_rkx ) then
      write(stderr,*) 'DOMAIN FILE : ', dsx/1000.0
      write(stderr,*) 'NAMELIST    : ', ds
      call die('Mismatch: DS in DOMAIN file /= DS in namelist')
    end if
    istatus = nf90_get_att(ncid, nf90_global, &
                           'latitude_of_projection_origin', iclat)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read attribute latitude_of_projection_origin')
    if (abs(real(iclat,rkx)-real(clat,rkx)) > 0.001_rkx ) then
      write(stderr,*) 'DOMAIN FILE : ', iclat
      write(stderr,*) 'NAMELIST    : ', clat
      call die('Mismatch: CLAT in DOMAIN file /= CLAT in namelist')
    end if
    istatus = nf90_get_att(ncid, nf90_global, &
                           'longitude_of_projection_origin', iclon)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read attribute longitude_of_projection_origin')
    if (abs(real(iclon,rkx)-real(clon,rkx)) > 0.001_rkx ) then
      write(stderr,*) 'DOMAIN FILE : ', iclon
      write(stderr,*) 'NAMELIST    : ', clon
      call die('Mismatch: CLON in DOMAIN file /= CLON in namelist')
    end if
    istatus = nf90_get_att(ncid, nf90_global, &
                           'index_of_projection_origin', icntr)
    if ( istatus /= nf90_noerr ) then
      write(stderr,*) 'WARNING: USING DOMAIN FILE FROM PREVIOUS MODEL VERSION'
      write(stderr,*) 'WARNING: ASSUMING PROJECTION CENTER IS DOMAIN CENTER'
    else
      if (abs(real(icntr(2),rkx)-real(cntri,rkx)) > 0.001_rkx ) then
        write(stderr,*) 'DOMAIN FILE : ', icntr(2)
        write(stderr,*) 'NAMELIST    : ', cntri
        call die('Mismatch: CNTRI in DOMAIN file /= CNTRI in namelist')
      end if
      if (abs(real(icntr(1),rkx)-real(cntrj,rkx)) > 0.001_rkx ) then
        write(stderr,*) 'DOMAIN FILE : ', icntr(1)
        write(stderr,*) 'NAMELIST    : ', cntrj
        call die('Mismatch: CNTRJ in DOMAIN file /= CNTRJ in namelist')
      end if
    end if
    istatus = nf90_get_att(ncid,nf90_global,'grid_factor',xcone)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read attribute grid_factor')
  end subroutine check_domain

end module mod_domain
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

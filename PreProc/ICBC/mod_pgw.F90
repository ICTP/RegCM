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

module mod_pgw

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_constants
  use mod_stdio
  use mod_memutil
  use mod_grid
  use mod_write
  use mod_vertint
  use mod_earth
  use mod_hgt
  use mod_humid
  use mod_projections
  use mod_vectutil
  use mod_message
  use mod_nchelper
  use mod_kdinterp
  use netcdf

  private

  public :: init_pgw , get_pgw , conclude_pgw

  integer :: ncid = -1
  integer :: ilon , jlat , klev

  real(rkx) , pointer , dimension(:) :: glat
  real(rkx) , pointer , dimension(:) :: glon
  real(rkx) , pointer , public , dimension(:) :: plevs

  real(rkx) , pointer , dimension(:,:,:) :: u2 , v2 , q2 , z2 , t2
  real(rkx) , pointer , dimension(:,:) :: ps2 , ts2

  real(rkx) , pointer , dimension(:,:,:) :: d3u , d3v , h3u , h3v

  type(h_interpolator) :: cross_hint , udot_hint , vdot_hint

  character(len=256) :: pgwfile

  integer , parameter :: ndelta = 7
  character(len=10) , dimension(ndelta) :: varname
  integer , dimension(ndelta) :: ivardelta
  data varname /'delta_ta ' , 'delta_zg ' , &
                'delta_hus' , 'delta_ua ' , &
                'delta_va ' , 'delta_ps ' , &
                'delta_ts'/

  contains

  subroutine init_pgw(filename)
    implicit none
    character(len=*) , intent(in) :: filename
    integer :: istatus , idimid , ivarid , kkrec

    pgwfile = filename
    istatus = nf90_open(pgwfile,nf90_nowrite,ncid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error open file '//trim(pgwfile))
    istatus = nf90_inq_dimid(ncid,'lat',idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Missing lat dimension in file '//trim(pgwfile))
    istatus = nf90_inquire_dimension(ncid,idimid,len=jlat)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Error reading lat dimelen in file '//trim(pgwfile))
    istatus = nf90_inq_dimid(ncid,'lon',idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Missing lon dimension in file '//trim(pgwfile))
    istatus = nf90_inquire_dimension(ncid,idimid,len=ilon)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Error reading lon dimelen in file '//trim(pgwfile))
    istatus = nf90_inq_dimid(ncid,'plev',idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
            'Missing plev dimension in file '//trim(pgwfile))
    istatus = nf90_inquire_dimension(ncid,idimid,len=klev)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Error reading plev dimelen in file '//trim(pgwfile))

    call getmem1d(plevs,1,klev,'pgw:plevs')
    call getmem1d(glat,1,jlat,'pgw:glat')
    call getmem1d(glon,1,ilon,'pgw:glon')

    istatus = nf90_inq_varid(ncid,'lat',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Missing lat variable in file '//trim(pgwfile))
    istatus = nf90_get_var(ncid,ivarid,glat)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Error reading lat variable in file '//trim(pgwfile))
    istatus = nf90_inq_varid(ncid,'lon',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Missing lon variable in file '//trim(pgwfile))
    istatus = nf90_get_var(ncid,ivarid,glon)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Error reading lon variable in file '//trim(pgwfile))
    istatus = nf90_inq_varid(ncid,'plev',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
            'Missing plev variable in file '//trim(pgwfile))
    istatus = nf90_get_var(ncid,ivarid,plevs)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Error reading plev variable in file '//trim(pgwfile))

    call h_interpolator_create(cross_hint,glat,glon,xlat,xlon)
    if ( idynamic == 3 ) then
      call h_interpolator_create(udot_hint,glat,glon,ulat,ulon)
      call h_interpolator_create(vdot_hint,glat,glon,vlat,vlon)
    else
      call h_interpolator_create(udot_hint,glat,glon,dlat,dlon)
    end if

    call getmem3d(u2,1,ilon,1,jlat,1,klev,'pgw:u2')
    call getmem3d(v2,1,ilon,1,jlat,1,klev,'pgw:v2')
    call getmem3d(t2,1,ilon,1,jlat,1,klev,'pgw:t2')
    call getmem3d(q2,1,ilon,1,jlat,1,klev,'pgw:q2')
    call getmem3d(z2,1,ilon,1,jlat,1,klev,'pgw:z2')
    call getmem2d(ps2,1,ilon,1,jlat,'pgw:ps2')
    call getmem2d(ts2,1,ilon,1,jlat,'pgw:ts2')

    call getmem3d(h3u,1,jx,1,iy,1,klev,'pgw:h3u')
    call getmem3d(h3v,1,jx,1,iy,1,klev,'pgw:h3v')
    if ( idynamic == 3 ) then
      call getmem3d(d3u,1,jx,1,iy,1,klev,'pgw:d3u')
      call getmem3d(d3v,1,jx,1,iy,1,klev,'pgw:d3v')
    end if

    do kkrec = 1 , ndelta
      istatus = nf90_inq_varid(ncid,varname(kkrec),ivardelta(kkrec))
      call checkncerr(istatus,__FILE__,__LINE__, &
        'Error find var '//varname(kkrec))
    end do
  end subroutine init_pgw

  subroutine get_pgw(month)
    implicit none
    integer(ik4) , intent(in) :: month
    integer(ik4) :: istatus , it
    integer(ik4) , dimension(4) :: icount , istart

    it = month
    istart(1) = 1
    icount(1) = ilon
    istart(2) = 1
    icount(2) = jlat
    istart(3) = 1
    icount(3) = klev
    istart(4) = it
    icount(4) = 1
    istatus = nf90_get_var(ncid,ivardelta(1),t2,istart,icount)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var '//varname(1))
    istatus = nf90_get_var(ncid,ivardelta(2),z2,istart,icount)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var '//varname(2))
    istatus = nf90_get_var(ncid,ivardelta(3),q2,istart,icount)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var '//varname(3))
    istatus = nf90_get_var(ncid,ivardelta(4),u2,istart,icount)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var '//varname(4))
    istatus = nf90_get_var(ncid,ivardelta(5),v2,istart,icount)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var '//varname(5))
    istart(3) = 1
    icount(3) = 1
    istart(4) = it
    icount(4) = 1
    istatus = nf90_get_var(ncid,ivardelta(6),ps2,istart,icount)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var '//varname(6))
    istatus = nf90_get_var(ncid,ivardelta(7),ts2,istart,icount)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var '//varname(7))
    write (stdout,*) 'READ IN fields for month:' , month
    call h_interpolate_cont(cross_hint,t2,t4)
    call h_interpolate_cont(cross_hint,q2,q4)
    call h_interpolate_cont(cross_hint,z2,z4)
    if ( idynamic == 3 ) then
      call h_interpolate_cont(udot_hint,u2,u4)
      call h_interpolate_cont(udot_hint,v2,d3v)
      call h_interpolate_cont(vdot_hint,u2,d3u)
      call h_interpolate_cont(vdot_hint,v2,v4)
      call pju%wind_rotate(u4,d3v)
      call pjv%wind_rotate(d3u,v4)
    else
      call h_interpolate_cont(udot_hint,u2,u4)
      call h_interpolate_cont(udot_hint,v2,v4)
      call pjd%wind_rotate(u4,v4)
    end if
    call h_interpolate_cont(cross_hint,ts2,ts4)
    call h_interpolate_cont(cross_hint,ps2,ps4)
    q4 = q4/(d_one-q4)
  end subroutine get_pgw

  subroutine conclude_pgw
    implicit none
    integer :: istatus
    if ( ncid >= 0 ) then
      istatus = nf90_close(ncid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error close file '//pgwfile)
    end if
    call h_interpolator_destroy(cross_hint)
    call h_interpolator_destroy(udot_hint)
    if ( idynamic == 3 ) then
      call h_interpolator_destroy(vdot_hint)
    end if
  end subroutine conclude_pgw

end module mod_pgw
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

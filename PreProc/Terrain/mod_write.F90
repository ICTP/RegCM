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

module mod_write

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_nchelper
  use netcdf

  private

  public :: write_domain

  contains
!
  subroutine write_domain(fname,lsub,lndfudge,texfudge,lakfudge,ntype,sigma, &
                          xlat,xlon,dlat,dlon,xmap,dmap,coriol,mask,    &
                          htgrid,lndout,snowam,dpth,texout,frac_tex)
    implicit none
    character (len=*) , intent(in) :: fname
    logical , intent(in) :: lsub , lndfudge , texfudge , lakfudge
    integer(ik4) , intent(in) :: ntype
    real(rk4) , dimension(:) , pointer , intent(in) :: sigma
    real(rk4) , dimension(:,:) , pointer , intent(in) :: xlat , xlon
    real(rk4) , dimension(:,:) , pointer , intent(in) :: dlat , dlon
    real(rk4) , dimension(:,:) , pointer , intent(in) :: xmap , dmap , coriol
    real(rk4) , dimension(:,:) , pointer , intent(in) :: mask
    real(rk4) , dimension(:,:) , pointer , intent(in) :: htgrid , lndout
    real(rk4) , dimension(:,:) , pointer , intent(in) :: snowam
    real(rk4) , dimension(:,:) , pointer , intent(in) :: dpth
    real(rk4) , dimension(:,:) , pointer , intent(in) :: texout
    real(rk4) , dimension(:,:,:) , pointer , intent(in) :: frac_tex

    integer(ik4) :: ncid , incstat
    integer(ik4) :: nx , ny , nz , ipnt
    integer(ik4) , dimension(4) :: idims
    integer(ik4) , dimension(2) :: ihvar
    integer(ik4) , dimension(2) :: izvar
    integer(ik4) , dimension(16) :: ivar
    character(len=128) :: cdum
    real(rk4) :: hptop
    real(rk4) , pointer , dimension(:) :: yiy
    real(rk4) , pointer , dimension(:) :: xjx

    hptop = real(ptop * 10.0D0)

    ! here XLAT is transposed !!! (but why not change ?)
    nx = ubound(xlat,2)
    ny = ubound(xlat,1)
    nz = ubound(sigma,1)

    call createfile_withname(fname,ncid)
    call add_common_global_params(ncid,'terrain',lsub)

    ! Terrain related global parameters
    call cdumlogical(cdum,smthbdy)
    incstat = nf90_put_att(ncid, nf90_global, 'boundary_smoothing', cdum)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding global boundary_smoothing')
    incstat = nf90_put_att(ncid, nf90_global,                     &
                 'minimum_h2o_pct_for_water', h2opct)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding global min_h2o_pct_for_wat')
    call cdumlogical(cdum,h2ohgt)
    incstat = nf90_put_att(ncid, nf90_global,'h2o_hgt_over__water',cdum)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding global min_h2o_pct_for_wat')
    incstat = nf90_put_att(ncid, nf90_global, &
                 'middle_dataset_resolution_in_minutes', ntype)
    call checkncerr(incstat,__FILE__,__LINE__,'Error adding global ntype')
    call cdumlogical(cdum,lndfudge)
    incstat = nf90_put_att(ncid, nf90_global,'landuse_fudging', cdum)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding global landuse_fudging')
    if ( lakedpth ) then
      call cdumlogical(cdum,lakfudge)
      incstat = nf90_put_att(ncid, nf90_global,'lake_fudging', cdum)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error adding global lake_fudging')
    end if
    if ( ltexture ) then
      call cdumlogical(cdum,texfudge)
      incstat = nf90_put_att(ncid, nf90_global, 'texture_fudging', cdum)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error adding global texture_fudging')
    end if

    ipnt = 1
    call define_basic_dimensions(ncid,nx,ny,nz,ipnt,idims)
    if ( ltexture ) then
      call add_dimension(ncid,'ntex',ntex,ipnt,idims)
    end if

    call define_horizontal_coord(ncid,nx,ny,xjx,yiy,idims,ihvar)
    call define_vertical_coord(ncid,idims,izvar)

    ipnt = 1
    call define_cross_geolocation_coord(ncid,idims,ipnt,ivar)
    call define_dot_geolocation_coord(ncid,idims,ipnt,ivar)
    call define_topo_and_mask(ncid,idims,ipnt,ivar)
    call define_landuse(ncid,idims,ipnt,ivar)
    call define_mapfactor_and_coriolis(ncid,idims,ipnt,ivar)
    if ( .false. ) then
      call define_initial_snow(ncid,idims,ipnt,ivar)
    end if
    if ( lakedpth ) then
      call define_lakedepth(ncid,idims,ipnt,ivar)
    end if
    if ( ltexture ) then
      call define_textures(ncid,idims,ipnt,ivar,4)
    end if
!
    incstat = nf90_enddef(ncid)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error End Definitions NetCDF output')
!
    call write_vertical_coord(ncid,sigma,hptop,izvar)
    call write_horizontal_coord(ncid,xjx,yiy,ihvar)
    ipnt = 1
    call write_var2d_static(ncid,'xlat',xlat,ipnt,ivar,do_transpose)
    call write_var2d_static(ncid,'xlon',xlon,ipnt,ivar,do_transpose)
    call write_var2d_static(ncid,'dlat',dlat,ipnt,ivar,do_transpose)
    call write_var2d_static(ncid,'dlon',dlon,ipnt,ivar,do_transpose)
    call write_var2d_static(ncid,'topo',htgrid,ipnt,ivar,do_transpose)
    call write_var2d_static(ncid,'mask',mask,ipnt,ivar,do_transpose)
    call write_var2d_static(ncid,'landuse',lndout,ipnt,ivar,do_transpose)
    call write_var2d_static(ncid,'xmap',xmap,ipnt,ivar,do_transpose)
    call write_var2d_static(ncid,'dmap',dmap,ipnt,ivar,do_transpose)
    call write_var2d_static(ncid,'coriol',coriol,ipnt,ivar,do_transpose)
    if ( .false. ) then
      call write_var2d_static(ncid,'snowam',snowam,ipnt,ivar,do_transpose)
    end if
    if (lakedpth) then
      call write_var2d_static(ncid,'dhlake',dpth,ipnt,ivar,do_transpose)
    endif
    if ( ltexture ) then
      call write_var2d_static(ncid,'texout',texout,ipnt,ivar,do_transpose)
      call write_var3d_static(ncid,'frac_tex',frac_tex,ipnt,ivar,do_transpose)
    end if

    incstat = nf90_close(ncid)
    call checkncerr(incstat,__FILE__,__LINE__, &
               ('Error closing NetCDF output '//trim(fname)))

  end subroutine write_domain

end module mod_write

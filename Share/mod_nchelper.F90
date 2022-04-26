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

module mod_nchelper

  use netcdf
  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_constants
  use mod_memutil
  use mod_dynparam
  use mod_message

  implicit none

  private

  public :: cdumlogical

  public :: openfile_withname
  public :: createfile_withname
  public :: closefile
  public :: add_common_global_params
  public :: define_basic_dimensions
  public :: define_horizontal_coord
  public :: define_vertical_coord
  public :: define_cross_geolocation_coord
  public :: write_vertical_coord
  public :: write_vertical_coord_zita
  public :: write_horizontal_coord
  public :: check_dims
  public :: check_var
  public :: ncd_inqdim
  public :: add_dimension
  public :: add_variable
  public :: add_attribute
  public :: get_attribute
  public :: check_dimlen
  public :: checkncerr

#ifdef SINGLE_PRECISION_REAL
  integer(ik4) , public :: regcm_vartype = nf90_real
#else
  integer(ik4) , public :: regcm_vartype = nf90_double
#endif

  interface read_var1d_static
    module procedure read_var1d_static_double_fix
    module procedure read_var1d_static_single_fix
    module procedure read_var1d_static_integer_fix
    module procedure read_var1d_static_double
    module procedure read_var1d_static_single
    module procedure read_var1d_static_integer
    module procedure read_var1d_static_text
  end interface read_var1d_static

  interface read_var2d_static
    module procedure read_var2d_static_double_fix
    module procedure read_var2d_static_single_fix
    module procedure read_var2d_static_integer_fix
    module procedure read_var2d_static_double
    module procedure read_var2d_static_single
    module procedure read_var2d_static_integer
  end interface read_var2d_static

  interface read_var3d_static
    module procedure read_var3d_static_double_fix
    module procedure read_var3d_static_single_fix
    module procedure read_var3d_static_integer_fix
    module procedure read_var3d_static_double
    module procedure read_var3d_static_single
    module procedure read_var3d_static_integer
  end interface read_var3d_static

  interface write_var1d_static
    module procedure write_var1d_static_double
    module procedure write_var1d_static_single
    module procedure write_var1d_static_integer
    module procedure write_var1d_static_text
  end interface write_var1d_static

  interface get_attribute
    module procedure get_attribute_char
    module procedure get_attribute_int4
    module procedure get_attribute_real4
    module procedure get_attribute_real8
  end interface get_attribute

  public :: read_var1d_static
  public :: read_var2d_static
  public :: read_var3d_static

  public :: write_var1d_static
  public :: write_var2d_static
  public :: write_var3d_static

  integer(ik4) :: incstat

  contains
!
  subroutine cdumlogical(cdum,yesno)
    implicit none
    character(len=*) , intent(out) :: cdum
    logical , intent(in) :: yesno
    if (yesno) then
      write (cdum,'(a)') 'Yes'
    else
      write (cdum,'(a)') 'No'
    end if
  end subroutine cdumlogical

  subroutine add_common_global_params(ncid,prgname,lsub)
    implicit none

    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: prgname
    logical :: lsub

    character(len=256) :: history
    real(rk4) , dimension(2) :: trlat
    integer(ik4) , dimension(8) :: tvals

    incstat = nf90_put_att(ncid, nf90_global, 'title',  &
               'ICTP Regional Climatic model V4')
    call checkncerr(incstat,__FILE__,__LINE__,'Error adding global title')
    incstat = nf90_put_att(ncid, nf90_global, 'institution','ICTP')
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding global institution')
    incstat = nf90_put_att(ncid, nf90_global, 'source', &
               'RegCM Model simulation '//trim(prgname)//' output')
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding global source')
    incstat = nf90_put_att(ncid, nf90_global, 'Conventions','CF-1.4')
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding global Conventions')
    call date_and_time(values=tvals)
    write (history,'(i0.4,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a,i0.2,a)')   &
         tvals(1) , '-' , tvals(2) , '-' , tvals(3) , ' ' ,         &
         tvals(5) , ':' , tvals(6) , ':' , tvals(7) ,               &
         ' : Created by RegCM '//trim(prgname)
    incstat = nf90_put_att(ncid, nf90_global, 'history', history)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding global history')
    incstat = nf90_put_att(ncid, nf90_global, 'references', &
               'https://github.com/ICTP/RegCM')
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding global references')
    incstat = nf90_put_att(ncid, nf90_global, 'model_revision',GIT_VER)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding global institution')
    incstat = nf90_put_att(ncid, nf90_global, 'experiment',domname)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding global experiment')
    incstat = nf90_put_att(ncid, nf90_global, 'projection',iproj)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding global projection')
    if ( lsub ) then
      incstat = nf90_put_att(ncid, nf90_global, &
                             'grid_size_in_meters', (ds*1000.0)/real(nsg,rkx))
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error adding global gridsize')
      incstat = nf90_put_att(ncid, nf90_global, 'model_subgrid', 'Yes');
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error adding global subgrid flag')
    else
      incstat = nf90_put_att(ncid, nf90_global,'grid_size_in_meters', ds*1000.0)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error adding global gridsize')
    end if
    incstat = nf90_put_att(ncid, nf90_global, &
                 'latitude_of_projection_origin', clat)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding global clat')
    incstat = nf90_put_att(ncid, nf90_global,                     &
                 'longitude_of_projection_origin', clon)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding global clon')
    if (iproj == 'ROTMER') then
      incstat = nf90_put_att(ncid, nf90_global,'grid_north_pole_latitude', plat)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error adding global plat')
      incstat = nf90_put_att(ncid, nf90_global, &
                             'grid_north_pole_longitude', plon)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error adding global plon')
    else if (iproj == 'LAMCON') then
      trlat(1) = real(truelatl)
      trlat(2) = real(truelath)
      incstat = nf90_put_att(ncid, nf90_global,'standard_parallel', trlat)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error adding global truelat')
    end if
    incstat = nf90_put_att(ncid, nf90_global, 'grid_factor', xcone)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding global grid_factor')
  end subroutine add_common_global_params

  subroutine define_horizontal_coord(ncid,nx,ny,xjx,yiy,idims,ihvar)
    implicit none
    integer(ik4) , intent(in) :: ncid
    integer(ik4) , intent(in) :: nx , ny
    integer(ik4) , intent(in) , dimension(:) :: idims
    integer(ik4) , dimension(2) , intent(out) :: ihvar
    real(rk4) , pointer , dimension(:) , intent(inout) :: xjx
    real(rk4) , pointer , dimension(:) , intent(inout) :: yiy
    integer(ik4) :: i , j

    call getmem1d(yiy,1,ny,'mod_write:yiy')
    call getmem1d(xjx,1,nx,'mod_write:xjx')
    yiy(1) = -real((real(iy-1,rkx)/d_two) * ds,rk4)
    xjx(1) = -real((real(jx-1,rkx)/d_two) * ds,rk4)
    do i = 2 , ny
      yiy(i) = real(real(yiy(i-1),rkx)+ds,rk4)
    end do
    do j = 2 , nx
      xjx(j) = real(real(xjx(j-1),rkx)+ds,rk4)
    end do
    incstat = nf90_def_var(ncid, 'jx', nf90_float, idims(1), ihvar(1))
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding variable jx')
    incstat = nf90_put_att(ncid, ihvar(1), 'standard_name',       &
                           'projection_x_coordinate')
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding jx standard_name')
    incstat = nf90_put_att(ncid, ihvar(1), 'long_name',           &
                           'x-coordinate in Cartesian system')
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding jx long_name')
    incstat = nf90_put_att(ncid, ihvar(1), 'units', 'km')
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding jx units')
    incstat = nf90_def_var(ncid, 'iy', nf90_float, idims(2), ihvar(2))
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding variable iy')
    incstat = nf90_put_att(ncid, ihvar(2), 'standard_name',       &
                           'projection_y_coordinate')
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding iy standard_name')
    incstat = nf90_put_att(ncid, ihvar(2), 'long_name',           &
                           'y-coordinate in Cartesian system')
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding iy long_name')
    incstat = nf90_put_att(ncid, ihvar(2), 'units', 'km')
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding iy units')
  end subroutine define_horizontal_coord

  subroutine define_vertical_coord(ncid,idims,izvar)
    implicit none
    integer(ik4) , intent(in) :: ncid
    integer(ik4) , intent(in) , dimension(:) :: idims
    integer(ik4) , intent(out) , dimension(3) :: izvar

    incstat = nf90_def_var(ncid, 'kz', nf90_double, idims(3), izvar(1))
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding variable kz')
    if ( idynamic < 3 ) then
      incstat = nf90_put_att(ncid, izvar(1), 'standard_name', &
                             'atmosphere_sigma_coordinate')
    else
      incstat = nf90_put_att(ncid, izvar(1), 'standard_name', &
                             'atmosphere_hybrid_height_coordinate')
    end if
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding sigma standard_name')
    incstat = nf90_put_att(ncid, izvar(1), 'long_name', &
                           'Sigma at model layer midpoints')
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding sigma long_name')
    incstat = nf90_put_att(ncid, izvar(1), 'units', '1')
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding sigma units')
    incstat = nf90_put_att(ncid, izvar(1), 'axis', 'Z')
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding sigma axis')
    incstat = nf90_put_att(ncid, izvar(1), 'positive', 'down')
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding sigma positive')
    if ( idynamic == 3 ) then
      incstat = nf90_put_att(ncid, izvar(1), 'formula', &
             'z(k,j,i) = a(k) + b(k) * topo(j,i)')
      call checkncerr(incstat,__FILE__,__LINE__, &
                             'Error adding zeta formula')
    else if ( idynamic == 2 ) then
      incstat = nf90_put_att(ncid, izvar(1), 'formula', &
             'p(n,k,j,i) = ptop + sigma(k)*(p0(j,i)-ptop) + ppa(n,k,j,i)')
      call checkncerr(incstat,__FILE__,__LINE__, &
                             'Error adding sigma formula')
    else
      incstat = nf90_put_att(ncid, izvar(1), 'formula_terms', &
                             'kz: sigma ps: ps ptop: ptop')
      call checkncerr(incstat,__FILE__,__LINE__, &
                             'Error adding sigma formula_terms')
    end if
    if ( idynamic < 3 ) then
      incstat = nf90_def_var(ncid, 'ptop', nf90_double, varid=izvar(2))
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error adding variable ptop')
      incstat = nf90_put_att(ncid, izvar(2), 'standard_name', &
                             'air_pressure')
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error adding ptop standard_name')
      incstat = nf90_put_att(ncid, izvar(2), 'long_name', &
                             'Pressure at model top')
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error adding ptop long_name')
      incstat = nf90_put_att(ncid, izvar(2), 'units', 'hPa')
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error adding ptop units')
    else
      incstat = nf90_def_var(ncid, 'a', nf90_double, idims(3), izvar(2))
      incstat = nf90_put_att(ncid, izvar(2), 'long_name', &
                             'vertical coordinate formula term a(k)')
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error adding a long_name')
      incstat = nf90_def_var(ncid, 'b', nf90_double, idims(3), izvar(3))
      incstat = nf90_put_att(ncid, izvar(3), 'long_name', &
                             'vertical coordinate formula term b(k)')
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error adding b long_name')
    end if
  end subroutine define_vertical_coord

  subroutine define_cross_geolocation_coord(ncid,idims,ipnt,ivar)
    implicit none
    integer(ik4) , intent(in) :: ncid
    integer(ik4) , dimension(:) , intent(in) :: idims
    integer(ik4) , intent(inout) :: ipnt
    integer(ik4) , dimension(:) , intent(out) :: ivar

    incstat = nf90_def_var(ncid, 'xlat', nf90_double, idims(1:2), ivar(ipnt))
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding variable xlat')
#ifdef NETCDF4_HDF5
#if defined (NETCDF4_COMPRESS)
    incstat = nf90_def_var_deflate(ncid, ivar(ipnt), 1, 1, deflate_level)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error setting deflate on xlat')
#endif
#endif
    incstat = nf90_put_att(ncid, ivar(ipnt), 'standard_name', 'latitude')
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding xlat standard_name')
    incstat = nf90_put_att(ncid, ivar(ipnt), 'long_name', &
                           'Latitude at cross points')
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding xlat long_name')
    incstat = nf90_put_att(ncid, ivar(ipnt), 'units', 'degrees_north')
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding xlat units')
    ipnt = ipnt+1
    incstat = nf90_def_var(ncid, 'xlon', nf90_double, idims(1:2), ivar(ipnt))
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding variable xlon')
#ifdef NETCDF4_HDF5
#if defined (NETCDF4_COMPRESS)
    incstat = nf90_def_var_deflate(ncid, ivar(ipnt), 1, 1, deflate_level)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error setting deflate on xlon')
#endif
#endif
    incstat = nf90_put_att(ncid, ivar(ipnt), 'standard_name', 'longitude')
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding xlon standard_name')
    incstat = nf90_put_att(ncid, ivar(ipnt), 'long_name', &
                           'Longitude at cross points')
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding xlon long_name')
    incstat = nf90_put_att(ncid, ivar(ipnt), 'units', 'degrees_east')
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding xlon units')
    ipnt = ipnt + 1
  end subroutine define_cross_geolocation_coord

  subroutine write_vertical_coord(ncid,sigma,ptop,izvar)
    implicit none
    integer(ik4) , intent(in) :: ncid
    real(rk4) , dimension(:) , intent(in) :: sigma
    real(rk4) , intent(in) :: ptop
    integer(ik4) , intent(in) , dimension(3) :: izvar
    incstat = nf90_put_var(ncid, izvar(1), sigma)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error variable sigma write')
    incstat = nf90_put_var(ncid, izvar(2), ptop)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error variable ptop write')
  end subroutine write_vertical_coord

  subroutine write_vertical_coord_zita(ncid,sigma,a,b,izvar)
    implicit none
    integer(ik4) , intent(in) :: ncid
    real(rk4) , dimension(:) , intent(in) :: sigma
    real(rk4) , dimension(:) , intent(in) :: a
    real(rk4) , dimension(:) , intent(in) :: b
    integer(ik4) , intent(in) , dimension(3) :: izvar
    incstat = nf90_put_var(ncid, izvar(1), sigma)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error variable sigma write')
    incstat = nf90_put_var(ncid, izvar(2), a)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error variable a write')
    incstat = nf90_put_var(ncid, izvar(3), b)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error variable b write')
  end subroutine write_vertical_coord_zita

  subroutine write_horizontal_coord(ncid,xjx,yiy,ihvar)
    implicit none
    integer(ik4) , intent(in) :: ncid
    real(rk4) , dimension(:) , intent(in) :: xjx , yiy
    integer(ik4) , intent(in) , dimension(2) :: ihvar
    incstat = nf90_put_var(ncid, ihvar(1), xjx)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error variable jx write')
    incstat = nf90_put_var(ncid, ihvar(2), yiy)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error variable iy write')
  end subroutine write_horizontal_coord

  subroutine write_var1d_static_single(ncid,vnam,values,ipnt,ivar)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) :: vnam
    real(rk4) , dimension(:) , intent(in) :: values
    integer(ik4) , intent(inout) :: ipnt
    integer(ik4) , intent(in) , dimension(:) :: ivar
    incstat = nf90_put_var(ncid, ivar(ipnt), values)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error variable '//vnam//' write')
    ipnt = ipnt + 1
    if ( debug_level > 2 ) then
      incstat = nf90_sync(ncid)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error variable '//vnam//' sync')
    end if
  end subroutine write_var1d_static_single

  subroutine write_var1d_static_double(ncid,vnam,values,ipnt,ivar)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) :: vnam
    real(rk8) , dimension(:) , intent(in) :: values
    integer(ik4) , intent(inout) :: ipnt
    integer(ik4) , intent(in) , dimension(:) :: ivar
    incstat = nf90_put_var(ncid, ivar(ipnt), values)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error variable '//vnam//' write')
    ipnt = ipnt + 1
    if ( debug_level > 2 ) then
      incstat = nf90_sync(ncid)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error variable '//vnam//' sync')
    end if
  end subroutine write_var1d_static_double

  subroutine write_var1d_static_integer(ncid,vnam,values,ipnt,ivar)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) :: vnam
    integer(ik4) , dimension(:) , intent(in) :: values
    integer(ik4) , intent(inout) :: ipnt
    integer(ik4) , intent(in) , dimension(:) :: ivar
    incstat = nf90_put_var(ncid, ivar(ipnt), values)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error variable '//vnam//' write')
    ipnt = ipnt + 1
    if ( debug_level > 2 ) then
      incstat = nf90_sync(ncid)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error variable '//vnam//' sync')
    end if
  end subroutine write_var1d_static_integer

  subroutine write_var1d_static_text(ncid,vnam,values,ipnt,ivar)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) :: vnam
    character(len=*) , intent(in) :: values
    integer(ik4) , intent(inout) :: ipnt
    integer(ik4) , intent(in) , dimension(:) :: ivar
    incstat = nf90_put_var(ncid, ivar(ipnt), values)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error variable '//vnam//' write')
    ipnt = ipnt + 1
    if ( debug_level > 2 ) then
      incstat = nf90_sync(ncid)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error variable '//vnam//' sync')
    end if
  end subroutine write_var1d_static_text

  subroutine write_var2d_static(ncid,vnam,values,ipnt,ivar)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) :: vnam
    real(rk4) , dimension(:,:) , intent(in) :: values
    integer(ik4) , intent(inout) :: ipnt
    integer(ik4) , intent(in) , dimension(:) :: ivar
    incstat = nf90_put_var(ncid, ivar(ipnt), values)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error variable '//vnam//' write')
    ipnt = ipnt + 1
    if ( debug_level > 2 ) then
      incstat = nf90_sync(ncid)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error variable '//vnam//' sync')
    end if
  end subroutine write_var2d_static

  subroutine write_var3d_static(ncid,vnam,values,ipnt,ivar)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    real(rk4) , dimension(:,:,:) , intent(in) :: values
    integer(ik4) , intent(inout) :: ipnt
    integer(ik4) , intent(in) , dimension(:) :: ivar
    integer(ik4) , dimension(3) :: istart
    integer(ik4) , dimension(3) :: icount
    istart(1) = 1
    istart(2) = 1
    icount(1) = ubound(values,1)
    icount(2) = ubound(values,2)
    incstat = nf90_put_var(ncid,ivar(ipnt),values)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error variable '//vnam//' write')
    if ( debug_level > 2 ) then
      incstat = nf90_sync(ncid)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error variable '//vnam//' sync')
    end if
    ipnt = ipnt + 1
  end subroutine write_var3d_static

  subroutine read_var1d_static_text(ncid,vnam,values)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    character(len=*) , dimension(:) :: values
    integer(ik4) :: ivarid
    incstat = nf90_inq_varid(ncid, vnam, ivarid)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error search '//vnam)
    incstat = nf90_get_var(ncid, ivarid, values)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error read '//vnam)
  end subroutine read_var1d_static_text

  subroutine read_var1d_static_single(ncid,vnam,values,lerror)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    real(rk4) , pointer , dimension(:) :: values
    logical , intent(inout) , optional :: lerror
    integer(ik4) :: ivarid
    incstat = nf90_inq_varid(ncid, vnam, ivarid)
    if ( incstat /= nf90_noerr ) then
      if ( present(lerror) ) then
        if ( lerror ) then
          lerror = .false.
          return
        end if
      end if
    end if
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error search '//vnam)
    incstat = nf90_get_var(ncid, ivarid, values)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error read '//vnam)
  end subroutine read_var1d_static_single

  subroutine read_var1d_static_double(ncid,vnam,values,lerror)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    real(rk8) , pointer , dimension(:) :: values
    logical , intent(inout) , optional :: lerror
    integer(ik4) :: ivarid
    incstat = nf90_inq_varid(ncid, vnam, ivarid)
    if ( incstat /= nf90_noerr ) then
      if ( present(lerror) ) then
        if ( lerror ) then
          lerror = .false.
          return
        end if
      end if
    end if
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error search '//vnam)
    incstat = nf90_get_var(ncid, ivarid, values)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error read '//vnam)
  end subroutine read_var1d_static_double

  subroutine read_var1d_static_integer(ncid,vnam,values,lerror)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    integer(ik4) , pointer , dimension(:) :: values
    logical , intent(inout) , optional :: lerror
    integer(ik4) :: ivarid
    incstat = nf90_inq_varid(ncid, vnam, ivarid)
    if ( incstat /= nf90_noerr ) then
      if ( present(lerror) ) then
        if ( lerror ) then
          lerror = .false.
          return
        end if
      end if
    end if
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error search '//vnam)
    incstat = nf90_get_var(ncid, ivarid, values)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error read '//vnam)
  end subroutine read_var1d_static_integer

  subroutine read_var1d_static_double_fix(ncid,vnam,n,values,lerror)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    integer(ik4) , intent(in) :: n
    real(rk8) , dimension(n) :: values
    logical , optional , intent(inout) :: lerror
    integer(ik4) :: ivarid
    incstat = nf90_inq_varid(ncid, vnam, ivarid)
    if ( incstat /= nf90_noerr ) then
      if ( present(lerror) ) then
        if ( lerror ) then
          lerror = .false.
          return
        end if
      end if
    end if
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error search '//vnam)
    incstat = nf90_get_var(ncid, ivarid, values)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error read '//vnam)
  end subroutine read_var1d_static_double_fix

  subroutine read_var1d_static_single_fix(ncid,vnam,n,values,lerror)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    integer(ik4) , intent(in) :: n
    real(rk4) , dimension(n) :: values
    logical , optional , intent(inout) :: lerror
    integer(ik4) :: ivarid
    incstat = nf90_inq_varid(ncid, vnam, ivarid)
    if ( incstat /= nf90_noerr ) then
      if ( present(lerror) ) then
        if ( lerror ) then
          lerror = .false.
          return
        end if
      end if
    end if
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error search '//vnam)
    incstat = nf90_get_var(ncid, ivarid, values)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error read '//vnam)
  end subroutine read_var1d_static_single_fix

  subroutine read_var1d_static_integer_fix(ncid,vnam,n,values,lerror)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    integer(ik4) , intent(in) :: n
    integer(ik4) , dimension(n) :: values
    logical , optional , intent(inout) :: lerror
    integer(ik4) :: ivarid
    incstat = nf90_inq_varid(ncid, vnam, ivarid)
    if ( incstat /= nf90_noerr ) then
      if ( present(lerror) ) then
        if ( lerror ) then
          lerror = .false.
          return
        end if
      end if
    end if
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error search '//vnam)
    incstat = nf90_get_var(ncid, ivarid, values)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error read '//vnam)
  end subroutine read_var1d_static_integer_fix

  subroutine read_var2d_static_double(ncid,vnam,values,lerror,istart,icount)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    real(rk8) , dimension(:,:) :: values
    logical , intent(inout) , optional :: lerror
    integer(ik4) , dimension(:) , optional :: istart , icount
    integer(ik4) :: ivarid
    if ( present(lerror) ) then
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      if ( incstat /= nf90_noerr .and. lerror ) then
        lerror = .false.
        return
      end if
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    else
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    end if
    if ( present(istart) .and. present(icount) ) then
      incstat = nf90_get_var(ncid, ivarid, values, istart, icount)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    else
      incstat = nf90_get_var(ncid, ivarid, values)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    end if
  end subroutine read_var2d_static_double

  subroutine read_var2d_static_single(ncid,vnam,values,lerror,istart,icount)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    real(rk4) , dimension(:,:) :: values
    logical , intent(inout) , optional :: lerror
    integer(ik4) , dimension(:) , optional :: istart , icount
    integer(ik4) :: ivarid
    if ( present(lerror) ) then
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      if ( incstat /= nf90_noerr .and. lerror ) then
        lerror = .false.
        return
      end if
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    else
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    end if
    if ( present(istart) .and. present(icount) ) then
      incstat = nf90_get_var(ncid, ivarid, values, istart, icount)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    else
      incstat = nf90_get_var(ncid, ivarid, values)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    end if
  end subroutine read_var2d_static_single

  subroutine read_var2d_static_integer(ncid,vnam,values,lerror,istart,icount)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    integer(ik4) , dimension(:,:) :: values
    logical , intent(inout) , optional :: lerror
    integer(ik4) , dimension(:) , optional :: istart , icount
    integer(ik4) :: ivarid
    if ( present(lerror) ) then
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      if ( incstat /= nf90_noerr .and. lerror ) then
        lerror = .false.
        return
      end if
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    else
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    end if
    if ( present(istart) .and. present(icount) ) then
      incstat = nf90_get_var(ncid, ivarid, values, istart, icount)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    else
      incstat = nf90_get_var(ncid, ivarid, values)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    end if
  end subroutine read_var2d_static_integer

  subroutine read_var2d_static_double_fix(ncid,vnam,n,m,values, &
                                          lerror,istart,icount)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    integer(ik4) , intent(in) :: n , m
    real(rk8) , dimension(n,m) :: values
    logical , intent(inout) , optional :: lerror
    integer(ik4) , dimension(:) , optional :: istart , icount
    integer(ik4) :: ivarid
    if ( present(lerror) ) then
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      if ( incstat /= nf90_noerr .and. lerror ) then
        lerror = .false.
        return
      end if
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    else
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    end if
    if ( present(istart) .and. present(icount) ) then
      incstat = nf90_get_var(ncid, ivarid, values, istart, icount)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    else
      incstat = nf90_get_var(ncid, ivarid, values)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    end if
  end subroutine read_var2d_static_double_fix

  subroutine read_var2d_static_single_fix(ncid,vnam,n,m,values, &
                                          lerror,istart,icount)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    integer(ik4) , intent(in) :: n , m
    real(rk4) , dimension(n,m) :: values
    logical , intent(inout) , optional :: lerror
    integer(ik4) , dimension(:) , optional :: istart , icount
    integer(ik4) :: ivarid
    if ( present(lerror) ) then
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      if ( incstat /= nf90_noerr .and. lerror ) then
        lerror = .false.
        return
      end if
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    else
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    end if
    if ( present(istart) .and. present(icount) ) then
      incstat = nf90_get_var(ncid, ivarid, values, istart, icount)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    else
      incstat = nf90_get_var(ncid, ivarid, values)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    end if
  end subroutine read_var2d_static_single_fix

  subroutine read_var2d_static_integer_fix(ncid,vnam,n,m,values, &
                                          lerror,istart,icount)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    integer(ik4) , intent(in) :: n , m
    integer(ik4) , dimension(n,m) :: values
    logical , intent(inout) , optional :: lerror
    integer(ik4) , dimension(:) , optional :: istart , icount
    integer(ik4) :: ivarid
    if ( present(lerror) ) then
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      if ( incstat /= nf90_noerr .and. lerror ) then
        lerror = .false.
        return
      end if
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    else
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    end if
    if ( present(istart) .and. present(icount) ) then
      incstat = nf90_get_var(ncid, ivarid, values, istart, icount)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    else
      incstat = nf90_get_var(ncid, ivarid, values)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    end if
  end subroutine read_var2d_static_integer_fix

  subroutine read_var3d_static_double(ncid,vnam,values,lerror,istart,icount)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    real(rk8) , dimension(:,:,:) :: values
    logical , intent(inout) , optional :: lerror
    integer(ik4) , dimension(:) , optional :: istart , icount
    integer(ik4) :: ivarid
    if ( present(lerror) ) then
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      if ( incstat /= nf90_noerr .and. lerror ) then
        lerror = .false.
        return
      end if
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    else
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    end if
    if ( present(istart) .and. present(icount) ) then
      incstat = nf90_get_var(ncid, ivarid, values, istart, icount)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    else
      incstat = nf90_get_var(ncid, ivarid, values)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    end if
  end subroutine read_var3d_static_double

  subroutine read_var3d_static_single(ncid,vnam,values,lerror,istart,icount)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    real(rk4) , dimension(:,:,:) :: values
    logical , intent(inout) , optional :: lerror
    integer(ik4) , dimension(:) , optional :: istart , icount
    integer(ik4) :: ivarid
    if ( present(lerror) ) then
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      if ( incstat /= nf90_noerr .and. lerror ) then
        lerror = .false.
        return
      end if
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    else
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    end if
    if ( present(istart) .and. present(icount) ) then
      incstat = nf90_get_var(ncid, ivarid, values, istart, icount)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    else
      incstat = nf90_get_var(ncid, ivarid, values)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    end if
  end subroutine read_var3d_static_single

  subroutine read_var3d_static_integer(ncid,vnam,values,lerror,istart,icount)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    integer(ik4) , dimension(:,:,:) :: values
    logical , intent(inout) , optional :: lerror
    integer(ik4) , dimension(:) , optional :: istart , icount
    integer(ik4) :: ivarid
    if ( present(lerror) ) then
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      if ( incstat /= nf90_noerr .and. lerror ) then
        lerror = .false.
        return
      end if
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    else
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    end if
    if ( present(istart) .and. present(icount) ) then
      incstat = nf90_get_var(ncid, ivarid, values, istart, icount)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    else
      incstat = nf90_get_var(ncid, ivarid, values)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    end if
  end subroutine read_var3d_static_integer

  subroutine read_var3d_static_double_fix(ncid,vnam,n,m,l,values, &
                                          lerror,istart,icount)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    integer , intent(in) :: n , m , l
    real(rk8) , dimension(n,m,l) :: values
    logical , intent(inout) , optional :: lerror
    integer(ik4) , dimension(:) , optional :: istart , icount
    integer(ik4) :: ivarid
    if ( present(lerror) ) then
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      if ( incstat /= nf90_noerr .and. lerror ) then
        lerror = .false.
        return
      end if
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    else
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    end if
    if ( present(istart) .and. present(icount) ) then
      incstat = nf90_get_var(ncid, ivarid, values, istart, icount)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    else
      incstat = nf90_get_var(ncid, ivarid, values)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    end if
  end subroutine read_var3d_static_double_fix

  subroutine read_var3d_static_single_fix(ncid,vnam,n,m,l,values, &
                                          lerror,istart,icount)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    integer , intent(in) :: n , m , l
    real(rk4) , dimension(n,m,l) :: values
    logical , intent(inout) , optional :: lerror
    integer(ik4) , dimension(:) , optional :: istart , icount
    integer(ik4) :: ivarid
    if ( present(lerror) ) then
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      if ( incstat /= nf90_noerr .and. lerror ) then
        lerror = .false.
        return
      end if
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    else
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    end if
    if ( present(istart) .and. present(icount) ) then
      incstat = nf90_get_var(ncid, ivarid, values, istart, icount)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    else
      incstat = nf90_get_var(ncid, ivarid, values)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    end if
  end subroutine read_var3d_static_single_fix

  subroutine read_var3d_static_integer_fix(ncid,vnam,n,m,l,values, &
                                           lerror,istart,icount)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vnam
    integer , intent(in) :: n , m , l
    integer(ik4) , dimension(n,m,l) :: values
    logical , intent(inout) , optional :: lerror
    integer(ik4) , dimension(:) , optional :: istart , icount
    integer(ik4) :: ivarid
    if ( present(lerror) ) then
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      if ( incstat /= nf90_noerr .and. lerror ) then
        lerror = .false.
        return
      end if
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    else
      incstat = nf90_inq_varid(ncid, vnam, ivarid)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error search '//vnam)
    end if
    if ( present(istart) .and. present(icount) ) then
      incstat = nf90_get_var(ncid, ivarid, values, istart, icount)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    else
      incstat = nf90_get_var(ncid, ivarid, values)
      call checkncerr(incstat,__FILE__,__LINE__, &
                      'Error read '//vnam)
    end if
  end subroutine read_var3d_static_integer_fix

  subroutine define_basic_dimensions(ncid,nx,ny,nz,ipnt,idims)
    implicit none
    integer(ik4) , intent(in) :: ncid
    integer(ik4) , intent(in) :: nx , ny , nz
    integer(ik4) , intent(inout) , dimension(:) :: idims
    integer(ik4) , intent(inout) :: ipnt
    incstat = nf90_def_dim(ncid, 'jx', nx, idims(ipnt))
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding dimension jx')
    incstat = nf90_def_dim(ncid, 'iy', ny, idims(ipnt+1))
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding dimension iy')
    incstat = nf90_def_dim(ncid, 'kz', nz, idims(ipnt+2))
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding dimension kz')
    ipnt = ipnt + 3
  end subroutine define_basic_dimensions

  subroutine add_dimension(ncid,dnam,nd,ipnt,idims)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: dnam
    integer(ik4) , intent(in) :: nd
    integer(ik4) , intent(inout) , dimension(:) :: idims
    integer(ik4) , intent(inout) :: ipnt
    integer(ik4) :: incstat
    incstat = nd
    if ( nd == -1 ) incstat = nf90_unlimited
    incstat = nf90_def_dim(ncid, dnam, incstat, idims(ipnt))
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding dimension '//dnam)
    ipnt = ipnt + 1
  end subroutine add_dimension

  subroutine createfile_withname(fname,ncid)
    implicit none
    character(len=*) , intent(in) :: fname
    integer(ik4) , intent(out) :: ncid
    incstat = nf90_create(fname, iomode, ncid)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error creating NetCDF output '//trim(fname))
  end subroutine createfile_withname

  subroutine openfile_withname(fname,ncid)
    implicit none
    character(len=*) , intent(in) :: fname
    integer(ik4) , intent(out) :: ncid
    incstat = nf90_open(fname, nf90_nowrite, ncid)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error open NetCDF input '//trim(fname))
  end subroutine openfile_withname

  subroutine ncd_inqdim(ncid,dname,dlen,lerror,lexist)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: dname
    integer(ik4) , intent(out) , optional :: dlen
    logical , intent(in) , optional :: lerror
    logical , intent(out) , optional :: lexist
    integer(ik4) :: istatus
    integer(ik4) :: idimid
    istatus = nf90_inq_dimid(ncid, dname, idimid)
    if ( istatus /= nf90_noerr ) then
      if ( present(lexist) ) then
        lexist = .false.
        return
      end if
      if ( present(lerror) ) then
        if ( lerror ) then
          if ( present(dlen) ) then
            dlen = -1
          end if
          return
        end if
      end if
    else
      if ( present(lexist) ) then
        lexist = .true.
        return
      end if
    end if
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error search dimension '//dname)
    if ( present(dlen) ) then
      istatus = nf90_inquire_dimension(ncid, idimid, len=dlen)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read dimension '//dname)
    end if
  end subroutine ncd_inqdim

  subroutine add_attribute(ncid,aname,aval,ivar)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: aname
    character(len=*) , intent(in) :: aval
    integer(ik4) , intent(in) , optional :: ivar
    integer :: istat
    if ( present(ivar) ) then
      istat = nf90_put_att(ncid,ivar,aname,aval)
    else
      istat = nf90_put_att(ncid,nf90_global,aname,aval)
    end if
    call checkncerr(istat,__FILE__,__LINE__, &
                    'Error adding attribute '//aname)
  end subroutine add_attribute

  subroutine get_attribute_char(ncid,aname,aval,ivar)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: aname
    character(len=*) , intent(out) :: aval
    integer(ik4) , intent(in) , optional :: ivar
    integer :: istat
    if ( present(ivar) ) then
      istat = nf90_get_att(ncid,ivar,aname,aval)
    else
      istat = nf90_get_att(ncid,nf90_global,aname,aval)
    end if
    call checkncerr(istat,__FILE__,__LINE__, &
                    'Error reading attribute '//aname)
  end subroutine get_attribute_char

  subroutine get_attribute_real4(ncid,aname,aval,ivar)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: aname
    real(rk4) , intent(out) :: aval
    integer(ik4) , intent(in) , optional :: ivar
    integer :: istat
    if ( present(ivar) ) then
      istat = nf90_get_att(ncid,ivar,aname,aval)
    else
      istat = nf90_get_att(ncid,nf90_global,aname,aval)
    end if
    call checkncerr(istat,__FILE__,__LINE__, &
                    'Error reading attribute '//aname)
  end subroutine get_attribute_real4

  subroutine get_attribute_real8(ncid,aname,aval,ivar)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: aname
    real(rk8) , intent(out) :: aval
    integer(ik4) , intent(in) , optional :: ivar
    integer :: istat
    if ( present(ivar) ) then
      istat = nf90_get_att(ncid,ivar,aname,aval)
    else
      istat = nf90_get_att(ncid,nf90_global,aname,aval)
    end if
    call checkncerr(istat,__FILE__,__LINE__, &
                    'Error reading attribute '//aname)
  end subroutine get_attribute_real8

  subroutine get_attribute_int4(ncid,aname,aval,ivar)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: aname
    integer(ik4) , intent(out) :: aval
    integer(ik4) , intent(in) , optional :: ivar
    integer :: istat
    if ( present(ivar) ) then
      istat = nf90_get_att(ncid,ivar,aname,aval)
    else
      istat = nf90_get_att(ncid,nf90_global,aname,aval)
    end if
    call checkncerr(istat,__FILE__,__LINE__, &
                    'Error reading attribute '//aname)
  end subroutine get_attribute_int4

  logical function check_dimlen(ncid,dname,ival)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: dname
    integer(ik4) , intent(in) :: ival
    integer(ik4) :: idimid , dlen , istatus
    check_dimlen = .false.
    istatus = nf90_inq_dimid(ncid, dname, idimid)
    if ( istatus /= nf90_noerr ) return
    istatus = nf90_inquire_dimension(ncid, idimid, len=dlen)
    if ( istatus /= nf90_noerr ) return
    if ( dlen /= ival ) return
    check_dimlen = .true.
  end function check_dimlen

  subroutine check_dims(ncid)
    implicit none
    integer(ik4) , intent(in) :: ncid
    integer(ik4) :: istatus
    integer(ik4) :: idimid
    integer(ik4) :: iyy , jxx , kzz
    istatus = nf90_inq_dimid(ncid, 'jx', idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error search dimension JX')
    istatus = nf90_inquire_dimension(ncid, idimid, len=jxx)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read dimension JX')
    if ( jx /= jxx ) then
      write(stderr,*) 'DOMAIN FILE : ', jxx
      write(stderr,*) 'NAMELIST    : ', jx
      call die('Mismatch: JX in DOMAIN file /= JX in namelist')
    end if
    istatus = nf90_inq_dimid(ncid, 'iy', idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error search dimension IY')
    istatus = nf90_inquire_dimension(ncid, idimid, len=iyy)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read dimension IY')
    if ( iy /= iyy ) then
      write(stderr,*) 'DOMAIN FILE : ', iyy
      write(stderr,*) 'NAMELIST    : ', iy
      call die('Mismatch: IY in DOMAIN file /= IY in namelist')
    end if
    istatus = nf90_inq_dimid(ncid, 'kz', idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error search dimension KZ')
    istatus = nf90_inquire_dimension(ncid, idimid, len=kzz)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read dimension KZ')
    if ( kz /= kzz ) then
      write(stderr,*) 'DOMAIN FILE : ', kzz
      write(stderr,*) 'NAMELIST    : ', kz
      call die('Mismatch: KZ in DOMAIN file /= KZ in namelist')
    end if
  end subroutine check_dims

  subroutine add_variable(ncid,varname,long_name,units,idims,ipnt,ivars)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: varname
    character(len=*) , intent(in) :: long_name
    character(len=*) , intent(in) :: units
    integer(ik4) , dimension(:) , intent(in) :: idims
    integer(ik4) , intent(inout) :: ipnt
    integer(ik4) , dimension(:) , intent(inout) :: ivars
    integer(ik4) :: incstat
    incstat = nf90_def_var(ncid, varname, regcm_vartype, idims, ivars(ipnt))
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error adding variable '//varname)
#ifdef NETCDF4_HDF5
#if defined (NETCDF4_COMPRESS)
    incstat = nf90_def_var_deflate(ncid, ivars(ipnt), 1, 1, deflate_level)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error setting deflate on xlat')
#endif
#endif
    incstat = nf90_put_att(ncid, ivars(ipnt), 'long_name',long_name)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error long_name to '//varname)
    incstat = nf90_put_att(ncid, ivars(ipnt), 'units', units)
    call checkncerr(incstat,__FILE__,__LINE__, &
                    'Error units to '//varname)
    ipnt = ipnt + 1
  end subroutine add_variable

  subroutine check_var(ncid,vname,lerror)
    implicit none
    integer(ik4) , intent(in) :: ncid
    character(len=*) , intent(in) :: vname
    logical , intent(inout) , optional :: lerror
    integer(ik4) :: istatus
    integer(ik4) :: ivarid
    istatus = nf90_inq_varid(ncid,vname,ivarid)
    if ( istatus /= nf90_noerr ) then
      if ( present(lerror) ) then
        if ( lerror ) then
          lerror = .false.
          return
        end if
      end if
    end if
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error search '//vname)
  end subroutine check_var

  subroutine closefile(ncid)
    implicit none
    integer(ik4) , intent(in) :: ncid
    integer(ik4) :: istatus
    istatus = nf90_close(ncid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error closing file')
  end subroutine closefile

  subroutine checkncerr(ival,filename,line,arg)
    implicit none
    integer(ik4) , intent(in) :: ival , line
    character(len=8) :: cline
    character(*) , intent(in) :: filename , arg
    if ( ival /= nf90_noerr ) then
      write (cline,'(i8)') line
      write (stderr,*) nf90_strerror(ival)
      call die(filename,trim(cline)//':'//arg,ival)
    end if
  end subroutine checkncerr
!
end module mod_nchelper
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

! Unit test of Share/mod_nchelper.F90: a file is written with the module
! (dimensions, coordinates, variables, attributes) in the current
! directory and read back, including netCDF-4 string attributes written
! by other tools, missing variables and dimensions with lerror and lexist,
! and the check of the domain dimensions.

module test_nchelper_c
  use, intrinsic :: iso_c_binding
  implicit none
  interface
    integer(c_int) function nc_put_att_string(ncid,varid,name,n,op) bind(c)
      import :: c_int, c_char, c_size_t, c_ptr
      integer(c_int), value :: ncid, varid
      character(kind=c_char), dimension(*), intent(in) :: name
      integer(c_size_t), value :: n
      type(c_ptr), dimension(*), intent(in) :: op
    end function nc_put_att_string
  end interface
  contains
  ! A string attribute as netCDF4-python writes them; varid is the
  ! Fortran one (nf90_global = 0).
  integer function put_string_att(ncid,varid,name,val) result(istat)
    implicit none
    integer, intent(in) :: ncid, varid
    character(len=*), intent(in) :: name, val
    character(kind=c_char,len=len(val)+1), target :: cval
    type(c_ptr), dimension(1) :: pp
    cval = val//c_null_char
    pp(1) = c_loc(cval)
    istat = nc_put_att_string(ncid,varid-1,name//c_null_char, &
                              1_c_size_t,pp)
  end function put_string_att
end module test_nchelper_c

program test_nchelper
  use netcdf
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_memutil
  use mod_nchelper
  use test_nchelper_c
  use share_test_utils
  implicit none
  character(len=*), parameter :: fname = 'test_nchelper.nc'
  integer(ik4), parameter :: nx = 6, ny = 4, nz = 3
  integer(ik4) :: ncid, ipnt, dlen, istat, ivid
  integer(ik4), dimension(4) :: idims
  integer(ik4), dimension(4) :: ivar
  integer(ik4), dimension(2) :: ihvar
  integer(ik4), dimension(3) :: izvar
  real(rk4), pointer, contiguous, dimension(:) :: xjx => null()
  real(rk4), pointer, contiguous, dimension(:) :: yiy => null()
  real(rk4), dimension(nx,ny) :: v2
  real(rk8), dimension(nx,ny) :: r2
  real(rk4), dimension(nx,ny,nz) :: v3
  real(rk8), dimension(nx,ny,nz) :: r3
  real(rk4), dimension(nz) :: sig
  real(rk8), dimension(nz) :: rs
  real(rk4) :: ra4
  real(rk8) :: ra8
  integer(ik4) :: ia
  character(len=64) :: ca
  logical :: lerror, lexist
  integer(ik4) :: i, j, k

  call memory_init
  jx = nx
  iy = ny
  kz = nz
  ds = 25.0_rkx
  nsg = 3
  idynamic = 1
  iomode = ior(nf90_clobber,nf90_netcdf4)

  do k = 1, nz
    sig(k) = 0.25_rk4*real(k,rk4)
    do j = 1, ny
      do i = 1, nx
        v3(i,j,k) = real(i+10*j+100*k,rk4)
      end do
    end do
  end do
  v2 = v3(:,:,1)

  call createfile_withname(fname,ncid)
  ipnt = 1
  call define_basic_dimensions(ncid,nx,ny,nz,ipnt,idims)
  call add_dimension(ncid,'time',-1,ipnt,idims)
  call define_horizontal_coord(ncid,nx,ny,xjx,yiy,idims,ihvar)
  call define_vertical_coord(ncid,idims,izvar)
  ipnt = 1
  call add_variable(ncid,'v2','two dimensional','K',idims(1:2),ipnt,ivar)
  call add_variable(ncid,'v3','three dimensional','K',idims(1:3), &
                    ipnt,ivar)
  call add_attribute(ncid,'note','global text')
  call add_attribute(ncid,'comment','variable text',ivar(1))
  istat = nf90_put_att(ncid,nf90_global,'count',42)
  istat = nf90_put_att(ncid,nf90_global,'scale',2.5_rk8)
  call check(put_string_att(ncid,nf90_global,'sglobal','string global') &
             == nf90_noerr,'write a global string attribute')
  call check(put_string_att(ncid,ivar(2),'svar','string variable') &
             == nf90_noerr,'write a variable string attribute')
  istat = nf90_enddef(ncid)
  call write_horizontal_coord(ncid,xjx,yiy,ihvar)
  call write_vertical_coord(ncid,sig,5.0_rk4,izvar)
  ipnt = 1
  call write_var2d_static(ncid,'v2',v2,ipnt,ivar)
  call write_var3d_static(ncid,'v3',v3,ipnt,ivar)
  call check(ipnt == 3,'write_var*_static advance the pointer')
  call closefile(ncid)

  call openfile_withname(fname,ncid)
  call check_dims(ncid)
  call check(.true.,'domain dimensions match')
  call ncd_inqdim(ncid,'kz',dlen)
  call check(dlen == nz,'ncd_inqdim length')
  call ncd_inqdim(ncid,'jx',dlen,lexist=lexist)
  call check(lexist .and. dlen == nx,'ncd_inqdim with lexist')
  call ncd_inqdim(ncid,'nothere',dlen,lexist=lexist)
  call check(.not. lexist .and. dlen == -1,'ncd_inqdim, missing')
  call ncd_inqdim(ncid,'nothere',dlen,lerror=.true.)
  call check(dlen == -1,'ncd_inqdim, missing with lerror')
  call check(check_dimlen(ncid,'iy',ny) .and. &
             .not. check_dimlen(ncid,'iy',ny+1),'check_dimlen')

  call read_var2d_static(ncid,'v2',r2)
  call check(all(abs(r2-real(v2,rk8)) < 1.0e-6_rk8),'read 2-D variable')
  call read_var3d_static(ncid,'v3',nx,ny,nz,r3)
  call check(all(abs(r3-real(v3,rk8)) < 1.0e-6_rk8),'read 3-D variable')
  call read_var3d_static(ncid,'v3',r3(1:2,1:2,1:1),istart=[2,3,2], &
                         icount=[2,2,1])
  call check(all(abs(r3(1:2,1:2,1)-real(v3(2:3,3:4,2),rk8)) &
             < 1.0e-6_rk8),'read part of a 3-D variable')
  call read_var1d_static(ncid,'kz',nz,rs)
  call check(all(abs(rs-real(sig,rk8)) < 1.0e-6_rk8), &
    'read the vertical coordinate')
  lerror = .true.
  call read_var2d_static(ncid,'missing',r2,lerror)
  call check(.not. lerror,'missing variable with lerror')
  lerror = .true.
  call check_var(ncid,'v3',lerror)
  call check(lerror,'check_var, existing variable')
  call check_var(ncid,'missing',lerror)
  call check(.not. lerror,'check_var, missing variable')

  istat = nf90_inq_varid(ncid,'kz',ivid)
  call get_attribute(ncid,'formula_terms',ca,ivid)
  call check(ca == 'sigma: kz ps: ps ptop: ptop', &
    'CF formula_terms of the sigma coordinate')

  call read_var1d_static(ncid,'jx',nx,rs(1:1))
  call check(abs(rs(1)+2.5_rk8*ds) < 1.0e-4_rk8, &
    'x coordinate, centred on the domain')

  call get_attribute(ncid,'note',ca)
  call check(ca == 'global text','global text attribute')
  call get_attribute(ncid,'comment',ca,ivar(1))
  call check(ca == 'variable text','variable text attribute')
  call get_attribute(ncid,'count',ia)
  call check(ia == 42,'integer attribute')
  call get_attribute(ncid,'scale',ra8)
  call get_attribute(ncid,'scale',ra4)
  call check(ra8 == 2.5_rk8 .and. ra4 == 2.5_rk4,'real attributes')
  call get_attribute(ncid,'sglobal',ca)
  call check(ca == 'string global','global string attribute')
  call get_attribute(ncid,'svar',ca,ivar(2))
  call check(ca == 'string variable','variable string attribute')
  istat = nf90_inq_varid(ncid,'v3',ivid)
  istat = nf90_myget_att_text(ncid,ivid,'svar',ca)
  call check(istat == nf90_noerr .and. ca == 'string variable', &
    'nf90_myget_att_text')
  istat = nf90_myget_att_text(ncid,ivid,'units',ca)
  call check(istat == nf90_noerr .and. ca == 'K', &
    'nf90_myget_att_text, text attribute')
  call closefile(ncid)

  ! The coordinates of a subgrid file use the subgrid distance.
  call createfile_withname(fname,ncid)
  ipnt = 1
  call define_basic_dimensions(ncid,nx*nsg,ny*nsg,nz,ipnt,idims)
  call define_horizontal_coord(ncid,nx*nsg,ny*nsg,xjx,yiy,idims,ihvar)
  call check(abs(xjx(2)-xjx(1)-real(ds,rk4)/3.0_rk4) < 1.0e-4_rk4 .and. &
             abs(xjx(1)+xjx(nx*nsg)) < 1.0e-4_rk4 .and. &
             abs(yiy(1)+yiy(ny*nsg)) < 1.0e-4_rk4, &
    'subgrid coordinates')
  call closefile(ncid)

  call memory_destroy
  call test_end('test_nchelper')
end program test_nchelper

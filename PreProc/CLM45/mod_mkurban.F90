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
module mod_mkurban
  use mod_realkinds
  use mod_intkinds
  use mod_dynparam
  use mod_grid
  use mod_getwindow
  use mod_bilinear
  use mod_nchelper
  use mod_memutil
  use netcdf

  implicit none

  private

  public :: mkurban_base , mkurban_param
  public :: ip2d , ip3d , ip4d

  character(len=16) , parameter :: latdim = 'lat'
  character(len=16) , parameter :: londim = 'lon'
  character(len=16) , parameter :: densdim = 'density_class'
  character(len=16) , parameter :: regiondim = 'region'
  character(len=16) , parameter :: levurbdim = 'nlevurb'
  character(len=16) , parameter :: solardim = 'numsolar'
  character(len=16) , parameter :: raddim = 'numrad'
  character(len=16) , parameter :: latvar = 'LAT'
  character(len=16) , parameter :: lonvar = 'LON'
  character(len=16) , parameter :: varname = 'PCT_URBAN'
  character(len=16) , parameter :: maskname = 'LANDMASK'
  character(len=16) , parameter :: regionname = 'REGION_ID'

  integer(ik4) , public , parameter :: npu2d = 14
  integer(ik4) , public , parameter :: npu3d = 6
  integer(ik4) , public , parameter :: npu4d = 4
  integer , parameter :: nparam = npu2d + npu3d + npu4d

  character(len=16) , dimension(npu2d) , public , parameter :: parm2d = &
    (/'CANYON_HWR      ', 'EM_IMPROAD      ', 'EM_PERROAD      ', &
      'EM_ROOF         ', 'EM_WALL         ', 'HT_ROOF         ', &
      'NLEV_IMPROAD    ', 'THICK_ROOF      ', 'THICK_WALL      ', &
      'T_BUILDING_MAX  ', 'T_BUILDING_MIN  ', 'WIND_HGT_CANYON ', &
      'WTLUNIT_ROOF    ', 'WTROAD_PERV     '/)
  character(len=36) , dimension(npu2d) , public , parameter :: lngn2d = &
    (/'canyon height to width ratio        ', &
      'emissivity of impervious road       ', &
      'emissivity of pervious road         ', &
      'emissivity of roof                  ', &
      'emissivity of wall                  ', &
      'height of roof                      ', &
      'number of impervious road layers    ', &
      'thickness of roof                   ', &
      'thickness of wall                   ', &
      'maximum intern building temperature ', &
      'minimum intern building temperature ', &
      'height of wind in canyon            ', &
      'fraction of roof                    ', &
      'fraction of pervious road           '/)
  character(len=4) , dimension(npu2d) , public , parameter :: unit2d = &
    (/'1   ', '1   ', '1   ', '1   ', '1   ', 'm   ', &
      '1   ', 'm   ', 'm   ', 'K   ', 'K   ', 'm   ', &
      '1   ', '1   '/)
  character(len=16) , dimension(npu3d) , public , parameter :: parm3d = &
    (/'CV_IMPROAD      ', 'CV_ROOF         ', 'CV_WALL         ', &
      'TK_IMPROAD      ', 'TK_ROOF         ', 'TK_WALL         '/)
  character(len=36) , dimension(npu3d) , public , parameter :: lngn3d = &
    (/'vol heat capacity of impervious road', &
      'vol heat capacity of roof           ', &
      'vol heat capacity of wall           ', &
      'thermal conductivity of imperv road ', &
      'thermal conductivity of roof        ', &
      'thermal conductivity of wall        '/)
  character(len=8) , dimension(npu3d) , public , parameter :: unit3d = &
    (/'J/m^3*K ' , 'J/m^3*K ', 'J/m^3*K ', &
      'W/m*K   ' , 'W/m*K   ', 'W/m*K   '/)
  character(len=16) , dimension(npu4d) , public , parameter :: parm4d = &
    (/'ALB_IMPROAD     ', 'ALB_PERROAD     ', 'ALB_ROOF        ', &
      'ALB_WALL        '/)
  character(len=36) , dimension(npu4d) , public , parameter :: lngn4d = &
    (/'albedo of impervious road           ', &
      'albedo of pervious road             ', &
      'albedo of roof                      ', &
      'albedo of wall                      '/)
  character(len=2) , dimension(npu4d) , public , parameter :: unit4d = &
    (/'1 ', '1 ', '1 ', '1 '/)
  character(len=16) , dimension(nparam) :: parmname

  integer :: nreg , ndens

  integer(ik4) , dimension(nparam) , parameter :: parmdim = &
    (/3,3,3,3,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,5,5,5,5/)
  real(rkx) :: vmin = 0.0_rkx
  real(rkx) :: vmisdat = -9999.0_rkx

  contains

  subroutine mkurban_base(urbanfile,urban)
    implicit none
    character(len=*) , intent(in) :: urbanfile
    real(rkx) , dimension(:,:,:) , intent(out) :: urban
    integer(ik4) :: nlat , nlon
    integer(ik4) :: idimid , ivarid , ivarmask , ncid
    integer(ik4) , dimension(3) :: istart , icount
    integer(ik4) :: istatus , i , ip , ipt , li , lo
    real(rkx) , dimension(:,:) , allocatable :: rmask
    real(rkx) , dimension(:,:,:) , allocatable :: rvar
    real(rkx) , dimension(:) , allocatable :: glat , glon , rlat , rlon
    type(global_domain) :: domain

    character(len=256) :: inpfile

    ipt = 1
    do ip = 1 , npu2d
      parmname(ipt) = parm2d(ip)
      ipt = ipt + 1
    end do
    do ip = 1 , npu3d
      parmname(ipt) = parm3d(ip)
      ipt = ipt + 1
    end do
    do ip = 1 , npu4d
      parmname(ipt) = parm4d(ip)
      ipt = ipt + 1
    end do
    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//urbanfile
    istatus = nf90_open(inpfile,nf90_nowrite,ncid)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Cannot open file '//trim(inpfile))

    istatus = nf90_inq_dimid(ncid,latdim,idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find dimension lat in file '//trim(inpfile))
    istatus = nf90_inquire_dimension(ncid,idimid,len=nlat)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read dimension lat in file '//trim(inpfile))

    allocate(glat(nlat))

    istatus = nf90_inq_varid(ncid,latvar,ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable lat/LAT in file '//trim(inpfile))
    istatus = nf90_get_var(ncid,ivarid,glat)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read variable lat in file '//trim(inpfile))

    istatus = nf90_inq_dimid(ncid,londim,idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find dimension lon in file '//trim(inpfile))
    istatus = nf90_inquire_dimension(ncid,idimid,len=nlon)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read dimension lon in file '//trim(inpfile))

    allocate(glon(nlon))

    istatus = nf90_inq_varid(ncid,lonvar,ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable lon in file '//trim(inpfile))
    istatus = nf90_get_var(ncid,ivarid,glon)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read variable lon in file '//trim(inpfile))
    ! Put longitudes in -180 - 180 range
    where ( glon >  180.0_rkx )
      glon = glon - 360.0_rkx
    end where
    where ( glon < -180.0_rkx )
      glon = glon + 360.0_rkx
    end where

    istatus = nf90_inq_dimid(ncid,densdim,idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find dimension density_class in file '//trim(inpfile))
    istatus = nf90_inquire_dimension(ncid,idimid,len=ndens)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read dimension density_class in file '//trim(inpfile))

    istatus = nf90_inq_dimid(ncid,regiondim,idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find dimension region in file '//trim(inpfile))
    istatus = nf90_inquire_dimension(ncid,idimid,len=nreg)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read dimension region in file '//trim(inpfile))

    call get_window(glat,glon,domain)

    allocate(rvar(sum(domain%ni),domain%nj,ndens))
    allocate(rmask(sum(domain%ni),domain%nj))
    allocate(rlon(sum(domain%ni)))
    allocate(rlat(domain%nj))

    istatus = nf90_inq_varid(ncid,varname,ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable urban in file '//trim(inpfile))
    istatus = nf90_inq_varid(ncid,maskname,ivarmask)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable landmask in file '//trim(inpfile))

    li = 1
    do i = 1 , domain%ntiles
      istart(1) = domain%igstart(i)
      icount(1) = domain%ni(i)
      istart(2) = domain%jgstart
      icount(2) = domain%nj
      istart(3) = 1
      icount(3) = ndens
      lo = li+domain%ni(i)-1
      istatus = nf90_get_var(ncid,ivarid,rvar(li:lo,:,:),istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
        'Cannot read variable urban from file '//trim(inpfile))
      istatus = nf90_get_var(ncid,ivarmask,rmask(li:lo,:), &
        istart(1:2),icount(1:2))
      call checkncerr(istatus,__FILE__,__LINE__, &
        'Cannot read variable mask from file '//trim(inpfile))
      rlon(li:lo) = glon(domain%igstart(i):domain%igstop(i))
      li = li + domain%ni(i)
    end do
    rlat = glat(domain%jgstart:domain%jgstop)

    call bilinear(rvar,rmask,rlon,rlat,urban,xlon,xlat,vmin,vmisdat)

    deallocate(glat,glon,rlat,rlon,rvar,rmask)
  end subroutine mkurban_base

  subroutine mkurban_param(urbanfile,urban3d,urban4d,urban5d)
    implicit none
    character(len=*) , intent(in) :: urbanfile
    real(rkx) , dimension(:,:,:,:) , intent(out) :: urban3d
    real(rkx) , dimension(:,:,:,:,:) , intent(out) :: urban4d
    real(rkx) , dimension(:,:,:,:,:,:) , intent(out) :: urban5d
    integer(ik4) :: nlat , nlon
    integer(ik4) :: idimid , ivarid , ivarmask , ivarreg , ncid
    integer(ik4) , dimension(5) :: istart , icount
    integer(ik4) :: istatus , i , j , ic , il , ir , is , n , li , lo
    integer(ik4) :: nurb , nrad , nsol
    real(rkx) , dimension(:,:,:) , allocatable :: rvar3d
    real(rkx) , dimension(:,:,:,:) , allocatable :: rvar4d
    real(rkx) , dimension(:,:,:,:,:) , allocatable :: rvar5d
    integer(ik4) :: i4 , i5 , i6
    real(rkx) , dimension(:,:) , allocatable :: mread
    integer(ik4) , dimension(:,:) , allocatable :: region
    real(rkx) , dimension(:,:) , allocatable :: rmask
    real(rkx) , dimension(:) , allocatable :: glat , glon , rlat , rlon
    type(global_domain) :: domain

    character(len=256) :: inpfile

    inpfile = trim(inpglob)//pthsep//'CLM45'// &
                             pthsep//'surface'//pthsep//urbanfile
    istatus = nf90_open(inpfile,nf90_nowrite,ncid)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Cannot open file '//trim(inpfile))

    istatus = nf90_inq_dimid(ncid,latdim,idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find dimension lat in file '//trim(inpfile))
    istatus = nf90_inquire_dimension(ncid,idimid,len=nlat)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read dimension lat in file '//trim(inpfile))

    allocate(glat(nlat))

    istatus = nf90_inq_varid(ncid,latvar,ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable lat/LAT in file '//trim(inpfile))
    istatus = nf90_get_var(ncid,ivarid,glat)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read variable lat in file '//trim(inpfile))

    istatus = nf90_inq_dimid(ncid,londim,idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find dimension lon in file '//trim(inpfile))
    istatus = nf90_inquire_dimension(ncid,idimid,len=nlon)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read dimension lon in file '//trim(inpfile))

    allocate(glon(nlon))

    istatus = nf90_inq_varid(ncid,lonvar,ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable lon in file '//trim(inpfile))
    istatus = nf90_get_var(ncid,ivarid,glon)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read variable lon in file '//trim(inpfile))

    ! Put longitudes in -180 - 180 range
    where ( glon >  180.0_rkx )
      glon = glon - 360.0_rkx
    end where
    where ( glon < -180.0_rkx )
      glon = glon + 360.0_rkx
    end where

    istatus = nf90_inq_dimid(ncid,levurbdim,idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find dimension nlevurb in file '//trim(inpfile))
    istatus = nf90_inquire_dimension(ncid,idimid,len=nurb)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read dimension nlevurb in file '//trim(inpfile))

    istatus = nf90_inq_dimid(ncid,solardim,idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find dimension numsolar in file '//trim(inpfile))
    istatus = nf90_inquire_dimension(ncid,idimid,len=nsol)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read dimension numsolar in file '//trim(inpfile))

    istatus = nf90_inq_dimid(ncid,raddim,idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find dimension numrad in file '//trim(inpfile))
    istatus = nf90_inquire_dimension(ncid,idimid,len=nrad)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot read dimension numrad in file '//trim(inpfile))

    call get_window(glat,glon,domain)

    allocate(region(sum(domain%ni),domain%nj))
    allocate(mread(ndens,0:nreg))
    allocate(rvar3d(sum(domain%ni),domain%nj,ndens))
    allocate(rvar4d(sum(domain%ni),domain%nj,nurb,ndens))
    allocate(rvar5d(sum(domain%ni),domain%nj,nrad,nsol,ndens))
    allocate(rmask(sum(domain%ni),domain%nj))
    allocate(rlon(sum(domain%ni)))
    allocate(rlat(domain%nj))

    mread(:,0) = vmisdat

    istatus = nf90_inq_varid(ncid,maskname,ivarmask)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable landmask in file '//trim(inpfile))

    istatus = nf90_inq_varid(ncid,regionname,ivarreg)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Cannot find variable region in file '//trim(inpfile))

    li = 1
    do i = 1 , domain%ntiles
      istart(1) = domain%igstart(i)
      icount(1) = domain%ni(i)
      istart(2) = domain%jgstart
      icount(2) = domain%nj
      lo = li+domain%ni(i)-1
      istatus = nf90_get_var(ncid,ivarmask,rmask(li:lo,:), &
                             istart(1:2),icount(1:2))
      call checkncerr(istatus,__FILE__,__LINE__, &
        'Cannot read variable mask from file '//trim(inpfile))
      istatus = nf90_get_var(ncid,ivarreg,region(li:lo,:), &
                             istart(1:2),icount(1:2))
      call checkncerr(istatus,__FILE__,__LINE__, &
        'Cannot read variable mask from file '//trim(inpfile))
      rlon(li:lo) = glon(domain%igstart(i):domain%igstop(i))
      li = li + domain%ni(i)
    end do
    rlat = glat(domain%jgstart:domain%jgstop)

    istart(1) = 1
    icount(1) = ndens
    istart(2) = 1
    icount(2) = nreg
    do n = 1 , nparam
      istatus = nf90_inq_varid(ncid,parmname(n),ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
        'Cannot find variable '//trim(parmname(n))// &
        ' in file '//trim(inpfile))
      select case (parmdim(n))
        case (3)
          istatus = nf90_get_var(ncid,ivarid,mread(:,1:ndens), &
                          istart(1:2),icount(1:2))
          call checkncerr(istatus,__FILE__,__LINE__,'Cannot read variable '// &
            trim(parmname(n))//' from file '//trim(inpfile))
          do ic = 1 , ndens
            do i = 1 , size(rvar3d,2)
              do j = 1 , size(rvar3d,1)
                rvar3d(j,i,ic) = mread(ic,region(j,i))
              end do
            end do
            where ( rvar3d(:,:,ic) < 0.0 )
              rvar3d(:,:,ic) = minval(mread(ic,1:nreg))
            end where
          end do
          i4 = ip2d(parmname(n))
          call bilinear(rvar3d,rmask,rlon,rlat,urban3d(:,:,:,i4), &
                   xlon,xlat,vmin,vmisdat)
          do ic = 1, ndens
            where ( urban3d(:,:,ic,i4) < 0.0 )
              urban3d(:,:,ic,i4) = minval(rvar3d(:,:,ic))
            end where
          end do
        case (4)
          do il = 1 , nurb
            istart(3) = il
            icount(3) = 1
            istatus = nf90_get_var(ncid,ivarid,mread(:,1:ndens), &
                            istart(1:3),icount(1:3))
            call checkncerr(istatus,__FILE__,__LINE__, &
              'Cannot read variable '// &
              trim(parmname(n))//' from file '//trim(inpfile))
            do ic = 1, ndens
              do i = 1 , size(rvar4d,2)
                do j = 1 , size(rvar4d,1)
                  rvar4d(j,i,il,ic) = mread(ic,region(j,i))
                end do
              end do
              where ( rvar4d(:,:,il,ic) < 0.0 )
                rvar4d(:,:,il,ic) = minval(mread(ic,1:nreg))
              end where
            end do
          end do
          i5 = ip3d(parmname(n))
          call bilinear(rvar4d,rmask,rlon,rlat,urban4d(:,:,:,:,i5), &
                   xlon,xlat,vmin,vmisdat)
          do il = 1 , nurb
            do ic = 1, ndens
              where ( urban4d(:,:,il,ic,i5) < 0.0 )
                urban4d(:,:,il,ic,i5) = minval(rvar4d(:,:,il,ic))
              end where
            end do
          end do
        case (5)
          do is = 1 , nsol
            istart(4) = is
            icount(4) = 1
            do ir = 1 , nrad
              istart(3) = ir
              icount(3) = 1
              istatus = nf90_get_var(ncid,ivarid,mread(:,1:ndens), &
                              istart(1:4),icount(1:4))
              call checkncerr(istatus,__FILE__,__LINE__, &
                'Cannot read variable '// &
                trim(parmname(n))//' from file '//trim(inpfile))
              do ic = 1 , ndens
                do i = 1 , size(rvar5d,2)
                  do j = 1 , size(rvar5d,1)
                    rvar5d(j,i,ir,is,ic) = mread(ic,region(j,i))
                  end do
                end do
                where ( rvar5d(:,:,ir,is,ic) < 0.0 )
                  rvar5d(:,:,ir,is,ic) = minval(mread(ic,1:nreg))
                end where
              end do
            end do
          end do
          i6 = ip4d(parmname(n))
          call bilinear(rvar5d,rmask,rlon,rlat,urban5d(:,:,:,:,:,i6), &
                   xlon,xlat,vmin,vmisdat)
          do is = 1 , nsol
            do ir = 1 , nrad
              do ic = 1 , ndens
                where ( urban5d(:,:,ir,is,ic,i6) < 0.0 )
                  urban5d(:,:,ir,is,ic,i6) = minval(rvar5d(:,:,ir,is,ic))
                end where
              end do
            end do
          end do
        case default
          call die(__FILE__,'Variable dimension not implemented',__LINE__)
      end select
    end do

    deallocate(glat,glon,rlat,rlon,region,mread,rvar3d,rvar4d,rvar5d,rmask)
  end subroutine mkurban_param

  integer(ik4) function ip2d(pname) result(ip)
    implicit none
    character(len=*) :: pname
    do ip = 1 , npu2d
      if ( pname == parm2d(ip) ) then
        return
      end if
    end do
    call die(__FILE__,'Variable '//pname//' NOT FOUND',__LINE__)
  end function ip2d

  integer(ik4) function ip3d(pname) result(ip)
    implicit none
    character(len=*) :: pname
    do ip = 1 , npu3d
      if ( pname == parm3d(ip) ) then
        return
      end if
    end do
    call die(__FILE__,'Variable '//pname//' NOT FOUND',__LINE__)
  end function ip3d

  integer(ik4) function ip4d(pname) result(ip)
    implicit none
    character(len=*) :: pname
    do ip = 1 , npu4d
      if ( pname == parm4d(ip) ) then
        return
      end if
    end do
    call die(__FILE__,'Variable '//pname//' NOT FOUND',__LINE__)
  end function ip4d

end module mod_mkurban
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

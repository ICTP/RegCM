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

! Unit test of Share/mod_domain.F90: domain files are written with
! mod_ncstream as terrain writes them (hydrostatic core, MOLOCH, subgrid)
! in the current directory, checked with check_domain (also as a file of
! surface fields for another sigma core) and read back with the three
! forms of read_domain.

program test_domain
  use netcdf
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_memutil
  use mod_ncstream
  use mod_nchelper, only : openfile_withname, closefile
  use mod_domain
  use share_test_utils
  implicit none
  integer(ik4), parameter :: nj = 12, ni = 10, nk = 5
  real(rkx), pointer, contiguous, dimension(:) :: sigma => null()
  real(rkx), pointer, contiguous, dimension(:) :: sig => null()
  real(rkx), pointer, contiguous, dimension(:,:) :: rlat => null()
  real(rkx), pointer, contiguous, dimension(:,:) :: rlon => null()
  real(rkx), pointer, contiguous, dimension(:,:) :: rht => null()
  real(rkx), pointer, contiguous, dimension(:,:) :: rdlat => null()
  real(rkx), pointer, contiguous, dimension(:,:) :: rulat => null()
  real(rkx), pointer, contiguous, dimension(:,:,:) :: rbnds => null()
  real(rkx), target, dimension(nj,ni) :: xlat, xlon, topo, mask, lu, cor
  real(rkx), target, dimension(nj,ni) :: ones
  real(rkx), target, dimension(nj*3,ni*3) :: sxlat, sxlon, stopo, smask
  real(rkx), target, dimension(4,nj,ni) :: bnds
  real(rkx), pointer, contiguous, dimension(:,:) :: sl => null()
  integer(ik4) :: ncid, i, j

  call memory_init
  jx = nj
  iy = ni
  kz = nk
  kzp1 = nk + 1
  ds = 25.0_rkx
  nsg = 3
  domname = 'test'
  calendar = 'gregorian'
  ical = 1
  iproj = 'LAMCON'
  xcone = 0.71_rkx
  clat = 43.0_rkx
  clon = 15.0_rkx
  cntrj = real(nj+1,rkx)*0.5_rkx
  cntri = real(ni+1,rkx)*0.5_rkx
  ptop = 5.0_rkx
  truelatl = 30.0_rkx
  truelath = 60.0_rkx
  allocate(sigma(nk+1))
  sigma = [0.0_rkx, 0.3_rkx, 0.6_rkx, 0.8_rkx, 0.95_rkx, 1.0_rkx]
  do i = 1, ni
    do j = 1, nj
      xlat(j,i) = 40.0_rkx + 0.25_rkx*real(i,rkx)
      xlon(j,i) = 10.0_rkx + 0.25_rkx*real(j,rkx)
      topo(j,i) = real(10*j+i,rkx)
      mask(j,i) = real(mod(i+j,2),rkx)
      lu(j,i) = real(1+mod(j,5),rkx)
      cor(j,i) = 1.0e-4_rkx
      ones(j,i) = 1.0_rkx
      bnds(:,j,i) = xlat(j,i) + [-0.1_rkx, -0.1_rkx, 0.1_rkx, 0.1_rkx]
    end do
  end do
  do i = 1, 3*ni
    do j = 1, 3*nj
      sxlat(j,i) = 40.0_rkx + real(i,rkx)/12.0_rkx
      sxlon(j,i) = 10.0_rkx + real(j,rkx)/12.0_rkx
      stopo(j,i) = real(j+i,rkx)
      smask(j,i) = 1.0_rkx
    end do
  end do

  ! Hydrostatic core.
  idynamic = 1
  call write_domain('test_domain.nc',.false.)
  call openfile_withname('test_domain.nc',ncid)
  call check_domain(ncid)
  call check(.true.,'check_domain, hydrostatic domain')
  xcone = 0.0_rkx
  call read_domain(ncid)
  call check(abs(xcone-0.71_rkx) < 1.0e-6_rkx, &
    'check_domain sets the cone factor from the file')
  call check(all(mddom_io%xlat == r4(xlat)) .and. &
             all(mddom_io%xlon == r4(xlon)),'read_domain: coordinates')
  call check(all(mddom_io%ht == r4(topo)) .and. &
             all(mddom_io%mask == mask) .and. &
             all(mddom_io%lndcat == lu),'read_domain: surface fields')
  call check(all(mddom_io%dlat == r4(xlat)) .and. &
             all(mddom_io%msfx == 1.0_rkx),'read_domain: dot points')
  call check(all(mddom_io%xlat_bnds == r4(bnds)), &
    'read_domain: cell bounds')
  call check(all(abs(mddom_io%sigma-sigma) < 1.0e-12_rkx), &
    'read_domain: sigma')
  call check(all(mddom_io%hlake == 0.0_rkx), &
    'read_domain: missing lake depth is zero')
  call getmem(sig,1,nk+1,'sig')
  call getmem(rlat,1,nj,1,ni,'rlat')
  call getmem(rlon,1,nj,1,ni,'rlon')
  call getmem(rht,1,nj,1,ni,'rht')
  call getmem(rdlat,1,nj,1,ni,'rdlat')
  call getmem(rbnds,1,4,1,nj,1,ni,'rbnds')
  call read_domain(ncid,sig,rlat,rlon,dlat=rdlat,ht=rht, &
                   xlat_bnds=rbnds)
  call check(all(rlat == r4(xlat)) .and. all(rht == r4(topo)) .and. &
             all(rdlat == r4(xlat)) .and. all(rbnds == r4(bnds)), &
    'read_domain with optional arrays')
  ! A file of surface fields only needs the same vertical coordinate.
  idynamic = 2
  call check_domain(ncid,lsurface=.true.)
  call check(.true.,'surface file of the hydrostatic core, idynamic = 2')
  idynamic = 1
  call closefile(ncid)

  ! MOLOCH: staggered u and v points, no ptop.
  idynamic = 3
  call write_domain('test_domain_moloch.nc',.false.)
  call openfile_withname('test_domain_moloch.nc',ncid)
  call getmem(rulat,1,nj,1,ni,'rulat')
  call read_domain(ncid,sig,rlat,rlon,ulat=rulat,ht=rht)
  call check(all(rulat == r4(xlat)) .and. ptop == 0.0_rkx, &
    'MOLOCH domain')
  call closefile(ncid)

  ! Subgrid file.
  idynamic = 1
  ptop = 5.0_rkx
  call write_domain('test_domain_sub.nc',.true.)
  call openfile_withname('test_domain_sub.nc',ncid)
  call relmem(rlat)
  call relmem(rht)
  call getmem(rlat,1,3*nj,1,3*ni,'rlat')
  call getmem(rht,1,3*nj,1,3*ni,'rht')
  call read_domain(ncid,sig,rlat,ht=rht,lsubgrid=.true.)
  call check(all(rlat == r4(sxlat)) .and. all(rht == r4(stopo)), &
    'subgrid domain')
  call closefile(ncid)

  call memory_destroy
  call test_end('test_domain')

  contains

  ! The values as stored in the file (single precision).
  elemental real(rkx) function r4(x)
    implicit none
    real(rkx), intent(in) :: x
    r4 = real(real(x,rk4),rkx)
  end function r4

  subroutine write_domain(fname,lsub)
    implicit none
    character(len=*), intent(in) :: fname
    logical, intent(in) :: lsub
    type(nc_output_stream) :: ncout
    type(ncoutstream_params) :: opar
    type(ncvariable2d_mixed), dimension(12) :: v2
    type(ncvariable3d_mixed), dimension(4) :: vb
    character(len=9), dimension(4), parameter :: bnames = &
      [character(len=9) :: 'xlat_bnds', 'xlon_bnds', 'dlat_bnds', 'dlon_bnds']
    character(len=8), dimension(12) :: names
    character(len=8), dimension(4), parameter :: subnames = &
      [character(len=8) :: 'xlat', 'xlon', 'topo', 'mask']
    integer(ik4) :: n, nv
    if ( idynamic == 3 ) then
      names = [character(len=8) :: 'xlat', 'xlon', 'ulat', 'ulon', &
               'vlat', 'vlon', 'topo', 'mask', 'landuse', 'coriol', &
               'xmap', 'snowam']
    else
      names = [character(len=8) :: 'xlat', 'xlon', 'dlat', 'dlon', &
               'xmap', 'dmap', 'topo', 'mask', 'landuse', 'coriol', &
               'snowam', 'areacella']
    end if
    nv = 12
    if ( lsub ) nv = 4
    opar%fname = fname
    opar%pname = 'terrain'
    opar%l_bound = .true.
    opar%l_subgrid = lsub
    opar%l_full_sigma = .true.
    call outstream_setup(ncout,opar)
    do n = 1, nv
      if ( lsub ) then
        v2(n)%vname = subnames(n)
      else
        v2(n)%vname = names(n)
      end if
      v2(n)%vunit = '1'
      v2(n)%long_name = v2(n)%vname
      v2(n)%standard_name = ''
      v2(n)%lrecords = .false.
      call outstream_addvar(ncout,v2(n))
    end do
    if ( .not. lsub ) then
      do n = 1, 4
        vb(n)%vname = bnames(n)
        vb(n)%vunit = 'degrees'
        vb(n)%long_name = ''
        vb(n)%standard_name = ''
        vb(n)%axis = 'nxy'
        vb(n)%lrecords = .false.
        call outstream_addvar(ncout,vb(n))
      end do
    end if
    call outstream_enable(ncout,sigma)
    do n = 1, nv
      if ( lsub ) then
        select case (n)
          case (1)
            sl => sxlat
          case (2)
            sl => sxlon
          case (3)
            sl => stopo
          case default
            sl => smask
        end select
      else
        select case (trim(names(n)))
          case ('xlat','dlat','ulat','vlat')
            sl => xlat
          case ('xlon','dlon','ulon','vlon')
            sl => xlon
          case ('topo')
            sl => topo
          case ('mask')
            sl => mask
          case ('landuse')
            sl => lu
          case ('coriol')
            sl => cor
          case ('xmap','dmap')
            sl => ones
          case default
            sl => cor
        end select
      end if
      v2(n)%rval => sl
      call outstream_writevar(ncout,v2(n))
    end do
    if ( .not. lsub ) then
      do n = 1, 4
        vb(n)%rval => bnds
        call outstream_writevar(ncout,vb(n))
      end do
    end if
    call outstream_dispose(ncout)
    call outstream_async_finalize( )
  end subroutine write_domain

end program test_domain

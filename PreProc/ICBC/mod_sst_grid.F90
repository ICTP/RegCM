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

module mod_sst_grid

  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_dynparam
  use mod_memutil
  use mod_message
  use mod_nchelper
  use mod_ncstream_types
  use mod_ncstream
  use mod_domain
  implicit none

  private

  real(rkx), public, pointer, contiguous, dimension(:,:) :: sstmm, icemm
  real(rkx), public, pointer, contiguous, dimension(:,:) :: xlat, xlon
  real(rkx), pointer, contiguous, dimension(:,:,:) :: xlat_bnds, xlon_bnds
  real(rkx), pointer, contiguous, dimension(:,:) :: topo, mask, lndcat
  real(rkx), pointer, contiguous, dimension(:,:) :: sstsrc
  real(rkx), pointer, contiguous, dimension(:) :: sigma

  type(nc_output_stream), save :: ncout
  integer(ik4), parameter :: nvar2d = 4
  integer(ik4), parameter :: nvar_bnds = 2   ! xlon_bnds, xlat_bnds
  type(ncvariable2d_mixed), save, dimension(nvar2d) :: v2dvar_base
  type(ncvariable3d_mixed), save, dimension(nvar_bnds) :: v3dvar_bnds
  type(ncvariable2d_mixed), save :: v2dvar_sst

  public :: init_grid, read_domain_info, setup_outvars, open_sstfile, &
            close_sstfile, writerec

  contains

  subroutine init_grid
    implicit none
    call getmem(sstmm,1,jx,1,iy,'mod_sst_grid:sstmm')
    call getmem(icemm,1,jx,1,iy,'mod_sst_grid:icemm')
    call getmem(xlat,1,jx,1,iy,'mod_sst_grid:xlat')
    call getmem(xlon,1,jx,1,iy,'mod_sst_grid:xlon')
    call getmem(xlat_bnds,1,4,1,jx,1,iy,'mod_sst_grid:xlat_bnds')
    call getmem(xlon_bnds,1,4,1,jx,1,iy,'mod_sst_grid:xlon_bnds')
    call getmem(topo,1,jx,1,iy,'mod_sst_grid:topo')
    call getmem(mask,1,jx,1,iy,'mod_sst_grid:mask')
    call getmem(sstsrc,1,jx,1,iy,'mod_sst_grid:sstsrc')
    call getmem(sigma,1,kzp1,'mod_sst_grid:sigma')
  end subroutine init_grid

  subroutine read_domain_info(terfile)
    implicit none
    character(len=256), intent(in) :: terfile
    integer(ik4) :: incin
    call openfile_withname(terfile,incin)
    ! Use the type reader so that mddom_io (including bound arrays) is fully
    ! populated; then alias the fields we need into local pointers.
    call read_domain(incin)
    sigma    => mddom_io%sigma
    xlat     => mddom_io%xlat
    xlon     => mddom_io%xlon
    xlat_bnds => mddom_io%xlat_bnds
    xlon_bnds => mddom_io%xlon_bnds
    topo     => mddom_io%ht
    mask     => mddom_io%mask
    lndcat   => mddom_io%lndcat
    call closefile(incin)
  end subroutine read_domain_info

  subroutine setup_outvars
    implicit none
    v2dvar_base(1)%vname = 'xlon'
    v2dvar_base(1)%vunit = 'degrees_east'
    v2dvar_base(1)%long_name = 'Longitude on Cross Points'
    v2dvar_base(1)%standard_name = 'longitude'
    v2dvar_base(2)%vname = 'xlat'
    v2dvar_base(2)%vunit = 'degrees_north'
    v2dvar_base(2)%long_name = 'Latitude on Cross Points'
    v2dvar_base(2)%standard_name = 'latitude'
    v2dvar_base(3)%vname = 'mask'
    v2dvar_base(3)%vunit = '1'
    v2dvar_base(3)%long_name = 'Land Mask'
    v2dvar_base(3)%standard_name = 'land_binary_mask'
    v2dvar_base(4)%vname = 'topo'
    v2dvar_base(4)%vunit = 'm'
    v2dvar_base(4)%long_name = 'Surface Model Elevation'
    v2dvar_base(4)%standard_name = 'surface_altitude'

    v3dvar_bnds(1)%vname = 'xlon_bnds'
    v3dvar_bnds(1)%vunit = 'degrees_east'
    v3dvar_bnds(1)%long_name = 'Cell corner longitudes on Cross Points'
    v3dvar_bnds(1)%standard_name = 'longitude'
    v3dvar_bnds(1)%axis = 'nxy'
    v3dvar_bnds(2)%vname = 'xlat_bnds'
    v3dvar_bnds(2)%vunit = 'degrees_north'
    v3dvar_bnds(2)%long_name = 'Cell corner latitudes on Cross Points'
    v3dvar_bnds(2)%standard_name = 'latitude'
    v3dvar_bnds(2)%axis = 'nxy'
    v3dvar_bnds(1)%j1 = -1
    v3dvar_bnds(1)%j2 = -1
    v3dvar_bnds(1)%i1 = -1
    v3dvar_bnds(1)%i2 = -1
    v3dvar_bnds(1)%k1 = -1
    v3dvar_bnds(1)%k2 = -1
    v3dvar_bnds(2)%j1 = -1
    v3dvar_bnds(2)%j2 = -1
    v3dvar_bnds(2)%i1 = -1
    v3dvar_bnds(2)%i2 = -1
    v3dvar_bnds(2)%k1 = -1
    v3dvar_bnds(2)%k2 = -1

    v2dvar_sst%vname = 'sst'
    v2dvar_sst%vunit = 'K'
    v2dvar_sst%long_name = 'Sea Surface Temperature'
    v2dvar_sst%standard_name = 'sea_surface_temperature'
    v2dvar_sst%lfillvalue = .true.
    v2dvar_sst%rmissval = -9999.0
    v2dvar_sst%lrecords = .true.
  end subroutine setup_outvars

  subroutine open_sstfile(idate1)
    implicit none
    type(rcm_time_and_date), intent(in) :: idate1
    character(len=256) :: sstname
    type(ncoutstream_params) :: opar
    integer(ik4) :: ivar
    sstname = trim(dirglob)//pthsep//trim(domname)//'_SST.nc'
    opar%fname = sstname
    opar%pname = 'sst'
    opar%zero_date = idate1
    opar%l_bound = .true.
    opar%l_full_sigma = .true.
    call outstream_setup(ncout,opar)
    call outstream_addatt(ncout,ncattribute_string('sst_source',ssttyp))
    v2dvar_base(1)%rval => xlon
    v2dvar_base(2)%rval => xlat
    v2dvar_base(3)%rval => mask
    v2dvar_base(4)%rval => topo
    v3dvar_bnds(1)%rval => xlon_bnds
    v3dvar_bnds(2)%rval => xlat_bnds
    do ivar = 1, nvar2d
      call outstream_addvar(ncout,v2dvar_base(ivar))
    end do
    do ivar = 1, nvar_bnds
      call outstream_addvar(ncout,v3dvar_bnds(ivar))
    end do
    ! CF 'bounds' attributes on coordinate variables
    call outstream_addvaratt(ncout,v2dvar_base(1), &
      ncattribute_string('bounds','xlon_bnds'))
    call outstream_addvaratt(ncout,v2dvar_base(2), &
      ncattribute_string('bounds','xlat_bnds'))
    call outstream_addvar(ncout,v2dvar_sst)
    v2dvar_sst%rval => sstmm
    call outstream_enable(ncout,sigma)
    do ivar = 1, nvar2d
      call outstream_writevar(ncout,v2dvar_base(ivar))
    end do
    do ivar = 1, nvar_bnds
      call outstream_writevar(ncout,v3dvar_bnds(ivar))
    end do
  end subroutine open_sstfile

  subroutine close_sstfile
    implicit none
    call outstream_dispose(ncout)
  end subroutine close_sstfile

  subroutine writerec(idate)
    implicit none
    type(rcm_time_and_date), intent(in) :: idate
    call fill_ocean_sst(idate)
    call outstream_addrec(ncout,idate)
    call outstream_writevar(ncout,v2dvar_sst)
  end subroutine writerec
  !
  ! Ocean points (landuse 15) left without SST by the interpolation get
  ! the inverse distance weighted mean of the valid points in the first
  ! ring around them holding any, floored at the freezing point of sea
  ! water. Only values from the source data are used, never filled ones,
  ! also from land points. All points other than ocean are then set
  ! missing, so the file shows where icbc uses the SST.
  !
  subroutine fill_ocean_sst(idate)
    implicit none
    type(rcm_time_and_date), intent(in) :: idate
    integer(ik4) :: i, j, ii, jj, nr, maxn, nfill, nleft
    real(rkx) :: icet, wt, wtsum, ssum
    logical :: found, lany

    select case (ssttyp)
      case ('EIN15','EIN75','EIXXX')
        icet = 271.46_rkx
      case default
        icet = 271.35_rkx
    end select
    sstsrc(:,:) = sstmm(:,:)
    lany = any(sstsrc > -900.0_rkx)
    maxn = max(min(jx,iy)/4,1)
    nfill = 0
    nleft = 0
    do i = 1, iy
      do j = 1, jx
        if ( lndcat(j,i) < 14.9_rkx .or. lndcat(j,i) > 15.1_rkx ) cycle
        if ( sstsrc(j,i) > -900.0_rkx ) cycle
        found = .false.
        do nr = 1, merge(maxn,0,lany)
          ssum = 0.0_rkx
          wtsum = 0.0_rkx
          do ii = max(i-nr,1), min(i+nr,iy)
            do jj = max(j-nr,1), min(j+nr,jx)
              if ( max(abs(ii-i),abs(jj-j)) /= nr ) cycle
              if ( sstsrc(jj,ii) > -900.0_rkx ) then
                wt = 1.0_rkx/sqrt(real(ii-i,rkx)**2+real(jj-j,rkx)**2)
                ssum = ssum + max(icet,sstsrc(jj,ii))*wt
                wtsum = wtsum + wt
              end if
            end do
          end do
          if ( wtsum > 0.0_rkx ) then
            sstmm(j,i) = ssum/wtsum
            found = .true.
            exit
          end if
        end do
        if ( found ) then
          nfill = nfill + 1
        else
          nleft = nleft + 1
        end if
      end do
    end do
    ! Only ocean points receive the SST in icbc: mark all others missing
    where ( lndcat < 14.9_rkx .or. lndcat > 15.1_rkx )
      sstmm = -9999.0_rkx
    end where
    if ( nfill > 0 ) then
      write(stdout,'(a,i8,a,a)') ' Filled SST at ', nfill, &
        ' ocean points for ', tochar(idate)
    end if
    if ( nleft > 0 ) then
      write(stderr,'(a,i8,a,a)') ' WARNING: no SST at ', nleft, &
        ' ocean points for ', tochar(idate)
    end if
  end subroutine fill_ocean_sst

end module mod_sst_grid
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

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

module mod_rdldtr

  use mod_intkinds
  use mod_realkinds
  use mod_stdio
  use mod_constants
  use mod_memutil
  use mod_nchelper
  use mod_message
  use mod_earth
  use mod_kdinterp

  private

  public :: globalfile , gfopen , gfclose , gfread , bestaround

  public :: read_ncglob

  real(rk8) , allocatable , dimension(:) :: glat
  real(rk8) , allocatable , dimension(:) :: glon

  type globalfile
    integer(ik4) :: ncid
    type(global_domain) :: gdomain
    type(h_interpolator) :: hint
    logical :: lmask = .false.
    integer(ik4) , dimension(:,:) , allocatable :: mask
  end type globalfile

  integer(ik4) :: ncid , istatus
  type(global_domain) :: gdomain

  interface read_ncglob
    module procedure read_ncglob2d
    module procedure read_ncglob3d
    module procedure read_ncglob2d3d
  end interface read_ncglob

  interface gfread
    module procedure gfread_2d
    module procedure gfread_2di
    module procedure gfread_2d_landuse
    module procedure gfread_2d3d
    module procedure gfread_3d
    module procedure gfread_3d_lookup
    module procedure gfread_4d
    module procedure gfread_4d_lookup
    module procedure gfread_5d
    module procedure gfread_5d_lookup
  end interface gfread

  interface bestaround
    module procedure bestaround2d
    module procedure bestaround3d
  end interface

  contains
!
!
! Read a netcdf file with a global dataset on a regular latitude/longitude grid
!
!  cfile = name of the input file
!  cvar = name of the var in the file
!  xlat , xlon = grid latitudes and longitudes of requested grid
!  iband = band flag
!  iores = output wanted resolution in minutes
!          iores == 0 then no resampling is performed
!  imeth = flag for selecting resampling method
!     default = nearest
!          1  = mean
!          2  = median
!          3  = most present
!          4  = mean of 50% central median filtered
!  In output :
!  grdlnma,grdlnmn : Minimum and maximum in longitude of the read grid
!  grdltma,grdltmn : Minimum and maximum in latitude of the read grid
!  nlatin,nlonin   : Size of the read grid
!  values : allocated space by the sub containing data: the caller is
!           in charge of the deallocate
!
  subroutine read_ncglob2d(cfile,cvar,iores,imeth,iband,xlat,xlon, &
                           grdlnma,grdlnmn,grdltma,grdltmn,        &
                           nlatin,nlonin,values)
    use netcdf
    implicit none
    character(len=*) , intent(in) :: cfile , cvar
    integer(ik4) , intent(in) :: iband , iores , imeth
    real(rkx) , dimension(:,:) , intent(in) :: xlat , xlon
    real(rk8) , intent(out) :: grdlnma , grdlnmn , grdltma , grdltmn
    integer(ik4) , intent(out) :: nlatin , nlonin
    real(rkx) , dimension(:,:) , intent(inout) , pointer :: values
    integer(ik4) :: nlat , nlon , iti , itf , itile , ivar
    integer(ik4) :: i , j , inpsec , iopsec , ifrac
    integer(ik4) , dimension(2) :: istart , icount
    real(rk8) :: deltalat , deltalon
    real(rkx) , allocatable , dimension(:,:) :: readbuf

#ifdef DEBUG
    write(stdout,*) 'Opening '//trim(cfile)
#endif
    istatus = nf90_open(cfile, nf90_nowrite, ncid)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error: '//trim(cfile))

    call read_geolocation(cfile)

    istatus = nf90_inq_varid(ncid, cvar, ivar)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')

    call get_window(glat,glon,xlat,xlon,iband,gdomain)

    nlat = gdomain%nj
    nlon = sum(gdomain%ni)

    inpsec = int(abs(glat(2)-glat(1))*3600.0_rk8)
    iopsec = max(int(real(iores,rk8)*60.0_rk8),inpsec)
    ifrac = max(iopsec/inpsec,1)
#ifdef DEBUG
    write(stderr,*) 'INPSEC = ', inpsec
    write(stderr,*) 'IOPSEC = ', iopsec
    write(stderr,*) 'IFRAC  = ', ifrac
#endif

    if ( ifrac > 1 ) then
      nlatin = (nlat/ifrac) - 1
      nlonin = (nlon/ifrac) - 1
      deltalat = (real(iopsec,rk8)/2.0_rk8) / 3600.0_rk8
      deltalon = (real(iopsec,rk8)/2.0_rk8) / 3600.0_rk8
    else
      nlatin = nlat
      nlonin = nlon
      deltalat = 0.0_rk8
      deltalon = 0.0_rk8
    end if

    grdlnmn = glon(gdomain%igstart(1)) + deltalon
    grdlnma = grdlnmn + real((nlonin-1)*iopsec,rk8) / 3600.0_rk8
    grdltmn = glat(gdomain%jgstart) + deltalat
    grdltma = grdltmn + real((nlatin-1)*iopsec,rk8) / 3600.0_rk8

    deallocate(glat)
    deallocate(glon)

    call getmem2d(values,1,nlonin,1,nlatin,'rdldtr:values')
    allocate(readbuf(nlon,nlat))

#ifdef DEBUG
    write(stderr,*) 'BOUNDS IN LAT: ', grdltmn , grdltma
    write(stderr,*) 'BOUNDS IN LON: ', grdlnmn , grdlnma
#endif

#ifdef DEBUG
    write(stderr,*) 'WILL READ ', nlon , 'x', nlat, ' points'
    write(stderr,*) 'WILL GIVE ', nlonin , 'x', nlatin, ' points'
#endif

    iti = 1
    readbuf = -1000000000
    do itile = 1 , gdomain%ntiles
      istart(1) = gdomain%igstart(itile)
      icount(1) = gdomain%ni(itile)
      istart(2) = gdomain%jgstart
      icount(2) = gdomain%nj
      itf = iti + gdomain%ni(itile) - 1
      istatus = nf90_get_var(ncid,ivar,readbuf(iti:itf,:),istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')
      iti = itf + 1
    end do

    istatus = nf90_close(ncid)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')

    ! Fix Poles for interpolations
    if ( istart(2) == 1 ) then
#ifdef DEBUG
      write (stdout,*) 'Correcting South pole.'
#endif
      readbuf(:,1) = readbuf(:,2)
    else if ( istart(2)+icount(2)-1 == gdomain%global_nj ) then
#ifdef DEBUG
      write (stdout,*) 'Correcting North pole.'
#endif
      readbuf(:,nlat) = readbuf(:,nlat-1)
    end if

    if ( ifrac > 1 ) then
      call resampling(ifrac,imeth,iband,nlat,nlon,nlatin,nlonin,readbuf,values)
    else
      do i = 1 , nlatin
        do j = 1 , nlonin
          values(j,i) = readbuf(j,i)
        end do
      end do
    end if
    deallocate(readbuf)
#ifdef DEBUG
    write(stdout,'(a)') ' Done.'
#endif
  end subroutine read_ncglob2d

  subroutine read_ncglob3d(cfile,cvar,iores,imeth,iband,xlat,xlon, &
                           grdlnma,grdlnmn,grdltma,grdltmn,        &
                           nlatin,nlonin,values)
    use netcdf
    implicit none
    character(len=*) , intent(in) :: cfile , cvar
    integer(ik4) , intent(in) :: iband , iores , imeth
    real(rkx) , dimension(:,:) , intent(in) :: xlat , xlon
    real(rk8) , intent(out) :: grdlnma , grdlnmn , grdltma , grdltmn
    integer(ik4) , intent(out) :: nlatin , nlonin
    real(rkx) , dimension(:,:,:) , intent(inout) , pointer :: values
    integer(ik4) :: nlat , nlon , iti , itf , itile , ivar
    integer(ik4) :: i , j , n , inpsec , iopsec , ifrac , nd
    integer(ik4) , dimension(3) :: idims , istart , icount
    real(rk8) :: deltalat , deltalon
    real(rkx) , allocatable , dimension(:,:,:) :: readbuf

#ifdef DEBUG
    write(stdout,*) 'Opening '//trim(cfile)
#endif
    istatus = nf90_open(cfile, nf90_nowrite, ncid)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error: '//trim(cfile))

    call read_geolocation(cfile)

    istatus = nf90_inq_varid(ncid, cvar, ivar)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')

    istatus = nf90_inquire_variable(ncid,ivar,ndims=nd)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')
    if ( nd /= 3 ) then
      call die('rdldtr','Variable '//trim(cvar)// &
                        ' is not 3d in file '//trim(cfile),1)
    end if
    istatus = nf90_inquire_variable(ncid,ivar,dimids=idims)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')

    istatus = nf90_inquire_dimension(ncid,idims(3),len=nd)
    call checkncerr(istatus,__FILE__,__LINE__, &
         'Error reading lon dimelen in file '//trim(cfile))

    call get_window(glat,glon,xlat,xlon,iband,gdomain)

    nlat = gdomain%nj
    nlon = sum(gdomain%ni)

    inpsec = int(abs(glat(2)-glat(1))*3600.0_rk8)
    iopsec = max(int(real(iores,rk8)*60.0_rk8),inpsec)
    ifrac = max(iopsec/inpsec,1)
#ifdef DEBUG
    write(stderr,*) 'INPSEC = ', inpsec
    write(stderr,*) 'IOPSEC = ', iopsec
    write(stderr,*) 'IFRAC  = ', ifrac
#endif

    if ( ifrac > 1 ) then
      nlatin = (nlat/ifrac) - 1
      nlonin = (nlon/ifrac) - 1
      deltalat = (real(iopsec,rk8)/2.0_rk8) / 3600.0_rk8
      deltalon = (real(iopsec,rk8)/2.0_rk8) / 3600.0_rk8
    else
      nlatin = nlat
      nlonin = nlon
      deltalat = 0.0_rk8
      deltalon = 0.0_rk8
    end if

    grdlnmn = glon(gdomain%igstart(1)) + deltalon
    grdlnma = grdlnmn + real((nlonin-1)*iopsec,rk8) / 3600.0_rk8
    grdltmn = glat(gdomain%jgstart) + deltalat
    grdltma = grdltmn + real((nlatin-1)*iopsec,rk8) / 3600.0_rk8

    deallocate(glat)
    deallocate(glon)

    call getmem3d(values,1,nlonin,1,nlatin,1,nd,'rdldtr:values')
    allocate(readbuf(nlon,nlat,nd))

#ifdef DEBUG
    write(stderr,*) 'BOUNDS IN LAT: ', grdltmn , grdltma
    write(stderr,*) 'BOUNDS IN LON: ', grdlnmn , grdlnma
#endif

#ifdef DEBUG
    write(stderr,*) 'WILL READ ', nlon , 'x', nlat, ' points'
    write(stderr,*) 'WILL GIVE ', nlonin , 'x', nlatin, ' points'
#endif

    iti = 1
    readbuf = -1000000000
    do itile = 1 , gdomain%ntiles
      istart(1) = gdomain%igstart(itile)
      icount(1) = gdomain%ni(itile)
      istart(2) = gdomain%jgstart
      icount(2) = gdomain%nj
      istart(3) = 1
      icount(3) = nd
      itf = iti + gdomain%ni(itile) - 1
      istatus = nf90_get_var(ncid,ivar,readbuf(iti:itf,:,:),istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')
      iti = itf + 1
    end do

    istatus = nf90_close(ncid)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')

    ! Fix Poles for interpolations
    if ( istart(2) == 1 ) then
#ifdef DEBUG
      write (stdout,*) 'Correcting South pole.'
#endif
      readbuf(:,1,:) = readbuf(:,2,:)
    else if ( istart(2)+icount(2)-1 == gdomain%global_nj ) then
#ifdef DEBUG
      write (stdout,*) 'Correcting North pole.'
#endif
      readbuf(:,nlat,:) = readbuf(:,nlat-1,:)
    end if

    if ( ifrac > 1 ) then
      do n = 1 , nd
        call resampling(ifrac,imeth,iband,nlat,nlon,nlatin,nlonin, &
                        readbuf(:,:,n),values(:,:,n))
      end do
    else
      do n = 1 , nd
        do i = 1 , nlatin
          do j = 1 , nlonin
            values(j,i,n) = readbuf(j,i,n)
          end do
        end do
      end do
    end if
    deallocate(readbuf)
#ifdef DEBUG
    write(stdout,'(a)') ' Done.'
#endif
  end subroutine read_ncglob3d

  subroutine read_ncglob2d3d(cfile,cvar,iores,imeth,iband,xlat,xlon, &
                             grdlnma,grdlnmn,grdltma,grdltmn,        &
                             nlatin,nlonin,values,isel)
    use netcdf
    implicit none
    character(len=*) , intent(in) :: cfile , cvar
    integer(ik4) , intent(in) :: iband , iores , imeth
    real(rkx) , dimension(:,:) , intent(in) :: xlat , xlon
    real(rk8) , intent(out) :: grdlnma , grdlnmn , grdltma , grdltmn
    integer(ik4) , intent(out) :: nlatin , nlonin
    real(rkx) , dimension(:,:) , intent(inout) , pointer :: values
    integer(ik4) , intent(in) :: isel
    integer(ik4) :: nlat , nlon , iti , itf , itile , ivar
    integer(ik4) :: i , j , inpsec , iopsec , ifrac , nd
    integer(ik4) , dimension(3) :: istart , icount
    real(rk8) :: deltalat , deltalon
    real(rkx) , allocatable , dimension(:,:) :: readbuf

#ifdef DEBUG
    write(stdout,*) 'Opening '//trim(cfile)
#endif
    istatus = nf90_open(cfile, nf90_nowrite, ncid)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error: '//trim(cfile))

    call read_geolocation(cfile)

    istatus = nf90_inq_varid(ncid, cvar, ivar)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')

    istatus = nf90_inquire_variable(ncid,ivar,ndims=nd)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')
    if ( nd /= 3 ) then
      call die('rdldtr','Variable '//trim(cvar)// &
                        ' is not 3d in file '//trim(cfile),1)
    end if

    call get_window(glat,glon,xlat,xlon,iband,gdomain)

    nlat = gdomain%nj
    nlon = sum(gdomain%ni)

    inpsec = int(abs(glat(2)-glat(1))*3600.0_rk8)
    iopsec = max(int(real(iores,rk8)*60.0_rk8),inpsec)
    ifrac = max(iopsec/inpsec,1)
#ifdef DEBUG
    write(stderr,*) 'INPSEC = ', inpsec
    write(stderr,*) 'IOPSEC = ', iopsec
    write(stderr,*) 'IFRAC  = ', ifrac
#endif

    if ( ifrac > 1 ) then
      nlatin = (nlat/ifrac) - 1
      nlonin = (nlon/ifrac) - 1
      deltalat = (real(iopsec,rk8)/2.0_rk8) / 3600.0_rk8
      deltalon = (real(iopsec,rk8)/2.0_rk8) / 3600.0_rk8
    else
      nlatin = nlat
      nlonin = nlon
      deltalat = 0.0_rk8
      deltalon = 0.0_rk8
    end if

    grdlnmn = glon(gdomain%igstart(1)) + deltalon
    grdlnma = grdlnmn + real((nlonin-1)*iopsec,rk8) / 3600.0_rk8
    grdltmn = glat(gdomain%jgstart) + deltalat
    grdltma = grdltmn + real((nlatin-1)*iopsec,rk8) / 3600.0_rk8

    deallocate(glat)
    deallocate(glon)

    call getmem2d(values,1,nlonin,1,nlatin,'rdldtr:values')
    allocate(readbuf(nlon,nlat))

#ifdef DEBUG
    write(stderr,*) 'BOUNDS IN LAT: ', grdltmn , grdltma
    write(stderr,*) 'BOUNDS IN LON: ', grdlnmn , grdlnma
#endif

#ifdef DEBUG
    write(stderr,*) 'WILL READ ', nlon , 'x', nlat, ' points'
    write(stderr,*) 'WILL GIVE ', nlonin , 'x', nlatin, ' points'
#endif

    iti = 1
    readbuf = -1000000000
    do itile = 1 , gdomain%ntiles
      istart(1) = gdomain%igstart(itile)
      icount(1) = gdomain%ni(itile)
      istart(2) = gdomain%jgstart
      icount(2) = gdomain%nj
      istart(3) = isel
      icount(3) = 1
      itf = iti + gdomain%ni(itile) - 1
      istatus = nf90_get_var(ncid,ivar,readbuf(iti:itf,:),istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')
      iti = itf + 1
    end do

    istatus = nf90_close(ncid)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')

    ! Fix Poles for interpolations
    if ( istart(2) == 1 ) then
#ifdef DEBUG
      write (stdout,*) 'Correcting South pole.'
#endif
      readbuf(:,1) = readbuf(:,2)
    else if ( istart(2)+icount(2)-1 == gdomain%global_nj ) then
#ifdef DEBUG
      write (stdout,*) 'Correcting North pole.'
#endif
      readbuf(:,nlat) = readbuf(:,nlat-1)
    end if

    if ( ifrac > 1 ) then
      call resampling(ifrac,imeth,iband,nlat,nlon,nlatin,nlonin,readbuf,values)
    else
      do i = 1 , nlatin
        do j = 1 , nlonin
          values(j,i) = readbuf(j,i)
        end do
      end do
    end if
    deallocate(readbuf)
#ifdef DEBUG
    write(stdout,'(a)') ' Done.'
#endif
  end subroutine read_ncglob2d3d

  subroutine read_geolocation(cfile)
    use netcdf
    implicit none
    character(len=*) , intent(in) :: cfile
    integer(ik4) :: idimid , idvar
    integer(ik4) :: jlat , ilon

    istatus = nf90_inq_dimid(ncid,'lat',idimid)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_dimid(ncid,'latitude',idimid)
      if ( istatus /= nf90_noerr ) then
        istatus = nf90_inq_dimid(ncid,'LAT',idimid)
        if ( istatus /= nf90_noerr ) then
          istatus = nf90_inq_dimid(ncid,'LATITUDE',idimid)
          if ( istatus /= nf90_noerr ) then
            istatus = nf90_inq_dimid(ncid,'lsmlat',idimid)
          end if
        end if
      end if
    end if
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Missing lat dimension in file '//trim(cfile))
    istatus = nf90_inquire_dimension(ncid,idimid,len=jlat)
    call checkncerr(istatus,__FILE__,__LINE__, &
         'Error reading lat dimelen in file '//trim(cfile))

    istatus = nf90_inq_dimid(ncid,'lon',idimid)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_dimid(ncid,'longitude',idimid)
      if ( istatus /= nf90_noerr ) then
        istatus = nf90_inq_dimid(ncid,'LON',idimid)
        if ( istatus /= nf90_noerr ) then
          istatus = nf90_inq_dimid(ncid,'LONGITUDE',idimid)
          if ( istatus /= nf90_noerr ) then
            istatus = nf90_inq_dimid(ncid,'lsmlon',idimid)
          end if
        end if
      end if
    end if
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Missing lon dimension in file '//trim(cfile))
    istatus = nf90_inquire_dimension(ncid,idimid,len=ilon)
    call checkncerr(istatus,__FILE__,__LINE__, &
         'Error reading lon dimelen in file '//trim(cfile))

    allocate(glat(jlat))
    allocate(glon(ilon))

    istatus = nf90_inq_varid(ncid,'lat',idvar)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_varid(ncid,'latitude',idvar)
      if ( istatus /= nf90_noerr ) then
        istatus = nf90_inq_varid(ncid,'LAT',idvar)
        if ( istatus /= nf90_noerr ) then
          istatus = nf90_inq_varid(ncid,'LATITUDE',idvar)
          if ( istatus /= nf90_noerr ) then
            istatus = nf90_inq_varid(ncid,'lsmlat',idvar)
          end if
        end if
      end if
    end if
    call checkncerr(istatus,__FILE__,__LINE__, &
        'Missing lat variable in file '//trim(cfile))
    istatus = nf90_get_var(ncid,idvar,glat)
    call checkncerr(istatus,__FILE__,__LINE__, &
        'Error reading lat variable in file '//trim(cfile))

    istatus = nf90_inq_varid(ncid,'lon',idvar)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_varid(ncid,'longitude',idvar)
      if ( istatus /= nf90_noerr ) then
        istatus = nf90_inq_varid(ncid,'LON',idvar)
        if ( istatus /= nf90_noerr ) then
          istatus = nf90_inq_varid(ncid,'LONGITUDE',idvar)
          if ( istatus /= nf90_noerr ) then
            istatus = nf90_inq_varid(ncid,'lsmlon',idvar)
          end if
        end if
      end if
    end if
    call checkncerr(istatus,__FILE__,__LINE__, &
        'Missing lon variable in file '//trim(cfile))
    istatus = nf90_get_var(ncid,idvar,glon)
    call checkncerr(istatus,__FILE__,__LINE__, &
        'Error reading lon variable in file '//trim(cfile))
  end subroutine read_geolocation

  subroutine resampling(ifrac,imeth,iband,nlat,nlon,nlatin,nlonin, &
                        readbuf,values)
    implicit none
    integer(ik4) , intent(in) :: ifrac , imeth , iband
    integer(ik4) , intent(in) :: nlat , nlon
    integer(ik4) , intent(in) :: nlatin , nlonin
    real(rkx) , dimension(nlon,nlat) , intent(in) :: readbuf
    real(rkx) , dimension(nlonin,nlatin) , intent(out) :: values
    integer(ik4) :: nfrac , np , i , j , ib , jb , iprint
    real(rkx) , allocatable , dimension(:) :: copybuf

    write(stdout,'(a)',advance='no') ' Resampling'
    iprint = max(nlatin/20,5)
    nfrac = (ifrac+1)*(ifrac+1)
    allocate(copybuf(nfrac))
    select case (imeth)
      case (1)
        do i = 1 , nlatin
          if (mod(i,iprint) == 0) write(stdout,'(a)',advance='no') '.'
          ib = (i-1)*ifrac+1
          do j = 1 , nlonin
            jb = (j-1)*ifrac+1
            call fillbuf(copybuf,readbuf,nlon,nlat,jb,ib,ifrac+1,iband)
            values(j,i) = sum(copybuf)/real(size(copybuf),rkx)
          end do
        end do
      case (2)
        do i = 1 , nlatin
          if (mod(i,iprint) == 0) write(stdout,'(a)',advance='no') '.'
          ib = (i-1)*ifrac+1
          do j = 1 , nlonin
            jb = (j-1)*ifrac+1
            call fillbuf(copybuf,readbuf,nlon,nlat,jb,ib,ifrac+1,iband)
            call qsort(copybuf)
            !values(j,i) = 0.5*(copybuf(nfrac/2)+copybuf(nfrac/2+1))
            values(j,i) = copybuf(max(nfrac/2,1))
          end do
        end do
      case (3)
        do i = 1 , nlatin
          if (mod(i,iprint) == 0) write(stdout,'(a)',advance='no') '.'
          ib = (i-1)*ifrac+1
          do j = 1 , nlonin
            jb = (j-1)*ifrac+1
            call fillbuf(copybuf,readbuf,nlon,nlat,jb,ib,ifrac+1,iband)
            values(j,i) = real(mpindex(copybuf),rkx)
          end do
        end do
      case (4)
        do i = 1 , nlatin
          if (mod(i,iprint) == 0) write(stdout,'(a)',advance='no') '.'
          ib = (i-1)*ifrac+1
          do j = 1 , nlonin
            jb = (j-1)*ifrac+1
            call fillbuf(copybuf,readbuf,nlon,nlat,jb,ib,ifrac+1,iband)
            call qsort(copybuf)
            np = (ifrac*ifrac)/4
            values(j,i) = sum(copybuf(1+np:size(copybuf)-np+1)) / &
                           real(size(copybuf)-2*np,rkx)
          end do
        end do
      case default
        do i = 1 , nlatin
        ib = (i-1)*ifrac+1 + ifrac/2
          if (mod(i,iprint) == 0) write(stdout,'(a)',advance='no') '.'
          do j = 1 , nlonin
            jb = (j-1)*ifrac+1 + ifrac/2
            values(j,i) = readbuf(jb,ib)
          end do
        end do
    end select
    deallocate(copybuf)
    write(stdout,'(a)',advance='no') new_line('a')
  end subroutine resampling

  subroutine fillbuf(copybuf,readbuf,ni,nj,i,j,isize,iband)
    implicit none
    integer(ik4) , intent(in) :: ni , nj , isize , iband
    real(rkx) , dimension(isize*isize) , intent(out) :: copybuf
    real(rkx) , dimension(ni,nj) , intent(in) :: readbuf
    integer(ik4) , intent(in) :: i , j
    integer(ik4) :: imin , imax , jmin , jmax , icnt , jcnt , ip
    integer(ik4) :: ib , jb

    imin = i
    imax = i + isize
    jmin = j
    jmax = j + isize

    ip = 1
    do icnt = 1 , isize
      do jcnt = 1 , isize
        ib = imin+icnt-1
        jb = jmin+jcnt-1
        if ( iband == 1 ) then
          if (ib < 1) then
            ib = ni - ib
          else if (ib > ni) then
            ib = ib - ni
          end if
        end if
        ib = min(max(ib,1),ni)
        jb = min(max(jb,1),nj)
        copybuf(ip) = readbuf(ib,jb)
        ip = ip + 1
      end do
    end do
  end subroutine fillbuf

  pure integer(ik4) function mpindex(x) result(res)
    implicit none
    real(rkx) , dimension(:) , intent(in) :: x
    integer(ik4) , dimension(32) :: cnt
    integer(ik4) :: i
    cnt = 0
    do i = 1 , 32
      cnt(i) = count(int(x) == i)
    end do
    res = maxloc(cnt,1,cnt>0)
  end function mpindex

  recursive subroutine qsort(a)
    implicit none
    real(rkx) , dimension(:) , intent(in out) :: a
    integer(ik4) :: np , isplit

    np = size(a)
    if (np > 1) then
     call partition(a, isplit)
     call qsort(a(:isplit-1))
     call qsort(a(isplit:))
    end if
  end subroutine qsort

  subroutine partition(a, marker)
    implicit none
    real(rkx) , dimension(:) , intent(inout) :: a
    integer(ik4) , intent(out) :: marker
    integer(ik4) :: np , left , right
    real(rkx) :: temp , pivot

    np = size(a)
    pivot = (a(1) + a(np))/2.0_rkx
    left = 0
    right = np + 1

    do while (left < right)
      right = right - 1
      do while (a(right) > pivot)
        right = right-1
      end do
      left = left + 1
      do while (a(left) < pivot)
        left = left + 1
      end do
      if (left < right) then
        temp = a(left)
        a(left) = a(right)
        a(right) = temp
      end if
    end do
    if (left == right) then
      marker = left + 1
    else
      marker = left
    end if
  end subroutine partition

  subroutine gfopen(gfile,cfile,xlat,xlon,ds,roi,iband)
    use netcdf
    implicit none
    type(globalfile) , intent(out) :: gfile
    character(len=*) , intent(in) :: cfile
    real(rkx) , dimension(:,:) , intent(in) :: xlat , xlon
    real(rkx) , intent(in) :: ds , roi
    integer(ik4) , intent(in) :: iband
    integer(ik4) :: nlat , nlon , n , i , j , js , itf , iti , itile
    real(rk8) , dimension(:) , allocatable :: glat , glon
    real(rkx) , dimension(:) , allocatable :: rglat , rglon
    integer(ik4) :: idimid , idvar
    integer(ik4) :: jlat , ilon
    integer(ik4) , dimension(2) :: istart , icount

#ifdef DEBUG
    write(stdout,*) 'Opening '//trim(cfile)
#endif
    istatus = nf90_open(cfile, nf90_nowrite, gfile%ncid)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error: '//trim(cfile))

    istatus = nf90_inq_dimid(gfile%ncid,'lat',idimid)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_dimid(gfile%ncid,'latitude',idimid)
      if ( istatus /= nf90_noerr ) then
        istatus = nf90_inq_dimid(gfile%ncid,'LAT',idimid)
        if ( istatus /= nf90_noerr ) then
          istatus = nf90_inq_dimid(gfile%ncid,'LATITUDE',idimid)
          if ( istatus /= nf90_noerr ) then
            istatus = nf90_inq_dimid(gfile%ncid,'lsmlat',idimid)
          end if
        end if
      end if
    end if
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Missing lat dimension in file '//trim(cfile))
    istatus = nf90_inquire_dimension(gfile%ncid,idimid,len=jlat)
    call checkncerr(istatus,__FILE__,__LINE__, &
         'Error reading lat dimelen in file '//trim(cfile))

    istatus = nf90_inq_dimid(gfile%ncid,'lon',idimid)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_dimid(gfile%ncid,'longitude',idimid)
      if ( istatus /= nf90_noerr ) then
        istatus = nf90_inq_dimid(gfile%ncid,'LON',idimid)
        if ( istatus /= nf90_noerr ) then
          istatus = nf90_inq_dimid(gfile%ncid,'LONGITUDE',idimid)
          if ( istatus /= nf90_noerr ) then
            istatus = nf90_inq_dimid(gfile%ncid,'lsmlon',idimid)
          end if
        end if
      end if
    end if
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Missing lon dimension in file '//trim(cfile))
    istatus = nf90_inquire_dimension(gfile%ncid,idimid,len=ilon)
    call checkncerr(istatus,__FILE__,__LINE__, &
         'Error reading lon dimelen in file '//trim(cfile))

    allocate(glat(jlat))
    allocate(glon(ilon))

    istatus = nf90_inq_varid(gfile%ncid,'lat',idvar)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_varid(gfile%ncid,'latitude',idvar)
      if ( istatus /= nf90_noerr ) then
        istatus = nf90_inq_varid(gfile%ncid,'LAT',idvar)
        if ( istatus /= nf90_noerr ) then
          istatus = nf90_inq_varid(gfile%ncid,'LATITUDE',idvar)
          if ( istatus /= nf90_noerr ) then
            istatus = nf90_inq_varid(gfile%ncid,'lsmlat',idvar)
          end if
        end if
      end if
    end if
    call checkncerr(istatus,__FILE__,__LINE__, &
        'Missing lat variable in file '//trim(cfile))
    istatus = nf90_get_var(gfile%ncid,idvar,glat)
    call checkncerr(istatus,__FILE__,__LINE__, &
        'Error reading lat variable in file '//trim(cfile))

    istatus = nf90_inq_varid(gfile%ncid,'lon',idvar)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_varid(gfile%ncid,'longitude',idvar)
      if ( istatus /= nf90_noerr ) then
        istatus = nf90_inq_varid(gfile%ncid,'LON',idvar)
        if ( istatus /= nf90_noerr ) then
          istatus = nf90_inq_varid(gfile%ncid,'LONGITUDE',idvar)
          if ( istatus /= nf90_noerr ) then
            istatus = nf90_inq_varid(gfile%ncid,'lsmlon',idvar)
          end if
        end if
      end if
    end if
    call checkncerr(istatus,__FILE__,__LINE__, &
        'Missing lon variable in file '//trim(cfile))
    istatus = nf90_get_var(gfile%ncid,idvar,glon)
    call checkncerr(istatus,__FILE__,__LINE__, &
        'Error reading lon variable in file '//trim(cfile))

    call get_window(glat,glon,xlat,xlon,iband,gfile%gdomain)

    nlat = gfile%gdomain%nj
    nlon = sum(gfile%gdomain%ni)
    allocate(rglat(nlat))
    allocate(rglon(nlon))
    rglat = 1e-20_rkx
    rglon = 1e-20_rkx
    do i = 1 , nlat
      rglat(i) = real(glat(gfile%gdomain%jgstart+i-1),rkx)
    end do
    js = 1
    do n = 1 , gfile%gdomain%ntiles
      do j = js , js + gfile%gdomain%ni(n) - 1
        rglon(j) = real(glon(gfile%gdomain%igstart(n) + j - js),rkx)
      end do
      js = js + gfile%gdomain%ni(n)
    end do
    deallocate(glat)
    deallocate(glon)
    call h_interpolator_create(gfile%hint,rglat,rglon,xlat,xlon,ds,roi)
    deallocate(rglat)
    deallocate(rglon)

    istatus = nf90_inq_varid(gfile%ncid,'LANDMASK',idvar)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_varid(gfile%ncid,'landmask',idvar)
      if ( istatus /= nf90_noerr ) then
        istatus = nf90_inq_varid(gfile%ncid,'MASK',idvar)
        if ( istatus /= nf90_noerr ) then
          istatus = nf90_inq_varid(gfile%ncid,'mask',idvar)
        end if
      end if
    end if
    if ( istatus == nf90_noerr ) then
      gfile%lmask = .true.
      allocate(gfile%mask(nlon,nlat))
      iti = 1
      gfile%mask = -10000
      do itile = 1 , gfile%gdomain%ntiles
        istart(1) = gfile%gdomain%igstart(itile)
        icount(1) = gfile%gdomain%ni(itile)
        istart(2) = gfile%gdomain%jgstart
        icount(2) = gfile%gdomain%nj
        itf = iti + gfile%gdomain%ni(itile) - 1
        istatus = nf90_get_var(gfile%ncid,idvar, &
                               gfile%mask(iti:itf,:),istart,icount)
        call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')
        iti = itf + 1
      end do
    else
      gfile%lmask = .false.
    end if
  end subroutine gfopen

  subroutine gfread_2di(gfile,vname,var,idef)
    use netcdf
    implicit none
    type(globalfile) , intent(in) :: gfile
    character(len=*) , intent(in) :: vname
    integer(ik4) , intent(in) :: idef
    integer(ik4) , dimension(:,:) , intent(out) :: var
    integer(ik4) , dimension(:,:) , allocatable :: vread
    integer(ik4) :: nlat , nlon , itile , ivar , iti , itf , j , i
    integer(ik4) , dimension(2) :: istart , icount

    nlat = gfile%gdomain%nj
    nlon = sum(gfile%gdomain%ni)
    allocate(vread(nlon,nlat))

    istatus = nf90_inq_varid(gfile%ncid, vname, ivar)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')

    iti = 1
    vread = -1000000000
    do itile = 1 , gfile%gdomain%ntiles
      istart(1) = gfile%gdomain%igstart(itile)
      icount(1) = gfile%gdomain%ni(itile)
      istart(2) = gfile%gdomain%jgstart
      icount(2) = gfile%gdomain%nj
      itf = iti + gfile%gdomain%ni(itile) - 1
      istatus = nf90_get_var(gfile%ncid,ivar,vread(iti:itf,:),istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')
      iti = itf + 1
    end do
    if ( gfile%lmask ) then
      do i = 1 , nlat
        do j = 1 , nlon
          if ( gfile%mask(j,i) == 0 ) then
            vread(j,i) = idef
          else
            if ( vread(j,i) < 0 ) vread(j,i) = idef
          end if
        end do
      end do
    end if
    call h_interpolate_class(gfile%hint,vread,var)
    deallocate(vread)
  end subroutine gfread_2di

  subroutine gfread_2d_landuse(gfile,vname,var,iw,h2opct,rdef)
    use netcdf
    implicit none
    type(globalfile) , intent(in) :: gfile
    character(len=*) , intent(in) :: vname
    integer(ik4) , intent(in) :: iw
    real(rkx) , intent(in) :: h2opct , rdef
    real(rkx) , dimension(:,:) , intent(out) :: var
    integer(ik4) , dimension(:,:) , allocatable :: vread
    integer(ik4) :: nlat , nlon , itile , ivar , iti , itf
    integer(ik4) , dimension(2) :: istart , icount
    real(rkx) :: unused

    nlat = gfile%gdomain%nj
    nlon = sum(gfile%gdomain%ni)
    allocate(vread(nlon,nlat))
    unused = rdef

    istatus = nf90_inq_varid(gfile%ncid, vname, ivar)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')

    iti = 1
    vread = -1000000000
    do itile = 1 , gfile%gdomain%ntiles
      istart(1) = gfile%gdomain%igstart(itile)
      icount(1) = gfile%gdomain%ni(itile)
      istart(2) = gfile%gdomain%jgstart
      icount(2) = gfile%gdomain%nj
      itf = iti + gfile%gdomain%ni(itile) - 1
      istatus = nf90_get_var(gfile%ncid,ivar,vread(iti:itf,:),istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')
      iti = itf + 1
    end do
    call h_interpolate_class(gfile%hint,vread,var,iw,h2opct)
    deallocate(vread)
  end subroutine gfread_2d_landuse

  subroutine gfread_2d(gfile,vname,var,rdef)
    use netcdf
    implicit none
    type(globalfile) , intent(in) :: gfile
    character(len=*) , intent(in) :: vname
    real(rkx) , intent(in) :: rdef
    real(rkx) , dimension(:,:) , intent(out) :: var
    real(rkx) , dimension(:,:) , allocatable :: vread
    integer(ik4) :: nlat , nlon , itile , ivar , iti , itf , i , j
    integer(ik4) , dimension(2) :: istart , icount

    nlat = gfile%gdomain%nj
    nlon = sum(gfile%gdomain%ni)
    allocate(vread(nlon,nlat))

    istatus = nf90_inq_varid(gfile%ncid, vname, ivar)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')

    iti = 1
    vread = -1000000000
    do itile = 1 , gfile%gdomain%ntiles
      istart(1) = gfile%gdomain%igstart(itile)
      icount(1) = gfile%gdomain%ni(itile)
      istart(2) = gfile%gdomain%jgstart
      icount(2) = gfile%gdomain%nj
      itf = iti + gfile%gdomain%ni(itile) - 1
      istatus = nf90_get_var(gfile%ncid,ivar,vread(iti:itf,:),istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')
      iti = itf + 1
    end do

    if ( gfile%lmask ) then
      do i = 1 , nlat
        do j = 1 , nlon
          if ( gfile%mask(j,i) == 0 ) then
            vread(j,i) = h_missing_value
          else
            if ( vread(j,i) < 0.0_rkx ) vread(j,i) = rdef
          end if
        end do
      end do
    else
      do i = 1 , nlat
        do j = 1 , nlon
          if ( vread(j,i) < 0.0_rkx ) vread(j,i) = rdef
        end do
      end do
    end if

    call h_interpolate_cont(gfile%hint,vread,var)
    deallocate(vread)
  end subroutine gfread_2d

  subroutine gfread_2d3d(gfile,vname,var,isel,rdef)
    use netcdf
    implicit none
    type(globalfile) , intent(in) :: gfile
    character(len=*) , intent(in) :: vname
    integer(ik4) , intent(in) :: isel
    real(rkx) , intent(in) :: rdef
    real(rkx) , dimension(:,:) , intent(out) :: var
    real(rkx) , dimension(:,:) , allocatable :: vread
    integer(ik4) :: nlat , nlon , itile , ivar , iti , itf , nd , i , j
    integer(ik4) , dimension(3) :: idims , istart , icount

    if ( isel < 0 ) then
      call die('rdldtr','Variable '//trim(vname)// &
                        ' requested slice not in file')
    end if

    nlat = gfile%gdomain%nj
    nlon = sum(gfile%gdomain%ni)
    allocate(vread(nlon,nlat))

    istatus = nf90_inq_varid(gfile%ncid, vname, ivar)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')
    istatus = nf90_inquire_variable(gfile%ncid,ivar,dimids=idims)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')
    istatus = nf90_inquire_dimension(gfile%ncid,idims(3),len=nd)
    call checkncerr(istatus,__FILE__,__LINE__, &
         'Error reading dimelens for variable '//trim(vname))

    if ( isel > nd ) then
      call die('rdldtr','Variable '//trim(vname)// &
                        ' requested slice not in file')
    end if

    iti = 1
    vread = -1000000000
    do itile = 1 , gfile%gdomain%ntiles
      istart(1) = gfile%gdomain%igstart(itile)
      icount(1) = gfile%gdomain%ni(itile)
      istart(2) = gfile%gdomain%jgstart
      icount(2) = gfile%gdomain%nj
      istart(3) = isel
      icount(3) = 1
      itf = iti + gfile%gdomain%ni(itile) - 1
      istatus = nf90_get_var(gfile%ncid,ivar,vread(iti:itf,:),istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')
      iti = itf + 1
    end do

    if ( gfile%lmask ) then
      do i = 1 , nlat
        do j = 1 , nlon
          if ( gfile%mask(j,i) == 0 ) then
            vread(j,i) = h_missing_value
          else
            if ( vread(j,i) < 0.0_rkx ) vread(j,i) = rdef
          end if
        end do
      end do
    else
      do i = 1 , nlat
        do j = 1 , nlon
          if ( vread(j,i) < 0.0_rkx ) vread(j,i) = rdef
        end do
      end do
    end if

    call h_interpolate_cont(gfile%hint,vread,var)
    deallocate(vread)
  end subroutine gfread_2d3d

  subroutine gfread_3d(gfile,vname,var,rdef)
    use netcdf
    implicit none
    type(globalfile) , intent(in) :: gfile
    character(len=*) , intent(in) :: vname
    real(rkx) , intent(in) :: rdef
    real(rkx) , dimension(:,:,:) , intent(out) :: var
    real(rkx) , dimension(:,:,:) , allocatable :: vread
    integer(ik4) :: nlat , nlon , itile , ivar , iti , itf , nd , i , j , n
    integer(ik4) , dimension(3) :: idims , istart , icount

    istatus = nf90_inq_varid(gfile%ncid, vname, ivar)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')
    istatus = nf90_inquire_variable(gfile%ncid,ivar,dimids=idims)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')
    istatus = nf90_inquire_dimension(gfile%ncid,idims(3),len=nd)
    call checkncerr(istatus,__FILE__,__LINE__, &
         'Error reading dimelens for variable '//trim(vname))

    nlat = gfile%gdomain%nj
    nlon = sum(gfile%gdomain%ni)
    allocate(vread(nlon,nlat,nd))

    iti = 1
    vread = -1000000000
    do itile = 1 , gfile%gdomain%ntiles
      istart(1) = gfile%gdomain%igstart(itile)
      icount(1) = gfile%gdomain%ni(itile)
      istart(2) = gfile%gdomain%jgstart
      icount(2) = gfile%gdomain%nj
      istart(3) = 1
      icount(3) = nd
      itf = iti + gfile%gdomain%ni(itile) - 1
      istatus = nf90_get_var(gfile%ncid,ivar,vread(iti:itf,:,:),istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')
      iti = itf + 1
    end do

    if ( gfile%lmask ) then
      do n = 1 , nd
        do i = 1 , nlat
          do j = 1 , nlon
            if ( gfile%mask(j,i) == 0 ) then
              vread(j,i,n) = h_missing_value
            else
              if ( vread(j,i,n) < 0.0_rkx ) vread(j,i,n) = rdef
            end if
          end do
        end do
      end do
    else
      do n = 1 , nd
        do i = 1 , nlat
          do j = 1 , nlon
            if ( vread(j,i,n) < 0.0_rkx ) vread(j,i,n) = rdef
          end do
        end do
      end do
    end if

    call h_interpolate_cont(gfile%hint,vread,var)
    deallocate(vread)
  end subroutine gfread_3d

  subroutine gfread_3d_lookup(gfile,vname,lkdim,lkvar,var,lrev,rdef)
    use netcdf
    implicit none
    type(globalfile) , intent(in) :: gfile
    character(len=*) , intent(in) :: vname
    character(len=*) , intent(in) :: lkdim
    character(len=*) , intent(in) :: lkvar
    logical , intent(in) :: lrev
    real(rkx) , intent(in) :: rdef
    real(rkx) , dimension(:,:,:) , intent(out) :: var
    real(rkx) , dimension(:,:,:) , allocatable :: vread
    real(rkx) , dimension(:,:) , allocatable :: temp
    integer(ik4) , dimension(:,:) , allocatable :: mpu
    integer(ik4) :: nlat , nlon , nlev , nmpu , i , j , n
    integer(ik4) :: itile , ivar , ilook , impud , iti , itf
    integer(ik4) , dimension(2) :: istart , icount

    nlat = gfile%gdomain%nj
    nlon = sum(gfile%gdomain%ni)
    nlev = size(var,3)
    allocate(mpu(nlon,nlat))
    allocate(vread(nlon,nlat,nlev))

    istatus = nf90_inq_varid(gfile%ncid, vname, ivar)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')
    istatus = nf90_inq_varid(gfile%ncid, lkvar, ilook)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')

    istatus = nf90_inq_dimid(gfile%ncid,lkdim,impud)
    call checkncerr(istatus,__FILE__,__LINE__, &
        'Cannot find dimension '//trim(lkdim)//' in file')

    istatus = nf90_inquire_dimension(gfile%ncid,impud,len=nmpu)
    call checkncerr(istatus,__FILE__,__LINE__, &
        'Cannot read dimension '//trim(lkdim)//' in file')

    if ( lrev ) then
      allocate(temp(nlev,nmpu))
    else
      allocate(temp(nmpu,nlev))
    end if

    istatus = nf90_get_var(gfile%ncid,ivar,temp)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')

    iti = 1
    vread = -1000000000
    do itile = 1 , gfile%gdomain%ntiles
      istart(1) = gfile%gdomain%igstart(itile)
      icount(1) = gfile%gdomain%ni(itile)
      istart(2) = gfile%gdomain%jgstart
      icount(2) = gfile%gdomain%nj
      itf = iti + gfile%gdomain%ni(itile) - 1
      istatus = nf90_get_var(gfile%ncid,ilook,mpu(iti:itf,:),istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')
      iti = itf + 1
    end do

    if ( lrev ) then
      do n = 1 , nlev
        do i = 1 , nlat
          do j = 1 , nlon
            if ( mpu(j,i) > 0 .and. mpu(j,i) <= nmpu ) then
              vread(j,i,n) = temp(n,mpu(j,i))
            else
              vread(j,i,n) = h_missing_value
            end if
          end do
        end do
      end do
    else
      do n = 1 , nlev
        do i = 1 , nlat
          do j = 1 , nlon
            if ( mpu(j,i) > 0 .and. mpu(j,i) <= nmpu ) then
              vread(j,i,n) = temp(mpu(j,i),n)
            else
              vread(j,i,n) = h_missing_value
            end if
          end do
        end do
      end do
    end if

    if ( gfile%lmask ) then
      do n = 1 , nlev
        do i = 1 , nlat
          do j = 1 , nlon
            if ( gfile%mask(j,i) == 0 ) then
              vread(j,i,n) = h_missing_value
            else
              if ( vread(j,i,n) < 0.0_rkx ) vread(j,i,n) = rdef
            end if
          end do
        end do
      end do
    else
      do n = 1 , nlev
        do i = 1 , nlat
          do j = 1 , nlon
            if ( vread(j,i,n) < 0.0_rkx ) vread(j,i,n) = rdef
          end do
        end do
      end do
    end if

    call h_interpolate_cont(gfile%hint,vread,var)
    deallocate(vread)
    deallocate(temp)
    deallocate(mpu)
  end subroutine gfread_3d_lookup

  subroutine gfread_4d(gfile,vname,var,rdef)
    use netcdf
    implicit none
    type(globalfile) , intent(in) :: gfile
    character(len=*) , intent(in) :: vname
    real(rkx) , intent(in) :: rdef
    real(rkx) , dimension(:,:,:,:) , intent(out) :: var
    real(rkx) , dimension(:,:,:,:) , allocatable :: vread
    integer(ik4) :: nlat , nlon , itile , ivar , iti , itf
    integer(ik4) :: nd , ne , i , j , m , n
    integer(ik4) , dimension(4) :: idims , istart , icount

    istatus = nf90_inq_varid(gfile%ncid, vname, ivar)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')
    istatus = nf90_inquire_variable(gfile%ncid,ivar,dimids=idims)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')
    istatus = nf90_inquire_dimension(gfile%ncid,idims(3),len=nd)
    call checkncerr(istatus,__FILE__,__LINE__, &
         'Error reading dimelens for variable '//trim(vname))
    istatus = nf90_inquire_dimension(gfile%ncid,idims(4),len=ne)
    call checkncerr(istatus,__FILE__,__LINE__, &
         'Error reading dimelens for variable '//trim(vname))

    nlat = gfile%gdomain%nj
    nlon = sum(gfile%gdomain%ni)
    allocate(vread(nlon,nlat,nd,ne))

    iti = 1
    vread = -1000000000
    do itile = 1 , gfile%gdomain%ntiles
      istart(1) = gfile%gdomain%igstart(itile)
      icount(1) = gfile%gdomain%ni(itile)
      istart(2) = gfile%gdomain%jgstart
      icount(2) = gfile%gdomain%nj
      istart(3) = 1
      icount(3) = nd
      istart(4) = 1
      icount(4) = ne
      itf = iti + gfile%gdomain%ni(itile) - 1
      istatus = nf90_get_var(gfile%ncid,ivar,vread(iti:itf,:,:,:),istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')
      iti = itf + 1
    end do

    if ( gfile%lmask ) then
      do m = 1 , ne
        do n = 1 , nd
          do i = 1 , nlat
            do j = 1 , nlon
              if ( gfile%mask(j,i) == 0 ) then
                vread(j,i,n,m) = h_missing_value
              else
                if ( vread(j,i,n,m) < 0.0_rkx ) vread(j,i,n,m) = rdef
              end if
            end do
          end do
        end do
      end do
    else
      do m = 1 , ne
        do n = 1 , nd
          do i = 1 , nlat
            do j = 1 , nlon
              if ( vread(j,i,n,m) < 0.0_rkx ) vread(j,i,n,m) = rdef
            end do
          end do
        end do
      end do
    end if

    call h_interpolate_cont(gfile%hint,vread,var)
    deallocate(vread)
  end subroutine gfread_4d

  subroutine gfread_4d_lookup(gfile,vname,lkdim,lkvar,var,lrev,rdef)
    use netcdf
    implicit none
    type(globalfile) , intent(in) :: gfile
    character(len=*) , intent(in) :: vname
    character(len=*) , intent(in) :: lkdim
    character(len=*) , intent(in) :: lkvar
    logical , intent(in) :: lrev
    real(rkx) , intent(in) :: rdef
    real(rkx) , dimension(:,:,:,:) , intent(out) :: var
    real(rkx) , dimension(:,:,:,:) , allocatable :: vread
    real(rkx) , dimension(:,:,:) , allocatable :: temp
    integer(ik4) , dimension(:,:) , allocatable :: mpu
    integer(ik4) :: nlat , nlon , nlev1 , nlev2 , nmpu
    integer(ik4) :: i , j , n1 , n2
    integer(ik4) :: itile , ivar , ilook , impud , iti , itf
    integer(ik4) , dimension(2) :: istart , icount

    nlat = gfile%gdomain%nj
    nlon = sum(gfile%gdomain%ni)
    nlev1 = size(var,3)
    nlev2 = size(var,4)
    allocate(mpu(nlon,nlat))
    allocate(vread(nlon,nlat,nlev1,nlev2))

    istatus = nf90_inq_varid(gfile%ncid, vname, ivar)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')
    istatus = nf90_inq_varid(gfile%ncid, lkvar, ilook)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')

    istatus = nf90_inq_dimid(gfile%ncid,lkdim,impud)
    call checkncerr(istatus,__FILE__,__LINE__, &
        'Cannot find dimension '//trim(lkdim)//' in file')

    istatus = nf90_inquire_dimension(gfile%ncid,impud,len=nmpu)
    call checkncerr(istatus,__FILE__,__LINE__, &
        'Cannot read dimension '//trim(lkdim)//' in file')

    if ( lrev ) then
      allocate(temp(nlev2,nmpu,nlev1))
    else
      allocate(temp(nlev1,nmpu,nlev2))
    end if

    istatus = nf90_get_var(gfile%ncid,ivar,temp)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')

    iti = 1
    vread = -1000000000
    do itile = 1 , gfile%gdomain%ntiles
      istart(1) = gfile%gdomain%igstart(itile)
      icount(1) = gfile%gdomain%ni(itile)
      istart(2) = gfile%gdomain%jgstart
      icount(2) = gfile%gdomain%nj
      itf = iti + gfile%gdomain%ni(itile) - 1
      istatus = nf90_get_var(gfile%ncid,ilook,mpu(iti:itf,:),istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')
      iti = itf + 1
    end do

    if ( lrev ) then
      do n2 = 1 , nlev2
        do n1 = 1 , nlev1
          do i = 1 , nlat
            do j = 1 , nlon
              if ( mpu(j,i) > 0 .and. mpu(j,i) <= nmpu ) then
                vread(j,i,n1,n2) = temp(n2,mpu(j,i),n1)
              else
                vread(j,i,n1,n2) = h_missing_value
              end if
            end do
          end do
        end do
      end do
    else
      do n2 = 1 , nlev2
        do n1 = 1 , nlev1
          do i = 1 , nlat
            do j = 1 , nlon
              if ( mpu(j,i) > 0 .and. mpu(j,i) <= nmpu ) then
                vread(j,i,n1,n2) = temp(n1,mpu(j,i),n2)
              else
                vread(j,i,n1,n2) = h_missing_value
              end if
            end do
          end do
        end do
      end do
    end if

    if ( gfile%lmask ) then
      do n2 = 1 , nlev2
        do n1 = 1 , nlev1
          do i = 1 , nlat
            do j = 1 , nlon
              if ( gfile%mask(j,i) == 0 ) then
                vread(j,i,n1,n2) = h_missing_value
              else
                if ( vread(j,i,n1,n2) < 0.0_rkx ) vread(j,i,n1,n2) = rdef
              end if
            end do
          end do
        end do
      end do
    else
      do n2 = 1 , nlev2
        do n1 = 1 , nlev1
          do i = 1 , nlat
            do j = 1 , nlon
              if ( vread(j,i,n1,n2) < 0.0_rkx ) vread(j,i,n1,n2) = rdef
            end do
          end do
        end do
      end do
    end if

    call h_interpolate_cont(gfile%hint,vread,var)
    deallocate(vread)
    deallocate(temp)
    deallocate(mpu)
  end subroutine gfread_4d_lookup

  subroutine gfread_5d(gfile,vname,var,rdef)
    use netcdf
    implicit none
    type(globalfile) , intent(in) :: gfile
    character(len=*) , intent(in) :: vname
    real(rkx) , intent(in) :: rdef
    real(rkx) , dimension(:,:,:,:,:) , intent(out) :: var
    real(rkx) , dimension(:,:,:,:,:) , allocatable :: vread
    integer(ik4) :: nlat , nlon , itile , ivar , iti , itf , nd , ne , nf
    integer(ik4) :: i , j , l , m , n
    integer(ik4) , dimension(5) :: idims , istart , icount

    istatus = nf90_inq_varid(gfile%ncid, vname, ivar)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')
    istatus = nf90_inquire_variable(gfile%ncid,ivar,dimids=idims)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')
    istatus = nf90_inquire_dimension(gfile%ncid,idims(3),len=nd)
    call checkncerr(istatus,__FILE__,__LINE__, &
         'Error reading dimelens for variable '//trim(vname))
    istatus = nf90_inquire_dimension(gfile%ncid,idims(4),len=ne)
    call checkncerr(istatus,__FILE__,__LINE__, &
         'Error reading dimelens for variable '//trim(vname))
    istatus = nf90_inquire_dimension(gfile%ncid,idims(5),len=nf)
    call checkncerr(istatus,__FILE__,__LINE__, &
         'Error reading dimelens for variable '//trim(vname))

    nlat = gfile%gdomain%nj
    nlon = sum(gfile%gdomain%ni)
    allocate(vread(nlon,nlat,nd,ne,nf))

    iti = 1
    vread = -1000000000
    do itile = 1 , gfile%gdomain%ntiles
      istart(1) = gfile%gdomain%igstart(itile)
      icount(1) = gfile%gdomain%ni(itile)
      istart(2) = gfile%gdomain%jgstart
      icount(2) = gfile%gdomain%nj
      istart(3) = 1
      icount(3) = nd
      istart(4) = 1
      icount(4) = ne
      istart(5) = 1
      icount(5) = nf
      itf = iti + gfile%gdomain%ni(itile) - 1
      istatus = nf90_get_var(gfile%ncid,ivar, &
                             vread(iti:itf,:,:,:,:),istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')
      iti = itf + 1
    end do

    if ( gfile%lmask ) then
      do l = 1 , nd
        do m = 1 , ne
          do n = 1 , nf
            do i = 1 , nlat
              do j = 1 , nlon
                if ( gfile%mask(j,i) == 0 ) then
                  vread(j,i,n,m,l) = h_missing_value
                else
                  if ( vread(j,i,n,m,l) < 0.0_rkx ) vread(j,i,n,m,l) = rdef
                end if
              end do
            end do
          end do
        end do
      end do
    else
      do l = 1 , nd
        do m = 1 , ne
          do n = 1 , nf
            do i = 1 , nlat
              do j = 1 , nlon
                if ( vread(j,i,n,m,l) < 0.0_rkx ) vread(j,i,n,m,l) = rdef
              end do
            end do
          end do
        end do
      end do
    end if

    call h_interpolate_cont(gfile%hint,vread,var)
    deallocate(vread)
  end subroutine gfread_5d

  subroutine gfread_5d_lookup(gfile,vname,lkdim,lkvar,var,lrev,rdef)
    use netcdf
    implicit none
    type(globalfile) , intent(in) :: gfile
    character(len=*) , intent(in) :: vname
    character(len=*) , intent(in) :: lkdim
    character(len=*) , intent(in) :: lkvar
    logical , intent(in) :: lrev
    real(rkx) , intent(in) :: rdef
    real(rkx) , dimension(:,:,:,:,:) , intent(out) :: var
    real(rkx) , dimension(:,:,:,:,:) , allocatable :: vread
    real(rkx) , dimension(:,:,:,:) , allocatable :: temp
    integer(ik4) , dimension(:,:) , allocatable :: mpu
    integer(ik4) :: nlat , nlon , nlev1 , nlev2 , nlev3 , nmpu
    integer(ik4) :: i , j , n1 , n2 , n3
    integer(ik4) :: itile , ivar , ilook , impud , iti , itf
    integer(ik4) , dimension(2) :: istart , icount
    integer(ik4) , dimension(4) :: vdims , dimlen

    nlat = gfile%gdomain%nj
    nlon = sum(gfile%gdomain%ni)
    nlev1 = size(var,3)
    nlev2 = size(var,4)
    nlev3 = size(var,5)
    allocate(mpu(nlon,nlat))
    allocate(vread(nlon,nlat,nlev1,nlev2,nlev3))

    istatus = nf90_inq_varid(gfile%ncid, vname, ivar)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')
    istatus = nf90_inquire_variable(gfile%ncid,ivar,dimids=vdims)
    call checkncerr(istatus,__FILE__,__LINE__, &
        'Cannot read variable dimension '//trim(vname)//' in file')
    do i = 1 , 4
      istatus = nf90_inquire_dimension(gfile%ncid,vdims(i),len=dimlen(i))
    end do

    istatus = nf90_inq_varid(gfile%ncid, lkvar, ilook)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')

    istatus = nf90_inq_dimid(gfile%ncid,lkdim,impud)
    call checkncerr(istatus,__FILE__,__LINE__, &
        'Cannot find dimension '//trim(lkdim)//' in file')

    istatus = nf90_inquire_dimension(gfile%ncid,impud,len=nmpu)
    call checkncerr(istatus,__FILE__,__LINE__, &
        'Cannot read dimension '//trim(lkdim)//' in file')

    if ( lrev ) then
      allocate(temp(nlev3,nmpu,nlev2,nlev1))
    else
      allocate(temp(nlev1,nmpu,nlev2,nlev3))
    end if

    istatus = nf90_get_var(gfile%ncid,ivar,temp)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')

    iti = 1
    vread = -1000000000
    do itile = 1 , gfile%gdomain%ntiles
      istart(1) = gfile%gdomain%igstart(itile)
      icount(1) = gfile%gdomain%ni(itile)
      istart(2) = gfile%gdomain%jgstart
      icount(2) = gfile%gdomain%nj
      itf = iti + gfile%gdomain%ni(itile) - 1
      istatus = nf90_get_var(gfile%ncid,ilook,mpu(iti:itf,:),istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')
      iti = itf + 1
    end do

    if ( lrev ) then
      do n3 = 1 , nlev3
        do n2 = 1 , nlev2
          do n1 = 1 , nlev1
            do i = 1 , nlat
              do j = 1 , nlon
                if ( mpu(j,i) > 0 .and. mpu(j,i) <= nmpu ) then
                  vread(j,i,n1,n2,n3) = temp(n3,mpu(j,i),n2,n1)
                else
                  vread(j,i,n1,n2,n3) = h_missing_value
                end if
              end do
            end do
          end do
        end do
      end do
    else
      do n3 = 1 , nlev3
        do n2 = 1 , nlev2
          do n1 = 1 , nlev1
            do i = 1 , nlat
              do j = 1 , nlon
                if ( mpu(j,i) > 0 .and. mpu(j,i) <= nmpu ) then
                  vread(j,i,n1,n2,n3) = temp(n1,mpu(j,i),n2,n3)
                else
                  vread(j,i,n1,n2,n3) = h_missing_value
                end if
              end do
            end do
          end do
        end do
      end do
    end if

    if ( gfile%lmask ) then
      do n3 = 1 , nlev3
        do n2 = 1 , nlev2
          do n1 = 1 , nlev1
            do i = 1 , nlat
              do j = 1 , nlon
                if ( gfile%mask(j,i) < 0 ) then
                  vread(j,i,n1,n2,n3) = h_missing_value
                else
                  if ( vread(j,i,n1,n2,n3) < 0.0_rkx ) then
                    vread(j,i,n1,n2,n3) = rdef
                  end if
                end if
              end do
            end do
          end do
        end do
      end do
    else
      do n3 = 1 , nlev3
        do n2 = 1 , nlev2
          do n1 = 1 , nlev1
            do i = 1 , nlat
              do j = 1 , nlon
                if ( vread(j,i,n1,n2,n3) < 0.0_rkx ) then
                  vread(j,i,n1,n2,n3) = rdef
                end if
              end do
            end do
          end do
        end do
      end do
    end if

    call h_interpolate_cont(gfile%hint,vread,var)
    deallocate(vread)
    deallocate(temp)
    deallocate(mpu)
  end subroutine gfread_5d_lookup

  subroutine gfclose(gfile)
    use netcdf
    implicit none
    type(globalfile) , intent(inout) :: gfile
    integer(ik4) :: istatus
    call h_interpolator_destroy(gfile%hint)
    if ( gfile%ncid > 0 ) then
      istatus = nf90_close(gfile%ncid)
      call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')
    end if
    if ( gfile%lmask ) then
      deallocate(gfile%mask)
    end if
  end subroutine gfclose

  subroutine bestaround2d(grid,i,j)
    implicit none
    integer(ik4) , intent(in) :: i , j
    real(rkx) , dimension(:,:) , intent(inout) :: grid
    integer(ik4) :: ii , jj , js , is , il
    integer(ik4) :: iis , jjs , iie , jje
    integer , parameter :: maxil = 10

    jjs = lbound(grid,1)
    iis = lbound(grid,2)
    jje = ubound(grid,1)
    iie = ubound(grid,2)

    il = 1
    do
      do ii = i - il , i + il
        do jj = j - il , j + il
          is = ii
          js = jj
          if ( js < jjs ) js = jjs-js
          if ( js > jje ) js = 2*jje-js
          if ( is < iis ) is = iis-is
          if ( is > iie ) is = 2*iie - is
          if ( grid(js,is) > h_missing_value ) then
            grid(j,i) = grid(js,is)
            return
          end if
        end do
      end do
      il = il + 1
      if ( il == maxil ) then
        grid(j,i) = h_missing_value
        return
      end if
    end do
  end subroutine bestaround2d

  subroutine bestaround3d(grid,i,j)
    implicit none
    integer(ik4) , intent(in) :: i , j
    real(rkx) , dimension(:,:,:) , intent(inout) :: grid
    integer(ik4) :: ii , jj , js , is , n , il
    integer(ik4) :: iis , jjs , iie , jje , kks , kke
    integer , parameter :: maxil = 10

    jjs = lbound(grid,1)
    iis = lbound(grid,2)
    kks = lbound(grid,3)
    jje = ubound(grid,1)
    iie = ubound(grid,2)
    kke = ubound(grid,3)

    il = 1
    do
      do ii = i - il , i + il
        do jj = j - il , j + il
          is = ii
          js = jj
          if ( js < jjs ) js = jjs-js
          if ( js > jje ) js = 2*jje - js
          if ( is < iis ) is = iis-is
          if ( is > iie ) is = 2*iie - is
          if ( grid(js,is,kks) > h_missing_value ) then
            do n = kks , kke
              grid(j,i,n) = grid(js,is,n)
            end do
            return
          end if
        end do
      end do
      il = il + 1
      if ( il == maxil ) then
        do n = kks , kke
          grid(j,i,n) = h_missing_value
        end do
        return
      end if
    end do
  end subroutine bestaround3d

end module mod_rdldtr
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

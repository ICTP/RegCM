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
  use mod_interp

  private

  public :: read_ncglob , longitude_circle

  real(rkx) , allocatable , dimension(:) :: glat
  real(rkx) , allocatable , dimension(:) :: glon

  integer(ik4) :: ncid , istatus
  type(global_domain) :: gdomain

  interface read_ncglob
    module procedure read_ncglob2d
    module procedure read_ncglob3d
    module procedure read_ncglob2d3d
  end interface read_ncglob

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
    real(rkx) , intent(out) :: grdlnma , grdlnmn , grdltma , grdltmn
    integer(ik4) , intent(out) :: nlatin , nlonin
    real(rkx) , dimension(:,:) , intent(out) , pointer :: values
    integer(ik4) :: nlat , nlon , iti , itf , itile , ivar
    integer(ik4) :: i , j , inpsec , iopsec , ifrac
    integer(ik4) , dimension(2) :: istart , icount
    real(rkx) :: delta
    real(rkx) , allocatable , dimension(:,:) :: readbuf

    write(stdout,*) 'Opening '//trim(cfile)
    istatus = nf90_open(cfile, nf90_nowrite, ncid)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')

    call read_geolocation(cfile)

    istatus = nf90_inq_varid(ncid, cvar, ivar)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')

    call get_window(glat,glon,xlat,xlon,iband,gdomain)

    nlat = gdomain%nj
    nlon = sum(gdomain%ni)

    inpsec = int(abs(glat(2)-glat(1))*3600.0_rkx)
    iopsec = max(int(dble(iores*60.0_rkx)),1)
    ifrac = max(iopsec/inpsec,1)
#ifdef DEBUG
    write(stderr,*) 'INPSEC = ', inpsec
    write(stderr,*) 'IOPSEC = ', iopsec
    write(stderr,*) 'IFRAC  = ', ifrac
#endif

    nlatin = (nlat/ifrac)
    nlonin = (nlon/ifrac)

    delta = 0.0_rkx
    if ( ifrac > 1 ) then
      delta = (dble(iopsec)/2.0_rkx) / 3600.0_rkx
    end if

    grdlnmn = glon(gdomain%igstart(1)) + delta
    grdlnma = glon(gdomain%igstart(gdomain%ntiles) + &
                   gdomain%ni(gdomain%ntiles)-1) - delta
    grdltmn = glat(gdomain%jgstart) + delta
    grdltma = glat(gdomain%jgstart+gdomain%nj-1) - delta
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
      write (stdout,*) 'Correcting South pole.'
      readbuf(:,1) = readbuf(:,2)
    else if ( istart(2)+icount(2)-1 == gdomain%global_nj ) then
      write (stdout,*) 'Correcting North pole.'
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
    write(stdout,'(a)') ' Done.'
  end subroutine read_ncglob2d

  subroutine read_ncglob3d(cfile,cvar,iores,imeth,iband,xlat,xlon, &
                           grdlnma,grdlnmn,grdltma,grdltmn,        &
                           nlatin,nlonin,values)
    use netcdf
    implicit none
    character(len=*) , intent(in) :: cfile , cvar
    integer(ik4) , intent(in) :: iband , iores , imeth
    real(rkx) , dimension(:,:) , intent(in) :: xlat , xlon
    real(rkx) , intent(out) :: grdlnma , grdlnmn , grdltma , grdltmn
    integer(ik4) , intent(out) :: nlatin , nlonin
    real(rkx) , dimension(:,:,:) , intent(out) , pointer :: values
    integer(ik4) :: nlat , nlon , iti , itf , itile , ivar
    integer(ik4) :: i , j , n , inpsec , iopsec , ifrac , nd
    integer(ik4) , dimension(3) :: idims , istart , icount
    real(rkx) :: delta
    real(rkx) , allocatable , dimension(:,:,:) :: readbuf

    write(stdout,*) 'Opening '//trim(cfile)
    istatus = nf90_open(cfile, nf90_nowrite, ncid)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')

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

    inpsec = int(abs(glat(2)-glat(1))*3600.0_rkx)
    iopsec = max(int(dble(iores*60.0_rkx)),1)
    ifrac = max(iopsec/inpsec,1)
#ifdef DEBUG
    write(stderr,*) 'INPSEC = ', inpsec
    write(stderr,*) 'IOPSEC = ', iopsec
    write(stderr,*) 'IFRAC  = ', ifrac
#endif

    nlatin = (nlat/ifrac)
    nlonin = (nlon/ifrac)

    delta = 0.0_rkx
    if ( ifrac > 1 ) then
      delta = (dble(iopsec)/2.0_rkx) / 3600.0_rkx
    end if

    grdlnmn = glon(gdomain%igstart(1)) + delta
    grdlnma = glon(gdomain%igstart(gdomain%ntiles) + &
                   gdomain%ni(gdomain%ntiles)-1) - delta
    grdltmn = glat(gdomain%jgstart) + delta
    grdltma = glat(gdomain%jgstart+gdomain%nj-1) - delta
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
      write (stdout,*) 'Correcting South pole.'
      readbuf(:,1,:) = readbuf(:,2,:)
    else if ( istart(2)+icount(2)-1 == gdomain%global_nj ) then
      write (stdout,*) 'Correcting North pole.'
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
    write(stdout,'(a)') ' Done.'
  end subroutine read_ncglob3d

  subroutine read_ncglob2d3d(cfile,cvar,iores,imeth,iband,xlat,xlon, &
                             grdlnma,grdlnmn,grdltma,grdltmn,        &
                             nlatin,nlonin,values,isel)
    use netcdf
    implicit none
    character(len=*) , intent(in) :: cfile , cvar
    integer(ik4) , intent(in) :: iband , iores , imeth
    real(rkx) , dimension(:,:) , intent(in) :: xlat , xlon
    real(rkx) , intent(out) :: grdlnma , grdlnmn , grdltma , grdltmn
    integer(ik4) , intent(out) :: nlatin , nlonin
    real(rkx) , dimension(:,:) , intent(out) , pointer :: values
    integer(ik4) , intent(in) :: isel
    integer(ik4) :: nlat , nlon , iti , itf , itile , ivar
    integer(ik4) :: i , j , inpsec , iopsec , ifrac , nd
    integer(ik4) , dimension(3) :: idims , istart , icount
    real(rkx) :: delta
    real(rkx) , allocatable , dimension(:,:) :: readbuf

    write(stdout,*) 'Opening '//trim(cfile)
    istatus = nf90_open(cfile, nf90_nowrite, ncid)
    call checkncerr(istatus,__FILE__,__LINE__,'NetCDF Error')

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

    inpsec = int(abs(glat(2)-glat(1))*3600.0_rkx)
    iopsec = max(int(dble(iores*60.0_rkx)),1)
    ifrac = max(iopsec/inpsec,1)
#ifdef DEBUG
    write(stderr,*) 'INPSEC = ', inpsec
    write(stderr,*) 'IOPSEC = ', iopsec
    write(stderr,*) 'IFRAC  = ', ifrac
#endif

    nlatin = (nlat/ifrac)
    nlonin = (nlon/ifrac)

    delta = 0.0_rkx
    if ( ifrac > 1 ) then
      delta = (dble(iopsec)/2.0_rkx) / 3600.0_rkx
    end if

    grdlnmn = glon(gdomain%igstart(1)) + delta
    grdlnma = glon(gdomain%igstart(gdomain%ntiles) + &
                   gdomain%ni(gdomain%ntiles)-1) - delta
    grdltmn = glat(gdomain%jgstart) + delta
    grdltma = glat(gdomain%jgstart+gdomain%nj-1) - delta
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
      write (stdout,*) 'Correcting South pole.'
      readbuf(:,1) = readbuf(:,2)
    else if ( istart(2)+icount(2)-1 == gdomain%global_nj ) then
      write (stdout,*) 'Correcting North pole.'
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
    write(stdout,'(a)') ' Done.'
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
    integer(ik4) :: nfrac , np , i , j
    real(rkx) , allocatable , dimension(:) :: copybuf

    write(stdout,'(a)',advance='no') ' Resampling'
    nfrac = ifrac*ifrac
    allocate(copybuf(nfrac))
    select case (imeth)
      case (1)
        do i = 1 , nlatin
          if (mod(i,10) == 0) write(stdout,'(a)',advance='no') '.'
          do j = 1 , nlonin
            call fillbuf(copybuf,readbuf,nlon,nlat,(j-1)*ifrac+1,&
                         (i-1)*ifrac+1,ifrac,iband)
            values(j,i) = sum(copybuf)/real(size(copybuf),rkx)
          end do
        end do
      case (2)
        do i = 1 , nlatin
          if (mod(i,10) == 0) write(stdout,'(a)',advance='no') '.'
          do j = 1 , nlonin
            call fillbuf(copybuf,readbuf,nlon,nlat,(j-1)*ifrac+1,&
                         (i-1)*ifrac+1,ifrac,iband)
            call qsort(copybuf)
            !values(j,i) = 0.5*(copybuf(nfrac/2)+copybuf(nfrac/2+1))
            values(j,i) = copybuf(max(nfrac/2,1))
          end do
        end do
      case (3)
        do i = 1 , nlatin
          if (mod(i,10) == 0) write(stdout,'(a)',advance='no') '.'
          do j = 1 , nlonin
            call fillbuf(copybuf,readbuf,nlon,nlat,(j-1)*ifrac+1,&
                         (i-1)*ifrac+1,ifrac,iband)
            values(j,i) = real(mpindex(copybuf),rkx)
          end do
        end do
      case (4)
        do i = 1 , nlatin
          if (mod(i,10) == 0) write(stdout,'(a)',advance='no') '.'
          do j = 1 , nlonin
            call fillbuf(copybuf,readbuf,nlon,nlat,(j-1)*ifrac+1,&
                         (i-1)*ifrac+1,ifrac,iband)
            call qsort(copybuf)
            np = (ifrac*ifrac)/4
            values(j,i) = sum(copybuf(1+np:size(copybuf)-np+1)) / &
                           real(size(copybuf)-2*np,rkx)
          end do
        end do
      case default
        do i = 1 , nlatin
          if (mod(i,10) == 0) write(stdout,'(a)',advance='no') '.'
          do j = 1 , nlonin
            values(j,i) = readbuf((j-1)*ifrac+1,(i-1)*ifrac+1)
          end do
        end do
    end select
    deallocate(copybuf)
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
    pivot = (a(1) + a(np))/2.0E0
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

  real(rkx) function longitude_circle(lat) result(er)
    implicit none
    real(rkx) , intent(in) :: lat
    er = d_two * mathpi * erkm * cos(lat*degrad)
  end function longitude_circle

end module mod_rdldtr
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

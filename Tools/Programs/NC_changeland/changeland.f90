program changeland
  use netcdf
  implicit none
  integer :: ncid
  integer :: istatus
  character(len=256) :: arg1
  integer :: jx , iy
  integer :: jxdimid , iydimid , ivarid

  real(4) , pointer , dimension(:,:) :: xlat , xlon , landuse
  logical , pointer , dimension(:,:) :: mask

  call getarg(1,arg1)

  istatus = nf90_open(arg1,nf90_write,ncid)
  if ( istatus /= nf90_noerr) then
    print *, nf90_strerror(istatus)
    stop
  end if

  istatus = nf90_inq_dimid(ncid, "jx", jxdimid)
  if ( istatus /= nf90_noerr) then
    print *, nf90_strerror(istatus)
    stop
  end if
  istatus = nf90_inquire_dimension(ncid, jxdimid, len=jx)
  if ( istatus /= nf90_noerr) then
    print *, nf90_strerror(istatus)
    stop
  end if
  istatus = nf90_inq_dimid(ncid, "iy", iydimid)
  if ( istatus /= nf90_noerr) then
    print *, nf90_strerror(istatus)
    stop
  end if
  istatus = nf90_inquire_dimension(ncid, iydimid, len=iy)
  if ( istatus /= nf90_noerr) then
    print *, nf90_strerror(istatus)
    stop
  end if

  allocate(xlat(jx,iy))
  allocate(xlon(jx,iy))
  allocate(landuse(jx,iy))
  allocate(mask(jx,iy))

  istatus = nf90_inq_varid(ncid, "xlat", ivarid)
  if ( istatus /= nf90_noerr) then
    print *, nf90_strerror(istatus)
    stop
  end if
  istatus = nf90_get_var(ncid, ivarid, xlat)
  if ( istatus /= nf90_noerr) then
    print *, nf90_strerror(istatus)
    stop
  end if

  istatus = nf90_inq_varid(ncid, "xlon", ivarid)
  if ( istatus /= nf90_noerr) then
    print *, nf90_strerror(istatus)
    stop
  end if
  istatus = nf90_get_var(ncid, ivarid, xlon)
  if ( istatus /= nf90_noerr) then
    print *, nf90_strerror(istatus)
    stop
  end if

  istatus = nf90_inq_varid(ncid, "landuse", ivarid)
  if ( istatus /= nf90_noerr) then
    print *, nf90_strerror(istatus)
    stop
  end if
  istatus = nf90_get_var(ncid, ivarid, landuse)
  if ( istatus /= nf90_noerr) then
    print *, nf90_strerror(istatus)
    stop
  end if

  mask(:,:) = .false.
  where ( (xlat > 45.0 .and. xlat < 48.0) .and. &
          (xlon > 13.0 .and. xlon < 15.0) )
    mask = .true.
  end where

  where ( mask .and. landuse == 5 )
    landuse = 1
  end where

  istatus = nf90_put_var(ncid, ivarid, landuse)
  if ( istatus /= nf90_noerr) then
    print *, nf90_strerror(istatus)
    stop
  end if

  istatus = nf90_close(ncid)
  if ( istatus /= nf90_noerr) then
    print *, nf90_strerror(istatus)
    stop
  end if

end program changeland

program changeland
  use netcdf
  implicit none
  integer :: ncid
  integer :: istatus
  character(len=256) :: arg1
  integer :: jx , iy
  integer :: jxdimid , iydimid , ivarid

  real(4) , pointer , dimension(:,:) :: xlat , xlon , var
  real(4) :: mean_elevation

  integer :: ip , i , j , ii , jj
  integer , parameter :: nlakes = 4
  integer , dimension(nlakes) :: jxlak
  integer , dimension(nlakes) :: iylak
  real(4) , dimension(nlakes) :: lak_dpth
  real(4) , dimension(nlakes) :: elevation

  data jxlak /161,156,158,155/
  data iylak /186,180,181,179/
  !data jxlak /13,12,13,12/
  !data iylak /18,18,19,19/
  data lak_dpth /310.0,145.0,82.0,71.0/
  data elevation /372.0,231.0,447.0,373.0/
  !data elevation /372.0,372.0,372.0,372.0/

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
  allocate(var(jx,iy))

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
  istatus = nf90_get_var(ncid, ivarid, var)
  if ( istatus /= nf90_noerr) then
    print *, nf90_strerror(istatus)
    stop
  end if
  do ip = 1 , nlakes
    var(jxlak(ip),iylak(ip)) = 14.0
  end do
  istatus = nf90_put_var(ncid, ivarid, var)
  if ( istatus /= nf90_noerr) then
    print *, nf90_strerror(istatus)
    stop
  end if

  istatus = nf90_inq_varid(ncid, "mask", ivarid)
  if ( istatus /= nf90_noerr) then
    print *, nf90_strerror(istatus)
    stop
  end if
  istatus = nf90_get_var(ncid, ivarid, var)
  if ( istatus /= nf90_noerr) then
    print *, nf90_strerror(istatus)
    stop
  end if
  do ip = 1 , nlakes
    var(jxlak(ip),iylak(ip)) = 0.0
  end do
  istatus = nf90_put_var(ncid, ivarid, var)
  if ( istatus /= nf90_noerr) then
    print *, nf90_strerror(istatus)
    stop
  end if

  istatus = nf90_inq_varid(ncid, "dhlake", ivarid)
  if ( istatus /= nf90_noerr) then
    print *, nf90_strerror(istatus)
    stop
  end if
  istatus = nf90_get_var(ncid, ivarid, var)
  if ( istatus /= nf90_noerr) then
    print *, nf90_strerror(istatus)
    stop
  end if
  do ip = 1 , nlakes
    var(jxlak(ip),iylak(ip)) = lak_dpth(ip)
  end do
  istatus = nf90_put_var(ncid, ivarid, var)
  if ( istatus /= nf90_noerr) then
    print *, nf90_strerror(istatus)
    stop
  end if

  istatus = nf90_inq_varid(ncid, "topo", ivarid)
  if ( istatus /= nf90_noerr) then
    print *, nf90_strerror(istatus)
    stop
  end if
  istatus = nf90_get_var(ncid, ivarid, var)
  if ( istatus /= nf90_noerr) then
    print *, nf90_strerror(istatus)
    stop
  end if
  do ip = 1 , nlakes
    var(jxlak(ip),iylak(ip)) = elevation(ip)
  end do
  do ip = 1 , nlakes
    ii = iylak(ip) - 3
    jj = jxlak(ip) - 3
    var(jxlak(ip),iylak(ip)) = elevation(ip)
    mean_elevation = elevation(ip) * 5.0
    do i = 1 , 5
      do j = 1 , 5
        mean_elevation = mean_elevation + var(jj+j,ii+i)
      end do
    end do
    mean_elevation = mean_elevation / 30.0
    do i = 1 , 5
      do j = 1 , 5
        var(jj+j,ii+i) = (var(jj+j,ii+i) + mean_elevation*2.0) / 3.0
      end do
    end do
  end do
  do ip = 1 , nlakes
    var(jxlak(ip),iylak(ip)) = elevation(ip)
  end do
  istatus = nf90_put_var(ncid, ivarid, var)
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

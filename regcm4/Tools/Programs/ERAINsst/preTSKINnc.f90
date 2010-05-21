      program fgennc
      include 'netcdf.inc'
! error status return
      integer  iret
! netCDF id
      integer  ncid
! dimension ids
      integer  longitude_dim
      integer  latitude_dim
      integer  time_dim
! dimension lengths
      integer  longitude_len
      integer  latitude_len
      integer  time_len
      parameter (longitude_len = 240)
      parameter (latitude_len = 121)
      parameter (time_len = NF_UNLIMITED)
! variable ids
      integer  longitude_id
      integer  latitude_id
      integer  time_id
      integer  skt_id
! rank (number of dimensions) for each variable
      integer  longitude_rank
      integer  latitude_rank
      integer  time_rank
      integer  skt_rank
      parameter (longitude_rank = 1)
      parameter (latitude_rank = 1)
      parameter (time_rank = 1)
      parameter (skt_rank = 3)
! variable shapes
      integer  longitude_dims(longitude_rank)
      integer  latitude_dims(latitude_rank)
      integer  time_dims(time_rank)
      integer  skt_dims(skt_rank)
! data variables
      real  longitude(longitude_len)
      real  latitude(latitude_len)
! attribute vectors
      integer  int2val(1)
      double precision  doubleval(1)
! enter define mode
      iret = nf_create('tskinERAIN.1989-2009.00.nc', NF_CLOBBER, ncid)
      call check_err(iret)
! define dimensions
      iret = nf_def_dim(ncid, 'longitude', longitude_len, longitude_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'latitude', latitude_len, latitude_dim)
      call check_err(iret)
      iret = nf_def_dim(ncid, 'time', NF_UNLIMITED, time_dim)
      call check_err(iret)
! define variables
      longitude_dims(1) = longitude_dim
      iret = nf_def_var(ncid, 'longitude', NF_REAL, longitude_rank, longitude_dims, longitude_id)
      call check_err(iret)
      latitude_dims(1) = latitude_dim
      iret = nf_def_var(ncid, 'latitude', NF_REAL, latitude_rank, latitude_dims, latitude_id)
      call check_err(iret)
      time_dims(1) = time_dim
      iret = nf_def_var(ncid, 'time', NF_INT, time_rank, time_dims, time_id)
      call check_err(iret)
      skt_dims(3) = time_dim
      skt_dims(2) = latitude_dim
      skt_dims(1) = longitude_dim
      iret = nf_def_var(ncid, 'skt', NF_INT2, skt_rank, skt_dims, skt_id)
      call check_err(iret)
! assign attributes
      iret = nf_put_att_text(ncid, longitude_id, 'units', 12, 'degrees_east')
      call check_err(iret)
      iret = nf_put_att_text(ncid, longitude_id, 'long_name', 9, 'longitude')
      call check_err(iret)
      iret = nf_put_att_text(ncid, latitude_id, 'units', 13, 'degrees_north')
      call check_err(iret)
      iret = nf_put_att_text(ncid, latitude_id, 'long_name', 8, 'latitude')
      call check_err(iret)
      iret = nf_put_att_text(ncid, time_id, 'units', 32, 'hours since 1900-01-01 00:00:0.0')
      call check_err(iret)
      iret = nf_put_att_text(ncid, time_id, 'long_name', 4, 'time')
      call check_err(iret)
      doubleval(1) = 0.00060759467941234814
      iret = nf_put_att_double(ncid, skt_id, 'scale_factor', NF_DOUBLE, 1, doubleval)
      call check_err(iret)
      doubleval(1) = 288.1350402832031
      iret = nf_put_att_double(ncid, skt_id, 'add_offset', NF_DOUBLE, 1, doubleval)
      call check_err(iret)
      int2val(1) = -32767
      iret = nf_put_att_int(ncid, skt_id, '_FillValue', NF_INT2, 1, int2val)
      call check_err(iret)
      int2val(1) = -32767
      iret = nf_put_att_int(ncid, skt_id, 'missing_value', NF_INT2, 1, int2val)
      call check_err(iret)
      iret = nf_put_att_text(ncid, skt_id, 'units', 1, 'K')
      call check_err(iret)
      iret = nf_put_att_text(ncid, skt_id, 'long_name', 16, 'Skin temperature')
      call check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'Conventions', 6, 'CF-1.0')
      call check_err(iret)
      iret = nf_put_att_text(ncid, NF_GLOBAL, 'history', 43, '2009-09-14 10:45:23 GMT by mars2netcdf-0.92')
      call check_err(iret)
! leave define mode
      iret = nf_enddef(ncid)
      call check_err(iret)
! store longitude
      do i=1,longitude_len
         longitude(i) = float(i-1)*1.5
      enddo
      iret = nf_put_var_real(ncid, longitude_id, longitude)
      call check_err(iret)
! store latitude
      do i=1,latitude_len
         latitude(i) = 90.-float(i-1)*1.5
      enddo
      iret = nf_put_var_real(ncid, latitude_id, latitude)
      call check_err(iret)
       
! Write record variables
      call writerecs(ncid,time_id,skt_id)
       
      iret = nf_close(ncid)
      call check_err(iret)
      end
       
      subroutine writerecs(ncid,time_id,skt_id)
       
! netCDF id
      integer  ncid
! variable ids
      integer  time_id
      integer  skt_id
       
      include 'netcdf.inc'
! error status return
      integer  iret
       
! netCDF dimension sizes for dimensions used with record variables
      integer  longitude_len
      parameter (longitude_len = 240)
      integer  latitude_len
      parameter (latitude_len = 121)
       
! rank (number of dimensions) for each variable
      integer  time_rank
      integer  skt_rank
      parameter (time_rank = 1)
      parameter (skt_rank = 3)
! starts and counts for array sections of record variables
      integer  time_start(time_rank), time_count(time_rank)
      integer  skt_start(skt_rank), skt_count(skt_rank)
       
! data variables
       
      integer  time_nr
      parameter (time_nr = 29824)
      integer  time(time_nr)
       
      integer  skt_nr
      parameter (skt_nr = 29824)
      integer*2  skt(longitude_len, latitude_len, skt_nr)
       
      open(21,file='tskinI2.dat',form='unformatted',recl=240*121*2,access='direct')
      do i=1,29824
         time(i) = i*6+780162
         read(21,rec=i) skt(:,:,i)
      enddo
       
! store time
      time_start(1) = 1
      time_count(1) = time_nr
      iret = nf_put_vara_int(ncid, time_id, time_start, time_count, time)
      call check_err(iret)
! store skt
      skt_start(1) = 1
      skt_start(2) = 1
      skt_start(3) = 1
      skt_count(1) = longitude_len
      skt_count(2) = latitude_len
      skt_count(3) = skt_nr
      iret = nf_put_vara_int2(ncid, skt_id, skt_start, skt_count, skt)
      call check_err(iret)
       
      end
       
      subroutine check_err(iret)
      integer iret
      include 'netcdf.inc'
      if (iret .ne. NF_NOERR) then
      print *, nf_strerror(iret)
      stop
      endif
      end

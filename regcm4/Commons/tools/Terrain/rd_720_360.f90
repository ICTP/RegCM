 !-------------------------------------------------------
 !  
 !  rd_netcdf_model.f90
 !    This file is a fortran template file designed to read the given
 !     netCDF file 'GTOPO30_30MIN.CDF' into memory.
 !  
 !  History: 
 !  Date       Name          Action
 !-------------------------------------------------------
 !  
 !  ?? Oct 93  B. Schwartz   Created.
 !  24 Aug 94  E. Boer       Modified dimension output.
 !  29 Jul 97  S. Rupert     Standardized and reinstated required 
 !                            include of netcdf.inc. 
 !  30 Jul 87  S. Rupert     Added usage message and command line inputs.
 !  3  Apr 2003 H. Yan       Change netcdf2 fortran code to netcdf3.5 fortran 90
 !                           code and modified attributes maxmum dimesions

 !-------------------------------------------------------
 !   Do not forget to include the -I path_to_netcdf_
 !   includes in your compile statement Required includes.
 !!! Also note: need 'netcdf.lib' or 'netcdfs.lib' when link 
       include 'netcdf.inc'

 !     Define Variables.
 !     Variable ids run sequentially from 1 to nvars=            5
 !     nvars=            5   ! number of variables
      integer,parameter :: nrec=            1   ! change this to generalize
      integer*4  ncid, status    ! file control
      integer*4 recdim   ! record dimension
 !-------------------------------------------------------------
 !     Below   5 variables is the data in netCDF file
      real*4         :: lon(          720 )
      real*4         :: lat(          360 )
      real*8         :: time(nrec)
      integer*2      ::  HT( 720, 360,nrec)
      integer*2      ::  HTSD( 720, 360,nrec)
      real*4         ::  a(720,360)
 !     above   5 variables is the data in netCDF file
 !-------------------------------------------------------------
      integer*4   :: start(10)
      integer*4   :: count(10)
      integer     :: dimids(10)! allow up to 10 dimensions
      integer     :: dimid, xtype
      character(len=31) :: dummy

 ! Open netCDF file.
       status=nf_open('GTOPO30_30MIN.CDF',nf_nowrite,ncid)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)


 !----------------------------------------------------
 !   Retrieve data for Variable 'lon'
      status=nf_inq_var(ncid,   1,dummy,xtype,ndim,dimids,natts)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_real(ncid,   1,start,count,lon)

 !----------------------------------------------------
 !   Retrieve data for Variable 'lat'
      status=nf_inq_var(ncid,   2,dummy,xtype,ndim,dimids,natts)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_real(ncid,   2,start,count,lat)

 !----------------------------------------------------
 !   Retrieve data for Variable 'time'
      status=nf_inq_var(ncid,   3,dummy,xtype,ndim,dimids,natts)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_double(ncid,   3,start,count,time)

 !----------------------------------------------------
 !   Retrieve data for Variable 'HT'
      status=nf_inq_var(ncid,   4,dummy,xtype,ndim,dimids,natts)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_int2(ncid,   4,start,count,HT)

 !----------------------------------------------------
 !   Retrieve data for Variable 'HTSD'
      status=nf_inq_var(ncid,   5,dummy,xtype,ndim,dimids,natts)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_int2(ncid,   5,start,count,HTSD)

 !----------------------------------------------------
 !  Begin writing statements to use the data.


 !     Here write your own code please!
       open(20,file='HT_pal.dat',form='unformatted',recl=720*360*4,access='direct')
       do j=1,360
       do i=1,720
          a(i,j) = HT(i,j,1)
       enddo
       enddo
       write(20,rec=1) a
       do j=1,360
       do i=1,720
          a(i,j) = HTSD(i,j,1)
       enddo
       enddo
       write(20,rec=2) a
       close(20)




 !----------------------------------------------------
 !  End Program

         end

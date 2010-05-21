 !-------------------------------------------------------
 !  
 !  rd_netcdf_model.f90
 !    This file is a fortran template file designed to read the given
 !     netCDF file 'tskin.1989-2009.12.nc' into memory.
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
 !     Variable ids run sequentially from 1 to nvars=            4
 !     nvars=            4   ! number of variables
      integer,parameter :: nrec=         7456   ! change this to generalize
      integer*4  ncid, status    ! file control
      integer*4 recdim   ! record dimension
 !-------------------------------------------------------------
 !     Below   4 variables is the data in netCDF file
      real*4         :: longitude(240)
      real*4         :: latitude(121)
      integer*4      :: time(nrec)
      integer*2      ::  skt(240,121,nrec)
      real*4         ::  tskin(240,121)
 !     above   4 variables is the data in netCDF file
 !-------------------------------------------------------------
      integer*4   :: start(10)
      integer*4   :: count(10)
      integer     :: dimids(10)! allow up to 10 dimensions
      integer     :: dimid, xtype
      character(len=31) :: dummy
      real*8  xscale,xadd
      integer :: i,j,k

 ! Open netCDF file.
       status=nf_open('/media/disk/RCM_DATA/ERAIN150/tskin.1989-2009.12.nc',nf_nowrite,ncid)
       xscale=0.0026407060784414
       xadd  =257.712584513648
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)


 !----------------------------------------------------
 !   Retrieve data for Variable 'longitude'
      status=nf_inq_var(ncid,   1,dummy,xtype,ndim,dimids,natts)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_real(ncid,   1,start,count,longitude)

 !----------------------------------------------------
 !   Retrieve data for Variable 'latitude'
      status=nf_inq_var(ncid,   2,dummy,xtype,ndim,dimids,natts)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_real(ncid,   2,start,count,latitude)

 !----------------------------------------------------
 !   Retrieve data for Variable 'time'
      status=nf_inq_var(ncid,   3,dummy,xtype,ndim,dimids,natts)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_int(ncid,   3,start,count,time)

 !----------------------------------------------------
 !   Retrieve data for Variable 'skt'
      status=nf_inq_var(ncid,   4,dummy,xtype,ndim,dimids,natts)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      do j=1,ndim
      status=nf_inq_dim(ncid,dimids(j),dummy,len)
           if ( status/=nf_noerr ) write (*,*) nf_strerror(status)
      start(j)=1 ; count(j)=len
      end do
      status=nf_get_vara_int2(ncid,   4,start,count,skt)

 !----------------------------------------------------
 !  Begin writing statements to use the data.
      open(20,file='tskin_12.dat',form='unformatted',recl=240*121*4,access='direct')
      do n=1,nrec
         do j=1,121
         do i=1,240
            if(skt(i,j,n).eq.-32767) then
               tskin(i,j) = -9999.0
            else
               tskin(i,j) = skt(i,j,n)*xscale+xadd
            endif
         enddo
         enddo
         write(20,rec=n) tskin
      enddo
      close(20)

 !     Here write your own code please!



 !----------------------------------------------------
 !  End Program

         end


program trend
  use netcdf

  implicit none

  ! The variable to search
  character(len=16) , dimension(2) , parameter :: timenames = &
                 (/'time','Time'/)
  character(len=16) , dimension(5) , parameter :: latnames = &
                 (/'xlat    ','lat     ','latitude','LAT     ','iy      '/)
  character(len=16) , dimension(5) , parameter :: lonnames = &
                 (/'xlon     ','lon      ','longitude','LON      ','jx       '/)

  character(len=256) :: progname
  character(len=256) :: inputfile
  character(len=256) :: outputfile
  character(len=32) :: varname
  character(len=32) :: dname , vname , aname
  integer :: ncid , ncout , ixdimid , iydimid , itimid
  integer :: idlat , idlon , idlatout , idlonout , ndims
  integer :: istatus
  integer , dimension(3) :: odims , istart , icount
  integer :: nvar , ivar , natt , iatt , idvar , idtime
  integer :: ivarouta , ivaroutb , ivaroutc
  integer :: i , nx , ny , nt , ipoint , npoint , it
  logical :: coord_2d = .false. , has_coord = .false.
  real , allocatable , dimension(:,:) :: abetrend , reshaped , coord
  real , allocatable , dimension(:) :: coord1d , a , b , r , atime , x , y
  real , allocatable , dimension(:,:,:) :: inpindx

  call get_command_argument(0,value=progname)
  call get_command_argument(1,value=inputfile)
  call get_command_argument(2,value=outputfile)
  call get_command_argument(3,value=varname)

  istatus = nf90_open(inputfile,nf90_nowrite,ncid)
  if ( istatus /= nf90_noerr ) then
    write(0,*) 'Cannot open input file : ',trim(inputfile)
    stop
  end if
  istatus = nf90_create(outputfile,nf90_clobber,ncout)
  if ( istatus /= nf90_noerr ) then
    write(0,*) 'Cannot open output file : ',trim(outputfile)
    stop
  end if

  istatus = nf90_inq_dimid(ncid,'jx',ixdimid)
  if ( istatus /= nf90_noerr ) then
    istatus = nf90_inq_dimid(ncid,'x',ixdimid)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_dimid(ncid,'lon',ixdimid)
      if ( istatus /= nf90_noerr ) then
        istatus = nf90_inq_dimid(ncid,'LON',ixdimid)
        if ( istatus /= nf90_noerr ) then
          istatus = nf90_inq_dimid(ncid,'longitude',ixdimid)
          if ( istatus /= nf90_noerr ) then
            ! Give up
            write(0,*) 'Cannot find any coded WE dimension in file'
            stop
          end if
        end if
      end if
    end if
  end if
  istatus = nf90_inq_dimid(ncid,'iy',iydimid)
  if ( istatus /= nf90_noerr ) then
    istatus = nf90_inq_dimid(ncid,'y',iydimid)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_dimid(ncid,'lat',iydimid)
      if ( istatus /= nf90_noerr ) then
        istatus = nf90_inq_dimid(ncid,'LAT',iydimid)
        if ( istatus /= nf90_noerr ) then
          istatus = nf90_inq_dimid(ncid,'latitude',iydimid)
          if ( istatus /= nf90_noerr ) then
            ! Give up
            write(0,*) 'Cannot find any coded SN dimension in file'
            stop
          end if
        end if
      end if
    end if
  end if

  istatus = nf90_inquire_dimension(ncid,ixdimid,name=dname,len=nx)
  if ( istatus /= nf90_noerr ) then
    write(0,*) 'Cannot read dimension WE from file'
    stop
  end if
  istatus = nf90_def_dim(ncout,dname,nx,odims(1))
  if ( istatus /= nf90_noerr ) then
    write(0,*) 'Cannot copy dimension WE to outfile'
    stop
  end if
  istatus = nf90_inquire_dimension(ncid,iydimid,name=dname,len=ny)
  if ( istatus /= nf90_noerr ) then
    write(0,*) 'Cannot read dimension SN from file'
    stop
  end if
  istatus = nf90_def_dim(ncout,dname,ny,odims(2))
  if ( istatus /= nf90_noerr ) then
    write(0,*) 'Cannot copy dimension SN to outfile'
    stop
  end if

  istatus = nf90_inq_dimid(ncid,'time',itimid)
  if ( istatus /= nf90_noerr ) then
    istatus = nf90_inq_dimid(ncid,'TIME',itimid)
    if ( istatus /= nf90_noerr ) then
      write(0,*) 'Cannot find any coded time dimension in file'
      stop
    end if
  end if
  istatus = nf90_inquire_dimension(ncid,itimid,name=dname,len=nt)
  if ( istatus /= nf90_noerr ) then
    write(0,*) 'Cannot read dimension time from file'
    stop
  end if

  istatus = nf90_inquire(ncid,nVariables=nvar)
  if ( istatus /= nf90_noerr ) then
    write(0,*) 'Cannot inquire input file : ',trim(inputfile)
    stop
  end if

  idvar = -1
  varloop: &
  do ivar = 1 , nvar
    istatus = nf90_inquire_variable(ncid,ivar,name=vname,natts=natt,ndims=ndims)
    if ( istatus /= nf90_noerr ) then
      write(0,*) 'Cannot inquire variable ',ivar
      stop
    end if
    if ( vname == varname ) then
      idvar = ivar
    end if
    do i = 1 , size(latnames)
      if ( vname == latnames(i) ) then
        idlat = ivar
        if ( ndims == 2 ) then
          coord_2d = .true.
          istatus = nf90_def_var(ncout,vname,nf90_real,odims(1:2),idlatout)
          if ( istatus /= nf90_noerr ) then
            write(0,*) 'Cannot define variable '//trim(vname)
            stop
          end if
        else
          istatus = nf90_def_var(ncout,vname,nf90_real,odims(2:2),idlatout)
          if ( istatus /= nf90_noerr ) then
            write(0,*) 'Cannot define variable '//trim(vname)
            stop
          end if
        end if
        do iatt = 1 , natt
          istatus = nf90_inq_attname(ncid,ivar,iatt,aname)
          if ( istatus /= nf90_noerr ) then
            write(0,*) 'Cannot inquire attribute ',iatt
            stop
          end if
          istatus = nf90_copy_att(ncid,ivar,aname,ncout,idlatout)
          if ( istatus /= nf90_noerr ) then
            write(0,*) 'Cannot copy attribute '//trim(aname)
            stop
          end if
        end do
        has_coord = .true.
        cycle varloop
      end if
    end do
    do i = 1 , size(lonnames)
      if ( vname == lonnames(i) ) then
        idlon = ivar
        if ( ndims == 2 ) then
          coord_2d = .true.
          istatus = nf90_def_var(ncout,vname,nf90_real,odims(1:2),idlonout)
          if ( istatus /= nf90_noerr ) then
            write(0,*) 'Cannot define variable '//trim(vname)
            stop
          end if
        else
          istatus = nf90_def_var(ncout,vname,nf90_real,odims(1:1),idlonout)
          if ( istatus /= nf90_noerr ) then
            write(0,*) 'Cannot define variable '//trim(vname)
            stop
          end if
        end if
        do iatt = 1 , natt
          istatus = nf90_inq_attname(ncid,ivar,iatt,aname)
          if ( istatus /= nf90_noerr ) then
            write(0,*) 'Cannot inquire attribute ',iatt
            stop
          end if
          istatus = nf90_copy_att(ncid,ivar,aname,ncout,idlonout)
          if ( istatus /= nf90_noerr ) then
            write(0,*) 'Cannot copy attribute '//trim(aname)
            stop
          end if
        end do
        has_coord = .true.
        cycle varloop
      end if
    end do
    do i = 1 , size(timenames)
      if ( vname == timenames(i) ) then
        idtime = ivar
      end if
      cycle varloop
    end do
  end do varloop

  if ( idvar == -1 ) then
    write(0,*) 'I have not found the variable ',trim(varname), &
               ' in the input file.'
    write(0,*) 'Usage : ',trim(progname),' infile.nc outfile.nc varname'
    stop
  end if

  istatus = nf90_def_var(ncout,'trend_a',nf90_real,odims(1:2),ivarouta)
  if ( istatus /= nf90_noerr ) then
    write(0,*) 'Cannot define variable trend'
    stop
  end if
  istatus = nf90_put_att(ncout,ivarouta,'long_name', &
               'Trend for var '//trim(varname)//' a coefficient in var=a+bt')
  if ( istatus /= nf90_noerr ) then
    write(0,*) 'Cannot add atrribute to trend variable'
    stop
  end if
  istatus = nf90_copy_att(ncid,idvar,'units',ncout,ivarouta)
  if ( istatus /= nf90_noerr ) then
    write(0,*) 'Cannot copy attribute units to trend variable'
    stop
  end if
  istatus = nf90_def_var(ncout,'trend_b',nf90_real,odims(1:2),ivaroutb)
  if ( istatus /= nf90_noerr ) then
    write(0,*) 'Cannot define variable trend'
    stop
  end if
  istatus = nf90_put_att(ncout,ivaroutb,'long_name', &
               'Trend for var '//trim(varname)//' b coefficient in var=a+bt')
  if ( istatus /= nf90_noerr ) then
    write(0,*) 'Cannot add atrribute to trend variable'
    stop
  end if
  istatus = nf90_put_att(ncout,ivaroutb,'units','year-1')
  if ( istatus /= nf90_noerr ) then
    write(0,*) 'Cannot add atrribute to trend variable'
    stop
  end if
  istatus = nf90_def_var(ncout,'trend_r',nf90_real,odims(1:2),ivaroutc)
  if ( istatus /= nf90_noerr ) then
    write(0,*) 'Cannot define variable trend'
    stop
  end if
  istatus = nf90_put_att(ncout,ivaroutc,'long_name', &
               'Trend correlation coefficient for var '//trim(varname))
  if ( istatus /= nf90_noerr ) then
    write(0,*) 'Cannot add atrribute to trend variable'
    stop
  end if
  istatus = nf90_put_att(ncout,ivaroutc,'units','1')
  if ( istatus /= nf90_noerr ) then
    write(0,*) 'Cannot add atrribute to trend variable'
    stop
  end if
  if ( has_coord .and. coord_2d ) then
    istatus = nf90_put_att(ncout,ivarouta,'coordinates','xlat xlon')
    if ( istatus /= nf90_noerr ) then
      write(0,*) 'Cannot add atrribute to trend variable'
      stop
    end if
    istatus = nf90_put_att(ncout,ivaroutb,'coordinates','xlat xlon')
    if ( istatus /= nf90_noerr ) then
      write(0,*) 'Cannot add atrribute to trend variable'
      stop
    end if
    istatus = nf90_put_att(ncout,ivaroutc,'coordinates','xlat xlon')
    if ( istatus /= nf90_noerr ) then
      write(0,*) 'Cannot add atrribute to trend variable'
      stop
    end if
  end if

  istatus = nf90_enddef(ncout)
  if ( istatus /= nf90_noerr ) then
    write(0,*) 'Cannot finalixe output file'
    stop
  end if

  if ( has_coord ) then
    if ( coord_2d ) then
      allocate(coord(nx,ny))
      istatus = nf90_get_var(ncid,idlat,coord)
      if ( istatus /= nf90_noerr ) then
        write(0,*) 'Cannot read lat coord'
        stop
      end if
      istatus = nf90_put_var(ncout,idlatout,coord)
      if ( istatus /= nf90_noerr ) then
        write(0,*) 'Cannot write lat coord'
        stop
      end if
      istatus = nf90_get_var(ncid,idlon,coord)
      if ( istatus /= nf90_noerr ) then
        write(0,*) 'Cannot read lon coord'
        stop
      end if
      istatus = nf90_put_var(ncout,idlonout,coord)
      if ( istatus /= nf90_noerr ) then
        write(0,*) 'Cannot write lon coord'
        stop
      end if
      deallocate(coord)
    else
      allocate(coord1d(ny))
      istatus = nf90_get_var(ncid,idlat,coord1d)
      if ( istatus /= nf90_noerr ) then
        write(0,*) 'Cannot read lat coord'
        stop
      end if
      istatus = nf90_put_var(ncout,idlatout,coord1d)
      if ( istatus /= nf90_noerr ) then
        write(0,*) 'Cannot write lat coord'
        stop
      end if
      deallocate(coord1d)
      allocate(coord1d(nx))
      istatus = nf90_get_var(ncid,idlon,coord1d)
      if ( istatus /= nf90_noerr ) then
        write(0,*) 'Cannot read lon coord'
        stop
      end if
      istatus = nf90_put_var(ncout,idlonout,coord1d)
      if ( istatus /= nf90_noerr ) then
        write(0,*) 'Cannot write lon coord'
        stop
      end if
      deallocate(coord1d)
    end if
  end if

  npoint = nx*ny
  allocate(inpindx(nx,ny,nt))
  allocate(reshaped(nt,npoint))
  allocate(atime(nt))
  allocate(x(nt))
  allocate(y(nt))
  allocate(a(npoint))
  allocate(b(npoint))
  allocate(r(npoint))
  allocate(abetrend(nx,ny))

  write (6,*) 'Read ',nt,' ',trim(varname),' steps'
  istatus = nf90_get_var(ncid,idvar,inpindx)
  if ( istatus /= nf90_noerr ) then
    write(0,*) 'Error reading ',nt,' steps'
    stop
  end if

  do it = 1 , nt
    atime(it) = real(it)
    reshaped(it,:) = reshape(inpindx(:,:,it),(/npoint/))
  end do

  do it = 1 , nt
    x(it) = atime(it)
  end do
  do ipoint = 1 , npoint
    do it = 1 , nt
      y(it) = reshaped(it,ipoint)
    end do
    call linear_regression_coefficients(x,y,nt,a(ipoint),b(ipoint),r(ipoint))
  end do

  istart(1) = 1
  istart(2) = 1
  icount(1) = nx
  icount(2) = ny
  abetrend(:,:) = reshape(a,(/nx,ny/))
  istatus = nf90_put_var(ncout,ivarouta,abetrend,istart,icount)
  if ( istatus /= nf90_noerr ) then
    write(0,*) 'Cannot write duration variable'
    stop
  end if
  abetrend(:,:) = reshape(b,(/nx,ny/))
  istatus = nf90_put_var(ncout,ivaroutb,abetrend,istart,icount)
  if ( istatus /= nf90_noerr ) then
    write(0,*) 'Cannot write duration variable'
    stop
  end if
  abetrend(:,:) = reshape(r,(/nx,ny/))
  istatus = nf90_put_var(ncout,ivaroutc,abetrend,istart,icount)
  if ( istatus /= nf90_noerr ) then
    write(0,*) 'Cannot write duration variable'
    stop
  end if

  istatus = nf90_close(ncid)
  if ( istatus /= nf90_noerr ) then
    write(0,*) 'Cannot close input file : ',trim(inputfile)
    stop
  end if
  istatus = nf90_close(ncout)
  if ( istatus /= nf90_noerr ) then
    write(0,*) 'Cannot close output file : ',trim(outputfile)
    stop
  end if

  write(6,*) 'Done.'

  contains

   !
   ! Coefficients of linear regression
   !
   subroutine linear_regression_coefficients(x,y,n,a,b,r)
     implicit none
     real , dimension(:) , intent(in) :: x , y
     integer , intent(in) :: n
     real , intent(out) :: a , b , r
     real(8) , allocatable , dimension(:) :: yd , xd
     real(8) :: nn , sx , sy , sxx , syy , sxy , alpha , beta , corr
     allocate(xd(n),yd(n))
     xd = x(1:n)
     yd = y(1:n)
     nn = dble(n)
     sx = sum(xd)
     sy = sum(yd)
     sxx = sum(xd*xd)
     syy = sum(yd*yd)
     sxy = sum(xd*yd)
     beta = (nn*sxy-sx*sy)/(nn*sxx-sx*sx)
     alpha = sy/nn - beta*sx/nn
     corr = (nn*sxy-sx*sy)/sqrt((nn*sxx-sx*sx)*(nn*syy-sy*sy))
     a = real(alpha)
     b = real(beta)
     r = real(corr)
     deallocate(xd,yd)
   end subroutine linear_regression_coefficients

end program trend
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

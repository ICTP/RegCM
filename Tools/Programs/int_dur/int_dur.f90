
program int_dur
  use netcdf
  use mod_date

  implicit none

  ! The variable to search
  character(len=16) , dimension(2) , parameter :: prenames =  &
                 (/'pre ','pre2'/)
  character(len=16) , dimension(2) , parameter :: timenames = &
                 (/'time','Time'/)
  character(len=16) , dimension(5) , parameter :: latnames = &
                 (/'xlat    ','lat     ','latitude','LAT     ','iy      '/)
  character(len=16) , dimension(5) , parameter :: lonnames = &
                 (/'xlon     ','lon      ','longitude','LON      ','jx       '/)

  character(len=256) :: inputfile
  character(len=256) :: outputfile
  character(len=32) :: dname , vname , aname
  character(len=64) :: time_unit , time_calendar, time_unitout
  integer :: ncid , ncout , ixdimid , iydimid , itimid
  integer :: idlat , idlon , idlatout , idlonout , ndims
  integer :: istatus
  integer , dimension(3) :: odims , istart , icount
  integer :: nvar , ivar , natt , iatt , idpre , idtime
  integer :: idpreout , iddurout , idtimeout
  integer :: i , it , iy , jx , nx , ny , nt
  logical :: coord_2d = .false. , has_coord = .false.
  real , allocatable , dimension(:,:) :: meanpre , duration , coord
  real , allocatable , dimension(:) :: coord1d
  real , allocatable , dimension(:,:,:) :: preslice
  real(8) :: pcount , sumpre , idursum , inumsum , percent , lastpercent
  real(8) , dimension(1) :: xtime
  integer :: iyear , nyear , ndays , nprocessed , maxdays
  type(rcm_time_and_date) :: idate1 , idate2 , idatecheck
  type(rcm_time_interval) :: tdif
  integer :: year1 , month1 , day1 , hour1
  integer :: year2 , month2 , day2 , hour2
  integer :: maxdpy , daypy , dstart

  call get_command_argument(1,value=inputfile)
  call get_command_argument(2,value=outputfile)

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
  istatus = nf90_def_dim(ncout,dname,nf90_unlimited,odims(3))
  if ( istatus /= nf90_noerr ) then
    write(0,*) 'Cannot copy dimension time to outfile'
    stop
  end if

  istatus = nf90_inquire(ncid,nVariables=nvar)
  if ( istatus /= nf90_noerr ) then
    write(0,*) 'Cannot inquire input file : ',trim(inputfile)
    stop
  end if

  varloop: &
  do ivar = 1 , nvar
    istatus = nf90_inquire_variable(ncid,ivar,name=vname,natts=natt,ndims=ndims)
    if ( istatus /= nf90_noerr ) then
      write(0,*) 'Cannot inquire variable ',ivar
      stop
    end if
    do i = 1 , size(prenames)
      if ( vname == prenames(i) ) then
        idpre = ivar
        istatus = nf90_def_var(ncout,vname,nf90_real,odims,idpreout)
        if ( istatus /= nf90_noerr ) then
          write(0,*) 'Cannot define variable '//trim(vname)
          stop
        end if
        do iatt = 1 , natt
          istatus = nf90_inq_attname(ncid,ivar,iatt,aname)
          if ( istatus /= nf90_noerr ) then
            write(0,*) 'Cannot inquire attribute ',iatt
            stop
          end if
          istatus = nf90_copy_att(ncid,ivar,aname,ncout,idpreout)
          if ( istatus /= nf90_noerr ) then
            write(0,*) 'Cannot copy attribute '//trim(aname)
            stop
          end if
        end do
        cycle varloop
      end if
    end do
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
        istatus = nf90_def_var(ncout,vname,nf90_double,odims(3:3),idtimeout)
        if ( istatus /= nf90_noerr ) then
          write(0,*) 'Cannot define variable time'
          stop
        end if
        att_time_loop: &
        do iatt = 1 , natt
          istatus = nf90_inq_attname(ncid,ivar,iatt,aname)
          if ( istatus /= nf90_noerr ) then
            write(0,*) 'Cannot inquire attribute ',iatt
            stop
          end if
          if ( aname == 'units' ) then
            istatus = nf90_get_att(ncid,ivar,'units',time_unit)
            if ( istatus /= nf90_noerr ) then
              write(0,*) 'Cannot read time units'
              stop
            end if
            cycle att_time_loop
          else if ( aname == 'calendar' ) then
            istatus = nf90_get_att(ncid,ivar,'calendar',time_calendar)
            if ( istatus /= nf90_noerr ) then
              write(0,*) 'Cannot read time calendar'
              stop
            end if
          end if
          if ( aname == 'bounds' ) cycle att_time_loop
          istatus = nf90_copy_att(ncid,ivar,aname,ncout,idtimeout)
          if ( istatus /= nf90_noerr ) then
            write(0,*) 'Cannot copy attribute '//trim(aname)
            stop
          end if
        end do att_time_loop
      end if
      cycle varloop
    end do
  end do varloop

  istatus = nf90_def_var(ncout,'duration',nf90_real,odims,iddurout)
  if ( istatus /= nf90_noerr ) then
    write(0,*) 'Cannot define variable duration'
    stop
  end if
  istatus = nf90_put_att(ncout,iddurout,'long_name','Rain period duration')
  if ( istatus /= nf90_noerr ) then
    write(0,*) 'Cannot add atrribute to duration variable'
    stop
  end if
  istatus = nf90_put_att(ncout,iddurout,'units','days')
  if ( istatus /= nf90_noerr ) then
    write(0,*) 'Cannot add atrribute to duration variable'
    stop
  end if
  istatus = nf90_put_att(ncout,iddurout,'coordinates','xlat xlon')
  if ( istatus /= nf90_noerr ) then
    write(0,*) 'Cannot add atrribute to duration variable'
    stop
  end if

  ! Calculate how many years we have to proces
  
  istart(1) = 1
  icount(1) = 1
  istatus = nf90_get_var(ncid,idtime,xtime,istart(1:1),icount(1:1))
  if ( istatus /= nf90_noerr ) then
    write(0,*) 'Cannot read time variable'
    stop
  end if
  idate1 = timeval2date(xtime(1), time_unit, time_calendar)
  call split_idate(idate1,year1,month1,day1,hour1)
  istart(1) = nt
  icount(1) = 1
  istatus = nf90_get_var(ncid,idtime,xtime,istart(1:1),icount(1:1))
  if ( istatus /= nf90_noerr ) then
    write(0,*) 'Cannot read time variable'
    stop
  end if
  idate2 = timeval2date(xtime(1), time_unit, time_calendar)
  call split_idate(idate2,year2,month2,day2,hour2)
 
  nyear = year2 - year1 + 1
  write(6,*) 'I have a total of ',nyear,' years in input file.'
  select case (idate1%calendar)
    case (noleap)
      maxdpy = 365+1 ! This for duration at 31 dec.
      daypy = 365
    case (y360)
      maxdpy = 360+1
      daypy = 360
    case default
      maxdpy = 366+1
      daypy = -1
  end select

  write(time_unitout,'(a,i4,a)') 'years since ',year1,'-06-15'

  istatus = nf90_put_att(ncout,idtimeout,'units',time_unitout)
  if ( istatus /= nf90_noerr ) then
    write(0,*) 'Cannot add time units'
    stop
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

  allocate(meanpre(nx,ny))
  allocate(duration(nx,ny))
  allocate(preslice(nx,ny,maxdpy))
 
  ! Check first date is Jan 01
  dstart = 1
  if ( month1 /= 1 .and. day1 /= 1 ) then
    idatecheck = year1*1000000+10112
    call setcal(idatecheck,idate1)
    tdif = idate1-idatecheck
    dstart = int(tohours(tdif)/24.0D0)
  end if

  nprocessed = 1
  do iyear = year1 , year2
    lastpercent = 0.0D0
    if ( idate1%calendar == gregorian ) then
      daypy = 365
      if ( lleap(iyear) ) daypy = 366
    end if
    ndays = daypy-dstart+1+1
    if ( nprocessed + ndays > nt ) ndays = nt-nprocessed+1
    istart(1) = 1
    istart(2) = 1
    istart(3) = nprocessed
    icount(1) = nx
    icount(2) = ny
    icount(3) = ndays
    write (6,*) 'Read ',ndays,' pre steps from ',nprocessed
    istatus = nf90_get_var(ncid,idpre,preslice,istart,icount)
    if ( istatus /= nf90_noerr ) then
      write(0,*) 'Error reading ',ndays,' pre steps from ',nprocessed
      stop
    end if
    meanpre = 0.0
    maxdays = ndays-1
    if ( iyear == year2 ) maxdays = ndays
    duration = real(maxdays)
    write(6,*) 'Data read for year :',iyear,', starting computation.'
    do iy = 1 , ny
      do jx = 1 , nx
        pcount = 0.0D0
        sumpre = 0.0D0
        do it = 1 , maxdays
          if ( preslice(jx,iy,it) >= 1.0 ) then
            sumpre = sumpre + dble(preslice(jx,iy,it))
            pcount = pcount + 1.0D0
          end if
        end do
        if ( pcount > 0 ) then
          meanpre(jx,iy) = real(sumpre/pcount)
        end if
        idursum = 0.0D0
        inumsum = 0.0D0
        do it = 1 , maxdays
          if ( preslice(jx,iy,it) < 1.0 ) then
            idursum = idursum + 1.0D0
          end if
        end do
        do it = 1 , maxdays
          if ( preslice(jx,iy,it) < 1.0 .and. preslice(jx,iy,it+1) >= 1.0 ) then
            inumsum = inumsum + 1.0D0
          end if
        end do
        if ( inumsum > 0.0 ) then
          duration(jx,iy) = real(idursum/inumsum)
        end if
        percent = 100.0D0*(dble((iy-1)*nx+jx)/dble(nx*ny))
        if ( percent > lastpercent+9.999D0 ) then
          write(6,*) 'Point ',jx , iy,' done ',int(percent),'%'
          lastpercent = percent
        end if
      end do
    end do
    write(6,*) 'Year processing completed.'

    istart(1) = iyear-year1+1
    icount(1) = 1
    xtime(1) = istart(1)-1
    istatus = nf90_put_var(ncout,idtimeout,xtime,istart(1:1),icount(1:1))
    if ( istatus /= nf90_noerr ) then
      write(0,*) 'Cannot write time variable'
      stop
    end if
    istart(1) = 1
    istart(2) = 1
    istart(3) = iyear-year1+1
    icount(1) = nx
    icount(2) = ny
    icount(3) = 1
    istatus = nf90_put_var(ncout,idpreout,meanpre,istart,icount)
    if ( istatus /= nf90_noerr ) then
      write(0,*) 'Cannot write pre variable'
      stop
    end if
    istatus = nf90_put_var(ncout,iddurout,duration,istart,icount)
    if ( istatus /= nf90_noerr ) then
      write(0,*) 'Cannot write duration variable'
      stop
    end if
    dstart = 1
    nprocessed = nprocessed + ndays - 2
  end do

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

end program int_dur

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

module mod_sst_1deg

  use mod_intkinds
  use mod_realkinds
  use mod_memutil
  use mod_stdio
  use mod_dynparam
  use mod_sst_grid
  use mod_interp
  use mod_nchelper
  use mod_message
  use netcdf

  private

  public :: sst_1deg

  integer(ik4) :: year , month , day , hour

  integer(ik4) :: istatus , latid , lonid
  integer(ik4) :: ilon , jlat
  real(rkx) , pointer , dimension(:) :: lati
  real(rkx) , pointer , dimension(:) :: loni
  real(rkx) , pointer , dimension(:,:) :: sst , ice
  integer(2) , pointer , dimension(:,:) :: work , work1

  contains
  !
  !
  ! Comments on dataset sources and location:
  !
  ! GISST2.3b    UKMO SST (Rayner et al 1996), 1 degree
  !              from UKMO DATA archive (http://www.badc.rl.ac.uk/)
  !              and reformed as direct-accessed binary GrADS format
  !              in file GISST_187101_200209
  !          ML = 1 is-179.5; ML = 2 is-178.5; => ML = 360 is 179.5E
  !          NL = 1 is -89.5; NL = 2 is -88.5; => NL = 180 is  89.5
  !              see the GrADS control file for details.
  !
  ! OISST        from CAC Optimal Interpolation dataset.
  !              in the original netCDF format.
  !          ML = 1 is   0.5; ML = 2 is   1.5; => ML = 360 is 359.5E
  !          NL = 1 is -89.5; NL = 2 is -88.5; => NL = 180 is  89.5
  !
  ! OI2ST        both SST and SeaIce in the original netCDF format.
  !
  ! OI_WK        weekly OISST in the original netCDF format.
  !
  ! OI2WK        weekly OISST and SeaIce in the original netCDF format.
  !
  subroutine sst_1deg
    implicit none
    real(rk4) , dimension(360,180) :: gisst
    integer(ik4) :: i , j , k , iwk , nrec
    integer(ik4) :: nsteps
    integer :: gireclen
    type(rcm_time_and_date) :: idate , idateo , idatef , idatem , irefd
    character(len=256) :: inpfile
    logical :: there

    if ( ssttyp == 'GISST' ) then
      ilon = 360
      jlat = 180
      call getmem1d(lati,1,jlat,'mod_sst_1deg:lati')
      call getmem1d(loni,1,ilon,'mod_sst_1deg:loni')
      call getmem2d(sst,1,ilon,1,jlat,'mod_sst_1deg:sst')
      if ( globidate1 < 1947121512 .or. globidate2 > 2002091512 ) then
        write (stderr,*) 'GISST data required are not available'
        write (stderr,*) 'IDATE1, IDATE2 = ' , globidate1 , globidate2
        call die('sst_1deg')
      end if
      ! SET UP LONGITUDES AND LATITUDES FOR SST DATA
      do i = 1 , ilon
        loni(i) = .5 + float(i-1)
      end do
      do j = 1 , jlat
        lati(j) = -89.5 + 1.*float(j-1)
      end do
      inpfile = trim(inpglob)//'/SST/GISST_194712_200209'
      inquire (file=inpfile,exist=there)
      if ( .not. there ) then
        call die('sst_1deg','GISST_194712_200209 is not available'//  &
                 ' under '//trim(inpglob)//'/SST/',1)
      end if
      inquire(iolength=gireclen) gisst
      write (stdout,*) 'OPEN ', trim(inpfile)
      open(121,file=inpfile,form='unformatted', &
           access='direct',recl=gireclen,action='read',status='old')
    else if ( ssttyp == 'OISST' .or. ssttyp == 'OI_NC' .or. &
              ssttyp == 'OI2ST' ) then
      if ( globidate1 < 1981121512 .or. globidate2 < 1981121512 ) then
        write (stderr,*) 'OISST data required are not available'
        write (stderr,*) 'IDATE1, IDATE2 = ' , globidate1 , globidate2
        call die('sst_1deg')
      end if
    else if ( ssttyp == 'OI_WK' .or. ssttyp == 'OI2WK' ) then
      if ( globidate1 < 1981110100 .or. globidate2 < 1981110106 ) then
        write (stderr,*) 'OI_WK (or OI2WK) data are not available'
        write (stderr,*) 'IDATE1, IDATE2 = ' , globidate1 , globidate2
        call die('sst_1deg')
      end if
    else
      write (stderr,*) 'PLEASE SET right SSTTYP in regcm.in'
      write (stderr,*) 'Supported are GISST OISST OI_NC OI2ST OI_WK OI2WK'
      call die('sst_1deg')
    end if

    ! Montly dataset
    if ( ssttyp /= 'OI_WK' .and. ssttyp /= 'OI2WK' ) then
      idateo = monfirst(globidate1)
      if (lfhomonth(globidate1)) then
        idateo = prevmon(globidate1)
      end if
      idatef = monfirst(globidate2)
      if (idatef < globidate2) then
        idatef = nextmon(idatef)
      end if
      nsteps = imondiff(idatef,idateo) + 1
      idatem = monmiddle(idateo)
      call open_sstfile(idatem)
    else
      ! Weekly dataset
      idateo = ifdoweek(globidate1)
      idatef = ildoweek(globidate2)
      nsteps = iwkdiff(idatef,idateo) + 1
      call open_sstfile(idateo)
    end if

    write (stdout,*) 'NSTEPS     : ' , nsteps

    ! ****** OISST SST DATA, 1 Deg data, AVAILABLE FROM 12/1981 TO PRESENT
    ! ****** GISST SST DATA, 1 Deg data, AVAILABLE FROM 12/1947 TO 9/2002

    idate = idateo

    if ( ssttyp /= 'OI_WK' .and. ssttyp /= 'OI2WK' ) then

      do k = 1 , nsteps

        call split_idate(idate,year,month,day,hour)

        if ( ssttyp == 'GISST' ) then
          nrec = (year-1947)*12 + month - 11
          read (121,rec=nrec) gisst
          sst = gisst
        else if ( ssttyp == 'OISST' .or. ssttyp == 'OI_NC' .or. &
                  ssttyp == 'OI2ST') then
          inpfile=trim(inpglob)//'/SST/sst.mnmean.nc'
          call sst_mn(idate,idateo,inpfile)
          if ( ssttyp == 'OI2ST' ) then
            inpfile=trim(inpglob)//'/SST/icec.mnmean.nc'
            call ice_mn(idate,idateo,inpfile)
          end if
        end if

        call bilinx(sstmm,sst,xlon,xlat,loni,lati,ilon,jlat,jx,iy)
        if ( ssttyp == 'OI2ST' ) then
          call bilinx(icemm,ice,xlon,xlat,loni,lati,ilon,jlat,jx,iy)
        end if

        write (stdout,*) 'XLON,XLAT,SST = ', xlon(1,1), xlat(1,1), sstmm(1,1)

        do j = 1 , jx
          do i = 1 , iy
            if ( sstmm(j,i) > -100. ) then
              sstmm(j,i) = sstmm(j,i) + 273.15
            else
              sstmm(j,i) = -9999.
            end if
          end do
        end do

        call writerec(idate)

        write (stdout,*) 'WRITTEN OUT SST DATA : ' , tochar(idate)

        idate = nextmon(idate)
        idatem = monmiddle(idate)

      end do

    else
      ! Weekly data

      do k = 1 , nsteps

        if ( idate < 1989123100 ) then
          irefd = 1981110100
          inpfile=trim(inpglob)//'/SST/sst.wkmean.1981-1989.nc'
          iwk = iwkdiff(idate,irefd)
        else
          irefd = 1989123100
          inpfile=trim(inpglob)//'/SST/sst.wkmean.1990-present.nc'
          iwk = iwkdiff(idate,irefd)
        end if

        call sst_wk(idate,iwk,inpfile)
        call bilinx(sstmm,sst,xlon,xlat,loni,lati,ilon,jlat,jx,iy)

        if ( ssttyp == 'OI2WK') then
          if ( idate < 19891231 ) then
            inpfile=trim(inpglob)//'/SST/icec.wkmean.1981-1989.nc'
          else
            inpfile=trim(inpglob)//'/SST/icec.wkmean.1990-present.nc'
          end if
          call ice_wk(idate,iwk,inpfile)
          call bilinx(icemm,ice,xlon,xlat,loni,lati,ilon,jlat,jx,iy)
        end if

        do i = 1 , iy
          do j = 1 , jx
            if ( sstmm(j,i) > -100. ) then
              sstmm(j,i) = sstmm(j,i) + 273.15
            else
              sstmm(j,i) = -9999.
            end if
          end do
        end do

        call writerec(idate)

        write (stdout,*) 'WRITTEN OUT SST DATA : ' , tochar(idate)
        idate = nextwk(idate)
      end do
    end if
  end subroutine sst_1deg

  subroutine sst_mn(idate,idate0,pathaddname)
    implicit none
    type(rcm_time_and_date) :: idate , idate0
    character(len=256) , intent(in) :: pathaddname
    integer(ik4) :: i , it , j
    character(len=5) :: varname
    integer(ik4) , dimension(3) , save :: icount , istart
    integer(ik4) , save :: inet , ivar
    real(rkx) , save :: xadd , xscale
    !
    ! This is the latitude, longitude dimension of the grid to be read.
    ! This corresponds to the lat and lon dimension variables in the
    ! netCDF file.
    !
    ! The data are packed into short integers (INTEGER*2).  The array
    ! work will be used to hold the packed integers. The array 'sst'
    ! will contain the unpacked data.
    !
    ! DATA ARRAY AND WORK ARRAY
    !
    data varname/'sst'/

    if ( idate == idate0 ) then
      istatus = nf90_open(pathaddname,nf90_nowrite,inet)
      call checkncerr(istatus,__FILE__,__LINE__, &
              'Error open '//trim(pathaddname))
      istatus = nf90_inq_varid(inet,varname,ivar)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find '//trim(varname))
      istatus = nf90_get_att(inet,ivar,'scale_factor',xscale)
      call checkncerr(istatus,__FILE__,__LINE__, &
              'Error find attribute scale_factor')
      istatus = nf90_get_att(inet,ivar,'add_offset',xadd)
      call checkncerr(istatus,__FILE__,__LINE__, &
              'Error find attribute add_offset')
      istatus = nf90_inq_dimid(inet,'lat',latid)
      if ( istatus /= nf90_noerr ) then
        istatus = nf90_inq_dimid(inet,'latitude',latid)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error find dim lat')
      end if
      istatus = nf90_inq_dimid(inet,'lon',lonid)
      if ( istatus /= nf90_noerr ) then
        istatus = nf90_inq_dimid(inet,'longitude',lonid)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error find dim lon')
      end if
      istatus = nf90_inquire_dimension(inet,latid,len=jlat)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error inquire dim lat')
      istatus = nf90_inquire_dimension(inet,lonid,len=ilon)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error inquire dim lon')
      istatus = nf90_inq_varid(inet,'lat',latid)
      if ( istatus /= nf90_noerr ) then
        istatus = nf90_inq_varid(inet,'latitude',latid)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error find var lat')
      end if
      istatus = nf90_inq_varid(inet,'lon',lonid)
      if ( istatus /= nf90_noerr ) then
        istatus = nf90_inq_varid(inet,'longitude',lonid)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error find var lon')
      end if
      call getmem1d(lati,1,jlat,'mod_sst_1deg:lati')
      call getmem1d(loni,1,ilon,'mod_sst_1deg:loni')
      call getmem2d(sst,1,ilon,1,jlat,'mod_sst_1deg:sst')
      call getmem2d(work,1,ilon,1,jlat,'mod_sst_1deg:work')
      istatus = nf90_get_var(inet,latid,lati)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var lat')
      istatus = nf90_get_var(inet,lonid,loni)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var lon')
      istart(1) = 1
      istart(2) = 1
      icount(1) = ilon
      icount(2) = jlat
    end if

    it = (year-1981)*12 + month - 11

    istart(3) = it
    icount(3) = 1

    istatus = nf90_get_var(inet,ivar,work,istart,icount)
    call checkncerr(istatus,__FILE__,__LINE__, &
            'Error read variable '//trim(varname))

    do j = 1 , jlat
      do i = 1 , ilon
        if ( work(i,j) == 32767 ) then
           sst(i,j) = -9999.
        else
           sst(i,j) = real(dble(work(i,j))*xscale + xadd)
        end if
      end do
    end do
  end subroutine sst_mn

  subroutine ice_mn(idate,idate0,pathaddname)
    implicit none
    type(rcm_time_and_date) :: idate , idate0
    character(len=256) , intent(in) :: pathaddname
    integer(ik4) :: i , it , j
    character(len=5) :: varname
    integer(ik4) , dimension(3) , save :: icount , istart
    integer(ik4) , save :: inet , ivar
    real(rkx) , save :: xadd , xscale
    !
    ! This is the latitude, longitude dimension of the grid to be read.
    ! This corresponds to the lat and lon dimension variables in the
    ! netCDF file.
    !
    ! The data are packed into short integers (INTEGER*2).  The array
    ! work will be used to hold the packed integers. The array 'sst'
    ! will contain the unpacked data.
    !
    ! DATA ARRAY AND WORK ARRAY
    !
    data varname/'icec'/
    if ( idate == idate0 ) then
      istatus = nf90_open(pathaddname,nf90_nowrite,inet)
      call checkncerr(istatus,__FILE__,__LINE__, &
              'Error open file '//trim(pathaddname))
      istatus = nf90_inq_varid(inet,varname,ivar)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find var '//varname)
      istatus = nf90_get_att(inet,ivar,'scale_factor',xscale)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find att scale_factor')
      istatus = nf90_get_att(inet,ivar,'add_offset',xadd)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find att add_offset')
      call getmem2d(ice,1,ilon,1,jlat,'mod_sst_1deg:ice')
      istart(1) = 1
      istart(2) = 1
      icount(1) = ilon
      icount(2) = jlat
    end if

    it = (year-1981)*12 + month - 11

    istart(3) = it
    icount(3) = 1
    istatus = nf90_get_var(inet,ivar,work,istart,icount)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var '//varname)

    do j = 1 , jlat
      do i = 1 , ilon
        if ( work(i,j) == 32767 ) then
           ice(i,j) = -9999.
        else
           ice(i,j) = real(dble(work(i,j))*xscale + xadd)
        end if
      end do
    end do
  end subroutine ice_mn

  subroutine sst_wk(idate,kkk,pathaddname)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    integer(ik4) , intent(in) :: kkk
    character(len=256) , intent(in) :: pathaddname
    integer(ik4) :: i , j
    character(len=3) :: varname
    integer(ik4) , dimension(3) , save :: icount , istart
    integer(ik4) , save :: inet , ivar
    real(rkx) , save :: xadd , xscale
    character(len=256) , save :: usename
    !
    ! This is the latitude, longitude dimension of the grid to be read.
    ! This corresponds to the lat and lon dimension variables in the
    ! netCDF file.
    !
    ! The data are packed into short integers (INTEGER*2).  The array
    ! work will be used to hold the packed integers. The array 'sst'
    ! will contain the unpacked data.
    !
    ! DATA ARRAY AND WORK ARRAY
    !
    data varname/'sst'/
    data usename/'none'/
    data inet/-1/
    if ( pathaddname /= usename ) then
      if (inet >= 0) then
        istatus = nf90_close(inet)
      end if
      istatus = nf90_open(pathaddname,nf90_nowrite,inet)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error open file '//trim(pathaddname))
      istatus = nf90_inq_varid(inet,varname,ivar)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find var '//varname)
      istatus = nf90_get_att(inet,ivar,'scale_factor',xscale)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find att scale_factor')
      istatus = nf90_get_att(inet,ivar,'add_offset',xadd)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find att add_offset')
      istatus = nf90_inq_dimid(inet,'lat',latid)
      if ( istatus /= nf90_noerr ) then
        istatus = nf90_inq_dimid(inet,'latitude',latid)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error find dim lat')
      end if
      istatus = nf90_inq_dimid(inet,'lon',lonid)
      if ( istatus /= nf90_noerr ) then
        istatus = nf90_inq_dimid(inet,'longitude',lonid)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error find dim lon')
      end if
      istatus = nf90_inquire_dimension(inet,latid,len=jlat)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error inquire dim lat')
      istatus = nf90_inquire_dimension(inet,lonid,len=ilon)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error inquire dim lon')
      istatus = nf90_inq_varid(inet,'lat',latid)
      if ( istatus /= nf90_noerr ) then
        istatus = nf90_inq_varid(inet,'latitude',latid)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error find var lat')
      end if
      istatus = nf90_inq_varid(inet,'lon',lonid)
      if ( istatus /= nf90_noerr ) then
        istatus = nf90_inq_varid(inet,'longitude',lonid)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error find var lon')
      end if
      call getmem1d(lati,1,jlat,'mod_sst_1deg:lati')
      call getmem1d(loni,1,ilon,'mod_sst_1deg:loni')
      call getmem2d(sst,1,ilon,1,jlat,'mod_sst_1deg:sst')
      call getmem2d(work,1,ilon,1,jlat,'mod_sst_1deg:work')
      call getmem2d(work1,1,ilon,1,jlat,'mod_sst_1deg:work1')
      istatus = nf90_get_var(inet,latid,lati)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var lat')
      istatus = nf90_get_var(inet,lonid,loni)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var lon')
      istart(1) = 1
      istart(2) = 1
      icount(1) = ilon
      icount(2) = jlat
      usename = pathaddname
    end if

    istart(3) = kkk
    icount(3) = 1
    istatus = nf90_get_var(inet,ivar,work,istart,icount)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var '//varname)
    if (idate < 1989123100) then
      istart(3) = kkk-1
      icount(3) = 1
      istatus = nf90_get_var(inet,ivar,work1,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//varname)
    end if

    do j = 1 , jlat
      do i = 1 , ilon
        if ( work(i,j) == 32767 ) then
           sst(i,j) = -9999.
        else
           sst(i,j) = real(dble(work(i,j))*xscale + xadd)
        end if
      end do
    end do

    if (idate < 1989123100) then
      do j = 1 , jlat
        do i = 1 , ilon
          if ( work1(i,j) == 32767 ) then
            sst(i,j) = -9999.
          else
             sst(i,j) = (sst(i,j)+real(dble(work1(i,j))*xscale+xadd))*0.5
          end if
        end do
      end do
    end if
  end subroutine sst_wk

  subroutine ice_wk(idate,kkk,pathaddname)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    integer(ik4) , intent(in) :: kkk
    character(len=256) , intent(in) :: pathaddname
    integer(ik4) :: i , j
    character(len=4) :: varname
    integer(2) , dimension(ilon,jlat) :: work , work1
    integer(ik4) , dimension(3) , save :: icount , istart
    integer(ik4) , save :: inet , ivar
    real(rkx) , save :: xadd , xscale
    character(len=256) , save :: usename
    !
    ! This is the latitude, longitude dimension of the grid to be read.
    ! This corresponds to the lat and lon dimension variables in the
    ! netCDF file.
    !
    ! The data are packed into short integers (INTEGER*2).  The array
    ! work will be used to hold the packed integers. The array 'sst'
    ! will contain the unpacked data.
    !
    ! DATA ARRAY AND WORK ARRAY
    !
    data varname/'icec'/
    data usename/'none'/
    data inet/-1/

    if ( pathaddname /= usename ) then
      if (inet >= 0) then
        istatus = nf90_close(inet)
      end if
      istatus = nf90_open(pathaddname,nf90_nowrite,inet)
      call checkncerr(istatus,__FILE__,__LINE__, &
              'Error open file '//trim(pathaddname))
      istatus = nf90_inq_varid(inet,varname,ivar)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find var '//varname)
      istatus = nf90_get_att(inet,ivar,'scale_factor',xscale)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find att scale_factor')
      istatus = nf90_get_att(inet,ivar,'add_offset',xadd)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find att add_offset')
      call getmem2d(ice,1,ilon,1,jlat,'mod_sst_1deg:ice')
      istart(1) = 1
      istart(2) = 1
      icount(1) = ilon
      icount(2) = jlat
      usename = pathaddname
    end if

    istart(3) = kkk
    icount(3) = 1
    istatus = nf90_get_var(inet,ivar,work,istart,icount)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read var '//varname)
    if (idate < 1989123100) then
      istart(3) = kkk-1
      icount(3) = 1
      istatus = nf90_get_var(inet,ivar,work1,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//varname)
    end if

    do j = 1 , jlat
      do i = 1 , ilon
        if ( work(i,j) == 32767 ) then
           ice(i,j) = -9999.
        else
           ice(i,j) = real(dble(work(i,j))*xscale + xadd)
        end if
      end do
    end do

    if (idate < 1989123100) then
      do j = 1 , jlat
        do i = 1 , ilon
          if ( work1(i,j) == 32767 ) then
             ice(i,j) = -9999.
          else
             ice(i,j) = (ice(i,j)+real(dble(work1(i,j))*xscale+xadd))*0.5
          end if
        end do
      end do
    end if
  end subroutine ice_wk

end module mod_sst_1deg
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

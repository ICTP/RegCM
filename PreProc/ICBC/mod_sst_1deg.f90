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

  use m_realkinds
  use m_die
  use m_stdio
  use m_zeit

  contains

  subroutine sst_1deg

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Comments on dataset sources and location:                          c
!                                                                    c
! GISST2.3b    UKMO SST (Rayner et al 1996), 1 degree                c
!              from UKMO DATA archive (http://www.badc.rl.ac.uk/)    c
!              and reformed as direct-accessed binary GrADS format   c
!              in file GISST_187101_200209                           c
!          ML = 1 is-179.5; ML = 2 is-178.5; => ML = 360 is 179.5E   c
!          NL = 1 is -89.5; NL = 2 is -88.5; => NL = 180 is  89.5    c
!              see the GrADS control file for details.               c
!                                                                    c
! OISST        from CAC Optimal Interpolation dataset.               c
!              in the original netCDF format.                        c
!          ML = 1 is   0.5; ML = 2 is   1.5; => ML = 360 is 359.5E   c
!          NL = 1 is -89.5; NL = 2 is -88.5; => NL = 180 is  89.5    c
!                                                                    c
! OI2ST        both SST and SeaIce in the original netCDF format.    c
!                                                                    c
! OI_WK        weekly OISST in the original netCDF format.           c
!                                                                    c
! OI2WK        weekly OISST and SeaIce in the original netCDF format.c
!                                                                    c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  use mod_sst_grid
  use mod_date
  use mod_interp , only : bilinx

  implicit none
!
  integer , parameter :: ilon = 360 , jlat = 180
!
  real(sp) , dimension(ilon,jlat) :: sst , ice
  integer :: i , j , k , iwk , iv , ludom , lumax , nrec
  integer :: idate , idateo , idatef , idatem , nsteps
  integer :: nyear , nmo , nday , nho
  real(sp) , dimension(jlat) :: lati
  real(sp) , dimension(ilon) :: loni
  integer , dimension(25) :: lund
  character(256) :: inpfile
  logical :: there
!
  call zeit_ci('sst_1deg')
  if ( ssttyp == 'GISST' ) then
    if ( globidate1 < 1947121512 .or. globidate2 > 2002091512 ) then
      write (stderr,*) 'GISST data required are not available'
      write (stderr,*) 'IDATE1, IDATE2 = ' , globidate1 , globidate2
      call die('sst_1deg')
    end if
    inquire (file=trim(inpglob)//'/SST/GISST_194712_200209',exist=there)
    if ( .not. there ) then
      call die('sst_1deg','GISST_194712_200209 is not available'//  &
               ' under '//trim(inpglob)//'/SST/',1)
    end if
    open (11,file=trim(inpglob)//'/SST/GISST_194712_200209',        &
          form='unformatted',recl=ilon*jlat*ibyte,access='direct',  &
          status = 'old')
  else if ( ssttyp == 'OISST' .or. ssttyp == 'OI_NC' .or. &
            ssttyp == 'OI2ST' ) then
    if ( globidate1 < 1981121512 .or. globidate2 < 1981121512 ) then
      write (stderr,*) 'OISST data required are not available'
      write (stderr,*) 'IDATE1, IDATE2 = ' , globidate1 , globidate2
      call die('sst_1deg')
    end if
    inquire (file=trim(inpglob)//'/SST/sst.mnmean.nc',exist=there)
    if ( .not.there ) then
      call die('sst_1deg','sst.mnmean.nc is not available'// &
               ' under '//trim(inpglob)//'/SST/',1)
    end if
    if ( ssttyp == 'OI2ST' ) then
      inquire (file=trim(inpglob)//'/SST/icec.mnmean.nc',exist=there)
      if ( .not. there ) then
        call die('sst_1deg','icec.mnmean.nc is not available'// &
                 ' under '//trim(inpglob)//'/SST/',1)
      end if
    end if
  else if ( ssttyp == 'OI_WK' .or. ssttyp == 'OI2WK' ) then
    if ( globidate1 < 1981110100 .or. globidate2 < 1981110106 ) then
      write (stderr,*) 'OI_WK (or OI2WK) data are not available'
      write (stderr,*) 'IDATE1, IDATE2 = ' , globidate1 , globidate2
      call die('sst_1deg')
    end if
    inquire (file=trim(inpglob)//'/SST/sst.wkmean.1981-1989.nc',exist=there)
    if ( .not.there ) then
      call die('sst_1deg','sst.wkmean.1981-1989.nc is not '// &
               'available under '//trim(inpglob)//'/SST/',1)
    end if
    inquire (file=trim(inpglob)//'/SST/sst.wkmean.1990-present.nc', &
             exist=there)
    if ( .not.there ) then
      call die('sst_1deg','sst.wkmean.1990-present.nc is not '// &
               'available under '//trim(inpglob)//'/SST/',1)
    end if
    if ( ssttyp == 'OI2WK' ) then
      inquire (file=trim(inpglob)//'/SST/icec.wkmean.1981-1989.nc', &
              exist=there)
      if ( .not.there ) then
        call die('sst_1deg','icec.wkmean.1981-1989.nc is not '// &
                 'available under '//trim(inpglob)//'/SST/',1)
      end if
      inquire (file=trim(inpglob)//'/SST/icec.wkmean.1990-present.nc', &
               exist=there)
      if ( .not.there ) then
        call die('sst_1deg','icec.wkmean.1990-present.nc is not '// &
                 'available under '//trim(inpglob)//'/SST/',1)
      end if
    end if
  else
    write (stderr,*) 'PLEASE SET right SSTTYP in regcm.in'
    write (stderr,*) 'Supported are GISST OISST OI_NC OI2ST OI_WK OI2WK'
    call die('sst_1deg')
  end if

  ! Montly dataset
  if ( ssttyp /= 'OI_WK' .and. ssttyp /= 'OI2WK' ) then
    idateo = imonfirst(globidate1)
    if (lfhomonth(globidate1)) then
      idateo = iprevmon(globidate1)
    end if
    idatef = imonfirst(globidate2)
    if (idatef < globidate2) then
      idatef = inextmon(idatef)
    end if
    nsteps = imondiff(idatef,idateo) + 1
    idatem = imonmiddle(idateo)
    call open_sstfile(idatem)
  ! Weekly dataset
  else
    idateo = ifodweek(globidate1)
    if (idateo > globidate1) then
      idateo = iprevwk(idateo)
    end if
    idatef = ifodweek(globidate2)
    if (idatef < globidate2) then
      idatef = inextwk(idatef)
    end if
    nsteps = iwkdiff(idatef,idateo) + 1
    call open_sstfile(idateo)
  end if


  ! SET UP LONGITUDES AND LATITUDES FOR SST DATA
  do i = 1 , ilon
    loni(i) = .5 + float(i-1)
  end do
  do j = 1 , jlat
    lati(j) = -89.5 + 1.*float(j-1)
  end do
 
!       ****** OISST SST DATA, 1 Deg data, AVAILABLE FROM 12/1981 TO
!       PRESENT ****** GISST SST DATA, 1 Deg data, AVAILABLE FROM
!       12/1947 TO 9/2002

  idate = idateo

  if ( ssttyp /= 'OI_WK' .and. ssttyp /= 'OI2WK' ) then

    do k = 1 , nsteps

      call split_idate(idate, nyear, nmo, nday, nho)

      if ( ssttyp == 'GISST' ) then
        nrec = (nyear-1947)*12 + nmo - 11
        read (11,rec=nrec) sst
      else if ( ssttyp == 'OISST' .or. ssttyp == 'OI_NC' .or. &
                ssttyp == 'OI2ST') then
        inpfile=trim(inpglob)//'/SST/sst.mnmean.nc'
        call sst_mn(idate,idateo,ilon,jlat,sst,inpfile)
        if ( ssttyp == 'OI2ST' ) then
          inpfile=trim(inpglob)//'/SST/icec.mnmean.nc'
          call ice_mn(idate,idateo,ilon,jlat,ice,inpfile)
        end if
      end if
 
      call bilinx(sst,sstmm,xlon,xlat,loni,lati,ilon,jlat,iy,jx,1)
      if ( ssttyp == 'OI2ST' ) then
        call bilinx(ice,icemm,xlon,xlat,loni,lati,ilon,jlat,iy,jx,1)
      end if

      write (stdout,*) 'XLON,XLAT,SST = ', xlon(1,1), xlat(1,1), sstmm(1,1)
 
      do j = 1 , jx
        do i = 1 , iy
          if ( sstmm(i,j) < -5000. .and. &
               (lu(i,j) > 13.5 .and. lu(i,j) < 15.5) ) then
            do iv = 1 , 20
              lund(iv) = 0
            end do
            lund(nint(lu(i-1,j-1))) = lund(nint(lu(i-1,j-1))) + 2
            lund(nint(lu(i-1,j))) = lund(nint(lu(i-1,j))) + 3
            lund(nint(lu(i-1,j+1))) = lund(nint(lu(i-1,j+1))) + 2
            lund(nint(lu(i,j-1))) = lund(nint(lu(i,j-1))) + 3
            lund(nint(lu(i,j+1))) = lund(nint(lu(i,j+1))) + 3
            lund(nint(lu(i+1,j-1))) = lund(nint(lu(i+1,j-1))) + 2
            lund(nint(lu(i+1,j))) = lund(nint(lu(i+1,j))) + 3
            lund(nint(lu(i+1,j+1))) = lund(nint(lu(i+1,j+1))) + 2
            ludom = 18
            lumax = 0
            do iv = 1 , 20
              if ( iv <= 13 .or. iv >= 16 ) then
                if ( lund(iv) > lumax ) then
                  ludom = iv
                  lumax = lund(iv)
                end if
              end if
            end do
            lu(i,j) = float(ludom)
            write (stdout,*) ludom , sstmm(i,j)
          end if
          if ( sstmm(i,j) > -100. ) then
            sstmm(i,j) = sstmm(i,j) + 273.15
          else
            sstmm(i,j) = -9999.
          end if
        end do
      end do
 
!         ******           WRITE OUT SST DATA ON MM4 GRID
      if ( ssttyp /= 'OI2ST' ) then
        call writerec(idatem,.false.)
      else
        call writerec(idatem,.true.)
      end if

      write (stdout,*) 'WRITTEN OUT SST DATA : ' , idate

      idate = inextmon(idate)
      idatem = imonmiddle(idate)

    end do

  ! Weekly data
  else

    do k = 1 , nsteps

      call split_idate(idate, nyear, nmo, nday, nho)

      if ( idate < 1989123100 ) then
        inpfile=trim(inpglob)//'/SST/sst.wkmean.1981-1989.nc'
        iwk = iwkdiff(idate,1981110100) + 1
      else
        inpfile=trim(inpglob)//'/SST/sst.wkmean.1990-present.nc'
        iwk = iwkdiff(idate,1989123100) + 1
      end if

      call sst_wk(idate,iwk,ilon,jlat,sst,inpfile)
      call bilinx(sst,sstmm,xlon,xlat,loni,lati,ilon,jlat,iy,jx,1)
 
      if ( ssttyp == 'OI2WK') then
        if ( idate < 19891231 ) then
          inpfile=trim(inpglob)//'/SST/icec.wkmean.1981-1989.nc'
        else
          inpfile=trim(inpglob)//'/SST/icec.wkmean.1990-present.nc'
        end if
        call ice_wk(idate,iwk,ilon,jlat,ice,inpfile)
        call bilinx(ice,icemm,xlon,xlat,loni,lati,ilon,jlat,iy,jx,1)
      end if 

      write (stdout,*) 'XLON,XLAT,SST = ', xlon(1,1), xlat(1,1), sstmm(1,1)

      do j = 1 , jx
        do i = 1 , iy
          if ( sstmm(i,j) > -100. ) then
            sstmm(i,j) = sstmm(i,j) + 273.15
          else
            sstmm(i,j) = -9999.
          end if
        end do
      end do
 
!         ******           WRITE OUT SST DATA ON MM4 GRID
      if (ssttyp == 'OI_WK') then
        call writerec(idate,.false.)
      else
        call writerec(idate,.true.)
      endif

      write (stdout,*) 'WRITTEN OUT SST DATA : ' , idate

      idate = inextwk(idate)

    end do
  end if
  call zeit_co('sst_1deg')

  end subroutine sst_1deg
!
!-----------------------------------------------------------------------
!
  subroutine sst_mn(idate,idate0,ilon,jlat,sst,pathaddname)
  use netcdf
  use mod_date , only : split_idate
  implicit none
!
  integer :: idate , idate0 , ilon , jlat
  character(256) :: pathaddname
  intent (in) idate , idate0 , ilon , jlat , pathaddname
  real(sp) , dimension(ilon,jlat) :: sst
  intent (out) :: sst
!
  integer :: i , it , j , month , n , nday , nhour , nyear
  logical :: there
  character(5) :: varname
  integer(2) , dimension(ilon,jlat) :: work
  integer :: istatus
!
  integer , dimension(10) , save :: icount , istart
  integer , save :: inet , ivar
  real(dp) , save :: xadd , xscale
!
!     This is the latitude, longitude dimension of the grid to be read.
!     This corresponds to the lat and lon dimension variables in the
!     netCDF file.
!
!     The data are packed into short integers (INTEGER*2).  The array
!     work will be used to hold the packed integers. The array 'sst'
!     will contain the unpacked data.
!
!     DATA ARRAY AND WORK ARRAY
!
  data varname/'sst'/
!
  call zeit_ci('sst_mn')
  if ( idate == idate0 ) then
    inquire (file=pathaddname,exist=there)
    if ( .not.there ) then
      call die('sst_mn',trim(pathaddname)//' is not available',1)
    end if
    istatus = nf90_open(pathaddname,nf90_nowrite,inet)
    if ( istatus /= nf90_noerr ) then
      call die('sst_mn',trim(pathaddname)//' open error',1, &
              nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_inq_varid(inet,varname,ivar)
    if ( istatus /= nf90_noerr ) then
      call die('sst_mn',trim(pathaddname)//':'//trim(varname)// &
               ' error',1,nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_get_att(inet,ivar,'scale_factor',xscale)
    if ( istatus /= nf90_noerr ) then
      call die('sst_mn',trim(pathaddname)//':'//trim(varname)// &
               ':scale_factor',1,nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_get_att(inet,ivar,'add_offset',xadd)
    if ( istatus /= nf90_noerr ) then
      call die('sst_mn',trim(pathaddname)//':'//trim(varname)// &
               ':add_offset',1,nf90_strerror(istatus),istatus)
    end if
    istart(1) = 1
    istart(2) = 1
    icount(1) = ilon
    icount(2) = jlat
    do n = 4 , 10
      istart(n) = 0
      icount(n) = 0
    end do
  end if
 
  call split_idate(idate, nyear, month, nday, nhour)
  it = (nyear-1981)*12 + month - 11
 
  istart(3) = it
  icount(3) = 1

  istatus = nf90_get_var(inet,ivar,work,istart,icount)
  if ( istatus /= nf90_noerr ) then
    write (stderr,*) istart
    write (stderr,*) icount
    call die('sst_mn',trim(pathaddname)//':'//trim(varname)// &
             ':read',1,nf90_strerror(istatus),istatus)
  end if
!
  do j = 1 , jlat
    do i = 1 , ilon
      if ( work(i,j) == 32767 ) then
         sst(i,jlat+1-j) = -9999.
      else
         sst(i,jlat+1-j) = real(dble(work(i,j))*xscale + xadd)
      end if
    end do
  end do
  call zeit_co('sst_mn')
!
  end subroutine sst_mn
!
!-----------------------------------------------------------------------
!
  subroutine ice_mn(idate,idate0,ilon,jlat,ice,pathaddname)
  use netcdf
  use mod_date , only : split_idate
  implicit none
!
  integer :: idate , idate0 , ilon , jlat
  character(256) :: pathaddname
  intent (in) idate , idate0 , ilon , jlat , pathaddname
  real(sp) , dimension(ilon,jlat) :: ice
  intent (out) :: ice
!
  integer :: i , it , j , month , n , nday , nhour , nyear
  logical :: there
  character(5) :: varname
  integer(2) , dimension(ilon,jlat) :: work
  integer :: istatus
!
  integer , dimension(10) , save :: icount , istart
  integer , save :: inet , ivar
  real(dp) , save :: xadd , xscale
!
!     This is the latitude, longitude dimension of the grid to be read.
!     This corresponds to the lat and lon dimension variables in the
!     netCDF file.
!
!     The data are packed into short integers (INTEGER*2).  The array
!     work will be used to hold the packed integers. The array 'sst'
!     will contain the unpacked data.
!
!     DATA ARRAY AND WORK ARRAY
!
  data varname/'icec'/
!
  call zeit_ci('ice_mn')
  if ( idate == idate0 ) then
    inquire (file=pathaddname,exist=there)
    if ( .not.there ) then
      call die('ice_mn',trim(pathaddname)//' is not available',1)
    end if
    istatus = nf90_open(pathaddname,nf90_nowrite,inet)
    if ( istatus /= nf90_noerr ) then
      call die('ice_mn',trim(pathaddname)//' open error',1, &
              nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_inq_varid(inet,varname,ivar)
    if ( istatus /= nf90_noerr ) then
      call die('ice_mn',trim(pathaddname)//':'//trim(varname)// &
               ' error',1,nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_get_att(inet,ivar,'scale_factor',xscale)
    if ( istatus /= nf90_noerr ) then
      call die('ice_mn',trim(pathaddname)//':'//trim(varname)// &
               ':scale_factor',1,nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_get_att(inet,ivar,'add_offset',xadd)
    if ( istatus /= nf90_noerr ) then
      call die('ice_mn',trim(pathaddname)//':'//trim(varname)// &
               ':add_offset',1,nf90_strerror(istatus),istatus)
    end if
    istart(1) = 1
    istart(2) = 1
    icount(1) = ilon
    icount(2) = jlat
    do n = 4 , 10
      istart(n) = 0
      icount(n) = 0
    end do
  end if
 
  call split_idate(idate, nyear, month, nday, nhour)
  it = (nyear-1981)*12 + month - 11
 
  istart(3) = it
  icount(3) = 1
  istatus = nf90_get_var(inet,ivar,work,istart,icount)
  if ( istatus /= nf90_noerr ) then
    write (stderr,*) istart
    write (stderr,*) icount
    call die('ice_mn',trim(pathaddname)//':'//trim(varname)// &
             ':read',1,nf90_strerror(istatus),istatus)
  end if
!
  do j = 1 , jlat
    do i = 1 , ilon
      if ( work(i,j) == 32767 ) then
         ice(i,jlat+1-j) = -9999.
      else
         ice(i,jlat+1-j) = real(dble(work(i,j))*xscale + xadd)
      end if
    end do
  end do
  call zeit_co('ice_mn')
!
  end subroutine ice_mn
!
!-----------------------------------------------------------------------
!
  subroutine sst_wk(idate,kkk,ilon,jlat,sst,pathaddname)
  use netcdf
  implicit none
!
  integer :: idate , kkk , ilon , jlat
  character(256) :: pathaddname
  intent (in) idate , kkk , ilon , jlat , pathaddname
  real(sp) , dimension(ilon,jlat) :: sst
  intent (out) :: sst
!
  integer :: i , j , n
  logical :: there
  character(3) :: varname
  integer :: istatus
  integer(2) , dimension(ilon,jlat) :: work , work1
!
  integer , dimension(10) , save :: icount , istart
  integer , save :: inet , ivar
  real(dp) , save :: xadd , xscale
  character(256) , save :: usename
!
!     This is the latitude, longitude dimension of the grid to be read.
!     This corresponds to the lat and lon dimension variables in the
!     netCDF file.
!
!     The data are packed into short integers (INTEGER*2).  The array
!     work will be used to hold the packed integers. The array 'sst'
!     will contain the unpacked data.
!
!     DATA ARRAY AND WORK ARRAY
!
  data varname/'sst'/
  data usename/'none'/
  data inet/-1/
!
  call zeit_ci('sst_wk')
  if ( pathaddname /= usename ) then
    if (inet >= 0) then
      istatus = nf90_close(inet)
    end if
    inquire (file=pathaddname,exist=there)
    if ( .not.there ) then
      call die('sst_wk',trim(pathaddname)//' is not available',1)
    end if
    istatus = nf90_open(pathaddname,nf90_nowrite,inet)
    if ( istatus /= nf90_noerr ) then
      call die('sst_wk',trim(pathaddname)//' open error',1, &
              nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_inq_varid(inet,varname,ivar)
    if ( istatus /= nf90_noerr ) then
      call die('sst_wk',trim(pathaddname)//':'//trim(varname)// &
               ' error',1,nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_get_att(inet,ivar,'scale_factor',xscale)
    if ( istatus /= nf90_noerr ) then
      call die('sst_wk',trim(pathaddname)//':'//trim(varname)// &
               ':scale_factor',1,nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_get_att(inet,ivar,'add_offset',xadd)
    if ( istatus /= nf90_noerr ) then
      call die('sst_wk',trim(pathaddname)//':'//trim(varname)// &
               ':add_offset',1,nf90_strerror(istatus),istatus)
    end if
    istart(1) = 1
    istart(2) = 1
    icount(1) = ilon
    icount(2) = jlat
    do n = 4 , 10
      istart(n) = 0
      icount(n) = 0
    end do
    usename = pathaddname
  end if

  istart(3) = kkk
  icount(3) = 1
  istatus = nf90_get_var(inet,ivar,work,istart,icount)
  if ( istatus /= nf90_noerr ) then
    write (stderr,*) istart
    write (stderr,*) icount
    call die('sst_wk',trim(pathaddname)//':'//trim(varname)// &
             ':read',1,nf90_strerror(istatus),istatus)
  end if
  if (idate < 1989123100) then
    istart(3) = kkk-1
    icount(3) = 1
    istatus = nf90_get_var(inet,ivar,work1,istart,icount)
    if ( istatus /= nf90_noerr ) then
      write (stderr,*) istart
      write (stderr,*) icount
      call die('sst_wk',trim(pathaddname)//':'//trim(varname)// &
               ':read',1,nf90_strerror(istatus),istatus)
    end if
  end if

  do j = 1 , jlat
    do i = 1 , ilon
      if ( work(i,j) == 32767 ) then
         sst(i,jlat+1-j) = -9999.
      else
         sst(i,jlat+1-j) = real(dble(work(i,j))*xscale + xadd)
      end if
    end do
  end do

  if (idate < 1989123100) then
    do j = 1 , jlat
      do i = 1 , ilon
        if ( work1(i,j) == 32767 ) then
           sst(i,jlat+1-j) = -9999.
        else
           sst(i,jlat+1-j) = &
              (sst(i,jlat+1-j)+real(dble(work1(i,j))*xscale+xadd))*0.5
        end if
      end do
    end do
  end if
  call zeit_co('sst_wk')

  end subroutine sst_wk
!
!-----------------------------------------------------------------------
!
  subroutine ice_wk(idate,kkk,ilon,jlat,ice,pathaddname)
  use netcdf
  implicit none
!
  integer :: idate , kkk , ilon , jlat
  character(256) :: pathaddname
  intent (in) idate , kkk , ilon , jlat , pathaddname
  real(sp) , dimension(ilon,jlat) :: ice
  intent (out) :: ice
!
  integer :: i , j , n
  logical :: there
  character(4) :: varname
  integer(2) , dimension(ilon,jlat) :: work , work1
  integer :: istatus
!
  integer , dimension(10) , save :: icount , istart
  integer , save :: inet , ivar
  real(dp) , save :: xadd , xscale
  character(256) , save :: usename
!
!     This is the latitude, longitude dimension of the grid to be read.
!     This corresponds to the lat and lon dimension variables in the
!     netCDF file.
!
!     The data are packed into short integers (INTEGER*2).  The array
!     work will be used to hold the packed integers. The array 'sst'
!     will contain the unpacked data.
!
!     DATA ARRAY AND WORK ARRAY
!
  data varname/'icec'/
  data usename/'none'/
  data inet/-1/
!
  call zeit_ci('ice_wk')
  if ( pathaddname /= usename ) then
    if (inet >= 0) then
      istatus = nf90_close(inet)
    end if
    inquire (file=pathaddname,exist=there)
    if ( .not.there ) then
      call die('ice_wk',trim(pathaddname)//' is not available',1)
    end if
    istatus = nf90_open(pathaddname,nf90_nowrite,inet)
    if ( istatus /= nf90_noerr ) then
      call die('ice_wk',trim(pathaddname)//' open error',1, &
              nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_inq_varid(inet,varname,ivar)
    if ( istatus /= nf90_noerr ) then
      call die('ice_wk',trim(pathaddname)//':'//trim(varname)// &
               ' error',1,nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_get_att(inet,ivar,'scale_factor',xscale)
    if ( istatus /= nf90_noerr ) then
      call die('ice_wk',trim(pathaddname)//':'//trim(varname)// &
               ':scale_factor',1,nf90_strerror(istatus),istatus)
    end if
    istatus = nf90_get_att(inet,ivar,'add_offset',xadd)
    if ( istatus /= nf90_noerr ) then
      call die('ice_wk',trim(pathaddname)//':'//trim(varname)// &
               ':add_offset',1,nf90_strerror(istatus),istatus)
    end if
    istart(1) = 1
    istart(2) = 1
    icount(1) = ilon
    icount(2) = jlat
    do n = 4 , 10
      istart(n) = 0
      icount(n) = 0
    end do
    usename = pathaddname
  end if

  istart(3) = kkk
  icount(3) = 1
  istatus = nf90_get_var(inet,ivar,work,istart,icount)
  if ( istatus /= nf90_noerr ) then
    write (stderr,*) istart
    write (stderr,*) icount
    call die('ice_wk',trim(pathaddname)//':'//trim(varname)// &
             ':read',1,nf90_strerror(istatus),istatus)
  end if
  if (idate < 1989123100) then
    istart(3) = kkk-1
    icount(3) = 1
    istatus = nf90_get_var(inet,ivar,work1,istart,icount)
    if ( istatus /= nf90_noerr ) then
      write (stderr,*) istart
      write (stderr,*) icount
      call die('ice_wk',trim(pathaddname)//':'//trim(varname)// &
               ':read',1,nf90_strerror(istatus),istatus)
    end if
  end if

  do j = 1 , jlat
    do i = 1 , ilon
      if ( work(i,j) == 32767 ) then
         ice(i,jlat+1-j) = -9999.
      else
         ice(i,jlat+1-j) = real(dble(work(i,j))*xscale + xadd)
      end if
    end do
  end do

  if (idate < 1989123100) then
    do j = 1 , jlat
      do i = 1 , ilon
        if ( work1(i,j) == 32767 ) then
           ice(i,jlat+1-j) = -9999.
        else
           ice(i,jlat+1-j) = &
                 (ice(i,jlat+1-j)+real(dble(work1(i,j))*xscale+xadd))*0.5
        end if
      end do
    end do
  end if
  call zeit_co('ice_wk')

  end subroutine ice_wk
!
end module mod_sst_1deg

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

module mod_ein

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_stdio
  use mod_memutil
  use mod_grid
  use mod_write
  use mod_interp
  use mod_vertint
  use mod_hgt
  use mod_humid
  use mod_mksst
  use mod_uvrot
  use mod_vectutil
  use mod_message
  use mod_nchelper
  use netcdf

  private

  integer(ik4) :: jlat , ilon , klev , timlen

  real(rkx) , pointer , dimension(:,:,:) :: b3
  real(rkx) , pointer , dimension(:,:,:) :: d3

  real(rkx) , pointer :: u3(:,:,:) , v3(:,:,:)
  real(rkx) , pointer :: h3(:,:,:) , q3(:,:,:) , t3(:,:,:)
  real(rkx) , pointer :: uvar(:,:,:) , vvar(:,:,:)
  real(rkx) , pointer :: hvar(:,:,:) , rhvar(:,:,:) , tvar(:,:,:)

  real(rkx) , pointer , dimension(:,:,:) :: b2
  real(rkx) , pointer , dimension(:,:,:) :: d2
  real(rkx) , pointer , dimension(:) :: glat
  real(rkx) , pointer , dimension(:) :: grev
  real(rkx) , pointer , dimension(:) :: glon
  !real(rkx) , pointer , dimension(:,:) :: glat2
  !real(rkx) , pointer , dimension(:,:) :: glon2
  integer(ik4) , pointer , dimension(:) :: plevs
  real(rkx) , pointer , dimension(:) :: sigma1 , sigmar
  real(rkx) :: pss
  integer(2) , pointer , dimension(:,:,:) :: work

  integer(ik4) , dimension(5,4) :: inet5
  integer(ik4) , dimension(5,4) :: ivar5
  real(rkx) , dimension(5,4) :: xoff , xscl
  type(rcm_time_and_date) , pointer , dimension(:) :: itimes
  integer(ik4) , pointer , dimension(:) :: xtimes

  type(global_domain) :: gdomain

  public :: getein , headerein

  contains

  subroutine getein(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    !
    ! Read data at idate
    !
    call ein6hour(dattyp,idate,globidate1)
    write (stdout,*) 'READ IN fields at DATE:' , tochar(idate)
    !
    ! Horizontal interpolation of both the scalar and vector fields
    !
    call bilinx(b3,b2,xlon,xlat,glon,glat,ilon,jlat,jx,iy,klev*3)
    call bilinx(d3,d2,xlon,xlat,glon,glat,ilon,jlat,jx,iy,klev*2)
    !call cressmcr(b3,b2,xlon,xlat,glon2,glat2,jx,iy,ilon,jlat,klev,3)
    !call cressmdt(d3,d2,dlon,dlat,glon2,glat2,jx,iy,ilon,jlat,klev,2)
    !
    ! Rotate u-v fields after horizontal interpolation
    !
    call uvrot4(u3,v3,dlon,dlat,clon,clat,xcone,jx,iy,klev,plon,plat,iproj)
    !
    ! Invert vertical order top -> bottom for RegCM convention
    !
    call top2btm(t3,jx,iy,klev)
    call top2btm(q3,jx,iy,klev)
    call top2btm(h3,jx,iy,klev)
    call top2btm(u3,jx,iy,klev)
    call top2btm(v3,jx,iy,klev)
    !
    ! Vertical interpolation
    ! New calculation of p* on rcm topography.
    !
    call intgtb(pa,za,tlayer,topogm,t3,h3,pss,sigmar,jx,iy,klev)
    call intpsn(ps4,topogm,pa,za,tlayer,ptop,jx,iy)
    call crs2dot(pd4,ps4,jx,iy,i_band)
    !
    ! Interpolation from pressure levels
    !
    call intv3(ts4,t3,ps4,pss,sigmar,ptop,jx,iy,klev)
    call readsst(ts4,idate)
    !
    ! Interpolate U, V, T, and Q.
    !
    call intv1(u4,u3,pd4,sigmah,pss,sigmar,ptop,jx,iy,kz,klev)
    call intv1(v4,v3,pd4,sigmah,pss,sigmar,ptop,jx,iy,kz,klev)
    call intv2(t4,t3,ps4,sigmah,pss,sigmar,ptop,jx,iy,kz,klev)
    call intv1(q4,q3,ps4,sigmah,pss,sigmar,ptop,jx,iy,kz,klev)
    !
    ! Get back to mixing ratio
    !
    call rh2mxr(t4,q4,ps4,ptop,sigmah,jx,iy,kz)
  end subroutine getein
!
!-----------------------------------------------------------------------
!
  subroutine ein6hour(dattyp,idate,idate0)
    implicit none
    character(len=5) , intent(in) :: dattyp
    type(rcm_time_and_date) , intent(in) :: idate , idate0
    integer(ik4) :: i , inet , it , j , k4 , kkrec , istatus , ivar
    integer(ik4) :: timid
    character(len=64) :: inname
    character(len=256) :: pathaddname
    character(len=1) , dimension(5) :: varname
    character(len=4) , dimension(5) :: fname
    character(len=4) , dimension(4) :: hname
    character(len=64) :: cunit , ccal
    real(rkx) :: xadd , xscale
    integer(ik4) , dimension(4) :: icount , istart
    integer(ik4) :: year , month , day , hour , monthp1
    integer(ik4) , save :: lastmonth , lastyear
    type(rcm_time_and_date) :: xdate
    type(rcm_time_interval) :: tdif
    !
    ! This is the latitude, longitude dimension of the grid to be read.
    ! This corresponds to the lat and lon dimension variables in the
    ! netCDF file.
    ! The data are packed into short integers (INTEGER*2).  The array
    ! work will be used to hold the packed integers.  The array 'x'
    ! will contain the unpacked data.
    !
    data varname /'t' , 'z' , 'r' , 'u' , 'v'/
    data fname   /'air','hgt','rhum','uwnd','vwnd'/
    data hname   /'.00.','.06.','.12.','.18.'/

    k4 = 1
    call split_idate(idate,year,month,day,hour)

    if ( dattyp == 'EIXXX' ) then
      if ( idate == idate0 .or. month /= lastmonth ) then
        lastmonth = month
        do kkrec = 1 , 5
          monthp1 = month+1
          if ( monthp1 == 13 ) monthp1 = 1
          write(inname,'(a,i0.2,a,i0.2,a)') &
             varname(kkrec)//'_xxxx',month,'0100-xxxx',monthp1,'0100.nc'
          pathaddname = trim(inpglob)//pthsep//'ERAIN_MEAN'//&
                        pthsep//'XXXX'//pthsep//inname
          istatus = nf90_open(pathaddname,nf90_nowrite,inet5(kkrec,1))
          call checkncerr(istatus,__FILE__,__LINE__, &
            'Error open file '//trim(pathaddname))
          istatus = nf90_inq_varid(inet5(kkrec,1),varname(kkrec), &
                                   ivar5(kkrec,1))
          call checkncerr(istatus,__FILE__,__LINE__, &
            'Error find var '//varname(kkrec))
          istatus = nf90_get_att(inet5(kkrec,1),ivar5(kkrec,1), &
                   'scale_factor',xscl(kkrec,1))
          call checkncerr(istatus,__FILE__,__LINE__, &
            'Error find att scale_factor')
          istatus = nf90_get_att(inet5(kkrec,1),ivar5(kkrec,1),  &
                   'add_offset',xoff(kkrec,1))
          call checkncerr(istatus,__FILE__,__LINE__, &
            'Error find att add_offset')
          write (stdout,*) inet5(kkrec,1) , trim(pathaddname) ,   &
                           xscl(kkrec,1) , xoff(kkrec,1)
          if ( kkrec == 1 ) then
            istatus = nf90_inq_dimid(inet5(1,1),'time',timid)
            call checkncerr(istatus,__FILE__,__LINE__, &
                            'Error find dim time')
            istatus = nf90_inquire_dimension(inet5(1,1),timid,len=timlen)
            call checkncerr(istatus,__FILE__,__LINE__, &
                            'Error inquire time')
            istatus = nf90_inq_varid(inet5(1,1),'time',timid)
            if ( istatus /= nf90_noerr ) then
              istatus = nf90_inq_varid(inet5(1,1),'date',timid)
              call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error find var time/date')
            end if
            cunit = 'hours since 1900-01-01 00:00:00'
            ccal = 'noleap'
            call getmem1d(itimes,1,timlen,'mod_ein:itimes')
            call getmem1d(xtimes,1,timlen,'mod_ein:xtimes')
            istatus = nf90_get_var(inet5(1,1),timid,xtimes)
            call checkncerr(istatus,__FILE__,__LINE__, &
                            'Error read time')
            do it = 1 , timlen
              itimes(it) = timeval2date(real(xtimes(it),rkx),cunit,ccal)
            end do
          end if
        end do
      end if
      xdate = 1979000000 + month*10000+day*100+hour
      call setcal(xdate,'noleap')
      tdif = xdate - itimes(1)
      it = nint(tohours(tdif))/6 + 1
    else
      if ( idate == idate0 .or. year /= lastyear ) then
        lastyear = year
        do k4 = 1 , 4
          do kkrec = 1 , 5
            write(inname,'(i4,a,a,i4,a)') &
              year, pthsep, trim(fname(kkrec))//'.', year, hname(k4)//'nc'
            pathaddname = trim(inpglob)//pthsep//dattyp//pthsep//inname
            istatus = nf90_open(pathaddname,nf90_nowrite,inet5(kkrec,k4))
            call checkncerr(istatus,__FILE__,__LINE__, &
              'Error open file '//trim(pathaddname))
            istatus = nf90_inq_varid(inet5(kkrec,k4),varname(kkrec), &
                                     ivar5(kkrec,k4))
            call checkncerr(istatus,__FILE__,__LINE__, &
              'Error find var '//varname(kkrec))
            istatus = nf90_get_att(inet5(kkrec,k4),ivar5(kkrec,k4), &
                     'scale_factor',xscl(kkrec,k4))
            call checkncerr(istatus,__FILE__,__LINE__, &
              'Error find att scale_factor')
            istatus = nf90_get_att(inet5(kkrec,k4),ivar5(kkrec,k4),  &
                     'add_offset',xoff(kkrec,k4))
            call checkncerr(istatus,__FILE__,__LINE__, &
              'Error find att add_offset')
            write (stdout,*) inet5(kkrec,k4) , trim(pathaddname) ,   &
                             xscl(kkrec,k4) , xoff(kkrec,k4)
            if ( k4 == 1 .and. kkrec == 1 ) then
              istatus = nf90_inq_dimid(inet5(1,1),'time',timid)
              call checkncerr(istatus,__FILE__,__LINE__, &
                              'Error find dim time')
              istatus = nf90_inquire_dimension(inet5(1,1),timid,len=timlen)
              call checkncerr(istatus,__FILE__,__LINE__, &
                              'Error inquire time')
              istatus = nf90_inq_varid(inet5(1,1),'time',timid)
              call checkncerr(istatus,__FILE__,__LINE__, &
                              'Error find var time')
              istatus = nf90_get_att(inet5(1,1),timid,'units',cunit)
              call checkncerr(istatus,__FILE__,__LINE__, &
                              'Error read time units')
              ccal = 'gregorian'
              call getmem1d(itimes,1,timlen,'mod_ein:itimes')
              call getmem1d(xtimes,1,timlen,'mod_ein:xtimes')
              istatus = nf90_get_var(inet5(1,1),timid,xtimes)
              call checkncerr(istatus,__FILE__,__LINE__, &
                              'Error read time')
              do it = 1 , timlen
                itimes(it) = timeval2date(real(xtimes(it),rkx),cunit,ccal)
              end do
            end if
          end do
        end do
      end if
      k4 = hour/6 + 1
      tdif = idate - itimes(1)
      it = nint(tohours(tdif))/24 + 1
    end if

    istart(3) = 1
    icount(3) = klev
    istart(4) = it
    icount(4) = 1

    do kkrec = 1 , 5
      inet = inet5(kkrec,k4)
      ivar = ivar5(kkrec,k4)
      xscale = xscl(kkrec,k4)
      xadd = xoff(kkrec,k4)
      call getwork(kkrec)
      if ( kkrec == 1 ) then
        do j = 1 , jlat
          do i = 1 , ilon
            tvar(i,j,:) = real(real(work(i,j,:),rkx)*xscale+xadd,rkx)
          end do
        end do
      else if ( kkrec == 2 ) then
        do j = 1 , jlat
          do i = 1 , ilon
            hvar(i,j,:) = real(real(work(i,j,:),rkx) * &
                      xscale+xadd,rkx)/9.80616_rk4
          end do
        end do
      else if ( kkrec == 3 ) then
        do j = 1 , jlat
          do i = 1 , ilon
            rhvar(i,j,:) = &
              max(real(real(work(i,j,:),rkx)*xscale+xadd,rkx)*0.01_rkx,0.0_rkx)
          end do
        end do
      else if ( kkrec == 4 ) then
        do j = 1 , jlat
          do i = 1 , ilon
            uvar(i,j,:) = real(real(work(i,j,:),rkx)*xscale+xadd,rkx)
          end do
        end do
      else if ( kkrec == 5 ) then
        do j = 1 , jlat
          do i = 1 , ilon
            vvar(i,j,:) = real(real(work(i,j,:),rkx)*xscale+xadd,rkx)
          end do
        end do
      end if
    end do

    contains

      subroutine getwork(irec)
        implicit none
        integer(ik4) , intent(in) :: irec
        integer(ik4) :: itile , iti , itf
        iti = 1
        do itile = 1 , gdomain%ntiles
          istart(1) = gdomain%igstart(itile)
          icount(1) = gdomain%ni(itile)
          ! Latitudes are reversed in original file
          istart(2) = gdomain%jgstart
          icount(2) = gdomain%nj
          itf = iti + gdomain%ni(itile) - 1
          istatus = nf90_get_var(inet,ivar,work(iti:itf,:,:),istart,icount)
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error read var '//varname(irec))
          iti = iti + gdomain%ni(itile)
        end do
      end subroutine getwork
  end subroutine ein6hour

  subroutine headerein
    implicit none
    integer(ik4) :: i , j , k , kr
    integer(ik4) :: year , month , monthp1 , day , hour
    character(len=256) :: pathaddname
    integer(ik4) :: istatus , ncid , ivarid , idimid
    character(len=64) :: inname

    call split_idate(globidate1,year,month,day,hour)
    if ( dattyp == 'EIXXX' ) then
      monthp1 = month+1
      if ( monthp1 == 13 ) monthp1 = 1
      write(inname,'(a,i0.2,a,i0.2,a)') &
         't_xxxx',month,'0100-xxxx',monthp1,'0100.nc'
      pathaddname = trim(inpglob)//pthsep//'ERAIN_MEAN'//&
          pthsep//'XXXX'//pthsep//inname
      istatus = nf90_open(pathaddname,nf90_nowrite,ncid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error open file '//trim(pathaddname))
    else
      write(inname,'(i4,a,a,i4,a)') &
        year, pthsep, 'air.', year, '.00.nc'
      pathaddname = trim(inpglob)//pthsep//dattyp//pthsep//inname
      istatus = nf90_open(pathaddname,nf90_nowrite,ncid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error open file '//trim(pathaddname))
    end if
    istatus = nf90_inq_dimid(ncid,'latitude',idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Missing latitude dimension in file '//trim(pathaddname))
    istatus = nf90_inquire_dimension(ncid,idimid,len=jlat)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Error reading latitude dimelen in file '//trim(pathaddname))
    istatus = nf90_inq_dimid(ncid,'longitude',idimid)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Missing longitude dimension in file '//trim(pathaddname))
    istatus = nf90_inquire_dimension(ncid,idimid,len=ilon)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Error reading longitude dimelen in file '//trim(pathaddname))
    istatus = nf90_inq_dimid(ncid,'levelist',idimid)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_dimid(ncid,'level',idimid)
      call checkncerr(istatus,__FILE__,__LINE__, &
            'Missing level/levelist dimension in file '//trim(pathaddname))
    end if
    istatus = nf90_inquire_dimension(ncid,idimid,len=klev)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Error reading levelist dimelen in file '//trim(pathaddname))
    !
    ! Allocate working space
    !
    call getmem1d(plevs,1,klev,'mod_ein:plevs')
    call getmem1d(glat,1,jlat,'mod_ein:glat')
    call getmem1d(glon,1,ilon,'mod_ein:glon')
    call getmem1d(grev,1,max(jlat,ilon),'mod_ein:grev')
    call getmem1d(sigma1,1,klev,'mod_ein:sigma1')
    call getmem1d(sigmar,1,klev,'mod_ein:sigmar')
    call getmem3d(b3,1,jx,1,iy,1,klev*3,'mod_ein:b3')
    call getmem3d(d3,1,jx,1,iy,1,klev*2,'mod_ein:d3')

    istatus = nf90_inq_varid(ncid,'latitude',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Missing latitude variable in file '//trim(pathaddname))
    istatus = nf90_get_var(ncid,ivarid,glat)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Error reading latitude variable in file '//trim(pathaddname))
    istatus = nf90_inq_varid(ncid,'longitude',ivarid)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Missing longitude variable in file '//trim(pathaddname))
    istatus = nf90_get_var(ncid,ivarid,glon)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Error reading longitude variable in file '//trim(pathaddname))
    istatus = nf90_inq_varid(ncid,'levelist',ivarid)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_varid(ncid,'level',ivarid)
      call checkncerr(istatus,__FILE__,__LINE__, &
            'Missing level/levelist variable in file '//trim(pathaddname))
    end if
    istatus = nf90_get_var(ncid,ivarid,plevs)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Error reading levelist variable in file '//trim(pathaddname))
    istatus = nf90_close(ncid)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Error close file '//trim(pathaddname))
    sigmar(:) = real(plevs(:),rkx)/plevs(klev)
    pss = plevs(klev)/10.0_rkx ! mb -> cb
    !
    ! CHANGE ORDER OF VERTICAL INDEXES FOR PRESSURE LEVELS
    !
    do k = 1 , klev
      kr = klev - k + 1
      sigma1(k) = sigmar(kr)
    end do
    !
    ! Find window to read
    !
    call get_window(glat,glon,xlat,xlon,i_band,gdomain)

    grev(1:jlat) = glat
    jlat = gdomain%nj
    call getmem1d(glat,1,jlat,'mod_ein:glat')
    glat = grev(gdomain%jgstart:gdomain%jgstop)
    grev(1:ilon) = glon
    ilon = sum(gdomain%ni)
    call getmem1d(glon,1,ilon,'mod_ein:glon')
    glon(1:gdomain%ni(1)) = grev(gdomain%igstart(1):gdomain%igstop(1))
    if ( gdomain%ntiles == 2 ) then
      glon(gdomain%ni(1)+1:ilon) = grev(gdomain%igstart(2):gdomain%igstop(2))
    end if

    !call getmem2d(glat2,1,ilon,1,jlat,'mod_ein:glat2')
    !call getmem2d(glon2,1,ilon,1,jlat,'mod_ein:glon2')

    !do i = 1 , ilon
    !  glat2(i,:) = glat
    !end do
    !do j = 1 , jlat
    !  glon2(:,j) = glon
    !end do

    call getmem3d(b2,1,ilon,1,jlat,1,klev*3,'mod_ein:b2')
    call getmem3d(d2,1,ilon,1,jlat,1,klev*2,'mod_ein:d2')
    call getmem3d(work,1,ilon,1,jlat,1,klev,'mod_ein:work')
    !
    ! Set up pointers
    !
    u3 => d3(:,:,1:klev)
    v3 => d3(:,:,klev+1:2*klev)
    t3 => b3(:,:,1:klev)
    h3 => b3(:,:,klev+1:2*klev)
    q3 => b3(:,:,2*klev+1:3*klev)
    uvar => d2(:,:,1:klev)
    vvar => d2(:,:,klev+1:2*klev)
    tvar => b2(:,:,1:klev)
    hvar => b2(:,:,klev+1:2*klev)
    rhvar => b2(:,:,2*klev+1:3*klev)
  end subroutine headerein

end module mod_ein
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

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
  use mod_date
  use mod_constants
  use mod_memutil
  use mod_grid
  use mod_write
  use mod_vertint
  use mod_earth
  use mod_hgt
  use mod_humid
  use mod_mksst
  use mod_projections
  use mod_vectutil
  use mod_message
  use mod_nchelper
  use mod_kdinterp
  use netcdf

  private

  integer(ik4) :: jlat , ilon , klev , timlen

  real(rkx) , pointer , dimension(:,:,:) :: b3
  real(rkx) , pointer , dimension(:,:,:) :: d3
  real(rkx) , pointer , dimension(:,:,:) :: d3u , d3v

  real(rkx) , pointer :: u3(:,:,:) , v3(:,:,:)
  real(rkx) , pointer :: u3v(:,:,:) , v3u(:,:,:)
  real(rkx) , pointer :: h3(:,:,:) , q3(:,:,:) , t3(:,:,:)
  real(rkx) , pointer :: h3u(:,:,:) , h3v(:,:,:)
  real(rkx) , pointer :: uvar(:,:,:) , vvar(:,:,:)
  real(rkx) , pointer :: hvar(:,:,:) , qvar(:,:,:) , tvar(:,:,:)
  real(rkx) , pointer :: topou(:,:) , topov(:,:)

  real(rkx) , pointer , dimension(:,:,:) :: b2
  real(rkx) , pointer , dimension(:,:,:) :: d2
  real(rkx) , pointer , dimension(:) :: glat
  real(rkx) , pointer , dimension(:) :: grev
  real(rkx) , pointer , dimension(:) :: glon
  real(rkx) , pointer , dimension(:) :: plevs
  real(rkx) , pointer , dimension(:) :: sigmar
  real(rkx) :: pss , pst
  integer(2) , pointer , dimension(:,:,:) :: work

  integer , parameter :: numvars = 7
  integer , parameter :: numruns = 4

  integer(ik4) , dimension(numvars,numruns) :: ncfile
  integer(ik4) , dimension(numvars,numruns) :: varid
  real(rkx) , dimension(numvars,numruns) :: xoff , xscl
  type(rcm_time_and_date) , pointer , dimension(:) :: itimes
  integer(ik4) , pointer , dimension(:) :: xtimes
  logical :: lqas = .false.

  type(global_domain) :: gdomain
  type(h_interpolator) :: cross_hint , udot_hint , vdot_hint

  public :: init_ein , get_ein , conclude_ein

  contains

  subroutine init_ein
    implicit none
    integer(ik4) :: k
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
    call getmem1d(sigmar,1,klev,'mod_ein:sigmar')
    call getmem3d(b3,1,jx,1,iy,1,klev*3,'mod_ein:b3')
    if ( idynamic == 3 ) then
      call getmem3d(d3u,1,jx,1,iy,1,klev*2,'mod_ein:d3u')
      call getmem3d(d3v,1,jx,1,iy,1,klev*2,'mod_ein:d3v')
      call getmem3d(h3u,1,jx,1,iy,1,klev,'mod_ein:h3u')
      call getmem3d(h3v,1,jx,1,iy,1,klev,'mod_ein:h3v')
      call getmem2d(topou,1,jx,1,iy,'mod_ein:topou')
      call getmem2d(topov,1,jx,1,iy,'mod_ein:topov')
    else
      call getmem3d(d3,1,jx,1,iy,1,klev*2,'mod_ein:d3')
    end if

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
    do k = 1 , klev
      sigmar(k) = (plevs(klev-k+1)-plevs(1))/(plevs(klev)-plevs(1))
    end do
    pss = (plevs(klev)-plevs(1))/10.0_rkx ! mb -> cb
    pst = plevs(1)/10.0_rkx ! mb -> cb
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

    call h_interpolator_create(cross_hint,glat,glon,xlat,xlon)
    if ( idynamic == 3 ) then
      call h_interpolator_create(udot_hint,glat,glon,ulat,ulon)
      call h_interpolator_create(vdot_hint,glat,glon,vlat,vlon)
    else
      call h_interpolator_create(udot_hint,glat,glon,dlat,dlon)
    end if

    call getmem3d(b2,1,ilon,1,jlat,1,klev*3,'mod_ein:b2')
    call getmem3d(d2,1,ilon,1,jlat,1,klev*2,'mod_ein:d2')
    call getmem3d(work,1,ilon,1,jlat,1,klev,'mod_ein:work')
    !
    ! Set up pointers
    !
    if ( idynamic == 3 ) then
      u3 => d3u(:,:,1:klev)
      v3u => d3u(:,:,klev+1:2*klev)
      u3v => d3v(:,:,1:klev)
      v3 => d3v(:,:,klev+1:2*klev)
    else
      u3 => d3(:,:,1:klev)
      v3 => d3(:,:,klev+1:2*klev)
    end if
    t3 => b3(:,:,1:klev)
    h3 => b3(:,:,klev+1:2*klev)
    q3 => b3(:,:,2*klev+1:3*klev)
    uvar => d2(:,:,1:klev)
    vvar => d2(:,:,klev+1:2*klev)
    tvar => b2(:,:,1:klev)
    hvar => b2(:,:,klev+1:2*klev)
    qvar => b2(:,:,2*klev+1:3*klev)
    if ( idynamic == 3 ) then
      call ucrs2dot(zud4,z0,jx,iy,kz,i_band)
      call vcrs2dot(zvd4,z0,jx,iy,kz,i_crm)
      call ucrs2dot(topou,topogm,jx,iy,i_band)
      call vcrs2dot(topov,topogm,jx,iy,i_crm)
    end if
  end subroutine init_ein

  subroutine get_ein(idate)
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
    call h_interpolate_cont(cross_hint,b2,b3)
    if ( idynamic == 3 ) then
      call h_interpolate_cont(udot_hint,d2,d3u)
      call h_interpolate_cont(vdot_hint,d2,d3v)
    else
      call h_interpolate_cont(udot_hint,d2,d3)
    end if
    !
    ! Rotate u-v fields after horizontal interpolation
    !
    if ( idynamic == 3 ) then
      call pju%wind_rotate(u3,v3u)
      call pju%wind_rotate(u3v,v3)
    else
      call pjd%wind_rotate(u3,v3)
    end if
    !
    ! Invert vertical order top -> bottom for RegCM convention
    !
!$OMP SECTIONS
!$OMP SECTION
    call top2btm(t3)
!$OMP SECTION
    call top2btm(q3)
!$OMP SECTION
    call top2btm(h3)
!$OMP SECTION
    call top2btm(u3)
!$OMP SECTION
    call top2btm(v3)
!$OMP END SECTIONS
    !
    ! Vertical interpolation
    ! New calculation of p* on rcm topography.
    !
    if ( idynamic == 3 ) then
      call ucrs2dot(h3u,h3,jx,iy,klev,i_band)
      call vcrs2dot(h3v,h3,jx,iy,klev,i_crm)
      call intzps(ps4,topogm,t3,h3,pss,sigmar,pst, &
                  xlat,yeardayfrac(idate),jx,iy,klev)
      call intz3(ts4,t3,h3,topogm,jx,iy,klev,0.6_rkx,0.5_rkx,0.85_rkx)
    else
      call intgtb(pa,za,tlayer,topogm,t3,h3,pss,sigmar,pst,jx,iy,klev)
      call intpsn(ps4,topogm,pa,za,tlayer,ptop,jx,iy)
      call crs2dot(pd4,ps4,jx,iy,i_band,i_crm)
      call intv3(ts4,t3,ps4,pss,sigmar,ptop,pst,jx,iy,klev)
    end if

    call readsst(ts4,idate)
    !
    ! Interpolate U, V, T, and Q.
    !
    if ( idynamic == 3 ) then
!$OMP SECTIONS
!$OMP SECTION
      call intz1(u4,u3,zud4,h3u,topou,jx,iy,kz,klev,0.6_rkx,0.2_rkx,0.2_rkx)
!$OMP SECTION
      call intz1(v4,v3,zvd4,h3v,topov,jx,iy,kz,klev,0.6_rkx,0.2_rkx,0.2_rkx)
!$OMP SECTION
      call intz1(t4,t3,z0,h3,topogm,jx,iy,kz,klev,0.6_rkx,0.5_rkx,0.85_rkx)
!$OMP SECTION
      call intz1(q4,q3,z0,h3,topogm,jx,iy,kz,klev,0.7_rkx,0.4_rkx,0.7_rkx)
!$OMP END SECTIONS
    else
!$OMP SECTIONS
!$OMP SECTION
      call intv1(u4,u3,pd4,sigmah,pss,sigmar,ptop,pst,jx,iy,kz,klev,1)
!$OMP SECTION
      call intv1(v4,v3,pd4,sigmah,pss,sigmar,ptop,pst,jx,iy,kz,klev,1)
!$OMP SECTION
      call intv2(t4,t3,ps4,sigmah,pss,sigmar,ptop,pst,jx,iy,kz,klev)
!$OMP SECTION
      call intv1(q4,q3,ps4,sigmah,pss,sigmar,ptop,pst,jx,iy,kz,klev,1)
!$OMP END SECTIONS
    end if
    ! Get from RHUM to mixing ratio
    if ( lqas ) then
      q4 = d_10**q4
    else
      call rh2mxr(t4,q4,ps4,ptop,sigmah,jx,iy,kz)
    end if
  end subroutine get_ein

  subroutine ein6hour(dattyp,idate,idate0)
    implicit none
    character(len=5) , intent(in) :: dattyp
    type(rcm_time_and_date) , intent(in) :: idate , idate0
    integer(ik4) :: i , ifile , it , j , irun , inet , istatus , ivar
    integer(ik4) :: timid
    character(len=64) :: inname
    character(len=256) :: pathaddname
    character(len=4) , dimension(7) :: varname
    character(len=4) , dimension(7) :: fname
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
    data varname /'t' , 'z' , 'r' , 'u' , 'v' , 'q', 'clwc'/
    data fname   /'air','hgt','rhum','uwnd','vwnd', 'qas', 'qcs'/
    data hname   /'.00.','.06.','.12.','.18.'/

    irun = 1
    call split_idate(idate,year,month,day,hour)

    if ( dattyp == 'EIXXX' ) then
      if ( idate == idate0 .or. month /= lastmonth ) then
        lastmonth = month
        do inet = 1 , 5
          monthp1 = month+1
          if ( monthp1 == 13 ) monthp1 = 1
          write(inname,'(a,i0.2,a,i0.2,a)') &
             varname(inet)//'_xxxx',month,'0100-xxxx',monthp1,'0100.nc'
          pathaddname = trim(inpglob)//pthsep//'ERAIN_MEAN'//&
                        pthsep//'XXXX'//pthsep//inname
          istatus = nf90_open(pathaddname,nf90_nowrite,ncfile(inet,1))
          call checkncerr(istatus,__FILE__,__LINE__, &
            'Error open file '//trim(pathaddname))
          istatus = nf90_inq_varid(ncfile(inet,1),varname(inet), &
                                   varid(inet,1))
          call checkncerr(istatus,__FILE__,__LINE__, &
            'Error find var '//varname(inet))
          istatus = nf90_get_att(ncfile(inet,1),varid(inet,1), &
                   'scale_factor',xscl(inet,1))
          call checkncerr(istatus,__FILE__,__LINE__, &
            'Error find att scale_factor')
          istatus = nf90_get_att(ncfile(inet,1),varid(inet,1),  &
                   'add_offset',xoff(inet,1))
          call checkncerr(istatus,__FILE__,__LINE__, &
            'Error find att add_offset')
          write (stdout,*) ncfile(inet,1) , trim(pathaddname) ,   &
                           xscl(inet,1) , xoff(inet,1)
          if ( inet == 1 ) then
            istatus = nf90_inq_dimid(ncfile(1,1),'time',timid)
            call checkncerr(istatus,__FILE__,__LINE__, &
                            'Error find dim time')
            istatus = nf90_inquire_dimension(ncfile(1,1),timid,len=timlen)
            call checkncerr(istatus,__FILE__,__LINE__, &
                            'Error inquire time')
            istatus = nf90_inq_varid(ncfile(1,1),'time',timid)
            if ( istatus /= nf90_noerr ) then
              istatus = nf90_inq_varid(ncfile(1,1),'date',timid)
              call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error find var time/date')
            end if
            cunit = 'hours since 1900-01-01 00:00:00'
            ccal = 'noleap'
            call getmem1d(itimes,1,timlen,'mod_ein:itimes')
            call getmem1d(xtimes,1,timlen,'mod_ein:xtimes')
            istatus = nf90_get_var(ncfile(1,1),timid,xtimes)
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
        do irun = 1 , 4
          do inet = 1 , 5
            if ( inet == 3 ) then
              write(inname,'(i4,a,a,i4,a)') &
                year, pthsep, trim(fname(inet))//'.', year, hname(irun)//'nc'
              pathaddname = trim(inpglob)//pthsep//dattyp//pthsep//inname
              istatus = nf90_open(pathaddname,nf90_nowrite,ncfile(inet,irun))
              if ( istatus /= nf90_noerr ) then
                write(inname,'(i4,a,a,i4,a)') &
                  year, pthsep, trim(fname(6))//'.', year, hname(irun)//'nc'
                pathaddname = trim(inpglob)//pthsep//dattyp//pthsep//inname
                istatus = nf90_open(pathaddname,nf90_nowrite,ncfile(inet,irun))
                lqas = .true.
              end if
            else
              write(inname,'(i4,a,a,i4,a)') &
                  year, pthsep, trim(fname(inet))//'.', year, hname(irun)//'nc'
              pathaddname = trim(inpglob)//pthsep//dattyp//pthsep//inname
              istatus = nf90_open(pathaddname,nf90_nowrite,ncfile(inet,irun))
            end if
            call checkncerr(istatus,__FILE__,__LINE__, &
              'Error open file '//trim(pathaddname))
            if ( inet == 3 .and. lqas ) then
              istatus = nf90_inq_varid(ncfile(inet,irun),varname(6), &
                                       varid(inet,irun))
            else
              istatus = nf90_inq_varid(ncfile(inet,irun),varname(inet), &
                                       varid(inet,irun))
            end if
            call checkncerr(istatus,__FILE__,__LINE__, &
              'Error find var '//varname(inet))
            istatus = nf90_get_att(ncfile(inet,irun),varid(inet,irun), &
                     'scale_factor',xscl(inet,irun))
            call checkncerr(istatus,__FILE__,__LINE__, &
              'Error find att scale_factor')
            istatus = nf90_get_att(ncfile(inet,irun),varid(inet,irun),  &
                     'add_offset',xoff(inet,irun))
            call checkncerr(istatus,__FILE__,__LINE__, &
              'Error find att add_offset')
            write (stdout,*) ncfile(inet,irun) , trim(pathaddname) ,   &
                             xscl(inet,irun) , xoff(inet,irun)
            if ( irun == 1 .and. inet == 1 ) then
              istatus = nf90_inq_dimid(ncfile(1,1),'time',timid)
              call checkncerr(istatus,__FILE__,__LINE__, &
                              'Error find dim time')
              istatus = nf90_inquire_dimension(ncfile(1,1),timid,len=timlen)
              call checkncerr(istatus,__FILE__,__LINE__, &
                              'Error inquire time')
              istatus = nf90_inq_varid(ncfile(1,1),'time',timid)
              call checkncerr(istatus,__FILE__,__LINE__, &
                              'Error find var time')
              istatus = nf90_get_att(ncfile(1,1),timid,'units',cunit)
              call checkncerr(istatus,__FILE__,__LINE__, &
                              'Error read time units')
              ccal = 'gregorian'
              call getmem1d(itimes,1,timlen,'mod_ein:itimes')
              call getmem1d(xtimes,1,timlen,'mod_ein:xtimes')
              istatus = nf90_get_var(ncfile(1,1),timid,xtimes)
              call checkncerr(istatus,__FILE__,__LINE__, &
                              'Error read time')
              do it = 1 , timlen
                itimes(it) = timeval2date(real(xtimes(it),rkx),cunit,ccal)
              end do
            end if
          end do
        end do
      end if
      irun = hour/6 + 1
      tdif = idate - itimes(1)
      it = nint(tohours(tdif))/24 + 1
    end if

    istart(3) = 1
    icount(3) = klev
    istart(4) = it
    icount(4) = 1

    do inet = 1 , 5
      ifile = ncfile(inet,irun)
      ivar = varid(inet,irun)
      xscale = xscl(inet,irun)
      xadd = xoff(inet,irun)
      call getwork(inet)
      if ( inet == 1 ) then
        do j = 1 , jlat
          do i = 1 , ilon
            tvar(i,j,:) = real(real(work(i,j,:),rkx)*xscale+xadd,rkx)
          end do
        end do
      else if ( inet == 2 ) then
        do j = 1 , jlat
          do i = 1 , ilon
            hvar(i,j,:) = real(real(work(i,j,:),rkx) * &
                      xscale+xadd,rkx)/9.80616_rk4
          end do
        end do
      else if ( inet == 3 ) then
        do j = 1 , jlat
          do i = 1 , ilon
            qvar(i,j,:) = &
              max(real(real(work(i,j,:),rkx)*xscale+xadd,rkx),0.0_rkx)
          end do
        end do
        if ( lqas ) then
          call sph2mxr(qvar,ilon,jlat,klev)
          qvar = log10(qvar)
        else
          qvar = qvar * 0.01_rkx
        end if
      else if ( inet == 4 ) then
        do j = 1 , jlat
          do i = 1 , ilon
            uvar(i,j,:) = real(real(work(i,j,:),rkx)*xscale+xadd,rkx)
          end do
        end do
      else if ( inet == 5 ) then
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
          istatus = nf90_get_var(ifile,ivar,work(iti:itf,:,:),istart,icount)
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error read var '//varname(irec))
          iti = iti + gdomain%ni(itile)
        end do
      end subroutine getwork
  end subroutine ein6hour

  subroutine conclude_ein
    implicit none
    call h_interpolator_destroy(cross_hint)
    call h_interpolator_destroy(udot_hint)
    if ( idynamic == 3 ) then
      call h_interpolator_destroy(vdot_hint)
    end if
  end subroutine conclude_ein

end module mod_ein
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

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

module mod_era5rda

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_stdio
  use mod_memutil
  use mod_grid
  use mod_date
  use mod_constants
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

  integer(ik4) :: jlat, ilon, klev, timlen

  real(rkx), pointer, contiguous, dimension(:,:,:) :: b3
  real(rkx), pointer, contiguous, dimension(:,:,:) :: d3
  real(rkx), pointer, contiguous, dimension(:,:,:) :: d3u
  real(rkx), pointer, contiguous, dimension(:,:,:) :: d3v
  real(rkx), pointer, contiguous, dimension(:,:,:) :: b2
  real(rkx), pointer, contiguous, dimension(:,:,:) :: d2

  real(rkx), pointer, contiguous, dimension(:,:,:) :: u3, v3, h3, q3, t3
  real(rkx), pointer, contiguous, dimension(:,:,:) :: u3v, v3u, h3u, h3v
  real(rkx), pointer, contiguous, dimension(:,:,:) :: uvar, vvar
  real(rkx), pointer, contiguous, dimension(:,:,:) :: hvar, qvar, tvar
  real(rkx), pointer, contiguous, dimension(:,:) :: prvar, topou, topov

  real(rkx), pointer, contiguous, dimension(:) :: glat
  real(rkx), pointer, contiguous, dimension(:) :: grev
  real(rkx), pointer, contiguous, dimension(:) :: glon
  real(rkx), pointer, contiguous, dimension(:) :: plevs
  real(rkx), pointer, contiguous, dimension(:) :: sigmar
  real(rkx) :: pss, pst
  real(rkx), pointer, contiguous, dimension(:,:,:) :: work
  real(rkx), pointer, contiguous, dimension(:,:) :: iwork

  integer(ik4), dimension(5) :: inet5
  integer(ik4), dimension(5) :: ivar5
  real(rkx), dimension(5) :: xoff, xscl
  type(rcm_time_and_date), pointer, contiguous, dimension(:) :: itimes
  integer(ik4), pointer, contiguous, dimension(:) :: xtimes

  type(global_domain) :: gdomain
  type(h_interpolator) :: cross_hint, udot_hint, vdot_hint

  public :: init_era5rda, get_era5rda, conclude_era5rda

  contains

  subroutine init_era5rda
    implicit none
    integer(ik4) :: k
    integer(ik4) :: year, month, day, hour
    character(len=256) :: pathaddname
    integer(ik4) :: istatus, ncid, ivarid, idimid
    character(len=64) :: inname

    ! set the subdirectory path to 3D data: ds633.0/e5.oper.an.pl
    character(len=21) :: subdir_plev = 'ds633.0' // pthsep // 'e5.oper.an.pl'
    character(len=6) :: yyyymm

    call split_idate(globidate1,year,month,day,hour)

    ! store the year and month in a string of format YYYYMM
    write(yyyymm,'(i4.4,i2.2)') year, month

    ! set the file name to be read
    !(format is e5.oper.an.pl.128_129_z.ll025sc.1979010100_1979010123.nc)
    write(inname,'(a,i0.4,i0.2,i0.2,a,i0.4,i0.2,i0.2,a)') &
        'e5.oper.an.pl.128_129_z.ll025sc.', year, month, day,'00_', &
         year, month, day, '23.nc'

    ! set the path to the file to be read
    pathaddname = trim(inpglob)//pthsep//trim(dattyp)//pthsep//subdir_plev//pthsep//yyyymm//pthsep//inname

    istatus = nf90_open(pathaddname,nf90_nowrite,ncid)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error open file '//trim(pathaddname))
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
          'Error reading longitude dimlen in file '//trim(pathaddname))
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
    call getmem1d(plevs,1,klev,'mod_era5:plevs')
    call getmem1d(glat,1,jlat,'mod_era5:glat')
    call getmem1d(glon,1,ilon,'mod_era5:glon')
    call getmem1d(grev,1,max(jlat,ilon),'mod_era5:grev')
    call getmem1d(sigmar,1,klev,'mod_era5:sigmar')
    call getmem3d(b3,1,jx,1,iy,1,klev*3,'mod_era5:b3')
    if ( idynamic == 3 ) then
      call getmem3d(d3u,1,jx,1,iy,1,klev*2,'mod_era5:d3u')
      call getmem3d(d3v,1,jx,1,iy,1,klev*2,'mod_era5:d3v')
      call getmem3d(h3u,1,jx,1,iy,1,klev,'mod_era5:h3u')
      call getmem3d(h3v,1,jx,1,iy,1,klev,'mod_era5:h3v')
      call getmem2d(topou,1,jx,1,iy,'mod_era5:topou')
      call getmem2d(topov,1,jx,1,iy,'mod_era5:topov')
    else
      call getmem3d(d3,1,jx,1,iy,1,klev*2,'mod_era5:d3')
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
    do k = 1, klev
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
    call getmem1d(glat,1,jlat,'mod_era5:glat')
    glat = grev(gdomain%jgstart:gdomain%jgstop)
    grev(1:ilon) = glon
    ilon = sum(gdomain%ni)
    call getmem1d(glon,1,ilon,'mod_era5:glon')
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

    call getmem3d(b2,1,ilon,1,jlat,1,klev*3,'mod_era5:b2')
    call getmem3d(d2,1,ilon,1,jlat,1,klev*2,'mod_era5:d2')
    call getmem3d(work,1,ilon,1,jlat,1,klev,'mod_era5:work')
    !
    ! Set up pointer
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
  end subroutine init_era5rda

  subroutine get_era5rda(idate)
    implicit none
    type(rcm_time_and_date), intent(in) :: idate
    !
    ! Read data at idate
    !
    call read_era5rda(dattyp,idate,globidate1)
    write (stdout,*) 'READ IN fields at DATE:', tochar(idate)
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
      call pjv%wind_rotate(u3v,v3)
    else
      call pjd%wind_rotate(u3,v3)
    end if
    !
    ! Invert vertical order, set BOTTOM -> TOP
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
    q4 = d_10**q4
  end subroutine get_era5rda

  subroutine read_era5rda(dattyp,idate,idate0)
    implicit none
    character(len=5), intent(in) :: dattyp
    type(rcm_time_and_date), intent(in) :: idate, idate0
    integer(ik4) :: i, inet, it, j, kkrec, istatus, ivar
    integer(ik4) :: timid
    character(len=64) :: inname
    character(len=256) :: pathaddname
    character(len=1), dimension(5) :: varname
    character(len=33), dimension(5) :: fname
    character(len=64) :: cunit, ccal
    real(rkx) :: xadd, xscale
    integer(ik4), dimension(4) :: icount, istart
    integer(ik4) :: year, month, day, hour
    integer(ik4), save :: lastday
    type(rcm_time_interval) :: tdif
    character(len=21) :: subdir_plev = 'ds633.0' // pthsep // 'e5.oper.an.pl'
    character(len=6) :: yyyymm
    !
    ! This is the latitude, longitude dimension of the grid to be read.
    ! This corresponds to the lat and lon dimension variables in the
    ! netCDF file.
    ! The data are packed into short integers (INTEGER*2).  The array
    ! work will be used to hold the packed integers.  The array 'x'
    ! will contain the unpacked data.
    !
    data varname /    &
                  'T',&
                  'Z',&
                  'Q',&
                  'U',&
                  'V' &
                 /
    data fname   /    &
                  'e5.oper.an.pl.128_130_t.ll025sc.',&
                  'e5.oper.an.pl.128_129_z.ll025sc.',&
                  'e5.oper.an.pl.128_133_q.ll025sc.',&
                  'e5.oper.an.pl.128_131_u.ll025uv.',&
                  'e5.oper.an.pl.128_132_v.ll025uv.' &
                 /

    call split_idate(idate,year,month,day,hour)

    ! store the year and month in a string of format YYYYMM
    write(yyyymm,'(i4.4,i2.2)') year, month

    ! determine whether this is a new day; open the next file if so
    if ( idate == idate0 .or. day /= lastday ) then
      lastday = day
      if ( idate /= idate0 ) then
        do kkrec = 1, 5
          istatus = nf90_close(inet5(kkrec))
          call checkncerr(istatus,__FILE__,__LINE__, &
              'Error close file')
        end do
      end if
      do kkrec = 1, 5
        ! set the file name to be read
        ! (format is e5.oper.an.pl.128_129_z.ll025sc.1979010100_1979010123.nc)
        write(inname,'(a,i0.4,i0.2,i0.2,a,i0.4,i0.2,i0.2,a)') &
        trim(fname(kkrec)), year, month, day,'00_', &
        year, month, day, '23.nc'

        ! set the path to the file to be read
        pathaddname = trim(inpglob)//pthsep//trim(dattyp)//pthsep//subdir_plev//pthsep//yyyymm//pthsep//inname

        ! open the file
        istatus = nf90_open(pathaddname,nf90_nowrite,inet5(kkrec))
        call checkncerr(istatus,__FILE__,__LINE__, &
          'Error open file '//trim(pathaddname))

        ! get the current variable
        istatus = nf90_inq_varid(inet5(kkrec),varname(kkrec), &
                                 ivar5(kkrec))
        call checkncerr(istatus,__FILE__,__LINE__, &
          'Error find var '//varname(kkrec))

        ! get the scale factor
        istatus = nf90_get_att(inet5(kkrec),ivar5(kkrec), &
                 'scale_factor',xscl(kkrec))
        ! assume a scale factor of 1 if it's missing
        if ( istatus /= 0 ) then
            xscl(kkrec) = 1.0
        end if

        ! get the offset
        istatus = nf90_get_att(inet5(kkrec),ivar5(kkrec),  &
                   'add_offset',xoff(kkrec))
        ! assume an offset of 0 if it's missing
        if ( istatus /= 0 ) then
            xoff(kkrec) = 0.0
        end if
        
        write (stdout,*) inet5(kkrec), trim(pathaddname),   &
                         xscl(kkrec), xoff(kkrec)
        ! if this is the first variable, get the time and metadata
        if ( kkrec == 1 ) then
          ! get the time dimension/variable
          istatus = nf90_inq_dimid(inet5(1),'time',timid)
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error find dim time')
          istatus = nf90_inquire_dimension(inet5(1),timid,len=timlen)
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error inquire time')
          istatus = nf90_inq_varid(inet5(1),'time',timid)
          call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find var time')

          ! get the units and calendar
          istatus = nf90_get_att(inet5(1),timid,'units',cunit)
          call checkncerr(istatus,__FILE__,__LINE__, &
                              'Error read time units')
          istatus = nf90_get_att(inet5(1),timid,'calendar',ccal)
          call checkncerr(istatus,__FILE__,__LINE__, &
                              'Error read time units')
          call getmem1d(itimes,1,timlen,'mod_era5:itimes')
          call getmem1d(xtimes,1,timlen,'mod_era5:xtimes')
          istatus = nf90_get_var(inet5(1),timid,xtimes)
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error read time')
          ! convert the time to a date
          do it = 1, timlen
            itimes(it) = timeval2date(real(xtimes(it),rkx),cunit,ccal)
          end do
        end if
      end do
    end if
    ! determine the time index of the current step
    it = hour + 1

    ! set the start and count for the read
    istart(3) = 1
    icount(3) = klev
    istart(4) = it
    icount(4) = 1

    do kkrec = 1, 5
      inet = inet5(kkrec)
      ivar = ivar5(kkrec)
      xscale = xscl(kkrec)
      xadd = xoff(kkrec)
      call getwork(kkrec)
      ! read temperature
      if ( kkrec == 1 ) then
        do j = 1, jlat
          do i = 1, ilon
            tvar(i,j,:) = real(real(work(i,j,:),rkx)*xscale+xadd,rkx)
          end do
        end do
      ! read geopotential
      else if ( kkrec == 2 ) then
        do j = 1, jlat
          do i = 1, ilon
            hvar(i,j,:) = real(real(work(i,j,:),rkx) * &
                      xscale+xadd,rkx)/9.80616_rk4
          end do
        end do
      ! read specific humidity
      else if ( kkrec == 3 ) then
        do j = 1, jlat
          do i = 1, ilon
            qvar(i,j,:) = &
              max(real(real(work(i,j,:),rkx)*xscale+xadd,rkx),0.0_rkx)
          end do
        end do
        call sph2mxr(qvar,ilon,jlat,klev)
        qvar = log10(max(qvar,dlowval))
      ! read u-wind
      else if ( kkrec == 4 ) then
        do j = 1, jlat
          do i = 1, ilon
            uvar(i,j,:) = real(real(work(i,j,:),rkx)*xscale+xadd,rkx)
          end do
        end do
      ! read v-wind
      else if ( kkrec == 5 ) then
        do j = 1, jlat
          do i = 1, ilon
            vvar(i,j,:) = real(real(work(i,j,:),rkx)*xscale+xadd,rkx)
          end do
        end do
      end if
    end do

    contains

      subroutine getwork(irec)
        implicit none
        integer(ik4), intent(in) :: irec
        integer(ik4) :: itile, iti, itf
        iti = 1
        do itile = 1, gdomain%ntiles
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
  end subroutine read_era5rda

  subroutine conclude_era5rda
    implicit none
    call h_interpolator_destroy(cross_hint)
    call h_interpolator_destroy(udot_hint)
    if ( idynamic == 3 ) then
      call h_interpolator_destroy(vdot_hint)
    end if
  end subroutine conclude_era5rda


end module mod_era5rda
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

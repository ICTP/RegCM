!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_era5

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
  real(rkx), pointer, contiguous, dimension(:,:,:) :: uvar, vvar, hvar, qvar, tvar
  real(rkx), pointer, contiguous, dimension(:,:) :: prvar, topou, topov

  real(rkx), pointer, contiguous, dimension(:) :: glat
  real(rkx), pointer, contiguous, dimension(:) :: grev
  real(rkx), pointer, contiguous, dimension(:) :: glon
  real(rkx), pointer, contiguous, dimension(:) :: plevs
  real(rkx), pointer, contiguous, dimension(:) :: sigmar
  real(rkx) :: pss, pst
  integer(2), pointer, contiguous, dimension(:,:,:) :: iwork3
  real(rkx), pointer, contiguous, dimension(:,:,:) :: rwork3
  integer(2), pointer, contiguous, dimension(:,:) :: iwork
  real(rkx), pointer, contiguous, dimension(:,:) :: rwork

  integer(ik4), dimension(5) :: inet5
  integer(ik4), dimension(5) :: ivar5
  real(rkx), dimension(5) :: xoff, xscl
  type(rcm_time_and_date), pointer, contiguous, dimension(:) :: itimes
  integer(ik4), pointer, contiguous, dimension(:) :: xtimes

  type(global_domain) :: gdomain
  type(h_interpolator) :: cross_hint, udot_hint, vdot_hint

  public :: init_era5, get_era5, conclude_era5
  public :: init_era5h, get_era5h, conclude_era5h

  contains

  subroutine init_era5h
    implicit none
    integer(ik4) :: year, month, day, hour
    character(len=256) :: pathaddname
    integer(ik4) :: istatus, ncid, ivarid, idimid
    character(len=64) :: inname

    call split_idate(globidate1,year,month,day,hour)
    if ( dattyp == 'ERAXX' ) then
      write(inname,'(a,i0.2,a)') 'pr_XXXX_', month,'.nc'
      pathaddname = trim(inpglob)//pthsep//'ERA5_MEAN'//pthsep// &
         'hourly'//pthsep//inname
    else
      write(inname,'(a,i0.4,a,i0.2,a)') 'pr_', year, '_', month,'.nc'
      pathaddname = trim(inpglob)//pthsep//'ERA5'//pthsep// &
         'hourly'//pthsep//inname
    end if
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
          'Error reading longitude dimelen in file '//trim(pathaddname))
    !
    ! Allocate working space
    !
    call getmem1d(glat,1,jlat,'mod_era5:glat')
    call getmem1d(glon,1,ilon,'mod_era5:glon')
    call getmem1d(grev,1,max(jlat,ilon),'mod_era5:grev')
    call getmem3d(b3,1,jx,1,iy,1,3,'mod_era5:b3')

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
    istatus = nf90_close(ncid)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Error close file '//trim(pathaddname))
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
    call getmem3d(b2,1,ilon,1,jlat,1,3,'mod_era5:b2')
    call getmem2d(prvar,1,ilon,1,jlat,'mod_era5:prvar')
    call getmem2d(iwork,1,ilon,1,jlat,'mod_era5:iwork')
    call getmem2d(rwork,1,ilon,1,jlat,'mod_era5:rwork')
  end subroutine init_era5h

  subroutine init_era5
    implicit none
    integer(ik4) :: k
    integer(ik4) :: year, month, day, hour
    character(len=256) :: pathaddname
    integer(ik4) :: istatus, ncid, ivarid, idimid
    character(len=64) :: inname

    call split_idate(globidate1,year,month,day,hour)
    if ( dattyp == 'ERAXX' ) then
      write(inname,'(a,a,a,a,a,i0.2,a)') &
          'XXXX', pthsep, 'geop_', 'XXXX', '_', month,'.nc'
      pathaddname = trim(inpglob)//pthsep//'ERA5_MEAN'//pthsep//inname
    else
      write(inname,'(i4,a,a,i0.4,a,i0.2,a)') &
      year, pthsep, 'geop_', year, '_', month,'.nc'
      pathaddname = trim(inpglob)//pthsep//dattyp(1:4)//pthsep//inname
    end if
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
          'Error reading longitude dimelen in file '//trim(pathaddname))
    istatus = nf90_inq_dimid(ncid,'levelist',idimid)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_dimid(ncid,'level',idimid)
      if ( istatus /= nf90_noerr ) then
        istatus = nf90_inq_dimid(ncid,'pressure_level',idimid)
        call checkncerr(istatus,__FILE__,__LINE__, &
              'Missing level/levelist dimension in file '//trim(pathaddname))
      end if
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
      if ( istatus /= nf90_noerr ) then
        istatus = nf90_inq_varid(ncid,'pressure_level',ivarid)
        call checkncerr(istatus,__FILE__,__LINE__, &
            'Missing level/levelist variable in file '//trim(pathaddname))
      end if
    end if
    istatus = nf90_get_var(ncid,ivarid,plevs)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Error reading levelist variable in file '//trim(pathaddname))
    istatus = nf90_close(ncid)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Error close file '//trim(pathaddname))
    if ( plevs(1) > plevs(klev) ) then
      do k = 1, klev
        sigmar(k) = (plevs(k)-plevs(klev))/(plevs(1)-plevs(klev))
      end do
      pss = (plevs(1)-plevs(klev))/10.0_rkx ! mb -> cb
      pst = plevs(klev)/10.0_rkx ! mb -> cb
    else
      do k = 1, klev
        sigmar(k) = (plevs(klev-k+1)-plevs(1))/(plevs(klev)-plevs(1))
      end do
      pss = (plevs(klev)-plevs(1))/10.0_rkx ! mb -> cb
      pst = plevs(1)/10.0_rkx ! mb -> cb
    end if
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
    call getmem3d(iwork3,1,ilon,1,jlat,1,klev,'mod_era5:iwork3')
    call getmem3d(rwork3,1,ilon,1,jlat,1,klev,'mod_era5:rwork3')
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
  end subroutine init_era5

  subroutine get_era5h(idate)
    implicit none
    type(rcm_time_and_date), intent(in) :: idate
    call era5hour(idate,globidate1)
    write (stdout,*) 'READ IN fields at DATE:', tochar(idate)
    call h_interpolate_nn(cross_hint,prvar,pr)
    call h_interpolate_cont(cross_hint,b2,b3)
    pr(:,:) = pr(:,:) / secph * 0.001_rkx
    ssr(:,:) = b3(:,:,1) / secph
    strd(:,:) = b3(:,:,2) / secph
    clt(:,:) = b3(:,:,3)
  end subroutine get_era5h

  subroutine get_era5(idate)
    implicit none
    type(rcm_time_and_date), intent(in) :: idate
    !
    ! Read data at idate
    !
    call era56hour(dattyp,idate,globidate1)
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
    if ( h3(1,1,1) > h3(1,1,klev) ) then
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
    end if
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
  end subroutine get_era5

  subroutine era5hour(idate,idate0)
    implicit none
    type(rcm_time_and_date), intent(in) :: idate, idate0
    integer(ik4) :: i, inet, it, j, kkrec, istatus, ivar
    integer(ik4) :: timid
    character(len=64) :: inname
    character(len=256) :: pathaddname
    character(len=4), dimension(4) :: varname
    character(len=4), dimension(4) :: fname
    character(len=64) :: cunit, ccal
    real(rkx) :: xadd, xscale
    integer(ik4), dimension(3) :: icount, istart
    integer(ik4) :: year, month, day, hour
    integer(ik4), save :: lastmonth
    type(rcm_time_interval) :: tdif
    logical, save :: has_offset = .false.
    !
    ! This is the latitude, longitude dimension of the grid to be read.
    ! This corresponds to the lat and lon dimension variables in the
    ! netCDF file.
    ! The data are packed into short integers (INTEGER*2).  The array
    ! work will be used to hold the packed integers.  The array 'x'
    ! will contain the unpacked data.
    !
    data varname /'tp  ', 'ssr ', 'strd', 'tcc '/
    data fname   /'pr  ', 'ssr ', 'strd', 'clt '/

    call split_idate(idate,year,month,day,hour)

    if ( idate == idate0 .or. month /= lastmonth ) then
      lastmonth = month
      if ( idate /= idate0 ) then
        do kkrec = 1, 4
          istatus = nf90_close(inet5(kkrec))
          call checkncerr(istatus,__FILE__,__LINE__, &
              'Error close file')
        end do
      end if
      do kkrec = 1, 4
        if ( dattyp == 'ERAXX' ) then
          write(inname,'(a,a,i0.2,a)') &
            trim(fname(kkrec)), '_XXXX_', month,'.nc'
          pathaddname = trim(inpglob)//pthsep//'ERA5_MEAN'//pthsep// &
              pthsep//'hourly'//pthsep//inname
        else
          write(inname,'(a,a,i0.4,a,i0.2,a)') &
            trim(fname(kkrec)), '_', year, '_', month,'.nc'
          pathaddname = trim(inpglob)//pthsep//'ERA5'//pthsep// &
              pthsep//'hourly'//pthsep//inname
        end if
        istatus = nf90_open(pathaddname,nf90_nowrite,inet5(kkrec))
        call checkncerr(istatus,__FILE__,__LINE__, &
          'Error open file '//trim(pathaddname))
        istatus = nf90_inq_varid(inet5(kkrec),trim(varname(kkrec)), &
                                 ivar5(kkrec))
        call checkncerr(istatus,__FILE__,__LINE__, &
          'Error find var '//varname(kkrec))
        istatus = nf90_get_att(inet5(kkrec),ivar5(kkrec), &
                 'scale_factor',xscl(kkrec))
        if ( istatus == nf90_noerr ) then
          istatus = nf90_get_att(inet5(kkrec),ivar5(kkrec),  &
                     'add_offset',xoff(kkrec))
          call checkncerr(istatus,__FILE__,__LINE__, &
              'Error find att add_offset')
          has_offset = .true.
          write (stdout,*) inet5(kkrec), trim(pathaddname),   &
                           xscl(kkrec), xoff(kkrec)
        else
          has_offset = .false.
          write (stdout,*) inet5(kkrec), trim(pathaddname)
        end if
        if ( kkrec == 1 ) then
          if ( dattyp == 'ERAXX' ) then
            call getmem1d(itimes,1,1,'mod_era5:itimes')
            itimes(1) = year*1000000 + month*10000+100
            call setcal(itimes(1),'noleap')
          else
            istatus = nf90_inq_dimid(inet5(1),'time',timid)
            if ( istatus /= nf90_noerr ) then
              istatus = nf90_inq_dimid(inet5(1),'valid_time',timid)
              call checkncerr(istatus,__FILE__,__LINE__, &
                              'Error find dim time')
            end if
            istatus = nf90_inquire_dimension(inet5(1),timid,len=timlen)
            call checkncerr(istatus,__FILE__,__LINE__, &
                            'Error inquire time')
            istatus = nf90_inq_varid(inet5(1),'time',timid)
            if ( istatus /= nf90_noerr ) then
              istatus = nf90_inq_varid(inet5(1),'valid_time',timid)
              if ( istatus /= nf90_noerr ) then
                istatus = nf90_inq_varid(inet5(1),'date',timid)
                call checkncerr(istatus,__FILE__,__LINE__, &
                                'Error find var time/date')
              end if
            end if
            istatus = nf90_myget_att_text(inet5(1),timid,'units',cunit)
            call checkncerr(istatus,__FILE__,__LINE__, &
                            'Error read time units')
            istatus = nf90_myget_att_text(inet5(1),timid,'calendar',ccal)
            call checkncerr(istatus,__FILE__,__LINE__, &
                            'Error read time units')
            call getmem1d(itimes,1,timlen,'mod_era5:itimes')
            call getmem1d(xtimes,1,timlen,'mod_era5:xtimes')
            istatus = nf90_get_var(inet5(1),timid,xtimes)
            call checkncerr(istatus,__FILE__,__LINE__, &
                            'Error read time')
            do it = 1, timlen
              itimes(it) = timeval2date(real(xtimes(it),rkx),cunit,ccal)
            end do
          end if
        end if
      end do
    end if

    tdif = idate - itimes(1)
    it = nint(tohours(tdif)) + 1

    istart(3) = it
    icount(3) = 1

    do kkrec = 1, 4
      inet = inet5(kkrec)
      ivar = ivar5(kkrec)
      xscale = xscl(kkrec)
      xadd = xoff(kkrec)
      if ( has_offset ) then
        call igetwork(kkrec)
        if ( kkrec == 1 ) then
          do j = 1, jlat
            do i = 1, ilon
              prvar(i,j) = real(real(iwork(i,j),rkx)*xscale+xadd,rkx)
            end do
          end do
        else
          do j = 1, jlat
            do i = 1, ilon
              b2(i,j,kkrec-1) = real(real(iwork(i,j),rkx)*xscale+xadd,rkx)
            end do
          end do
        end if
      else
        call rgetwork(kkrec)
        if ( kkrec == 1 ) then
          do j = 1, jlat
            do i = 1, ilon
              prvar(i,j) = rwork(i,j)
            end do
          end do
        else
          do j = 1, jlat
            do i = 1, ilon
              b2(i,j,kkrec-1) = rwork(i,j)
            end do
          end do
        end if
      end if
    end do

    contains

      subroutine igetwork(irec)
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
          istatus = nf90_get_var(inet,ivar,iwork(iti:itf,:),istart,icount)
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error read var '//varname(irec))
          iti = iti + gdomain%ni(itile)
        end do
      end subroutine igetwork

      subroutine rgetwork(irec)
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
          istatus = nf90_get_var(inet,ivar,rwork(iti:itf,:),istart,icount)
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error read var '//varname(irec))
          iti = iti + gdomain%ni(itile)
        end do
      end subroutine rgetwork

  end subroutine era5hour

  subroutine era56hour(dattyp,idate,idate0)
    implicit none
    character(len=5), intent(in) :: dattyp
    type(rcm_time_and_date), intent(in) :: idate, idate0
    integer(ik4) :: i, inet, it, j, kkrec, istatus, ivar
    integer(ik4) :: timid
    character(len=64) :: inname
    character(len=256) :: pathaddname
    character(len=1), dimension(5) :: varname
    character(len=4), dimension(5) :: fname
    character(len=64) :: cunit, ccal
    real(rkx) :: xadd, xscale
    integer(ik4), dimension(4) :: icount, istart
    integer(ik4) :: year, month, day, hour
    integer(ik4), save :: lastmonth
    type(rcm_time_interval) :: tdif
    logical, save :: has_offset = .false.
    !
    ! This is the latitude, longitude dimension of the grid to be read.
    ! This corresponds to the lat and lon dimension variables in the
    ! netCDF file.
    ! The data are packed into short integers (INTEGER*2).  The array
    ! work will be used to hold the packed integers.  The array 'x'
    ! will contain the unpacked data.
    !
    data varname /'t', 'z', 'q', 'u', 'v'/
    data fname   /'tatm','geop','qhum','uwnd','vwnd'/

    call split_idate(idate,year,month,day,hour)

    if ( idate == idate0 .or. month /= lastmonth ) then
      lastmonth = month
      if ( idate /= idate0 ) then
        do kkrec = 1, 5
          istatus = nf90_close(inet5(kkrec))
          call checkncerr(istatus,__FILE__,__LINE__, &
              'Error close file')
        end do
      end if
      do kkrec = 1, 5
        if ( dattyp == 'ERAXX' ) then
          write(inname,'(a,a,a,a,a,a,i0.2,a)') &
          'XXXX', pthsep, fname(kkrec), '_', 'XXXX', '_', month,'.nc'
          pathaddname = trim(inpglob)//pthsep//'ERA5_MEAN'//pthsep//inname
        else
          write(inname,'(i4,a,a,a,i0.4,a,i0.2,a)') &
          year, pthsep, fname(kkrec), '_', year, '_', month,'.nc'
          pathaddname = trim(inpglob)//pthsep//dattyp(1:4)//pthsep//inname
        end if
        istatus = nf90_open(pathaddname,nf90_nowrite,inet5(kkrec))
        call checkncerr(istatus,__FILE__,__LINE__, &
          'Error open file '//trim(pathaddname))
        istatus = nf90_inq_varid(inet5(kkrec),varname(kkrec), &
                                 ivar5(kkrec))
        call checkncerr(istatus,__FILE__,__LINE__, &
          'Error find var '//varname(kkrec))
        istatus = nf90_get_att(inet5(kkrec),ivar5(kkrec), &
                 'scale_factor',xscl(kkrec))
        if ( istatus == nf90_noerr ) then
          istatus = nf90_get_att(inet5(kkrec),ivar5(kkrec),  &
                     'add_offset',xoff(kkrec))
          call checkncerr(istatus,__FILE__,__LINE__, &
              'Error find att add_offset')
          has_offset = .true.
          write (stdout,*) inet5(kkrec), trim(pathaddname),   &
                           xscl(kkrec), xoff(kkrec)
        else
          has_offset = .false.
          write (stdout,*) inet5(kkrec), trim(pathaddname)
        end if
        if ( kkrec == 1 ) then
          if ( dattyp == 'ERAXX' ) then
            call getmem1d(itimes,1,1,'mod_era5:itimes')
            itimes(1) = year*1000000 + month*10000+100
            call setcal(itimes(1),'noleap')
          else
            istatus = nf90_inq_dimid(inet5(1),'time',timid)
            if ( istatus /= nf90_noerr ) then
              istatus = nf90_inq_dimid(inet5(1),'valid_time',timid)
              call checkncerr(istatus,__FILE__,__LINE__, &
                              'Error find dim time')
            end if
            istatus = nf90_inquire_dimension(inet5(1),timid,len=timlen)
            call checkncerr(istatus,__FILE__,__LINE__, &
                            'Error inquire time')
            istatus = nf90_inq_varid(inet5(1),'time',timid)
            if ( istatus /= nf90_noerr ) then
              istatus = nf90_inq_varid(inet5(1),'valid_time',timid)
              if ( istatus /= nf90_noerr ) then
                istatus = nf90_inq_varid(inet5(1),'date',timid)
                call checkncerr(istatus,__FILE__,__LINE__, &
                                'Error find var time/date')
              end if
            end if
            istatus = nf90_myget_att_text(inet5(1),timid,'units',cunit)
            call checkncerr(istatus,__FILE__,__LINE__, &
                                'Error read time units')
            istatus = nf90_myget_att_text(inet5(1),timid,'calendar',ccal)
            call checkncerr(istatus,__FILE__,__LINE__, &
                                'Error read time units')
            call getmem1d(itimes,1,timlen,'mod_era5:itimes')
            call getmem1d(xtimes,1,timlen,'mod_era5:xtimes')
            istatus = nf90_get_var(inet5(1),timid,xtimes)
            call checkncerr(istatus,__FILE__,__LINE__, &
                            'Error read time')
            do it = 1, timlen
              itimes(it) = timeval2date(real(xtimes(it),rkx),cunit,ccal)
            end do
          end if
        end if
      end do
    end if
    tdif = idate - itimes(1)
    it = nint(tohours(tdif))/ibdyfrq + 1

    istart(3) = 1
    icount(3) = klev
    istart(4) = it
    icount(4) = 1

    if ( has_offset ) then
      do kkrec = 1, 5
        inet = inet5(kkrec)
        ivar = ivar5(kkrec)
        xscale = xscl(kkrec)
        xadd = xoff(kkrec)
        call igetwork(kkrec)
        if ( kkrec == 1 ) then
          do j = 1, jlat
            do i = 1, ilon
              tvar(i,j,:) = real(real(iwork3(i,j,:),rkx)*xscale+xadd,rkx)
            end do
          end do
        else if ( kkrec == 2 ) then
          do j = 1, jlat
            do i = 1, ilon
              hvar(i,j,:) = real(real(iwork3(i,j,:),rkx) * &
                        xscale+xadd,rkx)/9.80616_rk4
            end do
          end do
        else if ( kkrec == 3 ) then
          do j = 1, jlat
            do i = 1, ilon
              qvar(i,j,:) = &
                max(real(real(iwork3(i,j,:),rkx)*xscale+xadd,rkx),0.0_rkx)
            end do
          end do
          call sph2mxr(qvar,ilon,jlat,klev)
          qvar = log10(max(qvar,dlowval))
        else if ( kkrec == 4 ) then
          do j = 1, jlat
            do i = 1, ilon
              uvar(i,j,:) = real(real(iwork3(i,j,:),rkx)*xscale+xadd,rkx)
            end do
          end do
        else if ( kkrec == 5 ) then
          do j = 1, jlat
            do i = 1, ilon
              vvar(i,j,:) = real(real(iwork3(i,j,:),rkx)*xscale+xadd,rkx)
            end do
          end do
        end if
      end do
    else
      do kkrec = 1, 5
        inet = inet5(kkrec)
        ivar = ivar5(kkrec)
        call rgetwork(kkrec)
        if ( kkrec == 1 ) then
          tvar = rwork3
        else if ( kkrec == 2 ) then
          hvar = rwork3/9.80616_rkx
        else if ( kkrec == 3 ) then
          qvar = max(rwork3,0.0_rkx)
          call sph2mxr(qvar,ilon,jlat,klev)
          qvar = log10(max(qvar,dlowval))
        else if ( kkrec == 4 ) then
          uvar = rwork3
        else if ( kkrec == 5 ) then
          vvar = rwork3
        end if
      end do
    end if

    contains

      subroutine igetwork(irec)
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
          istatus = nf90_get_var(inet,ivar,iwork3(iti:itf,:,:),istart,icount)
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error read var '//varname(irec))
          iti = iti + gdomain%ni(itile)
        end do
      end subroutine igetwork

      subroutine rgetwork(irec)
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
          istatus = nf90_get_var(inet,ivar,rwork3(iti:itf,:,:),istart,icount)
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error read var '//varname(irec))
          iti = iti + gdomain%ni(itile)
        end do
      end subroutine rgetwork

  end subroutine era56hour

  subroutine conclude_era5
    implicit none
    call h_interpolator_destroy(cross_hint)
    call h_interpolator_destroy(udot_hint)
    if ( idynamic == 3 ) then
      call h_interpolator_destroy(vdot_hint)
    end if
  end subroutine conclude_era5

  subroutine conclude_era5h
    implicit none
    call h_interpolator_destroy(cross_hint)
  end subroutine conclude_era5h

end module mod_era5
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

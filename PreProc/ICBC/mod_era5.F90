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

  use mod_intkinds, only : ik2, ik4
  use mod_realkinds, only : rk8, rkx
  use mod_dynparam, only : dattyp, dayspy, globidate1, i_band, i_crm, &
    ibdyfrq, idynamic, inpglob, iy, jx, kz, pthsep, ptop
  use mod_stdio, only : stderr, stdout
  use mod_message, only : die
  use mod_memutil, only : getmem
  use mod_grid, only : dlat, dlon, pa, pjd, pju, pjv, sigmah, tlayer, topogm, &
    topou, topov, ulat, ulon, vlat, vlon, xlat, xlon, z0, za, zetau, zetav
  use mod_date, only : rcm_time_and_date, rcm_time_interval, setcal, &
    split_idate, timeval2date, tochar, tohours, uhrs, yeardayfrac, &
    assignment(=), operator(+), operator(-), operator(/=), operator(<), &
    operator(<=), operator(==), operator(>), operator(>=)
  use mod_constants, only : d_10, dlowval, secph
  use mod_write, only : clt, pd4, pr, ps4, q4, ssr, strd, t4, ts4, u4, v4
  use mod_vertint, only : intgtb, intpsn, intv1, intv2, intv3, intz1, intz3, &
    intzps
  use mod_earth, only : get_window, global_domain
  use mod_humid, only : sph2mxr
  use mod_mksst, only : readsst
  use mod_vectutil, only : crs2dot, top2btm, ucrs2dot, vcrs2dot
  use mod_nchelper, only : checkncerr, nf90_myget_att_text
  use mod_kdinterp, only : h_interpolate_cont, h_interpolate_nn, &
    h_interpolator, h_interpolator_create, h_interpolator_destroy
  use netcdf
  implicit none

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
  real(rkx), pointer, contiguous, dimension(:,:,:) :: uvar, vvar, hvar
  real(rkx), pointer, contiguous, dimension(:,:,:) :: qvar, tvar
  real(rkx), pointer, contiguous, dimension(:,:) :: prvar

  real(rkx), pointer, contiguous, dimension(:) :: glat
  real(rkx), pointer, contiguous, dimension(:) :: grev
  real(rkx), pointer, contiguous, dimension(:) :: glon
  real(rkx), pointer, contiguous, dimension(:) :: plevs
  real(rkx), pointer, contiguous, dimension(:) :: sigmar
  real(rkx) :: pss, pst
  integer(ik2), pointer, contiguous, dimension(:,:,:) :: iwork3
  real(rkx), pointer, contiguous, dimension(:,:,:) :: rwork3
  integer(ik2), pointer, contiguous, dimension(:,:) :: iwork
  real(rkx), pointer, contiguous, dimension(:,:) :: rwork

  type(rcm_time_and_date), pointer, contiguous, dimension(:) :: itimes
  real(rk8), pointer, contiguous, dimension(:) :: xtimes

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
    call getmem(glat,1,jlat,'mod_era5:glat')
    call getmem(glon,1,ilon,'mod_era5:glon')
    call getmem(grev,1,max(jlat,ilon),'mod_era5:grev')
    call getmem(b3,1,jx,1,iy,1,3,'mod_era5:b3')

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
    call getmem(glat,1,jlat,'mod_era5:glat')
    glat = grev(gdomain%jgstart:gdomain%jgstop)
    grev(1:ilon) = glon
    ilon = sum(gdomain%ni)
    call getmem(glon,1,ilon,'mod_era5:glon')
    glon(1:gdomain%ni(1)) = grev(gdomain%igstart(1):gdomain%igstop(1))
    if ( gdomain%ntiles == 2 ) then
      glon(gdomain%ni(1)+1:ilon) = grev(gdomain%igstart(2):gdomain%igstop(2))
    end if

    call h_interpolator_create(cross_hint,glat,glon,xlat,xlon)
    call getmem(b2,1,ilon,1,jlat,1,3,'mod_era5:b2')
    call getmem(prvar,1,ilon,1,jlat,'mod_era5:prvar')
    call getmem(iwork,1,ilon,1,jlat,'mod_era5:iwork')
    call getmem(rwork,1,ilon,1,jlat,'mod_era5:rwork')
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
    call getmem(plevs,1,klev,'mod_era5:plevs')
    call getmem(glat,1,jlat,'mod_era5:glat')
    call getmem(glon,1,ilon,'mod_era5:glon')
    call getmem(grev,1,max(jlat,ilon),'mod_era5:grev')
    call getmem(sigmar,1,klev,'mod_era5:sigmar')
    call getmem(b3,1,jx,1,iy,1,klev*3,'mod_era5:b3')
    if ( idynamic == 3 ) then
      call getmem(d3u,1,jx,1,iy,1,klev*2,'mod_era5:d3u')
      call getmem(d3v,1,jx,1,iy,1,klev*2,'mod_era5:d3v')
      call getmem(h3u,1,jx,1,iy,1,klev,'mod_era5:h3u')
      call getmem(h3v,1,jx,1,iy,1,klev,'mod_era5:h3v')
    else
      call getmem(d3,1,jx,1,iy,1,klev*2,'mod_era5:d3')
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
    call getmem(glat,1,jlat,'mod_era5:glat')
    glat = grev(gdomain%jgstart:gdomain%jgstop)
    grev(1:ilon) = glon
    ilon = sum(gdomain%ni)
    call getmem(glon,1,ilon,'mod_era5:glon')
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

    call getmem(b2,1,ilon,1,jlat,1,klev*3,'mod_era5:b2')
    call getmem(d2,1,ilon,1,jlat,1,klev*2,'mod_era5:d2')
    call getmem(iwork3,1,ilon,1,jlat,1,klev,'mod_era5:iwork3')
    call getmem(rwork3,1,ilon,1,jlat,1,klev,'mod_era5:rwork3')
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
  end subroutine init_era5

  subroutine get_era5h(idate)
    implicit none
    type(rcm_time_and_date), intent(in) :: idate
    type(rcm_time_interval) :: onehour
    ! ERA5 accumulations are stamped at the end of their hour: the fluxes
    ! for the hour starting at idate are in the record of idate + 1 hour.
    onehour = rcm_time_interval(1,uhrs)
    call era5hour(idate+onehour,globidate1+onehour)
    write (stdout,*) 'READ IN fields at DATE:', tochar(idate)
    ! Precipitation by nearest neighbour, a choice: on a grid much finer
    ! than ERA5 (0.25 deg) it shows the ERA5 cells as blocks, but keeps
    ! the rain edges sharp. The continuous interpolation would smooth the
    ! blocks and spread the rain edges into drizzle.
    call h_interpolate_nn(cross_hint,prvar,pr)
    call h_interpolate_cont(cross_hint,b2,b3)
    ! m of water in the hour -> kg m-2 s-1
    ! Clip the round-off negatives of the interpolation
    pr(:,:) = max(pr(:,:) * 1000.0_rkx / secph,0.0_rkx)
    ! ssr is the ERA5 NET surface shortwave: clmsa recovers the downward
    ! flux with the CLM albedo. The ERA5 downward flux (ssrd) would be
    ! better, and should replace it if it becomes available for download.
    ssr(:,:) = max(b3(:,:,1) / secph,0.0_rkx)
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
                  xlat,yeardayfrac(idate),dayspy,jx,iy,klev)
      call intz3(ts4,t3,h3,topogm,0.6_rkx,0.5_rkx,0.85_rkx,jx,iy,klev)
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
      call intz1(u4,u3,zetau,h3u,topou,jx,iy,kz,klev,0.6_rkx,0.2_rkx,0.2_rkx)
!$OMP SECTION
      call intz1(v4,v3,zetav,h3v,topov,jx,iy,kz,klev,0.6_rkx,0.2_rkx,0.2_rkx)
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
    integer(ik4) :: inet, it, kkrec, istatus, ivar
    character(len=64) :: inname
    character(len=256) :: pathaddname
    character(len=4), dimension(4) :: varname
    character(len=4), dimension(4) :: fname
    integer(ik4), dimension(3) :: icount, istart
    integer(ik4) :: year, month, day, hour
    integer(ik4), save :: lastmonth
    type(rcm_time_interval) :: tdif
    ! Each file may be packed (short with scale_factor and add_offset)
    ! or not, independently of the others.
    logical, save, dimension(4) :: has_offset = .false.
    real(rkx), save, dimension(4) :: xoff, xscl
    integer(ik4), save, dimension(4) :: inet5
    integer(ik4), save, dimension(4) :: ivar5

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
              'hourly'//pthsep//inname
        else
          write(inname,'(a,a,i0.4,a,i0.2,a)') &
            trim(fname(kkrec)), '_', year, '_', month,'.nc'
          pathaddname = trim(inpglob)//pthsep//'ERA5'//pthsep// &
              'hourly'//pthsep//inname
        end if
        call open_era5_file(pathaddname,varname(kkrec),inet5(kkrec), &
               ivar5(kkrec),has_offset(kkrec),xscl(kkrec),xoff(kkrec))
        if ( kkrec == 1 ) then
          call read_era5_times(inet5(1),year,month)
        end if
      end do
    end if

    it = record_index(idate,1)

    istart(3) = it
    icount(3) = 1

    do kkrec = 1, 4
      inet = inet5(kkrec)
      ivar = ivar5(kkrec)
      if ( has_offset(kkrec) ) then
        call igetwork(kkrec)
        rwork = real(iwork,rkx)*xscl(kkrec) + xoff(kkrec)
      else
        call rgetwork(kkrec)
      end if
      if ( kkrec == 1 ) then
        prvar = rwork
      else
        b2(:,:,kkrec-1) = rwork
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
    integer(ik4) :: inet, it, kkrec, istatus, ivar
    character(len=64) :: inname
    character(len=256) :: pathaddname
    character(len=1), dimension(5) :: varname
    character(len=4), dimension(5) :: fname
    integer(ik4), save, dimension(5) :: inet5
    integer(ik4), save, dimension(5) :: ivar5
    ! Each file may be packed (short with scale_factor and add_offset)
    ! or not, independently of the others.
    logical, save, dimension(5) :: has_offset = .false.
    real(rkx), save, dimension(5) :: xoff, xscl
    integer(ik4), dimension(4) :: icount, istart
    integer(ik4) :: year, month, day, hour
    integer(ik4), save :: lastmonth
    ! Standard gravity of the ERA5 geopotential
    real(rkx), parameter :: era5_g = 9.80665_rkx

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
        call open_era5_file(pathaddname,varname(kkrec),inet5(kkrec), &
               ivar5(kkrec),has_offset(kkrec),xscl(kkrec),xoff(kkrec))
        if ( kkrec == 1 ) then
          call read_era5_times(inet5(1),year,month)
        end if
      end do
    end if

    it = record_index(idate,ibdyfrq)

    istart(3) = 1
    icount(3) = klev
    istart(4) = it
    icount(4) = 1

    do kkrec = 1, 5
      inet = inet5(kkrec)
      ivar = ivar5(kkrec)
      if ( has_offset(kkrec) ) then
        call igetwork(kkrec)
        rwork3 = real(iwork3,rkx)*xscl(kkrec) + xoff(kkrec)
      else
        call rgetwork(kkrec)
      end if
      select case (kkrec)
        case (1)
          tvar = rwork3
        case (2)
          hvar = rwork3/era5_g
        case (3)
          qvar = max(rwork3,0.0_rkx)
          call sph2mxr(qvar,ilon,jlat,klev)
          qvar = log10(max(qvar,dlowval))
        case (4)
          uvar = rwork3
        case (5)
          vvar = rwork3
      end select
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
  !
  ! Open an ERA5 file, find the variable and its packing: a short
  ! variable has a scale_factor and, optionally, an add_offset.
  !
  subroutine open_era5_file(pathaddname,vname,ncid,ivar,packed,xscl,xoff)
    implicit none
    character(len=*), intent(in) :: pathaddname, vname
    integer(ik4), intent(out) :: ncid, ivar
    logical, intent(out) :: packed
    real(rkx), intent(out) :: xscl, xoff
    integer(ik4) :: istatus
    istatus = nf90_open(pathaddname,nf90_nowrite,ncid)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Error open file '//trim(pathaddname))
    istatus = nf90_inq_varid(ncid,trim(vname),ivar)
    call checkncerr(istatus,__FILE__,__LINE__, &
      'Error find var '//trim(vname)//' in '//trim(pathaddname))
    istatus = nf90_get_att(ncid,ivar,'scale_factor',xscl)
    packed = ( istatus == nf90_noerr )
    if ( packed ) then
      istatus = nf90_get_att(ncid,ivar,'add_offset',xoff)
      if ( istatus /= nf90_noerr ) xoff = 0.0_rkx
      write (stdout,*) ncid, trim(pathaddname), xscl, xoff
    else
      xscl = 1.0_rkx
      xoff = 0.0_rkx
      write (stdout,*) ncid, trim(pathaddname)
    end if
  end subroutine open_era5_file
  !
  ! Read the times of a monthly ERA5 file. The climatology (ERAXX)
  ! has no time axis: its records start at the first of the month.
  !
  subroutine read_era5_times(ncid,year,month)
    implicit none
    integer(ik4), intent(in) :: ncid, year, month
    integer(ik4) :: istatus, timid, it
    character(len=64) :: cunit, ccal
    if ( dattyp == 'ERAXX' ) then
      timlen = 0
      call getmem(itimes,1,1,'mod_era5:itimes')
      itimes(1) = year*1000000 + month*10000+100
      call setcal(itimes(1),'noleap')
      return
    end if
    istatus = nf90_inq_dimid(ncid,'time',timid)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_dimid(ncid,'valid_time',timid)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find dim time')
    end if
    istatus = nf90_inquire_dimension(ncid,timid,len=timlen)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error inquire time')
    istatus = nf90_inq_varid(ncid,'time',timid)
    if ( istatus /= nf90_noerr ) then
      istatus = nf90_inq_varid(ncid,'valid_time',timid)
      if ( istatus /= nf90_noerr ) then
        istatus = nf90_inq_varid(ncid,'date',timid)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error find var time/date')
      end if
    end if
    istatus = nf90_myget_att_text(ncid,timid,'units',cunit)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read time units')
    istatus = nf90_myget_att_text(ncid,timid,'calendar',ccal)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read time calendar')
    call getmem(itimes,1,timlen,'mod_era5:itimes')
    call getmem(xtimes,1,timlen,'mod_era5:xtimes')
    istatus = nf90_get_var(ncid,timid,xtimes)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error read time')
    do it = 1, timlen
      itimes(it) = timeval2date(xtimes(it),cunit,ccal)
    end do
  end subroutine read_era5_times
  !
  ! Record of idate in the open files. The step is the one of the file
  ! times, and the time of the record must be idate. The climatology
  ! (ERAXX) has no times: its step is the nominal one (hours).
  !
  integer(ik4) function record_index(idate,nominal_step) result(it)
    implicit none
    type(rcm_time_and_date), intent(in) :: idate
    integer(ik4), intent(in) :: nominal_step
    real(rk8) :: step
    if ( timlen == 0 ) then
      it = nint(tohours(idate - itimes(1)))/nominal_step + 1
      return
    end if
    if ( timlen > 1 ) then
      step = tohours(itimes(2) - itimes(1))
    else
      step = 1.0_rk8
    end if
    it = nint(tohours(idate - itimes(1))/step) + 1
    if ( it < 1 .or. it > timlen ) then
      write (stderr,*) 'Date ', tochar(idate), ' outside the ERA5 file: ', &
        tochar(itimes(1)), ' to ', tochar(itimes(timlen))
      call die('mod_era5','Missing date in ERA5 file',1)
    end if
    if ( itimes(it) /= idate ) then
      write (stderr,*) 'Date ', tochar(idate), ' requested, record ', it, &
        ' of the ERA5 file is at ', tochar(itimes(it))
      call die('mod_era5','Wrong time step in ERA5 file',1)
    end if
  end function record_index

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

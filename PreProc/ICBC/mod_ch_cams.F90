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

!
! Chemistry (CH) and aerosol (AE) boundary conditions from the CAMS global
! reanalysis (EAC4) on pressure levels, chemtyp = 'CAMSR'. One file per
! variable and month: <inpglob>/CAMS/<year>/<NAME>_<year>_<month>.nc, with
! the variable names of the ADS (aermr01..aermr11, so2, go3, co, ...),
! mass mixing ratios in kg/kg, packed or not.
!
! The fields are interpolated horizontally on the model grid and in the
! vertical, linearly in log(p), to the pressure of the model sigma levels
! p = (ps - ptop) * sigma + ptop, with the surface pressure of the ICBC
! file at the same time. This is the vertical coordinate of all the
! chemistry boundary files: the model interpolates them to its levels at
! run time (non hydrostatic and MOLOCH cores).
!
! The IFS-AER aerosols are converted to the RegCM species (Remy et al.,
! 2019, GMD 12, 4627):
!  - sea salt is given at 80 % relative humidity: the mass is divided by
!    4.3 to obtain the dry mass, and the bin radii by 1.99, giving dry
!    diameter bins 0.03-0.5, 0.5-5 and 5-20 micron;
!  - organic matter is divided by the OM:OC ratio of 1.8 used in IFS-AER
!    to obtain the organic carbon of RegCM;
!  - the bins are redistributed on the RegCM bins (dry diameter: dust
!    0.01-1, 1-2.5, 2.5-5, 5-20 micron; sea salt 0.05-1, 1-10 micron)
!    with the IFS-AER size distributions (dust lognormal, mass median
!    diameter 2.45 micron, sigma 2; sea salt coarse mode, mass median dry
!    diameter 8.5 micron, sigma 2). Mass outside the RegCM bins is lost.
!
module mod_cams

  use mod_intkinds, only : ik4
  use mod_realkinds, only : rk8, rkx
  use mod_dynparam, only : chemtyp, dirglob, domname, globidate1, i_band, &
    inpglob, iy, jx, kz, pthsep, ptop
  use mod_stdio, only : stdout
  use mod_memutil, only : getmem
  use mod_grid, only : sigmah, xlat, xlon
  use mod_date, only : rcm_time_and_date, split_idate, timeval2date, &
    tochar, tochar10, monfirst, lsamemonth, assignment(=), operator(==)
  use mod_wrtoxd, only : aev4, chv4, write_ae_icbc, write_ch_icbc
  use mod_earth, only : get_window, global_domain
  use mod_message, only : die
  use mod_nchelper, only : checkncerr
  use mod_kdinterp, only : h_interpolate_cont, h_interpolator, &
    h_interpolator_create, h_interpolator_destroy
  use mod_ch_param, only : ae_bchb, ae_bchl, ae_dust1, ae_dust2, ae_dust3, &
    ae_dust4, ae_ochb, ae_ochl, ae_so2, ae_so4, ae_sslt1, ae_sslt2, cb_c2h6, &
    cb_ch4, cb_co, cb_h2o2, cb_hcho, cb_hno3, cb_isop, cb_no, cb_no2, cb_o3, &
    cb_pan, cb_par, cb_so2
  use netcdf
  implicit none

  private

  ! Sea salt: wet (80 % RH) to dry mass
  real(rkx), parameter :: sslt_wet2dry = 1.0_rkx/4.3_rkx
  ! Organic matter to organic carbon
  real(rkx), parameter :: om2oc = 1.0_rkx/1.8_rkx
  ! Fractions of the CAMS coarse dust bin (diameter 1.8-40 micron) in the
  ! RegCM bins 2 (1-2.5), 3 (2.5-5) and 4 (5-20 micron)
  real(rkx), parameter :: fdu3_2 = 0.27_rkx
  real(rkx), parameter :: fdu3_3 = 0.50_rkx
  real(rkx), parameter :: fdu3_4 = 0.23_rkx
  ! Fractions of the CAMS sea salt bin 2 (dry diameter 0.5-5 micron) in
  ! the RegCM bins 1 (0.05-1) and 2 (1-10 micron), and of the bin 3 (5-20
  ! micron) in the RegCM bin 2
  real(rkx), parameter :: fss2_1 = 0.10_rkx
  real(rkx), parameter :: fss2_2 = 0.90_rkx
  real(rkx), parameter :: fss3_2 = 0.55_rkx

  integer(ik4), parameter :: ngas = 13
  integer(ik4), parameter :: naer = 12

  ! ADS variable names and file names. The order is the one of the
  ! assignments in cams_to_regcm.
  character(len=8), dimension(ngas), parameter :: gasname = &
    [ character(len=8) :: 'go3', 'co', 'c2h6', 'hcho', 'c3h8', 'no', &
      'no2', 'h2o2', 'c5h8', 'so2', 'hno3', 'ch4', 'pan' ]
  character(len=8), dimension(ngas), parameter :: gasfile = &
    [ character(len=8) :: 'O3', 'CO', 'ETH', 'CH2O', 'C3H8', 'NO', &
      'NO2', 'H2O2', 'ISOP', 'SO2', 'HNO3', 'CH4', 'PAN' ]
  character(len=8), dimension(naer), parameter :: aername = &
    [ character(len=8) :: 'aermr01', 'aermr02', 'aermr03', 'aermr04', &
      'aermr05', 'aermr06', 'aermr07', 'aermr08', 'aermr09', 'aermr10', &
      'aermr11', 'so2' ]
  character(len=8), dimension(naer), parameter :: aerfile = &
    [ character(len=8) :: 'SSLT1', 'SSLT2', 'SSLT3', 'DUST1', 'DUST2', &
      'DUST3', 'OCHL', 'OCHB', 'BCHL', 'BCHB', 'SO4', 'SO2' ]

  integer(ik4) :: nsp, jlat, ilon, klev, timlen, ncicbc, ivarps
  character(len=8), dimension(max(ngas,naer)) :: vname, fname
  integer(ik4), dimension(max(ngas,naer)) :: ncid, ivar
  real(rkx), dimension(max(ngas,naer)) :: xscl, xoff
  logical :: lopen = .false.
  logical :: lpsopen = .false.
  type(rcm_time_and_date) :: filedate, psdate
  type(rcm_time_and_date), pointer, contiguous, dimension(:) :: itimes

  real(rkx), pointer, contiguous, dimension(:,:,:) :: b2, b3
  real(rkx), pointer, contiguous, dimension(:,:,:) :: work
  real(rkx), pointer, contiguous, dimension(:,:) :: ps3
  real(rkx), pointer, contiguous, dimension(:,:,:,:) :: xv
  real(rkx), pointer, contiguous, dimension(:) :: glat, glon, grev
  real(rk8), dimension(:), allocatable :: plog
  ! Pressure levels of the input, top to bottom (Pa)
  real(rkx), pointer, contiguous, dimension(:) :: plevs
  logical :: lflip

  type(global_domain) :: gdomain
  type(h_interpolator) :: cross_hint

  public :: init_cams, get_cams, conclude_cams

  contains

  subroutine init_cams(typ)
    implicit none
    character(len=2), intent(in) :: typ
    integer(ik4) :: istatus, inc, ivarid, idimid
    character(len=256) :: pathaddname
    character(len=64) :: cunit

    select case ( typ )
      case ( 'CH' )
        nsp = ngas
        vname(1:nsp) = gasname
        fname(1:nsp) = gasfile
      case ( 'AE' )
        nsp = naer
        vname(1:nsp) = aername
        fname(1:nsp) = aerfile
      case default
        call die('mod_ch_cams','Unknown type: support only CH and AE')
    end select
    lopen = .false.

    pathaddname = cams_file(globidate1,fname(1))
    istatus = nf90_open(pathaddname,nf90_nowrite,inc)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error open file '//trim(pathaddname))
    call get_dim(inc,'latitude','',jlat)
    call get_dim(inc,'longitude','',ilon)
    call get_dim(inc,'pressure_level','level',klev)
    call getmem(glat,1,jlat,'mod_cams:glat')
    call getmem(glon,1,ilon,'mod_cams:glon')
    call getmem(grev,1,max(jlat,ilon),'mod_cams:grev')
    call getmem(plevs,1,klev,'mod_cams:plevs')
    call get_coord(inc,'latitude','',glat)
    call get_coord(inc,'longitude','',glon)
    call get_coord(inc,'pressure_level','level',plevs)
    ! Pressure levels in hPa (ADS) or Pa
    istatus = nf90_inq_varid(inc,'pressure_level',ivarid)
    if ( istatus /= nf90_noerr ) istatus = nf90_inq_varid(inc,'level',ivarid)
    cunit = 'hPa'
    istatus = nf90_get_att(inc,ivarid,'units',cunit)
    if ( cunit(1:2) /= 'Pa' ) plevs = plevs * 100.0_rkx
    istatus = nf90_close(inc)
    call checkncerr(istatus,__FILE__,__LINE__, &
                    'Error close file '//trim(pathaddname))
    ! Vertical interpolation from top to bottom
    lflip = ( plevs(1) > plevs(klev) )
    if ( lflip ) plevs = plevs(klev:1:-1)
    if ( allocated(plog) ) deallocate(plog)
    allocate(plog(klev))
    plog = log(real(plevs,rk8))
    !
    ! Find window to read
    !
    call get_window(glat,glon,xlat,xlon,i_band,gdomain)
    grev(1:jlat) = glat
    jlat = gdomain%nj
    call getmem(glat,1,jlat,'mod_cams:glat')
    glat = grev(gdomain%jgstart:gdomain%jgstop)
    grev(1:ilon) = glon
    ilon = sum(gdomain%ni)
    call getmem(glon,1,ilon,'mod_cams:glon')
    glon(1:gdomain%ni(1)) = grev(gdomain%igstart(1):gdomain%igstop(1))
    if ( gdomain%ntiles == 2 ) then
      glon(gdomain%ni(1)+1:ilon) = grev(gdomain%igstart(2):gdomain%igstop(2))
    end if
    call h_interpolator_create(cross_hint,glat,glon,xlat,xlon)

    call getmem(b2,1,ilon,1,jlat,1,klev*nsp,'mod_cams:b2')
    call getmem(b3,1,jx,1,iy,1,klev*nsp,'mod_cams:b3')
    call getmem(work,1,ilon,1,jlat,1,klev,'mod_cams:work')
    call getmem(ps3,1,jx,1,iy,'mod_cams:ps3')
    call getmem(xv,1,jx,1,iy,1,kz,1,nsp,'mod_cams:xv')

  contains

    subroutine get_dim(nc,dname,altname,dlen)
      implicit none
      integer(ik4), intent(in) :: nc
      character(len=*), intent(in) :: dname, altname
      integer(ik4), intent(out) :: dlen
      istatus = nf90_inq_dimid(nc,dname,idimid)
      if ( istatus /= nf90_noerr .and. len(altname) > 0 ) then
        istatus = nf90_inq_dimid(nc,altname,idimid)
      end if
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Missing dimension '//dname//' in CAMS file')
      istatus = nf90_inquire_dimension(nc,idimid,len=dlen)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read dimension '//dname//' in CAMS file')
    end subroutine get_dim

    subroutine get_coord(nc,cname,altname,val)
      implicit none
      integer(ik4), intent(in) :: nc
      character(len=*), intent(in) :: cname, altname
      real(rkx), dimension(:), intent(out) :: val
      istatus = nf90_inq_varid(nc,cname,ivarid)
      if ( istatus /= nf90_noerr .and. len(altname) > 0 ) then
        istatus = nf90_inq_varid(nc,altname,ivarid)
      end if
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Missing variable '//cname//' in CAMS file')
      istatus = nf90_get_var(nc,ivarid,val)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read variable '//cname//' in CAMS file')
    end subroutine get_coord

  end subroutine init_cams

  subroutine get_cams(idate,typ)
    implicit none
    type(rcm_time_and_date), intent(in) :: idate
    character(len=2), intent(in) :: typ
    integer(ik4) :: n

    call read_cams(idate)
    write (stdout,*) 'READ IN CAMS '//typ//' fields at DATE:', tochar(idate)
    call h_interpolate_cont(cross_hint,b2,b3)
    call read_icbc_ps(idate)
    do n = 1, nsp
      call vertint(b3(:,:,(n-1)*klev+1:n*klev),xv(:,:,:,n))
    end do
    call cams_to_regcm(typ)
    if ( typ == 'CH' ) call write_ch_icbc(idate)
    if ( typ == 'AE' ) call write_ae_icbc(idate)
  end subroutine get_cams

  subroutine cams_to_regcm(typ)
    implicit none
    character(len=2), intent(in) :: typ
    if ( typ == 'CH' ) then
      chv4 = 0.0_rkx
      chv4(:,:,:,cb_o3) = xv(:,:,:,1)
      chv4(:,:,:,cb_co) = xv(:,:,:,2)
      chv4(:,:,:,cb_c2h6) = xv(:,:,:,3)
      chv4(:,:,:,cb_hcho) = xv(:,:,:,4)
      ! functional groups: three PAR for each propane
      chv4(:,:,:,cb_par) = 3.0_rkx * xv(:,:,:,5)
      chv4(:,:,:,cb_no) = xv(:,:,:,6)
      chv4(:,:,:,cb_no2) = xv(:,:,:,7)
      chv4(:,:,:,cb_h2o2) = xv(:,:,:,8)
      chv4(:,:,:,cb_isop) = xv(:,:,:,9)
      chv4(:,:,:,cb_so2) = xv(:,:,:,10)
      chv4(:,:,:,cb_hno3) = xv(:,:,:,11)
      chv4(:,:,:,cb_ch4) = xv(:,:,:,12)
      chv4(:,:,:,cb_pan) = xv(:,:,:,13)
    else
      ! Sea salt (1-3), dry mass
      aev4(:,:,:,ae_sslt1) = sslt_wet2dry * &
        (xv(:,:,:,1) + fss2_1*xv(:,:,:,2))
      aev4(:,:,:,ae_sslt2) = sslt_wet2dry * &
        (fss2_2*xv(:,:,:,2) + fss3_2*xv(:,:,:,3))
      ! Dust (4-6)
      aev4(:,:,:,ae_dust1) = xv(:,:,:,4)
      aev4(:,:,:,ae_dust2) = xv(:,:,:,5) + fdu3_2*xv(:,:,:,6)
      aev4(:,:,:,ae_dust3) = fdu3_3*xv(:,:,:,6)
      aev4(:,:,:,ae_dust4) = fdu3_4*xv(:,:,:,6)
      ! Organic matter (7 hydrophilic, 8 hydrophobic) to organic carbon
      aev4(:,:,:,ae_ochl) = om2oc * xv(:,:,:,7)
      aev4(:,:,:,ae_ochb) = om2oc * xv(:,:,:,8)
      ! Black carbon (9 hydrophilic, 10 hydrophobic)
      aev4(:,:,:,ae_bchl) = xv(:,:,:,9)
      aev4(:,:,:,ae_bchb) = xv(:,:,:,10)
      ! Sulfate and sulfur dioxide
      aev4(:,:,:,ae_so4) = xv(:,:,:,11)
      aev4(:,:,:,ae_so2) = xv(:,:,:,12)
    end if
  end subroutine cams_to_regcm
  !
  ! Linear interpolation in log(p) from the input pressure levels to the
  ! model sigma levels, constant above the top and below the lowest level.
  !
  subroutine vertint(f3,fz)
    implicit none
    real(rkx), dimension(:,:,:), intent(in) :: f3
    real(rkx), dimension(:,:,:), intent(out) :: fz
    real(rkx), dimension(klev) :: f
    integer(ik4) :: i, j, k, kk
    real(rk8) :: lp, w
    do i = 1, iy
      do j = 1, jx
        if ( lflip ) then
          f = f3(j,i,klev:1:-1)
        else
          f = f3(j,i,:)
        end if
        kk = 2
        do k = 1, kz
          lp = log((real(ps3(j,i),rk8)-ptop*1000.0_rk8)*sigmah(k) + &
                   ptop*1000.0_rk8)
          if ( lp <= plog(1) ) then
            fz(j,i,k) = f(1)
          else if ( lp >= plog(klev) ) then
            fz(j,i,k) = f(klev)
          else
            do while ( plog(kk) < lp )
              kk = kk + 1
            end do
            ! weight of the lower level kk, 1 at p = plevs(kk)
            w = (lp-plog(kk-1))/(plog(kk)-plog(kk-1))
            fz(j,i,k) = max(real(w*f(kk)+(1.0_rk8-w)*f(kk-1),rkx),0.0_rkx)
          end if
        end do
      end do
    end do
  end subroutine vertint

  subroutine read_cams(idate)
    implicit none
    type(rcm_time_and_date), intent(in) :: idate
    integer(ik4) :: n, it, istatus, itile, iti, itf, timid
    integer(ik4), dimension(4) :: istart, icount
    character(len=256) :: pathaddname
    character(len=64) :: cunit, ccal
    real(rk8), dimension(:), allocatable :: xtimes

    if ( .not. lopen .or. .not. lsamemonth(idate,filedate) ) then
      if ( lopen ) then
        do n = 1, nsp
          istatus = nf90_close(ncid(n))
          call checkncerr(istatus,__FILE__,__LINE__,'Error close file')
        end do
      end if
      filedate = idate
      do n = 1, nsp
        pathaddname = cams_file(idate,fname(n))
        istatus = nf90_open(pathaddname,nf90_nowrite,ncid(n))
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error open file '//trim(pathaddname))
        istatus = nf90_inq_varid(ncid(n),vname(n),ivar(n))
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error find var '//trim(vname(n))// &
                        ' in '//trim(pathaddname))
        ! Unpacked data have no scale factor and offset
        istatus = nf90_get_att(ncid(n),ivar(n),'scale_factor',xscl(n))
        if ( istatus /= nf90_noerr ) xscl(n) = 1.0_rkx
        istatus = nf90_get_att(ncid(n),ivar(n),'add_offset',xoff(n))
        if ( istatus /= nf90_noerr ) xoff(n) = 0.0_rkx
        write (stdout,*) 'Opened ', trim(pathaddname)
      end do
      istatus = nf90_inq_dimid(ncid(1),'valid_time',timid)
      if ( istatus /= nf90_noerr ) then
        istatus = nf90_inq_dimid(ncid(1),'time',timid)
      end if
      call checkncerr(istatus,__FILE__,__LINE__,'Error find dim time')
      istatus = nf90_inquire_dimension(ncid(1),timid,len=timlen)
      call checkncerr(istatus,__FILE__,__LINE__,'Error inquire time')
      istatus = nf90_inq_varid(ncid(1),'valid_time',timid)
      if ( istatus /= nf90_noerr ) then
        istatus = nf90_inq_varid(ncid(1),'time',timid)
      end if
      call checkncerr(istatus,__FILE__,__LINE__,'Error find var time')
      istatus = nf90_get_att(ncid(1),timid,'units',cunit)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read time units')
      istatus = nf90_get_att(ncid(1),timid,'calendar',ccal)
      if ( istatus /= nf90_noerr ) ccal = 'gregorian'
      allocate(xtimes(timlen))
      istatus = nf90_get_var(ncid(1),timid,xtimes)
      call checkncerr(istatus,__FILE__,__LINE__,'Error read time')
      call getmem(itimes,1,timlen,'mod_cams:itimes')
      do it = 1, timlen
        itimes(it) = timeval2date(xtimes(it),cunit,ccal)
      end do
      deallocate(xtimes)
      lopen = .true.
    end if
    it = 0
    do n = 1, timlen
      if ( itimes(n) == idate ) then
        it = n
        exit
      end if
    end do
    if ( it == 0 ) then
      call die('mod_ch_cams','Date '//tochar(idate)//' not in CAMS file '// &
               trim(cams_file(idate,fname(1))))
    end if
    istart(3) = 1
    icount(3) = klev
    istart(4) = it
    icount(4) = 1
    do n = 1, nsp
      iti = 1
      do itile = 1, gdomain%ntiles
        istart(1) = gdomain%igstart(itile)
        icount(1) = gdomain%ni(itile)
        istart(2) = gdomain%jgstart
        icount(2) = gdomain%nj
        itf = iti + gdomain%ni(itile) - 1
        istatus = nf90_get_var(ncid(n),ivar(n),work(iti:itf,:,:), &
                               istart,icount)
        call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error read var '//trim(vname(n)))
        iti = iti + gdomain%ni(itile)
      end do
      b2(:,:,(n-1)*klev+1:n*klev) = work*xscl(n) + xoff(n)
    end do
  end subroutine read_cams
  !
  ! Surface pressure (Pa) of the ICBC file of the domain at idate
  !
  subroutine read_icbc_ps(idate)
    implicit none
    type(rcm_time_and_date), intent(in) :: idate
    type(rcm_time_and_date) :: fdate, tdate
    integer(ik4) :: istatus, timid, n, it, nt
    integer(ik4), dimension(3) :: istart, icount
    character(len=256) :: icbcname
    character(len=64) :: cunit, ccal
    real(rk8), dimension(:), allocatable :: xtimes

    if ( lsamemonth(idate,globidate1) ) then
      fdate = globidate1
    else
      fdate = monfirst(idate)
    end if
    if ( .not. lpsopen .or. .not. (fdate == psdate) ) then
      if ( lpsopen ) then
        istatus = nf90_close(ncicbc)
        call checkncerr(istatus,__FILE__,__LINE__,'Error close ICBC file')
      end if
      write (icbcname,'(a,a,a,a,a,a)') trim(dirglob), pthsep, &
              trim(domname), '_ICBC.', trim(tochar10(fdate)), '.nc'
      istatus = nf90_open(icbcname,nf90_nowrite,ncicbc)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error open ICBC file '//trim(icbcname))
      istatus = nf90_inq_varid(ncicbc,'ps',ivarps)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error find var ps in ICBC file '//trim(icbcname))
      psdate = fdate
      lpsopen = .true.
    end if
    istatus = nf90_inq_dimid(ncicbc,'time',timid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find ICBC time dim')
    istatus = nf90_inquire_dimension(ncicbc,timid,len=nt)
    call checkncerr(istatus,__FILE__,__LINE__,'Error inquire ICBC time')
    istatus = nf90_inq_varid(ncicbc,'time',timid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find ICBC time var')
    istatus = nf90_get_att(ncicbc,timid,'units',cunit)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read ICBC time units')
    istatus = nf90_get_att(ncicbc,timid,'calendar',ccal)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read ICBC calendar')
    allocate(xtimes(nt))
    istatus = nf90_get_var(ncicbc,timid,xtimes)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read ICBC time')
    it = 0
    do n = 1, nt
      tdate = timeval2date(xtimes(n),cunit,ccal)
      if ( tdate == idate ) then
        it = n
        exit
      end if
    end do
    deallocate(xtimes)
    if ( it == 0 ) then
      call die('mod_ch_cams','Date '//tochar(idate)//' not in ICBC file')
    end if
    istart = [ 1, 1, it ]
    icount = [ jx, iy, 1 ]
    istatus = nf90_get_var(ncicbc,ivarps,ps3,istart,icount)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read ICBC ps')
    ! hPa to Pa
    ps3 = ps3 * 100.0_rkx
  end subroutine read_icbc_ps

  character(len=256) function cams_file(idate,name)
    implicit none
    type(rcm_time_and_date), intent(in) :: idate
    character(len=*), intent(in) :: name
    integer(ik4) :: year, month, day, hour
    character(len=64) :: inname
    call split_idate(idate,year,month,day,hour)
    write(inname,'(i4,a,a,a,i0.4,a,i0.2,a)') &
      year, pthsep, trim(name), '_', year, '_', month, '.nc'
    cams_file = trim(inpglob)//pthsep//chemtyp(1:4)//pthsep//inname
  end function cams_file

  subroutine conclude_cams
    implicit none
    integer(ik4) :: n, istatus
    if ( lopen ) then
      do n = 1, nsp
        istatus = nf90_close(ncid(n))
      end do
      lopen = .false.
    end if
    if ( lpsopen ) then
      istatus = nf90_close(ncicbc)
      lpsopen = .false.
    end if
    call h_interpolator_destroy(cross_hint)
  end subroutine conclude_cams

end module mod_cams
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

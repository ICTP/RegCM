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

module mod_cams

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_stdio
  use mod_memutil
  use mod_grid
  use mod_date
  use mod_constants
  use mod_wrtoxd
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
  use mod_ch_param
  use netcdf

  private

  integer(ik4) :: jlat , ilon , klev , timlen
  real(rkx) , pointer , dimension(:,:,:) :: b3 , b3a
  real(rkx) , pointer , dimension(:,:,:) :: b2 , b2a

  real(rkx) , pointer , dimension(:,:,:) :: hz , hza , o3 , co , &
    dust1 , dust2 , dust3 , eth , xch2o , h2o2 , bchl , ochl , bchb , ochb
  real(rkx) , pointer , dimension(:,:,:) :: oh , isop , hno3 , no , &
    pan , c3h8 , sslt1 , sslt2 , sslt3 , so4 , so2 , no2 , ch4 , aso2
  real(rkx) , pointer , dimension(:,:,:) :: hzvar , hzavar , o3var , &
    covar , dust1var , dust2var , dust3var , ethvar , xch2ovar, &
    h2o2var , bchlvar , rochlvar , bchbvar , ochbvar , ochlvar
  real(rkx) , pointer , dimension(:,:,:) :: ohvar , isopvar , hno3var , &
    novar , panvar , c3h8var , sslt1var , sslt2var , sslt3var , so4var , &
    so2var , no2var , ch4var , aso2var
  real(rkx) , pointer , dimension(:,:) :: prvar , topou , topov
  real(rkx) , pointer , dimension(:) :: glat
  real(rkx) , pointer , dimension(:) :: grev
  real(rkx) , pointer , dimension(:) :: glon
  real(rkx) , pointer , dimension(:) :: plevs
  real(rkx) , pointer , dimension(:) :: sigmar
  real(rkx) :: pss , pst
  integer(2) , pointer , dimension(:,:,:) :: work
  integer(2) , pointer , dimension(:,:) :: iwork

  integer(ik4) , dimension(20) :: inet5 !care the first dimension must be
                                        !bigger than the number of species
                                        !for gas or aer
  integer(ik4) , dimension(20) :: ivar5
  real(rkx) , dimension(20) :: xoff , xscl
  type(rcm_time_and_date) , pointer , dimension(:) :: itimes
  integer(ik4) , pointer , dimension(:) :: xtimes

  type(global_domain) :: gdomain
  type(h_interpolator) :: cross_hint , udot_hint , vdot_hint

  public :: init_cams , get_cams, conclude_cams

  contains

  subroutine init_cams(typ)

    implicit none
    character(len=2) :: typ
    integer(ik4) :: k
    integer(ik4) :: year , month , day , hour
    character(len=256) :: pathaddname
    integer(ik4) :: istatus , ncid , ivarid , idimid
    character(len=64) :: inname

    if ( idynamic < 3 ) &
      call die('mod_ch_cams','CAMSR only works for idynamic=3')

    call split_idate(globidate1,year,month,day,hour)
    write(inname,'(i4,a,a,i0.4,a,i0.2,a)') &
      year, pthsep, 'geop_', year, '_', month,'.nc'
    pathaddname = trim(inpglob)//pthsep//chemtyp(1:4)//pthsep//inname
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
      call checkncerr(istatus,__FILE__,__LINE__, &
            'Missing level/levelist dimension in file '//trim(pathaddname))
    end if
    istatus = nf90_inquire_dimension(ncid,idimid,len=klev)
    call checkncerr(istatus,__FILE__,__LINE__, &
          'Error reading levelist dimelen in file '//trim(pathaddname))
    !
    ! Allocate working space
    !
    call getmem1d(plevs,1,klev,'mod_cams:plevs')
    call getmem1d(glat,1,jlat,'mod_cams:glat')
    call getmem1d(glon,1,ilon,'mod_cams:glon')
    call getmem1d(grev,1,max(jlat,ilon),'mod_cams:grev')
    call getmem1d(sigmar,1,klev,'mod_cams:sigmar')

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
    call getmem1d(glat,1,jlat,'mod_cams:glat')
    glat = grev(gdomain%jgstart:gdomain%jgstop)
    grev(1:ilon) = glon
    ilon = sum(gdomain%ni)
    call getmem1d(glon,1,ilon,'mod_cams:glon')
    glon(1:gdomain%ni(1)) = grev(gdomain%igstart(1):gdomain%igstop(1))
    if ( gdomain%ntiles == 2 ) then
      glon(gdomain%ni(1)+1:ilon) = grev(gdomain%igstart(2):gdomain%igstop(2))
    end if

    call h_interpolator_create(cross_hint,glat,glon,xlat,xlon)

    if ( typ == 'CH' ) then
      call getmem3d(b3,1,jx,1,iy,1,klev*15,'mod_cams:b3')
      call getmem3d(b2,1,ilon,1,jlat,1,klev*15,'mod_cams:b2')
      call getmem3d(work,1,ilon,1,jlat,1,klev,'mod_cams:work')
      !
      ! Set up pointers
      !
      hz    => b3(:,:,1:klev)
      o3    => b3(:,:,klev+1:2*klev)
      co    => b3(:,:,2*klev+1:3*klev)
      eth   => b3(:,:,3*klev+1:4*klev)
      xch2o => b3(:,:,4*klev+1:5*klev)
      h2o2  => b3(:,:,5*klev+1:6*klev)
      oh    => b3(:,:,6*klev+1:7*klev)
      isop  => b3(:,:,7*klev+1:8*klev)
      hno3  => b3(:,:,8*klev+1:9*klev)
      no    => b3(:,:,9*klev+1:10*klev)
      pan   => b3(:,:,10*klev+1:11*klev)
      c3h8  => b3(:,:,11*klev+1:12*klev)
      so2   => b3(:,:,12*klev+1:13*klev)
      no2   => b3(:,:,13*klev+1:14*klev)
      ch4   => b3(:,:,14*klev+1:15*klev)

      hzvar    => b2(:,:,1:klev)
      o3var    => b2(:,:,klev+1:2*klev)
      covar    => b2(:,:,2*klev+1:3*klev)
      ethvar   => b2(:,:,3*klev+1:4*klev)
      xch2ovar => b2(:,:,4*klev+1:5*klev)
      h2o2var  => b2(:,:,5*klev+1:6*klev)
      ohvar    => b2(:,:,6*klev+1:7*klev)
      isopvar  => b2(:,:,7*klev+1:8*klev)
      hno3var  => b2(:,:,8*klev+1:9*klev)
      novar    => b2(:,:,9*klev+1:10*klev)
      panvar   => b2(:,:,10*klev+1:11*klev)
      c3h8var  => b2(:,:,11*klev+1:12*klev)
      so2var   => b2(:,:,12*klev+1:13*klev)
      no2var   => b2(:,:,13*klev+1:14*klev)
      ch4var   => b2(:,:,14*klev+1:15*klev)

    else if ( typ =='AE' ) then
      call getmem3d(b3a,1,jx,1,iy,1,klev*15,'mod_cams:b3a')
      call getmem3d(b2a,1,ilon,1,jlat,1,klev*15,'mod_cams:b2a')
      call getmem3d(work,1,ilon,1,jlat,1,klev,'mod_cams:work')
      !
      ! Set up pointers
      !
      hza   => b3a(:,:,1:klev)
      dust1 => b3a(:,:,klev+1:2*klev)
      dust2 => b3a(:,:,2*klev+1:3*klev)
      dust3 => b3a(:,:,3*klev+1:4*klev)
      bchl  => b3a(:,:,4*klev+1:5*klev)
      ochl  => b3a(:,:,5*klev+1:6*klev)
      bchb  => b3a(:,:,6*klev+1:7*klev)
      ochb  => b3a(:,:,7*klev+1:8*klev)
      sslt1 => b3a(:,:,8*klev+1:9*klev)
      sslt2 => b3a(:,:,9*klev+1:10*klev)
      sslt3 => b3a(:,:,10*klev+1:11*klev)
      so4   => b3a(:,:,11*klev+1:12*klev)
      aso2  => b3a(:,:,12*klev+1:13*klev)

      hzavar   => b2a(:,:,1:klev)
      dust1var => b2a(:,:,1*klev+1:2*klev)
      dust2var => b2a(:,:,2*klev+1:3*klev)
      dust3var => b2a(:,:,3*klev+1:4*klev)
      bchlvar  => b2a(:,:,4*klev+1:5*klev)
      ochlvar  => b2a(:,:,5*klev+1:6*klev)
      bchbvar  => b2a(:,:,6*klev+1:7*klev)
      ochbvar  => b2a(:,:,7*klev+1:8*klev)
      sslt1var => b2a(:,:,8*klev+1:9*klev)
      sslt2var => b2a(:,:,9*klev+1:10*klev)
      sslt3var => b2a(:,:,10*klev+1:11*klev)
      so4var   => b2a(:,:,11*klev+1:12*klev)
      aso2var  => b2a(:,:,12*klev+1:13*klev)
    end if
  end subroutine init_cams

  subroutine get_cams(idate,typ)
    implicit none
    character(len=2), intent(in) :: typ
    type(rcm_time_and_date) , intent(in) :: idate
    !
    ! Read data at idate
    !
    if ( typ == 'CH' ) then

      call cams6hour(dattyp,idate,globidate1,'CH')
      write (stdout,*) 'READ IN CAMS fields at DATE:' , tochar(idate)
      !
      ! Horizontal interpolation of both the scalar and vector fields
      !
      call h_interpolate_cont(cross_hint,b2,b3)
      !
      ! Invert vertical order, set BOTTOM -> TOP
!$OMP SECTIONS
!$OMP SECTION
      call top2btm(hz)
!$OMP SECTION
      call top2btm(o3)
!$OMP SECTION
      call top2btm(co)
!$OMP SECTION
      call top2btm(eth )
!$OMP SECTION
      call top2btm(xch2o)
!$OMP SECTION
      call top2btm(h2o2 )
!$OMP SECTION
      call top2btm(oh)
!$OMP SECTION
      call top2btm(isop)
!$OMP SECTION
      call top2btm (hno3)
!$OMP SECTION
      call top2btm(no)
!$OMP SECTION
      call top2btm(pan)
!$OMP SECTION
      call top2btm(c3h8 )
!$OMP SECTION
      call top2btm(so2)
!$OMP SECTION
      call top2btm(no2)
!$OMP SECTION
      call top2btm(ch4)
!$OMP END SECTIONS
    else if ( typ == 'AE' ) then
      call cams6hour(dattyp,idate,globidate1,'AE')
      write (stdout,*) 'READ IN CAMS AER fields at DATE:' , tochar(idate)
      !
      ! Horizontal interpolation of both the scalar and vector fields
      !
      call h_interpolate_cont(cross_hint,b2a,b3a)
      !
!$OMP SECTIONS
!$OMP SECTION
      call top2btm(hza)
!$OMP SECTION
      call top2btm(dust1)
!$OMP SECTION
      call top2btm(dust2)
!$OMP SECTION
      call top2btm(dust3)
!$OMP SECTION
      call top2btm(bchl)
!$OMP SECTION
      call top2btm(ochl)
!$OMP SECTION
      call top2btm(bchb)
!$OMP SECTION
      call top2btm(ochb)
!$OMP SECTION
      call top2btm(sslt1)
!$OMP SECTION
      call top2btm(sslt2)
!$OMP SECTION
      call top2btm(sslt3)
!$OMP SECTION
      call top2btm(so4)
!$OMP SECTION
      call top2btm(aso2)
!$OMP END SECTIONS
    end if
    ! New calculation of p* on rcm topography.
    !
    ! Interpolate chemicals.
    !
    if ( typ == 'CH' ) then
!$OMP SECTIONS
!$OMP SECTION
      call intz1(chv4(:,:,:,cb_O3),o3,z0,hz,topogm, &
         jx,iy,kz,klev,0.7_rkx,0.4_rkx,0.7_rkx)
!$OMP SECTION
      call intz1(chv4(:,:,:,cb_CO),co,z0,hz,topogm, &
         jx,iy,kz,klev,0.7_rkx,0.4_rkx,0.7_rkx)
!$OMP SECTION
      call intz1(chv4(:,:,:,cb_C2H6),eth,z0,hz,topogm, &
          jx,iy,kz,klev,0.7_rkx,0.4_rkx,0.7_rkx)
!$OMP SECTION
      call intz1(chv4(:,:,:,cb_HCHO),xch2o,z0,hz,topogm, &
          jx,iy,kz,klev,0.7_rkx,0.4_rkx,0.7_rkx)
!$OMP SECTION
      call intz1(chv4(:,:,:,cb_H2O2),h2o2,z0,hz,topogm, &
          jx,iy,kz,klev,0.7_rkx,0.4_rkx,0.7_rkx)
!$OMP SECTION
      call intz1(chv4(:,:,:,cb_ISOP),isop,z0,hz,topogm, &
          jx,iy,kz,klev,0.7_rkx,0.4_rkx,0.7_rkx)
!$OMP SECTION
      call intz1(chv4(:,:,:,cb_HNO3),hno3,z0,hz,topogm, &
          jx,iy,kz,klev,0.7_rkx,0.4_rkx,0.7_rkx)
!$OMP SECTION
      call intz1(chv4(:,:,:,cb_NO),no,z0,hz,topogm, &
          jx,iy,kz,klev,0.7_rkx,0.4_rkx,0.7_rkx)
!$OMP SECTION
      call intz1(chv4(:,:,:,cb_PAN),pan,z0,hz,topogm, &
          jx,iy,kz,klev,0.7_rkx,0.4_rkx,0.7_rkx)
!$OMP SECTION
      call intz1(chv4(:,:,:,cb_PAR),c3h8,z0,hz,topogm, &
          jx,iy,kz,klev,0.7_rkx,0.4_rkx,0.7_rkx)
!$OMP SECTION
      call intz1(chv4(:,:,:,cb_SO2),so2,z0,hz,topogm, &
          jx,iy,kz,klev,0.7_rkx,0.4_rkx,0.7_rkx)
!$OMP SECTION
      call intz1(chv4(:,:,:,cb_NO2),no2,z0,hz,topogm, &
          jx,iy,kz,klev,0.7_rkx,0.4_rkx,0.7_rkx)
!$OMP SECTION
      call intz1(chv4(:,:,:,cb_CH4),ch4,z0,hz,topogm, &
          jx,iy,kz,klev,0.7_rkx,0.4_rkx,0.7_rkx)
!$OMP END SECTIONS
      ! treat special case of functional groups
      chv4(:,:,:,cb_PAR) = 3._rkx* chv4(:,:,:,cb_PAR)
      ! consider putting OH

    else if ( typ == 'AE' ) then
!$OMP SECTIONS
      call intz1(aev4(:,:,:,ae_dust1),dust1,z0,hza,topogm, &
          jx,iy,kz,klev,0.7_rkx,0.4_rkx,0.7_rkx)
!$OMP SECTION
      call intz1(aev4(:,:,:,ae_dust2),dust2,z0,hza,topogm, &
          jx,iy,kz,klev,0.7_rkx,0.4_rkx,0.7_rkx)
!$OMP SECTION
      call intz1(aev4(:,:,:,ae_dust3),dust3,z0,hza,topogm, &
          jx,iy,kz,klev,0.7_rkx,0.4_rkx,0.7_rkx)
!$OMP SECTION
      call intz1(aev4(:,:,:,ae_bchl),bchl,z0,hza,topogm, &
          jx,iy,kz,klev,0.7_rkx,0.4_rkx,0.7_rkx)
!$OMP SECTION
      call intz1(aev4(:,:,:,ae_ochl),ochl,z0,hza,topogm, &
          jx,iy,kz,klev,0.7_rkx,0.4_rkx,0.7_rkx)
!$OMP SECTION
      call intz1(aev4(:,:,:,ae_bchb),bchb,z0,hza,topogm, &
          jx,iy,kz,klev,0.7_rkx,0.4_rkx,0.7_rkx)
!$OMP SECTION
      call intz1(aev4(:,:,:,ae_ochb),ochb,z0,hza,topogm, &
          jx,iy,kz,klev,0.7_rkx,0.4_rkx,0.7_rkx)
!$OMP SECTION
      call intz1(aev4(:,:,:,ae_sslt1),sslt1,z0,hza,topogm, &
          jx,iy,kz,klev,0.7_rkx,0.4_rkx,0.7_rkx)
!$OMP SECTION
      call intz1(aev4(:,:,:,ae_sslt2),sslt2,z0,hza,topogm, &
          jx,iy,kz,klev,0.7_rkx,0.4_rkx,0.7_rkx)
!$OMP SECTION
      call intz1(aev4(:,:,:,ae_so4),so4,z0,hza,topogm, &
          jx,iy,kz,klev,0.7_rkx,0.4_rkx,0.7_rkx)
!$OMP SECTION
      call intz1(aev4(:,:,:,ae_so2),aso2,z0,hza,topogm, &
          jx,iy,kz,klev,0.7_rkx,0.4_rkx,0.7_rkx)
!$OMP END SECTIONS

      ! Perform mass redistribution from cams bins to regcm bins for dust and
      ! seasalt. The weights are roughly estimated from bin size and dust
      ! volume distribution between 0 and 20um

      aev4(:,:,:,ae_dust1) = aev4(:,:,:,ae_dust1) + aev4(:,:,:,ae_dust2)
      aev4(:,:,:,ae_dust2) = 0.4_rkx * aev4(:,:,:,ae_dust3)
      aev4(:,:,:,ae_dust4) = 0.2_rkx * aev4(:,:,:,ae_dust3) ! 4 before 3 !!
      aev4(:,:,:,ae_dust3) = 0.4_rkx * aev4(:,:,:,ae_dust3)

      aev4(:,:,:,ae_sslt1) = aev4(:,:,:,ae_sslt1) + &
                             0.2_rkx * aev4(:,:,:,ae_sslt2)
      aev4(:,:,:,ae_sslt2) = 0.8_rkx * aev4(:,:,:,ae_sslt2)
      ! did not use sslt3 consider using a third ssbin ...

    end if

    !write out chv4
    if ( typ=='CH' ) call  write_ch_icbc(idate)
    if ( typ=='AE' ) call  write_ae_icbc(idate)

  end subroutine get_cams

  subroutine cams6hour(dattyp,idate,idate0,typ)
    implicit none
    character(len=5) , intent(in) :: dattyp
    character(len=2), intent(in) ::typ
    type(rcm_time_and_date) , intent(in) :: idate , idate0
    integer(ik4) :: i , inet , it , j , kkrec , istatus , ivar
    integer(ik4) :: timid
    character(len=64) :: inname
    character(len=256) :: pathaddname
    character(len=7) , dimension(14) :: gasname , gfname
    character(len=7) , dimension(13) :: aername , afname
    ! make sure that dim is large enough
    character(len=7) , dimension(15) ::varname
    character(len=7) , dimension(15) ::fname

    character(len=64) :: cunit , ccal
    real(rkx) :: xadd , xscale
    integer(ik4) , dimension(4) :: icount , istart
    integer(ik4) :: year , month , day , hour , nsp
    integer(ik4) , save :: lastmonth
    type(rcm_time_interval) :: tdif

    data gasname /'z','go3','co','c2h6','hcho','c3h8','no', &
                  'no2','h2o2','c5h8' ,'so2','hno3','ch4','pan' /
    data gfname  /'geop','O3','CO','ETH','CH2O','C3H8','NO', &
                  'NO2','H2O2','ISOP','SO2','HNO3','CH4','PAN'/
    data aername /'z','aermr04','aermr05','aermr06','aermr09',&
                  'aermr07','aermr10','aermr08','aermr01', &
                  'aermr02' ,'aermr03','aermr11','so2' /
    data afname  /'geop','DUST1','DUST2','DUST3','BCHL','OCHL', &
                  'BCHB','OCHB','SSLT1','SSLT2','SSLT3','SO4','SO2'/

    !
    ! This is the latitude, longitude dimension of the grid to be read.
    ! This corresponds to the lat and lon dimension variables in the
    ! netCDF file.
    ! The data are packed into short integers (INTEGER*2).  The array
    ! work will be used to hold the packed integers.  The array 'x'
    ! will contain the unpacked data.
    !
    if (typ == 'CH') then
      nsp = 14
      varname(1:nsp) = gasname(1:nsp)
      fname(1:nsp) = gfname(1:nsp)
    else if (typ == 'AE') then
      nsp = 13
      varname(1:nsp) = aername(1:nsp)
      fname(1:nsp) = afname(1:nsp)
    end if

    call split_idate(idate,year,month,day,hour)

    if ( idate == idate0 .or. month /= lastmonth ) then
      lastmonth = month
      if ( idate /= idate0 ) then
        do kkrec = 1 , nsp ! boulce sur nbre variables gaz
          istatus = nf90_close(inet5(kkrec))
          call checkncerr(istatus,__FILE__,__LINE__, &
             'Error close file')
        end do
      end if
      do kkrec = 1 , nsp ! boulce sur nbre variables gaz
        !verifie le path pour fichier CAMS !!
        write(inname,'(i4,a,a,a,i0.4,a,i0.2,a)') &
        year, pthsep, trim(fname(kkrec)), '_', year, '_', month,'.nc'
        pathaddname = trim(inpglob)//pthsep//chemtyp(1:4)//pthsep//inname
        istatus = nf90_open(pathaddname,nf90_nowrite,inet5(kkrec))

        call checkncerr(istatus,__FILE__,__LINE__, &
          'Error open file '//trim(pathaddname))
        istatus = nf90_inq_varid(inet5(kkrec),varname(kkrec), &
                                 ivar5(kkrec))
        call checkncerr(istatus,__FILE__,__LINE__, &
          'Error find var '//varname(kkrec))
        istatus = nf90_get_att(inet5(kkrec),ivar5(kkrec), &
                 'scale_factor',xscl(kkrec))
        call checkncerr(istatus,__FILE__,__LINE__, &
            'Error find att scale_factor')
        istatus = nf90_get_att(inet5(kkrec),ivar5(kkrec),  &
                   'add_offset',xoff(kkrec))
        call checkncerr(istatus,__FILE__,__LINE__, &
            'Error find att add_offset')
        write (stdout,*) inet5(kkrec) , trim(pathaddname) ,   &
                         xscl(kkrec) , xoff(kkrec)
        if ( kkrec == 1 ) then
          istatus = nf90_inq_dimid(inet5(1),'time',timid)
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error find dim time')
          istatus = nf90_inquire_dimension(inet5(1),timid,len=timlen)
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error inquire time')
          istatus = nf90_inq_varid(inet5(1),'time',timid)
          if ( istatus /= nf90_noerr ) then
            istatus = nf90_inq_varid(inet5(1),'date',timid)
            call checkncerr(istatus,__FILE__,__LINE__, &
                        'Error find var time/date')
          end if
          istatus = nf90_get_att(inet5(1),timid,'units',cunit)
          call checkncerr(istatus,__FILE__,__LINE__, &
                              'Error read time units')
          istatus = nf90_get_att(inet5(1),timid,'calendar',ccal)
          call checkncerr(istatus,__FILE__,__LINE__, &
                              'Error read time units')
          call getmem1d(itimes,1,timlen,'mod_cams:itimes')
          call getmem1d(xtimes,1,timlen,'mod_cams:xtimes')
          istatus = nf90_get_var(inet5(1),timid,xtimes)
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error read time')
          do it = 1 , timlen
            itimes(it) = timeval2date(real(xtimes(it),rkx),cunit,ccal)
          end do
        end if
      end do
    end if
    tdif = idate - itimes(1)
    it = nint(tohours(tdif))/6 + 1

    istart(3) = 1
    icount(3) = klev
    istart(4) = it
    icount(4) = 1

    if ( typ == 'CH' ) then
      do kkrec = 1 , nsp ! loop on number variables
        inet = inet5(kkrec)
        ivar = ivar5(kkrec)
        xscale = xscl(kkrec)
        xadd = xoff(kkrec)
        call getwork(kkrec)
        !care the order must match varname data - could be more elegant,
        if ( kkrec == 1 ) hzvar(1:ilon,1:jlat,:) = &
          real(real(work(1:ilon,1:jlat,:),rkx)*xscale+xadd,rkx)/9.80616_rk4
        if ( kkrec == 2 ) o3var(1:ilon,1:jlat,:) = &
          real(real(work(1:ilon,1:jlat,:),rkx)*xscale+xadd,rkx)
        if ( kkrec == 3 ) covar(1:ilon,1:jlat,:) = &
          real(real(work(1:ilon,1:jlat,:),rkx)*xscale+xadd,rkx)
        if ( kkrec == 4 ) ethvar(1:ilon,1:jlat,:) = &
          real(real(work(1:ilon,1:jlat,:),rkx)*xscale+xadd,rkx)
        if ( kkrec == 5 ) xch2ovar(1:ilon,1:jlat,:) = &
          real(real(work(1:ilon,1:jlat,:),rkx)*xscale+xadd,rkx)
        if ( kkrec == 6 ) c3h8var(1:ilon,1:jlat,:) = &
          real(real(work(1:ilon,1:jlat,:),rkx)*xscale+xadd,rkx)
        if ( kkrec == 7 ) novar(1:ilon,1:jlat,:) = &
          real(real(work(1:ilon,1:jlat,:),rkx)*xscale+xadd,rkx)
        if ( kkrec == 8 ) no2var(1:ilon,1:jlat,:) = &
          real(real(work(1:ilon,1:jlat,:),rkx)*xscale+xadd,rkx)
        if ( kkrec == 9 ) h2o2var(1:ilon,1:jlat,:) = &
          real(real(work(1:ilon,1:jlat,:),rkx)*xscale+xadd,rkx)
        if ( kkrec == 10 ) isopvar(1:ilon,1:jlat,:) = &
          real(real(work(1:ilon,1:jlat,:),rkx)*xscale+xadd,rkx)
        if ( kkrec == 11 ) so2var(1:ilon,1:jlat,:) = &
          real(real(work(1:ilon,1:jlat,:),rkx)*xscale+xadd,rkx)
        if ( kkrec == 12 ) hno3var(1:ilon,1:jlat,:) = &
          real(real(work(1:ilon,1:jlat,:),rkx)*xscale+xadd,rkx)
        if ( kkrec == 13 ) ch4var(1:ilon,1:jlat,:) = &
          real(real(work(1:ilon,1:jlat,:),rkx)*xscale+xadd,rkx)
        if ( kkrec == 14 ) panvar(1:ilon,1:jlat,:) = &
          real(real(work(1:ilon,1:jlat,:),rkx)*xscale+xadd,rkx)
      end do
    else if ( typ == 'AE' ) then
      do kkrec = 1 , nsp ! loop on number variables
        inet = inet5(kkrec)
        ivar = ivar5(kkrec)
        xscale = xscl(kkrec)
        xadd = xoff(kkrec)
        call getwork(kkrec)
        !care the order must match varname data - could be more elegant,
        if ( kkrec == 1 ) hzavar(1:ilon,1:jlat,:) = &
          real(real(work(1:ilon,1:jlat,:),rkx)*xscale+xadd,rkx) / 9.80616_rk4
        if ( kkrec == 2 ) dust1var(1:ilon,1:jlat,:) = &
          real(real(work(1:ilon,1:jlat,:),rkx)*xscale+xadd,rkx)
        if ( kkrec == 3 ) dust2var(1:ilon,1:jlat,:) = &
          real(real(work(1:ilon,1:jlat,:),rkx)*xscale+xadd,rkx)
        if ( kkrec == 4 ) dust3var(1:ilon,1:jlat,:) = &
          real(real(work(1:ilon,1:jlat,:),rkx)*xscale+xadd,rkx)
        if ( kkrec == 5 ) bchlvar(1:ilon,1:jlat,:) = &
          real(real(work(1:ilon,1:jlat,:),rkx)*xscale+xadd,rkx)
        if ( kkrec == 6 ) ochlvar(1:ilon,1:jlat,:) = &
          real(real(work(1:ilon,1:jlat,:),rkx)*xscale+xadd,rkx)
        if ( kkrec == 7 ) bchbvar(1:ilon,1:jlat,:) = &
          real(real(work(1:ilon,1:jlat,:),rkx)*xscale+xadd,rkx)
        if ( kkrec == 8 ) ochbvar(1:ilon,1:jlat,:) = &
          real(real(work(1:ilon,1:jlat,:),rkx)*xscale+xadd,rkx)
        if ( kkrec == 9 ) sslt1var(1:ilon,1:jlat,:) = &
          real(real(work(1:ilon,1:jlat,:),rkx)*xscale+xadd,rkx)
        if ( kkrec == 10 )sslt2var(1:ilon,1:jlat,:) = &
          real(real(work(1:ilon,1:jlat,:),rkx)*xscale+xadd,rkx)
        if ( kkrec == 11 )sslt3var(1:ilon,1:jlat,:) = &
          real(real(work(1:ilon,1:jlat,:),rkx)*xscale+xadd,rkx)
        if ( kkrec == 12 )so4var(1:ilon,1:jlat,:) = &
          real(real(work(1:ilon,1:jlat,:),rkx)*xscale+xadd,rkx)
        if ( kkrec == 13 )aso2var(1:ilon,1:jlat,:) = &
          real(real(work(1:ilon,1:jlat,:),rkx)*xscale+xadd,rkx)
      end do
    end if

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
  end subroutine cams6hour

  subroutine conclude_cams
    implicit none
    call h_interpolator_destroy(cross_hint)
    call h_interpolator_destroy(udot_hint)
    if ( idynamic == 3 ) then
      call h_interpolator_destroy(vdot_hint)
    end if
  end subroutine conclude_cams

end module mod_cams
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

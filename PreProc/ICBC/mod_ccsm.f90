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

module mod_ccsm

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx!
!
! This is a package of subroutines to read CCSM T85 and T42 L26 data in
! NETCDF format and to prepare Initial and boundary conditions for RegCM3.
! Both Global and Window of CCSM data are acceptable.
! Written By Moetasim Ashfaq Dec-2005 @ PURDUE.EDU
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx!
!
!  SUBROUTINE CAM85
!  Read unpacked CCSM NETCDF T85 L26 (six hourly) data and save into data
!    arrays. Each varaible is read in seperate monthly data files.
!
!  SUBROUTINE CAM42
!  Read unpacked CCSM NETCDF T42 L26 (six hourly) data and save into data
!     arrays. Each varaible is read in seperate yearly data files.
!
!  SUBROUTINE HEADER_CAM85 & HEADER_CAM42
!  Define pressure levels, Ak and Bk coeffcients and global grid dimensions
!     In CCSM, the vertical coordinate is a hybrid sigma-pressure system
!     The pressure is defined as:
!     P = Ak*Po+Bk*PS(i,j)
!     All 3D fields required for ICBC are at mid-points, so
!     Ak refers to hyam, the hybrid A coefficient at midpoints, and
!     Bk refers to hybm, the hybrid B coefficient at midpoints
!     Po = 1000mb
!
!  SUBROUTINE GET_CAM85
!  Main subroutine to read data arrays from CAM85 and prepare ICBCs at RCM Grid
!
!  SUBROUTINE GET_CAM42
!  Main subroutine to read data arrays from CAM42 and prepare ICBCs at RCM Grid
!
!  SUBROUTINE INITDATE3
!  Initialize CCSM 365 days calendar (No leap years)
!
!  SUBROUTINE CAMCLNDR
!  Subroutine for SST preparation with No leap year
!
!  SUBROUTINE HANDLE_ERR
!  Handle error for NETCDF calls
!
!  SUBROUTINES CALLED FROM ICBC.f
!     1) INTLIN
!     2) INTLOG
!     3) HUMIF1FV
!     3) BILINX2CR
!     4) BILINX2DT
!     5) TOP2BTM
!     6) UVROT4
!     7) INTGTB
!     8) P1P2
!     8) INTPSN
!     9) INTV3
!     10) MKSST
!     10) HUMID2FV
!     10) INTV1
!     11) WRITEF
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx!
!  DATA PREPARATION
!  Dataset required to use this code can be preapred using NCO utilitiles
!     such as NCKS, NCRCAT etc.
!  Prepare:
!     Monthly data files for CAM85
!     Yearly data files for CAM42
!  For example:
!     To extract global data of CAM42 for Specific Humidity
!     ncks -v Q input.nc ccsm.shum.nyear.nc    ,and
!     to extract a subset (window) of CAM85 data for Specific Humidity
!     ncks -d lat,min,max -d lon,min,max -v Q input.nc ccsm.shumJAN.nyear.nc
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx!
!  NAMING CONVENTION (Global Data Files) CAM85
!  (MONTH) = JAN/FEB/MAR/APR/MAY/JUN/JUL/AUG/SEP/OCT/NOV/DEC
!
!  ccsm.air(MONTH).nyear.nc     for 'T'     (Temperature)
!  ccsm.hgt(MONTH).nyear.nc     for 'Z3'    (Geopotential Height)
!  ccsm.shum(MONTH).nyear.nc    for 'Q'     (Specific Humidity)
!  ccsm.uwnd(MONTH).nyear.nc    for 'U'     (Zonal Wind)
!  ccsm.vwnd(MONTH).nyear.nc    for 'V'     (Meridonial Wind)
!  ccsm.pres(MONTH).nyear.nc    for 'PS'    (Surface Pressure)
!
!  PATH /DATA/CAM85/NYEAR/
!  ccsm_ht.nc      for 'PHIS'  (Surface Geopotential-static field)
!
!  PATH /DATA/CAM85/
!
!  NAMING CONVENTION (Global Data Files) CAM42
!  ccsm.air.nyear        for 'T'   (Temperature)
!  ccsm.hgt.nyear.nc     for 'Z3'  (Geopotential Height)
!  ccsm.shum.nyear.nc    for 'Q'   (Specific Humidity)
!  ccsm.uwnd.nyear.nc    for 'U'   (Zonal Wind)
!  ccsm.vwnd.nyear.nc    for 'V'   (Meridonial Wind)
!  ccsm.pres.nyear.nc    for 'PS'  (Surface Pressure)
!
!  PATH /DATA/CAM42/NYEAR/
!  ccsm_ht.nc      for 'PHIS'  (Surface Geopotential-static field)
!
!  PATH /DATA/CAM42/
!
!
!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx!

  use m_realkinds
  use m_stdio
  use m_die
  use m_mall
  use m_zeit
  use mod_dynparam
  use mod_constants
  use mod_grid
  use mod_write
  use mod_interp
  use mod_vertint
  use mod_hgt
  use mod_humid
  use mod_mksst
  use mod_uvrot
  use mod_vectutil

  implicit none

  private

  integer :: nlon , nlat , klev

  integer , parameter :: npl = 18

  ! Whole space
  real(sp) , allocatable , target , dimension(:,:,:) :: b2
  real(sp) , allocatable , target , dimension(:,:,:) :: d2
  real(sp) , allocatable , target , dimension(:,:,:) :: b3
  real(sp) , allocatable , target , dimension(:,:,:) :: d3

  real(sp) , pointer :: u3(:,:,:) , v3(:,:,:)
  real(sp) , pointer :: h3(:,:,:) , q3(:,:,:) , t3(:,:,:)
  real(sp) , pointer :: up(:,:,:) , vp(:,:,:)
  real(sp) , pointer :: hp(:,:,:) , qp(:,:,:) , tp(:,:,:)

  ! Input space
  real(sp) :: p0
  real(sp) , allocatable , dimension(:,:) :: psvar , zsvar
  real(sp) , allocatable , dimension(:) :: ak , bk
  real(sp) , allocatable , dimension(:) :: glat
  real(sp) , allocatable , dimension(:) :: glon
  real(sp) , allocatable , dimension(:,:,:) :: hvar , qvar , tvar , &
                                               uvar , vvar , pp3d
  integer :: timlen
  type(rcm_time_and_date) , allocatable , dimension(:) :: itimes
  real(sp) , dimension(npl) :: pplev , sigmar

  ! Shared by netcdf I/O routines
  type(rcm_time_and_date) :: ilastdate
  integer , dimension(4) :: icount , istart
  integer , dimension(6) :: inet6 , ivar6

  public :: get_ccsm , headccsm , footerccsm

  contains
!
  subroutine headccsm
    use netcdf
!
    implicit none
!
    integer :: istatus , ivar1 , inet1 , ilat , ilon , ihyam , ihybm , ip0 , k
    integer :: lonid , latid , ilevid
    character(256) :: pathaddname
    real(8) :: dp0
!
    call zeit_ci('headccsm')
    pathaddname = trim(inpglob)//'/CCSM/ccsm_ht.nc'
    istatus = nf90_open(pathaddname,nf90_nowrite,inet1)
    if ( istatus /= nf90_noerr ) call handle_err(istatus)

    istatus = nf90_inq_dimid(inet1,'lon',lonid)
    if ( istatus /= nf90_noerr ) call handle_err(istatus)
    istatus = nf90_inquire_dimension(inet1,lonid,len=nlon)
    if ( istatus /= nf90_noerr ) call handle_err(istatus)
    istatus = nf90_inq_dimid(inet1,'lat',latid)
    if ( istatus /= nf90_noerr ) call handle_err(istatus)
    istatus = nf90_inquire_dimension(inet1,latid,len=nlat)
    if ( istatus /= nf90_noerr ) call handle_err(istatus)
    istatus = nf90_inq_dimid(inet1,'lev',ilevid)
    if ( istatus /= nf90_noerr ) call handle_err(istatus)
    istatus = nf90_inquire_dimension(inet1,ilevid,len=klev)
    if ( istatus /= nf90_noerr ) call handle_err(istatus)

    istatus = nf90_inq_varid(inet1,'lat',ilat)
    if ( istatus /= nf90_noerr ) call handle_err(istatus)
    istatus = nf90_inq_varid(inet1,'lon',ilon)
    if ( istatus /= nf90_noerr ) call handle_err(istatus)
    istatus = nf90_inq_varid(inet1,'hyam',ihyam)
    if ( istatus /= nf90_noerr ) call handle_err(istatus)
    istatus = nf90_inq_varid(inet1,'hybm',ihybm)
    if ( istatus /= nf90_noerr ) call handle_err(istatus)
    istatus = nf90_inq_varid(inet1,'PHIS',ivar1)
    if ( istatus /= nf90_noerr ) call handle_err(istatus)
    istatus = nf90_inq_varid(inet1,'P0',ip0)
    if ( istatus/=nf90_noerr ) call handle_err(istatus)

    ! Input layer and pressure interpolated values

    allocate(glat(nlat),stat=istatus)
    if (istatus /= 0) call die('mod_ccsm','Allocation error on glat',1)
    call mall_mci(glat,'mod_ccsm')
    allocate(glon(nlon),stat=istatus)
    if (istatus /= 0) call die('mod_ccsm','Allocation error on glon',1)
    call mall_mci(glon,'mod_ccsm')
    allocate(zsvar(nlon,nlat),stat=istatus)
    if (istatus /= 0) call die('mod_ccsm','Allocation error on zsvar',1)
    call mall_mci(zsvar,'mod_ccsm')
    allocate(psvar(nlon,nlat),stat=istatus)
    if (istatus /= 0) call die('mod_ccsm','Allocation error on psvar',1)
    call mall_mci(psvar,'mod_ccsm')
    allocate(qvar(nlon,nlat,klev),stat=istatus)
    if (istatus /= 0) call die('mod_ccsm','Allocation error on qsvar',1)
    call mall_mci(qvar,'mod_ccsm')
    allocate(tvar(nlon,nlat,klev),stat=istatus)
    if (istatus /= 0) call die('mod_ccsm','Allocation error on tvar',1)
    call mall_mci(tvar,'mod_ccsm')
    allocate(hvar(nlon,nlat,klev),stat=istatus)
    if (istatus /= 0) call die('mod_ccsm','Allocation error on hvar',1)
    call mall_mci(hvar,'mod_ccsm')
    allocate(uvar(nlon,nlat,klev),stat=istatus)
    if (istatus /= 0) call die('mod_ccsm','Allocation error on uvar',1)
    call mall_mci(uvar,'mod_ccsm')
    allocate(vvar(nlon,nlat,klev),stat=istatus)
    if (istatus /= 0) call die('mod_ccsm','Allocation error on vvar',1)
    call mall_mci(vvar,'mod_ccsm')
    allocate(pp3d(nlon,nlat,klev),stat=istatus)
    if (istatus /= 0) call die('mod_ccsm','Allocation error on pp3d',1)
    call mall_mci(pp3d,'mod_ccsm')
 
    allocate(ak(klev),stat=istatus)
    if (istatus /= 0) call die('mod_ccsm','Allocation error on ak',1)
    call mall_mci(ak,'mod_ccsm')
    allocate(bk(klev),stat=istatus)
    if (istatus /= 0) call die('mod_ccsm','Allocation error on bk',1)
    call mall_mci(bk,'mod_ccsm')

    istatus = nf90_get_var(inet1,ilat,glat)
    if ( istatus /= nf90_noerr ) call handle_err(istatus)
    istatus = nf90_get_var(inet1,ilon,glon)
    if ( istatus /= nf90_noerr ) call handle_err(istatus)
    istatus = nf90_get_var(inet1,ihyam,ak)
    if ( istatus /= nf90_noerr ) call handle_err(istatus)
    istatus = nf90_get_var(inet1,ihybm,bk)
    if ( istatus /= nf90_noerr ) call handle_err(istatus)
    istatus = nf90_get_var(inet1,ip0,dp0)
    if ( istatus/=nf90_noerr ) call handle_err(istatus)
    p0 = real(dp0)

    icount(1) = nlon
    icount(2) = nlat
    icount(3) = 1
    istart(1) = 1
    istart(2) = 1
    istart(3) = 1
    istatus = nf90_get_var(inet1,ivar1,zsvar,istart(1:3),icount(1:3))
    if ( istatus /= nf90_noerr ) call handle_err(istatus)
    zsvar = zsvar/real(egrav)
    where (zsvar < 0.0) zsvar = 0.0

    write (stdout,*) 'Read in Static fields from ',trim(pathaddname),' .'

    istatus = nf90_close(inet1)
    if ( istatus /= nf90_noerr ) call handle_err(istatus)

    pplev(1) = 30.
    pplev(2) = 50.
    pplev(3) = 70.
    pplev(4) = 100.
    pplev(5) = 150.
    pplev(6) = 200.
    pplev(7) = 250.
    pplev(8) = 300.
    pplev(9) = 350.
    pplev(10) = 420.
    pplev(11) = 500.
    pplev(12) = 600.
    pplev(13) = 700.
    pplev(14) = 780.
    pplev(15) = 850.
    pplev(16) = 920.
    pplev(17) = 960.
    pplev(18) = 1000.
!
    do k = 1 , npl
      sigmar(k) = pplev(k)*0.001
    end do
 
    allocate(b3(jx,iy,npl*3),stat=istatus)
    if (istatus /= 0) call die('mod_ccsm','Allocation error on b3',1)
    call mall_mci(b3,'mod_ccsm')
    allocate(d3(jx,iy,npl*2),stat=istatus)
    if (istatus /= 0) call die('mod_ccsm','Allocation error on d3',1)
    call mall_mci(d3,'mod_ccsm')
    allocate(b2(nlon,nlat,npl*3),stat=istatus)
    if (istatus /= 0) call die('mod_ccsm','Allocation error on b2',1)
    call mall_mci(b2,'mod_ccsm')
    allocate(d2(nlon,nlat,npl*2),stat=istatus)
    if (istatus /= 0) call die('mod_ccsm','Allocation error on d2',1)
    call mall_mci(d2,'mod_ccsm')

!       Set up pointers
 
    up => d2(:,:,1:npl)
    vp => d2(:,:,npl+1:2*npl)
    tp => b2(:,:,1:npl)
    hp => b2(:,:,npl+1:2*npl)
    qp => b2(:,:,2*npl+1:3*npl)

    u3 => d3(:,:,1:npl)
    v3 => d3(:,:,npl+1:2*npl)
    t3 => b3(:,:,1:npl)
    h3 => b3(:,:,npl+1:2*npl)
    q3 => b3(:,:,2*npl+1:3*npl)

    call zeit_co('headccsm')
  end subroutine headccsm
!
!-----------------------------------------------------------------------
! 
  subroutine get_ccsm(idate)

    use netcdf

    implicit none
!
    type(rcm_time_and_date) , intent(in) :: idate
!

    call zeit_ci('get_ccsm')
    call readccsm(idate)

    write (stdout,*) 'Read in fields at Date: ', idate%tostring()
 
    call height(hp,hvar,tvar,psvar,pp3d,zsvar,nlon,nlat,klev,pplev,npl)
 
    call humid1fv(tvar,qvar,pp3d,nlon,nlat,klev)

    call intlin(up,uvar,psvar,pp3d,nlon,nlat,klev,pplev,npl)
    call intlin(vp,vvar,psvar,pp3d,nlon,nlat,klev,pplev,npl)
    call intlog(tp,tvar,psvar,pp3d,nlon,nlat,klev,pplev,npl)
    call intlin(qp,qvar,psvar,pp3d,nlon,nlat,klev,pplev,npl)
 
    call bilinx2(b3,b2,xlon,xlat,glon,glat,nlon,nlat,jx,iy,npl*3)
    call bilinx2(d3,d2,dlon,dlat,glon,glat,nlon,nlat,jx,iy,npl*2)

    call uvrot4(u3,v3,dlon,dlat,clon,clat,grdfac,jx,iy,npl,plon,plat,iproj)
 
    call top2btm(t3,jx,iy,npl)
    call top2btm(q3,jx,iy,npl)
    call top2btm(h3,jx,iy,npl)
    call top2btm(u3,jx,iy,npl)
    call top2btm(v3,jx,iy,npl)
 
    call intgtb(pa,za,tlayer,topogm,t3,h3,sigmar,jx,iy,npl)
    call intpsn(ps4,topogm,pa,za,tlayer,ptop,jx,iy)
 
    if(i_band.eq.1) then
       call p1p2_band(b3pd,ps4,jx,iy)
    else
       call p1p2(b3pd,ps4,jx,iy)
    endif
 
    call intv3(ts4,t3,ps4,sigmar,ptop,jx,iy,npl)
    call readsst(ts4,idate)
    call intv1(u4,u3,b3pd,sigma2,sigmar,ptop,jx,iy,kz,npl)
    call intv1(v4,v3,b3pd,sigma2,sigmar,ptop,jx,iy,kz,npl)
    call intv2(t4,t3,ps4,sigma2,sigmar,ptop,jx,iy,kz,npl)
 
    call intv1(q4,q3,ps4,sigma2,sigmar,ptop,jx,iy,kz,npl)
    call humid2(t4,q4,ps4,ptop,sigma2,jx,iy,kz)
 
    call hydrost(h4,t4,topogm,ps4,ptop,sigmaf,sigma2,dsigma,jx,iy,kz)
 
    call zeit_co('get_ccsm')
  end subroutine get_ccsm
!
!-----------------------------------------------------------------------
! 
  subroutine readccsm(idate)
!
    use netcdf
!
    implicit none
!
    type(rcm_time_and_date) , intent(in) :: idate
!
    integer :: istatus
    integer :: i , it , j , k , kkrec , timid
    integer :: inet , ivar
    character(25) :: inname
    character(256) :: pathaddname
    character(2) , dimension(6) :: varname
    real(dp) , allocatable , dimension(:) :: xtimes
    character(3) , dimension(12) :: mname
    character(64) :: cunit , ccal
    logical :: lfound , lfirst
!
    data mname   /'JAN','FEB','MAR','APR','MAY','JUN', &
                  'JUL','AUG','SEP','OCT','NOV','DEC'/
    data varname /'T' , 'Z3' , 'Q' , 'U' , 'V' , 'PS'/
    data lfirst  /.true./
!
    call zeit_ci('readccsm')
!
    if ( lfirst .or. .not. lsamemonth(idate,ilastdate) ) then
      do kkrec = 1 , 6
        if ( kkrec == 1 ) then
          write (inname,99001) idate%year, 'air', mname(idate%month), idate%year
        else if ( kkrec == 2 ) then
          write (inname,99001) idate%year, 'hgt', mname(idate%month), idate%year
        else if ( kkrec == 3 ) then
          write (inname,99002) idate%year, 'shum', mname(idate%month), idate%year
        else if ( kkrec == 4 ) then
          write (inname,99002) idate%year, 'uwnd', mname(idate%month), idate%year
        else if ( kkrec == 5 ) then
          write (inname,99002) idate%year, 'vwnd', mname(idate%month), idate%year
        else if ( kkrec == 6 ) then
          write (inname,99002) idate%year, 'pres', mname(idate%month), idate%year
        end if
 
        pathaddname = trim(inpglob)//'/CCSM/'//inname
 
        istatus = nf90_open(pathaddname,nf90_nowrite,inet6(kkrec))
        if ( istatus /= nf90_noerr ) call handle_err(istatus)
        istatus = nf90_inq_varid(inet6(kkrec),varname(kkrec), ivar6(kkrec))
        if ( istatus /= nf90_noerr ) call handle_err(istatus)
        write (stdout,*) inet6(kkrec), trim(pathaddname), ' : ', varname(kkrec)
        if ( kkrec == 1 ) then
          istatus = nf90_inq_dimid(inet6(kkrec),'time',timid)
          if ( istatus /= nf90_noerr ) call handle_err(istatus)
          istatus = nf90_inquire_dimension(inet6(kkrec),timid, len=timlen)
          if ( istatus /= nf90_noerr ) call handle_err(istatus)
          istatus = nf90_inq_varid(inet6(kkrec),'time',timid)
          if ( istatus /= nf90_noerr ) call handle_err(istatus)
          istatus = nf90_get_att(inet6(kkrec),timid,'units',cunit)
          if ( istatus /= nf90_noerr ) call handle_err(istatus)
          cunit = '-'//trim(cunit)//' GMT-'
          istatus = nf90_get_att(inet6(kkrec),timid,'calendar',ccal)
          if ( istatus /= nf90_noerr ) call handle_err(istatus)

          if (allocated(itimes)) then
            deallocate(itimes)
          end if
          allocate(xtimes(timlen),stat=istatus)
          if (istatus /= 0) call die('mod_ccsm','Allocation error on xtimes',1)
          call mall_mci(xtimes,'mod_ccsm')
          allocate(itimes(timlen),stat=istatus)
          if (istatus /= 0) call die('mod_ccsm','Allocation error on itimes',1)
          istatus = nf90_get_var(inet6(kkrec),timid,xtimes)
          if ( istatus /= nf90_noerr ) call handle_err(istatus)
          do it = 1 , timlen
            itimes(it) = timeval2date(xtimes(it)*24.0,cunit,ccal)
          end do
          deallocate(xtimes)
          call mall_mco(xtimes,'mod_ccsm')
        end if
      end do
      if (lfirst) lfirst = .false.
    end if

    do kkrec = 1 , 6

      lfound = .false.
      do it = 1 , timlen
        if (itimes(it) == idate) then
          lfound = .true.
          exit
        end if
      end do
 
      if ( .not. lfound ) then
        write (stderr,*) idate%tostring(), ' not found in ', trim(pathaddname)
        write (stderr,*) 'Extremes are : ', itimes(1)%tostring(), &
                         '-', itimes(timlen)%tostring()
        call die('readccsm')
      end if

      inet = inet6(kkrec)
      ivar = ivar6(kkrec)

      if ( kkrec == 6 ) then
        icount(1) = nlon
        icount(2) = nlat
        icount(3) = 1
        istart(1) = 1
        istart(2) = 1
        istart(3) = it
        istatus = nf90_get_var(inet,ivar,psvar,istart(1:3),icount(1:3))
        if ( istatus /= nf90_noerr ) call handle_err(istatus)
      else
        icount(1) = nlon
        icount(2) = nlat
        icount(3) = klev
        icount(4) = 1
        istart(1) = 1
        istart(2) = 1
        istart(3) = 1
        istart(4) = it
        if ( kkrec == 1 ) then
          istatus = nf90_get_var(inet,ivar,tvar,istart,icount)
          if ( istatus /= nf90_noerr ) call handle_err(istatus)
        else if ( kkrec == 2 ) then
          istatus = nf90_get_var(inet,ivar,hvar,istart,icount)
          if ( istatus /= nf90_noerr ) call handle_err(istatus)
        else if ( kkrec == 3 ) then
          istatus = nf90_get_var(inet,ivar,qvar,istart,icount)
          if ( istatus /= nf90_noerr ) call handle_err(istatus)
        else if ( kkrec == 4 ) then
          istatus = nf90_get_var(inet,ivar,uvar,istart,icount)
          if ( istatus /= nf90_noerr ) call handle_err(istatus)
        else if ( kkrec == 5 ) then
          istatus = nf90_get_var(inet,ivar,vvar,istart,icount)
          if ( istatus /= nf90_noerr ) call handle_err(istatus)
        end if
      end if
    end do

    do k = 1 , klev
      do j = 1 , nlat
        do i = 1 , nlon
          pp3d(i,j,k) = (ak(k)*p0 + bk(k)*psvar(i,j))*0.01
        end do
      end do
    end do
    do j = 1 , nlat
      do i = 1 , nlon
        psvar(i,j) = psvar(i,j)*0.01
      end do
    end do
 
    ilastdate = idate
    call zeit_co('readccsm')

99001   format (i4,'/','ccsm.',a3,a3,'.',i4,'.nc')
99002   format (i4,'/','ccsm.',a4,a3,'.',i4,'.nc')

  end subroutine readccsm
!
!-----------------------------------------------------------------------
!
  subroutine footerccsm
    implicit none
    call mall_mco(glat,'mod_ccsm')
    deallocate(glat)
    call mall_mco(glon,'mod_ccsm')
    deallocate(glon)
    call mall_mco(zsvar,'mod_ccsm')
    deallocate(zsvar)
    call mall_mco(psvar,'mod_ccsm')
    deallocate(psvar)
    call mall_mco(tvar,'mod_ccsm')
    deallocate(tvar)
    call mall_mco(hvar,'mod_ccsm')
    deallocate(hvar)
    call mall_mco(qvar,'mod_ccsm')
    deallocate(qvar)
    call mall_mco(uvar,'mod_ccsm')
    deallocate(uvar)
    call mall_mco(vvar,'mod_ccsm')
    deallocate(vvar)
    call mall_mco(pp3d,'mod_ccsm')
    deallocate(pp3d)
    call mall_mco(ak,'mod_ccsm')
    deallocate(ak)
    call mall_mco(bk,'mod_ccsm')
    deallocate(bk)
    call mall_mco(b3,'mod_ccsm')
    deallocate(b3)
    call mall_mco(d3,'mod_ccsm')
    deallocate(d3)
    call mall_mco(b2,'mod_ccsm')
    deallocate(b2)
    call mall_mco(d2,'mod_ccsm')
    deallocate(d2)
    deallocate(itimes)
  end subroutine footerccsm
!
!-----------------------------------------------------------------------
!
!     Error Handler for NETCDF Calls
!
  subroutine handle_err(istatus)
    use netcdf
    implicit none
!
    integer :: istatus
    intent (in) :: istatus
!
    write(stderr,*) nf90_strerror(istatus)
    call die('handle_err','Netcdf Error',1)

  end subroutine handle_err
!
end module mod_ccsm

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
  use mod_message

  implicit none

  private

  integer :: nlon , nlat , klev

  integer , parameter :: npl = 18

  ! Whole space
  real(sp) , pointer , dimension(:,:,:) :: b2
  real(sp) , pointer , dimension(:,:,:) :: d2
  real(sp) , pointer , dimension(:,:,:) :: b3
  real(sp) , pointer , dimension(:,:,:) :: d3

  real(sp) , pointer :: u3(:,:,:) , v3(:,:,:)
  real(sp) , pointer :: h3(:,:,:) , q3(:,:,:) , t3(:,:,:)
  real(sp) , pointer :: up(:,:,:) , vp(:,:,:)
  real(sp) , pointer :: hp(:,:,:) , qp(:,:,:) , tp(:,:,:)

  ! Input space
  real(sp) :: p0
  real(sp) , pointer , dimension(:,:) :: psvar , zsvar
  real(sp) , pointer , dimension(:) :: ak , bk
  real(sp) , pointer , dimension(:) :: glat
  real(sp) , pointer , dimension(:) :: glon
  real(sp) , pointer , dimension(:,:,:) :: hvar , qvar , tvar , &
                                               uvar , vvar , pp3d
  integer :: timlen
  type(rcm_time_and_date) , pointer , dimension(:) :: itimes
  real(dp) , pointer , dimension(:) :: xtimes
  real(sp) , dimension(npl) :: pplev , sigmar

  ! Shared by netcdf I/O routines
  type(rcm_time_and_date) :: ilastdate
  integer , dimension(4) :: icount , istart
  integer , dimension(6) :: inet6 , ivar6

  public :: get_ccsm , headccsm

  character(256) :: pathaddname

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
    call checkncerr(istatus,__FILE__,__LINE__,'Error open file '//trim(pathaddname))

    istatus = nf90_inq_dimid(inet1,'lon',lonid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find dim lon')
    istatus = nf90_inquire_dimension(inet1,lonid,len=nlon)
    call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dim lon')
    istatus = nf90_inq_dimid(inet1,'lat',latid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find dim lat')
    istatus = nf90_inquire_dimension(inet1,latid,len=nlat)
    call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dim lat')
    istatus = nf90_inq_dimid(inet1,'lev',ilevid)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find dim lev')
    istatus = nf90_inquire_dimension(inet1,ilevid,len=klev)
    call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dim lev')

    istatus = nf90_inq_varid(inet1,'lat',ilat)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find var lat')
    istatus = nf90_inq_varid(inet1,'lon',ilon)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find var lon')
    istatus = nf90_inq_varid(inet1,'hyam',ihyam)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find var hyam')
    istatus = nf90_inq_varid(inet1,'hybm',ihybm)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find var hybm')
    istatus = nf90_inq_varid(inet1,'PHIS',ivar1)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find var PHIS')
    istatus = nf90_inq_varid(inet1,'P0',ip0)
    call checkncerr(istatus,__FILE__,__LINE__,'Error find var P0')

    ! Input layer and pressure interpolated values

    call getmem1d(glat,1,nlat,'mod_ccsm:glat')
    call getmem1d(glon,1,nlon,'mod_ccsm:glon')
    call getmem2d(zsvar,1,nlon,1,nlat,'mod_ccsm:zsvar')
    call getmem2d(psvar,1,nlon,1,nlat,'mod_ccsm:psvar')
    call getmem3d(qvar,1,nlon,1,nlat,1,klev,'mod_ccsm:qvar')
    call getmem3d(tvar,1,nlon,1,nlat,1,klev,'mod_ccsm:tvar')
    call getmem3d(hvar,1,nlon,1,nlat,1,klev,'mod_ccsm:hvar')
    call getmem3d(uvar,1,nlon,1,nlat,1,klev,'mod_ccsm:uvar')
    call getmem3d(vvar,1,nlon,1,nlat,1,klev,'mod_ccsm:vvar')
    call getmem3d(pp3d,1,nlon,1,nlat,1,klev,'mod_ccsm:pp3d')
    call getmem1d(ak,1,klev,'mod_ccsm:ak')
    call getmem1d(bk,1,klev,'mod_ccsm:bk')
 
    istatus = nf90_get_var(inet1,ilat,glat)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var lat')
    istatus = nf90_get_var(inet1,ilon,glon)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var lon')
    istatus = nf90_get_var(inet1,ihyam,ak)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var hyam')
    istatus = nf90_get_var(inet1,ihybm,bk)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var hybm')
    istatus = nf90_get_var(inet1,ip0,dp0)
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var P0')
    p0 = real(dp0)

    icount(1) = nlon
    icount(2) = nlat
    icount(3) = 1
    istart(1) = 1
    istart(2) = 1
    istart(3) = 1
    istatus = nf90_get_var(inet1,ivar1,zsvar,istart(1:3),icount(1:3))
    call checkncerr(istatus,__FILE__,__LINE__,'Error read var PHIS')
    zsvar = zsvar/real(egrav)
    where (zsvar < 0.0) zsvar = 0.0

    write (stdout,*) 'Read in Static fields from ',trim(pathaddname),' .'

    istatus = nf90_close(inet1)
    call checkncerr(istatus,__FILE__,__LINE__,'Error close file')

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
 
    call getmem3d(b3,1,jx,1,iy,1,npl*3,'mod_ccsm:b3')
    call getmem3d(d3,1,jx,1,iy,1,npl*2,'mod_ccsm:d3')
    call getmem3d(b2,1,nlon,1,nlat,1,npl*3,'mod_ccsm:b2')
    call getmem3d(d2,1,nlon,1,nlat,1,npl*2,'mod_ccsm:d2')

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
    character(2) , dimension(6) :: varname
    character(3) , dimension(12) :: mname
    character(64) :: cunit , ccal
    type(rcm_time_interval) :: tdif
    logical :: lfirst
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
        call checkncerr(istatus,__FILE__,__LINE__,'Error open file '//trim(pathaddname))
        istatus = nf90_inq_varid(inet6(kkrec),varname(kkrec), ivar6(kkrec))
        call checkncerr(istatus,__FILE__,__LINE__,'Error find var '//varname(kkrec))
        write (stdout,*) inet6(kkrec), trim(pathaddname), ' : ', varname(kkrec)
        if ( kkrec == 1 ) then
          istatus = nf90_inq_dimid(inet6(kkrec),'time',timid)
          call checkncerr(istatus,__FILE__,__LINE__,'Error find dim time')
          istatus = nf90_inquire_dimension(inet6(kkrec),timid, len=timlen)
          call checkncerr(istatus,__FILE__,__LINE__,'Error inquire dim time')
          istatus = nf90_inq_varid(inet6(kkrec),'time',timid)
          call checkncerr(istatus,__FILE__,__LINE__,'Error find var time')
          istatus = nf90_get_att(inet6(kkrec),timid,'units',cunit)
          call checkncerr(istatus,__FILE__,__LINE__,'Error find att units')
          istatus = nf90_get_att(inet6(kkrec),timid,'calendar',ccal)
          call checkncerr(istatus,__FILE__,__LINE__,'Error find att calendar')

          call getmem1d(xtimes,1,timlen,'mod_ccsm:xtimes')
          call getmem1d(itimes,1,timlen,'mod_ccsm:itimes')
          istatus = nf90_get_var(inet6(kkrec),timid,xtimes)
          call checkncerr(istatus,__FILE__,__LINE__,'Error read var time')
          do it = 1 , timlen
            itimes(it) = timeval2date(xtimes(it),cunit,ccal)
          end do
        end if
      end do
      if (lfirst) lfirst = .false.
    end if

    tdif = idate - itimes(1)
    it = idnint(tdif%hours())/6 + 1

    do kkrec = 1 , 6

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
        call checkncerr(istatus,__FILE__,__LINE__,'Error read var ps')
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
          call checkncerr(istatus,__FILE__,__LINE__,'Error read var t')
        else if ( kkrec == 2 ) then
          istatus = nf90_get_var(inet,ivar,hvar,istart,icount)
          call checkncerr(istatus,__FILE__,__LINE__,'Error read var h')
        else if ( kkrec == 3 ) then
          istatus = nf90_get_var(inet,ivar,qvar,istart,icount)
          call checkncerr(istatus,__FILE__,__LINE__,'Error read var q')
        else if ( kkrec == 4 ) then
          istatus = nf90_get_var(inet,ivar,uvar,istart,icount)
          call checkncerr(istatus,__FILE__,__LINE__,'Error read var u')
        else if ( kkrec == 5 ) then
          istatus = nf90_get_var(inet,ivar,vvar,istart,icount)
          call checkncerr(istatus,__FILE__,__LINE__,'Error read var v')
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
end module mod_ccsm

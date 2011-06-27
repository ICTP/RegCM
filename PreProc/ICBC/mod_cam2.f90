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

module mod_cam2

!xxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxxx!
!
! This is a package of subroutines to read CAM2 data in
! NETCDF format and to prepare Initial and boundary conditions for RegCM.
! Both Global and Window of CAM2 data are acceptable.
! Written By Moetasim Ashfaq Dec-2005 @ PURDUE.EDU
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
  integer , dimension(4) :: icount , istart
  integer :: inet
  integer , dimension(6) :: ivar

  public :: get_cam2 , headcam2 , footercam2

  character(256) :: pathaddname
  type(rcm_time_and_date) :: refdate
  type(rcm_time_and_date) :: filedate
  logical :: ixfile

  data inet /-1/
  data ixfile /.false./
  data refdate /rcm_time_and_date(noleap,1989,12,27,0,0,0,32845,0)/

  contains
!
  subroutine headcam2
    use netcdf
!
    implicit none
!
    integer :: istatus , ivar1 , inet1 , ilat , ilon , ihyam , ihybm , ip0 , k
    integer :: lonid , latid , ilevid
    character(256) :: pathaddname
    real(8) :: dp0
!
    call zeit_ci('headcam2')
    pathaddname = trim(inpglob)//'/CAM2/USGS-gtopo30_0.9x1.25_remap_c051027.nc'
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
    if (istatus /= 0) call die('mod_cam2','Allocation error on glat',1)
    call mall_mci(glat,'mod_cam2')
    allocate(glon(nlon),stat=istatus)
    if (istatus /= 0) call die('mod_cam2','Allocation error on glon',1)
    call mall_mci(glon,'mod_cam2')
    allocate(zsvar(nlon,nlat),stat=istatus)
    if (istatus /= 0) call die('mod_cam2','Allocation error on zsvar',1)
    call mall_mci(zsvar,'mod_cam2')
    allocate(psvar(nlon,nlat),stat=istatus)
    if (istatus /= 0) call die('mod_cam2','Allocation error on psvar',1)
    call mall_mci(psvar,'mod_cam2')
    allocate(qvar(nlon,nlat,klev),stat=istatus)
    if (istatus /= 0) call die('mod_cam2','Allocation error on qsvar',1)
    call mall_mci(qvar,'mod_cam2')
    allocate(tvar(nlon,nlat,klev),stat=istatus)
    if (istatus /= 0) call die('mod_cam2','Allocation error on tvar',1)
    call mall_mci(tvar,'mod_cam2')
    allocate(hvar(nlon,nlat,klev),stat=istatus)
    if (istatus /= 0) call die('mod_cam2','Allocation error on hvar',1)
    call mall_mci(hvar,'mod_cam2')
    allocate(uvar(nlon,nlat,klev),stat=istatus)
    if (istatus /= 0) call die('mod_cam2','Allocation error on uvar',1)
    call mall_mci(uvar,'mod_cam2')
    allocate(vvar(nlon,nlat,klev),stat=istatus)
    if (istatus /= 0) call die('mod_cam2','Allocation error on vvar',1)
    call mall_mci(vvar,'mod_cam2')
    allocate(pp3d(nlon,nlat,klev),stat=istatus)
    if (istatus /= 0) call die('mod_cam2','Allocation error on pp3d',1)
    call mall_mci(pp3d,'mod_cam2')
 
    allocate(ak(klev),stat=istatus)
    if (istatus /= 0) call die('mod_cam2','Allocation error on ak',1)
    call mall_mci(ak,'mod_cam2')
    allocate(bk(klev),stat=istatus)
    if (istatus /= 0) call die('mod_cam2','Allocation error on bk',1)
    call mall_mci(bk,'mod_cam2')

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
    if (istatus /= 0) call die('mod_cam2','Allocation error on b3',1)
    call mall_mci(b3,'mod_cam2')
    allocate(d3(jx,iy,npl*2),stat=istatus)
    if (istatus /= 0) call die('mod_cam2','Allocation error on d3',1)
    call mall_mci(d3,'mod_cam2')
    allocate(b2(nlon,nlat,npl*3),stat=istatus)
    if (istatus /= 0) call die('mod_cam2','Allocation error on b2',1)
    call mall_mci(b2,'mod_cam2')
    allocate(d2(nlon,nlat,npl*2),stat=istatus)
    if (istatus /= 0) call die('mod_cam2','Allocation error on d2',1)
    call mall_mci(d2,'mod_cam2')

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

    call zeit_co('headcam2')
  end subroutine headcam2
!
!-----------------------------------------------------------------------
! 
  subroutine get_cam2(idate)

    use netcdf

    implicit none
!
    type(rcm_time_and_date) , intent(in) :: idate
!

    call zeit_ci('get_cam2')
    call readcam2(idate)

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
 
    call zeit_co('get_cam2')
  end subroutine get_cam2
!
!-----------------------------------------------------------------------
! 
  subroutine readcam2(idate)
!
    use netcdf
!
    implicit none
!
    type(rcm_time_and_date) , intent(in) :: idate
!
    integer :: istatus
    integer :: i , it , j , k , timid
    character(64) :: inname
    character(2) , dimension(6) :: varname
    real(dp) , allocatable , dimension(:) :: xtimes

    integer :: nscur , kkrec
    character(64) :: cunit , ccal
    type(rcm_time_interval) :: tdif
!
    data varname /'T' , 'Z3' , 'Q' , 'U' , 'V' , 'PS'/
!
    call zeit_ci('readcam2')
!
    if ( idate < itimes(1) .or. idate > itimes(timlen) ) then
      tdif = rcm_time_interval(7,uday)
      if (inet > 0) then
        istatus = nf90_close(inet)
        if ( istatus /= nf90_noerr ) call handle_err(istatus)
        filedate = filedate + tdif
        ixfile = .not. ixfile
      else
        filedate = refdate
        do while (idate <= filedate)
          filedate = filedate + tdif
          ixfile = .not. ixfile
        end do
      end if

      nscur = 0
      if (ixfile) nscur = 43200
      write (inname,99001) filedate%year, filedate%month, filedate%day, nscur
 
      pathaddname = trim(inpglob)//'/CAM2/'//inname
      istatus = nf90_open(pathaddname,nf90_nowrite,inet)
      if ( istatus /= nf90_noerr ) call handle_err(istatus)
      istatus = nf90_inq_dimid(inet,'time',timid)
      if ( istatus /= nf90_noerr ) call handle_err(istatus)
      istatus = nf90_inquire_dimension(inet,timid, len=timlen)
      if ( istatus /= nf90_noerr ) call handle_err(istatus)
      istatus = nf90_inq_varid(inet,'time',timid)
      if ( istatus /= nf90_noerr ) call handle_err(istatus)
      istatus = nf90_get_att(inet,timid,'units',cunit)
      if ( istatus /= nf90_noerr ) call handle_err(istatus)
      istatus = nf90_get_att(inet,timid,'calendar',ccal)
      if ( istatus /= nf90_noerr ) call handle_err(istatus)
      if (allocated(itimes)) then
        deallocate(itimes)
      end if
      allocate(xtimes(timlen),stat=istatus)
      if (istatus /= 0) call die('mod_cam2','Allocation error on xtimes',1)
      call mall_mci(xtimes,'mod_cam2')
      allocate(itimes(timlen),stat=istatus)
      if (istatus /= 0) call die('mod_cam2','Allocation error on itimes',1)
      istatus = nf90_get_var(inet,timid,xtimes)
      if ( istatus /= nf90_noerr ) call handle_err(istatus)
      do it = 1 , timlen
        itimes(it) = timeval2date(xtimes(it),cunit,ccal)
      end do
      deallocate(xtimes)
      call mall_mco(xtimes,'mod_cam2')
      do kkrec = 1 , 6
        istatus = nf90_inq_varid(inet,varname(kkrec), ivar(kkrec))
        if ( istatus /= nf90_noerr ) call handle_err(istatus)
        write (stdout,*) inet, trim(pathaddname), ' : ', varname(kkrec)
      end do
    end if

    tdif = idate - itimes(1)
    it = idnint(tdif%hours())/6 + 1

    icount(1) = nlon
    icount(2) = nlat
    icount(3) = 1
    istart(1) = 1
    istart(2) = 1
    istart(3) = it
    istatus = nf90_get_var(inet,ivar(6),psvar,istart(1:3),icount(1:3))
    if ( istatus /= nf90_noerr ) call handle_err(istatus)
    icount(1) = nlon
    icount(2) = nlat
    icount(3) = klev
    icount(4) = 1
    istart(1) = 1
    istart(2) = 1
    istart(3) = 1
    istart(4) = it
    istatus = nf90_get_var(inet,ivar(1),tvar,istart,icount)
    if ( istatus /= nf90_noerr ) call handle_err(istatus)
    istatus = nf90_get_var(inet,ivar(2),hvar,istart,icount)
    if ( istatus /= nf90_noerr ) call handle_err(istatus)
    istatus = nf90_get_var(inet,ivar(3),qvar,istart,icount)
    if ( istatus /= nf90_noerr ) call handle_err(istatus)
    istatus = nf90_get_var(inet,ivar(4),uvar,istart,icount)
    if ( istatus /= nf90_noerr ) call handle_err(istatus)
    istatus = nf90_get_var(inet,ivar(5),vvar,istart,icount)
    if ( istatus /= nf90_noerr ) call handle_err(istatus)

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
 
    call zeit_co('readcam2')

99001   format ('sococa.ts1.r1.cam2.h1.',i4,'-',i2,'-',i2,'-',i5,'.nc')

  end subroutine readcam2
!
!-----------------------------------------------------------------------
!
  subroutine footercam2
    implicit none
    call mall_mco(glat,'mod_cam2')
    deallocate(glat)
    call mall_mco(glon,'mod_cam2')
    deallocate(glon)
    call mall_mco(zsvar,'mod_cam2')
    deallocate(zsvar)
    call mall_mco(psvar,'mod_cam2')
    deallocate(psvar)
    call mall_mco(tvar,'mod_cam2')
    deallocate(tvar)
    call mall_mco(hvar,'mod_cam2')
    deallocate(hvar)
    call mall_mco(qvar,'mod_cam2')
    deallocate(qvar)
    call mall_mco(uvar,'mod_cam2')
    deallocate(uvar)
    call mall_mco(vvar,'mod_cam2')
    deallocate(vvar)
    call mall_mco(pp3d,'mod_cam2')
    deallocate(pp3d)
    call mall_mco(ak,'mod_cam2')
    deallocate(ak)
    call mall_mco(bk,'mod_cam2')
    deallocate(bk)
    call mall_mco(b3,'mod_cam2')
    deallocate(b3)
    call mall_mco(d3,'mod_cam2')
    deallocate(d3)
    call mall_mco(b2,'mod_cam2')
    deallocate(b2)
    call mall_mco(d2,'mod_cam2')
    deallocate(d2)
    deallocate(itimes)
  end subroutine footercam2
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
end module mod_cam2

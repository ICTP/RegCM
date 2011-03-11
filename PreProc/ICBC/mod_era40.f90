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

module mod_era40

  use mod_dynparam
  use m_realkinds
  use m_die
  use m_stdio
  use m_zeit
  use m_mall

  private

  integer , parameter :: klev = 23 , jlat = 73 , ilon = 144

  real(sp) , target , dimension(ilon,jlat,klev*3) :: b2
  real(sp) , target , dimension(ilon,jlat,klev*2) :: d2
  real(sp) , target , dimension(ilon,jlat,4*3+1) :: s2
  real(sp) , allocatable , target , dimension(:,:,:) :: b3
  real(sp) , allocatable , target , dimension(:,:,:) :: d3
  real(sp) , allocatable , target , dimension(:,:,:) :: s3

  real(sp) , dimension(ilon,jlat,klev) :: wvar

  real(sp) , pointer :: u3(:,:,:) , v3(:,:,:)
  real(sp) , pointer :: h3(:,:,:) , q3(:,:,:) , t3(:,:,:)
  real(sp) , pointer :: uvar(:,:,:) , vvar(:,:,:)
  real(sp) , pointer :: hvar(:,:,:) , rhvar(:,:,:) , tvar(:,:,:)
  real(sp) , pointer , dimension(:,:,:) :: qsoil , tsice , tsoil
  real(sp) , pointer , dimension(:,:) :: snw
  real(sp) , pointer , dimension(:,:,:) :: qs3 , ti3 , ts3
  real(sp) , pointer , dimension(:,:) :: snow

  real(sp) , dimension(jlat) :: glat
  real(sp) , dimension(ilon) :: glon
  real(sp) , dimension(klev) :: sigma1 , sigmar

  public :: getera40 , headerera , footerera

  contains

  subroutine getera40(idate)
  use mod_grid
  use mod_write
  use mod_interp , only : bilinx2
  use mod_vertint
  use mod_hgt
  use mod_humid
  use mod_mksst
  use mod_uvrot
  use mod_vectutil
  implicit none
!
  integer :: idate
!
!     D      BEGIN LOOP OVER NTIMES
!
  call zeit_ci('getera40')
  call era6hour(dattyp,idate,globidate1)
  write (stdout,*) 'READ IN fields at DATE:' , idate
!
!     HORIZONTAL INTERPOLATION OF BOTH THE SCALAR AND VECTOR FIELDS
!
  call bilinx2(b3,b2,xlon,xlat,glon,glat,ilon,jlat,jx,iy,klev*3)
  call bilinx2(d3,d2,dlon,dlat,glon,glat,ilon,jlat,jx,iy,klev*2)
!
!     ROTATE U-V FIELDS AFTER HORIZONTAL INTERPOLATION
!
  call uvrot4(u3,v3,dlon,dlat,clon,clat,grdfac,jx,iy,klev,plon,plat,iproj)
!
!     X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!     X X
!     V E R T I C A L   I N T E R P O L A T I O N
!
!     X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!     X X
!HH:  CHANGE THE VERTICAL ORDER.
  call top2btm(t3,jx,iy,klev)
  call top2btm(q3,jx,iy,klev)
  call top2btm(h3,jx,iy,klev)
  call top2btm(u3,jx,iy,klev)
  call top2btm(v3,jx,iy,klev)
!HH:OVER
!
!     ******           NEW CALCULATION OF P* ON RCM TOPOGRAPHY.
  call intgtb(pa,za,tlayer,topogm,t3,h3,sigmar,jx,iy,klev)
 
  call intpsn(ps4,topogm,pa,za,tlayer,ptop,jx,iy)
  if(i_band == 1) then
     call p1p2_band(b3pd,ps4,jx,iy)
  else
     call p1p2(b3pd,ps4,jx,iy)
  endif
 
!
!     F0  DETERMINE SURFACE TEMPS ON RCM TOPOGRAPHY.
!     INTERPOLATION FROM PRESSURE LEVELS AS IN INTV2
  call intv3(ts4,t3,ps4,sigmar,ptop,jx,iy,klev)
 
  call readsst(ts4,idate)

!     F2  DETERMINE P* AND HEIGHT.
!
!     F3  INTERPOLATE U, V, T, AND Q.
  call intv1(u4,u3,b3pd,sigma2,sigmar,ptop,jx,iy,kz,klev)
  call intv1(v4,v3,b3pd,sigma2,sigmar,ptop,jx,iy,kz,klev)
!
  call intv2(t4,t3,ps4,sigma2,sigmar,ptop,jx,iy,kz,klev)
 
  call intv1(q4,q3,ps4,sigma2,sigmar,ptop,jx,iy,kz,klev)
  call humid2(t4,q4,ps4,ptop,sigma2,jx,iy,kz)
!
!     F4  DETERMINE H
  call hydrost(h4,t4,topogm,ps4,ptop,sigmaf,sigma2,dsigma,jx,iy,kz)
!
  call zeit_co('getera40')
  end subroutine getera40

  subroutine era6hour(dattyp,idate,idate0)
  use netcdf
  implicit none
!
  character(5) :: dattyp
  integer :: idate , idate0
  intent (in) dattyp , idate , idate0
!
  integer :: i , inet , it , j , k , k4 , kkrec , month , nday , &
             nhour , nyear , istatus
  character(24) :: inname
  character(256) :: pathaddname
!     character(5) , dimension(3,4) :: sarname
!     character(2) :: snownm
  logical :: there
!     character(5) , dimension(6) :: varname
  integer(2) , dimension(ilon,jlat,klev) :: work
  real(dp) :: xadd , xscale

  integer , dimension(10) , save :: icount , istart
  real(dp) , dimension(5,4) , save :: xoff , xscl
  integer , dimension(5,4) , save :: inet6
!
!     This is the latitude, longitude dimension of the grid to be read.
!     This corresponds to the lat and lon dimension variables in the
!     netCDF file.
!
!     The data are packed into short integers (INTEGER*2).  The array
!     work will be used to hold the packed integers.  The array 'x'
!     will contain the unpacked data.
!
!     DATA ARRAY AND WORK ARRAY
!
!     data varname/'air' , 'hgt' , 'rhum' , 'uwnd' , 'vwnd' , 'omega'/
!     data sarname/'swvl1' , 'istl1' , 'stl1' , 'swvl2' , 'istl2' ,  &
!         &'stl2' , 'swvl3' , 'istl3' , 'stl3' , 'swvl4' , 'istl4' , &
!         &'stl4'/
!     data snownm/'sd'/
!
!     Below in the ncopen call is the file name of the netCDF file.
!     You may want to add code to read in the file name and the
!     variable name.
!     OPEN FILE AND GET FILES ID AND VARIABLE ID(S)
!
!bxq
  if ( idate < 1957090100 .or. idate > 2002083118 ) then
    call die('getera40', 'ERA40 datasets is just available from'// &
             ' 1957090100 to 2002083118', 1)
  end if
 
  nyear = idate/1000000
  month = idate/10000 - nyear*100
  nday = idate/100 - nyear*10000 - month*100
  nhour = idate - nyear*1000000 - month*10000 - nday*100
  if ( idate == idate0 .or.                                           &
       (mod(idate,100000) == 10100 .and. mod(idate,1000000) /= 110100) ) then
    do k4 = 1 , 4
      do kkrec = 1 , 5
        if ( kkrec == 1 ) then
          if ( k4 == 1 ) then
            write (inname,99001) nyear , 'air.' , nyear
          else if ( k4 == 2 ) then
            write (inname,99002) nyear , 'air.' , nyear
          else if ( k4 == 3 ) then
            write (inname,99003) nyear , 'air.' , nyear
          else if ( k4 == 4 ) then
            write (inname,99004) nyear , 'air.' , nyear
          else
          end if
        else if ( kkrec == 2 ) then
          if ( k4 == 1 ) then
            write (inname,99001) nyear , 'hgt.' , nyear
          else if ( k4 == 2 ) then
            write (inname,99002) nyear , 'hgt.' , nyear
          else if ( k4 == 3 ) then
            write (inname,99003) nyear , 'hgt.' , nyear
          else if ( k4 == 4 ) then
            write (inname,99004) nyear , 'hgt.' , nyear
          else
          end if
        else if ( kkrec == 3 ) then
          if ( k4 == 1 ) then
            write (inname,99005) nyear , 'rhum.' , nyear
          else if ( k4 == 2 ) then
            write (inname,99006) nyear , 'rhum.' , nyear
          else if ( k4 == 3 ) then
            write (inname,99007) nyear , 'rhum.' , nyear
          else if ( k4 == 4 ) then
            write (inname,99008) nyear , 'rhum.' , nyear
          else
          end if
        else if ( kkrec == 4 ) then
          if ( k4 == 1 ) then
            write (inname,99005) nyear , 'uwnd.' , nyear
          else if ( k4 == 2 ) then
            write (inname,99006) nyear , 'uwnd.' , nyear
          else if ( k4 == 3 ) then
            write (inname,99007) nyear , 'uwnd.' , nyear
          else if ( k4 == 4 ) then
            write (inname,99008) nyear , 'uwnd.' , nyear
          else
          end if
        else if ( kkrec == 5 ) then
          if ( k4 == 1 ) then
            write (inname,99005) nyear , 'vwnd.' , nyear
          else if ( k4 == 2 ) then
            write (inname,99006) nyear , 'vwnd.' , nyear
          else if ( k4 == 3 ) then
            write (inname,99007) nyear , 'vwnd.' , nyear
          else if ( k4 == 4 ) then
            write (inname,99008) nyear , 'vwnd.' , nyear
          else
          end if
        else if ( kkrec == 6 ) then
          if ( k4 == 1 ) then
            write (inname,99009) nyear , 'omega.' , nyear
          else if ( k4 == 2 ) then
            write (inname,99010) nyear , 'omega.' , nyear
          else if ( k4 == 3 ) then
            write (inname,99011) nyear , 'omega.' , nyear
          else if ( k4 == 4 ) then
            write (inname,99012) nyear , 'omega.' , nyear
          else
          end if
        else
        end if
 
        pathaddname = trim(inpglob)//dattyp//'/'//inname
        inquire (file=pathaddname,exist=there)
        if ( .not.there ) then
          call die('getera40', trim(pathaddname)//' is not available', 1)
        end if
        istatus = nf90_open(pathaddname,nf90_nowrite,inet6(kkrec,k4))
        if ( istatus /= nf90_noerr ) then
          call die('getera40', trim(pathaddname), 1, &
                   nf90_strerror(istatus), istatus)
        end if
        istatus = nf90_get_att(inet6(kkrec,k4),5,'scale_factor',xscl(kkrec,k4))
        if ( istatus /= nf90_noerr ) then
          call die('getera40', trim(pathaddname)//': scale_factor', 1, &
                   nf90_strerror(istatus), istatus)
        end if
        istatus = nf90_get_att(inet6(kkrec,k4),5,'add_offset',xoff(kkrec,k4))
        if ( istatus /= nf90_noerr ) then
          call die('getera40', trim(pathaddname)//': add_offset', 1, &
                   nf90_strerror(istatus), istatus)
        end if
        write (stdout,*) inet6(kkrec,k4) , trim(pathaddname) ,  &
                    xscl(kkrec,k4) , xoff(kkrec,k4)
      end do
    end do
 
  end if
 
  k4 = nhour/6 + 1
  it = nday
  if ( month == 2 ) it = it + 31
  if ( month == 3 ) it = it + 59
  if ( month == 4 ) it = it + 90
  if ( month == 5 ) it = it + 120
  if ( month == 6 ) it = it + 151
  if ( month == 7 ) it = it + 181
  if ( month == 8 ) it = it + 212
  if ( month == 9 ) it = it + 243
  if ( month == 10 ) it = it + 273
  if ( month == 11 ) it = it + 304
  if ( month == 12 ) it = it + 334
  if ( mod(nyear,4) == 0 .and. month > 2 ) it = it + 1
  if ( mod(nyear,100) == 0 .and. month > 2 ) it = it - 1
  if ( mod(nyear,400) == 0 .and. month > 2 ) it = it + 1
  do k = 1 , 4
    istart(k) = 1
  end do
  do k = 5 , 10
    istart(k) = 0
    icount(k) = 0
  end do
  icount(1) = ilon
  icount(2) = jlat
  icount(3) = klev
  icount(4) = 365
  if ( mod(nyear,4) == 0 ) icount(4) = 366
  if ( mod(nyear,100) == 0 ) icount(4) = 365
  if ( mod(nyear,400) == 0 ) icount(4) = 366
  if ( nyear == 2002 ) icount(4) = 243
  if ( nyear == 1957 ) icount(4) = 122
  if ( nyear == 1957 ) it = it - 243
  istart(4) = it
  icount(4) = 1
!bxq_
  do kkrec = 1 , 5
    inet = inet6(kkrec,k4)
    istatus = nf90_get_var(inet,5,work,istart,icount)
    if ( istatus /= nf90_noerr ) then
      call die('getera40', trim(pathaddname)//': getvar', 1, &
               nf90_strerror(istatus), istatus)
    end if
    xscale = xscl(kkrec,k4)
    xadd = xoff(kkrec,k4)
    if ( kkrec == 1 ) then
      do k = 1 , klev
        do j = 1 , jlat
          do i = 1 , ilon
            tvar(i,jlat+1-j,k) = real(dble(work(i,j,k))*xscale + xadd)
          end do
        end do
      end do
    else if ( kkrec == 2 ) then
      do k = 1 , klev
        do j = 1 , jlat
          do i = 1 , ilon
            hvar(i,jlat+1-j,k) = real(dble(work(i,j,k))*xscale + xadd)
            hvar(i,jlat+1-j,k) = hvar(i,jlat+1-j,k)/9.80616
          end do
        end do
      end do
    else if ( kkrec == 3 ) then
      do k = 1 , klev
        do j = 1 , jlat
          do i = 1 , ilon
            rhvar(i,jlat+1-j,k) = real(dble(work(i,j,k))*xscale + xadd)
            rhvar(i,jlat+1-j,k) = rhvar(i,jlat+1-j,k)*0.01
!               RHvar(i,jlat+1-j,k) = amax1(RHvar(i,jlat+1-j,k),1.05)
          end do
        end do
      end do
    else if ( kkrec == 4 ) then
      do k = 1 , klev
        do j = 1 , jlat
          do i = 1 , ilon
            uvar(i,jlat+1-j,k) = real(dble(work(i,j,k))*xscale + xadd)
          end do
        end do
      end do
    else if ( kkrec == 5 ) then
      do k = 1 , klev
        do j = 1 , jlat
          do i = 1 , ilon
            vvar(i,jlat+1-j,k) = real(dble(work(i,j,k))*xscale + xadd)
          end do
        end do
      end do
    else if ( kkrec == 6 ) then
      do k = 1 , klev
        do j = 1 , jlat
          do i = 1 , ilon
            wvar(i,jlat+1-j,k) = real(dble(work(i,j,k))*xscale + xadd)
          end do
        end do
      end do
    end if
  end do
 
99001 format (i4,'/',a4,i4,'.00.nc')
99002 format (i4,'/',a4,i4,'.06.nc')
99003 format (i4,'/',a4,i4,'.12.nc')
99004 format (i4,'/',a4,i4,'.18.nc')
99005 format (i4,'/',a5,i4,'.00.nc')
99006 format (i4,'/',a5,i4,'.06.nc')
99007 format (i4,'/',a5,i4,'.12.nc')
99008 format (i4,'/',a5,i4,'.18.nc')
99009 format (i4,'/',a6,i4,'.00.nc')
99010 format (i4,'/',a6,i4,'.06.nc')
99011 format (i4,'/',a6,i4,'.12.nc')
99012 format (i4,'/',a6,i4,'.18.nc')
!
  end subroutine era6hour

  subroutine headerera
  implicit none
!
  integer :: i , j , k , kr , ierr
!
  sigmar(1) = .001
  sigmar(2) = .002
  sigmar(3) = .003
  sigmar(4) = .005
  sigmar(5) = .007
  sigmar(6) = .01
  sigmar(7) = .02
  sigmar(8) = .03
  sigmar(9) = .05
  sigmar(10) = .07
  sigmar(11) = .1
  sigmar(12) = .15
  sigmar(13) = .2
  sigmar(14) = .25
  sigmar(15) = .3
  sigmar(16) = .4
  sigmar(17) = .5
  sigmar(18) = .6
  sigmar(19) = .7
  sigmar(20) = .775
  sigmar(21) = .85
  sigmar(22) = .925
  sigmar(23) = 1.00
!
!     INITIAL GLOBAL GRID-POINT LONGITUDE & LATITUDE
!
  do i = 1 , ilon
    glon(i) = float(i-1)*2.5
  end do
  do j = 1 , jlat
    glat(j) = -90.0 + float(j-1)*2.5
  end do
!HH:OVER
!     CHANGE ORDER OF VERTICAL INDEXES FOR PRESSURE LEVELS
!
  do k = 1 , klev
    kr = klev - k + 1
    sigma1(k) = sigmar(kr)
  end do
 
  allocate(b3(jx,iy,klev*3), stat=ierr)
  if (ierr /= 0) call die('headerera','allocate b3',ierr)
  call mall_mci(b3,'mod_era40')
  allocate(d3(jx,iy,klev*2), stat=ierr)
  if (ierr /= 0) call die('headerera','allocate d3',ierr)
  call mall_mci(d3,'mod_era40')
  allocate(s3(jx,iy,4*3+1), stat=ierr)
  if (ierr /= 0) call die('headerera','allocate s3',ierr)
  call mall_mci(s3,'mod_era40')

!     Set up pointers

  u3 => d3(:,:,1:klev)
  v3 => d3(:,:,klev+1:2*klev)
  t3 => b3(:,:,1:klev)
  h3 => b3(:,:,klev+1:2*klev)
  q3 => b3(:,:,2*klev+1:3*klev)
  qs3 => s3(:,:,1:4)
  ti3 => s3(:,:,5:8)
  ts3 => s3(:,:,9:12)
  snow => s3(:,:,13)
  uvar => d2(:,:,1:klev)
  vvar => d2(:,:,klev+1:2*klev)
  tvar => b2(:,:,1:klev)
  hvar => b2(:,:,klev+1:2*klev)
  rhvar => b2(:,:,2*klev+1:3*klev)
  qsoil => s2(:,:,1:4)
  tsice => s2(:,:,5:8)
  tsoil => s2(:,:,9:12)
  snw => s2(:,:,13)

  end subroutine headerera

  subroutine footerera
    implicit none
    call mall_mco(d3,'mod_era40')
    deallocate(d3)
    call mall_mco(b3,'mod_era40')
    deallocate(b3)
    call mall_mco(s3,'mod_era40')
    deallocate(s3)
  end subroutine footerera

end module mod_era40

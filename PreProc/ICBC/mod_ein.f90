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

  use mod_dynparam
  use m_realkinds
  use m_stdio
  use m_die
  use m_zeit
  use mod_memutil
  use mod_grid
  use mod_write
  use mod_interp
  use mod_vertint
  use mod_hgt
  use mod_humid
  use mod_mksst
  use mod_uvrot
  use mod_vectutil

  private

  integer :: inlev , klev , jlat , ilon , imindat , imaxdat

  real(sp) , pointer , dimension(:,:,:) :: b3
  real(sp) , pointer , dimension(:,:,:) :: d3

  real(sp) , pointer :: u3(:,:,:) , v3(:,:,:)
  real(sp) , pointer :: h3(:,:,:) , q3(:,:,:) , t3(:,:,:)
  real(sp) , pointer :: uvar(:,:,:) , vvar(:,:,:)
  real(sp) , pointer :: hvar(:,:,:) , rhvar(:,:,:) , tvar(:,:,:)

  real(sp) :: xres
  real(sp) , pointer , dimension(:,:,:) :: b2
  real(sp) , pointer , dimension(:,:,:) :: d2
  real(sp) , pointer , dimension(:) :: glat
  real(sp) , pointer , dimension(:) :: glon
  real(sp) , pointer , dimension(:) :: sigma1 , sigmar
  integer(2) , pointer , dimension(:,:,:) :: work

  integer , dimension(5,4) :: inet5
  integer , dimension(5,4) :: ivar5
  real(dp) , dimension(5,4) :: xoff , xscl

  public :: getein , headerein

  contains

  subroutine getein(idate)
  implicit none
!
  type(rcm_time_and_date) , intent(in) :: idate
!
!     READ DATA AT IDATE
!
  call zeit_ci('get')

  call ein6hour(dattyp,idate,globidate1)

  write (stdout,*) 'READ IN fields at DATE:' , idate%tostring()
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
!
  call intgtb(pa,za,tlayer,topogm,t3,h3,sigmar,jx,iy,klev)
 
  call intpsn(ps4,topogm,pa,za,tlayer,ptop,jx,iy)
  if(i_band == 1) then
     call p1p2_band(b3pd,ps4,jx,iy)
  else
     call p1p2(b3pd,ps4,jx,iy)
  endif
 
!     INTERPOLATION FROM PRESSURE LEVELS
  call intv3(ts4,t3,ps4,sigmar,ptop,jx,iy,klev)

  call readsst(ts4,idate)
 
!     F3  INTERPOLATE U, V, T, AND Q.

  call intv1(u4,u3,b3pd,sigma2,sigmar,ptop,jx,iy,kz,klev)
  call intv1(v4,v3,b3pd,sigma2,sigmar,ptop,jx,iy,kz,klev)
  call intv2(t4,t3,ps4,sigma2,sigmar,ptop,jx,iy,kz,klev)
  call intv1(q4,q3,ps4,sigma2,sigmar,ptop,jx,iy,kz,klev)
! 
  call humid2(t4,q4,ps4,ptop,sigma2,jx,iy,kz)
!
!     F4  DETERMINE H
  call hydrost(h4,t4,topogm,ps4,ptop,sigmaf,sigma2,dsigma,jx,iy,kz)
!
  call zeit_co('get')
!
  end subroutine getein
!
!-----------------------------------------------------------------------
!
  subroutine ein6hour(dattyp,idate,idate0)
  use netcdf
  implicit none
!
  character(5) , intent(in) :: dattyp
  type(rcm_time_and_date) , intent(in) :: idate , idate0
!
  integer :: i , inet , it , j , k , k4 , kkrec , istatus , ivar
  character(24) :: inname
  character(256) :: pathaddname
  character(1) , dimension(5) :: varname
  real(dp) :: xadd , xscale
  integer , dimension(10) :: icount , istart
!
!     This is the latitude, longitude dimension of the grid to be read.
!     This corresponds to the lat and lon dimension variables in the
!     netCDF file.
!     The data are packed into short integers (INTEGER*2).  The array
!     work will be used to hold the packed integers.  The array 'x'
!     will contain the unpacked data.
!
  data varname/'t' , 'z' , 'r' , 'u' , 'v'/
!
  call zeit_ci('readein')
!
  if ( idate < imindat .or. idate > imaxdat ) then
    write (stderr, *) 'EIN data for resolution ',xres,' degrees ', &
                  'are available only from ',imindat,' up to ', imaxdat
    call die('ein6hour','EIN dataset unavailable',1)
  end if
 
  if ( idate == idate0 .or. (lfdoyear(idate) .and. lmidnight(idate))) then
    do k4 = 1 , 4
      do kkrec = 1 , 5
        if ( kkrec == 1 ) then
          if ( k4 == 1 ) then
            write (inname,99001) idate%year , 'air.' , idate%year
          else if ( k4 == 2 ) then
            write (inname,99002) idate%year , 'air.' , idate%year
          else if ( k4 == 3 ) then
            write (inname,99003) idate%year , 'air.' , idate%year
          else if ( k4 == 4 ) then
            write (inname,99004) idate%year , 'air.' , idate%year
          end if
        else if ( kkrec == 2 ) then
          if ( k4 == 1 ) then
            write (inname,99001) idate%year , 'hgt.' , idate%year
          else if ( k4 == 2 ) then
            write (inname,99002) idate%year , 'hgt.' , idate%year
          else if ( k4 == 3 ) then
            write (inname,99003) idate%year , 'hgt.' , idate%year
          else if ( k4 == 4 ) then
            write (inname,99004) idate%year , 'hgt.' , idate%year
          end if
        else if ( kkrec == 3 ) then
          if ( k4 == 1 ) then
            write (inname,99005) idate%year , 'rhum.' , idate%year
          else if ( k4 == 2 ) then
            write (inname,99006) idate%year , 'rhum.' , idate%year
          else if ( k4 == 3 ) then
            write (inname,99007) idate%year , 'rhum.' , idate%year
          else if ( k4 == 4 ) then
            write (inname,99008) idate%year , 'rhum.' , idate%year
          end if
        else if ( kkrec == 4 ) then
          if ( k4 == 1 ) then
            write (inname,99005) idate%year , 'uwnd.' , idate%year
          else if ( k4 == 2 ) then
            write (inname,99006) idate%year , 'uwnd.' , idate%year
          else if ( k4 == 3 ) then
            write (inname,99007) idate%year , 'uwnd.' , idate%year
          else if ( k4 == 4 ) then
            write (inname,99008) idate%year , 'uwnd.' , idate%year
          end if
        else if ( kkrec == 5 ) then
          if ( k4 == 1 ) then
            write (inname,99005) idate%year , 'vwnd.' , idate%year
          else if ( k4 == 2 ) then
            write (inname,99006) idate%year , 'vwnd.' , idate%year
          else if ( k4 == 3 ) then
            write (inname,99007) idate%year , 'vwnd.' , idate%year
          else if ( k4 == 4 ) then
            write (inname,99008) idate%year , 'vwnd.' , idate%year
          end if
        end if
 
        pathaddname = trim(inpglob)//pthsep//dattyp//pthsep//inname
        istatus = nf90_open(pathaddname,nf90_nowrite,inet5(kkrec,k4))
        if ( istatus /= nf90_noerr) then
          call die('ein6hour',trim(pathaddname)//':open',1, &
                    nf90_strerror(istatus),istatus)
        end if
        istatus = nf90_inq_varid(inet5(kkrec,k4),varname(kkrec),ivar5(kkrec,k4))
        if ( istatus /= nf90_noerr) then
          call die('ein6hour',trim(pathaddname)//':'//varname(kkrec),1, &
                   nf90_strerror(istatus),istatus)
        end if
        istatus = nf90_get_att(inet5(kkrec,k4),ivar5(kkrec,k4), &
                 'scale_factor',xscl(kkrec,k4))
        if ( istatus /= nf90_noerr) then
          call die('ein6hour',trim(pathaddname)//':'// &
                   varname(kkrec)//':scale_factor',1,  &
                   nf90_strerror(istatus),istatus)
        end if
        istatus = nf90_get_att(inet5(kkrec,k4),ivar5(kkrec,k4),  &
                 'add_offset',xoff(kkrec,k4))
        if ( istatus /= nf90_noerr) then
          call die('ein6hour',trim(pathaddname)//':'// &
                   varname(kkrec)//':add_offset',1,    &
                   nf90_strerror(istatus),istatus)
        end if
        write (stdout,*) inet5(kkrec,k4) , trim(pathaddname) ,   &
                         xscl(kkrec,k4) , xoff(kkrec,k4)
      end do
    end do
 
  end if
 
  k4 = idate%hour/6 + 1
  it = idate%day
  if ( idate%month == 2 ) it = it + 31
  if ( idate%month == 3 ) it = it + 59
  if ( idate%month == 4 ) it = it + 90
  if ( idate%month == 5 ) it = it + 120
  if ( idate%month == 6 ) it = it + 151
  if ( idate%month == 7 ) it = it + 181
  if ( idate%month == 8 ) it = it + 212
  if ( idate%month == 9 ) it = it + 243
  if ( idate%month == 10 ) it = it + 273
  if ( idate%month == 11 ) it = it + 304
  if ( idate%month == 12 ) it = it + 334
  if ( mod(idate%year,4) == 0 .and. idate%month > 2 ) it = it + 1
  if ( mod(idate%year,100) == 0 .and. idate%month > 2 ) it = it - 1
  if ( mod(idate%year,400) == 0 .and. idate%month > 2 ) it = it + 1
  do k = 1 , 4
    istart(k) = 1
  end do
  do k = 5 , 10
    istart(k) = 0
    icount(k) = 0
  end do
  icount(1) = ilon
  icount(2) = jlat
  icount(3) = inlev
  icount(4) = 365
  if ( mod(idate%year,4) == 0 ) icount(4) = 366
  if ( mod(idate%year,100) == 0 ) icount(4) = 365
  if ( mod(idate%year,400) == 0 ) icount(4) = 366
  istart(4) = it
  icount(4) = 1
!bxq_
  do kkrec = 1 , 5
    inet = inet5(kkrec,k4)
    ivar = ivar5(kkrec,k4)
    istatus = nf90_get_var(inet,ivar,work,istart,icount)
    if (istatus /= nf90_noerr) then
      call die('ein6hour',trim(pathaddname)//':'// &
               varname(kkrec)//':readvar',1,       &
               nf90_strerror(istatus),istatus)
    end if
    xscale = xscl(kkrec,k4)
    xadd = xoff(kkrec,k4)
    if ( kkrec == 1 ) then
      do j = 1 , jlat
        do i = 1 , ilon
!             Tvar(i,jlat+1-j,k) = work(i,j,k)*xscale+xadd
          tvar(i,jlat+1-j,1) = real(dble(work(i,j,1))*xscale+xadd)
          tvar(i,jlat+1-j,2) = real(dble(work(i,j,2))*xscale+xadd)
          tvar(i,jlat+1-j,3) = real(dble(work(i,j,3))*xscale+xadd)
          tvar(i,jlat+1-j,4) = real(dble(work(i,j,4))*xscale+xadd)
          tvar(i,jlat+1-j,5) = real(dble(work(i,j,5))*xscale+xadd)
          tvar(i,jlat+1-j,6) = real(dble(work(i,j,6))*xscale+xadd)
          tvar(i,jlat+1-j,7) = real(dble(work(i,j,7))*xscale+xadd)
          tvar(i,jlat+1-j,8) = real(dble(work(i,j,8))*xscale+xadd)
          tvar(i,jlat+1-j,9) = real(dble(work(i,j,9))*xscale+xadd)
          tvar(i,jlat+1-j,10) = real(dble(work(i,j,10))*xscale+xadd)
          tvar(i,jlat+1-j,11) = real(dble(work(i,j,11))*xscale+xadd)
          tvar(i,jlat+1-j,12) = real(dble(work(i,j,13))*xscale+xadd)
          tvar(i,jlat+1-j,13) = real(dble(work(i,j,15))*xscale+xadd)
          tvar(i,jlat+1-j,14) = real(dble(work(i,j,17))*xscale+xadd)
          tvar(i,jlat+1-j,15) = real(dble(work(i,j,18))*xscale+xadd)
          tvar(i,jlat+1-j,16) = real(dble(work(i,j,20))*xscale+xadd)
          tvar(i,jlat+1-j,17) = real(dble(work(i,j,22))*xscale+xadd)
          tvar(i,jlat+1-j,18) = real(dble(work(i,j,24))*xscale+xadd)
          tvar(i,jlat+1-j,19) = real(dble(work(i,j,26))*xscale+xadd)
          tvar(i,jlat+1-j,20) = real(dble(work(i,j,28))*xscale+xadd)
          tvar(i,jlat+1-j,21) = real(dble(work(i,j,31))*xscale+xadd)
          tvar(i,jlat+1-j,22) = real(dble(work(i,j,34))*xscale+xadd)
          tvar(i,jlat+1-j,23) = real(dble(work(i,j,37))*xscale+xadd)
        end do
      end do
    else if ( kkrec == 2 ) then
      do j = 1 , jlat
        do i = 1 , ilon
!             Hvar(i,jlat+1-j,k) = work(i,j,k)*xscale+xadd
          hvar(i,jlat+1-j,1) = real(dble(work(i,j,1))*xscale+xadd)
          hvar(i,jlat+1-j,2) = real(dble(work(i,j,2))*xscale+xadd)
          hvar(i,jlat+1-j,3) = real(dble(work(i,j,3))*xscale+xadd)
          hvar(i,jlat+1-j,4) = real(dble(work(i,j,4))*xscale+xadd)
          hvar(i,jlat+1-j,5) = real(dble(work(i,j,5))*xscale+xadd)
          hvar(i,jlat+1-j,6) = real(dble(work(i,j,6))*xscale+xadd)
          hvar(i,jlat+1-j,7) = real(dble(work(i,j,7))*xscale+xadd)
          hvar(i,jlat+1-j,8) = real(dble(work(i,j,8))*xscale+xadd)
          hvar(i,jlat+1-j,9) = real(dble(work(i,j,9))*xscale+xadd)
          hvar(i,jlat+1-j,10) = real(dble(work(i,j,10))*xscale+xadd)
          hvar(i,jlat+1-j,11) = real(dble(work(i,j,11))*xscale+xadd)
          hvar(i,jlat+1-j,12) = real(dble(work(i,j,13))*xscale+xadd)
          hvar(i,jlat+1-j,13) = real(dble(work(i,j,15))*xscale+xadd)
          hvar(i,jlat+1-j,14) = real(dble(work(i,j,17))*xscale+xadd)
          hvar(i,jlat+1-j,15) = real(dble(work(i,j,18))*xscale+xadd)
          hvar(i,jlat+1-j,16) = real(dble(work(i,j,20))*xscale+xadd)
          hvar(i,jlat+1-j,17) = real(dble(work(i,j,22))*xscale+xadd)
          hvar(i,jlat+1-j,18) = real(dble(work(i,j,24))*xscale+xadd)
          hvar(i,jlat+1-j,19) = real(dble(work(i,j,26))*xscale+xadd)
          hvar(i,jlat+1-j,20) = real(dble(work(i,j,28))*xscale+xadd)
          hvar(i,jlat+1-j,21) = real(dble(work(i,j,31))*xscale+xadd)
          hvar(i,jlat+1-j,22) = real(dble(work(i,j,34))*xscale+xadd)
          hvar(i,jlat+1-j,23) = real(dble(work(i,j,37))*xscale+xadd)
        end do
      end do
      do k = 1 , klev
        do j = 1 , jlat
          do i = 1 , ilon
            hvar(i,j,k) = hvar(i,j,k)/9.80616
          end do
        end do
      end do
    else if ( kkrec == 3 ) then
      do j = 1 , jlat
        do i = 1 , ilon
!             RHvar(i,jlat+1-j,k) = work(i,j,k)*xscale+xadd
          rhvar(i,jlat+1-j,1) = real(dble(work(i,j,1))*xscale+xadd)
          rhvar(i,jlat+1-j,2) = real(dble(work(i,j,2))*xscale+xadd)
          rhvar(i,jlat+1-j,3) = real(dble(work(i,j,3))*xscale+xadd)
          rhvar(i,jlat+1-j,4) = real(dble(work(i,j,4))*xscale+xadd)
          rhvar(i,jlat+1-j,5) = real(dble(work(i,j,5))*xscale+xadd)
          rhvar(i,jlat+1-j,6) = real(dble(work(i,j,6))*xscale+xadd)
          rhvar(i,jlat+1-j,7) = real(dble(work(i,j,7))*xscale+xadd)
          rhvar(i,jlat+1-j,8) = real(dble(work(i,j,8))*xscale+xadd)
          rhvar(i,jlat+1-j,9) = real(dble(work(i,j,9))*xscale+xadd)
          rhvar(i,jlat+1-j,10) = real(dble(work(i,j,10))*xscale+xadd)
          rhvar(i,jlat+1-j,11) = real(dble(work(i,j,11))*xscale+xadd)
          rhvar(i,jlat+1-j,12) = real(dble(work(i,j,13))*xscale+xadd)
          rhvar(i,jlat+1-j,13) = real(dble(work(i,j,15))*xscale+xadd)
          rhvar(i,jlat+1-j,14) = real(dble(work(i,j,17))*xscale+xadd)
          rhvar(i,jlat+1-j,15) = real(dble(work(i,j,18))*xscale+xadd)
          rhvar(i,jlat+1-j,16) = real(dble(work(i,j,20))*xscale+xadd)
          rhvar(i,jlat+1-j,17) = real(dble(work(i,j,22))*xscale+xadd)
          rhvar(i,jlat+1-j,18) = real(dble(work(i,j,24))*xscale+xadd)
          rhvar(i,jlat+1-j,19) = real(dble(work(i,j,26))*xscale+xadd)
          rhvar(i,jlat+1-j,20) = real(dble(work(i,j,28))*xscale+xadd)
          rhvar(i,jlat+1-j,21) = real(dble(work(i,j,31))*xscale+xadd)
          rhvar(i,jlat+1-j,22) = real(dble(work(i,j,34))*xscale+xadd)
          rhvar(i,jlat+1-j,23) = real(dble(work(i,j,37))*xscale+xadd)
        end do
      end do
      do k = 1 , 23
        do j = 1 , jlat
          do i = 1 , ilon
            rhvar(i,j,k) = amax1(rhvar(i,j,k)*0.01,0.00)
          end do
        end do
      end do
    else if ( kkrec == 4 ) then
      do j = 1 , jlat
        do i = 1 , ilon
!             Uvar(i,jlat+1-j,k) = work(i,j,k)*xscale+xadd
          uvar(i,jlat+1-j,1) = real(dble(work(i,j,1))*xscale+xadd)
          uvar(i,jlat+1-j,2) = real(dble(work(i,j,2))*xscale+xadd)
          uvar(i,jlat+1-j,3) = real(dble(work(i,j,3))*xscale+xadd)
          uvar(i,jlat+1-j,4) = real(dble(work(i,j,4))*xscale+xadd)
          uvar(i,jlat+1-j,5) = real(dble(work(i,j,5))*xscale+xadd)
          uvar(i,jlat+1-j,6) = real(dble(work(i,j,6))*xscale+xadd)
          uvar(i,jlat+1-j,7) = real(dble(work(i,j,7))*xscale+xadd)
          uvar(i,jlat+1-j,8) = real(dble(work(i,j,8))*xscale+xadd)
          uvar(i,jlat+1-j,9) = real(dble(work(i,j,9))*xscale+xadd)
          uvar(i,jlat+1-j,10) = real(dble(work(i,j,10))*xscale+xadd)
          uvar(i,jlat+1-j,11) = real(dble(work(i,j,11))*xscale+xadd)
          uvar(i,jlat+1-j,12) = real(dble(work(i,j,13))*xscale+xadd)
          uvar(i,jlat+1-j,13) = real(dble(work(i,j,15))*xscale+xadd)
          uvar(i,jlat+1-j,14) = real(dble(work(i,j,17))*xscale+xadd)
          uvar(i,jlat+1-j,15) = real(dble(work(i,j,18))*xscale+xadd)
          uvar(i,jlat+1-j,16) = real(dble(work(i,j,20))*xscale+xadd)
          uvar(i,jlat+1-j,17) = real(dble(work(i,j,22))*xscale+xadd)
          uvar(i,jlat+1-j,18) = real(dble(work(i,j,24))*xscale+xadd)
          uvar(i,jlat+1-j,19) = real(dble(work(i,j,26))*xscale+xadd)
          uvar(i,jlat+1-j,20) = real(dble(work(i,j,28))*xscale+xadd)
          uvar(i,jlat+1-j,21) = real(dble(work(i,j,31))*xscale+xadd)
          uvar(i,jlat+1-j,22) = real(dble(work(i,j,34))*xscale+xadd)
          uvar(i,jlat+1-j,23) = real(dble(work(i,j,37))*xscale+xadd)
        end do
      end do
    else if ( kkrec == 5 ) then
      do j = 1 , jlat
        do i = 1 , ilon
!             Vvar(i,jlat+1-j,k) = work(i,j,k)*xscale+xadd
          vvar(i,jlat+1-j,1) = real(dble(work(i,j,1))*xscale+xadd)
          vvar(i,jlat+1-j,2) = real(dble(work(i,j,2))*xscale+xadd)
          vvar(i,jlat+1-j,3) = real(dble(work(i,j,3))*xscale+xadd)
          vvar(i,jlat+1-j,4) = real(dble(work(i,j,4))*xscale+xadd)
          vvar(i,jlat+1-j,5) = real(dble(work(i,j,5))*xscale+xadd)
          vvar(i,jlat+1-j,6) = real(dble(work(i,j,6))*xscale+xadd)
          vvar(i,jlat+1-j,7) = real(dble(work(i,j,7))*xscale+xadd)
          vvar(i,jlat+1-j,8) = real(dble(work(i,j,8))*xscale+xadd)
          vvar(i,jlat+1-j,9) = real(dble(work(i,j,9))*xscale+xadd)
          vvar(i,jlat+1-j,10) = real(dble(work(i,j,10))*xscale+xadd)
          vvar(i,jlat+1-j,11) = real(dble(work(i,j,11))*xscale+xadd)
          vvar(i,jlat+1-j,12) = real(dble(work(i,j,13))*xscale+xadd)
          vvar(i,jlat+1-j,13) = real(dble(work(i,j,15))*xscale+xadd)
          vvar(i,jlat+1-j,14) = real(dble(work(i,j,17))*xscale+xadd)
          vvar(i,jlat+1-j,15) = real(dble(work(i,j,18))*xscale+xadd)
          vvar(i,jlat+1-j,16) = real(dble(work(i,j,20))*xscale+xadd)
          vvar(i,jlat+1-j,17) = real(dble(work(i,j,22))*xscale+xadd)
          vvar(i,jlat+1-j,18) = real(dble(work(i,j,24))*xscale+xadd)
          vvar(i,jlat+1-j,19) = real(dble(work(i,j,26))*xscale+xadd)
          vvar(i,jlat+1-j,20) = real(dble(work(i,j,28))*xscale+xadd)
          vvar(i,jlat+1-j,21) = real(dble(work(i,j,31))*xscale+xadd)
          vvar(i,jlat+1-j,22) = real(dble(work(i,j,34))*xscale+xadd)
          vvar(i,jlat+1-j,23) = real(dble(work(i,j,37))*xscale+xadd)
        end do
      end do
    else
    end if
  end do

  call zeit_co('readein')

99001 format (i4,'/',a4,i4,'.00.nc')
99002 format (i4,'/',a4,i4,'.06.nc')
99003 format (i4,'/',a4,i4,'.12.nc')
99004 format (i4,'/',a4,i4,'.18.nc')
99005 format (i4,'/',a5,i4,'.00.nc')
99006 format (i4,'/',a5,i4,'.06.nc')
99007 format (i4,'/',a5,i4,'.12.nc')
99008 format (i4,'/',a5,i4,'.18.nc')
!
  end subroutine ein6hour
!
!-----------------------------------------------------------------------
!
  subroutine headerein(ires)
  implicit none
!
  integer , intent(in) :: ires
  integer :: i , j , k , kr

  klev = 23
  inlev = 37
  select case (ires)
    case (15)
      jlat = 121
      ilon = 240
      xres = 1.50
      imindat = 1989010100
      imaxdat = 2010033118
    case (25)
      jlat = 73
      ilon = 144
      xres = 2.50
      imindat = 1989010100
      imaxdat = 1998123118
    case (75)
      jlat = 241
      ilon = 480
      xres = 0.750
      imindat = 1989010100
      imaxdat = 2007123118
    case default
      call die('headerein','Unsupported resolution',1)
  end select
!
!     Allocate working space
!
  call getmem3d(b2,1,ilon,1,jlat,1,klev*3,'mod_ein:b2')
  call getmem3d(d2,1,ilon,1,jlat,1,klev*2,'mod_ein:d2')
  call getmem1d(glat,1,jlat,'mod_ein:glat')
  call getmem1d(glon,1,ilon,'mod_ein:glon')
  call getmem1d(sigma1,1,klev,'mod_ein:sigma1')
  call getmem1d(sigmar,1,klev,'mod_ein:sigmar')
  call getmem3d(b3,1,jx,1,iy,1,klev*3,'mod_ein:b3')
  call getmem3d(d3,1,jx,1,iy,1,klev*2,'mod_ein:d3')
  call getmem3d(work,1,ilon,1,jlat,1,inlev,'mod_ein:work')

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
    glon(i) = float(i-1)*xres
  end do
  do j = 1 , jlat
    glat(j) = -90.0 + float(j-1)*xres
  end do
!
!     CHANGE ORDER OF VERTICAL INDEXES FOR PRESSURE LEVELS
!
  do k = 1 , klev
    kr = klev - k + 1
    sigma1(k) = sigmar(kr)
  end do
 
!     Set up pointers

  u3 => d3(:,:,1:klev)
  v3 => d3(:,:,klev+1:2*klev)
  t3 => b3(:,:,1:klev)
  h3 => b3(:,:,klev+1:2*klev)
  q3 => b3(:,:,2*klev+1:3*klev)
  uvar => d2(:,:,1:klev)
  vvar => d2(:,:,klev+1:2*klev)
  tvar => b2(:,:,1:klev)
  hvar => b2(:,:,klev+1:2*klev)
  rhvar => b2(:,:,2*klev+1:3*klev)

  end subroutine headerein

end module mod_ein

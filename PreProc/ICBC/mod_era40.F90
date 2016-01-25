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

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_stdio
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
  use mod_message
  use mod_nchelper
  use netcdf

  private

  integer(ik4) , parameter :: klev = 23 , jlat = 73 , ilon = 144

  real(rk8) , target , dimension(ilon,jlat,klev*3) :: b2
  real(rk8) , target , dimension(ilon,jlat,klev*2) :: d2
  real(rk8) , target , dimension(ilon,jlat,4*3+1) :: s2
  real(rk8) , pointer , dimension(:,:,:) :: b3
  real(rk8) , pointer , dimension(:,:,:) :: d3
  real(rk8) , pointer , dimension(:,:,:) :: s3

  real(rk8) , dimension(ilon,jlat,klev) :: wvar

  real(rk8) , pointer :: u3(:,:,:) , v3(:,:,:)
  real(rk8) , pointer :: h3(:,:,:) , q3(:,:,:) , t3(:,:,:)
  real(rk8) , pointer :: uvar(:,:,:) , vvar(:,:,:)
  real(rk8) , pointer :: hvar(:,:,:) , rhvar(:,:,:) , tvar(:,:,:)
  real(rk8) , pointer , dimension(:,:,:) :: qsoil , tsice , tsoil
  real(rk8) , pointer , dimension(:,:) :: snw
  real(rk8) , pointer , dimension(:,:,:) :: qs3 , ti3 , ts3
  real(rk8) , pointer , dimension(:,:) :: snow

  real(rk8) , dimension(jlat) :: glat
  real(rk8) , dimension(ilon) :: glon
  real(rk8) , dimension(klev) :: sigma1 , sigmar

  public :: getera40 , headerera

  contains

  subroutine getera40(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate

    call era6hour(dattyp,idate,globidate1)
    write (stdout,*) 'READ IN fields at DATE:' , tochar(idate)
    !
    ! Horizontal interpolation of both the scalar and vector fields
    !
    call bilinx2(b3,b2,xlon,xlat,glon,glat,ilon,jlat,jx,iy,klev*3)
    call bilinx2(d3,d2,dlon,dlat,glon,glat,ilon,jlat,jx,iy,klev*2)
    !
    ! Rotate U-V fields after horizontal interpolation
    !
    call uvrot4(u3,v3,dlon,dlat,clon,clat,xcone,jx,iy,klev,plon,plat,iproj)
    !
    ! Vertical interpolation
    !
    call top2btm(t3,jx,iy,klev)
    call top2btm(q3,jx,iy,klev)
    call top2btm(h3,jx,iy,klev)
    call top2btm(u3,jx,iy,klev)
    call top2btm(v3,jx,iy,klev)
    !
    ! New calculation of P* on rcm topography.
    call intgtb(pa,za,tlayer,topogm,t3,h3,sigmar,jx,iy,klev)

    call intpsn(ps4,topogm,pa,za,tlayer,ptop,jx,iy)
    call crs2dot(pd4,ps4,jx,iy,i_band)
    !
    ! Determine surface temps on rcm topography.
    ! Interpolation from pressure levels
    call intv3(ts4,t3,ps4,sigmar,ptop,jx,iy,klev)

    call readsst(ts4,idate)

    ! Interpolate U, V, T, and Q.
    call intv1(u4,u3,pd4,sigmah,sigmar,ptop,jx,iy,kz,klev)
    call intv1(v4,v3,pd4,sigmah,sigmar,ptop,jx,iy,kz,klev)
    call intv2(t4,t3,ps4,sigmah,sigmar,ptop,jx,iy,kz,klev)
    call intv1(q4,q3,ps4,sigmah,sigmar,ptop,jx,iy,kz,klev)
    call humid2(t4,q4,ps4,ptop,sigmah,jx,iy,kz)
  end subroutine getera40

  subroutine era6hour(dattyp,idate,idate0)
    implicit none

    character(len=5) , intent(in) :: dattyp
    type(rcm_time_and_date) , intent(in) :: idate , idate0

    integer(ik4) :: i , inet , ivar , it , j , k , k4 , kkrec , istatus
    integer(ik4) :: year , month , day , hour
    character(len=24) :: inname
    character(len=256) :: pathaddname
    integer(2) , dimension(ilon,jlat,klev) :: work
    real(rk8) :: xadd , xscale

    integer(ik4) , dimension(10) , save :: icount , istart
    real(rk8) , dimension(5,4) , save :: xoff , xscl
    integer(ik4) , dimension(5,4) , save :: inet6
    integer(ik4) , dimension(5,4) , save :: ivar6
    character(len=5) , dimension(5) :: varname
    data varname/'t' , 'z' , 'r' , 'u' , 'v'/

    call split_idate(idate,year,month,day,hour)

    if ( idate < 1957090100 .or. idate > 2002083118 ) then
      call die('getera40', 'ERA40 datasets is just available from'// &
               ' 1957090100 to 2002083118', 1)
    end if

    if ( idate == idate0 .or. (lfdoyear(idate) .and. lmidnight(idate))) then
      do k4 = 1 , 4
        do kkrec = 1 , 5
          if ( kkrec == 1 ) then
            if ( k4 == 1 ) then
              write (inname,99001) year , 'air.' , year
            else if ( k4 == 2 ) then
              write (inname,99002) year , 'air.' , year
            else if ( k4 == 3 ) then
              write (inname,99003) year , 'air.' , year
            else if ( k4 == 4 ) then
              write (inname,99004) year , 'air.' , year
            end if
          else if ( kkrec == 2 ) then
            if ( k4 == 1 ) then
              write (inname,99001) year , 'hgt.' , year
            else if ( k4 == 2 ) then
              write (inname,99002) year , 'hgt.' , year
            else if ( k4 == 3 ) then
              write (inname,99003) year , 'hgt.' , year
            else if ( k4 == 4 ) then
              write (inname,99004) year , 'hgt.' , year
            end if
          else if ( kkrec == 3 ) then
            if ( k4 == 1 ) then
              write (inname,99005) year , 'rhum.' , year
            else if ( k4 == 2 ) then
              write (inname,99006) year , 'rhum.' , year
            else if ( k4 == 3 ) then
              write (inname,99007) year , 'rhum.' , year
            else if ( k4 == 4 ) then
              write (inname,99008) year , 'rhum.' , year
            end if
          else if ( kkrec == 4 ) then
            if ( k4 == 1 ) then
              write (inname,99005) year , 'uwnd.' , year
            else if ( k4 == 2 ) then
              write (inname,99006) year , 'uwnd.' , year
            else if ( k4 == 3 ) then
              write (inname,99007) year , 'uwnd.' , year
            else if ( k4 == 4 ) then
              write (inname,99008) year , 'uwnd.' , year
            end if
          else if ( kkrec == 5 ) then
            if ( k4 == 1 ) then
              write (inname,99005) year , 'vwnd.' , year
            else if ( k4 == 2 ) then
              write (inname,99006) year , 'vwnd.' , year
            else if ( k4 == 3 ) then
              write (inname,99007) year , 'vwnd.' , year
            else if ( k4 == 4 ) then
              write (inname,99008) year , 'vwnd.' , year
            end if
          end if

          pathaddname = trim(inpglob)//pthsep//dattyp//pthsep//inname
          istatus = nf90_open(pathaddname,nf90_nowrite,inet6(kkrec,k4))
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error open '//trim(pathaddname))
          istatus = nf90_inq_varid(inet6(kkrec,k4), &
                                   varname(kkrec),ivar6(kkrec,k4))
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error find var '//varname(kkrec))
          istatus = nf90_get_att(inet6(kkrec,k4),ivar6(kkrec,k4), &
                                 'scale_factor',xscl(kkrec,k4))
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error find att scale_factor')
          istatus = nf90_get_att(inet6(kkrec,k4),ivar6(kkrec,k4), &
                                 'add_offset',xoff(kkrec,k4))
          call checkncerr(istatus,__FILE__,__LINE__, &
                          'Error find att add_offset')
          write (stdout,*) inet6(kkrec,k4) , trim(pathaddname) ,  &
                      xscl(kkrec,k4) , xoff(kkrec,k4)
        end do
      end do

    end if

    k4 = hour/6 + 1
    it = day
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
    if ( mod(year,4) == 0 .and. month > 2 ) it = it + 1
    if ( mod(year,100) == 0 .and. month > 2 ) it = it - 1
    if ( mod(year,400) == 0 .and. month > 2 ) it = it + 1
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
    if ( mod(year,4) == 0 ) icount(4) = 366
    if ( mod(year,100) == 0 ) icount(4) = 365
    if ( mod(year,400) == 0 ) icount(4) = 366
    if ( year == 2002 ) icount(4) = 243
    if ( year == 1957 ) icount(4) = 122
    if ( year == 1957 ) it = it - 243
    istart(4) = it
    icount(4) = 1

    do kkrec = 1 , 5
      inet = inet6(kkrec,k4)
      ivar = ivar6(kkrec,k4)
      istatus = nf90_get_var(inet,ivar,work,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
                      'Error read var '//varname(kkrec))
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

  end subroutine era6hour

  subroutine headerera
    implicit none
    integer(ik4) :: i , j , k , kr

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
    ! Initial global grid-point longitude & latitude
    !
    do i = 1 , ilon
      glon(i) = float(i-1)*2.5
    end do
    do j = 1 , jlat
      glat(j) = -90.0 + float(j-1)*2.5
    end do
    !
    ! Change order of vertical indexes for pressure levels
    !
    do k = 1 , klev
      kr = klev - k + 1
      sigma1(k) = sigmar(kr)
    end do

    call getmem3d(b3,1,jx,1,iy,1,klev*3,'mod_era40:b3')
    call getmem3d(d3,1,jx,1,iy,1,klev*2,'mod_era40:d3')
    call getmem3d(s3,1,jx,1,iy,1,4*3+1,'mod_era40:s3')

    ! Set up pointers

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

end module mod_era40

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

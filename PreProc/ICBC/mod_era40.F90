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
  use mod_kdinterp
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

  real(rkx) , target , dimension(ilon,jlat,klev*3) :: b2
  real(rkx) , target , dimension(ilon,jlat,klev*2) :: d2
  real(rkx) , pointer , dimension(:,:,:) :: b3
  real(rkx) , pointer , dimension(:,:,:) :: d3
  real(rkx) , pointer , dimension(:,:,:) :: d3u , d3v

  real(rkx) , dimension(ilon,jlat,klev) :: wvar

  real(rkx) , pointer :: u3(:,:,:) , v3(:,:,:)
  real(rkx) , pointer :: u3v(:,:,:) , v3u(:,:,:)
  real(rkx) , pointer :: h3(:,:,:) , q3(:,:,:) , t3(:,:,:)
  real(rkx) , pointer :: uvar(:,:,:) , vvar(:,:,:)
  real(rkx) , pointer :: hvar(:,:,:) , rhvar(:,:,:) , tvar(:,:,:)

  real(rkx) , dimension(jlat) :: glat
  real(rkx) , dimension(ilon) :: glon
  real(rkx) , dimension(klev) :: sigma1 , sigmar
  real(rkx) , parameter :: pss = 100.0_rkx

  type(h_interpolator) :: cross_hint , udot_hint , vdot_hint

  public :: get_era40 , init_era40 , conclude_era40

  contains

  subroutine init_era40
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
    if ( idynamic == 3 ) then
      call getmem3d(d3u,1,jx,1,iy,1,klev*2,'mod_era40:d3u')
      call getmem3d(d3v,1,jx,1,iy,1,klev*2,'mod_era40:d3v')
    else
      call getmem3d(d3,1,jx,1,iy,1,klev*2,'mod_era40:d3')
    end if

    call h_interpolator_create(cross_hint,glat,glon,xlat,xlon)
    if ( idynamic == 3 ) then
      call h_interpolator_create(udot_hint,glat,glon,ulat,ulon)
      call h_interpolator_create(vdot_hint,glat,glon,vlat,vlon)
    else
      call h_interpolator_create(udot_hint,glat,glon,dlat,dlon)
    end if

    ! Set up pointers
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
    rhvar => b2(:,:,2*klev+1:3*klev)

  end subroutine init_era40

  subroutine get_era40(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate

    call era6hour(dattyp,idate,globidate1)
    write (stdout,*) 'READ IN fields at DATE:' , tochar(idate)
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
    ! Rotate U-V fields after horizontal interpolation
    !
    if ( idynamic == 3 ) then
      call uvrot4(u3,v3u,ulon,ulat,clon,clat,xcone,jx,iy,klev,plon,plat,iproj)
      call uvrot4(u3v,v3,vlon,vlat,clon,clat,xcone,jx,iy,klev,plon,plat,iproj)
    else
      call uvrot4(u3,v3,dlon,dlat,clon,clat,xcone,jx,iy,klev,plon,plat,iproj)
    end if
    !
    ! Vertical interpolation
    !
!$OMP SECTIONS
!$OMP SECTION
    call top2btm(t3,jx,iy,klev)
!$OMP SECTION
    call top2btm(q3,jx,iy,klev)
!$OMP SECTION
    call top2btm(h3,jx,iy,klev)
!$OMP SECTION
    call top2btm(u3,jx,iy,klev)
!$OMP SECTION
    call top2btm(v3,jx,iy,klev)
!$OMP END SECTIONS
    !
    ! New calculation of P* on rcm topography.
    !
    if ( idynamic == 3 ) then
      call ucrs2dot(zud4,z0,jx,iy,kz,i_band)
      call vcrs2dot(zvd4,z0,jx,iy,kz,i_crm)
      call intzps(ps4,topogm,t3,h3,pss,sigmar,jx,iy,klev)
      call intz3(ts4,t3,h3,topogm,jx,iy,klev,0.6_rkx,0.85_rkx,0.5_rkx)
    else
      call intgtb(pa,za,tlayer,topogm,t3,h3,pss,sigmar,jx,iy,klev)
      call intpsn(ps4,topogm,pa,za,tlayer,ptop,jx,iy)
      call crs2dot(pd4,ps4,jx,iy,i_band,i_crm)
      call intv3(ts4,t3,ps4,pss,sigmar,ptop,jx,iy,klev)
    end if

    call readsst(ts4,idate)
    !
    ! Interpolate U, V, T, and Q.
    !
    if ( idynamic == 3 ) then
!$OMP SECTIONS
!$OMP SECTION
      call intz1(u4,u3,zud4,h3,topogm,jx,iy,kz,klev,0.6_rkx,0.2_rkx,0.2_rkx)
!$OMP SECTION
      call intz1(v4,v3,zvd4,h3,topogm,jx,iy,kz,klev,0.6_rkx,0.2_rkx,0.2_rkx)
!$OMP SECTION
      call intz1(t4,t3,z0,h3,topogm,jx,iy,kz,klev,0.6_rkx,0.85_rkx,0.5_rkx)
!$OMP SECTION
      call intz1(q4,q3,z0,h3,topogm,jx,iy,kz,klev,0.7_rkx,0.7_rkx,0.4_rkx)
!$OMP END SECTIONS
    else
!$OMP SECTIONS
!$OMP SECTION
      call intv1(u4,u3,pd4,sigmah,pss,sigmar,ptop,jx,iy,kz,klev,1)
!$OMP SECTION
      call intv1(v4,v3,pd4,sigmah,pss,sigmar,ptop,jx,iy,kz,klev,1)
!$OMP SECTION
      call intv2(t4,t3,ps4,sigmah,pss,sigmar,ptop,jx,iy,kz,klev)
!$OMP SECTION
      call intv1(q4,q3,ps4,sigmah,pss,sigmar,ptop,jx,iy,kz,klev,1)
!$OMP END SECTIONS
    end if
    call rh2mxr(t4,q4,ps4,ptop,sigmah,jx,iy,kz)
  end subroutine get_era40

  subroutine era6hour(dattyp,idate,idate0)
    implicit none

    character(len=5) , intent(in) :: dattyp
    type(rcm_time_and_date) , intent(in) :: idate , idate0

    integer(ik4) :: i , inet , ivar , it , j , k , k4 , kkrec , istatus
    integer(ik4) :: year , month , day , hour
    character(len=24) :: inname
    character(len=256) :: pathaddname
    integer(2) , dimension(ilon,jlat,klev) :: work
    real(rkx) :: xadd , xscale

    integer(ik4) , dimension(10) , save :: icount , istart
    real(rkx) , dimension(5,4) , save :: xoff , xscl
    integer(ik4) , dimension(5,4) , save :: inet6
    integer(ik4) , dimension(5,4) , save :: ivar6
    character(len=5) , dimension(5) :: varname
    character (len=*) , parameter :: f99001 = '(i4,"/",a4,i4,".00.nc")'
    character (len=*) , parameter :: f99002 = '(i4,"/",a4,i4,".06.nc")'
    character (len=*) , parameter :: f99003 = '(i4,"/",a4,i4,".12.nc")'
    character (len=*) , parameter :: f99004 = '(i4,"/",a4,i4,".18.nc")'
    character (len=*) , parameter :: f99005 = '(i4,"/",a5,i4,".00.nc")'
    character (len=*) , parameter :: f99006 = '(i4,"/",a5,i4,".06.nc")'
    character (len=*) , parameter :: f99007 = '(i4,"/",a5,i4,".12.nc")'
    character (len=*) , parameter :: f99008 = '(i4,"/",a5,i4,".18.nc")'
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
              write (inname,f99001) year , 'air.' , year
            else if ( k4 == 2 ) then
              write (inname,f99002) year , 'air.' , year
            else if ( k4 == 3 ) then
              write (inname,f99003) year , 'air.' , year
            else if ( k4 == 4 ) then
              write (inname,f99004) year , 'air.' , year
            end if
          else if ( kkrec == 2 ) then
            if ( k4 == 1 ) then
              write (inname,f99001) year , 'hgt.' , year
            else if ( k4 == 2 ) then
              write (inname,f99002) year , 'hgt.' , year
            else if ( k4 == 3 ) then
              write (inname,f99003) year , 'hgt.' , year
            else if ( k4 == 4 ) then
              write (inname,f99004) year , 'hgt.' , year
            end if
          else if ( kkrec == 3 ) then
            if ( k4 == 1 ) then
              write (inname,f99005) year , 'rhum.' , year
            else if ( k4 == 2 ) then
              write (inname,f99006) year , 'rhum.' , year
            else if ( k4 == 3 ) then
              write (inname,f99007) year , 'rhum.' , year
            else if ( k4 == 4 ) then
              write (inname,f99008) year , 'rhum.' , year
            end if
          else if ( kkrec == 4 ) then
            if ( k4 == 1 ) then
              write (inname,f99005) year , 'uwnd.' , year
            else if ( k4 == 2 ) then
              write (inname,f99006) year , 'uwnd.' , year
            else if ( k4 == 3 ) then
              write (inname,f99007) year , 'uwnd.' , year
            else if ( k4 == 4 ) then
              write (inname,f99008) year , 'uwnd.' , year
            end if
          else if ( kkrec == 5 ) then
            if ( k4 == 1 ) then
              write (inname,f99005) year , 'vwnd.' , year
            else if ( k4 == 2 ) then
              write (inname,f99006) year , 'vwnd.' , year
            else if ( k4 == 3 ) then
              write (inname,f99007) year , 'vwnd.' , year
            else if ( k4 == 4 ) then
              write (inname,f99008) year , 'vwnd.' , year
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
  end subroutine era6hour

  subroutine conclude_era40
    implicit none
    call h_interpolator_destroy(cross_hint)
    call h_interpolator_destroy(udot_hint)
    if ( idynamic == 3 ) then
      call h_interpolator_destroy(vdot_hint)
    end if
  end subroutine conclude_era40

end module mod_era40

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

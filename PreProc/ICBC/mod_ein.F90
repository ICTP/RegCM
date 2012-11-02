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

  integer(ik4) , parameter :: klev = 37

  real(rk8) , dimension(klev) :: plevs

  integer(ik4) :: jlat , ilon

  real(rk8) , pointer , dimension(:,:,:) :: b3
  real(rk8) , pointer , dimension(:,:,:) :: d3

  real(rk8) , pointer :: u3(:,:,:) , v3(:,:,:)
  real(rk8) , pointer :: h3(:,:,:) , q3(:,:,:) , t3(:,:,:)
  real(rk8) , pointer :: uvar(:,:,:) , vvar(:,:,:)
  real(rk8) , pointer :: hvar(:,:,:) , rhvar(:,:,:) , tvar(:,:,:)

  real(rk8) :: xres
  real(rk8) , pointer , dimension(:,:,:) :: b2
  real(rk8) , pointer , dimension(:,:,:) :: d2
  real(rk8) , pointer , dimension(:) :: glat
  real(rk8) , pointer , dimension(:) :: glon
  real(rk8) , pointer , dimension(:) :: sigma1 , sigmar
  integer(2) , pointer , dimension(:,:,:) :: work

  integer(ik4) , dimension(5,4) :: inet5
  integer(ik4) , dimension(5,4) :: ivar5
  real(rk8) , dimension(5,4) :: xoff , xscl

  public :: getein , headerein

  data plevs / 1.0, 2.0, 3.0, 5.0, 7.0, 10.0, 20.0, 30.0, 50.0, 70.0,  &
               100.0, 125.0, 150.0, 175.0, 200.0, 225.0, 250.0, 300.0, &
               350.0, 400.0, 450.0, 500.0, 550.0, 600.0, 650.0, 700.0, &
               750.0, 775.0, 800.0, 825.0, 850.0, 875.0, 900.0, 925.0, &
               950.0, 975.0, 1000.0 /
  contains

  subroutine getein(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    !
    ! READ DATA AT IDATE
    !
    call ein6hour(dattyp,idate,globidate1)
    write (stdout,*) 'READ IN fields at DATE:' , tochar(idate)
    !
    ! HORIZONTAL INTERPOLATION OF BOTH THE SCALAR AND VECTOR FIELDS
    !
    call bilinx2(b3,b2,xlon,xlat,glon,glat,ilon,jlat,jx,iy,klev*3)
    call bilinx2(d3,d2,dlon,dlat,glon,glat,ilon,jlat,jx,iy,klev*2)
    !
    ! ROTATE U-V FIELDS AFTER HORIZONTAL INTERPOLATION
    !
    call uvrot4(u3,v3,dlon,dlat,clon,clat,xcone,jx,iy,klev,plon,plat,iproj)
    !
    ! X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
    !
    ! V E R T I C A L   I N T E R P O L A T I O N
    !
    ! X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
    !
    !HH:  CHANGE THE VERTICAL ORDER.
    call top2btm(t3,jx,iy,klev)
    call top2btm(q3,jx,iy,klev)
    call top2btm(h3,jx,iy,klev)
    call top2btm(u3,jx,iy,klev)
    call top2btm(v3,jx,iy,klev)
    !HH:OVER
    !
    ! NEW CALCULATION OF P* ON RCM TOPOGRAPHY.
    !
    call intgtb(pa,za,tlayer,topogm,t3,h3,sigmar,jx,iy,klev)
    call intpsn(ps4,topogm,pa,za,tlayer,ptop,jx,iy)
    if ( i_band == 1 ) then
      call p1p2_band(b3pd,ps4,jx,iy)
    else
      call p1p2(b3pd,ps4,jx,iy)
    end if
    ! 
    ! INTERPOLATION FROM PRESSURE LEVELS
    !
    call intv3(ts4,t3,ps4,sigmar,ptop,jx,iy,klev)

    call readsst(ts4,idate)
    ! 
    ! F3  INTERPOLATE U, V, T, AND Q.
    !
    call intv1(u4,u3,b3pd,sigma2,sigmar,ptop,jx,iy,kz,klev)
    call intv1(v4,v3,b3pd,sigma2,sigmar,ptop,jx,iy,kz,klev)
    call intv2(t4,t3,ps4,sigma2,sigmar,ptop,jx,iy,kz,klev)
    call intv1(q4,q3,ps4,sigma2,sigmar,ptop,jx,iy,kz,klev)
    ! 
    call humid2(t4,q4,ps4,ptop,sigma2,jx,iy,kz)
    !
    ! F4  DETERMINE H
    !
    call hydrost(h4,t4,topogm,ps4,ptop,sigmaf,sigma2,dsigma,jx,iy,kz)
!
  end subroutine getein
!
!-----------------------------------------------------------------------
!
  subroutine ein6hour(dattyp,idate,idate0)
    implicit none
!
    character(len=5) , intent(in) :: dattyp
    type(rcm_time_and_date) , intent(in) :: idate , idate0
!
    integer(ik4) :: i , inet , it , j , k , k4 , kkrec , istatus , ivar
    character(len=64) :: inname
    character(len=256) :: pathaddname
    character(len=1) , dimension(5) :: varname
    character(len=4) , dimension(5) :: fname
    character(len=4) , dimension(4) :: hname
    real(rk8) :: xadd , xscale
    integer(ik4) , dimension(4) :: icount , istart
    integer(ik4) :: year , month , day , hour , monthp1
    integer(ik4) , save :: lastmonth , lastyear
    !
    ! This is the latitude, longitude dimension of the grid to be read.
    ! This corresponds to the lat and lon dimension variables in the
    ! netCDF file.
    ! The data are packed into short integers (INTEGER*2).  The array
    ! work will be used to hold the packed integers.  The array 'x'
    ! will contain the unpacked data.
    !
    data varname /'t' , 'z' , 'r' , 'u' , 'v'/
    data fname   /'air','hgt','rhum','uwnd','vwnd'/
    data hname   /'.00.','.06.','.12.','.18.'/
!
    k4 = 1
    call split_idate(idate,year,month,day,hour)

    if ( dattyp == 'EIXXX' ) then
      if ( idate == idate0 .or. month /= lastmonth ) then
        lastmonth = month
        do kkrec = 1 , 5
          monthp1 = month+1
          if ( monthp1 == 13 ) monthp1 = 1
          write(inname,'(a,i0.2,a,i0.2,a)') &
             varname(kkrec)//'_xxxx',month,'0100-xxxx',monthp1,'0100.nc'
          pathaddname = trim(inpglob)//pthsep//'ERAIN_MEAN'//&
                        pthsep//'XXXX'//pthsep//inname
          istatus = nf90_open(pathaddname,nf90_nowrite,inet5(kkrec,1))
          call checkncerr(istatus,__FILE__,__LINE__, &
            'Error open file '//trim(pathaddname))
          istatus = nf90_inq_varid(inet5(kkrec,1),varname(kkrec), &
                                   ivar5(kkrec,1))
          call checkncerr(istatus,__FILE__,__LINE__, &
            'Error find var '//varname(kkrec))
          istatus = nf90_get_att(inet5(kkrec,1),ivar5(kkrec,1), &
                   'scale_factor',xscl(kkrec,1))
          call checkncerr(istatus,__FILE__,__LINE__, &
            'Error find att scale_factor')
          istatus = nf90_get_att(inet5(kkrec,1),ivar5(kkrec,1),  &
                   'add_offset',xoff(kkrec,1))
          call checkncerr(istatus,__FILE__,__LINE__, &
            'Error find att add_offset')
          write (stdout,*) inet5(kkrec,1) , trim(pathaddname) ,   &
                           xscl(kkrec,1) , xoff(kkrec,1)
        end do
      end if
      it = (day-1)*4 + hour/6 + 1
    else
      if ( idate == idate0 .or. year /= lastyear ) then
        lastyear = year
        do k4 = 1 , 4
          do kkrec = 1 , 5
            write(inname,'(i4,a,a,i4,a)') &
              year, pthsep, trim(fname(kkrec))//'.', year, hname(k4)//'nc'
            pathaddname = trim(inpglob)//pthsep//dattyp//pthsep//inname
            istatus = nf90_open(pathaddname,nf90_nowrite,inet5(kkrec,k4))
            call checkncerr(istatus,__FILE__,__LINE__, &
              'Error open file '//trim(pathaddname))
            istatus = nf90_inq_varid(inet5(kkrec,k4),varname(kkrec), &
                                     ivar5(kkrec,k4))
            call checkncerr(istatus,__FILE__,__LINE__, &
              'Error find var '//varname(kkrec))
            istatus = nf90_get_att(inet5(kkrec,k4),ivar5(kkrec,k4), &
                     'scale_factor',xscl(kkrec,k4))
            call checkncerr(istatus,__FILE__,__LINE__, &
              'Error find att scale_factor')
            istatus = nf90_get_att(inet5(kkrec,k4),ivar5(kkrec,k4),  &
                     'add_offset',xoff(kkrec,k4))
            call checkncerr(istatus,__FILE__,__LINE__, &
              'Error find att add_offset')
            write (stdout,*) inet5(kkrec,k4) , trim(pathaddname) ,   &
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
    end if
 
    do k = 1 , 4
      istart(k) = 1
    end do
    icount(1) = ilon
    icount(2) = jlat
    icount(3) = klev
    istart(4) = it
    icount(4) = 1

    do kkrec = 1 , 5
      inet = inet5(kkrec,k4)
      ivar = ivar5(kkrec,k4)
      istatus = nf90_get_var(inet,ivar,work,istart,icount)
      call checkncerr(istatus,__FILE__,__LINE__, &
        'Error read var '//varname(kkrec))
      xscale = xscl(kkrec,k4)
      xadd = xoff(kkrec,k4)
      if ( kkrec == 1 ) then
        do j = 1 , jlat
          do i = 1 , ilon
            tvar(i,jlat+1-j,:) = real(dble(work(i,j,:))*xscale+xadd)
          end do
        end do
      else if ( kkrec == 2 ) then
        do j = 1 , jlat
          do i = 1 , ilon
            hvar(i,jlat+1-j,:) = real(dble(work(i,j,:))*xscale+xadd)/9.80616
          end do
        end do
      else if ( kkrec == 3 ) then
        do j = 1 , jlat
          do i = 1 , ilon
            rhvar(i,jlat+1-j,:) = &
                amax1(real(dble(work(i,j,:))*xscale+xadd)*0.01,0.0)
          end do
        end do
      else if ( kkrec == 4 ) then
        do j = 1 , jlat
          do i = 1 , ilon
            uvar(i,jlat+1-j,:) = real(dble(work(i,j,:))*xscale+xadd)
          end do
        end do
      else if ( kkrec == 5 ) then
        do j = 1 , jlat
          do i = 1 , ilon
            vvar(i,jlat+1-j,:) = real(dble(work(i,j,:))*xscale+xadd)
          end do
        end do
      end if
    end do
  end subroutine ein6hour
!
!-----------------------------------------------------------------------
!
  subroutine headerein(ires)
    implicit none
!
    integer(ik4) , intent(in) :: ires
    integer(ik4) :: i , j , k , kr

    select case (ires)
      case (15)
        jlat = 121
        ilon = 240
        xres = 1.50
      case (25)
        jlat = 73
        ilon = 144
        xres = 2.50
      case (75)
        jlat = 241
        ilon = 480
        xres = 0.750
      case default
        call die('headerein','Unsupported resolution',1)
    end select
    !
    ! Allocate working space
    !
    call getmem3d(b2,1,ilon,1,jlat,1,klev*3,'mod_ein:b2')
    call getmem3d(d2,1,ilon,1,jlat,1,klev*2,'mod_ein:d2')
    call getmem1d(glat,1,jlat,'mod_ein:glat')
    call getmem1d(glon,1,ilon,'mod_ein:glon')
    call getmem1d(sigma1,1,klev,'mod_ein:sigma1')
    call getmem1d(sigmar,1,klev,'mod_ein:sigmar')
    call getmem3d(b3,1,jx,1,iy,1,klev*3,'mod_ein:b3')
    call getmem3d(d3,1,jx,1,iy,1,klev*2,'mod_ein:d3')
    call getmem3d(work,1,ilon,1,jlat,1,klev,'mod_ein:work')

    sigmar(:) = plevs(:)/1000.0
    !
    ! INITIAL GLOBAL GRID-POINT LONGITUDE & LATITUDE
    !
    do i = 1 , ilon
      glon(i) = float(i-1)*xres
    end do
    do j = 1 , jlat
      glat(j) = -90.0 + float(j-1)*xres
    end do
    !
    ! CHANGE ORDER OF VERTICAL INDEXES FOR PRESSURE LEVELS
    !
    do k = 1 , klev
      kr = klev - k + 1
      sigma1(k) = sigmar(kr)
    end do
    ! 
    ! Set up pointers
    !
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

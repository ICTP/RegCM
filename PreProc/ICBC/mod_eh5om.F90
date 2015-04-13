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

module mod_eh5om

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_stdio
  use mod_grid
  use mod_memutil
  use mod_write
  use mod_interp
  use mod_vertint
  use mod_hgt
  use mod_humid
  use mod_mksst
  use mod_uvrot
  use mod_vectutil
  use mod_message
  use netcdf

  private

  integer(ik4) , parameter :: klev = 17 , jlat = 96 , ilon = 192
  integer(ik4) , parameter :: mlev = 31

  real(rk8) , target , dimension(ilon,jlat,klev*3) :: b2
  real(rk8) , target , dimension(ilon,jlat,klev*2) :: d2
  real(rk8) , pointer , dimension(:,:,:) :: b3
  real(rk8) , pointer , dimension(:,:,:) :: d3

  real(rk8) , pointer :: u3(:,:,:) , v3(:,:,:)
  real(rk8) , pointer :: h3(:,:,:) , q3(:,:,:) , t3(:,:,:)
  real(rk8) , pointer :: uvar(:,:,:) , vvar(:,:,:)
  real(rk8) , pointer :: hvar(:,:,:) , rhvar(:,:,:) , tvar(:,:,:)

  real(rk8) , dimension(jlat) :: glat
  real(rk8) , dimension(ilon) :: glon
  real(rk8) , dimension(klev) :: sigma1 , sigmar

  real(rk8) , dimension(mlev+1) :: hyai , hybi
  real(rk8) , dimension(mlev) :: hyam , hybm

  integer(4) , dimension(10) :: icount , istart

  public :: geteh5om , headermpi

  contains

  subroutine geteh5om(idate)
    implicit none
    type(rcm_time_and_date) , intent(in) :: idate
    character(len=3) , dimension(12) :: chmon
    character(len=21) :: finm
    integer(ik4) :: i , i2 , ii , j , j2 , k , nrec , numx , numy
    integer(2) , dimension(ilon,jlat) :: itmp
    logical :: there
    real(rk8) :: offset , xscale
    character(len=4) , dimension(100) :: yr_a2
    character(len=4) , dimension(61) :: yr_rf
    character(len=4) :: namepart
    character(len=3) :: codepart
    integer(ik4) :: year , month , day , hour
    integer :: iehlen

    data yr_rf/'1941' , '1942' , '1943' , '1944' , '1945' , '1946' ,  &
               '1947' , '1948' , '1949' , '1950' , '1951' , '1952' ,  &
               '1953' , '1954' , '1955' , '1956' , '1957' , '1958' ,  &
               '1959' , '1960' , '1961' , '1962' , '1963' , '1964' ,  &
               '1965' , '1966' , '1967' , '1968' , '1969' , '1970' ,  &
               '1971' , '1972' , '1973' , '1974' , '1975' , '1976' ,  &
               '1977' , '1978' , '1979' , '1980' , '1981' , '1982' ,  &
               '1983' , '1984' , '1985' , '1986' , '1987' , '1988' ,  &
               '1989' , '1990' , '1991' , '1992' , '1993' , '1994' ,  &
               '1995' , '1996' , '1997' , '1998' , '1999' , '2000' ,  &
               '2001'/
    data yr_a2/'2001' , '2002' , '2003' , '2004' , '2005' , '2006' ,  &
               '2007' , '2008' , '2009' , '2010' , '2011' , '2012' ,  &
               '2013' , '2014' , '2015' , '2016' , '2017' , '2018' ,  &
               '2019' , '2020' , '2021' , '2022' , '2023' , '2024' ,  &
               '2025' , '2026' , '2027' , '2028' , '2029' , '2030' ,  &
               '2031' , '2032' , '2033' , '2034' , '2035' , '2036' ,  &
               '2037' , '2038' , '2039' , '2040' , '2041' , '2042' ,  &
               '2043' , '2044' , '2045' , '2046' , '2047' , '2048' ,  &
               '2049' , '2050' , '2051' , '2052' , '2053' , '2054' ,  &
               '2055' , '2056' , '2057' , '2058' , '2059' , '2060' ,  &
               '2061' , '2062' , '2063' , '2064' , '2065' , '2066' ,  &
               '2067' , '2068' , '2069' , '2070' , '2071' , '2072' ,  &
               '2073' , '2074' , '2075' , '2076' , '2077' , '2078' ,  &
               '2079' , '2080' , '2081' , '2082' , '2083' , '2084' ,  &
               '2085' , '2086' , '2087' , '2088' , '2089' , '2090' ,  &
               '2091' , '2092' , '2093' , '2094' , '2095' , '2096' ,  &
               '2097' , '2098' , '2099' , '2100'/
    data chmon/'JAN' , 'FEB' , 'MAR' , 'APR' , 'MAY' , 'JUN' , 'JUL' ,&
               'AUG' , 'SEP' , 'OCT' , 'NOV' , 'DEC'/

    call split_idate(idate,year,month,day,hour)

    if ( dattyp(1:2) == 'EH' ) then
      if ( idate < 1941010106 ) then
        call die('geteh5om','EH dataset is only available from 1941010106',1)
      end if
      if ( idate >= 2101010100 ) then
        call die('geteh5om','EH dataset is only available up to 2100123118',1)
      end if
    end if
    if (abs(lon1-lon0) < 1D-30 .and. abs(lat1-lat0) < 1D-30) then
      write (stdout, *) 'Assuming You have global dataset EH5OM'
      lon0 = 0.0D0
      lon1 = 358.125D0
      lat0 = -89.0625D0
      lat1 = 89.0625D0
    end if
    numx = nint((lon1-lon0)/1.875) + 1
    numy = nint((lat1-lat0)/1.875) + 1

    if ( dattyp == 'EHA1B' ) then
      if ( numx /= ilon .or. numy /= jlat ) then
        namepart = 'EH'
      else
        namepart = 'Eg'
      end if
      codepart = 'A1B'
    else
      if ( numx /= ilon .or. numy /= jlat ) then
        namepart = 'EH_'
      else
        namepart = 'EHg'
      end if
      codepart = dattyp(4:5)
    end if
    ! Any Normal month, inside the month
    if ( day /= 1 .or. hour /= 0 ) then
      if ( year < 2001 ) then
        finm = 'RF/'//yr_rf(year-1940)//'/'//trim(namepart)//'RF'//  &
                 yr_rf(year-1940)//chmon(month)
      else
        finm = trim(codepart)//'/'//yr_a2(year-2000)//'/'//&
               trim(namepart)//trim(codepart)//  &
               yr_a2(year-2000)//chmon(month)
      end if
    else if ( month /= 1 ) then
      ! Day 1 and hour 0 , but NOT January
      !   (First time is in previous month file)
      if ( year < 2001 ) then
        finm = 'RF/'//yr_rf(year-1940)//'/'//trim(namepart)//'RF'//    &
               yr_rf(year-1940)//chmon(month-1)
      else
        finm = trim(codepart)//'/'//yr_a2(year-2000)//'/'//&
               trim(namepart)//trim(codepart)//  &
               yr_a2(year-2000)//chmon(month-1)
      end if
    else ! First of january at 00:00:00
      if ( year < 2001 ) then
          finm = 'RF/'//yr_rf(year-1940)//'/'//trim(namepart)//'RF'//   &
                 yr_rf(year-1940)//chmon(12)
      else if ( year == 2001 ) then
        finm = 'RF/'//'2000'//'/'//'EHgRF'//'2000'//chmon(12)
      else
        finm = trim(codepart)//'/'//yr_a2(year-2000)//'/'//&
               trim(namepart)//trim(codepart)//  &
               yr_a2(year-2000)//chmon(12)
      end if
    end if
    do k = 1 , klev*3
      do j = 1 , jlat
        do i = 1 , ilon
          b2(i,j,k) = -9999.
        end do
      end do
    end do
    do k = 1 , klev*2
      do j = 1 , jlat
        do i = 1 , ilon
          d2(i,j,k) = -9999.
        end do
      end do
    end do

    inquire (file=trim(inpglob)//'/EH5OM/'//finm,exist=there)
    if ( .not.there ) then
      call die('geteh5om',trim(inpglob)//'/EH5OM/'//finm// &
               ' is not available',1)
    end if
    inquire(iolength=iehlen) offset , xscale , itmp
    open (63,file=trim(inpglob)//'/EH5OM/'//finm,form='unformatted',  &
          recl=iehlen,access='direct',action='read',status='old')
    if ( day /= 1 .or. hour /= 0 ) then
      nrec = ((day-1)*4+hour/6-1)*(klev*5)
    else if ( month == 1 .or. month == 2 .or. &
              month == 4 .or. month == 6 .or. &
              month == 8 .or. month == 9 .or. month == 11 ) then
      nrec = (mlev*4-1)*(klev*5)
    else if ( month == 5 .or. month == 7 .or. &
              month == 10 .or. month == 12 ) then
      nrec = (30*4-1)*(klev*5)
    else if ( mod(year,4) == 0 .and. year /= 2100 ) then
      nrec = (29*4-1)*(klev*5)
    else
      nrec = (28*4-1)*(klev*5)
    end if

    do k = klev , 1 , -1
      nrec = nrec + 1
      read (63,rec=nrec) offset , xscale , itmp
      do j = nint((lat0+.9375)/1.875) , nint((lat1+.9375)/1.875)
        do i = nint(lon0/1.875) , nint(lon1/1.875)
          ii = i + 1
          if ( ii <= 0 ) ii = ii + ilon
          if ( ii > ilon ) ii = ii - ilon
          i2 = i - nint(lon0/1.875) + 1
          j2 = j - nint((lat0+.9375)/1.875) + 1
          if ( numx == ilon .and. numy == jlat ) then
            hvar(ii,49-j,k) = real(dble(itmp(i2,j2))*xscale + offset)
          else
            hvar(ii,j+jlat/2,k) = real(dble(itmp(i2,j2))*xscale + offset)
          end if
        end do
      end do
    end do
    do k = klev , 1 , -1
      nrec = nrec + 1
      read (63,rec=nrec) offset , xscale , itmp
      do j = nint((lat0+.9375)/1.875) , nint((lat1+.9375)/1.875)
        do i = nint(lon0/1.875) , nint(lon1/1.875)
          ii = i + 1
          if ( ii <= 0 ) ii = ii + ilon
          if ( ii > ilon ) ii = ii - ilon
          i2 = i - nint(lon0/1.875) + 1
          j2 = j - nint((lat0+.9375)/1.875) + 1
          if ( numx == ilon .and. numy == jlat ) then
            rhvar(ii,49-j,k) = real(dmin1(dmax1(dble(itmp(i2,j2))*xscale+ &
                                      offset,0.D0),1.D0))
          else
            rhvar(ii,j+jlat/2,k) = real(dmin1(dmax1(dble(itmp(i2,j2))*xscale+ &
                                          offset,0.D0),1.D0))
          end if
        end do
      end do
    end do
    do k = klev , 1 , -1
      nrec = nrec + 1
      read (63,rec=nrec) offset , xscale , itmp
      do j = nint((lat0+.9375)/1.875) , nint((lat1+.9375)/1.875)
        do i = nint(lon0/1.875) , nint(lon1/1.875)
          ii = i + 1
          if ( ii <= 0 ) ii = ii + ilon
          if ( ii > ilon ) ii = ii - ilon
          i2 = i - nint(lon0/1.875) + 1
          j2 = j - nint((lat0+.9375)/1.875) + 1
          if ( numx == ilon .and. numy == jlat ) then
            tvar(ii,49-j,k) = real(dble(itmp(i2,j2))*xscale + offset)
          else
            tvar(ii,j+jlat/2,k) = real(dble(itmp(i2,j2))*xscale + offset)
          end if
        end do
      end do
    end do
    do k = klev , 1 , -1
      nrec = nrec + 1
      read (63,rec=nrec) offset , xscale , itmp
      do j = nint((lat0+.9375)/1.875) , nint((lat1+.9375)/1.875)
        do i = nint(lon0/1.875) , nint(lon1/1.875)
          ii = i + 1
          if ( ii <= 0 ) ii = ii + ilon
          if ( ii > ilon ) ii = ii - ilon
          i2 = i - nint(lon0/1.875) + 1
          j2 = j - nint((lat0+.9375)/1.875) + 1
          if ( numx == ilon .and. numy == jlat ) then
            uvar(ii,49-j,k) = real(dble(itmp(i2,j2))*xscale + offset)
          else
            uvar(ii,j+jlat/2,k) = real(dble(itmp(i2,j2))*xscale + offset)
          end if
        end do
      end do
    end do
    do k = klev , 1 , -1
      nrec = nrec + 1
      read (63,rec=nrec) offset , xscale , itmp
      do j = nint((lat0+.9375)/1.875) , nint((lat1+.9375)/1.875)
        do i = nint(lon0/1.875) , nint(lon1/1.875)
          ii = i + 1
          if ( ii <= 0 ) ii = ii + ilon
          if ( ii > ilon ) ii = ii - ilon
          i2 = i - nint(lon0/1.875) + 1
          j2 = j - nint((lat0+.9375)/1.875) + 1
          if ( numx == ilon .and. numy == jlat ) then
            vvar(ii,49-j,k) = real(dble(itmp(i2,j2))*xscale + offset)
          else
            vvar(ii,j+jlat/2,k) = real(dble(itmp(i2,j2))*xscale + offset)
          end if
        end do
      end do
    end do
    close (63)
    write (stdout,*) 'READ IN fields at DATE:' , tochar(idate)

    call bilinx2(b3,b2,xlon,xlat,glon,glat,ilon,jlat,jx,iy,klev*3)
    call bilinx2(d3,d2,dlon,dlat,glon,glat,ilon,jlat,jx,iy,klev*2)

    call uvrot4(u3,v3,dlon,dlat,clon,clat,xcone,jx,iy,klev,plon,plat,iproj)

    call top2btm(t3,jx,iy,klev)
    call top2btm(q3,jx,iy,klev)
    call top2btm(h3,jx,iy,klev)
    call top2btm(u3,jx,iy,klev)
    call top2btm(v3,jx,iy,klev)

    call intgtb(pa,za,tlayer,topogm,t3,h3,sigmar,jx,iy,klev)

    call intpsn(ps4,topogm,pa,za,tlayer,ptop,jx,iy)
    if(i_band == 1) then
       call p1p2_band(b3pd,ps4,jx,iy)
    else
       call p1p2(b3pd,ps4,jx,iy)
    endif

    call intv3(ts4,t3,ps4,sigmar,ptop,jx,iy,klev)
    call readsst(ts4,idate)

    call intv1(u4,u3,b3pd,sigmah,sigmar,ptop,jx,iy,kz,klev)
    call intv1(v4,v3,b3pd,sigmah,sigmar,ptop,jx,iy,kz,klev)

    call intv2(t4,t3,ps4,sigmah,sigmar,ptop,jx,iy,kz,klev)

    call intv1(q4,q3,ps4,sigmah,sigmar,ptop,jx,iy,kz,klev)
    call humid2(t4,q4,ps4,ptop,sigmah,jx,iy,kz)

    call hydrost(h4,t4,topogm,ps4,ptop,sigmah,jx,iy,kz)
  end subroutine geteh5om

  subroutine headermpi
    implicit none
    integer(ik4) :: i , k , kr

    glat(1) = -88.5719985961914
    glat(2) = -86.7229995727539
    glat(3) = -84.8619995117188
    glat(4) = -82.9990005493164
    glat(5) = -81.1350021362305
    glat(6) = -79.2710037231445
    glat(7) = -77.4059982299805
    glat(8) = -75.5410003662109
    glat(9) = -73.6760025024414
    glat(10) = -71.8109970092773
    glat(11) = -69.9459991455078
    glat(12) = -68.0810012817383
    glat(13) = -66.2160034179688
    glat(14) = -64.3509979248047
    glat(15) = -62.4860000610352
    glat(16) = -60.6199989318848
    glat(17) = -58.7550010681152
    glat(18) = -56.8899993896484
    glat(19) = -55.0250015258789
    glat(20) = -53.1599998474121
    glat(21) = -51.2939987182617
    glat(22) = -49.4290008544922
    glat(23) = -47.5639991760254
    glat(24) = -45.6990013122559
    glat(25) = -43.8330001831055
    glat(26) = -41.9679985046387
    glat(27) = -40.1030006408691
    glat(28) = -38.2379989624023
    glat(29) = -36.3720016479492
    glat(30) = -34.5069999694824
    glat(31) = -32.6419982910156
    glat(32) = -30.7770004272461
    glat(33) = -28.9109992980957
    glat(34) = -27.0459995269775
    glat(35) = -25.1809997558594
    glat(36) = -23.3159999847412
    glat(37) = -21.4500007629395
    glat(38) = -19.5849990844727
    glat(39) = -17.7199993133545
    glat(40) = -15.8549995422363
    glat(41) = -13.9890003204346
    glat(42) = -12.1239995956421
    glat(43) = -10.2589998245239
    glat(44) = -8.39400005340576
    glat(45) = -6.52799987792969
    glat(46) = -4.66300010681152
    glat(47) = -2.79800009727478
    glat(48) = -0.933000028133392
    glat(49) = 0.933000028133392
    glat(50) = 2.79800009727478
    glat(51) = 4.66300010681152
    glat(52) = 6.52799987792969
    glat(53) = 8.39400005340576
    glat(54) = 10.2589998245239
    glat(55) = 12.1239995956421
    glat(56) = 13.9890003204346
    glat(57) = 15.8549995422363
    glat(58) = 17.7199993133545
    glat(59) = 19.5849990844727
    glat(60) = 21.4500007629395
    glat(61) = 23.3159999847412
    glat(62) = 25.1809997558594
    glat(63) = 27.0459995269775
    glat(64) = 28.9109992980957
    glat(65) = 30.7770004272461
    glat(66) = 32.6419982910156
    glat(67) = 34.5069999694824
    glat(68) = 36.3720016479492
    glat(69) = 38.2379989624023
    glat(70) = 40.1030006408691
    glat(71) = 41.9679985046387
    glat(72) = 43.8330001831055
    glat(73) = 45.6990013122559
    glat(74) = 47.5639991760254
    glat(75) = 49.4290008544922
    glat(76) = 51.2939987182617
    glat(77) = 53.1599998474121
    glat(78) = 55.0250015258789
    glat(79) = 56.8899993896484
    glat(80) = 58.7550010681152
    glat(81) = 60.6199989318848
    glat(82) = 62.4860000610352
    glat(83) = 64.3509979248047
    glat(84) = 66.2160034179688
    glat(85) = 68.0810012817383
    glat(86) = 69.9459991455078
    glat(87) = 71.8109970092773
    glat(88) = 73.6760025024414
    glat(89) = 75.5410003662109
    glat(90) = 77.4059982299805
    glat(91) = 79.2710037231445
    glat(92) = 81.1350021362305
    glat(93) = 82.9990005493164
    glat(94) = 84.8619995117188
    glat(95) = 86.7229995727539
    glat(96) = 88.5719985961914

    sigmar(1) = .01
    sigmar(2) = .03
    sigmar(3) = .05
    sigmar(4) = .07
    sigmar(5) = .1
    sigmar(6) = .15
    sigmar(7) = .2
    sigmar(8) = .25
    sigmar(9) = .3
    sigmar(10) = .4
    sigmar(11) = .5
    sigmar(12) = .6
    sigmar(13) = .7
    sigmar(14) = .775
    sigmar(15) = .85
    sigmar(16) = .925
    sigmar(17) = 1.0

    do i = 1 , ilon
      glon(i) = float(i-1)*1.875
    end do

    do k = 1 , klev
      kr = klev - k + 1
      sigma1(k) = sigmar(kr)
    end do

    call getmem3d(b3,1,jx,1,iy,1,klev*3,'mod_eh5om:b3')
    call getmem3d(d3,1,jx,1,iy,1,klev*2,'mod_eh5om:d3')

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

  end subroutine headermpi

end module mod_eh5om
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

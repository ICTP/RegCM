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

      use mod_dynparam
      use m_realkinds
      use m_stdio
      use m_die
      use m_zeit
      use m_mall

      private

      integer , parameter :: klev = 17 , jlat = 96 , ilon = 192
      integer , parameter :: mlev = 31

      real(sp) , target , dimension(ilon,jlat,klev*3) :: b2
      real(sp) , target , dimension(ilon,jlat,klev*2) :: d2
      real(sp) , allocatable , target , dimension(:,:,:) :: b3
      real(sp) , allocatable , target , dimension(:,:,:) :: d3
      real(sp) , allocatable , dimension(:,:,:) :: sulfate3
      real(sp) , allocatable , dimension(:,:) :: pso4_3

      real(sp) , pointer :: u3(:,:,:) , v3(:,:,:)
      real(sp) , pointer :: h3(:,:,:) , q3(:,:,:) , t3(:,:,:)
      real(sp) , pointer :: uvar(:,:,:) , vvar(:,:,:)
      real(sp) , pointer :: hvar(:,:,:) , rhvar(:,:,:) , tvar(:,:,:)

      real(sp) , dimension(jlat) :: glat
      real(sp) , dimension(ilon) :: glon
      real(sp) , dimension(klev) :: sigma1 , sigmar

      real(sp) , dimension(mlev+1) :: hyai , hybi
      real(sp) , dimension(mlev) :: hyam , hybm
      real(sp) , dimension(ilon,jlat) :: pso4_0
      real(sp) , dimension(ilon,jlat,mlev,12) :: sulfate

      real(sp) , dimension(ilon,jlat) :: pso4_2
      real(sp) , dimension(ilon,jlat,mlev,2) :: sulfate1
      real(sp) , dimension(ilon,jlat,mlev) :: sulfate2

      integer(4) , dimension(10) :: icount , istart

      public :: geteh5om , headermpi , footermpi

      contains

      subroutine geteh5om(idate)
      use netcdf
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
! Dummy arguments
!
      integer :: idate
!
! Local variables
!
      character(3) , dimension(12) :: chmon
      character(21) :: finm , psnm
      character(256) :: fnso4
      integer :: i , i2 , ii , j , j2 , k , k0 , krec , l , month ,     &
               & nday , nhour , nrec , numx , numy , nyear
      integer(2) , dimension(ilon,jlat) :: itmp
      real(dp) :: offset , xscale
      real(sp) :: pmpi , pmpj , prcm
      logical :: there
      character(4) , dimension(100) :: yr_a2
      character(4) , dimension(61) :: yr_rf
      integer :: ncid , istatus

!
      data yr_rf/'1941' , '1942' , '1943' , '1944' , '1945' , '1946' ,  &
          &'1947' , '1948' , '1949' , '1950' , '1951' , '1952' ,        &
         & '1953' , '1954' , '1955' , '1956' , '1957' , '1958' ,        &
         & '1959' , '1960' , '1961' , '1962' , '1963' , '1964' ,        &
         & '1965' , '1966' , '1967' , '1968' , '1969' , '1970' ,        &
         & '1971' , '1972' , '1973' , '1974' , '1975' , '1976' ,        &
         & '1977' , '1978' , '1979' , '1980' , '1981' , '1982' ,        &
         & '1983' , '1984' , '1985' , '1986' , '1987' , '1988' ,        &
         & '1989' , '1990' , '1991' , '1992' , '1993' , '1994' ,        &
         & '1995' , '1996' , '1997' , '1998' , '1999' , '2000' , '2001'/
      data yr_a2/'2001' , '2002' , '2003' , '2004' , '2005' , '2006' ,  &
          &'2007' , '2008' , '2009' , '2010' , '2011' , '2012' ,        &
         & '2013' , '2014' , '2015' , '2016' , '2017' , '2018' ,        &
         & '2019' , '2020' , '2021' , '2022' , '2023' , '2024' ,        &
         & '2025' , '2026' , '2027' , '2028' , '2029' , '2030' ,        &
         & '2031' , '2032' , '2033' , '2034' , '2035' , '2036' ,        &
         & '2037' , '2038' , '2039' , '2040' , '2041' , '2042' ,        &
         & '2043' , '2044' , '2045' , '2046' , '2047' , '2048' ,        &
         & '2049' , '2050' , '2051' , '2052' , '2053' , '2054' ,        &
         & '2055' , '2056' , '2057' , '2058' , '2059' , '2060' ,        &
         & '2061' , '2062' , '2063' , '2064' , '2065' , '2066' ,        &
         & '2067' , '2068' , '2069' , '2070' , '2071' , '2072' ,        &
         & '2073' , '2074' , '2075' , '2076' , '2077' , '2078' ,        &
         & '2079' , '2080' , '2081' , '2082' , '2083' , '2084' ,        &
         & '2085' , '2086' , '2087' , '2088' , '2089' , '2090' ,        &
         & '2091' , '2092' , '2093' , '2094' , '2095' , '2096' ,        &
         & '2097' , '2098' , '2099' , '2100'/
      data chmon/'JAN' , 'FEB' , 'MAR' , 'APR' , 'MAY' , 'JUN' , 'JUL' ,&
          &'AUG' , 'SEP' , 'OCT' , 'NOV' , 'DEC'/
!
!     D      BEGIN LOOP OVER NTIMES
!
      call zeit_ci('get_eh5om')
      nyear = idate/1000000
                         !, SST2(JX,IY)
      month = idate/10000 - nyear*100
      nday = idate/100 - (idate/10000)*100
      nhour = mod(idate,100)
 
      if ( ssttyp=='EH5RF' ) then
        if ( idate<1941010106 ) then
          call die('geteh5om','EH5RF dataset is only available' // &
                   ' from 1941010106',1)
        end if
        if ( ehso4 ) then
          if ( idate<1950010100 ) then
            call die('geteh5om','EH5RF sulfate is only available' // &
                     ' from 1950010100',1)
          end if
        end if
        if ( idate>2001010100 ) then
          call die('geteh5om','EH5RF dataset is only available' // &
                   ' up to 2001010100',1)
        end if
      end if
      if ( ssttyp=='EH5A2' ) then
        if ( idate<2001010100 ) then
          call die('geteh5om','EH5A2 dataset is only available' // &
                   ' from 2001010100',1)
        end if
        if ( idate>=2101010100 ) then
          call die('geteh5om','EH5A2 dataset is only available' // &
                   ' up to 2100123118',1)
        end if
      end if
      if ( ssttyp=='EH5B1' ) then
        if ( idate<2001010100 ) then
          call die('geteh5om','EH5B1 dataset is only available' // &
                   ' from 2001010100',1)
        end if
        if ( idate>=2101010100 ) then
          call die('geteh5om','EH5B1 dataset is only available' // &
                   ' up to 2100123118',1)
        end if
      end if
      if ( ssttyp=='EHA1B' ) then
        if ( idate<2001010100 ) then
          call die('geteh5om','EHA1B dataset is only available' // &
                   ' from 2001010100',1)
        end if
        if ( idate>=2101010100 ) then
          call die('geteh5om','EHA1B dataset is only available' // &
                   ' up to 2100123118',1)
        end if
      end if
      if (abs(lon1-lon0)<1E-30 .and. abs(lat1-lat0)<1E-30) then
        write (stdout, *) 'Assuming You have global dataset EH5OM'
        lon0 = 0
        lon1 = 358.125
        lat0 = -89.0625
        lat1 = 89.0625
      end if
      numx = nint((lon1-lon0)/1.875) + 1
      numy = nint((lat1-lat0)/1.875) + 1
      if ( numx/=ilon .or. numy/=jlat ) then
        if ( nday/=1 .or. nhour/=0 ) then
          if ( ssttyp=='EH5RF' ) then
            finm = 'RF/'//yr_rf(nyear-1940)//'/'//'EH_RF'//             &
                 & yr_rf(nyear-1940)//chmon(month)
            if ( ehso4 ) psnm = 'RF/'//yr_rf(nyear-1940)//'/'//'EH_PS'//&
                              & yr_rf(nyear-1940)//chmon(month)
          else if ( ssttyp=='EH5A2' ) then
            finm = 'A2/'//yr_a2(nyear-2000)//'/'//'EH_A2'//             &
                 & yr_a2(nyear-2000)//chmon(month)
            if ( ehso4 ) psnm = 'A2/'//yr_a2(nyear-2000)//'/'//'EH_PS'//&
                              & yr_a2(nyear-2000)//chmon(month)
          else if ( ssttyp=='EH5B1' ) then
            finm = 'B1/'//yr_a2(nyear-2000)//'/'//'EH_B1'//             &
                 & yr_a2(nyear-2000)//chmon(month)
            if ( ehso4 ) psnm = 'B1/'//yr_a2(nyear-2000)//'/'//'EH_PS'//&
                              & yr_a2(nyear-2000)//chmon(month)
          else if ( ssttyp=='EHA1B' ) then
            finm = 'A1B/'//yr_a2(nyear-2000)//'/'//'E_A1B'//            &
                 & yr_a2(nyear-2000)//chmon(month)
            if ( ehso4 ) psnm = 'A1B/'//yr_a2(nyear-2000)               &
                              & //'/'//'EH_PS'//yr_a2(nyear-2000)       &
                              & //chmon(month)
          else
            call die('geteh5om','ERROR IN geteh5om',1)
          end if
        else if ( month/=1 ) then
          if ( ssttyp=='EH5RF' ) then
            finm = 'RF/'//yr_rf(nyear-1940)//'/'//'EH_RF'//             &
                 & yr_rf(nyear-1940)//chmon(month-1)
            if ( ehso4 ) psnm = 'RF/'//yr_rf(nyear-1940)//'/'//'EH_PS'//&
                              & yr_rf(nyear-1940)//chmon(month-1)
          else if ( ssttyp=='EH5A2' ) then
            finm = 'A2/'//yr_a2(nyear-2000)//'/'//'EH_A2'//             &
                 & yr_a2(nyear-2000)//chmon(month-1)
            if ( ehso4 ) psnm = 'A2/'//yr_a2(nyear-2000)//'/'//'EH_PS'//&
                              & yr_a2(nyear-2000)//chmon(month-1)
          else if ( ssttyp=='EH5B1' ) then
            finm = 'B1/'//yr_a2(nyear-2000)//'/'//'EH_B1'//             &
                 & yr_a2(nyear-2000)//chmon(month-1)
            if ( ehso4 ) psnm = 'B1/'//yr_a2(nyear-2000)//'/'//'EH_PS'//&
                              & yr_a2(nyear-2000)//chmon(month-1)
          else if ( ssttyp=='EHA1B' ) then
            finm = 'A1B/'//yr_a2(nyear-2000)//'/'//'E_A1B'//            &
                 & yr_a2(nyear-2000)//chmon(month-1)
            if ( ehso4 ) psnm = 'A1B/'//yr_a2(nyear-2000)               &
                              & //'/'//'EH_PS'//yr_a2(nyear-2000)       &
                              & //chmon(month-1)
          else
            call die('geteh5om','ERROR IN geteh5om',1)
          end if
        else if ( ssttyp=='EH5RF' ) then
          finm = 'RF/'//yr_rf(nyear-1941)//'/'//'EH_RF'//               &
               & yr_rf(nyear-1941)//chmon(12)
          if ( ehso4 ) psnm = 'RF/'//yr_rf(nyear-1941)//'/'//'EH_PS'//  &
                            & yr_rf(nyear-1941)//chmon(12)
        else if ( ssttyp=='EH5A2' ) then
          if ( nyear==2001 ) then
            finm = 'RF/'//yr_rf(nyear-1941)//'/'//'EH_RF'//             &
                 & yr_rf(nyear-1941)//chmon(12)
            if ( ehso4 ) psnm = 'RF/'//yr_rf(nyear-1941)//'/'//'EH_PS'//&
                              & yr_rf(nyear-1941)//chmon(12)
          else
            finm = 'A2/'//yr_a2(nyear-2001)//'/'//'EH_A2'//             &
                 & yr_a2(nyear-2001)//chmon(12)
            if ( ehso4 ) psnm = 'A2/'//yr_a2(nyear-2001)//'/'//'EH_PS'//&
                              & yr_a2(nyear-2001)//chmon(12)
          end if
        else if ( ssttyp=='EH5B1' ) then
          if ( nyear==2001 ) then
            finm = 'RF/'//yr_rf(nyear-1941)//'/'//'EH_RF'//             &
                 & yr_rf(nyear-1941)//chmon(12)
            if ( ehso4 ) psnm = 'RF/'//yr_rf(nyear-1941)//'/'//'EH_PS'//&
                              & yr_rf(nyear-1941)//chmon(12)
          else
            finm = 'B1/'//yr_a2(nyear-2001)//'/'//'EH_B1'//             &
                 & yr_a2(nyear-2001)//chmon(12)
            if ( ehso4 ) psnm = 'B1/'//yr_a2(nyear-2001)//'/'//'EH_PS'//&
                              & yr_a2(nyear-2001)//chmon(12)
          end if
        else if ( ssttyp=='EHA1B' ) then
          if ( nyear==2001 ) then
            finm = 'RF/'//yr_rf(nyear-1941)//'/'//'EH_RF'//             &
                 & yr_rf(nyear-1941)//chmon(12)
            if ( ehso4 ) psnm = 'RF/'//yr_rf(nyear-1941)//'/'//'EH_PS'//&
                              & yr_rf(nyear-1941)//chmon(12)
          else
            finm = 'A1B/'//yr_a2(nyear-2001)//'/'//'E_A1B'//            &
                 & yr_a2(nyear-2001)//chmon(12)
            if ( ehso4 ) psnm = 'A1B/'//yr_a2(nyear-2001)               &
                              & //'/'//'EH_PS'//yr_a2(nyear-2001)       &
                              & //chmon(12)
          end if
        else
          call die('geteh5om','ERROR IN geteh5om',1)
        end if
      else if ( nday/=1 .or. nhour/=0 ) then
        if ( ssttyp=='EH5RF' ) then
          finm = 'RF/'//yr_rf(nyear-1940)//'/'//'EHgRF'//               &
               & yr_rf(nyear-1940)//chmon(month)
          if ( ehso4 ) psnm = 'RF/'//yr_rf(nyear-1940)//'/'//'EHgPS'//  &
                            & yr_rf(nyear-1940)//chmon(month)
        else if ( ssttyp=='EH5A2' ) then
          finm = 'A2/'//yr_a2(nyear-2000)//'/'//'EHgA2'//               &
               & yr_a2(nyear-2000)//chmon(month)
          if ( ehso4 ) psnm = 'A2/'//yr_a2(nyear-2000)//'/'//'EHgPS'//  &
                            & yr_a2(nyear-2000)//chmon(month)
        else if ( ssttyp=='EH5B1' ) then
          finm = 'B1/'//yr_a2(nyear-2000)//'/'//'EHgB1'//               &
               & yr_a2(nyear-2000)//chmon(month)
          if ( ehso4 ) psnm = 'B1/'//yr_a2(nyear-2000)//'/'//'EHgPS'//  &
                            & yr_a2(nyear-2000)//chmon(month)
        else if ( ssttyp=='EHA1B' ) then
          finm = 'A1B/'//yr_a2(nyear-2000)//'/'//'EgA1B'//              &
               & yr_a2(nyear-2000)//chmon(month)
          if ( ehso4 ) psnm = 'A1B/'//yr_a2(nyear-2000)//'/'//'EHgPS'// &
                            & yr_a2(nyear-2000)//chmon(month)
        else
          call die('geteh5om','ERROR IN geteh5om',1)
        end if
      else if ( month/=1 ) then
        if ( ssttyp=='EH5RF' ) then
          finm = 'RF/'//yr_rf(nyear-1940)//'/'//'EHgRF'//               &
               & yr_rf(nyear-1940)//chmon(month-1)
          if ( ehso4 ) psnm = 'RF/'//yr_rf(nyear-1940)//'/'//'EHgPS'//  &
                            & yr_rf(nyear-1940)//chmon(month-1)
        else if ( ssttyp=='EH5A2' ) then
          finm = 'A2/'//yr_a2(nyear-2000)//'/'//'EHgA2'//               &
               & yr_a2(nyear-2000)//chmon(month-1)
          if ( ehso4 ) psnm = 'A2/'//yr_a2(nyear-2000)//'/'//'EHgPS'//  &
                            & yr_a2(nyear-2000)//chmon(month-1)
        else if ( ssttyp=='EH5B1' ) then
          finm = 'B1/'//yr_a2(nyear-2000)//'/'//'EHgB1'//               &
               & yr_a2(nyear-2000)//chmon(month-1)
          if ( ehso4 ) psnm = 'B1/'//yr_a2(nyear-2000)//'/'//'EHgPS'//  &
                            & yr_a2(nyear-2000)//chmon(month-1)
        else if ( ssttyp=='EHA1B' ) then
          finm = 'A1B/'//yr_a2(nyear-2000)//'/'//'EgA1B'//              &
               & yr_a2(nyear-2000)//chmon(month-1)
          if ( ehso4 ) psnm = 'A1B/'//yr_a2(nyear-2000)//'/'//'EHgPS'// &
                            & yr_a2(nyear-2000)//chmon(month-1)
        else
          call die('geteh5om','ERROR IN geteh5om',1)
        end if
      else if ( ssttyp=='EH5RF' ) then
        finm = 'RF/'//yr_rf(nyear-1941)//'/'//'EHgRF'//yr_rf(nyear-1941)&
             & //chmon(12)
        if ( ehso4 ) psnm = 'RF/'//yr_rf(nyear-1941)//'/'//'EHgPS'//    &
                          & yr_rf(nyear-1941)//chmon(12)
      else if ( ssttyp=='EH5A2' ) then
        if ( nyear==2001 ) then
          finm = 'RF/'//yr_rf(nyear-1941)//'/'//'EHgRF'//               &
               & yr_rf(nyear-1941)//chmon(12)
          if ( ehso4 ) psnm = 'RF/'//yr_rf(nyear-1941)//'/'//'EHgPS'//  &
                            & yr_rf(nyear-1941)//chmon(12)
        else
          finm = 'A2/'//yr_a2(nyear-2001)//'/'//'EHgA2'//               &
               & yr_a2(nyear-2001)//chmon(12)
          if ( ehso4 ) psnm = 'A2/'//yr_a2(nyear-2001)//'/'//'EHgPS'//  &
                            & yr_a2(nyear-2001)//chmon(12)
        end if
      else if ( ssttyp=='EH5B1' ) then
        if ( nyear==2001 ) then
          finm = 'RF/'//yr_rf(nyear-1941)//'/'//'EHgRF'//               &
               & yr_rf(nyear-1941)//chmon(12)
          if ( ehso4 ) psnm = 'RF/'//yr_rf(nyear-1941)//'/'//'EHgPS'//  &
                            & yr_rf(nyear-1941)//chmon(12)
        else
          finm = 'B1/'//yr_a2(nyear-2001)//'/'//'EHgB1'//               &
               & yr_a2(nyear-2001)//chmon(12)
          if ( ehso4 ) psnm = 'B1/'//yr_a2(nyear-2001)//'/'//'EHgPS'//  &
                            & yr_a2(nyear-2001)//chmon(12)
        end if
      else if ( ssttyp=='EHA1B' ) then
        if ( nyear==2001 ) then
          finm = 'RF/'//yr_rf(nyear-1941)//'/'//'EHgRF'//               &
               & yr_rf(nyear-1941)//chmon(12)
          if ( ehso4 ) psnm = 'RF/'//yr_rf(nyear-1941)//'/'//'EHgPS'//  &
                            & yr_rf(nyear-1941)//chmon(12)
        else
          finm = 'A1B/'//yr_a2(nyear-2001)//'/'//'EgA1B'//              &
               & yr_a2(nyear-2001)//chmon(12)
          if ( ehso4 ) psnm = 'A1B/'//yr_a2(nyear-2001)//'/'//'EHgPS'// &
                            & yr_a2(nyear-2001)//chmon(12)
        end if
      else
        call die('geteh5om','ERROR IN geteh5om',1)
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
      if ( ehso4 ) then
        do k = 1 , mlev
          do j = 1 , jlat
            do i = 1 , ilon
              sulfate2(i,j,k) = -9999.
            end do
          end do
        end do
        do j = 1 , jlat
          do i = 1 , ilon
            pso4_2(i,j) = -9999.
          end do
        end do
      end if
 
      inquire (file=trim(inpglob)//'/EH5OM/'//finm,exist=there)
      if ( .not.there ) then
        call die('geteh5om',trim(inpglob)//'/EH5OM/'//finm// &
                 ' is not available',1)
      end if
      open (63,file=trim(inpglob)//'/EH5OM/'//finm,form='unformatted',  &
          & recl=(numx*numy*2+16)/4*ibyte,access='direct')
      if ( ehso4 ) then
        inquire (file=trim(inpglob)//'/EH5OM/'//psnm,exist=there)
        if ( .not.there ) then
          call die('geteh5om',trim(inpglob)//'/EH5OM/'//psnm// &
                   ' is not available',1)
        end if
        open (62,file=trim(inpglob)//'/EH5OM/'//psnm,form='unformatted',&
            & recl=(numx*numy*2+16)/4*ibyte,access='direct')
      end if
      if ( nday/=1 .or. nhour/=0 ) then
        nrec = ((nday-1)*4+nhour/6-1)*(klev*5)
        if ( ehso4 ) krec = (nday-1)*4 + nhour/6
      else if ( month==1 .or. month==2 .or. month==4 .or. month==6 .or. &
              & month==8 .or. month==9 .or. month==11 ) then
        nrec = (mlev*4-1)*(klev*5)
        if ( ehso4 ) krec = mlev*4
      else if ( month==5 .or. month==7 .or. month==10 .or. month==12 )  &
              & then
        nrec = (30*4-1)*(klev*5)
        if ( ehso4 ) krec = 30*4
      else if ( mod(nyear,4)==0 .and. nyear/=2100 ) then
        nrec = (29*4-1)*(klev*5)
        if ( ehso4 ) krec = 29*4
      else
        nrec = (28*4-1)*(klev*5)
        if ( ehso4 ) krec = 28*4
      end if
      if ( ehso4 ) then
        read (62,rec=krec) offset , xscale ,                            &
                         & ((itmp(i,j),i=1,numx),j=1,numy)
        do j = nint((lat0+.9375)/1.875) , nint((lat1+.9375)/1.875)
          do i = nint(lon0/1.875) , nint(lon1/1.875)
            ii = i + 1
            if ( ii<=0 ) ii = ii + ilon
            if ( ii>ilon ) ii = ii - ilon
            i2 = i - nint(lon0/1.875) + 1
            j2 = j - nint((lat0+.9375)/1.875) + 1
            if ( numx==ilon .and. numy==jlat ) then
              pso4_2(ii,49-j) = itmp(i2,j2)*xscale + offset +           &
                              & pso4_0(ii,49-j)
            else
              pso4_2(ii,j+jlat/2) = itmp(i2,j2)*xscale + offset +       &
                                & pso4_0(ii,j+jlat/2)
            end if
          end do
        end do
        close (62)
      end if
 
      do k = klev , 1 , -1
        nrec = nrec + 1
        read (63,rec=nrec) offset , xscale ,                            &
                         & ((itmp(i,j),i=1,numx),j=1,numy)
        do j = nint((lat0+.9375)/1.875) , nint((lat1+.9375)/1.875)
          do i = nint(lon0/1.875) , nint(lon1/1.875)
            ii = i + 1
            if ( ii<=0 ) ii = ii + ilon
            if ( ii>ilon ) ii = ii - ilon
            i2 = i - nint(lon0/1.875) + 1
            j2 = j - nint((lat0+.9375)/1.875) + 1
            if ( numx==ilon .and. numy==jlat ) then
              hvar(ii,49-j,k) = itmp(i2,j2)*xscale + offset
            else
              hvar(ii,j+jlat/2,k) = itmp(i2,j2)*xscale + offset
            end if
          end do
        end do
      end do
      do k = klev , 1 , -1
        nrec = nrec + 1
        read (63,rec=nrec) offset , xscale ,                            &
                         & ((itmp(i,j),i=1,numx),j=1,numy)
        do j = nint((lat0+.9375)/1.875) , nint((lat1+.9375)/1.875)
          do i = nint(lon0/1.875) , nint(lon1/1.875)
            ii = i + 1
            if ( ii<=0 ) ii = ii + ilon
            if ( ii>ilon ) ii = ii - ilon
            i2 = i - nint(lon0/1.875) + 1
            j2 = j - nint((lat0+.9375)/1.875) + 1
            if ( numx==ilon .and. numy==jlat ) then
              rhvar(ii,49-j,k) = dmin1(dmax1(itmp(i2,j2)*xscale+offset, &
                               & 0.D0),1.D0)
            else
              rhvar(ii,j+jlat/2,k) = dmin1(dmax1(itmp(i2,j2)*xscale+    &
                              & offset,0.D0),1.D0)
            end if
          end do
        end do
      end do
      do k = klev , 1 , -1
        nrec = nrec + 1
        read (63,rec=nrec) offset , xscale ,                            &
                         & ((itmp(i,j),i=1,numx),j=1,numy)
        do j = nint((lat0+.9375)/1.875) , nint((lat1+.9375)/1.875)
          do i = nint(lon0/1.875) , nint(lon1/1.875)
            ii = i + 1
            if ( ii<=0 ) ii = ii + ilon
            if ( ii>ilon ) ii = ii - ilon
            i2 = i - nint(lon0/1.875) + 1
            j2 = j - nint((lat0+.9375)/1.875) + 1
            if ( numx==ilon .and. numy==jlat ) then
              tvar(ii,49-j,k) = itmp(i2,j2)*xscale + offset
            else
              tvar(ii,j+jlat/2,k) = itmp(i2,j2)*xscale + offset
            end if
          end do
        end do
      end do
      do k = klev , 1 , -1
        nrec = nrec + 1
        read (63,rec=nrec) offset , xscale ,                            &
                         & ((itmp(i,j),i=1,numx),j=1,numy)
        do j = nint((lat0+.9375)/1.875) , nint((lat1+.9375)/1.875)
          do i = nint(lon0/1.875) , nint(lon1/1.875)
            ii = i + 1
            if ( ii<=0 ) ii = ii + ilon
            if ( ii>ilon ) ii = ii - ilon
            i2 = i - nint(lon0/1.875) + 1
            j2 = j - nint((lat0+.9375)/1.875) + 1
            if ( numx==ilon .and. numy==jlat ) then
              uvar(ii,49-j,k) = itmp(i2,j2)*xscale + offset
            else
              uvar(ii,j+jlat/2,k) = itmp(i2,j2)*xscale + offset
            end if
          end do
        end do
      end do
      do k = klev , 1 , -1
        nrec = nrec + 1
        read (63,rec=nrec) offset , xscale ,                            &
                         & ((itmp(i,j),i=1,numx),j=1,numy)
        do j = nint((lat0+.9375)/1.875) , nint((lat1+.9375)/1.875)
          do i = nint(lon0/1.875) , nint(lon1/1.875)
            ii = i + 1
            if ( ii<=0 ) ii = ii + ilon
            if ( ii>ilon ) ii = ii - ilon
            i2 = i - nint(lon0/1.875) + 1
            j2 = j - nint((lat0+.9375)/1.875) + 1
            if ( numx==ilon .and. numy==jlat ) then
              vvar(ii,49-j,k) = itmp(i2,j2)*xscale + offset
            else
              vvar(ii,j+jlat/2,k) = itmp(i2,j2)*xscale + offset
            end if
          end do
        end do
      end do
      close (63)
      write (stdout,*) 'READ IN fields at DATE:' , idate
!
!     HORIZONTAL INTERPOLATION OF BOTH THE SCALAR AND VECTOR FIELDS
!
      call bilinx2(b3,b2,xlon,xlat,glon,glat,ilon,jlat,jx,iy,klev*3)
      call bilinx2(d3,d2,dlon,dlat,glon,glat,ilon,jlat,jx,iy,klev*2)
!
!     ROTATE U-V FIELDS AFTER HORIZONTAL INTERPOLATION
!
      call uvrot4(u3,v3,dlon,dlat,clon,clat,grdfac,jx,iy,klev,plon,plat,&
                & iproj)
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
      if(i_band.eq.1) then
         call p1p2_band(b3pd,ps4,jx,iy)
      else
         call p1p2(b3pd,ps4,jx,iy)
      endif
 
!
!     F0  DETERMINE SURFACE TEMPS ON RCM TOPOGRAPHY.
!     INTERPOLATION FROM PRESSURE LEVELS AS IN INTV2
      call intv3(ts4,t3,ps4,sigmar,ptop,jx,iy,klev)
      call readsst(ts4,idate)
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
      if ( ehso4 ) then
        if ( ssttyp=='EH5RF' ) then
          fnso4 = trim(inpglob)//'/EH5OM/SO4/RF/T63L31_skg_'//       &
                   & yr_rf(nyear-1940)//'.nc'
          inquire (file=fnso4,exist=there)
          if ( .not.there ) then
            call die('geteh5om', fnso4//' is not available',1)
          end if
          istatus = nf90_open(fnso4,nf90_nowrite,ncid)
          if (istatus /= nf90_noerr) then
            call die('geteh5om', fnso4//':open',1, &
                     nf90_strerror(istatus),istatus)
          end if
          istatus = nf90_get_var(ncid,10,sulfate,istart,icount)
          if (istatus /= nf90_noerr) then
            call die('geteh5om', fnso4//':getvar',1, &
                     nf90_strerror(istatus),istatus)
          end if
          istatus = nf90_close(ncid)
          if (istatus /= nf90_noerr) then
            call die('geteh5om', fnso4//':close',1, &
                     nf90_strerror(istatus),istatus)
          end if
          if ( nyear==1950 .and. month==1 .and. nday<16 ) then
            do k = 1 , mlev
              do j = 1 , jlat
                do i = 1 , ilon
                  sulfate2(i,j,k) = sulfate(i,j,k,1)
                end do
              end do
            end do
          else if ( nyear==2000 .and. month==12 .and. nday>=16 ) then
            do k = 1 , mlev
              do j = 1 , jlat
                do i = 1 , ilon
                  sulfate2(i,j,k) = sulfate(i,j,k,12)
                end do
              end do
            end do
          else if ( (month==1 .or. month==3. .or. month==5 .or.         &
                  & month==7 .or. month==8 .or. month==10) .and.        &
                  & nday>=16 ) then
            if ( month==1 ) then
              do k = 1 , mlev
                do j = 1 , jlat
                  do i = 1 , ilon
                    sulfate2(i,j,k) = sulfate(i,j,k,month)              &
                                    & *(1.-float(nday-16)/30.)          &
                                    & + sulfate(i,j,k,month+1)          &
                                    & *(float(nday-16)/30.)
                  end do
                end do
              end do
            else if ( month==3 .or. month==5 .or. month==8 .or.         &
                    & month==10 ) then
              do k = 1 , mlev
                do j = 1 , jlat
                  do i = 1 , ilon
                    sulfate2(i,j,k) = sulfate(i,j,k,month)              &
                                    & *(1.-float(nday-16)/31.)          &
                                    & + sulfate(i,j,k,month+1)          &
                                    & *(float(nday-16)/31.)
                  end do
                end do
              end do
            else if ( month==7 ) then
              do k = 1 , mlev
                do j = 1 , jlat
                  do i = 1 , ilon
                    sulfate2(i,j,k) = sulfate(i,j,k,month)              &
                                    & *(1.-float(nday-16)/32.)          &
                                    & + sulfate(i,j,k,month+1)          &
                                    & *(float(nday-16)/32.)
                  end do
                end do
              end do
            else
            end if
          else if ( month==2 .and. nday>=15 ) then
            do k = 1 , mlev
              do j = 1 , jlat
                do i = 1 , ilon
                  sulfate2(i,j,k) = sulfate(i,j,k,month)                &
                                  & *(1.-float(nday-16)/30.)            &
                                  & + sulfate(i,j,k,month+1)            &
                                  & *(float(nday-16)/30.)
                end do
              end do
            end do
          else if ( (month==4. .or. month==6 .or. month==9 .or.         &
                  & month==11) .and. nday>=16 ) then
            do k = 1 , mlev
              do j = 1 , jlat
                do i = 1 , ilon
                  sulfate2(i,j,k) = sulfate(i,j,k,month)                &
                                  & *(1.-float(nday-16)/31.)            &
                                  & + sulfate(i,j,k,month+1)            &
                                  & *(float(nday-16)/31.)
                end do
              end do
            end do
          else if ( month==12. .and. nday>=16 ) then
            do k = 1 , mlev
              do j = 1 , jlat
                do i = 1 , ilon
                  sulfate1(i,j,k,1) = sulfate(i,j,k,12)
                end do
              end do
            end do
            fnso4 = trim(inpglob)//'/EH5OM/SO4/RF/T63L31_skg_'//     &
                     & yr_rf(nyear-1939)//'.nc'
            inquire (file=fnso4,exist=there)
            if ( .not.there ) then
              call die('geteh5om',fnso4//' is not available',1)
            end if
            istatus = nf90_open(fnso4,nf90_nowrite,ncid)
            if (istatus /= nf90_noerr) then
              call die('geteh5om', fnso4//':open',1, &
                       nf90_strerror(istatus),istatus)
            end if
            istatus = nf90_get_var(ncid,10,sulfate,istart,icount)
            if (istatus /= nf90_noerr) then
              call die('geteh5om', fnso4//':getvar',1, &
                       nf90_strerror(istatus),istatus)
            end if
            istatus = nf90_close(ncid)
            if (istatus /= nf90_noerr) then
              call die('geteh5om', fnso4//':close',1, &
                       nf90_strerror(istatus),istatus)
            end if
            do k = 1 , mlev
              do j = 1 , jlat
                do i = 1 , ilon
                  sulfate1(i,j,k,2) = sulfate(i,j,k,1)
                end do
              end do
            end do
            do k = 1 , mlev
              do j = 1 , jlat
                do i = 1 , ilon
                  sulfate2(i,j,k) = sulfate1(i,j,k,1)                   &
                                  & *(1.-float(nday-16)/32.)            &
                                  & + sulfate1(i,j,k,2)                 &
                                  & *(float(nday-16)/32.)
                end do
              end do
            end do
          else if ( (month==3 .or. month==5 .or. month==7 .or.          &
                  & month==8 .or. month==10 .or. month==12) .and.       &
                  & nday<16 ) then
            if ( month==3 ) then
              do k = 1 , mlev
                do j = 1 , jlat
                  do i = 1 , ilon
                    sulfate2(i,j,k) = sulfate(i,j,k,month-1)            &
                                    & *(float(16-nday)/30.)             &
                                    & + sulfate(i,j,k,month)            &
                                    & *(1.-float(16-nday)/30.)
                  end do
                end do
              end do
            else if ( month==5 .or. month==7 .or. month==10 .or.        &
                    & month==12 ) then
              do k = 1 , mlev
                do j = 1 , jlat
                  do i = 1 , ilon
                    sulfate2(i,j,k) = sulfate(i,j,k,month-1)            &
                                    & *(float(16-nday)/31.)             &
                                    & + sulfate(i,j,k,month)            &
                                    & *(1.-float(16-nday)/31.)
                  end do
                end do
              end do
            else if ( month==8 ) then
              do k = 1 , mlev
                do j = 1 , jlat
                  do i = 1 , ilon
                    sulfate2(i,j,k) = sulfate(i,j,k,month-1)            &
                                    & *(float(16-nday)/32.)             &
                                    & + sulfate(i,j,k,month)            &
                                    & *(1.-float(16-nday)/32.)
                  end do
                end do
              end do
            else
            end if
          else if ( month==2 .and. nday<15 ) then
            do k = 1 , mlev
              do j = 1 , jlat
                do i = 1 , ilon
                  sulfate2(i,j,k) = sulfate(i,j,k,month-1)              &
                                  & *(float(15-nday)/30.)               &
                                  & + sulfate(i,j,k,month)              &
                                  & *(1.-float(15-nday)/30.)
                end do
              end do
            end do
          else if ( (month==4. .or. month==6 .or. month==9 .or.         &
                  & month==11) .and. nday<16 ) then
            do k = 1 , mlev
              do j = 1 , jlat
                do i = 1 , ilon
                  sulfate2(i,j,k) = sulfate(i,j,k,month-1)              &
                                  & *(float(16-nday)/31.)               &
                                  & + sulfate(i,j,k,month)              &
                                  & *(1.-float(16-nday)/31.)
                end do
              end do
            end do
          else if ( month==1 .and. nday<16 ) then
            do k = 1 , mlev
              do j = 1 , jlat
                do i = 1 , ilon
                  sulfate1(i,j,k,2) = sulfate(i,j,k,1)
                end do
              end do
            end do
            fnso4 = trim(inpglob)//'/EH5OM/SO4/RF/T63L31_skg_'//     &
                     & yr_rf(nyear-1941)//'.nc'
            inquire (file=fnso4,exist=there)
            if ( .not.there ) then
              call die('geteh5om',fnso4//' is not available',1)
            end if
            istatus = nf90_open(fnso4,nf90_nowrite,ncid)
            if (istatus /= nf90_noerr) then
              call die('geteh5om', fnso4//':open',1, &
                       nf90_strerror(istatus),istatus)
            end if
            istatus = nf90_get_var(ncid,10,sulfate,istart,icount)
            if (istatus /= nf90_noerr) then
              call die('geteh5om', fnso4//':getvar',1, &
                       nf90_strerror(istatus),istatus)
            end if
            istatus = nf90_close(ncid)
            if (istatus /= nf90_noerr) then
              call die('geteh5om', fnso4//':close',1, &
                       nf90_strerror(istatus),istatus)
            end if
            do k = 1 , mlev
              do j = 1 , jlat
                do i = 1 , ilon
                  sulfate1(i,j,k,1) = sulfate(i,j,k,12)
                end do
              end do
            end do
            do k = 1 , mlev
              do j = 1 , jlat
                do i = 1 , ilon
                  sulfate2(i,j,k) = sulfate1(i,j,k,1)                   &
                                  & *(float(16-nday)/32.)               &
                                  & + sulfate1(i,j,k,2)                 &
                                  & *(1.-float(16-nday)/32.)
                end do
              end do
            end do
          else
          end if
        else
          if ( ssttyp=='EH5A2' ) then
            fnso4 = trim(inpglob)//'/EH5OM/SO4/A2/T63L31_skg_A2_'//  &
                     & yr_a2(nyear-2000)//'.nc'
          else if ( ssttyp=='EHA1B' ) then
            fnso4 = trim(inpglob)//                                 &
                 & '/EH5OM/SO4/A1B/T63L31_skg_A1B_'//                   &
                 & yr_a2(nyear-2000)//'.nc'
          else if ( ssttyp=='EH5B1' ) then
            fnso4 = trim(inpglob)//'/EH5OM/SO4/B1/T63L31_skg_B1_'//  &
                     & yr_a2(nyear-2000)//'.nc'
          end if
          inquire (file=fnso4,exist=there)
          if ( .not.there ) then
            call die('geteh5om',fnso4//' is not available',1)
          end if
          istatus = nf90_open(fnso4,nf90_nowrite,ncid)
          if (istatus /= nf90_noerr) then
            call die('geteh5om', fnso4//':open',1, &
                     nf90_strerror(istatus),istatus)
          end if
          istatus = nf90_get_var(ncid,10,sulfate,istart,icount)
          if (istatus /= nf90_noerr) then
            call die('geteh5om', fnso4//':getvar',1, &
                     nf90_strerror(istatus),istatus)
          end if
          istatus = nf90_close(ncid)
          if (istatus /= nf90_noerr) then
            call die('geteh5om', fnso4//':close',1, &
                     nf90_strerror(istatus),istatus)
          end if
          if ( nyear==2001 .and. month==1 .and. nday<16 ) then
            do k = 1 , mlev
              do j = 1 , jlat
                do i = 1 , ilon
                  sulfate2(i,j,k) = sulfate(i,j,k,1)
                end do
              end do
            end do
          else if ( nyear==2100 .and. month==12 .and. nday>=16 ) then
            do k = 1 , mlev
              do j = 1 , jlat
                do i = 1 , ilon
                  sulfate2(i,j,k) = sulfate(i,j,k,12)
                end do
              end do
            end do
          else if ( (month==1 .or. month==3. .or. month==5 .or.         &
                  & month==7 .or. month==8 .or. month==10) .and.        &
                  & nday>=16 ) then
            if ( month==1 ) then
              do k = 1 , mlev
                do j = 1 , jlat
                  do i = 1 , ilon
                    sulfate2(i,j,k) = sulfate(i,j,k,month)              &
                                    & *(1.-float(nday-16)/30.)          &
                                    & + sulfate(i,j,k,month+1)          &
                                    & *(float(nday-16)/30.)
                  end do
                end do
              end do
            else if ( month==3 .or. month==5 .or. month==8 .or.         &
                    & month==10 ) then
              do k = 1 , mlev
                do j = 1 , jlat
                  do i = 1 , ilon
                    sulfate2(i,j,k) = sulfate(i,j,k,month)              &
                                    & *(1.-float(nday-16)/31.)          &
                                    & + sulfate(i,j,k,month+1)          &
                                    & *(float(nday-16)/31.)
                  end do
                end do
              end do
            else if ( month==7 ) then
              do k = 1 , mlev
                do j = 1 , jlat
                  do i = 1 , ilon
                    sulfate2(i,j,k) = sulfate(i,j,k,month)              &
                                    & *(1.-float(nday-16)/32.)          &
                                    & + sulfate(i,j,k,month+1)          &
                                    & *(float(nday-16)/32.)
                  end do
                end do
              end do
            else
            end if
          else if ( month==2 .and. nday>=15 ) then
            do k = 1 , mlev
              do j = 1 , jlat
                do i = 1 , ilon
                  sulfate2(i,j,k) = sulfate(i,j,k,month)                &
                                  & *(1.-float(nday-16)/30.)            &
                                  & + sulfate(i,j,k,month+1)            &
                                  & *(float(nday-16)/30.)
                end do
              end do
            end do
          else if ( (month==4. .or. month==6 .or. month==9 .or.         &
                  & month==11) .and. nday>=16 ) then
            do k = 1 , mlev
              do j = 1 , jlat
                do i = 1 , ilon
                  sulfate2(i,j,k) = sulfate(i,j,k,month)                &
                                  & *(1.-float(nday-16)/31.)            &
                                  & + sulfate(i,j,k,month+1)            &
                                  & *(float(nday-16)/31.)
                end do
              end do
            end do
          else if ( month==12. .and. nday>=16 ) then
            do k = 1 , mlev
              do j = 1 , jlat
                do i = 1 , ilon
                  sulfate1(i,j,k,1) = sulfate(i,j,k,12)
                end do
              end do
            end do
            if ( ssttyp=='EH5A2' ) then
              fnso4 = trim(inpglob)//'/EH5OM/SO4/A2/T63L31_skg_A2_'//&
                       & yr_a2(nyear-1999)//'.nc'
            else if ( ssttyp=='EHA1B' ) then
              fnso4 = trim(inpglob)//                               &
                        & '/EH5OM/SO4/A1B/T63L31_skg_A1B_'//            &
                        & yr_a2(nyear-1999)//'.nc'
            else if ( ssttyp=='EH5B1' ) then
              fnso4 = trim(inpglob)//'/EH5OM/SO4/B1/T63L31_skg_B1_'//&
                       & yr_a2(nyear-1999)//'.nc'
            end if
            inquire (file=fnso4,exist=there)
            if ( .not.there ) then
              call die('geteh5om',fnso4//' is not available',1)
            end if
            istatus = nf90_open(fnso4,nf90_nowrite,ncid)
            if (istatus /= nf90_noerr) then
              call die('geteh5om', fnso4//':open',1, &
                       nf90_strerror(istatus),istatus)
            end if
            istatus = nf90_get_var(ncid,10,sulfate,istart,icount)
            if (istatus /= nf90_noerr) then
              call die('geteh5om', fnso4//':getvar',1, &
                       nf90_strerror(istatus),istatus)
            end if
            istatus = nf90_close(ncid)
            if (istatus /= nf90_noerr) then
              call die('geteh5om', fnso4//':close',1, &
                       nf90_strerror(istatus),istatus)
            end if
            do k = 1 , mlev
              do j = 1 , jlat
                do i = 1 , ilon
                  sulfate1(i,j,k,2) = sulfate(i,j,k,1)
                end do
              end do
            end do
            do k = 1 , mlev
              do j = 1 , jlat
                do i = 1 , ilon
                  sulfate2(i,j,k) = sulfate1(i,j,k,1)                   &
                                  & *(1.-float(nday-16)/32.)            &
                                  & + sulfate1(i,j,k,2)                 &
                                  & *(float(nday-16)/32.)
                end do
              end do
            end do
          else if ( (month==3 .or. month==5 .or. month==7 .or.          &
                  & month==8 .or. month==10 .or. month==12) .and.       &
                  & nday<16 ) then
            if ( month==3 ) then
              do k = 1 , mlev
                do j = 1 , jlat
                  do i = 1 , ilon
                    sulfate2(i,j,k) = sulfate(i,j,k,month-1)            &
                                    & *(float(16-nday)/30.)             &
                                    & + sulfate(i,j,k,month)            &
                                    & *(1.-float(16-nday)/30.)
                  end do
                end do
              end do
            else if ( month==5 .or. month==7 .or. month==10 .or.        &
                    & month==12 ) then
              do k = 1 , mlev
                do j = 1 , jlat
                  do i = 1 , ilon
                    sulfate2(i,j,k) = sulfate(i,j,k,month-1)            &
                                    & *(float(16-nday)/31.)             &
                                    & + sulfate(i,j,k,month)            &
                                    & *(1.-float(16-nday)/31.)
                  end do
                end do
              end do
            else if ( month==8 ) then
              do k = 1 , mlev
                do j = 1 , jlat
                  do i = 1 , ilon
                    sulfate2(i,j,k) = sulfate(i,j,k,month-1)            &
                                    & *(float(16-nday)/32.)             &
                                    & + sulfate(i,j,k,month)            &
                                    & *(1.-float(16-nday)/32.)
                  end do
                end do
              end do
            else
            end if
          else if ( month==2 .and. nday<15 ) then
            do k = 1 , mlev
              do j = 1 , jlat
                do i = 1 , ilon
                  sulfate2(i,j,k) = sulfate(i,j,k,month-1)              &
                                  & *(float(15-nday)/30.)               &
                                  & + sulfate(i,j,k,month)              &
                                  & *(1.-float(15-nday)/30.)
                end do
              end do
            end do
          else if ( (month==4. .or. month==6 .or. month==9 .or.         &
                  & month==11) .and. nday<16 ) then
            do k = 1 , mlev
              do j = 1 , jlat
                do i = 1 , ilon
                  sulfate2(i,j,k) = sulfate(i,j,k,month-1)              &
                                  & *(float(16-nday)/31.)               &
                                  & + sulfate(i,j,k,month)              &
                                  & *(1.-float(16-nday)/31.)
                end do
              end do
            end do
          else if ( month==1 .and. nday<16 ) then
            do k = 1 , mlev
              do j = 1 , jlat
                do i = 1 , ilon
                  sulfate1(i,j,k,2) = sulfate(i,j,k,1)
                end do
              end do
            end do
            if ( ssttyp=='EH5A2' ) then
              fnso4 = trim(inpglob)//'/EH5OM/SO4/A2/T63L31_skg_A2_'//&
                       & yr_a2(nyear-2001)//'.nc'
            else if ( ssttyp=='EHA1B' ) then
              fnso4 = trim(inpglob)//                               &
                        & '/EH5OM/SO4/A1B/T63L31_skg_A1B_'//            &
                        & yr_a2(nyear-2001)//'.nc'
            else if ( ssttyp=='EH5B1' ) then
              fnso4 = trim(inpglob)//'/EH5OM/SO4/B1/T63L31_skg_B1_'//&
                       & yr_a2(nyear-2001)//'.nc'
            end if
            inquire (file=fnso4,exist=there)
            if ( .not.there ) then
              call die('geteh5om',fnso4//' is not available',1)
            end if
            istatus = nf90_open(fnso4,nf90_nowrite,ncid)
            if (istatus /= nf90_noerr) then
              call die('geteh5om', fnso4//':open',1, &
                       nf90_strerror(istatus),istatus)
            end if
            istatus = nf90_get_var(ncid,10,sulfate,istart,icount)
            if (istatus /= nf90_noerr) then
              call die('geteh5om', fnso4//':getvar',1, &
                       nf90_strerror(istatus),istatus)
            end if
            istatus = nf90_close(ncid)
            if (istatus /= nf90_noerr) then
              call die('geteh5om', fnso4//':close',1, &
                       nf90_strerror(istatus),istatus)
            end if
            do k = 1 , mlev
              do j = 1 , jlat
                do i = 1 , ilon
                  sulfate1(i,j,k,1) = sulfate(i,j,k,12)
                end do
              end do
            end do
            do k = 1 , mlev
              do j = 1 , jlat
                do i = 1 , ilon
                  sulfate2(i,j,k) = sulfate1(i,j,k,1)                   &
                                  & *(float(16-nday)/32.)               &
                                  & + sulfate1(i,j,k,2)                 &
                                  & *(1.-float(16-nday)/32.)
                end do
              end do
            end do
          else
          end if
        end if
        call bilinx2(sulfate3,sulfate2,xlon,xlat,glon,glat,             &
                  &  ilon,jlat,jx,iy,mlev)
        call bilinx2(pso4_3,pso4_2,xlon,xlat,glon,glat,                 &
                  &  ilon,jlat,jx,iy,1)
        do i = 1 , iy
          do j = 1 , jx
            do l = 1 , kz
              prcm = ((ps4(j,i)-ptop)*sigma2(l)+ptop)*10.
              k0 = -1
              do k = mlev , 1 , -1
                pmpi = (pso4_3(j,i)*hybm(k)+hyam(k))*0.01
                k0 = k
                if ( prcm>pmpi ) exit
              end do
              if ( k0==mlev ) then
                pmpj = (pso4_3(j,i)*hybm(mlev-1)+hyam(mlev-1))*0.01
                pmpi = (pso4_3(j,i)*hybm(mlev)+hyam(mlev))*0.01
                sulfate4(j,i,l) = sulfate3(j,i,mlev)                    &
                                & + (sulfate3(j,i,mlev)                 &
                                & -sulfate3(j,i,mlev-1))*(prcm-pmpi)    &
                                & /(pmpi-pmpj)
              else if ( k0>=1 ) then
                pmpj = (pso4_3(j,i)*hybm(k0)+hyam(k0))*0.01
                pmpi = (pso4_3(j,i)*hybm(k0+1)+hyam(k0+1))*0.01
                sulfate4(j,i,l) = (sulfate3(j,i,k0+1)*(prcm-pmpj)+      &
                                & sulfate3(j,i,k0)*(pmpi-prcm))         &
                                & /(pmpi-pmpj)
              else
              end if
            end do
          end do
        end do
!
      end if

      call zeit_co('get_eh5om')
      end subroutine geteh5om

      subroutine headermpi(ehso4)
      implicit none
!
      logical :: ehso4
      intent (in) ehso4
!
      integer :: i , k , kr , ierr
      logical :: there
!
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
!
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
!
!     INITIAL GLOBAL GRID-POINT LONGITUDE & LATITUDE
!
      do i = 1 , ilon
        glon(i) = float(i-1)*1.875
      end do
!     CHANGE ORDER OF VERTICAL INDEXES FOR PRESSURE LEVELS
!
      do k = 1 , klev
        kr = klev - k + 1
        sigma1(k) = sigmar(kr)
      end do
 
      if ( ehso4 ) then
        hyai(1) = 0.
        hyai(2) = 2000.
        hyai(3) = 4000.
        hyai(4) = 6000.
        hyai(5) = 8000.
        hyai(6) = 9976.13671875
        hyai(7) = 11820.5400390625
        hyai(8) = 13431.3896484375
        hyai(9) = 14736.3603515625
        hyai(10) = 15689.2099609375
        hyai(11) = 16266.6103515625
        hyai(12) = 16465.
        hyai(13) = 16297.6201171875
        hyai(14) = 15791.599609375
        hyai(15) = 14985.26953125
        hyai(16) = 13925.51953125
        hyai(17) = 12665.2900390625
        hyai(18) = 11261.23046875
        hyai(19) = 9771.40625
        hyai(20) = 8253.2109375
        hyai(21) = 6761.33984375
        hyai(22) = 5345.9140625
        hyai(23) = 4050.71801757812
        hyai(24) = 2911.56909179688
        hyai(25) = 1954.80505371094
        hyai(26) = 1195.89001464844
        hyai(27) = 638.14892578125
        hyai(28) = 271.626495361328
        hyai(29) = 72.0635833740234
        hyai(30) = 0.
        hyai(31) = 0.
        hyai(32) = 0.
 
        hybi(1) = 0.
        hybi(2) = 0.
        hybi(3) = 0.
        hybi(4) = 0.
        hybi(5) = 0.
        hybi(6) = 0.000390858185710385
        hybi(7) = 0.0029197009280324
        hybi(8) = 0.00919413194060326
        hybi(9) = 0.0203191600739956
        hybi(10) = 0.0369748584926128
        hybi(11) = 0.0594876408576965
        hybi(12) = 0.0878949835896492
        hybi(13) = 0.122003600001335
        hybi(14) = 0.161441504955292
        hybi(15) = 0.205703303217888
        hybi(16) = 0.254188597202301
        hybi(17) = 0.306235402822495
        hybi(18) = 0.361144989728928
        hybi(19) = 0.418202310800552
        hybi(20) = 0.476688086986542
        hybi(21) = 0.535886585712433
        hybi(22) = 0.595084190368652
        hybi(23) = 0.65356457233429
        hybi(24) = 0.710594415664673
        hybi(25) = 0.765405178070068
        hybi(26) = 0.817166984081268
        hybi(27) = 0.86495578289032
        hybi(28) = 0.907715916633606
        hybi(29) = 0.944213211536407
        hybi(30) = 0.972985208034515
        hybi(31) = 0.992281496524811
        hybi(32) = 1.0
 
        hyam(1) = 1000.
        hyam(2) = 3000.
        hyam(3) = 5000.
        hyam(4) = 7000.
        hyam(5) = 8988.068359375
        hyam(6) = 10898.3383789062
        hyam(7) = 12625.96484375
        hyam(8) = 14083.875
        hyam(9) = 15212.78515625
        hyam(10) = 15977.91015625
        hyam(11) = 16365.8051757812
        hyam(12) = 16381.3100585938
        hyam(13) = 16044.6098632812
        hyam(14) = 15388.4345703125
        hyam(15) = 14455.39453125
        hyam(16) = 13295.4047851562
        hyam(17) = 11963.2602539062
        hyam(18) = 10516.318359375
        hyam(19) = 9012.30859375
        hyam(20) = 7507.275390625
        hyam(21) = 6053.626953125
        hyam(22) = 4698.31604003906
        hyam(23) = 3481.1435546875
        hyam(24) = 2433.18707275391
        hyam(25) = 1575.34753417969
        hyam(26) = 917.019470214844
        hyam(27) = 454.887710571289
        hyam(28) = 171.845039367676
        hyam(29) = 36.0317916870117
        hyam(30) = 0.
        hyam(31) = 0.
 
        hybm(1) = 0.
        hybm(2) = 0.
        hybm(3) = 0.
        hybm(4) = 0.
        hybm(5) = 0.000195429092855193
        hybm(6) = 0.00165527955687139
        hybm(7) = 0.00605691643431783
        hybm(8) = 0.0147566460072994
        hybm(9) = 0.0286470092833042
        hybm(10) = 0.0482312496751547
        hybm(11) = 0.0736913122236729
        hybm(12) = 0.104949291795492
        hybm(13) = 0.141722552478313
        hybm(14) = 0.18357240408659
        hybm(15) = 0.229945950210094
        hybm(16) = 0.280212000012398
        hybm(17) = 0.333690196275711
        hybm(18) = 0.38967365026474
        hybm(19) = 0.447445198893547
        hybm(20) = 0.506287336349487
        hybm(21) = 0.565485388040543
        hybm(22) = 0.624324381351471
        hybm(23) = 0.682079493999481
        hybm(24) = 0.737999796867371
        hybm(25) = 0.791286081075668
        hybm(26) = 0.841061383485794
        hybm(27) = 0.886335849761963
        hybm(28) = 0.925964564085007
        hybm(29) = 0.958599209785461
        hybm(30) = 0.982633352279663
        hybm(31) = 0.996140748262405
 
        istart(1) = 1
        istart(2) = 1
        istart(3) = 1
        istart(4) = 1
        istart(5) = 0
        istart(6) = 0
        istart(7) = 0
        istart(8) = 0
        istart(9) = 0
        istart(10) = 0
 
        icount(1) = ilon
        icount(2) = jlat
        icount(3) = mlev
        icount(4) = 12
        icount(5) = 0
        icount(6) = 0
        icount(7) = 0
        icount(8) = 0
        icount(9) = 0
        icount(10) = 0
 
        inquire (file=trim(inpglob)//'/EH5OM/EHgPS.dat',exist=there)
        if ( .not.there ) then
          call die('headermpi',trim(inpglob)//'/EH5OM/EHgPS.dat '// &
                   'is not available',1)
        end if
        open (30,file=trim(inpglob)//'/EH5OM/EHgPS.dat',                &
            & form='unformatted',recl=ilon*jlat*4,access='direct')
        read (30,rec=1) pso4_0
        close (30)
      end if
 
      allocate(b3(jx,iy,klev*3), stat=ierr)
      if (ierr /= 0) call die('headermpi','allocate b3',ierr)
      call mall_mci(b3,'mod_eh5om')
      allocate(d3(jx,iy,klev*2), stat=ierr)
      if (ierr /= 0) call die('headermpi','allocate d3',ierr)
      call mall_mci(d3,'mod_eh5om')
      if ( ehso4 ) then
        allocate(sulfate3(jx,iy,mlev), stat=ierr)
        if (ierr /= 0) call die('headermpi','allocate sulfate3',ierr)
        call mall_mci(sulfate3,'mod_eh5om')
        allocate(pso4_3(jx,iy), stat=ierr)
        if (ierr /= 0) call die('headermpi','allocate pso4_3',ierr)
        call mall_mci(pso4_3,'mod_eh5om')
      end if

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

      end subroutine headermpi

      subroutine footermpi(ehso4)
        implicit none
        logical , intent(in) :: ehso4
        call mall_mco(d3,'mod_eh5om')
        deallocate(d3)
        call mall_mco(b3,'mod_eh5om')
        deallocate(b3)
        if ( ehso4 ) then
          call mall_mco(sulfate3,'mod_eh5om')
          deallocate(sulfate3)
          call mall_mco(pso4_3,'mod_eh5om')
          deallocate(pso4_3)
        end if
      end subroutine footermpi

      end module mod_eh5om

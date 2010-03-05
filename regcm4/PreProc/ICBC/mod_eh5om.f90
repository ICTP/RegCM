      module mod_eh5om
      use mod_param
      implicit none

      integer , parameter :: klev = 17 , jlat = 96 , ilon = 192
      integer , parameter :: mlev = 31

      real , target , dimension(ilon,jlat,klev*3) :: b2
      real , target , dimension(ilon,jlat,klev*2) :: d2
      real , target , dimension(jx,iy,klev*3) :: b3
      real , target , dimension(jx,iy,klev*2) :: d3

      real , pointer :: u3(:,:,:) , v3(:,:,:)
      real , pointer :: h3(:,:,:) , q3(:,:,:) , t3(:,:,:)
      real , pointer :: uvar(:,:,:) , vvar(:,:,:)
      real , pointer :: hvar(:,:,:) , rhvar(:,:,:) , tvar(:,:,:)

      real , dimension(jlat) :: glat
      real , dimension(ilon) :: glon
      real , dimension(klev) :: sigma1 , sigmar

      real , dimension(mlev+1) :: hyai , hybi
      real , dimension(mlev) :: hyam , hybm
      real(4) , dimension(ilon,jlat) :: pso4_0
      real(4) , dimension(ilon,jlat,mlev,12) :: sulfate

      real , dimension(ilon,jlat) :: pso4_2
      real , dimension(jx,iy) :: pso4_3
      real , dimension(ilon,jlat,mlev,2) :: sulfate1
      real , dimension(ilon,jlat,mlev) :: sulfate2
      real , dimension(jx,iy,mlev) :: sulfate3
      real , dimension(jx,iy,kz) :: sulfate4

      contains

      subroutine geteh5om(idate)
      use netcdf
      use mod_grid
      use mod_write
      use mod_var4
      implicit none
!
! Dummy arguments
!
      integer :: idate
!
! Local variables
!
      real , dimension(jx,iy) :: b3pd , pa , sst1 , tlayer , za
      character(3) , dimension(12) :: chmon
      character(21) :: finm , psnm
      character(44) :: fnso4_a1b
      character(42) :: fnso4_a2 , fnso4_b1
      character(39) :: fnso4_rf
      integer :: i , i2 , ii , j , j2 , k , k0 , krec , l , month ,     &
               & nday , nhour , nrec , numx , numy , nyear
      integer(2) , dimension(192,96) :: itmp
      real(8) :: offset , xscale
      real :: pmpi , pmpj , prcm
      logical :: there
      character(4) , dimension(100) :: yr_a2
      character(4) , dimension(61) :: yr_rf
      integer(4) , dimension(10) :: icount , istart
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
!
!     D      BEGIN LOOP OVER NTIMES
!
      nyear = idate/1000000
                         !, SST2(JX,IY)
      month = idate/10000 - nyear*100
      nday = idate/100 - (idate/10000)*100
      nhour = mod(idate,100)
 
      if ( ssttyp=='EH5RF' ) then
        if ( idate<1941010106 ) then
          write (*,*) 'EH5RF dataset is just available from 1941010106'
          stop
        end if
        if ( ehso4 ) then
          if ( idate<1950010100 ) then
            write (*,*) 'the monthly EH5RF sulfate data are just' ,     &
                       &' available from 1950010100'
            stop
          end if
        end if
        if ( idate>2001010100 ) then
          write (*,*) 'EH5RF dataset is just available  to  2001010100'
          stop
        end if
      end if
      if ( ssttyp=='EH5A2' ) then
        if ( idate<2001010100 ) then
          write (*,*) 'EH5A2 dataset is just avaiable from 2001010100'
          stop
        end if
        if ( idate>=2101010100 ) then
          write (*,*) 'EH5A2 dataset is just avaiable  to  2100123118'
          stop
        end if
      end if
      if ( ssttyp=='EH5B1' ) then
        if ( idate<2001010100 ) then
          write (*,*) 'EH5B1 dataset is just avaiable from 2001010100'
          stop
        end if
        if ( idate>=2101010100 ) then
          write (*,*) 'EH5B1 dataset is just avaiable  to  2100123118'
          stop
        end if
      end if
      if ( ssttyp=='EHA1B' ) then
        if ( idate<2001010100 ) then
          write (*,*) 'EHA1B dataset is just avaiable from 2001010100'
          stop
        end if
        if ( idate>=2101010100 ) then
          write (*,*) 'EHA1B dataset is just avaiable  to  2100123118'
          stop
        end if
      end if
      numx = nint((lon1-lon0)/1.875) + 1
      numy = nint((lat1-lat0)/1.875) + 1
      if ( numx/=192 .or. numy/=96 ) then
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
            write (*,*) 'ERROR IN SSTTYP'
            stop
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
            write (*,*) 'ERROR IN SSTTYP'
            stop
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
          write (*,*) 'ERROR IN SSTTYP'
          stop
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
          write (*,*) 'ERROR IN SSTTYP'
          stop
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
          write (*,*) 'ERROR IN SSTTYP'
          stop
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
        write (*,*) 'ERROR IN SSTTYP'
        stop
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
 
      inquire (file='../DATA/EH5OM/'//finm,exist=there)
      if ( .not.there ) then
        write (*,*) '../DATA/EH5OM/'//finm , ' is not available'
        write (*,*) 'please copy EH5OM output under ../DATA/EH5OM/'
        stop
      end if
      open (63,file='../DATA/EH5OM/'//finm,form='unformatted',          &
          & recl=(numx*numy*2+16)/4*ibyte,access='direct')
      if ( ehso4 ) then
        inquire (file='../DATA/EH5OM/'//psnm,exist=there)
        if ( .not.there ) then
          write (*,*) '../DATA/EH5OM/'//psnm , ' is not available'
          write (*,*) 'please copy EH5OM output under ../DATA/EH5OM/'
          stop
        end if
        open (62,file='../DATA/EH5OM/'//psnm,form='unformatted',        &
            & recl=(numx*numy*2+16)/4*ibyte,access='direct')
      end if
      if ( nday/=1 .or. nhour/=0 ) then
        nrec = ((nday-1)*4+nhour/6-1)*(klev*5)
        if ( ehso4 ) krec = (nday-1)*4 + nhour/6
      else if ( month==1 .or. month==2 .or. month==4 .or. month==6 .or. &
              & month==8 .or. month==9 .or. month==11 ) then
        nrec = (31*4-1)*(klev*5)
        if ( ehso4 ) krec = 31*4
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
        read (62,rec=krec) offset , xscale ,                             &
                         & ((itmp(i,j),i=1,numx),j=1,numy)
        do j = nint((lat0+.9375)/1.875) , nint((lat1+.9375)/1.875)
          do i = nint(lon0/1.875) , nint(lon1/1.875)
            ii = i + 1
            if ( ii<=0 ) ii = ii + 192
            if ( ii>192 ) ii = ii - 192
            i2 = i - nint(lon0/1.875) + 1
            j2 = j - nint((lat0+.9375)/1.875) + 1
            if ( numx==192 .and. numy==96 ) then
              pso4_2(ii,49-j) = itmp(i2,j2)*xscale + offset +            &
                              & pso4_0(ii,49-j)
            else
              pso4_2(ii,j+96/2) = itmp(i2,j2)*xscale + offset +          &
                                & pso4_0(ii,j+96/2)
            end if
          end do
        end do
        close (62)
      end if
 
      do k = klev , 1 , -1
        nrec = nrec + 1
        read (63,rec=nrec) offset , xscale ,                             &
                         & ((itmp(i,j),i=1,numx),j=1,numy)
        do j = nint((lat0+.9375)/1.875) , nint((lat1+.9375)/1.875)
          do i = nint(lon0/1.875) , nint(lon1/1.875)
            ii = i + 1
            if ( ii<=0 ) ii = ii + 192
            if ( ii>192 ) ii = ii - 192
            i2 = i - nint(lon0/1.875) + 1
            j2 = j - nint((lat0+.9375)/1.875) + 1
            if ( numx==192 .and. numy==96 ) then
              hvar(ii,49-j,k) = itmp(i2,j2)*xscale + offset
            else
              hvar(ii,j+96/2,k) = itmp(i2,j2)*xscale + offset
            end if
          end do
        end do
      end do
      do k = klev , 1 , -1
        nrec = nrec + 1
        read (63,rec=nrec) offset , xscale ,                             &
                         & ((itmp(i,j),i=1,numx),j=1,numy)
        do j = nint((lat0+.9375)/1.875) , nint((lat1+.9375)/1.875)
          do i = nint(lon0/1.875) , nint(lon1/1.875)
            ii = i + 1
            if ( ii<=0 ) ii = ii + 192
            if ( ii>192 ) ii = ii - 192
            i2 = i - nint(lon0/1.875) + 1
            j2 = j - nint((lat0+.9375)/1.875) + 1
            if ( numx==192 .and. numy==96 ) then
              rhvar(ii,49-j,k) = dmin1(dmax1(itmp(i2,j2)*xscale+offset,  &
                               & 0.D0),1.D0)
            else
              rhvar(ii,j+96/2,k) = dmin1(dmax1(itmp(i2,j2)*xscale+offset,&
                                 & 0.D0),1.D0)
            end if
          end do
        end do
      end do
      do k = klev , 1 , -1
        nrec = nrec + 1
        read (63,rec=nrec) offset , xscale ,                             &
                         & ((itmp(i,j),i=1,numx),j=1,numy)
        do j = nint((lat0+.9375)/1.875) , nint((lat1+.9375)/1.875)
          do i = nint(lon0/1.875) , nint(lon1/1.875)
            ii = i + 1
            if ( ii<=0 ) ii = ii + 192
            if ( ii>192 ) ii = ii - 192
            i2 = i - nint(lon0/1.875) + 1
            j2 = j - nint((lat0+.9375)/1.875) + 1
            if ( numx==192 .and. numy==96 ) then
              tvar(ii,49-j,k) = itmp(i2,j2)*xscale + offset
            else
              tvar(ii,j+96/2,k) = itmp(i2,j2)*xscale + offset
            end if
          end do
        end do
      end do
      do k = klev , 1 , -1
        nrec = nrec + 1
        read (63,rec=nrec) offset , xscale ,                             &
                         & ((itmp(i,j),i=1,numx),j=1,numy)
        do j = nint((lat0+.9375)/1.875) , nint((lat1+.9375)/1.875)
          do i = nint(lon0/1.875) , nint(lon1/1.875)
            ii = i + 1
            if ( ii<=0 ) ii = ii + 192
            if ( ii>192 ) ii = ii - 192
            i2 = i - nint(lon0/1.875) + 1
            j2 = j - nint((lat0+.9375)/1.875) + 1
            if ( numx==192 .and. numy==96 ) then
              uvar(ii,49-j,k) = itmp(i2,j2)*xscale + offset
            else
              uvar(ii,j+96/2,k) = itmp(i2,j2)*xscale + offset
            end if
          end do
        end do
      end do
      do k = klev , 1 , -1
        nrec = nrec + 1
        read (63,rec=nrec) offset , xscale ,                             &
                         & ((itmp(i,j),i=1,numx),j=1,numy)
        do j = nint((lat0+.9375)/1.875) , nint((lat1+.9375)/1.875)
          do i = nint(lon0/1.875) , nint(lon1/1.875)
            ii = i + 1
            if ( ii<=0 ) ii = ii + 192
            if ( ii>192 ) ii = ii - 192
            i2 = i - nint(lon0/1.875) + 1
            j2 = j - nint((lat0+.9375)/1.875) + 1
            if ( numx==192 .and. numy==96 ) then
              vvar(ii,49-j,k) = itmp(i2,j2)*xscale + offset
            else
              vvar(ii,j+96/2,k) = itmp(i2,j2)*xscale + offset
            end if
          end do
        end do
      end do
      close (63)
      write (*,*) 'READ IN fields at DATE:' , idate
!
!     HORIZONTAL INTERPOLATION OF BOTH THE SCALAR AND VECTOR FIELDS
!
      call bilinx(b3,b2,xlon,xlat,glon,glat,ilon,jlat,jx,iy,klev*3)
      call bilinx(d3,d2,dlon,dlat,glon,glat,ilon,jlat,jx,iy,klev*2)
!
!     ROTATE U-V FIELDS AFTER HORIZONTAL INTERPOLATION
!
      call uvrot4(u3,v3,dlon,dlat,clon,clat,grdfac,jx,iy,klev,plon,plat,&
                & lgtype)
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
      call p1p2(b3pd,ps4,jx,iy)
 
!
!     F0  DETERMINE SURFACE TEMPS ON RCM TOPOGRAPHY.
!     INTERPOLATION FROM PRESSURE LEVELS AS IN INTV2
      call intv3(ts4,t3,ps4,sigmar,ptop,jx,iy,klev)
 
!     F1  CALCULATE SSTS FOR DATE FROM OBSERVED SSTS
!     PRINT *, 'INPUT DAY FOR SST DATA ACQUISITION:', IDATE
      call mksst3(ts4,sst1,topogm,xlandu,jx,iy,idate)
 
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
!     G   WRITE AN INITIAL FILE FOR THE RCM
      call writef(u4,v4,t4,q4,ps4,ts4,ptop,jx,iy,kz,idate)
      if ( ehso4 ) then
        if ( ssttyp=='EH5RF' ) then
          fnso4_rf = '../DATA/EH5OM/SO4/RF/T63L31_skg_'//               &
                   & yr_rf(nyear-1940)//'.nc'
          inquire (file=fnso4_rf,exist=there)
          if ( .not.there ) then
            write (*,*) fnso4_rf , ' is not available'
            stop
          end if
          istatus = nf90_open(fnso4_rf,nf90_nowrite,ncid)
          istatus = nf90_get_var(ncid,10,sulfate,istart,icount)
          istatus = nf90_close(ncid)
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
            fnso4_rf = '../DATA/EH5OM/SO4/RF/T63L31_skg_'//             &
                     & yr_rf(nyear-1939)//'.nc'
            inquire (file=fnso4_rf,exist=there)
            if ( .not.there ) then
              write (*,*) fnso4_rf , ' is not available'
              stop
            end if
            istatus = nf90_open(fnso4_rf,nf90_nowrite,ncid)
            istatus = nf90_get_var(ncid,10,sulfate,istart,icount)
            istatus = nf90_close(ncid)
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
            fnso4_rf = '../DATA/EH5OM/SO4/RF/T63L31_skg_'//             &
                     & yr_rf(nyear-1941)//'.nc'
            inquire (file=fnso4_rf,exist=there)
            if ( .not.there ) then
              write (*,*) fnso4_rf , ' is not available'
              stop
            end if
            istatus = nf90_open(fnso4_rf,nf90_nowrite,ncid)
            istatus = nf90_get_var(ncid,10,sulfate,istart,icount)
            istatus = nf90_close(ncid)
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
            fnso4_a2 = '../DATA/EH5OM/SO4/A2/T63L31_skg_A2_'//          &
                     & yr_a2(nyear-2000)//'.nc'
            inquire (file=fnso4_a2,exist=there)
            if ( .not.there ) then
              write (*,*) fnso4_a2 , ' is not available'
              stop
            end if
            istatus = nf90_open(fnso4_a2,nf90_nowrite,ncid)
          else if ( ssttyp=='EHA1B' ) then
            fnso4_a1b = '../DATA/EH5OM/SO4/A1B/T63L31_skg_A1B_'//       &
                      & yr_a2(nyear-2000)//'.nc'
            inquire (file=fnso4_a1b,exist=there)
            if ( .not.there ) then
              write (*,*) fnso4_a1b , ' is not available'
              stop
            end if
            istatus = nf90_open(fnso4_a1b,nf90_nowrite,ncid)
          else if ( ssttyp=='EH5B1' ) then
            fnso4_b1 = '../DATA/EH5OM/SO4/B1/T63L31_skg_B1_'//          &
                     & yr_a2(nyear-2000)//'.nc'
            inquire (file=fnso4_b1,exist=there)
            if ( .not.there ) then
              write (*,*) fnso4_b1 , ' is not available'
              stop
            end if
            istatus = nf90_open(fnso4_b1,nf90_nowrite,ncid)
          else
          end if
          istatus = nf90_get_var(ncid,10,sulfate,istart,icount)
          istatus = nf90_close(ncid)
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
              fnso4_a2 = '../DATA/EH5OM/SO4/A2/T63L31_skg_A2_'//        &
                       & yr_a2(nyear-1999)//'.nc'
              inquire (file=fnso4_a2,exist=there)
              if ( .not.there ) then
                write (*,*) fnso4_a2 , ' is not available'
                stop
              end if
              istatus = nf90_open(fnso4_a2,nf90_nowrite,ncid)
            else if ( ssttyp=='EHA1B' ) then
              fnso4_a1b = '../DATA/EH5OM/SO4/A1B/T63L31_skg_A1B_'//     &
                        & yr_a2(nyear-1999)//'.nc'
              inquire (file=fnso4_a1b,exist=there)
              if ( .not.there ) then
                write (*,*) fnso4_a1b , ' is not available'
                stop
              end if
              istatus = nf90_open(fnso4_a1b,nf90_nowrite,ncid)
            else if ( ssttyp=='EH5B1' ) then
              fnso4_b1 = '../DATA/EH5OM/SO4/B1/T63L31_skg_B1_'//        &
                       & yr_a2(nyear-1999)//'.nc'
              inquire (file=fnso4_b1,exist=there)
              if ( .not.there ) then
                write (*,*) fnso4_b1 , ' is not available'
                stop
              end if
              istatus = nf90_open(fnso4_b1,nf90_nowrite,ncid)
            else
            end if
            istatus = nf90_get_var(ncid,10,sulfate,istart,icount)
            istatus = nf90_close(ncid)
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
              fnso4_a2 = '../DATA/EH5OM/SO4/A2/T63L31_skg_A2_'//        &
                       & yr_a2(nyear-2001)//'.nc'
              inquire (file=fnso4_a2,exist=there)
              if ( .not.there ) then
                write (*,*) fnso4_a2 , ' is not available'
                stop
              end if
              istatus = nf90_open(fnso4_a2,nf90_nowrite,ncid)
            else if ( ssttyp=='EHA1B' ) then
              fnso4_a1b = '../DATA/EH5OM/SO4/A1B/T63L31_skg_A1B_'//     &
                        & yr_a2(nyear-2001)//'.nc'
              inquire (file=fnso4_a1b,exist=there)
              if ( .not.there ) then
                write (*,*) fnso4_a1b , ' is not available'
                stop
              end if
              istatus = nf90_open(fnso4_a1b,nf90_nowrite,ncid)
            else if ( ssttyp=='EH5B1' ) then
              fnso4_b1 = '../DATA/EH5OM/SO4/B1/T63L31_skg_B1_'//        &
                       & yr_a2(nyear-2001)//'.nc'
              inquire (file=fnso4_b1,exist=there)
              if ( .not.there ) then
                write (*,*) fnso4_b1 , ' is not available'
                stop
              end if
              istatus = nf90_open(fnso4_b1,nf90_nowrite,ncid)
            else
            end if
            istatus = nf90_get_var(ncid,10,sulfate,istart,icount)
            istatus = nf90_close(ncid)
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
        call bilinx(sulfate3,sulfate2,xlon,xlat,glon,glat,ilon,jlat,jx, &
                  & iy,mlev)
        call bilinx(pso4_3,pso4_2,xlon,xlat,glon,glat,ilon,jlat,jx,iy,1)
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
        do k = kz , 1 , -1
          noutrec = noutrec + 1
          write (64,rec=noutrec) ((sulfate4(j,i,k),j=1,jx),i=1,iy)
        end do
      end if
      end subroutine geteh5om

      end module mod_eh5om

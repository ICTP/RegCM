!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      module mod_fvgcm
      use mod_regcm_param , only : ix , jx , kx , ibyte , dattyp
      use mod_preproc_param

      implicit none

      private

      integer , parameter :: nlev2 = 18 , nlat2 = 181 , nlon2 = 288

      real , dimension(nlev2+1) :: ak , bk
      real , dimension(nlev2) :: pplev , sigma1 , sigmar
      real , dimension(nlat2) :: vlat
      real , dimension(nlon2) :: vlon

      real , target , dimension(nlon2,nlat2,nlev2*4+1) :: bb
      real , target , dimension(nlon2,nlat2,nlev2*3) :: b2
      real , target , dimension(nlon2,nlat2,nlev2*2) :: d2
      real , target , dimension(jx,ix,nlev2*3) :: b3
      real , target , dimension(jx,ix,nlev2*2) :: d3

      real , dimension(nlon2,nlat2) :: zs2
      real , dimension(nlon2,nlat2,nlev2) :: pp3d , z1

      real , dimension(jx,ix) :: b3pd
      real , dimension(jx,ix,nlev2) :: w3

      real , pointer , dimension(:,:) :: ps2
      real , pointer , dimension(:,:,:) :: q2 , t2 , u2 , v2
      real , pointer , dimension(:,:,:) :: tp , qp , hp
      real , pointer , dimension(:,:,:) :: up , vp
      real , pointer , dimension(:,:,:) :: t3 , q3 , h3
      real , pointer , dimension(:,:,:) :: u3 , v3

      public :: getfvgcm , headerfv

      contains

      subroutine getfvgcm(idate)
      use mod_grid
      use mod_write
      implicit none
!
! Dummy arguments
!
      integer :: idate
!
! Local variables
!
      character(3) , dimension(12) :: chmon
      character(20) :: finm , fips
      character(5) :: fn_a2 , fn_rf , pn_a2 , pn_rf
      integer :: i , i2 , ii , j , j2 , k , month , mrec , nday ,       &
               & nhour , nmop , nrec , numx , numy , nyear , nyrp
      integer(2) , dimension(288,181) :: itmp
      real(8) :: offset , xscale
      real , dimension(288,181) :: temp
      logical :: there
      real :: wt
      character(4) , dimension(30) :: yr_a2 , yr_rf
!
      data fn_rf/'FV_RF'/ , fn_a2/'FV_A2'/
      data pn_rf/'PS_RF'/ , pn_a2/'PS_A2'/
      data yr_rf/'1961' , '1962' , '1963' , '1964' , '1965' , '1966' ,  &
          &'1967' , '1968' , '1969' , '1970' , '1971' , '1972' ,        &
         & '1973' , '1974' , '1975' , '1976' , '1977' , '1978' ,        &
         & '1979' , '1980' , '1981' , '1982' , '1983' , '1984' ,        &
         & '1985' , '1986' , '1987' , '1988' , '1989' , '1990'/
      data yr_a2/'2071' , '2072' , '2073' , '2074' , '2075' , '2076' ,  &
          &'2077' , '2078' , '2079' , '2080' , '2081' , '2082' ,        &
         & '2083' , '2084' , '2085' , '2086' , '2087' , '2088' ,        &
         & '2089' , '2090' , '2091' , '2092' , '2093' , '2094' ,        &
         & '2095' , '2096' , '2097' , '2098' , '2099' , '2100'/
      data chmon/'JAN' , 'FEB' , 'MAR' , 'APR' , 'MAY' , 'JUN' , 'JUL' ,&
          &'AUG' , 'SEP' , 'OCT' , 'NOV' , 'DEC'/
!
!
      if ( idate==idate1 ) then
        numx = nint((lon1-lon0)/1.25) + 1
        numy = nint(lat1-lat0) + 1
        inquire (file='../DATA/FVGCM/HT_SRF',exist=there)
        if ( .not.there ) then
          write (*,*) '../DATA/FVGCM/HT_SRF is not available'
          stop
        end if
        open (61,file='../DATA/FVGCM/HT_SRF',form='unformatted',        &
            & recl=numx*numy*ibyte,access='direct')
        read (61,rec=1) ((temp(i,j),i=1,numx),j=1,numy)
        do j = nint(lat0) , nint(lat1)
          do i = nint(lon0/1.25) , nint(lon1/1.25)
            ii = i + 1
            if ( ii<=0 ) ii = ii + 288
            if ( ii>288 ) ii = ii - 288
            i2 = i - nint(lon0/1.25) + 1
            j2 = j - nint(lat0) + 1
            zs2(ii,j+91) = temp(i2,j2)/9.80616
          end do
        end do
        close (61)
      end if
 
      nyear = idate/1000000
      month = idate/10000 - nyear*100
      nday = idate/100 - (idate/10000)*100
      nhour = mod(idate,100)
 
      if ( nday/=1 .or. nhour/=0 ) then
        if ( ssttyp=='FV_RF' ) then
          finm = 'RF/'//yr_rf(nyear-1960)//'/'//fn_rf//yr_rf(nyear-1960)&
               & //chmon(month)
          fips = 'RF/'//yr_rf(nyear-1960)//'/'//pn_rf//yr_rf(nyear-1960)&
               & //chmon(month)
        else if ( ssttyp=='FV_A2' ) then
          finm = 'A2/'//yr_a2(nyear-2070)//'/'//fn_a2//yr_a2(nyear-2070)&
               & //chmon(month)
          fips = 'A2/'//yr_a2(nyear-2070)//'/'//pn_a2//yr_a2(nyear-2070)&
               & //chmon(month)
        else
          write (*,*) 'ERROR IN SSTTYP'
          stop
        end if
      else if ( month/=1 ) then
        if ( ssttyp=='FV_RF' ) then
          finm = 'RF/'//yr_rf(nyear-1960)//'/'//fn_rf//yr_rf(nyear-1960)&
               & //chmon(month-1)
          fips = 'RF/'//yr_rf(nyear-1960)//'/'//pn_rf//yr_rf(nyear-1960)&
               & //chmon(month-1)
        else if ( ssttyp=='FV_A2' ) then
          finm = 'A2/'//yr_a2(nyear-2070)//'/'//fn_a2//yr_a2(nyear-2070)&
               & //chmon(month-1)
          fips = 'A2/'//yr_a2(nyear-2070)//'/'//pn_a2//yr_a2(nyear-2070)&
               & //chmon(month-1)
        else
          write (*,*) 'ERROR IN SSTTYP'
          stop
        end if
      else if ( ssttyp=='FV_RF' ) then
        if ( nyear==1961 ) then
          write (*,*) 'Fields on 00z01jan1961 is not saved'
          write (*,*) 'Please run from 00z02jan1961'
          stop
        end if
        finm = 'RF/'//yr_rf(nyear-1961)//'/'//fn_rf//yr_rf(nyear-1961)  &
             & //chmon(12)
        fips = 'RF/'//yr_rf(nyear-1961)//'/'//pn_rf//yr_rf(nyear-1961)  &
             & //chmon(12)
      else if ( ssttyp=='FV_A2' ) then
        if ( nyear==2071 ) then
          write (*,*) 'Fields on 00z01jan2071 is not saved'
          write (*,*) 'Please run from 00z02jan2071'
          stop
        end if
        finm = 'A2/'//yr_a2(nyear-2071)//'/'//fn_a2//yr_a2(nyear-2071)  &
             & //chmon(12)
        fips = 'A2/'//yr_a2(nyear-2071)//'/'//pn_a2//yr_a2(nyear-2071)  &
             & //chmon(12)
      else
        write (*,*) 'ERROR IN SSTTYP'
        stop
      end if
      numx = nint((lon1-lon0)/1.25) + 1
      numy = nint(lat1-lat0) + 1
      do k = 1 , nlev2*4 + 1
        do j = 1 , nlat2
          do i = 1 , nlon2
            bb(i,j,k) = -9999.
          end do
        end do
      end do
      inquire (file='../DATA/FVGCM/'//finm,exist=there)
      if ( .not.there ) then
        write (*,*) '../DATA/FVGCM/'//finm , ' is not available'
        write (*,*) 'please copy FVGCM output under ../DATA/FVGCM/'
        stop
      end if
      open (63,file='../DATA/FVGCM/'//finm,form='unformatted',          &
          & recl=(numx*numy*2+16)/4*ibyte,access='direct')
      open (62,file='../DATA/FVGCM/'//fips,form='unformatted',          &
          & recl=numx*numy*ibyte,access='direct')
      if ( nday/=1 .or. nhour/=0 ) then
        nrec = ((nday-1)*4+nhour/6-1)*(nlev2*4)
        mrec = (nday-1)*4 + nhour/6 - 1
      else if ( month==1 .or. month==2 .or. month==4 .or. month==6 .or. &
              & month==8 .or. month==9 .or. month==11 ) then
        nrec = (31*4-1)*(nlev2*4)
        mrec = 31*4 - 1
      else if ( month==5 .or. month==7 .or. month==10 .or. month==12 )  &
              & then
        nrec = (30*4-1)*(nlev2*4)
        mrec = 30*4 - 1
      else if ( mod(nyear,4)==0 .and. nyear/=2100 ) then
        nrec = (29*4-1)*(nlev2*4)
        mrec = 29*4 - 1
      else
        nrec = (28*4-1)*(nlev2*4)
        mrec = 28*4 - 1
      end if
      mrec = mrec + 1
      read (62,rec=mrec) ((temp(i,j),i=1,numx),j=1,numy)
      do j = nint(lat0) , nint(lat1)
        do i = nint(lon0/1.25) , nint(lon1/1.25)
          ii = i + 1
          if ( ii<=0 ) ii = ii + 288
          if ( ii>288 ) ii = ii - 288
          i2 = i - nint(lon0/1.25) + 1
          j2 = j - nint(lat0) + 1
          ps2(ii,j+91) = temp(i2,j2)*0.01
        end do
      end do
      do k = 1 , nlev2
        nrec = nrec + 1
        read (63,rec=nrec) offset , xscale ,                             &
                         & ((itmp(i,j),i=1,numx),j=1,numy)
        do j = nint(lat0) , nint(lat1)
          do i = nint(lon0/1.25) , nint(lon1/1.25)
            ii = i + 1
            if ( ii<=0 ) ii = ii + 288
            if ( ii>288 ) ii = ii - 288
            i2 = i - nint(lon0/1.25) + 1
            j2 = j - nint(lat0) + 1
            u2(ii,j+91,k) = itmp(i2,j2)*xscale + offset
          end do
        end do
      end do
      do k = 1 , nlev2
        nrec = nrec + 1
        read (63,rec=nrec) offset , xscale ,                             &
                         & ((itmp(i,j),i=1,numx),j=1,numy)
        do j = nint(lat0) , nint(lat1)
          do i = nint(lon0/1.25) , nint(lon1/1.25)
            ii = i + 1
            if ( ii<=0 ) ii = ii + 288
            if ( ii>288 ) ii = ii - 288
            i2 = i - nint(lon0/1.25) + 1
            j2 = j - nint(lat0) + 1
            v2(ii,j+91,k) = itmp(i2,j2)*xscale + offset
          end do
        end do
      end do
      do k = 1 , nlev2
        nrec = nrec + 1
        read (63,rec=nrec) offset , xscale ,                             &
                         & ((itmp(i,j),i=1,numx),j=1,numy)
        do j = nint(lat0) , nint(lat1)
          do i = nint(lon0/1.25) , nint(lon1/1.25)
            ii = i + 1
            if ( ii<=0 ) ii = ii + 288
            if ( ii>288 ) ii = ii - 288
            i2 = i - nint(lon0/1.25) + 1
            j2 = j - nint(lat0) + 1
            t2(ii,j+91,k) = itmp(i2,j2)*xscale + offset
          end do
        end do
      end do
      do k = 1 , nlev2
        nrec = nrec + 1
        read (63,rec=nrec) offset , xscale ,                             &
                         & ((itmp(i,j),i=1,numx),j=1,numy)
        do j = nint(lat0) , nint(lat1)
          do i = nint(lon0/1.25) , nint(lon1/1.25)
            ii = i + 1
            if ( ii<=0 ) ii = ii + 288
            if ( ii>288 ) ii = ii - 288
            i2 = i - nint(lon0/1.25) + 1
            j2 = j - nint(lat0) + 1
            q2(ii,j+91,k) = itmp(i2,j2)*xscale + offset
          end do
        end do
      end do
      close (63)
      close (62)
      write (*,*) 'READ IN fields at DATE:' , idate , ' from ' , finm
      do k = 1 , nlev2
        do j = 1 , nlat2
          do i = 1 , nlon2
            if ( ps2(i,j)>-9995. ) then
              pp3d(i,j,k) = ps2(i,j)*0.5*(bk(k)+bk(k+1))                &
                          & + 0.5*(ak(k)+ak(k+1))
            else
              pp3d(i,j,k) = -9999.0
            end if
          end do
        end do
      end do
!
!     to calculate Heights on sigma surfaces.
      call htsig(t2,z1,pp3d,ps2,zs2,nlon2,nlat2,nlev2)
!
!     to interpolate H,U,V,T,Q and QC
!     1. For Heights
      call height(hp,z1,t2,ps2,pp3d,zs2,nlon2,nlat2,nlev2,pplev,nlev2)
!     2. For Zonal and Meridional Winds
      call intlin(up,u2,ps2,pp3d,nlon2,nlat2,nlev2,pplev,nlev2)
      call intlin(vp,v2,ps2,pp3d,nlon2,nlat2,nlev2,pplev,nlev2)
!     3. For Temperatures
      call intlog(tp,t2,ps2,pp3d,nlon2,nlat2,nlev2,pplev,nlev2)
!     4. For Moisture qva & qca
      call humid1fv(t2,q2,pp3d,nlon2,nlat2,nlev2)
      call intlin(qp,q2,ps2,pp3d,nlon2,nlat2,nlev2,pplev,nlev2)
!
!     HORIZONTAL INTERPOLATION OF BOTH THE SCALAR AND VECTOR FIELDS
!
      call bilinx2(b3,b2,xlon,xlat,vlon,vlat,nlon2,nlat2,jx,ix,nlev2*3)
      call bilinx2(d3,d2,dlon,dlat,vlon,vlat,nlon2,nlat2,jx,ix,nlev2*2)
!
!     ROTATE U-V FIELDS AFTER HORIZONTAL INTERPOLATION
!
      call uvrot4(u3,v3,dlon,dlat,clon,clat,grdfac,jx,ix,nlev2,plon,    &
                & plat,iproj)
!
!     X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!     X X
!     V E R T I C A L   I N T E R P O L A T I O N
!
!     X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X X
!     X X
!HH:  CHANGE THE VERTICAL ORDER.
      call top2btm(t3,jx,ix,nlev2)
      call top2btm(q3,jx,ix,nlev2)
      call top2btm(h3,jx,ix,nlev2)
      call top2btm(u3,jx,ix,nlev2)
      call top2btm(v3,jx,ix,nlev2)
!HH:OVER
!
!     ******           NEW CALCULATION OF P* ON RegCM TOPOGRAPHY.
      call intgtb(pa,za,tlayer,topogm,t3,h3,sigmar,jx,ix,nlev2)
 
      call intpsn(ps4,topogm,pa,za,tlayer,ptop,jx,ix)
      call p1p2(b3pd,ps4,jx,ix)
!
!     F0    DETERMINE SURFACE TEMPS ON RegCM TOPOGRAPHY.
!     INTERPOLATION FROM PRESSURE LEVELS AS IN INTV2
      call intv3(ts4,t3,ps4,sigmar,ptop,jx,ix,nlev2)
 
      if ( ssttyp/='OI_WK' .and. ssttyp/='OI2WK' ) then
!       F1    CALCULATE SSTS FOR DATE FROM OBSERVED SSTS
!       PRINT *, 'INPUT DAY FOR SST DATA ACQUISITION:', IDATE
        call julian(idate,nyrp,nmop,wt)
!
        if ( ssttyp=='OI2ST' ) then
          call mkssta(ts4,sst1,sst2,ice1,ice2,topogm,xlandu,jx,ix,nyrp, &
               &      nmop,wt)
        else
          call mksst(ts4,sst1,sst2,topogm,xlandu,jx,ix,nyrp,nmop,wt)
        end if
      else
        if ( ssttyp=='OI2WK' ) then
          call mksst2(ts4,sst1,sst2,ice1,ice2,topogm,xlandu,jx,ix,      &
               &      idate/100)
        else
          call mksst2(ts4,sst1,sst2,topogm,xlandu,jx,ix,idate/100)
        end if
      end if
 
!     F2     DETERMINE P* AND HEIGHT.
!
!     F3     INTERPOLATE U, V, T, AND Q.
      call intv1(u4,u3,b3pd,sigma2,sigmar,ptop,jx,ix,kx,nlev2)
      call intv1(v4,v3,b3pd,sigma2,sigmar,ptop,jx,ix,kx,nlev2)
!
      call intv2(t4,t3,ps4,sigma2,sigmar,ptop,jx,ix,kx,nlev2)
 
      call intv1(q4,q3,ps4,sigma2,sigmar,ptop,jx,ix,kx,nlev2)
      call humid2fv(t4,q4,ps4,ptop,sigma2,jx,ix,kx)
!
!     F4     DETERMINE H
      call hydrost(h4,t4,topogm,ps4,ptop,sigmaf,sigma2,dsigma,jx,ix,kx)
!
!     G      WRITE AN INITIAL FILE FOR THE RegCM
      call writef(ptop,idate)
!
      end subroutine getfvgcm

      subroutine headerfv
      implicit none
!
! Local variables
!
      integer :: i , j , k , kr
!
      do j = 1 , nlat2
        vlat(j) = float(j-1) - 90.0
      end do
!
      do i = 1 , nlon2
        vlon(i) = float(i-1)*1.25
      end do
 
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
 
      do k = 1 , nlev2
        sigmar(k) = pplev(k)*0.001
      end do
!HH:OVER
!     CHANGE ORDER OF VERTICAL INDEXES FOR PRESSURE LEVELS
!
      do k = 1 , nlev2
        kr = nlev2 - k + 1
        sigma1(k) = sigmar(kr)
      end do
!
      ak(1) = 2.9170
      ak(2) = 7.9292
      ak(3) = 21.5539
      ak(4) = 49.1834
      ak(5) = 83.1425
      ak(6) = 79.9308
      ak(7) = 75.7738
      ak(8) = 70.5752
      ak(9) = 64.2963
      ak(10) = 56.9838
      ak(11) = 48.7913
      ak(12) = 39.9895
      ak(13) = 30.9631
      ak(14) = 22.1902
      ak(15) = 14.2039
      ak(16) = 7.5413
      ak(17) = 2.6838
      ak(18) = 0.0
      ak(19) = 0.0
 
      bk(1) = 0.0
      bk(2) = 0.0
      bk(3) = 0.0
      bk(4) = 0.0
      bk(5) = 0.0
      bk(6) = 0.0380541
      bk(7) = 0.0873088
      bk(8) = 0.1489307
      bk(9) = 0.2232996
      bk(10) = 0.3099406
      bk(11) = 0.4070096
      bk(12) = 0.5112977
      bk(13) = 0.6182465
      bk(14) = 0.7221927
      bk(15) = 0.8168173
      bk(16) = 0.8957590
      bk(17) = 0.9533137
      bk(18) = 0.9851222
      bk(19) = 1.0
 
!     Set up pointers

      ps2 => bb(:,:,1)
      t2 => bb(:,:,2:nlev2+1)
      q2 => bb(:,:,nlev2+2:2*nlev2+1)
      u2 => bb(:,:,2*nlev2+2:3*nlev2+1)
      v2 => bb(:,:,3*nlev2+2:4*nlev2+1)
      tp => b2(:,:,1:nlev2)
      qp => b2(:,:,nlev2+1:2*nlev2)
      hp => b2(:,:,2*nlev2+1:3*nlev2)
      up => d2(:,:,1:nlev2)
      vp => d2(:,:,nlev2+1:2*nlev2)
      t3 => b3(:,:,1:nlev2)
      q3 => b3(:,:,nlev2+1:2*nlev2)
      h3 => b3(:,:,2*nlev2+1:3*nlev2)
      u3 => d3(:,:,1:nlev2)
      v3 => d3(:,:,nlev2+1:2*nlev2)

      end subroutine headerfv

      end module mod_fvgcm

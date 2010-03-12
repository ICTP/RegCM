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

      module mod_gfs11
      use mod_param
      implicit none

      integer , parameter :: klev = 26 , jlat = 181 , ilon = 360

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

      contains

      subroutine getgfs11(idate)
      use mod_grid
      use mod_write
      implicit none
!
! PARAMETER definitions
!
      real , parameter :: ptop = 5.0 , clat = 50.00 , clon = 13.00 ,    &
                        & plat = clat , plon = clon
      integer , parameter :: igrads = 1 , ibigend = 1
!
! Dummy arguments
!
      integer :: idate
!
! Local variables
!
      real , dimension(jx,iy) :: b3pd , pa , sst1 , sst2 , tlayer , za ,&
                    &            ice1 , ice2
      character(3) , dimension(12) :: chmon
      character(17) :: finm
      integer :: i , i2 , ii , j , j2 , k , month , nday , nhour ,      &
               & nmop , nrec , numx , numy , nyear , nyrp
      integer(2) , dimension(360,181) :: itmp
      real(8) :: offset , xscale
      logical :: there
      real :: wt
      character(4) , dimension(9) :: yrgfs
!
!     DIMENSION SURFACE TEMPERATURE ON RCM SURFACE; NOT GIVEN BY ECMWF
!     DATA READ FROM THE OUTPUT OF SST STEP
!
      data yrgfs/'2000' , '2001' , '2002' , '2003' , '2004' , '2005' ,  &
          &'2006' , '2007' , '2008'/
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
      month = idate/10000 - nyear*100
      nday = idate/100 - (idate/10000)*100
      nhour = mod(idate,100)
      if ( idate<2000010106 ) then
        write (*,*) 'GFS 1x1 datasets is just avaiable from 2000010106'
        stop
      end if
!     IF(IDATE.GT.2008060100) THEN
!     WRITE(*,*) 'GFS 1x1 datasets is just avaiable to 2008060100'
!     STOP
!     ENDIF
      numx = nint((lon1-lon0)/1.0) + 1
      numy = nint((lat1-lat0)/1.0) + 1
      if ( numx/=360 .or. numy/=181 ) then
        if ( nday/=1 .or. nhour/=0 ) then
          finm = yrgfs(nyear-1999)//'/'//'GFS__'//yrgfs(nyear-1999)     &
               & //chmon(month)
        else if ( month/=1 ) then
          finm = yrgfs(nyear-1999)//'/'//'GFS__'//yrgfs(nyear-1999)     &
               & //chmon(month-1)
        else
          finm = yrgfs(nyear-2000)//'/'//'GFS__'//yrgfs(nyear-2000)     &
               & //chmon(12)
        end if
      else if ( nday/=1 .or. nhour/=0 ) then
        finm = yrgfs(nyear-1999)//'/'//'GFSg_'//yrgfs(nyear-1999)       &
             & //chmon(month)
      else if ( month/=1 ) then
        finm = yrgfs(nyear-1999)//'/'//'GFSg_'//yrgfs(nyear-1999)       &
             & //chmon(month-1)
      else
        finm = yrgfs(nyear-2000)//'/'//'GFSg_'//yrgfs(nyear-2000)       &
             & //chmon(12)
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
      inquire (file='../DATA/GFS11/'//finm,exist=there)
      if ( .not.there ) then
        write (*,*) '../DATA/GFS11/'//finm , ' is not available'
        write (*,*) 'please copy GFS11 datasets under ../DATA/GFS11/'
        stop
      end if
      open (63,file='../DATA/GFS11/'//finm,form='unformatted',          &
          & recl=(numx*numy*2+16)/4*ibyte,access='direct')
      if ( nday/=1 .or. nhour/=0 ) then
        nrec = ((nday-1)*4+nhour/6-1)*127
      else if ( month==1 .or. month==2 .or. month==4 .or. month==6 .or. &
              & month==8 .or. month==9 .or. month==11 ) then
        nrec = (31*4-1)*127
      else if ( month==5 .or. month==7 .or. month==10 .or. month==12 )  &
              & then
        nrec = (30*4-1)*127
      else if ( mod(nyear,4)==0 ) then
        nrec = (29*4-1)*127
      else
        nrec = (28*4-1)*127
      end if
      do k = klev , 1 , -1
        nrec = nrec + 1
        read (63,rec=nrec) offset , xscale ,                             &
                         & ((itmp(i,j),i=1,numx),j=1,numy)
        do j = nint(lat0/1.0) , nint(lat1/1.0)
          do i = nint(lon0/1.0) , nint(lon1/1.0)
            ii = i + 1
            if ( ii<=0 ) ii = ii + 360
            if ( ii>360 ) ii = ii - 360
            i2 = i - nint(lon0/1.0) + 1
            j2 = j - nint(lat0/1.0) + 1
            if ( numx==360 .and. numy==181 ) then
              hvar(ii,91-j,k) = itmp(i2,j2)*xscale + offset
            else
              hvar(ii,j+90,k) = itmp(i2,j2)*xscale + offset
            end if
          end do
        end do
      end do
      do k = klev , 1 , -1
        nrec = nrec + 1
        read (63,rec=nrec) offset , xscale ,                             &
                         & ((itmp(i,j),i=1,numx),j=1,numy)
        do j = nint(lat0/1.0) , nint(lat1/1.0)
          do i = nint(lon0/1.0) , nint(lon1/1.0)
            ii = i + 1
            if ( ii<=0 ) ii = ii + 360
            if ( ii>360 ) ii = ii - 360
            i2 = i - nint(lon0/1.0) + 1
            j2 = j - nint(lat0/1.0) + 1
            if ( numx==360 .and. numy==181 ) then
              tvar(ii,91-j,k) = itmp(i2,j2)*xscale + offset
            else
              tvar(ii,j+90,k) = itmp(i2,j2)*xscale + offset
            end if
          end do
        end do
      end do
      do k = klev , 6 , -1
        nrec = nrec + 1
        read (63,rec=nrec) offset , xscale ,                             &
                         & ((itmp(i,j),i=1,numx),j=1,numy)
        do j = nint(lat0/1.0) , nint(lat1/1.0)
          do i = nint(lon0/1.0) , nint(lon1/1.0)
            ii = i + 1
            if ( ii<=0 ) ii = ii + 360
            if ( ii>360 ) ii = ii - 360
            i2 = i - nint(lon0/1.0) + 1
            j2 = j - nint(lat0/1.0) + 1
            if ( numx==360 .and. numy==181 ) then
              rhvar(ii,91-j,k) = itmp(i2,j2)*xscale + offset
            else
              rhvar(ii,j+90,k) = itmp(i2,j2)*xscale + offset
            end if
          end do
        end do
      end do
      do k = 5 , 1 , -1
        do j = 1 , jlat
          do i = 1 , ilon
            rhvar(i,j,k) = rhvar(i,j,k+1)
          end do
        end do
      end do
      do k = 1 , klev
        do j = 1 , jlat
          do i = 1 , ilon
            rhvar(i,j,k) = amax1(rhvar(i,j,k)*0.01,1.0)
          end do
        end do
      end do
      do k = klev , 1 , -1
        nrec = nrec + 1
        read (63,rec=nrec) offset , xscale ,                             &
                         & ((itmp(i,j),i=1,numx),j=1,numy)
        do j = nint(lat0/1.0) , nint(lat1/1.0)
          do i = nint(lon0/1.0) , nint(lon1/1.0)
            ii = i + 1
            if ( ii<=0 ) ii = ii + 360
            if ( ii>360 ) ii = ii - 360
            i2 = i - nint(lon0/1.0) + 1
            j2 = j - nint(lat0/1.0) + 1
            if ( numx==360 .and. numy==181 ) then
              uvar(ii,91-j,k) = itmp(i2,j2)*xscale + offset
            else
              uvar(ii,j+90,k) = itmp(i2,j2)*xscale + offset
            end if
          end do
        end do
      end do
      do k = klev , 1 , -1
        nrec = nrec + 1
        read (63,rec=nrec) offset , xscale ,                             &
                         & ((itmp(i,j),i=1,numx),j=1,numy)
        do j = nint(lat0/1.0) , nint(lat1/1.0)
          do i = nint(lon0/1.0) , nint(lon1/1.0)
            ii = i + 1
            if ( ii<=0 ) ii = ii + 360
            if ( ii>360 ) ii = ii - 360
            i2 = i - nint(lon0/1.0) + 1
            j2 = j - nint(lat0/1.0) + 1
            if ( numx==360 .and. numy==181 ) then
              vvar(ii,91-j,k) = itmp(i2,j2)*xscale + offset
            else
              vvar(ii,j+90,k) = itmp(i2,j2)*xscale + offset
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
      call p1p2(b3pd,ps4,jx,iy)
 
!
!     F0  DETERMINE SURFACE TEMPS ON RCM TOPOGRAPHY.
!     INTERPOLATION FROM PRESSURE LEVELS AS IN INTV2
      call intv3(ts4,t3,ps4,sigmar,ptop,jx,iy,klev)
 
      if ( ssttyp/='OI_WK' .and. ssttyp/='OI2WK' ) then
!       F1  CALCULATE SSTS FOR DATE FROM OBSERVED SSTS
        print * , 'INPUT DAY FOR SST DATA ACQUISITION:' , idate
        call julian(idate,nyrp,nmop,wt)
!
        if ( ssttyp=='OI2ST' ) then
          call mkssta(ts4,sst1,sst2,ice1,ice2,topogm,xlandu,jx,iy,nyrp, &
                  &   nmop,wt)
        else
          call mksst(ts4,sst1,sst2,topogm,xlandu,jx,iy,nyrp,nmop,wt)
        end if
      else
        if ( ssttyp=='OI2WK' ) then
          call mksst2a(ts4,sst1,sst2,ice1,ice2,topogm,xlandu,jx,iy,     &
              &        idate/100)
        else
          call mksst2(ts4,sst1,sst2,topogm,xlandu,jx,iy,idate/100)
        end if
      end if
 
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
      call writef(ptop,idate)
!
      end subroutine getgfs11
!
      subroutine headergfs
      implicit none
!
! Local variables
!
      integer :: i , j , k , kr
!
      sigmar(1) = .01
      sigmar(2) = .02
      sigmar(3) = .03
      sigmar(4) = .05
      sigmar(5) = .07
      sigmar(6) = .1
      sigmar(7) = .15
      sigmar(8) = .2
      sigmar(9) = .25
      sigmar(10) = .3
      sigmar(11) = .35
      sigmar(12) = .4
      sigmar(13) = .45
      sigmar(14) = .5
      sigmar(15) = .55
      sigmar(16) = .6
      sigmar(17) = .65
      sigmar(18) = .7
      sigmar(19) = .75
      sigmar(20) = .8
      sigmar(21) = .85
      sigmar(22) = .9
      sigmar(23) = .925
      sigmar(24) = .95
      sigmar(25) = .975
      sigmar(26) = 1.00
!
!     INITIAL GLOBAL GRID-POINT LONGITUDE & LATITUDE
!
      do i = 1 , ilon
        glon(i) = float(i-1)*1.0
      end do
      do j = 1 , jlat
        glat(j) = -90.0 + float(j-1)*1.0
      end do
!HH:OVER
!     CHANGE ORDER OF VERTICAL INDEXES FOR PRESSURE LEVELS
!
      do k = 1 , klev
        kr = klev - k + 1
        sigma1(k) = sigmar(kr)
      end do
 
      end subroutine headergfs

      end module mod_gfs11

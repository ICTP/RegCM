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

      module mod_era40
      use mod_param

      implicit none

      integer , parameter :: klev = 23 , jlat = 73 , ilon = 144

      real , target , dimension(ilon,jlat,klev*3) :: b2
      real , target , dimension(ilon,jlat,klev*2) :: d2
      real , target , dimension(ilon,jlat,4*3+1) :: s2
      real , target , dimension(jx,iy,klev*3) :: b3
      real , target , dimension(jx,iy,klev*2) :: d3
      real , target , dimension(jx,iy,4*3+1) :: s3

      real , dimension(ilon,jlat,klev) :: wvar

      real , pointer :: u3(:,:,:) , v3(:,:,:)
      real , pointer :: h3(:,:,:) , q3(:,:,:) , t3(:,:,:)
      real , pointer :: uvar(:,:,:) , vvar(:,:,:)
      real , pointer :: hvar(:,:,:) , rhvar(:,:,:) , tvar(:,:,:)
      real , pointer , dimension(:,:,:) :: qsoil , tsice , tsoil
      real , pointer , dimension(:,:) :: snw
      real , pointer , dimension(:,:,:) :: qs3 , ti3 , ts3
      real , pointer , dimension(:,:) :: snow

      real , dimension(jlat) :: glat
      real , dimension(ilon) :: glon
      real , dimension(klev) :: sigma1 , sigmar

      contains

      subroutine getera40(idate)
      use mod_grid
      use mod_var4
      use mod_write
      implicit none
!
! Dummy arguments
!
      integer :: idate
!
! Local variables
!
      real , dimension(jx,iy) :: b3pd , pa , sst1 , sst2 , tlayer , za ,&
                          &      ice1 , ice2
      integer :: nmop , nyrp
      real :: wt
!
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
!
!     D      BEGIN LOOP OVER NTIMES
!
      call era6hour(dattyp,lsmtyp,idate,idate1)
      write (*,*) 'READ IN fields at DATE:' , idate
!
!     HORIZONTAL INTERPOLATION OF BOTH THE SCALAR AND VECTOR FIELDS
!
      call bilinx(b3,b2,xlon,xlat,glon,glat,ilon,jlat,jx,iy,klev*3)
      call bilinx(d3,d2,dlon,dlat,glon,glat,ilon,jlat,jx,iy,klev*2)
      if ( lsmtyp=='USGS' ) call bilinx(s3,s2,xlon,xlat,glon,glat,ilon, &
                                      & jlat,jx,iy,4*3+1)
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
 
      if ( ssttyp/='OI_WK' .and. ssttyp/='OI2WK' ) then
!       F1  CALCULATE SSTS FOR DATE FROM OBSERVED SSTS
        print * , 'INPUT DAY FOR SST DATA ACQUISITION:' , idate
        call julian(idate,nyrp,nmop,wt)
!
        if ( ssttyp/='OI2ST' ) then
          call mkssta(ts4,sst1,sst2,ice1,ice2,topogm,xlandu,jx,iy,      &
                 &    nyrp,nmop,wt)
        else
          call mksst(ts4,sst1,sst2,topogm,xlandu,jx,iy,nyrp,nmop,wt)
        end if
      else
        if ( ssttyp/='OI2WK' ) then
          call mksst2a(ts4,sst1,sst2,ice1,ice2,topogm,xlandu,jx,iy,     &
                   &   idate/100)
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
      call writefs(u4,v4,t4,q4,ps4,ts4,qs3,ti3,ts3,snow,ptop,jx,iy,kz,  &
                 & idate,lsmtyp)
!
      end subroutine getera40

      subroutine era6hour(dattyp,lsmtyp,idate,idate0)
      use netcdf
      implicit none
!
! Dummy arguments
!
      character(5) :: dattyp
      integer :: idate , idate0
      character(4) :: lsmtyp
      intent (in) dattyp , idate , idate0 , lsmtyp
!
! Local variables
!
      integer :: i , inet , it , j , k , k4 , kkrec , l4 , month ,      &
               & nday , nhour , nyear , istatus
      character(24) :: inname
      character(64) :: pathaddname
      character(5) , dimension(3,4) :: sarname
      character(2) :: snownm
      logical :: there
      character(5) , dimension(6) :: varname
      integer(2) , dimension(ilon,jlat,klev) :: work
      integer(2) , dimension(ilon,jlat) :: work2d
      real(8) :: xadd , xscale
      integer , dimension(10) :: icount , istart
      real(8) , dimension(5,4) :: xoff , xscl
      real(8) , dimension(3,4,4) :: xoff_s , xscl_s
      real(8) , dimension(4) :: xoff_sn , xscl_sn
      integer , dimension(5,4) :: inet6
      integer , dimension(3,4,4) :: isnet3
      integer , dimension(4) :: isnow
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
!bxq
!bxq_
      data varname/'air' , 'hgt' , 'rhum' , 'uwnd' , 'vwnd' , 'omega'/
      data sarname/'swvl1' , 'istl1' , 'stl1' , 'swvl2' , 'istl2' ,     &
          &'stl2' , 'swvl3' , 'istl3' , 'stl3' , 'swvl4' , 'istl4' ,    &
          &'stl4'/
      data snownm/'sd'/
!
!     Below in the ncopen call is the file name of the netCDF file.
!     You may want to add code to read in the file name and the
!     variable name.
!     OPEN FILE AND GET FILES ID AND VARIABLE ID(S)
!
!bxq
      if ( idate<1957090100 .or. idate>2002083118 ) then
        write (*,*) 'ERA40 datasets is just available from' ,           &
                   &' 1957090100 to 2002083118'
        stop
      end if
 
      nyear = idate/1000000
      month = idate/10000 - nyear*100
      nday = idate/100 - nyear*10000 - month*100
      nhour = idate - nyear*1000000 - month*10000 - nday*100
      if ( idate==idate0 .or.                                           &
         & (mod(idate,100000)==10100 .and. mod(idate,1000000)/=110100) )&
         & then
        do k4 = 1 , 4
          do kkrec = 1 , 5
            if ( kkrec==1 ) then
              if ( k4==1 ) then
                write (inname,99001) nyear , 'air.' , nyear
              else if ( k4==2 ) then
                write (inname,99002) nyear , 'air.' , nyear
              else if ( k4==3 ) then
                write (inname,99003) nyear , 'air.' , nyear
              else if ( k4==4 ) then
                write (inname,99004) nyear , 'air.' , nyear
              else
              end if
            else if ( kkrec==2 ) then
              if ( k4==1 ) then
                write (inname,99001) nyear , 'hgt.' , nyear
              else if ( k4==2 ) then
                write (inname,99002) nyear , 'hgt.' , nyear
              else if ( k4==3 ) then
                write (inname,99003) nyear , 'hgt.' , nyear
              else if ( k4==4 ) then
                write (inname,99004) nyear , 'hgt.' , nyear
              else
              end if
            else if ( kkrec==3 ) then
              if ( k4==1 ) then
                write (inname,99005) nyear , 'rhum.' , nyear
              else if ( k4==2 ) then
                write (inname,99006) nyear , 'rhum.' , nyear
              else if ( k4==3 ) then
                write (inname,99007) nyear , 'rhum.' , nyear
              else if ( k4==4 ) then
                write (inname,99008) nyear , 'rhum.' , nyear
              else
              end if
            else if ( kkrec==4 ) then
              if ( k4==1 ) then
                write (inname,99005) nyear , 'uwnd.' , nyear
              else if ( k4==2 ) then
                write (inname,99006) nyear , 'uwnd.' , nyear
              else if ( k4==3 ) then
                write (inname,99007) nyear , 'uwnd.' , nyear
              else if ( k4==4 ) then
                write (inname,99008) nyear , 'uwnd.' , nyear
              else
              end if
            else if ( kkrec==5 ) then
              if ( k4==1 ) then
                write (inname,99005) nyear , 'vwnd.' , nyear
              else if ( k4==2 ) then
                write (inname,99006) nyear , 'vwnd.' , nyear
              else if ( k4==3 ) then
                write (inname,99007) nyear , 'vwnd.' , nyear
              else if ( k4==4 ) then
                write (inname,99008) nyear , 'vwnd.' , nyear
              else
              end if
            else if ( kkrec==6 ) then
              if ( k4==1 ) then
                write (inname,99009) nyear , 'omega.' , nyear
              else if ( k4==2 ) then
                write (inname,99010) nyear , 'omega.' , nyear
              else if ( k4==3 ) then
                write (inname,99011) nyear , 'omega.' , nyear
              else if ( k4==4 ) then
                write (inname,99012) nyear , 'omega.' , nyear
              else
              end if
            else
            end if
 
            pathaddname = '../DATA/'//dattyp//'/'//inname
            inquire (file=pathaddname,exist=there)
            if ( .not.there ) then
              print * , pathaddname , ' is not available'
              stop
            end if
            istatus = nf90_open(pathaddname,nf90_nowrite,               &
                   & inet6(kkrec,k4))
            istatus = nf90_get_att(inet6(kkrec,k4),5,'scale_factor',    &
                   & xscl(kkrec,k4))
            istatus = nf90_get_att(inet6(kkrec,k4),5,'add_offset',      &
                   & xoff(kkrec,k4))
            write (*,*) inet6(kkrec,k4) , pathaddname , xscl(kkrec,k4) ,&
                      & xoff(kkrec,k4)
          end do
        end do
 
        if ( lsmtyp=='USGS' ) then
          do k4 = 1 , 4
            do l4 = 1 , 4
              do kkrec = 1 , 3
                if ( kkrec==1 ) then
                  if ( k4==1 ) then
                    write (inname,99013) 'Qsoil_' , l4
                  else if ( k4==2 ) then
                    write (inname,99014) 'Qsoil_' , l4
                  else if ( k4==3 ) then
                    write (inname,99015) 'Qsoil_' , l4
                  else if ( k4==4 ) then
                    write (inname,99016) 'Qsoil_' , l4
                  else
                  end if
                else if ( kkrec==2 ) then
                  if ( k4==1 ) then
                    write (inname,99017) 'Tice_' , l4
                  else if ( k4==2 ) then
                    write (inname,99018) 'Tice_' , l4
                  else if ( k4==3 ) then
                    write (inname,99019) 'Tice_' , l4
                  else if ( k4==4 ) then
                    write (inname,99020) 'Tice_' , l4
                  else
                  end if
                else if ( kkrec==3 ) then
                  if ( k4==1 ) then
                    write (inname,99013) 'Tsoil_' , l4
                  else if ( k4==2 ) then
                    write (inname,99014) 'Tsoil_' , l4
                  else if ( k4==3 ) then
                    write (inname,99015) 'Tsoil_' , l4
                  else if ( k4==4 ) then
                    write (inname,99016) 'Tsoil_' , l4
                  else
                  end if
                else
                end if
 
                pathaddname = '../DATA/'//dattyp//'/0surface/'//inname
                inquire (file=pathaddname,exist=there)
                if ( .not.there ) then
                  print * , pathaddname , ' is not available'
                  stop
                end if
                istatus = nf90_open(pathaddname,nf90_nowrite,           &
                       & isnet3(kkrec,l4,k4))
                istatus = nf90_get_att(isnet3(kkrec,l4,k4),4,           &
                        &'scale_factor',xscl_s(kkrec,l4,k4))
                istatus = nf90_get_att(isnet3(kkrec,l4,k4),4,           &
                        &'add_offset',xoff_s(kkrec,l4,k4))
                write (*,*) isnet3(kkrec,l4,k4) , pathaddname ,         &
                          & xscl_s(kkrec,l4,k4) , xoff_s(kkrec,l4,k4)
              end do
            end do
 
            if ( k4==1 ) then
              write (inname,99021) 'snowdpth'
            else if ( k4==2 ) then
              write (inname,99022) 'snowdpth'
            else if ( k4==3 ) then
              write (inname,99023) 'snowdpth'
            else if ( k4==4 ) then
              write (inname,99024) 'snowdpth'
            else
            end if
            pathaddname = '../DATA/'//dattyp//'/0surface/'//inname
            inquire (file=pathaddname,exist=there)
            if ( .not.there ) then
              print * , pathaddname , ' is not available'
              stop
            end if
            istatus = nf90_open(pathaddname,nf90_nowrite,isnow(k4))
            istatus = nf90_get_att(isnow(k4),4,'scale_factor',          &
                   & xscl_sn(k4))
            istatus = nf90_get_att(isnow(k4),4,'add_offset',            &
                   & xoff_sn(k4))
            write (*,*) isnow(k4) , pathaddname , xscl_sn(k4) ,         &
                      & xoff_sn(k4)
 
          end do
        end if
      end if
 
      k4 = nhour/6 + 1
      it = nday
      if ( month==2 ) it = it + 31
      if ( month==3 ) it = it + 59
      if ( month==4 ) it = it + 90
      if ( month==5 ) it = it + 120
      if ( month==6 ) it = it + 151
      if ( month==7 ) it = it + 181
      if ( month==8 ) it = it + 212
      if ( month==9 ) it = it + 243
      if ( month==10 ) it = it + 273
      if ( month==11 ) it = it + 304
      if ( month==12 ) it = it + 334
      if ( mod(nyear,4)==0 .and. month>2 ) it = it + 1
      if ( mod(nyear,100)==0 .and. month>2 ) it = it - 1
      if ( mod(nyear,400)==0 .and. month>2 ) it = it + 1
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
      if ( mod(nyear,4)==0 ) icount(4) = 366
      if ( mod(nyear,100)==0 ) icount(4) = 365
      if ( mod(nyear,400)==0 ) icount(4) = 366
      if ( nyear==2002 ) icount(4) = 243
      if ( nyear==1957 ) icount(4) = 122
      if ( nyear==1957 ) it = it - 243
      istart(4) = it
      icount(4) = 1
!bxq_
      do kkrec = 1 , 5
        inet = inet6(kkrec,k4)
        istatus = nf90_get_var(inet,5,work,istart,icount)
        xscale = xscl(kkrec,k4)
        xadd = xoff(kkrec,k4)
        if ( kkrec==1 ) then
          do k = 1 , klev
            do j = 1 , jlat
              do i = 1 , ilon
                tvar(i,jlat+1-j,k) = work(i,j,k)*xscale + xadd
              end do
            end do
          end do
        else if ( kkrec==2 ) then
          do k = 1 , klev
            do j = 1 , jlat
              do i = 1 , ilon
                hvar(i,jlat+1-j,k) = work(i,j,k)*xscale + xadd
                hvar(i,jlat+1-j,k) = hvar(i,jlat+1-j,k)/9.80616
              end do
            end do
          end do
        else if ( kkrec==3 ) then
          do k = 1 , klev
            do j = 1 , jlat
              do i = 1 , ilon
                rhvar(i,jlat+1-j,k) = work(i,j,k)*xscale + xadd
                rhvar(i,jlat+1-j,k) = rhvar(i,jlat+1-j,k)*0.01
!               RHvar(i,jlat+1-j,k)=amax1(RHvar(i,jlat+1-j,k),1.05)
              end do
            end do
          end do
        else if ( kkrec==4 ) then
          do k = 1 , klev
            do j = 1 , jlat
              do i = 1 , ilon
                uvar(i,jlat+1-j,k) = work(i,j,k)*xscale + xadd
              end do
            end do
          end do
        else if ( kkrec==5 ) then
          do k = 1 , klev
            do j = 1 , jlat
              do i = 1 , ilon
                vvar(i,jlat+1-j,k) = work(i,j,k)*xscale + xadd
              end do
            end do
          end do
        else if ( kkrec==6 ) then
          do k = 1 , klev
            do j = 1 , jlat
              do i = 1 , ilon
                wvar(i,jlat+1-j,k) = work(i,j,k)*xscale + xadd
              end do
            end do
          end do
        else
        end if
        istatus = nf90_close(inet)
      end do
 
      if ( lsmtyp=='USGS' ) then
        k4 = nhour/6 + 1
        it = nday
        if ( month==2 ) it = it + 31
        if ( month==3 ) it = it + 59
        if ( month==4 ) it = it + 90
        if ( month==5 ) it = it + 120
        if ( month==6 ) it = it + 151
        if ( month==7 ) it = it + 181
        if ( month==8 ) it = it + 212
        if ( month==9 ) it = it + 243
        if ( month==10 ) it = it + 273
        if ( month==11 ) it = it + 304
        if ( month==12 ) it = it + 334
        if ( mod(nyear,4)==0 .and. month>2 ) it = it + 1
        if ( mod(nyear,100)==0 .and. month>2 ) it = it - 1
        if ( mod(nyear,400)==0 .and. month>2 ) it = it + 1
        do k = 1957 , nyear - 1
          it = it + 365
          if ( mod(k,4)==0 ) it = it + 1
        end do
        it = it - 243
 
        do k = 1 , 3
          istart(k) = 1
        end do
        do k = 4 , 10
          istart(k) = 0
          icount(k) = 0
        end do
        icount(1) = ilon
        icount(2) = jlat
        istart(3) = it
        icount(3) = 1
!bxq_
        do l4 = 1 , 4
          do kkrec = 1 , 3
            inet = isnet3(kkrec,l4,k4)
            istatus = nf90_get_var(inet,4,work2d,istart,icount)
            xscale = xscl_s(kkrec,l4,k4)
            xadd = xoff_s(kkrec,l4,k4)
            if ( kkrec==1 ) then
              do j = 1 , jlat
                do i = 1 , ilon
                  qsoil(i,jlat+1-j,l4) = work2d(i,j)*xscale + xadd
                end do
              end do
            else if ( kkrec==2 ) then
              do j = 1 , jlat
                do i = 1 , ilon
                  tsice(i,jlat+1-j,l4) = work2d(i,j)*xscale + xadd
                end do
              end do
            else if ( kkrec==3 ) then
              do j = 1 , jlat
                do i = 1 , ilon
                  tsoil(i,jlat+1-j,l4) = work2d(i,j)*xscale + xadd
                end do
              end do
            else
            end if
            istatus = nf90_close(inet)
          end do
        end do
        inet = isnow(k4)
        istatus = nf90_get_var(inet,4,work2d,istart,icount)
        xscale = xscl_sn(k4)
        xadd = xoff_sn(k4)
        do j = 1 , jlat
          do i = 1 , ilon
            snw(i,jlat+1-j) = work2d(i,j)*xscale + xadd
          end do
        end do
        istatus = nf90_close(inet)
      end if
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
99013 format (a6,i1,'L.00.nc')
99014 format (a6,i1,'L.06.nc')
99015 format (a6,i1,'L.12.nc')
99016 format (a6,i1,'L.18.nc')
99017 format (a5,i1,'L.00.nc')
99018 format (a5,i1,'L.06.nc')
99019 format (a5,i1,'L.12.nc')
99020 format (a5,i1,'L.18.nc')
99021 format (a8,'.00.nc')
99022 format (a8,'.06.nc')
99023 format (a8,'.12.nc')
99024 format (a8,'.18.nc')
!
      end subroutine era6hour

      subroutine headerera
      implicit none
!
! Local variables
!
      integer :: i , j , k , kr
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
 
      end subroutine headerera

      end module mod_era40

      module mod_ncep

      use mod_domain

      implicit none

      integer , parameter :: klev = 13 , jlat = 73 , ilon = 144
      real(4) , dimension(ilon,jlat) :: psvar
      real , dimension(jlat) :: glat
      real , dimension(ilon) :: glon
      real , dimension(klev) :: sigma1 , sigmar
      real :: psref
      real :: delx , grdfac
      character(6) :: lgtype

      real , dimension(ilon,jlat,klev) :: wvar

      real , target , dimension(ilon,jlat,klev*3) :: b2
      real , target , dimension(ilon,jlat,klev*2) :: d2
      real , target , dimension(jx,iy,klev*3) :: b3
      real , target , dimension(jx,iy,klev*2) :: d3
      
      real , pointer :: u3(:,:,:) , v3(:,:,:)
      real , pointer :: h3(:,:,:) , q3(:,:,:) , t3(:,:,:)
      real , pointer :: uvar(:,:,:) , vvar(:,:,:)
      real , pointer :: hvar(:,:,:) , rhvar(:,:,:) , tvar(:,:,:)

      contains

      subroutine getncep(idate)
      use mod_domain
      use mod_geo
      use mod_var4
      implicit none
!
! Dummy arguments
!
      integer :: idate
!
! Local variables
!
      real , dimension(jx,iy) :: b3pd , pa , sst1 , sst2 , tlayer , za
      integer :: nmop , nyrp
      real :: wt

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
      call cdc6hour(dattyp,idate,idate1)

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
 
      if ( ssttyp/='OI_WK' ) then
!       F1  CALCULATE SSTS FOR DATE FROM OBSERVED SSTS
        print * , 'INPUT DAY FOR SST DATA ACQUISITION:' , idate
        call julian(idate,nyrp,nmop,wt)
!
        call mksst(ts4,sst1,sst2,topogm,xlandu,jx,iy,nyrp,nmop,wt)
      else
        call mksst2(ts4,sst1,sst2,topogm,xlandu,jx,iy,idate/100)
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
      call writef(u4,v4,t4,q4,ps4,ts4,ptop,jx,iy,kz,idate)
!
      end subroutine getncep

      subroutine cdc6hour(dattyp,idate,idate0)
      use netcdf
      implicit none
!
! Dummy arguments
!
      character(5) :: dattyp
      integer :: idate , idate0
      intent (in) dattyp , idate , idate0
!
! Local variables
!
      integer :: i , ilev , inet , it , j , k , kkrec , m , month ,     &
               & nday , nhour , nlev , nyear , istatus
      character(21) :: inname
      character(35) :: pathaddname
      logical :: there
      character(5) , dimension(7) :: varname
      integer(2) , dimension(ilon,jlat,klev) :: work
      real(8) :: xadd , xscale
      integer , dimension(10) :: icount , istart
      integer , dimension(6) :: inet7 , ivar7
      real(8) , dimension(7) :: xoff , xscl
!
!     This is the latitude, longitude dimension of the grid to be read.
!     This corresponds to the lat and lon dimension variables in the
!     netCDF file.
!
!     The data are packed into short integers (INTEGER*2).  The array
!     work will be used to hold the packed integers.
!
!     DATA ARRAY AND WORK ARRAY
!bxq
!bxq_
      data varname/'air' , 'hgt' , 'rhum' , 'uwnd' , 'vwnd' , 'omega' , &
          &'pres'/
!
!     Below in the ncopen call is the file name of the netCDF file.
!     You may want to add code to read in the file name and the
!     variable name.
!     OPEN FILE AND GET FILES ID AND VARIABLE ID(S)
!
!bxq
      nyear = idate/1000000
      month = idate/10000 - nyear*100
      nday = idate/100 - nyear*10000 - month*100
      nhour = idate - nyear*1000000 - month*10000 - nday*100
!fix  do kkrec=1,7
      do kkrec = 1 , 5
        if ( dattyp=='NNRP1' ) then
          if ( kkrec==1 .or. kkrec==2 .or. kkrec==4 .or. kkrec==5 )     &
             & nlev = klev
          if ( kkrec==6 ) nlev = 12
          if ( kkrec==3 ) nlev = 8
          if ( kkrec==7 ) nlev = 0
        else if ( dattyp=='NNRP2' ) then
          if ( kkrec<=6 ) nlev = klev
          if ( kkrec==7 ) nlev = 0
        else
        end if
        if ( idate==idate0 .or.                                         &
           & (mod(idate,100000)==10100 .and. mod(idate,1000000)/=110100)&
           & ) then
          if ( kkrec==1 ) then
            write (inname,99001) nyear , 'air.' , nyear
          else if ( kkrec==2 ) then
            write (inname,99001) nyear , 'hgt.' , nyear
          else if ( kkrec==3 ) then
            write (inname,99002) nyear , 'rhum.' , nyear
          else if ( kkrec==4 ) then
            write (inname,99002) nyear , 'uwnd.' , nyear
          else if ( kkrec==5 ) then
            write (inname,99002) nyear , 'vwnd.' , nyear
          else if ( kkrec==6 ) then
            write (inname,99003) nyear , 'omega.' , nyear
          else if ( kkrec==7 ) then
            write (inname,99004) nyear , 'pres.sfc.' , nyear
          else
          end if
 
          if ( dattyp=='NNRP1' ) then
            pathaddname = '../DATA/NNRP1/'//inname
          else if ( dattyp=='NNRP2' ) then
            pathaddname = '../DATA/NNRP2/'//inname
          else
          end if
          inquire (file=pathaddname,exist=there)
          if ( .not.there ) then
            print * , pathaddname , ' is not available'
            stop
          end if
          istatus = nf90_open(pathaddname,nf90_nowrite,inet7(kkrec))
          istatus = nf90_get_att(inet7(kkrec),5,'scale_factor',         &
                 & xscl(kkrec))
          istatus = nf90_get_att(inet7(kkrec),5,'add_offset',           &
                 & xoff(kkrec))
          write (*,*) inet7(kkrec) , pathaddname , xscl(kkrec) ,        &
                    & xoff(kkrec)
        end if
 
        it = (nday-1)*4 + nhour/6 + 1
        if ( month==2 ) it = it + 31*4
        if ( month==3 ) it = it + 59*4
        if ( month==4 ) it = it + 90*4
        if ( month==5 ) it = it + 120*4
        if ( month==6 ) it = it + 151*4
        if ( month==7 ) it = it + 181*4
        if ( month==8 ) it = it + 212*4
        if ( month==9 ) it = it + 243*4
        if ( month==10 ) it = it + 273*4
        if ( month==11 ) it = it + 304*4
        if ( month==12 ) it = it + 334*4
        if ( mod(nyear,4)==0 .and. month>2 ) it = it + 4
        if ( mod(nyear,100)==0 .and. month>2 ) it = it - 4
        if ( mod(nyear,400)==0 .and. month>2 ) it = it + 4
!bxq_
        do m = 1 , 4
          istart(m) = 1
        end do
        do m = 5 , 10
          istart(m) = 0
          icount(m) = 0
        end do
        icount(1) = ilon
        icount(2) = jlat
        icount(4) = 1460
        if ( mod(nyear,4)==0 ) icount(4) = 1464
        if ( mod(nyear,100)==0 ) icount(4) = 1460
        if ( mod(nyear,400)==0 ) icount(4) = 1464
        istart(4) = it
        icount(4) = 1
        inet = inet7(kkrec)
        if ( nlev>0 ) then
          icount(3) = nlev
          istatus = nf90_get_var(inet,5,work,istart,icount)
          xscale = xscl(kkrec)
          xadd = xoff(kkrec)
          do ilev = 1 , nlev
            if ( kkrec==1 ) then
              do j = 1 , jlat
                do i = 1 , ilon
                  tvar(i,jlat+1-j,14-ilev) = work(i,j,ilev)             &
                  & *xscale + xadd
                end do
              end do
            else if ( kkrec==2 ) then
              do j = 1 , jlat
                do i = 1 , ilon
                  hvar(i,jlat+1-j,14-ilev) = work(i,j,ilev)             &
                  & *xscale + xadd
                end do
              end do
            else if ( kkrec==3 ) then
              do j = 1 , jlat
                do i = 1 , ilon
                  rhvar(i,jlat+1-j,14-ilev)                             &
                  & = dmin1((work(i,j,ilev)*xscale+xadd)*0.01,1.D0)
                end do
              end do
            else if ( kkrec==4 ) then
              do j = 1 , jlat
                do i = 1 , ilon
                  uvar(i,jlat+1-j,14-ilev) = work(i,j,ilev)             &
                  & *xscale + xadd
                end do
              end do
            else if ( kkrec==5 ) then
              do j = 1 , jlat
                do i = 1 , ilon
                  vvar(i,jlat+1-j,14-ilev) = work(i,j,ilev)             &
                  & *xscale + xadd
                end do
              end do
            else if ( kkrec==6 ) then
              do j = 1 , jlat
                do i = 1 , ilon
                  wvar(i,jlat+1-j,14-ilev) = work(i,j,ilev)             &
                  & *xscale + xadd
                end do
              end do
            else
            end if
          end do
        else if ( nlev==0 ) then
          icount(3) = 1
          istatus = nf90_get_var(inet,5,work,istart,icount)
          if ( kkrec==7 ) then
            do j = 1 , jlat
              do i = 1 , ilon
                psvar(i,jlat+1-j) = work(i,j,1)*xscale + xadd
              end do
            end do
          end if
        else
        end if
        if ( dattyp=='NNRP1' ) then
!         It's a pity that we have to nudge the values by the following
!         way
          do k = 5 , 1 , -1
            do j = 1 , jlat
              do i = 1 , ilon
                rhvar(i,j,k) = rhvar(i,j,k+1)
              end do
            end do
          end do
 
          do j = 1 , jlat
            do i = 1 , ilon
              wvar(i,j,1) = 0.0
            end do
          end do
        end if
      end do
99001 format (i4,'/',a4,i4,'.nc')
99002 format (i4,'/',a5,i4,'.nc')
99003 format (i4,'/',a6,i4,'.nc')
99004 format (i4,'/',a9,i4,'.nc')
!
      end subroutine cdc6hour

      end module mod_ncep

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

      module mod_ein15

      use mod_regcm_param , only : iy , jx , kz , ibyte , dattyp
      use mod_preproc_param

      implicit none

      private

      integer , parameter :: klev = 23 , jlat = 121 , ilon = 240

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

      integer , dimension(5,4) :: inet5
      integer , dimension(5,4) :: ivar5
      real(8) , dimension(5,4) :: xoff , xscl

      public :: getein15 , headerein15

      contains

      subroutine getein15(idate)
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
      integer :: nmop , nyrp
      real :: wt
!
!
!     READ DATA AT IDATE
!
      call ein156hour(dattyp,idate,idate1)

      write (*,*) 'READ IN fields at DATE:' , idate
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
!
      call intgtb(pa,za,tlayer,topogm,t3,h3,sigmar,jx,iy,klev)
 
      call intpsn(ps4,topogm,pa,za,tlayer,ptop,jx,iy)
      call p1p2(b3pd,ps4,jx,iy)
 
!     INTERPOLATION FROM PRESSURE LEVELS
      call intv3(ts4,t3,ps4,sigmar,ptop,jx,iy,klev)

      if ( ssttyp=='OI_WK' ) then
        call mksst2(ts4,sst1,sst2,topogm,xlandu,jx,iy,idate/100)
      else if ( ssttyp=='OI2WK' ) then
        call mksst2a(ts4,sst1,sst2,ice1,ice2,topogm,xlandu,jx,iy,       &
              &      idate/100)
      else if ( ssttyp=='ERSST' .or. ssttyp=='ERSKT' ) then
        call mksst3(ts4,sst1,topogm,xlandu,jx,iy,idate)
      else
!       F1  CALCULATE SSTS FOR DATE FROM OBSERVED SSTS
        print * , 'INPUT DAY FOR SST DATA ACQUISITION:' , idate
        call julian(idate,nyrp,nmop,wt)
!
        if ( ssttyp=='OI2ST' ) then
          call mkssta(ts4,sst1,sst2,ice1,ice2,topogm,xlandu,jx,iy,nyrp, &
                   &  nmop,wt)
        else
          call mksst(ts4,sst1,sst2,topogm,xlandu,jx,iy,nyrp,nmop,wt)
        end if 
      end if
 
!     F3  INTERPOLATE U, V, T, AND Q.

      call intv1(u4,u3,b3pd,sigma2,sigmar,ptop,jx,iy,kz,klev)
      call intv1(v4,v3,b3pd,sigma2,sigmar,ptop,jx,iy,kz,klev)
      call intv2(t4,t3,ps4,sigma2,sigmar,ptop,jx,iy,kz,klev)
      call intv1(q4,q3,ps4,sigma2,sigmar,ptop,jx,iy,kz,klev)
! 
      call humid2(t4,q4,ps4,ptop,sigma2,jx,iy,kz)
!
!     F4  DETERMINE H
      call hydrost(h4,t4,topogm,ps4,ptop,sigmaf,sigma2,dsigma,jx,iy,kz)
!
!     G   WRITE AN INITIAL FILE FOR THE RegCM
      call writef(ptop,idate)
!
      end subroutine getein15
!
!-----------------------------------------------------------------------
!
      subroutine ein156hour(dattyp,idate,idate0)
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
      integer :: i , inet , it , j , k , k4 , kkrec , month , nday ,    &
               & nhour , nyear , istatus , ivar
      character(24) :: inname
      character(38) :: pathaddname
      logical :: there
      character(1) , dimension(5) :: varname
      real(8) :: xadd , xscale
      integer , dimension(10) :: icount , istart
      integer(2) , dimension(ilon,jlat,37) :: work
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
      data varname/'t' , 'z' , 'r' , 'u' , 'v'/
!
!     Below in the ncopen call is the file name of the netCDF file.
!     You may want to add code to read in the file name and the
!     variable name.
!     OPEN FILE AND GET FILES ID AND VARIABLE ID(S)
!
!bxq
      if ( idate<1989010100 .or. idate>2009053118 ) then
        write (*,*) 'EIN15 datasets is just available from' ,           &
                   &' 1989010100 to 2009053118'
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
            else
            end if
 
            pathaddname = '../DATA/'//dattyp//'/'//inname
            inquire (file=pathaddname,exist=there)
            if ( .not.there ) then
              print * , pathaddname , ' is not available'
              stop
            end if
            istatus = nf90_open(pathaddname,nf90_nowrite,               &
                   & inet5(kkrec,k4))
            if ( istatus /= nf90_noerr) then
              write (*,*) 'Error opening ' , trim(pathaddname)
              stop
            end if
            istatus = nf90_inq_varid(inet5(kkrec,k4),varname(kkrec),    &
                   & ivar5(kkrec,k4))
            if ( istatus /= nf90_noerr) then
              write (*,*) 'Error searching var ' , varname(kkrec)
              stop
            end if
            istatus = nf90_get_att(inet5(kkrec,k4),ivar5(kkrec,k4),     &
                   & 'scale_factor',xscl(kkrec,k4))
            if ( istatus /= nf90_noerr) then
              write (*,*) 'Error attribure scale_factor for var ' ,     &
                     &     varname(kkrec)
              stop
            end if
            istatus = nf90_get_att(inet5(kkrec,k4),ivar5(kkrec,k4),     &
                   & 'add_offset',xoff(kkrec,k4))
            if ( istatus /= nf90_noerr) then
              write (*,*) 'Error attribure add_offset for var ' ,     &
                     &     varname(kkrec)
              stop
            end if
            write (*,*) inet5(kkrec,k4) , pathaddname , xscl(kkrec,k4) ,&
                      & xoff(kkrec,k4)
          end do
        end do
 
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
      icount(3) = 37
      icount(4) = 365
      if ( mod(nyear,4)==0 ) icount(4) = 366
      if ( mod(nyear,100)==0 ) icount(4) = 365
      if ( mod(nyear,400)==0 ) icount(4) = 366
      istart(4) = it
      icount(4) = 1
!bxq_
      do kkrec = 1 , 5
        inet = inet5(kkrec,k4)
        ivar = ivar5(kkrec,k4)
        istatus = nf90_get_var(inet,ivar,work,istart,icount)
        if (istatus /= nf90_noerr) then
          write (*,*) 'Error reading variable ' , varname(kkrec) ,      &
              & ' at ' , istart , ' for ' , icount
          stop
        end if
        xscale = xscl(kkrec,k4)
        xadd = xoff(kkrec,k4)
        if ( kkrec==1 ) then
          do j = 1 , jlat
            do i = 1 , ilon
!             Tvar(i,jlat+1-j,k)=work(i,j,k)*xscale+xadd
              tvar(i,jlat+1-j,1) = work(i,j,1)*xscale + xadd
              tvar(i,jlat+1-j,2) = work(i,j,2)*xscale + xadd
              tvar(i,jlat+1-j,3) = work(i,j,3)*xscale + xadd
              tvar(i,jlat+1-j,4) = work(i,j,4)*xscale + xadd
              tvar(i,jlat+1-j,5) = work(i,j,5)*xscale + xadd
              tvar(i,jlat+1-j,6) = work(i,j,6)*xscale + xadd
              tvar(i,jlat+1-j,7) = work(i,j,7)*xscale + xadd
              tvar(i,jlat+1-j,8) = work(i,j,8)*xscale + xadd
              tvar(i,jlat+1-j,9) = work(i,j,9)*xscale + xadd
              tvar(i,jlat+1-j,10) = work(i,j,10)*xscale + xadd
              tvar(i,jlat+1-j,11) = work(i,j,11)*xscale + xadd
              tvar(i,jlat+1-j,12) = work(i,j,13)*xscale + xadd
              tvar(i,jlat+1-j,13) = work(i,j,15)*xscale + xadd
              tvar(i,jlat+1-j,14) = work(i,j,17)*xscale + xadd
              tvar(i,jlat+1-j,15) = work(i,j,18)*xscale + xadd
              tvar(i,jlat+1-j,16) = work(i,j,20)*xscale + xadd
              tvar(i,jlat+1-j,17) = work(i,j,22)*xscale + xadd
              tvar(i,jlat+1-j,18) = work(i,j,24)*xscale + xadd
              tvar(i,jlat+1-j,19) = work(i,j,26)*xscale + xadd
              tvar(i,jlat+1-j,20) = work(i,j,28)*xscale + xadd
              tvar(i,jlat+1-j,21) = work(i,j,31)*xscale + xadd
              tvar(i,jlat+1-j,22) = work(i,j,34)*xscale + xadd
              tvar(i,jlat+1-j,23) = work(i,j,37)*xscale + xadd
            end do
          end do
        else if ( kkrec==2 ) then
          do j = 1 , jlat
            do i = 1 , ilon
!             Hvar(i,jlat+1-j,k)=work(i,j,k)*xscale+xadd
              hvar(i,jlat+1-j,1) = work(i,j,1)*xscale + xadd
              hvar(i,jlat+1-j,2) = work(i,j,2)*xscale + xadd
              hvar(i,jlat+1-j,3) = work(i,j,3)*xscale + xadd
              hvar(i,jlat+1-j,4) = work(i,j,4)*xscale + xadd
              hvar(i,jlat+1-j,5) = work(i,j,5)*xscale + xadd
              hvar(i,jlat+1-j,6) = work(i,j,6)*xscale + xadd
              hvar(i,jlat+1-j,7) = work(i,j,7)*xscale + xadd
              hvar(i,jlat+1-j,8) = work(i,j,8)*xscale + xadd
              hvar(i,jlat+1-j,9) = work(i,j,9)*xscale + xadd
              hvar(i,jlat+1-j,10) = work(i,j,10)*xscale + xadd
              hvar(i,jlat+1-j,11) = work(i,j,11)*xscale + xadd
              hvar(i,jlat+1-j,12) = work(i,j,13)*xscale + xadd
              hvar(i,jlat+1-j,13) = work(i,j,15)*xscale + xadd
              hvar(i,jlat+1-j,14) = work(i,j,17)*xscale + xadd
              hvar(i,jlat+1-j,15) = work(i,j,18)*xscale + xadd
              hvar(i,jlat+1-j,16) = work(i,j,20)*xscale + xadd
              hvar(i,jlat+1-j,17) = work(i,j,22)*xscale + xadd
              hvar(i,jlat+1-j,18) = work(i,j,24)*xscale + xadd
              hvar(i,jlat+1-j,19) = work(i,j,26)*xscale + xadd
              hvar(i,jlat+1-j,20) = work(i,j,28)*xscale + xadd
              hvar(i,jlat+1-j,21) = work(i,j,31)*xscale + xadd
              hvar(i,jlat+1-j,22) = work(i,j,34)*xscale + xadd
              hvar(i,jlat+1-j,23) = work(i,j,37)*xscale + xadd
            end do
          end do
          do k = 1 , klev
            do j = 1 , jlat
              do i = 1 , ilon
                hvar(i,j,k) = hvar(i,j,k)/9.80616
              end do
            end do
          end do
        else if ( kkrec==3 ) then
          do j = 1 , jlat
            do i = 1 , ilon
!             RHvar(i,jlat+1-j,k)=work(i,j,k)*xscale+xadd
              rhvar(i,jlat+1-j,1) = work(i,j,1)*xscale + xadd
              rhvar(i,jlat+1-j,2) = work(i,j,2)*xscale + xadd
              rhvar(i,jlat+1-j,3) = work(i,j,3)*xscale + xadd
              rhvar(i,jlat+1-j,4) = work(i,j,4)*xscale + xadd
              rhvar(i,jlat+1-j,5) = work(i,j,5)*xscale + xadd
              rhvar(i,jlat+1-j,6) = work(i,j,6)*xscale + xadd
              rhvar(i,jlat+1-j,7) = work(i,j,7)*xscale + xadd
              rhvar(i,jlat+1-j,8) = work(i,j,8)*xscale + xadd
              rhvar(i,jlat+1-j,9) = work(i,j,9)*xscale + xadd
              rhvar(i,jlat+1-j,10) = work(i,j,10)*xscale + xadd
              rhvar(i,jlat+1-j,11) = work(i,j,11)*xscale + xadd
              rhvar(i,jlat+1-j,12) = work(i,j,13)*xscale + xadd
              rhvar(i,jlat+1-j,13) = work(i,j,15)*xscale + xadd
              rhvar(i,jlat+1-j,14) = work(i,j,17)*xscale + xadd
              rhvar(i,jlat+1-j,15) = work(i,j,18)*xscale + xadd
              rhvar(i,jlat+1-j,16) = work(i,j,20)*xscale + xadd
              rhvar(i,jlat+1-j,17) = work(i,j,22)*xscale + xadd
              rhvar(i,jlat+1-j,18) = work(i,j,24)*xscale + xadd
              rhvar(i,jlat+1-j,19) = work(i,j,26)*xscale + xadd
              rhvar(i,jlat+1-j,20) = work(i,j,28)*xscale + xadd
              rhvar(i,jlat+1-j,21) = work(i,j,31)*xscale + xadd
              rhvar(i,jlat+1-j,22) = work(i,j,34)*xscale + xadd
              rhvar(i,jlat+1-j,23) = work(i,j,37)*xscale + xadd
            end do
          end do
          do k = 1 , 23
            do j = 1 , jlat
              do i = 1 , ilon
                rhvar(i,j,k) = amax1(rhvar(i,j,k)*0.01,0.00)
              end do
            end do
          end do
        else if ( kkrec==4 ) then
          do j = 1 , jlat
            do i = 1 , ilon
!             Uvar(i,jlat+1-j,k)=work(i,j,k)*xscale+xadd
              uvar(i,jlat+1-j,1) = work(i,j,1)*xscale + xadd
              uvar(i,jlat+1-j,2) = work(i,j,2)*xscale + xadd
              uvar(i,jlat+1-j,3) = work(i,j,3)*xscale + xadd
              uvar(i,jlat+1-j,4) = work(i,j,4)*xscale + xadd
              uvar(i,jlat+1-j,5) = work(i,j,5)*xscale + xadd
              uvar(i,jlat+1-j,6) = work(i,j,6)*xscale + xadd
              uvar(i,jlat+1-j,7) = work(i,j,7)*xscale + xadd
              uvar(i,jlat+1-j,8) = work(i,j,8)*xscale + xadd
              uvar(i,jlat+1-j,9) = work(i,j,9)*xscale + xadd
              uvar(i,jlat+1-j,10) = work(i,j,10)*xscale + xadd
              uvar(i,jlat+1-j,11) = work(i,j,11)*xscale + xadd
              uvar(i,jlat+1-j,12) = work(i,j,13)*xscale + xadd
              uvar(i,jlat+1-j,13) = work(i,j,15)*xscale + xadd
              uvar(i,jlat+1-j,14) = work(i,j,17)*xscale + xadd
              uvar(i,jlat+1-j,15) = work(i,j,18)*xscale + xadd
              uvar(i,jlat+1-j,16) = work(i,j,20)*xscale + xadd
              uvar(i,jlat+1-j,17) = work(i,j,22)*xscale + xadd
              uvar(i,jlat+1-j,18) = work(i,j,24)*xscale + xadd
              uvar(i,jlat+1-j,19) = work(i,j,26)*xscale + xadd
              uvar(i,jlat+1-j,20) = work(i,j,28)*xscale + xadd
              uvar(i,jlat+1-j,21) = work(i,j,31)*xscale + xadd
              uvar(i,jlat+1-j,22) = work(i,j,34)*xscale + xadd
              uvar(i,jlat+1-j,23) = work(i,j,37)*xscale + xadd
            end do
          end do
        else if ( kkrec==5 ) then
          do j = 1 , jlat
            do i = 1 , ilon
!             Vvar(i,jlat+1-j,k) =work(i,j,k)*xscale+xadd
              vvar(i,jlat+1-j,1) = work(i,j,1)*xscale + xadd
              vvar(i,jlat+1-j,2) = work(i,j,2)*xscale + xadd
              vvar(i,jlat+1-j,3) = work(i,j,3)*xscale + xadd
              vvar(i,jlat+1-j,4) = work(i,j,4)*xscale + xadd
              vvar(i,jlat+1-j,5) = work(i,j,5)*xscale + xadd
              vvar(i,jlat+1-j,6) = work(i,j,6)*xscale + xadd
              vvar(i,jlat+1-j,7) = work(i,j,7)*xscale + xadd
              vvar(i,jlat+1-j,8) = work(i,j,8)*xscale + xadd
              vvar(i,jlat+1-j,9) = work(i,j,9)*xscale + xadd
              vvar(i,jlat+1-j,10) = work(i,j,10)*xscale + xadd
              vvar(i,jlat+1-j,11) = work(i,j,11)*xscale + xadd
              vvar(i,jlat+1-j,12) = work(i,j,13)*xscale + xadd
              vvar(i,jlat+1-j,13) = work(i,j,15)*xscale + xadd
              vvar(i,jlat+1-j,14) = work(i,j,17)*xscale + xadd
              vvar(i,jlat+1-j,15) = work(i,j,18)*xscale + xadd
              vvar(i,jlat+1-j,16) = work(i,j,20)*xscale + xadd
              vvar(i,jlat+1-j,17) = work(i,j,22)*xscale + xadd
              vvar(i,jlat+1-j,18) = work(i,j,24)*xscale + xadd
              vvar(i,jlat+1-j,19) = work(i,j,26)*xscale + xadd
              vvar(i,jlat+1-j,20) = work(i,j,28)*xscale + xadd
              vvar(i,jlat+1-j,21) = work(i,j,31)*xscale + xadd
              vvar(i,jlat+1-j,22) = work(i,j,34)*xscale + xadd
              vvar(i,jlat+1-j,23) = work(i,j,37)*xscale + xadd
            end do
          end do
        else
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
!
      end subroutine ein156hour
!
!-----------------------------------------------------------------------
!
      subroutine headerein15
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
        glon(i) = float(i-1)*1.50
      end do
      do j = 1 , jlat
        glat(j) = -90.0 + float(j-1)*1.50
      end do
!HH:OVER
!     CHANGE ORDER OF VERTICAL INDEXES FOR PRESSURE LEVELS
!
      do k = 1 , klev
        kr = klev - k + 1
        sigma1(k) = sigmar(kr)
      end do
 
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

      end subroutine headerein15

      end module mod_ein15

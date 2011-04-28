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

      module mod_sst_ersst

      contains

      subroutine sst_ersst

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! Comments on dataset sources and location:                          !
!                                                                    !
! ERAIN    ERAIN_SST is provided by ERA-Interim project              !
!          6 hourly frequncy, 1.5x1.5 degree resolution              !
!          from 1989010100 to 2009053118.                            !
!          'ERSST' for using the sea surface temperature;            !
!          'ERSKT' for using the skin temperature.                   !
!                                                                    !
!          ML= 1 is   0.0; ML= 2 is   1.5; => ML=240 is 358.5E       !
!          NL= 1 is  90.0; ML= 2 is  88.5; => ML=121 is -90.         !
!                                                                    !
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      use mod_sst_grid
      use mod_date
      use mod_interp , only : bilinx
      use netcdf

      implicit none
!
      integer , parameter :: ilon = 240 , jlat = 121
      integer , parameter :: idtbc = 6
!
      integer :: i , it , j , nday , nhour , nmo , nyear
      real(4) , dimension(jlat) :: lati
      real(4) , dimension(ilon) :: loni
      integer :: idate , ierastart , ierrec , nsteps
      real(4) , dimension(ilon,jlat) :: sst
      character(256) :: inpfile
!
      if ( ssttyp /= 'ERSST' .and. ssttyp /= 'ERSKT' ) then
        write (*,*) 'PLEASE SET right SSTTYP in regcm.in'
        write (*,*) 'Supported types are ERSST ERSKT'
        stop
      end if

      nsteps = idatediff(globidate2,globidate1)/idtbc + 1

      write (*,*) 'GLOBIDATE1 : ' , globidate1
      write (*,*) 'GLOBIDATE2 : ' , globidate2
      write (*,*) 'NSTEPS     : ' , nsteps

      call open_sstfile(globidate1)
 
!     ******    SET UP LONGITUDES AND LATITUDES FOR SST DATA
      do i = 1 , ilon
        loni(i) = float(i-1)*1.5
      end do
      do j = 1 , jlat
        lati(j) = -90. + 1.5*float(j-1)
      end do
 
      idate = globidate1
      ierastart = 1989010100
      do it = 1 , nsteps

        ierrec = idatediff(idate,ierastart)/idtbc+1

        if ( ssttyp=='ERSST' ) then
          inpfile = trim(inpglob)//'/SST/sstERAIN.1989-2009.nc'
          call sst_erain(ierrec,ilon,jlat,sst,inpfile)
        else if ( ssttyp=='ERSKT' ) then
          inpfile = trim(inpglob)//'/SST/tskinERAIN.1989-2009.nc'
          call skt_erain(ierrec,ilon,jlat,sst,inpfile)
        else
        end if

        call split_idate(idate, nyear, nmo, nday, nhour)
 
        call bilinx(sst,sstmm,xlon,xlat,loni,lati,ilon,jlat,iy,jx,1)
        print * , 'XLON,XLAT,SST=' , xlon(1,1) , xlat(1,1) , sstmm(1,1)
 
!       ******           WRITE OUT SST DATA ON MM4 GRID
        call writerec(idate,.false.)
        print * , 'WRITING OUT MM4 SST DATA:' , nmo , nyear

        call addhours(idate, idtbc)

      end do
 
      end subroutine sst_ersst
!
!-----------------------------------------------------------------------
!
      subroutine skt_erain(it,ilon,jlat,sst,pathaddname)
      use netcdf
      implicit none
!
      integer :: it , ilon , jlat
      character(256) :: pathaddname
      intent (in) it , ilon , jlat , pathaddname
      real(4) , dimension(ilon,jlat) :: sst
      intent (out) :: sst
!
      integer :: i , j , n
      character(5) :: varname
      integer(2) , dimension(ilon,jlat) :: work
      integer :: istatus
!
      integer , dimension(10) , save :: icount , istart
      integer , save :: inet , ivar
      real(8) , save :: xadd , xscale , xmiss
!
!     This is the latitude, longitude dimension of the grid to be read.
!     This corresponds to the lat and lon dimension variables in the
!     netCDF file.
!
!     The data are packed into short integers (INTEGER*2).  The array
!     work will be used to hold the packed integers. The array 'sst'
!     will contain the unpacked data.
!
!     DATA ARRAY AND WORK ARRAY
!
      data varname/'skt'/
!
      if ( it==1 ) then
        istatus = nf90_open(pathaddname,nf90_nowrite,inet)
        if ( istatus/=nf90_noerr ) then
          write ( 6,* ) 'Cannot open input file ', trim(pathaddname)
          stop 'INPUT FILE OPEN ERROR'
        end if
        istatus = nf90_inq_varid(inet,varname,ivar)
        if ( istatus/=nf90_noerr ) then
          write ( 6,* ) 'Cannot find variable ', varname,               &
               &        ' in input file ', trim(pathaddname)
          stop 'INPUT FILE ERROR'
        end if
        istatus = nf90_get_att(inet,ivar,'scale_factor',xscale)
        istatus = nf90_get_att(inet,ivar,'add_offset',xadd)
        istatus = nf90_get_att(inet,ivar,'_FillValue',xmiss)
        istart(1) = 1
        istart(2) = 1
        icount(1) = 240
        icount(2) = 121
        do n = 4 , 10
          istart(n) = 0
          icount(n) = 0
        end do
      end if
! 
      istart(3) = it
      icount(3) = 1
      istatus = nf90_get_var(inet,ivar,work,istart,icount)
      if ( istatus/=nf90_noerr ) then
        write ( 6,* ) 'Cannot get ', varname, ' from file'
        write ( 6,* ) nf90_strerror(istatus)
        stop 'ERROR READ SST'
      end if
!
      do j = 1 , jlat
        do i = 1 , ilon
          if (work(i,j)/=xmiss) then
            sst(i,jlat+1-j) = real(dble(work(i,j))*xscale + xadd)
          else
            sst(i,jlat+1-j)=-9999.0
          end if
        end do
      end do
!
      end subroutine skt_erain
!
!-----------------------------------------------------------------------
!
      subroutine sst_erain(it,ilon,jlat,sst,pathaddname)
      use netcdf
      implicit none
!
      integer :: it , ilon , jlat
      character(256) :: pathaddname
      intent (in) it , ilon , jlat , pathaddname
      real(4) , dimension(ilon,jlat) :: sst
      intent (out) :: sst
!
      integer :: i , j , n
      character(5) :: varname
      integer(2) , dimension(ilon,jlat) :: work
      integer :: istatus
!
      integer , save :: inet , ivar
      real(8) , save :: xadd , xscale , xmiss
      integer , dimension(10) , save :: icount , istart
!
!     This is the latitude, longitude dimension of the grid to be read.
!     This corresponds to the lat and lon dimension variables in the
!     netCDF file.
!
!     The data are packed into short integers (INTEGER*2).  The array
!     work will be used to hold the packed integers. The array 'sst'
!     will contain the unpacked data.
!
!     DATA ARRAY AND WORK ARRAY
!
      data varname/'sst'/
!
      if ( it==1 ) then
        istatus = nf90_open(pathaddname,nf90_nowrite,inet)
        if ( istatus/=nf90_noerr ) then
          write ( 6,* ) 'Cannot open input file ', trim(pathaddname)
          stop 'INPUT FILE OPEN ERROR'
        end if
        istatus = nf90_inq_varid(inet,varname,ivar)
        if ( istatus/=nf90_noerr ) then
          write ( 6,* ) 'Cannot find variable ', varname,               &
               &        ' in input file ', trim(pathaddname)
          stop 'INPUT FILE ERROR'
        end if
        istatus = nf90_get_att(inet,ivar,'scale_factor',xscale)
        istatus = nf90_get_att(inet,ivar,'add_offset',xadd)
        istatus = nf90_get_att(inet,ivar,'_FillValue',xmiss)
        istart(1) = 1
        istart(2) = 1
        icount(1) = 240
        icount(2) = 121
        do n = 4 , 10
          istart(n) = 0
          icount(n) = 0
        end do
      end if
!
      istart(3) = it
      icount(3) = 1
      istatus = nf90_get_var(inet,ivar,work,istart,icount)
      if ( istatus/=nf90_noerr ) then
        write ( 6,* ) 'Cannot get ', varname, ' from file'
        write ( 6,* ) nf90_strerror(istatus)
        stop 'ERROR READ SST'
      end if
!
      do j = 1 , jlat
        do i = 1 , ilon
          if (work(i,j)/=xmiss) then
            sst(i,jlat+1-j) = real(dble(work(i,j))*xscale + xadd)
          else
            sst(i,jlat+1-j)=-9999.0
          end if
        end do
      end do
!
      end subroutine sst_erain
!
      end module mod_sst_ersst

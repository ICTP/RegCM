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

      module mod_sst_ccsm

      use m_realkinds
      use m_die
      use m_stdio
      use m_mall
      use m_zeit

      contains

      subroutine sst_ccsm
!
!*******************************************************************************
!
! This is a package of subroutines to read CCSM SST 1x1 degree data in
! NETCDF format and interpolate SST at RCM grid
! Written By Moetasim Ashfaq Dec-2005 @ PURDUE.EDU
!
!******************************************************************************
!******************************************************************************
!       DATA PREPARATION
!       Dataset required to use this code can be preapred using NCO utilitiles
!       such as NCKS, NCRCAT etc.
!       We need top level of TEMP for SSTs which can be extarcted as following:
!       ncks -v time,TEMP -d z_t,0 input.nc output.nc
!       Files can be further concatenated using 'ncrcat'
!       Finally, the POP grid can be converted into lat/lon grid at 1x1 degree
!       resolution  using PopLatLon function in NCL
!******************************************************************************
!     NAMING CONVENTION (Global Data File)
!       ccsm_mn.sst.nc  for TEMP
!     PATH /DATA/SST/
!
!******************************************************************************
      use mod_sst_grid
      use mod_date
      use mod_interp , only : bilinx
      use m_die
      use m_mall
      use m_zeit

      implicit none
!
! PARAMETER definitions
!
      integer , parameter :: ilon = 360 , jlat = 180
!
! Local variables
!
      real(sp) , dimension(jlat) :: glat
      real(sp) , dimension(ilon) :: glon
      real(sp) , dimension(ilon,jlat) :: sst
      integer :: idate
      integer :: i , idatef , idateo , j , k , ludom , lumax , &
               & nday , nmo , nyear , nho , iv , nsteps
      integer , dimension(20) :: lund
      character(256) :: inpfile
!
      call zeit_ci('sst_ccsm')
      do i = 1 , ilon
        glon(i) = 0.5 + float(i-1)
      end do
      do j = 1 , jlat
        glat(j) = -89.5 + 1.*float(j-1)
      end do
 
      idateo = imonfirst(globidate1)
      if (lfhomonth(globidate1)) then
        idateo = iprevmon(globidate1)
      end if
      idatef = imonfirst(globidate2)
      if (idatef < globidate2) then
        idatef = inextmon(idatef)
      end if
      nsteps = imondiff(idatef,idateo) + 1
 
      call open_sstfile(idateo)
 
      idate = idateo
      do k = 1 , nsteps

        call split_idate(idate, nyear, nmo, nday, nho)
 
        inpfile = trim(inpglob)//'/SST/ccsm_mn.sst.nc'

        call ccsm_sst(idate,idateo,ilon,jlat,sst,inpfile)
        call bilinx(sst,sstmm,xlon,xlat,glon,glat,ilon,jlat,iy,jx,1)

        write (stdout,*) &
          'XLON,XLAT,SST=' , xlon(1,1) , xlat(1,1) , sstmm(1,1)+273.
 
        do j = 1 , jx
          do i = 1 , iy
            if ( sstmm(i,j)<-5000 .and.                                 &
               & (lu(i,j)>13.5 .and. lu(i,j)<15.5) ) then
              do iv = 1 , 20
                lund(iv) = 0.0
              end do
              lund(nint(lu(i-1,j-1))) = lund(nint(lu(i-1,j-1))) + 2
              lund(nint(lu(i-1,j))) = lund(nint(lu(i-1,j))) + 3
              lund(nint(lu(i-1,j+1))) = lund(nint(lu(i-1,j+1))) + 2
              lund(nint(lu(i,j-1))) = lund(nint(lu(i,j-1))) + 3
              lund(nint(lu(i,j+1))) = lund(nint(lu(i,j+1))) + 3
              lund(nint(lu(i+1,j-1))) = lund(nint(lu(i+1,j-1))) + 2
              lund(nint(lu(i+1,j))) = lund(nint(lu(i+1,j))) + 3
              lund(nint(lu(i+1,j+1))) = lund(nint(lu(i+1,j+1))) + 2
              ludom = 18
              lumax = 0
              do iv = 1 , 20
                if ( iv<=13 .or. iv>=16 ) then
                  if ( lund(iv)>lumax ) then
                    ludom = k
                    lumax = lund(iv)
                  end if
                end if
              end do
              lu(i,j) = float(ludom)
              write (stdout,*) ludom , sstmm(i,j)
            end if
            if ( sstmm(i,j)>-100. ) then
              sstmm(i,j) = sstmm(i,j) + 273.15
            else
              sstmm(i,j) = -9999.
            end if
          end do
        end do

        call writerec(idate,.false.)
        write (stdout,*) 'WRITEN OUT SST DATA : ' , idate

        idate = inextmon(idate)

      end do
 
      call zeit_co('sst_ccsm')

      end subroutine sst_ccsm
!
!-----------------------------------------------------------------------
!
!     Subroutine to read required records from SST data file
!
      subroutine ccsm_sst(idate,idate0,ilon,jlat,sst,pathaddname)

      use netcdf

      implicit none

      integer , intent (in) :: idate , idate0
      integer , intent (in) :: ilon , jlat
      character(len=256) ,intent(in) :: pathaddname
      real(sp) , dimension(ilon, jlat) , intent(out) :: sst

      integer, dimension(12) :: ndays
      character(len=4), dimension(2) :: varname
      integer, allocatable ::  work1(:)
      real(sp) , dimension (ilon , jlat) :: work2
      real(sp) :: imisng

      integer :: nyear , month
      integer :: inet1
      integer, dimension(10) :: istart , icount , istartt , icountt
      integer, dimension(2) :: ivar2
      integer :: it , icode , i , j , npos , nrec
      integer :: latid , lonid , timid
      integer :: latlen , lonlen , timlen
      logical :: there
      integer :: istatus

      data ndays/31,59,90,120,151,181,212,243,273,304,334,365/
      data varname/'time','TEMP'/
      
      call zeit_ci('read_ccsm')
      nyear = idate/1000000
      month = idate/10000 - nyear*100
      
      if (idate == idate0) then
         inquire(file=pathaddname,exist=there)
         if (.not.there) then
           call die('ccsm_sst',trim(pathaddname)//' is not available',1)
         endif
         istatus = nf90_open(pathaddname,nf90_nowrite,inet1)
         if ( istatus/=nf90_noerr ) then
           call die('ccsm_sst','Error opening '//trim(pathaddname),1, &
                    nf90_strerror(istatus),istatus)
         end if
         
         write(stdout,*) inet1 , trim(pathaddname) , icode
      endif  
!     GET DIMENSION IDs
      istatus = nf90_inq_dimid(inet1,'lat',latid)
      if ( istatus/=nf90_noerr ) then
        call die('ccsm_sst','Error '//trim(pathaddname)//' dim lat',1,  &
                 nf90_strerror(istatus),istatus)
      end if
      istatus = nf90_inq_dimid(inet1,'lon',lonid)
      if ( istatus/=nf90_noerr ) then
        call die('ccsm_sst','Error '//trim(pathaddname)//' dim lon',1,  &
                 nf90_strerror(istatus),istatus)
      end if
      istatus = nf90_inq_dimid(inet1,'time',timid)
      if ( istatus/=nf90_noerr ) then
        call die('ccsm_sst','Error '//trim(pathaddname)//' dim time',1, &
                 nf90_strerror(istatus),istatus)
      end if

!     GET DIMENSION LENGTHS
      istatus = nf90_inquire_dimension(inet1,latid,len=latlen)
      if ( istatus/=nf90_noerr ) then
        call die('ccsm_sst','Error '//trim(pathaddname)//' dim lat',1,  &
                 nf90_strerror(istatus),istatus)
      end if
      istatus = nf90_inquire_dimension(inet1,lonid,len=lonlen)
      if ( istatus/=nf90_noerr ) then
        call die('ccsm_sst','Error '//trim(pathaddname)//' dim lon',1,  &
                 nf90_strerror(istatus),istatus)
      end if
      istatus = nf90_inquire_dimension(inet1,timid,len=timlen)
      if ( istatus/=nf90_noerr ) then
        call die('ccsm_sst','Error '//trim(pathaddname)//' dim time',1, &
                 nf90_strerror(istatus),istatus)
      end if
      allocate(work1(timlen))
      if (istatus /= 0) call die('ccsm_sst','allocate work1',istatus)
      call mall_mci(work1,'mod_sst_ccsm')
      
!     MAKE SURE THAT SST DATA IS AT 1X1 DEGREE
      if(latlen /= jlat .or. lonlen /= ilon) then
        write (stderr,*) 'DIMENSIONS DO NOT MATCH'
        write (stderr,*) 'No. of LON in SST file =',lonlen
        write (stderr,*) 'No. of LON in 1x1 degree gloabl grid =',ilon
        write (stderr,*) 'No. of LAT in SST file =',latlen
        write (stderr,*) 'No. of LON in 1x1 degree gloabl grid =',jlat
        call die('ccsm_sst')
      endif
!     GET VARIABLE IDs
      istatus = nf90_inq_varid(inet1,varname(1),ivar2(1))
      if ( istatus/=nf90_noerr ) then
        call die('ccsm_sst','Error '//trim(pathaddname)//':'// &
                 varname(1),1,nf90_strerror(istatus),istatus)
      end if
      istatus = nf90_inq_varid(inet1,varname(2),ivar2(2))
      if ( istatus/=nf90_noerr ) then
        call die('ccsm_sst','Error '//trim(pathaddname)//':'// &
                 varname(2),1,nf90_strerror(istatus),istatus)
      end if
!     GET MISSING DATA VALUE
      istatus = nf90_get_att(inet1,ivar2(2),'_FillValue',imisng)
      if ( istatus/=nf90_noerr ) then
        call die('ccsm_sst','Error '//trim(pathaddname)//':'// &
           varname(2)//':_FillValue',1,nf90_strerror(istatus),istatus)
      end if
!     GET TIME VALUES
      istartt(1) = 1
      icountt(1) = timlen
      istatus = nf90_get_var(inet1,ivar2(1),work1,istartt,icountt)
      if ( istatus/=nf90_noerr ) then
        call die('ccsm_sst','Error '//trim(pathaddname)//':'// &
            varname(1)//' read',1,nf90_strerror(istatus),istatus)
      end if
      
!     CHECK FOR THE REQUIRED RECORD IN DATA FILE  
      npos = (nyear - 1000) * 365 + ndays(month)
      if (npos < work1(1) .or. npos > work1(timlen)) then
        write (stderr,*) 'Error in finding SST data for', &
                         (idate-100)/10000
        write (stderr,*) 'Required NREC = ', npos
        call die('ccsm_sst')
      end if

      i = 1
      do
        if (work1(i) == npos) then
          nrec=i
          exit
        end if
        i = i + 1
      end do
 
      it = nrec
      icount(1) = ilon
      icount(2) = jlat
      icount(3) = 1
      istart(1) = 1
      istart(2) = 1
      istart(3) = 1
      istart(3) = it
      istatus = nf90_get_var(inet1,ivar2(2),work2,istart,icount)
      if ( istatus/=nf90_noerr ) then
        call die('ccsm_sst','Error '//trim(pathaddname)//':'// &
            varname(2)//' read',1,nf90_strerror(istatus),istatus)
      end if
      do j = 1 , jlat
         do i = 1 , ilon
            if (work2(i,j) > (imisng+10) .and. work2(i,j) < 10000.0) then
               sst(i,j) = work2(i,j)
            else
               sst(i,j) = -9999.
            end if
         end do
      end do

      deallocate(work1)
      call mall_mco(work1,'mod_sst_ccsm')

      call zeit_co('read_ccsm')

      end subroutine ccsm_sst
!
      end module mod_sst_ccsm

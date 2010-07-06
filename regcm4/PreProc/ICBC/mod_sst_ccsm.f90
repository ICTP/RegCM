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

      implicit none
!
! PARAMETER definitions
!
      integer , parameter :: ilon = 360 , jlat = 180
!
! Local variables
!
      real(4) , dimension(jlat) :: glat
      real(4) , dimension(ilon) :: glon
      real(4) , dimension(ilon,jlat) :: sst
      integer :: idate
      integer :: i , idatef , idateo , j , k , ludom , lumax , &
               & nday , nmo , nyear , nho , iv , nsteps
      integer , dimension(20) :: lund
      character(256) :: inpfile
!
      do i = 1 , ilon
        glon(i) = 0.5 + float(i-1)
      end do
      do j = 1 , jlat
        glat(j) = -89.5 + 1.*float(j-1)
      end do
 
      idateo = imonfirst(globidate1)
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

        print * , 'XLON,XLAT,SST=' , xlon(1,1) , xlat(1,1) , sstmm(1,1) &
            & + 273.
 
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
              print * , ludom , sstmm(i,j)
            end if
            if ( sstmm(i,j)>-100. ) then
              sstmm(i,j) = sstmm(i,j) + 273.15
            else
              sstmm(i,j) = -9999.
            end if
          end do
        end do

        call writerec(idate,.false.)
        print * , 'WRITEN OUT SST DATA : ' , idate

        idate = inextmon(idate)

      end do
 
      return

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
      real(4) , dimension(ilon, jlat) , intent(out) :: sst

      integer, dimension(12) :: ndays
      character(len=4), dimension(2) :: varname
      integer, allocatable ::  work1(:)
      real(4) , dimension (ilon , jlat) :: work2
      real(4) :: imisng

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
      
      nyear = idate/1000000
      month = idate/10000 - nyear*100
      
      if (idate == idate0) then
         inquire(file=pathaddname,exist=there)
         if (.not.there) then
            print *, trim(pathaddname),' is not available'
            stop
         endif
         istatus = nf90_open(pathaddname,nf90_nowrite,inet1)
         if ( istatus/=nf90_noerr ) then
           write ( 6,*) 'Error opening ', trim(pathaddname)
           stop 'ERROR OPEN FILE'
         end if
         
         write(*,*) inet1 , trim(pathaddname) , icode
      endif  
!     GET DIMENSION IDs
      istatus = nf90_inq_dimid(inet1,'lat',latid)
      istatus = nf90_inq_dimid(inet1,'lon',lonid)
      istatus = nf90_inq_dimid(inet1,'time',timid)

!     GET DIMENSION LENGTHS
      istatus = nf90_inquire_dimension(inet1,latid,len=latlen)
      istatus = nf90_inquire_dimension(inet1,lonid,len=lonlen)
      istatus = nf90_inquire_dimension(inet1,timid,len=timlen)
      allocate(work1(timlen))
      
!     MAKE SURE THAT SST DATA IS AT 1X1 DEGREE
      if(latlen /= jlat .or. lonlen /= ilon) then
         print*,'DIMENSIONS DO NOT MATCH'
         print*,'No. of LON in SST file =',lonlen
         print*,'No. of LON in 1x1 degree gloabl grid =',ilon
         print*,'No. of LAT in SST file =',latlen
         print*,'No. of LON in 1x1 degree gloabl grid =',jlat
         STOP   'Check SST data file' 
      endif
!     GET VARIABLE IDs
      istatus = nf90_inq_varid(inet1,varname(1),ivar2(1))
      istatus = nf90_inq_varid(inet1,varname(2),ivar2(2))
!     GET MISSING DATA VALUE
      istatus = nf90_get_att(inet1,ivar2(2),'_FillValue',imisng)
!     GET TIME VALUES
      istartt(1) = 1
      icountt(1) = timlen
      istatus = nf90_get_var(inet1,ivar2(1),work1,istartt,icountt)
      
!     CHECK FOR THE REQUIRED RECORD IN DATA FILE  
      npos = (nyear - 1000) * 365 + ndays(month)
      i = 0
      print *,  npos
 10   continue
      i = i + 1
      if (npos < work1(i) .or. npos > work1(timlen)) then
         print *, 'Error in finding SST data for',(idate-100)/10000
         print *, 'Required NREC=',npos
         stop    'Check SST data file' 
      else if (work1(i) == npos) then
         nrec=i
         go to 20
      end if
      go to 10
 
 20   it = nrec
      icount(1) = ilon
      icount(2) = jlat
      icount(3) = 1
      istart(1) = 1
      istart(2) = 1
      istart(3) = 1
      istart(3) = it
      istatus = nf90_get_var(inet1,ivar2(2),work2,istart,icount)
      do j = 1 , jlat
         do i = 1 , ilon
            if (work2(i,j) > (imisng+10) .and. work2(i,j) < 10000.0) then
               sst(i,j) = work2(i,j)
            else
               sst(i,j) = -9999.
            end if
         end do
      end do

      end subroutine ccsm_sst
!
      end module mod_sst_ccsm

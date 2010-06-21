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

      module mod_sst_1deg

      contains

      subroutine sst_1deg

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Comments on dataset sources and location:                          c
!                                                                    c
! GISST2.3b    UKMO SST (Rayner et al 1996), 1 degree                c
!              from UKMO DATA archive (http://www.badc.rl.ac.uk/)    c
!              and reformed as direct-accessed binary GrADS format   c
!              in file GISST_187101_200209                           c
!              ML= 1 is-179.5; ML= 2 is-178.5; => ML=360 is 179.5E   c
!              NL= 1 is -89.5; NL= 2 is -88.5; => NL=180 is  89.5    c
!              see the GrADS control file for details.               c
!                                                                    c
! OISST        from CAC Optimal Interpolation dataset.               c
!              in the original netCDF format.                        c
!              ML= 1 is   0.5; ML= 2 is   1.5; => ML=360 is 359.5E   c
!              NL= 1 is -89.5; NL= 2 is -88.5; => NL=180 is  89.5    c
!                                                                    c
! OI2ST        both SST and SeaIce in the original netCDF format.    c
!                                                                    c
! OI_WK        weekly OISST in the original netCDF format.           c
!                                                                    c
! OI2WK        weekly OISST and SeaIce in the original netCDF format.c
!                                                                    c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      use mod_sst_grid
      use mod_interp , only : bilinx
      use mod_printl

      implicit none
!
! PARAMETER definitions
!
      integer , parameter :: ilon = 360 , jlat = 180
!
! Local variables
!
      real(4) , dimension(ilon,jlat) :: sst , ice
      integer :: idate , idate0 , kend , kstart
      integer , dimension(427+1097) :: wkday
      integer :: i , idatef , idateo , j , k , ludom , lumax , mrec ,   &
               & nday , nmo , nrec , nyear , nsteps
      real(4) , dimension(jlat) :: lati
      real(4) , dimension(ilon) :: loni
      integer , dimension(25) :: lund
      character(256) :: sstfile , inpfile
      logical :: there
!
      kstart = 0
      kend = 0
!
      if ( ssttyp=='GISST' ) then
        if ( globidate1<1947121512 .or. globidate2>2002091512 ) then
          print * , 'GISST data required are not available'
          print * , 'IDATE1, IDATE2 = ' , globidate1 , globidate2
          stop
        end if
        open (11,file=trim(inpglob)//'/SST/GISST_194712_200209',        &
             &form='unformatted',recl=ilon*jlat*ibyte,access='direct',  &
            & status='old',err=100)
      else if ( ssttyp=='OISST' .or. ssttyp=='OI_NC' .or.               &
            &   ssttyp=='OI2ST' ) then
        if ( globidate1<1981121512 .or. globidate2<1981121512 ) then
          print * , 'OISST data required are not available'
          print * , 'IDATE1, IDATE2 = ' , globidate1 , globidate2
          stop
        end if
        inquire (file=trim(inpglob)//'/SST/sst.mnmean.nc',exist=there)
        if ( .not.there ) print * , 'sst.mnmean.nc is not available' ,  &
                               &' under ',trim(inpglob),'/SST/'
        if ( ssttyp=='OI2ST' ) then
          inquire (file=trim(inpglob)//'/SST/icec.mnmean.nc',           &
              &    exist=there)
          if ( .not. there )                                            &
            & print *, 'icec.mnmean.nc is not available' ,              &
                   &' under ',trim(inpglob),'/SST/'
        end if
      else if ( ssttyp=='OI_WK' .or. ssttyp=='OI2WK' ) then
        if ( globidate1<1981110100 .or. globidate2<1981110106 ) then
          print * , 'OI_WK (or OI2WK) data required are not available'
          print * , 'IDATE1, IDATE2 = ' , globidate1 , globidate2
          stop
        end if
        inquire (file=trim(inpglob)//'/SST/sst.wkmean.1981-1989.nc',    &
             &   exist=there)
        if ( .not.there ) print * ,                                     &
                             &'sst.wkmean.1981-1989.nc is not available'&
                            & , ' under ',trim(inpglob),'/SST/'
        inquire (file=trim(inpglob)//'/SST/sst.wkmean.1990-present.nc', &
               & exist=there)
        if ( .not.there ) print * ,                                     &
                          &'sst.wkmean.1990-present.nc is not available'&
                         & , ' under ',trim(inpglob),'/SST/'
        call headwk(wkday)
        if ( ssttyp=='OI2WK' ) then
          inquire (file=trim(inpglob)//'/SST/icec.wkmean.1981-1989.nc', &
                & exist=there)
          if ( .not.there ) then
            print *, 'icec.wkmean.1981-1989.nc is not available',       &
                &    ' under ',trim(inpglob),'/SST/'
          end if
          inquire (file=trim(inpglob)//                                 &
                 '/SST/icec.wkmean.1990-present.nc',exist=there)
          if ( .not.there ) then
            print *, 'icec.wkmean.1990-present.nc is not available',    &
                &    ' under ',trim(inpglob),'/SST/'
          end if
        end if
      else
        write (*,*) 'PLEASE SET right SSTTYP in regcm.in'
        write (*,*) 'Supported are GISST OISST OI_NC OI2ST OI_WK OI2WK'
        stop
      end if

      write (sstfile,99001) trim(dirglob), pthsep, trim(domname),       &
             & '_SST.RCM'
      open (21,file=sstfile,form='unformatted',status='replace')
      write (sstfile,99001) trim(dirglob), pthsep, trim(domname),       &
          & '_RCM_SST.dat'
      open (25,file=sstfile,status='unknown',form='unformatted',        &
          & recl=iy*jx*ibyte,access='direct')

      if ( igrads==1 ) then
        write (sstfile,99001) trim(dirglob), pthsep, trim(domname),     &
           &  '_RCM_SST.ctl'
        open (31,file=sstfile,status='replace')
        write (31,'(a,a,a)') 'dset ^',trim(domname),'_RCM_SST.dat'
      end if
      if ( ssttyp=='OI2ST' .or. ssttyp=='OI2WK' ) then
        write (sstfile,99001) trim(dirglob), pthsep, trim(domname),     &
           &  '_RCM_ICE.dat'
        open (26,file=sstfile,status='unknown',form='unformatted',      &
          & recl=iy*jx*ibyte,access='direct')
        if ( igrads==1 ) then
          write (sstfile,99001) trim(dirglob), pthsep, trim(domname),   &
            & '_RCM_ICE.ctl'
          open ( 32,file=sstfile,status='replace')
          write (32,'(a,a,a)') 'dset ^',trim(domname),'_RCM_ICE.dat'
        end if
      end if
!#####
      if ( ssttyp/='OI_WK' .and. ssttyp/='OI2WK' ) then
!#####
        idate = globidate1/10000
        if ( idate-(idate/100)*100==1 ) then
          idate = idate - 89
        else
          idate = idate - 1
        end if
        idateo = idate
        idate0 = idateo*10000 + 100
        idate = globidate2/10000
        if ( idate-(idate/100)*100==12 ) then
          idate = idate + 89
        else
          idate = idate + 1
        end if
        idatef = idate
        nsteps = (idatef/100-idateo/100)*12 + (idatef-(idatef/100)*100) &
               & - (idateo-(idateo/100)*100) + 1
        write (*,*) 'GLOBIDATE1 : ' , globidate1
        write (*,*) 'GLOBIDATE2 : ' , globidate2
        write (*,*) 'NSTEPS     : ' , nsteps

        call setup_sstfile(idateo,nsteps)
!#####
      else
!#####
        idate = globidate1/100
        do k = 427 + 1097 , 1 , -1
          if ( wkday(k)<=idate ) then
            kstart = k
            exit
          end if
        end do
        idate = globidate2/100
        do k = 1 , 427 + 1097
          if ( wkday(k)>idate ) then
            kend = k
            exit
          end if
        end do
        idateo = wkday(kstart)
        idate0 = wkday(kstart)*100
        idatef = wkday(kend)
        print * , globidate1 , globidate2 , idateo , idatef ,           &
                &  kend - kstart + 1
        call setup_sst_ice_file(idateo,kend-kstart+1)
!#####
      end if
!#####
      mrec = 0
 
!     ******    SET UP LONGITUDES AND LATITUDES FOR SST DATA
      do i = 1 , ilon
        loni(i) = .5 + float(i-1)
      end do
      do j = 1 , jlat
        lati(j) = -89.5 + 1.*float(j-1)
      end do
 
!#####
      if ( ssttyp/='OI_WK' .and. ssttyp/='OI2WK' ) then
!#####
!       ****** OISST SST DATA, 1 Deg data, AVAILABLE FROM 12/1981 TO
!       PRESENT ****** GISST SST DATA, 1 Deg data, AVAILABLE FROM
!       12/1947 TO 9/2002
        idate = idateo
        do while ( idate<=idatef )
          nyear = idate/100
          nmo = idate - nyear*100
          if ( ssttyp=='GISST' ) then
            nrec = (nyear-1947)*12 + nmo - 11
            read (11,rec=nrec) sst
          else if ( ssttyp=='OISST' .or. ssttyp=='OI_NC' .or.           &
                    ssttyp=='OI2ST') then
            write (*,*) idate*10000 + 100 , idate0
            inpfile = trim(inpglob)//'/SST/sst.mnmean.nc'
            call sst_mn(idate*10000+100,idate0,ilon,jlat,sst,inpfile)
            if ( ssttyp=='OI2ST' ) then
              inpfile = trim(inpglob)//'/SST/icec.mnmean.nc'
              call ice_mn(idate*10000+100,idate0,ilon,jlat,ice,inpfile)
            end if
          else
          end if
 
!         ******           PRINT OUT DATA AS A CHECK
          if ( nmo==1 ) call printl(sst,ilon,jlat)
 
          call bilinx(sst,sstmm,xlon,xlat,loni,lati,ilon,jlat,iy,jx,1)
          print * , 'XLON,XLAT,SST=' , xlon(1,1) , xlat(1,1) ,          &
              & sstmm(1,1)
 
          if ( ssttyp=='OI2ST' ) call bilinx(ice,icemm,xlon,xlat,       &
              & loni,lati,ilon,jlat,iy,jx,1)
          do j = 1 , jx
            do i = 1 , iy
              if ( sstmm(i,j)<-5000. .and. (lu(i,j)>13.5 .and.          &
                &  lu(i,j)<15.5) ) then
                do k = 1 , 20
                  lund(k) = 0.0
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
                do k = 1 , 20
                  if ( k<=13 .or. k>=16 ) then
                    if ( lund(k)>lumax ) then
                      ludom = k
                      lumax = lund(k)
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
 
!         ******           WRITE OUT SST DATA ON MM4 GRID
          nday = 1
          if ( ssttyp/='OI2ST' ) then
            write (21) nday , nmo , nyear , sstmm
          else
            write (21) nday , nmo , nyear , sstmm , icemm
          end if
          print * , 'WRITING OUT MM4 SST DATA:' , nmo , nyear
          idate = idate + 1
          if ( nmo==12 ) idate = idate + 88
          mrec = mrec + 1
          write (25,rec=mrec) ((sstmm(i,j),j=1,jx),i=1,iy)
          if ( ssttyp=='OI2ST' )                                        &
           &    write (26,rec=mrec) ((icemm(i,j),j=1,jx),i=1,iy)
        end do
        write (10,rec=4) ((lu(i,j),j=1,jx),i=1,iy)
!#####
      else
!#####
        do k = kstart , kend
          idate = wkday(k)
          nyear = idate/10000
          nmo = idate/100 - nyear*100
          nday = mod(idate,100)
          write (*,*) idate*100 , idate0 , k
          if ( idate<19891231 ) then
            inpfile = trim(inpglob)//'/SST/sst.wkmean.1981-1989.nc'
          else
            inpfile = trim(inpglob)//'/SST/sst.wkmean.1990-present.nc'
          end if
          call sst_wk(idate*100,idate0,k,ilon,jlat,sst,inpfile)
 
          call bilinx(sst,sstmm,xlon,xlat,loni,lati,ilon,jlat,iy,jx,1)
          print * , 'XLON,XLAT,SST=' , xlon(1,1) , xlat(1,1) ,          &
              & sstmm(1,1)
          if ( ssttyp=='OI2WK') then
            if ( idate<19891231 ) then
              inpfile = trim(inpglob)//'/SST/icec.wkmean.1981-1989.nc'
            else
             inpfile = trim(inpglob)//'/SST/icec.wkmean.1990-present.nc'
            end if
            call ice_wk(idate*100,idate0,k,ilon,jlat,sst,inpfile)
            call bilinx(ice,icemm,xlon,xlat,loni,lati,ilon,jlat,iy,jx,1)
          end if 

          do j = 1 , jx
            do i = 1 , iy
              if ( sstmm(i,j)>-100. ) then
                sstmm(i,j) = sstmm(i,j) + 273.15
              else
                sstmm(i,j) = -9999.
              end if
            end do
          end do
 
!         ******           WRITE OUT SST DATA ON MM4 GRID
          if(ssttyp.eq.'OI_WK') then
             write (21) nday , nmo , nyear , sstmm
          else
             write (21) nday , nmo , nyear , sstmm,icemm
          endif
          idate = wkday(k)
          print * , 'WRITING OUT MM4 SST DATA:' , nmo , nyear , idate
          mrec = mrec + 1
          write (25,rec=mrec) ((sstmm(i,j),j=1,jx),i=1,iy)
          if(ssttyp.eq.'OI2WK') &
          write (26,rec=mrec) ((icemm(i,j),j=1,jx),i=1,iy)
        end do
 
!#####
      end if
!#####
      return

 100  continue
      print * , 'ERROR OPENING GISST FILE'
      stop '4810 IN PROGRAM RDSST'
 200  continue
      print * , 'ERROR OPENING DOMAIN HEADER FILE'
      stop '4830 IN PROGRAM RDSST'

99001 format (a,a,a,a)
      end subroutine sst_1deg
!
!-----------------------------------------------------------------------
!
      subroutine headwk(wkday)
      implicit none
!
! Dummy arguments
!
      integer , intent(out) , dimension(427+1097) :: wkday
!
! Local variables
!
      integer :: i , mday , month , myear
!
      wkday(1) = 19811029
      do i = 2 , 427
        wkday(i) = wkday(i-1) + 7
        myear = wkday(i)/10000
        month = wkday(i)/100 - myear*100
        mday = mod(wkday(i),10000) - month*100
        if ( month==1 .or. month==3 .or. month==5 .or. month==7 .or.    &
           & month==8 .or. month==10 ) then
          if ( mday>31 ) then
            mday = mday - 31
            month = month + 1
          end if
        else if ( month==12 ) then
          if ( mday>31 ) then
            mday = mday - 31
            month = 1
            myear = myear + 1
          end if
        else if ( month==4 .or. month==6 .or. month==9 .or. month==11 ) &
                & then
          if ( mday>30 ) then
            mday = mday - 30
            month = month + 1
          end if
        else if ( mod(myear,4)/=0 ) then
          if ( mday>28 ) then
            mday = mday - 28
            month = month + 1
          end if
        else if ( mod(myear,400)==0 ) then
          if ( mday>29 ) then
            mday = mday - 29
            month = month + 1
          end if
        else if ( mod(myear,100)==0 ) then
          if ( mday>28 ) then
            mday = mday - 28
            month = month + 1
          end if
        else if ( mday>29 ) then
          mday = mday - 29
          month = month + 1
        else
        end if
        wkday(i) = myear*10000 + month*100 + mday
      end do
!
      wkday(428) = 19891231
      do i = 429 , 427 + 1097
        wkday(i) = wkday(i-1) + 7
        myear = wkday(i)/10000
        month = wkday(i)/100 - myear*100
        mday = mod(wkday(i),10000) - month*100
        if ( month==1 .or. month==3 .or. month==5 .or. month==7 .or.    &
           & month==8 .or. month==10 ) then
          if ( mday>31 ) then
            mday = mday - 31
            month = month + 1
          end if
        else if ( month==12 ) then
          if ( mday>31 ) then
            mday = mday - 31
            month = 1
            myear = myear + 1
          end if
        else if ( month==4 .or. month==6 .or. month==9 .or. month==11 ) &
                & then
          if ( mday>30 ) then
            mday = mday - 30
            month = month + 1
          end if
        else if ( mod(myear,4)/=0 ) then
          if ( mday>28 ) then
            mday = mday - 28
            month = month + 1
          end if
        else if ( mod(myear,400)==0 ) then
          if ( mday>29 ) then
            mday = mday - 29
            month = month + 1
          end if
        else if ( mod(myear,100)==0 ) then
          if ( mday>28 ) then
            mday = mday - 28
            month = month + 1
          end if
        else if ( mday>29 ) then
          mday = mday - 29
          month = month + 1
        else
        end if
        wkday(i) = myear*10000 + month*100 + mday
      end do
!
      end subroutine headwk
!
!-----------------------------------------------------------------------
!
      subroutine sst_mn(idate,idate0,ilon,jlat,sst,pathaddname)
      use netcdf
      implicit none
!
! Dummy arguments
!
      integer :: idate , idate0 , ilon , jlat
      character(256) :: pathaddname
      intent (in) idate , idate0 , ilon , jlat , pathaddname
      real(4) , dimension(ilon,jlat) :: sst
      intent (out) :: sst
!
! Local variables
!
      integer :: i , it , j , month , n , nday , nhour , nyear
      logical :: there
      character(5) :: varname
      integer(2) , dimension(ilon,jlat) :: work
      integer :: istatus
!
      integer , dimension(10) , save :: icount , istart
      integer , save :: inet , ivar
      real(8) , save :: xadd , xscale
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
      if ( idate==idate0 ) then
        inquire (file=pathaddname,exist=there)
        if ( .not.there ) then
          write (*,*) trim(pathaddname) , ' is not available'
          stop
        end if
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
        istart(1) = 1
        istart(2) = 1
        icount(1) = ilon
        icount(2) = jlat
        do n = 4 , 10
          istart(n) = 0
          icount(n) = 0
        end do
      end if
 
!bxq
      nyear = idate/1000000
      month = idate/10000 - nyear*100
      nday = idate/100 - nyear*10000 - month*100
      nhour = idate - nyear*1000000 - month*10000 - nday*100
 
      it = (nyear-1981)*12 + month - 11
 
      istart(3) = it
      icount(3) = 1
      istatus = nf90_get_var(inet,ivar,work,istart,icount)
      if ( istatus/=nf90_noerr ) then
        write ( 6,* ) 'Cannot get ', varname, ' from file'
        write ( 6,* ) istart
        write ( 6,* ) icount
        write ( 6,* ) nf90_strerror(istatus)
        stop 'ERROR READ SST'
      end if
!bxq_
!
      do j = 1 , jlat
        do i = 1 , ilon
          if ( work(i,j)==32767 ) then
             sst(i,jlat+1-j) = -9999.
          else
             sst(i,jlat+1-j) = work(i,j)*xscale + xadd
          end if
        end do
      end do
!
      end subroutine sst_mn
!
!-----------------------------------------------------------------------
!
      subroutine ice_mn(idate,idate0,ilon,jlat,ice,pathaddname)
      use netcdf
      implicit none
!
! Dummy arguments
!
      integer :: idate , idate0 , ilon , jlat
      character(256) :: pathaddname
      intent (in) idate , idate0 , ilon , jlat , pathaddname
      real(4) , dimension(ilon,jlat) :: ice
      intent (out) :: ice
!
! Local variables
!
      integer :: i , it , j , month , n , nday , nhour , nyear
      logical :: there
      character(5) :: varname
      integer(2) , dimension(ilon,jlat) :: work
      integer :: istatus
!
      integer , dimension(10) , save :: icount , istart
      integer , save :: inet , ivar
      real(8) , save :: xadd , xscale
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
      data varname/'icec'/
!
      if ( idate==idate0 ) then
        inquire (file=pathaddname,exist=there)
        if ( .not.there ) then
          write (*,*) trim(pathaddname) , ' is not available'
          stop
        end if
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
        istart(1) = 1
        istart(2) = 1
        icount(1) = ilon
        icount(2) = jlat
        do n = 4 , 10
          istart(n) = 0
          icount(n) = 0
        end do
      end if
 
!bxq
      nyear = idate/1000000
      month = idate/10000 - nyear*100
      nday = idate/100 - nyear*10000 - month*100
      nhour = idate - nyear*1000000 - month*10000 - nday*100
 
      it = (nyear-1981)*12 + month - 11
 
      istart(3) = it
      icount(3) = 1
      istatus = nf90_get_var(inet,ivar,work,istart,icount)
      if ( istatus/=nf90_noerr ) then
        write ( 6,* ) 'Cannot get ', varname, ' from file'
        write ( 6,* ) istart
        write ( 6,* ) icount
        write ( 6,* ) nf90_strerror(istatus)
        stop 'ERROR READ SST'
      end if
!bxq_
!
      do j = 1 , jlat
        do i = 1 , ilon
          if ( work(i,j)==32767 ) then
             ice(i,jlat+1-j) = -9999.
          else
             ice(i,jlat+1-j) = work(i,j)*xscale + xadd
          end if
        end do
      end do
!
      end subroutine ice_mn
!
!-----------------------------------------------------------------------
!
      subroutine sst_wk(idate,idate0,kkk,ilon,jlat,sst,pathaddname)
      use netcdf
      implicit none
!
! Dummy arguments
!
      integer :: idate , idate0 , kkk , ilon , jlat
      character(256) :: pathaddname
      intent (in) idate , idate0 , kkk , ilon , jlat , pathaddname
      real(4) , dimension(ilon,jlat) :: sst
      intent (out) :: sst
!
! Local variables
!
      integer :: i , it , j , month , n , nday , nhour , nyear
      logical :: there
      character(3) :: varname
      integer :: istatus
      integer(2) , dimension(ilon,jlat) :: work
!
      integer , dimension(10) , save :: icount , istart
      integer , save :: inet , ivar
      real(8) , save :: xadd , xscale
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
      if ( idate==idate0 ) then
        inquire (file=pathaddname,exist=there)
        if ( .not.there ) then
          write (*,*) trim(pathaddname) , ' is not available'
          stop
        end if
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
        istart(1) = 1
        istart(2) = 1
        icount(1) = ilon
        icount(2) = jlat
        do n = 4 , 10
          istart(n) = 0
          icount(n) = 0
        end do
      end if
      if ( idate0<1989123100 .and. idate==1989123100 ) then
        istatus = nf90_close(inet)
        inquire (file=pathaddname,exist=there)
        if ( .not.there ) then
          write (*,*) trim(pathaddname) , ' is not available'
          stop
        end if
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
        istart(1) = 1
        istart(2) = 1
        icount(1) = ilon
        icount(2) = jlat
        do n = 4 , 10
          istart(n) = 0
          icount(n) = 0
        end do
      end if
!bxq
      nyear = idate/1000000
      month = idate/10000 - nyear*100
      nday = idate/100 - nyear*10000 - month*100
      nhour = idate - nyear*1000000 - month*10000 - nday*100
 
      it = kkk
      if ( kkk>427 ) it = kkk - 427
 
      istart(3) = it
      icount(3) = 1
      istatus = nf90_get_var(inet,ivar,work,istart,icount)
      if ( istatus/=nf90_noerr ) then
        write ( 6,* ) 'Cannot get ', varname, ' from file'
        write ( 6,* ) istart
        write ( 6,* ) icount
        write ( 6,* ) nf90_strerror(istatus)
        stop 'ERROR READ SST'
      end if
!bxq_
      do j = 1 , jlat
        do i = 1 , ilon
          if ( work(i,j)==32767 ) then
             sst(i,jlat+1-j) = -9999.
          else
             sst(i,jlat+1-j) = work(i,j)*xscale + xadd
          end if
        end do
      end do

      end subroutine sst_wk
!
!-----------------------------------------------------------------------
!
      subroutine ice_wk(idate,idate0,kkk,ilon,jlat,ice,pathaddname)
      use netcdf
      implicit none
!
! Dummy arguments
!
      integer :: idate , idate0 , kkk , ilon , jlat
      character(256) :: pathaddname
      intent (in) idate , idate0 , kkk , ilon , jlat , pathaddname
      real(4) , dimension(ilon,jlat) :: ice
      intent (out) :: ice
!
! Local variables
!
      integer :: i , it , j , month , n , nday , nhour , nyear
      logical :: there
      character(4) :: varname
      integer(2) , dimension(ilon,jlat) :: work
      integer :: istatus
!
      integer , dimension(10) , save :: icount , istart
      integer , save :: inet , ivar
      real(8) , save :: xadd , xscale
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
      data varname/'icec'/
!
      if ( idate==idate0 ) then
        inquire (file=pathaddname,exist=there)
        if ( .not.there ) then
          write (*,*) trim(pathaddname) , ' is not available'
          stop
        end if
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
        istart(1) = 1
        istart(2) = 1
        icount(1) = ilon
        icount(2) = jlat
        do n = 4 , 10
          istart(n) = 0
          icount(n) = 0
        end do
      end if
      if ( idate0<1989123100 .and. idate==1989123100 ) then
        istatus = nf90_close(inet)
        inquire (file=pathaddname,exist=there)
        if ( .not.there ) then
          write (*,*) trim(pathaddname) , ' is not available'
          stop
        end if
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
        istart(1) = 1
        istart(2) = 1
        icount(1) = ilon
        icount(2) = jlat
        do n = 4 , 10
          istart(n) = 0
          icount(n) = 0
        end do
      end if
!bxq
      nyear = idate/1000000
      month = idate/10000 - nyear*100
      nday = idate/100 - nyear*10000 - month*100
      nhour = idate - nyear*1000000 - month*10000 - nday*100
 
      it = kkk
      if ( kkk>427 ) it = kkk - 427
 
      istart(3) = it
      icount(3) = 1
      istatus = nf90_get_var(inet,ivar,work,istart,icount)
      if ( istatus/=nf90_noerr ) then
        write ( 6,* ) 'Cannot get ', varname, ' from file'
        write ( 6,* ) istart
        write ( 6,* ) icount
        write ( 6,* ) nf90_strerror(istatus)
        stop 'ERROR READ SST'
      end if
!bxq_
      do j = 1 , jlat
        do i = 1 , ilon
          if ( work(i,j)==32767 ) then
             ice(i,jlat+1-j) = -9999.
          else
             ice(i,jlat+1-j) = work(i,j)*xscale + xadd
          end if
        end do
      end do

      end subroutine ice_wk
!
      end module mod_sst_1deg

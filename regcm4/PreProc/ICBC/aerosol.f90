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

      program aerosol

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Comments on dataset sources and location:                          c
!                                                                    c
! EDGAR                                                              c
!                                                                    c
! LIOUSSE96                                                          c
!                                                                    c
! BOND                                                               c
!                                                                    c
! GEIA                                                               c
!                                                                    c
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      use mod_dynparam
      use netcdf

      implicit none
!
! PARAMETER definitions
!
      integer , parameter :: ilon = 360 , jlat = 180
!
! Local variables
!
      logical :: there
      integer :: i , j , nrec , ierr
      integer :: iyy , jxx
      integer :: istatus , incin , idimid , ivarid
      real(4) , dimension(jlat) :: lati
      real(4) , dimension(ilon) :: loni
      real(4) , dimension(ilon,jlat) :: aer2
      real(4) , allocatable , dimension(:,:) :: aermm , xlat , xlon
      real(4) , allocatable , dimension(:,:) :: finmat
      character(256) :: namelistfile, prgname , terfile , aerofile
!
!     Read input global namelist
!
      call getarg(0, prgname)
      call getarg(1, namelistfile)
      call initparam(namelistfile, ierr)
      if ( ierr/=0 ) then
        write ( 6, * ) 'Parameter initialization not completed'
        write ( 6, * ) 'Usage : '
        write ( 6, * ) '          ', trim(prgname), ' regcm.in'
        write ( 6, * ) ' '
        write ( 6, * ) 'Check argument and namelist syntax'
        stop
      end if
!
      allocate(aermm(iy,jx))
      allocate(xlat(iy,jx))
      allocate(xlon(iy,jx))
      allocate(finmat(jx,iy))

      inquire (file=trim(inpglob)//'/AERGLOB/AEROSOL.dat',exist=there)
      if ( .not.there ) print * , 'AEROSOL.dat is not available' ,      &
                             &' under ',trim(inpglob),'/AERGLOB/'
      open (11,file=trim(inpglob)//'/AERGLOB/AEROSOL.dat',              &
          & form='unformatted',recl=ilon*jlat*ibyte,access='direct',    &
          & status='old',err=100)

      write (aerofile,99001)                                            &
        & trim(dirglob), pthsep, trim(domname), '_AERO.dat'
      open (25,file=aerofile,form='unformatted',                        &
          & recl=iy*jx*ibyte,access='direct',status='replace')
      if ( igrads==1 ) then
        write (aerofile,99001)                                          &
          & trim(dirglob), pthsep, trim(domname), '_AERO.ctl'
        open (31,file=aerofile,status='replace')
        write (31,'(a,a,a)') 'dset ^',trim(domname),'_AERO.dat'
      end if
 
      terfile = trim(dirter)//pthsep//trim(domname)//'_DOMAIN000.nc'
      istatus = nf90_open(terfile, nf90_nowrite, incin)
      if ( istatus /= nf90_noerr) then
        write (6,*) 'Error Opening Domain file ', trim(terfile)
        write (6,*) nf90_strerror(istatus)
        stop
      end if

      istatus = nf90_inq_dimid(incin, "iy", idimid)
      if (istatus /= nf90_noerr) then
        write (6,*) 'Dimension iy missing'
        write (6,*) nf90_strerror(istatus)
        stop
      end if
      istatus = nf90_inquire_dimension(incin, idimid, len=iyy)
      if (istatus /= nf90_noerr) then
        write (6,*) 'Error dimension iy'
        write (6,*) nf90_strerror(istatus)
        stop
      end if
      istatus = nf90_inq_dimid(incin, "jx", idimid)
      if (istatus /= nf90_noerr) then
        write (6,*) 'Dimension jx missing'
        write (6,*) nf90_strerror(istatus)
        stop
      end if
      istatus = nf90_inquire_dimension(incin, idimid, len=jxx)
      if (istatus /= nf90_noerr) then
        write (6,*) 'Error dimension jx'
        write (6,*) nf90_strerror(istatus)
        stop
      end if
      if ( iyy/=iy .or. jxx/=jx ) then
        print * , 'IMPROPER DIMENSION SPECIFICATION'
        print * , '  namelist   : ' , iy , jx
        print * , '  DOMAIN     : ' , iyy , jxx
        stop 'Dimensions mismatch'
      end if
      istatus = nf90_inq_varid(incin, "xlat", ivarid)
      if (istatus /= nf90_noerr) then
        write (6,*) 'Error xlat variable undefined'
        write (6,*) nf90_strerror(istatus)
        stop
       end if
      istatus = nf90_get_var(incin, ivarid, finmat)
      if (istatus /= nf90_noerr) then
        write (6,*) 'Error reading xlat variable'
        write (6,*) nf90_strerror(istatus)
        stop
      end if
      xlat = transpose(finmat)
      istatus = nf90_inq_varid(incin, "xlon", ivarid)
      if (istatus /= nf90_noerr) then
        write (6,*) 'Error xlon variable undefined'
        write (6,*) nf90_strerror(istatus)
        stop
      end if
      istatus = nf90_get_var(incin, ivarid, finmat)
      if (istatus /= nf90_noerr) then
        write (6,*) 'Error reading xlon variable'
        write (6,*) nf90_strerror(istatus)
        stop
      end if
      xlon = transpose(finmat)
      istatus = nf90_close(incin)
      if (istatus /= nf90_noerr) then
        write (6,*) 'Error closing Domain file ', trim(terfile)
        write (6,*) nf90_strerror(istatus)
        stop
      end if
!
      deallocate(finmat)
!
      call gridml(xlon,xlat,iy,jx,ibyte)
!
 
!     ******    SET UP LONGITUDES AND LATITUDES FOR AEROSOL DATA
      do i = 1 , ilon
        loni(i) = -179.5 + float(i-1)
      end do
      do j = 1 , jlat
        lati(j) = -89.5 + 1.*float(j-1)
      end do
 
!     ****** ALL AEROSOL DATA, 1 Deg data, Climate value
      do nrec = 1 , 39
        read (11,rec=nrec) aer2
 
        call bilinx(aer2,aermm,xlon,xlat,loni,lati,ilon,jlat,iy,jx,1)
 
!       ******           WRITE OUT AEROSOL DATA ON RegCM GRID
        write (25,rec=nrec) ((aermm(i,j),j=1,jx),i=1,iy)
      end do
 
      deallocate(aermm)
      deallocate(xlat)
      deallocate(xlon)

      print *, 'Successfully built aerosol data for domain'
      stop

 100  continue
      print * , 'ERROR OPENING AEROSOL FILE'
      stop '4810 IN PROGRAM AEROSOL'
99001 format (a,a,a,a)
      end program aerosol
!
!-----------------------------------------------------------------------
!
      subroutine gridml(xlon,xlat)
      use mod_dynparam
      implicit none
!
! Dummy arguments
!
      real(4) , dimension(iy,jx) :: xlat , xlon
      intent (in) xlat , xlon
!
! Local variables
!
      real(4) :: alatmax , alatmin , alonmax , alonmin , centeri ,      &
            & centerj , dsinm , rlatinc , rloninc
      character(3) , dimension(12) :: cmonth
      integer :: i , j , month , nx , ny , period
!
      data cmonth/'jan' , 'feb' , 'mar' , 'apr' , 'may' , 'jun' ,       &
         & 'jul' , 'aug' , 'sep' , 'oct' , 'nov' , 'dec'/
!
      alatmin = 999999.
      alatmax = -999999.
      alonmin = 999999.
      alonmax = -999999.
      nx = 0
      ny = 0
      dsinm = ds* 1000.0
!
      if ( igrads==1 ) then
        write (31,'(a)')                                                &
                     &'title AEROSOL fields for RegCM domain, kg/m2/sec'
        if ( ibigend==1 ) then
          write (31,'(a)') 'options big_endian'
        else
          write (31,'(a)') 'options little_endian'
        end if
        write (31,'(a)') 'undef -9999.'
        if ( iproj=='LAMCON' .or. iproj=='ROTMER' ) then
          do j = 1 , jx
            if ( xlat(1,j)<alatmin ) alatmin = xlat(1,j)
            if ( xlat(iy,j)>alatmax ) alatmax = xlat(iy,j)
          end do
          do i = 1 , iy
            do j = 1 , jx
              if ( clon>=0.0 ) then
                if ( xlon(i,j)>=0.0 ) then
                  alonmin = amin1(alonmin,xlon(i,j))
                  alonmax = amax1(alonmax,xlon(i,j))
                else if ( abs(clon-xlon(i,j))<abs(clon-(xlon(i,j)+360.))&
                        & ) then
                  alonmin = amin1(alonmin,xlon(i,j))
                  alonmax = amax1(alonmax,xlon(i,j))
                else
                  alonmin = amin1(alonmin,xlon(i,j)+360.)
                  alonmax = amax1(alonmax,xlon(i,j)+360.)
                end if
              else if ( xlon(i,j)<0.0 ) then
                alonmin = amin1(alonmin,xlon(i,j))
                alonmax = amax1(alonmax,xlon(i,j))
              else if ( abs(clon-xlon(i,j))<abs(clon-(xlon(i,j)-360.)) )&
                      & then
                alonmin = amin1(alonmin,xlon(i,j))
                alonmax = amax1(alonmax,xlon(i,j))
              else
                alonmin = amin1(alonmin,xlon(i,j)-360.)
                alonmax = amax1(alonmax,xlon(i,j)-360.)
              end if
            end do
          end do
          rlatinc = ds/111./2.
          rloninc = ds/111./2.
          ny = 2 + nint(abs(alatmax-alatmin)/rlatinc)
          nx = 1 + nint(abs((alonmax-alonmin)/rloninc))
 
          centerj = jx/2.
          centeri = iy/2.
        end if
        if ( iproj=='LAMCON' ) then        ! Lambert projection
          write (31,99001) jx , iy , clat , clon , centerj , centeri ,  &
                         & truelatl , truelath , clon , dsinm , dsinm
          write (31,99002) nx + 2 , alonmin - rloninc , rloninc
          write (31,99003) ny + 2 , alatmin - rlatinc , rlatinc
        else if ( iproj=='POLSTR' ) then   !
        else if ( iproj=='NORMER' ) then
          write (31,99004) jx , xlon(1,1) , xlon(1,2) - xlon(1,1)
          write (31,99005) iy
          write (31,99006) (xlat(i,1),i=1,iy)
        else if ( iproj=='ROTMER' ) then
          write (*,*) 'Note that rotated Mercartor (ROTMER)' ,          &
                     &' projections are not supported by GrADS.'
          write (*,*) '  Although not exact, the eta.u projection' ,    &
                     &' in GrADS is somewhat similar.'
          write (*,*) ' FERRET, however, does support this projection.'
          write (31,99007) jx , iy , plon , plat , ds/111. ,            &
                         & ds/111.*.95238
          write (31,99002) nx + 2 , alonmin - rloninc , rloninc
          write (31,99003) ny + 2 , alatmin - rlatinc , rlatinc
        else
          write (*,*) 'Are you sure your map projection is correct ?'
          stop
        end if
        write (31,99008) 1 , 1000.
        month = 1
        period = 1
        write (31,99009) period , cmonth(month) , 2001
        write (31,99010) 39
        write (31,99011) 'so2   ' ,                                     &
                        &'Anthropogenic SO2 emission, EDGAR       '
        write (31,99011) 'bc    ' ,                                     &
                        &'Anthropogenic Black Carbon (BC), EDGAR  '
        write (31,99011) 'oc    ' ,                                     &
                        &'Anthropogenic Organic Carbon (OC), EDGAR'
        write (31,99011) 'so201 ' ,                                     &
                        &'Biomass SO2 emission, EDGAR, January    '
        write (31,99011) 'so202 ' ,                                     &
                        &'Biomass SO2 emission, EDGAR, February   '
        write (31,99011) 'so203 ' ,                                     &
                        &'Biomass SO2 emission, EDGAR, March      '
        write (31,99011) 'so204 ' ,                                     &
                        &'Biomass SO2 emission, EDGAR, April      '
        write (31,99011) 'so205 ' ,                                     &
                        &'Biomass SO2 emission, EDGAR, May        '
        write (31,99011) 'so206 ' ,                                     &
                        &'Biomass SO2 emission, EDGAR, June       '
        write (31,99011) 'so207 ' ,                                     &
                        &'Biomass SO2 emission, EDGAR, July       '
        write (31,99011) 'so208 ' ,                                     &
                        &'Biomass SO2 emission, EDGAR, August     '
        write (31,99011) 'so209 ' ,                                     &
                        &'Biomass SO2 emission, EDGAR, September  '
        write (31,99011) 'so210 ' ,                                     &
                        &'Biomass SO2 emission, EDGAR, October    '
        write (31,99011) 'so211 ' ,                                     &
                        &'Biomass SO2 emission, EDGAR, November   '
        write (31,99011) 'so212 ' ,                                     &
                        &'Biomass SO2 emission, EDGAR, December   '
        write (31,99011) 'bc_01 ' ,                                     &
                        &'Biomass BC emission, LIOUSSE, January   '
        write (31,99011) 'bc_02 ' ,                                     &
                        &'Biomass BC emission, LIOUSSE, February  '
        write (31,99011) 'bc_03 ' ,                                     &
                        &'Biomass BC emission, LIOUSSE, March     '
        write (31,99011) 'bc_04 ' ,                                     &
                        &'Biomass BC emission, LIOUSSE, April     '
        write (31,99011) 'bc_05 ' ,                                     &
                        &'Biomass BC emission, LIOUSSE, May       '
        write (31,99011) 'bc_06 ' ,                                     &
                        &'Biomass BC emission, LIOUSSE, June      '
        write (31,99011) 'bc_07 ' ,                                     &
                        &'Biomass BC emission, LIOUSSE, July      '
        write (31,99011) 'bc_08 ' ,                                     &
                        &'Biomass BC emission, LIOUSSE, August    '
        write (31,99011) 'bc_09 ' ,                                     &
                        &'Biomass BC emission, LIOUSSE, September '
        write (31,99011) 'bc_10 ' ,                                     &
                        &'Biomass BC emission, LIOUSSE, October   '
        write (31,99011) 'bc_11 ' ,                                     &
                        &'Biomass BC emission, LIOUSSE, November  '
        write (31,99011) 'bc_12 ' ,                                     &
                        &'Biomass BC emission, LIOUSSE, December  '
        write (31,99011) 'oc_01 ' ,                                     &
                        &'Biomass OC emission, LIOUSSE, January   '
        write (31,99011) 'oc_02 ' ,                                     &
                        &'Biomass OC emission, LIOUSSE, February  '
        write (31,99011) 'oc_03 ' ,                                     &
                        &'Biomass OC emission, LIOUSSE, March     '
        write (31,99011) 'oc_04 ' ,                                     &
                        &'Biomass OC emission, LIOUSSE, April     '
        write (31,99011) 'oc_05 ' ,                                     &
                        &'Biomass OC emission, LIOUSSE, May       '
        write (31,99011) 'oc_06 ' ,                                     &
                        &'Biomass OC emission, LIOUSSE, June      '
        write (31,99011) 'oc_07 ' ,                                     &
                        &'Biomass OC emission, LIOUSSE, July      '
        write (31,99011) 'oc_08 ' ,                                     &
                        &'Biomass OC emission, LIOUSSE, August    '
        write (31,99011) 'oc_09 ' ,                                     &
                        &'Biomass OC emission, LIOUSSE, September '
        write (31,99011) 'oc_10 ' ,                                     &
                        &'Biomass OC emission, LIOUSSE, October   '
        write (31,99011) 'oc_11 ' ,                                     &
                        &'Biomass OC emission, LIOUSSE, November  '
        write (31,99011) 'oc_12 ' ,                                     &
                        &'Biomass OC emission, LIOUSSE, December  '
 
        write (31,'(a)') 'endvars'
        close (31)
      end if
99001 format ('pdef ',i4,1x,i4,1x,'lcc',7(1x,f7.2),1x,2(f7.0,1x))
99002 format ('xdef ',i4,' linear ',f7.2,1x,f7.4)
99003 format ('ydef ',i4,' linear ',f7.2,1x,f7.4)
99004 format ('xdef ',i3,' linear ',f9.4,' ',f9.4)
99005 format ('ydef ',i3,' levels')
99006 format (10F7.2)
99007 format ('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
99008 format ('zdef ',i1,' levels ',f7.2)
99009 format ('tdef ',i4,' linear 00z16',a3,i4,' 1mo')
99010 format ('vars ',i2)
99011 format (a6,'0 99 ',a40)
!
      end subroutine gridml
!
!
!
      subroutine bilinx(fin,fout,lono,lato,loni,lati,nloni,nlati,iy,jx, &
                      & nflds)
      implicit none
!
! Dummy arguments
!
      integer :: iy , jx , nflds , nlati , nloni
      real(4) , dimension(nloni,nlati,nflds) :: fin
      real(4) , dimension(nlati) :: lati
      real(4) , dimension(iy,jx) :: lato , lono
      real(4) , dimension(nloni) :: loni
      real(4) , dimension(iy,jx,nflds) :: fout
      intent (in) fin , iy , jx , lati , lato , loni , lono , nflds ,   &
                & nlati , nloni
      intent (out) fout
!
! Local variables
!
      real(4) :: bas , lon180 , p , q , xsum , xind , yind
      integer :: i , ip , ipp1 , j , jq , jqp1 , l
!
!
!     PERFORMING BI-LINEAR INTERPOLATION USING 4 GRID POINTS FROM A
!     BIGGER RECTANGULAR GRID TO A GRID DESCRIBED BY XLONS AND XLATS OF
!     GRID2. A POINT ON GRID2 IS TRAPPED WITHIN FOUR GRID POINTS ON
!     GRID4.THE GRID POINTS ARE ALWAYS TO THE NORTH AND EAST OF THE
!     TRAPPED POINT.. THE ALGORITHM COMPUTES THE FRACTIONAL DISTANCES
!     IN BOTH X AND Y DIRECTION OF THE TRAPPED GRID POINT AND USES THE
!     INFORMATION AS WEIGHTING FACTORS IN THE INTERPOLATION.
!     THERE IS ONE LESS ROW AND COLUMN WHEN THE SCALAR FIELDS ARE
!     INTERPOLATED BECAUSE XLATS AND XLONS ARE NOT DEFINED FOR
!     THE CROSS POINTS IN THE MM4 MODEL.
!
!     IN(NLONI,NLATI,NFLDS)  IS THE INPUT FIELD ON REGULAR LAT/LON GRID.
!     OUT(NLATO,NLONO,NFLDS) IS THE OUTPUT FIELD ON LAMBERT CONFORMAL
!     GRID. LONI.....LONGITUDE VALUES IN DEGREES OF THE LAT-LON GRID.
!     LATI.....LATITUDE VALUES IN DEGREES OF THE LAT-LON GRID.
!     P.........EAST-WEST WEIGHTING FACTOR.
!     Q.........NORTH-SOUTH WEIGHTING FACTOR.
!     IP........GRID POINT LOCATION IN EAST-WEST OF TRAPPED GRID POINT.
!     IQ........GRID POINT LOCATION IN NORTH-SOUTH OF TRAPPED GRID
 
!     POINT.
 
      do j = 1 , jx
        do i = 1 , iy
 
          yind = (((lato(i,j)-lati(1))/(lati(nlati)-lati(1)))           &
               & *float(nlati-1)) + 1.
          jq = int(yind)
          jq = max0(jq,1)
          jqp1 = min0(jq+1,nlati)
          q = yind - jq
 
          lon180 = lono(i,j)
          if ( lono(i,j)<-180. ) lon180 = lono(i,j) + 360.
          if ( lono(i,j)>180. ) lon180 = lono(i,j) - 360.
          xind = (((lon180-loni(1))/(loni(nloni)-loni(1)))              &
               & *float(nloni-1)) + 1.
          ip = int(xind)
          ip = max0(ip,1)
          ipp1 = min0(ip+1,nloni)
          p = xind - ip
 
          do l = 1 , nflds
            xsum = 0.0
            bas = 0.0
            if ( fin(ip,jq,l)<-9990.0 .and. fin(ipp1,jq,l)<-9990.0 .and.&
               & fin(ipp1,jqp1,l)<-9990.0 .and. fin(ip,jqp1,l)<-9990.0 )&
               & then
              fout(i,j,l) = -9999.
            else
              if ( fin(ip,jq,l)>-9990.0 ) then
                xsum = xsum + (1.-q)*(1.-p)*fin(ip,jq,l)
                bas = bas + (1.-q)*(1.-p)
              end if
              if ( fin(ipp1,jq,l)>-9990.0 ) then
                xsum = xsum + (1.-q)*p*fin(ipp1,jq,l)
                bas = bas + (1.-q)*p
              end if
              if ( fin(ipp1,jqp1,l)>-9990.0 ) then
                xsum = xsum + q*p*fin(ipp1,jqp1,l)
                bas = bas + q*p
              end if
              if ( fin(ip,jqp1,l)>-9990.0 ) then
                xsum = xsum + q*(1.-p)*fin(ip,jqp1,l)
                bas = bas + q*(1.-p)
              end if
              fout(i,j,l) = xsum/bas
            end if
          end do
        end do
 
      end do
 
      end subroutine bilinx

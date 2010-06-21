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

      module mod_sst_grid
      use netcdf
      use mod_dynparam
      implicit none

      real(4) , allocatable , dimension(:,:) :: lu , sstmm , icemm ,    &
                                  &             xlat , xlon , finmat

      contains

      subroutine init_grid
        implicit none
        allocate(lu(iy,jx))
        allocate(sstmm(iy,jx))
        allocate(icemm(iy,jx))
        allocate(xlat(iy,jx))
        allocate(xlon(iy,jx))
        allocate(finmat(jx,iy))
      end subroutine init_grid

      subroutine free_grid
        implicit none
        deallocate(lu)
        deallocate(sstmm)
        deallocate(icemm)
        deallocate(xlat)
        deallocate(xlon)
        deallocate(finmat)
      end subroutine free_grid

      subroutine read_domain(terfile)
        implicit none
        character(256) :: terfile
        intent(in) :: terfile

        integer :: iyy , jxx
        integer :: istatus , incin , idimid , ivarid

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

        istatus = nf90_inq_varid(incin, "landuse", ivarid)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error landuse variable undefined'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
         istatus = nf90_get_var(incin, ivarid, finmat)
        if (istatus /= nf90_noerr) then
          write (6,*) 'Error reading landuse variable'
          write (6,*) nf90_strerror(istatus)
          stop
        end if
        lu = transpose(finmat)
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

      end subroutine read_domain

      subroutine setup_sstfile(idate1,nsteps)

      implicit none
!
! Dummy arguments
!
      integer :: idate1 , nsteps
      intent (in) idate1 , nsteps
!
! Local variables
!
      real(4) :: alatmax , alatmin , alonmax , alonmin , centeri ,      &
            & centerj , rlatinc , rloninc , dsinm
      character(3) , dimension(12) :: cmonth
      integer :: i , j , nx , ny
      integer :: nyear , nmonth , nday , nhour
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
      dsinm = ds * 1000.0
      nyear = idate1/1000000
      nmonth = (idate1-nyear*1000000)/10000
      nday = (idate1-nyear*1000000-nmonth*10000)/100
      nhour = mod(idate1,100)
!
      if ( igrads==1 ) then
        write (31,'(a)') 'title SST fields for RegCM domain'
        if ( ibigend==1 ) then
          write (31,'(a)') 'options big_endian'
        else
          write (31,'(a)') 'options little_endian'
        end if
        write (31,'(a)') 'undef -9999.'
        if ( ssttyp=='OI2ST' .or. ssttyp=='OI2WK' ) then
           write (32,'(a)') 'title SeaIce fields for RegCM domain'
           if ( ibigend.eq.1 ) then
             write (32,'(a)') 'options big_endian'
           else
             write (32,'(a)') 'options little_endian'
           end if
           write (32,'(a)') 'undef -9999.'
        end if
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
          if ( ssttyp=='OI2ST' .or. ssttyp=='OI2WK' ) then
            write (32,99001) jx , iy , clat , clon , centerj , centeri ,&
                        & truelatl , truelath , clon , dsinm, dsinm
            write (32,99002) nx + 2 , alonmin - rloninc , rloninc
            write (32,99003) ny + 2 , alatmin - rlatinc , rlatinc
          end if
        else if ( iproj=='POLSTR' ) then   !
        else if ( iproj=='NORMER' ) then
          write (31,99004) jx , xlon(1,1) , xlon(1,2) - xlon(1,1)
          write (31,99005) iy
          write (31,99006) (xlat(i,1),i=1,iy)
          if ( ssttyp=='OI2ST' .or. ssttyp=='OI2WK' ) then
            write (32,99004) jx , xlon(1,1) , xlon(1,2) - xlon(1,1)
            write (32,99005) iy
            write (32,99006) (xlat(i,1),i=1,iy)
          end if 
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
          if ( ssttyp=='OI2ST' .or. ssttyp=='OI2WK' ) then
            write (32,99007) jx , iy , plon , plat , ds/111. ,          &
                           & ds/111.*.95238
            write (32,99002) nx + 2 , alonmin - rloninc , rloninc
            write (32,99003) ny + 2 , alatmin - rlatinc , rlatinc
          end if
        else
          write (*,*) 'Are you sure your map projection is correct ?'
          stop
        end if
        if ( ssttyp=='OI_ST' .or. ssttyp=='OI_WK' .or.                  &
          &  ssttyp=='FV_RF' .or. ssttyp=='FV_A2' .or.                  &
          &  ssttyp=='FV_B2') then
          write (31,99008) 1 , 1000.
          write (31,99009) nsteps , cmonth(nmonth) , idate1/100
          write (31,99011) 1
          write (31,99012) 'sst ' , 'Sea Surface Temperature    '
          write (31,'(a)') 'endvars'
          close (31)
        else
          write (31,99008) 1 , 1000.
          write (31,99010) nsteps , nhour , nday , cmonth(nmonth) , nyear
          write (31,99011) 1
          write (31,99012) 'sst ' , 'Sea Surface Temperature    '
          write (31,'(a)') 'endvars'
          close (31)
        end if
        if ( ssttyp=='OI2ST' .or. ssttyp=='OI2WK' ) then
          write (32,99008) 1 , 1000.
          write (32,99009) nsteps , cmonth(nmonth) , idate1/100
          write (32,99011) 1
          write (32,99012) 'ice ' , 'Sea Ice fraction           '
          write (32,'(a)') 'endvars'
          close (32)
        end if
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
99010 format ('tdef ',i6,' linear ',i2,'z',i0.2,a3,i4,' 6hr')
99011 format ('vars ',i1)
99012 format (a4,'0 99 ',a26)
!
      end subroutine setup_sstfile

      subroutine setup_sst_ice_file(idate1,inumber)
      implicit none
!
! Dummy arguments
!
      integer :: idate1 , inumber
      intent (in) idate1 , inumber
!
! Local variables
!
      real(4) :: alatmax , alatmin , alonmax , alonmin , centeri ,      &
            & centerj , rlatinc , rloninc , dsinm
      character(2) , dimension(31) :: cday
      character(3) , dimension(12) :: cmonth
      integer :: day , i , j , month , nx , ny
!
      data cmonth/'jan' , 'feb' , 'mar' , 'apr' , 'may' , 'jun' ,       &
         & 'jul' , 'aug' , 'sep' , 'oct' , 'nov' , 'dec'/
      data cday/'01' , '02' , '03' , '04' , '05' , '06' , '07' , '08' , &
         & '09' , '10' , '11' , '12' , '13' , '14' , '15' , '16' ,      &
         & '17' , '18' , '19' , '20' , '21' , '22' , '23' , '24' ,      &
         & '25' , '26' , '27' , '28' , '29' , '30' , '31'/
!
      alatmin = 999999.
      alatmax = -999999.
      alonmin = 999999.
      alonmax = -999999.
      nx = 0
      ny = 0
      dsinm = ds*1000.0
!
      if ( igrads==1 ) then
        write (31,'(a)') 'title SST fields for RegCM domain'
        if ( ibigend==1 ) then
          write (31,'(a)') 'options big_endian'
        else
          write (31,'(a)') 'options little_endian'
        end if
        write (31,'(a)') 'undef -9999.'
        if ( ssttyp=='OI2ST' .or. ssttyp=='OI2WK' ) then
           write (32,'(a)') 'title SeaIce fields for RegCM domain'
           if ( ibigend.eq.1 ) then
             write (32,'(a)') 'options big_endian'
           else
             write (32,'(a)') 'options little_endian'
           end if
           write (32,'(a)') 'undef -9999.'
        end if
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
          if ( ssttyp=='OI2ST' .or. ssttyp=='OI2WK' ) then
            write (32,99001) jx , iy , clat , clon , centerj , centeri ,&
                  & truelatl , truelath , clon , dsinm , dsinm
            write (32,99002) nx + 2 , alonmin - rloninc , rloninc
            write (32,99003) ny + 2 , alatmin - rlatinc , rlatinc
          end if
        else if ( iproj=='POLSTR' ) then   !
        else if ( iproj=='NORMER' ) then
          write (31,99004) jx , xlon(1,1) , xlon(1,2) - xlon(1,1)
          write (31,99005) iy
          write (31,99006) (xlat(i,1),i=1,iy)
          if ( ssttyp=='OI2ST' .or. ssttyp=='OI2WK' ) then
            write (32,99004) jx , xlon(1,1) , xlon(1,2) - xlon(1,1)
            write (32,99005) iy
            write (32,99006) (xlat(i,1),i=1,iy)
          end if 
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
          if ( ssttyp=='OI2ST' .or. ssttyp=='OI2WK' ) then
            write (32,99007) jx , iy , plon , plat , ds/111. ,          &
                           & ds/111.*.95238
            write (32,99002) nx + 2 , alonmin - rloninc , rloninc
            write (32,99003) ny + 2 , alatmin - rlatinc , rlatinc
          end if
        else
          write (*,*) 'Are you sure your map projection is correct ?'
          stop
        end if
        write (31,99008) 1 , 1000.
        month = idate1/100 - (idate1/10000)*100
        day = mod(idate1,100)
        write (31,99009) inumber , cday(day) , cmonth(month) ,          &
                       & idate1/10000
        write (31,99010) 1
        write (31,99011) 'sst ' , 'surface elevation          '
        write (31,'(a)') 'endvars'
        close (31)
        if ( ssttyp=='OI2ST' .or. ssttyp=='OI2WK' ) then
          write (32,99008) 1 , 1000.
          write (32,99009) inumber , cday(day) , cmonth(month) ,        &
                       & idate1/10000
          write (32,99010) 1
          write (32,99011) 'ice ' , 'Sea Ice fraction           '
          write (32,'(a)') 'endvars'
          close (32)
        end if
      end if
99001 format ('pdef ',i4,1x,i4,1x,'lcc',7(1x,f7.2),1x,2(f7.0,1x))
99002 format ('xdef ',i4,' linear ',f7.2,1x,f7.4)
99003 format ('ydef ',i4,' linear ',f7.2,1x,f7.4)
99004 format ('xdef ',i3,' linear ',f9.4,' ',f9.4)
99005 format ('ydef ',i3,' levels')
99006 format (10F7.2)
99007 format ('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
99008 format ('zdef ',i1,' levels ',f7.2)
99009 format ('tdef ',i4,' linear 00z',a2,a3,i4,' 7dy')
99010 format ('vars ',i1)
99011 format (a4,'0 99 ',a26)
!
      end subroutine setup_sst_ice_file

      end module mod_sst_grid

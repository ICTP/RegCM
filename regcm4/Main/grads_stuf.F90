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
 
      subroutine gradsbat(ctlname)
      use mod_dynparam
      use mod_message , only : fatal
      use mod_date
      use mod_param1
      use mod_param2
      use mod_param3 , only : r8pt , a
!#ifdef INTEL
!  include 'ifport.f90'
!#endif
      use mod_main
#ifdef MPP1
      use mod_mppio
#endif
      implicit none
!
! Dummy arguments
!
      character(18) :: ctlname
!
! Local variables
!
      real(4) :: alatmax , alatmin , alonmax , alonmin , centeri ,      &
               & centerj , rlatinc , rloninc
      character(2) , dimension(31) :: cday
      character(3) , dimension(12) :: cmonth
      integer :: i , ifrq , j , jbend , mnend , month , myear , nbase , &
               & nday , nhour , nnumb , nx , ny
!
      data cday/'01' , '02' , '03' , '04' , '05' , '06' , '07' , '08' , &
          &'09' , '10' , '11' , '12' , '13' , '14' , '15' , '16' ,      &
         & '17' , '18' , '19' , '20' , '21' , '22' , '23' , '24' ,      &
         & '25' , '26' , '27' , '28' , '29' , '30' , '31'/
      data cmonth/'jan' , 'feb' , 'mar' , 'apr' , 'may' , 'jun' ,       &
         & 'jul' , 'aug' , 'sep' , 'oct' , 'nov' , 'dec'/
!
      alatmin = 999999.
      alatmax = -999999.
      alonmin = 999999.
      alonmax = -999999.
      nx = 0
      ny = 0
      centerj = (jxm2)/2.
      centeri = (iym2)/2.

      open (31,file=trim(dirout)//pthsep//ctlname,status='replace')
      write (31,99001) ctlname(1:14)
      write (31,99002)
      if ( ibigend.eq.1 ) then
        write (31,99003)
      else if ( ibigend.eq.0 ) then
        write (31,99004)
      else
      end if
      if ( iotyp.eq.2 ) write (31,'(a)') 'options sequential'
      write (31,99005)
      if ( iproj.eq.'LAMCON' .or. iproj.eq.'ROTMER' ) then
        do j = 2 , jx-1
#ifdef MPP1
          if ( xlat_io(2,j).lt.alatmin ) alatmin = xlat_io(2,j)
          if ( xlat_io(iy-1,j).gt.alatmax ) alatmax = xlat_io(iy-1,j)
#else
          if ( xlat(2,j).lt.alatmin ) alatmin = xlat(2,j)
          if ( xlat(iy-1,j).gt.alatmax ) alatmax = xlat(iy-1,j)
#endif
        end do
        do i = 2 , iy-1
          do j = 2 , jx-1
#ifdef MPP1
            if ( clon.ge.0.0 ) then
              if ( xlong_io(i,j).ge.0.0 ) then
                alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
                alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
              else if ( abs(clon-xlong_io(i,j))                         &
                      & .lt.abs(clon-(xlong_io(i,j)+360.)) ) then
                alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
                alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
              else
                alonmin = amin1(alonmin,sngl(xlong_io(i,j))+360.)
                alonmax = amax1(alonmax,sngl(xlong_io(i,j))+360.)
              end if
            else if ( xlong_io(i,j).lt.0.0 ) then
              alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
              alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
            else if ( abs(clon-xlong_io(i,j))                           &
                    & .lt.abs(clon-(xlong_io(i,j)-360.)) ) then
              alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
              alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
            else
              alonmin = amin1(alonmin,sngl(xlong_io(i,j))-360.)
              alonmax = amax1(alonmax,sngl(xlong_io(i,j))-360.)
            end if
#else
            if ( clon.ge.0.0 ) then
              if ( xlong(i,j).ge.0.0 ) then
                alonmin = amin1(alonmin,sngl(xlong(i,j)))
                alonmax = amax1(alonmax,sngl(xlong(i,j)))
              else if ( abs(clon-xlong(i,j))                            &
                      & .lt.abs(clon-(xlong(i,j)+360.)) ) then
                alonmin = amin1(alonmin,sngl(xlong(i,j)))
                alonmax = amax1(alonmax,sngl(xlong(i,j)))
              else
                alonmin = amin1(alonmin,sngl(xlong(i,j))+360.)
                alonmax = amax1(alonmax,sngl(xlong(i,j))+360.)
              end if
            else if ( xlong(i,j).lt.0.0 ) then
              alonmin = amin1(alonmin,sngl(xlong(i,j)))
              alonmax = amax1(alonmax,sngl(xlong(i,j)))
            else if ( abs(clon-xlong(i,j))                              &
                    & .lt.abs(clon-(xlong(i,j)-360.)) ) then
              alonmin = amin1(alonmin,sngl(xlong(i,j)))
              alonmax = amax1(alonmax,sngl(xlong(i,j)))
            else
              alonmin = amin1(alonmin,sngl(xlong(i,j))-360.)
              alonmax = amax1(alonmax,sngl(xlong(i,j))-360.)
            end if
#endif
          end do
        end do
        rlatinc = dx*0.001/111./2.
        rloninc = dx*0.001/111./2.
        ny = 2 + nint(abs(alatmax-alatmin)/rlatinc)
        nx = 1 + nint(abs((alonmax-alonmin)/rloninc))
      end if
      if ( iotyp.eq.1 ) then
        if ( iproj.eq.'LAMCON' ) then   ! Lambert projection
          write (31,99006) jxm2 , iym2 , clat , clon , centerj ,        &
                         & centeri , truelatl , truelath , clon , dx ,&
                         & dx
          write (31,99007) nx + 2 , alonmin - rloninc , rloninc
          write (31,99008) ny + 2 , alatmin - rlatinc , rlatinc
        else if ( iproj.eq.'POLSTR' ) then
                                        !
        else if ( iproj.eq.'NORMER' ) then
#ifdef MPP1
          write (31,99009) jxm2 , xlong_io(2,2) , xlong_io(2,3)         &
                         & - xlong_io(2,2)
          write (31,99010) iym2
          write (31,99011) (xlat_io(i,2),i=2,iym1)
#else
          write (31,99009) jxm2 , xlong(2,2) , xlong(2,3) - xlong(2,2)
          write (31,99010) iym2
          write (31,99011) (xlat(i,2),i=2,iym1)
#endif
        else if ( iproj.eq.'ROTMER' ) then
          write (*,*) 'Note that rotated Mercartor (ROTMER)' ,          &
                     &' projections are not supported by GrADS.'
          write (*,*) '  Although not exact, the eta.u projection' ,    &
                     &' in GrADS is somewhat similar.'
          write (*,*) ' FERRET, however, does support this projection.'
          write (31,99012) jxm2 , iym2 , plon , plat ,                  &
                         & dx/111000. , dx/111000.*.95238
          write (31,99007) nx + 2 , alonmin - rloninc , rloninc
          write (31,99008) ny + 2 , alatmin - rlatinc , rlatinc
        else
          call fatal(__FILE__,__LINE__,'INVALID MAP PROJECTION')
        end if
      else if ( iotyp.eq.2 ) then
#ifdef MPP1
        write (31,99009) jxm2 , xlong_io(2,2) , xlong_io(2,3)           &
                       & - xlong_io(2,2)
        write (31,99013) iym2 , xlat_io(2,2) , xlat_io(3,2)             &
                       & - xlat_io(2,2)
#else
        write (31,99009) jxm2 , xlong(2,2) , xlong(2,3) - xlong(2,2)
        write (31,99013) iym2 , xlat(2,2) , xlat(3,2) - xlat(2,2)
#endif
      else
      end if
      write (31,99014) (1013.25-r8pt*10.)*a(kz) + r8pt*10.
      myear = ldatez/1000000
      month = (ldatez-myear*1000000)/10000
      nday = (ldatez-myear*1000000-month*10000)/100
      nhour = mod(ldatez,100)
      call finddate(nbase,ldatez)
      if ( month.eq.12 ) then
        call finddate(mnend,myear*1000000+1010100)
      else
        call finddate(mnend,myear*1000000+month*10000+10100)
      end if
      call finddate(jbend,idate2)
      if ( ldatez.eq.idate0 ) then
        nnumb = (ibdyfrq/batfrq+0.00001)*(min0(jbend,mnend)-nbase) + 1
      else
        nnumb = (ibdyfrq/batfrq+0.00001)*(min0(jbend,mnend)-nbase)
      end if
      ifrq = batfrq + 0.00001
      if ( ldatez.eq.idate0 ) then
        write (31,99015) nnumb , nhour , cday(nday) , cmonth(month) ,   &
                       & myear , ifrq
      else
        write (31,99015) nnumb , nhour + ifrq , cday(nday) ,            &
                       & cmonth(month) , myear , ifrq
      end if
      if ( iotyp.eq.2 ) write (31,'(a)') 'theader 4'
      write (31,99016) 21 + 6
      if ( iproj.eq.'LAMCON' .and. iotyp.eq.1 ) then   ! Lambert projection
        write (31,99018) 'u10m    ' ,                                   &
                        &'westerly  wind at 10m (m/s)          '
        write (31,99019) 'v10m    ' ,                                   &
                        &'southerly wind at 10m (m/s)          '
      else
        write (31,99017) 'u10m    ' ,                                   &
                        &'westerly  wind at 10m (m/s)          '
        write (31,99017) 'v10m    ' ,                                   &
                        &'southerly wind at 10m (m/s)          '
      end if
      write (31,99017) 'uvdrag  ' ,                                     &
                      &'surface drag stress                  '
      write (31,99017) 'tg      ' ,                                     &
                      &'ground temperature (degree)          '
      write (31,99017) 'tlef    ' ,                                     &
                      &'temperature of foliage               '
      write (31,99017) 't2m     ' ,                                     &
                      &'air temperature at 2m (K)            '
      write (31,99017) 'q2m     ' ,                                     &
                      &'water vapor mixing ratio at 2m(kg/kg)'
      write (31,99017) 'ssw     ' ,                                     &
                      &'upper layer soil water               '
      write (31,99017) 'rsw     ' ,                                     &
                      &'root zone soil water                 '
      write (31,99017) 'tpr     ' ,                                     &
                      &'total precipitation (mm/day)         '
      write (31,99017) 'evp     ' ,                                     &
                      &'evapotranspiration (mm/day)          '
      write (31,99017) 'runoff  ' ,                                     &
                      &'surface runoff (mm/day)              '
      write (31,99017) 'scv     ' ,                                     &
                      &'total snow amount (mm)               '
      write (31,99017) 'sena    ' ,                                     &
                      &'sensible heat flux (W/m2)            '
      write (31,99017) 'flw     ' ,                                     &
                      &'net infrared energy flux (W/m2)      '
      write (31,99017) 'fsw     ' ,                                     &
                      &'net absorbed solar energy flux (W/m2)'
      write (31,99017) 'flwd    ' ,                                     &
                      &'downward infrared energy flux (W/m2) '
      write (31,99017) 'sina    ' ,                                     &
                      &'incident solar energy flux (W/m2)    '
      write (31,99017) 'prcv    ' ,                                     &
                      &'convective precipitation (mm/day)    '
      write (31,99017) 'psb     ' ,                                     &
                      &'surface pressure (hPa)               '
      write (31,99017) 'zpbl    ' ,                                     &
                      &'PBL layer height                     '
      write (31,99017) 'tgmax   ' ,                                     &
                      &'maximum ground temperature (K)       '
      write (31,99017) 'tgmin   ' ,                                     &
                      &'minimum ground temperature (K)       '
      write (31,99017) 't2max   ' ,                                     &
                      &'maximum 2m air temperature (K)       '
      write (31,99017) 't2min   ' ,                                     &
                      &'minimum 2m air temperature (K)       '
      write (31,99017) 'w10max  ' ,                                     &
                      &'maximum 10m wind speed (m/s)         '
      write (31,99017) 'ps_min  ' ,                                     &
                      &'minimum surface pressure (hPa)       '
      write (31,99020)
      close (31)
99001 format ('dset ^',a14)
99002 format ('title RegCM normal output variables')
99003 format ('options big_endian')
99004 format ('options little_endian')
99005 format ('undef -1.e34')
99006 format ('pdef ',i4,1x,i4,1x,'lccr',7(1x,f7.2),1x,2(f7.0,1x))
99007 format ('xdef ',i4,' linear ',f7.2,1x,f7.4)
99008 format ('ydef ',i4,' linear ',f7.2,1x,f7.4)
99009 format ('xdef ',i3,' linear ',f9.4,' ',f9.4)
99010 format ('ydef ',i3,' levels')
99011 format (10F7.2)
99012 format ('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
99013 format ('ydef ',i3,' linear ',f9.4,' ',f9.4)
99014 format ('zdef 1',' levels ',f7.2)
99015 format ('tdef ',i4,' linear ',i2,'z',a2,a3,i4,' ',i2,'hr')
99016 format ('vars ',i2)
99017 format (a8,'0 99 ',a36)
99018 format (a8,'0 33,105 ',a36)
99019 format (a8,'0 34,105 ',a36)
99020 format ('endvars')
      end subroutine gradsbat
!
      subroutine gradssub(ctlname)
      use mod_dynparam
      use mod_message , only : fatal
      use mod_date
      use mod_param1
      use mod_param2
      use mod_param3
      use mod_iunits
      use mod_main
#ifdef MPP1
      use mod_mppio
#endif
      implicit none
!
! Dummy arguments
!
      character(18) :: ctlname
!
! Local variables
!
      real(4) :: alatmax , alatmin , alonmax , alonmin , centeri ,      &
               & centerj , rlatinc , rloninc
      character(2) , dimension(31) :: cday
      character(3) , dimension(12) :: cmonth
      integer :: i , ifrq , j , jbend , mnend , month , myear , nbase , &
               & nday , nhour , nnumb , nx , ny
#ifdef MPP1
      real(4) , dimension(jxsg,iysg) :: xlat_s_io , xlon_s_io
#else
      real(4) , dimension(jxsg,iysg) :: xlat_s , xlon_s
#endif
!
      data cday/'01' , '02' , '03' , '04' , '05' , '06' , '07' , '08' , &
          &'09' , '10' , '11' , '12' , '13' , '14' , '15' , '16' ,      &
         & '17' , '18' , '19' , '20' , '21' , '22' , '23' , '24' ,      &
         & '25' , '26' , '27' , '28' , '29' , '30' , '31'/
      data cmonth/'jan' , 'feb' , 'mar' , 'apr' , 'may' , 'jun' ,       &
         & 'jul' , 'aug' , 'sep' , 'oct' , 'nov' , 'dec'/
!
      alatmin = 999999.
      alatmax = -999999.
      alonmin = 999999.
      alonmax = -999999.
      nx = 0
      ny = 0
      centerj = (jxm2)/2.
      centeri = (iym2)/2.

      open (31,file=trim(dirout)//pthsep//ctlname,status='replace')
      write (31,99001) ctlname(1:14)
      write (31,99002)
      if ( ibigend.eq.1 ) then
        write (31,99003)
      else if ( ibigend.eq.0 ) then
        write (31,99004)
      else
      end if
      if ( iotyp.eq.2 ) write (31,'(a)') 'options sequential'
      write (31,99005)
      if ( iproj.eq.'LAMCON' .or. iproj.eq.'ROTMER' ) then
        do j = 2 , jx-1
#ifdef MPP1
          if ( xlat_io(2,j).lt.alatmin ) alatmin = xlat_io(2,j)
          if ( xlat_io(iy-1,j).gt.alatmax ) alatmax = xlat_io(iy-1,j)
#else
          if ( xlat(2,j).lt.alatmin ) alatmin = xlat(2,j)
          if ( xlat(iy-1,j).gt.alatmax ) alatmax = xlat(iy-1,j)
#endif
        end do
        do i = 2 , iy-1
          do j = 2 , jx-1
#ifdef MPP1
            if ( clon.ge.0.0 ) then
              if ( xlong_io(i,j).ge.0.0 ) then
                alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
                alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
              else if ( abs(clon-xlong_io(i,j))                         &
                      & .lt.abs(clon-(xlong_io(i,j)+360.)) ) then
                alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
                alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
              else
                alonmin = amin1(alonmin,sngl(xlong_io(i,j))+360.)
                alonmax = amax1(alonmax,sngl(xlong_io(i,j))+360.)
              end if
            else if ( xlong_io(i,j).lt.0.0 ) then
              alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
              alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
            else if ( abs(clon-xlong_io(i,j))                           &
                    & .lt.abs(clon-(xlong_io(i,j)-360.)) ) then
              alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
              alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
            else
              alonmin = amin1(alonmin,sngl(xlong_io(i,j))-360.)
              alonmax = amax1(alonmax,sngl(xlong_io(i,j))-360.)
            end if
#else
            if ( clon.ge.0.0 ) then
              if ( xlong(i,j).ge.0.0 ) then
                alonmin = amin1(alonmin,sngl(xlong(i,j)))
                alonmax = amax1(alonmax,sngl(xlong(i,j)))
              else if ( abs(clon-xlong(i,j))                            &
                      & .lt.abs(clon-(xlong(i,j)+360.)) ) then
                alonmin = amin1(alonmin,sngl(xlong(i,j)))
                alonmax = amax1(alonmax,sngl(xlong(i,j)))
              else
                alonmin = amin1(alonmin,sngl(xlong(i,j))+360.)
                alonmax = amax1(alonmax,sngl(xlong(i,j))+360.)
              end if
            else if ( xlong(i,j).lt.0.0 ) then
              alonmin = amin1(alonmin,sngl(xlong(i,j)))
              alonmax = amax1(alonmax,sngl(xlong(i,j)))
            else if ( abs(clon-xlong(i,j))                              &
                    & .lt.abs(clon-(xlong(i,j)-360.)) ) then
              alonmin = amin1(alonmin,sngl(xlong(i,j)))
              alonmax = amax1(alonmax,sngl(xlong(i,j)))
            else
              alonmin = amin1(alonmin,sngl(xlong(i,j))-360.)
              alonmax = amax1(alonmax,sngl(xlong(i,j))-360.)
            end if
#endif
          end do
        end do
        rlatinc = dx*0.001/111./2.
        rloninc = dx*0.001/111./2.
        ny = 2 + nint(abs(alatmax-alatmin)/rlatinc)
        nx = 1 + nint(abs((alonmax-alonmin)/rloninc))
      end if
      if ( iotyp.eq.1 ) then
        if ( iproj.eq.'LAMCON' ) then   ! Lambert projection
          write (31,99006) jxm2sg , iym2sg , clat , clon ,              &
                         & centerj*nsg , centeri*nsg , truelatl ,       &
                         & truelath , clon , dx/nsg , dx/nsg
          write (31,99007) nx + 2 , alonmin - rloninc , rloninc
          write (31,99008) ny + 2 , alatmin - rlatinc , rlatinc
        else if ( iproj.eq.'POLSTR' ) then
                                        !
        else if ( iproj.eq.'NORMER' ) then
#ifdef MPP1
          read (iutin1,rec=5) xlat_s_io
          read (iutin1,rec=6) xlon_s_io
          write (31,99009) jxm2sg , xlon_s_io(nsg,nsg) ,                &
                         & xlon_s_io(nsg+1,nsg) - xlon_s_io(nsg,nsg)
          write (31,99010) iym2sg
          write (31,99011) (xlat_s_io(nsg+1,i),i=nsg+1,iym1sg)
#else
          read (iutin1,rec=5) xlat_s
          read (iutin1,rec=6) xlon_s
          write (31,99009) jxm2sg , xlon_s(nsg,nsg) ,                   &
                         & xlon_s(nsg+1,nsg) - xlon_s(nsg,nsg)
          write (31,99010) iym2sg
          write (31,99011) (xlat_s(nsg+1,i),i=nsg+1,iym1sg)
#endif
        else if ( iproj.eq.'ROTMER' ) then
          write (*,*) 'Note that rotated Mercartor (ROTMER)' ,          &
                     &' projections are not supported by GrADS.'
          write (*,*) '  Although not exact, the eta.u projection' ,    &
                     &' in GrADS is somewhat similar.'
          write (*,*) ' FERRET, however, does support this projection.'
          write (31,99012) jxm2sg , iym2sg , plon , plat ,              &
                         & dx/111000./nsg , dx/111000.*.95238/nsg
          write (31,99007) nx + 2 , alonmin - rloninc , rloninc
          write (31,99008) ny + 2 , alatmin - rlatinc , rlatinc
        else
          call fatal(__FILE__,__LINE__,'INVALID MAP PROJECTION')
        end if
      else if ( iotyp.eq.2 ) then
#ifdef MPP1
        write (31,99009) jxm2sg , xlong_io(2,2) ,                       &
                       & (xlong_io(2,3)-xlong_io(2,2))/nsg
        write (31,99013) iym2sg , xlat_io(2,2) ,                        &
                       & (xlat_io(3,2)-xlat_io(2,2))/nsg
#else
        write (31,99009) jxm2 , xlong(2,2) , xlong(2,3) - xlong(2,2)
        write (31,99013) iym2 , xlat(2,2) , xlat(3,2) - xlat(2,2)
#endif
      else
      end if
      write (31,99014) 1 , (1013.25-r8pt*10.)*a(kz) + r8pt*10.
      myear = ldatez/1000000
      month = (ldatez-myear*1000000)/10000
      nday = (ldatez-myear*1000000-month*10000)/100
      nhour = mod(ldatez,100)
      call finddate(nbase,ldatez)
      if ( month.eq.12 ) then
        call finddate(mnend,myear*1000000+1010100)
      else
        call finddate(mnend,myear*1000000+month*10000+10100)
      end if
      call finddate(jbend,idate2)
      if ( ldatez.eq.idate0 ) then
        nnumb = (ibdyfrq/batfrq+0.00001)*(min0(jbend,mnend)-nbase) + 1
      else
        nnumb = (ibdyfrq/batfrq+0.00001)*(min0(jbend,mnend)-nbase)
      end if
      ifrq = batfrq + 0.00001
      if ( ldatez.eq.idate0 ) then
        write (31,99015) nnumb , nhour , cday(nday) , cmonth(month) ,   &
                       & myear , ifrq
      else
        write (31,99015) nnumb , nhour + ifrq , cday(nday) ,            &
                       & cmonth(month) , myear , ifrq
      end if
      if ( iotyp.eq.2 ) write (31,'(a)') 'theader 4'
      write (31,99016) 16
      if ( iproj.eq.'LAMCON' .and. iotyp.eq.1 ) then   ! Lambert projection
        write (31,99018) 'u10m    ' ,                                   &
                        &'westerly  wind at 10m (m/s)          '
        write (31,99019) 'v10m    ' ,                                   &
                        &'southerly wind at 10m (m/s)          '
      else
        write (31,99017) 'u10m    ' ,                                   &
                        &'westerly  wind at 10m (m/s)          '
        write (31,99017) 'v10m    ' ,                                   &
                        &'southerly wind at 10m (m/s)          '
      end if
      write (31,99017) 'uvdrag  ' ,                                     &
                      &'surface drag stress                  '
      write (31,99017) 'tg      ' ,                                     &
                      &'ground temperature (degree)          '
      write (31,99017) 'tlef    ' ,                                     &
                      &'temperature of foliage               '
      write (31,99017) 't2m     ' ,                                     &
                      &'air temperature at 2m (K)            '
      write (31,99017) 'q2m     ' ,                                     &
                      &'water vapor mixing ratio at 2m(kg/kg)'
      write (31,99017) 'ssw     ' ,                                     &
                      &'upper layer soil water               '
      write (31,99017) 'rsw     ' ,                                     &
                      &'root zone soil water                 '
      write (31,99017) 'tpr     ' ,                                     &
                      &'total precipitation (mm/day)         '
      write (31,99017) 'evp     ' ,                                     &
                      &'evapotranspiration (mm/day)          '
      write (31,99017) 'runoff  ' ,                                     &
                      &'surface runoff (mm/day)              '
      write (31,99017) 'scv     ' ,                                     &
                      &'total snow amount (mm)               '
      write (31,99017) 'sena    ' ,                                     &
                      &'sensible heat flux (W/m2)            '
      write (31,99017) 'prcv    ' ,                                     &
                      &'convective precipitation (mm/day)    '
      write (31,99017) 'ps      ' ,                                     &
                      &'surface pressure (hPa)               '
      write (31,99020)
      close (31)
99001 format ('dset ^',a14)
99002 format ('title RegCM normal output variables')
99003 format ('options big_endian')
99004 format ('options little_endian')
99005 format ('undef -1.e34')
99006 format ('pdef ',i4,1x,i4,1x,'lccr',7(1x,f7.2),1x,2(f7.0,1x))
99007 format ('xdef ',i4,' linear ',f7.2,1x,f7.4)
99008 format ('ydef ',i4,' linear ',f7.2,1x,f7.4)
99009 format ('xdef ',i3,' linear ',f9.4,' ',f9.4)
99010 format ('ydef ',i3,' levels')
99011 format (10F7.2)
99012 format ('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
99013 format ('ydef ',i3,' linear ',f9.4,' ',f9.4)
99014 format ('zdef ',i2,' levels ',30F7.2)
99015 format ('tdef ',i4,' linear ',i2,'z',a2,a3,i4,' ',i2,'hr')
99016 format ('vars ',i2)
99017 format (a8,'0 99 ',a36)
99018 format (a8,'0 33,105 ',a36)
99019 format (a8,'0 34,105 ',a36)
99020 format ('endvars')
      end subroutine gradssub
!
      subroutine gradschem(ctlname)
      use mod_dynparam
      use mod_message , only : fatal
      use mod_date
      use mod_param1
      use mod_param2
      use mod_param3
      use mod_main
#ifdef MPP1
      use mod_mppio
#endif
      implicit none
!
! Dummy arguments
!
      character(18) :: ctlname
!
! Local variables
!
      real(4) :: alatmax , alatmin , alonmax , alonmin , centeri ,      &
               & centerj , rlatinc , rloninc
      character(2) , dimension(31) :: cday
      character(3) , dimension(12) :: cmonth
      integer :: i , ifrq , itr , j , jbend , k , mnend , month ,       &
               & myear , nbase , nday , nhour , nnumb , nx , ny
!
      data cday/'01' , '02' , '03' , '04' , '05' , '06' , '07' , '08' , &
          &'09' , '10' , '11' , '12' , '13' , '14' , '15' , '16' ,      &
         & '17' , '18' , '19' , '20' , '21' , '22' , '23' , '24' ,      &
         & '25' , '26' , '27' , '28' , '29' , '30' , '31'/
      data cmonth/'jan' , 'feb' , 'mar' , 'apr' , 'may' , 'jun' ,       &
         & 'jul' , 'aug' , 'sep' , 'oct' , 'nov' , 'dec'/
!
      alatmin = 999999.
      alatmax = -999999.
      alonmin = 999999.
      alonmax = -999999.
      nx = 0
      ny = 0
      centerj = (jxm2)/2.
      centeri = (iym2)/2.

      open (31,file=trim(dirout)//pthsep//ctlname,status='replace')
      write (31,99001) ctlname(1:14)
      write (31,99002)
      if ( ibigend.eq.1 ) then
        write (31,99003)
      else if ( ibigend.eq.0 ) then
        write (31,99004)
      else
      end if
      if ( iotyp.eq.2 ) write (31,'(a)') 'options sequential'
      write (31,99005)
      if ( iproj.eq.'LAMCON' .or. iproj.eq.'ROTMER' ) then
        do j = 2 , jx-1
#ifdef MPP1
          if ( xlat_io(2,j).lt.alatmin ) alatmin = xlat_io(2,j)
          if ( xlat_io(iy-1,j).gt.alatmax ) alatmax = xlat_io(iy-1,j)
#else
          if ( xlat(2,j).lt.alatmin ) alatmin = xlat(2,j)
          if ( xlat(iy-1,j).gt.alatmax ) alatmax = xlat(iy-1,j)
#endif
        end do
        do i = 2 , iy-1
          do j = 2 , jx-1
#ifdef MPP1
            if ( clon.ge.0.0 ) then
              if ( xlong_io(i,j).ge.0.0 ) then
                alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
                alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
              else if ( abs(clon-xlong_io(i,j))                         &
                      & .lt.abs(clon-(xlong_io(i,j)+360.)) ) then
                alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
                alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
              else
                alonmin = amin1(alonmin,sngl(xlong_io(i,j))+360.)
                alonmax = amax1(alonmax,sngl(xlong_io(i,j))+360.)
              end if
            else if ( xlong_io(i,j).lt.0.0 ) then
              alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
              alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
            else if ( abs(clon-xlong_io(i,j))                           &
                    & .lt.abs(clon-(xlong_io(i,j)-360.)) ) then
              alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
              alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
            else
              alonmin = amin1(alonmin,sngl(xlong_io(i,j))-360.)
              alonmax = amax1(alonmax,sngl(xlong_io(i,j))-360.)
            end if
#else
            if ( clon.ge.0.0 ) then
              if ( xlong(i,j).ge.0.0 ) then
                alonmin = amin1(alonmin,sngl(xlong(i,j)))
                alonmax = amax1(alonmax,sngl(xlong(i,j)))
              else if ( abs(clon-xlong(i,j))                            &
                      & .lt.abs(clon-(xlong(i,j)+360.)) ) then
                alonmin = amin1(alonmin,sngl(xlong(i,j)))
                alonmax = amax1(alonmax,sngl(xlong(i,j)))
              else
                alonmin = amin1(alonmin,sngl(xlong(i,j))+360.)
                alonmax = amax1(alonmax,sngl(xlong(i,j))+360.)
              end if
            else if ( xlong(i,j).lt.0.0 ) then
              alonmin = amin1(alonmin,sngl(xlong(i,j)))
              alonmax = amax1(alonmax,sngl(xlong(i,j)))
            else if ( abs(clon-xlong(i,j))                              &
                    & .lt.abs(clon-(xlong(i,j)-360.)) ) then
              alonmin = amin1(alonmin,sngl(xlong(i,j)))
              alonmax = amax1(alonmax,sngl(xlong(i,j)))
            else
              alonmin = amin1(alonmin,sngl(xlong(i,j))-360.)
              alonmin = amin1(alonmin,sngl(xlong(i,j))-360.)
              alonmax = amax1(alonmax,sngl(xlong(i,j))-360.)
            end if
#endif
          end do
        end do
        rlatinc = dx*0.001/111./2.
        rloninc = dx*0.001/111./2.
        ny = 2 + nint(abs(alatmax-alatmin)/rlatinc)
        nx = 1 + nint(abs((alonmax-alonmin)/rloninc))
      end if
      if ( iotyp.eq.1 ) then
        if ( iproj.eq.'LAMCON' ) then   ! Lambert projection
          write (31,99006) jxm2 , iym2 , clat , clon , centerj ,        &
                         & centeri , truelatl , truelath , clon , dx ,&
                         & dx
          write (31,99007) nx + 2 , alonmin - rloninc , rloninc
          write (31,99008) ny + 2 , alatmin - rlatinc , rlatinc
        else if ( iproj.eq.'POLSTR' ) then
                                        !
        else if ( iproj.eq.'NORMER' ) then
#ifdef MPP1
          write (31,99009) jxm2 , xlong_io(2,2) , xlong_io(2,3)         &
                         & - xlong_io(2,2)
          write (31,99010) iym2
          write (31,99011) (xlat_io(i,2),i=2,iym1)
#else
          write (31,99009) jxm2 , xlong(2,2) , xlong(2,3) - xlong(2,2)
          write (31,99010) iym2
          write (31,99011) (xlat(i,2),i=2,iym1)
#endif
        else if ( iproj.eq.'ROTMER' ) then
          write (*,*) 'Note that rotated Mercartor (ROTMER)' ,          &
                     &' projections are not supported by GrADS.'
          write (*,*) '  Although not exact, the eta.u projection' ,    &
                     &' in GrADS is somewhat similar.'
          write (*,*) ' FERRET, however, does support this projection.'
          write (31,99012) jxm2 , iym2 , plon , plat ,                  &
                         & dx/111000. , dx/111000.*.95238
          write (31,99007) nx + 2 , alonmin - rloninc , rloninc
          write (31,99008) ny + 2 , alatmin - rlatinc , rlatinc
        else
          call fatal(__FILE__,__LINE__,'INVALID MAP PROJECTION')
        end if
      else if ( iotyp.eq.2 ) then
#ifdef MPP1
        write (31,99009) jxm2 , xlong_io(2,2) , xlong_io(2,3)           &
                       & - xlong_io(2,2)
        write (31,99013) iym2 , xlat_io(2,2) , xlat_io(3,2)             &
                       & - xlat_io(2,2)
#else
        write (31,99009) jxm2 , xlong(2,2) , xlong(2,3) - xlong(2,2)
        write (31,99013) iym2 , xlat(2,2) , xlat(3,2) - xlat(2,2)
#endif
      else
      end if
      write (31,99014) kz , ((1013.25-r8pt*10.)*a(k)+r8pt*10.,k=kz,1,-1)
      myear = ldatez/1000000
      month = (ldatez-myear*1000000)/10000
      nday = (ldatez-myear*1000000-month*10000)/100
      nhour = mod(ldatez,100)
      call finddate(nbase,ldatez)
      if ( month.eq.12 ) then
        call finddate(mnend,myear*1000000+1010100)
      else
        call finddate(mnend,myear*1000000+month*10000+10100)
      end if
      call finddate(jbend,idate2)
      if ( ldatez.eq.idate0 ) then
        nnumb = (ibdyfrq/tapfrq+0.00001)*(min0(jbend,mnend)-nbase) + 1
      else
        nnumb = (ibdyfrq/tapfrq+0.00001)*(min0(jbend,mnend)-nbase)
      end if
      ifrq = chemfrq + 0.00001
      if ( ldatez.eq.idate0 ) then
        write (31,99015) nnumb , nhour , cday(nday) , cmonth(month) ,   &
                       & myear , ifrq
      else
        write (31,99015) nnumb , nhour + ifrq , cday(nday) ,            &
                       & cmonth(month) , myear , ifrq
      end if
      if ( iotyp.eq.2 ) write (31,'(a)') 'theader 4'
      write (31,99016) ntr + 3 + 7*ntr + 2 + 1
 
      do itr = 1 , ntr
        if ( itr.lt.10 ) then
          write (31,99018) 'trac' , itr , kz ,                          &
                          &'tracer mix. rat  (Kg/Kg)'
        else
          write (31,99019) 'trac' , itr , kz ,                          &
                          &'tracer mix. rat  (Kg/Kg)'
        end if
      end do
      write (31,99018) 'aext' , 8 , kz , 'aer mix. ext. coef      '
      write (31,99018) 'assa' , 8 , kz , 'aer mix. sin. scat. alb '
      write (31,99018) 'agfu' , 8 , kz , 'aer mix. ass. par       '
 
      do itr = 1 , ntr
        if ( itr.lt.10 ) then
          write (31,99020) 'colb__tr' , itr , 'columnburden inst(mg/m2)'
          write (31,99020) 'wdlsc_tr' , itr , 'wet dep lgscale(mg/m2/d)'
          write (31,99020) 'wdcvc_tr' , itr , 'wet dep convect(mg/m2/d)'
          write (31,99020) 'sdrdp_tr' , itr , 'surf dry depos.(mg/m2/d)'
          write (31,99020) 'xgasc_tr' , itr , 'chem gas conv. (mg/m2/d)'
          write (31,99020) 'xaquc_tr' , itr , 'chem aqu conv. (mg/m2/d)'
          write (31,99020) 'emiss_tr' , itr , 'surf emission  (mg/m2/d)'
        else
          write (31,99021) 'colb__tr' , itr , 'columnburden inst(mg/m2)'
          write (31,99021) 'wdlsc_tr' , itr , 'wet dep lgscale(mg/m2/d)'
          write (31,99021) 'wdcvc_tr' , itr , 'wet dep convect(mg/m2/d)'
          write (31,99021) 'sdrdp_tr' , itr , 'surf dry depos.(mg/m2/d)'
          write (31,99021) 'xgasc_tr' , itr , 'chem gas conv. (mg/m2/d)'
          write (31,99021) 'xaquc_tr' , itr , 'chem aqu conv. (mg/m2/d)'
          write (31,99021) 'emiss_tr' , itr , 'surf emission  (mg/m2/d)'
        end if
      end do
      write (31,99017) 'acstoarf' , ' TOArad forcing av.(W/m2)'
      write (31,99017) 'acstsrrf' , ' SRFrad forcing av.(W/m2)'
      write (31,99017) 'psa' , ' Surface Pressure (hPa)'
 
      write (31,99022)
      close (31)
99001 format ('dset ^',a14)
99002 format ('title RegCM chemistry/tracor variables')
99003 format ('options big_endian')
99004 format ('options little_endian')
99005 format ('undef -1.e34')
99006 format ('pdef ',i4,1x,i4,1x,'lcc',7(1x,f7.2),1x,2(f7.0,1x))
99007 format ('xdef ',i4,' linear ',f7.2,1x,f7.4)
99008 format ('ydef ',i4,' linear ',f7.2,1x,f7.4)
99009 format ('xdef ',i3,' linear ',f9.4,' ',f9.4)
99010 format ('ydef ',i3,' levels')
99011 format (10F7.2)
99012 format ('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
99013 format ('ydef ',i3,' linear ',f9.4,' ',f9.4)
99014 format ('zdef ',i2,' levels ',30F7.2)
99015 format ('tdef ',i4,' linear ',i2,'z',a2,a3,i4,' ',i2,'hr')
99016 format ('vars ',i2)
99017 format (a8,' 0 99 ',a26)
99018 format (a4,i1,' ',i2,' 0 ',a26)
99019 format (a4,i2,' ',i2,' 0 ',a26)
99020 format (a8,i1,' 0 99 ',a26)
99021 format (a8,i2,' 0 99 ',a26)
99022 format ('endvars')
      end subroutine gradschem
!
      subroutine gradsctl(ctlname)
      use mod_dynparam
      use mod_message , only : fatal
      use mod_date
      use mod_param1
      use mod_param2
      use mod_main
#ifdef MPP1
      use mod_mppio
#endif
      implicit none
!
! Dummy arguments
!
      character(12) :: ctlname
!
! Local variables
!
      real(4) :: alatmax , alatmin , alonmax , alonmin , centeri ,      &
               & centerj , rlatinc , rloninc
      integer :: i , j , nx , ny
!
      alatmin = 999999.
      alatmax = -999999.
      alonmin = 999999.
      alonmax = -999999.
      nx = 0
      ny = 0
      centerj = (jxm2)/2.
      centeri = (iym2)/2.

      open (31,file=trim(dirout)//pthsep//ctlname,status='replace')
      write (31,99001)
      write (31,99002)
      if ( ibigend.eq.1 ) then
        write (31,99003)
      else if ( ibigend.eq.0 ) then
        write (31,99004)
      else
      end if
      write (31,99005)
      if ( iproj.eq.'LAMCON' .or. iproj.eq.'ROTMER' ) then
        do j = 2 , jx-1
#ifdef MPP1
          if ( xlat_io(2,j).lt.alatmin ) alatmin = xlat_io(2,j)
          if ( xlat_io(iy-1,j).gt.alatmax ) alatmax = xlat_io(iy-1,j)
#else
          if ( xlat(2,j).lt.alatmin ) alatmin = xlat(2,j)
          if ( xlat(iy-1,j).gt.alatmax ) alatmax = xlat(iy-1,j)
#endif
        end do
        do i = 2 , iy-1
          do j = 2 , jx-1
#ifdef MPP1
            if ( clon.ge.0.0 ) then
              if ( xlong_io(i,j).ge.0.0 ) then
                alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
                alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
              else if ( abs(clon-xlong_io(i,j))                         &
                      & .lt.abs(clon-(xlong_io(i,j)+360.)) ) then
                alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
                alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
              else
                alonmin = amin1(alonmin,sngl(xlong_io(i,j))+360.)
                alonmax = amax1(alonmax,sngl(xlong_io(i,j))+360.)
              end if
            else if ( xlong_io(i,j).lt.0.0 ) then
              alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
              alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
            else if ( abs(clon-xlong_io(i,j))                           &
                    & .lt.abs(clon-(xlong_io(i,j)-360.)) ) then
              alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
              alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
            else
              alonmin = amin1(alonmin,sngl(xlong_io(i,j))-360.)
              alonmax = amax1(alonmax,sngl(xlong_io(i,j))-360.)
            end if
#else
            if ( clon.ge.0.0 ) then
              if ( xlong(i,j).ge.0.0 ) then
                alonmin = amin1(alonmin,sngl(xlong(i,j)))
                alonmax = amax1(alonmax,sngl(xlong(i,j)))
              else if ( abs(clon-xlong(i,j))                            &
                      & .lt.abs(clon-(xlong(i,j)+360.)) ) then
                alonmin = amin1(alonmin,sngl(xlong(i,j)))
                alonmax = amax1(alonmax,sngl(xlong(i,j)))
              else
                alonmin = amin1(alonmin,sngl(xlong(i,j))+360.)
                alonmax = amax1(alonmax,sngl(xlong(i,j))+360.)
              end if
            else if ( xlong(i,j).lt.0.0 ) then
              alonmin = amin1(alonmin,sngl(xlong(i,j)))
              alonmax = amax1(alonmax,sngl(xlong(i,j)))
            else if ( abs(clon-xlong(i,j))                              &
                    & .lt.abs(clon-(xlong(i,j)-360.)) ) then
              alonmin = amin1(alonmin,sngl(xlong(i,j)))
              alonmax = amax1(alonmax,sngl(xlong(i,j)))
            else
              alonmin = amin1(alonmin,sngl(xlong(i,j))-360.)
              alonmax = amax1(alonmax,sngl(xlong(i,j))-360.)
            end if
#endif
          end do
        end do
        rlatinc = dx*0.001/111./2.
        rloninc = dx*0.001/111./2.
        ny = 2 + nint(abs(alatmax-alatmin)/rlatinc)
        nx = 1 + nint(abs((alonmax-alonmin)/rloninc))
      end if
      if ( iproj.eq.'LAMCON' ) then     ! Lambert projection
        write (31,99006) jxm2 , iym2 , clat , clon , centerj ,          &
                       & centeri , truelatl , truelath , clon , dx ,  &
                       & dx
        write (31,99007) nx + 2 , alonmin - rloninc , rloninc
        write (31,99008) ny + 2 , alatmin - rlatinc , rlatinc
      else if ( iproj.eq.'POLSTR' ) then !
      else if ( iproj.eq.'NORMER' ) then
#ifdef MPP1
        write (31,99009) jxm2 , xlong_io(2,2) , xlong_io(2,3)           &
                       & - xlong_io(2,2)
        write (31,99010) iym2
        write (31,99011) (xlat_io(i,2),i=2,iym1)
#else
        write (31,99009) jxm2 , xlong(2,2) , xlong(2,3) - xlong(2,2)
        write (31,99010) iym2
        write (31,99011) (xlat(i,2),i=2,iym1)
#endif
      else if ( iproj.eq.'ROTMER' ) then
        write (*,*) 'Note that rotated Mercartor (ROTMER)' ,            &
                   &' projections are not supported by GrADS.'
        write (*,*) '  Although not exact, the eta.u projection' ,      &
                   &' in GrADS is somewhat similar.'
        write (*,*) ' FERRET, however, does support this projection.'
       write (31,99012) jxm2 , iym2 , plon , plat , dx/111000. ,      &
                       & dx/111000.*.95238
        write (31,99007) nx + 2 , alonmin - rloninc , rloninc
        write (31,99008) ny + 2 , alatmin - rlatinc , rlatinc

      else
        call fatal(__FILE__,__LINE__,'INVALID MAP PROJECTION')
      end if
      write (31,99013) 1 , 1000.
      write (31,99014) 1
      write (31,99015) 11
      write (31,99016) 'head    ' , 'header information         '
      write (31,99016) 'ht      ' , 'surface elevation          '
      write (31,99016) 'htsd    ' , 'surface elevation std dev  '
      write (31,99016) 'veg2d   ' , 'vegetation type in BATS    '
      write (31,99016) 'landuse ' , 'surface landuse mod_type       '
      write (31,99016) 'xlat    ' , 'latitude  of cross points  '
      write (31,99016) 'xlong   ' , 'longitude of cross points  '
      write (31,99016) 'xmap    ' , 'map factors of cross points'
      write (31,99016) 'dmap    ' , 'map factors of dot points  '
      write (31,99016) 'coriol  ' , 'coriol force               '
      write (31,99016) 'mask    ' , 'land/sea mask              '
      write (31,99017)
      close (31)
99001 format ('dset ^OUT_HEAD')
99002 format ('title RegCM domain information')
99003 format ('options big_endian')
99004 format ('options little_endian')
99005 format ('undef -1.e34')
99006 format ('pdef ',i4,1x,i4,1x,'lcc',7(1x,f7.2),1x,2(f7.0,1x))
99007 format ('xdef ',i4,' linear ',f7.2,1x,f7.4)
99008 format ('ydef ',i4,' linear ',f7.2,1x,f7.4)
99009 format ('xdef ',i3,' linear ',f9.4,' ',f9.4)
99010 format ('ydef ',i3,' levels')
99011 format (10F7.2)
99012 format ('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
99013 format ('zdef ',i1,' levels ',f7.2)
99014 format ('tdef ',i1,' linear 00z01Jan2001 1mo')
99015 format ('vars ',i2)
99016 format (a8,'0 99 ',a26)
99017 format ('endvars')
      end subroutine gradsctl
!
      subroutine gradsout(ctlname)
      use mod_dynparam
      use mod_message , only : fatal
      use mod_date
      use mod_param1
      use mod_param2
      use mod_param3
      use mod_main
#ifdef MPP1
      use mod_mppio
#endif
      implicit none
!
! Dummy arguments
!
      character(18) :: ctlname
!
! Local variables
!
      real(4) :: alatmax , alatmin , alonmax , alonmin , centeri ,      &
               & centerj , rlatinc , rloninc
      character(2) , dimension(31) :: cday
      character(3) , dimension(12) :: cmonth
      integer :: i , ifrq , j , jbend , k , mnend , month , myear ,     &
               & nbase , nday , nhour , nnumb , nx , ny
!
      data cday/'01' , '02' , '03' , '04' , '05' , '06' , '07' , '08' , &
          &'09' , '10' , '11' , '12' , '13' , '14' , '15' , '16' ,      &
         & '17' , '18' , '19' , '20' , '21' , '22' , '23' , '24' ,      &
         & '25' , '26' , '27' , '28' , '29' , '30' , '31'/
      data cmonth/'jan' , 'feb' , 'mar' , 'apr' , 'may' , 'jun' ,       &
         & 'jul' , 'aug' , 'sep' , 'oct' , 'nov' , 'dec'/
!
      alatmin = 999999.
      alatmax = -999999.
      alonmin = 999999.
      alonmax = -999999.
      nx = 0
      ny = 0
      centerj = (jxm2)/2.
      centeri = (iym2)/2.

      open (31,file=trim(dirout)//pthsep//ctlname,status='replace')
      write (31,99001) ctlname(1:14)
      write (31,99002)
      if ( ibigend.eq.1 ) then
        write (31,99003)
      else if ( ibigend.eq.0 ) then
        write (31,99004)
      else
      end if
      if ( iotyp.eq.2 ) write (31,'(a)') 'options sequential'
      write (31,99005)
      if ( iproj.eq.'LAMCON' .or. iproj.eq.'ROTMER' ) then
        do j = 2 , jx-1
#ifdef MPP1
          if ( xlat_io(2,j).lt.alatmin ) alatmin = xlat_io(2,j)
          if ( xlat_io(iy-1,j).gt.alatmax ) alatmax = xlat_io(iy-1,j)
#else
          if ( xlat(2,j).lt.alatmin ) alatmin = xlat(2,j)
          if ( xlat(iy-1,j).gt.alatmax ) alatmax = xlat(iy-1,j)
#endif
        end do
        do i = 2 , iy-1
          do j = 2 , jx-1
#ifdef MPP1
            if ( clon.ge.0.0 ) then
              if ( xlong_io(i,j).ge.0.0 ) then
                alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
                alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
              else if ( abs(clon-xlong_io(i,j))                         &
                      & .lt.abs(clon-(xlong_io(i,j)+360.)) ) then
                alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
                alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
              else
                alonmin = amin1(alonmin,sngl(xlong_io(i,j))+360.)
                alonmax = amax1(alonmax,sngl(xlong_io(i,j))+360.)
              end if
            else if ( xlong_io(i,j).lt.0.0 ) then
              alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
              alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
            else if ( abs(clon-xlong_io(i,j))                           &
                    & .lt.abs(clon-(xlong_io(i,j)-360.)) ) then
              alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
              alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
            else
              alonmin = amin1(alonmin,sngl(xlong_io(i,j))-360.)
              alonmax = amax1(alonmax,sngl(xlong_io(i,j))-360.)
            end if
#else
            if ( clon.ge.0.0 ) then
              if ( xlong(i,j).ge.0.0 ) then
                alonmin = amin1(alonmin,sngl(xlong(i,j)))
                alonmax = amax1(alonmax,sngl(xlong(i,j)))
              else if ( abs(clon-xlong(i,j))                            &
                      & .lt.abs(clon-(xlong(i,j)+360.)) ) then
                alonmin = amin1(alonmin,sngl(xlong(i,j)))
                alonmax = amax1(alonmax,sngl(xlong(i,j)))
              else
                alonmin = amin1(alonmin,sngl(xlong(i,j))+360.)
                alonmax = amax1(alonmax,sngl(xlong(i,j))+360.)
              end if
            else if ( xlong(i,j).lt.0.0 ) then
              alonmin = amin1(alonmin,sngl(xlong(i,j)))
              alonmax = amax1(alonmax,sngl(xlong(i,j)))
            else if ( abs(clon-xlong(i,j))                              &
                    & .lt.abs(clon-(xlong(i,j)-360.)) ) then
              alonmin = amin1(alonmin,sngl(xlong(i,j)))
              alonmax = amax1(alonmax,sngl(xlong(i,j)))
            else
              alonmin = amin1(alonmin,sngl(xlong(i,j))-360.)
              alonmax = amax1(alonmax,sngl(xlong(i,j))-360.)
            end if
#endif
          end do
        end do
        rlatinc = dx*0.001/111./2.
        rloninc = dx*0.001/111./2.
        ny = 2 + nint(abs(alatmax-alatmin)/rlatinc)
        nx = 1 + nint(abs((alonmax-alonmin)/rloninc))
      end if
      if ( iotyp.eq.1 ) then
        if ( iproj.eq.'LAMCON' ) then   ! Lambert projection
          write (31,99006) jxm2 , iym2 , clat , clon , centerj ,        &
                         & centeri , truelatl , truelath , clon , dx ,&
                         & dx
          write (31,99007) nx + 2 , alonmin - rloninc , rloninc
          write (31,99008) ny + 2 , alatmin - rlatinc , rlatinc
        else if ( iproj.eq.'POLSTR' ) then
                                        !
        else if ( iproj.eq.'NORMER' ) then
#ifdef MPP1
          write (31,99009) jxm2 , xlong_io(2,2) , xlong_io(2,3)         &
                         & - xlong_io(2,2)
          write (31,99010) iym2
          write (31,99011) (xlat_io(i,2),i=2,iym1)
#else
          write (31,99009) jxm2 , xlong(2,2) , xlong(2,3) - xlong(2,2)
          write (31,99010) iym2
          write (31,99011) (xlat(i,2),i=2,iym1)
#endif
        else if ( iproj.eq.'ROTMER' ) then
          write (*,*) 'Note that rotated Mercartor (ROTMER)' ,          &
                     &' projections are not supported by GrADS.'
          write (*,*) '  Although not exact, the eta.u projection' ,    &
                     &' in GrADS is somewhat similar.'
          write (*,*) ' FERRET, however, does support this projection.'
          write (31,99012) jxm2 , iym2 , plon , plat ,                  &
                         & dx/111000. , dx/111000.*.95238
          write (31,99007) nx + 2 , alonmin - rloninc , rloninc
          write (31,99008) ny + 2 , alatmin - rlatinc , rlatinc
        else
          call fatal(__FILE__,__LINE__,'INVALID MAP PROJECTION')
        end if
      else if ( iotyp.eq.2 ) then
#ifdef MPP1
        write (31,99009) jxm2 , xlong_io(2,2) , xlong_io(2,3)           &
                       & - xlong_io(2,2)
        write (31,99013) iym2 , xlat_io(2,2) , xlat_io(3,2)             &
                       & - xlat_io(2,2)
#else
        write (31,99009) jxm2 , xlong(2,2) , xlong(2,3) - xlong(2,2)
        write (31,99013) iym2 , xlat(2,2) , xlat(3,2) - xlat(2,2)
#endif
      else
      end if
      write (31,99014) kz , ((1013.25-r8pt*10.)*a(k)+r8pt*10.,k=kz,1,-1)
      myear = ldatez/1000000
      month = (ldatez-myear*1000000)/10000
      nday = (ldatez-myear*1000000-month*10000)/100
      nhour = mod(ldatez,100)
      call finddate(nbase,ldatez)
      if ( month.eq.12 ) then
        call finddate(mnend,myear*1000000+1010100)
      else
        call finddate(mnend,myear*1000000+month*10000+10100)
      end if
      call finddate(jbend,idate2)
      if ( ldatez.eq.idate0 ) then
        nnumb = (ibdyfrq/tapfrq+0.00001)*(min0(jbend,mnend)-nbase) + 1
      else
        nnumb = (ibdyfrq/tapfrq+0.00001)*(min0(jbend,mnend)-nbase)
      end if
      ifrq = tapfrq + 0.00001
      if ( ldatez.eq.idate0 ) then
        write (31,99015) nnumb , nhour , cday(nday) , cmonth(month) ,   &
                       & myear , ifrq
      else
        write (31,99015) nnumb , nhour + ifrq , cday(nday) ,            &
                       & cmonth(month) , myear , ifrq
      end if
      if ( iotyp.eq.2 ) write (31,'(a)') 'theader 4'
      write (31,99016) 10 + 1
      if ( iproj.eq.'LAMCON' .and. iotyp.eq.1 ) then   ! Lambert projection
        write (31,99019) 'u       ' , kz , 'westerly wind (m/s)        '
        write (31,99020) 'v       ' , kz , 'southerly wind (m/s)       '
      else
        write (31,99018) 'u       ' , kz , 'westerly wind (m/s)        '
        write (31,99018) 'v       ' , kz , 'southerly wind (m/s)       '
      end if
      write (31,99018) 'w       ' , kz , 'omega (hPa/s)   p-velocity '
      write (31,99018) 't       ' , kz , 'air temperature (degree)   '
      write (31,99018) 'qv      ' , kz , 'water vapor mixing ratio   '
      write (31,99018) 'qc      ' , kz , 'cloud water mixing ratio   '
      write (31,99017) 'psa     ' , 'surface pressure (hPa)     '
      write (31,99017) 'tpr     ' , 'total precipitation(mm/day)'
      write (31,99017) 'tgb     ' , 'lower groud temp. in BATS  '
      write (31,99017) 'swt     ' , 'total soil water in mm H2O '
      write (31,99017) 'rno     ' , 'accumulated infiltration   '
      write (31,99021)
      close (31)
99001 format ('dset ^',a14)
99002 format ('title RegCM normal output variables')
99003 format ('options big_endian')
99004 format ('options little_endian')
99005 format ('undef -1.e34')
99006 format ('pdef ',i4,1x,i4,1x,'lccr',7(1x,f7.2),1x,2(f7.0,1x))
99007 format ('xdef ',i4,' linear ',f7.2,1x,f7.4)
99008 format ('ydef ',i4,' linear ',f7.2,1x,f7.4)
99009 format ('xdef ',i3,' linear ',f9.4,' ',f9.4)
99010 format ('ydef ',i3,' levels')
99011 format (10F7.2)
99012 format ('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
99013 format ('ydef ',i3,' linear ',f9.4,' ',f9.4)
99014 format ('zdef ',i2,' levels ',30F7.2)
99015 format ('tdef ',i4,' linear ',i2,'z',a2,a3,i4,' ',i2,'hr')
99016 format ('vars ',i2)
99017 format (a8,'0 99 ',a26)
99018 format (a8,i2,' 0 ',a26)
99019 format (a8,i2,' 33,100 ',a36)
99020 format (a8,i2,' 34,100 ',a36)
99021 format ('endvars')
      end subroutine gradsout
!
      subroutine gradsrad(ctlname)
      use mod_dynparam
      use mod_message , only : fatal
      use mod_date
      use mod_param1
      use mod_param2
      use mod_param3
      use mod_main
#ifdef MPP1
      use mod_mppio
#endif
      implicit none
!
! Dummy arguments
!
      character(18) :: ctlname
!
! Local variables
!
      real(4) :: alatmax , alatmin , alonmax , alonmin , centeri ,      &
               & centerj , rlatinc , rloninc
      character(2) , dimension(31) :: cday
      character(3) , dimension(12) :: cmonth
      integer :: i , ifrq , j , jbend , k , mnend , month , myear ,     &
               & nbase , nday , nhour , nnumb , nx , ny
!
      data cday/'01' , '02' , '03' , '04' , '05' , '06' , '07' , '08' , &
          &'09' , '10' , '11' , '12' , '13' , '14' , '15' , '16' ,      &
         & '17' , '18' , '19' , '20' , '21' , '22' , '23' , '24' ,      &
         & '25' , '26' , '27' , '28' , '29' , '30' , '31'/
      data cmonth/'jan' , 'feb' , 'mar' , 'apr' , 'may' , 'jun' ,       &
         & 'jul' , 'aug' , 'sep' , 'oct' , 'nov' , 'dec'/
!
      alatmin = 999999.
      alatmax = -999999.
      alonmin = 999999.
      alonmax = -999999.
      nx = 0
      ny = 0
      centerj = (jxm2)/2.
      centeri = (iym2)/2.

      open (31,file=trim(dirout)//pthsep//ctlname,status='replace')
      write (31,99001) ctlname(1:14)
      write (31,99002)
      if ( ibigend.eq.1 ) then
        write (31,99003)
      else if ( ibigend.eq.0 ) then
        write (31,99004)
      else
      end if
      if ( iotyp.eq.2 ) write (31,'(a)') 'options sequential'
      write (31,99005)
      if ( iproj.eq.'LAMCON' .or. iproj.eq.'ROTMER' ) then
        do j = 2 , jx-1
#ifdef MPP1
          if ( xlat_io(2,j).lt.alatmin ) alatmin = xlat_io(2,j)
          if ( xlat_io(iy-1,j).gt.alatmax ) alatmax = xlat_io(iy-1,j)
#else
          if ( xlat(2,j).lt.alatmin ) alatmin = xlat(2,j)
          if ( xlat(iy-1,j).gt.alatmax ) alatmax = xlat(iy-1,j)
#endif
        end do
        do i = 2 , iy-1
          do j = 2 , jx-1
#ifdef MPP1
            if ( clon.ge.0.0 ) then
              if ( xlong_io(i,j).ge.0.0 ) then
                alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
                alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
              else if ( abs(clon-xlong_io(i,j))                         &
                      & .lt.abs(clon-(xlong_io(i,j)+360.)) ) then
                alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
                alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
              else
                alonmin = amin1(alonmin,sngl(xlong_io(i,j))+360.)
                alonmax = amax1(alonmax,sngl(xlong_io(i,j))+360.)
              end if
            else if ( xlong_io(i,j).lt.0.0 ) then
              alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
              alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
            else if ( abs(clon-xlong_io(i,j))                           &
                    & .lt.abs(clon-(xlong_io(i,j)-360.)) ) then
              alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
              alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
            else
              alonmin = amin1(alonmin,sngl(xlong_io(i,j))-360.)
              alonmax = amax1(alonmax,sngl(xlong_io(i,j))-360.)
            end if
#else
            if ( clon.ge.0.0 ) then
              if ( xlong(i,j).ge.0.0 ) then
                alonmin = amin1(alonmin,sngl(xlong(i,j)))
                alonmax = amax1(alonmax,sngl(xlong(i,j)))
              else if ( abs(clon-xlong(i,j))                            &
                      & .lt.abs(clon-(xlong(i,j)+360.)) ) then
                alonmin = amin1(alonmin,sngl(xlong(i,j)))
                alonmax = amax1(alonmax,sngl(xlong(i,j)))
              else
                alonmin = amin1(alonmin,sngl(xlong(i,j))+360.)
                alonmax = amax1(alonmax,sngl(xlong(i,j))+360.)
              end if
            else if ( xlong(i,j).lt.0.0 ) then
              alonmin = amin1(alonmin,sngl(xlong(i,j)))
              alonmax = amax1(alonmax,sngl(xlong(i,j)))
            else if ( abs(clon-xlong(i,j))                              &
                    & .lt.abs(clon-(xlong(i,j)-360.)) ) then
              alonmin = amin1(alonmin,sngl(xlong(i,j)))
              alonmax = amax1(alonmax,sngl(xlong(i,j)))
            else
              alonmin = amin1(alonmin,sngl(xlong(i,j))-360.)
              alonmax = amax1(alonmax,sngl(xlong(i,j))-360.)
            end if
#endif
          end do
        end do
        rlatinc = dx*0.001/111./2.
        rloninc = dx*0.001/111./2.
        ny = 2 + nint(abs(alatmax-alatmin)/rlatinc)
        nx = 1 + nint(abs((alonmax-alonmin)/rloninc))
      end if
      if ( iotyp.eq.1 ) then
        if ( iproj.eq.'LAMCON' ) then   ! Lambert projection
          write (31,99006) jxm2 , iym2 , clat , clon , centerj ,        &
                         & centeri , truelatl , truelath , clon , dx ,&
                         & dx
          write (31,99007) nx + 2 , alonmin - rloninc , rloninc
          write (31,99008) ny + 2 , alatmin - rlatinc , rlatinc
        else if ( iproj.eq.'POLSTR' ) then
                                        !
        else if ( iproj.eq.'NORMER' ) then
#ifdef MPP1
          write (31,99009) jxm2 , xlong_io(2,2) , xlong_io(2,3)         &
                         & - xlong_io(2,2)
          write (31,99010) iym2
          write (31,99011) (xlat_io(i,2),i=2,iym1)
#else
          write (31,99009) jxm2 , xlong(2,2) , xlong(2,3) - xlong(2,2)
          write (31,99010) iym2
          write (31,99011) (xlat(i,2),i=2,iym1)
#endif
        else if ( iproj.eq.'ROTMER' ) then
          write (*,*) 'Note that rotated Mercartor (ROTMER)' ,          &
                     &' projections are not supported by GrADS.'
          write (*,*) '  Although not exact, the eta.u projection' ,    &
                     &' in GrADS is somewhat similar.'
          write (*,*) ' FERRET, however, does support this projection.'
          write (31,99012) jxm2 , iym2 , plon , plat ,                  &
                         & dx/111000. , dx/111000.*.95238
          write (31,99007) nx + 2 , alonmin - rloninc , rloninc
          write (31,99008) ny + 2 , alatmin - rlatinc , rlatinc
        else
          call fatal(__FILE__,__LINE__,'INVALID MAP PROJECTION')
        end if
      else if ( iotyp.eq.2 ) then
#ifdef MPP1
        write (31,99009) jxm2 , xlong_io(2,2) , xlong_io(2,3)           &
                       & - xlong_io(2,2)
        write (31,99013) iym2 , xlat_io(2,2) , xlat_io(3,2)             &
                       & - xlat_io(2,2)
#else
        write (31,99009) jxm2 , xlong(2,2) , xlong(2,3) - xlong(2,2)
        write (31,99013) iym2 , xlat(2,2) , xlat(3,2) - xlat(2,2)
#endif
      else
      end if
      write (31,99014) kz , ((1013.25-r8pt*10.)*a(k)+r8pt*10.,k=kz,1,-1)
      myear = ldatez/1000000
      month = (ldatez-myear*1000000)/10000
      nday = (ldatez-myear*1000000-month*10000)/100
      nhour = mod(ldatez,100)
      call finddate(nbase,ldatez)
      if ( month.eq.12 ) then
        call finddate(mnend,myear*1000000+1010100)
      else
        call finddate(mnend,myear*1000000+month*10000+10100)
      end if
      call finddate(jbend,idate2)
      if ( ldatez.eq.idate0 ) then
        nnumb = (ibdyfrq/radisp+0.00001)*(min0(jbend,mnend)-nbase) + 1
      else
        nnumb = (ibdyfrq/radisp+0.00001)*(min0(jbend,mnend)-nbase)
      end if
      ifrq = radisp + 0.00001
      if ( ldatez.eq.idate0 ) then
        write (31,99015) nnumb , nhour , cday(nday) , cmonth(month) ,   &
                       & myear , ifrq
      else
        write (31,99015) nnumb , nhour + ifrq , cday(nday) ,            &
                       & cmonth(month) , myear , ifrq
      end if
      if ( iotyp.eq.2 ) write (31,'(a)') 'theader 4'
      write (31,99016) 14
      write (31,99018) 'cld   ' , kz ,                                  &
                      &'cloud fractional cover               '
      write (31,99018) 'clwp  ' , kz ,                                  &
                      &'cloud liquid water path              '
      write (31,99018) 'qrs   ' , kz ,                                  &
                      &'solar heating rate                   '
      write (31,99018) 'qrl   ' , kz ,                                  &
                      &'longwave cooling rate                '
      write (31,99017) 'frsa  ' ,                                       &
                      &'surface absorbed solar flux          '
      write (31,99017) 'frla  ' ,                                       &
                      &'longwave cooling of surface          '
      write (31,99017) 'clrst ' ,                                       &
                      &'clearsky total column abs solar flux '
      write (31,99017) 'clrss ' ,                                       &
                      &'clearsky surface absorbed solar flux '
      write (31,99017) 'clrlt ' ,                                       &
                      &'clearsky net upward LW flux at TOA   '
      write (31,99017) 'clrls ' ,                                       &
                      &'clearsky LW cooling at surface (W/m2)'
      write (31,99017) 'solin ' ,                                       &
                      &'instantaneous incident solar (W/m2)  '
      write (31,99017) 'sabtp ' ,                                       &
                      &'total column absorbed solar flux W/m2'
      write (31,99017) 'firtp ' ,                                       &
                      &'net upward LW flux at TOA (W/m2)     '
      write (31,99017) 'psa   ' ,                                       &
                      &'Surface pressure (hPa)               '
      write (31,99019)
      close (31)
99001 format ('dset ^',a14)
99002 format ('title RegCM normal output variables')
99003 format ('options big_endian')
99004 format ('options little_endian')
99005 format ('undef -1.34')
99006 format ('pdef ',i4,1x,i4,1x,'lcc',7(1x,f7.2),1x,2(f7.0,1x))
99007 format ('xdef ',i4,' linear ',f7.2,1x,f7.4)
99008 format ('ydef ',i4,' linear ',f7.2,1x,f7.4)
99009 format ('xdef ',i3,' linear ',f9.4,' ',f9.4)
99010 format ('ydef ',i3,' levels')
99011 format (10F7.2)
99012 format ('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
99013 format ('ydef ',i3,' linear ',f9.4,' ',f9.4)
99014 format ('zdef ',i2,' levels ',30F7.2)
99015 format ('tdef ',i4,' linear ',i2,'z',a2,a3,i4,' ',i2,'hr')
99016 format ('vars ',i2)
99017 format (a6,'0 99 ',a36)
99018 format (a6,i2,' 0 ',a36)
99019 format ('endvars')
      end subroutine gradsrad

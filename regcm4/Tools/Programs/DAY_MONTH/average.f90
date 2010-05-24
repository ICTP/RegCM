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

! This code treats the daily & monthly mean of all RegCM Input/Output
! (ICBC, ATM, RAD, SRF, SUB and CHE) files, written by Xunqiang Bi, ESP/ICTP
!
      program average
      implicit none
      integer iy,jx,kz,ntr,nsg,ibyte
      logical day_ICBC,day_ATM,day_RAD,day_SRF,day_SUB,day_CHE
      logical mon_ICBC,mon_ATM,mon_RAD,mon_SRF,mon_SUB,mon_CHE
      integer idate0,idate1,idate2
      character*80 Path_Input,Path_Output
      character*10 DomainName

      namelist /shareparam/iy,jx,kz,nsg,ntr,ibyte &
                          ,Path_Input,DomainName,Path_Output
      namelist /dateparam/ idate0,idate1,idate2
      namelist /averageparam/ &
                day_ICBC,day_ATM,day_RAD,day_SRF,day_SUB,day_CHE &
               ,mon_ICBC,mon_ATM,mon_RAD,mon_SRF,mon_SUB,mon_CHE

      read(*,shareparam)
      read(*,dateparam)
      read(*,averageparam)

      if(day_ICBC) call day2('ICBC',iy,jx,kz,ibyte, Path_Input &
                            ,DomainName,idate0,idate1,idate2)
      if(day_ATM)call day('ATM',iy,jx,kz,ntr,ibyte, Path_Output &
                          ,idate0,idate1,idate2)
      if(day_RAD)call day('RAD',iy,jx,kz,ntr,ibyte, Path_Output &
                          ,idate0,idate1,idate2)
      if(day_SRF)call day('SRF',iy,jx,kz,ntr,ibyte, Path_Output &
                          ,idate0,idate1,idate2)
      if(nsg.gt.1.and.day_SUB) &
                call day3('SUB',iy,jx,kz,nsg,ibyte, Path_Input &
                         ,DomainName,Path_Output,idate0,idate1,idate2)
      if(day_CHE)call day('CHE',iy,jx,kz,ntr,ibyte, Path_Output &
                         ,idate0,idate1,idate2)

      if(mon_ICBC) call mon2('ICBC',iy,jx,kz,ibyte, Path_Input &
                            ,DomainName,idate0,idate1,idate2)
      if(mon_ATM)call mon('ATM',iy,jx,kz,ntr,ibyte, Path_Output &
                         ,idate0,idate1,idate2)
      if(mon_RAD)call mon('RAD',iy,jx,kz,ntr,ibyte, Path_Output &
                         ,idate0,idate1,idate2)
      if(mon_SRF)call mon('SRF',iy,jx,kz,ntr,ibyte, Path_Output &
                         ,idate0,idate1,idate2)
      if(nsg.gt.1.and.mon_SUB) &
                call mon3('SUB',iy,jx,kz,nsg,ibyte, Path_Input &
                         ,DomainName,Path_Output,idate0,idate1,idate2)
      if(mon_CHE)call mon('CHE',iy,jx,kz,ntr,ibyte, Path_Output &
                         ,idate0,idate1,idate2)

      stop
      end

      subroutine day(filename,iy,jx,kz,ntr,ibyte, Path_Output &
                    ,idate0,idate1,idate2)
      implicit none
      character*3 filename
      character*80 Path_Output
      integer iy,jx,kz,ntr,ibyte,idate0,idate1,idate2
      integer iiy,jjx,kkz
      integer mdate0,ibltyp,icup,ipptls,iboudy
      real*4  truelatL,truelatH
      real*4  dxsp,ptsp,clat,clon,plat,plon
      real*4  dto,dtb,dtr,dtc,dt
      integer iotyp
      character*6 iproj
      real*4, allocatable ::  sigma(:)
      character*4 :: chy(1941:2100)
      data chy/'1941','1942','1943','1944','1945','1946','1947','1948',&
        '1949','1950','1951','1952','1953','1954','1955','1956','1957',&
        '1958','1959','1960','1961','1962','1963','1964','1965','1966',&
        '1967','1968','1969','1970','1971','1972','1973','1974','1975',&
        '1976','1977','1978','1979','1980','1981','1982','1983','1984',&
        '1985','1986','1987','1988','1989','1990','1991','1992','1993',&
        '1994','1995','1996','1997','1998','1999','2000','2001','2002',&
        '2003','2004','2005','2006','2007','2008','2009','2010','2011',&
        '2012','2013','2014','2015','2016','2017','2018','2019','2020',&
        '2021','2022','2023','2024','2025','2026','2027','2028','2029',&
        '2030','2031','2032','2033','2034','2035','2036','2037','2038',&
        '2039','2040','2041','2042','2043','2044','2045','2046','2047',&
        '2048','2049','2050','2051','2052','2053','2054','2055','2056',&
        '2057','2058','2059','2060','2061','2062','2063','2064','2065',&
        '2066','2067','2068','2069','2070','2071','2072','2073','2074',&
        '2075','2076','2077','2078','2079','2080','2081','2082','2083',&
        '2084','2085','2086','2087','2088','2089','2090','2091','2092',&
        '2093','2094','2095','2096','2097','2098','2099','2100'/
      character*2 cday(31)
      data cday/'01','02','03','04','05','06','07','08','09','10', &
                '11','12','13','14','15','16','17','18','19','20', &
                '21','22','23','24','25','26','27','28','29','30','31'/
      character*2 chm(12)
      data chm/'01','02','03','04','05','06','07','08','09','10', &
               '11','12'/
      character*3 chmc(12)
      data chmc/'jan','feb','mar','apr','may','jun'  &
               ,'jul','aug','sep','oct','nov','dec'/
      character*14 filein
      character*12 fileout
      integer ntype,nfile,nyear,month,n_slice,mrec,nrec
      logical there
      real*4  alatmin,alatmax,alonmin,alonmax,rlatinc,rloninc
      real*4  centerj,centeri
      integer ny,nx
      integer i,j,k,l,m,n,itr
      integer number,nvar,mvar,n_month,nrecord,nday,mday
      integer idatex
      real*4, allocatable ::  c(:,:,:),b(:,:),xlat(:,:),xlon(:,:)

      allocate(sigma(kz+1))
      allocate(b(jx-2,iy-2))
      allocate(xlat(jx-2,iy-2))
      allocate(xlon(jx-2,iy-2))

      open(10,file=trim(Path_Output)//'OUT_HEAD',form='unformatted' &
             ,recl=(jx-2)*(iy-2)*ibyte,access='direct')
      read(10,rec=1) mdate0,ibltyp,icup,ipptls,iboudy  &
                    ,iiy,jjx,kkz,(sigma(k),k=1,kz+1)   &
                    ,dxsp,ptsp,clat,clon,plat,plon     &
                    ,iproj,dto,dtb,dtr,dtc,iotyp,truelatL,truelatH
      if(iiy.ne.iy.or.jjx.ne.jx.or.kkz.ne.kz) then
         write(*,*) 'iy,jx,kz in parameter = ',iy,jx,kz
         write(*,*) 'iy,jx,kz in OUT_HEAD ',iiy,jjx,kkz
         write(*,*) 'They are not consistent'
         stop
      endif
      read(10,rec=6) xlat
      read(10,rec=7) xlon
      close(10)

      if(filename.eq.'ATM') then
         dt = dto
         nvar = kz*6+5
      else if(filename.eq.'SRF') then
         dt = dtb
         nvar = 27
      else if(filename.eq.'RAD') then
         dt = dtr
         nvar = kz*4+9
      else if(filename.eq.'CHE') then
         dt = dtc
         nvar = ntr*kz+kz*3+ntr*7+3
!        nvar = ntr*kz+kz*3+ntr*7+2      ! for RegCM3, one record less
      else
         write(*,*) 'filename is not correct'
         stop
      endif
      allocate(c(jx-2,iy-2,nvar))
      mvar = nvar
      if(filename.eq.'SRF') mvar = mvar -6

      if(idate1.gt.1940010100) then        ! original file
         n_month = (idate2/1000000-idate1/1000000)*12 &
                 + (mod(idate2/10000,100)-mod(idate1/10000,100))
         if(mod(idate2,10000).gt.0100) n_month = n_month+1
         ntype = 0
      else
         write(*,*) 'date should be in 10 digit'
         stop
      endif

      do nfile=1,n_month
         nyear = idate1/1000000
         month = idate1/10000-nyear*100 +nfile-1
         nyear = nyear + (month-1)/12
         month = mod(month,12)
         if(month.eq.0) month = 12

         if(nfile.eq.1.or.month.eq.1) then
         write(*,*) 'Calculate the daily mean of ',filename,nyear,month
         else
         write(*,*) '                            ',filename,nyear,month
         endif

         nrec = 0
         if(month.eq.1.or.month.eq.3.or.month.eq.5.or.  &
            month.eq.7.or.month.eq.8.or.month.eq.10.or. &
            month.eq.12) then
            nrecord = 31
         else if(month.eq.4.or.month.eq.6.or.month.eq.9.or. &
                 month.eq.11) then
            nrecord = 30
         else
            nrecord = 28
            if(mod(nyear,4).eq.0) nrecord = 29
            if(mod(nyear,100).eq.0) nrecord = 28
            if(mod(nyear,400).eq.0) nrecord = 29
         endif
         if(nfile.eq.n_month.and.mod(idate2/100-1,100).ne.0) &
            nrecord = min(nrecord,mod(idate2/100-1,100))
         if(nfile.eq.1.and.idate0.eq.idate1) nrec = nvar
         n_slice = 24/nint(dt)

         if(filename.eq.'ATM') then
            filein = 'ATM.'//chy(nyear)//chm(month)//'0100'
            fileout= 'ATM.'//chy(nyear)//chm(month)//'01'
         else if(filename.eq.'SRF') then
            filein = 'SRF.'//chy(nyear)//chm(month)//'0100'
            fileout= 'SRF.'//chy(nyear)//chm(month)//'01'
         else if(filename.eq.'RAD') then
            filein = 'RAD.'//chy(nyear)//chm(month)//'0100'
            fileout= 'RAD.'//chy(nyear)//chm(month)//'01'
         else if(filename.eq.'CHE') then
            filein = 'CHE.'//chy(nyear)//chm(month)//'0100'
            fileout= 'CHE.'//chy(nyear)//chm(month)//'01'
         endif

         inquire(file=trim(Path_Output)//fileout//'.ctl',exist=there)
         if(there) then
           open(31,file=trim(Path_Output)//fileout//'.ctl'  &
                  ,status='replace')
         else
           open(31,file=trim(Path_Output)//fileout//'.ctl',status='new')
         endif
         write(31,10) fileout
  10     format('dset ^',A12)
         write(31,20)
  20     format('title RegCM daily output variables')
         write(31,30)
  30     format('options big_endian')
         write(31,50)
  50     format('undef -1.e34')
         if(iproj.eq.'LAMCON'.or.iproj.eq.'ROTMER') then
            alatmin= 999999.
            alatmax=-999999.
            do j=1,jx-2
               if(xlat(j,1).lt.alatmin) alatmin=xlat(j,1)
               if(xlat(j,iy-2).gt.alatmax) alatmax=xlat(j,iy-2)
            enddo
            alonmin= 999999.
            alonmax=-999999.
            do i=1,iy-2
            do j=1,jx-2
               if(clon.ge.0.0) then
                  if(xlon(j,i).ge.0.0) then
                     alonmin = amin1(alonmin,xlon(j,i))
                     alonmax = amax1(alonmax,xlon(j,i))
                  else
                     if(abs(clon-xlon(j,i)).lt. &
                        abs(clon-(xlon(j,i)+360.))) then
                        alonmin = amin1(alonmin,xlon(j,i))
                        alonmax = amax1(alonmax,xlon(j,i))
                     else
                        alonmin = amin1(alonmin,xlon(j,i)+360.)
                        alonmax = amax1(alonmax,xlon(j,i)+360.)
                     endif
                  endif
               else
                  if(xlon(j,i).lt.0.0) then
                     alonmin = amin1(alonmin,xlon(j,i))
                     alonmax = amax1(alonmax,xlon(j,i))
                  else
                     if(abs(clon-xlon(j,i).lt. &
                        abs(clon-xlon(j,i)-360.))) then
                        alonmin = amin1(alonmin,xlon(j,i))
                        alonmax = amax1(alonmax,xlon(j,i))
                     else
                        alonmin = amin1(alonmin,xlon(j,i)-360.)
                        alonmax = amax1(alonmax,xlon(j,i)-360.)
                     endif
                  endif
               endif
            enddo
            enddo
            rlatinc=dxsp/111./2.
            rloninc=dxsp/111./2.
            ny=2+nint(abs(alatmax-alatmin)/rlatinc)
            nx=1+nint(abs((alonmax-alonmin)/rloninc))
            centerj=(jx-2)/2.
            centeri=(iy-2)/2.
         endif

         if(iproj.eq.'LAMCON') then        ! Lambert projection
            write(31,100) jx-2,iy-2,clat,clon,centerj,centeri, &
                          truelatL,truelatH,clon,dxsp*1000.,dxsp*1000.
 100  format('pdef ',i4,1x,i4,1x,'lccr',7(1x,f7.2),1x,2(f7.0,1x))
            write(31,110) nx+2,alonmin-rloninc,rloninc
 110  format('xdef ',i4,' linear ',f7.2,1x,f7.4)
            write(31,120) ny+2,alatmin-rlatinc,rlatinc
 120  format('ydef ',i4,' linear ',f7.2,1x,f7.4)
         elseif(iproj.eq.'POLSTR') then    !
         elseif(iproj.eq.'NORMER') then
            write(31,200)  jx-2,xlon(1,1),xlon(2,1)-xlon(1,1)
 200  format('xdef ',I3,' linear ',f9.4,' ',f9.4)
            write(31,210) iy-2
 210  format('ydef ',I3,' levels')
            write(31,220) (xlat(1,i),i=1,iy-2)
 220  format(10f7.2)
         elseif(iproj.eq.'ROTMER') then
            write(*,*) 'Note that rotated Mercartor (ROTMER)' &
                   ,' projections are not supported by GrADS.'
            write(*,*) '  Although not exact, the eta.u projection' &
                   ,' in GrADS is somewhat similar.'
            write(*,*) ' FERRET, however, does support this projection.'
            write(31,230) jx-2,iy-2,plon,plat,dxsp/111. &
                                          ,dxsp/111.*.95238
 230  format('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
            write(31,110) nx+2,alonmin-rloninc,rloninc
            write(31,120) ny+2,alatmin-rlatinc,rlatinc
         else
            write(*,*) 'Are you sure your map projection is right ?'
            stop
         endif
         if(filename.eq.'SRF') then
      write(31,300) (1013.25-ptsp*10.)*(sigma(1)+sigma(2))*.5+ptsp*10.
         else
      write(31,318) kz,((1013.25-ptsp*10.)*(sigma(k)+sigma(k+1))*.5 &
                       +ptsp*10.,k=1,kz)
         endif
 300  format('zdef 1',' levels ',f7.2)
 318  format('zdef ',I2,' levels ',30f7.2)
         if(nfile.eq.1) then
            nday=mod(idate0,10000)/100
            nrecord = nrecord+1-nday 
            write(31,400)nrecord,cday(nday),chmc(month),nyear
         else if(nfile.eq.n_month) then
            nday=mod(idate2,10000)/100
            if(nday.ne.1) nrecord = nday
            write(31,400)nrecord,cday(1),chmc(month),nyear
         else
            write(31,400)nrecord,cday(1),chmc(month),nyear
         endif
 400  format('tdef ',I4,' linear ',A2,A3,I4,' ','1dy')
      if(filename.eq.'ATM') then
         write(31,500) 6+5
      else if(filename.eq.'RAD') then
         write(31,500) 4+9
      else if(filename.eq.'SRF') then
         write(31,500) nvar
      else if(filename.eq.'CHE') then
         write(31,500) ntr*8+5
      endif
 500  format('vars ',I3)
 600  format(A8,'0 99 ',A36)
 611  format(A8,'0 33,105 ',A36)
 612  format(A8,'0 34,105 ',A36)
      if(filename.eq.'SRF') then
       if(iproj.eq.'LAMCON') then        ! Lambert projection
        write(31,611) 'u10m    ','westerly  wind at 10m (m/s)          '
        write(31,612) 'v10m    ','southerly wind at 10m (m/s)          '
       else
        write(31,600) 'u10m    ','westerly  wind at 10m (m/s)          '
        write(31,600) 'v10m    ','southerly wind at 10m (m/s)          '
       endif
        write(31,600) 'uvdrag  ','surface drag stress                  '
        write(31,600) 'tg      ','ground temperature (degree)          '
        write(31,600) 'tlef    ','temperature of foliage               '
        write(31,600) 't2m     ','air temperature at 2m (K)            '
        write(31,600) 'q2m     ','water vapor mixing ratio at 2m(kg/kg)'
        write(31,600) 'ssw     ','upper layer soil water               '
        write(31,600) 'rsw     ','root zone soil water                 '
        write(31,600) 'tpr     ','total precipitation (mm/day)         '
        write(31,600) 'evp     ','evapotranspiration (mm/day)          '
        write(31,600) 'runoff  ','surface runoff (mm/day)              '
        write(31,600) 'scv     ','total snow amount                    '
        write(31,600) 'sena    ','sensible heat flux (W/m2)            '
        write(31,600) 'flw     ','net infrared energy flux (W/m2)      '
        write(31,600) 'fsw     ','net absorbed solar energy flux (W/m2)'
        write(31,600) 'flwd    ','downward infrared energy flux (W/m2) '
        write(31,600) 'sina    ','incident solar energy flux (W/m2)    '
        write(31,600) 'prcv    ','convective precipitation (mm/day)    '
        write(31,600) 'psb     ','surface pressure (hPa)               '
        write(31,600) 'zpbl    ','PBL layer height                     '
        write(31,600) 'tgmax   ','maximum ground temperature (K)       '
        write(31,600) 'tgmin   ','minimum ground temperature (K)       '
        write(31,600) 't2max   ','maximum 2m air temperature (K)       '
        write(31,600) 't2min   ','minimum 2m air temperature (K)       '
        write(31,600) 'w10max  ','maximum 10m wind speed (m/s)         '
        write(31,600) 'ps_min  ','minimum surface pressure (hPa)       '
      else if(filename.eq.'ATM') then
 660  format(A8,I2,' 0 ',A36)
 661  format(A8,I2,' 33,100 ',A36)
 662  format(A8,I2,' 34,100 ',A36)
        if(iproj.eq.'LAMCON') then
           write(31,661) 'u       ',kz,'westerly wind (m/s)         '
           write(31,662) 'v       ',kz,'southerly wind (m/s)        '
        else
           write(31,660) 'u       ',kz,'westerly wind (m/s)         '
           write(31,660) 'v       ',kz,'southerly wind (m/s)        '
        endif
        write(31,660) 'w       ',kz,'omega (hPa/s)   p-velocity  '
        write(31,660) 't       ',kz,'air temperature (degree)    '
        write(31,660) 'qv      ',kz,'water vapor mixing ratio    '
        write(31,660) 'qc      ',kz,'cloud water mixing ratio    '
        write(31,600) 'psa     ',   'surface pressure (hPa)      '
        write(31,600) 'tpr     ',   'total precipitation(mm/day) '
        write(31,600) 'tgb     ',   'lower groud temp. in BATS   '
        write(31,600) 'swt     ',   'total soil water in mm H2O  '
        write(31,600) 'rno     ',   'accumulated infiltration    '
      else if(filename.eq.'RAD') then
        write(31,660)'cld   ',kz,'cloud fractional cover               '
        write(31,660)'clwp  ',kz,'cloud liquid water path              '
        write(31,660)'qrs   ',kz,'solar heating rate                   '
        write(31,660)'qrl   ',kz,'longwave cooling rate                '
        write(31,600)'frsa  ',   'surface absorbed solar flux          '
        write(31,600)'frla  ',   'longwave cooling of surface          '
        write(31,600)'clrst ',   'clearsky total column abs solar flux '
        write(31,600)'clrss ',   'clearsky surface absorbed solar flux '
        write(31,600)'clrlt ',   'clearsky net upward LW flux at TOA   '
        write(31,600)'clrls ',   'clearsky LW cooling at surface (W/m2)'
        write(31,600)'solin ',   'instantaneous incident solar (W/m2)  '
        write(31,600)'sabtp ',   'total column absorbed solar flux W/m2'
        write(31,600)'firtp ',   'net upward LW flux at TOA (W/m2)     '
      else if(filename.eq.'CHE') then
       do itr=1,ntr
         if(itr.lt.10) then
           write(31,650) 'trac',itr,kz, 'tracer mix. rat  (Kg/Kg)  '
         else
           write(31,651) 'trac',itr,kz, 'tracer mix. rat  (Kg/Kg)  '
         endif
       enddo
 650  format(A4,I1,' ',I2,' 0 ',A36)
 651  format(A4,I2,' ',I2,' 0 ',A36)
 655  format(A8,I1,' 0 99 ',A36)
 656  format(A8,I2,' 0 99 ',A36)
       write(31,650) 'aext',8,kz, 'aer mix. ext. coef      '
       write(31,650) 'assa',8,kz, 'aer mix. sin. scat. alb '
       write(31,650) 'agfu',8,kz, 'aer mix. ass. par       '
       do itr=1,ntr
        if(itr.lt.10) then
          write(31,655) 'colb__tr',itr,'columnburden inst(mg/m2)'
          write(31,655) 'wdlsc_tr',itr,'wet dep lgscale(mg/m2/d)'
          write(31,655) 'wdcvc_tr',itr,'wet dep convect(mg/m2/d)'
          write(31,655) 'sdrdp_tr',itr,'surf dry depos.(mg/m2/d)'
          write(31,655) 'xgasc_tr',itr,'chem gas conv. (mg/m2/d)'
          write(31,655) 'xaquc_tr',itr,'chem aqu conv. (mg/m2/d)'
          write(31,655) 'emiss_tr',itr,'surf emission  (mg/m2/d)'
        else
          write(31,656) 'colb__tr',itr,'columnburden inst(mg/m2)'
          write(31,656) 'wdlsc_tr',itr,'wet dep lgscale(mg/m2/d)'
          write(31,656) 'wdcvc_tr',itr,'wet dep convect(mg/m2/d)'
          write(31,656) 'sdrdp_tr',itr,'surf dry depos.(mg/m2/d)'
          write(31,656) 'xgasc_tr',itr,'chem gas conv. (mg/m2/d)'
          write(31,656) 'xaquc_tr',itr,'chem aqu conv. (mg/m2/d)'
          write(31,656) 'emiss_tr',itr,'surf emission  (mg/m2/d)'
        endif
       enddo
       write(31,600) 'acstoarf',' TOArad forcing av.(W/m2)            '
       write(31,600) 'acstsrrf',' SRFrad forcing av.(W/m2)            '
       write(31,600) 'psb     ','surface pressure (hPa)               '
      endif
      write(31,700)
 700  format('endvars')
      close(31)

         if(iotyp.eq.1) then
            open(10,file=trim(Path_Output)//filein,form='unformatted' &
                   ,recl=(iy-2)*(jx-2)*ibyte,access='direct')
         else if(iotyp.eq.2) then
            open(10,file=trim(Path_Output)//filein,form='unformatted')
         endif
         open(20,file=trim(Path_Output)//fileout,form='unformatted' &
                ,recl=(iy-2)*(jx-2)*ibyte,access='direct')
         mrec = 0
         do nday=1,nrecord
            do n=1,mvar
               do i=1,iy-2
               do j=1,jx-2
                  c(j,i,n) = 0.0
               enddo
               enddo
            enddo
            if(iotyp.eq.2) read(10) idatex
            if(filename.eq.'SRF') then
               do i=1,iy-2
               do j=1,jx-2
                  c(j,i,nvar-5)= -1.e20
                  c(j,i,nvar-4)=  1.e20
                  c(j,i,nvar-3)= -1.e20
                  c(j,i,nvar-2)=  1.e20
                  c(j,i,nvar-1)= -1.e20
                  c(j,i,nvar)  =  1.e20
               enddo
               enddo
            endif
            do mday = 1,n_slice
               do n=1,mvar
                  if(iotyp.eq.1) then
                     nrec=nrec+1
                     read(10,rec=nrec) b
                  else if(iotyp.eq.2) then
                     read(10) b
                  endif
                  do i=1,iy-2
                  do j=1,jx-2
                     c(j,i,n) = c(j,i,n)+b(j,i)
                  enddo
                  enddo
               enddo
               if(filename.eq.'SRF') then
                  if(iotyp.eq.1) then
                     nrec=nrec+1
                     read(10,rec=nrec) b
                  else if(iotyp.eq.2) then
                     read(10) b
                  endif
                  do i=1,iy-2
                  do j=1,jx-2
                     if(c(j,i,nvar-5).lt.b(j,i)) &
                        c(j,i,nvar-5) =  b(j,i)
                  enddo
                  enddo
                  if(iotyp.eq.1) then
                     nrec=nrec+1
                     read(10,rec=nrec) b
                  else if(iotyp.eq.2) then
                     read(10) b
                  endif
                  do i=1,iy-2
                  do j=1,jx-2
                     if(c(j,i,nvar-4).gt.b(j,i)) &
                        c(j,i,nvar-4) =  b(j,i)
                  enddo
                  enddo
                  if(iotyp.eq.1) then
                     nrec=nrec+1
                     read(10,rec=nrec) b
                  else if(iotyp.eq.2) then
                     read(10) b
                  endif
                  do i=1,iy-2
                  do j=1,jx-2
                     if(c(j,i,nvar-3).lt.b(j,i)) &
                        c(j,i,nvar-3) =  b(j,i)
                  enddo
                  enddo
                  if(iotyp.eq.1) then
                     nrec=nrec+1
                     read(10,rec=nrec) b
                  else if(iotyp.eq.2) then
                     read(10) b
                  endif
                  do i=1,iy-2
                  do j=1,jx-2
                     if(c(j,i,nvar-2).gt.b(j,i)) &
                        c(j,i,nvar-2) =  b(j,i)
                  enddo
                  enddo
                  if(iotyp.eq.1) then
                     nrec=nrec+1
                     read(10,rec=nrec) b
                  else if(iotyp.eq.2) then
                     read(10) b
                  endif
                  do i=1,iy-2
                  do j=1,jx-2
                     if(c(j,i,nvar-1).lt.b(j,i)) &
                        c(j,i,nvar-1) =  b(j,i)
                  enddo
                  enddo
                  if(iotyp.eq.1) then
                     nrec=nrec+1
                     read(10,rec=nrec) b
                  else if(iotyp.eq.2) then
                     read(10) b
                  endif
                  do i=1,iy-2
                  do j=1,jx-2
                     if(c(j,i,nvar).gt.b(j,i)) &
                        c(j,i,nvar) =  b(j,i)
                  enddo
                  enddo
               endif
            enddo
            do n=1,mvar
               do i=1,iy-2
               do j=1,jx-2
                  c(j,i,n) = c(j,i,n)/float(n_slice)
               enddo
               enddo
               mrec=mrec+1
               write(20,rec=mrec) ((c(j,i,n),j=1,jx-2),i=1,iy-2)
            enddo
            if(filename.eq.'SRF') then
               do n=nvar-5,nvar
                  mrec=mrec+1
                  write(20,rec=mrec) ((c(j,i,n),j=1,jx-2),i=1,iy-2)
               enddo
            endif
         enddo
         close(10)
         close(20)
      enddo
      deallocate(sigma)
      deallocate(b)
      deallocate(xlat)
      deallocate(xlon)
      deallocate(c)
      return
      end

      subroutine day2(filename,iy,jx,kz,ibyte, Path_Input,DomainName &
                     ,idate0,idate1,idate2)
      implicit none
      character*4 filename
      character*80 Path_Input
      character*10 DomainName
      integer iy,jx,kz,ntr,ibyte,idate0,idate1,idate2
      integer iiy,jjx,kkz
      real*4  truelatL,truelatH
      real*4  dsinm,ptop,clat,clon,plat,plon,GRDFAC
      integer igrads,ibigend
      integer iotyp
      character*6 iproj
      real*4, allocatable ::  sigma(:)
      character*4 :: chy(1941:2100)
      data chy/'1941','1942','1943','1944','1945','1946','1947','1948',&
        '1949','1950','1951','1952','1953','1954','1955','1956','1957',&
        '1958','1959','1960','1961','1962','1963','1964','1965','1966',&
        '1967','1968','1969','1970','1971','1972','1973','1974','1975',&
        '1976','1977','1978','1979','1980','1981','1982','1983','1984',&
        '1985','1986','1987','1988','1989','1990','1991','1992','1993',&
        '1994','1995','1996','1997','1998','1999','2000','2001','2002',&
        '2003','2004','2005','2006','2007','2008','2009','2010','2011',&
        '2012','2013','2014','2015','2016','2017','2018','2019','2020',&
        '2021','2022','2023','2024','2025','2026','2027','2028','2029',&
        '2030','2031','2032','2033','2034','2035','2036','2037','2038',&
        '2039','2040','2041','2042','2043','2044','2045','2046','2047',&
        '2048','2049','2050','2051','2052','2053','2054','2055','2056',&
        '2057','2058','2059','2060','2061','2062','2063','2064','2065',&
        '2066','2067','2068','2069','2070','2071','2072','2073','2074',&
        '2075','2076','2077','2078','2079','2080','2081','2082','2083',&
        '2084','2085','2086','2087','2088','2089','2090','2091','2092',&
        '2093','2094','2095','2096','2097','2098','2099','2100'/
      character*2 cday(31)
      data cday/'01','02','03','04','05','06','07','08','09','10', &
                '11','12','13','14','15','16','17','18','19','20', &
                '21','22','23','24','25','26','27','28','29','30','31'/
      character*2 chm(12)
      data chm/'01','02','03','04','05','06','07','08','09','10', &
               '11','12'/
      character*3 chmc(12)
      data chmc/'jan','feb','mar','apr','may','jun'  &
               ,'jul','aug','sep','oct','nov','dec'/
      character*15 filein
      character*13 fileout
      integer ntype,nfile,nyear,month,n_slice,mrec,nrec
      logical there
      real*4  alatmin,alatmax,alonmin,alonmax,rlatinc,rloninc
      real*4  centerj,centeri
      integer ny,nx
      integer i,j,k,l,m,n,itr
      integer number,nvar,n_month,nrecord,nday,mday
      integer idatex
      real*4, allocatable ::  c(:,:,:),b(:,:),xlat(:,:),xlon(:,:)
      real*4, allocatable ::  a(:,:)

      allocate(sigma(kz+1))
      allocate(b(jx,iy))

      allocate(xlat(jx-2,iy-2))
      allocate(xlon(jx-2,iy-2))

      open(10,file=trim(Path_Input)//trim(DomainName)//'.INFO'  &
             ,form='unformatted',recl=jx*iy*ibyte,access='direct')
      read(10,rec=1) iiy,jjx,kkz,dsinm,clat,clon,plat,plon,GRDFAC  &
                    ,iproj,(sigma(k),k=1,kz+1),ptop,igrads,ibigend &
                    ,truelatL,truelatH
      if(iiy.ne.iy.or.jjx.ne.jx.or.kkz.ne.kz) then
         write(*,*) 'iy,jx,kz in parameter = ',iy,jx,kz
         write(*,*) 'iy,jx,kz in DOMAIN.INFO ',iiy,jjx,kkz
         write(*,*) 'They are not consistent'
         stop
      endif
      read(10,rec=5) b
      do i=1,iy-2
      do j=1,jx-2
         xlat(j,i) = b(j+1,i+1)
      enddo
      enddo
      read(10,rec=6) b
      do i=1,iy-2
      do j=1,jx-2
         xlon(j,i) = b(j+1,i+1)
      enddo
      enddo
      close(10)

      if(filename.eq.'ICBC') then
         nvar = kz*4+2
      else
         write(*,*) 'filename is not correct'
         stop
      endif
      allocate(c(jx-2,iy-2,nvar))

      if(idate1.gt.1940010100) then        ! original file
         n_month = (idate2/1000000-idate1/1000000)*12 &
                 + (mod(idate2/10000,100)-mod(idate1/10000,100))
         if(mod(idate2,10000).gt.0100) n_month = n_month+1
         ntype = 0
      else
         write(*,*) 'date should be in 10 digit'
         stop
      endif

      do nfile=1,n_month
         nyear = idate1/1000000
         month = idate1/10000-nyear*100 +nfile-1
         nyear = nyear + (month-1)/12
         month = mod(month,12)
         if(month.eq.0) month = 12

         if(nfile.eq.1.or.month.eq.1) then
         write(*,*) 'Calculate the daily mean of ',filename,nyear,month
         else
         write(*,*) '                            ',filename,nyear,month
         endif

         if(month.eq.1.or.month.eq.3.or.month.eq.5.or.  &
            month.eq.7.or.month.eq.8.or.month.eq.10.or. &
            month.eq.12) then
            nrecord = 31
         else if(month.eq.4.or.month.eq.6.or.month.eq.9.or. &
                 month.eq.11) then
            nrecord = 30
         else
            nrecord = 28
            if(mod(nyear,4).eq.0) nrecord = 29
            if(mod(nyear,100).eq.0) nrecord = 28
            if(mod(nyear,400).eq.0) nrecord = 29
         endif
         if(nfile.eq.n_month.and.mod(idate2/100-1,100).ne.0) &
            nrecord = min(nrecord,mod(idate2/100-1,100))
         nrec = nvar+1
         n_slice = 24/6

         if(filename.eq.'ICBC') then
            filein = '_ICBC'//chy(nyear)//chm(month)//'0100'
            fileout= '_ICBC'//chy(nyear)//chm(month)//'01'
         endif

         inquire(file= &
        trim(Path_Input)//trim(DomainName)//fileout//'.ctl',exist=there)
         if(there) then
            open(31,file= &
                 trim(Path_Input)//trim(DomainName)//fileout//'.ctl' &
                                                      ,status='replace')
         else
            open(31,file= &
                 trim(Path_Input)//trim(DomainName)//fileout//'.ctl' &
                                                      ,status='new')
         endif
         write(31,10) '^'//trim(DomainName)//fileout
  10     format('dset ',A24)
         write(31,20)
  20     format('title RegCM daily input variables')
         write(31,30)
  30     format('options big_endian')
         write(31,50)
  50     format('undef -1.e34')
         if(iproj.eq.'LAMCON'.or.iproj.eq.'ROTMER') then
            alatmin= 999999.
            alatmax=-999999.
            do j=1,jx-2
               if(xlat(j,1).lt.alatmin) alatmin=xlat(j,1)
               if(xlat(j,iy-2).gt.alatmax) alatmax=xlat(j,iy-2)
            enddo
            alonmin= 999999.
            alonmax=-999999.
            do i=1,iy-2
            do j=1,jx-2
               if(clon.ge.0.0) then
                  if(xlon(j,i).ge.0.0) then
                     alonmin = amin1(alonmin,xlon(j,i))
                     alonmax = amax1(alonmax,xlon(j,i))
                  else
                     if(abs(clon-xlon(j,i)).lt. &
                        abs(clon-(xlon(j,i)+360.))) then
                        alonmin = amin1(alonmin,xlon(j,i))
                        alonmax = amax1(alonmax,xlon(j,i))
                     else
                        alonmin = amin1(alonmin,xlon(j,i)+360.)
                        alonmax = amax1(alonmax,xlon(j,i)+360.)
                     endif
                  endif
               else
                  if(xlon(j,i).lt.0.0) then
                     alonmin = amin1(alonmin,xlon(j,i))
                     alonmax = amax1(alonmax,xlon(j,i))
                  else
                     if(abs(clon-xlon(j,i).lt. &
                        abs(clon-xlon(j,i)-360.))) then
                        alonmin = amin1(alonmin,xlon(j,i))
                        alonmax = amax1(alonmax,xlon(j,i))
                     else
                        alonmin = amin1(alonmin,xlon(j,i)-360.)
                        alonmax = amax1(alonmax,xlon(j,i)-360.)
                     endif
                  endif
               endif
            enddo
            enddo
            rlatinc=dsinm*0.001/111./2.
            rloninc=dsinm*0.001/111./2.
            ny=2+nint(abs(alatmax-alatmin)/rlatinc)
            nx=1+nint(abs((alonmax-alonmin)/rloninc))
            centerj=(jx-2)/2.
            centeri=(iy-2)/2.
         endif

         if(iproj.eq.'LAMCON') then        ! Lambert projection
            write(31,100) jx-2,iy-2,clat,clon,centerj,centeri, &
                          truelatL,truelatH,clon,dsinm,dsinm
 100  format('pdef ',i4,1x,i4,1x,'lccr',7(1x,f7.2),1x,2(f7.0,1x))
            write(31,110) nx+2,alonmin-rloninc,rloninc
 110  format('xdef ',i4,' linear ',f7.2,1x,f7.4)
            write(31,120) ny+2,alatmin-rlatinc,rlatinc
 120  format('ydef ',i4,' linear ',f7.2,1x,f7.4)
         elseif(iproj.eq.'POLSTR') then    !
         elseif(iproj.eq.'NORMER') then
            write(31,200)  jx-2,xlon(1,1),xlon(2,1)-xlon(1,1)
 200  format('xdef ',I3,' linear ',f9.4,' ',f9.4)
            write(31,210) iy-2
 210  format('ydef ',I3,' levels')
            write(31,220) (xlat(1,i),i=1,iy-2)
 220  format(10f7.2)
         elseif(iproj.eq.'ROTMER') then
            write(*,*) 'Note that rotated Mercartor (ROTMER)' &
                   ,' projections are not supported by GrADS.'
            write(*,*) '  Although not exact, the eta.u projection' &
                   ,' in GrADS is somewhat similar.'
            write(*,*) ' FERRET, however, does support this projection.'
            write(31,230) jx-2,iy-2,plon,plat,dsinm/111000. &
                                          ,dsinm/111000.*.95238
 230  format('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
            write(31,110) nx+2,alonmin-rloninc,rloninc
            write(31,120) ny+2,alatmin-rlatinc,rlatinc
         else
            write(*,*) 'Are you sure your map projection is right ?'
            stop
         endif
      write(31,318) kz,((1013.25-ptop*10.)*(sigma(k)+sigma(k+1))*.5 &
                       +ptop*10.,k=kz,1,-1)
 318  format('zdef ',I2,' levels ',30f7.2)
         if(nfile.eq.1) then
            nday=mod(idate0,10000)/100
            nrecord = nrecord+1-nday 
            write(31,400)nrecord,cday(nday),chmc(month),nyear
         else if(nfile.eq.n_month) then
            nday=mod(idate2,10000)/100
            if(nday.ne.1) nrecord = nday
            write(31,400)nrecord,cday(1),chmc(month),nyear
         else
            write(31,400)nrecord,cday(1),chmc(month),nyear
         endif
 400  format('tdef ',I4,' linear ',A2,A3,I4,' ','1dy')
      write(31,500) 4+2
 500  format('vars ',I3)
 600  format(A8,'0 99 ',A36)
 611  format(A8,'0 33,105 ',A36)
 612  format(A8,'0 34,105 ',A36)
      if(filename.eq.'ICBC') then
 660  format(A8,I2,' 0 ',A36)
 661  format(A8,I2,' 33,100 ',A36)
 662  format(A8,I2,' 34,100 ',A36)
        if(iproj.eq.'LAMCON') then
           write(31,661) 'u       ',kz,'westerly wind (m/s)         '
           write(31,662) 'v       ',kz,'southerly wind (m/s)        '
        else
           write(31,660) 'u       ',kz,'westerly wind (m/s)         '
           write(31,660) 'v       ',kz,'southerly wind (m/s)        '
        endif
        write(31,660) 't       ',kz,'air temperature (degree)    '
        write(31,660) 'q       ',kz,'specific moisture (kg/kg)   '
        write(31,600) 'psa     ',   'surface pressure (hPa)      '
        write(31,600) 'ts      ',   'surface air temperature     '
      endif
      write(31,700)
 700  format('endvars')
      close(31)

         open(10,file=trim(Path_Input)//trim(DomainName)//filein &
                ,form='unformatted',recl=iy*jx*ibyte,access='direct')
         open(20,file=trim(Path_Input)//trim(DomainName)//fileout  &
          ,form='unformatted',recl=(iy-2)*(jx-2)*ibyte,access='direct')
         mrec = 0
         do nday=1,nrecord
            do n=1,nvar
               do i=1,iy-2
               do j=1,jx-2
                  c(j,i,n) = 0.0
               enddo
               enddo
            enddo
            do mday = 1,n_slice
               nrec=nrec+1
               read(10,rec=nrec) idatex
               do n=1,nvar
                  nrec=nrec+1
                  read(10,rec=nrec) b
                  do i=1,iy-2
                  do j=1,jx-2
                     c(j,i,n) = c(j,i,n)+b(j+1,i+1)
                  enddo
                  enddo
               enddo
            enddo
            do n=1,nvar
               do i=1,iy-2
               do j=1,jx-2
                  c(j,i,n) = c(j,i,n)/float(n_slice)
               enddo
               enddo
               mrec=mrec+1
               write(20,rec=mrec) ((c(j,i,n),j=1,jx-2),i=1,iy-2)
            enddo
         enddo
         close(10)
         close(20)
      enddo
      deallocate(sigma)
      deallocate(b)
      deallocate(xlat)
      deallocate(xlon)
      deallocate(c)
      return
      end

      subroutine day3(filename,iy,jx,kz,nsg,ibyte, Path_Input &
                     ,DomainName,Path_Output ,idate0,idate1,idate2)
      implicit none
      character*3 filename
      character*80 Path_Input,Path_Output
      character*10 DomainName
      integer iy,jx,kz,nsg,ibyte,idate0,idate1,idate2
      integer iiy,jjx,kkz
      integer mdate0,ibltyp,icup,ipptls,iboudy
      real*4  truelatL,truelatH
      real*4  dxsp,ptsp,clat,clon,plat,plon
      real*4  dto,dtb,dtr,dtc,dt
      integer iotyp
      character*6 iproj
      real*4, allocatable ::  sigma(:)
      character*4 :: chy(1941:2100)
      data chy/'1941','1942','1943','1944','1945','1946','1947','1948',&
        '1949','1950','1951','1952','1953','1954','1955','1956','1957',&
        '1958','1959','1960','1961','1962','1963','1964','1965','1966',&
        '1967','1968','1969','1970','1971','1972','1973','1974','1975',&
        '1976','1977','1978','1979','1980','1981','1982','1983','1984',&
        '1985','1986','1987','1988','1989','1990','1991','1992','1993',&
        '1994','1995','1996','1997','1998','1999','2000','2001','2002',&
        '2003','2004','2005','2006','2007','2008','2009','2010','2011',&
        '2012','2013','2014','2015','2016','2017','2018','2019','2020',&
        '2021','2022','2023','2024','2025','2026','2027','2028','2029',&
        '2030','2031','2032','2033','2034','2035','2036','2037','2038',&
        '2039','2040','2041','2042','2043','2044','2045','2046','2047',&
        '2048','2049','2050','2051','2052','2053','2054','2055','2056',&
        '2057','2058','2059','2060','2061','2062','2063','2064','2065',&
        '2066','2067','2068','2069','2070','2071','2072','2073','2074',&
        '2075','2076','2077','2078','2079','2080','2081','2082','2083',&
        '2084','2085','2086','2087','2088','2089','2090','2091','2092',&
        '2093','2094','2095','2096','2097','2098','2099','2100'/
      character*2 cday(31)
      data cday/'01','02','03','04','05','06','07','08','09','10', &
                '11','12','13','14','15','16','17','18','19','20', &
                '21','22','23','24','25','26','27','28','29','30','31'/
      character*2 chm(12)
      data chm/'01','02','03','04','05','06','07','08','09','10', &
               '11','12'/
      character*3 chmc(12)
      data chmc/'jan','feb','mar','apr','may','jun'  &
               ,'jul','aug','sep','oct','nov','dec'/
      character*14 filein
      character*12 fileout
      integer ntype,nfile,nyear,month,n_slice,mrec,nrec
      logical there
      real*4  alatmin,alatmax,alonmin,alonmax,rlatinc,rloninc
      real*4  centerj,centeri
      integer ny,nx
      integer i,j,k,l,m,n,itr
      integer number,nvar,mvar,n_month,nrecord,nday,mday
      integer idatex
      real*4, allocatable ::  c(:,:,:),b(:,:),xlat(:,:),xlon(:,:)

      allocate(sigma(kz+1))
      allocate(b((jx-2)*nsg,(iy-2)*nsg))
      allocate(xlat(jx-2,iy-2))
      allocate(xlon(jx-2,iy-2))

      open(10,file=trim(Path_Output)//'OUT_HEAD',form='unformatted' &
             ,recl=(jx-2)*(iy-2)*ibyte,access='direct')
      read(10,rec=1) mdate0,ibltyp,icup,ipptls,iboudy  &
                    ,iiy,jjx,kkz,(sigma(k),k=1,kz+1)   &
                    ,dxsp,ptsp,clat,clon,plat,plon     &
                    ,iproj,dto,dtb,dtr,dtc,iotyp,truelatL,truelatH
      if(iiy.ne.iy.or.jjx.ne.jx.or.kkz.ne.kz) then
         write(*,*) 'iy,jx,kz in parameter = ',iy,jx,kz
         write(*,*) 'iy,jx,kz in OUT_HEAD ',iiy,jjx,kkz
         write(*,*) 'They are not consistent'
         stop
      endif
      read(10,rec=6) xlat
      read(10,rec=7) xlon
      close(10)

      if(filename.eq.'SUB') then
         dt = dtb
         nvar = 16
      else
         write(*,*) 'filename is not correct'
         stop
      endif
      allocate(c((jx-2)*nsg,(iy-2)*nsg,nvar))

      if(idate1.gt.1940010100) then        ! original file
         n_month = (idate2/1000000-idate1/1000000)*12 &
                 + (mod(idate2/10000,100)-mod(idate1/10000,100))
         if(mod(idate2,10000).gt.0100) n_month = n_month+1
         ntype = 0
      else
         write(*,*) 'date should be in 10 digit'
         stop
      endif

      do nfile=1,n_month
         nyear = idate1/1000000
         month = idate1/10000-nyear*100 +nfile-1
         nyear = nyear + (month-1)/12
         month = mod(month,12)
         if(month.eq.0) month = 12

         if(nfile.eq.1.or.month.eq.1) then
         write(*,*) 'Calculate the daily mean of ',filename,nyear,month
         else
         write(*,*) '                            ',filename,nyear,month
         endif

         nrec = 0
         if(month.eq.1.or.month.eq.3.or.month.eq.5.or.  &
            month.eq.7.or.month.eq.8.or.month.eq.10.or. &
            month.eq.12) then
            nrecord = 31
         else if(month.eq.4.or.month.eq.6.or.month.eq.9.or. &
                 month.eq.11) then
            nrecord = 30
         else
            nrecord = 28
            if(mod(nyear,4).eq.0) nrecord = 29
            if(mod(nyear,100).eq.0) nrecord = 28
            if(mod(nyear,400).eq.0) nrecord = 29
         endif
         if(nfile.eq.n_month.and.mod(idate2/100-1,100).ne.0) &
            nrecord = min(nrecord,mod(idate2/100-1,100))
         if(nfile.eq.1.and.idate0.eq.idate1) nrec = nvar
         n_slice = 24/nint(dt)

         if(filename.eq.'SUB') then
            filein = 'SUB.'//chy(nyear)//chm(month)//'0100'
            fileout= 'SUB.'//chy(nyear)//chm(month)//'01'
         endif

         inquire(file=trim(Path_Output)//fileout//'.ctl',exist=there)
         if(there) then
           open(31,file=trim(Path_Output)//fileout//'.ctl' &
                                                      ,status='replace')
         else
           open(31,file=trim(Path_Output)//fileout//'.ctl',status='new')
         endif
         write(31,10) fileout
  10     format('dset ^',A12)
         write(31,20)
  20     format('title RegCM daily output variables')
         write(31,30)
  30     format('options big_endian')
         write(31,50)
  50     format('undef -1.e34')
         if(iproj.eq.'LAMCON'.or.iproj.eq.'ROTMER') then
            alatmin= 999999.
            alatmax=-999999.
            do j=1,jx-2
               if(xlat(j,1).lt.alatmin) alatmin=xlat(j,1)
               if(xlat(j,iy-2).gt.alatmax) alatmax=xlat(j,iy-2)
            enddo
            alonmin= 999999.
            alonmax=-999999.
            do i=1,iy-2
            do j=1,jx-2
               if(clon.ge.0.0) then
                  if(xlon(j,i).ge.0.0) then
                     alonmin = amin1(alonmin,xlon(j,i))
                     alonmax = amax1(alonmax,xlon(j,i))
                  else
                     if(abs(clon-xlon(j,i)).lt. &
                        abs(clon-(xlon(j,i)+360.))) then
                        alonmin = amin1(alonmin,xlon(j,i))
                        alonmax = amax1(alonmax,xlon(j,i))
                     else
                        alonmin = amin1(alonmin,xlon(j,i)+360.)
                        alonmax = amax1(alonmax,xlon(j,i)+360.)
                     endif
                  endif
               else
                  if(xlon(j,i).lt.0.0) then
                     alonmin = amin1(alonmin,xlon(j,i))
                     alonmax = amax1(alonmax,xlon(j,i))
                  else
                     if(abs(clon-xlon(j,i).lt. &
                        abs(clon-xlon(j,i)-360.))) then
                        alonmin = amin1(alonmin,xlon(j,i))
                        alonmax = amax1(alonmax,xlon(j,i))
                     else
                        alonmin = amin1(alonmin,xlon(j,i)-360.)
                        alonmax = amax1(alonmax,xlon(j,i)-360.)
                     endif
                  endif
               endif
            enddo
            enddo
            rlatinc=dxsp/111./2.
            rloninc=dxsp/111./2.
            ny=2+nint(abs(alatmax-alatmin)/rlatinc)
            nx=1+nint(abs((alonmax-alonmin)/rloninc))
            centerj=(jx-2)/2.
            centeri=(iy-2)/2.
         endif

         if(iproj.eq.'LAMCON') then        ! Lambert projection
            write(31,100) jx-2,iy-2,clat,clon,centerj,centeri, &
                          truelatL,truelatH,clon,dxsp*1000.,dxsp*1000.
 100  format('pdef ',i4,1x,i4,1x,'lccr',7(1x,f7.2),1x,2(f7.0,1x))
            write(31,110) nx+2,alonmin-rloninc,rloninc
 110  format('xdef ',i4,' linear ',f7.2,1x,f7.4)
            write(31,120) ny+2,alatmin-rlatinc,rlatinc
 120  format('ydef ',i4,' linear ',f7.2,1x,f7.4)
         elseif(iproj.eq.'POLSTR') then    !
         elseif(iproj.eq.'NORMER') then
            write(31,200)  jx-2,xlon(1,1),xlon(2,1)-xlon(1,1)
 200  format('xdef ',I3,' linear ',f9.4,' ',f9.4)
            write(31,210) iy-2
 210  format('ydef ',I3,' levels')
            write(31,220) (xlat(1,i),i=1,iy-2)
 220  format(10f7.2)
         elseif(iproj.eq.'ROTMER') then
            write(*,*) 'Note that rotated Mercartor (ROTMER)' &
                   ,' projections are not supported by GrADS.'
            write(*,*) '  Although not exact, the eta.u projection' &
                   ,' in GrADS is somewhat similar.'
            write(*,*) ' FERRET, however, does support this projection.'
            write(31,230) jx-2,iy-2,plon,plat,dxsp/111. &
                                          ,dxsp/111.*.95238
 230  format('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
            write(31,110) nx+2,alonmin-rloninc,rloninc
            write(31,120) ny+2,alatmin-rlatinc,rlatinc
         else
            write(*,*) 'Are you sure your map projection is right ?'
            stop
         endif
         if(filename.eq.'SRF') then
      write(31,300) (1013.25-ptsp*10.)*(sigma(1)+sigma(2))*.5+ptsp*10.
         else
      write(31,318) kz,((1013.25-ptsp*10.)*(sigma(k)+sigma(k+1))*.5 &
                       +ptsp*10.,k=1,kz)
         endif
 300  format('zdef 1',' levels ',f7.2)
 318  format('zdef ',I2,' levels ',30f7.2)
         if(nfile.eq.1) then
            nday=mod(idate0,10000)/100
            nrecord = nrecord+1-nday 
            write(31,400)nrecord,cday(nday),chmc(month),nyear
         else if(nfile.eq.n_month) then
            nday=mod(idate2,10000)/100
            if(nday.ne.1) nrecord = nday
            write(31,400)nrecord,cday(1),chmc(month),nyear
         else
            write(31,400)nrecord,cday(1),chmc(month),nyear
         endif
 400  format('tdef ',I4,' linear ',A2,A3,I4,' ','1dy')
         write(31,500) nvar
 500  format('vars ',I3)
 600  format(A8,'0 99 ',A36)
 611  format(A8,'0 33,105 ',A36)
 612  format(A8,'0 34,105 ',A36)
         if(filename.eq.'SRF') then
          if(iproj.eq.'LAMCON') then        ! Lambert projection
        write(31,611) 'u10m    ','westerly  wind at 10m (m/s)          '
        write(31,612) 'v10m    ','southerly wind at 10m (m/s)          '
          else
        write(31,600) 'u10m    ','westerly  wind at 10m (m/s)          '
        write(31,600) 'v10m    ','southerly wind at 10m (m/s)          '
          endif
        write(31,600) 'uvdrag  ','surface drag stress                  '
        write(31,600) 'tg      ','ground temperature (degree)          '
        write(31,600) 'tlef    ','temperature of foliage               '
        write(31,600) 't2m     ','air temperature at 2m (K)            '
        write(31,600) 'q2m     ','water vapor mixing ratio at 2m(kg/kg)'
        write(31,600) 'ssw     ','upper layer soil water               '
        write(31,600) 'rsw     ','root zone soil water                 '
        write(31,600) 'tpr     ','total precipitation (mm/day)         '
        write(31,600) 'evp     ','evapotranspiration (mm/day)          '
        write(31,600) 'runoff  ','surface runoff (mm/day)              '
        write(31,600) 'scv     ','total snow amount                    '
        write(31,600) 'sena    ','sensible heat flux (W/m2)            '
        write(31,600) 'flw     ','net infrared energy flux (W/m2)      '
        write(31,600) 'fsw     ','net absorbed solar energy flux (W/m2)'
        write(31,600) 'flwd    ','downward infrared energy flux (W/m2) '
        write(31,600) 'sina    ','incident solar energy flux (W/m2)    '
        write(31,600) 'prcv    ','convective precipitation (mm/day)    '
        write(31,600) 'psb     ','surface pressure (hPa)               '
        write(31,600) 'zpbl    ','PBL layer height                     '
        write(31,600) 'tgmax   ','maximum ground temperature (K)       '
        write(31,600) 'tgmin   ','minimum ground temperature (K)       '
        write(31,600) 't2max   ','maximum 2m air temperature (K)       '
        write(31,600) 't2min   ','minimum 2m air temperature (K)       '
        write(31,600) 'w10max  ','maximum 10m wind speed (m/s)         '
        write(31,600) 'ps_min  ','minimum surface pressure (hPa)       '
         endif
         write(31,700)
 700  format('endvars')
         close(31)

         if(iotyp.eq.1) then
            open(10,file=trim(Path_Output)//filein,form='unformatted' &
                   ,recl=(iy-2)*nsg*(jx-2)*nsg*ibyte,access='direct')
         else if(iotyp.eq.2) then
            open(10,file=trim(Path_Output)//filein,form='unformatted')
         endif
         open(20,file=fileout,form='unformatted' &
                ,recl=(iy-2)*nsg*(jx-2)*nsg*ibyte,access='direct')
         mrec = 0
         do nday=1,nrecord
            do n=1,mvar
               do i=1,(iy-2)*nsg
               do j=1,(jx-2)*nsg
                  c(j,i,n) = 0.0
               enddo
               enddo
            enddo
            if(iotyp.eq.2) read(10) idatex
            do mday = 1,n_slice
               do n=1,mvar
                  if(iotyp.eq.1) then
                     nrec=nrec+1
                     read(10,rec=nrec) b
                  else if(iotyp.eq.2) then
                     read(10) b
                  endif
                  do i=1,(iy-2)*nsg
                  do j=1,(jx-2)*nsg
                     c(j,i,n) = c(j,i,n)+b(j,i)
                  enddo
                  enddo
               enddo
            enddo
            do n=1,mvar
               do i=1,(iy-2)*nsg
               do j=1,(jx-2)*nsg
                  c(j,i,n) = c(j,i,n)/float(n_slice)
               enddo
               enddo
               mrec=mrec+1
               write(20,rec=mrec) ((c(j,i,n),j=1,(jx-2)*nsg)  &
                                            ,i=1,(iy-2)*nsg)
            enddo
         enddo
         close(10)
         close(20)
      enddo
      deallocate(sigma)
      deallocate(b)
      deallocate(xlat)
      deallocate(xlon)
      deallocate(c)
      return
      end

      subroutine mon(filename,iy,jx,kz,ntr,ibyte, Path_Output &
                    ,idate0,idate1,idate2)
      implicit none
      character*3 filename
      character*80 Path_Output
      integer iy,jx,kz,ntr,ibyte,idate0,idate1,idate2
      integer iiy,jjx,kkz
      integer mdate0,ibltyp,icup,ipptls,iboudy
      real*4  truelatL,truelatH
      real*4  dxsp,ptsp,clat,clon,plat,plon
      real*4  dto,dtb,dtr,dtc
      integer iotyp
      character*6 iproj
      real*4, allocatable ::  sigma(:)
      character*4 :: chy(1941:2100)
      data chy/'1941','1942','1943','1944','1945','1946','1947','1948',&
        '1949','1950','1951','1952','1953','1954','1955','1956','1957',&
        '1958','1959','1960','1961','1962','1963','1964','1965','1966',&
        '1967','1968','1969','1970','1971','1972','1973','1974','1975',&
        '1976','1977','1978','1979','1980','1981','1982','1983','1984',&
        '1985','1986','1987','1988','1989','1990','1991','1992','1993',&
        '1994','1995','1996','1997','1998','1999','2000','2001','2002',&
        '2003','2004','2005','2006','2007','2008','2009','2010','2011',&
        '2012','2013','2014','2015','2016','2017','2018','2019','2020',&
        '2021','2022','2023','2024','2025','2026','2027','2028','2029',&
        '2030','2031','2032','2033','2034','2035','2036','2037','2038',&
        '2039','2040','2041','2042','2043','2044','2045','2046','2047',&
        '2048','2049','2050','2051','2052','2053','2054','2055','2056',&
        '2057','2058','2059','2060','2061','2062','2063','2064','2065',&
        '2066','2067','2068','2069','2070','2071','2072','2073','2074',&
        '2075','2076','2077','2078','2079','2080','2081','2082','2083',&
        '2084','2085','2086','2087','2088','2089','2090','2091','2092',&
        '2093','2094','2095','2096','2097','2098','2099','2100'/
      character*2 cday(31)
      data cday/'01','02','03','04','05','06','07','08','09','10', &
                '11','12','13','14','15','16','17','18','19','20', &
                '21','22','23','24','25','26','27','28','29','30','31'/
      character*2 chm(12)
      data chm/'01','02','03','04','05','06','07','08','09','10', &
               '11','12'/
      character*3 chmc(12)
      data chmc/'jan','feb','mar','apr','may','jun'  &
               ,'jul','aug','sep','oct','nov','dec'/
      character*12 filein
      character*7 fileout
      integer ntype,nfile,nyear,month,mrec,nrec
      integer i,j,k,l,m,n,itr
      integer number,nvar,n_month,nrecord,nday,mday
      logical there
      real*4  alatmin,alatmax,alonmin,alonmax,rlatinc,rloninc
      real*4  centerj,centeri
      integer ny,nx

      real*4, allocatable ::  c(:,:,:),b(:,:),xlat(:,:),xlon(:,:)

      allocate(sigma(kz+1))
      allocate(b(jx-2,iy-2))
      allocate(xlat(jx-2,iy-2))
      allocate(xlon(jx-2,iy-2))

      open(10,file=trim(Path_Output)//'OUT_HEAD',form='unformatted' &
             ,recl=(jx-2)*(iy-2)*ibyte,access='direct')
      read(10,rec=1) mdate0,ibltyp,icup,ipptls,iboudy  &
                    ,iiy,jjx,kkz,(sigma(k),k=1,kz+1)   &
                    ,dxsp,ptsp,clat,clon,plat,plon     &
                    ,iproj,dto,dtb,dtr,dtc,iotyp,truelatL,truelatH
      if(iiy.ne.iy.or.jjx.ne.jx.or.kkz.ne.kz) then
         write(*,*) 'iy,jx,kz in parameter = ',iy,jx,kz
         write(*,*) 'iy,jx,kz in OUT_HEAD ',iiy,jjx,kkz
         write(*,*) 'They are not consistent'
         stop
      endif
      read(10,rec=6) xlat
      read(10,rec=7) xlon
      close(10)

      if(filename.eq.'ATM') then
         nvar = kz*6+5
      else if(filename.eq.'SRF') then
         nvar = 27
      else if(filename.eq.'RAD') then
         nvar = kz*4+9
      else if(filename.eq.'CHE') then
         nvar = ntr*kz+kz*3+ntr*7+3
!        nvar = ntr*kz+kz*3+ntr*7+2      ! for RegCM3, one record less
      else
         write(*,*) 'filename is not correct'
         stop
      endif
      allocate(c(jx-2,iy-2,nvar))

      if(idate1.gt.1940010100) then        ! original file
         n_month = (idate2/1000000-idate1/1000000)*12 &
                 + (mod(idate2/10000,100)-mod(idate1/10000,100))
         if(mod(idate2,10000).gt.0100) n_month = n_month+1
         ntype = 0
      else if(idate1.gt.19400101) then     ! daily mean file
         n_month = (idate2/10000-idate1/10000)*12 &
                 + (mod(idate2/100,100)-mod(idate1/100,100))
         if(mod(idate2,100).gt.01) n_month = n_month+1
         ntype = 1
      else
         write(*,*) 'date should be in 10 digits, or 8 digits'
         stop
      endif

      if(filename.eq.'ATM') then
         fileout= 'ATM.mon'
      else if(filename.eq.'SRF') then
         fileout= 'SRF.mon'
      else if(filename.eq.'RAD') then
         fileout= 'RAD.mon'
      else if(filename.eq.'CHE') then
         fileout= 'CHE.mon'
      endif

      open(20,file=trim(Path_Output)//fileout,form='unformatted' &
             ,recl=(iy-2)*(jx-2)*ibyte,access='direct')
      mrec = 0

      inquire(file=trim(Path_Output)//fileout//'.ctl',exist=there)
      if(there) then
       open(31,file=trim(Path_Output)//fileout//'.ctl',status='replace')
      else
         open(31,file=trim(Path_Output)//fileout//'.ctl',status='new')
      endif
      write(31,10) fileout
  10  format('dset ^',A7)
      write(31,20)
  20  format('title RegCM monthly output variables')
      write(31,30)
  30  format('options big_endian')
      write(31,50)
  50  format('undef -1.e34')
      if(iproj.eq.'LAMCON'.or.iproj.eq.'ROTMER') then
         alatmin= 999999.
         alatmax=-999999.
         do j=1,jx-2
            if(xlat(j,1).lt.alatmin) alatmin=xlat(j,1)
            if(xlat(j,iy-2).gt.alatmax) alatmax=xlat(j,iy-2)
         enddo
         alonmin= 999999.
         alonmax=-999999.
         do i=1,iy-2
         do j=1,jx-2
            if(clon.ge.0.0) then
               if(xlon(j,i).ge.0.0) then
                  alonmin = amin1(alonmin,xlon(j,i))
                  alonmax = amax1(alonmax,xlon(j,i))
               else
                  if(abs(clon-xlon(j,i)).lt. &
                     abs(clon-(xlon(j,i)+360.))) then
                     alonmin = amin1(alonmin,xlon(j,i))
                     alonmax = amax1(alonmax,xlon(j,i))
                  else
                     alonmin = amin1(alonmin,xlon(j,i)+360.)
                     alonmax = amax1(alonmax,xlon(j,i)+360.)
                  endif
               endif
            else
               if(xlon(j,i).lt.0.0) then
                  alonmin = amin1(alonmin,xlon(j,i))
                  alonmax = amax1(alonmax,xlon(j,i))
               else
                  if(abs(clon-xlon(j,i).lt. &
                     abs(clon-xlon(j,i)-360.))) then
                     alonmin = amin1(alonmin,xlon(j,i))
                     alonmax = amax1(alonmax,xlon(j,i))
                  else
                     alonmin = amin1(alonmin,xlon(j,i)-360.)
                     alonmax = amax1(alonmax,xlon(j,i)-360.)
                  endif
               endif
            endif
         enddo
         enddo
         rlatinc=dxsp/111./2.
         rloninc=dxsp/111./2.
         ny=2+nint(abs(alatmax-alatmin)/rlatinc)
         nx=1+nint(abs((alonmax-alonmin)/rloninc))
         centerj=(jx-2)/2.
         centeri=(iy-2)/2.
      endif

      if(iproj.eq.'LAMCON') then        ! Lambert projection
         write(31,100) jx-2,iy-2,clat,clon,centerj,centeri, &
                       truelatL,truelatH,clon,dxsp*1000.,dxsp*1000.
 100  format('pdef ',i4,1x,i4,1x,'lccr',7(1x,f7.2),1x,2(f7.0,1x))
         write(31,110) nx+2,alonmin-rloninc,rloninc
 110  format('xdef ',i4,' linear ',f7.2,1x,f7.4)
         write(31,120) ny+2,alatmin-rlatinc,rlatinc
 120  format('ydef ',i4,' linear ',f7.2,1x,f7.4)
      elseif(iproj.eq.'POLSTR') then    !
      elseif(iproj.eq.'NORMER') then
         write(31,200)  jx-2,xlon(1,1),xlon(2,1)-xlon(1,1)
 200  format('xdef ',I3,' linear ',f9.4,' ',f9.4)
         write(31,210) iy-2
 210  format('ydef ',I3,' levels')
         write(31,220) (xlat(1,i),i=1,iy-2)
 220  format(10f7.2)
      elseif(iproj.eq.'ROTMER') then
         write(*,*) 'Note that rotated Mercartor (ROTMER)' &
                   ,' projections are not supported by GrADS.'
         write(*,*) '  Although not exact, the eta.u projection' &
                   ,' in GrADS is somewhat similar.'
         write(*,*) ' FERRET, however, does support this projection.'
         write(31,230) jx-2,iy-2,plon,plat,dxsp/111. &
                                          ,dxsp/111.*.95238
 230  format('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
         write(31,110) nx+2,alonmin-rloninc,rloninc
         write(31,120) ny+2,alatmin-rlatinc,rlatinc
      else
         write(*,*) 'Are you sure your map projection is right ?'
         stop
      endif
      if(filename.eq.'SRF'.or.filename.eq.'SUB') then
      write(31,300) (1013.25-ptsp*10.)*(sigma(1)+sigma(2))*.5+ptsp*10.
      else
      write(31,318) kz,((1013.25-ptsp*10.)*(sigma(k)+sigma(k+1))*.5 &
                       +ptsp*10.,k=1,kz)
      endif
 300  format('zdef 1',' levels ',f7.2)
 318  format('zdef ',I2,' levels ',30f7.2)
      nyear=idate0/1000000
      month=(idate0-nyear*1000000)/10000
      nday =(idate0-nyear*1000000-month*10000)/100
      write(31,400)n_month,cday(16),chmc(month),nyear
 400  format('tdef ',I4,' linear ',A2,A3,I4,' ','1mo')
      if(filename.eq.'ATM') then
         write(31,500) 6+5
      else if(filename.eq.'RAD') then
         write(31,500) 4+9
      else if(filename.eq.'SRF') then
         write(31,500) nvar
      else if(filename.eq.'CHE') then
         write(31,500) ntr*8+5
      endif
 500  format('vars ',I3)
 600  format(A8,'0 99 ',A36)
 611  format(A8,'0 33,105 ',A36)
 612  format(A8,'0 34,105 ',A36)
      if(filename.eq.'SRF') then
       if(iproj.eq.'LAMCON') then        ! Lambert projection
        write(31,611) 'u10m    ','westerly  wind at 10m (m/s)          '
        write(31,612) 'v10m    ','southerly wind at 10m (m/s)          '
       else
        write(31,600) 'u10m    ','westerly  wind at 10m (m/s)          '
        write(31,600) 'v10m    ','southerly wind at 10m (m/s)          '
       endif
        write(31,600) 'uvdrag  ','surface drag stress                  '
        write(31,600) 'tg      ','ground temperature (degree)          '
        write(31,600) 'tlef    ','temperature of foliage               '
        write(31,600) 't2m     ','air temperature at 2m (K)            '
        write(31,600) 'q2m     ','water vapor mixing ratio at 2m(kg/kg)'
        write(31,600) 'ssw     ','upper layer soil water               '
        write(31,600) 'rsw     ','root zone soil water                 '
        write(31,600) 'tpr     ','total precipitation (mm/day)         '
        write(31,600) 'evp     ','evapotranspiration (mm/day)          '
        write(31,600) 'runoff  ','surface runoff (mm/day)              '
        write(31,600) 'scv     ','total snow amount                    '
        write(31,600) 'sena    ','sensible heat flux (W/m2)            '
        write(31,600) 'flw     ','net infrared energy flux (W/m2)      '
        write(31,600) 'fsw     ','net absorbed solar energy flux (W/m2)'
        write(31,600) 'flwd    ','downward infrared energy flux (W/m2) '
        write(31,600) 'sina    ','incident solar energy flux (W/m2)    '
        write(31,600) 'prcv    ','convective precipitation (mm/day)    '
        write(31,600) 'psb     ','surface pressure (hPa)               '
        write(31,600) 'zpbl    ','PBL layer height                     '
        write(31,600) 'tgmax   ','maximum ground temperature (K)       '
        write(31,600) 'tgmin   ','minimum ground temperature (K)       '
        write(31,600) 't2max   ','maximum 2m air temperature (K)       '
        write(31,600) 't2min   ','minimum 2m air temperature (K)       '
        write(31,600) 'w10max  ','maximum 10m wind speed (m/s)         '
        write(31,600) 'ps_min  ','minimum surface pressure (hPa)       '
      else if(filename.eq.'ATM') then
 660  format(A8,I2,' 0 ',A36)
 661  format(A8,I2,' 33,100 ',A36)
 662  format(A8,I2,' 34,100 ',A36)
        if(iproj.eq.'LAMCON') then
           write(31,661) 'u       ',kz,'westerly wind (m/s)         '
           write(31,662) 'v       ',kz,'southerly wind (m/s)        '
        else
           write(31,660) 'u       ',kz,'westerly wind (m/s)         '
           write(31,660) 'v       ',kz,'southerly wind (m/s)        '
        endif
        write(31,660) 'w       ',kz,'omega (hPa/s)   p-velocity  '
        write(31,660) 't       ',kz,'air temperature (degree)    '
        write(31,660) 'qv      ',kz,'water vapor mixing ratio    '
        write(31,660) 'qc      ',kz,'cloud water mixing ratio    '
        write(31,600) 'psa     ',   'surface pressure (hPa)      '
        write(31,600) 'tpr     ',   'total precipitation(mm/day) '
        write(31,600) 'tgb     ',   'lower groud temp. in BATS   '
        write(31,600) 'swt     ',   'total soil water in mm H2O  '
        write(31,600) 'rno     ',   'accumulated infiltration    '
      else if(filename.eq.'RAD') then
        write(31,660)'cld   ',kz,'cloud fractional cover               '
        write(31,660)'clwp  ',kz,'cloud liquid water path              '
        write(31,660)'qrs   ',kz,'solar heating rate                   '
        write(31,660)'qrl   ',kz,'longwave cooling rate                '
        write(31,600)'frsa  ',   'surface absorbed solar flux          '
        write(31,600)'frla  ',   'longwave cooling of surface          '
        write(31,600)'clrst ',   'clearsky total column abs solar flux '
        write(31,600)'clrss ',   'clearsky surface absorbed solar flux '
        write(31,600)'clrlt ',   'clearsky net upward LW flux at TOA   '
        write(31,600)'clrls ',   'clearsky LW cooling at surface (W/m2)'
        write(31,600)'solin ',   'instantaneous incident solar (W/m2)  '
        write(31,600)'sabtp ',   'total column absorbed solar flux W/m2'
        write(31,600)'firtp ',   'net upward LW flux at TOA (W/m2)     '
      else if(filename.eq.'CHE') then
       do itr=1,ntr
         if(itr.lt.10) then
           write(31,650) 'trac',itr,kz, 'tracer mix. rat  (Kg/Kg)  '
         else
           write(31,651) 'trac',itr,kz, 'tracer mix. rat  (Kg/Kg)  '
         endif
       enddo
 650  format(A4,I1,' ',I2,' 0 ',A36)
 651  format(A4,I2,' ',I2,' 0 ',A36)
 655  format(A8,I1,' 0 99 ',A36)
 656  format(A8,I2,' 0 99 ',A36)
       write(31,650) 'aext',8,kz, 'aer mix. ext. coef      '
       write(31,650) 'assa',8,kz, 'aer mix. sin. scat. alb '
       write(31,650) 'agfu',8,kz, 'aer mix. ass. par       '
       do itr=1,ntr
        if(itr.lt.10) then
          write(31,655) 'colb__tr',itr,'columnburden inst(mg/m2)'
          write(31,655) 'wdlsc_tr',itr,'wet dep lgscale(mg/m2/d)'
          write(31,655) 'wdcvc_tr',itr,'wet dep convect(mg/m2/d)'
          write(31,655) 'sdrdp_tr',itr,'surf dry depos.(mg/m2/d)'
          write(31,655) 'xgasc_tr',itr,'chem gas conv. (mg/m2/d)'
          write(31,655) 'xaquc_tr',itr,'chem aqu conv. (mg/m2/d)'
          write(31,655) 'emiss_tr',itr,'surf emission  (mg/m2/d)'
        else
          write(31,656) 'colb__tr',itr,'columnburden inst(mg/m2)'
          write(31,656) 'wdlsc_tr',itr,'wet dep lgscale(mg/m2/d)'
          write(31,656) 'wdcvc_tr',itr,'wet dep convect(mg/m2/d)'
          write(31,656) 'sdrdp_tr',itr,'surf dry depos.(mg/m2/d)'
          write(31,656) 'xgasc_tr',itr,'chem gas conv. (mg/m2/d)'
          write(31,656) 'xaquc_tr',itr,'chem aqu conv. (mg/m2/d)'
          write(31,656) 'emiss_tr',itr,'surf emission  (mg/m2/d)'
        endif
       enddo
       write(31,600) 'acstoarf',' TOArad forcing av.(W/m2)            '
       write(31,600) 'acstsrrf',' SRFrad forcing av.(W/m2)            '
       write(31,600) 'psb     ','surface pressure (hPa)               '
      endif
      write(31,700)
 700  format('endvars')
      close(31)

      do nfile=1,n_month
         if(ntype.eq.0) then
            nyear = idate1/1000000
            month = idate1/10000-nyear*100 +nfile-1
         else if(ntype.eq.1) then
            nyear = idate1/10000
            month = idate1/100-nyear*100 +nfile-1
         endif

         nyear = nyear + (month-1)/12
         month = mod(month,12)
         if(month.eq.0) month = 12

         if(nfile.eq.1.or.month.eq.1) then
         write(*,*)'Calculate the monthly mean of ',filename,nyear,month
         else
         write(*,*)'                              ',filename,nyear,month
         endif

         if(month.eq.1.or.month.eq.3.or.month.eq.5.or.  &
            month.eq.7.or.month.eq.8.or.month.eq.10.or. &
            month.eq.12) then
            nrecord = 31
         else if(month.eq.4.or.month.eq.6.or.month.eq.9.or. &
                 month.eq.11) then
            nrecord = 30
         else
            nrecord = 28
            if(mod(nyear,4).eq.0) nrecord = 29
            if(mod(nyear,100).eq.0) nrecord = 28
            if(mod(nyear,400).eq.0) nrecord = 29
         endif
         if(ntype.eq.0) then
            if(nfile.eq.n_month.and.mod(idate2/100-1,100).ne.0) &
               nrecord = min(nrecord,mod(idate2/100-1,100))
         else if(ntype.eq.1) then
            if(nfile.eq.n_month.and.mod(idate2-1,100).ne.0) &
               nrecord = min(nrecord,mod(idate2-1,100)+1)
         endif

         if(filename.eq.'ATM') then
            filein = 'ATM.'//chy(nyear)//chm(month)//'01'
         else if(filename.eq.'SRF') then
            filein = 'SRF.'//chy(nyear)//chm(month)//'01'
         else if(filename.eq.'RAD') then
            filein = 'RAD.'//chy(nyear)//chm(month)//'01'
         else if(filename.eq.'CHE') then
            filein = 'CHE.'//chy(nyear)//chm(month)//'01'
         endif
         open(10,file=trim(Path_Output)//filein,form='unformatted' &
                ,recl=(iy-2)*(jx-2)*ibyte,access='direct')
         nrec = 0

         do n=1,nvar
            do i=1,iy-2
            do j=1,jx-2
               c(j,i,n) = 0.0
            enddo
            enddo
         enddo

         do nday=1,nrecord
            do n=1,nvar
               nrec=nrec+1
               read(10,rec=nrec) b
               do i=1,iy-2
               do j=1,jx-2
                  c(j,i,n) = c(j,i,n)+b(j,i)
               enddo
               enddo
            enddo
         enddo
         do n=1,nvar
            do i=1,iy-2
            do j=1,jx-2
               c(j,i,n) = c(j,i,n)/float(nrecord)
            enddo
            enddo
            mrec=mrec+1
            write(20,rec=mrec) ((c(j,i,n),j=1,jx-2),i=1,iy-2)
         enddo
         close(10)
      enddo
      close(20)
      deallocate(sigma)
      deallocate(b)
      deallocate(xlat)
      deallocate(xlon)
      deallocate(c)
      return
      end

      subroutine mon2(filename,iy,jx,kz,ibyte, Path_Input,DomainName &
                     ,idate0,idate1,idate2)
      implicit none
      character*4 filename
      character*80 Path_Input
      character*10 DomainName
      integer iy,jx,kz,ibyte,idate0,idate1,idate2
      integer iiy,jjx,kkz
      real*4  truelatL,truelatH
      real*4  dsinm,ptop,clat,clon,plat,plon,GRDFAC
      integer igrads,ibigend
      character*6 iproj
      real*4, allocatable ::  sigma(:)
      character*4 :: chy(1941:2100)
      data chy/'1941','1942','1943','1944','1945','1946','1947','1948',&
        '1949','1950','1951','1952','1953','1954','1955','1956','1957',&
        '1958','1959','1960','1961','1962','1963','1964','1965','1966',&
        '1967','1968','1969','1970','1971','1972','1973','1974','1975',&
        '1976','1977','1978','1979','1980','1981','1982','1983','1984',&
        '1985','1986','1987','1988','1989','1990','1991','1992','1993',&
        '1994','1995','1996','1997','1998','1999','2000','2001','2002',&
        '2003','2004','2005','2006','2007','2008','2009','2010','2011',&
        '2012','2013','2014','2015','2016','2017','2018','2019','2020',&
        '2021','2022','2023','2024','2025','2026','2027','2028','2029',&
        '2030','2031','2032','2033','2034','2035','2036','2037','2038',&
        '2039','2040','2041','2042','2043','2044','2045','2046','2047',&
        '2048','2049','2050','2051','2052','2053','2054','2055','2056',&
        '2057','2058','2059','2060','2061','2062','2063','2064','2065',&
        '2066','2067','2068','2069','2070','2071','2072','2073','2074',&
        '2075','2076','2077','2078','2079','2080','2081','2082','2083',&
        '2084','2085','2086','2087','2088','2089','2090','2091','2092',&
        '2093','2094','2095','2096','2097','2098','2099','2100'/
      character*2 cday(31)
      data cday/'01','02','03','04','05','06','07','08','09','10', &
                '11','12','13','14','15','16','17','18','19','20', &
                '21','22','23','24','25','26','27','28','29','30','31'/
      character*2 chm(12)
      data chm/'01','02','03','04','05','06','07','08','09','10', &
               '11','12'/
      character*3 chmc(12)
      data chmc/'jan','feb','mar','apr','may','jun'  &
               ,'jul','aug','sep','oct','nov','dec'/
      character*14 filein
      character*9 fileout
      integer ntype,nfile,nyear,month,mrec,nrec
      logical there
      real*4  alatmin,alatmax,alonmin,alonmax,rlatinc,rloninc
      real*4  centerj,centeri
      integer ny,nx
      integer i,j,k,l,m,n
      integer number,nvar,n_month,nrecord,nday,mday
      real*4, allocatable ::  c(:,:,:),b(:,:),xlat(:,:),xlon(:,:)

      allocate(sigma(kz+1))
      allocate(b(jx,iy))
      allocate(xlat(jx-2,iy-2))
      allocate(xlon(jx-2,iy-2))


      open(10,file=trim(Path_Input)//trim(DomainName)//'.INFO'  &
             ,form='unformatted',recl=jx*iy*ibyte,access='direct')
      read(10,rec=1) iiy,jjx,kkz,dsinm,clat,clon,plat,plon,GRDFAC  &
                    ,iproj,(sigma(k),k=1,kz+1),ptop,igrads,ibigend &
                    ,truelatL,truelatH
      if(iiy.ne.iy.or.jjx.ne.jx.or.kkz.ne.kz) then
         write(*,*) 'iy,jx,kz in parameter = ',iy,jx,kz
         write(*,*) 'iy,jx,kz in DOMAIN.INFO ',iiy,jjx,kkz
         write(*,*) 'They are not consistent'
         stop
      endif
      read(10,rec=5) b
      do i=1,iy-2
      do j=1,jx-2
         xlat(j,i) = b(j+1,i+1)
      enddo
      enddo
      read(10,rec=6) b
      do i=1,iy-2
      do j=1,jx-2
         xlon(j,i) = b(j+1,i+1)
      enddo
      enddo
      close(10)

      deallocate(b)
      allocate(b(jx-2,iy-2))

      if(filename.eq.'ICBC') then
         nvar = kz*4+2
      else
         write(*,*) 'filename is not correct'
         stop
      endif
      allocate(c(jx-2,iy-2,nvar))

      if(idate1.gt.1940010100) then        ! original file
         n_month = (idate2/1000000-idate1/1000000)*12 &
                 + (mod(idate2/10000,100)-mod(idate1/10000,100))
         if(mod(idate2,10000).gt.0100) n_month = n_month+1
         ntype = 0
      else if(idate1.gt.19400101) then     ! daily mean file
         n_month = (idate2/10000-idate1/10000)*12 &
                 + (mod(idate2/100,100)-mod(idate1/100,100))
         if(mod(idate2,100).gt.01) n_month = n_month+1
         ntype = 1
      else
         write(*,*) 'date should be in 10 digits, or 8 digits'
         stop
      endif

      if(filename.eq.'ICBC') then
         fileout= '_ICBC.mon'
      endif

      open(20,file=trim(Path_Input)//trim(DomainName)//fileout  &
          ,form='unformatted',recl=(iy-2)*(jx-2)*ibyte,access='direct')
      nrec = 0

      inquire(file=trim(Path_Input)//trim(DomainName)//fileout//'.ctl' &
             ,exist=there)
      if(there) then
         open(31,file= &
             trim(Path_Input)//trim(DomainName)//fileout//'.ctl' &
                                                   ,status='replace')
      else
         open(31,file= &
       trim(Path_Input)//trim(DomainName)//fileout//'.ctl',status='new')
      endif
      write(31,10) '^'//trim(DomainName)//fileout
  10  format('dset ',A20)
      write(31,20)
  20  format('title RegCM daily input variables')
      write(31,30)
  30  format('options big_endian')
      write(31,50)
  50  format('undef -1.e34')
      if(iproj.eq.'LAMCON'.or.iproj.eq.'ROTMER') then
         alatmin= 999999.
         alatmax=-999999.
         do j=1,jx-2
            if(xlat(j,1).lt.alatmin) alatmin=xlat(j,1)
            if(xlat(j,iy-2).gt.alatmax) alatmax=xlat(j,iy-2)
         enddo
         alonmin= 999999.
         alonmax=-999999.
         do i=1,iy-2
         do j=1,jx-2
            if(clon.ge.0.0) then
               if(xlon(j,i).ge.0.0) then
                  alonmin = amin1(alonmin,xlon(j,i))
                  alonmax = amax1(alonmax,xlon(j,i))
               else
                  if(abs(clon-xlon(j,i)).lt. &
                     abs(clon-(xlon(j,i)+360.))) then
                     alonmin = amin1(alonmin,xlon(j,i))
                     alonmax = amax1(alonmax,xlon(j,i))
                  else
                     alonmin = amin1(alonmin,xlon(j,i)+360.)
                     alonmax = amax1(alonmax,xlon(j,i)+360.)
                  endif
               endif
            else
               if(xlon(j,i).lt.0.0) then
                  alonmin = amin1(alonmin,xlon(j,i))
                  alonmax = amax1(alonmax,xlon(j,i))
               else
                  if(abs(clon-xlon(j,i).lt. &
                     abs(clon-xlon(j,i)-360.))) then
                     alonmin = amin1(alonmin,xlon(j,i))
                     alonmax = amax1(alonmax,xlon(j,i))
                  else
                     alonmin = amin1(alonmin,xlon(j,i)-360.)
                     alonmax = amax1(alonmax,xlon(j,i)-360.)
                  endif
               endif
            endif
         enddo
         enddo
         rlatinc=dsinm*0.001/111./2.
         rloninc=dsinm*0.001/111./2.
         ny=2+nint(abs(alatmax-alatmin)/rlatinc)
         nx=1+nint(abs((alonmax-alonmin)/rloninc))
         centerj=(jx-2)/2.
         centeri=(iy-2)/2.
      endif

      if(iproj.eq.'LAMCON') then        ! Lambert projection
         write(31,100) jx-2,iy-2,clat,clon,centerj,centeri, &
                       truelatL,truelatH,clon,dsinm,dsinm
 100  format('pdef ',i4,1x,i4,1x,'lccr',7(1x,f7.2),1x,2(f7.0,1x))
         write(31,110) nx+2,alonmin-rloninc,rloninc
 110  format('xdef ',i4,' linear ',f7.2,1x,f7.4)
         write(31,120) ny+2,alatmin-rlatinc,rlatinc
 120  format('ydef ',i4,' linear ',f7.2,1x,f7.4)
      elseif(iproj.eq.'POLSTR') then    !
      elseif(iproj.eq.'NORMER') then
         write(31,200)  jx-2,xlon(1,1),xlon(2,1)-xlon(1,1)
 200  format('xdef ',I3,' linear ',f9.4,' ',f9.4)
         write(31,210) iy-2
 210  format('ydef ',I3,' levels')
         write(31,220) (xlat(1,i),i=1,iy-2)
 220  format(10f7.2)
      elseif(iproj.eq.'ROTMER') then
         write(*,*) 'Note that rotated Mercartor (ROTMER)' &
                   ,' projections are not supported by GrADS.'
         write(*,*) '  Although not exact, the eta.u projection' &
                   ,' in GrADS is somewhat similar.'
         write(*,*) ' FERRET, however, does support this projection.'
         write(31,230) jx-2,iy-2,plon,plat,dsinm/111000. &
                                          ,dsinm/111000.*.95238
 230  format('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
         write(31,110) nx+2,alonmin-rloninc,rloninc
         write(31,120) ny+2,alatmin-rlatinc,rlatinc
      else
         write(*,*) 'Are you sure your map projection is right ?'
         stop
      endif
      write(31,318) kz,((1013.25-ptop*10.)*(sigma(k)+sigma(k+1))*.5 &
                       +ptop*10.,k=kz,1,-1)
 318  format('zdef ',I2,' levels ',30f7.2)
      write(31,400)n_month,cday(16),chmc(month),nyear
 400  format('tdef ',I4,' linear ',A2,A3,I4,' ','1mo')
      write(31,500) 4+2
 500  format('vars ',I3)
 600  format(A8,'0 99 ',A36)
 611  format(A8,'0 33,105 ',A36)
 612  format(A8,'0 34,105 ',A36)
      if(filename.eq.'ICBC') then
 660  format(A8,I2,' 0 ',A36)
 661  format(A8,I2,' 33,100 ',A36)
 662  format(A8,I2,' 34,100 ',A36)
        if(iproj.eq.'LAMCON') then
           write(31,661) 'u       ',kz,'westerly wind (m/s)         '
           write(31,662) 'v       ',kz,'southerly wind (m/s)        '
        else
           write(31,660) 'u       ',kz,'westerly wind (m/s)         '
           write(31,660) 'v       ',kz,'southerly wind (m/s)        '
        endif
        write(31,660) 't       ',kz,'air temperature (degree)    '
        write(31,660) 'q       ',kz,'specific moisture (kg/kg)   '
        write(31,600) 'psa     ',   'surface pressure (hPa)      '
        write(31,600) 'ts      ',   'surface air temperature     '
      endif
      write(31,700)
 700  format('endvars')
      close(31)

      do nfile=1,n_month
         nrec = 0
         if(ntype.eq.0) then
            nyear = idate1/1000000
            month = idate1/10000-nyear*100 +nfile-1
         else if(ntype.eq.1) then
            nyear = idate1/10000
            month = idate1/100-nyear*100 +nfile-1
         endif

         nyear = nyear + (month-1)/12
         month = mod(month,12)
         if(month.eq.0) month = 12

         if(month.eq.1.or.month.eq.3.or.month.eq.5.or.  &
            month.eq.7.or.month.eq.8.or.month.eq.10.or. &
            month.eq.12) then
            nrecord = 31
         else if(month.eq.4.or.month.eq.6.or.month.eq.9.or. &
                 month.eq.11) then
            nrecord = 30
         else
            nrecord = 28
            if(mod(nyear,4).eq.0) nrecord = 29
            if(mod(nyear,100).eq.0) nrecord = 28
            if(mod(nyear,400).eq.0) nrecord = 29
         endif
         if(ntype.eq.0) then
            if(nfile.eq.n_month.and.mod(idate2/100-1,100).ne.0) &
               nrecord = min(nrecord,mod(idate2/100-1,100))
         else if(ntype.eq.1) then
            if(nfile.eq.n_month.and.mod(idate2-1,100).ne.0) &
               nrecord = min(nrecord,mod(idate2-1,100)+1)
         endif

         if(filename.eq.'ICBC') then
            filein = '_ICBC'//chy(nyear)//chm(month)//'01'
         endif
         open(10,file=trim(Path_Input)//trim(DomainName)//filein &
           ,form='unformatted',recl=(iy-2)*(jx-2)*ibyte,access='direct')

         if(nfile.eq.1.or.month.eq.1) then
         write(*,*)'Calculate the monthly mean of ',filename,nyear,month
         else
         write(*,*)'                              ',filename,nyear,month
         endif

         do n=1,nvar
            do i=1,iy-2
            do j=1,jx-2
               c(j,i,n) = 0.0
            enddo
            enddo
         enddo

         do nday=1,nrecord
            do n=1,nvar
               nrec=nrec+1
               read(10,rec=nrec) b
               do i=1,iy-2
               do j=1,jx-2
                  c(j,i,n) = c(j,i,n)+b(j,i)
               enddo
               enddo
            enddo
         enddo
         do n=1,nvar
            do i=1,iy-2
            do j=1,jx-2
               c(j,i,n) = c(j,i,n)/float(nrecord)
            enddo
            enddo
            mrec=mrec+1
            write(20,rec=mrec) ((c(j,i,n),j=1,jx-2),i=1,iy-2)
         enddo
         close(10)
      enddo
      close(20)
      return
      end

      subroutine mon3(filename,iy,jx,kz,nsg, Path_Input,DomainName &
                     ,Path_Output,ibyte,idate0,idate1,idate2)
      implicit none
      character*3 filename
      character*80 Path_Input,Path_Output
      character*10 DomainName
      integer iy,jx,kz,nsg,ibyte,idate0,idate1,idate2
      integer iiy,jjx,kkz
      integer mdate0,ibltyp,icup,ipptls,iboudy
      real*4  truelatL,truelatH
      real*4  dxsp,ptsp,clat,clon,plat,plon
      real*4  dto,dtb,dtr,dtc
      integer iotyp
      character*6 iproj
      real*4, allocatable ::  sigma(:)
      character*4 :: chy(1941:2100)
      data chy/'1941','1942','1943','1944','1945','1946','1947','1948',&
        '1949','1950','1951','1952','1953','1954','1955','1956','1957',&
        '1958','1959','1960','1961','1962','1963','1964','1965','1966',&
        '1967','1968','1969','1970','1971','1972','1973','1974','1975',&
        '1976','1977','1978','1979','1980','1981','1982','1983','1984',&
        '1985','1986','1987','1988','1989','1990','1991','1992','1993',&
        '1994','1995','1996','1997','1998','1999','2000','2001','2002',&
        '2003','2004','2005','2006','2007','2008','2009','2010','2011',&
        '2012','2013','2014','2015','2016','2017','2018','2019','2020',&
        '2021','2022','2023','2024','2025','2026','2027','2028','2029',&
        '2030','2031','2032','2033','2034','2035','2036','2037','2038',&
        '2039','2040','2041','2042','2043','2044','2045','2046','2047',&
        '2048','2049','2050','2051','2052','2053','2054','2055','2056',&
        '2057','2058','2059','2060','2061','2062','2063','2064','2065',&
        '2066','2067','2068','2069','2070','2071','2072','2073','2074',&
        '2075','2076','2077','2078','2079','2080','2081','2082','2083',&
        '2084','2085','2086','2087','2088','2089','2090','2091','2092',&
        '2093','2094','2095','2096','2097','2098','2099','2100'/
      character*2 cday(31)
      data cday/'01','02','03','04','05','06','07','08','09','10', &
                '11','12','13','14','15','16','17','18','19','20', &
                '21','22','23','24','25','26','27','28','29','30','31'/
      character*2 chm(12)
      data chm/'01','02','03','04','05','06','07','08','09','10', &
               '11','12'/
      character*3 chmc(12)
      data chmc/'jan','feb','mar','apr','may','jun'  &
               ,'jul','aug','sep','oct','nov','dec'/
      character*12 filein
      character*7 fileout
      character*14 headfile
      integer ntype,nfile,nyear,month,mrec,nrec
      integer i,j,k,l,m,n,itr
      integer number,nvar,n_month,nrecord,nday,mday
      logical there
      real*4  alatmin,alatmax,alonmin,alonmax,rlatinc,rloninc
      real*4  centerj,centeri
      integer ny,nx

      real*4, allocatable ::  c(:,:,:),b(:,:),xlat(:,:),xlon(:,:)

      allocate(sigma(kz+1))
      allocate(b((jx-2)*nsg,(iy-2)*nsg))
      allocate(xlat((jx-2)*nsg,(iy-2)*nsg))
      allocate(xlon((jx-2)*nsg,(iy-2)*nsg))

      if(nsg.lt.10) then
         write(headfile,101) nsg
      else
         write(headfile,102) nsg
      endif
 101  format('.INFO',I1)
 102  format('.INFO',I2)

      open(10,file='OUT_HEAD',form='unformatted' &
             ,recl=(jx-2)*(iy-2)*ibyte,access='direct')
      read(10,rec=1) mdate0,ibltyp,icup,ipptls,iboudy  &
                    ,iiy,jjx,kkz,(sigma(k),k=1,kz+1)   &
                    ,dxsp,ptsp,clat,clon,plat,plon     &
                    ,iproj,dto,dtb,dtr,dtc,iotyp,truelatL,truelatH
      if(iiy.ne.iy.or.jjx.ne.jx.or.kkz.ne.kz) then
         write(*,*) 'iy,jx,kz in parameter = ',iy,jx,kz
         write(*,*) 'iy,jx,kz in OUT_HEAD ',iiy,jjx,kkz
         write(*,*) 'They are not consistent'
         stop
      endif
      read(10,rec=6) xlat
      read(10,rec=7) xlon
      close(10)

      if(filename.eq.'SUB') then
         nvar = 16
      else
         write(*,*) 'filename is not correct'
         stop
      endif
      allocate(c((jx-2)*nsg,(iy-2)*nsg,nvar))

      if(idate1.gt.1940010100) then        ! original file
         n_month = (idate2/1000000-idate1/1000000)*12 &
                 + (mod(idate2/10000,100)-mod(idate1/10000,100))
         if(mod(idate2,10000).gt.0100) n_month = n_month+1
         ntype = 0
      else if(idate1.gt.19400101) then     ! daily mean file
         n_month = (idate2/10000-idate1/10000)*12 &
                 + (mod(idate2/100,100)-mod(idate1/100,100))
         if(mod(idate2,100).gt.01) n_month = n_month+1
         ntype = 1
      else
         write(*,*) 'date should be in 10 digits, or 8 digits'
         stop
      endif

      if(filename.eq.'SUB') then
         fileout= 'SUB.mon'
      endif

      open(20,file=trim(Path_Output)//fileout,form='unformatted' &
             ,recl=(iy-2)*nsg*(jx-2)*nsg*ibyte,access='direct')
      mrec = 0

      inquire(file=trim(Path_Output)//fileout//'.ctl',exist=there)
      if(there) then
       open(31,file=trim(Path_Output)//fileout//'.ctl',status='replace')
      else
         open(31,file=trim(Path_Output)//fileout//'.ctl',status='new')
      endif
      write(31,10) fileout
  10  format('dset ^',A7)
      write(31,20)
  20  format('title RegCM monthly output variables')
      write(31,30)
  30  format('options big_endian')
      write(31,50)
  50  format('undef -1.e34')
      if(iproj.eq.'LAMCON'.or.iproj.eq.'ROTMER') then
         alatmin= 999999.
         alatmax=-999999.
         do j=1,jx-2
            if(xlat(j,1).lt.alatmin) alatmin=xlat(j,1)
            if(xlat(j,iy-2).gt.alatmax) alatmax=xlat(j,iy-2)
         enddo
         alonmin= 999999.
         alonmax=-999999.
         do i=1,iy-2
         do j=1,jx-2
            if(clon.ge.0.0) then
               if(xlon(j,i).ge.0.0) then
                  alonmin = amin1(alonmin,xlon(j,i))
                  alonmax = amax1(alonmax,xlon(j,i))
               else
                  if(abs(clon-xlon(j,i)).lt. &
                     abs(clon-(xlon(j,i)+360.))) then
                     alonmin = amin1(alonmin,xlon(j,i))
                     alonmax = amax1(alonmax,xlon(j,i))
                  else
                     alonmin = amin1(alonmin,xlon(j,i)+360.)
                     alonmax = amax1(alonmax,xlon(j,i)+360.)
                  endif
               endif
            else
               if(xlon(j,i).lt.0.0) then
                  alonmin = amin1(alonmin,xlon(j,i))
                  alonmax = amax1(alonmax,xlon(j,i))
               else
                  if(abs(clon-xlon(j,i).lt. &
                     abs(clon-xlon(j,i)-360.))) then
                     alonmin = amin1(alonmin,xlon(j,i))
                     alonmax = amax1(alonmax,xlon(j,i))
                  else
                     alonmin = amin1(alonmin,xlon(j,i)-360.)
                     alonmax = amax1(alonmax,xlon(j,i)-360.)
                  endif
               endif
            endif
         enddo
         enddo
         rlatinc=dxsp/111./2.
         rloninc=dxsp/111./2.
         ny=2+nint(abs(alatmax-alatmin)/rlatinc)
         nx=1+nint(abs((alonmax-alonmin)/rloninc))
         centerj=(jx-2)/2.
         centeri=(iy-2)/2.
      endif

      if(iproj.eq.'LAMCON') then        ! Lambert projection
         write(31,100) jx-2,iy-2,clat,clon,centerj,centeri, &
                       truelatL,truelatH,clon,dxsp*1000.,dxsp*1000.
 100  format('pdef ',i4,1x,i4,1x,'lccr',7(1x,f7.2),1x,2(f7.0,1x))
         write(31,110) nx+2,alonmin-rloninc,rloninc
 110  format('xdef ',i4,' linear ',f7.2,1x,f7.4)
         write(31,120) ny+2,alatmin-rlatinc,rlatinc
 120  format('ydef ',i4,' linear ',f7.2,1x,f7.4)
      elseif(iproj.eq.'POLSTR') then    !
      elseif(iproj.eq.'NORMER') then
         write(31,200)  jx-2,xlon(1,1),xlon(2,1)-xlon(1,1)
 200  format('xdef ',I3,' linear ',f9.4,' ',f9.4)
         write(31,210) iy-2
 210  format('ydef ',I3,' levels')
         write(31,220) (xlat(1,i),i=1,iy-2)
 220  format(10f7.2)
      elseif(iproj.eq.'ROTMER') then
         write(*,*) 'Note that rotated Mercartor (ROTMER)' &
                   ,' projections are not supported by GrADS.'
         write(*,*) '  Although not exact, the eta.u projection' &
                   ,' in GrADS is somewhat similar.'
         write(*,*) ' FERRET, however, does support this projection.'
         write(31,230) jx-2,iy-2,plon,plat,dxsp/111. &
                                          ,dxsp/111.*.95238
 230  format('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
         write(31,110) nx+2,alonmin-rloninc,rloninc
         write(31,120) ny+2,alatmin-rlatinc,rlatinc
      else
         write(*,*) 'Are you sure your map projection is right ?'
         stop
      endif
      write(31,300) (1013.25-ptsp*10.)*(sigma(1)+sigma(2))*.5+ptsp*10.
 300  format('zdef 1',' levels ',f7.2)
      nyear=idate0/1000000
      month=(idate0-nyear*1000000)/10000
      nday =(idate0-nyear*1000000-month*10000)/100
      write(31,400)n_month,cday(16),chmc(month),nyear
 400  format('tdef ',I4,' linear ',A2,A3,I4,' ','1mo')
      write(31,500) nvar
 500  format('vars ',I3)
 600  format(A8,'0 99 ',A36)
 611  format(A8,'0 33,105 ',A36)
 612  format(A8,'0 34,105 ',A36)
      if(filename.eq.'SUB') then
       if(iproj.eq.'LAMCON') then        ! Lambert projection
        write(31,611) 'u10m    ','westerly  wind at 10m (m/s)          '
        write(31,612) 'v10m    ','southerly wind at 10m (m/s)          '
       else
        write(31,600) 'u10m    ','westerly  wind at 10m (m/s)          '
        write(31,600) 'v10m    ','southerly wind at 10m (m/s)          '
       endif
        write(31,600) 'uvdrag  ','surface drag stress                  '
        write(31,600) 'tg      ','ground temperature (degree)          '
        write(31,600) 'tlef    ','temperature of foliage               '
        write(31,600) 't2m     ','air temperature at 2m (K)            '
        write(31,600) 'q2m     ','water vapor mixing ratio at 2m(kg/kg)'
        write(31,600) 'ssw     ','upper layer soil water               '
        write(31,600) 'rsw     ','root zone soil water                 '
        write(31,600) 'tpr     ','total precipitation (mm/day)         '
        write(31,600) 'evp     ','evapotranspiration (mm/day)          '
        write(31,600) 'runoff  ','surface runoff (mm/day)              '
        write(31,600) 'scv     ','total snow amount                    '
        write(31,600) 'sena    ','sensible heat flux (W/m2)            '
        write(31,600) 'prcv    ','convective precipitation (mm/day)    '
        write(31,600) 'psb     ','surface pressure (hPa)               '
        write(31,700)
 700  format('endvars')
      endif
      close(31)

      do nfile=1,n_month
         if(ntype.eq.0) then
            nyear = idate1/1000000
            month = idate1/10000-nyear*100 +nfile-1
         else if(ntype.eq.1) then
            nyear = idate1/10000
            month = idate1/100-nyear*100 +nfile-1
         endif

         nyear = nyear + (month-1)/12
         month = mod(month,12)
         if(month.eq.0) month = 12

         if(nfile.eq.1.or.month.eq.1) then
         write(*,*)'Calculate the monthly mean of ',filename,nyear,month
         else
         write(*,*)'                              ',filename,nyear,month
         endif

         if(month.eq.1.or.month.eq.3.or.month.eq.5.or.  &
            month.eq.7.or.month.eq.8.or.month.eq.10.or. &
            month.eq.12) then
            nrecord = 31
         else if(month.eq.4.or.month.eq.6.or.month.eq.9.or. &
                 month.eq.11) then
            nrecord = 30
         else
            nrecord = 28
            if(mod(nyear,4).eq.0) nrecord = 29
            if(mod(nyear,100).eq.0) nrecord = 28
            if(mod(nyear,400).eq.0) nrecord = 29
         endif
         if(ntype.eq.0) then
            if(nfile.eq.n_month.and.mod(idate2/100-1,100).ne.0) &
               nrecord = min(nrecord,mod(idate2/100-1,100))
         else if(ntype.eq.1) then
            if(nfile.eq.n_month.and.mod(idate2-1,100).ne.0) &
               nrecord = min(nrecord,mod(idate2-1,100)+1)
         endif

         if(filename.eq.'SUB') then
            filein = 'SUB.'//chy(nyear)//chm(month)//'01'
         endif
         open(10,file=trim(Path_Output)//filein,form='unformatted' &
                ,recl=(iy-2)*nsg*(jx-2)*nsg*ibyte,access='direct')
         nrec = 0

         do n=1,nvar
            do i=1,(iy-2)*nsg
            do j=1,(jx-2)*nsg
               c(j,i,n) = 0.0
            enddo
            enddo
         enddo

         do nday=1,nrecord
            do n=1,nvar
               nrec=nrec+1
               read(10,rec=nrec) b
               do i=1,(iy-2)*nsg
               do j=1,(jx-2)*nsg
                  c(j,i,n) = c(j,i,n)+b(j,i)
               enddo
               enddo
            enddo
         enddo
         do n=1,nvar
            do i=1,(iy-2)*nsg
            do j=1,(jx-2)*nsg
               c(j,i,n) = c(j,i,n)/float(nrecord)
            enddo
            enddo
            mrec=mrec+1
            write(20,rec=mrec)((c(j,i,n),j=1,(jx-2)*nsg),i=1,(iy-2)*nsg)
         enddo
         close(10)
      enddo
      close(20)
      deallocate(sigma)
      deallocate(b)
      deallocate(xlat)
      deallocate(xlon)
      deallocate(c)
      return
      end

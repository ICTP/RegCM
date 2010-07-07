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
      integer iy,jx,kz,nsg,ntr,i_band,ibyte,igrads
      logical day_ICBC,day_ATM,day_RAD,day_SRF,day_SUB,day_CHE
      logical day_ICBC_P,day_ATM_P
      logical mon_ICBC,mon_ATM,mon_RAD,mon_SRF,mon_SUB,mon_CHE
      logical mon_ICBC_P,mon_ATM_P
      integer idate0,idate1,idate2
      character*128 Path_Input,Path_Output
      character*20 DomainName
      integer len_path
      integer np
      real*4, save ::  plev(11)
      logical s2p_ICBC,s2p_ATM

      namelist /shareparam/iy,jx,kz,nsg,ntr,i_band,ibyte,igrads &
                          ,Path_Input,DomainName,Path_Output
      namelist /sigma2p_param/ np,plev,s2p_ICBC,s2p_ATM
      namelist /dateparam/ idate0,idate1,idate2
      namelist /averageparam/ &
                day_ICBC,day_ATM,day_RAD,day_SRF,day_SUB,day_CHE &
               ,day_ICBC_P,day_ATM_P &
               ,mon_ICBC,mon_ATM,mon_RAD,mon_SRF,mon_SUB,mon_CHE &
               ,mon_ICBC_P,mon_ATM_P

      i_band = 0
      read(*,shareparam)

      len_path = len(trim(Path_Input))
      if(Path_Input(len_path:len_path).ne.'/') &
         Path_Input=trim(Path_Input)//'/'
      len_path = len(trim(Path_Output))
      if(Path_Output(len_path:len_path).ne.'/') &
         Path_Output=trim(Path_Output)//'/'

      read(*,sigma2p_param)
      if(np.ne.11) then
         write(*,*) 'Number of pressure levels does not equal to 11'
         write(*,*) 'Please reset np, or reset plev in average.f90'
         stop
      endif

      read(*,dateparam)

      read(*,averageparam)

      if(day_ICBC) call day2('ICBC',iy,jx,kz,i_band,ibyte, Path_Input &
                         ,DomainName,idate0,idate1,idate2,igrads)
      if(day_ICBC_P) call day2p('ICBC_P',iy,jx,kz,np,plev,i_band,ibyte &
                     ,Path_Input,DomainName,idate0,idate1,idate2,igrads)

      if(day_ATM)call day('ATM',iy,jx,kz,ntr,i_band,ibyte, Path_Output &
                          ,idate0,idate1,idate2,igrads)
      if(day_RAD)call day('RAD',iy,jx,kz,ntr,i_band,ibyte, Path_Output &
                          ,idate0,idate1,idate2,igrads)
      if(day_SRF)call day('SRF',iy,jx,kz,ntr,i_band,ibyte, Path_Output &
                          ,idate0,idate1,idate2,igrads)
      if(nsg.gt.1.and.day_SUB) &
                call day3('SUB',iy,jx,kz,nsg,i_band,ibyte, Path_Input &
                    ,DomainName,Path_Output,idate0,idate1,idate2,igrads)
      if(day_CHE)call day('CHE',iy,jx,kz,ntr,i_band,ibyte, Path_Output &
                         ,idate0,idate1,idate2,igrads)

      if(day_ATM_P)call dayp('ATM_P',iy,jx,kz,np,plev,i_band,ibyte  &
                        ,Path_Output,idate0,idate1,idate2,igrads)

      if(mon_ICBC) call mon2('ICBC',iy,jx,kz,i_band,ibyte, Path_Input &
                            ,DomainName,idate1,idate2,igrads)
      if(mon_ICBC_P)call mon2p('ICBC_P',iy,jx,kz,np,plev,i_band,ibyte  &
                     ,Path_Input,DomainName,idate1,idate2,igrads)

      if(mon_ATM)call mon('ATM',iy,jx,kz,ntr,i_band,ibyte, Path_Output &
                         ,idate1,idate2,igrads)
      if(mon_RAD)call mon('RAD',iy,jx,kz,ntr,i_band,ibyte, Path_Output &
                         ,idate1,idate2,igrads)
      if(mon_SRF)call mon('SRF',iy,jx,kz,ntr,i_band,ibyte, Path_Output &
                         ,idate1,idate2,igrads)
      if(nsg.gt.1.and.mon_SUB) &
                call mon3('SUB',iy,jx,kz,nsg,i_band,ibyte, Path_Input &
                   ,DomainName,Path_Output,idate1,idate2,igrads)
      if(mon_CHE)call mon('CHE',iy,jx,kz,ntr,i_band,ibyte, Path_Output &
                         ,idate1,idate2,igrads)

      if(mon_ATM_P)call monp('ATM_P',iy,jx,kz,np,plev,i_band,ibyte  &
                           ,Path_Output,idate1,idate2,igrads)

      stop
      end

      subroutine day(filename,iy,jx,kz,ntr,i_band,ibyte, Path_Output &
                    ,idate0,idate1,idate2,igrads)
      implicit none
      character*3 filename
      character*128 Path_Output
      integer iy,jx,kz,ntr,i_band,ibyte,idate0,idate1,idate2,igrads
      integer iiy,jjx,kkz
      integer mdate0,ibltyp,icup,ipptls,iboudy
      real*4  truelatL,truelatH
      real*4  dxsp,ptsp,clat,clon,plat,plon
      real*4  dto,dtb,dtr,dtc,dt
      integer iotyp
      character*6 iproj
      real*4, allocatable ::  sigma(:)
      character*4 :: chy
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
      integer i,j,k,n,itr,jx_len
      integer nvar,mvar,n_month,nrecord,nday,mday
      integer idatex
      real*4, allocatable ::  c(:,:,:),b(:,:),xlat(:,:),xlon(:,:)

      if(i_band.eq.1) then
         jx_len = jx
      else
         jx_len = jx-2
      endif

      allocate(sigma(kz+1))
      allocate(b(jx_len,iy-2))
      allocate(xlat(jx_len,iy-2))
      allocate(xlon(jx_len,iy-2))

      inquire(file=trim(Path_Output)//'OUT_HEAD',exist=there)
      if(.not.there) then
         write(*,*) trim(Path_Output)//'OUT_HEAD',' is not avaiable'
         stop
      endif
      open(10,file=trim(Path_Output)//'OUT_HEAD',form='unformatted' &
             ,recl=jx_len*(iy-2)*ibyte,access='direct')
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
         nvar = kz*4+10
      else if(filename.eq.'CHE') then
         dt = dtc
         nvar = ntr*kz+kz*3+ntr*7+3
!        nvar = ntr*kz+kz*3+ntr*7+2      ! for RegCM3, one record less
      else
         write(*,*) 'filename is not correct'
         stop
      endif
      allocate(c(jx_len,iy-2,nvar))
      mvar = nvar
      if(filename.eq.'SRF') mvar = mvar -6

      if(idate1.ge.1000010100) then        ! original file
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
      write(*,*) 'Calculate the daily mean of    ',filename,nyear,month
         else
      write(*,*) '                               ',filename,nyear,month
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
         n_slice = nint(24./dt)

         write(chy,199) nyear
         if(filename.eq.'ATM') then
            filein = 'ATM.'//chy//chm(month)//'0100'
            fileout= 'ATM.'//chy//chm(month)//'01'
         else if(filename.eq.'SRF') then
            filein = 'SRF.'//chy//chm(month)//'0100'
            fileout= 'SRF.'//chy//chm(month)//'01'
         else if(filename.eq.'RAD') then
            filein = 'RAD.'//chy//chm(month)//'0100'
            fileout= 'RAD.'//chy//chm(month)//'01'
         else if(filename.eq.'CHE') then
            filein = 'CHE.'//chy//chm(month)//'0100'
            fileout= 'CHE.'//chy//chm(month)//'01'
         endif

         IF(igrads.eq.1) THEN
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
            do j=1,jx_len
               if(xlat(j,1).lt.alatmin) alatmin=xlat(j,1)
               if(xlat(j,iy-2).gt.alatmax) alatmax=xlat(j,iy-2)
            enddo
            alonmin= 999999.
            alonmax=-999999.
            do i=1,iy-2
            do j=1,jx_len
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
                     if(abs(clon-xlon(j,i)).lt. &
                        abs(clon-(xlon(j,i)-360.))) then
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
            centerj=jx_len/2.
            centeri=(iy-2)/2.
         endif

         if(iproj.eq.'LAMCON') then        ! Lambert projection
            write(31,100) jx_len,iy-2,clat,clon,centerj,centeri, &
                          truelatL,truelatH,clon,dxsp*1000.,dxsp*1000.
 100  format('pdef ',i4,1x,i4,1x,'lccr',7(1x,f7.2),1x,2(f7.0,1x))
            write(31,110) nx+2,alonmin-rloninc,rloninc
 110  format('xdef ',i4,' linear ',f7.2,1x,f7.4)
            write(31,120) ny+2,alatmin-rlatinc,rlatinc
 120  format('ydef ',i4,' linear ',f7.2,1x,f7.4)
         elseif(iproj.eq.'POLSTR') then    !
         elseif(iproj.eq.'NORMER') then
            write(31,200)  jx_len,xlon(1,1),xlon(2,1)-xlon(1,1)
 200  format('xdef ',I3,' linear ',f9.4,' ',f9.4)
            write(31,210) iy-2
 210  format('ydef ',I3,' levels')
            write(31,220) (xlat(1,i),i=1,iy-2)
 220  format(10f7.2)
         elseif(iproj.eq.'ROTMER') then
         if(nfile.eq.1.or.month.eq.1) then
            write(*,*) 'Note that rotated Mercartor (ROTMER)' &
                   ,' projections are not supported by GrADS.'
            write(*,*) '  Although not exact, the eta.u projection' &
                   ,' in GrADS is somewhat similar.'
            write(*,*) ' FERRET, however, does support this projection.'
         endif
            write(31,230) jx_len,iy-2,plon,plat,dxsp/111. &
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
         if(nfile.eq.1.and.idate1.eq.idate0) then
            nday=mod(idate1,10000)/100
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
         write(31,500) 4+10
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
        write(31,600) 'scv     ','snow amount (mm, water equivalent)   '
        write(31,600) 'sena    ','sensible heat flux (W/m2)            '
        write(31,600) 'flw     ','net infrared energy flux (W/m2)      '
        write(31,600) 'fsw     ','net absorbed solar energy flux (W/m2)'
        write(31,600) 'flwd    ','downward infrared energy flux (W/m2) '
        write(31,600) 'sina    ','incident solar energy flux (W/m2)    '
        write(31,600) 'prcv    ','convective precipitation (mm/day)    '
        write(31,600) 'ps      ','surface pressure (hPa)               '
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
        write(31,600) 'ps      ',   'surface pressure (hPa)      '
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
        write(31,600)'ps    ',   'surface pressure (hPa)               '
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
       write(31,600) 'ps      ','surface pressure (hPa)               '
      endif
      write(31,700)
 700  format('endvars')
      close(31)
         ENDIF

         inquire(file=trim(Path_Output)//filein,exist=there)
         if(.not.there) then
            write(*,*) trim(Path_Output)//filein,' is not avaiable'
            stop
         endif
         if(iotyp.eq.1) then
            open(10,file=trim(Path_Output)//filein,form='unformatted' &
                   ,recl=(iy-2)*jx_len*ibyte,access='direct')
         else if(iotyp.eq.2) then
            open(10,file=trim(Path_Output)//filein,form='unformatted')
         endif
         open(20,file=trim(Path_Output)//fileout,form='unformatted' &
                ,recl=(iy-2)*jx_len*ibyte,access='direct')
         mrec = 0
         do nday=1,nrecord
            do n=1,mvar
               do i=1,iy-2
               do j=1,jx_len
                  c(j,i,n) = 0.0
               enddo
               enddo
            enddo
            if(iotyp.eq.2) read(10) idatex
            if(filename.eq.'SRF') then
               do i=1,iy-2
               do j=1,jx_len
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
                  do j=1,jx_len
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
                  do j=1,jx_len
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
                  do j=1,jx_len
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
                  do j=1,jx_len
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
                  do j=1,jx_len
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
                  do j=1,jx_len
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
                  do j=1,jx_len
                     if(c(j,i,nvar).gt.b(j,i)) &
                        c(j,i,nvar) =  b(j,i)
                  enddo
                  enddo
               endif
            enddo
            do n=1,mvar
               do i=1,iy-2
               do j=1,jx_len
                  c(j,i,n) = c(j,i,n)/float(n_slice)
               enddo
               enddo
               mrec=mrec+1
               write(20,rec=mrec) ((c(j,i,n),j=1,jx_len),i=1,iy-2)
            enddo
            if(filename.eq.'SRF') then
               do n=nvar-5,nvar
                  mrec=mrec+1
                  write(20,rec=mrec) ((c(j,i,n),j=1,jx_len),i=1,iy-2)
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
 199  format(I4)
      return
      end

      subroutine dayp(filename,iy,jx,kz,np,plev,i_band,ibyte  &
                     ,Path_Output,idate0,idate1,idate2,igrads)
      implicit none
      character*5 filename
      character*128 Path_Output
      integer iy,jx,kz,np,i_band,ibyte,idate0,idate1,idate2,igrads
      real*4  plev(np)
      integer iiy,jjx,kkz
      integer mdate0,ibltyp,icup,ipptls,iboudy
      real*4  truelatL,truelatH
      real*4  dxsp,ptsp,clat,clon,plat,plon
      real*4  dto,dtb,dtr,dtc
      integer iotyp
      character*6 iproj
      real*4, allocatable ::  sigma(:)
      character*4 :: chy
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
      character*16 filein
      character*14 fileout
      integer ntype,nfile,nyear,month,n_slice,mrec,nrec
      logical there
      real*4  alatmin,alatmax,alonmin,alonmax,rlatinc,rloninc
      real*4  centerj,centeri
      integer ny,nx
      integer i,j,k,n,jx_len
      integer nvar,n_month,nrecord,nday,mday
      real*4, allocatable ::  c(:,:,:),b(:,:),xlat(:,:),xlon(:,:)

      if(i_band.eq.1) then
         jx_len = jx
      else
         jx_len = jx-2
      endif

      allocate(sigma(kz+1))
      allocate(b(jx_len,iy-2))
      allocate(xlat(jx_len,iy-2))
      allocate(xlon(jx_len,iy-2))

      inquire(file=trim(Path_Output)//'OUT_HEAD',exist=there)
      if(.not.there) then
         write(*,*) trim(Path_Output)//'OUT_HEAD',' is not avaiable'
         stop
      endif
      open(10,file=trim(Path_Output)//'OUT_HEAD',form='unformatted' &
             ,recl=jx_len*(iy-2)*ibyte,access='direct')
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

      if(filename.eq.'ATM_P') then
         nvar = np*8+6
      else
         write(*,*) 'filename is not correct'
         stop
      endif
      allocate(c(jx_len,iy-2,nvar))

      if(idate1.ge.1000010100) then        ! original file
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
         write(*,*) 'Calculate the daily mean of  ',filename,nyear,month
         else
         write(*,*) '                             ',filename,nyear,month
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
         n_slice = nint(24./dto)

         write(chy,199) nyear
         filein = 'ATM_P.'//chy//chm(month)//'0100'
         fileout= 'ATM_P.'//chy//chm(month)//'01'

         IF(igrads.eq.1) THEN
         inquire(file=trim(Path_Output)//fileout//'.ctl',exist=there)
         if(there) then
           open(31,file=trim(Path_Output)//fileout//'.ctl'  &
                  ,status='replace')
         else
           open(31,file=trim(Path_Output)//fileout//'.ctl',status='new')
         endif
         write(31,10) fileout
  10     format('dset ^',A14)
         write(31,20)
  20     format('title RegCM daily pressure level output variables')
         write(31,30)
  30     format('options big_endian')
         write(31,50)
  50     format('undef -1.e34')
         if(iproj.eq.'LAMCON'.or.iproj.eq.'ROTMER') then
            alatmin= 999999.
            alatmax=-999999.
            do j=1,jx_len
               if(xlat(j,1).lt.alatmin) alatmin=xlat(j,1)
               if(xlat(j,iy-2).gt.alatmax) alatmax=xlat(j,iy-2)
            enddo
            alonmin= 999999.
            alonmax=-999999.
            do i=1,iy-2
            do j=1,jx_len
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
                     if(abs(clon-xlon(j,i)).lt. &
                        abs(clon-(xlon(j,i)-360.))) then
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
            centerj=jx_len/2.
            centeri=(iy-2)/2.
         endif

         if(iproj.eq.'LAMCON') then        ! Lambert projection
            write(31,100) jx_len,iy-2,clat,clon,centerj,centeri, &
                          truelatL,truelatH,clon,dxsp*1000.,dxsp*1000.
 100  format('pdef ',i4,1x,i4,1x,'lccr',7(1x,f7.2),1x,2(f7.0,1x))
            write(31,110) nx+2,alonmin-rloninc,rloninc
 110  format('xdef ',i4,' linear ',f7.2,1x,f7.4)
            write(31,120) ny+2,alatmin-rlatinc,rlatinc
 120  format('ydef ',i4,' linear ',f7.2,1x,f7.4)
         elseif(iproj.eq.'POLSTR') then    !
         elseif(iproj.eq.'NORMER') then
            write(31,200)  jx_len,xlon(1,1),xlon(2,1)-xlon(1,1)
 200  format('xdef ',I3,' linear ',f9.4,' ',f9.4)
            write(31,210) iy-2
 210  format('ydef ',I3,' levels')
            write(31,220) (xlat(1,i),i=1,iy-2)
 220  format(10f7.2)
         elseif(iproj.eq.'ROTMER') then
         if(nfile.eq.1.or.month.eq.1) then
            write(*,*) 'Note that rotated Mercartor (ROTMER)' &
                   ,' projections are not supported by GrADS.'
            write(*,*) '  Although not exact, the eta.u projection' &
                   ,' in GrADS is somewhat similar.'
            write(*,*) ' FERRET, however, does support this projection.'
         endif
            write(31,230) jx_len,iy-2,plon,plat,dxsp/111. &
                                          ,dxsp/111.*.95238
 230  format('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
            write(31,110) nx+2,alonmin-rloninc,rloninc
            write(31,120) ny+2,alatmin-rlatinc,rlatinc
         else
            write(*,*) 'Are you sure your map projection is right ?'
            stop
         endif
      write(31,318) np,(plev(k),k=1,np)
 300  format('zdef 1',' levels ',f7.2)
 318  format('zdef ',I2,' levels ',30f7.2)
         if(nfile.eq.1.and.idate1.eq.idate0) then
            nday=mod(idate1,10000)/100
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
         write(31,500) 8+6
 500  format('vars ',I3)
 600  format(A8,'0 99 ',A36)
 660  format(A8,I2,' 0 ',A36)
 661  format(A8,I2,' 33,100 ',A36)
 662  format(A8,I2,' 34,100 ',A36)
        if(iproj.eq.'LAMCON') then
           write(31,661) 'u       ',np,'westerly wind (m/s)         '
           write(31,662) 'v       ',np,'southerly wind (m/s)        '
        else
           write(31,660) 'u       ',np,'westerly wind (m/s)         '
           write(31,660) 'v       ',np,'southerly wind (m/s)        '
        endif
        write(31,660) 'omega   ',np,'omega (hPa/s)   p-velocity  '
        write(31,660) 'h       ',np,'geopotential height (m)     '
        write(31,660) 't       ',np,'air temperature (deg, K)    '
        write(31,660) 'rh      ',np,'relative moisture (%)       '
        write(31,660) 'qv      ',np,'water vapor mixing ratio    '
        write(31,660) 'qc      ',np,'cloud water mixing ratio    '
        write(31,600) 'ps      ',   'surface pressure (hPa)      '
        write(31,600) 'slp     ',   'sea level pressure (hPa)    '
        write(31,600) 'tpr     ',   'total precipitation(mm/day) '
        write(31,600) 'tgb     ',   'lower groud temp. in BATS   '
        write(31,600) 'swt     ',   'total soil water in mm H2O  '
        write(31,600) 'rno     ',   'accumulated infiltration    '
        write(31,700)
 700  format('endvars')
        close(31)
         ENDIF

         inquire(file=trim(Path_Output)//filein,exist=there)
         if(.not.there) then
            write(*,*) trim(Path_Output)//filein,' is not avaiable'
            stop
         endif
         open(10,file=trim(Path_Output)//filein,form='unformatted' &
                ,recl=(iy-2)*jx_len*ibyte,access='direct')
         open(20,file=trim(Path_Output)//fileout,form='unformatted' &
                ,recl=(iy-2)*jx_len*ibyte,access='direct')
         mrec = 0
         do nday=1,nrecord
            do n=1,nvar
               do i=1,iy-2
               do j=1,jx_len
                  c(j,i,n) = 0.0
               enddo
               enddo
            enddo
            do mday = 1,n_slice
               do n=1,nvar
                  nrec=nrec+1
                  read(10,rec=nrec) b
                  do i=1,iy-2
                  do j=1,jx_len
                     c(j,i,n) = c(j,i,n)+b(j,i)
                  enddo
                  enddo
               enddo
            enddo
            do n=1,nvar
               do i=1,iy-2
               do j=1,jx_len
                  c(j,i,n) = c(j,i,n)/float(n_slice)
               enddo
               enddo
               mrec=mrec+1
               write(20,rec=mrec) ((c(j,i,n),j=1,jx_len),i=1,iy-2)
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
 199  format(I4)
      return
      end

      subroutine day2(filename,iy,jx,kz,i_band,ibyte, Path_Input &
                     ,DomainName,idate0,idate1,idate2,igrads)
      implicit none
      character*4 filename
      character*128 Path_Input
      character*20 DomainName
      integer iy,jx,kz,i_band,ibyte,idate0,idate1,idate2,igrads
      integer iiy,jjx,kkz
      real*4  truelatL,truelatH
      real*4  dsinm,ptop,clat,clon,plat,plon,GRDFAC
      integer jgrads,ibigend
      character*6 iproj
      real*4, allocatable ::  sigma(:)
      character*4 :: chy
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
      integer i,j,k,n,jx_len
      integer nvar,n_month,nrecord,nday,mday
      integer idatex
      real*4, allocatable ::  c(:,:,:),b(:,:),xlat(:,:),xlon(:,:)

      if(i_band.eq.1) then
         jx_len = jx
      else
         jx_len = jx-2
      endif

      allocate(sigma(kz+1))
      allocate(b(jx,iy))

      allocate(xlat(jx_len,iy-2))
      allocate(xlon(jx_len,iy-2))

      inquire(file=trim(Path_Input)//trim(DomainName)//'.INFO' &
             ,exist=there)
      if(.not.there) then
         write(*,*) trim(Path_Input)//trim(DomainName)//'.INFO' &
                   ,' is not avaiable'
         stop
      endif
      open(10,file=trim(Path_Input)//trim(DomainName)//'.INFO'  &
             ,form='unformatted',recl=jx*iy*ibyte,access='direct')
      read(10,rec=1) iiy,jjx,kkz,dsinm,clat,clon,plat,plon,GRDFAC  &
                    ,iproj,(sigma(k),k=1,kz+1),ptop,jgrads,ibigend &
                    ,truelatL,truelatH
      if(iiy.ne.iy.or.jjx.ne.jx.or.kkz.ne.kz) then
         write(*,*) 'iy,jx,kz in parameter = ',iy,jx,kz
         write(*,*) 'iy,jx,kz in DOMAIN.INFO ',iiy,jjx,kkz
         write(*,*) 'They are not consistent'
         stop
      endif
      read(10,rec=5) b
      do i=1,iy-2
      do j=1,jx_len
         if(i_band.eq.1) then
            xlat(j,i) = b(j,i+1)
         else
            xlat(j,i) = b(j+1,i+1)
         endif
      enddo
      enddo
      read(10,rec=6) b
      do i=1,iy-2
      do j=1,jx_len
         if(i_band.eq.1) then
            xlon(j,i) = b(j,i+1)
         else
            xlon(j,i) = b(j+1,i+1)
         endif
      enddo
      enddo
      close(10)

      if(filename.eq.'ICBC') then
         nvar = kz*4+2
      else
         write(*,*) 'filename is not correct'
         stop
      endif
      allocate(c(jx_len,iy-2,nvar))

      if(idate1.ge.1000010100) then        ! original file
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
       write(*,*) 'Calculate the daily mean of   ',filename,nyear,month
         else
       write(*,*) '                              ',filename,nyear,month
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

         write(chy,199) nyear
         if(filename.eq.'ICBC') then
            filein = 'ICBC'//chy//chm(month)//'0100'
            fileout= 'ICBC'//chy//chm(month)//'01'
         endif

         IF(igrads.eq.1) THEN
         inquire(file= &
        trim(Path_Input)//trim(DomainName)//'_'//fileout//'.ctl'  &
                ,exist=there)
         if(there) then
            open(31,file= &
         trim(Path_Input)//trim(DomainName)//'_'//fileout//'.ctl' &
                                                  ,status='replace')
         else
            open(31,file= &
         trim(Path_Input)//trim(DomainName)//'_'//fileout//'.ctl' &
                                                      ,status='new')
         endif
         write(31,10) '^'//trim(DomainName)//'_'//fileout
  10     format('dset ',A34)
         write(31,20)
  20     format('title RegCM daily input variables')
         write(31,30)
  30     format('options big_endian')
         write(31,50)
  50     format('undef -1.e34')
         if(iproj.eq.'LAMCON'.or.iproj.eq.'ROTMER') then
            alatmin= 999999.
            alatmax=-999999.
            do j=1,jx_len
               if(xlat(j,1).lt.alatmin) alatmin=xlat(j,1)
               if(xlat(j,iy-2).gt.alatmax) alatmax=xlat(j,iy-2)
            enddo
            alonmin= 999999.
            alonmax=-999999.
            do i=1,iy-2
            do j=1,jx_len
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
                     if(abs(clon-xlon(j,i)).lt. &
                        abs(clon-(xlon(j,i)-360.))) then
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
            centerj=jx_len/2.
            centeri=(iy-2)/2.
         endif

         if(iproj.eq.'LAMCON') then        ! Lambert projection
            write(31,100) jx_len,iy-2,clat,clon,centerj,centeri, &
                          truelatL,truelatH,clon,dsinm,dsinm
 100  format('pdef ',i4,1x,i4,1x,'lccr',7(1x,f7.2),1x,2(f7.0,1x))
            write(31,110) nx+2,alonmin-rloninc,rloninc
 110  format('xdef ',i4,' linear ',f7.2,1x,f7.4)
            write(31,120) ny+2,alatmin-rlatinc,rlatinc
 120  format('ydef ',i4,' linear ',f7.2,1x,f7.4)
         elseif(iproj.eq.'POLSTR') then    !
         elseif(iproj.eq.'NORMER') then
            write(31,200)  jx_len,xlon(1,1),xlon(2,1)-xlon(1,1)
 200  format('xdef ',I3,' linear ',f9.4,' ',f9.4)
            write(31,210) iy-2
 210  format('ydef ',I3,' levels')
            write(31,220) (xlat(1,i),i=1,iy-2)
 220  format(10f7.2)
         elseif(iproj.eq.'ROTMER') then
         if(nfile.eq.1.or.month.eq.1) then
            write(*,*) 'Note that rotated Mercartor (ROTMER)' &
                   ,' projections are not supported by GrADS.'
            write(*,*) '  Although not exact, the eta.u projection' &
                   ,' in GrADS is somewhat similar.'
            write(*,*) ' FERRET, however, does support this projection.'
         endif
            write(31,230) jx_len,iy-2,plon,plat,dsinm/111000. &
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
         if(nfile.eq.1.and.idate1.eq.idate0) then
            nday=mod(idate1,10000)/100
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
        write(31,660) 'qv      ',kz,'specific moisture (kg/kg)   '
        write(31,600) 'ps      ',   'surface pressure (hPa)      '
        write(31,600) 'ts      ',   'surface air temperature     '
      endif
      write(31,700)
 700  format('endvars')
      close(31)
         ENDIF

         inquire(file=trim(Path_Input)//trim(DomainName)//'_'//filein  &
             ,exist=there)
         if(.not.there) then
            write(*,*) trim(Path_Input)//trim(DomainName)//'_'//filein &
                      ,' is not avaiable'
            stop
         endif
         open(10,file=trim(Path_Input)//trim(DomainName)//'_'//filein  &
                ,form='unformatted',recl=iy*jx*ibyte,access='direct')
         open(20,file=trim(Path_Input)//trim(DomainName)//'_'//fileout &
          ,form='unformatted',recl=(iy-2)*jx_len*ibyte,access='direct')
         mrec = 0
         do nday=1,nrecord
            do n=1,nvar
               do i=1,iy-2
               do j=1,jx_len
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
                  do j=1,jx_len
                     if(i_band.eq.1) then
                        c(j,i,n) = c(j,i,n)+b(j,i+1)
                     else
                        c(j,i,n) = c(j,i,n)+b(j+1,i+1)
                     endif
                  enddo
                  enddo
               enddo
            enddo
            do n=1,nvar
               do i=1,iy-2
               do j=1,jx_len
                  c(j,i,n) = c(j,i,n)/float(n_slice)
               enddo
               enddo
               mrec=mrec+1
               write(20,rec=mrec) ((c(j,i,n),j=1,jx_len),i=1,iy-2)
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
 199  format(I4)
      return
      end

      subroutine day2p(filename,iy,jx,kz,np,plev,i_band,ibyte &
                    ,Path_Input,DomainName,idate0,idate1,idate2,igrads)
      implicit none
      character*6 filename
      character*128 Path_Input
      character*20 DomainName
      integer iy,jx,kz,np,i_band,ibyte,idate0,idate1,idate2,igrads
      real*4  plev(np)
      integer iiy,jjx,kkz
      real*4  truelatL,truelatH
      real*4  dsinm,ptop,clat,clon,plat,plon,GRDFAC
      integer jgrads,ibigend
      character*6 iproj
      real*4, allocatable ::  sigma(:)
      character*4 :: chy
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
      character*16 filein
      character*14 fileout
      integer ntype,nfile,nyear,month,n_slice,mrec,nrec
      logical there
      real*4  alatmin,alatmax,alonmin,alonmax,rlatinc,rloninc
      real*4  centerj,centeri
      integer ny,nx
      integer i,j,k,n,jx_len
      integer nvar,n_month,nrecord,nday,mday
      real*4, allocatable ::  c(:,:,:),b(:,:),xlat(:,:),xlon(:,:)

      if(i_band.eq.1) then
         jx_len = jx
      else
         jx_len = jx-2
      endif

      allocate(sigma(kz+1))
      allocate(b(jx,iy))

      allocate(xlat(jx_len,iy-2))
      allocate(xlon(jx_len,iy-2))

      inquire(file=trim(Path_Input)//trim(DomainName)//'.INFO' &
             ,exist=there)
      if(.not.there) then
         write(*,*) trim(Path_Input)//trim(DomainName)//'.INFO' &
                   ,' is not avaiable'
         stop
      endif
      open(10,file=trim(Path_Input)//trim(DomainName)//'.INFO'  &
             ,form='unformatted',recl=jx*iy*ibyte,access='direct')
      read(10,rec=1) iiy,jjx,kkz,dsinm,clat,clon,plat,plon,GRDFAC  &
                    ,iproj,(sigma(k),k=1,kz+1),ptop,jgrads,ibigend &
                    ,truelatL,truelatH
      if(iiy.ne.iy.or.jjx.ne.jx.or.kkz.ne.kz) then
         write(*,*) 'iy,jx,kz in parameter = ',iy,jx,kz
         write(*,*) 'iy,jx,kz in DOMAIN.INFO ',iiy,jjx,kkz
         write(*,*) 'They are not consistent'
         stop
      endif
      read(10,rec=5) b
      do i=1,iy-2
      do j=1,jx_len
         if(i_band.eq.1) then
            xlat(j,i) = b(j,i+1)
         else
            xlat(j,i) = b(j+1,i+1)
         endif
      enddo
      enddo
      read(10,rec=6) b
      do i=1,iy-2
      do j=1,jx_len
         if(i_band.eq.1) then
            xlon(j,i) = b(j,i+1)
         else
            xlon(j,i) = b(j+1,i+1)
         endif
      enddo
      enddo
      close(10)
      deallocate(b)
      allocate(b(jx_len,iy-2))

      if(filename.eq.'ICBC_P') then
         nvar = np*6+2
      else
         write(*,*) 'filename is not correct'
         stop
      endif
      allocate(c(jx_len,iy-2,nvar))

      if(idate1.ge.1000010100) then        ! original file
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
         nrec = nvar
         n_slice = 24/6

         write(chy,199) nyear
         if(filename.eq.'ICBC_P') then
            filein = 'ICBC_P'//chy//chm(month)//'0100'
            fileout= 'ICBC_P'//chy//chm(month)//'01'
         endif

         IF(igrads.eq.1) THEN
         inquire(file= &
        trim(Path_Input)//trim(DomainName)//'_'//fileout//'.ctl'  &
                ,exist=there)
         if(there) then
            open(31,file= &
         trim(Path_Input)//trim(DomainName)//'_'//fileout//'.ctl' &
                                                  ,status='replace')
         else
            open(31,file= &
         trim(Path_Input)//trim(DomainName)//'_'//fileout//'.ctl' &
                                                      ,status='new')
         endif
         write(31,10) '^'//trim(DomainName)//'_'//fileout
  10     format('dset ',A34)
         write(31,20)
  20     format('title RegCM daily pressure level ICBC variables')
         write(31,30)
  30     format('options big_endian')
         write(31,50)
  50     format('undef -1.e34')
         if(iproj.eq.'LAMCON'.or.iproj.eq.'ROTMER') then
            alatmin= 999999.
            alatmax=-999999.
            do j=1,jx_len
               if(xlat(j,1).lt.alatmin) alatmin=xlat(j,1)
               if(xlat(j,iy-2).gt.alatmax) alatmax=xlat(j,iy-2)
            enddo
            alonmin= 999999.
            alonmax=-999999.
            do i=1,iy-2
            do j=1,jx_len
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
                     if(abs(clon-xlon(j,i)).lt. &
                        abs(clon-(xlon(j,i)-360.))) then
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
            centerj=jx_len/2.
            centeri=(iy-2)/2.
         endif

         if(iproj.eq.'LAMCON') then        ! Lambert projection
            write(31,100) jx_len,iy-2,clat,clon,centerj,centeri, &
                          truelatL,truelatH,clon,dsinm,dsinm
 100  format('pdef ',i4,1x,i4,1x,'lccr',7(1x,f7.2),1x,2(f7.0,1x))
            write(31,110) nx+2,alonmin-rloninc,rloninc
 110  format('xdef ',i4,' linear ',f7.2,1x,f7.4)
            write(31,120) ny+2,alatmin-rlatinc,rlatinc
 120  format('ydef ',i4,' linear ',f7.2,1x,f7.4)
         elseif(iproj.eq.'POLSTR') then    !
         elseif(iproj.eq.'NORMER') then
            write(31,200)  jx_len,xlon(1,1),xlon(2,1)-xlon(1,1)
 200  format('xdef ',I3,' linear ',f9.4,' ',f9.4)
            write(31,210) iy-2
 210  format('ydef ',I3,' levels')
            write(31,220) (xlat(1,i),i=1,iy-2)
 220  format(10f7.2)
         elseif(iproj.eq.'ROTMER') then
         if(nfile.eq.1.or.month.eq.1) then
            write(*,*) 'Note that rotated Mercartor (ROTMER)' &
                   ,' projections are not supported by GrADS.'
            write(*,*) '  Although not exact, the eta.u projection' &
                   ,' in GrADS is somewhat similar.'
            write(*,*) ' FERRET, however, does support this projection.'
         endif
            write(31,230) jx_len,iy-2,plon,plat,dsinm/111000. &
                                          ,dsinm/111000.*.95238
 230  format('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
            write(31,110) nx+2,alonmin-rloninc,rloninc
            write(31,120) ny+2,alatmin-rlatinc,rlatinc
         else
            write(*,*) 'Are you sure your map projection is right ?'
            stop
         endif
      write(31,318) np,(plev(k),k=1,np)
 318  format('zdef ',I2,' levels ',30f7.2)
         if(nfile.eq.1.and.idate1.eq.idate0) then
            nday=mod(idate1,10000)/100
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
      write(31,500) 6+2
 500  format('vars ',I3)
 600  format(A8,'0 99 ',A36)
 611  format(A8,'0 33,105 ',A36)
 612  format(A8,'0 34,105 ',A36)
      if(filename.eq.'ICBC_P') then
 660  format(A8,I2,' 0 ',A36)
 661  format(A8,I2,' 33,100 ',A36)
 662  format(A8,I2,' 34,100 ',A36)
        if(iproj.eq.'LAMCON') then
           write(31,661) 'u       ',np,'westerly wind (m/s)         '
           write(31,662) 'v       ',np,'southerly wind (m/s)        '
        else
           write(31,660) 'u       ',np,'westerly wind (m/s)         '
           write(31,660) 'v       ',np,'southerly wind (m/s)        '
        endif
        write(31,660) 'h       ',np,'geopotential height (m)     '
        write(31,660) 't       ',np,'air temperature (degree)    '
        write(31,660) 'rh      ',np,'relative moisture (%)       '
        write(31,660) 'qv      ',np,'specific moisture (kg/kg)   '
        write(31,600) 'ps      ',   'surface pressure (hPa)      '
        write(31,600) 'slp     ',   'sea level pressure (hPa)    '
      endif
      write(31,700)
 700  format('endvars')
      close(31)
         ENDIF

         inquire(file=trim(Path_Input)//trim(DomainName)//'_'//filein  &
             ,exist=there)
         if(.not.there) then
            write(*,*) trim(Path_Input)//trim(DomainName)//'_'//filein &
                      ,' is not avaiable'
            stop
         endif
         open(10,file=trim(Path_Input)//trim(DomainName)//'_'//filein  &
          ,form='unformatted',recl=(iy-2)*jx_len*ibyte,access='direct')
         open(20,file=trim(Path_Input)//trim(DomainName)//'_'//fileout &
          ,form='unformatted',recl=(iy-2)*jx_len*ibyte,access='direct')
         mrec = 0
         do nday=1,nrecord
            do n=1,nvar
               do i=1,iy-2
               do j=1,jx_len
                  c(j,i,n) = 0.0
               enddo
               enddo
            enddo
            do mday = 1,n_slice
               do n=1,nvar
                  nrec=nrec+1
                  read(10,rec=nrec) b
                  do i=1,iy-2
                  do j=1,jx_len
                     c(j,i,n) = c(j,i,n)+b(j,i)
                  enddo
                  enddo
               enddo
            enddo
            do n=1,nvar
               do i=1,iy-2
               do j=1,jx_len
                  c(j,i,n) = c(j,i,n)/float(n_slice)
               enddo
               enddo
               mrec=mrec+1
               write(20,rec=mrec) ((c(j,i,n),j=1,jx_len),i=1,iy-2)
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
 199  format(I4)
      return
      end

      subroutine day3(filename,iy,jx,kz,nsg,i_band,ibyte, Path_Input &
                    ,DomainName,Path_Output,idate0,idate1,idate2,igrads)
      implicit none
      character*3 filename
      character*128 Path_Input,Path_Output
      character*20 DomainName
      integer iy,jx,kz,nsg,i_band,ibyte,idate0,idate1,idate2,igrads
      integer iiy,jjx,kkz
      integer mdate0,ibltyp,icup,ipptls,iboudy
      real*4  truelatL,truelatH
      real*4  dxsp,ptsp,clat,clon,plat,plon
      real*4  dto,dtb,dtr,dtc,dt
      integer iotyp
      character*6 iproj
      real*4, allocatable ::  sigma(:)
      character*4 :: chy
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
      character*3 chnsg
      character*14 filein
      character*12 fileout
      integer ntype,nfile,nyear,month,n_slice,mrec,nrec
      logical there
      real*4  alatmin,alatmax,alonmin,alonmax,rlatinc,rloninc
      real*4  centerj,centeri
      integer ny,nx
      integer i,j,k,n,jx_len
      integer nvar,n_month,nrecord,nday,mday
      integer idatex
      real*4, allocatable ::  c(:,:,:),b(:,:),xlat(:,:),xlon(:,:)
      real*4, allocatable ::  xlat_s(:,:),xlon_s(:,:),a(:,:)

      if(i_band.eq.1) then
         jx_len = jx
      else
         jx_len = jx-2
      endif

      allocate(sigma(kz+1))
      allocate(b(jx_len*nsg,(iy-2)*nsg))
      allocate(xlat(jx_len,iy-2))
      allocate(xlon(jx_len,iy-2))

      allocate(xlat_s(jx_len*nsg,(iy-2)*nsg))
      allocate(xlon_s(jx_len*nsg,(iy-2)*nsg))
      allocate(a(jx*nsg,iy*nsg))

      if(nsg.lt.10) then
         write(chnsg,71) nsg
      else if(nsg.lt.100) then
         write(chnsg,72) nsg
      else
         write(chnsg,73) nsg
      endif
 71   format('00',I1)
 72   format('0', I2)
 73   format(I3)
      inquire(file=trim(Path_Input)//trim(DomainName)//chnsg//'.INFO' &
                  ,exist=there)
      if(.not.there) then
         write(*,*) trim(Path_Input)//trim(DomainName)//chnsg//'.INFO' &
                   ,' is not avaiable'
         stop
      endif
      open(10,file=trim(Path_Input)//trim(DomainName)//chnsg//'.INFO' &
       ,form='unformatted',recl=jx*nsg*iy*nsg*ibyte,access='direct')
      read(10,rec=5) a
      do i=nsg+1,(iy-1)*nsg
      do j=nsg+1,(jx-1)*nsg
         xlat_s(j-nsg,i-nsg) = a(j,i)
      enddo
      enddo
      read(10,rec=6) a
      do i=nsg+1,(iy-1)*nsg
      do j=nsg+1,(jx-1)*nsg
         xlon_s(j-nsg,i-nsg) = a(j,i)
      enddo
      enddo
      close(10)
      deallocate(a)

      inquire(file=trim(Path_Output)//'OUT_HEAD',exist=there)
      if(.not.there) then
         write(*,*) trim(Path_Output)//'OUT_HEAD',' is not avaiable'
         stop
      endif
      open(10,file=trim(Path_Output)//'OUT_HEAD',form='unformatted' &
             ,recl=jx_len*(iy-2)*ibyte,access='direct')
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
      allocate(c(jx_len*nsg,(iy-2)*nsg,nvar))

      if(idate1.ge.1000010100) then        ! original file
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
      write(*,*) 'Calculate the daily mean of    ',filename,nyear,month
         else
      write(*,*) '                               ',filename,nyear,month
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

         write(chy,199) nyear
         if(filename.eq.'SUB') then
            filein = 'SUB.'//chy//chm(month)//'0100'
            fileout= 'SUB.'//chy//chm(month)//'01'
         endif

         IF(igrads.eq.1) THEN
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
            do j=1,jx_len
               if(xlat(j,1).lt.alatmin) alatmin=xlat(j,1)
               if(xlat(j,iy-2).gt.alatmax) alatmax=xlat(j,iy-2)
            enddo
            alonmin= 999999.
            alonmax=-999999.
            do i=1,iy-2
            do j=1,jx_len
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
                     if(abs(clon-xlon(j,i)).lt. &
                        abs(clon-(xlon(j,i)-360.))) then
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
            centerj=jx_len/2.
            centeri=(iy-2)/2.
         endif

         if(iproj.eq.'LAMCON') then        ! Lambert projection
            write(31,100) jx_len*nsg,(iy-2)*nsg,clat,clon, &
                          centerj*nsg,centeri*nsg, &
                 truelatL,truelatH,clon,dxsp/nsg*1000.,dxsp/nsg*1000.
 100  format('pdef ',i4,1x,i4,1x,'lccr',7(1x,f7.2),1x,2(f7.0,1x))
            write(31,110) nx+2,alonmin-rloninc,rloninc
 110  format('xdef ',i4,' linear ',f7.2,1x,f7.4)
            write(31,120) ny+2,alatmin-rlatinc,rlatinc
 120  format('ydef ',i4,' linear ',f7.2,1x,f7.4)
         elseif(iproj.eq.'POLSTR') then    !
         elseif(iproj.eq.'NORMER') then
            write(31,200) jx_len*nsg,xlon_s(1,1),xlon_s(2,1)-xlon_s(1,1)
 200  format('xdef ',I3,' linear ',f9.4,' ',f9.4)
            write(31,210) (iy-2)*nsg
 210  format('ydef ',I3,' levels')
            write(31,220) (xlat_s(1,i),i=1,(iy-2)*nsg)
 220  format(10f7.2)
         elseif(iproj.eq.'ROTMER') then
         if(nfile.eq.1.or.month.eq.1) then
            write(*,*) 'Note that rotated Mercartor (ROTMER)' &
                   ,' projections are not supported by GrADS.'
            write(*,*) '  Although not exact, the eta.u projection' &
                   ,' in GrADS is somewhat similar.'
            write(*,*) ' FERRET, however, does support this projection.'
         endif
            write(31,230) jx_len*nsg,(iy-2)*nsg,plon,plat, &
                          dxsp/111./nsg ,dxsp/111.*.95238/nsg
 230  format('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
            write(31,110) nx+2,alonmin-rloninc,rloninc
            write(31,120) ny+2,alatmin-rlatinc,rlatinc
         else
            write(*,*) 'Are you sure your map projection is right ?'
            stop
         endif
      write(31,300) (1013.25-ptsp*10.)*(sigma(1)+sigma(2))*.5+ptsp*10.
 300  format('zdef 1',' levels ',f7.2)
         if(nfile.eq.1.and.idate1.eq.idate0) then
            nday=mod(idate1,10000)/100
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
        write(31,600) 'scv     ','snow amount (mm, water equivalent)   '
        write(31,600) 'sena    ','sensible heat flux (W/m2)            '
        write(31,600) 'prcv    ','convective precipitation (mm/day)    '
        write(31,600) 'ps      ','surface pressure (hPa)               '
         write(31,700)
 700  format('endvars')
         close(31)
         ENDIF

         inquire(file=trim(Path_Output)//filein,exist=there)
         if(.not.there) then
            write(*,*) trim(Path_Output)//filein,' is not avaiable'
            stop
         endif
         if(iotyp.eq.1) then
            open(10,file=trim(Path_Output)//filein,form='unformatted' &
                   ,recl=(iy-2)*nsg*jx_len*nsg*ibyte,access='direct')
         else if(iotyp.eq.2) then
            open(10,file=trim(Path_Output)//filein,form='unformatted')
         endif
         open(20,file=trim(Path_Output)//fileout,form='unformatted' &
                ,recl=(iy-2)*nsg*jx_len*nsg*ibyte,access='direct')
         mrec = 0
         do nday=1,nrecord
            do n=1,nvar
               do i=1,(iy-2)*nsg
               do j=1,jx_len*nsg
                  c(j,i,n) = 0.0
               enddo
               enddo
            enddo
            if(iotyp.eq.2) read(10) idatex
            do mday = 1,n_slice
               do n=1,nvar
                  if(iotyp.eq.1) then
                     nrec=nrec+1
                     read(10,rec=nrec) b
                  else if(iotyp.eq.2) then
                     read(10) b
                  endif
                  do i=1,(iy-2)*nsg
                  do j=1,jx_len*nsg
                     c(j,i,n) = c(j,i,n)+b(j,i)
                  enddo
                  enddo
               enddo
            enddo
            do n=1,nvar
               do i=1,(iy-2)*nsg
               do j=1,jx_len*nsg
                  c(j,i,n) = c(j,i,n)/float(n_slice)
               enddo
               enddo
               mrec=mrec+1
               write(20,rec=mrec) ((c(j,i,n),j=1,jx_len*nsg)  &
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

      deallocate(xlat_s)
      deallocate(xlon_s)

 199  format(I4)
      return
      end

      subroutine mon(filename,iy,jx,kz,ntr,i_band,ibyte, Path_Output &
                    ,idate1,idate2,igrads)
      implicit none
      character*3 filename
      character*128 Path_Output
      integer iy,jx,kz,ntr,i_band,ibyte,idate1,idate2,igrads
      integer iiy,jjx,kkz
      integer mdate0,ibltyp,icup,ipptls,iboudy
      real*4  truelatL,truelatH
      real*4  dxsp,ptsp,clat,clon,plat,plon
      real*4  dto,dtb,dtr,dtc
      integer iotyp
      character*6 iproj
      real*4, allocatable ::  sigma(:)
      character*4 :: chy
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
      integer i,j,k,n,itr,jx_len
      integer nvar,n_month,nrecord,nday
      logical there
      real*4  alatmin,alatmax,alonmin,alonmax,rlatinc,rloninc
      real*4  centerj,centeri
      integer ny,nx

      real*4, allocatable ::  c(:,:,:),b(:,:),xlat(:,:),xlon(:,:)

      if(i_band.eq.1) then
         jx_len = jx
      else
         jx_len = jx-2
      endif

      allocate(sigma(kz+1))
      allocate(b(jx_len,iy-2))
      allocate(xlat(jx_len,iy-2))
      allocate(xlon(jx_len,iy-2))

      inquire(file=trim(Path_Output)//'OUT_HEAD',exist=there)
      if(.not.there) then
         write(*,*) trim(Path_Output)//'OUT_HEAD',' is not avaiable'
         stop
      endif
      open(10,file=trim(Path_Output)//'OUT_HEAD',form='unformatted' &
             ,recl=jx_len*(iy-2)*ibyte,access='direct')
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
         nvar = kz*4+10
      else if(filename.eq.'CHE') then
         nvar = ntr*kz+kz*3+ntr*7+3
!        nvar = ntr*kz+kz*3+ntr*7+2      ! for RegCM3, one record less
      else
         write(*,*) 'filename is not correct'
         stop
      endif
      allocate(c(jx_len,iy-2,nvar))

      if(idate1.ge.1000010100) then        ! original file
         n_month = (idate2/1000000-idate1/1000000)*12 &
                 + (mod(idate2/10000,100)-mod(idate1/10000,100))
         if(mod(idate2,10000).gt.0100) n_month = n_month+1
         ntype = 0
         month = mod(idate1/10000,100)
         nyear = idate1/1000000
      else if(idate1.ge.10000101) then     ! daily mean file
         n_month = (idate2/10000-idate1/10000)*12 &
                 + (mod(idate2/100,100)-mod(idate1/100,100))
         if(mod(idate2,100).gt.01) n_month = n_month+1
         ntype = 1
         month = mod(idate1/100,100)
         nyear = idate1/10000
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
             ,recl=(iy-2)*jx_len*ibyte,access='direct')
      mrec = 0

      IF(igrads.eq.1) THEN
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
         do j=1,jx_len
            if(xlat(j,1).lt.alatmin) alatmin=xlat(j,1)
            if(xlat(j,iy-2).gt.alatmax) alatmax=xlat(j,iy-2)
         enddo
         alonmin= 999999.
         alonmax=-999999.
         do i=1,iy-2
         do j=1,jx_len
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
                  if(abs(clon-xlon(j,i)).lt. &
                     abs(clon-(xlon(j,i)-360.))) then
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
         centerj=jx_len/2.
         centeri=(iy-2)/2.
      endif

      if(iproj.eq.'LAMCON') then        ! Lambert projection
         write(31,100) jx_len,iy-2,clat,clon,centerj,centeri, &
                       truelatL,truelatH,clon,dxsp*1000.,dxsp*1000.
 100  format('pdef ',i4,1x,i4,1x,'lccr',7(1x,f7.2),1x,2(f7.0,1x))
         write(31,110) nx+2,alonmin-rloninc,rloninc
 110  format('xdef ',i4,' linear ',f7.2,1x,f7.4)
         write(31,120) ny+2,alatmin-rlatinc,rlatinc
 120  format('ydef ',i4,' linear ',f7.2,1x,f7.4)
      elseif(iproj.eq.'POLSTR') then    !
      elseif(iproj.eq.'NORMER') then
         write(31,200)  jx_len,xlon(1,1),xlon(2,1)-xlon(1,1)
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
         write(31,230) jx_len,iy-2,plon,plat,dxsp/111. &
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
      if(ntype.eq.0) then
         nyear=idate1/1000000
         month=(idate1-nyear*1000000)/10000
      else if(ntype.eq.1) then
         nyear=idate1/10000
         month=(idate1-nyear*10000)/100
      endif
      write(31,400)n_month,cday(16),chmc(month),nyear
 400  format('tdef ',I4,' linear ',A2,A3,I4,' ','1mo')
      if(filename.eq.'ATM') then
         write(31,500) 6+5
      else if(filename.eq.'RAD') then
         write(31,500) 4+10
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
        write(31,600) 'scv     ','snow amount (mm, water equivalent)   '
        write(31,600) 'sena    ','sensible heat flux (W/m2)            '
        write(31,600) 'flw     ','net infrared energy flux (W/m2)      '
        write(31,600) 'fsw     ','net absorbed solar energy flux (W/m2)'
        write(31,600) 'flwd    ','downward infrared energy flux (W/m2) '
        write(31,600) 'sina    ','incident solar energy flux (W/m2)    '
        write(31,600) 'prcv    ','convective precipitation (mm/day)    '
        write(31,600) 'ps      ','surface pressure (hPa)               '
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
        write(31,600) 'ps      ',   'surface pressure (hPa)      '
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
        write(31,600)'ps    ',   'surface pressure (hPa)               '
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
       write(31,600) 'ps      ','surface pressure (hPa)               '
      endif
      write(31,700)
 700  format('endvars')
      close(31)
      ENDIF

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
      write(*,*)'Calculate the monthly mean of    ',filename,nyear,month
         else
      write(*,*)'                                 ',filename,nyear,month
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

         write(chy,199) nyear
         if(filename.eq.'ATM') then
            filein = 'ATM.'//chy//chm(month)//'01'
         else if(filename.eq.'SRF') then
            filein = 'SRF.'//chy//chm(month)//'01'
         else if(filename.eq.'RAD') then
            filein = 'RAD.'//chy//chm(month)//'01'
         else if(filename.eq.'CHE') then
            filein = 'CHE.'//chy//chm(month)//'01'
         endif
         inquire(file=trim(Path_Output)//filein,exist=there)
         if(.not.there) then
            write(*,*) trim(Path_Output)//filein,' is not avaiable'
            write(*,*) 'Note: Daily mean files are required.'
            stop
         endif
         open(10,file=trim(Path_Output)//filein,form='unformatted' &
                ,recl=(iy-2)*jx_len*ibyte,access='direct')
         nrec = 0

         do n=1,nvar
            do i=1,iy-2
            do j=1,jx_len
               c(j,i,n) = 0.0
            enddo
            enddo
         enddo

         do nday=1,nrecord
            do n=1,nvar
               nrec=nrec+1
               read(10,rec=nrec) b
               do i=1,iy-2
               do j=1,jx_len
                  c(j,i,n) = c(j,i,n)+b(j,i)
               enddo
               enddo
            enddo
         enddo
         do n=1,nvar
            do i=1,iy-2
            do j=1,jx_len
               c(j,i,n) = c(j,i,n)/float(nrecord)
            enddo
            enddo
            mrec=mrec+1
            write(20,rec=mrec) ((c(j,i,n),j=1,jx_len),i=1,iy-2)
         enddo
         close(10)
      enddo
      close(20)
      deallocate(sigma)
      deallocate(b)
      deallocate(xlat)
      deallocate(xlon)
      deallocate(c)
 199  format(I4)
      return
      end

      subroutine monp(filename,iy,jx,kz,np,plev,i_band,ibyte, Path_Output &
                    ,idate1,idate2,igrads)
      implicit none
      character*5 filename
      character*128 Path_Output
      integer iy,jx,kz,np,i_band,ibyte,idate1,idate2,igrads
      real*4  plev(np)
      integer iiy,jjx,kkz
      integer mdate0,ibltyp,icup,ipptls,iboudy
      real*4  truelatL,truelatH
      real*4  dxsp,ptsp,clat,clon,plat,plon
      real*4  dto,dtb,dtr,dtc
      integer iotyp
      character*6 iproj
      real*4, allocatable ::  sigma(:)
      character*4 :: chy
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
      integer i,j,k,n,jx_len
      integer nvar,n_month,nrecord,nday
      logical there
      real*4  alatmin,alatmax,alonmin,alonmax,rlatinc,rloninc
      real*4  centerj,centeri
      integer ny,nx

      real*4, allocatable ::  c(:,:,:),b(:,:),xlat(:,:),xlon(:,:)

      if(i_band.eq.1) then
         jx_len = jx
      else
         jx_len = jx-2
      endif

      allocate(sigma(kz+1))
      allocate(b(jx_len,iy-2))
      allocate(xlat(jx_len,iy-2))
      allocate(xlon(jx_len,iy-2))

      inquire(file=trim(Path_Output)//'OUT_HEAD',exist=there)
      if(.not.there) then
         write(*,*) trim(Path_Output)//'OUT_HEAD',' is not avaiable'
         stop
      endif
      open(10,file=trim(Path_Output)//'OUT_HEAD',form='unformatted' &
             ,recl=jx_len*(iy-2)*ibyte,access='direct')
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

      if(filename.eq.'ATM_P') then
         nvar = np*8+6
      else
         write(*,*) 'filename is not correct'
         stop
      endif
      allocate(c(jx_len,iy-2,nvar))

      if(idate1.ge.1000010100) then        ! original file
         n_month = (idate2/1000000-idate1/1000000)*12 &
                 + (mod(idate2/10000,100)-mod(idate1/10000,100))
         if(mod(idate2,10000).gt.0100) n_month = n_month+1
         ntype = 0
         month = mod(idate1/10000,100)
         nyear = idate1/1000000
      else if(idate1.ge.10000101) then     ! daily mean file
         n_month = (idate2/10000-idate1/10000)*12 &
                 + (mod(idate2/100,100)-mod(idate1/100,100))
         if(mod(idate2,100).gt.01) n_month = n_month+1
         ntype = 1
         month = mod(idate1/100,100)
         nyear = idate1/10000
      else
         write(*,*) 'date should be in 10 digits, or 8 digits'
         stop
      endif

      fileout= 'ATM_P.mon'

      open(20,file=trim(Path_Output)//fileout,form='unformatted' &
             ,recl=(iy-2)*jx_len*ibyte,access='direct')
      mrec = 0

      IF(igrads.eq.1) THEN
      inquire(file=trim(Path_Output)//fileout//'.ctl',exist=there)
      if(there) then
       open(31,file=trim(Path_Output)//fileout//'.ctl',status='replace')
      else
         open(31,file=trim(Path_Output)//fileout//'.ctl',status='new')
      endif
      write(31,10) fileout
  10  format('dset ^',A9)
      write(31,20)
  20  format('title RegCM monthly pressure level output variables')
      write(31,30)
  30  format('options big_endian')
      write(31,50)
  50  format('undef -1.e34')
      if(iproj.eq.'LAMCON'.or.iproj.eq.'ROTMER') then
         alatmin= 999999.
         alatmax=-999999.
         do j=1,jx_len
            if(xlat(j,1).lt.alatmin) alatmin=xlat(j,1)
            if(xlat(j,iy-2).gt.alatmax) alatmax=xlat(j,iy-2)
         enddo
         alonmin= 999999.
         alonmax=-999999.
         do i=1,iy-2
         do j=1,jx_len
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
                  if(abs(clon-xlon(j,i)).lt. &
                     abs(clon-(xlon(j,i)-360.))) then
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
         centerj=jx_len/2.
         centeri=(iy-2)/2.
      endif

      if(iproj.eq.'LAMCON') then        ! Lambert projection
         write(31,100) jx_len,iy-2,clat,clon,centerj,centeri, &
                       truelatL,truelatH,clon,dxsp*1000.,dxsp*1000.
 100  format('pdef ',i4,1x,i4,1x,'lccr',7(1x,f7.2),1x,2(f7.0,1x))
         write(31,110) nx+2,alonmin-rloninc,rloninc
 110  format('xdef ',i4,' linear ',f7.2,1x,f7.4)
         write(31,120) ny+2,alatmin-rlatinc,rlatinc
 120  format('ydef ',i4,' linear ',f7.2,1x,f7.4)
      elseif(iproj.eq.'POLSTR') then    !
      elseif(iproj.eq.'NORMER') then
         write(31,200)  jx_len,xlon(1,1),xlon(2,1)-xlon(1,1)
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
         write(31,230) jx_len,iy-2,plon,plat,dxsp/111. &
                                          ,dxsp/111.*.95238
 230  format('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
         write(31,110) nx+2,alonmin-rloninc,rloninc
         write(31,120) ny+2,alatmin-rlatinc,rlatinc
      else
         write(*,*) 'Are you sure your map projection is right ?'
         stop
      endif
      write(31,318) np,(plev(k),k=1,np)
 318  format('zdef ',I2,' levels ',30f7.2)
      if(ntype.eq.0) then
         nyear=idate1/1000000
         month=(idate1-nyear*1000000)/10000
      else if(ntype.eq.1) then
         nyear=idate1/10000
         month=(idate1-nyear*10000)/100
      endif
      write(31,400)n_month,cday(16),chmc(month),nyear
 400  format('tdef ',I4,' linear ',A2,A3,I4,' ','1mo')
      write(31,500) 8+6
 500  format('vars ',I3)
 600  format(A8,'0 99 ',A36)
 660  format(A8,I2,' 0 ',A36)
 661  format(A8,I2,' 33,100 ',A36)
 662  format(A8,I2,' 34,100 ',A36)
      if(iproj.eq.'LAMCON') then
         write(31,661) 'u       ',np,'westerly wind (m/s)         '
         write(31,662) 'v       ',np,'southerly wind (m/s)        '
      else
         write(31,660) 'u       ',np,'westerly wind (m/s)         '
         write(31,660) 'v       ',np,'southerly wind (m/s)        '
      endif
      write(31,660) 'omega   ',np,'omega (hPa/s)   p-velocity  '
      write(31,660) 'h       ',np,'geopotential height (m)     '
      write(31,660) 't       ',np,'air temperature (deg, K)    '
      write(31,660) 'rh      ',np,'relative moisture (%)       '
      write(31,660) 'qv      ',np,'water vapor mixing ratio    '
      write(31,660) 'qc      ',np,'cloud water mixing ratio    '
      write(31,600) 'ps      ',   'surface pressure (hPa)      '
      write(31,600) 'slp     ',   'sea level pressure (hPa)    '
      write(31,600) 'tpr     ',   'total precipitation(mm/day) '
      write(31,600) 'tgb     ',   'lower groud temp. in BATS   '
      write(31,600) 'swt     ',   'total soil water in mm H2O  '
      write(31,600) 'rno     ',   'accumulated infiltration    '
      write(31,700)
 700  format('endvars')
      close(31)
      ENDIF

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
       write(*,*)'Calculate the monthly mean of  ',filename,nyear,month
         else
       write(*,*)'                               ',filename,nyear,month
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

         write(chy,199) nyear
         filein = 'ATM_P.'//chy//chm(month)//'01'
         inquire(file=trim(Path_Output)//filein,exist=there)
         if(.not.there) then
            write(*,*) trim(Path_Output)//filein,' is not avaiable'
            write(*,*) 'Note: Daily mean files are required.'
            stop
         endif
         open(10,file=trim(Path_Output)//filein,form='unformatted' &
                ,recl=(iy-2)*jx_len*ibyte,access='direct')
         nrec = 0

         do n=1,nvar
            do i=1,iy-2
            do j=1,jx_len
               c(j,i,n) = 0.0
            enddo
            enddo
         enddo

         do nday=1,nrecord
            do n=1,nvar
               nrec=nrec+1
               read(10,rec=nrec) b
               do i=1,iy-2
               do j=1,jx_len
                  c(j,i,n) = c(j,i,n)+b(j,i)
               enddo
               enddo
            enddo
         enddo
         do n=1,nvar
            do i=1,iy-2
            do j=1,jx_len
               c(j,i,n) = c(j,i,n)/float(nrecord)
            enddo
            enddo
            mrec=mrec+1
            write(20,rec=mrec) ((c(j,i,n),j=1,jx_len),i=1,iy-2)
         enddo
         close(10)
      enddo
      close(20)
      deallocate(sigma)
      deallocate(b)
      deallocate(xlat)
      deallocate(xlon)
      deallocate(c)
 199  format(I4)
      return
      end

      subroutine mon2(filename,iy,jx,kz,i_band,ibyte, Path_Input,DomainName &
                     ,idate1,idate2,igrads)
      implicit none
      character*4 filename
      character*128 Path_Input
      character*20 DomainName
      integer iy,jx,kz,i_band,ibyte,idate1,idate2,igrads
      integer iiy,jjx,kkz
      real*4  truelatL,truelatH
      real*4  dsinm,ptop,clat,clon,plat,plon,GRDFAC
      integer jgrads,ibigend
      character*6 iproj
      real*4, allocatable ::  sigma(:)
      character*4 :: chy
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
      character*8 fileout
      integer ntype,nfile,nyear,month,mrec,nrec
      logical there
      real*4  alatmin,alatmax,alonmin,alonmax,rlatinc,rloninc
      real*4  centerj,centeri
      integer ny,nx
      integer i,j,k,n,jx_len
      integer nvar,n_month,nrecord,nday
      real*4, allocatable ::  c(:,:,:),b(:,:),xlat(:,:),xlon(:,:)

      if(i_band.eq.1) then
         jx_len = jx
      else
         jx_len = jx-2
      endif

      allocate(sigma(kz+1))
      allocate(b(jx,iy))
      allocate(xlat(jx_len,iy-2))
      allocate(xlon(jx_len,iy-2))


      inquire(file=trim(Path_Input)//trim(DomainName)//'.INFO'  &
             ,exist=there)
      if(.not.there) then
         write(*,*) trim(Path_Input)//trim(DomainName)//'.INFO'  &
                   ,' is not avaiable'
         stop
      endif
      open(10,file=trim(Path_Input)//trim(DomainName)//'.INFO'  &
             ,form='unformatted',recl=jx*iy*ibyte,access='direct')
      read(10,rec=1) iiy,jjx,kkz,dsinm,clat,clon,plat,plon,GRDFAC  &
                    ,iproj,(sigma(k),k=1,kz+1),ptop,jgrads,ibigend &
                    ,truelatL,truelatH
      if(iiy.ne.iy.or.jjx.ne.jx.or.kkz.ne.kz) then
         write(*,*) 'iy,jx,kz in parameter = ',iy,jx,kz
         write(*,*) 'iy,jx,kz in DOMAIN.INFO ',iiy,jjx,kkz
         write(*,*) 'They are not consistent'
         stop
      endif
      read(10,rec=5) b
      do i=1,iy-2
      do j=1,jx_len
         if(i_band.eq.1) then
            xlat(j,i) = b(j,i+1)
         else
            xlat(j,i) = b(j+1,i+1)
         endif
      enddo
      enddo
      read(10,rec=6) b
      do i=1,iy-2
      do j=1,jx_len
         if(i_band.eq.1) then
            xlon(j,i) = b(j,i+1)
         else
            xlon(j,i) = b(j+1,i+1)
         endif
      enddo
      enddo
      close(10)

      deallocate(b)
      allocate(b(jx_len,iy-2))

      if(filename.eq.'ICBC') then
         nvar = kz*4+2
      else
         write(*,*) 'filename is not correct'
         stop
      endif
      allocate(c(jx_len,iy-2,nvar))

      if(idate1.ge.1000010100) then        ! original file
         n_month = (idate2/1000000-idate1/1000000)*12 &
                 + (mod(idate2/10000,100)-mod(idate1/10000,100))
         if(mod(idate2,10000).gt.0100) n_month = n_month+1
         ntype = 0
         month = mod(idate1/10000,100)
         nyear = idate1/1000000
      else if(idate1.ge.10000101) then     ! daily mean file
         n_month = (idate2/10000-idate1/10000)*12 &
                 + (mod(idate2/100,100)-mod(idate1/100,100))
         if(mod(idate2,100).gt.01) n_month = n_month+1
         ntype = 1
         month = mod(idate1/100,100)
         nyear = idate1/10000
      else
         write(*,*) 'date should be in 10 digits, or 8 digits'
         stop
      endif

      if(filename.eq.'ICBC') then
         fileout= 'ICBC.mon'
      endif

      open(20,file=trim(Path_Input)//trim(DomainName)//'_'//fileout  &
          ,form='unformatted',recl=(iy-2)*jx_len*ibyte,access='direct')
      mrec = 0

      IF(igrads.eq.1) THEN
      inquire(file= &
            trim(Path_Input)//trim(DomainName)//'_'//fileout//'.ctl' &
             ,exist=there)
      if(there) then
         open(31,file= &
             trim(Path_Input)//trim(DomainName)//'_'//fileout//'.ctl' &
                                                   ,status='replace')
      else
         open(31,file= &
             trim(Path_Input)//trim(DomainName)//'_'//fileout//'.ctl' &
                                                       ,status='new')
      endif
      write(31,10) '^'//trim(DomainName)//'_'//fileout
  10  format('dset ',A30)
      write(31,20)
  20  format('title RegCM monthly mean model level ICBC variables')
      write(31,30)
  30  format('options big_endian')
      write(31,50)
  50  format('undef -1.e34')
      if(iproj.eq.'LAMCON'.or.iproj.eq.'ROTMER') then
         alatmin= 999999.
         alatmax=-999999.
         do j=1,jx_len
            if(xlat(j,1).lt.alatmin) alatmin=xlat(j,1)
            if(xlat(j,iy-2).gt.alatmax) alatmax=xlat(j,iy-2)
         enddo
         alonmin= 999999.
         alonmax=-999999.
         do i=1,iy-2
         do j=1,jx_len
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
                  if(abs(clon-xlon(j,i)).lt. &
                     abs(clon-(xlon(j,i)-360.))) then
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
         centerj=jx_len/2.
         centeri=(iy-2)/2.
      endif

      if(iproj.eq.'LAMCON') then        ! Lambert projection
         write(31,100) jx_len,iy-2,clat,clon,centerj,centeri, &
                       truelatL,truelatH,clon,dsinm,dsinm
 100  format('pdef ',i4,1x,i4,1x,'lccr',7(1x,f7.2),1x,2(f7.0,1x))
         write(31,110) nx+2,alonmin-rloninc,rloninc
 110  format('xdef ',i4,' linear ',f7.2,1x,f7.4)
         write(31,120) ny+2,alatmin-rlatinc,rlatinc
 120  format('ydef ',i4,' linear ',f7.2,1x,f7.4)
      elseif(iproj.eq.'POLSTR') then    !
      elseif(iproj.eq.'NORMER') then
         write(31,200)  jx_len,xlon(1,1),xlon(2,1)-xlon(1,1)
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
         write(31,230) jx_len,iy-2,plon,plat,dsinm/111000. &
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
      if(ntype.eq.0) then
         nyear=idate1/1000000
         month=(idate1-nyear*1000000)/10000
      else if(ntype.eq.1) then
         nyear=idate1/10000
         month=(idate1-nyear*10000)/100
      endif
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
        write(31,660) 'qv      ',kz,'specific moisture (kg/kg)   '
        write(31,600) 'ps      ',   'surface pressure (hPa)      '
        write(31,600) 'ts      ',   'surface air temperature     '
      endif
      write(31,700)
 700  format('endvars')
      close(31)
      ENDIF

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

         write(chy,199) nyear
         if(filename.eq.'ICBC') then
            filein = 'ICBC'//chy//chm(month)//'01'
         endif
         inquire(file=trim(Path_Input)//trim(DomainName)//'_'//filein  &
                ,exist=there)
         if(.not.there) then
            write(*,*) trim(Path_Input)//trim(DomainName)//'_'//filein &
                      ,' is not avaiable'
            write(*,*) 'Note: Daily mean files are required.'
            stop
         endif
         open(10,file=trim(Path_Input)//trim(DomainName)//'_'//filein &
           ,form='unformatted',recl=(iy-2)*jx_len*ibyte,access='direct')

         if(nfile.eq.1.or.month.eq.1) then
      write(*,*)'Calculate the monthly mean of   ',filename,nyear,month
         else
      write(*,*)'                                ',filename,nyear,month
         endif

         do n=1,nvar
            do i=1,iy-2
            do j=1,jx_len
               c(j,i,n) = 0.0
            enddo
            enddo
         enddo

         do nday=1,nrecord
            do n=1,nvar
               nrec=nrec+1
               read(10,rec=nrec) b
               do i=1,iy-2
               do j=1,jx_len
                  c(j,i,n) = c(j,i,n)+b(j,i)
               enddo
               enddo
            enddo
         enddo
         do n=1,nvar
            do i=1,iy-2
            do j=1,jx_len
               c(j,i,n) = c(j,i,n)/float(nrecord)
            enddo
            enddo
            mrec=mrec+1
            write(20,rec=mrec) ((c(j,i,n),j=1,jx_len),i=1,iy-2)
         enddo
         close(10)
      enddo
      close(20)
 199  format(I4)
      return
      end

      subroutine mon2p(filename,iy,jx,kz,np,plev,i_band,ibyte  &
                     ,Path_Input,DomainName,idate1,idate2,igrads)
      implicit none
      character*6 filename
      character*128 Path_Input
      character*20 DomainName
      integer iy,jx,kz,np,i_band,ibyte,idate1,idate2,igrads
      real*4  plev(np)
      integer iiy,jjx,kkz
      real*4  truelatL,truelatH
      real*4  dsinm,ptop,clat,clon,plat,plon,GRDFAC
      integer jgrads,ibigend
      character*6 iproj
      real*4, allocatable ::  sigma(:)
      character*4 :: chy
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
      character*10 fileout
      integer ntype,nfile,nyear,month,mrec,nrec
      logical there
      real*4  alatmin,alatmax,alonmin,alonmax,rlatinc,rloninc
      real*4  centerj,centeri
      integer ny,nx
      integer i,j,k,n,jx_len
      integer nvar,n_month,nrecord,nday
      real*4, allocatable ::  c(:,:,:),b(:,:),xlat(:,:),xlon(:,:)

      if(i_band.eq.1) then
         jx_len = jx
      else
         jx_len = jx-2
      endif

      allocate(sigma(kz+1))
      allocate(b(jx,iy))
      allocate(xlat(jx_len,iy-2))
      allocate(xlon(jx_len,iy-2))

      inquire(file=trim(Path_Input)//trim(DomainName)//'.INFO'  &
             ,exist=there)
      if(.not.there) then
         write(*,*) trim(Path_Input)//trim(DomainName)//'.INFO'  &
                   ,' is not avaiable'
         stop
      endif
      open(10,file=trim(Path_Input)//trim(DomainName)//'.INFO'  &
             ,form='unformatted',recl=jx*iy*ibyte,access='direct')
      read(10,rec=1) iiy,jjx,kkz,dsinm,clat,clon,plat,plon,GRDFAC  &
                    ,iproj,(sigma(k),k=1,kz+1),ptop,jgrads,ibigend &
                    ,truelatL,truelatH
      if(iiy.ne.iy.or.jjx.ne.jx.or.kkz.ne.kz) then
         write(*,*) 'iy,jx,kz in parameter = ',iy,jx,kz
         write(*,*) 'iy,jx,kz in DOMAIN.INFO ',iiy,jjx,kkz
         write(*,*) 'They are not consistent'
         stop
      endif
      read(10,rec=5) b
      do i=1,iy-2
      do j=1,jx_len
         if(i_band.eq.1) then
            xlat(j,i) = b(j,i+1)
         else
            xlat(j,i) = b(j+1,i+1)
         endif
      enddo
      enddo
      read(10,rec=6) b
      do i=1,iy-2
      do j=1,jx_len
         if(i_band.eq.1) then
            xlon(j,i) = b(j,i+1)
         else
            xlon(j,i) = b(j+1,i+1)
         endif
      enddo
      enddo
      close(10)

      deallocate(b)
      allocate(b(jx_len,iy-2))

      if(filename.eq.'ICBC_P') then
         nvar = np*6+2
      else
         write(*,*) 'filename is not correct'
         stop
      endif
      allocate(c(jx_len,iy-2,nvar))

      if(idate1.ge.1000010100) then        ! original file
         n_month = (idate2/1000000-idate1/1000000)*12 &
                 + (mod(idate2/10000,100)-mod(idate1/10000,100))
         if(mod(idate2,10000).gt.0100) n_month = n_month+1
         ntype = 0
         month = mod(idate1/10000,100)
         nyear = idate1/1000000
      else if(idate1.ge.10000101) then     ! daily mean file
         n_month = (idate2/10000-idate1/10000)*12 &
                 + (mod(idate2/100,100)-mod(idate1/100,100))
         if(mod(idate2,100).gt.01) n_month = n_month+1
         ntype = 1
         month = mod(idate1/100,100)
         nyear = idate1/10000
      else
         write(*,*) 'date should be in 10 digits, or 8 digits'
         stop
      endif

      fileout= 'ICBC_P.mon'

      open(20,file=trim(Path_Input)//trim(DomainName)//'_'//fileout  &
          ,form='unformatted',recl=(iy-2)*jx_len*ibyte,access='direct')
      mrec = 0

      IF(igrads.eq.1) THEN
      inquire(file= &
            trim(Path_Input)//trim(DomainName)//'_'//fileout//'.ctl' &
             ,exist=there)
      if(there) then
         open(31,file= &
             trim(Path_Input)//trim(DomainName)//'_'//fileout//'.ctl' &
                                                   ,status='replace')
      else
         open(31,file= &
             trim(Path_Input)//trim(DomainName)//'_'//fileout//'.ctl' &
                                                       ,status='new')
      endif
      write(31,10) '^'//trim(DomainName)//'_'//fileout
  10  format('dset ',A30)
      write(31,20)
  20  format('title RegCM monthly mean pressure level ICBC variables')
      write(31,30)
  30  format('options big_endian')
      write(31,50)
  50  format('undef -1.e34')
      if(iproj.eq.'LAMCON'.or.iproj.eq.'ROTMER') then
         alatmin= 999999.
         alatmax=-999999.
         do j=1,jx_len
            if(xlat(j,1).lt.alatmin) alatmin=xlat(j,1)
            if(xlat(j,iy-2).gt.alatmax) alatmax=xlat(j,iy-2)
         enddo
         alonmin= 999999.
         alonmax=-999999.
         do i=1,iy-2
         do j=1,jx_len
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
                  if(abs(clon-xlon(j,i)).lt. &
                     abs(clon-(xlon(j,i)-360.))) then
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
         centerj=jx_len/2.
         centeri=(iy-2)/2.
      endif

      if(iproj.eq.'LAMCON') then        ! Lambert projection
         write(31,100) jx_len,iy-2,clat,clon,centerj,centeri, &
                       truelatL,truelatH,clon,dsinm,dsinm
 100  format('pdef ',i4,1x,i4,1x,'lccr',7(1x,f7.2),1x,2(f7.0,1x))
         write(31,110) nx+2,alonmin-rloninc,rloninc
 110  format('xdef ',i4,' linear ',f7.2,1x,f7.4)
         write(31,120) ny+2,alatmin-rlatinc,rlatinc
 120  format('ydef ',i4,' linear ',f7.2,1x,f7.4)
      elseif(iproj.eq.'POLSTR') then    !
      elseif(iproj.eq.'NORMER') then
         write(31,200)  jx_len,xlon(1,1),xlon(2,1)-xlon(1,1)
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
         write(31,230) jx_len,iy-2,plon,plat,dsinm/111000. &
                                          ,dsinm/111000.*.95238
 230  format('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
         write(31,110) nx+2,alonmin-rloninc,rloninc
         write(31,120) ny+2,alatmin-rlatinc,rlatinc
      else
         write(*,*) 'Are you sure your map projection is right ?'
         stop
      endif
      write(31,318) np,(plev(k),k=1,np)
 318  format('zdef ',I2,' levels ',30f7.2)
      if(ntype.eq.0) then
         nyear=idate1/1000000
         month=(idate1-nyear*1000000)/10000
      else if(ntype.eq.1) then
         nyear=idate1/10000
         month=(idate1-nyear*10000)/100
      endif
      write(31,400)n_month,cday(16),chmc(month),nyear
 400  format('tdef ',I4,' linear ',A2,A3,I4,' ','1mo')
      write(31,500) 6+2
 500  format('vars ',I3)
 600  format(A8,'0 99 ',A36)
 611  format(A8,'0 33,105 ',A36)
 612  format(A8,'0 34,105 ',A36)
      if(filename.eq.'ICBC_P') then
 660  format(A8,I2,' 0 ',A36)
 661  format(A8,I2,' 33,100 ',A36)
 662  format(A8,I2,' 34,100 ',A36)
        if(iproj.eq.'LAMCON') then
           write(31,661) 'u       ',np,'westerly wind (m/s)         '
           write(31,662) 'v       ',np,'southerly wind (m/s)        '
        else
           write(31,660) 'u       ',np,'westerly wind (m/s)         '
           write(31,660) 'v       ',np,'southerly wind (m/s)        '
        endif
        write(31,660) 'h       ',np,'geopotential height (m)     '
        write(31,660) 't       ',np,'air temperature (degree)    '
        write(31,660) 'rh      ',np,'relative moisture (%)       '
        write(31,660) 'qv      ',np,'specific moisture (kg/kg)   '
        write(31,600) 'ps      ',   'surface pressure (hPa)      '
        write(31,600) 'slp     ',   'sea level pressure (hPa)    '
      endif
      write(31,700)
 700  format('endvars')
      close(31)
      ENDIF

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

         write(chy,199) nyear
         filein = 'ICBC_P'//chy//chm(month)//'01'
         inquire(file=trim(Path_Input)//trim(DomainName)//'_'//filein  &
                ,exist=there)
         if(.not.there) then
            write(*,*) trim(Path_Input)//trim(DomainName)//'_'//filein &
                      ,' is not avaiable'
            write(*,*) 'Note: Daily mean files are required.'
            stop
         endif
         open(10,file=trim(Path_Input)//trim(DomainName)//'_'//filein &
           ,form='unformatted',recl=(iy-2)*jx_len*ibyte,access='direct')

         if(nfile.eq.1.or.month.eq.1) then
         write(*,*)'Calculate the monthly mean of ',filename,nyear,month
         else
         write(*,*)'                              ',filename,nyear,month
         endif

         do n=1,nvar
            do i=1,iy-2
            do j=1,jx_len
               c(j,i,n) = 0.0
            enddo
            enddo
         enddo

         do nday=1,nrecord
            do n=1,nvar
               nrec=nrec+1
               read(10,rec=nrec) b
               do i=1,iy-2
               do j=1,jx_len
                  c(j,i,n) = c(j,i,n)+b(j,i)
               enddo
               enddo
            enddo
         enddo
         do n=1,nvar
            do i=1,iy-2
            do j=1,jx_len
               c(j,i,n) = c(j,i,n)/float(nrecord)
            enddo
            enddo
            mrec=mrec+1
            write(20,rec=mrec) ((c(j,i,n),j=1,jx_len),i=1,iy-2)
         enddo
         close(10)
      enddo
      close(20)
 199  format(I4)
      return
      end

      subroutine mon3(filename,iy,jx,kz,nsg,i_band,ibyte, Path_Input &
                   ,DomainName,Path_Output,idate1,idate2,igrads)
      implicit none
      character*3 filename
      character*128 Path_Input,Path_Output
      character*20 DomainName
      integer iy,jx,kz,nsg,i_band,ibyte,idate1,idate2,igrads
      integer iiy,jjx,kkz
      integer mdate0,ibltyp,icup,ipptls,iboudy
      real*4  truelatL,truelatH
      real*4  dxsp,ptsp,clat,clon,plat,plon
      real*4  dto,dtb,dtr,dtc
      integer iotyp
      character*6 iproj
      real*4, allocatable ::  sigma(:)
      character*4 :: chy
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
      character*3 chnsg
      character*12 filein
      character*7 fileout
      integer ntype,nfile,nyear,month,mrec,nrec
      integer i,j,k,n,jx_len
      integer nvar,n_month,nrecord,nday
      logical there
      real*4  alatmin,alatmax,alonmin,alonmax,rlatinc,rloninc
      real*4  centerj,centeri
      integer ny,nx

      real*4, allocatable ::  c(:,:,:),b(:,:),xlat(:,:),xlon(:,:)
      real*4, allocatable ::  xlat_s(:,:),xlon_s(:,:),a(:,:)

      if(i_band.eq.1) then
         jx_len = jx
      else
         jx_len = jx-2
      endif

      allocate(sigma(kz+1))
      allocate(b(jx_len*nsg,(iy-2)*nsg))
      allocate(xlat(jx_len,iy-2))
      allocate(xlon(jx_len,iy-2))

      allocate(xlat_s(jx_len*nsg,(iy-2)*nsg))
      allocate(xlon_s(jx_len*nsg,(iy-2)*nsg))
      allocate(a(jx*nsg,iy*nsg))

      if(nsg.lt.10) then
         write(chnsg,71) nsg
      else if(nsg.lt.100) then
         write(chnsg,72) nsg
      else
         write(chnsg,73) nsg
      endif
 71   format('00',I1)
 72   format('0', I2)
 73   format(I3)
      inquire(file=trim(Path_Input)//trim(DomainName)//chnsg//'.INFO' &
                  ,exist=there)
      if(.not.there) then
         write(*,*) trim(Path_Input)//trim(DomainName)//chnsg//'.INFO' &
                   ,' is not avaiable'
         stop
      endif
      open(10,file=trim(Path_Input)//trim(DomainName)//chnsg//'.INFO' &
       ,form='unformatted',recl=jx*nsg*iy*nsg*ibyte,access='direct')
      read(10,rec=5) a
      do i=nsg+1,(iy-1)*nsg
      do j=nsg+1,(jx-1)*nsg
         xlat_s(j-nsg,i-nsg) = a(j,i)
      enddo
      enddo
      read(10,rec=6) a
      do i=nsg+1,(iy-1)*nsg
      do j=nsg+1,(jx-1)*nsg
         xlon_s(j-nsg,i-nsg) = a(j,i)
      enddo
      enddo
      close(10)
      deallocate(a)

      inquire(file=trim(Path_Output)//'OUT_HEAD',exist=there)
      if(.not.there) then
         write(*,*) trim(Path_Output)//'OUT_HEAD',' is not avaiable'
         stop
      endif
      open(10,file=trim(Path_Output)//'OUT_HEAD',form='unformatted' &
             ,recl=jx_len*(iy-2)*ibyte,access='direct')
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
      allocate(c(jx_len*nsg,(iy-2)*nsg,nvar))

      if(idate1.ge.1000010100) then        ! original file
         n_month = (idate2/1000000-idate1/1000000)*12 &
                 + (mod(idate2/10000,100)-mod(idate1/10000,100))
         if(mod(idate2,10000).gt.0100) n_month = n_month+1
         ntype = 0
         month = mod(idate1/10000,100)
         nyear = idate1/1000000
      else if(idate1.ge.10000101) then     ! daily mean file
         n_month = (idate2/10000-idate1/10000)*12 &
                 + (mod(idate2/100,100)-mod(idate1/100,100))
         if(mod(idate2,100).gt.01) n_month = n_month+1
         ntype = 1
         month = mod(idate1/100,100)
         nyear = idate1/10000
      else
         write(*,*) 'date should be in 10 digits, or 8 digits'
         stop
      endif

      if(filename.eq.'SUB') then
         fileout= 'SUB.mon'
      endif

      open(20,file=trim(Path_Output)//fileout,form='unformatted' &
             ,recl=(iy-2)*nsg*jx_len*nsg*ibyte,access='direct')
      mrec = 0

      IF(igrads.eq.1) THEN
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
         do j=1,jx_len
            if(xlat(j,1).lt.alatmin) alatmin=xlat(j,1)
            if(xlat(j,iy-2).gt.alatmax) alatmax=xlat(j,iy-2)
         enddo
         alonmin= 999999.
         alonmax=-999999.
         do i=1,iy-2
         do j=1,jx_len
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
                  if(abs(clon-xlon(j,i)).lt. &
                     abs(clon-(xlon(j,i)-360.))) then
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
         centerj=jx_len/2.
         centeri=(iy-2)/2.
      endif

      if(iproj.eq.'LAMCON') then        ! Lambert projection
         write(31,100) jx_len*nsg,(iy-2)*nsg,clat,clon,  &
                       centerj*nsg,centeri*nsg, &
               truelatL,truelatH,clon,dxsp*1000./nsg,dxsp*1000./nsg
 100  format('pdef ',i4,1x,i4,1x,'lccr',7(1x,f7.2),1x,2(f7.0,1x))
         write(31,110) nx+2,alonmin-rloninc,rloninc
 110  format('xdef ',i4,' linear ',f7.2,1x,f7.4)
         write(31,120) ny+2,alatmin-rlatinc,rlatinc
 120  format('ydef ',i4,' linear ',f7.2,1x,f7.4)
      elseif(iproj.eq.'POLSTR') then    !
      elseif(iproj.eq.'NORMER') then
         write(31,200)  jx_len*nsg,xlon_s(1,1),xlon_s(2,1)-xlon_s(1,1)
 200  format('xdef ',I3,' linear ',f9.4,' ',f9.4)
         write(31,210) (iy-2)
 210  format('ydef ',I3,' levels')
         write(31,220) (xlat_s(1,i),i=1,(iy-2)*nsg)
 220  format(10f7.2)
      elseif(iproj.eq.'ROTMER') then
         write(*,*) 'Note that rotated Mercartor (ROTMER)' &
                   ,' projections are not supported by GrADS.'
         write(*,*) '  Although not exact, the eta.u projection' &
                   ,' in GrADS is somewhat similar.'
         write(*,*) ' FERRET, however, does support this projection.'
         write(31,230) jx_len*nsg,(iy-2)*nsg,plon,plat, &
                       dxsp/111./nsg ,dxsp/111.*.95238/nsg
 230  format('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
         write(31,110) nx+2,alonmin-rloninc,rloninc
         write(31,120) ny+2,alatmin-rlatinc,rlatinc
      else
         write(*,*) 'Are you sure your map projection is right ?'
         stop
      endif
      write(31,300) (1013.25-ptsp*10.)*(sigma(1)+sigma(2))*.5+ptsp*10.
 300  format('zdef 1',' levels ',f7.2)
      if(ntype.eq.0) then
         nyear=idate1/1000000
         month=(idate1-nyear*1000000)/10000
      else if(ntype.eq.1) then
         nyear=idate1/10000
         month=(idate1-nyear*10000)/100
      endif
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
        write(31,600) 'scv     ','snow amount (mm, water equivalent)   '
        write(31,600) 'sena    ','sensible heat flux (W/m2)            '
        write(31,600) 'prcv    ','convective precipitation (mm/day)    '
        write(31,600) 'ps      ','surface pressure (hPa)               '
        write(31,700)
 700  format('endvars')
      endif
      close(31)
      ENDIF

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
      write(*,*)'Calculate the monthly mean of    ',filename,nyear,month
         else
      write(*,*)'                                 ',filename,nyear,month
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

         write(chy,199) nyear
         if(filename.eq.'SUB') then
            filein = 'SUB.'//chy//chm(month)//'01'
         endif
         inquire(file=trim(Path_Output)//filein,exist=there)
         if(.not.there) then
            write(*,*) trim(Path_Output)//filein,' is not avaiable'
            write(*,*) 'Note: Daily mean files are required.'
            stop
         endif
         open(10,file=trim(Path_Output)//filein,form='unformatted' &
                ,recl=(iy-2)*nsg*jx_len*nsg*ibyte,access='direct')
         nrec = 0

         do n=1,nvar
            do i=1,(iy-2)*nsg
            do j=1,jx_len*nsg
               c(j,i,n) = 0.0
            enddo
            enddo
         enddo

         do nday=1,nrecord
            do n=1,nvar
               nrec=nrec+1
               read(10,rec=nrec) b
               do i=1,(iy-2)*nsg
               do j=1,jx_len*nsg
                  c(j,i,n) = c(j,i,n)+b(j,i)
               enddo
               enddo
            enddo
         enddo
         do n=1,nvar
            do i=1,(iy-2)*nsg
            do j=1,jx_len*nsg
               c(j,i,n) = c(j,i,n)/float(nrecord)
            enddo
            enddo
            mrec=mrec+1
            write(20,rec=mrec)((c(j,i,n),j=1,jx_len*nsg),i=1,(iy-2)*nsg)
         enddo
         close(10)
      enddo
      close(20)
      deallocate(sigma)
      deallocate(b)
      deallocate(xlat)
      deallocate(xlon)
      deallocate(c)

      deallocate(xlat_s)
      deallocate(xlon_s)
 199  format(I4)
      return
      end

C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C    This file is part of RegCM model.
C
C    RegCM model is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    RegCM model is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
C
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      subroutine GRADSBAT(CTLNAME)
      implicit none
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#     include "include/date.cb"
#     include "include/param1.cb"
#     include "include/param2.cb"
#     include "include/param3.cb"
#ifdef MPP1
#     include "include/io.cb"
      include 'mpiregcm.h'
      integer ierr
#else
#     include "include/main.cb"
#endif

      real(kind=4)  clat,clon,dxsp,ptsp,plat,plon
      character proj*6
      COMMON /GRADS/clat,clon,dxsp,ptsp,plat,plon,proj
      real(kind=4)  TRUELATL,TRUELATH
      COMMON /SAVEPAR/ TRUELATL,TRUELATH
      character*18 CTLNAME
      logical there
      integer i,j,ny,nx
      real(kind=4)  alatmin,alatmax,alonmin,alonmax,rlatinc,rloninc
      real(kind=4)  centerj,centeri
      integer myear,month,nday,nhour
      integer number,nbase,mnend,jbend,ifrq
      character*2 cday(31)
      data cday/'01','02','03','04','05','06','07','08','09','10',
     &          '11','12','13','14','15','16','17','18','19','20',
     &     '21','22','23','24','25','26','27','28','29','30','31'/
      character*3 cmonth(12)
      data cmonth/'jan','feb','mar','apr','may','jun',
     &            'jul','aug','sep','oct','nov','dec'/
      integer isystm,system
      external system
c
      inquire(file='output/'//CTLNAME,exist=there)
      if(there) isystm=system('/bin/rm output/'//CTLNAME)
      open(31,file='output/'//CTLNAME,status='new')
      write(31,10) CTLNAME(1:14)
  10  format('dset ^',A14)
      write(31,20)
  20  format('title RegCM normal output variables')
      if(ibintyp.eq.1) then
         write(31,30)
  30     format('options big_endian')
      else if(ibintyp.eq.2) then
         write(31,40)
  40     format('options little_endian')
      endif
      if(iotyp.eq.2) write(31,'(a)') 'options sequential'
      write(31,50)
  50  format('undef -1.e34')
      if(proj.eq.'LAMCON'.or.proj.eq.'ROTMER') then
         alatmin= 999999.
         alatmax=-999999.
#ifdef MPP1
         do j=1,mjx
            if(xlat_io(1 ,j).lt.alatmin)alatmin=xlat_io(1 ,j)
            if(xlat_io(ix,j).gt.alatmax)alatmax=xlat_io(ix,j)
#else
         do j=1,jx
            if(xlat(1 ,j).lt.alatmin)alatmin=xlat(1 ,j)
            if(xlat(ix,j).gt.alatmax)alatmax=xlat(ix,j)
#endif
         enddo
         alonmin= 999999.
         alonmax=-999999.
         do i=1,ix
#ifdef MPP1
         do j=1,mjx
#else
         do j=1,jx
#endif
            if(clon.ge.0.0) then
#ifdef MPP1
               if(xlong_io(i,j).ge.0.0) then
                  alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
                  alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
#else
               if(xlong(i,j).ge.0.0) then
                  alonmin = amin1(alonmin,sngl(xlong(i,j)))
                  alonmax = amax1(alonmax,sngl(xlong(i,j)))
#endif
               else
#ifdef MPP1
                  if(abs(clon-xlong_io(i,j)).lt.
     &               abs(clon-(xlong_io(i,j)+360.))) then
                     alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
                     alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
#else
                  if(abs(clon-xlong(i,j)).lt.
     &               abs(clon-(xlong(i,j)+360.))) then
                     alonmin = amin1(alonmin,sngl(xlong(i,j)))
                     alonmax = amax1(alonmax,sngl(xlong(i,j)))
#endif
                  else
#ifdef MPP1
                     alonmin = amin1(alonmin,sngl(xlong_io(i,j))+360.)
                     alonmax = amax1(alonmax,sngl(xlong_io(i,j))+360.)
#else
                     alonmin = amin1(alonmin,sngl(xlong(i,j))+360.)
                     alonmax = amax1(alonmax,sngl(xlong(i,j))+360.)
#endif
                  endif
               endif
            else
#ifdef MPP1
               if(xlong_io(i,j).lt.0.0) then
                  alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
                  alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
#else
               if(xlong(i,j).lt.0.0) then
                  alonmin = amin1(alonmin,sngl(xlong(i,j)))
                  alonmax = amax1(alonmax,sngl(xlong(i,j)))
#endif
               else
#ifdef MPP1
                  if(abs(clon-xlong_io(i,j)).lt.
     &               abs(clon-(xlong_io(i,j)-360.))) then
                     alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
                     alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
#else
                  if(abs(clon-xlong(i,j)).lt.
     &               abs(clon-(xlong(i,j)-360.))) then
                     alonmin = amin1(alonmin,sngl(xlong(i,j)))
                     alonmax = amax1(alonmax,sngl(xlong(i,j)))
#endif
                  else
#ifdef MPP1
                     alonmin = amin1(alonmin,sngl(xlong_io(i,j))-360.)
                     alonmax = amax1(alonmax,sngl(xlong_io(i,j))-360.)
#else
                     alonmin = amin1(alonmin,sngl(xlong(i,j))-360.)
                     alonmax = amax1(alonmax,sngl(xlong(i,j))-360.)
#endif
                  endif
               endif
            endif
         enddo
         enddo
         rlatinc=dxsp*0.001/111./2.
         rloninc=dxsp*0.001/111./2.
         ny=2+nint(abs(alatmax-alatmin)/rlatinc)
         nx=1+nint(abs((alonmax-alonmin)/rloninc))

#ifdef MPP1
         centerj=(mjx-2)/2.
#else
         centerj=(jx-2)/2.
#endif
         centeri=(ix-2)/2.
      endif
      if(iotyp.eq.1) then
      if(proj.eq.'LAMCON') then        ! Lambert projection
#ifdef MPP1
         write(31,100) mjx-2,ix-2,clat,clon,centerj,centeri,
     &                 truelatL,truelatH,clon,dxsp,dxsp
#else
         write(31,100) jx-2,ix-2,clat,clon,centerj,centeri,
     &                 truelatL,truelatH,clon,dxsp,dxsp
#endif
 100  format('pdef ',i4,1x,i4,1x,'lccr',7(1x,f7.2),1x,2(f7.0,1x))
         write(31,110) nx+2,alonmin-rloninc,rloninc
 110  format('xdef ',i4,' linear ',f7.2,1x,f7.4)
         write(31,120) ny+2,alatmin-rlatinc,rlatinc
 120  format('ydef ',i4,' linear ',f7.2,1x,f7.4)
      elseif(proj.eq.'POLSTR') then    !
      elseif(proj.eq.'NORMER') then
#ifdef MPP1
         write(31,200)  mjx-2,xlong_io(2,2),xlong_io(2,3)-xlong_io(2,2)
#else
         write(31,200)  jx-2,xlong(2,2),xlong(2,3)-xlong(2,2)
#endif
 200  format('xdef ',I3,' linear ',f9.4,' ',f9.4)
         write(31,210) ix-2
 210  format('ydef ',I3,' levels')
#ifdef MPP1
         write(31,220) (xlat_io(i,2),i=2,ix-1)
#else
         write(31,220) (xlat(i,2),i=2,ix-1)
#endif
 220  format(10f7.2)
      elseif(proj.eq.'ROTMER') then
         write(*,*) 'Note that rotated Mercartor (ROTMER)'
     &             ,' projections are not supported by GrADS.'
         write(*,*) '  Although not exact, the eta.u projection'
     &             ,' in GrADS is somewhat similar.'
         write(*,*) ' FERRET, however, does support this projection.'
#ifdef MPP1
         write(31,230) mjx-2,ix-2,plon,plat,dxsp/111000.
     &                                ,dxsp/111000.*.95238
#else
         write(31,230) jx-2,ix-2,plon,plat,dxsp/111000.
     &                                ,dxsp/111000.*.95238
#endif
 230  format('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
         write(31,110) nx+2,alonmin-rloninc,rloninc
         write(31,120) ny+2,alatmin-rlatinc,rlatinc
      else
         write(*,*) 'Are you sure your map projection is right ?'
#ifdef MPP1
         call MPI_ABORT(MPI_COMM_WORLD,ierr)
#endif
         stop
      endif
      else if(iotyp.eq.2) then
#ifdef MPP1
         write(31,200) mjx-2,xlong_io(2,2),xlong_io(2,3)-xlong_io(2,2)
         write(31,201) ix-2,xlat_io(2,2),xlat_io(3,2)-xlat_io(2,2)
#else
         write(31,200) jx-2,xlong(2,2),xlong(2,3)-xlong(2,2)
         write(31,201) ix-2,xlat(2,2),xlat(3,2)-xlat(2,2)
#endif
      endif
 201  format('ydef ',I3,' linear ',f9.4,' ',f9.4)
      write(31,300) (1013.25-ptop*10.)*a(kx)+ptop*10.
 300  format('zdef 1',' levels ',f7.2)
      myear=ldatez/1000000
      month=(ldatez-myear*1000000)/10000
      nday =(ldatez-myear*1000000-month*10000)/100
      nhour=MOD(ldatez,100)
      CALL FINDDATE(nbase,ldatez)
      if(month.eq.12) then
         CALL FINDDATE(mnend,myear*1000000+1010100)
      else
         CALL FINDDATE(mnend,myear*1000000+month*10000+10100)
      endif
      CALL FINDDATE(jbend,IDATE2)
      if(ldatez.EQ.IDATE0) then
         number=(ibdyfrq/batfrq+0.00001)*(min0(jbend,mnend)-nbase)+1
      else
         number=(ibdyfrq/batfrq+0.00001)*(min0(jbend,mnend)-nbase)
      endif
      ifrq=batfrq+0.00001
      if(ldatez.EQ.IDATE0) then
         write(31,400)number,nhour     ,cday(nday),cmonth(month)
     &               ,myear,ifrq
      else
         write(31,400)number,nhour+ifrq,cday(nday),cmonth(month)
     &               ,myear,ifrq
      endif
 400  format('tdef ',I4,' linear ',I2,'z',A2,A3,I4,' ',I2,'hr')
      if(iotyp.eq.2)  write(31,'(a)') 'theader 4'
      write(31,500) 21+6
 500  format('vars ',I2)
 600  format(A8,'0 99 ',A36)
 611  format(A8,'0 33,105 ',A36)
 612  format(A8,'0 34,105 ',A36)
      if(proj.eq.'LAMCON'.and.iotyp.eq.1) then        ! Lambert projection
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
      write(31,700)
 700  format('endvars')
      close(31)
      return
      end
      subroutine GRADSSUB(CTLNAME)
      implicit none
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#     include "include/date.cb"
#     include "include/iunits.cb"
#     include "include/param1.cb"
#     include "include/param2.cb"
#     include "include/param3.cb"
#ifdef MPP1
#     include "include/io.cb"
      include 'mpiregcm.h'
      INTEGER ierr
      real(kind=4)  xlat_s_io(mjx*NSG,ix*NSG),xlon_s_io(mjx*NSG,ix*NSG)
#else
#     include "include/main.cb"
      real(kind=4)  xlat_s(jx*NSG,ix*NSG),xlon_s(jx*NSG,ix*NSG)
#endif
      real(kind=4)  clat,clon,dxsp,ptsp,plat,plon
      character proj*6
      COMMON /GRADS/clat,clon,dxsp,ptsp,plat,plon,proj
      real(kind=4)  TRUELATL,TRUELATH
      COMMON /SAVEPAR/ TRUELATL,TRUELATH
      character*18 CTLNAME
      logical there
      integer i,j,ny,nx
      real(kind=4)  alatmin,alatmax,alonmin,alonmax,rlatinc,rloninc
      real(kind=4)  centerj,centeri
      integer myear,month,nday,nhour
      integer number,nbase,mnend,jbend,ifrq
      character*2 cday(31)
      data cday/'01','02','03','04','05','06','07','08','09','10',
     &          '11','12','13','14','15','16','17','18','19','20',
     &     '21','22','23','24','25','26','27','28','29','30','31'/
      character*3 cmonth(12)
      data cmonth/'jan','feb','mar','apr','may','jun',
     &            'jul','aug','sep','oct','nov','dec'/
      integer isystm,system
      external system
c
      inquire(file='output/'//CTLNAME,exist=there)
      if(there) isystm=system('/bin/rm output/'//CTLNAME)
      open(31,file='output/'//CTLNAME,status='new')
      write(31,10) CTLNAME(1:14)
  10  format('dset ^',A14)
      write(31,20)
  20  format('title RegCM normal output variables')
      if(ibintyp.eq.1) then
         write(31,30)
  30     format('options big_endian')
      else if(ibintyp.eq.2) then
         write(31,40)
  40     format('options little_endian')
      endif
      if(iotyp.eq.2) write(31,'(a)') 'options sequential'
      write(31,50)
  50  format('undef -1.e34')
      if(proj.eq.'LAMCON'.or.proj.eq.'ROTMER') then
         alatmin= 999999.
         alatmax=-999999.
#ifdef MPP1
         do j=1,mjx
            if(xlat_io(1 ,j).lt.alatmin)alatmin=xlat_io(1 ,j)
            if(xlat_io(ix,j).gt.alatmax)alatmax=xlat_io(ix,j)
#else
         do j=1,jx
            if(xlat(1 ,j).lt.alatmin)alatmin=xlat(1 ,j)
            if(xlat(ix,j).gt.alatmax)alatmax=xlat(ix,j)
#endif
         enddo
         alonmin= 999999.
         alonmax=-999999.
         do i=1,ix
#ifdef MPP1
         do j=1,mjx
#else
         do j=1,jx
#endif
            if(clon.ge.0.0) then
#ifdef MPP1
               if(xlong_io(i,j).ge.0.0) then
                  alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
                  alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
#else
               if(xlong(i,j).ge.0.0) then
                  alonmin = amin1(alonmin,sngl(xlong(i,j)))
                  alonmax = amax1(alonmax,sngl(xlong(i,j)))
#endif
               else
#ifdef MPP1
                  if(abs(clon-xlong_io(i,j)).lt.
     &               abs(clon-(xlong_io(i,j)+360.))) then
                     alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
                     alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
#else
                  if(abs(clon-xlong(i,j)).lt.
     &               abs(clon-(xlong(i,j)+360.))) then
                     alonmin = amin1(alonmin,sngl(xlong(i,j)))
                     alonmax = amax1(alonmax,sngl(xlong(i,j)))
#endif
                  else
#ifdef MPP1
                     alonmin = amin1(alonmin,sngl(xlong_io(i,j))+360.)
                     alonmax = amax1(alonmax,sngl(xlong_io(i,j))+360.)
#else
                     alonmin = amin1(alonmin,sngl(xlong(i,j))+360.)
                     alonmax = amax1(alonmax,sngl(xlong(i,j))+360.)
#endif
                  endif
               endif
            else
#ifdef MPP1
               if(xlong_io(i,j).lt.0.0) then
                  alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
                  alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
#else
               if(xlong(i,j).lt.0.0) then
                  alonmin = amin1(alonmin,sngl(xlong(i,j)))
                  alonmax = amax1(alonmax,sngl(xlong(i,j)))
#endif
               else
#ifdef MPP1
                  if(abs(clon-xlong_io(i,j)).lt.
     &               abs(clon-(xlong_io(i,j)-360.))) then
                     alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
                     alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
#else
                  if(abs(clon-xlong(i,j)).lt.
     &               abs(clon-(xlong(i,j)-360.))) then
                     alonmin = amin1(alonmin,sngl(xlong(i,j)))
                     alonmax = amax1(alonmax,sngl(xlong(i,j)))
#endif
                  else
#ifdef MPP1
                     alonmin = amin1(alonmin,sngl(xlong_io(i,j))-360.)
                     alonmax = amax1(alonmax,sngl(xlong_io(i,j))-360.)
#else
                     alonmin = amin1(alonmin,sngl(xlong(i,j))-360.)
                     alonmax = amax1(alonmax,sngl(xlong(i,j))-360.)
#endif
                  endif
               endif
            endif
         enddo
         enddo
         rlatinc=dxsp*0.001/111./2.
         rloninc=dxsp*0.001/111./2.
         ny=2+nint(abs(alatmax-alatmin)/rlatinc)
         nx=1+nint(abs((alonmax-alonmin)/rloninc))

#ifdef MPP1
         centerj=(mjx-2)/2.
#else
         centerj=(jx-2)/2.
#endif
         centeri=(ix-2)/2.
      endif
      if(iotyp.eq.1) then
      if(proj.eq.'LAMCON') then        ! Lambert projection
#ifdef MPP1
         write(31,100) (mjx-2)*NSG,(ix-2)*NSG,clat,clon,centerj*NSG,
     &        centeri*NSG,truelatL,truelatH,clon,dxsp/NSG,dxsp/NSG
#else
         write(31,100) (jx-2)*NSG,(ix-2)*NSG,clat,clon,centerj*NSG,
     &        centeri*NSG,truelatL,truelatH,clon,dxsp/NSG,dxsp/NSG
#endif
 100  format('pdef ',i4,1x,i4,1x,'lccr',7(1x,f7.2),1x,2(f7.0,1x))
         write(31,110) nx+2,alonmin-rloninc,rloninc
 110  format('xdef ',i4,' linear ',f7.2,1x,f7.4)
         write(31,120) ny+2,alatmin-rlatinc,rlatinc
 120  format('ydef ',i4,' linear ',f7.2,1x,f7.4)
      elseif(proj.eq.'POLSTR') then    !
      elseif(proj.eq.'NORMER') then
#ifdef MPP1
         read(iutin1,rec=5) xlat_s_io
         read(iutin1,rec=6) xlon_s_io
         write(31,200) (mjx-2)*NSG,xlon_s_io(NSG,NSG),
     &                  xlon_s_io(NSG+1,NSG)-xlon_s_io(NSG,NSG)
#else
         read(iutin1,rec=5) xlat_s
         read(iutin1,rec=6) xlon_s
         write(31,200) (jx-2)*NSG,xlon_s(NSG,NSG),
     &                            xlon_s(NSG+1,NSG)-xlon_s(NSG,NSG)
#endif
 200  format('xdef ',I3,' linear ',f9.4,' ',f9.4)
         write(31,210) (ix-2)*NSG
 210  format('ydef ',I3,' levels')
#ifdef MPP1
         write(31,220) (xlat_s_io(NSG+1,i),i=NSG+1,(ix-1)*NSG)
#else
         write(31,220) (xlat_s(NSG+1,i),i=NSG+1,(ix-1)*NSG)
#endif
 220  format(10f7.2)
      elseif(proj.eq.'ROTMER') then
         write(*,*) 'Note that rotated Mercartor (ROTMER)'
     &             ,' projections are not supported by GrADS.'
         write(*,*) '  Although not exact, the eta.u projection'
     &             ,' in GrADS is somewhat similar.'
         write(*,*) ' FERRET, however, does support this projection.'
#ifdef MPP1
         write(31,230) (mjx-2)*NSG,(ix-2)*NSG,plon,plat,
     &                 dxsp/111000./NSG,dxsp/111000.*.95238/NSG
#else
         write(31,230) (jx-2)*NSG,(ix-2)*NSG,plon,plat,
     &                 dxsp/111000./NSG,dxsp/111000.*.95238/NSG
#endif
 230  format('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
         write(31,110) nx+2,alonmin-rloninc,rloninc
         write(31,120) ny+2,alatmin-rlatinc,rlatinc
      else
         write(*,*) 'Are you sure your map projection is right ?'
#ifdef MPP1
         call MPI_ABORT(MPI_COMM_WORLD,ierr)
#endif
         stop
      endif
      else if(iotyp.eq.2) then
#ifdef MPP1
         write(31,200) (mjx-2)*NSG,xlong_io(2,2),
     &                          (xlong_io(2,3)-xlong_io(2,2))/NSG
         write(31,201) (ix-2)*NSG,xlat_io(2,2),
     &                          (xlat_io(3,2)-xlat_io(2,2))/NSG
#else
         write(31,200) (jx-2)*NSG,xlong(2,2),(xlong(2,3)-xlong(2,2))/NSG
         write(31,201) (ix-2)*NSG,xlat(2,2),(xlat(3,2)-xlat(2,2))/NSG
#endif
      endif
 201  format('ydef ',I3,' linear ',f9.4,' ',f9.4)
      write(31,300) 1,(1013.25-ptop*10.)*a(kx)+ptop*10.
 300  format('zdef ',I2,' levels ',30f7.2)
      myear=ldatez/1000000
      month=(ldatez-myear*1000000)/10000
      nday =(ldatez-myear*1000000-month*10000)/100
      nhour=MOD(ldatez,100)
      CALL FINDDATE(nbase,ldatez)
      if(month.eq.12) then
         CALL FINDDATE(mnend,myear*1000000+1010100)
      else
         CALL FINDDATE(mnend,myear*1000000+month*10000+10100)
      endif
      CALL FINDDATE(jbend,IDATE2)
      if(ldatez.EQ.IDATE0) then
         number=(ibdyfrq/batfrq+0.00001)*(min0(jbend,mnend)-nbase)+1
      else
         number=(ibdyfrq/batfrq+0.00001)*(min0(jbend,mnend)-nbase)
      endif
      ifrq=batfrq+0.00001
      if(ldatez.EQ.IDATE0) then
         write(31,400)number,nhour     ,cday(nday),cmonth(month)
     &               ,myear,ifrq
      else
         write(31,400)number,nhour+ifrq,cday(nday),cmonth(month)
     &               ,myear,ifrq
      endif
 400  format('tdef ',I4,' linear ',I2,'z',A2,A3,I4,' ',I2,'hr')
      if(iotyp.eq.2)  write(31,'(a)') 'theader 4'
      write(31,500) 16
 500  format('vars ',I2)
 600  format(A8,'0 99 ',A36)
 611  format(A8,'0 33,105 ',A36)
 612  format(A8,'0 34,105 ',A36)
      if(proj.eq.'LAMCON'.and.iotyp.eq.1) then        ! Lambert projection
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
      write(31,600) 'ps      ','surface pressure (hPa)               '
      write(31,700)
 700  format('endvars')
      close(31)
      return
      end
      subroutine GRADSCHEM(CTLNAME)
      implicit none
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#     include "include/date.cb"
#     include "include/param1.cb"
#     include "include/param2.cb"
#     include "include/param3.cb"
#     include "include/trachem.cb"
#ifdef MPP1
#     include "include/io.cb"
      include 'mpiregcm.h'
      INTEGER ierr
#else
#     include "include/main.cb"
#endif
      real(kind=4)  clat,clon,dxsp,ptsp,plat,plon
      character proj*6
      COMMON /GRADS/clat,clon,dxsp,ptsp,plat,plon,proj
      real(kind=4)  TRUELATL,TRUELATH
      COMMON /SAVEPAR/ TRUELATL,TRUELATH
      character*18 CTLNAME
      logical there
      integer i,j,itr,ny,nx
      real(kind=4)  alatmin,alatmax,alonmin,alonmax,rlatinc,rloninc
      real(kind=4)  centerj,centeri
      integer myear,month,nday,nhour,k
      integer number,nbase,mnend,jbend,ifrq
      character*2 cday(31)
      data cday/'01','02','03','04','05','06','07','08','09','10',
     &          '11','12','13','14','15','16','17','18','19','20',
     &     '21','22','23','24','25','26','27','28','29','30','31'/
      character*3 cmonth(12)
      data cmonth/'jan','feb','mar','apr','may','jun',
     &            'jul','aug','sep','oct','nov','dec'/
      integer isystm,system
      external system
c
      inquire(file='output/'//CTLNAME,exist=there)
      if(there) isystm=system('/bin/rm output/'//CTLNAME)
      open(31,file='output/'//CTLNAME,status='new')
      write(31,10) CTLNAME(1:14)
  10  format('dset ^',A14)
      write(31,20)
  20  format('title RegCM chemistry/tracor variables')
      if(ibintyp.eq.1) then
         write(31,30)
  30     format('options big_endian')
      elseif(ibintyp.eq.2) then
         write(31,40)
  40     format('options little_endian')
      endif
      if(iotyp.eq.2) write(31,'(a)') 'options sequential'
      write(31,50)
  50  format('undef -1.e34')
      if(proj.eq.'LAMCON'.or.proj.eq.'ROTMER') then
         alatmin= 999999.
         alatmax=-999999.
#ifdef MPP1
         do j=1,mjx
            if(xlat_io(1 ,j).lt.alatmin)alatmin=xlat_io(1 ,j)
            if(xlat_io(ix,j).gt.alatmax)alatmax=xlat_io(ix,j)
#else
         do j=1,jx
            if(xlat(1 ,j).lt.alatmin)alatmin=xlat(1 ,j)
            if(xlat(ix,j).gt.alatmax)alatmax=xlat(ix,j)
#endif
         enddo
         alonmin= 999999.
         alonmax=-999999.
         do i=1,ix
#ifdef MPP1
         do j=1,mjx
#else
         do j=1,jx
#endif
            if(clon.ge.0.0) then
#ifdef MPP1
               if(xlong_io(i,j).ge.0.0) then
                  alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
                  alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
#else
               if(xlong(i,j).ge.0.0) then
                  alonmin = amin1(alonmin,sngl(xlong(i,j)))
                  alonmax = amax1(alonmax,sngl(xlong(i,j)))
#endif
               else
#ifdef MPP1
                  if(abs(clon-xlong_io(i,j)).lt.
     &               abs(clon-(xlong_io(i,j)+360.))) then
                     alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
                     alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
#else
                  if(abs(clon-xlong(i,j)).lt.
     &               abs(clon-(xlong(i,j)+360.))) then
                     alonmin = amin1(alonmin,sngl(xlong(i,j)))
                     alonmax = amax1(alonmax,sngl(xlong(i,j)))
#endif
                  else
#ifdef MPP1
                     alonmin = amin1(alonmin,sngl(xlong_io(i,j))+360.)
                     alonmax = amax1(alonmax,sngl(xlong_io(i,j))+360.)
#else
                     alonmin = amin1(alonmin,sngl(xlong(i,j))+360.)
                     alonmax = amax1(alonmax,sngl(xlong(i,j))+360.)
#endif
                  endif
               endif
            else
#ifdef MPP1
               if(xlong_io(i,j).lt.0.0) then
                  alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
                  alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
#else
               if(xlong(i,j).lt.0.0) then
                  alonmin = amin1(alonmin,sngl(xlong(i,j)))
                  alonmax = amax1(alonmax,sngl(xlong(i,j)))
#endif
               else
#ifdef MPP1
                  if(abs(clon-xlong_io(i,j)).lt.
     &               abs(clon-(xlong_io(i,j)-360.))) then
                     alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
                     alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
#else
                  if(abs(clon-xlong(i,j)).lt.
     &               abs(clon-(xlong(i,j)-360.))) then
                     alonmin = amin1(alonmin,sngl(xlong(i,j)))
                     alonmax = amax1(alonmax,sngl(xlong(i,j)))
#endif
                  else
#ifdef MPP1
                     alonmin = amin1(alonmin,sngl(xlong_io(i,j))-360.)
                     alonmax = amax1(alonmax,sngl(xlong_io(i,j))-360.)
#else
                     alonmin = amin1(alonmin,sngl(xlong(i,j))-360.)
                     alonmax = amax1(alonmax,sngl(xlong(i,j))-360.)
#endif
                  endif
               endif
            endif
         enddo
         enddo
         rlatinc=dxsp*0.001/111./2.
         rloninc=dxsp*0.001/111./2.
         ny=2+nint(abs(alatmax-alatmin)/rlatinc)
         nx=1+nint(abs((alonmax-alonmin)/rloninc))

#ifdef MPP1
         centerj=(mjx-2)/2.
#else
         centerj=(jx-2)/2.
#endif
         centeri=(ix-2)/2.
      endif
      if(iotyp.eq.1) then
      if(proj.eq.'LAMCON') then        ! Lambert projection
#ifdef MPP1
         write(31,100) mjx-2,ix-2,clat,clon,centerj,centeri,
     &                 truelatL,truelatH,clon,dxsp,dxsp
#else
         write(31,100) jx-2,ix-2,clat,clon,centerj,centeri,
     &                 truelatL,truelatH,clon,dxsp,dxsp
#endif
 100  format('pdef ',i4,1x,i4,1x,'lcc',7(1x,f7.2),1x,2(f7.0,1x))
         write(31,110) nx+2,alonmin-rloninc,rloninc
 110  format('xdef ',i4,' linear ',f7.2,1x,f7.4)
         write(31,120) ny+2,alatmin-rlatinc,rlatinc
 120  format('ydef ',i4,' linear ',f7.2,1x,f7.4)
      elseif(proj.eq.'POLSTR') then    !
      elseif(proj.eq.'NORMER') then
#ifdef MPP1
         write(31,200)  mjx-2,xlong_io(2,2),xlong_io(2,3)-xlong_io(2,2)
#else
         write(31,200)  jx-2,xlong(2,2),xlong(2,3)-xlong(2,2)
#endif
 200  format('xdef ',I3,' linear ',f9.4,' ',f9.4)
         write(31,210) ix-2
 210  format('ydef ',I3,' levels')
#ifdef MPP1
         write(31,220) (xlat_io(i,2),i=2,ix-1)
#else
         write(31,220) (xlat(i,2),i=2,ix-1)
#endif
 220  format(10f7.2)
      elseif(proj.eq.'ROTMER') then
         write(*,*) 'Note that rotated Mercartor (ROTMER)'
     &             ,' projections are not supported by GrADS.'
         write(*,*) '  Although not exact, the eta.u projection'
     &             ,' in GrADS is somewhat similar.'
         write(*,*) ' FERRET, however, does support this projection.'
#ifdef MPP1
         write(31,230) mjx-2,ix-2,plon,plat,dxsp/111000.
     &                                ,dxsp/111000.*.95238
#else
         write(31,230) jx-2,ix-2,plon,plat,dxsp/111000.
     &                                ,dxsp/111000.*.95238
#endif
 230  format('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
         write(31,110) nx+2,alonmin-rloninc,rloninc
         write(31,120) ny+2,alatmin-rlatinc,rlatinc
      else
         write(*,*) 'Are you sure your map projection is right ?'
#ifdef MPP1
         call MPI_ABORT(MPI_COMM_WORLD,ierr)
#endif
         stop
      endif
      else if(iotyp.eq.2) then
#ifdef MPP1
         write(31,200) mjx-2,xlong_io(2,2),xlong_io(2,3)-xlong_io(2,2)
         write(31,201) ix-2,xlat_io(2,2),xlat_io(3,2)-xlat_io(2,2)
#else
         write(31,200) jx-2,xlong(2,2),xlong(2,3)-xlong(2,2)
         write(31,201) ix-2,xlat(2,2),xlat(3,2)-xlat(2,2)
#endif
      endif
 201  format('ydef ',I3,' linear ',f9.4,' ',f9.4)
      write(31,300) kx,((1013.25-ptop*10.)*a(k)+ptop*10.,k=kx,1,-1)
 300  format('zdef ',I2,' levels ',30f7.2)
      myear=ldatez/1000000
      month=(ldatez-myear*1000000)/10000
      nday =(ldatez-myear*1000000-month*10000)/100
      nhour=MOD(ldatez,100)
      CALL FINDDATE(nbase,ldatez)
      if(month.eq.12) then
         CALL FINDDATE(mnend,myear*1000000+1010100)
      else
         CALL FINDDATE(mnend,myear*1000000+month*10000+10100)
      endif
      CALL FINDDATE(jbend,IDATE2)
      if(ldatez.EQ.IDATE0) then
         number=(ibdyfrq/tapfrq+0.00001)*(min0(jbend,mnend)-nbase)+1
      else
         number=(ibdyfrq/tapfrq+0.00001)*(min0(jbend,mnend)-nbase)
      endif
      ifrq=chemfrq+0.00001
      if(ldatez.EQ.IDATE0) then
         write(31,400)number,nhour     ,cday(nday),cmonth(month)
     &               ,myear,ifrq
      else
         write(31,400)number,nhour+ifrq,cday(nday),cmonth(month)
     &               ,myear,ifrq
      endif
 400  format('tdef ',I4,' linear ',I2,'z',A2,A3,I4,' ',I2,'hr')
      if(iotyp.eq.2)  write(31,'(a)') 'theader 4'
      write(31,500) ntr +3 + 7*ntr +2
 500  format('vars ',I2)
 600  format(A8,' 0 99 ',A26)
 650  format(A4,I1,' ',I2,' 0 ',A26)
 651  format(A4,I2,' ',I2,' 0 ',A26)
 655  format(A8,I1,' 0 99 ',A26)
 656  format(A8,I2,' 0 99 ',A26)

      do itr=1,ntr
       if(itr.lt.10) then
         write(31,650) 'trac',itr,kx, 'tracer mix. rat  (Kg/Kg)'
       else
         write(31,651) 'trac',itr,kx, 'tracer mix. rat  (Kg/Kg)'
       endif
      end do  
      write(31,650) 'aext',8,kx, 'aer mix. ext. coef      '
      write(31,650) 'assa',8,kx, 'aer mix. sin. scat. alb '
      write(31,650) 'agfu',8,kx, 'aer mix. ass. par       '
      
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
      end do 
      write(31,600) 'acstoarf',' TOArad forcing av.(W/m2)'
      write(31,600) 'acstsrrf',' SRFrad forcing av.(W/m2)'

      write(31,700)
 700  format('endvars')
      close(31)
      return
      end
      subroutine GRADSCTL(CTLNAME)
      implicit none
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#     include "include/param2.cb"
#ifdef MPP1
#     include "include/io.cb"
      include 'mpiregcm.h'
      INTEGER ierr
#else
#     include "include/main.cb"
#endif
      real(kind=4)  clat,clon,dxsp,ptsp,plat,plon
      character proj*6
      COMMON /GRADS/clat,clon,dxsp,ptsp,plat,plon,proj
      real(kind=4)  TRUELATL,TRUELATH
      COMMON /SAVEPAR/ TRUELATL,TRUELATH
      character*12 CTLNAME
      logical there
      integer i,j,ny,nx
      real(kind=4)  alatmin,alatmax,alonmin,alonmax,rlatinc,rloninc
      real(kind=4)  centerj,centeri
      integer isystm,system
      external system
c
      inquire(file='output/'//CTLNAME,exist=there)
      if(there) isystm=system('/bin/rm output/'//CTLNAME)
      open(31,file='output/'//CTLNAME,status='new')
      write(31,10)
  10  format('dset ^OUT_HEAD')
      write(31,20)
  20  format('title RegCM domain information')
      if(ibintyp.eq.1) then
         write(31,30)
  30     format('options big_endian')
      elseif(ibintyp.eq.2) then
         write(31,40)
  40     format('options little_endian')
      endif
      write(31,50)
  50  format('undef -1.e34')
      if(proj.eq.'LAMCON'.or.proj.eq.'ROTMER') then
         alatmin= 999999.
         alatmax=-999999.
#ifdef MPP1
         do j=1,mjx
            if(xlat_io(1 ,j).lt.alatmin)alatmin=xlat_io(1 ,j)
            if(xlat_io(ix,j).gt.alatmax)alatmax=xlat_io(ix,j)
#else
         do j=1,jx
            if(xlat(1 ,j).lt.alatmin)alatmin=xlat(1 ,j)
            if(xlat(ix,j).gt.alatmax)alatmax=xlat(ix,j)
#endif
         enddo
         alonmin= 999999.
         alonmax=-999999.
         do i=1,ix
#ifdef MPP1
         do j=1,mjx
#else
         do j=1,jx
#endif
            if(clon.ge.0.0) then
#ifdef MPP1
               if(xlong_io(i,j).ge.0.0) then
                  alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
                  alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
#else
               if(xlong(i,j).ge.0.0) then
                  alonmin = amin1(alonmin,sngl(xlong(i,j)))
                  alonmax = amax1(alonmax,sngl(xlong(i,j)))
#endif
               else
#ifdef MPP1
                  if(abs(clon-xlong_io(i,j)).lt.
     &               abs(clon-(xlong_io(i,j)+360.))) then
                     alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
                     alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
#else
                  if(abs(clon-xlong(i,j)).lt.
     &               abs(clon-(xlong(i,j)+360.))) then
                     alonmin = amin1(alonmin,sngl(xlong(i,j)))
                     alonmax = amax1(alonmax,sngl(xlong(i,j)))
#endif
                  else
#ifdef MPP1
                     alonmin = amin1(alonmin,sngl(xlong_io(i,j))+360.)
                     alonmax = amax1(alonmax,sngl(xlong_io(i,j))+360.)
#else
                     alonmin = amin1(alonmin,sngl(xlong(i,j))+360.)
                     alonmax = amax1(alonmax,sngl(xlong(i,j))+360.)
#endif
                  endif
               endif
            else
#ifdef MPP1
               if(xlong_io(i,j).lt.0.0) then
                  alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
                  alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
#else
               if(xlong(i,j).lt.0.0) then
                  alonmin = amin1(alonmin,sngl(xlong(i,j)))
                  alonmax = amax1(alonmax,sngl(xlong(i,j)))
#endif
               else
#ifdef MPP1
                  if(abs(clon-xlong_io(i,j)).lt.
     &               abs(clon-(xlong_io(i,j)-360.))) then
                     alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
                     alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
#else
                  if(abs(clon-xlong(i,j)).lt.
     &               abs(clon-(xlong(i,j)-360.))) then
                     alonmin = amin1(alonmin,sngl(xlong(i,j)))
                     alonmax = amax1(alonmax,sngl(xlong(i,j)))
#endif
                  else
#ifdef MPP1
                     alonmin = amin1(alonmin,sngl(xlong_io(i,j))-360.)
                     alonmax = amax1(alonmax,sngl(xlong_io(i,j))-360.)
#else
                     alonmin = amin1(alonmin,sngl(xlong(i,j))-360.)
                     alonmax = amax1(alonmax,sngl(xlong(i,j))-360.)
#endif
                  endif
               endif
            endif
         enddo
         enddo
         rlatinc=dxsp*0.001/111./2.
         rloninc=dxsp*0.001/111./2.
         ny=2+nint(abs(alatmax-alatmin)/rlatinc)
         nx=1+nint(abs((alonmax-alonmin)/rloninc))

#ifdef MPP1
         centerj=(mjx-2)/2.
#else
         centerj=(jx-2)/2.
#endif
         centeri=(ix-2)/2.
      endif
      if(proj.eq.'LAMCON') then        ! Lambert projection
#ifdef MPP1
         write(31,100) mjx-2,ix-2,clat,clon,centerj,centeri,
     &                 truelatL,truelatH,clon,dxsp,dxsp
#else
         write(31,100) jx-2,ix-2,clat,clon,centerj,centeri,
     &                 truelatL,truelatH,clon,dxsp,dxsp
#endif
 100  format('pdef ',i4,1x,i4,1x,'lcc',7(1x,f7.2),1x,2(f7.0,1x))
         write(31,110) nx+2,alonmin-rloninc,rloninc
 110  format('xdef ',i4,' linear ',f7.2,1x,f7.4)
         write(31,120) ny+2,alatmin-rlatinc,rlatinc
 120  format('ydef ',i4,' linear ',f7.2,1x,f7.4)
      elseif(proj.eq.'POLSTR') then    !
      elseif(proj.eq.'NORMER') then
#ifdef MPP1
         write(31,200)  mjx-2,xlong_io(2,2),xlong_io(2,3)-xlong_io(2,2)
#else
         write(31,200)  jx-2,xlong(2,2),xlong(2,3)-xlong(2,2)
#endif
 200  format('xdef ',I3,' linear ',f9.4,' ',f9.4)
         write(31,210) ix-2
 210  format('ydef ',I3,' levels')
#ifdef MPP1
         write(31,220) (xlat_io(i,2),i=2,ix-1)
#else
         write(31,220) (xlat(i,2),i=2,ix-1)
#endif
 220  format(10f7.2)
      elseif(proj.eq.'ROTMER') then
         write(*,*) 'Note that rotated Mercartor (ROTMER)'
     &             ,' projections are not supported by GrADS.'
         write(*,*) '  Although not exact, the eta.u projection'
     &             ,' in GrADS is somewhat similar.'
         write(*,*) ' FERRET, however, does support this projection.'
#ifdef MPP1
         write(31,230) mjx-2,ix-2,plon,plat,dxsp/111000.
     &                                ,dxsp/111000.*.95238
#else
         write(31,230) jx-2,ix-2,plon,plat,dxsp/111000.
     &                                ,dxsp/111000.*.95238
#endif
 230  format('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
         write(31,110) nx+2,alonmin-rloninc,rloninc
         write(31,120) ny+2,alatmin-rlatinc,rlatinc
      else
         write(*,*) 'Are you sure your map projection is right ?'
#ifdef MPP1
         call MPI_ABORT(MPI_COMM_WORLD,ierr)
#endif
         stop
      endif
      write(31,300) 1,1000.
 300  format('zdef ',I1,' levels ',f7.2)
      write(31,400) 1
 400  format('tdef ',I1,' linear 00z01Jan2001 1mo')
      write(31,500) 11
 500  format('vars ',I2)
 600  format(A8,'0 99 ',A26)
      write(31,600) 'head    ','header information         '
      write(31,600) 'ht      ','surface elevation          '
      write(31,600) 'htsd    ','surface elevation std dev  '
      write(31,600) 'veg2d   ','vegetation type in BATS    '
      write(31,600) 'landuse ','surface landuse type       '
      write(31,600) 'xlat    ','latitude  of cross points  '
      write(31,600) 'xlong   ','longitude of cross points  '
      write(31,600) 'xmap    ','map factors of cross points'
      write(31,600) 'dmap    ','map factors of dot points  '
      write(31,600) 'coriol  ','coriol force               '
      write(31,600) 'mask    ','land/sea mask              '
      write(31,700)
 700  format('endvars')
      close(31)
      return
      end
      subroutine GRADSOUT(CTLNAME)
      implicit none
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#     include "include/date.cb"
#     include "include/param1.cb"
#     include "include/param2.cb"
#     include "include/param3.cb"
#ifdef MPP1
#     include "include/io.cb"
      include 'mpiregcm.h'
      INTEGER ierr
#else
#     include "include/main.cb"
#endif
      real(kind=4)  clat,clon,dxsp,ptsp,plat,plon
      character proj*6
      COMMON /GRADS/clat,clon,dxsp,ptsp,plat,plon,proj
      real(kind=4)  TRUELATL,TRUELATH
      COMMON /SAVEPAR/ TRUELATL,TRUELATH
      character*18 CTLNAME
      logical there
      integer i,j,ny,nx
      real(kind=4)  alatmin,alatmax,alonmin,alonmax,rlatinc,rloninc
      real(kind=4)  centerj,centeri
      integer myear,month,nday,nhour,k
      integer number,nbase,mnend,jbend,ifrq
      character*2 cday(31)
      data cday/'01','02','03','04','05','06','07','08','09','10',
     &          '11','12','13','14','15','16','17','18','19','20',
     &     '21','22','23','24','25','26','27','28','29','30','31'/
      character*3 cmonth(12)
      data cmonth/'jan','feb','mar','apr','may','jun',
     &            'jul','aug','sep','oct','nov','dec'/
      integer isystm,system
      external system
c
      inquire(file='output/'//CTLNAME,exist=there)
      if(there) isystm=system('/bin/rm output/'//CTLNAME)
      open(31,file='output/'//CTLNAME,status='new')
      write(31,10) CTLNAME(1:14)
  10  format('dset ^',A14)
      write(31,20)
  20  format('title RegCM normal output variables')
      if(ibintyp.eq.1) then
         write(31,30)
  30     format('options big_endian')
      elseif(ibintyp.eq.2) then
         write(31,40)
  40     format('options little_endian')
      endif
      if(iotyp.eq.2) write(31,'(a)') 'options sequential'
      write(31,50)
  50  format('undef -1.e34')
      if(proj.eq.'LAMCON'.or.proj.eq.'ROTMER') then
         alatmin= 999999.
         alatmax=-999999.
#ifdef MPP1
         do j=1,mjx
            if(xlat_io(1 ,j).lt.alatmin)alatmin=xlat_io(1 ,j)
            if(xlat_io(ix,j).gt.alatmax)alatmax=xlat_io(ix,j)
#else
         do j=1,jx
            if(xlat(1 ,j).lt.alatmin)alatmin=xlat(1 ,j)
            if(xlat(ix,j).gt.alatmax)alatmax=xlat(ix,j)
#endif
         enddo
         alonmin= 999999.
         alonmax=-999999.
         do i=1,ix
#ifdef MPP1
         do j=1,mjx
#else
         do j=1,jx
#endif
            if(clon.ge.0.0) then
#ifdef MPP1
               if(xlong_io(i,j).ge.0.0) then
                  alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
                  alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
#else
               if(xlong(i,j).ge.0.0) then
                  alonmin = amin1(alonmin,sngl(xlong(i,j)))
                  alonmax = amax1(alonmax,sngl(xlong(i,j)))
#endif
               else
#ifdef MPP1
                  if(abs(clon-xlong_io(i,j)).lt.
     &               abs(clon-(xlong_io(i,j)+360.))) then
                     alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
                     alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
#else
                  if(abs(clon-xlong(i,j)).lt.
     &               abs(clon-(xlong(i,j)+360.))) then
                     alonmin = amin1(alonmin,sngl(xlong(i,j)))
                     alonmax = amax1(alonmax,sngl(xlong(i,j)))
#endif
                  else
#ifdef MPP1
                     alonmin = amin1(alonmin,sngl(xlong_io(i,j))+360.)
                     alonmax = amax1(alonmax,sngl(xlong_io(i,j))+360.)
#else
                     alonmin = amin1(alonmin,sngl(xlong(i,j))+360.)
                     alonmax = amax1(alonmax,sngl(xlong(i,j))+360.)
#endif
                  endif
               endif
            else
#ifdef MPP1
               if(xlong_io(i,j).lt.0.0) then
                  alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
                  alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
#else
               if(xlong(i,j).lt.0.0) then
                  alonmin = amin1(alonmin,sngl(xlong(i,j)))
                  alonmax = amax1(alonmax,sngl(xlong(i,j)))
#endif
               else
#ifdef MPP1
                  if(abs(clon-xlong_io(i,j)).lt.
     &               abs(clon-(xlong_io(i,j)-360.))) then
                     alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
                     alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
#else
                  if(abs(clon-xlong(i,j)).lt.
     &               abs(clon-(xlong(i,j)-360.))) then
                     alonmin = amin1(alonmin,sngl(xlong(i,j)))
                     alonmax = amax1(alonmax,sngl(xlong(i,j)))
#endif
                  else
#ifdef MPP1
                     alonmin = amin1(alonmin,sngl(xlong_io(i,j))-360.)
                     alonmax = amax1(alonmax,sngl(xlong_io(i,j))-360.)
#else
                     alonmin = amin1(alonmin,sngl(xlong(i,j))-360.)
                     alonmax = amax1(alonmax,sngl(xlong(i,j))-360.)
#endif
                  endif
               endif
            endif
         enddo
         enddo
         rlatinc=dxsp*0.001/111./2.
         rloninc=dxsp*0.001/111./2.
         ny=2+nint(abs(alatmax-alatmin)/rlatinc)
         nx=1+nint(abs((alonmax-alonmin)/rloninc))

#ifdef MPP1
         centerj=(mjx-2)/2.
#else
         centerj=(jx-2)/2.
#endif
         centeri=(ix-2)/2.
      endif
      if(iotyp.eq.1) then
      if(proj.eq.'LAMCON') then        ! Lambert projection
#ifdef MPP1
         write(31,100) mjx-2,ix-2,clat,clon,centerj,centeri,
     &                 truelatL,truelatH,clon,dxsp,dxsp
#else
         write(31,100) jx-2,ix-2,clat,clon,centerj,centeri,
     &                 truelatL,truelatH,clon,dxsp,dxsp
#endif
 100  format('pdef ',i4,1x,i4,1x,'lccr',7(1x,f7.2),1x,2(f7.0,1x))
         write(31,110) nx+2,alonmin-rloninc,rloninc
 110  format('xdef ',i4,' linear ',f7.2,1x,f7.4)
         write(31,120) ny+2,alatmin-rlatinc,rlatinc
 120  format('ydef ',i4,' linear ',f7.2,1x,f7.4)
      elseif(proj.eq.'POLSTR') then    !
      elseif(proj.eq.'NORMER') then
#ifdef MPP1
         write(31,200)  mjx-2,xlong_io(2,2),xlong_io(2,3)-xlong_io(2,2)
#else
         write(31,200)  jx-2,xlong(2,2),xlong(2,3)-xlong(2,2)
#endif
 200  format('xdef ',I3,' linear ',f9.4,' ',f9.4)
         write(31,210) ix-2
 210  format('ydef ',I3,' levels')
#ifdef MPP1
         write(31,220) (xlat_io(i,2),i=2,ix-1)
#else
         write(31,220) (xlat(i,2),i=2,ix-1)
#endif
 220  format(10f7.2)
      elseif(proj.eq.'ROTMER') then
         write(*,*) 'Note that rotated Mercartor (ROTMER)'
     &             ,' projections are not supported by GrADS.'
         write(*,*) '  Although not exact, the eta.u projection'
     &             ,' in GrADS is somewhat similar.'
         write(*,*) ' FERRET, however, does support this projection.'
#ifdef MPP1
         write(31,230) mjx-2,ix-2,plon,plat,dxsp/111000.
     &                                ,dxsp/111000.*.95238
#else
         write(31,230) jx-2,ix-2,plon,plat,dxsp/111000.
     &                                ,dxsp/111000.*.95238
#endif
 230  format('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
         write(31,110) nx+2,alonmin-rloninc,rloninc
         write(31,120) ny+2,alatmin-rlatinc,rlatinc
      else
         write(*,*) 'Are you sure your map projection is right ?'
#ifdef MPP1
         call MPI_ABORT(MPI_COMM_WORLD,ierr)
#endif
         stop
      endif
      else if(iotyp.eq.2) then
#ifdef MPP1
         write(31,200) mjx-2,xlong_io(2,2),xlong_io(2,3)-xlong_io(2,2)
         write(31,201) ix-2,xlat_io(2,2),xlat_io(3,2)-xlat_io(2,2)
#else
         write(31,200) jx-2,xlong(2,2),xlong(2,3)-xlong(2,2)
         write(31,201) ix-2,xlat(2,2),xlat(3,2)-xlat(2,2)
#endif
      endif
 201  format('ydef ',I3,' linear ',f9.4,' ',f9.4)
      write(31,300) kx,((1013.25-ptop*10.)*a(k)+ptop*10.,k=kx,1,-1)
 300  format('zdef ',I2,' levels ',30f7.2)
      myear=ldatez/1000000
      month=(ldatez-myear*1000000)/10000
      nday =(ldatez-myear*1000000-month*10000)/100
      nhour=MOD(ldatez,100)
      CALL FINDDATE(nbase,ldatez)
      if(month.eq.12) then
         CALL FINDDATE(mnend,myear*1000000+1010100)
      else
         CALL FINDDATE(mnend,myear*1000000+month*10000+10100)
      endif
      CALL FINDDATE(jbend,IDATE2)
      if(ldatez.EQ.IDATE0) then
         number=(ibdyfrq/tapfrq+0.00001)*(min0(jbend,mnend)-nbase)+1
      else
         number=(ibdyfrq/tapfrq+0.00001)*(min0(jbend,mnend)-nbase)
      endif
      ifrq=tapfrq+0.00001
      if(ldatez.EQ.IDATE0) then
         write(31,400)number,nhour     ,cday(nday),cmonth(month)
     &               ,myear,ifrq
      else
         write(31,400)number,nhour+ifrq,cday(nday),cmonth(month)
     &               ,myear,ifrq
      endif
 400  format('tdef ',I4,' linear ',I2,'z',A2,A3,I4,' ',I2,'hr')
      if(iotyp.eq.2) write(31,'(a)') 'theader 4'
      write(31,500) 10+1
 500  format('vars ',I2)
 600  format(A8,'0 99 ',A26)
 650  format(A8,I2,' 0 ',A26)
 651  format(A8,I2,' 33,100 ',A36)
 652  format(A8,I2,' 34,100 ',A36)
      if(proj.eq.'LAMCON'.and.iotyp.eq.1) then        ! Lambert projection
         write(31,651) 'u       ',kx,'westerly wind (m/s)        '
         write(31,652) 'v       ',kx,'southerly wind (m/s)       '
      else
         write(31,650) 'u       ',kx,'westerly wind (m/s)        '
         write(31,650) 'v       ',kx,'southerly wind (m/s)       '
      endif
      write(31,650) 'w       ',kx,'omega (hPa/s)   p-velocity '
      write(31,650) 't       ',kx,'air temperature (degree)   '
      write(31,650) 'qv      ',kx,'water vapor mixing ratio   '
      write(31,650) 'qc      ',kx,'cloud water mixing ratio   '
      write(31,600) 'psa     ',   'surface pressure (hPa)     '
      write(31,600) 'tpr     ',   'total precipitation(mm/day)'
      write(31,600) 'tgb     ',   'lower groud temp. in BATS  '
      write(31,600) 'swt     ',   'total soil water in mm H2O '
      write(31,600) 'rno     ',   'accumulated infiltration   '
      write(31,700)
 700  format('endvars')
      close(31)
      return
      end
      subroutine GRADSRAD(CTLNAME)
      implicit none
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#     include "include/date.cb"
#     include "include/param1.cb"
#     include "include/param2.cb"
#     include "include/param3.cb"
#ifdef MPP1
#     include "include/io.cb"
      include 'mpiregcm.h'
      INTEGER ierr
#else
#     include "include/main.cb"
#endif
      real(kind=4)  clat,clon,dxsp,ptsp,plat,plon
      character proj*6
      COMMON /GRADS/clat,clon,dxsp,ptsp,plat,plon,proj
      real(kind=4)  TRUELATL,TRUELATH
      COMMON /SAVEPAR/ TRUELATL,TRUELATH
      character*18 CTLNAME
      logical there
      integer i,j,ny,nx
      real(kind=4)  alatmin,alatmax,alonmin,alonmax,rlatinc,rloninc
      real(kind=4)  centerj,centeri
      integer myear,month,nday,nhour,k
      integer number,nbase,mnend,jbend,ifrq
      character*2 cday(31)
      data cday/'01','02','03','04','05','06','07','08','09','10',
     &          '11','12','13','14','15','16','17','18','19','20',
     &     '21','22','23','24','25','26','27','28','29','30','31'/
      character*3 cmonth(12)
      data cmonth/'jan','feb','mar','apr','may','jun',
     &            'jul','aug','sep','oct','nov','dec'/
      integer isystm,system
      external system
c
      inquire(file='output/'//CTLNAME,exist=there)
      if(there) isystm=system('/bin/rm output/'//CTLNAME)
      open(31,file='output/'//CTLNAME,status='new')
      write(31,10) CTLNAME(1:14)
  10  format('dset ^',A14)
      write(31,20)
  20  format('title RegCM normal output variables')
      if (ibintyp.eq.1) then
         write(31,30)
  30     format('options big_endian')
      else if (ibintyp.eq.2) then
         write(31,40)
  40     format('options little_endian')
      endif
      if(iotyp.eq.2) write(31,'(a)') 'options sequential'
      write(31,50)
  50  format('undef -1.34')
      if(proj.eq.'LAMCON'.or.proj.eq.'ROTMER') then
         alatmin= 999999.
         alatmax=-999999.
#ifdef MPP1
         do j=1,mjx
            if(xlat_io(1 ,j).lt.alatmin)alatmin=xlat_io(1 ,j)
            if(xlat_io(ix,j).gt.alatmax)alatmax=xlat_io(ix,j)
#else
         do j=1,jx
            if(xlat(1 ,j).lt.alatmin)alatmin=xlat(1 ,j)
            if(xlat(ix,j).gt.alatmax)alatmax=xlat(ix,j)
#endif
         enddo
         alonmin= 999999.
         alonmax=-999999.
         do i=1,ix
#ifdef MPP1
         do j=1,mjx
#else
         do j=1,jx
#endif
            if(clon.ge.0.0) then
#ifdef MPP1
               if(xlong_io(i,j).ge.0.0) then
                  alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
                  alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
#else
               if(xlong(i,j).ge.0.0) then
                  alonmin = amin1(alonmin,sngl(xlong(i,j)))
                  alonmax = amax1(alonmax,sngl(xlong(i,j)))
#endif
               else
#ifdef MPP1
                  if(abs(clon-xlong_io(i,j)).lt.
     &               abs(clon-(xlong_io(i,j)+360.))) then
                     alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
                     alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
#else
                  if(abs(clon-xlong(i,j)).lt.
     &               abs(clon-(xlong(i,j)+360.))) then
                     alonmin = amin1(alonmin,sngl(xlong(i,j)))
                     alonmax = amax1(alonmax,sngl(xlong(i,j)))
#endif
                  else
#ifdef MPP1
                     alonmin = amin1(alonmin,sngl(xlong_io(i,j))+360.)
                     alonmax = amax1(alonmax,sngl(xlong_io(i,j))+360.)
#else
                     alonmin = amin1(alonmin,sngl(xlong(i,j))+360.)
                     alonmax = amax1(alonmax,sngl(xlong(i,j))+360.)
#endif
                  endif
               endif
            else
#ifdef MPP1
               if(xlong_io(i,j).lt.0.0) then
                  alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
                  alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
#else
               if(xlong(i,j).lt.0.0) then
                  alonmin = amin1(alonmin,sngl(xlong(i,j)))
                  alonmax = amax1(alonmax,sngl(xlong(i,j)))
#endif
               else
#ifdef MPP1
                  if(abs(clon-xlong_io(i,j)).lt.
     &               abs(clon-(xlong_io(i,j)-360.))) then
                     alonmin = amin1(alonmin,sngl(xlong_io(i,j)))
                     alonmax = amax1(alonmax,sngl(xlong_io(i,j)))
#else
                  if(abs(clon-xlong(i,j)).lt.
     &               abs(clon-(xlong(i,j)-360.))) then
                     alonmin = amin1(alonmin,sngl(xlong(i,j)))
                     alonmax = amax1(alonmax,sngl(xlong(i,j)))
#endif
                  else
#ifdef MPP1
                     alonmin = amin1(alonmin,sngl(xlong_io(i,j))-360.)
                     alonmax = amax1(alonmax,sngl(xlong_io(i,j))-360.)
#else
                     alonmin = amin1(alonmin,sngl(xlong(i,j))-360.)
                     alonmax = amax1(alonmax,sngl(xlong(i,j))-360.)
#endif
                  endif
               endif
            endif
         enddo
         enddo
         rlatinc=dxsp*0.001/111./2.
         rloninc=dxsp*0.001/111./2.
         ny=2+nint(abs(alatmax-alatmin)/rlatinc)
         nx=1+nint(abs((alonmax-alonmin)/rloninc))

#ifdef MPP1
         centerj=(mjx-2)/2.
#else
         centerj=(jx-2)/2.
#endif
         centeri=(ix-2)/2.
      endif
      if(iotyp.eq.1) then
      if(proj.eq.'LAMCON') then        ! Lambert projection
#ifdef MPP1
         write(31,100) mjx-2,ix-2,clat,clon,centerj,centeri,
     &                 truelatL,truelatH,clon,dxsp,dxsp
#else
         write(31,100) jx-2,ix-2,clat,clon,centerj,centeri,
     &                 truelatL,truelatH,clon,dxsp,dxsp
#endif
 100  format('pdef ',i4,1x,i4,1x,'lcc',7(1x,f7.2),1x,2(f7.0,1x))
         write(31,110) nx+2,alonmin-rloninc,rloninc
 110  format('xdef ',i4,' linear ',f7.2,1x,f7.4)
         write(31,120) ny+2,alatmin-rlatinc,rlatinc
 120  format('ydef ',i4,' linear ',f7.2,1x,f7.4)
      elseif(proj.eq.'POLSTR') then    !
      elseif(proj.eq.'NORMER') then
#ifdef MPP1
         write(31,200)  mjx-2,xlong_io(2,2),xlong_io(2,3)-xlong_io(2,2)
#else
         write(31,200)  jx-2,xlong(2,2),xlong(2,3)-xlong(2,2)
#endif
 200  format('xdef ',I3,' linear ',f9.4,' ',f9.4)
         write(31,210) ix-2
 210  format('ydef ',I3,' levels')
#ifdef MPP1
         write(31,220) (xlat_io(i,2),i=2,ix-1)
#else
         write(31,220) (xlat(i,2),i=2,ix-1)
#endif
 220  format(10f7.2)
      elseif(proj.eq.'ROTMER') then
         write(*,*) 'Note that rotated Mercartor (ROTMER)'
     &             ,' projections are not supported by GrADS.'
         write(*,*) '  Although not exact, the eta.u projection'
     &             ,' in GrADS is somewhat similar.'
         write(*,*) ' FERRET, however, does support this projection.'
#ifdef MPP1
         write(31,230) mjx-2,ix-2,plon,plat,dxsp/111000.
     &                                ,dxsp/111000.*.95238
#else
         write(31,230) jx-2,ix-2,plon,plat,dxsp/111000.
     &                                ,dxsp/111000.*.95238
#endif
 230  format('pdef ',i4,1x,i4,1x,'eta.u',2(1x,f7.3),2(1x,f9.5))
         write(31,110) nx+2,alonmin-rloninc,rloninc
         write(31,120) ny+2,alatmin-rlatinc,rlatinc
      else
         write(*,*) 'Are you sure your map projection is right ?'
#ifdef MPP1
         call MPI_ABORT(MPI_COMM_WORLD,ierr)
#endif
         stop
      endif
      else if(iotyp.eq.2) then
#ifdef MPP1
         write(31,200) mjx-2,xlong_io(2,2),xlong_io(2,3)-xlong_io(2,2)
         write(31,201) ix-2,xlat_io(2,2),xlat_io(3,2)-xlat_io(2,2)
#else
         write(31,200) jx-2,xlong(2,2),xlong(2,3)-xlong(2,2)
         write(31,201) ix-2,xlat(2,2),xlat(3,2)-xlat(2,2)
#endif
      endif
 201  format('ydef ',I3,' linear ',f9.4,' ',f9.4)
      write(31,300) kx,((1013.25-ptop*10.)*a(k)+ptop*10.,k=kx,1,-1)
 300  format('zdef ',I2,' levels ',30f7.2)
      myear=ldatez/1000000
      month=(ldatez-myear*1000000)/10000
      nday =(ldatez-myear*1000000-month*10000)/100
      nhour=MOD(ldatez,100)
      CALL FINDDATE(nbase,ldatez)
      if(month.eq.12) then
         CALL FINDDATE(mnend,myear*1000000+1010100)
      else
         CALL FINDDATE(mnend,myear*1000000+month*10000+10100)
      endif
      CALL FINDDATE(jbend,IDATE2)
      if(ldatez.EQ.IDATE0) then
         number=(ibdyfrq/radisp+0.00001)*(min0(jbend,mnend)-nbase)+1
      else
         number=(ibdyfrq/radisp+0.00001)*(min0(jbend,mnend)-nbase)
      endif
      ifrq=radisp+0.00001
      if(ldatez.EQ.IDATE0) then
         write(31,400)number,nhour     ,cday(nday),cmonth(month)
     &               ,myear,ifrq
      else
         write(31,400)number,nhour+ifrq,cday(nday),cmonth(month)
     &               ,myear,ifrq
      endif
 400  format('tdef ',I4,' linear ',I2,'z',A2,A3,I4,' ',I2,'hr')
      if(iotyp.eq.2) write(31,'(a)') 'theader 4'
      write(31,500) 13
 500  format('vars ',I2)
 600  format(A6,'0 99 ',A36)
 650  format(A6,I2,' 0 ',A36)
      write(31,650)'cld   ',kx,'cloud fractional cover               '
      write(31,650)'clwp  ',kx,'cloud liquid water path              '
      write(31,650)'qrs   ',kx,'solar heating rate                   '
      write(31,650)'qrl   ',kx,'longwave cooling rate                '
      write(31,600)'frsa  ',   'surface absorbed solar flux          '
      write(31,600)'frla  ',   'longwave cooling of surface          '
      write(31,600)'clrst ',   'clearsky total column abs solar flux '
      write(31,600)'clrss ',   'clearsky surface absorbed solar flux '
      write(31,600)'clrlt ',   'clearsky net upward LW flux at TOA   '
      write(31,600)'clrls ',   'clearsky LW cooling at surface (W/m2)'
      write(31,600)'solin ',   'instantaneous incident solar (W/m2)  '
      write(31,600)'sabtp ',   'total column absorbed solar flux W/m2'
      write(31,600)'firtp ',   'net upward LW flux at TOA (W/m2)     '
      write(31,700)
 700  format('endvars')
      close(31)
      return
      end

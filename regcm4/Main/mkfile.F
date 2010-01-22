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

      subroutine mkfile

      implicit none
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#     include "include/param1.cb"
#     include "include/param2.cb"
#     include "include/param3.cb"
#     include "include/iunits.cb"
#     include "include/date.cb"
#     include "include/bats.cb2"
#     include "include/trachem.cb"

      character itype*3, filatm*14, filsrf*14, filsub*14, filrad*14
     &        , fillak*14, filchem *14           
     &        , a78*78, a16*16, a22*22, a1*1, a2*2, a120*120
     &        , filnam1*40, filnam2*40
      integer n,nmo,nline,nfmax,idatepp,iym1,iym2
      parameter (nfmax=999)
      integer iyr(nfmax),imo(nfmax),idate1d(nfmax)
      real(kind=4)  dto,dtb,dtr,dtc
      logical there

      integer nrcout
      COMMON /GENREC/ nrcout
      integer nrcrad
      COMMON /RADREC/ nrcrad

Cchem2
      integer nrcchem
      COMMON /CHEREC/ nrcchem
Cchem2_ 
      integer isystm,system
      external system

      print*,' '
      print*,'******* OPENING NEW OUTPUT FILES:',idatex
      if (iftape) then
        close(iutdat)
        itype = 'ATM'
        write(filatm,199) itype,idatex
        if (iotyp.eq.1) then
#ifdef MPP1
          open(iutdat,file='output/'//filatm,status='unknown'
     &        ,form='unformatted',recl=(ix-2)*(mjx-2)*ibyte
     &        ,access='direct')
#else
          open(iutdat,file='output/'//filatm,status='unknown'
     &        ,form='unformatted',recl=(ix-2)*(jx-2)*ibyte
     &        ,access='direct')
#endif
          nrcout=0
        else if (iotyp.eq.2) then
          open(iutdat,file='output/'//filatm,status='unknown'
     &        ,form='unformatted')
        end if
        call GRADSOUT(filatm//'.ctl')
        print*,'OPENING NEW OUT FILE: output/',filatm
      end if

      if (ifbat) then
        close(iutbat)
        itype = 'SRF'
        write(filsrf,199) itype,idatex
        if (iotyp.eq.1) then
#ifdef MPP1
          open(iutbat,file='output/'//filsrf,status='unknown'
     &        ,form='unformatted',recl=(ix-2)*(mjx-2)*ibyte
     &        ,access='direct')
#else
          open(iutbat,file='output/'//filsrf,status='unknown'
     &        ,form='unformatted',recl=(ix-2)*(jx-2)*ibyte
     &        ,access='direct')
#endif
          nrcbat=0
        else if (iotyp.eq.2) then
          open(iutbat,file='output/'//filsrf,status='unknown'
     &        ,form='unformatted')
        end if
        call GRADSBAT(filsrf//'.ctl')
        print*,'OPENING NEW BAT FILE: output/',filsrf
      end if

      if (NSG.gt.1.and.ifsub) then
        close(iutsub)
        itype = 'SUB'
        write(filsub,199) itype,idatex
        if (iotyp.eq.1) then
#ifdef MPP1
          open(iutsub,file='output/'//filsub,status='unknown'
     &        ,form='unformatted',recl=(ix-2)*(mjx-2)*NNSG*ibyte
     &        ,access='direct')
#else
          open(iutsub,file='output/'//filsub,status='unknown'
     &        ,form='unformatted',recl=(ix-2)*(jx-2)*NNSG*ibyte
     &        ,access='direct')
#endif
          nrcsub=0
        else if (iotyp.eq.2) then
          open(iutsub,file='output/'//filsub,status='unknown'
     &        ,form='unformatted')
        end if
        call GRADSSUB(filsub//'.ctl')
        print*,'OPENING NEW SUB FILE: output/',filsub
      end if

      if (ifrad) then
        close(iutrad)
        itype = 'RAD'
        write(filrad,199) itype,idatex
        if (iotyp.eq.1) then
#ifdef MPP1
          open(iutrad,file='output/'//filrad,status='unknown'
     &        ,form='unformatted',recl=(ix-2)*(mjx-2)*ibyte
     &        ,access='direct')
#else
          open(iutrad,file='output/'//filrad,status='unknown'
     &        ,form='unformatted',recl=(ix-2)*(jx-2)*ibyte
     &        ,access='direct')
#endif
          nrcrad=0
        else if (iotyp.eq.2) then
          open(iutrad,file='output/'//filrad,status='unknown'
     &        ,form='unformatted')
        end if
        call GRADSRAD(filrad//'.ctl')
        print*,'OPENING NEW RAD FILE: output/',filrad
      end if

Cchem2 
      if (ichem .eq. 1) then
        if (ifchem) then
          close(iutchem)
          itype = 'CHE'
          write(filchem,199) itype,idatex
          if (iotyp.eq.1) then
#ifdef MPP1
            open(iutchem,file='output/'//filchem,status='unknown'
     &        ,form='unformatted',recl=(ix-2)*(mjx-2)*ibyte
     &        ,access='direct')
#else
            open(iutchem,file='output/'//filchem,status='unknown'
     &        ,form='unformatted',recl=(ix-2)*(jx-2)*ibyte
     &        ,access='direct')
#endif
            nrcchem=0
          else if (iotyp.eq.2) then
            open(iutchem,file='output/'//filchem,status='unknown'
     &           ,form='unformatted')
          end if
          call GRADSCHEM(filchem//'.ctl')
          print*,'OPENING NEW CHEM FILE: output/',filchem
        end if
      end if
Cchem2_
      if(lakemod.eq.1)then
        close(iutlak)
        itype = 'LAK'
        write(fillak,199) itype,idatex
        open(iutlak,file='output/'//fillak,status='unknown'
     &      ,form='unformatted')
        print*,'OPENING NEW LAK FILE: ',fillak
      end if

 199  format(a3,'.',i10)

      if (jyear.eq.jyear0.and.ktau.eq.0) then
        inquire(file='postproc.in',exist=there)
        if (there) then
          isystm=system('/bin/mv -f postproc.in postproc.in.bak')
        end if
        nmo = (idate2/1000000-idate0/1000000)*12
     &      + (idate2/10000-(idate2/1000000)*100)
     &      - (idate0/10000-(idate0/1000000)*100)
        nmo=min(nmo,nfmax)
        idatepp = idate0
        iyr(1) = idatepp/1000000
        imo(1) = (idatepp - iyr(1)*1000000)/10000
        idate1d(1) = idatepp
        do n=2,nmo
          idatepp = (idatepp/10000)*10000 + 10100
          iyr(n) = idatepp/1000000
          imo(n) = (idatepp - iyr(n)*1000000)/10000
          if (imo(n).gt.12) then
            iyr(n) = iyr(n) + 1
            imo(n) = 1
            idatepp = iyr(n)*1000000 + imo(n)*10000 + 100
          end if
          idate1d(n) = idatepp
        enddo
        if(nmo.le.1) then
           nmo=1
           idate1d(2) = idate0
        end if

c **** Write out postprocessing information
        OPEN(99,FILE='postproc.in',FORM='FORMATTED',STATUS='unknown')
        a1=char(39)
        write(99,101) idate1d(2)
 101    format(i10,',     ! idate0 = First date in File (yymmddhh)')
        write(99,102) idate1d(2)
 102    format(i10
     &   ,',     ! idate1 = Start date for averaging and re-writing')
        write(99,103) idate2
 103    format(i10
     &        ,',     ! idate2 = End date for averaging and re-writing')
        a78 = '2,              ! iotype: 1=I*2 NetCDF;'//
     &        ' 2=r*4 NetCDF; 3=grads; 4=vis5d'
        write(99,105) a78
        a78 = '.false.,        ! Write out header?'
        write(99,105) a78
        a78 = '.false.,        ! Write out all RegCM data b/twn'//
     &        ' idate1 & idate2?'
        write(99,105) a78
        a78 = '.false.,        ! Average RegCM data b/twn idate1'//
     &        ' & idate2?'
        write(99,105) a78
        a78 = '.false.,        ! Diurnali avg of RegCM data b/twn'//
     &        ' idate1 & idate2?'
        write(99,105) a78
        a78 = '.true.,         ! Continually average ATM data b/twn'//
     &        ' idate1 & idate2?'
        write(99,105) a78
        a78='-1.,            ! No. Days for continual averaging'//
     &       ' (-1=monthly;1=daily;5=5day)'
        write(99,105)a78
        if (nmo.le.1) then
          iym1 = iyr(nmo)*100+imo(nmo)
          write(a16,196) a1,iym1,a1
        else if (nmo.le.2) then
          iym1 = iyr(nmo)*100+imo(nmo)
          write(a16,196) a1,iym1,a1
        else
          iym1 = iyr(2)*100+imo(2)
          iym2 = iyr(nmo)*100+imo(nmo)
          write(a16,197) a1,iym1,iym2,a1
        end if
 196    format(a1,i6,a1,',')
 197    format(a1,i6,'_',i6,a1,',')
        a78 = a16//'! postproc output filename (not'//
     &        ' including type & ext)'
        write(99,105) a78

        a78=a1//'../Input'//a1//',     ! ICBC directory'
        write(99,105) a78
        a78=a1//'output'//a1//',       ! RegCM Output directory'
        write(99,105) a78
        a78=a1//'DOMAIN.INFO'//a1//
     &          ',  ! Domain Info Filename (from terrain)'
        write(99,105)a78
        a78=a1//'OUT_HEAD'//a1//',     ! Header File name (from RegCM)'
        write(99,105)a78
        if (nmo.le.1) then
          a2 = 'st'
          write(a22,198) a1,idate1d(1),a1,nmo,a2
          a78 = a22//' RegCM Output File Extension'
          write(99,105) a78
        else
          do n=2,nmo
            if (n.eq.2) then
              a2 = 'st'
            else if (n.eq.3) then
              a2 = 'nd'
            else if (n.eq.4) then
              a2 = 'rd'
            else
              a2 = 'th'
            end if
            write(a22,198) a1,idate1d(n),a1,n-1,a2
            a78 = a22//'RegCM Output File Extension'
            write(99,105) a78
          end do
        end if
 198    format(a1,i10,a1,',   !',i3,a2)

 105    format(a78)
        close(98)
        close(99)

        inquire(file='../PostProc',exist=there)
        if(.not.there) then
          isystm=system('mkdir ../PostProc')
        endif
        filnam1 = '../PostProc/postproc.param'
        filnam2 = '../PostProc/postproc.param.bak'
        inquire(file=filnam1,exist=there)
        if (there) then
          a120 = '/bin/mv -f '//filnam1//' '//filnam2
          isystm=system(a120)
        end if
        OPEN(99,FILE=filnam1,FORM='FORMATTED'
     &      ,STATUS='unknown')
        nline = 0
        dto = tapfrq
        dtb = batfrq
        dtr = radisp
        dtc = chemfrq
        a78 = 'cccc SET DOMAIN DIMENSIONS'
        write(99,105)a78
        a78 = 'cccc ny = number of north-south points'
        write(99,105)a78
        a78 = 'cccc nx = number of east-west points'
        write(99,105)a78
        a78 = 'cccc nz = number of vertical levels'
        write(99,105)a78
        a78 = 'cccc SET MODEL OUTPUT INTERVALS'
        write(99,105)a78
        a78 = 'cccc dtbc = ibdyfrq = ICBC output interval (hrs)'
        write(99,105)a78
        a78 = 'cccc dtout = atmfrq = atmospheric output interval (hrs)'
        write(99,105)a78
        a78 = 'cccc dtbat = srffrq = surface output interval (hrs)'
        write(99,105)a78
        a78 = 'cccc dtrad = radfrq = radiation output interval (hrs)'
        write(99,105)a78
        a78 = 'cccc DIRECT ACCESS BINARY SPECIFICATION'
        write(99,105)a78
        a78 = 'cccc ibyte = 4 for PGI, IFC; 1 for SUN, SGI, DEC, IBM'
        write(99,105)a78
        a78 = '      integer nxf,nyf,nz,nxs,nys'
        write(99,105)a78
        write(a78,182) ix
 182    format('      parameter (nyf =',i4,')')
        write(99,105)a78
#ifdef MPP1
        write(a78,183) mjx
#else
        write(a78,183) jx
#endif
 183    format('      parameter (nxf =',i4,')')
        write(99,105)a78
        write(a78,184) kx
 184    format('      parameter (nz  =',i4,')')
        write(99,105)a78
        write(a78,185) 1 ! nxs when SUBBATS si in
 185    format('      parameter (nxs =',i4,')')
        write(99,105)a78
        write(a78,186) 1 ! nys when SUBBATS si in
 186    format('      parameter (nys =',i4,')')
        write(99,105)a78
        a78 = '      integer ibyte'
        write(99,105)a78
        write(a78,187) ibyte
 187    format('      parameter (ibyte =',i4,')')
        write(99,105)a78
        a78 = '      real dtbc,dtout,dtbat,dtrad,dtche,dtsub'
        write(99,105)a78
        write(a78,188) ibdyfrq
 188    format('      parameter (dtbc  =',i7,'.00)')
        write(99,105)a78
        write(a78,189) dto
 189    format('      parameter (dtout =',f10.2,')')
        write(99,105)a78
        write(a78,190) dtb
 190    format('      parameter (dtbat =',f10.2,')')
        write(99,105)a78
        write(a78,191) dtr
 191    format('      parameter (dtrad =',f10.2,')')
        write(99,105)a78
        write(a78,192) dtc
 192    format('      parameter (dtche =',f10.2,')')
        write(99,105)a78
        write(a78,193) dtb
 193    format('      parameter (dtsub =',f10.2,')')
        write(99,105)a78
        close(99)
      end if

      return
      end

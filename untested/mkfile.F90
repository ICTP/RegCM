!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
      subroutine mkfile
 
      use mod_regcm_param
      use mod_param1
      use mod_param2
      use mod_iunits
      use mod_date
      use mod_bats
      use mod_mainchem
      use mod_outrad
      implicit none
!
! PARAMETER definitions
!
      integer , parameter :: nfmax = 999
!
! Local variables
!
      character(1) :: a1
      character(120) :: a120
      character(16) :: a16
      character(2) :: a2
      character(22) :: a22
      character(78) :: a78
      real(4) :: dtb , dtc , dto , dtr
      character(14) :: filatm , filchem , fillak , filrad , filsrf ,    &
                     & filsub
      character(40) :: filnam1 , filnam2
      integer , dimension(nfmax) :: idate1d , imo , iyr
      integer :: idatepp , isystm , iym1 , iym2 , n , nmo
      character(3) :: itype
      logical :: there
!
      print * , ' '
      print * , '******* OPENING NEW OUTPUT FILES:' , idatex
      if ( iftape ) then
        close (iutdat)
        itype = 'ATM'
        write (filatm,99001) itype , idatex
        if ( iotyp.eq.1 ) then
          open (iutdat,file='output/'//filatm,status='unknown',         &
               &form='unformatted',recl=(ixm2)*(jxm2)*ibyte,            &
               &access='direct')
          nrcout = 0
        else if ( iotyp.eq.2 ) then
          open (iutdat,file='output/'//filatm,status='unknown',         &
               &form='unformatted')
        else
        end if
        call gradsout(filatm//'.ctl')
        print * , 'OPENING NEW OUT FILE: output/' , filatm
      end if
 
      if ( ifbat ) then
        close (iutbat)
        itype = 'SRF'
        write (filsrf,99001) itype , idatex
        if ( iotyp.eq.1 ) then
          open (iutbat,file='output/'//filsrf,status='unknown',         &
               &form='unformatted',recl=(ixm2)*(jxm2)*ibyte,            &
               &access='direct')
          nrcbat = 0
        else if ( iotyp.eq.2 ) then
          open (iutbat,file='output/'//filsrf,status='unknown',         &
               &form='unformatted')
        else
        end if
        call gradsbat(filsrf//'.ctl')
        print * , 'OPENING NEW BAT FILE: output/' , filsrf
      end if
 
      if ( nsg.gt.1 .and. ifsub ) then
        close (iutsub)
        itype = 'SUB'
        write (filsub,99001) itype , idatex
        if ( iotyp.eq.1 ) then
          open (iutsub,file='output/'//filsub,status='unknown',         &
               &form='unformatted',recl=(ixm2)*(jxm2)*nnsg*ibyte,       &
               &access='direct')
          nrcsub = 0
        else if ( iotyp.eq.2 ) then
          open (iutsub,file='output/'//filsub,status='unknown',         &
               &form='unformatted')
        else
        end if
        call gradssub(filsub//'.ctl')
        print * , 'OPENING NEW SUB FILE: output/' , filsub
      end if
 
      if ( ifrad ) then
        close (iutrad)
        itype = 'RAD'
        write (filrad,99001) itype , idatex
        if ( iotyp.eq.1 ) then
          open (iutrad,file='output/'//filrad,status='unknown',         &
               &form='unformatted',recl=(ixm2)*(jxm2)*ibyte,            &
               &access='direct')
          nrcrad = 0
        else if ( iotyp.eq.2 ) then
          open (iutrad,file='output/'//filrad,status='unknown',         &
               &form='unformatted')
        else
        end if
        call gradsrad(filrad//'.ctl')
        print * , 'OPENING NEW RAD FILE: output/' , filrad
      end if
 
!chem2
      if ( ichem.eq.1 ) then
        if ( ifchem ) then
          close (iutchem)
          itype = 'CHE'
          write (filchem,99001) itype , idatex
          if ( iotyp.eq.1 ) then
            open (iutchem,file='output/'//filchem,status='unknown',     &
                 &form='unformatted',recl=(ixm2)*(jxm2)*ibyte,          &
                 &access='direct')
            nrcchem = 0
          else if ( iotyp.eq.2 ) then
            open (iutchem,file='output/'//filchem,status='unknown',     &
                 &form='unformatted')
          else
          end if
          call gradschem(filchem//'.ctl')
          print * , 'OPENING NEW CHEM FILE: output/' , filchem
        end if
      end if
!chem2_
      if ( lakemod.eq.1 ) then
        close (iutlak)
        itype = 'LAK'
        write (fillak,99001) itype , idatex
        open (iutlak,file='output/'//fillak,status='unknown',           &
             &form='unformatted')
        print * , 'OPENING NEW LAK FILE: ' , fillak
      end if
 
      if ( jyear.eq.jyear0 .and. ktau.eq.0 ) then
        inquire (file='postproc.in',exist=there)
        if ( there ) isystm = system(                                   &
                             &'/bin/mv -f postproc.in postproc.in.bak')
        nmo = (idate2/1000000-idate0/1000000)                           &
            & *12 + (idate2/10000-(idate2/1000000)*100)                 &
            & - (idate0/10000-(idate0/1000000)*100)
        nmo = min(nmo,nfmax)
        idatepp = idate0
        iyr(1) = idatepp/1000000
        imo(1) = (idatepp-iyr(1)*1000000)/10000
        idate1d(1) = idatepp
        do n = 2 , nmo
          idatepp = (idatepp/10000)*10000 + 10100
          iyr(n) = idatepp/1000000
          imo(n) = (idatepp-iyr(n)*1000000)/10000
          if ( imo(n).gt.12 ) then
            iyr(n) = iyr(n) + 1
            imo(n) = 1
            idatepp = iyr(n)*1000000 + imo(n)*10000 + 100
          end if
          idate1d(n) = idatepp
        end do
        if ( nmo.le.1 ) then
          nmo = 1
          idate1d(2) = idate0
        end if
 
!       **** Write out postprocessing information
        open (99,file='postproc.in',form='FORMATTED',status='unknown')
        a1 = char(39)
        write (99,99002) idate1d(2)
        write (99,99003) idate1d(2)
        write (99,99004) idate2
        a78 = '2,              ! iotype: 1=I*2 NetCDF;'//               &
             &' 2=r*4 NetCDF; 3=grads; 4=vis5d'
        write (99,99008) a78
        a78 = '.false.,        ! Write out header?'
        write (99,99008) a78
        a78 = '.false.,        ! Write out all RegCM data b/twn'//      &
             &' idate1 & idate2?'
        write (99,99008) a78
        a78 = '.false.,        ! Average RegCM data b/twn idate1'//     &
             &' & idate2?'
        write (99,99008) a78
        a78 = '.false.,        ! Diurnali avg of RegCM data b/twn'//    &
             &' idate1 & idate2?'
        write (99,99008) a78
        a78 = '.true.,         ! Continually average ATM data b/twn'//  &
             &' idate1 & idate2?'
        write (99,99008) a78
        a78 = '-1.,            ! No. Days for continual averaging'//    &
             &' (-1=monthly;1=daily;5=5day)'
        write (99,99008) a78
        if ( nmo.le.1 ) then
          iym1 = iyr(nmo)*100 + imo(nmo)
          write (a16,99005) a1 , iym1 , a1
        else if ( nmo.le.2 ) then
          iym1 = iyr(nmo)*100 + imo(nmo)
          write (a16,99005) a1 , iym1 , a1
        else
          iym1 = iyr(2)*100 + imo(2)
          iym2 = iyr(nmo)*100 + imo(nmo)
          write (a16,99006) a1 , iym1 , iym2 , a1
        end if
        a78 = a16//'! postproc output filename (not'//                  &
             &' including type & ext)'
        write (99,99008) a78
 
        a78 = a1//'../Input'//a1//',     ! ICBC directory'
        write (99,99008) a78
        a78 = a1//'output'//a1//',       ! RegCM Output directory'
        write (99,99008) a78
        a78 = a1//'DOMAIN.INFO'//a1//                                   &
             &',  ! Domain Info Filename (from terrain)'
        write (99,99008) a78
        a78 = a1//'OUT_HEAD'//a1//                                      &
             &',     ! Header File name (from RegCM)'
        write (99,99008) a78
        if ( nmo.le.1 ) then
          a2 = 'st'
          write (a22,99007) a1 , idate1d(1) , a1 , nmo , a2
          a78 = a22//' RegCM Output File Extension'
          write (99,99008) a78
        else
          do n = 2 , nmo
            if ( n.eq.2 ) then
              a2 = 'st'
            else if ( n.eq.3 ) then
              a2 = 'nd'
            else if ( n.eq.4 ) then
              a2 = 'rd'
            else
              a2 = 'th'
            end if
            write (a22,99007) a1 , idate1d(n) , a1 , n - 1 , a2
            a78 = a22//'RegCM Output File Extension'
            write (99,99008) a78
          end do
        end if
        close (98)
        close (99)
 
        inquire (file='../PostProc',exist=there)
        if ( .not.there ) isystm = system('mkdir ../PostProc')
        filnam1 = '../PostProc/postproc.param'
        filnam2 = '../PostProc/postproc.param.bak'
        inquire (file=filnam1,exist=there)
        if ( there ) then
          a120 = '/bin/mv -f '//filnam1//' '//filnam2
          isystm = system(a120)
        end if
        open (99,file=filnam1,form='FORMATTED',status='unknown')
        dto = tapfrq
        dtb = batfrq
        dtr = radisp
        dtc = chemfrq
        a78 = 'cccc SET DOMAIN DIMENSIONS'
        write (99,99008) a78
        a78 = 'cccc ny = number of north-south points'
        write (99,99008) a78
        a78 = 'cccc nx = number of east-west points'
        write (99,99008) a78
        a78 = 'cccc nz = number of vertical levels'
        write (99,99008) a78
        a78 = 'cccc SET MODEL OUTPUT INTERVALS'
        write (99,99008) a78
        a78 = 'cccc dtbc = ibdyfrq = ICBC output interval (hrs)'
        write (99,99008) a78
        a78 = 'cccc dtout = atmfrq = atmospheric output interval (hrs)'
        write (99,99008) a78
        a78 = 'cccc dtbat = srffrq = surface output interval (hrs)'
        write (99,99008) a78
        a78 = 'cccc dtrad = radfrq = radiation output interval (hrs)'
        write (99,99008) a78
        a78 = 'cccc DIRECT ACCESS BINARY SPECIFICATION'
        write (99,99008) a78
        a78 = 'cccc ibyte = 4 for PGI, IFC; 1 for SUN, SGI, DEC, IBM'
        write (99,99008) a78
        a78 = '      integer nxf,nyf,nz,nxs,nys'
        write (99,99008) a78
        write (a78,99009) ix
        write (99,99008) a78
        write (a78,99010) jx
        write (99,99008) a78
        write (a78,99011) kx
        write (99,99008) a78
        write (a78,99012) 1
                         ! nxs when SUBBATS si in
        write (99,99008) a78
        write (a78,99013) 1
                         ! nys when SUBBATS si in
        write (99,99008) a78
        a78 = '      integer ibyte'
        write (99,99008) a78
        write (a78,99014) ibyte
        write (99,99008) a78
        a78 = '      real dtbc,dtout,dtbat,dtrad,dtche,dtsub'
        write (99,99008) a78
        write (a78,99015) ibdyfrq
        write (99,99008) a78
        write (a78,99016) dto
        write (99,99008) a78
        write (a78,99017) dtb
        write (99,99008) a78
        write (a78,99018) dtr
        write (99,99008) a78
        write (a78,99019) dtc
        write (99,99008) a78
        write (a78,99020) dtb
        write (99,99008) a78
        close (99)
      end if
 
99001 format (a3,'.',i10)
99002 format (i10,',     ! idate0 = First date in File (yymmddhh)')
99003 format (i10,                                                      &
             &',     ! idate1 = Start date for averaging and re-writing'&
            & )
99004 format (i10,                                                      &
             &',     ! idate2 = End date for averaging and re-writing')
99005 format (a1,i6,a1,',')
99006 format (a1,i6,'_',i6,a1,',')
99007 format (a1,i10,a1,',   !',i3,a2)
 
99008 format (a78)
99009 format ('      parameter (nyf =',i4,')')
99010 format ('      parameter (nxf =',i4,')')
99011 format ('      parameter (nz  =',i4,')')
99012 format ('      parameter (nxs =',i4,')')
99013 format ('      parameter (nys =',i4,')')
99014 format ('      parameter (ibyte =',i4,')')
99015 format ('      parameter (dtbc  =',i7,'.00)')
99016 format ('      parameter (dtout =',f10.2,')')
99017 format ('      parameter (dtbat =',f10.2,')')
99018 format ('      parameter (dtrad =',f10.2,')')
99019 format ('      parameter (dtche =',f10.2,')')
99020 format ('      parameter (dtsub =',f10.2,')')
 
      end subroutine mkfile

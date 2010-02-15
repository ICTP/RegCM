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
 
      subroutine for_next
      use mod_regcm_param
      use mod_param1
      use mod_param2
      use mod_param3 , only : jxsex , kxout, ncld
      use mod_date
      use mod_convect
      use mod_pmoist
      use mod_trachem
      implicit none
!
! Local variables
!
      character(80) :: chname
      integer :: itr , nc
!
      open (99,file='regcm.in.new',form='FORMATTED',status='unknown')
 
      write (99,99002) '&restartparam'
      ifrest = .true.
      if ( ifrest ) then
        write (99,99006) 'ifrest  = .true. '
      else
        write (99,99006) 'ifrest  = .false.'
      end if
      write (99,99007) 'idate0  = ' , idate0
      write (99,99007) 'idate1  = ' , idate1
      write (99,99007) 'idate2  = ' , idate2
      write (99,99007) 'nslice  = ' , nslice
      write (99,99001) '/'
 
      write (99,99003) '&timeparam'
      write (99,99008) 'radfrq  = ' , radfrq
      write (99,99008) 'abemh   = ' , abemh
      write (99,99008) 'abatm   = ' , abatm
      write (99,99008) 'dt      = ' , dtmin*60.
      write (99,99007) 'ibdyfrq = ' , ibdyfrq
      write (99,99001) '/'
 
      write (99,99004) '&outparam'
      if ( ifsave ) then
        write (99,99006) 'ifsave  = .true. '
      else
        write (99,99006) 'ifsave  = .false.'
      end if
      write (99,99009) 'savfrq  = ' , savfrq
      if ( iftape ) then
        write (99,99006) 'iftape  = .true. '
      else
        write (99,99006) 'iftape  = .false.'
      end if
      write (99,99009) 'tapfrq  = ' , tapfrq
      if ( ifrad ) then
        write (99,99006) 'ifrad   = .true. '
      else
        write (99,99006) 'ifrad   = .false.'
      end if
      write (99,99009) 'radisp  = ' , radisp
      if ( ifbat ) then
        write (99,99006) 'ifbat   = .true. '
      else
        write (99,99006) 'ifbat   = .false.'
      end if
      if ( ifsub ) then
        write (99,99006) 'ifsub   = .true. '
      else
        write (99,99006) 'ifsub   = .false.'
      end if
      write (99,99009) 'batfrq  = ' , batfrq
      if ( ifprt ) then
        write (99,99006) 'ifprt   = .true. '
      else
        write (99,99006) 'ifprt   = .false.'
      end if
      write (99,99009) 'prtfrq  = ' , prtfrq
      write (99,99010) 'kxout   = ' , kxout
      write (99,99010) 'jxsex   = ' , jxsex
      write (99,99007) 'iotyp   = ' , iotyp
      write (99,99007) 'ibintyp = ' , ibintyp
      if ( ifchem ) then
        write (99,99006) 'ifchem  = .true. '
      else
        write (99,99006) 'ifchem  = .false.'
      end if
      write (99,99009) 'chemfrq = ' , chemfrq
      write (99,99001) '/'
 
      write (99,99002) '&physicsparam'
      write (99,99007) 'iboudy  = ' , iboudy
      write (99,99007) 'ibltyp  = ' , ibltyp
      write (99,99007) 'icup    = ' , icup
      write (99,99010) 'igcc    = ' , igcc
      write (99,99007) 'ipptls  = ' , ipptls
      write (99,99007) 'iocnflx = ' , iocnflx
      write (99,99007) 'ipgf    = ' , ipgf
      write (99,99007) 'iemiss  = ' , iemiss
      write (99,99007) 'lakemod = ' , lakemod
      write (99,99007) 'ichem   = ' , ichem
      write (99,99001) '/'
 
      write (99,99005) '&subexparam'
      if ( ncld.ne.1 ) write (99,99007) 'ncld    = ' , ncld
      if ( fcmax.ne.0.80 ) write (99,99011) 'fcmax   = ' , fcmax
      if ( qck1land.ne.0.0005 ) write (99,99012) 'qck1land= ' , qck1land
      if ( qck1oce.ne.0.0005 ) write (99,99012) 'qck1oce = ' , qck1oce
      if ( gulland.ne.0.4 ) write (99,99011) 'gulland = ' , gulland
      if ( guloce.ne.0.4 ) write (99,99011) 'guloce  = ' , guloce
      if ( rhmax.ne.1.01 ) write (99,99011) 'rhmax   = ' , rhmax
      if ( rh0oce.ne.0.90 ) write (99,99011) 'rh0oce  = ' , rh0oce
      if ( rh0land.ne.0.80 ) write (99,99011) 'rh0land = ' , rh0land
      if ( tc0.ne.238.0 ) write (99,99011) 'tc0     = ' , tc0
      if ( cevap.ne.0.00002 ) write (99,99012) 'cevap   = ' , cevap
      if ( caccr.ne.6.0 ) write (99,99011) 'caccr   = ' , caccr
      write (99,99001) '/'
 
      write (99,99005) '&grellparam'
      if ( shrmin.ne.0.25 ) write (99,99011) 'shrmin  = ' , shrmin
      if ( shrmax.ne.0.50 ) write (99,99011) 'shrmax  = ' , shrmax
      if ( edtmin.ne.0.25 ) write (99,99011) 'edtmin  = ' , edtmin
      if ( edtmax.ne.1.00 ) write (99,99011) 'edtmax  = ' , edtmax
      if ( edtmino.ne.0.00 ) write (99,99011) 'edtmino = ' , edtmino
      if ( edtmaxo.ne.1.00 ) write (99,99011) 'edtmaxo = ' , edtmaxo
      if ( edtminx.ne.0.25 ) write (99,99011) 'edtminx = ' , edtminx
      if ( edtmaxx.ne.1.00 ) write (99,99011) 'edtmaxx = ' , edtmaxx
      if ( pbcmax.ne.150. ) write (99,99011) 'pbcmax  = ' , pbcmax
      if ( mincld.ne.150. ) write (99,99011) 'mincld  = ' , mincld
      if ( htmin.ne.-250. ) write (99,99011) 'htmin   = ' , htmin
      if ( htmax.ne.500. ) write (99,99011) 'htmax   = ' , htmax
      if ( skbmax.ne.0.4 ) write (99,99011) 'skbmax  = ' , skbmax
      if ( dtauc.ne.30. ) write (99,99011) 'dtauc   = ' , dtauc
      write (99,99001) '/'
 
      write (99,99003) '&emanparam'
      if ( minsig.ne.0.95 ) write (99,99011) 'minsig  = ' , minsig
      if ( elcrit.ne.0.0011 ) write (99,99011) 'elcrit  = ' , elcrit
      if ( tlcrit.ne.-55.0 ) write (99,99011) 'tlcrit  = ' , tlcrit
      if ( entp.ne.1.5 ) write (99,99011) 'entp    = ' , entp
      if ( sigd.ne.0.05 ) write (99,99011) 'sigd    = ' , sigd
      if ( sigs.ne.0.12 ) write (99,99011) 'sigs    = ' , sigs
      if ( omtrain.ne.50.0 ) write (99,99011) 'omtrain = ' , omtrain
      if ( omtsnow.ne.5.5 ) write (99,99011) 'omtsnow = ' , omtsnow
      if ( coeffr.ne.1.0 ) write (99,99011) 'coeffr  = ' , coeffr
      if ( coeffs.ne.0.8 ) write (99,99011) 'coeffs  = ' , coeffs
      if ( cu.ne.0.7 ) write (99,99011) 'cu      = ' , cu
      if ( betae.ne.10.0 ) write (99,99011) 'betae   = ' , betae
      if ( dtmax.ne.0.9 ) write (99,99011) 'dtmax   = ' , dtmax
      if ( alphae.ne.0.2 ) write (99,99011) 'alphae  = ' , alphae
      if ( damp.ne.0.1 ) write (99,99011) 'damp    = ' , damp
      write (99,99001) '/'
 
      write (99,99005) '&chemparam'
      if ( ichremlsc.ne.1 ) write (99,99013) 'ichremlsc = ' , ichremlsc
      if ( ichremcvc.ne.1 ) write (99,99013) 'ichremcvc = ' , ichremcvc
      if ( ichdrdepo.ne.1 ) write (99,99013) 'ichdrdepo = ' , ichdrdepo
      if ( ichcumtra.ne.1 ) write (99,99013) 'ichcumtra = ' , ichcumtra
      write (99,99013) 'idirect   = ' , idirect
      nc = 0
 
      do itr = 1 , ntr
        if ( chtrname(itr).eq.'SO2' .or. chtrname(itr).eq.'SO4' ) then
          write (chname(nc+1:nc+6),99015) chtrname(itr)
          nc = nc + 6
        else if ( chtrname(itr).eq.'DUST' ) then
          write (chname(nc+1:nc+7),99016) chtrname(itr)
          nc = nc + 7
        else
          write (chname(nc+1:nc+8),99017) chtrname(itr)
          nc = nc + 8
        end if
      end do
      write (99,99014) 'chtrname  = ' , chname(1:nc)
      write (99,99018) 'chtrsol   = ' , chtrsol
      write (99,99019) 'chtrdpv   = ' , chtrdpv
      write (99,99018) 'dustbsiz  = ' , dustbsiz
      write (99,99001) '/'
 
      close (99)
99001 format (1x,a1)
99002 format (1x,a13)
99003 format (1x,a10)
99004 format (1x,a9)
99005 format (1x,a11)
99006 format (1x,a17,',')
99007 format (1x,a10,i10,',')
99008 format (1x,a10,f6.0,',')
99009 format (3x,a10,f6.0,',')
99010 format (3x,a10,i10,',')
99011 format (1x,a10,f10.3,',')
99012 format (1x,a10,e10.3,',')
99013 format (1x,a12,i4,',')
99014 format (1x,a12,a75)
99015 format ("'",a3,"',")
99016 format ("'",a4,"',")
99017 format ("'",a5,"',")
99018 format (1x,a12,10(f5.2,','))
99019 format (1x,a12,20(f7.5,','))
      end subroutine for_next

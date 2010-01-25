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

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c-----------------------------NOTICE------------------------------------
c
c            NCAR COMMUNITY CLIMATE MODEL, VERSION 3.0
c            COPYRIGHT (C) 1996
c            UNIVERSITY CORPORATION FOR ATMOSPHERIC RESEARCH
c            ALL RIGHTS RESERVED
c
c               ------------ ----- --- ---------- ------
c  ********** | distribution terms and conditions notice | ************
c               ------------ ----- --- ---------- ------
c
c (c) copyright 1996 university corporation for atmospheric research/
c national center for atmospheric research/
c climate and global dynamics division
c
c this software, the community climate model (ccm), version ccm3, was
c developed by the climate and global dynamics division (cgd) climate
c modeling section (cms) of the national center for atmospheric research
c (ncar), which is operated by the university corporation for
c atmospheric research (ucar) and sponsored by the national science
c foundation (nsf).
c
c access and use of this software shall impose the following obligations
c and understandings on the user.  the user is granted the right,
c without any fee or cost, to use, copy, modify, alter, enhance and
c distribute this software, and any derivative works thereof, and its
c supporting documentation for any purpose whatsoever, except commercial
c sales, provided that this entire notice appears in all copies of the
c software, derivative works and supporting documentation.  further, the
c user agrees to credit ucar/ncar/cgd in any publications that result
c from the use of this software or in any software package that includes
c this software.  the names ucar/ncar/cgd, however, may not be used in
c any advertising or publicity to endorse or promote any products or
c commercial entity unless specific written permission is obtained from
c ucar/ncar/cgd.
c
c the ccm3 materials are made available with the understanding that
c ucar/ncar/cgd is not obligated to provide (and will not provide) the
c user with any support, consulting, training, or assistance of any kind
c with regard to the use, operation and performance of this software, nor
c to provide the user with any updates, revisions, new versions, or "bug
c fixes."
c
c this software is provided by ucar/ncar/cgd "as is" and any express or
c implied warranties, including but not limited to, the implied
c warranties of merchantability and fitness for a particular purpose are
c disclaimed.  in no event shall ucar/ncar/cgd be liable for any
c special, indirect or consequential damages or any damages whatsoever,
c including but not limited to claims associated with the loss of data
c or profits, which may result from an action in contract, negligence or
c other tortious claim that arises out of or in connection with the
c access, use or performance of this software.
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
      subroutine colmod3(jslc)
c
c     ccm3 column radiation model (crm)
c     Jeffrey Kiehl, Bruce Briegleb, and Charlie Zender
c     May 1996
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     This radiation model has been implemented in regcm3
c     by Keiichi Nishizawa on leave from CRIEPI (Japan)
c     in April, 1997
c
c     NB:
c     The following comments have not been changed from the original
c     in the ccm3 column radiation model (crm)
c     Therefore, they are not necessarily appropriate
c     for the regcm3 radiation package
c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c
c     All routines except crm(), getdat(), radctl(), and four
c     dummy routines (readric(), writeric(), outfld(), and radozn()) are
c     included straight from CCM3 code. The purpose of the (non-dummy)
c     routines is:
c
c     crm(): the main() routine. This routine takes the place of tphysbc()
c     in the ccm3 calling tree. It places the required calls to the cloud
c     routines cldefr() and cldems() directly, rather than invoking cldint().
c
c     getdat(): reads the stdin file to parse the user-specified column
c     profile. overwrites the co2vmr previsouly set by radini() with the
c     user specified value. computes o3vmr (instead of in radozn()).
c     sets the calday variable in comtim.h used in zenith() and radinp().
c
c     radctl(): main radiation driving routine, same as ccm3 version except:
c     receives additional input variable o3vmr, and passes out
c     additional diagnostic radiation quantities and eccf usually local to
c     radctl(). does all cgs-->mks conversion for all fluxes.
c
c
c     The files that need to be either -included-, specified in the
c     compilation statement, or linked from a ccm3 library are...
c
c     aermix.F    fmrgrid.F    radclr.F   radinp.F    trcabn.F    whenfgt.F
c     albland.F   intmax.F     radclw.F   radoz2.F    trcems.F    whenflt.F
c     albocean.F  isrchfgt.F   radcsw.F   radtpl.F    trcmix.F    whenne.F
c     blkdat.F    isrchfle.F   radded.F   resetr.F    trcplk.F    zenith.F
c     cldefr.F    myhandler.F  radems.F   torgrid.F   trcpth.F
c     cldems.F    radabs.F     radini.F   trcab.F     wheneq.F
c
c
c     Standard input file is a text (ascii) file. At the end of the included
c     input file are notes describing more fully the input. The standard
c     output file is also a text (ascii) file.
c
c $Id: implicit.h,v 1.1.1.1 1995/02/09 23:26:52 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
      integer jslc
c
cKN   added below
c
C------------------------------Parameters-------------------------------
#     include "parameter.inc"
#     include "include/parrad.cb"
#     include "include/param1.cb"
#     include "include/param2.cb"
#     include "include/iunits.cb"
#     include "include/bats.cb"
#     include "include/bats.cb2"
c
C------------------------------Commons----------------------------------
#     include "include/comtim.cb"
#     include "include/comctl.cb"
C------------------------------Arguments--------------------------------
c
c     Fields specified by the user in getdat()
c
      integer ioro(plond)       ! land/ocean/sea ice flag
c
      real(kind=8)  clat(plond),                ! Current latitude (radians)
     $     cld(plond,plevp),    ! fractional cloud cover
     $     clwp(plond,plev),    ! cloud liquid water path
     $     coslat(plond),              ! cosine latitude
c
c     NB: o3mmr and o3vmr should be dimensioned (plond,plevr) if a different
c     size radiation grid is used. Clashes between prgrid.h and ptrrgrid.h
c     (they both define plngbuf) prevent us from dimensioning anything by
c     plevr in this top level crm() routine.
c
cKN   added below
     $     h2ommr(plond,plev),  ! water vapor mass mixing ratio
cKN   added above
     $     o3mmr(plond,plev),  ! Ozone mass mixing ratio
     $     o3vmr(plond,plev),  ! Ozone volume mixing ratio
     $     pilnm1(plond,plevp), ! natural log of pintm1
     $     pintm1(plond,plevp), ! model interface pressures
     $     pmidm1(plond,plev),  ! model level pressures
     $     pmlnm1(plond,plev),  ! natural log of pmidm1
cKN     $     ps(plond),           ! surface pressure
     $     qm1(plond,plev),     ! model level specific humidity
     $     tm1(plond,plev),     ! model level temperatures
     $     ts(plond)            ! surface air temperature
c
c     Fields computed from user input
c
      real(kind=8)  effcld(plond,plevp), ! effective cloud=cld*emis
     $     emis(plond,plev),    ! cloud emissivity
     $     fice(plond,plev),    ! fractional amount of ice
     $     rei(plond,plev),     ! ice particle size
     $     rel(plond,plev)      ! liquid effective drop size (microns)
      real(kind=8)  eccf,                ! earth/sun distance factor
     $     loctim(plond),       ! local time of solar computation
     $     srfrad(plond)        ! srf radiative heat flux
c
      real(kind=8)  asdir(plond),        ! albedo: shortwave, direct
     $     asdif(plond),        ! albedo: shortwave, diffuse
     $     aldir(plond),        ! albedo: longwave, direct
     $     aldif(plond)         ! albedo: longwave, diffuse
C
C     Output longwave arguments from radctl()
C
      real(kind=8)  flwds(plond),        ! Surface down longwave flux
     $     qrl(plond,plev)      ! Longwave cooling rate
C
C     Output shortwave arguments from radctl()
C
      real(kind=8)  fsns(plond),         ! Surface absorbed solar flux
     $     qrs(plond,plev),     ! Solar heating rate
     $     soll(plond),         ! Downward solar rad onto surface (lw direct)
     $     solld(plond),        ! Downward solar rad onto surface (lw diffuse)
     $     sols(plond),         ! Downward solar rad onto surface (sw direct)
     $     solsd(plond)         ! Downward solar rad onto surface (sw diffuse)
c
c     Additional CRM diagnostic output from radctl()
c
      real(kind=8)  flns(plond),         ! srf longwave cooling (up-dwn) flux
     $     flnsc(plond),        ! clr sky lw flx at srf (up-dwn)
     $     flnt(plond),         ! net outgoing lw flx at model top
     $     flntc(plond),        ! clr sky lw flx at model top
     $     fsnsc(plond),        ! clr sky surface abs solar flux
     $     fsnt(plond),         ! total column absorbed solar flux
     $     fsntc(plond),        ! clr sky total column abs solar flux
     $     solin(plond)         ! solar incident flux
c
CEES next 3 added, they are calculated in radcsw
      real(kind=8) fsnirt(plond)      ! Near-IR flux absorbed at toa
      real(kind=8) fsnrtc(plond)      ! Clear sky near-IR flux absorbed at toa
      real(kind=8) fsnirtsq(plond)    ! Near-IR flux absorbed at toa >= 0.7 microns
      real(kind=8) fsds(plond)        ! Flux Shortwave Downwelling Surface
      real(kind=8) alb(plond)
      real(kind=8) albc(plond)
c
      integer
     $     i,                   ! longitude index
     $     k                    ! level index
      integer n
c
c     Fundamental constants needed by radini()
c
      real(kind=8)  cpairx,              ! heat capacity dry air at constant prs (J/kg/K)
     $     epsilox,             ! ratio mean mol weight h2o to dry air
     $     gravx,               ! gravitational acceleration (m/s**2)
     $     stebolx              ! Sefan-Boltzmann constant (W/m**2/K**4)
      integer ii0,ii1,ii2
c
c     Set latitude index to jslc
c
      lat = jslc
c
c     Set parameters in common block comtim.h: nstep,dosw,dolw,doabsems
c
c     nstep  =  ktau
      dosw  =  .true.
      dolw  =  .true.
      doabsems = .true.
c
c     Set parameters in common block comctl.h: anncyc,dodiavg,iradsw,iradlw,iradae
c
      anncyc = .true.
      dodiavg = .false.
      iradsw = 1
      iradlw = 1
      iradae = 1
c
c     Set parameters required for call to radini():
c
      gravx   =   9.80616
      cpairx  =   1.00464e3
      epsilox =   0.622
      stebolx =   5.67e-8
c
c     Given these four constants in MKS, radini() will define their CGS
c     equivalents, as well as setting many radiation parameters stored
c     in common blocks. radini() must be called before getdat(), because
c     the co2 mixing ratio set (by the user) in getdat() should overwrite
c     the default CCM3 co2 mixing ratio set by radini().
c
      call radini(gravx,cpairx,epsilox,stebolx)
c
cKN      write(6,*) ' .... Read in profile data ....'
c
cKN   added below
c
c     NB: orography types are specified in the following
c
      do i=1,plon
        ii0=0
        ii1=0
        ii2=0
        do n=1,NNSG
         if (ldoc1d(n,i).gt.0.1.and.sice1d(n,i).eq.0.0) then
           ii1=ii1+1
         else if (sice1d(n,i).gt.0.0) then
           ii2=ii2+1
         else
           ii0=ii0+1
         end if
        enddo
        if(ii0.ge.ii1.and.ii0.ge.ii2) ioro(i)=0
        if(ii1.gt.ii0.and.ii1.ge.ii2) ioro(i)=1
        if(ii2.gt.ii0.and.ii2.gt.ii1) ioro(i)=2
      end do
c
cKN   added above
c
c     getdat() also sets calday (used in zenith() and radinp()).
c
cKN   modified below
      call getdat(jslc,h2ommr,
cKN   modified above
     $     clat,
     $     cld,
     $     clwp,
     $     coslat,
     $     loctim,
     $     o3mmr,
     $     o3vmr,
     $     pilnm1,
     $     pintm1,
     $     pmidm1,
     $     pmlnm1,
     $     ps,
     $     qm1,
     $     tm1,
     $     ts)
c
cKN   modified below
c
c     NB:
c     variable coszrs is not calculated here in zenith()
c     but passed from vecbats() to colmod3() as an argument
c     therefore, subroutine zenith() is not necessary in regcm3
c
c     Get coszrs: needed as input to albland(), albocean(), radctl()
c
cKN      call zenith (calday  ,dodiavg ,clat    ,coszrs  )
c
c     NB:
c     land and ocean albedos are calculated
c     not in albland() and albocean() located below
c     but in albedov() called from subroutine vecbats()
c     therefore, the following subroutines are not called here
c
c     Find the albedo for land points
c
cKN      call albland(lat     ,ioro    ,sndpth  ,coszrs  ,asdir   ,
cKN     $     aldir   ,asdif   ,aldif   )
c
c     Find the albedo for ocean/sea-ice points
c
cKN      call albocean(lat     ,ioro    ,sndpth  ,coszrs  ,asdir   ,
cKN     $     aldir   ,asdif   ,aldif   )
c
c     NB:
c     albedos are copied from common block /bat1d/ in 'bats.cb2',
c     because variable names for albedos are somewhat different
c
      do i=1,plon
         asdir(i) = aldirs(i)
         asdif(i) = aldifs(i)
         aldir(i) = aldirl(i)
         aldif(i) = aldifl(i)
      end do
c
cKN   modified above
C
C Cloud particle size and fraction of ice
C
      call cldefr(ioro, tm1, rel, rei, fice, ps, pmidm1)
C
C Cloud emissivity
C
      call cldems(clwp, fice, rei, emis)
C
C Effective cloud cover
C
      do k=1,plev
         do i=1,plon
            effcld(i,k) = cld(i,k)*emis(i,k)
         end do
      end do
C
C Cloud cover at surface interface always zero (for safety's sake)
C
      do i=1,plon
         effcld(i,plevp) = 0.
         cld(i,plevp)    = 0.
      end do
C
C     Main radiation driving routine.
C     NB: All fluxes returned from radctl() have already been converted to MKS.
C
      call radctl(jslc      ,clat    ,coslat  ,ts   ,
     $     pmidm1  ,pintm1  ,pmlnm1  ,pilnm1  ,tm1     ,
     $     qm1     ,cld     ,effcld  ,clwp    ,
     $     asdir   ,asdif   ,aldir   ,aldif   ,fsns    ,
     $     qrs     ,qrl     ,flwds   ,rel     ,
     $     rei     ,fice    ,sols    ,soll    ,solsd   ,
     $     solld   ,emiss1d ,
c++csz
     $     fsnt,fsntc,fsnsc,flnt,flns,flntc,flnsc,solin, ! output
Cadd by bixq
     &     alb,albc,fsds,fsnirt,fsnrtc,fsnirtsq,
Cadd_
     $     eccf,
     $     o3vmr) ! input
c--csz
c
c write out final results:
c
cKN   added below
c
c     subroutine radout() is not included in the ccm3 crm itself
c     but introduced from the former regcm2 radiation package
c     for the output of results from radiation calculations
c     this subroutine is used also for the coupling with bats
c
c     NB:
c     Names of some output variables (in MKS) have been changed
c     from those in the CCM2 radiation package.
c
      call radout(solin  ,fsnt   ,fsns   ,fsntc  ,fsnsc  ,
     $            qrs    ,flnt   ,flns   ,flntc  ,flnsc  ,
     $            qrl, flwds ,srfrad,
     $            sols,soll,solsd,solld,
Cadd by bixq
     &     alb,albc,fsds,fsnirt,fsnrtc,fsnirtsq,
Cadd_
     $            jslc, h2ommr, cld, clwp)
c
cKN   added above
c
      return
      end

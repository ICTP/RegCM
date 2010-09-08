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
 
      module mod_colmod3

      use mod_constants
      use mod_dynparam
      use mod_bats
      use mod_date
      use mod_radiation

      private

      public :: colmod3

      logical :: aeregen , aeres , anncyc , cpuchek , dodiavg , lbrnch ,&
               & ldebug , nlend , nlhst , nlres , ozncyc , sstcyc
      integer :: iradae , iradlw , iradsw , itsst , nsrest

      contains

      subroutine colmod3(jslc)

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!-----------------------------NOTICE------------------------------------
!
!            NCAR COMMUNITY CLIMATE MODEL, VERSION 3.0
!            COPYRIGHT (C) 1996
!            UNIVERSITY CORPORATION FOR ATMOSPHERIC RESEARCH
!            ALL RIGHTS RESERVED
!
!               ------------ ----- --- ---------- ------
!  ********** | distribution terms and conditions notice | ************
!               ------------ ----- --- ---------- ------
!
! (c) copyright 1996 university corporation for atmospheric research/
! national center for atmospheric research/
! climate and global dynamics division
!
! this software, the community climate model (ccm), version ccm3, was
! developed by the climate and global dynamics division (cgd) climate
! modeling section (cms) of the national center for atmospheric research
! (ncar), which is operated by the university corporation for
! atmospheric research (ucar) and sponsored by the national science
! foundation (nsf).
!
! access and use of this software shall impose the following obligations
! and understandings on the user.  the user is granted the right,
! without any fee or cost, to use, copy, modify, alter, enhance and
! distribute this software, and any derivative works thereof, and its
! supporting documentation for any purpose whatsoever, except commercial
! sales, provided that this entire notice appears in all copies of the
! software, derivative works and supporting documentation.  further, the
! user agrees to credit ucar/ncar/cgd in any publications that result
! from the use of this software or in any software package that includes
! this software.  the names ucar/ncar/cgd, however, may not be used in
! any advertising or publicity to endorse or promote any products or
! commercial entity unless specific written permission is obtained from
! ucar/ncar/cgd.
!
! the ccm3 materials are made available with the understanding that
! ucar/ncar/cgd is not obligated to provide (and will not provide) the
! user with any support, consulting, training, or assistance of any kind
! with regard to the use, operation and performance of this software, nor
! to provide the user with any updates, revisions, new versions, or "bug
! fixes."
!
! this software is provided by ucar/ncar/cgd "as is" and any express or
! implied warranties, including but not limited to, the implied
! warranties of merchantability and fitness for a particular purpose are
! disclaimed.  in no event shall ucar/ncar/cgd be liable for any
! special, indirect or consequential damages or any damages whatsoever,
! including but not limited to claims associated with the loss of data
! or profits, which may result from an action in contract, negligence or
! other tortious claim that arises out of or in connection with the
! access, use or performance of this software.
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     ccm3 column radiation model (crm)
!     Jeffrey Kiehl, Bruce Briegleb, and Charlie Zender
!     May 1996
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     This radiation model has been implemented in regcm3
!     by Keiichi Nishizawa on leave from CRIEPI (Japan)
!     in April, 1997
!
!     NB:
!     The following comments have not been changed from the original
!     in the ccm3 column radiation model (crm)
!     Therefore, they are not necessarily appropriate
!     for the regcm3 radiation package
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!     All routines except crm(), getdat(), radctl(), and four
!     dummy routines (readric(), writeric(), outfld(), and radozn()) are
!     included straight from CCM3 code. The purpose of the (non-dummy)
!     routines is:
!
!     crm(): the main() routine. This routine takes the place of tphysbc()
!     in the ccm3 calling tree. It places the required calls to the cloud
!     routines cldefr() and cldems() directly, rather than invoking cldint().
!
!     getdat(): reads the stdin file to parse the user-specified column
!     profile. overwrites the co2vmr previsouly set by radini() with the
!     user specified value. computes o3vmr (instead of in radozn()).
!     sets the calday variable in comtim.h used in zenith() and radinp().
!
!     radctl(): main radiation driving routine, same as ccm3 version except:
!     receives additional input variable o3vmr, and passes out
!     additional diagnostic radiation quantities and eccf usually local to
!     radctl(). does all cgs-->mks conversion for all fluxes.
!
!
!     The files that need to be either -included-, specified in the
!     compilation statement, or linked from a ccm3 library are...
!
!     aermix.F    fmrgrid.F    radclr.F   radinp.F    trcabn.F    whenfgt.F
!     albland.F   intmax.F     radclw.F   radoz2.F    trcems.F    whenflt.F
!     albocean.F  isrchfgt.F   radcsw.F   radtpl.F    trcmix.F    whenne.F
!     blkdat.F    isrchfle.F   radded.F   resetr.F    trcplk.F    zenith.F
!     cldefr.F    myhandler.F  radems.F   torgrid.F   trcpth.F
!     cldems.F    radabs.F     radini.F   trcab.F     wheneq.F
!
!
!     Standard input file is a text (ascii) file. At the end of the included
!     input file are notes describing more fully the input. The standard
!     output file is also a text (ascii) file.
!
      implicit none
!
! Dummy arguments
!
      integer :: jslc
!
! Local variables
!
      real(8) , dimension(iym1) :: alb , albc , aldif , aldir , asdif ,&
                                  & asdir , alat , coslat , flns ,      &
                                  & flnsc , flnt , flntc , flwds ,      &
                                  & fsds , fsnirt , fsnirtsq , fsnrtc , &
                                  & fsns , fsnsc , fsnt , fsntc ,       &
                                  & loctim , solin , soll , solld ,     &
                                  & sols , solsd , srfrad , ts
      real(8) , dimension(iym1,kzp1) :: cld , effcld , pilnm1 , pintm1
      real(8) , dimension(iym1,kz) :: clwp , emis , fice , h2ommr ,  &
           & o3mmr , o3vmr , pmidm1 , pmlnm1 , qm1 , qrl , qrs , rei ,  &
           & rel , tm1
      real(8) :: eccf
      integer :: i , ii0 , ii1 , ii2 , k , n
      integer , dimension(iym1) :: ioro
!
!KN   added below
!
!
!     Fields specified by the user in getdat()
!
! land/ocean/sea ice flag
! Current latitude (radians)
! fractional cloud cover
! cloud liquid water path
!
!     NB: o3mmr and o3vmr should be dimensioned (iym1,kz) if a
!     different size radiation grid is used. Clashes between prgrid.h
!     and ptrrgrid.h (they both define plngbuf) prevent us from
!     dimensioning anything by kz in this top level crm() routine.
!
! cosine latitude
! water vapor mass mixing ratio
! Ozone mass mixing ratio
! Ozone volume mixing ratio
! natural log of pintm1
! model interface pressures
! model level pressures
! natural log of pmidm1
! model level specific humidity
! model level temperatures
!
!     Fields computed from user input
!
! surface air temperature
! effective cloud=cld*emis
! cloud emissivity
! fractional amount of ice
! ice particle size
! liquid effective drop size (microns)
! earth/sun distance factor
! local time of solar computation
! srf radiative heat flux
! albedo: shortwave, direct
! albedo: shortwave, diffuse
! albedo: longwave, direct
!
!     Output longwave arguments from radctl()
!
! albedo: longwave, diffuse
! Surface down longwave flux
!
!     Output shortwave arguments from radctl()
!
! Longwave cooling rate
! Surface absorbed solar flux
! Solar heating rate
! Downward solar rad onto surface (lw direct)
! Downward solar rad onto surface (lw diffuse)
! Downward solar rad onto surface (sw direct)
!
!     Additional CRM diagnostic output from radctl()
!
! Downward solar rad onto surface (sw diffuse)
! srf longwave cooling (up-dwn) flux
! clr sky lw flx at srf (up-dwn)
! net outgoing lw flx at model top
! clr sky lw flx at model top
! clr sky surface abs solar flux
! total column absorbed solar flux
! clr sky total column abs solar flux
!
!EES  next 3 added, they are calculated in radcsw
! solar incident flux
! Near-IR flux absorbed at toa
! Clear sky near-IR flux absorbed at toa
! Near-IR flux absorbed at toa >= 0.7 microns
! Flux Shortwave Downwelling Surface
!
!     Fundamental constants needed by radini()
!
! heat capacity dry air at constant prs (J/kg/K)
! ratio mean mol weight h2o to dry air
! gravitational acceleration (m/s**2)
! Sefan-Boltzmann constant (W/m**2/K**4)
!
!     Set latitude index to jslc
!
      ilat = jslc
!
!     Set parameters in common block comtim : dosw,dolw,doabsems
!
      dosw = .true.
      dolw = .true.
      doabsems = .true.
!
!     Set parameters in common block comctl:
!     anncyc,dodiavg,iradsw,iradlw,iradae
!
      anncyc = .true.
      dodiavg = .false.
      iradsw = 1
      iradlw = 1
      iradae = 1
!
!     radini sets many radiation parameters
!     radini() must be called before getdat(), because
!     the co2 mixing ratio set (by the user) in getdat() should
!     overwrite the default CCM3 co2 mixing ratio set by radini().
!
      call radini
!
!     NB: orography types are specified in the following
!
      do i = 1 , iym1
        ii0 = 0
        ii1 = 0
        ii2 = 0
        do n = 1 , nnsg
#ifdef CLM
          if ( ocld2d(n,i,jslc).gt.0.1 .and. sice1d(n,i).eq.0.0 ) then
#else
          if ( ldoc1d(n,i).gt.0.1 .and. sice1d(n,i).eq.0.0 ) then
#endif
            ii1 = ii1 + 1
          else if ( sice1d(n,i).gt.0.0 ) then
            ii2 = ii2 + 1
          else
            ii0 = ii0 + 1
          end if
        end do
        if ( ii0.ge.ii1 .and. ii0.ge.ii2 ) ioro(i) = 0
        if ( ii1.gt.ii0 .and. ii1.ge.ii2 ) ioro(i) = 1
        if ( ii2.gt.ii0 .and. ii2.gt.ii1 ) ioro(i) = 2
      end do
!
!     getdat() also sets calday (used in zenith() and radinp()).
!
      call getdat(jslc,h2ommr,alat,cld,clwp,coslat,loctim,o3mmr,o3vmr,  &
                & pilnm1,pintm1,pmidm1,pmlnm1,ps,qm1,tm1,ts)
!
!     NB:
!     variable coszrs is not calculated here in zenith()
!     but passed from vecbats() to colmod3() as an argument
!     therefore, subroutine zenith() is not necessary in regcm3
!
!     Get coszrs: needed as input to albland(), albocean(), radctl()
!
!KN   call zenith (calday  ,dodiavg ,alat    ,coszrs  )
!
!     NB:
!     land and ocean albedos are calculated
!     not in albland() and albocean() located below
!     but in albedov() called from subroutine vecbats()
!     therefore, the following subroutines are not called here
!
!     Find the albedo for land points
!
!KN   call albland(ilat     ,ioro    ,sndpth  ,coszrs  ,asdir   ,
!KN   $     aldir   ,asdif   ,aldif   )
!
!     Find the albedo for ocean/sea-ice points
!
!KN   call albocean(ilat     ,ioro    ,sndpth  ,coszrs  ,asdir   ,
!KN   $     aldir   ,asdif   ,aldif   )
!
!     NB:
!     albedos are copied from module bats2,
!     because variable names for albedos are somewhat different
!
      do i = 1 , iym1
        asdir(i) = aldirs(i)
        asdif(i) = aldifs(i)
        aldir(i) = aldirl(i)
        aldif(i) = aldifl(i)
      end do
!
!KN   modified above
!
!     Cloud particle size and fraction of ice
!
      call cldefr(ioro,tm1,rel,rei,fice,ps,pmidm1)
!
!     Cloud emissivity
!
      call cldems(clwp,fice,rei,emis)
!
!     Effective cloud cover
!
      do k = 1 , kz
        do i = 1 , iym1
          effcld(i,k) = cld(i,k)*emis(i,k)
        end do
      end do
!
!     Cloud cover at surface interface always zero (for safety's sake)
!
      do i = 1 , iym1
        effcld(i,kzp1) = 0.
        cld(i,kzp1) = 0.
      end do
!
!     Main radiation driving routine.
!     NB: All fluxes returned from radctl() have already been converted
!     to MKS.
!++csz
!add  by bixq
!add_
      call radctl(jslc,alat,coslat,ts,pmidm1,pintm1,pmlnm1,pilnm1,tm1,  &
                & qm1,cld,effcld,clwp,asdir,asdif,aldir,aldif,fsns,qrs, &
                & qrl,flwds,rel,rei,fice,sols,soll,solsd,solld,emiss1d, &
                & fsnt,fsntc,fsnsc,flnt,flns,flntc,flnsc,solin,alb,albc,&
                & fsds,fsnirt,fsnrtc,fsnirtsq,eccf,o3vmr)
                  ! input
!--csz
!
!     write out final results:
!
!KN   added below
!
!     subroutine radout() is not included in the ccm3 crm itself
!     but introduced from the former regcm2 radiation package
!     for the output of results from radiation calculations
!     this subroutine is used also for the coupling with bats
!
!     NB:
!     Names of some output variables (in MKS) have been changed
!     from those in the CCM2 radiation package.
!
!add  by bixq
!add_
      call radout(solin,fsnt,fsns,fsntc,fsnsc,qrs,flnt,flns,flntc,flnsc,&
                & qrl,flwds,srfrad,sols,soll,solsd,solld,alb,albc,fsds, &
                & fsnirt,fsnrtc,fsnirtsq,jslc,h2ommr,cld,clwp)
!
!KN   added above
!
      end subroutine colmod3
!
      subroutine cldefr(ioro,t,rel,rei,fice,ps,pmid)
!-----------------------------------------------------------------------
!
! Compute cloud drop size
!
!---------------------------Code history--------------------------------
!
! Original version:  J. Kiehl, January 1993
!
!-----------------------------------------------------------------------
!
      implicit none
!
!------------------------------Arguments--------------------------------
!
!     Input arguments
!
! ioro   - nint(oro(i))
! t      - Temperature
! ps     - surface pressure
! pmid   - midpoint pressures
!
!     Output arguments
!
! rel    - liquid effective drop size (microns)
! rei    - ice effective drop size (microns)
! fice   - fractional ice content within cloud
! pirnge - nrmlzd pres range for ice particle changes
! picemn - normalized pressure below which rei=reimax
! rirnge - range of ice radii (reimax - 10 microns)
! reimax - maximum ice effective radius
! pnrml  - normalized pressure
! weight - coef. for determining rei as fn of P/PS
!
!
! Dummy arguments
!
      real(8) , dimension(iym1,kz) :: fice , pmid , rei , rel , t
      integer , dimension(iym1) :: ioro
      real(8) , dimension(iym1) :: ps
      intent (in) ioro , pmid , ps , t
      intent (out) fice , rei , rel
!
!---------------------------Local workspace-----------------------------
!
! i, k   - longitude, level indices
! rliq   - temporary liquid drop size
!
!-----------------------------------------------------------------------
!
! Local variables
!
      integer :: i , k
      real(8) :: picemn , pirnge , pnrml , reimax , rirnge , rliq ,     &
               & weight
!
      do k = 1 , kz
        do i = 1 , iym1
!
!         Define liquid drop size
!
          if ( ioro(i).ne.1 ) then
!
!           Effective liquid radius over ocean and sea ice
!
            rliq = 10.0
          else
!
!           Effective liquid radius over land
!
            rliq = 5.0 + 5.0*dmin1(1.D0,dmax1(0.D0,(263.16-t(i,k))*0.05)&
                 & )
          end if
!
          rel(i,k) = rliq
!fil
!         test radius = 10.0
!         rel(i,k) = 10.0
!fil
!+        rei(i,k) = 30.0
!
!         Determine rei as function of normalized pressure
!
          reimax = 30.0
          rirnge = 20.0
          pirnge = 0.4
          picemn = 0.4
!
          pnrml = pmid(i,k)/ps(i)
          weight = dmax1(dmin1((pnrml-picemn)/pirnge,1.D0),0.D0)
          rei(i,k) = reimax - rirnge*weight
!
!         Define fractional amount of cloud that is ice
!
!         if warmer than -10 degrees C then water phase
!
          if ( t(i,k).gt.263.16 ) fice(i,k) = 0.0
!
!         if colder than -10 degrees C but warmer than -30 C mixed phase
!
          if ( t(i,k).le.263.16 .and. t(i,k).ge.243.16 ) fice(i,k)      &
             & = (263.16-t(i,k))/20.0
!
!         if colder than -30 degrees C then ice phase
!
          if ( t(i,k).lt.243.16 ) fice(i,k) = 1.0
!
!         Turn off ice radiative properties by setting fice = 0.0
!
!fil      no-ice test
!         fice(i,k) = 0.0
!
        end do
      end do
!
      end subroutine cldefr
!
      subroutine cldems(clwp,fice,rei,emis)

!-----------------------------------------------------------------------
!
! Compute cloud emissivity using cloud liquid water path (g/m**2)
!
!---------------------------Code history--------------------------------
!
! Original version:  J. Kiehl
! Standardized:      J. Rosinski, June 1992
! Reviewed:          J. Hack, J. Kiehl, August 1992
!
!-----------------------------------------------------------------------
!
      implicit none
!
! PARAMETER definitions
!
!     longwave absorption coeff (m**2/g)
      real(8) , parameter :: kabsl = 0.090361
!
!------------------------------Arguments--------------------------------
!
!     Input arguments
!
! clwp    - cloud liquid water path (g/m**2)
! rei     - ice effective drop size (microns)
! fice    - fractional ice content within cloud
!
!     Output arguments
!
! emis    - cloud emissivity (fraction)
!
!-----------------------------------------------------------------------
!
! Dummy arguments
!
      real(8) , dimension(iym1,kz) :: clwp , emis , fice , rei
      intent (in) clwp , fice , rei
      intent (out) emis
!
! Local variables
!
!-----------------------------------------------------------------------
!
! i, k    - longitude, level indices
! kabs    - longwave absorption coeff (m**2/g)
! kabsi   - ice absorption coefficient
!
!-----------------------------------------------------------------------
!
      integer :: i , k
      real(8) :: kabs , kabsi
!
      do k = 1 , kz
        do i = 1 , iym1
          kabsi = 0.005 + 1./rei(i,k)
          kabs = kabsl*(1.-fice(i,k)) + kabsi*fice(i,k)
          emis(i,k) = 1. - dexp(-1.66*kabs*clwp(i,k))
        end do
      end do
!
      end subroutine cldems
!
      end module mod_colmod3

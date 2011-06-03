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
!
module mod_colmod3
!
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
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  use mod_runparams
  use mod_bats
  use mod_date
  use mod_radiation
  use mod_main
  use mod_rad
  use mod_outrad
!
  private
!
  public :: colmod3
  real(8) , parameter :: lowcld = 1.0D-30
!
  contains
!
  subroutine colmod3(j)
!
    implicit none
!
    integer :: j
!
    real(8) , dimension(iym1) :: alb , albc , aldif , aldir , asdif , &
                                  asdir , alat , coslat , flns ,      &
                                  flnsc , flnt , flntc , flwds ,      &
                                  fsds , fsnirt , fsnirtsq , fsnrtc , &
                                  fsns , fsnsc , fsnt , fsntc ,       &
                                  loctim , solin , soll , solld ,     &
                                  sols , solsd , srfrad , ts
    real(8) , dimension(iym1,kzp1) :: cld , effcld , pilnm1 , pintm1
    real(8) , dimension(iym1,kz) :: clwp , emis , fice , h2ommr ,  &
           o3mmr , o3vmr , pmidm1 , pmlnm1 , qm1 , qrl , qrs , rei ,  &
           rel , tm1
    integer , dimension(iym1) :: ioro
    real(8) :: eccf
    integer :: i , ii0 , ii1 , ii2 , k , n
!
!   Fields specified by the user in getdat()
!
! land/ocean/sea ice flag
! Current latitude (radians)
! fractional cloud cover
! cloud liquid water path
!
!   NB: o3mmr and o3vmr should be dimensioned (iym1,kz) if a
!   different size radiation grid is used. Clashes between prgrid.h
!   and ptrrgrid.h (they both define plngbuf) prevent us from
!   dimensioning anything by kz in this top level crm() routine.
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
!   Fields computed from user input
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
!   Output longwave arguments from radctl()
!
! albedo: longwave, diffuse
! Surface down longwave flux
!
!   Output shortwave arguments from radctl()
!
! Longwave cooling rate
! Surface absorbed solar flux
! Solar heating rate
! Downward solar rad onto surface (lw direct)
! Downward solar rad onto surface (lw diffuse)
! Downward solar rad onto surface (sw direct)
!
!   Additional CRM diagnostic output from radctl()
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
!   Fundamental constants needed by radini()
!
! heat capacity dry air at constant prs (J/kg/K)
! ratio mean mol weight h2o to dry air
! gravitational acceleration (m/s**2)
! Sefan-Boltzmann constant (W/m**2/K**4)
!
!   Reset all arrays
!
    alb(:) = d_zero
    albc(:) = d_zero
    aldif(:) = d_zero
    aldir(:) = d_zero
    asdif(:) = d_zero
    asdir(:) = d_zero
    alat(:) = d_zero
    coslat(:) = d_zero
    flns(:) = d_zero
    flnsc(:) = d_zero
    flnt(:) = d_zero
    flntc(:) = d_zero
    flwds(:) = d_zero
    fsds(:) = d_zero
    fsnirt(:) = d_zero
    fsnirtsq(:) = d_zero
    fsnrtc(:) = d_zero
    fsns(:) = d_zero
    fsnsc(:) = d_zero
    fsnt(:) = d_zero
    fsntc(:) = d_zero
    loctim(:) = d_zero
    solin(:) = d_zero
    soll(:) = d_zero
    solld(:) = d_zero
    sols(:) = d_zero
    solsd(:) = d_zero
    srfrad(:) = d_zero
    ts(:) = d_zero
    cld(:,:) = d_zero
    effcld(:,:) = d_zero
    pilnm1(:,:) = d_zero
    pintm1(:,:) = d_zero
    clwp(:,:) = d_zero
    emis(:,:) = d_zero
    fice(:,:) = d_zero
    h2ommr(:,:) = d_zero
    o3mmr(:,:) = d_zero
    o3vmr(:,:) = d_zero
    pmidm1(:,:) = d_zero
    pmlnm1(:,:) = d_zero
    qm1(:,:) = d_zero
    qrl(:,:) = d_zero
    qrs(:,:) = d_zero
    rei(:,:) = d_zero
    rel(:,:) = d_zero
    tm1(:,:) = d_zero
    ioro(:) = 0
!
!   Set latitude index to j
!
    ilat = j
!
!   Set parameters in common block comtim : dosw,dolw,doabsems
!
    dosw = .true.
    dolw = .true.
    doabsems = .true.
!
!   radini sets many radiation parameters
!   radini() must be called before getdat(), because
!   the co2 mixing ratio set (by the user) in getdat() should
!   overwrite the default CCM3 co2 mixing ratio set by radini().
!
    call radini
!
!   NB: orography types are specified in the following
!
    do i = 1 , iym1
      ii0 = 0
      ii1 = 0
      ii2 = 0
      do n = 1 , nnsg
        if ( ocld2d(n,i,j) == 2 ) then
          ii2 = ii2 + 1
        else if ( ocld2d(n,i,j) == 1 ) then
          ii1 = ii1 + 1
        else
          ii0 = ii0 + 1
        end if
      end do
      if ( ii0 >= ii1 .and. ii0 >= ii2 ) ioro(i) = 0
      if ( ii1 > ii0 .and. ii1 >= ii2 ) ioro(i) = 1
      if ( ii2 > ii0 .and. ii2 > ii1 ) ioro(i) = 2
    end do
!
!   getdat() also sets calday (used in zenith() and radinp()).
!
    call getdat(j,h2ommr,alat,cld,clwp,coslat,loctim,o3mmr,o3vmr,  &
                pilnm1,pintm1,pmidm1,pmlnm1,ps,qm1,tm1,ts)
!
!   NB:
!   variable coszrs is not calculated here in zenith()
!   but passed from vecbats() to colmod3() as an argument
!   therefore, subroutine zenith() is not necessary in regcm3
!
!   NB:
!   land and ocean albedos are calculated
!   not in albland() and albocean() located below
!   but in albedov() called from subroutine vecbats()
!   therefore, the following subroutines are not called here
!
!   NB:
!   albedos are copied from module bats2,
!   because variable names for albedos are somewhat different
!
    do i = 1 , iym1
      asdir(i) = aldirs(i)
      asdif(i) = aldifs(i)
      aldir(i) = aldirl(i)
      aldif(i) = aldifl(i)
    end do
!
!   Cloud particle size and fraction of ice
!
    call cldefr(ioro,tm1,rel,rei,fice,ps,pmidm1)
!
!   Cloud emissivity
!
    call cldems(clwp,fice,rei,emis)
!
!   Effective cloud cover
!
    do k = 1 , kz
      do i = 1 , iym1
        effcld(i,k) = cld(i,k)*emis(i,k)
      end do
    end do
!
!   Cloud cover at surface interface always zero (for safety's sake)
!
    do i = 1 , iym1
      effcld(i,kzp1) = d_zero
      cld(i,kzp1) = d_zero
    end do
!
!   Main radiation driving routine.
!   NB: All fluxes returned from radctl() have already been converted
!   to MKS.
!
    call radctl(j,alat,coslat,ts,pmidm1,pintm1,pmlnm1,pilnm1,tm1,  &
                qm1,cld,effcld,clwp,asdir,asdif,aldir,aldif,fsns,qrs, &
                qrl,flwds,rel,rei,fice,sols,soll,solsd,solld,emiss1d, &
                fsnt,fsntc,fsnsc,flnt,flns,flntc,flnsc,solin,alb,albc,&
                fsds,fsnirt,fsnrtc,fsnirtsq,eccf,o3vmr)
!
!   subroutine radout() is not included in the ccm3 crm itself
!   but introduced from the former regcm2 radiation package
!   for the output of results from radiation calculations
!   this subroutine is used also for the coupling with bats
!
!   NB:
!   Names of some output variables (in MKS) have been changed
!   from those in the CCM2 radiation package.
!
    call radout(solin,fsnt,fsns,fsntc,fsnsc,qrs,flnt,flns,flntc,flnsc,&
                qrl,flwds,srfrad,sols,soll,solsd,solld,alb,albc,fsds, &
                fsnirt,fsnrtc,fsnirtsq,j,h2ommr,cld,clwp)
!
  end subroutine colmod3
!
!-----------------------------------------------------------------------
!
! Compute cloud drop size
!
!---------------------------Code history--------------------------------
!
! Original version:  J. Kiehl, January 1993
!
!------------------------------Arguments--------------------------------
!
!     Input arguments
!
! ioro   - idnint(oro(i))
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
!-----------------------------------------------------------------------
!
  subroutine cldefr(ioro,t,rel,rei,fice,ps,pmid)
!
    implicit none
!
    real(8) , dimension(iym1,kz) :: fice , pmid , rei , rel , t
    integer , dimension(iym1) :: ioro
    real(8) , dimension(iym1) :: ps
    intent (in) ioro , pmid , ps , t
    intent (out) fice , rei , rel
!
    integer :: i , k
    real(8) :: picemn , pirnge , pnrml , reimax , rirnge , rliq , weight
!
    do k = 1 , kz
      do i = 1 , iym1
!
!       Define liquid drop size
!
        if ( ioro(i) /= 1 ) then
!
!         Effective liquid radius over ocean and sea ice
!
          rliq = d_10
        else
!
!         Effective liquid radius over land
!
          rliq = d_five+d_five*dmin1(d_one,dmax1(d_zero,(263.16D0-t(i,k))*0.05D0))
        end if
!
        rel(i,k) = rliq
!fil
!       test radius = d_10
!       rel(i,k) = d_10
!fil
!+      rei(i,k) = 30.0
!
!       Determine rei as function of normalized pressure
!
        reimax = 30.0D0
        rirnge = 20.0D0
        pirnge = 0.4D0
        picemn = 0.4D0
!
        pnrml = pmid(i,k)/ps(i)
        weight = dmax1(dmin1((pnrml-picemn)/pirnge,d_one),d_zero)
        rei(i,k) = reimax - rirnge*weight
!
!       Define fractional amount of cloud that is ice
!
!       if warmer than -10 degrees C then water phase
!
        if ( t(i,k) > 263.16D0 ) fice(i,k) = d_zero
!
!       if colder than -10 degrees C but warmer than -30 C mixed phase
!
        if ( t(i,k) <= 263.16D0 .and. t(i,k) >= 243.16D0 ) &
          fice(i,k) = (263.16D0-t(i,k))/20.0D0
!
!       if colder than -30 degrees C then ice phase
!
        if ( t(i,k) < 243.16D0 ) fice(i,k) = d_one
!
!       Turn off ice radiative properties by setting fice = 0.0
!
!fil    no-ice test
!       fice(i,k) = d_zero
!
      end do
    end do
!
  end subroutine cldefr
!
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
  subroutine cldems(clwp,fice,rei,emis)
!
    implicit none
!
    real(8) , dimension(iym1,kz) :: clwp , emis , fice , rei
    intent (in) clwp , fice , rei
    intent (out) emis
!
!   longwave absorption coeff (m**2/g)
!
    real(8) , parameter :: kabsl = 0.090361D0
!
! i, k    - longitude, level indices
! kabs    - longwave absorption coeff (m**2/g)
! kabsi   - ice absorption coefficient
!
    integer :: i , k
    real(8) :: kabs , kabsi
!
    do k = 1 , kz
      do i = 1 , iym1
        kabsi = 0.005D0 + d_one/rei(i,k)
        kabs = kabsl*(d_one-fice(i,k)) + kabsi*fice(i,k)
        emis(i,k) = d_one - dexp(-1.66D0*kabs*clwp(i,k))
      end do
    end do
!
  end subroutine cldems
!
!-----------------------------------------------------------------------
!
! interface routine for column model that both initializes
! certain constants and reads external data:
!
! o3 mass mixing ratios are read in, but the model also requires the
! path lengths; they are computed here
!
! also, from the cloud input (fraction and liquid water path), the
! cloud longwave emissivity must be computed; this is done here
!
!     output arguments
!
! model latitude in radians
! cloud fraction
! cloud liquid water path (g/m**2)
! cosine latitude
! local time of solar computation
! o3 mass mixing ratio
! o3 volume mixing ratio
! ln(pintm1)
! pressure at model interfaces
! pressure at model mid-levels
! ln(pmidm1)
! model surface pressure field
! moisture field
! atmospheric temperature
! surface (air)  temperature
!
!-----------------------------------------------------------------------
!
  subroutine getdat(j,h2ommr,alat,cld,clwp,coslat,loctim,o3mmr,  &
                    o3vmr,pilnm1,pintm1,pmidm1,pmlnm1,ps,qm1,tm1,ts)
!
    implicit none
!
    integer :: j
    real(8) , dimension(iym1) :: alat , coslat , loctim , ps , ts
    real(8) , dimension(iym1,kzp1) :: cld , pilnm1 , pintm1
    real(8) , dimension(iym1,kz) :: clwp , h2ommr , o3mmr , o3vmr , &
                                    pmidm1 , pmlnm1 , qm1 , tm1
    intent (in) j
    intent (out) clwp , coslat , loctim , o3vmr , pilnm1 , pmlnm1 , qm1 , ts
    intent (inout) alat , cld , h2ommr , o3mmr , pintm1 , pmidm1 , ps , tm1
!
    real(8) :: ccvtem , clwtem , vmmr
    real(8) , dimension(iym1,kz) :: deltaz
    integer :: i , k , kj , ncldm1
    real(8) , dimension(iym1) :: rlat
!
    real(8) , parameter :: amd = 28.9644D0
    real(8) , parameter :: amo = 48.0000D0
!
!   begin read of data:
!
!   surface pressure and scaled pressure, from which level pressures are computed
    do i = 1 , iym1
      ps(i) = (sps2%ps(i,j)+r8pt)*d_10
      do k = 1 , kz
        pmidm1(i,k) = (sps2%ps(i,j)*a(k)+r8pt)*d_10
!KN     sclpr(k) = pmidm1(i,k)/ps(i)
      end do
    end do
!
!   convert pressures from mb to pascals and define interface pressures:
!
    do i = 1 , iym1
      ps(i) = ps(i)*d_100
      do k = 1 , kz
        pmidm1(i,k) = pmidm1(i,k)*d_100
        pmlnm1(i,k) = dlog(pmidm1(i,k))
      end do
    end do
    do k = 1 , kzp1
      do i = 1 , iym1
        pintm1(i,k) = (sps2%ps(i,j)*sigma(k)+r8pt)*d_1000
        pilnm1(i,k) = dlog(pintm1(i,k))
      end do
    end do
!
!   air temperatures
!
    do k = 1 , kz
      do i = 1 , iym1
        tm1(i,k) = atm2%t(i,k,j)/sps2%ps(i,j)
      end do
    end do
!
!   surface air temperature
!
!
!   h2o mass mixing ratio
!
    do k = 1 , kz
      do i = 1 , iym1
        h2ommr(i,k) = dmax1(1.0D-7,atm2%qv(i,k,j)/sps2%ps(i,j))
        qm1(i,k) = h2ommr(i,k)
      end do
    end do
!
!   o3 mass mixing ratio
!
    do k = 1 , kz
      do i = 1 , iym1
        kj = kzp1 - k
        o3mmr(i,k) = o3prof(i,kj,j)
      end do
    end do
!
!   fractional cloud cover (dependent on relative humidity)
!
!   qc   = gary's mods for clouds/radiation tie-in to exmois
    do k = 1 , kz
      do i = 1 , iym1
   
        ccvtem = d_zero   !cqc mod
!KN     cldfrc(i,k)=dmax1(cldfra(i,k)*0.9999999,ccvtem)
        cld(i,k) = dmax1(cldfra(i,k)*0.9999999D0,ccvtem)
!KN     cldfrc(i,k)=dmin1(cldfrc(i,k),0.9999999)
        cld(i,k) = dmin1(cld(i,k),0.9999999D0)
!
!       implement here the new formula then multiply by 10e6
!qc     if (tm1(i,k) > t0max) clwtem=clwmax
!qc     if (tm1(i,k) >= t0st .and. tm1(i,k) <= t0max) then
!qc       clwtem=clw0st+((tm1(i,k)-t0st)/(t0max-t0st))**2*(clwmax-clw0st)
!qc     end if
!qc     if (tm1(i,k) >= t0min .and. tm1(i,k) < t0st)
!qc       clwtem=clw0st+(tm1(i,k)-t0st)/(t0min-t0st)*(clwmin-clw0st)
!qc     end if
!qc     if (tm1(i,k) < t0min) clwtem=clwmin
!qc     clwtem=clwtem*1.e6
!
!       convert liquid water content into liquid water path, i.e.
!       multiply b deltaz
        clwtem = cldlwc(i,k) !cqc mod
        deltaz(i,k) = rgas*tm1(i,k)*(pintm1(i,k+1) - &
                        pintm1(i,k))/(egrav*pmidm1(i,k))
        clwp(i,k) = clwtem*deltaz(i,k)
!KN     if (cldfrc(i,k) == 0.) clwp(i,k)=d_zero
        if ( dabs(cld(i,k)) < lowcld ) then
          cld(i,k) = d_zero
          clwp(i,k) = d_zero
        end if
      end do
    end do
   
!   only allow thin clouds (<0.25) above 400 mb (yhuang, 11/97)
!   do k = 1 , kz
!     do i = 1 , iym1
!       if ( pintm1(i,k+1) < 40000.0D0 ) then
!         cld(i,k) = dmin1(cld(i,k),0.25d0)
!       else
!         cld(i,k)=dmin1(cld(i,k),0.7d0)
!       end if
!     end do
!   end do
!
!   set cloud fractional cover at top model level = 0
    do i = 1 , iym1
      cld(i,1) = d_zero
      clwp(i,1) = d_zero
      cld(i,2) = d_zero       !yhuang, 8/97 two-level
      clwp(i,2) = d_zero
    end do
!
!   set cloud fractional cover at bottom (ncld) model levels = 0
!
    ncldm1 = ncld - 1
    do k = kz - ncldm1 , kz
      do i = 1 , iym1
!KN     cldfrc(i,k)=d_zero
        cld(i,k) = d_zero
        clwp(i,k) = d_zero
      end do
    end do
!
!   ground temperature
!
    do i = 1 , iym1
!     tg(i)=tgb(i,jlsc)
!     when using bats calculate an equivalent ground (skin)
!     temperature by averaging over vegetated and non-vegetated areas
!jsp  tg(i)=((d_one-vgfrac(i))*tgb(i,j)**4.+vgfrac(i)*tlef2d(i,j)**4.)**0.25
!jsp  tg(i)=sfsta%tgbb(i,j)
      ts(i) = sfsta%tgbb(i,j)
    end do
!
!   cloud cover at surface interface always zero
!
    do i = 1 , iym1
!KN   effcld(i,kzp1) = d_zero
!KN   cldfrc(i,kzp1) = d_zero
      cld(i,kzp1) = d_zero
    end do
!
!KN adopted from regcm2 above
!
!----------------------------------------------------------------------
!
    do i = 1 , iym1
!
      do k = 1 , kz
        if ( cld(i,k) > 0.999D0 ) cld(i,k) = 0.999D0
      end do
!
      rlat(i) = mddom%xlat(i,j)
      calday = dble(julday) + dble(nnnnnn-nstrt0)*xdfbdy + &
                  (xtime/minph+gmt)/houpd
!
      loctim(i) = (calday-dint(calday))*houpd
      alat(i) = rlat(i)*degrad
      coslat(i) = dcos(alat(i))
!
    end do
!
!   Convert ozone mass mixing ratio to ozone volume mixing ratio:
!
    vmmr = amo/amd
    do k = 1 , kz
      do i = 1 , iym1
!       o3mmr(i,k) = vmmr*o3vmr(i,k)
        o3vmr(i,k) = o3mmr(i,k)/vmmr
      end do
    end do
!
  end subroutine getdat
!
end module mod_colmod3

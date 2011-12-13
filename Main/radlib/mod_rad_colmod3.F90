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
module mod_rad_colmod3
!
  use mod_realkinds
  use mod_dynparam
  use mod_rad_radiation
  use mod_rad_common
  use mod_rad_outrad
  use mod_rrtmg_driver
!
  private
!
  public :: allocate_mod_rad_colmod3 , colmod3
!
  real(dp) , parameter :: lowcld = 1.0D-30
!
!   longwave absorption coeff (m**2/g)
!
  real(dp) , parameter :: kabsl = 0.090361D0
!
  real(dp) , pointer , dimension(:) :: alb , albc , alat , ptrop ,    &
    flns , flnsc , flnt , flntc , flwds , fsds ,  fsnirt , fsnirtsq , &
    fsnrtc , fsns , fsnsc , fsnt , fsntc , solin , soll , solld ,     &
    sols , solsd , srfrad , ps , ts , emsvt1 , totcf
  real(dp) , pointer , dimension(:,:) :: cld , effcld , pilnm1 , pintm1
  real(dp) , pointer , dimension(:,:) :: clwp , emis , fice , h2ommr , &
    o3mmr , o3vmr , pmidm1 , pmlnm1 , qm1 , qrl , qrs , rei , rel ,    &
    deltaz , tm1
  integer , pointer , dimension(:) :: ioro
!
  contains
!
    subroutine allocate_mod_rad_colmod3
      implicit none
      call getmem1d(alb,1,jxp,'colmod3:alb')
      call getmem1d(albc,1,jxp,'colmod3:albc')
      call getmem1d(alat,1,jxp,'colmod3:alat')
      call getmem1d(ptrop,1,jxp,'colmod3:ptrop')
      call getmem1d(flns,1,jxp,'colmod3:flns')
      call getmem1d(flnsc,1,jxp,'colmod3:flnsc')
      call getmem1d(flnt,1,jxp,'colmod3:flnt')
      call getmem1d(flntc,1,jxp,'colmod3:flntc')
      call getmem1d(flwds,1,jxp,'colmod3:flwds')
      call getmem1d(fsds,1,jxp,'colmod3:fsds')
      call getmem1d(fsnirt,1,jxp,'colmod3:fsnirt')
      call getmem1d(fsnirtsq,1,jxp,'colmod3:fsnirtsq')
      call getmem1d(fsnrtc,1,jxp,'colmod3:fsnrtc')
      call getmem1d(fsns,1,jxp,'colmod3:fsns')
      call getmem1d(fsnsc,1,jxp,'colmod3:fsnsc')
      call getmem1d(fsnt,1,jxp,'colmod3:fsnt')
      call getmem1d(fsntc,1,jxp,'colmod3:fsntc')
      call getmem1d(solin,1,jxp,'colmod3:solin')
      call getmem1d(soll,1,jxp,'colmod3:soll')
      call getmem1d(solld,1,jxp,'colmod3:solld')
      call getmem1d(sols,1,jxp,'colmod3:sols')
      call getmem1d(solsd,1,jxp,'colmod3:solsd')
      call getmem1d(totcf,1,jxp,'colmod3:totcf')
      call getmem1d(srfrad,1,jxp,'colmod3:srfrad')
      call getmem1d(ps,1,jxp,'colmod3:ps')
      call getmem1d(ts,1,jxp,'colmod3:ts')
      call getmem1d(emsvt1,1,jxp,'colmod3:emsvt1')

      call getmem2d(cld,1,jxp,1,kzp1,'colmod3:cld')
      call getmem2d(effcld,1,jxp,1,kzp1,'colmod3:effcld')
      call getmem2d(pilnm1,1,jxp,1,kzp1,'colmod3:pilnm1')
      call getmem2d(pintm1,1,jxp,1,kzp1,'colmod3:pintm1')

      call getmem2d(clwp,1,jxp,1,kz,'colmod3:clwp')
      call getmem2d(emis,1,jxp,1,kz,'colmod3:emis')
      call getmem2d(fice,1,jxp,1,kz,'colmod3:fice')
      call getmem2d(h2ommr,1,jxp,1,kz,'colmod3:h2ommr')
      call getmem2d(o3mmr,1,jxp,1,kz,'colmod3:o3mmr')
      call getmem2d(o3vmr,1,jxp,1,kz,'colmod3:o3vmr')
      call getmem2d(pmidm1,1,jxp,1,kz,'colmod3:pmidm1')
      call getmem2d(pmlnm1,1,jxp,1,kz,'colmod3:pmlnm1')
      call getmem2d(qm1,1,jxp,1,kz,'colmod3:qm1')
      call getmem2d(qrl,1,jxp,1,kz,'colmod3:qrl')
      call getmem2d(qrs,1,jxp,1,kz,'colmod3:qrs')
      call getmem2d(rei,1,jxp,1,kz,'colmod3:rei')
      call getmem2d(rel,1,jxp,1,kz,'colmod3:rel')
      call getmem2d(tm1,1,jxp,1,kz,'colmod3:tm1')
      call getmem2d(deltaz,1,jxp,1,kz,'colmod3:deltaz')

      call getmem1d(ioro,1,jxp,'colmod3:ioro')
    end subroutine allocate_mod_rad_colmod3
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
  subroutine colmod3(jstart,jend,istart,iend,iyear,eccf,lout,labsem)
!
    implicit none
!
    integer , intent(in) :: iyear
    logical , intent(in) :: lout , labsem
    integer , intent(in) :: jstart , jend , istart , iend
    real(dp) , intent(in) :: eccf
!
    integer :: i , j , k , n , jj0 , jj1 , jj2
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
! solar incident flux
! Near-IR flux absorbed at toa
! Clear sky near-IR flux absorbed at toa
! Near-IR flux absorbed at toa >= 0.7 microns
! Flux Shortwave Downwelling Surface
!
  do i = istart , iend
!
!     Reset all arrays
!
      alb(:) = d_zero
      albc(:) = d_zero
      alat(:) = d_zero
      ptrop(:) = d_zero
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
!     Set parameters dosw , dolw , doabsems
!
      dosw = .true.
      dolw = .true.
      doabsems = .true.
!
!     radini sets many radiation parameters
!     radini() must be called before getdat(), because
!     the co2 mixing ratio set (by the user) in getdat() should
!     overwrite the default CCM3 co2 mixing ratio set by radini().
!
      call radini(iyear)
!
!     NB: orography types are specified in the following
!
      do j = jstart , jend
        jj0 = 0
        jj1 = 0
        jj2 = 0
        do n = 1 , nnsg
          if ( lndocnicemsk(n,j,i) == 2 ) then
            jj2 = jj2 + 1
          else if ( lndocnicemsk(n,j,i) == 1 ) then
            jj1 = jj1 + 1
          else
            jj0 = jj0 + 1
          end if
        end do
        if ( jj0 >= jj1 .and. jj0 >= jj2 ) ioro(j) = 0
        if ( jj1 > jj0 .and. jj1 >= jj2 ) ioro(j) = 1
        if ( jj2 > jj0 .and. jj2 > jj1 ) ioro(j) = 2
      end do
!
      call getdat(jstart,jend,i)
!
!     NB:
!     land and ocean albedos are calculated
!     not in albland() and albocean() located below
!     but in albedov() called from subroutine vecbats()
!     therefore, the following subroutines are not called here
!
!     Cloud particle size and fraction of ice
!
      call cldefr(jstart,jend)
!
!     Cloud emissivity
!
      call cldems(jstart,jend)
!
!     Effective cloud cover
!
      do k = 1 , kz
        do j = jstart , jend
          effcld(j,k) = cld(j,k)*emis(j,k)
        end do
      end do
!
!     Cloud cover at surface interface always zero (for safety's sake)
!
      do j = jstart , jend
        effcld(j,kzp1) = d_zero
        cld(j,kzp1) = d_zero
      end do
!
!     Main radiation driving routine.
!     NB: All fluxes returned from radctl() have already been converted
!     to MKS.
!
      call radctl(jstart,jend,i,alat,ptrop,ts,pmidm1,pintm1,pmlnm1,       &
                  pilnm1,tm1,qm1,cld,effcld,clwp,fsns,qrs,qrl,flwds,rel,  &
                  rei,fice,sols,soll,solsd,solld,emsvt1,fsnt,fsntc,fsnsc, &
                  flnt,flns,flntc,flnsc,solin,alb,albc,fsds,fsnirt,       &
                  fsnrtc,fsnirtsq,totcf,eccf,o3vmr,labsem)
!
!     subroutine radout() is not included in the ccm3 crm itself
!     but introduced from the former regcm radiation package
!     for the output of results from radiation calculations
!     this subroutine is used also for the coupling with bats
!
!     NB:
!     Names of some output variables (in MKS) have been changed
!     from those in the CCM2 radiation package.
!
      call radout(jstart,jend,i,lout,solin,fsnt,fsns,fsntc,fsnsc,qrs, &
                  flnt,flns,flntc,flnsc,qrl,flwds,srfrad,sols,soll,   &
                  solsd,solld,alb,albc,fsds,fsnirt,fsnrtc,fsnirtsq,   &
                  totcf,h2ommr,cld,clwp)
    end do
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
! ioro   - idnint(oro(j))
! t      - Temperature
! ps     - surface pressure
! pmid   - midpoint pressures
!
!     Output arguments
!
! rel    - liquid effective drop size (microns)
! rei    - ice effective drop size (microns)
! fice   - fractional ice content within cloud
! pnrml  - normalized pressure
! weight - coef. for determining rei as fn of P/PS
!
!-----------------------------------------------------------------------
!
  subroutine cldefr(jstart,jend)
    implicit none
    integer , intent(in) :: jstart , jend
!
    integer :: j , k
    real(dp) :: pnrml , rliq , weight
!
!   reimax - maximum ice effective radius
    real(dp) , parameter :: reimax = 30.0D0
!   rirnge - range of ice radii (reimax - 10 microns)
    real(dp) , parameter :: rirnge = 20.0D0
!   pirnge - nrmlzd pres range for ice particle changes
    real(dp) , parameter :: pirnge = 0.4D0
!   picemn - normalized pressure below which rei=reimax
    real(dp) , parameter :: picemn = 0.4D0
!   Temperatures in K (263.16 , 243.16)
    real(dp) , parameter :: minus10 = wattp-d_10
    real(dp) , parameter :: minus30 = wattp-(d_three*d_10)
!
    do k = 1 , kz
      do j = jstart , jend
!
!       Define liquid drop size
!
        if ( ioro(j) /= 1 ) then
!
!         Effective liquid radius over ocean and sea ice
!
          rliq = d_10
        else
!
!         Effective liquid radius over land
!
          rliq = d_five+d_five* & 
                  dmin1(d_one,dmax1(d_zero,(minus10-tm1(j,k))*0.05D0))
        end if
!
        rel(j,k) = rliq
!fil
!       test radius = d_10
!       rel(j,k) = d_10
!fil
!+      rei(j,k) = 30.0
!
!       Determine rei as function of normalized pressure
!
        pnrml = pmidm1(j,k)/ps(j)
        weight = dmax1(dmin1((pnrml-picemn)/pirnge,d_one),d_zero)
        rei(j,k) = reimax - rirnge*weight
!
!       Define fractional amount of cloud that is ice
!
!       if warmer than -10 degrees C then water phase
!
        if ( tm1(j,k) > minus10 ) fice(j,k) = d_zero
!
!       if colder than -10 degrees C but warmer than -30 C mixed phase
!
        if ( tm1(j,k) <= minus10 .and. tm1(j,k) >= minus30 ) &
          fice(j,k) = (minus10-tm1(j,k))/20.0D0
!
!       if colder than -30 degrees C then ice phase
!
        if ( tm1(j,k) < minus30 ) fice(j,k) = d_one
!
!       Turn off ice radiative properties by setting fice = 0.0
!
!fil    no-ice test
!       fice(j,k) = d_zero
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
  subroutine cldems(jstart,jend)
    implicit none
    integer , intent(in) :: jstart , jend
!
!   kabs    - longwave absorption coeff (m**2/g)
!   kabsi   - ice absorption coefficient
!
    integer :: j , k
    real(dp) :: kabs , kabsi
!
    do k = 1 , kz
      do j = jstart , jend
        kabsi = 0.005D0 + d_one/rei(j,k)
        kabs = kabsl*(d_one-fice(j,k)) + kabsi*fice(j,k)
        emis(j,k) = d_one - dexp(-1.66D0*kabs*clwp(j,k))
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
  subroutine getdat(jstart,jend,i)
!
    implicit none
!
    integer , intent(in) :: jstart , jend , i
!
    integer :: j , k , kj , ncldm1
    real(dp) :: ccvtem , clwtem , vmmr
!
    real(dp) , parameter :: amd = 28.9644D0
    real(dp) , parameter :: amo = 48.0000D0
!
    if ( iemiss == 1 ) then
      do j = jstart , jend
        emsvt1(j) = emsvt(j,i)
      end do
    end if
!
!   surface pressure and scaled pressure, from which level are computed
!
    do j = jstart , jend
      ps(j) = (sfps(j,i)+ptp)*d_10
      do k = 1 , kz
        pmidm1(j,k) = (sfps(j,i)*hlev(k)+ptp)*d_10
      end do
    end do
!
!   convert pressures from mb to pascals and define interface pressures:
!
    do j = jstart , jend
      ps(j) = ps(j)*d_100
      do k = 1 , kz
        pmidm1(j,k) = pmidm1(j,k)*d_100
        pmlnm1(j,k) = dlog(pmidm1(j,k))
      end do
    end do
    do k = 1 , kzp1
      do j = jstart , jend
        pintm1(j,k) = (sfps(j,i)*flev(k)+ptp)*d_1000
        pilnm1(j,k) = dlog(pintm1(j,k))
      end do
    end do
!
!   air temperatures
!
    do k = 1 , kz
      do j = jstart , jend
        tm1(j,k) = tatms(j,i,k)
      end do
    end do
!
!   h2o mass mixing ratio
!
    do k = 1 , kz
      do j = jstart , jend
        h2ommr(j,k) = dmax1(1.0D-7,qvatms(j,i,k))
        qm1(j,k) = h2ommr(j,k)
      end do
    end do
!
!   o3 mass mixing ratio
!
    do k = 1 , kz
      do j = jstart , jend
        kj = kzp1 - k
        o3mmr(j,k) = o3prof(j,i,kj)
      end do
    end do
!
!   fractional cloud cover (dependent on relative humidity)
!
    do k = 1 , kz
      do j = jstart , jend
   
        ccvtem = d_zero   !cqc mod
        cld(j,k) = dmax1(cldfra(j,i,k)*0.9999999D0,ccvtem)
        cld(j,k) = dmin1(cld(j,k),0.9999999D0)
!
!qc     gary's mods for clouds/radiation tie-in to exmois
!qc     implement here the new formula then multiply by 10e6
!qc     if (tm1(j,k) > t0max) clwtem=clwmax
!qc     if (tm1(j,k) >= t0st .and. tm1(j,k) <= t0max) then
!qc       clwtem=clw0st+((tm1(j,k)-t0st)/(t0max-t0st))**2*(clwmax-clw0st)
!qc     end if
!qc     if (tm1(j,k) >= t0min .and. tm1(j,k) < t0st)
!qc       clwtem=clw0st+(tm1(j,k)-t0st)/(t0min-t0st)*(clwmin-clw0st)
!qc     end if
!qc     if (tm1(j,k) < t0min) clwtem=clwmin
!qc     clwtem=clwtem*1.e6
!
!       convert liquid water content into liquid water path, i.e.
!       multiply b deltaz
!
        clwtem = cldlwc(j,i,k) !cqc mod
        deltaz(j,k) = rgas*tm1(j,k)*(pintm1(j,k+1) - &
                        pintm1(j,k))/(egrav*pmidm1(j,k))
        clwp(j,k) = clwtem*deltaz(j,k)
        if ( dabs(cld(j,k)) < lowcld ) then
          cld(j,k) = d_zero
          clwp(j,k) = d_zero
        end if
      end do
    end do
!   
!   only allow thin clouds (<0.25) above 400 mb (yhuang, 11/97)
!   do k = 1 , kz
!     do j = jstart , jend
!       if ( pintm1(j,k+1) < 40000.0D0 ) then
!         cld(j,k) = dmin1(cld(j,k),0.25d0)
!       else
!         cld(j,k)=dmin1(cld(j,k),0.7d0)
!       end if
!     end do
!   end do
!
!   set cloud fractional cover at top model level = 0
    do j = jstart , jend
      cld(j,1) = d_zero
      cld(j,2) = d_zero       !yhuang, 8/97 two-level
      clwp(j,1) = d_zero
      clwp(j,2) = d_zero
    end do
!
!   set cloud fractional cover at bottom (ncld) model levels = 0
!
    ncldm1 = ncld - 1
    do k = kz - ncldm1 , kz
      do j = jstart , jend
        cld(j,k) = d_zero
        clwp(j,k) = d_zero
      end do
    end do
!
!   ground temperature
!
    do j = jstart , jend
!     tg(j)=tgb(i,j)
!     when using bats calculate an equivalent ground (skin)
!     temperature by averaging over vegetated and non-vegetated areas
!jsp  tg(j)=((d_one-vgfrac(j))*tgb(i,j)**4.+vgfrac(j)*tlef2d(i,j)**4.)**0.25
      ts(j) = tground(j,i)
    end do
!
!   cloud cover at surface interface always zero
!
    do j = jstart , jend
      cld(j,kzp1) = d_zero
    end do
!
!----------------------------------------------------------------------
!
    do j = jstart , jend
!
      do k = 1 , kz
        if ( cld(j,k) > 0.999D0 ) cld(j,k) = 0.999D0
      end do
!
      alat(j) = xlat(j,i)*degrad
      ! pressure of tropopause
      ptrop(j) = 250.0D2 - 150.0D2*dcos(alat(j))**d_two
!
    end do
!
!   Convert ozone mass mixing ratio to ozone volume mixing ratio:
!
    vmmr = amo/amd
    do k = 1 , kz
      do j = jstart , jend
        o3vmr(j,k) = o3mmr(j,k)/vmmr
      end do
    end do
!
  end subroutine getdat
!
end module mod_rad_colmod3

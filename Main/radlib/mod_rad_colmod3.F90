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
  use mod_service
  use mod_rad_radiation
  use mod_rad_common
  use mod_rad_outrad
  use mod_rrtmg_driver
  use mod_rad_aerosol
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
  real(dp) , parameter :: nearone  = 0.99D+00
!
  real(dp) , pointer , dimension(:) :: alb , albc , &
    flns , flnsc , flnt , flntc , flwds , fsds ,  fsnirt , fsnirtsq ,  &
    fsnrtc , fsns , fsnsc , fsnt , fsntc , solin , soll , solld ,      &
    sols , solsd , ps , ts , emsvt1 , totcf , totcl , totci , xptrop , &
    dlat , czen
  real(dp) , pointer , dimension(:) :: aeradfo , aeradfos
  real(dp) , pointer , dimension(:) :: aerlwfo , aerlwfos
  real(dp) , pointer , dimension(:) :: adirsw , adifsw , adirlw , adiflw
  real(dp) , pointer , dimension(:) :: asw , alw
  real(dp) , pointer , dimension(:) :: abv , sol
  real(dp) , pointer , dimension(:,:) :: cld , effcld , pilnm1 , pintm1
  real(dp) , pointer , dimension(:,:) :: clwp , emis , fice , h2ommr , &
    o3mmr , o3vmr , pmidm1 , pmlnm1 , qm1 , qrl , qrs , rei , rel ,    &
    deltaz , tm1 , rh1
  real(dp) , pointer , dimension(:,:,:) :: tauxcl , tauxci
  real(dp) , pointer , dimension(:,:,:) :: aermmr
  real(dp) , pointer , dimension(:,:,:) :: absgasnxt
  real(dp) , pointer , dimension(:,:,:) :: absgastot
  real(dp) , pointer , dimension(:,:) :: emsgastot
  logical , pointer , dimension(:) :: czengt0
  integer , pointer , dimension(:) :: ioro
!
  integer :: npr

  contains
!
    subroutine allocate_mod_rad_colmod3(ichem)
      implicit none
      integer , intent(in) :: ichem
      npr = (jci2-jci1+1)*(ici2-ici1+1)
      call getmem1d(alb,1,npr,'colmod3:alb')
      call getmem1d(albc,1,npr,'colmod3:albc')
      call getmem1d(flns,1,npr,'colmod3:flns')
      call getmem1d(flnsc,1,npr,'colmod3:flnsc')
      call getmem1d(flnt,1,npr,'colmod3:flnt')
      call getmem1d(flntc,1,npr,'colmod3:flntc')
      call getmem1d(flwds,1,npr,'colmod3:flwds')
      call getmem1d(fsds,1,npr,'colmod3:fsds')
      call getmem1d(fsnirt,1,npr,'colmod3:fsnirt')
      call getmem1d(fsnirtsq,1,npr,'colmod3:fsnirtsq')
      call getmem1d(fsnrtc,1,npr,'colmod3:fsnrtc')
      call getmem1d(fsns,1,npr,'colmod3:fsns')
      call getmem1d(fsnsc,1,npr,'colmod3:fsnsc')
      call getmem1d(fsnt,1,npr,'colmod3:fsnt')
      call getmem1d(fsntc,1,npr,'colmod3:fsntc')
      call getmem1d(solin,1,npr,'colmod3:solin')
      call getmem1d(soll,1,npr,'colmod3:soll')
      call getmem1d(solld,1,npr,'colmod3:solld')
      call getmem1d(sols,1,npr,'colmod3:sols')
      call getmem1d(solsd,1,npr,'colmod3:solsd')
      call getmem1d(totcf,1,npr,'colmod3:totcf')
      call getmem1d(totcl,1,npr,'colmod3:totcl')
      call getmem1d(totci,1,npr,'colmod3:totci')
      call getmem1d(ps,1,npr,'colmod3:ps')
      call getmem1d(ts,1,npr,'colmod3:ts')
      call getmem1d(emsvt1,1,npr,'colmod3:emsvt1')
      call getmem1d(xptrop,1,npr,'rad:xptrop')
      call getmem1d(dlat,1,npr,'rad:dlat')
      call getmem1d(adirsw,1,npr,'rad:adirsw')
      call getmem1d(adifsw,1,npr,'rad:adifsw')
      call getmem1d(adirlw,1,npr,'rad:adirlw')
      call getmem1d(adiflw,1,npr,'rad:adiflw')
      call getmem1d(asw,1,npr,'rad:asw')
      call getmem1d(alw,1,npr,'rad:alw')
      call getmem1d(abv,1,npr,'rad:abv')
      call getmem1d(sol,1,npr,'rad:sol')

      call getmem2d(cld,1,npr,1,kzp1,'colmod3:cld')
      call getmem2d(effcld,1,npr,1,kzp1,'colmod3:effcld')
      call getmem2d(pilnm1,1,npr,1,kzp1,'colmod3:pilnm1')
      call getmem2d(pintm1,1,npr,1,kzp1,'colmod3:pintm1')

      call getmem2d(rh1,1,npr,1,kz,'colmod3:rh1')
      call getmem2d(clwp,1,npr,1,kz,'colmod3:clwp')
      call getmem2d(emis,1,npr,1,kz,'colmod3:emis')
      call getmem2d(fice,1,npr,1,kz,'colmod3:fice')
      call getmem2d(h2ommr,1,npr,1,kz,'colmod3:h2ommr')
      call getmem2d(o3mmr,1,npr,1,kz,'colmod3:o3mmr')
      call getmem2d(o3vmr,1,npr,1,kz,'colmod3:o3vmr')
      call getmem2d(pmidm1,1,npr,1,kz,'colmod3:pmidm1')
      call getmem2d(pmlnm1,1,npr,1,kz,'colmod3:pmlnm1')
      call getmem2d(qm1,1,npr,1,kz,'colmod3:qm1')
      call getmem2d(qrl,1,npr,1,kz,'colmod3:qrl')
      call getmem2d(qrs,1,npr,1,kz,'colmod3:qrs')
      call getmem2d(rei,1,npr,1,kz,'colmod3:rei')
      call getmem2d(rel,1,npr,1,kz,'colmod3:rel')
      call getmem2d(tm1,1,npr,1,kz,'colmod3:tm1')
      call getmem2d(deltaz,1,npr,1,kz,'colmod3:deltaz')
      call getmem1d(aeradfo,1,npr,'colmod3:aeradfo')
      call getmem1d(aeradfos,1,npr,'colmod3:aeradfos')
      call getmem1d(aerlwfo,1,npr,'colmod3:aerlwfo')
      call getmem1d(aerlwfos,1,npr,'colmod3:aerlwfos')
      call getmem1d(czen,1,npr,'colmod3:czen')
      call getmem1d(czengt0,1,npr,'colmod3:czengt0')
      call getmem3d(absgasnxt,1,npr,1,kz,1,4,'colmod3:absgasnxt')
      call getmem3d(absgastot,1,npr,1,kzp1,1,kzp1,'colmod3:absgastot')
      call getmem2d(emsgastot,1,npr,1,kzp1,'colmod3:emsgastot')
      call getmem3d(tauxcl,1,npr,0,kz,1,nspi,'colmod3:tauxcl')
      call getmem3d(tauxci,1,npr,0,kz,1,nspi,'colmod3:tauxci')

      call getmem1d(ioro,1,npr,'colmod3:ioro')

      if ( ichem == 1 ) then
        call getmem3d(aermmr,1,npr,1,kz,1,ntr,'colmod3:aermmr')
      end if

      dosw = .true.
      dolw = .true.
      doabsems = .true.

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
  subroutine colmod3(iyear,eccf,lout,labsem)
!
    implicit none
!
    integer , intent(in) :: iyear
    logical , intent(in) :: lout , labsem
    real(dp) , intent(in) :: eccf
!
    integer :: i , j , n , jj0 , jj1 , jj2
    character (len=64) :: subroutine_name='colmod3'
    integer :: indx = 0
!
    call time_begin(subroutine_name,indx)
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
!   Reset all arrays
!
    czen(:) = d_zero
    czengt0(:) = .false.
    alb(:) = d_zero
    albc(:) = d_zero
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
    ioro(:) = -1
!
!   radini sets many radiation parameters
!   radini() must be called before getdat(), because
!   the co2 mixing ratio set (by the user) in getdat() should
!   overwrite the default CCM3 co2 mixing ratio set by radini().
!
    call radini(iyear)
!
!   NB: orography types are specified in the following
!
    n = 1
    do i = ici1 , ici2
      do j = jci1 , jci2
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
        if ( jj0 >= jj1 .and. jj0 >= jj2 ) ioro(n) = 0
        if ( jj1 > jj0 .and. jj1 >= jj2 ) ioro(n) = 1
        if ( jj2 > jj0 .and. jj2 > jj1 ) ioro(n) = 2
        n = n + 1
      end do
    end do
!
    call getdat
!
!   Cloud particle size and fraction of ice
!
    call cldefr
!
!   Cloud emissivity
!
    call cldems
!
!   Effective cloud cover
!
    effcld(:,1:kz) = cld(:,1:kz)*emis(:,1:kz)
!
!   Cloud cover at surface interface always zero (for safety's sake)
!
    effcld(:,kzp1) = d_zero
    cld(:,kzp1) = d_zero
!
!   Main radiation driving routine.
!   NB: All fluxes returned from radctl() have already been converted
!   to MKS.
!
    call radctl(1,npr,dlat,xptrop,ts,pmidm1,pintm1,pmlnm1,pilnm1,       &
                tm1,qm1,rh1,cld,effcld,clwp,aermmr,fsns,qrs,qrl,flwds,  &
                rel,rei,fice,sols,soll,solsd,solld,emsvt1,fsnt,fsntc,   &
                fsnsc,flnt,flns,flntc,flnsc,solin,alb,albc,fsds,fsnirt, &
                fsnrtc,fsnirtsq,totcf,eccf,o3vmr,czen,czengt0,adirsw,   &
                adifsw,adirlw,adiflw,asw,alw,abv,sol,aeradfo,aeradfos,  &
                aerlwfo,aerlwfos,absgasnxt,absgastot,emsgastot,tauxcl,  &
                tauxci,labsem)
!
!   Save gas emission/absorbtion
!
    n = 1
    do i = ici1 , ici2
      do j = jci1 , jci2
        gasabsnxt(j,i,:,:) = absgasnxt(n,:,:)
        gasabstot(j,i,:,:) = absgastot(n,:,:)
        gasemstot(j,i,:) = emsgastot(n,:)
        taucldsp(j,i,:,:) = tauxcl(n,:,:) + tauxci(n,:,:)
        n = n + 1
      end do
    end do
!
!   subroutine radout() is not included in the ccm3 crm itself
!   but introduced from the former regcm radiation package
!   for the output of results from radiation calculations
!   this subroutine is used also for the coupling with bats
!
!   NB:
!     Names of some output variables (in MKS) have been changed
!     from those in the CCM2 radiation package.
!
    call radout(1,npr,lout,solin,fsnt,fsns,fsntc,fsnsc,qrs,flnt,flns, &
                flntc,flnsc,qrl,flwds,sols,soll,solsd,solld,alb,albc, &
                fsds,fsnirt,fsnrtc,fsnirtsq,totcf,totcl,totci,h2ommr, &
                cld,clwp,abv,sol,aeradfo,aeradfos,aerlwfo,aerlwfos,   &
                tauxar3d,tauasc3d,gtota3d)
!
    call time_end(subroutine_name,indx)
  end subroutine colmod3
!
!-----------------------------------------------------------------------
!
! Compute cloud drop size
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
  subroutine cldefr
    implicit none
!
    integer :: n , k
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
    character (len=64) :: subroutine_name='cldefr'
    integer :: indx = 0
!
    call time_begin(subroutine_name,indx)
!
    totci(:) = d_zero
    do k = 1 , kz
      do n = 1 , npr
!
!       Define liquid drop size
!
        if ( ioro(n) /= 1 ) then
!
!         Effective liquid radius over ocean and sea ice
!
          rliq = d_10
        else
!
!         Effective liquid radius over land
!
          rliq = d_five+d_five* & 
                  dmin1(d_one,dmax1(d_zero,(minus10-tm1(n,k))*0.05D0))
        end if
!
        rel(n,k) = rliq
!fil
!       test radius = d_10
!       rel(n,k) = d_10
!fil
!+      rei(n,k) = 30.0
!
!       Determine rei as function of normalized pressure
!
        pnrml = pmidm1(n,k)/ps(n)
        weight = dmax1(dmin1((pnrml-picemn)/pirnge,d_one),d_zero)
        rei(n,k) = reimax - rirnge*weight
!
!       Define fractional amount of cloud that is ice
!
!       if warmer than -10 degrees C then water phase
!
        if ( tm1(n,k) > minus10 ) fice(n,k) = d_zero
!
!       if colder than -10 degrees C but warmer than -30 C mixed phase
!
        if ( tm1(n,k) <= minus10 .and. tm1(n,k) >= minus30 ) &
          fice(n,k) = (minus10-tm1(n,k))/20.0D0
!
!       if colder than -30 degrees C then ice phase
!
        if ( tm1(n,k) < minus30 ) fice(n,k) = d_one
!
!       Turn off ice radiative properties by setting fice = 0.0
!
!fil    no-ice test
!       fice(n,k) = d_zero
        totci(n) = totci(n) + clwp(n,k)*fice(n,k)*d_r1000
!
      end do
    end do
    call time_end(subroutine_name,indx)
!
  end subroutine cldefr
!
!-----------------------------------------------------------------------
!
! Compute cloud emissivity using cloud liquid water path (g/m**2)
!
!     Output arguments
!
! emis    - cloud emissivity (fraction)
!
!-----------------------------------------------------------------------
!
  subroutine cldems
    implicit none
!
!   kabs    - longwave absorption coeff (m**2/g)
!   kabsi   - ice absorption coefficient
!
    integer :: n , k
    real(dp) :: kabs , kabsi
    character (len=64) :: subroutine_name='cldems'
    integer :: indx = 0
!
    call time_begin(subroutine_name,indx)
!
    do k = 1 , kz
      do n = 1 , npr
        kabsi = 0.005D0 + d_one/rei(n,k)
        kabs = kabsl*(d_one-fice(n,k)) + kabsi*fice(n,k)
        emis(n,k) = d_one - dexp(-1.66D0*kabs*clwp(n,k))
      end do
    end do
!
    call time_end(subroutine_name,indx)
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
  subroutine getdat
!
    implicit none
!
    integer :: n , i , j , k , k2 , itr , krev , ncldm1
    real(dp) :: ccvtem , clwtem
!
    real(dp) , parameter :: amd = 28.9644D0
    real(dp) , parameter :: amo = 48.0000D0
    real(dp) , parameter :: vmmr = amo/amd
    logical , save :: ifirst
!
    character (len=64) :: subroutine_name='getdat'
    integer :: indx = 0
!
    data ifirst /.true./
!
    call time_begin(subroutine_name,indx)
!
!   Static informations
!
    if ( ifirst ) then
      dlat = reshape(dabs(xlat(jci1:jci2,ici1:ici2)),(/npr/))
      xptrop = reshape(ptrop(jci1:jci2,ici1:ici2),(/npr/))
      ifirst = .false.
    end if
!
    if ( iemiss == 1 ) then
      emsvt1 = reshape(emsvt(jci1:jci2,ici1:ici2),(/npr/))
    end if
!
!   Albedoes
!
    adirsw = reshape(swdiralb(jci1:jci2,ici1:ici2),(/npr/))
    adifsw = reshape(swdifalb(jci1:jci2,ici1:ici2),(/npr/))
    adirlw = reshape(lwdiralb(jci1:jci2,ici1:ici2),(/npr/))
    adiflw = reshape(lwdifalb(jci1:jci2,ici1:ici2),(/npr/))
    asw = reshape(swalb(jci1:jci2,ici1:ici2),(/npr/))
    alw = reshape(lwalb(jci1:jci2,ici1:ici2),(/npr/))
!
!   Sun elevation
!
    czen = reshape(coszen(jci1:jci2,ici1:ici2),(/npr/))
    where ( czen > d_zero )
      czengt0 = .true.
    end where

    do n = 1 , 4
      do k = 1 , kz
        absgasnxt(:,k,n) = reshape(gasabsnxt(jci1:jci2,ici1:ici2,k,n),(/npr/))
      end do
    end do
    do k = 1 , kzp1
      do k2 = 1 , kzp1
        absgastot(:,k2,k) = reshape(gasabstot(jci1:jci2,ici1:ici2,k2,k),(/npr/))
      end do
    end do
    do k = 1 , kzp1
      emsgastot(:,k) = reshape(gasemstot(jci1:jci2,ici1:ici2,k),(/npr/))
    end do
!
!   surface pressure and scaled pressure, from which level are computed
!
    ps = reshape(((sfps(jci1:jci2,ici1:ici2)+ptp)*d_1000),(/npr/))
!
!   convert pressures from mb to pascals and define interface pressures:
!
    do k = 1 , kz
      pmidm1(:,k) = reshape(((sfps(jci1:jci2,ici1:ici2)*hlev(k)+ptp)*d_1000), &
                            (/npr/))
    end do
    pmlnm1(:,:) = dlog(pmidm1(:,:))

    do k = 1 , kzp1
      pintm1(:,k) = reshape(((sfps(jci1:jci2,ici1:ici2)*flev(k)+ptp)*d_1000), &
                            (/npr/))
    end do
    pilnm1(:,:) = dlog(pintm1(:,:))
!
!   air temperature and relative humidity
!
    do k = 1 , kz
      tm1(:,k) = reshape(tatms(jci1:jci2,ici1:ici2,k),(/npr/))
      rh1(:,k) = reshape( &
         dmax1(dmin1(rhatms(jci1:jci2,ici1:ici2,k),nearone),d_zero),(/npr/))
    end do
!
!   h2o mass mixing ratio
!
    do k = 1 , kz
      h2ommr(:,k) = reshape(dmax1(1.0D-8,qvatms(jci1:jci2,ici1:ici2,k)),(/npr/))
    end do
    qm1(:,:) = h2ommr(:,:)
!
!   fractional cloud cover (dependent on relative humidity)
!
    totcl(:) = d_zero
    do k = 1 , kz
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
   
          ccvtem = d_zero   !cqc mod
          cld(n,k) = dmax1(cldfra(j,i,k)*0.9999999D0,ccvtem)
          cld(n,k) = dmin1(cld(n,k),0.9999999D0)
!
!qc       gary's mods for clouds/radiation tie-in to exmois
!qc       implement here the new formula then multiply by 10e6
!qc       if (tm1(n,k) > t0max) clwtem=clwmax
!qc       if (tm1(n,k) >= t0st .and. tm1(n,k) <= t0max) then
!qc         clwtem=clw0st+((tm1(n,k)-t0st)/(t0max-t0st))**2*(clwmax-clw0st)
!qc       end if
!qc       if (tm1(n,k) >= t0min .and. tm1(n,k) < t0st)
!qc         clwtem=clw0st+(tm1(n,k)-t0st)/(t0min-t0st)*(clwmin-clw0st)
!qc       end if
!qc       if (tm1(n,k) < t0min) clwtem=clwmin
!qc       clwtem=clwtem*1.e6
!
!         convert liquid water content into liquid water path, i.e.
!         multiply b deltaz
!
          clwtem = cldlwc(j,i,k) !cqc mod
          deltaz(n,k) = rgas*tm1(n,k)*(pintm1(n,k+1) - &
                          pintm1(n,k))/(egrav*pmidm1(n,k))
          clwp(n,k) = clwtem*deltaz(n,k)
          if ( dabs(cld(n,k)) < lowcld ) then
            cld(n,k) = d_zero
            clwp(n,k) = d_zero
          end if
          totcl(n) = totcl(n)+clwp(n,k)*d_r1000
          n = n + 1
        end do
      end do
    end do
!   
!   only allow thin clouds (<0.25) above 400 mb (yhuang, 11/97)
!   do k = 1 , kz
!     do n = 1 , npr
!       if ( pintm1(n,k+1) < 40000.0D0 ) then
!         cld(n,k) = dmin1(cld(n,k),0.25d0)
!       else
!         cld(n,k)=dmin1(cld(n,k),0.7d0)
!       end if
!     end do
!   end do
!
!   set cloud fractional cover at top model level = 0
!
    cld(:,1:2) = d_zero
    clwp(:,1:2) = d_zero
!
!   set cloud fractional cover at bottom (ncld) model levels = 0
!
    ncldm1 = ncld - 1
    do k = kz - ncldm1 , kz
      do n = 1 , npr
        cld(n,k) = d_zero
        clwp(n,k) = d_zero
      end do
    end do
!
!   cloud cover (at surface interface always zero)
!
    cld(:,kzp1) = d_zero
!
!   Give correct range if out of range
!
    do k = 1 , kz
      do n = 1 , npr
        if ( cld(n,k) > 0.999D0 ) cld(n,k) = 0.999D0
      end do
    end do
!
!   ground temperature
!
    ts = reshape(tground(jci1:jci2,ici1:ici2),(/npr/))
!
!   o3 mass and volume mixing ratios
!
    do k = 1 , kz
      krev = kzp1 - k
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          o3mmr(n,k) = o3prof(j,i,krev)
          n = n + 1
        end do
      end do
    end do
    o3vmr(:,:) = o3mmr(:,:)/vmmr
!
    if ( lchem ) then
      do itr = 1 , ntr
        do k = 1 , kz
          n = 1
          do i = ici1 , ici2
            do j = jci1 , jci2
              aermmr(n,k,itr) = chspmix(j,i,k,itr)/psfps(j,i)
              n = n + 1
            end do
          end do
        end do
      end do
    end if

    call time_end(subroutine_name,indx)
  end subroutine getdat
!
end module mod_rad_colmod3

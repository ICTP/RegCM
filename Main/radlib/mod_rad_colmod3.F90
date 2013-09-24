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
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_service
  use mod_runparams , only : iemiss , iqv
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
! Allowed range for cloud fraction
!
  real(rk8) , parameter :: lowcld = 0.00001D0
  real(rk8) , parameter :: hicld  = 0.99999D0
!
!   longwave absorption coeff (m**2/g)
!
  real(rk8) , parameter :: kabsl = 0.090361D0
  real(rk8) , parameter :: nearone  = 0.99D+00
!
  real(rk8) , pointer , dimension(:) :: alb , albc , &
    flns , flnsc , flnt , flntc , flwds , fsds ,  fsnirt , fsnirtsq ,  &
    fsnrtc , fsns , fsnsc , fsnt , fsntc , solin , soll , solld ,      &
    sols , solsd , ps , ts , emsvt1 , totcf , totcl , totci , xptrop , &
    dlat , czen
  real(rk8) , pointer , dimension(:) :: aeradfo , aeradfos
  real(rk8) , pointer , dimension(:) :: aerlwfo , aerlwfos
  real(rk8) , pointer , dimension(:) :: adirsw , adifsw , adirlw , adiflw
  real(rk8) , pointer , dimension(:) :: asw , alw
  real(rk8) , pointer , dimension(:) :: abv , sol
  real(rk8) , pointer , dimension(:,:) :: cld , effcld , pilnm1 , pintm1
  real(rk8) , pointer , dimension(:,:) :: clwp , emis , fice , h2ommr , &
    o3mmr , o3vmr , pmidm1 , pmlnm1 , qm1 , qrl , qrs , rei , rel ,    &
    deltaz , tm1 , rh1
  real(rk8) , pointer , dimension(:,:,:) :: tauxcl , tauxci
  real(rk8) , pointer , dimension(:,:,:) :: aermmr
  real(rk8) , pointer , dimension(:,:,:) :: absgasnxt
  real(rk8) , pointer , dimension(:,:,:) :: absgastot
  real(rk8) , pointer , dimension(:,:) :: emsgastot
  real(rk8) , pointer , dimension(:,:,:) :: outtaucl , outtauci
  logical , pointer , dimension(:) :: czengt0
  integer(ik4) , pointer , dimension(:) :: ioro
!
  integer(ik4) :: npr

  contains
!
    subroutine allocate_mod_rad_colmod3
      implicit none
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
      call getmem3d(outtaucl,1,npr,1,kzp1,1,4,'colmod3:outtaucl')
      call getmem3d(outtauci,1,npr,1,kzp1,1,4,'colmod3:outtauci')

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
    integer(ik4) , intent(in) :: iyear
    logical , intent(in) :: lout , labsem
    real(rk8) , intent(in) :: eccf
!
    integer(ik4) :: i , j , k , n , m , k2 , jj0 , jj1 , jj2
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'colmod3'
    integer(ik4) :: indx = 0
    call time_begin(subroutine_name,indx)
#endif
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
    ! radini sets many radiation parameters
    !
    call radini(iyear)
    !
    ! NB: orography types are specified in the following
    !
    m = 1
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
        if ( jj0 >= jj1 .and. jj0 >= jj2 ) ioro(m) = 0
        if ( jj1 >  jj0 .and. jj1 >= jj2 ) ioro(m) = 1
        if ( jj2 >  jj0 .and. jj2 >  jj1 ) ioro(m) = 2
        m = m + 1
      end do
    end do
    !
    ! Copy data from atmosphere
    !
    call getdat
    !
    ! Cloud particle size and fraction of ice
    !
    call cldefr
    !
    ! Cloud emissivity
    !
    call cldems
    !
    ! Effective cloud cover
    !
    effcld(:,1:kz) = cld(:,1:kz)*emis(:,1:kz)
    !
    ! Main radiation driving routine.
    ! NB: All fluxes returned from radctl() have already been converted to MKS.
    !
    call radctl(1,npr,dlat,xptrop,ts,pmidm1,pintm1,pmlnm1,pilnm1,       &
                tm1,qm1,rh1,cld,effcld,clwp,aermmr,fsns,qrs,qrl,flwds,  &
                rel,rei,fice,sols,soll,solsd,solld,emsvt1,fsnt,fsntc,   &
                fsnsc,flnt,flns,flntc,flnsc,solin,alb,albc,fsds,fsnirt, &
                fsnrtc,fsnirtsq,totcf,eccf,o3vmr,czen,czengt0,adirsw,   &
                adifsw,adirlw,adiflw,asw,alw,abv,sol,aeradfo,aeradfos,  &
                aerlwfo,aerlwfos,absgasnxt,absgastot,emsgastot,tauxcl,  &
                tauxci,outtaucl,outtauci,labsem)
    !
    ! Save gas emission/absorbtion
    !
    if ( labsem ) then
      do m = 1 , 4
        do k = 1 , kz
          n = 1
          do i = ici1 , ici2
            do j = jci1 , jci2
              gasabsnxt(j,i,k,m) = absgasnxt(n,k,m)
              n = n + 1
            end do
          end do
        end do
      end do
      do k = 1 , kzp1
        do k2 = 1 , kzp1
          n = 1
          do i = ici1 , ici2
            do j = jci1 , jci2
              gasabstot(j,i,k2,k) = absgastot(n,k2,k)
              n = n + 1
            end do
          end do
        end do
      end do
      do k = 1 , kzp1
        n = 1
        do i = ici1 , ici2
          do j = jci1 , jci2
            gasemstot(j,i,k) = emsgastot(n,k)
            n = n + 1
          end do
        end do
      end do
    end if
    if ( ichem == 1 ) then
      do m = 1 , nspi
        do k = 0 , kz
          n = 1
          do i = ici1 , ici2
            do j = jci1 , jci2
              taucldsp(j,i,k,m)  = tauxcl(n,k,m) + tauxci(n,k,m)
              n = n + 1
            end do
          end do
        end do
      end do
    end if
    !
    ! subroutine radout() copies back the data to RegCM for surface
    ! computations and output purposes.
    !
    call radout(lout,solin,fsnt,fsns,fsntc,fsnsc,qrs,flnt,flns,  &
                flntc,flnsc,qrl,flwds,sols,soll,solsd,solld,     &
                totcf,totcl,totci,cld,clwp,abv,sol,aeradfo,      &
                aeradfos,aerlwfo,aerlwfos,tauxar3d,tauasc3d,     &
                gtota3d,deltaz,outtaucl,outtauci)
#ifdef DEBUG
    call time_end(subroutine_name,indx)
#endif
  end subroutine colmod3
!
!-----------------------------------------------------------------------
!
! Compute cloud drop size
!
!-----------------------------------------------------------------------
!
  subroutine cldefr
    implicit none
!
    integer(ik4) :: n , k , nt
    real(rk8) :: pnrml , rliq , weight , rhoa , nc , aerc , lwc , kparam
    ! reimax - maximum ice effective radius
    real(rk8) , parameter :: reimax = 30.0D0
    ! rirnge - range of ice radii (reimax - 10 microns)
    real(rk8) , parameter :: rirnge = 20.0D0
    ! pirnge - nrmlzd pres range for ice particle changes
    real(rk8) , parameter :: pirnge = 0.4D0
    ! picemn - normalized pressure below which rei=reimax
    real(rk8) , parameter :: picemn = 0.4D0
    ! Temperatures in K (263.16 , 243.16)
    real(rk8) , parameter :: minus10 = wattp-d_10
    real(rk8) , parameter :: minus30 = wattp-(d_three*d_10)
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'cldefr'
    integer(ik4) :: indx = 0
    call time_begin(subroutine_name,indx)
#endif
!
    totci(:) = d_zero
    do k = 1 , kz
      do n = 1 , npr
        ! Define liquid drop size
        if ( ioro(n) /= 1 ) then
          ! Effective liquid radius over ocean and sea ice
          rliq = d_10
        else
          ! Effective liquid radius over land
          rliq = d_five+d_five* & 
                   dmin1(d_one,dmax1(d_zero,(minus10-tm1(n,k))*0.05D0))
        end if
        ! rel : liquid effective drop size (microns)
        rel(n,k) = rliq
!fil
!       test radius
!       rel(n,k) = d_10
!       rei(n,k) = 30.0
!fil
        ! Determine rei as function of normalized pressure
        pnrml = pmidm1(n,k)/ps(n)
        ! weight coef. for determining rei as fn of P/PS
        weight = dmax1(dmin1((pnrml-picemn)/pirnge,d_one),d_zero)
        ! rei : ice effective drop size (microns)
        rei(n,k) = reimax - rirnge*weight
        ! Define fractional amount of cloud that is ice
        ! if warmer than -10 degrees C then water phase
        if ( tm1(n,k) > minus10 ) then
          fice(n,k) = d_zero
        else if ( tm1(n,k) <= minus10 .and. tm1(n,k) >= minus30 ) then
        ! if colder than -10 degrees C but warmer than -30 C mixed phase
        ! fice : fractional ice content within cloud
          fice(n,k) = (minus10-tm1(n,k))/20.0D0
        !  if colder than -30 degrees C then ice phase
        else
          fice(n,k) = d_one
        end if
        !  Turn off ice radiative properties by setting fice = 0.0
!
!fil    no-ice test
!       fice(n,k) = d_zero
!fil
        totci(n) = totci(n) + clwp(n,k)*fice(n,k)*d_r1000
      end do
    end do

    !FAB : reintroduce simple sulfate indirect effect 
    ! from Qian  1999
    ! clwp is passed  in g/m2 
    if ( ichem == 1 .and. iindirect == 1 ) then
      do nt = 1 , ntr 
        if ( chtrname(nt) /= 'SO4' ) cycle
        do k = 1 , kz
          do n = 1 , npr
            rhoa = pmidm1(n,k)* amdk / (tm1(n,k) * 8.314D0)
            aerc = rhoa * aermmr(n,k,nt) * 1.D9 !microg/m3
            ! thershold of 0.1 microg/m3 for activation of indirect effect
            if ( aerc > 0.1D0 ) then
              ! ccn number concentration in cm-3 
              nc = 90.7D0 * aerc**0.45D0 + 23.D0
              ! kg/m3, already account fro cum and ls clouds      
              lwc = clwp(n,k) / deltaz(n,k) * d_r1000
              if ( lwc < 1.D-6 ) cycle
              if ( ioro(n) /= 1 ) then
                ! Martin et al.(1994) parameter over ocean and sea ice
                kparam = 0.80D0
              else
               ! and over land
               kparam = 0.67D0
              end if
              !finally modify effective radius 
              !(1.D6 to convert to rel to microm, 1D6 to vonvert nc in m-3)
              rel(n,k) = 1.D6 * ( d_three*lwc / &
                   (d_four*mathpi*rhoh2o*kparam*nc*1.D6 ) )**(d_one/d_three)
            end if
          end do
        end do
      end do
    end if
#ifdef DEBUG
    call time_end(subroutine_name,indx)
#endif
  end subroutine cldefr
!
!-----------------------------------------------------------------------
!
! Compute cloud emissivity using cloud liquid water path (g/m**2)
!
!-----------------------------------------------------------------------
!
  subroutine cldems
    implicit none
    integer(ik4) :: n , k
    real(rk8) :: kabs , kabsi
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'cldems'
    integer(ik4) :: indx = 0
    call time_begin(subroutine_name,indx)
#endif
    do k = 1 , kz
      do n = 1 , npr
        ! ice absorption coefficient
        kabsi = 0.005D0 + d_one/rei(n,k)
        ! longwave absorption coeff (m**2/g)
        kabs = kabsl*(d_one-fice(n,k)) + kabsi*fice(n,k)
        ! cloud emissivity (fraction)
        emis(n,k) = d_one - dexp(-1.66D0*kabs*clwp(n,k))
      end do
    end do
#ifdef DEBUG
    call time_end(subroutine_name,indx)
#endif
  end subroutine cldems
!
!-----------------------------------------------------------------------
!
! Interface routine for column model that initializes internal variables
! A copy if them is performed on local 1D arrays.
!
!-----------------------------------------------------------------------
!
  subroutine getdat
    implicit none
!
    integer(ik4) :: n , m , i , j , k , k2 , itr , krev , kmincld , kmaxcld
    real(rk8) , parameter :: amd = 28.9644D0
    real(rk8) , parameter :: amo = 48.0000D0
    real(rk8) , parameter :: vmmr = amo/amd
    logical , save :: ifirst = .true.
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'getdat'
    integer(ik4) :: indx = 0
    call time_begin(subroutine_name,indx)
#endif
    !
    ! Static informations
    !
    if ( ifirst ) then
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          dlat(n) = dabs(xlat(j,i))
          xptrop(n) = ptrop(j,i)
          n = n + 1
        end do
      end do
      if ( iemiss == 1 ) then
        n = 1
        do i = ici1 , ici2
          do j = jci1 , jci2
            emsvt1(n) = emsvt(j,i)
            n = n + 1
          end do
        end do
      else
        emsvt1(:) = 0.9995D0
      end if
      ifirst = .false.
    end if
    !
    ! Albedoes
    !
    n = 1
    do i = ici1 , ici2
      do j = jci1 , jci2
        adirsw(n) = swdiralb(j,i)
        adifsw(n) = swdifalb(j,i)
        adirlw(n) = lwdiralb(j,i)
        adiflw(n) = lwdifalb(j,i)
        asw(n)    = swalb(j,i)
        alw(n)    = lwalb(j,i)
        n = n + 1
      end do
    end do
    !
    ! Sun elevation
    !
    n = 1
    do i = ici1 , ici2
      do j = jci1 , jci2
        czen(n) = coszen(j,i)
        if ( czen(n) > d_zero ) czengt0(n) = .true.
        n = n + 1
      end do
    end do
    !
    ! Gas concentrations
    !
    do m = 1 , 4
      do k = 1 , kz
        n = 1
        do i = ici1 , ici2
          do j = jci1 , jci2
            absgasnxt(n,k,m) = gasabsnxt(j,i,k,m)
            n = n + 1
          end do
        end do
      end do
    end do
    do k = 1 , kzp1
      do k2 = 1 , kzp1
        n = 1
        do i = ici1 , ici2
          do j = jci1 , jci2
            absgastot(n,k2,k) = gasabstot(j,i,k2,k)
            n = n + 1
          end do
        end do
      end do
    end do
    do k = 1 , kzp1
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          emsgastot(n,k) = gasemstot(j,i,k)
          n = n + 1
        end do
      end do
    end do
    !
    ! Surface pressure and scaled pressure, from which level are computed
    !
    n = 1
    do i = ici1 , ici2
      do j = jci1 , jci2
        ps(n) = (sfps(j,i)+ptop)*d_1000
        n = n + 1
      end do
    end do
    !
    ! convert pressures from mb to pascals and define interface pressures:
    !
    do k = 1 , kz
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          pmidm1(n,k) = (sfps(j,i)*hsigma(k)+ptop)*d_1000
          n = n + 1
        end do
      end do
    end do
    pmlnm1(:,:) = dlog(pmidm1(:,:))
    do k = 1 , kzp1
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          pintm1(n,k) = (sfps(j,i)*sigma(k)+ptop)*d_1000
          n = n + 1
        end do
      end do
    end do
    pilnm1(:,:) = dlog(pintm1(:,:))
    !
    ! Air temperature and relative humidity
    !
    do k = 1 , kz
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          tm1(n,k) = tatms(j,i,k)
          rh1(n,k) = dmax1(dmin1(rhatms(j,i,k),nearone),d_zero)
          n = n + 1
        end do
      end do
    end do
    !
    ! H2O mass mixing ratio
    !
    do k = 1 , kz
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          h2ommr(n,k) = dmax1(minqx,qxatms(j,i,k,iqv))
          n = n + 1
        end do
      end do
    end do
    qm1(:,:) = h2ommr(:,:)
    !
    ! deltaz
    !
    do k = 1 , kz
      do n = 1 , npr
        deltaz(n,k) = rgas*tm1(n,k)*(pintm1(n,k+1) - &
                      pintm1(n,k))/(egrav*pmidm1(n,k))
      end do
    end do
    !
    ! Fractional cloud cover (dependent on relative humidity)
    ! Set cloud
    !   - NOT on the topmost two layers
    !   - Starting from ncld levels from the surface
    !
    totcl(:) = d_zero
    kmaxcld = 3
    kmincld = kz-ncld
    do k = kmaxcld , kmincld
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          cld(n,k) = dmax1(cldfra(j,i,k),lowcld)
          cld(n,k) = dmin1(cld(n,k),hicld)
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
          !
          ! Convert liquid water content into liquid water path, i.e.
          ! multiply b deltaz
          !
          clwp(n,k) = cldlwc(j,i,k)*deltaz(n,k)
          if ( cldfra(j,i,k) < lowcld ) then
            clwp(n,k) = d_zero
          end if
          totcl(n) = totcl(n)+clwp(n,k)*d_r1000
          n = n + 1
        end do
      end do
    end do
!
!   only allow thin clouds (<0.25) above 400 mb (yhuang, 11/97)
!
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
    !
    ! Ground temperature
    !
    n = 1
    do i = ici1 , ici2
      do j = jci1 , jci2
        ts(n) = tground(j,i)
        n = n + 1
      end do
    end do
    !
    ! O3 mass and volume mixing ratios
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
    ! Tracers mixing ratios
    !
    if ( ichem == 1 ) then
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
#ifdef DEBUG
    call time_end(subroutine_name,indx)
#endif
  end subroutine getdat
!
end module mod_rad_colmod3

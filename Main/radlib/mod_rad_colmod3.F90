!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_rad_colmod3
  use mod_intkinds, only : ik4
  use mod_realkinds, only : rkx
  use mod_constants, only : amd, amo3, mathpi
  use mod_constants, only : rgas, rhoh2o, tzero
  use mod_dynparam, only : jci1, jci2, ici1, ici2, kz, kzp1
  use mod_dynparam, only : ntr, nspi, myid
  use mod_memutil, only : getmem, getmem, getmem
  use mod_runparams, only : ipptls, eccf, ichem, idirect, iindirect
  use mod_runparams, only : iqv, iqc, iqi, ncld, chtrname
  use mod_runparams, only : cftotmax, iclimaaer
  use mod_rad_radiation, only : radctl, radtype
  use mod_rad_common, only : gasabsnxt, gasabstot, gasemstot
  use mod_rad_common, only : o3prof, taucldsp, dosw, dolw, doabsems
  use mod_rad_outrad, only : radout
  use mod_rad_aerosol, only : aermmr, tauxar3d, tauasc3d, gtota3d
  use mod_regcm_types, only : mod_2_rad, rad_2_mod
#ifdef DEBUG
  use mod_service, only : time_begin, time_end, dbgslen
#endif

  implicit none

  private

  public :: allocate_mod_rad_colmod3, colmod3

  type(radtype) :: rt

  contains

  subroutine allocate_mod_rad_colmod3
    implicit none
    integer(ik4) :: npr
    npr = (jci2-jci1+1)*(ici2-ici1+1)
    rt%n1 = 1
    rt%n2 = npr
    call getmem(rt%ioro,1,npr,'colmod3:ioro')
    call getmem(rt%dlat,1,npr,'rad:dlat')
    call getmem(rt%xptrop,1,npr,'rad:xptrop')
    call getmem(rt%ts,1,npr,'colmod3:ts')
    call getmem(rt%ps,1,npr,'colmod3:ps')
    call getmem(rt%totcl,1,npr,'colmod3:totcl')
    call getmem(rt%totci,1,npr,'colmod3:totci')
    call getmem(rt%totwv,1,npr,'colmod3:totwv')
    call getmem(rt%fsns,1,npr,'colmod3:fsns')
    call getmem(rt%flwds,1,npr,'colmod3:flwds')
    call getmem(rt%sols,1,npr,'colmod3:sols')
    call getmem(rt%soll,1,npr,'colmod3:soll')
    call getmem(rt%solsd,1,npr,'colmod3:solsd')
    call getmem(rt%solld,1,npr,'colmod3:solld')
    call getmem(rt%emiss,1,npr,'colmod3:emiss')
    call getmem(rt%fsnt,1,npr,'colmod3:fsnt')
    call getmem(rt%fsntc,1,npr,'colmod3:fsntc')
    call getmem(rt%fsnsc,1,npr,'colmod3:fsnsc')
    call getmem(rt%flnt,1,npr,'colmod3:flnt')
    call getmem(rt%lwout,1,npr,'colmod3:lwout')
    call getmem(rt%lwin,1,npr,'colmod3:lwin')
    call getmem(rt%flns,1,npr,'colmod3:flns')
    call getmem(rt%flntc,1,npr,'colmod3:flntc')
    call getmem(rt%flnsc,1,npr,'colmod3:flnsc')
    call getmem(rt%solin,1,npr,'colmod3:solin')
    call getmem(rt%solout,1,npr,'colmod3:solout')
    call getmem(rt%alb,1,npr,'colmod3:alb')
    call getmem(rt%albc,1,npr,'colmod3:albc')
    call getmem(rt%fsds,1,npr,'colmod3:fsds')
    call getmem(rt%fsnirt,1,npr,'colmod3:fsnirt')
    call getmem(rt%fsnrtc,1,npr,'colmod3:fsnrtc')
    call getmem(rt%fsnirtsq,1,npr,'colmod3:fsnirtsq')
    call getmem(rt%totcf,1,npr,'colmod3:totcf')
    call getmem(rt%czen,1,npr,'colmod3:czen')
    call getmem(rt%czengt0,1,npr,'colmod3:czengt0')
    call getmem(rt%adirsw,1,npr,'rad:adirsw')
    call getmem(rt%adifsw,1,npr,'rad:adifsw')
    call getmem(rt%adirlw,1,npr,'rad:adirlw')
    call getmem(rt%adiflw,1,npr,'rad:adiflw')
    call getmem(rt%asw,1,npr,'rad:asw')
    call getmem(rt%alw,1,npr,'rad:alw')
    call getmem(rt%abv,1,npr,'rad:abv')
    call getmem(rt%sol,1,npr,'rad:sol')
    call getmem(rt%aeradfo,1,npr,'colmod3:aeradfo')
    call getmem(rt%aeradfos,1,npr,'colmod3:aeradfos')
    call getmem(rt%aerlwfo,1,npr,'colmod3:aerlwfo')
    call getmem(rt%aerlwfos,1,npr,'colmod3:aerlwfos')
    ! Only for RRTM
    rt%asaeradfo => null( )
    rt%asaeradfos => null( )
    rt%asaerlwfo => null( )
    rt%asaerlwfos => null( )
    call getmem(rt%pmid,1,kz,1,npr,'colmod3:pmid')
    call getmem(rt%pint,1,kzp1,1,npr,'colmod3:pint')
    call getmem(rt%pmln,1,kz,1,npr,'colmod3:pmln')
    call getmem(rt%piln,1,kzp1,1,npr,'colmod3:piln')
    call getmem(rt%t,1,kz,1,npr,'colmod3:t')
    call getmem(rt%q,1,kz,1,npr,'colmod3:q')
    call getmem(rt%ql,1,kz,1,npr,'colmod3:ql')
    call getmem(rt%qi,1,kz,1,npr,'colmod3:qi')
    call getmem(rt%dz,1,kz,1,npr,'colmod3:dz')
    call getmem(rt%rh,1,kz,1,npr,'colmod3:rh')
    call getmem(rt%rho,1,kz,1,npr,'colmod3:rho')
    call getmem(rt%cld,1,kzp1,1,npr,'colmod3:cld')
    call getmem(rt%effcld,1,kzp1,1,npr,'colmod3:effcld')
    call getmem(rt%clwp,1,kz,1,npr,'colmod3:clwp')
    call getmem(rt%qrs,1,kz,1,npr,'colmod3:qrs')
    call getmem(rt%qrl,1,kz,1,npr,'colmod3:qrl')
    call getmem(rt%rel,1,kz,1,npr,'colmod3:rel')
    call getmem(rt%rei,1,kz,1,npr,'colmod3:rei')
    call getmem(rt%fice,1,kz,1,npr,'colmod3:fice')
    call getmem(rt%o3vmr,1,kz,1,npr,'colmod3:o3vmr')
    call getmem(rt%emsgastot,1,kzp1,1,npr,'colmod3:emsgastot')
    call getmem(rt%absgasnxt,1,kz,1,4,1,npr,'colmod3:absgasnxt')
    call getmem(rt%absgastot,1,kzp1,1,kzp1,1,npr,'colmod3:absgastot')
    call getmem(rt%tauxcl,0,kz,1,npr,1,nspi,'colmod3:tauxcl')
    call getmem(rt%tauxci,0,kz,1,npr,1,nspi,'colmod3:tauxci')
    call getmem(rt%outtaucl,1,kzp1,1,4,1,npr,'colmod3:outtaucl')
    call getmem(rt%outtauci,1,kzp1,1,4,1,npr,'colmod3:outtauci')
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
  subroutine colmod3(iyear,imonth,lout,labsem,m2r,r2m)
    !@acc use nvtx
    implicit none
    type(mod_2_rad), intent(in) :: m2r
    type(rad_2_mod), intent(inout) :: r2m
    integer(ik4), intent(in) :: iyear, imonth
    logical, intent(in) :: lout, labsem

    integer(ik4) :: n, m, i, j, k, k2, itr
    integer(ik4) :: ni, nj, npr
    real(rkx) :: totci_l, totcl_l, totwv_l
    real(rkx) :: nc, aerc, lwc, kparam
    real(rkx) :: kabs, kabsi, kabsl, cldemis, arg
    real(rkx) :: iwc, tempc, tcels, fsr, aiwc, biwc, desr
    real(rkx), dimension(kz,rt%n1:rt%n2) :: temp
    !real(rkx) :: tpara
    real(rkx), parameter :: minus20 = 253.15_rkx
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
    !   NB: o3vmr should be dimensioned (iym1,kz) if a
    !   different size radiation grid is used. Clashes between prgrid.h
    !   and ptrrgrid.h (they both define plngbuf) prevent us from
    !   dimensioning anything by kz in this top level crm() routine.
    !
    ! cosine latitude
    ! water vapor mass mixing ratio
    ! Ozone mass mixing ratio
    ! Ozone volume mixing ratio
    ! natural log of pint
    ! model interface pressures
    ! model level pressures
    ! natural log of pmid
    ! model level specific humidity
    ! model level temperatures
    !
    !   Fields computed from user input
    !
    ! surface air temperature
    ! cloud emissivity
    ! effective cloud = cloud fraction * cloud emissivity
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
    !@acc call nvtxStartRange("colmod3")
#if 0
    !$acc kernels
    rt%alb(:) = 0.0_rkx
    rt%albc(:) = 0.0_rkx
    rt%flns(:) = 0.0_rkx
    rt%flnsc(:) = 0.0_rkx
    rt%flnt(:) = 0.0_rkx
    rt%lwout(:) = 0.0_rkx
    rt%lwin(:) = 0.0_rkx
    rt%flntc(:) = 0.0_rkx
    rt%flwds(:) = 0.0_rkx
    rt%fsds(:) = 0.0_rkx
    rt%fsnirt(:) = 0.0_rkx
    rt%fsnirtsq(:) = 0.0_rkx
    rt%fsnrtc(:) = 0.0_rkx
    rt%fsns(:) = 0.0_rkx
    rt%fsnsc(:) = 0.0_rkx
    rt%fsnt(:) = 0.0_rkx
    rt%fsntc(:) = 0.0_rkx
    rt%solin(:) = 0.0_rkx
    rt%solout(:) = 0.0_rkx
    rt%soll(:) = 0.0_rkx
    rt%solld(:) = 0.0_rkx
    rt%sols(:) = 0.0_rkx
    rt%solsd(:) = 0.0_rkx
    rt%qrl(:,:) = 0.0_rkx
    rt%qrs(:,:) = 0.0_rkx
    !$acc end kernels
#else
    npr = (jci2-jci1+1)*(ici2-ici1+1)
    do concurrent ( n = 1:npr )
     rt%alb(n) = 0.0_rkx
     rt%albc(n) = 0.0_rkx
     rt%flns(n) = 0.0_rkx
     rt%flnsc(n) = 0.0_rkx
     rt%flnt(n) = 0.0_rkx
     rt%lwout(n) = 0.0_rkx
     rt%lwin(n) = 0.0_rkx
     rt%flntc(n) = 0.0_rkx
     rt%flwds(n) = 0.0_rkx
     rt%fsds(n) = 0.0_rkx
     rt%fsnirt(n) = 0.0_rkx
     rt%fsnirtsq(n) = 0.0_rkx
     rt%fsnrtc(n) = 0.0_rkx
     rt%fsns(n) = 0.0_rkx
     rt%fsnsc(n) = 0.0_rkx
     rt%fsnt(n) = 0.0_rkx
     rt%fsntc(n) = 0.0_rkx
     rt%solin(n) = 0.0_rkx
     rt%solout(n) = 0.0_rkx
     rt%soll(n) = 0.0_rkx
     rt%solld(n) = 0.0_rkx
     rt%sols(n) = 0.0_rkx
     rt%solsd(n) = 0.0_rkx
    end do
    do concurrent ( k = 1:kz, n = 1:npr)
     rt%qrl(k,n) = 0.0_rkx
     rt%qrs(k,n) = 0.0_rkx
    end do
#endif
    !
    ! Dimensions
    !
    ni = 1+(ici2-ici1)
    nj = 1+(jci2-jci1)

    do concurrent ( j = jci1:jci2, i = ici1:ici2 )
      n = (j-jci1)+(i-ici1)*nj+1
    !
    ! NB: orography types are specified in the following
    !
      rt%ioro(n) = m2r%ldmsk(j,i)
    !
    ! Copy data from atmosphere, compute cloud particle size and
    ! fraction of ice and LW emissivity
    !
    ! Static informations
    !
      rt%dlat(n) = m2r%xlat(j,i)
      rt%xptrop(n) = m2r%ptrop(j,i)
    !
    ! Albedoes and surface emissivity
    !
      rt%emiss(n)  = m2r%emiss(j,i)
      rt%adirsw(n) = m2r%aldirs(j,i)
      rt%adifsw(n) = m2r%aldifs(j,i)
      rt%adirlw(n) = m2r%aldirl(j,i)
      rt%adiflw(n) = m2r%aldifl(j,i)
      rt%asw(n)    = m2r%albvs(j,i)
      rt%alw(n)    = m2r%albvl(j,i)
    !
    ! Sun elevation
    !
      rt%czen(n) = m2r%coszrs(j,i)
    end do
    do concurrent ( n = rt%n1:rt%n2 )
      if ( rt%czen(n) < 1.e-3_rkx ) then
        rt%czen(n) = 0.0_rkx
      end if
      if ( rt%czen(n) > 0.0_rkx ) then
        rt%czengt0(n) = .true.
      else
        rt%czengt0(n) = .false.
      end if
    end do
    !
    ! Gas concentrations
    !
    do concurrent ( k = 1:kz, m = 1:4, j = jci1:jci2, i = ici1:ici2 )
      n = (j-jci1)+(i-ici1)*nj+1
      rt%absgasnxt(k,m,n) = gasabsnxt(j,i,k,m)
    end do
    do concurrent ( k2 = 1:kzp1, k = 1:kzp1, j = jci1:jci2, i = ici1:ici2 )
      n = (j-jci1)+(i-ici1)*nj+1
      rt%absgastot(k2,k,n) = gasabstot(j,i,k2,k)
    end do
    do concurrent ( k = 1:kzp1, j = jci1:jci2, i = ici1:ici2 )
      n = (j-jci1)+(i-ici1)*nj+1
      rt%emsgastot(k,n) = gasemstot(j,i,k)
    end do
    !
    ! Surface pressure and scaled pressure, from which level are computed
    !
    do concurrent ( j = jci1:jci2, i = ici1:ici2 )
      n = (j-jci1)+(i-ici1)*nj+1
      rt%ps(n) = m2r%psatms(j,i)
      rt%ts(n) = m2r%tg(j,i)
    end do
    !
    ! rearrange input
    !
    do concurrent ( k = 1:kz, j = jci1:jci2, i = ici1:ici2 )
      n = (j-jci1)+(i-ici1)*nj+1
      rt%pmid(k,n) = m2r%phatms(j,i,k)
      rt%pmln(k,n) = log(m2r%phatms(j,i,k))
      rt%t(k,n)    = m2r%tatms(j,i,k)
      rt%rh(k,n)   = m2r%rhatms(j,i,k)
      rt%q(k,n)    = m2r%qxatms(j,i,k,iqv)
      rt%ql(k,n)   = m2r%qxatms(j,i,k,iqc)
      rt%dz(k,n)   = m2r%deltaz(j,i,k)
    end do
    do concurrent ( k = 1:kzp1, j = jci1:jci2, i = ici1:ici2 )
      n = (j-jci1)+(i-ici1)*nj+1
      rt%pint(k,n) = m2r%pfatms(j,i,k)
      rt%piln(k,n) = log(m2r%pfatms(j,i,k))
    end do
    !
    ! Air temperature and relative humidity
    !
    do concurrent ( k = 1:kz, n = rt%n1:rt%n2 )
      rt%rho(k,n) = (rt%pmid(k,n))/(rt%t(k,n)*rgas)
    end do
    if ( ipptls > 1 ) then
      do concurrent ( k = 1:kz, j = jci1:jci2, i = ici1:ici2 )
        n = (j-jci1)+(i-ici1)*nj+1
        rt%qi(k,n) = m2r%qxatms(j,i,k,iqi)
      end do
      do concurrent ( k = 1:kz, n = rt%n1:rt%n2 )
        if ( rt%qi(k,n) > 1.0e-11_rkx ) then
          if ( rt%ql(k,n) > 1.0e-11_rkx ) then
            rt%fice(k,n) = rt%qi(k,n) / (rt%ql(k,n)+rt%qi(k,n))
          else
            rt%fice(k,n) = 1.0_rkx
          end if
        else
          rt%fice(k,n) = 0.0_rkx
        end if
      end do
    else
      do concurrent ( k = 1:kz,  n = rt%n1:rt%n2 )
        ! Define fractional amount of cloud that is ice
        ! if warmer than -10 degrees C then water phase
        if ( rt%t(k,n) >= tzero ) then
          rt%fice(k,n) = 0.0_rkx
        else if ( rt%t(k,n) <= minus20 ) then
          !  if colder than -20 degrees C then ice phase
           rt%fice(k,n) = 1.0_rkx
        else
          ! if colder than 0 degrees C but warmer than -20 C mixed phase
          ! fice : fractional ice content within cloud
          rt%fice(k,n) = (tzero-rt%t(k,n))/20.0_rkx
        end if
        rt%qi(k,n) = rt%ql(k,n) * rt%fice(k,n)
        rt%ql(k,n) = rt%ql(k,n) - rt%qi(k,n)
      end do
    end if
    !
    ! Fractional cloud cover (dependent on relative humidity)
    ! Set cloud
    !   - NOT on the topmost two layers
    !   - Starting from ncld levels from the surface
    !
    do concurrent ( k = 1:kz, j = jci1:jci2, i = ici1:ici2 )
      n = (j-jci1)+(i-ici1)*nj+1
      ! Convert liquid water content into liquid water path
      rt%clwp(k,n) = m2r%cldlwc(j,i,k)*m2r%deltaz(j,i,k)
      ! O3 mass and volume mixing ratios
      rt%o3vmr(k,n) = 0.5_rkx*(o3prof(j,i,k+1)+o3prof(j,i,k))*(amd/amo3)
      ! Set a minimum clwp in the cloud of of 0.01 mm
      if ( rt%clwp(k,n) > 0.01_rkx ) then
        temp(k,n) = m2r%cldfrc(j,i,k)
      else
        rt%clwp(k,n) = 0.0_rkx
        temp(k,n) = 0.0_rkx
      end if
    end do
    do concurrent ( n = rt%n1:rt%n2 )
      rt%cld(kzp1,n) = temp(kz,n)
    end do
    do concurrent( k = 2:kz, n = rt%n1:rt%n2 )
      ! Use Maximum Random Overlap assumption
      rt%cld(k,n) = temp(k-1,n)+temp(k,n) - (temp(k-1,n)*temp(k,n))
    end do
    do concurrent ( n = rt%n1:rt%n2 )
      rt%cld(1,n) = 0.0_rkx
    end do
    if ( ncld > 0 ) then
      do concurrent ( k = kzp1-ncld:kzp1,  n = rt%n1:rt%n2 )
        rt%cld(k,n) = 0.0_rkx
      end do
    end if
    !
    ! only allow thin clouds (<0.25) above 400 mb (yhuang, 11/97)
    !
    !   do n = rt%n1, rt%n2
    !     do k = 1, kz
    !       if ( rt%pint(k+1,n) < 40000.0_rkx ) then
    !         rt%cld(k,n) = min(rt%cld(k,n),0.25_rkx)
    !       else
    !         rt%cld(k,n) = min(rt%cld(k,n),cftotmax)
    !       end if
    !     end do
    !   end do
    !
    !
    ! Tracers mixing ratios
    !
    if ( ichem == 1 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz, itr = 1:ntr )
        n = (j-jci1)+(i-ici1)*nj+1
        aermmr(n,k,itr) = max(m2r%chiatms(j,i,k,itr),0.0_rkx)
      end do
    end if
    !
    ! Compute cloud drop size and emissivity
    !
    do concurrent( k = 1:kz, n = rt%n1:rt%n2 )
      ! Define liquid drop size
      ! rel : liquid effective drop size (microns)
      !
      ! Global distribution of Cloud Droplet Effective Radius from
      ! POLDER polarization measurements
      ! On average, droplets are 2 to 3 micron smaller overland than over
      ! the ocean, but the smaller droplets are found over highly
      ! polluted regions and in areas affected by smoke from biomass
      ! burning activity. Largest droplets are found in remote tropical
      ! oceans, away from polarized reflectance of liquid clouds.
      ! Toward the poles, for colder temperature, large droplets freeze
      ! and only the smaller ones are kept liquid.
      !
      ! tpara = min(1.0_rkx,max(0.0_rkx,(minus10-rt%t(k,n))*0.05_rkx))
      ! if ( rt%ioro(n) == 1 ) then
      !   rt%rel(k,n) = 6.0_rkx + 3.0_rkx * tpara
      ! else
      !   rt%rel(k,n) = 11.0_rkx
      ! end if
      !
      ! if ( rt%ioro(n) == 1 ) then
      !   rt%rel(k,n) = 8.5_rkx
      ! else
      !   rt%rel(k,n) = 11.0_rkx
      ! end if
      !
      ! Long-wave optical properties of water clouds and rain
      ! Savijarvi,Raisanene (Tellus 1998)
      !
      ! g/m3, already account for cum and ls clouds
      if ( temp(k,n) > 0.0_rkx ) then
        lwc = (rt%ql(k,n) * rt%rho(k,n) * 1000.0_rkx)/temp(k,n)
        if ( rt%ioro(n) == 1 ) then
          ! Effective liquid radius over land
          rt%rel(k,n) = min(4.0_rkx + 7.0_rkx*lwc,15.0_rkx)
        else
          ! Effective liquid radius over ocean and sea ice
          rt%rel(k,n) = min(5.5_rkx + 9.5_rkx*lwc,18.0_rkx)
        end if
        !
        ! Stengel, Fokke Meirink, Eliasson (2023)
        ! On the Temperature Dependence of the Cloud Ice Particle Effective
        ! Radius - A Satellite Perspective
        !
        iwc = (rt%qi(k,n) * rt%rho(k,n) * 1000.0_rkx)/temp(k,n)
        tempc = rt%t(k,n) - 83.15_rkx
        tcels = tempc - tzero
        fsr = 1.2351_rkx + 0.0105_rkx * tcels
        aiwc = 45.8966_rkx * iwc**0.2214_rkx
        biwc = 0.7957_rkx * iwc**0.2535_rkx
        desr = fsr*(aiwc+biwc*tempc)
        desr = max(30.0_rkx,min(155.0_rkx,desr))
        ! rei : ice effective drop size (microns)
        rt%rei(k,n) = 0.64952_rkx*desr
      else
        ! filler for no clouds
        rt%rel(k,n) = 8.5_rkx
        rt%rei(k,n) = 20.0_rkx
      end if
    end do
#if 0
    !$acc kernels
    rt%totci(rt%n1:rt%n2) = 0.0_rkx
    rt%totcl(rt%n1:rt%n2) = 0.0_rkx
    rt%totwv(rt%n1:rt%n2) = 0.0_rkx
    !$acc end kernels
#else
    do concurrent ( n = rt%n1:rt%n2 )
      rt%totci(n) = 0.0_rkx
      rt%totcl(n) = 0.0_rkx
      rt%totwv(n) = 0.0_rkx
    end do
#endif
#ifdef STDPAR_FIXED
    !$acc parallel loop gang vector
    do concurrent ( n = rt%n1:rt%n2 )
#else
    !$acc parallel loop gang vector
    do n = rt%n1, rt%n2
#endif
      totci_l = 0.0_rkx
      totcl_l = 0.0_rkx
      totwv_l = 0.0_rkx
      !$acc loop seq
      do k = 1, kz
        ! Turn off ice radiative properties by setting fice = 0.0
        ! rt%fice(k,n) = 0.0_rkx
        totcl_l = totcl_l + &
             (rt%clwp(k,n)*rt%cld(k,n)*(1.0_rkx-rt%fice(k,n)))*0.001_rkx
        totci_l = totci_l + &
             (rt%clwp(k,n)*rt%cld(k,n)*rt%fice(k,n))*0.001_rkx
        totwv_l = totwv_l + rt%q(k,n)*rt%rho(k,n)*rt%dz(k,n)
      end do
      rt%totci(n)= totci_l
      rt%totcl(n)= totcl_l
      rt%totwv(n)= totwv_l
    end do

    !FAB : reintroduce simple sulfate indirect effect
    ! from Qian  1999
    ! rt%clwp is passed in g/m2
    if ( (ichem == 1 .and. iindirect == 1) .or. iclimaaer == 1 ) then
      do itr = 1, ntr
        if ( chtrname(itr) /= 'SO4' ) cycle
        do concurrent ( k = 1:kz, n = rt%n1:rt%n2 )
          aerc = rt%rho(k,n) * aermmr(k,n,itr) * 1.e9_rkx !microg/m3
          ! thershold of 0.1 microg/m3 for activation of indirect effect
          if ( aerc > 0.1_rkx ) then
            ! ccn number concentration in cm-3
            nc = 90.7_rkx * aerc**0.45_rkx + 23.0_rkx
            ! kg/m3, already account for cum and ls clouds
            lwc = (rt%clwp(k,n) / rt%dz(k,n)) * 0.001_rkx
            if ( lwc < 1.e-6_rkx ) cycle
            if ( rt%ioro(n) == 1 ) then
              ! Martin et al.(1994) parameter over land
              kparam = 0.67_rkx
            else
              ! Martin et al.(1994) parameter over ocean and sea ice
              kparam = 0.80_rkx
            end if
            !finally modify effective radius
            !(1.e6 to convert to rel to microm,
            ! 1.e6 to convert nc in m-3)
            rt%rel(k,n) = 1.e6_rkx * ( 3.0_rkx*lwc / &
              (4.0_rkx*mathpi*rhoh2o*kparam*nc*1.e6_rkx ) )**(1.0_rkx/3.0_rkx)
          end if
        end do
      end do
    end if
    do concurrent ( k = 1:kz, n = rt%n1:rt%n2 )
      ! ice absorption coefficient
      kabsi = 0.005_rkx + (1.0_rkx/rt%rei(k,n))
      ! liquid water absorption coefficient for broadband long wave
      ! (Hannu Savijärvi & Petri Räisänen)
      kabsl = 0.31_rkx * exp(-0.08_rkx*rt%rel(k,n))
      ! longwave absorption coeff (m**2/g)
      kabs = kabsl*(1.0_rkx-rt%fice(k,n)) + (kabsi*rt%fice(k,n))
      ! cloud emissivity (fraction)
      arg = min(kabs*rt%clwp(k,n),25.0_rkx)
      cldemis = max(1.0_rkx - exp(-arg),0.00_rkx)
      ! Effective cloud cover
      rt%effcld(k,n) = rt%cld(k,n)*cldemis
    end do
    rt%eccf = real(eccf,rkx)
    rt%labsem = labsem
    !
    ! Main radiation driving routine.
    ! NB: All fluxes returned from radctl() have already been converted to MKS.
    !
    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !############################
    call radctl(rt,iyear,imonth)
    !############################
    !@@@@@@@@@@@@@@@@@@@@@@@@@@@@
    !
    ! Save gas emission/absorbtion
    !
    if ( labsem ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz, m = 1:4 )
        n = (j-jci1)+(i-ici1)*nj+1
        gasabsnxt(j,i,k,m) = rt%absgasnxt(k,m,n)
      end do
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k2 = 1:kzp1, k = 1:kzp1 )
        n = (j-jci1)+(i-ici1)*nj+1
        gasabstot(j,i,k2,k) = rt%absgastot(k2,k,n)
      end do
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kzp1 )
        n = (j-jci1)+(i-ici1)*nj+1
        gasemstot(j,i,k) = rt%emsgastot(k,n)
      end do
    end if
    if ( ichem == 1 .or. iclimaaer == 1 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 0:kz, m = 1:nspi )
        n = (j-jci1)+(i-ici1)*nj+1
        taucldsp(j,i,k,m)  = rt%tauxcl(k,n,m) + rt%tauxci(k,n,m)
      end do
    end if
    !
    ! subroutine radout() copies back the data to RegCM for surface
    ! computations and output purposes.
    !
    call radout(lout,rt%solin,rt%solout,rt%fsns,rt%fsntc,rt%fsnsc,    &
                rt%qrs,rt%lwout,rt%flns,rt%flntc,rt%flnsc,rt%qrl,     &
                rt%flwds,rt%sols,rt%soll,rt%solsd,rt%solld,rt%totcf,  &
                rt%totwv,rt%totcl,rt%totci,rt%cld,rt%clwp,rt%abv,     &
                rt%sol,rt%aeradfo,rt%aeradfos,rt%aerlwfo,rt%aerlwfos, &
                tauxar3d,tauasc3d,gtota3d,rt%dz,rt%o3vmr,rt%outtaucl, &
                rt%outtauci,rt%asaeradfo,rt%asaeradfos,rt%asaerlwfo,  &
                rt%asaerlwfos,r2m,m2r)
#ifdef DEBUG
    call time_end(subroutine_name,indx)
#endif

  !@acc call nvtxEndRange
  end subroutine colmod3

end module mod_rad_colmod3

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

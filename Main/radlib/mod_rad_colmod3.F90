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
  use mod_intkinds , only : ik4
  use mod_realkinds , only : rkx
  use mod_constants , only : amd , amdk , amo3 , mathpi
  use mod_constants , only : rgasmol , rhoh2o , tzero
  use mod_dynparam , only : jci1 , jci2 , ici1 , ici2 , kz , kzp1
  use mod_dynparam , only : ntr , nspi
  use mod_memutil , only : getmem1d , getmem2d , getmem3d
  use mod_runparams , only : ipptls , eccf , ichem , idirect , iindirect
  use mod_runparams , only : iqv , iqc , iqi , ncld , chtrname
  use mod_runparams , only : cftotmax , iclimaaer
  use mod_rad_radiation , only : radini , radctl , radtype
  use mod_rad_common , only : gasabsnxt , gasabstot , gasemstot
  use mod_rad_common , only : o3prof , taucldsp , dosw , dolw , doabsems
  use mod_rad_outrad , only : radout
  use mod_rad_aerosol , only : aermmr , tauxar3d , tauasc3d , gtota3d
  use mod_regcm_types , only : mod_2_rad , rad_2_mod
#ifdef DEBUG
  use mod_service , only : time_begin , time_end , dbgslen
#endif

  implicit none

  private

  public :: allocate_mod_rad_colmod3 , colmod3

  type(radtype) :: rt

  contains

  subroutine allocate_mod_rad_colmod3
    implicit none
    integer(ik4) :: npr
    npr = (jci2-jci1+1)*(ici2-ici1+1)
    rt%n1 = 1
    rt%n2 = npr
    call getmem1d(rt%alb,1,npr,'colmod3:alb')
    call getmem1d(rt%albc,1,npr,'colmod3:albc')
    call getmem1d(rt%flns,1,npr,'colmod3:flns')
    call getmem1d(rt%flnsc,1,npr,'colmod3:flnsc')
    call getmem1d(rt%flnt,1,npr,'colmod3:flnt')
    call getmem1d(rt%lwout,1,npr,'colmod3:lwout')
    call getmem1d(rt%lwin,1,npr,'colmod3:lwin')
    call getmem1d(rt%flntc,1,npr,'colmod3:flntc')
    call getmem1d(rt%flwds,1,npr,'colmod3:flwds')
    call getmem1d(rt%fsds,1,npr,'colmod3:fsds')
    call getmem1d(rt%fsnirt,1,npr,'colmod3:fsnirt')
    call getmem1d(rt%fsnirtsq,1,npr,'colmod3:fsnirtsq')
    call getmem1d(rt%fsnrtc,1,npr,'colmod3:fsnrtc')
    call getmem1d(rt%fsns,1,npr,'colmod3:fsns')
    call getmem1d(rt%fsnsc,1,npr,'colmod3:fsnsc')
    call getmem1d(rt%fsnt,1,npr,'colmod3:fsnt')
    call getmem1d(rt%fsntc,1,npr,'colmod3:fsntc')
    call getmem1d(rt%solin,1,npr,'colmod3:solin')
    call getmem1d(rt%solout,1,npr,'colmod3:solout')
    call getmem1d(rt%soll,1,npr,'colmod3:soll')
    call getmem1d(rt%solld,1,npr,'colmod3:solld')
    call getmem1d(rt%sols,1,npr,'colmod3:sols')
    call getmem1d(rt%solsd,1,npr,'colmod3:solsd')
    call getmem1d(rt%totcf,1,npr,'colmod3:totcf')
    call getmem1d(rt%totwv,1,npr,'colmod3:totwv')
    call getmem1d(rt%totcl,1,npr,'colmod3:totcl')
    call getmem1d(rt%totci,1,npr,'colmod3:totci')
    call getmem1d(rt%ps,1,npr,'colmod3:ps')
    call getmem1d(rt%ts,1,npr,'colmod3:ts')
    call getmem1d(rt%emiss,1,npr,'colmod3:emiss')
    call getmem1d(rt%xptrop,1,npr,'rad:xptrop')
    call getmem1d(rt%dlat,1,npr,'rad:dlat')
    call getmem1d(rt%adirsw,1,npr,'rad:adirsw')
    call getmem1d(rt%adifsw,1,npr,'rad:adifsw')
    call getmem1d(rt%adirlw,1,npr,'rad:adirlw')
    call getmem1d(rt%adiflw,1,npr,'rad:adiflw')
    call getmem1d(rt%asw,1,npr,'rad:asw')
    call getmem1d(rt%alw,1,npr,'rad:alw')
    call getmem1d(rt%abv,1,npr,'rad:abv')
    call getmem1d(rt%sol,1,npr,'rad:sol')

    call getmem2d(rt%cld,1,npr,1,kzp1,'colmod3:cld')
    call getmem2d(rt%effcld,1,npr,1,kzp1,'colmod3:effcld')
    call getmem2d(rt%piln,1,npr,1,kzp1,'colmod3:piln')
    call getmem2d(rt%pint,1,npr,1,kzp1,'colmod3:pint')

    call getmem2d(rt%rh,1,npr,1,kz,'colmod3:rh')
    call getmem2d(rt%rho,1,npr,1,kz,'colmod3:rh')
    call getmem2d(rt%clwp,1,npr,1,kz,'colmod3:clwp')
    call getmem2d(rt%fice,1,npr,1,kz,'colmod3:fice')
    call getmem2d(rt%o3vmr,1,npr,1,kz,'colmod3:o3vmr')
    call getmem2d(rt%pmid,1,npr,1,kz,'colmod3:pmid')
    call getmem2d(rt%pmln,1,npr,1,kz,'colmod3:pmln')
    call getmem2d(rt%q,1,npr,1,kz,'colmod3:q')
    call getmem2d(rt%ql,1,npr,1,kz,'colmod3:ql')
    call getmem2d(rt%qi,1,npr,1,kz,'colmod3:qi')
    call getmem2d(rt%qrl,1,npr,1,kz,'colmod3:qrl')
    call getmem2d(rt%qrs,1,npr,1,kz,'colmod3:qrs')
    call getmem2d(rt%rei,1,npr,1,kz,'colmod3:rei')
    call getmem2d(rt%rel,1,npr,1,kz,'colmod3:rel')
    call getmem2d(rt%t,1,npr,1,kz,'colmod3:t')
    call getmem2d(rt%dz,1,npr,1,kz,'colmod3:dz')
    call getmem1d(rt%aeradfo,1,npr,'colmod3:aeradfo')
    call getmem1d(rt%aeradfos,1,npr,'colmod3:aeradfos')
    call getmem1d(rt%aerlwfo,1,npr,'colmod3:aerlwfo')
    call getmem1d(rt%aerlwfos,1,npr,'colmod3:aerlwfos')
    call getmem1d(rt%czen,1,npr,'colmod3:czen')
    call getmem1d(rt%czengt0,1,npr,'colmod3:czengt0')
    call getmem3d(rt%absgasnxt,1,npr,1,kz,1,4,'colmod3:absgasnxt')
    call getmem3d(rt%absgastot,1,npr,1,kzp1,1,kzp1,'colmod3:absgastot')
    call getmem2d(rt%emsgastot,1,npr,1,kzp1,'colmod3:emsgastot')
    call getmem3d(rt%tauxcl,1,npr,0,kz,1,nspi,'colmod3:tauxcl')
    call getmem3d(rt%tauxci,1,npr,0,kz,1,nspi,'colmod3:tauxci')
    call getmem3d(rt%outtaucl,1,npr,1,kzp1,1,4,'colmod3:outtaucl')
    call getmem3d(rt%outtauci,1,npr,1,kzp1,1,4,'colmod3:outtauci')

    call getmem1d(rt%ioro,1,npr,'colmod3:ioro')

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
    implicit none
    type(mod_2_rad) , intent(in) :: m2r
    type(rad_2_mod) , intent(inout) :: r2m
    integer(ik4) , intent(in) :: iyear , imonth
    logical , intent(in) :: lout , labsem

    integer(ik4) :: n , m , i , j , k , k2 , itr , kmincld , kmaxcld
    real(rkx) :: nc , aerc , lwc , kparam
    real(rkx) :: kabs , kabsi , kabsl , cldemis , arg
    real(rkx) :: iwc , tempc , tcels , fsr , aiwc , biwc , desr
    !real(rkx) :: tpara
    real(rkx) , parameter :: minus20 = 253.15_rkx
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
    !
    ! NB: orography types are specified in the following
    !
    n = 1
    do i = ici1 , ici2
      do j = jci1 , jci2
        rt%ioro(n) = m2r%ldmsk(j,i)
        n = n + 1
      end do
    end do
    !
    ! Copy data from atmosphere, compute cloud particle size and
    ! fraction of ice and LW emissivity
    !
    ! Static informations
    !
    n = 1
    do i = ici1 , ici2
      do j = jci1 , jci2
        rt%dlat(n) = m2r%xlat(j,i)
        rt%xptrop(n) = m2r%ptrop(j,i)
        n = n + 1
      end do
    end do
    !
    ! radini sets many radiation parameters
    !
    call radini(rt%n1,rt%n2,iyear,imonth,rt%dlat)
    !
    ! Albedoes and surface emissivity
    !
    n = 1
    do i = ici1 , ici2
      do j = jci1 , jci2
        rt%emiss(n)  = m2r%emiss(j,i)
        rt%adirsw(n) = m2r%aldirs(j,i)
        rt%adifsw(n) = m2r%aldifs(j,i)
        rt%adirlw(n) = m2r%aldirl(j,i)
        rt%adiflw(n) = m2r%aldifl(j,i)
        rt%asw(n)    = m2r%albvs(j,i)
        rt%alw(n)    = m2r%albvl(j,i)
        n = n + 1
      end do
    end do
    !
    ! Sun elevation
    !
    n = 1
    do i = ici1 , ici2
      do j = jci1 , jci2
        rt%czen(n) = m2r%coszrs(j,i)
        n = n + 1
      end do
    end do
    do n = rt%n1 , rt%n2
      if ( rt%czen(n) < 1.e-3_rkx ) then
        rt%czen(n) = 0.0_rkx
      end if
    end do
    do n = rt%n1 , rt%n2
      if ( rt%czen(n) > 0.0_rkx ) then
        rt%czengt0(n) = .true.
      else
        rt%czengt0(n) = .false.
      end if
    end do
    !
    ! Gas concentrations
    !
    do m = 1 , 4
      do k = 1 , kz
        n = 1
        do i = ici1 , ici2
          do j = jci1 , jci2
            rt%absgasnxt(n,k,m) = gasabsnxt(j,i,k,m)
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
            rt%absgastot(n,k2,k) = gasabstot(j,i,k2,k)
            n = n + 1
          end do
        end do
      end do
    end do
    do k = 1 , kzp1
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          rt%emsgastot(n,k) = gasemstot(j,i,k)
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
        rt%ps(n) = m2r%psatms(j,i)
        rt%ts(n) = m2r%tg(j,i)
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
          rt%pmid(n,k) = m2r%phatms(j,i,k)
          rt%pmln(n,k) = log(rt%pmid(n,k))
          n = n + 1
        end do
      end do
    end do
    do k = 1 , kzp1
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          rt%pint(n,k) = m2r%pfatms(j,i,k)
          rt%piln(n,k) = log(rt%pint(n,k))
          n = n + 1
        end do
      end do
    end do
    !
    ! Air temperature and relative humidity
    !
    do k = 1 , kz
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          rt%t(n,k) = m2r%tatms(j,i,k)
          rt%rh(n,k) = m2r%rhatms(j,i,k)
          rt%rho(n,k) = (rt%pmid(n,k)*amdk)/(rt%t(n,k)*rgasmol)
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
          ! Specific humidity (h2o mass mix ratio)
          rt%q(n,k) = m2r%qxatms(j,i,k,iqv)/(1.0_rkx+m2r%qxatms(j,i,k,iqv))
          rt%ql(n,k) = m2r%qxatms(j,i,k,iqc)
          n = n + 1
        end do
      end do
    end do
    if ( ipptls > 1 ) then
      do k = 1 , kz
        n = 1
        do i = ici1 , ici2
          do j = jci1 , jci2
            rt%qi(n,k) = m2r%qxatms(j,i,k,iqi)
            n = n + 1
          end do
        end do
      end do
      do k = 1 , kz
        do n = rt%n1 , rt%n2
          if ( rt%qi(n,k) > 1.0e-11_rkx ) then
            if ( rt%ql(n,k) > 1.0e-11_rkx ) then
              rt%fice(n,k) = rt%qi(n,k) / (rt%ql(n,k)+rt%qi(n,k))
            else
              rt%fice(n,k) = 1.0_rkx
            end if
          else
            rt%fice(n,k) = 0.0_rkx
          end if
        end do
      end do
    else
      do k = 1 , kz
        do n = rt%n1 , rt%n2
          ! Define fractional amount of cloud that is ice
          ! if warmer than -10 degrees C then water phase
          if ( rt%t(n,k) >= tzero ) then
            rt%fice(n,k) = 0.0_rkx
          else if ( rt%t(n,k) <= minus20 ) then
            !  if colder than -20 degrees C then ice phase
            rt%fice(n,k) = 1.0_rkx
          else
            ! if colder than 0 degrees C but warmer than -20 C mixed phase
            ! fice : fractional ice content within cloud
            rt%fice(n,k) = (tzero-rt%t(n,k))/20.0_rkx
          end if
        end do
      end do
      do k = 1 , kz
        n = 1
        do i = ici1 , ici2
          do j = jci1 , jci2
            rt%qi(n,k) = rt%ql(n,k) * rt%fice(n,k)
            rt%ql(n,k) = rt%ql(n,k) - rt%qi(n,k)
          end do
        end do
      end do
    end if
    !
    ! deltaz
    !
    do k = 1 , kz
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          rt%dz(n,k) = m2r%deltaz(j,i,k)
          n = n + 1
        end do
      end do
    end do
    !
    ! Fractional cloud cover (dependent on relative humidity)
    ! Set cloud
    !   - NOT on the topmost two layers
    !   - Starting from ncld levels from the surface
    !
    rt%cld  = 0.0_rkx
    rt%clwp = 0.0_rkx
    kmaxcld = 3
    kmincld = kz - ncld
    do k = kmaxcld , kmincld
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          ! Convert liquid water content into liquid water path
          rt%clwp(n,k) = m2r%cldlwc(j,i,k)*m2r%deltaz(j,i,k)
          n = n + 1
        end do
      end do
    end do
    do k = kmaxcld , kmincld
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( rt%clwp(n,k) > 0.0_rkx ) then
            ! Use Maximum Random Overlap assumption
            rt%cld(n,k) = m2r%cldfrc(j,i,k-1)+m2r%cldfrc(j,i,k) - &
                         (m2r%cldfrc(j,i,k-1)*m2r%cldfrc(j,i,k))
            rt%cld(n,k) = min(rt%cld(n,k),cftotmax)
          end if
          n = n + 1
        end do
      end do
    end do
    !
    ! only allow thin clouds (<0.25) above 400 mb (yhuang, 11/97)
    !
    !   do k = 1 , kz
    !     do n = rt%n1 , rt%n2
    !       if ( rt%pint(n,k+1) < 40000.0_rkx ) then
    !         rt%cld(n,k) = min(rt%cld(n,k),0.25_rkx)
    !       else
    !         rt%cld(n,k) = min(rt%cld(n,k),cftotmax)
    !       end if
    !     end do
    !   end do
    !
    ! O3 mass and volume mixing ratios
    !
    do k = 1 , kz
      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          rt%o3vmr(n,k) = 0.5_rkx*(o3prof(j,i,k+1)+o3prof(j,i,k))*(amd/amo3)
          n = n + 1
        end do
      end do
    end do
    !
    ! Tracers mixing ratios
    !
    if ( ichem == 1 ) then
      do itr = 1 , ntr
        do k = 1 , kz
          n = 1
          do i = ici1 , ici2
            do j = jci1 , jci2
              aermmr(n,k,itr) = max(m2r%chiatms(j,i,k,itr),0.0_rkx)
              n = n + 1
            end do
          end do
        end do
      end do
    end if
    !
    ! Compute cloud drop size and emissivity
    !
    rt%totci(rt%n1:rt%n2) = 0.0_rkx
    rt%totcl(rt%n1:rt%n2) = 0.0_rkx
    rt%totwv(rt%n1:rt%n2) = 0.0_rkx
    do k = 1 , kz
      do n = rt%n1 , rt%n2
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
        ! tpara = min(1.0_rkx,max(0.0_rkx,(minus10-rt%t(n,k))*0.05_rkx))
        ! if ( rt%ioro(n) == 1 ) then
        !   rt%rel(n,k) = 6.0_rkx + 3.0_rkx * tpara
        ! else
        !   rt%rel(n,k) = 11.0_rkx
        ! end if
        !
        ! if ( rt%ioro(n) == 1 ) then
        !   rt%rel(n,k) = 8.5_rkx
        ! else
        !   rt%rel(n,k) = 11.0_rkx
        ! end if
        !
        ! Long-wave optical properties of water clouds and rain
        ! Savijarvi,Raisanene (Tellus 1998)
        !
        ! g/m3, already account for cum and ls clouds
        lwc = (rt%clwp(n,k) / rt%dz(n,k))
        if ( rt%ioro(n) == 1 ) then
          ! Effective liquid radius over land
          rt%rel(n,k) = min(4.0_rkx + 7.0_rkx*lwc,15.0_rkx)
        else
          ! Effective liquid radius over ocean and sea ice
          rt%rel(n,k) = min(5.5_rkx + 9.5_rkx*lwc,18.0_rkx)
        end if
        !
        ! Stengel, Fokke Meirink, Eliasson (2023)
        ! On the Temperature Dependence of the Cloud Ice Particle Effective
        ! Radius - A Satellite Perspective
        !
        iwc = rt%qi(n,k) * rt%rho(n,k) * 1000.0_rkx
        tempc = rt%t(n,k) - 83.15_rkx
        tcels = tempc - tzero
        fsr = 1.2351_rkx + 0.0105_rkx * tcels
        aiwc = 45.8966_rkx * iwc**0.2214_rkx
        biwc = 0.7957_rkx * iwc**0.2535_rkx
        desr = fsr*(aiwc+biwc*tempc)
        desr = max(30.0_rkx,min(155.0_rkx,desr))
        ! rei : ice effective drop size (microns)
        rt%rei(n,k) = 0.64952_rkx*desr
      end do
    end do
    do k = 1 , kz
      do n = rt%n1 , rt%n2
        ! Turn off ice radiative properties by setting fice = 0.0
        ! rt%fice(n,k) = 0.0_rkx
        rt%totcl(n) = rt%totcl(n) + &
             (rt%clwp(n,k)*rt%cld(n,k)*(1.0_rkx-rt%fice(n,k)))*0.001_rkx
        rt%totci(n) = rt%totci(n) + &
             (rt%clwp(n,k)*rt%cld(n,k)*rt%fice(n,k))*0.001_rkx
        rt%totwv(n) = rt%totwv(n) + rt%q(n,k)*rt%rho(n,k)*rt%dz(n,k)
      end do
    end do

    !FAB : reintroduce simple sulfate indirect effect
    ! from Qian  1999
    ! rt%clwp is passed in g/m2
    if ( (ichem == 1 .and. iindirect == 1) .or. iclimaaer == 1 ) then
      do itr = 1 , ntr
        if ( chtrname(itr) /= 'SO4' ) cycle
        do k = 1 , kz
          do n = rt%n1 , rt%n2
            aerc = rt%rho(n,k) * aermmr(n,k,itr) * 1.e9_rkx !microg/m3
            ! thershold of 0.1 microg/m3 for activation of indirect effect
            if ( aerc > 0.1_rkx ) then
              ! ccn number concentration in cm-3
              nc = 90.7_rkx * aerc**0.45_rkx + 23.0_rkx
              ! kg/m3, already account for cum and ls clouds
              lwc = (rt%clwp(n,k) / rt%dz(n,k)) * 0.001_rkx
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
              rt%rel(n,k) = 1.e6_rkx * ( 3.0_rkx*lwc / &
                (4.0_rkx*mathpi*rhoh2o*kparam*nc*1.e6_rkx ) )**(1.0_rkx/3.0_rkx)
            end if
          end do
        end do
      end do
    end if
    do k = 1 , kz
      do n = rt%n1 , rt%n2
        ! ice absorption coefficient
        kabsi = 0.005_rkx + (1.0_rkx/rt%rei(n,k))
        ! liquid water absorption coefficient for broadband long wave
        ! (Hannu Savijärvi & Petri Räisänen)
        kabsl = 0.31_rkx * exp(-0.08_rkx*rt%rel(n,k))
        ! longwave absorption coeff (m**2/g)
        kabs = kabsl*(1.0_rkx-rt%fice(n,k)) + (kabsi*rt%fice(n,k))
        ! cloud emissivity (fraction)
        arg = min(kabs*rt%clwp(n,k),25.0_rkx)
        cldemis = max(1.0_rkx - exp(-arg),0.00_rkx)
        ! Effective cloud cover
        rt%effcld(n,k) = rt%cld(n,k)*cldemis
      end do
    end do
    !
    ! Main radiation driving routine.
    ! NB: All fluxes returned from radctl() have already been converted to MKS.
    !
    rt%eccf = real(eccf,rkx)
    rt%labsem = labsem
    call radctl(rt)
    !
    ! Save gas emission/absorbtion
    !
    if ( labsem ) then
      do m = 1 , 4
        do k = 1 , kz
          n = 1
          do i = ici1 , ici2
            do j = jci1 , jci2
              gasabsnxt(j,i,k,m) = rt%absgasnxt(n,k,m)
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
              gasabstot(j,i,k2,k) = rt%absgastot(n,k2,k)
              n = n + 1
            end do
          end do
        end do
      end do
      do k = 1 , kzp1
        n = 1
        do i = ici1 , ici2
          do j = jci1 , jci2
            gasemstot(j,i,k) = rt%emsgastot(n,k)
            n = n + 1
          end do
        end do
      end do
    end if
    if ( ichem == 1 .or. iclimaaer == 1 ) then
      do m = 1 , nspi
        do k = 0 , kz
          n = 1
          do i = ici1 , ici2
            do j = jci1 , jci2
              taucldsp(j,i,k,m)  = rt%tauxcl(n,k,m) + rt%tauxci(n,k,m)
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

  end subroutine colmod3

end module mod_rad_colmod3

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

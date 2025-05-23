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

module mod_bats_common
!
! Storage for Surface (BATS and shared by CLM) variables
!
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_runparams, only : ichem, iemiss, syncro_srf, replacemoist
  use mod_runparams, only : rhmax, rhmin
  use mod_runparams, only : rcmtimer
  use mod_mppparam
  use mod_mpmessage
  use mod_constants
  use mod_service
  use mod_stdio
  use mod_bats_param
  use mod_bats_internal
  use mod_bats_bndry
  use mod_bats_leaftemp
  use mod_bats_drag
  use mod_bats_albedo
  use mod_regcm_types

  implicit none

  private

  public :: dtbat
  public :: allocate_mod_bats_internal
  public :: interf, initbats, vecbats, albedobats

  real(rkx), public :: rdnnsg

  real(rkx), public, pointer, contiguous, dimension(:,:,:) :: xqs

  contains

  !
  ! provides initial fields to boundary subroutine
  ! units are si
  !
  subroutine initbats(lm,lms)
    implicit none
    type(lm_exchange), intent(inout) :: lm
    type(lm_state), intent(inout) :: lms
    integer(ik4) :: i, j, n
    logical, parameter :: snowhack = .false.
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'initbats'
    integer(ik4), save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    if ( nsg > 1 ) then
      call getmem3d(xqs,1,nnsg,jci1,jci2,ici1,ici2,'bats:xqs')
    end if
    if ( rcmtimer%start( ) ) then
      call c2l_gs(lndcomm,lm%ht,ht)
      call c2l_ss(lndcomm,lm%ht1,hts)
      call c2l_ss(lndcomm,lm%iveg1,lveg)
      call c2l_ss(lndcomm,lm%itex1,ltex)
      call c2l_ss(lndcomm,lm%xlat1,lat)
      call c2l_gs(lndcomm,lm%tg,tgrd)
      call c2l_gs(lndcomm,lm%snowam,sncv)
      call c2l_gs(lndcomm,lm%zencos,czenith)
      ! Reset if invalid
      where ( sncv > 1.0e+10_rkx ) sncv = d_zero
      ! Remove water -> no-data
      where( ltex == 14 )
        ltex = 17
      end where
      tgbrd = tgrd
      tlef  = tgrd
      taf   = tgrd
      if ( snowhack ) then
        where ( tgrd < 263.0_rkx .and. hts > 500.0_rkx * egrav )
          sncv = hts*d_r100*regrav*(1.0_rkx - &
                 (min(1.0_rkx,max((tgrd-253.0_rkx)/20.0_rkx,0.0_rkx))))
        end where
      end if
      if ( lsmoist ) then
        call c2l_gs(lndcomm,lm%smoist,ssw)
        do i = ilndbeg, ilndend
          ! Initialize surface soil moisture
          ssw(i) = depuv(lveg(i))*xmopor(ltex(i))*ssw(i)
          rsw(i) = deprv(lveg(i))*xmopor(ltex(i))*slmo(lveg(i))
          tsw(i) = deptv(lveg(i))*xmopor(ltex(i))*slmo(lveg(i))
          gwet(i) = d_half
        end do
      else
        do i = ilndbeg, ilndend
          ! Initialize soil moisture in the 3 layers
          ssw(i) = depuv(lveg(i))*xmopor(ltex(i))*slmo(lveg(i))
          rsw(i) = deprv(lveg(i))*xmopor(ltex(i))*slmo(lveg(i))
          tsw(i) = deptv(lveg(i))*xmopor(ltex(i))*slmo(lveg(i))
          gwet(i) = d_half
        end do
      end if
      if ( replacemoist ) then
        if ( myid == italk ) then
          write(stdout,*) 'Initializing moisture from DOMAIN file'
        end if
        ! Replace surface and root soil moisture
        call l2c_ss(lndcomm,ssw,lms%ssw)
        do i = ici1, ici2
          do j = jci1, jci2
            if (lm%rmoist(j,i,1) < 1.0e+10_rkx ) then
              do n = 1, nsg
                lms%ssw(n,:,:) = lm%rmoist(:,:,1)
              end do
            end if
          end do
        end do
        call c2l_ss(lndcomm,lms%ssw,ssw)
        call l2c_ss(lndcomm,rsw,lms%rsw)
        do i = ici1, ici2
          do j = jci1, jci2
            if (lm%rmoist(j,i,2) < 1.0e+10_rkx ) then
              do n = 1, nsg
                lms%rsw(n,:,:) = lm%rmoist(:,:,2)
              end do
            end if
          end do
        end do
        call c2l_ss(lndcomm,lms%rsw,rsw)
        call l2c_ss(lndcomm,tsw,lms%tsw)
        do i = ici1, ici2
          do j = jci1, jci2
            if (lm%rmoist(j,i,3) < 1.0e+10_rkx ) then
              do n = 1, nsg
                lms%tsw(n,:,:) = lm%rmoist(:,:,3)
              end do
            end if
          end do
        end do
        call c2l_ss(lndcomm,lms%tsw,tsw)
      end if
      !
      ! Calculate emission coefficients
      !
      call fseas(tgbrd,aseas)
      if ( iemiss == 1 ) then
        do i = ilndbeg, ilndend
          emiss(i) = lndemiss(lveg(i)) - seasemi(lveg(i)) * aseas(i)
        end do
      else
        do i = ilndbeg, ilndend
          emiss(i) = lnd_sfcemiss
        end do
      end if
    else
      call c2l_gs(lndcomm,lm%ht,ht)
      call c2l_ss(lndcomm,lm%ht1,hts)
      call c2l_ss(lndcomm,lm%iveg1,lveg)
      call c2l_ss(lndcomm,lm%itex1,ltex)
      call c2l_ss(lndcomm,lm%xlat1,lat)
      call c2l_ss(lndcomm,lms%tgrd,tgrd)
      call c2l_ss(lndcomm,lms%tgbrd,tgbrd)
      call c2l_ss(lndcomm,lms%gwet,gwet)
      call c2l_ss(lndcomm,lms%sncv,sncv)
      call c2l_ss(lndcomm,lms%snag,snag)
      call c2l_ss(lndcomm,lms%ldew,ldew)
      call c2l_ss(lndcomm,lms%taf,taf)
      call c2l_ss(lndcomm,lms%emisv,emiss)
      call c2l_ss(lndcomm,lms%tlef,tlef)
      call c2l_ss(lndcomm,lms%ssw,ssw)
      call c2l_ss(lndcomm,lms%rsw,rsw)
      call c2l_ss(lndcomm,lms%tsw,tsw)
      call c2l_gs(lndcomm,lm%zencos,czenith)
      ! Remove water -> no-data
      where( ltex == 14 )
        ltex = 17
      end where
      call fseas(tgbrd,aseas)
    end if
    ! Prepare for albedo computation
    do i = ilndbeg, ilndend
      lncl(i) = mfcv(lveg(i)) - seasf(lveg(i))*aseas(i)
    end do
    call depth
    dzh = hts - ht
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine initbats
  !
  !=======================================================================
  !  based on: bats version 1e          copyright 18 august 1989
  !=======================================================================
  !      written by robert e. dickinson and patrick j. kennedy with
  ! contributions from klaus bluemel, filippo giorgi, and ann henderson-
  ! sellers.  version 1b (1 dec 1988) greatly modifies versions 1 and 1a,
  ! for example:  general cleanup and stabilization of iterations;
  ! three-layer soil moisture; improvement and bug fixes by bluemel
  ! including an improved canopy model; improvement in soil temperature
  ! force-restore method for overlying snow; stomatal resistance
  ! dependence on solar radiation made dependent on lai; inclusion of a
  ! non-zero displacement level; elimination of soil water movement
  ! in frozen ground; fixing of the carbon routine;
  ! and derivation of leaf iteration scheme.
  !
  !      version 1c (1 march 1989) modifies version 1b by:
  ! adding a dependence of stomatal resistance on vapor pressure deficit,
  ! allowing for solar zenith dependence of incident light in
  ! calculating dependence of stomatal resistance on light;
  ! further streamlining the iterative calculation of leaf temperature;
  ! and fixing a bug in the soil moisture calculation of bats1b.
  !
  !      version 1d (1 june 1989) fixes two bugs introduced into 1c,
  ! in subroutines lftemp and tgrund, and includes minor changes
  ! suggested by k. bluemel.  it redefines output fluxes of sensible and
  ! latent heat to agree exactly with those used for soil energy balance
  ! to ensure conservation of energy.  some other cleanup is done also.
  !
  !      verson 1e (18 august 1989) makes changes consistent with 1d to
  ! insure water conservation.  latent heat of sublimation now applies
  ! to all soil terms where snow-covered or below freezing.  two minor
  ! bugs from misplaced lines of code are fixed.  leaf temperature
  ! calculation reformulated to use average longwave rather than over
  ! bare soil as done previously, to better match ccm1 tie-in.
  !
  !      vector bats (1989) was developed from bats1e by gary bates.
  ! it is designed such that one call to bats perfoms calculations at
  ! a specified arbitrary number of gridpoints instead of just one point
  ! as in the standard version of bats1e.  this allows a more efficient
  ! vectorization of the code.  vector bats was coupled to the mm4
  ! in februrary 1992 by gary bates.
  !
  !                  main drive program
  !
  ! the core model begins with subroutine bndry, which is driven by
  ! this routine.
  !
  ! *********************************************************************
  !
  !
  !                      * ts (lowest model layer "midpt" at 75m)
  !
  !
  !                         i-----------i
  !                         i         ======= taf (temp air in leaves)
  !                         i     tlef  i
  !                         i-----------i
  !                              i  i
  !          ta (anemom)         i  i
  !          is actually a ss    i  i
  !             i   i temp (4')  i  i
  !  -----------i---i------------i--i------------------------------------
  !                                                  tg
  ! -------------------------------------------------------------------
  !                                                    tgb
  ! ********************************************************************
  !  **  type1-crop
  !  **  type2-short grass                  ****************************
  !  **  type3-evergreen needle leaf tree   * soils parameters are a fn of
  !  **  type4-deciduous ""  """ "" "" "    * soil colour -- from 1(light)
  !  **  type5- "" """"  broadleaf tree     *                  to 8(dark)
  !  **  type6- evergreen   """  """"       * soil texture -- from 1(sand)
  !  **  type7-tall grass                   *     thru 6(loam) to 12(clay)
  !  **  type8-desert                       *
  !  **  type9-tundra                       * soil drainage - 7(free),
  !  **  type10-irrig crop                  *      8(poor), 9(impeded)
  !  **  type11-semi-desert                 *
  !  **  type12-ice                         *  (drainage to be used yet)
  !  **  type13-bog or marsh                *
  !  **  type14-(inland water)              *
  !  **  type15-(sea)                       *****************************
  !  **  type16-evgr shrub
  !  **  type17-decid shrub
  !  **  type18-mixed tree
  !
  !
  !  flow from this driver and order of subroutines is:
  !
  !   initbats
  !   vecbats ==> soilbc
  !                bndry   ==>   drag  ==> dragdn  ==> depth
  !                             satur
  !                            vcover
  !                              drip
  !                            lftemp   =====================> stomat
  !                              drip                          frawat
  !                               co2 ===> carbon             frawat
  !                            tseice                           root
  !                            tgrund                          satur
  !                              snow                         lfdrag
  !                             water                         condch
  !                                                           condcq
  !                                                            deriv
  !
  !=======================================================================
  !  si version  - water fluxes are generally calculated in kg/m**2/s.
  !  1000 kg/m**2/s = 1 m/s  - converted to energy units for display:
  !                            1 kg/m**2/s = 2.5 x 10.e6  watts/m**2.
  !  note also  1 kg/m**2/s = 1 mm/m**2/s so fluxes can be thought of
  !                            in mm/s.
  !=======================================================================
  !
  subroutine vecbats(lm,lms)
    implicit none
    type(lm_exchange), intent(inout) :: lm
    type(lm_state), intent(inout) :: lms
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'vecbats'
    integer(ik4), save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !----------------------------
    ! Excange from model to BATS
    !----------------------------
    call interf(lm,lms,1)
    ! Calculate surface fluxes and hydrology budgets
    call soilbc
    call bndry
    !----------------------------
    ! Excange from BATS to model
    !----------------------------
    call interf(lm,lms,2)
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine vecbats
  !
  !  this subroutine interfaces regcm and bats variables
  !
  !  ivers = 1,   regcm --> bats
  !  ivers = 2,   bats --> regcm
  !
  subroutine interf(lm,lms,ivers)
    implicit none
    type(lm_exchange), intent(inout) :: lm
    type(lm_state), intent(inout) :: lms
    integer(ik4), intent (in) :: ivers
    real(rkx) :: facb, facs, fact, factuv, facv, fracb,  &
                 fracs, fracv, rh0, solvt, xqs0, xqsdif, wspd
    integer(ik4) :: i, j, n
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'interf'
    integer(ik4), save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    if ( ivers == 1 ) then ! regcm2d --> bats

      call c2l_gs(lndcomm,lm%hgt,zh)
      call c2l_gs(lndcomm,lm%sfps,p0)
      call c2l_gs(lndcomm,lm%sfta,ts0)
      call c2l_gs(lndcomm,lm%qvatm,qs0)
      call c2l_gs(lndcomm,lm%swdir,swd0)
      call c2l_gs(lndcomm,lm%swdif,swf0)
      call c2l_gs(lndcomm,lm%zencos,czenith)
      call c2l_gs(lndcomm,lm%ncprate,ncpr0)
      call c2l_gs(lndcomm,lm%cprate,cpr0)
      call c2l_gs(lndcomm,lm%hfx,sent)
      call c2l_gs(lndcomm,lm%qfx,evpr)
      call c2l_gs(lndcomm,lm%uatm,usw)
      call c2l_gs(lndcomm,lm%vatm,vsw)
      call c2l_gs(lndcomm,lm%solar,swsi)
      call c2l_gs(lndcomm,lm%rswf,swflx)
      call c2l_gs(lndcomm,lm%rlwf,lwflx)
      call c2l_gs(lndcomm,lm%vegswab,abswveg)

      z1log  = log(zh)
      z2fra  = log(zh*d_half)
      z10fra = log(zh*d_r10)
      zlglnd = log(zh/zlnd)
      zlgsno = log(zh/zsno)

      call fseas(tgbrd,aseas)
      do i = ilndbeg, ilndend
        lncl(i) = mfcv(lveg(i)) - seasf(lveg(i))*aseas(i)
      end do
      call depth

      do i = ilndbeg, ilndend
        xqs0 = pfwsat(ts0(i),p0(i))
        rh0 = min(max(qs0(i)/xqs0,rhmin),rhmax)
        solvt = swd0(i) + swf0(i)
        sts(i) = ts0(i)-lrate*regrav*dzh(i)
        sfcp(i) = p0(i)*(sts(i)/ts0(i))
        xqs0 = pfwsat(sts(i),sfcp(i))
        qs(i) = max(rh0*xqs0,d_zero)
        ! Move to specific humidity
        qs(i) = qs(i)/(d_one+qs(i))
        rhs(i) = sfcp(i)/(rgas*sts(i))
        ! Average over the period
        prcp(i) = (ncpr0(i) + cpr0(i))*syncro_srf%rw
        !
        ! quantities stored on 2d surface array for bats use only
        !
        zlgveg(i) = log(zh(i)/rough(lveg(i)))
        zlgdis(i) = log((zh(i)-displa(lveg(i)))/rough(lveg(i)))
        if ( solvt > d_zero ) then
          fracd(i) = swf0(i)/solvt
        else
          fracd(i) = 0.2_rkx
        end if
      end do

      ! Smooth out big qs differences
      if ( nsg > 1 ) then
        call l2c_ss(lndcomm,qs,xqs)
        do i = ici1, ici2
          do j = jci1, jci2
            xqsdif = d_zero
            ! Coarse scale specific humidity
            xqs0 = lm%qvatm(j,i)/(d_one+lm%qvatm(j,i))
            do n = 1, nnsg
              if ( lm%ldmsk1(n,j,i) == 1 ) then
                xqsdif = xqsdif + (xqs(n,j,i)-xqs0)
              end if
            end do
            ! Mean difference from subgrid values and coarse grid
            xqsdif = xqsdif*rdnnsg
            do n = 1, nnsg
              xqs(n,j,i) = max(xqs(n,j,i)-xqsdif,d_zero)
            end do
          end do
        end do
        call c2l_ss(lndcomm,xqs,qs)
      end if

    else if ( ivers == 2 ) then ! bats --> regcm2d

      call fseas(tgbrd,aseas)
      do i = ilndbeg, ilndend
        lncl(i) = mfcv(lveg(i)) - seasf(lveg(i))*aseas(i)
      end do
      call depth

      if ( iemiss == 1 ) then
        do i = ilndbeg, ilndend
          fracs = lncl(i)*wt(i) + (d_one-lncl(i))*scvk(i)
          emiss(i) = (lndemiss(lveg(i))-seasemi(lveg(i))*aseas(i)) * &
                  (d_one-fracs) + 0.992_rkx*fracs
        end do
      end if
      call l2c_ss(lndcomm,tgbrd,lms%tgbrd)
      call l2c_ss(lndcomm,tlef,lms%tlef)
      call l2c_ss(lndcomm,tgrd,lms%tgrd)
      call l2c_ss(lndcomm,gwet,lms%gwet)
      call l2c_ss(lndcomm,snag,lms%snag)
      call l2c_ss(lndcomm,ldew,lms%ldew)
      call l2c_ss(lndcomm,sncv,lms%sncv)
      call l2c_ss(lndcomm,ssw,lms%ssw)
      call l2c_ss(lndcomm,rsw,lms%rsw)
      call l2c_ss(lndcomm,tsw,lms%tsw)
      call l2c_ss(lndcomm,emiss,lms%emisv)
      call l2c_ss(lndcomm,taf,lms%taf)
      call l2c_ss(lndcomm,sfcp,lms%sfcp)
      call l2c_ss(lndcomm,sent,lms%sent)
      call l2c_ss(lndcomm,evpr,lms%evpr)
      call l2c_ss(lndcomm,prcp,lms%prcp)
      call l2c_ss(lndcomm,trnof,lms%trnof)
      call l2c_ss(lndcomm,srnof,lms%srnof)
      call l2c_ss(lndcomm,drag,lms%drag)
      call l2c_ss(lndcomm,sm,lms%snwm)
      call l2c_ss(lndcomm,rib,lms%br)
      call l2c_ss(lndcomm,cgrnds,lms%rah1)
      call l2c_ss(lndcomm,cdrx,lms%ram1)

      do i = ici1, ici2
        do j = jci1, jci2
          do n = 1, nnsg
            if ( lm%ldmsk1(n,j,i) == 1 ) then
              wspd = max(sqrt(lm%uatm(j,i)*lm%uatm(j,i) + &
                              lm%vatm(j,i)*lm%vatm(j,i)),0.1_rkx)
              lms%rah1(n,j,i) = d_one/lms%rah1(n,j,i)
              lms%ram1(n,j,i) = wspd/lms%ram1(n,j,i)
            end if
          end do
        end do
      end do

      do i = ilndbeg, ilndend
        fracv = sigf(i)
        fracb = (d_one-lncl(i))*(d_one-scvk(i))
        fracs = lncl(i)*wt(i) + (d_one-lncl(i))*scvk(i)
        zo(i) = fracv * rough(lveg(i)) + fracb * zlnd + fracs * zsno
        facv = z2fra(i)/zlgveg(i)
        facb = z2fra(i)/zlglnd(i)
        facs = z2fra(i)/zlgsno(i)
        fact = fracv*facv + fracb*facb + fracs*facs
        qs(i) = qs(i) - delq(i)*fact
        sts(i) = sts(i) - delt(i)*fact
        facv = z10fra(i)/zlgveg(i)
        facb = z10fra(i)/zlglnd(i)
        facs = z10fra(i)/zlgsno(i)
        factuv = fracv*facv + fracb*facb + fracs*facs
        usw(i) = usw(i)*(d_one-factuv)
        vsw(i) = vsw(i)*(d_one-factuv)
        tgbb(i) = ((d_one-lncl(i))*tgrd(i)**4+(lncl(i))*tlef(i)**4)**d_rfour
      end do

      call l2c_ss(lndcomm,usw,lms%u10m)
      call l2c_ss(lndcomm,vsw,lms%v10m)
      call l2c_ss(lndcomm,qs,lms%q2m)
      call l2c_ss(lndcomm,sts,lms%t2m)
      call l2c_ss(lndcomm,tgbb,lms%tgbb)
      call l2c_ss(lndcomm,zo,lms%zo)

      if ( ichem == 1 ) then
        call l2c_ss(lndcomm,delt,lms%deltat)
        call l2c_ss(lndcomm,delq,lms%deltaq)
        call l2c_ss(lndcomm,xlai,lms%xlai)
        call l2c_ss(lndcomm,wt,lms%wt)
        call l2c_ss(lndcomm,scvk,lms%scvk)
        call l2c_ss(lndcomm,sigf,lms%sigf)
        call l2c_ss(lndcomm,lncl,lms%lncl)
      end if
    end if ! Versus of the interface (1,2)

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif

    contains

#include <pfesat.inc>
#include <pfwsat.inc>

  end subroutine interf

  subroutine albedobats(lm,lms)
    implicit none
    type(lm_exchange), intent(inout) :: lm
    type(lm_state), intent(inout) :: lms
    call albedo
    call l2c_ss(lndcomm,swal,lms%swalb)
    call l2c_ss(lndcomm,lwal,lms%lwalb)
    call l2c_ss(lndcomm,swdiral,lms%swdiralb)
    call l2c_ss(lndcomm,lwdiral,lms%lwdiralb)
    call l2c_ss(lndcomm,swdifal,lms%swdifalb)
    call l2c_ss(lndcomm,lwdifal,lms%lwdifalb)
  end subroutine albedobats

end module mod_bats_common
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

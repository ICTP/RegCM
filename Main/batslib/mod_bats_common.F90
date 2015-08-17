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

module mod_bats_common
!
! Storage for Surface (BATS and shared by CLM) variables
!
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_runparams , only : ichem , iemiss , rtsrf , ktau , replacemoist
  use mod_mppparam
  use mod_mpmessage
  use mod_constants
  use mod_service
  use mod_bats_param
  use mod_bats_internal
  use mod_bats_bndry
  use mod_bats_leaftemp
  use mod_bats_albedo
  use mod_regcm_types

  implicit none

  private

  public :: dtbat
  public :: ldesseas
  public :: allocate_mod_bats_internal
  public :: interf , initbats , vecbats , albedobats

  real(rk8) , public :: rdnnsg

  real(rk8) , public , pointer , dimension(:,:,:) :: xqs

  contains

  !
  ! provides initial fields to boundary subroutine
  ! units are si
  !
  subroutine initbats(lm,lms)
    implicit none
    type(lm_exchange) , intent(inout) :: lm
    type(lm_state) , intent(inout) :: lms
    integer(ik4) :: i , itex
    logical , parameter :: snowhack = .false.
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'initbats'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    if ( nsg > 1 ) then
      call getmem3d(xqs,1,nnsg,jci1,jci2,ici1,ici2,'bats:xqs')
    end if
    if ( ktau == 0 ) then
      call c2l_gs(lndcomm,lm%ht,ht)
      call c2l_ss(lndcomm,lm%ht1,hts)
      call c2l_ss(lndcomm,lm%iveg1,lveg)
      call c2l_ss(lndcomm,lm%xlat1,lat)
      call c2l_gs(lndcomm,lm%tground2,tgrd)
      call c2l_gs(lndcomm,lm%snowam,sncv)

      tgbrd = tgrd
      tlef  = tgrd
      taf   = tgrd
      if ( snowhack ) then
        where ( tgrd < 263.0D0 )
          sncv = hts/10.0D0
        end where
      end if

      if ( lsmoist ) then
        call c2l_gs(lndcomm,lm%smoist,gwet)
        do i = ilndbeg , ilndend
          itex  = iexsol(lveg(i))
          ! Initialize soil moisture in the 3 layers
          tsw(i) = gwet(i)*deptv(lveg(i))
          rsw(i) = gwet(i)*deprv(lveg(i))
          ssw(i) = gwet(i)*depuv(lveg(i))
          gwet(i) = d_half
        end do
      else
        do i = ilndbeg , ilndend
          itex  = iexsol(lveg(i))
          ! Initialize soil moisture in the 3 layers
          tsw(i) = deptv(lveg(i))*xmopor(itex)*slmo(lveg(i))
          rsw(i) = deprv(lveg(i))*xmopor(itex)*slmo(lveg(i))
          ssw(i) = depuv(lveg(i))*xmopor(itex)*slmo(lveg(i))
          gwet(i) = d_half
        end do
      end if
      if ( replacemoist ) then
        if ( myid == italk ) then
          write(stdout,*) 'Initializing moisture from DOMAIN file'
        end if
        ! Replace soil moisture
        lm%smoist = lm%rmoist(:,:,2)
        call c2l_gs(lndcomm,lm%smoist,rsw)
        lm%smoist = lm%rmoist(:,:,1)
        call c2l_gs(lndcomm,lm%smoist,ssw)
        tsw = rsw + ssw
      end if
      !
      ! Calculate emission coefficients
      !
      call fseas(tgbrd,aseas)
      if ( iemiss == 1 ) then
        do i = ilndbeg , ilndend
          emiss(i) = lndemiss(lveg(i)) - seasemi(lveg(i)) * aseas(i)
        end do
      else
        do i = ilndbeg , ilndend
          emiss(i) = lnd_sfcemiss
        end do
      end if
    else
      call c2l_gs(lndcomm,lm%ht,ht)
      call c2l_ss(lndcomm,lm%ht1,hts)
      call c2l_ss(lndcomm,lm%iveg1,lveg)
      call c2l_ss(lndcomm,lm%xlat1,lat)
      call c2l_ss(lndcomm,lms%tgrd,tgrd)
      call c2l_ss(lndcomm,lms%tgbrd,tgbrd)
      call c2l_ss(lndcomm,lms%sncv,sncv)
      call c2l_ss(lndcomm,lms%scvk,scvk)
      call c2l_ss(lndcomm,lms%gwet,gwet)
      call c2l_ss(lndcomm,lms%snag,snag)
      call c2l_ss(lndcomm,lms%scvk,scvk)
      call c2l_ss(lndcomm,lms%ldew,ldew)
      call c2l_ss(lndcomm,lms%taf,taf)
      call c2l_ss(lndcomm,lms%emisv,emiss)
      call c2l_ss(lndcomm,lms%tlef,tlef)
      call c2l_ss(lndcomm,lms%ssw,ssw)
      call c2l_ss(lndcomm,lms%rsw,rsw)
      call c2l_ss(lndcomm,lms%tsw,tsw)
      call c2l_gs(lndcomm,lm%zencos,czenith)
      call fseas(tgbrd,aseas)
    end if
    do i = ilndbeg , ilndend
      lncl(i) = mfcv(lveg(i)) - seasf(lveg(i))*aseas(i)
    end do
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
    type(lm_exchange) , intent(inout) :: lm
    type(lm_state) , intent(inout) :: lms
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'vecbats'
    integer(ik4) , save :: idindx = 0
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
  !  ivers = 1 ,   regcm --> bats
  !  ivers = 2 ,   bats --> regcm
  !
  subroutine interf(lm,lms,ivers)
    implicit none
    type(lm_exchange) , intent(inout) :: lm
    type(lm_state) , intent(inout) :: lms
    integer(ik4) , intent (in) :: ivers
    real(rk8) :: facb , facs , fact , factuv , facv , fracb ,  &
                 fracs , fracv , rh0 , solvt , xqs0 , xqsdif
    integer(ik4) :: i , j , n
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'interf'
    integer(ik4) , save :: idindx = 0
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

      z1log  = dlog(zh)
      z2fra  = dlog(zh*d_half)
      z10fra = dlog(zh*d_r10)
      zlglnd = dlog(zh/zlnd)
      zlgsno = dlog(zh/zsno)

      call fseas(tgbrd,aseas)

      do i = ilndbeg , ilndend
        xqs0 = pfqsat(ts0(i),p0(i))
        rh0 = max(qs0(i)/xqs0,d_zero)
        solvt = swd0(i) + swf0(i)
        sts(i) = ts0(i)-lrate*regrav*dzh(i)
        sfcp(i) = p0(i)*(sts(i)/ts0(i))
        xqs0 = pfqsat(sts(i),sfcp(i))
        qs(i) = max(rh0*xqs0,d_zero)
        ! Move to specific humidity
        qs(i) = qs(i)/(d_one+qs(i))
        rhs(i) = sfcp(i)/(rgas*sts(i))
        ! Average over the period
        prcp(i) = (ncpr0(i) + cpr0(i))*rtsrf
        !
        ! quantities stored on 2d surface array for bats use only
        !
        lncl(i) = mfcv(lveg(i)) - seasf(lveg(i))*aseas(i)
        zlgveg(i) = log(zh(i)/rough(lveg(i)))
        zlgdis(i) = log((zh(i)-displa(lveg(i)))/rough(lveg(i)))
        if ( solvt > d_zero ) then
          fracd(i) = swf0(i)/solvt
        else
          fracd(i) = 0.2D0
        end if
      end do

      ! Smooth out big qs differences
      if ( nsg > 1 ) then
        call l2c_ss(lndcomm,qs,xqs)
        do i = ici1 , ici2
          do j = jci1 , jci2
            xqsdif = d_zero
            ! Coarse scale specific humidity
            xqs0 = lm%qvatm(j,i)/(d_one+lm%qvatm(j,i))
            do n = 1 , nnsg
              if ( lm%ldmsk1(n,j,i) == 1 ) then
                xqsdif = xqsdif + (xqs(n,j,i)-xqs0)
              end if
            end do
            ! Mean difference from subgrid values and coarse grid
            xqsdif = xqsdif/rdnnsg
            do n = 1 , nnsg
              xqs(n,j,i) = max(xqs(n,j,i)-xqsdif,d_zero)
            end do
          end do
        end do
        call c2l_ss(lndcomm,xqs,qs)
      end if

    else if ( ivers == 2 ) then ! bats --> regcm2d

      call fseas(tgbrd,aseas)
      if ( iemiss == 1 ) then
        do i = ilndbeg , ilndend
          fracs = lncl(i)*wt(i) + (d_one-lncl(i))*scvk(i)
          emiss(i) = (lndemiss(lveg(i))-seasemi(lveg(i))*aseas(i)) * &
                  (d_one-fracs) + 0.992*fracs
        end do
      end if
      call l2c_ss(lndcomm,tlef,lms%tlef)
      call l2c_ss(lndcomm,tgrd,lms%tgrd)
      call l2c_ss(lndcomm,tgbrd,lms%tgbrd)
      call l2c_ss(lndcomm,gwet,lms%gwet)
      call l2c_ss(lndcomm,ldew,lms%ldew)
      call l2c_ss(lndcomm,snag,lms%snag)
      call l2c_ss(lndcomm,sncv,lms%sncv)
      call l2c_ss(lndcomm,scvk,lms%scvk)
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

      do i = ilndbeg , ilndend
        fracv = sigf(i)
        fracb = (d_one-lncl(i))*(d_one-scvk(i))
        fracs = lncl(i)*wt(i) + (d_one-lncl(i))*scvk(i)
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

      if ( ichem == 1 ) then
        call l2c_ss(lndcomm,delt,lms%deltat)
        call l2c_ss(lndcomm,delq,lms%deltaq)
        call l2c_ss(lndcomm,sigf,lms%sigf)
        call l2c_ss(lndcomm,lncl,lms%lncl)
        call l2c_ss(lndcomm,wt,lms%wt)
        call l2c_ss(lndcomm,xlai,lms%xlai)
      end if
    end if ! Versus of the interface (1,2)

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine interf

  subroutine albedobats(lm,lms)
    implicit none
    type(lm_exchange) , intent(inout) :: lm
    type(lm_state) , intent(inout) :: lms
    call interf(lm,lms,1)
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

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

module mod_cu_em
  !
  ! Kerry Emanuel MIT Convective scheme
  !
  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_memutil
  use mod_runparams , only : alphae , betae , coeffr , coeffs , cu ,    &
    damp , dtmax , entp , minorig , omtrain , omtsnow , sigd , sigs ,   &
    tlcrit , iqv , ichem , clfrcv , ichcumtra , icup , dt
  use mod_runparams , only : k2_const , kfac_shal , kfac_deep
  use mod_cu_common
  use mod_service
  use mod_regcm_types

  implicit none

  private

  public :: allocate_mod_cu_em , cupemandrv

  real(rkx) , parameter :: mincbmf = 1.0e-30_rkx
  real(rkx) , parameter :: delt0 = 300.0_rkx

  ! Latent heats
  real(rkx) , parameter :: clq = 2500.0_rkx

  real(rkx) , public , pointer , dimension(:,:) :: cbmf2d
  real(rkx) , public , pointer , dimension(:,:) :: elcrit2d
  real(rkx) , public , pointer , dimension(:,:) :: epmax2d

  integer :: ncp , nap

  real(rkx) , pointer , dimension(:) :: cbmf , pret , qprime , &
    tprime , wd , elcrit , epmax
  real(rkx) , pointer , dimension(:,:) :: fq , ft , fu , fv , pcup , &
    qcup , qscup , tcup , ucup , vcup , zcup , cldfra , phcup , ppcp
  real(rkx) , pointer , dimension(:,:,:) :: ftra , tra
  integer(ik4) , pointer , dimension(:) :: iflag , kbase , ktop , &
    imap , jmap

  contains

  subroutine allocate_mod_cu_em
    implicit none

    ncp = (jci2-jci1+1) * (ici2-ici1+1)

    call getmem2d(cbmf2d,jci1,jci2,ici1,ici2,'emanuel:cbmf2d')
    call getmem2d(elcrit2d,jci1,jci2,ici1,ici2,'emanuel:elcrit2d')
    call getmem2d(epmax2d,jci1,jci2,ici1,ici2,'emanuel:epmax2d')

    call getmem1d(imap,1,ncp,'emanuel:imap')
    call getmem1d(jmap,1,ncp,'emanuel:jmap')
    call getmem1d(cbmf,1,ncp,'emanuel:cbmf')
    call getmem1d(pret,1,ncp,'emanuel:pret')
    call getmem1d(qprime,1,ncp,'emanuel:qprime')
    call getmem1d(tprime,1,ncp,'emanuel:tprime')
    call getmem1d(wd,1,ncp,'emanuel:wd')
    call getmem1d(elcrit,1,ncp,'emanuel:elcrit')
    call getmem1d(epmax,1,ncp,'emanuel:epmax')
    call getmem1d(iflag,1,ncp,'emanuel:iflag')
    call getmem1d(kbase,1,ncp,'emanuel:kbase')
    call getmem1d(ktop,1,ncp,'emanuel:ktop')
    call getmem2d(fq,1,ncp,1,kz,'emanuel:fq')
    call getmem2d(ft,1,ncp,1,kz,'emanuel:ft')
    call getmem2d(fu,1,ncp,1,kz,'emanuel:fu')
    call getmem2d(fv,1,ncp,1,kz,'emanuel:fv')
    call getmem2d(pcup,1,ncp,1,kz,'emanuel:pcup')
    call getmem2d(phcup,1,ncp,1,kzp1,'emanuel:phcup')
    call getmem2d(qcup,1,ncp,1,kz,'emanuel:qcup')
    call getmem2d(qscup,1,ncp,1,kz,'emanuel:qscup')
    call getmem2d(tcup,1,ncp,1,kz,'emanuel:tcup')
    call getmem2d(ucup,1,ncp,1,kz,'emanuel:ucup')
    call getmem2d(vcup,1,ncp,1,kz,'emanuel:vcup')
    call getmem2d(zcup,1,ncp,1,kz,'emanuel:zcup')
    call getmem2d(cldfra,1,ncp,1,kz,'emanuel:cldfra')
    call getmem2d(ppcp,1,ncp,1,kz,'emanuel:ppcp')
    if ( ichem == 1 ) then
      call getmem3d(ftra,1,ncp,1,kz,1,ntr,'emanuel:ftra')
      call getmem3d(tra,1,ncp,1,kz,1,ntr,'emanuel:tra')
    end if
  end subroutine allocate_mod_cu_em
  !
  ! **********************************************
  ! **** Driver for Emanuel Convection Scheme ****
  ! **********************************************
  !
  subroutine cupemandrv(m2c)
    implicit none
    type(mod_2_cum) , intent(in) :: m2c
    integer(ik4) :: i , j , k , n , kk
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'cupemandrv'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    nap = 0
    do i = ici1 , ici2
      do j = jci1 , jci2
        if ( cuscheme(j,i) == 4 ) then
          nap = nap + 1
          imap(nap) = i
          jmap(nap) = j
        end if
      end do
    end do

    if ( nap == 0 ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    do k = 1 , kz
      do n = 1 , nap
        i = imap(n)
        j = jmap(n)
        kk = kzp1 - k
        tcup(n,k) = m2c%tas(j,i,kk)                                   ! [K]
        ! Model wants specific humidities and pressures in mb
        qcup(n,k) = m2c%qxas(j,i,kk,iqv)/(d_one+m2c%qxas(j,i,kk,iqv)) ! [kg/kg]
        qscup(n,k) = pfqsat(m2c%tas(j,i,kk),m2c%pas(j,i,kk))
        ucup(n,k) = m2c%uas(j,i,kk)                                   ! [m/s]
        vcup(n,k) = m2c%vas(j,i,kk)                                   ! [m/s]
        zcup(n,k) = m2c%zas(j,i,kk)                                   ! [m/s]
        pcup(n,k) = m2c%pas(j,i,kk)*d_r100                            ! [hPa]
        fu(n,k) = d_zero
        fv(n,k) = d_zero
        ft(n,k) = d_zero
        fq(n,k) = d_zero
        cldfra(n,k) = d_zero
        ppcp(n,k) = d_zero
      end do
    end do

    do k = 1 , kzp1
      do n = 1 , nap
        i = imap(n)
        j = jmap(n)
        kk = kzp1 - k + 1
        phcup(n,k) = m2c%pasf(j,i,kk) * d_r100 ! [hPa]
      end do
    end do

    do n = 1 , nap
      i = imap(n)
      j = jmap(n)
      elcrit(n) = elcrit2d(j,i)
      epmax(n) = epmax2d(j,i)
      ! Past history
      cbmf(n) = cbmf2d(j,i) ! [(kg/m**2)/s]
    end do

    do n = 1 , nap
      pret(n) = d_zero
      tprime(n) = d_zero
      qprime(n) = d_zero
      iflag(n) = 0
      ktop(n) = -1
      kbase(n) = -1
    end do

    if ( ichem == 1 ) then
      do n = 1 , nap
        i = imap(n)
        j = jmap(n)
        do k = 1 , kz
          kk = kzp1 - k
          tra(n,k,:) = m2c%chias(j,i,kk,:)   ! [kg/kg]
          ftra(n,k,:) = d_zero
        end do
      end do
    end if

    call cupeman(nap,kz,ntr,tcup,qcup,qscup,ucup,vcup,tra,pcup,phcup, &
                 iflag,ft,fq,fu,fv,ftra,pret,ppcp,wd,tprime,qprime,   &
                 cbmf,cldfra,kbase,ktop,elcrit,epmax)

    do n = 1 , nap
      i = imap(n)
      j = jmap(n)
      cbmf2d(j,i) = cbmf(n)
    end do

    do n = 1 , nap
      ! iflag=0: No moist convection; atmosphere stable or surface
      !          temperature < 250K or surface humidity is negative.
      ! iflag=1: Moist convection occurs.
      ! iflag=2: No moist convection: lifted condensation level above 200 mb.
      ! iflag=3: No moist convection: cloud base higher than level kzm2.
      ! iflag=4: Moist convection occurs, but CFL condition on the
      !          subsidence warming is violated. (Does not terminate scheme.)
      if ( iflag(n) == 1 .or. iflag(n) == 4 ) then
        i = imap(n)
        j = jmap(n)
        ! Tendencies
        do k = 1 , kz
          kk = kzp1 - k
          cu_tten(j,i,kk) = ft(n,k)
          ! Move specific humidity tendency to mixing ratio tendency
          cu_qten(j,i,kk,iqv) = fq(n,k)/(d_one-qcup(n,k))**2
          ! There is a bit of an inconsistency here...  The wind
          ! tendencies from convection are on cross points, but the
          ! model wants them on dot points.
          cu_uten(j,i,kk) = fu(n,k)
          cu_vten(j,i,kk) = fv(n,k)
        end do

        if ( ktop(n) > 0 .and. kbase(n) > 0 ) then
          cu_ktop(j,i) = kzp1-ktop(n)
          cu_kbot(j,i) = kzp1-kbase(n)
        end if
        do k = 1 , kz
          kk = kzp1 - k
          cu_cldfrc(j,i,k) = cldfra(n,kk)
        end do

        ! Precipitation
        if ( pret(n) > dlowval ) then
          ! The order top/bottom for regcm is reversed.
          cu_prate(j,i) = cu_prate(j,i) + pret(n)
          total_precip_points = total_precip_points + 1
        end if
      end if
    end do

    ! Tracer tendency
    if ( ichem == 1 .and. ichcumtra == 1 .and. &
         (.not. any( icup == 2 .or. icup == 6 )) ) then
      do n = 1 , nap
        if ( iflag(n) == 1 .or. iflag(n) == 4 ) then
          i = imap(n)
          j = jmap(n)
          do k = 1 , kz
            kk = kzp1 - k
            cu_chiten(j,i,kk,:) = ftra(n,k,:)
          end do
          ! Build for chemistry 3d table of constant precipitation rate
          ! from the surface to the top of the convection
          do k = 1 , ktop(n)-1
            kk = kzp1 - k
            cu_convpr(j,i,kk) = ppcp(n,k)
          end do
        end if
      end do
    end if

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
    contains
#include <pfesat.inc>
#include <pfqsat.inc>
  end subroutine cupemandrv
!
!**************************************************************************
!****                    subroutine cupeman (formerly convect)        *****
!****                          version 4.3C                           *****
!****                          20 May, 2002                           *****
!****                          Kerry Emanuel                          *****
!**************************************************************************
!
!-----------------------------------------------------------------------------
!    *** on input:      ***
!
!     t:   array of absolute temperature (k) of dimension nd, with first
!           index corresponding to lowest model level.
!
!     q:   array of specific humidity (gm/gm) of dimension nd, with first
!            index corresponding to lowest model level. must be defined
!            at same grid levels as t.
!
!     qs:  array of saturation specific humidity of dimension nd, with first
!            index corresponding to lowest model level. must be defined
!            at same grid levels as t.
!
!     u:   array of zonal wind velocity (m/s) of dimension nd, witth first
!            index corresponding with the lowest model level. defined at
!            same levels as t.
!
!     v:   same as u but for meridional velocity.
!
!     tra: array of passive tracer mixing ratio, of dimensions (nd,ntra),
!            where ntra is the number of different tracers. if no
!            convective tracer transport is needed, define a dummy
!            input array of dimension (nd,1). tracers are defined at
!            same vertical levels as t.
!
!     p:   array of pressure (mb) of dimension nd, with first
!            index corresponding to lowest model level. must be defined
!            at same grid levels as t.
!
!     ph:  array of pressure (mb) of dimension nd+1, with first index
!            corresponding to lowest level. these pressures are defined at
!            levels intermediate between those of p, t, q and qs. the first
!            value of ph should be greater than (i.e. at a lower level than)
!            the first value of the array p.
!
!     nd:  the dimension of the arrays t,q,qs,p,ph,ft and fq
!
!     ntra:the number of different tracers. if no tracer transport
!            is needed, set this equal to 1. (on most compilers, setting
!            ntra to 0 will bypass tracer calculation, saving some cpu.)
!
!----------------------------------------------------------------------------
!    ***   on output:         ***
!
!     iflag: an output integer whose value denotes the following:
!
!                value                        interpretation
!                -----                        --------------
!                  0               no moist convection; atmosphere is not
!                                  unstable, or surface temperature is less
!                                  than 250 k or surface specific humidity
!                                  is non-positive.
!
!                  1               moist convection occurs.
!
!                  2               no moist convection: lifted condensation
!                                  level is above the 200 mb level.
!
!                  3               no moist convection: cloud base is higher
!                                  then the level nl-1.
!
!                  4               moist convection occurs, but a cfl condition
!                                  on the subsidence warming is violated. this
!                                  does not cause mod_the scheme to terminate.
!
!     ft:   array of temperature tendency (k/s) of dimension nd, defined at same
!             grid levels as t, q, qs and p.
!
!     fq:   array of specific humidity tendencies ((gm/gm)/s) of dimension nd,
!             defined at same grid levels as t, q, qs and p.
!
!     fu:   array of forcing of zonal velocity (m/s^2) of dimension nd,
!             defined at same grid levels as t.
!
!     fv:   same as fu, but for forcing of meridional velocity.
!
!     ftra: array of forcing of tracer content, in tracer mixing ratio per
!             second, defined at same levels as t. dimensioned (nd,ntra).
!
!     precip: scalar convective precipitation rate (mm/s).
!
!     wd:    a convective downdraft velocity scale. for use in surface
!             flux parameterizations. see convect.ps file for details.
!
!     tprime: a convective downdraft temperature perturbation scale (k).
!              for use in surface flux parameterizations. see convect.ps
!              file for details.
!
!     qprime: a convective downdraft specific humidity
!              perturbation scale (gm/gm).
!              for use in surface flux parameterizations. see convect.ps
!              file for details.
!
!     cbmf:   the cloud base mass flux ((kg/m**2)/s). this scalar value must
!              be stored by the calling program and returned to convect at
!              its next call. that is, the value of cbmf must be "remembered"
!              by the calling program between calls to convect.
!
!------------------------------------------------------------------------------
!
!    ***  the parameter na should in general be greater than   ***
!    ***                or equal to  nd + 1                    ***
!
!------------------------------------------------------------------------------
! modifications for regcm:
!   1. units for precipitation were change from mm/day to mm/s
!   2. the thermodynamic constants were made consistent with those
!      of regcm.
!   3. dependence of latent heat of vaporization on temperature
!      removed because it is neglected in regcm. this is done by
!      setting cpv equal to cl and setting wlhv to the regcm value.
!   4. added cloud base (icb) and cloud top (ict) to the output to
!      compute the cloud fraction and cloud liquid water content.
!   5. each variable is now explicitly declared.  that is, the
!      "implicit none" option was added.
!   6. the value minorig is increased because the thickness of the
!      lowest layer(s) is(are) too small. if the thickness is too
!      small for the given timestep, the mass of the layer is likely
!      to be evacuated.
!   7. a maximum value to the cloud base mass flux has been added.
!   8. Cloud fraction computation from massflux has been added.
!   9. Remove dry adiabatic adjustment: regcm has a pbl scheme.
!
  subroutine cupeman(np,nd,ntra,t,q,qs,u,v,tra,p,ph,iflag,ft,fq,fu,fv,  &
                     ftra,precip,ppcp,wd,tprime,qprime,cbmf,cldfra,kcb, &
                     kct,elcrit,epmax)
    implicit none
    integer , intent(in) :: np , nd , ntra
    real(rkx) , pointer , dimension(:) , intent(inout) :: cbmf
    real(rkx) , pointer , dimension(:) , intent(in) :: elcrit , epmax
    real(rkx) , pointer , dimension(:,:) , intent(in) :: p , q , qs , t , u , v
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: tra
    real(rkx) , pointer , dimension(:,:) , intent(in) :: ph
    real(rkx) , pointer , dimension(:) , intent(inout) :: precip
    real(rkx) , pointer , dimension(:) , intent(inout) :: qprime , tprime , wd
    integer(ik4) , pointer , dimension(:) , intent(inout) :: iflag , kcb , kct
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: fq , ft , fu , fv
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: cldfra , ppcp
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: ftra

    real(rkx) :: ad , afac , ahmax , ahmin , alt , altem , am ,          &
                 anum , asij , awat , b6 , bf2 , bsum , by , byp , c6 ,  &
                 cape , capem , cbmfold , chi , coeff , cpinv , cwat ,   &
                 damps , dbo , dbosum , defrac , dei , delm , delp ,     &
                 denom , dhdp , dpinv , dtma , dtmnx , dtpbl , elacrit , &
                 ents , fac , fqold , frac , ftold , ftraold , fuold ,   &
                 fvold , plcl , qp1 , qsm , qstm , qti , rat , revap ,   &
                 rh , scrit , sigt , sjmax , sjmin , smid , smin ,       &
                 stemp , tca , traav , tvaplcl , tvpplcl , tvx ,  tvy ,  &
                 uav , vav , wdtrain , amp1
    real(rkx) , dimension(nd+1) :: clw , cpn , ep , evap , gz , h , hm , &
                                   hp , lv , lvcp , m , mp , qp , sigp , &
                                   tp , tv , tvp , up , vp , water , wt
    integer(ik4) :: n , nl , i , ihmin , icb , ict , ict1 , j , jtt , k , nk
    integer(ik4) , dimension(nd+1) :: nent
    real(rkx) , dimension(nd+1,nd+1) :: elij , ment , qent , sij , uent , vent
    real(rkx) , dimension(nd+1,nd+1,ntra) :: traent
    real(rkx) , dimension(nd+1,ntra) :: trap
    logical :: chemcutran
    !
    !   specify switches
    !
    !   minorig: lowest level from which convection may originate
    !            (should be first model level at which t is defined
    !            for models using bulk pbl schemes; otherwise, it should
    !            be the first model level at which t is defined above
    !            the surface layer)
    !

    nl = nd - 1

    chemcutran = ( ichem == 1 .and. ichcumtra > 0 )

    pointloop: &
    do  n = 1 , np
      !
      ! Calculate arrays of geopotential, heat capacity and static energy
      !
      gz(1) = d_zero
      cpn(1) = cpd*(d_one-q(n,1)) + cpv*q(n,1)
      h(1) = t(n,1)*cpn(1)
      lv(1) = wlh(t(n,1))
      hm(1) = lv(1)*q(n,1)
      tv(1) = t(n,1)*(d_one+q(n,1)*rgowi-q(n,1))
      ahmin = 1.0e12_rkx
      ihmin = nl
      do i = 2 , nl + 1
        tvx = t(n,i)*(d_one+q(n,i)*rgowi-q(n,i))
        tvy = t(n,i-1)*(d_one+q(n,i-1)*rgowi-q(n,i-1))
        gz(i) = gz(i-1) + (rgas*d_half)*(tvx+tvy)*(p(n,i-1)-p(n,i))/ph(n,i)
        cpn(i) = cpd*(d_one-q(n,i)) + cpv*q(n,i)
        h(i) = t(n,i)*cpn(i) + gz(i)
        lv(i) = wlh(t(n,i))
        hm(i) = (cpd*(d_one-q(n,i))+cpv*q(n,i))*(t(n,i)-t(n,1))+ &
                 lv(i)*q(n,i)+gz(i)
        tv(i) = t(n,i)*(d_one+q(n,i)*rgowi-q(n,i))
        !
        ! Find level of minimum moist static energy
        !
        if ( i >= minorig .and. hm(i) < ahmin .and. hm(i) < hm(i-1) ) then
          ahmin = hm(i)
          ihmin = i
        end if
      end do
      ihmin = min(ihmin,nl-1)
      !
      ! Find the model level below the level of minimum moist
      ! static energy that has the maximum value of moist static energy
      !
      ahmax = d_zero
      nk = minorig
      do i = minorig , ihmin
        if ( hm(i) > ahmax ) then
          nk = i
          ahmax = hm(i)
        end if
      end do
      !
      ! Check whether parcel level temperature and specific humidity are
      ! reasonable. Skip convection if hm increases monotonically upward
      !
      if ( t(n,nk) < 250.0_rkx .or. q(n,nk) <= d_zero .or. &
        ihmin == (nl-1) ) then
        iflag(n) = 0
        cbmf(n) = d_zero
        cycle pointloop
      end if
      !
      ! Calculate lifted condensation level of air at parcel origin level
      ! (within 0.2% of formula of bolton, mon. wea. rev.,1980)
      !
      rh = max(min(q(n,nk)/qs(n,nk),1.1_rkx),0.1_rkx)
      chi = t(n,nk)/(1669.0_rkx-122.0_rkx*rh-t(n,nk))
      plcl = p(n,nk)*(rh**chi)
      if ( plcl < 200.0_rkx .or. plcl >= 2000.0_rkx ) then
        iflag(n) = 2
        cbmf(n) = d_zero
        cycle pointloop
      end if
      !
      ! Calculate first level above lcl (=icb)
      !
      icb = nl - 1
      do i = nk + 1 , nl
        if ( p(n,i) < plcl ) icb = min(icb,i)
      end do
      if ( icb >= (nl-1) ) then
        iflag(n) = 3
        cbmf(n) = d_zero
        cycle pointloop
      end if
      !
      ! Find temperature up through icb and test for instability
      !
      ! Subroutine tlift calculates part of the lifted parcel virtual
      ! temperature, the actual temperature and the adiabatic
      ! liquid water content
      !
      call tlift(n,p,t,q,qs,gz,icb,nk,tvp,tp,clw,nd,nl,1)
      do i = nk , icb
        tvp(i) = tvp(i) - tp(i)*q(n,nk)
      end do
      !
      ! If there was no convection at last time step and parcel
      ! is stable at icb then skip rest of calculation
      !
      if ( abs(cbmf(n)) < mincbmf .and. tvp(icb) <= (tv(icb)-dtmax) ) then
        iflag(n) = 0
        cycle pointloop
      end if
      !
      ! if this point is reached, moist convective adjustment is necessary
      !
      if ( iflag(n) /= 4 ) iflag = 1
      !
      ! Find the rest of the lifted parcel temperatures
      !
      call tlift(n,p,t,q,qs,gz,icb,nk,tvp,tp,clw,nd,nl,2)
      !
      ! Set the precipitation efficiencies and the fraction of
      ! precipitation falling outside of cloud
      ! these may be functions of tp(i), p(i) and clw(i)
      !
      do i = 1 , nk
        ep(i) = d_zero
        sigp(i) = sigs
      end do
      do i = nk + 1 , nl
        tca = tp(i) - tzero
        if ( tca >= d_zero ) then
          elacrit = elcrit(n)
        else
          elacrit = elcrit(n)*(d_one-tca/tlcrit)
        end if
        elacrit = max(elacrit,d_zero)
        ep(i) = epmax(n)*(d_one-elacrit/max(clw(i),1.0e-8_rkx))
        ep(i) = max(ep(i),d_zero)
        ep(i) = min(ep(i),epmax(n))
        sigp(i) = sigs
      end do
      !
      ! Calculate virtual temperature and lifted parcel
      ! virtual temperature
      !
      do i = icb + 1 , nl
        tvp(i) = tvp(i) - tp(i)*q(n,nk)
      end do
      tvp(nl+1) = tvp(nl) - (gz(nl+1)-gz(nl))*rcpd
      !
      ! now initialize various arrays used in the computations
      !
      do i = 1 , nl + 1
        hp(i) = h(i)
        nent(i) = 0
        water(i) = d_zero
        evap(i) = d_zero
        wt(i) = omtsnow
        mp(i) = d_zero
        m(i) = d_zero
        lvcp(i) = lv(i)/cpn(i)
        do j = 1 , nl + 1
          qent(i,j) = q(n,j)
          elij(i,j) = d_zero
          ment(i,j) = d_zero
          sij(i,j) = d_zero
          uent(i,j) = u(n,j)
          vent(i,j) = v(n,j)
          if ( chemcutran ) then
            do k = 1 , ntra
              traent(i,j,k) = tra(n,j,k)
            end do
          end if
        end do
      end do
      qp(1) = q(n,1)
      up(1) = u(n,1)
      vp(1) = v(n,1)
      if ( chemcutran ) then
        do i = 1 , ntra
          trap(1,i) = tra(n,1,i)
        end do
      end if
      do i = 2 , nl + 1
        qp(i) = q(n,i-1)
        up(i) = u(n,i-1)
        vp(i) = v(n,i-1)
        if ( chemcutran ) then
          do j = 1 , ntra
            trap(i,j) = tra(n,i-1,j)
          end do
        end if
      end do
      !
      ! Find the first model level (ict1) above the parcel's highest level
      ! of neutral buoyancy and the highest level of positive cape (ict)
      !
      cape = d_zero
      capem = d_zero
      ict = icb + 1
      ict1 = ict
      byp = d_zero
      do i = icb + 1 , nl - 1
        by = (tvp(i)-tv(i))*(ph(n,i)-ph(n,i+1))/p(n,i)
        cape = cape + by
        if ( by >= d_zero ) ict1 = i + 1
        if ( cape > d_zero ) then
          ict = i + 1
          byp = (tvp(i+1)-tv(i+1))*(ph(n,i+1)-ph(n,i+2))/p(n,i+1)
          capem = cape
        end if
      end do
      ict = max(ict,ict1)
      cape = capem + byp
      defrac = capem - cape
      defrac = max(defrac,0.001_rkx)
      frac = -cape/defrac
      frac = min(frac,d_one)
      frac = max(frac,d_zero)
      !
      ! Calculate liquid water static energy of lifted parcel
      !
      do i = icb , ict
        hp(i) = h(nk) + (lv(i)+(cpd-cpw)*t(n,i))*ep(i)*clw(i)
      end do
      !
      ! Calculate cloud base mass flux and rates of mixing, m(i),
      ! at each model level
      !
      dbosum = d_zero
      !
      ! Interpolate difference between lifted parcel and
      ! environmental temperatures to lifted condensation level
      !
      tvpplcl = tvp(icb-1) - rgas*tvp(icb-1)*(p(n,icb-1)-plcl) / &
                                             (cpn(icb-1)*p(n,icb-1))
      tvaplcl = tv(icb)+(tvp(icb)-tvp(icb+1))*(plcl-p(n,icb)) / &
                                             (p(n,icb)-p(n,icb+1))
      dtpbl = d_zero
      do i = nk , icb - 1
        dtpbl = dtpbl + (tvp(i)-tv(i))*(ph(n,i)-ph(n,i+1))
      end do
      dtpbl = dtpbl/(ph(n,nk)-ph(n,icb))
      dtmnx = tvpplcl - tvaplcl + dtmax + dtpbl
      dtma = dtmnx
      !
      ! Adjust cloud base mass flux
      !
      cbmfold = cbmf(n)
      damps = damp*dt/delt0
      cbmf(n) = (d_one-damps)*cbmf(n) + 0.1_rkx*alphae*dtma
      cbmf(n) = max(cbmf(n),d_zero)
      !
      ! If cloud base mass flux is zero, skip rest of calculation
      !
      if ( abs(cbmf(n)) < mincbmf .and. abs(cbmfold) < mincbmf ) then
        iflag(n) = 0
        cycle pointloop
      end if
      !
      ! calculate rates of mixing, m(i)
      !
      m(icb) = d_zero
      do i = icb + 1 , ict
        k = min(i,ict1)
        dbo = abs(tv(k)-tvp(k)) + entp*0.02_rkx*(ph(n,k)-ph(n,k+1))
        dbosum = dbosum + dbo
        m(i) = cbmf(n)*dbo
      end do
      do i = icb + 1 , ict
        m(i) = m(i)/dbosum
      end do
      !
      ! Calculate entrained air mass flux (ment),
      ! total water mixing ratio (qent),
      ! total condensed water (elij), and mixing fraction (sij)
      !
      do i = icb + 1 , ict
        qti = q(n,nk) - ep(i)*clw(i)
        do j = icb , ict
          bf2 = d_one + lv(j)*lv(j)*qs(n,j)/(rwat*t(n,j)*t(n,j)*cpd)
          anum = h(j) - hp(i) + (cpv-cpd)*t(n,j)*(qti-q(n,j))
          denom = h(i) - hp(i) + (cpd-cpv)*(q(n,i)-qti)*t(n,j)
          dei = denom
          if ( abs(dei) < 0.01_rkx ) dei = 0.01_rkx
          sij(i,j) = anum/dei
          sij(i,i) = d_one
          altem = sij(i,j)*q(n,i) + (d_one-sij(i,j))*qti - qs(n,j)
          altem = altem/bf2
          cwat = clw(j)*(d_one-ep(j))
          stemp = sij(i,j)
          if ( (stemp < d_zero .or. stemp > d_one .or. &
                altem > cwat) .and. j > i ) then
            anum = anum - lv(j)*(qti-qs(n,j)-cwat*bf2)
            denom = denom + lv(j)*(q(n,i)-qti)
            if ( abs(denom) < 0.01_rkx ) denom = 0.01_rkx
            sij(i,j) = anum/denom
            altem = sij(i,j)*q(n,i) + (d_one-sij(i,j))*qti - qs(n,j)
            altem = altem - (bf2-d_one)*cwat
          end if
          if ( sij(i,j) > d_zero .and. sij(i,j) < 0.9_rkx ) then
            qent(i,j) = sij(i,j)*q(n,i) + (d_one-sij(i,j))*qti
            uent(i,j) = sij(i,j)*u(n,i) + (d_one-sij(i,j))*u(n,nk)
            vent(i,j) = sij(i,j)*v(n,i) + (d_one-sij(i,j))*v(n,nk)
            if ( chemcutran ) then
              do k = 1 , ntra
                traent(i,j,k) = sij(i,j)*tra(n,i,k)+(d_one-sij(i,j))*tra(n,nk,k)
              end do
            end if
            elij(i,j) = altem
            elij(i,j) = max(d_zero,elij(i,j))
            ment(i,j) = m(i)/(d_one-sij(i,j))
            nent(i) = nent(i) + 1
          end if
          sij(i,j) = max(d_zero,sij(i,j))
          sij(i,j) = min(d_one,sij(i,j))
        end do
        !
        ! If no air can entrain at level i assume that updraft detrains
        ! at that level and calculate detrained air flux and properties
        !
        if ( nent(i) == 0 ) then
          ment(i,i) = m(i)
          qent(i,i) = q(n,nk) - ep(i)*clw(i)
          uent(i,i) = u(n,nk)
          vent(i,i) = v(n,nk)
          if ( chemcutran ) then
            do j = 1 , ntra
              traent(i,i,j) = tra(n,nk,j)
            end do
          end if
          elij(i,i) = clw(i)
          sij(i,i) = d_one
        end if
      end do
      sij(ict,ict) = d_one
      !
      ! Normalize entrained air mass fluxes to represent equal
      ! probabilities of mixing
      !
      do i = icb + 1 , ict
        if ( nent(i) /= 0 ) then
          qp1 = q(n,nk) - ep(i)*clw(i)
          anum = h(i) - hp(i) - lv(i)*(qp1-qs(n,i))
          denom = h(i) - hp(i) + lv(i)*(q(n,i)-qp1)
          if ( abs(denom) < 0.01_rkx ) denom = 0.01_rkx
          scrit = anum/denom
          alt = qp1 - qs(n,i) + scrit*(q(n,i)-qp1)
          if ( alt < d_zero ) scrit = d_one
          scrit = max(scrit,d_zero)
          asij = d_zero
          smin = d_one
          do j = icb , ict
            if ( sij(i,j) > d_zero .and. sij(i,j) < 0.9_rkx ) then
              if ( j > i ) then
                smid = min(sij(i,j),scrit)
                sjmax = smid
                sjmin = smid
                if ( smid < smin .and. sij(i,j+1) < smid ) then
                  smin = smid
                  sjmax = min(sij(i,j+1),sij(i,j),scrit)
                  sjmin = max(sij(i,j-1),sij(i,j))
                  sjmin = min(sjmin,scrit)
                end if
              else
                sjmax = max(sij(i,j+1),scrit)
                smid = max(sij(i,j),scrit)
                sjmin = d_zero
                if ( j > 1 ) sjmin = sij(i,j-1)
                sjmin = max(sjmin,scrit)
              end if
              delp = abs(sjmax-smid)
              delm = abs(sjmin-smid)
              asij = asij + (delp+delm)*(ph(n,j)-ph(n,j+1))
              ment(i,j) = ment(i,j)*(delp+delm)*(ph(n,j)-ph(n,j+1))
            end if
          end do
          asij = max(1.0e-21_rkx,asij)
          asij = d_one/asij
          do j = icb , ict
            ment(i,j) = ment(i,j)*asij
          end do
          bsum = d_zero
          do j = icb , ict
            bsum = bsum + ment(i,j)
          end do
          if ( bsum < 1.0e-18_rkx ) then
            nent(i) = 0
            ment(i,i) = m(i)
            qent(i,i) = q(n,nk) - ep(i)*clw(i)
            uent(i,i) = u(n,nk)
            vent(i,i) = v(n,nk)
            if ( chemcutran ) then
              do j = 1 , ntra
                traent(i,i,j) = tra(n,nk,j)
              end do
            end if
            elij(i,i) = clw(i)
            sij(i,i) = d_one
          end if
        end if
      end do
      !
      ! Check whether ep(ict)=0, if so, skip precipitating
      ! downdraft calculation
      !
      if ( ep(ict) >= 0.0001_rkx ) then
        !
        ! Integrate liquid water equation to find condensed water
        ! and condensed water flux
        !
        jtt = 2
        !
        ! Begin downdraft loop
        !
        do i = ict , 1 , -1
          !
          ! Calculate detrained precipitation
          !
          wdtrain = egrav*ep(i)*m(i)*clw(i)
          if ( i > 1 ) then
            do j = 1 , i - 1
              awat = elij(j,i) - (d_one-ep(i))*clw(i)
              awat = max(d_zero,awat)
              wdtrain = wdtrain + egrav*awat*ment(j,i)
            end do
          end if
          !
          ! Find rain water and evaporation using provisional
          ! estimates of qp(i)and qp(i-1)
          ! Value of terminal velocity and coefficient of evaporation for snow
          !
          !
          ! Value of terminal velocity and coefficient of evaporation for rain
          !
          if ( t(n,i) > tzero) then
            coeff = coeffr
            wt(i) = omtrain
          else
            coeff = coeffs
            wt(i) = omtsnow
          end if
          qsm = d_half*(q(n,i)+qp(i+1))
          afac = coeff*ph(n,i)*(qs(n,i)-qsm) / &
                      (1.0e4_rkx+2.0e3_rkx*ph(n,i)*qs(n,i))
          afac = max(afac,d_zero)
          sigt = sigp(i)
          sigt = max(d_zero,sigt)
          sigt = min(d_one,sigt)
          b6 = d_100*(ph(n,i)-ph(n,i+1))*sigt*afac/wt(i)
          c6 = (water(i+1)*wt(i+1)+wdtrain/sigd)/wt(i)
          revap = d_half*(-b6+sqrt(b6*b6+d_four*c6))
          evap(i) = sigt*afac*revap
          water(i) = revap*revap
          ppcp(n,i) = wt(i)*sigd*water(i)*regrav ! mm/s
          !
          ! Calculate precipitating downdraft mass flux under
          ! hydrostatic approximation
          !
          if ( i /= 1 ) then
            dhdp = (h(i)-h(i-1))/(p(n,i-1)-p(n,i))
            dhdp = max(dhdp,d_10)
            mp(i) = d_100*regrav*lv(i)*sigd*evap(i)/dhdp
            mp(i) = max(mp(i),d_zero)
            !
            ! Add small amount of inertia to downdraft
            !
            fac = 20.0_rkx/(ph(n,i-1)-ph(n,i))
            mp(i) = (fac*mp(i+1)+mp(i))/(d_one+fac)
            !
            ! Force mp to decrease linearly to zero
            ! between about 950 mb and the surface
            !
            if ( p(n,i) > (0.949_rkx*p(n,1)) ) then
              jtt = max(jtt,i)
              mp(i) = mp(jtt)*(p(n,1)-p(n,i))/(p(n,1)-p(n,jtt))
            end if
          end if
          !
          ! Find mixing ratio of precipitating downdraft
          !
          if ( i /= ict ) then
            if ( i == 1 ) then
              qstm = qs(n,1)
            else
              qstm = qs(n,i-1)
            end if
            if ( mp(i) > mp(i+1) ) then
              rat = mp(i+1)/mp(i)
              qp(i) = qp(i+1)*rat + q(n,i)*(d_one-rat) + &
                    d_100*regrav*sigd*(ph(n,i)-ph(n,i+1))*(evap(i)/mp(i))
              up(i) = up(i+1)*rat + u(n,i)*(d_one-rat)
              vp(i) = vp(i+1)*rat + v(n,i)*(d_one-rat)
              if ( chemcutran ) then
                do j = 1 , ntra
                  trap(i,j) = trap(i+1,j)*rat + trap(i,j)*(d_one-rat)
                end do
              end if
            else if ( mp(i+1) > d_zero ) then
              qp(i) = (gz(i+1)-gz(i)+qp(i+1)*(lv(i+1)+t(n,i+1)*(cpv-cpd)) + &
                       cpd*(t(n,i+1)-t(n,i)))/(lv(i)+t(n,i)*(cpv-cpd))
              up(i) = up(i+1)
              vp(i) = vp(i+1)
              if ( chemcutran ) then
                do j = 1 , ntra
                  trap(i,j) = trap(i+1,j)
                end do
              end if
            end if
            qp(i) = min(qp(i),qstm)
            qp(i) = max(qp(i),d_zero)
          end if
        end do
        !
        ! Calculate surface precipitation in mm/s
        !
        precip(n) = precip(n) + wt(1)*sigd*water(1)*regrav ! mm/s
      end if
      !
      ! Calculate downdraft velocity scale and surface temperature and
      ! water vapor fluctuations
      !
      wd = betae*abs(mp(icb))*0.01_rkx*rgas*t(n,icb)/(sigd*p(n,icb))
      qprime = d_half*(qp(1)-q(n,1))
      tprime = wlhv*qprime*rcpd
      !
      ! Calculate tendencies of lowest level potential temperature
      ! and mixing ratio
      !
      dpinv = 0.01_rkx/(ph(n,1)-ph(n,2))
      am = d_zero
      if ( nk == 1 ) then
        do k = 2 , ict
          am = am + m(k)
        end do
      end if
      if ( (d_two*egrav*dpinv*am) >= d_one/dt ) iflag = 4
      ft(n,1) = ft(n,1) + egrav*dpinv*am*(t(n,2)-t(n,1)+(gz(2)-gz(1))/cpn(1))
      ft(n,1) = ft(n,1) - lvcp(1)*sigd*evap(1)
      ft(n,1) = ft(n,1) + sigd*wt(2)*(cpw-cpd)*water(2)* &
                (t(n,2)-t(n,1))*dpinv/cpn(1)
      fq(n,1) = fq(n,1) + egrav*mp(2)*(qp(2)-q(n,1))*dpinv + sigd*evap(1)
      fq(n,1) = fq(n,1) + egrav*am*(q(n,2)-q(n,1))*dpinv
      fu(n,1) = fu(n,1) + egrav*dpinv*(mp(2)*(up(2)-u(n,1))+am*(u(n,2)-u(n,1)))
      fv(n,1) = fv(n,1) + egrav*dpinv*(mp(2)*(vp(2)-v(n,1))+am*(v(n,2)-v(n,1)))
      if ( chemcutran ) then
        do j = 1 , ntra
          ftra(n,1,j) = ftra(n,1,j) + egrav*dpinv * &
                 (mp(2)*(trap(2,j)-tra(n,1,j)) + am*(tra(n,2,j)-tra(n,1,j)))
        end do
      end if
      do j = 2 , ict
        fq(n,1) = fq(n,1) + egrav*dpinv*ment(j,1)*(qent(j,1)-q(n,1))
        fu(n,1) = fu(n,1) + egrav*dpinv*ment(j,1)*(uent(j,1)-u(n,1))
        fv(n,1) = fv(n,1) + egrav*dpinv*ment(j,1)*(vent(j,1)-v(n,1))
        if ( chemcutran ) then
          do k = 1 , ntra
            ftra(n,1,k) = ftra(n,1,k) + &
                    egrav*dpinv*ment(j,1)*(traent(j,1,k)-tra(n,1,k))
          end do
        end if
      end do
      !
      ! Calculate tendencies of potential temperature and mixing ratio
      ! at levels above the lowest level
      ! First find the net saturated updraft and downdraft mass fluxes
      ! through each level
      !
      do i = 2 , ict
        dpinv = 0.01_rkx/(ph(n,i)-ph(n,i+1))
        cpinv = d_one/cpn(i)
        ad = d_zero
        amp1 = d_zero
        if ( i >= nk ) then
          do k = i + 1 , ict + 1
            amp1 = amp1 + m(k)
          end do
        end if
        do k = 1 , i
          do j = i + 1 , ict + 1
            amp1 = amp1 + ment(k,j)
          end do
        end do
        if ( (d_two*egrav*dpinv*amp1) >= d_one/dt ) iflag = 4
        do k = 1 , i - 1
          do j = i , ict
            ad = ad + ment(j,k)
          end do
        end do
        ft(n,i) = ft(n,i) + egrav*dpinv*(amp1*(t(n,i+1)-t(n,i)+ &
              (gz(i+1)-gz(i))*cpinv)-ad*(t(n,i)-t(n,i-1)+ &
              (gz(i)-gz(i-1))*cpinv)) - sigd*lvcp(i)*evap(i)
        ft(n,i) = ft(n,i) + egrav*dpinv*ment(i,i) * &
              (hp(i)-h(i)+t(n,i)*(cpv-cpd)*(q(n,i)-qent(i,i)))*cpinv
        ft(n,i) = ft(n,i) + sigd*wt(i+1)*(cpw-cpd)*water(i+1) * &
              (t(n,i+1)-t(n,i))*dpinv*cpinv
        fq(n,i) = fq(n,i) + egrav*dpinv * &
              (amp1*(q(n,i+1)-q(n,i))-ad*(q(n,i)-q(n,i-1)))
        fu(n,i) = fu(n,i) + egrav*dpinv * &
              (amp1*(u(n,i+1)-u(n,i))-ad*(u(n,i)-u(n,i-1)))
        fv(n,i) = fv(n,i) + egrav*dpinv * &
              (amp1*(v(n,i+1)-v(n,i))-ad*(v(n,i)-v(n,i-1)))
        if ( chemcutran ) then
          do k = 1 , ntra
            ftra(n,i,k) = ftra(n,i,k) + egrav*dpinv * &
               (amp1*(tra(n,i+1,k)-tra(n,i,k))-ad*(tra(n,i,k)-tra(n,i-1,k)))
          end do
        end if
        do k = 1 , i - 1
          awat = elij(k,i) - (d_one-ep(i))*clw(i)
          awat = max(awat,d_zero)
          fq(n,i) = fq(n,i) + egrav*dpinv*ment(k,i)*(qent(k,i)-awat-q(n,i))
          fu(n,i) = fu(n,i) + egrav*dpinv*ment(k,i)*(uent(k,i)-u(n,i))
          fv(n,i) = fv(n,i) + egrav*dpinv*ment(k,i)*(vent(k,i)-v(n,i))
          if ( chemcutran ) then
            do j = 1 , ntra
              ftra(n,i,j) = ftra(n,i,j) + &
                egrav*dpinv*ment(k,i)*(traent(k,i,j)-tra(n,i,j))
            end do
          end if
        end do
        do k = i , ict
          fq(n,i) = fq(n,i) + egrav*dpinv*ment(k,i)*(qent(k,i)-q(n,i))
          fu(n,i) = fu(n,i) + egrav*dpinv*ment(k,i)*(uent(k,i)-u(n,i))
          fv(n,i) = fv(n,i) + egrav*dpinv*ment(k,i)*(vent(k,i)-v(n,i))
          if ( chemcutran ) then
            do j = 1 , ntra
              ftra(n,i,j) = ftra(n,i,j) + &
                egrav*dpinv*ment(k,i)*(traent(k,i,j)-tra(n,i,j))
            end do
          end if
        end do
        fq(n,i) = fq(n,i) + sigd*evap(i) + egrav * &
                (mp(i+1)*(qp(i+1)-q(n,i)) - mp(i)*(qp(i)-q(n,i-1)))*dpinv
        fu(n,i) = fu(n,i) + egrav * &
                (mp(i+1)*(up(i+1)-u(n,i)) - mp(i)*(up(i)-u(n,i-1)))*dpinv
        fv(n,i) = fv(n,i) + egrav * &
                (mp(i+1)*(vp(i+1)-v(n,i)) - mp(i)*(vp(i)-v(n,i-1)))*dpinv
        if ( chemcutran ) then
          do j = 1 , ntra
            ftra(n,i,j) = ftra(n,i,j) + egrav*dpinv * &
             (mp(i+1)*(trap(i+1,j)-tra(n,i,j)) - mp(i)*(trap(i,j)-tra(n,i-1,j)))
          end do
        end if
      end do
      !
      ! Adjust tendencies at top of convection layer to reflect
      ! actual position of the level zero cape
      !
      fqold = fq(n,ict)
      fq(n,ict) = fq(n,ict)*(d_one-frac)
      fq(n,ict-1) = fq(n,ict-1) + frac*fqold * &
                  ((ph(n,ict)-ph(n,ict+1))/(ph(n,ict-1)-ph(n,ict))) * &
                  lv(ict)/lv(ict-1)
      ftold = ft(n,ict)
      ft(n,ict) = ft(n,ict)*(d_one-frac)
      ft(n,ict-1) = ft(n,ict-1) + frac*ftold * &
                  ((ph(n,ict)-ph(n,ict+1))/(ph(n,ict-1)-ph(n,ict))) * &
                  cpn(ict)/cpn(ict-1)
      fuold = fu(n,ict)
      fu(n,ict) = fu(n,ict)*(d_one-frac)
      fu(n,ict-1) = fu(n,ict-1) + frac*fuold * &
                  ((ph(n,ict)-ph(n,ict+1))/(ph(n,ict-1)-ph(n,ict)))
      fvold = fv(n,ict)
      fv(n,ict) = fv(n,ict)*(d_one-frac)
      fv(n,ict-1) = fv(n,ict-1) + frac*fvold * &
                  ((ph(n,ict)-ph(n,ict+1))/(ph(n,ict-1)-ph(n,ict)))
      if ( chemcutran ) then
        do k = 1 , ntra
          ftraold = ftra(n,ict,k)
          ftra(n,ict,k) = ftra(n,ict,k)*(d_one-frac)
          ftra(n,ict-1,k) = ftra(n,ict-1,k) + frac*ftraold * &
                          (ph(n,ict)-ph(n,ict+1))/(ph(n,ict-1)-ph(n,ict))
        end do
      end if
      !
      ! Very slightly adjust tendencies to force exact
      ! enthalpy, momentum and tracer conservation
      !
      ents = d_zero
      uav = d_zero
      vav = d_zero
      do i = 1 , ict
        ents = ents + (cpn(i)*ft(n,i)+lv(i)*fq(n,i))*(ph(n,i)-ph(n,i+1))
        uav = uav + fu(n,i)*(ph(n,i)-ph(n,i+1))
        vav = vav + fv(n,i)*(ph(n,i)-ph(n,i+1))
      end do
      ents = ents/(ph(n,1)-ph(n,ict+1))
      uav = uav/(ph(n,1)-ph(n,ict+1))
      vav = vav/(ph(n,1)-ph(n,ict+1))
      do i = 1 , ict
        ft(n,i) = ft(n,i) - ents/cpn(i)
        fu(n,i) = (d_one-cu)*(fu(n,i)-uav)
        fv(n,i) = (d_one-cu)*(fv(n,i)-vav)
      end do
      if ( chemcutran ) then
        do k = 1 , ntra
          traav = d_zero
          do i = 1 , ict
            traav = traav + ftra(n,i,k)*(ph(n,i)-ph(n,i+1))
          end do
          traav = traav/(ph(n,1)-ph(n,ict+1))
          do i = 1 , ict
            ftra(n,i,k) = ftra(n,i,k) - traav
          end do
        end do
      end if
      !
      ! Xu, K.-M., and S. K. Krueger:
      !   Evaluation of cloudiness parameterizations using a cumulus
      !   ensemble model, Mon. Wea. Rev., 119, 342-367, 1991.
      !
      ! Identify Deep concection if cloud depth is > 2000m
      !
      if ( zcup(n,ict) - zcup(n,icb) >= 2000.0_rkx ) then
        do i = icb , ict
          cldfra(n,i) = kfac_deep*log(d_one+(k2_const*d_half*(m(i)+m(i+1))))
          cldfra(n,i) = min(max(0.01_rkx,cldfra(n,i)),0.6_rkx)
        end do
      else
        do i = icb , ict
          cldfra(n,i) = kfac_shal*log(d_one+(k2_const*d_half*(m(i)+m(i+1))))
          cldfra(n,i) = min(max(0.01_rkx,cldfra(n,i)),0.2_rkx)
        end do
      end if
      kcb(n) = icb
      kct(n) = ict
    end do pointloop

    contains

#include <wlh.inc>
#include <pfesat.inc>
#include <pfqsat.inc>

      !
      ! Calculate lifting level temperature
      !
      subroutine tlift(n,p,t,q,qs,gz,icb,nk,tvp,tpk,clw,nd,nl,kk)
        implicit none
        integer , intent(in) :: n , nd , nk , nl , kk
        integer(ik4) , intent(in) :: icb
        real(rkx) , dimension(:,:) , intent(in) :: p , q , qs , t
        real(rkx) , dimension(nd) , intent(in) :: gz
        real(rkx) , dimension(nd) , intent(inout) :: clw , tpk , tvp
        real(rkx) :: ah0 , ahg , alv , cpinv , cpp , qg , rg , s , tg
        real(rkx) :: ppa
        !real(rkx) :: tc , es , denom
        integer(ik4) :: i , j , nsb , nst
        !
        ! calculate certain parcel quantities, including static energy
        !
        ah0 = (cpd*(d_one-q(n,nk))+cpv*q(n,nk))*t(n,nk) + q(n,nk) * &
              wlh(t(n,nk)) + gz(nk)
        cpp = cpd*(d_one-q(n,nk)) + cpv*q(n,nk)
        cpinv = d_one/cpp
        if ( kk == 1 ) then
          !
          ! calculate lifted parcel quantities below cloud base
          !
          do i = 1 , icb - 1
            clw(i) = d_zero
          end do
          do i = nk , icb - 1
            tpk(i) = t(n,nk) - (gz(i)-gz(nk))*cpinv
            tvp(i) = tpk(i)*(d_one+q(n,nk)*rgowi)
          end do
        end if
        !
        ! find lifted parcel quantities above cloud base
        !
        nst = icb
        nsb = icb
        if ( kk == 2 ) then
          nst = nl
          nsb = icb + 1
        end if
        do i = nsb , nst
          tg = t(n,i)
          qg = qs(n,i)
          alv = wlh(tg)
          do j = 1 , 2
            s = d_one/(cpd + alv*alv*qg/(rwat*t(n,i)*t(n,i)))
            ahg = cpd*tg + (cpv-cpd)*q(n,nk)*t(n,i) + alv*qg + gz(i)
            tg = tg + s*(ah0-ahg)
            !tg = max(tg + s*(ah0-ahg),35.0_rkx)
            !tc = tg - tzero
            !denom = 243.5_rkx + tc
            !if ( tc >= 0.0_rkx ) then
            !  es = 6.112_rkx * exp(17.67_rkx*tc/denom)
            !else
            !  es = exp(23.33086_rkx-6111.72784_rkx/tg+0.15215_rkx*log(tg))
            !end if
            !qg = rgow * es/(p(n,i) - es * (1.0_rkx-rgow))
            ppa = p(n,i)*100.0_rkx
            qg = pfqsat(tg,ppa)
          end do
          tpk(i) = (ah0-(cpv-cpd)*q(n,nk)*t(n,i)-gz(i)-alv*qg)*rcpd
          clw(i) = q(n,nk) - qg
          clw(i) = max(d_zero,clw(i))
          rg = qg/(d_one-q(n,nk))
          tvp(i) = tpk(i)*(d_one+rg*rgowi)
        end do
      end subroutine tlift

  end subroutine cupeman

end module mod_cu_em
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

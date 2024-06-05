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
!
!    Documentation:
!
!    Authors: code loosely based on toy model developed by Tompkins
!        implemented and further developed by Nogherotto in REGCM in 2013-2015
!
!    Code implements a 5 phase prognostic cloud microphysics scheme with a
!    diagnostic treatment of cloud fraction
!
!    - The cloud cover is parameterized using the diagnostic scheme of
!      Xu and Randall 96
!
!    - Microphysics were originally based on Tiedtke 93 and Tompkins et al. 2007
!      but with a number of modifications to the melting, collection,
!      autoconversion
!
!    - Solver is a simple implicit solver, implemented layer by layer
!      Original tests in REGCM4 published in Nogherotto et al. 2016
!
!    Modifications:
!      201704: Tompkins and Nogherotto
!       a) repetition of code removed for layer k=1 and tidy up of
!          duplicate/redundant arrays
!       b) Bug fixes to autoconversion terms for all parameterization scheme
!          (CCOVER missing)
!       c) Sedimentation term moved before BF-type process and loss sink
!          included in the first guess to improve supercooled liquid water
!          at top of mixed phase clouds (no need for cloud top distance
!          deposition "fudge".
!       d) Factor 2 fudge removed from deposition term
!       e) A number of clean up terms added to remove tiny cloud/precipiation
!          amounts
!       f) Collection fudge factor removed from autoconversion and replaced by
!          explicit parameterization from Khairoutdinov and Kogan [2000].
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_micro_nogtom
  use mod_realkinds
  use mod_dynparam
  use mod_stdio
  use mod_mpmessage
  use mod_memutil
  use mod_service
  use mod_regcm_types
  use mod_constants , only : d_zero , d_one , d_half , d_two , d_1000
  use mod_constants , only : dlowval , mathpi
  use mod_constants , only : tzero , rtice , rtwat_rtice_r
  use mod_constants , only : c5alvcp , c5alscp , rhoh2o , rovcp
  use mod_constants , only : wlhfocp , wlhsocp , wlhvocp
  use mod_constants , only : rwat , wlhs , wlhv
  use mod_constants , only : c5les , c5ies , c3ies , c3les , c4les , c4ies
  use mod_constants , only : c2es , ep1
  use mod_constants , only : egrav , regrav , ep1
  use mod_runparams , only : nqx
  use mod_runparams , only : iqqv => iqv !vapor
  use mod_runparams , only : iqql => iqc !liquid
  use mod_runparams , only : iqqr => iqr !rain
  use mod_runparams , only : iqqi => iqi !ice
  use mod_runparams , only : iqqs => iqs !snow
  use mod_runparams , only : dt , rdt
  use mod_runparams , only : ipptls , ichem , iaerosol , iindirect , rcrit
  use mod_runparams , only : budget_compute , nssopt , iautoconv
  use mod_runparams , only : auto_rate_khair , auto_rate_kessl , &
                             auto_rate_klepi , rcldiff
  use mod_runparams , only : rkconv , skconv , rcovpmin , rpecons

#ifdef DEBUG
  use mod_runparams , only : stats
#endif

  implicit none

  private

  logical , parameter :: lmicro = .true.

  ! critical autoconversion
  real(rkx) , parameter :: rlcritsnow = 4.e-5_rkx

  real(rkx) , parameter :: auto_expon_khair = 1.47_rkx
  real(rkx) , parameter :: rldcp = d_one/wlhfocp  ! Cp/Lf
  ! 1/autoconversion time scale (s)
  real(rkx) , parameter :: autocrit_kessl = 5.e-4_rkx
  real(rkx) , parameter :: rclcrit_land = 5.e-4_rkx
  real(rkx) , parameter :: rclcrit_sea = 3.e-4_rkx
  real(rkx) , parameter :: rhcrit_lnd = 0.80_rkx
  real(rkx) , parameter :: rhcrit_sea = 0.90_rkx
  real(rkx) , parameter :: rprc1 = 3.e2_rkx  ! in Sundqvist = 300
  real(rkx) , parameter :: siglow = 0.8_rkx
  ! Cloud fraction threshold that defines cloud top
  real(rkx) , parameter :: cldtopcf = 0.1_rkx
  ! Fraction of deposition rate in cloud top layer
  real(rkx) , parameter :: depliqrefrate = 0.1_rkx
  ! Depth of supercooled liquid water layer (m)
  real(rkx) , parameter :: depliqrefdepth = 500.0_rkx
  ! max threshold rh for evaporation for a precip coverage of zero
  real(rkx) , parameter :: rprecrhmax = 0.7_rkx
  ! evaporation rate coefficient Numerical fit to wet bulb temperature
  !real(rkx) , parameter :: tw1 = 1329.31_rkx
  !real(rkx) , parameter :: tw2 = 0.0074615_rkx
  !real(rkx) , parameter :: tw3 = 0.85e5_rkx
  !real(rkx) , parameter :: tw4 = 40.637_rkx
  !real(rkx) , parameter :: tw5 = 275.0_rkx
  !real(rkx) , parameter :: rtaumel = 1.1880e4_rkx
  ! temperature homogeneous freezing
  real(rkx) , parameter :: thomo = 235.16_rkx  ! -38.00 Celsius
  ! initial mass of ice particle
  real(rkx) , parameter :: iceinit = 1.e-12_rkx
  real(rkx) , parameter :: rkoop1 = 2.583_rkx
  real(rkx) , parameter :: rkoop2 = 0.48116e-2_rkx ! 1/207.8
  !------------------------------------------------
  real(rkx) , parameter :: ciden13 = 8.87_rkx      ! ice density 700**0.333
  real(rkx) , parameter :: airconduct = 2.4e-2_rkx ! conductivity of air
  real(rkx) , parameter :: spherefac = (4.0_rkx/3.0_rkx)*mathpi

  public :: allocate_mod_nogtom , init_nogtom , nogtom

  ! Total water and enthalpy budget diagnostics variables
  ! marker for water phase of each species
  ! 0 = vapour, 1 = liquid, 2 = ice
  integer(ik4) , pointer , dimension(:) :: iphase
  ! marks melting linkage for ice categories
  ! ice->liquid, snow->rain
  integer(ik4) , pointer , dimension(:) :: imelt
  logical , pointer , dimension(:) :: lfall

  real(rkx) , pointer , dimension(:,:,:):: sumh0 , sumq0
  real(rkx) , pointer , dimension(:,:,:) :: sumh1 , sumq1
  real(rkx) , pointer , dimension(:,:) :: errorq , errorh
  real(rkx) , pointer , dimension(:,:,:):: tentkp
  real(rkx) , pointer , dimension(:,:,:,:) :: tenqkp
  ! distance from the top of the cloud
  real(rkx) , pointer , dimension(:,:,:) :: cldtopdist
  ! Mass variables
  ! Microphysics
  real(rkx) , pointer , dimension(:,:,:) :: dqsatdt
  ! for convection detrainment source and subsidence source/sink terms
  real(rkx) , pointer , dimension(:) :: convsrce
  real(rkx) , pointer , dimension(:,:,:) :: eewmt
  ! fluxes convergence of species
  real(rkx) , pointer , dimension(:,:,:) :: qliq

  real(rkx) , pointer , dimension(:,:,:) :: eew
  ! ice water saturation
  real(rkx) , pointer , dimension(:,:,:) :: qsice
  ! diagnostic mixed phase RH
  real(rkx) , pointer , dimension(:,:,:) :: qsmix
  ! Storage for eeliq , eeice
  real(rkx) , pointer , dimension(:,:,:) :: eeliq
  real(rkx) , pointer , dimension(:,:,:) :: eeice
  ! water/ice saturation mixing ratio
  real(rkx) , pointer , dimension(:,:,:) :: eeliqt
  ! liq+rain sedim flux
  real(rkx) , pointer , dimension(:,:,:) :: pfplsl
  ! ice+snow sedim flux
  real(rkx) , pointer , dimension(:,:,:) :: pfplsn
  ! decoupled temperature tendency
  real(rkx) , pointer , dimension(:,:,:) :: ttendc
  ! critical factors
  real(rkx) , pointer , dimension(:,:) :: xlcrit
  real(rkx) , pointer , dimension(:,:) :: rhcrit
  ! fall speeds of three categories
  real(rkx) , pointer , dimension(:) :: vqx
  ! decoupled mixing ratios tendency
  real(rkx) , pointer , dimension(:,:,:,:) :: qxtendc
  ! j,i,n ! generalized precipitation flux
  real(rkx) , pointer , dimension(:,:,:,:) :: pfplsx
  real(rkx) , pointer, dimension(:,:,:,:) :: qx
  real(rkx) , pointer, dimension(:,:,:) :: tx
  real(rkx) , pointer, dimension(:,:,:) :: pf
  real(rkx) , pointer, dimension(:,:,:) :: qdetr
  real(rkx) , pointer, dimension(:,:,:) :: delz
  real(rkx) , pointer, dimension(:,:,:) :: ph
  real(rkx) , pointer, dimension(:,:,:) :: rho
  real(rkx) , pointer, dimension(:,:,:) :: ccn
  real(rkx) , pointer, dimension(:,:,:) :: pverv
  real(rkx) , pointer, dimension(:,:,:) :: heatrt
  real(rkx) , pointer, dimension(:,:,:) :: fcc
  real(rkx) , pointer, dimension(:,:,:) :: remrat
  ! saturation mixing ratio with respect to water
  real(rkx) , pointer, dimension(:,:,:) :: qsliq
  ! koop
  ! se T < 0 la nuvola si forma o quando q e' maggiore della liquid
  ! water saturation minima, oppure se e' maggiore del mixing ratio
  ! wrt ice critica a cui inizia l'homogeneaous ice nucleation
  ! At temperatures below 0 degC new cloud forms in any non-cloudy part
  ! of the grid box where the humidity exceeds either the minimum of
  ! the liquid water saturation specific humidity (qsl), or the
  ! critical vapour saturation mixing ratio with respect to ice at
  ! which homogeneous ice nucleation initiates
  ! empirical fit given by Karcher and Lohmann (2002) which is a
  ! function of temperature and ranges from 45% supersaturation at
  ! T = 235 K to 67% at T = 190 K.
  ! At temperatures warmer than -38 degC the cloud formation over a
  ! timestep results entirely in liquid cloud,
  ! i.e. koop = eeliq/eeice, mentre per T < -38 koop = RHhomo
  ! while below this threshold the liquid water or aqueous sulphate
  ! solutes are assumed to freeze instantaneously and the process is
  ! a source for cloud ice.
  ! koop modifies the ice saturation mixing ratio for homogeneous
  ! nucleation
  real(rkx) , pointer, dimension(:,:,:) :: koop
  ! Delta pressure
  real(rkx) , pointer, dimension(:,:,:) :: dpfs

  real(rkx) , parameter :: zerocf = 0.01_rkx
  real(rkx) , parameter :: onecf  = 0.99_rkx

  real(rkx) , parameter :: activqx = 1.0e-12_rkx
  real(rkx) , parameter :: verylowqx = 1.0e-20_rkx
  real(rkx) , parameter :: activcf = 2.0_rkx*zerocf
  real(rkx) , parameter :: maxsat  = 0.5_rkx

  contains

  subroutine allocate_mod_nogtom
    implicit none
    call getmem1d(vqx,1,nqx,'cmicro:vqx')
    call getmem1d(imelt,1,nqx,'cmicro:imelt')
    call getmem1d(lfall,1,nqx,'cmicro:lfall')
    call getmem1d(iphase,1,nqx,'cmicro:iphase')
    call getmem3d(qliq,1,kz,jci1,jci2,ici1,ici2,'cmicro:qliq')
    call getmem3d(eewmt,1,kz,jci1,jci2,ici1,ici2,'cmicro:eewmt')
    call getmem3d(qsmix,1,kz,jci1,jci2,ici1,ici2,'cmicro:qsmix')
    call getmem3d(ttendc,1,kz,jci1,jci2,ici1,ici2,'cmicro:ttendc')
    call getmem1d(convsrce,1,nqx,'cmicro:convsrce')
    call getmem3d(eew,1,kz,jci1,jci2,ici1,ici2,'cmicro:eew')
    call getmem3d(qsice,1,kz,jci1,jci2,ici1,ici2,'cmicro:qsice')
    call getmem4d(qx,1,nqx,1,kz,jci1,jci2,ici1,ici2,'cmicro:qx')
    call getmem3d(tx,1,kz,jci1,jci2,ici1,ici2,'cmicro:tx')
    call getmem3d(qdetr,1,kz,jci1,jci2,ici1,ici2,'cmicro:qdetr')
    call getmem3d(delz,1,kz,jci1,jci2,ici1,ici2,'cmicro:delz')
    call getmem3d(ph,1,kz,jci1,jci2,ici1,ici2,'cmicro:ph')
    call getmem3d(rho,1,kz,jci1,jci2,ici1,ici2,'cmicro:rho')
    call getmem3d(ccn,1,kz,jci1,jci2,ici1,ici2,'cmicro:ccn')
    call getmem3d(pverv,1,kz,jci1,jci2,ici1,ici2,'cmicro:pverv')
    call getmem3d(heatrt,1,kz,jci1,jci2,ici1,ici2,'cmicro:heatrt')
    call getmem3d(fcc,1,kz,jci1,jci2,ici1,ici2,'cmicro:fcc')
    if ( ichem == 1 ) then
      call getmem3d(remrat,1,kz,jci1,jci2,ici1,ici2,'cmicro:remrat')
    end if
    call getmem3d(pf,1,kzp1,jci1,jci2,ici1,ici2,'cmicro:pf')
    call getmem3d(qsliq,1,kz,jci1,jci2,ici1,ici2,'cmicro:qsliq')
    call getmem3d(cldtopdist,1,kz,jci1,jci2,ici1,ici2,'cmicro:cldtopdist')
    call getmem3d(dqsatdt,jci1,jci2,ici1,ici2,1,kz,'cmicro:dqsatdt')
    call getmem3d(pfplsl,1,kzp1,jci1,jci2,ici1,ici2,'cmicro:pfplsl')
    call getmem3d(pfplsn,1,kzp1,jci1,jci2,ici1,ici2,'cmicro:pfplsn')
    call getmem3d(koop,1,kz,jci1,jci2,ici1,ici2,'cmicro:koop')
    call getmem2d(xlcrit,jci1,jci2,ici1,ici2,'cmicro:xlcrit')
    call getmem2d(rhcrit,jci1,jci2,ici1,ici2,'cmicro:rhcrit')
    call getmem3d(eeliq,1,kz,jci1,jci2,ici1,ici2,'cmicro:eeliq')
    call getmem3d(eeice,1,kz,jci1,jci2,ici1,ici2,'cmicro:eeice')
    call getmem3d(eeliqt,1,kz,jci1,jci2,ici1,ici2,'cmicro:eeliqt')
    call getmem4d(qxtendc,1,nqx,1,kz,jci1,jci2,ici1,ici2,'cmicro:qxtendc')
    call getmem4d(pfplsx,1,nqx,1,kzp1,jci1,jci2,ici1,ici2,'cmicro:pfplsx')
    call getmem3d(dpfs,1,kz,jci1,jci2,ici1,ici2,'cmicro:dpfs')
    if ( budget_compute ) then
      call getmem3d(sumq0,1,kz,jci1,jci2,ici1,ici2,'cmicro:sumq0')
      call getmem3d(sumh0,1,kz,jci1,jci2,ici1,ici2,'cmicro:sumh0')
      call getmem3d(sumq1,1,kz,jci1,jci2,ici1,ici2,'cmicro:sumq1')
      call getmem3d(sumh1,1,kz,jci1,jci2,ici1,ici2,'cmicro:sumh1')
      call getmem2d(errorq,jci1,jci2,ici1,ici2,'cmicro:errorq')
      call getmem2d(errorh,jci1,jci2,ici1,ici2,'cmicro:errorh')
      call getmem4d(tenqkp,1,nqx,1,kz,jci1,jci2,ici1,ici2,'cmicro:tenqkp')
      call getmem3d(tentkp,1,kz,jci1,jci2,ici1,ici2,'cmicro:tentkp')
    end if
  end subroutine allocate_mod_nogtom

  subroutine init_nogtom(ldmsk)
    use mod_runparams , only : vfqr , vfqi , vfqs
    implicit none
    integer(ik4) , pointer , dimension(:,:) , intent(in) :: ldmsk
    integer(ik4) :: i , j , n
    ! Define species phase, 0 = vapour, 1 = liquid, 2 = ice
    iphase(iqqv) = 0
    iphase(iqql) = 1
    iphase(iqqi) = 2
    iphase(iqqr) = 1
    iphase(iqqs) = 2

    ! Set up melting/freezing index,
    ! if an ice category melts/freezes, where does it go?

    imelt(iqqv) = -99
    imelt(iqql) = iqqi
    imelt(iqqi) = iqql
    imelt(iqqr) = iqqs
    imelt(iqqs) = iqqr

    ! Set the fall velocities
    vqx(iqqv) = d_zero ! * sqrt(QX(JL,JK,IQV))
    vqx(iqql) = d_zero ! * sqrt(QX(JL,JK,IQL))
    vqx(iqqi) = vfqi   !0.15_rkx * sqrt(QX(JL,JK,IQI))
    vqx(iqqr) = vfqr   !4.0_rkx  * sqrt(QX(JL,JK,IQR))
    vqx(iqqs) = vfqs   !1.0_rkx  * sqrt(QX(JL,JK,IQS))

    ! Set lfall
    do n = 1 , nqx
      if ( vqx(n) > d_zero ) then
        lfall(n) = .true. !falling species
      else
        lfall(n) = .false.
      end if
    end do

    ! modify autoconversion threshold dependent on:
    ! land (polluted, high ccn, smaller droplets, higher threshold)
    ! sea  (clean, low ccn, larger droplets, lower threshold)
    do concurrent ( j = jci1:jci2, i = ici1:ici2 )
      if ( ldmsk(j,i) == 1 ) then ! landmask =1 land
        xlcrit(j,i) = rclcrit_land ! landrclcrit_land = 5.e-4
        rhcrit(j,i) = rhcrit_lnd
      else
        xlcrit(j,i) = rclcrit_sea  ! oceanrclcrit_sea  = 3.e-4
        rhcrit(j,i) = rhcrit_sea
      end if
    end do
  end subroutine init_nogtom

#ifdef DEBUG
  subroutine nogtom(mo2mc,mc2mo,ngs)
    implicit none
    type(nogtom_stats) , intent(inout) :: ngs
#else
  subroutine nogtom(mo2mc,mc2mo)
    implicit none
#endif
    type(mod_2_micro) , intent(in) :: mo2mc
    type(micro_2_mod) , intent(out) :: mc2mo

    integer(ik4) :: i , j , k , kk , n , m , jn , jo
    real(rkx) :: tnew , dp , qe , tmpl , tmpi
    real(rkx) :: qprev , hprev , zdelta , prainx , psnowx
    ! for sedimentation source/sink terms
    real(rkx) , dimension(nqx) :: fallsrce , fallsink
    ! n x n matrix storing the LHS of implicit solver
    real(rkx) , dimension(nqx,nqx) :: qlhs
    ! explicit sources and sinks "q s exp"=q source explicit
    real(rkx) , dimension(nqx,nqx) :: qsexp
    ! implicit sources and sinks "q s imp"=q source/sink implicit
    real(rkx) , dimension(nqx,nqx) :: qsimp
    ! Initial values
    real(rkx) , dimension(nqx) :: qx0 , qxfg , qxn
    real(rkx) , dimension(nqx) :: ratio , sinksum
    ! array for sorting explicit terms
    ! integer(ik4) , dimension(nqx) :: iorder
    real(rkx) :: tk , tc , dens , pbot , snowp , rainp , supsat , subsat
    real(rkx) :: totcond , qliqfrac , qicefrac , qicetot
    real(rkx) :: gdp , dtgdp , rdtgdp , alpha1 , ldefr
    real(rkx) :: facl , faci , facw , corr , acond , zdl , infactor
    real(rkx) :: alfaw , phases , qexc , rhc , zsig , preclr , arg
    real(rkx) :: rexplicit , xlcondlim , tmpa , zrh , beta , beta1
    real(rkx) :: cond , dtdp , cdmax , tdiff , cons1 , qvnew
    ! local real constants for evaporation
    real(rkx) :: dpr , denom , dpevap , evapi , evapl , excess
    real(rkx) :: dqsmixdt , dqsicedt , dqsliqdt
    real(rkx) :: corqsliq , corqsice , corqsmix , evaplimmix
    real(rkx) :: ql_incld , qi_incld , qli_incld , ldifdt , sink
    ! Cloud coverage and clearsky portion
    real(rkx) :: covptot , covpclr
    ! real(rkx) :: botm , rm
    real(rkx) :: qold , tcond , dqs , chng , chngmax , icenuclei
    real(rkx) :: qpretot , fluxq
    ! constants for deposition process
    real(rkx) :: vpice , vpliq , xadd , xbdd , cvds , &
                 qice0 , qinew , rainaut , snowaut
    ! constants for condensation and turbulent mixing erosion of clouds
    real(rkx) :: dpmxdt , wtot , dtdiab , dtforc , &
                 qp , qsat , cond1 , levap , leros
    real(rkx) :: qsmixv , ccover , lccover , rain , rainh
    real(rkx) :: precip , cfpr , acrit
    logical :: lccn , lerror , ldetr , lconden , lactiv , locast
    logical :: ltkgt0 , ltklt0 , ltkgthomo , lcloud
    logical , dimension(nqx,nqx) :: lind2
    integer(ik4) :: ii , jj , ll , imax , nn
    real(rkx) :: aamax , dum , xsum , swap
    real(rkx) , dimension(nqx) :: vv
    integer(ik4) , dimension(nqx) :: indx

#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'microphys'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    lccn = ( ichem == 1 .and. iaerosol == 1 .and. iindirect == 2 )

    if ( idynamic == 3 ) then
      do concurrent ( n = 1:nqx, k = 1:kz, j = jci1:jci2, i = ici1:ici2 )
        qxtendc(n,k,j,i) = mc2mo%qxten(j,i,k,n)
      end do
      do concurrent ( k = 1:kz, j = jci1:jci2, i = ici1:ici2 )
        ttendc(k,j,i) = mc2mo%tten(j,i,k)
      end do
    else
      ! Decouple tendencies
      do concurrent ( n = 1:nqx, k = 1:kz, j = jci1:jci2, i = ici1:ici2 )
        qxtendc(n,k,j,i) = mc2mo%qxten(j,i,k,n) / mo2mc%psb(j,i)
      end do
      do concurrent ( k = 1:kz, j = jci1:jci2, i = ici1:ici2 )
        ttendc(k,j,i) = mc2mo%tten(j,i,k) / mo2mc%psb(j,i)
      end do
    end if

    ! Define the initial array qx and reset total precipitation variables
    do concurrent ( n = 1:nqx, k = 1:kz, j = jci1:jci2, i = ici1:ici2 )
      qx(n,k,j,i) = mo2mc%qxx(j,i,k,n)
      pfplsx(n,k,j,i) = d_zero
    end do

    ! Copy with transpose all the other variables
    do concurrent ( k = 1:kz, j = jci1:jci2, i = ici1:ici2 )
      tx(k,j,i) = mo2mc%t(j,i,k)
      ph(k,j,i) = mo2mc%phs(j,i,k)
      qdetr(k,j,i) = mo2mc%qdetr(j,i,k)
      delz(k,j,i) = mo2mc%delz(j,i,k)
      rho(k,j,i) = mo2mc%rho(j,i,k)
      pverv(k,j,i) = mo2mc%pverv(j,i,k)
      heatrt(k,j,i) = mo2mc%heatrt(j,i,k)
      fcc(k,j,i) = mo2mc%cldf(j,i,k)
    end do

    if ( lccn ) then
      do concurrent ( k = 1:kz, j = jci1:jci2, i = ici1:ici2 )
        ccn(k,j,i) = mo2mc%ccn(j,i,k)
      end do
    end if

    do concurrent ( k = 1:kzp1, j = jci1:jci2, i = ici1:ici2 )
      pf(k,j,i) = mo2mc%pfs(j,i,k)
    end do

    ! Delta pressure
    do concurrent ( k = 1:kz, j = jci1:jci2, i = ici1:ici2 )
      dpfs(k,j,i) = pf(k+1,j,i)-pf(k,j,i)
    end do

    !-----------------------------------
    ! initialization for cloud variables
    ! -------------------------------------
    ! Define qliq the function for mixed phase
    !     PHASE is calculated to distinguish the three cases:
    !     PHASE = 1            water phase
    !     PHASE = 0            ice phase
    !     0 < PHASE < 1        mixed phase
    ! Define pressure at full levels
    ! pf = Pressure on fuLL levels (Pa)
    ! Define a new array for detrainment

    do concurrent ( k = 1:kz, j = jci1:jci2, i = ici1:ici2 )
      qliq(k,j,i) = max(min(d_one,((max(rtice,min(tzero, &
                        tx(k,j,i)))-rtice)*rtwat_rtice_r)**2),d_zero)
    end do

    ! Compute supersaturations
    do concurrent ( k = 1:kz, j = jci1:jci2, i = ici1:ici2 )
      eeliq(k,j,i) = c2es*exp(c3les*((tx(k,j,i)-tzero)/(tx(k,j,i)-c4les)))
      eeice(k,j,i) = c2es*exp(c3ies*((tx(k,j,i)-tzero)/(tx(k,j,i)-c4ies)))
      koop(k,j,i) = min(rkoop1-rkoop2*tx(k,j,i),eeliq(k,j,i)/eeice(k,j,i))
    end do

    !-------------------------------------
    ! Initial enthalpy and total water diagnostics
    !-------------------------------------
    !
    ! Starting budget if requested
    !
    if ( budget_compute ) then
      ! Record the tendencies
      do concurrent ( n = 1:nqx, k = 1:kz, j = jci1:jci2, i = ici1:ici2 )
        tenqkp(n,k,j,i) = qxtendc(n,k,j,i)
      end do
      do concurrent ( k = 1:kz, j = jci1:jci2, i = ici1:ici2 )
        tentkp(k,j,i) = ttendc(k,j,i)
        sumq0(k,j,i)  = d_zero
        sumh0(k,j,i)  = d_zero
      end do
      ! initialize the flux arrays
#ifndef __GFORTRAN__
      do concurrent ( j = jci1:jci2, i = ici1:ici2 ) &
        local(tnew , dp , qe , tmpl , tmpi , alfaw)
#else
      do i = ici1 , ici2
        do j = jci1 , jci2
#endif
          tnew = tx(1,j,i)
          dp = dpfs(1,j,i)
          qe = qdetr(1,j,i)
          tmpl = qx(iqql,1,j,i)+qx(iqqr,1,j,i)
          tmpi = qx(iqqi,1,j,i)+qx(iqqs,1,j,i)
          tnew = tnew - wlhvocp*tmpl - wlhsocp*tmpi
          sumq0(1,j,i) = sumq0(1,j,i)+(tmpl+tmpi+qx(iqqv,1,j,i))*dp*regrav
          ! Detrained water treated here
          if ( lmicro .and. abs(qe) > activqx ) then
            sumq0(1,j,i) = sumq0(1,j,i) + qe*dp*regrav
            alfaw = qliq(1,j,i)
            tnew = tnew-(wlhvocp*alfaw+wlhsocp*(d_one-alfaw))*qe
          end if
          sumh0(1,j,i) = sumh0(1,j,i) + dp*tnew
#ifdef __GFORTRAN__
        end do
#endif
      end do
#ifndef __GFORTRAN__
      do concurrent ( j = jci1:jci2, i = ici1:ici2 ) &
        local(tnew,tmpi,alfaw,tmpl,qe,dp,qprev,hprev,k)
#else
      do i = ici1 , ici2
        do j = jci1 , jci2
#endif
          qprev = sumq0(1,j,i)
          hprev = sumh0(1,j,i)
          do k = 2 , kz
            tnew = tx(k,j,i)
            dp = dpfs(k,j,i)
            qe = qdetr(k,j,i)
            sumq0(k,j,i) = qprev ! total water
            sumh0(k,j,i) = hprev ! liquid water temperature
            tmpl = qx(iqql,k,j,i)+qx(iqqr,k,j,i)
            tmpi = qx(iqqi,k,j,i)+qx(iqqs,k,j,i)
            tnew = tnew - wlhvocp*tmpl - wlhsocp*tmpi
            sumq0(k,j,i) = sumq0(k,j,i)+(tmpl+tmpi+qx(iqqv,k,j,i))*dp*regrav
            ! Detrained water treated here
            if ( lmicro .and. abs(qe) > activqx ) then
              sumq0(k,j,i) = sumq0(k,j,i) + qe*dp*regrav
              alfaw = qliq(k,j,i)
              tnew = tnew-(wlhvocp*alfaw+wlhsocp*(d_one-alfaw))*qe
            end if
            sumh0(k,j,i) = sumh0(k,j,i) + dp*tnew
            qprev = sumq0(k,j,i)
            hprev = sumh0(k,j,i)
          end do
#ifdef __GFORTRAN__
        end do
#endif
      end do
      do concurrent ( k = 1:kz, j = jci1:jci2, i = ici1:ici2 )
        sumh0(k,j,i) = sumh0(k,j,i)/pf(k+1,j,i)
      end do
    end if ! budget_compute

    ! -------------------------------
    ! Define saturation values
    !---------------------------
#ifndef __GFORTRAN__
    do concurrent ( k = 1:kz, j = jci1:jci2, i = ici1:ici2 ) &
      local(zdelta , phases)
#else
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
#endif
          ! zdelta = 1 if t > tzero
          ! zdelta = 0 if t < tzero
          zdelta = max(d_zero,sign(d_one,tx(k,j,i)-tzero))
          !---------------------------------------------
          ! mixed phase saturation
          !--------------------------------------------
          phases = qliq(k,j,i)
          eewmt(k,j,i) = eeliq(k,j,i)*phases + eeice(k,j,i)*(d_one-phases)
          eewmt(k,j,i) = min(eewmt(k,j,i)/ph(k,j,i),maxsat)
          qsmix(k,j,i) = eewmt(k,j,i)
          ! ep1 = rwat/rgas - d_one
          qsmix(k,j,i) = qsmix(k,j,i)/(d_one-ep1*qsmix(k,j,i))
          !--------------------------------------------
          ! ice saturation T < 273K
          ! liquid water saturation for T > 273K
          !--------------------------------------------
          eew(k,j,i) = (zdelta*eeliq(k,j,i) + &
               (d_one-zdelta)*eeice(k,j,i))/ph(k,j,i)
          eew(k,j,i) = min(eew(k,j,i),maxsat)
          !ice water saturation
          qsice(k,j,i) = min(eeice(k,j,i)/ph(k,j,i),maxsat)
          qsice(k,j,i) = qsice(k,j,i)/(d_one-ep1*qsice(k,j,i))
          !----------------------------------
          ! liquid water saturation
          !----------------------------------
          !eeliq is the saturation vapor pressure es(T)
          !the saturation mixing ratio is ws = es(T)/p *0.622
          !ws = ws/(-(d_one/eps - d_one)*ws)
          eeliqt(k,j,i) = min(eeliq(k,j,i)/ph(k,j,i),maxsat)
          qsliq(k,j,i) = eeliqt(k,j,i)
          qsliq(k,j,i) = qsliq(k,j,i)/(d_one-ep1*qsliq(k,j,i))
#ifdef __GFORTRAN__
        end do
      end do
#endif
    end do

    !--------------------------------ADEED BY RITA
    ! Calculate distance from cloud top
    ! defined by cloudy layer below a layer with cloud frac <0.01
    !--------------------------------------------------------------

    do concurrent ( k = 1:kz, j = jci1:jci2, i = ici1:ici2 )
      cldtopdist(k,j,i) = d_zero
    end do
#ifndef __GFORTRAN__
    do concurrent ( j = jci1:jci2, i = ici1:ici2 ) local(k,kk)
#else
    do i = ici1 , ici2
      do j = jci1 , jci2
#endif
        do k = 2 , kz
          do kk = 2 , k
            if ( fcc(kk-1,j,i) > cldtopcf .and. &
                 fcc(kk,j,i)  <= cldtopcf ) then
              cldtopdist(k,j,i) = cldtopdist(k,j,i) + delz(kk,j,i)
            end if
          end do
        end do
#ifdef __GFORTRAN__
      end do
#endif
    end do

#ifdef DEBUG
    if ( stats ) then
      ngs%statssupw(:,:,:) = d_zero
      ngs%statssupc(:,:,:) = d_zero
      ngs%statserosw(:,:,:) = d_zero
      ngs%statserosc(:,:,:) = d_zero
      ngs%statsdetrw(:,:,:) = d_zero
      ngs%statsdetrc(:,:,:) = d_zero
      ngs%statsevapw(:,:,:) = d_zero
      ngs%statsevapc(:,:,:) = d_zero
      ngs%statscond1w(:,:,:) = d_zero
      ngs%statscond1c(:,:,:) = d_zero
      ngs%statsdepos(:,:,:) = d_zero
      ngs%statsmelt(:,:,:) = d_zero
      ngs%statsfrz(:,:,:) = d_zero
      ngs%statsrainev(:,:,:) = d_zero
      ngs%statssnowev(:,:,:) = d_zero
      ngs%statsautocvw(:,:,:) = d_zero
      ngs%statsautocvc(:,:,:) = d_zero
    end if
#endif
    !
    !----------------------------------------------------------------------
    !                       START OF VERTICAL LOOP
    !----------------------------------------------------------------------
    !
    ! Loop over points
    !
#ifndef __GFORTRAN__
    do concurrent ( j = jci1:jci2, i = ici1:ici2 ) &
      local(fallsrce,fallsink,ldefr,qlhs,qsexp,qsimp,qx0,qxfg,qxn,lind2, &
            ratio,sinksum,tk,tc,dens,pbot,snowp,rainp,supsat,subsat,     &
            totcond,qliqfrac,qicefrac,qicetot,dp,gdp,dtgdp,rdtgdp,alpha1,&
            facl,faci,facw,corr,acond,zdl,infactor,alfaw,phases,qexc,    &
            rhc,zsig,qe,preclr,arg,rexplicit,xlcondlim,tmpa,zrh,beta,    &
            beta1,cond,dtdp,cdmax,tdiff,cons1,dpr,denom,dpevap,evapi,    &
            evapl,excess,dqsmixdt,dqsicedt,dqsliqdt,corqsliq,corqsice,   &
            corqsmix,evaplimmix,ql_incld,qi_incld,qli_incld,ldifdt,sink, &
            covptot,covpclr,qold,tcond,dqs,chng,chngmax,icenuclei,       &
            acrit,precip,cfpr,qpretot,fluxq,vpice,vpliq,xadd,xbdd,cvds,  &
            qice0,qinew,rainaut,snowaut,dpmxdt,wtot,dtdiab,dtforc,qp,    &
            qsat,cond1,levap,leros,qsmixv,ccover,lccover,k,n,m,jn,jo,    &
            ldetr,lconden,lactiv,locast,ltkgt0,ltklt0,ltkgthomo,lcloud,  &
            ii,jj,kk,ll,imax,nn,aamax,dum,xsum,swap,vv,indx)
#else
    do i = ici1 , ici2
      do j = jci1 , jci2
#endif

        pbot = pf(kzp1,j,i)
        covptot = d_zero
        covpclr = d_zero

        ! Loop over levels

        do k = 1 , kz

          supsat      = d_zero
          subsat      = d_zero
          ldefr       = d_zero
          do n = 1 , nqx
            fallsrce(n) = d_zero
            fallsink(n) = d_zero
            convsrce(n) = d_zero
          end do

          !-------------------------------------------------------
          ! SOURCE/SINK array for implicit and explicit terms
          !-------------------------------------------------------
          !
          ! a POSITIVE value entered into the arrays is a...
          !
          !             Source of this variable
          !             |
          !             |   Sink of this variable
          !             |   |
          !             V   V
          ! qsexp/imp:q(iqA,iqB)
          !
          ! Thus if qsexp/imp(iql,iqv) = qk where qk > 0 then this is
          ! a source of iql and a sink of iqv
          !
          ! put 'magic' source terms such as qdetr from
          ! detrainment into explicit source/sink array diagnognal
          ! qsexp(iql,iql) = qdetr
          !--------------------------------------------------------
          ! Define the microphysics
          ! the matrix will be sparse is this a problem ?
          ! (X,Y) means a sink of X and a source of Y
          ! for the implementation I will use flexible pointers
          ! such that it will be written (iqr,iqg) to indicate graupel to rain
          ! and the parametrization can have different variables switched on
          ! and off.
          ! each of these is a parametrization for a microphysical process.
          !--------------------------------------------------------
          !
          do n = 1 , nqx
            do m = 1 , nqx
              qsexp(m,n)  = d_zero
              qsimp(m,n)  = d_zero
            end do
          end do
          !
          !---------------------------------
          ! First guess microphysics
          !---------------------------------
          do n = 1 , nqx
            qx0(n)  = qx(n,k,j,i)
            qxfg(n) = qx0(n)
          end do

          ldetr = ( abs(qdetr(k,j,i)) > activqx )
          totcond = qxfg(iqql)+qxfg(iqqi)
          lconden = ( qxfg(iqql) > activqx .and. qxfg(iqqi) > activqx )
          if ( lconden ) then
            qliqfrac = qxfg(iqql)/totcond
            qicefrac = d_one-qliqfrac
          else
            if ( qxfg(iqql) > activqx ) then
              qliqfrac = d_one
              qicefrac = d_zero
            else if ( qxfg(iqqi) > activqx ) then
              qliqfrac = d_zero
              qicefrac = d_one
            else
              qliqfrac = d_zero
              qicefrac = d_zero
            end if
          end if

          qicetot = d_zero
          do n = 1 , nqx
            if ( iphase(n) == 2 ) then
              qicetot = qicetot + qxfg(n)
            end if
          end do

          dp       = dpfs(k,j,i)
          tk       = tx(k,j,i)
          tc       = tk - tzero
          dens     = rho(k,j,i)
          qsmixv   = qsmix(k,j,i)
          ccover   = fcc(k,j,i)
          ccover   = min(max(ccover,zerocf),onecf)

          if ( k == 1 ) then
            lccover = d_zero
            rainp   = d_zero
            snowp   = d_zero
          else
            lccover = fcc(k-1,j,i)
            lccover = min(max(lccover,zerocf),onecf)
            rainp   = pfplsx(iqqr,k,j,i)
            snowp   = pfplsx(iqqs,k,j,i)
          end if

          ltkgt0    = ( tk > tzero )
          ltklt0    = ( .not. ltkgt0 )
          ltkgthomo = ( tk > thomo )
          lcloud    = ( ccover > activcf )
          locast    = ( ccover >= onecf )

          ! Derived variables needed
          gdp = egrav/dp       ! g/dp  =(1/m)
          dtgdp = dt*gdp       ! (dt*g)/dp =(dt/m)
          rdtgdp = d_one/dtgdp ! dp/(gdt)=m/dt  [Kg/m2/s]
          !------------------------------------
          ! calculate dqs/dT
          !------------------------------------
          ! liquid
          facw     = c5les/((tk - c4les)**2)
          corr     = d_one/(d_one - ep1*eeliqt(k,j,i))
          dqsliqdt = facw*corr*qsliq(k,j,i)
          corqsliq = d_one + wlhvocp*dqsliqdt
          ! ice
          faci     = c5ies/((tk - c4ies)**2)
          corr     = d_one/(d_one - ep1*eew(k,j,i))
          dqsicedt = faci*corr*qsice(k,j,i)
          corqsice = d_one + wlhsocp*dqsicedt
          ! diagnostic mixed
          alfaw    = qliq(k,j,i)
          facl     = alfaw*facw + (d_one - alfaw)*faci
          corr     = d_one/(d_one - ep1*eewmt(k,j,i))
          dqsmixdt = facl*corr*qsmixv
          corqsmix = d_one/(d_one + eldcpm(tk)*dqsmixdt)
          !--------------------------------
          ! evaporation/sublimation limits
          !--------------------------------
          evaplimmix = max((qsmixv-qxfg(iqqv))*corqsmix,d_zero)

          !--------------------------------
          ! in-cloud consensate amount
          !--------------------------------
          tmpa = d_one/ccover
          ql_incld = qxfg(iqql)*tmpa
          qi_incld = qxfg(iqqi)*tmpa
          qli_incld  = ql_incld + qi_incld

          !------------------------------------------------------------------
          !  MICROPHYSICS START HERE
          !------------------------------------------------------------------

          !------------------------------------------------------------------
          ! Turn on/off microphysics
          !------------------------------------------------------------------

          if ( lmicro ) then

            !-------------------------------------------------------
            !  FALL SOURCE
            !-------------------------------------------------------
            qpretot = d_zero
            if ( k > 1 ) then
              do n = 1 , nqx
                if ( lfall(n) ) then
                  ! Source from layer above
                  fallsrce(n) = pfplsx(n,k,j,i)*dtgdp
                  qsexp(n,n) = qsexp(n,n) + fallsrce(n)
                  qxfg(n) = qxfg(n) + fallsrce(n)
                  qpretot = qpretot + qxfg(n)
                endif
              end do
            else
              do n = 1 , nqx
                if ( lfall(n) ) then
                  qpretot = qpretot + qxfg(n)
                end if
              end do
            end if

            !------------------------------------------------
            ! Evaporate very small amounts of liquid and ice
            !------------------------------------------------

            if ( qx0(iqql) < verylowqx ) then
              qsexp(iqqv,iqql) = qsexp(iqqv,iqql) + qx0(iqql)
              qsexp(iqql,iqqv) = qsexp(iqql,iqqv) - qx0(iqql)
              qxfg(iqql) = qxfg(iqql) - qx0(iqql)
              qxfg(iqqv) = qxfg(iqqv) + qx0(iqql)
            end if
            if ( qx0(iqqi) < verylowqx ) then
              qsexp(iqqv,iqqi) = qsexp(iqqv,iqqi) + qx0(iqqi)
              qsexp(iqqi,iqqv) = qsexp(iqqi,iqqv) - qx0(iqqi)
              qxfg(iqqi) = qxfg(iqqi) - qx0(iqqi)
              qxfg(iqqv) = qxfg(iqqv) + qx0(iqqi)
            end if

            !------------------------------------------------------------------
            !  SEDIMENTATION/FALLING OF *ALL* MICROPHYSICAL SPECIES
            !
            !     now that rain and snow species are prognostic
            !     the precipitation flux can be defined directly level
            !     by level
            !     There is no vertical memory required from the flux
            !     variable
            !
            !     *AMT* moved sedimentation before the deposition and
            !     included sink in first guess in order to account for
            !     supercooled water enhancement at cloud top
            !
            !------------------------------------------------------------------
            do n = 1 , nqx
              if ( lfall(n) ) then
                ! Sink to next layer, constant fall speed
                ! *AMT* now included in first guess.
                sink = vqx(n) * dens * dtgdp
                fallsink(n) = sink
                qxfg(n) = qxfg(n)/(d_one+sink)
              end if  !lfall
            end do ! n

            !-----------------------------------------------------------------
            !  ICE SUPERSATURATION ADJUSTMENT
            !-----------------------------------------------------------------
            ! Note that the supersaturation adjustment is made with respect to
            ! liquid saturation:  when T > 0C
            ! ice saturation:     when T < 0C
            !                     with an adjustment made to allow for ice
            !                     supersaturation in the clear sky
            ! Note also that the KOOP factor automatically clips the
            ! supersaturation to a maximum set by the liquid water saturation
            ! mixing ratio
            ! important for temperatures near to but below 0C
            ! qv_max = qs * (fcc + (1-fcc) *RH_homo ) if T < 0C
            ! qv_max = qs                             if T > 0C
            !-----------------------------------------------------------------
            !-----------------------------------
            ! Supersaturation limit (from Koop)
            !-----------------------------------
            if ( nssopt == 0 )  then
              facl = d_one
            else
              if ( ltkgt0 ) then
                facl = d_one
              else
                facl = ccover + koop(k,j,i)*(d_one-ccover)
              end if
            end if

            !-----------------------------------------------------------------
            ! Calculate supersaturation wrt Koop including dqs/dT
            ! correction factor
            !-----------------------------------------------------------------
            ! Here the supersaturation is turned into liquid water
            ! However, if the temperature is below the threshold for homogeneous
            ! freezing then the supersaturation is turned instantly to ice.
            ! Moreover the RH is clipped to the limit of
            ! qv_max = qs * (fcc + (1-fcc) *RH_homo )
            !------------------------------------------------------------------
            supsat = max((qxfg(iqqv)-facl*qsmixv)*corqsmix,d_zero)
            ! e < esi, because for e > esi ice still present
            subsat = min((qxfg(iqqv)-facl*qsmixv)*corqsmix,d_zero)
            if ( supsat > dlowval ) then
              if ( ltkgthomo ) then
                ! turn supersaturation into liquid water
                qsexp(iqql,iqqv) = qsexp(iqql,iqqv) + supsat
                qsexp(iqqv,iqql) = qsexp(iqqv,iqql) - supsat
                qxfg(iqql) = qxfg(iqql) + supsat
                qxfg(iqqv) = qxfg(iqqv) - supsat
#ifdef DEBUG
                if ( stats ) then
                  ngs%statssupw(k,j,i) = ngs%statssupw(k,j,i) + supsat
                end if
#endif
              else if  ( ltklt0 ) then
                ! turn supersaturation into ice water
                qsexp(iqqi,iqqv) = qsexp(iqqi,iqqv) + supsat
                qsexp(iqqv,iqqi) = qsexp(iqqv,iqqi) - supsat
                qxfg(iqqi) = qxfg(iqqi) + supsat
                qxfg(iqqv) = qxfg(iqqv) - supsat
#ifdef DEBUG
                if ( stats ) then
                  ngs%statssupc(k,j,i) = ngs%statssupc(k,j,i) - supsat
                end if
#endif
              end if
            else
              if ( subsat < d_zero .and. lconden .and. .not. lcloud ) then
                ! turn subsaturation into vapor, where there is no cloud
                excess = totcond + subsat
                if ( excess < d_zero ) then
                  if ( ltkgthomo ) then
                    evapl = max(-qxfg(iqql),subsat)
                    qsexp(iqqv,iqql) = qsexp(iqqv,iqql) - evapl
                    qsexp(iqql,iqqv) = qsexp(iqql,iqqv) + evapl
                    qxfg(iqql) = qxfg(iqql) + evapl
                    qxfg(iqqv) = qxfg(iqqv) - evapl
                  else if  ( ltklt0 ) then
                    evapi = max(-qxfg(iqqi),subsat)
                    ! turn subsaturation into vapour
                    qsexp(iqqv,iqqi) = qsexp(iqqv,iqqi) - evapi
                    qsexp(iqqi,iqqv) = qsexp(iqqi,iqqv) + evapi
                    qxfg(iqqi) = qxfg(iqqi) + evapi
                    qxfg(iqqv) = qxfg(iqqv) - evapi
                  end if
                end if
              end if
            end if
            !
            !call addpath(iqql,iqqv,supsatl,qsexp,qsimp,d_zero,qxfg)
            !call addpath(iqqi,iqqv,supsati,qsexp,qsimp,d_zero,qxfg)
            !
            !-------------------------------------------------------
            ! source/sink array for implicit and explicit terms
            !-------------------------------------------------------
            !
            ! a positive value is:
            !
            !        Source   Sink of this variable
            !             |   |
            !             V   V
            ! QSEXP/IMP:q(IQa,IQb)
            !
            ! Thus if QSEXP/IMP(IQL,IQV) = K where K > 0 then this is
            ! a source of IQL and a sink of IQV
            !
            ! put external source terms in the diagonal entries
            !--------------------------------------------------------

            !------------------------------------------------------------------
            ! convective detrainment
            !------------------------------------------------------------------
            if ( ldetr ) then
              !qice = 1 if T < 250, qice = 0 if T > 273
              qe = qdetr(k,j,i)
              alfaw = qliq(k,j,i)
              convsrce(iqql) = alfaw*qe
              convsrce(iqqi) = (d_one-alfaw)*qe
              qsexp(iqql,iqql) = qsexp(iqql,iqql) + convsrce(iqql)
              qsexp(iqqi,iqqi) = qsexp(iqqi,iqqi) + convsrce(iqqi)
              qxfg(iqql) = qxfg(iqql) + convsrce(iqql)
              qxfg(iqqi) = qxfg(iqqi) + convsrce(iqqi)
#ifdef DEBUG
              if ( stats ) then
                ngs%statsdetrw(k,j,i) = convsrce(iqql)
                ngs%statsdetrc(k,j,i) = convsrce(iqqi)
              end if
#endif
            end if

            !---------------------------------------
            ! EROSION OF CLOUDS BY TURBULENT MIXING
            !--------------------------------------
            ! rcldiff  : Diffusion coefficient for evaporation by turbulent
            ! mixing (IBID., EQU. 30) rcldiff = 1.0e-6_rkx
            ldifdt = rcldiff*dt
            !Increase by factor of 5 for convective points
            if ( lconden ) then
              leros = ccover * ldifdt * max(qsmixv-qxfg(iqqv),d_zero)
              leros = min(leros,evaplimmix)
              leros = min(leros,totcond)
              facl = qliqfrac*leros
              faci = qicefrac*leros
              qsexp(iqql,iqqv) = qsexp(iqql,iqqv) - facl
              qsexp(iqqv,iqql) = qsexp(iqqv,iqql) + facl
              qsexp(iqqi,iqqv) = qsexp(iqqi,iqqv) - faci
              qsexp(iqqv,iqqi) = qsexp(iqqv,iqqi) + faci
              qxfg(iqql) = qxfg(iqql) - facl
              qxfg(iqqi) = qxfg(iqqi) - faci
#ifdef DEBUG
              if ( stats ) then
                ngs%statserosw(k,j,i) = qliqfrac*leros
                ngs%statserosc(k,j,i) = qicefrac*leros
              end if
#endif
            end if

            !------------------------------------------------------------------
            ! condensation/evaporation due to dqsat/dt
            !------------------------------------------------------------------
            ! calculate dqs/dt and use to calculate the cloud source
            ! note that old diagnostic mix phased qsat is retained for moment
            !------------------------------------------------------------------
            dtdp   = rovcp*tk/ph(k,j,i)
            dpmxdt = dp*rdt
            wtot   = pverv(k,j,i)
            wtot   = min(dpmxdt,max(-dpmxdt,wtot))
            dtdiab = min(dpmxdt*dtdp, &
                     max(-dpmxdt*dtdp,heatrt(k,j,i)))*dt+wlhfocp*ldefr
            ! ldefr = 0
            ! note: ldefr should be set to the difference between the mixed
            ! phase functions in the convection and cloud scheme, and
            ! for now we set it to zero and the functions are the same.
            ! In RegCM not all convection schemes provide such info.
            dtforc = dtdp*wtot*dt + dtdiab
            qold   = qsmixv
            tcond  = tk + dtforc
            tcond  = max(tcond,160.0_rkx)
            ! the goal is to produce dqs = qsmix - qold, where qsmix is
            ! reduced because of the condensation. so that dqs is negative?
            qp = d_one/ph(k,j,i)
            phases = max(min(d_one,((max(rtice,min(tzero, &
                       tcond))-rtice)*rtwat_rtice_r)**2),d_zero)
            ! saturation mixing ratio ws
            qsat = eewm(tcond,phases) * qp
            qsat = min(qsat,maxsat)          ! ws < 0.5        WHY?
            corr  = d_one/(d_one-ep1*qsat)
            qsat = qsat*corr
            cond = (qsmixv-qsat)/(d_one + qsat*edem(tcond,phases))
            tcond = tcond + eldcpm(tcond)*cond
            phases = max(min(d_one,((max(rtice,min(tzero, &
                       tcond))-rtice)*rtwat_rtice_r)**2),d_zero)
            qsmixv = qsmixv - cond
            qsat = eewm(tcond,phases) * qp
            qsat = min(qsat,maxsat)
            corr = d_one/(d_one-ep1*qsat)
            qsat = qsat*corr
            cond1 = (qsmixv-qsat)/(d_one + qsat*edem(tcond,phases))
            tcond = tcond + eldcpm(tcond)*cond1
            qsmixv = qsmixv - cond1
            dqs = qsmixv - qold
            qsmixv = qold

            !----------------------------------------------------------------
            ! dqs > 0:  evaporation of clouds
            !----------------------------------------------------------------
            ! erosion term is explicit in for cloud liquid
            ! changed to be uniform distribution in cloud region
            ! previous function based on delta distribution in cloud:
            if ( dqs > d_zero ) then
              !levap = C*min( dqs/dt , (qi+ql)/C )
              levap = ccover*min(dqs,qli_incld)
              levap = min(levap,evaplimmix)
              levap = min(levap,max(qsmixv-qxfg(iqqv),d_zero))
              facl = qliqfrac*levap
              faci = qicefrac*levap
              qsexp(iqqv,iqql) = qsexp(iqqv,iqql) + facl
              qsexp(iqql,iqqv) = qsexp(iqql,iqqv) - facl
              qsexp(iqqv,iqqi) = qsexp(iqqv,iqqi) + faci
              qsexp(iqqi,iqqv) = qsexp(iqqi,iqqv) - faci
              qxfg(iqql) = qxfg(iqql) - facl
              qxfg(iqqi) = qxfg(iqqi) - faci
#ifdef DEBUG
              if ( stats ) then
                ngs%statsevapw(k,j,i) = qliqfrac*levap
                ngs%statsevapc(k,j,i) = qicefrac*levap
              end if
#endif
            !-----------------------------------------------------------------
            ! dqs < 0: formation of clouds
            !-----------------------------------------------------------------
            else if ( dqs < d_zero ) then
              ! (1) increase of cloud water in existing clouds
              if ( lcloud ) then
                ! new limiter
                chng = -dqs
                ! old limiter
                !  (significantly improves upper tropospheric humidity rms)
                if ( locast ) then
                  corr = d_one/(d_one-ep1*qsmixv)
                  cdmax = (qxfg(iqqv)-qsmixv)/(d_one+corr*qsmixv*edem(tk,alfaw))
                else
                  cdmax = (qxfg(iqqv)-ccover*qsmixv)/ccover
                end if
                chng = min(chng,cdmax)
                chng = ccover*chng
                chng = max(chng,d_zero)
                !-------------------------------------------------------------
                ! All increase goes into liquid unless so cold cloud
                ! homogeneously freezes
                ! include new liquid formation in first guess value, otherwise
                ! liquid remains at cold temperatures until next timestep.
                !-------------------------------------------------------------
                if ( ltkgthomo ) then
                  qsexp(iqql,iqqv) = qsexp(iqql,iqqv) + chng
                  qsexp(iqqv,iqql) = qsexp(iqqv,iqql) - chng
                  qxfg(iqql) = qxfg(iqql) + chng
                  qxfg(iqqv) = qxfg(iqqv) - chng
#ifdef DEBUG
                  if ( stats ) then
                    ngs%statscond1w(k,j,i) = chng
                  end if
#endif
                else if ( ltklt0 ) then
                  qsexp(iqqi,iqqv) = qsexp(iqqi,iqqv) + chng
                  qsexp(iqqv,iqqi) = qsexp(iqqv,iqqi) - chng
                  qxfg(iqqi) = qxfg(iqqi) + chng
                  qxfg(iqqv) = qxfg(iqqv) - chng
#ifdef DEBUG
                  if ( stats ) then
                    ngs%statscond1c(k,j,i) = chng
                  end if
#endif
                end if
              else
                ! (2) generation of new clouds (dc/dt>0)
                qexc = 0.0_rkx ! Make compiler happy
                select case (nssopt)
                  case (0,1)
                    qexc = max((qxfg(iqqv)-ccover*qsmixv) / &
                        (d_one-ccover),d_zero)
                  case (2) ! Khairoutdinov and Kogan (2000)
                    qexc = qxfg(iqqv)
                  case (3) ! Kessler(1969)
                    qexc = qxfg(iqqv)/totcond
                end select
                rhc = rhcrit(j,i)
                zsig = ph(k,j,i)/pbot
                if ( zsig > siglow ) then
                  ! increase RHcrit to 1.0 towards the surface (sigma>0.8)
                  rhc = rhc + (d_one-rhc)*((zsig-siglow)/(d_one-siglow))**2
                end if
                ! supersaturation options
                if ( ltkgt0 .or. nssopt == 0 ) then
                  ! no ice supersaturation allowed
                  facl = d_one
                else
                  ! ice supersaturation
                  facl = koop(k,j,i)
                end if
                if ( qexc >= rhc*qsmixv*facl .and. qexc < qsmixv*facl ) then
                  ! note: not **2 on 1-a term if qe is used.
                  ! added correction term fac to numerator 15/03/2010
                  acond = -(d_one-ccover)*facl*dqs / &
                          max(d_two*(facl*qsmixv-qexc),dlowval)
                  acond = min(acond,d_one-ccover) ! put the limiter back
                  ! linear term:
                  ! added correction term fac 15/03/2010
                  chng = -facl*dqs*d_half*acond !mine linear
                  ! new limiter formulation
                  ! qsice(k,j,i)-qexc) /
                  tmpa = d_one-ccover
                  zdl = d_two*(facl*qsmixv-qexc) / tmpa
                  ! added correction term fac 15/03/2010
                  if ( facl*dqs < -zdl ) then
                    ! qsice(k,j,i)+qvnow
                    xlcondlim = (ccover-d_one)*facl*dqs-facl*qsmixv+qxfg(iqqv)
                    chng = min(chng,xlcondlim)
                  end if
                  chng = max(chng,d_zero)
                  !-------------------------------------------------------------
                  ! all increase goes into liquid unless so cold cloud
                  ! homogeneously freezes
                  ! include new liquid formation in first guess value, otherwise
                  ! liquid remains at cold temperatures until next timestep.
                  !-------------------------------------------------------------
                  if ( ltkgthomo ) then
                    chng = min(chng,qxfg(iqql))
                    qsexp(iqql,iqqv) = qsexp(iqql,iqqv) + chng
                    qsexp(iqqv,iqql) = qsexp(iqqv,iqql) - chng
                    qxfg(iqql) = qxfg(iqql) + chng
                    qxfg(iqqv) = qxfg(iqqv) - chng
                  else
                    ! homogeneous freezing
                    chng = min(chng,qxfg(iqqi))
                    qsexp(iqqi,iqqv) = qsexp(iqqi,iqqv) + chng
                    qsexp(iqqv,iqqi) = qsexp(iqqv,iqqi) - chng
                    qxfg(iqqi) = qxfg(iqqi) + chng
                    qxfg(iqqv) = qxfg(iqqv) - chng
                  end if
#ifdef DEBUG
                  if ( stats ) then
                    ngs%statscond1c(k,j,i) = ngs%statscond1c(k,j,i) + chng
                  end if
#endif
                end if
              end if
            end if

            !------------------------------------------------------------------
            ! DEPOSITION:
            ! Growth of ice by vapour deposition
            ! and fudged ice contact nucleation included here.
            !
            !------------------------------------------------------------------
            ! Following Rotstayn et al. 2001 and Meyers et al. 1992
            !
            ! clouds are exactly saturated with
            ! respect to liquid water (well mixed), (or koop)
            !
            ! Growth considered as sink of liquid water
            !
            ! Bergeron-Findeisen adjustment not required.
            !
            ! Can not treat if liquid not present as would require
            ! additional variable to model in-cloud vapour mixing ratio
            !
            ! *AMT* 03/2017 removed factor 2, and cloud top reduction
            ! introduce enhancement due to contact nucleation when
            ! collisions occurs between liquid and ice crystals
            ! By considering sedimentation first and including the
            ! implicit loss term in the first guess of ice.
            !--------------------------------------------------------------
            lactiv = qx0(iqql) > activqx .and. ltklt0
            if ( lactiv ) then
              vpice = eeice(k,j,i) !saturation vapor pressure wrt ice
              vpliq = eeliq(k,j,i) !saturation vapor pressure wrt liq
              ! Meyers et al 1992
              icenuclei = d_1000*exp(12.96_rkx * &
                          ((vpliq-vpice)/vpice)-0.639_rkx)

              !---------------------------------------------------------
              ! *AMT* contact nucleation fudge factor
              ! Note this refers to contact between liquid and ice
              ! crystals
              ! not contact nucleation by contact with heterogeneous
              ! nuclei
              ! process acts as 1/liqfrac , when liqfrac=1, no speed up
              ! this is the max(activqx,qliqfrac) factor...
              !---------------------------------------------------------

              xadd  = wlhs*(wlhs/(rwat*tk)-d_one)/(airconduct*tk)
              xbdd  = rwat*tk*ph(k,j,i)/(2.21_rkx*vpice)
              cvds = 7.8_rkx * (icenuclei/dens)**0.666_rkx * &
                     (vpliq-vpice)/(ciden13*(xadd+xbdd)*vpice)
              cvds = max(cvds,d_zero)

              !---------------------------------------------------
              ! iceinit = 1.e-12 is initial mass of ice particle
              !           used if no ice present to start process
              !---------------------------------------------------
              qice0 = max(qi_incld, icenuclei*iceinit/dens)

              !-----------------------------------------------------
              ! new value of ice mixing ratio
              ! Note: eqn 8 in Rotstayn et al. (2000) is incorrect
              !-----------------------------------------------------
              qinew = (0.666_rkx*cvds*dt+qice0**0.666_rkx)**1.5_rkx
              qinew = max(qinew,d_zero)

              !-------------------------------------------------------
              ! grid-mean deposition rate:
              ! Use of CCOVER assumes that clouds are completely well
              ! mixed
              !-------------------------------------------------------
              chng = ccover*(qinew-qice0)
              !re-added by Rita 3/2/2022
              infactor = min(icenuclei/15000.0_rkx,d_one)
              chng = chng*min(infactor + (d_one-infactor)* &
                  (depliqrefrate+cldtopdist(k,j,i)/depliqrefdepth),d_one)
              chng = min(chng,qxfg(iqql))

              !-------------------------------------------------------------
              ! limit deposition to liquid water amount
              ! can't treat vapour in ice-only cloud without extra
              ! prognostic variable
              !-------------------------------------------------------------
              chng = max(chng,d_zero)

              !--------------
              ! add to matrix
              !--------------
              qsexp(iqqi,iqql) = qsexp(iqqi,iqql) + chng
              qsexp(iqql,iqqi) = qsexp(iqql,iqqi) - chng
              qxfg(iqql) = qxfg(iqql) - chng
              qxfg(iqqi) = qxfg(iqqi) + chng
#ifdef DEBUG
              if ( stats ) then
                ngs%statsdepos(k,j,i) = chng
              end if
#endif
            end if

            tmpa = d_one/ccover
            ql_incld = qxfg(iqql)*tmpa
            qi_incld = qxfg(iqqi)*tmpa

            !---------------------------------------------------------------
            ! Precip cover overlap using MAX-RAN Overlap
            ! Since precipitation is now prognostic we must
            !   1) apply an arbitrary minimum coverage (0.3) if precip>0
            !   2) abandon the 2-flux clr/cld treatment
            !   3) Thus, since we have no memory of the clear sky precip
            !      fraction, we mimic the previous method by reducing
            !      COVPTOT(JL), which has the memory, proportionally with
            !      the precip evaporation rate, taking cloud fraction
            !      into account
            !   #3 above leads to much smoother vertical profiles of
            !   precipitation fraction than the Klein-Jakob scheme which
            !   monotonically increases precip fraction and then resets
            !   it to zero in a step function once clear-sky precip reaches
            !   zero.
            !   Maximum overlap for clouds in adjacent levels and random
            !   overlap for clouds separated by clear levels.
            !---------------------------------------------------------------
            if ( qpretot > d_zero ) then
              covptot = d_one - ((d_one-covptot) * &
                  (d_one - max(ccover,lccover))/(d_one-lccover))
              covptot = max(covptot,rcovpmin)
              covpclr = max(covptot-ccover,d_zero)
            else
              covptot = d_zero ! no flux - reset cover
              covpclr = d_zero ! no flux - reset cover
            end if
            ! clear sky proportion

            !---------------------------------------------------------------
            !   WARM PHASE AUTOCONVERSION
            !---------------------------------------------------------------
            rainaut = 0.0_rkx ! Make compiler happy
            if ( ql_incld > d_zero ) then
              select case (iautoconv)
                case (1) ! Klein & Pincus (2000)
                  rainaut = dt*auto_rate_klepi*(ql_incld**(2.3_rkx))
                  qsimp(iqql,iqqv) = d_zero
                  qsimp(iqqr,iqql) = qsimp(iqqr,iqql) + rainaut
                  qsexp(iqqr,iqql) = d_zero
                case (2) ! Khairoutdinov and Kogan (2000)
                  rainaut = dt*auto_rate_khair*(ql_incld**(auto_expon_khair))
                  qsimp(iqql,iqqv) = d_zero
                  qsimp(iqqr,iqql) = qsimp(iqqr,iqql) + rainaut
                case (3) ! Kessler(1969)
                  rainaut = dt*auto_rate_kessl*autocrit_kessl
                  qsimp(iqql,iqqv) = d_zero
                  qsexp(iqqr,iqql) = qsexp(iqqr,iqql) - rainaut
                  qsexp(iqql,iqqr) = qsexp(iqql,iqqr) + rainaut
                  qsimp(iqqr,iqql) = qsimp(iqqr,iqql) + rainaut
                case (4) ! Sundqvist
                  !alpha1 = min(rkconv*dt,ql_incld)
                  alpha1 = rkconv*dt
                  acrit = xlcrit(j,i)
                  if ( lccn ) then
                    if ( ccn(k,j,i) > 0._rkx ) then
                      ! aerosol second indirect effect on autoconversion
                      ! threshold, rcrit is a critical cloud radius for cloud
                      ! water undergoing autoconversion
                      ! ccn = number of ccn /m3
                      acrit = ccn(k,j,i)*spherefac * &
                          ((rcrit*1e-6_rkx)**3)*rhoh2o
                    end if
                  end if
                  !-------------------------------------------------------
                  ! parameters for cloud collection by rain and snow.
                  ! note that with new prognostic variable it is now
                  ! possible to replace this with an explicit collection
                  ! parametrization to be replaced by
                  ! Khairoutdinov and Kogan [2000]:
                  !-------------------------------------------------------
                  if ( covptot > d_zero ) then
                    precip = (rainp+snowp)/covptot
                    cfpr = d_one + rprc1*sqrt(max(precip,d_zero))
                    alpha1 = alpha1*cfpr
                    acrit = acrit/cfpr
                  end if

                  ! security for exp for some compilers
                  arg = (ql_incld/acrit)**2
                  if ( arg < 25.0_rkx ) then
                    rainaut = alpha1*(d_one - exp(-arg))
                  else
                    rainaut = alpha1
                  end if
                  ! clean up
                  qsimp(iqql,iqqv) = d_zero
                  if ( ltkgt0 ) then
                    qsimp(iqqr,iqql) = qsimp(iqqr,iqql) + rainaut
                  else
                    !-----------------------
                    ! rain freezes instantly
                    !-----------------------
                    qsimp(iqqs,iqql) = qsimp(iqqs,iqql) + rainaut
                  end if
              end select
#ifdef DEBUG
              if ( stats ) then
                if ( ltkgt0 ) then
                  ngs%statsautocvw(k,j,i) = ngs%statsautocvw(k,j,i) + rainaut
                else
                  ngs%statsautocvc(k,j,i) = ngs%statsautocvc(k,j,i) + rainaut
                end if
              end if
#endif
! save the precip production for chem. wet. dep.
              if ( ichem == 1 )  then
                remrat(k,j,i) = rainaut/dt
              end if
            end if ! appreciable liquid cloud

            !------------
            ! Cold clouds
            !------------
            if ( ltklt0 ) then
              ! Snow Autoconversion rate follow Lin et al. 1983
              if ( qi_incld > d_zero ) then
                alpha1 = dt*skconv*exp(0.025_rkx*tc)
                arg = (qi_incld/rlcritsnow)**2
                if ( arg < 25.0_rkx ) then
                  snowaut = alpha1 * (d_one - exp(-arg))
                else
                  snowaut = alpha1
                end if
                qsimp(iqqs,iqqi) = qsimp(iqqs,iqqi) + snowaut
#ifdef DEBUG
                if ( stats ) then
                  ngs%statsautocvc(k,j,i) = ngs%statsautocvc(k,j,i) + snowaut
                end if
#endif
              end if
            else
              !---------------------------------------------------------------
              !                         MELTING
              !---------------------------------------------------------------
              ! The melting of ice and snow are treated explicitly.
              ! First water and ice saturation are found
              !---------------------------------------------
              ! ice saturation T < 273K
              ! liquid water saturation for T > 273K
              !---------------------------------------------
              qicetot = qxfg(iqqi)+qxfg(iqqs)
              if ( qicetot > d_zero ) then
                ! Calculate subsaturation
                ! qsice(k,j,i)-qxfg(iqqv),d_zero)
                subsat = max(qsmixv-qxfg(iqqv),d_zero)
                ! Calculate difference between dry-bulb (t)  and the temperature
                ! at which the wet-bulb = 0degC
                ! Melting only occurs if the wet-bulb temperature >0
                ! i.e. warming of ice particle due to melting > cooling
                ! due to evaporation.
                ! The wet-bulb temperature is used in order to account for the
                ! thermal (cooling) ect of evaporation on the melting process
                ! in sub-saturated air. The evaporation counteracts the latent
                ! heating due to melting and allows snow particles to survive
                ! to slightly warmer temperatures when the relative
                ! humidity of the air is low. The wet-bulb temperature is
                ! approximated as in the scheme described by
                ! Wilson and Ballard(1999): Tw = Td-(qs-q)(A+B(p-c)-D(Td-E))
                ! tdiff = tc - subsat * &
                !     (tw1+tw2*(ph(k,j,i)-tw3)-tw4*(tk-tw5))
                tdiff = tc
                ! Ensure CONS1 is positive so that MELTMAX = 0 if TDMTW0 < 0
                ! cons1 = abs(dt*(d_one + d_half*tdiff)/rtaumel)
                ! cons1 = dt/rtaumel
                cons1 = d_one
                chngmax = max(tdiff*cons1*rldcp,d_zero)
                if ( chngmax > d_zero ) then
                  ! Loop over frozen hydrometeors (iphase == 2 (ice, snow))
                  do n = 1, nqx
                    if ( iphase(n) == 2 ) then
                      m = imelt(n) ! imelt(iqqi)=iqql, imelt(iqqs)=iqqr
                      if ( m < 0 ) cycle
                      phases = qxfg(n)/qicetot
                      chng = min(qxfg(n),phases*chngmax)
                      chng = max(chng,d_zero)
                      ! n = iqqi,iqqs; m = iqql,iqqr
                      qsexp(m,n) =  qsexp(m,n) + chng
                      qsexp(n,m) =  qsexp(n,m) - chng
                      qxfg(n) =  qxfg(n) - chng
                      qxfg(m) =  qxfg(m) + chng
#ifdef DEBUG
                      if ( stats ) then
                        ngs%statsmelt(k,j,i) = ngs%statsmelt(k,j,i) + chng
                      end if
#endif
                    end if
                  end do
                end if
              end if
            end if

            !------------------------------------------------------------!
            !                         FREEZING                           !
            !------------------------------------------------------------!

            ! Freezing of rain.
            ! All rain freezes in a timestep if the temperature is below 0 C
            ! calculate sublimation latent heat

            chngmax = max((tzero-tk)*rldcp,d_zero)
            if ( chngmax > d_zero .and. qx0(iqqr) > activqx ) then
              chng = min(qxfg(iqqr),chngmax)
              chng = max(chng,d_zero)
              qsexp(iqqs,iqqr) = qsexp(iqqs,iqqr) + chng
              qsexp(iqqr,iqqs) = qsexp(iqqr,iqqs) - chng
              qxfg(iqqs) = qxfg(iqqs) + chng
              qxfg(iqqr) = qxfg(iqqr) - chng
#ifdef DEBUG
              if ( stats ) then
                ngs%statsfrz(k,j,i) = chng
              end if
#endif
            end if

            !-------------------
            ! Freezing of liquid
            !-------------------

            chngmax = max((thomo-tk)*rldcp,d_zero)
            if ( chngmax > d_zero .and. qx0(iqql) > activqx ) then
              chng = min(qxfg(iqql),chngmax)
              chng = max(chng,d_zero)
              qsexp(iqqi,iqql) = qsexp(iqqi,iqql) + chng
              qsexp(iqql,iqqi) = qsexp(iqql,iqqi) - chng
              qxfg(iqql) = qxfg(iqql) - chng
              qxfg(iqqi) = qxfg(iqqi) + chng
#ifdef DEBUG
              if ( stats ) then
                ngs%statsfrz(k,j,i) = ngs%statsfrz(k,j,i) + chng
              end if
#endif
            end if
            !---------------------------------------------------------------
            ! evaporation - follows Jakob and Klein MWR 2000, with mods from
            !               Tompkins
            !------------------------------------------------------------
            ! recalculate qpretot since melting term may have changed it
            ! rprecrhmax is the threshold for the clear-sky RH that
            ! can be reached by evaporation of precipitation. This assumption
            ! is done to prevent the gridbox saturating due to the evaporation
            ! of precipitation occuring in a portion of the grid
            !------------------------------------------------------------
            qpretot = d_zero
            do n = 1 , nqx
              if ( lfall(n) ) then
                qpretot = qpretot + qxfg(n)
              end if
            end do

            zrh = rprecrhmax + (d_one-rprecrhmax)*covpclr/(d_one-ccover)
            zrh = min(max(zrh,rprecrhmax),d_one)

            ! This is a critical relative humidity that is used to limit
            ! moist environment to prevent the gridbox saturating when
            ! only part of the gridbox has evaporating precipitation
            qe = (qxfg(iqqv) - ccover*qsliq(k,j,i)) / (d_one-ccover)
            !---------------------------------------------
            ! humidity in moistest covpclr part of domain
            !---------------------------------------------
            qe = max(min(qe,qsliq(k,j,i)),d_zero)
            lactiv = covpclr > d_zero .and. &
                     covptot > d_zero .and. &
                     qpretot > d_zero .and.      &
                     qx0(iqqr) > activqx .and.   &
                     qe < zrh*qsliq(k,j,i)
            if ( lactiv ) then
              ! note: units of preclr and qpretot differ
              !       qpretot is a mixing ratio (hence "q" in name)
              !       preclr is a rain flux
              preclr = qpretot*covpclr/(covptot*dtgdp)
              !--------------------------------------
              ! actual microphysics formula in beta
              !--------------------------------------
              ! sensitivity test showed multiply rain evap rate by 0.5
              beta1 = sqrt(ph(k,j,i)/pbot)/5.09e-3_rkx*preclr/covpclr
              if ( beta1 > d_zero ) then
                beta = d_half*egrav*rpecons*(beta1)**0.5777_rkx
                denom = d_one + beta*dt*corqsliq
                dpr = covpclr * beta * (qsliq(k,j,i)-qe)/denom*dp*regrav
                dpevap = dpr*dtgdp

                !---------------------------------------------------------
                ! add evaporation term to explicit sink.
                ! this has to be explicit since if treated in the implicit
                ! term evaporation can not reduce rain to zero and model
                ! produces small amounts of rainfall everywhere.
                !---------------------------------------------------------

                ! evaporate rain
                chng = min(dpevap,qxfg(iqqr))
                chng = max(chng,d_zero)
                !-------------------------------------------------------------
                ! reduce the total precip coverage proportional to evaporation
                !-------------------------------------------------------------
                covptot = covptot - max(d_zero, &
                           (covptot-ccover)*dpevap/qpretot)
                covptot = max(covptot,rcovpmin)
              else
                chng = qxfg(iqqr)
              end if
              qsexp(iqqv,iqqr) = qsexp(iqqv,iqqr) + chng
              qsexp(iqqr,iqqv) = qsexp(iqqr,iqqv) - chng
              qxfg(iqqr)       = qxfg(iqqr) - chng
              qxfg(iqqv)       = qxfg(iqqv) + chng
#ifdef DEBUG
              if ( stats ) then
                ngs%statsrainev(k,j,i) = chng
              end if
#endif
            end if

            ! snow
            qe = (qxfg(iqqv) - ccover*qsice(k,j,i)) / (d_one-ccover)
            !---------------------------------------------
            ! humidity in moistest covpclr part of domain
            !---------------------------------------------
            qe = max(min(qe,qsice(k,j,i)),d_zero)
            lactiv = covpclr > d_zero .and. &
                     covptot > d_zero .and. &
                     qpretot > d_zero .and.      &
                     qx0(iqqs) > activqx .and.   &
                     qe < zrh*qsice(k,j,i)
            if ( lactiv ) then
              ! note: units of preclr and qpretot differ
              !       qpretot is a mixing ratio (hence "q" in name)
              !       preclr is a rain flux
              preclr = qpretot*covpclr/(covptot*dtgdp)
              !--------------------------------------
              ! actual microphysics formula in beta
              !--------------------------------------
              beta1 = sqrt(ph(k,j,i)/pbot)/5.09e-3_rkx*preclr/covpclr
              if ( beta1 >= d_zero ) then
                beta = d_half*egrav*rpecons*(beta1)**0.5777_rkx
                denom = d_one + beta*dt*corqsice
                dpr = covpclr * beta * (qsice(k,j,i)-qe)/denom*dp*regrav
                dpevap = dpr*dtgdp

                ! sublimation of  snow
                ! AMT just evaporate all if snow is very small
                if ( qxfg(iqqs) < activqx ) dpevap = qxfg(iqqs)

                chng = min(dpevap,qxfg(iqqs))
                chng = max(chng,d_zero)
                !-------------------------------------------------------------
                ! reduce the total precip coverage proportional to evaporation
                !-------------------------------------------------------------
                covptot = covptot - &
                     max(d_zero,(covptot-ccover)*dpevap/qpretot)
                covptot = max(covptot,rcovpmin)
              else
                chng = qxfg(iqqs)
              end if
              qsexp(iqqv,iqqs) = qsexp(iqqv,iqqs) + chng
              qsexp(iqqs,iqqv) = qsexp(iqqs,iqqv) - chng
              qxfg(iqqs)       = qxfg(iqqs) - chng
              qxfg(iqqv)       = qxfg(iqqv) + chng
#ifdef DEBUG
              if ( stats ) then
                ngs%statssnowev(k,j,i) = chng
              end if
#endif
            end if

          end if ! lmicro
          !------------------------------------------------------------------
          !  MICROPHYSICS ENDS HERE
          !------------------------------------------------------------------

          !--------------------------------
          ! solver for the microphysics
          !--------------------------------
          ! Truncate sum of explicit sinks to size of bin
          ! this approach is inaccurate, but conserves -
          ! prob best can do with explicit (i.e. not implicit!) terms
          !----------------------------------------------------------
          do n = 1 , nqx
            sinksum(n) = d_zero
          end do
          do n = 1 , nqx
            do m = 1 , nqx
              lind2(m,n) = .false.
            end do
          end do
          !----------------------------
          ! collect sink terms and mark
          !----------------------------
          do jn = 1 , nqx
            do n = 1 , nqx
              sinksum(n) = sinksum(n) - qsexp(n,jn)
            end do
          end do
          !---------------------------------------
          ! calculate overshoot and scaling factor
          !---------------------------------------
          do n = 1 , nqx
            ratio(n) = max(qx0(n),verylowqx) / &
              max(sinksum(n),max(qx0(n),verylowqx))
          end do
          !--------------------------------------------------------
          ! now sort ratio to find out which species run out first
          !--------------------------------------------------------
          ! iorder = argsort(ratio)
          !--------------------------------------------
          ! scale the sink terms, in the correct order,
          ! recalculating the scale factor each time
          !--------------------------------------------
          do n = 1 , nqx
            sinksum(n) = d_zero
          end do
          !----------------
          ! recalculate sum
          !----------------
          do n = 1 , nqx
            do jn = 1 , nqx
              !jo = iorder(n)
              !lind2(jo,jn) = qsexp(jo,jn) < d_zero
              !sinksum(jo) = sinksum(jo) - qsexp(jo,jn)
              lind2(n,jn) = qsexp(n,jn) < d_zero
              sinksum(n) = sinksum(n) - qsexp(n,jn)
            end do
          end do
          !---------------------------
          ! recalculate scaling factor
          !---------------------------
          do n = 1 , nqx
            !jo = iorder(n)
            !ratio(jo) = max(qx0(jo),verylowqx) / &
            !   max(sinksum(jo),max(qx0(jo),verylowqx))
            ratio(n) = max(qx0(n),verylowqx) / &
               max(sinksum(n),max(qx0(n),verylowqx))
          end do
          !------
          ! scale
          !------
          do n = 1 , nqx
            do jn = 1 , nqx
              !jo = iorder(n)
              !if ( lind2(jo,jn) ) then
              !  qsexp(jo,jn) = qsexp(jo,jn)*ratio(jo)
              !  qsexp(jn,jo) = qsexp(jn,jo)*ratio(jo)
              !end if
              if ( lind2(n,jn) ) then
                qsexp(n,jn) = qsexp(n,jn)*ratio(n)
                qsexp(jn,n) = qsexp(jn,n)*ratio(n)
              end if
            end do
          end do

          ! SOLVE THE LINEAR SYSTEM

          ! Set the LHS of equation
          do n = 1 , nqx
            do jn = 1 , nqx
              ! Diagonals: microphysical sink terms+transport
              if ( jn == n ) then
                qlhs(jn,n) = d_one + fallsink(n)
                do jo = 1 , nqx
                  qlhs(jn,n) = qlhs(jn,n) + qsimp(jo,jn)
                end do
                ! Non-diagonals: microphysical source terms
              else
                ! Here is the delta T - missing from doc.
                qlhs(jn,n) = -qsimp(jn,n)
              end if
            end do
          end do

          ! Set the RHS of equation

          do n = 1 , nqx
            ! Sum the explicit source and sink
            rexplicit = d_zero
            do jn = 1 , nqx
              ! Positive, since summed over 2nd index
              rexplicit = rexplicit + qsexp(n,jn)
            end do
            qxn(n) = qx0(n) + rexplicit
          end do

          do n = 1 , nqx
            aamax = d_zero
            do jn = 1 , nqx
              if ( abs(qlhs(n,jn)) > aamax ) aamax = abs(qlhs(n,jn))
            end do
#ifdef DEBUG
            if ( aamax == d_zero ) then
              do nn = 1 , nqx
                write(stderr,'(a,i2,f20.9)') 'QX0 ', nn , qx0(nn)
                do ll = 1 , nqx
                  write(stderr,'(a,i2,i2,f20.9)') 'QLHS ', &
                      ll , nn , qlhs(ll,nn)
                end do
              end do
            end if
#endif
            vv(n) = d_one/aamax ! Save the scaling.
          end do
          !                                                Ux=y
          ! solve A x = b-------------> LU x = b---------> Ly=b
          !
          do n = 1 , nqx
            ! This is the loop over columns
            if ( n > 1 ) then
              do m = 1 , n - 1
                xsum = qlhs(m,n)
                do kk = 1 , m - 1
                  xsum = xsum - qlhs(m,kk)*qlhs(kk,n)
                end do
                qlhs(m,n) = xsum
              end do
            end if
            ! Initialize the search for largest pivot element.
            aamax = d_zero
            imax = n
            do m = n , nqx
              xsum = qlhs(m,n)
              if ( n > 1 ) then
                do kk = 1 , n - 1
                  xsum = xsum - qlhs(m,kk)*qlhs(kk,n)
                end do
                qlhs(m,n) = xsum
              end if
              dum = vv(m)*abs(xsum)   ! Figure of merit for the pivot.
              if ( dum >= aamax ) then
                ! better than the best so far
                imax = m
                aamax = dum
              end if
            end do
            if ( n /= imax ) then
              ! Do we need to interchange rows? yes, do so...
              ! D = -D !...and change the parity of D.
              do ii = 1 , nqx
                swap = qlhs(imax,ii)
                qlhs(imax,ii) = qlhs(n,ii)
                qlhs(n,ii) = swap
              end do
              vv(imax) = vv(n) ! Also interchange the scale factor.
            end if
            indx(n) = imax
            if ( n /= nqx ) then
              dum = d_one/max(qlhs(n,n),verylowqx)
              do m = n + 1 , nqx
                qlhs(m,n) = qlhs(m,n)*dum
              end do
            end if
          end do
          !
          ! Now solve the set of n linear equations A * X = B.
          ! B(1:N) is input as the right-hand side vector B,
          ! and is used to store solution after back-substitution.
          !
          ii = 0
          ! When ii is set to a positive value, it will become
          ! the index of the  first nonvanishing element of B.
          ! We now do the forward substitution, and the only new
          ! wrinkle is to unscramble the permutation as we go.
          do m = 1 , nqx
            ll = indx(m)
            xsum = qxn(ll)
            qxn(ll) = qxn(m)
            if ( ii == 0 ) then
              if ( abs(xsum) > verylowqx ) ii = m
            else
              do jj = ii , m - 1
                xsum = xsum - qlhs(m,jj)*qxn(jj)
              end do
            end if
            qxn(m) = xsum
          end do

          ! Now we do the backsubstitution
          do m = nqx , 1 , -1
            xsum = qxn(m)
            do jj = m + 1 , nqx
              xsum = xsum - qlhs(m,jj)*qxn(jj)
            end do
            ! Store a component of the solution vector qxn.
            qxn(m) = xsum/qlhs(m,m)
          end do

          !-------------------------------------------------------------------
          !  Precipitation/sedimentation fluxes to next level
          !  diagnostic precipitation fluxes
          !  It is this scaled flux that must be used for source to next layer
          !-------------------------------------------------------------------
          do n = 1 , nqx
            chng = qxn(n) - qx0(n)
            pfplsx(n,k+1,j,i) = fallsink(n)*qxn(n)*rdtgdp
            ! Generalized precipitation flux
            ! this will be the source for the k
            ! Calculate fluxes in and out of box for conservation of TL
            fluxq = convsrce(n) + fallsrce(n) - fallsink(n)*qxn(n)
            ! Calculate the water variables tendencies
            qxtendc(n,k,j,i) = qxtendc(n,k,j,i) + chng*rdt
            ! Calculate the temperature tendencies
            if ( iphase(n) == 1 ) then
              ttendc(k,j,i) = ttendc(k,j,i)+wlhvocp*(chng-fluxq)*rdt
            else if ( iphase(n) == 2 ) then
              ttendc(k,j,i) = ttendc(k,j,i)+wlhsocp*(chng-fluxq)*rdt
            end if
          end do
        end do  ! kz : end of vertical loop
#ifdef __GFORTRAN__
      end do
#endif
    end do      ! jx, iy : end of latitude-longitude loop

    if ( idynamic == 3 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz, n = 1:nqx )
        mc2mo%qxten(j,i,k,n) = qxtendc(n,k,j,i)
      end do
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        mc2mo%tten(j,i,k) = ttendc(k,j,i)
      end do
    else
      !
      ! Couple tendencies with pressure
      !
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz, n = 1:nqx )
        mc2mo%qxten(j,i,k,n) = qxtendc(n,k,j,i)*mo2mc%psb(j,i)
      end do
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        mc2mo%tten(j,i,k) = ttendc(k,j,i)*mo2mc%psb(j,i)
      end do
    end if
    !
    !-------------------------------------
    ! Final enthalpy and total water diagnostics
    !-------------------------------------
    if ( budget_compute ) then

      ! Initialize the flux arrays
      do concurrent ( k = 1:kz, j = jci1:jci2, i = ici1:ici2 )
        sumh1(k,j,i) = d_zero
        sumq1(k,j,i) = d_zero
      end do
      do concurrent ( j = jci1:jci2, i = ici1:ici2 )
        errorq(j,i) = d_zero
        errorh(j,i) = d_zero
      end do
#ifndef __GFORTRAN__
      do concurrent ( j = jci1:jci2, i = ici1:ici2 ) &
        local(dp,tnew,qvnew,tmpl,tmpi)
#else
      do i = ici1 , ici2
        do j = jci1 , jci2
#endif
          dp = dpfs(1,j,i)
          tnew = tx(1,j,i)+dt*(ttendc(1,j,i)-tentkp(1,j,i))
          qvnew = qx(iqqv,1,j,i)+dt*(qxtendc(iqqv,1,j,i)-tenqkp(iqqv,1,j,i))
          tmpl = qx(iqql,1,j,i)+dt*(qxtendc(iqql,1,j,i)-tenqkp(iqql,1,j,i))+&
                 qx(iqqr,1,j,i)+dt*(qxtendc(iqqr,1,j,i)-tenqkp(iqqr,1,j,i))
          tmpi = qx(iqqi,1,j,i)+dt*(qxtendc(iqqi,1,j,i)-tenqkp(iqqi,1,j,i))+&
                 qx(iqqs,1,j,i)+dt*(qxtendc(iqqs,1,j,i)-tenqkp(iqqs,1,j,i))
          tnew = tnew - wlhvocp*tmpl - wlhsocp*tmpi
          sumq1(1,j,i) = sumq1(1,j,i) + (tmpl + tmpi + qvnew)*dp*regrav
          sumh1(1,j,i) = sumh1(1,j,i) + dp*tnew
#ifdef __GFORTRAN__
        end do
#endif
      end do
#ifndef __GFORTRAN__
      do concurrent ( j = jci1:jci2, i = ici1:ici2 ) &
        local(dp,tnew,qvnew,tmpl,tmpi,k)
#else
      do i = ici1 , ici2
        do j = jci1 , jci2
#endif
          do k = 2 , kz
            dp = dpfs(k,j,i)
            tnew = tx(k,j,i)+dt*(ttendc(k,j,i)-tentkp(k,j,i))
            qvnew = qx(iqqv,k,j,i)+dt*(qxtendc(iqqv,k,j,i)-tenqkp(iqqv,k,j,i))
            sumq1(k,j,i) = sumq1(k-1,j,i)
            sumh1(k,j,i) = sumh1(k-1,j,i)
            tmpl = qx(iqql,k,j,i)+dt*(qxtendc(iqql,k,j,i)-tenqkp(iqql,k,j,i))+&
                   qx(iqqr,k,j,i)+dt*(qxtendc(iqqr,k,j,i)-tenqkp(iqqr,k,j,i))
            tmpi = qx(iqqi,k,j,i)+dt*(qxtendc(iqqi,k,j,i)-tenqkp(iqqi,k,j,i))+&
                   qx(iqqs,k,j,i)+dt*(qxtendc(iqqs,k,j,i)-tenqkp(iqqs,k,j,i))
            tnew = tnew - wlhvocp*tmpl - wlhsocp*tmpi
            sumq1(k,j,i) = sumq1(k,j,i) + (tmpl + tmpi + qvnew)*dp*regrav
            sumh1(k,j,i) = sumh1(k,j,i) + dp*tnew
          end do
#ifdef __GFORTRAN__
        end do
#endif
      end do
#ifndef __GFORTRAN__
      do concurrent ( j = jci1:jci2, i = ici1:ici2 ) &
        local(dp,dtgdp,rain,rainh,k,n)
#else
      do i = ici1 , ici2
        do j = jci1 , jci2
#endif
          do k = 1 , kz
            dp = dpfs(k,j,i)
            dtgdp = dt*egrav/dp
            rain = d_zero
            rainh = d_zero
            do n = 1 , nqx
              rain = rain + dt*pfplsx(n,k+1,j,i)
              if ( iphase(n) == 1 ) then
                rainh = rainh+wlhvocp*dtgdp*pfplsx(n,k+1,j,i)*dp
              else if ( iphase(n) == 2 ) then
                rainh = rainh+wlhsocp*dtgdp*pfplsx(n,k+1,j,i)*dp
              end if
            end do
            sumq1(k,j,i) = sumq1(k,j,i) + rain
            sumh1(k,j,i) = sumh1(k,j,i) - rainh
          end do
#ifdef __GFORTRAN__
        end do
#endif
      end do
      do concurrent ( k = 1:kz, j = jci1:jci2, i = ici1:ici2 )
        sumh1(k,j,i) = sumh1(k,j,i) / pf(k+1,j,i)
      end do
#ifndef __GFORTRAN__
      do concurrent ( j = jci1:jci2, i = ici1:ici2 ) local(k)
#else
      do i = ici1 , ici2
        do j = jci1 , jci2
#endif
          do k = 1 , kz
            errorq(j,i) = errorq(j,i) + (sumq1(k,j,i)-sumq0(k,j,i))
            errorh(j,i) = errorh(j,i) + (sumh1(k,j,i)-sumh0(k,j,i))
          end do
#ifdef __GFORTRAN__
        end do
#endif
      end do

      lerror = .false.
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( abs(errorq(j,i)) > 1.e-10_rkx .or. &
               abs(errorh(j,i)) > 1.e-10_rkx) then
            if ( abs(errorq(j,i)) > 1.e-10_rkx ) then
              write(stderr,*) 'WATER NON CONSERVED AT '
              write(stderr,*) 'J = ',j
              write(stderr,*) 'I = ',i
              write(stderr,*) 'ERROR IS : ',errorq(j,i)
            end if
            if ( abs(errorh(j,i)) > 1.e-10_rkx ) then
              write(stderr,*) 'ENTHALPY NON CONSERVED AT '
              write(stderr,*) 'J = ',j
              write(stderr,*) 'I = ',i
              write(stderr,*) 'ERROR IS : ',errorh(j,i)
            end if
            lerror = .true.
          end if
        end do
      end do
      if ( lerror ) then
        call fatal(__FILE__,__LINE__, &
                'TOTAL WATER OR ENTHALPY NOT CONSERVED')
      end if
    end if ! budget_compute

    ! Sum fluxes over the levels
    ! Initialize fluxes
    do concurrent ( k = 1:kzp1, j = jci1:jci2, i = ici1:ici2 )
      pfplsl(k,j,i) = d_zero
      pfplsn(k,j,i) = d_zero
    end do

    !--------------------------------------------------------------------
    ! Copy general precip arrays back into FP arrays
    ! Add rain and liquid fluxes, ice and snow fluxes
    !--------------------------------------------------------------------

    ! Rain+liquid, snow+ice
    ! for each level k = 1 , kzp1, sum of the same phase elements
    do n = 1 , nqx
      do concurrent ( k = 1:kzp1, j = jci1:jci2, i = ici1:ici2 )
        if ( iphase(n) == 1 ) then
          pfplsl(k,j,i) = pfplsl(k,j,i) + pfplsx(n,k,j,i)
        else if ( iphase(n) == 2 ) then
          pfplsn(k,j,i) = pfplsn(k,j,i) + pfplsx(n,k,j,i)
        end if
      end do
    end do
    !
    if ( ichem == 1 ) then
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        mc2mo%rainls(j,i,k) = pfplsl(k+1,j,i)
        ! save the 3D precip for chemical washout
        mc2mo%rembc(j,i,k) =  mc2mo%rainls(j,i,k)
        mc2mo%remrat(j,i,k) = remrat(k,j,i)
      end do
    end if
    !--------------------------------------------------------------
    ! Convert the accumlated precipitation to appropriate units for
    ! the surface physics and the output sum up through the levels
    !--------------------------------------------------------------
#ifndef __GFORTRAN__
    do concurrent ( j = jci1:jci2, i = ici1:ici2 ) local(prainx,psnowx)
#else
    do i = ici1 , ici2
      do j = jci1 , jci2
#endif
        prainx = pfplsl(kzp1,j,i)*dt
        psnowx = pfplsn(kzp1,j,i)*dt
        if ( prainx > d_zero ) then
          mc2mo%rainnc(j,i) = mc2mo%rainnc(j,i) + prainx
          mc2mo%lsmrnc(j,i) = mc2mo%lsmrnc(j,i) + pfplsl(kzp1,j,i)
          mc2mo%trrate(j,i) = pfplsl(kzp1,j,i)
        end if
        if ( psnowx > d_zero ) then
          mc2mo%rainnc(j,i) = mc2mo%rainnc(j,i) + psnowx
          mc2mo%lsmrnc(j,i) = mc2mo%lsmrnc(j,i) + pfplsn(kzp1,j,i)
          mc2mo%trrate(j,i) = mc2mo%trrate(j,i) + pfplsn(kzp1,j,i)
          mc2mo%snownc(j,i) = mc2mo%snownc(j,i) + pfplsn(kzp1,j,i)
        end if
#ifdef __GFORTRAN__
      end do
#endif
    end do

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif

    contains

    pure real(rkx) function edem(t,phase)
!$acc routine seq
      implicit none
      real(rkx) , intent(in):: t , phase
      edem = phase * c5alvcp * (d_one/(t-c4les)**2) + &
               (d_one - phase) * c5alscp * (d_one/(t-c4ies)**2)
    end function edem

    pure real(rkx) function eldcpm(t)
!$acc routine seq
      implicit none
      real(rkx) , intent(in):: t
      real(rkx) :: phase
      phase = max(min(d_one,((max(rtice,min(tzero,t))-rtice)* &
                              rtwat_rtice_r)**2),d_zero)
      eldcpm = phase*wlhvocp + (d_one-phase)*wlhsocp
    end function eldcpm

    pure real(rkx) function eewm(t,phase)
!$acc routine seq
      implicit none
      real(rkx) , intent(in) :: t , phase
      real(rkx) :: eliq , eice
      eliq = c2es*exp(c3les*((t-tzero)/(t-c4les)))
      eice = c2es*exp(c3ies*((t-tzero)/(t-c4ies)))
      eewm = phase * eliq + (d_one-phase) * eice
    end function eewm

   ! subroutine addpath(src,snk,proc,zsqa,zsqb,beta,fg)
   !   implicit none
   !   real(rkx) , pointer , intent(inout) , dimension(:,:) :: zsqa , zsqb
   !   real(rkx) , pointer , intent(inout) , dimension(:) :: fg
   !   real(rkx) , intent(in) :: proc
   !   integer(ik4) , intent(in) :: src , snk
   !   real(rkx) , intent(in) :: beta
   !   zsqa(src,snk) = zsqa(src,snk) + (d_one-beta)*proc
   !   zsqa(snk,src) = zsqa(snk,src) - (d_one-beta)*proc
   !   fg(src) = fg(src) + (d_one-beta)*proc
   !   fg(snk) = fg(snk) - (d_one-beta)*proc
   !   zsqb(src,snk) = zsqb(src,snk) + beta*proc
   ! end subroutine addpath

   ! pure function argsort(a) result(b)
   !   implicit none
   !   real(rk8) , intent(in) :: a(:)
   !   integer(ik4) , dimension(size(a)) :: b
   !   integer(ik4) :: n , i , imin , temp1
   !   real(rk8) :: temp2
   !   real(rk8) , dimension(size(a)) :: a2
   !   a2 = a
   !   n = size(a)
   !   do i = 1 , n
   !     b(i) = i
   !   end do
   !   if ( n == 1 ) return
   !   do i = 1 , n-1
   !     imin = minloc(a2(i:),1) + i - 1
   !     if ( imin /= i ) then
   !       temp2 = a2(i)
   !       a2(i) = a2(imin)
   !       a2(imin) = temp2
   !       temp1 = b(i)
   !       b(i) = b(imin)
   !       b(imin) = temp1
   !     end if
   !   end do
   ! end function argsort

  end subroutine nogtom

end module mod_micro_nogtom

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

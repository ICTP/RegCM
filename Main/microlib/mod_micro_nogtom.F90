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
  use mod_regcm_types
  use mod_constants
  use mod_service
  use mod_runparams , only : nqx
  use mod_runparams , only : iqqv => iqv !vapor
  use mod_runparams , only : iqql => iqc !liquid
  use mod_runparams , only : iqqr => iqr !rain
  use mod_runparams , only : iqqi => iqi !ice
  use mod_runparams , only : iqqs => iqs !snow
  use mod_runparams , only : sigma
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
  real(rkx) , parameter :: rhcrit_lnd = 0.80_rkx
  real(rkx) , parameter :: rhcrit_sea = 0.90_rkx
  !------------------------------------------------
  real(rkx) , parameter :: ciden13 = 8.87_rkx      ! ice density 700**0.333
  real(rkx) , parameter :: airconduct = 2.4e-2_rkx ! conductivity of air

  public :: allocate_mod_nogtom , init_nogtom , nogtom

  ! Total water and enthalpy budget diagnostics variables
  ! marker for water phase of each species
  ! 0 = vapour, 1 = liquid, 2 = ice
  integer(ik4) , pointer , dimension(:) :: iphase
  ! marks melting linkage for ice categories
  ! ice->liquid, snow->rain
  integer(ik4) , pointer , dimension(:) :: imelt
  ! array for sorting explicit terms
  integer(ik4) , pointer , dimension(:) :: iorder
  logical , pointer , dimension(:) :: lfall
  logical , pointer , dimension(:,:) :: lind2

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
  ! for sedimentation source/sink terms
  real(rkx) , pointer , dimension(:) :: fallsrce
  real(rkx) , pointer , dimension(:) :: fallsink
  ! for convection detrainment source and subsidence source/sink terms
  real(rkx) , pointer , dimension(:) :: convsrce
  real(rkx) , pointer , dimension(:,:,:) :: eewmt
  ! fluxes convergence of species
  real(rkx) , pointer , dimension(:,:,:) :: qliq

  real(rkx) , pointer , dimension(:) :: ratio
  real(rkx) , pointer , dimension(:) :: sinksum
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
  ! Flux of liquid
  real(rkx) , pointer , dimension(:,:,:) :: pfsqlf
  ! Flux of ice
  real(rkx) , pointer , dimension(:,:,:) :: pfsqif
  ! decoupled temperature tendency
  real(rkx) , pointer , dimension(:,:,:) :: ttendc
  ! critical factors
  real(rkx) , pointer , dimension(:,:) :: xlcrit
  real(rkx) , pointer , dimension(:,:) :: rhcrit
  ! Cloud coverage and clearsky portion
  real(rkx) , pointer , dimension(:,:) :: covptot , covpclr
  ! fall speeds of three categories
  real(rkx) , pointer , dimension(:) :: vqx
  ! n x n matrix storing the LHS of implicit solver
  real(rkx) , pointer , dimension(:,:) :: qlhs
  ! explicit sources and sinks "q s exp"=q source explicit
  real(rkx) , pointer , dimension(:,:) :: qsexp
  ! implicit sources and sinks "q s imp"=q source/sink implicit
  real(rkx) , pointer , dimension(:,:) :: qsimp
  ! decoupled mixing ratios tendency
  real(rkx) , pointer , dimension(:,:,:,:) :: qxtendc
  ! j,i,n ! generalized precipitation flux
  real(rkx) , pointer , dimension(:,:,:,:) :: pfplsx
  real(rkx) , pointer, dimension(:,:,:,:) :: qx
  real(rkx) , pointer, dimension(:,:,:) :: tx
  ! Initial values
  real(rkx) , pointer, dimension(:) :: qx0
  real(rkx) , pointer, dimension(:) :: qxfg
  ! new values for qxx at time+1
  real(rkx) , pointer, dimension(:) :: qxn
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

  integer(ik4) , pointer , dimension(:) :: indx
  real(rkx) , pointer , dimension(:) :: vv

  real(rkx) , parameter :: activqx = 1.0e-12_rkx
  real(rkx) , parameter :: zerocf = 0.0001_rkx
  real(rkx) , parameter :: onecf  = 0.9999_rkx

  abstract interface
    subroutine voidsub
      implicit none
    end subroutine voidsub
  end interface

  contains

  subroutine allocate_mod_nogtom
    implicit none
    call getmem1d(vqx,1,nqx,'cmicro:vqx')
    call getmem1d(indx,1,nqx,'cmicro:indx')
    call getmem1d(vv,1,nqx,'cmicro:vv')
    call getmem1d(imelt,1,nqx,'cmicro:imelt')
    call getmem1d(lfall,1,nqx,'cmicro:lfall')
    call getmem1d(iphase,1,nqx,'cmicro:iphase')
    call getmem3d(qliq,jci1,jci2,ici1,ici2,1,kzp1,'cmicro:qliq')
    call getmem3d(eewmt,jci1,jci2,ici1,ici2,1,kz,'cmicro:eewmt')
    call getmem3d(qsmix,jci1,jci2,ici1,ici2,1,kz,'cmicro:qsmix')
    call getmem1d(iorder,1,nqx,'cmicro:iorder')
    call getmem3d(ttendc,jci1,jci2,ici1,ici2,1,kz,'cmicro:ttendc')
    call getmem1d(convsrce,1,nqx,'cmicro:convsrce')
    call getmem3d(eew,jci1,jci2,ici1,ici2,1,kz,'cmicro:eew')
    call getmem3d(qsice,jci1,jci2,ici1,ici2,1,kz,'cmicro:qsice')
    call getmem4d(qx,1,nqx,jci1,jci2,ici1,ici2,1,kz,'cmicro:qx')
    call getmem3d(tx,jci1,jci2,ici1,ici2,1,kz,'cmicro:tx')
    call getmem3d(qsliq,jci1,jci2,ici1,ici2,1,kz,'cmicro:qsliq')
    call getmem1d(fallsink,1,nqx,'cmicro:fallsink')
    call getmem1d(fallsrce,1,nqx,'cmicro:fallsrce')
    call getmem1d(ratio,1,nqx,'cmicro:ratio')
    call getmem1d(sinksum,1,nqx,'cmicro:sinksum')
    call getmem3d(cldtopdist,jci1,jci2,ici1,ici2,1,kz,'cmicro:cldtopdist')
    call getmem3d(dqsatdt,jci1,jci2,ici1,ici2,1,kz,'cmicro:dqsatdt')
    call getmem3d(pfplsl,jci1,jci2,ici1,ici2,1,kzp1,'cmicro:pfplsl')
    call getmem3d(pfplsn,jci1,jci2,ici1,ici2,1,kzp1,'cmicro:pfplsn')
    call getmem3d(pfsqlf,jci1,jci2,ici1,ici2,1,kzp1,'cmicro:pfsqlf')
    call getmem3d(pfsqif,jci1,jci2,ici1,ici2,1,kzp1,'cmicro:pfsqif')
    call getmem3d(koop,jci1,jci2,ici1,ici2,1,kz,'cmicro:koop')
    call getmem2d(xlcrit,jci1,jci2,ici1,ici2,'cmicro:xlcrit')
    call getmem2d(rhcrit,jci1,jci2,ici1,ici2,'cmicro:rhcrit')
    call getmem2d(covptot,jci1,jci2,ici1,ici2,'cmicro:covptot')
    call getmem2d(covpclr,jci1,jci2,ici1,ici2,'cmicro:covpclr')
    call getmem3d(eeliq,jci1,jci2,ici1,ici2,1,kz,'cmicro:eeliq')
    call getmem3d(eeice,jci1,jci2,ici1,ici2,1,kz,'cmicro:eeice')
    call getmem3d(eeliqt,jci1,jci2,ici1,ici2,1,kz,'cmicro:eeliqt')
    call getmem4d(qxtendc,1,nqx,jci1,jci2,ici1,ici2,1,kz,'cmicro:qxtendc')
    call getmem1d(qx0,1,nqx,'cmicro:qx0')
    call getmem1d(qxfg,1,nqx,'cmicro:qxfg')
    call getmem1d(qxn,1,nqx,'cmicro:qxn')
    call getmem2d(qlhs,1,nqx,1,nqx,'cmicro:qlhs')
    call getmem2d(qsexp,1,nqx,1,nqx,'cmicro:qsexp')
    call getmem2d(qsimp,1,nqx,1,nqx,'cmicro:qsimp')
    call getmem2d(lind2,1,nqx,1,nqx,'cmicro:lind2')
    call getmem4d(pfplsx,1,nqx,jci1,jci2,ici1,ici2,1,kzp1,'cmicro:pfplsx')
    call getmem3d(dpfs,jci1,jci2,ici1,ici2,1,kz,'cmicro:dpfs')
    if ( budget_compute ) then
      call getmem3d(sumq0,jci1,jci2,ici1,ici2,1,kz,'cmicro:sumq0')
      call getmem3d(sumh0,jci1,jci2,ici1,ici2,1,kz,'cmicro:sumh0')
      call getmem3d(sumq1,jci1,jci2,ici1,ici2,1,kz,'cmicro:sumq1')
      call getmem3d(sumh1,jci1,jci2,ici1,ici2,1,kz,'cmicro:sumh1')
      call getmem2d(errorq,jci1,jci2,ici1,ici2,'cmicro:errorq')
      call getmem2d(errorh,jci1,jci2,ici1,ici2,'cmicro:errorh')
      call getmem4d(tenqkp,1,nqx,jci1,jci2,ici1,ici2,1,kz,'cmicro:tenqkp')
      call getmem3d(tentkp,jci1,jci2,ici1,ici2,1,kz,'cmicro:tentkp')
    end if
  end subroutine allocate_mod_nogtom

  subroutine init_nogtom(ldmsk)
    use mod_runparams , only : vfqr , vfqi , vfqs
    implicit none
    integer , pointer , dimension(:,:) , intent(in) :: ldmsk
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
    do i = ici1 , ici2
      do j = jci1 , jci2
        if ( ldmsk(j,i) == 1 ) then ! landmask =1 land
          xlcrit(j,i) = rclcrit_land ! landrclcrit_land = 5.e-4
          rhcrit(j,i) = rhcrit_lnd
        else
          xlcrit(j,i) = rclcrit_sea  ! oceanrclcrit_sea  = 3.e-4
          rhcrit(j,i) = rhcrit_sea
        end if
      end do
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
    logical :: lactiv , ltkgt0 , ltklt0 , ltkgthomo , lcloud
    logical :: locast , lconden , lccn , lerror
    logical :: ldetr
    real(rkx) :: rexplicit , xlcondlim
    real(rkx) :: facl , faci , facw , corr , gdp , acond , zdl , infactor
    real(rkx) :: alfaw , phases , zdelta , tmpl , qexc , rhc , zsig , &
                 tmpi , tnew , qvnew , qe , rain , rainh , preclr , arg
    real(rkx) :: totcond ! total condensate liquid+ice
    ! total rain frac: fractional occurence of precipitation (%)
    ! for condensation
    ! ice nuclei concentration
    ! local real variables for autoconversion rate constants
    real(rkx) :: alpha1 ! coefficient autoconversion cold cloud
    real(rkx) :: tmpa
    ! real(rkx) :: zqadj
    real(rkx) :: zrh
    real(rkx) :: beta , beta1
    ! local variables for condensation
    real(rkx) :: cond , dtdp , cdmax
    ! local variables for melting
    real(rkx) :: tdiff
    real(rkx) :: cons1
    ! constant for converting the fluxes unit measures
    real(rkx) :: prainx , psnowx
    ! local real constants for evaporation
    real(rkx) :: dpr , denom , dpevap , evapi , evapl , excess
    real(rkx) :: dqsmixdt , dqsicedt , dqsliqdt
    real(rkx) :: dp , dtgdp , rdtgdp
    real(rkx) :: corqsliq , corqsice , corqsmix , evaplimmix
    real(rkx) :: ql_incld , qi_incld , qli_incld
    real(rkx) :: supsat , subsat
    real(rkx) :: ldifdt , sink
    ! real(rkx) :: botm , rm
    real(rkx) :: qold , tcond , dqs
    real(rkx) :: chng , chngmax
    real(rkx) :: icenuclei
    real(rkx) :: qpretot
    real(rkx) :: qicetot
    real(rkx) :: ldefr
    real(rkx) :: critauto
    real(rkx) :: qliqfrac
    real(rkx) :: qicefrac
    real(rkx) :: fluxq
    ! constants for deposition process
    real(rkx) :: vpice , vpliq , xadd , xbdd , cvds , &
                 qice0 , qinew , rainaut , snowaut
    ! constants for condensation and turbulent mixing erosion of clouds
    real(rkx) :: dpmxdt , wtot , dtdiab , dtforc , &
                 qp , qsat , cond1 , levap , leros
    real(rkx) :: sqmix , ccover , lccover
    real(rkx) :: tk , tc , dens , pbot , ccn
    real(rkx) :: snowp , rainp

#ifndef __PGI
    procedure (voidsub) , pointer :: selautoconv => null()
    procedure (voidsub) , pointer :: selnss => null()
#endif

#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'microphys'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    lccn = ( ichem == 1 .and. iaerosol == 1 .and. iindirect == 2 )

#ifndef __PGI
    !---------------------------------------------------------------
    !                         AUTOCONVERSION
    !---------------------------------------------------------------
    ! Warm clouds
    select case (iautoconv)
      case (1) ! Klein & Pincus (2000)
        selautoconv => klein_and_pincus
      case (2) ! Khairoutdinov and Kogan (2000)
        selautoconv => khairoutdinov_and_kogan
      case (3) ! Kessler(1969)
        selautoconv => kessler
      case (4) ! Sundqvist
        selautoconv => sundqvist
      case default
        call fatal(__FILE__,__LINE__,'UNKNOWN AUTOCONVERSION SCHEME')
    end select
    select case(nssopt)
      case(0,1)
        selnss => nss_tompkins
      case(2)
        selnss => nss_lohmann_and_karcher
      case(3)
        selnss => nss_gierens
      case default
        call fatal(__FILE__,__LINE__, 'NSSOPT IN CLOUD MUST BE IN RANGE 0-3')
    end select
#endif

    if ( idynamic == 3 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            do n = 1 , nqx
              qxtendc(n,j,i,k) = mc2mo%qxten(j,i,k,n)
            end do
          end do
        end do
      end do
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            ttendc(j,i,k) = mc2mo%tten(j,i,k)
          end do
        end do
      end do
    else
      ! Decouple tendencies
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            do n = 1 , nqx
              qxtendc(n,j,i,k) = mc2mo%qxten(j,i,k,n) / mo2mc%psb(j,i)
            end do
          end do
        end do
      end do
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            ttendc(j,i,k) = mc2mo%tten(j,i,k) / mo2mc%psb(j,i)
          end do
        end do
      end do
    end if

    ! Define the initial array qx
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          do n = 1 , nqx
            qx(n,j,i,k) = mo2mc%qxx(j,i,k,n)
          end do
        end do
      end do
    end do

    ! Define the initial array qx
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          tx(j,i,k) = mo2mc%t(j,i,k)
        end do
      end do
    end do

    ! Delta pressure
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          dpfs(j,i,k) = mo2mc%pfs(j,i,k+1)-mo2mc%pfs(j,i,k)
        end do
      end do
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

    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          qliq(j,i,k) = max(min(d_one,((max(rtice,min(tzero, &
                        tx(j,i,k)))-rtice)*rtwat_rtice_r)**2),d_zero)
        end do
      end do
    end do

    ! Reset total precipitation variables
    pfplsx(:,:,:,:) = d_zero

    ! Compute supersaturations
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          eeliq(j,i,k) = c2es*exp(c3les*((tx(j,i,k)-tzero)/(tx(j,i,k)-c4les)))
          eeice(j,i,k) = c2es*exp(c3ies*((tx(j,i,k)-tzero)/(tx(j,i,k)-c4ies)))
          koop(j,i,k) = min(rkoop1-rkoop2*tx(j,i,k),eeliq(j,i,k)/eeice(j,i,k))
        end do
      end do
    end do

    !-------------------------------------
    ! Initial enthalpy and total water diagnostics
    !-------------------------------------
    !
    ! Starting budget if requested
    !
    if ( budget_compute ) then

      ! Reset arrays
      tentkp(:,:,:)  = d_zero
      tenqkp(:,:,:,:) = d_zero

      ! Record the tendencies
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            do n = 1 , nqx
              tenqkp(n,j,i,k) = qxtendc(n,j,i,k)
            end do
          end do
        end do
      end do
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            tentkp(j,i,k) = ttendc(j,i,k)
          end do
        end do
      end do

      ! initialize the flux arrays
      sumq0(:,:,:)     = d_zero
      sumh0(:,:,:)     = d_zero

      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            tnew = tx(j,i,k)
            dp = dpfs(j,i,k)
            qe = mo2mc%qdetr(j,i,k)

            if ( k > 1 ) then
              sumq0(j,i,k) = sumq0(j,i,k-1) ! total water
              sumh0(j,i,k) = sumh0(j,i,k-1) ! liquid water temperature
            end if

            tmpl = qx(iqql,j,i,k)+qx(iqqr,j,i,k)
            tmpi = qx(iqqi,j,i,k)+qx(iqqs,j,i,k)
            tnew = tnew - wlhvocp*tmpl - wlhsocp*tmpi
            sumq0(j,i,k) = sumq0(j,i,k)+(tmpl+tmpi+qx(iqqv,j,i,k))*dp*regrav

            ! Detrained water treated here
            if ( lmicro .and. abs(qe) > activqx ) then
              sumq0(j,i,k) = sumq0(j,i,k) + qe*dp*regrav
              alfaw = qliq(j,i,k)
              tnew = tnew-(wlhvocp*alfaw+wlhsocp*(d_one-alfaw))*qe
            end if
            sumh0(j,i,k) = sumh0(j,i,k) + dp*tnew
          end do
        end do
      end do
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            sumh0(j,i,k) = sumh0(j,i,k)/mo2mc%pfs(j,i,k+1)
          end do
        end do
      end do
    end if ! budget_compute

    ! -------------------------------
    ! Define saturation values
    !---------------------------
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          ! zdelta = 1 if t > tzero
          ! zdelta = 0 if t < tzero
          zdelta = max(d_zero,sign(d_one,tx(j,i,k)-tzero))
          !---------------------------------------------
          ! mixed phase saturation
          !--------------------------------------------
          phases = qliq(j,i,k)
          eewmt(j,i,k) = eeliq(j,i,k)*phases + eeice(j,i,k)*(d_one-phases)
          eewmt(j,i,k) = min(eewmt(j,i,k)/mo2mc%phs(j,i,k),d_half)
          qsmix(j,i,k) = eewmt(j,i,k)
          ! ep1 = rwat/rgas - d_one
          qsmix(j,i,k) = qsmix(j,i,k)/(d_one-ep1*qsmix(j,i,k))
          !--------------------------------------------
          ! ice saturation T < 273K
          ! liquid water saturation for T > 273K
          !--------------------------------------------
          eew(j,i,k) = (zdelta*eeliq(j,i,k) + &
               (d_one-zdelta)*eeice(j,i,k))/mo2mc%phs(j,i,k)
          eew(j,i,k) = min(d_half,eew(j,i,k))
          !ice water saturation
          qsice(j,i,k) = min(eeice(j,i,k)/mo2mc%phs(j,i,k),d_half)
          qsice(j,i,k) = qsice(j,i,k)/(d_one-ep1*qsice(j,i,k))
          !----------------------------------
          ! liquid water saturation
          !----------------------------------
          !eeliq is the saturation vapor pressure es(T)
          !the saturation mixing ratio is ws = es(T)/p *0.622
          !ws = ws/(-(d_one/eps - d_one)*ws)
          eeliqt(j,i,k) = min(eeliq(j,i,k)/mo2mc%phs(j,i,k),d_half)
          qsliq(j,i,k) = eeliqt(j,i,k)
          qsliq(j,i,k) = qsliq(j,i,k)/(d_one-ep1*qsliq(j,i,k))
        end do
      end do
    end do

    !--------------------------------ADEED BY RITA
    ! Calculate distance from cloud top
    ! defined by cloudy layer below a layer with cloud frac <0.01
    !--------------------------------------------------------------

    cldtopdist(:,:,:) = d_zero
    do k = 2 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          do kk = 2 , k
            if ( mc2mo%fcc(j,i,kk-1) > cldtopcf .and. &
                 mc2mo%fcc(j,i,kk)  <= cldtopcf ) then
              cldtopdist(j,i,k) = cldtopdist(j,i,k) + mo2mc%delz(j,i,kk)
            end if
          end do
        end do
      end do
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
    !                       INITIALIZE STORAGE
    !----------------------------------------------------------------------
    !
    covptot(:,:) = d_zero
    covpclr(:,:) = d_zero
    !
    !----------------------------------------------------------------------
    !                       START OF VERTICAL LOOP
    !----------------------------------------------------------------------
    !
    ! Loop over levels and points
    !
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2

          supsat      = d_zero
          subsat      = d_zero
          fallsrce(:) = d_zero
          fallsink(:) = d_zero
          convsrce(:) = d_zero
          ldefr       = d_zero

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
          qsexp(:,:)  = d_zero
          qsimp(:,:)  = d_zero
          !
          !---------------------------------
          ! First guess microphysics
          !---------------------------------
          do n = 1 , nqx
            qx0(n)  = qx(n,j,i,k)
            qxfg(n) = qx0(n)
          end do

          ldetr = ( abs(mo2mc%qdetr(j,i,k)) > activqx )
          totcond = qxfg(iqql)+qxfg(iqqi)
          lconden = ( totcond > 2.0*activqx )
          if ( lconden ) then
            qliqfrac = qxfg(iqql)/totcond
            qicefrac = d_one-qliqfrac
          else
            qliqfrac = d_zero
            qicefrac = d_zero
          end if

          qicetot = d_zero
          do n = 1 , nqx
            if ( iphase(n) == 2 ) then
              qicetot = qicetot + qxfg(n)
            end if
          end do

          critauto = xlcrit(j,i)
          pbot     = mo2mc%pfs(j,i,kzp1)
          dp       = dpfs(j,i,k)
          tk       = tx(j,i,k)
          tc       = tk - tzero
          dens     = mo2mc%rho(j,i,k)
          sqmix    = qsmix(j,i,k)
          ccover   = mc2mo%fcc(j,i,k)
          ccover   = min(max(ccover,zerocf),onecf)

          if ( k == 1 ) then
            lccover = d_zero
            rainp   = d_zero
            snowp   = d_zero
          else
            lccover = mc2mo%fcc(j,i,k-1)
            lccover = min(max(lccover,zerocf),onecf)
            rainp   = pfplsx(iqqr,j,i,k)
            snowp   = pfplsx(iqqs,j,i,k)
          end if

          if ( lccn ) ccn = mo2mc%ccn(j,i,k)

          ltkgt0    = ( tk > tzero )
          ltklt0    = ( .not. ltkgt0 )
          ltkgthomo = ( tk > thomo )
          lcloud    = ( ccover > zerocf )
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
          corr     = d_one/(d_one - ep1*eeliqt(j,i,k))
          dqsliqdt = facw*corr*qsliq(j,i,k)
          corqsliq = d_one + wlhvocp*dqsliqdt
          ! ice
          faci     = c5ies/((tk - c4ies)**2)
          corr     = d_one/(d_one - ep1*eew(j,i,k))
          dqsicedt = faci*corr*qsice(j,i,k)
          corqsice = d_one + wlhsocp*dqsicedt
          ! diagnostic mixed
          alfaw    = qliq(j,i,k)
          facl     = alfaw*facw + (d_one - alfaw)*faci
          corr     = d_one/(d_one - ep1*eewmt(j,i,k))
          dqsmixdt = facl*corr*sqmix
          corqsmix = d_one/(d_one + eldcpm(tk)*dqsmixdt)
          !--------------------------------
          ! evaporation/sublimation limits
          !--------------------------------
          evaplimmix = max((sqmix-qxfg(iqqv))*corqsmix,d_zero)

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
                  fallsrce(n) = pfplsx(n,j,i,k)*dtgdp
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

            if ( qx0(iqql) < activqx ) then
              qsexp(iqqv,iqql) =  qx0(iqql)
              qsexp(iqql,iqqv) = -qx0(iqql)
              qxfg(iqql) = qxfg(iqql) - qx0(iqql)
              qxfg(iqqv) = qxfg(iqql) + qx0(iqql)
            end if
            if ( qx0(iqqi) < activqx ) then
              qsexp(iqqv,iqqi) =  qx0(iqqi)
              qsexp(iqqi,iqqv) = -qx0(iqqi)
              qxfg(iqqi) = qxfg(iqqi) - qx0(iqqi)
              qxfg(iqqv) = qxfg(iqqi) + qx0(iqqi)
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

            !-----------..........--------------------------------------------
            !  ICE SUPERSATURATION ADJUSTMENT
            !-..........------------------------------------------------------
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
                facl = ccover + koop(j,i,k)*(d_one-ccover)
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
            supsat = max((qxfg(iqqv)-facl*sqmix)*corqsmix,d_zero)
            ! e < esi, because for e > esi ice still present
            subsat = min((qxfg(iqqv)-facl*sqmix)*corqsmix,d_zero)
            if ( supsat > dlowval ) then
              if ( ltkgthomo ) then
                ! turn supersaturation into liquid water
                qsexp(iqql,iqqv) = qsexp(iqql,iqqv) + supsat
                qsexp(iqqv,iqql) = qsexp(iqqv,iqql) - supsat
                qxfg(iqql) = qxfg(iqql) + supsat
                qxfg(iqqv) = qxfg(iqqv) - supsat
#ifdef DEBUG
                if ( stats ) then
                  ngs%statssupw(j,i,k) = ngs%statssupw(j,i,k) + supsat
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
                  ngs%statssupc(j,i,k) = ngs%statssupc(j,i,k) - evapi
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
              qe = mo2mc%qdetr(j,i,k)
              alfaw = qliq(j,i,k)
              convsrce(iqql) = alfaw*qe
              convsrce(iqqi) = (d_one-alfaw)*qe
              qsexp(iqql,iqql) = qsexp(iqql,iqql) + convsrce(iqql)
              qsexp(iqqi,iqqi) = qsexp(iqqi,iqqi) + convsrce(iqqi)
              qxfg(iqql) = qxfg(iqql) + convsrce(iqql)
              qxfg(iqqi) = qxfg(iqqi) + convsrce(iqqi)
#ifdef DEBUG
              if ( stats ) then
                ngs%statsdetrw(j,i,k) = convsrce(iqql)
                ngs%statsdetrc(j,i,k) = convsrce(iqqi)
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
              leros = ccover * ldifdt * max(sqmix-qxfg(iqqv),d_zero)
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
                ngs%statserosw(j,i,k) = qliqfrac*leros
                ngs%statserosc(j,i,k) = qicefrac*leros
              end if
#endif
            end if

            !------------------------------------------------------------------
            ! condensation/evaporation due to dqsat/dt
            !------------------------------------------------------------------
            ! calculate dqs/dt and use to calculate the cloud source
            ! note that old diagnostic mix phased qsat is retained for moment
            !------------------------------------------------------------------
            dtdp   = rovcp*tk/mo2mc%phs(j,i,k)
            dpmxdt = dp*rdt
            wtot   = mo2mc%pverv(j,i,k)
            wtot   = min(dpmxdt,max(-dpmxdt,wtot))
            dtdiab = min(dpmxdt*dtdp, &
                     max(-dpmxdt*dtdp,mo2mc%heatrt(j,i,k)))*dt+wlhfocp*ldefr
            ! ldefr = 0
            ! note: ldefr should be set to the difference between the mixed
            ! phase functions in the convection and cloud scheme, and
            ! for now we set it to zero and the functions are the same.
            ! In RegCM not all convection schemes provide such info.
            dtforc = dtdp*wtot*dt + dtdiab
            qold   = sqmix
            tcond  = tk + dtforc
            tcond  = max(tcond,160.0_rkx)
            ! the goal is to produce dqs = qsmix - qold, where qsmix is
            ! reduced because of the condensation. so that dqs is negative?
            qp = d_one/mo2mc%phs(j,i,k)
            phases = max(min(d_one,((max(rtice,min(tzero, &
                       tcond))-rtice)*rtwat_rtice_r)**2),d_zero)
            ! saturation mixing ratio ws
            qsat = eewm(tcond,phases) * qp
            qsat = min(d_half,qsat)          ! ws < 0.5        WHY?
            corr  = d_one/(d_one-ep1*qsat)
            qsat = qsat*corr
            cond = (sqmix-qsat)/(d_one + qsat*edem(tcond,phases))
            tcond = tcond + eldcpm(tcond)*cond
            phases = max(min(d_one,((max(rtice,min(tzero, &
                       tcond))-rtice)*rtwat_rtice_r)**2),d_zero)
            sqmix = sqmix - cond
            qsat = eewm(tcond,phases) * qp
            qsat = min(d_half,qsat)
            corr = d_one/(d_one-ep1*qsat)
            qsat = qsat*corr
            cond1 = (sqmix-qsat)/(d_one + qsat*edem(tcond,phases))
            tcond = tcond + eldcpm(tcond)*cond1
            sqmix = sqmix - cond1
            dqs = sqmix - qold
            sqmix = qold

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
              levap = min(levap,max(sqmix-qxfg(iqqv),d_zero))
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
                ngs%statsevapw(j,i,k) = qliqfrac*levap
                ngs%statsevapc(j,i,k) = qicefrac*levap
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
                  corr = d_one/(d_one-ep1*sqmix)
                  cdmax = (qxfg(iqqv)-sqmix)/(d_one+corr*sqmix*edem(tk,alfaw))
                else
                  cdmax = (qxfg(iqqv)-ccover*sqmix)/ccover
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
                    ngs%statscond1w(j,i,k) = chng
                  end if
#endif
                else if ( ltklt0 ) then
                  qsexp(iqqi,iqqv) = qsexp(iqqi,iqqv) + chng
                  qsexp(iqqv,iqqi) = qsexp(iqqv,iqqi) - chng
                  qxfg(iqqi) = qxfg(iqqi) + chng
                  qxfg(iqqv) = qxfg(iqqv) - chng
#ifdef DEBUG
                  if ( stats ) then
                    ngs%statscond1c(j,i,k) = chng
                  end if
#endif
                end if
              else
                ! (2) generation of new clouds (dc/dt>0)
#ifdef __PGI
                select case (nssopt)
                  case (0,1)
                    call nss_tompkins
                  case (2) ! Khairoutdinov and Kogan (2000)
                    call nss_lohmann_and_karcher
                  case (3) ! Kessler(1969)
                    call nss_gierens
                end select
#else
                call selnss
#endif
                rhc = rhcrit(j,i)
                zsig = mo2mc%phs(j,i,k)/pbot
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
                  facl = koop(j,i,k)
                end if
                if ( qexc >= rhc*sqmix*facl .and. qexc < sqmix*facl ) then
                  ! note: not **2 on 1-a term if qe is used.
                  ! added correction term fac to numerator 15/03/2010
                  acond = -(d_one-ccover)*facl*dqs / &
                          max(d_two*(facl*sqmix-qexc),dlowval)
                  acond = min(acond,d_one-ccover) ! put the limiter back
                  ! linear term:
                  ! added correction term fac 15/03/2010
                  chng = -facl*dqs*d_half*acond !mine linear
                  ! new limiter formulation
                  ! qsice(j,i,k)-qexc) /
                  tmpa = d_one-ccover
                  zdl = d_two*(facl*sqmix-qexc) / tmpa
                  ! added correction term fac 15/03/2010
                  if ( facl*dqs < -zdl ) then
                    ! qsice(j,i,k)+qvnow
                    xlcondlim = (ccover-d_one)*facl*dqs-facl*sqmix+qxfg(iqqv)
                    chng = min(chng,xlcondlim)
                  end if
                  chng = max(chng,d_zero)
                  if ( chng < activqx ) then
                    chng = d_zero
                  end if
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
                    ngs%statscond1c(j,i,k) = ngs%statscond1c(j,i,k) + chng
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
            lactiv = qxfg(iqql) > activqx .and. ltklt0
            if ( lactiv ) then
              vpice = eeice(j,i,k) !saturation vapor pressure wrt ice
              vpliq = eeliq(j,i,k) !saturation vapor pressure wrt liq
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
              xbdd  = rwat*tk*mo2mc%phs(j,i,k)/(2.21_rkx*vpice)
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
                  (depliqrefrate+cldtopdist(j,i,k)/depliqrefdepth),d_one)
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
                ngs%statsdepos(j,i,k) = chng
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
              covptot(j,i) = d_one - ((d_one-covptot(j,i)) * &
                  (d_one - max(ccover,lccover))/(d_one-lccover))
              covptot(j,i) = max(covptot(j,i),rcovpmin)
              covpclr(j,i) = max(covptot(j,i)-ccover,d_zero)
            else
              covptot(j,i) = d_zero ! no flux - reset cover
              covpclr(j,i) = d_zero ! no flux - reset cover
            end if
            ! clear sky proportion

            !---------------------------------------------------------------
            !   WARM PHASE AUTOCONVERSION
            !---------------------------------------------------------------
            if ( ql_incld > d_zero ) then
#ifdef __PGI
              select case (iautoconv)
                case (1) ! Klein & Pincus (2000)
                  call klein_and_pincus
                case (2) ! Khairoutdinov and Kogan (2000)
                  call khairoutdinov_and_kogan
                case (3) ! Kessler(1969)
                  call kessler
                case (4) ! Sundqvist
                  call sundqvist
              end select
#else
              call selautoconv
#endif
#ifdef DEBUG
              if ( stats ) then
                if ( ltkgt0 ) then
                  ngs%statsautocvw(j,i,k) = ngs%statsautocvw(j,i,k) + rainaut
                else
                  ngs%statsautocvc(j,i,k) = ngs%statsautocvc(j,i,k) + rainaut
                end if
              end if
#endif
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
                  ngs%statsautocvc(j,i,k) = ngs%statsautocvc(j,i,k) + snowaut
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
                ! qsice(j,i,k)-qxfg(iqqv),d_zero)
                subsat = max(sqmix-qxfg(iqqv),d_zero)
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
                !     (tw1+tw2*(mo2mc%phs(j,i,k)-tw3)-tw4*(tk-tw5))
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
                        ngs%statsmelt(j,i,k) = ngs%statsmelt(j,i,k) + chng
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
            if ( chngmax > d_zero .and. qxfg(iqqr) > activqx ) then
              chng = min(qxfg(iqqr),chngmax)
              chng = max(chng,d_zero)
              qsexp(iqqs,iqqr) = qsexp(iqqs,iqqr) + chng
              qsexp(iqqr,iqqs) = qsexp(iqqr,iqqs) - chng
              qxfg(iqqs) = qxfg(iqqs) + chng
              qxfg(iqqr) = qxfg(iqqr) - chng
#ifdef DEBUG
              if ( stats ) then
                ngs%statsfrz(j,i,k) = chng
              end if
#endif
            end if

            !-------------------
            ! Freezing of liquid
            !-------------------

            chngmax = max((thomo-tk)*rldcp,d_zero)
            if ( chngmax > d_zero .and. qxfg(iqql) > activqx ) then
              chng = min(qxfg(iqql),chngmax)
              chng = max(chng,d_zero)
              qsexp(iqqi,iqql) = qsexp(iqqi,iqql) + chng
              qsexp(iqql,iqqi) = qsexp(iqql,iqqi) - chng
              qxfg(iqql) = qxfg(iqql) - chng
              qxfg(iqqi) = qxfg(iqqi) + chng
#ifdef DEBUG
              if ( stats ) then
                ngs%statsfrz(j,i,k) = ngs%statsfrz(j,i,k) + chng
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

            zrh = rprecrhmax + (d_one-rprecrhmax)*covpclr(j,i)/(d_one-ccover)
            zrh = min(max(zrh,rprecrhmax),d_one)

            ! This is a critical relative humidity that is used to limit
            ! moist environment to prevent the gridbox saturating when
            ! only part of the gridbox has evaporating precipitation
            qe = (qxfg(iqqv) - ccover*qsliq(j,i,k)) / (d_one-ccover)
            !---------------------------------------------
            ! humidity in moistest covpclr part of domain
            !---------------------------------------------
            qe = max(min(qe,qsliq(j,i,k)),d_zero)
            lactiv = covpclr(j,i) > d_zero .and. &
                     covptot(j,i) > d_zero .and. &
                     qpretot > d_zero .and.      &
                     qxfg(iqqr) > activqx .and.  &
                     qe < zrh*qsliq(j,i,k)
            if ( lactiv ) then
              ! note: units of preclr and qpretot differ
              !       qpretot is a mixing ratio (hence "q" in name)
              !       preclr is a rain flux
              preclr = qpretot*covpclr(j,i)/(covptot(j,i)*dtgdp)
              !--------------------------------------
              ! actual microphysics formula in beta
              !--------------------------------------
              ! sensitivity test showed multiply rain evap rate by 0.5
              beta1 = sqrt(mo2mc%phs(j,i,k)/pbot) / &
                           5.09e-3_rkx*preclr/covpclr(j,i)
              if ( beta1 > d_zero ) then
                beta = d_half*egrav*rpecons*(beta1)**0.5777_rkx
                denom = d_one + beta*dt*corqsliq
                dpr = covpclr(j,i) * beta * (qsliq(j,i,k)-qe)/denom*dp*regrav
                dpevap = dpr*dtgdp

                ! AMT just evaporate all rain if the rainfall is very small
                if ( qxfg(iqqr) < activqx ) dpevap = qxfg(iqqr)

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
                covptot(j,i) = covptot(j,i) - max(d_zero, &
                           (covptot(j,i)-ccover)*dpevap/qpretot)
                covptot(j,i) = max(covptot(j,i),rcovpmin)
              else
                chng = qxfg(iqqr)
              end if
              qsexp(iqqv,iqqr) = qsexp(iqqv,iqqr) + chng
              qsexp(iqqr,iqqv) = qsexp(iqqr,iqqv) - chng
              qxfg(iqqr)       = qxfg(iqqr) - chng
              qxfg(iqqv)       = qxfg(iqqv) + chng
#ifdef DEBUG
              if ( stats ) then
                ngs%statsrainev(j,i,k) = chng
              end if
#endif
            end if

            ! snow
            qe = (qxfg(iqqv) - ccover*qsice(j,i,k)) / (d_one-ccover)
            !---------------------------------------------
            ! humidity in moistest covpclr part of domain
            !---------------------------------------------
            qe = max(min(qe,qsice(j,i,k)),d_zero)
            lactiv = covpclr(j,i) > d_zero .and. &
                     covptot(j,i) > d_zero .and. &
                     qpretot > d_zero .and.      &
                     qxfg(iqqs) > activqx .and.  &
                     qe < zrh*qsice(j,i,k)
            if ( lactiv ) then
              ! note: units of preclr and qpretot differ
              !       qpretot is a mixing ratio (hence "q" in name)
              !       preclr is a rain flux
              preclr = qpretot*covpclr(j,i)/(covptot(j,i)*dtgdp)
              !--------------------------------------
              ! actual microphysics formula in beta
              !--------------------------------------
              beta1 = sqrt(mo2mc%phs(j,i,k)/pbot) / &
                           5.09e-3_rkx*preclr/covpclr(j,i)
              if ( beta1 >= d_zero ) then
                beta = d_half*egrav*rpecons*(beta1)**0.5777_rkx
                denom = d_one + beta*dt*corqsice
                dpr = covpclr(j,i) * beta * (qsice(j,i,k)-qe)/denom*dp*regrav
                dpevap = dpr*dtgdp

                ! sublimation of  snow
                ! AMT just evaporate all if snow is very small
                if ( qxfg(iqqs) < activqx ) dpevap = qxfg(iqqs)

                chng = min(dpevap,qxfg(iqqs))
                chng = max(chng,d_zero)
                !-------------------------------------------------------------
                ! reduce the total precip coverage proportional to evaporation
                !-------------------------------------------------------------
                covptot(j,i) = covptot(j,i) - &
                     max(d_zero,(covptot(j,i)-ccover)*dpevap/qpretot)
                covptot(j,i) = max(covptot(j,i),rcovpmin)
              else
                chng = qxfg(iqqs)
              end if
              qsexp(iqqv,iqqs) = qsexp(iqqv,iqqs) + chng
              qsexp(iqqs,iqqv) = qsexp(iqqs,iqqv) - chng
              qxfg(iqqs)       = qxfg(iqqs) - chng
              qxfg(iqqv)       = qxfg(iqqv) + chng
#ifdef DEBUG
              if ( stats ) then
                ngs%statssnowev(j,i,k) = chng
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
          sinksum(:) = d_zero
          lind2(:,:) = .false.
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
            ratio(n) = max(qx0(n),activqx) / &
              max(sinksum(n),max(qx0(n),activqx))
          end do
          !--------------------------------------------------------
          ! now sort ratio to find out which species run out first
          !--------------------------------------------------------
          iorder = argsort(ratio)

          !--------------------------------------------
          ! scale the sink terms, in the correct order,
          ! recalculating the scale factor each time
          !--------------------------------------------
          sinksum(:) = d_zero
          !----------------
          ! recalculate sum
          !----------------
          do n = 1 , nqx
            do jn = 1 , nqx
              jo = iorder(n)
              lind2(jo,jn) = qsexp(jo,jn) < d_zero
              sinksum(jo) = sinksum(jo) - qsexp(jo,jn)
            end do
          end do
          !---------------------------
          ! recalculate scaling factor
          !---------------------------
          do n = 1 , nqx
            jo = iorder(n)
            ratio(jo) = max(qx0(jo),activqx) / &
               max(sinksum(jo),max(qx0(jo),activqx))
          end do
          !------
          ! scale
          !------
          do n = 1 , nqx
            do jn = 1 , nqx
              jo = iorder(n)
              if ( lind2(jo,jn) ) then
                qsexp(jo,jn) = qsexp(jo,jn)*ratio(jo)
                qsexp(jn,jo) = qsexp(jn,jo)*ratio(jo)
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

          call mysolve

          !-------------------------------------------------------------------
          !  Precipitation/sedimentation fluxes to next level
          !  diagnostic precipitation fluxes
          !  It is this scaled flux that must be used for source to next layer
          !-------------------------------------------------------------------
          do n = 1 , nqx
            ! Generalized precipitation flux
            ! this will be the source for the k
            pfplsx(n,j,i,k+1) = fallsink(n)*qxn(n)*rdtgdp
            ! Calculate fluxes in and out of box for conservation of TL
            fluxq = convsrce(n) + fallsrce(n) - fallsink(n)*qxn(n)
            ! Calculate the water variables tendencies
            chng = qxn(n) - qx0(n)
            qxtendc(n,j,i,k) = qxtendc(n,j,i,k) + chng*rdt
            ! Calculate the temperature tendencies
            if ( iphase(n) == 1 ) then
              ttendc(j,i,k) = ttendc(j,i,k)+wlhvocp*(chng-fluxq)*rdt
            else if ( iphase(n) == 2 ) then
              ttendc(j,i,k) = ttendc(j,i,k)+wlhsocp*(chng-fluxq)*rdt
            end if
          end do

        end do ! jx : end of longitude loop
      end do   ! iy : end of latitude loop
    end do     ! kz : end of vertical loop

    if ( idynamic == 3 ) then
      do n = 1 , nqx
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              mc2mo%qxten(j,i,k,n) = qxtendc(n,j,i,k)
            end do
          end do
        end do
      end do
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            mc2mo%tten(j,i,k) = ttendc(j,i,k)
          end do
        end do
      end do
    else
      !
      ! Couple tendencies with pressure
      !
      do n = 1 , nqx
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              mc2mo%qxten(j,i,k,n) = qxtendc(n,j,i,k)*mo2mc%psb(j,i)
            end do
          end do
        end do
      end do
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            mc2mo%tten(j,i,k) = ttendc(j,i,k)*mo2mc%psb(j,i)
          end do
        end do
      end do
    end if
    !
    ! calculate the rain out kin constat remrat (s-1) for chemical scavenging
    ! here taken as rain water tendency divided by cloud liquid water mixing
    ! ratio
    if ( ichem == 1 )  then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
             mc2mo%remrat(j,i,k) = qxtendc(iqqr,j,i,k) / qx(iqql,j,i,k)
          end do
        end do
      end do
    end if

    !-------------------------------------
    ! Final enthalpy and total water diagnostics
    !-------------------------------------
    if ( budget_compute ) then

      ! Initialize the flux arrays
      sumh1(:,:,:)     = d_zero
      sumq1(:,:,:)     = d_zero
      errorq(:,:)    = d_zero
      errorh(:,:)    = d_zero

      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            dp = dpfs(j,i,k)
            tnew = tx(j,i,k)+dt*(ttendc(j,i,k)-tentkp(j,i,k))
            qvnew = qx(iqqv,j,i,k)+dt*(qxtendc(iqqv,j,i,k)-tenqkp(iqqv,j,i,k))
            if ( k > 1 ) then
              sumq1(j,i,k) = sumq1(j,i,k-1)
              sumh1(j,i,k) = sumh1(j,i,k-1)
            end if
            tmpl = qx(iqql,j,i,k)+dt*(qxtendc(iqql,j,i,k)-tenqkp(iqql,j,i,k))+&
                   qx(iqqr,j,i,k)+dt*(qxtendc(iqqr,j,i,k)-tenqkp(iqqr,j,i,k))
            tmpi = qx(iqqi,j,i,k)+dt*(qxtendc(iqqi,j,i,k)-tenqkp(iqqi,j,i,k))+&
                   qx(iqqs,j,i,k)+dt*(qxtendc(iqqs,j,i,k)-tenqkp(iqqs,j,i,k))
            tnew = tnew - wlhvocp*tmpl - wlhsocp*tmpi
            sumq1(j,i,k) = sumq1(j,i,k) + (tmpl + tmpi + qvnew)*dp*regrav
            sumh1(j,i,k) = sumh1(j,i,k) + dp*tnew
          end do
        end do
      end do
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            dp = dpfs(j,i,k)
            dtgdp = dt*egrav/dp
            rain = d_zero
            rainh = d_zero
            do n = 1 , nqx
              rain = rain + dt*pfplsx(n,j,i,k+1)
              if ( iphase(n) == 1 ) then
                rainh = rainh+wlhvocp*dtgdp*pfplsx(n,j,i,k+1)*dp
              else if ( iphase(n) == 2 ) then
                rainh = rainh+wlhsocp*dtgdp*pfplsx(n,j,i,k+1)*dp
              end if
            end do
            sumq1(j,i,k) = sumq1(j,i,k) + rain
            sumh1(j,i,k) = sumh1(j,i,k) - rainh
          end do
        end do
      end do
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            sumh1(j,i,k) = sumh1(j,i,k) / mo2mc%pfs(j,i,k+1)
            errorq(j,i) = errorq(j,i) + (sumq1(j,i,k)-sumq0(j,i,k))
            errorh(j,i) = errorh(j,i) + (sumh1(j,i,k)-sumh0(j,i,k))
          end do
        end do
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
    pfplsl(:,:,:) = d_zero
    pfplsn(:,:,:) = d_zero
    mc2mo%rainls(:,:,:) = d_zero

    !--------------------------------------------------------------------
    ! Copy general precip arrays back into FP arrays
    ! Add rain and liquid fluxes, ice and snow fluxes
    !--------------------------------------------------------------------

    ! Rain+liquid, snow+ice
    ! for each level k = 1 , kz, sum of the same phase elements
    do k = 1 , kzp1
      do i = ici1 , ici2
        do j = jci1 , jci2
          do n = 1 , nqx
            if ( iphase(n) == 1 ) then
              pfplsl(j,i,k) = pfplsl(j,i,k) + pfplsx(n,j,i,k)
              mc2mo%rainls(j,i,k) = pfplsl(j,i,k)
            else if ( iphase(n) == 2 ) then
              pfplsn(j,i,k) = pfplsn(j,i,k) + pfplsx(n,j,i,k)
            end if
          end do
        end do
      end do
    end do
    !--------------------------------------------------------------
    ! Convert the accumlated precipitation to appropriate units for
    ! the surface physics and the output sum up through the levels
    !--------------------------------------------------------------
    do i = ici1 , ici2
      do j = jci1 , jci2
        prainx = pfplsl(j,i,kzp1)*dt
        psnowx = pfplsn(j,i,kzp1)*dt
        if ( prainx > d_zero ) then
          mc2mo%rainnc(j,i) = mc2mo%rainnc(j,i) + prainx   !mm
          mc2mo%lsmrnc(j,i) = mc2mo%lsmrnc(j,i) + pfplsl(j,i,kzp1)
          mc2mo%trrate(j,i) = pfplsl(j,i,kzp1)
        end if
        if ( psnowx > d_zero ) then
          mc2mo%snownc(j,i) = mc2mo%snownc(j,i) + psnowx
          mc2mo%lsmrnc(j,i) = mc2mo%lsmrnc(j,i) + pfplsn(j,i,kzp1)
          mc2mo%trrate(j,i) = pfplsn(j,i,kzp1)
        end if
      end do
    end do

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif

    contains

    pure real(rkx) function edem(t,phase)
      implicit none
      real(rkx) , intent(in):: t , phase
      edem = phase * c5alvcp * (d_one/(t-c4les)**2) + &
               (d_one - phase) * c5alscp * (d_one/(t-c4ies)**2)
    end function edem

    pure real(rkx) function eldcpm(t)
      implicit none
      real(rkx) , intent(in):: t
      real(rkx) :: phase
      phase = max(min(d_one,((max(rtice,min(tzero,t))-rtice)* &
                              rtwat_rtice_r)**2),d_zero)
      eldcpm = phase*wlhvocp + (d_one-phase)*wlhsocp
    end function eldcpm

    pure real(rkx) function eewm(t,phase)
      implicit none
      real(rkx) , intent(in) :: t , phase
      real(rkx) :: eliq , eice
      eliq = c2es*exp(c3les*((t-tzero)/(t-c4les)))
      eice = c2es*exp(c3ies*((t-tzero)/(t-c4ies)))
      eewm = phase * eliq + (d_one-phase) * eice
    end function eewm

    subroutine nss_tompkins
      implicit none
      qexc = max((qxfg(iqqv)-ccover*sqmix)/(d_one-ccover),d_zero)
    end subroutine nss_tompkins

    subroutine nss_lohmann_and_karcher
      implicit none
      qexc = qxfg(iqqv)
    end subroutine nss_lohmann_and_karcher

    subroutine nss_gierens
      implicit none
      qexc = qxfg(iqqv)/totcond
    end subroutine nss_gierens

    subroutine klein_and_pincus
      implicit none
      rainaut = dt*auto_rate_klepi*(ql_incld**(2.3_rkx))
      qsimp(iqql,iqqv) = d_zero
      qsimp(iqqr,iqql) = qsimp(iqqr,iqql) + rainaut
      qsexp(iqqr,iqql) = d_zero
    end subroutine klein_and_pincus

    subroutine khairoutdinov_and_kogan
      implicit none
      rainaut = dt*auto_rate_khair*(ql_incld**(auto_expon_khair))
      qsimp(iqql,iqqv) = d_zero
      qsimp(iqqr,iqql) = qsimp(iqqr,iqql) + rainaut
    end subroutine khairoutdinov_and_kogan

    subroutine kessler
      implicit none
      rainaut = dt*auto_rate_kessl*autocrit_kessl
      qsimp(iqql,iqqv) = d_zero
      qsexp(iqqr,iqql) = qsexp(iqqr,iqql) - rainaut
      qsexp(iqql,iqqr) = qsexp(iqql,iqqr) + rainaut
      qsimp(iqqr,iqql) = qsimp(iqqr,iqql) + rainaut
    end subroutine kessler

    subroutine sundqvist
      implicit none
      real(rkx) :: precip , cfpr , arg , acrit
      real(rkx) , parameter :: spherefac = (4.0_rkx/3.0_rkx)*mathpi
      !alpha1 = min(rkconv*dt,ql_incld)
      alpha1 = rkconv*dt
      acrit = critauto
      if ( lccn ) then
        if ( ccn > 0._rkx ) then
          ! aerosol second indirect effect on autoconversion
          ! threshold, rcrit is a critical cloud radius for cloud
          ! water undergoing autoconversion
          ! ccn = number of ccn /m3
          acrit = ccn*spherefac*((rcrit*1e-6_rkx)**3)*rhoh2o
        endif
      endif
      !-----------------------------------------------------------
      ! parameters for cloud collection by rain and snow.
      ! note that with new prognostic variable it is now possible
      ! to replace this with an explicit collection
      ! parametrization to be replaced by Khairoutdinov and Kogan [2000]:
      !-----------------------------------------------------------
      if ( covptot(j,i) > d_zero ) then
        precip = (rainp+snowp)/covptot(j,i)
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
    end subroutine sundqvist

    subroutine mysolve
      implicit none
      integer(ik4) :: ii , jj , kk , ll , imax , n , nn
      real(rkx) :: aamax , dum , xsum , swap
      ! find implicit scaling information
      do n = 1 , nqx
        aamax = d_zero
        do jn = 1 , nqx
          if ( abs(qlhs(n,jn)) > aamax ) aamax = abs(qlhs(n,jn))
        end do
        if ( aamax == d_zero ) then
          do nn = 1 , nqx
            write(stderr,'(a,i2,f20.9)') 'QX0 ', nn , qx0(nn)
            do ll = 1 , nqx
              write(stderr,'(a,i2,i2,f20.9)') 'QLHS ', ll , nn , qlhs(ll,nn)
            end do
          end do
          call fatal(__FILE__,__LINE__, &
                     'System does not have a solution. Cannot solve.')
        end if
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
          dum = d_one/max(qlhs(n,n),activqx)
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
          if ( abs(xsum) > activqx ) ii = m
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
    end subroutine mysolve

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

    pure function argsort(a) result(b)
      implicit none
      real(rk8) , intent(in) :: a(:)
      integer(ik4) , dimension(size(a)) :: b
      integer :: n , i , imin , temp1
      real(rk8) :: temp2
      real(rk8) , dimension(size(a)) :: a2
      a2 = a
      n = size(a)
      do i = 1 , n
        b(i) = i
      end do
      if ( n == 1 ) return
      do i = 1 , n-1
        imin = minloc(a2(i:),1) + i - 1
        if ( imin /= i ) then
          temp2 = a2(i)
          a2(i) = a2(imin)
          a2(imin) = temp2
          temp1 = b(i)
          b(i) = b(imin)
          b(imin) = temp1
        end if
      end do
    end function argsort

  end subroutine nogtom

end module mod_micro_nogtom

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

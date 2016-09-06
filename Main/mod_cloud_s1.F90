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

module mod_cloud_s1
  use mod_realkinds
  use mod_dynparam
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
  use mod_runparams , only : dt
  use mod_runparams , only : ipptls , ichem , iaerosol , iindirect , rcrit
  use mod_runparams , only : budget_compute , nssopt , iautoconv
  use mod_runparams , only : auto_rate_khair , auto_rate_kessl , &
                             auto_rate_klepi
  use mod_runparams , only : rsemi , rkconv , rcovpmin , rpecons
  use mod_runparams , only : ktau
  use mod_runparams , only : rtsrf

#ifdef DEBUG
  use mod_runparams , only : stats
#endif

  implicit none

  private

  logical , parameter :: lmicro = .true.

  ! critical autoconversion
  real(rkx) , parameter :: rcldiff = 3.e-4_rkx
  real(rkx) , parameter :: rlcritsnow = 3.e-5_rkx

  real(rkx) , parameter :: auto_expon_khair = 1.47_rkx
  real(rkx) , parameter :: rldcp = d_one/(wlhsocp-wlhvocp)  ! Cp/Lf
  ! 1/autoconversion time scale (s)
  !real(rkx) , parameter :: rkconv = d_one/6000.0_rkx
  real(rkx) , parameter :: autocrit_kessl = 5.e-4_rkx
  ! real(rkx) , parameter :: qs0 = 1.0e-3_rkx             !g g^(-1)
  ! real(rkx) , parameter :: rcovpmin = 0.1
  real(rkx) , parameter :: rclcrit_land = 5.e-4_rkx
  real(rkx) , parameter :: rclcrit_sea = 3.e-4_rkx
  real(rkx) , parameter :: rprc1 = 3.e2_rkx  ! in Sundqvist = 300
  real(rkx) , parameter :: ramid = 0.8_rkx
  ! real(rkx) , parameter :: rlmin = 1.D-8
  ! max threshold rh for evaporation for a precip coverage of zero
  real(rkx) , parameter :: rprecrhmax = 0.7_rkx
  ! evaporation rate coefficient Numerical fit to wet bulb temperature
  ! real(rkx) , parameter :: tw1 = 1329.31_rkx
  ! real(rkx) , parameter :: tw2 = 0.0074615_rkx
  ! real(rkx) , parameter :: tw3 = 0.85e5_rkx
  ! real(rkx) , parameter :: tw4 = 40.637_rkx
  ! real(rkx) , parameter :: tw5 = 275.0_rkx
  ! real(rkx) , parameter :: rtaumel = 1.1880e4_rkx
  ! temperature homogeneous freezing
  real(rkx) , parameter :: thomo = 235.16_rkx  ! -38.00 Celsius
  ! Cloud fraction threshold that defines cloud top
  real(rkx) , parameter :: cldtopcf = 0.1_rkx
  ! Fraction of deposition rate in cloud top layer
  real(rkx) , parameter :: depliqrefrate = 0.1_rkx
  ! Depth of supercooled liquid water layer (m)
  real(rkx) , parameter :: depliqrefdepth = 500.0_rkx
  ! initial mass of ice particle
  real(rkx) , parameter :: iceinit = 1.e-12_rkx

  public :: allocate_mod_cloud_s1 , init_cloud_s1 , microphys

  real(rkx) :: oneodt                                 ! 1/dt
  integer(ik4) , pointer , dimension(:,:) :: ldmsk    ! mddom
  real(rkx) , pointer , dimension(:,:) :: psb         ! sfc
  real(rkx) , pointer , dimension(:,:) :: rainnc      ! sfc
  real(rkx) , pointer , dimension(:,:) :: lsmrnc      ! sfc
  real(rkx) , pointer , dimension(:,:) :: snownc      ! sfc
  real(rkx) , pointer , dimension(:,:,:) :: pfcc      ! from atm
  real(rkx) , pointer , dimension(:,:,:) :: phs       ! from atms
  real(rkx) , pointer , dimension(:,:,:) :: pfs       ! from atms
  real(rkx) , pointer , dimension(:,:,:) :: t         ! from atms
  real(rkx) , pointer , dimension(:,:,:) :: rho       ! from atms
  real(rkx) , pointer , dimension(:,:,:) :: pverv     ! from atms
  real(rkx) , pointer , dimension(:,:,:,:) :: qxx     ! from atms
  real(rkx) , pointer , dimension(:,:,:) :: radheatrt ! radiation heat rate
  real(rkx) , pointer , dimension(:,:,:) :: tten      ! tendency of temperature
  real(rkx) , pointer , dimension(:,:,:) :: qdetr     ! conv. detr. water
  real(rkx) , pointer , dimension(:,:,:) :: rainls    ! Rain from here
  real(rkx) , pointer , dimension(:,:,:,:) :: qxten   ! tendency of qx

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
  logical , pointer , dimension(:) ::   lind1
  logical , pointer , dimension(:,:) :: lind3

  real(rkx) , pointer , dimension(:,:,:):: sumh0 , sumq0
  real(rkx) , pointer , dimension(:,:,:) :: sumh1 , sumq1
  real(rkx) , pointer , dimension(:,:,:) :: errorq , errorh
  real(rkx) , pointer , dimension(:,:,:):: tentkeep
  real(rkx) , pointer , dimension(:,:,:,:) :: tenkeep
  ! Mass variables
  ! Microphysics
  real(rkx) , pointer , dimension(:,:,:) :: dqsatdt
  real(rkx) , pointer , dimension(:,:,:) :: pccn
  ! for sedimentation source/sink terms
  real(rkx) , pointer , dimension(:) :: fallsink
  real(rkx) , pointer , dimension(:) :: fallsrce
  ! for convection detrainment source and subsidence source/sink terms
  real(rkx) , pointer , dimension(:) :: convsrce
  real(rkx) , pointer , dimension(:) :: convsink
  ! total rain frac: fractional occurence of precipitation (%)
  ! for condensation
  ! distance from the top of the cloud
  real(rkx) , pointer , dimension(:,:,:) :: cldtopdist
  ! ice nuclei concentration
  real(rkx) , pointer , dimension(:,:,:) :: eewmt
  real(rkx) , pointer , dimension(:,:,:) :: qliq
  real(rkx) , pointer , dimension(:,:,:) :: qliqfrac
  real(rkx) , pointer , dimension(:,:,:) :: qicefrac
  real(rkx) , pointer , dimension(:,:,:) :: qlt
  ! fluxes convergence of species
  real(rkx) , pointer , dimension(:) :: ratio
  real(rkx) , pointer , dimension(:) :: sinksum
  real(rkx) , pointer , dimension(:,:,:) :: eew
  ! ice water saturation
  real(rkx) , pointer , dimension(:,:,:) :: qsice
  ! diagnostic mixed phase RH
  real(rkx) , pointer , dimension(:,:,:) :: qsmix
  ! Storage for eeliq , eeice
  real(rkx) , pointer , dimension(:,:,:) :: st_eeliq
  real(rkx) , pointer , dimension(:,:,:) :: st_eeice
  ! water saturation mixing ratio
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
  ! detrainment from tiedtke scheme
  real(rkx) , pointer , dimension(:,:,:) :: xqdetr
  ! real(rkx) , pointer , dimension(:,:,:) :: xqdetr2
  ! fall speeds of three categories
  real(rkx) , pointer , dimension(:) :: vqx
  ! n x n matrix storing the LHS of implicit solver
  real(rkx) , pointer , dimension(:,:) :: qlhs
  ! explicit sources and sinks
  real(rkx) , pointer , dimension(:,:) :: solqa
  ! implicit sources and sinks
  real(rkx) , pointer , dimension(:,:) :: solqb
  ! decoupled mixing ratios tendency
  real(rkx) , pointer , dimension(:,:,:,:) :: qxtendc
  ! j,i,n ! generalized precipitation flux
  real(rkx) , pointer , dimension(:,:,:,:) :: pfplsx
  real(rkx) , public  , pointer, dimension(:,:,:,:) :: qx
  ! Initial values
  real(rkx) , public  , pointer, dimension(:)   :: qx0
  ! new values for qxx at time+1
  real(rkx) , public  , pointer, dimension(:)   :: qxn
  ! first guess values including precip
  real(rkx) , public  , pointer, dimension(:)   :: qxfg
  ! first guess value for cloud fraction
  real(rkx) , public  , pointer, dimension(:,:,:)   :: fccfg
  ! relative humidity
  real(rkx) , public  , pointer, dimension(:,:,:)   :: relh
  ! saturation mixing ratio with respect to water
  real(rkx) , public  , pointer, dimension(:,:,:)   :: qsliq
  ! Delta pressure
  real(rkx) , public  , pointer, dimension(:,:,:)   :: dpfs

#ifdef DEBUG
  ! statistic only if stat =.true.
  real(rkx) , public  , pointer, dimension(:,:,:)   :: statssupw
  real(rkx) , public  , pointer, dimension(:,:,:)   :: statssupc
  real(rkx) , public  , pointer, dimension(:,:,:)   :: statserosw
  real(rkx) , public  , pointer, dimension(:,:,:)   :: statserosc
  real(rkx) , public  , pointer, dimension(:,:,:)   :: statsdetrw
  real(rkx) , public  , pointer, dimension(:,:,:)   :: statsdetrc
  real(rkx) , public  , pointer, dimension(:,:,:)   :: statsevapw
  real(rkx) , public  , pointer, dimension(:,:,:)   :: statsevapc
  real(rkx) , public  , pointer, dimension(:,:,:)   :: statscond1w
  real(rkx) , public  , pointer, dimension(:,:,:)   :: statscond1c
  real(rkx) , public  , pointer, dimension(:,:,:)   :: statscond2w
  real(rkx) , public  , pointer, dimension(:,:,:)   :: statscond2c
  real(rkx) , public  , pointer, dimension(:,:,:)   :: statsdepos
  real(rkx) , public  , pointer, dimension(:,:,:)   :: statsmelt
  real(rkx) , public  , pointer, dimension(:,:,:)   :: statsfrz
  real(rkx) , public  , pointer, dimension(:,:,:)   :: statsrainev
  real(rkx) , public  , pointer, dimension(:,:,:)   :: statssnowev
#endif
  integer(ik4) , pointer , dimension(:) :: indx
  real(rkx) , pointer , dimension(:) :: vv

  real(rkx) , parameter :: activqx = 1.0e-6_rkx
  real(rkx) , parameter :: clfeps = 1.0e-6_rkx
  real(rkx) , parameter :: zerocf = lowcld + 0.1_rkx
  real(rkx) , parameter :: onecf  = hicld - 0.1_rkx

  contains

  subroutine allocate_mod_cloud_s1
    implicit none
    if ( ipptls /= 2 ) return

    call getmem1d(vqx,1,nqx,'cmicro:vqx')
    call getmem1d(indx,1,nqx,'cmicro:indx')
    call getmem1d(vv,1,nqx,'cmicro:vv')
    call getmem1d(imelt,1,nqx,'cmicro:imelt')
    call getmem1d(lfall,1,nqx,'cmicro:lfall')
    call getmem1d(iphase,1,nqx,'cmicro:iphase')
    call getmem3d(cldtopdist,jci1,jci2,ici1,ici2,1,kz,'cmicro:cldtopdist')
    call getmem3d(qliqfrac,jci1,jci2,ici1,ici2,1,kz,'cmicro:qliqfrac')
    call getmem3d(qicefrac,jci1,jci2,ici1,ici2,1,kz,'cmicro:qicefrac')
    call getmem3d(eewmt,jci1,jci2,ici1,ici2,1,kz,'cmicro:eewmt')
    call getmem3d(qsmix,jci1,jci2,ici1,ici2,1,kz,'cmicro:qsmix')
    call getmem3d(qlt,jci1,jci2,ici1,ici2,1,kzp1,'cmicro:qlt')
    call getmem1d(iorder,1,nqx,'cmicro:iorder')
    call getmem3d(ttendc,jci1,jci2,ici1,ici2,1,kz,'cmicro:ttendc')
    call getmem1d(convsrce,1,nqx,'cmicro:convsrce')
    call getmem1d(convsink,1,nqx,'cmicro:convsink')
    call getmem3d(eew,jci1,jci2,ici1,ici2,1,kz,'cmicro:eew')
    call getmem3d(qsice,jci1,jci2,ici1,ici2,1,kz,'cmicro:qsice')
    call getmem4d(qx,1,nqx,jci1,jci2,ici1,ici2,1,kz,'cmicro:qx')
    call getmem3d(qsliq,jci1,jci2,ici1,ici2,1,kz,'cmicro:qsliq')
    call getmem1d(fallsink,1,nqx,'cmicro:fallsink')
    call getmem1d(fallsrce,1,nqx,'cmicro:fallsrce')
    call getmem1d(ratio,1,nqx,'cmicro:ratio')
    call getmem1d(sinksum,1,nqx,'cmicro:sinksum')
    call getmem3d(dqsatdt,jci1,jci2,ici1,ici2,1,kz,'cmicro:dqsatdt')
    call getmem3d(pfplsl,jci1,jci2,ici1,ici2,1,kzp1,'cmicro:pfplsl')
    call getmem3d(pfplsn,jci1,jci2,ici1,ici2,1,kzp1,'cmicro:pfplsn')
    call getmem3d(pfsqlf,jci1,jci2,ici1,ici2,1,kzp1,'cmicro:pfsqlf')
    call getmem3d(pfsqif,jci1,jci2,ici1,ici2,1,kzp1,'cmicro:pfsqif')
    call getmem3d(qliq,jci1,jci2,ici1,ici2,1,kzp1,'cmicro:qliq')
    call getmem1d(qxfg,1,nqx,'cmicro:qxfg')
    call getmem3d(fccfg,jci1,jci2,ici1,ici2,1,kz,'cmicro:fccfg')
    call getmem1d(lind1,1,nqx,'cmicro:lind1')
    call getmem3d(xqdetr,jci1,jci2,ici1,ici2,1,kz,'cmicro:xqdetr')
    !xqdetr after evaporation
    !call getmem3d(xqdetr2,jci1,jci2,ici1,ici2,1,kz,'cmicro:xqdetr2')
    call getmem3d(st_eeliq,jci1,jci2,ici1,ici2,1,kz,'cmicro:st_eeliq')
    call getmem3d(st_eeice,jci1,jci2,ici1,ici2,1,kz,'cmicro:st_eeice')
    call getmem3d(eeliqt,jci1,jci2,ici1,ici2,1,kz,'cmicro:eeliqt')
    call getmem4d(qxtendc,1,nqx,jci1,jci2,ici1,ici2,1,kz,'cmicro:qxtendc')
    call getmem1d(qx0,1,nqx,'cmicro:qx0')
    call getmem1d(qxn,1,nqx,'cmicro:qxn')
    call getmem2d(qlhs,1,nqx,1,nqx,'cmicro:qlhs')
    call getmem2d(solqa,1,nqx,1,nqx,'cmicro:solqa')
    call getmem2d(solqb,1,nqx,1,nqx,'cmicro:solqb')
    call getmem2d(lind3,1,nqx,1,nqx,'cmicro:lind3')
    call getmem4d(pfplsx,1,nqx,jci1,jci2,ici1,ici2,1,kzp1,'cmicro:pfplsx')
    call getmem3d(dpfs,jci1,jci2,ici1,ici2,1,kz,'cmicro:dpfs')
    if ( budget_compute ) then
      call getmem3d(sumq0,jci1,jci2,ici1,ici2,1,kz,'cmicro:sumq0')
      call getmem3d(sumh0,jci1,jci2,ici1,ici2,1,kz,'cmicro:sumh0')
      call getmem3d(sumq1,jci1,jci2,ici1,ici2,1,kz,'cmicro:sumq1')
      call getmem3d(sumh1,jci1,jci2,ici1,ici2,1,kz,'cmicro:sumh1')
      call getmem3d(errorq,jci1,jci2,ici1,ici2,1,kz,'cmicro:errorq')
      call getmem3d(errorh,jci1,jci2,ici1,ici2,1,kz,'cmicro:errorh')
      call getmem4d(tenkeep,1,nqx,jci1,jci2,ici1,ici2,1,kz,'cmicro:tenkeep')
      call getmem3d(tentkeep,jci1,jci2,ici1,ici2,1,kz,'cmicro:tentkeep')
    end if
#ifdef DEBUG
    if ( stats ) then
      call getmem3d(statssupw,jci1,jci2,ici1,ici2,1,kz,'cmicro:statssupw')
      call getmem3d(statssupc,jci1,jci2,ici1,ici2,1,kz,'cmicro:statssupc')
      call getmem3d(statsdetrw,jci1,jci2,ici1,ici2,1,kz,'cmicro:statsdetrw')
      call getmem3d(statsdetrc,jci1,jci2,ici1,ici2,1,kz,'cmicro:statsdetrc')
      call getmem3d(statserosw,jci1,jci2,ici1,ici2,1,kz,'cmicro:statserosw')
      call getmem3d(statserosc,jci1,jci2,ici1,ici2,1,kz,'cmicro:statserosc')
      call getmem3d(statsevapw,jci1,jci2,ici1,ici2,1,kz,'cmicro:statsevapw')
      call getmem3d(statsevapc,jci1,jci2,ici1,ici2,1,kz,'cmicro:statsevapc')
      call getmem3d(statscond1w,jci1,jci2,ici1,ici2,1,kz,'cmicro:statscond1w')
      call getmem3d(statscond1c,jci1,jci2,ici1,ici2,1,kz,'cmicro:statscond1c')
      call getmem3d(statscond2w,jci1,jci2,ici1,ici2,1,kz,'cmicro:statscond2w')
      call getmem3d(statscond2c,jci1,jci2,ici1,ici2,1,kz,'cmicro:statscond2c')
      call getmem3d(statsdepos,jci1,jci2,ici1,ici2,1,kz,'cmicro:statsdepos')
      call getmem3d(statsmelt,jci1,jci2,ici1,ici2,1,kz,'cmicro:statsmelt')
      call getmem3d(statsfrz,jci1,jci2,ici1,ici2,1,kz,'cmicro:statsfrz')
      call getmem3d(statsrainev,jci1,jci2,ici1,ici2,1,kz,'cmicro:statsrainev')
      call getmem3d(statssnowev,jci1,jci2,ici1,ici2,1,kz,'cmicro:statssnowev')
    end if
#endif
  end subroutine allocate_mod_cloud_s1

  subroutine init_cloud_s1
    use mod_atm_interface
    use mod_runparams , only : vfqr , vfqi , vfqs
    implicit none
    integer(ik4) :: n
    call assignpnt(atms%pb3d,phs)
    call assignpnt(atms%pf3d,pfs)
    call assignpnt(atms%tb3d,t)
    call assignpnt(atms%wpx3d,pverv)
    call assignpnt(atms%qxb3d,qxx)
    call assignpnt(atms%rhob3d,rho)
    call assignpnt(atms%rhb3d,relh)
    call assignpnt(aten%qx,qxten)
    call assignpnt(aten%t,tten)
    call assignpnt(sfs%psb,psb)
    call assignpnt(sfs%rainnc,rainnc)
    call assignpnt(sfs%snownc,snownc)
    call assignpnt(fcc,pfcc)
    call assignpnt(pptnc,lsmrnc)
    call assignpnt(heatrt,radheatrt)
    call assignpnt(q_detr,qdetr)
    call assignpnt(rain_ls,rainls)
    call assignpnt(mddom%ldmsk,ldmsk)
    if ( ichem == 1 .and. iaerosol == 1 .and. iindirect == 2 ) then
      call assignpnt(ccn,pccn)
    end if

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

  end subroutine init_cloud_s1

  subroutine microphys
    implicit none
    integer(ik4) :: i , j , k , n , m , jn , jo
    logical :: lactiv , ltkgt0 , ltklt0 , ltkgthomo , lcloud , lnocloud
    logical :: locast , lnocast , lliq , lccn
    real(rkx) :: rexplicit
    real(rkx) :: facl , faci , facw , corr , koop , gdp
    real(rkx) :: alfaw , phases , qice , zdelta , tmpl , &
                 tmpi , tnew , qe , rain , preclr , arg
    ! local real variables for autoconversion rate constants
    real(rkx) :: alpha1 ! coefficient autoconversion cold cloud
    real(rkx) :: tmpa
    real(rkx) :: xlcrit
    ! real(rkx) :: zqadj
    real(rkx) :: zrh
    real(rkx) :: beta , beta1
    ! local variables for condensation
    real(rkx) :: cond , dtdp , cdmax , rhc , zsig , &
                 acond , zdl , xlcondlim
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
    real(rkx) :: liqcld , icecld , licld
    real(rkx) :: supsat , subsat
    real(rkx) :: ldifdt
    real(rkx) :: qold , told , tcond , dqs
    real(rkx) :: chng , chngmax
    real(rkx) :: qexc
    real(rkx) :: icenuclei
    real(rkx) :: qpretot , covptot , covpclr
    real(rkx) :: qicetot , fluxq
    real(rkx) :: ldefr
    ! real(rkx) :: gdph_r
    ! constants for deposition process
    real(rkx) :: vpice , vpliq , xadd , xbdd , cvds , &
                 qice0 , qinew , infactor , snowaut
    ! constants for condensation and turbulent mixing erosion of clouds
    real(rkx) :: dpmxdt , wtot , dtdiab , dtforc , &
                 qp , qsat , cond1 , levap , leros
    real(rkx) :: qvnow , qlnow , qinow , sqmix , ccover , lccover
    real(rkx) :: tk , tc , dens , pbot , totliq , ccn
    real(rkx) :: snowp , rainp

    procedure () , pointer :: selautoconv => null()
    procedure () , pointer :: selnss => null()

#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'microphys'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    lccn = ( ichem == 1 .and. iaerosol == 1 .and. iindirect == 2 )

    select case(nssopt)
      case(0,1)
        selnss => nss_tompkins
      case(2)
        selnss => nss_lohmann_and_karcher
      case(3)
        selnss => nss_gierens
      case default
        call fatal(__FILE__,__LINE__, &
                   'NSSOPT IN CLOUD MUST BE IN RANGE 0-3')
    end select
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

    oneodt = d_one/dt

    ! Set the default 1.e-14_rkx = d_zero
    ! Define the inizial array qx
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          do n = 1 , nqx
            if ( qxx(j,i,k,n) > minqx ) then
              qx(n,j,i,k) = qxx(j,i,k,n)
            else
              qx(n,j,i,k) = minqx
            end if
          end do
        end do
      end do
    end do

    ! Detrainment : [xqdetr] = kg/(m^2*s)
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( qdetr(j,i,k) > 1.0e-8_rkx ) then
            xqdetr(j,i,k) = qdetr(j,i,k)
          else
            xqdetr(j,i,k) = d_zero
          end if
        end do
      end do
    end do

    ! Delta pressure
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          dpfs(j,i,k) = pfs(j,i,k+1)-pfs(j,i,k)
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
                        t(j,i,k)))-rtice)*rtwat_rtice_r)**2),d_zero)
        end do
      end do
    end do

    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          fccfg(j,i,k) = min(hicld,max(pfcc(j,i,k),lowcld))
        end do
      end do
    end do

    !--------------------------------------------------------------
    ! Calculate distance from cloud top
    ! defined by cloudy layer below a layer with cloud frac <0.01
    ! DZ = DP(JL)/(RHO(JL)*RG)
    !--------------------------------------------------------------
    cldtopdist(:,:,:) = d_zero
    do k = 2 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( fccfg(j,i,k-1) > cldtopcf .and. &
               fccfg(j,i,k)  <= cldtopcf ) then
            cldtopdist(j,i,k) = cldtopdist(j,i,k) + &
                                dpfs(j,i,k)/(rho(j,i,k)*egrav)
          end if
        end do
      end do
    end do

    ! Reset total precipitation variables
    pfplsx(:,:,:,:) = d_zero

    ! Compute supersaturations

    call eeliq
    call eeice

    ! Decouple tendencies
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          do n = 1 , nqx
            qxtendc(n,j,i,k) = qxten(j,i,k,n)/psb(j,i)
          end do
        end do
      end do
    end do
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          ttendc(j,i,k) = tten(j,i,k)/psb(j,i)
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
      errorq(:,:,:)    = d_zero
      errorh(:,:,:)    = d_zero
      sumh0(:,:,:)     = d_zero
      sumq0(:,:,:)     = d_zero
      sumh1(:,:,:)     = d_zero
      sumq1(:,:,:)     = d_zero
      tentkeep(:,:,:)  = d_zero
      tenkeep(:,:,:,:) = d_zero

      ! Record the tendencies
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            do n = 1 , nqx
              tenkeep(n,j,i,k) = qxtendc(n,j,i,k)
            end do
          end do
        end do
      end do
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            tentkeep(j,i,k) = ttendc(j,i,k)
          end do
        end do
      end do

      ! initialize the flux arrays
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            ! [tnew] = K
            tnew = t(j,i,k) + dt * (ttendc(j,i,k)-tentkeep(j,i,k))
            if ( k > 1 ) then
              sumq0(j,i,k) = sumq0(j,i,k-1) ! total water
              sumh0(j,i,k) = sumh0(j,i,k-1) ! liquid water temperature
            end if
            do n = 1 , nqx
              if ( iphase(n) == 1 ) then
                tmpl = qx(n,j,i,k)+dt*(qxtendc(n,j,i,k)-tenkeep(n,j,i,k))
                tmpi = d_zero
              else if ( iphase(n) == 2 ) then
                tmpi = qx(n,j,i,k)+dt*(qxtendc(n,j,i,k)-tenkeep(n,j,i,k))
                tmpl = d_zero
              end if
              sumq0(j,i,k) = sumq0(j,i,k) + &
                (qx(n,j,i,k)+(qxtendc(n,j,i,k)-tenkeep(n,j,i,k))*dt)* &
                dpfs(j,i,k)*regrav
              tnew = tnew - wlhvocp*tmpl - wlhsocp*tmpi
              sumq0(j,i,k) = sumq0(j,i,k) + &
                (tmpl+tmpi)*dpfs(j,i,k)*regrav    !(kg/m^2)
            end do
            ! Detrained water treated here
            qe = xqdetr(j,i,k)*dt*egrav/dpfs(j,i,k) ! 1 ?
            if ( qe > activqx ) then
              sumq0(j,i,k) = sumq0(j,i,k)+xqdetr(j,i,k)*dt
              alfaw = qliq(j,i,k)
              tnew = tnew-(wlhvocp*alfaw+wlhsocp*(d_one-alfaw))*qe
            end if
            sumh0(j,i,k) = sumh0(j,i,k) + dpfs(j,i,k)*tnew
          end do
        end do
      end do
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            sumh0(j,i,k) = sumh0(j,i,k)/(pfs(j,i,k+1)-pfs(j,i,1))
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
          !---------------------------------------------
          ! mixed phase saturation
          !--------------------------------------------
          phases = qliq(j,i,k)
          eewmt(j,i,k) = min(eewm(t(j,i,k),phases)/phs(j,i,k),d_half)
          qsmix(j,i,k) = eewmt(j,i,k)
          ! ep1 = rwat/rgas - d_one
          qsmix(j,i,k) = qsmix(j,i,k)/(d_one-ep1*qsmix(j,i,k))
          !--------------------------------------------
          ! ice saturation T < 273K
          ! liquid water saturation for T > 273K
          !--------------------------------------------
          zdelta = delta(t(j,i,k))
          eew(j,i,k) = min((zdelta*st_eeliq(j,i,k) + &
               (d_one-zdelta)*st_eeice(j,i,k))/phs(j,i,k),d_half)
          eew(j,i,k) = min(d_half,eew(j,i,k))
          !qsi saturation mixing ratio with respect to ice
          !qsice(j,i,k) = eew(j,i,k)/(d_one-ep1*eew(j,i,k))
          !ice water saturation
          qsice(j,i,k) = min(st_eeice(j,i,k)/phs(j,i,k),d_half)
          qsice(j,i,k) = qsice(j,i,k)/(d_one-ep1*qsice(j,i,k))
          !----------------------------------
          ! liquid water saturation
          !----------------------------------
          !eeliq is the saturation vapor pressure es(T)
          !the saturation mixing ratio is ws = es(T)/p *0.622
          !ws = ws/(-(d_one/eps - d_one)*ws)
          eeliqt(j,i,k) = min(st_eeliq(j,i,k)/phs(j,i,k),d_half)
          qsliq(j,i,k) = eeliqt(j,i,k)
          qsliq(j,i,k) = qsliq(j,i,k)/(d_one-ep1*qsliq(j,i,k))
        end do
      end do
    end do

    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          !-------------------------------------------------------------------
          ! Calculate liq/ice fractions (no longer a diagnostic relationship)
          !-------------------------------------------------------------------
          qlt(j,i,k) = qx(iqql,j,i,k)+qx(iqqi,j,i,k)
          if ( qlt(j,i,k) > activqx ) then
            qliqfrac(j,i,k) = qx(iqql,j,i,k)/qlt(j,i,k)
            qicefrac(j,i,k) = d_one-qliqfrac(j,i,k)
          else
            qliqfrac(j,i,k) = d_zero
            qicefrac(j,i,k) = d_zero
          end if
        end do
      end do
    end do

#ifdef DEBUG
    if ( stats ) then
      statssupw(:,:,:) = d_zero
      statssupc(:,:,:) = d_zero
    end if
#endif

    do i = ici1 , ici2
      do j = jci1 , jci2

        do n = 1 , nqx
          qxfg(n) = qx(n,j,i,1)
          qx0(n) = qxfg(n)
        end do

        qvnow = qx0(iqqv)
        qlnow = qx0(iqql)
        qinow = qx0(iqqi)

        solqa(:,:)  = d_zero
        solqb(:,:)  = d_zero
        supsat      = d_zero
        subsat      = d_zero
        convsrce(:) = d_zero
        fallsrce(:) = d_zero
        fallsink(:) = d_zero
        convsink(:) = d_zero
        qpretot     = d_zero
        covptot     = d_zero
        covpclr     = d_zero
        qicetot     = d_zero
        ldefr       = d_zero

        !---------------------------------
        ! First guess microphysics
        !---------------------------------

        pbot  = pfs(j,i,kzp1)
        dp      = dpfs(j,i,1)
        tk      = t(j,i,1)
        tc      = tk - tzero
        dens    = rho(j,i,1)
        alfaw   = qliq(j,i,1)
        sqmix   = qsmix(j,i,1)
        ccover  = fccfg(j,i,1)
        lccover = d_zero
        totliq  = qlt(j,i,1)
        if ( lccn) ccn = pccn(j,i,1)
        rainp = d_zero
        snowp = d_zero

        ltkgt0    = ( tk > tzero )
        ltklt0    = ( .not. ltkgt0 )
        ltkgthomo = ( tk > thomo )
        lcloud    = ( ccover >= zerocf )
        lnocloud  = ( .not. lcloud )
        locast    = ( ccover >= onecf )
        lnocast   = ( .not. locast )
        lliq      = ( totliq > activqx )

        ! Derived variables needed
        gdp = egrav/dp        ! g/dp  =(1/m)
        dtgdp = dt*gdp        ! (dt*g)/dp =(dt/m)
        rdtgdp = dp*(d_one/(dt*egrav)) ! dp/(gdt)=m/dt  [Kg/m2/s]
        !------------------------------------
        ! calculate dqs/dT
        !------------------------------------
        ! liquid
        facw     = c5les/((tk - c4les)**2)
        corr     = d_one/(d_one - ep1*eeliqt(j,i,1))
        dqsliqdt = facw*corr*qsliq(j,i,1)
        corqsliq = d_one + wlhvocp*dqsliqdt
        ! ice
        faci     = c5ies/((tk - c4ies)**2)
        corr     = d_one/(d_one - ep1*eew(j,i,1))
        dqsicedt = faci*corr*qsice(j,i,1)
        corqsice = d_one + wlhsocp*dqsicedt
        ! diagnostic mixed
        facl     = alfaw*facw + (d_one - alfaw)*faci
        corr     = d_one/(d_one - ep1*eewmt(j,i,1))
        dqsmixdt = facl*corr*sqmix
        corqsmix = d_one/(d_one + eldcpm(tk)*dqsmixdt)
        !--------------------------------
        ! evaporation/sublimation limits
        !--------------------------------
        evaplimmix = max((sqmix-qvnow)*corqsmix,d_zero)
        !--------------------------------
        ! in-cloud consensate amount
        !--------------------------------
        tmpa = d_one/ccover
        liqcld = qlnow*tmpa
        icecld = qinow*tmpa
        licld  = liqcld + icecld

        !------------------------------------------------
        ! Evaporate very small amounts of liquid and ice
        !------------------------------------------------
        if ( qlnow < activqx ) then
          solqa(iqqv,iqql) =  qlnow
          solqa(iqql,iqqv) = -qlnow
          qxfg(iqql) = qxfg(iqql) - qlnow
        end if
        if ( qinow < activqx ) then
          solqa(iqqv,iqqi) =  qinow
          solqa(iqqi,iqqv) = -qinow
          qxfg(iqqi) = qxfg(iqqi) - qinow
        end if

        !------------------------------------------------------------------
        !  ICE SUPERSATURATION ADJUSTMENT
        !------------------------------------------------------------------
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
        !-------------------------------------------------------------------
        !-----------------------------------
        ! Supersaturation limit (from Koop)
        !-----------------------------------
        if ( nssopt == 0 )  then
          facl = d_one
        else
          if ( ltkgt0 ) then
            facl = d_one
          else
            koop = fkoop(tk,st_eeliq(j,i,1),st_eeice(j,i,1))
            facl = ccover + koop*(d_one-ccover)
          end if
        end if

        !-------------------------------------------------------------------
        ! Calculate supersaturation wrt Koop including dqs/dT
        ! correction factor
        !-------------------------------------------------------------------
        ! Here the supersaturation is turned into liquid water
        ! However, if the temperature is below the threshold for homogeneous
        ! freezing then the supersaturation is turned instantly to ice.
        ! Moreover the RH is clipped to the limit of
        ! qv_max = qs * (fcc + (1-fcc) *RH_homo )
        !--------------------------------------------------------------------
        ! qsmix?
        ! corqslsliq? shouldn't it be corqsmix?
        supsat = max((qvnow-facl*sqmix)/corqsliq,d_zero)
        ! e < esi, because for e > esi ice still present
        subsat = min((qvnow-facl*sqmix)/corqsliq,d_zero)

        if ( supsat > dlowval ) then
          if ( ltkgthomo ) then
            ! turn supersaturation into liquid water
            solqa(iqql,iqqv) = solqa(iqql,iqqv)+supsat
            solqa(iqqv,iqql) = solqa(iqqv,iqql)-supsat
            qxfg(iqql) = qxfg(iqql)+supsat
          else
            ! turn supersaturation into ice water
            solqa(iqqi,iqqv) = solqa(iqqi,iqqv)+supsat
            solqa(iqqv,iqqi) = solqa(iqqv,iqqi)-supsat
            qxfg(iqqi) = qxfg(iqqi)+supsat
          end if
        else
          if ( subsat < d_zero .and. lnocloud .and. lliq ) then
            ! turn subsaturation into vapor, where there is no cloud
            excess = totliq+subsat
            if ( excess < d_zero ) then
              if ( ltkgthomo ) then
                evapl = max(-totliq,-evaplimmix)!*oneodt
                solqa(iqqv,iqql) = solqa(iqqv,iqql)-evapl
                solqa(iqql,iqqv) = solqa(iqql,iqqv)+evapl
                qxfg(iqql) = qxfg(iqql)+evapl
              else
                evapi = max(-totliq,-evaplimmix)!*oneodt
                ! turn subsaturation into vapour
                solqa(iqqv,iqqi) = solqa(iqqv,iqqi)-evapi
                solqa(iqqi,iqqv) = solqa(iqqi,iqqv)+evapi
                qxfg(iqqi) = qxfg(iqqi)+evapi
              end if
            end if
          end if
        end if

#ifdef DEBUG
        if ( stats ) then
          if ( supsat > dlowval ) then
            if ( ltkgthomo ) then
              statssupw(j,i,1) = supsat
            else
              statssupc(j,i,1) = supsat
            end if
          else
            if ( subsat < d_zero .and. lnocloud .and. lliq ) then
              excess = totliq+subsat
              if ( excess < d_zero ) then
                if ( ltkgthomo ) then
                  evapl = max(-totliq,-evaplimmix)!*oneodt
                  statssupw(j,i,1) =  statssupw(j,i,1) + evapl
                else
                  evapi = max(-totliq,-evaplimmix)!*oneodt
                  statssupc(j,i,1) =  statssupc(j,i,1) - evapi
                end if
              end if
            end if
          end if
        end if
#endif
        !
        ! call addpath(iqql,iqqv,supsatl,solqa,solqb,d_zero,qxfg)
        ! call addpath(iqqi,iqqv,supsati,solqa,solqb,d_zero,qxfg)
        !
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
        ! SOLQA/B:q(IQa,IQb)
        !
        ! Thus if SOLQA/B(IQL,IQV) = K where K > 0 then this is
        ! a source of IQL and a sink of IQV
        !
        ! put 'magic' source terms such as LUDE from
        ! detrainment into explicit source/sink array diagnognal
        ! SOLQA(IQL,IQL)=LUDE
        !--------------------------------------------------------
        ! Define the microphysics
        ! the matrix will be sparse is this a problem ?
        ! (X,Y) means a sink of X and a source of Y
        ! for the implementation I will use flexible pointers
        ! such that it will be written (IQR,IQG) to indicate graupel to rain
        ! and the parametrization can have different variables switched on
        ! and off.
        ! each of these is a parametrization for a microphysical process.
        !------------------------------------------------------------------
        !                 DETRAINMENT FROM CONVECTION
        !------------------------------------------------------------------
        !chng = d_zero
        xqdetr(j,i,1) = xqdetr(j,i,1)*dtgdp  !kg/kg
        if ( xqdetr(j,i,1) > activqx ) then
          !qice = 1 if T < 250, qice = 0 if T > 273
          qice   = d_one-alfaw
          convsrce(iqql) = alfaw*xqdetr(j,i,1)
          convsrce(iqqi) = qice*xqdetr(j,i,1)
          solqa(iqql,iqql) = solqa(iqql,iqql) + convsrce(iqql)
          solqa(iqqi,iqqi) = solqa(iqqi,iqqi) + convsrce(iqqi)
          qxfg(iqql) = qxfg(iqql)+convsrce(iqql)
          qxfg(iqqi) = qxfg(iqqi)+convsrce(iqqi)
        end if

#ifdef DEBUG
        if ( stats ) then
          statsdetrw(j,i,1) = convsrce(iqql)
          statsdetrc(j,i,1) = convsrce(iqqi)
        end if
#endif

        !------------------------------------------------------------------
        ! Turn on/off microphysics
        !------------------------------------------------------------------

        if ( lmicro ) then

          !---------------------------------------
          ! EROSION OF CLOUDS BY TURBULENT MIXING
          !--------------------------------------
          ! rcldiff  : Diffusion coefficient for evaporation by turbulent
          ! mixing (IBID., EQU. 30) rcldiff = 3.0e-6_rkx
          ldifdt = rcldiff*dt
          if ( xqdetr(j,i,1) > d_zero ) ldifdt = 5.0_rkx*ldifdt
          !Increase by factor of 5 for convective points
          if ( lliq ) then
            leros = ccover * ldifdt * max(sqmix-qvnow,d_zero)
            leros = min(leros,evaplimmix)
            leros = min(leros,totliq)
            facl = qliqfrac(j,i,1)*leros
            faci = qicefrac(j,i,1)*leros
            solqa(iqql,iqqv) = solqa(iqql,iqqv) - facl
            solqa(iqqv,iqql) = solqa(iqqv,iqql) + facl
            solqa(iqqi,iqqv) = solqa(iqqi,iqqv) - faci
            solqa(iqqv,iqqi) = solqa(iqqv,iqqi) + faci
            qxfg(iqql) = qxfg(iqql) - facl
            qxfg(iqqi) = qxfg(iqqi) - faci
          end if

#ifdef DEBUG
          if ( stats ) then
            if ( lliq ) then
              statserosw(j,i,1) = qliqfrac(j,i,1)*leros
              statserosc(j,i,1) = qicefrac(j,i,1)*leros
            end if
          end if
#endif
          !------------------------------------------------------------------
          ! CONDENSATION/EVAPORATION DUE TO DQSAT/DT
          !------------------------------------------------------------------
          ! calculate dqs/dt
          ! Note: For the separate prognostic Qi and Ql, one would ideally use
          ! Qsat/DT wrt liquid/Koop here, since the physics is that new clouds
          ! forms by liquid droplets [liq] or when aqueous aerosols [Koop]
          ! form.
          ! These would then instant. freeze if T<-38C or lead to ice growth
          ! by deposition in warmer mixed phase clouds.  However, since we do
          ! not have a separate prognostic equation for in-cloud humidity or a
          ! statistical scheme approach in place, the depositional growth of
          ! ice in the mixed phase can not be modelled and we resort to
          ! supersaturation
          ! wrt ice instanteously converting to ice over one timestep
          ! (see Tompkins et al. QJRMS 2007 for details)
          ! Thus for the initial implementation the diagnostic mixed phase is
          ! retained for the moment, and the level of approximation noted.
          !------------------------------------------------------------------
          dtdp   = rovcp*tk/phs(j,i,1)
          dpmxdt = dp*oneodt
          wtot   = pverv(j,i,1)
          wtot   = min(dpmxdt,max(-dpmxdt,wtot))
          dtdiab = min(dpmxdt*dtdp,max(-dpmxdt*dtdp,radheatrt(j,i,1)))*dt + &
                        wlhfocp*ldefr
          ! ldefr = 0
          ! note: ldefr should be set to the difference between the mixed
          ! phase functions in the convection and cloud scheme, and
          ! for now we set it to zero and the functions are the same.
          ! In RegCM not all convection schemes provide such info.
          dtforc = dtdp*wtot*dt + dtdiab
          qold   = sqmix
          told   = tk
          tcond  = tk+dtforc
          tcond  = max(tcond,160.0_rkx)
          ! the goal is to produce dqs = qsmix - qold, where qsmix is
          ! reduced because of the condensation. so that dqs is negative?
          qp = d_one/phs(j,i,1)
          phases = max(min(d_one,((max(rtice,min(tzero, &
                   tcond))-rtice)*rtwat_rtice_r)**2),d_zero)
          ! saturation mixing ratio ws
          qsat = eewm(tcond,phases)*qp
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
          tcond = tcond+eldcpm(tcond)*cond1
          sqmix = sqmix - cond1
          dqs = sqmix - qold
          sqmix = qold
          tcond = told

          !----------------------------------------------------------------
          ! dqs > 0:  evaporation of clouds
          !----------------------------------------------------------------
          ! erosion term is linear in l
          ! changed to be uniform distribution in cloud region
          ! previous function based on delta distribution in cloud:
          if ( dqs > d_zero ) then
            !levap = C*min( dqs/dt , (qi+ql)/C )
            levap = ccover*min(dqs,licld)
            levap = min(levap,evaplimmix)
            levap = min(levap,max(sqmix-qvnow,d_zero))
            facl = qliqfrac(j,i,1)*levap
            faci = qicefrac(j,i,1)*levap
            solqa(iqqv,iqql) = solqa(iqqv,iqql) + facl
            solqa(iqql,iqqv) = solqa(iqql,iqqv) - facl
            solqa(iqqv,iqqi) = solqa(iqqv,iqqi) + faci
            solqa(iqqi,iqqv) = solqa(iqqi,iqqv) - faci
            qxfg(iqql) = qxfg(iqql) - facl
            qxfg(iqqi) = qxfg(iqqi) - faci
          end if

#ifdef DEBUG
          if ( stats ) then
            if ( dqs > d_zero ) then
              statsevapw(j,i,1) = qliqfrac(j,i,1)*levap
              statsevapc(j,i,1) = qicefrac(j,i,1)*levap
            end if
          end if
#endif

          !-----------------------------------------------------------------
          ! dqs < 0: formation of clouds
          !-----------------------------------------------------------------
          ! (1) increase of cloud water in existing clouds
          chng = d_zero

          if ( lcloud .and. dqs <= -activqx ) then
            ! new limiter
            chng = max(-dqs,d_zero)
            ! old limiter
            !  (significantly improves upper tropospheric humidity rms)
            phases = alfaw
            if ( locast ) then
              corr = d_one/(d_one-ep1*sqmix)
              cdmax = (qvnow-sqmix) / (d_one + corr*sqmix*edem(tk,phases))
            else
              cdmax = (qvnow-ccover*sqmix) / ccover
            end if
            chng = max(min(chng,cdmax),d_zero)
            ! end old limiter
            chng = ccover*chng
            if ( chng < activqx ) chng = d_zero
            !----------------------------------------------------------------
            ! All increase goes into liquid unless so cold cloud
            ! homogeneously freezes
            ! include new liquid formation in first guess value, otherwise
            ! liquid remains at cold temperatures until next timestep.
            !----------------------------------------------------------------
            if ( ltkgthomo ) then
              solqa(iqql,iqqv) = solqa(iqql,iqqv) + chng
              solqa(iqqv,iqql) = solqa(iqqv,iqql) - chng
              qxfg(iqql) = qxfg(iqql) + chng
            else
              solqa(iqqi,iqqv) = solqa(iqqi,iqqv) + chng
              solqa(iqqv,iqqi) = solqa(iqqv,iqqi) - chng
              qxfg(iqqi) = qxfg(iqqi) + chng
            end if
          end if

#ifdef DEBUG
          if ( stats ) then
            if ( ltkgthomo ) then
              statscond1w(j,i,1) = chng
            else
              statscond1c(j,i,1) = chng
            end if
          end if
#endif

          chng = d_zero
          qexc = d_zero

          ! (2) generation of new clouds (dc/dt>0)
          if ( dqs <= -activqx .and. lnocast ) then

            call selnss

            !---------------------------
            ! critical relative humidity
            !---------------------------
            ! *RAMID*   REAL    BASE VALUE FOR CALCULATION OF RELATIVE
            !                   HUMIDITY THRESHOLD FOR ONSET OF STRATIFORM
            !                   CONDENSATION (TIEDTKE, 1993, EQUATION 24)
            rhc = ramid !=0.8
            zsig = phs(j,i,1)/pbot
            ! increase RHcrit to 1.0 towards the surface (sigma>0.8)
            if ( zsig > ramid ) then
              rhc = ramid + (d_one-ramid)*((zsig-ramid)/0.2_rkx)**2
            end if
            !---------------------------
            ! supersaturation options
            !---------------------------
            if ( ltkgt0 .or. nssopt == 0 ) then
              ! no ice supersaturation allowed
              facl = d_one
            else
              ! ice supersaturation
              facl = fkoop(tk,st_eeliq(j,i,1),st_eeice(j,i,1))
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
              ! qsice(j,i,1)-qexc) /
              tmpa = d_one-ccover
              zdl = d_two*(facl*sqmix-qexc) / tmpa
              ! added correction term fac 15/03/2010
              if ( facl*dqs < -zdl ) then
                ! qsice(j,i,1)+qvnow
                xlcondlim = (ccover-d_one)*facl*dqs - &
                             facl*sqmix+qvnow
                chng = min(chng,xlcondlim)
              end if
              chng = max(chng,d_zero)
              if ( chng < activqx ) then
                chng = d_zero
                acond = d_zero
              end if
              !-------------------------------------------------------------
              ! all increase goes into liquid unless so cold cloud
              ! homogeneously freezes
              ! include new liquid formation in first guess value, otherwise
              ! liquid remains at cold temperatures until next timestep.
              !-------------------------------------------------------------
              if ( ltkgthomo ) then
                solqa(iqql,iqqv) = solqa(iqql,iqqv) + chng
                solqa(iqqv,iqql) = solqa(iqqv,iqql) - chng
                qxfg(iqql) = qxfg(iqql) + chng
              else ! homogeneous freezing
                solqa(iqqi,iqqv) = solqa(iqqi,iqqv) + chng
                solqa(iqqv,iqqi) = solqa(iqqv,iqqi) - chng
                qxfg(iqqi) = qxfg(iqqi) + chng
              end if
            end if
          end if

#ifdef DEBUG
          if ( stats ) then
            if ( ltkgthomo ) then
              statscond2w(j,i,1) = chng
            else
              statscond2c(j,i,1) = chng
            end if
          end if
#endif

          !------------------------------------------------------------------
          ! DEPOSITION: Growth of ice by vapour deposition
          !------------------------------------------------------------------
          ! Following Rotstayn et al. 2001:
          ! does not use the ice nuclei number from cloudaer.F90
          ! but rather a simple Meyers et al. 1992 form based on the
          ! supersaturation and assuming clouds are saturated with
          ! respect to liquid water (well mixed), (or Koop adjustment)
          ! Growth considered as sink of liquid water if present so
          ! Bergeron-Findeisen adjustment in autoconversion term no
          ! longer needed
          !-----------------------------------------------------------------

          chng = d_zero

          !--------------------------------------------------------------
          ! only treat depositional growth if liquid present. due to fact
          ! that can not model ice growth from vapour without additional
          ! in-cloud water vapour variable
          ! If water droplets are present, the ice crystals are in an
          ! environment supersaturated with respect to ice and grow by
          ! deposition, reducing the water vapour and leading to
          ! subsaturation with respect to water.
          ! Water droplets then evaporate and the process continues with
          ! ice growth until the water droplets are completely evaporated.
          ! Thus in mixed phase clouds, the deposition process acts as a
          ! sink of cloud liquid and a source of ice cloud.
          !--------------------------------------------------------------
          if ( ltklt0 .and. qxfg(iqql) > activqx ) then
            vpice = st_eeice(j,i,1) !saturation vapor pressure wrt ice
            vpliq = st_eeliq(j,i,1) !saturation vapor pressure wrt liq
            ! Meyers et al 1992
            icenuclei = d_1000*exp(12.96_rkx*((vpliq-vpice)/vpice)-0.639_rkx)
            !------------------------------------------------
            !   2.4e-2 is conductivity of air
            !   8.87 = 700**1/3 = density of ice to the third
            !------------------------------------------------
            xadd  = wlhs*(wlhs/(rwat*tk)-d_one)/(2.4e-2_rkx*tk)
            xbdd  = rwat*tk*phs(j,i,1)/(2.21_rkx*vpice)
            cvds = 7.8_rkx*(icenuclei/dens)** &
                   0.666_rkx*(vpliq-vpice)/(8.87_rkx*(xadd+xbdd)*vpice)
            !-----------------------------------------------------
            ! iceinit = 1.e-12 is initial mass of ice particle
            !-----------------------------------------------------
            qice0 = max(icecld, icenuclei*iceinit/dens)
            !------------------
            ! new value of ice condensate amount ( Rotstayn et al. (2000) )
            !------------------
            qice0 = max(qice0,d_zero)
            cvds = max(cvds,d_zero)
            qinew = (0.666_rkx*cvds*dt+qice0**0.666_rkx)**1.5_rkx
            !---------------------------
            ! grid-mean deposition rate:
            !---------------------------
            chng = max(ccover*(qinew-qice0),d_zero)*2.0_rkx
            ! above increased by factor of 2 to retain similar mixed
            ! phase liq as in diagnostic scheme
            !---------------------------------------------------------------
            ! limit deposition to liquid water amount
            ! if liquid is all frozen, ice would use up reservoir of water
            ! vapour in excess of ice saturation mixing ratio - however this
            ! can not be represented without a in-cloud humidity variable.
            ! using the grid-mean humidity would imply a large artificial
            ! horizontal flux from the clear sky to the cloudy area.
            ! we thus rely on the supersaturation check to clean up any
            ! remaining supersaturation
            !---------------------------------------------------------------
            ! limit to liquid water amount
            chng = min(chng,qxfg(iqql))
            !---------------------------------------------------------------
            ! at top of cloud, reduce deposition rate near cloud top to
            ! account for small scale turbulent processes, limited ice
            ! nucleation and ice fallout
            !---------------------------------------------------------------
            ! Fraction of deposition rate in cloud top layer
            ! depliqrefrate = 0.1_rkx
            ! Depth of supercooled liquid water layer (m)
            ! depliqrefdepth = 500.0_rkx
            infactor = min(icenuclei/15000.0_rkx, d_one)
            chng = chng*min(infactor + (d_one-infactor)* &
                   (depliqrefrate+cldtopdist(j,i,1)/depliqrefdepth),d_one)
            !--------------
            ! add to matrix
            !--------------
            solqa(iqqi,iqql) = solqa(iqqi,iqql) + chng
            solqa(iqql,iqqi) = solqa(iqql,iqqi) - chng
            qxfg(iqql) = qxfg(iqql) - chng
            qxfg(iqqi) = qxfg(iqqi) + chng
          end if

#ifdef DEBUG
          if ( stats ) then
            statsdepos(j,i,1) = chng
          end if
#endif

          tmpa = d_one/ccover
          liqcld = qxfg(iqql)*tmpa
          icecld = qxfg(iqqi)*tmpa

          !------------------------------------------------------------------
          !  SEDIMENTATION/FALLING OF *ALL* MICROPHYSICAL SPECIES
          !     now that rain and snow species are prognostic
          !     the precipitation flux can be defined directly level by level
          !     There is no vertical memory required from the flux variable
          !------------------------------------------------------------------
          do n = 1 , nqx
            if ( lfall(n) ) then
              ! Sink to next layer, constant fall speed
              fallsink(n) = dtgdp*vqx(n)*dens   !Kg/Kg
            end if  !lfall
          end do ! n

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
          if ( qpretot > dlowval ) then
            covptot = d_one - ((d_one-covptot) * &
                (d_one - max(ccover,lccover))/(d_one - min(lccover,hicld)))
            covptot = max(covptot,rcovpmin)   !rcovpmin = 0.1
            ! clear sky proportion
            covpclr = max(d_zero,covptot-ccover)
          else
            covptot = d_zero  ! no flux - reset cover
            covpclr = d_zero  ! reset clear sky proportion
          end if
          !---------------------------------------------------------------
          !                         AUTOCONVERSION
          !---------------------------------------------------------------
          ! Warm clouds
          if ( liqcld > activqx ) then
            call selautoconv
          end if

          ! Cold clouds
          ! Snow Autoconversion rate follow Lin et al. 1983
          if ( ltklt0 ) then
            if ( icecld > activqx ) then
              alpha1 = dt*1.0e-3_rkx*exp(0.025*tc)
              xlcrit = rlcritsnow
              arg = (icecld/xlcrit)**2
              if ( arg < 25.0_rkx ) then
                snowaut = alpha1 * (d_one - exp(-arg))
              else
                snowaut = alpha1
              end if
              solqb(iqqs,iqqi) = solqb(iqqs,iqqi) + snowaut
            end if
          end if

          !---------------------------------------------------------------
          !                         MELTING
          !---------------------------------------------------------------
          ! The melting of ice and snow are treated explicitly.
          ! First water and ice saturation are found
          !---------------------------------------------
          ! ice saturation T < 273K
          ! liquid water saturation for T > 273K
          !---------------------------------------------

          do n = 1 , nqx
            if ( iphase(n) == 2 ) then
              qicetot = qicetot + qxfg(n)
            end if
          end do

          chngmax = d_zero

          if ( qicetot > activqx .and. ltkgt0 ) then
            ! Calculate subsaturation
            ! qsice(j,i,1)-qvnow,d_zero)
            subsat = max(sqmix-qvnow,d_zero)
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
            tdiff = tc ! - subsat * &
            !    (tw1+tw2*(phs(j,i,1)-tw3)-tw4*(tk-tw5))
            ! Ensure CONS1 is positive so that MELTMAX = 0 if TDMTW0 < 0
            cons1 = d_one ! abs(dt*(d_one + d_half*tdiff)/rtaumel)
            chngmax = max(tdiff*cons1*rldcp,d_zero)
          end if

          ! Loop over frozen hydrometeors (iphase == 2 (ice, snow))

          chng = d_zero
          if ( chngmax > dlowval .and. qicetot > activqx ) then
            do n = 1, nqx
              if ( iphase(n) == 2 ) then
                m = imelt(n) ! imelt(iqqi)=iqql, imelt(iqqs)=iqqr
                if ( m < 0 ) cycle
                phases = qxfg(n)/qicetot
                chng = min(qxfg(n),phases*chngmax)
                chng = max(chng,d_zero)
                ! n = iqqi,iqqs; m = iqql,iqqr
                qxfg(n) =  qxfg(n) - chng
                qxfg(m) =  qxfg(m) + chng
                solqa(m,n) =  solqa(m,n) + chng
                solqa(n,m) =  solqa(n,m) - chng
              end if
            end do
          end if

#ifdef DEBUG
          if ( stats ) then
            statsmelt(j,i,1) = chng
          end if
#endif
          !------------------------------------------------------------!
          !                         FREEZING                           !
          !------------------------------------------------------------!

          ! Freezing of rain.
          ! All rain freezes in a timestep if the temperature is below 0 C
          ! calculate sublimation latent heat

          chngmax = max((tzero-tk)*rldcp,d_zero)
          chng = d_zero

          if ( chngmax > dlowval .and. qxfg(iqqr) > activqx ) then
            chng = min(qxfg(iqqr),chngmax)
            chng = max(chng,d_zero)
            solqa(iqqs,iqqr) = solqa(iqqs,iqqr) + chng
            solqa(iqqr,iqqs) = solqa(iqqr,iqqs) - chng
          end if

#ifdef DEBUG
          if ( stats ) then
            statsfrz(j,i,1) = chng
          end if
#endif

          ! Freezing of liquid

          chngmax = max((thomo-tk)*rldcp,d_zero)
          chng = d_zero

          if ( chngmax > dlowval .and. qxfg(iqql) > activqx ) then
            chng = min(qxfg(iqql),chngmax)
            chng = max(chng,d_zero)
            solqa(iqqi,iqql) = solqa(iqqi,iqql) + chng
            solqa(iqql,iqqi) = solqa(iqql,iqqi) - chng
            qxfg(iqql) = qxfg(iqql) - chng
            qxfg(iqqi) = qxfg(iqqi) + chng
          end if

#ifdef DEBUG
          if ( stats ) then
            statsfrz(j,i,1) = statsfrz(j,i,1) + chng
          end if
#endif

          !---------------------------------------------------------------
          !                         EVAPORATION
          !---------------------------------------------------------------
          !------------------------------------------------------------
          ! recalculate qpretot since melting term may have changed it
          ! rprecrhmax = 0.7 is the threshold for the clear-sky RH that
          ! can be reached by evaporation of precipitation. THis assumption
          ! is done to prevent the gridbox saturating due to the evaporation
          ! of precipitation occuring in a portion of the grid
          !------------------------------------------------------------

          ! if ( iphase(n) == 2 ) then
          ! qpretot = qpretot+qxfg(n)
          qpretot = qxfg(iqqs) + qxfg(iqqr)
          ! end if

          ! rain
          chng = d_zero
          zrh = rprecrhmax + (d_one-rprecrhmax)*covpclr/(d_one-ccover)
          zrh = min(max(zrh,rprecrhmax),d_one)
          ! This is a critical relative humidity that is used to limit
          ! moist environment to prevent the gridbox saturating when
          ! only part of the gridbox has evaporating precipitation
          qe = (qvnow-ccover*qsliq(j,i,1))/(d_one-ccover)
          !---------------------------------------------
          ! humidity in moistest covpclr part of domain
          !---------------------------------------------
          qe = max(d_zero,min(qe,qsliq(j,i,1)))
          lactiv = covpclr > dlowval .and. &
          !      qpretot > dlowval .and. &
                 qxfg(iqqr) > activqx .and. qe < zrh*qsliq(j,i,1)
          if ( lactiv ) then
            ! note: units of preclr and qpretot differ
            !       qpretot is a mixing ratio (hence "q" in name)
            !       preclr is a rain flux
            preclr = qpretot*covpclr/(max(dlowval,covptot)*dtgdp)
            !--------------------------------------
            ! actual microphysics formula in beta
            !--------------------------------------
            ! sensitivity test showed multiply rain evap rate by 0.5
            beta1 = sqrt(phs(j,i,1)/pbot)/5.09e-3_rkx*preclr / &
                     max(covpclr,dlowval)

            if ( beta1 >= d_zero ) then
               beta = egrav*rpecons*d_half*(beta1)**0.5777_rkx
               denom = d_one + beta*dt*corqsliq
               dpr = covpclr * beta * (qsliq(j,i,1)-qe)/denom*dp*regrav
               dpevap = dpr*dtgdp
              !---------------------------------------------------------
              ! add evaporation term to explicit sink.
              ! this has to be explicit since if treated in the implicit
              ! term evaporation can not reduce rain to zero and model
              ! produces small amounts of rainfall everywhere.
              !---------------------------------------------------------

              ! evaporate rain
              chng = min(dpevap,qxfg(iqqr))
              !-------------------------------------------------------------
              ! reduce the total precip coverage proportional to evaporation
              !-------------------------------------------------------------
              covptot = covptot - max(d_zero,(covptot-ccover)*dpevap/qpretot)
            else
              chng = qxfg(iqqr)
            end if
            solqa(iqqv,iqqr) = solqa(iqqv,iqqr) + chng
            solqa(iqqr,iqqv) = solqa(iqqr,iqqv) - chng
            qxfg(iqqr)       = qxfg(iqqr) - chng
          end if

#ifdef DEBUG
          if ( stats ) then
            statsrainev(j,i,1) = chng
          end if
#endif

          ! snow

          chng = d_zero
          zrh = rprecrhmax + (d_one-rprecrhmax) * covpclr/(d_one-ccover)
          zrh = min(max(zrh,rprecrhmax),d_one)
          qe = (qvnow-ccover*qsice(j,i,1))/(d_one-ccover)
          !---------------------------------------------
          ! humidity in moistest covpclr part of domain
          !---------------------------------------------
          qe = max(d_zero,min(qe,qsice(j,i,1)))
          lactiv = covpclr > dlowval .and. &
                 qxfg(iqqs) > activqx .and. &
                 qe < zrh*qsice(j,i,1)
          if ( lactiv ) then
            ! note: units of preclr and qpretot differ
            !       qpretot is a mixing ratio (hence "q" in name)
            !       preclr is a rain flux
            preclr = qpretot*covpclr/(max(dlowval,covptot)*dtgdp)
            !--------------------------------------
            ! actual microphysics formula in beta
            !--------------------------------------
            beta1 = sqrt(phs(j,i,1)/pbot) / &
                  5.09e-3_rkx*preclr/max(covpclr,dlowval)

            if ( beta1 >= d_zero ) then
              beta = egrav*rpecons*(beta1)**0.5777_rkx
              !rpecons = alpha1
              denom = d_one + beta*dt*corqsice
              dpr = covpclr * beta * &
                     (qsice(j,i,1)-qe)/denom*dp*regrav
              dpevap = dpr*dtgdp
              !---------------------------------------------------------
              ! add evaporation term to explicit sink.
              ! this has to be explicit since if treated in the implicit
              ! term evaporation can not reduce snow to zero and model
              ! produces small amounts of snowfall everywhere.
              !---------------------------------------------------------
              ! evaporate snow
              chng = min(dpevap,qxfg(iqqs))
              chng = max(chng,d_zero)
              !-------------------------------------------------------------
              ! reduce the total precip coverage proportional to evaporation
              !-------------------------------------------------------------
              covptot = covptot-max(d_zero,(covptot-ccover) * &
                              dpevap/qpretot)
            else
              chng = qxfg(iqqs)
            end if
            solqa(iqqv,iqqs) = solqa(iqqv,iqqs) + chng
            solqa(iqqs,iqqv) = solqa(iqqs,iqqv) - chng
            qxfg(iqqs)       = qxfg(iqqs) - chng
          end if

#ifdef DEBUG
          if ( stats ) then
            statssnowev(j,i,1) = chng
          end if
#endif

        end if ! lmicro

        !--------------------------------
        ! solver for the microphysics
        !--------------------------------
        ! Truncate sum of explicit sinks to size of bin
        ! this approach is inaccurate, but conserves -
        ! prob best can do with explicit (i.e. not implicit!) terms
        !----------------------------------------------------------
        sinksum(:) = d_zero
        lind3(:,:) = .false.
        !----------------------------
        ! collect sink terms and mark
        !----------------------------
        do jn = 1 , nqx
          do n = 1 , nqx
            sinksum(n) = sinksum(n) - solqa(n,jn)
          end do
        end do
        !---------------------------------------
        ! calculate overshoot and scaling factor
        !---------------------------------------
        do n = 1 , nqx
          ratio(n) = qx0(n)/max(sinksum(n),qx0(n))
        end do
        !--------------------------------------------------------
        ! now sort ratio to find out which species run out first
        !--------------------------------------------------------
        do n = 1 , nqx
          iorder(n) = -999
          lind1(n) = .true.
        end do
        do n = 1 , nqx
          do jn = 1 , nqx
            if ( lind1(jn) .and. ratio(jn) > dlowval ) then
              iorder(n) = jn
            end if
          end do
        end do
        do n = 1 , nqx
          if ( iorder(n) > 0 ) then
            lind1(iorder(n)) = .false.
          end if
        end do
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
            if ( jo > 0 ) then
              lind3(jo,jn) = solqa(jo,jn) < d_zero
              sinksum(jo) = sinksum(jo) - solqa(jo,jn)
            end if
          end do
        end do
        !---------------------------
        ! recalculate scaling factor
        !---------------------------
        do n = 1 , nqx
          jo = iorder(n)
          if ( jo > 0 ) then
            ratio(jo) = qx0(jo)/max(sinksum(jo),qx0(jo))
          end if
        end do
        !------
        ! scale
        !------
        do n = 1 , nqx
          do jn = 1 , nqx
            jo = iorder(n)
            if ( jo > 0 ) then
              if ( lind3(jo,jn) ) then
                solqa(jo,jn) = solqa(jo,jn)*ratio(jo)
                solqa(jn,jo) = solqa(jn,jo)*ratio(jo)
              end if
            end if
          end do
        end do

        ! SOLVE THE LINEAR SYSTEM

        ! Set the LHS of equation
        do n = 1 , nqx
          do jn = 1 , nqx
            ! Diagonals: microphysical sink terms+transport
            if ( jn == n ) then
              qlhs(jn,n) = d_one + rsemi*fallsink(n)
              do jo = 1 , nqx
                qlhs(jn,n) = qlhs(jn,n) + rsemi*solqb(jo,jn)
              end do
              ! Non-diagonals: microphysical source terms
            else
              ! Here is the delta T - missing from doc.
              qlhs(jn,n) = -rsemi*solqb(jn,n)
            end if
          end do
        end do

        ! Set the RHS of equation

        do n = 1 , nqx
          ! Sum the explicit source and sink
          rexplicit = d_zero
          do jn = 1 , nqx
            ! Positive, since summed over 2nd index
            rexplicit = rexplicit + solqa(n,jn)
            if ( jn /= n ) then
              rexplicit =  rexplicit - &
                      (d_one-rsemi)*qx0(n)*solqb(jn,n) + &
                      (d_one-rsemi)*qx0(jn)*solqb(n,jn)
            end if
            rexplicit = rexplicit - (d_one-rsemi)*qx0(n)*fallsink(n)
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
          pfplsx(n,j,i,2) = rsemi*fallsink(n) * qxn(n)*rdtgdp + &
                  (d_one-rsemi)*fallsink(n)*qx0(n)*rdtgdp ! kg/m2/s
          ! Calculate fluxes in and out of box for conservation of TL
          fluxq = convsrce(n) + fallsrce(n) - &
                  (fallsink(n)+convsink(n)) * qxn(n)
          ! Calculate the water variables tendencies
          qxtendc(n,j,i,1) = qxtendc(n,j,i,1) + (qxn(n)-qx0(n))*oneodt
          ! Calculate the temperature tendencies
          if ( iphase(n) == 1 ) then
            ttendc(j,i,1) = ttendc(j,i,1) + &
                            wlhvocp * (qxn(n)-qx0(n)-fluxq)*oneodt
          else if ( iphase(n) == 2 ) then
            ttendc(j,i,1) = ttendc(j,i,1) + &
                            wlhsocp * (qxn(n)-qx0(n)-fluxq)*oneodt
          end if
        end do

      end do ! jx : end of longitude loop
    end do   ! iy : end of latitude loop
    !
    !----------------------------------------------------------------------
    !                       START OF VERTICAL LOOP
    !----------------------------------------------------------------------
    !
    ! Loop over levels and points
    !
    do k = 2 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2

          !---------------------------------
          ! First guess microphysics
          !---------------------------------
          do n = 1 , nqx
            qxfg(n) = qx(n,j,i,k)
            qx0(n) = qxfg(n)
          end do

          qvnow = qx0(iqqv)
          qlnow = qx0(iqql)
          qinow = qx0(iqqi)

          solqa(:,:)  = d_zero
          solqb(:,:)  = d_zero
          supsat      = d_zero
          subsat      = d_zero
          convsrce(:) = d_zero
          fallsrce(:) = d_zero
          fallsink(:) = d_zero
          convsink(:) = d_zero
          qpretot     = d_zero
          covptot     = d_zero
          covpclr     = d_zero
          qicetot     = d_zero
          ldefr       = d_zero

          pbot    = pfs(j,i,kzp1)
          dp      = dpfs(j,i,k)
          tk      = t(j,i,k)
          tc      = tk - tzero
          dens    = rho(j,i,k)
          alfaw   = qliq(j,i,k)
          sqmix   = qsmix(j,i,k)
          ccover  = fccfg(j,i,k)
          lccover = fccfg(j,i,k-1)
          totliq  = qlt(j,i,k)
          if ( lccn) ccn = pccn(j,i,k)
          rainp = pfplsx(iqqr,j,i,k)
          snowp = pfplsx(iqqs,j,i,k)

          ltkgt0    = ( tk > tzero )
          ltklt0    = ( .not. ltkgt0 )
          ltkgthomo = ( tk > thomo )
          lcloud    = ( ccover > zerocf )
          lnocloud  = ( .not. lcloud )
          locast    = ( ccover >= onecf )
          lnocast   = ( .not. locast )
          lliq      = ( totliq > activqx )

          ! Derived variables needed
          gdp = egrav/dp        ! g/dp  =(1/m)
          dtgdp = dt*gdp        ! (dt*g)/dp =(dt/m)
          rdtgdp = dp*(d_one/(dt*egrav)) ! dp/(gdt)=m/dt  [Kg/m2/s]
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
          facl     = alfaw*facw + (d_one - alfaw)*faci
          corr     = d_one/(d_one - ep1*eewmt(j,i,k))
          dqsmixdt = facl*corr*sqmix
          corqsmix = d_one/(d_one + eldcpm(tk)*dqsmixdt)
          !--------------------------------
          ! evaporation/sublimation limits
          !--------------------------------
          evaplimmix = max((sqmix-qvnow) * corqsmix,d_zero)
          !--------------------------------
          ! in-cloud consensate amount
          !--------------------------------
          tmpa = d_one/ccover
          liqcld = qlnow*tmpa
          icecld = qinow*tmpa
          licld  = liqcld + icecld
          !------------------------------------------------
          ! Evaporate very small amounts of liquid and ice
          !------------------------------------------------
          if ( qlnow < activqx ) then
            solqa(iqqv,iqql) =  qlnow
            solqa(iqql,iqqv) = -qlnow
            qxfg(iqql) = qxfg(iqql) - qlnow
          end if
          if ( qinow < activqx ) then
            solqa(iqqv,iqqi) =  qinow
            solqa(iqqi,iqqv) = -qinow
            qxfg(iqqi) = qxfg(iqqi) - qinow
          end if

          !------------------------------------------------------------------
          !  ICE SUPERSATURATION ADJUSTMENT
          !------------------------------------------------------------------
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
          !-------------------------------------------------------------------
          !-----------------------------------
          ! Supersaturation limit (from Koop)
          !-----------------------------------
          if ( nssopt == 0 )  then
            facl = d_one
          else
            if ( ltkgt0 ) then
              facl = d_one
            else
              koop = fkoop(tk,st_eeliq(j,i,k),st_eeice(j,i,k))
              facl = ccover + koop*(d_one-ccover)
            end if
          end if

          !-------------------------------------------------------------------
          ! Calculate supersaturation wrt Koop including dqs/dT
          ! correction factor
          !-------------------------------------------------------------------
          ! Here the supersaturation is turned into liquid water
          ! However, if the temperature is below the threshold for homogeneous
          ! freezing then the supersaturation is turned instantly to ice.
          ! Moreover the RH is clipped to the limit of
          ! qv_max = qs * (fcc + (1-fcc) *RH_homo )
          !--------------------------------------------------------------------
          ! qsmix?
          ! corqslsliq? shouldn't it be corqsmix?
          supsat = max((qvnow-facl*sqmix)/corqsliq,d_zero)
          ! e < esi, because for e > esi ice still present
          subsat = min((qvnow-facl*sqmix)/corqsliq,d_zero)

          if ( supsat > dlowval ) then
            if ( ltkgthomo ) then
              ! turn supersaturation into liquid water
              solqa(iqql,iqqv) = solqa(iqql,iqqv)+supsat
              solqa(iqqv,iqql) = solqa(iqqv,iqql)-supsat
              qxfg(iqql) = qxfg(iqql)+supsat
            else
              ! turn supersaturation into ice water
              solqa(iqqi,iqqv) = solqa(iqqi,iqqv)+supsat
              solqa(iqqv,iqqi) = solqa(iqqv,iqqi)-supsat
              qxfg(iqqi) = qxfg(iqqi)+supsat
            end if
          else
            if ( subsat < d_zero .and. lnocloud .and. lliq ) then
              ! turn subsaturation into vapor, where there is no cloud
              excess = totliq+subsat
              if ( excess < d_zero ) then
                if ( ltkgthomo ) then
                  evapl = max(-totliq,-evaplimmix)!*oneodt
                  solqa(iqqv,iqql) = solqa(iqqv,iqql)-evapl
                  solqa(iqql,iqqv) = solqa(iqql,iqqv)+evapl
                  qxfg(iqql) = qxfg(iqql)+evapl
                else
                  evapi = max(-totliq,-evaplimmix)!*oneodt
                  ! turn subsaturation into vapour
                  solqa(iqqv,iqqi) = solqa(iqqv,iqqi)-evapi
                  solqa(iqqi,iqqv) = solqa(iqqi,iqqv)+evapi
                  qxfg(iqqi) = qxfg(iqqi)+evapi
                end if
              end if
            end if
          end if

#ifdef DEBUG
          if ( stats ) then
            if ( supsat > dlowval ) then
              if ( ltkgthomo ) then
                statssupw(j,i,k) = supsat
              else
                statssupc(j,i,k) = supsat
              end if
            else
              if ( subsat < d_zero .and. lnocloud .and. lliq ) then
                excess = totliq+subsat
                if ( excess < d_zero ) then
                  if ( ltkgthomo ) then
                    evapl = max(-totliq,-evaplimmix)!*oneodt
                    statssupw(j,i,k) =  statssupw(j,i,k) + evapl
                  else
                    evapi = max(-totliq,-evaplimmix)!*oneodt
                    statssupc(j,i,k) =  statssupc(j,i,k) - evapi
                  end if
                end if
              end if
            end if
          end if
#endif
          !
          ! call addpath(iqql,iqqv,supsatl,solqa,solqb,d_zero,qxfg)
          ! call addpath(iqqi,iqqv,supsati,solqa,solqb,d_zero,qxfg)
          !
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
          ! SOLQA/B:q(IQa,IQb)
          !
          ! Thus if SOLQA/B(IQL,IQV) = K where K > 0 then this is
          ! a source of IQL and a sink of IQV
          !
          ! put 'magic' source terms such as LUDE from
          ! detrainment into explicit source/sink array diagnognal
          ! SOLQA(IQL,IQL)=LUDE
          !--------------------------------------------------------
          ! Define the microphysics
          ! the matrix will be sparse is this a problem ?
          ! (X,Y) means a sink of X and a source of Y
          ! for the implementation I will use flexible pointers
          ! such that it will be written (IQR,IQG) to indicate graupel to rain
          ! and the parametrization can have different variables switched on
          ! and off.
          ! each of these is a parametrization for a microphysical process.
          !------------------------------------------------------------------
          !                 DETRAINMENT FROM CONVECTION
          !------------------------------------------------------------------
          !chng = d_zero
          if ( k < kz ) then
            xqdetr(j,i,k) = xqdetr(j,i,k)*dtgdp  !kg/kg
            if ( xqdetr(j,i,k) > activqx ) then
              !qice = 1 if T < 250, qice = 0 if T > 273
              qice   = d_one-alfaw
              convsrce(iqql) = alfaw*xqdetr(j,i,k)
              convsrce(iqqi) = qice*xqdetr(j,i,k)
              solqa(iqql,iqql) = solqa(iqql,iqql) + convsrce(iqql)
              solqa(iqqi,iqqi) = solqa(iqqi,iqqi) + convsrce(iqqi)
              qxfg(iqql) = qxfg(iqql)+convsrce(iqql)
              qxfg(iqqi) = qxfg(iqqi)+convsrce(iqqi)
            end if
          end if

#ifdef DEBUG
          if ( stats ) then
            if ( k < kz ) then
              statsdetrw(j,i,k) = convsrce(iqql)
              statsdetrc(j,i,k) = convsrce(iqqi)
            end if
          end if
#endif

          !------------------------------------------------------------------
          ! Turn on/off microphysics
          !------------------------------------------------------------------

          if ( lmicro ) then

            !---------------------------------------
            ! EROSION OF CLOUDS BY TURBULENT MIXING
            !--------------------------------------
            ! rcldiff  : Diffusion coefficient for evaporation by turbulent
            ! mixing (IBID., EQU. 30) rcldiff = 3.0e-6_rkx
            ldifdt = rcldiff*dt
            if ( xqdetr(j,i,k) > d_zero ) ldifdt = 5.0_rkx*ldifdt
            !Increase by factor of 5 for convective points
            if ( lliq ) then
              leros = ccover * ldifdt * max(sqmix-qvnow,d_zero)
              leros = min(leros,evaplimmix)
              leros = min(leros,totliq)
              facl = qliqfrac(j,i,k)*leros
              faci = qicefrac(j,i,k)*leros
              solqa(iqql,iqqv) = solqa(iqql,iqqv) - facl
              solqa(iqqv,iqql) = solqa(iqqv,iqql) + facl
              solqa(iqqi,iqqv) = solqa(iqqi,iqqv) - faci
              solqa(iqqv,iqqi) = solqa(iqqv,iqqi) + faci
              qxfg(iqql) = qxfg(iqql) - facl
              qxfg(iqqi) = qxfg(iqqi) - faci
            end if

#ifdef DEBUG
            if ( stats ) then
              if ( lliq ) then
                statserosw(j,i,k) = qliqfrac(j,i,k)*leros
                statserosc(j,i,k) = qicefrac(j,i,k)*leros
              end if
            end if
#endif
            !------------------------------------------------------------------
            ! CONDENSATION/EVAPORATION DUE TO DQSAT/DT
            !------------------------------------------------------------------
            ! calculate dqs/dt
            ! Note: For the separate prognostic Qi and Ql, one would ideally use
            ! Qsat/DT wrt liquid/Koop here, since the physics is that new clouds
            ! forms by liquid droplets [liq] or when aqueous aerosols [Koop]
            ! form.
            ! These would then instant. freeze if T<-38C or lead to ice growth
            ! by deposition in warmer mixed phase clouds.  However, since we do
            ! not have a separate prognostic equation for in-cloud humidity or a
            ! statistical scheme approach in place, the depositional growth of
            ! ice in the mixed phase can not be modelled and we resort to
            ! supersaturation
            ! wrt ice instanteously converting to ice over one timestep
            ! (see Tompkins et al. QJRMS 2007 for details)
            ! Thus for the initial implementation the diagnostic mixed phase is
            ! retained for the moment, and the level of approximation noted.
            !------------------------------------------------------------------
            dtdp   = rovcp*tk/phs(j,i,k)
            dpmxdt = dp*oneodt
            wtot   = pverv(j,i,k)
            wtot   = min(dpmxdt,max(-dpmxdt,wtot))
            dtdiab = min(dpmxdt*dtdp,max(-dpmxdt*dtdp,radheatrt(j,i,k)))*dt + &
                          wlhfocp*ldefr
            ! ldefr = 0
            ! note: ldefr should be set to the difference between the mixed
            ! phase functions in the convection and cloud scheme, and
            ! for now we set it to zero and the functions are the same.
            ! In RegCM not all convection schemes provide such info.
            dtforc = dtdp*wtot*dt + dtdiab
            qold   = sqmix
            told   = tk
            tcond  = tk+dtforc
            tcond  = max(tcond,160.0_rkx)
            ! the goal is to produce dqs = qsmix - qold, where qsmix is
            ! reduced because of the condensation. so that dqs is negative?
            qp = d_one/phs(j,i,k)
            phases = max(min(d_one,((max(rtice,min(tzero, &
                       tcond))-rtice)*rtwat_rtice_r)**2),d_zero)
            ! saturation mixing ratio ws
            qsat = eewm(tcond,phases)*qp
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
            tcond = tcond+eldcpm(tcond)*cond1
            sqmix = sqmix - cond1
            dqs = sqmix - qold
            sqmix = qold
            tcond = told

            !----------------------------------------------------------------
            ! dqs > 0:  evaporation of clouds
            !----------------------------------------------------------------
            ! erosion term is linear in l
            ! changed to be uniform distribution in cloud region
            ! previous function based on delta distribution in cloud:
            if ( dqs > d_zero ) then
              !levap = C*min( dqs/dt , (qi+ql)/C )
              levap = ccover*min(dqs,licld)
              levap = min(levap,evaplimmix)
              levap = min(levap,max(sqmix-qvnow,d_zero))
              facl = qliqfrac(j,i,k)*levap
              faci = qicefrac(j,i,k)*levap
              solqa(iqqv,iqql) = solqa(iqqv,iqql) + facl
              solqa(iqql,iqqv) = solqa(iqql,iqqv) - facl
              solqa(iqqv,iqqi) = solqa(iqqv,iqqi) + faci
              solqa(iqqi,iqqv) = solqa(iqqi,iqqv) - faci
              qxfg(iqql) = qxfg(iqql) - facl
              qxfg(iqqi) = qxfg(iqqi) - faci
            end if

#ifdef DEBUG
            if ( stats ) then
              if ( dqs > d_zero ) then
                statsevapw(j,i,k) = qliqfrac(j,i,k)*levap
                statsevapc(j,i,k) = qicefrac(j,i,k)*levap
              end if
            end if
#endif

            !-----------------------------------------------------------------
            ! dqs < 0: formation of clouds
            !-----------------------------------------------------------------
            ! (1) increase of cloud water in existing clouds
            chng = d_zero

            if ( lcloud .and. dqs <= -activqx ) then
              ! new limiter
              chng = max(-dqs,d_zero)
              ! old limiter
              !  (significantly improves upper tropospheric humidity rms)
              phases = alfaw
              if ( locast ) then
                corr = d_one/(d_one-ep1*sqmix)
                cdmax = (qvnow-sqmix)/(d_one + corr*sqmix*edem(tk,phases))
              else
                cdmax = (qvnow-ccover*sqmix)/ccover
              end if
              chng = max(min(chng,cdmax),d_zero)
              ! end old limiter
              chng = ccover*chng
              if ( chng < activqx ) chng = d_zero
              !----------------------------------------------------------------
              ! All increase goes into liquid unless so cold cloud
              ! homogeneously freezes
              ! include new liquid formation in first guess value, otherwise
              ! liquid remains at cold temperatures until next timestep.
              !----------------------------------------------------------------
              if ( ltkgthomo ) then
                solqa(iqql,iqqv) = solqa(iqql,iqqv) + chng
                solqa(iqqv,iqql) = solqa(iqqv,iqql) - chng
                qxfg(iqql) = qxfg(iqql) + chng
              else
                solqa(iqqi,iqqv) = solqa(iqqi,iqqv) + chng
                solqa(iqqv,iqqi) = solqa(iqqv,iqqi) - chng
                qxfg(iqqi) = qxfg(iqqi) + chng
              end if
            end if

#ifdef DEBUG
            if ( stats ) then
              if ( ltkgthomo ) then
                statscond1w(j,i,k) = chng
              else
                statscond1c(j,i,k) = chng
              end if
            end if
#endif

            chng = d_zero
            qexc = d_zero

            ! (2) generation of new clouds (dc/dt>0)
            if ( dqs <= -activqx .and. lnocast ) then

              call selnss

              !---------------------------
              ! critical relative humidity
              !---------------------------
              ! *RAMID*   REAL    BASE VALUE FOR CALCULATION OF RELATIVE
              !                   HUMIDITY THRESHOLD FOR ONSET OF STRATIFORM
              !                   CONDENSATION (TIEDTKE, 1993, EQUATION 24)
              rhc = ramid !=0.8
              zsig = phs(j,i,k)/pbot
              ! increase RHcrit to 1.0 towards the surface (sigma>0.8)
              if ( zsig > ramid ) then
                rhc = ramid+(d_one-ramid)*((zsig-ramid)/0.2_rkx)**2
              end if
              !---------------------------
              ! supersaturation options
              !---------------------------
              if ( ltkgt0 .or. nssopt == 0 ) then
                ! no ice supersaturation allowed
                facl = d_one
              else
                ! ice supersaturation
                facl = fkoop(tk,st_eeliq(j,i,k),st_eeice(j,i,k))
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
                  xlcondlim = (ccover-d_one)*facl*dqs - &
                               facl*sqmix+qvnow
                  chng = min(chng,xlcondlim)
                end if
                chng = max(chng,d_zero)
                if ( chng < activqx ) then
                  chng = d_zero
                  acond = d_zero
                end if
                !-------------------------------------------------------------
                ! all increase goes into liquid unless so cold cloud
                ! homogeneously freezes
                ! include new liquid formation in first guess value, otherwise
                ! liquid remains at cold temperatures until next timestep.
                !-------------------------------------------------------------
                if ( ltkgthomo ) then
                  solqa(iqql,iqqv) = solqa(iqql,iqqv) + chng
                  solqa(iqqv,iqql) = solqa(iqqv,iqql) - chng
                  qxfg(iqql) = qxfg(iqql) + chng
                else ! homogeneous freezing
                  solqa(iqqi,iqqv) = solqa(iqqi,iqqv) + chng
                  solqa(iqqv,iqqi) = solqa(iqqv,iqqi) - chng
                  qxfg(iqqi) = qxfg(iqqi) + chng
                end if
              end if
            end if

#ifdef DEBUG
            if ( stats ) then
              if ( ltkgthomo ) then
                statscond2w(j,i,k) = chng
              else
                statscond2c(j,i,k) = chng
              end if
            end if
#endif

            !------------------------------------------------------------------
            ! DEPOSITION: Growth of ice by vapour deposition
            !------------------------------------------------------------------
            ! Following Rotstayn et al. 2001:
            ! does not use the ice nuclei number from cloudaer.F90
            ! but rather a simple Meyers et al. 1992 form based on the
            ! supersaturation and assuming clouds are saturated with
            ! respect to liquid water (well mixed), (or Koop adjustment)
            ! Growth considered as sink of liquid water if present so
            ! Bergeron-Findeisen adjustment in autoconversion term no
            ! longer needed
            !-----------------------------------------------------------------

            chng = d_zero

            !--------------------------------------------------------------
            ! only treat depositional growth if liquid present. due to fact
            ! that can not model ice growth from vapour without additional
            ! in-cloud water vapour variable
            ! If water droplets are present, the ice crystals are in an
            ! environment supersaturated with respect to ice and grow by
            ! deposition, reducing the water vapour and leading to
            ! subsaturation with respect to water.
            ! Water droplets then evaporate and the process continues with
            ! ice growth until the water droplets are completely evaporated.
            ! Thus in mixed phase clouds, the deposition process acts as a
            ! sink of cloud liquid and a source of ice cloud.
            !--------------------------------------------------------------
            if ( ltklt0 .and. qxfg(iqql) > activqx ) then
              vpice = st_eeice(j,i,k) !saturation vapor pressure wrt ice
              vpliq = st_eeliq(j,i,k) !saturation vapor pressure wrt liq
              ! Meyers et al 1992
              icenuclei = d_1000*exp(12.96_rkx*((vpliq-vpice)/vpice)-0.639_rkx)
              !------------------------------------------------
              !   2.4e-2 is conductivity of air
              !   8.87 = 700**1/3 = density of ice to the third
              !------------------------------------------------
              xadd  = wlhs*(wlhs/(rwat*tk)-d_one)/(2.4e-2_rkx*tk)
              xbdd  = rwat*tk*phs(j,i,k)/(2.21_rkx*vpice)
              cvds = 7.8_rkx*(icenuclei/dens)** &
                     0.666_rkx*(vpliq-vpice)/(8.87_rkx*(xadd+xbdd)*vpice)
              !-----------------------------------------------------
              ! iceinit = 1.e-12 is initial mass of ice particle
              !-----------------------------------------------------
              qice0 = max(icecld, icenuclei*iceinit/dens)
              !------------------
              ! new value of ice condensate amount ( Rotstayn et al. (2000) )
              !------------------
              qice0 = max(qice0,d_zero)
              cvds = max(cvds,d_zero)
              qinew = (0.666_rkx*cvds*dt+qice0**0.666_rkx)**1.5_rkx
              !---------------------------
              ! grid-mean deposition rate:
              !---------------------------
              chng = max(ccover*(qinew-qice0),d_zero)*2.0_rkx
              ! above increased by factor of 2 to retain similar mixed
              ! phase liq as in diagnostic scheme
              !---------------------------------------------------------------
              ! limit deposition to liquid water amount
              ! if liquid is all frozen, ice would use up reservoir of water
              ! vapour in excess of ice saturation mixing ratio - however this
              ! can not be represented without a in-cloud humidity variable.
              ! using the grid-mean humidity would imply a large artificial
              ! horizontal flux from the clear sky to the cloudy area.
              ! we thus rely on the supersaturation check to clean up any
              ! remaining supersaturation
              !---------------------------------------------------------------
              ! limit to liquid water amount
              chng = min(chng,qxfg(iqql))
              !---------------------------------------------------------------
              ! at top of cloud, reduce deposition rate near cloud top to
              ! account for small scale turbulent processes, limited ice
              ! nucleation and ice fallout
              !---------------------------------------------------------------
              ! Fraction of deposition rate in cloud top layer
              ! depliqrefrate = 0.1_rkx
              ! Depth of supercooled liquid water layer (m)
              ! depliqrefdepth = 500.0_rkx
              infactor = min(icenuclei/15000.0_rkx, d_one)
              chng = chng*min(infactor + (d_one-infactor)* &
                     (depliqrefrate+cldtopdist(j,i,k)/depliqrefdepth),d_one)
              !--------------
              ! add to matrix
              !--------------
              solqa(iqqi,iqql) = solqa(iqqi,iqql) + chng
              solqa(iqql,iqqi) = solqa(iqql,iqqi) - chng
              qxfg(iqql) = qxfg(iqql) - chng
              qxfg(iqqi) = qxfg(iqqi) + chng
            end if

#ifdef DEBUG
            if ( stats ) then
              statsdepos(j,i,k) = chng
            end if
#endif

            tmpa = d_one/ccover
            liqcld = qxfg(iqql)*tmpa
            icecld = qxfg(iqqi)*tmpa

            !------------------------------------------------------------------
            !  SEDIMENTATION/FALLING OF *ALL* MICROPHYSICAL SPECIES
            !     now that rain and snow species are prognostic
            !     the precipitation flux can be defined directly level by level
            !     There is no vertical memory required from the flux variable
            !------------------------------------------------------------------
            do n = 1 , nqx
              if ( lfall(n) ) then
                ! Source from layer above
                fallsrce(n) = pfplsx(n,j,i,k)*dtgdp
                solqa(n,n) = solqa(n,n) + fallsrce(n)
                qxfg(n) = qxfg(n) + fallsrce(n)
                qpretot = qpretot + qxfg(n)
                ! Sink to next layer, constant fall speed
                fallsink(n) = dtgdp*vqx(n)*dens   !Kg/Kg
              end if  !lfall
            end do ! n

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
            if ( qpretot > dlowval ) then
              covptot = d_one - ((d_one-covptot) * &
                              (d_one - max(ccover,lccover)) / &
                              (d_one - min(lccover,hicld)))
              covptot = max(covptot,rcovpmin)   !rcovpmin = 0.1
              ! clear sky proportion
              covpclr = max(d_zero,covptot-ccover)
            else
              covptot = d_zero  ! no flux - reset cover
              covpclr = d_zero  ! reset clear sky proportion
            end if
            !---------------------------------------------------------------
            !                         AUTOCONVERSION
            !---------------------------------------------------------------
            ! Warm clouds
            if ( liqcld > activqx ) then
              call selautoconv
            end if

            ! Cold clouds
            ! Snow Autoconversion rate follow Lin et al. 1983
            if ( ltklt0 ) then
              if ( icecld > activqx ) then
                alpha1 = dt*1.0e-3_rkx*exp(0.025*tc)
                xlcrit = rlcritsnow
                arg = (icecld/xlcrit)**2
                if ( arg < 25.0_rkx ) then
                  snowaut = alpha1 * (d_one - exp(-arg))
                else
                  snowaut = alpha1
                end if
                solqb(iqqs,iqqi) = solqb(iqqs,iqqi) + snowaut
              end if
            end if

            !---------------------------------------------------------------
            !                         MELTING
            !---------------------------------------------------------------
            ! The melting of ice and snow are treated explicitly.
            ! First water and ice saturation are found
            !---------------------------------------------
            ! ice saturation T < 273K
            ! liquid water saturation for T > 273K
            !---------------------------------------------

            do n = 1 , nqx
              if ( iphase(n) == 2 ) then
                qicetot = qicetot + qxfg(n)
              end if
            end do

            chngmax = d_zero

            if ( qicetot > activqx .and. ltkgt0 ) then
              ! Calculate subsaturation
              ! qsice(j,i,k)-qvnow,d_zero)
              subsat = max(sqmix-qvnow,d_zero)
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
              tdiff = tc ! - subsat * &
              !    (tw1+tw2*(phs(j,i,k)-tw3)-tw4*(tk-tw5))
              ! Ensure CONS1 is positive so that MELTMAX = 0 if TDMTW0 < 0
              cons1 = d_one ! abs(dt*(d_one + d_half*tdiff)/rtaumel)
              chngmax = max(tdiff*cons1*rldcp,d_zero)
            end if

            ! Loop over frozen hydrometeors (iphase == 2 (ice, snow))

            chng = d_zero
            if ( chngmax > dlowval .and. qicetot > activqx ) then
              do n = 1, nqx
                if ( iphase(n) == 2 ) then
                  m = imelt(n) ! imelt(iqqi)=iqql, imelt(iqqs)=iqqr
                  if ( m < 0 ) cycle
                  phases = qxfg(n)/qicetot
                  chng = min(qxfg(n),phases*chngmax)
                  chng = max(chng,d_zero)
                  ! n = iqqi,iqqs; m = iqql,iqqr
                  qxfg(n) =  qxfg(n) - chng
                  qxfg(m) =  qxfg(m) + chng
                  solqa(m,n) =  solqa(m,n) + chng
                  solqa(n,m) =  solqa(n,m) - chng
                end if
              end do
            end if

#ifdef DEBUG
            if ( stats ) then
              statsmelt(j,i,k) = chng
            end if
#endif
            !------------------------------------------------------------!
            !                         FREEZING                           !
            !------------------------------------------------------------!

            ! Freezing of rain.
            ! All rain freezes in a timestep if the temperature is below 0 C
            ! calculate sublimation latent heat

            chngmax = max((tzero-tk)*rldcp,d_zero)
            chng = d_zero

            if ( chngmax > dlowval .and. qxfg(iqqr) > activqx ) then
              chng = min(qxfg(iqqr),chngmax)
              chng = max(chng,d_zero)
              solqa(iqqs,iqqr) = solqa(iqqs,iqqr) + chng
              solqa(iqqr,iqqs) = solqa(iqqr,iqqs) - chng
            end if

#ifdef DEBUG
            if ( stats ) then
              statsfrz(j,i,k) = chng
            end if
#endif

            ! Freezing of liquid

            chngmax = max((thomo-tk)*rldcp,d_zero)
            chng = d_zero

            if ( chngmax > dlowval .and. qxfg(iqql) > activqx ) then
              chng = min(qxfg(iqql),chngmax)
              chng = max(chng,d_zero)
              solqa(iqqi,iqql) = solqa(iqqi,iqql) + chng
              solqa(iqql,iqqi) = solqa(iqql,iqqi) - chng
              qxfg(iqql) = qxfg(iqql) - chng
              qxfg(iqqi) = qxfg(iqqi) + chng
            end if

#ifdef DEBUG
            if ( stats ) then
              statsfrz(j,i,k) = statsfrz(j,i,k) + chng
            end if
#endif

            !---------------------------------------------------------------
            !                         EVAPORATION
            !---------------------------------------------------------------
            !------------------------------------------------------------
            ! recalculate qpretot since melting term may have changed it
            ! rprecrhmax = 0.7 is the threshold for the clear-sky RH that
            ! can be reached by evaporation of precipitation. THis assumption
            ! is done to prevent the gridbox saturating due to the evaporation
            ! of precipitation occuring in a portion of the grid
            !------------------------------------------------------------

            ! if ( iphase(n) == 2 ) then
            ! qpretot = qpretot+qxfg(n)
            qpretot = qxfg(iqqs) + qxfg(iqqr)
            ! end if

            ! rain
            chng = d_zero
            zrh = rprecrhmax + &
              (d_one-rprecrhmax)*covpclr/(d_one-ccover)
            zrh = min(max(zrh,rprecrhmax),d_one)
            ! This is a critical relative humidity that is used to limit
            ! moist environment to prevent the gridbox saturating when
            ! only part of the gridbox has evaporating precipitation
            qe = (qvnow-ccover*qsliq(j,i,k)) / (d_one-ccover)
            !---------------------------------------------
            ! humidity in moistest covpclr part of domain
            !---------------------------------------------
            qe = max(d_zero,min(qe,qsliq(j,i,k)))
            lactiv = covpclr > dlowval .and. &
            !      qpretot > dlowval .and. &
                   qxfg(iqqr) > activqx .and. qe < zrh*qsliq(j,i,k)
            if ( lactiv ) then
              ! note: units of preclr and qpretot differ
              !       qpretot is a mixing ratio (hence "q" in name)
              !       preclr is a rain flux
              preclr = qpretot*covpclr/(max(dlowval,covptot)*dtgdp)
              !--------------------------------------
              ! actual microphysics formula in beta
              !--------------------------------------
              ! sensitivity test showed multiply rain evap rate by 0.5
              beta1 = sqrt(phs(j,i,k)/pbot)/5.09e-3_rkx*preclr / &
                       max(covpclr,dlowval)

              if ( beta1 >= d_zero ) then
                 beta = egrav*rpecons*d_half*(beta1)**0.5777_rkx
                 denom = d_one + beta*dt*corqsliq
                 dpr = covpclr * beta * (qsliq(j,i,k)-qe)/denom*dp*regrav
                 dpevap = dpr*dtgdp
                !---------------------------------------------------------
                ! add evaporation term to explicit sink.
                ! this has to be explicit since if treated in the implicit
                ! term evaporation can not reduce rain to zero and model
                ! produces small amounts of rainfall everywhere.
                !---------------------------------------------------------

                ! evaporate rain
                chng = min(dpevap,qxfg(iqqr))
                !-------------------------------------------------------------
                ! reduce the total precip coverage proportional to evaporation
                !-------------------------------------------------------------
                covptot = covptot - max(d_zero, &
                           (covptot-ccover)*dpevap/qpretot)
              else
                chng = qxfg(iqqr)
              end if
              solqa(iqqv,iqqr) = solqa(iqqv,iqqr) + chng
              solqa(iqqr,iqqv) = solqa(iqqr,iqqv) - chng
              qxfg(iqqr)       = qxfg(iqqr) - chng
            end if

#ifdef DEBUG
            if ( stats ) then
              statsrainev(j,i,k) = chng
            end if
#endif

            ! snow

            chng = d_zero

            zrh = rprecrhmax + (d_one-rprecrhmax) * covpclr/(d_one-ccover)
            zrh = min(max(zrh,rprecrhmax),d_one)
            qe = (qvnow - ccover*qsice(j,i,k))/(d_one-ccover)
            !---------------------------------------------
            ! humidity in moistest covpclr part of domain
            !---------------------------------------------
            qe = max(d_zero,min(qe,qsice(j,i,k)))
            lactiv = covpclr > dlowval .and. &
                   qxfg(iqqs) > activqx .and. &
                   qe < zrh*qsice(j,i,k)
            if ( lactiv ) then
              ! note: units of preclr and qpretot differ
              !       qpretot is a mixing ratio (hence "q" in name)
              !       preclr is a rain flux
              preclr = qpretot*covpclr/(max(dlowval,covptot)*dtgdp)
              !--------------------------------------
              ! actual microphysics formula in beta
              !--------------------------------------
               beta1 = sqrt(phs(j,i,k)/pbot) / &
                    5.09e-3_rkx*preclr/max(covpclr,dlowval)

              if ( beta1 >= d_zero ) then
                beta = egrav*rpecons*(beta1)**0.5777_rkx
                !rpecons = alpha1
                denom = d_one + beta*dt*corqsice
                dpr = covpclr * beta * &
                       (qsice(j,i,k)-qe)/denom*dp*regrav
                dpevap = dpr*dtgdp
                !---------------------------------------------------------
                ! add evaporation term to explicit sink.
                ! this has to be explicit since if treated in the implicit
                ! term evaporation can not reduce snow to zero and model
                ! produces small amounts of snowfall everywhere.
                !---------------------------------------------------------
                ! evaporate snow
                chng = min(dpevap,qxfg(iqqs))
                chng = max(chng,d_zero)
                !-------------------------------------------------------------
                ! reduce the total precip coverage proportional to evaporation
                !-------------------------------------------------------------
                covptot = covptot-max(d_zero,(covptot-ccover) * &
                                dpevap/qpretot)
              else
                chng = qxfg(iqqs)
              end if
              solqa(iqqv,iqqs) = solqa(iqqv,iqqs) + chng
              solqa(iqqs,iqqv) = solqa(iqqs,iqqv) - chng
              qxfg(iqqs)       = qxfg(iqqs) - chng
            end if

#ifdef DEBUG
            if ( stats ) then
              statssnowev(j,i,k) = chng
            end if
#endif

          end if ! lmicro

          !--------------------------------
          ! solver for the microphysics
          !--------------------------------
          ! Truncate sum of explicit sinks to size of bin
          ! this approach is inaccurate, but conserves -
          ! prob best can do with explicit (i.e. not implicit!) terms
          !----------------------------------------------------------
          sinksum(:) = d_zero
          lind3(:,:) = .false.
          !----------------------------
          ! collect sink terms and mark
          !----------------------------
          do jn = 1 , nqx
            do n = 1 , nqx
              sinksum(n) = sinksum(n) - solqa(n,jn)
            end do
          end do
          !---------------------------------------
          ! calculate overshoot and scaling factor
          !---------------------------------------
          do n = 1 , nqx
            ratio(n) = qx0(n)/max(sinksum(n),qx0(n))
          end do
          !--------------------------------------------------------
          ! now sort ratio to find out which species run out first
          !--------------------------------------------------------
          do n = 1 , nqx
            iorder(n) = -999
            lind1(n) = .true.
          end do
          do n = 1 , nqx
            do jn = 1 , nqx
              if ( lind1(jn) .and. ratio(jn) > dlowval ) then
                iorder(n) = jn
              end if
            end do
          end do
          do n = 1 , nqx
            if ( iorder(n) > 0 ) then
              lind1(iorder(n)) = .false.
            end if
          end do
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
              if ( jo > 0 ) then
                lind3(jo,jn) = solqa(jo,jn) < d_zero
                sinksum(jo) = sinksum(jo) - solqa(jo,jn)
              end if
            end do
          end do
          !---------------------------
          ! recalculate scaling factor
          !---------------------------
          do n = 1 , nqx
            jo = iorder(n)
            if ( jo > 0 ) then
              ratio(jo) = qx0(jo)/max(sinksum(jo),qx0(jo))
            end if
          end do
          !------
          ! scale
          !------
          do n = 1 , nqx
            do jn = 1 , nqx
              jo = iorder(n)
              if ( jo > 0 ) then
                if ( lind3(jo,jn) ) then
                  solqa(jo,jn) = solqa(jo,jn)*ratio(jo)
                  solqa(jn,jo) = solqa(jn,jo)*ratio(jo)
                end if
              end if
            end do
          end do

          ! SOLVE THE LINEAR SYSTEM

          ! Set the LHS of equation
          do n = 1 , nqx
            do jn = 1 , nqx
              ! Diagonals: microphysical sink terms+transport
              if ( jn == n ) then
                qlhs(jn,n) = d_one + rsemi*fallsink(n)
                do jo = 1 , nqx
                  qlhs(jn,n) = qlhs(jn,n) + rsemi*solqb(jo,jn)
                end do
                ! Non-diagonals: microphysical source terms
              else
                ! Here is the delta T - missing from doc.
                qlhs(jn,n) = -rsemi*solqb(jn,n)
              end if
            end do
          end do

          ! Set the RHS of equation

          do n = 1 , nqx
            ! Sum the explicit source and sink
            rexplicit = d_zero
            do jn = 1 , nqx
              ! Positive, since summed over 2nd index
              rexplicit = rexplicit + solqa(n,jn)
              if ( jn /= n ) then
                rexplicit =  rexplicit - &
                        (d_one-rsemi)*qx0(n)*solqb(jn,n) + &
                        (d_one-rsemi)*qx0(jn)*solqb(n,jn)
              end if
              rexplicit = rexplicit - (d_one-rsemi)*qx0(n)*fallsink(n)
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
            pfplsx(n,j,i,k+1) = rsemi*fallsink(n) * qxn(n)*rdtgdp + &
                    (d_one-rsemi)*fallsink(n)*qx0(n)*rdtgdp ! kg/m2/s
            ! Calculate fluxes in and out of box for conservation of TL
            fluxq = convsrce(n) + fallsrce(n) - &
                    (fallsink(n)+convsink(n)) * qxn(n)
            ! Calculate the water variables tendencies
            qxtendc(n,j,i,k) = qxtendc(n,j,i,k) + (qxn(n)-qx0(n))*oneodt
            ! Calculate the temperature tendencies
            if ( iphase(n) == 1 ) then
              ttendc(j,i,k) = ttendc(j,i,k) + &
                       wlhvocp * (qxn(n)-qx0(n)-fluxq)*oneodt
            else if ( iphase(n) == 2 ) then
              ttendc(j,i,k) = ttendc(j,i,k) + &
                       wlhsocp * (qxn(n)-qx0(n)-fluxq)*oneodt
            end if
          end do

        end do ! jx : end of longitude loop
      end do   ! iy : end of latitude loop
    end do     ! kz : end of vertical loop
    !
    ! Couple tendencies with pressure
    !
    do n = 1 , nqx
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            qxten(j,i,k,n) = qxtendc(n,j,i,k)*psb(j,i)
          end do
        end do
      end do
    end do
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          tten(j,i,k) = ttendc(j,i,k)*psb(j,i)
        end do
      end do
    end do

    !-------------------------------------
    ! Final enthalpy and total water diagnostics
    !-------------------------------------
    if ( budget_compute ) then
      ! Initialize the flux arrays
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            tnew = t(j,i,k)+dt*(ttendc(j,i,k)-tentkeep(j,i,k))
            if ( k == 1 ) then
              sumq1(j,i,k) = d_zero ! total water
              sumh1(j,i,k) = d_zero ! liquid water temperature
            else
              sumq1(j,i,k) = sumq1(j,i,k-1)
              sumh1(j,i,k) = sumh1(j,i,k-1)
            end if
            ! cld vars
            do n = 1 , nqx
              if ( iphase(n) == 1 ) then
                tnew = tnew-wlhvocp*(qx(n,j,i,k)+ &
                        (qxtendc(n,j,i,k)-tenkeep(n,j,i,k))*dt)
              else if ( iphase(n) == 2 ) then
                tnew = tnew-wlhsocp*(qx(n,j,i,k)+ &
                        (qxtendc(n,j,i,k)-tenkeep(n,j,i,k))*dt)
              end if
              sumq1(j,i,k) = sumq1(j,i,k) + &
                (qx(n,j,i,k)+(qxtendc(n,j,i,k)-tenkeep(n,j,i,k))*dt)* &
                dpfs(j,i,k)*regrav
            end do
            sumh1(j,i,k) = sumh1(j,i,k)+dpfs(j,i,k)*tnew
            rain = d_zero
            do n = 1 , nqx
            rain = rain + dt*pfplsx(n,j,i,k+1)
          end do
          errorq(j,i,k) = sumq1(j,i,k)+rain-sumq0(j,i,k)
        end do
      end do
    end do
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          dtgdp = dt*egrav/dpfs(j,i,k)
          rain = d_zero
          do n = 1 , nqx
            if ( iphase(n) == 1 ) then
              rain = rain+wlhvocp*dtgdp*pfplsx(n,j,i,k+1)*dpfs(j,i,k)
            else if ( iphase(n) == 2 ) then
              rain = rain+wlhsocp*dtgdp*pfplsx(n,j,i,k+1)*dpfs(j,i,k)
            end if
          end do
          sumh1(j,i,k) = (sumh1(j,i,k)-rain)/(pfs(j,i,k+1)-pfs(j,i,1))
          errorh(j,i,k) = sumh1(j,i,k)-sumh0(j,i,k)
        end do
      end do
    end do

    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( abs(errorq(j,i,kz)) > 1.e-12_rkx .or. &
               abs(errorh(j,i,kz)) > 1.e-12_rkx) then
            if ( abs(errorq(j,i,kz)) > 1.e-12_rkx ) then
              write(stderr,*) 'WATER NON CONSERVED AT '
              write(stderr,*) 'J = ',j
              write(stderr,*) 'I = ',i
              write(stderr,*) 'K = ',k
              write(stderr,*) 'ERROR IS : ',errorq(j,i,kz)
            end if
            if ( abs(errorh(j,i,kz)) > 1.e-12_rkx ) then
              write(stderr,*) 'ENTHALPY NON CONSERVED AT '
              write(stderr,*) 'J = ',j
              write(stderr,*) 'I = ',i
              write(stderr,*) 'K = ',k
              write(stderr,*) 'ERROR IS : ',errorh(j,i,kz)
            end if
            call fatal(__FILE__,__LINE__, &
              'TOTAL WATER OR ENTHALPY NOT CONSERVED')
          end if
        end do
      end do
    end do
  end if ! budget_compute

  ! Sum fluxes over the levels
  ! Initialize fluxes
  pfplsl(:,:,:) = d_zero
  pfplsn(:,:,:) = d_zero
  rainls(:,:,:) = d_zero

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
            rainls(j,i,k) = pfplsl(j,i,k)
          else if ( iphase(n) == 2 ) then
            pfplsn(j,i,k) = pfplsn(j,i,k)+ pfplsx(n,j,i,k)
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
      if ( prainx > dlowval ) then
        rainnc(j,i) =  rainnc(j,i) + prainx   !mm
        lsmrnc(j,i) =  lsmrnc(j,i) + pfplsl(j,i,kzp1)
      end if
      if ( psnowx > dlowval ) then
        snownc(j,i) = snownc(j,i) + psnowx
        lsmrnc(j,i) =  lsmrnc(j,i) + pfplsn(j,i,kzp1)
      end if
    end do
  end do

#ifdef DEBUG
  call time_end(subroutine_name,idindx)
#endif

  contains

    pure real(rkx) function delta(t)
      !delta = 1 if t > tzero
      !delta = 0 if t < tzero
      implicit none
      real(rkx) , intent(in):: t
      delta = max(d_zero,sign(d_one,t-tzero))
    end function delta

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
      eldcpm = phase*wlhvocp+(d_one-phase)*wlhsocp
    end function eldcpm

    subroutine eeliq ! = 0.622*esw  !Teten's formula
      implicit none
      integer(ik4) :: i , j , k
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            st_eeliq(j,i,k) = c2es*exp(c3les*((t(j,i,k)-tzero) / &
                                              (t(j,i,k)-c4les)))
          end do
        end do
      end do
    end subroutine eeliq

    subroutine eeice ! = 0.622*esi
      implicit none
      integer(ik4) :: i , j , k
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            st_eeice(j,i,k) = c2es*exp(c3ies*((t(j,i,k)-tzero) / &
                                              (t(j,i,k)-c4ies)))
          end do
        end do
      end do
    end subroutine eeice

    pure real(rkx) function eewm(t,phase)
      implicit none
      real(rkx) , intent(in):: t , phase
      real(rkx) :: eliq , eice
      eliq = c2es*exp(c3les*((t-tzero)/(t-c4les)))
      eice = c2es*exp(c3ies*((t-tzero)/(t-c4ies)))
      eewm = phase * eliq + (d_one-phase) * eice
    end function eewm

    pure real(rkx) function fkoop(t,eeliq,eeice)
      implicit none
      ! se T < 0 la nuvola si forma o quando q e' maggiore della liquid
      ! water saturation minima, oppure se e' maggiore del mixing ratio
      ! wrt ice critica a cui inizia l'homogeneaous ice nucleation
      ! At temperatures below 0◦C new cloud forms in any non-cloudy part
      ! of the grid box where the humidity exceeds either the minimum of
      ! the liquid water saturation specific humidity (qsl), or the
      ! critical vapour saturation mixing ratio with respect to ice at
      ! which homogeneous ice nucleation initiates
      ! empirical fit given by Karcher and Lohmann (2002) which is a
      ! function of temperature and ranges from 45% supersaturation at
      ! T = 235 K to 67% at T = 190 K.
      ! At temperatures warmer than -38 degC the cloud formation over a
      ! timestep results entirely in liquid cloud,
      ! i.e. fkoop = eeliq/eeice, mentre per T < -38 fkoop = RHhomo
      ! while below this threshold the liquid water or aqueous sulphate
      ! solutes are assumed to freeze instantaneously and the process is
      ! a source for cloud ice.
      ! fkoop modifies the ice saturation mixing ratio for homogeneous
      ! nucleation
      real(rkx) , parameter :: rkoop1 = 2.583_rkx
      real(rkx) , parameter :: rkoop2 = 0.48116e-2_rkx ! 1/207.8
      real(rkx) , intent(in) :: t , eeliq , eeice
      fkoop = min(rkoop1-rkoop2*t,eeliq/eeice)
    end function fkoop

    subroutine klein_and_pincus
      implicit none
      solqb(iqql,iqqv) = d_zero
      solqb(iqqr,iqql) = solqb(iqqr,iqql) + &
                         dt*auto_rate_klepi * (qlnow**(2.3_rkx))
      solqa(iqqr,iqql) = d_zero
    end subroutine klein_and_pincus

    subroutine khairoutdinov_and_kogan
      implicit none
      solqb(iqql,iqqv) = d_zero
      solqb(iqqr,iqql) = solqb(iqqr,iqql) + &
                         dt*auto_rate_khair*(qlnow**(auto_expon_khair))
    end subroutine khairoutdinov_and_kogan

    subroutine kessler
      implicit none
      solqb(iqql,iqqv) = d_zero
      solqa(iqqr,iqql) = solqa(iqqr,iqql) - &
                         auto_rate_kessl*autocrit_kessl
      solqa(iqql,iqqr) = solqa(iqql,iqqr) + &
                         auto_rate_kessl*autocrit_kessl
      solqb(iqqr,iqql) = solqb(iqqr,iqql) + dt*auto_rate_kessl
    end subroutine kessler

    subroutine sundqvist
      implicit none
      real(rkx) :: precip , cfpr , rainaut , arg
      solqb(iqql,iqqv) = d_zero
      alpha1 = rkconv*dt
      ! modify autoconversion threshold dependent on:
      ! land (polluted, high ccn, smaller droplets, higher
      !       threshold)
      ! sea  (clean, low ccn, larger droplets, lower threshold)
      if ( ldmsk(j,i) == 0 ) then  ! landmask =0 land, =1 ocean
        ! THRESHOLD VALUE FOR RAIN AUTOCONVERSION OVER LAND
        xlcrit = rclcrit_land ! landrclcrit_land = 5.e-4
      else
        xlcrit = rclcrit_sea  ! oceanrclcrit_sea  = 3.e-4
      end if
      if ( lccn ) then
        if ( ccn > 0._rkx ) then
          ! aerosol second indirect effect on autoconversion
          ! threshold, rcrit is a critical cloud radius for cloud
          ! water undergoing autoconversion
          ! ccn = number of ccn /m3
          xlcrit = ccn*(4.0_rkx/3.0_rkx)*mathpi * &
                   ((rcrit*1e-6_rkx)**3)*rhoh2o
        endif
      endif
      !-----------------------------------------------------------
      ! parameters for cloud collection by rain and snow.
      ! note that with new prognostic variable it is now possible
      ! to replace this with an explicit collection
      ! parametrization
      !-----------------------------------------------------------
      precip = (rainp+snowp)/max(dlowval,covptot)
      cfpr = d_one + rprc1*sqrt(max(precip,d_zero))
      alpha1 = alpha1*cfpr
      xlcrit = xlcrit/max(cfpr,dlowval)
      arg = (liqcld/xlcrit)**2
      ! security for exp for some compilers
      if ( arg < 25.0_rkx ) then
        rainaut = alpha1*(d_one - exp(-arg))
      else
        rainaut = alpha1
      end if
      !-----------------------
      ! rain freezes instantly
      !-----------------------
      if ( ltklt0 ) then
        solqb(iqqs,iqql) = solqb(iqqs,iqql)+rainaut
      else
        solqb(iqqr,iqql) = solqb(iqqr,iqql)+rainaut
      end if
    end subroutine sundqvist

    subroutine nss_tompkins
      implicit none
      qexc = max((qvnow-ccover*sqmix)/(d_one-ccover),d_zero)
    end subroutine nss_tompkins

    subroutine nss_lohmann_and_karcher
      implicit none
      qexc = qvnow
    end subroutine nss_lohmann_and_karcher

    subroutine nss_gierens
      implicit none
      qexc = qvnow/totliq
    end subroutine nss_gierens

    subroutine mysolve
      implicit none
      integer(ik4) :: ii , jj , kk , ll , imax
      real(rkx) :: aamax , dum , xsum , swap
      ! find implicit scaling information
      do n = 1 , nqx
        aamax = d_zero
        do jn = 1 , nqx
          if ( abs(qlhs(n,jn)) > aamax ) aamax = abs(qlhs(n,jn))
        end do
        if ( aamax < dlowval ) then
          call fatal(__FILE__,__LINE__, &
                     'SINGULAR MATRIX')
        end if ! Singular matrix
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
        ! Initialize for the search for largest pivot element.
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
        dum = d_one/qlhs(n,n)
        if ( n /= nqx ) then
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
!
!  subroutine addpath(src,snk,proc,zsqa,zsqb,beta,fg)
!    implicit none
!    real(rkx) , pointer , intent(inout) , dimension(:,:) :: zsqa , zsqb
!    real(rkx) , pointer , intent(inout) , dimension(:) :: fg
!    real(rkx) , intent(in) :: proc
!    integer(ik4) , intent(in) :: src , snk
!    real(rkx) , intent(in) :: beta
!    zsqa(src,snk) = zsqa(src,snk) + (d_one-beta)*proc
!    zsqa(snk,src) = zsqa(snk,src) - (d_one-beta)*proc
!    fg(src) = fg(src) + (d_one-beta)*proc
!    fg(snk) = fg(snk) - (d_one-beta)*proc
!    zsqb(src,snk) = zsqb(src,snk) + beta*proc
!  end subroutine addpath
!
  end subroutine microphys

end module mod_cloud_s1

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

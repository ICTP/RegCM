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
  use mod_runparams , only : stats , budget_compute , nssopt , iautoconv
  use mod_runparams , only : auto_rate_khair , auto_rate_kessl , &
                             auto_rate_klepi
  use mod_runparams , only : rsemi , rkconv , rcovpmin , rpecons
  use mod_runparams , only : ktau
  use mod_runparams , only : rtsrf

  implicit none

  private

  logical , parameter :: lmicro = .true.
  logical , parameter :: fscheme = .true.

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
  real(rkx) , parameter :: cldtopcf = d_r100
  ! Fraction of deposition rate in cloud top layer
  real(rkx) , parameter :: depliqrefrate = d_r10
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
  integer(ik4) , pointer , dimension(:,:,:) :: iorder
  logical , pointer , dimension(:) :: lfall
  logical , pointer , dimension(:,:,:) ::   lind1
  logical , pointer , dimension(:,:,:,:) :: lind3

  real(rkx) , pointer , dimension(:,:,:):: sumh0 , sumq0
  real(rkx) , pointer , dimension(:,:,:) :: sumh1 , sumq1
  real(rkx) , pointer , dimension(:,:,:) :: errorq , errorh
  real(rkx) , pointer , dimension(:,:,:):: tentkeep
  real(rkx) , pointer , dimension(:,:,:,:) :: tenkeep
  ! Mass variables
  real(rkx) , pointer , dimension(:,:) :: dp     ! dp
  real(rkx) , pointer , dimension(:,:) :: dtgdp  ! dt * g/dp
  real(rkx) , pointer , dimension(:,:) :: rdtgdp ! dp / (dt * g)  [Kg/(m*s)]
  ! Microphysics
  real(rkx) , pointer , dimension(:,:) :: qexc
  real(rkx) , pointer , dimension(:,:) :: chng
  real(rkx) , pointer , dimension(:,:) :: chngmax
  real(rkx) , pointer , dimension(:,:) :: qicetot
  real(rkx) , pointer , dimension(:,:) :: prcflxw
  real(rkx) , pointer , dimension(:,:) :: prcflxc
  real(rkx) , pointer , dimension(:,:,:) :: dqsatdt
  real(rkx) , pointer , dimension(:,:,:) :: pccn
  ! for sedimentation source/sink terms
  real(rkx) , pointer , dimension(:,:,:) :: fallsink
  real(rkx) , pointer , dimension(:,:,:) :: fallsrce
  ! for convection detrainment source and subsidence source/sink terms
  real(rkx) , pointer , dimension(:,:) :: corqsice
  real(rkx) , pointer , dimension(:,:) :: corqsliq
  real(rkx) , pointer , dimension(:,:) :: corqsmix
  real(rkx) , pointer , dimension(:,:) :: evaplimmix
  real(rkx) , pointer , dimension(:,:,:) :: convsrce
  real(rkx) , pointer , dimension(:,:,:) :: convsink
  ! total rain frac: fractional occurence of precipitation (%)
  real(rkx) , pointer , dimension(:,:) :: covptot
  ! for condensation
  real(rkx) , pointer , dimension(:,:) :: covpclr
  real(rkx) , pointer , dimension(:,:) :: qpretot
  real(rkx) , pointer , dimension(:,:) :: liqcld
  real(rkx) , pointer , dimension(:,:) :: icecld
  real(rkx) , pointer , dimension(:,:) :: supsat
  real(rkx) , pointer , dimension(:,:) :: subsat
  real(rkx) , pointer , dimension(:,:) :: licld
  real(rkx) , pointer , dimension(:,:) :: ldefr
  real(rkx) , pointer , dimension(:,:) :: qold
  real(rkx) , pointer , dimension(:,:) :: told
  real(rkx) , pointer , dimension(:,:) :: dqs
  real(rkx) , pointer , dimension(:,:,:) :: tcond
  ! distance from the top of the cloud
  real(rkx) , pointer , dimension(:,:) :: cldtopdist
  ! ice nuclei concentration
  real(rkx) , pointer , dimension(:,:) :: icenuclei
  real(rkx) , pointer , dimension(:,:,:) :: eewmt
  real(rkx) , pointer , dimension(:,:,:) :: qliq
  real(rkx) , pointer , dimension(:,:,:) :: qliqfrac
  real(rkx) , pointer , dimension(:,:,:) :: qicefrac
  real(rkx) , pointer , dimension(:,:,:) :: qlt
  ! fluxes convergence of species
  real(rkx) , pointer , dimension(:,:,:) :: fluxq
  real(rkx) , pointer , dimension(:,:,:) :: ratio
  real(rkx) , pointer , dimension(:,:,:) :: sinksum
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
  real(rkx) , pointer , dimension(:,:,:,:) :: solqa
  ! implicit sources and sinks
  real(rkx) , pointer , dimension(:,:,:,:) :: solqb
  ! decoupled mixing ratios tendency
  real(rkx) , pointer , dimension(:,:,:,:) :: qxtendc
  ! j,i,n ! generalized precipitation flux
  real(rkx) , pointer , dimension(:,:,:,:) :: pfplsx
  real(rkx) , public  , pointer, dimension(:,:,:,:) :: qx0
  ! new values for qxx at time+1
  real(rkx) , public  , pointer, dimension(:)   :: qxn
  ! first guess values including precip
  real(rkx) , public  , pointer, dimension(:,:,:)   :: qxfg
  ! first guess value for cloud fraction
  real(rkx) , public  , pointer, dimension(:,:,:)   :: fccfg
  ! relative humidity
  real(rkx) , public  , pointer, dimension(:,:,:)   :: relh
  ! saturation mixing ratio with respect to water
  real(rkx) , public  , pointer, dimension(:,:,:)   :: qsliq
  ! turbulent erosion rate
  real(rkx) , pointer , dimension(:,:) :: ldifdt

  ! statistic only if stas =.true.
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
  integer(ik4) , pointer , dimension(:) :: indx
  real(rkx) , pointer , dimension(:) :: vv

!  interface addpath
!    module procedure addpath_array
!    module procedure addpath_real
!  end interface

  real(rkx) , parameter :: clfeps = 1.0e-6_rkx
  real(rkx) , parameter :: zerocf = lowcld - clfeps
  real(rkx) , parameter :: onecf  = hicld + clfeps

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
    call getmem2d(ldefr,jci1,jci2,ici1,ici2,'cmicro:ldefr')
    call getmem2d(qold,jci1,jci2,ici1,ici2,'cmicro:qold')
    call getmem2d(told,jci1,jci2,ici1,ici2,'cmicro:told')
    call getmem2d(dqs,jci1,jci2,ici1,ici2,'cmicro:dqs')
    call getmem2d(qicetot,jci1,jci2,ici1,ici2,'cmicro:qicetot')
    call getmem2d(qexc,jci1,jci2,ici1,ici2,'cmicro:qexc')
    call getmem2d(chng,jci1,jci2,ici1,ici2,'cmicro:chng')
    call getmem2d(chngmax,jci1,jci2,ici1,ici2,'cmicro:chngmax')
    call getmem2d(dp,jci1,jci2,ici1,ici2,'cmicro:dp')
    call getmem2d(dtgdp,jci1,jci2,ici1,ici2,'cmicro:dtgdp')
    call getmem2d(rdtgdp,jci1,jci2,ici1,ici2,'cmicro:rdtgdp')
    call getmem2d(prcflxw,jci1,jci2,ici1,ici2,'cmicro:prcflxw')
    call getmem2d(prcflxc,jci1,jci2,ici1,ici2,'cmicro:prcflxc')
    call getmem2d(covptot,jci1,jci2,ici1,ici2,'cmicro:covptot')
    call getmem2d(covpclr,jci1,jci2,ici1,ici2,'cmicro:covpclr')
    call getmem2d(qpretot,jci1,jci2,ici1,ici2,'cmicro:qpretot')
    call getmem2d(cldtopdist,jci1,jci2,ici1,ici2,'cmicro:cldtopdist')
    call getmem2d(icenuclei,jci1,jci2,ici1,ici2,'cmicro:icenuclei')
    call getmem2d(licld,jci1,jci2,ici1,ici2,'cmicro:licld')
    call getmem2d(corqsmix,jci1,jci2,ici1,ici2,'cmicro:corqsmix')
    call getmem2d(evaplimmix,jci1,jci2,ici1,ici2,'cmicro:evaplimmix')
    call getmem2d(liqcld,jci1,jci2,ici1,ici2,'cmicro:liqcld')
    call getmem2d(icecld,jci1,jci2,ici1,ici2,'cmicro:icecld')
    call getmem2d(supsat,jci1,jci2,ici1,ici2,'cmicro:supsat')
    call getmem2d(subsat,jci1,jci2,ici1,ici2,'cmicro:subsat')
    call getmem2d(corqsice,jci1,jci2,ici1,ici2,'cmicro:corqsice')
    call getmem2d(corqsliq,jci1,jci2,ici1,ici2,'cmicro:corqsliq')
    call getmem2d(ldifdt,jci1,jci2,ici1,ici2,'cmicro:ldifdt')
    call getmem3d(tcond,jci1,jci2,ici1,ici2,1,kz,'cmicro:tcond')
    call getmem3d(qliqfrac,jci1,jci2,ici1,ici2,1,kz,'cmicro:qliqfrac')
    call getmem3d(qicefrac,jci1,jci2,ici1,ici2,1,kz,'cmicro:qicefrac')
    call getmem3d(eewmt,jci1,jci2,ici1,ici2,1,kz,'cmicro:eewmt')
    call getmem3d(qsmix,jci1,jci2,ici1,ici2,1,kz,'cmicro:qsmix')
    call getmem3d(qlt,jci1,jci2,ici1,ici2,1,kz+1,'cmicro:qlt')
    call getmem3d(iorder,1,nqx,jci1,jci2,ici1,ici2,'cmicro:iorder')
    call getmem3d(ttendc,jci1,jci2,ici1,ici2,1,kz,'cmicro:ttendc')
    call getmem3d(convsrce,1,nqx,jci1,jci2,ici1,ici2,'cmicro:convsrce')
    call getmem3d(convsink,1,nqx,jci1,jci2,ici1,ici2,'cmicro:convsink')
    call getmem3d(eew,jci1,jci2,ici1,ici2,1,kz,'cmicro:eew')
    call getmem3d(qsice,jci1,jci2,ici1,ici2,1,kz,'cmicro:qsice')
    call getmem4d(qx0,1,nqx,jci1,jci2,ici1,ici2,1,kz,'cmicro:qx0')
    call getmem3d(qsliq,jci1,jci2,ici1,ici2,1,kz,'cmicro:qsliq')
    call getmem3d(fallsink,1,nqx,jci1,jci2,ici1,ici2,'cmicro:fallsink')
    call getmem3d(fallsrce,1,nqx,jci1,jci2,ici1,ici2,'cmicro:fallsrce')
    call getmem3d(fluxq,1,nqx,jci1,jci2,ici1,ici2,'cmicro:fluxq')
    call getmem3d(ratio,1,nqx,jci1,jci2,ici1,ici2,'cmicro:ratio')
    call getmem3d(sinksum,1,nqx,jci1,jci2,ici1,ici2,'cmicro:sinksum')
    call getmem3d(dqsatdt,jci1,jci2,ici1,ici2,1,kz,'cmicro:dqsatdt')
    call getmem3d(pfplsl,jci1,jci2,ici1,ici2,1,kz+1,'cmicro:pfplsl')
    call getmem3d(pfplsn,jci1,jci2,ici1,ici2,1,kz+1,'cmicro:pfplsn')
    call getmem3d(pfsqlf,jci1,jci2,ici1,ici2,1,kz+1,'cmicro:pfsqlf')
    call getmem3d(pfsqif,jci1,jci2,ici1,ici2,1,kz+1,'cmicro:pfsqif')
    call getmem3d(qliq,jci1,jci2,ici1,ici2,1,kz+1,'cmicro:qliq')
    call getmem3d(qxfg,1,nqx,jci1,jci2,ici1,ici2,'cmicro:qxfg')
    call getmem3d(fccfg,jci1,jci2,ici1,ici2,1,kz,'cmicro:fccfg')
    call getmem3d(lind1,1,nqx,jci1,jci2,ici1,ici2,'cmicro:lind1')
    call getmem3d(xqdetr,jci1,jci2,ici1,ici2,1,kz,'cmicro:xqdetr')
    !xqdetr after evaporation
    !call getmem3d(xqdetr2,jci1,jci2,ici1,ici2,1,kz,'cmicro:xqdetr2')
    call getmem3d(st_eeliq,jci1,jci2,ici1,ici2,1,kz,'cmicro:st_eeliq')
    call getmem3d(st_eeice,jci1,jci2,ici1,ici2,1,kz,'cmicro:st_eeice')
    call getmem3d(eeliqt,jci1,jci2,ici1,ici2,1,kz,'cmicro:eeliqt')
    call getmem4d(qxtendc,1,nqx,jci1,jci2,ici1,ici2,1,kz,'cmicro:qxtendc')
    call getmem1d(qxn,1,nqx,'cmicro:qxn')
    call getmem2d(qlhs,1,nqx,1,nqx,'cmicro:qlhs')
    call getmem4d(solqa,1,nqx,1,nqx,jci1,jci2,ici1,ici2,'cmicro:solqa')
    call getmem4d(solqb,1,nqx,1,nqx,jci1,jci2,ici1,ici2,'cmicro:solqb')
    call getmem4d(lind3,1,nqx,1,nqx,jci1,jci2,ici1,ici2,'cmicro:lind3')
    call getmem4d(pfplsx,1,nqx,jci1,jci2,ici1,ici2,1,kz+1,'cmicro:pfplsx')
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
  end subroutine allocate_mod_cloud_s1

  subroutine init_cloud_s1
    use mod_atm_interface
    use mod_runparams , only : vfqr , vfqi , vfqs , nssopt
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
    call assignpnt(ccn,pccn)

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

    ! Sanity check
    if ( nssopt < 0 .or. nssopt > 3 ) then
      call fatal(__FILE__,__LINE__,'NSSOPT IN CLOUD MUST BE IN RANGE 0-3')
    end if
  end subroutine init_cloud_s1

  subroutine microphys
    implicit none
    integer(ik4) :: i , j , k , n , m , jn , jo
    integer(ik4) :: ii , jj , kk , ll , imax
    logical :: lactiv
    real(rkx) :: rexplicit
    real(rkx) :: fac , faci , facw , corr , koop , gdp
    real(rkx) :: alfaw , phases , qice , zdelta , tmpl , &
                 tmpi , tnew , qe , rain , preclr , arg
    real(rkx) :: aamax , dum , xsum , swap
    ! local real variables for autoconversion rate constants
    real(rkx) :: alpha1 ! coefficient autoconversion cold cloud
    real(rkx) :: tmpa
    real(rkx) :: cfpr
    real(rkx) :: xlcrit
    real(rkx) :: precip
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
    ! real(rkx) :: gdph_r
    ! constants for deposition process
    real(rkx) :: vpice , vpliq , xadd , xbdd , cvds , &
                 qice0 , qinew , infactor , rainaut , snowaut
    ! constants for condensation and turbulent mixing erosion of clouds
    real(rkx) :: dpmxdt , wtot , dtdiab , dtforc , &
                 qp , qsat , cond1 , levap , leros
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'microphys'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    if ( .not. fscheme ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if

    oneodt = d_one/dt

    ! Set the default 1.e-14_rkx = d_zero
    ! Define the inizial array qx0
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          do n = 1 , nqx
            if ( qxx(j,i,k,n) > minqx ) then
              qx0(n,j,i,k) = qxx(j,i,k,n)
            else
              qx0(n,j,i,k) = minqx
            end if
          end do
        end do
      end do
    end do

    ! Detrainment : [xqdetr] = kg/(m^2*s)
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( qdetr(j,i,k) > minqx ) then
            xqdetr(j,i,k) = qdetr(j,i,k)
          else
            xqdetr(j,i,k) = d_zero
          end if
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
          qliq(j,i,k) = phase(t(j,i,k))
        end do
      end do
    end do

    ! Reset variables
    covptot(:,:)    = d_zero
    covpclr(:,:)    = d_zero
    cldtopdist(:,:) = d_zero
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
                tmpl = qx0(n,j,i,k)+dt*(qxtendc(n,j,i,k)-tenkeep(n,j,i,k))
                tmpi = d_zero
              else if ( iphase(n) == 2 ) then
                tmpi = qx0(n,j,i,k)+dt*(qxtendc(n,j,i,k)-tenkeep(n,j,i,k))
                tmpl = d_zero
              end if
              sumq0(j,i,k) = sumq0(j,i,k) + &
                (qx0(n,j,i,k)+(qxtendc(n,j,i,k)-tenkeep(n,j,i,k))*dt)* &
                (pfs(j,i,k+1)-pfs(j,i,k))*regrav
              tnew = tnew - wlhvocp*tmpl - wlhsocp*tmpi
              sumq0(j,i,k) = sumq0(j,i,k) + &
                (tmpl+tmpi)*(pfs(j,i,k+1)-pfs(j,i,k))*regrav    !(kg/m^2)
            end do
            ! Detrained water treated here
            qe = xqdetr(j,i,k)*dt*egrav/(pfs(j,i,k+1)-pfs(j,i,k)) ! 1 ?
            if ( qe > minqq ) then
              sumq0(j,i,k) = sumq0(j,i,k)+xqdetr(j,i,k)*dt
              alfaw = qliq(j,i,k)
              tnew = tnew-(wlhvocp*alfaw+wlhsocp*(d_one-alfaw))*qe
            end if
            sumh0(j,i,k) = sumh0(j,i,k)+(pfs(j,i,k+1)-pfs(j,i,k))*tnew
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
          qlt(j,i,k) = qx0(iqql,j,i,k)+qx0(iqqi,j,i,k)
          if ( qlt(j,i,k) > minqq ) then
            qliqfrac(j,i,k) = qx0(iqql,j,i,k)/qlt(j,i,k)
            qicefrac(j,i,k) = d_one-qliqfrac(j,i,k)
          else
            qliqfrac(j,i,k) = d_zero
            qicefrac(j,i,k) = d_zero
          end if
        end do
      end do
    end do

    !----------------------------------------------------------------------
    !                       START OF VERTICAL LOOP
    !----------------------------------------------------------------------
    !
    ! Loop over levels
    !
    do k = 1 , kz

      !---------------------------------
      ! First guess microphysics
      !---------------------------------
      do i = ici1 , ici2
        do j = jci1 , jci2
          do n = 1 , nqx
            qxfg(n,j,i) = qx0(n,j,i,k)
          end do
        end do
      end do

      do i = ici1 , ici2
        do j = jci1 , jci2
          fccfg(j,i,k) = min(hicld,max(pfcc(j,i,k),lowcld))
        end do
      end do

      ! Reset matrix so missing pathways are set
      solqb(:,:,:,:)  = d_zero
      solqa(:,:,:,:)  = d_zero
      fallsrce(:,:,:) = d_zero
      fallsink(:,:,:) = d_zero
      convsrce(:,:,:) = d_zero
      convsink(:,:,:) = d_zero
      ratio(:,:,:)    = d_zero
      qicetot(:,:)    = d_zero
      supsat(:,:)     = d_zero
      subsat(:,:)     = d_zero
      licld(:,:)      = d_zero
      ldefr(:,:)      = d_zero

      ! Set j,i arrays to zero
      qpretot(:,:) = d_zero
      ! xqdetr2(:,:,:) = d_zero

      ! Set stats variables to zero
      if ( stats ) then
        statssupw(:,:,:) = d_zero
        statssupc(:,:,:) = d_zero
      end if

      ! Derived variables needed
      do i = ici1, ici2
        do j = jci1, jci2
          dp(j,i) = pfs(j,i,k+1) - pfs(j,i,k)      ! dp
          gdp = egrav/dp(j,i)                      ! g/dp  =(1/m)
          dtgdp(j,i) = dt*gdp                      ! (dt*g)/dp =(dt/m)
          rdtgdp(j,i) = dp(j,i)*(d_one/(dt*egrav)) ! dp/(gdt)=m/dt  [Kg/m2/s]
          !------------------------------------
          ! calculate dqs/dT
          !------------------------------------
          ! liquid
          facw     = c5les/((t(j,i,k) - c4les)**2)
          corr      = d_one/(d_one - ep1*eeliqt(j,i,k))
          dqsliqdt = facw*corr*qsliq(j,i,k)
          corqsliq(j,i) = d_one + wlhvocp*dqsliqdt
          ! ice
          faci     = c5ies/((t(j,i,k) - c4ies)**2)
          corr      = d_one/(d_one - ep1*eew(j,i,k))
          dqsicedt = faci*corr*qsice(j,i,k)
          corqsice(j,i) = d_one + wlhsocp*dqsicedt
          ! diagnostic mixed
          alfaw    = qliq(j,i,k)
          fac      = alfaw*facw + (d_one - alfaw)*faci
          corr      = d_one/(d_one - ep1*eewmt(j,i,k))
          dqsmixdt = fac*corr*qsmix(j,i,k)
          corqsmix(j,i) = d_one/(d_one + eldcpm(t(j,i,k))*dqsmixdt)
          !--------------------------------
          ! evaporation/sublimation limits
          !--------------------------------
          evaplimmix(j,i) = max((qsmix(j,i,k) - qx0(iqqv,j,i,k)) * &
                                  corqsmix(j,i),d_zero)
          !--------------------------------
          ! in-cloud consensate amount
          !--------------------------------
          tmpa = d_one/fccfg(j,i,k)
          liqcld(j,i) = qx0(iqql,j,i,k)*tmpa
          icecld(j,i) = qx0(iqqi,j,i,k)*tmpa
          licld(j,i)  = liqcld(j,i) + icecld(j,i)
        end do
      end do
      !------------------------------------------------
      ! Evaporate very small amounts of liquid and ice
      !------------------------------------------------
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( qx0(iqql,j,i,k) < minqq ) then
            solqa(iqqv,iqql,j,i) =  qx0(iqql,j,i,k)
            solqa(iqql,iqqv,j,i) = -qx0(iqql,j,i,k)
            qxfg(iqql,j,i) = qxfg(iqql,j,i) - qx0(iqql,j,i,k)
          end if
          if ( qx0(iqqi,j,i,k) < minqq ) then
            solqa(iqqv,iqqi,j,i) =  qx0(iqqi,j,i,k)
            solqa(iqqi,iqqv,j,i) = -qx0(iqqi,j,i,k)
            qxfg(iqqi,j,i) = qxfg(iqqi,j,i) - qx0(iqqi,j,i,k)
          end if
        end do
      end do

      !---------------------------------------------------------------------
      !  ICE SUPERSATURATION ADJUSTMENT
      !---------------------------------------------------------------------
      ! Note that the supersaturation adjustment is made with respect to
      ! liquid saturation:  when T > 0C
      ! ice saturation:     when T < 0C
      !                     with an adjustment made to allow for ice
      !                     supersaturation in the clear sky
      ! Note also that the KOOP factor automatically clips the supersaturation
      ! to a maximum set by the liquid water saturation mixing ratio
      ! important for temperatures near to but below 0C
      !qv_max = qs * (fcc + (1-fcc) *RH_homo ) if T < 0C
      !qv_max = qs                             if T > 0C
      !-----------------------------------------------------------------------
      do i = ici1 , ici2
        do j = jci1 , jci2
          !-----------------------------------
          ! Supersaturation limit (from Koop)
          !-----------------------------------
          if ( nssopt == 0 )  then
            fac  = d_one
          else
            if ( t(j,i,k) >= tzero ) then
              fac  = d_one
            else
              koop = fkoop(t(j,i,k),st_eeliq(j,i,k),st_eeice(j,i,k))
              fac  = fccfg(j,i,k) + koop*(d_one-fccfg(j,i,k))
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
          supsat(j,i) = max((qx0(iqqv,j,i,k)-fac*qsmix(j,i,k)) / &
                              corqsliq(j,i),d_zero)
          ! e < esi, because for e > esi ice still present
          subsat(j,i) = min((qx0(iqqv,j,i,k)-fac*qsmix(j,i,k))/  &
                              corqsliq(j,i),d_zero)

          if ( supsat(j,i) > dlowval ) then
            if ( t(j,i,k) > thomo ) then
              ! turn supersaturation into liquid water
              solqa(iqql,iqqv,j,i) = solqa(iqql,iqqv,j,i)+supsat(j,i)
              solqa(iqqv,iqql,j,i) = solqa(iqqv,iqql,j,i)-supsat(j,i)
              qxfg(iqql,j,i) = qxfg(iqql,j,i)+supsat(j,i)
            else
              ! turn supersaturation into ice water
              solqa(iqqi,iqqv,j,i) = solqa(iqqi,iqqv,j,i)+supsat(j,i)
              solqa(iqqv,iqqi,j,i) = solqa(iqqv,iqqi,j,i)-supsat(j,i)
              qxfg(iqqi,j,i) = qxfg(iqqi,j,i)+supsat(j,i)
            end if
          else
            if ( subsat(j,i) < d_zero .and. &
                 fccfg(j,i,k) < zerocf .and. &
                 qlt(j,i,k) > minqx ) then
              ! turn subsaturation into vapor, where there is no cloud
              excess = qlt(j,i,k)+subsat(j,i)
              if ( excess < d_zero ) then
                if ( t(j,i,k) > thomo ) then
                  evapl = max(-qlt(j,i,k),-evaplimmix(j,i))!*oneodt
                  solqa(iqqv,iqql,j,i) = solqa(iqqv,iqql,j,i)-evapl
                  solqa(iqql,iqqv,j,i) = solqa(iqql,iqqv,j,i)+evapl
                  qxfg(iqql,j,i) = qxfg(iqql,j,i)+evapl
                else
                  evapi = max(-qlt(j,i,k),-evaplimmix(j,i))!*oneodt
                  ! turn subsaturation into vapour
                  solqa(iqqv,iqqi,j,i) = solqa(iqqv,iqqi,j,i)-evapi
                  solqa(iqqi,iqqv,j,i) = solqa(iqqi,iqqv,j,i)+evapi
                  qxfg(iqqi,j,i) = qxfg(iqqi,j,i)+evapi
                end if
              end if
            end if
          end if
        end do
      end do

      if ( stats ) then
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( supsat(j,i) > dlowval ) then
              if ( t(j,i,k) > thomo ) then
                statssupw(j,i,k) = supsat(j,i)
              else
                statssupc(j,i,k) = supsat(j,i)
              end if
            else
              if ( subsat(j,i) < d_zero .and. &
                   fccfg(j,i,k) < zerocf .and. &
                   qlt(j,i,k) > minqx ) then
                excess = qlt(j,i,k)+subsat(j,i)
                if ( excess < d_zero ) then
                  if ( t(j,i,k) > thomo ) then
                    evapl = max(-qlt(j,i,k),-evaplimmix(j,i))!*oneodt
                    statssupw(j,i,k) =  statssupw(j,i,k) + evapl
                  else
                    evapi = max(-qlt(j,i,k),-evaplimmix(j,i))!*oneodt
                    statssupc(j,i,k) =  statssupc(j,i,k) - evapi
                  end if
                end if
              end if
            end if
          end do
        end do
      end if

      ! call addpath_real(iqql,iqqv,supsatl,solqa,solqb,d_zero,qxfg)
      ! call addpath_real(iqqi,iqqv,supsati,solqa,solqb,d_zero,qxfg)
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
      !chng(:,:) = d_zero
      if ( k < kz .and. k >= 1 ) then
        do i = ici1 , ici2
          do j = jci1 , jci2
            xqdetr(j,i,k) = xqdetr(j,i,k)*dtgdp(j,i)  !kg/kg
            if ( xqdetr(j,i,k) > minqq ) then
              !qice = 1 if T < 250, qice = 0 if T > 273
              alfaw = qliq(j,i,k)
              qice   = d_one-alfaw
              convsrce(iqql,j,i) = alfaw*xqdetr(j,i,k)
              convsrce(iqqi,j,i) = qice*xqdetr(j,i,k)
              solqa(iqql,iqql,j,i) = solqa(iqql,iqql,j,i) + &
                                      convsrce(iqql,j,i)
              solqa(iqqi,iqqi,j,i) = solqa(iqqi,iqqi,j,i) + &
                                      convsrce(iqqi,j,i)
              qxfg(iqql,j,i) = qxfg(iqql,j,i)+convsrce(iqql,j,i)
              qxfg(iqqi,j,i) = qxfg(iqqi,j,i)+convsrce(iqqi,j,i)

              ! adjustment of the cloud fraction, that increases when qdetr
              ! adds liquid to the scheme
              ! qe = (qx0(iqqv,j,i,k)-fccfg(j,i,k)*qsmix(j,i,k)) / &
              !        qsice(j,i,k)) / (d_one-fccfg(j,i,k))
              ! qe = max(d_zero,min(qe,qsmix(j,i,k))) !qsice(j,i,k)))
              ! define the new cloud
              ! fccfg(j,i,k) = (relh(j,i,k)**0.25_rkx)* &
              !                (d_one-dexp((-100.0_rkx*(qxfg(iqql,j,i) + &
              !                 qxfg(iqqi,j,i))/qx0(iqqi,j,i,k)
              !                sqrt((d_one-relh(j,i,k))*qsice(j,i,k)))))
              ! fccfg(j,i,k) = dmin1(dmax1(fccfg(j,i,k),0.01_rkx),0.99_rkx)
              ! chng(j,i) = min(xqdetr(j,i,k),(fccfg(j,i,k) - &
              !                            fccfg(j,i,k))*(qsice(j,i,k)-qe))
              ! chng(j,i) = max(chng(j,i),d_zero)
              ! xqdetr2(j,i,k) = xqdetr2(j,i,k) - chng(j,i)
              ! if (xqdetr2(j,i,k) > d_zero) then
              !   fccfg(j,i,k) = fccfg(j,i,k)
              ! else
              !   fccfg(j,i,k) = pfcc(j,i,k)
              ! end if
              ! solqa(iqqv,iqql,j,i) = solqa(iqqv,iqql,j,i) + chng(j,i)
              ! solqa(iqql,iqqv,j,i) = solqa(iqql,iqqv,j,i) - chng(j,i)
              ! qxfg(iqql,j,i) = qxfg(iqql,j,i) - chng(j,i)
            end if
          end do
        end do
      end if

      if ( stats ) then
        if ( k < kz .and. k >= 1 ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              statsdetrw(j,i,k) = convsrce(iqql,j,i)
              statsdetrc(j,i,k) = convsrce(iqqi,j,i)
            end do
          end do
        end if
      end if


      !------------------------------------------------------------------
      ! Turn on/off microphysics
      !------------------------------------------------------------------

      if ( lmicro ) then
        !---------------------------------------
        ! EROSION OF CLOUDS BY TURBULENT MIXING
        !--------------------------------------
        ! rcldiff  : Diffusion coefficient for evaporation by turbulent
        ! mixing (IBID., EQU. 30) rcldiff = 3.0e-6_rkx
        do i = ici1 , ici2
          do j = jci1 , jci2
            ldifdt(j,i) = rcldiff*dt
            if ( xqdetr(j,i,k) > d_zero ) ldifdt(j,i) = 5.0_rkx*ldifdt(j,i)
          end do
        end do
        !Increase by factor of 5 for convective points
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( qlt(j,i,k) > minqx ) then
              leros = fccfg(j,i,k) * ldifdt(j,i) * &
                      max(qsmix(j,i,k)-qx0(iqqv,j,i,k),d_zero)
              leros = min(leros,evaplimmix(j,i))
              leros = min(leros,qlt(j,i,k))
              fac  = qliqfrac(j,i,k)*leros
              faci = qicefrac(j,i,k)*leros
              solqa(iqql,iqqv,j,i) = solqa(iqql,iqqv,j,i) - fac
              solqa(iqqv,iqql,j,i) = solqa(iqqv,iqql,j,i) + fac
              solqa(iqqi,iqqv,j,i) = solqa(iqqi,iqqv,j,i) - faci
              solqa(iqqv,iqqi,j,i) = solqa(iqqv,iqqi,j,i) + faci
              qxfg(iqql,j,i) = qxfg(iqql,j,i) - fac
              qxfg(iqqi,j,i) = qxfg(iqqi,j,i) - faci
            end if
          end do
        end do

        if ( stats ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              if ( qlt(j,i,k) > minqx ) then
                leros = fccfg(j,i,k)*ldifdt(j,i) * &
                        max(qsmix(j,i,k)-qx0(iqqv,j,i,k),d_zero)
                leros = min(leros,evaplimmix(j,i))
                leros = min(leros,qlt(j,i,k))
                statserosw(j,i,k) = qliqfrac(j,i,k)*leros
                statserosc(j,i,k) = qicefrac(j,i,k)*leros
              end if
            end do
          end do
        end if

        !----------------------------------------------------------------------
        ! CONDENSATION/EVAPORATION DUE TO DQSAT/DT
        !----------------------------------------------------------------------
        !  calculate dqs/dt
        !  Note: For the separate prognostic Qi and Ql, one would ideally use
        !  Qsat/DT wrt liquid/Koop here, since the physics is that new clouds
        !  forms by liquid droplets [liq] or when aqueous aerosols [Koop] form.
        !  These would then instant. freeze if T<-38C or lead to ice growth
        !  by deposition in warmer mixed phase clouds.  However, since we do
        !  not have a separate prognostic equation for in-cloud humidity or a
        !  statistical scheme approach in place, the depositional growth of ice
        !  in the mixed phase can not be modelled and we resort to
        !  supersaturation
        !  wrt ice instanteously converting to ice over one timestep
        !  (see Tompkins et al. QJRMS 2007 for details)
        !  Thus for the initial implementation the diagnostic mixed phase is
        !  retained for the moment, and the level of approximation noted.
        !----------------------------------------------------------------------
        do i = ici1 , ici2
          do j = jci1 , jci2
            dtdp   = rovcp*t(j,i,k)/phs(j,i,k)
            dpmxdt = dp(j,i)*oneodt
            wtot   = pverv(j,i,k)
            wtot   = min(dpmxdt,max(-dpmxdt,wtot))
            dtdiab = min(dpmxdt*dtdp,max(-dpmxdt*dtdp,radheatrt(j,i,k)))*dt + &
                          wlhfocp*ldefr(j,i)
            ! ldefr = 0
            ! note: ldefr should be set to the difference between the mixed
            ! phase functions in the convection and cloud scheme, and
            ! for now we set it to zero and the functions are the same.
            ! In RegCM not all convection schemes provide such info.
            dtforc       = dtdp*wtot*dt + dtdiab
            qold(j,i)    = qsmix(j,i,k)
            told(j,i)    = t(j,i,k)
            tcond(j,i,k) = t(j,i,k)+dtforc
            tcond(j,i,k) = max(tcond(j,i,k),160.0_rkx)
          end do
        end do
        !this loop's goal is to produce dqs = qsmix - qold, where qsmix is
        !reduced because of the condensation. so that dqs is negative?
        do i = ici1 , ici2
          do j = jci1 , jci2
            qp = d_one/phs(j,i,k)
            phases = phase(tcond(j,i,k))
            ! saturation mixing ratio ws
            qsat = eewm(tcond(j,i,k),phases)*qp
            qsat = min(d_half,qsat)          ! ws < 0.5        WHY?
            corr  = d_one/(d_one-ep1*qsat)
            qsat = qsat*corr
            cond = (qsmix(j,i,k)-qsat) / &
                    (d_one + qsat*edem(tcond(j,i,k),phases))
            tcond(j,i,k) = tcond(j,i,k)+eldcpm(tcond(j,i,k))*cond

            phases = phase(tcond(j,i,k))
            qsmix(j,i,k) = qsmix(j,i,k)-cond
            qsat = eewm(tcond(j,i,k),phases)*qp
            qsat = min(d_half,qsat)
            corr = d_one/(d_one-ep1*qsat)
            qsat = qsat*corr
            cond1 = (qsmix(j,i,k)-qsat)  / &
                     (d_one + qsat*edem(tcond(j,i,k),phases))
            tcond(j,i,k) = tcond(j,i,k)+eldcpm(tcond(j,i,k))*cond1
            qsmix(j,i,k) = qsmix(j,i,k)-cond1
          end do
        end do
        do i = ici1 , ici2
          do j = jci1 , jci2
            dqs(j,i)     = qsmix(j,i,k)-qold(j,i)
            qsmix(j,i,k) = qold(j,i)
            tcond(j,i,k) = told(j,i)
          end do
        end do

        !----------------------------------------------------------------
        ! dqs(jl) > 0:  evaporation of clouds
        !----------------------------------------------------------------
        ! erosion term is linear in l
        ! changed to be uniform distribution in cloud region
        do i = ici1 , ici2
          do j = jci1 , jci2
            ! previous function based on delta distribution in cloud:
            if ( dqs(j,i) > d_zero ) then
              !levap = C*min( dqs/dt , (qi+ql)/C )
              levap = fccfg(j,i,k)*min(dqs(j,i),licld(j,i))
              levap = min(levap,evaplimmix(j,i))
              levap = min(levap,max(qsmix(j,i,k)-qx0(iqqv,j,i,k),d_zero))
              fac = qliqfrac(j,i,k)*levap
              faci = qicefrac(j,i,k)*levap
              solqa(iqqv,iqql,j,i) = solqa(iqqv,iqql,j,i) + fac
              solqa(iqql,iqqv,j,i) = solqa(iqql,iqqv,j,i) - fac
              solqa(iqqv,iqqi,j,i) = solqa(iqqv,iqqi,j,i) + faci
              solqa(iqqi,iqqv,j,i) = solqa(iqqi,iqqv,j,i) - faci
              qxfg(iqql,j,i) = qxfg(iqql,j,i) - fac
              qxfg(iqqi,j,i) = qxfg(iqqi,j,i) - faci
            end if
          end do
        end do
        if ( stats ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              if ( dqs(j,i) > d_zero ) then
                levap = fccfg(j,i,k)*min(dqs(j,i),licld(j,i))
                levap = min(levap,evaplimmix(j,i))
                levap = min(levap,max(qsmix(j,i,k)-qx0(iqqv,j,i,k),d_zero))
                statsevapw(j,i,k) = qliqfrac(j,i,k)*levap
                statsevapc(j,i,k) = qicefrac(j,i,k)*levap
              end if
            end do
          end do
        end if

        !----------------------------------------------------------------------
        ! dqs(j,i) < 0: formation of clouds
        !----------------------------------------------------------------------
        ! (1) increase of cloud water in existing clouds
        chng(:,:) = d_zero

        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( fccfg(j,i,k) > zerocf .and. dqs(j,i) <= -minqq ) then
              ! new limiter
              chng(j,i) = max(-dqs(j,i),d_zero)
              ! old limiter
              !  (significantly improves upper tropospheric humidity rms)
              phases = qliq(j,i,k)
              if ( fccfg(j,i,k) > onecf ) then
                corr = d_one/(d_one-ep1*qsmix(j,i,k))
                cdmax = (qx0(iqqv,j,i,k)-qsmix(j,i,k)) / &
                         (d_one + corr*qsmix(j,i,k)*edem(t(j,i,k),phases))
              else
                cdmax = (qx0(iqqv,j,i,k)-fccfg(j,i,k)*qsmix(j,i,k)) / &
                        fccfg(j,i,k)
              end if
              chng(j,i) = max(min(chng(j,i),cdmax),d_zero)
              ! end old limiter
              chng(j,i) = fccfg(j,i,k)*chng(j,i)
              if ( chng(j,i) < minqq ) chng(j,i) = d_zero
              !----------------------------------------------------------------
              ! All increase goes into liquid unless so cold cloud
              ! homogeneously freezes
              ! include new liquid formation in first guess value, otherwise
              ! liquid remains at cold temperatures until next timestep.
              !----------------------------------------------------------------
              if ( t(j,i,k) > thomo ) then
                solqa(iqql,iqqv,j,i) = solqa(iqql,iqqv,j,i) + chng(j,i)
                solqa(iqqv,iqql,j,i) = solqa(iqqv,iqql,j,i) - chng(j,i)
                qxfg(iqql,j,i) = qxfg(iqql,j,i) + chng(j,i)
              else
                solqa(iqqi,iqqv,j,i) = solqa(iqqi,iqqv,j,i) + chng(j,i)
                solqa(iqqv,iqqi,j,i) = solqa(iqqv,iqqi,j,i) - chng(j,i)
                qxfg(iqqi,j,i) = qxfg(iqqi,j,i) + chng(j,i)
              end if
            end if
          end do
        end do

        if ( stats ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              if ( t(j,i,k) > thomo ) then
                statscond1w(j,i,k) = chng(j,i)
              else
                statscond1c(j,i,k) = chng(j,i)
              end if
            end do
          end do
        end if

        chng(:,:) = d_zero
        qexc(:,:) = d_zero

        if ( nssopt == 0 .or. nssopt == 1 ) then ! Tompkins
          do i = ici1 , ici2
            do j = jci1 , jci2
              if ( dqs(j,i) <= -minqq .and. fccfg(j,i,k) < onecf ) then
                qexc(j,i) = max((qx0(iqqv,j,i,k) -  &
                        fccfg(j,i,k)*qsmix(j,i,k)) / &
                        (d_one-fccfg(j,i,k)),d_zero)
              end if
            end do
          end do
        else if ( nssopt == 2 ) then ! Lohmann and Karcher
          do i = ici1 , ici2
            do j = jci1 , jci2
              if ( dqs(j,i) <= -minqq .and. fccfg(j,i,k) < onecf ) then
                qexc(j,i) = qx0(iqqv,j,i,k)
              end if
            end do
          end do
        else if ( nssopt == 3 ) then ! Gierens
          do i = ici1 , ici2
            do j = jci1 , jci2
              if ( dqs(j,i) <= -minqq .and. fccfg(j,i,k) < onecf ) then
                qexc(j,i) = qx0(iqqv,j,i,k)/qlt(j,i,k)
              end if
            end do
          end do
        end if

        ! (2) generation of new clouds (dc/dt>0)
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( dqs(j,i) <= -minqq .and. fccfg(j,i,k) < onecf ) then
              !---------------------------
              ! critical relative humidity
              !---------------------------
              ! *RAMID*   REAL    BASE VALUE FOR CALCULATION OF RELATIVE
              !                   HUMIDITY THRESHOLD FOR ONSET OF STRATIFORM
              !                   CONDENSATION (TIEDTKE, 1993, EQUATION 24)
              rhc = ramid !=0.8
              zsig = phs(j,i,k)/pfs(j,i,kz+1)
              ! increase RHcrit to 1.0 towards the surface (sigma>0.8)
              if ( zsig > 0.8_rkx ) then
                rhc = ramid+(d_one-ramid)*((zsig-0.8_rkx)/0.2_rkx)**2
              end if
              !---------------------------
              ! supersaturation options
              !---------------------------
              if ( t(j,i,k) >= tzero .or. nssopt == 0 ) then
                ! no ice supersaturation allowed
                fac = d_one
              else
                ! ice supersaturation
                fac = fkoop(t(j,i,k),st_eeliq(j,i,k),st_eeice(j,i,k))
              end if
              if ( qexc(j,i) >= rhc*qsmix(j,i,k)*fac .and. &
                   qexc(j,i) < qsmix(j,i,k)*fac ) then
                ! note: not **2 on 1-a term if qe is used.
                ! added correction term fac to numerator 15/03/2010
                acond = -(d_one-fccfg(j,i,k))*fac*dqs(j,i) / &
                          max(d_two*(fac*qsmix(j,i,k)-qexc(j,i)),dlowval)
                acond = min(acond,d_one-fccfg(j,i,k)) ! put the limiter back
                ! linear term:
                ! added correction term fac 15/03/2010
                chng(j,i) = -fac*dqs(j,i)*d_half*acond !mine linear
                ! new limiter formulation
                ! qsice(j,i,k)-qexc(j,i)) /
                tmpa = d_one-fccfg(j,i,k)
                zdl = d_two*(fac*qsmix(j,i,k)-qexc(j,i)) / tmpa
                ! added correction term fac 15/03/2010
                if ( fac*dqs(j,i) < -zdl ) then
                  ! qsice(j,i,k)+qx0(iqqv,j,i,k)
                  xlcondlim = (fccfg(j,i,k)-d_one)*fac*dqs(j,i) - &
                               fac*qsmix(j,i,k)+qx0(iqqv,j,i,k)
                  chng(j,i) = min(chng(j,i),xlcondlim)
                end if
                chng(j,i) = max(chng(j,i),d_zero)
                if ( chng(j,i) < minqq ) then
                  chng(j,i) = d_zero
                  acond     = d_zero
                end if
                !-------------------------------------------------------------
                ! all increase goes into liquid unless so cold cloud
                ! homogeneously freezes
                ! include new liquid formation in first guess value, otherwise
                ! liquid remains at cold temperatures until next timestep.
                !-------------------------------------------------------------
                if ( t(j,i,k) > thomo ) then
                  solqa(iqql,iqqv,j,i) = solqa(iqql,iqqv,j,i) + chng(j,i)
                  solqa(iqqv,iqql,j,i) = solqa(iqqv,iqql,j,i) - chng(j,i)
                  qxfg(iqql,j,i) = qxfg(iqql,j,i) + chng(j,i)
                else ! homogeneous freezing
                  solqa(iqqi,iqqv,j,i) = solqa(iqqi,iqqv,j,i) + chng(j,i)
                  solqa(iqqv,iqqi,j,i) = solqa(iqqv,iqqi,j,i) - chng(j,i)
                  qxfg(iqqi,j,i) = qxfg(iqqi,j,i) + chng(j,i)
                end if
              end if
            end if
          end do
        end do

        if ( stats ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              if ( t(j,i,k) > thomo ) then
                statscond2w(j,i,k) = chng(j,i)
              else
                statscond2c(j,i,k) = chng(j,i)
              end if
            end do
          end do
        end if

        !----------------------------------------------------------------------
        ! DEPOSITION: Growth of ice by vapour deposition
        !----------------------------------------------------------------------
        ! Following Rotstayn et al. 2001:
        ! does not use the ice nuclei number from cloudaer.F90
        ! but rather a simple Meyers et al. 1992 form based on the
        ! supersaturation and assuming clouds are saturated with
        ! respect to liquid water (well mixed), (or Koop adjustment)
        ! Growth considered as sink of liquid water if present so
        ! Bergeron-Findeisen adjustment in autoconversion term no longer needed
        !----------------------------------------------------------------------

        chng(:,:) = d_zero

        do i = ici1 , ici2
          do j = jci1 , jci2
            !--------------------------------------------------------------
            ! Calculate distance from cloud top
            ! defined by cloudy layer below a layer with cloud frac <0.01
            ! DZ = DP(JL)/(RHO(JL)*RG)
            !--------------------------------------------------------------
            if ( k > 1 ) then
              if ( fccfg(j,i,k-1) < cldtopcf .and. &
                   fccfg(j,i,k)  >= cldtopcf ) then
                cldtopdist(j,i) = d_zero
              else
                cldtopdist(j,i) = cldtopdist(j,i) + dp(j,i)/(rho(j,i,k)*egrav)
              end if
            end if
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
            if ( t(j,i,k) < tzero .and. qxfg(iqql,j,i) > minqq ) then
              vpice = st_eeice(j,i,k) !saturation vapor pressure wrt ice
              vpliq = st_eeliq(j,i,k) !saturation vapor pressure wrt liq
              ! Meyers et al 1992
              icenuclei(j,i) = d_1000 * &
                exp(12.96_rkx*((vpliq-vpice)/vpice)-0.639_rkx)
              !------------------------------------------------
              !   2.4e-2 is conductivity of air
              !   8.87 = 700**1/3 = density of ice to the third
              !------------------------------------------------
              xadd  = wlhs*(wlhs/(rwat*t(j,i,k))-d_one)/(2.4e-2_rkx*t(j,i,k))
              xbdd  = rwat*t(j,i,k)*phs(j,i,k)/(2.21_rkx*vpice)
              cvds = 7.8_rkx*(icenuclei(j,i)/rho(j,i,k))** &
                       0.666_rkx*(vpliq-vpice)/(8.87_rkx*(xadd+xbdd)*vpice)
              !-----------------------------------------------------
              ! iceinit = 1.e-12 is initial mass of ice particle
              !-----------------------------------------------------
              qice0 = max(icecld(j,i), icenuclei(j,i)*iceinit/rho(j,i,k))
              !------------------
              ! new value of ice condensate amount ( Rotstayn et al. (2000) )
              !------------------
              qice0 = max(qice0,d_zero)
              cvds = max(cvds,d_zero)
              qinew = (0.666_rkx*cvds*dt+qice0**0.666_rkx)**1.5_rkx
              !---------------------------
              ! grid-mean deposition rate:
              !---------------------------
              chng(j,i) = max(fccfg(j,i,k)*(qinew-qice0),d_zero)*2.0_rkx
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
              chng(j,i) = min(chng(j,i),qxfg(iqql,j,i))
              !---------------------------------------------------------------
              ! at top of cloud, reduce deposition rate near cloud top to
              ! account for small scale turbulent processes, limited ice
              ! nucleation and ice fallout
              !---------------------------------------------------------------
              ! Fraction of deposition rate in cloud top layer
              ! depliqrefrate  = d_r10
              ! Depth of supercooled liquid water layer (m)
              ! depliqrefdepth = 500.0_rkx
              infactor = min(icenuclei(j,i)/15000.0_rkx, d_one)
              chng(j,i) = chng(j,i)*min(infactor + (d_one-infactor)* &
                     (depliqrefrate+cldtopdist(j,i)/depliqrefdepth),d_one)
              !--------------
              ! add to matrix
              !--------------
              solqa(iqqi,iqql,j,i) = solqa(iqqi,iqql,j,i) + chng(j,i)
              solqa(iqql,iqqi,j,i) = solqa(iqql,iqqi,j,i) - chng(j,i)
              qxfg(iqql,j,i) = qxfg(iqql,j,i) - chng(j,i)
              qxfg(iqqi,j,i) = qxfg(iqqi,j,i) + chng(j,i)
            end if
          end do
        end do

        if ( stats ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              statsdepos(j,i,k) = chng(j,i)
            end do
          end do
        end if

        do i = ici1 , ici2
          do j = jci1 , jci2
            tmpa = d_one/fccfg(j,i,k)
            liqcld(j,i) = qxfg(iqql,j,i)*tmpa
            icecld(j,i) = qxfg(iqqi,j,i)*tmpa
          end do
        end do

        !------------------------------------------------------------------
        !  SEDIMENTATION/FALLING OF *ALL* MICROPHYSICAL SPECIES
        !     now that rain and snow species are prognostic
        !     the precipitation flux can be defined directly level by level
        !     There is no vertical memory required from the flux variable
        !------------------------------------------------------------------
        do i = ici1 , ici2
          do j = jci1 , jci2
            do n = 1 , nqx
              if ( lfall(n) ) then
                ! Source from layer above
                if ( k > 1 ) then
                  fallsrce(n,j,i) = pfplsx(n,j,i,k)*dtgdp(j,i)
                  solqa(n,n,j,i) = solqa(n,n,j,i) + fallsrce(n,j,i)
                  qxfg(n,j,i) = qxfg(n,j,i) + fallsrce(n,j,i)
                  qpretot(j,i)= qpretot(j,i) + qxfg(n,j,i)
                end if
                ! Sink to next layer, constant fall speed
                fallsink(n,j,i) = dtgdp(j,i)*vqx(n)*rho(j,i,k)   !Kg/Kg
              end if  !lfall
            end do
          end do
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
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( qpretot(j,i) > dlowval ) then
              covptot(j,i) = d_one - ((d_one-covptot(j,i)) * &
                              (d_one - max(fccfg(j,i,k),fccfg(j,i,k-1))) / &
                              (d_one - min(fccfg(j,i,k-1),hicld)))
              covptot(j,i) = max(covptot(j,i),rcovpmin)   !rcovpmin = 0.1
              ! clear sky proportion
              covpclr(j,i) = max(d_zero,covptot(j,i)-fccfg(j,i,k))
            else
              covptot(j,i) = d_zero  ! no flux - reset cover
              covpclr(j,i) = d_zero  ! reset clear sky proportion
            end if
          end do
        end do
        !---------------------------------------------------------------
        !                         AUTOCONVERSION
        !---------------------------------------------------------------
        ! Warm clouds
        select case (iautoconv)
          case (1) ! Klein & Pincus (2000)
            do i = ici1 , ici2
              do j = jci1 , jci2
                if ( liqcld(j,i) > minqx ) then
                  solqb(iqql,iqqv,j,i) = d_zero
                  solqb(iqqr,iqql,j,i) = solqb(iqqr,iqql,j,i) + &
                    dt*auto_rate_klepi * (qx0(iqql,j,i,k)**(2.3_rkx))
                  solqa(iqqr,iqql,j,i) = d_zero
                end if
              end do
            end do
          case (2) ! Khairoutdinov and Kogan (2000)
            do i = ici1 , ici2
              do j = jci1 , jci2
                if ( liqcld(j,i) > minqx ) then
                  solqb(iqql,iqqv,j,i) = d_zero
                  solqb(iqqr,iqql,j,i) = solqb(iqqr,iqql,j,i) + &
                    dt*auto_rate_khair*(qx0(iqql,j,i,k)**(auto_expon_khair))
                end if
              end do
            end do
          case (3) ! Kessler(1969)
            do i = ici1 , ici2
              do j = jci1 , jci2
                if ( liqcld(j,i) > minqx ) then
                  solqb(iqql,iqqv,j,i) = d_zero
                  solqa(iqqr,iqql,j,i) = solqa(iqqr,iqql,j,i) - &
                               auto_rate_kessl*autocrit_kessl
                  solqa(iqql,iqqr,j,i) = solqa(iqql,iqqr,j,i) + &
                               auto_rate_kessl*autocrit_kessl
                  solqb(iqqr,iqql,j,i) = solqb(iqqr,iqql,j,i) + &
                               dt*auto_rate_kessl
                end if
              end do
            end do
          case (4) ! Sundqvist
            do i = ici1 , ici2
              do j = jci1 , jci2
                if ( liqcld(j,i) > minqx ) then
                  solqb(iqql,iqqv,j,i) = d_zero
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
                  if ( ichem == 1 .and. &
                       iaerosol == 1 .and. &
                       iindirect == 2 ) then
                    if ( pccn(j,i,k) > 0._rkx ) then
                      ! aerosol second indirect effect on autoconversion
                      ! threshold, rcrit is a critical cloud radius for cloud
                      ! water undergoing autoconversion
                      ! pccn = number of ccn /m3
                      xlcrit = pccn(j,i,k)*(4.0_rkx/3.0_rkx)*mathpi * &
                               ((rcrit*1e-6_rkx)**3)*rhoh2o
                    endif
                  endif
                  !-----------------------------------------------------------
                  ! parameters for cloud collection by rain and snow.
                  ! note that with new prognostic variable it is now possible
                  ! to replace this with an explicit collection
                  ! parametrization
                  !-----------------------------------------------------------
                  precip = (pfplsx(iqqs,j,i,k)+pfplsx(iqqr,j,i,k)) / &
                             max(dlowval,covptot(j,i))
                  cfpr = d_one + rprc1*sqrt(max(precip,d_zero))
                  alpha1 = alpha1*cfpr
                  xlcrit = xlcrit/max(cfpr,dlowval)
                  ! security for exp for some compilers
                  if ( (liqcld(j,i)/xlcrit)**2 < 25.0_rkx ) then
                    rainaut = alpha1*(d_one - exp(-(liqcld(j,i)/xlcrit)**2))
                  else
                    rainaut = alpha1
                  end if

                  !-----------------------
                  ! rain freezes instantly
                  !-----------------------
                  if ( t(j,i,k) <= tzero ) then
                    solqb(iqqs,iqql,j,i) = solqb(iqqs,iqql,j,i)+rainaut
                  else
                    solqb(iqqr,iqql,j,i) = solqb(iqqr,iqql,j,i)+rainaut
                  end if
                end if
              end do
            end do
        end select

        ! Cold clouds
        ! Snow Autoconversion rate follow Lin et al. 1983
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( t(j,i,k) <= tzero ) then
              if ( icecld(j,i) > minqx ) then
                alpha1 = dt*1.0e-3_rkx*exp(0.025*(t(j,i,k)-tzero))
                xlcrit = rlcritsnow
                arg = (icecld(j,i)/xlcrit)**2
                if ( arg < 25.0_rkx ) then
                  snowaut = alpha1 * (d_one - exp(-arg))
                else
                  snowaut = alpha1
                end if
                solqb(iqqs,iqqi,j,i) = solqb(iqqs,iqqi,j,i) + snowaut
              end if
            end if
          end do
        end do

        !---------------------------------------------------------------
        !                         MELTING
        !---------------------------------------------------------------
        ! The melting of ice and snow are treated explicitly.
        ! First water and ice saturation are found
        !---------------------------------------------
        ! ice saturation T < 273K
        ! liquid water saturation for T > 273K
        !---------------------------------------------

        do i = ici1, ici2
          do j = jci1, jci2
            do n = 1 , nqx
              if ( iphase(n) == 2 ) then
                qicetot(j,i) = qicetot(j,i) + qxfg(n,j,i)
              end if
            end do
          end do
        end do

        chngmax(:,:) = d_zero

        do i = ici1, ici2
          do j = jci1, jci2
            if ( qicetot(j,i) > minqx .and. t(j,i,k) > tzero ) then
              ! Calculate subsaturation
              ! qsice(j,i,k)-qx0(iqqv,j,i,k),d_zero)
              subsat(j,i) = max(qsmix(j,i,k)-qx0(iqqv,j,i,k),d_zero)
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
              tdiff = t(j,i,k)-tzero ! - subsat * &
              !    (tw1+tw2*(phs(j,i,k)-tw3)-tw4*(t(j,i,k)-tw5))
              ! Ensure CONS1 is positive so that MELTMAX = 0 if TDMTW0 < 0
              cons1 = d_one ! abs(dt*(d_one + d_half*tdiff)/rtaumel)
              chngmax(j,i) = max(tdiff*cons1*rldcp,d_zero)
            end if
          end do
        end do

        ! Loop over frozen hydrometeors (iphase == 2 (ice, snow))

        chng(:,:) = d_zero
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( chngmax(j,i) > dlowval .and. qicetot(j,i) > minqx ) then
              do n = 1, nqx
                if ( iphase(n) == 2 ) then
                  m = imelt(n) ! imelt(iqqi)=iqql, imelt(iqqs)=iqqr
                  if ( m < 0 ) cycle
                  phases = qxfg(n,j,i)/qicetot(j,i)
                  chng(j,i) = min(qxfg(n,j,i),phases*chngmax(j,i))
                  chng(j,i) = max(chng(j,i),d_zero)
                  ! n = iqqi,iqqs; m = iqql,iqqr
                  qxfg(n,j,i) =  qxfg(n,j,i) - chng(j,i)
                  qxfg(m,j,i) =  qxfg(m,j,i) + chng(j,i)
                  solqa(m,n,j,i) =  solqa(m,n,j,i) + chng(j,i)
                  solqa(n,m,j,i) =  solqa(n,m,j,i) - chng(j,i)
                end if
              end do
            end if
          end do
        end do

        if ( stats ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              statsmelt(j,i,k) = chng(j,i)
            end do
          end do
        end if

        !------------------------------------------------------------!
        !                         FREEZING                           !
        !------------------------------------------------------------!

        ! Freezing of rain.
        ! All rain freezes in a timestep if the temperature is below 0 C
        ! calculate sublimation latent heat

        do i = ici1 , ici2
          do j = jci1 , jci2
            chngmax(j,i) = max((tzero-t(j,i,k))*rldcp,d_zero)
          end do
        end do

        chng(:,:) = d_zero

        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( chngmax(j,i) > dlowval .and. qxfg(iqqr,j,i) > minqx ) then
              chng(j,i) = min(qxfg(iqqr,j,i),chngmax(j,i))
              chng(j,i) = max(chng(j,i),d_zero)
              solqa(iqqs,iqqr,j,i) = solqa(iqqs,iqqr,j,i) + chng(j,i)
              solqa(iqqr,iqqs,j,i) = solqa(iqqr,iqqs,j,i) - chng(j,i)
            end if
          end do
        end do

        if ( stats ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              statsfrz(j,i,k) = chng(j,i)
            end do
          end do
        end if

        ! Freezing of liquid

        do i = ici1 , ici2
          do j = jci1 , jci2
            chngmax(j,i) = max((thomo-t(j,i,k))*rldcp,d_zero)
          end do
        end do

        chng(:,:) = d_zero

        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( chngmax(j,i) > dlowval .and. qxfg(iqql,j,i) > minqx ) then
              chng(j,i) = min(qxfg(iqql,j,i),chngmax(j,i))
              chng(j,i) = max(chng(j,i),d_zero)
              solqa(iqqi,iqql,j,i) = solqa(iqqi,iqql,j,i) + chng(j,i)
              solqa(iqql,iqqi,j,i) = solqa(iqql,iqqi,j,i) - chng(j,i)
              qxfg(iqql,j,i) = qxfg(iqql,j,i) - chng(j,i)
              qxfg(iqqi,j,i) = qxfg(iqqi,j,i) + chng(j,i)
            end if
          end do
        end do

        if ( stats ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              statsfrz(j,i,k) = statsfrz(j,i,k) + chng(j,i)
            end do
          end do
        end if

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

        do i = ici1 , ici2
          do j = jci1 , jci2
           ! if ( iphase(n) == 2 ) then
           !   qpretot(j,i) = qpretot(j,i)+qxfg(n,j,i)
           qpretot(j,i) = qxfg(iqqs,j,i) + qxfg(iqqr,j,i)
           ! end if
          end do
        end do

        ! rain
        chng(:,:) = d_zero
        do i = ici1 , ici2
          do j = jci1 , jci2
            zrh = rprecrhmax + &
              (d_one-rprecrhmax)*covpclr(j,i)/(d_one-fccfg(j,i,k))
            zrh = min(max(zrh,rprecrhmax),d_one)
            ! This is a critical relative humidity that is used to limit
            ! moist environment to prevent the gridbox saturating when
            ! only part of the gridbox has evaporating precipitation
            qe = (qx0(iqqv,j,i,k)-fccfg(j,i,k)*qsliq(j,i,k)) / &
                      (d_one-fccfg(j,i,k))
            !---------------------------------------------
            ! humidity in moistest covpclr part of domain
            !---------------------------------------------
            qe = max(d_zero,min(qe,qsliq(j,i,k)))
            lactiv = covpclr(j,i) > dlowval .and. &
            !      qpretot(jl) > dlowval .and. &
                   qxfg(iqqr,j,i) > minqq .and. qe < zrh*qsliq(j,i,k)
            if ( lactiv ) then
              ! note: units of preclr and qpretot differ
              !       qpretot is a mixing ratio (hence "q" in name)
              !       preclr is a rain flux
              preclr = qpretot(j,i)*covpclr(j,i) / &
                          (max(dlowval,covptot(j,i))*dtgdp(j,i))
              !--------------------------------------
              ! actual microphysics formula in beta
              !--------------------------------------
              ! sensitivity test showed multiply rain evap rate by 0.5
              beta1 = sqrt(phs(j,i,k)/pfs(j,i,kz+1))/5.09e-3_rkx*preclr / &
                       max(covpclr(j,i),dlowval)

              if ( beta1 >= d_zero ) then
                 beta = egrav*rpecons*d_half*(beta1)**0.5777_rkx
                 denom = d_one + beta*dt*corqsliq(j,i)
                 dpr = covpclr(j,i) * beta * &
                   (qsliq(j,i,k)-qe)/denom*dp(j,i)*regrav
                 dpevap = dpr*dtgdp(j,i)
                !---------------------------------------------------------
                ! add evaporation term to explicit sink.
                ! this has to be explicit since if treated in the implicit
                ! term evaporation can not reduce rain to zero and model
                ! produces small amounts of rainfall everywhere.
                !---------------------------------------------------------

                ! evaporate rain
                chng(j,i) = min(dpevap,qxfg(iqqr,j,i))
                !-------------------------------------------------------------
                ! reduce the total precip coverage proportional to evaporation
                !-------------------------------------------------------------
                covptot(j,i) = covptot(j,i) - max(d_zero, &
                           (covptot(j,i)-fccfg(j,i,k))*dpevap/qpretot(j,i))
              else
                chng(j,i) = qxfg(iqqr,j,i)
              end if
              solqa(iqqv,iqqr,j,i) = solqa(iqqv,iqqr,j,i) + chng(j,i)
              solqa(iqqr,iqqv,j,i) = solqa(iqqr,iqqv,j,i) - chng(j,i)
              qxfg(iqqr,j,i)       = qxfg(iqqr,j,i) - chng(j,i)
            end if
          end do
        end do

        if ( stats ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              statsrainev(j,i,k) = chng(j,i)
            end do
          end do
        end if

        ! snow

        chng(:,:) = d_zero

        do i = ici1 , ici2
          do j = jci1 , jci2
            zrh = rprecrhmax + (d_one-rprecrhmax) * &
                   covpclr(j,i)/(d_one-fccfg(j,i,k))
            zrh = min(max(zrh,rprecrhmax),d_one)
            qe = (qx0(iqqv,j,i,k)-fccfg(j,i,k) * &
                   qsice(j,i,k))/(d_one-fccfg(j,i,k))
            !---------------------------------------------
            ! humidity in moistest covpclr part of domain
            !---------------------------------------------
            qe = max(d_zero,min(qe,qsice(j,i,k)))
            lactiv = covpclr(j,i) > dlowval .and. &
                   qxfg(iqqs,j,i) > minqq .and. &
                   qe < zrh*qsice(j,i,k)
            if ( lactiv ) then
              ! note: units of preclr and qpretot differ
              !       qpretot is a mixing ratio (hence "q" in name)
              !       preclr is a rain flux
              preclr = qpretot(j,i)*covpclr(j,i) / &
                    (max(dlowval,covptot(j,i))*dtgdp(j,i))
              !--------------------------------------
              ! actual microphysics formula in beta
              !--------------------------------------
               beta1 = sqrt(phs(j,i,k)/pfs(j,i,kz+1)) / &
                    5.09e-3_rkx*preclr/max(covpclr(j,i),dlowval)

              if ( beta1 >= d_zero ) then
                beta = egrav*rpecons*(beta1)**0.5777_rkx
                !rpecons = alpha1
                denom = d_one + beta*dt*corqsice(j,i)
                dpr = covpclr(j,i) * beta * &
                       (qsice(j,i,k)-qe)/denom*dp(j,i)*regrav
                dpevap = dpr*dtgdp(j,i)
                !---------------------------------------------------------
                ! add evaporation term to explicit sink.
                ! this has to be explicit since if treated in the implicit
                ! term evaporation can not reduce snow to zero and model
                ! produces small amounts of snowfall everywhere.
                !---------------------------------------------------------
                ! evaporate snow
                chng(j,i) = min(dpevap,qxfg(iqqs,j,i))
                chng(j,i) = max(chng(j,i),d_zero)
                !-------------------------------------------------------------
                ! reduce the total precip coverage proportional to evaporation
                !-------------------------------------------------------------
                covptot(j,i) = covptot(j,i)-max(d_zero, &
                               (covptot(j,i)-fccfg(j,i,k)) * &
                                dpevap/qpretot(j,i))
              else
                chng(j,i) = qxfg(iqqs,j,i)
              end if
              solqa(iqqv,iqqs,j,i) = solqa(iqqv,iqqs,j,i) + chng(j,i)
              solqa(iqqs,iqqv,j,i) = solqa(iqqs,iqqv,j,i) - chng(j,i)
              qxfg(iqqs,j,i)       = qxfg(iqqs,j,i) - chng(j,i)
            end if
          end do
        end do

        if ( stats ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              statssnowev(j,i,k) = chng(j,i)
            end do
          end do
        end if

      end if !lmicro

      !--------------------------------
      ! solver for the microphysics
      !--------------------------------
      ! Truncate sum of explicit sinks to size of bin
      ! this approach is inaccurate, but conserves -
      ! prob best can do with explicit (i.e. not implicit!) terms
      !----------------------------------------------------------
      sinksum(:,:,:) = d_zero
      lind3(:,:,:,:) = .false.
      !----------------------------
      ! collect sink terms and mark
      !----------------------------
      do i = ici1 , ici2
        do j = jci1 , jci2
          do jn = 1 , nqx
            do n = 1 , nqx
              sinksum(n,j,i) = sinksum(n,j,i) - solqa(n,jn,j,i)
            end do
          end do
        end do
      end do
      !---------------------------------------
      ! calculate overshoot and scaling factor
      !---------------------------------------
      do i = ici1 , ici2
        do j = jci1 , jci2
          do n = 1 , nqx
            ratio(n,j,i) = qx0(n,j,i,k)/max(sinksum(n,j,i),qx0(n,j,i,k))
          end do
        end do
      end do
      !--------------------------------------------------------
      ! now sort ratio to find out which species run out first
      !--------------------------------------------------------
      do i = ici1 , ici2
        do j = jci1 , jci2
          do n = 1 , nqx
            iorder(n,j,i) = -999
            lind1(n,j,i) = .true.
          end do
        end do
      end do
      do i = ici1 , ici2
        do j = jci1 , jci2
          do n = 1 , nqx
            do jn = 1 , nqx
              if ( lind1(jn,j,i) .and. ratio(jn,j,i) > dlowval ) then
                iorder(n,j,i) = jn
              end if
            end do
          end do
        end do
      end do
      do i = ici1 , ici2
        do j = jci1 , jci2
          do n = 1 , nqx
            if ( iorder(n,j,i) > 0 ) then
              lind1(iorder(n,j,i),j,i) = .false.
            end if
          end do
        end do
      end do
      !--------------------------------------------
      ! scale the sink terms, in the correct order,
      ! recalculating the scale factor each time
      !--------------------------------------------
      sinksum(:,:,:) = d_zero
      !----------------
      ! recalculate sum
      !----------------
      do i = ici1 , ici2
        do j = jci1 , jci2
          do n = 1 , nqx
            do jn = 1 , nqx
              jo = iorder(n,j,i)
              if ( jo > 0 ) then
                lind3(jo,jn,j,i) = solqa(jo,jn,j,i) < d_zero
                sinksum(jo,j,i) = sinksum(jo,j,i) - solqa(jo,jn,j,i)
              end if
            end do
          end do
        end do
      end do
      !---------------------------
      ! recalculate scaling factor
      !---------------------------
      do i = ici1 , ici2
        do j = jci1 , jci2
          do n = 1 , nqx
            jo = iorder(n,j,i)
            if ( jo > 0 ) then
              ratio(jo,j,i) = qx0(jo,j,i,k)/max(sinksum(jo,j,i),qx0(jo,j,i,k))
            end if
          end do
        end do
      end do
      !------
      ! scale
      !------
      do i = ici1 , ici2
        do j = jci1 , jci2
          do n = 1 , nqx
            do jn = 1 , nqx
              jo = iorder(n,j,i)
              if ( jo > 0 ) then
                if ( lind3(jo,jn,j,i) ) then
                  solqa(jo,jn,j,i) = solqa(jo,jn,j,i)*ratio(jo,j,i)
                  solqa(jn,jo,j,i) = solqa(jn,jo,j,i)*ratio(jo,j,i)
                end if
              end if
            end do
          end do
        end do
      end do

      ! FOREACH POINT, SOLVE A LINEAR SYSTEM

      do i = ici1 , ici2
        do j = jci1 , jci2

          ! Set the LHS of equation
          do n = 1 , nqx
            aamax = d_zero
            do jn = 1 , nqx
              ! Diagonals: microphysical sink terms+transport
              if ( jn == n ) then
                qlhs(jn,n) = d_one + rsemi*fallsink(n,j,i)
                do jo = 1 , nqx
                  qlhs(jn,n) = qlhs(jn,n) + rsemi*solqb(jo,jn,j,i)
                end do
                ! Non-diagonals: microphysical source terms
              else
                ! Here is the delta T - missing from doc.
                qlhs(jn,n) = -rsemi*solqb(jn,n,j,i)
              end if
            end do
          end do

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

          ! Set the RHS of equation

          do n = 1 , nqx
            ! Sum the explicit source and sink
            rexplicit = d_zero
            do jn = 1 , nqx
              ! Positive, since summed over 2nd index
              rexplicit = rexplicit + solqa(n,jn,j,i)
              if ( jn /= n ) then
                rexplicit =  rexplicit - &
                        (d_one-rsemi)*qx0(n,j,i,k)*solqb(jn,n,j,i) + &
                        (d_one-rsemi)*qx0(jn,j,i,k)*solqb(n,jn,j,i)
              end if
              rexplicit = rexplicit - &
                    (d_one-rsemi)*qx0(n,j,i,k)*fallsink(n,j,i)
            end do
            qxn(n) = qx0(n,j,i,k) + rexplicit
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
              if ( abs(xsum) > minqq ) ii = m
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
            ! Generalized precipitation flux
            ! this will be the source for the k
            pfplsx(n,j,i,k+1) = rsemi*fallsink(n,j,i) * &
              qxn(n)*rdtgdp(j,i) + (d_one-rsemi)*fallsink(n,j,i)* &
              qx0(n,j,i,k)*rdtgdp(j,i) ! kg/m2/s
            ! Calculate fluxes in and out of box for conservation of TL
            fluxq(n,j,i) = convsrce(n,j,i)+ fallsrce(n,j,i) - &
                         (fallsink(n,j,i)+convsink(n,j,i))*qxn(n)
            ! Calculate the water variables tendencies
            qxtendc(n,j,i,k) = qxtendc(n,j,i,k) + &
                                (qxn(n)-qx0(n,j,i,k))*oneodt
            ! Calculate the temperature tendencies
            if ( iphase(n) == 1 ) then
              ttendc(j,i,k) = ttendc(j,i,k) + &
                              wlhvocp*(qxn(n)-qx0(n,j,i,k) - &
                               fluxq(n,j,i))*oneodt
            else if ( iphase(n) == 2 ) then
              ttendc(j,i,k) = ttendc(j,i,k) + &
                              wlhsocp*(qxn(n)-qx0(n,j,i,k) - &
                              fluxq(n,j,i))*oneodt
            end if
          end do

        end do
      end do

      ! Couple tendencies with pressure
      do n = 1 , nqx
        do i = ici1 , ici2
          do j = jci1 , jci2
            qxten(j,i,k,n) =  qxtendc(n,j,i,k)*psb(j,i)
          end do
        end do
      end do
      do i = ici1 , ici2
        do j = jci1 , jci2
          tten(j,i,k) = ttendc(j,i,k)*psb(j,i)
        end do
      end do

    end do   ! kz : end of vertical loop

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
                tnew = tnew-wlhvocp*(qx0(n,j,i,k)+ &
                        (qxtendc(n,j,i,k)-tenkeep(n,j,i,k))*dt)
              else if ( iphase(n) == 2 ) then
                tnew = tnew-wlhsocp*(qx0(n,j,i,k)+ &
                        (qxtendc(n,j,i,k)-tenkeep(n,j,i,k))*dt)
              end if
              sumq1(j,i,k) = sumq1(j,i,k) + &
                (qx0(n,j,i,k)+(qxtendc(n,j,i,k)-tenkeep(n,j,i,k))*dt)* &
                (pfs(j,i,k+1)-pfs(j,i,k))*regrav
            end do
            sumh1(j,i,k) = sumh1(j,i,k)+(pfs(j,i,k+1)-pfs(j,i,k))*tnew
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
          dtgdp(j,i) = dt*egrav/(pfs(j,i,k+1)-pfs(j,i,k))
          rain = d_zero
          do n = 1 , nqx
            if ( iphase(n) == 1 ) then
              rain = rain+wlhvocp*dtgdp(j,i)*pfplsx(n,j,i,k+1)* & !k+1?
                       (pfs(j,i,k+1)-pfs(j,i,k))
            else if ( iphase(n) == 2 ) then
              rain = rain+wlhsocp*dtgdp(j,i)*pfplsx(n,j,i,k+1)* &
                      (pfs(j,i,k+1)-pfs(j,i,k))
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
  prcflxw(:,:) = d_zero
  prcflxc(:,:) = d_zero
  pfplsl(:,:,:) = d_zero
  pfplsn(:,:,:) = d_zero
  rainls(:,:,:) = d_zero

  !--------------------------------------------------------------------
  ! Copy general precip arrays back into FP arrays
  ! Add rain and liquid fluxes, ice and snow fluxes
  !--------------------------------------------------------------------

  ! Rain+liquid, snow+ice
  ! for each level k = 1 , kz, sum of the same phase elements
  do k = 1 , kz+1
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
      prainx = pfplsl(j,i,kz+1)*dt
      psnowx = pfplsn(j,i,kz+1)*dt
      if ( prainx > dlowval ) then
        rainnc(j,i) =  rainnc(j,i) + prainx   !mm
        lsmrnc(j,i) =  lsmrnc(j,i) + pfplsl(j,i,kz+1)
      end if
      if ( psnowx > dlowval ) then
        snownc(j,i) = snownc(j,i) + psnowx
        lsmrnc(j,i) =  lsmrnc(j,i) + pfplsn(j,i,kz+1)
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

    pure real(rkx) function phase(t)
      !phase = 1      if t > tzero
      !phase = 0      if t < rtice
      !0<phase < 1    if rtice < t < tzero
      implicit none
      real(rkx) , intent(in):: t
      phase = max(min(d_one,((max(rtice,min(tzero,t))-rtice)* &
                              rtwat_rtice_r)**2),d_zero)
    end function phase

    pure real(rkx) function edem(t,phase)
      implicit none
      real(rkx) , intent(in):: t , phase
      edem = phase * c5alvcp * (d_one/(t-c4les)**2) + &
               (d_one - phase) * c5alscp * (d_one/(t-c4ies)**2)
    end function edem

    pure real(rkx) function eldcpm(t)
      implicit none
      real(rkx) , intent(in):: t
      eldcpm = phase(t)*wlhvocp+(d_one-phase(t))*wlhsocp
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

end subroutine microphys

!  subroutine addpath_array(src,snk,proc2,zsqa,zsqb,beta,fg,j,i)
!    implicit none
!    real(rkx) , pointer , intent(inout) , dimension(:,:,:,:) :: zsqa , zsqb
!    real(rkx) , pointer , intent(inout) , dimension(:,:,:) :: fg
!    real(rkx) , pointer , intent(in) , dimension(:,:)  :: proc2
!    integer(ik4) , intent(in) :: src, snk
!    integer(ik4), intent(in)  :: i , j
!    real(rkx) , intent(in) :: beta
!    zsqa(j,i,src,snk) = zsqa(j,i,src,snk) + (d_one-beta)*proc2(j,i)
!    zsqa(j,i,snk,src) = zsqa(j,i,snk,src) - (d_one-beta)*proc2(j,i)
!    fg(j,i,src) = fg(j,i,src) + (d_one-beta)*proc2(j,i)
!    fg(j,i,snk) = fg(j,i,snk) - (d_one-beta)*proc2(j,i)
!    zsqb(j,i,src,snk) = zsqb(j,i,src,snk) + beta*proc2(j,i)
!  end subroutine addpath_array
!
!  subroutine addpath_real(src,snk,proc,zsqa,zsqb,beta,fg)
!    implicit none
!    real(rkx) , pointer , intent(inout) , dimension(:,:,:,:) :: zsqa , zsqb
!    real(rkx) , pointer , intent(inout) , dimension(:,:,:) :: fg
!    real(rkx) , intent(in) :: proc
!    integer(ik4) , intent(in) :: src , snk
!    integer(ik4) :: i , j
!    real(rkx) , intent(in) :: beta
!    do i = ici1 , ici2
!      do j = jci1 , jci2
!        zsqa(j,i,src,snk) = zsqa(j,i,src,snk) + (d_one-beta)*proc
!        zsqa(j,i,snk,src) = zsqa(j,i,snk,src) - (d_one-beta)*proc
!        fg(j,i,src) = fg(j,i,src) + (d_one-beta)*proc
!        fg(j,i,snk) = fg(j,i,snk) - (d_one-beta)*proc
!        zsqb(j,i,src,snk) = zsqb(j,i,src,snk) + beta*proc
!      end do
!    end do
!  end subroutine addpath_real

end module mod_cloud_s1

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

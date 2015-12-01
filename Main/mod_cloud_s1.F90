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

  use mod_cloud_variables
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
  use mod_runparams , only : ipptls
  use mod_runparams , only : stats , budget_compute , nssopt , kautoconv
  use mod_runparams , only : zauto_rate_khair , zauto_rate_kessl , &
                             zauto_rate_klepi
  use mod_runparams , only : rkconv , rcovpmin , rpecons
  use mod_runparams , only : ktau
  use mod_runparams , only : rtsrf

  implicit none

  private

  public :: allocate_mod_cloud_s1 , init_cloud_s1 , microphys

  real(rk8) :: zqtmst                                  ! 1/dt
  real(rk8) , pointer , dimension(:,:)   :: psb , rainnc, lsmrnc, snownc
  real(rk8) , pointer , dimension(:,:,:) :: pfcc       ! from atm
  integer(ik4) , pointer , dimension(:,:) :: ldmsk
  real(rk8) , pointer , dimension(:,:,:) :: paph       ! from atms
  real(rk8) , pointer , dimension(:,:,:) :: papf       ! from atms
  real(rk8) , pointer , dimension(:,:,:) :: zt         ! from atms
  real(rk8) , pointer , dimension(:,:,:) :: rho        ! from atms
  real(rk8) , pointer , dimension(:,:,:) :: pverv      ! from atms
  real(rk8) , pointer , dimension(:,:,:,:) :: zqxx     ! from atms
  real(rk8) , pointer , dimension(:,:,:) :: radheatrt  ! radiation heat rate
  real(rk8) , pointer , dimension(:,:,:) :: ztten      ! tendency of temperature
  real(rk8) , pointer , dimension(:,:,:) :: qdetr      ! conv. detr. water
  real(rk8) , pointer , dimension(:,:,:) :: rainls     ! Rain from here
  real(rk8) , pointer , dimension(:,:,:,:) :: zqxten   ! tendency of zqx

  ! Total water and enthalpy budget diagnostics variables
  ! marker for water phase of each species
  ! 0 = vapour, 1 = liquid, 2 = ice
  integer(ik4) , pointer , dimension(:) :: kphase
  ! marks melting linkage for ice categories
  ! ice->liquid, snow->rain
  integer(ik4) , pointer , dimension(:) :: imelt
  ! array for sorting explicit terms
  integer(ik4) , pointer , dimension(:,:,:) :: iorder
  logical , pointer , dimension(:) :: llfall
  logical , pointer , dimension(:,:,:) ::   llindex1
  logical , pointer , dimension(:,:,:,:) :: llindex3

  real(rk8) , pointer , dimension(:,:,:):: zsumh0 , zsumq0
  real(rk8) , pointer , dimension(:,:,:) :: zsumh1 , zsumq1
  real(rk8) , pointer , dimension(:,:,:) :: zerrorq , zerrorh
  real(rk8) , pointer , dimension(:,:,:):: ztentkeep
  real(rk8) , pointer , dimension(:,:,:,:) :: ztenkeep
  ! Mass variables
  real(rk8) , pointer , dimension(:,:) :: zdp     ! dp
  real(rk8) , pointer , dimension(:,:) :: zgdp    ! g/dp
  real(rk8) , pointer , dimension(:,:) :: zdtgdp  ! dt * g/dp
  real(rk8) , pointer , dimension(:,:) :: zrdtgdp ! dp / (dt * g)  [Kg/(m*s)]
  ! Microphysics
  real(rk8) , pointer , dimension(:,:) :: zchng
  real(rk8) , pointer , dimension(:,:) :: zchngmax
  real(rk8) , pointer , dimension(:,:) :: zicetot
  real(rk8) , pointer , dimension(:,:) :: prcflxw
  real(rk8) , pointer , dimension(:,:) :: prcflxc
  real(rk8) , pointer , dimension(:,:,:) :: dqsatdt
  ! for sedimentation source/sink terms
  real(rk8) , pointer , dimension(:,:,:) :: zfallsink
  real(rk8) , pointer , dimension(:,:,:) :: zfallsrce
  ! for convection detrainment source and subsidence source/sink terms
  real(rk8) , pointer , dimension(:,:) :: zcorqsice
  real(rk8) , pointer , dimension(:,:) :: zcorqsliq
  real(rk8) , pointer , dimension(:,:) :: zcorqsmix
  real(rk8) , pointer , dimension(:,:) :: zevaplimmix
  real(rk8) , pointer , dimension(:,:) :: zsnowaut
  real(rk8) , pointer , dimension(:,:,:) :: zconvsrce
  real(rk8) , pointer , dimension(:,:,:) :: zconvsink
  ! total rain frac: fractional occurence of precipitation (%)
  real(rk8) , pointer , dimension(:,:) :: zcovptot
  ! for condensation
  real(rk8) , pointer , dimension(:,:) :: koop
  real(rk8) , pointer , dimension(:,:) :: zcovpclr
  real(rk8) , pointer , dimension(:,:) :: zqpretot
  real(rk8) , pointer , dimension(:,:) :: zliqcld
  real(rk8) , pointer , dimension(:,:) :: zicecld
  real(rk8) , pointer , dimension(:,:) :: zsupsat
  real(rk8) , pointer , dimension(:,:) :: zsubsat
  real(rk8) , pointer , dimension(:,:) :: zmin
  real(rk8) , pointer , dimension(:,:) :: zlicld
  real(rk8) , pointer , dimension(:,:) :: zldefr
  real(rk8) , pointer , dimension(:,:) :: zqold
  real(rk8) , pointer , dimension(:,:) :: ztold
  real(rk8) , pointer , dimension(:,:) :: zdqs
  real(rk8) , pointer , dimension(:,:,:) :: ztcond
  ! tropopause
  real(rk8) , pointer , dimension(:,:) :: zpapfd
  real(rk8) , pointer , dimension(:,:) :: ztrpaus
  ! distance from the top of the cloud
  real(rk8) , pointer , dimension(:,:) :: zcldtopdist
  ! ice nuclei concentration
  real(rk8) , pointer , dimension(:,:) :: zicenuclei
  real(rk8) , pointer , dimension(:,:,:) :: zfoeewmt
  real(rk8) , pointer , dimension(:,:,:) :: zliq
  real(rk8) , pointer , dimension(:,:,:) :: zliqfrac
  real(rk8) , pointer , dimension(:,:,:) :: zicefrac
  real(rk8) , pointer , dimension(:,:,:) :: zli
  ! fluxes convergence of species
  real(rk8) , pointer , dimension(:,:,:) :: zfluxq
  real(rk8) , pointer , dimension(:,:,:) :: zratio
  real(rk8) , pointer , dimension(:,:,:) :: zsinksum
  real(rk8) , pointer , dimension(:,:,:) :: zfoeew
  ! ice water saturation
  real(rk8) , pointer , dimension(:,:,:) :: zqsice
  ! diagnostic mixed phase RH
  real(rk8) , pointer , dimension(:,:,:) :: zqsmix
  ! melting
  real(rk8) , pointer , dimension(:,:,:) :: zfoeeliq
  ! water saturation mixing ratio
  real(rk8) , pointer , dimension(:,:,:) :: zfoeeliqt
  ! ice saturation mixing ratio
  real(rk8) , pointer , dimension(:,:,:) :: zfoeeicet
  ! liq+rain sedim flux
  real(rk8) , pointer , dimension(:,:,:) :: zpfplsl
  ! ice+snow sedim flux
  real(rk8) , pointer , dimension(:,:,:) :: zpfplsn
  ! Flux of liquid
  real(rk8) , pointer , dimension(:,:,:) :: pfsqlf
  ! Flux of ice
  real(rk8) , pointer , dimension(:,:,:) :: pfsqif
  ! decoupled temperature tendency
  real(rk8) , pointer , dimension(:,:,:) :: zttendc
  ! detrainment from tiedtke scheme
  real(rk8) , pointer , dimension(:,:,:) :: zqdetr
  ! real(rk8) , pointer , dimension(:,:,:) :: zqdetr2
  ! fall speeds of three categories
  real(rk8) , pointer , dimension(:) :: zvqx
  ! n x n matrix storing the LHS of implicit solver
  real(rk8) , pointer , dimension(:,:,:,:) :: zqlhs
  ! explicit sources and sinks
  real(rk8) , pointer , dimension(:,:,:,:) :: zsolqa
  ! implicit sources and sinks
  real(rk8) , pointer , dimension(:,:,:,:) :: zsolqb
  ! decoupled mixing ratios tendency
  real(rk8) , pointer , dimension(:,:,:,:) :: zqxtendc
  ! j,i,n ! generalized precipitation flux
  real(rk8) , pointer , dimension(:,:,:,:) :: zpfplsx
  real(rk8) , public  , pointer, dimension(:,:,:,:) :: zqx0
  ! new values for zqxx at time+1
  real(rk8) , public  , pointer, dimension(:,:,:)   :: zqxn
  ! first guess values including precip
  real(rk8) , public  , pointer, dimension(:,:,:)   :: zqxfg
  ! first guess value for cloud fraction
  real(rk8) , public  , pointer, dimension(:,:,:)   :: fccfg
  ! relative humidity
  real(rk8) , public  , pointer, dimension(:,:,:)   :: relh
  ! saturation mixing ratio with respect to water
  real(rk8) , public  , pointer, dimension(:,:,:)   :: zqsliq
  ! turbulent erosion rate
  real(rk8) , pointer , dimension(:,:) :: zldifdt

  ! statistic only if stas =.true.
  real(rk8) , public  , pointer, dimension(:,:,:)   :: statssupw
  real(rk8) , public  , pointer, dimension(:,:,:)   :: statssupc
  real(rk8) , public  , pointer, dimension(:,:,:)   :: statserosw
  real(rk8) , public  , pointer, dimension(:,:,:)   :: statserosc
  real(rk8) , public  , pointer, dimension(:,:,:)   :: statsdetrw
  real(rk8) , public  , pointer, dimension(:,:,:)   :: statsdetrc
  real(rk8) , public  , pointer, dimension(:,:,:)   :: statsevapw
  real(rk8) , public  , pointer, dimension(:,:,:)   :: statsevapc
  real(rk8) , public  , pointer, dimension(:,:,:)   :: statscond1w
  real(rk8) , public  , pointer, dimension(:,:,:)   :: statscond1c
  real(rk8) , public  , pointer, dimension(:,:,:)   :: statscond2w
  real(rk8) , public  , pointer, dimension(:,:,:)   :: statscond2c
  real(rk8) , public  , pointer, dimension(:,:,:)   :: statsdepos
  real(rk8) , public  , pointer, dimension(:,:,:)   :: statsmelt
  real(rk8) , public  , pointer, dimension(:,:,:)   :: statsfrz
  real(rk8) , public  , pointer, dimension(:,:,:)   :: statsrainev
  real(rk8) , public  , pointer, dimension(:,:,:)   :: statssnowev
#ifdef USE_LAPACK
  integer(ik4) , pointer , dimension(:) :: ipivot
#endif

  real(rk8) , parameter :: zminqx = minqq ! 1.0D-8

  abstract interface
    pure function myzqe(q,m,f,l)
      use mod_constants , only : rk8
      implicit none
      real(rk8) :: myzqe
      real(rk8) , intent(in) :: q , m , f , l
    end function myzqe
  end interface

  procedure(myzqe) , pointer :: p_zqe => null()

!  interface addpath
!    module procedure addpath_array
!    module procedure addpath_real
!  end interface

  real(rk8) , parameter :: clfeps = 1.0D-6
  real(rk8) , parameter :: zerocf = lowcld - clfeps
  real(rk8) , parameter :: onecf  = hicld + clfeps
  contains

  subroutine allocate_mod_cloud_s1
    implicit none
    if ( ipptls /= 2 ) return

    call getmem1d(zvqx,1,nqx,'cmicro:zvqx')
    call getmem1d(imelt,1,nqx,'cmicro:imelt')
    call getmem1d(llfall,1,nqx,'cmicro:llfall')
    call getmem1d(kphase,1,nqx,'cmicro:kphase')
    call getmem2d(zldefr,jci1,jci2,ici1,ici2,'cmicro:zldefr')
    call getmem2d(zqold,jci1,jci2,ici1,ici2,'cmicro:zqold')
    call getmem2d(ztold,jci1,jci2,ici1,ici2,'cmicro:ztold')
    call getmem2d(zdqs,jci1,jci2,ici1,ici2,'cmicro:zdqs')
    call getmem2d(koop,jci1,jci2,ici1,ici2,'cmicro:koop')
    call getmem2d(zicetot,jci1,jci2,ici1,ici2,'cmicro:zicetot')
    call getmem2d(zchng,jci1,jci2,ici1,ici2,'cmicro:zchng')
    call getmem2d(zchngmax,jci1,jci2,ici1,ici2,'cmicro:zchngmax')
    call getmem2d(zdp,jci1,jci2,ici1,ici2,'cmicro:zdp')
    call getmem2d(zgdp,jci1,jci2,ici1,ici2,'cmicro:zgdp')
    call getmem2d(zdtgdp,jci1,jci2,ici1,ici2,'cmicro:zdtgdp')
    call getmem2d(zrdtgdp,jci1,jci2,ici1,ici2,'cmicro:zrdtgdp')
    call getmem2d(prcflxw,jci1,jci2,ici1,ici2,'cmicro:prcflxw')
    call getmem2d(prcflxc,jci1,jci2,ici1,ici2,'cmicro:prcflxc')
    call getmem2d(zcovptot,jci1,jci2,ici1,ici2,'cmicro:zcovptot')
    call getmem2d(zcovpclr,jci1,jci2,ici1,ici2,'cmicro:zcovpclr')
    call getmem2d(zqpretot,jci1,jci2,ici1,ici2,'cmicro:zqpretot')
    call getmem2d(zsnowaut,jci1,jci2,ici1,ici2,'cmicro:zsnowaut')
    call getmem2d(zcldtopdist,jci1,jci2,ici1,ici2,'cmicro:zcldtopdist')
    call getmem2d(zicenuclei,jci1,jci2,ici1,ici2,'cmicro:zicenuclei')
    call getmem2d(zlicld,jci1,jci2,ici1,ici2,'cmicro:zlicld')
    call getmem2d(zcorqsmix,jci1,jci2,ici1,ici2,'cmicro:zcorqsmix')
    call getmem2d(zevaplimmix,jci1,jci2,ici1,ici2,'cmicro:zevaplimmix')
    call getmem2d(zliqcld,jci1,jci2,ici1,ici2,'cmicro:zliqcld')
    call getmem2d(zicecld,jci1,jci2,ici1,ici2,'cmicro:zicecld')
    call getmem2d(zsupsat,jci1,jci2,ici1,ici2,'cmicro:zsupsat')
    call getmem2d(zsubsat,jci1,jci2,ici1,ici2,'cmicro:zsubsat')
    call getmem2d(zmin,jci1,jci2,ici1,ici2,'cmicro:zmin')
    call getmem2d(zcorqsice,jci1,jci2,ici1,ici2,'cmicro:zcorqsice')
    call getmem2d(zcorqsliq,jci1,jci2,ici1,ici2,'cmicro:zcorqsliq')
    call getmem2d(ztrpaus,jci1,jci2,ici1,ici2,'cmicro:ztrpaus')
    call getmem2d(zpapfd,jci1,jci2,ici1,ici2,'cmicro:zpapfd')
    call getmem2d(zldifdt,jci1,jci2,ici1,ici2,'cmicro:zldifdt')
    call getmem3d(ztcond,jci1,jci2,ici1,ici2,1,kz,'cmicro:ztcond')
    call getmem3d(zliqfrac,jci1,jci2,ici1,ici2,1,kz,'cmicro:zliqfrac')
    call getmem3d(zicefrac,jci1,jci2,ici1,ici2,1,kz,'cmicro:zicefrac')
    call getmem3d(zfoeewmt,jci1,jci2,ici1,ici2,1,kz,'cmicro:zfoeewmt')
    call getmem3d(zqsmix,jci1,jci2,ici1,ici2,1,kz,'cmicro:zqsmix')
    call getmem3d(zli,jci1,jci2,ici1,ici2,1,kz+1,'cmicro:zli')
    call getmem3d(iorder,jci1,jci2,ici1,ici2,1,nqx,'cmicro:iorder')
    call getmem3d(zttendc,jci1,jci2,ici1,ici2,1,kz,'cmicro:zttendc')
    call getmem3d(zconvsrce,jci1,jci2,ici1,ici2,1,nqx,'cmicro:zconvsrce')
    call getmem3d(zconvsink,jci1,jci2,ici1,ici2,1,nqx,'cmicro:zconvsink')
    call getmem3d(zfoeew,jci1,jci2,ici1,ici2,1,kz,'cmicro:zfoeew')
    call getmem3d(zfoeeliq,jci1,jci2,ici1,ici2,1,kz,'cmicro:zfoeeliq')
    call getmem3d(zqsice,jci1,jci2,ici1,ici2,1,kz,'cmicro:zqsice')
    call getmem4d(zqx0,jci1,jci2,ici1,ici2,1,kz,1,nqx,'cmicro:zqx0')
    call getmem3d(zqsliq,jci1,jci2,ici1,ici2,1,kz,'cmicro:zqsliq')
    call getmem3d(zfallsink,jci1,jci2,ici1,ici2,1,nqx,'cmicro:zfallsink')
    call getmem3d(zfallsrce,jci1,jci2,ici1,ici2,1,nqx,'cmicro:zfallsrce')
    call getmem3d(zfluxq,jci1,jci2,ici1,ici2,1,nqx,'cmicro:zfluxq')
    call getmem3d(zratio,jci1,jci2,ici1,ici2,1,nqx,'cmicro:zratio')
    call getmem3d(zsinksum,jci1,jci2,ici1,ici2,1,nqx,'cmicro:zsinksum')
    call getmem3d(dqsatdt,jci1,jci2,ici1,ici2,1,kz,'cmicro:dqsatdt')
    call getmem3d(zpfplsl,jci1,jci2,ici1,ici2,1,kz+1,'cmicro:zpfplsl')
    call getmem3d(zpfplsn,jci1,jci2,ici1,ici2,1,kz+1,'cmicro:zpfplsn')
    call getmem3d(pfsqlf,jci1,jci2,ici1,ici2,1,kz+1,'cmicro:pfsqlf')
    call getmem3d(pfsqif,jci1,jci2,ici1,ici2,1,kz+1,'cmicro:pfsqif')
    call getmem3d(zliq,jci1,jci2,ici1,ici2,1,kz+1,'cmicro:zliq')
    call getmem3d(zqxfg,jci1,jci2,ici1,ici2,1,nqx,'cmicro:zqxfg')
    call getmem3d(fccfg,jci1,jci2,ici1,ici2,1,kz,'cmicro:fccfg')
    call getmem3d(llindex1,jci1,jci2,ici1,ici2,1,nqx,'cmicro:llindex1')
    call getmem3d(zqdetr,jci1,jci2,ici1,ici2,1,kz,'cmicro:zqdetr')
    !zqdetr after evaporation
    !call getmem3d(zqdetr2,jci1,jci2,ici1,ici2,1,kz,'cmicro:zqdetr2')
    call getmem3d(zfoeeliqt,jci1,jci2,ici1,ici2,1,kz,'cmicro:zfoeeliqt')
    call getmem3d(zfoeeicet,jci1,jci2,ici1,ici2,1,kz,'cmicro:zfoeeicet')
    call getmem4d(zqxtendc,jci1,jci2,ici1,ici2,1,kz,1,nqx,'cmicro:zqxtendc')
    call getmem3d(zqxn,1,nqx,jci1,jci2,ici1,ici2,'cmicro:zqxn')
    call getmem4d(zqlhs,1,nqx,1,nqx,jci1,jci2,ici1,ici2,'cmicro:zqlhs')
    call getmem4d(zsolqa,jci1,jci2,ici1,ici2,1,nqx,1,nqx,'cmicro:zsolqa')
    call getmem4d(zsolqb,jci1,jci2,ici1,ici2,1,nqx,1,nqx,'cmicro:zsolqb')
    call getmem4d(llindex3,jci1,jci2,ici1,ici2,1,nqx,1,nqx,'cmicro:llindex3')
    call getmem4d(zpfplsx,jci1,jci2,ici1,ici2,1,kz+1,1,nqx,'cmicro:zpfplsx')
#ifdef USE_LAPACK
    call getmem1d(ipivot,1,nqx,'cmicro:ipivot')
#endif
    if ( budget_compute ) then
      call getmem3d(zsumq0,jci1,jci2,ici1,ici2,1,kz,'cmicro:zsumq0')
      call getmem3d(zsumh0,jci1,jci2,ici1,ici2,1,kz,'cmicro:zsumh0')
      call getmem3d(zsumq1,jci1,jci2,ici1,ici2,1,kz,'cmicro:zsumq1')
      call getmem3d(zsumh1,jci1,jci2,ici1,ici2,1,kz,'cmicro:zsumh1')
      call getmem3d(zerrorq,jci1,jci2,ici1,ici2,1,kz,'cmicro:zerrorq')
      call getmem3d(zerrorh,jci1,jci2,ici1,ici2,1,kz,'cmicro:zerrorh')
      call getmem4d(ztenkeep,jci1,jci2,ici1,ici2,1,kz,1,nqx,'cmicro:ztenkeep')
      call getmem3d(ztentkeep,jci1,jci2,ici1,ici2,1,kz,'cmicro:ztentkeep')
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
    use mod_runparams , only : vqxr , vqxi , vqxs , nssopt
    implicit none
    integer(ik4) :: n
    call assignpnt(atms%pb3d,paph)
    call assignpnt(atms%pf3d,papf)
    call assignpnt(atms%tb3d,zt)
    call assignpnt(atms%wpx3d,pverv)
    call assignpnt(atms%qxb3d,zqxx)
    call assignpnt(atms%rhob3d,rho)
    call assignpnt(atms%rhb3d,relh)
    call assignpnt(aten%qx,zqxten)
    call assignpnt(aten%t,ztten)
    call assignpnt(sfs%psb,psb)
    call assignpnt(sfs%rainnc,rainnc)
    call assignpnt(sfs%snownc,snownc)
    call assignpnt(fcc,pfcc)
    call assignpnt(pptnc,lsmrnc)
    call assignpnt(heatrt,radheatrt)
    call assignpnt(q_detr,qdetr)
    call assignpnt(rain_ls,rainls)
    call assignpnt(mddom%ldmsk,ldmsk)

    ! Define species phase, 0 = vapour, 1 = liquid, 2 = ice
    kphase(iqqv) = 0
    kphase(iqql) = 1
    kphase(iqqr) = 1
    kphase(iqqi) = 2
    kphase(iqqs) = 2

    ! Set up melting/freezing index,
    ! if an ice category melts/freezes, where does it go?

    imelt(iqqv) = -99
    imelt(iqql) = iqqi
    imelt(iqqr) = iqqs
    imelt(iqqi) = iqql
    imelt(iqqs) = iqqr

    ! Set the fall velocities
    zvqx(iqqv) = d_zero ! * sqrt(ZQX(JL,JK,IQV))
    zvqx(iqql) = d_zero ! * sqrt(ZQX(JL,JK,IQL))
    zvqx(iqqr) = vqxr   !4.0D0  * sqrt(ZQX(JL,JK,IQR))
    zvqx(iqqi) = vqxi   !0.15D0 * sqrt(ZQX(JL,JK,IQI))
    zvqx(iqqs) = vqxs   !1.0D0  * sqrt(ZQX(JL,JK,IQS))

    ! Set llfall
    do n = 1 , nqx
      if ( zvqx(n) > d_zero ) then
        llfall(n) = .true. !falling species
      else
        llfall(n) = .false.
      end if
    end do

    ! Sanity check
    if ( nssopt < 0 .or. nssopt > 3 ) then
      call fatal(__FILE__,__LINE__,'NSSOPT IN CLOUD MUST BE IN RANGE 0-3')
    end if

    select case(nssopt)
      case (0)
        p_zqe => noscheme_zqe
      case (1)
        p_zqe => tompkins_zqe
      case (2)
        p_zqe => lohmann_karcher_zqe
      case (3)
        p_zqe => gierens_zqe
    end select
  end subroutine init_cloud_s1

  subroutine microphys
    implicit none
    integer(ik4) :: i , j , k , n , m , jn , jo
    logical :: llo1
    real(rk8) :: zexplicit
    real(rk8) :: zfac , zfaci , zfacw , zcor , zfokoop
    real(rk8) :: zalfaw , zphases , zice , zdelta , ztmpl , &
                 ztmpi , ztnew , zqe , zrain , zpreclr , zarg
    ! local real variables for autoconversion rate constants
    real(rk8) :: alpha1 ! coefficient autoconversion cold cloud
    real(rk8) :: ztmpa
    real(rk8) :: zcfpr
    real(rk8) :: zlcrit
    real(rk8) :: zprecip
    ! real(rk8) :: zqadj
    real(rk8) :: zzrh
    real(rk8) :: zbeta , zbeta1
    ! local variables for condensation
    real(rk8) :: zcond , zdtdp , zcdmax , zrhc , zsig , &
                 zacond , zzdl , zlcondlim
    ! local variables for melting
    real(rk8) :: ztdiff
    real(rk8) :: zcons1
    ! constant for converting the fluxes unit measures
    real(rk8) :: prainx , psnowx
    ! local real constants for evaporation
    real(rk8) :: zdpr , zdenom , zdpevap , zevapi , zevapl , excess
    real(rk8) :: zdqsmixdt , zdqsicedt , zdqsliqdt
    ! real(rk8) :: zgdph_r
    ! constants for deposition process
    real(rk8) :: zvpice , zvpliq , zadd , zbdd , zcvds , &
                 zice0 , zinew , zinfactor , zrainaut , zsnowaut
    ! constants for condensation
    real(rk8) :: zdpmxdt , zwtot , zzzdt , zdtdiab , zdtforc , &
                 zqp , zqsat , zcond1 , zlevap
    ! constants for turbulent mixing erosion of clouds
    real(rk8) :: ze , zleros
#ifdef USE_LAPACK
    integer :: ires
#endif
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'microphys'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    if ( .not. fscheme ) return

    ! Set the default 1.D-14 = d_zero
    ! Define the inizial array zqx0
    do n = 1 , nqx
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( zqxx(j,i,k,n) > zepsec ) then
              zqx0(j,i,k,n) = zqxx(j,i,k,n)
            else
              zqx0(j,i,k,n) = zepsec
            end if
          end do
        end do
      end do
    end do

    ! Detrainment : [zqdetr] = kg/(m^2*s)
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( qdetr(j,i,k) > zepsec ) then
            zqdetr(j,i,k) = qdetr(j,i,k)
          else
            zqdetr(j,i,k) = d_zero
          end if
        end do
      end do
    end do

    !-----------------------------------
    ! initialization for cloud variables
    ! -------------------------------------
    ! Define zliq the function for mixed phase
    !     PHASES is calculated to distinguish the three cases:
    !     PHASES = 1            water phase
    !     PHASES = 0            ice phase
    !     0 < PHASES < 1        mixed phase
    ! Define pressure at full levels (half levels for ECMWF)
    ! PAPH = PROVISIONAL PRESSURE ON HALF LEVELS            (Pa)
    ! Define a new array for detrainment
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          zliq(j,i,k) = phases(zt(j,i,k))    ! zliq = zfoealfa in IFS
        end do
      end do
    end do

    ! Reset variables
    zcovptot(:,:)    = d_zero
    zcovpclr(:,:)    = d_zero
    zcldtopdist(:,:) = d_zero
    zpfplsx(:,:,:,:) = d_zero

    ! Decouple tendencies
    do n = 1 , nqx
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            zqxtendc(j,i,k,n) = zqxten(j,i,k,n)/psb(j,i)
          end do
        end do
      end do
    end do
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          zttendc(j,i,k) = ztten(j,i,k)/psb(j,i)
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
      zerrorq(:,:,:)       = d_zero
      zerrorh(:,:,:)       = d_zero
      zsumh0(:,:,:)        = d_zero
      zsumq0(:,:,:)        = d_zero
      zsumh1(:,:,:)        = d_zero
      zsumq1(:,:,:)        = d_zero
      ztentkeep(:,:,:)     = d_zero
      ztenkeep(:,:,:,:)    = d_zero

      ! Record the tendencies
      do n = 1 , nqx
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              ztenkeep(j,i,k,n) = zqxtendc(j,i,k,n)
            end do
          end do
        end do
      end do
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            ztentkeep(j,i,k) = zttendc(j,i,k)
          end do
        end do
      end do

      ! initialize the flux arrays
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            ! [ztnew] = K
            ztnew = zt(j,i,k)+dt*(zttendc(j,i,k)-ztentkeep(j,i,k))
            if ( k > 1 ) then
              zsumq0(j,i,k) = zsumq0(j,i,k-1) ! total water
              zsumh0(j,i,k) = zsumh0(j,i,k-1) ! liquid water temperature
            end if
            do n = 1 , nqx
              if ( kphase(n) == 1 ) then
                ztmpl = zqx0(j,i,k,n)+dt*(zqxtendc(j,i,k,n)-ztenkeep(j,i,k,n))
                ztmpi = d_zero
              else if ( kphase(n) == 2 ) then
                ztmpi = zqx0(j,i,k,n)+dt*(zqxtendc(j,i,k,n)-ztenkeep(j,i,k,n))
                ztmpl = d_zero
              end if
              zsumq0(j,i,k) = zsumq0(j,i,k) + &
                (zqx0(j,i,k,n)+(zqxtendc(j,i,k,n)-ztenkeep(j,i,k,n))*dt)* &
                (papf(j,i,k+1)-papf(j,i,k))*regrav
            end do
            ztnew = ztnew - wlhvocp*ztmpl - wlhsocp*ztmpi
            zsumq0(j,i,k) = zsumq0(j,i,k) + &
              (ztmpl+ztmpi)*(papf(j,i,k+1)-papf(j,i,k))*regrav    !(kg/m^2)
            ! Detrained water treated here
            zqe = zqdetr(j,i,k)*dt*egrav/(papf(j,i,k+1)-papf(j,i,k)) ! 1 ?
            if ( zqe > zminqx ) then
              zsumq0(j,i,k) = zsumq0(j,i,k)+zqdetr(j,i,k)*dt
              zalfaw = zliq(j,i,k)
              ztnew = ztnew-(wlhvocp*zalfaw+wlhsocp*(d_one-zalfaw))*zqe
            end if
            zsumh0(j,i,k) = zsumh0(j,i,k)+(papf(j,i,k+1)-papf(j,i,k))*ztnew
          end do
        end do
      end do
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            zsumh0(j,i,k) = zsumh0(j,i,k)/(papf(j,i,k+1)-papf(j,i,1))
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
          zphases = zliq(j,i,k)
          zfoeewmt(j,i,k) = min(foeewm(zt(j,i,k),zphases)/paph(j,i,k),d_half)
          zqsmix(j,i,k) = zfoeewmt(j,i,k)
          ! ep1 = rwat/rgas - d_one
          zqsmix(j,i,k) = zqsmix(j,i,k)/(d_one-ep1*zqsmix(j,i,k))
          !--------------------------------------------
          ! ice saturation T < 273K
          ! liquid water saturation for T > 273K
          !--------------------------------------------
          zdelta = delta(zt(j,i,k))
          zfoeew(j,i,k) = min((zdelta*foeeliq(zt(j,i,k)) + &
               (d_one-zdelta)*foeeice(zt(j,i,k)))/paph(j,i,k),d_half)
          zfoeew(j,i,k) = min(d_half,zfoeew(j,i,k))
          !qsi saturation mixing ratio with respect to ice
          !zqsice(j,i,k) = zfoeew(j,i,k)/(d_one-ep1*zfoeew(j,i,k))
          !ice water saturation
          zfoeeicet(j,i,k) = min(foeeice(zt(j,i,k))/paph(j,i,k),d_half)
          zqsice(j,i,k) = zfoeeicet(j,i,k)
          zqsice(j,i,k) = zqsice(j,i,k)/(d_one-ep1*zqsice(j,i,k))
          !----------------------------------
          ! liquid water saturation
          !----------------------------------
          !foeeliq is the saturation vapor pressure es(T)
          !the saturation mixing ratio is ws = es(T)/p *0.622
          !ws = ws/(-(d_one/eps - d_one)*ws)
          zfoeeliqt(j,i,k) = min(foeeliq(zt(j,i,k))/paph(j,i,k),d_half)
          zqsliq(j,i,k) = zfoeeliqt(j,i,k)
          zqsliq(j,i,k) = zqsliq(j,i,k)/(d_one-ep1*zqsliq(j,i,k))
        end do
      end do
    end do

    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          !-------------------------------------------------------------------
          ! Calculate liq/ice fractions (no longer a diagnostic relationship)
          !-------------------------------------------------------------------
          zli(j,i,k) = zqx0(j,i,k,iqql)+zqx0(j,i,k,iqqi)
          if ( zli(j,i,k) > zminqx ) then
            zliqfrac(j,i,k) = zqx0(j,i,k,iqql)/zli(j,i,k)
            zicefrac(j,i,k) = d_one-zliqfrac(j,i,k)
          else
            zliqfrac(j,i,k) = d_zero
            zicefrac(j,i,k) = d_zero
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
      do n = 1 , nqx
        do i = ici1 , ici2
          do j = jci1 , jci2
            zqxfg(j,i,n) = zqx0(j,i,k,n)
          end do
        end do
      end do

      do i = ici1 , ici2
        do j = jci1 , jci2
          fccfg(j,i,k) = min(hicld,max(pfcc(j,i,k),lowcld))
        end do
      end do

      ! Reset matrix so missing pathways are set
      zsolqb(:,:,:,:)  = d_zero
      zsolqa(:,:,:,:)  = d_zero
      zfallsrce(:,:,:) = d_zero
      zfallsink(:,:,:) = d_zero
      zconvsrce(:,:,:) = d_zero
      zconvsink(:,:,:) = d_zero
      zratio(:,:,:)    = d_zero
      zicetot(:,:)     = d_zero
      zsupsat(:,:)     = d_zero
      zsubsat(:,:)     = d_zero
      zlicld(:,:)      = d_zero
      zldefr(:,:)      = d_zero
      zqxn(:,:,:)      = d_zero
      zqlhs(:,:,:,:)   = d_zero

      ! Set j,i arrays to zero
      zqpretot(:,:)  = d_zero
      ! zqdetr2(:,:,:) = d_zero

      ! Set stats variables to zero
      if ( stats ) then
        statssupw(:,:,:) = d_zero
        statssupc(:,:,:) = d_zero
      end if

      ! Derived variables needed
      do i = ici1, ici2
        do j = jci1, jci2
          zdp(j,i) = papf(j,i,k+1) - papf(j,i,k)     !dp
          zgdp(j,i) = egrav/zdp(j,i)                 !g/dp  =(1/m)
          zdtgdp(j,i) = dt*zgdp(j,i)                 !(dt*g)/dp =(dt/m)
          zrdtgdp(j,i) = zdp(j,i)*(d_one/(dt*egrav)) !dp/(gdt)=m/dt  [Kg/m2/s]
          zqtmst = d_one/dt                          !1/dt
          !------------------------------------
          ! calculate dqs/dT
          !------------------------------------
          ! liquid
          zfacw     = r5les/((zt(j,i,k) - r4les)**2)
          zcor      = d_one/(d_one - ep1*zfoeeliqt(j,i,k))
          zdqsliqdt = zfacw*zcor*zqsliq(j,i,k)
          zcorqsliq(j,i) = d_one + wlhvocp*zdqsliqdt
          ! ice
          zfaci     = r5ies/((zt(j,i,k) - r4ies)**2)
          zcor      = d_one/(d_one - ep1*zfoeew(j,i,k))
          zdqsicedt = zfaci*zcor*zqsice(j,i,k)
          zcorqsice(j,i) = d_one + wlhsocp*zdqsicedt
          ! diagnostic mixed
          zalfaw    = zliq(j,i,k)
          zfac      = zalfaw*zfacw + (d_one - zalfaw)*zfaci
          zcor      = d_one/(d_one - ep1*zfoeewmt(j,i,k))
          zdqsmixdt = zfac*zcor*zqsmix(j,i,k)
          zcorqsmix(j,i) = d_one/(d_one + foeldcpm(zt(j,i,k))*zdqsmixdt)
          !--------------------------------
          ! evaporation/sublimation limits
          !--------------------------------
          zevaplimmix(j,i) = max((zqsmix(j,i,k) - zqx0(j,i,k,iqqv)) * &
                                  zcorqsmix(j,i),d_zero)
          !--------------------------------
          ! in-cloud consensate amount
          !--------------------------------
          ztmpa = d_one/fccfg(j,i,k)
          zliqcld(j,i) = zqx0(j,i,k,iqql)*ztmpa
          zicecld(j,i) = zqx0(j,i,k,iqqi)*ztmpa
          zlicld(j,i)  = zliqcld(j,i) + zicecld(j,i)
        end do
      end do
      !------------------------------------------------
      ! Evaporate very small amounts of liquid and ice
      !------------------------------------------------
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( zqx0(j,i,k,iqql) < zminqx ) then
            zsolqa(j,i,iqqv,iqql) =  zqx0(j,i,k,iqql)
            zsolqa(j,i,iqql,iqqv) = -zqx0(j,i,k,iqql)
            zqxfg(j,i,iqql) = zqxfg(j,i,iqql) - zqx0(j,i,k,iqql)
          end if
          if ( zqx0(j,i,k,iqqi) < zminqx ) then
            zsolqa(j,i,iqqv,iqqi) =  zqx0(j,i,k,iqqi)
            zsolqa(j,i,iqqi,iqqv) = -zqx0(j,i,k,iqqi)
            zqxfg(j,i,iqqi) = zqxfg(j,i,iqqi) - zqx0(j,i,k,iqqi)
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
            zfac  = d_one
          else
            if ( zt(j,i,k) >= tzero ) then
              zfac  = d_one
            else
              zfokoop = fokoop(zt(j,i,k),foeeliq(zt(j,i,k)),foeeice(zt(j,i,k)))
              zfac  = fccfg(j,i,k) + zfokoop*(d_one-fccfg(j,i,k))
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
          ! zqsmix?
          ! zcorqslsliq? shouldn't it be zcorqsmix?
          zsupsat(j,i) = max((zqx0(j,i,k,iqqv)-zfac*zqsmix(j,i,k)) / &
                              zcorqsliq(j,i),d_zero)
          ! e < esi, because for e > esi ice still present
          zsubsat(j,i) = min((zqx0(j,i,k,iqqv)-zfac*zqsmix(j,i,k))/  &
                              zcorqsliq(j,i),d_zero)

          if ( zsupsat(j,i) > zepsec ) then
            if ( zt(j,i,k) > thomo ) then
              ! turn supersaturation into liquid water
              zsolqa(j,i,iqql,iqqv) = zsolqa(j,i,iqql,iqqv)+zsupsat(j,i)
              zsolqa(j,i,iqqv,iqql) = zsolqa(j,i,iqqv,iqql)-zsupsat(j,i)
              zqxfg(j,i,iqql) = zqxfg(j,i,iqql)+zsupsat(j,i)
            else
              ! turn supersaturation into ice water
              zsolqa(j,i,iqqi,iqqv) = zsolqa(j,i,iqqi,iqqv)+zsupsat(j,i)
              zsolqa(j,i,iqqv,iqqi) = zsolqa(j,i,iqqv,iqqi)-zsupsat(j,i)
              zqxfg(j,i,iqqi) = zqxfg(j,i,iqqi)+zsupsat(j,i)
            end if
          else
            if ( zsubsat(j,i) < d_zero .and. &
                 fccfg(j,i,k) < zerocf .and. &
                 zli(j,i,k) > zepsec ) then
              ! turn subsaturation into vapor, where there is no cloud
              excess = zli(j,i,k)+zsubsat(j,i)
              if ( excess < d_zero ) then
                if ( zt(j,i,k) > thomo ) then
                  zevapl = max(-zli(j,i,k),-zevaplimmix(j,i))!*zqtmst
                  zsolqa(j,i,iqqv,iqql) = zsolqa(j,i,iqqv,iqql)-zevapl
                  zsolqa(j,i,iqql,iqqv) = zsolqa(j,i,iqql,iqqv)+zevapl
                  zqxfg(j,i,iqql) = zqxfg(j,i,iqql)+zevapl
                else
                  zevapi = max(-zli(j,i,k),-zevaplimmix(j,i))!*zqtmst
                  ! turn subsaturation into vapour
                  zsolqa(j,i,iqqv,iqqi) = zsolqa(j,i,iqqv,iqqi)-zevapi
                  zsolqa(j,i,iqqi,iqqv) = zsolqa(j,i,iqqi,iqqv)+zevapi
                  zqxfg(j,i,iqqi) = zqxfg(j,i,iqqi)+zevapi
                end if
              end if
            end if
          end if
        end do
      end do

      if ( stats ) then
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( zsupsat(j,i) > zepsec ) then
              if ( zt(j,i,k) > thomo ) then
                statssupw(j,i,k) = zsupsat(j,i)
              else
                statssupc(j,i,k) = zsupsat(j,i)
              end if
            else
              if ( zsubsat(j,i) < d_zero .and. &
                   fccfg(j,i,k) < zerocf .and. &
                   zli(j,i,k) > zepsec ) then
                excess = zli(j,i,k)+zsubsat(j,i)
                if ( excess < d_zero ) then
                  if ( zt(j,i,k) > thomo ) then
                    zevapl = max(-zli(j,i,k),-zevaplimmix(j,i))!*zqtmst
                    statssupw(j,i,k) =  statssupw(j,i,k) + zevapl
                  else
                    zevapi = max(-zli(j,i,k),-zevaplimmix(j,i))!*zqtmst
                    statssupc(j,i,k) =  statssupc(j,i,k) - zevapi
                  end if
                end if
              end if
            end if
          end do
        end do
      end if

      ! call addpath_real(iqql,iqqv,zsupsatl,zsolqa,zsolqb,d_zero,zqxfg)
      ! call addpath_real(iqqi,iqqv,zsupsati,zsolqa,zsolqb,d_zero,zqxfg)
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
      ! ZSOLQA/B:q(IQa,IQb)
      !
      ! Thus if ZSOLQA/B(IQL,IQV) = K where K > 0 then this is
      ! a source of IQL and a sink of IQV
      !
      ! put 'magic' source terms such as PLUDE from
      ! detrainment into explicit source/sink array diagnognal
      ! ZSOLQA(IQL,IQL)=PLUDE
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
      ! zchng(:,:) = d_zero
      if ( k < kz .and. k >= 1 ) then
        do i = ici1 , ici2
          do j = jci1 , jci2
            zqdetr(j,i,k) = zqdetr(j,i,k)*zdtgdp(j,i)  !kg/kg
            if ( zqdetr(j,i,k) > zminqx ) then
              !zice = 1 if T < 250, zice = 0 if T > 273
              zalfaw = zliq(j,i,k)
              zice   = d_one-zalfaw
              zconvsrce(j,i,iqql) = zalfaw*zqdetr(j,i,k)
              zconvsrce(j,i,iqqi) = zice*zqdetr(j,i,k)
              zsolqa(j,i,iqql,iqql) = zsolqa(j,i,iqql,iqql) + &
                                      zconvsrce(j,i,iqql)
              zsolqa(j,i,iqqi,iqqi) = zsolqa(j,i,iqqi,iqqi) + &
                                      zconvsrce(j,i,iqqi)
              zqxfg(j,i,iqql) = zqxfg(j,i,iqql)+zconvsrce(j,i,iqql)
              zqxfg(j,i,iqqi) = zqxfg(j,i,iqqi)+zconvsrce(j,i,iqqi)

              ! adjustment of the cloud fraction, that increases when qdetr
              ! adds liquid to the scheme
              ! zqe = (zqx0(j,i,k,iqqv)-fccfg(j,i,k)*zqsmix(j,i,k)) / &
              !        zqsice(j,i,k)) / (d_one-fccfg(j,i,k))
              ! zqe = max(d_zero,min(zqe,zqsmix(j,i,k))) !zqsice(j,i,k)))
              ! define the new cloud
              ! fccfg(j,i,k) = (relh(j,i,k)**0.25D0)* &
              !                (d_one-dexp((-100.0D0*(zqxfg(j,i,iqql) + &
              !                 zqxfg(j,i,iqqi))/qx3(j,i,k,iqi)
              !                sqrt((d_one-relh(j,i,k))*zqsice(j,i,k)))))
              ! fccfg(j,i,k) = dmin1(dmax1(fccfg(j,i,k),0.01D0),0.99D0)
              ! zchng(j,i) = min(zqdetr(j,i,k),(fccfg(j,i,k) - &
              !                            fccfg(j,i,k))*(zqsice(j,i,k)-zqe))
              ! zchng(j,i) = max(zchng(j,i),d_zero)
              ! zqdetr2(j,i,k) = zqdetr2(j,i,k) - zchng(j,i)
              ! if (zqdetr2(j,i,k) > d_zero) then
              !   fccfg(j,i,k) = fccfg(j,i,k)
              ! else
              !   fccfg(j,i,k) = pfcc(j,i,k)
              ! end if
              ! zsolqa(j,i,iqqv,iqql) = zsolqa(j,i,iqqv,iqql) + zchng(j,i)
              ! zsolqa(j,i,iqql,iqqv) = zsolqa(j,i,iqql,iqqv) - zchng(j,i)
              ! zqxfg(j,i,iqql) = zqxfg(j,i,iqql) - zchng(j,i)
            end if
          end do
        end do
      end if

      if ( stats ) then
        if ( k < kz .and. k >= 1 ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              statsdetrw(j,i,k) = zconvsrce(j,i,iqql)
              statsdetrc(j,i,k) = zconvsrce(j,i,iqqi)
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
        ! mixing (IBID., EQU. 30) rcldiff = 3.0D-06
        do i = ici1 , ici2
          do j = jci1 , jci2
            zldifdt(j,i) = rcldiff*dt
            if ( zqdetr(j,i,k) > d_zero ) zldifdt(j,i) = 5.0D0*zldifdt(j,i)
          end do
        end do
        !Increase by factor of 5 for convective points
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( zli(j,i,k) > zepsec ) then
              ze = zldifdt(j,i)*max(zqsmix(j,i,k)-zqx0(j,i,k,iqqv),d_zero)
              zleros = fccfg(j,i,k)*ze
              zleros = min(zleros,zevaplimmix(j,i))
              zleros = min(zleros,zli(j,i,k))
              zfac  = zliqfrac(j,i,k)*zleros
              zfaci = zicefrac(j,i,k)*zleros
              zsolqa(j,i,iqql,iqqv) = zsolqa(j,i,iqql,iqqv) - zfac
              zsolqa(j,i,iqqv,iqql) = zsolqa(j,i,iqqv,iqql) + zfac
              zsolqa(j,i,iqqi,iqqv) = zsolqa(j,i,iqqi,iqqv) - zfaci
              zsolqa(j,i,iqqv,iqqi) = zsolqa(j,i,iqqv,iqqi) + zfaci
              zqxfg(j,i,iqql) = zqxfg(j,i,iqql) - zfac
              zqxfg(j,i,iqqi) = zqxfg(j,i,iqqi) - zfaci
            end if
          end do
        end do

        if ( stats ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              if ( zli(j,i,k) > zepsec ) then
                ze = zldifdt(j,i)*max(zqsmix(j,i,k)-zqx0(j,i,k,iqqv),d_zero)
                zleros = fccfg(j,i,k)*ze
                zleros = min(zleros,zevaplimmix(j,i))
                zleros = min(zleros,zli(j,i,k))
                statserosw(j,i,k) = zliqfrac(j,i,k)*zleros
                statserosc(j,i,k) = zicefrac(j,i,k)*zleros
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
            zdtdp   = rovcp*zt(j,i,k)/paph(j,i,k)
            zdpmxdt = zdp(j,i)*zqtmst
            zwtot   = pverv(j,i,k)
            zwtot   = min(zdpmxdt,max(-zdpmxdt,zwtot))
            zzzdt   = radheatrt(j,i,k)
            zdtdiab = min(zdpmxdt*zdtdp,max(-zdpmxdt*zdtdp,zzzdt))*dt + &
                          wlhfocp*zldefr(j,i)
            ! zldefr = 0
            ! note: zldefr should be set to the difference between the mixed
            ! phase functions in the convection and cloud scheme, but this is
            ! not calculated in the IFS, so is zero and the functions must be
            ! the same. We should check what happens in RegCM

            zdtforc       = zdtdp*zwtot*dt + zdtdiab
            zqold(j,i)    = zqsmix(j,i,k)
            ztold(j,i)    = zt(j,i,k)
            ztcond(j,i,k) = zt(j,i,k)+zdtforc
            ztcond(j,i,k) = max(ztcond(j,i,k),160.0D0)
          end do
        end do
        !this loop's goal is to produce zdqs = zqsmix - zqold, where zqsmix is
        !reduced because of the condensation. so that zdqs is negative?
        do i = ici1 , ici2
          do j = jci1 , jci2
            zqp = d_one/paph(j,i,k)
            zphases = phases(ztcond(j,i,k))
            ! saturation mixing ratio ws
            zqsat = foeewm(ztcond(j,i,k),zphases)*zqp
            zqsat = min(d_half,zqsat)          ! ws < 0.5        WHY?
            zcor  = d_one/(d_one-ep1*zqsat)
            zqsat = zqsat*zcor
            zcond = (zqsmix(j,i,k)-zqsat) / &
                    (d_one + zqsat*foedem(ztcond(j,i,k),zphases))
            ztcond(j,i,k) = ztcond(j,i,k)+foeldcpm(ztcond(j,i,k))*zcond

            zphases = phases(ztcond(j,i,k))
            zqsmix(j,i,k) = zqsmix(j,i,k)-zcond
            zqsat = foeewm(ztcond(j,i,k),zphases)*zqp
            zqsat = min(d_half,zqsat)
            zcor = d_one/(d_one-ep1*zqsat)
            zqsat = zqsat*zcor
            zcond1 = (zqsmix(j,i,k)-zqsat)  / &
                     (d_one + zqsat*foedem(ztcond(j,i,k),zphases))
            ztcond(j,i,k) = ztcond(j,i,k)+foeldcpm(ztcond(j,i,k))*zcond1
            zqsmix(j,i,k) = zqsmix(j,i,k)-zcond1
          end do
        end do
        do i = ici1 , ici2
          do j = jci1 , jci2
            zdqs(j,i)     = zqsmix(j,i,k)-zqold(j,i)
            zqsmix(j,i,k) = zqold(j,i)
            ztcond(j,i,k) = ztold(j,i)
          end do
        end do

        !----------------------------------------------------------------
        ! zdqs(jl) > 0:  evaporation of clouds
        !----------------------------------------------------------------
        ! erosion term is linear in l
        ! changed to be uniform distribution in cloud region
        do i = ici1 , ici2
          do j = jci1 , jci2
            ! previous function based on delta distribution in cloud:
            if ( zdqs(j,i) > d_zero ) then
              !zlevap = C*min( dqs/dt , (qi+ql)/C )
              zlevap = fccfg(j,i,k)*min(zdqs(j,i),zlicld(j,i))
              zlevap = min(zlevap,zevaplimmix(j,i))
              zlevap = min(zlevap,max(zqsmix(j,i,k)-zqx0(j,i,k,iqqv),d_zero))
              zfac = zliqfrac(j,i,k)*zlevap
              zfaci = zicefrac(j,i,k)*zlevap
              zsolqa(j,i,iqqv,iqql) = zsolqa(j,i,iqqv,iqql) + zfac
              zsolqa(j,i,iqql,iqqv) = zsolqa(j,i,iqql,iqqv) - zfac
              zsolqa(j,i,iqqv,iqqi) = zsolqa(j,i,iqqv,iqqi) + zfaci
              zsolqa(j,i,iqqi,iqqv) = zsolqa(j,i,iqqi,iqqv) - zfaci
              zqxfg(j,i,iqql) = zqxfg(j,i,iqql) - zfac
              zqxfg(j,i,iqqi) = zqxfg(j,i,iqqi) - zfaci
            end if
          end do
        end do
        if ( stats ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              if ( zdqs(j,i) > d_zero ) then
                zlevap = fccfg(j,i,k)*min(zdqs(j,i),zlicld(j,i))
                zlevap = min(zlevap,zevaplimmix(j,i))
                zlevap = min(zlevap,max(zqsmix(j,i,k)-zqx0(j,i,k,iqqv),d_zero))
                statsevapw(j,i,k) = zliqfrac(j,i,k)*zlevap
                statsevapc(j,i,k) = zicefrac(j,i,k)*zlevap
              end if
            end do
          end do
        end if

        !----------------------------------------------------------------------
        ! zdqs(j,i) < 0: formation of clouds
        !----------------------------------------------------------------------
        ! (1) increase of cloud water in existing clouds
        zchng(:,:) = d_zero

        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( fccfg(j,i,k) > zerocf .and. &
                 zdqs(j,i) <= -zminqx ) then
              ! new limiter
              zchng(j,i) = max(-zdqs(j,i),d_zero)
              ! old limiter
              !  (significantly improves upper tropospheric humidity rms)
              zphases = zliq(j,i,k)
              if ( fccfg(j,i,k) > onecf ) then
                zcor = d_one/(d_one-ep1*zqsmix(j,i,k))
                zcdmax = (zqx0(j,i,k,iqqv)-zqsmix(j,i,k)) / &
                         (d_one + zcor*zqsmix(j,i,k)*foedem(zt(j,i,k),zphases))
              else
                zcdmax = (zqx0(j,i,k,iqqv)-fccfg(j,i,k)*zqsmix(j,i,k)) / &
                        fccfg(j,i,k)
              end if
              zchng(j,i) = max(min(zchng(j,i),zcdmax),d_zero)
              ! end old limiter
              zchng(j,i) = fccfg(j,i,k)*zchng(j,i)
              if ( zchng(j,i) < zminqx ) zchng(j,i) = d_zero
              !----------------------------------------------------------------
              ! All increase goes into liquid unless so cold cloud
              ! homogeneously freezes
              ! include new liquid formation in first guess value, otherwise
              ! liquid remains at cold temperatures until next timestep.
              !----------------------------------------------------------------
              if ( zt(j,i,k) > thomo ) then
                zsolqa(j,i,iqql,iqqv) = zsolqa(j,i,iqql,iqqv) + zchng(j,i)
                zsolqa(j,i,iqqv,iqql) = zsolqa(j,i,iqqv,iqql) - zchng(j,i)
                zqxfg(j,i,iqql) = zqxfg(j,i,iqql) + zchng(j,i)
              else
                zsolqa(j,i,iqqi,iqqv) = zsolqa(j,i,iqqi,iqqv) + zchng(j,i)
                zsolqa(j,i,iqqv,iqqi) = zsolqa(j,i,iqqv,iqqi) - zchng(j,i)
                zqxfg(j,i,iqqi) = zqxfg(j,i,iqqi) + zchng(j,i)
              end if
            end if
          end do
        end do

        if ( stats ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              if ( zt(j,i,k) > thomo ) then
                statscond1w(j,i,k) = zchng(j,i)
              else
                statscond1c(j,i,k) = zchng(j,i)
              end if
            end do
          end do
        end if

        zchng(:,:) = d_zero

        ! (2) generation of new clouds (dc/dt>0)
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( zdqs(j,i) <= -zminqx .and. &
                 fccfg(j,i,k) < onecf ) then
              !---------------------------
              ! critical relative humidity
              !---------------------------
              ! *RAMID*   REAL    BASE VALUE FOR CALCULATION OF RELATIVE
              !                   HUMIDITY THRESHOLD FOR ONSET OF STRATIFORM
              !                   CONDENSATION (TIEDTKE, 1993, EQUATION 24)
              zrhc = ramid !=0.8
              zsig = paph(j,i,k)/papf(j,i,kz+1)
              ! increase RHcrit to 1.0 towards the surface (sigma>0.8)
              if ( zsig > 0.8D0 ) then
                zrhc = ramid+(d_one-ramid)*((zsig-0.8D0)/0.2D0)**2
              end if
              !---------------------------
              ! supersaturation options
              !---------------------------
              if ( zt(j,i,k) >= tzero .or. nssopt == 0 ) then
                ! no ice supersaturation allowed
                zfac = d_one
              else
                ! ice supersaturation
                zfac = fokoop(zt(j,i,k),foeeliq(zt(j,i,k)),foeeice(zt(j,i,k)))
              end if
              zqe = p_zqe(zqx0(j,i,k,iqqv),fccfg(j,i,k), &
                          zqsmix(j,i,k),zli(j,i,k))
              if ( zqe >= zrhc*zqsmix(j,i,k)*zfac .and. &
                   zqe < zqsmix(j,i,k)*zfac ) then
                ! note: not **2 on 1-a term if zqe is used.
                ! added correction term zfac to numerator 15/03/2010
                zacond = -(d_one-fccfg(j,i,k))*zfac*zdqs(j,i) / &
                          max(d_two*(zfac*zqsmix(j,i,k)-zqe),zepsec)
                zacond = min(zacond,d_one-fccfg(j,i,k)) ! put the limiter back
                ! linear term:
                ! added correction term zfac 15/03/2010
                zchng(j,i) = -zfac*zdqs(j,i)*d_half*zacond !mine linear
                ! new limiter formulation
                ! zqsice(j,i,k)-zqe) / &
                ztmpa = d_one-fccfg(j,i,k)
                zzdl = d_two*(zfac*zqsmix(j,i,k)-zqe) / ztmpa
                ! added correction term zfac 15/03/2010
                if ( zfac*zdqs(j,i) < -zzdl ) then
                  ! zqsice(j,i,k)+zqx0(j,i,k,iqqv)
                  zlcondlim = (fccfg(j,i,k)-d_one)*zfac*zdqs(j,i)- &
                               zfac*zqsmix(j,i,k)+zqx0(j,i,k,iqqv)
                  zchng(j,i) = min(zchng(j,i),zlcondlim)
                end if
                zchng(j,i) = max(zchng(j,i),d_zero)
                if ( zchng(j,i) < zminqx ) then
                  zchng(j,i) = d_zero
                  zacond     = d_zero
                end if
                !-------------------------------------------------------------
                ! all increase goes into liquid unless so cold cloud
                ! homogeneously freezes
                ! include new liquid formation in first guess value, otherwise
                ! liquid remains at cold temperatures until next timestep.
                !-------------------------------------------------------------
                if ( zt(j,i,k) > thomo ) then
                  zsolqa(j,i,iqql,iqqv) = zsolqa(j,i,iqql,iqqv) + zchng(j,i)
                  zsolqa(j,i,iqqv,iqql) = zsolqa(j,i,iqqv,iqql) - zchng(j,i)
                  zqxfg(j,i,iqql) = zqxfg(j,i,iqql) + zchng(j,i)
                else ! homogeneous freezing
                  zsolqa(j,i,iqqi,iqqv) = zsolqa(j,i,iqqi,iqqv) + zchng(j,i)
                  zsolqa(j,i,iqqv,iqqi) = zsolqa(j,i,iqqv,iqqi) - zchng(j,i)
                  zqxfg(j,i,iqqi) = zqxfg(j,i,iqqi) + zchng(j,i)
                end if
              end if
            end if
          end do
        end do

        if ( stats ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              if ( zt(j,i,k) > thomo ) then
                statscond2w(j,i,k) = zchng(j,i)
              else
                statscond2c(j,i,k) = zchng(j,i)
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
        zchng(:,:) = d_zero

        do i = ici1 , ici2
          do j = jci1 , jci2
            !--------------------------------------------------------------
            ! Calculate distance from cloud top
            ! defined by cloudy layer below a layer with cloud frac <0.01
            ! ZDZ = ZDP(JL)/(ZRHO(JL)*RG)
            !--------------------------------------------------------------
            if ( k > 1 ) then
              if ( fccfg(j,i,k-1) < rcldtopcf .and. &
                   fccfg(j,i,k)  >= rcldtopcf ) then
                zcldtopdist(j,i) = d_zero
              else
                zcldtopdist(j,i) = zcldtopdist(j,i) + &
                             zdp(j,i)/(rho(j,i,k)*egrav)
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
            if ( zt(j,i,k) < tzero .and. zqxfg(j,i,iqql) > zminqx ) then
              zvpice = foeeice(zt(j,i,k)) !saturation vapor pressure wrt ice
              zvpliq = foeeliq(zt(j,i,k)) !saturation vapor pressure wrt liq
              ! Meyers et al 1992
              zicenuclei(j,i) = d_1000 * &
                exp(12.96D0*((zvpliq-zvpice)/zvpice)-0.639D0)
              !------------------------------------------------
              !   2.4e-2 is conductivity of air
              !   8.87 = 700**1/3 = density of ice to the third
              !------------------------------------------------
              zadd  = wlhs*(wlhs/(rwat*zt(j,i,k))-d_one)/(2.4D-2*zt(j,i,k))
              zbdd  = rwat*zt(j,i,k)*paph(j,i,k)/(2.21D0*zvpice)
              zcvds = 7.8D0*(zicenuclei(j,i)/rho(j,i,k))** &
                       0.666D0*(zvpliq-zvpice)/(8.87D0*(zadd+zbdd)*zvpice)
              !-----------------------------------------------------
              ! riceinit = 1.e-12_jprb is initial mass of ice particle
              !-----------------------------------------------------
              zice0 = max(zicecld(j,i), zicenuclei(j,i)*riceinit/rho(j,i,k))
              !------------------
              ! new value of ice condensate amount ( Rotstayn et al. (2000) )
              !------------------
              zice0 = max(zice0,d_zero)
              zcvds = max(zcvds,d_zero)
              zinew = (0.666D0*zcvds*dt+zice0**0.666D0)**1.5D0
              !---------------------------
              ! grid-mean deposition rate:
              !---------------------------
              zchng(j,i) = max(fccfg(j,i,k)*(zinew-zice0),d_zero)*2.0D0
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
              zchng(j,i) = min(zchng(j,i),zqxfg(j,i,iqql))
              !---------------------------------------------------------------
              ! at top of cloud, reduce deposition rate near cloud top to
              ! account for small scale turbulent processes, limited ice
              ! nucleation and ice fallout
              !---------------------------------------------------------------
              ! Fraction of deposition rate in cloud top layer
              ! rdepliqrefrate  = d_r10
              ! Depth of supercooled liquid water layer (m)
              ! rdepliqrefdepth = 500.0D0
              zinfactor = min(zicenuclei(j,i)/15000.0D0, d_one)
              zchng(j,i) = zchng(j,i)*min(zinfactor + (d_one-zinfactor)* &
                     (rdepliqrefrate+zcldtopdist(j,i)/rdepliqrefdepth),d_one)
              !--------------
              ! add to matrix
              !--------------
              zsolqa(j,i,iqqi,iqql) = zsolqa(j,i,iqqi,iqql) + zchng(j,i)
              zsolqa(j,i,iqql,iqqi) = zsolqa(j,i,iqql,iqqi) - zchng(j,i)
              zqxfg(j,i,iqqi) = zqxfg(j,i,iqqi) + zchng(j,i)
              zqxfg(j,i,iqql) = zqxfg(j,i,iqql) - zchng(j,i)
            end if
          end do
        end do

        if ( stats ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              statsdepos(j,i,k) = zchng(j,i)
            end do
          end do
        end if

        do i = ici1 , ici2
          do j = jci1 , jci2
            ztmpa = d_one/fccfg(j,i,k)
            zliqcld(j,i) = zqxfg(j,i,iqql)*ztmpa
            zicecld(j,i) = zqxfg(j,i,iqqi)*ztmpa
          end do
        end do
        !------------------------------------------------------------------
        !  SEDIMENTATION/FALLING OF *ALL* MICROPHYSICAL SPECIES
        !     now that rain and snow species are prognostic
        !     the precipitation flux can be defined directly level by level
        !     There is no vertical memory required from the flux variable
        !------------------------------------------------------------------
        do n = 1 , nqx
          if ( llfall(n) ) then
            do i = ici1 , ici2
              do j = jci1 , jci2
                ! Source from layer above
                if ( k > 1 ) then
                  zfallsrce(j,i,n) = zpfplsx(j,i,k,n)*zdtgdp(j,i)
                  zsolqa(j,i,n,n) = zsolqa(j,i,n,n) + zfallsrce(j,i,n)
                  zqxfg(j,i,n) = zqxfg(j,i,n) + zfallsrce(j,i,n)
                  zqpretot(j,i)= zqpretot(j,i) + zqxfg(j,i,n)
                end if
                ! Sink to next layer, constant fall speed
                zfallsink(j,i,n) = zdtgdp(j,i)*zvqx(n)*rho(j,i,k)   !Kg/Kg
              end do
            end do
          end if  !llfall
        end do ! n

        !---------------------------------------------------------------
        ! Precip cover overlap using MAX-RAN Overlap
        ! Since precipitation is now prognostic we must
        !   1) apply an arbitrary minimum coverage (0.3) if precip>0
        !   2) abandon the 2-flux clr/cld treatment
        !   3) Thus, since we have no memory of the clear sky precip
        !      fraction, we mimic the previous method by reducing
        !      ZCOVPTOT(JL), which has the memory, proportionally with
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
            if ( zqpretot(j,i) > zepsec ) then
              zcovptot(j,i) = d_one - ((d_one-zcovptot(j,i))*&
                              (d_one - max(fccfg(j,i,k),fccfg(j,i,k-1)))/&
                              (d_one - min(fccfg(j,i,k-1),hicld)))
              zcovptot(j,i) = max(zcovptot(j,i),rcovpmin)   !rcovpmin = 0.1
              ! clear sky proportion
              zcovpclr(j,i) = max(d_zero,zcovptot(j,i)-fccfg(j,i,k))
            else
              zcovptot(j,i) = d_zero  ! no flux - reset cover
              zcovpclr(j,i) = d_zero  ! reset clear sky proportion
            end if
          end do
        end do
        !---------------------------------------------------------------
        !                         AUTOCONVERSION
        !---------------------------------------------------------------

        ! Warm clouds

        select case (kautoconv)
          case (1) ! Klein & Pincus (2000)
            do i = ici1 , ici2
              do j = jci1 , jci2
                if ( zliqcld(j,i) > zepsec ) then
                  zsolqb(j,i,iqql,iqqv) = d_zero
                  zsolqb(j,i,iqqr,iqql) = zsolqb(j,i,iqqr,iqql) + &
                    dt*zauto_rate_klepi * (zqx0(j,i,k,iqql)**(2.3D0))
                  zsolqa(j,i,iqqr,iqql) = d_zero
                end if
              end do
            end do
          case (2) ! Khairoutdinov and Kogan (2000)
            do i = ici1 , ici2
              do j = jci1 , jci2
                if ( zliqcld(j,i) > zepsec ) then
                  zsolqb(j,i,iqql,iqqv) = d_zero
                  zsolqb(j,i,iqqr,iqql) = zsolqb(j,i,iqqr,iqql) + &
                    dt*zauto_rate_khair*(zqx0(j,i,k,iqql)**(zauto_expon_khair))
                end if
              end do
            end do
          case (3) ! Kessler(1969)
            do i = ici1 , ici2
              do j = jci1 , jci2
                if ( zliqcld(j,i) > zepsec ) then
                  zsolqb(j,i,iqql,iqqv) = d_zero
                  if ( zqx0(j,i,k,iqql) > zautocrit_kessl ) then
                    zsolqa(j,i,iqqr,iqql) = zsolqa(j,i,iqqr,iqql) - &
                                  zauto_rate_kessl*zautocrit_kessl
                    zsolqa(j,i,iqql,iqqr) = zsolqa(j,i,iqql,iqqr) + &
                                  zauto_rate_kessl*zautocrit_kessl
                    zsolqb(j,i,iqqr,iqql) = zsolqb(j,i,iqqr,iqql) + &
                                  dt*zauto_rate_kessl
                  end if
                end if
              end do
            end do
          case (4) ! Sundqvist
            do i = ici1 , ici2
              do j = jci1 , jci2
                if ( zliqcld(j,i) > zepsec ) then
                  zsolqb(j,i,iqql,iqqv) = d_zero
                  alpha1 = rkconv*dt
                  ! modify autoconversion threshold dependent on:
                  ! land (polluted, high ccn, smaller droplets, higher
                  !       threshold)
                  ! sea  (clean, low ccn, larger droplets, lower threshold)
                  if ( ldmsk(j,i) == 0 ) then  ! landmask =0 land, =1 ocean
                    ! THRESHOLD VALUE FOR RAIN AUTOCONVERSION OVER LAND
                    zlcrit = rclcrit_land ! landrclcrit_land = 5.e-4_jprb
                  else
                    zlcrit = rclcrit_sea  ! oceanrclcrit_sea  = 3.e-4_jprb
                  end if
                  !-----------------------------------------------------------
                  ! parameters for cloud collection by rain and snow.
                  ! note that with new prognostic variable it is now possible
                  ! to replace this with an explicit collection parametrization
                  !-----------------------------------------------------------
                  zprecip = (zpfplsx(j,i,k,iqqs)+zpfplsx(j,i,k,iqqr)) / &
                             max(zepsec,zcovptot(j,i))
                  zcfpr = d_one + rprc1*sqrt(max(zprecip,d_zero))
                  alpha1 = alpha1*zcfpr
                  zlcrit = zlcrit/max(zcfpr,zepsec)
                  ! security for exp for some compilers
                  if ( zliqcld(j,i)/zlcrit < 25.0D0 ) then
                    zrainaut = alpha1*(d_one - exp(-(zliqcld(j,i)/zlcrit)**2))
                  else
                    zrainaut = alpha1
                  end if
                  !-----------------------
                  ! rain freezes instantly
                  !-----------------------
                  if ( zt(j,i,k) <= tzero ) then
                    zsolqb(j,i,iqqs,iqql) = zsolqb(j,i,iqqs,iqql)+zrainaut
                  else
                    zsolqb(j,i,iqqr,iqql) = zsolqb(j,i,iqqr,iqql)+zrainaut
                  end if
                end if
              end do
            end do
        end select

        ! Cold clouds
        ! Snow Autoconversion rate follow Lin et al. 1983
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( zt(j,i,k) <= tzero ) then
              if ( zicecld(j,i) > zepsec ) then
                alpha1 = dt*1.0D-3*exp(0.025*(zt(j,i,k)-tzero))
                zlcrit = rlcritsnow
                zarg = (zicecld(j,i)/zlcrit)**2
                if ( zarg < 25.0D0 ) then
                  zsnowaut = alpha1 * (d_one - exp(-zarg))
                else
                  zsnowaut = alpha1
                end if
                zsolqb(j,i,iqqs,iqqi) = zsolqb(j,i,iqqs,iqqi) + zsnowaut
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

        do n = 1 , nqx
          if ( kphase(n) == 2 ) then
            do i = ici1, ici2
              do j = jci1, jci2
                zicetot(j,i) = zicetot(j,i) + zqxfg(j,i,n)
              end do
            end do
          end if
        end do

        zchngmax(:,:) = d_zero

        do i = ici1, ici2
          do j = jci1, jci2
            if ( zicetot(j,i) > zepsec .and. zt(j,i,k) > tzero ) then
              ! Calculate subsaturation
              ! zqsice(j,i,k)-zqx0(j,i,k,iqqv),d_zero)
              zsubsat(j,i) = max(zqsmix(j,i,k)-zqx0(j,i,k,iqqv),d_zero)
              ! Calculate difference between dry-bulb (zt)  and the temperature
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
              ztdiff = zt(j,i,k)-tzero ! - zsubsat * &
              !    (ztw1+ztw2*(paph(j,i,k)-ztw3)-ztw4*(zt(j,i,k)-ztw5))
              ! Ensure ZCONS1 is positive so that ZMELTMAX = 0 if ZTDMTW0 < 0
              zcons1 = d_one ! abs(dt*(d_one + d_half*ztdiff)/rtaumel)
              zchngmax(j,i) = max(ztdiff*zcons1*zrldcp,d_zero)
            end if
          end do
        end do

        ! Loop over frozen hydrometeors (kphase == 2 (ice, snow))

        zchng(:,:) = d_zero
        do n = 1, nqx
          if ( kphase(n) == 2 ) then
            m = imelt(n) ! imelt(iqqi)=iqql, imelt(iqqs)=iqqr
            if ( m < 0 ) cycle
            do i = ici1 , ici2
              do j = jci1 , jci2
                if ( zchngmax(j,i) > zepsec .and. &
                     zicetot(j,i) > zepsec ) then
                  zphases = zqxfg(j,i,n)/zicetot(j,i)
                  zchng(j,i) = min(zqxfg(j,i,n),zphases*zchngmax(j,i))
                  zchng(j,i) = max(zchng(j,i),d_zero)
                  ! n = iqqi,iqqs; m = iqql,iqqr
                  zqxfg(j,i,n) =  zqxfg(j,i,n) - zchng(j,i)
                  zqxfg(j,i,m) =  zqxfg(j,i,m) + zchng(j,i)
                  zsolqa(j,i,m,n) =  zsolqa(j,i,m,n) + zchng(j,i)
                  zsolqa(j,i,n,m) =  zsolqa(j,i,n,m) - zchng(j,i)
                end if
              end do
            end do
          end if
        end do

        if ( stats ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              statsmelt(j,i,k) = zchng(j,i)
            end do
          end do
        end if

        !------------------------------------------------------------!
        !                         FREEZING                           !
        !------------------------------------------------------------!
        !Freezing of rain.
        !All rain freezes in a timestep if the temperature is below 0 C
        !calculate sublimation latent heat

        do i = ici1 , ici2
          do j = jci1 , jci2
            zchngmax(j,i) = max((tzero-zt(j,i,k))*zrldcp,d_zero)
          end do
        end do

        zchng(:,:) = d_zero

        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( zchngmax(j,i) > zepsec .and. zqxfg(j,i,iqqr) > zepsec ) then
              zchng(j,i) = min(zqxfg(j,i,iqqr),zchngmax(j,i))
              zchng(j,i) = max(zchng(j,i),d_zero)
              zsolqa(j,i,iqqs,iqqr) = zsolqa(j,i,iqqs,iqqr) + zchng(j,i)
              zsolqa(j,i,iqqr,iqqs) = zsolqa(j,i,iqqr,iqqs) - zchng(j,i)
            end if
          end do
        end do

        if ( stats ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              statsfrz(j,i,k) = zchng(j,i)
            end do
          end do
        end if

        ! Freezing of liquid

        do i = ici1 , ici2
          do j = jci1 , jci2
            zchngmax(j,i) = max((thomo-zt(j,i,k))*zrldcp,d_zero)
          end do
        end do

        zchng(:,:) = d_zero

        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( zchngmax(j,i) > zepsec .and. zqxfg(j,i,iqql) > zepsec ) then
              zchng(j,i) = min(zqxfg(j,i,iqql),zchngmax(j,i))
              zchng(j,i) = max(zchng(j,i),d_zero)
              zsolqa(j,i,iqqi,iqql) = zsolqa(j,i,iqqi,iqql) + zchng(j,i)
              zsolqa(j,i,iqql,iqqi) = zsolqa(j,i,iqql,iqqi) - zchng(j,i)
              zqxfg(j,i,iqqi) = zqxfg(j,i,iqqi) + zchng(j,i)
              zqxfg(j,i,iqql) = zqxfg(j,i,iqql) - zchng(j,i)
            end if
          end do
        end do

        if ( stats ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              statsfrz(j,i,k) = statsfrz(j,i,k) + zchng(j,i)
            end do
          end do
        end if
        !---------------------------------------------------------------
        !                         EVAPORATION
        !---------------------------------------------------------------
        !------------------------------------------------------------
        ! recalculate zqpretot since melting term may have changed it
        ! rprecrhmax = 0.7 is the threshold for the clear-sky RH that
        ! can be reached by evaporation of precipitation. THis assumption
        ! is done to prevent the gridbox saturating due to the evaporation
        ! of precipitation occuring in a portion of the grid
        !------------------------------------------------------------
        !do n = 1 , nqx
        do i = ici1 , ici2
          do j = jci1 , jci2
           ! if ( kphase(n) == 2 ) then
           !   zqpretot(j,i) = zqpretot(j,i)+zqxfg(j,i,n)
           zqpretot(j,i) = zqxfg(j,i,iqqs) + zqxfg(j,i,iqqr)
           ! end if
          end do
        end do
        !end do

        ! rain
        zchng(:,:) = d_zero
        do i = ici1 , ici2
          do j = jci1 , jci2
            zzrh = rprecrhmax + &
              (d_one-rprecrhmax)*zcovpclr(j,i)/(d_one-fccfg(j,i,k))
            zzrh = min(max(zzrh,rprecrhmax),d_one)
            ! This is a critical relative humidity that is used to limit
            ! moist environment to prevent the gridbox saturating when
            ! only part of the gridbox has evaporating precipitation
            zqe = (zqx0(j,i,k,iqqv)-fccfg(j,i,k)*zqsliq(j,i,k)) / &
                      (d_one-fccfg(j,i,k))
            !---------------------------------------------
            ! humidity in moistest zcovpclr part of domain
            !---------------------------------------------
            zqe = max(d_zero,min(zqe,zqsliq(j,i,k)))
            llo1 = zcovpclr(j,i) > zepsec .and. &
            !      zqpretot(jl) > zepsec .and. &
                   zqxfg(j,i,iqqr) > zminqx .and. &
                   zqe < zzrh*zqsliq(j,i,k)
            if ( llo1 ) then
              ! note: units of zpreclr and zqpretot differ
              !       zqpretot is a mixing ratio (hence "q" in name)
              !       zpreclr is a rain flux
              zpreclr = zqpretot(j,i)*zcovpclr(j,i)/(zcovptot(j,i)*zdtgdp(j,i))
              !--------------------------------------
              ! actual microphysics formula in zbeta
              !--------------------------------------
              ! sensitivity test showed multiply rain evap rate by 0.5
              zbeta1 = sqrt(paph(j,i,k)/papf(j,i,kz+1))/5.09D-3*zpreclr / &
                       max(zcovpclr(j,i),zepsec)

              if ( zbeta1 >= d_zero ) then
                 zbeta = egrav*rpecons*d_half*(zbeta1)**0.5777D0
                 zdenom = d_one + zbeta*dt*zcorqsliq(j,i)
                 zdpr = zcovpclr(j,i) * zbeta * &
                   (zqsliq(j,i,k)-zqe)/zdenom*zdp(j,i)*regrav
                 zdpevap = zdpr*zdtgdp(j,i)
                !---------------------------------------------------------
                ! add evaporation term to explicit sink.
                ! this has to be explicit since if treated in the implicit
                ! term evaporation can not reduce rain to zero and model
                ! produces small amounts of rainfall everywhere.
                !---------------------------------------------------------

                ! evaporate rain
                zchng(j,i) = min(zdpevap,zqxfg(j,i,iqqr))
                !-------------------------------------------------------------
                ! reduce the total precip coverage proportional to evaporation
                !-------------------------------------------------------------
                zcovptot(j,i) = zcovptot(j,i) - max(d_zero, &
                           (zcovptot(j,i)-fccfg(j,i,k))*zdpevap/zqpretot(j,i))
              else
                zchng(j,i) = zqxfg(j,i,iqqr)
              end if
              zsolqa(j,i,iqqv,iqqr) = zsolqa(j,i,iqqv,iqqr) + zchng(j,i)
              zsolqa(j,i,iqqr,iqqv) = zsolqa(j,i,iqqr,iqqv) - zchng(j,i)
              zqxfg(j,i,iqqr)       = zqxfg(j,i,iqqr) - zchng(j,i)
            end if
          end do
        end do

        if ( stats ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              statsrainev(j,i,k) = zchng(j,i)
            end do
          end do
        end if

        ! snow

        zchng(:,:) = d_zero

        do i = ici1 , ici2
          do j = jci1 , jci2
            zzrh = rprecrhmax + (d_one-rprecrhmax) * &
                   zcovpclr(j,i)/(d_one-fccfg(j,i,k))
            zzrh = min(max(zzrh,rprecrhmax),d_one)
            zqe = (zqx0(j,i,k,iqqv)-fccfg(j,i,k) * &
                   zqsice(j,i,k))/(d_one-fccfg(j,i,k))
            !---------------------------------------------
            ! humidity in moistest zcovpclr part of domain
            !---------------------------------------------
            zqe = max(d_zero,min(zqe,zqsice(j,i,k)))
            llo1 = zcovpclr(j,i) > zepsec .and. &
                   zqxfg(j,i,iqqs) > zminqx .and. &
                   zqe < zzrh*zqsice(j,i,k)
            if ( llo1 ) then
              ! note: units of zpreclr and zqpretot differ
              !       zqpretot is a mixing ratio (hence "q" in name)
              !       zpreclr is a rain flux
              zpreclr = zqpretot(j,i)*zcovpclr(j,i)/(zcovptot(j,i)*zdtgdp(j,i))
              !--------------------------------------
              ! actual microphysics formula in zbeta
              !--------------------------------------
               zbeta1 = sqrt(paph(j,i,k)/&
                        papf(j,i,kz+1))/5.09D-3*zpreclr/&
                        max(zcovpclr(j,i),zepsec)

              if ( zbeta1 >= d_zero ) then
                zbeta = egrav*rpecons*(zbeta1)**0.5777D0
                !rpecons = alpha1
                zdenom = d_one + zbeta*dt*zcorqsice(j,i)
                zdpr = zcovpclr(j,i) * zbeta * &
                       (zqsice(j,i,k)-zqe)/zdenom*zdp(j,i)*regrav
                zdpevap = zdpr*zdtgdp(j,i)
                !---------------------------------------------------------
                ! add evaporation term to explicit sink.
                ! this has to be explicit since if treated in the implicit
                ! term evaporation can not reduce snow to zero and model
                ! produces small amounts of snowfall everywhere.
                !---------------------------------------------------------
                ! evaporate snow
                zchng(j,i) = min(zdpevap,zqxfg(j,i,iqqs))
                zchng(j,i) = max(zchng(j,i),d_zero)
                !-------------------------------------------------------------
                ! reduce the total precip coverage proportional to evaporation
                !-------------------------------------------------------------
                zcovptot(j,i) = zcovptot(j,i)-max(d_zero, &
                               (zcovptot(j,i)-fccfg(j,i,k)) * &
                                zdpevap/zqpretot(j,i))
              else
                zchng(j,i) = zqxfg(j,i,iqqs)
              end if
              zsolqa(j,i,iqqv,iqqs) = zsolqa(j,i,iqqv,iqqs) + zchng(j,i)
              zsolqa(j,i,iqqs,iqqv) = zsolqa(j,i,iqqs,iqqv) - zchng(j,i)
              zqxfg(j,i,iqqs)       = zqxfg(j,i,iqqs) - zchng(j,i)
            end if
          end do
        end do

        if ( stats ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              statssnowev(j,i,k) = zchng(j,i)
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
      zsinksum(:,:,:) = d_zero
      llindex3(:,:,:,:) = .false.
      !----------------------------
      ! collect sink terms and mark
      !----------------------------
      do jn = 1 , nqx
        do n = 1 , nqx
          do i = ici1 , ici2
            do j = jci1 , jci2
              zsinksum(j,i,n) = zsinksum(j,i,n) - zsolqa(j,i,n,jn)
            end do
          end do
        end do
      end do
      !---------------------------------------
      ! calculate overshoot and scaling factor
      !---------------------------------------
      do n = 1 , nqx
        do i = ici1 , ici2
          do j = jci1 , jci2
            zratio(j,i,n) = zqx0(j,i,k,n)/max(zsinksum(j,i,n),zqx0(j,i,k,n))
          end do
        end do
      end do
      !--------------------------------------------------------
      ! now sort zratio to find out which species run out first
      !--------------------------------------------------------
      do n = 1 , nqx
        do i = ici1 , ici2
          do j = jci1 , jci2
            iorder(j,i,n) = -999
          end do
        end do
      end do
      do jn = 1 , nqx
        do i = ici1 , ici2
          do j = jci1 , jci2
            llindex1(j,i,jn) = .true.
          end do
        end do
      end do
      do n = 1 , nqx
        do i = ici1 , ici2
          do j = jci1 , jci2
            zmin(j,i) = 1.E32
          end do
        end do
        do jn = 1 , nqx
          do i = ici1 , ici2
            do j = jci1 , jci2
              if ( llindex1(j,i,jn) .and. zratio(j,i,jn) < zmin(j,i) ) then
                iorder(j,i,n) = jn
                zmin(j,i) = zratio(j,i,jn)
              end if
            end do
          end do
        end do
        do i = ici1 , ici2
          do j = jci1 , jci2
            llindex1(j,i,iorder(j,i,n)) = .false.
          end do
        end do
      end do
      !--------------------------------------------
      ! scale the sink terms, in the correct order,
      ! recalculating the scale factor each time
      !--------------------------------------------
      zsinksum(:,:,:) = d_zero
      !----------------
      ! recalculate sum
      !----------------
      do n = 1 , nqx
        do jn = 1 , nqx
          do i = ici1 , ici2
            do j = jci1 , jci2
              jo = iorder(j,i,n)
              llindex3(j,i,jo,jn) = zsolqa(j,i,jo,jn) < d_zero
              zsinksum(j,i,jo) = zsinksum(j,i,jo) - zsolqa(j,i,jo,jn)
            end do
          end do
        end do
        !---------------------------
        ! recalculate scaling factor
        !---------------------------
        do i = ici1 , ici2
          do j = jci1 , jci2
            jo = iorder(j,i,n)
            zratio(j,i,jo) = zqx0(j,i,k,jo)/max(zsinksum(j,i,jo),zqx0(j,i,k,jo))
          end do
        end do
        !------
        ! scale
        !------
        do jn = 1 , nqx
          do i = ici1 , ici2
            do j = jci1 , jci2
              jo = iorder(j,i,n)
              if ( llindex3(j,i,jo,jn) ) then
                zsolqa(j,i,jo,jn) = zsolqa(j,i,jo,jn)*zratio(j,i,jo)
                zsolqa(j,i,jn,jo) = zsolqa(j,i,jn,jo)*zratio(j,i,jo)
              end if
            end do
          end do
        end do
      end do

      ! Set the LHS of equation

      do n = 1 , nqx
        do jn = 1 , nqx
          do i = ici1 , ici2
            do j = jci1 , jci2
              ! Diagonals: microphysical sink terms+transport
              if ( jn == n ) then
                zqlhs(jn,n,j,i) = d_one + ksemi*zfallsink(j,i,n)
                do jo = 1 , nqx
                  zqlhs(jn,n,j,i) = zqlhs(jn,n,j,i) + ksemi*zsolqb(j,i,jo,jn)
                end do
                ! Non-diagonals: microphysical source terms
              else
                ! Here is the delta T - missing from doc.
                zqlhs(jn,n,j,i) = -ksemi*zsolqb(j,i,jn,n)
              end if
            end do
          end do
        end do
      end do

      ! Set the RHS of equation

      do n = 1 , nqx
        do i = ici1 , ici2
          do j = jci1 , jci2
            ! Sum the explicit source and sink
            zexplicit = d_zero
            do jn = 1 , nqx
              ! Positive, since summed over 2nd index
              zexplicit = zexplicit + zsolqa(j,i,n,jn)
              if ( jn /= n ) then
                zexplicit =  zexplicit - &
                        (d_one-ksemi)*zqx0(j,i,k,n)*zsolqb(j,i,jn,n) + &
                        (d_one-ksemi)*zqx0(j,i,k,jn)*zsolqb(j,i,n,jn)
              end if
              zexplicit = zexplicit - &
                    (d_one-ksemi)*zqx0(j,i,k,n)*zfallsink(j,i,n)
            end do
            zqxn(n,j,i) = zqx0(j,i,k,n) + zexplicit
          end do
        end do
      end do
      ! Solve by LU decomposition: dummy test
      ! zqlhs(:,:,:) = 0.0
      ! do jm = 1,nqx
      !   zqlhs(1,jm,jm)=1.0
      ! end do
      ! zqlhs(1,3,3)=zqlhs(1,3,3)+0.47
      ! zqlhs(1,5,3)=zqlhs(1,5,3)-0.47
#ifdef USE_LAPACK
      do i = ici1 , ici2
        do j = jci1 , jci2
          call dgesv(nqx,1,zqlhs(:,:,j,i),nqx,ipivot,zqxn(:,j,i),nqx,ires)
          if ( ires /= 0 ) then
            write(stderr,*) 'Error from lapack subroutine DGESV'
            if ( ires < 0 ) then
              write(stderr,*) 'Argument ',ires,' has illegal value'
            else if ( ires > 0 ) then
              write(stderr,*) 'U(',ires,',',ires,') is zero: Singularity'
            end if
            call fatal(__FILE__,__LINE__,'LAPACK DGESV ERROR')
          end if
        end do
      end do
#else
      ! Internally coded solution
      call mysolve(zqlhs,zqxn)
#endif
      !------------------------------------------------------------------------
      !  Precipitation/sedimentation fluxes to next level
      !     diagnostic precipitation fluxes
      !     It is this scaled flux that must be used for source to next layer
      !------------------------------------------------------------------------
      ! Generalized precipitation flux
      do n = 1 , nqx
        do i = ici1 , ici2
          do j = jci1 , jci2
            ! this will be the source for the k
            zpfplsx(j,i,k+1,n) = ksemi*zfallsink(j,i,n) * &
              zqxn(n,j,i)*zrdtgdp(j,i) + (d_one-ksemi)*zfallsink(j,i,n)* &
              zqx0(j,i,k,n)*zrdtgdp(j,i) ! kg/m2/s
           end do
        end do
      end do
      ! Calculate fluxes in and out of box for conservation of TL
      do n = 1 , nqx
        do i = ici1 , ici2
          do j = jci1 , jci2
            zfluxq(j,i,n) = zconvsrce(j,i,n)+ zfallsrce(j,i,n) - &
                        (zfallsink(j,i,n)+zconvsink(j,i,n))*zqxn(n,j,i)
          end do
        end do
      end do
      ! Calculate the water variables tendencies
      do n = 1 , nqx
        do i = ici1 , ici2
          do j = jci1 , jci2
            zqxtendc(j,i,k,n) = zqxtendc(j,i,k,n) + &
              (zqxn(n,j,i)-zqx0(j,i,k,n))*zqtmst
          end do
        end do
      end do
     ! Calculate the temperature tendencies
      do n = 1 , nqx
        if ( kphase(n) == 1 ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
                zttendc(j,i,k) = zttendc(j,i,k) + &
                                 wlhvocp*(zqxn(n,j,i)-zqx0(j,i,k,n) - &
                                 zfluxq(j,i,n))*zqtmst
            end do
          end do
        else if ( kphase(n) == 2 ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              zttendc(j,i,k) = zttendc(j,i,k) + &
                               wlhsocp*(zqxn(n,j,i)-zqx0(j,i,k,n) - &
                               zfluxq(j,i,n))*zqtmst
            end do
          end do
        end if
      end do
      ! Couple tendencies with pressure
      do n = 1 , nqx
        do i = ici1 , ici2
          do j = jci1 , jci2
            zqxten(j,i,k,n) =  zqxtendc(j,i,k,n)*psb(j,i)
          end do
        end do
      end do
      do i = ici1 , ici2
        do j = jci1 , jci2
          ztten(j,i,k) = zttendc(j,i,k)*psb(j,i)
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
            ztnew = zt(j,i,k)+dt*(zttendc(j,i,k)-ztentkeep(j,i,k))
            if ( k == 1 ) then
              zsumq1(j,i,k) = d_zero ! total water
              zsumh1(j,i,k) = d_zero ! liquid water temperature
            else
              zsumq1(j,i,k) = zsumq1(j,i,k-1)
              zsumh1(j,i,k) = zsumh1(j,i,k-1)
            end if
            ! cld vars
            do n = 1 , nqx
              if ( kphase(n) == 1 ) then
                ztnew = ztnew-wlhvocp*(zqx0(j,i,k,n)+ &
                        (zqxtendc(j,i,k,n)-ztenkeep(j,i,k,n))*dt)
              else if ( kphase(n) == 2 ) then
                ztnew = ztnew-wlhsocp*(zqx0(j,i,k,n)+ &
                        (zqxtendc(j,i,k,n)-ztenkeep(j,i,k,n))*dt)
              end if
              zsumq1(j,i,k) = zsumq1(j,i,k) + &
                (zqx0(j,i,k,n)+(zqxtendc(j,i,k,n)-ztenkeep(j,i,k,n))*dt)* &
                (papf(j,i,k+1)-papf(j,i,k))*regrav
            end do
            zsumh1(j,i,k) = zsumh1(j,i,k)+(papf(j,i,k+1)-papf(j,i,k))*ztnew
            zrain = d_zero
            do n = 1 , nqx
              zrain = zrain+dt*zpfplsx(j,i,k+1,n)
            end do
            zerrorq(j,i,k) = zsumq1(j,i,k)+zrain-zsumq0(j,i,k)
          end do
        end do
      end do
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            zdtgdp(j,i) = dt*egrav/(papf(j,i,k+1)-papf(j,i,k))
            zrain = d_zero
            do n = 1 , nqx
              if ( kphase(n) == 1 ) then
                zrain = zrain+wlhvocp*zdtgdp(j,i)*zpfplsx(j,i,k+1,n)* & !k+1?
                         (papf(j,i,k+1)-papf(j,i,k))
              else if ( kphase(n) == 2 ) then
                zrain = zrain+wlhsocp*zdtgdp(j,i)*zpfplsx(j,i,k+1,n)* &
                        (papf(j,i,k+1)-papf(j,i,k))
              end if
            end do
            zsumh1(j,i,k) = (zsumh1(j,i,k)-zrain)/(papf(j,i,k+1)-papf(j,i,1))
            zerrorh(j,i,k) = zsumh1(j,i,k)-zsumh0(j,i,k)
          end do
        end do
      end do

      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( abs(zerrorq(j,i,kz)) > 1.D-12 .or. &
                 abs(zerrorh(j,i,kz)) > 1.D-12) then
              if ( abs(zerrorq(j,i,kz)) > 1.D-12 ) then
                write(stderr,*) 'WATER NON CONSERVED AT '
                write(stderr,*) 'J = ',j
                write(stderr,*) 'I = ',i
                write(stderr,*) 'K = ',k
                write(stderr,*) 'ERROR IS : ',zerrorq(j,i,kz)
              end if
              if ( abs(zerrorh(j,i,kz)) > 1.D-12 ) then
                write(stderr,*) 'ENTHALPY NON CONSERVED AT '
                write(stderr,*) 'J = ',j
                write(stderr,*) 'I = ',i
                write(stderr,*) 'K = ',k
                write(stderr,*) 'ERROR IS : ',zerrorh(j,i,kz)
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
    zpfplsl(:,:,:) = d_zero
    zpfplsn(:,:,:) = d_zero
    rainls(:,:,:) = d_zero

    !--------------------------------------------------------------------
    ! Copy general precip arrays back into PFP arrays for GRIB archiving
    ! Add rain and liquid fluxes, ice and snow fluxes
    !--------------------------------------------------------------------
    !Rain+liquid, snow+ice
    !for each level k = 1 , kz, sum of the same phase elements
    do n = 1 , nqx
      if ( kphase(n) == 1 ) then
        do k = 1 , kz+1
          do i = ici1 , ici2
            do j = jci1 , jci2
              zpfplsl(j,i,k) = zpfplsl(j,i,k) + zpfplsx(j,i,k,n)
              rainls(j,i,k) = zpfplsl(j,i,k)
            end do
          end do
        end do
      else if ( kphase(n) == 2 ) then
        do k = 1 , kz+1
          do i = ici1 , ici2
            do j = jci1 , jci2
              zpfplsn(j,i,k) = zpfplsn(j,i,k)+ zpfplsx(j,i,k,n)
            end do
          end do
        end do
      end if
    end do

    !--------------------------------------------------------------------
    !Convert the accumlated precipitation to appropriate units for
    !the surface physics and the output sum up through the levels
    !--------------------------------------------------------------------
    ! sum over the levels
    do i = ici1 , ici2
      do j = jci1 , jci2
        prainx = zpfplsl(j,i,kz+1)*dt
        psnowx = zpfplsn(j,i,kz+1)*dt
        if ( prainx > dlowval ) then
          rainnc(j,i) =  rainnc(j,i) + prainx   !mm
          lsmrnc(j,i) =  lsmrnc(j,i) + zpfplsl(j,i,kz+1)
        end if
        if ( psnowx > dlowval ) then
          snownc(j,i) = snownc(j,i) + psnowx
          lsmrnc(j,i) =  lsmrnc(j,i) + zpfplsn(j,i,kz+1)
        end if
      end do
     end do

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif

    contains

      pure real(rk8) function delta(zt)
        !delta = 1 if zt > tzero
        !delta = 0 if zt < tzero
        implicit none
        real(rk8) , intent(in):: zt
        delta = max(d_zero,sign(d_one,zt-tzero))
      end function delta

      pure real(rk8) function phases(zt)
        !phases = 1      if zt > tzero
        !phases = 0      if zt < rtice
        !0<phases < 1    if rtice < zt < tzero
        implicit none
        real(rk8) , intent(in):: zt
        phases = max(min(d_one,((max(rtice,min(tzero,zt))-rtice)* &
                                rtwat_rtice_r)**2),d_zero)
      end function phases

      pure real(rk8) function foedem(zt,phase)
        implicit none
        real(rk8) , intent(in):: zt , phase
        foedem = phase * r5alvcp * (d_one/(zt-r4les)**2) + &
                 (d_one - phase) * r5alscp * (d_one/(zt-r4ies)**2)
      end function foedem

      pure real(rk8) function foeldcpm(zt)
        implicit none
        real(rk8) , intent(in):: zt
        foeldcpm = phases(zt)*wlhvocp+(d_one-phases(zt))*wlhsocp
      end function foeldcpm

      pure real(rk8) function foeeliq(zt) ! = 0.622*esw  !Teten's formula
        implicit none
        real(rk8) , intent(in):: zt
        foeeliq = r2es*exp(r3les*(zt-tzero)/(zt-r4les))
      end function foeeliq

      pure real(rk8) function foeeice(zt) ! = 0.622*esi
        implicit none
        real(rk8) , intent(in):: zt
        foeeice = r2es*exp(r3ies*((zt-tzero)/(zt-r4ies)))
      end function foeeice

      pure real(rk8) function foeewm(zt,phase)
        implicit none
        real(rk8) , intent(in):: zt , phase
        foeewm = phase * foeeliq(zt) + (d_one-phase) * foeeice(zt)
      end function foeewm

      pure real(rk8) function fokoop(zt,foeeliq,foeeice)
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
        ! i.e. fokoop = foeeliq/foeeice, mentre per T < -38 fokoop = RHhomo
        ! while below this threshold the liquid water or aqueous sulphate
        ! solutes are assumed to freeze instantaneously and the process is
        ! a source for cloud ice.
        ! fokoop modifies the ice saturation mixing ratio for homogeneous
        ! nucleation
        real(rk8) , intent(in) :: zt , foeeliq , foeeice
        fokoop = min(rkoop1-rkoop2*zt,foeeliq/foeeice)
      end function fokoop

      subroutine mysolve(aam,bbm)
        implicit none
        real(rk8) , pointer , intent(inout) , dimension(:,:,:,:) :: aam
        real(rk8) , pointer , intent(inout) , dimension(:,:,:) :: bbm
        integer(ik4) :: i , j , k , ii , jj , ll , imax , m , n
        real(rk8) :: aamax , dum , xsum , swap
        integer(ik4) , dimension(nqx) :: indx
        real(rk8) , dimension(nqx) :: vv

        do i = ici1 , ici2
          do j = jci1 , jci2
            !                                                Ux=y
            ! solve A x = b-------------> LU x = b---------> Ly=b
            !
            ! LOOP OVER ROWS TO GET THE IMPLICIT SCALING INFORMATION.
            !
            do m = 1 , nqx
              aamax = d_zero
              do n = 1 , nqx
                if ( abs(aam(m,n,j,i)) > aamax ) aamax = abs(aam(m,n,j,i))
              end do
              if ( aamax < dlowval ) then
                call fatal(__FILE__,__LINE__,'SINGULAR MATRIX')
              end if ! SINGULAR MATRIX
              vv(m) = d_one/aamax ! SAVE THE SCALING.
            end do
            do n = 1 , nqx
              ! THIS IS THE LOOP OVER COLUMNS OF CROUT S METHOD.
              if ( n > 1 ) then
                do m = 1 , n - 1
                  ! THIS IS EQUATION (2.3.12) EXCEPT FOR I = J.
                  xsum = aam(m,n,j,i)
                  do k = 1 , m - 1
                    xsum = xsum - aam(m,k,j,i)*aam(k,n,j,i)
                  end do
                  aam(m,n,j,i) = xsum
                end do
              end if
              ! INITIALIZE FOR THE SEARCH FOR LARGEST PIVOT ELEMENT.
              aamax = d_zero
              imax = n
              ! THIS IS I = J OF EQUATION (2.3.12)
              do m = n , nqx
                ! AND I = J+1. . .N OF EQUATION (2.3.13).
                xsum = aam(m,n,j,i)
                if ( n > 1 ) then
                  do k = 1 , n - 1
                    xsum = xsum - aam(m,k,j,i)*aam(k,n,j,i)
                  end do
                  aam(m,n,j,i) = xsum
                end if
                dum = vv(m)*abs(xsum)   ! FIGURE OF MERIT FOR THE PIVOT.
                if ( dum >= aamax ) then ! IS IT BETTER THAN THE BEST SO FAR?
                  imax = m
                  aamax = dum
                end if
              end do
              if ( n /= imax ) then
                ! DO WE NEED TO INTERCHANGE ROWS? YES, DO SO...
                ! D = -D !...AND CHANGE THE PARITY OF D.
                do ii = 1 , nqx
                  swap = aam(imax,ii,j,i)
                  aam(imax,ii,j,i) = aam(n,ii,j,i)
                  aam(n,ii,j,i) = swap
                end do
                vv(imax) = vv(n) ! ALSO INTERCHANGE THE SCALE FACTOR.
              end if
              indx(n) = imax
              if ( abs(aam(n,n,j,i)) < dlowval ) aam(n,n,j,i) = dlowval
              dum = d_one/aam(n,n,j,i)
              if ( n /= nqx ) then
                do m = n + 1 , nqx
                  aam(m,n,j,i) = aam(m,n,j,i)*dum
                end do
              end if
            end do
            if ( abs(aam(nqx,nqx,j,i)) < dlowval ) then
              aam(nqx,nqx,j,i) = dlowval
            end if
            !
            ! SOLVE THE SET OF N LINEAR EQUATIONS A * X = B.
            ! HERE A IS INPUT, NOT AS THE MATRIX A BUT RATHER AS
            ! ITS LU DECOMPOSITION, DETERMINED BY THE ROUTINE LUDCMP.
            ! INDX IS INPUT AS THE PERMUTATION VECTOR RETURNED BY LUDCMP.
            ! B(1:N) IS INPUT AS THE RIGHT-HAND SIDE VECTOR B,
            ! AND RETURNS WITH THE SOLUTION VECTOR X. A, N, NP,
            ! AND INDX ARE NOT MODIFIED BY THIS ROUTINE AND CAN BE
            ! LEFT IN PLACE FOR SUCCESSIVE CALLS WITH DI
            !
            ii = 0
            ! WHEN II IS SET TO A POSITIVE VALUE, IT WILL BECOME
            ! THE INDEX OF THE  FIRST NONVANISHING ELEMENT OF B.
            ! WE NOW DO THE FORWARD SUBSTITUTION, EQUATION (2.3.6).
            ! THE ONLY NEW WRINKLE IS TO UNSCRAMBLE THE PERMUTATION AS WE GO.
            do m = 1 , nqx
              ll = indx(m)
              xsum = bbm(ll,j,i)
              bbm(ll,j,i) = bbm(m,j,i)
              if ( ii == 0 ) then
                if ( abs(xsum) > dlowval ) ii = m
              else
                do jj = ii , m - 1
                  xsum = xsum - aam(m,jj,j,i)*bbm(jj,j,i)
                end do
              end if
              bbm(m,j,i) = xsum
            end do
            ! NOW WE DO THE BACKSUBSTITUTION, EQUATION (2.3.7)
            do m = nqx , 1 , -1
              xsum = bbm(m,j,i)
              do jj = m + 1 , nqx
                xsum = xsum - aam(m,jj,j,i)*bbm(jj,j,i)
              end do
              ! STORE A COMPONENT OF THE SOLUTION VECTOR X.
              bbm(m,j,i) = xsum/aam(m,m,j,i)
            end do
          end do
        end do
      end subroutine mysolve

  end subroutine microphys

  pure real(rk8) function noscheme_zqe(q,m,f,l) result(zqe)
    implicit none
    real(rk8) , intent(in) :: q , m , f , l
    zqe = max((q-m*f)/(d_one-f),d_zero)
  end function noscheme_zqe

  pure real(rk8) function tompkins_zqe(q,m,f,l) result(zqe)
    implicit none
    real(rk8) , intent(in) :: q , m , f , l
    zqe = max((q-m*f)/(d_one-f),d_zero)
  end function tompkins_zqe

  pure real(rk8) function lohmann_karcher_zqe(q,m,f,l) result(zqe)
    implicit none
    real(rk8) , intent(in) :: q , m , f , l
    zqe = q
  end function lohmann_karcher_zqe

  pure real(rk8) function gierens_zqe(q,m,f,l) result(zqe)
    implicit none
    real(rk8) , intent(in) :: q , m , f , l
    zqe = q/l
  end function gierens_zqe

!  subroutine addpath_array(src,snk,proc2,zsqa,zsqb,beta,fg,j,i)
!    implicit none
!    real(rk8) , pointer , intent(inout) , dimension(:,:,:,:) :: zsqa , zsqb
!    real(rk8) , pointer , intent(inout) , dimension(:,:,:) :: fg
!    real(rk8) , pointer , intent(in) , dimension(:,:)  :: proc2
!    integer(ik4) , intent(in) :: src, snk
!    integer(ik4), intent(in)  :: i , j
!    real(rk8) , intent(in) :: beta
!    zsqa(j,i,src,snk) = zsqa(j,i,src,snk) + (d_one-beta)*proc2(j,i)
!    zsqa(j,i,snk,src) = zsqa(j,i,snk,src) - (d_one-beta)*proc2(j,i)
!    fg(j,i,src) = fg(j,i,src) + (d_one-beta)*proc2(j,i)
!    fg(j,i,snk) = fg(j,i,snk) - (d_one-beta)*proc2(j,i)
!    zsqb(j,i,src,snk) = zsqb(j,i,src,snk) + beta*proc2(j,i)
!  end subroutine addpath_array
!
!  subroutine addpath_real(src,snk,proc,zsqa,zsqb,beta,fg)
!    implicit none
!    real(rk8) , pointer , intent(inout) , dimension(:,:,:,:) :: zsqa , zsqb
!    real(rk8) , pointer , intent(inout) , dimension(:,:,:) :: fg
!    real(rk8) , intent(in) :: proc
!    integer(ik4) , intent(in) :: src , snk
!    integer(ik4) :: i , j
!    real(rk8) , intent(in) :: beta
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

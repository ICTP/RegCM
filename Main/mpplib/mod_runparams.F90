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

module mod_runparams

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_date
  use mod_dynparam
  use mod_mpmessage
  use mod_memutil
  use mod_timer
  use mod_spline

  implicit none

  private

  character(len=256) , public :: namelistfile , prgname

  logical , public :: carb_aging_control = .false.

  integer(ik4) , public :: nqx , iqfrst , iqlst
  integer(ik4) , public , parameter :: iqv = 1
  integer(ik4) , public , parameter :: iqc = 2
  integer(ik4) , public , parameter :: iqr = 3
  integer(ik4) , public , parameter :: iqi = 4
  integer(ik4) , public , parameter :: iqs = 5

  integer(ik4) , public , parameter :: number_of_prognostic_components = 3
  integer(ik4) , public , parameter :: pc_total      = 1
  integer(ik4) , public , parameter :: pc_dynamic    = 2
  integer(ik4) , public , parameter :: pc_physic     = 3

  type(rcm_time_and_date) , save , public :: idate0 , idate1 , idate2

  type(rcm_timer) , save , public , pointer :: rcmtimer

  type(rcm_alarm) , save , public , pointer :: alarm_hour
  type(rcm_alarm) , save , public , pointer :: alarm_day
  type(rcm_alarm) , save , public , pointer :: alarm_out_nwf
  type(rcm_alarm) , save , public , pointer :: alarm_out_sav
  type(rcm_alarm) , save , public , pointer :: alarm_out_atm
  type(rcm_alarm) , save , public , pointer :: alarm_out_rad
  type(rcm_alarm) , save , public , pointer :: alarm_out_srf
  type(rcm_alarm) , save , public , pointer :: alarm_out_shf
  type(rcm_alarm) , save , public , pointer :: alarm_out_sts
  type(rcm_alarm) , save , public , pointer :: alarm_out_che
  type(rcm_alarm) , save , public , pointer :: alarm_out_lak
  type(rcm_alarm) , save , public , pointer :: alarm_out_opt
  type(rcm_alarm) , save , public , pointer :: alarm_out_sub

  type(rcm_alarm) , save , public , pointer :: alarm_in_bdy

  type(rcm_syncro) , save , public , pointer :: syncro_dbg
  type(rcm_syncro) , save , public , pointer :: syncro_rep
  type(rcm_syncro) , save , public , pointer :: syncro_srf
  type(rcm_syncro) , save , public , pointer :: syncro_rad
  type(rcm_syncro) , save , public , pointer :: syncro_radfor
  type(rcm_syncro) , save , public , pointer :: syncro_emi
  type(rcm_syncro) , save , public , pointer :: syncro_cum
  type(rcm_syncro) , save , public , pointer :: syncro_che

  type(rcm_syncro) , save , public , pointer :: syncro_cpl

  ! Orbital paramters
  real(rk8) , public :: year_offset
  real(rk8) , public :: eccen
  real(rk8) , public :: obliq
  real(rk8) , public :: mvelp
  real(rk8) , public :: obliqr
  real(rk8) , public :: lambm0
  real(rk8) , public :: mvelpp
  real(rk8) , public :: eccf

  type(rcm_time_and_date) , save , public :: bdydate1 , bdydate2
  type(rcm_time_and_date) , save , public :: somdate1 , somdate2

  type(rcm_time_interval) , save , public :: intbdy
  type(rcm_time_interval) , save , public :: intsom

  real(rk8) , public :: declin
  real(rkx) , public :: xbctime
  real(rkx) , public :: calday , twodt
  real(rk8) , public :: xslabtime

  real(rkx) , public :: solcon , scon

  logical , public :: uvrotate

  ! Step counters to activate surface and radiation schemes
  real(rkx) , public :: rnsrf_for_srffrq , rnsrf_for_day , &
     rnsrf_for_lakfrq , rnsrf_for_subfrq , rnrad_for_optfrq , &
     rnrad_for_srffrq , rnrad_for_radfrq
  ! Step of surface scheme in one atmosphere I/O interval
  real(rkx) , public :: rnsrf_for_atmfrq
  ! One over seconds in one surface I/O interval
  real(rkx) , public :: rsrffrq_sec
  ! Model base timestep in seconds
  real(rkx) , public :: dtsec
  !
  ! Run with idealized conditions
  integer(ik4) , public :: irceideal
  ! Cumulus scheme index
  integer(ik4) , public :: icup_lnd
  integer(ik4) , public :: icup_ocn
  integer(ik4) , dimension(2) , public :: icup
  ! Closure index for Grell
  integer(ik4) , public :: igcc
  ! Cumulus cloud model
  integer(ik4) , public :: icumcloud
  ! Boundary layer index
  integer(ik4) , public :: ibltyp
  ! Diffusion scheme
  integer(ik4) , public :: idiffu
  integer , public :: idif
  ! Lake model activation index
  integer(ik4) , public :: lakemod
  ! Diurnal cycle SST index
  integer(ik4) , public :: idcsst
  integer(ik4) , public :: ipcpcool
  integer(ik4) , public :: iwhitecap
  ! Sea Ice scheme index
  integer(ik4) , public :: iseaice
  real(rkx) , public :: icetriggert
  ! Seasonal albedo for desert index
  integer(ik4) , public :: idesseas
  ! Ocean model switch indexes
  integer(ik4) , public :: iocnrough , iocnflx , iocncpl , iocnzoq
  ! Wave model switch indexes
  integer(ik4) , public :: iwavcpl
  ! COP switch indexes
  integer(ik4) , public :: icopcpl
  ! Radiation switch controls
  integer(ik4) , public :: idirect , iindirect , iemiss , isolconst , ifixsolar
  integer(ik4) , public :: isnowdark
  integer(ik4) , public :: ichecold
  ! Fixed solar constant for ifixsolar = 1
  real(rkx) , public :: fixedsolarval
  ! Semi-Langrangian advection for tracers
  integer(ik4) , public :: isladvec
  integer(ik4) , public :: iqmsl
  ! Use for large scale lwp the same algo used for convective lwp
  integer(ik4) , public :: iconvlwp
  ! Cloud fraction scheme
  integer(ik4) , public :: icldfrac
  ! Marine stratocumulus considered
  integer(ik4) , public :: icldmstrat
  ! vqrang is the range limit on vqflx.
  real(rkx) , public :: vqrang = 5.0e-4_rkx
  ! Aerosol effects
  real(rkx) , public :: rcrit
  real(rkx) , public :: coef_ccn
  real(rkx) , public :: abulk

  character(len=8) , public :: scenario
  integer(ik4) , public :: ghg_year_const

  ! Moloch
  logical , parameter , public :: moloch_do_test_1 = .false.
  logical , parameter , public :: moloch_do_test_2 = .false.
  real(rkx) , public :: mo_dzita , mo_anu2
  real(rkx) , public :: mo_wmax , mo_cflhmax , mo_cflsmax
  logical , public :: mo_filterpai
  integer(ik4) , public :: mo_nzfilt
  integer(ik4) , public :: mo_nadv
  integer(ik4) , public :: mo_nsound

  real(rkx) , public :: dt , dt2 , dtsq , dtcb , dtbdys , rdt
  real(rkx) , public :: dx , dx2 , dx4 , dx8 , dx16 , dxsq
  real(rkx) , public :: rdx , rdxsq
  real(rkx) , public :: dtsrf , dtabem , dtrad , dtcum , dtche
  real(rkx) , public :: cpldt , zomax , ustarmax
  real(rkx) , public :: ckh
  integer(ik4) , public :: diffu_hgtf

  ! Set this to zero to remove dynamical dependency of diffusion
  real(rkx) , public :: adyndif = d_one

#ifdef CLM45
  integer(ik4) , public :: sfbcread
  logical , public :: rcm_megan_enabled
#endif
  integer(ik4) , public :: iboudy , ichem , ipgf , ipptls
  ! usefull flags for chemistry
  integer(ik4) , public :: iaerosol , igaschem , ioxclim , iisoropia
  character(len=6) , pointer , dimension(:) , public :: chtrname
  integer(ik4) , public :: nbin = 1
  integer(ik4) , public :: nmine = 1

  logical , public :: do_parallel_netcdf_in , do_parallel_netcdf_out
  logical , public :: do_parallel_save
  logical , public :: ifrest , doing_restart , lsync , chechgact

  real(rkx) , pointer , dimension(:) , public :: dtau , dtsplit
  real(rkx) , pointer , dimension(:) , public :: hsigma , dsigma , qcon
  real(rkx) , pointer , dimension(:) , public :: sigma , zita , zitah
  real(rkx) , pointer , dimension(:) , public :: ffilt , ak , bk
  real(rkx) , pointer , dimension(:,:) , public :: twt

  real(rkx) , public :: clfrcv ! Cloud fractional cover for convective precip

  ! Moisture from previous run

  logical , public :: replacemoist = .false.

  ! Number od split exp modes

  integer(ik4) , public :: nsplit
  logical , public :: lstand

  ! Asselin filter parameter

  real(rkx) , public :: gnu1
  real(rkx) , public :: gnu2

  ! Advection tuning

  logical , public :: upstream_mode
  real(rkx) , public :: uoffc

  logical , public :: stability_enhance
  real(rkx) , public :: t_extrema
  real(rkx) , public :: q_rel_extrema

  ! Non hydrostatic core parameters

  ! Upper radiative BC for non-hydrostatic core
  real(rkx) , public :: base_state_ts0
  integer(ik4) , public :: ifupr
  ! BET parameter in sound waves removal
  real(rkx) , public :: nhbet
  ! XKD parameter in sound waves removal
  real(rkx) , public :: nhxkd
  ! Raleigh dumping activation
  integer(ik4) , public :: ifrayd
  integer(ik4) , public :: itopnudge
  integer(ik4) , public :: rayndamp
  real(rkx) , public :: rayalpha0
  real(rkx) , public :: rayzd
  real(rkx) , public :: rayhd

  ! Grell cumulus scheme parameters

  real(rkx) , public :: mincld
  real(rkx) , public :: skbmax
  real(rkx) , public :: shrmax_ocn
  real(rkx) , public :: shrmin_ocn
  real(rkx) , public :: edtmax_ocn
  real(rkx) , public :: edtmin_ocn
  real(rkx) , public :: edtmaxo_ocn
  real(rkx) , public :: edtmino_ocn
  real(rkx) , public :: edtmaxx_ocn
  real(rkx) , public :: edtminx_ocn
  real(rkx) , public :: shrmax
  real(rkx) , public :: shrmin
  real(rkx) , public :: edtmax
  real(rkx) , public :: edtmin
  real(rkx) , public :: edtmaxo
  real(rkx) , public :: edtmino
  real(rkx) , public :: edtmaxx
  real(rkx) , public :: edtminx
  real(rkx) , public :: dtauc
  real(rkx) , public :: pbcmax
  real(rkx) , public :: htmax
  real(rkx) , public :: htmin

  ! Emanuel MIT cumulus scheme parameters

  real(rkx) , public :: alphae
  real(rkx) , public :: betae
  real(rkx) , public :: coeffr
  real(rkx) , public :: coeffs
  real(rkx) , public :: cu
  real(rkx) , public :: damp
  real(rkx) , public :: dtmax
  real(rkx) , public :: elcrit_ocn
  real(rkx) , public :: elcrit_lnd
  real(rkx) , public :: entp
  real(rkx) , public :: omtrain
  real(rkx) , public :: omtsnow
  real(rkx) , public :: sigd
  real(rkx) , public :: sigs
  real(rkx) , public :: tlcrit
  real(rkx) , public :: epmax_ocn
  real(rkx) , public :: epmax_lnd
  integer(ik4) , public :: minorig

  ! Tiedtke cumulus scheme parameters

  logical , public :: lmfpen    = .true.  ! penetrative conv is switched on
  logical , public :: lmfmid    = .true.  ! midlevel conv is switched on
  logical , public :: lmfdd     = .true.  ! cumulus downdraft is switched on
  logical , public :: lepcld    = .true.  ! prognostic cloud scheme is on
  logical , public :: lmfdudv   = .true.  ! cumulus friction is switched on
  logical , public :: lmfscv    = .true.  ! shallow convection is switched on
  logical , public :: lmfuvdis  = .true.  ! use kinetic energy dissipation
  logical , public :: lmftrac   = .true.  ! chemical tracer transport is on
  logical , public :: lmfsmooth = .false. ! smoot of mass fluxes for tracers
  logical , public :: lmfwstar  = .false. ! Grant w* closure for shallow conv

  integer(ik4) , public :: iconv

  real(rkx) , public :: entrdd  ! entrainment rate for cumulus downdrafts

  ! ICONV 1, 2, 3
  real(rkx) , public :: entrpen_lnd ! entrainment rate for penetrative conv.
  real(rkx) , public :: entrpen_ocn ! entrainment rate for penetrative conv.
  real(rkx) , public :: entrscv ! entrainment rate for shallow convection
  real(rkx) , public :: entrmid ! entrainment rate for midlevel convection
  real(rkx) , public :: cprcon  ! coefficients for determining conversion
                                ! from cloud water to rain
  real(rkx) , public :: entrmax ! Max entrainment

  ! ICONV 4 paramters
  real(rkx) , public :: cmtcape   ! CAPE removal timestep
  real(rkx) , public :: rcuc_lnd  ! Convective cloud cover for rain evporation
  real(rkx) , public :: rcuc_ocn  ! Convective cloud cover for rain evporation
  real(rkx) , public :: rcpec_lnd ! Coefficient for rain evaporation below cloud
  real(rkx) , public :: rcpec_ocn ! Coefficient for rain evaporation below cloud
  real(rkx) , public :: rhebc_lnd ! Critical relative humidity below
                                  ! cloud at which evaporation starts for land
  real(rkx) , public :: rhebc_ocn ! Critical relative humidity below
                                  ! cloud at which evaporation starts for ocean
  real(rkx) , public :: rprc_lnd  ! coeff for conversion from cloud water
  real(rkx) , public :: rprc_ocn  ! coeff for conversion from cloud water
  real(rkx) , public :: detrpen_lnd ! Detrainment rate for penetrative conv
  real(rkx) , public :: detrpen_ocn ! Detrainment rate for penetrative conv
  real(rkx) , public :: entshalp  ! shallow entrainment factor for entrorg

  ! Kain-Fritsch parameter

  real(rkx) , public :: kf_entrate
  real(rkx) , public :: kf_convrate
  real(rkx) , public :: kf_min_pef
  real(rkx) , public :: kf_max_pef
  real(rkx) , public :: kf_dpp
  real(rkx) , public :: kf_min_dtcape
  real(rkx) , public :: kf_max_dtcape
  real(rkx) , public :: kf_tkemax
  real(rkx) , public :: kf_wthreshold

  ! Tweak Global data

  integer(ik4) , public :: itweak
  integer(ik4) , public :: itweak_sst
  integer(ik4) , public :: itweak_temperature
  integer(ik4) , public :: itweak_solar_irradiance
  integer(ik4) , public :: itweak_greenhouse_gases
  real(rkx) , public :: sst_tweak
  real(rkx) , public :: temperature_tweak
  real(rkx) , public :: solar_tweak
  real(rkx) , dimension(5) , public :: gas_tweak_factors

  ! RRTM scheme parameters

  integer(ik4) , public :: irrtm
  integer(ik4) , public :: irrtm_cldov
  integer(ik4) , public :: irrtm_sw_opcliq
  integer(ik4) , public :: irrtm_sw_opcice
  integer(ik4) , public :: inflgsw
  integer(ik4) , public :: iceflgsw
  integer(ik4) , public :: liqflgsw
  integer(ik4) , public :: icld
  integer(ik4) , public :: irng
  integer(ik4) , public :: imcica
  integer(ik4) , public :: inflglw
  integer(ik4) , public :: iceflglw
  integer(ik4) , public :: liqflglw
  integer(ik4) , public :: nradfo

  ! Radiation schemes common parametrs

  integer(ik4) , public :: iclimao3
  integer(ik4) , public :: iclimaaer
  character(len=256) , public :: radclimpath

  ! UW PBL parameters

  integer(ik4) , public :: iuwvadv
  real(rkx) , public :: rstbl
  real(rkx) , public :: atwo
  real(rkx) , public  :: czero
  real(rkx) , public  :: nuk

  ! Holtslag PBL parameters

  real(rkx) , public :: ricr_ocn
  real(rkx) , public :: ricr_lnd
  real(rkx) , public :: zhnew_fac
  integer(ik4) , public :: ifaholtth10
  integer(ik4) , public :: ifaholt
  integer(ik4) , public :: holtth10iter

  ! Chemistry nameliste option

  character(len=8) , public :: chemsimtype
  integer(ik4) , public :: ichcumtra
  integer(ik4) , public :: ichdrdepo
  integer(ik4) , public :: ichremcvc
  integer(ik4) , public :: ichremlsc
  integer(ik4) , public :: ichsursrc
  integer(ik4) , public :: ichsolver
  integer(ik4) , public :: ismoke ! marijuana
  integer(ik4) , public :: ichdustemd
  integer(ik4) , public :: ichdustparam
  integer(ik4) , public :: ichdiag
  integer(ik4) , public :: ichebdy
  integer(ik4) , public :: ichjphcld
  integer(ik4) , public :: ichbion
  integer(ik4) , public :: ichlinox
  real(rkx) , public :: rdstemfac
  real(rkx) , public :: rocemfac

  ! chemistry species indices that are used not only in chemlib but also in
  ! other interface ( e.g CLM4.5)/ other species are delcared in
  ! chemlib/mod_che_indices
  integer(ik4) , public , parameter :: nchlmax = 3
  integer(ik4) , public :: nochl , nbchl
  integer(ik4) , public , dimension(nchlmax) :: iochl , ibchl
  integer(ik4) , public :: ibchb , iochb , ianh4 , iano3 , iisop , &
                           ich4 , ism1 , ism2 , ino

  ! Cloud control parameters

  integer(ik4) , public :: ncld   ! # of bottom model levels with no clouds
  logical , public :: lsrfhack    ! Surface radiation hack
  logical , public :: larcticcorr ! Vavrus and Valiser Arctic correction
  real(rkx) , public :: rhmin     ! Minimum value for relative humidity
  real(rkx) , public :: rhmax     ! Maximum value for relative humidity
  real(rkx) , public :: rh0land   ! Trigger value for RH for land (SUBEX)
  real(rkx) , public :: rh0oce    ! Trigger value for RH for ocean (SUBEX)
  real(rkx) , public :: tc0       ! Temperature at which cloud are high (SUBEX)
  real(rkx) , public :: cllwcv    ! Cloud liquid water content (Kuo, Grell)
  real(rkx) , public :: clfrcvmax ! Maximum non-scaled convetive cloud fraction
  real(rkx) , public :: cftotmax  ! Maximum total cloud fraction for radiation
  real(rkx) , public :: kfac_shal ! Convective cloud fraction factor for
                                  ! empirical relationship with mass flux
  real(rkx) , public :: kfac_deep ! Convective cloud fraction factor for
                                  ! empirical relationship with mass flux
  real(rkx) , public :: k2_const  ! Convective cloud fraction factor for
                                  ! empirical relationship with mass flux

  ! Large scale SUBEX parameters

  real(rkx) , public :: qck1land
  real(rkx) , public :: qck1oce
  real(rkx) , public :: gulland
  real(rkx) , public :: guloce
  real(rkx) , public :: caccrlnd
  real(rkx) , public :: cevaplnd
  real(rkx) , public :: caccroce
  real(rkx) , public :: cevapoce
  real(rkx) , public :: conf     ! Condensation efficiency.

  !
  ! SLAB ocean parameters
  !
  integer(ik4) , public :: islab_ocean
  logical , public :: ifslaboc = .false.
  logical , public :: do_qflux_adj
  logical , public :: do_restore_sst
  ! TK mod; restoring time scale for SST in days
  real(rkx) , public :: sst_restore_timescale
  real(rkx) , public :: mixed_layer_depth
  integer(ik4) , dimension(mpy) , public :: stepcount

#ifdef CLM
  ! CLM options
  integer(ik4) , public :: imask
  integer(ik4) , public :: ilawrence_albedo
  real(rkx) , public :: clmfrq
  character(len=256) , public :: dirclm
#endif

  !
  ! New microphys parameters
  !
  ! Total water and enthalpy budget on/off
  logical , public :: stats
  logical , public :: budget_compute
  ! Super saturation option
  integer , public :: nssopt
  ! Choose the autoconversion paramaterization
  integer , public :: iautoconv
  ! Fall speed values
  real(rkx) , public :: vfqr
  real(rkx) , public :: vfqi
  real(rkx) , public :: vfqs
  ! autoconversion values
  real(rkx) , public :: auto_rate_khair
  real(rkx) , public :: auto_rate_kessl
  real(rkx) , public :: auto_rate_klepi
  real(rkx) , public :: rkconv
  real(rkx) , public :: skconv
  real(rkx) , public :: rcldiff
  ! limit values
  real(rkx) , public :: rcovpmin
  ! evaporation values
  real(rkx) , public :: rpecons
  !
  ! Debug limit values for tendencies
  !
  real(rkx) , public :: temp_tend_maxval
  real(rkx) , public :: wind_tend_maxval

  ! option for writing tendency diagnostic
  integer(ik4) , public :: idiag
  integer(ik4) , public :: icosp

  data doing_restart /.false./

  public :: allocate_mod_runparams
  public :: iswater , isocean , islake
  public :: exponential_nudging

  contains

  subroutine allocate_mod_runparams
    implicit none
    if ( idynamic < 3 ) then
      if ( idynamic == 1 ) then
        call getmem1d(dtau,1,nsplit,'mod_runparams:dtau')
        call getmem1d(dtsplit,1,nsplit,'mod_runparams:dtsplit')
      end if
    else
      call getmem1d(zita,1,kzp1,'mod_runparams:zita')
      call getmem1d(zitah,1,kz,'mod_runparams:zitah')
      call getmem1d(ffilt,1,kz,'mod_runparams:ffilt')
      call getmem1d(ak,1,kz,'mod_runparams:ak')
      call getmem1d(bk,1,kz,'mod_runparams:bk')
    end if
    call getmem1d(dsigma,1,kz,'mod_runparams:dsigma')
    call getmem1d(hsigma,1,kz,'mod_runparams:hsigma')
    call getmem1d(sigma,1,kzp1,'mod_runparams:sigma')
    call getmem1d(qcon,1,kz,'mod_runparams:qcon')
    call getmem2d(twt,1,kz,1,2,'mod_runparams:twt')
  end subroutine allocate_mod_runparams

  logical function iswater(a)
    real(rkx) , intent(in) :: a
    iswater = .false.
    if (a > 13.5_rkx .and. a < 15.5_rkx) iswater = .true.
  end function

  logical function isocean(a)
    real(rkx) , intent(in) :: a
    isocean = .false.
    if (a > 14.5_rkx .and. a < 15.5_rkx) isocean = .true.
  end function

  logical function islake(a)
    real(rkx) , intent(in) :: a
    islake = .false.
    if (a > 13.5_rkx .and. a < 14.5_rkx) islake = .true.
  end function

  subroutine exponential_nudging(nudge)
    implicit none
    real(rkx) , dimension(kz) , intent(out) :: nudge
    real(rkx) , dimension(3) :: ncin
    real(rkx) , dimension(3) :: zcin
    real(rkx) , dimension(3) :: ycin
    data ycin /0.0_rkx, 0.0_rkx, 0.0_rkx/
    ncin(1) = high_nudge
    ncin(2) = medium_nudge
    ncin(3) = low_nudge
    if ( idynamic == 3 ) then
      zcin(1) = 0.0_rkx
      zcin(2) = 0.5_rkx
      zcin(3) = 1.0_rkx
    else
      zcin(1) = sigma(1)
      zcin(2) = (sigma(kzp1)-sigma(1))*0.5_rkx
      zcin(3) = sigma(kzp1)
    end if
    call spline1d(3,zcin,ncin,ycin,kz,hsigma,nudge)
    if ( myid == 0 ) then
      call vprntv(nudge,kz,'Nudging coefficient profile')
    end if
  end subroutine exponential_nudging

end module mod_runparams
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

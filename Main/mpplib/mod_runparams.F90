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
  use mod_date
  use mod_dynparam
  use mod_memutil

  implicit none

  private

  character(len=256) , public :: namelistfile , prgname

  integer(ik4) , public :: nqx , iqfrst , iqlst
  integer(ik4) , public , parameter :: iqv = 1
  integer(ik4) , public , parameter :: iqc = 2
  integer(ik4) , public , parameter :: iqi = 3
  integer(ik4) , public , parameter :: iqr = 4
  integer(ik4) , public , parameter :: iqs = 5

  type(rcm_time_and_date) , save , public :: idate0 , idate1 , idate2

  type(rcm_time_and_date) , save , public :: idatex
  integer(ik4) , public :: xyear , xmonth , xday , xhour

  ! Orbital paramters
  real(rkx) , public :: eccen
  real(rkx) , public :: obliqr
  real(rkx) , public :: lambm0
  real(rkx) , public :: mvelpp
  real(rkx) , public :: eccf

  type(rcm_time_and_date) , save , public :: bdydate1 , bdydate2
  type(rcm_time_and_date) , save , public :: somdate1 , somdate2

  type(rcm_time_interval) , save , public :: intmdl
  type(rcm_time_interval) , save , public :: intbdy
  type(rcm_time_interval) , save , public :: intsom

  real(rkx) , public :: declin , deltmx
  real(rkx) , public :: xbctime
  real(rkx) , public :: calday , twodt
  real(rk8) , public :: xslabtime

  real(rkx) , public :: solcon , scon

  ! Step counter. Is zero at idate0, always increasing, never reset.
  integer(ik8) , public :: ktau
  ! Final number of step for THIS run
  integer(ik8) , public :: mtau
  ! How many steps for an hour (updates date fields Y m d H)
  integer(ik8) , public :: khour
  ! How many steps for a day (updates date fields Y m d)
  integer(ik8) , public :: kday
  ! Output k values for I/O operations.
  integer(ik8) , public :: katm , krad , kche , ksav , kdbg , kbdy , &
                  ksrf , ksub , klak , krep
  ! Seconds counter in between boundary conditions read
  integer(ik8) , public :: nbdytime
  ! Step counters to activate surface and radiation schemes
  integer(ik8) , public :: ntsrf , ntrad , ntcum , ntabem , ntche , ntcpl
  real(rkx) , public :: rtsrf , rtrad , rnsrf_for_srffrq , rnsrf_for_day , &
               rnsrf_for_lakfrq , rnsrf_for_subfrq , rnrad_for_chem , &
               rnrad_for_radfrq
  real(rkx) , public :: afdout , cfdout
  ! Step of surface scheme in one atmosphere I/O interval
  real(rkx) , public :: rsrf_in_atm
  ! One over seconds in one surface I/O interval
  real(rkx) , public :: rsrffrq_sec
  ! Model timestep in seconds (real and integer)
  integer(ik8) , public :: ntsec
  real(rkx) , public :: dtsec
  ! Internal count for how many SRF outputs per day
  integer(ik8) , public :: ksts , kstsoff
  !
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
  ! Lake model activation index
  integer(ik4) , public :: lakemod
  ! Diurnal cycle SST index
  integer(ik4) , public :: idcsst
  ! Sea Ice scheme index
  integer(ik4) , public :: iseaice
  ! Seasonal albedo for desert index
  integer(ik4) , public :: idesseas
  ! Ocean model switch indexes
  integer(ik4) , public :: iocnrough , iocnflx , iocncpl , iocnzoq
  ! Wave model switch indexes
  integer(ik4) , public :: iwavcpl
  ! Radiation switch controls
  integer(ik4) , public :: idirect , iindirect , iemiss , isolconst
  ! Semi-Langrangian advection for tracers
  integer(ik4) , public :: isladvec
  integer(ik4) , public :: iqmsl
  ! Use for large scale lwp the same algo used for convective lwp
  integer(ik4) , public :: iconvlwp
  ! Cloud fraction scheme
  integer(ik4) , public :: icldfrac
  ! Marine stratocumulus considered
  integer(ik4) , public :: icldmstrat
  ! Upper radiative BC for non-hydrostatic core
  integer(ik4) , public :: ifupr
  ! BET parameter in sound waves removal
  real(rkx) , public :: nhbet
  ! XKD parameter in sound waves removal
  real(rkx) , public :: nhxkd
  ! vqrang is the range limit on vqflx.
  real(rkx) , public :: vqrang = 5.0e-4_rkx
  ! Aerosol effects
  real(rkx) , public :: rcrit
  real(rkx) , public :: coef_ccn
  real(rkx) , public :: abulk

  character(len=8) , public :: scenario

  real(rkx) , public :: dt , dt2 , dtsq , dtcb , dtbdys , rdt
  real(rkx) , public :: dx , dx2 , dx4 , dx8 , dx16 , dxsq
  real(rkx) , public :: rdxsq
  real(rkx) , public :: dtsrf , dtabem , dtrad , dtcum , dtche
  real(rkx) , public :: cpldt , zomax , ustarmax
  real(rkx) , public :: ckh
  integer(ik4) , public :: diffu_hgtf

  integer(ik4) , public :: iboudy , ichem , ipgf , ipptls
  ! usefull flags for chemistry
  integer(ik4) , public :: iaerosol , igaschem , ioxclim , iisoropia
  character(len=6) , pointer , dimension(:) , public :: chtrname
  integer(ik4) , public :: nbin = 1
  integer(ik4) , public :: nmine = 1

  logical , public :: do_parallel_netcdf_in , do_parallel_netcdf_out
  logical , public :: ifrest , doing_restart , lsync

  integer(ik4) , public :: kchi , kclo , kcmd

  real(rkx) , pointer , dimension(:) , public :: dtau
  real(rkx) , pointer , dimension(:) , public :: hsigma , dsigma , qcon
  real(rkx) , pointer , dimension(:) , public :: sigma
  real(rkx) , pointer , dimension(:) , public :: anudgh , anudgf
  real(rkx) , pointer , dimension(:,:) , public :: twt

  real(rkx) , public :: clfrcv ! Cloud fractional cover for convective precip

  ! Moisture from previous run

  logical , public :: replacemoist = .false.

  ! Non hydrostatic core parameters

  real(rkx) , public :: logp_lrate

  ! Grell cumulus scheme parameters

  real(rkx) , public :: gcr0
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
  real(rkx) , public :: minsig
  real(rkx) , public :: omtrain
  real(rkx) , public :: omtsnow
  real(rkx) , public :: sigd
  real(rkx) , public :: sigs
  real(rkx) , public :: tlcrit
  real(rkx) , public :: epmax_ocn
  real(rkx) , public :: epmax_lnd
  integer(ik4) , public :: minorig

  ! Tiedtke cumulus scheme parameters

  integer(ik4) , public :: iconv

  real(rkx) , public :: entrdd  ! entrainment rate for cumulus downdrafts

  ! ICONV 1, 2, 3
  real(rkx) , public :: entrpen ! entrainment rate for penetrative convection
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
  real(rkx) , public :: detrpen   ! Detrainment rate for penetrative convection
  real(rkx) , public :: entshalp  ! shallow entrainment factor for entrorg

  ! Kain-Fritsch parameter

  real(rkx) , public :: kf_entrate
  real(rkx) , public :: kf_min_pef
  real(rkx) , public :: kf_max_pef
  real(rkx) , public :: kf_dpp
  real(rkx) , public :: kf_min_dtcape
  real(rkx) , public :: kf_max_dtcape
  real(rkx) , public :: kf_tkemax
  !

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
  integer(ik4) , public :: ichdiag
  integer(ik4) , public :: ichebdy
  integer(ik4) , public :: ichjphcld
  integer(ik4) , public :: ichbion
  real(rkx) , public :: rdstemfac

  ! chemistry species indices that are used not only in chemlib but also in
  ! other interface ( e.g CLM4.5)/ other species are delcared in
  ! chemlib/mod_che_indices
  integer(ik4) , public :: ibchb , ibchl , iochl , iochb , ianh4 , &
                           iano3 , iisop , ich4 , ism1 , ism2 , ino

  ! Large scale SUBEX parameters

  integer(ik4) , public :: ncld ! # of bottom model levels with no clouds
  logical , public :: lsrfhack ! Surface radiation hack
  real(rkx) , public :: qck1land
  real(rkx) , public :: qck1oce
  real(rkx) , public :: gulland
  real(rkx) , public :: guloce
  real(rkx) , public :: rhmin
  real(rkx) , public :: rhmax
  real(rkx) , public :: rh0land
  real(rkx) , public :: rh0oce
  real(rkx) , public :: caccrlnd
  real(rkx) , public :: cevaplnd
  real(rkx) , public :: caccroce
  real(rkx) , public :: cevapoce
  real(rkx) , public :: tc0
  real(rkx) , public :: cllwcv   ! Cloud liquid water content
                                 ! for convective precip.
  real(rkx) , public :: clfrcvmax
  real(rkx) , public :: cftotmax ! Maximum total cloud fraction for radiation
  real(rkx) , public :: conf     ! Condensation threshold.

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
  ! Implict option
  real(rkx) , public :: rsemi
  ! Fall speed values
  real(rkx) , public :: vfqr
  real(rkx) , public :: vfqi
  real(rkx) , public :: vfqs
  ! autoconversion values
  real(rkx) , public :: auto_rate_khair
  real(rkx) , public :: auto_rate_kessl
  real(rkx) , public :: auto_rate_klepi
  real(rkx) , public :: rkconv
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

  public :: allocate_mod_runparams , iswater , isocean , islake

  contains

  subroutine allocate_mod_runparams
    implicit none
    call getmem1d(hsigma,1,kz,'mod_runparams:hsigma')
    call getmem1d(dsigma,1,kz,'mod_runparams:dsigma')
    call getmem1d(qcon,1,kz,'mod_runparams:qcon')
    call getmem1d(sigma,1,kzp1,'mod_runparams:sigma')
    call getmem2d(twt,1,kz,1,2,'mod_runparams:twt')
    call getmem1d(anudgh,1,kz,'mod_runparams:anudgh')
    call getmem1d(anudgf,1,kzp1,'mod_runparams:anudgf')
    call getmem1d(dtau,1,nsplit,'mod_runparams:nsplit')
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

end module mod_runparams
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

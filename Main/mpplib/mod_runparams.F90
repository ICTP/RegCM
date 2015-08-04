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
  integer(ik4) , public , parameter :: iqr = 3
  integer(ik4) , public , parameter :: iqi = 4
  integer(ik4) , public , parameter :: iqs = 5

  type(rcm_time_and_date) , save , public :: idate0 , idate1 , idate2

  type(rcm_time_and_date) , save , public :: idatex
  integer(ik4) , public :: xyear , xmonth , xday , xhour

  ! Orbital paramters
  real(rk8) , public :: eccen
  real(rk8) , public :: obliqr
  real(rk8) , public :: lambm0
  real(rk8) , public :: mvelpp
  real(rk8) , public :: eccf

  type(rcm_time_and_date) , save , public :: bdydate1 , bdydate2
  type(rcm_time_and_date) , save , public :: somdate1 , somdate2

  type(rcm_time_interval) , save , public :: intmdl
  type(rcm_time_interval) , save , public :: intbdy
  type(rcm_time_interval) , save , public :: intsom

  real(rk8) , public :: declin , deltmx
  real(rk8) , public :: xbctime , xslabtime
  real(rk8) , public :: calday , twodt

  real(rk8) , public :: solcon , scon

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
  integer(ik8) , public :: ntsrf , ntrad , ntabem , ntcpl
  real(rk8) , public :: rtsrf , rtrad , rnsrf_for_srffrq , rnsrf_for_day , &
               rnsrf_for_lakfrq , rnsrf_for_subfrq , rnrad_for_chem , &
               rnrad_for_radfrq
  real(rk8) , public :: afdout , cfdout
  ! Step of surface scheme in one atmosphere I/O interval
  real(rk8) , public :: rsrf_in_atm
  ! One over seconds in one surface I/O interval
  real(rk8) , public :: rsrffrq_sec
  ! Model timestep in seconds (real and integer)
  integer(ik8) , public :: ntsec
  real(rk8) , public :: dtsec
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
  ! Radiation switch controls
  integer(ik4) , public :: idirect , iindirect , iemiss , isolconst
  ! Semi-Langrangian advection for tracers
  integer(ik4) , public :: isladvec
  ! Convective LWP scheme
  integer(ik4) , public :: iconvlwp
  ! Upper radiative BC for non-hydrostatic core
  integer(ik4) , public :: ifupr

  character(len=8) , public :: scenario

  real(rk8) , public :: dt , dt2 , dtsq , dtcb , dtbdys , rdt
  real(rk8) , public :: dx , dx2 , dx4 , dx8 , dx16 , dxsq
  real(rk8) , public :: c200 , rdxsq , dtsrf , dtabem , dtrad , cpldt
  real(rk8) , public :: xkhmax , xkhz

  integer(ik4) , public :: iboudy , ichem , ipgf , ipptls
  ! usefull flags for chemistry
  integer(ik4) , public :: iaerosol , igaschem , ioxclim , iisoropia
  character(len=6) , pointer , dimension(:) , public :: chtrname
  integer(ik4) , public :: nbin = 1
!
  logical , public :: do_parallel_netcdf_in , do_parallel_netcdf_out
  logical , public :: ifrest , rfstrt , doing_restart , lsync

  integer(ik4) , public :: kchi , kclo , kcmd , cpldbglevel
!
  real(rk8) , public :: akht1 , akht2

  real(rk8) , pointer , dimension(:) , public :: dtau
  real(rk8) , pointer , dimension(:) , public :: hsigma , dsigma , qcon
  real(rk8) , pointer , dimension(:) , public :: sigma
  real(rk8) , pointer , dimension(:) , public :: anudgh , anudgf
  real(rk8) , pointer , dimension(:,:) , public :: twt

  real(rk8) , public :: clfrcv ! Cloud fractional cover for convective precip

  ! Moisture from previous run

  logical , public :: replacemoist = .false.

  ! Non hydrostatic core parameters

  real(rk8) , public :: logp_lrate

  ! Grell cumulus scheme parameters

  real(rk8) , public :: mincld
  real(rk8) , public :: skbmax
  real(rk8) , public :: shrmax_ocn
  real(rk8) , public :: shrmin_ocn
  real(rk8) , public :: edtmax_ocn
  real(rk8) , public :: edtmin_ocn
  real(rk8) , public :: edtmaxo_ocn
  real(rk8) , public :: edtmino_ocn
  real(rk8) , public :: edtmaxx_ocn
  real(rk8) , public :: edtminx_ocn
  real(rk8) , public :: shrmax
  real(rk8) , public :: shrmin
  real(rk8) , public :: edtmax
  real(rk8) , public :: edtmin
  real(rk8) , public :: edtmaxo
  real(rk8) , public :: edtmino
  real(rk8) , public :: edtmaxx
  real(rk8) , public :: edtminx
  real(rk8) , public :: dtauc
  real(rk8) , public :: pbcmax
  real(rk8) , public :: htmax
  real(rk8) , public :: htmin

  ! Emanuel MIT cumulus scheme parameters

  real(rk8) , public :: alphae
  real(rk8) , public :: betae
  real(rk8) , public :: coeffr
  real(rk8) , public :: coeffs
  real(rk8) , public :: cu
  real(rk8) , public :: damp
  real(rk8) , public :: dtmax
  real(rk8) , public :: elcrit_ocn
  real(rk8) , public :: elcrit_lnd
  real(rk8) , public :: entp
  real(rk8) , public :: minsig
  real(rk8) , public :: omtrain
  real(rk8) , public :: omtsnow
  real(rk8) , public :: sigd
  real(rk8) , public :: sigs
  real(rk8) , public :: tlcrit
  real(rk8) , public :: epmax_ocn
  real(rk8) , public :: epmax_lnd
  integer(ik4) , public :: minorig

  ! Tiedtke cumulus scheme parameters

  integer(ik4) , public :: iconv

  real(rk8) , public :: entrdd  ! entrainment rate for cumulus downdrafts

  ! ICONV 1, 2, 3
  real(rk8) , public :: entrpen ! entrainment rate for penetrative convection
  real(rk8) , public :: entrscv ! entrainment rate for shallow convection
  real(rk8) , public :: entrmid ! entrainment rate for midlevel convection
  real(rk8) , public :: cprcon  ! coefficients for determining conversion
                                ! from cloud water to rain
  real(rk8) , public :: entrmax ! Max entrainment

  ! ICONV 4 paramters
  real(rk8) , public :: rcuc_lnd  ! Convective cloud cover for rain evporation
  real(rk8) , public :: rcuc_ocn  ! Convective cloud cover for rain evporation
  real(rk8) , public :: rcpec_lnd ! Coefficient for rain evaporation below cloud
  real(rk8) , public :: rcpec_ocn ! Coefficient for rain evaporation below cloud
  real(rk8) , public :: rhebc_lnd ! Critical relative humidity below
                                  ! cloud at which evaporation starts for land
  real(rk8) , public :: rhebc_ocn ! Critical relative humidity below
                                  ! cloud at which evaporation starts for ocean
  real(rk8) , public :: rprc_lnd  ! coeff for conversion from cloud water
  real(rk8) , public :: rprc_ocn  ! coeff for conversion from cloud water
  real(rk8) , public :: detrpen   ! Detrainment rate for penetrative convection
  real(rk8) , public :: entshalp  ! shallow entrainment factor for entrorg

  ! Kain-Fritsch parameter

  real(rk8) , public :: kf_entrate
  integer(ik4) , public :: kf_trigger
  !

  ! Tweak Global data

  integer(ik4) , public :: itweak
  integer(ik4) , public :: itweak_sst
  integer(ik4) , public :: itweak_temperature
  integer(ik4) , public :: itweak_solar_irradiance
  integer(ik4) , public :: itweak_greenhouse_gases
  real(rk8) , public :: sst_tweak
  real(rk8) , public :: temperature_tweak
  real(rk8) , public :: solar_tweak
  real(rk8) , dimension(5) , public :: gas_tweak_factors

  ! RRTM scheme parameters

  integer(ik4) , public :: irrtm
  integer(ik4) , public :: irrtm_cldov
  integer(ik4) , public :: irrtm_sw_opcliq
  integer(ik4) , public :: irrtm_sw_opcice
  integer(ik4) , public :: inflgsw
  integer(ik4) , public :: iceflgsw
  integer(ik4) , public :: liqflgsw
  integer(ik4) , public :: icld
  integer(ik4) , public :: idrv
  integer(ik4) , public :: irng
  integer(ik4) , public :: inflglw
  integer(ik4) , public :: iceflglw
  integer(ik4) , public :: liqflglw

  ! Radiation schemes common parametrs

  integer(ik4) , public :: iclimao3

  ! UW PBL parameters

  integer(ik4) , public :: iuwvadv
  real(rk8) , public :: rstbl
  real(rk8) , public :: atwo
  real(rk8) , public  :: czero

  ! Holtslag PBL parameters

  real(rk8) , public :: ricr_ocn
  real(rk8) , public :: ricr_lnd
  real(rk8) , public :: zhnew_fac
  integer(ik4) , public :: ifaholtth10
  integer(ik4) , public :: ifaholtmax
  integer(ik4) , public :: ifaholtmin

  ! Chemistry nameliste option

  character(len=8) , public :: chemsimtype
  integer(ik4) , public :: ichcumtra
  integer(ik4) , public :: ichdrdepo
  integer(ik4) , public :: ichremcvc
  integer(ik4) , public :: ichremlsc
  integer(ik4) , public :: ichsursrc
  integer(ik4) , public :: ichsolver
  integer(ik4) , public :: ichdustemd
  integer(ik4) , public :: ichdiag
  integer(ik4) , public :: ichebdy
  integer(ik4) , public :: ichjphcld
  integer(ik4) , public :: ichbion
  real(rk8) , public :: rdstemfac

  ! chemistry species indices that are used not only in chemlib but also in
  ! other interface ( e.g CLM4.5)/ other species are delcared in
  ! chemlib/mod_che_indices
  integer(ik4) , public :: ibchb,ibchl,iochl,iochb,ianh4,iano3

  ! Large scale SUBEX parameters

  integer(ik4) , public :: ncld ! # of bottom model levels with no clouds
  logical , public :: lsrfhack ! Surface radiation hack
  real(rk8) , public :: qck1land
  real(rk8) , public :: qck1oce
  real(rk8) , public :: gulland
  real(rk8) , public :: guloce
  real(rk8) , public :: rhmax
  real(rk8) , public :: rh0land
  real(rk8) , public :: rh0oce
  real(rk8) , public :: caccrlnd
  real(rk8) , public :: cevaplnd
  real(rk8) , public :: caccroce
  real(rk8) , public :: cevapoce
  real(rk8) , public :: tc0
  real(rk8) , public :: cllwcv   ! Cloud liquid water content
                                 ! for convective precip.
  real(rk8) , public :: clfrcvmax
  real(rk8) , public :: cftotmax ! Maximum total cloud fraction for radiation
  real(rk8) , public :: conf     ! Condensation threshold.

  !
  ! SLAB ocean parameters
  !
  integer(ik4) , public :: islab_ocean
  logical , public :: ifslaboc = .false.
  logical , public :: do_qflux_adj
  logical , public :: do_restore_sst
  ! TK mod; restoring time scale for SST in days
  real(rk8) , public :: sst_restore_timescale
  real(rk8) , public :: mixed_layer_depth
  integer(ik4) , dimension(mpy) , public :: stepcount

#ifdef CLM
  ! CLM options
  integer(ik4) , public :: imask
  integer(ik4) , public :: ilawrence_albedo
  real(rk8) , public :: clmfrq
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
  integer , public :: kautoconv
  ! Implict option
  real(rk8) , public :: ksemi
  ! Fall speed values
  real(rk8) , public :: vqxr
  real(rk8) , public :: vqxi
  real(rk8) , public :: vqxs
  ! autoconversion values
  real(rk8) , public :: zauto_rate_khair
  real(rk8) , public :: zauto_rate_kessl
  real(rk8) , public :: zauto_rate_klepi
  real(rk8) , public :: rkconv
  ! limit values
  real(rk8) , public :: rcovpmin
  ! evaporation values
  real(rk8) , public :: rpecons
  !
  ! Debug limit values for tendencies
  !
  real(rk8) , public :: temp_tend_maxval
  real(rk8) , public :: wind_tend_maxval

  ! option for writing tendency diagnostic
  integer(ik4) , public :: idiag

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
    real(rk8) , intent(in) :: a
    iswater = .false.
    if (a > 13.5D0 .and. a < 15.5D0) iswater = .true.
  end function

  logical function isocean(a)
    real(rk8) , intent(in) :: a
    isocean = .false.
    if (a > 14.5D0 .and. a < 15.5D0) isocean = .true.
  end function

  logical function islake(a)
    real(rk8) , intent(in) :: a
    islake = .false.
    if (a > 13.5D0 .and. a < 14.5D0) islake = .true.
  end function

end module mod_runparams
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

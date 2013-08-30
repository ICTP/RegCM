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

  character(len=256) :: namelistfile , prgname

  integer(ik4) , public :: nqx , iqfrst , iqlst
  integer(ik4) , public , parameter :: iqv = 1
  integer(ik4) , public , parameter :: iqc = 2
  integer(ik4) , public , parameter :: iqr = 3
  integer(ik4) , public , parameter :: iqi = 4
  integer(ik4) , public , parameter :: iqs = 5

  type(rcm_time_and_date) , save :: idate0 , idate1 , idate2

  type(rcm_time_and_date) , save :: idatex
  integer(ik4) :: xyear , xmonth , xday , xhour

  type(rcm_time_and_date) , save :: bdydate1 , bdydate2
  type(rcm_time_and_date) , save :: somdate1 , somdate2

  type(rcm_time_interval) , save :: intmdl
  type(rcm_time_interval) , save :: intbdy
  type(rcm_time_interval) , save :: intsom

  real(rk8) :: declin , deltmx
  real(rk8) :: xbctime , xslabtime
  real(rk8) :: calday , twodt

  real(rk8) :: solcon , scon

  ! Step counter. Is zero at idate0, always increasing, never reset.
  integer(ik8) :: ktau
  ! Final number of step for THIS run
  integer(ik8) :: mtau
  ! How many steps for an hour (updates date fields Y m d H)
  integer(ik8) :: khour
  ! How many steps for a day (updates date fields Y m d)
  integer(ik8) :: kday
  ! Output k values for I/O operations.
  integer(ik8) :: katm , krad , kche , ksav , kdbg , kbdy , &
                  ksrf , ksub , klak , krep
  ! Seconds counter in between boundary conditions read
  integer(ik8) :: nbdytime
  ! Step counters to activate surface and radiation schemes
  integer(ik8) :: ntsrf , ntrad
  real(rk8) :: rtsrf , rtrad , rnsrf_for_srffrq , rnsrf_for_day , &
               rnsrf_for_lakfrq , rnsrf_for_subfrq , rnrad_for_chem , &
               rnrad_for_radfrq
  real(rk8) :: afdout , cfdout
  ! Step of surface scheme in one atmosphere I/O interval
  real(rk8) :: rsrf_in_atm
  ! One over seconds in one surface I/O interval
  real(rk8) :: rsrffrq_sec
  ! Model timestep in seconds (real and integer)
  integer(ik8) :: ntsec
  real(rk8) :: dtsec
  ! Internal count for how many SRF outputs per day
  integer(ik8) :: ksts , kstsoff
  !
  ! Cumulus scheme index
  integer(ik4) :: icup
  ! Closure index for Grell
  integer(ik4) :: igcc
  ! Cumulus cloud model
  integer(ik4) :: icumcloud
  ! Boundary layer index
  integer(ik4) :: ibltyp
  ! Lake model activation index
  integer(ik4) :: lakemod
  ! Diurnal cycle SST index
  integer(ik4) :: idcsst
  ! Sea Ice scheme index
  integer(ik4) :: iseaice
  ! Seasonal albedo for desert index
  integer(ik4) :: idesseas
  ! Ocean model switch indexes
  integer(ik4) :: iocnrough , iocnflx , iocncpl
  ! Radiation switch controls
  integer(ik4) :: idirect , iindirect , iemiss , isolconst
  ! Semi-Langrangian advection for tracers
  integer(ik4) :: isladvec
!
  character(len=8) :: scenario
!
  real(rk8) :: dt , dt2 , dtsq , dtcb , dtbdys , rdt
  real(rk8) :: dx , dx2 , dx4 , dx8 , dx16 , dxsq
  real(rk8) :: c200 , rdxsq , dtsrf , dtabem , dtrad , cpldt
  real(rk8) :: xkhmax , xkhz

  integer(ik4) :: iboudy , ichem , ipgf , ipptls
  ! usefull flags for chemistry
  integer(ik4) :: iaerosol , igaschem , ioxclim
  character(len=6) , pointer , dimension(:) :: chtrname
!
  logical :: do_parallel_netcdf_in , do_parallel_netcdf_out
  logical :: ifrest , rfstrt , doing_restart , lsync

  integer(ik4) :: kchi , kclo , kcmd , cpldbglevel
!
  real(rk8) :: akht1 , akht2

  real(rk8) , pointer , dimension(:) :: dtau
  real(rk8) , pointer , dimension(:) :: hsigma , dsigma , qcon
  real(rk8) , pointer , dimension(:) :: sigma
  real(rk8) , pointer , dimension(:,:) :: twt
  real(rk8) , pointer , dimension(:) :: anudg

  real(rk8) :: clfrcv ! Cloud fractional cover for convective precip
  real(rk8) :: cllwcv ! Cloud liquid water content for convective precip.

  ! Grell cumulus scheme parameters

  real(rk8) :: gulland
  real(rk8) :: guloce
  real(rk8) :: mincld
  real(rk8) :: qck1land
  real(rk8) :: qck1oce
  real(rk8) :: rh0land
  real(rk8) :: rh0oce
  real(rk8) :: skbmax
  real(rk8) :: clfrcvmax
  real(rk8) :: shrmax_ocn
  real(rk8) :: shrmin_ocn
  real(rk8) :: edtmax_ocn
  real(rk8) :: edtmin_ocn
  real(rk8) :: edtmaxo_ocn
  real(rk8) :: edtmino_ocn
  real(rk8) :: edtmaxx_ocn
  real(rk8) :: edtminx_ocn
  real(rk8) :: shrmax
  real(rk8) :: shrmin
  real(rk8) :: edtmax
  real(rk8) :: edtmin
  real(rk8) :: edtmaxo
  real(rk8) :: edtmino
  real(rk8) :: edtmaxx
  real(rk8) :: edtminx
  real(rk8) :: dtauc
  real(rk8) :: pbcmax
  real(rk8) :: htmax
  real(rk8) :: htmin

  ! Emanuel MIT cumulus scheme parameters

  real(rk8) :: alphae
  real(rk8) :: betae
  real(rk8) :: coeffr
  real(rk8) :: coeffs
  real(rk8) :: cu
  real(rk8) :: damp
  real(rk8) :: dtmax
  real(rk8) :: elcrit_ocn
  real(rk8) :: elcrit_lnd
  real(rk8) :: entp
  real(rk8) :: minsig
  real(rk8) :: omtrain
  real(rk8) :: omtsnow
  real(rk8) :: sigd
  real(rk8) :: sigs
  real(rk8) :: tlcrit
  real(rk8) :: epmax_ocn
  real(rk8) :: epmax_lnd
  integer(ik4) :: minorig
 
  ! Tiedtke cumulus scheme parameters

  real(rk8) :: entrpen      ! entrainment rate for penetrative convection
  real(rk8) :: entrscv      ! entrainment rate for shallow convection
  real(rk8) :: entrmid      ! entrainment rate for midlevel convection
  real(rk8) :: entrdd       ! entrainment rate for cumulus downdrafts
  real(rk8) :: cmfctop      ! relat. cloud massflux at level above nonbuoyanc
  real(rk8) :: cmtcape      ! CAPE adjustment timescale parameter
  real(rk8) :: zdlev        ! Restrict rainfall up to this elevation
  real(rk8) :: cmfcmax      ! maximum massflux value allowed for
  real(rk8) :: cmfcmin      ! minimum massflux value (for safety)
  real(rk8) :: cmfdeps      ! fractional massflux for downdrafts at lfs
  real(rk8) :: rhcdd        ! relative saturation in downdrafts
  real(rk8) :: cprcon       ! coefficients for determining conversion
                            ! from cloud water to rain
  real(rk8) :: ctrigger     ! coefficients for triggering convection
  real(rk8) :: cmcptop      ! Top pressure for midlevel convection

  integer(ik4) :: iconv

  logical :: lmfpen    !  true if penetrative convection is switched on
  logical :: lmfscv    !  true if shallow convection is switched on
  logical :: lmfmid    !  true if midlevel convection is switched on
  logical :: lmfdd     !  true if cumulus downdraft is switched on
  logical :: lmfdudv   !  true if cumulus friction is switched on

  ! RRTM scheme parameters

  integer(ik4) :: irrtm
  integer(ik4) :: irrtm_cldov
  integer(ik4) :: irrtm_sw_opcliq
  integer(ik4) :: irrtm_sw_opcice
  integer(ik4) :: inflgsw
  integer(ik4) :: iceflgsw
  integer(ik4) :: liqflgsw
  integer(ik4) :: icld
  integer(ik4) :: idrv
  integer(ik4) :: irng
  integer(ik4) :: inflglw
  integer(ik4) :: iceflglw
  integer(ik4) :: liqflglw

  ! Radiation schemes common parametrs

  real(rk8) , public :: cftotmax ! Maximum total cloud fraction for radiation
  integer(ik4) :: iclimao3
  integer(ik4) :: ncld ! # of bottom model levels with no clouds

  ! UW PBL parameters

  integer(ik4) :: iuwvadv
  real(rk8) :: rstbl
  real(rk8) :: atwo

  ! Holtslag PBL parameters

  real(rk8) :: ricr_ocn
  real(rk8) :: ricr_lnd
  real(rk8) :: zhnew_fac
  integer(ik4) :: ifaholtth10
  integer(ik4) :: ifaholtmax
  integer(ik4) :: ifaholtmin

  ! Chemistry nameliste option

  character(len=8) :: chemsimtype 
  integer(ik4) :: ichcumtra
  integer(ik4) :: ichdrdepo
  integer(ik4) :: ichremcvc
  integer(ik4) :: ichremlsc
  integer(ik4) :: ichsursrc
  integer(ik4) :: ichsolver
  integer(ik4) :: ichdustemd
  integer(ik4) :: ichdiag
  integer(ik4) :: ichebdy
  real(rk8) :: rdstemfac

  ! Large scale SUBEX parameters

  real(rk8) :: caccr
  real(rk8) :: cevap
  real(rk8) :: rhmax
  real(rk8) :: tc0
  real(rk8) :: conf
  integer(ik4) :: iconvlwp

  !
  ! SLAB ocean parameters
  !
  integer(ik4) :: islab_ocean
  logical :: ifslaboc = .false.
  logical :: do_qflux_adj
  logical :: do_restore_sst
  ! TK mod; restoring time scale for SST in days
  real(rk8) :: sst_restore_timescale 
  real(rk8) :: mixed_layer_depth 
  integer(ik4) , dimension(mpy) :: stepcount

  ! CLM options
#ifdef CLM
  integer(ik4) :: imask
  integer(ik4) :: ilawrence_albedo
  real(rk8) :: clmfrq
  character(len=256) :: dirclm
#endif

  !
  ! New microphys parameters
  !
  ! Total water and enthalpy budget on/off
  logical :: budget_compute
  ! Super saturation option
  integer :: nssopt
  ! Choose the autoconversion paramaterization
  integer :: kautoconv
  ! Implict option
  real(rk8) :: ksemi

  !
  ! Debug limit values for tendencies
  !
  real(rk8) :: temp_tend_maxval
  real(rk8) :: wind_tend_maxval

  ! option for writing tendency diagnostic
  integer(ik4) :: idiag

  data doing_restart /.false./

  contains

  subroutine allocate_mod_runparams
    implicit none
    call getmem1d(hsigma,1,kz,'mod_runparams:hsigma')
    call getmem1d(dsigma,1,kz,'mod_runparams:dsigma')
    call getmem1d(qcon,1,kz,'mod_runparams:qcon')
    call getmem1d(sigma,1,kzp1,'mod_runparams:sigma')
    call getmem1d(anudg,1,kz,'mod_runparams:anudg')
    call getmem2d(twt,1,kz,1,2,'mod_runparams:twt')
    call getmem1d(dtau,1,nsplit,'mod_runparams:nsplit')
  end subroutine allocate_mod_runparams
!
  logical function iswater(a)
    real(rk8) , intent(in) :: a
    iswater = .false.
    if (a > 13.5D0 .and. a < 15.5D0) iswater = .true.
  end function
!
  logical function isocean(a)
    real(rk8) , intent(in) :: a
    isocean = .false.
    if (a > 14.5D0 .and. a < 15.5D0) isocean = .true.
  end function
!
  logical function islake(a)
    real(rk8) , intent(in) :: a
    islake = .false.
    if (a > 13.5D0 .and. a < 14.5D0) islake = .true.
  end function
!
  subroutine cdumlbcs(cdum)
    implicit none
    character(len=*) , intent(out) :: cdum
    select case (iboudy)
      case(0)
       write (cdum,'(a)') 'Fixed'
      case(1)
       write (cdum,'(a)') 'Relaxation, linear technique'
      case(2)
       write (cdum,'(a)') 'Time-dependent'
      case(3)
       write (cdum,'(a)') 'Time and inflow/outflow dependent'
      case(4)
       write (cdum,'(a)') 'Sponge (Perkey & Kreitzberg, MWR 1976)'
      case(5)
       write (cdum,'(a)') 'Relaxation, exponential technique'
      case default 
       write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdumlbcs

  subroutine cdumcums(cdum)
    implicit none
    character(len=*) , intent(out) :: cdum
    select case (icup)
      case(1)
       write (cdum,'(a)') 'Kuo'
      case(2)
       write (cdum,'(a)') 'Grell'
      case(3)
       write (cdum,'(a)') 'Betts-Miller (1986)'
      case(4)
       write (cdum,'(a)') 'Emanuel (1991)'
      case(5)
        write (cdum,'(a)') 'Tiedtke (1986)'
      case(96)
        write (cdum,'(a)') 'Tiedtke (1986) over ocean, Grell over land'
      case(97)
        write (cdum,'(a)') 'Tiedtke (1986) over ocean, Emanuel (1991) over land'
      case(98)
        write (cdum,'(a)') 'Grell over ocean, Emanuel (1991) over land'
      case(99)
        write (cdum,'(a)') 'Emanuel (1991) over ocean, Grell over land'
      case default 
        write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdumcums

  subroutine cdumcumcl(cdum)
    implicit none
    character(len=*) , intent(out) :: cdum
    select case (igcc)
      case(1)
        write (cdum,'(a)') 'Arakawa & Schubert (1974)'
      case(2)
        write (cdum,'(a)') 'Fritsch & Chappell (1980)'
      case default 
        write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdumcumcl

  subroutine cdumpbl(cdum)
    implicit none
    character(len=*) , intent(out) :: cdum
    select case (ibltyp)
      case(0)
        write (cdum,'(a)') 'Frictionless'
      case(1)
        write (cdum,'(a)') 'Holtslag PBL (Holtslag, 1990)'
      case(2)
        write (cdum,'(a)') 'UW PBL (Bretherton and McCaa, 2004)'
      case(99)
        write (cdum,'(a)') 'Holtslag PBL, with UW in diag. mode'
      case default 
        write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdumpbl

  subroutine cdummoist(cdum)
    implicit none
    character(len=*) , intent(out) :: cdum
    select case (ipptls)
      case(1)
        write (cdum,'(a)') 'Explicit moisture (SUBEX; Pal et al 2000)'
      case(2)
        write (cdum,'(a)') 'Explicit moisture Nogherotto/Tompkins'
      case default 
        write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdummoist

  subroutine cdumocnflx(cdum)
    implicit none
    character(len=*) , intent(out) :: cdum
    select case (iocnflx)
      case(1)
        write (cdum,'(a)') 'Use BATS1e Monin-Obukhov'
      case(2)
        write (cdum,'(a)') 'Zeng et al (1998)'
      case default 
        write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdumocnflx

  subroutine cdumpgfs(cdum)
    implicit none
    character(len=*) , intent(out) :: cdum
    select case (ipgf)
      case(0)
        write (cdum,'(a)') 'Use full fields'
      case(1)
        write (cdum,'(a)') 'Hydrostatic deduction with perturbation temperature'
      case default 
        write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdumpgfs

  subroutine cdumemiss(cdum)
    implicit none
    character(len=*) , intent(out) :: cdum
    select case (iemiss)
      case(0)
        write (cdum,'(a)') 'No'
      case(1)
        write (cdum,'(a)') 'Yes'
      case default 
        write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdumemiss

  subroutine cdumlakes(cdum)
    implicit none
    character(len=*) , intent(out) :: cdum
    select case (lakemod)
      case(0)
        write (cdum,'(a)') 'No'
      case(1)
        write (cdum,'(a)') 'Yes'
      case default 
        write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdumlakes

  subroutine cdumchems(cdum)
    implicit none
    character(len=*) , intent(out) :: cdum
    select case (ichem)
      case(0)
        write (cdum,'(a)') 'Not active'
      case(1)
        write (cdum,'(a)') 'Active'
      case default 
        write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdumchems

  subroutine cdumdcsst(cdum)
    implicit none
    character(len=*) , intent(out) :: cdum
    select case (idcsst)
      case(0)
        write (cdum,'(a)') 'Not active'
      case(1)
        write (cdum,'(a)') 'Active'
      case default 
        write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdumdcsst

  subroutine cdumseaice(cdum)
    implicit none
    character(len=*) , intent(out) :: cdum
    select case (iseaice)
      case(0)
        write (cdum,'(a)') 'Not active'
      case(1)
        write (cdum,'(a)') 'Active'
      case default 
        write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdumseaice

  subroutine cdumdesseas(cdum)
    implicit none
    character(len=*) , intent(out) :: cdum
    select case (idesseas)
      case(0)
        write (cdum,'(a)') 'Not active'
      case(1)
        write (cdum,'(a)') 'Active'
      case default 
        write (cdum,'(a)') 'Unknown or not specified'
    end select
  end subroutine cdumdesseas

end module mod_runparams

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

module mod_rad_aerosol

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_runparams
  use mod_constants
  use mod_memutil
  use mod_mpmessage
  use mod_stdatm
  use mod_rad_common
  use mod_regcm_types
  use mod_mppparam
  use mod_nhinterp
  use mod_kdinterp
  use mod_interp
  use mod_zita
  use mod_date
  use mod_stdio
  use mod_vertint
  use parrrsw , only : nbndsw
  use parrrtm , only : nbndlw
  use netcdf
  use mo_simple_plumes
  implicit none

  private

  public :: tauxar3d , tauasc3d , gtota3d , ftota3d
  public :: tauxar , tauasc , gtota , ftota
  public :: aermmr , aertrlw , tauxar3d_lw
  public :: allocate_mod_rad_aerosol , aeroppt
  public :: read_aerclima , close_aerclima
  public :: init_aeroppdata , read_aeroppdata
  public :: cmip6_plume_profile
  public :: aerclima_ntr , aerclima_nbin
  !
  character(len=256) :: macv2sp_hist , macv2sp_scen
  real(rk8) , pointer , dimension(:) :: lambdaw
  real(rk8) , pointer , dimension(:) :: latr4 , lonr4 , altr4
  real(rk8) , pointer , dimension(:,:) :: z , dz
  real(rk8) , pointer , dimension(:) :: dnovrnr4
  real(rk8) , pointer , dimension(:,:,:) :: extprofr4
  real(rk8) , pointer , dimension(:,:,:) :: ssaprofr4
  real(rk8) , pointer , dimension(:,:,:) :: asyprofr4
  real(rkx) , pointer , dimension(:) :: lat , lon
  real(rkx) , pointer , dimension(:,:) :: alon , alat
  real(rkx) , pointer , dimension(:,:,:,:) :: extprof
  real(rkx) , pointer , dimension(:,:,:,:) :: ssaprof
  real(rkx) , pointer , dimension(:,:,:,:) :: asyprof
  real(rkx) , pointer , dimension(:,:,:) :: zpr3d , zdzr3d
  real(rkx) , pointer , dimension(:,:,:) :: rdvar
  real(rkx) , pointer , dimension(:,:,:) :: hzivar
  real(rkx) , pointer , dimension(:,:,:) :: plvar , p1 , p2
  real(rkx) , pointer , dimension(:,:,:,:) :: plext1 , plext2
  real(rkx) , pointer , dimension(:,:,:,:) :: plssa1 , plssa2
  real(rkx) , pointer , dimension(:,:,:,:) :: plasy1 , plasy2
  real(rkx) , pointer , dimension(:,:,:,:) :: pldp1 , pldp2 , pl1 , pl2
  real(rkx) , pointer , dimension(:,:,:) :: sgvar
  real(rkx) , pointer , dimension(:,:,:,:) :: sgext1 , sgext2
  real(rkx) , pointer , dimension(:,:,:,:) :: sgssa1 , sgssa2
  real(rkx) , pointer , dimension(:,:,:,:) :: sgasy1 , sgasy2
  type(h_interpolator) :: hint
  integer(ik4) :: ncid = -1
  integer(ik4) :: clnlon , clnlat , clnlev
  !
  integer(ik4) , parameter :: ncoefs = 5  ! Number of coefficients
  integer(ik4) , parameter :: nwav = 19
  integer(ik4) , parameter :: nih = 8

  logical :: aerclima_init = .false.
  integer(ik4) , parameter :: aerclima_ntr = 12
  integer(ik4) , parameter :: aerclima_nbin = 4

  integer(ik4) , parameter :: nacwb  = 5 ! number of waveband MERRA aerosol clim
                                         ! might be completed in the future
  integer(ik4) :: ncaec , ncstatus , naetime
  integer(ik4) , dimension(aerclima_ntr) :: ncaevar
  type(rcm_time_and_date) , dimension(:) , allocatable :: aetime
  type(rcm_time_and_date) :: d1 , d2
  type(rcm_time_interval) :: aefreq
  real(rkx) :: aerfreq
  real(rkx) , pointer , dimension(:,:,:,:) :: aerm1 , aerm2
  real(rkx) , pointer , dimension(:,:,:) :: aerio

  real(rkx) , parameter :: d10e5  = 1.0e+5_rkx
  real(rkx) , parameter :: d10e4  = 1.0e+4_rkx
  real(rkx) , parameter :: minimum_aerosol = 1.0e-14_rkx
  real(rkx) , parameter :: minimum_utaer   = 1.0e-10_rkx
  real(rkx) , parameter :: minimum_waer   = 1.0e-30_rkx
  real(rkx) , parameter :: minimum_gaer   = 1.0e-20_rkx
  real(rkx) , parameter :: fiveothree  = d_five/d_three
  !
  ! Sulfate param for standard scheme / works only with rad standard
  ! (Brieglieb et al.)
  real(rkx) , pointer , dimension(:,:,:) :: aermmr

  ! optical properties for dust and org / for rrtm and standard scheme

  real(rkx) , dimension(nbndsw) :: gsbc_hb_rrtm , gsbc_hl_rrtm , &
    gsoc_hb_rrtm , gsoc_hl_rrtm , ksbc_hb_rrtm , ksbc_hl_rrtm ,  &
    ksoc_hb_rrtm , ksoc_hl_rrtm , wsbc_hb_rrtm , wsbc_hl_rrtm ,  &
    wsoc_hb_rrtm , wsoc_hl_rrtm

  real(rkx) , dimension(nspi,2) :: gssslt , kssslt , wssslt
  !
  ! Depth
  !
  real(rkx) , pointer , dimension(:,:) :: path
  !
  ! Aerosol optical properties (for the mixing)
  !
  real(rkx) , pointer , dimension(:) :: gsbc_hb , gsbc_hl , gsoc_hb , &
            gsoc_hl , ksbc_hb , ksbc_hl , ksoc_hb , ksoc_hl ,         &
            wsbc_hb , wsbc_hl , wsoc_hb , wsoc_hl , gssm1 , gssm2 ,   &
            kssm1 , kssm2 , wssm1 , wssm2

  real(rkx) , pointer,  dimension(:,:) :: gsdust , ksdust , wsdust , ksdust_lw

  real(rkx) , pointer , dimension(:,:,:) :: ftota3d , gtota3d , &
    tauasc3d , tauxar3d, tauxar3d_lw
  real(rkx) , pointer , dimension(:,:) :: ftota , gtota , tauasc , tauxar
  !
  ! Work arrays for aeroppt (aerosol individual optical properties SW)
  !
  ! uaer, tauxar  - Aerosol radiative properties (local arrays)
  ! wa            - Aerosol single scattering albedo
  ! ga            - Aerosol asimmetry parameter
  ! fa            - Aerosol forward scattered fraction
  ! utaer, tauaer - Total column aerosol extinction
  ! waer          - Aerosol single scattering albedo
  ! gaer          - Aerosol asymmetry parameter
  ! faer          - Aerosol forward scattered fraction
  !
  !
  real(rkx) , pointer , dimension(:,:) :: aermtot , aervtot
  real(rkx) , pointer , dimension(:,:,:) :: fa , ga , tx , uaer , wa
  real(rkx) , pointer , dimension(:,:) :: faer , gaer , tauaer , utaer , waer
  integer(ik4) :: npoints , nj
  integer(ik4) :: nband
  !
  ! Aersol LW optical properties
  !
  real(rkx) , pointer , dimension(:,:,:) ::  aertrlw
  !
  !--------------------------------------------------------------------------
  !                  DATA SECTION
  !--------------------------------------------------------------------------
  !
  ! kscoef  - specific extinction (m2/g)
  ! wscoef  - single partical albedo
  ! gscoef  - asymmetry parameter
  ! ksbase  - specific extinction (m2/g) base
  ! wsbase  - single partical albedo base
  ! gsbase  - asymmetry parameter base
  ! ksdust  - specific extinction (m2/g) dust
  ! wsdust  - single partical albedo dust
  ! gsdust  - asymmetry parameter dust
  !
  !
  ! DUST OP data base for external mixing : maximum of 4 bin for the
  ! momeent , determined from Zender et al.
  !
  ! DATA section for optical properties relative to RRTM
  ! based on of line calculation considering the Kok et al., 2011 distribution
  !

  real(rkx) , dimension(nspi) , parameter :: ksbase = [            &
       5.206e+0_rkx , 5.206e+0_rkx , 5.206e+0_rkx , 5.206e+0_rkx , &
       5.206e+0_rkx , 5.206e+0_rkx , 5.206e+0_rkx , 3.203e+0_rkx , &
       3.203e+0_rkx , 1.302e+0_rkx , 5.992e-1_rkx , 2.948e-1_rkx , &
       1.475e-1_rkx , 7.387e-2_rkx , 1.683e-1_rkx , 2.655e-1_rkx , &
       5.770e-2_rkx , 2.290e-1_rkx , 2.270e-1_rkx ]

  real(rkx) , dimension(nspi) , parameter :: wsbase = [            &
       7.371e-8_rkx , 7.371e-8_rkx , 7.371e-8_rkx , 7.371e-8_rkx , &
       7.371e-8_rkx , 7.371e-8_rkx , 7.371e-8_rkx , 6.583e-8_rkx , &
       6.583e-8_rkx , 3.656e-6_rkx , 4.919e-5_rkx , 3.539e-3_rkx , &
       2.855e-2_rkx , 2.126e-1_rkx , 8.433e-1_rkx , 9.653e-1_rkx , &
       6.198e-1_rkx , 9.642e-1_rkx , 9.699e-1_rkx ]

  real(rkx) , dimension(nspi) , parameter :: gsbase = [            &
       6.899e-1_rkx , 6.899e-1_rkx , 6.899e-1_rkx , 6.899e-1_rkx , &
       6.899e-1_rkx , 6.899e-1_rkx , 6.899e-1_rkx , 6.632e-1_rkx , &
       6.632e-1_rkx , 5.912e-1_rkx , 5.111e-1_rkx , 4.269e-1_rkx , &
       3.321e-1_rkx , 2.197e-1_rkx , 1.305e-1_rkx , 7.356e-2_rkx , &
       1.602e-1_rkx , 6.883e-2_rkx , 6.304e-2_rkx ]

  real(rkx) , dimension(nspi) , parameter :: ksbc_hb_stand = [ &
       20.7830_rkx , 17.2120_rkx , 15.8640_rkx , 15.0530_rkx , &
       14.3040_rkx , 13.6130_rkx , 11.9660_rkx ,  6.5782_rkx , &
        4.3961_rkx ,  4.3800_rkx ,  2.1100_rkx ,  2.1100_rkx , &
        2.1100_rkx ,  2.1100_rkx ,  1.4000_rkx ,  1.4000_rkx , &
        1.4000_rkx ,  1.4000_rkx ,  1.4000_rkx ]

  real(rkx) , dimension(nspi) , parameter :: wsbc_hb_stand = [     &
       0.245240_rkx , 0.209620_rkx , 0.195000_rkx , 0.185860_rkx , &
       0.177190_rkx , 0.168970_rkx , 0.148490_rkx , 0.071748_rkx , &
       0.037536_rkx , 0.089000_rkx , 0.025000_rkx , 0.025000_rkx , &
       0.025000_rkx , 0.025000_rkx , 0.009000_rkx , 0.009000_rkx , &
       0.009000_rkx , 0.009000_rkx , 0.009000_rkx ]

  real(rkx) , dimension(nspi) , parameter :: gsbc_hb_stand = [     &
       0.213870_rkx , 0.171540_rkx , 0.155640_rkx , 0.146100_rkx , &
       0.137320_rkx , 0.129230_rkx , 0.110060_rkx , 0.049175_rkx , &
       0.026638_rkx , 0.220000_rkx , 0.123000_rkx , 0.123000_rkx , &
       0.123000_rkx , 0.123000_rkx , 0.073000_rkx , 0.073000_rkx , &
       0.073000_rkx , 0.073000_rkx , 0.073000_rkx ]

  real(rkx) , dimension(nspi) , parameter :: ksbc_hl_stand = [ &
       14.8510_rkx , 14.2580_rkx , 13.9430_rkx , 13.7240_rkx , &
       13.5070_rkx , 13.2950_rkx , 12.7220_rkx ,  9.4434_rkx , &
        6.9653_rkx ,  4.3800_rkx ,  2.1100_rkx ,  2.1100_rkx , &
        2.1100_rkx ,  2.1100_rkx ,  1.4000_rkx ,  1.4000_rkx , &
        1.4000_rkx ,  1.4000_rkx ,  1.4000_rkx ]

  real(rkx) , dimension(nspi) , parameter :: wsbc_hl_stand = [ &
       0.46081_rkx , 0.44933_rkx , 0.44397_rkx , 0.44065_rkx , &
       0.43737_rkx , 0.43394_rkx , 0.42346_rkx , 0.35913_rkx , &
       0.29579_rkx , 0.08900_rkx , 0.02500_rkx , 0.02500_rkx , &
       0.02500_rkx , 0.02500_rkx , 0.00900_rkx , 0.00900_rkx , &
       0.00900_rkx , 0.00900_rkx , 0.00900_rkx ]

  real(rkx) , dimension(nspi) , parameter :: gsbc_hl_stand = [ &
       0.69038_rkx , 0.65449_rkx , 0.63711_rkx , 0.62542_rkx , &
       0.61407_rkx , 0.60319_rkx , 0.57467_rkx , 0.42050_rkx , &
       0.30660_rkx , 0.22000_rkx , 0.12300_rkx , 0.12300_rkx , &
       0.12300_rkx , 0.12300_rkx , 0.07300_rkx , 0.07300_rkx , &
       0.07300_rkx , 0.07300_rkx , 0.07300_rkx ]

  real(rkx) , dimension(nspi) , parameter :: ksoc_hb_stand = [         &
       6.0584e+0_rkx , 6.0654e+0_rkx , 6.1179e+0_rkx , 6.0102e+0_rkx , &
       5.8000e+0_rkx , 5.6957e+0_rkx , 5.6494e+0_rkx , 4.3283e+0_rkx , &
       3.1485e+0_rkx , 1.3020e+0_rkx , 5.9920e-1_rkx , 2.9480e-1_rkx , &
       1.4750e-1_rkx , 7.3870e-2_rkx , 1.6830e-1_rkx , 2.6550e-1_rkx , &
       5.7700e-2_rkx , 2.2900e-1_rkx , 2.2700e-1_rkx ]

  real(rkx) , dimension(nspi) , parameter :: wsoc_hb_stand = [ &
       0.91735_rkx , 0.92365_rkx , 0.92941_rkx , 0.93067_rkx , &
       0.93311_rkx , 0.93766_rkx , 0.94042_rkx , 0.95343_rkx , &
       0.95480_rkx , 0.70566_rkx , 0.70566_rkx , 0.70566_rkx , &
       0.70566_rkx , 0.70566_rkx , 0.70566_rkx , 0.70566_rkx , &
       0.70566_rkx , 0.70566_rkx , 0.70566_rkx ]

  real(rkx) , dimension(nspi) , parameter :: gsoc_hb_stand = [             &
       0.67489e+0_rkx , 0.67003e+0_rkx , 0.67725e+0_rkx , 0.65487e+0_rkx , &
       0.65117e+0_rkx , 0.66116e+0_rkx , 0.64547e+0_rkx , 0.60033e+0_rkx , &
       0.55389e+0_rkx , 5.91200e-1_rkx , 5.11100e-1_rkx , 4.26900e-1_rkx , &
       3.32100e-1_rkx , 2.19700e-1_rkx , 1.30500e-1_rkx , 7.35600e-2_rkx , &
       1.60200e-1_rkx , 6.88300e-2_rkx , 6.30400e-2_rkx ]

  real(rkx) , dimension(nspi) , parameter :: ksoc_hl_stand = [         &
       3.5430e+0_rkx , 3.6230e+0_rkx , 3.7155e+0_rkx , 3.7120e+0_rkx , &
       3.6451e+0_rkx , 3.6376e+0_rkx , 3.6844e+0_rkx , 3.4588e+0_rkx , &
       2.9846e+0_rkx , 1.3020e+0_rkx , 5.9920e-1_rkx , 2.9480e-1_rkx , &
       1.4750e-1_rkx , 7.3870e-2_rkx , 1.6830e-1_rkx , 2.6550e-1_rkx , &
       5.7700e-2_rkx , 2.2900e-1_rkx , 2.2700e-1_rkx ]

  real(rkx) , dimension(nspi) , parameter :: wsoc_hl_stand = [ &
       0.87931_rkx , 0.88292_rkx , 0.89214_rkx , 0.89631_rkx , &
       0.89996_rkx , 0.90540_rkx , 0.90805_rkx , 0.93423_rkx , &
       0.95012_rkx , 0.85546_rkx , 0.85546_rkx , 0.85546_rkx , &
       0.85546_rkx , 0.85546_rkx , 0.85546_rkx , 0.85546_rkx , &
       0.85546_rkx , 0.85546_rkx , 0.85546_rkx ]

  real(rkx) , dimension(nspi) , parameter :: gsoc_hl_stand = [             &
       0.73126e+0_rkx , 0.71089e+0_rkx , 0.72042e+0_rkx , 0.69924e+0_rkx , &
       0.69908e+0_rkx , 0.70696e+0_rkx , 0.68479e+0_rkx , 0.64879e+0_rkx , &
       0.63433e+0_rkx , 5.91200e-1_rkx , 5.11100e-1_rkx , 4.26900e-1_rkx , &
       3.32100e-1_rkx , 2.19700e-1_rkx , 1.30500e-1_rkx , 7.35600e-2_rkx , &
       1.60200e-1_rkx , 6.88300e-2_rkx , 6.30400e-2_rkx ]

  real(rkx) , dimension(8) , parameter :: rhp = [   &
        0.00_rkx , 0.50_rkx , 0.70_rkx , 0.80_rkx , &
        0.90_rkx , 0.95_rkx , 0.98_rkx , 0.99_rkx ]

  ! sea salt oppt  param for standard scheme

  real(rkx) , dimension(nbndsw) , parameter :: kssm1_rrtm = [       &
      0.0600_rkx,  0.0840_rkx,  0.1110_rkx, 0.1420_rkx, 0.1940_rkx, &
      0.3010_rkx,  0.4110_rkx,  0.9340_rkx, 2.3400_rkx, 4.8800_rkx, &
      8.5250_rkx, 12.8970_rkx, 15.9180_rkx, 0.0281_rkx ]

  real(rkx) , dimension(nbndsw) , parameter :: gssm1_rrtm = [     &
      0.0180_rkx, 0.0270_rkx, 0.0390_rkx, 0.0510_rkx, 0.0680_rkx, &
      0.0990_rkx, 0.1260_rkx, 0.2160_rkx, 0.3830_rkx, 0.5310_rkx, &
      0.6290_rkx, 0.6930_rkx, 0.7270_rkx, 0.0046_rkx ]

  real(rkx) , dimension(nbndsw) , parameter :: wssm1_rrtm = [     &
      0.0710_rkx, 0.1230_rkx, 0.1880_rkx, 0.2550_rkx, 0.3400_rkx, &
      0.4610_rkx, 0.5410_rkx, 0.6710_rkx, 0.7970_rkx, 0.8480_rkx, &
      0.8750_rkx, 0.8840_rkx, 0.8750_rkx, 0.0114_rkx ]

  real(rkx) , dimension(nbndsw) , parameter :: kssm2_rrtm = [       &
      0.0450_rkx,  0.0650_rkx,  0.0930_rkx, 0.1290_rkx, 0.1920_rkx, &
      0.3340_rkx,  0.4810_rkx,  1.1640_rkx, 2.9190_rkx, 5.7520_rkx, &
      9.2420_rkx, 15.5280_rkx, 13.7550_rkx, 0.0181_rkx ]

  real(rkx) , dimension(nbndsw) , parameter :: gssm2_rrtm = [     &
      0.0258_rkx, 0.0390_rkx, 0.0550_rkx, 0.0720_rkx, 0.0970_rkx, &
      0.1410_rkx, 0.1790_rkx, 0.2950_rkx, 0.4790_rkx, 0.5990_rkx, &
      0.6710_rkx, 0.7130_rkx, 0.7210_rkx, 0.0070_rkx ]

  real(rkx) , dimension(nbndsw) , parameter :: wssm2_rrtm = [     &
      0.1670_rkx, 0.2680_rkx, 0.3750_rkx, 0.4680_rkx, 0.5650_rkx, &
      0.6760_rkx, 0.7380_rkx, 0.8180_rkx, 0.8870_rkx, 0.9120_rkx, &
      0.9230_rkx, 0.9210_rkx, 0.9040_rkx, 0.0290_rkx ]

  character(len=6) , dimension(aerclima_ntr) , parameter :: aerclima_chtr = &
     [ 'BC_HB ' , 'BC_HL ' , 'OC_HB ' , 'OC_HL ' , 'SO2   ' , 'SO4   ' ,    &
       'SSLT01' , 'SSLT02' , 'DUST01' , 'DUST02' , 'DUST03' , 'DUST04' ]

  real(rkx) , dimension(nspi,ncoefs) , parameter :: kscoef = reshape(   &
  [ 0.1126E+02_rkx , 0.1126E+02_rkx , 0.1126E+02_rkx , 0.1126E+02_rkx , &
    0.1126E+02_rkx , 0.1126E+02_rkx , 0.1126E+02_rkx , 0.1124E+02_rkx , &
    0.1124E+02_rkx , 0.1222E+02_rkx , 0.1357E+02_rkx , 0.1557E+02_rkx , &
    0.1758E+02_rkx , 0.1597E+02_rkx , 0.2107E+02_rkx , -.2424E+00_rkx , &
    0.2535E+02_rkx , -.1545E+00_rkx , 0.8835E+00_rkx , -.2502E+00_rkx , &
    -.2502E+00_rkx , -.2502E+00_rkx , -.2502E+00_rkx , -.2502E+00_rkx , &
    -.2502E+00_rkx , -.2502E+00_rkx , -.3040E+00_rkx , -.3040E+00_rkx , &
    -.3770E+00_rkx , -.4190E+00_rkx , -.4353E+00_rkx , -.4389E+00_rkx , &
    -.4337E+00_rkx , -.3041E+00_rkx , -.1770E+00_rkx , -.2270E+00_rkx , &
    -.1661E+00_rkx , -.1590E+00_rkx , -.1087E+01_rkx , -.1087E+01_rkx , &
    -.1087E+01_rkx , -.1087E+01_rkx , -.1087E+01_rkx , -.1087E+01_rkx , &
    -.1087E+01_rkx , -.1088E+01_rkx , -.1088E+01_rkx , -.1089E+01_rkx , &
    -.1087E+01_rkx , -.1083E+01_rkx , -.1078E+01_rkx , -.1073E+01_rkx , &
    -.1067E+01_rkx , -.1032E+01_rkx , -.1052E+01_rkx , -.1030E+01_rkx , &
    -.1029E+01_rkx , -.1794E+03_rkx , -.1794E+03_rkx , -.1794E+03_rkx , &
    -.1794E+03_rkx , -.1794E+03_rkx , -.1794E+03_rkx , -.1794E+03_rkx , &
    -.1776E+03_rkx , -.1776E+03_rkx , -.1898E+03_rkx , -.2070E+03_rkx , &
    -.2382E+03_rkx , -.2716E+03_rkx , -.2510E+03_rkx , -.2494E+03_rkx , &
    0.1469E+00_rkx , -.2528E+03_rkx , -.4698E-03_rkx , -.2838E+02_rkx , &
    0.1556E+02_rkx , 0.1556E+02_rkx , 0.1556E+02_rkx , 0.1556E+02_rkx , &
    0.1556E+02_rkx , 0.1556E+02_rkx , 0.1556E+02_rkx , 0.1537E+02_rkx , &
    0.1537E+02_rkx , 0.1504E+02_rkx , 0.1478E+02_rkx , 0.1486E+02_rkx , &
    0.1505E+02_rkx , 0.1527E+02_rkx , 0.1166E+02_rkx , 0.1947E+01_rkx , &
    0.9888E+01_rkx , 0.7275E-01_rkx , 0.2734E+02_rkx ], [nspi,ncoefs])

  real(rkx) , dimension(nspi,ncoefs) , parameter :: wscoef = reshape(   &
  [ 0.2492E+01_rkx , 0.2492E+01_rkx , 0.2492E+01_rkx , 0.2492E+01_rkx , &
    0.2492E+01_rkx , 0.2492E+01_rkx , 0.2492E+01_rkx , 0.1139E+01_rkx , &
    0.1139E+01_rkx , 0.1848E+01_rkx , 0.5459E+01_rkx , 0.1187E+01_rkx , &
    -.3640E+01_rkx , -.5634E+01_rkx , 0.1826E+00_rkx , 0.2164E+01_rkx , &
    0.2268E+00_rkx , 0.2178E+01_rkx , 0.1713E+01_rkx , -.5210E-01_rkx , &
    -.5210E-01_rkx , -.5210E-01_rkx , -.5210E-01_rkx , -.5210E-01_rkx , &
    -.5210E-01_rkx , -.5210E-01_rkx , -.1110E-01_rkx , -.1110E-01_rkx , &
    -.3920E-03_rkx , 0.9357E+00_rkx , 0.2241E+00_rkx , 0.2552E+00_rkx , &
    0.2068E+00_rkx , 0.6588E-01_rkx , 0.1194E+00_rkx , 0.3266E-01_rkx , &
    0.1151E+00_rkx , 0.9166E-01_rkx , -.1036E+01_rkx , -.1036E+01_rkx , &
    -.1036E+01_rkx , -.1036E+01_rkx , -.1036E+01_rkx , -.1036E+01_rkx , &
    -.1036E+01_rkx , -.1011E+01_rkx , -.1011E+01_rkx , -.9924E+00_rkx , &
    -.1626E+01_rkx , -.1226E+01_rkx , -.1168E+01_rkx , -.1122E+01_rkx , &
    -.1098E+01_rkx , -.1044E+01_rkx , -.1064E+01_rkx , -.1042E+01_rkx , &
    -.1039E+01_rkx , -.4398E+02_rkx , -.4398E+02_rkx , -.4398E+02_rkx , &
    -.4398E+02_rkx , -.4398E+02_rkx , -.4398E+02_rkx , -.4398E+02_rkx , &
    -.7754E+01_rkx , -.7754E+01_rkx , -.1607E+01_rkx , -.5282E+01_rkx , &
    0.1442E+02_rkx , 0.4458E+02_rkx , 0.7528E+02_rkx , -.1996E-01_rkx , &
    -.3221E+02_rkx , -.2677E-01_rkx , -.3325E+02_rkx , -.2660E+02_rkx , &
    0.1724E+02_rkx , 0.1724E+02_rkx , 0.1724E+02_rkx , 0.1724E+02_rkx , &
    0.1724E+02_rkx , 0.1724E+02_rkx , 0.1724E+02_rkx , 0.6737E+01_rkx , &
    0.6737E+01_rkx , 0.8587E+00_rkx , 0.1066E+01_rkx , -.1402E+02_rkx , &
    0.1152E+02_rkx , 0.1290E+02_rkx , 0.1618E+00_rkx , 0.1564E+02_rkx , &
    0.1309E+00_rkx , 0.1600E+02_rkx , 0.1629E+02_rkx ], [nspi,ncoefs])

  real(rkx) , dimension(nspi,ncoefs) , parameter :: gscoef = reshape(   &
  [ -.9874E+00_rkx , -.9874E+00_rkx , -.9874E+00_rkx , -.9874E+00_rkx , &
    -.9874E+00_rkx , -.9874E+00_rkx , -.9874E+00_rkx , -.3666E+00_rkx , &
    -.3666E+00_rkx , 0.5824E+00_rkx , 0.1238E+01_rkx , 0.2299E+01_rkx , &
    0.3037E+01_rkx , 0.4683E+01_rkx , 0.3842E+01_rkx , 0.3237E+01_rkx , &
    0.4181E+01_rkx , 0.3378E+01_rkx , 0.3943E+01_rkx , -.3033E+02_rkx , &
    -.3033E+02_rkx , -.3033E+02_rkx , -.3033E+02_rkx , -.3033E+02_rkx , &
    -.3033E+02_rkx , -.3033E+02_rkx , -.1319E+01_rkx , -.1319E+01_rkx , &
    -.1875E+00_rkx , -.1550E+00_rkx , -.1686E+00_rkx , -.1447E+00_rkx , &
    -.2307E+00_rkx , -.6301E+00_rkx , -.4530E+00_rkx , -.4140E+00_rkx , &
    -.4334E+00_rkx , -.3952E+00_rkx , -.2138E+02_rkx , -.2138E+02_rkx , &
    -.2138E+02_rkx , -.2138E+02_rkx , -.2138E+02_rkx , -.2138E+02_rkx , &
    -.2138E+02_rkx , -.3311E+01_rkx , -.3311E+01_rkx , -.1567E+01_rkx , &
    -.1368E+01_rkx , -.1304E+01_rkx , -.1223E+01_rkx , -.1241E+01_rkx , &
    -.1367E+01_rkx , -.1204E+01_rkx , -.1284E+01_rkx , -.1188E+01_rkx , &
    -.1170E+01_rkx , -.2265E+01_rkx , -.2265E+01_rkx , -.2265E+01_rkx , &
    -.2265E+01_rkx , -.2265E+01_rkx , -.2265E+01_rkx , -.2265E+01_rkx , &
    -.2821E-01_rkx , -.2821E-01_rkx , -.4402E+01_rkx , -.1127E+02_rkx , &
    -.2677E+02_rkx , -.2609E+02_rkx , -.4312E+02_rkx , -.4144E+02_rkx , &
    -.3234E+02_rkx , -.4489E+02_rkx , -.3664E+02_rkx , -.4415E+02_rkx , &
    0.5238E+01_rkx , 0.5238E+01_rkx , 0.5238E+01_rkx , 0.5238E+01_rkx , &
    0.5238E+01_rkx , 0.5238E+01_rkx , 0.5238E+01_rkx , 0.8844E+00_rkx , &
    0.8844E+00_rkx , 0.6268E+01_rkx , 0.8334E+01_rkx , 0.1101E+02_rkx , &
    0.8267E+01_rkx , 0.8838E+01_rkx , 0.9620E+01_rkx , 0.8946E+01_rkx , &
    0.9950E+01_rkx , 0.9786E+01_rkx , 0.1031E+02_rkx ], [nspi,ncoefs])

  real(rkx) , dimension(nspi,4) , parameter :: ksdust_stand = reshape(  &
  [ 0.188010E+01_rkx, 0.202540E+01_rkx, 0.195470E+01_rkx, 0.189960E+01_rkx, &
    0.179460E+01_rkx, 0.171490E+01_rkx, 0.154310E+01_rkx, 0.244820E+01_rkx, &
    0.310670E+01_rkx, 0.503910E+00_rkx, 0.503910E+00_rkx, 0.503910E+00_rkx, &
    0.503910E+00_rkx, 0.503910E+00_rkx, 0.503910E+00_rkx, 0.503910E+00_rkx, &
    0.503910E+00_rkx, 0.503910E+00_rkx, 0.503910E+00_rkx, 0.760170E+00_rkx, &
    0.783780E+00_rkx, 0.763890E+00_rkx, 0.749160E+00_rkx, 0.740440E+00_rkx, &
    0.754010E+00_rkx, 0.823220E+00_rkx, 0.856800E+00_rkx, 0.744880E+00_rkx, &
    0.672450E+00_rkx, 0.672450E+00_rkx, 0.672450E+00_rkx, 0.672450E+00_rkx, &
    0.672450E+00_rkx, 0.672450E+00_rkx, 0.672450E+00_rkx, 0.672450E+00_rkx, &
    0.672450E+00_rkx, 0.672450E+00_rkx, 0.366810E+00_rkx, 0.368450E+00_rkx, &
    0.371190E+00_rkx, 0.372200E+00_rkx, 0.369840E+00_rkx, 0.368920E+00_rkx, &
    0.375050E+00_rkx, 0.380780E+00_rkx, 0.436880E+00_rkx, 0.536050E+00_rkx, &
    0.536050E+00_rkx, 0.536050E+00_rkx, 0.536050E+00_rkx, 0.536050E+00_rkx, &
    0.536050E+00_rkx, 0.536050E+00_rkx, 0.536050E+00_rkx, 0.536050E+00_rkx, &
    0.536050E+00_rkx, 0.169330E+00_rkx, 0.170020E+00_rkx, 0.170320E+00_rkx, &
    0.170520E+00_rkx, 0.170720E+00_rkx, 0.170900E+00_rkx, 0.171430E+00_rkx, &
    0.173960E+00_rkx, 0.181040E+00_rkx, 0.205990E+00_rkx, 0.205990E+00_rkx, &
    0.205990E+00_rkx, 0.205990E+00_rkx, 0.205990E+00_rkx, 0.205990E+00_rkx, &
    0.205990E+00_rkx, 0.205990E+00_rkx, 0.205990E+00_rkx, 0.205990E+00_rkx ], &
     [nspi,4])

  real(rkx) , dimension(nspi,4) , parameter :: wsdust_stand = reshape(  &
  [ 0.643280E+00_rkx, 0.677570E+00_rkx, 0.673160E+00_rkx, 0.662450E+00_rkx, &
    0.681320E+00_rkx, 0.679600E+00_rkx, 0.726790E+00_rkx, 0.947300E+00_rkx, &
    0.975360E+00_rkx, 0.895680E+00_rkx, 0.895680E+00_rkx, 0.895680E+00_rkx, &
    0.895680E+00_rkx, 0.895680E+00_rkx, 0.895680E+00_rkx, 0.895680E+00_rkx, &
    0.895680E+00_rkx, 0.895680E+00_rkx, 0.895680E+00_rkx, 0.551960E+00_rkx, &
    0.569090E+00_rkx, 0.560270E+00_rkx, 0.553380E+00_rkx, 0.574400E+00_rkx, &
    0.584670E+00_rkx, 0.687440E+00_rkx, 0.881810E+00_rkx, 0.891610E+00_rkx, &
    0.963220E+00_rkx, 0.963220E+00_rkx, 0.963220E+00_rkx, 0.963220E+00_rkx, &
    0.963220E+00_rkx, 0.963220E+00_rkx, 0.963220E+00_rkx, 0.963220E+00_rkx, &
    0.963220E+00_rkx, 0.963220E+00_rkx, 0.537480E+00_rkx, 0.536390E+00_rkx, &
    0.538750E+00_rkx, 0.539470E+00_rkx, 0.543280E+00_rkx, 0.542420E+00_rkx, &
    0.585640E+00_rkx, 0.807610E+00_rkx, 0.863780E+00_rkx, 0.950080E+00_rkx, &
    0.950080E+00_rkx, 0.950080E+00_rkx, 0.950080E+00_rkx, 0.950080E+00_rkx, &
    0.950080E+00_rkx, 0.950080E+00_rkx, 0.950080E+00_rkx, 0.950080E+00_rkx, &
    0.950080E+00_rkx, 0.543420E+00_rkx, 0.542320E+00_rkx, 0.541810E+00_rkx, &
    0.541490E+00_rkx, 0.541430E+00_rkx, 0.541130E+00_rkx, 0.545760E+00_rkx, &
    0.704550E+00_rkx, 0.758000E+00_rkx, 0.892930E+00_rkx, 0.892930E+00_rkx, &
    0.892930E+00_rkx, 0.892930E+00_rkx, 0.892930E+00_rkx, 0.892930E+00_rkx, &
    0.892930E+00_rkx, 0.892930E+00_rkx, 0.892930E+00_rkx, 0.892930E+00_rkx ], &
     [nspi,4])

  real(rkx) , dimension(nspi,4) , parameter :: gsdust_stand = reshape(  &
  [ 0.871140E+00_rkx, 0.861270E+00_rkx, 0.838000E+00_rkx, 0.817600E+00_rkx, &
    0.770880E+00_rkx, 0.739250E+00_rkx, 0.606950E+00_rkx, 0.643930E+00_rkx, &
    0.747600E+00_rkx, 0.267610E+00_rkx, 0.267610E+00_rkx, 0.267610E+00_rkx, &
    0.267610E+00_rkx, 0.267610E+00_rkx, 0.267610E+00_rkx, 0.267610E+00_rkx, &
    0.267610E+00_rkx, 0.267610E+00_rkx, 0.267610E+00_rkx, 0.925560E+00_rkx, &
    0.921000E+00_rkx, 0.911940E+00_rkx, 0.904420E+00_rkx, 0.885170E+00_rkx, &
    0.883640E+00_rkx, 0.860860E+00_rkx, 0.764570E+00_rkx, 0.620410E+00_rkx, &
    0.560450E+00_rkx, 0.560450E+00_rkx, 0.560450E+00_rkx, 0.560450E+00_rkx, &
    0.560450E+00_rkx, 0.560450E+00_rkx, 0.560450E+00_rkx, 0.560450E+00_rkx, &
    0.560450E+00_rkx, 0.560450E+00_rkx, 0.945420E+00_rkx, 0.943550E+00_rkx, &
    0.943040E+00_rkx, 0.942390E+00_rkx, 0.937100E+00_rkx, 0.935480E+00_rkx, &
    0.918240E+00_rkx, 0.813300E+00_rkx, 0.806650E+00_rkx, 0.688970E+00_rkx, &
    0.688970E+00_rkx, 0.688970E+00_rkx, 0.688970E+00_rkx, 0.688970E+00_rkx, &
    0.688970E+00_rkx, 0.688970E+00_rkx, 0.688970E+00_rkx, 0.688970E+00_rkx, &
    0.688970E+00_rkx, 0.948310E+00_rkx, 0.948130E+00_rkx, 0.948030E+00_rkx, &
    0.947960E+00_rkx, 0.947750E+00_rkx, 0.947630E+00_rkx, 0.944730E+00_rkx, &
    0.877840E+00_rkx, 0.859740E+00_rkx, 0.681740E+00_rkx, 0.681740E+00_rkx, &
    0.681740E+00_rkx, 0.681740E+00_rkx, 0.681740E+00_rkx, 0.681740E+00_rkx, &
    0.681740E+00_rkx, 0.681740E+00_rkx, 0.681740E+00_rkx, 0.681740E+00_rkx ], &
     [nspi,4])

  real(rkx) , dimension(nspi,12) , parameter :: ksdust12_stand = reshape(&
  [ 0.582605E+01_rkx, 0.392185E+01_rkx, 0.334761E+01_rkx, 0.302282E+01_rkx, &
    0.273157E+01_rkx, 0.246940E+01_rkx, 0.188028E+01_rkx, 0.547130E+00_rkx, &
    0.113816E+00_rkx, 0.143850E-01_rkx, 0.143301E-01_rkx, 0.143301E-01_rkx, &
    0.143301E-01_rkx, 0.143301E-01_rkx, 0.142754E-01_rkx, 0.142754E-01_rkx, &
    0.361802E-02_rkx, 0.178064E-02_rkx, 0.178064E-02_rkx, 0.530982E+01_rkx, &
    0.544157E+01_rkx, 0.544018E+01_rkx, 0.543838E+01_rkx, 0.539746E+01_rkx, &
    0.537280E+01_rkx, 0.527157E+01_rkx, 0.357883E+01_rkx, 0.195724E+01_rkx, &
    0.238501E+00_rkx, 0.237798E+00_rkx, 0.237798E+00_rkx, 0.237798E+00_rkx, &
    0.237798E+00_rkx, 0.237098E+00_rkx, 0.237098E+00_rkx, 0.196700E-01_rkx, &
    0.431951E-02_rkx, 0.431951E-02_rkx, 0.144072E+01_rkx, 0.142496E+01_rkx, &
    0.141674E+01_rkx, 0.141367E+01_rkx, 0.141976E+01_rkx, 0.142070E+01_rkx, &
    0.144707E+01_rkx, 0.181699E+01_rkx, 0.210074E+01_rkx, 0.581777E+00_rkx, &
    0.581505E+00_rkx, 0.581505E+00_rkx, 0.581505E+00_rkx, 0.581505E+00_rkx, &
    0.581236E+00_rkx, 0.581236E+00_rkx, 0.240837E+00_rkx, 0.518355E-01_rkx, &
    0.518355E-01_rkx, 0.646088E+00_rkx, 0.651380E+00_rkx, 0.650973E+00_rkx, &
    0.652048E+00_rkx, 0.655378E+00_rkx, 0.660716E+00_rkx, 0.672827E+00_rkx, &
    0.707262E+00_rkx, 0.713763E+00_rkx, 0.612041E+00_rkx, 0.611995E+00_rkx, &
    0.611995E+00_rkx, 0.611995E+00_rkx, 0.611995E+00_rkx, 0.611932E+00_rkx, &
    0.611932E+00_rkx, 0.690188E+00_rkx, 0.220056E+00_rkx, 0.220056E+00_rkx, &
    0.406936E+00_rkx, 0.409548E+00_rkx, 0.410955E+00_rkx, 0.412576E+00_rkx, &
    0.413786E+00_rkx, 0.414107E+00_rkx, 0.415358E+00_rkx, 0.428308E+00_rkx, &
    0.429890E+00_rkx, 0.534167E+00_rkx, 0.534208E+00_rkx, 0.534208E+00_rkx, &
    0.534208E+00_rkx, 0.534208E+00_rkx, 0.534250E+00_rkx, 0.534250E+00_rkx, &
    0.742462E+00_rkx, 0.425718E+00_rkx, 0.425718E+00_rkx, 0.289033E+00_rkx, &
    0.290722E+00_rkx, 0.291565E+00_rkx, 0.291860E+00_rkx, 0.292314E+00_rkx, &
    0.293162E+00_rkx, 0.294211E+00_rkx, 0.302019E+00_rkx, 0.297812E+00_rkx, &
    0.429108E+00_rkx, 0.429167E+00_rkx, 0.429167E+00_rkx, 0.429167E+00_rkx, &
    0.429167E+00_rkx, 0.429172E+00_rkx, 0.429172E+00_rkx, 0.511023E+00_rkx, &
    0.493497E+00_rkx, 0.493497E+00_rkx, 0.233752E+00_rkx, 0.234956E+00_rkx, &
    0.235470E+00_rkx, 0.235824E+00_rkx, 0.236243E+00_rkx, 0.236530E+00_rkx, &
    0.237435E+00_rkx, 0.242130E+00_rkx, 0.248246E+00_rkx, 0.343648E+00_rkx, &
    0.343706E+00_rkx, 0.343706E+00_rkx, 0.343706E+00_rkx, 0.343706E+00_rkx, &
    0.343731E+00_rkx, 0.343731E+00_rkx, 0.291703E+00_rkx, 0.462575E+00_rkx, &
    0.462575E+00_rkx, 0.184485E+00_rkx, 0.185304E+00_rkx, 0.185664E+00_rkx, &
    0.185914E+00_rkx, 0.186151E+00_rkx, 0.186366E+00_rkx, 0.187021E+00_rkx, &
    0.190466E+00_rkx, 0.195955E+00_rkx, 0.249103E+00_rkx, 0.249093E+00_rkx, &
    0.249093E+00_rkx, 0.249093E+00_rkx, 0.249093E+00_rkx, 0.249061E+00_rkx, &
    0.249061E+00_rkx, 0.173126E+00_rkx, 0.340392E+00_rkx, 0.340392E+00_rkx, &
    0.125724E+00_rkx, 0.126168E+00_rkx, 0.126366E+00_rkx, 0.126495E+00_rkx, &
    0.126625E+00_rkx, 0.126751E+00_rkx, 0.127096E+00_rkx, 0.129002E+00_rkx, &
    0.130571E+00_rkx, 0.144168E+00_rkx, 0.144179E+00_rkx, 0.144179E+00_rkx, &
    0.144179E+00_rkx, 0.144179E+00_rkx, 0.144196E+00_rkx, 0.144196E+00_rkx, &
    0.156436E+00_rkx, 0.143665E+00_rkx, 0.143665E+00_rkx, 0.730034E-01_rkx, &
    0.731844E-01_rkx, 0.732652E-01_rkx, 0.733181E-01_rkx, 0.733705E-01_rkx, &
    0.734223E-01_rkx, 0.735618E-01_rkx, 0.743120E-01_rkx, 0.750355E-01_rkx, &
    0.834307E-01_rkx, 0.834394E-01_rkx, 0.834394E-01_rkx, 0.834394E-01_rkx, &
    0.834394E-01_rkx, 0.834493E-01_rkx, 0.834493E-01_rkx, 0.823404E-01_rkx, &
    0.915974E-01_rkx, 0.915974E-01_rkx, 0.100000E-19_rkx, 0.100000E-19_rkx, &
    0.100000E-19_rkx, 0.100000E-19_rkx, 0.434019E-01_rkx, 0.434237E-01_rkx, &
    0.434826E-01_rkx, 0.438195E-01_rkx, 0.441201E-01_rkx, 0.471212E-01_rkx, &
    0.471195E-01_rkx, 0.471195E-01_rkx, 0.471195E-01_rkx, 0.471195E-01_rkx, &
    0.471169E-01_rkx, 0.471169E-01_rkx, 0.455453E-01_rkx, 0.441742E-01_rkx, &
    0.441742E-01_rkx, 0.100000E-19_rkx, 0.100000E-19_rkx, 0.100000E-19_rkx, &
    0.100000E-19_rkx, 0.100000E-19_rkx, 0.100000E-19_rkx, 0.100000E-19_rkx, &
    0.100000E-19_rkx, 0.282240E-01_rkx, 0.290008E-01_rkx, 0.290007E-01_rkx, &
    0.290007E-01_rkx, 0.290007E-01_rkx, 0.290007E-01_rkx, 0.290011E-01_rkx, &
    0.290011E-01_rkx, 0.295979E-01_rkx, 0.300372E-01_rkx, 0.300372E-01_rkx ], &
     [nspi,12])

  real(rkx) , dimension(nspi,12) , parameter :: wsdust12_stand = reshape(&
  [ 0.803780E+00_rkx, 0.799301E+00_rkx, 0.797433E+00_rkx, 0.795686E+00_rkx, &
    0.793633E+00_rkx, 0.791376E+00_rkx, 0.791539E+00_rkx, 0.816759E+00_rkx, &
    0.833502E+00_rkx, 0.208783E+00_rkx, 0.208637E+00_rkx, 0.208637E+00_rkx, &
    0.208637E+00_rkx, 0.208637E+00_rkx, 0.208492E+00_rkx, 0.208492E+00_rkx, &
    0.914436E-01_rkx, 0.293192E-01_rkx, 0.293192E-01_rkx, 0.736897E+00_rkx, &
    0.794829E+00_rkx, 0.820288E+00_rkx, 0.833319E+00_rkx, 0.847253E+00_rkx, &
    0.860061E+00_rkx, 0.885778E+00_rkx, 0.952989E+00_rkx, 0.978783E+00_rkx, &
    0.706520E+00_rkx, 0.706464E+00_rkx, 0.706464E+00_rkx, 0.706464E+00_rkx, &
    0.706464E+00_rkx, 0.706408E+00_rkx, 0.706408E+00_rkx, 0.718841E+00_rkx, &
    0.495389E+00_rkx, 0.495389E+00_rkx, 0.578261E+00_rkx, 0.600564E+00_rkx, &
    0.613124E+00_rkx, 0.623304E+00_rkx, 0.633117E+00_rkx, 0.643199E+00_rkx, &
    0.682406E+00_rkx, 0.876231E+00_rkx, 0.971464E+00_rkx, 0.957579E+00_rkx, &
    0.957584E+00_rkx, 0.957584E+00_rkx, 0.957584E+00_rkx, 0.957584E+00_rkx, &
    0.957590E+00_rkx, 0.957590E+00_rkx, 0.969097E+00_rkx, 0.928109E+00_rkx, &
    0.928109E+00_rkx, 0.536202E+00_rkx, 0.545354E+00_rkx, 0.550490E+00_rkx, &
    0.555428E+00_rkx, 0.562435E+00_rkx, 0.571515E+00_rkx, 0.605722E+00_rkx, &
    0.807373E+00_rkx, 0.937056E+00_rkx, 0.979975E+00_rkx, 0.979988E+00_rkx, &
    0.979988E+00_rkx, 0.979988E+00_rkx, 0.979988E+00_rkx, 0.979998E+00_rkx, &
    0.979998E+00_rkx, 0.988435E+00_rkx, 0.984462E+00_rkx, 0.984462E+00_rkx, &
    0.534536E+00_rkx, 0.535369E+00_rkx, 0.537324E+00_rkx, 0.540009E+00_rkx, &
    0.542884E+00_rkx, 0.545557E+00_rkx, 0.561734E+00_rkx, 0.759934E+00_rkx, &
    0.909948E+00_rkx, 0.973860E+00_rkx, 0.973911E+00_rkx, 0.973911E+00_rkx, &
    0.973911E+00_rkx, 0.973911E+00_rkx, 0.973929E+00_rkx, 0.973929E+00_rkx, &
    0.988239E+00_rkx, 0.990160E+00_rkx, 0.990160E+00_rkx, 0.537743E+00_rkx, &
    0.536572E+00_rkx, 0.536678E+00_rkx, 0.536671E+00_rkx, 0.537174E+00_rkx, &
    0.538667E+00_rkx, 0.546270E+00_rkx, 0.723814E+00_rkx, 0.885800E+00_rkx, &
    0.964989E+00_rkx, 0.965026E+00_rkx, 0.965026E+00_rkx, 0.965026E+00_rkx, &
    0.965026E+00_rkx, 0.965149E+00_rkx, 0.965149E+00_rkx, 0.982017E+00_rkx, &
    0.990949E+00_rkx, 0.990949E+00_rkx, 0.540048E+00_rkx, 0.538641E+00_rkx, &
    0.538149E+00_rkx, 0.537987E+00_rkx, 0.538106E+00_rkx, 0.538213E+00_rkx, &
    0.542012E+00_rkx, 0.700613E+00_rkx, 0.870473E+00_rkx, 0.957201E+00_rkx, &
    0.957130E+00_rkx, 0.957130E+00_rkx, 0.957130E+00_rkx, 0.957130E+00_rkx, &
    0.957123E+00_rkx, 0.957123E+00_rkx, 0.967379E+00_rkx, 0.989883E+00_rkx, &
    0.989883E+00_rkx, 0.542302E+00_rkx, 0.541019E+00_rkx, 0.540465E+00_rkx, &
    0.540171E+00_rkx, 0.539910E+00_rkx, 0.539690E+00_rkx, 0.540825E+00_rkx, &
    0.675834E+00_rkx, 0.849436E+00_rkx, 0.946153E+00_rkx, 0.946201E+00_rkx, &
    0.946201E+00_rkx, 0.946201E+00_rkx, 0.946201E+00_rkx, 0.946308E+00_rkx, &
    0.946308E+00_rkx, 0.946405E+00_rkx, 0.985447E+00_rkx, 0.985447E+00_rkx, &
    0.545056E+00_rkx, 0.544150E+00_rkx, 0.543734E+00_rkx, 0.543458E+00_rkx, &
    0.543194E+00_rkx, 0.542933E+00_rkx, 0.542557E+00_rkx, 0.638271E+00_rkx, &
    0.800642E+00_rkx, 0.924151E+00_rkx, 0.924183E+00_rkx, 0.924183E+00_rkx, &
    0.924183E+00_rkx, 0.924183E+00_rkx, 0.924198E+00_rkx, 0.924198E+00_rkx, &
    0.941366E+00_rkx, 0.963209E+00_rkx, 0.963209E+00_rkx, 0.547463E+00_rkx, &
    0.546961E+00_rkx, 0.546725E+00_rkx, 0.546566E+00_rkx, 0.546407E+00_rkx, &
    0.546247E+00_rkx, 0.545805E+00_rkx, 0.595887E+00_rkx, 0.726096E+00_rkx, &
    0.890923E+00_rkx, 0.890979E+00_rkx, 0.890979E+00_rkx, 0.890979E+00_rkx, &
    0.890979E+00_rkx, 0.891029E+00_rkx, 0.891029E+00_rkx, 0.912332E+00_rkx, &
    0.947218E+00_rkx, 0.947218E+00_rkx, 0.500000E+00_rkx, 0.500000E+00_rkx, &
    0.500000E+00_rkx, 0.500000E+00_rkx, 0.548104E+00_rkx, 0.548018E+00_rkx, &
    0.547773E+00_rkx, 0.568268E+00_rkx, 0.648635E+00_rkx, 0.846069E+00_rkx, &
    0.846132E+00_rkx, 0.846132E+00_rkx, 0.846132E+00_rkx, 0.846132E+00_rkx, &
    0.846194E+00_rkx, 0.846194E+00_rkx, 0.866683E+00_rkx, 0.914608E+00_rkx, &
    0.914608E+00_rkx, 0.500000E+00_rkx, 0.500000E+00_rkx, 0.500000E+00_rkx, &
    0.500000E+00_rkx, 0.500000E+00_rkx, 0.500000E+00_rkx, 0.500000E+00_rkx, &
    0.500000E+00_rkx, 0.594727E+00_rkx, 0.795945E+00_rkx, 0.796001E+00_rkx, &
    0.796001E+00_rkx, 0.796001E+00_rkx, 0.796001E+00_rkx, 0.796056E+00_rkx, &
    0.796056E+00_rkx, 0.822320E+00_rkx, 0.888785E+00_rkx, 0.888785E+00_rkx ], &
     [nspi,12])

  real(rkx) , dimension(nspi,12) , parameter :: gsdust12_stand = reshape(&
  [ 0.546546E+00_rkx, 0.415561E+00_rkx, 0.361926E+00_rkx, 0.331551E+00_rkx, &
    0.304818E+00_rkx, 0.281297E+00_rkx, 0.231473E+00_rkx, 0.108432E+00_rkx, &
    0.528125E-01_rkx, 0.910449E-02_rkx, 0.908684E-02_rkx, 0.908684E-02_rkx, &
    0.908684E-02_rkx, 0.908684E-02_rkx, 0.906926E-02_rkx, 0.906926E-02_rkx, &
    0.319155E-02_rkx, 0.126936E-02_rkx, 0.126936E-02_rkx, 0.719875E+00_rkx, &
    0.720841E+00_rkx, 0.726308E+00_rkx, 0.721714E+00_rkx, 0.716225E+00_rkx, &
    0.716554E+00_rkx, 0.718668E+00_rkx, 0.628714E+00_rkx, 0.506675E+00_rkx, &
    0.101909E+00_rkx, 0.101787E+00_rkx, 0.101787E+00_rkx, 0.101787E+00_rkx, &
    0.101787E+00_rkx, 0.101666E+00_rkx, 0.101666E+00_rkx, 0.386167E-01_rkx, &
    0.154572E-01_rkx, 0.154572E-01_rkx, 0.886115E+00_rkx, 0.850893E+00_rkx, &
    0.831440E+00_rkx, 0.822061E+00_rkx, 0.810844E+00_rkx, 0.795068E+00_rkx, &
    0.766308E+00_rkx, 0.666552E+00_rkx, 0.650917E+00_rkx, 0.345008E+00_rkx, &
    0.344985E+00_rkx, 0.344985E+00_rkx, 0.344985E+00_rkx, 0.344985E+00_rkx, &
    0.344961E+00_rkx, 0.344961E+00_rkx, 0.312532E+00_rkx, 0.116778E+00_rkx, &
    0.116778E+00_rkx, 0.936182E+00_rkx, 0.927581E+00_rkx, 0.921453E+00_rkx, &
    0.917115E+00_rkx, 0.912670E+00_rkx, 0.908398E+00_rkx, 0.891809E+00_rkx, &
    0.791368E+00_rkx, 0.686444E+00_rkx, 0.586731E+00_rkx, 0.586594E+00_rkx, &
    0.586594E+00_rkx, 0.586594E+00_rkx, 0.586594E+00_rkx, 0.586444E+00_rkx, &
    0.586444E+00_rkx, 0.654273E+00_rkx, 0.451257E+00_rkx, 0.451257E+00_rkx, &
    0.944853E+00_rkx, 0.942077E+00_rkx, 0.940080E+00_rkx, 0.938567E+00_rkx, &
    0.936667E+00_rkx, 0.934233E+00_rkx, 0.923936E+00_rkx, 0.827524E+00_rkx, &
    0.756116E+00_rkx, 0.670752E+00_rkx, 0.670696E+00_rkx, 0.670696E+00_rkx, &
    0.670696E+00_rkx, 0.670696E+00_rkx, 0.670656E+00_rkx, 0.670656E+00_rkx, &
    0.729345E+00_rkx, 0.645941E+00_rkx, 0.645941E+00_rkx, 0.947046E+00_rkx, &
    0.946199E+00_rkx, 0.945564E+00_rkx, 0.944893E+00_rkx, 0.944115E+00_rkx, &
    0.943280E+00_rkx, 0.938002E+00_rkx, 0.857309E+00_rkx, 0.771482E+00_rkx, &
    0.691530E+00_rkx, 0.691508E+00_rkx, 0.691508E+00_rkx, 0.691508E+00_rkx, &
    0.691508E+00_rkx, 0.691547E+00_rkx, 0.691547E+00_rkx, 0.694262E+00_rkx, &
    0.722575E+00_rkx, 0.722575E+00_rkx, 0.947651E+00_rkx, 0.947249E+00_rkx, &
    0.946937E+00_rkx, 0.946655E+00_rkx, 0.946317E+00_rkx, 0.945823E+00_rkx, &
    0.942902E+00_rkx, 0.871887E+00_rkx, 0.802731E+00_rkx, 0.688915E+00_rkx, &
    0.688807E+00_rkx, 0.688807E+00_rkx, 0.688807E+00_rkx, 0.688807E+00_rkx, &
    0.688747E+00_rkx, 0.688747E+00_rkx, 0.571069E+00_rkx, 0.737631E+00_rkx, &
    0.737631E+00_rkx, 0.948064E+00_rkx, 0.947862E+00_rkx, 0.947729E+00_rkx, &
    0.947619E+00_rkx, 0.947471E+00_rkx, 0.947269E+00_rkx, 0.945948E+00_rkx, &
    0.887280E+00_rkx, 0.819214E+00_rkx, 0.680801E+00_rkx, 0.680825E+00_rkx, &
    0.680825E+00_rkx, 0.680825E+00_rkx, 0.680825E+00_rkx, 0.680876E+00_rkx, &
    0.680876E+00_rkx, 0.511984E+00_rkx, 0.698317E+00_rkx, 0.698317E+00_rkx, &
    0.948455E+00_rkx, 0.948376E+00_rkx, 0.948332E+00_rkx, 0.948297E+00_rkx, &
    0.948258E+00_rkx, 0.948210E+00_rkx, 0.947904E+00_rkx, 0.907566E+00_rkx, &
    0.844326E+00_rkx, 0.694185E+00_rkx, 0.694119E+00_rkx, 0.694119E+00_rkx, &
    0.694119E+00_rkx, 0.694119E+00_rkx, 0.694035E+00_rkx, 0.694035E+00_rkx, &
    0.742725E+00_rkx, 0.561306E+00_rkx, 0.561306E+00_rkx, 0.948708E+00_rkx, &
    0.948709E+00_rkx, 0.948704E+00_rkx, 0.948700E+00_rkx, 0.948695E+00_rkx, &
    0.948688E+00_rkx, 0.948668E+00_rkx, 0.927961E+00_rkx, 0.879297E+00_rkx, &
    0.786935E+00_rkx, 0.786899E+00_rkx, 0.786899E+00_rkx, 0.786899E+00_rkx, &
    0.786899E+00_rkx, 0.786855E+00_rkx, 0.786855E+00_rkx, 0.768741E+00_rkx, &
    0.759343E+00_rkx, 0.759343E+00_rkx, 0.990000E+00_rkx, 0.990000E+00_rkx, &
    0.990000E+00_rkx, 0.990000E+00_rkx, 0.948823E+00_rkx, 0.948829E+00_rkx, &
    0.948848E+00_rkx, 0.940308E+00_rkx, 0.910074E+00_rkx, 0.816865E+00_rkx, &
    0.816810E+00_rkx, 0.816810E+00_rkx, 0.816810E+00_rkx, 0.816810E+00_rkx, &
    0.816755E+00_rkx, 0.816755E+00_rkx, 0.804314E+00_rkx, 0.744671E+00_rkx, &
    0.744671E+00_rkx, 0.990000E+00_rkx, 0.990000E+00_rkx, 0.990000E+00_rkx, &
    0.990000E+00_rkx, 0.990000E+00_rkx, 0.990000E+00_rkx, 0.990000E+00_rkx, &
    0.990000E+00_rkx, 0.930485E+00_rkx, 0.842174E+00_rkx, 0.842163E+00_rkx, &
    0.842163E+00_rkx, 0.842163E+00_rkx, 0.842163E+00_rkx, 0.842159E+00_rkx, &
    0.842159E+00_rkx, 0.840731E+00_rkx, 0.801060E+00_rkx, 0.801060E+00_rkx ], &
     [nspi,12])

  real(rkx) , dimension(nwav,2,nih) , parameter :: ksslt = reshape( &
  [ 0.110670E+01_rkx, 0.113470E+01_rkx, 0.114490E+01_rkx, 0.115360E+01_rkx, &
    0.116140E+01_rkx, 0.116810E+01_rkx, 0.118630E+01_rkx, 0.126130E+01_rkx, &
    0.128830E+01_rkx, 0.111740E+01_rkx, 0.107340E+01_rkx, 0.101450E+01_rkx, &
    0.937800E+00_rkx, 0.887000E+00_rkx, 0.857100E+00_rkx, 0.772600E+00_rkx, &
    0.551900E+00_rkx, 0.236900E+00_rkx, 0.236900E+00_rkx, 0.243500E+00_rkx, &
    0.244900E+00_rkx, 0.245300E+00_rkx, 0.245700E+00_rkx, 0.246100E+00_rkx, &
    0.246400E+00_rkx, 0.247300E+00_rkx, 0.251200E+00_rkx, 0.256100E+00_rkx, &
    0.264800E+00_rkx, 0.267600E+00_rkx, 0.270700E+00_rkx, 0.274000E+00_rkx, &
    0.276100E+00_rkx, 0.277500E+00_rkx, 0.281700E+00_rkx, 0.296200E+00_rkx, &
    0.318800E+00_rkx, 0.318800E+00_rkx, 0.276050E+01_rkx, 0.279270E+01_rkx, &
    0.280620E+01_rkx, 0.281920E+01_rkx, 0.283250E+01_rkx, 0.284560E+01_rkx, &
    0.288990E+01_rkx, 0.308680E+01_rkx, 0.320500E+01_rkx, 0.299350E+01_rkx, &
    0.291440E+01_rkx, 0.283310E+01_rkx, 0.269720E+01_rkx, 0.260410E+01_rkx, &
    0.253200E+01_rkx, 0.231690E+01_rkx, 0.174510E+01_rkx, 0.960700E+00_rkx, &
    0.960700E+00_rkx, 0.625000E+00_rkx, 0.627400E+00_rkx, 0.628200E+00_rkx, &
    0.629000E+00_rkx, 0.629700E+00_rkx, 0.630400E+00_rkx, 0.632200E+00_rkx, &
    0.640200E+00_rkx, 0.648300E+00_rkx, 0.665400E+00_rkx, 0.670000E+00_rkx, &
    0.674600E+00_rkx, 0.680400E+00_rkx, 0.684400E+00_rkx, 0.687400E+00_rkx, &
    0.695800E+00_rkx, 0.714900E+00_rkx, 0.779500E+00_rkx, 0.779500E+00_rkx, &
    0.343580E+01_rkx, 0.348660E+01_rkx, 0.350620E+01_rkx, 0.352400E+01_rkx, &
    0.354090E+01_rkx, 0.355630E+01_rkx, 0.360340E+01_rkx, 0.384490E+01_rkx, &
    0.401160E+01_rkx, 0.383520E+01_rkx, 0.375150E+01_rkx, 0.367460E+01_rkx, &
    0.353150E+01_rkx, 0.343520E+01_rkx, 0.335000E+01_rkx, 0.309100E+01_rkx, &
    0.236570E+01_rkx, 0.139570E+01_rkx, 0.139570E+01_rkx, 0.789700E+00_rkx, &
    0.792100E+00_rkx, 0.793000E+00_rkx, 0.793700E+00_rkx, 0.794400E+00_rkx, &
    0.795100E+00_rkx, 0.797000E+00_rkx, 0.806000E+00_rkx, 0.814900E+00_rkx, &
    0.835600E+00_rkx, 0.840900E+00_rkx, 0.846300E+00_rkx, 0.852400E+00_rkx, &
    0.857000E+00_rkx, 0.860300E+00_rkx, 0.869800E+00_rkx, 0.890800E+00_rkx, &
    0.969800E+00_rkx, 0.969800E+00_rkx, 0.413430E+01_rkx, 0.418210E+01_rkx, &
    0.420170E+01_rkx, 0.421980E+01_rkx, 0.423730E+01_rkx, 0.425330E+01_rkx, &
    0.430420E+01_rkx, 0.459390E+01_rkx, 0.480120E+01_rkx, 0.468220E+01_rkx, &
    0.459660E+01_rkx, 0.452830E+01_rkx, 0.438560E+01_rkx, 0.429240E+01_rkx, &
    0.419710E+01_rkx, 0.390190E+01_rkx, 0.303000E+01_rkx, 0.189030E+01_rkx, &
    0.189030E+01_rkx, 0.952200E+00_rkx, 0.956900E+00_rkx, 0.958500E+00_rkx, &
    0.959800E+00_rkx, 0.960800E+00_rkx, 0.961600E+00_rkx, 0.963200E+00_rkx, &
    0.973700E+00_rkx, 0.983000E+00_rkx, 0.100670E+01_rkx, 0.101270E+01_rkx, &
    0.101870E+01_rkx, 0.102600E+01_rkx, 0.103140E+01_rkx, 0.103510E+01_rkx, &
    0.104550E+01_rkx, 0.106800E+01_rkx, 0.115860E+01_rkx, 0.115860E+01_rkx, &
    0.584240E+01_rkx, 0.588300E+01_rkx, 0.590140E+01_rkx, 0.591970E+01_rkx, &
    0.593840E+01_rkx, 0.595670E+01_rkx, 0.602000E+01_rkx, 0.638780E+01_rkx, &
    0.669990E+01_rkx, 0.673060E+01_rkx, 0.665350E+01_rkx, 0.661700E+01_rkx, &
    0.649700E+01_rkx, 0.642870E+01_rkx, 0.631860E+01_rkx, 0.596400E+01_rkx, &
    0.477770E+01_rkx, 0.329640E+01_rkx, 0.329640E+01_rkx, 0.136500E+01_rkx, &
    0.136820E+01_rkx, 0.136940E+01_rkx, 0.137050E+01_rkx, 0.137140E+01_rkx, &
    0.137230E+01_rkx, 0.137490E+01_rkx, 0.138870E+01_rkx, 0.140120E+01_rkx, &
    0.143000E+01_rkx, 0.143720E+01_rkx, 0.144470E+01_rkx, 0.145260E+01_rkx, &
    0.145840E+01_rkx, 0.146280E+01_rkx, 0.147540E+01_rkx, 0.150540E+01_rkx, &
    0.161800E+01_rkx, 0.161800E+01_rkx, 0.847370E+01_rkx, 0.857490E+01_rkx, &
    0.861230E+01_rkx, 0.864390E+01_rkx, 0.867140E+01_rkx, 0.869340E+01_rkx, &
    0.874920E+01_rkx, 0.921250E+01_rkx, 0.961140E+01_rkx, 0.994990E+01_rkx, &
    0.990040E+01_rkx, 0.992670E+01_rkx, 0.987340E+01_rkx, 0.986950E+01_rkx, &
    0.975570E+01_rkx, 0.936830E+01_rkx, 0.780820E+01_rkx, 0.596950E+01_rkx, &
    0.596950E+01_rkx, 0.203340E+01_rkx, 0.203540E+01_rkx, 0.203640E+01_rkx, &
    0.203730E+01_rkx, 0.203810E+01_rkx, 0.203880E+01_rkx, 0.204110E+01_rkx, &
    0.206220E+01_rkx, 0.207770E+01_rkx, 0.211370E+01_rkx, 0.212310E+01_rkx, &
    0.213290E+01_rkx, 0.214430E+01_rkx, 0.215220E+01_rkx, 0.215760E+01_rkx, &
    0.217230E+01_rkx, 0.220780E+01_rkx, 0.234700E+01_rkx, 0.234700E+01_rkx, &
    0.146884E+02_rkx, 0.148200E+02_rkx, 0.148680E+02_rkx, 0.149089E+02_rkx, &
    0.149452E+02_rkx, 0.149753E+02_rkx, 0.150547E+02_rkx, 0.156303E+02_rkx, &
    0.161627E+02_rkx, 0.170787E+02_rkx, 0.171301E+02_rkx, 0.172957E+02_rkx, &
    0.174347E+02_rkx, 0.176012E+02_rkx, 0.175259E+02_rkx, 0.172179E+02_rkx, &
    0.152499E+02_rkx, 0.132278E+02_rkx, 0.132278E+02_rkx, 0.360330E+01_rkx, &
    0.360900E+01_rkx, 0.361130E+01_rkx, 0.361330E+01_rkx, 0.361520E+01_rkx, &
    0.361700E+01_rkx, 0.362240E+01_rkx, 0.364850E+01_rkx, 0.366570E+01_rkx, &
    0.372400E+01_rkx, 0.373710E+01_rkx, 0.375280E+01_rkx, 0.376940E+01_rkx, &
    0.378100E+01_rkx, 0.378800E+01_rkx, 0.380740E+01_rkx, 0.385640E+01_rkx, &
    0.403860E+01_rkx, 0.403860E+01_rkx, 0.224888E+02_rkx, 0.226277E+02_rkx, &
    0.226797E+02_rkx, 0.227257E+02_rkx, 0.227683E+02_rkx, 0.228059E+02_rkx, &
    0.229152E+02_rkx, 0.236106E+02_rkx, 0.242436E+02_rkx, 0.257450E+02_rkx, &
    0.259205E+02_rkx, 0.262202E+02_rkx, 0.265736E+02_rkx, 0.269240E+02_rkx, &
    0.269291E+02_rkx, 0.268394E+02_rkx, 0.248584E+02_rkx, 0.233536E+02_rkx, &
    0.233536E+02_rkx, 0.562270E+01_rkx, 0.561380E+01_rkx, 0.561090E+01_rkx, &
    0.560950E+01_rkx, 0.560980E+01_rkx, 0.561170E+01_rkx, 0.562400E+01_rkx, &
    0.566380E+01_rkx, 0.569990E+01_rkx, 0.576950E+01_rkx, 0.578650E+01_rkx, &
    0.580620E+01_rkx, 0.582770E+01_rkx, 0.584240E+01_rkx, 0.585120E+01_rkx, &
    0.587590E+01_rkx, 0.594260E+01_rkx, 0.616670E+01_rkx, 0.616670E+01_rkx ], &
     [nwav,2,nih])

  real(rkx) , dimension(nwav,2,nih) , parameter :: wsslt = reshape( &
  [ 0.999700E+00_rkx, 0.999800E+00_rkx, 0.999800E+00_rkx, 0.999900E+00_rkx, &
    0.999900E+00_rkx, 0.999900E+00_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, &
    0.100000E+01_rkx, 0.997800E+00_rkx, 0.997000E+00_rkx, 0.994800E+00_rkx, &
    0.992600E+00_rkx, 0.992500E+00_rkx, 0.991100E+00_rkx, 0.986900E+00_rkx, &
    0.956000E+00_rkx, 0.989700E+00_rkx, 0.989700E+00_rkx, 0.997800E+00_rkx, &
    0.998700E+00_rkx, 0.999000E+00_rkx, 0.999200E+00_rkx, 0.999400E+00_rkx, &
    0.999500E+00_rkx, 0.999800E+00_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, &
    0.987600E+00_rkx, 0.984500E+00_rkx, 0.975900E+00_rkx, 0.967100E+00_rkx, &
    0.965200E+00_rkx, 0.960800E+00_rkx, 0.947900E+00_rkx, 0.855100E+00_rkx, &
    0.978200E+00_rkx, 0.978200E+00_rkx, 0.999800E+00_rkx, 0.999900E+00_rkx, &
    0.100000E+01_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, &
    0.100000E+01_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, 0.998900E+00_rkx, &
    0.996100E+00_rkx, 0.975800E+00_rkx, 0.960600E+00_rkx, 0.968600E+00_rkx, &
    0.967300E+00_rkx, 0.959700E+00_rkx, 0.706400E+00_rkx, 0.933700E+00_rkx, &
    0.933700E+00_rkx, 0.999300E+00_rkx, 0.999600E+00_rkx, 0.999600E+00_rkx, &
    0.999700E+00_rkx, 0.999800E+00_rkx, 0.999800E+00_rkx, 0.999900E+00_rkx, &
    0.100000E+01_rkx, 0.100000E+01_rkx, 0.991000E+00_rkx, 0.983800E+00_rkx, &
    0.960400E+00_rkx, 0.943400E+00_rkx, 0.949000E+00_rkx, 0.945900E+00_rkx, &
    0.934700E+00_rkx, 0.717500E+00_rkx, 0.854700E+00_rkx, 0.854700E+00_rkx, &
    0.999800E+00_rkx, 0.999900E+00_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, &
    0.100000E+01_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, &
    0.100000E+01_rkx, 0.999000E+00_rkx, 0.996000E+00_rkx, 0.975100E+00_rkx, &
    0.959600E+00_rkx, 0.967900E+00_rkx, 0.966700E+00_rkx, 0.959200E+00_rkx, &
    0.701600E+00_rkx, 0.929400E+00_rkx, 0.929400E+00_rkx, 0.999400E+00_rkx, &
    0.999700E+00_rkx, 0.999700E+00_rkx, 0.999800E+00_rkx, 0.999900E+00_rkx, &
    0.999900E+00_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, &
    0.991000E+00_rkx, 0.983200E+00_rkx, 0.959800E+00_rkx, 0.942800E+00_rkx, &
    0.948200E+00_rkx, 0.945100E+00_rkx, 0.933700E+00_rkx, 0.717700E+00_rkx, &
    0.833600E+00_rkx, 0.833600E+00_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, &
    0.100000E+01_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, &
    0.100000E+01_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, 0.999000E+00_rkx, &
    0.995900E+00_rkx, 0.974800E+00_rkx, 0.959200E+00_rkx, 0.967600E+00_rkx, &
    0.966400E+00_rkx, 0.959100E+00_rkx, 0.701000E+00_rkx, 0.927200E+00_rkx, &
    0.927200E+00_rkx, 0.999400E+00_rkx, 0.999700E+00_rkx, 0.999700E+00_rkx, &
    0.999800E+00_rkx, 0.999900E+00_rkx, 0.999900E+00_rkx, 0.100000E+01_rkx, &
    0.100000E+01_rkx, 0.100000E+01_rkx, 0.990700E+00_rkx, 0.982500E+00_rkx, &
    0.959200E+00_rkx, 0.942100E+00_rkx, 0.947300E+00_rkx, 0.944100E+00_rkx, &
    0.932300E+00_rkx, 0.717200E+00_rkx, 0.817200E+00_rkx, 0.817200E+00_rkx, &
    0.100000E+01_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, &
    0.100000E+01_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, &
    0.100000E+01_rkx, 0.999000E+00_rkx, 0.995800E+00_rkx, 0.974400E+00_rkx, &
    0.958800E+00_rkx, 0.967300E+00_rkx, 0.966200E+00_rkx, 0.959200E+00_rkx, &
    0.703500E+00_rkx, 0.924600E+00_rkx, 0.924600E+00_rkx, 0.999600E+00_rkx, &
    0.999700E+00_rkx, 0.999800E+00_rkx, 0.999800E+00_rkx, 0.999900E+00_rkx, &
    0.999900E+00_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, &
    0.990100E+00_rkx, 0.981000E+00_rkx, 0.957900E+00_rkx, 0.940500E+00_rkx, &
    0.945200E+00_rkx, 0.941700E+00_rkx, 0.928800E+00_rkx, 0.715000E+00_rkx, &
    0.787600E+00_rkx, 0.787600E+00_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, &
    0.100000E+01_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, &
    0.100000E+01_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, 0.998900E+00_rkx, &
    0.995500E+00_rkx, 0.974000E+00_rkx, 0.958600E+00_rkx, 0.967000E+00_rkx, &
    0.966000E+00_rkx, 0.959300E+00_rkx, 0.708700E+00_rkx, 0.922400E+00_rkx, &
    0.922400E+00_rkx, 0.999700E+00_rkx, 0.999800E+00_rkx, 0.999800E+00_rkx, &
    0.999900E+00_rkx, 0.999900E+00_rkx, 0.999900E+00_rkx, 0.100000E+01_rkx, &
    0.100000E+01_rkx, 0.100000E+01_rkx, 0.988900E+00_rkx, 0.978800E+00_rkx, &
    0.956100E+00_rkx, 0.937900E+00_rkx, 0.941800E+00_rkx, 0.937700E+00_rkx, &
    0.923200E+00_rkx, 0.710800E+00_rkx, 0.755700E+00_rkx, 0.755700E+00_rkx, &
    0.100000E+01_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, &
    0.100000E+01_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, &
    0.100000E+01_rkx, 0.998600E+00_rkx, 0.995000E+00_rkx, 0.973300E+00_rkx, &
    0.958000E+00_rkx, 0.966300E+00_rkx, 0.965300E+00_rkx, 0.959000E+00_rkx, &
    0.716700E+00_rkx, 0.917900E+00_rkx, 0.917900E+00_rkx, 0.999800E+00_rkx, &
    0.999900E+00_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, &
    0.100000E+01_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, &
    0.986600E+00_rkx, 0.975300E+00_rkx, 0.952800E+00_rkx, 0.933500E+00_rkx, &
    0.935900E+00_rkx, 0.930800E+00_rkx, 0.913200E+00_rkx, 0.701600E+00_rkx, &
    0.712900E+00_rkx, 0.712900E+00_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, &
    0.100000E+01_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, &
    0.100000E+01_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, 0.998300E+00_rkx, &
    0.994300E+00_rkx, 0.972400E+00_rkx, 0.957200E+00_rkx, 0.965400E+00_rkx, &
    0.964400E+00_rkx, 0.958200E+00_rkx, 0.722100E+00_rkx, 0.911800E+00_rkx, &
    0.911800E+00_rkx, 0.999800E+00_rkx, 0.999900E+00_rkx, 0.100000E+01_rkx, &
    0.100000E+01_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, 0.100000E+01_rkx, &
    0.100000E+01_rkx, 0.100000E+01_rkx, 0.984600E+00_rkx, 0.972100E+00_rkx, &
    0.949800E+00_rkx, 0.929300E+00_rkx, 0.930500E+00_rkx, 0.924400E+00_rkx, &
    0.904200E+00_rkx, 0.692800E+00_rkx, 0.682700E+00_rkx, 0.682700E+00_rkx ], &
     [nwav,2,nih])

  real(rkx) , dimension(nwav,2,nih) , parameter :: gsslt = reshape( &
  [ 0.730800E+00_rkx, 0.718200E+00_rkx, 0.714200E+00_rkx, 0.710900E+00_rkx, &
    0.708200E+00_rkx, 0.706000E+00_rkx, 0.701200E+00_rkx, 0.695100E+00_rkx, &
    0.695900E+00_rkx, 0.701400E+00_rkx, 0.699400E+00_rkx, 0.697100E+00_rkx, &
    0.694700E+00_rkx, 0.696000E+00_rkx, 0.695300E+00_rkx, 0.692200E+00_rkx, &
    0.640500E+00_rkx, 0.594900E+00_rkx, 0.594900E+00_rkx, 0.808400E+00_rkx, &
    0.810400E+00_rkx, 0.810800E+00_rkx, 0.810900E+00_rkx, 0.810600E+00_rkx, &
    0.810000E+00_rkx, 0.806300E+00_rkx, 0.799500E+00_rkx, 0.790200E+00_rkx, &
    0.777600E+00_rkx, 0.773900E+00_rkx, 0.769200E+00_rkx, 0.764300E+00_rkx, &
    0.761800E+00_rkx, 0.761300E+00_rkx, 0.759800E+00_rkx, 0.741800E+00_rkx, &
    0.705900E+00_rkx, 0.705900E+00_rkx, 0.785300E+00_rkx, 0.784900E+00_rkx, &
    0.784600E+00_rkx, 0.784000E+00_rkx, 0.783100E+00_rkx, 0.782000E+00_rkx, &
    0.777000E+00_rkx, 0.771800E+00_rkx, 0.773100E+00_rkx, 0.778600E+00_rkx, &
    0.777700E+00_rkx, 0.777600E+00_rkx, 0.778900E+00_rkx, 0.781100E+00_rkx, &
    0.782200E+00_rkx, 0.784700E+00_rkx, 0.775900E+00_rkx, 0.705300E+00_rkx, &
    0.705300E+00_rkx, 0.844900E+00_rkx, 0.846300E+00_rkx, 0.846800E+00_rkx, &
    0.847200E+00_rkx, 0.847600E+00_rkx, 0.848000E+00_rkx, 0.849300E+00_rkx, &
    0.848900E+00_rkx, 0.843900E+00_rkx, 0.834500E+00_rkx, 0.833300E+00_rkx, &
    0.835900E+00_rkx, 0.836400E+00_rkx, 0.832600E+00_rkx, 0.832800E+00_rkx, &
    0.834300E+00_rkx, 0.886900E+00_rkx, 0.806500E+00_rkx, 0.806500E+00_rkx, &
    0.803000E+00_rkx, 0.796500E+00_rkx, 0.794400E+00_rkx, 0.792600E+00_rkx, &
    0.791200E+00_rkx, 0.790000E+00_rkx, 0.787500E+00_rkx, 0.780100E+00_rkx, &
    0.779900E+00_rkx, 0.786500E+00_rkx, 0.785800E+00_rkx, 0.786100E+00_rkx, &
    0.788100E+00_rkx, 0.790600E+00_rkx, 0.792000E+00_rkx, 0.795700E+00_rkx, &
    0.793900E+00_rkx, 0.724500E+00_rkx, 0.724500E+00_rkx, 0.848800E+00_rkx, &
    0.847900E+00_rkx, 0.847900E+00_rkx, 0.848100E+00_rkx, 0.848500E+00_rkx, &
    0.849000E+00_rkx, 0.851800E+00_rkx, 0.853800E+00_rkx, 0.848400E+00_rkx, &
    0.841800E+00_rkx, 0.840500E+00_rkx, 0.843200E+00_rkx, 0.844100E+00_rkx, &
    0.840400E+00_rkx, 0.840700E+00_rkx, 0.842300E+00_rkx, 0.893700E+00_rkx, &
    0.818400E+00_rkx, 0.818400E+00_rkx, 0.804400E+00_rkx, 0.802400E+00_rkx, &
    0.801500E+00_rkx, 0.800700E+00_rkx, 0.799800E+00_rkx, 0.799000E+00_rkx, &
    0.796200E+00_rkx, 0.786000E+00_rkx, 0.785100E+00_rkx, 0.790400E+00_rkx, &
    0.790100E+00_rkx, 0.790800E+00_rkx, 0.793400E+00_rkx, 0.796100E+00_rkx, &
    0.797800E+00_rkx, 0.802600E+00_rkx, 0.806200E+00_rkx, 0.739000E+00_rkx, &
    0.739000E+00_rkx, 0.835900E+00_rkx, 0.850000E+00_rkx, 0.854000E+00_rkx, &
    0.856700E+00_rkx, 0.858400E+00_rkx, 0.859000E+00_rkx, 0.855900E+00_rkx, &
    0.855700E+00_rkx, 0.854800E+00_rkx, 0.847400E+00_rkx, 0.846700E+00_rkx, &
    0.849900E+00_rkx, 0.850600E+00_rkx, 0.847000E+00_rkx, 0.847200E+00_rkx, &
    0.848800E+00_rkx, 0.898400E+00_rkx, 0.827100E+00_rkx, 0.827100E+00_rkx, &
    0.817100E+00_rkx, 0.815500E+00_rkx, 0.814700E+00_rkx, 0.813800E+00_rkx, &
    0.812900E+00_rkx, 0.812000E+00_rkx, 0.808700E+00_rkx, 0.795500E+00_rkx, &
    0.790600E+00_rkx, 0.795900E+00_rkx, 0.795900E+00_rkx, 0.797400E+00_rkx, &
    0.800800E+00_rkx, 0.803500E+00_rkx, 0.805600E+00_rkx, 0.811900E+00_rkx, &
    0.826800E+00_rkx, 0.760400E+00_rkx, 0.760400E+00_rkx, 0.844500E+00_rkx, &
    0.850400E+00_rkx, 0.852300E+00_rkx, 0.853700E+00_rkx, 0.854600E+00_rkx, &
    0.855000E+00_rkx, 0.854200E+00_rkx, 0.861000E+00_rkx, 0.856800E+00_rkx, &
    0.855600E+00_rkx, 0.855000E+00_rkx, 0.858300E+00_rkx, 0.860200E+00_rkx, &
    0.857300E+00_rkx, 0.857700E+00_rkx, 0.859600E+00_rkx, 0.905800E+00_rkx, &
    0.841600E+00_rkx, 0.841600E+00_rkx, 0.829600E+00_rkx, 0.826100E+00_rkx, &
    0.824700E+00_rkx, 0.823600E+00_rkx, 0.822700E+00_rkx, 0.822000E+00_rkx, &
    0.820500E+00_rkx, 0.804400E+00_rkx, 0.798400E+00_rkx, 0.799500E+00_rkx, &
    0.799800E+00_rkx, 0.801700E+00_rkx, 0.805800E+00_rkx, 0.808200E+00_rkx, &
    0.810700E+00_rkx, 0.818300E+00_rkx, 0.844700E+00_rkx, 0.780000E+00_rkx, &
    0.780000E+00_rkx, 0.839700E+00_rkx, 0.842900E+00_rkx, 0.844300E+00_rkx, &
    0.845600E+00_rkx, 0.846800E+00_rkx, 0.848000E+00_rkx, 0.851900E+00_rkx, &
    0.864600E+00_rkx, 0.863300E+00_rkx, 0.862300E+00_rkx, 0.862200E+00_rkx, &
    0.866200E+00_rkx, 0.868300E+00_rkx, 0.866000E+00_rkx, 0.866500E+00_rkx, &
    0.868800E+00_rkx, 0.913700E+00_rkx, 0.858300E+00_rkx, 0.858300E+00_rkx, &
    0.837900E+00_rkx, 0.838000E+00_rkx, 0.837800E+00_rkx, 0.837400E+00_rkx, &
    0.836800E+00_rkx, 0.836000E+00_rkx, 0.832400E+00_rkx, 0.817400E+00_rkx, &
    0.808000E+00_rkx, 0.804500E+00_rkx, 0.804000E+00_rkx, 0.805600E+00_rkx, &
    0.809700E+00_rkx, 0.811200E+00_rkx, 0.813900E+00_rkx, 0.822600E+00_rkx, &
    0.863700E+00_rkx, 0.798200E+00_rkx, 0.798200E+00_rkx, 0.839200E+00_rkx, &
    0.837900E+00_rkx, 0.838100E+00_rkx, 0.838800E+00_rkx, 0.840200E+00_rkx, &
    0.842100E+00_rkx, 0.850900E+00_rkx, 0.864000E+00_rkx, 0.870300E+00_rkx, &
    0.867000E+00_rkx, 0.868500E+00_rkx, 0.873000E+00_rkx, 0.876400E+00_rkx, &
    0.875400E+00_rkx, 0.876500E+00_rkx, 0.880100E+00_rkx, 0.924100E+00_rkx, &
    0.882100E+00_rkx, 0.882100E+00_rkx, 0.844600E+00_rkx, 0.845000E+00_rkx, &
    0.845000E+00_rkx, 0.844800E+00_rkx, 0.844500E+00_rkx, 0.844000E+00_rkx, &
    0.841600E+00_rkx, 0.828000E+00_rkx, 0.817600E+00_rkx, 0.809000E+00_rkx, &
    0.808000E+00_rkx, 0.809200E+00_rkx, 0.812700E+00_rkx, 0.812900E+00_rkx, &
    0.815500E+00_rkx, 0.824200E+00_rkx, 0.874500E+00_rkx, 0.808300E+00_rkx, &
    0.808300E+00_rkx, 0.807500E+00_rkx, 0.820000E+00_rkx, 0.824600E+00_rkx, &
    0.828500E+00_rkx, 0.832000E+00_rkx, 0.835000E+00_rkx, 0.843300E+00_rkx, &
    0.862700E+00_rkx, 0.870200E+00_rkx, 0.870300E+00_rkx, 0.872400E+00_rkx, &
    0.877300E+00_rkx, 0.881500E+00_rkx, 0.881400E+00_rkx, 0.882800E+00_rkx, &
    0.887300E+00_rkx, 0.930800E+00_rkx, 0.900600E+00_rkx, 0.900600E+00_rkx ], &
     [nwav,2,nih])

  real(rkx) , dimension(nbndsw,4) , parameter :: ksdust_rrtm = reshape(&
  [ 0.324733E-01_rkx, 0.733946E-01_rkx, 0.142493E+00_rkx, 0.223231E+00_rkx, &
    0.354726E+00_rkx, 0.643902E+00_rkx, 0.962607E+00_rkx, 0.179510E+01_rkx, &
    0.307554E+01_rkx, 0.352915E+01_rkx, 0.275913E+01_rkx, 0.187476E+01_rkx, &
    0.212181E+01_rkx, 0.392232E-02_rkx, 0.345723E+00_rkx, 0.627518E+00_rkx, &
    0.824609E+00_rkx, 0.107062E+01_rkx, 0.119388E+01_rkx, 0.135678E+01_rkx, &
    0.137925E+01_rkx, 0.968465E+00_rkx, 0.658609E+00_rkx, 0.863619E+00_rkx, &
    0.718698E+00_rkx, 0.744859E+00_rkx, 0.726664E+00_rkx, 0.474096E-01_rkx, &
    0.611556E+00_rkx, 0.635682E+00_rkx, 0.549972E+00_rkx, 0.439040E+00_rkx, &
    0.320780E+00_rkx, 0.317542E+00_rkx, 0.408074E+00_rkx, 0.385324E+00_rkx, &
    0.378327E+00_rkx, 0.339740E+00_rkx, 0.341137E+00_rkx, 0.334836E+00_rkx, &
    0.329921E+00_rkx, 0.185288E+00_rkx, 0.131353E+00_rkx, 0.185898E+00_rkx, &
    0.170663E+00_rkx, 0.135025E+00_rkx, 0.157492E+00_rkx, 0.151201E+00_rkx, &
    0.154931E+00_rkx, 0.150441E+00_rkx, 0.145725E+00_rkx, 0.143190E+00_rkx, &
    0.141600E+00_rkx, 0.140409E+00_rkx, 0.139311E+00_rkx, 0.221978E+00_rkx ], &
     [nbndsw,4])

  real(rkx) , dimension(nbndsw,4) , parameter :: wsdust_rrtm = reshape(&
  [ 0.907714E+00_rkx, 0.941418E+00_rkx, 0.959479E+00_rkx, 0.968132E+00_rkx, &
    0.974402E+00_rkx, 0.979381E+00_rkx, 0.981505E+00_rkx, 0.985477E+00_rkx, &
    0.987068E+00_rkx, 0.966195E+00_rkx, 0.868209E+00_rkx, 0.669915E+00_rkx, &
    0.621792E+00_rkx, 0.544753E+00_rkx, 0.984963E+00_rkx, 0.987864E+00_rkx, &
    0.988337E+00_rkx, 0.988170E+00_rkx, 0.988054E+00_rkx, 0.985916E+00_rkx, &
    0.983429E+00_rkx, 0.966435E+00_rkx, 0.933441E+00_rkx, 0.877245E+00_rkx, &
    0.696429E+00_rkx, 0.591606E+00_rkx, 0.543167E+00_rkx, 0.935908E+00_rkx, &
    0.989606E+00_rkx, 0.986320E+00_rkx, 0.979395E+00_rkx, 0.968969E+00_rkx, &
    0.951572E+00_rkx, 0.937219E+00_rkx, 0.947774E+00_rkx, 0.926920E+00_rkx, &
    0.909271E+00_rkx, 0.787770E+00_rkx, 0.611392E+00_rkx, 0.543835E+00_rkx, &
    0.536133E+00_rkx, 0.985617E+00_rkx, 0.948398E+00_rkx, 0.951511E+00_rkx, &
    0.936367E+00_rkx, 0.912397E+00_rkx, 0.917815E+00_rkx, 0.899879E+00_rkx, &
    0.891498E+00_rkx, 0.866052E+00_rkx, 0.824030E+00_rkx, 0.679679E+00_rkx, &
    0.551325E+00_rkx, 0.542113E+00_rkx, 0.544148E+00_rkx, 0.987402E+00_rkx ], &
     [nbndsw,4])

  real(rkx) , dimension(nbndsw,4) , parameter :: gsdust_rrtm = reshape(&
  [ 0.714045E-01_rkx, 0.109877E+00_rkx, 0.158598E+00_rkx, 0.207594E+00_rkx, &
    0.284977E+00_rkx, 0.453583E+00_rkx, 0.584279E+00_rkx, 0.636315E+00_rkx, &
    0.723353E+00_rkx, 0.737234E+00_rkx, 0.684414E+00_rkx, 0.627575E+00_rkx, &
    0.817597E+00_rkx, 0.184466E-01_rkx, 0.557701E+00_rkx, 0.618631E+00_rkx, &
    0.659187E+00_rkx, 0.723095E+00_rkx, 0.721698E+00_rkx, 0.736860E+00_rkx, &
    0.725706E+00_rkx, 0.612740E+00_rkx, 0.547276E+00_rkx, 0.779661E+00_rkx, &
    0.821398E+00_rkx, 0.895405E+00_rkx, 0.930697E+00_rkx, 0.134666E+00_rkx, &
    0.729551E+00_rkx, 0.727320E+00_rkx, 0.681214E+00_rkx, 0.600556E+00_rkx, &
    0.501661E+00_rkx, 0.556412E+00_rkx, 0.727370E+00_rkx, 0.755553E+00_rkx, &
    0.779724E+00_rkx, 0.817743E+00_rkx, 0.904363E+00_rkx, 0.939012E+00_rkx, &
    0.946158E+00_rkx, 0.455520E+00_rkx, 0.525760E+00_rkx, 0.739899E+00_rkx, &
    0.768289E+00_rkx, 0.716044E+00_rkx, 0.763721E+00_rkx, 0.765208E+00_rkx, &
    0.792828E+00_rkx, 0.809833E+00_rkx, 0.833725E+00_rkx, 0.888829E+00_rkx, &
    0.941787E+00_rkx, 0.947955E+00_rkx, 0.948354E+00_rkx, 0.680178E+00_rkx ], &
     [nbndsw,4])

  real(rkx) , dimension(nbndsw,12) , parameter :: ksdust12_rrtm = reshape(&
  [ 0.272948E-02_rkx, 0.396960E-02_rkx, 0.558669E-02_rkx, 0.726846E-02_rkx, &
    0.100291E-01_rkx, 0.169911E-01_rkx, 0.249299E-01_rkx, 0.609426E-01_rkx, &
    0.201171E+00_rkx, 0.622058E+00_rkx, 0.177007E+01_rkx, 0.385513E+01_rkx, &
    0.761601E+01_rkx, 0.997750E-03_rkx, 0.105097E-01_rkx, 0.223419E-01_rkx, &
    0.428785E-01_rkx, 0.686183E-01_rkx, 0.117021E+00_rkx, 0.247432E+00_rkx, &
    0.393107E+00_rkx, 0.891504E+00_rkx, 0.239732E+01_rkx, 0.383522E+01_rkx, &
    0.517408E+01_rkx, 0.504019E+01_rkx, 0.365543E+01_rkx, 0.176373E-02_rkx, &
    0.128690E+00_rkx, 0.257957E+00_rkx, 0.433757E+00_rkx, 0.646740E+00_rkx, &
    0.961487E+00_rkx, 0.132737E+01_rkx, 0.174189E+01_rkx, 0.203876E+01_rkx, &
    0.186163E+01_rkx, 0.117175E+01_rkx, 0.125537E+01_rkx, 0.130027E+01_rkx, &
    0.119341E+01_rkx, 0.143158E-01_rkx, 0.458140E+00_rkx, 0.675648E+00_rkx, &
    0.939034E+00_rkx, 0.102583E+01_rkx, 0.116843E+01_rkx, 0.113745E+01_rkx, &
    0.105356E+01_rkx, 0.682357E+00_rkx, 0.722278E+00_rkx, 0.687013E+00_rkx, &
    0.672982E+00_rkx, 0.630313E+00_rkx, 0.621575E+00_rkx, 0.647129E-01_rkx, &
    0.625384E+00_rkx, 0.747434E+00_rkx, 0.762785E+00_rkx, 0.683835E+00_rkx, &
    0.561863E+00_rkx, 0.372238E+00_rkx, 0.333502E+00_rkx, 0.462954E+00_rkx, &
    0.393826E+00_rkx, 0.427299E+00_rkx, 0.410053E+00_rkx, 0.401259E+00_rkx, &
    0.394794E+00_rkx, 0.147315E+00_rkx, 0.567730E+00_rkx, 0.512739E+00_rkx, &
    0.376847E+00_rkx, 0.285331E+00_rkx, 0.253167E+00_rkx, 0.347236E+00_rkx, &
    0.396214E+00_rkx, 0.320372E+00_rkx, 0.296730E+00_rkx, 0.306369E+00_rkx, &
    0.292894E+00_rkx, 0.291139E+00_rkx, 0.287308E+00_rkx, 0.211338E+00_rkx, &
    0.420692E+00_rkx, 0.288169E+00_rkx, 0.209533E+00_rkx, 0.219521E+00_rkx, &
    0.286868E+00_rkx, 0.288531E+00_rkx, 0.223144E+00_rkx, 0.247584E+00_rkx, &
    0.249322E+00_rkx, 0.242542E+00_rkx, 0.236534E+00_rkx, 0.233807E+00_rkx, &
    0.231411E+00_rkx, 0.242352E+00_rkx, 0.220517E+00_rkx, 0.164365E+00_rkx, &
    0.215889E+00_rkx, 0.252894E+00_rkx, 0.216731E+00_rkx, 0.189873E+00_rkx, &
    0.226244E+00_rkx, 0.195196E+00_rkx, 0.192359E+00_rkx, 0.189225E+00_rkx, &
    0.185151E+00_rkx, 0.183256E+00_rkx, 0.181592E+00_rkx, 0.251536E+00_rkx, &
    0.146855E+00_rkx, 0.149056E+00_rkx, 0.116729E+00_rkx, 0.140686E+00_rkx, &
    0.133201E+00_rkx, 0.129320E+00_rkx, 0.118460E+00_rkx, 0.125126E+00_rkx, &
    0.122366E+00_rkx, 0.121452E+00_rkx, 0.119775E+00_rkx, 0.118858E+00_rkx, &
    0.118014E+00_rkx, 0.194438E+00_rkx, 0.805498E-01_rkx, 0.784490E-01_rkx, &
    0.828706E-01_rkx, 0.736198E-01_rkx, 0.783249E-01_rkx, 0.761732E-01_rkx, &
    0.790690E-01_rkx, 0.748942E-01_rkx, 0.737752E-01_rkx, 0.732786E-01_rkx, &
    0.726544E-01_rkx, 0.722256E-01_rkx, 0.718499E-01_rkx, 0.100233E+00_rkx, &
    0.474997E-01_rkx, 0.471086E-01_rkx, 0.462397E-01_rkx, 0.462776E-01_rkx, &
    0.455693E-01_rkx, 0.452760E-01_rkx, 0.442922E-01_rkx, 0.445527E-01_rkx, &
    0.440550E-01_rkx, 0.437889E-01_rkx, 0.435189E-01_rkx, 0.433345E-01_rkx, &
    0.431715E-01_rkx, 0.490044E-01_rkx, 0.300359E-01_rkx, 0.293560E-01_rkx, &
    0.292904E-01_rkx, 0.290049E-01_rkx, 0.289160E-01_rkx, 0.288026E-01_rkx, &
    0.285342E-01_rkx, 0.284688E-01_rkx, 0.282308E-01_rkx, 0.280936E-01_rkx, &
    0.279676E-01_rkx, 0.100000E-19_rkx, 0.100000E-19_rkx, 0.317719E-01_rkx ], &
     [nbndsw,12])

  real(rkx) , dimension(nbndsw,12) , parameter :: wsdust12_rrtm = reshape(&
  [ 0.109340E+00_rkx, 0.178006E+00_rkx, 0.259037E+00_rkx, 0.330266E+00_rkx, &
    0.419588E+00_rkx, 0.556481E+00_rkx, 0.648767E+00_rkx, 0.774388E+00_rkx, &
    0.902079E+00_rkx, 0.901246E+00_rkx, 0.869690E+00_rkx, 0.842400E+00_rkx, &
    0.848520E+00_rkx, 0.209253E-01_rkx, 0.744437E+00_rkx, 0.832865E+00_rkx, &
    0.886310E+00_rkx, 0.914146E+00_rkx, 0.936441E+00_rkx, 0.958262E+00_rkx, &
    0.968363E+00_rkx, 0.976674E+00_rkx, 0.984543E+00_rkx, 0.972366E+00_rkx, &
    0.932847E+00_rkx, 0.858172E+00_rkx, 0.709509E+00_rkx, 0.291017E+00_rkx, &
    0.971476E+00_rkx, 0.978937E+00_rkx, 0.981759E+00_rkx, 0.983361E+00_rkx, &
    0.986080E+00_rkx, 0.987344E+00_rkx, 0.987373E+00_rkx, 0.986210E+00_rkx, &
    0.977098E+00_rkx, 0.891438E+00_rkx, 0.763443E+00_rkx, 0.657181E+00_rkx, &
    0.570206E+00_rkx, 0.811386E+00_rkx, 0.987375E+00_rkx, 0.988943E+00_rkx, &
    0.988627E+00_rkx, 0.988473E+00_rkx, 0.987200E+00_rkx, 0.983070E+00_rkx, &
    0.976232E+00_rkx, 0.953073E+00_rkx, 0.941579E+00_rkx, 0.850295E+00_rkx, &
    0.701149E+00_rkx, 0.572085E+00_rkx, 0.535622E+00_rkx, 0.953685E+00_rkx, &
    0.990051E+00_rkx, 0.988544E+00_rkx, 0.985287E+00_rkx, 0.980679E+00_rkx, &
    0.971644E+00_rkx, 0.946494E+00_rkx, 0.937328E+00_rkx, 0.936816E+00_rkx, &
    0.907323E+00_rkx, 0.812766E+00_rkx, 0.631386E+00_rkx, 0.548832E+00_rkx, &
    0.534532E+00_rkx, 0.981216E+00_rkx, 0.988374E+00_rkx, 0.982510E+00_rkx, &
    0.969910E+00_rkx, 0.952761E+00_rkx, 0.938151E+00_rkx, 0.944305E+00_rkx, &
    0.942950E+00_rkx, 0.918098E+00_rkx, 0.890279E+00_rkx, 0.776625E+00_rkx, &
    0.592917E+00_rkx, 0.541316E+00_rkx, 0.537518E+00_rkx, 0.987701E+00_rkx, &
    0.983892E+00_rkx, 0.967867E+00_rkx, 0.945027E+00_rkx, 0.939744E+00_rkx, &
    0.946501E+00_rkx, 0.934754E+00_rkx, 0.907787E+00_rkx, 0.904123E+00_rkx, &
    0.876244E+00_rkx, 0.748010E+00_rkx, 0.576659E+00_rkx, 0.539378E+00_rkx, &
    0.539774E+00_rkx, 0.989537E+00_rkx, 0.967989E+00_rkx, 0.945076E+00_rkx, &
    0.947820E+00_rkx, 0.949006E+00_rkx, 0.932188E+00_rkx, 0.911477E+00_rkx, &
    0.920349E+00_rkx, 0.887840E+00_rkx, 0.853588E+00_rkx, 0.716036E+00_rkx, &
    0.561551E+00_rkx, 0.540216E+00_rkx, 0.542083E+00_rkx, 0.990139E+00_rkx, &
    0.954207E+00_rkx, 0.941136E+00_rkx, 0.914237E+00_rkx, 0.922952E+00_rkx, &
    0.908616E+00_rkx, 0.889347E+00_rkx, 0.865349E+00_rkx, 0.848355E+00_rkx, &
    0.802472E+00_rkx, 0.659625E+00_rkx, 0.547911E+00_rkx, 0.543356E+00_rkx, &
    0.545190E+00_rkx, 0.983821E+00_rkx, 0.928060E+00_rkx, 0.910361E+00_rkx, &
    0.899370E+00_rkx, 0.873460E+00_rkx, 0.867375E+00_rkx, 0.840595E+00_rkx, &
    0.830818E+00_rkx, 0.788857E+00_rkx, 0.732884E+00_rkx, 0.605375E+00_rkx, &
    0.545667E+00_rkx, 0.546259E+00_rkx, 0.547380E+00_rkx, 0.969371E+00_rkx, &
    0.896477E+00_rkx, 0.872626E+00_rkx, 0.847013E+00_rkx, 0.830465E+00_rkx, &
    0.808301E+00_rkx, 0.777415E+00_rkx, 0.750775E+00_rkx, 0.715187E+00_rkx, &
    0.657070E+00_rkx, 0.569713E+00_rkx, 0.547216E+00_rkx, 0.547985E+00_rkx, &
    0.548590E+00_rkx, 0.948524E+00_rkx, 0.858330E+00_rkx, 0.823501E+00_rkx, &
    0.793639E+00_rkx, 0.770237E+00_rkx, 0.746269E+00_rkx, 0.711804E+00_rkx, &
    0.686180E+00_rkx, 0.649858E+00_rkx, 0.601391E+00_rkx, 0.554799E+00_rkx, &
    0.548324E+00_rkx, 0.500000E+00_rkx, 0.500000E+00_rkx, 0.928404E+00_rkx ], &
     [nbndsw,12])

  real(rkx) , dimension(nbndsw,12) , parameter :: gsdust12_rrtm = reshape(&
  [ 0.346379E-02_rkx, 0.532454E-02_rkx, 0.760988E-02_rkx, 0.978847E-02_rkx, &
    0.129334E-01_rkx, 0.192973E-01_rkx, 0.252011E-01_rkx, 0.419456E-01_rkx, &
    0.828151E-01_rkx, 0.147716E+00_rkx, 0.284132E+00_rkx, 0.504929E+00_rkx, &
    0.632938E+00_rkx, 0.887657E-03_rkx, 0.299321E-01_rkx, 0.458652E-01_rkx, &
    0.654008E-01_rkx, 0.840745E-01_rkx, 0.111347E+00_rkx, 0.168701E+00_rkx, &
    0.225066E+00_rkx, 0.402610E+00_rkx, 0.619089E+00_rkx, 0.685494E+00_rkx, &
    0.744738E+00_rkx, 0.752464E+00_rkx, 0.696838E+00_rkx, 0.771129E-02_rkx, &
    0.202337E+00_rkx, 0.330338E+00_rkx, 0.495852E+00_rkx, 0.601693E+00_rkx, &
    0.619352E+00_rkx, 0.658747E+00_rkx, 0.726645E+00_rkx, 0.729282E+00_rkx, &
    0.685932E+00_rkx, 0.544204E+00_rkx, 0.733322E+00_rkx, 0.865095E+00_rkx, &
    0.893021E+00_rkx, 0.506849E-01_rkx, 0.616252E+00_rkx, 0.641233E+00_rkx, &
    0.722164E+00_rkx, 0.718575E+00_rkx, 0.738646E+00_rkx, 0.718152E+00_rkx, &
    0.679960E+00_rkx, 0.552694E+00_rkx, 0.702482E+00_rkx, 0.779828E+00_rkx, &
    0.850731E+00_rkx, 0.905884E+00_rkx, 0.934684E+00_rkx, 0.180875E+00_rkx, &
    0.721069E+00_rkx, 0.735060E+00_rkx, 0.724757E+00_rkx, 0.695022E+00_rkx, &
    0.630267E+00_rkx, 0.490052E+00_rkx, 0.529683E+00_rkx, 0.713001E+00_rkx, &
    0.739448E+00_rkx, 0.807462E+00_rkx, 0.892127E+00_rkx, 0.933183E+00_rkx, &
    0.944583E+00_rkx, 0.366911E+00_rkx, 0.735882E+00_rkx, 0.699022E+00_rkx, &
    0.603182E+00_rkx, 0.489171E+00_rkx, 0.494111E+00_rkx, 0.702517E+00_rkx, &
    0.776974E+00_rkx, 0.758849E+00_rkx, 0.766561E+00_rkx, 0.840093E+00_rkx, &
    0.912016E+00_rkx, 0.942013E+00_rkx, 0.946880E+00_rkx, 0.523451E+00_rkx, &
    0.700282E+00_rkx, 0.571750E+00_rkx, 0.473265E+00_rkx, 0.558598E+00_rkx, &
    0.714231E+00_rkx, 0.771167E+00_rkx, 0.717543E+00_rkx, 0.760124E+00_rkx, &
    0.802612E+00_rkx, 0.853897E+00_rkx, 0.924673E+00_rkx, 0.945047E+00_rkx, &
    0.947578E+00_rkx, 0.619055E+00_rkx, 0.557231E+00_rkx, 0.498460E+00_rkx, &
    0.686063E+00_rkx, 0.772881E+00_rkx, 0.765081E+00_rkx, 0.740547E+00_rkx, &
    0.793729E+00_rkx, 0.782627E+00_rkx, 0.815096E+00_rkx, 0.872429E+00_rkx, &
    0.934774E+00_rkx, 0.946979E+00_rkx, 0.948038E+00_rkx, 0.674429E+00_rkx, &
    0.699890E+00_rkx, 0.767045E+00_rkx, 0.721469E+00_rkx, 0.775748E+00_rkx, &
    0.769894E+00_rkx, 0.791408E+00_rkx, 0.784479E+00_rkx, 0.817454E+00_rkx, &
    0.844005E+00_rkx, 0.899659E+00_rkx, 0.944573E+00_rkx, 0.948268E+00_rkx, &
    0.948490E+00_rkx, 0.673365E+00_rkx, 0.750332E+00_rkx, 0.760635E+00_rkx, &
    0.801505E+00_rkx, 0.787758E+00_rkx, 0.811382E+00_rkx, 0.824420E+00_rkx, &
    0.842833E+00_rkx, 0.851547E+00_rkx, 0.876250E+00_rkx, 0.923888E+00_rkx, &
    0.948195E+00_rkx, 0.948695E+00_rkx, 0.948714E+00_rkx, 0.639275E+00_rkx, &
    0.794881E+00_rkx, 0.810746E+00_rkx, 0.824644E+00_rkx, 0.835258E+00_rkx, &
    0.844376E+00_rkx, 0.859087E+00_rkx, 0.866250E+00_rkx, 0.884484E+00_rkx, &
    0.906868E+00_rkx, 0.939399E+00_rkx, 0.948839E+00_rkx, 0.948832E+00_rkx, &
    0.948767E+00_rkx, 0.680069E+00_rkx, 0.823712E+00_rkx, 0.837959E+00_rkx, &
    0.853110E+00_rkx, 0.861751E+00_rkx, 0.872182E+00_rkx, 0.886589E+00_rkx, &
    0.895115E+00_rkx, 0.909941E+00_rkx, 0.928012E+00_rkx, 0.945923E+00_rkx, &
    0.948915E+00_rkx, 0.990000E+00_rkx, 0.990000E+00_rkx, 0.761202E+00_rkx ], &
     [nbndsw,12])

  real(rkx) , dimension(nbndlw,4) , parameter :: ksdust_rrtm_lw = reshape(&
  [ 0.903338E-07_rkx, 0.529696E-04_rkx, 0.103134E-03_rkx, 0.114509E-03_rkx, &
    0.242895E-03_rkx, 0.578706E-03_rkx, 0.118229E-02_rkx, 0.156680E-02_rkx, &
    0.711396E-03_rkx, 0.142405E-02_rkx, 0.254564E-02_rkx, 0.656414E-02_rkx, &
    0.965619E-02_rkx, 0.118215E-01_rkx, 0.147898E-01_rkx, 0.263245E-01_rkx, &
    0.166503E-05_rkx, 0.102266E-02_rkx, 0.198371E-02_rkx, 0.217008E-02_rkx, &
    0.457206E-02_rkx, 0.113272E-01_rkx, 0.231581E-01_rkx, 0.295729E-01_rkx, &
    0.117724E-01_rkx, 0.246294E-01_rkx, 0.413630E-01_rkx, 0.100878E+00_rkx, &
    0.144428E+00_rkx, 0.163034E+00_rkx, 0.185000E+00_rkx, 0.262229E+00_rkx, &
    0.139793E-04_rkx, 0.102906E-01_rkx, 0.188747E-01_rkx, 0.189807E-01_rkx, &
    0.387656E-01_rkx, 0.888989E-01_rkx, 0.135854E+00_rkx, 0.132997E+00_rkx, &
    0.478602E-01_rkx, 0.997541E-01_rkx, 0.142638E+00_rkx, 0.292038E+00_rkx, &
    0.371668E+00_rkx, 0.392363E+00_rkx, 0.400920E+00_rkx, 0.406941E+00_rkx, &
    0.123597E-03_rkx, 0.909172E-01_rkx, 0.867894E-01_rkx, 0.734485E-01_rkx, &
    0.120615E+00_rkx, 0.142348E+00_rkx, 0.996127E-01_rkx, 0.900509E-01_rkx, &
    0.683209E-01_rkx, 0.119819E+00_rkx, 0.125996E+00_rkx, 0.122892E+00_rkx, &
    0.106664E+00_rkx, 0.976398E-01_rkx, 0.913925E-01_rkx, 0.763635E-01_rkx ],&
    [nbndlw,4])

  real(rkx) , dimension(nbndlw,12) , parameter :: ksdust12_rrtm_lw = reshape(&
  [ 0.943320E-09_rkx, 0.549405E-06_rkx, 0.106979E-05_rkx, 0.119006E-05_rkx, &
    0.252349E-05_rkx, 0.596784E-05_rkx, 0.121514E-04_rkx, 0.161361E-04_rkx, &
    0.748894E-05_rkx, 0.148517E-04_rkx, 0.266658E-04_rkx, 0.692718E-04_rkx, &
    0.992866E-04_rkx, 0.121841E-03_rkx, 0.152949E-03_rkx, 0.274053E-03_rkx, &
    0.242747E-07_rkx, 0.141745E-04_rkx, 0.276003E-04_rkx, 0.306811E-04_rkx, &
    0.650768E-04_rkx, 0.154341E-03_rkx, 0.314708E-03_rkx, 0.417667E-03_rkx, &
    0.192223E-03_rkx, 0.382566E-03_rkx, 0.686078E-03_rkx, 0.177937E-02_rkx, &
    0.257410E-02_rkx, 0.315727E-02_rkx, 0.396055E-02_rkx, 0.709637E-02_rkx, &
    0.412358E-06_rkx, 0.244980E-03_rkx, 0.476715E-03_rkx, 0.527223E-03_rkx, &
    0.111702E-02_rkx, 0.269697E-02_rkx, 0.553147E-02_rkx, 0.727898E-02_rkx, &
    0.317285E-02_rkx, 0.645799E-02_rkx, 0.113876E-01_rkx, 0.290183E-01_rkx, &
    0.438245E-01_rkx, 0.528732E-01_rkx, 0.648010E-01_rkx, 0.107272E+00_rkx, &
    0.249960E-05_rkx, 0.156266E-02_rkx, 0.302346E-02_rkx, 0.328672E-02_rkx, &
    0.690405E-02_rkx, 0.172493E-01_rkx, 0.349294E-01_rkx, 0.436492E-01_rkx, &
    0.167118E-01_rkx, 0.351449E-01_rkx, 0.571044E-01_rkx, 0.133955E+00_rkx, &
    0.186049E+00_rkx, 0.207183E+00_rkx, 0.234767E+00_rkx, 0.326611E+00_rkx, &
    0.867040E-05_rkx, 0.598433E-02_rkx, 0.112907E-01_rkx, 0.117650E-01_rkx, &
    0.243645E-01_rkx, 0.596038E-01_rkx, 0.104488E+00_rkx, 0.111867E+00_rkx, &
    0.390610E-01_rkx, 0.800462E-01_rkx, 0.116591E+00_rkx, 0.270079E+00_rkx, &
    0.342407E+00_rkx, 0.353059E+00_rkx, 0.373278E+00_rkx, 0.423478E+00_rkx, &
    0.200439E-04_rkx, 0.156573E-01_rkx, 0.276508E-01_rkx, 0.267074E-01_rkx, &
    0.534150E-01_rkx, 0.113530E+00_rkx, 0.152373E+00_rkx, 0.139808E+00_rkx, &
    0.542006E-01_rkx, 0.115857E+00_rkx, 0.156714E+00_rkx, 0.310977E+00_rkx, &
    0.374026E+00_rkx, 0.373059E+00_rkx, 0.382935E+00_rkx, 0.371831E+00_rkx, &
    0.348134E-04_rkx, 0.299434E-01_rkx, 0.474603E-01_rkx, 0.419379E-01_rkx, &
    0.802288E-01_rkx, 0.149212E+00_rkx, 0.154806E+00_rkx, 0.135289E+00_rkx, &
    0.633148E-01_rkx, 0.131083E+00_rkx, 0.168116E+00_rkx, 0.291058E+00_rkx, &
    0.337175E+00_rkx, 0.326846E+00_rkx, 0.320989E+00_rkx, 0.275473E+00_rkx, &
    0.639697E-04_rkx, 0.576057E-01_rkx, 0.734485E-01_rkx, 0.592327E-01_rkx, &
    0.112439E+00_rkx, 0.161819E+00_rkx, 0.136628E+00_rkx, 0.118753E+00_rkx, &
    0.691587E-01_rkx, 0.135583E+00_rkx, 0.159980E+00_rkx, 0.225479E+00_rkx, &
    0.228502E+00_rkx, 0.215787E+00_rkx, 0.207929E+00_rkx, 0.150134E+00_rkx, &
    0.180706E-03_rkx, 0.957461E-01_rkx, 0.813079E-01_rkx, 0.756272E-01_rkx, &
    0.117072E+00_rkx, 0.113926E+00_rkx, 0.782340E-01_rkx, 0.734595E-01_rkx, &
    0.640591E-01_rkx, 0.100870E+00_rkx, 0.971746E-01_rkx, 0.758939E-01_rkx, &
    0.660027E-01_rkx, 0.647597E-01_rkx, 0.650166E-01_rkx, 0.701898E-01_rkx, &
    0.362353E-03_rkx, 0.527134E-01_rkx, 0.500051E-01_rkx, 0.556937E-01_rkx, &
    0.553413E-01_rkx, 0.420387E-01_rkx, 0.433355E-01_rkx, 0.431202E-01_rkx, &
    0.416369E-01_rkx, 0.419620E-01_rkx, 0.396648E-01_rkx, 0.471680E-01_rkx, &
    0.475567E-01_rkx, 0.454902E-01_rkx, 0.432298E-01_rkx, 0.408139E-01_rkx, &
    0.648472E-03_rkx, 0.281645E-01_rkx, 0.266968E-01_rkx, 0.253841E-01_rkx, &
    0.255754E-01_rkx, 0.270152E-01_rkx, 0.264115E-01_rkx, 0.261331E-01_rkx, &
    0.236985E-01_rkx, 0.246062E-01_rkx, 0.244858E-01_rkx, 0.251570E-01_rkx, &
    0.245834E-01_rkx, 0.244727E-01_rkx, 0.246475E-01_rkx, 0.246670E-01_rkx, &
    0.100623E-02_rkx, 0.183919E-01_rkx, 0.172178E-01_rkx, 0.162879E-01_rkx, &
    0.165092E-01_rkx, 0.168667E-01_rkx, 0.170195E-01_rkx, 0.168656E-01_rkx, &
    0.153070E-01_rkx, 0.157179E-01_rkx, 0.157377E-01_rkx, 0.159894E-01_rkx, &
    0.160022E-01_rkx, 0.158424E-01_rkx, 0.157455E-01_rkx, 0.157207E-01_rkx ],&
    [nbndlw,12])

  data ncaec / -1 /

  contains

    subroutine allocate_mod_rad_aerosol
      implicit none
      integer(ik4) :: n , k , kk , kk1

      if ( irrtm == 1 ) then
        nband = nbndsw
      else
        nband = nspi
      end if

      nj = (jci2-jci1+1)
      npoints = nj*(ici2-ici1+1)

      if ( ichem == 1 .or. iclimaaer == 1 ) then
        call getmem3d(aermmr,1,npoints,1,kz,1,ntr,'aerosol:aermmr')
        if ( iclimaaer == 1 ) then
          call getmem4d(aerm1,jce1,jce2,ice1,ice2,1,kz,1,ntr,'aerosol:aerm1')
          call getmem4d(aerm2,jce1,jce2,ice1,ice2,1,kz,1,ntr,'aerosol:aerm2')
          if ( .not. do_parallel_netcdf_in ) then
            if ( myid == iocpu ) then
              call getmem3d(aerio,jdot1,jdot2,idot1,idot2,1,kz,'aerosol:aerio')
            end if
          end if
        end if
      end if

      call getmem1d(gsbc_hb,1,nband,'aerosol:gsbc_hb')
      call getmem1d(gsbc_hl,1,nband,'aerosol:gsbc_hl')
      call getmem1d(gsoc_hb,1,nband,'aerosol:gsoc_hb')
      call getmem1d(gsoc_hl,1,nband,'aerosol:gsoc_hl')
      call getmem1d(gssm1,1,nband,'aerosol:gssm1')
      call getmem1d(gssm2,1,nband,'aerosol:gssm2')

      call getmem1d(ksbc_hb,1,nband,'aerosol:ksbc_hb')
      call getmem1d(ksbc_hl,1,nband,'aerosol:ksbc_hl')
      call getmem1d(ksoc_hb,1,nband,'aerosol:ksoc_hb')
      call getmem1d(ksoc_hl,1,nband,'aerosol:ksoc_hl')
      call getmem1d(kssm1,1,nband,'aerosol:kssm1')
      call getmem1d(kssm2,1,nband,'aerosol:kssm2')

      call getmem1d(wsbc_hb,1,nband,'aerosol:wsbc_hb')
      call getmem1d(wsbc_hl,1,nband,'aerosol:wsbc_hl')
      call getmem1d(wsoc_hb,1,nband,'aerosol:wsoc_hl')
      call getmem1d(wsoc_hl,1,nband,'aerosol:wsoc_hl')
      call getmem1d(wssm1,1,nband,'aerosol:wssm1')
      call getmem1d(wssm2,1,nband,'aerosol:wssm2')

      call getmem2d(gsdust,1,nband,1,nbin,'aerosol:gsdust')
      call getmem2d(ksdust,1,nband,1,nbin,'aerosol:ksdust')
      call getmem2d(wsdust,1,nband,1,nbin,'aerosol:wsdust')
      ! op propert lw for rrtm
      call getmem2d(ksdust_lw,1,nbndlw,1,nbin,'aerosol:ksdust_lw')

      call getmem2d(path,1,npoints,1,kz,'aerosol:path')
      call getmem3d(ftota3d,1,npoints,0,kz,1,nband,'aerosol:ftota3d')

      ! these variables are defined on full rad grid including hat
      if ( irrtm == 1 ) then
        call getmem3d(gtota3d,1,npoints,0,kth,1,nband,'aerosol:gtota3d')
        call getmem3d(tauasc3d,1,npoints,0,kth,1,nband,'aerosol:tauasc3d')
        call getmem3d(tauxar3d,1,npoints,0,kth,1,nband,'aerosol:tauxar3d')
        call getmem3d(tauxar3d_lw,1,npoints, &
                      0,kth,1,nbndlw,'aerosol:tauxar3d_lw')
      else ! standard scheme has one extra strato level at k = 0
        call getmem3d(gtota3d,1,npoints,0,kz,1,nband,'aerosol:gtota3d')
        call getmem3d(tauasc3d,1,npoints,0,kz,1,nband,'aerosol:tauasc3d')
        call getmem3d(tauxar3d,1,npoints,0,kz,1,nband,'aerosol:tauxar3d')
      end if

      call getmem2d(ftota,1,npoints,1,nband,'aerosol:ftota')
      call getmem2d(gtota,1,npoints,1,nband,'aerosol:gtota')
      call getmem2d(tauasc,1,npoints,1,nband,'aerosol:tauasc')
      call getmem2d(tauxar,1,npoints,1,nspi,'aerosol:tauxar')
      call getmem2d(aermtot,1,npoints,1,kz,'aerosol:aermtot')
      call getmem2d(aervtot,1,npoints,1,kz,'aerosol:aervtot')
      call getmem3d(aertrlw,1,npoints,1,kzp1,1,kzp1,'aerosol:aertrlw')
      if ( ichem == 1 .or. iclimaaer == 1 ) then
        call getmem3d(fa,1,npoints,0,kz,1,ntr,'aerosol:fa')
        call getmem3d(ga,1,npoints,0,kz,1,ntr,'aerosol:ga')
        call getmem3d(tx,1,npoints,0,kz,1,ntr,'aerosol:tx')
        call getmem3d(uaer,1,npoints,0,kz,1,ntr,'aerosol:uaer')
        call getmem3d(wa,1,npoints,0,kz,1,ntr,'aerosol:wa')
        call getmem2d(faer,1,npoints,1,ntr,'aerosol:faer')
        call getmem2d(gaer,1,npoints,1,ntr,'aerosol:gaer')
        call getmem2d(tauaer,1,npoints,1,ntr,'aerosol:tauaer')
        call getmem2d(utaer,1,npoints,1,ntr,'aerosol:utaer')
        call getmem2d(waer,1,npoints,1,ntr,'aerosol:waer')
      end if

      if ( iclimaaer == 2 ) then
        if ( myid == iocpu ) then
          call getmem2d(alon,jcross1,jcross2,icross1,icross2,'aerosol:alon')
          call getmem2d(alat,jcross1,jcross2,icross1,icross2,'aerosol:alat')
        end if
        ! FAB note that prof are always determined on kth level,
        ! even with standard scheme
        call getmem4d(extprof,jci1,jci2,ici1,ici2,1,kth,1,nacwb,'rad:extprof')
        call getmem4d(asyprof,jci1,jci2,ici1,ici2,1,kth,1,nacwb,'rad:asyprof')
        call getmem4d(ssaprof,jci1,jci2,ici1,ici2,1,kth,1,nacwb,'rad:ssaprof')
      else if ( iclimaaer == 3 ) then
        macv2sp_hist = trim(inpglob)//pthsep//'CMIP6'//pthsep// &
              'AEROSOL'//pthsep//'MACv2.0-SP_v1.nc'
        select case (scenario)
          case ('ssp119', 'SSP119')
            macv2sp_scen = trim(inpglob)//pthsep//'CMIP6'//pthsep// &
              'AEROSOL'//pthsep//'MACv2.0-SP_IMAGE-SSP1-19-SPA1.nc'
          case ('ssp126', 'SSP126')
            macv2sp_scen = trim(inpglob)//pthsep//'CMIP6'//pthsep// &
              'AEROSOL'//pthsep//'MACv2.0-SP_IMAGE-SSP1-26-SPA1.nc'
          case ('ssp245', 'SSP245')
            macv2sp_scen = trim(inpglob)//pthsep//'CMIP6'//pthsep// &
              'AEROSOL'//pthsep//'MACv2.0-SP_MESSAGE-GLOBIOM-SSP2-45-SPA2.nc'
          case ('ssp370', 'SSP370')
            macv2sp_scen = trim(inpglob)//pthsep//'CMIP6'//pthsep// &
              'AEROSOL'//pthsep//'MACv2.0-SP_AIM-SSP3-Ref-SPA0.nc'
          case ('ssp434', 'SSP434')
            macv2sp_scen = trim(inpglob)//pthsep//'CMIP6'//pthsep// &
              'AEROSOL'//pthsep//'MACv2.0-SP_GCAM4-SSP4-34-SPA4.nc'
          case ('ssp460', 'SSP460')
            macv2sp_scen = trim(inpglob)//pthsep//'CMIP6'//pthsep// &
              'AEROSOL'//pthsep//'MACv2.0-SP_GCAM4-SSP4-60-SPA4.nc'
          case ('ssp534', 'SSP534')
            macv2sp_scen = trim(inpglob)//pthsep//'CMIP6'//pthsep// &
              'AEROSOL'//pthsep//'MACv2.0-SP_REMIND-MAGPIE-SSP5-34-OS.nc'
          case ('ssp585', 'SSP585')
            macv2sp_scen = trim(inpglob)//pthsep//'CMIP6'//pthsep// &
              'AEROSOL'//pthsep//'MACv2.0-SP_REMIND-MAGPIE-SSP5-Ref.nc'
          case default
            macv2sp_scen = macv2sp_hist
        end select
        call getmem1d(dnovrnr4,1,npoints,'rad:dnovrnr4')
        call getmem3d(extprofr4,1,npoints,1,kth,1,nband,'rad:extprofr4')
        call getmem3d(asyprofr4,1,npoints,1,kth,1,nband,'rad:asyprofr4')
        call getmem3d(ssaprofr4,1,npoints,1,kth,1,nband,'rad:ssaprofr4')
        call getmem1d(lambdaw,1,nband,'rad:lambdaw')
        call getmem1d(latr4,1,npoints,'aerosol:latr4')
        call getmem1d(lonr4,1,npoints,'aerosol:lonr4')
        call getmem1d(altr4,1,npoints,'aerosol:altr4')
        call getmem2d(z,1,npoints,1,kth,'aerosol:z')
        call getmem2d(dz,1,npoints,1,kth,'aerosol:dz')
        do k = 1, kth-kz
          kk = n_prehlev-k+1
          kk1 = n_hreflev-k+1
          do n = 1, npoints
            z(:,k) = stdhlevh(kk)*d_1000
            dz(:,k) = (stdhlevf(kk1)-stdhlevf(kk1-1))*d_1000
          end do
        end do
      end if

      ! initialise aerosol properties in function of the radiation
      ! and dust scheme option
      if ( irrtm == 1 ) then
        gsbc_hb = gsbc_hb_rrtm
        gsbc_hl = gsbc_hl_rrtm
        gsoc_hb = gsoc_hb_rrtm
        gsoc_hl = gsoc_hl_rrtm
        ksbc_hb = ksbc_hb_rrtm
        ksbc_hl = ksbc_hl_rrtm
        ksoc_hb = ksoc_hb_rrtm
        ksoc_hl = ksoc_hl_rrtm
        wsbc_hb = wsbc_hb_rrtm
        wsbc_hl = wsbc_hl_rrtm
        wsoc_hb = wsoc_hb_rrtm
        wsoc_hl = wsoc_hl_rrtm
        if ( ismoke > 0) then
          gssm1 = gssm1_rrtm
          kssm1 = kssm1_rrtm
          wssm1 = wssm1_rrtm
          gssm2 = gssm2_rrtm
          kssm2 = kssm2_rrtm
          wssm2 = wssm2_rrtm
        end if
        if ( nbin == 4 ) then
          gsdust =  gsdust_rrtm
          ksdust =  ksdust_rrtm
          wsdust =  wsdust_rrtm
          ksdust_lw = ksdust_rrtm_lw
        else if ( nbin == 12 ) then
          gsdust =  gsdust12_rrtm
          ksdust =  ksdust12_rrtm
          wsdust =  wsdust12_rrtm
          ksdust_lw = ksdust12_rrtm_lw
        end if
      else if ( irrtm == 0 ) then
        gsbc_hb = gsbc_hb_stand
        gsbc_hl = gsbc_hl_stand
        gsoc_hb = gsoc_hb_stand
        gsoc_hl = gsoc_hl_stand
        ksbc_hb = ksbc_hb_stand
        ksbc_hl = ksbc_hl_stand
        ksoc_hb = ksoc_hb_stand
        ksoc_hl = ksoc_hl_stand
        wsbc_hb = wsbc_hb_stand
        wsbc_hl = wsbc_hl_stand
        wsoc_hb = wsoc_hb_stand
        wsoc_hl = wsoc_hl_stand
        if ( nbin == 4 ) then
          gsdust =  gsdust_stand
          ksdust =  ksdust_stand
          wsdust =  wsdust_stand
        else if (nbin == 12) then
          gsdust =  gsdust12_stand
          ksdust =  ksdust12_stand
          wsdust =  wsdust12_stand
        end if
      end if
    end subroutine allocate_mod_rad_aerosol

    subroutine init_aerclima
      implicit none
      integer(ik4) :: itr
      type (rcm_time_and_date) :: aedate

      allocate(chtrname(ntr))
      do itr = 1 , ntr
        chtrname(itr) = aerclima_chtr(itr)
      end do

      aedate = idate1
      if ( aedate /= globidate1 ) then
        aedate = monfirst(aedate)
      end if

      call open_aerclima(aedate)

      aefreq%ival = ibdyfrq
      aefreq%iunit = uhrs
      aerfreq = real(ibdyfrq,rkx)
      d1 = idate1
      d2 = idate1
    end subroutine init_aerclima

    subroutine read_aerclima(idatex,m2r)
      implicit none
      type (rcm_time_and_date) , intent(in) :: idatex
      type(mod_2_rad) , intent(in) :: m2r
      type(rcm_time_interval) :: tdif
      real(rkx) :: w1 , w2 , step
      integer(ik4) :: i , j , k , n , ib

      if ( .not. aerclima_init ) then
        call init_aerclima( )
        aerclima_init = .true.
      end if
      if ( d1 == d2 ) then
        call doread(d1,m2r,1,aerm1)
        d2 = d1 + aefreq
        call doread(d2,m2r,2,aerm2)
        if ( myid == italk ) then
          write(stdout, *) 'Ready AEROSOL data from ', &
              tochar10(d1),' to ',tochar10(d2)
        end if
        step = d_zero
      else
        step = dtsec/3600.0_rkx
      end if

      if ( idatex >= d2 ) then
        aerm2 = aerm1
        d1 = d2
        d2 = d1 + aefreq
        if ( d2 > aetime(naetime) ) then
          call open_aerclima(d2)
        end if
        call doread(d2,m2r,2,aerm2)
        if ( myid == italk ) then
          write(stdout, *) 'Ready AEROSOL data from ', &
              tochar10(d1),' to ',tochar10(d2)
        end if
      end if

      tdif = d2 - idatex
      w1 = (real(tohours(tdif),rkx)-step)/aerfreq
      w2 = d_one - w1
      do n = 1 , ntr
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              ib = (j-jci1)+(i-ici1)*nj+1
              aermmr(ib,k,n) = max(w1*aerm1(j,i,k,n) + w2*aerm2(j,i,k,n),d_zero)
            end do
          end do
        end do
      end do
    end subroutine read_aerclima

    integer(ik4) function findrec(idate) result(irec)
      implicit none
      type (rcm_time_and_date) , intent(in) :: idate
      type(rcm_time_interval) :: tdif

      irec = -1
      if ( idate < aetime(1) ) return
      if ( idate > aetime(naetime) ) return
      tdif = idate - aetime(1)
      irec = (nint(tohours(tdif))/ibdyfrq)+1
    end function findrec

    subroutine doread(idate,m2r,step,aerm)
      implicit none
      type (rcm_time_and_date) , intent(in) :: idate
      type(mod_2_rad) , intent(in) :: m2r
      integer(ik4) , intent(in) :: step
      real(rkx) , dimension(:,:,:,:) , pointer , intent(inout) :: aerm
      real(rkx) , dimension(:,:,:) , pointer :: pnt
      integer(ik4) :: n , irec
      integer(ik4) , dimension(4) :: istart , icount

      irec = findrec(idate)
      if ( irec < 0 ) then
        write(stderr,* ) 'Searching for date '//tochar(idate)// &
          ' in AE file : NOT FOUND'
        write(stderr,* ) 'Available dates : from '//tochar10(aetime(1))// &
          ' to '//tochar10(aetime(naetime))
        call fatal(__FILE__,__LINE__,'NO AEROSOL DATA')
      end if

      if ( .not. do_parallel_netcdf_in ) then
        if ( myid /= iocpu ) then
          do n = 1 , ntr
            call assignpnt(aerm,pnt,n)
            call grid_distribute(aerio,pnt,jce1,jce2,ice1,ice2,1,kz)
          end do
          return
        else
          istart(1) = 1
          istart(2) = 1
          istart(3) = 1
          istart(4) = irec
          icount(1) = jdot2-jdot1 + 1
          icount(2) = idot2-idot1 + 1
          icount(3) = kz
          icount(4) = 1
          do n = 1 , ntr
            call assignpnt(aerm,pnt,n)
            ncstatus = nf90_get_var(ncaec,ncaevar(n),aerio,istart,icount)
            call check_ok(__FILE__,__LINE__, &
               'Error reading variable '//trim(aerclima_chtr(n))// &
               ' in AE file','AEBC FILE')
            call grid_distribute(aerio,pnt,jce1,jce2,ice1,ice2,1,kz)
          end do
        end if
      else
        istart(1) = jce1
        istart(2) = ice1
        istart(3) = 1
        istart(4) = irec
        icount(1) = (jce2-jce1)+1
        icount(2) = (ice2-ice1)+1
        icount(3) = kz
        icount(4) = 1
        do n = 1 , ntr
          call assignpnt(aerm,pnt,n)
          ncstatus = nf90_get_var(ncaec,ncaevar(n),pnt,istart,icount)
          call check_ok(__FILE__,__LINE__, &
             'Error reading variable '//trim(aerclima_chtr(n))// &
             ' in AE file','AEBC FILE')
        end do
      end if

      if ( idynamic == 2 ) then
        if ( step == 1 ) then
          call nhinterp(ice1,ice2,jce1,jce2,kz,ntr, &
                        hsigma,sigma,aerm,m2r%btv0,m2r%bps0,m2r%ps0)
        else if ( step == 2 ) then
          call nhinterp(ice1,ice2,jce1,jce2,kz,ntr, &
                        hsigma,sigma,aerm,m2r%btv1,m2r%bps1,m2r%ps0)
        end if
      end if
      if ( idynamic == 3 ) then
        if ( step == 1 ) then
          call zita_interp(jce1,jce2,ice1,ice2,kz,ntr, &
                           aerm,m2r%za,m2r%btv0,hsigma,m2r%bps0,0)
        else if ( step == 2 ) then
          call zita_interp(jce1,jce2,ice1,ice2,kz,ntr, &
                           aerm,m2r%za,m2r%btv1,hsigma,m2r%bps1,0)
        end if
      end if

    end subroutine doread

    subroutine open_aerclima(idate)
      implicit none
      type (rcm_time_and_date) , intent(in) :: idate
      integer(ik4) :: nae , idtime , n
      character(len=256) :: aefile
      real(rkx) , allocatable , dimension(:) :: xtime
      character(len=64) :: timeunits , timecal
      character(len=11) :: ctime

      if ( .not. do_parallel_netcdf_in ) then
        if ( myid /= iocpu ) then
          call bcast(naetime)
          if ( allocated(aetime) ) then
            deallocate(aetime)
          end if
          allocate(aetime(naetime))
          call bcast(aetime)
          return
        end if
      end if

      call close_aerclima

      write (ctime, '(a)') tochar10(idate)
      aefile = trim(dirglob)//pthsep//trim(domname)// &
               '_AEBC.'//trim(ctime)//'.nc'

      if ( myid == italk ) then
        write (stdout, *) 'Opening aerosol file : '//trim(aefile)
      end if

      ncstatus = nf90_open(aefile,nf90_nowrite,ncaec)
      call check_ok(__FILE__,__LINE__, &
        'Error Open AE file '//trim(aefile),'AEBC FILE')
      do nae = 1 , ntr
        ncstatus = nf90_inq_varid(ncaec,aerclima_chtr(nae),ncaevar(nae))
        call check_ok(__FILE__,__LINE__, &
          'Error searching variable '//trim(aerclima_chtr(nae))// &
          ' in AE file '//trim(aefile),'AEBC FILE')
      end do

      ncstatus = nf90_inq_dimid(ncaec,'time',idtime)
      call check_ok(__FILE__,__LINE__, &
         'Error searching dimension time in AE file '//trim(aefile),'AEBC FILE')
      ncstatus = nf90_inquire_dimension(ncaec,idtime,len=naetime)
      call check_ok(__FILE__,__LINE__, &
         'Error reading dimension time in AE file '//trim(aefile),'AEBC FILE')
      allocate(xtime(naetime))
      allocate(aetime(naetime))
      ncstatus = nf90_inq_varid(ncaec,'time',idtime)
      call check_ok(__FILE__,__LINE__, &
         'Error searching variable time in AE file '//trim(aefile),'AEBC FILE')
      ncstatus = nf90_get_att(ncaec,idtime,'units',timeunits)
      call check_ok(__FILE__,__LINE__, &
         'variable time missing attribute units '// &
         'in AE file '//trim(aefile),'AEBC FILE')
      ncstatus = nf90_get_att(ncaec,idtime,'calendar',timecal)
      call check_ok(__FILE__,__LINE__, &
         'variable time missing attribute calendar '// &
         'in AE file '//trim(aefile),'AEBC FILE')
      ncstatus = nf90_get_var(ncaec,idtime,xtime)
      call check_ok(__FILE__,__LINE__, &
         'Error reading variable time in AE file '//trim(aefile),'AEBC FILE')
      do n = 1 , naetime
        aetime(n) = timeval2date(xtime(n),timeunits,timecal)
      end do
      deallocate(xtime)
      if ( .not. do_parallel_netcdf_in ) then
        call bcast(naetime)
        call bcast(aetime)
      end if
    end subroutine open_aerclima

    subroutine close_aerclima
      implicit none
      if ( .not. do_parallel_netcdf_in ) then
        if ( myid /= iocpu ) return
      end if
      if ( ncaec > 0 ) then
        if ( allocated(aetime) ) then
          deallocate(aetime)
        end if
        ncstatus = nf90_close(ncaec)
        call check_ok(__FILE__,__LINE__, &
               'Error Close AE file','AEBC FILE')
        ncaec = -1
      end if
    end subroutine close_aerclima

    subroutine check_ok(f,l,m1,mf)
      implicit none
      character(*) , intent(in) :: f, m1 , mf
      integer(ik4) , intent(in) :: l
      if ( ncstatus /= nf90_noerr ) then
        write (stderr,*) trim(m1)
        write (stderr,*) nf90_strerror(ncstatus)
        call fatal(f,l,trim(mf))
      end if
    end subroutine check_ok
    !
    !!FAB :  Read directly climatologies of optical properties
    !! (ext,ssa,asy)

    subroutine read_aeroppdata(idatex,m2r)
      implicit none
      type (rcm_time_and_date) , intent(in) :: idatex
      type(mod_2_rad) , intent(in) :: m2r
      logical , save :: lfirst = .true.
      logical :: dointerp
      real(rkx) , dimension(kth) :: opprnt
      real(rkx) :: xfac1 , xfac2 , odist
      type (rcm_time_and_date) :: imonmidd
      integer(ik4) :: iyear , imon , iday , ihour
      integer(ik4) :: i , j , k , kk , im1 , iy1 , im2 , iy2 , wn
      integer(ik4) , save :: ism , isy
      type (rcm_time_and_date) :: iref1 , iref2
      type (rcm_time_interval) :: tdif
      data ism /-1/
      data isy /-1/

      call split_idate(idatex,iyear,imon,iday,ihour)
      imonmidd = monmiddle(idatex)

      if (iyear < 1979 .and. iyear > 2020) then
        write (stderr,*) &
          'NO CLIMATIC AEROPP DATA AVAILABLE FOR ', iyear*100+imon
        return
      end if

      if ( lfirst ) then
        call grid_collect(m2r%xlon,alon,jce1,jce2,ice1,ice2)
        call grid_collect(m2r%xlat,alat,jce1,jce2,ice1,ice2)
        if ( myid == iocpu ) then
          call getfile(iyear,imon,ncid,3) ! open just clim vis for latlon
                                           ! reading
          call getmem1d(lat,1,clnlat,'aeropp:lat')
          call getmem1d(lon,1,clnlon,'aeropp:lon')
          call getmem3d(rdvar,1,clnlon,1,clnlat,1,clnlev, 'aerosol:rdvar')
          call getmem3d(hzivar,1,njcross,1,nicross,1,clnlev,'aerosol:hziext1')
          call init_aeroppdata(ncid,lat,lon)
          call h_interpolator_create(hint,lat,lon,alat,alon)
          call bcast(clnlev)
        else
          call bcast(clnlev)
        end if
        call getmem4d(plext1,jci1,jci2,ici1,ici2,1,clnlev, &
                      1,nacwb,'aerosol:plext1')
        call getmem4d(plssa1,jci1,jci2,ici1,ici2,1,clnlev, &
                      1,nacwb,'aerosol:plssa1')
        call getmem4d(plasy1,jci1,jci2,ici1,ici2,1,clnlev, &
                      1,nacwb,'aerosol:plasy1')
        call getmem4d(pldp1,jci1,jci2,ici1,ici2,1,clnlev, &
                      1,nacwb,'aerosol:pldp1')
        call getmem4d(pl1,jci1,jci2,ici1,ici2,1,clnlev, &
                      1,nacwb,'aerosol:pl1')
        call getmem4d(plext2,jci1,jci2,ici1,ici2,1,clnlev, &
                      1,nacwb,'aerosol:plext2')
        call getmem4d(plssa2,jci1,jci2,ici1,ici2,1,clnlev, &
                      1,nacwb,'aerosol:plssa2')
        call getmem4d(plasy2,jci1,jci2,ici1,ici2,1,clnlev, &
                      1,nacwb,'aerosol:plasy2')
        call getmem4d(pldp2,jci1,jci2,ici1,ici2,1,clnlev, &
                      1,nacwb,'aerosol:pldp2')
        call getmem4d(pl2,jci1,jci2,ici1,ici2,1,clnlev, &
                      1,nacwb,'aerosol:pl2')

        call getmem4d(sgext1,jci1,jci2,ici1,ici2,1,kth,1,nacwb,'aerosol:sgext1')
        call getmem4d(sgext2,jci1,jci2,ici1,ici2,1,kth,1,nacwb,'aerosol:sgext2')
        call getmem4d(sgssa1,jci1,jci2,ici1,ici2,1,kth,1,nacwb,'aerosol:sgssa1')
        call getmem4d(sgssa2,jci1,jci2,ici1,ici2,1,kth,1,nacwb,'aerosol:sgssa2')
        call getmem4d(sgasy1,jci1,jci2,ici1,ici2,1,kth,1,nacwb,'aerosol:sgasy1')
        call getmem4d(sgasy2,jci1,jci2,ici1,ici2,1,kth,1,nacwb,'aerosol:sgasy2')
        call getmem3d(zpr3d,jci1,jci2,ici1,ici2,1,kth,'aerosol:zpr3d')
        call getmem3d(zdzr3d,jci1,jci2,ici1,ici2,1,kth,'aerosol:zdz3d')
        !
        ! RRTMG radiative hat (up to kth)
        !
        do k = 1, kth-kz
          kk = n_prehlev-k+1
          zpr3d(:,:,k) = stdplevh(kk)*d_100
          kk = n_hreflev-k+1
          zdzr3d(:,:,k) = (stdhlevf(kk)-stdhlevf(kk-1))*d_1000
        end do
      end if
      im1 = imon
      iy1 = iyear
      im2 = imon
      iy2 = iyear
      if ( idatex > imonmidd ) then
        call inextmon(iy2,im2)
        iref1 = imonmidd
        iref2 = monmiddle(nextmon(idatex))
      else
        call iprevmon(iy1,im1)
        iref1 = monmiddle(prevmon(idatex))
        iref2 = imonmidd
      end if
      dointerp = .false.
      if ( ism /= im1 .or. isy /= iy1 ) then
        ism = im1
        isy = iy1
        dointerp = .true.
      end if
      if ( dointerp ) then
        if ( myid == iocpu ) then
          if ( myid == italk ) then
            write (stdout,*) 'Reading EXT,SSA,ASY Data...'
          end if
          if ( lfirst ) then
            do wn = 1 , nacwb
              call getfile(iy1,im1,ncid,wn)
              call readvar3d(ncid,'EXTTOT',rdvar)
              call remove_nans(rdvar,d_zero)
              call h_interpolate_cont(hint,rdvar,hzivar)
              call assignpnt(plext1,plvar,wn)
              call grid_distribute(hzivar,plvar,jci1,jci2,ici1,ici2,1,clnlev)
              call readvar3d(ncid,'SSATOT',rdvar)
              call remove_nans(rdvar,d_one)
              call h_interpolate_cont(hint,rdvar,hzivar)
              call assignpnt(plssa1,plvar,wn)
              call grid_distribute(hzivar,plvar,jci1,jci2,ici1,ici2,1,clnlev)
              call readvar3d(ncid,'GTOT',rdvar)
              call remove_nans(rdvar,0.2_rkx)
              call h_interpolate_cont(hint,rdvar,hzivar)
              call assignpnt(plasy1,plvar,wn)
              call grid_distribute(hzivar,plvar,jci1,jci2,ici1,ici2,1,clnlev)
              call readvar3d(ncid,'DELP',rdvar)
              call h_interpolate_cont(hint,rdvar,hzivar)
              call assignpnt(pldp1,plvar,wn)
              call grid_distribute(hzivar,plvar,jci1,jci2,ici1,ici2,1,clnlev)

              call getfile(iy2,im2,ncid,wn)
              call readvar3d(ncid,'EXTTOT',rdvar)
              call remove_nans(rdvar,d_zero)
              call h_interpolate_cont(hint,rdvar,hzivar)
              call assignpnt(plext2,plvar,wn)
              call grid_distribute(hzivar,plvar,jci1,jci2,ici1,ici2,1,clnlev)
              call readvar3d(ncid,'SSATOT',rdvar)
              call remove_nans(rdvar,d_one)
              call h_interpolate_cont(hint,rdvar,hzivar)
              call assignpnt(plssa2,plvar,wn)
              call grid_distribute(hzivar,plvar,jci1,jci2,ici1,ici2,1,clnlev)
              call readvar3d(ncid,'GTOT',rdvar)
              call remove_nans(rdvar,0.2_rkx)
              call h_interpolate_cont(hint,rdvar,hzivar)
              call assignpnt(plasy2,plvar,wn)
              call grid_distribute(hzivar,plvar,jci1,jci2,ici1,ici2,1,clnlev)
              call readvar3d(ncid,'DELP',rdvar)
              call h_interpolate_cont(hint,rdvar,hzivar)
              call assignpnt(pldp2,plvar,wn)
              call grid_distribute(hzivar,plvar,jci1,jci2,ici1,ici2,1,clnlev)
            end do ! clim wav band loop
          else
            plext1 = plext2
            plssa1 = plssa2
            plasy1 = plasy2
            do wn = 1 , nacwb
              call getfile(iy2,im2,ncid,wn)
              call readvar3d(ncid,'EXTTOT',rdvar)
              call remove_nans(rdvar,d_zero)
              call h_interpolate_cont(hint,rdvar,hzivar)
              call assignpnt(plext2,plvar,wn)
              call grid_distribute(hzivar,plvar,jci1,jci2,ici1,ici2,1,clnlev)
              call readvar3d(ncid,'SSATOT',rdvar)
              call remove_nans(rdvar,d_one)
              call h_interpolate_cont(hint,rdvar,hzivar)
              call assignpnt(plssa2,plvar,wn)
              call grid_distribute(hzivar,plvar,jci1,jci2,ici1,ici2,1,clnlev)
              call readvar3d(ncid,'GTOT',rdvar)
              call remove_nans(rdvar,0.2_rkx)
              call h_interpolate_cont(hint,rdvar,hzivar)
              call assignpnt(plasy2,plvar,wn)
              call grid_distribute(hzivar,plvar,jci1,jci2,ici1,ici2,1,clnlev)
              call readvar3d(ncid,'DELP',rdvar)
              call h_interpolate_cont(hint,rdvar,hzivar)
              call assignpnt(pldp2,plvar,wn)
              call grid_distribute(hzivar,plvar,jci1,jci2,ici1,ici2,1,clnlev)
            end do ! clim wav band loop
          end if
        else
          if ( lfirst ) then
            do wn = 1 , nacwb
              call assignpnt(plext1,plvar,wn)
              call grid_distribute(hzivar,plvar,jci1,jci2,ici1,ici2,1,clnlev)
              call assignpnt(plssa1,plvar,wn)
              call grid_distribute(hzivar,plvar,jci1,jci2,ici1,ici2,1,clnlev)
              call assignpnt(plasy1,plvar,wn)
              call grid_distribute(hzivar,plvar,jci1,jci2,ici1,ici2,1,clnlev)
              call assignpnt(pldp1,plvar,wn)
              call grid_distribute(hzivar,plvar,jci1,jci2,ici1,ici2,1,clnlev)
              call assignpnt(plext2,plvar,wn)
              call grid_distribute(hzivar,plvar,jci1,jci2,ici1,ici2,1,clnlev)
              call assignpnt(plssa2,plvar,wn)
              call grid_distribute(hzivar,plvar,jci1,jci2,ici1,ici2,1,clnlev)
              call assignpnt(plasy2,plvar,wn)
              call grid_distribute(hzivar,plvar,jci1,jci2,ici1,ici2,1,clnlev)
              call assignpnt(pldp2,plvar,wn)
              call grid_distribute(hzivar,plvar,jci1,jci2,ici1,ici2,1,clnlev)
            end do
          else
            plext1 = plext2
            plssa1 = plssa2
            plasy1 = plasy2
            do wn = 1 , nacwb
              call assignpnt(plext2,plvar,wn)
              call grid_distribute(hzivar,plvar,jci1,jci2,ici1,ici2,1,clnlev)
              call assignpnt(plssa2,plvar,wn)
              call grid_distribute(hzivar,plvar,jci1,jci2,ici1,ici2,1,clnlev)
              call assignpnt(plasy2,plvar,wn)
              call grid_distribute(hzivar,plvar,jci1,jci2,ici1,ici2,1,clnlev)
              call assignpnt(pldp2,plvar,wn)
              call grid_distribute(hzivar,plvar,jci1,jci2,ici1,ici2,1,clnlev)
            end do
          end if
        end if
        ! The pressure at the MERRA top is a fixed constant:
        !
        !                 PTOP = 0.01 hPa (1 Pa)
        !
        ! Pressures at edges should be computed by summing the DELP
        ! pressure thickness starting at PTOP.
        ! A representative pressure for the layer can then be obtained
        ! from these.
        ! Here just use linear av for now, should be improved !!
        ! MERRA grid is top down
        if ( lfirst ) then
          do wn = 1 , nacwb
            pl1(:,:,1,wn)  = d_one + pldp1(:,:,1,wn)*d_half
            do k = 2 , clnlev
              pl1(:,:,k,wn) = pl1(:,:,k-1,wn) + &
                (pldp1(:,:,k-1,wn) + pldp1(:,:,k,wn))*d_half
            end do
          end do
          do wn = 1 , nacwb
            pl2(:,:,1,wn)  = d_one + pldp2(:,:,1,wn)*d_half
            do k = 2 , clnlev
              pl2(:,:,k,wn) = pl2(:,:,k-1,wn) + &
                (pldp2(:,:,k-1,wn) + pldp2(:,:,k,wn))*d_half
            end do
          end do
        else
          pl1 = pl2
          do wn = 1 , nacwb
            pl2(:,:,1,wn)  = d_one + pldp2(:,:,1,wn)*d_half
            do k = 2 , clnlev
              pl2(:,:,k,wn) = pl2(:,:,k-1,wn) + &
                (pldp2(:,:,k-1,wn) + pldp2(:,:,k,wn))*d_half
            end do
          end do
        end if
      end if ! end of interp
      !
      ! VERTICAL Interpolation
      !
      ! Model pressure levels (in Pa)
      !
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kz )
        zpr3d(j,i,kth-kz+k) = m2r%phatms(j,i,k)
        zdzr3d(j,i,kth-kz+k) = m2r%deltaz(j,i,k)
      end do
      do wn = 1 , nacwb
        call assignpnt(pl1,p1,wn)
        call assignpnt(pl2,p2,wn)
        call assignpnt(plext1,plvar,wn)
        call assignpnt(sgext1,sgvar,wn)
        call intp1(sgvar,plvar,zpr3d,p1,jci1,jci2,ici1,ici2, &
                   kth,clnlev,0.7_rkx,0.7_rkx,0.4_rkx)
        call assignpnt(plext2,plvar,wn)
        call assignpnt(sgext2,sgvar,wn)
        call intp1(sgvar,plvar,zpr3d,p2,jci1,jci2,ici1,ici2, &
                   kth,clnlev,0.7_rkx,0.7_rkx,0.4_rkx)
        call assignpnt(plssa1,plvar,wn)
        call assignpnt(sgssa1,sgvar,wn)
        call intp1(sgvar,plvar,zpr3d,p1,jci1,jci2,ici1,ici2, &
                   kth,clnlev,0.7_rkx,0.7_rkx,0.4_rkx)
        call assignpnt(plssa2,plvar,wn)
        call assignpnt(sgssa2,sgvar,wn)
        call intp1(sgvar,plvar,zpr3d,p2,jci1,jci2,ici1,ici2, &
                   kth,clnlev,0.7_rkx,0.7_rkx,0.4_rkx)
        call assignpnt(plasy1,plvar,wn)
        call assignpnt(sgasy1,sgvar,wn)
        call intp1(sgvar,plvar,zpr3d,p1,jci1,jci2,ici1,ici2, &
                   kth,clnlev,0.7_rkx,0.7_rkx,0.4_rkx)
        call assignpnt(plasy2,plvar,wn)
        call assignpnt(sgasy2,sgvar,wn)
        call intp1(sgvar,plvar,zpr3d,p2,jci1,jci2,ici1,ici2, &
                   kth,clnlev,0.7_rkx,0.7_rkx,0.4_rkx)
      end do
      tdif = idatex-iref1
      xfac1 = real(tohours(tdif),rkx)
      tdif = idatex-iref2
      xfac2 = real(tohours(tdif),rkx)
      odist = xfac1 - xfac2
      xfac1 = xfac1 / odist
      xfac2 = d_one - xfac1
      !
      ! Important :  radiation schemes expect AOD per layer, calculated
      ! from extinction
      !
      do concurrent ( j = jci1:jci2, i = ici1:ici2, k = 1:kth , wn = 1:nacwb )
        extprof(j,i,k,wn) = (sgext1(j,i,k,wn)*xfac2 + &
                             sgext2(j,i,k,wn)*xfac1) * &
                             zdzr3d(j,i,k)
        ssaprof(j,i,k,wn) = (sgssa1(j,i,k,wn)*xfac2 + &
                             sgssa2(j,i,k,wn)*xfac1)
        asyprof(j,i,k,wn) = (sgasy1(j,i,k,wn)*xfac2 + &
                             sgasy2(j,i,k,wn)*xfac1)
      end do
      if ( myid == italk .and. dointerp ) then
        do k = 1 , kth
          opprnt(k) = extprof(jci1,ici1,k,3) !!vis band in MERRA clim
        end do
        call vprntv(opprnt,kth,'Updated VIS ext profile')
        do k = 1 , kth
          opprnt(k) = ssaprof(jci1,ici1,k,3)
        end do
        call vprntv(opprnt,kth,'Updated VIS ssa profile')
        do k = 1 , kth
          opprnt(k) = asyprof(jci1,ici1,k,3)
        end do
        call vprntv(opprnt,kth,'Updated VIS asy profile')
      end if
      lfirst = .false.
    end subroutine read_aeroppdata

    subroutine inextmon(iyear,imon)
      implicit none
      integer(ik4) , intent(inout) :: iyear , imon
      imon = imon + 1
      if ( imon > 12 ) then
        imon = 1
        iyear = iyear + 1
      end if
    end subroutine inextmon

    subroutine iprevmon(iyear,imon)
      implicit none
      integer(ik4) , intent(inout) :: iyear , imon
      imon = imon - 1
      if ( imon < 1 ) then
        imon = 12
        iyear = iyear - 1
      end if
    end subroutine iprevmon

    subroutine init_aeroppdata(ncid,lat,lon)
      implicit none
      integer(ik4) , intent(in) :: ncid
      real(rkx) , intent(inout) , dimension(:) :: lat , lon
      call readvar1d(ncid,'lat',lat)
      call readvar1d(ncid,'lon',lon)
    end subroutine init_aeroppdata

    subroutine readvar1d(ncid,vname,val)
      implicit none
      integer(ik4) , intent(in) :: ncid
      character(len=*) , intent(in) :: vname
      real(rkx) , intent(out) , dimension(:) :: val
      integer(ik4) :: icvar , iret
      iret = nf90_inq_varid(ncid,vname,icvar)
      if ( iret /= nf90_noerr ) then
        write (stderr, *) nf90_strerror(iret)
        call fatal(__FILE__,__LINE__,'CANNOT READ FROM AEROPP FILE')
      end if
      iret = nf90_get_var(ncid,icvar,val)
      if ( iret /= nf90_noerr ) then
        write (stderr, *) nf90_strerror(iret)
        call fatal(__FILE__,__LINE__,'CANNOT READ FROM AEROPP FILE')
      end if
    end subroutine readvar1d

    subroutine readvar3d(ncid,vname,val)
      implicit none
      integer(ik4) , intent(in) :: ncid
      character(len=*) , intent(in) :: vname
      real(rkx) , intent(out) , dimension(:,:,:) :: val
      integer(ik4) , save :: icvar
      integer(ik4) , save , dimension(3) :: istart , icount
      integer(ik4) :: iret
      data icvar /-1/
      data istart  /  1 ,  1 ,  1/

      icount(1) = clnlon
      icount(2) = clnlat
      icount(3) = clnlev
      iret = nf90_inq_varid(ncid,vname,icvar)
      if ( iret /= nf90_noerr ) then
        write (stderr, *) nf90_strerror(iret)
        call fatal(__FILE__,__LINE__,'CANNOT READ FROM AEROPPCLIM FILE')
      end if
      iret = nf90_get_var(ncid,icvar,val,istart,icount)
      if ( iret /= nf90_noerr ) then
        write (stderr, *) nf90_strerror(iret)
        call fatal(__FILE__,__LINE__,'CANNOT READ FROM AEROPPCLIM FILE')
      end if
    end subroutine readvar3d
    !
    ! SUBROUTINE AEROPPT
    !
    subroutine aeroppt(rh,pint,n1,n2)
      implicit none
      !
      ! Interface pressure, relative humidity
      !
      integer(ik4) , intent(in) :: n1 , n2
      real(rkx) , intent(in) , pointer , dimension(:,:) :: pint
      real(rkx) , intent(in) , pointer , dimension(:,:) :: rh

      integer(ik4) :: n , l , ibin , jbin , itr , k1 , k2 , ns
      integer(ik4) :: j , i , k , kk , visband
      real(rkx) :: uaerdust , qabslw , rh0
      real(rkx) , dimension(13) :: wavn ! central sw wave number for rrtm
      real(rkx) , dimension(4) :: wavncl ! central sw wave number for clim data
      real(rkx) , dimension(4) :: src
      real(rkx) , dimension(13) :: dest
      !-
      ! Aerosol forced by Optical Properties Climatology
      ! distinguish between standard scheme and RRTM scheme
      ! Rq extprof has been scaled for layer height before
      ! ( represents here the layer AOD )
      ! Directly force aerosol properties passed to radiation driver and exit

      if ( iclimaaer == 2 ) then
        aertrlw(:,:,:) = d_one
        if ( irrtm == 1 ) then
          do ns = 1 , 13
            wavn(ns) = 0.5_rkx*(wavnm2(ns) + wavnm1(ns))
          end do
          !
          ! wavncl are the climatological spectral bands 3,6,10,13
          ! MIGHT EVOLVE IN FUTURE
          !
          wavncl(1)  = wavn(3)
          wavncl(2)  = wavn(6)
          wavncl(3)  = wavn(10)
          wavncl(4)  = wavn(13)
          do k = 1 , kth
            do i = ici1 , ici2
              do j = jci1 , jci2
                n = (j-jci1)+(i-ici1)*nj+1
                ! Spectral interpolation for all RRTM band
                ! Special band 14 is left aside. Use constant extrapolation
                ! between clim data wn points ( check interp1d code in Share)
                ! to avoid negative values at wn = 1
                do ns = 1 , 4
                  src(ns) = extprof(j,i,k,ns)
                end do
                call interp1d(wavncl,src,wavn,dest,1._rkx,1._rkx,0._rkx)
                do ns = 1 , 13
                  tauxar3d(n,k,ns) = min(max(dest(ns),0._rkx),1.0_rkx)
                end do
                do ns = 1 , 4
                  src(ns) = ssaprof(j,i,k,ns)
                end do
                call interp1d(wavncl,src,wavn,dest,1._rkx,1._rkx,0._rkx)
                do ns = 1 , 13
                  tauasc3d(n,k,ns) = min(max(dest(ns),0._rkx),0.999_rkx)
                end do
                do ns = 1 , 4
                  src(ns) = asyprof(j,i,k,ns)
                end do
                call interp1d(wavncl,src,wavn,dest,1._rkx,1._rkx,0._rkx)
                do ns = 1 , 13
                  gtota3d(n,k,ns) = min(max(dest(ns),0._rkx),1.0_rkx)
                end do
                ! The LW prop assumed to be spectrally constant for now
                ! MIGHT CHANGE IN FUTURE
                do ns = 1 , nbndlw
                  tauxar3d_lw(n,k,ns) = extprof(j,i,k,5)
                end do
              end do
            end do
          end do
        else if ( irrtm == 0 ) then
          visband = 8
          ns = visband
          tauxar(ns,:) = d_zero
          tauasc(ns,:) = d_zero
          gtota(ns,:) = d_zero
          ftota(ns,:) = d_zero
          ! adapt the clim vert grid (1 to kth) to the standard
          ! rad grid ( 0 to kz) first Treat the top radiative layer
          ! tauxar3d(ns,0,n)
          if ( kth > kz ) then
            do i = ici1 , ici2
              do j = jci1 , jci2
                !
                ! index 3 correponds to clim vis band
                ! MIGHT CHANGE IN FUTURE
                !
                n = (j-jci1)+(i-ici1)*nj+1
                tauxar3d(ns,0,n) = max(sum(extprof(j,i,1:kth-kz,3)),0.0_rkx)
                tauasc3d(ns,0,n) = tauxar3d(ns,0,n) * &
                  max(sum(ssaprof(j,i,1:kth-kz,3))/real(kth-kz,rkx),0.0_rkx)
                gtota3d(ns,n,0) = max(sum(asyprof(j,i,1:kth-kz,3)) / &
                         real(kth-kz,rkx),0.0_rkx)
                ftota3d(ns,0,n) = tauasc3d(ns,0,n) * gtota3d(ns,0,n)**2
                gtota3d(ns,0,n) = tauasc3d(ns,0,n) * gtota3d(ns,0,n)
              end do
            end do
          end if
          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
                ! already scaled for layer height
                ! grid is top down
                ! FAB : index 2 is the vis band in MERRA aerclim
                ! MIGHT CHANGE
                n = (j-jci1)+(i-ici1)*nj+1
                tauxar3d(ns,k,n) = max(extprof(j,i,kth-kz+k,3),0.0_rkx)
                tauasc3d(ns,k,n) = max(ssaprof(j,i,kth-kz+k,3),0.0_rkx)
                gtota3d(ns,k,n) = max(asyprof(j,i,kth-kz+k,3),0.0_rkx)
                ! here the standard scheme expect layer scaled quantity
                tauasc3d(ns,k,n) = tauasc3d(ns,k,n) * tauxar3d(ns,k,n)
                ftota3d(ns,k,n) = gtota3d(ns,k,n)**2 * tauasc3d(ns,k,n)
                gtota3d(ns,k,n) = gtota3d(ns,k,n) * tauasc3d(ns,k,n)
                tauxar(ns,n) = tauxar(ns,n) + tauxar3d(ns,k,n)
                tauasc(ns,n) = tauasc(ns,n) + tauasc3d(ns,k,n)
                gtota(ns,n)  = gtota(ns,n)  + gtota3d(ns,k,n)
                ftota(ns,n)  = ftota(ns,n)  + ftota3d(ns,k,n)
              end do
            end do
          end do
        end if
        return  ! important
      else if ( iclimaaer == 3 ) then
        ! Limit only to SW (?)
        aertrlw(:,:,:) = d_one
        if ( irrtm == 1 ) then
          do n = 1 , nband
            do k = 1 , kth
              do i = 1 , npoints
                tauxar3d(i,k,n) = extprofr4(i,k,n)
                tauasc3d(i,k,n) = ssaprofr4(i,k,n)
                gtota3d(i,k,n)  = asyprofr4(i,k,n)
              end do
            end do
          end do
        else if ( irrtm == 0 ) then
          do n = 1 , nband
            do k = 1 , kth
              do i = 1 , npoints
                ! already scaled for layer height
                tauxar3d(n,k,i) = extprofr4(i,k,n)
                ! here the standard scheme expect layer scaled quantity
                tauasc3d(n,k,i) = ssaprofr4(i,k,n) * tauxar3d(n,k,i)
                gtota3d(n,k,i)  = asyprofr4(i,k,n) * &
                       ssaprofr4(i,k,n) * tauxar3d(n,k,i)
                ftota3d(n,k,i)  = asyprofr4(i,k,n)**2 * &
                       ssaprofr4(i,k,n) * tauxar3d(n,k,i)

                ! define also tauxar for std scheme clear sky diagnostics
                tauxar(n,i) = tauxar(n,i) + tauxar3d(n,k,i)
                tauasc(n,i) = tauasc(n,i) + tauasc3d(n,k,i)
                ftota(n,i) =  ftota(n,i)  + ftota3d(n,k,i)
                gtota(n,i) =  gtota(n,i)  + gtota3d(n,k,i)
              end do
            end do
          end do
        end if
        return  ! important
      end if
      ! Exit if you don't want aerosol optical properties calculated from
      ! concentrations (either intractive aerosol, or prescribed concentration
      ! climatology)

      if ( ichem /= 1 .and. iclimaaer /= 1 ) then
        tauxar(:,:)     = d_zero
        tauasc(:,:)     = d_zero
        gtota(:,:)      = d_zero
        ftota(:,:)      = d_zero
        tauxar3d(:,:,:) = d_zero
        tauasc3d(:,:,:) = d_zero
        gtota3d(:,:,:)  = d_zero
        ftota3d(:,:,:)  = d_zero
        aertrlw (:,:,:) = d_one
        if ( irrtm == 1 ) then
          tauxar3d_lw(:,:,:) = d_zero
        end if
        return
      end if
      !
      !Calculate aerosol properties passed to radiation from concentrations
      !(interactive or prescribed from climatology)
      !
      tx = d_zero
      wa = d_zero
      ga = d_zero
      fa = d_zero
      !
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      !
      !   Melange externe
      !
      !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
      !
      !   Spectral loop
      !
      do concurrent ( n = n1:n2, k = 1:kz )
        path(n,k) = (pint(n,k+1)-pint(n,k))*regravgts
      end do

      do ns = 1 , nband
        tauxar(:,ns) = d_zero
        tauasc(:,ns) = d_zero
        gtota(:,ns) = d_zero
        ftota(:,ns) = d_zero

        tauxar3d(:,:,ns) = d_zero
        tauasc3d(:,:,ns) = d_zero
        gtota3d(:,:,ns) = d_zero
        ftota3d(:,:,ns) = d_zero

        uaer(:,:,:) = d_zero
        tx(:,:,:) = d_zero
        wa(:,:,:) = d_zero
        ga(:,:,:) = d_zero
        fa(:,:,:) = d_zero
        utaer(:,:) = d_zero
        tauaer(:,:) = d_zero
        waer(:,:) = d_zero
        gaer(:,:) = d_zero
        faer(:,:) = d_zero
        !
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !
        ! calculate optical properties of each aerosol component
        ! correct for rh if needed
        !
        !CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
        !
        ibin = 0
        jbin = 0
        do itr = 1 , ntr
          if ( chtrname(itr) == 'XXXXX') then
            continue
          end if
          if ( chtrname(itr)(1:4) == 'DUST') then
            ibin = ibin + 1
            do k = 1 , kz
              do n = n1 , n2
                uaer(n,k,itr) = aermmr(n,k,itr)*path(n,k)
                tx(n,k,itr) = d10e5*uaer(n,k,itr)*ksdust(ns,ibin)
                wa(n,k,itr) = wsdust(ns,ibin)
                ga(n,k,itr) = gsdust(ns,ibin)
                fa(n,k,itr) = gsdust(ns,ibin)*gsdust(ns,ibin)
              end do
            end do
          else if ( chtrname(itr)(1:3) == 'SO4' .or.  &
                    chtrname(itr)(1:5) == 'H2SO4'.or. &
                    chtrname(itr)(1:4) == 'ANO3' .or. &
                    chtrname(itr)(1:4) == 'ANH4' ) then
            do k = 1 , kz
              do n = n1 , n2
                rh0 = min(0.97_rkx,max(0.0_rkx,rh(n,k)))
                ! maximum limit for effect on sulfate extinction
                uaer(n,k,itr) = aermmr(n,k,itr)*path(n,k)
                tx(n,k,itr) = d10e5*uaer(n,k,itr) *    &
                  ksbase(ns)*exp(kscoef(ns,1) +        &
                     kscoef(ns,2)/(rh0+kscoef(ns,3)) + &
                     kscoef(ns,4)/(rh0+kscoef(ns,5)))
                wa(n,k,itr) = d_one - wsbase(ns) * exp(wscoef(ns,1) + &
                  wscoef(ns,2) / (rh0+wscoef(ns,3)) +             &
                  wscoef(ns,4) / (rh0+wscoef(ns,5)))
                ga(n,k,itr) = gsbase(ns) * exp(gscoef(ns,1) + &
                  gscoef(ns,2) / (rh0+gscoef(ns,3)) +     &
                  gscoef(ns,4) / (rh0+gscoef(ns,5)))
                fa(n,k,itr) = ga(n,k,itr)*ga(n,k,itr)
              end do
            end do
          else if ( chtrname(itr)(1:4) == 'SSLT' ) then
            jbin = jbin+1
            do k = 1 , kz
              do n = n1 , n2
                rh0 = min(0.99_rkx,max(0.0_rkx,rh(n,k)))
                do l = 1 , 7
                  if ( rh0 > rhp(l) .and. rh0 <= rhp(l+1) ) then
                    ! FAB : test according to li et al., ksslt cannot exceed 1.3
                    ! quick fix for now, update parameterisation to LI et al,
                    ! ACP 2008 in a near future
                    kssslt(ns,jbin) = min(ksslt(ns,jbin,l),1.2_rkx)
                    gssslt(ns,jbin) = gsslt(ns,jbin,l)
                    wssslt(ns,jbin) = wsslt(ns,jbin,l)
                  end if
                end do
                uaer(n,k,itr) = aermmr(n,k,itr)*path(n,k)
                tx(n,k,itr) = d10e5*uaer(n,k,itr)*kssslt(ns,jbin)
                wa(n,k,itr) = wssslt(ns,jbin)
                ga(n,k,itr) = gssslt(ns,jbin)
                fa(n,k,itr) = gssslt(ns,jbin)*gssslt(ns,jbin)
              end do
            end do
          else if ( chtrname(itr)(1:5) == 'OC_HL' ) then
            do k = 1 , kz
              do n = n1 , n2
                uaer(n,k,itr) = aermmr(n,k,itr)*path(n,k)
                rh0 = min(0.99_rkx,max(0.0_rkx,rh(n,k)))
                ! Humidity effect !
                tx(n,k,itr) = d10e5*uaer(n,k,itr)*ksoc_hl(ns) * &
                              (d_one-rh0)**(-0.25_rkx)
                wa(n,k,itr) = wsoc_hl(ns)
                ga(n,k,itr) = gsoc_hl(ns)
                fa(n,k,itr) = ga(n,k,itr)*ga(n,k,itr)
              end do
            end do
          else if ( chtrname(itr)(1:5) == 'BC_HL' ) then
            do k = 1 , kz
              do n = n1 , n2
                uaer(n,k,itr) = aermmr(n,k,itr)*path(n,k)
                rh0 = min(0.99_rkx,max(0.0_rkx,rh(n,k)))
                ! Humidity effect !
                tx(n,k,itr) = d10e5*uaer(n,k,itr)*ksbc_hl(ns) * &
                              (d_one-rh0)**(-0.20_rkx)
                wa(n,k,itr) = wsbc_hl(ns)
                ga(n,k,itr) = gsbc_hl(ns)
                fa(n,k,itr) = ga(n,k,itr)*ga(n,k,itr)
              end do
            end do
          else if ( chtrname(itr)(1:5) == 'OC_HB' ) then
            do k = 1 , kz
              do n = n1 , n2
                uaer(n,k,itr) = aermmr(n,k,itr)*path(n,k)
                tx(n,k,itr) = d10e5*uaer(n,k,itr)*ksoc_hb(ns)
                wa(n,k,itr) = wsoc_hb(ns)
                ga(n,k,itr) = gsoc_hb(ns)
                fa(n,k,itr) = gsoc_hb(ns)*gsoc_hb(ns)
              end do
            end do
          else if ( chtrname(itr)(1:5) == 'BC_HB' ) then
            do k = 1 , kz
              do n = n1 , n2
                uaer(n,k,itr) = aermmr(n,k,itr)*path(n,k)
                ! Absorbing aerosols (soot type)
                tx(n,k,itr) = d10e5*uaer(n,k,itr)*ksbc_hb(ns)
                wa(n,k,itr) = wsbc_hb(ns)
                ga(n,k,itr) = gsbc_hb(ns)
                fa(n,k,itr) = gsbc_hb(ns)*gsbc_hb(ns)
              end do
            end do
          else if ( chtrname(itr)(1:3) == 'SM1' ) then
            do k = 1 , kz
              do n = n1 , n2
                uaer(n,k,itr) = aermmr(n,k,itr)*path(n,k)
                rh0 = min(0.99_rkx,max(0.0_rkx,rh(n,k)))
                tx(n,k,itr) = d10e5*uaer(n,k,itr)*kssm1(ns) * &
                              (d_one-rh0)**(-0.15_rkx)
                wa(n,k,itr) = wssm1(ns)
                ga(n,k,itr) = gssm1(ns)
                fa(n,k,itr) = gssm1(ns)*gssm1(ns)
              end do
            end do
          else if ( chtrname(itr)(1:3) == 'SM2' ) then
            do k = 1 , kz
              do n = n1 , n2
                uaer(n,k,itr) = aermmr(n,k,itr)*path(n,k)
                rh0 = min(0.99_rkx,max(0.0_rkx,rh(n,k)))
                tx(n,k,itr) = d10e5*uaer(n,k,itr)*kssm2(ns) * &
                              (d_one-rh0)**(-0.25_rkx)
                wa(n,k,itr) = wssm2(ns)
                ga(n,k,itr) = gssm2(ns)
                fa(n,k,itr) = gssm2(ns)*gssm2(ns)
              end do
            end do
          end if
        end do ! end tracer loop
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! calculate optical properties of the aerosol mixture
        !             passed to radiation scheme
        !ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
        ! optical properties for the clear sky diagnostic standard scheme
        do itr = 1 , ntr
          do k = 1 , kz
            do n = n1 , n2
              utaer(n,itr) = utaer(n,itr) + uaer(n,k,itr)
              tauaer(n,itr) = tauaer(n,itr) + tx(n,k,itr)
              waer(n,itr) = waer(n,itr) + wa(n,k,itr)*uaer(n,k,itr)
              gaer(n,itr) = gaer(n,itr) + ga(n,k,itr)*uaer(n,k,itr)
              faer(n,itr) = faer(n,itr) + fa(n,k,itr)*uaer(n,k,itr)
            end do
          end do
        end do
        do itr = 1 , ntr
          do n = n1 , n2
            if ( utaer(n,itr) <= minimum_utaer ) utaer(n,itr) = minimum_utaer
            waer(n,itr) = waer(n,itr)/utaer(n,itr)
            gaer(n,itr) = gaer(n,itr)/utaer(n,itr)
            faer(n,itr) = faer(n,itr)/utaer(n,itr)
          end do
        end do
        !
        ! Calculate the EXTERNAL Mixing of aerosols
        ! melange externe
        !
        ! only for climatic feedback allowed

        if ( irrtm == 0 ) then
          do itr = 1 , ntr
            do k = 0 , kz
              do n = n1 , n2
                tauxar3d(n,k,ns) = tauxar3d(n,k,ns) + tx(n,k,itr)
                tauasc3d(n,k,ns) = tauasc3d(n,k,ns) + tx(n,k,itr)*wa(n,k,itr)
                gtota3d(n,k,ns) = gtota3d(n,k,ns) + ga(n,k,itr) * &
                                  tx(n,k,itr)*wa(n,k,itr)
                ftota3d(n,k,ns) = ftota3d(n,k,ns) + fa(n,k,itr) * &
                                  tx(n,k,itr)*wa(n,k,itr)
              end do
            end do
          end do
          do k = 0 , kz
            do n = n1 , n2
              !consider a minimal extinction and reflectivity background
              if ( tauxar3d(n,k,ns) < 1.E-10_rkx ) then
                tauxar3d(n,k,ns) = 1.E-10_rkx
                tauasc3d(n,k,ns) = 0.999999_rkx * tauxar3d(n,k,ns)
                gtota3d(n,k,ns) = 0.5_rkx * tauasc3d(n,k,ns)
                ftota3d(n,k,ns) = 0.5_rkx * gtota3d(n,k,ns)
              end if
            end do
          end do
          !
          ! Clear sky (always calcuated if ichdir >=1 for
          ! diagnostic radiative forcing)
          !
          do itr = 1 , ntr
            do n = n1 , n2
              tauxar(n,ns) = tauxar(n,ns) + tauaer(n,itr)
              if (waer(n,itr) > minimum_waer) then
                tauasc(n,ns) = tauasc(n,ns) + tauaer(n,itr)*waer(n,itr)
              end if
              if (gaer(n,itr) > minimum_gaer .and.  &
                waer(n,itr) > minimum_gaer) then
                gtota(n,ns) = gtota(n,ns) + gaer(n,itr) * &
                                tauaer(n,itr)*waer(n,itr)
                ftota(n,ns) = ftota(n,ns) + faer(n,itr) * &
                                tauaer(n,itr)*waer(n,itr)
              end if
            end do
          end do
          ! in the case RRTM expect the layer extinction, and effective
          ! SSA and asym relative to the mixture
        else if ( irrtm == 1 ) then
          do itr = 1 , ntr
            do k = 1 , kz
              do n = n1 , n2
                kk = kth -kz + k
                tauxar3d(n,kk,ns) = tauxar3d(n,kk,ns) + tx(n,k,itr)
                tauasc3d(n,kk,ns) = tauasc3d(n,kk,ns) + tx(n,k,itr)*wa(n,k,itr)
                gtota3d(n,kk,ns) = gtota3d(n,kk,ns) + ga(n,k,itr) * &
                                  tx(n,k,itr)*wa(n,k,itr)
              end do
            end do
          end do
          do k = kth -kz+1 , kth
            do n = n1 , n2
              !consider a minimal extinction background
              if ( tauxar3d(n,k,ns) > 1.E-10_rkx ) then
                tauasc3d(n,k,ns) = tauasc3d(n,k,ns) / tauxar3d(n,k,ns)
                gtota3d(n,k,ns) = gtota3d(n,k,ns) / &
                  (tauasc3d(n,k,ns)*tauxar3d(n,k,ns))
              else
                tauxar3d(n,k,ns) = 1.E-10_rkx
                tauasc3d(n,k,ns) = 0.999999_rkx
                gtota3d(n,k,ns) = 0.5_rkx
              end if
            end do
          end do
          ! radiative hat
          tauxar3d(n1:n2,0:kth -kz,ns) = 1.E-10_rkx
          tauasc3d(n1:n2,0:kth -kz,ns) = 0.999999_rkx
          gtota3d(n1:n2 ,0:kth -kz,ns) = 0.5_rkx
        end if
      end do ! end spectral loop

      ! DUST LW emissivity
      !
      if ( irrtm == 0 ) then
        ! qabslw = absorption coeff between k1 and  k2 (m2.g-1) in the LW :
        qabslw = d_r10
        ! initialisation à 1 = perfect transmittivity
        aertrlw (:,:,:) = d_one
        !
        do itr = 1 , ntr
          if ( chtrname(itr)(1:4) == 'DUST' ) then
            do k1 = 1 , kzp1
              do k2 = 1 , kzp1
                do n = n1 , n2
                  if ( k1 == k2 ) aertrlw(n,k1,k2) = d_one
                  ! aerosol path btw k1 and k2 flux level
                  uaerdust = d_zero
                  if ( k1<k2 ) then
                    uaerdust =  uaerdust + d10e5 * sum(uaer(n,k1:k2-1,itr))
                    aertrlw(n,k1,k2) = exp(-fiveothree * qabslw * uaerdust)
                  else if ( k1>k2 ) then
                    uaerdust =  uaerdust + d10e5 * sum(uaer(n,k2:k1-1,itr))
                    aertrlw(n,k1,k2) = exp(-fiveothree * qabslw * uaerdust)
                  end if
                end do
              end do
            end do
          end if
        end do
      else if ( irrtm == 1 ) then
        ! in this case use directly the LW extinction.
        do ns = 1 , nbndlw
          tauxar3d_lw(:,:,ns) = 1.0E-10_rkx
          ibin = 0
          do itr = 1 , ntr
            if ( chtrname(itr)(1:4) == 'DUST') then
              ibin = ibin + 1
              do k = 1 , kz
                do n = n1 , n2
                  kk = kth - kz +k
                  uaer(n,k,itr) = aermmr(n,k,itr)*path(n,k)
                  tx(n,k,itr) = d10e5*uaer(n,k,itr)*ksdust_lw(ns,ibin)
                  ! add the extinction for every bins
                  tauxar3d_lw(n,kk,ns) = tauxar3d_lw(n,kk,ns) + tx(n,k,itr)
                end do
              end do
            end if
          end do
        end do
      end if
    end subroutine aeroppt

    subroutine cmip6_plume_profile(x,m2r)
      implicit none
      type (rcm_time_and_date) , intent(in) :: x
      type(mod_2_rad) , intent(in) :: m2r
      logical , save :: lfirst = .true.
      real(rkx) , dimension(kth) :: opprnt
      integer(ik4) :: ibin , i , j , k , kk , n
      integer(ik4) :: iy , im , id
      integer(ik4) , save :: idlast = -1
      real(rk8) :: year_fr

      if ( lfirst ) then
        if ( irrtm == 0 ) then
          do n = 1 , nband
            lambdaw(n) = (wavmin(n)+wavmax(n))*d_half*d_1000
          end do
        else if ( irrtm == 1 ) then
          do n = 1 , nband
            ! wavenumber is in cm-1 , convert to wavelenght in nm for mac-v2
            lambdaw(n) = (1.e7_rkx/wavnm1(n)+1.e7_rkx/wavnm2(n))*d_half
          end do
        end if
        ibin = 1
        do i = ici1 , ici2
          do j = jci1 , jci2
            altr4(ibin) = m2r%ht(j,i)*regrav
            latr4(ibin) = m2r%xlat(j,i)
            lonr4(ibin) = m2r%xlon(j,i)
            ibin = ibin + 1
          end do
        end do
        lfirst = .false.
      end if
      call split_idate(x,iy,im,id)
      if ( id /= idlast ) then
        do k = 1 , kz
          ibin = 1
          kk = kth-k+1
          do i = ici1 , ici2
            do j = jci1 , jci2
              z(ibin,kk)  = m2r%za(j,i,kz-k+1)
              dz(ibin,kk) = m2r%deltaz(j,i,kz-k+1)
              ibin = ibin + 1
            end do
          end do
        end do
        year_fr = real(iy) + real(yeardayfrac(x))/real(yeardays(iy,x%calendar))
        do n = 1 , nband
          call sp_aop_profile(macv2sp_hist,macv2sp_scen, &
                              kth,npoints,lambdaw(n), &
                              altr4,lonr4,latr4,year_fr,z,dz, &
                              dnovrnr4,extprofr4(:,:,n), &
                              ssaprofr4(:,:,n),asyprofr4(:,:,n))
        end do
        if ( myid == italk ) then
          write(stdout,*) 'Updating aerosol optical properties...'
          do k = 1 , kth
            opprnt(k) = extprofr4(1,k,10)
          end do
          call vprntv(opprnt,kth,'Updated VIS ext profile')
          do k = 1 , kth
            opprnt(k) = ssaprofr4(1,k,10)
          end do
          call vprntv(opprnt,kth,'Updated VIS ssa profile')
          do k = 1 , kth
            opprnt(k) = asyprofr4(1,k,10)
          end do
          call vprntv(opprnt,kth,'Updated VIS asy profile')
        end if
        idlast = id
      end if
    end subroutine cmip6_plume_profile

    subroutine getfile(year,month,ncid,wbclim)
      implicit none
      integer(ik4) , intent(in) :: year,month,wbclim
      integer(ik4) , intent(inout) ::  ncid
      character(len=256) :: infile
      integer(ik4) :: iret , idimid
      character(len=5) :: filnum
      if (wbclim == 1) filnum = 'wb3.'
      if (wbclim == 2) filnum = 'wb6.'
      if (wbclim == 3) filnum = 'wb10.'
      if (wbclim == 4) filnum = 'wb13.'
      if (wbclim == 5) filnum = 'wb19.' !

      write(infile,'(A,A,I4,I0.2,A)') &
        trim(radclimpath)//pthsep//'MERRA2_OPPMONTH_', &
        trim(filnum),year,month,'.nc'
      if ( ncid < 0 ) then
        iret = nf90_open(infile,nf90_nowrite,ncid)
        if ( iret /= nf90_noerr ) then
          write (stderr, *) nf90_strerror(iret), trim(infile)
          call fatal(__FILE__,__LINE__,'CANNOT OPEN AEROSOL OP.PROP CLIM FILE')
        end if
      else
        iret = nf90_close(ncid)
        if ( iret /= nf90_noerr ) then
          write (stderr, *) nf90_strerror(iret), trim(infile)
          call fatal(__FILE__,__LINE__,'CANNOT CLOSE FILE')
        end if
        iret = nf90_open(infile,nf90_nowrite,ncid)
        if ( iret /= nf90_noerr ) then
          write (stderr, *) nf90_strerror(iret), trim(infile)
          call fatal(__FILE__,__LINE__, &
                     'CANNOT OPEN AEROSOL OP.PROP CLIM FILE')
        end if
        if ( myid == italk ) then
          write(stdout,*) 'AEROPP file open ', trim(infile)
        end if
      end if
      ncstatus = nf90_inq_dimid(ncid,'lev',idimid)
      call check_ok(__FILE__,__LINE__, &
         'Error searching dimension lev in file '//trim(infile),'OPP FILE')
      ncstatus = nf90_inquire_dimension(ncid,idimid,len=clnlev)
      call check_ok(__FILE__,__LINE__, &
         'Error reading dimension lev in file '//trim(infile),'OPP FILE')
      ncstatus = nf90_inq_dimid(ncid,'lon',idimid)
      call check_ok(__FILE__,__LINE__, &
         'Error searching dimension lon in file '//trim(infile),'OPP FILE')
      ncstatus = nf90_inquire_dimension(ncid,idimid,len=clnlon)
      call check_ok(__FILE__,__LINE__, &
         'Error reading dimension lon in file '//trim(infile),'OPP FILE')
      ncstatus = nf90_inq_dimid(ncid,'lat',idimid)
      call check_ok(__FILE__,__LINE__, &
         'Error searching dimension lat in file '//trim(infile),'OPP FILE')
      ncstatus = nf90_inquire_dimension(ncid,idimid,len=clnlat)
      call check_ok(__FILE__,__LINE__, &
         'Error reading dimension lat in file '//trim(infile),'OPP FILE')
    end subroutine getfile

    subroutine remove_nans(val,set)
      implicit none
      real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: val
      real(rkx) , intent(in) :: set
      where ( is_nan(val) )
        val = set
      end where
    end subroutine remove_nans

end module mod_rad_aerosol

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

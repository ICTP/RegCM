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
  public :: allocate_mod_rad_aerosol , aermix , aeroppt
  public :: init_aerclima , read_aerclima , close_aerclima
  public :: init_aeroppdata , read_aeroppdata
  public :: cmip6_plume_profile
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
  real(rkx) , pointer , dimension(:,:,:) :: ext1 , ext2
  real(rkx) , pointer , dimension(:,:,:) :: extprof
  real(rkx) , pointer , dimension(:,:,:) :: ssaprof
  real(rkx) , pointer , dimension(:,:,:) :: asyprof
  real(rkx) , pointer , dimension(:,:,:) :: ssa1 , ssa2 , asy1 , asy2
  real(rkx) , pointer , dimension(:,:,:) :: ext , ssa , asy , zp3d , zdz3d
  real(rkx) , pointer , dimension(:,:,:) :: zpr3d , zdzr3d
  real(rkx) , pointer , dimension(:,:,:) :: yext , yssa , yasy, ydelp , yphcl
  real(rkx) , pointer , dimension(:,:,:) :: xext1 , xext2
  real(rkx) , pointer , dimension(:,:,:) :: xssa1 , xssa2
  real(rkx) , pointer , dimension(:,:,:) :: xasy1 , xasy2
  real(rkx) , pointer , dimension(:,:,:) :: xdelp1 , xdelp2
  type(h_interpolator) :: hint
  integer(ik4) :: ncid = -1
  integer(ik4) :: clnlon , clnlat , clnlev
  !
  integer(ik4) , parameter :: ncoefs = 5  ! Number of coefficients
  integer(ik4) , parameter :: nwav = 19
  integer(ik4) , parameter :: nih = 8

  integer(ik4) , parameter :: aerclima_ntr = 12
  integer(ik4) , parameter :: aerclima_nbin = 4

  character(len=6) , dimension(aerclima_ntr) :: aerclima_chtr
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
  integer(ik4) , private :: ii , jj ! coefficient index

  ! Sulfate param for standard scheme / works only with rad standard
  ! (Brieglieb et al.)
  real(rkx) , pointer , dimension(:,:,:) :: aermmr
  real(rkx) , dimension(nspi) :: gsbase  , ksbase , wsbase
  real(rkx) , dimension(nspi,ncoefs) :: gscoef , kscoef , wscoef

  ! optical properties for dust and org / for rrtm and standard scheme
  real(rkx) , dimension(nspi) ::  gsbc_hb_stand , gsbc_hl_stand ,   &
    gsoc_hb_stand , gsoc_hl_stand , ksbc_hb_stand , ksbc_hl_stand , &
    ksoc_hb_stand , ksoc_hl_stand , wsbc_hb_stand , wsbc_hl_stand , &
    wsoc_hb_stand , wsoc_hl_stand

  real(rkx) , dimension(nbndsw) :: gsbc_hb_rrtm , gsbc_hl_rrtm , &
    gsoc_hb_rrtm , gsoc_hl_rrtm , ksbc_hb_rrtm , ksbc_hl_rrtm ,  &
    ksoc_hb_rrtm , ksoc_hl_rrtm , wsbc_hb_rrtm , wsbc_hl_rrtm ,  &
    wsoc_hb_rrtm , wsoc_hl_rrtm , gssm1_rrtm , gssm2_rrtm ,      &
    kssm1_rrtm , kssm2_rrtm , wssm1_rrtm , wssm2_rrtm

  real(rkx) , dimension(nspi,4) :: gsdust_stand , ksdust_stand , wsdust_stand
  real(rkx), dimension(nspi,12) :: gsdust12_stand  , ksdust12_stand , &
    wsdust12_stand

  real(rkx) , dimension(nbndsw,4) :: gsdust_rrtm , ksdust_rrtm  , wsdust_rrtm
  real(rkx), dimension(nbndsw,12) :: gsdust12_rrtm , ksdust12_rrtm , &
    wsdust12_rrtm
  real(rkx), dimension(nbndlw,4) ::  ksdust_rrtm_lw
  real(rkx), dimension(nbndlw,12) ::  ksdust12_rrtm_lw

  ! sea salt oppt  param for standard scheme
  real(rkx) , dimension(nwav,2,nih) :: ksslt , wsslt , gsslt
  real(rkx) , dimension(nspi,2) :: gssslt , kssslt , wssslt
  !
  real(rkx) , dimension(8) :: rhp
  !
  ! Depth
  !
  real(rkx) , pointer , dimension(:,:) :: path
  !
  ! Background aerosol mass mixing ratio
  !
  real(rkx) , pointer , dimension(:,:) :: aermmb
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
  integer(ik4) :: ll , mm , nn
  integer(ik4) :: npoints
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
  data ksbase / 5.206e+0_rkx , 5.206e+0_rkx , 5.206e+0_rkx , 5.206e+0_rkx , &
                5.206e+0_rkx , 5.206e+0_rkx , 5.206e+0_rkx , 3.203e+0_rkx , &
                3.203e+0_rkx , 1.302e+0_rkx , 5.992e-1_rkx , 2.948e-1_rkx , &
                1.475e-1_rkx , 7.387e-2_rkx , 1.683e-1_rkx , 2.655e-1_rkx , &
                5.770e-2_rkx , 2.290e-1_rkx , 2.270e-1_rkx /

  data ((kscoef(ii,jj),jj=1,ncoefs),ii=1,nspi)/         1.126e+1_rkx , &
       -2.502e-1_rkx , -1.087e+0_rkx , -1.794e+2_rkx ,  1.556e+1_rkx , &
        1.126e+1_rkx , -2.502e-1_rkx , -1.087e+0_rkx , -1.794e+2_rkx , &
        1.556e+1_rkx ,  1.126e+1_rkx , -2.502e-1_rkx , -1.087e+0_rkx , &
       -1.794e+2_rkx ,  1.556e+1_rkx ,  1.126e+1_rkx , -2.502e-1_rkx , &
       -1.087e+0_rkx , -1.794e+2_rkx ,  1.556e+1_rkx ,  1.126e+1_rkx , &
       -2.502e-1_rkx , -1.087e+0_rkx , -1.794e+2_rkx ,  1.556e+1_rkx , &
        1.126e+1_rkx , -2.502e-1_rkx , -1.087e+0_rkx , -1.794e+2_rkx , &
        1.556e+1_rkx ,  1.126e+1_rkx , -2.502e-1_rkx , -1.087e+0_rkx , &
       -1.794e+2_rkx ,  1.556e+1_rkx ,  1.124e+1_rkx , -3.040e-1_rkx , &
       -1.088e+0_rkx , -1.776e+2_rkx ,  1.537e+1_rkx ,  1.124e+1_rkx , &
       -3.040e-1_rkx , -1.088e+0_rkx , -1.776e+2_rkx ,  1.537e+1_rkx , &
        1.222e+1_rkx , -3.770e-1_rkx , -1.089e+0_rkx , -1.898e+2_rkx , &
        1.504e+1_rkx ,  1.357e+1_rkx , -4.190e-1_rkx , -1.087e+0_rkx , &
       -2.070e+2_rkx ,  1.478e+1_rkx ,  1.557e+1_rkx , -4.353e-1_rkx , &
       -1.083e+0_rkx , -2.382e+2_rkx ,  1.486e+1_rkx ,  1.758e+1_rkx , &
       -4.389e-1_rkx , -1.078e+0_rkx , -2.716e+2_rkx ,  1.505e+1_rkx , &
        1.597e+1_rkx , -4.337e-1_rkx , -1.073e+0_rkx , -2.510e+2_rkx , &
        1.527e+1_rkx ,  2.107e+1_rkx , -3.041e-1_rkx , -1.067e+0_rkx , &
       -2.494e+2_rkx ,  1.166e+1_rkx , -2.424e-1_rkx , -1.770e-1_rkx , &
       -1.032e+0_rkx ,  1.469e-1_rkx ,  1.947e+0_rkx ,  2.535e+1_rkx , &
       -2.270e-1_rkx , -1.052e+0_rkx , -2.528e+2_rkx ,  9.888e+0_rkx , &
       -1.545e-1_rkx , -1.661e-1_rkx , -1.030e+0_rkx , -4.698e-4_rkx , &
        7.275e-2_rkx ,  8.835e-1_rkx , -1.590e-1_rkx , -1.029e+0_rkx , &
       -2.838e+1_rkx ,  2.734e+1_rkx/

  data wsbase / 7.371e-8_rkx , 7.371e-8_rkx , 7.371e-8_rkx , 7.371e-8_rkx , &
                7.371e-8_rkx , 7.371e-8_rkx , 7.371e-8_rkx , 6.583e-8_rkx , &
                6.583e-8_rkx , 3.656e-6_rkx , 4.919e-5_rkx , 3.539e-3_rkx , &
                2.855e-2_rkx , 2.126e-1_rkx , 8.433e-1_rkx , 9.653e-1_rkx , &
                6.198e-1_rkx , 9.642e-1_rkx , 9.699e-1_rkx/

  data ((wscoef(ii,jj),jj=1,ncoefs),ii=1,nspi)   /      2.492e+0_rkx , &
       -5.210e-2_rkx , -1.036e+0_rkx , -4.398e+1_rkx ,  1.724e+1_rkx , &
        2.492e+0_rkx , -5.210e-2_rkx , -1.036e+0_rkx , -4.398e+1_rkx , &
        1.724e+1_rkx ,  2.492e+0_rkx , -5.210e-2_rkx , -1.036e+0_rkx , &
       -4.398e+1_rkx ,  1.724e+1_rkx ,  2.492e+0_rkx , -5.210e-2_rkx , &
       -1.036e+0_rkx , -4.398e+1_rkx ,  1.724e+1_rkx ,  2.492e+0_rkx , &
       -5.210e-2_rkx , -1.036e+0_rkx , -4.398e+1_rkx ,  1.724e+1_rkx , &
        2.492e+0_rkx , -5.210e-2_rkx , -1.036e+0_rkx , -4.398e+1_rkx , &
        1.724e+1_rkx ,  2.492e+0_rkx , -5.210e-2_rkx , -1.036e+0_rkx , &
       -4.398e+1_rkx ,  1.724e+1_rkx ,  1.139e+0_rkx , -1.110e-2_rkx , &
       -1.011e+0_rkx , -7.754e+0_rkx ,  6.737e+0_rkx ,  1.139e+0_rkx , &
       -1.110e-2_rkx , -1.011e+0_rkx , -7.754e+0_rkx ,  6.737e+0_rkx , &
        1.848e+0_rkx , -3.920e-4_rkx , -9.924e-1_rkx , -1.607e+0_rkx , &
        8.587e-1_rkx ,  5.459e+0_rkx ,  9.357e-1_rkx , -1.626e+0_rkx , &
       -5.282e+0_rkx ,  1.066e+0_rkx ,  1.187e+0_rkx ,  2.241e-1_rkx , &
       -1.226e+0_rkx ,  1.442e+1_rkx , -1.402e+1_rkx , -3.640e+0_rkx , &
        2.552e-1_rkx , -1.168e+0_rkx ,  4.458e+1_rkx ,  1.152e+1_rkx , &
       -5.634e+0_rkx ,  2.068e-1_rkx , -1.122e+0_rkx ,  7.528e+1_rkx , &
        1.290e+1_rkx ,  1.826e-1_rkx ,  6.588e-2_rkx , -1.098e+0_rkx , &
       -1.996e-2_rkx ,  1.618e-1_rkx ,  2.164e+0_rkx ,  1.194e-1_rkx , &
       -1.044e+0_rkx , -3.221e+1_rkx ,  1.564e+1_rkx ,  2.268e-1_rkx , &
        3.266e-2_rkx , -1.064e+0_rkx , -2.677e-2_rkx ,  1.309e-1_rkx , &
        2.178e+0_rkx ,  1.151e-1_rkx , -1.042e+0_rkx , -3.325e+1_rkx , &
        1.600e+1_rkx ,  1.713e+0_rkx ,  9.166e-2_rkx , -1.039e+0_rkx , &
       -2.660e+1_rkx ,  1.629e+1_rkx/

  data gsbase / 6.899e-1_rkx , 6.899e-1_rkx , 6.899e-1_rkx , 6.899e-1_rkx , &
                6.899e-1_rkx , 6.899e-1_rkx , 6.899e-1_rkx , 6.632e-1_rkx , &
                6.632e-1_rkx , 5.912e-1_rkx , 5.111e-1_rkx , 4.269e-1_rkx , &
                3.321e-1_rkx , 2.197e-1_rkx , 1.305e-1_rkx , 7.356e-2_rkx , &
                1.602e-1_rkx , 6.883e-2_rkx , 6.304e-2_rkx /

  data ((gscoef(ii,jj),jj=1,ncoefs),ii=1,nspi) /       -9.874e-1_rkx , &
       -3.033e+1_rkx , -2.138e+1_rkx , -2.265e+0_rkx ,  5.238e+0_rkx , &
       -9.874e-1_rkx , -3.033e+1_rkx , -2.138e+1_rkx , -2.265e+0_rkx , &
        5.238e+0_rkx , -9.874e-1_rkx , -3.033e+1_rkx , -2.138e+1_rkx , &
       -2.265e+0_rkx ,  5.238e+0_rkx , -9.874e-1_rkx , -3.033e+1_rkx , &
       -2.138e+1_rkx , -2.265e+0_rkx ,  5.238e+0_rkx , -9.874e-1_rkx , &
       -3.033e+1_rkx , -2.138e+1_rkx , -2.265e+0_rkx ,  5.238e+0_rkx , &
       -9.874e-1_rkx , -3.033e+1_rkx , -2.138e+1_rkx , -2.265e+0_rkx , &
        5.238e+0_rkx , -9.874e-1_rkx , -3.033e+1_rkx , -2.138e+1_rkx , &
       -2.265e+0_rkx ,  5.238e+0_rkx , -3.666e-1_rkx , -1.319e+0_rkx , &
       -3.311e+0_rkx , -2.821e-2_rkx ,  8.844e-1_rkx , -3.666e-1_rkx , &
       -1.319e+0_rkx , -3.311e+0_rkx , -2.821e-2_rkx ,  8.844e-1_rkx , &
        5.824e-1_rkx , -1.875e-1_rkx , -1.567e+0_rkx , -4.402e+0_rkx , &
        6.268e+0_rkx ,  1.238e+0_rkx , -1.550e-1_rkx , -1.368e+0_rkx , &
       -1.127e+1_rkx ,  8.334e+0_rkx ,  2.299e+0_rkx , -1.686e-1_rkx , &
       -1.304e+0_rkx , -2.677e+1_rkx ,  1.101e+1_rkx ,  3.037e+0_rkx , &
       -1.447e-1_rkx , -1.223e+0_rkx , -2.609e+1_rkx ,  8.267e+0_rkx , &
        4.683e+0_rkx , -2.307e-1_rkx , -1.241e+0_rkx , -4.312e+1_rkx , &
        8.838e+0_rkx ,  3.842e+0_rkx , -6.301e-1_rkx , -1.367e+0_rkx , &
       -4.144e+1_rkx ,  9.620e+0_rkx ,  3.237e+0_rkx , -4.530e-1_rkx , &
       -1.204e+0_rkx , -3.234e+1_rkx ,  8.946e+0_rkx ,  4.181e+0_rkx , &
       -4.140e-1_rkx , -1.284e+0_rkx , -4.489e+1_rkx ,  9.950e+0_rkx , &
        3.378e+0_rkx , -4.334e-1_rkx , -1.188e+0_rkx , -3.664e+1_rkx , &
        9.786e+0_rkx ,  3.943e+0_rkx , -3.952e-1_rkx , -1.170e+0_rkx , &
       -4.415e+1_rkx ,  1.031e+1_rkx/

  data ksbc_hb_stand /20.7830_rkx , 17.2120_rkx , 15.8640_rkx , 15.0530_rkx , &
                      14.3040_rkx , 13.6130_rkx , 11.9660_rkx ,  6.5782_rkx , &
                       4.3961_rkx ,  4.3800_rkx ,  2.1100_rkx ,  2.1100_rkx , &
                       2.1100_rkx ,  2.1100_rkx ,  1.4000_rkx ,  1.4000_rkx , &
                       1.4000_rkx ,  1.4000_rkx ,  1.4000_rkx/

  data wsbc_hb_stand &
     / 0.245240_rkx , 0.209620_rkx , 0.195000_rkx , 0.185860_rkx , &
       0.177190_rkx , 0.168970_rkx , 0.148490_rkx , 0.071748_rkx , &
       0.037536_rkx , 0.089000_rkx , 0.025000_rkx , 0.025000_rkx , &
       0.025000_rkx , 0.025000_rkx , 0.009000_rkx , 0.009000_rkx , &
       0.009000_rkx , 0.009000_rkx , 0.009000_rkx/

  data gsbc_hb_stand &
     / 0.213870_rkx , 0.171540_rkx , 0.155640_rkx , 0.146100_rkx , &
       0.137320_rkx , 0.129230_rkx , 0.110060_rkx , 0.049175_rkx , &
       0.026638_rkx , 0.220000_rkx , 0.123000_rkx , 0.123000_rkx , &
       0.123000_rkx , 0.123000_rkx , 0.073000_rkx , 0.073000_rkx , &
       0.073000_rkx , 0.073000_rkx , 0.073000_rkx/

  data ksbc_hl_stand /14.8510_rkx , 14.2580_rkx , 13.9430_rkx , 13.7240_rkx , &
                      13.5070_rkx , 13.2950_rkx , 12.7220_rkx ,  9.4434_rkx , &
                       6.9653_rkx ,  4.3800_rkx ,  2.1100_rkx ,  2.1100_rkx , &
                       2.1100_rkx ,  2.1100_rkx ,  1.4000_rkx ,  1.4000_rkx , &
                       1.4000_rkx ,  1.4000_rkx ,  1.4000_rkx/

  data wsbc_hl_stand /0.46081_rkx , 0.44933_rkx , 0.44397_rkx , 0.44065_rkx , &
                      0.43737_rkx , 0.43394_rkx , 0.42346_rkx , 0.35913_rkx , &
                      0.29579_rkx , 0.08900_rkx , 0.02500_rkx , 0.02500_rkx , &
                      0.02500_rkx , 0.02500_rkx , 0.00900_rkx , 0.00900_rkx , &
                      0.00900_rkx , 0.00900_rkx , 0.00900_rkx/

  data gsbc_hl_stand /0.69038_rkx , 0.65449_rkx , 0.63711_rkx , 0.62542_rkx , &
                      0.61407_rkx , 0.60319_rkx , 0.57467_rkx , 0.42050_rkx , &
                      0.30660_rkx , 0.22000_rkx , 0.12300_rkx , 0.12300_rkx , &
                      0.12300_rkx , 0.12300_rkx , 0.07300_rkx , 0.07300_rkx , &
                      0.07300_rkx , 0.07300_rkx , 0.07300_rkx/

  data ksoc_hb_stand / &
    6.0584e+0_rkx , 6.0654e+0_rkx , 6.1179e+0_rkx , 6.0102e+0_rkx , &
    5.8000e+0_rkx , 5.6957e+0_rkx , 5.6494e+0_rkx , 4.3283e+0_rkx , &
    3.1485e+0_rkx , 1.3020e+0_rkx , 5.9920e-1_rkx , 2.9480e-1_rkx , &
    1.4750e-1_rkx , 7.3870e-2_rkx , 1.6830e-1_rkx , 2.6550e-1_rkx , &
    5.7700e-2_rkx , 2.2900e-1_rkx , 2.2700e-1_rkx/

  data wsoc_hb_stand /0.91735_rkx , 0.92365_rkx , 0.92941_rkx , 0.93067_rkx , &
                      0.93311_rkx , 0.93766_rkx , 0.94042_rkx , 0.95343_rkx , &
                      0.95480_rkx , 0.70566_rkx , 0.70566_rkx , 0.70566_rkx , &
                      0.70566_rkx , 0.70566_rkx , 0.70566_rkx , 0.70566_rkx , &
                      0.70566_rkx , 0.70566_rkx , 0.70566_rkx/

  data gsoc_hb_stand / &
    0.67489e+0_rkx , 0.67003e+0_rkx , 0.67725e+0_rkx , 0.65487e+0_rkx , &
    0.65117e+0_rkx , 0.66116e+0_rkx , 0.64547e+0_rkx , 0.60033e+0_rkx , &
    0.55389e+0_rkx , 5.91200e-1_rkx , 5.11100e-1_rkx , 4.26900e-1_rkx , &
    3.32100e-1_rkx , 2.19700e-1_rkx , 1.30500e-1_rkx , 7.35600e-2_rkx , &
    1.60200e-1_rkx , 6.88300e-2_rkx , 6.30400e-2_rkx/

  data ksoc_hl_stand / &
    3.5430e+0_rkx , 3.6230e+0_rkx , 3.7155e+0_rkx , 3.7120e+0_rkx , &
    3.6451e+0_rkx , 3.6376e+0_rkx , 3.6844e+0_rkx , 3.4588e+0_rkx , &
    2.9846e+0_rkx , 1.3020e+0_rkx , 5.9920e-1_rkx , 2.9480e-1_rkx , &
    1.4750e-1_rkx , 7.3870e-2_rkx , 1.6830e-1_rkx , 2.6550e-1_rkx , &
    5.7700e-2_rkx , 2.2900e-1_rkx , 2.2700e-1_rkx/

  data wsoc_hl_stand  /0.87931_rkx , 0.88292_rkx , 0.89214_rkx , &
                       0.89631_rkx , 0.89996_rkx , 0.90540_rkx , &
                       0.90805_rkx , 0.93423_rkx , 0.95012_rkx , &
                       0.85546_rkx , 0.85546_rkx , 0.85546_rkx , &
                       0.85546_rkx , 0.85546_rkx , 0.85546_rkx , &
                       0.85546_rkx , 0.85546_rkx , 0.85546_rkx , &
                       0.85546_rkx/

  data gsoc_hl_stand  / &
    0.73126e+0_rkx , 0.71089e+0_rkx , 0.72042e+0_rkx , 0.69924e+0_rkx , &
    0.69908e+0_rkx , 0.70696e+0_rkx , 0.68479e+0_rkx , 0.64879e+0_rkx , &
    0.63433e+0_rkx , 5.91200e-1_rkx , 5.11100e-1_rkx , 4.26900e-1_rkx , &
    3.32100e-1_rkx , 2.19700e-1_rkx , 1.30500e-1_rkx , 7.35600e-2_rkx , &
    1.60200e-1_rkx , 6.88300e-2_rkx , 6.30400e-2_rkx/
  !
  ! DUST OP data base for external mixing : maximum of 4 bin for the
  ! momeent , determined from Zender et al.
  !
  data ((ksdust_stand (ii,jj),jj=1,4),ii=1,nspi)/ 1.88010_rkx , 0.76017_rkx , &
        0.36681_rkx , 0.16933_rkx , 2.02540_rkx , 0.78378_rkx , 0.36845_rkx , &
        0.17002_rkx , 1.95470_rkx , 0.76389_rkx , 0.37119_rkx , 0.17032_rkx , &
        1.89960_rkx , 0.74916_rkx , 0.37220_rkx , 0.17052_rkx , 1.79460_rkx , &
        0.74044_rkx , 0.36984_rkx , 0.17072_rkx , 1.71490_rkx , 0.75401_rkx , &
        0.36892_rkx , 0.17090_rkx , 1.54310_rkx , 0.82322_rkx , 0.37505_rkx , &
        0.17143_rkx , 2.44820_rkx , 0.85680_rkx , 0.38078_rkx , 0.17396_rkx , &
        3.10670_rkx , 0.74488_rkx , 0.43688_rkx , 0.18104_rkx , 0.50391_rkx , &
        0.67245_rkx , 0.53605_rkx , 0.20599_rkx , 0.50391_rkx , 0.67245_rkx , &
        0.53605_rkx , 0.20599_rkx , 0.50391_rkx , 0.67245_rkx , 0.53605_rkx , &
        0.20599_rkx , 0.50391_rkx , 0.67245_rkx , 0.53605_rkx , 0.20599_rkx , &
        0.50391_rkx , 0.67245_rkx , 0.53605_rkx , 0.20599_rkx , 0.50391_rkx , &
        0.67245_rkx , 0.53605_rkx , 0.20599_rkx , 0.50391_rkx , 0.67245_rkx , &
        0.53605_rkx , 0.20599_rkx , 0.50391_rkx , 0.67245_rkx , 0.53605_rkx , &
        0.20599_rkx , 0.50391_rkx , 0.67245_rkx , 0.53605_rkx , 0.20599_rkx , &
        0.50391_rkx , 0.67245_rkx , 0.53605_rkx , 0.20599_rkx/

  data ((wsdust_stand (ii,jj),jj=1,4),ii=1,nspi)/ 0.64328_rkx , 0.55196_rkx , &
        0.53748_rkx , 0.54342_rkx , 0.67757_rkx , 0.56909_rkx , 0.53639_rkx , &
        0.54232_rkx , 0.67316_rkx , 0.56027_rkx , 0.53875_rkx , 0.54181_rkx , &
        0.66245_rkx , 0.55338_rkx , 0.53947_rkx , 0.54149_rkx , 0.68132_rkx , &
        0.57440_rkx , 0.54328_rkx , 0.54143_rkx , 0.67960_rkx , 0.58467_rkx , &
        0.54242_rkx , 0.54113_rkx , 0.72679_rkx , 0.68744_rkx , 0.58564_rkx , &
        0.54576_rkx , 0.94730_rkx , 0.88181_rkx , 0.80761_rkx , 0.70455_rkx , &
        0.97536_rkx , 0.89161_rkx , 0.86378_rkx , 0.75800_rkx , 0.89568_rkx , &
        0.96322_rkx , 0.95008_rkx , 0.89293_rkx , 0.89568_rkx , 0.96322_rkx , &
        0.95008_rkx , 0.89293_rkx , 0.89568_rkx , 0.96322_rkx , 0.95008_rkx , &
        0.89293_rkx , 0.89568_rkx , 0.96322_rkx , 0.95008_rkx , 0.89293_rkx , &
        0.89568_rkx , 0.96322_rkx , 0.95008_rkx , 0.89293_rkx , 0.89568_rkx , &
        0.96322_rkx , 0.95008_rkx , 0.89293_rkx , 0.89568_rkx , 0.96322_rkx , &
        0.95008_rkx , 0.89293_rkx , 0.89568_rkx , 0.96322_rkx , 0.95008_rkx , &
        0.89293_rkx , 0.89568_rkx , 0.96322_rkx , 0.95008_rkx , 0.89293_rkx , &
        0.89568_rkx , 0.96322_rkx , 0.95008_rkx , 0.89293_rkx/

  data ((gsdust_stand (ii,jj),jj=1,4),ii=1,nspi)/ 0.87114_rkx , 0.92556_rkx , &
        0.94542_rkx , 0.94831_rkx , 0.86127_rkx , 0.92100_rkx , 0.94355_rkx , &
        0.94813_rkx , 0.83800_rkx , 0.91194_rkx , 0.94304_rkx , 0.94803_rkx , &
        0.81760_rkx , 0.90442_rkx , 0.94239_rkx , 0.94796_rkx , 0.77088_rkx , &
        0.88517_rkx , 0.93710_rkx , 0.94775_rkx , 0.73925_rkx , 0.88364_rkx , &
        0.93548_rkx , 0.94763_rkx , 0.60695_rkx , 0.86086_rkx , 0.91824_rkx , &
        0.94473_rkx , 0.64393_rkx , 0.76457_rkx , 0.81330_rkx , 0.87784_rkx , &
        0.74760_rkx , 0.62041_rkx , 0.80665_rkx , 0.85974_rkx , 0.26761_rkx , &
        0.56045_rkx , 0.68897_rkx , 0.68174_rkx , 0.26761_rkx , 0.56045_rkx , &
        0.68897_rkx , 0.68174_rkx , 0.26761_rkx , 0.56045_rkx , 0.68897_rkx , &
        0.68174_rkx , 0.26761_rkx , 0.56045_rkx , 0.68897_rkx , 0.68174_rkx , &
        0.26761_rkx , 0.56045_rkx , 0.68897_rkx , 0.68174_rkx , 0.26761_rkx , &
        0.56045_rkx , 0.68897_rkx , 0.68174_rkx , 0.26761_rkx , 0.56045_rkx , &
        0.68897_rkx , 0.68174_rkx , 0.26761_rkx , 0.56045_rkx , 0.68897_rkx , &
        0.68174_rkx , 0.26761_rkx , 0.56045_rkx , 0.68897_rkx , 0.68174_rkx , &
        0.26761_rkx , 0.56045_rkx , 0.68897_rkx , 0.68174_rkx/

  data ((ksdust12_stand (ii,jj),jj=1,12),ii=1,nspi) / &
      5.82605e+0_rkx, 5.30982e+0_rkx, 1.44072e+0_rkx, 6.46088e-1_rkx, &
      4.06936e-1_rkx, 2.89033e-1_rkx, 2.33752e-1_rkx, 1.84485e-1_rkx, &
      1.25724e-1_rkx, 7.30034e-2_rkx, 1.0000e-20_rkx, 1.0000e-20_rkx, &
      3.92185e+0_rkx, 5.44157e+0_rkx, 1.42496e+0_rkx, 6.51380e-1_rkx, &
      4.09548e-1_rkx, 2.90722e-1_rkx, 2.34956e-1_rkx, 1.85304e-1_rkx, &
      1.26168e-1_rkx, 7.31844e-2_rkx, 1.0000e-20_rkx, 1.0000e-20_rkx, &
      3.34761e+0_rkx, 5.44018e+0_rkx, 1.41674e+0_rkx, 6.50973e-1_rkx, &
      4.10955e-1_rkx, 2.91565e-1_rkx, 2.35470e-1_rkx, 1.85664e-1_rkx, &
      1.26366e-1_rkx, 7.32652e-2_rkx, 1.0000e-20_rkx, 1.0000e-20_rkx, &
      3.02282e+0_rkx, 5.43838e+0_rkx, 1.41367e+0_rkx, 6.52048e-1_rkx, &
      4.12576e-1_rkx, 2.91860e-1_rkx, 2.35824e-1_rkx, 1.85914e-1_rkx, &
      1.26495e-1_rkx, 7.33181e-2_rkx, 1.0000e-20_rkx, 1.0000e-20_rkx, &
      2.73157e+0_rkx, 5.39746e+0_rkx, 1.41976e+0_rkx, 6.55378e-1_rkx, &
      4.13786e-1_rkx, 2.92314e-1_rkx, 2.36243e-1_rkx, 1.86151e-1_rkx, &
      1.26625e-1_rkx, 7.33705e-2_rkx, 4.34019e-2_rkx, 1.0000e-20_rkx, &
      2.46940e+0_rkx, 5.37280e+0_rkx, 1.42070e+0_rkx, 6.60716e-1_rkx, &
      4.14107e-1_rkx, 2.93162e-1_rkx, 2.36530e-1_rkx, 1.86366e-1_rkx, &
      1.26751e-1_rkx, 7.34223e-2_rkx, 4.34237e-2_rkx, 1.0000e-20_rkx, &
      1.88028e+0_rkx, 5.27157e+0_rkx, 1.44707e+0_rkx, 6.72827e-1_rkx, &
      4.15358e-1_rkx, 2.94211e-1_rkx, 2.37435e-1_rkx, 1.87021e-1_rkx, &
      1.27096e-1_rkx, 7.35618e-2_rkx, 4.34826e-2_rkx, 1.0000e-20_rkx, &
      5.47130e-1_rkx, 3.57883e+0_rkx, 1.81699e+0_rkx, 7.07262e-1_rkx, &
      4.28308e-1_rkx, 3.02019e-1_rkx, 2.42130e-1_rkx, 1.90466e-1_rkx, &
      1.29002e-1_rkx, 7.43120e-2_rkx, 4.38195e-2_rkx, 1.0000e-20_rkx, &
      1.13816e-1_rkx, 1.95724e+0_rkx, 2.10074e+0_rkx, 7.13763e-1_rkx, &
      4.29890e-1_rkx, 2.97812e-1_rkx, 2.48246e-1_rkx, 1.95955e-1_rkx, &
      1.30571e-1_rkx, 7.50355e-2_rkx, 4.41201e-2_rkx, 2.82240e-2_rkx, &
      1.43850e-2_rkx, 2.38501e-1_rkx, 5.81777e-1_rkx, 6.12041e-1_rkx, &
      5.34167e-1_rkx, 4.29108e-1_rkx, 3.43648e-1_rkx, 2.49103e-1_rkx, &
      1.44168e-1_rkx, 8.34307e-2_rkx, 4.71212e-2_rkx, 2.90008e-2_rkx, &
      1.43301e-2_rkx, 2.37798e-1_rkx, 5.81505e-1_rkx, 6.11995e-1_rkx, &
      5.34208e-1_rkx, 4.29167e-1_rkx, 3.43706e-1_rkx, 2.49093e-1_rkx, &
      1.44179e-1_rkx, 8.34394e-2_rkx, 4.71195e-2_rkx, 2.90007e-2_rkx, &
      1.43301e-2_rkx, 2.37798e-1_rkx, 5.81505e-1_rkx, 6.11995e-1_rkx, &
      5.34208e-1_rkx, 4.29167e-1_rkx, 3.43706e-1_rkx, 2.49093e-1_rkx, &
      1.44179e-1_rkx, 8.34394e-2_rkx, 4.71195e-2_rkx, 2.90007e-2_rkx, &
      1.43301e-2_rkx, 2.37798e-1_rkx, 5.81505e-1_rkx, 6.11995e-1_rkx, &
      5.34208e-1_rkx, 4.29167e-1_rkx, 3.43706e-1_rkx, 2.49093e-1_rkx, &
      1.44179e-1_rkx, 8.34394e-2_rkx, 4.71195e-2_rkx, 2.90007e-2_rkx, &
      1.43301e-2_rkx, 2.37798e-1_rkx, 5.81505e-1_rkx, 6.11995e-1_rkx, &
      5.34208e-1_rkx, 4.29167e-1_rkx, 3.43706e-1_rkx, 2.49093e-1_rkx, &
      1.44179e-1_rkx, 8.34394e-2_rkx, 4.71195e-2_rkx, 2.90007e-2_rkx, &
      1.42754e-2_rkx, 2.37098e-1_rkx, 5.81236e-1_rkx, 6.11932e-1_rkx, &
      5.34250e-1_rkx, 4.29172e-1_rkx, 3.43731e-1_rkx, 2.49061e-1_rkx, &
      1.44196e-1_rkx, 8.34493e-2_rkx, 4.71169e-2_rkx, 2.90011e-2_rkx, &
      1.42754e-2_rkx, 2.37098e-1_rkx, 5.81236e-1_rkx, 6.11932e-1_rkx, &
      5.34250e-1_rkx, 4.29172e-1_rkx, 3.43731e-1_rkx, 2.49061e-1_rkx, &
      1.44196e-1_rkx, 8.34493e-2_rkx, 4.71169e-2_rkx, 2.90011e-2_rkx, &
      3.61802e-3_rkx, 1.96700e-2_rkx, 2.40837e-1_rkx, 6.90188e-1_rkx, &
      7.42462e-1_rkx, 5.11023e-1_rkx, 2.91703e-1_rkx, 1.73126e-1_rkx, &
      1.56436e-1_rkx, 8.23404e-2_rkx, 4.55453e-2_rkx, 2.95979e-2_rkx, &
      1.78064e-3_rkx, 4.31951e-3_rkx, 5.18355e-2_rkx, 2.20056e-1_rkx, &
      4.25718e-1_rkx, 4.93497e-1_rkx, 4.62575e-1_rkx, 3.40392e-1_rkx, &
      1.43665e-1_rkx, 9.15974e-2_rkx, 4.41742e-2_rkx, 3.00372e-2_rkx, &
      1.78064e-3_rkx, 4.31951e-3_rkx, 5.18355e-2_rkx, 2.20056e-1_rkx, &
      4.25718e-1_rkx, 4.93497e-1_rkx, 4.62575e-1_rkx, 3.40392e-1_rkx, &
      1.43665e-1_rkx, 9.15974e-2_rkx, 4.41742e-2_rkx, 3.00372e-2_rkx /

  data ((wsdust12_stand (ii,jj),jj=1,12),ii=1,nspi) / &
      8.03780e-1_rkx, 7.36897e-1_rkx, 5.78261e-1_rkx, 5.36202e-1_rkx, &
      5.34536e-1_rkx, 5.37743e-1_rkx, 5.40048e-1_rkx, 5.42302e-1_rkx, &
      5.45056e-1_rkx, 5.47463e-1_rkx, 5.00000e-1_rkx, 5.00000e-1_rkx, &
      7.99301e-1_rkx, 7.94829e-1_rkx, 6.00564e-1_rkx, 5.45354e-1_rkx, &
      5.35369e-1_rkx, 5.36572e-1_rkx, 5.38641e-1_rkx, 5.41019e-1_rkx, &
      5.44150e-1_rkx, 5.46961e-1_rkx, 5.00000e-1_rkx, 5.00000e-1_rkx, &
      7.97433e-1_rkx, 8.20288e-1_rkx, 6.13124e-1_rkx, 5.50490e-1_rkx, &
      5.37324e-1_rkx, 5.36678e-1_rkx, 5.38149e-1_rkx, 5.40465e-1_rkx, &
      5.43734e-1_rkx, 5.46725e-1_rkx, 5.00000e-1_rkx, 5.00000e-1_rkx, &
      7.95686e-1_rkx, 8.33319e-1_rkx, 6.23304e-1_rkx, 5.55428e-1_rkx, &
      5.40009e-1_rkx, 5.36671e-1_rkx, 5.37987e-1_rkx, 5.40171e-1_rkx, &
      5.43458e-1_rkx, 5.46566e-1_rkx, 5.00000e-1_rkx, 5.00000e-1_rkx, &
      7.93633e-1_rkx, 8.47253e-1_rkx, 6.33117e-1_rkx, 5.62435e-1_rkx, &
      5.42884e-1_rkx, 5.37174e-1_rkx, 5.38106e-1_rkx, 5.39910e-1_rkx, &
      5.43194e-1_rkx, 5.46407e-1_rkx, 5.48104e-1_rkx, 5.00000e-1_rkx, &
      7.91376e-1_rkx, 8.60061e-1_rkx, 6.43199e-1_rkx, 5.71515e-1_rkx, &
      5.45557e-1_rkx, 5.38667e-1_rkx, 5.38213e-1_rkx, 5.39690e-1_rkx, &
      5.42933e-1_rkx, 5.46247e-1_rkx, 5.48018e-1_rkx, 5.00000e-1_rkx, &
      7.91539e-1_rkx, 8.85778e-1_rkx, 6.82406e-1_rkx, 6.05722e-1_rkx, &
      5.61734e-1_rkx, 5.46270e-1_rkx, 5.42012e-1_rkx, 5.40825e-1_rkx, &
      5.42557e-1_rkx, 5.45805e-1_rkx, 5.47773e-1_rkx, 5.00000e-1_rkx, &
      8.16759e-1_rkx, 9.52989e-1_rkx, 8.76231e-1_rkx, 8.07373e-1_rkx, &
      7.59934e-1_rkx, 7.23814e-1_rkx, 7.00613e-1_rkx, 6.75834e-1_rkx, &
      6.38271e-1_rkx, 5.95887e-1_rkx, 5.68268e-1_rkx, 5.00000e-1_rkx, &
      8.33502e-1_rkx, 9.78783e-1_rkx, 9.71464e-1_rkx, 9.37056e-1_rkx, &
      9.09948e-1_rkx, 8.85800e-1_rkx, 8.70473e-1_rkx, 8.49436e-1_rkx, &
      8.00642e-1_rkx, 7.26096e-1_rkx, 6.48635e-1_rkx, 5.94727e-1_rkx, &
      2.08783e-1_rkx, 7.06520e-1_rkx, 9.57579e-1_rkx, 9.79975e-1_rkx, &
      9.73860e-1_rkx, 9.64989e-1_rkx, 9.57201e-1_rkx, 9.46153e-1_rkx, &
      9.24151e-1_rkx, 8.90923e-1_rkx, 8.46069e-1_rkx, 7.95945e-1_rkx, &
      2.08637e-1_rkx, 7.06464e-1_rkx, 9.57584e-1_rkx, 9.79988e-1_rkx, &
      9.73911e-1_rkx, 9.65026e-1_rkx, 9.57130e-1_rkx, 9.46201e-1_rkx, &
      9.24183e-1_rkx, 8.90979e-1_rkx, 8.46132e-1_rkx, 7.96001e-1_rkx, &
      2.08637e-1_rkx, 7.06464e-1_rkx, 9.57584e-1_rkx, 9.79988e-1_rkx, &
      9.73911e-1_rkx, 9.65026e-1_rkx, 9.57130e-1_rkx, 9.46201e-1_rkx, &
      9.24183e-1_rkx, 8.90979e-1_rkx, 8.46132e-1_rkx, 7.96001e-1_rkx, &
      2.08637e-1_rkx, 7.06464e-1_rkx, 9.57584e-1_rkx, 9.79988e-1_rkx, &
      9.73911e-1_rkx, 9.65026e-1_rkx, 9.57130e-1_rkx, 9.46201e-1_rkx, &
      9.24183e-1_rkx, 8.90979e-1_rkx, 8.46132e-1_rkx, 7.96001e-1_rkx, &
      2.08637e-1_rkx, 7.06464e-1_rkx, 9.57584e-1_rkx, 9.79988e-1_rkx, &
      9.73911e-1_rkx, 9.65026e-1_rkx, 9.57130e-1_rkx, 9.46201e-1_rkx, &
      9.24183e-1_rkx, 8.90979e-1_rkx, 8.46132e-1_rkx, 7.96001e-1_rkx, &
      2.08492e-1_rkx, 7.06408e-1_rkx, 9.57590e-1_rkx, 9.79998e-1_rkx, &
      9.73929e-1_rkx, 9.65149e-1_rkx, 9.57123e-1_rkx, 9.46308e-1_rkx, &
      9.24198e-1_rkx, 8.91029e-1_rkx, 8.46194e-1_rkx, 7.96056e-1_rkx, &
      2.08492e-1_rkx, 7.06408e-1_rkx, 9.57590e-1_rkx, 9.79998e-1_rkx, &
      9.73929e-1_rkx, 9.65149e-1_rkx, 9.57123e-1_rkx, 9.46308e-1_rkx, &
      9.24198e-1_rkx, 8.91029e-1_rkx, 8.46194e-1_rkx, 7.96056e-1_rkx, &
      9.14436e-2_rkx, 7.18841e-1_rkx, 9.69097e-1_rkx, 9.88435e-1_rkx, &
      9.88239e-1_rkx, 9.82017e-1_rkx, 9.67379e-1_rkx, 9.46405e-1_rkx, &
      9.41366e-1_rkx, 9.12332e-1_rkx, 8.66683e-1_rkx, 8.22320e-1_rkx, &
      2.93192e-2_rkx, 4.95389e-1_rkx, 9.28109e-1_rkx, 9.84462e-1_rkx, &
      9.90160e-1_rkx, 9.90949e-1_rkx, 9.89883e-1_rkx, 9.85447e-1_rkx, &
      9.63209e-1_rkx, 9.47218e-1_rkx, 9.14608e-1_rkx, 8.88785e-1_rkx, &
      2.93192e-2_rkx, 4.95389e-1_rkx, 9.28109e-1_rkx, 9.84462e-1_rkx, &
      9.90160e-1_rkx, 9.90949e-1_rkx, 9.89883e-1_rkx, 9.85447e-1_rkx, &
      9.63209e-1_rkx, 9.47218e-1_rkx, 9.14608e-1_rkx, 8.88785e-1_rkx /

  data ((gsdust12_stand (ii,jj),jj=1,12),ii=1,nspi) / &
      5.46546e-1_rkx, 7.19875e-1_rkx, 8.86115e-1_rkx, 9.36182e-1_rkx, &
      9.44853e-1_rkx, 9.47046e-1_rkx, 9.47651e-1_rkx, 9.48064e-1_rkx, &
      9.48455e-1_rkx, 9.48708e-1_rkx, 9.90000e-1_rkx, 9.90000e-1_rkx, &
      4.15561e-1_rkx, 7.20841e-1_rkx, 8.50893e-1_rkx, 9.27581e-1_rkx, &
      9.42077e-1_rkx, 9.46199e-1_rkx, 9.47249e-1_rkx, 9.47862e-1_rkx, &
      9.48376e-1_rkx, 9.48709e-1_rkx, 9.90000e-1_rkx, 9.90000e-1_rkx, &
      3.61926e-1_rkx, 7.26308e-1_rkx, 8.31440e-1_rkx, 9.21453e-1_rkx, &
      9.40080e-1_rkx, 9.45564e-1_rkx, 9.46937e-1_rkx, 9.47729e-1_rkx, &
      9.48332e-1_rkx, 9.48704e-1_rkx, 9.90000e-1_rkx, 9.90000e-1_rkx, &
      3.31551e-1_rkx, 7.21714e-1_rkx, 8.22061e-1_rkx, 9.17115e-1_rkx, &
      9.38567e-1_rkx, 9.44893e-1_rkx, 9.46655e-1_rkx, 9.47619e-1_rkx, &
      9.48297e-1_rkx, 9.48700e-1_rkx, 9.90000e-1_rkx, 9.90000e-1_rkx, &
      3.04818e-1_rkx, 7.16225e-1_rkx, 8.10844e-1_rkx, 9.12670e-1_rkx, &
      9.36667e-1_rkx, 9.44115e-1_rkx, 9.46317e-1_rkx, 9.47471e-1_rkx, &
      9.48258e-1_rkx, 9.48695e-1_rkx, 9.48823e-1_rkx, 9.90000e-1_rkx, &
      2.81297e-1_rkx, 7.16554e-1_rkx, 7.95068e-1_rkx, 9.08398e-1_rkx, &
      9.34233e-1_rkx, 9.43280e-1_rkx, 9.45823e-1_rkx, 9.47269e-1_rkx, &
      9.48210e-1_rkx, 9.48688e-1_rkx, 9.48829e-1_rkx, 9.90000e-1_rkx, &
      2.31473e-1_rkx, 7.18668e-1_rkx, 7.66308e-1_rkx, 8.91809e-1_rkx, &
      9.23936e-1_rkx, 9.38002e-1_rkx, 9.42902e-1_rkx, 9.45948e-1_rkx, &
      9.47904e-1_rkx, 9.48668e-1_rkx, 9.48848e-1_rkx, 9.90000e-1_rkx, &
      1.08432e-1_rkx, 6.28714e-1_rkx, 6.66552e-1_rkx, 7.91368e-1_rkx, &
      8.27524e-1_rkx, 8.57309e-1_rkx, 8.71887e-1_rkx, 8.87280e-1_rkx, &
      9.07566e-1_rkx, 9.27961e-1_rkx, 9.40308e-1_rkx, 9.90000e-1_rkx, &
      5.28125e-2_rkx, 5.06675e-1_rkx, 6.50917e-1_rkx, 6.86444e-1_rkx, &
      7.56116e-1_rkx, 7.71482e-1_rkx, 8.02731e-1_rkx, 8.19214e-1_rkx, &
      8.44326e-1_rkx, 8.79297e-1_rkx, 9.10074e-1_rkx, 9.30485e-1_rkx, &
      9.10449e-3_rkx, 1.01909e-1_rkx, 3.45008e-1_rkx, 5.86731e-1_rkx, &
      6.70752e-1_rkx, 6.91530e-1_rkx, 6.88915e-1_rkx, 6.80801e-1_rkx, &
      6.94185e-1_rkx, 7.86935e-1_rkx, 8.16865e-1_rkx, 8.42174e-1_rkx, &
      9.08684e-3_rkx, 1.01787e-1_rkx, 3.44985e-1_rkx, 5.86594e-1_rkx, &
      6.70696e-1_rkx, 6.91508e-1_rkx, 6.88807e-1_rkx, 6.80825e-1_rkx, &
      6.94119e-1_rkx, 7.86899e-1_rkx, 8.16810e-1_rkx, 8.42163e-1_rkx, &
      9.08684e-3_rkx, 1.01787e-1_rkx, 3.44985e-1_rkx, 5.86594e-1_rkx, &
      6.70696e-1_rkx, 6.91508e-1_rkx, 6.88807e-1_rkx, 6.80825e-1_rkx, &
      6.94119e-1_rkx, 7.86899e-1_rkx, 8.16810e-1_rkx, 8.42163e-1_rkx, &
      9.08684e-3_rkx, 1.01787e-1_rkx, 3.44985e-1_rkx, 5.86594e-1_rkx, &
      6.70696e-1_rkx, 6.91508e-1_rkx, 6.88807e-1_rkx, 6.80825e-1_rkx, &
      6.94119e-1_rkx, 7.86899e-1_rkx, 8.16810e-1_rkx, 8.42163e-1_rkx, &
      9.08684e-3_rkx, 1.01787e-1_rkx, 3.44985e-1_rkx, 5.86594e-1_rkx, &
      6.70696e-1_rkx, 6.91508e-1_rkx, 6.88807e-1_rkx, 6.80825e-1_rkx, &
      6.94119e-1_rkx, 7.86899e-1_rkx, 8.16810e-1_rkx, 8.42163e-1_rkx, &
      9.06926e-3_rkx, 1.01666e-1_rkx, 3.44961e-1_rkx, 5.86444e-1_rkx, &
      6.70656e-1_rkx, 6.91547e-1_rkx, 6.88747e-1_rkx, 6.80876e-1_rkx, &
      6.94035e-1_rkx, 7.86855e-1_rkx, 8.16755e-1_rkx, 8.42159e-1_rkx, &
      9.06926e-3_rkx, 1.01666e-1_rkx, 3.44961e-1_rkx, 5.86444e-1_rkx, &
      6.70656e-1_rkx, 6.91547e-1_rkx, 6.88747e-1_rkx, 6.80876e-1_rkx, &
      6.94035e-1_rkx, 7.86855e-1_rkx, 8.16755e-1_rkx, 8.42159e-1_rkx, &
      3.19155e-3_rkx, 3.86167e-2_rkx, 3.12532e-1_rkx, 6.54273e-1_rkx, &
      7.29345e-1_rkx, 6.94262e-1_rkx, 5.71069e-1_rkx, 5.11984e-1_rkx, &
      7.42725e-1_rkx, 7.68741e-1_rkx, 8.04314e-1_rkx, 8.40731e-1_rkx, &
      1.26936e-3_rkx, 1.54572e-2_rkx, 1.16778e-1_rkx, 4.51257e-1_rkx, &
      6.45941e-1_rkx, 7.22575e-1_rkx, 7.37631e-1_rkx, 6.98317e-1_rkx, &
      5.61306e-1_rkx, 7.59343e-1_rkx, 7.44671e-1_rkx, 8.01060e-1_rkx, &
      1.26936e-3_rkx, 1.54572e-2_rkx, 1.16778e-1_rkx, 4.51257e-1_rkx, &
      6.45941e-1_rkx, 7.22575e-1_rkx, 7.37631e-1_rkx, 6.98317e-1_rkx, &
      5.61306e-1_rkx, 7.59343e-1_rkx, 7.44671e-1_rkx, 8.01060e-1_rkx /

  data (((ksslt(nn,ll,mm),ll=1,2),nn=1,nwav),mm=1,nih) / &
    1.10670_rkx , 0.24350_rkx , 1.13470_rkx , 0.24490_rkx , 1.14490_rkx , &
    0.24530_rkx , 1.15360_rkx , 0.24570_rkx , 1.16140_rkx , 0.24610_rkx , &
    1.16810_rkx , 0.24640_rkx , 1.18630_rkx , 0.24730_rkx , 1.26130_rkx , &
    0.25120_rkx , 1.28830_rkx , 0.25610_rkx , 1.11740_rkx , 0.26480_rkx , &
    1.07340_rkx , 0.26760_rkx , 1.01450_rkx , 0.27070_rkx , 0.93780_rkx , &
    0.27400_rkx , 0.88700_rkx , 0.27610_rkx , 0.85710_rkx , 0.27750_rkx , &
    0.77260_rkx , 0.28170_rkx , 0.55190_rkx , 0.29620_rkx , 0.23690_rkx , &
    0.31880_rkx , 0.23690_rkx , 0.31880_rkx , 2.76050_rkx , 0.62500_rkx , &
    2.79270_rkx , 0.62740_rkx , 2.80620_rkx , 0.62820_rkx , 2.81920_rkx , &
    0.62900_rkx , 2.83250_rkx , 0.62970_rkx , 2.84560_rkx , 0.63040_rkx , &
    2.88990_rkx , 0.63220_rkx , 3.08680_rkx , 0.64020_rkx , 3.20500_rkx , &
    0.64830_rkx , 2.99350_rkx , 0.66540_rkx , 2.91440_rkx , 0.67000_rkx , &
    2.83310_rkx , 0.67460_rkx , 2.69720_rkx , 0.68040_rkx , 2.60410_rkx , &
    0.68440_rkx , 2.53200_rkx , 0.68740_rkx , 2.31690_rkx , 0.69580_rkx , &
    1.74510_rkx , 0.71490_rkx , 0.96070_rkx , 0.77950_rkx , 0.96070_rkx , &
    0.77950_rkx , 3.43580_rkx , 0.78970_rkx , 3.48660_rkx , 0.79210_rkx , &
    3.50620_rkx , 0.79300_rkx , 3.52400_rkx , 0.79370_rkx , 3.54090_rkx , &
    0.79440_rkx , 3.55630_rkx , 0.79510_rkx , 3.60340_rkx , 0.79700_rkx , &
    3.84490_rkx , 0.80600_rkx , 4.01160_rkx , 0.81490_rkx , 3.83520_rkx , &
    0.83560_rkx , 3.75150_rkx , 0.84090_rkx , 3.67460_rkx , 0.84630_rkx , &
    3.53150_rkx , 0.85240_rkx , 3.43520_rkx , 0.85700_rkx , 3.35000_rkx , &
    0.86030_rkx , 3.09100_rkx , 0.86980_rkx , 2.36570_rkx , 0.89080_rkx , &
    1.39570_rkx , 0.96980_rkx , 1.39570_rkx , 0.96980_rkx , 4.13430_rkx , &
    0.95220_rkx , 4.18210_rkx , 0.95690_rkx , 4.20170_rkx , 0.95850_rkx , &
    4.21980_rkx , 0.95980_rkx , 4.23730_rkx , 0.96080_rkx , 4.25330_rkx , &
    0.96160_rkx , 4.30420_rkx , 0.96320_rkx , 4.59390_rkx , 0.97370_rkx , &
    4.80120_rkx , 0.98300_rkx , 4.68220_rkx , 1.00670_rkx , 4.59660_rkx , &
    1.01270_rkx , 4.52830_rkx , 1.01870_rkx , 4.38560_rkx , 1.02600_rkx , &
    4.29240_rkx , 1.03140_rkx , 4.19710_rkx , 1.03510_rkx , 3.90190_rkx , &
    1.04550_rkx , 3.03000_rkx , 1.06800_rkx , 1.89030_rkx , 1.15860_rkx , &
    1.89030_rkx , 1.15860_rkx , 5.84240_rkx , 1.36500_rkx , 5.88300_rkx , &
    1.36820_rkx , 5.90140_rkx , 1.36940_rkx , 5.91970_rkx , 1.37050_rkx , &
    5.93840_rkx , 1.37140_rkx , 5.95670_rkx , 1.37230_rkx , 6.02000_rkx , &
    1.37490_rkx , 6.38780_rkx , 1.38870_rkx , 6.69990_rkx , 1.40120_rkx , &
    6.73060_rkx , 1.43000_rkx , 6.65350_rkx , 1.43720_rkx , 6.61700_rkx , &
    1.44470_rkx , 6.49700_rkx , 1.45260_rkx , 6.42870_rkx , 1.45840_rkx , &
    6.31860_rkx , 1.46280_rkx , 5.96400_rkx , 1.47540_rkx , 4.77770_rkx , &
    1.50540_rkx , 3.29640_rkx , 1.61800_rkx , 3.29640_rkx , 1.61800_rkx , &
    8.47370_rkx , 2.03340_rkx , 8.57490_rkx , 2.03540_rkx , 8.61230_rkx , &
    2.03640_rkx , 8.64390_rkx , 2.03730_rkx , 8.67140_rkx , 2.03810_rkx , &
    8.69340_rkx , 2.03880_rkx , 8.74920_rkx , 2.04110_rkx , 9.21250_rkx , &
    2.06220_rkx , 9.61140_rkx , 2.07770_rkx , 9.94990_rkx , 2.11370_rkx , &
    9.90040_rkx , 2.12310_rkx , 9.92670_rkx , 2.13290_rkx , 9.87340_rkx , &
    2.14430_rkx , 9.86950_rkx , 2.15220_rkx , 9.75570_rkx , 2.15760_rkx , &
    9.36830_rkx , 2.17230_rkx , 7.80820_rkx , 2.20780_rkx , 5.96950_rkx , &
    2.34700_rkx , 5.96950_rkx , 2.34700_rkx ,14.68840_rkx , 3.60330_rkx , &
   14.82000_rkx , 3.60900_rkx ,14.86800_rkx , 3.61130_rkx ,14.90890_rkx , &
    3.61330_rkx ,14.94520_rkx , 3.61520_rkx ,14.97530_rkx , 3.61700_rkx , &
   15.05470_rkx , 3.62240_rkx ,15.63030_rkx , 3.64850_rkx ,16.16270_rkx , &
    3.66570_rkx ,17.07870_rkx , 3.72400_rkx ,17.13010_rkx , 3.73710_rkx , &
   17.29570_rkx , 3.75280_rkx ,17.43470_rkx , 3.76940_rkx ,17.60120_rkx , &
    3.78100_rkx ,17.52590_rkx , 3.78800_rkx ,17.21790_rkx , 3.80740_rkx , &
   15.24990_rkx , 3.85640_rkx ,13.22780_rkx , 4.03860_rkx ,13.22780_rkx , &
    4.03860_rkx ,22.48880_rkx , 5.62270_rkx ,22.62770_rkx , 5.61380_rkx , &
   22.67970_rkx , 5.61090_rkx ,22.72570_rkx , 5.60950_rkx ,22.76830_rkx , &
    5.60980_rkx ,22.80590_rkx , 5.61170_rkx ,22.91520_rkx , 5.62400_rkx , &
   23.61060_rkx , 5.66380_rkx ,24.24360_rkx , 5.69990_rkx ,25.74500_rkx , &
    5.76950_rkx ,25.92050_rkx , 5.78650_rkx ,26.22020_rkx , 5.80620_rkx , &
   26.57360_rkx , 5.82770_rkx ,26.92400_rkx , 5.84240_rkx ,26.92910_rkx , &
    5.85120_rkx ,26.83940_rkx , 5.87590_rkx ,24.85840_rkx , 5.94260_rkx , &
   23.35360_rkx , 6.16670_rkx ,23.35360_rkx , 6.16670_rkx /

  data (((wsslt(nn,ll,mm),ll=1,2),nn=1,nwav),mm=1,nih) / &
    0.99970_rkx , 0.99780_rkx , 0.99980_rkx , 0.99870_rkx , 0.99980_rkx , &
    0.99900_rkx , 0.99990_rkx , 0.99920_rkx , 0.99990_rkx , 0.99940_rkx , &
    0.99990_rkx , 0.99950_rkx , 1.00000_rkx , 0.99980_rkx , 1.00000_rkx , &
    1.00000_rkx , 1.00000_rkx , 1.00000_rkx , 0.99780_rkx , 0.98760_rkx , &
    0.99700_rkx , 0.98450_rkx , 0.99480_rkx , 0.97590_rkx , 0.99260_rkx , &
    0.96710_rkx , 0.99250_rkx , 0.96520_rkx , 0.99110_rkx , 0.96080_rkx , &
    0.98690_rkx , 0.94790_rkx , 0.95600_rkx , 0.85510_rkx , 0.98970_rkx , &
    0.97820_rkx , 0.98970_rkx , 0.97820_rkx , 0.99980_rkx , 0.99930_rkx , &
    0.99990_rkx , 0.99960_rkx , 1.00000_rkx , 0.99960_rkx , 1.00000_rkx , &
    0.99970_rkx , 1.00000_rkx , 0.99980_rkx , 1.00000_rkx , 0.99980_rkx , &
    1.00000_rkx , 0.99990_rkx , 1.00000_rkx , 1.00000_rkx , 1.00000_rkx , &
    1.00000_rkx , 0.99890_rkx , 0.99100_rkx , 0.99610_rkx , 0.98380_rkx , &
    0.97580_rkx , 0.96040_rkx , 0.96060_rkx , 0.94340_rkx , 0.96860_rkx , &
    0.94900_rkx , 0.96730_rkx , 0.94590_rkx , 0.95970_rkx , 0.93470_rkx , &
    0.70640_rkx , 0.71750_rkx , 0.93370_rkx , 0.85470_rkx , 0.93370_rkx , &
    0.85470_rkx , 0.99980_rkx , 0.99940_rkx , 0.99990_rkx , 0.99970_rkx , &
    1.00000_rkx , 0.99970_rkx , 1.00000_rkx , 0.99980_rkx , 1.00000_rkx , &
    0.99990_rkx , 1.00000_rkx , 0.99990_rkx , 1.00000_rkx , 1.00000_rkx , &
    1.00000_rkx , 1.00000_rkx , 1.00000_rkx , 1.00000_rkx , 0.99900_rkx , &
    0.99100_rkx , 0.99600_rkx , 0.98320_rkx , 0.97510_rkx , 0.95980_rkx , &
    0.95960_rkx , 0.94280_rkx , 0.96790_rkx , 0.94820_rkx , 0.96670_rkx , &
    0.94510_rkx , 0.95920_rkx , 0.93370_rkx , 0.70160_rkx , 0.71770_rkx , &
    0.92940_rkx , 0.83360_rkx , 0.92940_rkx , 0.83360_rkx , 1.00000_rkx , &
    0.99940_rkx , 1.00000_rkx , 0.99970_rkx , 1.00000_rkx , 0.99970_rkx , &
    1.00000_rkx , 0.99980_rkx , 1.00000_rkx , 0.99990_rkx , 1.00000_rkx , &
    0.99990_rkx , 1.00000_rkx , 1.00000_rkx , 1.00000_rkx , 1.00000_rkx , &
    1.00000_rkx , 1.00000_rkx , 0.99900_rkx , 0.99070_rkx , 0.99590_rkx , &
    0.98250_rkx , 0.97480_rkx , 0.95920_rkx , 0.95920_rkx , 0.94210_rkx , &
    0.96760_rkx , 0.94730_rkx , 0.96640_rkx , 0.94410_rkx , 0.95910_rkx , &
    0.93230_rkx , 0.70100_rkx , 0.71720_rkx , 0.92720_rkx , 0.81720_rkx , &
    0.92720_rkx , 0.81720_rkx , 1.00000_rkx , 0.99960_rkx , 1.00000_rkx , &
    0.99970_rkx , 1.00000_rkx , 0.99980_rkx , 1.00000_rkx , 0.99980_rkx , &
    1.00000_rkx , 0.99990_rkx , 1.00000_rkx , 0.99990_rkx , 1.00000_rkx , &
    1.00000_rkx , 1.00000_rkx , 1.00000_rkx , 1.00000_rkx , 1.00000_rkx , &
    0.99900_rkx , 0.99010_rkx , 0.99580_rkx , 0.98100_rkx , 0.97440_rkx , &
    0.95790_rkx , 0.95880_rkx , 0.94050_rkx , 0.96730_rkx , 0.94520_rkx , &
    0.96620_rkx , 0.94170_rkx , 0.95920_rkx , 0.92880_rkx , 0.70350_rkx , &
    0.71500_rkx , 0.92460_rkx , 0.78760_rkx , 0.92460_rkx , 0.78760_rkx , &
    1.00000_rkx , 0.99970_rkx , 1.00000_rkx , 0.99980_rkx , 1.00000_rkx , &
    0.99980_rkx , 1.00000_rkx , 0.99990_rkx , 1.00000_rkx , 0.99990_rkx , &
    1.00000_rkx , 0.99990_rkx , 1.00000_rkx , 1.00000_rkx , 1.00000_rkx , &
    1.00000_rkx , 1.00000_rkx , 1.00000_rkx , 0.99890_rkx , 0.98890_rkx , &
    0.99550_rkx , 0.97880_rkx , 0.97400_rkx , 0.95610_rkx , 0.95860_rkx , &
    0.93790_rkx , 0.96700_rkx , 0.94180_rkx , 0.96600_rkx , 0.93770_rkx , &
    0.95930_rkx , 0.92320_rkx , 0.70870_rkx , 0.71080_rkx , 0.92240_rkx , &
    0.75570_rkx , 0.92240_rkx , 0.75570_rkx , 1.00000_rkx , 0.99980_rkx , &
    1.00000_rkx , 0.99990_rkx , 1.00000_rkx , 1.00000_rkx , 1.00000_rkx , &
    1.00000_rkx , 1.00000_rkx , 1.00000_rkx , 1.00000_rkx , 1.00000_rkx , &
    1.00000_rkx , 1.00000_rkx , 1.00000_rkx , 1.00000_rkx , 1.00000_rkx , &
    1.00000_rkx , 0.99860_rkx , 0.98660_rkx , 0.99500_rkx , 0.97530_rkx , &
    0.97330_rkx , 0.95280_rkx , 0.95800_rkx , 0.93350_rkx , 0.96630_rkx , &
    0.93590_rkx , 0.96530_rkx , 0.93080_rkx , 0.95900_rkx , 0.91320_rkx , &
    0.71670_rkx , 0.70160_rkx , 0.91790_rkx , 0.71290_rkx , 0.91790_rkx , &
    0.71290_rkx , 1.00000_rkx , 0.99980_rkx , 1.00000_rkx , 0.99990_rkx , &
    1.00000_rkx , 1.00000_rkx , 1.00000_rkx , 1.00000_rkx , 1.00000_rkx , &
    1.00000_rkx , 1.00000_rkx , 1.00000_rkx , 1.00000_rkx , 1.00000_rkx , &
    1.00000_rkx , 1.00000_rkx , 1.00000_rkx , 1.00000_rkx , 0.99830_rkx , &
    0.98460_rkx , 0.99430_rkx , 0.97210_rkx , 0.97240_rkx , 0.94980_rkx , &
    0.95720_rkx , 0.92930_rkx , 0.96540_rkx , 0.93050_rkx , 0.96440_rkx , &
    0.92440_rkx , 0.95820_rkx , 0.90420_rkx , 0.72210_rkx , 0.69280_rkx , &
    0.91180_rkx , 0.68270_rkx , 0.91180_rkx , 0.68270_rkx/

  data (((gsslt(nn,ll,mm),ll=1,2), nn=1,nwav),mm=1,nih) / &
    0.73080_rkx , 0.80840_rkx , 0.71820_rkx , 0.81040_rkx , 0.71420_rkx , &
    0.81080_rkx , 0.71090_rkx , 0.81090_rkx , 0.70820_rkx , 0.81060_rkx , &
    0.70600_rkx , 0.81000_rkx , 0.70120_rkx , 0.80630_rkx , 0.69510_rkx , &
    0.79950_rkx , 0.69590_rkx , 0.79020_rkx , 0.70140_rkx , 0.77760_rkx , &
    0.69940_rkx , 0.77390_rkx , 0.69710_rkx , 0.76920_rkx , 0.69470_rkx , &
    0.76430_rkx , 0.69600_rkx , 0.76180_rkx , 0.69530_rkx , 0.76130_rkx , &
    0.69220_rkx , 0.75980_rkx , 0.64050_rkx , 0.74180_rkx , 0.59490_rkx , &
    0.70590_rkx , 0.59490_rkx , 0.70590_rkx , 0.78530_rkx , 0.84490_rkx , &
    0.78490_rkx , 0.84630_rkx , 0.78460_rkx , 0.84680_rkx , 0.78400_rkx , &
    0.84720_rkx , 0.78310_rkx , 0.84760_rkx , 0.78200_rkx , 0.84800_rkx , &
    0.77700_rkx , 0.84930_rkx , 0.77180_rkx , 0.84890_rkx , 0.77310_rkx , &
    0.84390_rkx , 0.77860_rkx , 0.83450_rkx , 0.77770_rkx , 0.83330_rkx , &
    0.77760_rkx , 0.83590_rkx , 0.77890_rkx , 0.83640_rkx , 0.78110_rkx , &
    0.83260_rkx , 0.78220_rkx , 0.83280_rkx , 0.78470_rkx , 0.83430_rkx , &
    0.77590_rkx , 0.88690_rkx , 0.70530_rkx , 0.80650_rkx , 0.70530_rkx , &
    0.80650_rkx , 0.80300_rkx , 0.84880_rkx , 0.79650_rkx , 0.84790_rkx , &
    0.79440_rkx , 0.84790_rkx , 0.79260_rkx , 0.84810_rkx , 0.79120_rkx , &
    0.84850_rkx , 0.79000_rkx , 0.84900_rkx , 0.78750_rkx , 0.85180_rkx , &
    0.78010_rkx , 0.85380_rkx , 0.77990_rkx , 0.84840_rkx , 0.78650_rkx , &
    0.84180_rkx , 0.78580_rkx , 0.84050_rkx , 0.78610_rkx , 0.84320_rkx , &
    0.78810_rkx , 0.84410_rkx , 0.79060_rkx , 0.84040_rkx , 0.79200_rkx , &
    0.84070_rkx , 0.79570_rkx , 0.84230_rkx , 0.79390_rkx , 0.89370_rkx , &
    0.72450_rkx , 0.81840_rkx , 0.72450_rkx , 0.81840_rkx , 0.80440_rkx , &
    0.83590_rkx , 0.80240_rkx , 0.85000_rkx , 0.80150_rkx , 0.85400_rkx , &
    0.80070_rkx , 0.85670_rkx , 0.79980_rkx , 0.85840_rkx , 0.79900_rkx , &
    0.85900_rkx , 0.79620_rkx , 0.85590_rkx , 0.78600_rkx , 0.85570_rkx , &
    0.78510_rkx , 0.85480_rkx , 0.79040_rkx , 0.84740_rkx , 0.79010_rkx , &
    0.84670_rkx , 0.79080_rkx , 0.84990_rkx , 0.79340_rkx , 0.85060_rkx , &
    0.79610_rkx , 0.84700_rkx , 0.79780_rkx , 0.84720_rkx , 0.80260_rkx , &
    0.84880_rkx , 0.80620_rkx , 0.89840_rkx , 0.73900_rkx , 0.82710_rkx , &
    0.73900_rkx , 0.82710_rkx , 0.81710_rkx , 0.84450_rkx , 0.81550_rkx , &
    0.85040_rkx , 0.81470_rkx , 0.85230_rkx , 0.81380_rkx , 0.85370_rkx , &
    0.81290_rkx , 0.85460_rkx , 0.81200_rkx , 0.85500_rkx , 0.80870_rkx , &
    0.85420_rkx , 0.79550_rkx , 0.86100_rkx , 0.79060_rkx , 0.85680_rkx , &
    0.79590_rkx , 0.85560_rkx , 0.79590_rkx , 0.85500_rkx , 0.79740_rkx , &
    0.85830_rkx , 0.80080_rkx , 0.86020_rkx , 0.80350_rkx , 0.85730_rkx , &
    0.80560_rkx , 0.85770_rkx , 0.81190_rkx , 0.85960_rkx , 0.82680_rkx , &
    0.90580_rkx , 0.76040_rkx , 0.84160_rkx , 0.76040_rkx , 0.84160_rkx , &
    0.82960_rkx , 0.83970_rkx , 0.82610_rkx , 0.84290_rkx , 0.82470_rkx , &
    0.84430_rkx , 0.82360_rkx , 0.84560_rkx , 0.82270_rkx , 0.84680_rkx , &
    0.82200_rkx , 0.84800_rkx , 0.82050_rkx , 0.85190_rkx , 0.80440_rkx , &
    0.86460_rkx , 0.79840_rkx , 0.86330_rkx , 0.79950_rkx , 0.86230_rkx , &
    0.79980_rkx , 0.86220_rkx , 0.80170_rkx , 0.86620_rkx , 0.80580_rkx , &
    0.86830_rkx , 0.80820_rkx , 0.86600_rkx , 0.81070_rkx , 0.86650_rkx , &
    0.81830_rkx , 0.86880_rkx , 0.84470_rkx , 0.91370_rkx , 0.78000_rkx , &
    0.85830_rkx , 0.78000_rkx , 0.85830_rkx , 0.83790_rkx , 0.83920_rkx , &
    0.83800_rkx , 0.83790_rkx , 0.83780_rkx , 0.83810_rkx , 0.83740_rkx , &
    0.83880_rkx , 0.83680_rkx , 0.84020_rkx , 0.83600_rkx , 0.84210_rkx , &
    0.83240_rkx , 0.85090_rkx , 0.81740_rkx , 0.86400_rkx , 0.80800_rkx , &
    0.87030_rkx , 0.80450_rkx , 0.86700_rkx , 0.80400_rkx , 0.86850_rkx , &
    0.80560_rkx , 0.87300_rkx , 0.80970_rkx , 0.87640_rkx , 0.81120_rkx , &
    0.87540_rkx , 0.81390_rkx , 0.87650_rkx , 0.82260_rkx , 0.88010_rkx , &
    0.86370_rkx , 0.92410_rkx , 0.79820_rkx , 0.88210_rkx , 0.79820_rkx , &
    0.88210_rkx , 0.84460_rkx , 0.80750_rkx , 0.84500_rkx , 0.82000_rkx , &
    0.84500_rkx , 0.82460_rkx , 0.84480_rkx , 0.82850_rkx , 0.84450_rkx , &
    0.83200_rkx , 0.84400_rkx , 0.83500_rkx , 0.84160_rkx , 0.84330_rkx , &
    0.82800_rkx , 0.86270_rkx , 0.81760_rkx , 0.87020_rkx , 0.80900_rkx , &
    0.87030_rkx , 0.80800_rkx , 0.87240_rkx , 0.80920_rkx , 0.87730_rkx , &
    0.81270_rkx , 0.88150_rkx , 0.81290_rkx , 0.88140_rkx , 0.81550_rkx , &
    0.88280_rkx , 0.82420_rkx , 0.88730_rkx , 0.87450_rkx , 0.93080_rkx , &
    0.80830_rkx , 0.90060_rkx , 0.80830_rkx , 0.90060_rkx/

  data rhp /0.0_rkx,0.5_rkx,0.7_rkx,0.8_rkx,0.9_rkx,0.95_rkx,0.98_rkx,0.99_rkx/

  ! DATA section for optical properties relative to RRTM
  ! based on of line calculation considering the Kok et al., 2011 distribution

  data ((ksdust_rrtm (ii,jj),jj=1,4),ii=1,nbndsw) /     &
      3.24733e-2_rkx, 3.45723e-1_rkx, 6.11556e-1_rkx, 1.31353e-1_rkx, &
      7.33946e-2_rkx, 6.27518e-1_rkx, 6.35682e-1_rkx, 1.85898e-1_rkx, &
      1.42493e-1_rkx, 8.24609e-1_rkx, 5.49972e-1_rkx, 1.70663e-1_rkx, &
      2.23231e-1_rkx, 1.07062e+0_rkx, 4.39040e-1_rkx, 1.35025e-1_rkx, &
      3.54726e-1_rkx, 1.19388e+0_rkx, 3.20780e-1_rkx, 1.57492e-1_rkx, &
      6.43902e-1_rkx, 1.35678e+0_rkx, 3.17542e-1_rkx, 1.51201e-1_rkx, &
      9.62607e-1_rkx, 1.37925e+0_rkx, 4.08074e-1_rkx, 1.54931e-1_rkx, &
      1.79510e+0_rkx, 9.68465e-1_rkx, 3.85324e-1_rkx, 1.50441e-1_rkx, &
      3.07554e+0_rkx, 6.58609e-1_rkx, 3.78327e-1_rkx, 1.45725e-1_rkx, &
      3.52915e+0_rkx, 8.63619e-1_rkx, 3.39740e-1_rkx, 1.43190e-1_rkx, &
      2.75913e+0_rkx, 7.18698e-1_rkx, 3.41137e-1_rkx, 1.41600e-1_rkx, &
      1.87476e+0_rkx, 7.44859e-1_rkx, 3.34836e-1_rkx, 1.40409e-1_rkx, &
      2.12181e+0_rkx, 7.26664e-1_rkx, 3.29921e-1_rkx, 1.39311e-1_rkx, &
      3.92232e-3_rkx, 4.74096e-2_rkx, 1.85288e-1_rkx, 2.21978e-1_rkx /

  data ((wsdust_rrtm (ii,jj),jj=1,4),ii=1,nbndsw) / &
       9.07714e-1_rkx, 9.84963e-1_rkx, 9.89606e-1_rkx, 9.48398e-1_rkx, &
       9.41418e-1_rkx, 9.87864e-1_rkx, 9.86320e-1_rkx, 9.51511e-1_rkx, &
       9.59479e-1_rkx, 9.88337e-1_rkx, 9.79395e-1_rkx, 9.36367e-1_rkx, &
       9.68132e-1_rkx, 9.88170e-1_rkx, 9.68969e-1_rkx, 9.12397e-1_rkx, &
       9.74402e-1_rkx, 9.88054e-1_rkx, 9.51572e-1_rkx, 9.17815e-1_rkx, &
       9.79381e-1_rkx, 9.85916e-1_rkx, 9.37219e-1_rkx, 8.99879e-1_rkx, &
       9.81505e-1_rkx, 9.83429e-1_rkx, 9.47774e-1_rkx, 8.91498e-1_rkx, &
       9.85477e-1_rkx, 9.66435e-1_rkx, 9.26920e-1_rkx, 8.66052e-1_rkx, &
       9.87068e-1_rkx, 9.33441e-1_rkx, 9.09271e-1_rkx, 8.24030e-1_rkx, &
       9.66195e-1_rkx, 8.77245e-1_rkx, 7.87770e-1_rkx, 6.79679e-1_rkx, &
       8.68209e-1_rkx, 6.96429e-1_rkx, 6.11392e-1_rkx, 5.51325e-1_rkx, &
       6.69915e-1_rkx, 5.91606e-1_rkx, 5.43835e-1_rkx, 5.42113e-1_rkx, &
       6.21792e-1_rkx, 5.43167e-1_rkx, 5.36133e-1_rkx, 5.44148e-1_rkx, &
       5.44753e-1_rkx, 9.35908e-1_rkx, 9.85617e-1_rkx, 9.87402e-1_rkx /

  data ((gsdust_rrtm (ii,jj),jj=1,4),ii=1,nbndsw)  /  &
      7.14045e-2_rkx, 5.57701e-1_rkx, 7.29551e-1_rkx, 5.25760e-1_rkx, &
      1.09877e-1_rkx, 6.18631e-1_rkx, 7.27320e-1_rkx, 7.39899e-1_rkx, &
      1.58598e-1_rkx, 6.59187e-1_rkx, 6.81214e-1_rkx, 7.68289e-1_rkx, &
      2.07594e-1_rkx, 7.23095e-1_rkx, 6.00556e-1_rkx, 7.16044e-1_rkx, &
      2.84977e-1_rkx, 7.21698e-1_rkx, 5.01661e-1_rkx, 7.63721e-1_rkx, &
      4.53583e-1_rkx, 7.36860e-1_rkx, 5.56412e-1_rkx, 7.65208e-1_rkx, &
      5.84279e-1_rkx, 7.25706e-1_rkx, 7.27370e-1_rkx, 7.92828e-1_rkx, &
      6.36315e-1_rkx, 6.12740e-1_rkx, 7.55553e-1_rkx, 8.09833e-1_rkx, &
      7.23353e-1_rkx, 5.47276e-1_rkx, 7.79724e-1_rkx, 8.33725e-1_rkx, &
      7.37234e-1_rkx, 7.79661e-1_rkx, 8.17743e-1_rkx, 8.88829e-1_rkx, &
      6.84414e-1_rkx, 8.21398e-1_rkx, 9.04363e-1_rkx, 9.41787e-1_rkx, &
      6.27575e-1_rkx, 8.95405e-1_rkx, 9.39012e-1_rkx, 9.47955e-1_rkx, &
      8.17597e-1_rkx, 9.30697e-1_rkx, 9.46158e-1_rkx, 9.48354e-1_rkx, &
      1.84466e-2_rkx, 1.34666e-1_rkx, 4.55520e-1_rkx, 6.80178e-1_rkx /

  data ((ksdust12_rrtm (ii,jj),jj=1,12),ii=1,nbndsw) / &
      2.72948e-3_rkx, 1.05097e-2_rkx, 1.28690e-1_rkx, 4.58140e-1_rkx, &
      6.25384e-1_rkx, 5.67730e-1_rkx, 4.20692e-1_rkx, 2.20517e-1_rkx, &
      1.46855e-1_rkx, 8.05498e-2_rkx, 4.74997e-2_rkx, 3.00359e-2_rkx, &
      3.96960e-3_rkx, 2.23419e-2_rkx, 2.57957e-1_rkx, 6.75648e-1_rkx, &
      7.47434e-1_rkx, 5.12739e-1_rkx, 2.88169e-1_rkx, 1.64365e-1_rkx, &
      1.49056e-1_rkx, 7.84490e-2_rkx, 4.71086e-2_rkx, 2.93560e-2_rkx, &
      5.58669e-3_rkx, 4.28785e-2_rkx, 4.33757e-1_rkx, 9.39034e-1_rkx, &
      7.62785e-1_rkx, 3.76847e-1_rkx, 2.09533e-1_rkx, 2.15889e-1_rkx, &
      1.16729e-1_rkx, 8.28706e-2_rkx, 4.62397e-2_rkx, 2.92904e-2_rkx, &
      7.26846e-3_rkx, 6.86183e-2_rkx, 6.46740e-1_rkx, 1.02583e+0_rkx, &
      6.83835e-1_rkx, 2.85331e-1_rkx, 2.19521e-1_rkx, 2.52894e-1_rkx, &
      1.40686e-1_rkx, 7.36198e-2_rkx, 4.62776e-2_rkx, 2.90049e-2_rkx, &
      1.00291e-2_rkx, 1.17021e-1_rkx, 9.61487e-1_rkx, 1.16843e+0_rkx, &
      5.61863e-1_rkx, 2.53167e-1_rkx, 2.86868e-1_rkx, 2.16731e-1_rkx, &
      1.33201e-1_rkx, 7.83249e-2_rkx, 4.55693e-2_rkx, 2.89160e-2_rkx, &
      1.69911e-2_rkx, 2.47432e-1_rkx, 1.32737e+0_rkx, 1.13745e+0_rkx, &
      3.72238e-1_rkx, 3.47236e-1_rkx, 2.88531e-1_rkx, 1.89873e-1_rkx, &
      1.29320e-1_rkx, 7.61732e-2_rkx, 4.52760e-2_rkx, 2.88026e-2_rkx, &
      2.49299e-2_rkx, 3.93107e-1_rkx, 1.74189e+0_rkx, 1.05356e+0_rkx, &
      3.33502e-1_rkx, 3.96214e-1_rkx, 2.23144e-1_rkx, 2.26244e-1_rkx, &
      1.18460e-1_rkx, 7.90690e-2_rkx, 4.42922e-2_rkx, 2.85342e-2_rkx, &
      6.09426e-2_rkx, 8.91504e-1_rkx, 2.03876e+0_rkx, 6.82357e-1_rkx, &
      4.62954e-1_rkx, 3.20372e-1_rkx, 2.47584e-1_rkx, 1.95196e-1_rkx, &
      1.25126e-1_rkx, 7.48942e-2_rkx, 4.45527e-2_rkx, 2.84688e-2_rkx, &
      2.01171e-1_rkx, 2.39732e+0_rkx, 1.86163e+0_rkx, 7.22278e-1_rkx, &
      3.93826e-1_rkx, 2.96730e-1_rkx, 2.49322e-1_rkx, 1.92359e-1_rkx, &
      1.22366e-1_rkx, 7.37752e-2_rkx, 4.40550e-2_rkx, 2.82308e-2_rkx, &
      6.22058e-1_rkx, 3.83522e+0_rkx, 1.17175e+0_rkx, 6.87013e-1_rkx, &
      4.27299e-1_rkx, 3.06369e-1_rkx, 2.42542e-1_rkx, 1.89225e-1_rkx, &
      1.21452e-1_rkx, 7.32786e-2_rkx, 4.37889e-2_rkx, 2.80936e-2_rkx, &
      1.77007e+0_rkx, 5.17408e+0_rkx, 1.25537e+0_rkx, 6.72982e-1_rkx, &
      4.10053e-1_rkx, 2.92894e-1_rkx, 2.36534e-1_rkx, 1.85151e-1_rkx, &
      1.19775e-1_rkx, 7.26544e-2_rkx, 4.35189e-2_rkx, 2.79676e-2_rkx, &
      3.85513e+0_rkx, 5.04019e+0_rkx, 1.30027e+0_rkx, 6.30313e-1_rkx, &
      4.01259e-1_rkx, 2.91139e-1_rkx, 2.33807e-1_rkx, 1.83256e-1_rkx, &
      1.18858e-1_rkx, 7.22256e-2_rkx, 4.33345e-2_rkx, 1.0000e-20_rkx, &
      7.61601e+0_rkx, 3.65543e+0_rkx, 1.19341e+0_rkx, 6.21575e-1_rkx, &
      3.94794e-1_rkx, 2.87308e-1_rkx, 2.31411e-1_rkx, 1.81592e-1_rkx, &
      1.18014e-1_rkx, 7.18499e-2_rkx, 4.31715e-2_rkx, 1.0000e-20_rkx, &
      9.97750e-4_rkx, 1.76373e-3_rkx, 1.43158e-2_rkx, 6.47129e-2_rkx, &
      1.47315e-1_rkx, 2.11338e-1_rkx, 2.42352e-1_rkx, 2.51536e-1_rkx, &
      1.94438e-1_rkx, 1.00233e-1_rkx, 4.90044e-2_rkx, 3.17719e-2_rkx /

  data ((wsdust12_rrtm (ii,jj),jj=1,12),ii=1,nbndsw) / &
      1.09340e-1_rkx, 7.44437e-1_rkx, 9.71476e-1_rkx, 9.87375e-1_rkx, &
      9.90051e-1_rkx, 9.88374e-1_rkx, 9.83892e-1_rkx, 9.67989e-1_rkx, &
      9.54207e-1_rkx, 9.28060e-1_rkx, 8.96477e-1_rkx, 8.58330e-1_rkx, &
      1.78006e-1_rkx, 8.32865e-1_rkx, 9.78937e-1_rkx, 9.88943e-1_rkx, &
      9.88544e-1_rkx, 9.82510e-1_rkx, 9.67867e-1_rkx, 9.45076e-1_rkx, &
      9.41136e-1_rkx, 9.10361e-1_rkx, 8.72626e-1_rkx, 8.23501e-1_rkx, &
      2.59037e-1_rkx, 8.86310e-1_rkx, 9.81759e-1_rkx, 9.88627e-1_rkx, &
      9.85287e-1_rkx, 9.69910e-1_rkx, 9.45027e-1_rkx, 9.47820e-1_rkx, &
      9.14237e-1_rkx, 8.99370e-1_rkx, 8.47013e-1_rkx, 7.93639e-1_rkx, &
      3.30266e-1_rkx, 9.14146e-1_rkx, 9.83361e-1_rkx, 9.88473e-1_rkx, &
      9.80679e-1_rkx, 9.52761e-1_rkx, 9.39744e-1_rkx, 9.49006e-1_rkx, &
      9.22952e-1_rkx, 8.73460e-1_rkx, 8.30465e-1_rkx, 7.70237e-1_rkx, &
      4.19588e-1_rkx, 9.36441e-1_rkx, 9.86080e-1_rkx, 9.87200e-1_rkx, &
      9.71644e-1_rkx, 9.38151e-1_rkx, 9.46501e-1_rkx, 9.32188e-1_rkx, &
      9.08616e-1_rkx, 8.67375e-1_rkx, 8.08301e-1_rkx, 7.46269e-1_rkx, &
      5.56481e-1_rkx, 9.58262e-1_rkx, 9.87344e-1_rkx, 9.83070e-1_rkx, &
      9.46494e-1_rkx, 9.44305e-1_rkx, 9.34754e-1_rkx, 9.11477e-1_rkx, &
      8.89347e-1_rkx, 8.40595e-1_rkx, 7.77415e-1_rkx, 7.11804e-1_rkx, &
      6.48767e-1_rkx, 9.68363e-1_rkx, 9.87373e-1_rkx, 9.76232e-1_rkx, &
      9.37328e-1_rkx, 9.42950e-1_rkx, 9.07787e-1_rkx, 9.20349e-1_rkx, &
      8.65349e-1_rkx, 8.30818e-1_rkx, 7.50775e-1_rkx, 6.86180e-1_rkx, &
      7.74388e-1_rkx, 9.76674e-1_rkx, 9.86210e-1_rkx, 9.53073e-1_rkx, &
      9.36816e-1_rkx, 9.18098e-1_rkx, 9.04123e-1_rkx, 8.87840e-1_rkx, &
      8.48355e-1_rkx, 7.88857e-1_rkx, 7.15187e-1_rkx, 6.49858e-1_rkx, &
      9.02079e-1_rkx, 9.84543e-1_rkx, 9.77098e-1_rkx, 9.41579e-1_rkx, &
      9.07323e-1_rkx, 8.90279e-1_rkx, 8.76244e-1_rkx, 8.53588e-1_rkx, &
      8.02472e-1_rkx, 7.32884e-1_rkx, 6.57070e-1_rkx, 6.01391e-1_rkx, &
      9.01246e-1_rkx, 9.72366e-1_rkx, 8.91438e-1_rkx, 8.50295e-1_rkx, &
      8.12766e-1_rkx, 7.76625e-1_rkx, 7.48010e-1_rkx, 7.16036e-1_rkx, &
      6.59625e-1_rkx, 6.05375e-1_rkx, 5.69713e-1_rkx, 5.54799e-1_rkx, &
      8.69690e-1_rkx, 9.32847e-1_rkx, 7.63443e-1_rkx, 7.01149e-1_rkx, &
      6.31386e-1_rkx, 5.92917e-1_rkx, 5.76659e-1_rkx, 5.61551e-1_rkx, &
      5.47911e-1_rkx, 5.45667e-1_rkx, 5.47216e-1_rkx, 5.48324e-1_rkx, &
      8.42400e-1_rkx, 8.58172e-1_rkx, 6.57181e-1_rkx, 5.72085e-1_rkx, &
      5.48832e-1_rkx, 5.41316e-1_rkx, 5.39378e-1_rkx, 5.40216e-1_rkx, &
      5.43356e-1_rkx, 5.46259e-1_rkx, 5.47985e-1_rkx, 5.00000e-1_rkx, &
      8.48520e-1_rkx, 7.09509e-1_rkx, 5.70206e-1_rkx, 5.35622e-1_rkx, &
      5.34532e-1_rkx, 5.37518e-1_rkx, 5.39774e-1_rkx, 5.42083e-1_rkx, &
      5.45190e-1_rkx, 5.47380e-1_rkx, 5.48590e-1_rkx, 5.00000e-1_rkx, &
      2.09253e-2_rkx, 2.91017e-1_rkx, 8.11386e-1_rkx, 9.53685e-1_rkx, &
      9.81216e-1_rkx, 9.87701e-1_rkx, 9.89537e-1_rkx, 9.90139e-1_rkx, &
      9.83821e-1_rkx, 9.69371e-1_rkx, 9.48524e-1_rkx, 9.28404e-1_rkx /

  data ((gsdust12_rrtm (ii,jj),jj=1,12),ii=1,nbndsw) / &
      3.46379e-3_rkx, 2.99321e-2_rkx, 2.02337e-1_rkx, 6.16252e-1_rkx, &
      7.21069e-1_rkx, 7.35882e-1_rkx, 7.00282e-1_rkx, 5.57231e-1_rkx, &
      6.99890e-1_rkx, 7.50332e-1_rkx, 7.94881e-1_rkx, 8.23712e-1_rkx, &
      5.32454e-3_rkx, 4.58652e-2_rkx, 3.30338e-1_rkx, 6.41233e-1_rkx, &
      7.35060e-1_rkx, 6.99022e-1_rkx, 5.71750e-1_rkx, 4.98460e-1_rkx, &
      7.67045e-1_rkx, 7.60635e-1_rkx, 8.10746e-1_rkx, 8.37959e-1_rkx, &
      7.60988e-3_rkx, 6.54008e-2_rkx, 4.95852e-1_rkx, 7.22164e-1_rkx, &
      7.24757e-1_rkx, 6.03182e-1_rkx, 4.73265e-1_rkx, 6.86063e-1_rkx, &
      7.21469e-1_rkx, 8.01505e-1_rkx, 8.24644e-1_rkx, 8.53110e-1_rkx, &
      9.78847e-3_rkx, 8.40745e-2_rkx, 6.01693e-1_rkx, 7.18575e-1_rkx, &
      6.95022e-1_rkx, 4.89171e-1_rkx, 5.58598e-1_rkx, 7.72881e-1_rkx, &
      7.75748e-1_rkx, 7.87758e-1_rkx, 8.35258e-1_rkx, 8.61751e-1_rkx, &
      1.29334e-2_rkx, 1.11347e-1_rkx, 6.19352e-1_rkx, 7.38646e-1_rkx, &
      6.30267e-1_rkx, 4.94111e-1_rkx, 7.14231e-1_rkx, 7.65081e-1_rkx, &
      7.69894e-1_rkx, 8.11382e-1_rkx, 8.44376e-1_rkx, 8.72182e-1_rkx, &
      1.92973e-2_rkx, 1.68701e-1_rkx, 6.58747e-1_rkx, 7.18152e-1_rkx, &
      4.90052e-1_rkx, 7.02517e-1_rkx, 7.71167e-1_rkx, 7.40547e-1_rkx, &
      7.91408e-1_rkx, 8.24420e-1_rkx, 8.59087e-1_rkx, 8.86589e-1_rkx, &
      2.52011e-2_rkx, 2.25066e-1_rkx, 7.26645e-1_rkx, 6.79960e-1_rkx, &
      5.29683e-1_rkx, 7.76974e-1_rkx, 7.17543e-1_rkx, 7.93729e-1_rkx, &
      7.84479e-1_rkx, 8.42833e-1_rkx, 8.66250e-1_rkx, 8.95115e-1_rkx, &
      4.19456e-2_rkx, 4.02610e-1_rkx, 7.29282e-1_rkx, 5.52694e-1_rkx, &
      7.13001e-1_rkx, 7.58849e-1_rkx, 7.60124e-1_rkx, 7.82627e-1_rkx, &
      8.17454e-1_rkx, 8.51547e-1_rkx, 8.84484e-1_rkx, 9.09941e-1_rkx, &
      8.28151e-2_rkx, 6.19089e-1_rkx, 6.85932e-1_rkx, 7.02482e-1_rkx, &
      7.39448e-1_rkx, 7.66561e-1_rkx, 8.02612e-1_rkx, 8.15096e-1_rkx, &
      8.44005e-1_rkx, 8.76250e-1_rkx, 9.06868e-1_rkx, 9.28012e-1_rkx, &
      1.47716e-1_rkx, 6.85494e-1_rkx, 5.44204e-1_rkx, 7.79828e-1_rkx, &
      8.07462e-1_rkx, 8.40093e-1_rkx, 8.53897e-1_rkx, 8.72429e-1_rkx, &
      8.99659e-1_rkx, 9.23888e-1_rkx, 9.39399e-1_rkx, 9.45923e-1_rkx, &
      2.84132e-1_rkx, 7.44738e-1_rkx, 7.33322e-1_rkx, 8.50731e-1_rkx, &
      8.92127e-1_rkx, 9.12016e-1_rkx, 9.24673e-1_rkx, 9.34774e-1_rkx, &
      9.44573e-1_rkx, 9.48195e-1_rkx, 9.48839e-1_rkx, 9.48915e-1_rkx, &
      5.04929e-1_rkx, 7.52464e-1_rkx, 8.65095e-1_rkx, 9.05884e-1_rkx, &
      9.33183e-1_rkx, 9.42013e-1_rkx, 9.45047e-1_rkx, 9.46979e-1_rkx, &
      9.48268e-1_rkx, 9.48695e-1_rkx, 9.48832e-1_rkx, 9.90000e-1_rkx, &
      6.32938e-1_rkx, 6.96838e-1_rkx, 8.93021e-1_rkx, 9.34684e-1_rkx, &
      9.44583e-1_rkx, 9.46880e-1_rkx, 9.47578e-1_rkx, 9.48038e-1_rkx, &
      9.48490e-1_rkx, 9.48714e-1_rkx, 9.48767e-1_rkx, 9.90000e-1_rkx, &
      8.87657e-4_rkx, 7.71129e-3_rkx, 5.06849e-2_rkx, 1.80875e-1_rkx, &
      3.66911e-1_rkx, 5.23451e-1_rkx, 6.19055e-1_rkx, 6.74429e-1_rkx, &
      6.73365e-1_rkx, 6.39275e-1_rkx, 6.80069e-1_rkx, 7.61202e-1_rkx /

  data ((ksdust_rrtm_lw (ii,jj),jj=1,4),ii=1,nbndlw) / &
      9.03338e-8_rkx, 1.66503e-6_rkx, 1.39793e-5_rkx, 1.23597e-4_rkx, &
      5.29696e-5_rkx, 1.02266e-3_rkx, 1.02906e-2_rkx, 9.09172e-2_rkx, &
      1.03134e-4_rkx, 1.98371e-3_rkx, 1.88747e-2_rkx, 8.67894e-2_rkx, &
      1.14509e-4_rkx, 2.17008e-3_rkx, 1.89807e-2_rkx, 7.34485e-2_rkx, &
      2.42895e-4_rkx, 4.57206e-3_rkx, 3.87656e-2_rkx, 1.20615e-1_rkx, &
      5.78706e-4_rkx, 1.13272e-2_rkx, 8.88989e-2_rkx, 1.42348e-1_rkx, &
      1.18229e-3_rkx, 2.31581e-2_rkx, 1.35854e-1_rkx, 9.96127e-2_rkx, &
      1.56680e-3_rkx, 2.95729e-2_rkx, 1.32997e-1_rkx, 9.00509e-2_rkx, &
      7.11396e-4_rkx, 1.17724e-2_rkx, 4.78602e-2_rkx, 6.83209e-2_rkx, &
      1.42405e-3_rkx, 2.46294e-2_rkx, 9.97541e-2_rkx, 1.19819e-1_rkx, &
      2.54564e-3_rkx, 4.13630e-2_rkx, 1.42638e-1_rkx, 1.25996e-1_rkx, &
      6.56414e-3_rkx, 1.00878e-1_rkx, 2.92038e-1_rkx, 1.22892e-1_rkx, &
      9.65619e-3_rkx, 1.44428e-1_rkx, 3.71668e-1_rkx, 1.06664e-1_rkx, &
      1.18215e-2_rkx, 1.63034e-1_rkx, 3.92363e-1_rkx, 9.76398e-2_rkx, &
      1.47898e-2_rkx, 1.85000e-1_rkx, 4.00920e-1_rkx, 9.13925e-2_rkx, &
      2.63245e-2_rkx, 2.62229e-1_rkx, 4.06941e-1_rkx, 7.63635e-2_rkx /

  data ((ksdust12_rrtm_lw (ii,jj),jj=1,12),ii=1,nbndlw) / &
    9.4332e-10_rkx, 2.42747e-8_rkx, 4.12358e-7_rkx, 2.49960e-6_rkx, &
    8.67040e-6_rkx, 2.00439e-5_rkx, 3.48134e-5_rkx, 6.39697e-5_rkx, &
    1.80706e-4_rkx, 3.62353e-4_rkx, 6.48472e-4_rkx, 1.00623e-3_rkx, &
    5.49405e-7_rkx, 1.41745e-5_rkx, 2.44980e-4_rkx, 1.56266e-3_rkx, &
    5.98433e-3_rkx, 1.56573e-2_rkx, 2.99434e-2_rkx, 5.76057e-2_rkx, &
    9.57461e-2_rkx, 5.27134e-2_rkx, 2.81645e-2_rkx, 1.83919e-2_rkx, &
    1.06979e-6_rkx, 2.76003e-5_rkx, 4.76715e-4_rkx, 3.02346e-3_rkx, &
    1.12907e-2_rkx, 2.76508e-2_rkx, 4.74603e-2_rkx, 7.34485e-2_rkx, &
    8.13079e-2_rkx, 5.00051e-2_rkx, 2.66968e-2_rkx, 1.72178e-2_rkx, &
    1.19006e-6_rkx, 3.06811e-5_rkx, 5.27223e-4_rkx, 3.28672e-3_rkx, &
    1.17650e-2_rkx, 2.67074e-2_rkx, 4.19379e-2_rkx, 5.92327e-2_rkx, &
    7.56272e-2_rkx, 5.56937e-2_rkx, 2.53841e-2_rkx, 1.62879e-2_rkx, &
    2.52349e-6_rkx, 6.50768e-5_rkx, 1.11702e-3_rkx, 6.90405e-3_rkx, &
    2.43645e-2_rkx, 5.34150e-2_rkx, 8.02288e-2_rkx, 1.12439e-1_rkx, &
    1.17072e-1_rkx, 5.53413e-2_rkx, 2.55754e-2_rkx, 1.65092e-2_rkx, &
    5.96784e-6_rkx, 1.54341e-4_rkx, 2.69697e-3_rkx, 1.72493e-2_rkx, &
    5.96038e-2_rkx, 1.13530e-1_rkx, 1.49212e-1_rkx, 1.61819e-1_rkx, &
    1.13926e-1_rkx, 4.20387e-2_rkx, 2.70152e-2_rkx, 1.68667e-2_rkx, &
    1.21514e-5_rkx, 3.14708e-4_rkx, 5.53147e-3_rkx, 3.49294e-2_rkx, &
    1.04488e-1_rkx, 1.52373e-1_rkx, 1.54806e-1_rkx, 1.36628e-1_rkx, &
    7.82340e-2_rkx, 4.33355e-2_rkx, 2.64115e-2_rkx, 1.70195e-2_rkx, &
    1.61361e-5_rkx, 4.17667e-4_rkx, 7.27898e-3_rkx, 4.36492e-2_rkx, &
    1.11867e-1_rkx, 1.39808e-1_rkx, 1.35289e-1_rkx, 1.18753e-1_rkx, &
    7.34595e-2_rkx, 4.31202e-2_rkx, 2.61331e-2_rkx, 1.68656e-2_rkx, &
    7.48894e-6_rkx, 1.92223e-4_rkx, 3.17285e-3_rkx, 1.67118e-2_rkx, &
    3.90610e-2_rkx, 5.42006e-2_rkx, 6.33148e-2_rkx, 6.91587e-2_rkx, &
    6.40591e-2_rkx, 4.16369e-2_rkx, 2.36985e-2_rkx, 1.53070e-2_rkx, &
    1.48517e-5_rkx, 3.82566e-4_rkx, 6.45799e-3_rkx, 3.51449e-2_rkx, &
    8.00462e-2_rkx, 1.15857e-1_rkx, 1.31083e-1_rkx, 1.35583e-1_rkx, &
    1.00870e-1_rkx, 4.19620e-2_rkx, 2.46062e-2_rkx, 1.57179e-2_rkx, &
    2.66658e-5_rkx, 6.86078e-4_rkx, 1.13876e-2_rkx, 5.71044e-2_rkx, &
    1.16591e-1_rkx, 1.56714e-1_rkx, 1.68116e-1_rkx, 1.59980e-1_rkx, &
    9.71746e-2_rkx, 3.96648e-2_rkx, 2.44858e-2_rkx, 1.57377e-2_rkx, &
    6.92718e-5_rkx, 1.77937e-3_rkx, 2.90183e-2_rkx, 1.33955e-1_rkx, &
    2.70079e-1_rkx, 3.10977e-1_rkx, 2.91058e-1_rkx, 2.25479e-1_rkx, &
    7.58939e-2_rkx, 4.71680e-2_rkx, 2.51570e-2_rkx, 1.59894e-2_rkx, &
    9.92866e-5_rkx, 2.57410e-3_rkx, 4.38245e-2_rkx, 1.86049e-1_rkx, &
    3.42407e-1_rkx, 3.74026e-1_rkx, 3.37175e-1_rkx, 2.28502e-1_rkx, &
    6.60027e-2_rkx, 4.75567e-2_rkx, 2.45834e-2_rkx, 1.60022e-2_rkx, &
    1.21841e-4_rkx, 3.15727e-3_rkx, 5.28732e-2_rkx, 2.07183e-1_rkx, &
    3.53059e-1_rkx, 3.73059e-1_rkx, 3.26846e-1_rkx, 2.15787e-1_rkx, &
    6.47597e-2_rkx, 4.54902e-2_rkx, 2.44727e-2_rkx, 1.58424e-2_rkx, &
    1.52949e-4_rkx, 3.96055e-3_rkx, 6.48010e-2_rkx, 2.34767e-1_rkx, &
    3.73278e-1_rkx, 3.82935e-1_rkx, 3.20989e-1_rkx, 2.07929e-1_rkx, &
    6.50166e-2_rkx, 4.32298e-2_rkx, 2.46475e-2_rkx, 1.57455e-2_rkx, &
    2.74053e-4_rkx, 7.09637e-3_rkx, 1.07272e-1_rkx, 3.26611e-1_rkx, &
    4.23478e-1_rkx, 3.71831e-1_rkx, 2.75473e-1_rkx, 1.50134e-1_rkx, &
    7.01898e-2_rkx, 4.08139e-2_rkx, 2.46670e-2_rkx, 1.57207e-2_rkx /

  data kssm1_rrtm /0.06_rkx, 0.084_rkx, 0.111_rkx, 0.142_rkx, &
      0.194_rkx, 0.301_rkx, 0.411_rkx, 0.934_rkx, 2.34_rkx, 4.88_rkx, &
      8.525_rkx, 12.897_rkx, 15.918_rkx, 0.0281_rkx/

  data gssm1_rrtm /0.018_rkx, 0.027_rkx, 0.039_rkx, 0.051_rkx, &
      0.068_rkx, 0.099_rkx, 0.126_rkx, 0.216_rkx, &
      0.383_rkx, 0.531_rkx, 0.629_rkx, 0.693_rkx, &
      0.727_rkx, 0.0046_rkx/

  data wssm1_rrtm /0.071_rkx, 0.123_rkx, 0.188_rkx, 0.255_rkx, &
      0.340_rkx, 0.461_rkx, 0.541_rkx, 0.671_rkx, &
      0.797_rkx, 0.848_rkx, 0.875_rkx, 0.884_rkx, &
      0.875_rkx, 0.0114_rkx/

  data kssm2_rrtm /0.045_rkx, 0.065_rkx, 0.093_rkx, 0.129_rkx, &
      0.192_rkx, 0.334_rkx, 0.481_rkx, 1.164_rkx, 2.919_rkx, 5.752_rkx, &
      9.242_rkx, 15.528_rkx, 13.755_rkx, 0.0181_rkx/

  data gssm2_rrtm /0.0258_rkx, 0.039_rkx, 0.055_rkx, 0.072_rkx, &
      0.097_rkx, 0.141_rkx, 0.179_rkx, 0.295_rkx, &
      0.479_rkx, 0.599_rkx, 0.671_rkx, 0.713_rkx, &
      0.721_rkx, 0.007_rkx/

  data wssm2_rrtm /0.167_rkx, 0.268_rkx, 0.375_rkx, 0.468_rkx, &
      0.565_rkx, 0.676_rkx, 0.738_rkx, 0.818_rkx, &
      0.887_rkx, 0.912_rkx, 0.923_rkx, 0.921_rkx, &
      0.904_rkx, 0.029_rkx/

  data aerclima_chtr     / 'BC_HB ' , 'BC_HL ' , 'OC_HB ' , 'OC_HL ' , &
                           'SO2   ' , 'SO4   ' , 'SSLT01' , 'SSLT02' , &
                           'DUST01' , 'DUST02' , 'DUST03' , 'DUST04' /

  data ncaec / -1 /

  contains

    subroutine allocate_mod_rad_aerosol
      implicit none

      if ( irrtm == 1 ) then
        nband = nbndsw
      else
        nband = nspi
      end if

      npoints = (jci2-jci1+1)*(ici2-ici1+1)

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

      if ( iclimaaer == 2 ) then
        if ( myid == iocpu ) then
          call getmem2d(alon,jcross1,jcross2,icross1,icross2,'aerosol:alon')
          call getmem2d(alat,jcross1,jcross2,icross1,icross2,'aerosol:alat')
          call getmem3d(zp3d,jcross1,jcross2, &
                             icross1,icross2,1,kz,'aerosol:zp3d')
          call getmem3d(zdz3d,jcross1,jcross2, &
                             icross1,icross2,1,kz,'aerosol:zdz3d')
          ! FAB: now define radiative OP on kth radiative levels,
          ! including strato hat
          call getmem3d(zpr3d,jcross1,jcross2, &
                              icross1,icross2,1,kth,'aerosol:zpr3d')
          call getmem3d(zdzr3d,jcross1,jcross2, &
                              icross1,icross2,1,kth,'aerosol:zdz3d')
          call getmem3d(ext1,jcross1,jcross2, &
                             icross1,icross2,1,kth,'aerosol:ext1')
          call getmem3d(ext2,jcross1,jcross2, &
                             icross1,icross2,1,kth,'aerosol:ext2')
          call getmem3d(ext,jcross1,jcross2, &
                            icross1,icross2,1,kth,'aerosol:ext')
          call getmem3d(ssa1,jcross1,jcross2, &
                             icross1,icross2,1,kth,'aerosol:ssa1')
          call getmem3d(ssa2,jcross1,jcross2, &
                             icross1,icross2,1,kth,'aerosol:ssa2')
          call getmem3d(ssa,jcross1,jcross2, &
                            icross1,icross2,1,kth,'aerosol:ssa')
          call getmem3d(asy1,jcross1,jcross2, &
                             icross1,icross2,1,kth,'aerosol:asy1')
          call getmem3d(asy2,jcross1,jcross2, &
                             icross1,icross2,1,kth,'aerosol:asy2')
          call getmem3d(asy,jcross1,jcross2, &
                            icross1,icross2,1,kth,'aerosol:asy')
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
      call getmem2d(aermmb,1,npoints,1,kz,'aerosol:aermmb')
      call getmem3d(ftota3d,1,npoints,0,kz,1,nband,'aerosol:ftota3d')

      ! these variables are defined on full rad grid including hat
      if (irrtm == 1) then
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
        ! FAB note that prof are always determined on kth level,
        ! even with standard scheme
        call getmem3d(extprof,jci1,jci2,ici1,ici2,1,kth,'rad:extprof')
        call getmem3d(asyprof,jci1,jci2,ici1,ici2,1,kth,'rad:asyprof')
        call getmem3d(ssaprof,jci1,jci2,ici1,ici2,1,kth,'rad:ssaprof')
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
        call getmem3d(extprofr4,1,npoints,1,kz,1,nband,'rad:extprofr4')
        call getmem3d(asyprofr4,1,npoints,1,kz,1,nband,'rad:asyprofr4')
        call getmem3d(ssaprofr4,1,npoints,1,kz,1,nband,'rad:ssaprofr4')
        call getmem1d(lambdaw,1,nband,'rad:lambdaw')
        call getmem1d(latr4,1,npoints,'aerosol:latr4')
        call getmem1d(lonr4,1,npoints,'aerosol:lonr4')
        call getmem1d(altr4,1,npoints,'aerosol:altr4')
        call getmem2d(z,1,npoints,1,kz,'aerosol:z')
        call getmem2d(dz,1,npoints,1,kz,'aerosol:dz')
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

      if ( iclimaaer /= 1 ) return

      ntr = aerclima_ntr
      nbin = aerclima_nbin
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
          ib = 1
          do i = ici1 , ici2
            do j = jci1 , jci2
              aermmr(ib,k,n) = max(w1*aerm1(j,i,k,n) + w2*aerm2(j,i,k,n),d_zero)
              ib = ib + 1
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
        icount(1) = jce2
        icount(2) = ice2
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

      write (stdout, *) 'Opening aerosol file : '//trim(aefile)

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
      logical , save :: lfirst
      logical :: dointerp
      real(rkx) , dimension(kth) :: opprnt
      real(rkx) :: xfac1 , xfac2 , odist
      type (rcm_time_and_date) :: imonmidd
      integer(ik4) :: iyear , imon , iday , ihour
      integer(ik4) :: k , im1 , iy1 , im2 , iy2
      integer(ik4) , save :: ism , isy
      type (rcm_time_and_date) :: iref1 , iref2
      type (rcm_time_interval) :: tdif
      data lfirst /.true./
      data ism /-1/
      data isy /-1/

      call split_idate(idatex,iyear,imon,iday,ihour)
      imonmidd = monmiddle(idatex)

      if (iyear < 1979 .and. iyear > 2021) then
        write (stderr,*) 'NO CLIMATIC AEROPP DATA AVAILABLE FOR ',iyear*100+imon
        return
      end if

      if ( lfirst ) then
        call grid_collect(m2r%xlon,alon,jce1,jce2,ice1,ice2)
        call grid_collect(m2r%xlat,alat,jce1,jce2,ice1,ice2)
        if ( myid == iocpu ) then
          call getfile(iyear,imon,ncid)

          call getmem1d(lat,1,clnlat,'aeropp:lat')
          call getmem1d(lon,1,clnlon,'aeropp:lon')
          call getmem3d(xext1,1,clnlon, &
                              1,clnlat,1,clnlev, 'aerosol:xext1')
          call getmem3d(xext2,1,clnlon, &
                              1,clnlat,1,clnlev, 'aerosol:xext2')
          call getmem3d(xssa1,1,clnlon, &
                              1,clnlat,1,clnlev, 'aerosol:xssa1')
          call getmem3d(xssa2,1,clnlon, &
                              1,clnlat,1,clnlev, 'aerosol:xssa2')
          call getmem3d(xasy1,1,clnlon, &
                              1,clnlat,1,clnlev, 'aerosol:xasy1')
          call getmem3d(xasy2,1,clnlon, &
                              1,clnlat,1,clnlev, 'aerosol:xasy2')
          call getmem3d(xdelp1,1,clnlon, &
                               1,clnlat,1,clnlev, 'aerosol:xasy2')
          call getmem3d(xdelp2,1,clnlon, &
                               1,clnlat,1,clnlev, 'aerosol:xasy2')
          call getmem3d(yext,1,njcross,1,nicross,1,clnlev,':yext')
          call getmem3d(yssa,1,njcross,1,nicross,1,clnlev,':yssa')
          call getmem3d(yasy,1,njcross,1,nicross,1,clnlev,':yasy')
          call getmem3d(ydelp,1,njcross,1,nicross,1,clnlev,':ydelp')
          call getmem3d(yphcl,1,njcross,1,nicross,1,clnlev,':yphcl')
          call init_aeroppdata(ncid,lat,lon)
          call h_interpolator_create(hint,lat,lon,alat,alon)
        end if
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
        call grid_collect(m2r%phatms,zp3d,jci1,jci2,ici1,ici2,1,kz)
        call grid_collect(m2r%deltaz,zdz3d,jci1,jci2,ici1,ici2,1,kz)
        if ( myid == iocpu ) then
          ! First build model pressure levels (in Pa) including
          ! radiative hat (up to kth)
          do k = 1, kz
            zpr3d(1:njcross,1:nicross,kth-kz+k) = zp3d(1:njcross,1:nicross,k)
            zdzr3d(1:njcross,1:nicross,kth-kz+k) = zdz3d(1:njcross,1:nicross,k)
          end do
          do k = 1, kth -kz
            zpr3d(1:njcross,1:nicross,kth-kz-k+1) = stdplevh(kclimh+k-1)*d_100
            zdzr3d(1:njcross,1:nicross,kth-kz-k+1) = &
              (stdhlevf(kclimh +k)-stdhlevf(kclimh +k -1)) * d_1000
          end do

          write (stdout,*) 'Reading EXT.,SSA,ASY Data...'
          if ( lfirst ) then
            call getfile(iy1,im1,ncid)
            call readvar3d(ncid,'EXTTOT',xext1)
            call readvar3d(ncid,'SSATOT',xssa1)
            call readvar3d(ncid,'GTOT',xasy1)
            call readvar3d(ncid,'DELP',xdelp1)

            call getfile(iy2,im2,ncid)
            call readvar3d(ncid,'EXTTOT',xext2)
            call readvar3d(ncid,'SSATOT',xssa2)
            call readvar3d(ncid,'GTOT',xasy2)
            call readvar3d(ncid,'DELP',xdelp2)

            xext1  = max(xext1,d_zero)
            xssa1  = min(max(xssa1,d_zero),d_one)
            xasy1  = min(max(xasy1,d_zero),d_one)
            xdelp1 = max(xdelp1,d_zero)
            xext2  = max(xext2,d_zero)
            xssa2  = min(max(xssa2,d_zero),d_one)
            xasy2  = min(max(xasy2,d_zero),d_one)
            xdelp2 = max(xdelp2,d_zero)

            call h_interpolate_cont(hint,xext1,yext)
            call h_interpolate_cont(hint,xssa1,yssa)
            call h_interpolate_cont(hint,xasy1,yasy)
            call h_interpolate_cont(hint,xdelp1,ydelp)
            !
            !               VERTICAL Interpolation
            !
            ! The pressure at the MERRA top is a fixed constant:
            !
            !                 PTOP = 0.01 hPa (1 Pa)
            !
            ! Pressures at edges should be computed by summing the DELP
            ! pressure thickness starting at PTOP.
            ! A representative pressure for the layer can then be obtained
            ! from these. Here just use linear av for now, should be improved !!
            ! MERRA grid is top down
            yphcl(1:njcross,1:nicross,1)  = d_one + &
                ydelp(1:njcross,1:nicross,1)*d_half
            do k = 2 , clnlev
              yphcl(1:njcross,1:nicross,k) = &
                 yphcl(1:njcross,1:nicross,k-1)  + &
                 (ydelp(1:njcross,1:nicross,k-1) + &
                  ydelp(1:njcross,1:nicross,k))*d_half
            end do

            call intp1(ext1,yext,zpr3d,yphcl,njcross,nicross, &
                       kth,clnlev,0.7_rkx,0.7_rkx,0.4_rkx)
            call intp1(ssa1,yssa,zpr3d,yphcl,njcross,nicross, &
                       kth,clnlev,0.7_rkx,0.7_rkx,0.4_rkx)
            call intp1(asy1,yasy,zpr3d,yphcl,njcross,nicross, &
                       kth,clnlev,0.7_rkx,0.7_rkx,0.4_rkx)

            ! same for ext2
            call h_interpolate_cont(hint,xext2,yext)
            call h_interpolate_cont(hint,xssa2,yssa)
            call h_interpolate_cont(hint,xasy2,yasy)
            call h_interpolate_cont(hint,xdelp2,ydelp)
            yphcl(1:njcross,1:nicross,1) = d_one + &
                ydelp(1:njcross,1:nicross,1)*d_half
            do k = 2,clnlev
              yphcl(1:njcross,1:nicross,k) = &
                  yphcl(1:njcross,1:nicross,k-1)  + &
                  (ydelp(1:njcross,1:nicross,k-1) + &
                   ydelp(1:njcross,1:nicross,k))*d_half
            end do
            call intp1(ext2,yext,zpr3d,yphcl,njcross,nicross, &
                       kth,clnlev,0.7_rkx,0.7_rkx,0.4_rkx)
            call intp1(ssa2,yssa,zpr3d,yphcl,njcross,nicross, &
                       kth,clnlev,0.7_rkx,0.7_rkx,0.4_rkx)
            call intp1(asy2,yasy,zpr3d,yphcl,njcross,nicross, &
                       kth,clnlev,0.7_rkx,0.7_rkx,0.4_rkx)

          else  ! ( if not first call , just update ext2 )

            ext1 = ext2
            ssa1 = ssa2
            asy1 = asy2
            call getfile(iy2,im2,ncid)
            call readvar3d(ncid,'EXTTOT',xext2)
            call readvar3d(ncid,'SSATOT',xssa2)
            call readvar3d(ncid,'GTOT',xasy2)
            call readvar3d(ncid,'DELP',xdelp2)

            xext2 = max(xext2,d_zero)
            xssa2 = min(max(xssa2,d_zero),d_one)
            xasy2 = min(max(xasy2,d_zero),d_one)
            xdelp2= max(xdelp2,d_zero)

            call h_interpolate_cont(hint,xext2,yext)
            call h_interpolate_cont(hint,xssa2,yssa)
            call h_interpolate_cont(hint,xasy2,yasy)
            call h_interpolate_cont(hint,xdelp2,ydelp)

            yphcl(1:njcross,1:nicross,1) = d_one + &
                ydelp(1:njcross,1:nicross,1)*d_half
            do k = 2,clnlev
              yphcl(1:njcross,1:nicross,k) = &
                  yphcl(1:njcross,1:nicross,k-1)  + &
                  (ydelp(1:njcross,1:nicross,k-1) + &
                   ydelp(1:njcross,1:nicross,k))*d_half
            end do

            call intp1(ext2,yext,zpr3d,yphcl,njcross,nicross, &
                       kth,clnlev,0.7_rkx,0.7_rkx,0.4_rkx)
            call intp1(ssa2,yssa,zpr3d,yphcl,njcross,nicross, &
                       kth,clnlev,0.7_rkx,0.7_rkx,0.4_rkx)
            call intp1(asy2,yasy,zpr3d,yphcl,njcross,nicross, &
                       kth,clnlev,0.7_rkx,0.7_rkx,0.4_rkx)
          end if
        end if
      end if ! end of interp
      if ( myid == iocpu ) then
        ! FAB :  normally we should perform vertical interpolation here
        ! rather than in dointerp, since it depends on instantaneous model
        ! pressure field !
        tdif = idatex-iref1
        xfac1 = real(tohours(tdif),rkx)
        tdif = idatex-iref2
        xfac2 = real(tohours(tdif),rkx)
        odist = xfac1 - xfac2
        xfac1 = xfac1 / odist
        xfac2 = d_one - xfac1
        ! ! Important :  radiation schemes expect AOD per layer, calculated
        !   from extinction
        ext = (ext1*xfac2 + ext2*xfac1) * zdzr3d
        ssa = ssa1*xfac2 + ssa2*xfac1
        asy = asy1*xfac2 + asy2*xfac1

        where ( ext < 1.E-10_rkx ) ext = 1.E-10_rkx
        where ( ssa < 1.E-10_rkx ) ssa = 0.991_rkx
        where ( asy < 1.E-10_rkx ) ssa = 0.611_rkx
      end if
      call grid_distribute(ext,extprof,jci1,jci2,ici1,ici2,1,kth)
      call grid_distribute(ssa,ssaprof,jci1,jci2,ici1,ici2,1,kth)
      call grid_distribute(asy,asyprof,jci1,jci2,ici1,ici2,1,kth)
      if ( myid == italk .and. dointerp ) then
        do k = 1 , kth
          opprnt(k) = extprof(3,3,kth-k+1)
        end do
        call vprntv(opprnt,kth,'Updated ext profile at (3,3)')
        do k = 1 , kth
          opprnt(k) = ssaprof(3,3,kth-k+1)
        end do
        call vprntv(opprnt,kth,'Updated ssa profile at (3,3)')
        do k = 1 , kth
          opprnt(k) = asyprof(3,3,kth-k+1)
        end do
        call vprntv(opprnt,kth,'Updated asy profile at (3,3)')
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
      integer(ik4) , save , dimension(4) :: istart , icount
      integer(ik4) :: iret , irec
      data icvar /-1/
      data istart  /  1 ,  1 ,  1 ,  1/

      icount(1) = clnlon
      icount(2) = clnlat
      icount(3) = clnlev
      icount(4) = 1
      ! irec = ((iyear-cliyear)*12+imon)-1
      irec = 1
      iret = nf90_inq_varid(ncid,vname,icvar)
      if ( iret /= nf90_noerr ) then
        write (stderr, *) nf90_strerror(iret)
        call fatal(__FILE__,__LINE__,'CANNOT READ FROM AEROPPCLIM FILE')
      end if
      istart(4) = irec
      iret = nf90_get_var(ncid,icvar,val,istart,icount)
      if ( iret /= nf90_noerr ) then
        write (stderr, *) nf90_strerror(iret)
        call fatal(__FILE__,__LINE__,'CANNOT READ FROM AEROPPCLIM FILE')
      end if
    end subroutine readvar3d
    !
    !-----------------------------------------------------------------------
    !
    ! SUBROUTINE AERMIX
    !
    ! Set global mean tropospheric aerosol
    !
    ! Specify aerosol mixing ratio and compute relative humidity for later
    ! adjustment of aerosol optical properties. Aerosol mass mixing ratio
    ! is specified so that the column visible aerosol optical depth is a
    ! specified global number (tauvis). This means that the actual mixing
    ! ratio depends on pressure thickness of the lowest three atmospheric
    ! layers near the surface.
    !
    ! Optical properties and relative humidity parameterization are from:
    !
    ! J.T. Kiehl and B.P. Briegleb  "The Relative Roles of Sulfate Aerosols
    ! and Greenhouse Gases in Climate Forcing"  Science  260  pp311-314
    ! 16 April 1993
    !
    ! Visible (vis) here means 0.5-0.7 micro-meters
    ! Forward scattering fraction is taken as asymmetry parameter squared
    !
    !---------------------------Code history--------------------------------
    !
    ! Original version:  B. Briegleb  March 1995
    ! Standarized:       L. Buja,     Feb 1996
    ! Reviewed:          B. Briegleb, Mar 1996
    !
    !-----------------------------------------------------------------------
    !
    subroutine aermix(pint,n1,n2)
      implicit none
      integer(ik4) , intent(in) :: n1 , n2
      ! Radiation level interface pressures (dynes/cm2)
      real(rkx) , intent(in) , pointer , dimension(:,:) :: pint
      !
      !-----------------------------------------------------------------------
      !
      ! mxaerl - max nmbr aerosol levels counting up from surface
      ! tauvis - visible optical depth
      ! kaervs - visible extinction coefficiant of aerosol (m2/g)
      ! omgvis - visible omega0
      ! gvis   - visible forward scattering asymmetry parameter
      !
      !-----------------------------------------------------------------------
      !
      real(rkx) :: gvis , kaervs , omgvis , rhfac , tauvis
      integer(ik4) :: n , k , mxaerl

      data kaervs /5.3012_rkx/        ! multiplication factor for kaer
      data omgvis /0.999999_rkx/
      data gvis   /0.694889_rkx/
      data rhfac  /1.6718_rkx/        ! EES added for efficiency
      !
      !-----------------------------------------------------------------------
      !
      mxaerl = 4
      !
      !fil  tauvis = 0.01_rkx
      tauvis = 0.04_rkx
      !
      ! Set relative humidity and factor; then aerosol amount:
      !
      do k = 1 , kz
        do n = n1 , n2
          !
          ! Define background aerosol
          ! Find constant aerosol mass mixing ratio for specified levels
          ! in the column, converting units where appropriate
          ! for the moment no more used
          !
          if ( k >= kz + 1 - mxaerl ) then
            aermmb(n,k) = egravgts * tauvis / &
                         (d10e4*kaervs*rhfac*(d_one-omgvis*gvis*gvis) * &
                         (pint(n,kzp1)-pint(n,kzp1-mxaerl)))
          else
            aermmb(n,k) = d_zero
          end if
        end do
      end do
    end subroutine aermix
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
      integer(ik4) :: j , i , k , visband
      real(rkx) :: uaerdust , qabslw , rh0

      !-
      ! Aerosol forced by Optical Properties Climatology
      ! distinguish between standard scheme and RRTM scheme
      ! Rq extprof has been scaled for layer height before
      ! ( represents here the layer AOD )
      ! Directly force aerosol properties passed to radiation driver and exit

      if ( iclimaaer == 2 ) then
        aertrlw(:,:,:) = d_one
        if ( irrtm == 1 ) then
          visband = 10
           ! FAB try only the visible RRTM  now
          ns = visband
          do k = 1 , kth
            n = 1
            do i = ici1 , ici2
              do j = jci1 , jci2
                tauxar3d(n,k,ns) = extprof(j,i,k)  !already scaled, top down regcm
                tauasc3d(n,k,ns) = ssaprof(j,i,k)
                gtota3d(n,k,ns)  = asyprof(j,i,k)
                n = n + 1
              end do
            end do
          end do
        else if ( irrtm == 0 ) then
          visband = 8
          ns = visband
          tauxar(:,ns) = d_zero
          tauasc(:,ns) = d_zero
          gtota(:,ns) = d_zero
          ftota(:,ns) = d_zero
          ! adapt the clim vert grid (1 to kth) to the standard
          ! rad grid ( 0 to kz) first Treat the top radiative layer
          ! tauxar3d(n,0,ns)
          if ( kth > kz ) then
            n = 1
            do i = ici1 , ici2
              do j = jci1 , jci2
                tauxar3d(n,0,ns) = sum(extprof(j,i,1:kth-kz) )
                tauasc3d(n,0,ns) = sum(ssaprof(j,i,1:kth-kz)) / &
                                  real(kth-kz,rkx) * tauxar3d(n,0,ns)
                gtota3d(n,0,ns) = sum(asyprof(j,i,1:kth-kz)) / &
                                  real(kth-kz,rkx) * tauasc3d(n,0,ns)
                ftota3d(n,0,ns) = (sum(asyprof(j,i,1:kth-kz)) / &
                                  real(kth-kz,rkx))**2 * tauasc3d(n,0,ns)
                n = n +1
              end do
            end do
          end if
          do k = 1 , kz
            n = 1
            do i = ici1 , ici2
              do j = jci1 , jci2
                ! already scaled for layer height
                ! grid is top down
                tauxar3d(n,k,ns) = extprof(j,i,kth-kz+k)
                ! here the standard scheme expect layer scaled quantity
                tauasc3d(n,k,ns) = ssaprof(j,i,kth-kz+k) * tauxar3d(n,k,ns)
                gtota3d(n,k,ns)  = asyprof(j,i,kth-kz+k) * tauasc3d(n,k,ns)
                ftota3d(n,k,ns)  = asyprof(j,i,kth-kz+k)**2 * tauasc3d(n,k,ns)
                tauxar(n,ns) = tauxar(n,ns) + tauxar3d(n,k,ns)
                tauasc(n,ns) = tauasc(n,ns) + tauasc3d(n,k,ns)
                gtota(n,ns)  = gtota(n,ns)  + gtota3d(n,k,ns)
                ftota(n,ns)  = ftota(n,ns)  + ftota3d(n,k,ns)
                n = n + 1
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
            do k = 1 , kz
              do i = 1 , npoints
                tauxar3d(i,k,n) = extprofr4(i,k,n)
                tauasc3d(i,k,n) = ssaprofr4(i,k,n)
                gtota3d(i,k,n)  = asyprofr4(i,k,n)
              end do
            end do
          end do
        else if ( irrtm == 0 ) then
          do n = 1 , nband
            do k = 1 , kz
              do i = 1 , npoints
                ! already scaled for layer height
                tauxar3d(i,k,n) = extprofr4(i,k,n)
                ! here the standard scheme expect layer scaled quantity
                tauasc3d(i,k,n) = ssaprofr4(i,k,n) * tauxar3d(i,k,n)
                gtota3d(i,k,n)  = asyprofr4(i,k,n) * &
                       ssaprofr4(i,k,n) * tauxar3d(i,k,n)
                ftota3d(i,k,n)  = asyprofr4(i,k,n)**2 * &
                       ssaprofr4(i,k,n) * tauxar3d(i,k,n)

                ! define also tauxar for std scheme clear sky diagnostics
                tauxar(i,n) = tauxar(i,n) + tauxar3d(i,k,n)
                tauasc(i,n) = tauasc(i,n) + tauasc3d(i,k,n)
                ftota(i,n) =  ftota(i,n)  + ftota3d(i,k,n)
                gtota(i,n) =  gtota(i,n)  + gtota3d(i,k,n)
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
        tauxar(:,:) = d_zero
        tauasc(:,:) = d_zero
        gtota(:,:) = d_zero
        ftota(:,:) = d_zero
        tauxar3d(:,:,:) = d_zero
        tauasc3d(:,:,:) = d_zero
        gtota3d(:,:,:) = d_zero
        ftota3d(:,:,:) = d_zero
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
      do k = 1 , kz
        do n = n1 , n2
          path(n,k) = (pint(n,k+1)-pint(n,k))*regravgts
        end do
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
          else if ( chtrname(itr)(1:3) == 'SO4' .or.    &
                    chtrname(itr)(1:5) == 'H2SO4'.or.    &
                    chtrname(itr)(1:4) == 'ANO3' .or.   &
                    chtrname(itr)(1:4) == 'ANH4' ) then
            do k = 1 , kz
              do n = n1 , n2
                rh0 = min(0.97_rkx,max(d_zero,rh(n,k)))
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
                rh0 = min(0.99_rkx,max(d_zero,rh(n,k)))
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
                rh0 = min(0.99_rkx,max(d_zero,rh(n,k)))
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
                rh0 = min(0.99_rkx,max(d_zero,rh(n,k)))
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
                rh0 = min(0.99_rkx,max(d_zero,rh(n,k)))
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
                rh0 = min(0.99_rkx,max(d_zero,rh(n,k)))
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
        if ( irrtm == 1 ) then
          do k = 0 , kz
            do n = n1 , n2
              if ( tauxar3d(n,k,ns) > d_zero ) then
                tauasc3d(n,k,ns) = tauasc3d(n,k,ns) / tauxar3d(n,k,ns)
                gtota3d(n,k,ns) = gtota3d(n,k,ns) / &
                  (tauasc3d(n,k,ns)*tauxar3d(n,k,ns))
              else
                tauasc3d(n,k,ns) = 0.999999_rkx
                gtota3d(n,k,ns) = 0.5_rkx
              end if
            end do
          end do
        end if
      end do ! end spectral loop
      !
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
          tauxar3d_lw(:,:,ns) = d_zero
          ibin = 0
          do itr = 1 , ntr
            if ( chtrname(itr)(1:4) == 'DUST') then
              ibin = ibin + 1
              do k = 1 , kz
                do n = n1 , n2
                  uaer(n,k,itr) = aermmr(n,k,itr)*path(n,k)
                  tx(n,k,itr) = d10e5*uaer(n,k,itr)*ksdust_lw(ns,ibin)
                  ! add the extinction for every bins
                  tauxar3d_lw(n,k,ns) = tauxar3d_lw(n,k,ns) + tx(n,k,itr)
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
      integer(ik4) :: ibin , i , j , k , n
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
          do i = ici1 , ici2
            do j = jci1 , jci2
              z(ibin,k)  = m2r%za(j,i,k)
              dz(ibin,k) = m2r%deltaz(j,i,k)
              ibin = ibin + 1
            end do
          end do
        end do
        year_fr = real(iy) + real(yeardayfrac(x))/real(yeardays(iy,x%calendar))
        do n = 1 , nband
          call sp_aop_profile(macv2sp_hist,macv2sp_scen, &
                              kz,npoints,lambdaw(n), &
                              altr4,lonr4,latr4,year_fr,z,dz, &
                              dnovrnr4,extprofr4(:,:,n), &
                              ssaprofr4(:,:,n),asyprofr4(:,:,n))
        end do
        if ( myid == italk ) then
          write(stdout,*) 'Updating aerosol optical properties...'
        end if
        idlast = id
      end if
    end subroutine cmip6_plume_profile

    subroutine getfile(year,month,ncid)
      implicit none
      integer(ik4) , intent(in) :: year,month
      integer(ik4) , intent(inout) ::  ncid
      character(len=256) :: infile
      integer(ik4) :: iret , idimid
      write(infile,'(A,I4,I0.2,A)') &
        trim(radclimpath)//pthsep//'MERRA2_OPPMONTH_wb10.',year,month,'.nc'
      if ( ncid < 0 ) then
        iret = nf90_open(infile,nf90_nowrite,ncid)
        if ( iret /= nf90_noerr ) then
          write (stderr, *) nf90_strerror(iret), trim(infile)
          call fatal(__FILE__,__LINE__,'CANNOT OPEN AEROSOL OP.PROP CLIM FILE')
        else
          write(stdout,*) 'AEROPP file open ', trim(infile)
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
        write(stdout,*) 'AEROPP file open ', trim(infile)
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

end module mod_rad_aerosol

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

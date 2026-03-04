!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_stdatm

  use mod_intkinds
  use mod_realkinds
  use mod_constants

  implicit none

  private

  integer(ik4), public, parameter :: itropical = 1
  integer(ik4), public, parameter :: imidlatsummer = 2
  integer(ik4), public, parameter :: imidlatwinter = 3
  integer(ik4), public, parameter :: ipolarsummer = 4
  integer(ik4), public, parameter :: ipolarwinter = 5

  integer(ik4), public, parameter :: istdatm_hgtkm = 1  ! HGT in km
  integer(ik4), public, parameter :: istdatm_prsmb = 2  ! PRESS in hPa
  integer(ik4), public, parameter :: istdatm_tempk = 3  ! TEMP in K
  integer(ik4), public, parameter :: istdatm_airdn = 4  ! RHO in g/m^3
  integer(ik4), public, parameter :: istdatm_qdens = 5  ! Q in g/m^3
  integer(ik4), public, parameter :: istdatm_ozone = 6  ! Ozone in g/m^3

  integer(ik4), public, parameter :: n_atmzones = 5
  integer(ik4), public, parameter :: n_atmparms = 6
  integer(ik4), public, parameter :: n_atmlevls = 31
  integer(ik4), public, parameter :: n_preflev  = 31
  integer(ik4), public, parameter :: n_prehlev  = 30
  integer(ik4), public, parameter :: n_hreflev  = 31
  integer(ik4), public, parameter :: n_hrehlev  = 30
  real(rk8), dimension(n_hrehlev), parameter :: stdhlevh = &
  [    0.5_rk8,  1.5_rk8,  2.5_rk8,  3.5_rk8,  4.5_rk8,  5.5_rk8, &
       6.5_rk8,  7.5_rk8,  8.5_rk8,  9.5_rk8, 10.5_rk8, 11.5_rk8, &
      12.5_rk8, 13.5_rk8, 14.5_rk8, 15.5_rk8, 16.5_rk8, 17.5_rk8, &
      18.5_rk8, 19.5_rk8, 20.5_rk8, 21.5_rk8, 22.5_rk8, 23.5_rk8, &
      24.5_rk8, 27.5_rk8, 32.5_rk8, 37.5_rk8, 42.5_rk8, 47.5_rk8 ]
  real(rk8), dimension(n_hreflev), parameter :: stdhlevf = &
  [    0.0_rk8,  1.0_rk8,  2.0_rk8,  3.0_rk8,  4.0_rk8,  5.0_rk8, &
       6.0_rk8,  7.0_rk8,  8.0_rk8,  9.0_rk8, 10.0_rk8, 11.0_rk8, &
      12.0_rk8, 13.0_rk8, 14.0_rk8, 15.0_rk8, 16.0_rk8, 17.0_rk8, &
      18.0_rk8, 19.0_rk8, 20.0_rk8, 21.0_rk8, 22.0_rk8, 23.0_rk8, &
      24.0_rk8, 25.0_rk8, 30.0_rk8, 35.0_rk8, 40.0_rk8, 45.0_rk8, &
      50.0_rk8 ]
  real(rk8), dimension(n_prehlev), parameter :: stdplevh = &
  [  950.0_rk8, 850.0_rk8, 750.0_rk8, 675.0_rk8, 600.0_rk8, 525.0_rk8, &
     475.0_rk8, 415.0_rk8, 350.0_rk8, 300.0_rk8, 265.0_rk8, 215.0_rk8, &
     200.0_rk8, 165.0_rk8, 135.0_rk8, 115.0_rk8,  95.0_rk8,  85.0_rk8, &
      70.0_rk8,  60.0_rk8,  50.0_rk8,  42.5_rk8,  37.5_rk8,  32.5_rk8, &
      27.5_rk8,  17.5_rk8,  10.0_rk8,   5.0_rk8,   2.0_rk8,   1.0_rk8 ]
  real(rk8), dimension(n_preflev), parameter :: stdplevf = &
  [ 1013.0_rk8, 900.0_rk8, 800.0_rk8, 700.0_rk8, 650.0_rk8, 550.0_rk8, &
     500.0_rk8, 450.0_rk8, 375.0_rk8, 325.0_rk8, 275.0_rk8, 250.0_rk8, &
     225.0_rk8, 175.0_rk8, 150.0_rk8, 125.0_rk8, 100.0_rk8,  90.0_rk8, &
      75.0_rk8,  65.0_rk8,  55.0_rk8,  45.0_rk8,  40.0_rk8,  35.0_rk8, &
      30.0_rk8,  25.0_rk8,  12.5_rk8,   7.5_rk8,   2.5_rk8,   1.5_rk8, &
       0.5_rk8 ]

  real(rk8), dimension(n_atmparms,n_atmlevls,n_atmzones), &
    parameter :: stdatm = reshape( &
  [ 0.0_rk8,1013.00_rk8,300.00_rk8,0.1167e+4_rk8,0.190e+2_rk8,0.560e-4_rk8, &
    1.0_rk8, 904.00_rk8,294.00_rk8,0.1064e+4_rk8,0.130e+2_rk8,0.560e-4_rk8, &
    2.0_rk8, 805.00_rk8,288.00_rk8,0.9689e+3_rk8,0.930e+1_rk8,0.540e-4_rk8, &
    3.0_rk8, 715.00_rk8,284.00_rk8,0.8756e+3_rk8,0.470e+1_rk8,0.510e-4_rk8, &
    4.0_rk8, 633.00_rk8,277.00_rk8,0.7951e+3_rk8,0.220e+1_rk8,0.470e-4_rk8, &
    5.0_rk8, 559.00_rk8,270.00_rk8,0.7199e+3_rk8,0.150e+1_rk8,0.450e-4_rk8, &
    6.0_rk8, 492.00_rk8,264.00_rk8,0.6501e+3_rk8,0.850e+0_rk8,0.430e-4_rk8, &
    7.0_rk8, 432.00_rk8,257.00_rk8,0.5855e+3_rk8,0.470e+0_rk8,0.410e-4_rk8, &
    8.0_rk8, 378.00_rk8,250.00_rk8,0.5258e+3_rk8,0.250e+0_rk8,0.390e-4_rk8, &
    9.0_rk8, 329.00_rk8,244.00_rk8,0.4708e+3_rk8,0.120e+0_rk8,0.390e-4_rk8, &
   10.0_rk8, 286.00_rk8,237.00_rk8,0.4202e+3_rk8,0.500e-1_rk8,0.390e-4_rk8, &
   11.0_rk8, 247.00_rk8,230.00_rk8,0.3740e+3_rk8,0.170e-1_rk8,0.410e-4_rk8, &
   12.0_rk8, 213.00_rk8,224.00_rk8,0.3316e+3_rk8,0.600e-2_rk8,0.430e-4_rk8, &
   13.0_rk8, 182.00_rk8,217.00_rk8,0.2929e+3_rk8,0.180e-2_rk8,0.450e-4_rk8, &
   14.0_rk8, 156.00_rk8,210.00_rk8,0.2578e+3_rk8,0.100e-2_rk8,0.450e-4_rk8, &
   15.0_rk8, 132.00_rk8,204.00_rk8,0.2260e+3_rk8,0.760e-3_rk8,0.470e-4_rk8, &
   16.0_rk8, 111.00_rk8,197.00_rk8,0.1972e+3_rk8,0.640e-3_rk8,0.470e-4_rk8, &
   17.0_rk8,  93.70_rk8,195.00_rk8,0.1676e+3_rk8,0.560e-3_rk8,0.690e-4_rk8, &
   18.0_rk8,  78.90_rk8,199.00_rk8,0.1382e+3_rk8,0.500e-3_rk8,0.900e-4_rk8, &
   19.0_rk8,  66.60_rk8,203.00_rk8,0.1145e+3_rk8,0.490e-3_rk8,0.140e-3_rk8, &
   20.0_rk8,  56.50_rk8,207.00_rk8,0.9515e+2_rk8,0.450e-3_rk8,0.190e-3_rk8, &
   21.0_rk8,  48.00_rk8,211.00_rk8,0.7938e+2_rk8,0.510e-3_rk8,0.240e-3_rk8, &
   22.0_rk8,  40.90_rk8,215.00_rk8,0.6645e+2_rk8,0.510e-3_rk8,0.280e-3_rk8, &
   23.0_rk8,  35.00_rk8,217.00_rk8,0.5618e+2_rk8,0.540e-3_rk8,0.320e-3_rk8, &
   24.0_rk8,  30.00_rk8,219.00_rk8,0.4763e+2_rk8,0.600e-3_rk8,0.340e-3_rk8, &
   25.0_rk8,  25.70_rk8,221.00_rk8,0.4045e+2_rk8,0.670e-3_rk8,0.340e-3_rk8, &
   30.0_rk8,  12.20_rk8,232.00_rk8,0.1831e+2_rk8,0.360e-3_rk8,0.240e-3_rk8, &
   35.0_rk8,   6.00_rk8,243.00_rk8,0.8600e+1_rk8,0.110e-3_rk8,0.920e-4_rk8, &
   40.0_rk8,   3.05_rk8,254.00_rk8,0.4181e+1_rk8,0.430e-4_rk8,0.410e-4_rk8, &
   45.0_rk8,   1.59_rk8,265.00_rk8,0.2097e+1_rk8,0.190e-4_rk8,0.130e-4_rk8, &
   50.0_rk8,   0.85_rk8,270.00_rk8,0.1101e+1_rk8,0.630e-5_rk8,0.430e-5_rk8, &
    0.0_rk8,1013.00_rk8,294.00_rk8,0.1191e+4_rk8,0.140e+2_rk8,0.600e-4_rk8, &
    1.0_rk8, 902.00_rk8,290.00_rk8,0.1080e+4_rk8,0.930e+1_rk8,0.600e-4_rk8, &
    2.0_rk8, 802.00_rk8,285.00_rk8,0.9757e+3_rk8,0.590e+1_rk8,0.600e-4_rk8, &
    3.0_rk8, 710.00_rk8,279.00_rk8,0.8846e+3_rk8,0.330e+1_rk8,0.620e-4_rk8, &
    4.0_rk8, 628.00_rk8,273.00_rk8,0.7998e+3_rk8,0.190e+1_rk8,0.640e-4_rk8, &
    5.0_rk8, 554.00_rk8,267.00_rk8,0.7211e+3_rk8,0.100e+1_rk8,0.660e-4_rk8, &
    6.0_rk8, 487.00_rk8,261.00_rk8,0.6487e+3_rk8,0.610e+0_rk8,0.690e-4_rk8, &
    7.0_rk8, 426.00_rk8,255.00_rk8,0.5830e+3_rk8,0.370e+0_rk8,0.750e-4_rk8, &
    8.0_rk8, 372.00_rk8,248.00_rk8,0.5225e+3_rk8,0.210e+0_rk8,0.790e-4_rk8, &
    9.0_rk8, 324.00_rk8,242.00_rk8,0.4669e+3_rk8,0.120e+0_rk8,0.860e-4_rk8, &
   10.0_rk8, 281.00_rk8,235.00_rk8,0.4159e+3_rk8,0.640e-1_rk8,0.900e-4_rk8, &
   11.0_rk8, 243.00_rk8,229.00_rk8,0.3693e+3_rk8,0.220e-1_rk8,0.110e-3_rk8, &
   12.0_rk8, 209.00_rk8,222.00_rk8,0.3269e+3_rk8,0.600e-2_rk8,0.120e-3_rk8, &
   13.0_rk8, 179.00_rk8,216.00_rk8,0.2882e+3_rk8,0.180e-2_rk8,0.150e-3_rk8, &
   14.0_rk8, 153.00_rk8,216.00_rk8,0.2464e+3_rk8,0.100e-2_rk8,0.180e-3_rk8, &
   15.0_rk8, 130.00_rk8,216.00_rk8,0.2104e+3_rk8,0.760e-3_rk8,0.190e-3_rk8, &
   16.0_rk8, 111.00_rk8,216.00_rk8,0.1797e+3_rk8,0.640e-3_rk8,0.210e-3_rk8, &
   17.0_rk8,  95.00_rk8,216.00_rk8,0.1535e+3_rk8,0.560e-3_rk8,0.240e-3_rk8, &
   18.0_rk8,  81.20_rk8,216.00_rk8,0.1305e+3_rk8,0.500e-3_rk8,0.280e-3_rk8, &
   19.0_rk8,  69.50_rk8,217.00_rk8,0.1110e+3_rk8,0.490e-3_rk8,0.320e-3_rk8, &
   20.0_rk8,  59.50_rk8,218.00_rk8,0.9453e+2_rk8,0.450e-3_rk8,0.340e-3_rk8, &
   21.0_rk8,  51.00_rk8,219.00_rk8,0.8056e+2_rk8,0.510e-3_rk8,0.360e-3_rk8, &
   22.0_rk8,  43.70_rk8,220.00_rk8,0.6872e+2_rk8,0.510e-3_rk8,0.360e-3_rk8, &
   23.0_rk8,  37.60_rk8,222.00_rk8,0.5867e+2_rk8,0.540e-3_rk8,0.340e-3_rk8, &
   24.0_rk8,  32.20_rk8,223.00_rk8,0.5014e+2_rk8,0.600e-3_rk8,0.320e-3_rk8, &
   25.0_rk8,  27.70_rk8,224.00_rk8,0.4288e+2_rk8,0.670e-3_rk8,0.300e-3_rk8, &
   30.0_rk8,  13.20_rk8,234.00_rk8,0.1322e+2_rk8,0.360e-3_rk8,0.200e-3_rk8, &
   35.0_rk8,   6.52_rk8,245.00_rk8,0.6519e+1_rk8,0.110e-3_rk8,0.920e-4_rk8, &
   40.0_rk8,   3.33_rk8,258.00_rk8,0.3330e+1_rk8,0.430e-4_rk8,0.410e-4_rk8, &
   45.0_rk8,   1.76_rk8,270.60_rk8,0.1757e+1_rk8,0.190e-4_rk8,0.130e-4_rk8, &
   50.0_rk8,   0.95_rk8,276.00_rk8,0.9512e+0_rk8,0.630e-5_rk8,0.430e-5_rk8, &
    0.0_rk8,1018.00_rk8,272.20_rk8,0.1301e+4_rk8,0.350e+1_rk8,0.600e-4_rk8, &
    1.0_rk8, 897.30_rk8,268.70_rk8,0.1162e+4_rk8,0.250e+1_rk8,0.540e-4_rk8, &
    2.0_rk8, 789.70_rk8,265.20_rk8,0.1037e+4_rk8,0.180e+1_rk8,0.490e-4_rk8, &
    3.0_rk8, 693.80_rk8,261.70_rk8,0.9230e+3_rk8,0.120e+1_rk8,0.490e-4_rk8, &
    4.0_rk8, 608.10_rk8,255.70_rk8,0.8282e+3_rk8,0.660e+0_rk8,0.490e-4_rk8, &
    5.0_rk8, 531.30_rk8,249.70_rk8,0.7411e+3_rk8,0.380e+0_rk8,0.580e-4_rk8, &
    6.0_rk8, 462.70_rk8,243.70_rk8,0.6614e+3_rk8,0.210e+0_rk8,0.640e-4_rk8, &
    7.0_rk8, 401.60_rk8,237.70_rk8,0.5886e+3_rk8,0.850e-1_rk8,0.770e-4_rk8, &
    8.0_rk8, 347.30_rk8,231.70_rk8,0.5222e+3_rk8,0.350e-1_rk8,0.900e-4_rk8, &
    9.0_rk8, 299.20_rk8,225.70_rk8,0.4619e+3_rk8,0.160e-1_rk8,0.120e-3_rk8, &
   10.0_rk8, 256.80_rk8,219.70_rk8,0.4072e+3_rk8,0.750e-2_rk8,0.160e-3_rk8, &
   11.0_rk8, 219.90_rk8,219.20_rk8,0.3496e+3_rk8,0.690e-2_rk8,0.210e-3_rk8, &
   12.0_rk8, 188.20_rk8,218.70_rk8,0.2999e+3_rk8,0.600e-2_rk8,0.260e-3_rk8, &
   13.0_rk8, 161.00_rk8,218.20_rk8,0.2572e+3_rk8,0.180e-2_rk8,0.300e-3_rk8, &
   14.0_rk8, 137.80_rk8,217.70_rk8,0.2206e+3_rk8,0.100e-2_rk8,0.320e-3_rk8, &
   15.0_rk8, 117.80_rk8,217.20_rk8,0.1890e+3_rk8,0.760e-3_rk8,0.340e-3_rk8, &
   16.0_rk8, 100.70_rk8,216.70_rk8,0.1620e+3_rk8,0.640e-3_rk8,0.360e-3_rk8, &
   17.0_rk8,  86.10_rk8,216.20_rk8,0.1388e+3_rk8,0.560e-3_rk8,0.390e-3_rk8, &
   18.0_rk8,  73.50_rk8,215.70_rk8,0.1188e+3_rk8,0.500e-3_rk8,0.410e-3_rk8, &
   19.0_rk8,  62.80_rk8,215.20_rk8,0.1017e+3_rk8,0.490e-3_rk8,0.430e-3_rk8, &
   20.0_rk8,  53.70_rk8,215.20_rk8,0.8690e+2_rk8,0.450e-3_rk8,0.450e-3_rk8, &
   21.0_rk8,  45.80_rk8,215.20_rk8,0.7421e+2_rk8,0.510e-3_rk8,0.430e-3_rk8, &
   22.0_rk8,  39.10_rk8,215.20_rk8,0.6338e+2_rk8,0.510e-3_rk8,0.430e-3_rk8, &
   23.0_rk8,  33.40_rk8,215.20_rk8,0.5415e+2_rk8,0.540e-3_rk8,0.390e-3_rk8, &
   24.0_rk8,  28.60_rk8,215.20_rk8,0.4624e+2_rk8,0.600e-3_rk8,0.360e-3_rk8, &
   25.0_rk8,  24.30_rk8,215.20_rk8,0.3950e+2_rk8,0.670e-3_rk8,0.340e-3_rk8, &
   30.0_rk8,  11.10_rk8,217.40_rk8,0.1783e+2_rk8,0.360e-3_rk8,0.190e-3_rk8, &
   35.0_rk8,   5.18_rk8,227.80_rk8,0.7924e+1_rk8,0.110e-3_rk8,0.920e-4_rk8, &
   40.0_rk8,   2.53_rk8,243.20_rk8,0.3625e+1_rk8,0.430e-4_rk8,0.410e-4_rk8, &
   45.0_rk8,   1.29_rk8,258.50_rk8,0.1741e+1_rk8,0.190e-4_rk8,0.130e-4_rk8, &
   50.0_rk8,   0.68_rk8,265.70_rk8,0.8954e+0_rk8,0.630e-5_rk8,0.430e-5_rk8, &
    0.0_rk8,1010.00_rk8,287.00_rk8,0.1220e+4_rk8,0.910e+1_rk8,0.490e-4_rk8, &
    1.0_rk8, 896.00_rk8,282.00_rk8,0.1110e+4_rk8,0.600e+1_rk8,0.540e-4_rk8, &
    2.0_rk8, 792.90_rk8,276.00_rk8,0.9971e+3_rk8,0.420e+1_rk8,0.560e-4_rk8, &
    3.0_rk8, 700.00_rk8,271.00_rk8,0.8985e+3_rk8,0.270e+1_rk8,0.580e-4_rk8, &
    4.0_rk8, 616.00_rk8,266.00_rk8,0.8077e+3_rk8,0.170e+1_rk8,0.600e-4_rk8, &
    5.0_rk8, 541.00_rk8,260.00_rk8,0.7244e+3_rk8,0.100e+1_rk8,0.640e-4_rk8, &
    6.0_rk8, 473.00_rk8,253.00_rk8,0.6519e+3_rk8,0.540e+0_rk8,0.710e-4_rk8, &
    7.0_rk8, 413.00_rk8,246.00_rk8,0.5849e+3_rk8,0.290e+0_rk8,0.750e-4_rk8, &
    8.0_rk8, 359.00_rk8,239.00_rk8,0.5231e+3_rk8,0.130e-1_rk8,0.790e-4_rk8, &
    9.0_rk8, 310.70_rk8,232.00_rk8,0.4663e+3_rk8,0.420e-1_rk8,0.110e-3_rk8, &
   10.0_rk8, 267.70_rk8,225.00_rk8,0.4142e+3_rk8,0.150e-1_rk8,0.130e-3_rk8, &
   11.0_rk8, 230.00_rk8,225.00_rk8,0.3559e+3_rk8,0.940e-2_rk8,0.180e-3_rk8, &
   12.0_rk8, 197.70_rk8,225.00_rk8,0.3059e+3_rk8,0.600e-2_rk8,0.210e-3_rk8, &
   13.0_rk8, 170.00_rk8,225.00_rk8,0.2630e+3_rk8,0.180e-2_rk8,0.260e-3_rk8, &
   14.0_rk8, 146.00_rk8,225.00_rk8,0.2260e+3_rk8,0.100e-2_rk8,0.280e-3_rk8, &
   15.0_rk8, 125.00_rk8,225.00_rk8,0.1943e+3_rk8,0.760e-3_rk8,0.320e-3_rk8, &
   16.0_rk8, 108.00_rk8,225.00_rk8,0.1671e+3_rk8,0.640e-3_rk8,0.340e-3_rk8, &
   17.0_rk8,  92.80_rk8,225.00_rk8,0.1436e+3_rk8,0.560e-3_rk8,0.390e-3_rk8, &
   18.0_rk8,  79.80_rk8,225.00_rk8,0.1235e+3_rk8,0.500e-3_rk8,0.410e-3_rk8, &
   19.0_rk8,  68.60_rk8,225.00_rk8,0.1062e+3_rk8,0.490e-3_rk8,0.410e-3_rk8, &
   20.0_rk8,  58.90_rk8,225.00_rk8,0.9128e+2_rk8,0.450e-3_rk8,0.390e-3_rk8, &
   21.0_rk8,  50.70_rk8,225.00_rk8,0.7849e+2_rk8,0.510e-3_rk8,0.360e-3_rk8, &
   22.0_rk8,  43.60_rk8,225.00_rk8,0.6750e+2_rk8,0.510e-3_rk8,0.320e-3_rk8, &
   23.0_rk8,  37.50_rk8,225.00_rk8,0.5805e+2_rk8,0.540e-3_rk8,0.300e-3_rk8, &
   24.0_rk8,  32.27_rk8,226.00_rk8,0.4963e+2_rk8,0.600e-3_rk8,0.280e-3_rk8, &
   25.0_rk8,  27.80_rk8,228.00_rk8,0.4247e+2_rk8,0.670e-3_rk8,0.260e-3_rk8, &
   30.0_rk8,  13.40_rk8,235.00_rk8,0.1338e+2_rk8,0.360e-3_rk8,0.140e-3_rk8, &
   35.0_rk8,   6.61_rk8,247.00_rk8,0.6614e+1_rk8,0.110e-3_rk8,0.920e-4_rk8, &
   40.0_rk8,   3.40_rk8,262.00_rk8,0.3404e+1_rk8,0.430e-4_rk8,0.410e-4_rk8, &
   45.0_rk8,   1.81_rk8,274.00_rk8,0.1817e+1_rk8,0.190e-4_rk8,0.130e-4_rk8, &
   50.0_rk8,   0.99_rk8,277.00_rk8,0.9868e+0_rk8,0.630e-5_rk8,0.430e-5_rk8, &
    0.0_rk8,1013.00_rk8,257.10_rk8,0.1372e+4_rk8,0.120e+1_rk8,0.410e-4_rk8, &
    1.0_rk8, 887.80_rk8,259.10_rk8,0.1193e+4_rk8,0.120e+1_rk8,0.410e-4_rk8, &
    2.0_rk8, 777.50_rk8,255.90_rk8,0.1058e+4_rk8,0.940e+0_rk8,0.410e-4_rk8, &
    3.0_rk8, 679.80_rk8,252.70_rk8,0.9366e+3_rk8,0.680e+0_rk8,0.430e-4_rk8, &
    4.0_rk8, 593.20_rk8,247.70_rk8,0.8339e+3_rk8,0.410e+0_rk8,0.450e-4_rk8, &
    5.0_rk8, 515.80_rk8,240.90_rk8,0.7457e+3_rk8,0.200e+0_rk8,0.470e-4_rk8, &
    6.0_rk8, 446.70_rk8,234.10_rk8,0.6646e+3_rk8,0.980e-1_rk8,0.490e-4_rk8, &
    7.0_rk8, 385.30_rk8,227.30_rk8,0.5904e+3_rk8,0.540e-1_rk8,0.710e-4_rk8, &
    8.0_rk8, 330.80_rk8,220.60_rk8,0.5226e+3_rk8,0.110e-1_rk8,0.900e-4_rk8, &
    9.0_rk8, 282.90_rk8,217.20_rk8,0.4538e+3_rk8,0.840e-2_rk8,0.160e-3_rk8, &
   10.0_rk8, 241.80_rk8,217.20_rk8,0.3879e+3_rk8,0.550e-2_rk8,0.240e-3_rk8, &
   11.0_rk8, 206.70_rk8,217.20_rk8,0.3315e+3_rk8,0.380e-2_rk8,0.320e-3_rk8, &
   12.0_rk8, 176.60_rk8,217.20_rk8,0.2834e+3_rk8,0.260e-2_rk8,0.430e-3_rk8, &
   13.0_rk8, 151.00_rk8,217.20_rk8,0.2422e+3_rk8,0.180e-2_rk8,0.470e-3_rk8, &
   14.0_rk8, 129.10_rk8,217.20_rk8,0.2071e+3_rk8,0.100e-2_rk8,0.490e-3_rk8, &
   15.0_rk8, 110.30_rk8,217.20_rk8,0.1770e+3_rk8,0.760e-3_rk8,0.560e-3_rk8, &
   16.0_rk8,  94.31_rk8,216.60_rk8,0.1517e+3_rk8,0.640e-3_rk8,0.620e-3_rk8, &
   17.0_rk8,  80.58_rk8,216.00_rk8,0.1300e+3_rk8,0.560e-3_rk8,0.620e-3_rk8, &
   18.0_rk8,  68.82_rk8,215.40_rk8,0.1113e+3_rk8,0.500e-3_rk8,0.620e-3_rk8, &
   19.0_rk8,  58.75_rk8,214.80_rk8,0.9529e+2_rk8,0.490e-3_rk8,0.600e-3_rk8, &
   20.0_rk8,  50.14_rk8,214.10_rk8,0.8155e+2_rk8,0.450e-3_rk8,0.560e-3_rk8, &
   21.0_rk8,  42.77_rk8,213.60_rk8,0.6976e+2_rk8,0.510e-3_rk8,0.510e-3_rk8, &
   22.0_rk8,  36.47_rk8,213.00_rk8,0.5966e+2_rk8,0.510e-3_rk8,0.470e-3_rk8, &
   23.0_rk8,  31.09_rk8,212.40_rk8,0.5100e+2_rk8,0.540e-3_rk8,0.430e-3_rk8, &
   24.0_rk8,  26.49_rk8,211.80_rk8,0.4358e+2_rk8,0.600e-3_rk8,0.360e-3_rk8, &
   25.0_rk8,  22.56_rk8,211.20_rk8,0.3722e+2_rk8,0.670e-3_rk8,0.320e-3_rk8, &
   30.0_rk8,  10.20_rk8,216.00_rk8,0.1645e+2_rk8,0.360e-3_rk8,0.150e-3_rk8, &
   35.0_rk8,   4.70_rk8,222.20_rk8,0.7368e+1_rk8,0.110e-3_rk8,0.920e-4_rk8, &
   40.0_rk8,   2.24_rk8,234.70_rk8,0.3330e+1_rk8,0.430e-4_rk8,0.410e-4_rk8, &
   45.0_rk8,   1.11_rk8,247.00_rk8,0.1569e+1_rk8,0.190e-4_rk8,0.130e-4_rk8, &
   50.0_rk8,   0.57_rk8,259.00_rk8,0.7682e+0_rk8,0.630e-5_rk8,0.430e-5_rk8 ], &
   [n_atmparms,n_atmlevls,n_atmzones])

 !
 ! For OpenACC:
 ! "declare create" is a data region defined as having the same scope as the
 ! scoping unit in which it's used.  Primarily used for module and global
 ! variables. If using module data directly (i.e. not passed in as an argument)
 ! in a device subroutine, it's required.
 !

!!!$acc declare create(stdhlevh,stdhlevf,stdplevh,stdplevf,stdatm)

  interface stdatm_val
    module procedure stdatm_val_seasonal_r4
    module procedure stdatm_val_noseason_r4
    module procedure stdatm_val_seasonal_r8
    module procedure stdatm_val_noseason_r8
  end interface stdatm_val

  interface stdlrate
    module procedure stdlrate_r4
    module procedure stdlrate_r8
  end interface stdlrate

  public :: stdplevf, stdplevh, stdhlevh, stdhlevf, stdatm_val, stdlrate

!-------------------------------------------------------------------------------
!
!  Standard TROpical ATMosphere
!
!  Standard Mid-Latitudes Summer ATMosphere
!
!  Standard Mid-Latitudes Winter ATMosphere
!
!  Standard POlar Summer ATMosphere
!
!  Standard POlar Winter ATMosphere
!
   contains

     pure real(rk4) function stdatm_val_noseason_r4(lat,plev,ival)
!$acc routine seq
       implicit none
       real(rk4), intent(in) :: lat
       real(rk4), intent(in) :: plev
       integer(ik4), intent(in) :: ival
       stdatm_val_noseason_r4 = real(stdatm_val_noseason_r8( &
         real(lat,rk8),real(plev,rk8),ival),rk4)
     end function stdatm_val_noseason_r4

     pure real(rk8) function stdatm_val_noseason_r8(lat,plev,ival)
!$acc routine seq
       implicit none
       real(rk8), intent(in) :: lat
       real(rk8), intent(in) :: plev
       integer(ik4), intent(in) :: ival
       integer(ik4) :: kp1, kp2
       real(rk8) :: wtp1, wtp2, wtl2, wtl1, wts1, wts2
       real(rk8) :: vs1, vs2

       wts1 = 0.5_rk8
       wts2 = 1.0_rk8-wts1
       if ( abs(lat) >= 45.0_rk8 ) then
         wtl1 = ((90.0_rk8-abs(lat))/45.0_rk8)**3
         wtl2 = 1.0_rk8-wtl1
         kp1 = find_klev(plev,ipolarwinter)
         kp2 = find_klev(plev,ipolarsummer)
         wtp1 = plev_wgt(kp1,plev,ipolarwinter)
         wtp2 = plev_wgt(kp1,plev,ipolarsummer)
         vs1 = stdatm(ival,kp1,ipolarwinter)*wtp1 + &
               stdatm(ival,kp1+1,ipolarwinter)*(1.0_rk8-wtp1)
         vs2 = stdatm(ival,kp1,ipolarsummer)*wtp1 + &
               stdatm(ival,kp1+1,ipolarsummer)*(1.0_rk8-wtp1)
         stdatm_val_noseason_r8 = vs1*wts1+vs2*wts2
         kp1 = find_klev(plev,imidlatwinter)
         kp2 = find_klev(plev,imidlatsummer)
         wtp1 = plev_wgt(kp1,plev,imidlatwinter)
         wtp2 = plev_wgt(kp1,plev,imidlatsummer)
         vs1 = stdatm(ival,kp1,imidlatwinter)*wtp1 + &
               stdatm(ival,kp1+1,imidlatwinter)*(1.0_rk8-wtp1)
         vs2 = stdatm(ival,kp1,imidlatsummer)*wtp1 + &
               stdatm(ival,kp1+1,imidlatsummer)*(1.0_rk8-wtp1)
         stdatm_val_noseason_r8 = wtl2*stdatm_val_noseason_r8 + &
                                  wtl1*(vs1*wts1+vs2*wts2)
       else
         wtl1 = (abs(lat)/45.0_rk8)**2
         wtl2 = 1.0_rk8-wtl1
         kp1 = find_klev(plev,itropical)
         wtp1 = plev_wgt(kp1,plev,itropical)
         stdatm_val_noseason_r8 = stdatm(ival,kp1,itropical)*wtp1 + &
                      stdatm(ival,kp1+1,itropical)*(1.0_rk8-wtp1)
         kp1 = find_klev(plev,imidlatwinter)
         kp2 = find_klev(plev,imidlatsummer)
         wtp1 = plev_wgt(kp1,plev,imidlatwinter)
         wtp2 = plev_wgt(kp1,plev,imidlatsummer)
         vs1 = stdatm(ival,kp1,imidlatwinter)*wtp1 + &
               stdatm(ival,kp1+1,imidlatwinter)*(1.0_rk8-wtp1)
         vs2 = stdatm(ival,kp2,imidlatsummer)*wtp2 + &
               stdatm(ival,kp2+1,imidlatsummer)*(1.0_rk8-wtp2)
         stdatm_val_noseason_r8 = wtl2*stdatm_val_noseason_r8 + &
                                  wtl1*(vs1*wts1+vs2*wts2)
       end if
     end function stdatm_val_noseason_r8

     pure real(rk4) function stdatm_val_seasonal_r4(jday,dayspy,lat,plev,ival)
!$acc routine seq
       implicit none
       real(rk4), intent(in) :: jday
       real(rk4), intent(in) :: dayspy
       real(rk4), intent(in) :: lat
       real(rk4), intent(in) :: plev
       integer(ik4), intent(in) :: ival
       stdatm_val_seasonal_r4 = real(stdatm_val_seasonal_r8( &
         real(jday,rk8), real(dayspy,rk8), &
         real(lat,rk8), real(plev,rk8),ival), rk4)
     end function stdatm_val_seasonal_r4

     pure real(rk8) function stdatm_val_seasonal_r8(jday,dayspy,lat,plev,ival)
!$acc routine seq
       implicit none
       real(rk8), intent(in) :: jday
       real(rk8), intent(in) :: dayspy
       real(rk8), intent(in) :: lat
       real(rk8), intent(in) :: plev
       integer(ik4), intent(in) :: ival
       integer(ik4) :: kp1, kp2
       real(rk8) :: wtl1, wtl2, wts1, wts2, wtp1, wtp2
       real(rk8) :: vs1, vs2

       stdatm_val_seasonal_r8 = dmissval

       if ( abs(lat) >= 45.0_rk8 ) then
         wtl1 = ((90.0_rk8-abs(lat))/45.0_rk8)**3
         wtl2 = 1.0_rk8-wtl1
         wts1 = winter_wgt(jday,dayspy)
         wts2 = 1.0_rk8-wts1
         kp1 = find_klev(plev,ipolarwinter)
         kp2 = find_klev(plev,ipolarsummer)
         wtp1 = plev_wgt(kp1,plev,ipolarwinter)
         wtp2 = plev_wgt(kp1,plev,ipolarsummer)
         vs1 = stdatm(ival,kp1,ipolarwinter)*wtp1 + &
               stdatm(ival,kp1+1,ipolarwinter)*(1.0_rk8-wtp1)
         vs2 = stdatm(ival,kp1,ipolarsummer)*wtp2 + &
               stdatm(ival,kp1+1,ipolarsummer)*(1.0_rk8-wtp2)
         stdatm_val_seasonal_r8 = vs1*wts1+vs2*wts2
         kp1 = find_klev(plev,imidlatwinter)
         kp2 = find_klev(plev,imidlatsummer)
         wtp1 = plev_wgt(kp1,plev,imidlatwinter)
         wtp2 = plev_wgt(kp1,plev,imidlatsummer)
         vs1 = stdatm(ival,kp1,imidlatwinter)*wtp1 + &
               stdatm(ival,kp1+1,imidlatwinter)*(1.0_rk8-wtp1)
         vs2 = stdatm(ival,kp1,imidlatsummer)*wtp2 + &
               stdatm(ival,kp1+1,imidlatsummer)*(1.0_rk8-wtp2)
         stdatm_val_seasonal_r8 = wtl2*stdatm_val_seasonal_r8 + &
                                  wtl1*(vs1*wts1+vs2*wts2)
       else
         wtl1 = (abs(lat)/45.0_rk8)**2
         wtl2 = 1.0_rk8-wtl1
         wts1 = winter_wgt(jday,dayspy)
         wts2 = 1.0_rk8-wts1
         kp1 = find_klev(plev,itropical)
         wtp1 = plev_wgt(kp1,plev,itropical)
         stdatm_val_seasonal_r8 = stdatm(ival,kp1,itropical)*wtp1 + &
                      stdatm(ival,kp1+1,itropical)*(1.0_rk8-wtp1)
         kp1 = find_klev(plev,imidlatwinter)
         kp2 = find_klev(plev,imidlatsummer)
         wtp1 = plev_wgt(kp1,plev,imidlatwinter)
         wtp2 = plev_wgt(kp1,plev,imidlatsummer)
         vs1 = stdatm(ival,kp1,imidlatwinter)*wtp1 + &
               stdatm(ival,kp1+1,imidlatwinter)*(1.0_rk8-wtp1)
         vs2 = stdatm(ival,kp2,imidlatsummer)*wtp2 + &
               stdatm(ival,kp2+1,imidlatsummer)*(1.0_rk8-wtp2)
         stdatm_val_seasonal_r8 = wtl2*stdatm_val_seasonal_r8 + &
                                  wtl1*(vs1*wts1+vs2*wts2)
       end if
     end function stdatm_val_seasonal_r8

     pure real(rk8) function winter_wgt(jday,dayspy)
!$acc routine seq
       implicit none
       real(rk8), intent(in) :: jday, dayspy
       real(rk8) :: dis, half_dayspy, sixteenth_dayspy
       half_dayspy = dayspy * 0.5_rk8
       sixteenth_dayspy = dayspy * 0.0625_rk8
       dis = ((half_dayspy-jday-sixteenth_dayspy+1.0_rk8)/dayspy)*mathpi
       winter_wgt = real(sin(dis)**2,rk8)
     end function winter_wgt

     pure integer(ik4) function find_klev(plev,izone)
!$acc routine seq
       implicit none
       integer(ik4), intent(in) :: izone
       real(rk8), intent(in) :: plev
       integer(ik4) :: k
       find_klev = 1
       do k = 2, n_atmlevls
         find_klev = k-1
         if ( plev > stdatm(istdatm_prsmb,k,izone) ) exit
       end do
     end function find_klev

     pure integer(ik4) function find_zlev(z,izone)
!$acc routine seq
       implicit none
       integer(ik4), intent(in) :: izone
       real(rk8), intent(in) :: z
       integer(ik4) :: k
       find_zlev = 1
       do k = 2, n_atmlevls
         find_zlev = k-1
         if ( z < 1000.0_rk8 * stdatm(istdatm_hgtkm,k,izone) ) exit
       end do
     end function find_zlev

     pure real(rk8) function plev_wgt(k,plev,izone)
!$acc routine seq
       implicit none
       integer(ik4), intent(in) :: k, izone
       real(rk8), intent(in) :: plev
       if ( plev >= stdatm(istdatm_prsmb,k,izone) ) then
         plev_wgt = 1.0_rk8
       else if ( plev <= stdatm(istdatm_prsmb,k+1,izone) ) then
         plev_wgt = 0.0_rk8
       else
         plev_wgt = log(plev/stdatm(istdatm_prsmb,k+1,izone)) / &
             log(stdatm(istdatm_prsmb,k,izone)/stdatm(istdatm_prsmb,k+1,izone))
       end if
     end function plev_wgt

     pure real(rk8) function stdlrate_r4(jday,dayspy,lat)
!$acc routine seq
       implicit none
       real(rk4), intent(in) :: jday
       real(rk4), intent(in) :: dayspy
       real(rk4), intent(in) :: lat
       stdlrate_r4 = real(stdlrate_r8( &
         real(jday,rk8),real(dayspy,rk8),real(lat,rk8)),rk4)
     end function stdlrate_r4

     pure real(rk8) function stdlrate_r8(jday,dayspy,lat)
!$acc routine seq
       implicit none
       real(rk8), intent(in) :: jday
       real(rk8), intent(in) :: dayspy
       real(rk8), intent(in) :: lat
       real(rk8) :: wts1, wts2, wtp1, wtp2, vs1, vs2, vz1, vz2
       real(rk8) :: wtl1, wtl2, vs3, vs4
       integer(ik4) :: kp1, kp2
       real(rk8) :: polar, midlat, tropical
       if ( abs(lat) >= 45.0_rk8 ) then
         wtl1 = ((90.0_rk8-abs(lat))/45.0_rk8)**3
         wtl2 = 1.0_rk8-wtl1
         wts1 = winter_wgt(jday,dayspy)
         wts2 = 1.0_rk8-wts1
         kp1 = find_klev(300.0_rk8,ipolarwinter)
         kp2 = find_klev(300.0_rk8,ipolarsummer)
         wtp1 = plev_wgt(kp1,300.0_rk8,ipolarwinter)
         wtp2 = plev_wgt(kp1,300.0_rk8,ipolarsummer)
         vs1 = stdatm(istdatm_tempk,kp1,ipolarwinter)*wtp1 + &
               stdatm(istdatm_tempk,kp1+1,ipolarwinter)*(1.0_rk8-wtp1)
         vs2 = stdatm(istdatm_tempk,kp2,ipolarsummer)*wtp2 + &
               stdatm(istdatm_tempk,kp2+1,ipolarsummer)*(1.0_rk8-wtp2)
         vs3 = stdatm(istdatm_tempk,1,ipolarwinter)
         vs4 = stdatm(istdatm_tempk,1,ipolarsummer)
         vz1 = stdatm(istdatm_hgtkm,kp1,ipolarwinter)*wtp1 + &
               stdatm(istdatm_hgtkm,kp1+1,ipolarwinter)*(1.0_rk8-wtp1)
         vz2 = stdatm(istdatm_hgtkm,kp2,ipolarsummer)*wtp2 + &
               stdatm(istdatm_hgtkm,kp2+1,ipolarsummer)*(1.0_rk8-wtp2)
         polar = ((vs1*wts1+vs2*wts2) - (vs3*wts1+vs4*wts2)) / &
                 (d_1000*(vz1*wts1+vz2*wts2))
         kp1 = find_klev(300.0_rk8,imidlatwinter)
         kp2 = find_klev(300.0_rk8,imidlatsummer)
         wtp1 = plev_wgt(kp1,300.0_rk8,imidlatwinter)
         wtp2 = plev_wgt(kp1,300.0_rk8,imidlatsummer)
         vs1 = stdatm(istdatm_tempk,kp1,imidlatwinter)*wtp1 + &
               stdatm(istdatm_tempk,kp1+1,imidlatwinter)*(1.0_rk8-wtp1)
         vs2 = stdatm(istdatm_tempk,kp2,imidlatsummer)*wtp2 + &
               stdatm(istdatm_tempk,kp2+1,imidlatsummer)*(1.0_rk8-wtp2)
         vs3 = stdatm(istdatm_tempk,1,imidlatwinter)
         vs4 = stdatm(istdatm_tempk,1,imidlatsummer)
         vz1 = stdatm(istdatm_hgtkm,kp1,imidlatwinter)*wtp1 + &
               stdatm(istdatm_hgtkm,kp1+1,imidlatwinter)*(1.0_rk8-wtp1)
         vz2 = stdatm(istdatm_hgtkm,kp2,imidlatsummer)*wtp2 + &
               stdatm(istdatm_hgtkm,kp2+1,imidlatsummer)*(1.0_rk8-wtp2)
         midlat = ((vs1*wts1+vs2*wts2) - (vs3*wts1+vs4*wts2)) / &
                  (d_1000*(vz1*wts1+vz2*wts2))
         stdlrate_r8 = wtl2 * polar + wtl1 * midlat
       else
         wtl1 = (abs(lat)/45.0_rk8)**2
         wtl2 = 1.0_rk8-wtl1
         kp1 = find_klev(300.0_rk8,itropical)
         wtp1 = plev_wgt(kp1,300.0_rk8,itropical)
         vs1 = stdatm(istdatm_tempk,kp1,itropical)*wtp1 + &
               stdatm(istdatm_tempk,kp1+1,itropical)*(1.0_rk8-wtp1)
         vs2 = stdatm(istdatm_tempk,1,itropical)
         vz1 = stdatm(istdatm_hgtkm,kp1,itropical)*wtp1 + &
               stdatm(istdatm_hgtkm,kp1+1,itropical)*(1.0_rk8-wtp1)
         tropical = (vs1 - vs2)/(d_1000*vz1)
         wts1 = winter_wgt(jday,dayspy)
         wts2 = 1.0_rk8-wts1
         kp1 = find_klev(300.0_rk8,imidlatwinter)
         kp2 = find_klev(300.0_rk8,imidlatsummer)
         wtp1 = plev_wgt(kp1,300.0_rk8,imidlatwinter)
         wtp2 = plev_wgt(kp1,300.0_rk8,imidlatsummer)
         vs1 = stdatm(istdatm_tempk,kp1,imidlatwinter)*wtp1 + &
               stdatm(istdatm_tempk,kp1+1,imidlatwinter)*(1.0_rk8-wtp1)
         vs2 = stdatm(istdatm_tempk,kp1,imidlatsummer)*wtp2 + &
               stdatm(istdatm_tempk,kp1+1,imidlatsummer)*(1.0_rk8-wtp2)
         vs3 = stdatm(istdatm_tempk,1,imidlatwinter)
         vs4 = stdatm(istdatm_tempk,1,imidlatsummer)
         vz1 = stdatm(istdatm_hgtkm,kp1,imidlatwinter)*wtp1 + &
               stdatm(istdatm_hgtkm,kp1+1,imidlatwinter)*(1.0_rk8-wtp1)
         vz2 = stdatm(istdatm_hgtkm,kp1,imidlatsummer)*wtp2 + &
               stdatm(istdatm_hgtkm,kp1+1,imidlatsummer)*(1.0_rk8-wtp2)
         midlat = ((vs1*wts1+vs2*wts2) - (vs3*wts1+vs4*wts2)) / &
                  (d_1000*(vz1*wts1+vz2*wts2))
         stdlrate_r8 = wtl2 * midlat + wtl1 * tropical
       end if
    end function stdlrate_r8

end module mod_stdatm
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

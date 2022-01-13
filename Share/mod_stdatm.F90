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

module mod_stdatm

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam

  implicit none

  private

  integer(ik4) , public , parameter :: itropical = 1
  integer(ik4) , public , parameter :: imidlatsummer = 2
  integer(ik4) , public , parameter :: imidlatwinter = 3
  integer(ik4) , public , parameter :: ipolarsummer = 4
  integer(ik4) , public , parameter :: ipolarwinter = 5

  integer(ik4) , public , parameter :: istdatm_hgtkm = 1  ! HGT in km
  integer(ik4) , public , parameter :: istdatm_prsmb = 2  ! PRESS in hPa
  integer(ik4) , public , parameter :: istdatm_tempk = 3  ! TEMP in K
  integer(ik4) , public , parameter :: istdatm_airdn = 4  ! RHO in kg/m^3
  integer(ik4) , public , parameter :: istdatm_qdens = 5  ! Q in kg/m^3
  integer(ik4) , public , parameter :: istdatm_ozone = 6  ! Ozone in kg/m^3

  integer(ik4) , public , parameter :: n_atmzones = 5
  integer(ik4) , public , parameter :: n_atmparms = 6
  integer(ik4) , public , parameter :: n_atmlevls = 31
  integer(ik4) , public , parameter :: n_preflev  = 31
  integer(ik4) , public , parameter :: n_prehlev  = 30
  integer(ik4) , public , parameter :: n_hreflev  = 31
  integer(ik4) , public , parameter :: n_hrehlev  = 30
  integer(ik4) :: ip , il , iz
  real(rkx) , dimension(n_atmparms,n_atmlevls,n_atmzones) :: stdatm
  real(rkx) , dimension(n_hrehlev) :: stdhlevh
  real(rkx) , dimension(n_hreflev) :: stdhlevf
  real(rkx) , dimension(n_prehlev) :: stdplevh
  real(rkx) , dimension(n_preflev) :: stdplevf

  interface stdatm_val
    module procedure stdatm_val_seasonal
    module procedure stdatm_val_noseason
  end interface stdatm_val

  public :: stdplevf , stdplevh , stdhlevh , stdhlevf , stdatm_val , stdlrate

!-------------------------------------------------------------------------------
!
!  Standard reference pressure levels
!
  data stdplevf / &
    1013.0_rkx , 900.0_rkx , 800.0_rkx , 700.0_rkx , 650.0_rkx , 550.0_rkx , &
     500.0_rkx , 450.0_rkx , 375.0_rkx , 325.0_rkx , 275.0_rkx , 250.0_rkx , &
     225.0_rkx , 175.0_rkx , 150.0_rkx , 125.0_rkx , 100.0_rkx ,  90.0_rkx , &
      75.0_rkx ,  65.0_rkx ,  55.0_rkx ,  45.0_rkx ,  40.0_rkx ,  35.0_rkx , &
      30.0_rkx ,  25.0_rkx ,  12.5_rkx ,   7.5_rkx ,   2.5_rkx ,   1.5_rkx , &
       0.5_rkx /
  data stdplevh / &
     950.0_rkx , 850.0_rkx , 750.0_rkx , 675.0_rkx , 600.0_rkx , 525.0_rkx , &
     475.0_rkx , 415.0_rkx , 350.0_rkx , 300.0_rkx , 265.0_rkx , 215.0_rkx , &
     200.0_rkx , 165.0_rkx , 135.0_rkx , 115.0_rkx ,  95.0_rkx ,  85.0_rkx , &
      70.0_rkx ,  60.0_rkx ,  50.0_rkx ,  42.5_rkx ,  37.5_rkx ,  32.5_rkx , &
      27.5_rkx ,  17.5_rkx ,  10.0_rkx ,   5.0_rkx ,   2.0_rkx ,   1.0_rkx /
  data stdhlevf / &
       0.0_rkx ,  1.0_rkx ,  2.0_rkx ,  3.0_rkx ,  4.0_rkx ,  5.0_rkx , &
       6.0_rkx ,  7.0_rkx ,  8.0_rkx ,  9.0_rkx , 10.0_rkx , 11.0_rkx , &
      12.0_rkx , 13.0_rkx , 14.0_rkx , 15.0_rkx , 16.0_rkx , 17.0_rkx , &
      18.0_rkx , 19.0_rkx , 20.0_rkx , 21.0_rkx , 22.0_rkx , 23.0_rkx , &
      24.0_rkx , 25.0_rkx , 30.0_rkx , 35.0_rkx , 40.0_rkx , 45.0_rkx , &
      50.0_rkx /
  data stdhlevh / &
       0.5_rkx ,  1.5_rkx ,  2.5_rkx ,  3.5_rkx ,  4.5_rkx ,  5.5_rkx , &
       6.5_rkx ,  7.5_rkx ,  8.5_rkx ,  9.5_rkx , 10.5_rkx , 11.5_rkx , &
      12.5_rkx , 13.5_rkx , 14.5_rkx , 15.5_rkx , 16.5_rkx , 17.5_rkx , &
      18.5_rkx , 19.5_rkx , 20.5_rkx , 21.5_rkx , 22.5_rkx , 23.5_rkx , &
      24.5_rkx , 27.5_rkx , 32.5_rkx , 37.5_rkx , 42.5_rkx , 47.5_rkx /
!
!  Standard TROpical ATMosphere
!
  data (((stdatm(ip,il,iz),ip=1,6),il=1,31),iz=1,1) &
  / 0.0_rkx,1013.00_rkx,300.00_rkx,0.116700e+4_rkx,0.190000e+2_rkx,0.560000e-4_rkx , &
    1.0_rkx, 904.00_rkx,294.00_rkx,0.106400e+4_rkx,0.130000e+2_rkx,0.560000e-4_rkx , &
    2.0_rkx, 805.00_rkx,288.00_rkx,0.968900e+3_rkx,0.930000e+1_rkx,0.540000e-4_rkx , &
    3.0_rkx, 715.00_rkx,284.00_rkx,0.875600e+3_rkx,0.470000e+1_rkx,0.510000e-4_rkx , &
    4.0_rkx, 633.00_rkx,277.00_rkx,0.795100e+3_rkx,0.220000e+1_rkx,0.470000e-4_rkx , &
    5.0_rkx, 559.00_rkx,270.00_rkx,0.719900e+3_rkx,0.150000e+1_rkx,0.450000e-4_rkx , &
    6.0_rkx, 492.00_rkx,264.00_rkx,0.650100e+3_rkx,0.850000e+0_rkx,0.430000e-4_rkx , &
    7.0_rkx, 432.00_rkx,257.00_rkx,0.585500e+3_rkx,0.470000e+0_rkx,0.410000e-4_rkx , &
    8.0_rkx, 378.00_rkx,250.00_rkx,0.525800e+3_rkx,0.250000e+0_rkx,0.390000e-4_rkx , &
    9.0_rkx, 329.00_rkx,244.00_rkx,0.470800e+3_rkx,0.120000e+0_rkx,0.390000e-4_rkx , &
   10.0_rkx, 286.00_rkx,237.00_rkx,0.420200e+3_rkx,0.500000e-1_rkx,0.390000e-4_rkx , &
   11.0_rkx, 247.00_rkx,230.00_rkx,0.374000e+3_rkx,0.170000e-1_rkx,0.410000e-4_rkx , &
   12.0_rkx, 213.00_rkx,224.00_rkx,0.331600e+3_rkx,0.600000e-2_rkx,0.430000e-4_rkx , &
   13.0_rkx, 182.00_rkx,217.00_rkx,0.292900e+3_rkx,0.180000e-2_rkx,0.450000e-4_rkx , &
   14.0_rkx, 156.00_rkx,210.00_rkx,0.257800e+3_rkx,0.100000e-2_rkx,0.450000e-4_rkx , &
   15.0_rkx, 132.00_rkx,204.00_rkx,0.226000e+3_rkx,0.760000e-3_rkx,0.470000e-4_rkx , &
   16.0_rkx, 111.00_rkx,197.00_rkx,0.197200e+3_rkx,0.640000e-3_rkx,0.470000e-4_rkx , &
   17.0_rkx,  93.70_rkx,195.00_rkx,0.167600e+3_rkx,0.560000e-3_rkx,0.690000e-4_rkx , &
   18.0_rkx,  78.90_rkx,199.00_rkx,0.138200e+3_rkx,0.500000e-3_rkx,0.900000e-4_rkx , &
   19.0_rkx,  66.60_rkx,203.00_rkx,0.114500e+3_rkx,0.490000e-3_rkx,0.140000e-3_rkx , &
   20.0_rkx,  56.50_rkx,207.00_rkx,0.951500e+2_rkx,0.450000e-3_rkx,0.190000e-3_rkx , &
   21.0_rkx,  48.00_rkx,211.00_rkx,0.793800e+2_rkx,0.510000e-3_rkx,0.240000e-3_rkx , &
   22.0_rkx,  40.90_rkx,215.00_rkx,0.664500e+2_rkx,0.510000e-3_rkx,0.280000e-3_rkx , &
   23.0_rkx,  35.00_rkx,217.00_rkx,0.561800e+2_rkx,0.540000e-3_rkx,0.320000e-3_rkx , &
   24.0_rkx,  30.00_rkx,219.00_rkx,0.476300e+2_rkx,0.600000e-3_rkx,0.340000e-3_rkx , &
   25.0_rkx,  25.70_rkx,221.00_rkx,0.404500e+2_rkx,0.670000e-3_rkx,0.340000e-3_rkx , &
   30.0_rkx,  12.20_rkx,232.00_rkx,0.183100e+2_rkx,0.360000e-3_rkx,0.240000e-3_rkx , &
   35.0_rkx,   6.00_rkx,243.00_rkx,0.860000e+1_rkx,0.110000e-3_rkx,0.920000e-4_rkx , &
   40.0_rkx,   3.05_rkx,254.00_rkx,0.418100e+1_rkx,0.430000e-4_rkx,0.410000e-4_rkx , &
   45.0_rkx,   1.59_rkx,265.00_rkx,0.209700e+1_rkx,0.190000e-4_rkx,0.130000e-4_rkx , &
   50.0_rkx,   0.85_rkx,270.00_rkx,0.110100e+1_rkx,0.630000e-5_rkx,0.430000e-5_rkx /
!
!  Standard Mid-Latitudes Summer ATMosphere
!
  data (((stdatm(ip,il,iz),ip=1,6),il=1,31),iz=2,2) &
  / 0.0_rkx,1013.00_rkx,294.00_rkx,0.119100e+4_rkx,0.140000e+2_rkx,0.600000e-4_rkx , &
    1.0_rkx, 902.00_rkx,290.00_rkx,0.108000e+4_rkx,0.930000e+1_rkx,0.600000e-4_rkx , &
    2.0_rkx, 802.00_rkx,285.00_rkx,0.975700e+3_rkx,0.590000e+1_rkx,0.600000e-4_rkx , &
    3.0_rkx, 710.00_rkx,279.00_rkx,0.884600e+3_rkx,0.330000e+1_rkx,0.620000e-4_rkx , &
    4.0_rkx, 628.00_rkx,273.00_rkx,0.799800e+3_rkx,0.190000e+1_rkx,0.640000e-4_rkx , &
    5.0_rkx, 554.00_rkx,267.00_rkx,0.721100e+3_rkx,0.100000e+1_rkx,0.660000e-4_rkx , &
    6.0_rkx, 487.00_rkx,261.00_rkx,0.648700e+3_rkx,0.610000e+0_rkx,0.690000e-4_rkx , &
    7.0_rkx, 426.00_rkx,255.00_rkx,0.583000e+3_rkx,0.370000e+0_rkx,0.750000e-4_rkx , &
    8.0_rkx, 372.00_rkx,248.00_rkx,0.522500e+3_rkx,0.210000e+0_rkx,0.790000e-4_rkx , &
    9.0_rkx, 324.00_rkx,242.00_rkx,0.466900e+3_rkx,0.120000e+0_rkx,0.860000e-4_rkx , &
   10.0_rkx, 281.00_rkx,235.00_rkx,0.415900e+3_rkx,0.640000e-1_rkx,0.900000e-4_rkx , &
   11.0_rkx, 243.00_rkx,229.00_rkx,0.369300e+3_rkx,0.220000e-1_rkx,0.110000e-3_rkx , &
   12.0_rkx, 209.00_rkx,222.00_rkx,0.326900e+3_rkx,0.600000e-2_rkx,0.120000e-3_rkx , &
   13.0_rkx, 179.00_rkx,216.00_rkx,0.288200e+3_rkx,0.180000e-2_rkx,0.150000e-3_rkx , &
   14.0_rkx, 153.00_rkx,216.00_rkx,0.246400e+3_rkx,0.100000e-2_rkx,0.180000e-3_rkx , &
   15.0_rkx, 130.00_rkx,216.00_rkx,0.210400e+3_rkx,0.760000e-3_rkx,0.190000e-3_rkx , &
   16.0_rkx, 111.00_rkx,216.00_rkx,0.179700e+3_rkx,0.640000e-3_rkx,0.210000e-3_rkx , &
   17.0_rkx,  95.00_rkx,216.00_rkx,0.153500e+3_rkx,0.560000e-3_rkx,0.240000e-3_rkx , &
   18.0_rkx,  81.20_rkx,216.00_rkx,0.130500e+3_rkx,0.500000e-3_rkx,0.280000e-3_rkx , &
   19.0_rkx,  69.50_rkx,217.00_rkx,0.111000e+3_rkx,0.490000e-3_rkx,0.320000e-3_rkx , &
   20.0_rkx,  59.50_rkx,218.00_rkx,0.945300e+2_rkx,0.450000e-3_rkx,0.340000e-3_rkx , &
   21.0_rkx,  51.00_rkx,219.00_rkx,0.805600e+2_rkx,0.510000e-3_rkx,0.360000e-3_rkx , &
   22.0_rkx,  43.70_rkx,220.00_rkx,0.687200e+2_rkx,0.510000e-3_rkx,0.360000e-3_rkx , &
   23.0_rkx,  37.60_rkx,222.00_rkx,0.586700e+2_rkx,0.540000e-3_rkx,0.340000e-3_rkx , &
   24.0_rkx,  32.20_rkx,223.00_rkx,0.501400e+2_rkx,0.600000e-3_rkx,0.320000e-3_rkx , &
   25.0_rkx,  27.70_rkx,224.00_rkx,0.428800e+2_rkx,0.670000e-3_rkx,0.300000e-3_rkx , &
   30.0_rkx,  13.20_rkx,234.00_rkx,0.132200e+2_rkx,0.360000e-3_rkx,0.200000e-3_rkx , &
   35.0_rkx,   6.52_rkx,245.00_rkx,0.651900e+1_rkx,0.110000e-3_rkx,0.920000e-4_rkx , &
   40.0_rkx,   3.33_rkx,258.00_rkx,0.333000e+1_rkx,0.430000e-4_rkx,0.410000e-4_rkx , &
   45.0_rkx,   1.76_rkx,270.60_rkx,0.175700e+1_rkx,0.190000e-4_rkx,0.130000e-4_rkx , &
   50.0_rkx,   0.95_rkx,276.00_rkx,0.951200e+0_rkx,0.630000e-5_rkx,0.430000e-5_rkx /
!
!  Standard Mid-Latitudes Winter ATMosphere
!
  data (((stdatm(ip,il,iz),ip=1,6),il=1,31),iz=3,3) &
  / 0.0_rkx,1018.00_rkx,272.20_rkx,0.130100e+4_rkx,0.350000e+1_rkx,0.600000e-4_rkx , &
    1.0_rkx, 897.30_rkx,268.70_rkx,0.116200e+4_rkx,0.250000e+1_rkx,0.540000e-4_rkx , &
    2.0_rkx, 789.70_rkx,265.20_rkx,0.103700e+4_rkx,0.180000e+1_rkx,0.490000e-4_rkx , &
    3.0_rkx, 693.80_rkx,261.70_rkx,0.923000e+3_rkx,0.120000e+1_rkx,0.490000e-4_rkx , &
    4.0_rkx, 608.10_rkx,255.70_rkx,0.828200e+3_rkx,0.660000e+0_rkx,0.490000e-4_rkx , &
    5.0_rkx, 531.30_rkx,249.70_rkx,0.741100e+3_rkx,0.380000e+0_rkx,0.580000e-4_rkx , &
    6.0_rkx, 462.70_rkx,243.70_rkx,0.661400e+3_rkx,0.210000e+0_rkx,0.640000e-4_rkx , &
    7.0_rkx, 401.60_rkx,237.70_rkx,0.588600e+3_rkx,0.850000e-1_rkx,0.770000e-4_rkx , &
    8.0_rkx, 347.30_rkx,231.70_rkx,0.522200e+3_rkx,0.350000e-1_rkx,0.900000e-4_rkx , &
    9.0_rkx, 299.20_rkx,225.70_rkx,0.461900e+3_rkx,0.160000e-1_rkx,0.120000e-3_rkx , &
   10.0_rkx, 256.80_rkx,219.70_rkx,0.407200e+3_rkx,0.750000e-2_rkx,0.160000e-3_rkx , &
   11.0_rkx, 219.90_rkx,219.20_rkx,0.349600e+3_rkx,0.690000e-2_rkx,0.210000e-3_rkx , &
   12.0_rkx, 188.20_rkx,218.70_rkx,0.299900e+3_rkx,0.600000e-2_rkx,0.260000e-3_rkx , &
   13.0_rkx, 161.00_rkx,218.20_rkx,0.257200e+3_rkx,0.180000e-2_rkx,0.300000e-3_rkx , &
   14.0_rkx, 137.80_rkx,217.70_rkx,0.220600e+3_rkx,0.100000e-2_rkx,0.320000e-3_rkx , &
   15.0_rkx, 117.80_rkx,217.20_rkx,0.189000e+3_rkx,0.760000e-3_rkx,0.340000e-3_rkx , &
   16.0_rkx, 100.70_rkx,216.70_rkx,0.162000e+3_rkx,0.640000e-3_rkx,0.360000e-3_rkx , &
   17.0_rkx,  86.10_rkx,216.20_rkx,0.138800e+3_rkx,0.560000e-3_rkx,0.390000e-3_rkx , &
   18.0_rkx,  73.50_rkx,215.70_rkx,0.118800e+3_rkx,0.500000e-3_rkx,0.410000e-3_rkx , &
   19.0_rkx,  62.80_rkx,215.20_rkx,0.101700e+3_rkx,0.490000e-3_rkx,0.430000e-3_rkx , &
   20.0_rkx,  53.70_rkx,215.20_rkx,0.869000e+2_rkx,0.450000e-3_rkx,0.450000e-3_rkx , &
   21.0_rkx,  45.80_rkx,215.20_rkx,0.742100e+2_rkx,0.510000e-3_rkx,0.430000e-3_rkx , &
   22.0_rkx,  39.10_rkx,215.20_rkx,0.633800e+2_rkx,0.510000e-3_rkx,0.430000e-3_rkx , &
   23.0_rkx,  33.40_rkx,215.20_rkx,0.541500e+2_rkx,0.540000e-3_rkx,0.390000e-3_rkx , &
   24.0_rkx,  28.60_rkx,215.20_rkx,0.462400e+2_rkx,0.600000e-3_rkx,0.360000e-3_rkx , &
   25.0_rkx,  24.30_rkx,215.20_rkx,0.395000e+2_rkx,0.670000e-3_rkx,0.340000e-3_rkx , &
   30.0_rkx,  11.10_rkx,217.40_rkx,0.178300e+2_rkx,0.360000e-3_rkx,0.190000e-3_rkx , &
   35.0_rkx,   5.18_rkx,227.80_rkx,0.792400e+1_rkx,0.110000e-3_rkx,0.920000e-4_rkx , &
   40.0_rkx,   2.53_rkx,243.20_rkx,0.362500e+1_rkx,0.430000e-4_rkx,0.410000e-4_rkx , &
   45.0_rkx,   1.29_rkx,258.50_rkx,0.174100e+1_rkx,0.190000e-4_rkx,0.130000e-4_rkx , &
   50.0_rkx,   0.68_rkx,265.70_rkx,0.895400e+0_rkx,0.630000e-5_rkx,0.430000e-5_rkx /
!
!  Standard POlar Summer ATMosphere
!
  data (((stdatm(ip,il,iz),ip=1,6),il=1,31),iz=4,4) &
  / 0.0_rkx,1010.00_rkx,287.00_rkx,0.122000e+4_rkx,0.910000e+1_rkx,0.490000e-4_rkx , &
    1.0_rkx, 896.00_rkx,282.00_rkx,0.111000e+4_rkx,0.600000e+1_rkx,0.540000e-4_rkx , &
    2.0_rkx, 792.90_rkx,276.00_rkx,0.997100e+3_rkx,0.420000e+1_rkx,0.560000e-4_rkx , &
    3.0_rkx, 700.00_rkx,271.00_rkx,0.898500e+3_rkx,0.270000e+1_rkx,0.580000e-4_rkx , &
    4.0_rkx, 616.00_rkx,266.00_rkx,0.807700e+3_rkx,0.170000e+1_rkx,0.600000e-4_rkx , &
    5.0_rkx, 541.00_rkx,260.00_rkx,0.724400e+3_rkx,0.100000e+1_rkx,0.640000e-4_rkx , &
    6.0_rkx, 473.00_rkx,253.00_rkx,0.651900e+3_rkx,0.540000e+0_rkx,0.710000e-4_rkx , &
    7.0_rkx, 413.00_rkx,246.00_rkx,0.584900e+3_rkx,0.290000e+0_rkx,0.750000e-4_rkx , &
    8.0_rkx, 359.00_rkx,239.00_rkx,0.523100e+3_rkx,0.130000e-1_rkx,0.790000e-4_rkx , &
    9.0_rkx, 310.70_rkx,232.00_rkx,0.466300e+3_rkx,0.420000e-1_rkx,0.110000e-3_rkx , &
   10.0_rkx, 267.70_rkx,225.00_rkx,0.414200e+3_rkx,0.150000e-1_rkx,0.130000e-3_rkx , &
   11.0_rkx, 230.00_rkx,225.00_rkx,0.355900e+3_rkx,0.940000e-2_rkx,0.180000e-3_rkx , &
   12.0_rkx, 197.70_rkx,225.00_rkx,0.305900e+3_rkx,0.600000e-2_rkx,0.210000e-3_rkx , &
   13.0_rkx, 170.00_rkx,225.00_rkx,0.263000e+3_rkx,0.180000e-2_rkx,0.260000e-3_rkx , &
   14.0_rkx, 146.00_rkx,225.00_rkx,0.226000e+3_rkx,0.100000e-2_rkx,0.280000e-3_rkx , &
   15.0_rkx, 125.00_rkx,225.00_rkx,0.194300e+3_rkx,0.760000e-3_rkx,0.320000e-3_rkx , &
   16.0_rkx, 108.00_rkx,225.00_rkx,0.167100e+3_rkx,0.640000e-3_rkx,0.340000e-3_rkx , &
   17.0_rkx,  92.80_rkx,225.00_rkx,0.143600e+3_rkx,0.560000e-3_rkx,0.390000e-3_rkx , &
   18.0_rkx,  79.80_rkx,225.00_rkx,0.123500e+3_rkx,0.500000e-3_rkx,0.410000e-3_rkx , &
   19.0_rkx,  68.60_rkx,225.00_rkx,0.106200e+3_rkx,0.490000e-3_rkx,0.410000e-3_rkx , &
   20.0_rkx,  58.90_rkx,225.00_rkx,0.912800e+2_rkx,0.450000e-3_rkx,0.390000e-3_rkx , &
   21.0_rkx,  50.70_rkx,225.00_rkx,0.784900e+2_rkx,0.510000e-3_rkx,0.360000e-3_rkx , &
   22.0_rkx,  43.60_rkx,225.00_rkx,0.675000e+2_rkx,0.510000e-3_rkx,0.320000e-3_rkx , &
   23.0_rkx,  37.50_rkx,225.00_rkx,0.580500e+2_rkx,0.540000e-3_rkx,0.300000e-3_rkx , &
   24.0_rkx,  32.27_rkx,226.00_rkx,0.496300e+2_rkx,0.600000e-3_rkx,0.280000e-3_rkx , &
   25.0_rkx,  27.80_rkx,228.00_rkx,0.424700e+2_rkx,0.670000e-3_rkx,0.260000e-3_rkx , &
   30.0_rkx,  13.40_rkx,235.00_rkx,0.133800e+2_rkx,0.360000e-3_rkx,0.140000e-3_rkx , &
   35.0_rkx,   6.61_rkx,247.00_rkx,0.661400e+1_rkx,0.110000e-3_rkx,0.920000e-4_rkx , &
   40.0_rkx,   3.40_rkx,262.00_rkx,0.340400e+1_rkx,0.430000e-4_rkx,0.410000e-4_rkx , &
   45.0_rkx,   1.81_rkx,274.00_rkx,0.181700e+1_rkx,0.190000e-4_rkx,0.130000e-4_rkx , &
   50.0_rkx,   0.99_rkx,277.00_rkx,0.986800e+0_rkx,0.630000e-5_rkx,0.430000e-5_rkx /
!
!  Standard POlar Winter ATMosphere
!
  data (((stdatm(ip,il,iz),ip=1,6),il=1,31),iz=5,5) &
  / 0.0_rkx,1013.00_rkx,257.10_rkx,0.137200e+4_rkx,0.120000e+1_rkx,0.410000e-4_rkx , &
    1.0_rkx, 887.80_rkx,259.10_rkx,0.119300e+4_rkx,0.120000e+1_rkx,0.410000e-4_rkx , &
    2.0_rkx, 777.50_rkx,255.90_rkx,0.105800e+4_rkx,0.940000e+0_rkx,0.410000e-4_rkx , &
    3.0_rkx, 679.80_rkx,252.70_rkx,0.936600e+3_rkx,0.680000e+0_rkx,0.430000e-4_rkx , &
    4.0_rkx, 593.20_rkx,247.70_rkx,0.833900e+3_rkx,0.410000e+0_rkx,0.450000e-4_rkx , &
    5.0_rkx, 515.80_rkx,240.90_rkx,0.745700e+3_rkx,0.200000e+0_rkx,0.470000e-4_rkx , &
    6.0_rkx, 446.70_rkx,234.10_rkx,0.664600e+3_rkx,0.980000e-1_rkx,0.490000e-4_rkx , &
    7.0_rkx, 385.30_rkx,227.30_rkx,0.590400e+3_rkx,0.540000e-1_rkx,0.710000e-4_rkx , &
    8.0_rkx, 330.80_rkx,220.60_rkx,0.522600e+3_rkx,0.110000e-1_rkx,0.900000e-4_rkx , &
    9.0_rkx, 282.90_rkx,217.20_rkx,0.453800e+3_rkx,0.840000e-2_rkx,0.160000e-3_rkx , &
   10.0_rkx, 241.80_rkx,217.20_rkx,0.387900e+3_rkx,0.550000e-2_rkx,0.240000e-3_rkx , &
   11.0_rkx, 206.70_rkx,217.20_rkx,0.331500e+3_rkx,0.380000e-2_rkx,0.320000e-3_rkx , &
   12.0_rkx, 176.60_rkx,217.20_rkx,0.283400e+3_rkx,0.260000e-2_rkx,0.430000e-3_rkx , &
   13.0_rkx, 151.00_rkx,217.20_rkx,0.242200e+3_rkx,0.180000e-2_rkx,0.470000e-3_rkx , &
   14.0_rkx, 129.10_rkx,217.20_rkx,0.207100e+3_rkx,0.100000e-2_rkx,0.490000e-3_rkx , &
   15.0_rkx, 110.30_rkx,217.20_rkx,0.177000e+3_rkx,0.760000e-3_rkx,0.560000e-3_rkx , &
   16.0_rkx,  94.31_rkx,216.60_rkx,0.151700e+3_rkx,0.640000e-3_rkx,0.620000e-3_rkx , &
   17.0_rkx,  80.58_rkx,216.00_rkx,0.130000e+3_rkx,0.560000e-3_rkx,0.620000e-3_rkx , &
   18.0_rkx,  68.82_rkx,215.40_rkx,0.111300e+3_rkx,0.500000e-3_rkx,0.620000e-3_rkx , &
   19.0_rkx,  58.75_rkx,214.80_rkx,0.952900e+2_rkx,0.490000e-3_rkx,0.600000e-3_rkx , &
   20.0_rkx,  50.14_rkx,214.10_rkx,0.815500e+2_rkx,0.450000e-3_rkx,0.560000e-3_rkx , &
   21.0_rkx,  42.77_rkx,213.60_rkx,0.697600e+2_rkx,0.510000e-3_rkx,0.510000e-3_rkx , &
   22.0_rkx,  36.47_rkx,213.00_rkx,0.596600e+2_rkx,0.510000e-3_rkx,0.470000e-3_rkx , &
   23.0_rkx,  31.09_rkx,212.40_rkx,0.510000e+2_rkx,0.540000e-3_rkx,0.430000e-3_rkx , &
   24.0_rkx,  26.49_rkx,211.80_rkx,0.435800e+2_rkx,0.600000e-3_rkx,0.360000e-3_rkx , &
   25.0_rkx,  22.56_rkx,211.20_rkx,0.372200e+2_rkx,0.670000e-3_rkx,0.320000e-3_rkx , &
   30.0_rkx,  10.20_rkx,216.00_rkx,0.164500e+2_rkx,0.360000e-3_rkx,0.150000e-3_rkx , &
   35.0_rkx,   4.70_rkx,222.20_rkx,0.736800e+1_rkx,0.110000e-3_rkx,0.920000e-4_rkx , &
   40.0_rkx,   2.24_rkx,234.70_rkx,0.333000e+1_rkx,0.430000e-4_rkx,0.410000e-4_rkx , &
   45.0_rkx,   1.11_rkx,247.00_rkx,0.156900e+1_rkx,0.190000e-4_rkx,0.130000e-4_rkx , &
   50.0_rkx,   0.57_rkx,259.00_rkx,0.768200e+0_rkx,0.630000e-5_rkx,0.430000e-5_rkx /
!
   contains

     real(rkx) function stdatm_val_noseason(lat,plev,ival)
       implicit none
       real(rkx) , intent(in) :: lat
       real(rkx) , intent(in) :: plev
       integer(ik4) , intent(in) :: ival
       integer(ik4) :: kp1 , kp2
       real(rkx) :: wtp1 , wtp2 , wtl2 , wtl1 , wts1 , wts2
       real(rkx) :: vs1 , vs2

       wts1 = 0.5_rkx
       wts2 = 1.0_rkx-wts1
       if ( abs(lat) >= 45.0_rkx ) then
         wtl1 = ((90.0_rkx-abs(lat))/45.0_rkx)**3
         wtl2 = 1.0_rkx-wtl1
         kp1 = find_klev(plev,ipolarwinter)
         kp2 = find_klev(plev,ipolarsummer)
         wtp1 = plev_wgt(kp1,plev,ipolarwinter)
         wtp2 = plev_wgt(kp1,plev,ipolarsummer)
         vs1 = stdatm(ival,kp1,ipolarwinter)*wtp1 + &
               stdatm(ival,kp1+1,ipolarwinter)*(d_one-wtp1)
         vs2 = stdatm(ival,kp1,ipolarsummer)*wtp1 + &
               stdatm(ival,kp1+1,ipolarsummer)*(d_one-wtp1)
         stdatm_val_noseason = vs1*wts1+vs2*wts2
         kp1 = find_klev(plev,imidlatwinter)
         kp2 = find_klev(plev,imidlatsummer)
         wtp1 = plev_wgt(kp1,plev,imidlatwinter)
         wtp2 = plev_wgt(kp1,plev,imidlatsummer)
         vs1 = stdatm(ival,kp1,imidlatwinter)*wtp1 + &
               stdatm(ival,kp1+1,imidlatwinter)*(d_one-wtp1)
         vs2 = stdatm(ival,kp1,imidlatsummer)*wtp1 + &
               stdatm(ival,kp1+1,imidlatsummer)*(d_one-wtp1)
         stdatm_val_noseason = wtl2*stdatm_val_noseason + &
                               wtl1*(vs1*wts1+vs2*wts2)
       else
         wtl1 = (abs(lat)/45.0_rkx)**2
         wtl2 = 1.0_rkx-wtl1
         kp1 = find_klev(plev,itropical)
         wtp1 = plev_wgt(kp1,plev,itropical)
         stdatm_val_noseason = stdatm(ival,kp1,itropical)*wtp1 + &
                      stdatm(ival,kp1+1,itropical)*(d_one-wtp1)
         kp1 = find_klev(plev,imidlatwinter)
         kp2 = find_klev(plev,imidlatsummer)
         wtp1 = plev_wgt(kp1,plev,imidlatwinter)
         wtp2 = plev_wgt(kp1,plev,imidlatsummer)
         vs1 = stdatm(ival,kp1,imidlatwinter)*wtp1 + &
               stdatm(ival,kp1+1,imidlatwinter)*(d_one-wtp1)
         vs2 = stdatm(ival,kp1,imidlatsummer)*wtp1 + &
               stdatm(ival,kp1+1,imidlatsummer)*(d_one-wtp1)
         stdatm_val_noseason = wtl2*stdatm_val_noseason + &
                               wtl1*(vs1*wts1+vs2*wts2)
       end if
     end function stdatm_val_noseason

     real(rkx) function stdatm_val_seasonal(jday,lat,plev,ival)
       implicit none
       real(rk8) , intent(in) :: jday
       real(rkx) , intent(in) :: lat
       real(rkx) , intent(in) :: plev
       integer(ik4) , intent(in) :: ival
       integer(ik4) :: kp1 , kp2
       real(rkx) :: wtl1 , wtl2 , wts1 , wts2 , wtp1 , wtp2
       real(rkx) :: vs1 , vs2

       stdatm_val_seasonal = dmissval

       if ( abs(lat) >= 45.0_rkx ) then
         wtl1 = ((90.0_rkx-abs(lat))/45.0_rkx)**3
         wtl2 = 1.0_rkx-wtl1
         wts1 = winter_wgt(jday)
         wts2 = 1.0_rkx-wts1
         kp1 = find_klev(plev,ipolarwinter)
         kp2 = find_klev(plev,ipolarsummer)
         wtp1 = plev_wgt(kp1,plev,ipolarwinter)
         wtp2 = plev_wgt(kp1,plev,ipolarsummer)
         vs1 = stdatm(ival,kp1,ipolarwinter)*wtp1 + &
               stdatm(ival,kp1+1,ipolarwinter)*(d_one-wtp1)
         vs2 = stdatm(ival,kp1,ipolarsummer)*wtp1 + &
               stdatm(ival,kp1+1,ipolarsummer)*(d_one-wtp1)
         stdatm_val_seasonal = vs1*wts1+vs2*wts2
         kp1 = find_klev(plev,imidlatwinter)
         kp2 = find_klev(plev,imidlatsummer)
         wtp1 = plev_wgt(kp1,plev,imidlatwinter)
         wtp2 = plev_wgt(kp1,plev,imidlatsummer)
         vs1 = stdatm(ival,kp1,imidlatwinter)*wtp1 + &
               stdatm(ival,kp1+1,imidlatwinter)*(d_one-wtp1)
         vs2 = stdatm(ival,kp1,imidlatsummer)*wtp1 + &
               stdatm(ival,kp1+1,imidlatsummer)*(d_one-wtp1)
         stdatm_val_seasonal = wtl2*stdatm_val_seasonal + &
                               wtl1*(vs1*wts1+vs2*wts2)
       else
         wtl1 = (abs(lat)/45.0_rkx)**2
         wtl2 = 1.0_rkx-wtl1
         wts1 = winter_wgt(jday)
         wts2 = 1.0_rkx-wts1
         kp1 = find_klev(plev,itropical)
         wtp1 = plev_wgt(kp1,plev,itropical)
         stdatm_val_seasonal = stdatm(ival,kp1,itropical)*wtp1 + &
                      stdatm(ival,kp1+1,itropical)*(d_one-wtp1)
         kp1 = find_klev(plev,imidlatwinter)
         kp2 = find_klev(plev,imidlatsummer)
         wtp1 = plev_wgt(kp1,plev,imidlatwinter)
         wtp2 = plev_wgt(kp1,plev,imidlatsummer)
         vs1 = stdatm(ival,kp1,imidlatwinter)*wtp1 + &
               stdatm(ival,kp1+1,imidlatwinter)*(d_one-wtp1)
         vs2 = stdatm(ival,kp1,imidlatsummer)*wtp1 + &
               stdatm(ival,kp1+1,imidlatsummer)*(d_one-wtp1)
         stdatm_val_seasonal = wtl2*stdatm_val_seasonal + &
                               wtl1*(vs1*wts1+vs2*wts2)
       end if
     end function stdatm_val_seasonal

     real(rkx) function winter_wgt(jday)
       implicit none
       real(rkx) , intent(in) :: jday
       real(rk8) :: dis
       dis = ((half_dayspy-jday-sixteenth_dayspy+d_one)/dayspy)*mathpi
       winter_wgt = real(sin(dis)**2,rkx)
     end function winter_wgt

     integer(ik4) function find_klev(plev,izone)
       implicit none
       integer(ik4) , intent(in) :: izone
       real(rkx) , intent(in) :: plev
       integer(ik4) :: k
       find_klev = 1
       do k = 2 , n_atmlevls
         find_klev = k-1
         if ( plev > stdatm(istdatm_prsmb,k,izone) ) exit
       end do
     end function find_klev

     integer(ik4) function find_zlev(z,izone)
       implicit none
       integer(ik4) , intent(in) :: izone
       real(rkx) , intent(in) :: z
       integer(ik4) :: k
       find_zlev = 1
       do k = 2 , n_atmlevls
         find_zlev = k-1
         if ( z < d_1000 * stdatm(istdatm_hgtkm,k,izone) ) exit
       end do
     end function find_zlev

     real(rkx) function plev_wgt(k,plev,izone)
       implicit none
       integer(ik4) , intent(in) :: k , izone
       real(rkx) , intent(in) :: plev
       if ( plev >= stdatm(istdatm_prsmb,k,izone) ) then
         plev_wgt = d_one
       else if ( plev <= stdatm(istdatm_prsmb,k+1,izone) ) then
         plev_wgt = d_zero
       else
         plev_wgt = log(plev/stdatm(istdatm_prsmb,k+1,izone)) / &
             log(stdatm(istdatm_prsmb,k,izone)/stdatm(istdatm_prsmb,k+1,izone))
       end if
     end function plev_wgt

     real(rkx) function stdlrate(jday,lat)
       implicit none
       real(rkx) , intent(in) :: jday
       real(rkx) , intent(in) :: lat
       real(rkx) :: wts1 , wts2 , wtp1 , wtp2 , vs1 , vs2 , vz1 , vz2
       real(rkx) :: wtl1 , wtl2 , vs3 , vs4
       integer(ik4) :: kp1 , kp2
       real(rkx) :: polar , midlat , tropical
       if ( abs(lat) >= 45.0_rkx ) then
         wtl1 = ((90.0_rkx-abs(lat))/45.0_rkx)**3
         wtl2 = 1.0_rkx-wtl1
         wts1 = winter_wgt(jday)
         wts2 = 1.0_rkx-wts1
         kp1 = find_klev(300.0_rkx,ipolarwinter)
         kp2 = find_klev(300.0_rkx,ipolarsummer)
         wtp1 = plev_wgt(kp1,300.0_rkx,ipolarwinter)
         wtp2 = plev_wgt(kp1,300.0_rkx,ipolarsummer)
         vs1 = stdatm(istdatm_tempk,kp1,ipolarwinter)*wtp1 + &
               stdatm(istdatm_tempk,kp1+1,ipolarwinter)*(d_one-wtp1)
         vs2 = stdatm(istdatm_tempk,kp1,ipolarsummer)*wtp1 + &
               stdatm(istdatm_tempk,kp1+1,ipolarsummer)*(d_one-wtp1)
         vs3 = stdatm(istdatm_tempk,1,ipolarwinter)
         vs4 = stdatm(istdatm_tempk,1,ipolarsummer)
         vz1 = stdatm(istdatm_hgtkm,kp1,ipolarwinter)*wtp1 + &
               stdatm(istdatm_hgtkm,kp1+1,ipolarwinter)*(d_one-wtp1)
         vz2 = stdatm(istdatm_hgtkm,kp1,ipolarsummer)*wtp1 + &
               stdatm(istdatm_hgtkm,kp1+1,ipolarsummer)*(d_one-wtp1)
         polar = ((vs1*wts1+vs2*wts2) - (vs3*wts1+vs4*wts2)) / &
                 (d_1000*(vz1*wts1+vz2*wts2))
         kp1 = find_klev(300.0_rkx,imidlatwinter)
         kp2 = find_klev(300.0_rkx,imidlatsummer)
         wtp1 = plev_wgt(kp1,300.0_rkx,imidlatwinter)
         wtp2 = plev_wgt(kp1,300.0_rkx,imidlatsummer)
         vs1 = stdatm(istdatm_tempk,kp1,imidlatwinter)*wtp1 + &
               stdatm(istdatm_tempk,kp1+1,imidlatwinter)*(d_one-wtp1)
         vs2 = stdatm(istdatm_tempk,kp1,imidlatsummer)*wtp1 + &
               stdatm(istdatm_tempk,kp1+1,imidlatsummer)*(d_one-wtp1)
         vs3 = stdatm(istdatm_tempk,1,imidlatwinter)
         vs4 = stdatm(istdatm_tempk,1,imidlatsummer)
         vz1 = stdatm(istdatm_hgtkm,kp1,imidlatwinter)*wtp1 + &
               stdatm(istdatm_hgtkm,kp1+1,imidlatwinter)*(d_one-wtp1)
         vz2 = stdatm(istdatm_hgtkm,kp1,imidlatsummer)*wtp1 + &
               stdatm(istdatm_hgtkm,kp1+1,imidlatsummer)*(d_one-wtp1)
         midlat = ((vs1*wts1+vs2*wts2) - (vs3*wts1+vs4*wts2)) / &
                  (d_1000*(vz1*wts1+vz2*wts2))
         stdlrate = wtl2 * polar + wtl1 * midlat
       else
         wtl1 = (abs(lat)/45.0_rkx)**2
         wtl2 = 1.0_rkx-wtl1
         kp1 = find_klev(300.0_rkx,itropical)
         wtp1 = plev_wgt(kp1,300.0_rkx,itropical)
         vs1 = stdatm(istdatm_tempk,kp1,itropical)*wtp1 + &
               stdatm(istdatm_tempk,kp1+1,itropical)*(d_one-wtp1)
         vs2 = stdatm(istdatm_tempk,1,itropical)
         vz1 = stdatm(istdatm_hgtkm,kp1,itropical)*wtp1 + &
               stdatm(istdatm_hgtkm,kp1+1,itropical)*(d_one-wtp1)
         tropical = (vs1 - vs2)/(d_1000*vz1)
         wts1 = winter_wgt(jday)
         wts2 = 1.0_rkx-wts1
         kp1 = find_klev(300.0_rkx,imidlatwinter)
         kp2 = find_klev(300.0_rkx,imidlatsummer)
         wtp1 = plev_wgt(kp1,300.0_rkx,imidlatwinter)
         wtp2 = plev_wgt(kp1,300.0_rkx,imidlatsummer)
         vs1 = stdatm(istdatm_tempk,kp1,imidlatwinter)*wtp1 + &
               stdatm(istdatm_tempk,kp1+1,imidlatwinter)*(d_one-wtp1)
         vs2 = stdatm(istdatm_tempk,kp1,imidlatsummer)*wtp1 + &
               stdatm(istdatm_tempk,kp1+1,imidlatsummer)*(d_one-wtp1)
         vs3 = stdatm(istdatm_tempk,1,imidlatwinter)
         vs4 = stdatm(istdatm_tempk,1,imidlatsummer)
         vz1 = stdatm(istdatm_hgtkm,kp1,imidlatwinter)*wtp1 + &
               stdatm(istdatm_hgtkm,kp1+1,imidlatwinter)*(d_one-wtp1)
         vz2 = stdatm(istdatm_hgtkm,kp1,imidlatsummer)*wtp1 + &
               stdatm(istdatm_hgtkm,kp1+1,imidlatsummer)*(d_one-wtp1)
         midlat = ((vs1*wts1+vs2*wts2) - (vs3*wts1+vs4*wts2)) / &
                  (d_1000*(vz1*wts1+vz2*wts2))
         stdlrate = wtl2 * midlat + wtl1 * tropical
       end if
    end function stdlrate

end module mod_stdatm
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

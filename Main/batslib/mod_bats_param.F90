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
!
module mod_bats_param
!
  use mod_intkinds
  use mod_realkinds
  !
  ! Landuse Categories (BATS 1e , Table 1 , pag 17 + added classes)
  !
  !      1 = Crop/mixed farming
  !      2 = Short grass
  !      3 = Evergreen needleleaf tree
  !      4 = Deciduous needleleaf tree
  !      5 = Deciduous broadleaf tree
  !      6 = Evergreen broadleaf tree
  !      7 = Tall grass
  !      8 = Desert
  !      9 = Tundra
  !     10 = Irrigated Crop
  !     11 = Semi-desert
  !     12 = Ice cap/glacier
  !     13 = Bog or marsh
  !     14 = Internal Water
  !     15 = Ocean
  !     16 = Evergreen shrub
  !     17 = Deciduous shrub
  !     18 = Mixed Woodland
  !     19 = Forest/Field mosaic      ( added class )
  !     20 = Water and Land mixture   ( added class )
  !     21 = Urban                    ( added class )
  !     22 = Sub-Urban                ( added class )
  !
  implicit none

  public

  real(rkx) , dimension(8) :: solour
  real(rkx) , dimension(22) :: albvgl , albvgs , crough , deprv ,     &
                             deptv , depuv , displa , fc , freza ,  &
                             frezu , rough , rsmin , sai , seasf ,  &
                             sqrtdi , mfcv , xla , xlai0 , rootf ,  &
                             slmo , lndemiss , seasemi
  real(rkx) , dimension(12) :: bee , skrat , xmofc , xmohyd , xmopor ,&
                             xmosuc , xmowil
  integer(ik4) , dimension(22) :: iexsol , kolsol
  ! 0.005 ccm specific
  ! 0.01  high
  ! 0.04  Antartic
  ! 0.02  Artic best
  real(rkx) , parameter :: aarea = 0.02_rkx
  real(rkx) , parameter :: minsigf = 0.001_rkx
  ! Seasonal crop cutoff
  logical , parameter :: lcrop_cutoff = .false.
  !
  ! slmo is initial surface moisture availability in fraction of one
  !
  data slmo / 0.50_rkx , 0.50_rkx , 0.50_rkx , 0.50_rkx , 0.50_rkx , &
              0.50_rkx , 0.50_rkx , 0.50_rkx , 0.50_rkx , 0.50_rkx , &
              0.50_rkx , 0.50_rkx , 0.50_rkx , 1.00_rkx , 1.00_rkx , &
              0.50_rkx , 0.50_rkx , 0.50_rkx , 0.50_rkx , 0.50_rkx , &
              0.50_rkx , 0.50_rkx /
  !
  ! mfcv is maximum fractional vegetation cover
  ! BATS 1e , Table 2, a, pag 21
  !
  data mfcv  / 0.85_rkx , 0.80_rkx , 0.80_rkx , 0.80_rkx , 0.80_rkx , 0.90_rkx , &
               0.80_rkx , 0.00_rkx , 0.60_rkx , 0.80_rkx , 0.10_rkx , 0.00_rkx , &
               0.80_rkx , 0.00_rkx , 0.00_rkx , 0.80_rkx , 0.80_rkx , 0.80_rkx , &
               0.80_rkx , 0.80_rkx , 0.05_rkx , 0.40_rkx /
  !
  ! seasf is the difference between mfcv and cover at temperature of 269k
  ! BATS 1e , Table 2, b, pag 21
  !
  data seasf / 0.60_rkx , 0.10_rkx , 0.10_rkx , 0.30_rkx , 0.30_rkx , 0.50_rkx , &
               0.30_rkx , 0.00_rkx , 0.20_rkx , 0.60_rkx , 0.10_rkx , 0.00_rkx , &
               0.40_rkx , 0.00_rkx , 0.00_rkx , 0.20_rkx , 0.30_rkx , 0.20_rkx , &
               0.40_rkx , 0.40_rkx , 0.05_rkx , 0.15_rkx /
  !
  ! rough is an aerodynamic roughness length (m) =approx 0.1*veg
  ! height also used snow masking depth in subrout albedo
  ! BATS 1e , Table 2, c, pag 21
  !
  !data rough / 0.0600_rkx , 0.0200_rkx , 1.0000_rkx , 1.0000_rkx , 0.8000_rkx , &
  !             2.0000_rkx , 0.1000_rkx , 0.0500_rkx , 0.0400_rkx , 0.0600_rkx , &
  !             0.1000_rkx , 0.0100_rkx , 0.0300_rkx , 0.0024_rkx , 0.0024_rkx , &
  !             0.1000_rkx , 0.1000_rkx , 0.8000_rkx , 0.3000_rkx , 0.3000_rkx , &
  !             1.5000_rkx , 1.0000_rkx /
  ! Following Stulls "Meteorology for Scientists and Engineers",
  ! Davenport-Wieringa roughness-length classifcations.
  data rough / 0.1000_rkx , 0.0300_rkx , 1.0000_rkx , 1.0000_rkx , 1.0000_rkx , &
               1.0000_rkx , 0.3000_rkx , 0.0050_rkx , 0.0300_rkx , 0.1000_rkx , &
               0.0300_rkx , 0.0050_rkx , 0.1000_rkx , 0.0002_rkx , 0.0004_rkx , &
               0.2500_rkx , 0.1000_rkx , 1.0000_rkx , 0.5000_rkx , 0.3000_rkx , &
               2.0000_rkx , 1.0000_rkx /
  !
  ! displacement height (meter)
  ! if great parts of veg. are covered by snow, use displa=0
  ! because then the new displa-theory is not valid
  !
  data displa / 0.00_rkx , 0.00_rkx , 9.00_rkx , 9.00_rkx , 0.00_rkx ,18.00_rkx , &
                0.00_rkx , 0.00_rkx , 0.00_rkx , 0.00_rkx , 0.00_rkx , 0.00_rkx , &
                0.00_rkx , 0.00_rkx , 0.00_rkx , 0.00_rkx , 0.00_rkx , 0.00_rkx , &
                0.00_rkx , 0.00_rkx , 6.00_rkx , 2.50_rkx /
  !
  ! minimum stomatal resistance (s/m)
  ! Model is expecially sensible to this, tune it carefully
  ! increase in rsmin will lead to a decrease in evapotranspration
  ! BATS 1e , Table 2, i, pag 21
  !
  !data rsmin / 120.0_rkx , 200.0_rkx , 200.0_rkx , 200.0_rkx , 200.0_rkx , 150.0_rkx ,&
  !             200.0_rkx , 200.0_rkx , 200.0_rkx , 200.0_rkx , 200.0_rkx , 200.0_rkx ,&
  !             200.0_rkx , 200.0_rkx , 200.0_rkx , 200.0_rkx , 200.0_rkx , 200.0_rkx ,&
  !             200.0_rkx , 200.0_rkx , 120.0_rkx , 120.0_rkx /
  ! Modified by Laura Mariotti
  data rsmin /  45.0_rkx , 200.0_rkx , 200.0_rkx , 200.0_rkx , 200.0_rkx ,  50.0_rkx , &
               200.0_rkx , 200.0_rkx , 200.0_rkx , 200.0_rkx , 200.0_rkx , 200.0_rkx , &
               200.0_rkx , 200.0_rkx , 200.0_rkx , 200.0_rkx , 200.0_rkx , 200.0_rkx , &
               200.0_rkx , 200.0_rkx , 120.0_rkx ,  60.0_rkx /
  !
  ! maximum leaf area index (ratio unit cover per unit ground)
  ! BATS 1e , Table 2, j, pag 21
  !
  !data xla / 6.0_rkx , 2.0_rkx , 6.0_rkx , 6.0_rkx , 6.0_rkx , 6.0_rkx , &
  !           6.0_rkx , 0.0_rkx , 6.0_rkx , 6.0_rkx , 6.0_rkx , 0.0_rkx , &
  !           6.0_rkx , 0.0_rkx , 0.0_rkx , 6.0_rkx , 6.0_rkx , 6.0_rkx , &
  !           6.0_rkx , 6.0_rkx , 1.0_rkx , 2.0_rkx /
  data xla / 6.0_rkx , 2.0_rkx , 6.0_rkx , 6.0_rkx , 6.0_rkx , 6.0_rkx , &
             3.0_rkx , 0.0_rkx , 2.0_rkx , 4.0_rkx , 1.0_rkx , 0.0_rkx , &
             4.0_rkx , 0.0_rkx , 0.0_rkx , 4.0_rkx , 4.0_rkx , 5.0_rkx , &
             4.0_rkx , 1.0_rkx , 1.0_rkx , 2.0_rkx /
  !
  ! minimum leaf area index **lai depends on temp as veg cover
  ! BATS 1e , Table 2, k, pag 21
  !
  !data xlai0 / 0.5_rkx , 0.5_rkx , 5.0_rkx , 1.0_rkx , 1.0_rkx , 5.0_rkx , &
  !             0.5_rkx , 0.0_rkx , 0.5_rkx , 0.5_rkx , 0.5_rkx , 0.0_rkx , &
  !             0.5_rkx , 0.0_rkx , 0.0_rkx , 5.0_rkx , 1.0_rkx , 3.0_rkx , &
  !             0.5_rkx , 0.5_rkx , 0.5_rkx , 0.5_rkx /
  data xlai0 / 0.5_rkx , 0.5_rkx , 5.0_rkx , 1.0_rkx , 1.0_rkx , 5.0_rkx , &
               1.0_rkx , 0.0_rkx , 0.5_rkx , 2.0_rkx , 0.5_rkx , 0.0_rkx , &
               2.0_rkx , 0.0_rkx , 0.0_rkx , 3.0_rkx , 1.0_rkx , 3.0_rkx , &
               0.5_rkx , 1.0_rkx , 0.5_rkx , 1.0_rkx /
  !
  ! stem and dead matter area index (projected area of non-transpiring sfcs)
  ! BATS 1e , Table 2, l, pag 21
  !
  data sai / 0.5_rkx , 4.0_rkx , 2.0_rkx , 2.0_rkx , 2.0_rkx , 2.0_rkx , &
             2.0_rkx , 0.5_rkx , 0.5_rkx , 2.0_rkx , 2.0_rkx , 2.0_rkx , &
             2.0_rkx , 2.0_rkx , 2.0_rkx , 2.0_rkx , 2.0_rkx , 2.0_rkx , &
             2.0_rkx , 2.0_rkx , 0.5_rkx , 0.5_rkx /
  !
  ! inverse square root of leaf dimension - used for calculating fluxes
  ! from foliage (m**0.5)
  ! BATS 1e , Table 2, m, pag 21
  !
  data sqrtdi / 10.0_rkx , 5.0_rkx , 5.0_rkx , 5.0_rkx , 5.0_rkx , 5.0_rkx , &
                 5.0_rkx , 5.0_rkx , 5.0_rkx , 5.0_rkx , 5.0_rkx , 5.0_rkx , &
                 5.0_rkx , 5.0_rkx , 5.0_rkx , 5.0_rkx , 5.0_rkx , 5.0_rkx , &
                 5.0_rkx , 5.0_rkx , 5.0_rkx , 5.0_rkx /
  !
  ! fc = light dependence of stomatal resistance (m^2/W)
  ! BATS 1e , Table 2, n, pag 21
  !
  data fc / 0.02_rkx , 0.02_rkx , 0.06_rkx , 0.06_rkx , 0.06_rkx , 0.06_rkx , &
            0.02_rkx , 0.02_rkx , 0.02_rkx , 0.02_rkx , 0.02_rkx , 0.02_rkx , &
            0.02_rkx , 0.02_rkx , 0.02_rkx , 0.02_rkx , 0.02_rkx , 0.06_rkx , &
            0.02_rkx , 0.02_rkx , 0.02_rkx , 0.02_rkx /
  !
  ! depuv is depth of upper soil layer (mm)
  ! BATS 1e , Table 2, e, pag 21
  !
  data depuv / 100.0_rkx , 100.0_rkx , 100.0_rkx , 100.0_rkx , 100.0_rkx , &
               100.0_rkx , 100.0_rkx , 100.0_rkx , 100.0_rkx , 100.0_rkx , &
               100.0_rkx , 100.0_rkx , 100.0_rkx , 100.0_rkx , 100.0_rkx , &
               100.0_rkx , 100.0_rkx , 100.0_rkx , 100.0_rkx , 100.0_rkx , &
               100.0_rkx , 100.0_rkx /
  !
  ! deprv is depth of root zone (mm)
  ! BATS 1e , Table 2, d, pag 21
  !
  !data deprv / 1000.0_rkx , 1000.0_rkx , 1500.0_rkx , 1500.0_rkx , 2000.0_rkx , &
  !             1500.0_rkx , 1000.0_rkx , 1000.0_rkx , 1000.0_rkx , 1000.0_rkx , &
  !             1000.0_rkx , 1000.0_rkx , 1000.0_rkx , 1000.0_rkx , 1000.0_rkx , &
  !             1000.0_rkx , 1000.0_rkx , 2000.0_rkx , 1000.0_rkx , 1000.0_rkx , &
  !             1000.0_rkx , 1000.0_rkx /
  data deprv / 1000.0_rkx , 1000.0_rkx , 1500.0_rkx , 1500.0_rkx , 2000.0_rkx , &
               1500.0_rkx , 1000.0_rkx , 1000.0_rkx , 1000.0_rkx , 1000.0_rkx , &
               1000.0_rkx , 1000.0_rkx , 1000.0_rkx , 1000.0_rkx , 1000.0_rkx , &
               1000.0_rkx , 1000.0_rkx , 2000.0_rkx , 2000.0_rkx , 2000.0_rkx , &
               1000.0_rkx , 1000.0_rkx /
  !
  ! deptv is depth of total soil (mm)
  !
  data deptv / 3000.0_rkx , 3000.0_rkx , 3000.0_rkx , 3000.0_rkx , 3000.0_rkx , &
               3000.0_rkx , 3000.0_rkx , 3000.0_rkx , 3000.0_rkx , 3000.0_rkx , &
               3000.0_rkx , 3000.0_rkx , 3000.0_rkx , 3000.0_rkx , 3000.0_rkx , &
               3000.0_rkx , 3000.0_rkx , 3000.0_rkx , 3000.0_rkx , 3000.0_rkx , &
               3000.0_rkx , 3000.0_rkx /
  !
  ! kolsol is soil color type (see subr. albedo)
  !
  ! 1 -> Light colour
  ! 2 -> Dark colour
  !
  !data kolsol / 6 , 4 , 5 , 5 , 5 , 5 , 5 , 3 , 4 , 4 , 2 , 1 , 6 ,   &
  !              6 , 6 , 5 , 4 , 5 , 5 , 5 , 4 , 4 /
  data kolsol / 6 , 4 , 5 , 5 , 5 , 5 , 5 , 1 , 4 , 4 , 2 , 1 , 6 ,   &
                6 , 6 , 5 , 4 , 5 , 5 , 5 , 4 , 4 /
  !
  ! Saturated soil albedo for different soil coloures for short wave rad.
  ! BATS 1e , Table 3, II , b, first line, pag 27
  !
  !data solour / 0.12_rkx , 0.11_rkx , 0.10_rkx , 0.09_rkx , 0.08_rkx , &
  !              0.07_rkx , 0.06_rkx , 0.05_rkx /
  data solour / 0.16_rkx , 0.15_rkx , 0.10_rkx , 0.09_rkx , 0.08_rkx , &
                0.07_rkx , 0.06_rkx , 0.05_rkx /
  !
  ! iexsol is soil texture type (see subr soilbc)
  !
  ! 1  -> Sand
  ! 12 -> Clay
  !
  !data iexsol / 6 , 6 , 6 , 6 , 7 , 8 , 6 , 3 , 6 , 6 , 5 , 12 , 6 ,  &
  !              6 , 6 , 6 , 5 , 6 , 6 , 6 , 12 , 8 /
  data iexsol / 6 , 4 , 6 , 6 , 7 , 8 , 5 , 1 , 6 , 6 , 2 , 12 , 9 ,  &
                0 , 0 , 6 , 5 , 6 , 6 , 8 , 12 , 8 /
  !
  ! xmopor is fraction of soil that is voids
  ! BATS 1e , Table 3, I , a, pag 27
  !
  !data xmopor / 0.33_rkx , 0.36_rkx , 0.39_rkx , 0.42_rkx , 0.45_rkx , 0.48_rkx , &
  !              0.51_rkx , 0.54_rkx , 0.57_rkx , 0.60_rkx , 0.63_rkx , 0.66_rkx /
  data xmopor / 0.13_rkx , 0.26_rkx , 0.39_rkx , 0.42_rkx , 0.45_rkx , 0.48_rkx , &
                0.51_rkx , 0.54_rkx , 0.57_rkx , 0.60_rkx , 0.63_rkx , 0.66_rkx /
  !
  ! xmosuc is the minimum soil suction (mm)
  ! BATS 1e , Table 3, I , b, pag 27
  !
  !data xmosuc / 30.0_rkx ,  30.0_rkx ,  30.0_rkx , 200.0_rkx , 200.0_rkx , 200.0_rkx ,&
  !             200.0_rkx , 200.0_rkx , 200.0_rkx , 200.0_rkx , 200.0_rkx , 200.0_rkx /
  data xmosuc / 10.0_rkx ,  20.0_rkx ,  30.0_rkx , 200.0_rkx , 200.0_rkx , 200.0_rkx , &
               200.0_rkx , 200.0_rkx , 200.0_rkx , 200.0_rkx , 200.0_rkx , 200.0_rkx /
  !
  ! xmohyd is the saturated hydraulic conductivity (mm/s)
  ! BATS 1e , Table 3, I , c, pag 27
  !
  data xmohyd / 0.0200_rkx , 0.0800_rkx , 0.0320_rkx , 0.0130_rkx , 0.0089_rkx , &
                0.0063_rkx , 0.0045_rkx , 0.0032_rkx , 0.0022_rkx , 0.0016_rkx , &
                0.0011_rkx , 0.0008_rkx /
  !
  ! xmowilt is fraction of water content at which permanent wilting occurs
  ! (transpiration ceases)
  ! BATS 1e , Table 3, I , f, pag 27
  !
  !data xmowil / 0.088_rkx , 0.119_rkx , 0.151_rkx , 0.266_rkx , 0.300_rkx , &
  !              0.332_rkx , 0.378_rkx , 0.419_rkx , 0.455_rkx , 0.487_rkx , &
  !              0.516_rkx , 0.542_rkx /
  data xmowil / 0.095_rkx , 0.128_rkx , 0.161_rkx , 0.266_rkx , 0.300_rkx , &
                0.332_rkx , 0.378_rkx , 0.419_rkx , 0.455_rkx , 0.487_rkx , &
                0.516_rkx , 0.542_rkx /
  !
  ! Field capacity
  !
  data xmofc / 0.404_rkx , 0.477_rkx , 0.547_rkx , 0.614_rkx , 0.653_rkx ,  &
               0.688_rkx , 0.728_rkx , 0.763_rkx , 0.794_rkx , 0.820_rkx ,  &
               0.845_rkx , 0.866_rkx /
  !
  ! bee is the clapp and hornbereger "b" parameter
  ! BATS 1e , Table 3, I , e, pag 27
  !
  data bee / 3.5_rkx , 4.0_rkx , 4.5_rkx , 5.0_rkx , 5.5_rkx , 6.0_rkx , 6.8_rkx ,  &
             7.6_rkx , 8.4_rkx , 9.2_rkx ,10.0_rkx ,10.8_rkx /
  !
  ! bskrat is ratio of soil thermal conduc. to that of loam
  ! a function of texture
  ! BATS 1e , Table 3, I , d, pag 27
  !
  data skrat / 1.70_rkx , 1.50_rkx , 1.30_rkx , 1.20_rkx , 1.10_rkx , 1.00_rkx , &
               0.95_rkx , 0.90_rkx , 0.85_rkx , 0.80_rkx , 0.75_rkx , 0.70_rkx /
  !
  ! albvgs is vegetation albedo for wavelengths < 0.7 microns data
  ! BATS 1e , Table 2, I , g, pag 21
  !
  !data albvgs / 0.10_rkx , 0.10_rkx , 0.05_rkx , 0.05_rkx , 0.08_rkx , 0.04_rkx , &
  !              0.08_rkx , 0.20_rkx , 0.10_rkx , 0.08_rkx , 0.17_rkx , 0.80_rkx , &
  !              0.06_rkx , 0.07_rkx , 0.07_rkx , 0.05_rkx , 0.08_rkx , 0.06_rkx , &
  !              0.06_rkx , 0.06_rkx , 0.02_rkx , 0.06_rkx /
  data albvgs / 0.10_rkx , 0.10_rkx , 0.04_rkx , 0.04_rkx , 0.06_rkx , 0.04_rkx , &
                0.08_rkx , 0.20_rkx , 0.10_rkx , 0.08_rkx , 0.17_rkx , 0.80_rkx , &
                0.06_rkx , 0.07_rkx , 0.07_rkx , 0.05_rkx , 0.08_rkx , 0.05_rkx , &
                0.06_rkx , 0.06_rkx , 0.02_rkx , 0.06_rkx /
  !
  ! albvgl is vegetation albedo for wavelengths > 0.7 microns data
  ! BATS 1e , Table 2, I , h, pag 21
  !
  !data albvgl / 0.30_rkx , 0.30_rkx , 0.23_rkx , 0.23_rkx , 0.28_rkx , 0.20_rkx , &
  !              0.30_rkx , 0.40_rkx , 0.30_rkx , 0.28_rkx , 0.34_rkx , 0.60_rkx , &
  !              0.18_rkx , 0.20_rkx , 0.20_rkx , 0.23_rkx , 0.28_rkx , 0.24_rkx , &
  !              0.18_rkx , 0.18_rkx , 0.15_rkx , 0.18_rkx /
  data albvgl / 0.30_rkx , 0.30_rkx , 0.20_rkx , 0.20_rkx , 0.26_rkx , 0.20_rkx , &
                0.30_rkx , 0.40_rkx , 0.30_rkx , 0.28_rkx , 0.34_rkx , 0.60_rkx , &
                0.18_rkx , 0.20_rkx , 0.20_rkx , 0.23_rkx , 0.28_rkx , 0.23_rkx , &
                0.18_rkx , 0.18_rkx , 0.15_rkx , 0.18_rkx /
  !
  ! Fraction of water extracted by upper layer roots (saturated)
  ! BATS 1e , Table 2, I , h, pag 21
  !
  data rootf / 0.30_rkx , 0.80_rkx , 0.67_rkx , 0.67_rkx , 0.50_rkx , 0.80_rkx , &
               0.80_rkx , 0.90_rkx , 0.90_rkx , 0.30_rkx , 0.80_rkx , 0.50_rkx , &
               0.50_rkx , 0.50_rkx , 0.50_rkx , 0.50_rkx , 0.50_rkx , 0.50_rkx , &
               0.50_rkx , 0.50_rkx , 0.90_rkx , 0.50_rkx /
  !
  ! Emissivity coefficients (used in iemiss = 1)
  !
  data lndemiss / 0.983_rkx, 0.983_rkx, 0.983_rkx, 0.987_rkx, 0.981_rkx, 0.981_rkx, &
                  0.983_rkx, 0.965_rkx, 0.987_rkx, 0.985_rkx, 0.970_rkx, 0.993_rkx, &
                  0.992_rkx, 0.992_rkx, 0.992_rkx, 0.983_rkx, 0.972_rkx, 0.983_rkx, &
                  0.981_rkx, 0.991_rkx, 0.970_rkx, 0.972_rkx /
  !
  ! Seasonal variations (used in iemiss = 1)
  !
  data seasemi / 0.005_rkx, 0.002_rkx, 0.000_rkx, 0.004_rkx, 0.004_rkx, 0.000_rkx, &
                 0.002_rkx, 0.000_rkx, 0.000_rkx, 0.002_rkx, 0.000_rkx, 0.000_rkx, &
                 0.000_rkx, 0.000_rkx, 0.000_rkx, 0.000_rkx, 0.004_rkx, 0.002_rkx, &
                 0.004_rkx, 0.000_rkx, 0.000_rkx, 0.001_rkx /
  !
  !      1 = Crop/mixed farming
  !      2 = Short grass
  !      3 = Evergreen needleleaf tree
  !      4 = Deciduous needleleaf tree
  !      5 = Deciduous broadleaf tree
  !      6 = Evergreen broadleaf tree
  !      7 = Tall grass
  !      8 = Desert
  !      9 = Tundra
  !     10 = Irrigated Crop
  !     11 = Semi-desert
  !     12 = Ice cap/glacier
  !     13 = Bog or marsh
  !     14 = Internal Water
  !     15 = Ocean
  !     16 = Evergreen shrub
  !     17 = Deciduous shrub
  !     18 = Mixed Woodland
  !     19 = Forest/Field mosaic
  !     20 = Water and Land mixture
  !     21 = Urban
  !     22 = Sub-Urban
  !

end module mod_bats_param
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

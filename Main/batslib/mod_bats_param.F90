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

  real(rk8) , dimension(8) :: solour
  real(rk8) , dimension(22) :: albvgl , albvgs , crough , deprv ,     &
                             deptv , depuv , displa , fc , freza ,  &
                             frezu , rough , rsmin , sai , seasf ,  &
                             sqrtdi , mfcv , xla , xlai0 , rootf ,  &
                             slmo , lndemiss , seasemi
  real(rk8) , dimension(12) :: bee , skrat , xmofc , xmohyd , xmopor ,&
                             xmosuc , xmowil
  integer(ik4) , dimension(22) :: iexsol , kolsol
  ! 0.005 ccm specific
  ! 0.01  high
  ! 0.04  Antartic
  ! 0.02  Artic best
  real(rk8) , parameter :: aarea = 0.02D0
  real(rk8) , parameter :: minsigf = 0.001D+00
  !
  ! Saturated soil albedo for different soil coloures
  ! BATS 1e , Table 3, II , pag 27
  !
  !data solour / 0.12D0 , 0.11D0 , 0.10D0 , 0.09D0 , 0.08D0 ,        &
  !              0.07D0 , 0.06D0 , 0.05D0 /
  data solour / 0.16D0 , 0.15D0 , 0.10D0 , 0.09D0 , 0.08D0 ,        &
                0.07D0 , 0.06D0 , 0.05D0 /
  !
  ! mfcv is maximum fractional vegetation cover
  ! BATS 1e , Table 2, a, pag 21
  !
  !data mfcv  / 0.85D0 , 0.80D0 , 0.80D0 , 0.80D0 , 0.80D0 , 0.90D0 , &
  !             0.80D0 , 0.00D0 , 0.60D0 , 0.80D0 , 0.10D0 , 0.00D0 , &
  !             0.80D0 , 0.00D0 , 0.00D0 , 0.80D0 , 0.80D0 , 0.80D0 , &
  !             0.80D0 , 0.80D0 , 0.05D0 , 0.40D0 /
  data mfcv  / 0.85D0 , 0.80D0 , 0.80D0 , 0.80D0 , 0.80D0 , 0.90D0 , &
               0.80D0 , 0.00D0 , 0.60D0 , 0.80D0 , 0.35D0 , 0.00D0 , &
               0.80D0 , 0.00D0 , 0.00D0 , 0.80D0 , 0.80D0 , 0.80D0 , &
               0.80D0 , 0.80D0 , 0.05D0 , 0.40D0 /
  !
  ! seasf is the difference between mfcv and cover at temperature of 269k
  ! BATS 1e , Table 2, b, pag 21
  !
  data seasf / 0.60D0 , 0.10D0 , 0.10D0 , 0.30D0 , 0.30D0 , 0.50D0 , &
               0.30D0 , 0.00D0 , 0.20D0 , 0.60D0 , 0.10D0 , 0.00D0 , &
               0.40D0 , 0.00D0 , 0.00D0 , 0.20D0 , 0.30D0 , 0.20D0 , &
               0.40D0 , 0.40D0 , 0.05D0 , 0.15D0 /
  !
  ! rough is an aerodynamic roughness length (m) =approx 0.1*veg
  ! height also used snow masking depth in subrout albedo
  ! BATS 1e , Table 2, c, pag 21
  !
  !data rough / 0.0600D0 , 0.0200D0 , 1.0000D0 , 1.0000D0 , 0.8000D0 , &
  !             2.0000D0 , 0.1000D0 , 0.0500D0 , 0.0400D0 , 0.0600D0 , &
  !             0.1000D0 , 0.0100D0 , 0.0300D0 , 0.0024D0 , 0.0024D0 , &
  !             0.1000D0 , 0.1000D0 , 0.8000D0 , 0.3000D0 , 0.3000D0 , &
  !             1.5000D0 , 0.4000D0 /
  data rough / 0.0800D0 , 0.0500D0 , 1.0000D0 , 1.0000D0 , 0.8000D0 , &
               2.0000D0 , 0.1000D0 , 0.0500D0 , 0.0400D0 , 0.0600D0 , &
               0.1000D0 , 0.0100D0 , 0.0300D0 , 0.0004D0 , 0.0004D0 , &
               0.1000D0 , 0.1000D0 , 0.8000D0 , 0.3000D0 , 0.3000D0 , &
               1.5000D0 , 0.4000D0 /
  !
  ! displacement height (meter)
  ! if great parts of veg. are covered by snow, use displa=0
  ! because then the new displa-theory is not valid
  !
  data displa / 0.00D0 , 0.00D0 , 9.00D0 , 9.00D0 , 0.00D0 ,18.00D0 , &
                0.00D0 , 0.00D0 , 0.00D0 , 0.00D0 , 0.00D0 , 0.00D0 , &
                0.00D0 , 0.00D0 , 0.00D0 , 0.00D0 , 0.00D0 , 0.00D0 , &
                0.00D0 , 0.00D0 , 6.00D0 , 2.50D0 /
  !
  ! minimum stomatal resistance (s/m)
  ! data rsmin/153.0,4*200.0,150.0,14*200.0/   ! shuttleworth
  ! data rsmin/120.0,4*200.0,150.0,14*200.0/   ! bats1e numbers
  !Sara
  !     data rsmin/45.,60.,2*80.,120.,2*60.,200.,80.,45.,150.,200.,45.
  !     &          ,2*200.,80.,120.,100.,2*120./
  ! (increase in rsmin will lead to a decrease in evapotranspration)
  ! Model is expecially sensible to this, tune it carefully
  ! BATS 1e , Table 2, i, pag 21
  !
  !data rsmin / 120.0D0 , 200.0D0 , 200.0D0 , 200.0D0 , 200.0D0 , 150.0D0 , &
  !             200.0D0 , 200.0D0 , 200.0D0 , 200.0D0 , 200.0D0 , 200.0D0 , &
  !             200.0D0 , 200.0D0 , 200.0D0 , 200.0D0 , 200.0D0 , 200.0D0 , &
  !             200.0D0 , 200.0D0 , 120.0D0 ,  60.0D0 /
  data rsmin /  45.0D0 , 200.0D0 , 200.0D0 , 200.0D0 , 200.0D0 ,  50.0D0 , &
               200.0D0 , 200.0D0 , 200.0D0 , 200.0D0 , 200.0D0 , 200.0D0 , &
               200.0D0 , 200.0D0 , 200.0D0 , 200.0D0 , 200.0D0 , 200.0D0 , &
               200.0D0 , 200.0D0 , 120.0D0 ,  60.0D0 /
  !
  ! maximum leaf area index (ratio unit cover per unit ground)
  ! BATS 1e , Table 2, j, pag 21
  !
  !data xla / 6.0D0 , 2.0D0 , 6.0D0 , 6.0D0 , 6.0D0 , 6.0D0 , &
  !           6.0D0 , 0.0D0 , 6.0D0 , 6.0D0 , 6.0D0 , 0.0D0 , &
  !           6.0D0 , 0.0D0 , 0.0D0 , 6.0D0 , 6.0D0 , 6.0D0 , &
  !           6.0D0 , 6.0D0 , 1.0D0, 2.0D0 /
  data xla / 6.0D0 , 2.0D0 , 6.0D0 , 6.0D0 , 6.0D0 , 6.0D0 , &
             3.0D0 , 0.0D0 , 2.0D0 , 4.0D0 , 1.0D0 , 0.0D0 , &
             4.0D0 , 0.0D0 , 0.0D0 , 4.0D0 , 4.0D0 , 5.0D0 , &
             4.0D0 , 1.0D0 , 1.0D0 , 2.0D0 /
  !
  ! minimum leaf area index **lai depends on temp as veg cover
  ! BATS 1e , Table 2, k, pag 21
  !
  !data xlai0 / 0.5D0 , 0.5D0 , 5.0D0 , 1.0D0 , 1.0D0 , 5.0D0 , &
  !             0.5D0 , 0.0D0 , 0.5D0 , 0.5D0 , 0.5D0 , 0.0D0 , &
  !             0.5D0 , 0.0D0 , 0.0D0 , 5.0D0 , 1.0D0 , 3.0D0 , &
  !             0.5D0 , 0.5D0 , 0.5D0 , 0.5D0 /
  data xlai0 / 0.5D0 , 0.5D0 , 5.0D0 , 1.0D0 , 1.0D0 , 5.0D0 , &
               1.0D0 , 0.0D0 , 0.5D0 , 2.0D0 , 0.5D0 , 0.0D0 , &
               2.0D0 , 0.0D0 , 0.0D0 , 3.0D0 , 1.0D0 , 3.0D0 , &
               0.5D0 , 1.0D0 , 0.5D0 , 1.0D0 /
  !
  ! stem and dead matter area index (projected area of non-transpiring sfcs)
  ! BATS 1e , Table 2, l, pag 21
  !
  data sai / 0.5D0 , 4.0D0 , 2.0D0 , 2.0D0 , 2.0D0 , 2.0D0 , &
             2.0D0 , 0.5D0 , 0.5D0 , 2.0D0 , 2.0D0 , 2.0D0 , &
             2.0D0 , 2.0D0 , 2.0D0 , 2.0D0 , 2.0D0 , 2.0D0 , &
             2.0D0 , 2.0D0 , 0.5D0 , 0.5D0 /
  !
  ! inverse square root of leaf dimension - used for calculating fluxes
  ! from foliage (m**0.5)
  ! BATS 1e , Table 2, m, pag 21
  !
  data sqrtdi / 10.0D0 , 5.0D0 , 5.0D0 , 5.0D0 , 5.0D0 , 5.0D0 , &
                 5.0D0 , 5.0D0 , 5.0D0 , 5.0D0 , 5.0D0 , 5.0D0 , &
                 5.0D0 , 5.0D0 , 5.0D0 , 5.0D0 , 5.0D0 , 5.0D0 , &
                 5.0D0 , 5.0D0 , 5.0D0 , 5.0D0 /
  !
  ! fc = light dependence of stomatal resistance (m^2/W)
  ! BATS 1e , Table 2, n, pag 21
  !
  data fc / 0.02D0 , 0.02D0 , 0.06D0 , 0.06D0 , 0.06D0 , 0.06D0 , &
            0.02D0 , 0.02D0 , 0.02D0 , 0.02D0 , 0.02D0 , 0.02D0 , &
            0.02D0 , 0.02D0 , 0.02D0 , 0.02D0 , 0.02D0 , 0.06D0 , &
            0.02D0 , 0.02D0 , 0.02D0 , 0.02D0 /
  !
  ! depuv is depth of upper soil layer (mm)
  ! BATS 1e , Table 2, e, pag 21
  !
  data depuv / 100.0D0 , 100.0D0 , 100.0D0 , 100.0D0 , 100.0D0 , &
               100.0D0 , 100.0D0 , 100.0D0 , 100.0D0 , 100.0D0 , &
               100.0D0 , 100.0D0 , 100.0D0 , 100.0D0 , 100.0D0 , &
               100.0D0 , 100.0D0 , 100.0D0 , 100.0D0 , 100.0D0 , &
               100.0D0 , 100.0D0 /
  !
  ! deprv is depth of root zone (mm)
  ! BATS 1e , Table 2, d, pag 21
  !
  data deprv / 1000.0D0 , 1000.0D0 , 1500.0D0 , 1500.0D0 , 2000.0D0 , &
               1500.0D0 , 1000.0D0 , 1000.0D0 , 1000.0D0 , 1000.0D0 , &
               1000.0D0 , 1000.0D0 , 1000.0D0 , 1000.0D0 , 1000.0D0 , &
               1000.0D0 , 1000.0D0 , 2000.0D0 , 1000.0D0 , 1000.0D0 , &
               1000.0D0 , 1000.0D0 /
  !
  ! deptv is depth of total soil (mm)
  !
  data deptv / 3000.0D0 , 3000.0D0 , 3000.0D0 , 3000.0D0 , 3000.0D0 , &
               3000.0D0 , 3000.0D0 , 3000.0D0 , 3000.0D0 , 3000.0D0 , &
               3000.0D0 , 3000.0D0 , 3000.0D0 , 3000.0D0 , 3000.0D0 , &
               3000.0D0 , 3000.0D0 , 3000.0D0 , 3000.0D0 , 3000.0D0 , &
               3000.0D0 , 3000.0D0 /
  !
  ! iexsol is soil texture type (see subr soilbc)
  !
  !data iexsol / 6 , 6 , 6 , 6 , 7 , 8 , 6 , 3 , 6 , 6 , 5 , 12 , 6 ,  &
  !              6 , 6 , 6 , 5 , 6 , 6 , 6 , 12 , 8 /
  data iexsol / 6 , 6 , 6 , 6 , 7 , 8 , 6 , 1 , 6 , 6 , 5 , 12 , 6 ,  &
                6 , 6 , 6 , 5 , 6 , 6 , 6 , 12 , 8 /
  !
  ! kolsol is soil color type (see subr. albedo)
  !
  !data kolsol / 6 , 4 , 5 , 5 , 5 , 5 , 5 , 3 , 4 , 4 , 2 , 1 , 6 ,   &
  !              6 , 6 , 5 , 4 , 5 , 5 , 5 , 4 , 4 /
  data kolsol / 6 , 4 , 5 , 5 , 5 , 5 , 5 , 1 , 4 , 4 , 2 , 1 , 6 ,   &
                6 , 6 , 5 , 4 , 5 , 5 , 5 , 4 , 4 /
  !
  ! slmo is initial surface moisture availability in fraction of one
  !
  data slmo / 0.65D0 , 0.45D0 , 0.60D0 , 0.60D0 , 0.65D0 , &
              0.65D0 , 0.55D0 , 0.10D0 , 0.90D0 , 0.80D0 , &
              0.20D0 , 0.90D0 , 0.90D0 , 1.00D0 , 1.00D0 , &
              0.50D0 , 0.50D0 , 0.65D0 , 0.60D0 , 0.60D0 , &
              0.50D0 , 0.60D0 /
  !
  ! xmopor is fraction of soil that is voids
  ! BATS 1e , Table 3, I , a, pag 27
  !
  data xmopor / 0.33D0 , 0.36D0 , 0.39D0 , 0.42D0 , 0.45D0 , 0.48D0 , &
                0.51D0 , 0.54D0 , 0.57D0 , 0.60D0 , 0.63D0 , 0.66D0 /
  !
  ! xmosuc is the minimum soil suction (mm)
  ! BATS 1e , Table 3, I , b, pag 27
  !
  data xmosuc / 30.0D0 ,  30.0D0 ,  30.0D0 , 200.0D0 , 200.0D0 , 200.0D0 , &
               200.0D0 , 200.0D0 , 200.0D0 , 200.0D0 , 200.0D0 , 200.0D0 /
  !
  ! xmohyd is the max. hydraulic conductivity (mm/s)
  ! BATS 1e , Table 3, I , c, pag 27
  !
  data xmohyd / 0.2000D0 , 0.0800D0 , 0.0320D0 , 0.0130D0 , 0.0089D0 , &
                0.0063D0 , 0.0045D0 , 0.0032D0 , 0.0022D0 , 0.0016D0 , &
                0.0011D0 , 0.0008D0 /
  !
  ! xmowilt is fraction of water content at which permanent wilting occurs
  ! (transpiration ceases)
  ! BATS 1e , Table 3, I , f, pag 27
  !
  !data xmowil / 0.088D0 , 0.119D0 , 0.151D0 , 0.266D0 , 0.300D0 , &
  !              0.332D0 , 0.378D0 , 0.419D0 , 0.455D0 , 0.487D0 , &
  !              0.516D0 , 0.542D0 /
  data xmowil / 0.095D0 , 0.128D0 , 0.161D0 , 0.266D0 , 0.300D0 , &
                0.332D0 , 0.378D0 , 0.419D0 , 0.455D0 , 0.487D0 , &
                0.516D0 , 0.542D0 /
  !
  ! Field capacity (Lara Kueppers modif for irrigated agricultural lands)
  ! see Kueppers et al. (2008)
  !
  !data xmofc / 0.404D0 , 0.477D0 , 0.547D0 , 0.614D0 , 0.653D0 ,      &
  !             0.688D0 , 0.728D0 , 0.763D0 , 0.794D0 , 0.820D0 ,      &
  !             0.845D0 , 0.866D0 /
  !
  ! bee is the clapp and hornbereger "b" parameter
  ! BATS 1e , Table 3, I , e, pag 27
  !
  data bee / 3.5D0 , 4.0D0 , 4.5D0 , 5.0D0 , 5.5D0 , 6.0D0 , 6.8D0 ,  &
             7.6D0 , 8.4D0 , 9.2D0 ,10.0D0 ,10.8D0 /
  !
  ! bskrat is ratio of soil thermal conduc. to that of loam
  ! a function of texture
  ! BATS 1e , Table 3, I , d, pag 27
  !
  data skrat / 1.70D0 , 1.50D0 , 1.30D0 , 1.20D0 , 1.10D0 , 1.00D0 , &
               0.95D0 , 0.90D0 , 0.85D0 , 0.80D0 , 0.75D0 , 0.70D0 /
  !
  ! albvgs is vegetation albedo for wavelengths < 0.7 microns data
  ! BATS 1e , Table 2, I , g, pag 21
  !
  !data albvgs / 0.10D0 , 0.10D0 , 0.05D0 , 0.05D0 , 0.08D0 , 0.04D0 , &
  !              0.08D0 , 0.20D0 , 0.10D0 , 0.08D0 , 0.17D0 , 0.80D0 , &
  !              0.06D0 , 0.07D0 , 0.07D0 , 0.05D0 , 0.08D0 , 0.06D0 , &
  !              0.06D0 , 0.06D0 , 0.02D0 , 0.06D0 /
  data albvgs / 0.10D0 , 0.10D0 , 0.04D0 , 0.04D0 , 0.06D0 , 0.04D0 , &
                0.08D0 , 0.20D0 , 0.10D0 , 0.08D0 , 0.17D0 , 0.80D0 , &
                0.06D0 , 0.07D0 , 0.07D0 , 0.05D0 , 0.08D0 , 0.05D0 , &
                0.06D0 , 0.06D0 , 0.02D0 , 0.06D0 /
  !
  ! albvgl is vegetation albedo for wavelengths > 0.7 microns data
  ! BATS 1e , Table 2, I , h, pag 21
  !
  !data albvgl / 0.30D0 , 0.30D0 , 0.23D0 , 0.23D0 , 0.28D0 , 0.20D0 , &
  !              0.30D0 , 0.40D0 , 0.30D0 , 0.28D0 , 0.34D0 , 0.60D0 , &
  !              0.18D0 , 0.20D0 , 0.20D0 , 0.23D0 , 0.28D0 , 0.24D0 , &
  !              0.18D0 , 0.18D0 , 0.15D0 , 0.18D0 /
  data albvgl / 0.30D0 , 0.30D0 , 0.20D0 , 0.20D0 , 0.26D0 , 0.20D0 , &
                0.30D0 , 0.40D0 , 0.30D0 , 0.28D0 , 0.34D0 , 0.60D0 , &
                0.18D0 , 0.20D0 , 0.20D0 , 0.23D0 , 0.28D0 , 0.23D0 , &
                0.18D0 , 0.18D0 , 0.15D0 , 0.18D0 /
  !
  ! Fraction of water extracted by upper layer roots (saturated)
  ! BATS 1e , Table 2, I , h, pag 21
  !
  data rootf / 0.30D0 , 0.80D0 , 0.67D0 , 0.67D0 , 0.50D0 , 0.80D0 , &
               0.80D0 , 0.90D0 , 0.90D0 , 0.30D0 , 0.80D0 , 0.50D0 , &
               0.50D0 , 0.50D0 , 0.50D0 , 0.50D0 , 0.50D0 , 0.50D0 , &
               0.50D0 , 0.50D0 , 0.90D0 , 0.50D0 /
  !
  ! Emissivity coefficients (used in iemiss = 1)
  !
  data lndemiss / 0.983D0, 0.983D0, 0.983D0, 0.987D0, 0.981D0, 0.981D0, &
                  0.983D0, 0.965D0, 0.987D0, 0.985D0, 0.970D0, 0.993D0, &
                  0.992D0, 0.992D0, 0.992D0, 0.983D0, 0.972D0, 0.983D0, &
                  0.981D0, 0.991D0, 0.970D0, 0.972D0 /
  !
  ! Seasonal variations (used in iemiss = 1)
  !
  data seasemi / 0.005D0, 0.002D0, 0.000D0, 0.004D0, 0.004D0, 0.000D0, &
                 0.002D0, 0.000D0, 0.000D0, 0.002D0, 0.000D0, 0.000D0, &
                 0.000D0, 0.000D0, 0.000D0, 0.000D0, 0.004D0, 0.002D0, &
                 0.004D0, 0.000D0, 0.000D0, 0.001D0 /
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

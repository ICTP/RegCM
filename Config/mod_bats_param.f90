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
      implicit none
!
      real(8) , dimension(8) :: solour
      real(8) , dimension(22) :: albvgl , albvgs , crough , deprv ,     &
                               & deptv , depuv , displa , fc , freza ,  &
                               & frezu , rough , rsmin , sai , seasf ,  &
                               & sqrtdi , vegc , xla , xlai0 , rootf ,  &
                               & aest , arye
      real(8) , dimension(12) :: bee , skrat , xmofc , xmohyd , xmopor ,&
                               & xmosuc , xmowil
      integer , dimension(22) :: iexsol , kolsol
!
!------------------ DATA SECTION ----------------------------------------
!

!*    soil albedo for different coloured
     data solour / 0.16D0 , 0.15D0 , 0.10D0 , 0.09D0 , 0.08D0 ,        &
                   0.07D0 , 0.06D0 , 0.05D0 /

!*    vegc is maximum fractional cover of vegetation
      data vegc /0.85D0 , 0.8D0 , 0.8D0 , 0.8D0 , 0.8D0  , 0.9D0 ,      &
              &  0.8D0  , 0.0D0 , 0.6D0 , 0.8D0 , 0.35D0 , 0.0D0 ,      &
              &  0.8D0  , 2*0.0D0 , 5*0.8D0, 0.05D0, 0.40D0 /
!*    seasf is the difference between vegc and fractional cover at 269k
      data seasf /0.6D0 , 0.1D0 , 0.1D0 , 0.3D0 , 0.3D0 , 0.5D0 ,       &
              &   0.3D0 , 0.0D0 , 0.2D0 , 0.6D0 , 0.1D0 , 0.0D0 ,       &
              &   0.4D0 , 0.0D0 , 0.0D0 , 0.2D0 , 0.3D0 , 0.2D0 ,       &
              &   2*0.4D0, 0.05D0, 0.15D0 /
!*    rough is an aerodynamic roughness length (m) =approx 0.1*veg
!*    height also used snow masking depth in subrout albedo
      data rough /0.08D0 , 0.05D0 , 2*1.0D0 , 0.8D0 , 2.0D0  , 0.1D0  , &
              &   0.05D0 , 0.04D0 , 0.06D0 ,  0.1D0 , 0.01D0 , 0.03D0 , &
              &   2*0.0004D0 , 2*0.1D0 , 0.8D0 , 2*0.3D0, 1.5D0, 0.40D0 /
!     ******      displacement height (meter)
!     ******      if great parts of veg. are covered by snow, use
!     displa=0 ******      because mod_then the new displa-theory is not
!     valid
      data displa/0.0D0 , 0.0D0 , 9.0D0 , 9.0D0 , 0.0D0 , 18.0D0 ,      &
              &   14*0.0D0, 6.0D0, 2.5D0/
!     ******      min stomatl resistance (s/m)
!cc   data rsmin/153.0,4*200.0,150.0,14*200.0/   ! shuttleworth
!     data rsmin/120.0,4*200.0,150.0,14*200.0/   ! bats1e numbers
!Sara
!     data rsmin/45.,60.,2*80.,120.,2*60.,200.,80.,45.,150.,200.,45.
!     &          ,2*200.,80.,120.,100.,2*120./
      data rsmin /5*200.0D0,50.0D0,14*200.0D0, 120.0D0, 60.0D0/
!Sara_
!     ******      max leaf area index (ratio unit cover per unit ground)
!     ORIGINAL
!     data xla/6.,2.,5*6.,0.,3*6.,0.,6.,2*0.,5*6./
!     Laura 21/04/08
      data xla/4.0D0 , 2.0D0 , 4*6.0D0 , 3.0D0 , 0.0D0 , 2.0D0 , 4.0D0 ,&
            &  1.0D0 , 0.0D0 , 4.0D0 , 2*0.0D0 , 4.0D0 , 4.0D0 , 5.0D0 ,&
            &  4.0D0 , 1.0D0, 1.0D0, 2.0D0 /
!     ******      min leaf area index **lai depends on temp as veg cover
!     ORIGINAL
!     data xlai0/0.5,0.5,5.,2*1.,5.,0.5,0.,3*0.5,0.,0.5,2*0.,5.,1.,3.,
!     &2*0.5/
!     Laura 21/04/08
      data xlai0/2*0.5D0 , 5.0D0 , 2*1.0D0 , 5.0D0 , 1.0D0 , 0.0D0 ,    &
            &    0.5D0 , 2.0D0 , 0.5D0 , 0.0D0 , 2.0D0 , 2*0.0D0 ,      &
            &    3.0D0 , 1.0D0 , 3.0D0 , 0.5D0 , 1.0D0, 0.5D0, 1.0D0/

!     ******      stem area index (projected area of non-transpiring
!     sfcs)
      data sai/0.5D0 , 4.0D0 , 5*2.0D0 , 2*0.5D0 , 11*2.0D0, 2*0.5D0/
!     ******      inverse square root of leaf dimension - used for
!     ******      calculating fluxes from foliage
      data sqrtdi/10.0D0 , 19*5.0D0, 2*5.0D0/
!     ******      fc = light dependence of stomatal resistance
      data fc/0.02D0 , 0.02D0 , 4*0.06D0 , 11*0.02D0 , 0.06D0 ,         &
          &   2*0.02D0, 2*0.02D0/

!     ******      depuv is depth of upper soil layer (mm)
!     ******      deprv is depth of root zone (mm)
!     ******      deptv is depth of total soil (mm)
      data depuv/22*100.0D0/
      data deprv/2*1000.0D0 , 2*1500.0D0 , 2000.0D0 , 1500.0D0 ,        &
           &     11*1000.D0 , 2000.D0 , 2*2000.D0, 2*1000.0D0/
      data deptv/20*3000.D0, 2*3000.0D0/
!     ******      iexsol is soil texture type (see subr soilbc)
!     ORIGINAL
!     data iexsol/6,6,6,6,7,8,6,3,6,6,5,12,6,6,6,6,5,6,6,6/
!     Laura  04/04/08 changed soil texture for desert: 3->1
      data iexsol/6 , 6 , 6 , 6 , 7 , 8 , 6 , 1 , 6 , 6 , 5 , 12 , 6 ,  &
              &   6 , 6 , 6 , 5 , 6 , 6 , 6 , 12 , 8/
!     ******      kolsol is soil color type (see subr. albedo)
!     Dec. 15, 2008
!     data kolsol/5,3,4,4,4,4,4,1,3,3,2,1,5,5,5,4,3,4,4,4/
      data kolsol/6 , 4 , 5 , 5 , 5 , 5 , 5 , 1 , 4 , 4 , 2 , 1 , 6 ,   &
             &    6 , 6 , 5 , 4 , 5 , 5 , 5 , 4 , 4/
!     Dec. 15, 2008_
!     ******      xmopor is fraction of soil that is voids
      data xmopor/0.33D0 , 0.36D0 , 0.39D0 , 0.42D0 , 0.45D0 , 0.48D0 , &
             &    0.51D0 , 0.54D0 , 0.57D0 , 0.6D0 , 0.63D0 , 0.66D0/
!     ******      xmosuc is the minimum soil suction (mm)
      data xmosuc/3*30.0D0 , 9*200.D0/
!     ******      xmohyd is the max. hydraulic conductivity (mm/s)
      data xmohyd/0.20D-0 , 0.80D-1 , 0.32D-1 , 0.13D-1 , 0.89D-2 ,     &
         & 0.63D-2 , 0.45D-2 , 0.32D-2 , 0.22D-2 , 0.16D-2 , 0.11D-2 ,  &
         & 0.80D-3/
!     ******      xmowilt is fraction of water content at which
!     permanent wilting occurs
      data xmowil/0.095D0 , 0.128D0 , 0.161D0 , 0.266D0 , 0.3D0 ,       &
         &        0.332D0 , 0.378D0 , 0.419D0 , 0.455D0 , 0.487D0 ,     &
         &        0.516D0 , 0.542D0 /
      data xmofc/0.404D0 , 0.477D0 , 0.547D0 , 0.614D0 , 0.653D0 ,      &
         &       0.688D0 , 0.728D0 , 0.763D0 , 0.794D0 , 0.820D0 ,      &
         &       0.845D0 , 0.866D0/
!     ******      bee is the clapp and hornbereger "b" parameter
      data bee/3.5D0 , 4.0D0 , 4.5D0 , 5.0D0 , 5.5D0 , 6.0D0 , 6.8D0 ,  &
         &     7.6D0 , 8.4D0 , 9.2D0 , 10.0D0 , 10.8D0/
!     ******      bskrat is ratio of soil thermal conduc. to that of
!     loam - a function of texture
      data skrat/1.7D0 , 1.5D0 , 1.3D0 , 1.2D0 , 1.1D0 , 1.0D0 , .95D0 ,&
         &       0.9D0 , 0.85D0 , 0.80D0 , 0.75D0 , 0.7D0/

!     Dec. 15, 2008
!     ******      albvgs is vegetation albedo for wavelengths < 0.7
!     microns data
!     albvgs/.1,.1,.05,.05,.08,.04,.08,.2,.1,.08,.17,.8,.06,2*.07,
!     &.05,.08,.06,2*0.06/
      data albvgs/.1D0 , .1D0 , .04D0 , .04D0 , .06D0 , .04D0 , .08D0 , &
          &       .2D0 , .1D0 , .08D0 , .17D0 , .8D0 , .06D0 , 2*.07D0 ,&
          &       .05D0 , .08D0 , .05D0 , 2*0.06D0, 0.02D0, 0.06D0/
!     ******      albvgl is vegetation albedo for wavelengths > 0.7
!     microns data
!     albvgl/.3,.3,.23,.23,.28,.20,.30,.4,.3,.28,.34,.6,.18,2*.2,
!     &.23,.28,.24,2*.18/
      data albvgl/.3D0 , .3D0 , .20D0 , .20D0 , .26D0 , .20D0 , .30D0 , &
         &        .4D0 , .3D0 , .28D0 , .34D0 , .6D0 , .18D0 , 2*.2D0 , &
         &        .23D0 , .28D0 , .23D0 , 2*.18D0, 0.15D0, 0.18D0/
!     Dec. 15, 2008_

      data rootf/0.30D0 , 0.80D0 , 0.67D0 , 0.67D0 , 0.50D0 , 0.80D0 , &
                 0.80D0 , 0.90D0 , 0.90D0 , 0.30D0 , 0.80D0 , 9*0.50D0,&
                 0.90D0 , 0.50D0/
!
! Stokes parameters
!
      data aest     /0.80D0 , 0.80D0 , 0.8D0 , 0.8D0 , 1.2D0 , 1.20D0 , &
           2.0D0 , 1.5D0 ,  1.5D0 , 2.0D0 , 15.0D0 , 15.0D0 , 1.5D0 ,   &
           1.5D0 , 1.5D0 , 15.0D0 , 1.2D0 , 1.2D0 , 1.2D0 , 1.2D0 ,     &
           1.2D0 , 1.2D0 /
!
      data arye     /0.5D0 , 5.0D0 , 0.5D0 , 5.0D0 , 1.0D0 , 1.0D0 ,    &
         0.0001D0 , 5.0D0 , 10.0D0 , 10.0D0 , 0.0001D0 , 0.0001D0 ,     &
         0.56D0 , 0.56D0 , 0.56D0 , 0.56D0 ,  0.56D0 , 0.56D0 , 0.56D0 ,&
         0.56D0 , 1.0D0 , 1.0D0/
!
      end module mod_bats_param

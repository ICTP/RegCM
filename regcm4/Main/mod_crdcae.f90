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
!    MDRCHANTABILITY or FITNDSS FOR A PARTICULAR PURPOSD.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      module mod_crdcae
      implicit none
!
      real(8) , dimension(2) :: a1 , a2 , b1 , b2 , realk , st
      real(8) , dimension(4) :: c1 , c2 , c3 , c4 , c5 , c6 , c7
      real(8) :: c10 , c11 , c12 , c13 , c14 , c15 , c16 , c17 , c18 ,  &
               & c19 , c20 , c21 , c22 , c23 , c24 , c25 , c26 , c27 ,  &
               & c28 , c29 , c30 , c31 , c8 , c9 , cfa1 , fc1 , fwc1 ,  &
               & fwc2 , fwcoef
      real(8) :: co2vmr
      real(8) , dimension(3,4) :: coefa , coefc , coefe
      real(8) , dimension(4,4) :: coefb , coefd
      real(8) , dimension(6,2) :: coeff , coefi
      real(8) , dimension(2,4) :: coefg , coefh
      real(8) , dimension(3,2) :: coefj , coefk
!
!     H2O DMISSIVITY AND ABSORTIVITY CODFFICIDNTS
!
      data coefa/1.01400D+00 , 6.41695D-03 , 2.85787D-05 , 1.01320D+00 ,&
         & 6.86400D-03 , 2.96961D-05 , 1.02920D+00 , 1.01680D-02 ,      &
         & 5.30226D-05 , 1.02743D+00 , 9.85113D-03 , 5.00233D-05/
!
      data coefb/8.85675D+00 , -3.51620D-02 , 2.38653D-04 ,             &
         & -1.71439D-06 , 5.73841D+00 , -1.91919D-02 , 1.65993D-04 ,    &
         & -1.54665D-06 , 6.64034D+00 , 1.56651D-02 , -9.73357D-05 ,    &
         & 0.0 , 7.09281D+00 , 1.40056D-02 , -1.15774D-04 , 0.0/
!
      data coefc/9.90127D-01 , 1.22475D-03 , 4.90135D-06 , 9.89753D-01 ,&
         & 1.97081D-03 , 3.42046D-06 , 9.75230D-01 , 1.03341D-03 , 0.0 ,&
         & 9.77366D-01 , 8.60014D-04 , 0.0/
!
      data coefd/7.03047D-01 , -2.63501D-03 , -1.57023D-06 , 0.0 ,      &
         & 5.29269D-01 , -3.14754D-03 , 4.39595D-06 , 0.0 ,             &
         & 7.88193D-02 , 1.31290D-03 , 4.25827D-06 , -1.23982D-08 ,     &
         & 1.62744D-01 , 2.22847D-03 , 2.60102D-06 , -4.30133D-08/
!
      data coefe/3.93137D-02 , -4.34341D-05 , 3.74545D-07 ,             &
         & 3.67785D-02 , -3.10794D-05 , 2.94436D-07 , 7.42500D-02 ,     &
         & 3.97397D-05 , 0.0 , 7.52859D-02 , 4.18073D-05 , 0.0/
!
      data coeff/2.2037D-01 , 1.39719D-03 , -7.32011D-06 ,              &
         & -1.40262D-08 , 2.13638D-10 , -2.35955D-13 , 3.07431D-01 ,    &
         & 8.27225D-04 , -1.30067D-05 , 3.49847D-08 , 2.07835D-10 ,     &
         & -1.98937D-12/
!
      data coefg/9.04489D+00 , -9.56499D-03 , 1.80898D+01 ,             &
         & -1.91300D-02 , 8.72239D+00 , -9.53359D-03 , 1.74448D+01 ,    &
         & -1.90672D-02/
!
      data coefh/5.46557D+01 , -7.30387D-02 , 1.09311D+02 ,             &
         & -1.46077D-01 , 5.11479D+01 , -6.82615D-02 , 1.02296D+02 ,    &
         & -1.36523D-01/
!
      data coefi/3.31654D-01 , -2.86103D-04 , -7.87860D-06 ,            &
         & 5.88187D-08 , -1.25340D-10 , -1.37731D-12 , 3.14365D-01 ,    &
         & -1.33872D-03 , -2.15585D-06 , 6.07798D-08 , -3.45612D-10 ,   &
         & -9.34139D-15/
!
      data coefj/2.82096D-02 , 2.47836D-04 , 1.16904D-06 , 9.27379D-02 ,&
         & 8.04454D-04 , 6.88844D-06/
!
      data coefk/2.48852D-01 , 2.09667D-03 , 2.60377D-06 , 1.03594D+00 ,&
         & 6.58620D-03 , 4.04456D-06/
!
!
!     Narrow band data for H2O
!     200CM data for 800-1000 CM-1 and 1000-1200 CM-1.
!
      data realk/0.18967069430426D-04 , 0.70172244841851D-04/
      data st/0.31930234492350D-03 , 0.97907319939060D-03/
      data a1/0.28775403075736D-01 , 0.23236701470511D-01/
      data a2/ -0.57966222388131D-04 , -0.95105504388411D-04/
      data b1/0.29927771523756D-01 , 0.21737073577293D-01/
      data b2/ -0.86322071248593D-04 , -0.78543550629536D-04/

      end module mod_crdcae

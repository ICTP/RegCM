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

module physics_msis

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_spline

  implicit none

  private

  character(4) , public , dimension(2) :: chname , istime
  character(4) , public , dimension(3) :: isdate
  real(rkx) , public , dimension(25) :: sav
  integer(ik4) , public :: imr = 0
  real(rkx) , public :: tlb , s , db04 , db16 , db28 , db32 , db40 , &
              db48 , db01 , za , t0 , z0 , xg0 , rl , dd , db14 , tr12

  public :: gtd7d

  public :: gtd7 , gts7 , ghp7
  public :: tselec , tretrv

  integer(ik4) :: isw
  real(rkx) , dimension(25) :: sw , swc
  real(rkx) :: dm01 , dm04 , dm14 , dm16 , dm28 , dm32 , dm40
  real(rkx) , dimension(10,8) :: pdm
  real(rkx) , dimension(10) :: ptm
  real(rkx) :: apd , apdf , c2tloc , c3tloc , ctloc , day , df , dfa ,   &
              s2tloc , s3tloc , stloc
  real(rkx) , dimension(10) :: pavgm
  real(rkx) , dimension(4) :: apt
  real(rkx) , dimension(2) :: tgn1 , tgn2 , tgn3
  real(rkx) , dimension(5) :: tn1 , tn3
  real(rkx) , dimension(4) :: tn2
  integer(ik4) :: iyr
  real(rkx) , dimension(9,4) :: plg

  real(rkx) , dimension(150) , target :: pt      ! pt1 , pt2 , pt3
  real(rkx) , dimension(150,9) , target :: pd    ! pa1 , pa2 , pa3
                                                ! pb1 , pb2 , pb3
                                                ! pc1 , pc2 , pc3
                                                ! pd1 , pd2 , pd3
                                                ! pe1 , pe2 , pe3
                                                ! pf1 , pf2 , pf3
                                                ! pg1 , pg2 , pg3
                                                ! ph1 , ph2 , ph3
                                                ! pi1 , pi2 , pi3
  real(rkx) , dimension(150) , target :: ps      ! pj1 , pj2 , pj3
  real(rkx) , dimension(25,2) , target :: pdl    ! pk1
  real(rkx) , dimension(100,4) , target :: ptl   ! pl1 , pl2
                                                ! pm1 , pm2
                                                ! pn1 , pn2
                                                ! po1 , po2
  real(rkx) , dimension(100,10) , target :: pma  ! pp1 , pp2
                                                ! pq1 , pq2
                                                ! pr1 , pr2
                                                ! ps1 , ps2
                                                ! pu1 , pu2
                                                ! pv1 , pv2
                                                ! pw1 , pw2
                                                ! px1 , px2
                                                ! py1 , py2
                                                ! pz1 , pz2
  real(rkx) , dimension(100) , target :: sam     ! paa1 , paa2

  real(rkx) :: gsurf , re
  real(rkx) :: gb , rout , tinf
  real(rkx) , dimension(15) :: t
  ! Gas constant in J/K/mol
  real(rkx) , parameter :: r100gas = rgasmol*d_100
  real(rkx) , parameter :: nearzero = 0.000001_rkx
  real(rkx) , parameter :: dr = 1.72142e-2_rkx
  integer(ik4) , parameter :: mn1 = 5
  integer(ik4) , parameter :: mn2 = 4
  integer(ik4) , parameter :: mn3 = 5
  real(rkx) , dimension(mn1) :: zn1
  real(rkx) , dimension(mn2) :: zn2
  real(rkx) , dimension(mn3) :: zn3
  integer(ik4) , dimension(11) :: mt
  real(rkx) , dimension(9) :: alph
  real(rkx) , dimension(8) :: altl

  ! MSISE-00 01-FEB-02
  data isdate/'01-F' , 'EB-0' , '2   '/
  data istime/'15:4' , '9:27'/
  data chname/'MSIS' , 'E-00'/
  !
  data zn1 /120.0_rkx , 110.0_rkx , 100.0_rkx , 90.0_rkx , 72.50_rkx/
  data zn2 /72.5_rkx , 55.0_rkx , 45.0_rkx , 32.50_rkx/
  data zn3 /32.5_rkx , 20.0_rkx , 15.0_rkx , 10.0_rkx , 0.0_rkx/
  !
  data mt /48 , 0 , 4 , 16 , 28 , 32 , 40 , 1 , 49 , 14 , 17/
  data altl /200.0_rkx , 300.0_rkx , 160.0_rkx , 250.0_rkx , &
             240.0_rkx , 450.0_rkx , 320.0_rkx , 450.0_rkx/
  data alph / -0.380_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx , &
                0.170_rkx , 0.0_rkx , -0.380_rkx , 0.0_rkx , 0.0_rkx/

  ! TEMPERATURE
  data pt(1:50) / &
            9.86573e-1_rkx ,  1.62228e-2_rkx ,  1.55270e-2_rkx , -1.04323e-1_rkx , &
           -3.75801e-3_rkx , -1.18538e-3_rkx , -1.24043e-1_rkx ,  4.56820e-3_rkx , &
            8.76018e-3_rkx , -1.36235e-1_rkx , -3.52427e-2_rkx ,  8.84181e-3_rkx , &
           -5.92127e-3_rkx , -8.61650e+0_rkx ,  0.00000e+0_rkx ,  1.28492e-2_rkx , &
            0.00000e+0_rkx ,  1.30096e+2_rkx ,  1.04567e-2_rkx ,  1.65686e-3_rkx , &
           -5.53887e-6_rkx ,  2.97810e-3_rkx ,  0.00000e+0_rkx ,  5.13122e-3_rkx , &
            8.66784e-2_rkx ,  1.58727e-1_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx , -7.27026e-6_rkx ,  0.00000e+0_rkx ,  6.74494e+0_rkx , &
            4.93933e-3_rkx ,  2.21656e-3_rkx ,  2.50802e-3_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx , -2.08841e-2_rkx , -1.79873e+0_rkx ,  1.45103e-3_rkx , &
            2.81769e-4_rkx , -1.44703e-3_rkx , -5.16394e-5_rkx ,  8.47001e-2_rkx , &
            1.70147e-1_rkx ,  5.72562e-3_rkx ,  5.07493e-5_rkx ,  4.36148e-3_rkx , &
            1.17863e-4_rkx ,  4.74364e-3_rkx/
  data pt(51:100) / &
            6.61278e-3_rkx ,  4.34292e-5_rkx ,  1.44373e-3_rkx ,  2.41470e-5_rkx , &
            2.84426e-3_rkx ,  8.56560e-4_rkx ,  2.04028e-3_rkx ,  0.00000e+0_rkx , &
           -3.15994e+3_rkx , -2.46423e-3_rkx ,  1.13843e-3_rkx ,  4.20512e-4_rkx , &
            0.00000e+0_rkx , -9.77214e+1_rkx ,  6.77794e-3_rkx ,  5.27499e-3_rkx , &
            1.14936e-3_rkx ,  0.00000e+0_rkx , -6.61311e-3_rkx , -1.84255e-2_rkx , &
           -1.96259e-2_rkx ,  2.98618e+4_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  6.44574e+2_rkx ,  8.84668e-4_rkx ,  5.05066e-4_rkx , &
            0.00000e+0_rkx ,  4.02881e+3_rkx , -1.89503e-3_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  8.21407e-4_rkx ,  2.06780e-3_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           -1.20410e-2_rkx , -3.63963e-3_rkx ,  9.92070e-5_rkx , -1.15284e-4_rkx , &
           -6.33059e-5_rkx , -6.05545e-1_rkx ,  8.34218e-3_rkx , -9.13036e+1_rkx , &
            3.71042e-4_rkx ,  0.00000e+0_rkx/
  data pt(101:150) / &
            4.19000e-4_rkx ,  2.70928e-3_rkx ,  3.31507e-3_rkx , -4.44508e-3_rkx , &
           -4.96334e-3_rkx , -1.60449e-3_rkx ,  3.95119e-3_rkx ,  2.48924e-3_rkx , &
            5.09815e-4_rkx ,  4.05302e-3_rkx ,  2.24076e-3_rkx ,  0.00000e+0_rkx , &
            6.84256e-3_rkx ,  4.66354e-4_rkx ,  0.00000e+0_rkx , -3.68328e-4_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx , -1.46870e+2_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  1.09501e-3_rkx ,  4.65156e-4_rkx ,  5.62583e-4_rkx , &
            3.21596e+0_rkx ,  6.43168e-4_rkx ,  3.14860e-3_rkx ,  3.40738e-3_rkx , &
            1.78481e-3_rkx ,  9.62532e-4_rkx ,  5.58171e-4_rkx ,  3.43731e+0_rkx , &
           -2.33195e-1_rkx ,  5.10289e-4_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           -9.25347e+4_rkx ,  0.00000e+0_rkx , -1.99639e-3_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx/
  ! HE DENSITY
  data pd(1:50,1) / &
            1.09979e+0_rkx , -4.88060e-2_rkx , -1.97501e-1_rkx ,                &
           -9.10280e-2_rkx , -6.96558e-3_rkx ,  2.42136e-2_rkx ,  3.91333e-1_rkx , &
           -7.20068e-3_rkx , -3.22718e-2_rkx ,  1.41508e+0_rkx ,  1.68194e-1_rkx , &
            1.85282e-2_rkx ,  1.09384e-1_rkx , -7.24282e+0_rkx ,  0.00000e+0_rkx , &
            2.96377e-1_rkx , -4.97210e-2_rkx ,  1.04114e+2_rkx , -8.61108e-2_rkx , &
           -7.29177e-4_rkx ,  1.48998e-6_rkx ,  1.08629e-3_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  8.31090e-2_rkx ,  1.12818e-1_rkx , -5.75005e-2_rkx , &
           -1.29919e-2_rkx , -1.78849e-2_rkx , -2.86343e-6_rkx ,  0.00000e+0_rkx , &
           -1.51187e+2_rkx , -6.65902e-3_rkx ,  0.00000e+0_rkx , -2.02069e-3_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  4.32264e-2_rkx , -2.80444e+1_rkx , &
           -3.26789e-3_rkx ,  2.47461e-3_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            9.82100e-2_rkx ,  1.22714e-1_rkx , -3.96450e-2_rkx ,  0.00000e+0_rkx , &
           -2.76489e-3_rkx ,  0.00000e+0_rkx ,  1.87723e-3_rkx/
  data pd(51:100,1) / &
           -8.09813e-3_rkx ,  4.34428e-5_rkx , -7.70932e-3_rkx ,                &
            0.00000e+0_rkx , -2.28894e-3_rkx , -5.69070e-3_rkx , -5.22193e-3_rkx , &
            6.00692e-3_rkx , -7.80434e+3_rkx , -3.48336e-3_rkx , -6.38362e-3_rkx , &
           -1.82190e-3_rkx ,  0.00000e+0_rkx , -7.58976e+1_rkx , -2.17875e-2_rkx , &
           -1.72524e-2_rkx , -9.06287e-3_rkx ,  0.00000e+0_rkx ,  2.44725e-2_rkx , &
            8.66040e-2_rkx ,  1.05712e-1_rkx ,  3.02543e+4_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx , -6.01364e+3_rkx , -5.64668e-3_rkx , &
           -2.54157e-3_rkx ,  0.00000e+0_rkx ,  3.15611e+2_rkx , -5.69158e-3_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx , -4.47216e-3_rkx , -4.49523e-3_rkx , &
            4.64428e-3_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  4.51236e-2_rkx ,  2.46520e-2_rkx ,  6.17794e-3_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx , -3.62944e-1_rkx , -4.80022e-2_rkx , &
           -7.57230e+1_rkx , -1.99656e-3_rkx ,  0.00000e+0_rkx/
  data pd(101:150,1) / &
           -5.18780e-3_rkx , -1.73990e-2_rkx , -9.03485e-3_rkx ,                &
            7.48465e-3_rkx ,  1.53267e-2_rkx ,  1.06296e-2_rkx ,  1.18655e-2_rkx , &
            2.55569e-3_rkx ,  1.69020e-3_rkx ,  3.51936e-2_rkx , -1.81242e-2_rkx , &
            0.00000e+0_rkx , -1.00529e-1_rkx , -5.10574e-3_rkx ,  0.00000e+0_rkx , &
            2.10228e-3_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , -1.73255e+2_rkx , &
            5.07833e-1_rkx , -2.41408e-1_rkx ,  8.75414e-3_rkx ,  2.77527e-3_rkx , &
           -8.90353e-5_rkx , -5.25148e+0_rkx , -5.83899e-3_rkx , -2.09122e-2_rkx , &
           -9.63530e-3_rkx ,  9.77164e-3_rkx ,  4.07051e-3_rkx ,  2.53555e-4_rkx , &
           -5.52875e+0_rkx , -3.55993e-1_rkx , -2.49231e-3_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  2.86026e+1_rkx ,  0.00000e+0_rkx ,  3.42722e-4_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx/
  ! O DENSITY
  data pd(1:50,2) / &
            1.02315e+0_rkx , -1.59710e-1_rkx , -1.06630e-1_rkx ,                &
           -1.77074e-2_rkx , -4.42726e-3_rkx ,  3.44803e-2_rkx ,  4.45613e-2_rkx , &
           -3.33751e-2_rkx , -5.73598e-2_rkx ,  3.50360e-1_rkx ,  6.33053e-2_rkx , &
            2.16221e-2_rkx ,  5.42577e-2_rkx , -5.74193e+0_rkx ,  0.00000e+0_rkx , &
            1.90891e-1_rkx , -1.39194e-2_rkx ,  1.01102e+2_rkx ,  8.16363e-2_rkx , &
            1.33717e-4_rkx ,  6.54403e-6_rkx ,  3.10295e-3_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  5.38205e-2_rkx ,  1.23910e-1_rkx , -1.39831e-2_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx , -3.95915e-6_rkx ,  0.00000e+0_rkx , &
           -7.14651e-1_rkx , -5.01027e-3_rkx ,  0.00000e+0_rkx , -3.24756e-3_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  4.42173e-2_rkx , -1.31598e+1_rkx , &
           -3.15626e-3_rkx ,  1.24574e-3_rkx , -1.47626e-3_rkx , -1.55461e-3_rkx , &
            6.40682e-2_rkx ,  1.34898e-1_rkx , -2.42415e-2_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  6.13666e-4_rkx/
  data pd(51:100,2) / &
           -5.40373e-3_rkx ,  2.61635e-5_rkx , -3.33012e-3_rkx ,                &
            0.00000e+0_rkx , -3.08101e-3_rkx , -2.42679e-3_rkx , -3.36086e-3_rkx , &
            0.00000e+0_rkx , -1.18979e+3_rkx , -5.04738e-2_rkx , -2.61547e-3_rkx , &
           -1.03132e-3_rkx ,  1.91583e-4_rkx , -8.38132e+1_rkx , -1.40517e-2_rkx , &
           -1.14167e-2_rkx , -4.08012e-3_rkx ,  1.73522e-4_rkx , -1.39644e-2_rkx , &
           -6.64128e-2_rkx , -6.85152e-2_rkx , -1.34414e+4_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  6.07916e+2_rkx , -4.12220e-3_rkx , &
           -2.20996e-3_rkx ,  0.00000e+0_rkx ,  1.70277e+3_rkx , -4.63015e-3_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx , -2.25360e-3_rkx , -2.96204e-3_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  3.92786e-2_rkx ,  1.31186e-2_rkx , -1.78086e-3_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx , -3.90083e-1_rkx , -2.84741e-2_rkx , &
           -7.78400e+1_rkx , -1.02601e-3_rkx ,  0.00000e+0_rkx/
  data pd(101:150,2) / &
           -7.26485e-4_rkx , -5.42181e-3_rkx , -5.59305e-3_rkx ,                &
            1.22825e-2_rkx ,  1.23868e-2_rkx ,  6.68835e-3_rkx , -1.03303e-2_rkx , &
           -9.51903e-3_rkx ,  2.70021e-4_rkx , -2.57084e-2_rkx , -1.32430e-2_rkx , &
            0.00000e+0_rkx , -3.81000e-2_rkx , -3.16810e-3_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx , -9.05762e-4_rkx , -2.14590e-3_rkx , &
           -1.17824e-3_rkx ,  3.66732e+0_rkx , -3.79729e-4_rkx , -6.13966e-3_rkx , &
           -5.09082e-3_rkx , -1.96332e-3_rkx , -3.08280e-3_rkx , -9.75222e-4_rkx , &
            4.03315e+0_rkx , -2.52710e-1_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx/
  ! N2 DENSITY
  data pd(1:50,3) / &
            1.16112e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  3.33725e-2_rkx , &
            0.00000e+0_rkx ,  3.48637e-2_rkx , -5.44368e-3_rkx ,  0.00000e+0_rkx , &
           -6.73940e-2_rkx ,  1.74754e-1_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  1.74712e+2_rkx ,  0.00000e+0_rkx ,  1.26733e-1_rkx , &
            0.00000e+0_rkx ,  1.03154e+2_rkx ,  5.52075e-2_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  8.13525e-4_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            8.66784e-2_rkx ,  1.58727e-1_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , -2.50482e+1_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , -2.48894e-3_rkx , &
            6.16053e-4_rkx , -5.79716e-4_rkx ,  2.95482e-3_rkx ,  8.47001e-2_rkx , &
            1.70147e-1_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx/
  data pd(51:100,3) / &
           0.00000e+0_rkx ,  2.47425e-5_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx/
  data pd(101:150,3) / &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx/
  ! TLB
  data pd(1:50,4) / &
           9.44846e-1_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , -3.08617e-2_rkx , &
           0.00000e+0_rkx , -2.44019e-2_rkx ,  6.48607e-3_rkx ,  0.00000e+0_rkx , &
           3.08181e-2_rkx ,  4.59392e-2_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  1.74712e+2_rkx ,  0.00000e+0_rkx ,  2.13260e-2_rkx , &
           0.00000e+0_rkx , -3.56958e+2_rkx ,  0.00000e+0_rkx ,  1.82278e-4_rkx , &
           0.00000e+0_rkx ,  3.07472e-4_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           8.66784e-2_rkx ,  1.58727e-1_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           3.83054e-3_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , -1.93065e-3_rkx , &
          -1.45090e-3_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx , -1.23493e-3_rkx ,  1.36736e-3_rkx ,  8.47001e-2_rkx , &
           1.70147e-1_rkx ,  3.71469e-3_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx/
  data pd(51:100,4) / &
           5.10250e-3_rkx ,  2.47425e-5_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  3.68756e-3_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx/
  data pd(101:150,4) / &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx/
  ! O2 DENSITY
  data pd(1:50,5) / &
           1.35580e+0_rkx ,  1.44816e-1_rkx ,  0.00000e+0_rkx ,  6.07767e-2_rkx , &
           0.00000e+0_rkx ,  2.94777e-2_rkx ,  7.46900e-2_rkx ,  0.00000e+0_rkx , &
          -9.23822e-2_rkx ,  8.57342e-2_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  2.38636e+1_rkx ,  0.00000e+0_rkx ,  7.71653e-2_rkx , &
           0.00000e+0_rkx ,  8.18751e+1_rkx ,  1.87736e-2_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  1.49667e-2_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           8.66784e-2_rkx ,  1.58727e-1_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , -3.67874e+2_rkx , &
           5.48158e-3_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  8.47001e-2_rkx , &
           1.70147e-1_rkx ,  1.22631e-2_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx/
  data pd(51:100,5) / &
           8.17187e-3_rkx ,  3.71617e-5_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx , -2.10826e-3_rkx , -3.13640e-3_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
          -7.35742e-2_rkx , -5.00266e-2_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  1.94965e-2_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx/
  data pd(101:150,5) / &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx/
  ! AR DENSITY
  data pd(1:50,6) / &
           1.04761e+0_rkx ,  2.00165e-1_rkx ,  2.37697e-1_rkx ,  3.68552e-2_rkx , &
           0.00000e+0_rkx ,  3.57202e-2_rkx , -2.14075e-1_rkx ,  0.00000e+0_rkx , &
          -1.08018e-1_rkx , -3.73981e-1_rkx ,  0.00000e+0_rkx ,  3.10022e-2_rkx , &
          -1.16305e-3_rkx , -2.07596e+1_rkx ,  0.00000e+0_rkx ,  8.64502e-2_rkx , &
           0.00000e+0_rkx ,  9.74908e+1_rkx ,  5.16707e-2_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           8.66784e-2_rkx ,  1.58727e-1_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  3.46193e+2_rkx , &
           1.34297e-2_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , -3.48509e-3_rkx , &
          -1.54689e-4_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  8.47001e-2_rkx , &
           1.70147e-1_rkx ,  1.47753e-2_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx/
  data pd(51:100,6) / &
           1.89320e-2_rkx ,  3.68181e-5_rkx ,  1.32570e-2_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  3.59719e-3_rkx ,  7.44328e-3_rkx , -1.00023e-3_rkx , &
          -6.50528e+3_rkx ,  0.00000e+0_rkx ,  1.03485e-2_rkx , -1.00983e-3_rkx , &
          -4.06916e-3_rkx , -6.60864e+1_rkx , -1.71533e-2_rkx ,  1.10605e-2_rkx , &
           1.20300e-2_rkx , -5.20034e-3_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx , -2.62769e+3_rkx ,  7.13755e-3_rkx ,  4.17999e-3_rkx , &
           0.00000e+0_rkx ,  1.25910e+4_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx , -2.23595e-3_rkx ,  4.60217e-3_rkx ,  5.71794e-3_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
          -3.18353e-2_rkx , -2.35526e-2_rkx , -1.36189e-2_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  2.03522e-2_rkx , -6.67837e+1_rkx , &
          -1.09724e-3_rkx ,  0.00000e+0_rkx/
  data pd(101:150,6) / &
           -1.38821e-2_rkx ,  1.60468e-2_rkx ,  0.00000e+0_rkx ,                &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  1.51574e-2_rkx , -5.44470e-4_rkx , &
            0.00000e+0_rkx ,  7.28224e-2_rkx ,  6.59413e-2_rkx ,  0.00000e+0_rkx , &
           -5.15692e-3_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , -3.70367e+3_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  1.36131e-2_rkx ,  5.38153e-3_rkx , &
            0.00000e+0_rkx ,  4.76285e+0_rkx , -1.75677e-2_rkx ,  2.26301e-2_rkx , &
            0.00000e+0_rkx ,  1.76631e-2_rkx ,  4.77162e-3_rkx ,  0.00000e+0_rkx , &
            5.39354e+0_rkx ,  0.00000e+0_rkx , -7.51710e-3_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx , -8.82736e+1_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx/
  ! H DENSITY
  data pd(1:50,7) / &
           1.26376e+0_rkx , -2.14304e-1_rkx , -1.49984e-1_rkx ,  2.30404e-1_rkx , &
           2.98237e-2_rkx ,  2.68673e-2_rkx ,  2.96228e-1_rkx ,  2.21900e-2_rkx , &
          -2.07655e-2_rkx ,  4.52506e-1_rkx ,  1.20105e-1_rkx ,  3.24420e-2_rkx , &
           4.24816e-2_rkx , -9.14313e+0_rkx ,  0.00000e+0_rkx ,  2.47178e-2_rkx , &
          -2.88229e-2_rkx ,  8.12805e+1_rkx ,  5.10380e-2_rkx , -5.80611e-3_rkx , &
           2.51236e-5_rkx , -1.24083e-2_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           8.66784e-2_rkx ,  1.58727e-1_rkx , -3.48190e-2_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  2.89885e-5_rkx ,  0.00000e+0_rkx ,  1.53595e+2_rkx , &
          -1.68604e-2_rkx ,  0.00000e+0_rkx ,  1.01015e-2_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  2.84552e-4_rkx , &
          -1.22181e-3_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  8.47001e-2_rkx , &
           1.70147e-1_rkx , -1.04927e-2_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx , -5.91313e-3_rkx/
  data pd(51:100,7) / &
           -2.30501e-2_rkx ,  3.14758e-5_rkx ,  0.00000e+0_rkx ,                &
            0.00000e+0_rkx ,  1.26956e-2_rkx ,  8.35489e-3_rkx ,  3.10513e-4_rkx , &
            0.00000e+0_rkx ,  3.42119e+3_rkx , -2.45017e-3_rkx , -4.27154e-4_rkx , &
            5.45152e-4_rkx ,  1.89896e-3_rkx ,  2.89121e+1_rkx , -6.49973e-3_rkx , &
           -1.93855e-2_rkx , -1.48492e-2_rkx ,  0.00000e+0_rkx , -5.10576e-2_rkx , &
            7.87306e-2_rkx ,  9.51981e-2_rkx , -1.49422e+4_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  2.65503e+2_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  6.37110e-3_rkx ,  3.24789e-4_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  6.14274e-2_rkx ,  1.00376e-2_rkx , -8.41083e-4_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , -1.27099e-2_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx/
  data pd(101:150,7) / &
           -3.94077e-3_rkx , -1.28601e-2_rkx , -7.97616e-3_rkx ,                &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx , -6.71465e-3_rkx , -1.69799e-3_rkx , &
            1.93772e-3_rkx ,  3.81140e+0_rkx , -7.79290e-3_rkx , -1.82589e-2_rkx , &
           -1.25860e-2_rkx , -1.04311e-2_rkx , -3.02465e-3_rkx ,  2.43063e-3_rkx , &
            3.63237e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx/
  ! N DENSITY
  data pd(1:50,8) / &
           7.09557e+1_rkx , -3.26740e-1_rkx ,  0.00000e+0_rkx , -5.16829e-1_rkx , &
          -1.71664e-3_rkx ,  9.09310e-2_rkx , -6.71500e-1_rkx , -1.47771e-1_rkx , &
          -9.27471e-2_rkx , -2.30862e-1_rkx , -1.56410e-1_rkx ,  1.34455e-2_rkx , &
          -1.19717e-1_rkx ,  2.52151e+0_rkx ,  0.00000e+0_rkx , -2.41582e-1_rkx , &
           5.92939e-2_rkx ,  4.39756e+0_rkx ,  9.15280e-2_rkx ,  4.41292e-3_rkx , &
           0.00000e+0_rkx ,  8.66807e-3_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           8.66784e-2_rkx ,  1.58727e-1_rkx ,  9.74701e-2_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  6.70217e+1_rkx , &
          -1.31660e-3_rkx ,  0.00000e+0_rkx , -1.65317e-2_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  8.50247e-2_rkx ,  2.77428e+1_rkx ,  4.98658e-3_rkx , &
           6.15115e-3_rkx ,  9.50156e-3_rkx , -2.12723e-2_rkx ,  8.47001e-2_rkx , &
           1.70147e-1_rkx , -2.38645e-2_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  1.37380e-3_rkx/
  data pd(51:100,8) / &
           -8.41918e-3_rkx ,  2.80145e-5_rkx ,  7.12383e-3_rkx ,                &
            0.00000e+0_rkx , -1.66209e-2_rkx ,  1.03533e-4_rkx , -1.68898e-2_rkx , &
            0.00000e+0_rkx ,  3.64526e+3_rkx ,  0.00000e+0_rkx ,  6.54077e-3_rkx , &
            3.69130e-4_rkx ,  9.94419e-4_rkx ,  8.42803e+1_rkx , -1.16124e-2_rkx , &
           -7.74414e-3_rkx , -1.68844e-3_rkx ,  1.42809e-3_rkx , -1.92955e-3_rkx , &
            1.17225e-1_rkx , -2.41512e-2_rkx ,  1.50521e+4_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  1.60261e+3_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx , -3.54403e-4_rkx , -1.87270e-2_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  2.76439e-2_rkx ,  6.43207e-3_rkx , -3.54300e-2_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , -2.80221e-2_rkx , &
            8.11228e+1_rkx , -6.75255e-4_rkx ,  0.00000e+0_rkx/
  data pd(101:150,8) / &
           -1.05162e-2_rkx , -3.48292e-3_rkx , -6.97321e-3_rkx ,                &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx , -1.45546e-3_rkx , -1.31970e-2_rkx , &
           -3.57751e-3_rkx , -1.09021e+0_rkx , -1.50181e-2_rkx , -7.12841e-3_rkx , &
           -6.64590e-3_rkx , -3.52610e-3_rkx , -1.87773e-2_rkx , -2.22432e-3_rkx , &
           -3.93895e-1_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx/
  ! HOT O DENSITY
  data pd(1:50,9) / &
           6.04050e-2_rkx ,  1.57034e+0_rkx ,  2.99387e-2_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx , -1.51018e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx , -8.61650e+0_rkx ,  1.26454e-2_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  5.50878e-3_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           8.66784e-2_rkx ,  1.58727e-1_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           6.23881e-2_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  8.47001e-2_rkx , &
           1.70147e-1_rkx , -9.45934e-2_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx/
  data pd(51:100,9) / &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx/
  data pd(101:150,9) / &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx/
  ! S PARAM
  data ps(1:50) / &
           9.56827e-1_rkx ,  6.20637e-2_rkx ,  3.18433e-2_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  3.94900e-2_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
          -9.24882e-3_rkx , -7.94023e-3_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  1.74712e+2_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  2.74677e-3_rkx ,  0.00000e+0_rkx ,  1.54951e-2_rkx , &
           8.66784e-2_rkx ,  1.58727e-1_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx , -6.99007e-4_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  1.24362e-2_rkx , -5.28756e-3_rkx ,  8.47001e-2_rkx , &
           1.70147e-1_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx/
  data ps(51:100) / &
           0.00000e+0_rkx ,  2.47425e-5_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx/
  data ps(101:150) / &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx/
  ! TURBO
  data pdl/1.09930e+0_rkx ,  3.90631e+0_rkx ,  3.07165e+0_rkx ,  9.86161e-1_rkx , &
           1.63536e+1_rkx ,  4.63830e+0_rkx ,  1.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  1.28840e+0_rkx ,  3.10302e-2_rkx , &
           1.18339e-1_rkx ,  1.00000e+0_rkx ,  7.00000e-1_rkx ,  1.15020e+0_rkx , &
           3.44689e+0_rkx ,  1.28840e+0_rkx ,  1.00000e+0_rkx ,  1.08738e+0_rkx , &
           1.22947e+0_rkx ,  1.10016e+0_rkx ,  7.34129e-1_rkx ,  1.15241e+0_rkx , &
           2.22784e+0_rkx ,  7.95046e-1_rkx ,  4.01612e+0_rkx ,  4.47749e+0_rkx , &
           1.23435e+2_rkx , -7.60535e-2_rkx ,  1.68986e-6_rkx ,  7.44294e-1_rkx , &
           1.03604e+0_rkx ,  1.72783e+2_rkx ,  1.15020e+0_rkx ,  3.44689e+0_rkx , &
          -7.46230e-1_rkx ,  9.49154e-1_rkx/

  ! LOWER BOUNDARY
  data ptm/1.04130e+3_rkx ,  3.86000e+2_rkx ,  1.95000e+2_rkx ,  1.66728e+1_rkx , &
           2.13000e+2_rkx ,  1.20000e+2_rkx ,  2.40000e+2_rkx ,  1.87000e+2_rkx , &
          -2.00000e+0_rkx ,  0.00000e+0_rkx/
!
  data pdm/2.45600e+7_rkx ,  6.71072e-6_rkx ,  1.00000e+2_rkx ,  0.00000e+0_rkx , &
           1.10000e+2_rkx ,  1.00000e+1_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  8.5940e+10_rkx ,  1.00000e+0_rkx , &
           1.05000e+2_rkx , -8.00000e+0_rkx ,  1.10000e+2_rkx ,  1.00000e+1_rkx , &
           9.00000e+1_rkx ,  2.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           2.8100e+11_rkx ,  0.00000e+0_rkx ,  1.05000e+2_rkx ,  2.80000e+1_rkx , &
           2.89500e+1_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  3.3000e+10_rkx ,  2.68270e-1_rkx , &
           1.05000e+2_rkx ,  1.00000e+0_rkx ,  1.10000e+2_rkx ,  1.00000e+1_rkx , &
           1.10000e+2_rkx , -1.00000e+1_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           1.33000e+9_rkx ,  1.19615e-2_rkx ,  1.05000e+2_rkx ,  0.00000e+0_rkx , &
           1.10000e+2_rkx ,  1.00000e+1_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  1.76100e+5_rkx ,  1.00000e+0_rkx , &
           9.50000e+1_rkx , -8.00000e+0_rkx ,  1.10000e+2_rkx ,  1.00000e+1_rkx , &
           9.00000e+1_rkx ,  2.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           1.00000e+7_rkx ,  1.00000e+0_rkx ,  1.05000e+2_rkx , -8.00000e+0_rkx , &
           1.10000e+2_rkx ,  1.00000e+1_rkx ,  9.00000e+1_rkx ,  2.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  1.00000e+6_rkx ,  1.00000e+0_rkx , &
           1.05000e+2_rkx , -8.00000e+0_rkx ,  5.50000e+2_rkx ,  7.60000e+1_rkx , &
           9.00000e+1_rkx ,  2.00000e+0_rkx ,  0.00000e+0_rkx ,  4.00000e+3_rkx/
  ! TN1(2)
  data ptl(1:50,1) / &
           1.00858e+0_rkx ,  4.56011e-2_rkx , -2.22972e-2_rkx , -5.44388e-2_rkx , &
           5.23136e-4_rkx , -1.88849e-2_rkx ,  5.23707e-2_rkx , -9.43646e-3_rkx , &
           6.31707e-3_rkx , -7.80460e-2_rkx , -4.88430e-2_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx , -7.60250e+0_rkx ,  0.00000e+0_rkx , -1.44635e-2_rkx , &
          -1.76843e-2_rkx , -1.21517e+2_rkx ,  2.85647e-2_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  6.31792e-4_rkx ,  0.00000e+0_rkx ,  5.77197e-3_rkx , &
           8.66784e-2_rkx ,  1.58727e-1_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , -8.90272e+3_rkx , &
           3.30611e-3_rkx ,  3.02172e-3_rkx ,  0.00000e+0_rkx , -2.13673e-3_rkx , &
          -3.20910e-4_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  2.76034e-3_rkx , &
           2.82487e-3_rkx , -2.97592e-4_rkx , -4.21534e-3_rkx ,  8.47001e-2_rkx , &
           1.70147e-1_rkx ,  8.96456e-3_rkx ,  0.00000e+0_rkx , -1.08596e-2_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx/
  data ptl(51:100,1) / &
           5.57917e-3_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  9.65405e-3_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  2.00000e+0_rkx/
  ! TN1(3)
  data ptl(1:50,2) / &
           9.39664e-1_rkx ,  8.56514e-2_rkx , -6.79989e-3_rkx ,  2.65929e-2_rkx , &
          -4.74283e-3_rkx ,  1.21855e-2_rkx , -2.14905e-2_rkx ,  6.49651e-3_rkx , &
          -2.05477e-2_rkx , -4.24952e-2_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  1.19148e+1_rkx ,  0.00000e+0_rkx ,  1.18777e-2_rkx , &
          -7.28230e-2_rkx , -8.15965e+1_rkx ,  1.73887e-2_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx , -1.44691e-2_rkx ,  2.80259e-4_rkx , &
           8.66784e-2_rkx ,  1.58727e-1_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  2.16584e+2_rkx , &
           3.18713e-3_rkx ,  7.37479e-3_rkx ,  0.00000e+0_rkx , -2.55018e-3_rkx , &
          -3.92806e-3_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , -2.89757e-3_rkx , &
          -1.33549e-3_rkx ,  1.02661e-3_rkx ,  3.53775e-4_rkx ,  8.47001e-2_rkx , &
           1.70147e-1_rkx , -9.17497e-3_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx/
  data ptl(51:100,2) / &
           3.56082e-3_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx , -1.00902e-2_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  2.00000e+0_rkx/
  ! TN1(4)
  data ptl(1:50,3) / &
           9.85982e-1_rkx , -4.55435e-2_rkx ,  1.21106e-2_rkx ,  2.04127e-2_rkx , &
          -2.40836e-3_rkx ,  1.11383e-2_rkx , -4.51926e-2_rkx ,  1.35074e-2_rkx , &
          -6.54139e-3_rkx ,  1.15275e-1_rkx ,  1.28247e-1_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx , -5.30705e+0_rkx ,  0.00000e+0_rkx , -3.79332e-2_rkx , &
          -6.24741e-2_rkx ,  7.71062e-1_rkx ,  2.96315e-2_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  6.81051e-3_rkx , -4.34767e-3_rkx , &
           8.66784e-2_rkx ,  1.58727e-1_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  1.07003e+1_rkx , &
          -2.76907e-3_rkx ,  4.32474e-4_rkx ,  0.00000e+0_rkx ,  1.31497e-3_rkx , &
          -6.47517e-4_rkx ,  0.00000e+0_rkx , -2.20621e+1_rkx , -1.10804e-3_rkx , &
          -8.09338e-4_rkx ,  4.18184e-4_rkx ,  4.29650e-3_rkx ,  8.47001e-2_rkx , &
           1.70147e-1_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx/
  data ptl(51:100,3) / &
           -4.04337e-3_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,                &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , -9.52550e-4_rkx , &
            8.56253e-4_rkx ,  4.33114e-4_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  1.21223e-3_rkx ,  2.38694e-4_rkx ,  9.15245e-4_rkx , &
            1.28385e-3_rkx ,  8.67668e-4_rkx , -5.61425e-6_rkx ,  1.04445e+0_rkx , &
            3.41112e+1_rkx ,  0.00000e+0_rkx , -8.40704e-1_rkx , -2.39639e+2_rkx , &
            7.06668e-1_rkx , -2.05873e+1_rkx , -3.63696e-1_rkx ,  2.39245e+1_rkx , &
            0.00000e+0_rkx , -1.06657e-3_rkx , -7.67292e-4_rkx ,  1.54534e-4_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  2.00000e+0_rkx/
  ! TN1(5) TN2(1)
  data ptl(1:50,4) / &
           1.00320e+0_rkx ,  3.83501e-2_rkx , -2.38983e-3_rkx ,  2.83950e-3_rkx , &
           4.20956e-3_rkx ,  5.86619e-4_rkx ,  2.19054e-2_rkx , -1.00946e-2_rkx , &
          -3.50259e-3_rkx ,  4.17392e-2_rkx , -8.44404e-3_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  4.96949e+0_rkx ,  0.00000e+0_rkx , -7.06478e-3_rkx , &
          -1.46494e-2_rkx ,  3.13258e+1_rkx , -1.86493e-3_rkx ,  0.00000e+0_rkx , &
          -1.67499e-2_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  5.12686e-4_rkx , &
           8.66784e-2_rkx ,  1.58727e-1_rkx , -4.64167e-3_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  4.37353e-3_rkx , -1.99069e+2_rkx , &
           0.00000e+0_rkx , -5.34884e-3_rkx ,  0.00000e+0_rkx ,  1.62458e-3_rkx , &
           2.93016e-3_rkx ,  2.67926e-3_rkx ,  5.90449e+2_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx , -1.17266e-3_rkx , -3.58890e-4_rkx ,  8.47001e-2_rkx , &
           1.70147e-1_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  1.38673e-2_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx/
  data ptl(51:100,4) / &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  1.60571e-3_rkx ,  6.28078e-4_rkx , &
           5.05469e-5_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
          -1.57829e-3_rkx , -4.00855e-4_rkx ,  5.04077e-5_rkx , -1.39001e-3_rkx , &
          -2.33406e-3_rkx , -4.81197e-4_rkx ,  1.46758e+0_rkx ,  6.20332e+0_rkx , &
           0.00000e+0_rkx ,  3.66476e-1_rkx , -6.19760e+1_rkx ,  3.09198e-1_rkx , &
          -1.98999e+1_rkx ,  0.00000e+0_rkx , -3.29933e+2_rkx ,  0.00000e+0_rkx , &
          -1.10080e-3_rkx , -9.39310e-5_rkx ,  1.39638e-4_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  2.00000e+0_rkx/
  ! TN2(2)
  data pma(1:50,1) / &
           9.81637e-1_rkx , -1.41317e-3_rkx ,  3.87323e-2_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx , -3.58707e-2_rkx , -8.63658e-3_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx , -2.02226e+0_rkx ,  0.00000e+0_rkx , -8.69424e-3_rkx , &
          -1.91397e-2_rkx ,  8.76779e+1_rkx ,  4.52188e-3_rkx ,  0.00000e+0_rkx , &
           2.23760e-2_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx , -7.07572e-3_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx , -4.11210e-3_rkx ,  3.50060e+1_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx , -8.36657e-3_rkx ,  1.61347e+1_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , -1.45130e-2_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx/
  data pma(51:100,1) / &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  1.24152e-3_rkx ,  6.43365e-4_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           1.33255e-3_rkx ,  2.42657e-3_rkx ,  1.60666e-3_rkx , -1.85728e-3_rkx , &
          -1.46874e-3_rkx , -4.79163e-6_rkx ,  1.22464e+0_rkx ,  3.53510e+1_rkx , &
           0.00000e+0_rkx ,  4.49223e-1_rkx , -4.77466e+1_rkx ,  4.70681e-1_rkx , &
           8.41861e+0_rkx , -2.88198e-1_rkx ,  1.67854e+2_rkx ,  0.00000e+0_rkx , &
           7.11493e-4_rkx ,  6.05601e-4_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  2.00000e+0_rkx/
  ! TN2(3)
  data pma(1:50,2) / &
           1.00422e+0_rkx , -7.11212e-3_rkx ,  5.24480e-3_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx , -5.28914e-2_rkx , -2.41301e-2_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx , -2.12219e+1_rkx , -1.03830e-2_rkx , -3.28077e-3_rkx , &
           1.65727e-2_rkx ,  1.68564e+0_rkx , -6.68154e-3_rkx ,  0.00000e+0_rkx , &
           1.45155e-2_rkx ,  0.00000e+0_rkx ,  8.42365e-3_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx , -4.34645e-3_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  2.16780e-2_rkx ,  0.00000e+0_rkx , -1.38459e+2_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  7.04573e-3_rkx , -4.73204e+1_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  1.08767e-2_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx/
  data pma(51:100,2) / &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx , -8.08279e-3_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  5.21769e-4_rkx , -2.27387e-4_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           3.26769e-3_rkx ,  3.16901e-3_rkx ,  4.60316e-4_rkx , -1.01431e-4_rkx , &
           1.02131e-3_rkx ,  9.96601e-4_rkx ,  1.25707e+0_rkx ,  2.50114e+1_rkx , &
           0.00000e+0_rkx ,  4.24472e-1_rkx , -2.77655e+1_rkx ,  3.44625e-1_rkx , &
           2.75412e+1_rkx ,  0.00000e+0_rkx ,  7.94251e+2_rkx ,  0.00000e+0_rkx , &
           2.45835e-3_rkx ,  1.38871e-3_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  2.00000e+0_rkx/
  ! TN2(4) TN3(1)
  data pma(1:50,3) / &
           1.01890e+0_rkx , -2.46603e-2_rkx ,  1.00078e-2_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx , -6.70977e-2_rkx , -4.02286e-2_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx , -2.29466e+1_rkx , -7.47019e-3_rkx ,  2.26580e-3_rkx , &
           2.63931e-2_rkx ,  3.72625e+1_rkx , -6.39041e-3_rkx ,  0.00000e+0_rkx , &
           9.58383e-3_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx , -1.85291e-3_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  1.39717e+2_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  9.19771e-3_rkx , -3.69121e+2_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , -1.57067e-2_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx/
  data pma(51:100,3) / &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx , -7.07265e-3_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx , -2.92953e-3_rkx , -2.77739e-3_rkx , &
          -4.40092e-4_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           2.47280e-3_rkx ,  2.95035e-4_rkx , -1.81246e-3_rkx ,  2.81945e-3_rkx , &
           4.27296e-3_rkx ,  9.78863e-4_rkx ,  1.40545e+0_rkx , -6.19173e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx , -7.93632e+1_rkx ,  4.44643e-1_rkx , &
          -4.03085e+2_rkx ,  0.00000e+0_rkx ,  1.15603e+1_rkx ,  0.00000e+0_rkx , &
           2.25068e-3_rkx ,  8.48557e-4_rkx , -2.98493e-4_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  2.00000e+0_rkx/
  ! TN3(2)
  data pma(1:50,4) / &
           9.75801e-1_rkx ,  3.80680e-2_rkx , -3.05198e-2_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  3.85575e-2_rkx ,  5.04057e-2_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx , -1.76046e+2_rkx ,  1.44594e-2_rkx , -1.48297e-3_rkx , &
          -3.68560e-3_rkx ,  3.02185e+1_rkx , -3.23338e-3_rkx ,  0.00000e+0_rkx , &
           1.53569e-2_rkx ,  0.00000e+0_rkx , -1.15558e-2_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  4.89620e-3_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx , -1.00616e-2_rkx , -8.21324e-3_rkx , -1.57757e+2_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  6.63564e-3_rkx ,  4.58410e+1_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , -2.51280e-2_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx/
  data pma(51:100,4) / &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  9.91215e-3_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx , -8.73148e-4_rkx , -1.29648e-3_rkx , &
          -7.32026e-5_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
          -4.68110e-3_rkx , -4.66003e-3_rkx , -1.31567e-3_rkx , -7.39390e-4_rkx , &
           6.32499e-4_rkx , -4.65588e-4_rkx , -1.29785e+0_rkx , -1.57139e+2_rkx , &
           0.00000e+0_rkx ,  2.58350e-1_rkx , -3.69453e+1_rkx ,  4.10672e-1_rkx , &
           9.78196e+0_rkx , -1.52064e-1_rkx , -3.85084e+3_rkx ,  0.00000e+0_rkx , &
          -8.52706e-4_rkx , -1.40945e-3_rkx , -7.26786e-4_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  2.00000e+0_rkx/
  ! TN3(3)
  data pma(1:50,5) / &
           9.60722e-1_rkx ,  7.03757e-2_rkx , -3.00266e-2_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  2.22671e-2_rkx ,  4.10423e-2_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx , -1.63070e+2_rkx ,  1.06073e-2_rkx ,  5.40747e-4_rkx , &
           7.79481e-3_rkx ,  1.44908e+2_rkx ,  1.51484e-4_rkx ,  0.00000e+0_rkx , &
           1.97547e-2_rkx ,  0.00000e+0_rkx , -1.41844e-2_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  5.77884e-3_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  9.74319e-3_rkx ,  0.00000e+0_rkx , -2.88015e+3_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx , -4.44902e-3_rkx , -2.92760e+1_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  2.34419e-2_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx/
  data pma(51:100,5) / &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  5.36685e-3_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx , -4.65325e-4_rkx , -5.50628e-4_rkx , &
           3.31465e-4_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
          -2.06179e-3_rkx , -3.08575e-3_rkx , -7.93589e-4_rkx , -1.08629e-4_rkx , &
           5.95511e-4_rkx , -9.05050e-4_rkx ,  1.18997e+0_rkx ,  4.15924e+1_rkx , &
           0.00000e+0_rkx , -4.72064e-1_rkx , -9.47150e+2_rkx ,  3.98723e-1_rkx , &
           1.98304e+1_rkx ,  0.00000e+0_rkx ,  3.73219e+3_rkx ,  0.00000e+0_rkx , &
          -1.50040e-3_rkx , -1.14933e-3_rkx , -1.56769e-4_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  2.00000e+0_rkx/
  ! TN3(4)
  data pma(1:50,6) / &
           1.03123e+0_rkx , -7.05124e-2_rkx ,  8.71615e-3_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx , -3.82621e-2_rkx , -9.80975e-3_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  2.89286e+1_rkx ,  9.57341e-3_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  8.66153e+1_rkx ,  7.91938e-4_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  4.68917e-3_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  7.86638e-3_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  9.90827e-3_rkx ,  0.00000e+0_rkx ,  6.55573e+1_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx , -4.00200e+1_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  7.07457e-3_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx/
  data pma(51:100,6) / &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  5.72268e-3_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx , -2.04970e-4_rkx ,  1.21560e-3_rkx , &
          -8.05579e-6_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
          -2.49941e-3_rkx , -4.57256e-4_rkx , -1.59311e-4_rkx ,  2.96481e-4_rkx , &
          -1.77318e-3_rkx , -6.37918e-4_rkx ,  1.02395e+0_rkx ,  1.28172e+1_rkx , &
           0.00000e+0_rkx ,  1.49903e-1_rkx , -2.63818e+1_rkx ,  0.00000e+0_rkx , &
           4.70628e+1_rkx , -2.22139e-1_rkx ,  4.82292e-2_rkx ,  0.00000e+0_rkx , &
          -8.67075e-4_rkx , -5.86479e-4_rkx ,  5.32462e-4_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  2.00000e+0_rkx/
  ! TN3(5) SURFACE TEMP TSL
  data pma(1:50,7) / &
           1.00828e+0_rkx , -9.10404e-2_rkx , -2.26549e-2_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx , -2.32420e-2_rkx , -9.08925e-3_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  3.36105e+1_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx , -1.24957e+1_rkx , -5.87939e-3_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  2.79765e+1_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  2.01237e+3_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , -1.75553e-2_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx/
  data pma(51:100,7) / &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  3.29699e-3_rkx ,  1.26659e-3_rkx , &
           2.68402e-4_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           1.17894e-3_rkx ,  1.48746e-3_rkx ,  1.06478e-4_rkx ,  1.34743e-4_rkx , &
          -2.20939e-3_rkx , -6.23523e-4_rkx ,  6.36539e-1_rkx ,  1.13621e+1_rkx , &
           0.00000e+0_rkx , -3.93777e-1_rkx ,  2.38687e+3_rkx ,  0.00000e+0_rkx , &
           6.61865e+2_rkx , -1.21434e-1_rkx ,  9.27608e+0_rkx ,  0.00000e+0_rkx , &
           1.68478e-4_rkx ,  1.24892e-3_rkx ,  1.71345e-3_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  2.00000e+0_rkx/
  ! TGN3(2) SURFACE GRAD TSLG
  data pma(1:50,8) / &
           1.57293e+0_rkx , -6.78400e-1_rkx ,  6.47500e-1_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx , -7.62974e-2_rkx , -3.60423e-1_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  1.28358e+2_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  4.68038e+1_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx , -1.67898e-1_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  2.90994e+4_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  3.15706e+1_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           0.00000e+0_rkx ,  0.00000e+0_rkx/
  data pma(51:100,8) / &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  &
           0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  &
           0.00000e+0_rkx ,  2.00000e+0_rkx/
  ! TGN2(1) TGN1(2)
  data pma(1:50,9) / &
            8.60028e-1_rkx ,  3.77052e-1_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx , -1.17570e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  7.77757e-3_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  1.01024e+2_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  6.54251e+2_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx/
  data pma(51:100,9) / &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
           -1.56959e-2_rkx ,  1.91001e-2_rkx ,  3.15971e-2_rkx ,  1.00982e-2_rkx , &
           -6.71565e-3_rkx ,  2.57693e-3_rkx ,  1.38692e+0_rkx ,  2.82132e-1_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  3.81511e+2_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  2.00000e+0_rkx/
  ! TGN3(1) TGN2(2)
  data pma(1:50,10) / &
            1.06029e+0_rkx , -5.25231e-2_rkx ,  3.73034e-1_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  3.31072e-2_rkx , -3.88409e-1_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx , -1.65295e+2_rkx , -2.13801e-1_rkx , -4.38916e-2_rkx , &
           -3.22716e-1_rkx , -8.82393e+1_rkx ,  1.18458e-1_rkx ,  0.00000e+0_rkx , &
           -4.35863e-1_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx , -1.19782e-1_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  2.62229e+1_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx , -5.37443e+1_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , -4.55788e-1_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx/
  data pma(51:100,10) / &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  3.84009e-2_rkx ,  3.96733e-2_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            5.05494e-2_rkx ,  7.39617e-2_rkx ,  1.92200e-2_rkx , -8.46151e-3_rkx , &
           -1.34244e-2_rkx ,  1.96338e-2_rkx ,  1.50421e+0_rkx ,  1.88368e+1_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx , -5.13114e+1_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            5.11923e-2_rkx ,  3.61225e-2_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx ,  0.00000e+0_rkx , &
            0.00000e+0_rkx ,  2.00000e+0_rkx/
  ! SEMIANNUAL MULT SAM
  data sam(1:50) / &
            1.00000e+0_rkx , 1.00000e+0_rkx , 1.00000e+0_rkx , 1.00000e+0_rkx , &
            1.00000e+0_rkx , 1.00000e+0_rkx , 1.00000e+0_rkx , 1.00000e+0_rkx , &
            1.00000e+0_rkx , 1.00000e+0_rkx , 1.00000e+0_rkx , 1.00000e+0_rkx , &
            1.00000e+0_rkx , 1.00000e+0_rkx , 1.00000e+0_rkx , 1.00000e+0_rkx , &
            1.00000e+0_rkx , 1.00000e+0_rkx , 1.00000e+0_rkx , 1.00000e+0_rkx , &
            1.00000e+0_rkx , 1.00000e+0_rkx , 1.00000e+0_rkx , 1.00000e+0_rkx , &
            1.00000e+0_rkx , 1.00000e+0_rkx , 1.00000e+0_rkx , 1.00000e+0_rkx , &
            1.00000e+0_rkx , 1.00000e+0_rkx , 1.00000e+0_rkx , 1.00000e+0_rkx , &
            1.00000e+0_rkx , 1.00000e+0_rkx , 1.00000e+0_rkx , 1.00000e+0_rkx , &
            1.00000e+0_rkx , 1.00000e+0_rkx , 1.00000e+0_rkx , 1.00000e+0_rkx , &
            1.00000e+0_rkx , 1.00000e+0_rkx , 1.00000e+0_rkx , 1.00000e+0_rkx , &
            1.00000e+0_rkx , 1.00000e+0_rkx , 1.00000e+0_rkx , 1.00000e+0_rkx , &
            1.00000e+0_rkx , 1.00000e+0_rkx/
  data sam(51:100) / &
            0.00000e+0_rkx , 0.00000e+0_rkx , 0.00000e+0_rkx , 0.00000e+0_rkx , &
            0.00000e+0_rkx , 0.00000e+0_rkx , 0.00000e+0_rkx , 0.00000e+0_rkx , &
            0.00000e+0_rkx , 0.00000e+0_rkx , 0.00000e+0_rkx , 0.00000e+0_rkx , &
            0.00000e+0_rkx , 0.00000e+0_rkx , 0.00000e+0_rkx , 0.00000e+0_rkx , &
            0.00000e+0_rkx , 0.00000e+0_rkx , 0.00000e+0_rkx , 0.00000e+0_rkx , &
            0.00000e+0_rkx , 0.00000e+0_rkx , 0.00000e+0_rkx , 0.00000e+0_rkx , &
            0.00000e+0_rkx , 0.00000e+0_rkx , 0.00000e+0_rkx , 0.00000e+0_rkx , &
            0.00000e+0_rkx , 0.00000e+0_rkx , 0.00000e+0_rkx , 0.00000e+0_rkx , &
            0.00000e+0_rkx , 0.00000e+0_rkx , 0.00000e+0_rkx , 0.00000e+0_rkx , &
            0.00000e+0_rkx , 0.00000e+0_rkx , 0.00000e+0_rkx , 0.00000e+0_rkx , &
            0.00000e+0_rkx , 0.00000e+0_rkx , 0.00000e+0_rkx , 0.00000e+0_rkx , &
            0.00000e+0_rkx , 0.00000e+0_rkx , 0.00000e+0_rkx , 0.00000e+0_rkx , &
            0.00000e+0_rkx , 0.00000e+0_rkx/
  ! MIDDLE ATMOSPHERE AVERAGES
  data pavgm /2.61000e+2_rkx ,  2.64000e+2_rkx ,  2.29000e+2_rkx ,  2.17000e+2_rkx , &
              2.17000e+2_rkx ,  2.23000e+2_rkx ,  2.86760e+2_rkx , -2.93940e+0_rkx , &
              2.50000e+0_rkx ,  0.00000e+0_rkx/
  contains
!
!   Chemistry/dissociation correction for msis models
!
    real(rkx) function ccor(alt,r,h1,zh)
      implicit none
      ! alt - Altitude
      ! r   - Target ratio
      ! h1  - Transition scale length
      ! zh  - Altitude of 1/2 r
      real(rkx) , intent(in) :: alt , h1 , r , zh
      real(rkx) , parameter :: echeck = 70.0_rkx
      real(rkx) :: e , ex
      e = (alt-zh)/h1
      if ( e > echeck ) then
        ccor = d_zero
      else if ( e < -echeck ) then
        ccor = r
      else
        ex = exp(e)
        ccor = r/(d_one+ex)
      end if
      ccor = exp(ccor)
    end function ccor
!
!   O&O2 Chemistry/dissociation correction for msis models
!
    real(rkx) function ccor2(alt,r,h1,zh,h2)
      implicit none
      ! alt - Altitude
      ! r   - Target ratio
      ! h1  - Transition scale length
      ! zh  - Altitude of 1/2 r
      real(rkx) , intent(in) :: alt , h1 , h2 , r , zh
      real(rkx) :: e1 , e2 , ex1 , ex2
      real(rkx) , parameter :: echeck = 70.0_rkx
!
      e1 = (alt-zh)/h1
      e2 = (alt-zh)/h2
      if ( e1 > echeck .or. e2 > echeck ) then
        ccor2 = d_zero
      else if ( e1 < -echeck .and. e2 < -echeck ) then
        ccor2 = r
      else
        ex1 = exp(e1)
        ex2 = exp(e2)
        ccor2 = r/(d_one+d_half*(ex1+ex2))
      end if
      ccor2 = exp(ccor2)
    end function ccor2
!
! Calculate Temperature and Density Profiles for lower atmos.
!
    real(rkx) function densm(alt,d0,xm,tz)
      implicit none
!
      real(rkx) , intent(in) :: alt , d0 , xm
      real(rkx) , intent(out) :: tz
!
      real(rkx) :: expl , gamm , glb , t1 , t2 , x , y , yd1 , &
                 yd2 , yi , z , z1 , z2 , zg , zgdif
      integer(ik4) :: k , mn
      real(rkx) , dimension(10) :: xs , y2out , ys

      densm = d0
      if ( alt <= zn2(1) ) then
        ! STRATOSPHERE/MESOSPHERE TEMPERATURE
        z = max(alt,zn2(mn2))
        mn = mn2
        z1 = zn2(1)
        z2 = zn2(mn)
        t1 = tn2(1)
        t2 = tn2(mn)
        zg = zeta(z,z1)
        zgdif = zeta(z2,z1)
        ! Set up spline nodes
        do k = 1 , mn
          xs(k) = zeta(zn2(k),z1)/zgdif
          ys(k) = d_one/tn2(k)
        end do
        yd1 = -tgn2(1)/(t1*t1)*zgdif
        yd2 = -tgn2(2)/(t2*t2)*zgdif*((re+z2)/(re+z1))**2
        ! Calculate spline coefficients
        call spline(xs(1:mn),ys(1:mn),yd1,yd2,y2out(1:mn))
        x = zg/zgdif
        call splint(xs(1:mn),ys(1:mn),y2out(1:mn),x,y)
        ! Temperature at altitude
        tz = d_one/y
        if ( abs(xm) > nearzero ) then
          ! CALCULATE STRATOSPHERE/MESOSPHERE DENSITY
          glb = gsurf/(d_one+z1/re)**2
          gamm = xm*glb*zgdif/r100gas
          ! Integrate temperature profile
          call splini(xs(1:mn),ys(1:mn),y2out(1:mn),x,yi)
          expl = gamm*yi
          if ( expl > 50.0_rkx ) expl = 50.0_rkx
          ! Density at altitude
          densm = densm*(t1/tz)*exp(-expl)
        end if
        if ( alt <= zn3(1) ) then
          ! TROPOSPHERE/STRATOSPHERE TEMPERATURE
          z = alt
          mn = mn3
          z1 = zn3(1)
          z2 = zn3(mn)
          t1 = tn3(1)
          t2 = tn3(mn)
          zg = zeta(z,z1)
          zgdif = zeta(z2,z1)
          ! Set up spline nodes
          do k = 1 , mn
            xs(k) = zeta(zn3(k),z1)/zgdif
            ys(k) = d_one/tn3(k)
          end do
          yd1 = -tgn3(1)/(t1*t1)*zgdif
          yd2 = -tgn3(2)/(t2*t2)*zgdif*((re+z2)/(re+z1))**2
          ! Calculate spline coefficients
          call spline(xs(1:mn),ys(1:mn),yd1,yd2,y2out(1:mn))
          x = zg/zgdif
          call splint(xs(1:mn),ys(1:mn),y2out(1:mn),x,y)
          ! temperature at altitude
          tz = d_one/y
          if ( abs(xm) > nearzero ) then
            ! CALCULATE TROPOSPHERIC/STRATOSPHERE DENSITY
            glb = gsurf/(d_one+z1/re)**2
            gamm = xm*glb*zgdif/r100gas
            ! Integrate temperature profile
            call splini(xs(1:mn),ys(1:mn),y2out(1:mn),x,yi)
            expl = gamm*yi
            if ( expl > 50.0_rkx ) expl = 50.0_rkx
            ! Density at altitude
            densm = densm*(t1/tz)*exp(-expl)
          end if
        end if
      end if
      if ( abs(xm) < nearzero ) densm = tz
      contains
        real(rkx) function zeta(zz,zl)
          implicit none
          real(rkx) , intent(in) :: zz , zl
          zeta = (zz-zl)*(re+zl)/(re+zz)
        end function zeta
    end function densm
!
!   Calculate Temperature and Density Profiles for MSIS models
!   New lower thermo polynomial 10/30/89
!
    real(rkx) function densu(alt,dlb,t1,t2,xm,xalph,tz,zlb,s2)
      implicit none
!
      real(rkx) , intent(in) :: xalph , alt , t1 , t2 , dlb , s2 , xm , zlb
      real(rkx) , intent(out) :: tz
!
      real(rkx) :: densa , dta , expl , gamm , gammo , glb , &
                 tt1 , tt2 , ta , tt , x , y , yd1 , yd2 , yi , z ,   &
                 z1 , z2 , za , zg , zg2 , zgdif
      integer(ik4) :: k , mn
      real(rkx) , dimension(5) :: xs , y2out , ys
!
      densu = d_one
!     Joining altitude of Bates and spline
      za = zn1(1)
      z = max(alt,za)
!     Geopotential altitude difference from ZLB
      zg2 = zeta(z,zlb)
!     Bates temperature
      tt = t1 - (t1-t2)*exp(-s2*zg2)
      ta = tt
      tz = tt
      densu = tz
      if ( alt < za ) then
        ! CALCULATE TEMPERATURE BELOW ZA
        ! Temperature gradient at ZA from Bates profile
        dta = (t1-ta)*s2*((re+zlb)/(re+za))**2
        tgn1(1) = dta
        tn1(1) = ta
        z = max(alt,zn1(mn1))
        mn = mn1
        z1 = zn1(1)
        z2 = zn1(mn)
        tt1 = tn1(1)
        tt2 = tn1(mn)
        ! Geopotental difference from Z1
        zg = zeta(z,z1)
        zgdif = zeta(z2,z1)
        ! Set up spline nodes
        do k = 1 , mn
          xs(k) = zeta(zn1(k),z1)/zgdif
          ys(k) = d_one/tn1(k)
        end do
        ! End node derivatives
        yd1 = -tgn1(1)/(tt1*tt1)*zgdif
        yd2 = -tgn1(2)/(tt2*tt2)*zgdif*((re+z2)/(re+z1))**2
        ! Calculate spline coefficients
        call spline(xs,ys,yd1,yd2,y2out)
        x = zg/zgdif
        call splint(xs,ys,y2out,x,y)
        ! temperature at altitude
        tz = d_one/y
        densu = tz
      end if
      if ( abs(xm) > nearzero ) then
        ! CALCULATE DENSITY ABOVE ZA
        glb = gsurf/(d_one+zlb/re)**2
        gammo = xm*glb/(s2*r100gas*t1)
        expl = exp(-s2*gammo*zg2)
        if ( expl > 50.0_rkx .or. tt <= d_zero ) expl = 50.0_rkx
        ! Density at altitude
        densa = dlb*(t2/tt)**(d_one+xalph+gammo)*expl
        densu = densa
        if ( alt<za ) then
          ! CALCULATE DENSITY BELOW ZA
          glb = gsurf/(d_one+z1/re)**2
          gamm = xm*glb*zgdif/r100gas
          ! integrate spline temperatures
          call splini(xs,ys,y2out,x,yi)
          expl = gamm*yi
          if ( expl > 50.0_rkx .or. tz <= d_zero ) expl = 50.0_rkx
          ! Density at altitude
          densu = densu*(tt1/tz)**(d_one+xalph)*exp(-expl)
        end if
      end if
      contains
        real(rkx) function zeta(zz,zl)
          implicit none
          real(rkx) , intent(in) :: zz , zl
          zeta = (zz-zl)*(re+zl)/(re+zz)
        end function zeta
    end function densu
!
!   Turbopause correction for msis models
!
    real(rkx) function dnet(dd,dm,zhm,xmm,xm)
      implicit none
!
!     DD - diffusive density
!     DM - full mixed density
!     ZHM - transition scale length
!     XMM - full mixed molecular weight
!     XM  - species molecular weight
!     DNET - combined density
!
      real(rkx) , intent(in) :: dm , xm , xmm , zhm
      real(rkx) , intent(inout) :: dd
!
      real(rkx) :: a , ylog
!
      a = zhm/(xmm-xm)
      if ( dm <= d_zero .or. dd <= d_zero ) then
        if ( abs(dd) < nearzero .and. abs(dm) < nearzero ) dd = d_one
        if ( abs(dm) < nearzero ) then
          dnet = dd
        else if ( abs(dd) < nearzero ) then
          dnet = dm
        end if
      else
        ylog = a*log(dm/dd)
        if ( ylog < -d_10 ) then
          dnet = dd
        else if ( ylog > d_10 ) then
          dnet = dm
        else
          dnet = dd*(d_one+exp(ylog))**(d_one/a)
        end if
      end if
    end function dnet
!
!   SET SWITCHES
!   SW FOR MAIN TERMS, SWC FOR CROSS TERMS
!
!   TO TURN ON AND OFF PARTICULAR VARIATIONS CALL TSELEC(SV),
!   WHERE SV IS A 25 ELEMENT ARRAY CONTAINING 0. FOR OFF, 1.
!   FOR ON, OR 2. FOR MAIN EFFECTS OFF BUT CROSS TERMS ON
!
!   To get current values of SW: CALL TRETRV(SW)
!
    subroutine tselec(sv)
      implicit none
!
      real(rkx) , dimension(25) , intent(in) :: sv
!
      integer(ik4) :: i

      do i = 1 , 25
        sav(i) = sv(i)
        sw(i) = mod(sv(i),d_two)
        if ( abs(sv(i)-d_one) < nearzero .or. &
             abs(sv(i)-d_two) < nearzero ) then
          swc(i) = d_one
        else
          swc(i) = d_zero
        end if
      end do
      isw = 64999
      return
    end subroutine tselec
!
    subroutine tretrv(svv)
      implicit none
      real(rkx) , dimension(25) , intent(out) :: svv
      integer(ik4) :: i
      do i = 1 , 25
        svv(i) = sav(i)
      end do
    end subroutine tretrv
!
!   CALCULATE LATITUDE VARIABLE GRAVITY (GV) AND EFFECTIVE RADIUS (REFF)
!
    subroutine glatf(lat)
      implicit none
      real(rkx) , intent(in) :: lat
      real(rkx)  :: c2
      c2 = cos(d_two*degrad*lat)
      gsurf = egrav*100.0_rkx*(d_one-0.0026373_rkx*c2)
      re = d_two*gsurf/(3.085462e-6_rkx+2.27e-9_rkx*c2)*1.e-5_rkx
    end subroutine glatf
!
    real(rkx) function glob7s(glong,p)
      implicit none
!
      real(rkx) , intent(in) :: glong
      real(rkx) , dimension(:) , intent(out) :: p
!
      real(rkx) :: t82 , tt , t71 , t72 , t81
      integer(ik4) :: i
      real(rkx) , dimension(14) :: t
!
!     VERSION OF GLOBE FOR LOWER ATMOSPHERE 10/26/99
!
      real(rkx) , parameter :: pset = d_two
      real(rkx) , save :: dayl
      real(rkx) , save :: p32 , p18 , p14 , p39
      real(rkx) , save :: cd14 , cd18 , cd32 , cd39
      data dayl /-1.0_rkx/
      data p32  /-1000.0_rkx/
      data p18  /-1000.0_rkx/
      data p14  /-1000.0_rkx/
      data p39  /-1000.0_rkx/

      t(:) = d_zero

      if ( abs(p(100)) < nearzero ) p(100) = pset

      if ( abs(day-dayl) > nearzero .or. abs(p32-p(32)) > nearzero ) then
        cd32 = cos(dr*(day-p(32)))
      end if
      if ( abs(day-dayl) > nearzero .or. abs(p18-p(18)) > nearzero ) then
        cd18 = cos(d_two*dr*(day-p(18)))
      end if
      if ( abs(day-dayl) > nearzero .or. abs(p14-p(14)) > nearzero ) then
        cd14 = cos(dr*(day-p(14)))
      end if
      if ( abs(day-dayl) > nearzero .or. abs(p39-p(39)) > nearzero ) then
        cd39 = cos(d_two*dr*(day-p(39)))
      end if

      dayl = day
      p32 = p(32)
      p18 = p(18)
      p14 = p(14)
      p39 = p(39)
!
      t(1) = p(22)*dfa
      ! TIME INDEPENDENT
      t(2) = p(2)*plg(3,1) + p(3)*plg(5,1) + p(23)*plg(7,1) + &
             p(27)*plg(2,1) + p(15)*plg(4,1) + p(60)*plg(6,1)
      ! SYMMETRICAL ANNUAL
      t(3) = (p(19)+p(48)*plg(3,1)+p(30)*plg(5,1))*cd32
      ! SYMMETRICAL SEMIANNUAL
      t(4) = (p(16)+p(17)*plg(3,1)+p(31)*plg(5,1))*cd18
      ! ASYMMETRICAL ANNUAL
      t(5) = (p(10)*plg(2,1)+p(11)*plg(4,1)+p(21)*plg(6,1))*cd14
      ! ASYMMETRICAL SEMIANNUAL
      t(6) = (p(38)*plg(2,1))*cd39
      ! DIURNAL
      if ( abs(sw(7)) > d_zero ) then
        t71 = p(12)*plg(3,2)*cd14*swc(5)
        t72 = p(13)*plg(3,2)*cd14*swc(5)
        t(7) = ((p(4)*plg(2,2)+p(5)*plg(4,2)+t71)*ctloc + &
                (p(7)*plg(2,2)+p(8)*plg(4,2)+t72)*stloc)
      end if
      ! SEMIDIURNAL
      if ( abs(sw(8)) > d_zero ) then
        t81 = (p(24)*plg(4,3)+p(36)*plg(6,3))*cd14*swc(5)
        t82 = (p(34)*plg(4,3)+p(37)*plg(6,3))*cd14*swc(5)
        t(8) = ((p(6)*plg(3,3)+p(42)*plg(5,3)+t81)*c2tloc + &
                (p(9)*plg(3,3)+p(43)*plg(5,3)+t82)*s2tloc)
      end if
      ! TERDIURNAL
      if ( abs(sw(14)) > d_zero ) then
        t(14) = p(40)*plg(4,4)*s3tloc + p(41)*plg(4,4)*c3tloc
      end if
      ! MAGNETIC ACTIVITY
      if ( abs(sw(9)) > d_zero ) then
        if ( sw(9) > d_zero ) t(9) = apdf*(p(33)+p(46)*plg(3,1)*swc(2))
        if ( sw(9) < d_zero ) t(9) = (p(51)*apt(1)+p(97)*plg(3,1)*apt(1)*swc(2))
      end if
      ! LONGITUDINAL
      if ( abs(sw(10)) > d_zero .and. abs(sw(11)) > d_zero .and. &
           glong > -1000.0_rkx ) then
        t(11) = (d_one+plg(2,1)*(p(81)*swc(5)*cos(dr*(day-p(82))) +      &
                               p(86)*swc(6)*cos(d_two*dr*(day-p(87)))) + &
                               p(84)*swc(3)*cos(dr*(day-p(85))) +      &
                               p(88)*swc(4)*cos(d_two*dr*(day-p(89)))) * &
                ((p(65)*plg(3,2)+p(66)*plg(5,2)+p(67)*plg(7,2)+         &
                  p(75)*plg(2,2)+p(76)*plg(4,2)+p(77)*plg(6,2)) *       &
                  cos(degrad*glong) +                                   &
                 (p(91)*plg(3,2)+p(92)*plg(5,2)+p(93)*plg(7,2)+         &
                  p(78)*plg(2,2)+p(79)*plg(4,2)+p(80)*plg(6,2)) *       &
                  sin(degrad*glong))
      end if
      tt = d_zero
      do i = 1 , 14
        tt = tt + abs(sw(i))*t(i)
      end do
      glob7s = tt
    end function glob7s
!
!     FIND ALTITUDE OF PRESSURE SURFACE (PRESS) FROM GTD7 INPUT:
!     IYD - YEAR AND DAY AS YYDDD
!     SEC - UT(SEC)
!     GLAT - GEODETIC LATITUDE(DEG)
!     GLONG - GEODETIC LONGITUDE(DEG)
!     STL - LOCAL APPARENT SOLAR TIME(HRS)
!     F107A - 3 MONTH AVERAGE OF F10.7 FLUX
!     F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY
!     AP - MAGNETIC INDEX(DAILY) OR WHEN SW(9)=-1. :
!     - ARRAY CONTAINING:
!     (1) DAILY AP
!     (2) 3 HR AP INDEX FOR CURRENT TIME
!     (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
!     (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
!     (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
!     (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
!     TO CURRENT TIME
!     (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 59 HRS PRIOR
!     TO CURRENT TIME
!     PRESS - PRESSURE LEVEL(MB)
!     OUTPUT:
!     ALT - ALTITUDE(KM)
!     D(1) - HE NUMBER DENSITY(CM-3)
!     D(2) - O NUMBER DENSITY(CM-3)
!     D(3) - N2 NUMBER DENSITY(CM-3)
!     D(4) - O2 NUMBER DENSITY(CM-3)
!     D(5) - AR NUMBER DENSITY(CM-3)
!     D(6) - TOTAL MASS DENSITY(GM/CM3)
!     D(7) - H NUMBER DENSITY(CM-3)
!     D(8) - N NUMBER DENSITY(CM-3)
!     D(9) - HOT O NUMBER DENSITY(CM-3)
!     T(1) - EXOSPHERIC TEMPERATURE
!     T(2) - TEMPERATURE AT ALT
!
    subroutine ghp7(iyd,sec,alt,glat,glong,stl,f107a,f107,ap,d,t,press)
      implicit none
!
      real(rkx) , intent(in) :: f107 , f107a , glat , glong , press , sec , stl
      real(rkx) , intent(out) :: alt
      integer(ik4) , intent(in) :: iyd
      real(rkx) , intent(in) , dimension(7) :: ap
      real(rkx) , intent(inout) , dimension(9) :: d
      real(rkx) , intent(inout) , dimension(2) :: t
!
      real(rkx) :: ca , cd , cl , cl2 , diff , g , p , pl , sh , &
                 xm , xn , z , zi
      integer(ik4) :: iday , l

      real(rkx) , parameter :: bm = 1.3806e-19_rkx
      real(rkx) , parameter :: test = 0.00043_rkx
      integer(ik4) , parameter :: ltest = 12
!
      pl = log10(press)
!
!     Initial altitude estimate
      if ( pl >= -5.0_rkx ) then
        if ( pl > 2.5_rkx ) zi = 18.06_rkx*(3.00_rkx-pl)
        if ( pl > 0.75_rkx .and. pl <= 2.5_rkx ) zi = 14.98_rkx*(3.08_rkx-pl)
        if ( pl > -1.0_rkx .and. pl <= 0.75_rkx ) zi = 17.8_rkx*(2.72_rkx-pl)
        if ( pl > -2.0_rkx .and. pl <= -1.0_rkx ) zi = 14.28_rkx*(3.64_rkx-pl)
        if ( pl > -4.0_rkx .and. pl <= -2.0_rkx ) zi = 12.72_rkx*(4.32_rkx-pl)
        if ( pl <= -4.0_rkx ) zi = 25.3_rkx*(.11_rkx-pl)
        iday = mod(iyd,1000)
        cl = glat / 90.0_rkx
        cl2 = cl*cl
        if ( iday < 182 ) cd = d_one - iday/91.25_rkx
        if ( iday >= 182 ) cd = iday / 91.25_rkx - 3.0_rkx
        ca = d_zero
        if ( pl > -1.11_rkx .and. pl <= -0.23_rkx ) ca = d_one
        if ( pl > -0.23_rkx ) ca = (2.79_rkx-pl)/(2.79_rkx+0.23_rkx)
        if ( pl <= -1.11_rkx .and. pl > -3.0_rkx ) &
          ca = (-2.93_rkx-pl)/(-2.93_rkx+1.11_rkx)
        z = zi - 4.87_rkx*cl*cd*ca - 1.64_rkx*cl2*ca + 0.31_rkx*ca*cl
      end if
      if ( pl < -5.0_rkx ) z = 22.0_rkx*(pl+4.0_rkx)**2 + 110.0_rkx
      ! ITERATION LOOP
      l = 0
      do
        l = l + 1
        call gtd7(iyd,sec,z,glat,glong,stl,f107a,f107,ap,48,d,t)
        xn = d(1) + d(2) + d(3) + d(4) + d(5) + d(7) + d(8)
        p = bm*xn*t(2)
        if ( imr == 1 ) p = p*1.e-6_rkx
        diff = pl - log10(p)
        if ( abs(diff) < test .or. l == ltest ) then
          ! Non converging , should not happen
          alt = z
          exit
        else
          xm = d(6)/xn/1.66e-24_rkx
          if ( imr == 1 ) xm = xm*1.e3_rkx
          g = gsurf / (1.0_rkx+z/re)**2
          sh = r100gas*t(2)/(xm*g)
          ! New altitude estimate using scale height
          if ( l < 6 ) then
            z = z - sh*diff*2.302_rkx
          else
            z = z - sh*diff
          end if
          cycle
        end if
      end do
    end subroutine ghp7
!
    real(rkx) function globe7(yrd,sec,lat,long,tloc,f107a,f107,ap,p)
      implicit none
!
      real(rkx) , intent(in) :: f107 , f107a , lat , long , sec , tloc , yrd
      real(rkx) , dimension(:) , intent(in) :: ap
      real(rkx) , intent(out) , dimension(:) :: p
!
      real(rkx) :: c , c2 , c4 , exp1 , f1 , f2 , p44 , p45 , s , s2 ,    &
                  t71 , t72 , t81 , t82 , tix
      integer(ik4) :: i
      real(rkx) , save , dimension(25) :: sv = d_one
!
!     CALCULATE G(L) FUNCTION
!
!     Upper Thermosphere Parameters
      integer(ik4) , parameter :: nsw = 14
      real(rkx) , parameter :: hr = 0.2618_rkx
      real(rkx) , parameter :: sr = 7.2722e-5_rkx
      real(rkx) , save :: xl  = 1000.0_rkx
      real(rkx) , save :: tll = 1000.0_rkx
      real(rkx) , save :: sw9 = 1.0_rkx
      real(rkx) , save :: dayl = -1.0_rkx
      real(rkx) , save :: p14 = -1000.0_rkx
      real(rkx) , save :: p18 = -1000.0_rkx
      real(rkx) , save :: p32 = -1000.0_rkx
      real(rkx) , save :: p39 = -1000.0_rkx
      real(rkx) , save :: cd14 , cd18 , cd32 , cd39

      if ( isw /= 64999 ) call tselec(sv)
!
      t(:) = d_zero
!
      if ( sw(9) > d_zero ) sw9 =  d_one
      if ( sw(9) < d_zero ) sw9 = -d_one

      iyr = int(yrd*d_r1000)
      day = yrd - iyr*d_1000
!     Eq. A22 (remainder of code)
      if ( abs(xl-lat) > nearzero ) then
!       CALCULATE LEGENDRE POLYNOMIALS
        c = sin(lat*degrad)
        s = cos(lat*degrad)
        c2 = c*c
        c4 = c2*c2
        s2 = s*s
        plg(2,1) = c
        plg(3,1) = 0.5_rkx*(3.0_rkx*c2-d_one)
        plg(4,1) = 0.5_rkx*(5.0_rkx*c*c2-3.0_rkx*c)
        plg(5,1) = (35.0_rkx*c4-30.0_rkx*c2+3.0_rkx)/8.0_rkx
        plg(6,1) = (63.0_rkx*c2*c2*c-70.0_rkx*c2*c+15.0_rkx*c)/8.0_rkx
        plg(7,1) = (11.0_rkx*c*plg(6,1)-5.0_rkx*plg(5,1))/6.0_rkx
!       PLG(8,1) = (13.0_rkx*C*PLG(7,1) - 6.0_rkx*PLG(6,1))/7.0_rkx
        plg(2,2) = s
        plg(3,2) = 3.0_rkx*c*s
        plg(4,2) = 1.5_rkx*(5.0_rkx*c2-d_one)*s
        plg(5,2) = 2.5_rkx*(7.0_rkx*c2*c-3.0_rkx*c)*s
        plg(6,2) = 1.875_rkx*(21.0_rkx*c4-14.0_rkx*c2+d_one)*s
        plg(7,2) = (11.0_rkx*c*plg(6,2)-6.0_rkx*plg(5,2))/5.0_rkx
!       PLG(8,2) = (13.0_rkx*C*PLG(7,2)-7.0_rkx*PLG(6,2))/6.0_rkx
!       PLG(9,2) = (15.0_rkx*C*PLG(8,2)-8.0_rkx*PLG(7,2))/7.0_rkx
        plg(3,3) = 3.0_rkx*s2
        plg(4,3) = 15.0_rkx*s2*c
        plg(5,3) = 7.5_rkx*(7.0_rkx*c2-d_one)*s2
        plg(6,3) = 3.0_rkx*c*plg(5,3) - d_two*plg(4,3)
        plg(7,3) = (11.0_rkx*c*plg(6,3)-7.0_rkx*plg(5,3))/d_four
        plg(8,3) = (13.0_rkx*c*plg(7,3)-8.0_rkx*plg(6,3))/5.0_rkx
        plg(4,4) = 15.0_rkx*s2*s
        plg(5,4) = 105.0_rkx*s2*s*c
        plg(6,4) = (9.0_rkx*c*plg(5,4)-7.0_rkx*plg(4,4))/d_two
        plg(7,4) = (11.0_rkx*c*plg(6,4)-8.0_rkx*plg(5,4))/3.0_rkx
        xl = lat
      end if
      if ( abs(tll-tloc) > nearzero ) then
        if ( abs(sw(7)) > d_zero .or. abs(sw(8)) > d_zero .or. &
             abs(sw(14)) > d_zero ) then
          stloc = sin(hr*tloc)
          ctloc = cos(hr*tloc)
          s2tloc = sin(d_two*hr*tloc)
          c2tloc = cos(d_two*hr*tloc)
          s3tloc = sin(3.0_rkx*hr*tloc)
          c3tloc = cos(3.0_rkx*hr*tloc)
          tll = tloc
        end if
      end if
      if ( abs(day-dayl) > nearzero .or. abs(p(14)-p14) > nearzero ) &
        cd14 = cos(dr*(day-p(14)))
      if ( abs(day-dayl) > nearzero .or. abs(p(18)-p18) > nearzero ) &
        cd18 = cos(d_two*dr*(day-p(18)))
      if ( abs(day-dayl) > nearzero .or. abs(p(32)-p32) > nearzero ) &
        cd32 = cos(dr*(day-p(32)))
      if ( abs(day-dayl) > nearzero .or. abs(p(39)-p39) > nearzero ) &
        cd39 = cos(d_two*dr*(day-p(39)))
      dayl = day
      p14 = p(14)
      p18 = p(18)
      p32 = p(32)
      p39 = p(39)
      ! F10.7 EFFECT
      df = f107 - f107a
      dfa = f107a - 150.0_rkx
      t(1) = p(20)*df*(d_one+p(60)*dfa) + p(21)*df*df + p(22)*dfa+p(30)*dfa*dfa
      f1 = d_one + (p(48)*dfa+p(20)*df+p(21)*df*df)*swc(1)
      f2 = d_one + (p(50)*dfa+p(20)*df+p(21)*df*df)*swc(1)
      ! TIME INDEPENDENT
      t(2) = (p(2)*plg(3,1)+p(3)*plg(5,1)+p(23)*plg(7,1)) + &
             (p(15)*plg(3,1))*dfa*swc(1) + p(27)*plg(2,1)
      ! SYMMETRICAL ANNUAL
      t(3) = (p(19))*cd32
      ! SYMMETRICAL SEMIANNUAL
      t(4) = (p(16)+p(17)*plg(3,1))*cd18
      ! ASYMMETRICAL ANNUAL
      t(5) = f1*(p(10)*plg(2,1)+p(11)*plg(4,1))*cd14
      ! ASYMMETRICAL SEMIANNUAL
      t(6) = p(38)*plg(2,1)*cd39
      ! DIURNAL
      if ( abs(sw(7)) > d_zero ) then
        t71 = (p(12)*plg(3,2))*cd14*swc(5)
        t72 = (p(13)*plg(3,2))*cd14*swc(5)
        t(7) = f2*((p(4)*plg(2,2)+p(5)*plg(4,2)+p(28)*plg(6,2)+t71)*ctloc + &
                   (p(7)*plg(2,2)+p(8)*plg(4,2)+p(29)*plg(6,2)+t72)*stloc)
      end if
      ! SEMIDIURNAL
      if ( abs(sw(8)) > d_zero ) then
        t81 = (p(24)*plg(4,3)+p(36)*plg(6,3))*cd14*swc(5)
        t82 = (p(34)*plg(4,3)+p(37)*plg(6,3))*cd14*swc(5)
        t(8) = f2*((p(6)*plg(3,3)+p(42)*plg(5,3)+t81)*c2tloc + &
                   (p(9)*plg(3,3)+p(43)*plg(5,3)+t82)*s2tloc)
      end if
      ! TERDIURNAL
      if ( abs(sw(14)) > d_zero ) then
        t(14) = f2*((p(40)*plg(4,4) + &
              (p(94)*plg(5,4)+p(47)*plg(7,4))*cd14*swc(5))*s3tloc + &
              (p(41)*plg(4,4)+(p(95)*plg(5,4) + &
               p(49)*plg(7,4))*cd14*swc(5))*c3tloc)
      end if
      ! MAGNETIC ACTIVITY BASED ON DAILY AP
      if ( abs(sw9+d_one) > nearzero ) then
        apd = (ap(1)-d_four)
        p44 = p(44)
        p45 = p(45)
        if ( p44 < d_zero ) p44 = 1.e-5_rkx
        apdf = apd + (p45-d_one)*(apd+(exp(-p44*apd)-d_one)/p44)
        if ( abs(sw(9)) > nearzero ) then
          t(9) = apdf*(p(33)+p(46)*plg(3,1)+p(35)*plg(5,1) +           &
                      (p(101)*plg(2,1)+p(102)*plg(4,1) +               &
                       p(103)*plg(6,1))*cd14*swc(5)+(p(122)*plg(2,2) + &
                       p(123)*plg(4,2)+p(124)*plg(6,2))*swc(7)*        &
                       cos(hr*(tloc-p(125))))
        end if
      else if ( abs(p(52)) > nearzero ) then
        exp1 = exp(-10800.0_rkx*abs(p(52))/(d_one+p(139)*(45.0_rkx-abs(lat))))
        if ( exp1 > 0.99999_rkx ) exp1 = 0.99999_rkx
        if ( p(25) < 1.e-4_rkx ) p(25) = 1.e-4_rkx
        apt(1) = sg0(exp1)
!       APT(2)=SG2(EXP1)
!       APT(3)=SG0(EXP2)
!       APT(4)=SG2(EXP2)
        if ( abs(sw(9)) > d_zero ) then
          t(9) = apt(1)*(p(51)+p(97)*plg(3,1)+p(55)*plg(5,1) +   &
                        (p(126)*plg(2,1)+p(127)*plg(4,1) +       &
                         p(128)*plg(6,1))*cd14*swc(5) +          &
                        (p(129)*plg(2,2)+p(130)*plg(4,2) +       &
                         p(131)*plg(6,2))*swc(7)*cos(hr*(tloc-p(132))))
        end if
      end if
      if ( abs(sw(10)) > d_zero .and. long > -1000.0_rkx ) then
        ! LONGITUDINAL
        if ( abs(sw(11)) > d_zero ) then
          t(11) = (d_one+p(81)*dfa*swc(1)) *                               &
                    ((p(65)*plg(3,2)+p(66)*plg(5,2)+p(67)*plg(7,2) +     &
                      p(104)*plg(2,2)+p(105)*plg(4,2)+p(106)*plg(6,2) +  &
              swc(5)*(p(110)*plg(2,2)+p(111)*plg(4,2)+p(112)*plg(6,2)) * &
                  cd14)*cos(degrad*long) +                             &
                     (p(91)*plg(3,2)+p(92)*plg(5,2)+p(93)*plg(7,2) +     &
                      p(107)*plg(2,2)+p(108)*plg(4,2)+p(109)*plg(6,2) +  &
              swc(5)*(p(113)*plg(2,2)+p(114)*plg(4,2)+p(115)*plg(6,2)) * &
                  cd14)*sin(degrad*long))
        end if
        ! UT AND MIXED UT,LONGITUDE
        if ( abs(sw(12)) > nearzero ) then
          t(12) = (d_one+p(96)*plg(2,1))*(d_one+p(82)*dfa*swc(1)) *   &
                  (d_one+p(120)*plg(2,1)*swc(5)*cd14)  *            &
                ((p(69)*plg(2,1)+p(70)*plg(4,1)+p(71)*plg(6,1)) * &
                cos(sr*(sec-p(72))))
          t(12) = t(12) + swc(11) *                              &
                (p(77)*plg(4,3)+p(78)*plg(6,3)+p(79)*plg(8,3)) * &
                cos(sr*(sec-p(80))+d_two*degrad*long)*(d_one+p(138)*dfa*swc(1))
        end if
        ! UT,LONGITUDE MAGNETIC ACTIVITY
        if ( abs(sw(13)) > nearzero ) then
          if ( abs(sw9+d_one) > nearzero ) then
            t(13) = apdf*swc(11)*(d_one+p(121)*plg(2,1))                 &
                   *((p(61)*plg(3,2)+p(62)*plg(5,2)+p(63)*plg(7,2))    &
                   *cos(degrad*(long-p(64)))) + apdf*swc(11)*swc(5)  &
                   *(p(116)*plg(2,2)+p(117)*plg(4,2)+p(118)*plg(6,2))  &
                   *cd14*cos(degrad*(long-p(119))) + apdf*swc(12)    &
                   *(p(84)*plg(2,1)+p(85)*plg(4,1)+p(86)*plg(6,1))     &
                   *cos(sr*(sec-p(76)))
          else if ( p(52)/=0 ) then
            t(13) = apt(1)*swc(11)*(d_one+p(133)*plg(2,1))                &
                   *((p(53)*plg(3,2)+p(99)*plg(5,2)+p(68)*plg(7,2))     &
                   *cos(degrad*(long-p(98)))) + apt(1)*swc(11)*swc(5) &
                   *(p(134)*plg(2,2)+p(135)*plg(4,2)+p(136)*plg(6,2))   &
                   *cd14*cos(degrad*(long-p(137))) + apt(1)*swc(12)   &
                   *(p(56)*plg(2,1)+p(57)*plg(4,1)+p(58)*plg(6,1))      &
                   *cos(sr*(sec-p(59)))
          end if
        end if
      end if

      ! PARMS NOT USED: 83, 90,100,140-150
      tix = p(31)
      do i = 1 , nsw
        tix = tix + abs(sw(i))*t(i)
      end do
      globe7 = tix
      contains

!     3hr Magnetic activity functions
!     Eq. A24d
      real(rkx) function g0(a)
        implicit none
        real(rkx) , intent(in) :: a
        g0 = (a-d_four+(p(26)-d_one)*(a-d_four + &
             (exp(-abs(p(25))*(a-d_four))-d_one)/abs(p(25))))
      end function g0
!     Eq. A24c
      real(rkx) function sumex(ex)
        implicit none
        real(rkx) , intent(in) :: ex
        sumex = d_one + (d_one-ex**19)/(d_one-ex)*sqrt(ex)
      end function sumex
!     Eq. A24a
      real(rkx) function sg0(ex)
        implicit none
        real(rkx) , intent(in) :: ex
        sg0 = (g0(ap(2))+(g0(ap(3))*ex+g0(ap(4))*ex*ex + &
               g0(ap(5))*ex**3+(g0(ap(6))*ex**4 +    &
               g0(ap(7))*ex**12)*(d_one-ex**8)/(d_one-ex)))/sumex(ex)
      end function sg0
    end function globe7
!
    subroutine gtd7d(iyd,sec,alt,glat,glong,stl,f107a,f107,ap,mass,d,t)
      implicit none
!
      real(rkx) , intent(in) :: alt , f107 , f107a , glat , glong , sec , stl
      integer(ik4) :: iyd , mass
      real(rkx) , intent(in) , dimension(7) :: ap
      real(rkx) , intent(inout) , dimension(2) :: t
      real(rkx) , intent(inout) , dimension(9) :: d
!
!     NRLMSISE-00
!     -----------
!     This subroutine provides Effective Total Mass Density for
!     output D(6) which includes contributions from "anomalous
!     oxygen" which can affect satellite drag above 500 km.  This
!     subroutine is part of the distribution package for the
!     Neutral Atmosphere Empirical Model from the surface to lower
!     exosphere.  See subroutine GTD7 for more extensive comments.
!
!     INPUT VARIABLES:
!     IYD - YEAR AND DAY AS YYDDD (day of year from 1 to 365 (or 366))
!     (Year ignored in current model)
!     SEC - UT(SEC)
!     ALT - ALTITUDE(KM)
!     GLAT - GEODETIC LATITUDE(DEG)
!     GLONG - GEODETIC LONGITUDE(DEG)
!     STL - LOCAL APPARENT SOLAR TIME(HRS; see Note below)
!     F107A - 81 day AVERAGE OF F10.7 FLUX (centered on day DDD)
!     F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY
!     AP - MAGNETIC INDEX(DAILY) OR WHEN SW(9)=-1. :
!     - ARRAY CONTAINING:
!     (1) DAILY AP
!     (2) 3 HR AP INDEX FOR CURRENT TIME
!     (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
!     (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
!     (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
!     (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
!     TO CURRENT TIME
!     (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS PRIOR
!     TO CURRENT TIME
!     MASS - MASS NUMBER (ONLY DENSITY FOR SELECTED GAS IS
!     CALCULATED.  MASS 0 IS TEMPERATURE.  MASS 48 FOR ALL.
!     MASS 17 IS Anomalous O ONLY.)
!
!     NOTES ON INPUT VARIABLES:
!     UT, Local Time, and Longitude are used independently in the
!     model and are not of equal importance for every situation.
!     For the most physically realistic calculation these three
!     variables should be consistent (STL=SEC/3600+GLONG/15).
!     The Equation of Time departures from the above formula
!     for apparent local time can be included if available but
!     are of minor importance.
!
!     F107 and F107A values used to generate the model correspond
!     to the 10.7 cm radio flux at the actual distance of the Earth
!     from the Sun rather than the radio flux at 1 AU.
!
!     OUTPUT VARIABLES:
!     D(1) - HE NUMBER DENSITY(CM-3)
!     D(2) - O NUMBER DENSITY(CM-3)
!     D(3) - N2 NUMBER DENSITY(CM-3)
!     D(4) - O2 NUMBER DENSITY(CM-3)
!     D(5) - AR NUMBER DENSITY(CM-3)
!     D(6) - TOTAL MASS DENSITY(GM/CM3) [includes anomalous oxygen]
!     D(7) - H NUMBER DENSITY(CM-3)
!     D(8) - N NUMBER DENSITY(CM-3)
!     D(9) - Anomalous oxygen NUMBER DENSITY(CM-3)
!     T(1) - EXOSPHERIC TEMPERATURE
!     T(2) - TEMPERATURE AT ALT
!
      call gtd7(iyd,sec,alt,glat,glong,stl,f107a,f107,ap,mass,d,t)
!     TOTAL MASS DENSITY
      if ( mass == 48 ) then
        d(6) = 1.66e-24_rkx*(d_four*d(1)+16.0_rkx*d(2)+28.0_rkx*d(3) + &
                            32.0_rkx*d(4)+40.0_rkx*d(5)+d(7) +     &
                            14.0_rkx*d(8)+16.0_rkx*d(9))
        if ( imr == 1 ) d(6) = d(6)/1000.0_rkx
      end if
    end subroutine gtd7d
!
    subroutine gtd7(iyd,sec,alt,glat,glong,stl,f107a,f107,ap,mass,d,t)
      implicit none
!
      real(rkx) , intent(in) :: alt , f107 , f107a , glat , glong , sec , stl
      integer(ik4) , intent(in) :: iyd , mass
      real(rkx) , intent(in) , dimension(7) :: ap
      real(rkx) , intent(out) , dimension(9) :: d
      real(rkx) , intent(out) , dimension(2) :: t
!
      real(rkx) :: altt , dmc , dmr , dm28m , dz28 , tz , v1 , xlat , xmm
      real(rkx) , dimension(9) :: ds
      integer(ik4) :: j , mss
      real(rkx) , dimension(2) :: ts
!
      real(rkx) , dimension(25) , save :: sv
      integer(ik4) , save :: mssl
      real(rkx) , save :: alast
      real(rkx) , parameter :: zmix = 62.5_rkx
!
      data alast /99999.0_rkx/
      data mssl /-999/
      data sv /25*d_one/
!
!     NRLMSISE-00
!     -----------
!     Neutral Atmosphere Empirical Model from the surface to lower
!     exosphere
!
!     NEW FEATURES:
!     *Extensive satellite drag database used in model generation
!     *Revised O2 (and O) in lower thermosphere
!     *Additional nonlinear solar activity term
!     *"ANOMALOUS OXYGEN" NUMBER DENSITY, OUTPUT D(9)
!     At high altitudes (> 500 km), hot atomic oxygen or ionized
!     oxygen can become appreciable for some ranges of subroutine
!     inputs, thereby affecting drag on satellites and debris. We
!     group these species under the term "anomalous oxygen," since
!     their individual variations are not presently separable with
!     the drag data used to define this model component.
!
!     SUBROUTINES FOR SPECIAL OUTPUTS:
!
!     HIGH ALTITUDE DRAG: EFFECTIVE TOTAL MASS DENSITY
!     (SUBROUTINE GTD7D, OUTPUT D(6))
!     For atmospheric drag calculations at altitudes above 500 km,
!     call SUBROUTINE GTD7D to compute the "effective total mass
!     density" by including contributions from "anomalous oxygen."
!     See "NOTES ON OUTPUT VARIABLES" below on D(6).
!
!     PRESSURE GRID (SUBROUTINE GHP7)
!     See subroutine GHP7 to specify outputs at a pressure level
!     rather than at an altitude.
!
!     OUTPUT IN M-3 and KG/M3:   CALL METERS(.TRUE.)
!
!     INPUT VARIABLES:
!     IYD - YEAR AND DAY AS YYDDD (day of year from 1 to 365 (or 366))
!     (Year ignored in current model)
!     SEC - UT(SEC)
!     ALT - ALTITUDE(KM)
!     GLAT - GEODETIC LATITUDE(DEG)
!     GLONG - GEODETIC LONGITUDE(DEG)
!     STL - LOCAL APPARENT SOLAR TIME(HRS; see Note below)
!     F107A - 81 day AVERAGE OF F10.7 FLUX (centered on day DDD)
!     F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY
!     AP - MAGNETIC INDEX(DAILY) OR WHEN SW(9)=-1. :
!     - ARRAY CONTAINING:
!     (1) DAILY AP
!     (2) 3 HR AP INDEX FOR CURRENT TIME
!     (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
!     (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
!     (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
!     (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
!     TO CURRENT TIME
!     (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS PRIOR
!     TO CURRENT TIME
!     MASS - MASS NUMBER (ONLY DENSITY FOR SELECTED GAS IS
!     CALCULATED.  MASS 0 IS TEMPERATURE.  MASS 48 FOR ALL.
!     MASS 17 IS Anomalous O ONLY.)
!
!     NOTES ON INPUT VARIABLES:
!     UT, Local Time, and Longitude are used independently in the
!     model and are not of equal importance for every situation.
!     For the most physically realistic calculation these three
!     variables should be consistent (STL=SEC/3600+GLONG/15).
!     The Equation of Time departures from the above formula
!     for apparent local time can be included if available but
!     are of minor importance.
!
!     F107 and F107A values used to generate the model correspond
!     to the 10.7 cm radio flux at the actual distance of the Earth
!     from the Sun rather than the radio flux at 1 AU. The following
!     site provides both classes of values:
!     ftp://ftp.ngdc.noaa.gov/STP/SOLAR_DATA/SOLAR_RADIO/FLUX/
!
!     F107, F107A, and AP effects are neither large nor well
!     established below 80 km and these parameters should be set to
!     150., 150., and 4. respectively.
!
!     OUTPUT VARIABLES:
!     D(1) - HE NUMBER DENSITY(CM-3)
!     D(2) - O NUMBER DENSITY(CM-3)
!     D(3) - N2 NUMBER DENSITY(CM-3)
!     D(4) - O2 NUMBER DENSITY(CM-3)
!     D(5) - AR NUMBER DENSITY(CM-3)
!     D(6) - TOTAL MASS DENSITY(GM/CM3)
!     D(7) - H NUMBER DENSITY(CM-3)
!     D(8) - N NUMBER DENSITY(CM-3)
!     D(9) - Anomalous oxygen NUMBER DENSITY(CM-3)
!     T(1) - EXOSPHERIC TEMPERATURE
!     T(2) - TEMPERATURE AT ALT
!
!     NOTES ON OUTPUT VARIABLES:
!     TO GET OUTPUT IN M-3 and KG/M3:   CALL METERS(.TRUE.)
!
!     O, H, and N are set to zero below 72.5 km
!
!     T(1), Exospheric temperature, is set to global average for
!     altitudes below 120 km. The 120 km gradient is left at global
!     average value for altitudes below 72 km.
!
!     D(6), TOTAL MASS DENSITY, is NOT the same for subroutines GTD7
!     and GTD7D
!
!     SUBROUTINE GTD7 -- D(6) is the sum of the mass densities of the
!     species labeled by indices 1-5 and 7-8 in output variable D.
!     This includes He, O, N2, O2, Ar, H, and N but does NOT include
!     anomalous oxygen (species index 9).
!
!     SUBROUTINE GTD7D -- D(6) is the "effective total mass density
!     for drag" and is the sum of the mass densities of all species
!     in this model, INCLUDING anomalous oxygen.
!
!     SWITCHES: The following is for test and special purposes:
!
!     TO TURN ON AND OFF PARTICULAR VARIATIONS CALL TSELEC(SW),
!     WHERE SW IS A 25 ELEMENT ARRAY CONTAINING 0. FOR OFF, 1.
!     FOR ON, OR 2. FOR MAIN EFFECTS OFF BUT CROSS TERMS ON
!     FOR THE FOLLOWING VARIATIONS
!     1 - F10.7 EFFECT ON MEAN  2 - TIME INDEPENDENT
!     3 - SYMMETRICAL ANNUAL    4 - SYMMETRICAL SEMIANNUAL
!     5 - ASYMMETRICAL ANNUAL   6 - ASYMMETRICAL SEMIANNUAL
!     7 - DIURNAL               8 - SEMIDIURNAL
!     9 - DAILY AP             10 - ALL UT/LONG EFFECTS
!     11 - LONGITUDINAL         12 - UT AND MIXED UT/LONG
!     13 - MIXED AP/UT/LONG     14 - TERDIURNAL
!     15 - DEPARTURES FROM DIFFUSIVE EQUILIBRIUM
!     16 - ALL TINF VAR         17 - ALL TLB VAR
!     18 - ALL TN1 VAR           19 - ALL S VAR
!     20 - ALL TN2 VAR           21 - ALL NLB VAR
!     22 - ALL TN3 VAR           23 - TURBO SCALE HEIGHT VAR
!
!     To get current values of SW: CALL TRETRV(SW)
!
      if ( isw /= 64999 ) call tselec(sv)
!
      ! Test for changed input
      v1 = vtst7(iyd,sec,glat,glong,stl,f107a,f107,ap,1)
      !
      ! Latitude variation of gravity (none for SW(2)=0)
      xlat = glat
      if ( abs(sw(2)) < nearzero ) xlat = 45.0_rkx
      call glatf(xlat)
!
      xmm = pdm(5,3)
!
      ! THERMOSPHERE/MESOSPHERE (above ZN2(1))
      altt = max(alt,zn2(1))
      mss = mass
      ! Only calculate N2 in thermosphere if alt in mixed region
      if ( alt < zmix .and. mass > 0 ) mss = 28
      ! Only calculate thermosphere if input parameters changed
      ! or altitude above ZN2(1) in mesosphere
      if ( abs(v1-d_one) < nearzero .or. alt > zn2(1) .or. &
           alast > zn2(1) .or. mss /= mssl ) then
        call gts7(iyd,sec,altt,glat,glong,stl,f107a,f107,ap,mss,ds,ts)
        dm28m = dm28
!       metric adjustment
        if ( imr == 1 ) dm28m = dm28*1.e6_rkx
        mssl = mss
      end if
      t(1) = ts(1)
      t(2) = ts(2)
      if ( alt >= zn2(1) ) then
        do j = 1 , 9
          d(j) = ds(j)
        end do
        return
      end if
      !
      ! LOWER MESOSPHERE/UPPER STRATOSPHERE [between ZN3(1) and ZN2(1)]
      ! Temperature at nodes and gradients at end nodes
      ! Inverse temperature a linear function of spherical harmonics
      ! Only calculate nodes if input changed
      !
      if ( abs(v1-d_one) < nearzero .or. alast >= zn2(1) ) then
        tgn2(1) = tgn1(2)
        tn2(1) = tn1(5)
        tn2(2) = pma(1,1)*pavgm(1)/(d_one-sw(20)*glob7s(glong,pma(:,1)))
        tn2(3) = pma(1,2)*pavgm(2)/(d_one-sw(20)*glob7s(glong,pma(:,2)))
        tn2(4) = pma(1,3)*pavgm(3)/(d_one-sw(20)*sw(22)*glob7s(glong,pma(:,3)))
        tgn2(2) = pavgm(9)*pma(1,10)                                    &
                 *(d_one+sw(20)*sw(22)*glob7s(glong,pma(:,10)))*tn2(4)*tn2(4)   &
                 /(pma(1,3)*pavgm(3))**2
        tn3(1) = tn2(4)
      end if
      !
      ! LOWER STRATOSPHERE AND TROPOSPHERE [below ZN3(1)]
      ! Temperature at nodes and gradients at end nodes
      ! Inverse temperature a linear function of spherical harmonics
      ! Only calculate nodes if input changed
      !
      if ( alt < zn3(1) ) then
        if ( abs(v1-d_one) < nearzero .or. alast >= zn3(1) ) then
          tgn3(1) = tgn2(2)
          tn3(2) = pma(1,4)*pavgm(4)/(d_one-sw(22)*glob7s(glong,pma(:,4)))
          tn3(3) = pma(1,5)*pavgm(5)/(d_one-sw(22)*glob7s(glong,pma(:,5)))
          tn3(4) = pma(1,6)*pavgm(6)/(d_one-sw(22)*glob7s(glong,pma(:,6)))
          tn3(5) = pma(1,7)*pavgm(7)/(d_one-sw(22)*glob7s(glong,pma(:,7)))
          tgn3(2) = pma(1,8)*pavgm(8)*(d_one+sw(22)*glob7s(glong,pma(:,8)))  &
                   *tn3(5)*tn3(5)/(pma(1,7)*pavgm(7))**2
        end if
      end if

      if ( mass == 0 ) then
        dd = densm(alt,d_one,d_zero,tz)
        t(2) = tz
      else
        ! LINEAR TRANSITION TO FULL MIXING BELOW ZN2(1)
        dmc = d_zero
        if ( alt > zmix ) dmc = d_one - (zn2(1)-alt)/(zn2(1)-zmix)
        dz28 = ds(3)
        ! N2 DENSITY
        dmr = ds(3)/dm28m - d_one
        d(3) = densm(alt,dm28m,xmm,tz)
        d(3) = d(3)*(d_one+dmr*dmc)
        ! HE DENSITY
        d(1) = d_zero
        if ( mass == 4 .or. mass == 48 ) then
          dmr = ds(1)/(dz28*pdm(2,1)) - d_one
          d(1) = d(3)*pdm(2,1)*(d_one+dmr*dmc)
        end if
        ! O DENSITY
        d(2) = d_zero
        d(9) = d_zero
        ! O2 DENSITY
        d(4) = d_zero
        if ( mass == 32 .or. mass == 48 ) then
          dmr = ds(4)/(dz28*pdm(2,4)) - d_one
          d(4) = d(3)*pdm(2,4)*(d_one+dmr*dmc)
        end if
        ! AR DENSITY
        d(5) = d_zero
        if ( mass == 40 .or. mass == 48 ) then
          dmr = ds(5)/(dz28*pdm(2,5)) - d_one
          d(5) = d(3)*pdm(2,5)*(d_one+dmr*dmc)
        end if
        ! HYDROGEN DENSITY
        d(7) = d_zero
        ! ATOMIC NITROGEN DENSITY
        d(8) = d_zero
        !
        ! TOTAL MASS DENSITY
        !
        if ( mass == 48 ) then
          d(6) = 1.66e-24_rkx*(d_four*d(1)+16.0_rkx*d(2)+28.0_rkx*d(3) + &
                              32.0_rkx*d(4)+40.0_rkx*d(5)+d(7)+14.0_rkx*d(8))
          if ( imr == 1 ) d(6) = d(6)/1000.0_rkx
        end if
        t(2) = tz
      end if
      alast = alt
    end subroutine gtd7
!
    subroutine gts7(iyd,sec,alt,glat,glong,stl,f107a,f107,ap,mass,d,t)
      implicit none
!
      real(rkx) , intent(in) :: alt , f107 , f107a , glat , glong , sec , stl
      integer(ik4) , intent(in) :: iyd , mass
      real(rkx) , intent(in) , dimension(:) :: ap
      real(rkx) , intent(inout) , dimension(9) :: d
      real(rkx) , intent(inout) , dimension(2) :: t
!
      real(rkx) :: b01 , b04 , b14 , b16 , b28 , b32 , b40 ,     &
              day , db16h , ddum , g1 , g14 , g16 , g16h , g28 ,        &
              g32 , g4 , g40 , hc01 , hc04 , hc14 , hc16 , hc216 ,      &
              hc32 , hc40 , hcc01 , hcc14 , hcc16 , hcc232 , hcc32 ,    &
              rc01 , rc14 , rc16 , rc32 , t2 , tho , tz , v2 ,   &
              xmd , xmm , yrd , z , zc01 , zc04 , zc14 , zc16 , zc32 ,  &
              zc40 , zcc01 , zcc14 , zcc16 , zcc32 , zh01 , zh04 ,      &
              zh14 , zh16 , zh28 , zh32 , zh40 , zhf , zhm01 , zhm04 ,  &
              zhm14 , zhm16 , zhm28 , zhm32 , zhm40 , zmho , zsho , zsht
      integer(ik4) :: i , j
      real(rkx) , save :: alast
!
!     Thermospheric portion of NRLMSISE-00
!     See GTD7 for more extensive comments
!
!     OUTPUT IN M-3 and KG/M3:   CALL METERS(.TRUE.)
!
!     INPUT VARIABLES:
!     IYD - YEAR AND DAY AS YYDDD (day of year from 1 to 365 (or 366))
!     (Year ignored in current model)
!     SEC - UT(SEC)
!     ALT - ALTITUDE(KM) (>72.5 km)
!     GLAT - GEODETIC LATITUDE(DEG)
!     GLONG - GEODETIC LONGITUDE(DEG)
!     STL - LOCAL APPARENT SOLAR TIME(HRS; see Note below)
!     F107A - 81 day AVERAGE OF F10.7 FLUX (centered on day DDD)
!     F107 - DAILY F10.7 FLUX FOR PREVIOUS DAY
!     AP - MAGNETIC INDEX(DAILY) OR WHEN SW(9)=-1. :
!     - ARRAY CONTAINING:
!     (1) DAILY AP
!     (2) 3 HR AP INDEX FOR CURRENT TIME
!     (3) 3 HR AP INDEX FOR 3 HRS BEFORE CURRENT TIME
!     (4) 3 HR AP INDEX FOR 6 HRS BEFORE CURRENT TIME
!     (5) 3 HR AP INDEX FOR 9 HRS BEFORE CURRENT TIME
!     (6) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 12 TO 33 HRS PRIOR
!     TO CURRENT TIME
!     (7) AVERAGE OF EIGHT 3 HR AP INDICIES FROM 36 TO 57 HRS PRIOR
!     TO CURRENT TIME
!     MASS - MASS NUMBER (ONLY DENSITY FOR SELECTED GAS IS
!     CALCULATED.  MASS 0 IS TEMPERATURE.  MASS 48 FOR ALL.
!     MASS 17 IS Anomalous O ONLY.)
!
!     NOTES ON INPUT VARIABLES:
!     UT, Local Time, and Longitude are used independently in the
!     model and are not of equal importance for every situation.
!     For the most physically realistic calculation these three
!     variables should be consistent (STL=SEC/3600+GLONG/15).
!     The Equation of Time departures from the above formula
!     for apparent local time can be included if available but
!     are of minor importance.
!
!     F107 and F107A values used to generate the model correspond
!     to the 10.7 cm radio flux at the actual distance of the Earth
!     from the Sun rather than the radio flux at 1 AU. The following
!     site provides both classes of values:
!     ftp://ftp.ngdc.noaa.gov/STP/SOLAR_DATA/SOLAR_RADIO/FLUX/
!
!     F107, F107A, and AP effects are neither large nor well
!     established below 80 km and these parameters should be set to
!     150., 150., and 4. respectively.
!
!     OUTPUT VARIABLES:
!     D(1) - HE NUMBER DENSITY(CM-3)
!     D(2) - O NUMBER DENSITY(CM-3)
!     D(3) - N2 NUMBER DENSITY(CM-3)
!     D(4) - O2 NUMBER DENSITY(CM-3)
!     D(5) - AR NUMBER DENSITY(CM-3)
!     D(6) - TOTAL MASS DENSITY(GM/CM3) [Anomalous O NOT included]
!     D(7) - H NUMBER DENSITY(CM-3)
!     D(8) - N NUMBER DENSITY(CM-3)
!     D(9) - Anomalous oxygen NUMBER DENSITY(CM-3)
!     T(1) - EXOSPHERIC TEMPERATURE
!     T(2) - TEMPERATURE AT ALT
!
      data alast / -999.0_rkx/
!
      ! Test for changed input
      v2 = vtst7(iyd,sec,glat,glong,stl,f107a,f107,ap,2)
!
      yrd = iyd
      za = pdl(16,2)
      zn1(1) = za
      do j = 1 , 9
        d(j) = d_zero
      end do
      ! TINF VARIATIONS NOT IMPORTANT BELOW ZA OR ZN1(1)
      if ( alt > zn1(1) ) then
        if ( abs(v2-d_one) < nearzero .or. alast <= zn1(1) ) then
          tinf = ptm(1)*pt(1)*(d_one+sw(16)* &
             globe7(yrd,sec,glat,glong,stl,f107a,f107,ap,pt))
        end if
      else
        tinf = ptm(1)*pt(1)
      end if
      t(1) = tinf
      ! GRADIENT VARIATIONS NOT IMPORTANT BELOW ZN1(5)
      if ( alt > zn1(5) ) then
        if ( abs(v2-d_one) < nearzero .or. alast <= zn1(5) ) then
          xg0 = ptm(4)*ps(1)*(d_one+sw(19) * &
             globe7(yrd,sec,glat,glong,stl,f107a,f107,ap,ps))
        end if
      else
        xg0 = ptm(4)*ps(1)
      end if
      ! Calculate these temperatures only if input changed
      if ( abs(v2-d_one) < nearzero .or. alt < 300.0_rkx ) then
        tlb = ptm(2)*(d_one+sw(17)* &
          globe7(yrd,sec,glat,glong,stl,f107a,f107,ap,pd(:,4)))*pd(1,4)
      end if
      s = xg0/(tinf-tlb)
      ! Lower thermosphere temp variations not significant for
      ! density above 300 km
      if ( alt >= 300.0_rkx ) then
        tn1(2) = ptm(7)*ptl(1,1)
        tn1(3) = ptm(3)*ptl(1,2)
        tn1(4) = ptm(8)*ptl(1,3)
        tn1(5) = ptm(5)*ptl(1,4)
        tgn1(2) = ptm(9)*pma(1,9)*tn1(5)*tn1(5)/(ptm(5)*ptl(1,4))**2
      else if ( abs(v2-d_one) < nearzero .or. alast >= 300.0_rkx ) then
        tn1(2) = ptm(7)*ptl(1,1)/(d_one-sw(18)*glob7s(glong,ptl(:,1)))
        tn1(3) = ptm(3)*ptl(1,2)/(d_one-sw(18)*glob7s(glong,ptl(:,2)))
        tn1(4) = ptm(8)*ptl(1,3)/(d_one-sw(18)*glob7s(glong,ptl(:,3)))
        tn1(5) = ptm(5)*ptl(1,4)/(d_one-sw(18)*sw(20)*glob7s(glong,ptl(:,4)))
        tgn1(2) = ptm(9)*pma(1,9)*(d_one+sw(18)*sw(20)*glob7s(glong,pma(:,9)))   &
                & *tn1(5)*tn1(5)/(ptm(5)*ptl(1,4))**2
      end if
!
      z0 = zn1(4)
      t0 = tn1(4)
      tr12 = d_one
!
      if ( mass == 0 ) go to 700
!
      ! N2 variation factor at Zlb
      g28 = sw(21)*globe7(yrd,sec,glat,glong,stl,f107a,f107,ap,pd(:,3))
      day = mod(yrd,1000.0_rkx)
      ! VARIATION OF TURBOPAUSE HEIGHT
      zhf = pdl(25,2)*(d_one+sw(5) * &
            pdl(25,1)*sin(degrad*glat)*cos(dr*(day-pt(14))))
      yrd = iyd
      t(1) = tinf
      xmm = pdm(5,3)
      z = alt
!
      do j = 1 , 11
        if ( mass == mt(j) ) go to 100
      end do
!
      go to 800
!
 100  continue
!
      if ( z <= altl(6) .or. mass == 28 .or. mass == 48 ) then
        !
        !       **** N2 DENSITY ****
        !
        ! Diffusive density at Zlb
        db28 = pdm(1,3)*exp(g28)*pd(1,3)
        ! Diffusive density at Alt
        d(3) = densu(z,db28,tinf,tlb,28.0_rkx,alph(3),t(2),ptm(6),s)
        dd = d(3)
        ! Turbopause
        zh28 = pdm(3,3)*zhf
        zhm28 = pdm(4,3)*pdl(6,2)
        xmd = 28.0_rkx - xmm
        ! Mixed density at Zlb
        b28 = densu(zh28,db28,tinf,tlb,xmd,alph(3)-d_one,tz,ptm(6),s)
        if ( z <= altl(3) .and. abs(sw(15)) > d_zero ) then
          ! Mixed density at Alt
          dm28 = densu(z,b28,tinf,tlb,xmm,alph(3),tz,ptm(6),s)
          ! Net density at Alt
          d(3) = dnet(d(3),dm28,zhm28,xmm,28.0_rkx)
        end if
      end if
      if ( j == 2 ) go to 700
      if ( j == 4 .or. j == 9 ) then
      else if ( j == 5 ) then
        go to 800
      else if ( j == 6 ) then
        go to 200
      else if ( j == 7 ) then
        go to 300
      else if ( j == 8 ) then
        go to 400
      else if ( j == 10 ) then
        go to 500
      else if ( j == 11 ) then
        go to 600
      else
        !
        !       **** HE DENSITY ****
        !
        ! Density variation factor at Zlb
        g4 = sw(21)*globe7(yrd,sec,glat,glong,stl,f107a,f107,ap,pd(:,1))
        ! Diffusive density at Zlb
        db04 = pdm(1,1)*exp(g4)*pd(1,1)
        ! Diffusive density at Alt
        d(1) = densu(z,db04,tinf,tlb,d_four,alph(1),t(2),ptm(6),s)
        dd = d(1)
        if ( z <= altl(1) .and. abs(sw(15)) > d_zero ) then
          ! Turbopause
          zh04 = pdm(3,1)
          ! Mixed density at Zlb
          b04 = densu(zh04,db04,tinf,tlb,d_four-xmm,alph(1)-d_one,t(2),ptm(6),s)
          ! Mixed density at Alt
          dm04 = densu(z,b04,tinf,tlb,xmm,d_zero,t(2),ptm(6),s)
          zhm04 = zhm28
          ! Net density at Alt
          d(1) = dnet(d(1),dm04,zhm04,xmm,d_four)
          ! Correction to specified mixing ratio at ground
          rl = log(b28*pdm(2,1)/b04)
          zc04 = pdm(5,1)*pdl(1,2)
          hc04 = pdm(6,1)*pdl(2,2)
          ! Net density corrected at Alt
          d(1) = d(1)*ccor(z,rl,hc04,zc04)
        end if
        if ( mass /= 48 ) go to 800
      end if
      !
      !     **** O DENSITY ****
      !
      ! Density variation factor at Zlb
      g16 = sw(21)*globe7(yrd,sec,glat,glong,stl,f107a,f107,ap,pd(:,2))
      ! Diffusive density at Zlb
      db16 = pdm(1,2)*exp(g16)*pd(1,2)
      ! Diffusive density at Alt
      d(2) = densu(z,db16,tinf,tlb,16.0_rkx,alph(2),t(2),ptm(6),s)
      dd = d(2)
      if ( z <= altl(2) .and. abs(sw(15)) > d_zero ) then
        ! Corrected from PDM(3,1) to PDM(3,2)  12/2/85
        ! Turbopause
        zh16 = pdm(3,2)
        ! Mixed density at Zlb
        b16 = densu(zh16,db16,tinf,tlb,16.0_rkx-xmm,alph(2)-d_one,t(2),ptm(6),s)
        ! Mixed density at Alt
        dm16 = densu(z,b16,tinf,tlb,xmm,d_zero,t(2),ptm(6),s)
        zhm16 = zhm28
        ! Net density at Alt
        d(2) = dnet(d(2),dm16,zhm16,xmm,16.0_rkx)
        ! 3/16/99 Change form to match O2 departure from diff equil near
        ! 150 km and add dependence on F10.7
        ! RL=ALOG(B28*PDM(2,2)*ABS(PDL(17,2))/B16)
        rl = pdm(2,2)*pdl(17,2)*(d_one+sw(1)*pdl(24,1)*(f107a-150.0_rkx))
        hc16 = pdm(6,2)*pdl(4,2)
        zc16 = pdm(5,2)*pdl(3,2)
        hc216 = pdm(6,2)*pdl(5,2)
        d(2) = d(2)*ccor2(z,rl,hc16,zc16,hc216)
        ! Chemistry correction
        hcc16 = pdm(8,2)*pdl(14,2)
        zcc16 = pdm(7,2)*pdl(13,2)
        rc16 = pdm(4,2)*pdl(15,2)
        ! Net density corrected at Alt
        d(2) = d(2)*ccor(z,rc16,hcc16,zcc16)
      end if
      if ( mass /= 48 .and. mass /= 49 ) go to 800
      !
      !     **** O2 DENSITY ****
      !
      ! Density variation factor at Zlb
 200  continue
      g32 = sw(21)*globe7(yrd,sec,glat,glong,stl,f107a,f107,ap,pd(:,5))
      ! Diffusive density at Zlb
      db32 = pdm(1,4)*exp(g32)*pd(1,5)
      ! Diffusive density at Alt
      d(4) = densu(z,db32,tinf,tlb,32.0_rkx,alph(4),t(2),ptm(6),s)
      if ( mass == 49 ) then
        dd = dd + d_two*d(4)
      else
        dd = d(4)
      end if
      if ( abs(sw(15)) > d_zero ) then
        if ( z <= altl(4) ) then
          ! Turbopause
          zh32 = pdm(3,4)
          ! Mixed density at Zlb
          b32 = densu(zh32,db32,tinf,tlb,32.0_rkx-xmm,alph(4)-d_one, t(2),ptm(6),s)
          ! Mixed density at Alt
          dm32 = densu(z,b32,tinf,tlb,xmm,d_zero,t(2),ptm(6),s)
          zhm32 = zhm28
          ! Net density at Alt
          d(4) = dnet(d(4),dm32,zhm32,xmm,32.0_rkx)
          ! Correction to specified mixing ratio at ground
          rl = log(b28*pdm(2,4)/b32)
          hc32 = pdm(6,4)*pdl(8,2)
          zc32 = pdm(5,4)*pdl(7,2)
          d(4) = d(4)*ccor(z,rl,hc32,zc32)
        end if
        ! Correction for general departure from diffusive equilibrium
        ! above Zlb
        hcc32 = pdm(8,4)*pdl(23,2)
        hcc232 = pdm(8,4)*pdl(23,1)
        zcc32 = pdm(7,4)*pdl(22,2)
        rc32 = pdm(4,4)*pdl(24,2)*(d_one+sw(1)*pdl(24,1)*(f107a-150.0_rkx))
        ! Net density corrected at Alt
        d(4) = d(4)*ccor2(z,rc32,hcc32,zcc32,hcc232)
      end if
      if ( mass /= 48 ) go to 800
      !
      !     **** AR DENSITY ****
      !
      ! Density variation factor at Zlb
 300  continue
      g40 = sw(21)*globe7(yrd,sec,glat,glong,stl,f107a,f107,ap,pd(:,6))
      ! Diffusive density at Zlb
      db40 = pdm(1,5)*exp(g40)*pd(1,6)
      ! Diffusive density at Alt
      d(5) = densu(z,db40,tinf,tlb,40.0_rkx,alph(5),t(2),ptm(6),s)
      dd = d(5)
      if ( z <= altl(5) .and. abs(sw(15)) > d_zero ) then
        ! Turbopause
        zh40 = pdm(3,5)
        ! Mixed density at Zlb
        b40 = densu(zh40,db40,tinf,tlb,40.0_rkx-xmm,alph(5)-d_one,t(2),ptm(6),s)
        ! Mixed density at Alt
        dm40 = densu(z,b40,tinf,tlb,xmm,d_zero,t(2),ptm(6),s)
        zhm40 = zhm28
        ! Net density at Alt
        d(5) = dnet(d(5),dm40,zhm40,xmm,40.0_rkx)
        ! Correction to specified mixing ratio at ground
        rl = log(b28*pdm(2,5)/b40)
        hc40 = pdm(6,5)*pdl(10,2)
        zc40 = pdm(5,5)*pdl(9,2)
        ! Net density corrected at Alt
        d(5) = d(5)*ccor(z,rl,hc40,zc40)
      end if
      if ( mass /= 48 ) go to 800
      !
      !     **** HYDROGEN DENSITY ****
      !
      ! Density variation factor at Zlb
 400  continue
      g1 = sw(21)*globe7(yrd,sec,glat,glong,stl,f107a,f107,ap,pd(:,7))
      ! Diffusive density at Zlb
      db01 = pdm(1,6)*exp(g1)*pd(1,7)
      ! Diffusive density at Alt
      d(7) = densu(z,db01,tinf,tlb,d_one,alph(7),t(2),ptm(6),s)
      dd = d(7)
      if ( z <= altl(7) .and. abs(sw(15)) > d_zero ) then
        ! Turbopause
        zh01 = pdm(3,6)
        ! Mixed density at Zlb
        b01 = densu(zh01,db01,tinf,tlb,d_one-xmm,alph(7)-d_one,t(2),ptm(6),s)
        ! Mixed density at Alt
        dm01 = densu(z,b01,tinf,tlb,xmm,d_zero,t(2),ptm(6),s)
        zhm01 = zhm28
        ! Net density at Alt
        d(7) = dnet(d(7),dm01,zhm01,xmm,d_one)
        ! Correction to specified mixing ratio at ground
        rl = log(b28*pdm(2,6)*abs(pdl(18,2))/b01)
        hc01 = pdm(6,6)*pdl(12,2)
        zc01 = pdm(5,6)*pdl(11,2)
        d(7) = d(7)*ccor(z,rl,hc01,zc01)
        ! Chemistry correction
        hcc01 = pdm(8,6)*pdl(20,2)
        zcc01 = pdm(7,6)*pdl(19,2)
        rc01 = pdm(4,6)*pdl(21,2)
        ! Net density corrected at Alt
        d(7) = d(7)*ccor(z,rc01,hcc01,zcc01)
      end if
      if ( mass /= 48 ) go to 800
      !
      !     **** ATOMIC NITROGEN DENSITY ****
      !
      ! Density variation factor at Zlb
 500  continue
      g14 = sw(21)*globe7(yrd,sec,glat,glong,stl,f107a,f107,ap,pd(:,8))
      ! Diffusive density at Zlb
      db14 = pdm(1,7)*exp(g14)*pd(1,8)
      ! Diffusive density at Alt
      d(8) = densu(z,db14,tinf,tlb,14.0_rkx,alph(8),t(2),ptm(6),s)
      dd = d(8)
      if ( z <= altl(8) .and. abs(sw(15)) > d_zero ) then
        ! Turbopause
        zh14 = pdm(3,7)
        ! Mixed density at Zlb
        b14 = densu(zh14,db14,tinf,tlb,14.0_rkx-xmm,alph(8)-d_one,t(2),ptm(6),s)
        ! Mixed density at Alt
        dm14 = densu(z,b14,tinf,tlb,xmm,d_zero,t(2),ptm(6),s)
        zhm14 = zhm28
        ! Net density at Alt
        d(8) = dnet(d(8),dm14,zhm14,xmm,14.0_rkx)
        ! Correction to specified mixing ratio at ground
        rl = log(b28*pdm(2,7)*abs(pdl(3,1))/b14)
        hc14 = pdm(6,7)*pdl(2,1)
        zc14 = pdm(5,7)*pdl(1,1)
        d(8) = d(8)*ccor(z,rl,hc14,zc14)
        ! Chemistry correction
        hcc14 = pdm(8,7)*pdl(5,1)
        zcc14 = pdm(7,7)*pdl(4,1)
        rc14 = pdm(4,7)*pdl(6,1)
        ! Net density corrected at Alt
        d(8) = d(8)*ccor(z,rc14,hcc14,zcc14)
      end if
      if ( mass /= 48 ) go to 800
      !
      !     **** Anomalous OXYGEN DENSITY ****
      !
 600  continue
      g16h = sw(21)*globe7(yrd,sec,glat,glong,stl,f107a,f107,ap,pd(:,9))
      db16h = pdm(1,8)*exp(g16h)*pd(1,9)
      tho = pdm(10,8)*pdl(7,1)
      dd = densu(z,db16h,tho,tho,16.0_rkx,alph(9),t2,ptm(6),s)
      zsht = pdm(6,8)
      zmho = pdm(5,8)
      zsho = scalh(zmho,16.0_rkx,tho)
      d(9) = dd*exp(-zsht/zsho*(exp(-(z-zmho)/zsht)-d_one))
      if ( mass == 48 ) then
        !
        ! TOTAL MASS DENSITY
        !
        d(6) = 1.66e-24_rkx*(d_four*d(1)+16.0_rkx*d(2)+28.0_rkx*d(3) + &
                            32.0_rkx*d(4)+40.0_rkx*d(5)+d(7)+14.0_rkx*d(8))
        db48 = 1.66e-24_rkx*(d_four*db04+16.0_rkx*db16+28.0_rkx*db28+32.0_rkx*db32 + &
                         40.0_rkx*db40+db01+14.0_rkx*db14)
      end if
      go to 800
      ! TEMPERATURE AT ALTITUDE
 700  continue
      z = abs(alt)
      ddum = densu(z,d_one,tinf,tlb,d_zero,d_zero,t(2),ptm(6),s)
      ! ADJUST DENSITIES FROM CGS TO KGM
 800  continue
      if ( imr == 1 ) then
        do i = 1 , 9
          d(i) = d(i)*1.e6_rkx
        end do
        d(6) = d(6)/1000.0_rkx
      end if
      alast = alt
      return

      contains
      !
      !   Calculate scale height (km)
      !
      real(rkx) function scalh(alt,xm,temp)
        implicit none
        real(rkx) , intent(in) :: alt , temp , xm
        real(rkx) :: g
        g = gsurf/(d_one+alt/re)**2
        scalh = r100gas*temp/(g*xm)
      end function scalh
    end subroutine gts7
!
    real(rkx) function vtst7(iyd,sec,glat,glong,stl,f107a,f107,ap,ic)
      implicit none
!
      real(rkx) , intent(in) :: f107 , f107a , glat , glong , sec , stl
      integer(ik4) , intent(in) :: ic , iyd
      real(rkx) , intent(in) , dimension(7) :: ap
!
      real(rkx) , dimension(7,2) , save :: apl
      real(rkx) , dimension(2) , save :: fal , fl , glatl , gll , secl , stll
      integer(ik4) , dimension(2) , save :: iydl
      real(rkx) , dimension(25,2) , save :: swcl , swl
      integer(ik4) :: i
      logical :: ldc
!
!     Test if geophysical variables or switches changed and save
!     Return 0 if unchanged and 1 if changed
!
      data iydl  /2* -999/
      data secl  /2* -999.0_rkx/
      data glatl /2* -999.0_rkx/
      data gll   /2* -999.0_rkx/
      data stll  /2* -999.0_rkx/
      data fal   /2* -999.0_rkx/
      data fl    /2* -999.0_rkx/
      data apl   /14* -999.0_rkx/
      data swl   /50* -999.0_rkx/
      data swcl  /50* -999.0_rkx/
      vtst7 = d_zero
      ldc = .false.

      if ( iyd /= iydl(ic) ) then
        ldc = .true.
      else
        if ( abs(sec-secl(ic)) > nearzero ) then
          ldc = .true.
        else
          if ( abs(glat-glatl(ic)) > nearzero ) then
            ldc = .true.
          else
            if ( abs(glong-gll(ic)) > nearzero ) then
              ldc = .true.
            else
              if ( abs(stl-stll(ic)) > nearzero ) then
                ldc = .true.
              else
                if ( abs(f107a-fal(ic)) > nearzero ) then
                  ldc = .true.
                else
                  if ( abs(f107-fl(ic)) > nearzero ) then
                    ldc = .true.
                  else
                    do i = 1 , 7
                      if ( abs(ap(i)-apl(i,ic)) > nearzero ) then
                        ldc = .true.
                        exit
                      end if
                    end do
                    do i = 1 , 25
                      if ( abs(sw(i)-swl(i,ic)) > nearzero ) then
                        ldc = .true.
                        exit
                      end if
                      if ( abs(swc(i)-swcl(i,ic)) > nearzero ) then
                        ldc = .true.
                        exit
                      end if
                    end do
                  end if
                end if
              end if
            end if
          end if
        end if
      end if

      if ( ldc ) then
        vtst7 = d_one
        iydl(ic) = iyd
        secl(ic) = sec
        glatl(ic) = glat
        gll(ic) = glong
        stll(ic) = stl
        fal(ic) = f107a
        fl(ic) = f107
        do i = 1 , 7
          apl(i,ic) = ap(i)
        end do
        do i = 1 , 25
          swl(i,ic) = sw(i)
          swcl(i,ic) = swc(i)
        end do
      end if

    end function vtst7
!
end module physics_msis
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

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
 
module mod_che_drydep
!
! Chemical and aerosol surface emission and dry deposition
!
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_che_common
  use mod_che_dust
  use mod_mpmessage
  use mod_service
  use mod_che_ncio
  use mod_che_mppio
  use mod_che_indices


      private
!
      public :: drydep_aero, drydep_gas, a1,a2,a3,c1,c2,c3,c4,aa1,aa2,aa3
!
!     Dynamic Viscosity Parameters
!
      real(8) , parameter :: a1 = 145.8
      real(8) , parameter :: a2 = 1.5
      real(8) , parameter :: a3 = 110.4
!
!     Molecular Free Path calculation parameters
!
      real(8) , parameter :: c1 = 6.54E-8
      real(8) , parameter :: c2 = 1.818E-5
      real(8) , parameter :: c3 = 1.013E5
      real(8) , parameter :: c4 = 293.15
!
!     Cunningham slip correction factor parameters
!
      real(8) , parameter :: aa1 = 1.257
      real(8) , parameter :: aa2 = 0.4
      real(8) , parameter :: aa3 = 1.1


!
! Only one cover type per grid cell for now 

      integer, parameter :: luc = 1
!    
      integer, parameter :: ngasd = 31  ! number of gas taken into account by drydep scheme

      real(8), parameter :: rainthr = 0.1 ! threshold of rainfall intensity to activate water covered canopy option



!DATA section for the Zhang drydep scheme
!
!the Zhang scheme uses its own LAI . first index is consistent with BATS LU type, 15 is for the number of month 12 + 3 for interp.
!the tables are supposed to be consistant with BATS types determined by ivegcov (be careffull with ethe ocean option) 
!Rq: It would be better to have intercative directly :LAI and roughness from BATS and CLM 
!Rq2: Since ivegcov is defined even when clm is activated, the dry dep schem could in principle be used with CLM . 
!BUT , there is also the option of activating CLM PFT level dydep scheme.  


       real(8) lai(20,15)

       real(8) z01(20),z02(20)

       data (lai(1,kk), kk = 1, 15)/            &
           0.1  , 0.1  , 0.1  , 0.5  , 1.0  ,    &
           2.0  , 3.0  , 3.5  , 4.0  , 0.1  , &
           0.1  , 0.1  , 0.1  , 0.1  , 4.0  /

      data (lai(2,kk), kk = 1, 15)/ &
           1.0  , 1.0  , 1.0  , 1.0  , 1.0  , &
           1.0  , 1.0  , 1.0  , 1.0  , 1.0  , &
           1.0  , 1.0  , 1.0  , 1.0  , 1.0  /

      data (lai(3,kk), kk = 1, 15)/ &
           5.0  , 5.0  , 5.0  , 5.0  , 5.0  , &
           5.0  , 5.0  , 5.0  , 5.0  , 5.0  , &
           5.0  , 5.0  , 5.0  , 5.0  , 5.0  /

       data (lai(4,kk), kk = 1, 15)/         &
            0.1  , 0.1  , 0.5  , 1.0  , 2.0  , &
            4.0  , 5.0  , 5.0  , 4.0  , 2.0  , &
            1.0  , 0.1  , 0.1  , 0.1  , 5.0  /

       data (lai(5,kk), kk = 1, 15)/ &
            0.1  , 0.1  , 0.5  , 1.0  , 2.0  , &
            4.0  , 5.0  , 5.0  , 4.0  , 2.0  , &
            1.0  , 0.1  , 0.1  , 0.1  , 5.0  /

       data (lai(6,kk), kk = 1, 15)/ &
           6.0  , 6.0  , 6.0  , 6.0  , 6.0  , &
           6.0  , 6.0  , 6.0  , 6.0  , 6.0  , &
           6.0  , 6.0  , 6.0  , 6.0  , 6.0  /

       data (lai(7,kk), kk = 1, 15)/ &
            0.5  , 0.5  , 0.5  , 0.5  , 0.5  , &
            0.5  , 1.0  , 2.0  , 2.0  , 1.5  , &
            1.0  , 1.0  , 0.5  , 0.5  , 2.0  /
 
       data (lai(8,kk), kk = 1, 15)/ &
           0.0  , 0.0  , 0.0  , 0.0  , 0.0  , &
           0.0  , 0.0  , 0.0  , 0.0  , 0.0  , &
           0.0  , 0.0  , 0.0  , 0.0  , 0.0  /
    

       data (lai(9,kk), kk = 1, 15)/ &
            1.0  , 1.0  , 0.5  , 0.1  , 0.1  , &
            0.1  , 0.1  , 1.0  , 2.0  , 1.5  , &
            1.5  , 1.0  , 1.0  , 0.1  , 2.0  /

       data (lai(10,kk), kk = 1, 15)/ &
           1.0  , 1.0  , 1.0  , 1.0  , 1.0  , &
           1.0  , 1.0  , 1.0  , 1.0  , 1.0  , &
           1.0  , 1.0  , 1.0  , 1.0  , 1.0  /


       data (lai(11,kk), kk = 1, 15)/ &
           0.0  , 0.0  , 0.0  , 0.0  , 0.0  , &
           0.0  , 0.0  , 0.0  , 0.0  , 0.0  , &
           0.0  , 0.0  , 0.0  , 0.0  , 0.0  /

       data (lai(12,kk), kk = 1, 15)/ &
           0.0  , 0.0  , 0.0  , 0.0  , 0.0  , &
           0.0  , 0.0  , 0.0  , 0.0  , 0.0  , &
           0.0  , 0.0  , 0.0  , 0.0  , 0.0  /
   
       data (lai(13,kk), kk = 1, 15)/ &
           4.0  , 4.0  , 4.0  , 4.0  , 4.0  , &
           4.0  , 4.0  , 4.0  , 4.0  , 4.0  , &
           4.0  , 4.0  , 4.0  , 4.0  , 4.0  /

       data (lai(14,kk), kk = 1, 15)/ &
           0.0  , 0.0  , 0.0  , 0.0  , 0.0  , &
           0.0  , 0.0  , 0.0  , 0.0  , 0.0  , &
           0.0  , 0.0  , 0.0  , 0.0  , 0.0  /

       data (lai(15,kk), kk = 1, 15)/ &
           0.0  , 0.0  , 0.0  , 0.0  , 0.0  , &
           0.0  , 0.0  , 0.0  , 0.0  , 0.0  , &
           0.0  , 0.0  , 0.0  , 0.0  , 0.0  /

       data (lai(16,kk), kk = 1, 15)/ &
           3.0  , 3.0  , 3.0  , 3.0  , 3.0  , &
           3.0  , 3.0  , 3.0  , 3.0  , 3.0  , &
           3.0  , 3.0  , 3.0  , 3.0  , 3.0  /

       data (lai(17,kk), kk = 1, 15)/ &
           3.0  , 3.0  , 3.0  , 3.0  , 3.0  , &
           3.0  , 3.0  , 3.0  , 3.0  , 3.0  , &
           3.0  , 3.0  , 3.0  , 3.0  , 3.0  /


      data (lai(18,kk), kk = 1, 15)/ &
           3.0  , 3.0  , 3.0  , 4.0  , 4.5  ,&
           5.0  , 5.0  , 5.0  , 4.0  , 3.0  ,&
           3.0  , 3.0  , 3.0  , 3.0  , 5.0  /

      data (lai(19,kk), kk = 1, 15)/&
           3.0  , 3.0  , 3.0  , 4.0  , 4.5  ,&
           5.0  , 5.0  , 5.0  , 4.0  , 3.0  ,&
           3.0  , 3.0  , 3.0  , 3.0  , 5.0  /

      data (lai(20,kk), kk = 1, 15)/&
           0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,&
           0.0  , 0.0  , 0.0  , 0.0  , 0.0  ,&
           0.0  , 0.0  , 0.0  , 0.0  , 0.0  /

      data  z01/0.02, 0.04, 0.90, 0.40, 0.40, 2.00,&
               0.02, 0.04, 0.03, 0.05, 0.04, 0.01,&
               0.10, 0.00, 0.00, 0.20, 0.20, 0.60,&
               0.60, 0.00 /

      data  z02/0.10, 0.04, 0.90, 0.90, 1.00, 2.00,&
               0.10, 0.04, 0.03, 0.05, 0.04, 0.01,&
               0.10, 0.00, 0.00, 0.20, 0.20, 0.90,&
               0.90, 0.00 /
       

! Zhang stomatal resistance parameters

       real(8)::tmin(20),tmax(20)
       real(8)::rsminz(20),brs(20)
       real(8)::topt(20),bvpd(20)
       real(8)::psi1(20),psi2(20)
       real(8)::rac1(20),rac2(20)
       real(8)::rgo(20),rcutdO(20)
       real(8)::rcutwO(20),rcutdS(20)
       real(8)::rgs(20),sdmax(20)  
       real(8)::mw(31),rm(31)
       real(8)::alphaz(31),betaz(31) 


       data tmin/5.0, 5.0,  -5.0, -5.0, 0.0,    &
                0.0, 5.0,-999.0, -5.0, 5.0,     &
             -999.0,-999.0, 0.0,-999.0,-999.0,  &
                0.0,   0.0,-3.0,   0.0,-999.0 /

       data tmax/45.0, 40.0,  40.0, 40.0, 45.0, &
                45.0, 45.0,-999.0, 40.0, 45.0, & 
              -999.0,-999.0, 45.0,-999.0,-999.0,&
                45.0, 45.0,  42.0,  45.0,-999.0 /

       data rsminz/120.0, 150.0, 250.0, 250.0, 150.0,&
                 150.0, 100.0,-999.0, 150.0, 150.0, &
                -999.0,-999.0, 150.0,-999.0,-999.0, &
                 150.0, 250.0, 150.0, 150.0,-999.0 /

       data brs/40.0,  50.0,  44.0, 44.0, 43.0,   &
               40.0,  20.0,-999.0, 25.0, 40.0,    &
             -999.0,-999.0, 40.0,-999.0,-999.0,   &
               40.0,  44.0, 44.0,  43.0,-999.0 /

       data topt/27.0, 30.0,  15.0, 15.0, 27.0,  &
                 30.0, 25.0,-999.0, 20.0, 25.0,  &
              -999.0,-999.0, 20.0,-999.0,-999.0, &
                30.0, 25.0,  21.0,  25.0,-999.0 /

      data bvpd/0.00, 0.0,  0.31, 0.31, 0.36,    &
                0.27, 0.0,-999.0, 0.24,  0.0,    &    
              -999.0,-999.0,0.27,-999.0,-999.0,  &
                0.27, 0.27, 0.34, 0.31,-999.0 /

       data psi1/ -1.5,  -1.5,  -2.0,  -2.0,  -1.9, &
                 -1.0,  -1.5,-999.0,   0.0,  -1.5, &
               -999.0,-999.0,  -1.5,-999.0,-999.0, &
                 -2.0,  -2.0,  -2.0,  -2.0,-999.0 /

       data psi2/-2.5,  -2.5,  -2.5,  -2.5,  -2.5, &
                -5.0,  -2.5,-999.0,  -1.5,  -2.5, &
              -999.0,-999.0,  -2.5,-999.0,-999.0, &
                -4.0, -3.5,   -2.5,  -3.0,-999.0 /

      data rac1/10.00, 20.0, 100.0, 60.00, 40.0,&
               250.0, 10.0, 0.000, 40.00, 20.0, &
               0.000, 0.00, 20.00, 00.00, 0.00, &
               60.00, 40.0, 100.0, 100.0, 0.00 /

      data rac2/40.00, 20.0, 100.0, 100.0, 40.0, &
               250.0, 40.0, 0.000, 40.00, 20.0, &
               0.000, 0.00, 20.00, 0.000, 0.00, &
               60.00, 40.0, 100.0, 100.0, 0.00 /

       data rcutdO /4000.0, 4000.0, 4000.0, 4000.0, 6000.0, &
                   6000.0, 4000.0, -999.0, 8000.0, 4000.0,  &
                   -999.0, -999.0, 5000.0, -999.0, -999.0,&
                   6000.0, 5000.0, 4000.0, 4000.0, -999.0 /

       data rcutwO /200.00, 200.00, 200.00, 200.00, 400.00, &
                   400.00, 200.00, -999.0, 400.00, 200.00, &
                   -999.0, -999.0, 300.00, -999.0, -999.0, &
                   400.00, 300.00, 200.00, 200.00, -999.0 /

       data rgO    /200.00, 200.00, 200.00, 200.00, 200.00, &
                   200.00, 200.00, 500.00, 500.00, 500.00, &
                   500.00, 2000.0, 500.00, 2000.0, 2000.0, &
                   200.00, 200.00, 200.00, 200.00, 2000.0 /


       data rcutds /1500.0, 1000.0, 2000.0, 2000.0, 2500.0, &
                   2500.0, 1000.0, -999.0, 2000.0, 2000.0, &
                   -999.0, -999.0, 1500.0, -999.0, -999.0, &
                   2000.0, 2000.0, 2500.0, 2500.0, -999.0 /

       data rgs    /200.0, 200.0, 200.0, 200.0, 200.0, &
                   100.0, 200.0, 700.0, 300.0,  50.0, &
                   700.0,  70.0,  50.0,  20.0,  20.0, &
                   200.0, 200.0, 200.0, 200.0,  20.0 /

       data sdmax/10.0,  5.0, 200.0,   1.1, 200.0, &
                400.0, 20.0,   2.0,   2.0,  10.0, &
                  2.0,  1.0,  10.0,-999.0,-999.0, &
                 50.0, 50.0, 200.0, 200.0,-999.0 /
    
!   *****************************************************************
!   * Gas Properties (Total 31 species)                          ****
!   * Mesophyll resistance RM, scaling factors ALPHAZ and BETAZ,   ****
!   * molecular weight.                                          ****
!   *****************************************************************
! parameters are given for the folloewing species in the zhang scheme
! SO2    H2SO4   NO2     O3     H2O2      HONO           HNO4            NH3 
!                 PAN           PPN       APAN           MPAN   
!                 HCHO          MCHO      PALD           C4A  
!                 C7A           ACHO      MVK            MACR 
!                 MGLY          MOH       ETOH           POH 
!                 CRES          FORM      ACAC           ROOH   
!                 ONIT          INIT 

       data   rm / 0.   , 0.   , 0.   ,  0.  ,  0.  , &
                 0.   , 0.   , 0.   ,  0.  ,  0.  ,  &
                 0.   , 0.   , 0.   ,  0.  ,  100.,&
                 100. , 100. , 100. ,  100.,  0.  , &
                 100. , 0.   , 0.   ,  0.  ,  0.  , &
                 0.   , 0.   , 0.   ,  0.  ,  100., &
                 100.     /

       data  alphaz /  1.   , 1.   , 0.   ,  0.  ,  1.  , &
                     10.  , 2.   , 5.   ,  1.  ,  0.  , &
                      0.   , 0.   , 0.   ,  0.8 ,  0.  ,  &
                      0.   , 0.   , 0.   ,  0.  ,  0.  , &
                      0.   , 0.01 , 0.6  ,  0.6 ,  0.4 ,  &
                      0.01 , 2.   , 1.5  ,  0.1 ,  0.  , &
                      0.       /
 
       data   betaz /  0.   , 1.   , 0.8  ,  1.  ,  1.  , &
                    10.  , 2.   , 5.   ,  0.0 ,  0.6 , &
                    0.6  , 0.8  , 0.3  ,  0.2 ,  0.05, &
                    0.05 , 0.05 , 0.05 ,  0.05,  0.05, &
                    0.05 , 0.   , 0.1  ,  0.  ,  0.  , &
                    0.   , 0.   , 0.   ,  0.8 ,  0.5 , &
                    0.5      /


       data  mw / 64.0  , 98.0  , 46.0  ,  48.0 ,  34.0 , &
                 63.0  , 47.0  , 79.0  ,  17.0 ,  121.0, &
                 135.0 , 183.0 , 147.0 ,  30.0 ,  44.0 , &
                 58.0  , 72.0  , 128.0 ,  106.0,  70.0 , &
                 70.0  , 72.0  , 32.0  ,  46.0 ,  60.0 , &
                 104.0 , 46.0  , 60.0  ,  48.0 ,  77.0 , &
                 147.0     /
              
     contains

      subroutine drydep_aero (j,mbin,indsp,rhop,ivegcov,throw,      &
                        & roarow,shj,pressg,temp2,sutemp,srad,rh10,     &
                        & wind10,zeff,trsize,pdepv,ddepv)
!
      implicit none
!
      integer :: j
      integer :: mbin
      integer , dimension(mbin) :: indsp
      integer , dimension(iy) :: ivegcov

      real(8) , dimension(iy) :: pressg , rh10 , srad , sutemp ,       &
                                & temp2 , wind10 , zeff
      real(8) , dimension(iy,kz) :: roarow , throw
      real(8) , dimension(kz) :: shj
      real(8) , dimension(mbin,2) :: trsize
      real(8)  :: rhop 

! output table to be passed out. Care dimension is ntr  

      real(8) , dimension(iy,kz,ntr) :: pdepv
      real(8) , dimension(iy,ntr) :: ddepv

      intent (in) j, indsp, ivegcov , mbin ,       &
                & pressg , rh10 , roarow , shj , srad , sutemp , temp2 ,&
                & throw , trsize , wind10 , zeff, rhop
      intent (inout) pdepv
!
      real(8) :: amfp , amob , asq , ch , cm , cun , dtemp , dthv , eb ,&
               & eim , ein , es , fh , fm , frx1 , kui , logratio ,     &
               & mol , pre , prii , priiv , psit , psiu , ptemp2 , qs , &
               & r1 , ratioz , aa , rib , st , tbar , thstar , tsv ,    &
               & tsw , ustarsq , utstar , vp , vptemp , wvpm , ww , x , &
               & y , z , z0water , zdl , zl
      real(8) , dimension(iy,kz) :: amu
      real(8) , dimension(iy) :: anu , schm , zz0
      real(8) , dimension(iy,kz,isize) :: cfac , pdepvsub ,   pdiff ,  &
           & rhsize , taurel

      integer :: i , jc , k , kcov , l , lev , n , tot,ib
      real(8) , dimension(iy,luc) :: ra , ustar , vegcover
      real(8) , dimension(iy,luc,isize) :: rs
       real(8), dimension(iy,kz) :: wk, settend
      real(8) , parameter :: z10 = 10.0
      real(8) , dimension(isize) :: avesize
      character (len=64) :: subroutine_name='chdrydep'
      integer :: idindx=0
!
      call time_begin(subroutine_name,idindx)

    
      do n = 1 , isize
        avesize(n) = (aerosize(1,n)+aerosize(2,n))/2.0
      end do
      
!======================================================================
!     ********************************************************
!     *   aerosize - dry radius                           ****
!     *   rhop  - density for each aerosol type           ****
!     ********************************************************
      do n = 1 , isize
        do l = 1 , kz
          do i = 2 , iym2
 
!           ********************************************************
!           *  aerosol gravitational settling velocity          ****
!           *  and diffusion coefficient                        ****
!           *                                                   ****
!           * air's dynamic viscosity                           ****
!           ********************************************************
 
            amu(i,l) = a1*1.E-8*throw(i,l)**a2/(throw(i,l)+a3)
 
!           . . . . mid layer pressure in [pascal].
            pre = pressg(i)*shj(l)
!           ********************************************************
!           * mean molecular free path.                         ****
!           *     k.v. beard [1976], j atm. sci., 33            ****
!           ********************************************************
 
            amfp = c1*(amu(i,l)/c2)*(c3/pre)*(throw(i,l)/c4)**(1./2.)
            prii = 2./9.*egrav/amu(i,l)
            priiv = prii*(rhop-roarow(i,l))
 
!           ********************************************************
!           * cunningham slip correction factor and             ****
!           * relaxation time = vg/grav.                        ****
!           ********************************************************
 
            cfac(i,l,n) = 1. + amfp/avesize(n)                          &
                        & *(aa1+aa2*exp(-aa3*avesize(n)/amfp))
            taurel(i,l,n) = dmax1(priiv*avesize(n)**2*cfac(i,l,n)*regrav, &
                          & 0.D0)
 
!           ********************************************************
!           * stokes friction                                  *****
!           pdepvsub(i,l,n) ' sellting dep. velocity = '
!           ********************************************************
 
            pdepvsub(i,l,n) = taurel(i,l,n)*egrav
          end do
        end do
      end do

! find aerodynamic resistance

       call aerodyresis(zeff,wind10,temp2      &
     &   , sutemp,rh10,srad,ivegcov,ustar,ra )


 
!     *****************************************************
!     * the schmidt number is the ratio of the         ****
!     * kinematic viscosity of air to the particle     ****
!     * brownian diffusivity ===> sc=v/d               ****
!     *****************************************************
 
      do n = 1 , isize
        do l = 1 , kz
          do i = 2 , iym2
 
!           *****************************************************
!           * for now we will not consider the humidity      ****
!           * impact so we will set the variable frx1=1.0    ****
!           * i.e. only dry particles                        ****
!           *****************************************************
 
            frx1 = 1.0
            rhsize(i,l,n) = avesize(n)*frx1
            anu(i) = amu(i,l)/roarow(i,l)
            amob = 6.*mathpi*amu(i,l)*rhsize(i,l,n)/cfac(i,l,n)
            pdiff(i,l,n) = boltzk*throw(i,l)/amob
            schm(i) = anu(i)/pdiff(i,l,n)
 
!           ******************************************************
!           * for brownian diffusion, there is evidence that  ****
!           * its fromula depend on schmidt number as :       ****
!           * eb= schm x c^gama                               ****
!           * where gama is efficiency factor and its value   ****
!           * between 1/2 and 2/3 with larger values for      ****
!           * rougher surfaces                                ****
!           * ****************************************************
 
          end do
          if ( l.eq.kz ) then
            do k = 1 , luc ! luc  = 1 for the moment
              do i = 2 , iym2

!======================================================================
!     find the right table index for the cell cover ( ocean and lake
!     are 0 in the ivegcov and 14-15 in the table )

              if(ivegcov(i)==0) then
                  kcov = 14
              else if (ivegcov(i)>20) then
                  kcov = 20
              else
                  kcov = ivegcov(i)
              end if

 
!               ******************************************************
!               * the parameter governing impaction processes is *****
!               * the stokes number,st, which has the form of    *****
!               * 1) st = vg x u* /g a for vegetated surefaces   *****
!               *    (slinn, 1982)                               *****
!               * 2) st = vg x u*2/anu for smothed surfaces or   *****
!               *    surfaces with bluff roughness elements      *****
!               *    (giorgi,1988)                               *****
!               ******************************************************
 
                st = taurel(i,l,n)*ustar(i,k)*ustar(i,k)/anu(i)
 
                eb = schm(i)**(-0.666667)
!               eim=(st/(st+aest(k)))**2
                eim = (st/(st+aest(kcov)))**2
 
 
                eim = dmin1(eim,0.6D0)
                ein = 0.0
!               if (arye(k) .gt. 0.0001) then
!               ein = (1000.0*2.0*avesize(n)/arye(k))**1.5
!               end if
 
                if ( arye(kcov).gt.0.0001 )                             &
                   & ein = (1000.0*2.0*avesize(n)/arye(kcov))**1.5
 
                ein = dmin1(ein,0.5D0)
 
!               *****************************************************
!               * partickes larger than 5 micro may rebounded   *****
!               * after hitting a surface,this process may be   *****
!               * included by modifying the total collection    *****
!               * by the factor of r1, which represents the     *****
!               * fraction of particles sticking to the surface *****
!               * slinn (1982) suggested the following:         *****
!               * r = exp (- st^0.2)                            *****
!               *****************************************************
 
!               r1 = exp (-st**0.5)
                r1 = dmax1(0.5D0,exp(-st**0.5))
!               if (k .ge. 11 .and. r1 .lt. 0.5 ) r1=0.5
                if ( kcov.ge.11 .and. r1.lt.0.5 ) r1 = 0.5
                if ( r1.lt.0.4 ) r1 = 0.4
 
!               ***************************************************
!               * calculation of rs: the surface resistance   *****
!               * which depends on the collection efficiency  *****
!               * of the surface and is determined by the     *****
!               * various deposition processes                *****
!               ***************************************************
 
!               rs= 1.0/ustar(i,k)/(eb+eim+ein)/r1
                rs(i,k,n) = 1.0/3.0/ustar(i,k)/(eb+eim+ein)/r1
              end do
            end do
          end if
        end do
      end do
!======================================================================
 

!     average settling and deposition velocities on bin
! care we use pdepv and ddpv table that are dimensionned to ntr and not mbin ! 
!
  
      do ib = 1 , mbin
        tot = 0
         pdepv(:,:,indsp(ib)) = 0.
         ddepv(:,indsp(ib))  =0.
        do n = 1 , isize
          if ( avesize(n)*1.E6.ge.trsize(ib,1) .and. avesize(n)          &
             & *1.E6.lt.trsize(ib,2) ) then
 
              do i = 2 , iym2
                pdepv(i,:,indsp(ib)) = pdepv(i,:,indsp(ib)) + pdepvsub(i,:,n)
! agregate the dry deposition velocity, remember one cover per grid cell for now
! the dry deposition deposition velocity must accound also for the settling vrlocity at kz
! simple form now add the vs
                ddepv(i,indsp(ib)) =  ddepv(i,indsp(ib)) + ( 1.0/(ra(i,1)+rs(i,1,n)) +  pdepvsub(i,kz,n))

              end do          
            tot = tot + 1
          end if
        end do
        if ( tot.gt.0 ) then
        
            do i = 2 , iym2
              pdepv(i,:,indsp(ib)) = pdepv(i,:,indsp(ib))/tot
              ddepv(i,indsp(ib)) = ddepv(i,indsp(ib))/tot 
           end do       
        end if
      end do
! 
! Finally update the emission and settling tendencies for dust and sea salt 
!

          do ib = 1,mbin
!         deposition
          do k = 2 , kz
            do i = 2 , iym2
              wk(i,k) = (1./cpsb(i,j))*          &
                      & (ctwt(k,1)*chib(i,k,j,indsp(ib))+ &
                      &  ctwt(k,2)*chib(i,k-1,j,indsp(ib)))
            end do
          end do

          do i = 2 , iym2
            do k = 2 , kz - 1
                        ! do not apply to the first level
              settend(i,k) = (wk(i,k+1)*pdepv(i,k+1,indsp(ib)) - &
                              wk(i,k)*pdepv(i,k,indsp(ib)))*     &
                              egrav*1.E-3/cdsigma(k)

              chiten(i,k,j,indsp(ib)) = chiten(i,k,j,indsp(ib)) - settend(i,k)
            end do
!
! at first level include surface drydep velocity to calculte the divergence
            settend(i,kz) =   (chib(i,kz,j,indsp(ib))/ cpsb(i,j) * ddepv(i,indsp(ib))   &  
                              -wk(i,kz)*pdepv(i,kz,indsp(ib)) ) &
                              *egrav*1.E-3/cdsigma(kz)

            chiten(i,kz,j,indsp(ib)) = chiten(i,kz,j,indsp(ib)) - settend(i,kz)
 
!           dignoctic for dry deposition
            remdrd(i,j,indsp(ib) ) = remdrd(i,j,indsp(ib)) + settend(i,kz)*dtche/2.

!dry dep velocity diagnostic in m.s-1  (sum setlling  + drydep v. , accumulated between two outputs time step) 
            drydepv(i,j,indsp(ib)) =  drydepv(i,j,indsp(ib)) +  ddepv(i,indsp(ib)) + pdepv(i,kz,indsp(ib))

          end do
          end do





      call time_end(subroutine_name,idindx)
      end subroutine drydep_aero
!
! 
!*************************************************************************************************************
!*************************************************************************************************************
      subroutine drydep_gas (j, ivegcov ,       &
                 rh10, srad , tsurf , prec, temp10 ,  &
                 wind10 , zeff, drydepvg)


#ifdef CLM
       use clm_drydep, only : c2rvdep
#endif

       use mod_che_indices


      implicit none
!

      integer :: j   
      integer , dimension(iy) :: ivegcov
   
      real(8) , dimension(iy) :: rh10 , srad , tsurf ,  prec,     &
                                & temp10 , wind10 , zeff

      real(8), dimension(iy,ntr) :: drydepvg

      intent (in) j, ivegcov ,     &
                 rh10, srad , tsurf , temp10 , &
                 wind10 
      intent(inout) zeff
      intent(out) drydepvg


      integer :: n,i,k,im,iday_m,kcov
      real (8), dimension(iy,luc) :: ustar, resa
      real (8), dimension(ngasd, iy, luc) :: resb, resc
      real (8), dimension(ngasd, iy, luc) :: vdg
      real(8), dimension(iy) :: ddrem
      real(8) , dimension(iy)::lai_f,laimin,laimax, snow
      real(8) ::Kd

#ifndef CLM

! Different options for LAI and roughness 
! for the moment read from 
   
        im = INT(ccalday / 30.5 ) + 1
        iday_m =ccalday - INT((im-1)*30.5+0.5)

      if (iday_m .eq. 0) THEN
           im = im - 1
           iday_m =ccalday - (im - 1)*30.5
      end if 
     
       do i = 2, iym2
 
      if(ivegcov(i)==0) then
           kcov = 14
       else if (ivegcov(i)>20) then
           kcov = 20
       else
           kcov = ivegcov(i)
       end if
       
       lai_f(i)  = lai(kcov ,im) &
                    + iday_m / 30.5*(lai(kcov,im + 1 )-lai(kcov,im))    
       
       laimin(i) = lai(kcov,14)
       laimax(i) = lai(kcov,15)

       end do 


       call aerodyresis(zeff,wind10,temp10      &
     &   , tsurf,rh10,srad,ivegcov,ustar,resa )


 !      call zenitm(coszrs,iy,j)

       snow(:) =0. 


       call stomtresis(lai_f,laimin,laimax, ivegcov,ngasd,    &
            ustar,prec,snow,srad,                               &  
            tsurf,temp10,rh10,czen(j,:),resc,resb)


! now calculate the dry deposition velocities and select it according to the gasphase mechanism
!vdg in m.s-1
       vdg(:,:,:)=0.      
       do i=2,iym2
       do n=1,ngasd
       vdg(n,i,:) = 1./ (resa(i,:)+ resb(n,i,:)+ resc (n,i,:))
       end do
       end do
! this part depends on the chem mechanism
! for CBMZ , we can certainly improve this.


       drydepvg = 0.

       drydepvg(:,iso2)  =  vdg(1,:,1)
       drydepvg(:,iso4)  =  vdg(2,:,1)
       drydepvg(:,ino2)  =  vdg(3,:,1)!*0.5
       drydepvg(:,io3)   =  vdg(4,:,1)!*0.5
       drydepvg(:,ih2o2) =  vdg(5,:,1)!*0.5
       drydepvg(:,ihno3) =  vdg(6,:,1)!*0.5
       drydepvg(:,ipan)  =  vdg(10,:,1)!*0.5 
       drydepvg(:,ihcho) =  vdg(14,:,1)!*0.5
       drydepvg(:,iald2) =  vdg(15,:,1)!*0.5
       drydepvg(:,imoh)  =  vdg(23,:,1)!*0.5





! Finally : gas phase dry dep tendency calculation 

          do i=2,iym2
          do n = 1, ntr
            Kd =  drydepvg(i,n) / cdzq(j,i,kz) !Kd removal rate in s-1
             ddrem(i)  =  chib(i,kz,j,n) * (1 -   dexp(-Kd*dtche)) / dtche ! dry dep removal tendency (+)
!update chiten
             chiten(i,kz,j,n) = chiten(i,kz,j,n) - ddrem(i)
!drydep flux diagnostic (accumulated between two outputs time step) 
             remdrd(i,j,n) = remdrd(i,j,n) + ddrem(i) * dtche / 2 
!dry dep velocity diagnostic in m.s-1  (accumulated between two outputs time step) 
             drydepv(i,j,n) =  drydepv(i,j,n) + drydepvg(i,n) 
          end do
          end do


#endif


! if CLM is used then use directly the clm dry dep module.
#ifdef CLM

       jj = j + (jxp*myid)

          do i=2,iym2
          do n = 1, ntr

! use clm dry deposition velocity
            Kd =  c2rvdep(jj,i,itr) / dzq(j,i,kz) !Kd removal rate in s-1
             ddrem(i)  =  chib(i,kz,j,n) * (1 -   dexp(-Kd*dtche)) / dtche ! dry dep removal tendency (+)
!update chiten
             chiten(i,kz,j,n) = chiten(i,kz,j,n) - ddrem(i)
!drydep flux diagnostic (accumulated between two outputs time step) 
             remdrd(i,j,n) = remdrd(i,j,n) + ddrem(i) * dtche / 2 
!
!          drydepv(i,j,n) =  drydepv(i,j,n) + drydepvg(i,n) 
          end do
          end do
          end if
#endif



      
      end subroutine drydep_gas

!
!%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


       subroutine aerodyresis(zeff,wind10,temp2,         &
                             sutemp,rh10,srad,ivegcov,ustar,ra)

       implicit none
     
       
         
       integer ::ivegcov(iy)
       
       real(8) :: temp2(iy),wind10(iy),rh10(iy)
       real(8) :: sutemp(iy),srad(iy)
       real(8) :: zeff(iy)
       
       intent (in) zeff,wind10,temp2,         &
                   sutemp,rh10,srad,ivegcov

       real(8) :: ustar(iy,luc),ra(iy,luc)
       intent (out) ustar, ra

! local variables
       integer i,j       
       real(8) :: delz,vp,tsv
       real(8) :: z,zl,ww
       real(8) :: ptemp2,es,qs
       real(8) :: wvpm,vptemp,tsw,mol
       real(8) :: z0water,dthv,cun,zdl
       real(8) :: psiu,psit,x,y
       real(8) :: thstar,rib,dtemp,tbar
       real(8) :: ustarsq,utstar,kui
       real(8) :: ratioz,logratio,asq
       real(8) :: aa,cm,ch,fm,fh    
       real(8) :: zz0(iy)
       real(8) , parameter :: z10 = 10.0
  
!======================================================================
!     ****************************************************
!     * ra : is the aerodynamic resistance above the  ****
!     *      canopy and it is function in u* and      ****
!     *      z0: roughness length and the stability   ****
!     *      function                                 ****
!     * mol  - monin obukhov length (m) - calculated  ****
!     *           for each land use category          ****
!     * ptemp2 -potential temperature at z2  (deg. k) ****
!     * temp2 - temperature at 10m. (deg k)           ****
!     * z10   - 10 m.                                 ****
!     * sutemp -surface temperature (deg k)           ****
!     * srad   -solar irradiance at the ground(w/m**2)****
!     * rh10  - relative humidity of air at 10m.      ****
!     *           (0.0-1.0)                           ****
!     * stdpmb - sea level pressure (mb)               ****
!     ****************************************************
      do j = 1 , luc
        do i = 2 , iym2
          ww = dmax1(wind10(i),1.0D0)
          zz0(i) = zeff(i)
 
! ***************************************************************
!         *  potential temperature at z2  (deg. k)                 
!         *****
 
! ***************************************************************
          ptemp2 = temp2(i) + z10*0.0098
 
! ***************************************************************
!         * for calculations over water compute values of critical 
!         ***** * profile variables: l and ustar                       
!         ***** *           ******begin for water***                   
!         *****
 
! ***************************************************************
          if ( ivegcov(i).eq.0 ) then
 
! **************************************************************
!           *  vp  - vapour pressure at z2                          
!           ***** *  wvpm- water vapour mixing ratio at  z2            
!           ***** *  vptemp- virtual potential temperature at z2 (deg.
!           k)  *****
! **************************************************************
            es = 6.108*exp(17.27*(temp2(i)-tzero)/(temp2(i)-35.86))
            vp = rh10(i)*es
            wvpm = ep2*vp/(stdpmb-vp)
            vptemp = ptemp2*(1.0+0.61*wvpm)
 
! **************************************************************
!           *  assume rh10 at water surface is 100%                 
!           ***** *   vp = es(tsw-tzero) !sat. vap press at surface   
!           ***** *   saturated vapour pressure at surface             
!           ***** *   saturated mixing ratio at surface                
!           ***** *   tsv - virtual potential temperature at surface   
!           ***** *           (deg. k)                                 
!           *****
 
! **************************************************************
            tsw = sutemp(i)
            vp = 6.108*exp(17.27*(tsw-tzero)/(tsw-35.86))
            qs = ep2*vp/(stdpmb-vp)
            tsv = tsw*(1.+0.61*qs)
            z0water = 1.0E-4
! **************************************************************
!           * scalet  :  not required if  z2 = 10m                  
!           *****
! **************************************************************
            dthv = (vptemp-tsv)
! **************************************************************
!           * calculate drag coefficient cun with neutral condition 
!           ***** * assumption  garratt (1977)                         
!           *****
 
! **************************************************************
            cun = 7.5E-4 + 6.7E-5*ww
            mol = 9999.0
 
            if ( abs(dthv).gt.1.0E-6 )                                  &
               & mol = vptemp*cun**1.5*ww**2/(5.096E-3*dthv)
            if ( mol.gt.0. .and. mol.lt.5.0 ) mol = 5.0
            if ( mol.gt.-5.0 .and. mol.lt.0 ) mol = -5.0
            zdl = z10/mol
!
            if ( zdl.lt.0.0 ) then
 
! **************************************************************
!             *                        wind speed                     
!             *****
 
! **************************************************************
              x = (1.0-15.0*zdl)**0.25
              psiu = 2.*dlog(0.5*(1.0+x)) + dlog(0.5*(1.0+x*x))         &
                   & - 2.0*atan(x) + 0.5*mathpi
 
! **************************************************************
!             *                       pot temp                        
!             *****
! **************************************************************
              y = sqrt(1.-9.*zdl)
              psit = 2.*0.74*dlog((1+y)/2.0)
            else
              psiu = -4.7*zdl
              psit = psiu
            end if
            z0water = 0.000002*ww**2.5
!
            ustar(i,j) = vonkar*ww/(dlog(z10/z0water)-psiu)
            thstar = vonkar*(ptemp2-sutemp(i))                          &
                   & /(0.74*dlog(z10/z0water)-psit)
            zz0(i) = z0water
!
          else
 
! **************************************************************
!           * compute ustar and l for land use categories other than 
!           **** * water use louis method. !pkk 7/16/85, find bulk     
!           **** * richardson number.                                  
!           ****
 
! **************************************************************
            rib = egrav*z10*(ptemp2-sutemp(i))/(sutemp(i)*ww**2)
 
! ***************************************************************
!           * ensure that conditions over land are never stable when 
!           ***** * there is incoming solar radiatiom                  
!           *****
! ***************************************************************
            if ( srad(i).gt.0.0 .and. rib.gt.0.0 ) rib = 1.E-15
!
            dtemp = ptemp2 - sutemp(i)
            if ( dabs(dtemp).lt.1.E-10 ) dtemp = dsign(1.D-10,dtemp)
            tbar = 0.5*(ptemp2+sutemp(i))
!
            ratioz = z10/zz0(i)
            logratio = dlog(ratioz)
            asq = 0.16/(logratio**2)
!
            if ( rib.le.0.0 ) then
              aa = asq*9.4*sqrt(ratioz)
              cm = 7.4*aa
              ch = 5.3*aa
              fm = 1. - (9.4*rib/(1.+cm*sqrt(abs(rib))))
              fh = 1. - (9.4*rib/(1.+ch*sqrt(abs(rib))))
            else
              fm = 1./((1.+4.7*rib)**2)
              fh = fm
            end if
!
            ustarsq = asq*ww**2*fm
            utstar = asq*ww*dtemp*fh/0.74
            ustar(i,j) = sqrt(ustarsq)
            thstar = utstar/ustar(i,j)
!
            mol = tbar*ustarsq/(vonkar*egrav*thstar)
          end if
 
          kui = 1.0/(vonkar*ustar(i,j))
 
!         **************************************************************
!         * compute the values of  ra                            *******
!         **************************************************************
 
          z = z10
          zl = z/mol
 
          if ( zl.ge.0. ) then
            ra(i,j) = kui*(0.74*dlog(z/zz0(i))+4.7*zl)
          else
            ra(i,j) = kui*0.74*(dlog(z/zz0(i))-2.0*dlog((1+sqrt(1-9.*zl)&
                    & )*0.5))
          end if
          ra(i,j) = dmax1(ra(i,j),0.99D0)
          ra(i,j) = dmin1(ra(i,j),999.9D0)
        end do
      end do

      end subroutine aerodyresis


!************************************************************************
!***********************************************************************

        subroutine stomtresis(lai_f,laimin,laimax, ivegcov,igas,       &
                            ustar,prec,sd,srad,ts,t2,rh,coszen,rc,rb )


       implicit none


       integer i,l,k,kk,j
       integer im,iday_m
       integer jday

       integer kcov,n,igas,ig
      
       integer ivegcov(iy),veg(iy)
       

       real(8), dimension (iy) :: lai_f , laimin, laimax, coszen, srad,  &
                                   ts,rh,prec, sd, t2

       real(8), dimension (iy,luc) ::  ustar

       real(8) :: rb(igas,iy,luc),rc(igas,iy,luc)

       real(8) :: rst,wst,rac,rgs_f
       real(8) :: rdu,rdv,rgo_f
       real(8) :: rcuto_f,rcuts_f 
       real(8) :: ww1,ww2,ww3
       real(8) :: rdm,rdn,rv,rn
       real(8) :: ratio,sv,fv,fvv
       real(8) :: pardir, pardif 
       real(8) :: tmaxk , tmink
       real(8) :: pshad,psun,rshad,rsun
       real(8) :: gshad,gsun,fsun,fshad
       real(8) :: gspar, temps !C
       real(8) :: bt,gt,gw,ryx
       real(8) :: es,d0,gd,psi
       real(8) :: coedew,dq,usmin
       real(8) :: fsnow,rsnows

       real(8) :: dgas,di,vi
       real(8) :: dair,dh2o,dvh2o,rstom
       real(8) :: rcut,rg

       logical :: is_dew, is_rain
       
            

       dair= 0.369 * 29.0 + 6.29
       dh2o=0.369*18.+6.29
       is_rain = .false.
       is_dew = .false.

       do j=1,luc
       do i=2,iym2

       if (ivegcov(i) .eq. 0) then
           kcov = 14
       else
           kcov = ivegcov(i)
       end if



!        print*,'srad ====', srad(i)
!        print*,' ts  ====', ts(i)
!        print*,' coszen == ', coszen(i)

        tmaxk = tmax(kcov) + 273.15
        tmink = tmin(kcov) + 273.15
!             print*,' tmax, tmin ==== ', tmaxk, tmink      

!FAB
! initialise rst as undef
         rst= -999.
        
       if (srad(i)      .ge. 0.1      .and.    &
          ts(i)        .lt. tmaxk    .and.    &
          ts(i)        .gt. tmink    .and.    &
          lai_f(i)  .gt. 0.001    .and.    &
          coszen(i)    .gt. 0.001                     ) then


      !================================================================
      ! Calculate direct and diffuse PAR from solar radiation and
      ! solar zenith angle
      !================================================================

         rdu  = 600.0 * dexp(-0.185/coszen(i))*coszen(i)
         rdv  = 0.4 * (600.0 - rdu ) * coszen(i)
         ww1   = -dlog(coszen(i))/2.302585

!        print*,' ww1 = ', ww1

         ww2   = -1.195 + 0.4459 * ww1 - 0.0345 * ww1**2
         ww3   = 1320*10**ww2

!         print*,'ww= ', ww
         rdm  = (720.*dexp(-0.06/coszen(i))-ww3)*coszen(i)
!         print*,'ww3= ', ww3, rdm

         rdn  = 0.6 * (720.0 - rdm - ww3) * coszen(i)
         rv   = dmax1(0.1d0,  rdu + rdv)
         rn   = dmax1(0.01d0, rdm + rdn)
         ratio= dmin1(0.9d0,srad(i)/( rv + rn))
!         print*,'ratio= ', ratio, rdn, rv, rn
         sv   = ratio * rv                            ! Total PAR
         fv   = dmin1(0.99d0, (0.901 - ratio)/0.7)
!         print*,'sv  fv  = ', sv, fv
!         print*,'rv  xxxxx  = ', rv, (1.0 - fv**0.6667)

         fvv  = dmax1(0.01d0,rdu/rv*(1.0 - fv**0.6667))
!         print*,'fvv  = ', fvv  
                             !fraction of PAR in the direct beam
         pardir = fvv * sv
                             ! PAR from direct radiation
         pardif = sv - pardir
                             ! PAR from diffuse radiation

!         print*,'pardif========', pardif   

      !================================================================
      ! Calculate sunlit and shaded leaf area, PAR for sunlit and
      ! shaded leaves
      !===============================================================

      if (lai_f(i) .gt. 2.5 .and. srad(i) .gt. 200.0) then

         pshad=pardif * dexp(-0.5 * lai_f(i)**0.8)    &
     &        + 0.07 * pardir * (1.1- 0.1*lai_f(i))* dexp(-coszen(i))
         psun=pardir**0.8*0.5/coszen(i) + pshad

      else
         pshad = pardif * dexp(-0.5 * lai_f(i)**0.7)       &
     &        + 0.07 * pardir *(1.1- 0.1*lai_f(i)) * dexp(-coszen(i))
         psun = pardir * 0.5/coszen(i) + pshad
      end if
         
!       print*,'pshad   psun   ', pshad , psun 

         rshad= rsminz(kcov) + brs(kcov) * rsminz(kcov)/pshad
         rsun = rsminz(kcov) + brs(kcov) * rsminz(kcov)/psun
         gshad= 1.0/rshad
         gsun = 1.0/rsun
!       print*,'rshad  ----< ', rshad, rsun, ' >---------rsun'
!       print*,'gshad  ----< ', gshad, gsun, ' >-------- gsun'

      !================================================================
      ! Fsun, Fshade are the total sunlit and shaded leaf area
      ! index
      !================================================================

         fsun =2.0 * coszen(i)*(1.0 - dexp(-0.5*lai_f(i)/coszen(i)))
                                   ! Sunlit leaf area
         fshad=lai_f(i) - fsun
                                   ! Shaded leaf area
!           print*,'f, f ====',fshad,fsun
      !================================================================
      ! Stomatal conductance before including effects of
      ! temperature, vapor pressure defict and water stress.
      !================================================================

         gspar = fsun * gsun + fshad * gshad

      !================================================================
      ! function for temperature effect
      !================================================================

         temps = ts(i) - 273.156
!FAB         bt= (tmax(kcov) - topt(kcov))/(tmax(kcov) - tmin(kcov))
         bt= (tmax(kcov) - topt(kcov))/(topt(kcov) - tmin(kcov))
         gt= (tmax(kcov) - temps)/(tmax(kcov) - topt(kcov))
         gt= gt**bt
         gt= gt*(temps - tmin(kcov))/(topt(kcov) - tmin(kcov))

!         print*,'gt ==========',gt

      !================================================================
      ! function for vapor pressure deficit
      !================================================================

         es = 6.108*exp(17.27*(ts(i) - 273.156)/(ts(i) - 35.86))
         d0=  es *(1.0 - rh(i))/10.0           !kPa
         gd= 1.0 - bvpd(kcov) * d0

!         print*,'gd===',gd
      !================================================================
      ! function for water stress
      !================================================================

         psi = (-0.72 - 0.0013 * srad(i))
!        psi_s=(-0.395-0.043*(ts-273.15))*102.
         gw=(psi - psi2(kcov))/(psi1(kcov) - psi2(kcov))

!         print*,'gw==',gw
!TEST
!         gw=1
         if (gw .gt. 1.0) gw=1.0
         if (gw .lt. 0.1) gw=0.1
         if (gd .gt. 1.0) gd=1.0
         if (gd .lt. 0.1) gd=0.1

      !================================================================
      ! Stomatal resistance for water vapor
      !================================================================

         rst=1.0/(gspar * gt * gd * gw)
!         print*,'rst===',rst

      end if


              
       coedew = 0.1  ! for clear cloud
       es = 6.108*exp(17.27*(ts(i) - 273.156)/(ts(i) - 35.86))
       dq = 0.622/1000.0 * es *(1.0- rh(i))*1000.0    ! unit g/kg
       dq = dmax1(0.0001d0,dq)
       usmin=1.5/dq*coedew

!        print*,'prec===== ', prec(i)
!       print*,'usmin   ===  ', usmin
! what is the unit of precipitation threshold 
        if(ts(i) .gt. 273.156 .and. prec(i) .gt. rainthr) then
        is_rain = .true.

!        print*, 'rain==='
      else if (ts(i) .gt. 273.156 .and. ustar(i,j) .lt. usmin) then
        is_dew = .true.

!        print*, 'dew==='


      else
        is_rain = .false.
        is_dew = .false.

!         print*, 'NO dew, NO rain ==='


      end if

      !================================================================
      !Decide fraction of stomatal blocking due to wet conditions
      !================================================================

      wst = 0.0

      if ((is_dew .or. is_rain) .and. srad(i) .gt. 200.0) then

         wst = (srad(i) - 200.0)/800.0
         wst = dmin1(wst, 0.5d0)

      end if
      
      !================================================================
      ! In-canopy aerodynamic resistance
      !================================================================

        rac = rac1(kcov)+(lai_f(i)-   &
                  laimin(i))/(laimax(i)-laimin(i)+1.E-10)  &
                * (rac2(kcov)-rac1(kcov))

!        print*,'rac1 = ', rac

        rac = rac*lai_f(i)**0.25/ustar(i,j)/ustar(i,j)

!        print*,'rac2 = ', rac

      !================================================================
      ! Ground resistance for O3
      !================================================================

       if(ts(i) .lt. 272.15 .and. kcov .ne. 14 ) then

        rgo_f = dmin1( rgo(kcov)*2.0, rgo(kcov) *     &
                                 dexp(0.2*(272.15-ts(i))))

!       print*,'rgo_f1 =',rgo_f, ts(i)
      else
        rgo_f = rgo(kcov)

      end if
        
      !================================================================
      ! Ground resistance for SO2
      !================================================================

      if (kcov .eq. 12) then

         rgs_f = dmin1(rgs(kcov)*(275.15 - ts(i)), 500.d0)
         rgs_f = dmax1(rgs(kcov), 100.d0)

!         print*,'rgs_f ==== ', rgs_f

      else if (is_rain .and. kcov .ne. 14 ) then

         rgs_f = 50.0
!         print*,'rgs_f ==== ', rgs_f

      else if (is_dew .and. kcov .ne. 14 ) then

         rgs_f = 100.0
!         print*,'rgs_f ==== ', rgs_f


       else if(ts(i) .lt. 272.156 .and. kcov .ne. 14 ) then


         rgs_f = dmin1( rgs(kcov)*2.0, rgs(kcov) *     &
                                     dexp(0.2*(272.156 - ts(i))))
!         print*,'rgs_f ==== ', rgs_f

      else

        rgs_f =  rgs(kcov)
!        print*,'rgs_f ==== ', rgs_f


      end if

      !================================================================
      ! Cuticle resistance for O3 AND SO2
      !================================================================

      if (rcutdo(kcov) .le. -1) then
         rcuto_f = 1.E25
         rcuts_f = 1.E25

!         print*,'RCUT === ', rcuto_f,rcuts_f 


      else if (is_rain) then

        rcuto_f = rcutwo(kcov)/lai_f(i)**0.5/ustar(i,j)
        rcuts_f = 50.0/lai_f(i)**0.5/ustar(i,j)
        rcuts_f = dmax1(rcuts_f, 20.d0)
!        print*,'RCUT === ', rcuto_f,rcuts_f 


      else if (is_dew) then

        rcuto_f = rcutwo(kcov)/lai_f(i)**0.5/ustar(i,j)
        rcuts_f = 100.0/lai_f(i)**0.5/ustar(i,j)
        rcuts_f = dmax1 (rcuts_f, 20.d0)
!        print*,'RCUT === ', rcuto_f,rcuts_f 


      else if (ts(i) .lt. 272.156) then

        ryx = dexp(0.2 * (272.156 - ts(i) ))
        rcuto_f = rcutdo(kcov)/dexp(3.0 * rh(i))/     &
                 lai_f(i)**0.25/ustar(i,j)
        rcuts_f = rcutds(kcov)/dexp(3.0 * rh(i))/     &
                 lai_f(i)**0.25/ustar(i,j)
        rcuto_f = dmin1( rcuto_f * 2.0, rcuto_f * ryx )
        rcuts_f = dmin1( rcuts_f*2.0, rcuts_f * ryx )
        rcuto_f = dmax1 (rcuto_f,100.d0)
        rcuts_f = dmax1 (rcuts_f,100.d0)
!        print*,'RCUT === ', rcuto_f,rcuts_f 


      else
        rcuto_f=rcutdo(kcov)/exp(3.*rh(i))/lai_f(i)**0.25/ustar(i,j)
        rcuts_f=rcutds(kcov)/exp(3.*rh(i))/lai_f(i)**0.25/ustar(i,j)
        rcuto_f = dmax1 (rcuto_f, 100.d0)
        rcuts_f = dmax1 (rcuts_f, 100.d0)

!        print*,'RCUT === ', rcuto_f,rcuts_f 

      end if

      !================================================================
      ! If snow occurs, Rg and Rcut are adjusted by snow cover
      ! fraction
      !================================================================

         fsnow= sd(i)/sdmax(kcov)
         fsnow= dmin1(1.0d0, fsnow)   !snow cover fraction for leaves
!         print*,' fsnow=  ', fsnow 

      if (fsnow .gt. 0.0001 .and.kcov .ne. 20 .or.          &
                                 kcov .ne. 15 .or.          & 
                                 kcov .ne. 14 .or.          &
                                 kcov .ne. 12        ) then

        rsnows= dmin1(70.*(275.15-ts(i)), 500.d0)
        rsnows= dmax1(rsnows, 100.d0)
        rcuts_f=1.0/((1.0 - fsnow)/rcuts_f + fsnow/rsnows)
        rcuto_f=1.0/((1.0 - fsnow)/rcuto_f + fsnow/2000.0)
        fsnow  =dmin1(1.d0, fsnow*2.0)
                                        !snow cover fraction for ground
        rgs_f=1.0/((1.0 - fsnow)/rgs_f + fsnow/rsnows)
        rgo_f=1.0/((1.0 - fsnow)/rgo_f + fsnow/2000.0)

!        print*,'rsnows= ', rsnows, ' fsnow=  ', fsnow 
      end if


      !================================================================
      ! Calculate diffusivity for each gas species
      !================================================================

        do ig=1,igas

        dgas = 0.369 * mw(ig) + 6.29
        di=0.001*ts(i)**1.75*sqrt((29.0 + mw(ig))/mw(ig)/29.)
        di=di/1.0/(dair**0.3333 + dgas**0.3333)**2
        vi=145.8 * 1.e-4 * (ts(i) * 0.5 + t2(i) *0.5)**1.5/    &
                          (ts(i) * 0.5 + t2(i) *0.5 +110.4)

      
      !================================================================
      ! Calculate quasi-laminar resistance
      !================================================================

        rb(ig,i,j) =5.0/ustar(i,j) * (vi/di)**.666667

!        print*,'rb==', rb(ig,i,j)

      !================================================================
      ! Calculate stomatal resistance for each species from the ratio
      ! of  diffusity of water vapor to the gas species
      !================================================================

        dvh2o = 0.001*ts(i)**1.75 * sqrt((29.0 + 18.0 )/29.0/18.0)
        dvh2o = dvh2o/(dair**0.3333 + dh2o**0.3333)**2
        rstom = rst * dVh2o/di + rm(ig) 
! (rst <999) for bare surfaces)
      !================================================================
      ! Scale cuticle and ground resistances for each species
      !================================================================

        rcut = 1.0/(alphaz(ig)/rcuts_f+betaz(ig)/rcuto_f)
        rg   = 1.0/(alphaz(ig)/rgs_f+betaz(ig)/rgo_f)

      !================================================================
      ! Calculate total surface resistance
      !================================================================
!FAB
!account for zero stomatal resistance (rst and rstom are zero for bare surfaces)
!set wst to 1 also in that case (total stomatal blocking).             
      
        if(rst== -999.) wst = 1.0

!        rc(ig,i,j) = (1.0 - wst)/rstom + 1.0/(rg)+1.0/rcut
        rc(ig,i,j) = (1.0 - wst)/rstom + 1.0/(rac+rg)+1.0/rcut

        rc(ig,i,j) = dmax1(10.d0,1.0/rc(ig,i,j) )

       end do !igas

       end do !ilg
       end do !luc

      
      end subroutine stomtresis


end module mod_che_drydep	

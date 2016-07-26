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

module mod_che_dust
  !
  ! DUST module
  !
  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_mpmessage
  use mod_memutil
  use mod_mppparam
  use mod_che_common
  use mod_che_ncio
  use mod_che_mppio

  implicit none

  private

  real(rkx) , dimension(4,2) :: dustbsiz1
  real(rkx) , dimension(12,2) :: dustbsiz2

  ! Fix the actual dust aerosol bin size: diameter in microm

  !  data  dustbsiz1 / 0.01_rkx,  1.00_rkx,  2.50_rkx,  5.00_rkx,  1.00_rkx, &
  !                   2.50_rkx,  5.00_rkx, 20.00_rkx/

  data  dustbsiz1 / 0.09_rkx,  1.00_rkx,  2.50_rkx,  5.00_rkx,  1.00_rkx, &
                   2.50_rkx,  5.00_rkx, 63.00_rkx/


  ! new option defined from LISA optimized distribution
  data  dustbsiz2 / 0.09_rkx,  0.18_rkx,  0.60_rkx,  1.55_rkx,  2.50_rkx, &
                    3.75_rkx,  4.70_rkx,  5.70_rkx,  7.50_rkx, 14.50_rkx, &
                   26.00_rkx, 41.00_rkx,  0.18_rkx,  0.60_rkx,  1.55_rkx, &
                    2.50_rkx,  3.75_rkx,  4.70_rkx,  5.70_rkx,  7.50_rkx, &
                   14.50_rkx, 26.00_rkx, 41.00_rkx, 63.00_rkx /

  ! dust effective diameter

  real(rkx) , dimension(4)    ::  dustbed1
  real(rkx) , dimension(12)   ::  dustbed2

  ! has to be calculated from an assumed sub-bin distribution
  !  data dustbed1 /0.82_rkx , 1.8_rkx , 3.7_rkx , 12.5_rkx /

  !FAB
  ! if calculated from Kok distribution
  data dustbed1 /0.658184_rkx, 1.75093_rkx, 3.67936_rkx, 8.46347_rkx /

  ! 12 bins option calculated from Kok Distibution
  data dustbed2 / 0.14062217_rkx,  0.43004150_rkx,  1.10404692_rkx,  &
                  2.02166018_rkx,  3.10952699_rkx,  4.21185749_rkx,  &
                  5.18438211_rkx,  6.55182088_rkx,  9.96016755_rkx,  &
                 16.19150734_rkx, 26.74151275_rkx, 41.32554485_rkx /

  ! solubility of od dust aer for param of giorgi and chameides

  real(rkx) , dimension(4) ::  soldust1
  real(rkx) , dimension(12) ::  soldust2

  data  soldust1 /0.1_rkx , 0.1_rkx , 0.1_rkx , 0.1_rkx /
  data  soldust2 /0.1_rkx , 0.1_rkx , 0.1_rkx , 0.1_rkx, &
                  0.1_rkx , 0.1_rkx , 0.1_rkx , 0.1_rkx, &
                  0.1_rkx , 0.1_rkx , 0.1_rkx , 0.1_rkx /

  ! Basic dust aerosol density (ACE-2 ) in kg/m3
  real(rkx) , parameter :: rhodust = 2650.0_rkx

  integer(ik4) , parameter :: nsoil = 152
  integer(ik4) , parameter :: mode = 5
  integer(ik4) , parameter :: jsoilm = 1
  integer(ik4) , parameter :: jfs = 1
  integer(ik4) , parameter :: ust = 1
  integer(ik4) , parameter :: ndi = 6500

  !choice of emission distribution 1= alfaro/gomes
  !                                2 = Kok + Laurent et al.
  ! ichdustemd
  ! lognormal alfaro parameters
  ! define the aerosol distribution at the emission and the corresponding
  ! weighting factors in fuction of bin sizes
  !
  real(rkx) , parameter :: d1 = 1.5_rkx
  real(rkx) , parameter :: d2 = 6.7_rkx
  real(rkx) , parameter :: d3 = 14.2_rkx
  real(rkx) , parameter :: sigma1 = 1.7_rkx
  real(rkx) , parameter :: sigma2 = 1.2_rkx
  real(rkx) , parameter :: sigma3 = 1.5_rkx

  ! corresponding binding energies/
  real(rkx) , parameter  :: e1 = 3.61_rkx
  real(rkx) , parameter  :: e2 = 3.52_rkx
  real(rkx) , parameter  :: e3 = 3.46_rkx
  !
  ! parameter for alternative Kok emission distribution
  real(rkx) , parameter :: d = 3.4_rkx
  real(rkx) , parameter :: sigmas = 3.0_rkx
  ! Normalization constant
  real(rkx) , parameter :: cv = 12.62_rkx
  real(rkx) , parameter :: lambda = 12.0_rkx

  !FENNEC distribution parameters(Ryder et. al. 2013)
  real(rkx) , parameter :: d1F = 0.05_rkx
  real(rkx) , parameter :: d2F = 0.71_rkx
  real(rkx) , parameter :: d3F = 2.04_rkx
  real(rkx) , parameter :: d4F = 5.28_rkx

  real(rkx) , parameter :: sigma1F = 2.5_rkx
  real(rkx) , parameter :: sigma2F = 1.33_rkx
  real(rkx) , parameter :: sigma3F = 1.45_rkx
  real(rkx) , parameter :: sigma4F = 2.0_rkx

  real(rkx) , parameter :: N1 = 508.27_rkx
  real(rkx) , parameter :: N2 = 8.84_rkx
  real(rkx) , parameter :: N3 = 1.89_rkx
  real(rkx) , parameter :: N4 = 0.54_rkx
  ! soil variable, srel 2 d corresponds to the soil aggregate size distribution
  ! in each texture type.
  real(rkx) , pointer,  dimension(:) :: dustbed , soldust
  real(rkx) , pointer,  dimension(:) :: frac1 , frac2 , frac3 , frac
  real(rkx) , pointer , dimension(:,:,:) :: clay2row2 , sand2row2 , silt2row2
  real(rkx) , pointer , dimension(:,:) :: clayrow2 , sandrow2 , dustbsiz
  real(rkx) , pointer , dimension(:,:,:,:) :: srel2d
  real(rkx) , pointer , dimension(:,:,:) :: dustsotex
  !
  ! Mineralogy fraction of minerals in clay and silt categories
  real(rkx) , pointer , dimension(:,:,:) :: cminer,sminer
  real(rkx) , dimension (4):: cfrac,sfrac   !************jj -----Scanza et al.,2015
  data cfrac  /1._rkx, 0.97_rkx, 0.625_rkx, 0.429_rkx /  
  data sfrac / 0._rkx, 0.03_rkx, 0.375_rkx,  0.571_rkx /
  ! Name of variable changed ! SC. 06.10.2010
  real(rkx) , dimension(nsoil) :: dp_array

  public :: sandrow2
  public :: rhodust
  public :: soldust , dustbed , dustbsiz

  integer(ik4) :: ilg

  public :: allocate_mod_che_dust , inidust , sfflux , clm_dust_tend

#ifdef SINGLE_PRECISION_REAL
  real(rkx) , parameter :: mxarg = 15.0_rkx
#else
  real(rkx) , parameter :: mxarg = 25.0_rkx
#endif

  real(rkx) , parameter :: eps = 1.0e-7_rkx

  contains

    subroutine allocate_mod_che_dust
      implicit none
      if ( ichem == 1 ) then
        call getmem3d(dustsotex,jce1,jce2,ice1,ice2,1,nats,'che_dust:dustsotex')
        call getmem3d(clay2row2,ici1,ici2,1,nats,jci1,jci2,'che_dust:clay2row2')
        call getmem3d(sand2row2,ici1,ici2,1,nats,jci1,jci2,'che_dust:sand2row2')
        call getmem3d(silt2row2,ici1,ici2,1,nats,jci1,jci2,'che_dust:silt2row2')
        call getmem2d(clayrow2,ici1,ici2,jci1,jci2,'che_dust:clayrow2')
        call getmem2d(sandrow2,ici1,ici2,jci1,jci2,'che_dust:sandrow2')
        call getmem4d(srel2d,ici1,ici2,jci1,jci2,1,nsoil,1,nats, &
                      'che_dust:srel2d')
        call getmem2d(dustbsiz,1,nbin,1,2,'che_dust:dustbsiz')
        call getmem1d(dustbed,1,nbin,'che_dust:dustbed')
        call getmem1d(soldust,1,nbin,'che_dust:soldust')
        call getmem1d(frac1,1,nbin,'che_dust:frac1')
        call getmem1d(frac2,1,nbin,'che_dust:frac2')
        call getmem1d(frac3,1,nbin,'che_dust:frac3')
        call getmem1d(frac,1,nbin,'che_dust:frac')
        call getmem3d(cminer,jce1,jce2,ice1,ice2,1,nmine,'che_dust:cminer')
        call getmem3d(sminer,jce1,jce2,ice1,ice2,1,nmine,'che_dust:sminer')

    end if
      ilg = ici2-ici1+1
    end subroutine allocate_mod_che_dust
    !
    !  ***********************************************************
    !  * description of 12- soil categories                  *****
    !  *                                                     *****
    !  * i         cat                     sizing            *****
    !  * ------------------------------------------------    *****
    !  * 1         sand                   coarse             *****
    !  * 2         lomay sand             coarse             *****
    !  * 3         sand lomay             coarse-medium      *****
    !  * 4         silt loma              medium-fine        *****
    !  * 5         silt                   medium             *****
    !  * 6         loam                   fine               *****
    !  * 7         sandy clay loam        coarse-medium-fine *****
    !  * 8         silty clay loam        medium             *****
    !  * 9         clay loam              medium-fine        *****
    !  * 10        sandy clay             coarse-fine        *****
    !  * 11        silty clay             medium-fine        *****
    !  * 12        clay                   fine               *****
    !  ***********************************************************
    !
    subroutine inidust
      implicit none
      real(rkx) , dimension(nats) :: bcly , bslt , bsnd
      real(rkx) :: deldp , stotal , xk , xl , xm , xn
      integer(ik4) :: i , j , n , nm , ns , nt , itr
      real(rkx) , dimension(mode,nats) :: mmdd , pcentd , sigmad
      real(rkx) , dimension(mode,nats) :: mmd , pcent , sigma
      real(rkx) , dimension(iy,nsoil,nats) :: srel
      real(rkx) , dimension(nsoil) :: ss
      real(rkx) , dimension(ndi) :: di
      ! modif new distribution
      ! for each category, this is the percent of Coarse sand,
      ! Fine mode sand, silt , clay and salt ( cf Menut et al. ,2012)
      real(rkx) , dimension (mode,12) :: soiltexpc
      real(rkx) , dimension (mode)    :: texmmd , texstd

      logical :: rd_tex
      character(6) :: aerctl
      real(rkx) :: alogdi , amean1 , amean2 , amean3 , asigma1 , &
             asigma2 , asigma3 , totv1 , totv2 , totv3 , totv ,  &
             exp1 , exp2 , exp3
#ifdef __PGI
      real(rkx) , external :: erf
#endif
      !
      ! FAB update
      ! change type 1 and 2 and 3 to Laurent et al., 2008,
      ! marticorena et al., 1997 soil size parameter.
      ! change also the clay/sand/sil percentage (used to calculate the ratio
      ! of vertical/horizontal flux): This (bcly) is only effective when Kok
      ! size distribution is option enabled.
      ! Values are derived from Laurent et al. typical ranged adapted
      ! for our USDA texture types.
      !
      ! data bcly/0.00_rkx , 0.4e-2_rkx ,0.7e-2_rkx  , 0.7e-2_rkx , &
      !           0.4e-2_rkx , 1.e-2_rkx , 3.e-2_rkx , 3e-2_rkx ,   &
      !           5.e-2_rkx , 8.e-2_rkx , 8.e-2_rkx , 1.e-2_rkx/
      ! data bcly / 4.3e-2_rkx, 2.3e-2_rkx, 7.3e-2_rkx, 0.0_rkx,0.0_rkx, &
      !             0.0e-2_rkx,0.0_rkx, 0.0_rkx,0.0_rkx,0.0_rkx,0.0_rkx,0.0_rkx/
      data bcly / 6.0e-2_rkx , 2.3e-2_rkx , 7.3e-2_rkx , 0.0e-2_rkx , &
                  0.0e-2_rkx , 0.0e-2_rkx , 0.0e-2_rkx , 0.0e-2_rkx , &
                  0.0e-2_rkx , 0.0e-2_rkx , 0.0e-2_rkx , 0.0e-2_rkx /

      ! bsnd and bslt are not really used after /
      ! the data here are not consistent with clay.
      data bsnd / 0.90_rkx , 0.85_rkx , 0.80_rkx , 0.50_rkx , 0.45_rkx , &
                  0.35_rkx , 0.30_rkx , 0.30_rkx , 0.20_rkx , 0.65_rkx , &
                  0.60_rkx , 0.50_rkx /
      data bslt / 0.050_rkx , 0.050_rkx , 0.051_rkx , 0.350_rkx , 0.400_rkx , &
                  0.600_rkx , 0.650_rkx , 0.500_rkx , 0.050_rkx , 0.000_rkx , &
                  0.000_rkx , 0.000_rkx/

      data mmdd /690.0_rkx ,   0.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx, &
                 690.0_rkx , 210.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx, &
                 690.0_rkx , 210.0_rkx , 0.0_rkx , 0.0_rkx , 0.0_rkx, &
                 520.0_rkx , 100.0_rkx , 5.0_rkx , 0.0_rkx , 0.0_rkx, &
                 520.0_rkx ,  75.0_rkx , 2.5_rkx , 0.0_rkx , 0.0_rkx, &
                 520.0_rkx ,  75.0_rkx , 2.5_rkx , 0.0_rkx , 0.0_rkx, &
                 210.0_rkx ,  75.0_rkx , 2.5_rkx , 0.0_rkx , 0.0_rkx, &
                 210.0_rkx ,  50.0_rkx , 2.5_rkx , 0.0_rkx , 0.0_rkx, &
                 125.0_rkx ,  50.0_rkx , 1.0_rkx , 0.0_rkx , 0.0_rkx, &
                 100.0_rkx ,  10.0_rkx , 1.0_rkx , 0.0_rkx , 0.0_rkx, &
                 100.0_rkx ,  10.0_rkx , 0.5_rkx , 0.0_rkx , 0.0_rkx, &
                 100.0_rkx ,  10.0_rkx , 0.5_rkx , 0.0_rkx , 0.0_rkx /

      data sigmad / 1.6_rkx , 1.8_rkx , 1.8_rkx , 0.0_rkx , 0.0_rkx,  &
                    1.6_rkx , 1.8_rkx , 1.8_rkx , 0.0_rkx , 0.0_rkx,  &
                    1.6_rkx , 1.8_rkx , 1.8_rkx , 0.0_rkx , 0.0_rkx,  &
                    1.6_rkx , 1.7_rkx , 1.8_rkx , 0.0_rkx , 0.0_rkx,  &
                    1.6_rkx , 1.7_rkx , 1.8_rkx , 0.0_rkx , 0.0_rkx,  &
                    1.6_rkx , 1.7_rkx , 1.8_rkx , 0.0_rkx , 0.0_rkx,  &
                    1.7_rkx , 1.7_rkx , 1.8_rkx , 0.0_rkx , 0.0_rkx,  &
                    1.7_rkx , 1.7_rkx , 1.8_rkx , 0.0_rkx , 0.0_rkx,  &
                    1.7_rkx , 1.7_rkx , 1.8_rkx , 0.0_rkx , 0.0_rkx,  &
                    1.8_rkx , 1.8_rkx , 1.8_rkx , 0.0_rkx , 0.0_rkx,  &
                    1.8_rkx , 1.8_rkx , 1.8_rkx , 0.0_rkx , 0.0_rkx,  &
                    1.8_rkx , 1.8_rkx , 1.8_rkx,  0.0_rkx , 0.0_rkx /

      data pcentd / 1.00_rkx , 0.00_rkx , 0.00_rkx , 0.0_rkx , 0.0_rkx, &
                    0.90_rkx , 0.10_rkx , 0.00_rkx , 0.0_rkx , 0.0_rkx, &
                    0.80_rkx , 0.20_rkx , 0.00_rkx , 0.0_rkx , 0.0_rkx, &
                    0.50_rkx , 0.35_rkx , 0.15_rkx , 0.0_rkx , 0.0_rkx, &
                    0.45_rkx , 0.40_rkx , 0.15_rkx , 0.0_rkx , 0.0_rkx, &
                    0.35_rkx , 0.50_rkx , 0.15_rkx , 0.0_rkx , 0.0_rkx, &
                    0.30_rkx , 0.50_rkx , 0.20_rkx , 0.0_rkx , 0.0_rkx, &
                    0.30_rkx , 0.50_rkx , 0.20_rkx , 0.0_rkx , 0.0_rkx, &
                    0.20_rkx , 0.50_rkx , 0.30_rkx , 0.0_rkx , 0.0_rkx, &
                    0.65_rkx , 0.00_rkx , 0.35_rkx , 0.0_rkx , 0.0_rkx, &
                    0.60_rkx , 0.00_rkx , 0.40_rkx , 0.0_rkx , 0.0_rkx, &
                    0.50_rkx , 0.00_rkx , 0.50_rkx , 0.0_rkx , 0.0_rkx /
      !!
      !! new option
      !!
      data   soiltexpc / 0.46_rkx, 0.46_rkx, 0.05_rkx, 0.03_rkx, 0.0_rkx, &
                         0.41_rkx, 0.41_rkx, 0.18_rkx, 0.00_rkx, 0.0_rkx, &
                         0.29_rkx, 0.29_rkx, 0.32_rkx, 0.10_rkx, 0.0_rkx, &
                         0.00_rkx, 0.17_rkx, 0.70_rkx, 0.13_rkx, 0.0_rkx, &
                         0.00_rkx, 0.10_rkx, 0.85_rkx, 0.05_rkx, 0.0_rkx, &
                         0.00_rkx, 0.43_rkx, 0.39_rkx, 0.18_rkx, 0.0_rkx, &
                         0.29_rkx, 0.29_rkx, 0.15_rkx, 0.27_rkx, 0.0_rkx, &
                         0.00_rkx, 0.10_rkx, 0.56_rkx, 0.34_rkx, 0.0_rkx, &
                         0.00_rkx, 0.32_rkx, 0.34_rkx, 0.34_rkx, 0.0_rkx, &
                         0.00_rkx, 0.52_rkx, 0.06_rkx, 0.42_rkx, 0.0_rkx, &
                         0.00_rkx, 0.06_rkx, 0.47_rkx, 0.47_rkx, 0.0_rkx, &
                         0.00_rkx, 0.22_rkx, 0.20_rkx, 0.58_rkx, 0.0_rkx/

      data  texmmd  / 690.0_rkx, 210.0_rkx, 125.0_rkx,2.0_rkx, 520.0_rkx /
      data  texstd  / 1.6_rkx,   1.6_rkx,   1.8_rkx,  2.0_rkx, 1.50_rkx /

      mmd = d_zero
      sigma = d_zero
      pcent = d_zero
      if ( ichdustemd == 1 ) then
        mmd = mmdd
        sigma = sigmad
        pcent = pcentd
      else if ( ichdustemd == 2 ) then
        do nm = 1 , mode
          mmd(nm,:) = texmmd(nm)
          sigma(nm,:) = texstd(nm)
          pcent(nm,:) = soiltexpc(nm,:)
        end do
        where ( abs(soiltexpc(:,:)) < dlowval )
          mmd (:,:) = d_zero
          sigma(:,:) = d_zero
        end where
      end if

      rd_tex = .false.
      do itr = 1 , ntr
        aerctl = chtrname(itr)
        if ( aerctl(1:2) == 'DU' ) then
          rd_tex = .true.
          exit
        end if
      end do

      call bcast(rd_tex)

      if ( rd_tex ) then
        call read_texture(nats,dustsotex)
      end if

      ! read mineral fractions
      if ( imine(1,1) > 0 ) then
        call read_miner(nmine,cminer, sminer)
      end if

      clay2row2 = d_zero
      clayrow2  = d_zero
      sand2row2 = d_zero
      silt2row2 = d_zero
      srel2d    = d_zero

      do j = jci1 , jci2
        do i = ici1 , ici2
          do nt = 1 , nats
            clay2row2(i,nt,j) = bcly(nt)*d_100
            sand2row2(i,nt,j) = bsnd(nt)*d_100
            silt2row2(i,nt,j) = bslt(nt)*d_100
            sandrow2(i,j) = sandrow2(i,j) + dustsotex(j,i,nt)* sand2row2(i,nt,j)
            clayrow2(i,j) = clayrow2(i,j) + dustsotex(j,i,nt)*clay2row2(i,nt,j)
          end do
        end do
      end do ! end j loop

      dp_array(1) = 0.0001_rkx  !cm
      do ns = 2 , nsoil
        dp_array(ns) = dp_array(ns-1)*exp(0.0460517018598807_rkx)
        deldp = dp_array(ns) - dp_array(ns-1)
      end do

      di(1) = 0.01_rkx !microm
      do ns = 2 , ndi
        di(ns) = di(ns-1) + 0.01_rkx
      end do

      do j = jci1 , jci2
        srel(:,:,:) = d_zero
        do i = ici1 , ici2
          do nt = 1 , nats
            ss(:) = d_zero
            stotal = d_zero
            if ( sand2row2(i,nt,j) > d_zero ) then
              do ns = 1 , nsoil          !soil size segregatoin no
                do nm = 1 , mode       !soil mode = 5
                  if ( (pcent(nm,nt) > eps) .and.                    &
                       (sigma(nm,nt) /= d_zero) ) then
                    xk = pcent(nm,nt)/(sqrt(twopi)*log(sigma(nm,nt)))
                    xl = ((log(dp_array(ns))- &
                           log(mmd(nm,nt)*1.0e-4_rkx))**2) / &
                           (d_two*(log(sigma(nm,nt)))**2)
                    if ( xl > mxarg ) then
                      xm = d_zero
                    else
                      xm = xk*exp(-xl)
                    end if
                  else
                    xm = d_zero
                  end if
                  xn = rhodust*twot*(dp_array(ns)*d_half)
                  deldp = 0.0460517018598807_rkx
                  ! dp_array(2)-dp_array(1) ss(nsoil)
                  ss(ns) = ss(ns) + (xm*deldp/xn)
                end do
                stotal = stotal + ss(ns)
              end do
              do ns = 1 , nsoil
                if ( stotal > d_zero ) srel(i,ns,nt) = ss(ns)/stotal
              end do
            end if
          end do ! soil types
        end do

        do nt = 1 , nats
          do i = ici1 , ici2
            do ns = 1 , nsoil
              srel2d(i,j,ns,nt) = srel(i,ns,nt)
            end do
          end do
        end do

      end do  ! end J loop
      !
      ! Finally calculate the emission stribution weights in function of
      ! distribution parameters:(Alfaro, Kok)
      !
      if ( size(idust) == size(dustbsiz1,1)   ) then
        dustbsiz(:,:) = dustbsiz1(:,:)
        dustbed(:) = dustbed1(:)
        chtrsol(idust(:))=  soldust1(:)
      elseif ( size(idust) == size(dustbsiz2,1) ) then
        dustbsiz(:,:) = dustbsiz2(:,:)
        dustbed(:) = dustbed2(:)
        chtrsol(idust(:))=  soldust2(:)
      end if

      if ( ichdustemd == 1 ) then
        totv1 = d_zero
        totv2 = d_zero
        totv3 = d_zero
        amean1 = log10(d1)
        amean2 = log10(d2)
        amean3 = log10(d3)
        asigma1 = log10(sigma1)
        asigma2 = log10(sigma2)
        asigma3 = log10(sigma3)
        do ns = 1, ndi
          alogdi = log10(di(ns))
          exp1 = (alogdi-amean1)**2/(d_two*asigma1**2)
          exp2 = (alogdi-amean2)**2/(d_two*asigma2**2)
          exp3 = (alogdi-amean3)**2/(d_two*asigma3**2)
          do n = 1 , nbin
            if ( di(ns) > dustbsiz(n,1) .and. &
                 di(ns) <= dustbsiz(n,2) ) then
              ! the independant variable is diameter so going from
              ! dV/log10D to dV/dD implies a factor 1/(2.303)D
              if ( exp1 < mxarg ) then
                frac1(n) = frac1(n) + (d_one/di(ns)) * exp(-exp1)
              end if
              if ( exp2 < mxarg ) then
                frac2(n) = frac2(n) + (d_one/di(ns)) * exp(-exp2)
              end if
              if ( exp3 < mxarg ) then
                frac3(n) = frac3(n) + (d_one/di(ns)) * exp(-exp3)
              end if
            end if
          end do
          if ( exp1 < mxarg ) then
            totv1 = totv1 + (d_one/di(ns)) * exp(-exp1)
          end if
          if ( exp2 < mxarg ) then
            totv2 = totv2 + (d_one/di(ns)) * exp(-exp2)
          end if
          if ( exp3 < mxarg ) then
            totv3 = totv3 + (d_one/di(ns)) * exp(-exp3)
          end if
        end do
        frac1(:) = frac1(:) / totv1
        frac2(:) = frac2(:) / totv2
        frac3(:) = frac3(:) / totv3
      else if ( ichdustemd >= 2 ) then
        ! calculate the bin mass fraction to the total mass from Kok et al.
        ! distribution ( mass distribution).
        ! the independant variable is diameter so going from
        ! dV/dlnD to dV/dD implies a factor 1/D
        frac = 0._rkx
        totv = 0._rkx
        do ns = 1, ndi
          do n = 1, nbin
             if ( di(ns) > dustbsiz(n,1) .and. di(ns) <= dustbsiz(n,2) ) then
                frac(n) = frac(n) + d_one/cv * &
                  (d_one+erf(log(di(ns)/d)/sqrt(d_two)/ &
                  log(sigmas)))*exp(-(di(ns)/lambda)**3)  !see Kok (2011)
             end if
           end do
           totv = totv + d_one / cv * (d_one+erf(log(di(ns)/d)/sqrt(d_two)/ &
                  log(sigmas)))*exp(-(di(ns)/lambda)**3)
        end do
        frac(:) = frac(:) / totv
      end if

    end subroutine inidust
    !
    !   *****************************************************************
    !   * calculate of ustar01(d) using iversen and white (1982)     ****
    !   * for smoth surface:                                         ****
    !   * coded by :                                                 ****
    !   * dum    : particle diameter [um]                            ****
    !   * ustar0 : threshold frication velocity [m/s]                ****
    !   *****************************************************************
    !
    real(rkx) function ustart01(rhodust,dum,rhair)
      implicit none
      real(rkx) , intent(in) :: dum , rhair , rhodust
      real(rkx) , parameter :: a2 = 0.129_rkx
      real(rkx) , parameter :: c1 = 0.006_rkx
      real(rkx) , parameter :: c2 = 1.928_rkx
      real(rkx) , parameter :: c3 = 0.0858_rkx
      real(rkx) , parameter :: c4 = -0.0617_rkx
      real(rkx) , parameter :: c5 = 2.5_rkx
      real(rkx) , parameter :: y1 = 1331.647_rkx
      real(rkx) , parameter :: y2 = 1.561228_rkx
      real(rkx) , parameter :: y3 = 0.38194_rkx
      real(rkx) :: dm , rep , term , term1 , term2

      dm = dum  !* 1.0e-4      ! cm
      rep = y1*(dm**y2) + y3
      term1 = sqrt(d_one+(c1/(rhodust*egrav*0.1_rkx*(dm**c5))))
      term2 = sqrt(rhodust*egrav*d_100*dm/rhair)
      term = term1*term2
      ustart01 = cvmgt(a2*term*(d_one-c3*exp(c4*(rep-d_10))),  &
                 a2*term/sqrt(c2*(rep**0.092_rkx)-d_one),rep > d_10)
    contains

      real(rkx) function cvmgt(val1,val2,cond)
        implicit none
        logical , intent(in) :: cond
        real(rkx) , intent(in) :: val1 , val2
        if ( cond ) then
          cvmgt = val1
        else
          cvmgt = val2
        end if
      end function cvmgt

    end function ustart01
    !
    !   *****************************************************************
    !   *                                                            ****
    !   * y. shao, 13 june 2000                                      ****
    !   * calculate ustar0(d) using shao and lu (2000) for uncovered ****
    !   * dry surface                                                ****
    !   * dum:    particle diameter                   [um]           ****
    !   * ustar0: threshold friction velocity       [cm/s]           ****
    !   *****************************************************************
    !
    real(rkx) function ustart0(rhodust,dum,rhoa)
      implicit none
      real(rkx) , intent(in) :: dum , rhoa , rhodust
      real(rkx) , parameter :: agamma = 3.0e-4_rkx , f = 0.0123_rkx
      real(rkx) :: dm , sigma
      sigma = rhodust/rhoa
      dm = dum*1.0e-2_rkx
      ustart0 = f*(sigma*egrav*dm+agamma/(rhoa*dm))
      ustart0 = sqrt(ustart0)
      ustart0 = ustart0*d_100
    end function ustart0
    !
    !   **********************************************************
    !   *  dust emission scheme                             ******
    !   *                                                   ******
    !   * this scheme based on marticorena and bergametti,  ******
    !   * 1995; gong et al.,(2003); alfaro et al.,(1997)    ******
    !   * Zakey et al., 2006                                ******
    !   **********************************************************
    !
    subroutine sfflux(jloop,ivegcov,vegfrac,ustarnd,z0,soilw, &
                      surfwd,roarow,trsize,rsfrow)
      implicit none
      integer(ik4) , intent(in) :: jloop
      integer(ik4) , intent(in) , dimension(ici1:ici2) ::  ivegcov
      real(rkx) , intent(in) , dimension(ici1:ici2) :: roarow , soilw , &
                surfwd , vegfrac , z0 , ustarnd
      real(rkx) , intent(out) , dimension(ici1:ici2,nbin) :: rsfrow
      real(rkx) , intent(in) , dimension(nbin,2) :: trsize
      real(rkx) , dimension(ilg) :: xclayrow , xroarow , xsoilw , &
                xsurfwd , xvegfrac , xz0 , xustarnd , xsnowfrac
      real(rkx) , dimension(ilg,nbin) :: xrsfrow
      real(rkx) , dimension(ilg,nats) :: xftex , xalphaprop
      real(rkx) , dimension(ilg,nsoil,nats) :: xsrel2d
      integer(ik4) :: i , ieff , n , ns

      rsfrow = d_zero
      ! effective emitter cell ( depending on ivegcov)
      xvegfrac = d_zero
      xsnowfrac = d_zero
      xftex = d_zero
      xsoilw = d_zero
      xsurfwd = d_zero
      xz0 = d_zero
      xclayrow = d_zero
      xroarow = d_zero
      xsrel2d = d_zero
      xustarnd=d_zero
      xrsfrow = d_zero
      xalphaprop = d_zero

      ieff = 0
      do i = ici1 , ici2
        if (ivegcov(i) == 8 .or. ivegcov(i) == 11) then
          ieff = ieff + 1
          xvegfrac(ieff) = vegfrac(i)
          xsnowfrac(ieff) =  csfracs2d(jloop,i)
          xsoilw(ieff) = soilw(i)
          xsurfwd(ieff) = surfwd(i)
          xz0(ieff) = z0(i)
          xroarow(ieff) = roarow(i)
          xustarnd(ieff) = ustarnd(i)
          xclayrow(ieff) = clayrow2(i,jloop)
          do n = 1 , nats
            xftex(ieff,n) = dustsotex(jloop,i,n)
            if ( ichdustemd == 2 ) then
              if ( clay2row2(i,n,jloop) <= 20 ) then
                xalphaprop(ieff,n) = d_10**(0.134_rkx * &
                              clay2row2(i,n,jloop)-6.0_rkx)
!                             clay2row2(i,n,jloop)-6.0_rkx)*0.035_rkx
              else
                xalphaprop(ieff,n) = d_10**(-0.1_rkx * &
                              clay2row2(i,n,jloop)-6.0_rkx)
!                             clay2row2(i,n,jloop)-1.2_rkx)*0.035_rkx
              end if
            end if
            do  ns = 1 , nsoil
              xsrel2d(ieff,ns,n) = srel2d(i,jloop,ns,n)
            end do
          end do
        end if
      end do

      if ( ieff > 0 ) then
        call dust_module(1,ieff,trsize,xsoilw,xvegfrac,xsnowfrac,xsurfwd, &
                         xftex,xclayrow,xroarow,xalphaprop,xz0, &
                         xsrel2d,xustarnd,xrsfrow)
      end if

      ! put back the dust flux on the right grid

      ieff = 1
      do i = ici1 , ici2
        if  (ivegcov(i) == 8 .or. ivegcov(i) == 11) then
          do n = 1 , nbin
            rsfrow(i,n) = xrsfrow(ieff,n)
            if ( ichdrdepo == 1 ) then
              chiten(jloop,i,kz,idust(n)) = chiten(jloop,i,kz,idust(n)) + &
                   rsfrow(i,n)*egrav/(dsigma(kz)*1.e3_rkx)
            elseif ( ichdrdepo == 2 ) then
              ! pass the flux to BL scheme
              chifxuw(jloop,i,idust(n)) = chifxuw(jloop,i,idust(n)) + &
                   rsfrow(i,n)
            end if
            ! diagnostic source (accumulated)
            cemtrac(jloop,i,idust(n)) = cemtrac(jloop,i,idust(n)) + &
                     rsfrow(i,n)* cfdout

             if ( ichdiag == 1 ) then
             cemisdiag(jloop,i,idust(n)) = cemisdiag(jloop,i,idust(n)) + &
                       rsfrow(i,n) / &
                       (cdzq(jloop,i,kz)*crhob3d(jloop,i,kz)) * cfdout
             end if
          end do
          ieff = ieff + 1
        end if
      end do
      ! Mineralogy flux option
      ! introduce mineralogy here, implicit loop on mineral types
      if ( imine(1,1) > 0 ) then
        do n = 1 , nbin 
          do i = ici1 , ici2
            chiten(jloop,i,kz,imine(n,:)) = &
                   chiten(jloop,i,kz,imine(n,:)) + &
                   rsfrow(i,n)* ( cminer(jloop,i,:)*cfrac(n) + &
                                  sminer(jloop,i,:)*sfrac(n) ) &
                   *egrav / (dsigma(kz)*1.e3_rkx)


            cemtrac(jloop,i,imine(n,:)) = cemtrac(jloop,i,imine(n,:)) + &
                     rsfrow(i,n)* (cminer(jloop,i,:) *cfrac(n)  + &
                                   sminer(jloop,i,:) *sfrac(n)) &
                     * cfdout

          end do
        end do
      end if
    end subroutine sfflux

    subroutine dust_module(il1,il2,trsize,soilw,vegfrac,snowfrac,surfwd,ftex, &
                           clayrow,roarow,alphaprop,z0,srel,ustarnd,rsfrow)
      implicit none
      integer(ik4) :: il1 , il2
      real(rkx) , dimension(ilg) :: clayrow , roarow , soilw , surfwd ,   &
                                   vegfrac , z0 , ustarnd , snowfrac
      real(rkx) , dimension(ilg,nbin) :: rsfrow
      real(rkx) , dimension(ilg,nats) :: ftex , alphaprop
      real(rkx) , dimension(ilg,nsoil,nats) :: srel
      real(rkx) , dimension(nbin,2) :: trsize
      intent (in) clayrow , soilw , surfwd , z0 , ustarnd , ftex

      real(rkx) , dimension(ilg) :: alamda , hc , rc , srl , wprim
      real(rkx) :: arc1 , arc2 , cly1 , cly2 , tempd , &
          ustarns , uth , utmin , ustarfw
      integer(ik4) :: i
      real(rkx) , dimension(ilg) :: ustar
      real(rkx) , dimension(ilg,nsoil) :: utheff

      real(rkx) , parameter :: umin = 15.0_rkx
      real(rkx) , parameter :: xz = 0.25_rkx
      real(rkx) , parameter :: br = 202.0_rkx
      real(rkx) , parameter :: ym = 0.16_rkx
      real(rkx) , parameter :: sigr = 1.45_rkx
      real(rkx) , parameter :: z0s = 3.0e-3_rkx
      real(rkx) , parameter :: x = d_10

      do i = il1 , il2

        srl(i) = z0(i)*d_100
        rc(i) = d_one

        if ( jfs == 0 ) then
          ! * raupach et al. (1993)
          if ( vegfrac(i) < d_one ) then
            alamda(i) = xz*(log(d_one-vegfrac(i)))*(-d_one)
            arc1 = sigr*ym*alamda(i)
            arc2 = br*ym*alamda(i)
            if ( arc1 <= d_one .and. arc2 <= d_one ) then
              rc(i) = (sqrt(d_one-arc1)*sqrt(d_one+arc2))
            end if
          end if
        else if ( jfs == 1 ) then
          ! Marticorena et al., 1997: correction factor for non
          ! erodible elements
          rc(i) = d_one - (log(0.50e-2_rkx/z0s) / &
                          (log(0.35_rkx*(x/z0s)**0.8_rkx)))
        end if
        ! threshold velocity correction for soil humidity hc
        if ( jsoilm == 0 ) then
          if ( soilw(i) < d_zero ) then
            write(stderr,*) 'hc, rc = ' , soilw(i) , ' less than zero'
            call fatal(__FILE__,__LINE__,'NEGATIVE SOILW')
          else if ( soilw(i) < 0.03_rkx ) then
            hc(i) = exp(22.7_rkx*soilw(i))
          else if ( soilw(i) >= 0.03_rkx ) then
            hc(i) = exp(95.3_rkx*soilw(i)-2.029_rkx)
          else
            hc(i) = d_one
          end if
        else if ( jsoilm == 1 ) then
          cly1 = clayrow(i)
          cly2 = cly1*cly1
          wprim(i) = 0.0014_rkx*cly2 + 0.17_rkx*cly1
          tempd =  max(0.00001_rkx,soilw(i)*d_100 -wprim(i))
!         print*,'humidity',i,cly1,soilw(i)*100,wprim(i),tempd
          if ( soilw(i)*d_100 > wprim(i) ) then
            hc(i) = sqrt(d_one+1.21_rkx*tempd**0.68_rkx)
!           print*,'hc',i,hc(i)
          else
            hc(i) = d_one
          end if
          ! no soil humidity correction facor if jsoilm > 1
        else
          hc(i)=d_one
        end if
        ! * total correction factor for both hc and rc
        rc(i) = rc(i)/hc(i)
        ! * computation of the wind friction velocity
        ! * accounting for the increase of the roughness length
        ! * due to the saltation layer (gillette etal. jgr 103,
        ! * no. d6, p6203-6209, 1998
        ustarfw = (vonkar*100.0_rkx*surfwd(i))/(log(1000.0_rkx/srl(i)))
        ustarns = ustarnd(i)*d_100 !cm.s-1
        utmin = (umin/(d_100*vonkar*rc(i)))*log(d_1000/srl(i))
        if ( surfwd(i) >= utmin ) then
          ustar(i) = ustarns + 0.3_rkx*(surfwd(i)-utmin)*(surfwd(i)-utmin)
        else
          ustar(i) = ustarns
        end if
      end do       ! end i loop

      call uthefft(il1,il2,ust,nsoil,roarow,utheff,rhodust)

      call emission(il1,il2,rhodust,ftex,alphaprop,uth,roarow,rc,utheff, &
                    ustar,srel,rsfrow,vegfrac,snowfrac)

    end subroutine dust_module

    subroutine uthefft(il1,il2,ust,nsoil,roarow,utheff,rhodust)
      implicit none
      integer(ik4) :: il1 , il2 , nsoil , ust
      real(rkx) :: rhodust
      real(rkx) , dimension(ilg) :: roarow
      real(rkx) , dimension(ilg,nsoil) :: utheff
      intent (in) il1 , il2 , nsoil , ust
      intent (out) utheff
      integer(ik4) :: n , i
      do n = 1 , nsoil
        do i = il1 , il2
          if ( ust == 1 ) utheff(i,n) = ustart0(rhodust,dp_array(n),roarow(i))
          if ( ust == 1 ) utheff(i,n) = ustart01(rhodust,dp_array(n),roarow(i))
        end do
      end do
    end subroutine uthefft

    subroutine emission(il1,il2,rhodust,ftex,alphaprop,uth,roarow,rc, &
                        utheff,ustar,srel,rsfrow,vegfrac,snowfrac)
      implicit none
      integer(ik4) :: il1 , il2
      real(rkx) :: rhodust , uth
      real(rkx) , dimension(ilg) :: rc , ustar, roarow , vegfrac , snowfrac
      real(rkx) , dimension(ilg,nbin) :: rsfrow
      real(rkx) , dimension(ilg,nats) :: ftex , alphaprop
      real(rkx) , dimension(ilg,nsoil,nats) :: srel
      real(rkx) , dimension(ilg,nsoil) :: utheff
      intent (in)  il1 , il2 , rc , rhodust , roarow , srel ,  &
                   ustar , utheff , vegfrac, ftex
      intent (inout) rsfrow , uth
      real(rkx) :: beta , p1 , p2 , p3 , dec , ec , fdp1 , fdp2
      real(rkx) , dimension(ilg,nats) :: fsoil , fsoil1 , fsoil2 , fsoil3
      integer(ik4) :: i , k , n , nt , ns

      real(rkx), dimension(ilg,nbin,nats):: rsfrowt

      data beta  /16300.0_rkx/

      !
      ! Put rdstemfac consistent with soil parameters and Laurent et al., 08
      !
      ! rdstemfac = d_one
      !
      p1 = d_zero
      p2 = d_zero
      p3 = d_zero
      fsoil(:,:) = d_zero
      fsoil1(:,:) = d_zero
      fsoil2(:,:) = d_zero
      fsoil3(:,:) = d_zero

      do nt = 1 , nats
         do i = il1 , il2
           if (ftex(i,nt) < 1.e-10_rkx) cycle
           do ns = 1 , nsoil
            if ( rc(i) > d_zero .and. ustar(i) /= d_zero ) then
              uth = utheff(i,ns)/(rc(i)*ustar(i))
              if ( uth <= d_one ) then
                fdp1 = ustar(i)**3*(d_one-uth*uth)
                fdp2 = (d_one+uth)*rdstemfac*(1.0e-5_rkx)*roarow(i)*regrav
                if ( fdp2 <= d_zero ) fdp2 = d_zero
                ! FAB: with subgrid soil texture, the aggregation of vertical
                ! fluxes per texture type at the grid cell level is done in
                ! fine.
                ! fsoil(k) = srel(k,j,i)*fdp1*fdp2*aeffect*beffect
                ! FAB
                if ( ichdustemd == 1 ) then
                  fsoil(i,nt) = srel(i,ns,nt)*fdp1*fdp2
                  ! size-distributed kinetic energy flux(per texture type)
                  dec = fsoil(i,nt)*beta
                  ! individual kinetic energy for an aggregate of size dp (
                  ! g cm2 s-2) cf alfaro (dp) is in cm
                  ec = (mathpi/12.0_rkx)*rhodust*1.0e-3_rkx*(dp_array(ns)**3)* &
                        (20.0_rkx*ustar(i))**2
                  if ( ec > e1 ) then
                    p1 = (ec-e1)/(ec-e3)
                    p2 = (d_one-p1)*(ec-e2)/(ec-e3)
                    p3 = d_one - p1 - p2
                  else if ( ec > e2 .and. ec <= e1 ) then
                    p1 = d_zero
                    p2 = (ec-e2)/(ec-e3)
                    p3 = d_one - p2
                  else if ( ec > e3 .and. ec <= e2 ) then
                    p1 = d_zero
                    p2 = d_zero
                    p3 = d_one
                  else if ( ec <= e3 ) then
                    p1 = d_zero
                    p2 = d_zero
                    p3 = d_zero
                  end if
                  fsoil1(i,nt) = fsoil1(i,nt) + 1.0e-2_rkx*p1*(dec/e1)* &
                            (mathpi/6.0_rkx)*rhodust*((d1*1.0e-4_rkx)**3)
                  fsoil2(i,nt) = fsoil2(i,nt) + 1.0e-2_rkx*p2*(dec/e2)* &
                            (mathpi/6.0_rkx)*rhodust*((d2*1.0e-4_rkx)**3)
                  fsoil3(i,nt) = fsoil3(i,nt) + 1.0e-2_rkx*p3*(dec/e3)* &
                            (mathpi/6.0_rkx)*rhodust*((d3*1.0e-4_rkx)**3)
                else if ( ichdustemd == 2 ) then
                  fsoil(i,nt) = fsoil(i,nt) + alphaprop(i,nt)* &
                                srel(i,ns,nt)*fdp1*fdp2
                end if
              end if
            end if
          end do
        end do
      end do
      !
      ! calculate fluxes for each of transport bins
      !
      rsfrowt(:,:,:) = d_zero
      if ( ichdustemd == 1 ) then
        do nt = 1 , nats
          do n = 1 , nbin
            do i = il1 , il2
              rsfrowt (i,n,nt) = fsoil1(i,nt)*frac1(n) + &
                                 fsoil2(i,nt)*frac2(n) + &
                                 fsoil3(i,nt)*frac3(n)
            end do
          end do
        end do
      else if ( ichdustemd == 2 ) then
        do nt = 1 , nats
          do n = 1 , nbin
            do i = il1 , il2
              rsfrowt (i,n,nt) = fsoil(i,nt)*frac(n)
            end do
          end do
        end do
      end if
      !
      ! Finally, aggregation of the dust flux at the grid cell level =
      ! weighted sum over soil texture
      ! weighting by grid cell veg fraction, snow fraction,
      ! EBL = erodibility factor, to be introduced )
      ! f = d_zero
      ! aeffect = (1-f)*(1-vegfrac(k))
      ! beffect = 0.01*fland(k,i)*sand2row(k,i)
      ! Fab : fland is equal to 1 with bats
      ! the fraction of sand ( coarse particles) is intrinsically contained
      ! in dsrel soil aggregate distribution
      ! there is no need to multipky by sand2row.
      !
      rsfrow(:,:) = d_zero
      do k = 1 , nbin
        do nt = 1 , nats
          do i = il1 , il2
            rsfrow(i,k) =  rsfrow(i,k) + rsfrowt(i,k,nt)*ftex(i,nt) * &
                           (d_one - vegfrac(i))*(d_one - snowfrac(i))
            ! * EBL(i)
            ! * (1-snowfrac)
          end do
        end do
      end do
    end subroutine emission

    subroutine clm_dust_tend
      implicit none
      integer(ik4) :: i,j,n
      real(rkx) , pointer , dimension(:,:) :: sumdflux
      real(rkx) :: cdsfrq
      ! Update dust tendency with dust fluxes calculated in CLM
      ! here sump up the total flux from clm ( initially defined on 4 bins)
      ! and re-distribute it according to the selected dust emission size
      ! distrib frac(n) is the weight relatibe to bin n ..
      ! this insure consistency between emission distrbution and effective
      ! readius and optical properties
      ! use the same tuning erodibility factor rdstemfac than for
      ! standard scheme
      cdsfrq = real(ntsrf/kche,rkx)
      allocate(sumdflux(jci1:jci2,ici1:ici2))
      sumdflux = d_zero
#ifdef CLM45
      sumdflux = sum(cdustflx_clm,3) * rdstemfac
#endif
      do j = jci1 , jci2
        do i = ici1 , ici2
           do n = 1 , nbin
            if ( ichdrdepo == 1 ) then
              chiten(j,i,kz,idust(n)) = chiten(j,i,kz,idust(n)) + &
                   sumdflux(j,i)*frac(n) * egrav/(dsigma(kz)*1.e3_rkx)
            else if ( ichdrdepo == 2 ) then
              ! pass the flux to BL scheme
              chifxuw(j,i,idust(n)) = chifxuw(j,i,idust(n)) + &
                  sumdflux(j,i) * frac(n)
            end if
            ! diagnostic source (accumulated)

            cemtrac(j,i,idust(n)) = cemtrac(j,i,idust(n)) + &
                    sumdflux(j,i)*frac(n) * cdsfrq
             if ( ichdiag == 1 ) then
               cemisdiag(j,i,idust(n)) = cemisdiag(j,i,idust(n)) + &
                                 sumdflux(j,i)*frac(n) / &
                                 (cdzq(j,i,kz)*crhob3d(j,i,kz)) * cdsfrq
             end if
          end do
        end do
      end do
      deallocate(sumdflux)
    end subroutine clm_dust_tend

end module mod_che_dust

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

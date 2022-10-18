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
#ifdef CLM45
  real(rkx) , pointer , dimension(:,:) :: sumdflux
#endif

  ! Fix the actual dust aerosol bin size: diameter in microm

  data  dustbsiz1 / 0.01_rkx,  1.00_rkx,  2.50_rkx,  5.00_rkx,  1.00_rkx, &
                    2.50_rkx,  5.00_rkx, 20.00_rkx/

  !data  dustbsiz1 / 0.09_rkx,  1.00_rkx,  2.50_rkx,  5.00_rkx,  1.00_rkx, &
  !                 2.50_rkx,  5.00_rkx, 63.00_rkx/


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
  ! data dustbed1 /0.82_rkx , 1.8_rkx , 3.7_rkx , 12.5_rkx /

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
  !integer(ik4) , parameter :: jsoilm = 1 ! TEST PLEASE REMOVE
  integer(ik4) , parameter :: jsoilm = 0

  integer(ik4) , parameter :: ust = 1
  integer(ik4) , parameter :: ndi_4 = 4000
  integer(ik4) , parameter :: ndi_12 = 6500

  !choice emission scheme and size distribution 1=  MB95  + alfaro/gomes
  !                                             2 = Kok et al, 2014, 2011 with thresholf friction vel based on MB95
  ! ichdustemd is set in regcm.in
  !
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
  real(rkx) , pointer , dimension(:) :: fclay, fsand
  real(rkx) , pointer , dimension(:,:) :: clayrow2 , dustbsiz
  real(rkx) , pointer , dimension(:,:,:,:) :: srel2d
  real(rkx) , pointer , dimension(:,:,:) :: dustsotex
  !
  ! erodibility : source function parameter
  real(rkx) , pointer , dimension(:,:) :: erodfc
  ! aeolian roughness from satellite obs
  real(rkx) , pointer , dimension(:,:,:) :: aez0

  ! Mineralogy fraction of minerals in clay and silt categories
  real(rkx) , pointer , dimension(:,:,:) :: cminer , sminer

  !************ realative to mineralogy option -----Scanza et al.,2015
  real(rkx) , dimension (4):: cfrac , sfrac
  data cfrac  /1._rkx, 0.97_rkx, 0.625_rkx, 0.429_rkx /
  data sfrac / 0._rkx, 0.03_rkx, 0.375_rkx, 0.571_rkx /
  ! Name of variable changed ! SC. 06.10.2010
  real(rkx) , dimension(nsoil) :: dp_array

  public :: rhodust
  public :: soldust , dustbed , dustbsiz
  public :: dustsotex, fsand
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
        call getmem3d(dustsotex,jci1,jci2,ici1,ici2,1,nats,'che_dust:dustsotex')
        call getmem1d(fclay,1,nats,'che_dust:fclay')
        call getmem1d(fsand,1,nats,'che_dust:fsand')
        call getmem2d(clayrow2,jci1,jci2,ici1,ici2,'che_dust:clayrow2')
        call getmem4d(srel2d,1,nsoil,1,nats, &
                      jci1,jci2,ici1,ici2,'che_dust:srel2d')
        call getmem2d(dustbsiz,1,nbin,1,2,'che_dust:dustbsiz')
        call getmem2d(erodfc,jci1,jci2,ici1,ici2,'che_dust:erodfc')
        call getmem3d(aez0,jci1,jci2,ici1,ici2,1,12,'che_dust:aez0')
        call getmem1d(dustbed,1,nbin,'che_dust:dustbed')
        call getmem1d(soldust,1,nbin,'che_dust:soldust')
        call getmem1d(frac1,1,nbin,'che_dust:frac1')
        call getmem1d(frac2,1,nbin,'che_dust:frac2')
        call getmem1d(frac3,1,nbin,'che_dust:frac3')
        call getmem1d(frac,1,nbin,'che_dust:frac')
        if ( nmine > 0 ) then
          call getmem3d(cminer,jci1,jci2,ici1,ici2,1,nmine,'che_dust:cminer')
          call getmem3d(sminer,jci1,jci2,ici1,ici2,1,nmine,'che_dust:sminer')
        end if
#ifdef CLM45
        if ( ichdustemd == 3 ) then
          call getmem2d(sumdflux,jci1,jci2,ici1,ici2,'che_dust:sumdflux')
        end if
#endif
      end if
      ilg = (ici2-ici1+1)*(jci2-jci1+1)
    end subroutine allocate_mod_che_dust
    !
    !  ***********************************************************
    !  * description of 12- soil categories                  *****
    !  *                                                     *****
    !  * i         cat                     sizing            *****
    !  * ------------------------------------------------    *****
    !  * 1         sand                   coarse             *****
    !  * 2         lomay sand             coarse             *****
    !  * 3         sand loamy             coarse-medium      *****
    !  * 4         silt loay              medium-fine        *****
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
      real(rkx) :: deldp , stotal , xk , xl , xm , xn
      integer(ik4) :: i , j , n , nm , ns , nt , itr , ndi
      real(rkx) , dimension(mode,nats) :: mmd , pcent , sigma
      real(rkx) , dimension(nsoil) :: ss
      real(rkx) , dimension(:) , allocatable :: di
      real(rkx) , dimension (mode,12) :: soiltexpc
      real(rkx) , dimension (mode)    :: texmmd , texstd

      logical :: rd_tex
      integer(ik4) , dimension(1) :: p
      character(6) :: aerctl
      real(rkx) :: alogdi , amean1 , amean2 , amean3 , asigma1 , &
             asigma2 , asigma3 , totv1 , totv2 , totv3 , totv ,  &
             exp1 , exp2 , exp3 , term
#ifdef __PGI
      !real(rkx) , external :: erf
#endif
      !
      !! Soil texture characterisation/
      !! Based on STAT-FAO soil data types and summarized in Menut et al. 2013
      !! gives the texture composition (percent of CS,FMS,Silt, Clay, Salt)  for each of the 12 categories
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
      ! soil distribution parameters for each of the texture component, also summarized in Menut et al. 2013
      data  texmmd  / 690.0_rkx, 210.0_rkx, 125.0_rkx,2.0_rkx, 520.0_rkx /
      data  texstd  / 1.6_rkx,   1.6_rkx,   1.8_rkx,  2.0_rkx, 1.50_rkx /

      ! specific table for clay component
      fclay(:) = soiltexpc(4,:)
      fsand(:) = soiltexpc(1,:) + soiltexpc(2,:) 

      mmd = d_zero
      sigma = d_zero
      pcent = d_zero
      if ( ichdustemd <= 2 ) then
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

! read 2D texture class composition
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

      if (ichdustparam == 1) then
         call read_dust_param(erodfc, aez0 )
      else
         erodfc(jci1:jci2,ici1:ici2) = d_one
         ! This is a default value (in cm) for aeolian roughness lenght
         aez0(jci1:jci2,ici1:ici2,:) = 1.e-2_rkx
      end if

      ! read mineral fractions
      if ( nmine > 0 ) then
        call read_miner(nmine,cminer,sminer)
      end if

! calculation of srel2d
      clayrow2  = d_zero
      srel2d    = d_zero

      do i = ici1 , ici2
        do j = jci1 , jci2
          do nt = 1 , nats
            ! grid level clay fraction in percent
            clayrow2(j,i) = clayrow2(j,i) + dustsotex(j,i,nt)*fclay(nt)*100_rkx
          end do
        end do
      end do ! end j loop

      dp_array(1) = 0.0001_rkx  !cm
      do ns = 2 , nsoil
        dp_array(ns) = dp_array(ns-1)*exp(0.0460517018598807_rkx)
        deldp = dp_array(ns) - dp_array(ns-1)
      end do

      if ( nbin == 4 ) then
        ndi = ndi_4
      else
        ndi = ndi_12
      end if
      allocate(di(ndi))
      di(1) = 0.01_rkx !microm
      do ns = 2 , ndi
        di(ns) = 0.01_rkx * real(ns-1,rkx)
      end do

      do i = ici1 , ici2
        do j = jci1 , jci2
          do nt = 1 , nats
            ss(:) = d_zero
            stotal = d_zero

              do ns = 1 , nsoil          !soil size segregatoin no
                do nm = 1 , mode       ! soil mode = 5
                  if ( (pcent(nm,nt) > eps) .and. (sigma(nm,nt) > eps) ) then
                    xk = pcent(nm,nt)/(sqrt(twopi)*log(sigma(nm,nt)))
                    xl = ((log(dp_array(ns))- &
                           log(mmd(nm,nt)*1.0e-4_rkx))**2) / &
                           (d_two*(log(sigma(nm,nt)))**2)
                    if ( xl > mxarg ) then
                      xm = d_zero
                    else
                      xm = xk*exp(-xl)
                    end if
                    xn = rhodust*twot*(dp_array(ns)*d_half)
                    deldp = 0.0460517018598807_rkx
                    ss(ns) = ss(ns) + (xm*deldp/xn)
                  end if
                end do
                stotal = stotal + ss(ns)
              end do
              do ns = 1 , nsoil
                if ( stotal > d_zero ) then
                  srel2d(ns,nt,j,i) = min(ss(ns)/stotal,d_one)
                end if
              end do

          end do ! soil types
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
      else if ( size(idust) == size(dustbsiz2,1) ) then
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
      else if ( ichdustemd <= 3 ) then
        ! calculate the bin mass fraction to the total mass from Kok et al.,2011
        ! distribution (mass distribution).
        ! the independant variable is diameter so going from
        ! dV/dlnD to dV/dD implies a factor 1/D
        frac = 0._rkx
        totv = 0._rkx
        do ns = 1 , ndi
          ! kok 2011
          term = d_one/cv * (d_one+erf(log(di(ns)/d)/sqrt(d_two)/ &
                     log(sigmas)))*exp(-(di(ns)/lambda)**3)
          do n = 1 , nbin
             if ( di(ns) > dustbsiz(n,1) .and. di(ns) <= dustbsiz(n,2) ) then
                frac(n) = frac(n) + term
             end if
           end do
           totv = totv + term
        end do
        frac(:) = frac(:) / totv
        if ( abs(sum(frac) - d_one) > epsilon(1.0) ) then
          p = maxloc(frac)
          n = p(1)
          frac(n) = frac(n) + (d_one-sum(frac))
        end if
        if ( abs(sum(frac) - d_one) > 0.005 ) then
          write(stderr,*) 'TOTFRAC = ', sum(frac)
          write(stderr,*) 'ERROR = ', sum(frac)-d_one
          call fatal(__FILE__,__LINE__, &
               'Weight normalization failed')
        end if
      end if

      deallocate(di)

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
      real(rkx) , parameter :: agamma = 3.0e-4_rkx
      real(rkx) , parameter :: f = 0.0123_rkx
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
    subroutine sfflux(lmonth,ivegcov,vegfrac,snowfrac,ustarnd,z0,soilw, &
                      surfwd,roarow,trsize)
      implicit none
      integer(ik4) , intent(in) :: lmonth
      integer(ik4) , intent(in) , dimension(jci1:jci2,ici1:ici2) :: ivegcov
      real(rkx) , intent(in) , dimension(jci1:jci2,ici1:ici2) :: roarow
      real(rkx) , intent(in) , dimension(jci1:jci2,ici1:ici2) :: soilw
      real(rkx) , intent(in) , dimension(jci1:jci2,ici1:ici2) :: surfwd
      real(rkx) , intent(in) , dimension(jci1:jci2,ici1:ici2) :: vegfrac
      real(rkx) , intent(in) , dimension(jci1:jci2,ici1:ici2) :: snowfrac
      real(rkx) , intent(in) , dimension(jci1:jci2,ici1:ici2) :: z0
      real(rkx) , intent(in) , dimension(jci1:jci2,ici1:ici2) :: ustarnd
      real(rkx) , intent(in) , dimension(nbin,2) :: trsize
      real(rkx) , dimension(ilg) :: xclayrow , xroarow , xsoilw , &
                xsurfwd , xvegfrac , xz0 , xaez0 , xustarnd , xsnowfrac
      real(rkx) , dimension(ilg,nbin) :: xrsfrow
      real(rkx) , dimension(ilg,nats) :: xftex
      real(rkx) , dimension(ilg,nsoil,nats) :: xsrel2d
      integer(ik4) :: i , j , ieff , n , ns , m

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
      xaez0 = d_zero

      ieff = 0
      do i = ici1 , ici2
        do j = jci1 , jci2
!          if ( ivegcov(j,i) == 8 .or. ivegcov(j,i) == 11 ) then
            if ( vegfrac(j,i) < 0.999_rkx .and. ivegcov(j,i) > 0 ) then
            ieff = ieff + 1
            xvegfrac(ieff) = vegfrac(j,i)
            xsnowfrac(ieff) =  snowfrac(j,i)
            xsoilw(ieff) = soilw(j,i)
            xsurfwd(ieff) = surfwd(j,i)
            xz0(ieff) = z0(j,i)
            xroarow(ieff) = roarow(j,i)
            xustarnd(ieff) = ustarnd(j,i)
            xclayrow(ieff) = clayrow2(j,i)
            xaez0(ieff) = aez0(j,i,lmonth)
            do n = 1 , nats
              xftex(ieff,n) = dustsotex(j,i,n)
              do  ns = 1 , nsoil
                xsrel2d(ieff,ns,n) = srel2d(ns,n,j,i)
              end do
            end do
          end if
        end do
      end do
      if ( ieff > 0 ) then
        call dust_module(1,ieff,trsize,xsoilw,xvegfrac,xsnowfrac,xsurfwd, &
                         xftex,xclayrow,xroarow,fclay,xz0,xaez0, &
                         xsrel2d,xustarnd,xrsfrow)
      end if

      ! put back the dust flux on the right grid

      do n = 1 , nbin
        ieff = 1
        do i = ici1 , ici2
          do j = jci1 , jci2
           ! if ( ivegcov(j,i) == 8 .or. ivegcov(j,i) == 11 ) then
             if (vegfrac(j,i) < 0.999_rkx .and. ivegcov(j,i) > 0) then
              if ( ichdrdepo == 1 ) then
                if ( idynamic == 3 ) then
                  chiten(j,i,kz,idust(n)) = chiten(j,i,kz,idust(n)) + &
                       xrsfrow(ieff,n)/(cdzq(j,i,kz)*crhob3d(j,i,kz))
                else
                  ! kg m-2 s-1 => kg/kg * s-1 * psb
                  chiten(j,i,kz,idust(n)) = chiten(j,i,kz,idust(n)) + &
                       xrsfrow(ieff,n)/(cdzq(j,i,kz)*crhob3d(j,i,kz))*cpsb(j,i)
                end if
              else if ( ichdrdepo == 2 ) then
                ! pass the flux to BL scheme
                chifxuw(j,i,idust(n)) = chifxuw(j,i,idust(n)) + &
                     xrsfrow(ieff,n)
              end if
              ! diagnostic source (accumulated)
              cemtrac(j,i,idust(n)) = cemtrac(j,i,idust(n)) + &
                       xrsfrow(ieff,n)* cfdout
              if ( ichdiag > 0 ) then
                cemisdiag(j,i,kz,idust(n)) = &
                     cemisdiag(j,i,kz,idust(n)) + xrsfrow(ieff,n) / &
                     (cdzq(j,i,kz)*crhob3d(j,i,kz)) * cfdout
              end if
              ieff = ieff + 1
            end if
          end do
        end do
      end do

      ! Mineralogy flux option
      ! introduce mineralogy here, implicit loop on mineral types
      ! mines only in desert ?

      if ( nmine > 0 ) then
        do m = 1 , nmine
          do n = 1 , nbin
            ieff = 1
            do i = ici1 , ici2
              do j = jci1 , jci2
                if ( ivegcov(j,i) == 8 .or. ivegcov(j,i) == 11 ) then
                  if ( idynamic == 3 ) then
                    chiten(j,i,kz,imine(n,m)) = &
                         chiten(j,i,kz,imine(n,m)) + xrsfrow(ieff,n) * &
                         (cminer(j,i,m)*cfrac(n) + sminer(j,i,m)*sfrac(n)) / &
                         (cdzq(j,i,kz)*crhob3d(j,i,kz))
                  else
                    chiten(j,i,kz,imine(n,m)) = &
                         chiten(j,i,kz,imine(n,m)) + xrsfrow(ieff,n) * &
                         (cminer(j,i,m)*cfrac(n) + sminer(j,i,m)*sfrac(n)) / &
                         (cdzq(j,i,kz)*crhob3d(j,i,kz)) * cpsb(j,i)
                  end if
                  cemtrac(j,i,imine(n,m)) = cemtrac(j,i,imine(n,m)) + &
                         xrsfrow(ieff,n)* (cminer(j,i,m) *cfrac(n)  + &
                                       sminer(j,i,m) *sfrac(n)) * cfdout
                  ieff = ieff + 1
                end if
              end do
            end do
          end do
        end do
      end if
    end subroutine sfflux

    subroutine dust_module(jl1,jl2,trsize,soilw,vegfrac,snowfrac,surfwd,ftex, &
                           clayrow,roarow,fclay,z0,aez0,srel,ustarnd,rsfrow)
      implicit none
      integer(ik4) :: jl1 , jl2
      real(rkx) , dimension(ilg) :: clayrow , roarow , soilw , surfwd ,   &
                            vegfrac , z0 , ustarnd , snowfrac , aez0
      real(rkx) , dimension(ilg,nbin) :: rsfrow
      real(rkx) , dimension(ilg,nats) :: ftex
      real(rkx), dimension(nats)      :: fclay
      real(rkx) , dimension(ilg,nsoil,nats) :: srel
      real(rkx) , dimension(nbin,2) :: trsize
      intent (in) clayrow , soilw , surfwd , z0 , ustarnd , ftex

      real(rkx) , dimension(ilg) :: hc , rc , srl , wprim
      real(rkx) :: cly1 , cly2 , tempd , ustarns , uth , utmin
      integer(ik4) :: j,n
      real(rkx) , dimension(ilg) :: ustar
      real(rkx) , dimension(ilg,nsoil) :: utheff

      real(rkx) , parameter :: umin = 15.0_rkx
      real(rkx) , parameter :: xz = 0.25_rkx
      real(rkx) , parameter :: br = 202.0_rkx
      real(rkx) , parameter :: ym = 0.16_rkx
      real(rkx) , parameter :: sigr = 1.45_rkx
      real(rkx) , parameter :: z0s = 1.0e-3_rkx
      real(rkx) , parameter :: x = d_10

     ! threshold friction velocity on dry smooth surface
      do n = 1 , nsoil
        do j = jl1 , jl2
          if ( ust == 0 ) utheff(j,n) = ustart0(rhodust,dp_array(n),roarow(j))
          if ( ust == 1 ) utheff(j,n) = ustart01(rhodust,dp_array(n),roarow(j))
       end do
      end do

      do j = jl1 , jl2
        srl(j) = z0(j)*d_100
        rc(j) = d_one
        ! Marticorena et al., 1997: thershold friction velocity correction
        ! factor for non erodible elements. aez0 , aeolian roughness
        if ( aez0(j) <= d_zero ) then
          rc(j) = d_zero
        else if ( aez0(j) > z0s ) then
          rc(j) = d_one - (log(aez0(j)/z0s) / &
                  (log(0.35_rkx*(x/z0s)**0.8_rkx)))
        else
          rc(j) = d_one
        end if
        ! threshold velocity correction for soil humidity hc
        ! based on Fécan et al. 1999. Grid level.
        if (jsoilm == 1) then
          cly1 = clayrow(j)
          cly2 = cly1*cly1
          wprim(j) = 0.0014_rkx*cly2 + 0.17_rkx*cly1
          tempd =  max(0.00001_rkx,soilw(j)*d_100 -wprim(j))
!         print*,'humidity',i,cly1,soilw(j)*100,wprim(j),tempd
          if ( soilw(j)*d_100 > wprim(j) ) then
            hc(j) = sqrt(d_one+1.21_rkx*tempd**0.68_rkx)
!           print*,'hc',i,hc(j)
          else
            hc(j) = d_one
          end if
          ! no soil humidity corretion
        else
          hc(j) = d_one
        end if
        ! * total correction factor for both hc and rc
        rc(j) = rc(j)/hc(j)

        ! * computation of the wind friction velocity
        ! * accounting for the increase of the roughness length
        ! * due to the saltation layer (gillette etal. jgr 103,
        ! * no. d6, p6203-6209, 1998
        ustarns = ustarnd(j)*d_100 !cm.s-1
        utmin = (umin/(d_100*vonkar*rc(j)))*log(d_1000/srl(j))
        if (surfwd(j) >= utmin ) then
          ustar(j) = ustarns + 0.3_rkx*(surfwd(j)-utmin)**2
        else
          ustar(j) = ustarns
        end if
      end do       ! end i loop
      call emission(jl1,jl2,rhodust,ftex,fclay,uth,roarow,rc,utheff, &
                    ustar,srel,rsfrow,vegfrac,snowfrac)

    end subroutine dust_module

    subroutine uthefft(jl1,jl2,ust,nsoil,roarow,utheff,rhodust)
      implicit none
      integer(ik4) :: jl1 , jl2 , nsoil , ust
      real(rkx) :: rhodust
      real(rkx) , dimension(ilg) :: roarow
      real(rkx) , dimension(ilg,nsoil) :: utheff
      intent (in) jl1 , jl2 , nsoil , ust
      intent (out) utheff
      integer(ik4) :: n , j
      do n = 1 , nsoil
        do j = jl1 , jl2
          if ( ust == 0 ) utheff(j,n) = ustart0(rhodust,dp_array(n),roarow(j))
          if ( ust == 1 ) utheff(j,n) = ustart01(rhodust,dp_array(n),roarow(j))
       end do
      end do
    end subroutine uthefft

    subroutine emission(jl1,jl2,rhodust,ftex,fclay,uth,roarow,rc, &
                        utheff,ustar,srel,rsfrow,vegfrac,snowfrac)
      implicit none
      integer(ik4) :: jl1 , jl2
      real(rkx) :: rhodust , uth
      real(rkx) , dimension(ilg) :: rc , ustar, roarow , vegfrac , snowfrac
      real(rkx) , dimension(ilg,nbin) :: rsfrow
      real(rkx) , dimension(ilg,nats) :: ftex
      real(rkx) , dimension(nats) :: fclay
      real(rkx) , dimension(ilg,nsoil,nats) :: srel
      real(rkx) , dimension(ilg,nsoil) :: utheff
      intent (in)  jl1 , jl2 , rc , rhodust , roarow , srel ,  &
                   ustar , utheff , vegfrac, ftex
      intent (inout) rsfrow , uth
      real(rkx) :: p1 , p2 , p3 , dec , ec , fdp1 , fdp2
      real(rkx) , dimension(ilg,nats) :: fsoil , fsoil1 , fsoil2 , fsoil3
      integer(ik4) :: j , k , n , nt , ns
      real(rkx), dimension(ilg,nbin,nats):: rsfrowt

      real(rkx) , parameter :: beta = 16300.0_rkx

      ! Parameters used in K 2014
      real(rkx) :: Cd,k1,k2,utheffc, usst, ustark
      !
      real(rkx) , parameter :: roa0 = 1.225_rkx
      real(rkx) , parameter :: usst0 = 0.16_rkx
      real(rkx) , parameter :: calph = 2.7_rkx
      real(rkx) , parameter :: Ce = 2.0_rkx
      real(rkx) , parameter :: Cd0 = 4.5e-5_rkx
      !

      fsoil(:,:) = d_zero
      if ( ichdustemd == 1 ) then

        p1 = d_zero
        p2 = d_zero
        p3 = d_zero

        fsoil1(:,:) = d_zero
        fsoil2(:,:) = d_zero
        fsoil3(:,:) = d_zero

        do nt = 1 , nats
          do j = jl1 , jl2
            if ( ftex(j,nt) < 1.e-10_rkx ) cycle
            do ns = 1 , nsoil
              if ( rc(j) > d_zero .and. ustar(j) > d_zero ) then
                uth = utheff(j,ns)/(rc(j)*ustar(j))
                if ( uth <= d_one ) then
                  fdp1 = ustar(j)**3*(d_one-uth*uth)
                  fdp2 = (d_one+uth)*rdstemfac*(1.0e-5_rkx)*roarow(j)*regrav
                  if ( fdp2 <= d_zero ) fdp2 = d_zero
                  ! FAB: with subgrid soil texture, the aggregation of vertical
                  ! fluxes per texture type at the grid cell level is done in
                  ! fine.

                  fsoil(j,nt) = srel(j,ns,nt)*fdp1*fdp2
                  ! size-distributed kinetic energy flux(per texture type)
                  dec = fsoil(j,nt)*beta
                  ! individual kinetic energy for an aggregate of size dp (
                  ! g cm2 s-2) cf alfaro (dp) is in cm
                  ec = (mathpi/12.0_rkx)*rhodust*1.0e-3_rkx * &
                      (dp_array(ns)**3)*(20.0_rkx*ustar(j))**2
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
                  fsoil1(j,nt) = fsoil1(j,nt) + 1.0e-2_rkx*p1*(dec/e1)* &
                              (mathpi/6.0_rkx)*rhodust*((d1*1.0e-4_rkx)**3)
                  fsoil2(j,nt) = fsoil2(j,nt) + 1.0e-2_rkx*p2*(dec/e2)* &
                              (mathpi/6.0_rkx)*rhodust*((d2*1.0e-4_rkx)**3)
                  fsoil3(j,nt) = fsoil3(j,nt) + 1.0e-2_rkx*p3*(dec/e3)* &
                              (mathpi/6.0_rkx)*rhodust*((d3*1.0e-4_rkx)**3)


                end if
              end if
            end do
          end do
        end do
      end if

      if ( ichdustemd == 2 ) then
        ! Kok et al., ACP, 2014 parameterisation
        ! NB : the integration on soil aggregate size distribution is
        ! made explicitely using MB95. This differs from clm5 approach
        do nt = 1 , nats
          do j = jl1 , jl2
            if ( ftex(j,nt) < 1.e-10_rkx ) cycle
            if ( rc(j) > d_zero .and. ustar(j) /= d_zero ) then
              do ns = 1 , nsoil
                  utheffc = utheff(j,ns)/ rc(j)
                  uth =  utheffc / ustar(j)
                  usst =  utheffc * (roarow(j)/roa0)**0.5
                  usst = usst / 100._rkx  ! usst in m.s-1
                  utheffc = utheffc / 100._rkx ! utheffc in m.s-1
                  ustark = ustar(j) / 100._rkx !
                  if ( uth <= d_one ) then
                    k1 = calph * (usst - usst0) / usst0
                    k2 = roarow(j) * (ustark**2 - utheffc**2)/ usst
                    ! Cd is equivalent to erodibilty
                    ! maybe output it in the future
                    Cd = Cd0 * exp(-Ce * (usst - usst0)/ usst0)
                    !
                    ! finally integrate over nsoil
                    ! note the clay fraction is used
                    fsoil(j,nt) = fsoil(j,nt) +  &
                                  srel(j,ns,nt)* &
                                  fclay(nt) * Cd * &
                                  k1 * uth ** k2 * rdstemfac
                  end if
              end do
            end if
          end do
        end do

      end if ! end Kok 2014

      ! calculate fluxes for each of transport bins
      !
      rsfrowt(:,:,:) = d_zero
      if ( ichdustemd == 1 ) then
        do nt = 1 , nats
          do n = 1 , nbin
            do j = jl1 , jl2
              rsfrowt(j,n,nt) = fsoil1(j,nt)*frac1(n) + &
                                fsoil2(j,nt)*frac2(n) + &
                                fsoil3(j,nt)*frac3(n)
            end do
          end do
        end do
      else if ( ichdustemd == 2 ) then
        do nt = 1 , nats
          do n = 1 , nbin
            do j = jl1 , jl2
              rsfrowt(j,n,nt) = fsoil(j,nt)*frac(n)
            end do
          end do
        end do
      end if
      !
      ! Finally, aggregation of the dust flux at the grid cell level =
      ! weighted sum over soil texture
      ! weighting by grid cell veg fraction, snow fraction,
      !
      ! f = d_zero
      !

      rsfrow(:,:) = d_zero
      do k = 1 , nbin
        do nt = 1 , nats
          do j = jl1 , jl2
#ifdef CLM45
            ! CLM45 sends fraction of ground emitting dust
            rsfrow(j,k) = rsfrow(j,k) + rsfrowt(j,k,nt)*ftex(j,nt) * &
                          (d_one - vegfrac(j))
#else
            rsfrow(j,k) = rsfrow(j,k) + rsfrowt(j,k,nt)*ftex(j,nt) * &
                          (d_one - vegfrac(j))*(d_one - snowfrac(j))
#endif

            ! * (1-snowfrac)
          end do
        end do
      end do
    end subroutine emission

    subroutine clm_dust_tend
      implicit none
#ifdef CLM45
      integer(ik4) :: i , j , n , ib
      ! real(rkx) :: cdsfrq
      ! Update dust tendency with dust fluxes calculated in CLM
      ! here sump up the total flux from clm ( initially defined on 4 bins)
      ! and re-distribute it according to the selected dust emission size
      ! distrib frac(n) is the weight relatibe to bin n ..
      ! this insure consistency between emission distrbution and effective
      ! readius and optical properties
      ! use the same tuning  rdstemfac and also add the erodibility source function (fixed to 1 in CLM)
      ! from external data.

      do i = ici1 , ici2
        do j = jci1 , jci2
          sumdflux(j,i) = sum(cdustflx_clm(j,i,:)) * rdstemfac * erodfc(j,i)
        end do
      end do

      do n = 1 , nbin
        ib = idust(n)
        if ( ichdrdepo == 1 ) then
          if ( idynamic == 3 ) then
            do i = ici1 , ici2
              do j = jci1 , jci2
                chiten(j,i,kz,ib) = chiten(j,i,kz,ib) + &
                     sumdflux(j,i)*frac(n)/(cdzq(j,i,kz)*crhob3d(j,i,kz))
              end do
            end do
          else
            do i = ici1 , ici2
              do j = jci1 , jci2
                chiten(j,i,kz,ib) = chiten(j,i,kz,ib) + &
                     sumdflux(j,i)*frac(n)/(cdzq(j,i,kz)*crhob3d(j,i,kz)) * &
                     cpsb(j,i)
              end do
            end do
          end if
        else if ( ichdrdepo == 2 ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              ! pass the flux to BL scheme
              chifxuw(j,i,ib) = chifxuw(j,i,ib) + sumdflux(j,i) * frac(n)
            end do
          end do
        end if
        do i = ici1 , ici2
          do j = jci1 , jci2
            ! diagnostic source (accumulated)
            ! cdsfrq = cfdout
            cemtrac(j,i,ib) = cemtrac(j,i,ib) + sumdflux(j,i)*frac(n) * cfdout
          end do
        end do
        if ( ichdiag > 0 ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              cemisdiag(j,i,kz,ib) = cemisdiag(j,i,kz,ib) + &
                                 sumdflux(j,i)*frac(n) / &
                                 (cdzq(j,i,kz)*crhob3d(j,i,kz))*cfdout
            end do
          end do
        end if
      end do
#endif
    end subroutine clm_dust_tend

end module mod_che_dust

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

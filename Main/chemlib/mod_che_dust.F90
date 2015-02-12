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

  real(rk8) , dimension(4,2) :: dustbsiz1
  real(rk8) , dimension(12,2) :: dustbsiz2

  ! Fix the actual dust aerosol bin size: diameter in microm

!  data  dustbsiz1 / 0.01D0,  1.00D0,  2.50D0,  5.00D0,  1.00D0, &
!                   2.50D0,  5.00D0, 20.00D0/

  data  dustbsiz1 / 0.09D0,  1.00D0,  2.50D0,  5.00D0,  1.00D0, &
                   2.50D0,  5.00D0, 63.00D0/


  ! new option defined from LISA optimized distribution
  data  dustbsiz2 /0.09D0,0.18D0, 0.60D0, 1.55D0, 2.50D0, 3.75D0,   &
                   4.70D0, 5.70D0, 7.50D0, 14.50D0, 26.0D0, 41.0D0, &
                   0.18D0,0.60D0, 1.55D0, 2.50D0, 3.75D0, 4.70D0,   &
                   5.70D0, 7.50D0, 14.50D0, 26.0D0, 41.0D0, 63.0D0 /

  ! dust effective diameter

  real(rk8) , dimension(4)    ::  dustbed1
  real(rk8) , dimension(12)   ::  dustbed2

  ! has to be calculated from an assumed sub-bin distribution
!  data dustbed1 /0.82D0 , 1.8D0 , 3.7D0 , 12.5D0 /

!FAB
! if calculated from Kok distribution 
  data dustbed1 /0.658184, 1.75093D0, 3.67936D0, 8.46347D0 /

  ! 12 bins option calculated from Kok Distibution
  data dustbed2 / 0.14062217D0, 0.4300415D0,  1.10404692D0,   &
                  2.02166018D0, 3.10952699D0, 4.21185749D0,   &
                  5.18438211D0, 6.55182088D0, 9.96016755D0,   &
                 16.19150734D0, 26.74151275D0,  41.32554485D0 /

  ! solubility of od dust aer for param of giorgi and chameides

  real(rk8) , dimension(4) ::  soldust1
  real(rk8) , dimension(12) ::  soldust2

  data  soldust1 /0.1D0 , 0.1D0 , 0.1D0 , 0.1D0 /
  data  soldust2 /0.1D0 , 0.1D0 , 0.1D0 , 0.1D0, &
                  0.1D0 , 0.1D0 , 0.1D0 , 0.1D0, &
                  0.1D0 , 0.1D0 , 0.1D0 , 0.1D0 /

  ! Basic dust aerosol density (ACE-2 ) in kg/m3
  real(rk8) , parameter :: rhodust = 2650.0D0

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
  real(rk8) , parameter :: d1 = 1.5D0
  real(rk8) , parameter :: d2 = 6.7D0
  real(rk8) , parameter :: d3 = 14.2D0
  real(rk8) , parameter :: sigma1 = 1.7D0
  real(rk8) , parameter :: sigma2 = 1.2D0
  real(rk8) , parameter :: sigma3 = 1.5D0

  ! corresponding binding energies/
  real(rk8) , parameter  :: e1 = 3.61D0
  real(rk8) , parameter  :: e2 = 3.52D0
  real(rk8) , parameter  :: e3 = 3.46D0
  !
  ! parameter for alternative Kok emission distribution
  real(rk8) , parameter :: d = 3.4D0
  real(rk8) , parameter :: sigmas = 3.0D0
  ! Normalization constant
  real(rk8) , parameter :: cv = 12.62D0
  real(rk8) , parameter :: lambda = 12.0D0

  !FENNEC distribution parameters(Ryder et. al. 2013)
  real(rk8) , parameter :: d1F = 0.05D0
  real(rk8) , parameter :: d2F = 0.71D0
  real(rk8) , parameter :: d3F = 2.04D0
  real(rk8) , parameter :: d4F = 5.28D0

  real(rk8) , parameter :: sigma1F = 2.5D0
  real(rk8) , parameter :: sigma2F = 1.33D0
  real(rk8) , parameter :: sigma3F = 1.45D0
  real(rk8) , parameter :: sigma4F = 2.0D0

  real(rk8) , parameter :: N1 = 508.27D0
  real(rk8) , parameter :: N2 = 8.84D0
  real(rk8) , parameter :: N3 = 1.89D0
  real(rk8) , parameter :: N4 = 0.54D0
  ! soil variable, srel 2 d corresponds to the soil aggregate size distribution
  ! in each texture type.
  real(rk8) , pointer,  dimension(:) :: dustbed , soldust
  real(rk8) , pointer,  dimension(:) :: frac1 , frac2 , frac3 , frac
  real(rk8) , pointer , dimension(:,:,:) :: clay2row2 , sand2row2 , silt2row2
  real(rk8) , pointer , dimension(:,:) :: clayrow2 , sandrow2 , dustbsiz
  real(rk8) , pointer , dimension(:,:,:,:) :: srel2d
  real(rk8) , pointer , dimension(:,:,:) :: dustsotex
  ! Name of variable changed ! SC. 06.10.2010
  real(rk8) , dimension(nsoil) :: dp_array

  public :: sandrow2
  public :: rhodust
  public :: soldust , dustbed , dustbsiz

  integer(ik4) :: ilg

  public :: allocate_mod_che_dust , inidust , sfflux

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
      real(rk8) , dimension(nats) :: bcly , bslt , bsnd
      real(rk8) :: deldp , eps , stotal , xk , xl , xm , xn
      integer(ik4) :: i , j , n , nm , ns , nt , itr
      real(rk8) , dimension(mode,nats) :: mmdd , pcentd , sigmad
      real(rk8) , dimension(mode,nats) :: mmd , pcent , sigma
      real(rk8) , dimension(iy,nsoil,nats) :: srel
      real(rk8) , dimension(nsoil) :: ss
      real(rk8) , dimension(ndi) :: di
      ! modif new distribution
      ! for each category, this is the percent of Coarse sand,
      ! Fine mode sand, silt , clay and salt ( cf Menut et al. ,2012)
      real(rk8) , dimension (mode,12) :: soiltexpc
      real(rk8) , dimension (mode)    :: texmmd , texstd

      logical :: rd_tex
      character(6) :: aerctl
      real(rk8) :: alogdi , amean1 , amean2 , amean3 , asigma1 , amean , &
             asigma , asigma2 , asigma3 , rwi , totv1 , totv2 , totv3 , totv
      ! Fennec
      real(rk8) :: amean1F , amean2F , amean3F , amean4F , asigma1F , &
                   asigma2F , asigma3F , asigma4F, totvF, mass_dist
#ifdef __PGI
      real(rk8) , external :: derf
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
      ! data bcly/0.00D0 , 0.4D-2 ,0.7D-2  , 0.7D-2 , 0.4D-2 , 1.D-2 , &
      !           3.D-2 , 3D-2 , 5.D-2 , 8.D-2 , 8.D-2 , 1.D-2/

      ! data bcly / 4.3D-2, 2.3D-2, 7.3D-2, 0.0D0,0.0D0,0.0D-2,0.0D0, &
      !             0.0D0,0.0D0,0.0D0,0.0D0,0.0D0/
      data bcly / 6.D-2, 2.3D-2, 7.3D-2, 0.0D0,0.0D0,0.0D-2,0.0D0, &
                  0.0D0,0.0D0,0.0D0,0.0D0,0.0D0/

      ! bsnd and bslt are not really used after /
      ! the data here are not consistent with clay.
      data bsnd/0.90D0 , 0.85D0 , 0.80D0 , 0.50D0 , 0.45D0 , 0.35D0 , &
                0.30D0 , 0.30D0 , 0.20D0 , 0.65D0 , 0.60D0 , 0.50D0/
      data bslt/0.05D0 , 0.05D0 , 0.051D0 , 0.35D0 , 0.40D0 , 0.60D0 , &
                0.65D0 , 0.50D0 , 0.05D0 , 0.00D0 , 0.00D0 , 0.00D0/

      data eps/1.0D-7/

      data mmdd/690.0D0 ,  0.0D0 ,   0.0D0 , 0.0D0 ,   0.0D0, &
               690.0D0 , 210.0D0 ,   0.0D0 , 0.0D0 ,   0.0D0, &
               690.0D0 , 210.0D0 ,   0.0D0 , 0.0D0 ,   0.0D0, &
               520.0D0 , 100.0D0 ,   5.0D0 , 0.0D0 ,   0.0D0, &
               520.0D0 ,  75.0D0 ,   2.5D0 , 0.0D0 ,   0.0D0, &
               520.0D0 ,  75.0D0 ,   2.5D0 , 0.0D0 ,   0.0D0, &
               210.0D0 ,  75.0D0 ,   2.5D0 , 0.0D0 ,   0.0D0, &
               210.0D0 ,  50.0D0 ,   2.5D0 , 0.0D0 ,   0.0D0, &
               125.0D0 ,  50.0D0 ,   1.0D0 , 0.0D0 ,   0.0D0, &
               100.0D0 ,  10.0D0 ,   1.0D0 , 0.0D0 ,   0.0D0, &
               100.0D0 ,  10.0D0 ,   0.5D0 , 0.0D0 ,   0.0D0, &
               100.0D0 ,  10.0D0 ,   0.5D0, 0.0D0 ,   0.0D0 /

      data sigmad/1.6D0 , 1.8D0 , 1.8D0 , 0.0D0 ,   0.0D0,  &
                 1.6D0 , 1.8D0 , 1.8D0 ,  0.0D0 ,   0.0D0,  &
                 1.6D0 , 1.8D0 , 1.8D0 , 0.0D0 ,   0.0D0,   &
                 1.6D0 , 1.7D0 , 1.8D0 , 0.0D0 ,   0.0D0,   &
                 1.6D0 , 1.7D0 , 1.8D0 ,  0.0D0 ,   0.0D0,  &
                 1.6D0 , 1.7D0 , 1.8D0 , 0.0D0 ,   0.0D0,   &
                 1.7D0 , 1.7D0 , 1.8D0 , 0.0D0 ,   0.0D0,   &
                 1.7D0 , 1.7D0 , 1.8D0 , 0.0D0 ,   0.0D0,   &
                 1.7D0 , 1.7D0 , 1.8D0 , 0.0D0 ,   0.0D0,   &
                 1.8D0 , 1.8D0 , 1.8D0 ,  0.0D0 ,   0.0D0,  &
                 1.8D0 , 1.8D0 , 1.8D0 , 0.0D0 ,   0.0D0,   &
                 1.8D0 , 1.8D0 , 1.8D0,  0.0D0 ,   0.0D0 /

       data pcentd/1.00D0 , 0.00D0 , 0.00D0 ,  0.0D0 ,   0.0D0, &
                  0.90D0 , 0.10D0 , 0.00D0 ,  0.0D0 ,   0.0D0,  &
                  0.80D0 , 0.20D0 , 0.00D0 , 0.0D0 ,   0.0D0,   &
                  0.50D0 , 0.35D0 , 0.15D0 ,  0.0D0 ,   0.0D0,  &
                  0.45D0 , 0.40D0 , 0.15D0 , 0.0D0 ,   0.0D0,   &
                  0.35D0 , 0.50D0 , 0.15D0 , 0.0D0 ,   0.0D0,   &
                  0.30D0 , 0.50D0 , 0.20D0 , 0.0D0 ,   0.0D0,   &
                  0.30D0 , 0.50D0 , 0.20D0 ,  0.0D0 ,   0.0D0,  &
                  0.20D0 , 0.50D0 , 0.30D0 , 0.0D0 ,   0.0D0,   &
                  0.65D0 , 0.00D0 , 0.35D0 ,  0.0D0 ,   0.0D0,  &
                  0.60D0 , 0.00D0 , 0.40D0 , 0.0D0 ,   0.0D0,   &
                  0.50D0 , 0.00D0 , 0.50D0, 0.0D0 ,   0.0D0 /
       !!
       !! new option
       !!
       data   soiltexpc / 0.46D0, 0.46D0, 0.05D0,  0.03D0, 0.0D0, &
                          0.41D0, 0.41D0, 0.18D0,  0.00D0, 0.0D0, &
                          0.29D0, 0.29D0, 0.32D0,  0.10D0, 0.0D0, &
                          0.00D0, 0.17D0, 0.70D0,  0.13D0, 0.0D0, &
                          0.00D0, 0.10D0, 0.85D0,  0.05D0, 0.0D0, &
                          0.00D0, 0.43D0, 0.39D0,  0.18D0, 0.0D0, &
                          0.29D0, 0.29D0, 0.15D0,  0.27D0, 0.0D0, &
                          0.00D0, 0.10D0, 0.56D0,  0.34D0, 0.0D0, &
                          0.00D0, 0.32D0, 0.34D0,  0.34D0, 0.0D0, &
                          0.00D0, 0.52D0, 0.06D0,  0.42D0, 0.0D0, &
                          0.00D0, 0.06D0, 0.47D0,  0.47D0, 0.0D0, &
                          0.00D0, 0.22D0, 0.20D0,  0.58D0, 0.0D0/

      data  texmmd  / 690.0D0, 210.0D0, 125.0D0,2.0D0, 520.0D0 /
      data  texstd  / 1.6D0,   1.6D0,   1.8D0,  2.0D0, 1.50D0 /

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
      where ( dabs(soiltexpc(:,:)) < dlowval )
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

      dp_array(1) = 0.0001D0  !cm
      do ns = 2 , nsoil
        dp_array(ns) = dp_array(ns-1)*dexp(0.0460517018598807D0)
        deldp = dp_array(ns) - dp_array(ns-1)
      end do

      di(1) = 0.01D0 !microm
      do ns = 2 , ndi
        di(ns) = di(ns-1) + 0.01D0
      end do

      do j = jci1 , jci2
        srel(:,:,:) = d_zero
        do i = ici1 , ici2
          do nt = 1 , nats
            ss(:) =d_zero
            stotal = d_zero
            if ( sand2row2(i,nt,j) > d_zero ) then
              do ns = 1 , nsoil          !soil size segregatoin no
                do nm = 1 , mode       !soil mode = 5
                  if ( (pcent(nm,nt) > eps) .and.                    &
                         (sigma(nm,nt) /= d_zero) ) then
                    xk = pcent(nm,nt)/(dsqrt(twopi)*dlog(sigma(nm,nt)))
                    xl = ((dlog(dp_array(ns))- &
                           dlog(mmd(nm,nt)*1.0D-4))**2) &
                           /(d_two*(dlog(sigma(nm,nt)))**2)
                    xm = xk*dexp(-xl)
                  else
                    xm = d_zero
                  end if
                  xn = rhodust*twot*(dp_array(ns)*d_half)
                  deldp = 0.0460517018598807D0
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
      dustbsiz(:,:) = dustbsiz1(:,:)
      dustbed(:) = dustbed1(:)
      soldust(:) = soldust1(:)
      if ( chemsimtype == 'DU12' ) then
        dustbsiz(:,:) = dustbsiz2(:,:)
        dustbed(:) = dustbed2(:)
        soldust(:) = soldust2(:)
      end if

      if ( ichdustemd == 1 ) then
        totv1 = d_zero
        totv2 = d_zero
        totv3 = d_zero
        amean1 = dlog10(d1)
        amean2 = dlog10(d2)
        amean3 = dlog10(d3)
        asigma1 = dlog10(sigma1)
        asigma2 = dlog10(sigma2)
        asigma3 = dlog10(sigma3)
        do ns = 1, ndi
          alogdi = dlog10(di(ns))
          do n = 1 , nbin
            if ( di(ns) > dustbsiz(n,1) .and. &
                 di(ns) <= dustbsiz(n,2) ) then
              ! the independant variable is diameter so going from
              ! dV/dlog10D to dV/dD implies a factor 1/(2.303)D
              frac1(n) = frac1(n) + (d_one/di(ns)) * &
                dexp(-(alogdi-amean1)**2/(d_two*asigma1**2))
              frac2(n) = frac2(n) + (d_one/di(ns)) * &
                dexp(-(alogdi-amean2)**2/(d_two*asigma2**2))
              frac3(n) = frac3(n) + (d_one/di(ns)) * &
                dexp(-(alogdi-amean3)**2/(d_two*asigma3**2))
            end if
          end do
          totv1 = totv1 + (d_one/di(ns)) * &
            dexp(-(alogdi-amean1)**2/(d_two*asigma1**2))
          totv2 = totv2 + (d_one/di(ns)) * &
            dexp(-(alogdi-amean2)**2/(d_two*asigma2**2))
          totv3 = totv3 + (d_one/di(ns)) * &
            dexp(-(alogdi-amean3)**2/(d_two*asigma3**2))
        end do
        frac1(:) = frac1(:) / totv1
        frac2(:) = frac2(:) / totv2
        frac3(:) = frac3(:) / totv3
      else if ( ichdustemd == 2 ) then
        ! calculate the bin mass fraction to the total mass from Kok et al.
        ! distribution ( mass distribution).
        ! the independant variable is diameter so going from
        ! dV/dlnD to dV/dD implies a factor 1/D
        frac = 0.D0
        totv = 0.D0
        do ns = 1, ndi
          do n = 1, nbin
             if ( di(ns) > dustbsiz(n,1) .and. di(ns) <= dustbsiz(n,2) ) then
                frac(n) = frac(n) + d_one/cv * &
                  (d_one+derf(log(di(ns)/d)/sqrt(d_two)/ &
                  log(sigmas)))*exp(-(di(ns)/lambda)**3)  !see Kok (2011)
             end if
           end do
           totv = totv + d_one / cv * (d_one+derf(log(di(ns)/d)/sqrt(d_two)/ &
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
    real(rk8) function ustart01(rhodust,dum,rhair)
      implicit none
!
      real(rk8) , parameter :: a2 = 0.129D0 , c1 = 0.006D0 , c2 = 1.928D0 , &
                              c3 = 0.0858D0 , c4 = -0.0617D0 , c5 = 2.5D0 ,&
                              y1 = 1331.647D0 , y2 = 1.561228D0 ,          &
                              y3 = 0.38194D0
!
      real(rk8) , intent(in) :: dum , rhair , rhodust
      real(rk8) :: dm , rep , term , term1 , term2
!
      dm = dum  !* 1.0e-4      ! cm
      rep = y1*(dm**y2) + y3
      term1 = dsqrt(d_one+(c1/(rhodust*egrav*0.1D0*(dm**c5))))
      term2 = dsqrt(rhodust*egrav*d_100*dm/rhair)
      term = term1*term2
      ustart01 = cvmgt(a2*term*(d_one-c3*dexp(c4*(rep-d_10))),  &
                 a2*term/dsqrt(c2*(rep**0.092D0)-d_one),rep > d_10)
    contains
!
      real(rk8) function cvmgt(val1,val2,cond)
        implicit none
        logical , intent(in) :: cond
        real(rk8) , intent(in) :: val1 , val2
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
    real(rk8) function ustart0(rhodust,dum,rhoa)
      implicit none
      real(rk8) , intent(in) :: dum , rhoa , rhodust
      real(rk8) , parameter :: agamma = 3.0D-4 , f = 0.0123D0
      real(rk8) :: dm , sigma
      sigma = rhodust/rhoa
      dm = dum*1.0D-2
      ustart0 = f*(sigma*egrav*dm+agamma/(rhoa*dm))
      ustart0 = dsqrt(ustart0)
      ustart0 = ustart0*d_100
    end function ustart0
!
!   **********************************************************
!   *  dust emission scheme                             ******
!   *                                                   ******
!   * this scheme based on marticorena and bergametti,  ******
!   * 1995; gong et al.,(2003); alfaro et al.,(1997)    ******
!   * Zakey et al., 2006                                                  ******
!   **********************************************************
!
    subroutine sfflux(jloop,ivegcov,vegfrac,ustarnd,z0,soilw, &
                      surfwd,roarow,trsize,rsfrow)
!
      implicit none
!
      integer(ik4) , intent(in) :: jloop
      integer(ik4) , intent(in) , dimension(ici1:ici2) ::  ivegcov
      real(rk8) , intent(in) , dimension(ici1:ici2) :: roarow , soilw , &
                surfwd , vegfrac , z0 , ustarnd
      real(rk8) , intent(out) , dimension(ici1:ici2,nbin) :: rsfrow
      real(rk8) , intent(in) , dimension(nbin,2) :: trsize
!
      real(rk8) , dimension(ilg) :: xclayrow , xroarow , xsoilw , &
                xsurfwd , xvegfrac , xz0 , xustarnd , xsnowfrac
      real(rk8) , dimension(ilg,nbin) :: xrsfrow
      real(rk8) , dimension(ilg,nats) :: xftex , xalphaprop
      real(rk8) , dimension(ilg,nsoil,nats) :: xsrel2d
      integer(ik4) :: i , ieff , n , ns
!
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
                xalphaprop(ieff,n) = d_10**(0.134D0 * &
                              clay2row2(i,n,jloop)-6.0D0)
!                             clay2row2(i,n,jloop)-6.0D0)*0.035D0
              else
                xalphaprop(ieff,n) = d_10**(-0.1D0 * &
                              clay2row2(i,n,jloop)-6.0D0)
!                             clay2row2(i,n,jloop)-1.2D0)*0.035D0
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
                   rsfrow(i,n)*egrav/(dsigma(kz)*1.D3)
            elseif ( ichdrdepo == 2 ) then
              ! pass the flux to BL scheme
              chifxuw(jloop,i,idust(n)) = chifxuw(jloop,i,idust(n)) + &
                   rsfrow(i,n)
            end if
            ! diagnostic source (accumulated)
            cemtrac(jloop,i,idust(n)) = cemtrac(jloop,i,idust(n)) + &
                     rsfrow(i,n)* cfdout

             if ( ichdiag == 1 ) then
             cemisdiag(jloop,i,kz,idust(n)) = cemisdiag(jloop,i,kz,idust(n)) + &
                       rsfrow(i,n)/ ( cdzq(jloop,i,kz)*crhob3d(jloop,i,kz)) * cfdout
             end if
          end do
          ieff = ieff + 1
        end if
      end do

    end subroutine sfflux
!
    subroutine dust_module(il1,il2,trsize,soilw,vegfrac,snowfrac,surfwd,ftex, &
                           clayrow,roarow,alphaprop,z0,srel,ustarnd,rsfrow)
      implicit none
!
      integer(ik4) :: il1 , il2
      real(rk8) , dimension(ilg) :: clayrow , roarow , soilw , surfwd ,   &
                                   vegfrac , z0 , ustarnd , snowfrac
      real(rk8) , dimension(ilg,nbin) :: rsfrow
      real(rk8) , dimension(ilg,nats) :: ftex , alphaprop
      real(rk8) , dimension(ilg,nsoil,nats) :: srel
      real(rk8) , dimension(nbin,2) :: trsize
      intent (in) clayrow , soilw , surfwd , z0 , ustarnd , ftex
!
      real(rk8) , dimension(ilg) :: alamda , hc , rc , srl , wprim
      real(rk8) :: arc1 , arc2 , br , cly1 , cly2 , sigr , tempd ,   &
          umin , ustarns , uth , utmin , x , xz , ym , z0s , ustarfw
      integer(ik4) :: i
      real(rk8) , dimension(ilg) :: ustar
      real(rk8) , dimension(ilg,nsoil) :: utheff
!
      data umin/15.0D0/
      data xz/0.25D0/ , br/202.0D0/ , ym/0.16D0/ , sigr/1.45D0/
      data z0s/3.0D-3/ , x/d_10/

      do i = il1 , il2

        srl(i) = z0(i)*d_100
        rc(i) = d_one

        if ( jfs == 0 ) then
          ! * raupach et al. (1993)
          if ( vegfrac(i) < d_one ) then
            alamda(i) = xz*(dlog(d_one-vegfrac(i)))*(-d_one)
            arc1 = sigr*ym*alamda(i)
            arc2 = br*ym*alamda(i)
            if ( arc1 <= d_one .and. arc2 <= d_one ) then
              rc(i) = (dsqrt(d_one-arc1)*dsqrt(d_one+arc2))
            end if
          end if
        else if ( jfs == 1 ) then
          ! Marticorena et al., 1997: correction factor for non
          ! erodible elements
          rc(i) = d_one - (dlog(0.5D-2/z0s)/(dlog(0.35D0*(x/z0s)**0.8D0)))
        end if
        ! threshold velocity correction for soil humidity hc
        if ( jsoilm == 0 ) then
          if ( soilw(i) < d_zero ) then
            write(stderr,*) 'hc, rc = ' , soilw(i) , ' less than zero'
            call fatal(__FILE__,__LINE__,'NEGATIVE SOILW')
          else if ( soilw(i) < 0.03D0 ) then
            hc(i) = dexp(22.7D0*soilw(i))
          else if ( soilw(i) >= 0.03D0 ) then
            hc(i) = dexp(95.3D0*soilw(i)-2.029D0)
          else
            hc(i) = d_one
          end if
        else if ( jsoilm == 1 ) then
          cly1 = clayrow(i)
          cly2 = cly1*cly1
          wprim(i) = 0.0014D0*cly2 + 0.17D0*cly1
          tempd =  dmax1(0.00001D0,soilw(i)*d_100 -wprim(i))
!         print*,'humidity',i,cly1,soilw(i)*100,wprim(i),tempd
          if ( soilw(i)*d_100 > wprim(i) ) then
            hc(i) = dsqrt(d_one+1.21D0*tempd**0.68D0)
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
        ustarfw = (vonkar*100.0D0*surfwd(i))/(dlog(1000.0D0/srl(i)))
        ustarns = ustarnd(i)*d_100 !cm.s-1
        utmin = (umin/(d_100*vonkar*rc(i)))*dlog(d_1000/srl(i))
        if ( surfwd(i) >= utmin ) then
          ustar(i) = ustarns + 0.3D0*(surfwd(i)-utmin)*(surfwd(i)-utmin)
        else
          ustar(i) = ustarns
        end if
      end do       ! end i loop

      call uthefft(il1,il2,ust,nsoil,roarow,utheff,rhodust)

      call emission(il1,il2,rhodust,ftex,alphaprop, uth,roarow,rc,utheff, &
                    ustar,srel,rsfrow,trsize,vegfrac,snowfrac)

    end subroutine dust_module
!
    subroutine uthefft(il1,il2,ust,nsoil,roarow,utheff,rhodust)
      implicit none
      integer(ik4) :: il1 , il2 , nsoil , ust
      real(rk8) :: rhodust
      real(rk8) , dimension(ilg) :: roarow
      real(rk8) , dimension(ilg,nsoil) :: utheff
      intent (in) il1 , il2 , nsoil , ust
      intent (out) utheff
      integer(ik4) :: n , i
      do n = 1 , nsoil
        do i = il1 , il2
          if ( ust == 0 ) utheff(i,n) = ustart0(rhodust,dp_array(n),roarow(i))
          if ( ust == 1 ) utheff(i,n) = ustart01(rhodust,dp_array(n),roarow(i))
        end do
      end do
    end subroutine uthefft
!
    subroutine emission(il1,il2,rhodust,ftex,alphaprop,uth,roarow,rc, &
                        utheff,ustar,srel,rsfrow,trsize,vegfrac,snowfrac)

      implicit none
!
      integer(ik4) :: il1 , il2
      real(rk8) :: rhodust , uth
      real(rk8) , dimension(ilg) :: rc ,ustar, roarow , vegfrac , snowfrac
      real(rk8) , dimension(ilg,nbin) :: rsfrow
      real(rk8) , dimension(ilg,nats) :: ftex , alphaprop
      real(rk8) , dimension(ilg,nsoil,nats) :: srel
      real(rk8) , dimension(nbin,2) :: trsize
      real(rk8) , dimension(ilg,nsoil) :: utheff
      intent (in)  il1 , il2 , rc , rhodust , roarow , srel ,  &
                   trsize , ustar , utheff , vegfrac, ftex
      intent (inout) rsfrow , uth
!
!     real(rk8) :: beffect
      real(rk8) :: beta , p1 , p2 , p3 , rwi , dec , ec , fdp1 , fdp2
      real(rk8) , dimension(ilg,nats) :: fsoil , fsoil1 , fsoil2 , fsoil3
      integer(ik4) :: i , k , n , nt , ns

      real(rk8), dimension(ilg,nbin,nats):: rsfrowt

      data beta  /16300.0D0/

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
           if (ftex(i,nt) < 1.D-10) cycle
           do ns = 1 , nsoil
            if ( rc(i) > d_zero .and. ustar(i) /= d_zero ) then
              uth = utheff(i,ns)/(rc(i)*ustar(i))
              if ( uth <= d_one ) then
                fdp1 = ustar(i)**3*(d_one-uth*uth)
                fdp2 = (d_one+uth)*rdstemfac*(1.0D-5)*roarow(i)*regrav
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
                  ec = (mathpi/12.0D0)*rhodust*1.0D-3*(dp_array(ns)**3)* &
                        (20.0D0*ustar(i))**2
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
                  fsoil1(i,nt) = fsoil1(i,nt) + 1.0D-2*p1*(dec/e1)* &
                            (mathpi/6.0D0)*rhodust*((d1*1.0D-04)**3)
                  fsoil2(i,nt) = fsoil2(i,nt) + 1.0D-2*p2*(dec/e2)* &
                            (mathpi/6.0D0)*rhodust*((d2*1.0D-04)**3)
                  fsoil3(i,nt) = fsoil3(i,nt) + 1.0D-2*p3*(dec/e3)* &
                            (mathpi/6.0D0)*rhodust*((d3*1.0D-04)**3)
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

end module mod_che_dust

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

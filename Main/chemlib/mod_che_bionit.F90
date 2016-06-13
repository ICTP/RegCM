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

module mod_che_bionit
  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_che_common
  use mod_che_dust
  use mod_che_ncio
  use mod_dynparam

  implicit none

  public :: allocate_mod_che_bionit, ini_bionit

  real(rkx) , pointer, dimension(:,:) :: nmanure , nfert , soilph

  !coefficients and weights derived from neural net analysis
  real(rkx), parameter :: xwgt0 =  0.561651794427011_rkx
  real(rkx), parameter :: xwgt1 = -0.48932473825312_rkx
  real(rkx), parameter :: xwgt2 = -0.73521035872982_rkx
  real(rkx), parameter :: xwgt3 =  0.506600069632212_rkx
  real(rkx), parameter :: xwgt4 = -0.784867014304196_rkx
  real(rkx), parameter :: xwgt5 = -0.283241716518431_rkx
  real(rkx), parameter :: xwgt6 =  0.132539461337082_rkx
  real(rkx), parameter :: xwgt7 = -0.00839661549597702_rkx
  real(rkx), parameter :: xwgt8 = -1.62075908632141_rkx
  real(rkx), parameter :: xwgt9 =  0.638173941854311_rkx
  real(rkx), parameter :: xwgt10 =  3.88469485689393_rkx
  real(rkx), parameter :: xwgt11 = -0.942985468044301_rkx
  real(rkx), parameter :: xwgt12 = -0.862455616914003_rkx
  real(rkx), parameter :: xwgt13 = -2.68040193699105_rkx
  real(rkx), parameter :: xwgt14 =  1.61126351888328_rkx
  real(rkx), parameter :: xwgt15 =  0.134088164903734_rkx
  real(rkx), parameter :: xwgt16 = -0.21261983875851_rkx
  real(rkx), parameter :: xwgt17 =  0.901773966331639_rkx
  real(rkx), parameter :: xwgt18 = -5.18779902340853_rkx
  real(rkx), parameter :: xwgt19 =  1.23132977162784_rkx
  real(rkx), parameter :: xwgt20 = -2.62451302093078_rkx
  real(rkx), parameter :: xwgt21 = -0.27778477919531_rkx
  real(rkx), parameter :: xwgt22 = 0.413060247967231_rkx
  real(rkx), parameter :: xwgt23 = -0.560462552556124_rkx
  real(rkx), parameter :: xwgt24 = 0.499562769416134_rkx
  real(rkx), parameter :: xwgt25 = -1.23876483956298_rkx
  real(rkx), parameter :: xwgt26 = -1.41295235373665_rkx
  real(rkx), parameter :: xwgt27 = -1.20659105237301_rkx

  real(rkx), parameter :: xcoef1 = -2.453992_rkx
  real(rkx), parameter :: xcoef2 =  0.142680_rkx
  real(rkx), parameter :: xcoef3 = -4.609693_rkx
  real(rkx), parameter :: xcoef4 =  0.115964_rkx
  real(rkx), parameter :: xcoef5 = -2.717366_rkx
  real(rkx), parameter :: xcoef6 =  0.163039_rkx
  real(rkx), parameter :: xcoef7 = -0.364632_rkx
  real(rkx), parameter :: xcoef8 =  5.577532_rkx
  real(rkx), parameter :: xcoef9 = -1.535199_rkx
  real(rkx), parameter :: xcoef10 = 0.054909_rkx
  real(rkx), parameter :: xcoef11 = -25.554238_rkx
  real(rkx), parameter :: xcoef12 =   3.158129_rkx
  real(rkx), parameter :: xcoef13 =  -1.182905_rkx
  real(rkx), parameter :: xcoef14 =   0.614317_rkx
  real(rkx), parameter :: xcoef15 =   3.403007_rkx
  real(rkx), parameter :: xcoef16s =  9.205080_rkx
  real(rkx), parameter :: xcoef16l =  2.205080_rkx

contains

  !Subroutine originally based on the SOILEMISNO_n code of C.Delon (for
  !  the ISBA model) ...adapted for RegCM4-chem by F.Tummon (July 2011)
  !
  ! Calculates NOx and NH3 emissions from soil manure and fertiliser
  ! application rates
  !  ...input N application rates from Potter et. al 2010
  !
  ! Parameterises NOx emissions as a function of:
  !     - soil temperature    |
  !     - 10m winds           |read-in from regcm
  !     - soil moisture       |
  !     - % sand in the soil  |
  !     - soil fertilisation rates  |
  !     - soil pH                   |values read-in from external files
  !
  ! Parameterisation from a neural network using data from obs stations:
  !       - Aurade (France)
  !       - Grignon (France)
  !       - Escompte (France)
  !       - Hombori (Niger)
  !       - Neisouma (Niger)

  subroutine allocate_mod_che_bionit
    implicit none
    if ( ichem == 1 ) then
      call getmem2d(nfert,jce1,jce2,ice1,ice2,'che_bionit:nfert')
      call getmem2d(nmanure,jce1,jce2,ice1,ice2,'che_bionit:nmanure')
      call getmem2d(soilph,jce1,jce2,ice1,ice2,'che_bionit:soilph')
    end if
  end subroutine allocate_mod_che_bionit

  subroutine ini_bionit
    implicit none
    if ( ichbion == 1 ) call read_bionem (nfert,nmanure,soilph)
  end subroutine ini_bionit

  subroutine soilnitro_emissions(j,ivegcov,wid10)
    implicit none
    integer(ik4) , intent(in) :: j
    real(rkx) , dimension(ici1:ici2) , intent(in) :: wid10
    integer(ik4) , dimension(ici1:ici2), intent(in) :: ivegcov
    ! local variables
    integer(ik4) :: i
    ! manure application rate (kg/m2/s)
!!$ real, dimension(iy,jx), intent(in) :: manrate
    ! fertiliser application rate (kg/m2/s)
!!$ real, dimension(iy,jx), intent(in) :: fertrate
    ! soilpH
!!$ real, dimension(iy,jx), intent(in) :: soilph

    real(rkx), dimension(ici1:ici2) :: &
         soiltemp_surf ,&    ! surface soil temperature (C)
         soiltemp_deep ,&    ! deep soil temperature (C)
         sandper       ,&    ! sand percentage (%)
         porewater     ,&    ! pore space water content (%)
         canred        ,&    ! canopy reduction factor
         windsp        ,&    ! wind speed (m/s)
         soilfert      ,&    ! soil fertilisation rate (kg/m2/hr)
         norm_sm       ,&    ! normalised soil moisture
         norm_wi       ,&    ! normalised wind speed
         norm_fe       ,&    ! normalised fertilisation rate
         norm_sd       ,&    ! normalised deep soil temp
         norm_ss       ,&    ! normalised surf soil temp
         norm_ph       ,&    ! normalised ph value
         norm_sa       ,&    ! normalised sand % content
         norm_no       ,&    ! normalised NO flux
         nsum1         ,&    ! normalised sum 1
         nsum2         ,&    ! normalised sum 2
         nsum3               ! normalised sum 3
    real (rkx), dimension(ici1:ici2) :: &
         man1d,     &
         fert1d,    &
         ph1d,      &
         lai_int,   &  ! 1D manure/fertiliser app. rate + pH
         totn1d,    &  ! 1D total N app. rate
         fracnox,   &  ! fraction emitted as NOx
         fracnh3,   &  ! fraction emitted as NH3
         noxflux       ! calculated soil NOx flux

    ! getting just 1D fert/man. rates
    ! this  come from external data
    man1d = nmanure(j,ici1:ici2)
    fert1d = nfert(j,ici1:ici2)
    ph1d = soilph(j,ici1:ici2)

    ! iFAB  ! put interactive LAI
    lai_int = cxlai2d(j,ici1:ici2)

    soiltemp_surf = d_zero
    soiltemp_deep = d_zero ! deep soil temperature (C)
    sandper       = d_zero
    porewater     = d_zero
    canred        = d_zero
    windsp        = d_zero
    soilfert      = d_zero
    norm_sm       = d_zero
    norm_wi       = d_zero
    norm_fe       = d_zero
    norm_sd       = d_zero
    norm_ss       = d_zero
    norm_ph       = d_zero
    norm_sa       = d_zero
    norm_no       = d_zero
    nsum1         = d_zero
    nsum2         = d_zero
    nsum3         = d_zero
    noxflux       = d_zero

    totn1d = man1d + fert1d

    do i = ici1 , ici2
      ! cycle on sea points
      if ( ivegcov(i) == 0 ) cycle
      ! getting the soil sand percentage, pH and fert rate values
      sandper(i) = sandrow2(i,j)
      ! calculating water-filled pore space from soil moisture
      porewater(i) = cssw2da(j,i)
      ! converting soil moisture from kg/m2 to (m3 water/m3 soil)
      ! divide first by depth of soil layer = 10cm = 0.1m : cdepuv = 100mm
      ! then divide by density of water = 1000kg/m3
      porewater(i) = (porewater(i)/(cdepuv(ivegcov(i)) * d_1000)) * d_r1000
      ! calculating water-filled pore space (%)
      ! coefficient of 0.45 derived from obs at Grignon(0.536),
      ! Hombori(0.4) and Escompte(0.43) = avg. 0.45
      ! in regcm/bats  this parameter would be cxmopor : consider replacing ?
      porewater(i) = (porewater(i) * 0.45_rkx) * d_100

      ! converting temperature from Kelvin to Celsius
      soiltemp_deep(i) = ctg(j,i)
      soiltemp_surf(i) = ctga(j,i)
      if ( soiltemp_deep(i) /= 0 ) then
        soiltemp_deep(i) = soiltemp_deep(i) - tzero
      end if
      if ( soiltemp_surf(i) /= 0 ) then
        soiltemp_surf(i) = soiltemp_surf(i) - tzero
      end if

      ! calculating what percentage volatilised N gets incorporated
      ! into NH3 and NOx

      fracnh3(i) = totn1d(i)*0.3_rkx
      fracnox(i) = totn1d(i)*0.7_rkx

      ! calculation of NOx flux from soil
      ! normalised centered entries

      norm_ss(i) = xcoef1  + xcoef2*soiltemp_surf(i)
      norm_sm(i) = xcoef3  + xcoef4*porewater(i)
      norm_sd(i) = xcoef5  + xcoef6*soiltemp_deep(i)
      norm_fe(i) = xcoef7  + xcoef8*fracnox(i)
      norm_sa(i) = xcoef9  + xcoef10*sandper(i)
      norm_ph(i) = xcoef11 + xcoef12*ph1d(i)
      norm_wi(i) = xcoef13 + xcoef14*wid10(i)

      ! weighted sums (coefficients from soil_nox_params)
      nsum1(i) = xwgt0 + xwgt1*norm_ss(i) &
            + xwgt2*norm_sm(i) + xwgt3*norm_sd(i) &
            + xwgt4*norm_fe(i) + xwgt5*norm_sa(i) &
            + xwgt6*norm_ph(i) + xwgt7*norm_wi(i)

      nsum2(i) = xwgt8 + xwgt9*norm_ss(i) &
            + xwgt10*norm_sm(i) + xwgt11*norm_sd(i) &
            + xwgt12*norm_fe(i) + xwgt13*norm_sa(i) &
            + xwgt14*norm_ph(i) + xwgt15*norm_wi(i)

      nsum3(i) = xwgt16 + xwgt17*norm_ss(i) &
            + xwgt18*norm_sm(i) + xwgt19*norm_sd(i)  &
            + xwgt20*norm_fe(i) + xwgt21*norm_sa(i) &
            + xwgt22*norm_ph(i) + xwgt23*norm_wi(i)

      ! hyperbolic tangent calculation
      norm_no(i) = xwgt24 + xwgt25*tanh(nsum1(i)) &
            + xwgt26*tanh(nsum2(i)) + xwgt27*tanh(nsum3(i))

      !flux calculation
      ! If sand > 50%, pulse effect, amplitude coefficient is maximum.
      ! If sand < 50%, amplitude coefficient is reduced
      !                to avoid strong emissions
      ! Sand conditions are correlated to pH values.
      if ( ph1d(i) >= d_six ) then
        noxflux(i) = xcoef15 + xcoef16s*norm_no(i)
      else if (ph1d(i) <= d_six) then
        noxflux(i) = xcoef15 + xcoef16l*norm_no(i)
      end if

      !avoiding negative fluxes
      if ( noxflux(i) < d_zero ) then
        noxflux(i)= d_zero
      end if

      ! converting the NO flux from gN/ha/d to kg/m2/s
      ! g to kg: /1000
      ! ha to m2: /100 /100
      ! d to s: /86400
      noxflux(i) = noxflux(i)/ (1000._rkx*100._rkx*100._rkx*86400._rkx)

      ! flux reduction because of canopy absorption
      if ( lai_int(i) > 1.9_rkx .and. lai_int(i) <  5._rkx ) then
        canred(i) = 0.5_rkx
      else if ( lai_int(i) > 5._rkx ) then
        canred(i) = 0.2_rkx
      else
        canred(i) = d_one
      end if

      noxflux(i) = noxflux(i)*canred(i)

    end do

    !if ( j == 25 ) write(stdout,*) maxval(noxflux)
    !
    ! update tendency for NO flux
    do i  = ici1 , ici2
      if ( ivegcov(i) == 0 ) cycle
      if ( ichdrdepo == 1 ) then
        chiten(j,i,kz,ino) = chiten(j,i,kz,ino) + &
               noxflux(i)*egrav/(dsigma(kz)*1.e3_rkx)
      else if ( ichdrdepo == 2 ) then
        ! pass the flux to BL scheme
        chifxuw(j,i,ino) = chifxuw(j,i,ino) + noxflux(i)
      end if
      ! diagnostic source (accumulated)
      cemtrac(j,i,ino) = cemtrac(j,i,ino) + noxflux(i)* cfdout

      if ( ichdiag == 1 ) then
        cemisdiag(j,i,ino) = cemisdiag(j,i,ino) + &
               noxflux(i)/ ( cdzq(j,i,kz)*crhob3d(j,i,kz)) * cfdout
      end if
    end do
  end subroutine soilnitro_emissions

end module mod_che_bionit
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

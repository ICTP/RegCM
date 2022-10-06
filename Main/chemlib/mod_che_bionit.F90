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

  ! manure application rate (kg/m2/s)
  ! fertiliser application rate (kg/m2/s)
  ! soilpH/
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
    if ( ichem == 1 .and. ichbion == 1 ) then
      call getmem2d(nfert,jci1,jci2,ici1,ici2,'che_bionit:nfert')
      call getmem2d(nmanure,jci1,jci2,ici1,ici2,'che_bionit:nmanure')
      call getmem2d(soilph,jci1,jci2,ici1,ici2,'che_bionit:soilph')
    end if
  end subroutine allocate_mod_che_bionit

  subroutine ini_bionit
    implicit none
    if ( ichbion == 1 ) call read_bionem (nfert,nmanure,soilph)
  end subroutine ini_bionit

  subroutine soilnitro_emissions(ivegcov,wid10)
    implicit none
    real(rkx) , dimension(jci1:jci2,ici1:ici2) , intent(in) :: wid10
    integer(ik4) , dimension(jci1:jci2,ici1:ici2), intent(in) :: ivegcov
    ! local variables
    integer(ik4) :: i , j, nt

    real(rkx) :: &
      canred        ,&    ! canopy reduction factor
      fracnox       ,&  ! fraction emitted as NOx
      fracnh3       ,&  ! fraction emitted as NH3
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
    real(rkx), dimension(jci1:jci2,ici1:ici2) :: &
         soiltemp_surf ,&    ! surface soil temperature (C)
         soiltemp_deep ,&    ! deep soil temperature (C)
         sandper       ,&    ! sand percentage (%)
         porewater     ,&    ! pore space water content (%)
         windsp        ,&    ! wind speed (m/s)
         soilfert      ,&    ! soil fertilisation rate (kg/m2/hr)
         lai_int       ,&    ! manure/fertiliser app. rate + pH
         noxflux       ,&    ! calculated soil NOx flux
         totn                ! total N app. rate

    ! convert from kg/ha/year to kg/ha/hr needed by the neural network

    totn = (nmanure + nfert)/(24._rkx * 365._rkx)

    ! iFAB  ! put interactive LAI
    lai_int = cxlai2d(jci1:jci2,ici1:ici2)

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

#ifndef CLM45
    do i = ici1 , ici2
      do j = jci1 , jci2
        ! cycle on sea points
        if ( ivegcov(j,i) == 0 ) cycle
        ! getting the soil sand percentage, pH and fert rate values
        ! ABSOLUMENT àCHANGER
        !sandper(j,i) = sandrow2(j,i)
        sandper(j,i)  = 30.
        ! calculating water-filled pore space from soil moisture
        porewater(j,i) = cssw2da(j,i)
        ! converting soil moisture from kg/m2 to (m3 water/m3 soil)
        ! divide first by depth of soil layer = 10cm = 0.1m : cdepuv = 100mm
        ! then divide by density of water = 1000kg/m3
        porewater(j,i) = (porewater(j,i) / &
                   (cdepuv(ivegcov(j,i)) * d_1000)) * d_r1000
        ! calculating water-filled pore space (%)
        ! coefficient of 0.45 derived from obs at Grignon(0.536),
        ! Hombori(0.4) and Escompte(0.43) = avg. 0.45
        ! in regcm/bats  this parameter would be cxmopor : consider replacing ?
        porewater(j,i) = (porewater(j,i) / 0.45_rkx) * d_100
        ! converting temperature from Kelvin to Celsius
        soiltemp_deep(j,i) = ctg(j,i)
        soiltemp_surf(j,i) = ctga(j,i)
        if ( soiltemp_deep(j,i) /= 0 ) then
          soiltemp_deep(j,i) = soiltemp_deep(j,i) - tzero
        end if
        if ( soiltemp_surf(j,i) /= 0 ) then
          soiltemp_surf(j,i) = soiltemp_surf(j,i) - tzero
        end if
      end do
    end do
#endif
#ifdef CLM45
    do i = ici1 , ici2
      do j = jci1 , jci2
        ! cycle on sea points
        if ( ivegcov(j,i) == 0 ) cycle
        ! getting the soil sand percentage, pH and fert rate values
        do nt = 1, nats
          sandper(j,i) = sandper(j,i)  + dustsotex(j,i,nt)*fsand(nt)*100_rkx 
        end do 
        ! calculating water-filled pore space from soil moisture
        ! csw_vol voluletric soil moist (m3/m3)
        ! here consider second soil level / to be perhaps tested
        porewater(j,i) = csw_vol(j,i,2)

        ! calculating water-filled pore space (%)
        ! coefficient of 0.45 derived from obs at Grignon(0.536),
        ! Hombori(0.4) and Escompte(0.43) = avg. 0.45
        ! in regcm/bats  this parameter would be cxmopor : consider replacing ?
        porewater(j,i) = (porewater(j,i) / 0.45_rkx) * d_100

        ! temperature profile
        ! converting temperature from Kelvin to Celsius
        soiltemp_deep(j,i) = ctsoi(j,i,3)
        soiltemp_surf(j,i) = ctsoi(j,i,1)

        if ( soiltemp_deep(j,i) /= 0 ) then
          soiltemp_deep(j,i) = soiltemp_deep(j,i) - tzero
        end if
        if ( soiltemp_surf(j,i) /= 0 ) then
          soiltemp_surf(j,i) = soiltemp_surf(j,i) - tzero
        end if
      end do
    end do

#endif

    do i = ici1 , ici2
      do j = jci1 , jci2

        ! calculating what percentage volatilised N gets incorporated
        ! into NH3 and NOx

        fracnh3 = totn(j,i)*0.3_rkx
        fracnox = totn(j,i)*0.7_rkx

        ! calculation of NOx flux from soil
        ! normalised centered entries

        norm_ss = xcoef1  + xcoef2*soiltemp_surf(j,i)
        norm_sm = xcoef3  + xcoef4*porewater(j,i)
        norm_sd = xcoef5  + xcoef6*soiltemp_deep(j,i)
        norm_fe = xcoef7  + xcoef8*fracnox
        norm_sa = xcoef9  + xcoef10*sandper(j,i)
        norm_ph = xcoef11 + xcoef12*soilph(j,i)
        norm_wi = xcoef13 + xcoef14*wid10(j,i)

        ! weighted sums (coefficients from soil_nox_params)
        nsum1 = xwgt0 + xwgt1*norm_ss &
              + xwgt2*norm_sm + xwgt3*norm_sd &
              + xwgt4*norm_fe + xwgt5*norm_sa &
              + xwgt6*norm_ph + xwgt7*norm_wi

        nsum2 = xwgt8 + xwgt9*norm_ss &
              + xwgt10*norm_sm + xwgt11*norm_sd &
              + xwgt12*norm_fe + xwgt13*norm_sa &
              + xwgt14*norm_ph + xwgt15*norm_wi

        nsum3 = xwgt16 + xwgt17*norm_ss &
              + xwgt18*norm_sm + xwgt19*norm_sd  &
              + xwgt20*norm_fe + xwgt21*norm_sa &
              + xwgt22*norm_ph + xwgt23*norm_wi

        ! hyperbolic tangent calculation
        norm_no = xwgt24 + xwgt25*tanh(nsum1) &
              + xwgt26*tanh(nsum2) + xwgt27*tanh(nsum3)

        !flux calculation
        ! If sand > 50%, pulse effect, amplitude coefficient is maximum.
        ! If sand < 50%, amplitude coefficient is reduced
        !                to avoid strong emissions
        ! Sand conditions are correlated to pH values.
        if ( soilph(j,i) >= 6._rkx ) then
          noxflux(j,i) = xcoef15 + xcoef16s*norm_no
        else
          noxflux(j,i) = xcoef15 + xcoef16l*norm_no
        end if

        ! avoiding negative fluxes
        if ( noxflux(j,i) < d_zero ) then
          noxflux(j,i)= d_zero
        end if

        ! converting the NO flux from gN/ha/d to kg/m2/s
        ! g to kg: /1000
        ! ha to m2: /100 /100
        ! d to s: /86400
        noxflux(j,i) = noxflux(j,i) * &
          (30.0_rkx/14.0_rkx)/ (1000.0_rkx*100.0_rkx*100.0_rkx*86400.0_rkx)

        ! flux reduction because of canopy absorption
        if ( lai_int(j,i) > 1.9_rkx .and. lai_int(j,i) <  5._rkx ) then
          canred = 0.5_rkx
        else if ( lai_int(j,i) > 5._rkx ) then
          canred = 0.2_rkx
        else
          canred = d_one
        end if
        noxflux(j,i) = noxflux(j,i)*canred
      end do
    end do

    !if ( j == 25 ) write(stdout,*) maxval(noxflux)
    !
    ! update tendency for NO flux
    if ( ichdrdepo == 1 ) then
      if ( idynamic == 3 ) then
        do i  = ici1 , ici2
          do j  = jci1 , jci2
            if ( ivegcov(j,i) == 0 ) cycle
            chiten(j,i,kz,ino) = chiten(j,i,kz,ino) + &
                   noxflux(j,i)/(cdzq(j,i,kz)*crhob3d(j,i,kz))
          end do
        end do
      else
        do i  = ici1 , ici2
          do j  = jci1 , jci2
            if ( ivegcov(j,i) == 0 ) cycle
            chiten(j,i,kz,ino) = chiten(j,i,kz,ino) + &
                   noxflux(j,i)*cpsb(j,i)/(cdzq(j,i,kz)*crhob3d(j,i,kz))
          end do
        end do
      end if
    else if ( ichdrdepo == 2 ) then
      do i  = ici1 , ici2
        do j  = jci1 , jci2
          if ( ivegcov(j,i) == 0 ) cycle
          ! pass the flux to BL scheme
          chifxuw(j,i,ino) = chifxuw(j,i,ino) + noxflux(j,i)
        end do
      end do
    end if

    do i  = ici1 , ici2
      do j  = jci1 , jci2
        ! diagnostic source (accumulated)
        cemtrac(j,i,ino) = cemtrac(j,i,ino) + noxflux(j,i)* cfdout
      end do
    end do

    if ( ichdiag > 0 ) then
      do i  = ici1 , ici2
        do j  = jci1 , jci2
          cemisdiag(j,i,kz,ino) = cemisdiag(j,i,kz,ino) + &
               noxflux(j,i)/(cdzq(j,i,kz)*crhob3d(j,i,kz)) * cfdout
        end do
      end do
    end if
  end subroutine soilnitro_emissions

end module mod_che_bionit
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

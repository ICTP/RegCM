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

module mod_che_linox
  !
  ! Tracer convective transport
  !
  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_che_common
  use mod_dynparam

  implicit none

  private

  public :: linox_em

  contains

  subroutine linox_em(ivegcov)
    implicit none
    integer(ik4) , dimension(jci1:jci2,ici1:ici2), intent(in) :: ivegcov
    real(rkx), dimension(kz) :: amass , pic , pcg
    real(rkx), dimension(jci1:jci2,ici1:ici2,kz) :: znox_prod_ic
    real(rkx), dimension(jci1:jci2,ici1:ici2,kz) :: znox_prod_cg
    real(rkx):: zsum_lmass_ic , zsum_lmass_cg , flashrate , &
         pic_rate , pcg_rate , cloud_top_height , zthickice , zbeta , fac
    integer(ik4) :: i , j , k , kcumzero , kcumm10

    znox_prod_ic = d_zero
    znox_prod_cg = d_zero

    do i = ici1 , ici2
      do j = jci1 , jci2
        ! cloud top heigh from convection scheme
        if ( kcumtop(j,i) <= 0 ) cycle

        pic = d_zero
        pcg = d_zero

        cloud_top_height  =  cza(j,i,kcumtop(j,i)) + &
                   0.5_rkx * cdzq(j,i,kcumtop(j,i))

        ! Lightning frequency for continentaland maritime clouds as a function
        ! of cloud top height : Price and Rind

        if ( ivegcov(j,i) > 0 ) then
          ! Flashes/min over land boxes
          flashrate   = 3.44e-5_rkx * &
              ( ( cloud_top_height * 1.e-3_rkx )**4.9_rkx  )
        else
          flashrate  = 6.4e-4_rkx  * &
              ( ( cloud_top_height  * 1e-3_rkx )**1.73_rkx )
        end if
        ! flash per second parameterization Price and Rind
        flashrate = flashrate/60.0_rkx
        ! apply a scale factor ( Price and Rind, 1994)
        fac  = dx**2 * (mathpi/180.0_rkx)**2 / &
            (earthrad**2 *cos(cxlat(j,i) * mathpi/180.0_rkx))
        flashrate = flashrate*0.97241_rkx*exp(0.048203_rkx*fac)

        ! IC / CG flash ratio in function of icy cloud thickness
        kcumzero = 0
        icelev:  do  k = kcumbot(j,i), kcumtop(j,i),-1
          if ( ctb3d(j,i,k) < tzero ) then
            kcumzero = k
            exit icelev
          end if
        end do icelev

        ! if convective cloud is not icy, cycle to the next i column
        if ( kcumzero == 0 .or. kcumzero == kcumtop(j,i) ) cycle

        kcumm10 = 0
        minus10:  do  k = kcumzero, kcumtop(j,i),-1
          if (ctb3d(j,i,k) < tzero - 10.0_rkx ) then
            kcumm10 = k
            exit minus10
          end if
        end do minus10

        zthickice = (cloud_top_height - cza(j,i,kcumzero) ) * 1.e-3_rkx

        zbeta = (((0.021_rkx*zthickice - 0.648_rkx) * &
            zthickice   +  7.493_rkx) * &
            zthickice  - 36.540_rkx) * &
            zthickice + 63.090_rkx
        ! 5.5km < zthickice < 14km
        zbeta = min( 48.7_rkx, max( 0.19_rkx, zbeta ) )
        !

        amass(:) = 0.0_rkx
        zsum_lmass_ic = 0.0_rkx
        zsum_lmass_cg = 0.0_rkx

        ! IC_NOx production rate

        pic_rate = (zbeta / (1.0_rkx+zbeta))*flashrate
        do k = kcumzero, kcumtop(j,i), -1
          amass(k) = crhob3d(j,i,k)* cdzq(j,i,k) * dx**2
          zsum_lmass_ic = zsum_lmass_ic + amass(k)
          pic(k) = pic_rate * 6.7e25_rkx * amass(k)
        end do
        pic = pic / zsum_lmass_ic
        !
        ! cg_nox production rate
        ! ic can happen without cg
        if ( kcumm10 > 0 ) then
          pcg_rate  = (1.0_rkx / (1.0_rkx+zbeta))*flashrate
          do k = kz, kcumm10, -1
            amass(k) = crhob3d(j,i,k)* cdzq(j,i,k) * dx**2
            zsum_lmass_cg = zsum_lmass_cg + amass(k)
            pcg(k) = pcg_rate * 6.7e26_rkx * amass(k)
          end do
          pcg = pcg / zsum_lmass_cg
        end if
        !
        ! Update units (molecules NO /s => kg.kg-1.s-1)
        ! update emission tendencies and diag
        !
        ! ------------------------------------------------
        !!
        where ( amass > 0.0_rkx )
          pic = pic * 30.e-3_rkx / (navgdr*amass)
          pcg = pcg * 30.e-3_rkx / (navgdr*amass)
        end where

        znox_prod_ic(j,i,:) = pic(:)
        znox_prod_cg(j,i,:) = pcg(:)

      end do ! en loop on i
    end do ! en loop on i

    if ( idynamic == 3 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            chiten(j,i,k,ino) = chiten(j,i,k,ino) + &
                (znox_prod_ic(j,i,k) + znox_prod_cg(j,i,k))
          end do
        end do
      end do
    else
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            chiten(j,i,k,ino) = chiten(j,i,k,ino) + &
                (znox_prod_ic(j,i,k) + znox_prod_cg(j,i,k)) * cpsb(j,i)
          end do
        end do
      end do
    end if

    if ( ichdiag > 0 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            cemisdiag(j,i,k,ino) = cemisdiag(j,i,k,ino ) + &
                (znox_prod_ic(j,i,k) +  znox_prod_cg(j,i,k)) * cfdout
          end do
        end do
      end do
    end if
  end subroutine linox_em

end module mod_che_linox
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

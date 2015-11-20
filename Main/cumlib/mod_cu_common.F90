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

module mod_cu_common

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_runparams , only : clfrcv , icumcloud , icup
  use mod_humid , only : clwfromt
  use mod_regcm_types

  implicit none

  public

  real(rk8) :: cevapu ! Raindrop evap rate coef [[(kg m-2 s-1)-1/2]/s]

  integer(ik4) , pointer , dimension(:,:) :: cuscheme ! which scheme to use
  real(rk8) , pointer , dimension(:,:,:) :: rain_cc
  integer(ik4) :: total_precip_points

  real(rk8) , dimension(10) :: cld_profile
  real(rk8) , dimension(10) :: fixed_cld_profile
  real(rk8) , dimension(10) :: rnum

  real(rk8) , parameter :: maxcloud_dp =  60000.0D0 ! In Pa
  logical , parameter :: addnoise = .false.

  contains

  subroutine init_mod_cumulus
    implicit none
    integer(ik4) , dimension(:) , allocatable:: iseed
    integer :: k , nseed
    real(rk4) :: cputime
    if ( icumcloud == 2 ) then
      !
      ! Free hand draw of a generic ten layer cumulus cloud shape.
      !
      fixed_cld_profile(1)  = 0.130D0
      fixed_cld_profile(2)  = 0.125D0
      fixed_cld_profile(3)  = 0.120D0
      fixed_cld_profile(4)  = 0.080D0
      fixed_cld_profile(5)  = 0.080D0
      fixed_cld_profile(6)  = 0.080D0
      fixed_cld_profile(7)  = 0.085D0
      fixed_cld_profile(8)  = 0.085D0
      fixed_cld_profile(9)  = 0.105D0
      fixed_cld_profile(10) = 0.110D0
      if ( addnoise ) then
        call random_seed(size=nseed)
        call cpu_time(cputime)
        allocate(iseed(nseed))
        iseed = int(cputime) + 37*(/(k-1,k=1,nseed)/)
        call random_seed(put=iseed)
        deallocate(iseed)
      else
        cld_profile = fixed_cld_profile
      end if
    end if
  end subroutine init_mod_cumulus

  subroutine model_cumulus_cloud(m2c,c2m)
    implicit none
    type(mod_2_cum) , intent(in) :: m2c
    type(cum_2_mod) , intent(inout) :: c2m
    real(rk8) :: akclth , scalep , scalef
    integer(ik4):: i , j , k , ktop , kbot , kclth , ikh
    scalef = (d_one-clfrcv)
    if ( icumcloud <= 1 ) then
      iloop1: &
      do i = ici1 , ici2
        jloop1: &
        do j = jci1 , jci2
          if ( cuscheme(j,i) /= 1 .and. cuscheme(j,i) /= 3 ) cycle jloop1
          ! The regcm model is top to bottom
          ktop = c2m%kcumtop(j,i)
          kbot = c2m%kcumbot(j,i)
          kclth = kbot - ktop + 1
          if ( kclth < 2 ) cycle jloop1
          akclth = d_one/dble(kclth)
          do k = ktop , kbot
            c2m%cldfrc(j,i,k) = d_one - scalef**akclth
          end do
        end do jloop1
      end do iloop1
    else if ( icumcloud == 2 ) then
      if ( addnoise ) then
        ! Put 25% noise level. Update cld_profile each time.
        call random_number(rnum)
        cld_profile = (0.75D0+(rnum/2.0D0))*fixed_cld_profile
      end if
      iloop3: &
      do i = ici1 , ici2
        jloop3: &
        do j = jci1 , jci2
          if ( cuscheme(j,i) /= 1 .and. cuscheme(j,i) /= 3 ) cycle jloop3
          ktop = c2m%kcumtop(j,i)
          kbot = c2m%kcumbot(j,i)
          kclth = kbot - ktop + 1
          if ( kclth < 2 ) cycle jloop3
          scalep = min((m2c%pas(j,i,kbot)-m2c%pas(j,i,ktop)) / &
                  maxcloud_dp,d_one)
          do k = ktop , kbot
            ikh = max(1,min(10,int((dble(k-ktop+1)/dble(kclth))*d_10)))
            c2m%cldfrc(j,i,k) = cld_profile(ikh)*clfrcv*scalep
          end do
        end do jloop3
      end do iloop3
    end if
    if ( icumcloud == 0 ) then
      do i = ici1 , ici2
        do j = jci1 , jci2
          ktop = c2m%kcumtop(j,i)
          kbot = c2m%kcumbot(j,i)
          if ( ktop > 1 .and. kbot > 1 ) then
            do k = ktop , kbot
              if ( cuscheme(j,i) == 4 ) then
                c2m%cldlwc(j,i,k) = 0.5D-4
              else
                c2m%cldlwc(j,i,k) = 0.3D-3
              end if
            end do
          end if
        end do
      end do
    else
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            if ( c2m%cldfrc(j,i,k) > 0.001D0 ) then
              c2m%cldlwc(j,i,k) = clwfromt(m2c%tas(j,i,k))
            else
              c2m%cldlwc(j,i,k) = d_zero
            end if
          end do
        end do
      end do
    end if
  end subroutine model_cumulus_cloud

end module mod_cu_common
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

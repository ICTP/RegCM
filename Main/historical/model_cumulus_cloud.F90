!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
! Historical, not compiled: the convective cloud of the Kuo and
! Betts-Miller schemes (cucloud of mod_cu_interface with icup 1 or 3),
! once in mod_cu_common. With icumcloud <= 1 the cloud fraction in the
! layers between cloud base and top is 1-(1-clfrcv)**(1/depth in
! levels); with icumcloud = 2 the ten layer profile below, scaled with
! the cloud depth (maxcloud_dp = 60000 Pa) and, with addnoise, a 25 %
! noise. Module data of mod_cu_common used here:
!
!  real(rkx), dimension(10) :: cld_profile, fixed_cld_profile, rnum
!  real(rkx), parameter :: maxcloud_dp =  60000.0_rkx ! In Pa
!  logical, parameter :: addnoise = .false.
!
! Profile set in init_mod_cumulus:
!
!    if ( icumcloud == 2 ) then
!      !
!      ! Free hand draw of a generic ten layer cumulus cloud shape.
!      !
!      fixed_cld_profile(1)  = 0.130_rkx
!      fixed_cld_profile(2)  = 0.125_rkx
!      fixed_cld_profile(3)  = 0.120_rkx
!      fixed_cld_profile(4)  = 0.080_rkx
!      fixed_cld_profile(5)  = 0.080_rkx
!      fixed_cld_profile(6)  = 0.080_rkx
!      fixed_cld_profile(7)  = 0.085_rkx
!      fixed_cld_profile(8)  = 0.085_rkx
!      fixed_cld_profile(9)  = 0.105_rkx
!      fixed_cld_profile(10) = 0.110_rkx
!      if ( .not. addnoise ) then
!        cld_profile = fixed_cld_profile
!      end if
!    end if
!
  subroutine model_cumulus_cloud(m2c)
    implicit none
    type(mod_2_cum), intent(in) :: m2c
    real(rkx) :: akclth, scalep, scalef
    integer(ik4):: i, j, k, ktop, kbot, kclth, ikh
    scalef = (d_one-clfrcv)
    if ( icumcloud <= 1 ) then
      do i = ici1, ici2
        do j = jci1, jci2
          ! The regcm model is top to bottom
          ktop = cu_ktop(j,i)
          kbot = cu_kbot(j,i)
          kclth = kbot - ktop + 1
          if ( kclth < 2 ) cycle
          akclth = d_one/real(kclth,rkx)
          do k = ktop, kbot
            cu_cldfrc(j,i,k) = d_one - scalef**akclth
          end do
        end do
      end do
    else if ( icumcloud == 2 ) then
      if ( addnoise ) then
        ! Put 25% noise level. Update cld_profile each time.
        call random_number(rnum)
        cld_profile = (0.75_rkx+(rnum/2.0_rkx))*fixed_cld_profile
      end if
      do i = ici1, ici2
        do j = jci1, jci2
          ktop = cu_ktop(j,i)
          kbot = cu_kbot(j,i)
          kclth = kbot - ktop + 1
          if ( kclth < 2 ) cycle
          scalep = min((m2c%pas(j,i,kbot)-m2c%pas(j,i,ktop)) / &
                  maxcloud_dp,d_one)
          do k = ktop, kbot
            ikh = max(1,min(10,int((real(k-ktop+1,rkx)/real(kclth,rkx))*d_10)))
            cu_cldfrc(j,i,k) = cld_profile(ikh)*clfrcv*scalep
          end do
        end do
      end do
    end if
  end subroutine model_cumulus_cloud

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

module mod_rad_tracer

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam , only : kz

  implicit none

  private

  public :: trcmix

  contains

  !-----------------------------------------------------------------------
  !
  ! Specify zonal mean mass mixing ratios of CH4, N2O, CFC11 and
  ! CFC12
  !          Code: J.T.Kiehl November 21, 1994
  !
  !-----------------------------------------------------------------------
  !
  !------------------------------input------------------------------------
  !
  ! pmid   - model pressures
  !
  !------------------------------output-----------------------------------
  !
  ! n2o    - nitrous oxide mass mixing ratio
  ! ch4    - methane mass mixing ratio
  ! cfc11  - cfc11 mass mixing ratio
  ! cfc12  - cfc12 mass mixing ratio
  !
  !-----------------------------------------------------------------------
  !
  ! This is now inline in CCM, and is used here by RRTM
  !
  pure subroutine trcmix(n1,n2,dlat,ptrop,pmid,n2o0,ch40,cfc110,cfc120, &
                         n2o,ch4,cfc11,cfc12)
    implicit none
    !
    ! dlat   - latitude in degrees
    ! xn2o   - pressure scale height for n2o
    ! xch4   - pressure scale height for ch4
    ! xcfc11 - pressure scale height for cfc11
    ! xcfc12 - pressure scale height for cfc12
    ! ptrop  - pressure level of tropopause
    ! pratio - pressure divided by ptrop
    !
    integer(ik4) , intent(in) :: n1 , n2
    real(rkx) , dimension(n1:n2) , intent(in) :: dlat , ptrop
    real(rkx) , dimension(n1:n2,kz) , intent(in) :: pmid
    real(rkx) , dimension(n1:n2) , intent(in) :: n2o0 , ch40 , cfc110 , cfc120
    real(rkx) , dimension(n1:n2,kz) , intent(out) :: cfc11 , cfc12 , ch4 , n2o
#ifndef RCEMIP
    real(rkx) :: alat
#endif
    real(rkx) :: pratio , xcfc11 , xcfc12 , xch4 , xn2o
    integer(ik4) :: n , k

#ifdef RCEMIP
    xn2o = 0.3478_rkx
    xch4 = 0.2353_rkx
    xcfc11 = 0.7273_rkx
    xcfc12 = 0.4000_rkx
#endif
    do k = 1 , kz
#ifdef STDPAR
#ifdef RCEMIP
      do concurrent ( n = n1:n2 ) local(pratio)
#else
      do concurrent ( n = n1:n2 ) local(pratio,alat)
#endif
#else
      do n = n1 , n2
#endif
#ifndef RCEMIP
        alat = abs(dlat(n))
        if ( alat <= 45.0_rkx ) then
          xn2o = 0.3478_rkx + 0.00116_rkx*alat
          xch4 = 0.2353_rkx
          xcfc11 = 0.7273_rkx + 0.00606_rkx*alat
          xcfc12 = 0.4000_rkx + 0.00222_rkx*alat
        else
          xn2o = 0.4000_rkx + 0.013333_rkx*(alat-45.0_rkx)
          xch4 = 0.2353_rkx + 0.0225489_rkx*(alat-45.0_rkx)
          xcfc11 = 1.00_rkx + 0.013333_rkx*(alat-45.0_rkx)
          xcfc12 = 0.50_rkx + 0.024444_rkx*(alat-45.0_rkx)
        end if
#endif
        !  set stratospheric scale height factor for gases
        if ( pmid(n,k) >= ptrop(n) ) then
          ch4(n,k) = ch40(n)
          n2o(n,k) = n2o0(n)
          cfc11(n,k) = cfc110(n)
          cfc12(n,k) = cfc120(n)
        else
          pratio = pmid(n,k)/ptrop(n)
          ch4(n,k) = ch40(n)*(pratio**xch4)
          n2o(n,k) = n2o0(n)*(pratio**xn2o)
          cfc11(n,k) = cfc110(n)*(pratio**xcfc11)
          cfc12(n,k) = cfc120(n)*(pratio**xcfc12)
        end if
      end do
    end do
  end subroutine trcmix

end module mod_rad_tracer
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

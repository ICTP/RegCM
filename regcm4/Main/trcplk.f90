!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of RegCM model.
!
!    RegCM model is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    RegCM model is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
      subroutine trcplk(tint,tlayr,tplnke,emplnk,abplnk1,abplnk2)
!
      use mod_regcm_param
      use mod_crdcon
      implicit none
!
!----------------------------------------------------------------------
!   Calculate Planck factors for absorptivity and emissivity of
!   CH4, N2O, CFC11 and CFC12
!
!-----------------------------------------------------------------------
!
! Input arguments
!
! tint    - interface temperatures
! tlayr   - k-1 level temperatures
! tplnke  - Top Layer temperature
!
! output arguments
!
! emplnk  - emissivity Planck factor
! abplnk1 - non-nearest layer Plack factor
! abplnk2 - nearest layer factor
!
! Dummy arguments
!
      real(8) , dimension(14,ixm1,kxp1) :: abplnk1 , abplnk2
      real(8) , dimension(14,ixm1) :: emplnk
      real(8) , dimension(ixm1,kxp1) :: tint , tlayr
      real(8) , dimension(ixm1) :: tplnke
      intent (in) tint , tlayr , tplnke
      intent (out) abplnk1 , abplnk2 , emplnk
!
! Local variables
!
! wvl   - wavelength index
! f1    - Planck function factor
! f2    -       "
! f3    -       "
!
      real(8) , dimension(14) :: f1 , f2 , f3
      integer :: i , k , wvl
!
      data f1/5.85713D8 , 7.94950D8 , 1.47009D9 , 1.40031D9 ,           &
         & 1.34853D8 , 1.05158D9 , 3.35370D8 , 3.99601D8 , 5.35994D8 ,  &
         & 8.42955D8 , 4.63682D8 , 5.18944D8 , 8.83202D8 , 1.03279D9/
                                    !        "
      data f2/2.02493D11 , 3.04286D11 , 6.90698D11 , 6.47333D11 ,       &
         & 2.85744D10 , 4.41862D11 , 9.62780D10 , 1.21618D11 ,          &
         & 1.79905D11 , 3.29029D11 , 1.48294D11 , 1.72315D11 ,          &
         & 3.50140D11 , 4.31364D11/
      data f3/1383.D0 , 1531.D0 , 1879.D0 , 1849.D0 , 848.D0 , 1681.D0 ,&
         & 1148.D0 , 1217.D0 , 1343.D0 , 1561.D0 , 1279.D0 , 1328.D0 ,  &
         & 1586.D0 , 1671.D0/
!
!     Calculate emissivity Planck factor
!
      do wvl = 1 , 14
        do i = 1 , ixm1
          emplnk(wvl,i) = f1(wvl)                                       &
                        & /(tplnke(i)**4.0*(dexp(f3(wvl)/tplnke(i))-1.0)&
                        & )
        end do
      end do
!
!     Calculate absorptivity Planck factor for tint and tlayr temperatures
!
      do wvl = 1 , 14
        do k = 1 , kxp1
          do i = 1 , ixm1
!           non-nearlest layer function
            abplnk1(wvl,i,k) = (f2(wvl)*dexp(f3(wvl)/tint(i,k)))        &
                             & /(tint(i,k)                              &
                             & **5.0*(dexp(f3(wvl)/tint(i,k))-1.0)**2.0)
!           nearest layer function
            abplnk2(wvl,i,k) = (f2(wvl)*dexp(f3(wvl)/tlayr(i,k)))       &
                             & /(tlayr(i,k)                             &
                             & **5.0*(dexp(f3(wvl)/tlayr(i,k))-1.0)     &
                             & **2.0)
          end do
        end do
      end do

      end subroutine trcplk

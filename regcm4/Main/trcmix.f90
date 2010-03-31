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
 
      subroutine trcmix(pmid,clat,coslat,n2o,ch4,cfc11,cfc12)
!
      use mod_regcm_param
      use mod_tracer
      implicit none
!
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
! clat   - current latitude in radians
! coslat - cosine of latitude
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
! Dummy arguments
!
      real(8) , dimension(iym1,kz) :: cfc11 , cfc12 , ch4 , n2o ,    &
           & pmid
      real(8) , dimension(iym1) :: clat , coslat
      intent (in) clat , coslat , pmid
      intent (out) cfc11 , cfc12 , ch4 , n2o
!
! Local variables
!
! i      - longitude loop index
! k      - level index
! dlat   - latitude in degrees
! xn2o   - pressure scale height for n2o
! xch4   - pressure scale height for ch4
! xcfc11 - pressure scale height for cfc11
! xcfc12 - pressure scale height for cfc12
! ptrop  - pressure level of tropopause
! pratio - pressure divided by ptrop
!
      real(8) :: dlat , pratio , xcfc11 , xcfc12 , xch4 , xn2o
      integer :: i , k
      real(8) , dimension(iym1) :: ptrop
!
      do i = 1 , iym1
!       set stratospheric scale height factor for gases
        dlat = dabs(57.2958*clat(i))
        if ( dlat.le.45.0 ) then
          xn2o = 0.3478 + 0.00116*dlat
          xch4 = 0.2353
          xcfc11 = 0.7273 + 0.00606*dlat
          xcfc12 = 0.4000 + 0.00222*dlat
        else
          xn2o = 0.4000 + 0.013333*(dlat-45)
          xch4 = 0.2353 + 0.0225489*(dlat-45)
          xcfc11 = 1.00 + 0.013333*(dlat-45)
          xcfc12 = 0.50 + 0.024444*(dlat-45)
        end if
!
!       pressure of tropopause
        ptrop(i) = 250.0E2 - 150.0E2*coslat(i)**2.0
!
      end do
!
      do k = 1 , kz
        do i = 1 , iym1
          if ( pmid(i,k).ge.ptrop(i) ) then
            ch4(i,k) = ch40
            n2o(i,k) = n2o0
            cfc11(i,k) = cfc110
            cfc12(i,k) = cfc120
          else
            pratio = pmid(i,k)/ptrop(i)
            ch4(i,k) = ch40*(pratio)**xch4
            n2o(i,k) = n2o0*(pratio)**xn2o
            cfc11(i,k) = cfc110*(pratio)**xcfc11
            cfc12(i,k) = cfc120*(pratio)**xcfc12
          end if
        end do
      end do

      end subroutine trcmix

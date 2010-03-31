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
 
      subroutine dragdn
!
!     returns neutral drag coefficient for grid square
!
!     zlnd = soil roughness length
!     zoce = ocean roughness length
!     zsno = snow roughness length
!     vonkar = von karman constant
!
!     frav = fraction of grid point covered by vegetation
!     fras = fraction of grid point covered by snow
!     frab = fraction of grid point covered by bare soil
!     cdb = neutral drag coeff over bare soil, ocean, sea ice
!     cds = neutral drag coeff over snow
!     cdv = neutral drag coeff over vegetation
!     cdrn = neutral drag coeff for momentum avgd over grid point
!
      use mod_regcm_param
      use mod_bats , only : displa , lveg , z1d , z1 , z1log , cdrn ,   &
                   & sigf , sice1d , ldoc1d , veg1d , scvk , rough , wt
      use mod_constants , only : zlnd , zoce , zsno , vonkar
      implicit none
!
! Local variables
!
      real(8) :: asigf , cdb , cds , cdv , frab , fras , frav
      integer :: n , i
!
!     ******           sea ice classified same as desert
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( lveg(n,i).le.0 .and. sice1d(n,i).gt.0. ) lveg(n,i) = 8
        end do
      end do
 
      call depth
 
      do i = 2 , iym1
        do n = 1 , nnsg
 
          z1(n,i) = z1d(n,i)
          z1log(n,i) = dlog(z1(n,i))
 
#ifdef SEAICE
           if ( ldoc1d(n,i).gt.1.5 ) then
             sigf(n,i) = 0.0
             cdrn(n,i) = ( vonkar / dlog( z1(n,i)/zlnd ) )**2
           else if ( ldoc1d(n,i).gt.0.5 ) then
#else
           if ( ldoc1d(n,i).gt.0.5 ) then
#endif
 
!           ******           drag coeff over land
            frav = sigf(n,i)
            asigf = veg1d(n,i)
            fras = asigf*wt(n,i) + (1.-asigf)*scvk(n,i)
            frab = (1.-asigf)*(1.-scvk(n,i))
            cdb = (vonkar/dlog(z1(n,i)/zlnd))**2
            cds = (vonkar/dlog(z1(n,i)/zsno))**2
            cdv = (vonkar/dlog((z1(n,i)-displa(lveg(n,i)))/             &
                & rough(lveg(n,i))))**2
            cdrn(n,i) = frav*cdv + frab*cdb + fras*cds
 
          else
!           ******           drag coeff over ocean
            sigf(n,i) = 0.0
            cdrn(n,i) = (vonkar/dlog(z1(n,i)/zoce))**2
          end if
        end do
      end do
 
      end subroutine dragdn

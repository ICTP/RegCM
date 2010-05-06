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
 
      subroutine depth
!
!          wt = fraction of vegetation covered by snow
!        sigf = fraction of veg. cover, excluding snow-covered veg.
!        scvk = fraction of soil covered by snow
!
!     scrat = snow depth (m) =  .001 snow depth (mm) / snow density
!     rhosw = ratio of snow density to density of h2o
!
!     height of vegetation assumed to be 10 x vegetation roughness ht
!     densi is defined the same here as in subr. albedo
!
!     wt now scaled so that runs betw. 0 & 1 and wt=0.5 when depth
!     of snow equals height of vegetation
!
 
      use mod_dynparam
      use mod_bats , only : lveg , rhosw , densi , scrat , wt , rough , &
                    & veg1d , sigf , scvk , ldoc1d , sag1d , scv1d
      implicit none
!
! Local variables
!
      real(8) :: age
      integer :: n , i
! 
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 ) then
            age = (1.-1./(1.+sag1d(n,i)))
            rhosw(n,i) = .10*(1.+3.*age)
            densi(n,i) = .01/(1.+3.*age)
            scrat(n,i) = scv1d(n,i)*densi(n,i)
            wt(n,i) = 1.0
            if ( lveg(n,i).gt.0 ) then
              wt(n,i) = 0.1*scrat(n,i)/rough(lveg(n,i))
              wt(n,i) = wt(n,i)/(1.+wt(n,i))
            end if
            sigf(n,i) = (1.-wt(n,i))*veg1d(n,i)
            scvk(n,i) = scrat(n,i)/(0.1+scrat(n,i))
          end if
        end do
      end do
 
      end subroutine depth

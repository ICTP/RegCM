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
 
      use regcm_param
      use bats
      implicit none
!
! Local variables
!
      real(8) :: age
      integer :: n , np
! 
      do np = np1 , npts
        do n = 1 , nnsg
          if ( ldoc1d(n,np).gt.0.5 ) then
            age = (1.-1./(1.+sag1d(n,np)))
            rhosw(n,np) = .10*(1.+3.*age)
            densi(n,np) = .01/(1.+3.*age)
            scrat(n,np) = scv1d(n,np)*densi(n,np)
            wt(n,np) = 1.0
            if ( lveg(n,np).gt.0 ) then
              wt(n,np) = 0.1*scrat(n,np)/rough(lveg(n,np))
              wt(n,np) = wt(n,np)/(1.+wt(n,np))
            end if
            sigf(n,np) = (1.-wt(n,np))*veg1d(n,np)
            scvk(n,np) = scrat(n,np)/(0.1+scrat(n,np))
          end if
        end do
      end do
 
      end subroutine depth

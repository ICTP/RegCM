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
 
      subroutine cldfrac(j)
 
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     This subroutine computes the fractional cloud coverage and      c
!     liquid water content (in cloud value).  Both are use in         c
!     radiation.                                                      c
!                                                                     c
!     The fractional coverage of large scale clouds is a function of  c
!     relative humidity, using the relationship of sundqvist et       c
!     al., 1989.  The relative humidity at which clouds begin to      c
!     form is lower over land than ocean, due to the greater number   c
!     of cloud condensation nucleii.                                  c
!                                                                     c
!     The fracional coverage of convective clouds is passed in from   c
!     the convection scheme.                                          c
!                                                                     c
!     The large-scale and convective clouds are combined as follows:  c
!     1) If the convective cloud fraction > large scale fraction, the c
!     convective fraction and water content are used (this occurs     c
!     infrequently).                                                  c
!     2) Otherwise, the cloud fraction equals the large-scale         c
!     fraction AND the water content is a weighted average of both    c
!     types.                                                          c
!                                                                     c
!     Note: the incloud water content (g/m3) is passed to radiation   c
!                                                                     c
!     See Pal et al (200) for more info.                              c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
 
      use mod_regcm_param
      use mod_pmoist
      use mod_rad
      use mod_slice
      implicit none
!
! Dummy arguments
!
      integer :: j
      intent (in) j
!
! Local variables
!
      real(8) :: exlwc , rh0adj
      integer :: i , k
!--------------------------------------------------------------------
!     1.  Determine large-scale cloud fraction
!--------------------------------------------------------------------
      do k = 1 , kz         ! Adjusted relative humidity threshold
        do i = 2 , iym2
          if ( tb3d(i,k,j).gt.tc0 ) then
            rh0adj = rh0(i,j)
          else ! high cloud (less subgrid variability)
            rh0adj = rhmax - (rhmax-rh0(i,j))                           &
                   & /(1.0+0.15*(tc0-tb3d(i,k,j)))
          end if
          if ( rhb3d(i,k,j).ge.rhmax ) then    ! full cloud cover
            fcc(i,k,j) = 1.
          else if ( rhb3d(i,k,j).le.rh0adj ) then
                                               ! no cloud cover
            fcc(i,k,j) = 0.
          else                                ! partial cloud cover
            fcc(i,k,j) = 1. - sqrt(1.-(rhb3d(i,k,j)-rh0adj)/(rhmax-     &
                       & rh0adj))
            fcc(i,k,j) = dmin1(dmax1(fcc(i,k,j),0.01D0),0.99D0)
          end if                      !  rh0 threshold
        end do
      end do
 
!--------------------------------------------------------------------
!     2.  Combine large-scale and convective fraction and liquid water
!     to be passed into radiation.
!--------------------------------------------------------------------
      do k = 1 , kz
        do i = 2 , iym2
          exlwc = qcb3d(i,k,j)*rhob3d(i,k,j)*1.E03
                                                  ! Cloud Water Volume
!         units:  = kg gq / kg dry air * kg dry air / m3 * 1000 = g qc
!         / m3
          cldlwc(i,k) = (cldfra(i,k)*cldlwc(i,k)*1.E3+fcc(i,k,j)*exlwc) &
                      & /dmax1(cldfra(i,k)+fcc(i,k,j),0.01D0)
          cldfra(i,k) = dmin1(dmax1(cldfra(i,k),fcc(i,k,j)),fcmax)
        end do
      end do
 
      end subroutine cldfrac

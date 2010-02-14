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
 
      subroutine frawat

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  **** determines fraction of foliage covered by water fwet, and
!  **** the fraction of foliage that is dry transpiring leaf fdry.
!  note: their defns differ - fwet is the fraction of all veg surfaces
!  which are wet because mod_stems can evaporate, fdry is the fraction
!  of lai which is dry because only leaves can transpire
!
!  ldew1d(np) is in kg/m**2/s
!  fwet   = ratio of dew to max value to 2/3 power
!           ( 2/3 power comes from deardorff (1978) )
!              ** keep fwet le 1.0 **
!  dewmxi = inverse of max allowed dew depth on leaf in mm
!
      use mod_regcm_param
      use mod_bats , only : npts , ldoc1d , sigf , ldew1d , fwet ,      &
                   & xlsai , fdry , vegt , xlai
      use mod_constants , only : dewmxi
      implicit none
!
! Local variables
!
      integer :: n , np
!
      do np = np1 , npts
        do n = 1 , nnsg
          if ( ldoc1d(n,np).gt.0.5 ) then
            if ( sigf(n,np).gt.0.001 ) then
              fwet(n,np) = 0.
              if ( ldew1d(n,np).gt.0. ) then
                fwet(n,np) = ((dewmxi/vegt(n,np))*ldew1d(n,np))         &
                           & **.666666666666
                fwet(n,np) = dmin1(fwet(n,np),1.D0)
              end if
              fdry(n,np) = (1.-fwet(n,np))*xlai(n,np)/xlsai(n,np)
            end if
          end if
        end do
      end do
 
      end subroutine frawat

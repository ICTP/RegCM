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
 
      subroutine frawat

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!  **** determines fraction of foliage covered by water fwet, and
!  **** the fraction of foliage that is dry transpiring leaf fdry.
!  note: their defns differ - fwet is the fraction of all veg surfaces
!  which are wet because mod_stems can evaporate, fdry is the fraction
!  of lai which is dry because only leaves can transpire
!
!  ldew1d(i) is in kg/m**2/s
!  fwet   = ratio of dew to max value to 2/3 power
!           ( 2/3 power comes from deardorff (1978) )
!              ** keep fwet le 1.0 **
!  dewmxi = inverse of max allowed dew depth on leaf in mm
!
      use mod_regcm_param
      use mod_bats , only : xlai , ldoc1d , sigf , ldew1d , fwet ,      &
                   & xlsai , fdry , vegt
      use mod_constants , only : dewmxi
      implicit none
!
! Local variables
!
      integer :: n , i
!
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 ) then
            if ( sigf(n,i).gt.0.001 ) then
              fwet(n,i) = 0.
              if ( ldew1d(n,i).gt.0. ) then
                fwet(n,i) = ((dewmxi/vegt(n,i))*ldew1d(n,i))**(2.0/3.0)
                fwet(n,i) = dmin1(fwet(n,i),1.D0)
              end if
              fdry(n,i) = (1.-fwet(n,i))*xlai(n,i)/xlsai(n,i)
            end if
          end if
        end do
      end do
 
      end subroutine frawat

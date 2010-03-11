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

      subroutine solar1clm(xtime)

!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
!     This subroutine calls clm routines to compute the following     c
!     orbital variables based on julian date:                         c
!                                                                     c
!     eccen   : orbital eccentricity                                  c
!     obliq   : obliquity in degrees                                  c
!     obliqr  : obliquity in radians                                  c
!     mvelp   : moving vernal equinox long                            c
!     mvelpp  : moving vernal equinox long of perihelion plus pi(rad) c
!     lambm0  : mean long of perihelion at vernal equinox (rad)       c
!     delta   : solar declination angle in rad                        c
!     eccf    : Earth-sun distance factor ((1/r)**2)                  c
!                                                                     c
!     Input Variables:                                                c
!     xtime   : forecast time in minutes.                             c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

!
      use shr_orb_mod
      use clm_varsur ,  only : numdays
      use shr_kind_mod , only : r8 => shr_kind_r8
!
      use mod_clm , only : r2ceccen , r2cobliqr , r2clambm0 ,           &
             &             r2cmvelpp , r2ceccf
      use mod_message , only : aline
      use mod_date , only : declin , julday , gmt , nnnnnn , nstrt0 ,   &
             &              idate1
      use mod_constants , only : degrad , dayspy

      implicit none
!
! Dummy arguments
!
      real(8) :: xtime
      intent (in) xtime
!
! Local variables
!
!----------------------------------------------------------------------
!
!     local variables for orbital params
!
! obliq     : obliquity in degrees
! mvelp     : moving vernal equinox long
! delta     : solar decl in radians
! iyear_ad  : Year to calculate orbit for
! log_print : Flags print of status/error
! decdeg    : solar decl in degrees
!
!----------------------------------------------------------------------
!
      real(8) :: calday , decdeg
      real(r8) :: delta , mvelp , obliq
      integer :: iyear_ad
      logical :: log_print
!
      
      log_print = .false.

      iyear_ad = idate1/1000000
      numdays = dayspy
 
!     Get eccen,obliq,mvelp,obliqr,lambm0,mvelpp
      call shr_orb_params(iyear_ad,r2ceccen,obliq,mvelp,r2cobliqr,      &
                        & r2clambm0,r2cmvelpp,log_print)
 
      calday = dble(julday) + (nnnnnn-nstrt0)/4. + (xtime/60.+gmt)/24.
 
!     Get delta,eccf
      call shr_orb_decl(calday,r2ceccen,r2cmvelpp,r2clambm0,r2cobliqr,  &
                      & delta,r2ceccf)
 
!     convert delta to degrees
      declin = delta
      decdeg = declin/degrad
 
!     abt rcm below
      write (aline, 99001) decdeg
      call say
99001 format (11x,'*** solar declination angle = ',f6.2,' degrees.')
 
      end subroutine solar1clm

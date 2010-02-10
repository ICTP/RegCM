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
 
      subroutine aerout(jslc,tauxar_mix,tauasc_mix,gtota_mix,aeradfo,   &
                      & aeradfos)
 
      use mod_regcm_param
      use mod_param2
      use mod_trachem
      use mod_aerosol , only : nspi
      implicit none
!
! Dummy arguments
!
      integer :: jslc
      real(8) , dimension(ix - 1) :: aeradfo , aeradfos
      real(8) , dimension(ix - 1,0:kx,nspi) :: gtota_mix , tauasc_mix ,&
           & tauxar_mix
      intent (in) aeradfo , aeradfos , gtota_mix , jslc , tauasc_mix ,  &
                & tauxar_mix
!
! Local variables
!
      integer :: i , k , ntim
! 
      do k = 1 , kx
        do i = 2 , ix - 1
#ifdef MPP1
          aerext(i-1,k,jslc) = tauxar_mix(i,k,8)
          aerssa(i-1,k,jslc) = tauasc_mix(i,k,8)
          aerasp(i-1,k,jslc) = gtota_mix(i,k,8)
#else
          aerext(i-1,k,jslc-1) = tauxar_mix(i,k,8)
          aerssa(i-1,k,jslc-1) = tauasc_mix(i,k,8)
          aerasp(i-1,k,jslc-1) = gtota_mix(i,k,8)
#endif
        end do
      end do
!
!     CARE :Average the radiative forcing between chem output time
!     steps (in hour) according to radfrq (in min), aertarf is reset to
!     0 at each chem output (cf output.f)
!
      ntim = 60*chemfrq/radfrq
!
!     aersol radative forcing (care cgs to mks after radiation scheme !)
!
      do i = 2 , ix - 1
#ifdef MPP1
        aertarf(i-1,jslc) = aertarf(i-1,jslc) + aeradfo(i)*1.E-3/ntim
        aersrrf(i-1,jslc) = aersrrf(i-1,jslc) + aeradfos(i)*1.E-3/ntim
#else
        aertarf(i-1,jslc-1) = aertarf(i-1,jslc-1) + aeradfo(i)          &
                            & *1.E-3/ntim
        aersrrf(i-1,jslc-1) = aersrrf(i-1,jslc-1) + aeradfos(i)         &
                            & *1.E-3/ntim
#endif
      end do
 
      end subroutine aerout

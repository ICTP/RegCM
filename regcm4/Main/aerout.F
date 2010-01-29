C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
C
C    This file is part of RegCM model.
C
C    RegCM model is free software: you can redistribute it and/or modify
C    it under the terms of the GNU General Public License as published by
C    the Free Software Foundation, either version 3 of the License, or
C    (at your option) any later version.
C
C    RegCM model is distributed in the hope that it will be useful,
C    but WITHOUT ANY WARRANTY; without even the implied warranty of
C    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
C    GNU General Public License for more details.
C
C    You should have received a copy of the GNU General Public License
C    along with RegCM model.  If not, see <http://www.gnu.org/licenses/>.
C
C::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

      subroutine aerout(jslc, tauxar_mix, tauasc_mix, gtota_mix,
     &                  aeradfo, aeradfos) 

      implicit none

#     include "parameter.inc"

#     include "include/param2.cb"
#     include "include/trachem.cb"
#     include "include/aeroppt.cb"
#     include "include/parrad.cb"

      integer jslc
!     aerosol extinction optical depth
      real(kind=8)  tauxar_mix(plond,0:plev,nspi)
      real(kind=8)  tauasc_mix(plond,0:plev,nspi)
      real(kind=8)  gtota_mix(plond,0:plev,nspi)
      real(kind=8)  aeradfo(plond), aeradfos(plond) 

      integer ntim, i, k

C      Visible cas ns=8

      do k = 1,kx
        do i = 2,ix-1
#ifdef MPP1
          aerext(i-1,k,jslc) = tauxar_mix(i,k,8)        
          aerssa(i-1,k,jslc) = tauasc_mix(i,k,8)
          aerasp(i-1,k,jslc) =  gtota_mix(i,k,8) 
#else
          aerext(i-1,k,jslc-1) = tauxar_mix(i,k,8)        
          aerssa(i-1,k,jslc-1) = tauasc_mix(i,k,8)
          aerasp(i-1,k,jslc-1) =  gtota_mix(i,k,8) 
#endif
        end do
      end do
C
C CARE :Average the radiative forcing between chem output time steps (in
C hour) according to radfrq (in min), aertarf is reset to 0 at each chem
C output (cf output.f)
C
      ntim = 60 * chemfrq / radfrq
C
C aersol radative forcing (care cgs to mks after radiation scheme !)        
C
      do i = 2, ix-1
#ifdef MPP1
        aertarf(i-1,jslc) = aertarf(i-1,jslc)+ aeradfo(i) *1.e-3/ntim
        aersrrf(i-1,jslc) = aersrrf(i-1,jslc)+aeradfos(i)*1.e-3/ntim
#else
        aertarf(i-1,jslc-1) = aertarf(i-1,jslc-1)+aeradfo(i)*1.e-3/ntim
        aersrrf(i-1,jslc-1) = aersrrf(i-1,jslc-1)+aeradfos(i)*1.e-3/ntim
#endif
      end do   

      return
      end subroutine aerout

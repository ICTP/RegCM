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

      subroutine radoz2(o3vmr   ,pint    ,plol    ,plos    )
C-----------------------------------------------------------------------
C
C Computes the path length integrals to the model interfaces given the
C ozone volume mixing ratio
C
C---------------------------Code history--------------------------------
C
C Original version:     CCM1
C Standardized:         J. Rosinski, June 1992
C Reviewed:             J. Kiehl, B. Briegleb, August 1992
C Mixing ratio version: Bruce Biegleb, September 1992
C
C-----------------------------------------------------------------------
c
c $Id: radoz2.F,v 1.1.1.1 1995/02/09 23:27:03 ccm2 Exp $
c $Author: ccm2 $
c
c
c $Id: implicit.h,v 1.1.1.1 1995/02/09 23:26:52 ccm2 Exp $
c $Author: ccm2 $
c
      implicit none
C-----------------------------------------------------------------------
c
#     include "parameter.inc"
c
#     include "include/parrad.cb"
c
C------------------------------Commons----------------------------------
#     include "include/comozp.cb"
C------------------------------Input arguments--------------------------
C
      real(kind=8)  o3vmr(plond,plev)   ! ozone volume mixing ratio
      real(kind=8)  pint(plond,plevp)   ! Model interface pressures
C
C----------------------------Output arguments---------------------------
C
      real(kind=8)  plol(plond,plevp),  ! Ozone prs weighted path length (cm)
     $     plos(plond,plevp)   ! Ozone path length (cm)
C
C---------------------------Local workspace-----------------------------
C
      integer   i,             ! longitude index
     $          k              ! level index
C
C-----------------------------------------------------------------------
C
C Evaluate the ozone path length integrals to interfaces;
C factors of .1 and .01 to convert pressures from cgs to mks:
C
C Bug fix, 24 May 1996:  the 0.5 and 0.25 factors removed.
C
      do i=1,plon
         plos(i,1) = 0.1*cplos*o3vmr(i,1)*pint(i,1)
         plol(i,1) = 0.01*cplol*o3vmr(i,1)*pint(i,1)*pint(i,1)
      end do
      do k=2,plevp
         do i=1,plon
            plos(i,k) = plos(i,k-1) + 0.1*cplos*o3vmr(i,k-1)*
     $                  (pint(i,k) - pint(i,k-1))
            plol(i,k) = plol(i,k-1) + 0.01*cplol*o3vmr(i,k-1)*
     $                  (pint(i,k)*pint(i,k)-pint(i,k-1)*pint(i,k-1))
         end do
      end do
C
      return
      end

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

cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      subroutine cldems(clwp, fice, rei, emis)
C-----------------------------------------------------------------------
C
C Compute cloud emissivity using cloud liquid water path (g/m**2)
C
C---------------------------Code history--------------------------------
C
C Original version:  J. Kiehl
C Standardized:      J. Rosinski, June 1992
C Reviewed:          J. Hack, J. Kiehl, August 1992
C
C-----------------------------------------------------------------------
c
c $Id: cldems.F,v 1.2 1995/02/17 21:28:09 jhack Exp $
c $Author: jhack $
c
      implicit none
c
C------------------------------Parameters-------------------------------
c
#     include "parameter.inc"
c
#     include "include/parrad.cb"
C
      real(kind=8)  kabs                   ! longwave absorption coeff (m**2/g)
C------------------------------Arguments--------------------------------
C
C Input arguments
C
      real(kind=8)  clwp(plond,plev),      ! cloud liquid water path (g/m**2)
     $     rei(plond,plev),       ! ice effective drop size (microns)
     $     fice(plond,plev)       ! fractional ice content within cloud
C
C Output arguments
C
      real(kind=8)  emis(plond,plev)       ! cloud emissivity (fraction)
C
C---------------------------Local workspace-----------------------------
C
      integer i,k                 ! longitude, level indices
      real(kind=8)  kabsl,                 ! longwave absorption coeff (m**2/g)
     $     kabsi                  ! ice absorption coefficient
      parameter (kabsl = 0.090361)
C
C-----------------------------------------------------------------------
C
      do k=1,plev
         do i=1,plon
            kabsi = 0.005 + 1./rei(i,k)
            kabs = kabsl*(1.-fice(i,k))+kabsi*fice(i,k)
            emis(i,k) = 1. - dexp(-1.66*kabs*clwp(i,k))
         end do
      end do
C
      return
      end

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

      subroutine outsrf

c ******    write bats fields to unit iutbat

      implicit none
#     include "parameter.inc"
#     include "include/param1.cb"
#     include "include/param2.cb"
#     include "include/iunits.cb"
#     include "include/date.cb"

      integer i,j,n

      integer numbat
      parameter( numbat=21+6 )
      integer nrcbat,nrcsub
      COMMON /BATUNIT/ nrcbat,nrcsub
#ifdef MPP1
      real(kind=4)  fbat_io(mjx-2,ix-2,numbat)
      COMMON /BATOUTio/fbat_io
      real(kind=4)  fbat(jxp,ix-2,numbat)
#else
      real(kind=4)  fbat(jx-2,ix-2,numbat)
#endif
      COMMON /BATOUT/fbat

c ****** check if at desired output time for bats variables
      write(*,*) 'BATS variables written at ',idatex,xtime
      if (iotyp.eq.2) write(iutbat) idatex
      do n=1,numbat
        if (iotyp.eq.1) then
          nrcbat=nrcbat+1
#ifdef MPP1
       write(iutbat,rec=nrcbat)((fbat_io(j,i,n),j=1,mjx-2),i=1,ix-2)
#else
          write(iutbat,rec=nrcbat)((fbat(j,i,n),j=1,jx-2),i=1,ix-2)
#endif
        else if (iotyp.eq.2) then
#ifdef MPP1
          write(iutbat) ((fbat_io(j,i,n),j=1,mjx-2),i=1,ix-2)
#else
          write(iutbat) ((fbat(j,i,n),j=1,jx-2),i=1,ix-2)
#endif
        end if
      end do

      return
      end subroutine outsrf

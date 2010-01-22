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

      subroutine radtap
c
      implicit none
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#     include "include/param1.cb"
#     include "include/param2.cb"
#     include "include/parrad.cb"
#     include "include/bats.cb"
#     include "include/bats.cb2"
#     include "include/rad.cb"
#     include "include/date.cb"
#     include "include/iunits.cb"

      integer nrad2d,nrad3d
      parameter(nrad2d=21,nrad3d=5)
#ifdef MPP1
      include 'mpiregcm.h'
      real(kind=4)  frad2d_io(mjx-2,ix-2,nrad2d)
      real(kind=4)  frad3d_io(mjx-2,ix-2,kx,nrad3d)
      COMMON /OUTRADio/ frad2d_io,frad3d_io
#else
      real(kind=4)  frad2d(jx-2,ix-2,nrad2d),frad3d(jx-2,ix-2,kx,nrad3d)
      COMMON /OUTRAD/ frad2d,frad3d
#endif
      integer nrcrad
      COMMON /RADREC/ nrcrad
C
C Radiation resolution and I/O parameters
C
c
      integer       i    ! latitude index
      integer       j    ! longitude index
      integer       k    ! level index
      integer n
c
c      real(kind=8)  cgsmks     ! conversion factor for fluxes from cgs to mks
c
c      data cgsmks / 1.e-3 /
c
c compute total radiative heating flux for the surface,
c converting units from cgs to mks:
c
      print *, 'Writing rad fields at ktau = ', ktau,idatex
      if (iotyp.eq.1) then
        do n=1,nrad3d
          if (n.ne.1) then ! skip h20mmr
            do k=kx,1,-1
              nrcrad=nrcrad+1
#ifdef MPP1
              write(iutrad,rec=nrcrad)
     &                     ((frad3d_io(j,i,k,n),j=1,mjx-2),i=1,ix-2)
#else
              write(iutrad,rec=nrcrad)
     &                     ((frad3d(j,i,k,n),j=1,jx-2),i=1,ix-2)
#endif
           enddo
          end if
        enddo
        do n=1,nrad2d
          if (n.lt.10) then ! skip everything from alb (10 to 21)
            nrcrad=nrcrad+1
#ifdef MPP1
            write(iutrad,rec=nrcrad)
     &                   ((frad2d_io(j,i,n),j=1,mjx-2),i=1,ix-2)
#else
            write(iutrad,rec=nrcrad)
     &                   ((frad2d(j,i,n),j=1,jx-2),i=1,ix-2)
#endif
          end if
        enddo
      else if (iotyp.eq.2) then
        write(iutrad) idatex
        do n=1,nrad3d
          if (n.ne.1) then ! skip h20mmr
            do k=kx,1,-1
#ifdef MPP1
              write(iutrad) ((frad3d_io(j,i,k,n),j=1,mjx-2),i=1,ix-2)
#else
              write(iutrad) ((frad3d(j,i,k,n),j=1,jx-2),i=1,ix-2)
#endif
            end do
          end if
        end do
        do n=1,nrad2d
          if (n.lt.10) then ! skip everything from alb (10 to 21)
#ifdef MPP1
            write(iutrad) ((frad2d_io(j,i,n),j=1,mjx-2),i=1,ix-2)
#else
            write(iutrad) ((frad2d(j,i,n),j=1,jx-2),i=1,ix-2)
#endif
          end if
        end do
      end if

      return
      end

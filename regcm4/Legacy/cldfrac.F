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

      subroutine cldfrac(j)

ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
c                                                                     c
c     This subroutine computes the fractional cloud coverage and      c
c     liquid water content (in cloud value).  Both are use in         c
c     radiation.                                                      c
c                                                                     c
c     The fractional coverage of large scale clouds is a function of  c
c     relative humidity, using the relationship of sundqvist et       c
c     al., 1989.  The relative humidity at which clouds begin to      c
c     form is lower over land than ocean, due to the greater number   c
c     of cloud condensation nucleii.                                  c
c                                                                     c
c     The fracional coverage of convective clouds is passed in from   c
c     the convection scheme.                                          c
c                                                                     c
c     The large-scale and convective clouds are combined as follows:  c
c     1) If the convective cloud fraction > large scale fraction, the c
c     convective fraction and water content are used (this occurs     c
c     infrequently).                                                  c
c     2) Otherwise, the cloud fraction equals the large-scale         c
c     fraction AND the water content is a weighted average of both    c
c     types.                                                          c
c                                                                     c
c     Note: the incloud water content (g/m3) is passed to radiation   c
c                                                                     c
c     See Pal et al (200) for more info.                              c
c                                                                     c
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

      implicit none
      integer j

#     include "parameter.inc"
#     include "include/pmoist.cb"
#     include "include/rad.cb"
#     include "include/slice.cb"

      integer i,k
      real(kind=8)  exlwc
      real(kind=8)  rh0adj  ! Adjusted relative humidity threshold


c--------------------------------------------------------------------
c    1.  Determine large-scale cloud fraction
c--------------------------------------------------------------------
      do k=1,kx
      do i=2,ilxm
        if (tb3d(i,k,j).gt.tc0) then
          rh0adj = rh0(i,j)
        else   ! high cloud (less subgrid variability)
          rh0adj = rhmax-(rhmax-rh0(i,j))/(1.0+0.15*(tc0-tb3d(i,k,j)))
        end if
        if (rhb3d(i,k,j).ge.rhmax) then        ! full cloud cover
          fcc(i,k,j) = 1.
        else if (rhb3d(i,k,j).le.rh0adj) then  ! no cloud cover
          fcc(i,k,j) = 0.
        else                                  ! partial cloud cover
          fcc(i,k,j) = 1.-sqrt(1.-(rhb3d(i,k,j)-rh0adj)/(rhmax-rh0adj))
          fcc(i,k,j) = dmin1(dmax1(fcc(i,k,j),0.01D0),0.99D0)
        end if                        !  rh0 threshold
      end do
      end do

c--------------------------------------------------------------------
c    2.  Combine large-scale and convective fraction and liquid water
c        to be passed into radiation.
c--------------------------------------------------------------------
      do k=1,kx
      do i=2,ilxm
        exlwc = qcb3d(i,k,j)*rhob3d(i,k,j)*1.e03  ! Cloud Water Volume
c units:  = kg gq / kg dry air * kg dry air / m3 * 1000 = g qc / m3
        cldlwc(i,k) = (cldfra(i,k)*cldlwc(i,k)*1.e3 + fcc(i,k,j)*exlwc)
     &              / dmax1(cldfra(i,k)+fcc(i,k,j),0.01d0)
        cldfra(i,k) = dmin1(dmax1(cldfra(i,k),fcc(i,k,j)),fcmax)
      end do
      end do

      return
      end subroutine cldfrac

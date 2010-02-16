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

      subroutine snow

C
C     update snow cover and snow age
C
C     three-part if block:
C       if snow cover < 0, then snow cover and snow age = 0
C       if antarctica, snow age = 0 (katabatic winds keep snow fresh)
C       if elsewhere, snow age follows given formulae
C
C        ps = snow precipitation rate
C     evaps = moisture flux from ground to atmosphere
C        sm = snow melt rate
C     sdrop = snow fallen from vegetation
C
C     aging of snow consists of three factors:
C           age1: snow crystal growth
C           age2: surface melting
C           age3: accumulation  of other particles, soot, etc., which
C                      is small in southern hemisphere
C
      implicit none
#     include "parameter.inc"
#     include "include/bats.cb"

      real(kind=8)  sold(NNSG,nbmax)
      real(kind=8)  age3,arg,age1,arg2,age2,tage,dela0,dela,dels,sge

      integer n, np

      age3 = 0.3

C=======================================================================
C        1.   partition soil evaporation and precipitation
C                    between water and snow
C=======================================================================
      do np = np1, npts
        do n = 1, NNSG
          if (ldoc1d(n,np).gt.0.5) then

            evapw(n,np) = fevpg(n,np)
            evaps(n,np) = scvk(n,np) * evapw(n,np)
            if (ldoc1d(n,np).gt.1.5) evaps(n,np) = fevpg(n,np)
            evapw(n,np) = (1.-scvk(n,np)) * evapw(n,np)
C
C ******                tm  is temperature of precipitation
            if (tm(n,np) .ge. c(67)) then
              pw(n,np) = prcp1d(n,np) * (1.-sigf(n,np))
              ps(n,np) = 0.0
            else
C ******                snowing
              pw(n,np) = 0.0
              ps(n,np) = prcp1d(n,np) * (1.-sigf(n,np))
            end if
          end if
        end do
      end do
C
C=======================================================================
C        2.   update snow cover
C=======================================================================
      do np = np1, npts
        do n = 1, NNSG
          if (ldoc1d(n,np).gt.0.5) then
            sold(n,np) = scv1d(n,np)
            scv1d(n,np) = scv1d(n,np)
     $            + c(4)*(ps(n,np)-evaps(n,np)-sm(n,np)) + sdrop(n,np)
            scv1d(n,np) = dmax1( scv1d(n,np), 0.d0 )
            sag1d(n,np) = dmax1( sag1d(n,np), 0.d0 )

C ******           snow cover except for antarctica
C=======================================================================
C        3.   increment non-dimensional "age" of snow;
C             10 mm snow restores surface to that of new snow.
C=======================================================================
            if (scv1d(n,np).gt.0.) then
              arg=5.e3*(1./c(67)-1./tg1d(n,np))
              age1=dexp(arg)
              arg2=dmin1(0.d0,10.*arg)
              age2=dexp(arg2)
              tage=age1+age2+age3
              dela0=1.e-6*c(4)
              dela=dela0*tage
              dels=0.1*dmax1( 0.d0, scv1d(n,np)-sold(n,np) )
              sge=(sag1d(n,np)+dela)*(1.0-dels)
              sag1d(n,np) = dmax1( 0.d0, sge )
            end if

C ******           antarctica
            if (scv1d(n,np).gt.800.) sag1d(n,np)=0.
          end if
        end do
      end do

      return
      end subroutine snow

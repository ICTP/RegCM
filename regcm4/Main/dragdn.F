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

      subroutine dragdn
C
C     returns neutral drag coefficient for grid square
C
C     zlnd = soil roughness length
C     zoce = ocean roughness length
C     zsno = snow roughness length
C     vonkar = von karman constant
C
C     frav = fraction of grid point covered by vegetation
C     fras = fraction of grid point covered by snow
C     frab = fraction of grid point covered by bare soil
C     cdb = neutral drag coeff over bare soil, ocean, sea ice
C     cds = neutral drag coeff over snow
C     cdv = neutral drag coeff over vegetation
C     cdrn = neutral drag coeff for momentum avgd over grid point
C

      implicit none
#ifdef MPP1
      include 'regcm.param2'
#else
      include 'regcm.param'
#endif
      include 'parame'
#     include "include/bats.cb"

      integer n,np
      real(kind=8)  frav,asigf,fras,frab,cdb,cds,cdv

C ******           sea ice classified same as desert
      do np = np1, npts
        do n = 1, NNSG
          if (lveg(n,np).le.0 .and. sice1d(n,np).gt.0.) lveg(n,np)=8
        end do
      end do

CAS put call to depth back  -- need this to be called with every BATS call
      call depth
CAS

      do np = np1, npts
        do n = 1, NNSG

C* moved from inside of if block (2/4/92)
           z1(n,np) = z1d(n,np)
           z1log(n,np) = dlog(z1(n,np))

Ccc        if (lveg(n,np).gt.0.and.(lveg(n,np).ne.14.and.
Ccc  $         lveg(n,np).ne.15)) then
           if ( ldoc1d(n,np).gt.0.5 )then

C ******           drag coeff over land
             frav = sigf(n,np)
C            sfac = 1. - fseas(tgb1d(n,np))
C            asigf  = veg1d(n,np) - seasf(lveg(n,np)) * sfac
             asigf  = veg1d(n,np)
             fras = asigf*wt(n,np) + (1.-asigf)*scvk(n,np)
             frab = (1.-asigf) * (1.-scvk(n,np))
             cdb  = ( vonkar / dlog( z1(n,np)/zlnd ) )**2
             cds  = ( vonkar / dlog( z1(n,np)/zsno ) )**2
             cdv  = ( vonkar / dlog((z1(n,np)-displa(lveg(n,np)))
     a        / rough(lveg(n,np))) )**2
             cdrn(n,np) = frav*cdv + frab*cdb + fras*cds

          else
C ******           drag coeff over ocean
            sigf(n,np) = 0.0
            cdrn(n,np) = ( vonkar / dlog(z1(n,np)/zoce) )**2
          end if
        end do
      end do

      return
      end subroutine dragdn

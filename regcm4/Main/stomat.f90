!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    ICTP RegCM is free software: you can redistribute it and/or modify
!    it under the terms of the GNU General Public License as published by
!    the Free Software Foundation, either version 3 of the License, or
!    (at your option) any later version.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
!    GNU General Public License for more details.
!
!    You should have received a copy of the GNU General Public License
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
 
      subroutine stomat

!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!     gives leaf stomatal resistance from environmental parameters
!             under conditions of no moisture stress
!
!     standard lai from xla=max & xlai0=min lai (set in bconst)
!     seasb = fseas(tgb1d(n,i) (set in bndry) is a seasonal
!             factor for reduced winter lai and root water uptake
!        fc = light sensitivity for crops and grasses and has inverse
!             radiation units (m**2/watt)
!      rlai = sum of leaf and stem area indices
!     rmax0 = 5000. s/m (maximum resistance)
!     radu & radl = visible light intensity in upper & lower canopy
!     ft & fb = the fractional intercepted photo-active radiation
!             per unit (leaf & stem) area in the top (upper) and
!             bottom (lower) canopies, respectively
!     radfi = average of upper and lower canopy light factors
!        rs = stomatal resistance = min.res. * rad.factor * leaf factor
!      trup = transmission of the upper canopy, assumed to be the same
!             for the lower canopy,i.e., trup=dexp(-0.5*g*rlai/czen),
!             where g = attenuation factor
!
!     documented in ncar tech note, dickinson et al., 1986
!     improved stomatal shading, dickinson, nov 88.
!
      use mod_regcm_param
      use mod_bats
      use mod_ictp01
      use mod_constants , only : rmax0
      implicit none
!
! Local variables
!
      real(8) :: difzen , g , radfi , seas , vpdf , x
      real(8) :: fseas
      real(8) , dimension(nnsg,iym1) :: fsol0 , fsold , radf , rmini ,  &
           & trup , trupd
      integer :: il , ilmax , n , i
      real(8) , dimension(10) :: rad , radd
!
!     ***** seasonal temperature factor
!
      fseas(x) = dmax1(0.0D0,1.-0.0016*dmax1(298.-x,0.D0)**2)
 
!     ***** g is average leaf crosssection per unit lai
!     ***** difzen is ave of inverse of cos of angle of diffuse vis
!     light ***** ilmax is number of canopy layers
!     ***** czen is cosine solar zenith angle for incident light
!     *****   (to spec from input data need a good treatment of diffuse
!     rad) ***** trup is transmission of direct beam light in one
!     canopy layer ***** trupd is transmission of diffuse light in one
!     canopy layer
      g = 0.5
      difzen = 2.0
      ilmax = 4
!*    delete fracd here to put in diffuse mod_radiation from ccm
!cc   fracd = difrat         !  from shuttleworth mods #2
 
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 ) then
            if ( sigf(n,i).gt.0.001 ) then
!             **********            zenith angle set in zenitm
              if ( czen(i).gt.0.001 ) then
                trup(n,i) = dexp(-g*rlai(n,i)/(ilmax*czen(i)))
                trupd(n,i) = dexp(-difzen*g*rlai(n,i)/(ilmax))
                fsold(n,i) = fracd(i)*solis(i)*fc(lveg(n,i))
                fsol0(n,i) = (1.-fracd(i))*solis(i)*fc(lveg(n,i))
                rmini(n,i) = rsmin(lveg(n,i))/rmax0
              end if
            end if
          end if
        end do
      end do
 
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 ) then
            if ( sigf(n,i).gt.0.001 ) then
              if ( czen(i).gt.0.001 ) then
                rad(1) = (1.-trup(n,i))*fsol0(n,i)*ilmax/rlai(n,i)
                radd(1) = (1.-trupd(n,i))*fsold(n,i)*ilmax/rlai(n,i)
                do il = 2 , ilmax
                  rad(il) = trup(n,i)*rad(il-1)
                  radd(il) = trupd(n,i)*radd(il-1)
                end do
                radfi = 0.
                do il = 1 , ilmax
                  radfi = radfi + (rad(il)+radd(il)+rmini(n,i))         &
                        & /(1.+rad(il)+radd(il))
                end do
                radf(n,i) = ilmax/radfi
              end if
            end if
          end if
        end do
      end do
 
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i).gt.0.5 ) then
            if ( sigf(n,i).gt.0.001 ) then
              if ( czen(i).gt.0.001 ) then
                vpdf = 1./dmax1(0.3D0,1.D0-vpdc(n,i)*0.025)
                seas = 1./(rmini(n,i)+fseas(tlef1d(n,i)))
                rs(n,i) = rsmin(lveg(n,i))*radf(n,i)*seas*vpdf
                rs(n,i) = dmin1(rs(n,i),rmax0)
              else
                rs(n,i) = rmax0
              end if
            end if
          end if
        end do
      end do
      end subroutine stomat

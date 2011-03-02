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
 
      module mod_drag
!
      use mod_constants
      use mod_runparams
      use mod_dynparam
      use mod_bats
!
      private
!
      public :: drag , depth
!
      contains
! 
!=======================================================================
!l  based on: bats version 1e          copyright 18 august 1989
!=======================================================================
!
!     *** determines surface transfer coeffs. at anemom. level from
!     *** lowest model level based on monin-obukov theory using
!     *** deardorff parameterization in terms of bulk richardson no.
! 
!     ****  a.  calculates neutral drag coefficient (cdrn) as a fn of
!     ****             underlying surface
! 
!     ****  b.  modifies cdrn as fn of bulk rich. no. of surface layer
!
!=======================================================================
!
      subroutine drag
! 
      implicit none
!
      real(8) , dimension(nnsg,iym1) :: cdrmin , rib , ribl , ribn
      real(8) :: dthdz , u1 , u2 , zatild
      integer :: n , i
!
!=======================================================================
!     1.   get neutral drag coefficient
!=======================================================================
      call dragdn
 
      do i = 2 , iym1
        do n = 1 , nnsg
!=======================================================================
!         2.   compute stability as bulk rich. no. = rin/rid =
!         ri(numerator)/ri(denominator)
!=======================================================================
          if ( lveg(n,i) /= 0 ) then
            zatild = (z1(n,i)-displa(lveg(n,i)))*sigf(n,i) + z1(n,i)    &
                   & *(1.0D0-sigf(n,i))
          else
            zatild = z1(n,i)
          end if
          ribn(n,i) = zatild*gti*(ts1d(n,i)-sigf(n,i)*taf1d(n,i)-       &
                     & (1.0D0-sigf(n,i))*tg1d(n,i))/ts1d(n,i)
!=======================================================================
!         2.1  compute the bulk richardson number;
!         first get avg winds to use for ri number by summing the
!         squares of horiz., vertical, and convective velocities
!=======================================================================
          if ( ribn(n,i) <= 0.0D0 ) then
            dthdz = (1.0D0-sigf(n,i))*tg1d(n,i) + sigf(n,i)*taf1d(n,i) &
                  & - ts1d(n,i)
            u1 = wtur + 2.0D0*dsqrt(dthdz)
            ribd(n,i) = us1d(i)**2.0D0 + vs1d(i)**2.0D0 + u1**2.0D0
          else
            u2 = wtur
            ribd(n,i) = us1d(i)**2.0D0 + vs1d(i)**2.0D0 + u2**2.0D0
          end if
          vspda(n,i) = dsqrt(ribd(n,i))
          if ( vspda(n,i) < 1.0D0 ) then
            vspda(n,i) = 1.0D0
            ribd(n,i) = 1.0D0
          end if
          rib(n,i) = ribn(n,i)/ribd(n,i)
!=======================================================================
!         3.   obtain drag coefficient as product of neutral value
!         and stability correction
!=======================================================================
!         ****   -0.4 < rib < 0.2   (deardorff, jgr, 1968, 2549-2557)
          if ( rib(n,i) < 0.0D0 ) then
            cdr(n,i) = cdrn(n,i)* &
                     (1.0D0+24.5D0*dsqrt(-cdrn(n,i)*rib(n,i)))
          else
            cdr(n,i) = cdrn(n,i)/(1.0D0+11.5D0*rib(n,i))
          end if
!         3.1  apply lower limit to drag coefficient value
          cdrmin(n,i) = dmax1(0.25D0*cdrn(n,i),6.D-4)
          if ( cdr(n,i) < cdrmin(n,i) ) cdr(n,i) = cdrmin(n,i)
          cdrx(n,i) = cdr(n,i)
 
        end do
      end do
 
!=======================================================================
!     4.   obtain drag coefficient over sea ice as weighted average
!     over ice and leads
!     warning! the lat test below (4.1-4.3) is model dependent!
!=======================================================================
 
!     4.1  test if northern or southern hemisphere
      do i = 2 , iym1
        do n = 1 , nnsg
                                                    ! check each point
!cc   if(lat(i) ==     1) aarea(i) = 0.005  ! ccm specific code
!cc   if(lat(i) ==     2) aarea(i) = 0.01
!cc   if(lat(i) >= nlat2) aarea(i) = 0.04   !  4.2  antarctic
          if ( ldoc1d(n,i) > 1.5D0 ) aarea(n,i) = 0.02D0
                                                    !  4.3  arctic
        end do
      end do
 
!     4.4  neutral cd over lead water
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i) > 1.5D0 ) then       !  check each point
            cdrn(n,i) = (vonkar/dlog(z1(n,i)/zoce))**2.0D0
 
!           4.5  drag coefficient over leads
            ribl(n,i) = (1.0D0-271.5D0/ts1d(n,i))*z1(n,i)*gti/ribd(n,i)
            if ( ribl(n,i) >= 0.0D0 ) then
              clead(n,i) = cdrn(n,i)/(1.0D0+11.5D0*ribl(n,i))
            else
              clead(n,i) = cdrn(n,i)*(1.0D0+24.5D0*dsqrt(-cdrn(n,i)* &
                    & ribl(n,i)))
            end if
 
!           4.6  calculate weighted avg of ice and lead drag
!           coefficients
            cdrx(n,i) = (1.0D0-aarea(n,i))*cdr(n,i) + &
                               aarea(n,i)*clead(n,i)
          end if
        end do
      end do
 
      end subroutine drag
!
!=======================================================================
! DRAGDN
!
!     returns neutral drag coefficient for grid square
!
!     zlnd = soil roughness length
!     zoce = ocean roughness length
!     zsno = snow roughness length
!     vonkar = von karman constant
!
!     frav = fraction of grid point covered by vegetation
!     fras = fraction of grid point covered by snow
!     frab = fraction of grid point covered by bare soil
!     cdb = neutral drag coeff over bare soil, ocean, sea ice
!     cds = neutral drag coeff over snow
!     cdv = neutral drag coeff over vegetation
!     cdrn = neutral drag coeff for momentum avgd over grid point
!
!=======================================================================
!
      subroutine dragdn
!
      implicit none
!
      real(8) :: asigf , cdb , cds , cdv , frab , fras , frav
      integer :: n , i
!
!     ******           sea ice classified same as desert
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( lveg(n,i) <= 0 .and. sice1d(n,i) > 0.0D0 ) lveg(n,i) = 8
        end do
      end do
 
      call depth
 
      do i = 2 , iym1
        do n = 1 , nnsg
 
          z1(n,i) = z1d(n,i)
          z1log(n,i) = dlog(z1(n,i))
 
           if ( ldoc1d(n,i) > 1.5D0 ) then
             sigf(n,i) = 0.0D0
             cdrn(n,i) = ( vonkar / dlog( z1(n,i)/zlnd ) )**2.0D0
           else if ( ldoc1d(n,i) > 0.5D0 ) then
!           ******           drag coeff over land
            frav = sigf(n,i)
            asigf = veg1d(n,i)
            fras = asigf*wt(n,i) + (1.0D0-asigf)*scvk(n,i)
            frab = (1.0D0-asigf)*(1.0D0-scvk(n,i))
            cdb = (vonkar/dlog(z1(n,i)/zlnd))**2.0D0
            cds = (vonkar/dlog(z1(n,i)/zsno))**2.0D0
            cdv = (vonkar/dlog((z1(n,i)-displa(lveg(n,i)))/             &
                & rough(lveg(n,i))))**2.0D0
            cdrn(n,i) = frav*cdv + frab*cdb + fras*cds
 
          else
!           ******           drag coeff over ocean
            sigf(n,i) = 0.0D0
            cdrn(n,i) = (vonkar/dlog(z1(n,i)/zoce))**2.0D0
          end if
        end do
      end do
 
      end subroutine dragdn
!
!=======================================================================
! SNOW DEPTH
!
!          wt = fraction of vegetation covered by snow
!        sigf = fraction of veg. cover, excluding snow-covered veg.
!        scvk = fraction of soil covered by snow
!
!     scrat = snow depth (m) =  .001 snow depth (mm) / snow density
!     rhosw = ratio of snow density to density of h2o
!
!     height of vegetation assumed to be 10 x vegetation roughness ht
!     densi is defined the same here as in subr. albedo
!
!     wt now scaled so that runs betw. 0 & 1 and wt=0.5 when depth
!     of snow equals height of vegetation
!
!=======================================================================
!
      subroutine depth
!
      implicit none
!
      real(8) :: age
      integer :: n , i
! 
      do i = 2 , iym1
        do n = 1 , nnsg
          if ( ldoc1d(n,i) > 0.5D0 ) then
            age = (1.0D0-1.0D0/(1.0D0+sag1d(n,i)))
            rhosw(n,i) = 0.10D0*(1.0D0+3.0D0*age)
            densi(n,i) = 0.01D0/(1.0D0+3.0D0*age)
            scrat(n,i) = scv1d(n,i)*densi(n,i)
            wt(n,i) = 1.0D0
            if ( lveg(n,i) > 0 ) then
              wt(n,i) = 0.1D0*scrat(n,i)/rough(lveg(n,i))
              wt(n,i) = wt(n,i)/(1.0D0+wt(n,i))
            end if
            sigf(n,i) = (1.0D0-wt(n,i))*veg1d(n,i)
            scvk(n,i) = scrat(n,i)/(0.1D0+scrat(n,i))
          end if
        end do
      end do
! 
      end subroutine depth
!
      end module mod_drag

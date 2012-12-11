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
 
module mod_bats_drag
!
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_bats_common
  use mod_bats_internal
!
  private
!
  public :: dragc , depth
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
  subroutine dragc
! 
  implicit none
!
  real(rk8) :: dthdz , u1 , u2 , zatild , cdrmin
  integer(ik4) :: n , i , j
  !
  !=======================================================================
  !     1.   get neutral drag coefficient
  !=======================================================================
  !
  call dragdn
 
  do i = ici1 , ici2
    do j = jci1 , jci2
      do n = 1 , nnsg
        !=======================================================================
        ! 2.  compute stability as bulk rich. no. = rin/rid =
        !     ri(numerator)/ri(denominator)
        !=======================================================================
        if ( ldmsk1(n,j,i) /= 0 ) then
          zatild = (zh(n,j,i)-displa(lveg(n,j,i)))*sigf(n,j,i) + &
                    zh(n,j,i)*(d_one-sigf(n,j,i))
          ribn(n,j,i) = zatild*egrav*(sts(n,j,i)-sigf(n,j,i)*taf(n,j,i)- &
                    (d_one-sigf(n,j,i))*tgrd(n,j,i))/sts(n,j,i)
        else
          ribn(n,j,i) = zh(n,j,i)*egrav*(sts(n,j,i)-tgrd(n,j,i))/sts(n,j,i)
        end if
        !=======================================================================
        ! 2.1  compute the bulk richardson number;
        !      first get avg winds to use for ri number by summing the
        !      squares of horiz., vertical, and convective velocities
        !=======================================================================
        if ( ribn(n,j,i) <= d_zero ) then
          dthdz = (d_one-sigf(n,j,i))*tgrd(n,j,i) + &
                   sigf(n,j,i)*taf(n,j,i) - sts(n,j,i)
          u1 = wtur + d_two*dsqrt(dthdz)
          ribd(n,j,i) = usw(j,i)**2 + vsw(j,i)**2 + u1**2
        else
          u2 = wtur
          ribd(n,j,i) = usw(j,i)**2 + vsw(j,i)**2 + u2**2
        end if
        vspda(n,j,i) = dsqrt(ribd(n,j,i))
        if ( vspda(n,j,i) < d_one ) then
          vspda(n,j,i) = d_one
          ribd(n,j,i) = d_one
        end if
        rib(n,j,i) = ribn(n,j,i)/ribd(n,j,i)
        !=======================================================================
        ! 3.   obtain drag coefficient as product of neutral value
        !      and stability correction
        !=======================================================================
        ! -0.4 < rib < 0.2   (deardorff, jgr, 1968, 2549-2557)
        if ( rib(n,j,i) < d_zero ) then
          cdr(n,j,i) = cdrn(n,j,i)* &
                   (d_one+24.5D0*dsqrt(-cdrn(n,j,i)*rib(n,j,i)))
        else
          cdr(n,j,i) = cdrn(n,j,i)/(d_one+11.5D0*rib(n,j,i))
        end if
        ! 3.1  apply lower limit to drag coefficient value
        !if ( ldmsk1(n,j,i) == 0 ) then
        !  cdrmin = dmin1(0.25D0*cdrn(n,j,i),6.0D-4)
        !else
          cdrmin = dmax1(0.25D0*cdrn(n,j,i),6.0D-4)
        !end if
        if ( cdr(n,j,i) < cdrmin ) cdr(n,j,i) = cdrmin
        cdrx(n,j,i) = cdr(n,j,i)
      end do
    end do
  end do

  !=======================================================================
  ! 4. obtain drag coefficient over sea ice as weighted average
  !    over ice and leads
  ! warning! the lat test below (4.1-4.3) is model dependent!
  !=======================================================================
 
  ! 4.1  neutral cd over lead water
  do i = ici1 , ici2
    do j = jci1 , jci2
      do n = 1 , nnsg
        if ( ldmsk1(n,j,i) == 2 ) then       !  check each point
          cdrn(n,j,i) = (vonkar/zlgsno(n,j,i))**2
          ! 4.1  drag coefficient over leads
          ribl(n,j,i) = (d_one-271.5D0/sts(n,j,i))* &
                        zh(n,j,i)*egrav/ribd(n,j,i)
          if ( ribl(n,j,i) >= d_zero ) then
            clead(n,j,i) = cdrn(n,j,i)/(d_one+11.5D0*ribl(n,j,i))
          else
            clead(n,j,i) = cdrn(n,j,i)*(d_one+24.5D0* &
                         dsqrt(-cdrn(n,j,i)*ribl(n,j,i)))
          end if
          ! 4.2  calculate weighted avg of ice and lead drag
          !      coefficients
          cdrx(n,j,i) = (d_one-aarea)*cdr(n,j,i) + aarea*clead(n,j,i)
        end if
      end do
    end do
  end do
  end subroutine dragc
!
!=======================================================================
!
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
  real(rk8) :: asigf , cdb , cds , cdv , frab , fras , frav
  integer(ik4) :: n , i , j
!
  call depth
!
  do i = ici1 , ici2
    do j = jci1 , jci2
      do n = 1 , nnsg
        if ( ldmsk1(n,j,i) == 2 ) then
          ! drag coeff over seaice
          sigf(n,j,i) = d_zero
          cdrn(n,j,i) = ( vonkar / zlglnd(n,j,i) )**2
        else if ( ldmsk1(n,j,i) == 1 ) then
          ! drag coeff over land
          frav = sigf(n,j,i)
          asigf = lncl(n,j,i)
          fras = asigf*wt(n,j,i) + (d_one-asigf)*scvk(n,j,i)
          frab = (d_one-asigf)*(d_one-scvk(n,j,i))
          cdb = (vonkar/zlglnd(n,j,i))**2
          cds = (vonkar/zlgsno(n,j,i))**2
          cdv = (vonkar/zlgdis(n,j,i))**2
          cdrn(n,j,i) = frav*cdv + frab*cdb + fras*cds
        else
          ! drag coeff over ocean
          sigf(n,j,i) = d_zero
          cdrn(n,j,i) = (vonkar/zlgocn(n,j,i))**2
        end if
      end do
    end do
  end do
  end subroutine dragdn
!
!=======================================================================
!
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
  real(rk8) :: age
  integer(ik4) :: n , i , j
! 
  do i = ici1 , ici2
    do j = jci1 , jci2
      do n = 1 , nnsg
        if ( ldmsk1(n,j,i) /= 0 ) then
          age = (d_one-d_one/(d_one+snag(n,j,i)))
          rhosw(n,j,i) = 0.10D0*(d_one+d_three*age)
          densi(n,j,i) = 0.01D0/(d_one+d_three*age)
          scrat(n,j,i) = sncv(n,j,i)*densi(n,j,i)
          wt(n,j,i) = d_one
          if ( ldmsk1(n,j,i) /= 2 ) then
            wt(n,j,i) = 0.1D0*scrat(n,j,i)/rough(lveg(n,j,i))
            wt(n,j,i) = wt(n,j,i)/(d_one+wt(n,j,i))
          end if
          sigf(n,j,i) = (d_one-wt(n,j,i))*lncl(n,j,i)
          scvk(n,j,i) = scrat(n,j,i)/(0.1D0+scrat(n,j,i))
        end if
      end do
    end do
  end do
  end subroutine depth
!
end module mod_bats_drag

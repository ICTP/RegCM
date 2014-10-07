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
  use mod_service
  use mod_bats_internal
  use mod_bats_param

  implicit none
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
    implicit none
    real(rk8) :: dthdz , u1 , ribn , zatild , cdrmin
    integer(ik4) :: i
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'dragc'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    !======================================
    !     1.   get neutral drag coefficient
    !======================================

    call dragdn

    do i = ilndbeg , ilndend
      !==================================================
      ! 2.  compute stability as bulk rich. no. = rin/rid
      !     ri(numerator)/ri(denominator)
      !==================================================
      zatild = (zh(i)-displa(lveg(i)))*sigf(i) + &
                zh(i)*(d_one-sigf(i))
      ribn = zatild*egrav*(sts(i)-sigf(i)*taf(i)- &
             (d_one-sigf(i))*tgrd(i))/sts(i)
      !===========================================================
      ! 2.1  compute the bulk richardson number;
      !    first get avg winds to use for ri number by summing the
      !    squares of horiz., vertical, and convective velocities
      !===========================================================
      if ( ribn <= d_zero ) then
        dthdz = (d_one-sigf(i))*tgrd(i) + sigf(i)*taf(i) - sts(i)
        u1 = wtur + d_two*dsqrt(dthdz)
      else
        u1 = wtur
      end if
      ribd(i) = usw(i)**2 + vsw(i)**2 + u1**2
      vspda(i) = dsqrt(ribd(i))
      if ( vspda(i) < d_one ) then
        vspda(i) = d_one
        ribd(i) = d_one
      end if
      rib(i) = ribn/ribd(i)
      !=========================================================
      ! 3.   obtain drag coefficient as product of neutral value
      !      and stability correction
      !=========================================================
      ! -0.4 < rib < 0.2   (deardorff, jgr, 1968, 2549-2557)
      if ( rib(i) < d_zero ) then
        cdr(i) = cdrn(i)*(d_one+24.5D0*dsqrt(-cdrn(i)*rib(i)))
      else
        cdr(i) = cdrn(i)/(d_one+11.5D0*rib(i))
      end if
      ! 3.1  apply lower limit to drag coefficient value
      cdrmin = min(0.25D0*cdrn(i),6.0D-4)
      !end if
      if ( cdr(i) < cdrmin ) cdr(i) = cdrmin
      cdrx(i) = cdr(i)
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
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
    implicit none
    real(rk8) :: asigf , cdb , cds , cdv , frab , fras , frav
    integer(ik4) :: i
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'dragdn'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    call depth
    !
    do i = ilndbeg , ilndend
      ! drag coeff over land
      frav = sigf(i)
      asigf = lncl(i)
      fras = asigf*wt(i) + (d_one-asigf)*scvk(i)
      frab = (d_one-asigf)*(d_one-scvk(i))
      cdb = (vonkar/zlglnd(i))**2
      cds = (vonkar/zlgsno(i))**2
      cdv = (vonkar/zlgdis(i))**2
      cdrn(i) = frav*cdv + frab*cdb + fras*cds
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
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
    implicit none
    real(rk8) :: age , densi
    integer(ik4) :: i
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'depth'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    do i = ilndbeg , ilndend
      age = (d_one-d_one/(d_one+snag(i)))
      densi = 0.01D0/(d_one+d_three*age)
      scrat(i) = sncv(i)*densi
      wt(i) = 0.1D0*scrat(i)/rough(lveg(i))
      wt(i) = wt(i)/(d_one+wt(i))
      sigf(i) = (d_one-wt(i))*lncl(i)
      scvk(i) = scrat(i)/(0.1D0+scrat(i))
      rhosw(i) = 0.10D0*(d_one+d_three*age)
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine depth
!
end module mod_bats_drag
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

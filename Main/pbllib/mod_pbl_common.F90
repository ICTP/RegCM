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

module mod_pbl_common
!
! Storage parameters and constants related to the boundary layer
!
  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_memutil
  use mod_regcm_types

  implicit none

  private

  real(rkx) , public , pointer , dimension(:,:) :: ricr

  integer(ik4) , public :: kmxpbl

  !
  ! Pointers to the TCM state variables
  !
  type(tcm_state) , public :: uwstate
  real(rkx) , public , pointer , dimension(:,:,:,:) :: chiuwten

  integer(ik4) , parameter :: npsi = 1000

  real(rkx) , dimension(:) , pointer :: psim_stab , psim_unstab
  real(rkx) , dimension(:) , pointer  :: psih_stab , psih_unstab

  public :: init_minisfcscheme , minisfcscheme

  contains

  subroutine init_minisfcscheme
    implicit none
    integer(ik4) :: i
    real(rkx) :: zolf
    call getmem1d(psim_stab,0,npsi,'mod_pbl_gfs:psim_stab')
    call getmem1d(psim_unstab,0,npsi,'mod_pbl_gfs:psim_unstab')
    call getmem1d(psih_stab,0,npsi,'mod_pbl_gfs:psih_stab')
    call getmem1d(psih_unstab,0,npsi,'mod_pbl_gfs:psih_unstab')
    !
    ! Precomputation
    !
    psim_stab(0) = d_zero
    psih_stab(0) = d_zero
    psim_unstab(0) = d_zero
    psih_unstab(0) = d_zero
    do i = 1 , npsi
      zolf = real(i,rkx)*0.01_rkx
      psim_stab(i) = psim_stable_full(zolf)
      psih_stab(i) = psih_stable_full(zolf)
      zolf = -zolf
      psim_unstab(i) = psim_unstable_full(zolf)
      psih_unstab(i) = psih_unstable_full(zolf)
    end do
  end subroutine init_minisfcscheme

  subroutine minisfcscheme(ta,pa,tha,za,ua,va,qa,rhoa, &
                           ps,hf,qf,tsk,udrag,hpbl,z0,dx,ldm, &
                           br,psih,psim)
    implicit none
    real(rkx) , intent(in) :: ta , pa , tha , za , ua , va , qa , rhoa
    real(rkx) , intent(in) :: ps , hf , qf , tsk , udrag , hpbl , z0 , dx
    integer(ik4) , intent(in) :: ldm
    real(rkx) , intent(out) :: br , psih , psim

    real(rkx) :: tvcon , rrhox , tskv , xp
    real(rkx) :: wspd0 , zo , thvx , dthvdz , gz1oz0
    real(rkx) :: zol0 , zolzz , dthvm
    real(rkx) :: fluxc , vconv , vsgd , wspd , zol

    tvcon = d_one + ep1*qa
    rrhox = (rgas*(ta*tvcon))/pa

    xp =  (p00/ps)**rovcp
    zo = min(z0,za)
    gz1oz0 = log((za+zo)/zo)
    thvx = tha*tvcon
    tskv = tsk*xp*tvcon
    dthvdz = thvx - tskv
    wspd0 = sqrt(ua*ua + va*va)
    vsgd = 0.32_rkx * (max(dx/5000.0_rkx-d_one,d_zero))**0.33_rkx
    if ( ldm > 0 ) then
      fluxc = max(hf/rhoa*rcpd + qf/rhoa*ep1*tskv,0.0_rkx)
      vconv = d_one*(egrav/tsk*hpbl*fluxc)**0.33_rkx
    else
      if ( -dthvdz >= d_zero ) then
        dthvm = -dthvdz
      else
        dthvm = d_zero
      end if
      vconv = sqrt(dthvm)
    end if
    wspd = sqrt(wspd0*wspd0+vconv*vconv+vsgd*vsgd)
    br = (egrav/(ta*tvcon))*za*dthvdz/(wspd*wspd)
    if ( br > d_zero ) then
      if ( br > 250.0_rkx ) then
        zol = zolri(250.0_rkx,za,zo)
      else
        zol = zolri(br,za,zo)
      end if
    else if ( br < d_zero ) then
      if ( udrag < 0.001_rkx ) then
        zol = br*gz1oz0
      else
        if ( br < -250.0_rkx ) then
          zol = zolri(-250.0_rkx,za,zo)
        else
          zol = zolri(br,za,zo)
        end if
      end if
    else
      zol = d_zero
    end if
    zolzz = zol * (za+zo)/za
    zol0 = zol * zo/za
    if ( br < d_zero ) then
      psim = psim_unstable(zolzz)-psim_unstable(zol0)
      psih = psih_unstable(zolzz)-psih_unstable(zol0)
    else
      psim = psim_stable(zolzz)-psim_stable(zol0)
      psih = psih_stable(zolzz)-psih_stable(zol0)
    end if
    psim = max(min(0.9_rkx*gz1oz0,psim),0.1_rkx*gz1oz0)
    psih = max(min(0.9_rkx*gz1oz0,psih),0.1_rkx*gz1oz0)

    contains

      pure real(rkx) function psim_stable(zolf)
        implicit none
        real(rkx) , intent(in) :: zolf
        integer(ik4) :: nzol
        real(rkx) :: rzol
        nzol = int(zolf*d_100)
        rzol = zolf*100.0_rkx - real(nzol,rkx)
        if ( nzol >= 0 .and. nzol + 1 <= npsi ) then
          psim_stable = psim_stab(nzol) + &
                  rzol*(psim_stab(nzol+1)-psim_stab(nzol))
        else
          if ( zolf >= d_zero ) then
            psim_stable = psim_stable_full(zolf)
          else
            psim_stable = d_zero
          end if
        end if
      end function psim_stable

      pure real(rkx) function psim_unstable(zolf)
        implicit none
        real(rkx) , intent(in) :: zolf
        integer(ik4) :: nzol
        real(rkx) :: rzol
        nzol = int(-zolf*d_100)
        rzol = -zolf*100.0_rkx - real(nzol,rkx)
        if ( nzol >= 0 .and. nzol + 1 <= npsi ) then
          psim_unstable = psim_unstab(nzol) + &
                  rzol*(psim_unstab(nzol+1)-psim_unstab(nzol))
        else
          if ( zolf > d_zero ) then
            psim_unstable = d_zero
          else
            psim_unstable = psim_unstable_full(zolf)
          end if
        end if
      end function psim_unstable

      pure real(rkx) function psih_stable(zolf)
        implicit none
        real(rkx) , intent(in) :: zolf
        integer(ik4) :: nzol
        real(rkx) :: rzol
        nzol = int(zolf*d_100)
        rzol = zolf*100.0_rkx - real(nzol,rkx)
        if ( nzol >= 0 .and. nzol + 1 <= npsi ) then
          psih_stable = psih_stab(nzol) + &
                  rzol*(psih_stab(nzol+1)-psih_stab(nzol))
        else
          if ( zolf > d_zero ) then
            psih_stable = psih_stable_full(zolf)
          else
            psih_stable = d_zero
          end if
        end if
      end function psih_stable

      pure real(rkx) function psih_unstable(zolf)
        implicit none
        real(rkx) , intent(in) :: zolf
        integer(ik4) :: nzol
        real(rkx) :: rzol
        nzol = int(-zolf*d_100)
        rzol = -zolf*100.0_rkx - real(nzol,rkx)
        if ( nzol >= 0 .and. nzol + 1 <= npsi ) then
          psih_unstable = psih_unstab(nzol) + &
                  rzol*(psih_unstab(nzol+1)-psih_unstab(nzol))
        else
          if ( zolf > d_zero ) then
            psih_unstable = d_zero
          else
            psih_unstable = psih_unstable_full(zolf)
          end if
        end if
      end function psih_unstable

      pure real(rkx) function zolri2(zol2,ri2,z,z0)
        implicit none
        real(rkx) , intent(in) :: zol2 , ri2 , z , z0
        real(rkx) :: zol20 , zol3
        real(rkx) :: psix2 , psih2
        zol20 = zol2 * z0/z
        zol3 = zol2 + zol20
        if ( ri2 < d_zero ) then
          psix2 = log((z+z0)/z0)-(psim_unstable(zol3)-psim_unstable(zol20))
          psih2 = log((z+z0)/z0)-(psih_unstable(zol3)-psih_unstable(zol20))
        else
          psix2 = log((z+z0)/z0)-(psim_stable(zol3)-psim_stable(zol20))
          psih2 = log((z+z0)/z0)-(psih_stable(zol3)-psih_stable(zol20))
        end if
        zolri2 = zol2 * psih2/psix2**2 - ri2
      end function zolri2

      pure real(rkx) function zolri(ri,z,z0)
        implicit none
        real(rkx) , intent(in) :: ri , z , z0
        real(rkx) :: x1 , x2
        real(rkx) :: fx1 , fx2
        if ( ri < d_zero ) then
          x1 = -5.0_rkx
          x2 = 0.0_rkx
        else
          x1 = 0.0_rkx
          x2 = 5.0_rkx
        end if
        fx1 = zolri2(x1,ri,z,z0)
        fx2 = zolri2(x2,ri,z,z0)
        zolri = fx1
        do
          if ( abs(fx2) < abs(fx1) ) then
            x1 = x1-fx1/(fx2-fx1)*(x2-x1)
            fx1 = zolri2(x1,ri,z,z0)
            zolri = x1
          else
            x2 = x2-fx2/(fx2-fx1)*(x2-x1)
            fx2 = zolri2(x2,ri,z,z0)
            zolri = x2
          end if
          if ( abs(x1-x2) < 0.01_rkx ) exit
        end do
      end function zolri

  end subroutine minisfcscheme

  pure real(rkx) function psim_stable_full(zolf)
    implicit none
    real(rkx) , intent(in) :: zolf
    psim_stable_full = -6.1_rkx*log(zolf + &
            (d_one+zolf**2.5_rkx)**(d_one/2.5_rkx))
  end function psim_stable_full

  pure real(rkx) function psih_stable_full(zolf)
    implicit none
    real(rkx) , intent(in) :: zolf
    psih_stable_full = -5.3_rkx*log(zolf + &
            (d_one+zolf**1.1_rkx)**(d_one/1.1_rkx))
  end function psih_stable_full

  pure real(rkx) function psim_unstable_full(zolf)
    implicit none
    real(rkx) , intent(in) :: zolf
    real(rkx) :: x , psimk , ym , psimc
    x = (d_one-16.0_rkx*zolf)**0.25_rkx
    psimk = 2.0_rkx * log(0.5_rkx*(d_one+x)) + &
            log(0.5_rkx*(d_one+x*x)) - &
            2.0_rkx * atan(x) + 2.0_rkx * atan(d_one)
    ym = (d_one - 10.0_rkx * zolf)**0.33_rkx
    psimc = (3.0_rkx/2.0_rkx) * log((ym**2+ym+d_one)/3.0_rkx) - &
            sqrt(3.0_rkx)*atan((2.0_rkx*ym+d_one)/sqrt(3.0_rkx)) + &
            4.0_rkx*atan(d_one)/sqrt(3.0_rkx)
    psim_unstable_full = (psimk + zolf**2 * psimc) / (d_one+zolf**2)
  end function psim_unstable_full

  pure real(rkx) function psih_unstable_full(zolf)
    implicit none
    real(rkx) , intent(in) :: zolf
    real(rkx) :: y , psihk , yh , psihc
    y = (d_one - 16.0_rkx*zolf)**0.5_rkx
    psihk = 2.0_rkx*log((d_one+y)/2.0_rkx)
    yh = (d_one-34.0_rkx*zolf)**0.33_rkx
    psihc = (3.0_rkx/2.0_rkx)*log((yh**2+yh+d_one)/3.0_rkx) - &
            sqrt(3.0_rkx)*atan((2.0_rkx*yh+d_one)/sqrt(3.0_rkx)) + &
            4.0_rkx*atan(d_one)/sqrt(3.0_rkx)
    psih_unstable_full = (psihk + zolf**2 * psihc) / (d_one+zolf**2)
  end function psih_unstable_full

end module mod_pbl_common

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

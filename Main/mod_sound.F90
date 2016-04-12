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

module mod_sound

  use mod_realkinds
  use mod_intkinds
  use mod_memutil
  use mod_dynparam
  use mod_runparams
  use mod_mpmessage
  use mod_constants
  use mod_stdio
  use mod_mppparam
  use mod_cu_interface , only : total_precip_points
  use mod_atm_interface

  implicit none

  private

  public :: allocate_mod_sound , init_sound , sound

  real(rk8) , pointer , dimension(:,:,:) :: aa
  real(rk8) , pointer , dimension(:,:,:) :: b
  real(rk8) , pointer , dimension(:,:,:) :: c
  real(rk8) , pointer , dimension(:,:,:) :: rhs
  real(rk8) , pointer , dimension(:,:,:) :: sigdot
  real(rk8) , pointer , dimension(:,:,:) :: wo
  real(rk8) , pointer , dimension(:,:,:) :: ca
  real(rk8) , pointer , dimension(:,:,:) :: g1 , g2
  real(rk8) , pointer , dimension(:,:,:) :: ptend , pxup , pyvp
  real(rk8) , pointer , dimension(:,:,:) :: ucrs , vcrs
  real(rk8) , pointer , dimension(:,:,:) :: tk
  real(rk8) , pointer , dimension(:,:,:) :: cc , cdd , cj
  real(rk8) , pointer , dimension(:,:,:) :: pi
  real(rk8) , pointer , dimension(:,:,:) :: e , f
  real(rk8) , pointer , dimension(:,:) :: astore
  real(rk8) , pointer , dimension(:,:) :: rpsb
  real(rk8) , pointer , dimension(:,:) :: wpval

  real(rk8) , dimension(-6:6) :: fi , fj
  real(rk8) , dimension(0:6) :: fk , fl

  !
  !  BET IS IKAWA BETA PARAMETER (0.=CENTERED, 1.=BACKWARD)
  !
  !  The parameter BET determines the time-weighting, where zero gives a
  !  time-centered average and positive values give a bias towards the
  !  future time step that can be used for acoustic damping. In practice,
  !  values of BET = 0.2 - 0.4 are used (MM5 manual, Sec. 2.5.1). This
  !  time-weighting is applied on w and pp:
  !
  real(rk8) , parameter :: bet = 0.4D0
  real(rk8) , parameter :: xkd = 0.1D0
  real(rk8) , parameter :: xgamma = d_one/(d_one-rovcp)
  real(rk8) :: cs , bp , bm , bpxbm , bpxbp
  real(rk8) :: dtsmax
  real(rk8) :: npts

  contains

#include "cpmf.inc"

  subroutine allocate_mod_sound
    implicit none
    call getmem3d(aa,jci1,jci2,ici1,ici2,1,kzp1,'sound:aa')
    call getmem3d(b,jci1,jci2,ici1,ici2,1,kzp1,'sound:b')
    call getmem3d(c,jci1,jci2,ici1,ici2,1,kzp1,'sound:c')
    call getmem3d(rhs,jci1,jci2,ici1,ici2,1,kzp1,'sound:rhs')
    call getmem3d(sigdot,jci1,jci2,ici1,ici2,1,kzp1,'sound:sigdot')
    call getmem3d(wo,jci1,jci2,ici1,ici2,1,kzp1,'sound:wo')
    call getmem3d(e,jci1,jci2,ici1,ici2,1,kzp1,'sound:e')
    call getmem3d(f,jci1,jci2,ici1,ici2,1,kzp1,'sound:f')
    call getmem3d(ca,jci1,jci2,ici1,ici2,1,kz,'sound:ca')
    call getmem3d(g1,jci1,jci2,ici1,ici2,1,kz,'sound:g1')
    call getmem3d(g2,jci1,jci2,ici1,ici2,1,kz,'sound:g2')
    call getmem3d(ptend,jci1,jci2,ici1,ici2,1,kz,'sound:ptend')
    call getmem3d(pxup,jci1,jci2,ici1,ici2,1,kz,'sound:pxup')
    call getmem3d(pyvp,jci1,jci2,ici1,ici2,1,kz,'sound:pyvp')
    call getmem3d(tk,jci1,jci2,ici1,ici2,1,kz,'sound:tk')
    call getmem3d(cc,jci1,jci2,ici1,ici2,1,kz,'sound:cc')
    call getmem3d(cdd,jci1,jci2,ici1,ici2,1,kz,'sound:cdd')
    call getmem3d(cj,jci1,jci2,ici1,ici2,1,kz,'sound:cj')
    call getmem3d(pi,jci1,jci2,ici1,ici2,1,kz,'sound:pi')
    call getmem2d(astore,jci1,jci2,ici1,ici2,'sound:astore')
    call getmem2d(wpval,jci1,jci2,ici1,ici2,'sound:wpval')
    call getmem2d(rpsb,jce1gb,jce2gb,ice1gb,ice2gb,'sound:rpsb')
    call getmem3d(ucrs,jci1,jci2,ici1,ici2,1,kz,'sound:ucrs')
    call getmem3d(vcrs,jci1,jci2,ici1,ici2,1,kz,'sound:vcrs')
  end subroutine allocate_mod_sound

  subroutine init_sound
    implicit none
    integer(ik4) :: i

    if ( ifupr == 1 ) then
      !
      ! DEFINE VALUES OF FK, FL, FI & FJ FOR UPPER RADIATIVE BC
      !
      do i = 1 , 5
        fk(i) = d_two
        fl(i) = d_two
      end do
      fk(0) = d_one
      fl(0) = d_one
      fk(6) = d_one
      fl(6) = d_one
      do i = -5 , 5
        fi(i) = d_one
        fj(i) = d_one
      end do
      fi(-6) = d_half
      fj(-6) = d_half
      fi(6) = d_half
      fj(6) = d_half
    end if
    cs = sqrt(xgamma*rgas*stdt)
    ! Calculate short time-step
    dtsmax = dx/cs/(d_one+xkd)
    bp = (d_one+bet)*d_half
    bm = (d_one-bet)*d_half
    bpxbp = bp*bp
    bpxbm = bp*bm
    npts = dble((nicross-2)*(njcross-2))
  end subroutine init_sound

  subroutine sound
    implicit none
    real(rk8) :: cddtmp ,  cfl , check , chh , cjtmp , cpm , denom , &
      dppdp0 , dpterm , dts , ppold , rho , rho0s , rofac , maxcfl , &
      rll , rkk , ri , rj
    integer(ik4) :: i , j , k , km1 , kp1 , istep , it , iconvec
    logical , save :: cfl_error = .false.
    character (len=32) :: appdat
    !
    ! Variables to implement upper radiative bc
    !
    real(rk8) :: abar , atot , dxmsfb , ensq , rhon , rhontot , xkeff , &
                 xkleff , xleff , xmsftot , xmsf
    real(rk8) :: loc_abar , loc_rhon , loc_xmsf
    integer(ik4) :: inn , jnn , ll , kk , nsi , nsj
    !
    ! HT IS G*(TERR. HT.)
    ! UTENS, VTENS, PPTENS AND WTENS ARE SUPPLIED TO THIS ROUTINE
    ! UTENS,VTENS=(ADVECTION+CORIOLIS+DIFFUSION)
    ! TTENS=G*W*(TLP/R/T0-1/CP)+DP`DT/RHO/CP+Q/CP (ADIABATIC+MEAN ADVECT
    !    DIABATIC)+ADVECTION+DIFFUSION
    ! PPTENS(I,K)=
    !   XGAMMA*PR0(I,K)*QQ(I,K)/CP/T0(I,K)(HEATING)
    !   +ADVECTION+DIFFUSION
    ! WTENS(I,K)=WTL*G*((T(I,K)-T0(I,K))/T0(I,K)-R/CP*PP3D(I,K)/PR0(I,K)
    !   - WTU*G*((T(I,K-1)-T0(I,K-1))/T0(I,K-1)-R/CP*PP3D(I,K-1)/PR0(I,
    !   (BUOYANCY)+ADVECTION+DIFFUSION
    !   WHERE WTL=(SIGMA(K)-SIGMA(K-1))/(SIGMA(K+1)-SIGMA(K-1))
    !   WTU=(SIGMA(K+1)-SIGMA(K))/(SIGMA(K+1)-SIGMA(K-1))
    !   TIME-STEPS(ISTEP)
    !
    !  pp(BET)=0.5(1+BET)*pp(t+1)+0.5(1-BET)*pp(t)
    !   w(BET)=0.5(1+BET)* w(t+1)+0.5(1-BET)* w(t)
    !
    ! DTL LONG TIME-STEP (XXB-XXC)
    if ( ktau == 0 ) then
      istep = 4
    else
      istep = max(int(dt/dtsmax),2)
      if ( ktau > 1 ) then
        istep = max(4,istep)
      end if
    end if
    dts = dt/dble(istep)
    !
    ! Calculate the loop boundaries
    !
    if ( ktau == 0 .or. ktau == 2 ) then
      if ( myid == italk ) write(stdout,'(a,i2,a,f7.2,a,i3,a,f6.3,a,f6.3)') &
            ' ktau = ' , ktau , ' : Short time step ' , dts , &
            ', nstep = ', istep , ', beta = ' , bet , ', xkd = ' , xkd
    end if
    !
    !  Premultiply the tendency arrays by dts
    !
    !  Calculate initial arrays for short timestep
    !  xxb stores filtered old xxa without xxc term
    !  no asselin filter on boundary
    !
    do i = ice1gb , ice2gb
      do j = jce1gb , jce2gb
        rpsb(j,i) = d_one/sfs%psb(j,i)
      end do
    end do
    do k = 1 , kz
      do i = ide1 , ide2
        do j = jde1 , jde2
          aten%u(j,i,k) = aten%u(j,i,k) * dts
          atmc%u(j,i,k) = atm2%u(j,i,k)/sfs%psdotb(j,i)
          atm2%u(j,i,k) = omuhf*atm1%u(j,i,k) + gnuhf*atm2%u(j,i,k)
          aten%v(j,i,k) = aten%v(j,i,k) * dts
          atmc%v(j,i,k) = atm2%v(j,i,k)/sfs%psdotb(j,i)
          atm2%v(j,i,k) = omuhf*atm1%v(j,i,k) + gnuhf*atm2%v(j,i,k)
        end do
      end do
    end do
    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          atmc%qx(j,i,k,iqv) = atm2%qx(j,i,k,iqv) * rpsb(j,i)
        end do
      end do
    end do
    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          aten%pp(j,i,k) = aten%pp(j,i,k) * dts
          atmc%pp(j,i,k) = atm2%pp(j,i,k) * rpsb(j,i)
          atm2%pp(j,i,k) = omuhf*atm1%pp(j,i,k) + gnuhf*atm2%pp(j,i,k)
        end do
      end do
    end do
    do k = 1 , kzp1
      do i = ice1 , ice2
        do j = jce1 , jce2
          aten%w(j,i,k) = aten%w(j,i,k) * dts
          atmc%w(j,i,k) = atm2%w(j,i,k) * rpsb(j,i)
          atm2%w(j,i,k) = omuhf*atm1%w(j,i,k) + gnuhf*atm2%w(j,i,k)
        end do
      end do
    end do
    !
    ! Time Step loop
    !
    do it = 1 , istep
      if ( it > 1 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              atmc%pp(j,i,k) = atmc%pp(j,i,k) + xkd*pi(j,i,k)
            end do
          end do
        end do
      end if
      do k = 1 , kz
        kp1 = min(kz,k+1)
        km1 = max(1,k-1)
        do i = ice1 , ice2
          do j = jce1 , jce2
            atmc%t(j,i,k) = (atmc%pp(j,i,km1)-atmc%pp(j,i,kp1)) / &
                            (atm0%pr(j,i,km1)-atm0%pr(j,i,kp1))
          end do
        end do
      end do
      call exchange(atmc%t,1,jce1,jce2,ice1,ice2,1,kz)
      call exchange(atmc%pp,1,jce1,jce2,ice1,ice2,1,kz)
      !
      ! Advance u and v
      !
      do k = 1 , kz
        do i = idi1 , idi2
          do j = jdi1 , jdi2
            ! Predict u and v
            rho    = d_rfour * (atm1%rho(j,i,k)   + atm1%rho(j-1,i,k) + &
                                atm1%rho(j,i-1,k) + atm1%rho(j-1,i-1,k))
            dppdp0 = d_rfour * (atmc%t(j,i,k)   + atmc%t(j-1,i,k) + &
                                atmc%t(j,i-1,k) + atmc%t(j-1,i-1,k))
            ! Divide by map scale factor
            chh = d_half * dts / (rho*dx) / mddom%msfd(j,i)
            !
            ! Nonhydrostatic model: pressure gradient term in sigma vertical
            ! coordinanate: 4th RHS term in Eqs. 2.2.1, 2.2.2, 2.2.9, 2.2.10,
            ! 2.3.3, 2.3.4 in the MM5 manual.
            !
            atmc%u(j,i,k) = atmc%u(j,i,k) - &
                      chh * (atmc%pp(j,i,k)   - atmc%pp(j-1,i,k)   + &
                             atmc%pp(j,i-1,k) - atmc%pp(j-1,i-1,k) - &
                      ( atm0%pr(j,i,k)   - atm0%pr(j-1,i,k)  + &
                        atm0%pr(j,i-1,k) - atm0%pr(j-1,i-1,k)) * dppdp0)
            atmc%v(j,i,k) = atmc%v(j,i,k) - &
                      chh * (atmc%pp(j,i,k)   - atmc%pp(j,i-1,k)   + &
                             atmc%pp(j-1,i,k) - atmc%pp(j-1,i-1,k) - &
                      ( atm0%pr(j,i,k)   - atm0%pr(j,i-1,k)  + &
                        atm0%pr(j-1,i,k) - atm0%pr(j-1,i-1,k)) * dppdp0)
          end do
        end do
      end do
      do k = 1 , kz
        do i = ide1 , ide2
          do j = jde1 , jde2
            atmc%u(j,i,k) = atmc%u(j,i,k) + aten%u(j,i,k)
            atmc%v(j,i,k) = atmc%v(j,i,k) + aten%v(j,i,k)
          end do
        end do
      end do
      call exchange(atmc%u,1,jde1,jde2,ide1,ide2,1,kz)
      call exchange(atmc%v,1,jde1,jde2,ide1,ide2,1,kz)
      if ( it > 1 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              atmc%pp(j,i,k) = atmc%pp(j,i,k) - xkd*pi(j,i,k)
            end do
          end do
        end do
      end if
      !
      !  Semi-implicit solution for w and p
      !
      !  Presure perturbation tendency: 4th, 5th and 6th RHS terms
      !  in Eq.2.2.4. 4th and 5th RHS terms in Eq.2.3.8
      !  Vertical momentum    tendency: 1st subterm and part of the
      !  3rd subterm of the 4th RHS term in Eq.2.2.3 and 2.2.11
      !  (see Hint below). This is joined into the 4th RHS term in Eq. 2.3.7.
      !  Hint: R=Cp-Cv, gamma=Cp/Cv -> 1/gamma+R/Cp=1
      !
      do k = 1 , kzp1
        do i = ici1 , ici2
          do j = jci1 , jci2
            wo(j,i,k) = atmc%w(j,i,k)
          end do
        end do
      end do
      !
      ! Vertical boundary conditions, w=v.dh/dy at bottom, lid at top
      !
      do i = ici1 , ici2
        do j = jci1 , jci2
          atmc%w(j,i,kzp1) = d_half * d_rfour * regrav * &
                     ((atmc%v(j,i+1,kz)   + atmc%v(j,i,kz) +    &
                       atmc%v(j+1,i+1,kz) + atmc%v(j+1,i,kz)) * &
                      ( mddom%ht(j,i+1) - mddom%ht(j,i-1) ) +   &
                      (atmc%u(j,i+1,kz)   + atmc%u(j,i,kz) +    &
                       atmc%u(j+1,i+1,kz) + atmc%u(j+1,i,kz)) * &
                      ( mddom%ht(j+1,i) - mddom%ht(j-1,i))) /   &
                      ( dx * mddom%msfx(j,i) )
          e(j,i,kz) = d_zero
          f(j,i,kz) = atmc%w(j,i,kzp1)
        end do
      end do
      do i = ici1 , ici2
        do j = jci1 , jci2
          cc(j,i,1)  = xgamma * atm1%pr(j,i,1) * dts/ (dx*mddom%msfx(j,i))
          cdd(j,i,1) = xgamma * atm1%pr(j,i,1) * atm0%rho(j,i,1) * &
                       egrav * dts / (atm0%ps(j,i)*dsigma(1))
          cj(j,i,1)  = d_half * atm0%rho(j,i,1) * egrav * dts
          pxup(j,i,1) = 0.0625D0 *                              &
                      ( atm0%pr(j+1,i,1) - atm0%pr(j-1,i,1) ) * &
                      ( atmc%u(j,i,1)   + atmc%u(j+1,i,1)   +   &
                        atmc%u(j,i+1,1) + atmc%u(j+1,i+1,1) -   &
                        atmc%u(j,i,2)   - atmc%u(j+1,i,2)   -   &
                        atmc%u(j,i+1,2) - atmc%u(j+1,i+1,2) ) / &
                      ( atm0%pr(j,i,1) - atm0%pr(j,i,2))
          pyvp(j,i,1) = 0.0625D0 *                              &
                      ( atm0%pr(j,i+1,1) - atm0%pr(j,i-1,1) ) * &
                      ( atmc%v(j,i,1)   + atmc%v(j+1,i,1)   +   &
                        atmc%v(j,i+1,1) + atmc%v(j+1,i+1,1) -   &
                        atmc%v(j,i,2)   - atmc%v(j+1,i,2)   -   &
                        atmc%v(j,i+1,2) - atmc%v(j+1,i+1,2) ) / &
                      ( atm0%pr(j,i,1) - atm0%pr(j,i,2) )
          !
          ! Zero gradient (free slip) b.c.s on v at top and bottom
          !
          ! IG: at the top (k=1), w(x,y,1)=0, dw(x,y,1)/dsigma=0 so
          ! 3rd and 4th LHS in Eq. 2.5.1.4 vanish.
          !
          ptend(j,i,1) = aten%pp(j,i,1) - d_half * cc(j,i,1) *       &
                       ( ( atmc%v(j,i+1,1)   * mddom%msfd(j,i+1)   - &
                           atmc%v(j,i,1)     * mddom%msfd(j,i)     + &
                           atmc%v(j+1,i+1,1) * mddom%msfd(j+1,i+1) - &
                           atmc%v(j+1,i,1)   * mddom%msfd(j+1,i)   + &
                           atmc%u(j+1,i,1)   * mddom%msfd(j+1,i)   - &
                           atmc%u(j,i,1)     * mddom%msfd(j,i)     + &
                           atmc%u(j+1,i+1,1) * mddom%msfd(j+1,i+1) - &
                           atmc%u(j,i+1,1)   * mddom%msfd(j,i+1) ) / &
                       mddom%msfx(j,i) - d_two * (pyvp(j,i,1) + pxup(j,i,1)) )
          tk(j,i,1) = (d_half * atm0%ps(j,i) * atm0%t(j,i,1)) / &
                      (xgamma * atm0%pr(j,i,1) * atm2%t(j,i,1) * rpsb(j,i))
        end do
      end do
      do k = 2 , kz
        kp1 = min(k+1,kz)
        km1 = k-1
        do i = ici1 , ici2
          do j = jci1 , jci2
            tk(j,i,k) = (d_half * atm0%ps(j,i) * atm0%t(j,i,k)) / &
                        (xgamma * atm0%pr(j,i,k) * atm2%t(j,i,k) * rpsb(j,i))
            rofac = (dsigma(km1)*atm0%rho(j,i,k) +    &
                     dsigma(k)  *atm0%rho(j,i,km1)) / &
                    (dsigma(km1)*atm1%rho(j,i,k) +    &
                     dsigma(k)  *atm1%rho(j,i,km1))
            !
            ! Set factors for differencing
            !
            cc(j,i,k)  = xgamma * atm1%pr(j,i,k) * dts / (dx*mddom%msfx(j,i))
            cdd(j,i,k) = xgamma * atm1%pr(j,i,k) * atm0%rho(j,i,k) * &
                         egrav * dts / (atm0%ps(j,i)*dsigma(k))
            cj(j,i,k) = d_half * atm0%rho(j,i,k) * egrav * dts
            ca(j,i,k) = egrav * dts / (atm0%pr(j,i,k)-atm0%pr(j,i,km1)) * rofac
            g1(j,i,k) = d_one - dsigma(km1) * tk(j,i,k)
            g2(j,i,k) = d_one + dsigma(k) * tk(j,i,km1)
            !
            ! Implicit w equation coefficient arrays and rhs (ikawa method)
            !
            c(j,i,k) = -ca(j,i,k) * (cdd(j,i,km1)-cj(j,i,km1))*g2(j,i,k)*bpxbp
            b(j,i,k) = d_one + ca(j,i,k) * ( g1(j,i,k) *      &
                       (cdd(j,i,k) - cj(j,i,k)) + g2(j,i,k) * &
                       (cdd(j,i,km1) + cj(j,i,km1)) ) * bpxbp
            aa(j,i,k) = -ca(j,i,k) * (cdd(j,i,k)+cj(j,i,k))*g1(j,i,k)*bpxbp
            pyvp(j,i,k) = 0.125D0 * (atm0%pr(j,i+1,k) - atm0%pr(j,i-1,k)) * &
                          ( atmc%v(j,i,km1)   + atmc%v(j+1,i,km1)   +       &
                            atmc%v(j,i+1,km1) + atmc%v(j+1,i+1,km1) -       &
                            atmc%v(j,i,kp1)   - atmc%v(j+1,i,kp1)   -       &
                            atmc%v(j,i+1,kp1) - atmc%v(j+1,i+1,kp1) ) /     &
                          ( atm0%pr(j,i,km1) - atm0%pr(j,i,kp1) )
            pxup(j,i,k) = 0.125D0 * (atm0%pr(j+1,i,k) - atm0%pr(j-1,i,k)) * &
                          ( atmc%u(j,i,km1)   + atmc%u(j+1,i,km1)   +       &
                            atmc%u(j,i+1,km1) + atmc%u(j+1,i+1,km1) -       &
                            atmc%u(j,i,kp1)   - atmc%u(j+1,i,kp1)   -       &
                            atmc%u(j,i+1,kp1) - atmc%u(j+1,i+1,kp1) ) /     &
                          ( atm0%pr(j,i,km1) - atm0%pr(j,i,kp1) )
          end do
        end do
      end do
      !
      ! Zero gradient (free slip) b.c.s on v at top and bottom
      !
      ! IG: at the bottom (k=kz), w(x,y,kz)=0, dw(x,y,kz)/dsigma=0 so
      ! 3rd and 4th LHS in Eq. 2.5.1.4 vanish.
      !
      do i = ici1 , ici2
        do j = jci1 , jci2
          pyvp(j,i,kz) = pyvp(j,i,kz)*d_half
          pxup(j,i,kz) = pxup(j,i,kz)*d_half
        end do
      end do
      do k = 2 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            !
            ! Nonhydrostatic model.
            ! Presure perturbation tendency: 5th RHS terms in Eq.2.3.8
            !
            ptend(j,i,k) = aten%pp(j,i,k) - d_half * cc(j,i,k) *         &
                           ( (atmc%v(j,i+1,k)   * mddom%msfd(j,i+1)   -  &
                              atmc%v(j,i,k)     * mddom%msfd(j,i)     +  &
                              atmc%v(j+1,i+1,k) * mddom%msfd(j+1,i+1) -  &
                              atmc%v(j+1,i,k)   * mddom%msfd(j+1,i)   +  &
                              atmc%u(j+1,i,k)   * mddom%msfd(j+1,i)   -  &
                              atmc%u(j,i,k)     * mddom%msfd(j,i)     +  &
                              atmc%u(j+1,i+1,k) * mddom%msfd(j+1,i+1) -  &
                              atmc%u(j,i+1,k)   * mddom%msfd(j,i+1) ) /  &
                          mddom%msfx(j,i) - &
                          d_two*( pyvp(j,i,k) + pxup(j,i,k) ) )
            rhs(j,i,k) = atmc%w(j,i,k) + &
                    aten%w(j,i,k) + ca(j,i,k) * ( bpxbm *   &
                     ( (cdd(j,i,k-1) - cj(j,i,k-1))*g2(j,i,k)*wo(j,i,k-1) - &
                     ( (cdd(j,i,k-1) + cj(j,i,k-1))*g2(j,i,k) +             &
                       (cdd(j,i,k) - cj(j,i,k))*g1(j,i,k) ) * wo(j,i,k) +   &
                       (cdd(j,i,k) + cj(j,i,k))*g1(j,i,k)*wo(j,i,k+1) ) +   &
                     ( atmc%pp(j,i,k)   * g1(j,i,k)   -                     &
                       atmc%pp(j,i,k-1) * g2(j,i,k) ) +                     &
                     ( g1(j,i,k)*ptend(j,i,k) -                             &
                       g2(j,i,k)*ptend(j,i,k-1) ) * bp )
          end do
        end do
      end do
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            pi(j,i,k) = atmc%pp(j,i,k)
            !
            ! Nonhydrostatic model.
            ! Presure perturbation tendency: 4th RHS term and last subterm
            ! in 5th RHS term in Eq. 2.3.8. Also, cf. Eq. 2.5.1.4
            !
            atmc%pp(j,i,k) = atmc%pp(j,i,k) + ptend(j,i,k) +         &
                          ( cj(j,i,k)  * (wo(j,i,k+1) + wo(j,i,k)) + &
                            cdd(j,i,k) * (wo(j,i,k+1) - wo(j,i,k)) ) * bm
          end do
        end do
      end do
      !
      ! Upward calculation of coefficients
      !
      do k = kz , 2 , -1
        do i = ici1 , ici2
          do j = jci1 , jci2
            denom = aa(j,i,k)*e(j,i,k) + b(j,i,k)
            e(j,i,k-1) = -c(j,i,k) / denom
            f(j,i,k-1) = (rhs(j,i,k) - f(j,i,k)*aa(j,i,k)) / denom
          end do
        end do
      end do
      !
      ! First, set upper boundary condition, either w=0 or radiation
      !
      do i = ici1 , ici2
        do j = jci1 , jci2
          wpval(j,i) = d_zero
        end do
      end do
      !
      ! Upper radiative BC, compute the wpval here as in 2.7
      !
      if ( ifupr == 1 ) then
        do i = ici1 , ici2
          do j = jci1 , jci2
            denom = (cdd(j,i,1) + cj(j,i,1)) * bp
            estore(j,i) = atmc%pp(j,i,1) + f(j,i,1) * denom
            astore(j,i) = denom * e(j,i,1) + (cj(j,i,1) - cdd(j,i,1)) * bp
          end do
        end do
        call grid_fill(estore,estore_g,jci1,jci2,ici1,ici2)
        !
        ! If first time through and upper radiation b.c`s are used
        ! Need to calc some coefficients
        !
        if ( ( ktau == 0 .or. mod(ktau,kday) == 0 ) .and. it == 1 ) then
          ! Calculating means for upper radiative boundary conditions
          if ( myid == italk ) then
            write(stdout,'(a,i8)') &
               ' Updating upper radiative BC coefficients at ktau =',ktau
          end if
          atot = d_zero
          rhontot = d_zero
          xmsftot = d_zero
          do i = ici1 , ici2
            do j = jci1 , jci2
              atot = atot + astore(j,i)
              ensq = egrav*egrav/cpd/(atm2%t(j,i,1) * rpsb(j,i))
              rhontot = rhontot + atm1%rho(j,i,1)*sqrt(ensq)
              xmsftot = xmsftot + mddom%msfx(j,i)
            end do
          end do
          loc_abar = atot
          loc_rhon = rhontot
          loc_xmsf = xmsftot
          call sumall(loc_abar,abar)
          call sumall(loc_rhon,rhon)
          call sumall(loc_xmsf,xmsf)
          abar = abar/npts
          rhon = rhon/npts
          xmsf = xmsf/npts
          dxmsfb = d_two/dx/xmsf
          tmask(:,:) = d_zero
          do ll = 0 , 6
            rll = dble(ll)
            do kk = 0 , 6
              rkk = dble(kk)
              xkeff = dxmsfb * sin(mathpi*rkk/12.0D0)*cos(mathpi*rll/12.0D0)
              xleff = dxmsfb * sin(mathpi*rll/12.0D0)*cos(mathpi*rkk/12.0D0)
              xkleff = sqrt(xkeff*xkeff + xleff*xleff)
              do i = -6 , 6
                ri = dble(i)
                do j = -6 , 6
                  rj = dble(j)
                  tmask(j,i) = tmask(j,i) + &
                               fi(i)*fj(j)*fk(kk)*fl(ll)/144.0D0 * &
                               cos(2.0D0*mathpi*rkk*ri/12.0D0) *   &
                               cos(2.0D0*mathpi*rll*rj/12.0D0) *   &
                               xkleff / (rhon-abar*xkleff)
                end do
              end do
            end do
          end do
          ! Finished initial coefficient compute (goes in SAV file)
        end if
        !
        ! Apply upper rad cond. (not in the lateral boundary)
        !
        do i = ici1 , ici2
          do j = jci1 , jci2
            do nsi = -6 , 6
              inn = min(max(i,2),nicross-2)
              do nsj = -6 , 6
                jnn = min(max(j,2),njcross-2)
                wpval(j,i) = wpval(j,i) + &
                       estore_g(jnn,inn)*tmask(nsj,nsi)*wtbdy_g(jnn,inn)
              end do
            end do
          end do
        end do
      end if
      !
      ! Finished calc of radiation w, apply whichever
      !
      do i = ici1 , ici2
        do j = jci1 , jci2
          atmc%w(j,i,1) = wpval(j,i)
        end do
      end do
      !
      ! Downward sweep calculation of w
      !
      do i = ici1 , ici2
        do j = jci1 , jci2
          atmc%w(j,i,2) = bet * atmc%w(j,i,1) + &
                (d_one - bet) * (e(j,i,1)*atmc%w(j,i,1)+f(j,i,1))
        end do
      end do
      do k = 2 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            atmc%w(j,i,k+1) = e(j,i,k)*atmc%w(j,i,k) + f(j,i,k)
          end do
        end do
      end do
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            ucrs(j,i,k) = atmc%u(j,i,k) * mddom%msfd(j,i) +     &
                          atmc%u(j,i+1,k) * mddom%msfd(j,i+1) + &
                          atmc%u(j+1,i,k) * mddom%msfd(j+1,i) + &
                          atmc%u(j+1,i+1,k) * mddom%msfd(j+1,i+1)
            vcrs(j,i,k) = atmc%v(j,i,k) * mddom%msfd(j,i) +     &
                          atmc%v(j,i+1,k) * mddom%msfd(j,i+1) + &
                          atmc%v(j+1,i,k) * mddom%msfd(j+1,i) + &
                          atmc%v(j+1,i+1,k) * mddom%msfd(j+1,i+1)
          end do
        end do
      end do
      cfl = d_zero
      do k = kz , 2 , -1
        do i = ici1 , ici2
          do j = jci1 , jci2
            rho0s = twt(k,1)*atm0%rho(j,i,k) + twt(k,2)*atm0%rho(j,i,k-1)
            sigdot(j,i,k) = -rho0s*egrav*atmc%w(j,i,k)/atm0%ps(j,i) -  &
               sigma(k) * ( dpsdxm(j,i) * ( twt(k,1)*ucrs(j,i,k) +     &
                                            twt(k,2)*ucrs(j,i,k-1) ) + &
                            dpsdym(j,i) * ( twt(k,1)*vcrs(j,i,k) +     &
                                            twt(k,2)*vcrs(j,i,k-1) ) )
            check = abs(sigdot(j,i,k)) * dt / (dsigma(k) + dsigma(k-1))
            cfl = max(check,cfl)
          end do
        end do
      end do
      if ( mod(ktau,krep) == 0 .and. ktau > 0 .and. it == istep ) then
        call sumall(total_precip_points,iconvec)
        call maxall(cfl,maxcfl)
        if ( myid == italk ) then
          appdat = tochar(idatex)
          write(stdout,'(a,a23,a,i16)') ' $$$ ', appdat , ', ktau   = ', ktau
          write(stdout,'(a,e12.5)') ' $$$ max value of CFL = ',maxcfl
          if ( any(icup > 0) ) then
            write(stdout,'(a,i7)') &
              ' $$$  no. of points with active convection = ', iconvec
          end if
        end if
      end if
      if ( cfl > d_one ) then
        do k = kz , 2 , -1
          do i = ici1 , ici2
            do j = jci1 , jci2
              cfl = abs(sigdot(j,i,k)) * dt / (dsigma(k)+dsigma(k-1))
              if ( cfl > d_one ) then
                write(stderr,99003) cfl , atmc%w(j,i,k) , i , j , k
    99003       format ('CFL>1: CFL = ',f12.4,' W = ',f12.4,'  I = ',i5, &
                        '  J = ',i5,'  K = ',i5 )
                if ( cfl_error ) then
                  call fatal(__FILE__,__LINE__,'CFL violation')
                end if
                cfl_error = .true.
              end if
            end do
          end do
        end do
      else
        cfl_error = .false.
      end if
      !
      ! Now compute the new pressure
      !
      do k =  1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            ppold = pi(j,i,k)
            cddtmp = xgamma * atm1%pr(j,i,k) * atm0%rho(j,i,k) * &
                     egrav * dts / (atm0%ps(j,i)*dsigma(k))
            cjtmp = atm0%rho(j,i,k) * egrav * dts * d_half
            atmc%pp(j,i,k) = atmc%pp(j,i,k) + &
                          ( cjtmp  * (atmc%w(j,i,k+1) + atmc%w(j,i,k)) + &
                            cddtmp * (atmc%w(j,i,k+1) - atmc%w(j,i,k)) ) * bp
            pi(j,i,k) = atmc%pp(j,i,k) - ppold - aten%pp(j,i,k)
            !
            ! Compute pressure dp`/dt correction to the temperature
            !
            cpm = cpmf(atmc%qx(j,i,k,iqv))
            dpterm = sfs%psb(j,i)*(atmc%pp(j,i,k)-ppold) / (cpm*atm1%rho(j,i,k))
            atm2%t(j,i,k) = atm2%t(j,i,k) + gnuhf*dpterm
            atm1%t(j,i,k) = atm1%t(j,i,k) + dpterm
          end do
        end do
      end do
      ! Fix bdy pn pp
      if ( ma%has_bdyleft ) then
        do i = ice1 , ice2
          atmc%pp(jce1,i,:) = atmc%pp(jce1,i,:) + aten%pp(jce1,i,:)
        end do
      end if
      if ( ma%has_bdyright ) then
        do i = ice1 , ice2
          atmc%pp(jce2,i,:) = atmc%pp(jce2,i,:) + aten%pp(jce2,i,:)
        end do
      end if
      if ( ma%has_bdybottom ) then
        do j = jce1 , jce2
          atmc%pp(j,ice1,:) = atmc%pp(j,ice1,:) + aten%pp(j,ice1,:)
        end do
      end if
      if ( ma%has_bdytop ) then
        do j = jce1 , jce2
          atmc%pp(j,ice2,:) = atmc%pp(j,ice2,:) + aten%pp(j,ice2,:)
        end do
      end if
      ! Zero gradient conditions on w
      if ( ma%has_bdyleft ) then
        do i = ici1 , ici2
          atmc%w(jce1,i,:) = atmc%w(jci1,i,:)
        end do
        if ( ma%has_bdytop ) then
          atmc%w(jce1,ice2,:) = atmc%w(jci1,ici2,:)
        end if
        if ( ma%has_bdybottom ) then
          atmc%w(jce1,ice1,:) = atmc%w(jci1,ici1,:)
        end if
      end if
      if ( ma%has_bdyright ) then
        do i = ici1 , ici2
          atmc%w(jce2,i,:) = atmc%w(jci2,i,:)
        end do
        if ( ma%has_bdytop ) then
          atmc%w(jce2,ice2,:) = atmc%w(jci2,ici2,:)
        end if
        if ( ma%has_bdybottom ) then
          atmc%w(jce2,ice1,:) = atmc%w(jci2,ici1,:)
        end if
      end if
      if ( ma%has_bdytop ) then
        do j = jci1 , jci2
          atmc%w(j,ice2,:) = atmc%w(j,ici2,:)
        end do
      end if
      if ( ma%has_bdybottom ) then
        do j = jci1 , jci2
          atmc%w(j,ice1,:) = atmc%w(j,ici1,:)
        end do
      end if
      ! End of time loop
    end do
    !
    ! Transfer xxa to xxb, new values to xxa and apply time filter
    !
    do k = 1 , kz
      do i = idi1 , idi2
        do j = jdi1 , jdi2
          atm1%u(j,i,k) = sfs%psdotb(j,i) * atmc%u(j,i,k)
          atm1%v(j,i,k) = sfs%psdotb(j,i) * atmc%v(j,i,k)
          atm2%u(j,i,k) = atm2%u(j,i,k) + gnuhf*atm1%u(j,i,k)
          atm2%v(j,i,k) = atm2%v(j,i,k) + gnuhf*atm1%v(j,i,k)
        end do
      end do
    end do
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          atm1%pp(j,i,k) = sfs%psb(j,i) * atmc%pp(j,i,k)
          atm2%pp(j,i,k) = atm2%pp(j,i,k) + gnuhf*atm1%pp(j,i,k)
        end do
      end do
    end do
    do k = 1 , kzp1
      do i = ici1 , ici2
        do j = jci1 , jci2
          atm1%w(j,i,k) = sfs%psb(j,i) * atmc%w(j,i,k)
          if ( abs(atm1%w(j,i,k)) < dlowval) then
            atm1%w(j,i,k) = sign(dlowval,atm1%w(j,i,k))
          end if
          atm2%w(j,i,k) = atm2%w(j,i,k) + gnuhf*atm1%w(j,i,k)
          if ( abs(atm2%w(j,i,k)) < dlowval) then
            atm2%w(j,i,k) = sign(dlowval,atm2%w(j,i,k))
          end if
        end do
      end do
    end do
  end subroutine sound

end module mod_sound

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

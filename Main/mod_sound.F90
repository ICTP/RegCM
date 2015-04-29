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
  use mod_dynparam
  use mod_runparams
  use mod_constants
  use mod_stdio
  use mod_mppparam
  use mod_cu_interface , only : total_precip_points
  use mod_atm_interface

  implicit none

  private

  public :: sound

  contains

  subroutine sound(dtl,ktau)
    implicit none

    real(rk8) , intent(in) :: dtl
    integer(ik8) , intent(in) :: ktau

    real(rk8) , dimension(jci1:jci2,ici1:ici2,1:kzp1) :: aa
    real(rk8) , dimension(jci1:jci2,ici1:ici2,1:kzp1) :: b
    real(rk8) , dimension(jci1:jci2,ici1:ici2,1:kzp1) :: c
    real(rk8) , dimension(jci1:jci2,ici1:ici2,1:kzp1) :: rhs
    real(rk8) , dimension(jci1:jci2,ici1:ici2,1:kzp1) :: sigdot
    real(rk8) , dimension(jci1:jci2,ici1:ici2,1:kzp1) :: wo
    real(rk8) :: bet , bm , bp , bpxbm , bpxbp , cddtmp ,  cfl , check , &
      chh , cjtmp , cpm , cs , denom , dppdp0 , dpterm , dts , dtsmax ,  &
      ppold , rho , rho0s , rofac , xgamma , xkd , sumcfl
    real(rk8) , dimension(jci1:jci2,ici1:ici2) :: astore , estore
    real(rk8) , dimension(jci1:jci2,ici1:ici2,1:kz) :: ca
    real(rk8) , dimension(jci1:jci2,ici1:ici2,1:kz) :: g1 , g2
    real(rk8) , dimension(jci1:jci2,ici1:ici2,1:kz) :: ptend , pxup , pyvp
    real(rk8) , dimension(jci1:jci2,ici1:ici2,1:kz) :: tk
    integer(ik4) :: i , istep , it , j , k , km1 , kp1 , iconvec
    real(rk8) , dimension(jci1:jci2,ici1:ici2,1:kz) :: cc , cdd , cj
    real(rk8) , dimension(jci1:jci2,ici1:ici2,1:kz) :: pp3d
    real(rk8) , dimension(jci1:jci2,ici1:ici2,1:kz) :: qv3d
    real(rk8) , dimension(jci1:jci2,ici1:ici2,1:kz) :: pi
    real(rk8) , dimension(jci1:jci2,ici1:ici2,1:kz) ::t3d
    real(rk8) , dimension(jdi1:jdi2,idi1:idi2,1:kz) ::u3d , v3d
    real(rk8) , dimension(jci1:jci2,ici1:ici2,1:kzp1) :: e , f
    real(rk8) , dimension(jce1:jce2,ice1:ice2,1:kzp1) :: w3d
    character (len=32) :: appdat
    !
    ! Graziano:
    !
    ! Variables to implement upper radiative bc, for now disabled.
    ! Need deeper knowledge of this stuff...
    !
    ! real(rk8) :: abar , atot , dxmsfb , ensq , rhon , rhontot , xkeff , &
    !   xkleff , xleff , xmsfbar , xmsftot
    ! real(rk8) , dimension(-6:6) :: fi , fj
    ! real(rk8) , dimension(0:6) :: fk , fl
    ! real(rk8) , dimension(jci1:jci2,ici1:ici2) :: wpval
    ! integer(ik4) :: icut , inn , jnn , ll , lp1 , mp1 , npts , nsi , nsj
    ! real(rk8) , dimension(-6:6,-6:6) :: tmask
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
    ! BET IS IKAWA BETA PARAMETER (0.=CENTERED, 1.=BACKWARD)
    !
    bet = 0.4D0
    bp = (d_one+bet)*d_half
    bm = (d_one-bet)*d_half
    bpxbp = bp*bp
    bpxbm = bp*bm
    xkd = 0.1D0

    xgamma = d_one/(d_one-rovcp)
    ! CALCULATE SHORT TIME-STEP
    cs = sqrt(xgamma*rgas*stdt)
    dtsmax = dx/cs/(d_one+xkd)
    ! DTL LONG TIME-STEP (XXB-XXC)
    istep = int(dtl/dtsmax) + 1
    if ( ktau >= 1 ) istep = max(4,istep)
    dts = dtl/istep
    !
    ! Calculate the loop boundaries
    !
!    if ( ktau == 0 ) then
!      write(stdout,*) 'SHORT TIME STEP ' , dts , istep , &
!                      ' BETA = ' , bet , ' XKD = ' , xkd
!      do j = -6 , 6
!        do i = -6 , 6
!          tmask(j,i) = d_zero
!        end do
!      end do
!      if ( ifupr == 1 ) then
!        !
!        ! DEFINE VALUES OF FK, FL, FI & FJ FOR UPPER RADIATIVE BC
!        !
!        do i = 1 , 5
!          fk(i) = d_two
!          fl(i) = d_two
!        end do
!        fk(0) = d_one
!        fl(0) = d_one
!        fk(6) = d_one
!        fl(6) = d_one
!        do i = -5 , 5
!          fi(i) = d_one
!          fj(i) = d_one
!        end do
!        fi(-6) = d_half
!        fj(-6) = d_half
!        fi(6) = d_half
!        fj(6) = d_half
!      end if
!    end if
    !
    !  Premultiply the tendency arrays by dts
    !
    !  Calculate initial arrays for short timestep
    !  xxb stores filtered old xxa without xxc term
    !  no asselin filter on boundary
    !
    do k = 1 , kz
      do i = idi1 , idi2
        do j = jdi1 , jdi2
          aten%u(j,i,k) = aten%u(j,i,k)*dts
          aten%v(j,i,k) = aten%v(j,i,k)*dts
          u3d(j,i,k)    = atm2%u(j,i,k)/sfs%psdotb(j,i)
          v3d(j,i,k)    = atm2%v(j,i,k)/sfs%psdotb(j,i)
          atm2%u(j,i,k) = omuhf*atm1%u(j,i,k)/mddom%msfd(j,i) + &
                          gnuhf*atm2%u(j,i,k)
          atm2%v(j,i,k) = omuhf*atm1%v(j,i,k)/mddom%msfd(j,i) + &
                          gnuhf*atm2%v(j,i,k)
        end do
      end do
    end do
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          aten%pp(j,i,k) = aten%pp(j,i,k)*dts
          qv3d(j,i,k)    = atm2%qx(j,i,k,iqv)/sfs%psa(j,i)
          pp3d(j,i,k)    = atm2%pp(j,i,k)/sfs%psa(j,i)
          atm2%pp(j,i,k) = omuhf*atm1%pp(j,i,k) + gnuhf*atm2%pp(j,i,k)
        end do
      end do
    end do
    do k = 1 , kzp1
      do i = ici1 , ici2
        do j = jci1 , jci2
          aten%w(j,i,k) = aten%w(j,i,k)*dts
          w3d(j,i,k)    = atm2%w(j,i,k)/sfs%psa(j,i)
          atm2%w(j,i,k) = omuhf*atm1%w(j,i,k) + gnuhf*atm2%w(j,i,k)
        end do
      end do
    end do
    !
    ! Time Step loop
    !
    do it = 1 , istep
      if ( it /= 1 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              pp3d(j,i,k) = pp3d(j,i,k) + xkd*pi(j,i,k)
            end do
          end do
        end do
      end if
      do k = 1 , kz
        kp1 = min(kz,k+1)
        km1 = max(1,k-1)
        do i = ici1 , ici2
          do j = jci1 , jci2
            t3d(j,i,k) = (pp3d(j,i,km1)-pp3d(j,i,kp1)) / &
                         (atm0%pr(j,i,km1)-atm0%pr(j,i,kp1))
          end do
        end do
      end do
      !
      ! Advance u and v
      !
      do k = 1 , kz
        do i = idii1 , idii2
          do j = jdii1 , jdii2
            ! Predict u and v
            rho    = d_rfour * (atm2%rho(j,i,k)   + atm2%rho(j,i-1,k) + &
                                atm2%rho(j-1,i,k) + atm2%rho(j-1,i-1,k))
            dppdp0 = d_rfour * (t3d(j,i,k)   + t3d(j,i-1,k) + &
                                t3d(j-1,i,k) + t3d(j-1,i-1,k))
            ! Divide by map scale factor
            chh = d_half * dts / (rho*dx) / mddom%msfd(j,i)
            !
            ! Nonhydrostatic model: pressure gradient term in sigma vertical
            ! coordinanate: 4th RHS term in Eqs. 2.2.1, 2.2.2, 2.2.9, 2.2.10,
            ! 2.3.3, 2.3.4 in the MM5 manual.
            !
            u3d(j,i,k) = u3d(j,i,k) - &
                      chh * (pp3d(j,i,k)   - pp3d(j,i-1,k)   + &
                             pp3d(j-1,i,k) - pp3d(j-1,i-1,k) - &
                      ( atm0%pr(j,i,k)   - atm0%pr(j,i-1,k)  + &
                        atm0%pr(j-1,i,k) - atm0%pr(j-1,i-1,k)) * dppdp0)
            v3d(j,i,k) = v3d(j,i,k) - &
                      chh * (pp3d(j,i,k)   - pp3d(j-1,i,k)   + &
                             pp3d(j,i-1,k) - pp3d(j-1,i-1,k) - &
                      ( atm0%pr(j,i,k)   - atm0%pr(j-1,i,k)  + &
                        atm0%pr(j,i-1,k) - atm0%pr(j-1,i-1,k))*dppdp0)
          end do
        end do
      end do
      do k = 1 , kz
        do i = idi1 , idi2
          do j = jdi1 , jdi2
            u3d(j,i,k) = u3d(j,i,k) + aten%u(j,i,k)
            v3d(j,i,k) = v3d(j,i,k) + aten%v(j,i,k)
          end do
        end do
      end do
      if ( it /= 1 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              pp3d(j,i,k) = pp3d(j,i,k) - xkd*pi(j,i,k)
            end do
          end do
        end do
      end if
      !
      !  Semi-implicit solution for w and p
      !
      !  Nonhydrostatic model.
      !  Presure perturbation tendency: 4th, 5th and 6th RHS terms
      !  in Eq.2.2.4. 4th and 5th RHS terms in Eq.2.3.8
      !  Vertical momentum    tendency: 1st subterm and part of the
      !  3rd subterm of the 4th RHS term in Eq.2.2.3 and 2.2.11
      !  (see Hint below). This is joined into the 4th RHS term in Eq. 2.3.7.
      !  Hint: R=Cp-Cv, gamma=Cp/Cv -> 1/gamma+R/Cp=1
      !
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            wo(j,i,k) = w3d(j,i,k)
          end do
        end do
      end do
      !
      ! Vertical boundary conditions, w=v.dh/dy at bottom, lid at top
      !
      do i = ici1 , ici2
        do j = jci1 , jci2
          w3d(j,i,kzp1) = d_rfour * ((v3d(j+1,i,kz) + v3d(j,i,kz) +          &
                                      v3d(j+1,i+1,kz) + v3d(j,i+1,kz)) *     &
                                     ( mddom%ht(j+1,i) - mddom%ht(j-1,i) ) + &
                                     (u3d(j+1,i,kz) + u3d(j,i,kz) +          &
                                      u3d(j+1,i+1,kz) + u3d(j,i+1,kz)) *     &
                                     ( mddom%ht(j,i+1) - mddom%ht(j,i-1))) / &
                           ( d_two * dx * mddom%msfx(j,i) * egrav )
          e(j,i,kz) = d_zero
          f(j,i,kz) = w3d(j,i,kzp1)
          cc(j,i,1) = xgamma * atm2%pr(j,i,1) * dts/ (dx*mddom%msfx(j,i))
          cdd(j,i,1) = xgamma * atm2%pr(j,i,1) * atm0%rho(j,i,1) * &
                       egrav * dts / (atm0%ps(j,i)*dsigma(1))
          cj(j,i,1) = atm0%rho(j,i,1) * egrav * dts / d_two
          pxup(j,i,1) = 0.0625D0 *                              &
                      ( atm0%pr(j,i+1,1) - atm0%pr(j,i-1,1) ) * &
                      ( u3d(j,i,1) + u3d(j+1,i,1) +             &
                        u3d(j,i+1,1) + u3d(j+1,i+1,1) -         &
                        u3d(j,i,2) - u3d(j+1,i,2) -             &
                        u3d(j,i+1,2) - u3d(j+1,i+1,2) ) /       &
                      ( atm0%pr(j,i,1) - atm0%pr(j,i,2))
          pyvp(j,i,1) = 0.0625D0 *                              &
                      ( atm0%pr(j+1,i,1) - atm0%pr(j-1,i,1) ) * &
                      ( v3d(j,i,1) + v3d(j+1,i,1) +             &
                        v3d(j,i+1,1) + v3d(j+1,i+1,1) -         &
                        v3d(j,i,2) - v3d(j+1,i,2) -             &
                        v3d(j,i+1,2) - v3d(j+1,i+1,2) ) /       &
                      ( atm0%pr(j,i,1) - atm0%pr(j,i,2) )
          !
          ! Zero gradient (free slip) b.c.s on v at top and bottom
          !
          ptend(j,i,1) = aten%pp(j,i,1) - d_half * cc(j,i,1) *      &
                       ( ( v3d(j+1,i,1) * mddom%msfd(j+1,i) -       &
                           v3d(j,i,1) * mddom%msfd(j,i) +           &
                           v3d(j+1,i+1,1) * mddom%msfd(j+1,i+1) -   &
                           v3d(j,i+1,1) * mddom%msfd(j,i+1) +       &
                           u3d(j,i+1,1) * mddom%msfd(j,i+1) -       &
                           u3d(j,i,1) * mddom%msfd(j,i) +           &
                           u3d(j+1,i+1,1) * mddom%msfd(j+1,i+1) -   &
                           u3d(j+1,i,1) * mddom%msfd(j+1,i) ) /     &
                       mddom%msfx(j,i) - d_two * (pyvp(j,i,1) + pxup(j,i,1)) )
          tk(j,i,1) = atm0%ps(j,i) * atm0%t(j,i,1) / &
                      (d_two * xgamma * atm0%pr(j,i,1) * &
                      atm2%t(j,i,1) / sfs%psa(j,i))
        end do
      end do
      do k = 2 , kz
        kp1 = min(k+1,kz)
        km1 = k-1
        do i = ici1 , ici2
          do j = jci1 , jci2
            tk(j,i,k) = atm0%ps(j,i) * atm0%t(j,i,k) / &
                        (d_two * xgamma * atm0%pr(j,i,k) * &
                        atm2%t(j,i,k) / sfs%psa(j,i))
            rofac = (dsigma(k-1)*atm0%rho(j,i,k) + &
                     dsigma(k)*atm0%rho(j,i,k-1)) / &
                    (dsigma(k-1)*atm2%rho(j,i,k) + &
                     dsigma(k)*atm2%rho(j,i,k-1))
            !
            ! Set factors for differencing
            !
            cc(j,i,k) = xgamma * atm2%pr(j,i,k) * dts / (dx*mddom%msfx(j,i))
            cdd(j,i,k) = xgamma * atm2%pr(j,i,k) * atm0%rho(j,i,k) * &
                         egrav * dts / (atm0%ps(j,i)*dsigma(k))
            cj(j,i,k) = atm0%rho(j,i,k) * egrav * dts/d_two
            ca(j,i,k) = egrav * dts / (atm0%pr(j,i,k)-atm0%pr(j,i,km1)) * rofac
            g1(j,i,k) = d_one - dsigma(km1) * tk(j,i,k)
            g2(j,i,k) = d_one + dsigma(k) * tk(j,i,km1)
            !
            ! Implicit w equation coefficient arrays and rhs (ikawa method)
            !
            c(j,i,k) = -ca(j,i,k) * (cdd(j,i,k-1)-cj(j,i,k-1))*g2(j,i,k)*bpxbp
            b(j,i,k) = d_one + ca(j,i,k) * ( g1(j,i,k) *      &
                       (cdd(j,i,k) - cj(j,i,k)) + g2(j,i,k) * &
                       (cdd(j,i,k-1) + cj(j,i,k-1)) ) * bpxbp
            aa(j,i,k) = -ca(j,i,k) * (cdd(j,i,k)+cj(j,i,k))*g1(j,i,k)*bpxbp
            pyvp(j,i,k) = 0.125D0 * (atm0%pr(j+1,i,k) - atm0%pr(j-1,i,k)) * &
                          ( v3d(j,i,km1) + v3d(j+1,i,km1) +         &
                            v3d(j,i+1,km1) + v3d(j+1,i+1,km1) -     &
                            v3d(j,i,kp1) - v3d(j+1,i,kp1) -         &
                            v3d(j,i+1,kp1) - v3d(j+1,i+1,kp1) ) /   &
                          ( atm0%pr(j,i,km1) - atm0%pr(j,i,kp1) )
            pxup(j,i,k) = 0.125D0 * (atm0%pr(j,i+1,k) - atm0%pr(j,i-1,k)) * &
                          ( u3d(j,i,km1) + u3d(j+1,i,km1) +         &
                            u3d(j,i+1,km1) + u3d(j+1,i+1,km1) -     &
                            u3d(j,i,kp1) - u3d(j+1,i,kp1) -         &
                            u3d(j,i+1,kp1) - u3d(j+1,i+1,kp1) ) /   &
                          ( atm0%pr(j,i,km1) - atm0%pr(j,i,kp1) )
          end do
        end do
      end do
      !
      ! Zero gradient (free slip) b.c.s on v at top and bottom
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
            ptend(j,i,k) = aten%pp(j,i,k) - d_half * cc(j,i,k) *      &
                           ( (v3d(j+1,i,k) * mddom%msfd(j+1,i) -      &
                              v3d(j,i,k) * mddom%msfd(j,i) +          &
                              v3d(j+1,i+1,k) * mddom%msfd(j+1,i+1) -  &
                              v3d(j,i+1,k) * mddom%msfd(j,i+1) +      &
                              u3d(j,i+1,k) * mddom%msfd(j,i+1) -      &
                              u3d(j,i,k) * mddom%msfd(j,i) +          &
                              u3d(j+1,i+1,k) * mddom%msfd(j+1,i+1) -  &
                              u3d(j+1,i,k) * mddom%msfd(j+1,i) ) /    &
                          mddom%msfx(j,i) - &
                          d_two*( pyvp(j,i,k) + pxup(j,i,k) ) )
            rhs(j,i,k) = w3d(j,i,k) + aten%w(j,i,k) + ca(j,i,k) * ( bpxbm *   &
                       ( (cdd(j,i,k-1) - cj(j,i,k-1))*g2(j,i,k)*wo(j,i,k-1) - &
                        ( (cdd(j,i,k-1) + cj(j,i,k-1))*g2(j,i,k) +            &
                          (cdd(j,i,k) - cj(j,i,k))*g1(j,i,k) ) * wo(j,i,k) +  &
                          (cdd(j,i,k) + cj(j,i,k))*g1(j,i,k)*wo(j,i,k+1) ) +  &
                       ( pp3d(j,i,k) * g1(j,i,k) -                            &
                         pp3d(j,i,k-1) * g2(j,i,k) ) +                        &
                       ( g1(j,i,k)*ptend(j,i,k) -                             &
                         g2(j,i,k)*ptend(j,i,k-1) ) * bp )
          end do
        end do
      end do
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            pi(j,i,k) = pp3d(j,i,k)
            !
            ! Nonhydrostatic model.
            ! Presure perturbation tendency: 4th RHS term and last subterm
            ! in 5th RHS term in Eq. 2.3.8
            !
            pp3d(j,i,k) = pp3d(j,i,k) + ptend(j,i,k) +              &
                          ( cj(j,i,k) * (wo(j,i,k+1) + wo(j,i,k)) + &
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
      do i = ici1 , ici2
        do j = jci1 , jci2
          denom = (cdd(j,i,1)+cj(j,i,1))*bp
          estore(j,i) = pp3d(j,i,1) + f(j,i,1)*denom
          astore(j,i) = denom*e(j,i,1) + (cj(j,i,1)-cdd(j,i,1))*bp
        end do
      end do
      !
      ! If first time through and upper radiation b.c`s are used
      ! Need to calc some coefficients
      !
!      if ( ifupr == 1 .and. ktau == 0 .and. it == 1 ) then
!        ! CALCULATING MEANS FOR UP. RAD. B.C.
!        atot = d_zero
!        rhontot = d_zero
!        xmsftot = d_zero
!        npts = 0
!        do i = ici1 , ici2
!          do j = jci1 , jci2
!            atot = atot + astore(j,i)
!            ensq = egrav*egrav/cpd/(atm2%t(j,i,1)/sfs%psa(j,i))
!            rhontot = rhontot + atm2%rho(j,i,1)*sqrt(ensq)
!            xmsftot = xmsftot + mddom%msfx(j,i)
!          end do
!        end do
!        npts = (jci2-jci1+1)*(ici2-ici1+1)
!        abar = atot/npts
!        rhon = rhontot/npts
!        xmsfbar = xmsftot/npts
!        dxmsfb = d_two/dx/xmsfbar
!        do ll = 0 , 6
!          do k = 0 , 6
!            xkeff = dxmsfb * sin(mathpi*k/12.0D0)*cos(mathpi*ll/12.0D0)
!            xleff = dxmsfb * sin(mathpi*ll/12.0D0)*cos(mathpi*k/12.0D0)
!            xkleff = sqrt(xkeff*xkeff + xleff*xleff)
!            do j = -6 , 6
!              do i = -6 , 6
!                tmask(j,i) = tmask(j,i) + &
!                             fi(i)*fj(j)*fk(k)*fl(ll)/144.0D0 * &
!                             cos(2.0D0*mathpi*k*i/12.0D0) *     &
!                             cos(2.0D0*mathpi*ll*j/12.0D0) *    &
!                             xkleff / (rhon-abar*xkleff)
!              end do
!            end do
!          end do
!        end do
!      end if
      !
      !  Finished initial coefficient compute
      !  Now do downward sweep for w
      !
      ! First, set upper boundary condition, either w=0 or radiation
      !
!      do i = ici1 , ici2
!        do j = jci1 , jci2
!          wpval(j,i) = d_zero
!        end do
!      end do
!      if ( ifupr == 1 ) then
!        !
!        ! Apply upper rad cond. no w3d(top) in lateral sponge
!        !
!        do i = ici1+3 , ici2-3
!          inn = insi(i,nsi)
!          do j = jci1+3 , jci2-3
!            do nsj = -6 , 6
!              jnn = jnsj(j,nsj)
!              do nsi = -6 , 6
!                wpval(j,i) = wpval(j,i) + &
!                      estore(inn,jnn)*tmask(nsi,nsj)*wtij(j,i)
!              end do
!            end do
!          end do
!        end do
!      end if
!      !
!      ! Finished calc of radiation w
!      !
!      do i = ici1 , ici2
!        do j = jci1 , jci2
!          w3d(j,i,1) = wpval(j,i)
!        end do
!      end do
      do i = ici1 , ici2
        do j = jci1 , jci2
          w3d(j,i,1) = d_zero
        end do
      end do
      !
      ! Downward calculation of w
      !
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            w3d(j,i,k+1) = e(j,i,k)*w3d(j,i,k) + f(j,i,k)
          end do
        end do
      end do
      cfl = d_zero
      do k = kz , 2 , -1
        do i = ici1 , ici2
          do j = jci1 , jci2
            rho0s = twt(k,1)*atm0%rho(j,i,k) + twt(k,2)*atm0%rho(j,i,k-1)
            sigdot(j,i,k) = -rho0s*egrav*w3d(j,i,k)/sfs%psa(j,i)*d_r1000 -   &
               sigma(k) * ( dpsdxm(j,i) * ( twt(k,1)*atms%ubx3d(j,i,k) +     &
                                            twt(k,2)*atms%ubx3d(j,i,k-1) ) + &
                            dpsdym(j,i) * ( twt(k,1)*atms%vbx3d(j,i,k) +     &
                                            twt(k,2)*atms%vbx3d(j,i,k-1) ) )
            check = abs(sigdot(j,i,k)) * dtl / (dsigma(k) + dsigma(k-1))
            cfl = max(check,cfl)
          end do
        end do
      end do
      if ( mod(ktau,krep) == 0 ) then
        call sumall(total_precip_points,iconvec)
        call sumall(cfl,sumcfl)
        if ( myid == italk ) then
          appdat = tochar(idatex)
          sumcfl = sumcfl/nproc
          write(stdout,'(a,a23,a,i16)') ' $$$ ', appdat , ', ktau   = ', ktau
          write(stdout,'(a,e12.5)') ' $$$ mean value of CFL = ',sumcfl
          write(stdout,'(a,i7)') ' $$$  no. of points w/convection = ', iconvec
        end if
      end if
      if ( cfl > d_one ) then
        do k = kz , 2 , -1
          do i = ici1 , ici2
            do j = jci1 , jci2
              cfl = abs(sigdot(j,i,k)) * dtl / (dsigma(k)+dsigma(k-1))
              if ( cfl > d_one ) then
                write(stderr,99003) cfl , w3d(j,i,k) , i , j , k
    99003       format ('CFL>1: CFL = ',f12.4,' W = ',f12.4,'  I = ',i5, &
                        '  J = ',i5,'  K = ',i5 )
              end if
            end do
          end do
        end do
      end if
      !
      ! Now compute the new pressure
      !
      do k =  1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            ppold = pi(j,i,k)
            cddtmp = xgamma * atm2%pr(j,i,k) * atm0%rho(j,i,k) * &
                     egrav * dts / (atm0%ps(j,i)*dsigma(k))
            cjtmp = atm0%rho(j,i,k) * egrav * dts/d_two
            pp3d(j,i,k) = pp3d(j,i,k) + &
                          ( cjtmp * (w3d(j,i,k+1) + w3d(j,i,k)) + &
                            cddtmp * (w3d(j,i,k+1) - w3d(j,i,k)) )*bp
            pi(j,i,k) = pp3d(j,i,k) - ppold - aten%pp(j,i,k)
            !
            ! Compute pressure dp`/dt correction to the temperature
            !
            cpm = cpd * (d_one + 0.856D0*qv3d(j,i,k))
            dpterm = sfs%psa(j,i)*(pp3d(j,i,k)-ppold) / (cpm*atm2%rho(j,i,k))
            atm2%t(j,i,k) = atm2%t(j,i,k) + gnuhf*dpterm
            atm1%t(j,i,k) = atm1%t(j,i,k) + dpterm
          end do
        end do
      end do
      !
      ! Zero gradient conditions on w,  specified on pp
      !
      do k =  1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            pp3d(j,i,k) = pp3d(j,i,k) + aten%pp(j,i,k)
          end do
        end do
      end do
      do k =  1 , kzp1
        do i = ici1 , ici2
          do j = jci1 , jci2
            w3d(j,i,k) = w3d(j,i,k) + aten%w(j,i,k)
          end do
        end do
      end do
      if ( ma%has_bdyleft ) then
        do k = 1 , kzp1
          do i = ice1 , ice2
            w3d(jce1,i,k) = w3d(jci1,i,k)
          end do
        end do
      end if
      if ( ma%has_bdyright ) then
        do k = 1 , kzp1
          do i = ice1 , ice2
            w3d(jce2,i,k) = w3d(jci2,i,k)
          end do
        end do
      end if
      if ( ma%has_bdybottom ) then
        do k = 1 , kzp1
          do j = jce1 , jce2
            w3d(j,ice1,k) = w3d(j,ici1,k)
          end do
        end do
      end if
      if ( ma%has_bdytop ) then
        do k = 1 , kzp1
          do j = jce1 , jce2
            w3d(j,ice2,k) = w3d(j,ici2,k)
          end do
        end do
      end if
      ! End of time loop
    end do
    !
    ! Transfer xxa to xxb, new values to xxa and apply time filter
    !
    do k = 1 , kz
      do i = idi1 , idi1
        do j = jdi1 , jdi1
          atm1%u(j,i,k) = sfs%psdotb(j,i)*u3d(j,i,k)
          atm1%v(j,i,k) = sfs%psdotb(j,i)*v3d(j,i,k)
          atm2%u(j,i,k) = atm2%u(j,i,k) + gnuhf*atm1%u(j,i,k)
          atm2%v(j,i,k) = atm2%v(j,i,k) + gnuhf*atm1%v(j,i,k)
        end do
      end do
    end do
    do k = 1 , kz
      do i = ici1 , ici1
        do j = jci1 , jci1
          atm1%pp(j,i,k) = sfs%psa(j,i)*pp3d(j,i,k)
          atm2%pp(j,i,k) = atm2%pp(j,i,k) + gnuhf*atm1%pp(j,i,k)
        end do
      end do
    end do
    do k = 1 , kzp1
      do i = ici1 , ici1
        do j = jci1 , jci1
          atm1%w(j,i,k) = sfs%psa(j,i)*w3d(j,i,k)
          atm2%w(j,i,k) = atm2%w(j,i,k) + gnuhf*atm1%w(j,i,k)
        end do
      end do
    end do
  end subroutine sound

end module mod_sound

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

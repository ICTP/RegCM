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

module mod_nhinterp

  use mod_realkinds
  use mod_intkinds
  use mod_constants
  use mod_stdio

  implicit none

  private

  public :: nhbase , nhinterp , nhpp , nhw

  contains
    !
    ! Compute the nonhydrostatic base state.
    !
    subroutine nhbase(i1,i2,j1,j2,kx,a,ter,ptop,p0,tlp,ts0,ps0,pr0,t0,rho0)
      implicit none
      integer(ik4) , intent(in) :: i1 , i2 , j1 , j2 , kx
      real(rk8) , intent(in) , dimension(:) :: a
      real(rk8) , intent(in) , dimension(:,:) :: ter
      real(rk8) , intent(in) :: ptop , p0 , tlp , ts0
      real(rk8) , intent(out) , dimension(:,:) :: ps0
      real(rk8) , intent(out) , dimension(:,:,:) :: pr0
      real(rk8) , intent(out) , dimension(:,:,:) :: t0
      real(rk8) , intent(out) , dimension(:,:,:) :: rho0
      integer(ik4) :: i , j , k
      real(rk8) :: ac , alnp , b , piso , ziso , pt
      !
      ! Define ps0 from terrain heights and t0 profile.
      !
      pt = ptop * d_1000
      do i = i1 , i2
        do j = j1 , j2
          ! Height is already multiplied by g.
          ac = ter(i,j) / (d_two * tlp * rgas)
          b = ts0 / tlp
          alnp = -b + sqrt(b*b - d_four * ac)
          ps0(j,i) = p0 * exp(alnp) - pt
          ! Define reference state temperature at model points.
          do k = 1 , kx
            pr0(j,i,k) = ps0(j,i) * a(k) + pt
            t0(j,i,k) = max(ts0 + tlp*log(pr0(j,i,k)/p0), tiso)
            rho0(j,i,k) = pr0(j,i,k) / rgas / t0(j,i,k)
          end do
        end do
      end do
    end subroutine nhbase
    !
    ! Interpolate the hydrostatic input to nonhydrostatic coordinate.
    !
    subroutine nhinterp(i1,i2,j1,j2,kxs,f,tv,ter,ps,ps0,sigma,intmeth,a, &
                        ptop,p0,tlp,ts0,tiso)
      implicit none
      integer(ik4) , intent(in) :: i1 , i2 , j1 , j2 , kxs , intmeth
      real(rk8) , intent(in) , dimension(:) :: a
      real(rk8) , intent(in) , dimension(:) :: sigma
      real(rk8) , intent(in) , dimension(:,:) :: ter
      real(rk8) , intent(in) , dimension(:,:,:) :: tv
      real(rk8) , intent(in) :: ptop , tlp , tiso , p0 , ts0
      real(rk8) , intent(in) , dimension(:,:) :: ps
      real(rk8) , intent(in) , dimension(:,:) :: ps0
      real(rk8) , intent(inout) , dimension(:,:,:) :: f
      integer(ik4) :: i , j , k , l , ll
      integer(ik4) :: im , ip , jm , jp
      real(rk8) :: fl , fu , pt , piso , pr0 , alnp , alnqvn
      real(rk8) :: tvav , ziso , zl , zu
      real(rk8) , dimension(j1:j2,i1:i2,1:kxs) :: fn
      real(rk8) , dimension(j1:j2,i1:i2,1:kxs) :: z , z0
      real(rk8) , dimension(1:kxs+1) :: z0q , zq

      piso = p0 * exp((tiso - ts0)/tlp)
      pt = ptop * d_1000
      do k = 1 , kxs
        do i = i1 , i2
          do j = j1 , j2
            pr0 = ps0(j,i) * a(k) + pt
            if ( pr0 < piso ) then
              alnp = log(piso/p0)
              ziso = -(rdry*tlp/(d_two*egrav)*alnp*alnp + rdry*ts0/egrav*alnp)
              z0(j,i,k) = ziso - rdry * tiso / egrav * log(pr0/piso)
            else
              alnp = log(pr0/p0)
              z0(j,i,k) = -(rdry*tlp/(d_two*egrav)*alnp*alnp + &
                rdry*ts0/egrav*alnp)
            end if
          end do
        end do
      end do
      !
      !  Calculate heights of input temperature sounding for interpolation
      !  to nonhydrostatic model levels.
      !
      zq(kxs+1) = ter(j,i)/egrav
      do k = kxs , 1 , -1
        do i = i1 , i2
          do j = j1 , j2
            zq(k) = zq(k+1) - rdry/egrav * tv(j,i,k) * &
              log((sigma(k)*ps(j,i)+pt)/(sigma(k+1)*ps(j,i)+pt))
            z(j,i,k) = d_half * (zq(k) + zq(k+1))
          end do
        end do
      end do
      !
      ! Interpolate from z to z0 levels.
      !
      do i = i1 , i2
        do j = j1 , j2
          do k = 1 , kxs
            do ll = 1 , kxs - 1
              l = ll
              if (z(j,i,l+1) .lt. z0(j,i,k)) exit
            end do
            zu = z(j,i,l)
            zl = z(j,i,l+1)
            if ( intmeth  == 1 ) then
              fu = f(j,i,l)
              fl = f(j,i,l+1)
              fn(j,i,k) = (fu * (z0(j,i,k) - zl ) + &
                           fl * (zu - z0(j,i,k))) / (zu - zl)
            else if ( intmeth == 2 ) then
              f(j,i,l) = max(f(j,i,l), 1.D-10 )
              f(j,i,l+1) = max(f(j,i,l+1), 1.D-10 )
              if ( z0(j,i,k) > zu ) then
                fn(j,i,k) = f(j,i,l)
              else
                fu = log(f(j,i,l  ))
                fl = log(f(j,i,l+1))
                alnqvn = (fu * (z0(j,i,k) - zl ) + &
                          fl * (zu - z0(j,i,k))) / (zu - zl)
                fn(j,i,k) = exp(alnqvn)
              end if
            end if
          end do
          do k = 1 , kxs
            f(j,i,k) = fn(j,i,k)
          end do
        end do
      end do
    end subroutine nhinterp
    !
    ! Compute the pressure perturbation pp.
    !
    subroutine nhpp(i1,i2,j1,j2,kxs,t,pr0,t0,tv,ps,ps0,sigma,pt,pp)
      implicit none
      integer(ik4) , intent(in) :: i1 , i2 , j1 , j2 , kxs
      integer(ik4) :: i , j , k
      real(rk8) , intent(in) , dimension(:) :: sigma
      real(rk8) , intent(in) , dimension(:,:,:) :: t , t0 , tv
      real(rk8) , intent(in) , dimension(:,:,:) :: pr0
      real(rk8) , intent(in) , dimension(:,:) :: ps , ps0
      real(rk8) , intent(in) :: pt
      real(rk8) , intent(out) , dimension(:,:,:) :: pp
      real(rk8) :: aa , bb , cc , check , checkl , checkr , delp0 , p0surf
      real(rk8) :: psp , tk , tkp1 , tvk , tvkp1 , tvpot , wtl , wtu
      !
      !  Calculate p' at bottom and integrate upwards (hydrostatic balance).
      !
      do i = i1 , i2
        do j = j1 , j2
          !
          !  Correct pressure perturbation field (mass) to be consistent
          !  with psa assuming density perturbation remains constant in
          !  lowest half-layer.   Start with pp at surface.
          !
          p0surf = ps0(j,i) + pt
          psp = ps(j,i) - ps0(i,j)
          delp0 = p0surf - pr0(j,i,kxs)
          tvk = tv(j,i,kxs)
          tk = t(j,i,kxs)
          tvpot = (tvk - t0(j,i,kxs)) / tk
          pp(j,i,kxs) = (tvpot*delp0 + psp) / (d_one + delp0/pr0(j,i,kxs))
          do k = kxs - 1 , 1 , -1
            tvkp1 = tvk
            tvk = tv(j,i,k)
            tkp1 = tk
            tk = t(j,i,k)
            wtl = (sigma(k+1) - sigma(k  )) / (sigma(k+2) - sigma(k))
            wtu = (sigma(k+2) - sigma(k+1)) / (sigma(k+2) - sigma(k))
            aa = egrav / (pr0(j,i,k+1) - pr0(j,i,k))
            bb = egrav * wtl / pr0(j,i,k+1) * t0(j,i,k+1) / tkp1
            cc = egrav * wtu / pr0(j,i,k  ) * t0(j,i,k  ) / tK
            tvpot = wtl * ((tvkp1 - t0(j,i,k+1)) / tkp1) + &
              wtu * ((tvk - t0(j,i,k  )) / tk)
            pp(j,i,k) = (egrav * tvpot + pp(j,i,k+1) * &
              (aa - bb)) / (aa + cc)
          end do
        end do
      end do
      !
      !  Do vertical gradient check
      !
      do k = 1 , kxs - 1
        do i = i1 , i2
          do j = j1 , j2
            wtl = (sigma(k+1) - sigma(k  )) / (sigma(k+2) - sigma(k  ))
            wtu = (sigma(k+2) - sigma(k+1)) / (sigma(k+2) - sigma(k  ))
            tvpot = wtl * (tv(j,i,k+1) - t0(j,i,k+1)) / &
               t(j,i,k+1) + wtu * (tv(j,i,k  ) - t0(j,i,k  )) / t(j,i,k  )
            checkl = (pp(j,i,k+1) - pp(j,i,k)) / (pr0(j,i,k+1) - pr0(j,i,k))
            checkr = tvpot - &
                 wtl * t0(j,i,k+1)/t(j,i,k+1) * pp(j,i,k+1)/pr0(j,i,k+1) - &
                 wtu * t0(j,i,k)  /t(j,i,k)   * pp(j,i,k)  /pr0(j,i,k)
            check = checkl + checkr
            if ( abs(check) > 1.D-2 ) then
              write(stderr,'(A,3I4,3g12.6)') &
                'NHPP vert gradient check ',i,j,k,check,checkl,checkr
            end if
          end do
        end do
      end do
    end subroutine nhpp
    !
    ! Compute the nonhydrostatic initial vertical velocity (w) from the
    ! horizontal wind fields (u and v).
    !
    subroutine nhw(i1,i2,j1,j2,kxs,u,v,tv,rho0,ps,psdot,ps0,xmsfx,ter,sigma, &
                   w,wtop,ds,a,dsigma,ptop,p0,tlp,ts0,tiso)
      implicit none
      integer(ik4) , intent(in) :: i1 , i2 , j1 , j2 , kxs
      real(rk8) , intent(in) , dimension(:) :: a , sigma , dsigma
      real(rk8) , intent(in) , dimension(:,:) :: ter , xmsfx
      real(rk8) , intent(in) , dimension(:,:) :: ps , ps0 , psdot
      real(rk8) , intent(in) , dimension(:,:,:) :: rho0 , tv
      real(rk8) , intent(in) , dimension(:,:,:) :: u , v
      real(rk8) , intent(in) :: ds , p0 , ptop , tiso , tlp , ts0
      real(rk8) , intent(out) , dimension(:,:,:) :: w
      real(rk8) , intent(out) , dimension(:,:) :: wtop
      integer(ik4) :: i , j , k , km
      integer(ik4) :: l , ll , lm , lp
      real(rk8) :: pt , alnpq , dx2 , omegal , omegau , ubar , vbar
      real(rk8) :: piso , pr0 , rh0 , ziso , zl , zu , rho
      real(rk8) , dimension(kxs+1) :: omega , omegan , qdt
      real(rk8) , dimension(kxs+1) :: z0q , zq
      real(rk8) , dimension(j1:j2,i1:i2,kxs+1) :: wtmp

      pt = ptop * d_1000
      dx2 = d_two * ds
      piso = p0 * exp( (tiso-ts0) / tlp )
      do i = i1 , i2
        do j = j1 , j2
          qdt(kxs+1) = d_zero
          z0q(kxs+1) = ter(j,i)/egrav
          if ( ps0(j,i) + pt < piso ) THEN
            write(stderr,'(a,f5.1,a,f6.1,a)') &
              'The chosen value of Tiso, ',tiso,' K occurs at ',piso,' hPa.'
            write(stderr,'(A)') &
              'This value of pressure is .GT. the reference surface pressure.'
          end if
          do k = 1 , kxs
            pr0 = ps0(j,i) * sigma(k) + pt
            if ( pr0 < piso ) then
              alnpq = log(piso/p0)
              ziso = -(rdry*tlp / (d_two*egrav) * alnpq * alnpq + &
                rgas * ts0 / egrav * alnpq)
              z0q(k) = ziso - rdry * tiso / egrav * log(pr0/piso)
            else
              alnpq = log((ps0(j,i) * sigma(k) + pt) / p0)
              z0q(k) = -(rdry*tlp / (d_two*egrav) * alnpq * alnpq + &
                rdry * ts0 / egrav * alnpq)
            end if
          end do
          zq(kxs+1) = ter(j,i)
          qdt(kxs+1) = d_zero
          do l = kxs , 1 , -1
            lp = min(l,kxs)
            lm = max(l-1,1)
            zq(l) = zq(l+1) - rdry / egrav * tv(j,i,l) * &
              log((sigma(l) *   ps(j,i) + pt ) / &
                  (sigma(l+1) * ps(j,i) + pt ))
            qdt(l) = qdt(l+1) + &
                    (u(j+1,i+1,l) * psdot(j+1,i+1) + &
                     u(j+1,i  ,l) * psdot(j+1,i  ) - &
                     u(j  ,i+1,l) * psdot(j  ,i+1) - &
                     u(j  ,i  ,l) * psdot(j  ,i  ) + &
                     v(j+1,i+1,l) * psdot(j+1,i+1) + &
                     v(j  ,i+1,l) * psdot(j  ,i+1) - &
                     v(j+1,i  ,l) * psdot(j+1,i  ) - &
                     v(j  ,i  ,l) * psdot(j  ,i  )) / &
                 dx2 * xmsfx(j,i) * xmsfx(j,i)  * dsigma(l) / ps(j,i)
            ubar = 0.125 * (u(j  ,i  ,lp) + u(j  ,i  ,lm) + &
                            u(j  ,i+1,lp) + u(j  ,i+1,lm) + &
                            u(j+1,i  ,lp) + u(j+1,i  ,lm) + &
                            u(j+1,i+1,lp) + u(j+1,i+1,lm))
            vbar = 0.125 * (v(j  ,i  ,lp) + v(j  ,i  ,lm) + &
                            v(j  ,i+1,lp) + v(j  ,i+1,lm) + &
                            v(j+1,i  ,lp) + v(j+1,i  ,lm) + &
                            v(j+1,i+1,lp) + v(j+1,i+1,lm))
            !  Calculate omega (msfx not inverted).
            omega(l) = ps(j,i) * qdt(l) + sigma(l) * &
                       ((ps(j+1,i  ) - ps(j-1,i  )) * ubar +  &
                        (ps(j  ,i+1) - ps(j  ,i-1)) * vbar) / dx2 * xmsfx(j,i)
          end do
          !
          !  Vertical velocity from interpolated omega, zero at top.
          !
          do k = 2 , kxs + 1
            km = max(k-1,1)
            if ( k == kxs+1 ) km = kxs - 1
            do ll = 1 , kxs
              l = ll
              if ( zq(l+1) < z0q(k) )  exit
            end do
            zu = zq(l)
            zl = zq(l+1)
            omegau = omega(l)
            omegal = omega(l+1)
            omegan(k) = (omegau * (z0q(k) - zl ) + &
              omegal * (zu - z0q(k))) / (zu - zl)
            !  W =~ -OMEGA/RHO0/G *1000*PS0/1000. (OMEGA IN CB)
            rho = (rho0(j,i,km) * (a(k) - sigma(k)) + &
                   rho0(j,i,k) * (sigma(k) - a(km) )) / (a(k) - a(km))
            wtmp(j,i,k) = - omegan(k) / rho / egrav
          end do
          wtmp(j,i,1) = d_zero
        end do
      end do
      do i = i1 , i2
        do j = j1 , j2
          wtop(j,i) = wtmp(j,i,1)
          do k = 2 , kxs + 1
            w(j,i,k-1) = wtmp(j,i,k)
          end do
        end do
      end do
    end subroutine nhw

end module mod_nhinterp

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

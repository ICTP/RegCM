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

  public :: nhsetup , nhbase , nhinterp , nhpp , nhw

  real(rk8) :: ptop = 50.0D0  ! Centibars
  real(rk8) :: piso           ! Pascal
  real(rk8) :: ptoppa         ! Pascal
  real(rk8) :: p0 = stdp      ! Pascal
  real(rk8) :: ts0 = stdt     ! Kelvin
  real(rk8) :: tlp = 50.0D0

  interface nhinterp
    module procedure nhinterp3d
    module procedure nhinterp4d
  end interface nhinterp

  contains

    subroutine nhsetup(ptp,p,ts,lp)
      implicit none
      real(rk8) , intent(in) :: ptp , p , ts , lp
      ptop = ptp
      p0 = p
      ts0 = ts
      tlp = lp
      ptoppa = ptop * d_1000
      piso = p0 * exp((tiso - ts0)/tlp)
    end subroutine nhsetup
    !
    ! Compute the nonhydrostatic base state.
    !
    subroutine nhbase(i1,i2,j1,j2,kx,a,ter,ps0,pr0,t0,rho0)
      implicit none
      integer(ik4) , intent(in) :: i1 , i2 , j1 , j2 , kx
      real(rk8) , pointer , intent(in) , dimension(:) :: a         ! Adim 0-1
      real(rk8) , pointer , intent(in) , dimension(:,:) :: ter     ! Meters
      real(rk8) , pointer , intent(out) , dimension(:,:) :: ps0    ! Pascal
      real(rk8) , pointer , intent(out) , dimension(:,:,:) :: pr0  ! Pascal
      real(rk8) , pointer , intent(out) , dimension(:,:,:) :: t0   ! Kelvin
      real(rk8) , pointer , intent(out) , dimension(:,:,:) :: rho0 ! kg/kg
      integer(ik4) :: i , j , k
      real(rk8) :: ac , alnp , b
      !
      ! Define ps0 from terrain heights and t0 profile.
      !
      do i = i1 , i2
        do j = j1 , j2
          ac = d_half * govr * ter(j,i) / tlp
          b = ts0 / tlp
          alnp = -b + sqrt(b*b - d_four * ac)
          ps0(j,i) = p0 * exp(alnp) - ptoppa
          ! Define reference state temperature at model points.
          do k = 1 , kx
            pr0(j,i,k) = ps0(j,i) * a(k) + ptoppa
            if ( pr0(j,i,k) > 5474.9 ) then
              t0(j,i,k) = max(ts0 + tlp*log(pr0(j,i,k)/p0), tiso)
            else
              if ( pr0(j,i,k) > 5474.9 ) then
                t0(j,i,k) = tiso + min(log(5474.9/pr0(j,i,k)),228.65D0)
              else
                t0(j,i,k) = 228.65D0 + &
                        2.0D0*min(log(5474.9/pr0(j,i,k)),271.15D0)
              end if
            end if
            rho0(j,i,k) = pr0(j,i,k) / rgas / t0(j,i,k)
          end do
        end do
      end do
    end subroutine nhbase
    !
    ! Interpolate the hydrostatic input to nonhydrostatic coordinate.
    !
    subroutine nhinterp3d(i1,i2,j1,j2,kxs,a,sigma,ter,f,tv,ps,ps0,intmeth)
      implicit none
      integer(ik4) , intent(in) :: i1 , i2 , j1 , j2 , kxs , intmeth
      real(rk8) , pointer , intent(in) , dimension(:) :: a
      real(rk8) , pointer , intent(in) , dimension(:,:) :: ter  ! Meters
      real(rk8) , pointer , intent(in) , dimension(:) :: sigma
      real(rk8) , pointer , intent(in) , dimension(:,:,:) :: tv
      real(rk8) , pointer , intent(in) , dimension(:,:) :: ps
      real(rk8) , pointer , intent(in) , dimension(:,:) :: ps0
      real(rk8) , pointer , intent(inout) , dimension(:,:,:) :: f
      integer(ik4) :: i , j , k , l , ll , ip , im , jp , jm
      real(rk8) :: fl , fu , pr0 , alnp , alnqvn
      real(rk8) :: ziso , zl , zu
      real(rk8) , dimension(1:kxs) :: fn
      real(rk8) , dimension(j1:j2,i1:i2,1:kxs) :: z , z0
      real(rk8) , dimension(1:kxs+1) :: zq
      !
      ! We expect ps and ps0 to be already interpolated on dot points
      !
      do k = 1 , kxs
        do i = i1 , i2
          do j = j1 , j2
            pr0 = ps0(j,i) * a(k) + ptoppa
            if ( pr0 < piso ) then
              alnp = log(piso/p0)
              ziso = -(d_half*rovg*tlp*alnp*alnp + rovg*ts0*alnp)
              z0(j,i,k) = ziso - rovg * tiso * log(pr0/piso)
            else
              alnp = log(pr0/p0)
              z0(j,i,k) = -(d_half*rovg*tlp*alnp*alnp + rovg*ts0*alnp)
            end if
          end do
        end do
      end do
      !
      !  Calculate heights of input temperature sounding for interpolation
      !  to nonhydrostatic model levels.
      !
      do i = i1 , i2
        do j = j1 , j2
          ip = min(i2-1,i)
          jp = min(j2-1,j)
          im = max(i1,i-1)
          jm = max(j1,j-1)
          zq(kxs+1) = ter(j,i)
          do k = kxs , 1 , -1
            zq(k) = zq(k+1) - rovg * tv(j,i,k) * &
              log((sigma(k)*ps(j,i)+ptop)/(sigma(k+1)*ps(j,i)+ptop))
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
              if (z(j,i,l+1) < z0(j,i,k)) exit
            end do
            zu = z(j,i,l)
            zl = z(j,i,l+1)
            if ( intmeth == 1 ) then
              fu = f(j,i,l)
              fl = f(j,i,l+1)
              fn(k) = (fu * (z0(j,i,k) - zl ) + &
                       fl * (zu - z0(j,i,k))) / (zu - zl)
            else if ( intmeth == 2 ) then
              f(j,i,l) = max(f(j,i,l), minqx)
              f(j,i,l+1) = max(f(j,i,l+1), minqx)
              if ( z0(j,i,k) > zu ) then
                fn(k) = f(j,i,l)
              else
                fu = log(f(j,i,l  ))
                fl = log(f(j,i,l+1))
                alnqvn = (fu * (z0(j,i,k) - zl ) + &
                          fl * (zu - z0(j,i,k))) / (zu - zl)
                fn(k) = exp(alnqvn)
              end if
            end if
          end do
          do k = 1 , kxs
            f(j,i,k) = fn(k)
          end do
        end do
      end do
    end subroutine nhinterp3d

    subroutine nhinterp4d(i1,i2,j1,j2,kxs,nn,a,sigma,geo,f,tv,ps,ps0)
      implicit none
      integer(ik4) , intent(in) :: i1 , i2 , j1 , j2 , kxs , nn
      real(rk8) , pointer , intent(in) , dimension(:) :: a
      real(rk8) , pointer , intent(in) , dimension(:,:) :: geo  ! Geopotential
      real(rk8) , pointer , intent(in) , dimension(:) :: sigma
      real(rk8) , pointer , intent(in) , dimension(:,:,:) :: tv
      real(rk8) , pointer , intent(in) , dimension(:,:) :: ps
      real(rk8) , pointer , intent(in) , dimension(:,:) :: ps0
      real(rk8) , pointer , intent(inout) , dimension(:,:,:,:) :: f
      integer(ik4) :: i , j , k , n , l , ll , ip , im , jp , jm
      real(rk8) :: fl , fu , pr0 , alnp , alnqvn
      real(rk8) :: ziso , zl , zu
      real(rk8) , dimension(1:kxs) :: fn
      real(rk8) , dimension(j1:j2,i1:i2,1:kxs) :: z , z0
      real(rk8) , dimension(1:kxs+1) :: zq
      !
      ! We expect ps and ps0 to be already interpolated on dot points
      !
      do k = 1 , kxs
        do i = i1 , i2
          do j = j1 , j2
            pr0 = ps0(j,i) * a(k) + ptoppa
            if ( pr0 < piso ) then
              alnp = log(piso/p0)
              ziso = -(d_half*rovg*tlp*alnp*alnp + rovg*ts0*alnp)
              z0(j,i,k) = ziso - rovg * tiso * log(pr0/piso)
            else
              alnp = log(pr0/p0)
              z0(j,i,k) = -(d_half*rovg*tlp*alnp*alnp + rovg*ts0*alnp)
            end if
          end do
        end do
      end do
      !
      !  Calculate heights of input temperature sounding for interpolation
      !  to nonhydrostatic model levels.
      !
      do i = i1 , i2
        do j = j1 , j2
          ip = min(i2-1,i)
          jp = min(j2-1,j)
          im = max(i1,i-1)
          jm = max(j1,j-1)
          zq(kxs+1) = geo(j,i)*regrav
          do k = kxs , 1 , -1
            zq(k) = zq(k+1) - rovg * tv(j,i,k) * &
              log((sigma(k)*ps(j,i)+ptop)/(sigma(k+1)*ps(j,i)+ptop))
            z(j,i,k) = d_half * (zq(k) + zq(k+1))
          end do
        end do
      end do
      !
      ! Interpolate from z to z0 levels.
      !
      do n = 1 , nn
        do i = i1 , i2
          do j = j1 , j2
            do k = 1 , kxs
              do ll = 1 , kxs - 1
                l = ll
                if (z(j,i,l+1) < z0(j,i,k)) exit
              end do
              zu = z(j,i,l)
              zl = z(j,i,l+1)
              f(j,i,l,  n) = max(f(j,i,l,  n), dlowval)
              f(j,i,l+1,n) = max(f(j,i,l+1,n), dlowval)
              if ( z0(j,i,k) > zu ) then
                fn(k) = f(j,i,l,n)
              else
                fu = log(f(j,i,l  ,n))
                fl = log(f(j,i,l+1,n))
                alnqvn = (fu * (z0(j,i,k) - zl ) + &
                          fl * (zu - z0(j,i,k))) / (zu - zl)
                fn(k) = max(exp(alnqvn),dlowval)
              end if
            end do
            do k = 1 , kxs
              f(j,i,k,n) = fn(k)
            end do
          end do
        end do
      end do
    end subroutine nhinterp4d
    !
    ! Compute the pressure perturbation pp.
    !
    subroutine nhpp(i1,i2,j1,j2,kxs,sigma,t,pr0,t0,tv,ps,ps0,pp)
      implicit none
      integer(ik4) , intent(in) :: i1 , i2 , j1 , j2 , kxs
      integer(ik4) :: i , j , k
      real(rk8) , pointer , intent(in) , dimension(:) :: sigma
      real(rk8) , pointer , intent(in) , dimension(:,:,:) :: t , t0 , tv
      real(rk8) , pointer , intent(in) , dimension(:,:,:) :: pr0
      real(rk8) , pointer , intent(in) , dimension(:,:) :: ps , ps0
      real(rk8) , pointer , intent(out) , dimension(:,:,:) :: pp
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
          p0surf = ps0(j,i) + ptoppa
          psp = ps(j,i) * d_1000 - ps0(j,i)
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
            cc = egrav * wtu / pr0(j,i,k  ) * t0(j,i,k  ) / tk
            tvpot = wtl * ((tvkp1 - t0(j,i,k+1)) / tkp1) + &
                    wtu * ((tvk   - t0(j,i,k  )) / tk  )
            pp(j,i,k) = (egrav * tvpot + pp(j,i,k+1) * (aa - bb)) / (aa + cc)
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
            tvpot = wtl * (tv(j,i,k+1) - t0(j,i,k+1)) / t(j,i,k+1) + &
                    wtu * (tv(j,i,k  ) - t0(j,i,k  )) / t(j,i,k  )
            checkl = (pp(j,i,k+1) - pp(j,i,k)) / (pr0(j,i,k+1) - pr0(j,i,k))
            checkr = tvpot - &
                 wtl * t0(j,i,k+1)/t(j,i,k+1) * pp(j,i,k+1)/pr0(j,i,k+1) - &
                 wtu * t0(j,i,k)  /t(j,i,k)   * pp(j,i,k)  /pr0(j,i,k)
            check = checkl + checkr
            if ( abs(check) > 1.D-2 ) then
              write(stderr,'(A,3I4,3g14.6)') &
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
    subroutine nhw(i1,i2,j1,j2,kxs,a,sigma,dsigma,ter,u,v,tv,rho0, &
                   ps,psdot,ps0,xmsfx,xmsfd,w,wtop,ds)
      implicit none
      integer(ik4) , intent(in) :: i1 , i2 , j1 , j2 , kxs
      real(rk8) , pointer , intent(in) , dimension(:) :: a , sigma , dsigma
      real(rk8) , pointer , intent(in) , dimension(:,:) :: ter  ! Meters
      real(rk8) , pointer , intent(in) , dimension(:,:) :: xmsfx , xmsfd
      real(rk8) , pointer , intent(in) , dimension(:,:) :: ps , ps0 , psdot
      real(rk8) , pointer , intent(in) , dimension(:,:,:) :: rho0 , tv
      real(rk8) , pointer , intent(in) , dimension(:,:,:) :: u , v
      real(rk8) , intent(in) :: ds                    ! Kilometers
      real(rk8) , pointer , intent(out) , dimension(:,:,:) :: w
      real(rk8) , pointer , intent(out) , dimension(:,:) :: wtop
      integer(ik4) :: i , j , k , km , kp
      integer(ik4) :: l , ll , lm , lp , ip , im , jp , jm
      real(rk8) :: alnpq , dx2 , omegal , omegau , ubar , vbar
      real(rk8) :: pr0 , ziso , zl , zu , rho , omegan
      real(rk8) , dimension(kxs+1) :: omega , qdt
      real(rk8) , dimension(kxs+1) :: z0q , zq
      real(rk8) , dimension(j1:j2,i1:i2,kxs+1) :: wtmp
      real(rk8) , dimension(j1:j2,i1:i2) :: pspa
      real(rk8) , dimension(j1:j2,i1:i2) :: psdotpa

      wtmp(:,:,:) = d_zero
      dx2 = d_two * ds
      psdotpa = psdot * d_1000
      pspa = ps * d_1000
      do i = i1 , i2
        do j = j1 , j2
          ip = min(i+1,i2)
          im = max(i-1,i1)
          jp = min(j+1,j2)
          jm = max(j-1,j1)
          if ( ps0(j,i) + ptoppa < piso ) THEN
            write(stderr,'(a,f5.1,a,f6.1,a)') &
              'The chosen value of Tiso, ',tiso,' K occurs at ',piso,' hPa.'
            write(stderr,'(A)') &
              'This value of pressure is greater than ref. surface pressure.'
          end if
          do k = 1 , kxs
            pr0 = ps0(j,i) * sigma(k) + ptoppa
            if ( pr0 < piso ) then
              alnpq = log(piso/p0)
              ziso = -(d_half*rovg*tlp*alnpq*alnpq + rovg*ts0*alnpq)
              z0q(k) = ziso - rovg*tiso*log(pr0/piso)
            else
              alnpq = log(pr0/p0)
              z0q(k) = -(d_half*rovg*tlp*alnpq*alnpq + rovg*ts0*alnpq)
            end if
          end do
          z0q(kxs+1) = ter(j,i)
          zq(kxs+1) = ter(j,i)
          do l = kxs , 1 , -1
            zq(l) = zq(l+1) - rovg*tv(j,i,l) * &
              log((sigma(l)   * pspa(j,i) + ptoppa) / &
                  (sigma(l+1) * pspa(j,i) + ptoppa))
          end do
          qdt(1) = d_zero
          do l = 2 , kxs
            qdt(l) = qdt(l-1) - &
                      (u(jp,ip,l) * psdotpa(jp,ip) / xmsfd(jp,ip) +  &
                       u(jp,i ,l) * psdotpa(jp,i ) / xmsfd(jp,i ) -  &
                       u(j ,ip,l) * psdotpa(j ,ip) / xmsfd(j ,ip) -  &
                       u(j ,i ,l) * psdotpa(j ,i ) / xmsfd(j ,i ) +  &
                       v(jp,ip,l) * psdotpa(jp,ip) / xmsfd(jp,ip) +  &
                       v(j ,ip,l) * psdotpa(j ,ip) / xmsfd(j ,ip) -  &
                       v(jp,i ,l) * psdotpa(jp,i ) / xmsfd(jp,i ) -  &
                       v(j ,i ,l) * psdotpa(j ,i ) / xmsfd(j ,i )) / &
                   dx2 * xmsfx(j,i) * xmsfx(j,i) * dsigma(l) / pspa(j,i)
          end do
          qdt(kxs+1) = d_zero
          omega(1) = d_zero
          do l = 2 , kxs+1
            lp = min(l,kxs)
            lm = max(l-1,1)
            ubar = 0.125D0 * (u(j ,i ,lp) + u(j ,i ,lm) + &
                              u(j ,ip,lp) + u(j ,ip,lm) + &
                              u(jp,i ,lp) + u(jp,i ,lm) + &
                              u(jp,ip,lp) + u(jp,ip,lm))
            vbar = 0.125D0 * (v(j ,i ,lp) + v(j ,i ,lm) + &
                              v(j ,ip,lp) + v(j ,ip,lm) + &
                              v(jp,i ,lp) + v(jp,i ,lm) + &
                              v(jp,ip,lp) + v(jp,ip,lm))
            ! Calculate omega
            omega(l) = pspa(j,i) * qdt(l) + sigma(l) * &
                       ((pspa(jp,i ) - pspa(jm,i )) * ubar +  &
                        (pspa(j ,ip) - pspa(j ,im)) * vbar) / dx2 * xmsfx(j,i)
          end do
          !
          !  Vertical velocity from interpolated omega, zero at top.
          !
          do k = 2 , kxs + 1
            kp = min(k,kxs)
            km = min(max(k-1,1),kxs-1)
            do ll = 2 , kxs + 1
              l = ll
              if ( zq(l) < z0q(k) ) exit
            end do
            zu = zq(l-1)
            zl = zq(l)
            omegau = omega(l-1)
            omegal = omega(l)
            omegan = (omegau * (z0q(k) - zl ) + &
                      omegal * (zu - z0q(k))) / (zu - zl)
            !  W =~ -OMEGA/RHO0/G *1000*PS0/1000. (OMEGA IN CB)
            rho = (rho0(j,i,km) * (a(kp) - sigma(k)) + &
                   rho0(j,i,kp) * (sigma(k) - a(km))) / (a(kp)-a(km))
            wtmp(j,i,k) = -omegan/rho * regrav
          end do
        end do
      end do
      wtop(j1:j2,i1:i2) = wtmp(j1:j2,i1:i2,1)
      do k = 2 , kxs + 1
        do i = i1 , i2
          do j = j1 , j2
            w(j,i,k-1) = wtmp(j,i,k)
          end do
        end do
      end do
    end subroutine nhw

end module mod_nhinterp

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

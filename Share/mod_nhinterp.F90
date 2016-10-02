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

  public :: nhsetup , temppres , nhbase , nhinterp , nhpp , nhw

  real(rkx) :: ptop = 50.0_rkx  ! Centibars
  real(rkx) :: piso             ! Pascal
  real(rkx) :: ptoppa           ! Pascal
  real(rkx) :: p0 = stdp        ! Pascal
  real(rkx) :: ts0 = stdt       ! Kelvin
  real(rkx) :: tlp = 50.0_rkx

  interface nhinterp
    module procedure nhinterp3d
    module procedure nhinterp4d
  end interface nhinterp

  contains

    subroutine nhsetup(ptp,p,ts,lp)
      implicit none
      real(rkx) , intent(in) :: ptp , p , ts , lp
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
    real(rkx) function temppres(pr) result(tp)
      implicit none
      real(rkx) , intent(in) :: pr
      if ( pr > 22632.0_rkx ) then
        tp = max(ts0 + tlp*log(pr/p0), tiso)
      else if ( pr > 5474.9_rkx ) then
        tp = tiso + min(log(5474.9_rkx/pr),228.65_rkx)
      else
        tp = 228.65_rkx + 2.0_rkx*min(log(5474.9_rkx/pr),271.15_rkx)
      end if
    end function temppres

    subroutine nhbase(i1,i2,j1,j2,kx,sigmah,ter,ps0,pr0,t0,rho0)
      implicit none
      integer(ik4) , intent(in) :: i1 , i2 , j1 , j2 , kx
      real(rkx) , pointer , intent(in) , dimension(:) :: sigmah    ! Adim 0-1
      real(rkx) , pointer , intent(in) , dimension(:,:) :: ter     ! Meters
      real(rkx) , pointer , intent(out) , dimension(:,:) :: ps0    ! Pascal
      real(rkx) , pointer , intent(out) , dimension(:,:,:) :: pr0  ! Pascal
      real(rkx) , pointer , intent(out) , dimension(:,:,:) :: t0   ! Kelvin
      real(rkx) , pointer , intent(out) , dimension(:,:,:) :: rho0 ! kg/kg
      integer(ik4) :: i , j , k
      real(rkx) :: ac , alnp , b
      !
      ! Define ps0 from terrain heights and t0 profile.
      !
      do i = i1 , i2
        do j = j1 , j2
          ac = d_half * govr * ter(j,i) / tlp
          b = ts0 / tlp
          alnp = -b + sqrt(b*b - d_four * ac)
          ps0(j,i) = p0 * exp(alnp) - ptoppa
          if ( ps0(j,i) + ptoppa < piso ) THEN
            write(stderr,'(a,f5.1,a,f6.1,a)') &
              'The chosen value of Tiso, ',tiso,' K occurs at ',piso,' hPa.'
            write(stderr,'(a,i4,a,i4)') &
              'The value of pressure in point (J,I) ',j,',',i
            write(stderr,'(a)') 'is lower than tropopause pressure.'
          end if
          ! Define reference state temperature at model points.
          do k = 1 , kx
            pr0(j,i,k) = ps0(j,i) * sigmah(k) + ptoppa
            t0(j,i,k) = temppres(pr0(j,i,k))
            rho0(j,i,k) = pr0(j,i,k) / rgas / t0(j,i,k)
          end do
        end do
      end do
    end subroutine nhbase
    !
    ! Interpolate the hydrostatic input to nonhydrostatic coordinate.
    !
    subroutine nhinterp3d(i1,i2,j1,j2,kxs,sigmah,sigma,ter,f,tv,ps,ps0,intmeth)
      implicit none
      integer(ik4) , intent(in) :: i1 , i2 , j1 , j2 , kxs , intmeth
      real(rkx) , pointer , intent(in) , dimension(:) :: sigmah
      real(rkx) , pointer , intent(in) , dimension(:,:) :: ter  ! Meters
      real(rkx) , pointer , intent(in) , dimension(:) :: sigma
      real(rkx) , pointer , intent(in) , dimension(:,:,:) :: tv
      real(rkx) , pointer , intent(in) , dimension(:,:) :: ps
      real(rkx) , pointer , intent(in) , dimension(:,:) :: ps0
      real(rkx) , pointer , intent(inout) , dimension(:,:,:) :: f
      integer(ik4) :: i , j , k , l , ll
      real(rkx) :: fl , fu , pr0 , alnp , alnqvn
      real(rkx) :: ziso , zl , zu , wu , wl
      real(rkx) , dimension(1:kxs) :: fn
      real(rkx) , dimension(j1:j2,i1:i2,1:kxs) :: z , z0
      real(rkx) , dimension(1:kxs+1) :: zq
      !
      ! We expect ps and ps0 to be already interpolated on dot points
      !
      do k = 1 , kxs
        do i = i1 , i2
          do j = j1 , j2
            pr0 = ps0(j,i) * sigmah(k) + ptoppa
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
              f(j,i,l) = max(f(j,i,l), minqq)
              f(j,i,l+1) = max(f(j,i,l+1), minqq)
              if ( z0(j,i,k) > zu ) then
                fn(k) = f(j,i,l)
              else
                fu = log(f(j,i,l  ))
                fl = log(f(j,i,l+1))
                wu = (z0(j,i,k) - zl) / (zu - zl)
                wl = d_one - wu
                alnqvn = fu * wu + fl * wl
                fn(k) = exp(alnqvn)
              end if
              if ( fn(k) < minqq ) fn(k) = minqq
            end if
          end do
          do k = 1 , kxs
            f(j,i,k) = fn(k)
          end do
        end do
      end do
    end subroutine nhinterp3d

    subroutine nhinterp4d(i1,i2,j1,j2,kxs,nn,sigmah,sigma,geo,f,tv,ps,ps0)
      implicit none
      integer(ik4) , intent(in) :: i1 , i2 , j1 , j2 , kxs , nn
      real(rkx) , pointer , intent(in) , dimension(:) :: sigmah
      real(rkx) , pointer , intent(in) , dimension(:,:) :: geo  ! Geopotential
      real(rkx) , pointer , intent(in) , dimension(:) :: sigma
      real(rkx) , pointer , intent(in) , dimension(:,:,:) :: tv
      real(rkx) , pointer , intent(in) , dimension(:,:) :: ps
      real(rkx) , pointer , intent(in) , dimension(:,:) :: ps0
      real(rkx) , pointer , intent(inout) , dimension(:,:,:,:) :: f
      integer(ik4) :: i , j , k , n , l , ll
      real(rkx) :: fl , fu , pr0 , alnp , alnqvn
      real(rkx) :: ziso , zl , zu , wl , wu
      real(rkx) , dimension(1:kxs) :: fn
      real(rkx) , dimension(j1:j2,i1:i2,1:kxs) :: z , z0
      real(rkx) , dimension(1:kxs+1) :: zq
      !
      ! We expect ps and ps0 to be already interpolated on dot points
      !
      do k = 1 , kxs
        do i = i1 , i2
          do j = j1 , j2
            pr0 = ps0(j,i) * sigmah(k) + ptoppa
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
              f(j,i,l,  n) = max(f(j,i,l,  n),mintr)
              f(j,i,l+1,n) = max(f(j,i,l+1,n),mintr)
              if ( z0(j,i,k) > zu ) then
                fn(k) = f(j,i,l,n)
              else
                fu = log(f(j,i,l  ,n))
                fl = log(f(j,i,l+1,n))
                wu = (z0(j,i,k) - zl) / (zu - zl)
                wl = d_one - wu
                alnqvn = fu * wu + fl * wl
                fn(k) = exp(alnqvn)
              end if
              if ( fn(k) < dlowval ) fn(k) = mintr
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
      real(rkx) , pointer , intent(in) , dimension(:) :: sigma
      real(rkx) , pointer , intent(in) , dimension(:,:,:) :: t , t0 , tv
      real(rkx) , pointer , intent(in) , dimension(:,:,:) :: pr0
      real(rkx) , pointer , intent(in) , dimension(:,:) :: ps , ps0
      real(rkx) , pointer , intent(out) , dimension(:,:,:) :: pp
      real(rkx) :: aa , bb , cc , check , checkl , checkr , delp0 , p0surf
      real(rkx) :: psp , tk , tkp1 , tvk , tvkp1 , tvpot , wtl , wtu
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
            wtu = d_one - wtl
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
            wtu = (sigma(k+2) - sigma(k+1)) / (sigma(k+2) - sigma(k  ))
            wtl = d_one - wtu
            tvpot = wtl * (tv(j,i,k+1) - t0(j,i,k+1)) / t(j,i,k+1) + &
                    wtu * (tv(j,i,k  ) - t0(j,i,k  )) / t(j,i,k  )
            checkl = (pp(j,i,k+1) - pp(j,i,k)) / (pr0(j,i,k+1) - pr0(j,i,k))
            checkr = tvpot - &
                 wtl * t0(j,i,k+1)/t(j,i,k+1) * pp(j,i,k+1)/pr0(j,i,k+1) - &
                 wtu * t0(j,i,k)  /t(j,i,k)   * pp(j,i,k)  /pr0(j,i,k)
            check = checkl + checkr
            if ( abs(check) > 1.e-2_rkx ) then
              write(stderr,'(A,3I4,3g14.6)') &
                'NHPP vert gradient check error at ',i,j,k,check,checkl,checkr
            end if
          end do
        end do
      end do
    end subroutine nhpp
    !
    ! Compute the nonhydrostatic initial vertical velocity (w) from the
    ! horizontal wind fields (u and v).
    !
    subroutine nhw(i1,i2,j1,j2,kxs,sigmah,sigma,dsigma,ter,u,v,tv,rho0, &
                   ps,psdot,ps0,xmsfx,xmsfd,w,wtop,ds,iband)
      implicit none
      integer(ik4) , intent(in) :: i1 , i2 , j1 , j2 , kxs , iband
      real(rkx) , pointer , intent(in) , dimension(:) :: sigmah , sigma , dsigma
      real(rkx) , pointer , intent(in) , dimension(:,:) :: ter  ! Meters
      real(rkx) , pointer , intent(in) , dimension(:,:) :: xmsfx , xmsfd
      real(rkx) , pointer , intent(in) , dimension(:,:) :: ps , ps0 , psdot
      real(rkx) , pointer , intent(in) , dimension(:,:,:) :: rho0 , tv
      real(rkx) , pointer , intent(in) , dimension(:,:,:) :: u , v
      real(rkx) , intent(in) :: ds                    ! Kilometers
      real(rkx) , pointer , intent(out) , dimension(:,:,:) :: w
      real(rkx) , pointer , intent(out) , dimension(:,:) :: wtop
      integer(ik4) :: i , j , k , km , kp
      integer(ik4) :: l , ll , lm , lp , ip , im , jp , jm
      integer(ik4) :: ipp , imm , jpp , jmm
      real(rkx) :: alnpq , dx2 , omegal , omegau , ubar , vbar , wu , wl
      real(rkx) :: ziso , zl , zu , rho , omegan , pten
      real(rkx) :: ua , ub , va , vb
      real(rkx) , dimension(kxs) :: mdv
      real(rkx) , dimension(kxs+1) :: omega , qdt
      real(rkx) , dimension(kxs+1) :: z0q , zq
      real(rkx) , dimension(j1:j2,i1:i2,kxs+1) :: wtmp
      real(rkx) , dimension(j1:j2,i1:i2,kxs) :: pr0 , logprp
      real(rkx) , dimension(j1:j2,i1:i2) :: pspa , rpspa , smt1 , smt2
      real(rkx) , dimension(j1:j2,i1:i2) :: dummy , dummy1
      real(rkx) , dimension(j1:j2,i1:i2) :: psdotpam

      wtmp(:,:,:) = d_zero
      dx2 = d_two * ds
      psdotpam = psdot * d_1000 / xmsfd
      pspa = ps * d_1000
      rpspa = d_one / pspa
      dummy = (xmsfx * xmsfx) / dx2
      dummy1 = xmsfx / dx2

      do k = 1 , kxs
        do i = i1 , i2
          do j = j1 , j2
            pr0(j,i,k) = ps0(j,i) * sigma(k) + ptoppa
            logprp(j,i,k) = log((sigma(k)   * pspa(j,i) + ptoppa) / &
                                (sigma(k+1) * pspa(j,i) + ptoppa))
          end do
        end do
      end do

      do i = i1 , i2
        do j = j1 , j2
          z0q(kxs+1) = ter(j,i)
          do k = 1 , kxs
            if ( pr0(j,i,k) < piso ) then
              alnpq = log(piso/p0)
              ziso = -(d_half*rovg*tlp*alnpq*alnpq + rovg*ts0*alnpq)
              z0q(k) = ziso - rovg*tiso*log(pr0(j,i,k)/piso)
            else
              alnpq = log(pr0(j,i,k)/p0)
              z0q(k) = -(d_half*rovg*tlp*alnpq*alnpq + rovg*ts0*alnpq)
            end if
          end do
          zq(kxs+1) = ter(j,i)
          do l = kxs , 1 , -1
            zq(l) = zq(l+1) - rovg*tv(j,i,l) * logprp(j,i,l)
          end do
          ip = min(i+1,i2)
          im = max(i-1,i1)
          ipp = min(i+2,i2)
          imm = max(i-2,i1)
          if ( iband /= 1 ) then
            jp = min(j+1,j2)
            jm = max(j-1,j1)
            jpp = min(j+2,j2)
            jmm = max(j-2,j1)
          else
            if ( j == j2-1 ) then
              jpp = j1
              jp = j2
            else if ( j == j2 ) then
              jpp = j1 + 1
              jp = j1
            else
              jp = j + 1
              jpp = j + 2
            end if
            if ( j == j1 + 1 ) then
              jmm = j2
              jm = j1
            else if ( j == j1 ) then
              jmm = j2 - 1
              jm = j2
            else
              jmm = j - 2
              jm = j - 1
            end if
          end if
          pten = d_zero
          mdv(:) = d_zero
          do l = 1 , kxs
            ua = d_rfour * ( u(j ,i ,l) * psdotpam(j ,i ) + &
                             u(jm,i ,l) * psdotpam(jm,i ) + &
                             u(jm,ip,l) * psdotpam(jm,ip) + &
                             u(j ,ip,l) * psdotpam(j ,ip) )
            ub = d_rfour * ( u(jp, i ,l) * psdotpam(jp ,i ) + &
                             u(jpp,i ,l) * psdotpam(jpp,i ) + &
                             u(jpp,ip,l) * psdotpam(jpp,ip) + &
                             u(jp ,ip,l) * psdotpam(jp, ip) )
            va = d_rfour * ( v(j ,i ,l) * psdotpam(j ,i ) + &
                             v(jp,i ,l) * psdotpam(jp,i ) + &
                             v(jp,im,l) * psdotpam(jp,im) + &
                             v(j ,im,l) * psdotpam(j ,im) )
            vb = d_rfour * ( v(j ,ip ,l) * psdotpam(j ,ip ) + &
                             v(jp,ip ,l) * psdotpam(jp,ip ) + &
                             v(jp,ipp,l) * psdotpam(jp,ipp) + &
                             v(j ,ipp,l) * psdotpam(j ,ipp) )
            mdv(l) = (ub - ua + vb - va) * dummy(j,i)
            pten = pten - mdv(l) * dsigma(l)
          end do
          qdt(1) = d_zero
          do l = 2 , kxs + 1
            qdt(l) = qdt(l-1) - (pten + mdv(l-1)) * dsigma(l-1) * rpspa(j,i)
          end do
          do l = 2 , kxs+1
            lp = min(l,kxs)
            lm = max(l-1,1)
            ubar = 0.125_rkx * (u(j ,i ,lp) + u(j ,i ,lm)  + &
                                u(j ,ip,lp) + u(j ,ip,lm)  + &
                                u(jp,i ,lp) + u(jp,i ,lm)  + &
                                u(jp,ip,lp) + u(jp,ip,lm))
            vbar = 0.125_rkx * (v(j ,i ,lp) + v(j ,i ,lm)  + &
                                v(j ,ip,lp) + v(j ,ip,lm)  + &
                                v(jp,i ,lp) + v(jp,i ,lm)  + &
                                v(jp,ip,lp) + v(jp,ip,lm))
            ! Calculate omega
            omega(l) = pspa(j,i) * qdt(l) + sigma(l) *  &
                    ((pspa(jp,i) - pspa(jm,i)) * ubar + &
                     (pspa(j,ip) - pspa(j,im)) * vbar) * dummy1(j,i)
          end do
          omega(1) = d_zero
          !
          !  Vertical velocity from interpolated omega, zero at top.
          !
          do k = 2 , kxs + 1
            kp = min(k,kxs)
            km = min(max(k-1,1),kxs-1)
            do ll = 1 , kxs
              l = ll
              if ( zq(l+1) < z0q(k) ) exit
            end do
            zu = zq(l)
            zl = zq(l+1)
            omegau = omega(l)
            omegal = omega(l+1)
            wu = (z0q(k) - zl) / (zu - zl)
            wl = d_one - wu
            omegan = omegau * wu + omegal * wl
            wu = (sigmah(kp) - sigma(k)) / (sigmah(kp)-sigmah(km))
            wl = d_one - wu
            rho = rho0(j,i,km) * wu + rho0(j,i,kp) * wl
            !  W =~ -OMEGA/RHO0/G *1000*PS0/1000. (OMEGA IN CB)
            wtmp(j,i,k) = -omegan/rho * regrav
          end do
          wtmp(j,i,1) = -omega(2)/rho0(j,i,1)*regrav
        end do
      end do
      wtop(j1:j2,i1:i2) = wtmp(j1:j2,i1:i2,1)
      do k = 2 , kxs + 1
        smt1(:,:) = wtmp(:,:,k)
        smt2(:,:) = wtmp(:,:,k)
        do i = i1 , i2
          do j = j1+2 , j2-2
            smt2(j,i) = 0.10_rkx*(d_four*smt1(j,i)+ &
                                  d_two*smt1(j+1,i)+d_two*smt1(j-1,i) + &
                                  smt1(j+2,i)+smt1(j-2,i))
          end do
          smt2(2,i) = d_rfour*(d_two*smt1(2,i)+smt1(1,i)+smt1(3,i))
          smt2(j2-1,i) = d_rfour*(d_two*smt1(j2-1,i)+smt1(j2,i)+smt1(j2-2,i))
        end do
        do i = i1+2 , i2-2
          do j = j1 , j2
            smt1(j,i) = 0.10_rkx*(d_four*smt2(j,i)+ &
                                  d_two*smt2(j,i+1)+d_two*smt2(j,i-1) * &
                                  smt2(j,i+2)+smt2(j,i-2))
          end do
        end do
        do j = j1 , j2
          smt1(j,2) = d_rfour*(d_two*smt2(j,2)+smt2(j,3)+smt2(j,1))
          smt1(j,i2-1) = d_rfour*(d_two*smt2(j,i2-1)+smt2(j,i2-2)+smt2(j,i2))
        end do
        w(:,:,k-1) = smt1(:,:)
      end do
      w(j2,:,:) = d_zero
      w(:,i2,:) = d_zero
      wtop(j2,:) = d_zero
      wtop(:,i2) = d_zero
    end subroutine nhw

end module mod_nhinterp

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

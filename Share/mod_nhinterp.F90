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
  use mod_stdatm

  implicit none

  private

  public :: base_state_temperature , nhsetup , nhbase , nhinterp , nhpp , nhw

  real(rkx) :: ptop = 5.0_rkx      ! Centibars
  real(rkx) :: ptoppa = 5000.0_rkx ! Pascal
  real(rkx) :: p0 = 101325.0_rkx   ! Pascal
  real(rkx) :: tlp = 47.70_rkx     ! [K/ln(Pa)]
  real(rkx) :: st0 = 288.15        ! K

  interface nhinterp
    module procedure nhinterp3d
    module procedure nhinterp4d
  end interface nhinterp

  contains

    real(rkx) function base_state_temperature(i1,i2,j1,j2,xlat)
      implicit none
      integer(ik4) , intent(in) :: i1 , i2 , j1 , j2
      real(rkx) , pointer , intent(in) , dimension(:,:) :: xlat    ! Latitude
      real(rkx) , dimension(j1:j2,i1:i2) :: ts0
      integer(ik4) :: i , j
      do i = i1 , i2
        do j = j1 , j2
          ts0(j,i) = stdatm_val(xlat(j,i),p0*d_r100,istdatm_tempk)
        end do
      end do
      base_state_temperature = d_half*(maxval(ts0)+minval(ts0))
      write(stdout,*) 'Computed reference surface temperature for domain: ', &
         base_state_temperature
    end function base_state_temperature

    subroutine nhsetup(ptp,pbase,lp,ts0)
      implicit none
      real(rkx) , intent(in) :: ptp , pbase , lp , ts0
      ptop = ptp
      p0 = pbase
      tlp = lp
      ptoppa = ptop * d_1000
      st0 = ts0
    end subroutine nhsetup
    !
    ! Compute the nonhydrostatic base state.
    !
    subroutine nhbase(i1,i2,j1,j2,kx,sig,ter,ps0,pr0,t0,rho0,z0)
      implicit none
      integer(ik4) , intent(in) :: i1 , i2 , j1 , j2 , kx
      real(rkx) , pointer , intent(in) , dimension(:) :: sig       ! Adim 0-1
      real(rkx) , pointer , intent(in) , dimension(:,:) :: ter     ! Meters
      real(rkx) , pointer , intent(inout) , dimension(:,:) :: ps0    ! Pascal
      real(rkx) , pointer , intent(inout) , dimension(:,:,:) :: pr0  ! Pascal
      real(rkx) , pointer , intent(inout) , dimension(:,:,:) :: t0   ! Kelvin
      real(rkx) , pointer , intent(inout) , dimension(:,:,:) :: rho0 ! kg/kg
      real(rkx) , pointer , intent(inout) , dimension(:,:,:) :: z0   ! m
      integer(ik4) :: i , j , k
      real(rkx) :: ac , alnp , b
      !
      ! Define ps0 from terrain heights and t0 profile.
      !
      do i = i1 , i2
        do j = j1 , j2
          ac = d_half * govr * ter(j,i) / tlp
          b = st0 / tlp
          alnp = -b + sqrt(b*b - d_four * ac)
          ps0(j,i) = p0 * exp(alnp) - ptoppa
          ! Define reference state temperature at model points.
          do k = 1 , kx
            pr0(j,i,k) = ps0(j,i) * sig(k) + ptoppa
            t0(j,i,k) = max(st0 + tlp * log(pr0(j,i,k) / p0),tiso)
            rho0(j,i,k) = pr0(j,i,k) / rgas / t0(j,i,k)
            alnp = log(pr0(j,i,k)/(ps0(j,i)+ptoppa))
            z0(j,i,k) = max(-(d_half*rovg*tlp*alnp*alnp + rovg*st0*alnp),d_zero)
          end do
        end do
      end do
    end subroutine nhbase
    !
    ! Interpolate the hydrostatic input to nonhydrostatic coordinate.
    !
    subroutine nhinterp3d(i1,i2,j1,j2,kxs,sigmah,sigma,f,tv,ps,ps0,intmeth)
      implicit none
      integer(ik4) , intent(in) :: i1 , i2 , j1 , j2 , kxs , intmeth
      real(rkx) , pointer , intent(in) , dimension(:) :: sigmah
      real(rkx) , pointer , intent(in) , dimension(:) :: sigma
      real(rkx) , pointer , intent(in) , dimension(:,:,:) :: tv
      real(rkx) , pointer , intent(in) , dimension(:,:) :: ps
      real(rkx) , pointer , intent(in) , dimension(:,:) :: ps0
      real(rkx) , pointer , intent(inout) , dimension(:,:,:) :: f
      integer(ik4) :: i , j , k , l , ll
      real(rkx) :: fl , fu , pr0 , alnqvn , alnp
      real(rkx) :: zl , zu , wu , wl
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
            alnp = log(pr0/(ps0(j,i)+ptoppa))
            z0(j,i,k) = - (d_half*rovg*tlp*alnp*alnp + rovg*st0*alnp)
          end do
        end do
      end do
      !
      !  Calculate heights of input temperature sounding for interpolation
      !  to nonhydrostatic model levels.
      !
      do i = i1 , i2
        do j = j1 , j2
          zq(kxs+1) = d_zero
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
      if ( intmeth == 1 ) then
        do i = i1 , i2
          do j = j1 , j2
            do k = 1 , kxs
              do ll = 1 , kxs - 1
                l = ll
                if (z(j,i,l+1) < z0(j,i,k)) exit
              end do
              zu = z(j,i,l)
              zl = z(j,i,l+1)
              fu = f(j,i,l)
              fl = f(j,i,l+1)
              fn(k) = (fu * (z0(j,i,k) - zl ) + &
                       fl * (zu - z0(j,i,k))) / (zu - zl)
            end do
            do k = 1 , kxs
              f(j,i,k) = fn(k)
            end do
          end do
        end do
      else
        do i = i1 , i2
          do j = j1 , j2
            do k = 1 , kxs
              do ll = 1 , kxs - 1
                l = ll
                if (z(j,i,l+1) < z0(j,i,k)) exit
              end do
              zu = z(j,i,l)
              zl = z(j,i,l+1)
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
            end do
            do k = 1 , kxs
              f(j,i,k) = fn(k)
            end do
          end do
        end do
      end if
    end subroutine nhinterp3d

    subroutine nhinterp4d(i1,i2,j1,j2,kxs,nn,sigmah,sigma,f,tv,ps,ps0)
      implicit none
      integer(ik4) , intent(in) :: i1 , i2 , j1 , j2 , kxs , nn
      real(rkx) , pointer , intent(in) , dimension(:) :: sigmah
      real(rkx) , pointer , intent(in) , dimension(:) :: sigma
      real(rkx) , pointer , intent(in) , dimension(:,:,:) :: tv
      real(rkx) , pointer , intent(in) , dimension(:,:) :: ps
      real(rkx) , pointer , intent(in) , dimension(:,:) :: ps0
      real(rkx) , pointer , intent(inout) , dimension(:,:,:,:) :: f
      integer(ik4) :: i , j , k , n , l , ll
      real(rkx) :: fl , fu , pr0 , alnqvn , alnp
      real(rkx) :: zl , zu , wl , wu
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
            alnp = log(pr0/(ps0(j,i)+ptoppa))
            z0(j,i,k) = - (d_half*rovg*tlp*alnp*alnp + rovg*st0*alnp)
          end do
        end do
      end do
      !
      !  Calculate heights of input temperature sounding for interpolation
      !  to nonhydrostatic model levels.
      !
      do i = i1 , i2
        do j = j1 , j2
          zq(kxs+1) = d_zero
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
      real(rkx) , pointer , intent(inout) , dimension(:,:,:) :: pp
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
            tvkp1 = tv(j,i,k+1)
            tvk = tv(j,i,k)
            tkp1 = t(j,i,k+1)
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
    subroutine nhw(i1,i2,j1,j2,kxs,sigma,dsigma,u,v,tv, &
                   ps,psdot,ps0,xmsfx,w,wtop,ds,iband,icrm)
      implicit none
      integer(ik4) , intent(in) :: i1 , i2 , j1 , j2 , kxs , iband , icrm
      real(rkx) , pointer , intent(in) , dimension(:) :: sigma , dsigma
      real(rkx) , pointer , intent(in) , dimension(:,:) :: xmsfx
      real(rkx) , pointer , intent(in) , dimension(:,:) :: ps , ps0 , psdot
      real(rkx) , pointer , intent(in) , dimension(:,:,:) :: tv
      real(rkx) , pointer , intent(in) , dimension(:,:,:) :: u , v
      real(rkx) , intent(in) :: ds                    ! Kilometers
      real(rkx) , pointer , intent(inout) , dimension(:,:,:) :: w
      real(rkx) , pointer , intent(inout) , dimension(:,:) :: wtop
      integer(ik4) :: i , j , k
      integer(ik4) :: l , ll , ip , im , jp , jm , lm , lp
      real(rkx) :: dx2 , omegal , omegau , ubar , vbar , wu , wl
      real(rkx) :: zl , zu , rho , omegan , alnp
      real(rkx) :: ua , ub , va , vb
      real(rkx) , dimension(kxs) :: mdv
      real(rkx) , dimension(kxs+1) :: qdt
      real(rkx) , dimension(j1:j2,i1:i2,kxs+1) :: z0 , z
      real(rkx) , dimension(j1:j2,i1:i2,kxs+1) :: wtmp
      real(rkx) , dimension(j1:j2,i1:i2,kxs+1) :: pr0 , t0 , omega
      real(rkx) , dimension(j1:j2,i1:i2) :: dummy , dummy1

      wtmp(:,:,:) = d_zero
      omega(:,:,:) = d_zero
      dx2 = d_two * ds
      dummy = (xmsfx * xmsfx) / dx2
      dummy1 = xmsfx / dx2

      !
      ! We expect ps and ps0 to be already interpolated on dot points
      !
      do k = 1 , kxs+1
        do i = i1 , i2
          do j = j1 , j2
            pr0(j,i,k) = ps0(j,i) * sigma(k) + ptoppa
            t0(j,i,k) = max(st0 + tlp * log(pr0(j,i,k) / p0),tiso)
            alnp = log(pr0(j,i,k)/(ps0(j,i)+ptoppa))
            z0(j,i,k) = - (d_half*rovg*tlp*alnp*alnp + rovg*st0*alnp)
          end do
        end do
      end do
      !
      !  Calculate heights of input temperature sounding for interpolation
      !  to nonhydrostatic model levels.
      !
      do i = i1 , i2
        do j = j1 , j2
          z(j,i,kxs+1) = d_zero
          do k = kxs , 1 , -1
            z(j,i,k) = z(j,i,k+1) - rovg * tv(j,i,k) * &
              log((sigma(k)*ps(j,i)+ptop)/(sigma(k+1)*ps(j,i)+ptop))
          end do
        end do
      end do
      do i = i1 , i2
        do j = j1 , j2
          if ( icrm /= 1 ) then
            ip = min(i+1,i2)
            im = max(i-1,i1)
          else
            if ( i == i2-1 ) then
              ip = i2
            else if ( i == i2 ) then
              ip = i1
            else
              ip = i + 1
            end if
            if ( i == i1+1 ) then
              im = i1
            else if ( i == i1 ) then
              im = i2
            else
              im = i - 1
            end if
          end if
          if ( iband /= 1 ) then
            jp = min(j+1,j2)
            jm = max(j-1,j1)
          else
            if ( j == j2-1 ) then
              jp = j2
            else if ( j == j2 ) then
              jp = j1
            else
              jp = j + 1
            end if
            if ( j == j1+1 ) then
              jm = j1
            else if ( j == j1 ) then
              jm = j2
            else
              jm = j - 1
            end if
          end if

          mdv(:) = d_zero
          do l = 1 , kxs
            ua = u(j ,i ,l) * psdot(j,i)  + &
                 u(j ,ip,l) * psdot(j,ip)
            ub = u(jp, i,l) * psdot(jp,i) + &
                 u(jp,ip,l) * psdot(jp,ip)
            va = v(j ,i ,l) * psdot(j,i)  + &
                 v(jp,i ,l) * psdot(jp,i)
            vb = v(j ,ip,l) * psdot(j,ip) + &
                 v(jp,ip,l) * psdot(jp,ip)
            mdv(l) = (ub - ua + vb - va) * dummy(j,i) / ps(j,i)
          end do
          qdt(kxs+1) = d_zero
          do l = kxs , 1 , -1
            qdt(l) = qdt(l+1) + mdv(l) * dsigma(l)
          end do
          do l = kxs+1 , 1 , -1
            lp = min(l,kxs)
            lm = max(l-1,1)
            if ( l == kxs+1 ) lm = kxs-1
            ubar = 0.125_rkx * (u(j ,i ,lm) + u(j ,ip,lm) + &
                                u(jp,i ,lm) + u(jp,ip,lm) + &
                                u(j ,i ,lp) + u(j ,ip,lp) + &
                                u(jp,i ,lp) + u(jp,ip,lp))
            vbar = 0.125_rkx * (v(j ,i ,lm) + v(j ,ip,lm) + &
                                v(jp,i ,lm) + v(jp,ip,lm) + &
                                v(j ,i ,lp) + v(j ,ip,lp) + &
                                v(jp,i ,lp) + v(jp,ip,lp))
            ! Calculate omega
            omega(j,i,l) = ps(j,i) * qdt(l) + sigma(l) *  &
                    ((ps(jp,i) - ps(jm,i)) * ubar + &
                     (ps(j,ip) - ps(j,im)) * vbar) * dummy1(j,i)
          end do
        end do
      end do
      !
      ! Remove signal from grid (need because interpolation in ATM)
      !
      call smtdsmt(omega,i1,i2,j1,j2,1,kxs+1)
      !
      ! Vertical velocity from interpolated omega
      !
      do k = 2 , kxs + 1
        do i = i1 , i2
          do j = j1 , j2
            do ll = 1 , kxs
              l = ll
              if (z(j,i,l+1) < z0(j,i,k)) exit
            end do
            zu = z(j,i,l)
            zl = z(j,i,l+1)
            omegau = omega(j,i,l)
            omegal = omega(j,i,l+1)
            wu = (z0(j,i,k) - zl) / (zu - zl)
            wl = d_one - wu
            omegan = omegau * wu + omegal * wl
            rho = pr0(j,i,k) / rgas / t0(j,i,k)
            ! W =~ -OMEGA/RHO0/G *1000*PS0/1000. (OMEGA IN CB)
            wtmp(j,i,k) = -d_1000 * omegan/rho * regrav
          end do
        end do
      end do
      wtmp(j1,:,:) = wtmp(j1+1,:,:)
      wtmp(j2-1,:,:) = wtmp(j2-2,:,:)
      wtmp(:,i1,:) = wtmp(:,i1+1,:)
      wtmp(:,i2-1,:) = wtmp(:,i2-2,:)
      wtmp(j2,:,:) = wtmp(j2-1,:,:)
      wtmp(:,i2,:) = wtmp(:,i2-1,:)
      wtop(j1:j2,i1:i2) = wtmp(j1:j2,i1:i2,1)
      w(j1:j2,i1:i2,1:kxs) = wtmp(j1:j2,i1:i2,2:kxs+1)
    end subroutine nhw

    subroutine smtdsmt(slab,i1,i2,j1,j2,k1,k2)
      implicit none
      integer(ik4) , intent(in) :: i1 , i2 , j1 , j2 , k1 , k2
      real(rkx) , intent(inout) , dimension(j1:j2,i1:i2,k1:k2) :: slab
      real(rkx) :: aplus , asv , cell
      integer(ik4) :: i , is , ie , j , js , je , k , kp , np
      real(rkx) , dimension(2) :: xnu
      integer(ik4) , parameter :: npass = 16
      !
      ! purpose: spatially smooth data in slab to dampen short
      ! wavelength components
      !
      ie = i2-1
      je = j2-1
      is = i1+1
      js = j1+1
      xnu(1) =  0.50_rkx
      xnu(2) = -0.52_rkx
      do k = k1 , k2
        do np = 1 , npass
          do kp = 1 , 2
            ! first smooth in the ni direction
            do i = i1 , i2
              asv = slab(j1,i,k)
              do j = js , je
                cell = slab(j,i,k)
                aplus = slab(j+1,i,k)
                slab(j,i,k) = cell + xnu(kp)*( (asv+aplus)/d_two - cell)
                asv = cell
              end do
            end do
            ! smooth in the nj direction
            do j = j1 , j2
              asv = slab(j,i1,k)
              do i = is , ie
                cell = slab(j,i,k)
                aplus = slab(j,i+1,k)
                slab(j,i,k) = cell + xnu(kp)*((asv+aplus)/d_two - cell)
                asv = cell
              end do
            end do
          end do
          slab(j1,:,k) = slab(j1+1,:,k)
          slab(:,i1,k) = slab(:,i1+1,k)
          slab(j2-1,:,k) = slab(j2-2,:,k)
          slab(:,i2-1,k) = slab(:,i2-2,k)
          slab(j2,:,k) = slab(j2-1,:,k)
          slab(:,i2,k) = slab(:,i2-1,k)
        end do
      end do
    end subroutine smtdsmt

end module mod_nhinterp

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

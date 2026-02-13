!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_pbl_gfs

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_memutil
  use mod_dynparam, only : kz, kzp1, kzm1, ntr, idynamic
  use mod_dynparam, only : ici1, ici2, jci1, jci2
  use mod_runparams, only : dt, nqx, ichem, iqv
  use mod_regcm_types, only : mod_2_pbl, pbl_2_mod

  implicit none

  private

  public :: init_pbl_gfs, pbl_gfs

  integer(ik4) :: iblp, ibnt

  ! Input to moninq

  real(rkx), dimension(:,:), pointer, contiguous :: uo, vo, to, xpfac
  real(rkx), dimension(:,:,:), pointer, contiguous :: qo
  real(rkx), dimension(:), pointer, contiguous :: psk, rbsoil
  real(rkx), dimension(:), pointer, contiguous :: fm, fh, spd1
  real(rkx), dimension(:,:), pointer, contiguous :: prsi, phii, z, dz
  real(rkx), dimension(:,:), pointer, contiguous :: del, prsl, prslk, phil, thraten
  real(rkx), dimension(:), pointer, contiguous :: heat, evap, stress, ustar, rrho

  ! Output from moninq

  real(rkx), dimension(:), pointer, contiguous :: hpbl
  integer(ik4), dimension(:), pointer, contiguous :: kpbl
  real(rkx), dimension(:,:), pointer, contiguous :: dv, du
  real(rkx), dimension(:,:,:), pointer, contiguous :: rtg
  real(rkx), dimension(:,:), pointer, contiguous :: tau

  contains

    subroutine init_pbl_gfs
      implicit none

      iblp = (jci2-jci1+1)*(ici2-ici1+1)
      ibnt = nqx+ntr

      call getmem(hpbl,1,iblp,'mod_pbl_gfs:hpbl')
      call getmem(kpbl,1,iblp,'mod_pbl_gfs:kpbl')
      call getmem(du,1,iblp,1,kz,'mod_pbl_gfs:du')
      call getmem(dv,1,iblp,1,kz,'mod_pbl_gfs:dv')
      call getmem(tau,1,iblp,1,kz,'mod_pbl_gfs:tau')
      call getmem(rtg,1,iblp,1,kz,1,ibnt,'mod_pbl_gfs:rtg')

      call getmem(uo,1,iblp,1,kz,'mod_pbl_gfs:uo')
      call getmem(vo,1,iblp,1,kz,'mod_pbl_gfs:vo')
      call getmem(to,1,iblp,1,kz,'mod_pbl_gfs:to')
      call getmem(qo,1,iblp,1,kz,1,ibnt,'mod_pbl_gfs:qo')
      call getmem(psk,1,iblp,'mod_pbl_gfs:psk')
      call getmem(rbsoil,1,iblp,'mod_pbl_gfs:rbsoil')
      call getmem(fm,1,iblp,'mod_pbl_gfs:fm')
      call getmem(fh,1,iblp,'mod_pbl_gfs:fh')
      call getmem(spd1,1,iblp,'mod_pbl_gfs:spd1')
      call getmem(prsi,1,iblp,1,kzp1,'mod_pbl_gfs:prsi')
      call getmem(phii,1,iblp,1,kzp1,'mod_pbl_gfs:phii')
      call getmem(del,1,iblp,1,kz,'mod_pbl_gfs:del')
      call getmem(prsl,1,iblp,1,kz,'mod_pbl_gfs:prsl')
      call getmem(prslk,1,iblp,1,kz,'mod_pbl_gfs:prslk')
      call getmem(phil,1,iblp,1,kz,'mod_pbl_gfs:phil')
      call getmem(thraten,1,iblp,1,kz,'mod_pbl_gfs:thraten')
      call getmem(z,1,iblp,1,kz,'mod_pbl_gfs:z')
      call getmem(dz,1,iblp,1,kz,'mod_pbl_gfs:dz')
      call getmem(heat,1,iblp,'mod_pbl_gfs:heat')
      call getmem(evap,1,iblp,'mod_pbl_gfs:evap')
      call getmem(rrho,1,iblp,'mod_pbl_gfs:rrho')
      call getmem(stress,1,iblp,'mod_pbl_gfs:stress')
      call getmem(ustar,1,iblp,'mod_pbl_gfs:ustar')
      call getmem(xpfac,jci1,jci2,ici1,ici2,'mod_pbl_gfs:xpfac')

    end subroutine init_pbl_gfs

    subroutine pbl_gfs(m2p,p2m)
      implicit none
      type(mod_2_pbl), intent(in) :: m2p
      type(pbl_2_mod), intent(inout) :: p2m

      integer(ik4) :: i, j, k, kk, km, n
      integer(ik4) :: iq, it, iit
      real(rkx) :: ps, ua, va

      n = 1
      do i = ici1, ici2
        do j = jci1, jci2
          ua = max(m2p%u10m(j,i),0.05_rkx)
          va = max(m2p%v10m(j,i),0.05_rkx)
          ps = m2p%patmf(j,i,kzp1)
          fm(n) = m2p%ram1(j,i)
          fh(n) = m2p%rah1(j,i)
          rbsoil(n) = m2p%br(j,i)
          psk(n) = (ps/p00)**rovcp
          ustar(n) = m2p%ustar(j,i)
          rrho(n) = 1.0_rkx/m2p%rhox2d(j,i)
          stress(n) = (ustar(n)*ustar(n))*m2p%rhox2d(j,i)
          heat(n) = m2p%hfx(j,i)*rcpd*rrho(n)
          evap(n) = m2p%qfx(j,i)
          spd1(n) = sqrt(ua*ua+va*va)
          prsi(n,1) = ps*d_r1000
          phii(n,1) = d_zero
          n = n + 1
        end do
      end do

      do k = 1, kz
        n = 1
        kk = kzp1 - k
        do i = ici1, ici2
          do j = jci1, jci2
            to(n,kk) = m2p%tatm(j,i,k)
            uo(n,kk) = m2p%uxatm(j,i,k)
            vo(n,kk) = m2p%vxatm(j,i,k)
            prsl(n,kk) = m2p%patm(j,i,k)*d_r1000
            prslk(n,kk) = (m2p%patm(j,i,k)/p00)**rovcp
            dz(n,kk) = m2p%dzq(j,i,k)
            z(n,kk) = m2p%za(j,i,k)
            thraten(n,kk) = m2p%heatrt(j,i,k)/prslk(n,kk)
            n = n + 1
          end do
        end do
      end do

      do k = 2, kz
        n = 1
        km = k - 1
        do i = ici1, ici2
          do j = jci1, jci2
            del(n,km) = prsl(n,km)/rovg*dz(n,km)/to(n,km)
            prsi(n,k) = prsi(n,km)-del(n,km)
            phii(n,k) = (z(n,k)-z(n,1))*egrav
            phil(n,km) = d_half*(z(n,k)+z(n,km)-d_two*z(n,1))*egrav
            n = n + 1
          end do
        end do
      end do

      n = 1
      do i = ici1, ici2
        do j = jci1, jci2
          del(n,kz) = del(n,kzm1)
          prsi(n,kzp1) = prsi(n,kz)-del(n,kzm1)
          phii(n,kzp1) = phii(n,kz)+dz(n,kz)*egrav
          phil(n,kz) = phii(n,kz)-phil(n,kzm1)+phii(n,kz)
          n = n + 1
        end do
      end do

      do iq = 1, nqx
        do k = 1, kz
          n = 1
          kk = kzp1 - k
          do i = ici1, ici2
            do j = jci1, jci2
              qo(n,kk,iq) = m2p%qxatm(j,i,k,iq)/(d_one + m2p%qxatm(j,i,k,iq))
              n = n + 1
            end do
          end do
        end do
      end do

      if ( ichem == 1 ) then
        do it = 1, ntr
          iit = it + nqx
          do k = 1, kz
            n = 1
            kk = kzp1 - k
            do i = ici1, ici2
              do j = jci1, jci2
                qo(n,kk,iit) = m2p%chib(j,i,k,it)/(d_one + m2p%chib(j,i,k,it))
                n = n + 1
              end do
            end do
          end do
        end do
      end if

      du(:,:) = d_zero
      dv(:,:) = d_zero
      tau(:,:) = d_zero
      rtg(:,:,:) = d_zero
      hpbl(:) = d_zero
      kpbl(:) = 1

      call moninq(iblp,kz,ibnt)

      if ( idynamic == 3 ) then
        xpfac = d_one
      else
        xpfac = m2p%psb(jci1:jci2,ici1:ici2)
      end if
      do k = 1, kz
        n = 1
        kk = kzp1 - k
        do i = ici1, ici2
          do j = jci1, jci2
            p2m%uxten(j,i,k) = du(n,kk)
            p2m%vxten(j,i,k) = dv(n,kk)
            p2m%tten(j,i,k) = p2m%tten(j,i,k) + tau(n,kk)*xpfac(j,i)
            n = n + 1
          end do
        end do
      end do

      do iq = 1, nqx
        do k = 1, kz
          n = 1
          kk = kzp1 - k
          do i = ici1, ici2
            do j = jci1, jci2
              p2m%qxten(j,i,k,iq) = p2m%qxten(j,i,k,iq) + &
                     rtg(n,kk,iq)/(d_one-qo(n,kk,iq))**2 * xpfac(j,i)
              n = n + 1
            end do
          end do
        end do
      end do

      if ( ichem == 1 ) then
        do it = 1, ntr
          iit = it + nqx
          do k = 1, kz
            n = 1
            kk = kzp1 - k
            do i = ici1, ici2
              do j = jci1, jci2
                p2m%chiten(j,i,k,it) = p2m%chiten(j,i,k,it) + &
                      rtg(n,kk,iit)/(d_one-qo(n,kk,iit)) * xpfac(j,i)
                n = n + 1
              end do
            end do
          end do
        end do
      end if

      n = 1
      do i = ici1, ici2
        do j = jci1, jci2
          p2m%kpbl(j,i) = kzp1 - kpbl(n)
          p2m%zpbl(j,i) = hpbl(n)
          n = n + 1
        end do
      end do

    end subroutine pbl_gfs

    subroutine moninq(im,km,ntrac)
      implicit none

      integer(ik4), intent(in) :: im, km, ntrac
      integer(ik4) :: i, is, k, kk, km1, kmpbl
      integer(ik4), dimension(im) :: lcld, icld, kcld, krad, kx1
      integer(ik4), dimension(im) :: kpblx, kinver

      real(rkx), dimension(im) :: phih, phim
      real(rkx), dimension(im) :: rbdn, rbup, beta
      real(rkx), dimension(im) :: wscale, thermal, wstar3
      real(rkx), dimension(im) :: hpblx

      real(rkx), dimension(im,km) :: thvx, thlvx
      real(rkx), dimension(im,km) :: qlx, thetae, qtx
      real(rkx), dimension(im,km-1) :: bf, radx
      real(rkx), dimension(im) :: govrth, hrad, cteit
      real(rkx), dimension(im) :: radmin, vrad, zd, zdd, thlvx1
      real(rkx), dimension(im) :: hgamt, hgamq, tx1, tx2

      real(rkx), dimension(im,km-1) :: rdzt, dktx, xkzo, xkzmo
      real(rkx), dimension(im,km+1) :: zi
      real(rkx), dimension(im,km-1) :: dku, dkt, cku, ckt, al, au
      real(rkx), dimension(im,km) :: zl, ad, a1, theta
      real(rkx), dimension(im,km*ntrac) :: a2
      real(rkx), dimension(im) :: prinv, rent
      logical, dimension(im) :: pblflg, sfcflg, scuflg, flg
      real(rkx) :: bvf2, dk, dsdz2, dsdzq, dsdzt
      real(rkx) :: dsig, dtodsd, dtodsu, dw2, hol, hol1
      real(rkx) :: prnum, qtend, rbint, rdt, rdz, ri, rl2
      real(rkx) :: sflux, shr2, spdk2, sri, tem, ti
      real(rkx) :: ttend, utend, vtend, zfac, vpert, zk
      real(rkx) :: tem1, tem2, ptem, ptem1, ptem2

      real(rkx), parameter :: gocp = egrav*rcpd
      real(rkx), parameter :: rlam = 30.0_rkx
      real(rkx), parameter :: vk = vonkar
      real(rkx), parameter :: prmin = 0.25_rkx
      real(rkx), parameter :: prmax = 4.0_rkx
      real(rkx), parameter :: dw2min = 0.0001_rkx
      real(rkx), parameter :: dkmin = 0.0_rkx
      real(rkx), parameter :: dkmax = 1000.0_rkx
      real(rkx), parameter :: rimin = -100.0_rkx
      real(rkx), parameter :: rbcr = 0.25_rkx
      real(rkx), parameter :: wfac = 7.0_rkx
      real(rkx), parameter :: cfac = 6.5_rkx
      real(rkx), parameter :: pfac = 2.0_rkx
      real(rkx), parameter :: sfcfrac = 0.1_rkx
      real(rkx), parameter :: qmin = 1.e-8_rkx
      real(rkx), parameter :: xkzm = 1.0_rkx
      real(rkx), parameter :: zfmin = 1.e-8_rkx
      real(rkx), parameter :: aphi5 = 5.0_rkx
      real(rkx), parameter :: aphi16 = 16.0_rkx
      real(rkx), parameter :: tdzmin = 1.e-3_rkx
      real(rkx), parameter :: qlmin = 1.e-12_rkx
      real(rkx), parameter :: h1 = 0.33333333_rkx
      real(rkx), parameter :: cldtime = 500.0_rkx
      real(rkx), parameter :: xkzmu = 3.0_rkx
      real(rkx), parameter :: xkzminv = 0.3_rkx
      real(rkx), parameter :: gamcrt = 3.0_rkx
      real(rkx), parameter :: gamcrq = 0.0_rkx
      real(rkx), parameter :: rlamun = 150.0_rkx
      real(rkx), parameter :: rentf1 = 0.2_rkx
      real(rkx), parameter :: rentf2 = 1.0_rkx
      real(rkx), parameter :: radfac = 0.85_rkx
      real(rkx), parameter :: zstblmax = 2500.0_rkx
      real(rkx), parameter :: qlcr=3.5e-5_rkx
      real(rkx), parameter :: actei = 0.7_rkx

      rdt   = d_one / dt
      km1   = km - 1
      kmpbl = km / 2
      do k = 1, km
        do i = 1, im
          zi(i,k) = phii(i,k) * regrav
          zl(i,k) = phil(i,k) * regrav
        end do
      end do
      do i = 1, im
        zi(i,km+1) = phii(i,km+1) * regrav
      end do
      do k = 1, km1
        do i = 1, im
          rdzt(i,k) = d_one / (zl(i,k+1) - zl(i,k))
        end do
      end do
      do i = 1, im
        kinver(i) = km
        kx1(i) = 1
        tx1(i) = 1.0 / prsi(i,1)
        tx2(i) = tx1(i)
      end do
      ! Vertical background diffusivity
      do k = 1, km1
        do i = 1, im
          xkzo(i,k) = d_zero
          xkzmo(i,k) = d_zero
          if ( k < kinver(i) ) then
            ptem      = prsi(i,k+1) * tx1(i)
            tem1      = d_one - ptem
            tem1      = min(tem1 * tem1 * 10.0_rkx,25.0_rkx)
            xkzo(i,k) = xkzm * min(d_one, exp(-tem1))
            ! Vertical background diffusivity for momentum
            if ( ptem >= 0.2_rkx ) then
              xkzmo(i,k) = xkzmu
              kx1(i)     = k + 1
            else
              if (k == kx1(i) .and. k > 1) tx2(i) = 1.0 / prsi(i,k)
              tem1 = d_one - prsi(i,k+1) * tx2(i)
              tem1 = min(tem1 * tem1 * 5.0_rkx,25.0_rkx)
              xkzmo(i,k) = xkzmu * min(d_one, exp(-tem1))
            end if
          end if
        end do
      end do
      ! Diffusivity in the inversion layer is set to be xkzminv (m^2/s)
      do k = 1, kmpbl
        do i = 1, im
          if ( zi(i,k+1) > 250.0_rkx ) then
            tem1 = (to(i,k+1)-to(i,k)) * rdzt(i,k)
            if ( tem1 > 1.e-5_rkx ) then
              xkzo(i,k) = min(xkzo(i,k),xkzminv)
            end if
          end if
        end do
      end do
      do i = 1, im
        hgamt(i) = d_zero
        hgamq(i) = d_zero
        wscale(i)= d_zero
        kpbl(i)  = 1
        kpblx(i) = 1
        hpbl(i)  = zi(i,1)
        hpblx(i) = zi(i,1)
        pblflg(i)= .true.
        sfcflg(i)= .true.
        if ( rbsoil(i) > d_zero ) sfcflg(i) = .false.
        scuflg(i)= .true.
        if ( scuflg(i) ) then
          radmin(i)= d_zero
          cteit(i) = d_zero
          rent(i)  = rentf1
          hrad(i)  = zi(i,1)
          krad(i)  = 1
          icld(i)  = 0
          lcld(i)  = km1
          kcld(i)  = km1
          zd(i)    = d_zero
        end if
      end do
      do k = 1, km
        do i = 1, im
          theta(i,k) = to(i,k) * psk(i) / prslk(i,k)
          qlx(i,k)   = max(qo(i,k,2),qlmin)
          qtx(i,k)   = max(qo(i,k,1),qmin)+qlx(i,k)
          ptem       = qlx(i,k)
          ptem1      = wlhvocp*max(qo(i,k,1),qmin)/to(i,k)
          thetae(i,k)= theta(i,k)*(d_one+ptem1)
          thvx(i,k)  = theta(i,k)*(d_one+ep1*max(qo(i,k,1),qmin)-ptem)
          ptem2      = theta(i,k)-wlhvocp*ptem
          thlvx(i,k) = ptem2*(d_one+ep1*qtx(i,k))
        end do
      end do
      do k = 1, km1
        do i = 1, im
          dku(i,k)  = d_zero
          dkt(i,k)  = d_zero
          dktx(i,k) = d_zero
          cku(i,k)  = d_zero
          ckt(i,k)  = d_zero
          tem       = zi(i,k+1)-zi(i,k)
          radx(i,k) = tem*thraten(i,k)
        end do
      end do
      do i = 1, im
        flg(i) = scuflg(i)
      end do
      do k = 1, km1
        do i = 1, im
          if ( flg(i) .and. zl(i,k) >= zstblmax ) then
            lcld(i) = k
            flg(i) = .false.
          end if
        end do
      end do
      ! Compute buoyancy flux
      do k = 1, km1
        do i = 1, im
          bf(i,k) = (thvx(i,k+1)-thvx(i,k))*rdzt(i,k)
        end do
      end do
      do i = 1, im
        govrth(i) = egrav/theta(i,1)
      end do
      do i = 1, im
        beta(i) = dt / (zi(i,2)-zi(i,1))
      end do
      do i = 1, im
        thermal(i) = thvx(i,1)
      end do
      ! Compute the first guess pbl height
      do i = 1, im
        flg(i) = .false.
        rbup(i) = rbsoil(i)
      end do
      do k = 2, kmpbl
        do i = 1, im
          if ( .not. flg(i) ) then
            rbdn(i) = rbup(i)
            spdk2   = max((uo(i,k)**2+vo(i,k)**2),1.0_rkx)
            rbup(i) = (thvx(i,k)-thermal(i)) * &
                      (egrav*zl(i,k)/thvx(i,1))/spdk2
            kpbl(i) = k
            flg(i)  = (rbup(i) > rbcr)
          end if
        end do
      end do
      do i = 1, im
        k = kpbl(i)
        if ( rbdn(i) >= rbcr ) then
          rbint = d_zero
        else if ( rbup(i) <= rbcr ) then
          rbint = d_one
        else
          rbint = (rbcr-rbdn(i))/(rbup(i)-rbdn(i))
        end if
        hpbl(i) = zl(i,k-1) + rbint*(zl(i,k)-zl(i,k-1))
        if ( hpbl(i) < zi(i,kpbl(i)) ) kpbl(i) = kpbl(i) - 1
        hpblx(i) = hpbl(i)
        kpblx(i) = kpbl(i)
      end do
      ! Compute similarity parameters
      do i = 1, im
        sflux = heat(i) + evap(i)*ep1*theta(i,1)*rrho(i)
        if ( sfcflg(i) .and. sflux > d_zero ) then
          hol = max(rbsoil(i)*fm(i)*fm(i)/fh(i),rimin)
          hol = min(hol,-zfmin)
          hol1 = hol*hpbl(i)/zl(i,1)*sfcfrac
          tem     = d_one / (d_one - aphi16*hol1)
          phih(i) = sqrt(tem)
          phim(i) = sqrt(phih(i))
          wstar3(i) = govrth(i)*sflux*hpbl(i)
          tem1      = ustar(i)**3
          wscale(i) = (tem1+wfac*vk*wstar3(i)*sfcfrac)**h1
          wscale(i) = min(wscale(i),ustar(i)*aphi16)
          wscale(i) = max(wscale(i),ustar(i)/aphi5)
        else
          pblflg(i) = .false.
        end if
      end do
      ! Compute counter-gradient mixing term for heat and moisture
      do i = 1, im
        if ( pblflg(i) ) then
          hgamt(i)  = min(cfac*heat(i)/wscale(i),gamcrt)
          hgamq(i)  = min(cfac*evap(i)/wscale(i),gamcrq)
          vpert     = hgamt(i) + hgamq(i)*ep1*theta(i,1)
          vpert     = min(vpert,gamcrt)
          thermal(i)= thermal(i)+max(vpert,d_zero)
          hgamt(i)  = max(hgamt(i),d_zero)
          hgamq(i)  = max(hgamq(i),d_zero)
        end if
      end do
      ! Enhance the pbl height by considering the thermal excess
      do i = 1, im
        flg(i)  = .true.
        if ( pblflg(i) ) then
          flg(i)  = .false.
          rbup(i) = rbsoil(i)
        end if
      end do
      do k = 2, kmpbl
        do i = 1, im
          if ( .not. flg(i) ) then
            rbdn(i) = rbup(i)
            spdk2   = max((uo(i,k)**2+vo(i,k)**2),1.0_rkx)
            rbup(i) = (thvx(i,k)-thermal(i)) * &
                      (egrav*zl(i,k)/thvx(i,1))/spdk2
            kpbl(i) = k
            flg(i)  = (rbup(i) > rbcr)
          end if
        end do
      end do
      do i = 1, im
        if ( pblflg(i) ) then
          k = kpbl(i)
          if ( rbdn(i) >= rbcr ) then
            rbint = d_zero
          else if ( rbup(i) <= rbcr ) then
            rbint = d_one
          else
            rbint = (rbcr-rbdn(i))/(rbup(i)-rbdn(i))
          end if
          hpbl(i) = zl(i,k-1) + rbint*(zl(i,k)-zl(i,k-1))
          if ( hpbl(i) < zi(i,kpbl(i)) ) kpbl(i) = kpbl(i) - 1
          if ( kpbl(i) >= 1 ) pblflg(i) = .false.
        end if
      end do
      ! Look for stratocumulus
      do i = 1, im
        flg(i) = scuflg(i)
      end do
      do k = kmpbl, 1, -1
        do i = 1, im
          if ( flg(i) .and. k <= lcld(i) ) then
            if ( qlx(i,k) >= qlcr ) then
              kcld(i) = k
              flg(i) = .false.
            end if
          end if
        end do
      end do
      do i = 1, im
        if ( scuflg(i) .and. kcld(i) == km1 ) scuflg(i) = .false.
      end do
      do i = 1, im
        flg(i) = scuflg(i)
      end do
      do k = kmpbl, 1, -1
        do i = 1, im
          if ( flg(i) .and. k <= kcld(i) ) then
            if ( qlx(i,k) >= qlcr ) then
              if ( radx(i,k) < radmin(i) ) then
                radmin(i) = radx(i,k)
                krad(i) = k
              end if
            else
              flg(i) = .false.
            end if
          end if
        end do
      end do
      do i = 1, im
        if ( scuflg(i) .and. krad(i) <= 1 ) scuflg(i) = .false.
        if ( scuflg(i) .and. radmin(i) >= d_zero ) scuflg(i) = .false.
      end do
      do i = 1, im
        flg(i) = scuflg(i)
      end do
      do k = kmpbl, 2, -1
        do i = 1, im
          if ( flg(i) .and. k <= krad(i) ) then
            if ( qlx(i,k) >= qlcr ) then
              icld(i) = icld(i)+1
            else
              flg(i) = .false.
            end if
          end if
        end do
      end do
      do i = 1, im
        if ( scuflg(i) .and. icld(i) < 1 ) scuflg(i) = .false.
      enddo
      do i = 1, im
        if ( scuflg(i) ) then
          hrad(i) = zi(i,krad(i)+1)
        end if
      end do
      do i = 1, im
        if ( scuflg(i) .and. hrad(i) < zi(i,2) ) scuflg(i) = .false.
      end do
      do i = 1, im
        if ( scuflg(i) ) then
          k    = krad(i)
          tem  = zi(i,k+1)-zi(i,k)
          tem1 = cldtime*radmin(i)/tem
          thlvx1(i) = thlvx(i,k)+tem1
        end if
      end do
      do i = 1, im
        flg(i) = scuflg(i)
      end do
      do k = kmpbl, 1, -1
        do i = 1, im
          if ( flg(i) .and. k <= krad(i) ) then
            if (thlvx1(i) <= thlvx(i,k) ) then
              tem = zi(i,k+1)-zi(i,k)
              zd(i) = zd(i)+tem
            else
              flg(i) = .false.
            end if
          end if
        end do
      end do
      do i = 1, im
        if ( scuflg(i) ) then
          kk = max(1, krad(i)+1-icld(i))
          zdd(i) = hrad(i)-zi(i,kk)
        end if
      end do
      do i = 1, im
        if ( scuflg(i) ) then
          zd(i) = max(zd(i),zdd(i))
          zd(i) = min(zd(i),hrad(i))
          tem   = govrth(i)*zd(i)*(-radmin(i))
          vrad(i) = tem**h1
        end if
      end do
      ! Compute inverse Prandtl number
      do i = 1, im
        if ( pblflg(i) ) then
          tem = phih(i)/phim(i)+cfac*vk*sfcfrac
          prinv(i) = d_one / tem
          prinv(i) = min(prinv(i),prmax)
          prinv(i) = max(prinv(i),prmin)
        end if
      end do
      ! Compute diffusion coefficients below pbl
      do k = 1, kmpbl
        do i = 1, im
          if ( pblflg(i) .and. k < kpbl(i) ) then
            zfac = max((d_one-zi(i,k+1)/(hpbl(i))), zfmin)
            tem = wscale(i)*vk*zi(i,k+1)*zfac**pfac
            dku(i,k) = xkzmo(i,k) + tem
            dkt(i,k) = xkzo(i,k) + tem * prinv(i)
            dku(i,k) = min(dku(i,k),dkmax)
            dkt(i,k) = min(dkt(i,k),dkmax)
            dktx(i,k)= dkt(i,k)
          end if
        end do
      end do
      ! Compute diffusion coefficients based on local scheme
      do i = 1, im
        if ( .not. pblflg(i) ) then
          kpbl(i) = 1
        end if
      end do
      do k = 1, km1
        do i = 1, im
          if ( k >= kpbl(i) ) then
            rdz  = rdzt(i,k)
            ti   = d_two/(to(i,k)+to(i,k+1))
            dw2  = (uo(i,k)-uo(i,k+1))**2 +(vo(i,k)-vo(i,k+1))**2
            shr2 = max(dw2,dw2min)*rdz*rdz
            bvf2 = egrav*bf(i,k)*ti
            ri   = max(bvf2/shr2,rimin)
            zk   = vk*zi(i,k+1)
            if (ri < d_zero ) then ! Unstable regime
              rl2      = zk*rlamun/(rlamun+zk)
              dk       = rl2*rl2*sqrt(shr2)
              sri      = sqrt(-ri)
              dku(i,k) = xkzmo(i,k) + &
                      dk*(d_one+8.0_rkx*(-ri)/(d_one+1.746_rkx*sri))
              dkt(i,k) = xkzo(i,k) + &
                      dk*(d_one+8.0_rkx*(-ri)/(d_one+1.286_rkx*sri))
            else ! Stable regime
              rl2      = zk*rlam/(rlam+zk)
              dk       = rl2*rl2*sqrt(shr2)
              tem1     = dk/(d_one+5.0_rkx*ri)**2
              if ( k >= kpblx(i) ) then
                prnum = d_one + 2.1_rkx*ri
                prnum = min(prnum,prmax)
              else
                prnum = d_one
              end if
              dkt(i,k) = xkzo(i,k) + tem1
              dku(i,k) = xkzmo(i,k) + tem1 * prnum
            end if
            dku(i,k) = min(dku(i,k),dkmax)
            dkt(i,k) = min(dkt(i,k),dkmax)
          end if
        end do
      end do
      ! Compute diffusion coefficients for cloud-top driven diffusion
      ! If the condition for cloud-top instability is met,
      ! increase entrainment flux at cloud top
      do i = 1, im
        if ( scuflg(i) ) then
          k = krad(i)
          tem = thetae(i,k) - thetae(i,k+1)
          tem1 = qtx(i,k) - qtx(i,k+1)
          if ( tem > d_zero .and. tem1 > d_zero ) then
            cteit(i)= cpd*tem/(wlhv*tem1)
            if ( cteit(i) > actei ) rent(i) = rentf2
          end if
        end if
      end do
      do i = 1, im
        if ( scuflg(i) ) then
          k = krad(i)
          tem1  = max(bf(i,k),tdzmin)
          ckt(i,k) = -rent(i)*radmin(i)/tem1
          cku(i,k) = ckt(i,k)
        end if
      end do
      do k = 1, kmpbl
        do i = 1, im
          if ( scuflg(i) .and. k < krad(i) ) then
            tem1 = hrad(i)-zd(i)
            tem2 = zi(i,k+1)-tem1
            if ( tem2 > d_zero ) then
              ptem = tem2/zd(i)
              if ( ptem >= d_one ) ptem = d_one
              ptem = tem2*ptem*sqrt(d_one-ptem)
              ckt(i,k) = radfac*vk*vrad(i)*ptem
              cku(i,k) = 0.75_rkx*ckt(i,k)
              ckt(i,k) = max(ckt(i,k),dkmin)
              ckt(i,k) = min(ckt(i,k),dkmax)
              cku(i,k) = max(cku(i,k),dkmin)
              cku(i,k) = min(cku(i,k),dkmax)
            end if
          end if
        end do
      end do
      do k = 1, kmpbl
        do i = 1, im
          if ( scuflg(i) ) then
            dkt(i,k) = dkt(i,k)+ckt(i,k)
            dku(i,k) = dku(i,k)+cku(i,k)
            dkt(i,k) = min(dkt(i,k),dkmax)
            dku(i,k) = min(dku(i,k),dkmax)
          end if
        end do
      end do
      ! Compute tridiagonal matrix elements for heat and moisture
      do i = 1, im
         ad(i,1) = d_one
         a1(i,1) = to(i,1)   + beta(i) * heat(i)
         a2(i,1) = qo(i,1,1) + beta(i) * evap(i)
      end do
      if ( ntrac >= 2 ) then
        do k = 2, ntrac
          is = (k-1) * km
          do i = 1, im
            a2(i,1+is) = qo(i,1,k)
          end do
        end do
      end if
      do k = 1, km1
        do i = 1, im
          dtodsd = dt/del(i,k)
          dtodsu = dt/del(i,k+1)
          dsig   = prsl(i,k)-prsl(i,k+1)
          rdz    = rdzt(i,k)
          tem1   = dsig * dkt(i,k) * rdz
          if ( pblflg(i) .and. k < kpbl(i) ) then
             ptem1 = dsig * dktx(i,k) * rdz
             tem   = d_one / hpbl(i)
             dsdzt = tem1 * gocp - ptem1*hgamt(i)*tem
             dsdzq = ptem1 * (-hgamq(i)*tem)
             a2(i,k)   = a2(i,k)+dtodsd*dsdzq
             a2(i,k+1) = qo(i,k+1,1)-dtodsu*dsdzq
          else
             dsdzt = tem1 * gocp
             a2(i,k+1) = qo(i,k+1,1)
          endif
          dsdz2     = tem1 * rdz
          au(i,k)   = -dtodsd*dsdz2
          al(i,k)   = -dtodsu*dsdz2
          ad(i,k)   = ad(i,k)-au(i,k)
          ad(i,k+1) = d_one-al(i,k)
          a1(i,k)   = a1(i,k)+dtodsd*dsdzt
          a1(i,k+1) = to(i,k+1)-dtodsu*dsdzt
        end do
      end do
      if ( ntrac >= 2 ) then
        do kk = 2, ntrac
          is = (kk-1) * km
          do k = 1, km1
            do i = 1, im
              a2(i,k+1+is) = qo(i,k+1,kk)
            end do
          end do
        end do
      end if
      ! Solve tridiagonal problem for heat and moisture
      call tridin(im,km,ntrac,al,ad,au,a1,a2,au,a1,a2)
      ! Recover tendencies of heat and moisture
      do k = 1,km
        do i = 1, im
          ttend      = (a1(i,k)-to(i,k))*rdt
          qtend      = (a2(i,k)-qo(i,k,1))*rdt
          tau(i,k)   = tau(i,k)+ttend
          rtg(i,k,1) = rtg(i,k,1)+qtend
        end do
      end do
      if ( ntrac >= 2 ) then
        do kk = 2, ntrac
          is = (kk-1) * km
          do k = 1, km
            do i = 1, im
              qtend = (a2(i,k+is)-qo(i,k,kk))*rdt
              rtg(i,k,kk) = rtg(i,k,kk)+qtend
            end do
          end do
        end do
      end if
      ! Compute tridiagonal matrix elements for momentum
      do i = 1, im
        ad(i,1) = d_one + beta(i) * stress(i) / spd1(i)
        a1(i,1) = uo(i,1)
        a2(i,1) = vo(i,1)
      end do
      do k = 1, km1
        do i = 1, im
          dtodsd = dt/del(i,k)
          dtodsu = dt/del(i,k+1)
          dsig   = prsl(i,k)-prsl(i,k+1)
          rdz    = rdzt(i,k)
          tem1   = dsig*dku(i,k)*rdz
          a1(i,k+1) = uo(i,k+1)
          a2(i,k+1) = vo(i,k+1)
          dsdz2     = tem1*rdz
          au(i,k)   = -dtodsd*dsdz2
          al(i,k)   = -dtodsu*dsdz2
          ad(i,k)   = ad(i,k)-au(i,k)
          ad(i,k+1) = d_one-al(i,k)
        end do
      end do
      ! Solve tridiagonal problem for momentum
      call tridi2(im,km,al,ad,au,a1,a2,au,a1,a2)
      ! Recover tendencies of momentum
      do k = 1, km
        do i = 1, im
          utend = (a1(i,k)-uo(i,k))*rdt
          vtend = (a2(i,k)-vo(i,k))*rdt
          du(i,k)  = du(i,k)+utend
          dv(i,k)  = dv(i,k)+vtend
        end do
      end do
      ! PBL height for diagnostic purpose
      do i = 1, im
        hpbl(i) = hpblx(i)
        kpbl(i) = kpblx(i)
      end do

    end subroutine moninq

    subroutine tridi2(l,n,cl,cm,cu,r1,r2,au,a1,a2)
      implicit none
      integer(ik4), intent(in) :: l, n
      real(kind=rkx), dimension(l,2:n), intent(in) :: CL
      real(kind=rkx), dimension(l,n), intent(in) :: CM
      real(kind=rkx), dimension(l,n-1), intent(in) :: CU
      real(kind=rkx), dimension(l,n), intent(in) :: R1
      real(kind=rkx), dimension(l,n), intent(in) :: R2
      real(kind=rkx), dimension(l,n-1), intent(inout) :: au
      real(kind=rkx), dimension(l,n), intent(inout) :: a1
      real(kind=rkx), dimension(l,n), intent(inout) :: a2
      real(kind=rkx) :: fk
      integer(ik4) :: k, i

      do i = 1, l
        fk      = d_one/cm(i,1)
        au(i,1) = fk*cu(i,1)
        a1(i,1) = fk*r1(i,1)
        a2(i,1) = fk*r2(i,1)
      end do
      do k = 2, n-1
        do i = 1, l
          fk      = d_one/(cm(i,k)-cl(i,k)*au(i,k-1))
          au(i,k) = fk*cu(i,k)
          a1(i,k) = fk*(r1(i,k)-cl(i,k)*a1(i,k-1))
          a2(i,k) = fk*(r2(i,k)-cl(i,k)*a2(i,k-1))
        end do
      end do
      do i = 1, l
        fk      = d_one/(cm(i,n)-cl(i,n)*au(i,n-1))
        a1(i,n) = fk*(r1(i,n)-cl(i,n)*a1(i,n-1))
        a2(i,n) = fk*(r2(i,n)-cl(i,n)*a2(i,n-1))
      end do
      do k = n-1, 1, -1
        do i = 1, l
          a1(i,k) = a1(i,k)-au(i,k)*a1(i,k+1)
          a2(i,k) = a2(i,k)-au(i,k)*a2(i,k+1)
        end do
      end do
    end subroutine tridi2

    subroutine tridin(l,n,nt,cl,cm,cu,r1,r2,au,a1,a2)
      implicit none
      integer(ik4), intent(in) :: l, n, nt
      real(kind=rkx), dimension(l,2:n), intent(in):: cl
      real(kind=rkx), dimension(l,n), intent(in) :: cm
      real(kind=rkx), dimension(l,n-1), intent(in) ::  cu
      real(kind=rkx), dimension(l,n), intent(in) :: r1
      real(kind=rkx), dimension(l,n*nt), intent(in) :: r2
      real(kind=rkx), dimension(l,n-1), intent(inout) :: au
      real(kind=rkx), dimension(l,n), intent(inout) :: a1
      real(kind=rkx), dimension(l,n*nt), intent(inout) :: a2
      real(kind=rkx), dimension(l,2:n-1) :: fkk
      real(kind=rkx), dimension(l) :: fk
      integer(ik4) :: is, k, kk, i

      do i = 1, l
        fk(i)   = d_one/cm(i,1)
        au(i,1) = fk(i)*cu(i,1)
        a1(i,1) = fk(i)*r1(i,1)
      end do
      do k = 1, nt
        is = (k-1) * n
        do i = 1, l
          a2(i,1+is) = fk(i) * r2(i,1+is)
        end do
      end do
      do k = 2, n-1
        do i = 1, l
          fkk(i,k) = d_one/(cm(i,k)-cl(i,k)*au(i,k-1))
          au(i,k)  = fkk(i,k)*cu(i,k)
          a1(i,k)  = fkk(i,k)*(r1(i,k)-cl(i,k)*a1(i,k-1))
        end do
      end do
      do kk = 1, nt
        is = (kk-1) * n
        do k = 2, n-1
          do i = 1, l
            a2(i,k+is) = fkk(i,k)*(r2(i,k+is)-cl(i,k)*a2(i,k+is-1))
          end do
        end do
      end do
      do i = 1, l
        fk(i)   = d_one/(cm(i,n)-cl(i,n)*au(i,n-1))
        a1(i,n) = fk(i)*(r1(i,n)-cl(i,n)*a1(i,n-1))
      end do
      do k = 1, nt
        is = (k-1) * n
        do i = 1, l
          a2(i,n+is) = fk(i)*(r2(i,n+is)-cl(i,n)*a2(i,n+is-1))
        end do
      end do
      do k = n-1, 1, -1
        do i = 1, l
          a1(i,k) = a1(i,k) - au(i,k)*a1(i,k+1)
        end do
      end do
      do kk = 1, nt
        is = (kk-1) * n
        do k = n-1, 1, -1
          do i = 1, l
            a2(i,k+is) = a2(i,k+is) - au(i,k)*a2(i,k+is+1)
          end do
        end do
      end do
    end subroutine tridin

    subroutine tridit(l,n,nt,cl,cm,cu,rt,au,at)
      implicit none
      integer(ik4), intent(in) :: l, n, nt
      real(kind=rkx), dimension(l,2:n), intent(in) :: cl
      real(kind=rkx), dimension(l,n), intent(in) :: cm
      real(kind=rkx), dimension(l,n-1), intent(in) :: cu
      real(kind=rkx), dimension(l,n*nt), intent(in) :: rt
      real(kind=rkx), dimension(l,n-1), intent(inout) :: au
      real(kind=rkx), dimension(l,n*nt), intent(inout) :: at
      real(kind=rkx), dimension(l,2:n-1) :: fkk
      real(kind=rkx), dimension(l) :: fk
      integer(ik4) :: is, k, kk, i

      do i = 1, l
        fk(i)   = d_one/cm(i,1)
        au(i,1) = fk(i)*cu(i,1)
      end do
      do k = 1, nt
        is = (k-1) * n
        do i = 1, l
          at(i,1+is) = fk(i) * rt(i,1+is)
        end do
      end do
      do k = 2, n-1
        do i = 1, l
          fkk(i,k) = d_one/(cm(i,k)-cl(i,k)*au(i,k-1))
          au(i,k)  = fkk(i,k)*cu(i,k)
        end do
      end do
      do kk = 1, nt
        is = (kk-1) * n
        do k = 2, n-1
          do i = 1, l
            at(i,k+is) = fkk(i,k)*(rt(i,k+is)-cl(i,k)*at(i,k+is-1))
          end do
        end do
      end do
      do i = 1, l
        fk(i)   = d_one/(cm(i,n)-cl(i,n)*au(i,n-1))
      end do
      do k = 1, nt
        is = (k-1) * n
        do i = 1, l
          at(i,n+is) = fk(i)*(rt(i,n+is)-cl(i,n)*at(i,n+is-1))
        end do
      end do
      do kk = 1, nt
        is = (kk-1) * n
        do k = n-1, 1, -1
          do i = 1, l
            at(i,k+is) = at(i,k+is) - au(i,k)*at(i,k+is+1)
          end do
        end do
      end do
    end subroutine tridit

end module mod_pbl_gfs

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

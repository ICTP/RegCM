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

module mod_pbl_gfs

  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_memutil
  use mod_dynparam , only : kz , kzp1 , kzm1 , ntr
  use mod_dynparam , only : ici1 , ici2 , jci1 , jci2
  use mod_runparams , only : dt , dx , nqx , ichem , iqv
  use mod_regcm_types , only : mod_2_pbl , pbl_2_mod

  public :: blgfs2011

  integer(ik4) :: iblp , ibnt

  ! Input to moninq

  real(rkx) , dimension(:,:) , pointer :: uo , vo , t1
  real(rkx) , dimension(:,:,:) , pointer :: q1
  real(rkx) , dimension(:) , pointer :: psk , rbsoil
  real(rkx) , dimension(:) , pointer :: fm , fh , spd1
  real(rkx) , dimension(:,:) , pointer :: prsi , phii , z , dz
  real(rkx) , dimension(:,:) , pointer :: del , prsl , prslk , phil , thraten
  real(rkx) , dimension(:,:) , pointer :: rcs
  real(rkx) , dimension(:) , pointer :: heat , evap , stress

  ! Output from moninq

  real(rkx) , dimension(:) , pointer :: hpbl
  integer(ik4) , dimension(:) , pointer :: kpbl
  real(rkx) , dimension(:,:) , pointer:: dv , du
  real(rkx) , dimension(:,:,:) , pointer :: rtg
  real(rkx) , dimension(:,:) , pointer :: tau

  integer(ik4) , parameter :: npsi = 1000

  real(rkx) , dimension(:) , pointer :: psim_stab , psim_unstab
  real(rkx) , dimension(:) , pointer  :: psih_stab , psih_unstab

  contains

    subroutine init_pbl_gfs
      implicit none
      integer(ik4) :: i
      real(rkx) :: zolf

      iblp = (jci2-jci1+1)*(ici2-ici1+1)
      ibnt = nqx+ntr

      call getmem1d(hpbl,1,iblp,'mod_pbl_gfs:hpbl')
      call getmem1d(kpbl,1,iblp,'mod_pbl_gfs:kpbl')
      call getmem2d(du,1,iblp,1,kz,'mod_pbl_gfs:du')
      call getmem2d(dv,1,iblp,1,kz,'mod_pbl_gfs:dv')
      call getmem2d(tau,1,iblp,1,kz,'mod_pbl_gfs:tau')
      call getmem3d(rtg,1,iblp,1,kz,1,ibnt,'mod_pbl_gfs:rtg')

      call getmem2d(rcs,1,iblp,1,kz,'mod_pbl_gfs:rcs')
      call getmem2d(uo,1,iblp,1,kz,'mod_pbl_gfs:uo')
      call getmem2d(vo,1,iblp,1,kz,'mod_pbl_gfs:vo')
      call getmem2d(t1,1,iblp,1,kz,'mod_pbl_gfs:t1')
      call getmem3d(q1,1,iblp,1,kz,1,ibnt,'mod_pbl_gfs:q1')
      call getmem1d(psk,1,iblp,'mod_pbl_gfs:psk')
      call getmem1d(rbsoil,1,iblp,'mod_pbl_gfs:rbsoil')
      call getmem1d(fm,1,iblp,'mod_pbl_gfs:fm')
      call getmem1d(fh,1,iblp,'mod_pbl_gfs:fh')
      call getmem1d(spd1,1,iblp,'mod_pbl_gfs:spd1')
      call getmem2d(prsi,1,iblp,1,kzp1,'mod_pbl_gfs:prsi')
      call getmem2d(phii,1,iblp,1,kzp1,'mod_pbl_gfs:phii')
      call getmem2d(del,1,iblp,1,kz,'mod_pbl_gfs:del')
      call getmem2d(prsl,1,iblp,1,kz,'mod_pbl_gfs:prsl')
      call getmem2d(prslk,1,iblp,1,kz,'mod_pbl_gfs:prslk')
      call getmem2d(phil,1,iblp,1,kz,'mod_pbl_gfs:phil')
      call getmem2d(thraten,1,iblp,1,kz,'mod_pbl_gfs:thraten')
      call getmem2d(z,1,iblp,1,kz,'mod_pbl_gfs:z')
      call getmem2d(dz,1,iblp,1,kz,'mod_pbl_gfs:dz')
      call getmem1d(heat,1,iblp,'mod_pbl_gfs:heat')
      call getmem1d(evap,1,iblp,'mod_pbl_gfs:evap')
      call getmem1d(stress,1,iblp,'mod_pbl_gfs:stress')

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

    end subroutine init_pbl_gfs

    subroutine pbl_gfs(m2p,p2m)
      implicit none
      type(mod_2_pbl) , intent(in) :: m2p
      type(pbl_2_mod) , intent(inout) :: p2m

      integer(ik4) :: i , j , k , kk , km , n
      integer(ik4) :: iq , it , iit
      real(rkx) :: tvcon , zo , gz1oz0 , thvx , tskv , dthvdz , br
      real(rkx) :: ps , wspd , psim , psih , dthvm , fluxc , tsk
      real(rkx) :: za , ta , qa , pa , ua , va , tha , rhoa
      real(rkx) :: rrhox , cpm , vconv , vsgd , hf , qf , xp
      real(rkx) :: wspd0 , zol , zol0 , zolzz
      real(rkx) :: tv1 , tv2 , u1 , u2

      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          za = m2p%za(j,i,kz)
          ta = m2p%tatm(j,i,kz)
          qa = m2p%qxatm(j,i,kz,iqv)
          ps = m2p%patmf(j,i,kzp1)
          pa = m2p%patm(j,i,kz)
          ua = m2p%uxatm(j,i,kz)
          va = m2p%vxatm(j,i,kz)
          tha = m2p%thatm(j,i,kz)
          rhoa = m2p%rhox2d(j,i)
          hf = m2p%hfx(j,i)
          qf = m2p%qfx(j,i)
          tsk = m2p%tsk(j,i)
          tvcon = d_one + ep1*qa
          rrhox = (rgas*(ta*tvcon))/pa
          cpm = cpd * (d_one + 0.8_rkx * qa)
          zo = min(max(m2p%zo(j,i),0.01_rkx),za)
          xp =  (p00/ps)**rovcp
          gz1oz0 = log((za+zo)/zo)
          thvx = tha*tvcon
          tskv = tsk*xp*tvcon
          dthvdz = thvx - tskv
          wspd0 = sqrt(ua*ua + va*va)
          vsgd = 0.32_rkx * (max(dx/5000.0_rkx-d_one,d_zero))**0.33_rkx
          if ( m2p%ldmsk(j,i) > 0 ) then
            fluxc = max(hf/rhoa*rcpd + qf/rhoa*ep1*tskv,0.0_rkx)
            vconv = d_one*(egrav/tsk*hpbl(n)*fluxc)**0.33_rkx
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
            if ( m2p%uvdrag(j,i) < 0.001_rkx ) then
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
          fm(n) = gz1oz0 - psim
          fh(n) = gz1oz0 - psih
          psk(n) = (ps/p00)**rovcp
          stress(n) = vonkar*vonkar*wspd0*wspd0/(fm(n)*fm(n))
          heat(n) = hf/cpm*rrhox
          evap(n) = qf*rrhox
          spd1(n) = wspd0
          rbsoil(n) = br
          hpbl(n) = p2m%zpbl(j,i)
          prsi(n,1) = ps*d_r1000
          phii(n,1) = d_zero
          n = n + 1
        end do
      end do

      do k = 1 , kz
        n = 1
        kk = kzp1 - k
        do i = ici1 , ici2
          do j = jci1 , jci2
            t1(n,kk) = m2p%tatm(j,i,k)
            uo(n,kk) = m2p%uxatm(j,i,k)
            vo(n,kk) = m2p%vxatm(j,i,k)
            prsl(n,kk) = m2p%patm(j,i,k)*d_r1000
            prslk(n,kk) = (m2p%patm(j,i,k)/p00)**rovcp
            dz(n,kk) = m2p%dzq(j,i,k)
            z(n,kk) = m2p%za(j,i,k)
            thraten(n,kk) = m2p%heatrt(j,i,k)/prslk(n,kk)
            rcs(n,kk) = 1.0_rkx
            n = n + 1
          end do
        end do
      end do

      do k = 2 , kz
        n = 1
        km = k - 1
        do i = ici1 , ici2
          do j = jci1 , jci2
            del(n,km) = prsl(n,km)/rovg*dz(n,km)/t1(n,km)
            prsi(n,k) = prsi(n,km)-del(n,km)
            phii(n,k) = (z(n,k)-z(n,1))*egrav
            phil(n,km) = d_half*(z(n,k)+z(n,km)-d_two*z(n,1))*egrav
            n = n + 1
          end do
        end do
      end do

      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          del(n,kz) = del(n,kzm1)
          prsi(n,kzp1) = prsi(n,kz)-del(n,kzm1)
          phii(n,kzp1) = phii(n,kz)+dz(n,kz)*egrav
          phil(n,kz) = phii(n,kz)-phil(n,kzm1)+phii(n,kz)
          n = n + 1
        end do
      end do

      do iq = 1 , nqx
        do k = 1 , kz
          n = 1
          kk = kzp1 - k
          do i = ici1 , ici2
            do j = jci1 , jci2
              q1(n,kk,iq) = m2p%qxatm(j,i,k,iq)/(d_one + m2p%qxatm(j,i,k,iq))
              n = n + 1
            end do
          end do
        end do
      end do

      if ( ichem == 1 ) then
        do it = 1 , ntr
          iit = it + nqx
          do k = 1 , kz
            n = 1
            kk = kzp1 - k
            do i = ici1 , ici2
              do j = jci1 , jci2
                q1(n,kk,iit) = m2p%chib(j,i,k,it)/(d_one + m2p%chib(j,i,k,it))
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

      do k = 1 , kz
        n = 1
        kk = kzp1 - k
        do i = ici1 , ici2
          do j = jci1 , jci2
            p2m%uxten(j,i,k) = du(n,kk)
            p2m%vxten(j,i,k) = dv(n,kk)
            p2m%tten(j,i,k) = p2m%tten(j,i,k) + tau(n,kk)*m2p%psb(j,i)
            n = n + 1
          end do
        end do
      end do

      do iq = 1 , nqx
        do k = 1 , kz
          n = 1
          kk = kzp1 - k
          do i = ici1 , ici2
            do j = jci1 , jci2
              p2m%qxten(j,i,k,iq) = p2m%qxten(j,i,k,iq) + &
                     rtg(n,kk,iq)/(d_one-q1(n,kk,iq))**2 * m2p%psb(j,i)
              n = n + 1
            end do
          end do
        end do
      end do

      if ( ichem == 1 ) then
        do it = 1 , ntr
          iit = it + nqx
          do k = 1 , kz
            n = 1
            kk = kzp1 - k
            do i = ici1 , ici2
              do j = jci1 , jci2
                p2m%chiten(j,i,k,it) = p2m%chiten(j,i,k,it) + &
                      rtg(n,kk,iit)/(d_one-q1(n,kk,iit)) * m2p%psb(j,i)
                n = n + 1
              end do
            end do
          end do
        end do
      end if

      n = 1
      do i = ici1 , ici2
        do j = jci1 , jci2
          p2m%kpbl(j,i) = kzp1 - kpbl(n)
          p2m%zpbl(j,i) = hpbl(n)
          n = n + 1
        end do
      end do

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

    end subroutine pbl_gfs

    subroutine moninq(im,km,ntrac)
      implicit none

      integer(ik4) , intent(in) :: im , km , ntrac
      integer(ik4) :: i , is , k , kk , km1 , kmpbl
      integer(ik4) , dimension(im) :: lcld , icld , kcld , krad
      integer(ik4) , dimension(im) :: kpblx

      real(rkx) , dimension(im) :: phih , phim
      real(rkx) , dimension(im) :: rbdn , rbup , beta
      real(rkx) , dimension(im) :: ustar , wscale , thermal , wstar3
      real(rkx) , dimension(im) :: hpblx

      real(rkx) , dimension(im,km) :: thvx , thlvx
      real(rkx) , dimension(im,km) :: qlx , thetae , qtx
      real(rkx) , dimension(im,km-1) :: bf , radx
      real(rkx) , dimension(im,km) :: u1 , v1
      real(rkx) , dimension(im) :: govrth , hrad , cteit
      real(rkx) , dimension(im) :: radmin , vrad , zd , zdd , thlvx1
      real(rkx) , dimension(im) :: hgamt , hgamq

      real(rkx) , dimension(im,km-1) :: rdzt , dktx , dkux , xkzo , xkzmo
      real(rkx) , dimension(im,km+1) :: zi
      real(rkx) , dimension(im,km-1) :: dku , dkt , cku , ckt , al , au
      real(rkx) , dimension(im,km) :: zl , ad , a1 , theta
      real(rkx) , dimension(im,km*ntrac) :: a2
      real(rkx) , dimension(im) :: prinv , rent
      logical , dimension(im) :: pblflg , sfcflg , scuflg , flg
      real(rkx) :: bvf2 , dk , dsdz2 , dsdzq , dsdzt
      real(rkx) :: dsig , dtodsd , dtodsu , dw2 , hol , hol1
      real(rkx) :: prnum , qtend , rbint , rdt , rdz , ri , rl2
      real(rkx) :: sflux , shr2 , spdk2 , sri , tem , ti
      real(rkx) :: ttend , utend , vtend , zfac , vpert , zk
      real(rkx) :: tem1 , tem2 , ptem , ptem1 , ptem2

      real(rkx) , parameter :: gocp = egrav*rcpd
      real(rkx) , parameter :: rlam = 30.0_rkx
      real(rkx) , parameter :: vk = vonkar
      real(rkx) , parameter :: prmin = 0.25_rkx
      real(rkx) , parameter :: prmax = 4.0_rkx
      real(rkx) , parameter :: dw2min = 0.0001_rkx
      real(rkx) , parameter :: dkmin = 0.0_rkx
      real(rkx) , parameter :: dkmax = 1000.0_rkx
      real(rkx) , parameter :: rimin = -100.0_rkx
      real(rkx) , parameter :: rbcr = 0.25_rkx
      real(rkx) , parameter :: wfac = 7.0_rkx
      real(rkx) , parameter :: cfac = 6.5_rkx
      real(rkx) , parameter :: pfac = 2.0_rkx
      real(rkx) , parameter :: sfcfrac = 0.1_rkx
      real(rkx) , parameter :: qmin = 1.e-8_rkx
      real(rkx) , parameter :: xkzm = 1.0_rkx
      real(rkx) , parameter :: zfmin = 1.e-8_rkx
      real(rkx) , parameter :: aphi5 = 5.0_rkx
      real(rkx) , parameter :: aphi16 = 16.0_rkx
      real(rkx) , parameter :: tdzmin = 1.e-3_rkx
      real(rkx) , parameter :: qlmin = 1.e-12_rkx
      real(rkx) , parameter :: h1 = 0.33333333_rkx
      real(rkx) , parameter :: cldtime = 500.0_rkx
      real(rkx) , parameter :: xkzmu = 3.0_rkx
      real(rkx) , parameter :: xkzminv = 0.3_rkx
      real(rkx) , parameter :: gamcrt = 3.0_rkx
      real(rkx) , parameter :: gamcrq = 0.0_rkx
      real(rkx) , parameter :: rlamun = 150.0_rkx
      real(rkx) , parameter :: rentf1 = 0.2_rkx
      real(rkx) , parameter :: rentf2 = 1.0_rkx
      real(rkx) , parameter :: radfac = 0.85_rkx
      real(rkx) , parameter :: zstblmax = 2500.0_rkx
      real(rkx) , parameter :: qlcr=3.5e-5_rkx
      real(rkx) , parameter :: actei = 0.7_rkx

      rdt   = d_one / dt
      km1   = km - 1
      kmpbl = km / 2
      do k = 1 , km
        do i = 1 , im
          zi(i,k) = phii(i,k) * regrav
          zl(i,k) = phil(i,k) * regrav
          u1(i,k) = uo(i,k) * rcs(i,k)
          v1(i,k) = vo(i,k) * rcs(i,k)
        end do
      end do
      do i = 1 , im
        zi(i,km+1) = phii(i,km+1) * regrav
      end do
      do k = 1 , km1
        do i = 1 , im
          rdzt(i,k) = d_one / (zl(i,k+1) - zl(i,k))
        end do
      end do
      ! Vertical background diffusivity
      do k = 1 , km1
        do i = 1 , im
          tem1      = d_one - prsi(i,k+1) / prsi(i,1)
          tem1      = min(tem1 * tem1 * 10.0_rkx,25.0_rkx)
          xkzo(i,k) = xkzm * min(d_one, exp(-tem1))
        end do
      end do
      ! Vertical background diffusivity for momentum
      do k = 1 , km1
        do i = 1 , im
          ptem = prsi(i,k+1) / prsi(i,1)
          if ( ptem >= 0.2_rkx ) then
            xkzmo(i,k) = xkzmu
            ptem1 = prsi(i,k+1)
          else
            tem1 = d_one - prsi(i,k+1) / ptem1
            tem1 = min(tem1 * tem1 * 5.0_rkx,25.0_rkx)
            xkzmo(i,k) = xkzmu * min(d_one, exp(-tem1))
          end if
        end do
      end do
      ! Diffusivity in the inversion layer is set to be xkzminv (m^2/s)
      do k = 1 , kmpbl
        do i = 1 , im
          if ( zi(i,k+1) > 250.0_rkx ) then
            tem1 = (t1(i,k+1)-t1(i,k)) * rdzt(i,k)
            if ( tem1 > 1.e-5_rkx ) then
              xkzo(i,k) = min(xkzo(i,k),xkzminv)
            end if
          end if
        end do
      end do
      do i = 1 , im
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
      do k = 1 , km
        do i = 1 , im
          theta(i,k) = t1(i,k) * psk(i) / prslk(i,k)
          qlx(i,k)   = max(q1(i,k,2),qlmin)
          qtx(i,k)   = max(q1(i,k,1),qmin)+qlx(i,k)
          ptem       = qlx(i,k)
          ptem1      = wlhv*max(q1(i,k,1),qmin)/(cpd*t1(i,k))
          thetae(i,k)= theta(i,k)*(d_one+ptem1)
          thvx(i,k)  = theta(i,k)*(d_one+ep1*max(q1(i,k,1),qmin)-ptem)
          ptem2      = theta(i,k)-(wlhv/cpd)*ptem
          thlvx(i,k) = ptem2*(d_one+ep1*qtx(i,k))
        end do
      end do
      do k = 1 , km1
        do i = 1 , im
          dku(i,k)  = d_zero
          dkt(i,k)  = d_zero
          dktx(i,k) = d_zero
          dkux(i,k) = d_zero
          cku(i,k)  = d_zero
          ckt(i,k)  = d_zero
          tem       = zi(i,k+1)-zi(i,k)
          radx(i,k) = tem*thraten(i,k)
        end do
      end do
      do i = 1 , im
        flg(i) = scuflg(i)
      end do
      do k = 1 , km1
        do i = 1 , im
          if ( flg(i) .and. zl(i,k) >= zstblmax ) then
            lcld(i) = k
            flg(i) = .false.
          end if
        end do
      end do
      ! Compute buoyancy flux
      do k = 1 , km1
        do i = 1 , im
          bf(i,k) = (thvx(i,k+1)-thvx(i,k))*rdzt(i,k)
        end do
      end do
      do i = 1 , im
        govrth(i) = egrav/theta(i,1)
      end do
      do i = 1 , im
        beta(i) = dt / (zi(i,2)-zi(i,1))
      end do
      do i = 1 , im
        ustar(i) = sqrt(stress(i))
        thermal(i) = thvx(i,1)
      end do
      ! Compute the first guess pbl height
      do i = 1 , im
        flg(i) = .false.
        rbup(i) = rbsoil(i)
      end do
      do k = 2 , kmpbl
        do i = 1 , im
          if ( .not. flg(i) ) then
            rbdn(i) = rbup(i)
            spdk2   = max((u1(i,k)**2+v1(i,k)**2),d_one)
            rbup(i) = (thvx(i,k)-thermal(i))*        &
                      (egrav*zl(i,k)/thvx(i,1))/spdk2
            kpbl(i) = k
            flg(i)  = (rbup(i) > rbcr)
          end if
        end do
      end do
      do i = 1 , im
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
      do i = 1 , im
        sflux = heat(i) + evap(i)*ep1*theta(i,1)
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
      do i = 1 , im
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
      do i = 1 , im
        flg(i)  = .true.
        if ( pblflg(i) ) then
          flg(i)  = .false.
          rbup(i) = rbsoil(i)
        end if
      end do
      do k = 2 , kmpbl
        do i = 1 , im
          if ( .not. flg(i) ) then
            rbdn(i) = rbup(i)
            spdk2   = max((u1(i,k)**2+v1(i,k)**2),d_one)
            rbup(i) = (thvx(i,k)-thermal(i)) * &
                      (egrav*zl(i,k)/thvx(i,1))/spdk2
            kpbl(i) = k
            flg(i)  = (rbup(i) > rbcr)
          end if
        end do
      end do
      do i = 1 , im
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
      do i = 1 , im
        flg(i) = scuflg(i)
      end do
      do k = kmpbl , 1 , -1
        do i = 1 , im
          if ( flg(i) .and. k <= lcld(i) ) then
            if ( qlx(i,k) >= qlcr ) then
              kcld(i) = k
              flg(i) = .false.
            end if
          end if
        end do
      end do
      do i = 1 , im
        if ( scuflg(i) .and. kcld(i) == km1 ) scuflg(i) = .false.
      end do
      do i = 1 , im
        flg(i) = scuflg(i)
      end do
      do k = kmpbl , 1 , -1
        do i = 1 , im
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
      do i = 1 , im
        if ( scuflg(i) .and. krad(i) <= 1 ) scuflg(i) = .false.
        if ( scuflg(i) .and. radmin(i) >= d_zero ) scuflg(i) = .false.
      end do
      do i = 1 , im
        flg(i) = scuflg(i)
      end do
      do k = kmpbl , 2 , -1
        do i = 1 , im
          if ( flg(i) .and. k <= krad(i) ) then
            if ( qlx(i,k) >= qlcr ) then
              icld(i) = icld(i)+1
            else
              flg(i) = .false.
            end if
          end if
        end do
      end do
      do i = 1 , im
        if ( scuflg(i) .and. icld(i) < 1 ) scuflg(i) = .false.
      enddo
      do i = 1 , im
        if ( scuflg(i) ) then
          hrad(i) = zi(i,krad(i)+1)
        end if
      end do
      do i = 1 , im
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
      do i = 1 , im
        flg(i) = scuflg(i)
      end do
      do k = kmpbl , 1 , -1
        do i = 1 , im
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
      do i = 1 , im
        if ( scuflg(i) ) then
          zd(i) = max(zd(i),zdd(i))
          zd(i) = min(zd(i),hrad(i))
          tem   = govrth(i)*zd(i)*(-radmin(i))
          vrad(i) = tem**h1
        end if
      end do
      ! Compute inverse Prandtl number
      do i = 1 , im
        if ( pblflg(i) ) then
          tem = phih(i)/phim(i)+cfac*vk*sfcfrac
          prinv(i) = d_one / tem
          prinv(i) = min(prinv(i),prmax)
          prinv(i) = max(prinv(i),prmin)
        end if
      end do
      ! Compute diffusion coefficients below pbl
      do k = 1 , kmpbl
        do i = 1 , im
          if ( pblflg(i) .and. k < kpbl(i) ) then
            zfac = max((d_one-zi(i,k+1)/(hpbl(i))), zfmin)
            tem = wscale(i)*vk*zi(i,k+1)*zfac**pfac
            dku(i,k) = xkzmo(i,k) + tem
            dkt(i,k) = xkzo(i,k) + tem * prinv(i)
            dku(i,k) = min(dku(i,k),dkmax)
            dkt(i,k) = min(dkt(i,k),dkmax)
            dktx(i,k)= dkt(i,k)
            dkux(i,k)= dku(i,k)
          end if
        end do
      end do
      ! Compute diffusion coefficients based on local scheme
      do i = 1 , im
        if ( .not. pblflg(i) ) then
          kpbl(i) = 1
        end if
      end do
      do k = 1 , km1
        do i = 1 , im
          if ( k >= kpbl(i) ) then
            rdz  = rdzt(i,k)
            ti   = d_two/(t1(i,k)+t1(i,k+1))
            dw2  = (u1(i,k)-u1(i,k+1))**2 +(v1(i,k)-v1(i,k+1))**2
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
      do i = 1 , im
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
      do k = 1 , kmpbl
        do i = 1 , im
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
      do k = 1 , kmpbl
        do i = 1 , im
          if ( scuflg(i) ) then
            dkt(i,k) = dkt(i,k)+ckt(i,k)
            dku(i,k) = dku(i,k)+cku(i,k)
            dkt(i,k) = min(dkt(i,k),dkmax)
            dku(i,k) = min(dku(i,k),dkmax)
          end if
        end do
      end do
      ! Compute tridiagonal matrix elements for heat and moisture
      do i = 1 , im
         ad(i,1) = d_one
         a1(i,1) = t1(i,1)   + beta(i) * heat(i)
         a2(i,1) = q1(i,1,1) + beta(i) * evap(i)
      end do
      if ( ntrac >= 2 ) then
        do k = 2 , ntrac
          is = (k-1) * km
          do i = 1 , im
            a2(i,1+is) = q1(i,1,k)
          end do
        end do
      end if
      do k = 1 , km1
        do i = 1 , im
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
             a2(i,k+1) = q1(i,k+1,1)-dtodsu*dsdzq
          else
             dsdzt = tem1 * gocp
             a2(i,k+1) = q1(i,k+1,1)
          endif
          dsdz2     = tem1 * rdz
          au(i,k)   = -dtodsd*dsdz2
          al(i,k)   = -dtodsu*dsdz2
          ad(i,k)   = ad(i,k)-au(i,k)
          ad(i,k+1) = d_one-al(i,k)
          a1(i,k)   = a1(i,k)+dtodsd*dsdzt
          a1(i,k+1) = t1(i,k+1)-dtodsu*dsdzt
        end do
      end do
      if ( ntrac >= 2 ) then
        do kk = 2 , ntrac
          is = (kk-1) * km
          do k = 1 , km1
            do i = 1 , im
              a2(i,k+1+is) = q1(i,k+1,kk)
            end do
          end do
        end do
      end if
      ! Solve tridiagonal problem for heat and moisture
      call tridin(im,km,ntrac,al,ad,au,a1,a2,au,a1,a2)
      ! Recover tendencies of heat and moisture
      do k = 1,km
        do i = 1 , im
          ttend      = (a1(i,k)-t1(i,k))*rdt
          qtend      = (a2(i,k)-q1(i,k,1))*rdt
          tau(i,k)   = tau(i,k)+ttend
          rtg(i,k,1) = rtg(i,k,1)+qtend
        end do
      end do
      if ( ntrac >= 2 ) then
        do kk = 2 , ntrac
          is = (kk-1) * km
          do k = 1 , km
            do i = 1 , im
              qtend = (a2(i,k+is)-q1(i,k,kk))*rdt
              rtg(i,k,kk) = rtg(i,k,kk)+qtend
            end do
          end do
        end do
      end if
      ! Compute tridiagonal matrix elements for momentum
      do i = 1 , im
        ad(i,1) = d_one + beta(i) * stress(i) / spd1(i)
        a1(i,1) = u1(i,1)
        a2(i,1) = v1(i,1)
      end do
      do k = 1 , km1
        do i = 1 , im
          dtodsd = dt/del(i,k)
          dtodsu = dt/del(i,k+1)
          dsig   = prsl(i,k)-prsl(i,k+1)
          rdz    = rdzt(i,k)
          tem1   = dsig*dku(i,k)*rdz
          a1(i,k+1) = u1(i,k+1)
          a2(i,k+1) = v1(i,k+1)
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
      do k = 1 , km
        do i = 1 , im
          ptem = d_one/rcs(i,k)
          utend = (a1(i,k)-u1(i,k))*rdt
          vtend = (a2(i,k)-v1(i,k))*rdt
          du(i,k)  = du(i,k)+utend*ptem
          dv(i,k)  = dv(i,k)+vtend*ptem
        end do
      end do
      ! PBL height for diagnostic purpose
      do i = 1 , im
        hpbl(i) = hpblx(i)
        kpbl(i) = kpblx(i)
      end do
    end subroutine moninq

    subroutine tridi2(l,n,cl,cm,cu,r1,r2,au,a1,a2)
      implicit none
      integer(ik4) , intent(in) :: l , n
      real(kind=rkx) , dimension(l,2:n) , intent(in) :: CL
      real(kind=rkx) , dimension(l,n) , intent(in) :: CM
      real(kind=rkx) , dimension(l,n-1) , intent(in) :: CU
      real(kind=rkx) , dimension(l,n) , intent(in) :: R1
      real(kind=rkx) , dimension(l,n) , intent(in) :: R2
      real(kind=rkx) , dimension(l,n-1) , intent(inout) :: au
      real(kind=rkx) , dimension(l,n) , intent(inout) :: a1
      real(kind=rkx) , dimension(l,n) , intent(inout) :: a2
      real(kind=rkx) :: fk
      integer(ik4) :: k , i

      do i = 1 , l
        fk      = d_one/cm(i,1)
        au(i,1) = fk*cu(i,1)
        a1(i,1) = fk*r1(i,1)
        a2(i,1) = fk*r2(i,1)
      end do
      do k = 2 , n-1
        do i = 1 , l
          fk      = d_one/(cm(i,k)-cl(i,k)*au(i,k-1))
          au(i,k) = fk*cu(i,k)
          a1(i,k) = fk*(r1(i,k)-cl(i,k)*a1(i,k-1))
          a2(i,k) = fk*(r2(i,k)-cl(i,k)*a2(i,k-1))
        end do
      end do
      do i = 1 , l
        fk      = d_one/(cm(i,n)-cl(i,n)*au(i,n-1))
        a1(i,n) = fk*(r1(i,n)-cl(i,n)*a1(i,n-1))
        a2(i,n) = fk*(r2(i,n)-cl(i,n)*a2(i,n-1))
      end do
      do k = n-1 , 1 , -1
        do i = 1 , l
          a1(i,k) = a1(i,k)-au(i,k)*a1(i,k+1)
          a2(i,k) = a2(i,k)-au(i,k)*a2(i,k+1)
        end do
      end do
    end subroutine tridi2

    subroutine tridin(l,n,nt,cl,cm,cu,r1,r2,au,a1,a2)
      implicit none
      integer(ik4) , intent(in) :: l , n , nt
      real(kind=rkx) , dimension(l,2:n) , intent(in):: cl
      real(kind=rkx) , dimension(l,n) , intent(in) :: cm
      real(kind=rkx) , dimension(l,n-1) , intent(in) ::  cu
      real(kind=rkx) , dimension(l,n) , intent(in) :: r1
      real(kind=rkx) , dimension(l,n*nt) , intent(in) :: r2
      real(kind=rkx) , dimension(l,n-1) , intent(inout) :: au
      real(kind=rkx) , dimension(l,n) , intent(inout) :: a1
      real(kind=rkx) , dimension(l,n*nt) , intent(inout) :: a2
      real(kind=rkx) , dimension(l,2:n-1) :: fkk
      real(kind=rkx) , dimension(l) :: fk
      integer(ik4) :: is , k , kk , i

      do i = 1 , l
        fk(i)   = d_one/cm(i,1)
        au(i,1) = fk(i)*cu(i,1)
        a1(i,1) = fk(i)*r1(i,1)
      end do
      do k = 1 , nt
        is = (k-1) * n
        do i = 1 , l
          a2(i,1+is) = fk(i) * r2(i,1+is)
        end do
      end do
      do k = 2 , n-1
        do i = 1 , l
          fkk(i,k) = d_one/(cm(i,k)-cl(i,k)*au(i,k-1))
          au(i,k)  = fkk(i,k)*cu(i,k)
          a1(i,k)  = fkk(i,k)*(r1(i,k)-cl(i,k)*a1(i,k-1))
        end do
      end do
      do kk = 1 , nt
        is = (kk-1) * n
        do k = 2 , n-1
          do i = 1 , l
            a2(i,k+is) = fkk(i,k)*(r2(i,k+is)-cl(i,k)*a2(i,k+is-1))
          end do
        end do
      end do
      do i = 1 , l
        fk(i)   = d_one/(cm(i,n)-cl(i,n)*au(i,n-1))
        a1(i,n) = fk(i)*(r1(i,n)-cl(i,n)*a1(i,n-1))
      end do
      do k = 1 , nt
        is = (k-1) * n
        do i = 1 , l
          a2(i,n+is) = fk(i)*(r2(i,n+is)-cl(i,n)*a2(i,n+is-1))
        end do
      end do
      do k = n-1 , 1 , -1
        do i = 1 , l
          a1(i,k) = a1(i,k) - au(i,k)*a1(i,k+1)
        end do
      end do
      do kk = 1 , nt
        is = (kk-1) * n
        do k = n-1 , 1 , -1
          do i = 1 , l
            a2(i,k+is) = a2(i,k+is) - au(i,k)*a2(i,k+is+1)
          end do
        end do
      end do
    end subroutine tridin

    subroutine tridit(l,n,nt,cl,cm,cu,rt,au,at)
      implicit none
      integer(ik4) , intent(in) :: l , n , nt
      real(kind=rkx) , dimension(l,2:n) , intent(in) :: cl
      real(kind=rkx) , dimension(l,n) , intent(in) :: cm
      real(kind=rkx) , dimension(l,n-1) , intent(in) :: cu
      real(kind=rkx) , dimension(l,n*nt) , intent(in) :: rt
      real(kind=rkx) , dimension(l,n-1) , intent(inout) :: au
      real(kind=rkx) , dimension(l,n*nt) , intent(inout) :: at
      real(kind=rkx) , dimension(l,2:n-1) :: fkk
      real(kind=rkx) , dimension(l) :: fk
      integer(ik4) :: is , k , kk , i

      do i = 1 , l
        fk(i)   = d_one/cm(i,1)
        au(i,1) = fk(i)*cu(i,1)
      end do
      do k = 1 , nt
        is = (k-1) * n
        do i = 1 , l
          at(i,1+is) = fk(i) * rt(i,1+is)
        end do
      end do
      do k = 2 , n-1
        do i = 1 , l
          fkk(i,k) = d_one/(cm(i,k)-cl(i,k)*au(i,k-1))
          au(i,k)  = fkk(i,k)*cu(i,k)
        end do
      end do
      do kk = 1 , nt
        is = (kk-1) * n
        do k = 2 , n-1
          do i = 1 , l
            at(i,k+is) = fkk(i,k)*(rt(i,k+is)-cl(i,k)*at(i,k+is-1))
          end do
        end do
      end do
      do i = 1 , l
        fk(i)   = d_one/(cm(i,n)-cl(i,n)*au(i,n-1))
      end do
      do k = 1 , nt
        is = (k-1) * n
        do i = 1 , l
          at(i,n+is) = fk(i)*(rt(i,n+is)-cl(i,n)*at(i,n+is-1))
        end do
      end do
      do kk = 1 , nt
        is = (kk-1) * n
        do k = n-1 , 1 , -1
          do i = 1 , l
            at(i,k+is) = at(i,k+is) - au(i,k)*at(i,k+is+1)
          end do
        end do
      end do
    end subroutine tridit

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

end module mod_pbl_gfs

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

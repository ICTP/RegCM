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

module mod_pbl_myj
  !
  !
  ! REFERENCES:  Janjic (1994), Mon. Wea. Rev., 122, 927-945.
  !              Janjic (2002), NCEP Office Note 437.
  !              Mellor and Yamada (1982), Rev. Geophys. Space Phys.
  !
  ! ABSTRACT:
  !     MYJ updates the turbulent kinetic energy with the production/
  !     dissipation term and the vertical diffusion term
  !     (using an implicit formulation) from Mellor-Yamada
  !     level 2.5 as extended by Janjic.  Exchange coefficients for
  !     the surface and for all layer interfaces are computed from
  !     Monin-Obukhov theory via the companion module mod_pbl_myj_sflayer,
  !     which provides a self-consistent Janjic Eta surface layer.
  !     The turbulent vertical exchange is then executed.
  !
  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_memutil
  use mod_regcm_types
  use mod_runparams
  use mod_pbl_myj_sflayer, only : myj_sfclayer

  implicit none

  private

  public :: init_myjpbl, myjpbl, myjtkemin

  ! Half of epsq2, the minimum q**2
  real(rkx), parameter :: myjtkemin = 0.1_rkx

  real(rkx), parameter :: elocp = eliwv/cpd
  real(rkx), parameter :: epsq2 = 0.2_rkx
  real(rkx), parameter :: eps1 = 1.e-12_rkx
  real(rkx), parameter :: eps2 = 0.0_rkx
  real(rkx), parameter :: epsl = 0.32_rkx
  real(rkx), parameter :: epsru = 1.e-7_rkx
  real(rkx), parameter :: epsrs = 1.e-7_rkx
  real(rkx), parameter :: epstrb = 1.e-24_rkx
  real(rkx), parameter :: epsa = 1.e-8_rkx
  real(rkx), parameter :: epsit = 1.e-4_rkx
  real(rkx), parameter :: fh = 1.01_rkx
  real(rkx), parameter :: alph = 0.30_rkx
  real(rkx), parameter :: beta = d_one/273.0_rkx
  real(rkx), parameter :: el0max = 1000.0_rkx
  real(rkx), parameter :: el0min = 1.0_rkx
  real(rkx), parameter :: elfc = 0.23_rkx*0.5_rkx
  real(rkx), parameter :: gam1 = 0.2222222222222222222_rkx
  real(rkx), parameter :: prt = 1.0_rkx
  real(rkx), parameter :: a1 = 0.659888514560862645_rkx
  real(rkx), parameter :: axx = 0.6574209922667784586_rkx
  real(rkx), parameter :: b1 = 11.87799326209552761_rkx
  real(rkx), parameter :: b2 = 7.226971804046074028_rkx
  real(rkx), parameter :: c1 = 0.000830955950095854396_rkx
  real(rkx), parameter :: elz0 = 0.0_rkx
  real(rkx), parameter :: esq = 5.0_rkx
  real(rkx), parameter :: btg = beta*egrav
  real(rkx), parameter :: esqhf = 0.5_rkx*5.0_rkx
  real(rkx), parameter :: rb1 = d_one/b1
  real(rkx), parameter :: adnh = 9.0_rkx*a1*axx*axx * &
               (12.0_rkx*a1+3.0_rkx*b2)*btg*btg
  real(rkx), parameter :: adnm = 18.0_rkx*a1*a1*axx*(b2-3.0_rkx*axx)*btg
  real(rkx), parameter :: anmh = -9.0_rkx*a1*axx*axx*btg*btg
  real(rkx), parameter :: anmm = -3.0_rkx*a1*axx * &
               (3.0_rkx*axx+3.0_rkx*b2*c1+18.0_rkx*a1*c1-b2)*btg
  real(rkx), parameter :: bdnh = 3.0_rkx*axx*(7.0_rkx*a1+b2)*btg
  real(rkx), parameter :: bdnm = 6.0_rkx*a1*a1
  real(rkx), parameter :: beqh = axx*b1*btg+3.0_rkx*axx*(7.0_rkx*a1+b2)*btg
  real(rkx), parameter :: beqm = -a1*b1*(1.0_rkx-3.0_rkx*c1)+6.0_rkx*a1*a1
  real(rkx), parameter :: bnmh = -axx*btg
  real(rkx), parameter :: bnmm = a1*(1.0_rkx-3.0_rkx*c1)
  real(rkx), parameter :: bshh = 9.0_rkx*a1*axx*axx*btg
  real(rkx), parameter :: bshm = 18.0_rkx*a1*a1*axx*c1
  real(rkx), parameter :: bsmh = -3.0_rkx*a1*axx * &
             (3.0_rkx*axx+3.0_rkx*b2*c1+12.0_rkx*a1*c1-b2)*btg
  real(rkx), parameter :: cesh = axx
  real(rkx), parameter :: cesm = a1*(1.0_rkx-3.0_rkx*c1)
  !
  ! Free term in the equilibrium equation for (L/Q)**2
  !
  real(rkx), parameter :: aeqh = 9.0_rkx*a1*axx*axx*b1*btg*btg +  &
              9.0_rkx*a1*axx*axx*(12.0_rkx*a1+3.0_rkx*b2)*btg*btg
  real(rkx), parameter :: aeqm = 3.0_rkx*a1*axx*b1*(3.0_rkx*axx + &
              3.0_rkx*b2*c1+18.0_rkx*a1*c1-b2)*btg+18.0_rkx*a1 *  &
              a1*axx*(b2-3.0_rkx*axx)*btg
  !
  ! Forbidden turbulence area
  !
  real(rkx), parameter :: requ = -aeqh/aeqm
  real(rkx), parameter :: epsgh = 1.e-9_rkx
  real(rkx), parameter :: epsgm = requ*epsgh
  !
  ! Near isotropy for shear turbulence, ww/q2 lower limit
  !
  real(rkx), parameter :: ubryl = (18.0_rkx*requ*a1*a1*axx*b2*c1*btg + &
               9.0_rkx*a1*axx*axx*b2*btg*btg)/(requ*adnm+adnh)
  real(rkx), parameter :: ubry = (d_one+epsrs)*ubryl
  real(rkx), parameter :: ubry3 = 3.0_rkx*ubry

  real(rkx), parameter :: aubh = 27.0_rkx*a1*axx*axx*b2*btg*btg-adnh*ubry3
  real(rkx), parameter :: aubm = 54.0_rkx*a1*a1*axx*b2*c1*btg-adnm*ubry3
  real(rkx), parameter :: bubh = (9.0_rkx*a1*axx + &
              3.0_rkx*axx*b2)*btg-bdnh*ubry3
  real(rkx), parameter :: bubm = 18.0_rkx*a1*a1*c1-bdnm*ubry3
  real(rkx), parameter :: cubr = d_one-ubry3
  real(rkx), parameter :: rcubr = d_one/cubr

  real(rkx), parameter :: countergrad = d_zero

  integer(ik4) :: nspec, ispb

  contains

  subroutine init_myjpbl
    implicit none
    if ( ipptls > 1 ) then
      nspec = 4
    else
      nspec = 3
    end if
    if ( ichem == 1 ) then
      ispb  = nspec
      nspec = nspec + ntr
    end if
  end subroutine init_myjpbl

  subroutine myjpbl(m2p,p2m)
    implicit none
    type(mod_2_pbl), intent(in) :: m2p
    type(pbl_2_mod), intent(inout) :: p2m
    !
    ! nspec is the number of mass species to be vertically mixed
    !
    integer(ik4) :: i, j, k, n, lmxl, nums
    real(rkx) :: akms_dens, dqdt, dtdif, dtdt, &
          dtturbl, qold, ratiomx,rdtturbl, thnew, thold, tx, exner, ct
    real(rkx), dimension(nspec) :: rkbot, sflux, cts, sz0
    real(rkx), dimension(kz) :: cwmk, pk, q2k, qk, thek , &
            tk, uk, vk, qcwk, qcik
    real(rkx), dimension(nspec,kz) :: species
    real(rkx), dimension(kzm1) :: akhk, akmk, el, gh, gm
    real(rkx), dimension(kzp1) :: zhk
    real(rkx), dimension(kz) :: rhok
    real(rkx), dimension(jci1:jci2,ici1:ici2,kz) :: ape, the, th, cwm
    real(rkx), dimension(jci1:jci2,ici1:ici2,kzm1) :: akh, akm
    real(rkx), dimension(jci1:jci2,ici1:ici2,kzp1) :: zint
    real(rkx) :: pfac
    ! 2-D fields from the surface layer call
    real(rkx), dimension(jci1:jci2,ici1:ici2) :: ustar

    dtturbl = dt
    rdtturbl = d_one/dtturbl
    dtdif = dtturbl
    ct = countergrad

    do concurrent ( j = jci1:jci2, i = ici1:ici2,  k = 1:kzm1 )
      akm(j,i,k) = d_zero
    end do

    do concurrent ( j = jci1:jci2, i = ici1:ici2,  k = 1:kzp1 )
      zint(j,i,k) = m2p%zq(j,i,k) + m2p%ht(j,i) * regrav
    end do

    do concurrent ( j = jci1:jci2, i = ici1:ici2,  k = 1:kz )
      exner = (m2p%patm(j,i,k)/p00)**rovcp
      ape(j,i,k) = d_one/exner
      tx = m2p%tatm(j,i,k)
      th(j,i,k) = tx * ape(j,i,k)
      cwm(j,i,k) = m2p%qxatm(j,i,k,iqc)
      if ( ipptls > 1 ) cwm(j,i,k) = cwm(j,i,k)+m2p%qxatm(j,i,k,iqi)
      the(j,i,k) = (cwm(j,i,k)*(-elocp/tx)+d_one)*th(j,i,k)
    end do

    !
    ! Janjic Eta surface layer (Janjic 1994; NCEP ON 437): ustar, the
    ! lower boundary condition of the TKE. The surface fluxes are the
    ! ones of the surface scheme (see main_integration).
    !
    call myj_sfclayer(m2p, ustar)

    setup_integration: &
    do i = ici1, ici2
      do j = jci1, jci2
        !
        ! Fill 1-d vertical arrays: myj scheme counts downward from
        ! the domain's top
        !
        do k = 1, kz
          zhk(k) = zint(j,i,k)
          tk(k) = m2p%tatm(j,i,k)
          thek(k) = the(j,i,k)
          ratiomx = m2p%qxatm(j,i,k,iqv)
          qk(k) = ratiomx/(d_one+ratiomx)
          cwmk(k) = cwm(j,i,k)
          pk(k) = m2p%patm(j,i,k)
          uk(k) = m2p%uxatm(j,i,k)
          vk(k) = m2p%vxatm(j,i,k)
          !
          ! TKE = 0.5*(q**2) ==> q**2 = 2.0*TKE
          !
          q2k(k) = d_two*m2p%tkests(j,i,k)
        end do
        zhk(kzp1) = zint(j,i,kzp1)
        !
        ! Find the mixing length
        !
        call mixlen(uk,vk,tk,thek,qk,cwmk,q2k,zhk,gm,gh,el, &
                    p2m%zpbl(j,i),p2m%kpbl(j,i),lmxl,ct)
        !
        ! Solve for the production/dissipation of the turbulent kinetic energy
        !
        call prodq2(dtsec,ustar(j,i),gm,gh,el,q2k)
        !
        ! Find the exchange coefficients in the free atmosphere
        !
        call difcof(gm,gh,el,q2k,zhk,akmk,akhk)
        !
        ! Counting downward from the top, the exchange coefficients akh
        ! are defined on the bottoms of the layers 1 to kzm1.
        !
        do k = 1, kzm1
          akh(j,i,k) = akhk(k)
          akm(j,i,k) = akmk(k)
        end do
        !
        ! Carry out the vertical diffusion of turbulent kinetic energy
        !
        call vdifq(dtsec,q2k,el,zhk)
        !
        ! Save the new TKE and mixing length.
        !
        do k = 1, kz
          q2k(k) = max(q2k(k),epsq2)
          p2m%tkepbl(j,i,k) = d_half*q2k(k)
        end do
      end do
    end do setup_integration

    main_integration: &
    do i = ici1, ici2
      do j = jci1, jci2
        !
        ! Fill 1-d vertical arrays: myj scheme counts downward
        ! from the domain's top
        !
        do k = 1, kz
          zhk(k) = zint(j,i,k)
          thek(k) = the(j,i,k)
          ratiomx = m2p%qxatm(j,i,k,iqv)
          qk(k) = ratiomx/(d_one+ratiomx)
          cwmk(k) = cwm(j,i,k)
          qcwk(k) = m2p%qxatm(j,i,k,iqc)
          species(1,k) = thek(k)
          species(2,k) = qk(k)
          species(3,k) = qcwk(k)
          if ( ipptls > 1 ) then
            qcik(k) = m2p%qxatm(j,i,k,iqi)
            species(4,k) = qcik(k)
          end if
          rhok(k) = m2p%patm(j,i,k)/(rgas*m2p%tatm(j,i,k)*  &
                    (d_one+ep1*qk(k)-cwmk(k)))
        end do
        zhk(kzp1) = zint(j,i,kzp1)
        if ( ichem == 1 ) then
          do n = 1, ntr
            do k = 1, kz
              species(ispb+n,k) = m2p%chib(j,i,k,n)
            end do
          end do
        end if
        !
        ! Counting downward from the top, the exchange coefficients akh
        ! are defined on the bottoms of the layers 1 to kz-1.
        ! These coefficients are also multiplied by the density at the
        ! bottom interface level.
        !
        do k = 1, kzm1
          akhk(k) = akh(j,i,k)*d_half*(rhok(k)+rhok(k+1))
        end do

        !
        ! Surface boundary conditions: the sensible heat and moisture
        ! fluxes of the surface scheme are applied as they are, so that
        ! energy and water are conserved. The tracers have no flux (the
        ! emissions are added by the chemistry) and, with ichdrdepo = 2,
        ! an implicit dry deposition sink rho*vd*chi.
        !
        do nums = 1, nspec
          sz0(nums) = d_zero
          rkbot(nums) = d_zero
          sflux(nums) = d_zero
        end do
        sflux(1) = m2p%hfx(j,i)*ape(j,i,kz)/cpd
        sflux(2) = m2p%qfx(j,i)
        if ( ichem == 1 .and. ichdrdepo == 2 ) then
          do n = 1, ntr
            rkbot(ispb+n) = rhok(kz)*m2p%drydepv(j,i,n)
          end do
        end if

        cts(1) = ct
        do nums = 2, nspec
          cts(nums) = d_zero
        end do
        !
        ! Carry out the vertical diffusion of temperature and water vapor
        !
        call vdifh(dtdif,p2m%kpbl(j,i),sz0,rkbot,sflux,cts, &
                   species,nspec,akhk,zhk,rhok)
        !
        ! Compute primary variable tendencies.
        ! pfac: sigma-pressure tendency factor for the hydrostatic core;
        !
        do k = 1, kz
          thek(k) = species(1,k)
          qk  (k) = species(2,k)
          qcwk(k) = species(3,k)
          cwmk(k) = qcwk(k)
          if ( ipptls > 1 ) then
            qcik(k) = species(4,k)
            cwmk(k) = cwmk(k) + qcik(k)
          end if
        end do

        if ( idynamic == 3 ) then
          pfac = d_one
        else
          pfac = m2p%psb(j,i)
        end if
        do k = 1, kz
          thold = th(j,i,k)
          thnew = thek(k) + cwmk(k)*elocp*ape(j,i,k)
          dtdt  = (thnew-thold)*rdtturbl
          qold  = m2p%qxatm(j,i,k,iqv)/(d_one+m2p%qxatm(j,i,k,iqv))
          dqdt  = (qk(k)-qold)*rdtturbl
          exner = (m2p%patm(j,i,k)/p00)**rovcp
          p2m%tten(j,i,k)     = p2m%tten(j,i,k) + dtdt * exner * pfac
          p2m%qxten(j,i,k,iqv) = p2m%qxten(j,i,k,iqv) + &
                     (dqdt/(d_one-qk(k))**2) * pfac
          p2m%qxten(j,i,k,iqc) = p2m%qxten(j,i,k,iqc) + &
                     (qcwk(k)-m2p%qxatm(j,i,k,iqc))*rdtturbl*pfac
          if ( ipptls > 1 ) then
            p2m%qxten(j,i,k,iqi) = p2m%qxten(j,i,k,iqi) + &
                     (qcik(k)-m2p%qxatm(j,i,k,iqi))*rdtturbl*pfac
          end if
        end do
        if ( ichem == 1 ) then
          do n = 1, ntr
            do k = 1, kz
              p2m%chiten(j,i,k,n) = p2m%chiten(j,i,k,n) + &
                   (species(ispb+n,k)-m2p%chib(j,i,k,n))*rdtturbl*pfac
            end do
          end do
          if ( ichdrdepo == 2 ) then
            ! Deposited amount, as accumulated by the Holtslag scheme
            do n = 1, ntr
              p2m%remdrd(j,i,n) = p2m%remdrd(j,i,n) + &
                  m2p%chib(j,i,kz,n)*m2p%drydepv(j,i,n)*dtturbl * &
                  d_half*m2p%rhox2d(j,i)*egrav / &
                  (m2p%patmf(j,i,kzp1)-m2p%patmf(j,i,kz))*pfac
            end do
          end if
        end if
        !
        ! Density-weight the momentum exchange coefficients and diffuse
        ! the velocity components.
        !
        do k = 1, kzm1
          akmk(k) = akm(j,i,k) * (rhok(k)+rhok(k+1)) * d_half
        end do

        ! The surface stress of the surface scheme, uvdrag*|U|, applied
        ! implicitly with no slip at the surface
        akms_dens = m2p%uvdrag(j,i)     ! [kg m^-2 s^-1]
        do k = 1, kz
          uk(k) = m2p%uxatm(j,i,k)
          vk(k) = m2p%vxatm(j,i,k)
        end do
        !
        ! Carry out the vertical diffusion of velocity components
        !
        call vdifv(dtdif,akms_dens,uk,vk,akmk,zhk,rhok)
        !
        ! Compute momentum tendencies
        !
        do k = 1, kz
          p2m%uten(j,i,k) = p2m%uten(j,i,k) + &
                             (uk(k)-m2p%uxatm(j,i,k))*rdtturbl
          p2m%vten(j,i,k) = p2m%vten(j,i,k) + &
                             (vk(k)-m2p%vxatm(j,i,k))*rdtturbl
        end do
      end do
    end do main_integration

  end subroutine myjpbl
  !
  ! Level 2.5 Mixing Length
  !
  subroutine mixlen(u,v,t,the,q,cwm,q2,z,gm,gh,el,pblh,lpbl,lmxl,ct)
    implicit none
    integer(ik4), intent(out) :: lmxl, lpbl
    real(rkx), dimension(kz), intent(in) :: cwm, q, q2, t, the, u, v
    real(rkx), dimension(kzp1), intent(in) :: z
    real(rkx), intent(out) :: pblh
    real(rkx), dimension(kzm1), intent(out) :: el, gh, gm
    real(rkx), intent(inout) :: ct
    integer(ik4) :: k, lpblm
    real(rkx) :: a, aden, b, bden, aubr, bubr, blmx, el0, &
      eloq2, ghl, gml, qol2st, qol2un, qdzl, rdz, sq,    &
      srel, szq, tem, thm, vkrmz
    real(rkx), dimension(kz) :: q1
    real(rkx), dimension(kzm1) :: dth, elm, rel
    logical :: lfound
    !
    ! Find the height of the pbl
    !
    !
    ! Scan upward (k decreasing = moving away from surface) to find the
    ! top of the turbulent layer.  Stop when q2 drops to the minimum
    ! value OR when the interface height AGL exceeds 6000 m.
    !
    lpbl = kz
    lfound = .false.
    do k = kzm1, 1, -1
      if ( q2(k) <= epsq2*fh .or. (z(k+1)-z(kzp1)) > 6000.0_rkx ) then
        lpbl = k
        lfound = .true.
        pblh = z(lpbl+1) - z(kzp1)
        exit
      end if
    end do
    if ( .not. lfound ) then
      lpbl = 1
      pblh = z(lpbl+1) - z(kzp1)
    end if
    !
    ! The height of the pbl
    !
    do k = 1, kz
      q1(k) = d_zero
    end do
    do k = 1, kzm1
      dth(k) = the(k)-the(k+1)
    enddo
    do k = kz-2, 1, -1
      if ( dth(k) > d_zero .and. dth(k+1) <= d_zero ) then
        dth(k) = dth(k)+ct
        exit
      end if
    end do
    ct = d_zero
    do k = 1, kzm1
      rdz = d_two/(z(k)-z(k+2))
      gml = ((u(k)-u(k+1))**2+(v(k)-v(k+1))**2)*rdz*rdz
      gm(k) = max(gml,epsgm)
      tem = (t(k)+t(k+1))*d_half
      thm = (the(k)+the(k+1))*d_half
      a = thm*ep1
      b = (elocp/tem-d_one-ep1)*thm
      ghl = (dth(k)*((q(k)+q(k+1)+cwm(k)+cwm(k+1))*(d_half*ep1)+d_one) + &
            (q(k)-q(k+1)+cwm(k)-cwm(k+1))*a + (cwm(k)-cwm(k+1))*b)*rdz
      if ( abs(ghl) <= epsgh ) ghl = epsgh
      gh(k) = ghl
    end do
    !
    ! Find maximum mixing lengths and the level of the pbl top
    !
    lmxl = kz
    do k = 1, kzm1
      gml = gm(k)
      ghl = gh(k)
      if ( ghl >= epsgh ) then
        if ( gml/ghl <= requ ) then
          elm(k) = epsl
          lmxl   = k
        else
          aubr   = (aubm*gml+aubh*ghl)*ghl
          bubr   = bubm*gml+bubh*ghl
          qol2st = (-d_half*bubr+sqrt(bubr*bubr*0.25_rkx-aubr*cubr))*rcubr
          eloq2  = d_one/qol2st
          elm(k) = max(sqrt(eloq2*q2(k)),epsl)
        endif
      else
        aden   = (adnm*gml+adnh*ghl)*ghl
        bden   = bdnm*gml+bdnh*ghl
        qol2un = -d_half*bden+sqrt(bden*bden*0.25_rkx-aden)
        eloq2  = d_one/(qol2un+epsru)
        elm(k) = max(sqrt(eloq2*q2(k)),epsl)
      end if
    end do
    if ( elm(kzm1) == epsl ) lmxl = kz
    !
    ! The height of the mixed layer
    !
    blmx  = z(lmxl)-z(kzp1)
    do k = lpbl, kz
      q1(k) = sqrt(q2(k))
    enddo
    szq = d_zero
    sq  = d_zero
    do k = 1, kzm1
      qdzl = (q1(k)+q1(k+1))*(z(k+1)-z(k+2))
      szq  = (z(k+1)+z(k+2)-z(kzp1)-z(kzp1))*qdzl+szq
      sq   = qdzl+sq
    end do
    !
    ! Computation of asymptotic l in blackadar formula.
    ! el0 is bounded above by the mixed-layer depth (blmx) to prevent
    ! unrealistically large mixing lengths in deep convective situations.
    !
    el0 = min(alph*szq*d_half/sq, el0max)
    if ( lmxl < kz ) el0 = min(el0, blmx) ! cap at mixed-layer depth
    el0 = max(el0, el0min)
    !
    ! Above the pbl top
    !
    lpblm = max(lpbl-1,1)
    do k = 1, lpblm
      el(k)  = min((z(k)-z(k+2))*elfc,elm(k))
      rel(k) = el(k)/elm(k)
    end do
    !
    ! Inside the pbl
    !
    if ( lpbl < kz ) then
      do k = lpbl, kzm1
        vkrmz = (z(k+1)-z(kzp1))*vonkar
        el(k) = min(vkrmz/(vkrmz/el0+d_one),elm(k))
        rel(k) = el(k)/elm(k)
      end do
    end if
    do k = lpbl+1, kz-2
      srel  = min(((rel(k-1)+rel(k+1))*d_half+rel(k))*d_half,rel(k))
      el(k) = max(srel*elm(k),epsl)
    end do
  end subroutine mixlen
  !
  ! Level 2.5 Q2 Production/Dissipation
  !
  subroutine prodq2(dtturbl,ustar,gm,gh,el,q2)
    implicit none
    real(rkx), intent(in) :: dtturbl, ustar
    real(rkx), dimension(kzm1), intent(in) :: gh, gm
    real(rkx), dimension(kzm1), intent(inout) :: el
    real(rkx), dimension(kz), intent(inout) :: q2
    integer(ik4) :: k
    real(rkx) :: aden, aequ, anum, arhs, bden, bequ, bnum, brhs, &
      cden, crhs, dloq1, eloq11, eloq12, eloq13, eloq21, eloq22, &
      eloq31, eloq32, eloq41, eloq42, eloq51, eloq52, eloqn,      &
      eqol2, ghl, gml, rden1, rden2, rhs2, rhsp1, rhsp2, rhst2

    main_integration: &
    do k = 1, kzm1
      gml = gm(k)
      ghl = gh(k)
      !
      ! Coefficients of the equilibrium equation
      !
      aequ = (aeqm*gml+aeqh*ghl)*ghl
      bequ = beqm*gml+beqh*ghl
      !
      ! Equilibrium solution for L/Q
      !
      eqol2 = -d_half*bequ+sqrt(bequ*bequ*0.25_rkx-aequ)
      !
      ! Is there production/dissipation ?
      !
      if ( (gml+ghl*ghl <= epstrb) .or. &
           (ghl >= epsgh .and. gml/ghl <= requ) .or. &
           (eqol2 <= eps2) ) then
        !
        ! No Turbulence
        !
        q2(k) = epsq2
        el(k) = epsl
      else
        !
        ! Turbulence
        !
        ! Coefficients of the terms in the numerator
        !
        anum = (anmm*gml+anmh*ghl)*ghl
        bnum = bnmm*gml+bnmh*ghl
        !
        ! Coefficients of the terms in the denominator
        !
        aden = (adnm*gml+adnh*ghl)*ghl
        bden = bdnm*gml+bdnh*ghl
        cden = d_one
        !
        ! Coefficients of the numerator of the linearized eq.
        !
        arhs = -(anum*bden-bnum*aden)*d_two
        brhs = - anum*d_four
        crhs = - bnum*d_two
        !
        ! Initial value of L/Q
        !
        dloq1 = el(k)/sqrt(q2(k))
        !
        ! First iteration for L/Q, RHS=0
        !
        eloq21 = d_one/eqol2
        eloq11 = sqrt(eloq21)
        eloq31 = eloq21*eloq11
        eloq41 = eloq21*eloq21
        eloq51 = eloq21*eloq31
        !
        ! 1./Denominator
        !
        rden1 = d_one/(aden*eloq41+bden*eloq21+cden)
        !
        ! d(RHS)/d(L/Q)
        !
        rhsp1 = (arhs*eloq51+brhs*eloq31+crhs*eloq11)*rden1*rden1
        !
        ! First-guess solution
        !
        eloq12 = eloq11+(dloq1-eloq11)*exp(rhsp1*dtturbl)
        eloq12 = max(eloq12,eps1)
        !
        ! Second iteration for L/Q
        !
        eloq22 = eloq12*eloq12
        eloq32 = eloq22*eloq12
        eloq42 = eloq22*eloq22
        eloq52 = eloq22*eloq32
        !
        ! 1./Denominator
        !
        rden2 = d_one/(aden*eloq42+bden*eloq22+cden)
        rhs2  = -(anum*eloq42+bnum*eloq22)*rden2+rb1
        rhsp2 = (arhs*eloq52+brhs*eloq32+crhs*eloq12)*rden2*rden2
        rhst2 = rhs2/rhsp2
        !
        ! Corrected solution
        !
        eloq13 = eloq12-rhst2+(rhst2+dloq1-eloq12)*exp(rhsp2*dtturbl)
        eloq13 = max(eloq13,eps1)
        !
        ! Two iterations is enough in most cases ...
        !
        eloqn = eloq13
        if ( eloqn > eps1 ) then
          q2(k) = el(k)*el(k)/(eloqn*eloqn)
          q2(k) = max(q2(k),epsq2)
          if ( q2(k) == epsq2 ) then
            el(k) = epsl
          end if
        else
          q2(k) = epsq2
          el(k) = epsl
        end if
        !
        ! End of Turbulent Branch
        !
      end if
      ! End of Production/dissipation Loop
    end do main_integration
    !
    ! Lower boundary condition for Q2
    !
    q2(kz) = max(b1**(d_two/d_three)*ustar*ustar,epsq2)
  end subroutine prodq2
  !
  ! Level 2.5 Diffusion Coefficients
  !
  subroutine difcof(gm,gh,el,q2,z,akm,akh)
    implicit none
    real(rkx), dimension(kz), intent(in) :: q2
    real(rkx), dimension(kzm1), intent(in) :: el, gh, gm
    real(rkx), dimension(kzp1), intent(in) :: z
    real(rkx), dimension(kzm1), intent(out) :: akh, akm
    integer(ik4) :: k
    real(rkx) :: aden, bden, besh, besm, cden, &
      ell, eloq2, eloq4, elqdz, esh, esm, ghl, gml, q1l, &
      rden, rdz

    do k = 1, kzm1
      ell   = el(k)
      eloq2 = ell*ell/q2(k)
      eloq4 = eloq2*eloq2
      gml   = gm(k)
      ghl   = gh(k)
      !
      ! Coefficients of the terms in the denominator
      !
      aden = (adnm*gml+adnh*ghl)*ghl
      bden = bdnm*gml+bdnh*ghl
      cden = d_one
      !
      ! Coefficients for the sm determinant
      !
      besm = bsmh*ghl
      !
      ! Coefficients for the sh determinant
      !
      besh = bshm*gml+bshh*ghl
      !
      !  1./Denominator
      !
      rden = d_one/(aden*eloq4+bden*eloq2+cden)
      !
      ! SM and SH
      !
      esm = (besm*eloq2+cesm)*rden
      esh = (besh*eloq2+cesh)*rden
      !
      ! Diffusion Coefficients
      !
      rdz = d_two/(z(k)-z(k+2))
      q1l = sqrt(q2(k))
      elqdz = ell*q1l*rdz
      akm(k) = elqdz*esm
      akh(k) = elqdz*esh
    end do
  end subroutine difcof
  !
  ! Vertical turbulent diffusion of q2 (tke)
  !
  subroutine vdifq(dtdif,q2,el,z)
    implicit none
    real(rkx) ,intent(in) :: dtdif
    real(rkx) ,dimension(kzm1),intent(in) :: el
    real(rkx) ,dimension(kzp1),intent(in) :: z

    real(rkx) ,dimension(kz),intent(inout) :: q2

    integer(ik4) :: k
    real(rkx) :: akqs, cf, dtozs, esqhf
    real(rkx), dimension(kzm2) :: akq, cm, cr, dtoz, rsq2

    esqhf = d_half*esq

    do k = 1, kzm2
      dtoz(k) = (dtdif+dtdif)/(z(k)-z(k+2))
      akq(k)  = sqrt((q2(k)+q2(k+1))*d_half)*(el(k)+el(k+1))* &
                esqhf/(z(k+1)-z(k+2))
      cr(k)   = -dtoz(k)*akq(k)
    end do

    cm(1)   = dtoz(1)*akq(1)+d_one
    rsq2(1) = q2(1)

    do k = 2, kz-2
      cf      = -dtoz(k)*akq(k-1)/cm(k-1)
      cm(k)   = -cr(k-1)*cf+(akq(k-1)+akq(k))*dtoz(k)+d_one
      rsq2(k) = -rsq2(k-1)*cf+q2(k)
    end do

    dtozs = (dtdif+dtdif)/(z(kzm1)-z(kzp1))
    akqs = sqrt((q2(kzm1)+q2(kz))*d_half)*(el(kzm1)+elz0) * &
           esqhf/(z(kz)-z(kzp1))
    cf = -dtozs*akq(kz-2)/cm(kz-2)

    q2(kzm1) = (dtozs*akqs*q2(kz)-rsq2(kz-2)*cf+q2(kzm1)) / &
                ((akq(kz-2)+akqs)*dtozs-cr(kz-2)*cf+d_one)

    do k = kz-2, 1, -1
      q2(k) = (-cr(k)*q2(k+1)+rsq2(k))/cm(k)
    end do
  end subroutine vdifq
  !
  ! Vertical diffusion of mass variables
  !
  ! At the surface, each species has the flux sflux plus the implicit
  ! exchange rkbot*(sz0 - species(kz)).
  !
  subroutine vdifh(dtdif,lpbl,sz0,rkbot,sflux,cts, &
                   species,ns,rkh,zhk,rho)
    implicit none
    integer(ik4), intent(in) :: lpbl, ns
    real(rkx), intent(in) :: dtdif
    real(rkx), dimension(nspec), intent(in) :: rkbot, sflux, cts, sz0
    real(rkx), dimension(kzm1), intent(in) :: rkh
    real(rkx), dimension(kz), intent(in) :: rho
    real(rkx), dimension(kzp1), intent(in) :: zhk
    real(rkx), dimension(nspec,kz), intent(inout) :: species

    integer(ik4) :: k, m
    real(rkx) :: cf, cmb, cmsb, dtozl, dtozs, rcml, rhok ,&
      rkhh, rkhz, rkss, rssb
    real(rkx), dimension(kzm1) :: cm, cr, dtoz
    real(rkx), dimension(ns,kzm1) :: rkct, rss

    do k = 1, kzm1
      dtoz(k) = dtdif/(zhk(k)-zhk(k+1))
      cr(k) = -dtoz(k)*rkh(k)
      if ( k < lpbl ) then
        do m = 1, ns
          rkct(m,k) = d_zero
        end do
      else
        rkhz = rkh(k)*(zhk(k)-zhk(k+2))
        do m = 1, ns
          rkct(m,k) = d_half*rkhz*cts(m)
        end do
      end if
    end do
    ! Top level
    rhok  = rho(1)
    cm(1) = dtoz(1)*rkh(1)+rhok
    do m = 1, ns
      rss(m,1) = -rkct(m,1)*dtoz(1)+species(m,1)*rhok
    end do
    ! Intermediate levels
    do k = 2, kzm1
      dtozl = dtoz(k)
      cf    = -dtozl*rkh(k-1)/cm(k-1)
      rhok  = rho(k)
      cm(k) = -cr(k-1)*cf+(rkh(k-1)+rkh(k))*dtozl+rhok
      do m = 1, ns
        rss(m,k) = -rss(m,k-1)*cf + &
                   (rkct(m,k-1)-rkct(m,k))*dtozl+species(m,k)*rhok
      end do
    end do
    ! Bottom level
    dtozs = dtdif/(zhk(kz)-zhk(kzp1))
    rkhh  = rkh(kzm1)
    cf    = -dtozs*rkhh/cm(kzm1)
    cmb   = cr(kzm1)*cf
    rhok  = rho(kz)
    do m = 1, ns
      rkss = rkbot(m)
      cmsb = -cmb+(rkhh+rkss)*dtozs+rhok
      rssb = -rss(m,kzm1)*cf+rkct(m,kzm1)*dtozs+species(m,kz)*rhok
      species(m,kz) = (dtozs*(rkss*sz0(m)+sflux(m))+rssb)/cmsb
    end do
    ! Backsubstitution
    do k = kzm1, 1, -1
      rcml = d_one/cm(k)
      do m = 1, ns
        species(m,k) = (-cr(k)*species(m,k+1)+rss(m,k))*rcml
      end do
    end do
  end subroutine vdifh
  !
  ! Vertical diffusion of velocity components
  !
  ! At the surface the stress is rkms*u, with no slip.
  !
  subroutine vdifv(dtdif,rkms,u,v,rkm,z,rho)
    implicit none
    real(rkx), intent(in) :: rkms, dtdif
    real(rkx), dimension(kzm1), intent(in) :: rkm
    real(rkx), dimension(kz), intent(in) :: rho
    real(rkx), dimension(kzp1), intent(in) :: z
    real(rkx), dimension(kz), intent(inout) :: u, v
    integer(ik4) :: k
    real(rkx) :: cf, dtozl, dtozs, rcml, rcmvb, rhok, rkmh
    real(rkx), dimension(kzm1) :: cm, cr, dtoz, rsu, rsv

    do k = 1, kzm1
      dtoz(k) = dtdif/(z(k)-z(k+1))
      cr(k)   = -dtoz(k)*rkm(k)
    end do

    rhok   = rho(1)
    cm(1)  = dtoz(1)*rkm(1)+rhok
    rsu(1) = u(1)*rhok
    rsv(1) = v(1)*rhok

    do k = 2, kzm1
      dtozl  = dtoz(k)
      cf     = -dtozl*rkm(k-1)/cm(k-1)
      rhok   = rho(k)
      cm(k)  = -cr(k-1)*cf+(rkm(k-1)+rkm(k))*dtozl+rhok
      rsu(k) = -rsu(k-1)*cf+u(k)*rhok
      rsv(k) = -rsv(k-1)*cf+v(k)*rhok
    end do

    dtozs  = dtdif/(z(kz)-z(kzp1))
    rkmh   = rkm(kzm1)
    cf     = -dtozs*rkmh/cm(kzm1)
    rhok   = rho(kz)
    rcmvb  = d_one/((rkmh+rkms)*dtozs-cr(kzm1)*cf+rhok)
    u(kz) = (-rsu(kzm1)*cf+u(kz)*rhok)*rcmvb
    v(kz) = (-rsv(kzm1)*cf+v(kz)*rhok)*rcmvb

    do k = kzm1, 1, -1
      rcml = d_one/cm(k)
      u(k) = (-cr(k)*u(k+1)+rsu(k))*rcml
      v(k) = (-cr(k)*v(k+1)+rsv(k))*rcml
    end do
  end subroutine vdifv

end module mod_pbl_myj

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

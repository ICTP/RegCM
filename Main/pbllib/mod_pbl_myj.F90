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

module mod_pbl_myj
  !
  !
  ! REFERENCES:  Janjic (2002), NCEP Office Note 437
  !              Mellor and Yamada (1982), Rev. Geophys. Space Phys.
  !
  ! ABSTRACT:
  !     MYJ updates the turbulent kinetic energy with the production/
  !     dissipation term and the vertical diffusion term
  !     (using an implicit formulation) from Mellor-Yamada
  !     level 2.5 as extended by Janjic.  Exchange coefficients for
  !     the surface and for all layer interfaces are computed from
  !     Monin-Obukhov theory.
  !     The turbulent vertical exchange is then executed.
  !
  use mod_intkinds
  use mod_realkinds
  use mod_constants
  use mod_dynparam
  use mod_runparams

  implicit none

  private

  public :: init_myjpbl , myjpbl

  ! iteration count for mixing length computation
  integer(ik4) , parameter :: itrmx = 5

  real(rkx) , parameter :: rlivwv = wlhs/wlhv
  real(rkx) , parameter :: elocp = 2.72e6_rkx/cpd
  real(rkx) , parameter :: epsq2 = 0.2_rkx
  real(rkx) , parameter :: eps1 = 1.e-12_rkx
  real(rkx) , parameter :: eps2 = 0.0_rkx
  real(rkx) , parameter :: epsl = 0.32_rkx
  real(rkx) , parameter :: epsru = 1.e-7_rkx
  real(rkx) , parameter :: epsrs = 1.e-7_rkx
  real(rkx) , parameter :: epstrb = 1.e-24_rkx
  real(rkx) , parameter :: epsa = 1.e-8_rkx
  real(rkx) , parameter :: epsit = 1.e-4_rkx
  real(rkx) , parameter :: epsu2 = 1.e-4_rkx
  real(rkx) , parameter :: epsust = 0.07_rkx
  real(rkx) , parameter :: fh = 1.01_rkx
  real(rkx) , parameter :: alph = 0.30_rkx
  real(rkx) , parameter :: beta = d_one/273.0_rkx
  real(rkx) , parameter :: el0max = 1000.0_rkx
  real(rkx) , parameter :: el0min = 1.0_rkx
  real(rkx) , parameter :: elfc = 0.23_rkx*0.5_rkx
  real(rkx) , parameter :: gam1 = 0.2222222222222222222_rkx
  real(rkx) , parameter :: prt = 1.0_rkx
  real(rkx) , parameter :: a1 = 0.659888514560862645_rkx
  real(rkx) , parameter :: a2 = 17.2693882_rkx
  real(rkx) , parameter :: a3 = 273.16_rkx
  real(rkx) , parameter :: a4 = 35.86_rkx
  real(rkx) , parameter :: a2x = 0.6574209922667784586_rkx
  real(rkx) , parameter :: b1 = 11.87799326209552761_rkx
  real(rkx) , parameter :: b2 = 7.226971804046074028_rkx
  real(rkx) , parameter :: c1 = 0.000830955950095854396_rkx
  real(rkx) , parameter :: a2s = 17.2693882_rkx
  real(rkx) , parameter :: a3s = 273.16_rkx
  real(rkx) , parameter :: a4s = 35.86_rkx
  real(rkx) , parameter :: elz0 = 0.0_rkx
  real(rkx) , parameter :: esq = 5.0_rkx
  real(rkx) , parameter :: excm = 0.001_rkx
  real(rkx) , parameter :: fhneu = 0.8_rkx
  real(rkx) , parameter :: glkbr = 10.0_rkx
  real(rkx) , parameter :: glkbs = 30.0_rkx
  real(rkx) , parameter :: qvisc = 2.1e-5_rkx
  real(rkx) , parameter :: rfc = 0.191_rkx
  real(rkx) , parameter :: ric = 0.505_rkx
  real(rkx) , parameter :: small = 0.35_rkx
  real(rkx) , parameter :: sqpr = 0.84_rkx
  real(rkx) , parameter :: sqsc = 0.84_rkx
  real(rkx) , parameter :: sqvisc = 258.2_rkx
  real(rkx) , parameter :: tvisc = 2.1e-5_rkx
  real(rkx) , parameter :: ustc = 0.7_rkx
  real(rkx) , parameter :: ustr = 0.225_rkx
  real(rkx) , parameter :: visc = 1.5e-5_rkx
  real(rkx) , parameter :: wold = 0.15_rkx
  real(rkx) , parameter :: wwst = 1.2_rkx
  real(rkx) , parameter :: ztmax = 1.0_rkx
  real(rkx) , parameter :: ztfc = 1.0_rkx
  real(rkx) , parameter :: ztmin = -5.0_rkx
  real(rkx) , parameter :: seafc = 0.98_rkx
  real(rkx) , parameter :: pq0sea = pq0*seafc
  real(rkx) , parameter :: btg = beta*egrav
  real(rkx) , parameter :: cziv = small*glkbs
  real(rkx) , parameter :: esqhf = 0.5_rkx*5.0_rkx
  real(rkx) , parameter :: grrs = glkbr/glkbs
  real(rkx) , parameter :: rb1 = d_one/b1
  real(rkx) , parameter :: rtvisc = d_one/tvisc
  real(rkx) , parameter :: rvisc = d_one/visc
  real(rkx) , parameter :: zqrzt = sqsc/sqpr
  real(rkx) , parameter :: adnh = 9.0_rkx*a1*a2x*a2x * &
               (12.0_rkx*a1+3.0_rkx*b2)*btg*btg
  real(rkx) , parameter :: adnm = 18.0_rkx*a1*a1*a2x*(b2-3.0_rkx*a2x)*btg
  real(rkx) , parameter :: anmh = -9.0_rkx*a1*a2x*a2x*btg*btg
  real(rkx) , parameter :: anmm = -3.0_rkx*a1*a2x * &
               (3.0_rkx*a2x+3.0_rkx*b2*c1+18.0_rkx*a1*c1-b2)*btg
  real(rkx) , parameter :: bdnh = 3.0_rkx*a2x*(7.0_rkx*a1+b2)*btg
  real(rkx) , parameter :: bdnm = 6.0_rkx*a1*a1
  real(rkx) , parameter :: beqh = a2x*b1*btg+3.0_rkx*a2x*(7.0_rkx*a1+b2)*btg
  real(rkx) , parameter :: beqm = -a1*b1*(1.0_rkx-3.0_rkx*c1)+6.0_rkx*a1*a1
  real(rkx) , parameter :: bnmh = -a2x*btg
  real(rkx) , parameter :: bnmm = a1*(1.0_rkx-3.0_rkx*c1)
  real(rkx) , parameter :: bshh = 9.0_rkx*a1*a2x*a2x*btg
  real(rkx) , parameter :: bshm = 18.0_rkx*a1*a1*a2x*c1
  real(rkx) , parameter :: bsmh = -3.0_rkx*a1*a2x * &
               (3.0_rkx*a2x+3.0_rkx*b2*c1+12.0_rkx*a1*c1-b2)*btg
  real(rkx) , parameter :: cesh = a2x
  real(rkx) , parameter :: cesm = a1*(1.0_rkx-3.0_rkx*c1)
  real(rkx) , parameter :: cnv = ep1*egrav/btg
  real(rkx) , parameter :: elfcs = vonkar*btg
  real(rkx) , parameter :: fzq1 = rtvisc*qvisc*zqrzt
  real(rkx) , parameter :: fzq2 = rtvisc*qvisc*zqrzt
  real(rkx) , parameter :: fzt1 = rvisc*tvisc*sqpr
  real(rkx) , parameter :: fzt2 = cziv*grrs*tvisc*sqpr
  real(rkx) , parameter :: fzu1 = cziv*visc
  real(rkx) , parameter :: rfac = ric/(fhneu*rfc*rfc)
  real(rkx) , parameter :: rqvisc = d_one/qvisc
  real(rkx) , parameter :: rric = d_one/ric
  real(rkx) , parameter :: ustfc = 0.018_rkx*regrav
  real(rkx) , parameter :: wnew = d_one-wold
  real(rkx) , parameter :: wwst2 = wwst*wwst
  !
  ! Free term in the equilibrium equation for (L/Q)**2
  !
  real(rkx) , parameter :: aeqh = 9.0_rkx*a1*a2x*a2x*b1*btg*btg +  &
               9.0_rkx*a1*a2x*a2x*(12.0_rkx*a1+3.0_rkx*b2)*btg*btg
  real(rkx) , parameter :: aeqm = 3.0_rkx*a1*a2x*b1*(3.0_rkx*a2x + &
               3.0_rkx*b2*c1+18.0_rkx*a1*c1-b2)*btg+18.0_rkx*a1 *  &
               a1*a2x*(b2-3.0_rkx*a2x)*btg
  !
  ! Forbidden turbulence area
  !
  real(rkx) , parameter :: requ = -aeqh/aeqm
  real(rkx) , parameter :: epsgh = 1.e-9_rkx
  real(rkx) , parameter :: epsgm = requ*epsgh
  !
  ! Near isotropy for shear turbulence, ww/q2 lower limit
  !
  real(rkx) , parameter :: ubryl = (18.0_rkx*requ*a1*a1*a2x*b2*c1*btg + &
                9.0_rkx*a1*a2x*a2x*b2*btg*btg)/(requ*adnm+adnh)
  real(rkx) , parameter :: ubry = (d_one+epsrs)*ubryl
  real(rkx) , parameter :: ubry3 = 3.0_rkx*ubry

  real(rkx) , parameter :: aubh = 27.0_rkx*a1*a2x*a2x*b2*btg*btg-adnh*ubry3
  real(rkx) , parameter :: aubm = 54.0_rkx*a1*a1*a2x*b2*c1*btg-adnm*ubry3
  real(rkx) , parameter :: bubh = (9.0_rkx*a1*a2x + &
               3.0_rkx*a2x*b2)*btg-bdnh*ubry3
  real(rkx) , parameter :: bubm = 18.0_rkx*a1*a1*c1-bdnm*ubry3
  real(rkx) , parameter :: cubr = d_one-ubry3
  real(rkx) , parameter :: rcubr = d_one/cubr

  integer(ik4) :: nspec

  contains

  subroutine init_myjpbl
    implicit none
    nspec = 1 + nqx + ntr
  end subroutine init_myjpbl

  subroutine myjpbl(dt,stepbl,ht,dz,                               &
                    pmid,pint,th,t,exner,qv,qcw,qci,qcs,qcr,qcg,   &
                    u,v,rho,tsk,qsfc,chklowq,thz0,qz0,uz0,vz0,     &
                    lowlyr,xland,sice,snow,                        &
                    tke_myj,exch_h,ustar,znt,el_myj,pblh,kpbl,ct,  &
                    akhs,akms,elflx,mixht,                         &
                    rublten,rvblten,rthblten,rqvblten,rqcblten,    &
                    rqiblten,rqsblten,rqrblten,rqgblten)
    implicit none

    real(rkx) , intent(in) :: dt
    integer(ik4) , intent(in) :: stepbl
    integer(ik4) , dimension(jci1:jci2,ici1:ici2) , intent(in) :: lowlyr
    integer(ik4) , dimension(jci1:jci2,ici1:ici2) , intent(out) :: kpbl
    real(rkx) , dimension(jci1:jci2,ici1:ici2) , intent(in) :: ht
    real(rkx) , dimension(jci1:jci2,ici1:ici2) , intent(in) :: sice
    real(rkx) , dimension(jci1:jci2,ici1:ici2) , intent(in) :: snow
    real(rkx) , dimension(jci1:jci2,ici1:ici2) , intent(in) :: tsk
    real(rkx) , dimension(jci1:jci2,ici1:ici2) , intent(in) :: xland
    real(rkx) , dimension(jci1:jci2,ici1:ici2,kz),intent(in) :: dz,exner   &
                                                         ,pmid,pint  &
                                                         ,rho        &
                                                         ,t,th,u,v
    real(rkx) , intent(in) , dimension(jci1:jci2,ici1:ici2,kz)::    &
                                              qv,qcw,qci,qcs,qcr,qcg
    real(rkx) , dimension(jci1:jci2,ici1:ici2),intent(out) :: pblh,mixht
    real(rkx) , dimension(jci1:jci2,ici1:ici2),intent(inout) :: akhs,akms
    real(rkx) , intent(out), dimension(jci1:jci2,ici1:ici2,kz) ::         &
                                             el_myj,rthblten         &
                                            ,rublten,rvblten
    real(rkx) , intent(out),dimension(jci1:jci2,ici1:ici2,kz)::   &
                                                 rqcblten,rqvblten   &
                                                ,rqiblten,rqsblten   &
                                                ,rqrblten,rqgblten
    real(rkx) , dimension(jci1:jci2,ici1:ici2),intent(inout) :: ct,qsfc,qz0     &
                                                    ,thz0,ustar      &
                                                    ,uz0,vz0,znt
    real(rkx) , dimension(jci1:jci2,ici1:ici2,kz)                          &
        ,intent(inout) ::                    exch_h,tke_myj
    real(rkx) , dimension(jci1:jci2,ici1:ici2),intent(in) :: chklowq,elflx
    !
    ! nspec is the number of mass species to be vertically mixed
    ! It is assumed that cloud water is always being passed into the scheme.
    !
    integer(ik4) :: i,j,k,kflip,llow,lmh,lmxl, nums,nsp
    integer(ik4) , dimension(jci1:jci2,ici1:ici2) :: lpbl
    real(rkx) :: akhs_dens,akms_dens,apex,deltaz,dqdt,dtdif,dtdt     &
           ,dtturbl,dudt,dvdt,exnsfc,psfc,ptop,qfc1,qlow,qold   &
           ,ratiomx,rdtturbl,rg,rwmsk,seamask,thnew,thold,tx    &
           ,ulow,vlow,wmsk
    real(rkx) , dimension(nspec) :: clow , cts , sz0
    real(rkx) , dimension(kz) :: cwmk,pk,q2k,qk,thek,tk,uk,vk     &
                              ,qcwk,qcik,qcsk,qcrk,qcgk
    real(rkx) , dimension(nspec,kz) :: species
    real(rkx) , dimension(kzm1) :: akhk , akmk , el , gh , gm
    real(rkx) , dimension(kzp1) :: zhk
    real(rkx) , dimension(jci1:jci2,ici1:ici2) :: thsk
    real(rkx) , dimension(kz) :: rhok
    real(rkx) , dimension(jci1:jci2,ici1:ici2,kz) :: ape,the, cwm
    real(rkx) , dimension(jci1:jci2,ici1:ici2,kzm1) :: akh , akm
    real(rkx) , dimension(jci1:jci2,ici1:ici2,kzm1) :: zint

    dtturbl = dt*stepbl
    rdtturbl = d_one/dtturbl
    dtdif = dtturbl

    do k = 1 , kzm1
      do i = ici1 , ici2
        do j = jci1 , jci2
          akm(j,i,k) = d_zero
        enddo
      enddo
    enddo

    do k = 1 , kzp1
      do i = ici1 , ici2
        do j = jci1 , jci2
          zint(j,i,k) = d_zero
        end do
      end do
    end do

    do i = ici1 , ici2
      do j = jci1 , jci2
        zint(j,i,kzp1) = ht(j,i) ! z at bottom of lowest sigma layer
      end do
    end do

    do k = kz , 1 , -1
      kflip = kzp1-k
      do i = ici1 , ici2
        do j = jci1 , jci2
          zint(j,i,k) = zint(j,i,k+1)+dz(j,i,kflip)
          apex = d_one/exner(j,i,k)
          ape(j,i,k) = apex
          tx = t(j,i,k)
          cwm(j,i,k) = qcw(j,i,k)
          if ( ipptls > 1 ) cwm(j,i,k) = cwm(j,i,k)+qci(j,i,k)
          the(j,i,k) = (cwm(j,i,k)*(-elocp/tx)+d_one)*th(j,i,k)
        end do
      end do
    end do

    el_myj(jci1:jci2,ici1:ici2,:) = d_zero

    setup_integration: &
    do i = ici1 , ici2
      do j = jci1 , jci2
        !
        ! Lowest layer above ground must be flipped
        !
        lmh = kz-lowlyr(j,i)+1
        ptop = pint(j,i,kzp1)
        psfc = pint(j,i,lowlyr(j,i))
        !
        ! Convert land mask (1 for sea; 0 for land)
        !
        seamask = xland(j,i)-d_one
        !
        ! Fill 1-d vertical arrays and flip direction since myj scheme
        ! counts downward from the domain's top
        !
        do k = kz , 1 , -1
          kflip = kzp1-k
          tk(k) = t(j,i,kflip)
          thek(k) = the(j,i,kflip)
          ratiomx = qv(j,i,kflip)
          qk(k) = ratiomx/(d_one+ratiomx)
          cwmk(k) = cwm(j,i,kflip)
          pk(k) = pmid(j,i,kflip)
          uk(k) = u(j,i,kflip)
          vk(k) = v(j,i,kflip)
          !
          ! TKE = 0.5*(q**2) ==> q**2 = 2.0*TKE
          !
          q2k(k) = d_two*tke_myj(j,i,kflip)
          !
          ! Compute the heights of the layer interfaces
          !
          zhk(k) = zint(j,i,k)
        end do
        zhk(kzp1) = ht(j,i) ! z at bottom of lowest sigma layer
        !
        ! Find the mixing length
        !
        call mixlen(lmh,uk,vk,tk,thek,qk,cwmk,q2k,zhk,gm,gh,el, &
                    pblh(j,i),lpbl(j,i),lmxl,ct(j,i),mixht(j,i))
        !
        ! Solve for the production/dissipation of the turbulent kinetic energy
        !
        call prodq2(lmh,dtturbl,ustar(j,i),gm,gh,el,q2k)
        !
        ! The model layer (counting upward) containing the top of the pbl
        !
        kpbl(j,i) = kz-lpbl(j,i)+1
        !
        ! Find the exchange coefficients in the free atmosphere
        !
        call difcof(lmh,lmxl,gm,gh,el,tk,q2k,zhk,akmk,akhk)
        !
        ! Counting downward from the top, the exchange coefficients akh
        ! are defined on the bottoms of the layers 1 to kzm1.
        ! Counting upward from the bottom, those same coefficients as
        ! exch_h are defined on the bottoms of the layers 1 to kzm1
        ! thus exch_h indices increase upward with k=1 at the ground.
        !
        do k = 1 , kzm1
          kflip = kz-k
          akh(j,i,k) = akhk(k)
          akm(j,i,k) = akmk(k)
          deltaz = d_half*(zhk(kflip)-zhk(kflip+2))
          exch_h(j,i,k) = akhk(kflip)*deltaz
        end do
        !
        ! Carry out the vertical diffusion of turbulent kinetic energy
        !
        call vdifq(lmh,dtdif,q2k,el,zhk)
        !
        ! Save the new TKE and mixing length.
        !
        do k = 1 , kz
          kflip = kzp1-k
          q2k(kflip) = max(q2k(kflip),epsq2)
          tke_myj(j,i,k) = d_half*q2k(kflip)
          ! el is not defined at kz (ground surface)
          if ( kflip < kz ) el_myj(j,i,k) = el(kflip)
        end do
      end do
    end do setup_integration
    !
    ! Convert surface sensible temperature to potential temperature.
    !
    do i = ici1 , ici2
      do j = jci1 , jci2
        psfc = pint(j,i,lowlyr(j,i))
        thsk(j,i) = tsk(j,i)*(p00/psfc)**rovcp
      end do
    end do

    main_integration: &
    do i = ici1 , ici2
      do j = jci1 , jci2
        !
        ! Fill 1-d vertical arrays and flip direction since myj scheme
        ! counts downward from the domain's top
        !
        do k = kz , 1 , -1
          kflip = kzp1-k
          thek(k) = the(j,i,kflip)
          ratiomx = qv(j,i,kflip)
          qk(k) = ratiomx/(d_one+ratiomx)
          cwmk(k) = cwm(j,i,kflip)
          qcwk(k) = qcw(j,i,kflip)
          species(1,k) = thek(k)
          species(2,k) = qk(k)
          species(3,k) = qcwk(k)
          nsp = 3
          if ( ipptls > 1 ) then
            nsp = nsp+1
            qcik(k) = qci(j,i,kflip)
            species(nsp,k) = qcik(k)
          end if
          zhk(k) = zint(j,i,k)
          rhok(k) = pmid(j,i,kflip)/(rgas*t(j,i,kflip)*  &
                    (d_one+ep1*qk(k)-cwmk(k)))
        end do
        !
        ! Counting downward from the top, the exchange coefficients akh
        ! are defined on the bottoms of the layers 1 to kz-1.
        ! These coefficients are also multiplied by the density at the
        ! bottom interface level.
        !
        do k = 1 , kzm1
          akhk(k) = akh(j,i,k)*d_half*(rhok(k)+rhok(k+1))
        end do

        zhk(kzp1) = zint(j,i,kzp1)

        seamask = xland(j,i)-d_one
        thz0(j,i) = (d_one-seamask)*thsk(j,i)+seamask*thz0(j,i)

        llow = lowlyr(j,i)
        akhs_dens = akhs(j,i)*rhok(kzp1-llow)

        if ( seamask < d_half ) then
          qfc1 = wlhv*chklowq(j,i)*akhs_dens
          if ( snow(j,i) > d_zero .or. sice(j,i) > 0.5_rkx ) then
            qfc1 = qfc1*rlivwv
          endif
          if ( qfc1 > d_zero ) then
            qlow = qk(kzp1-llow)
            qsfc(j,i) = qlow+elflx(j,i)/qfc1
          end if
        else
          psfc = pint(j,i,lowlyr(j,i))
          exnsfc = (p00/psfc)**rovcp
          qsfc(j,i) = pq0sea/psfc*exp(a2*(thsk(j,i)-a3*exnsfc) / &
            (thsk(j,i)-a4*exnsfc))
        end if

        qz0(j,i) = (d_one-seamask)*qsfc(j,i)+seamask*qz0(j,i)

        sz0(1) = thz0(j,i)
        sz0(2) = qz0(j,i)
        do nums = 3 , nspec
          sz0(nums) = d_zero
        end do

        clow(1) = d_one
        clow(2) = chklowq(j,i)
        do nums = 3 , nspec
          clow(nums) = d_zero
        end do

        cts(1) = ct(j,i)
        do nums = 2 , nspec
          cts(nums) = d_zero
        end do
        !
        ! Lowest layer above ground must be flipped
        !
        lmh = kz - lowlyr(j,i) + 1
        !
        ! Carry out the vertical diffusion of temperature and water vapor
        !
        call vdifh(dtdif,lmh,lpbl(j,i),sz0,akhs_dens,clow,cts, &
                   species,nspec,akhk,zhk,rhok(1))
        !
        ! Compute primary variable tendencies
        !
        do k = 1 , kz
          thek(k) = species(1,k)
          qk  (k) = species(2,k)
          qcwk(k) = species(3,k)
          cwmk(k) = qcwk(k)
          nsp = 3
          if ( ipptls > 1 ) then
            nsp = nsp+1
            qcik(k) = species(nsp,k)
            cwmk(k) = cwmk(k)+qcik(k)
          end if
        end do

        do k = 1 , kz
          kflip = kzp1-k
          thold = th(j,i,k)
          thnew = thek(kflip)+cwmk(kflip)*elocp*ape(j,i,k)
          dtdt = (thnew-thold)*rdtturbl
          qold = qv(j,i,k)/(d_one+qv(j,i,k))
          dqdt = (qk(kflip)-qold)*rdtturbl
          rthblten(j,i,k) = dtdt
          rqvblten(j,i,k) = dqdt/(d_one-qk(kflip))**2
          rqcblten(j,i,k) = (qcwk(kflip)-qcw(j,i,k))*rdtturbl
          if ( ipptls > 1 ) rqiblten(j,i,k) = (qcik(kflip)-qci(j,i,k))*rdtturbl
        end do

        psfc = 0.01_rkx*pint(j,i,lowlyr(j,i))

        !
        ! Fill 1-d vertical arrays and flip direction since myj scheme
        ! counts downward from the domain's top
        !
        do k = 1 , kz-1
          akmk(k) = akm(j,i,k)
          akmk(k) = akmk(k)*(rhok(k)+rhok(k+1))*d_half
        end do

        llow = lowlyr(j,i)
        akms_dens = akms(j,i)*rhok(kzp1-llow)
        do k = kz , 1 , -1
          kflip = kzp1-k
          uk(k) = u(j,i,kflip)
          vk(k) = v(j,i,kflip)
          zhk(k) = zint(j,i,k)
        end do
        zhk(kzp1) = zint(j,i,kzp1)
        !
        ! Carry out the vertical diffusion of velocity components
        !
        call vdifv(lmh,dtdif,uz0(j,i),vz0(j,i), &
                   akms_dens,uk,vk,akmk,zhk,rhok(1))
        !
        ! compute primary variable tendencies
        !
        do k = 1 , kz
          kflip = kzp1-k
          dudt = (uk(kflip)-u(j,i,k))*rdtturbl
          dvdt = (vk(kflip)-v(j,i,k))*rdtturbl
          rublten(j,i,k) = dudt
          rvblten(j,i,k) = dvdt
        end do
      end do
    end do main_integration
  end subroutine myjpbl
  !
  ! Level 2.5 Mixing Length
  !
  subroutine mixlen(lmh,u,v,t,the,q,cwm,q2,z,gm,gh,el,pblh,lpbl,lmxl,ct,mixht)
    implicit none
    integer(ik4) , intent(in) :: lmh
    integer(ik4) , intent(out) :: lmxl , lpbl
    real(rkx) , dimension(kz) , intent(in) :: cwm , q , q2 , t , the , u , v
    real(rkx) , dimension(kzp1) , intent(in) :: z
    real(rkx) , intent(out) :: pblh , mixht
    real(rkx) , dimension(kzm1) , intent(out) :: el , gh , gm
    real(rkx) , intent(inout) :: ct
    integer(ik4) :: k , lpblm
    real(rkx) :: a , aden , b , bden , aubr , bubr , blmx , el0 , &
      eloq2x , ghl , gml , qol2st , qol2un , qdzl , rdz , sq ,    &
      srel , szq , tem , thm , vkrmz
    real(rkx) , dimension(kz) :: q1
    real(rkx) , dimension(kzm1) :: dth , elm , rel
    logical :: lfound
    !
    ! Find the height of the pbl
    !
    lpbl = lmh
    lfound = .false.
    do k = lmh-1 , 1 , -1
      if ( q2(k) <= epsq2*fh ) then
        lpbl = k
        lfound = .true.
        pblh = z(lpbl+1)-z(lmh+1)
      end if
    end do
    lpbl = 1
    if ( .not. lfound ) pblh = z(lpbl+1)-z(lmh+1)
    !
    ! The height of the pbl
    !
 110  pblh = z(lpbl+1)-z(lmh+1)

    do k = 1 , lmh
      q1(k) = d_zero
    end do
    do k = 1 , lmh-1
      dth(k) = the(k)-the(k+1)
    enddo
    do k = lmh-2 , 1 , -1
      if ( dth(k) > d_zero .and. dth(k+1) <= d_zero ) then
        dth(k) = dth(k)+ct
        exit
      end if
    end do
    ct = d_zero
    do k = 1 , lmh-1
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
    lmxl = lmh
    do k = 1 , lmh-1
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
          eloq2x = d_one/qol2st
          elm(k) = max(sqrt(eloq2x*q2(k)),epsl)
        endif
      else
        aden   = (adnm*gml+adnh*ghl)*ghl
        bden   = bdnm*gml+bdnh*ghl
        qol2un = -d_half*bden+sqrt(bden*bden*0.25_rkx-aden)
        eloq2x = d_one/(qol2un+epsru)
        elm(k) = max(sqrt(eloq2x*q2(k)),epsl)
      end if
    end do
    if ( elm(lmh-1) == epsl ) lmxl = lmh
    !
    ! The height of the mixed layer
    !
    blmx  = z(lmxl)-z(lmh+1)
    mixht = blmx
    do k = lpbl , lmh
      q1(k) = sqrt(q2(k))
    enddo
    szq = d_zero
    sq  = d_zero
    do k = 1 , lmh-1
      qdzl = (q1(k)+q1(k+1))*(z(k+1)-z(k+2))
      szq  = (z(k+1)+z(k+2)-z(lmh+1)-z(lmh+1))*qdzl+szq
      sq   = qdzl+sq
    end do
    !
    ! Computation of asymptotic l in blackadar formula
    !
    el0 = min(alph*szq*d_half/sq,el0max)
    el0 = max(el0,el0min)
    !
    ! Above the pbl top
    !
    lpblm = max(lpbl-1,1)
    do k = 1 , lpblm
      el(k)  = min((z(k)-z(k+2))*elfc,elm(k))
      rel(k) = el(k)/elm(k)
    end do
    !
    ! Inside the pbl
    !
    if ( lpbl < lmh ) then
      do k = lpbl , lmh-1
        vkrmz = (z(k+1)-z(lmh+1))*vonkar
        el(k) = min(vkrmz/(vkrmz/el0+d_one),elm(k))
        rel(k) = el(k)/elm(k)
      end do
    end if
    do k = lpbl+1 , lmh-2
      srel  = min(((rel(k-1)+rel(k+1))*d_half+rel(k))*d_half,rel(k))
      el(k) = max(srel*elm(k),epsl)
    end do
  end subroutine mixlen
  !
  ! Level 2.5 Q2 Production/Dissipation
  !
  subroutine prodq2(lmh,dtturbl,ustar,gm,gh,el,q2)
    implicit none
    integer(ik4) , intent(in) :: lmh
    real(rkx) , intent(in) :: dtturbl , ustar
    real(rkx) , dimension(kzm1) , intent(in) :: gh , gm
    real(rkx) , dimension(kzm1) , intent(inout) :: el
    real(rkx) , dimension(kz) , intent(inout) :: q2
    integer(ik4) :: k
    real(rkx) :: aden , aequ , anum , arhs , bden , bequ , bnum , brhs , &
      cden , crhs , dloq1 , eloq11 , eloq12 , eloq13 , eloq21 , eloq22 , &
      eloq31 , eloq32 , eloq41 , eloq42 , eloq51 , eloq52 , eloqn ,      &
      eqol2 , ghl , gml , rden1 , rden2 , rhs2 , rhsp1 , rhsp2 , rhst2

    main_integration: &
    do k = 1 , lmh-1
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
    q2(lmh) = max(b1**(d_two/d_three)*ustar*ustar,epsq2)
  end subroutine prodq2
  !
  ! Level 2.5 Diffusion Coefficients
  !
  subroutine difcof(lmh,lmxl,gm,gh,el,t,q2,z,akm,akh)
    implicit none
    integer(ik4) , intent(in) :: lmh , lmxl
    real(rkx) , dimension(kz) , intent(in) :: q2 , t
    real(rkx) , dimension(kzm1) , intent(in) :: el , gh , gm
    real(rkx) , dimension(kzp1) , intent(in) :: z
    real(rkx) , dimension(kzm1) , intent(out) :: akh , akm
    integer(ik4) :: k , kinv
    real(rkx) :: aden , akmin , bden , besh , besm , cden , d2t , &
      ell , eloq2 , eloq4 , elqdz , esh , esm , ghl , gml , q1l , &
      rden , rdz

    do k = 1 , lmh-1
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
  subroutine vdifq(lmh,dtdif,q2,el,z)
    implicit none
    integer(ik4) , intent(in) :: lmh
    real(rkx) ,intent(in) :: dtdif
    real(rkx) ,dimension(kzm1),intent(in) :: el
    real(rkx) ,dimension(kzp1),intent(in) :: z

    real(rkx) ,dimension(kz),intent(inout) :: q2

    integer(ik4) :: k
    real(rkx) :: akqs , cf , dtozs , esqhf , rdz
    real(rkx) , dimension(kzm2) :: akq , cm , cr , dtoz , rsq2

    esqhf = d_half*esq

    do k = 1 , lmh-2
      dtoz(k) = (dtdif+dtdif)/(z(k)-z(k+2))
      akq(k)  = sqrt((q2(k)+q2(k+1))*d_half)*(el(k)+el(k+1))* &
                esqhf/(z(k+1)-z(k+2))
      cr(k)   = -dtoz(k)*akq(k)
    end do

    cm(1)   = dtoz(1)*akq(1)+d_one
    rsq2(1) = q2(1)

    do k = 2 , lmh-2
      cf      = -dtoz(k)*akq(k-1)/cm(k-1)
      cm(k)   = -cr(k-1)*cf+(akq(k-1)+akq(k))*dtoz(k)+d_one
      rsq2(k) = -rsq2(k-1)*cf+q2(k)
    end do

    dtozs = (dtdif+dtdif)/(z(lmh-1)-z(lmh+1))
    akqs = sqrt((q2(lmh-1)+q2(lmh))*d_half)*(el(lmh-1)+elz0) * &
           esqhf/(z(lmh)-z(lmh+1))
    cf = -dtozs*akq(lmh-2)/cm(lmh-2)

    q2(lmh-1) = (dtozs*akqs*q2(lmh)-rsq2(lmh-2)*cf+q2(lmh-1)) / &
                ((akq(lmh-2)+akqs)*dtozs-cr(lmh-2)*cf+d_one)

    do k = lmh-2 , 1 , -1
      q2(k) = (-cr(k)*q2(k+1)+rsq2(k))/cm(k)
    end do
  end subroutine vdifq
  !
  ! Vertical diffusion of mass variables
  !
  subroutine vdifh(dtdif,lmh,lpbl,sz0,rkhs,clow,cts, &
                   species,nspec,rkh,zhk,rho)
    implicit none
    integer(ik4) , intent(in) :: lmh , lpbl , nspec
    real(rkx) , intent(in) :: dtdif , rkhs
    real(rkx) , dimension(nspec) , intent(in) :: clow , cts , sz0
    real(rkx) , dimension(kz-1) , intent(in) :: rkh
    real(rkx) , dimension(kz) , intent(in) :: rho
    real(rkx) , dimension(kzp1) , intent(in) :: zhk
    real(rkx) , dimension(nspec,kz) , intent(inout) :: species

    integer(ik4) :: k , m
    real(rkx) :: cf , cmb , cmsb , dtozl , dtozs , rcml , rhok ,&
      rkhh , rkhz , rkss , rssb
    real(rkx) , dimension(kzm1) :: cm , cr , dtoz
    real(rkx) , dimension(nspec,kzm1) :: rkct , rss

    do k = 1 , lmh-1
      dtoz(k) = dtdif/(zhk(k)-zhk(k+1))
      cr(k) = -dtoz(k)*rkh(k)
      if ( k < lpbl ) then
        do m = 1 , nspec
          rkct(m,k) = d_zero
        end do
      else
        rkhz = rkh(k)*(zhk(k)-zhk(k+2))
        do m = 1 , nspec
          rkct(m,k) = d_half*rkhz*cts(m)
        end do
      end if
    end do
    ! Top level
    rhok  = rho(1)
    cm(1) = dtoz(1)*rkh(1)+rhok
    do m = 1 , nspec
      rss(m,1) = -rkct(m,1)*dtoz(1)+species(m,1)*rhok
    end do
    ! Intermediate levels
    do k = 2 , lmh-1
      dtozl = dtoz(k)
      cf    = -dtozl*rkh(k-1)/cm(k-1)
      rhok  = rho(k)
      cm(k) = -cr(k-1)*cf+(rkh(k-1)+rkh(k))*dtozl+rhok
      do m = 1 , nspec
        rss(m,k) = -rss(m,k-1)*cf + &
                   (rkct(m,k-1)-rkct(m,k))*dtozl+species(m,k)*rhok
      end do
    end do
    ! Bottom level
    dtozs = dtdif/(zhk(lmh)-zhk(lmh+1))
    rkhh  = rkh(lmh-1)
    cf    = -dtozs*rkhh/cm(lmh-1)
    cmb   = cr(lmh-1)*cf
    rhok  = rho(lmh)
    do m = 1 , nspec
      rkss = rkhs*clow(m)
      cmsb = -cmb+(rkhh+rkss)*dtozs+rhok
      rssb = -rss(m,lmh-1)*cf+rkct(m,lmh-1)*dtozs+species(m,lmh)*rhok
      species(m,lmh) = (dtozs*rkss*sz0(m)+rssb)/cmsb
    end do
    ! Backsubstitution
    do k = lmh-1 , 1 , -1
      rcml = d_one/cm(k)
      do m = 1 , nspec
        species(m,k) = (-cr(k)*species(m,k+1)+rss(m,k))*rcml
      end do
    end do
  end subroutine vdifh
  !
  ! Vertical diffusion of mass variables
  !
  subroutine vdifx(dtdif,lmh,rkhs,ct,dust1,dust2,rkh,z,rho)
    implicit none
    integer(ik4) , intent(in) :: lmh
    real(rkx) , intent(in) :: ct,dtdif,rkhs
    real(rkx) , dimension(kzm1) , intent(in) :: rkh
    real(rkx) , dimension(kz) , intent(in) :: rho
    real(rkx) , dimension(kzp1) , intent(in) :: z
    real(rkx) , dimension(kz) , intent(inout) :: dust1 , dust2
    integer(ik4) :: k
    real(rkx) :: cf , cmb , cmdb , cthf , dtozl , dtozs , rcml , &
      rkhh , rsd1b , rsd2b
    real(rkx) , dimension(kzm1) :: cm , cr , dtoz , rkct , rsd1 , rsd2

    cthf = d_half*ct

    do k = 1 , lmh-1
      dtoz(k) = dtdif/(z(k)-z(k+1))
      cr(k)   = -dtoz(k)*rkh(k)
      rkct(k) = rkh(k)*(z(k)-z(k+2))*cthf
    enddo

    cm(1) = dtoz(1)*rkh(1)+rho(1)

    rsd1(1) = dust1(1)*rho(1)
    rsd2(1) = dust2(1)*rho(1)

    do k = 2 , lmh-1
      dtozl   = dtoz(k)
      cf      = -dtozl*rkh(k-1)/cm(k-1)
      cm(k)   = -cr(k-1)*cf+(rkh(k-1)+rkh(k))*dtozl+rho(k)
      rsd1(k) = -rsd1(k-1)*cf+dust1(k)*rho(k)
      rsd2(k) = -rsd2(k-1)*cf+dust2(k)*rho(k)
    end do

    dtozs = dtdif/(z(lmh)-z(lmh+1))
    rkhh = rkh(lmh-1)

    cf = -dtozs*rkhh/cm(lmh-1)

    cmb  = cr(lmh-1)*cf
    cmdb = -cmb+rkhh*dtozs+rho(lmh)

    rsd1b = -rsd1(lmh-1)*cf+dust1(lmh)*rho(lmh)
    rsd2b = -rsd2(lmh-1)*cf+dust2(lmh)*rho(lmh)

    dust1(lmh) = rsd1b/cmdb
    dust2(lmh) = rsd2b/cmdb

    do k = lmh-1 , 1, -1
      rcml     = d_one/cm(k)
      dust1(k) = (-cr(k)*dust1(k+1)+rsd1(k))*rcml
      dust2(k) = (-cr(k)*dust2(k+1)+rsd2(k))*rcml
    end do
  end subroutine vdifx
  !
  ! Vertical diffusion of velocity components
  !
  subroutine vdifv(lmh,dtdif,uz0,vz0,rkms,u,v,rkm,z,rho)
    implicit none
    integer(ik4) , intent(in) :: lmh
    real(rkx) , intent(in) :: rkms , dtdif , uz0 , vz0
    real(rkx) , dimension(kzm1) , intent(in) :: rkm
    real(rkx) , dimension(kz) , intent(in) :: rho
    real(rkx) , dimension(kzp1) , intent(in) :: z
    real(rkx) , dimension(kz) , intent(inout) :: u , v
    integer(ik4) :: k
    real(rkx) :: cf , dtozak , dtozl , dtozs , rcml , rcmvb , rhok , rkmh
    real(rkx) , dimension(kzm1) :: cm , cr , dtoz , rsu , rsv

    do k = 1 , lmh-1
      dtoz(k) = dtdif/(z(k)-z(k+1))
      cr(k)   = -dtoz(k)*rkm(k)
    end do

    rhok   = rho(1)
    cm(1)  = dtoz(1)*rkm(1)+rhok
    rsu(1) = u(1)*rhok
    rsv(1) = v(1)*rhok

    do k = 2 , lmh-1
      dtozl  = dtoz(k)
      cf     = -dtozl*rkm(k-1)/cm(k-1)
      rhok   = rho(k)
      cm(k)  = -cr(k-1)*cf+(rkm(k-1)+rkm(k))*dtozl+rhok
      rsu(k) = -rsu(k-1)*cf+u(k)*rhok
      rsv(k) = -rsv(k-1)*cf+v(k)*rhok
    end do

    dtozs  = dtdif/(z(lmh)-z(lmh+1))
    rkmh   = rkm(lmh-1)
    cf     = -dtozs*rkmh/cm(lmh-1)
    rhok   = rho(lmh)
    rcmvb  = d_one/((rkmh+rkms)*dtozs-cr(lmh-1)*cf+rhok)
    dtozak = dtozs*rkms
    u(lmh) = (dtozak*uz0-rsu(lmh-1)*cf+u(lmh)*rhok)*rcmvb
    v(lmh) = (dtozak*vz0-rsv(lmh-1)*cf+v(lmh)*rhok)*rcmvb

    do k = lmh-1 , 1 , -1
      rcml = d_one/cm(k)
      u(k) = (-cr(k)*u(k+1)+rsu(k))*rcml
      v(k) = (-cr(k)*v(k+1)+rsv(k))*rcml
    end do
  end subroutine vdifv

end module mod_pbl_myj

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

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

module mod_split
  !
  ! Split explicit time integration
  !
  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_mppparam
  use mod_runparams
  use mod_constants
  use mod_stdio
  use mod_atm_interface
  use mod_vmodes
  use mod_bdycod
  use mod_atm_interface
  use mod_memutil
  use mod_service

  implicit none

  private

  public :: allocate_mod_split, spinit, splitf

  real(rkx), pointer, contiguous, dimension(:) :: aam => null( )
  real(rkx), pointer, contiguous, dimension(:) :: an => null( )
  real(rkx), pointer, contiguous, dimension(:,:) :: am => null( )
  real(rkx), pointer, contiguous, dimension(:,:) :: map => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: uuu => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: vvv => null( )

  real(rkx), pointer, contiguous, dimension(:,:,:) :: ddsum => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: dhsum => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:,:) :: deld => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:,:) :: delh => null( )
  real(rkx), pointer, contiguous, dimension(:,:,:) :: work => null( )
  real(rkx), pointer, contiguous, dimension(:,:) :: uu => null( )
  real(rkx), pointer, contiguous, dimension(:,:) :: vv => null( )
  real(rkx), pointer, contiguous, dimension(:,:) :: xdelh => null( )

  contains

  subroutine allocate_mod_split
    implicit none
    call getmem1d(aam,1,nsplit,'split:aam')
    call getmem2d(am,1,kz,1,nsplit,'split:am')
    call getmem2d(map,jce1,jce2,ice1,ice2,'split:map')
    call getmem1d(an,1,nsplit,'split:naam')
    call getmem3d(ddsum,jde1,jde2,ide1,ide2,1,nsplit,'split:ddsum')
    call getmem4d(deld,jde1,jde2,ide1,ide2,1,nsplit,1,3,'split:deld')
    call getmem4d(delh,jde1,jde2,ide1,ide2,1,nsplit,1,3,'split:delh')
    call getmem2d(xdelh,jde1ga,jde2,ide1ga,ide2,'split:xdelh')
    call getmem3d(dhsum,jde1ga,jde2,ide1ga,ide2,1,nsplit,'split:dhsum')
    call getmem3d(work,jdi1,jdi2,idi1,idi2,1,3,'split:work')
    call getmem2d(uu,jdi1,jdi2ga,idi1,idi2ga,'split:uu')
    call getmem2d(vv,jdi1,jdi2ga,idi1,idi2ga,'split:vv')
    call getmem3d(uuu,jde1,jde2ga,ide1,ide2ga,1,kz,'split:uuu')
    call getmem3d(vvv,jde1,jde2ga,ide1,ide2ga,1,kz,'split:vvv')
  end subroutine allocate_mod_split
  !
  ! Intial computation of vertical modes.
  !
  subroutine spinit
    implicit none
    real(rkx) :: eps1, fac, pdlog
    integer(ik4) :: i, j, k, l, n, ns
    real(rkx) :: rnpts, lxps, ltbark, lxms, xmsf, rdx2
    real(rkx) :: eps
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'spinit'
    integer(ik4), save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    ! Compute m.
    !
    rdx2 = d_one/dx2
    do ns = 1, nsplit
      aam(ns) = nint(dtsec/dtau(ns))
    end do
    !
    ! lstand = .true. if standard atmosphere t to be used (ignore input
    ! tbarh and xps in that case). otherwise, xps and tbarh must
    ! be defined on input.
    !
    call allocate_mod_vmodes
    !
    do concurrent ( j = jce1:jce2, i = ice1:ice2 )
      map(j,i) = d_one / (mddom%msfx(j,i)*mddom%msfx(j,i))
    end do
    !
    if ( .not. lstand ) then
      !
      ! compute xps and tbarh for use in vmodes.
      !
      rnpts = d_one/real((nicross*njcross),rkx)
      lxps = d_zero
      lxms = d_zero
      do i = ice1, ice2
        do j = jce1, jce2
          lxms = lxms + rnpts * map(j,i)
          lxps = lxps + rnpts * sfs%psa(j,i) * map(j,i)
        end do
      end do
      call sumall(lxps,xps)
      call sumall(lxms,xmsf)
      xps = xps / xmsf
      do k = 1, kz
        ltbark = d_zero
        do i = ice1, ice2
          do j = jce1, jce2
            ltbark = ltbark + rnpts * map(j,i) * atm1%t(j,i,k)/sfs%psa(j,i)
          end do
        end do
        call sumall(ltbark,tbarh(k))
        tbarh(k) = tbarh(k)/xmsf
      end do
    end if
    !
    ! Compute vertical modes.
    !
    call vmodes
    !
    ! Compute am and an.
    !
    do n = 1, nsplit
      an(n) = d_zero
      do l = 1, kz
        an(n) = an(n) + dsigma(l)*zmatx(l,n)
      end do
      do k = 1, kz
        am(k,n) = d_zero
        tau(n,k) = d_zero
      end do
      do l = 1, kz
        do k = 1, kz
          am(k,n) = am(k,n) + a0(k,l)*zmatx(l,n)
          tau(n,k) = tau(n,k) + rgas*zmatxr(n,l)*hydros(l,k)
        end do
      end do
      do k = 1, kzp1
        varpa1(n,k) = d_zero
      end do
      do l = 1, kz
        do k = 1, kzp1
          varpa1(n,k) = varpa1(n,k) + rgas*zmatxr(n,l)*hydroc(l,k)
        end do
      end do
    end do
    !
    ! Multiply am, an and zmatx by factor.
    !
    do l = 1, nsplit
      fac = d_two*dtsec/(d_two*aam(l)+d_one)
      an(l) = an(l)*fac
      do k = 1, kz
        zmatx(k,l) = zmatx(k,l)*fac
        am(k,l) = am(k,l)*fac
      end do
      if ( myid == italk ) then
        write(stdout,'(a,i4,a,f11.4,a,f11.4)') &
          ' Split : ',l,' => aam :',aam(l),', fac :',fac
      end if
    end do
    !
    ! If a restart run, do not recalculate the hstor/dstor
    !
    if ( ifrest ) then
#ifdef DEBUG
      call time_end(subroutine_name,idindx)
#endif
      return
    end if
    !
    ! Divergence manipulations (0)
    !
    ! compute divergence z from u and v
    ! ( u must be pstar * u ; similarly for v )
    ! ( note: map scale factors have been inverted in model (init) )
    !
    do concurrent ( j = jde1:jde2, i = ide1:ide2, k = 1:kz )
      uuu(j,i,k) = atm2%u(j,i,k)*mddom%msfd(j,i)
      vvv(j,i,k) = atm2%v(j,i,k)*mddom%msfd(j,i)
    end do

    call exchange_rt(uuu,1,jde1,jde2,ide1,ide2,1,kz)
    call exchange_rt(vvv,1,jde1,jde2,ide1,ide2,1,kz)

    do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz, l = 1:nsplit )
      dstor(j,i,l) = dstor(j,i,l) + zmatxr(l,k) * map(j,i) * rdx2 * &
                 (-uuu(j,i+1,k)+uuu(j+1,i+1,k)-uuu(j,i,k)+uuu(j+1,i,k) + &
                   vvv(j,i+1,k)+vvv(j+1,i+1,k)-vvv(j,i,k)-vvv(j+1,i,k))
    end do
    !
    ! Geopotential manipulations
    !
    do l = 1, nsplit
      pdlog = varpa1(l,kzp1)*log(sigmah(kzp1)*pd+ptop)
      eps1 = varpa1(l,kzp1)*sigmah(kzp1)/(sigmah(kzp1)*pd+ptop)
      do concurrent ( j = jce1:jce2, i = ice1:ice2 )
        eps = eps1*(sfs%psb(j,i)-pd)
        hstor(j,i,l) = pdlog + eps
      end do

      do k = 1, kz
        pdlog = varpa1(l,k)*log(sigmah(k)*pd+ptop)
        eps1 = varpa1(l,k)*sigmah(k)/(sigmah(k)*pd+ptop)
        do concurrent ( j = jce1:jce2, i = ice1:ice2 )
          eps = eps1*(sfs%psb(j,i)-pd)
          hstor(j,i,l) = hstor(j,i,l) + pdlog + &
                         tau(l,k)*atm2%t(j,i,k)/sfs%psb(j,i) + eps
        end do
      end do
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine spinit
  !
  ! Compute deld, delh, integrate in time and add correction terms appropriately
  !
  subroutine splitf
    implicit none
    real(rkx) :: rdx2, eps1, gnuam, gnuan, gnuzm, pdlog
    real(rkx) :: eps, fac, x, y
    integer(ik4) :: i, j, k, l, n
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'splitf'
    integer(ik4), save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    rdx2 = d_one/dx2
    do concurrent ( j = jde1:jde2, i = ide1:ide2, n = 1:nsplit, l = 1:3 )
      deld(j,i,n,l) = d_zero
      delh(j,i,n,l) = d_zero
    end do
    !
    ! compute pressure on dot grid
    !
    call exchange(sfs%psa,1,jce1,jce2,ice1,ice2)
    call psc2psd(sfs%psa,sfs%psdota)
    !
    ! get deld(0), delh(0) from storage
    !
    do concurrent ( j = jde1:jde2, i = ide1:ide2, n = 1:nsplit )
      deld(j,i,n,1) = dstor(j,i,n)
      delh(j,i,n,1) = hstor(j,i,n)
    end do
    !
    ! Divergence manipulations (f)
    !
    do concurrent ( j = jde1:jde2, i = ide1:ide2, k = 1:kz )
      uuu(j,i,k) = atm1%u(j,i,k)*mddom%msfd(j,i)
      vvv(j,i,k) = atm1%v(j,i,k)*mddom%msfd(j,i)
    end do

    call exchange_rt(uuu,1,jde1,jde2,ide1,ide2,1,kz)
    call exchange_rt(vvv,1,jde1,jde2,ide1,ide2,1,kz)

    do l = 1, nsplit
      do concurrent ( j = jde1:jde2, i = ide1:ide2 )
        deld(j,i,l,3) = d_zero
      end do

      do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
        deld(j,i,l,3) = deld(j,i,l,3) + zmatxr(l,k) * rdx2 * map(j,i) * &
               (-uuu(j,i+1,k)+uuu(j+1,i+1,k)-uuu(j,i,k)+uuu(j+1,i,k) + &
                 vvv(j,i+1,k)+vvv(j+1,i+1,k)-vvv(j,i,k)-vvv(j+1,i,k))
      end do
    end do

    do concurrent ( j = jde1:jde2, i = ide1:ide2, n = 1:nsplit )
      deld(j,i,n,3) = deld(j,i,n,3) - deld(j,i,n,1)
    end do
    !
    ! Divergence manipulations (0)
    !
    do concurrent ( j = jde1:jde2, i = ide1:ide2, k = 1:kz )
      uuu(j,i,k) = atm2%u(j,i,k)*mddom%msfd(j,i)
      vvv(j,i,k) = atm2%v(j,i,k)*mddom%msfd(j,i)
    end do

    call exchange_rt(uuu,1,jde1,jde2,ide1,ide2,1,kz)
    call exchange_rt(vvv,1,jde1,jde2,ide1,ide2,1,kz)

    do l = 1, nsplit
      do concurrent ( j = jde1:jde2, i = ide1:ide2 )
        deld(j,i,l,2) = d_zero
      end do
      do concurrent ( j = jce1:jce2, i = ice1:ice2, k = 1:kz )
        deld(j,i,l,2) = deld(j,i,l,2) + zmatxr(l,k) * rdx2 * map(j,i) * &
              (-uuu(j,i+1,k)+uuu(j+1,i+1,k)-uuu(j,i,k)+uuu(j+1,i,k) + &
                vvv(j,i+1,k)+vvv(j+1,i+1,k)-vvv(j,i,k)-vvv(j+1,i,k))
      end do
    end do

    do concurrent ( j = jde1:jde2, i = ide1:ide2, n = 1:nsplit )
      deld(j,i,n,1) = deld(j,i,n,1) - deld(j,i,n,2)
    end do
    !
    ! Geopotential manipulations (f)
    !
    do l = 1, nsplit
      pdlog = varpa1(l,kzp1)*log(sigmah(kzp1)*pd+ptop)
      eps1 = varpa1(l,kzp1)*sigmah(kzp1)/(sigmah(kzp1)*pd+ptop)
      do concurrent ( j = jce1:jce2, i = ice1:ice2 )
        eps = eps1*(sfs%psa(j,i)-pd)
        delh(j,i,l,3) = pdlog + eps
      end do
      do k = 1, kz
        pdlog = varpa1(l,k)*log(sigmah(k)*pd+ptop)
        eps1 = varpa1(l,k)*sigmah(k)/(sigmah(k)*pd+ptop)
        do concurrent ( j = jce1:jce2, i = ice1:ice2 )
          eps = eps1*(sfs%psa(j,i)-pd)
          delh(j,i,l,3) = delh(j,i,l,3) + pdlog +  &
              tau(l,k)*atm1%t(j,i,k)/sfs%psa(j,i) + eps
        end do
      end do
    end do

    do concurrent ( j = jde1:jde2, i = ide1:ide2, n = 1:nsplit )
      delh(j,i,n,3) = delh(j,i,n,3) - delh(j,i,n,1)
    end do
    !
    ! Geopotential manipulations (0)
    !
    do l = 1, nsplit
      pdlog = varpa1(l,kzp1)*log(sigmah(kzp1)*pd+ptop)
      eps1 = varpa1(l,kzp1)*sigmah(kzp1)/(sigmah(kzp1)*pd+ptop)
      do concurrent ( j = jce1:jce2, i = ice1:ice2 )
        eps = eps1*(sfs%psb(j,i)-pd)
        delh(j,i,l,2) = pdlog + eps
      end do
      do k = 1, kz
        pdlog = varpa1(l,k)*log(sigmah(k)*pd+ptop)
        eps1 = varpa1(l,k)*sigmah(k)/(sigmah(k)*pd+ptop)
        do concurrent ( j = jce1:jce2, i = ice1:ice2 )
          eps = eps1*(sfs%psb(j,i)-pd)
          delh(j,i,l,2) = delh(j,i,l,2) + pdlog +  &
                   tau(l,k)*atm2%t(j,i,k)/sfs%psb(j,i) + eps
        end do
      end do
    end do

    do concurrent ( j = jde1:jde2, i = ide1:ide2, n = 1:nsplit )
      delh(j,i,n,1) = delh(j,i,n,1) - delh(j,i,n,2)
    end do
    !
    ! put deld(0), delh(0) into storage
    !
    do concurrent ( j = jde1:jde2, i = ide1:ide2, n = 1:nsplit )
      dstor(j,i,n) = deld(j,i,n,2)
      hstor(j,i,n) = delh(j,i,n,2)
    end do
    !
    ! split explicit time integration
    !
    call spstep
    !
    ! Add corrections to t and p;  u and v
    !
    do l = 1, nsplit
      gnuan = gnu1*an(l)
      do concurrent ( j = jci1:jci2, i = ici1:ici2 )
        sfs%psa(j,i) = sfs%psa(j,i) - an(l)*ddsum(j,i,l)
        sfs%psb(j,i) = sfs%psb(j,i) - gnuan*ddsum(j,i,l)
      end do
    end do
    do l = 1, nsplit
      do k = 1, kz
        gnuam = gnu1*am(k,l)
        do concurrent ( j = jci1:jci2, i = ici1:ici2 )
          atm1%t(j,i,k) = atm1%t(j,i,k) + am(k,l)*ddsum(j,i,l)
          atm2%t(j,i,k) = atm2%t(j,i,k) + gnuam*ddsum(j,i,l)
        end do
      end do
    end do

    call exchange_lb(dhsum,1,jde1,jde2,ide1,ide2,1,nsplit)

    do l = 1, nsplit
      do k = 1, kz
        gnuzm = gnu1*zmatx(k,l)
        do concurrent ( j = jdi1:jdi2, i = idi1:idi2 )
          fac = sfs%psdota(j,i)/(dx2*mddom%msfd(j,i))
          x = fac*(dhsum(j,i,l)+dhsum(j,i-1,l) - &
                   dhsum(j-1,i,l)-dhsum(j-1,i-1,l))
          y = fac*(dhsum(j,i,l)-dhsum(j,i-1,l) + &
                   dhsum(j-1,i,l)-dhsum(j-1,i-1,l))
          atm1%u(j,i,k) = atm1%u(j,i,k) - zmatx(k,l)*x
          atm1%v(j,i,k) = atm1%v(j,i,k) - zmatx(k,l)*y
          atm2%u(j,i,k) = atm2%u(j,i,k) - gnuzm*x
          atm2%v(j,i,k) = atm2%v(j,i,k) - gnuzm*y
        end do
      end do
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine splitf

  subroutine spstep
    implicit none
    real(rkx) :: rdx2, dtau2, fac
    integer(ik4) :: i, j, m2, n, n0, n1, n2, ns, nw
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'spstep'
    integer(ik4), save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    rdx2 = d_one/dx2

    do concurrent ( j = jde1:jde2, i = ide1:ide2, n = 1:nsplit )
      ddsum(j,i,n) = d_zero
      dhsum(j,i,n) = d_zero
    end do

    do ns = 1, nsplit
      n0 = 1
      n1 = 2
      n2 = n0
      m2 = int(aam(ns))*2
      dtau2 = dtau(ns)*d_two
      !
      ! below follows Madala (1987)
      !
      do concurrent ( j = jce1:jce2, i = ice1:ice2 )
        ddsum(j,i,ns) = deld(j,i,ns,n0)
        dhsum(j,i,ns) = delh(j,i,ns,n0)
      end do
      !
      ! first step, use forward scheme
      !
      ! compute gradient of delh;  output = (work1,work2)
      !
      xdelh(jde1:jde2,ide1:ide2) = delh(jde1:jde2,ide1:ide2,ns,n0)
      call exchange_lb(xdelh,1,jde1,jde2,ide1,ide2)
      do concurrent ( j = jdi1:jdi2, i = idi1:idi2 )
        fac = dx2*mddom%msfx(j,i)
        work(j,i,1) = (xdelh(j,i)  +xdelh(j,i-1) - &
                       xdelh(j-1,i)-xdelh(j-1,i-1))/fac
        work(j,i,2) = (xdelh(j,i)  +xdelh(j-1,i) - &
                       xdelh(j,i-1)-xdelh(j-1,i-1))/fac
      end do

      do concurrent ( j = jdi1:jdi2, i = idi1:idi2, nw = 1:2 )
        work(j,i,nw) = work(j,i,nw)*sfs%psdota(j,i)
      end do
      !
      ! compute divergence z from u and v
      ! ( u must be pstar * u ; similarly for v )
      ! ( note: map scale factors have been inverted in model (init) )
      !
      do concurrent ( j = jdi1:jdi2, i = idi1:idi2 )
        uu(j,i) = work(j,i,1)*mddom%msfd(j,i)
        vv(j,i) = work(j,i,2)*mddom%msfd(j,i)
      end do

      call exchange_rt(uu,1,jdi1,jdi2,idi1,idi2)
      call exchange_rt(vv,1,jdi1,jdi2,idi1,idi2)

      do concurrent ( j = jci1:jci2, i = ici1:ici2 )
        work(j,i,3) = rdx2 * map(j,i) * &
                (-uu(j,i+1)+uu(j+1,i+1)-uu(j,i)+uu(j+1,i) &
                 +vv(j,i+1)+vv(j+1,i+1)-vv(j,i)-vv(j+1,i))
      end do

      do concurrent ( j = jci1:jci2, i = ici1:ici2 )
        deld(j,i,ns,n1) = deld(j,i,ns,n0) - dtau(ns)*work(j,i,3) + &
                          deld(j,i,ns,3)/m2
        delh(j,i,ns,n1) = delh(j,i,ns,n0) - dtau(ns)*hbar(ns) * &
                          deld(j,i,ns,n0)/sfs%psa(j,i)+delh(j,i,ns,3)/m2
      end do
      !
      ! not in Madala (1987)
      !
      fac = (aam(ns)-d_one)/aam(ns)
      if ( ma%has_bdyleft ) then
        do i = ici1, ici2
          delh(jce1,i,ns,n1) = delh(jce1,i,ns,n0)*fac
        end do
      end if
      if ( ma%has_bdyright ) then
        do i = ici1, ici2
          delh(jce2,i,ns,n1) = delh(jce2,i,ns,n0)*fac
        end do
      end if
      if ( ma%has_bdybottom ) then
        do j = jce1, jce2
          delh(j,ice1,ns,n1) = delh(j,ice1,ns,n0)*fac
        end do
      end if
      if ( ma%has_bdytop ) then
        do j = jce1, jce2
          delh(j,ice2,ns,n1) = delh(j,ice2,ns,n0)*fac
        end do
      end if

      do concurrent ( j = jce1:jce2, i = ice1:ice2 )
        ddsum(j,i,ns) = ddsum(j,i,ns) + deld(j,i,ns,n1)
        dhsum(j,i,ns) = dhsum(j,i,ns) + delh(j,i,ns,n1)
      end do
      !
      ! subsequent steps, use leapfrog scheme
      !
      do n = 2, m2
        !
        ! compute gradient of delh;  output = (work1,work2)
        !
        xdelh(jde1:jde2,ide1:ide2) = delh(jde1:jde2,ide1:ide2,ns,n1)
        call exchange_lb(xdelh,1,jde1,jde2,ide1,ide2)
        do concurrent ( j = jdi1:jdi2, i = idi1:idi2 )
          fac = dx2*mddom%msfx(j,i)
          work(j,i,1) = (xdelh(j,i)+xdelh(j,i-1)- &
                         xdelh(j-1,i)-xdelh(j-1,i-1))/fac
          work(j,i,2) = (xdelh(j,i)+xdelh(j-1,i)- &
                         xdelh(j,i-1)-xdelh(j-1,i-1))/fac
        end do

        do concurrent ( j = jdi1:jdi2, i = idi1:idi2, nw = 1:2 )
          work(j,i,nw) = work(j,i,nw)*sfs%psdota(j,i)
        end do
        !
        ! compute divergence z from u and v
        ! ( u must be pstar * u ; similarly for v )
        ! ( note: map scale factors have been inverted in model (init) )
        !
        do concurrent ( j = jdi1:jdi2, i = idi1:idi2 )
          uu(j,i) = work(j,i,1)*mddom%msfd(j,i)
          vv(j,i) = work(j,i,2)*mddom%msfd(j,i)
        end do

        call exchange_rt(uu,1,jdi1,jdi2,idi1,idi2)
        call exchange_rt(vv,1,jdi1,jdi2,idi1,idi2)

        do concurrent ( j = jci1:jci2, i = ici1:ici2 )
          work(j,i,3) = rdx2 * map(j,i) * &
                  (-uu(j,i+1)+uu(j+1,i+1)-uu(j,i)+uu(j+1,i) + &
                    vv(j,i+1)+vv(j+1,i+1)-vv(j,i)-vv(j+1,i))
        end do

        do concurrent ( j = jci1:jci2, i = ici1:ici2 )
          deld(j,i,ns,n2) = deld(j,i,ns,n0) - dtau2*work(j,i,3) + &
                            deld(j,i,ns,3)/aam(ns)
          delh(j,i,ns,n2) = delh(j,i,ns,n0) - dtau2*hbar(ns) * &
                            deld(j,i,ns,n1)/sfs%psa(j,i) +     &
                            delh(j,i,ns,3)/aam(ns)
        end do
        !
        ! not in Madala (1987)
        !
        if ( ma%has_bdyleft ) then
          do i = ici1, ici2
            delh(jce1,i,ns,n2) = d_two*delh(jce1,i,ns,n1)-delh(jce1,i,ns,n0)
          end do
        end if
        if ( ma%has_bdyright ) then
          do i = ici1, ici2
            delh(jce2,i,ns,n2) = d_two*delh(jce2,i,ns,n1)-delh(jce2,i,ns,n0)
          end do
        end if
        if ( ma%has_bdybottom ) then
          do j = jce1, jce2
            delh(j,ice1,ns,n2) = d_two*delh(j,ice1,ns,n1)-delh(j,ice1,ns,n0)
          end do
        end if
        if ( ma%has_bdytop ) then
          do j = jce1, jce2
            delh(j,ice2,ns,n2) = d_two*delh(j,ice2,ns,n1)-delh(j,ice2,ns,n0)
          end do
        end if

        do concurrent ( j = jce1:jce2, i = ice1:ice2 )
          ddsum(j,i,ns) = ddsum(j,i,ns) + deld(j,i,ns,n2)
          dhsum(j,i,ns) = dhsum(j,i,ns) + delh(j,i,ns,n2)
        end do

        n0 = n1
        n1 = n2
        n2 = n0
      end do

    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
  end subroutine spstep

end module mod_split
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

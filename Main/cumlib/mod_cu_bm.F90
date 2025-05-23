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

module mod_cu_bm

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_constants
  use mod_memutil
  use mod_service
  use mod_cu_common
  use mod_runparams, only : iqv, dt, dtsec, ichem
  use mod_runparams, only : rhmin, rhmax
  use mod_regcm_types

!*****************************************************************
!
! Betts Miller Cumulus Convection Scheme
!
! Convective adjustment for deep or shallow convection
! Modified by Jack Kain of Penn State to replace the look-up table
!  by calculations.
!
!*****************************************************************
!                                                                *
!  references:                                                   *
!                                                                *
!  Betts, A.K., 1986:  A new convective adjustment scheme.       *
!    Part I: Observational and theoretical basis.  Quart. J. R.  *
!    Met. Soc., 112, 677-691.                                    *
!                                                                *
!  Betts, A.K., and M.J. Miller, 1986:  A new convective         *
!    adjustment scheme.  Part II: Single column tests using      *
!    gate wave, bomex, atex and arctic air mass data sets.       *
!    Quart. J. R. Met. Soc., 112, 693-709.                       *
!                                                                *
!  N.B.  Part of the code is scalar.  In global models           *
!  convection occurs in less than 30/100 points.  With           *
!  simulataneous vector processing for both deep and shallow     *
!  convection, there would be a lot of redundant vector          *
!  computations.  If vector processing is 10 times faster        *
!  than scalar, one might hope that the cpu time will be about   *
!  the same for both scalar and vector code.                     *
!                                                                *
!*****************************************************************
!
  implicit none

  private

  real(rkx) :: mxarg
  real(rkx), pointer, contiguous, dimension(:,:) :: cldefi
  real(rkx), pointer, contiguous, dimension(:,:,:) :: ape, q, qqmod, t
  real(rkx), pointer, contiguous, dimension(:,:,:) :: tmod, tref, z0
  real(rkx), pointer, contiguous, dimension(:) :: apek, apesk, difq, dift, &
         fpk, pdp, pk, psk, qk, qrefk, qsatk, therk, thsk,  &
         trefk, thvref, tk, tds
  real(rkx), pointer, contiguous, dimension(:,:) :: cldhgt, dsp0, dspb, dspt, &
         psp, xsm, thbt, thesp, ths, tthbt, tthes
  integer(ik4), pointer, contiguous, dimension(:,:) :: ifbuoy, ml
  integer(ik4), pointer, contiguous, dimension(:) :: kdp, ideep, ishal, jdeep, jshal

  public :: bmpara, allocate_mod_cu_bm
  public :: cldefi

  contains

  subroutine allocate_mod_cu_bm
    implicit none
    integer(ik4) :: intall
    call getmem2d(cldefi,jci1,jci2,ici1,ici2,'cu_bm:cldefi')
    call getmem1d(apek,1,kz,'cu_bm:apek')
    call getmem1d(tds,1,kz,'cu_bm:tds')
    call getmem1d(apesk,1,kz,'cu_bm:apesk')
    call getmem1d(difq,1,kz,'cu_bm:difq')
    call getmem1d(dift,1,kz,'cu_bm:dift')
    call getmem1d(fpk,1,kz,'cu_bm:fpk')
    call getmem1d(pdp,1,kz,'cu_bm:pdp')
    call getmem1d(pk,1,kz,'cu_bm:pk')
    call getmem1d(psk,1,kz,'cu_bm:psk')
    call getmem1d(qk,1,kz,'cu_bm:qk')
    call getmem1d(qrefk,1,kz,'cu_bm:qrefk')
    call getmem1d(qsatk,1,kz,'cu_bm:qsatk')
    call getmem1d(therk,1,kz,'cu_bm:therk')
    call getmem1d(thsk,1,kz,'cu_bm:thsk')
    call getmem1d(thvref,1,kz,'cu_bm:thrvref')
    call getmem1d(tk,1,kz,'cu_bm:tk')
    call getmem1d(trefk,1,kz,'cu_bm:trefk')
    call getmem1d(kdp,1,kz,'cu_bm:kdp')
    call getmem3d(ape,jci1,jci2,ici1,ici2,1,kz,'cu_bm:ape')
    call getmem3d(q,jci1,jci2,ici1,ici2,1,kz,'cu_bm:q')
    call getmem3d(qqmod,jci1,jci2,ici1,ici2,1,kz,'cu_bm:qqmod')
    call getmem3d(t,jci1,jci2,ici1,ici2,1,kz,'cu_bm:t')
    call getmem3d(tmod,jci1,jci2,ici1,ici2,1,kz,'cu_bm:tmod')
    call getmem3d(tref,jci1,jci2,ici1,ici2,1,kz,'cu_bm:tref')
    call getmem3d(z0,jci1,jci2,ici1,ici2,1,kz,'cu_bm:z0')
    call getmem2d(ifbuoy,jci1,jci2,ici1,ici2,'cu_bm:ifbuoy')
    call getmem2d(ml,jci1,jci2,ici1,ici2,'cu_bm:ml')
    call getmem2d(cldhgt,jci1,jci2,ici1,ici2,'cu_bm:cldhgt')
    call getmem2d(dsp0,jci1,jci2,ici1,ici2,'cu_bm:dsp0')
    call getmem2d(dspb,jci1,jci2,ici1,ici2,'cu_bm:dspb')
    call getmem2d(dspt,jci1,jci2,ici1,ici2,'cu_bm:dspt')
    call getmem2d(psp,jci1,jci2,ici1,ici2,'cu_bm:psp')
    call getmem2d(xsm,jci1,jci2,ici1,ici2,'cu_bm:xsm')
    call getmem2d(thbt,jci1,jci2,ici1,ici2,'cu_bm:thbt')
    call getmem2d(thesp,jci1,jci2,ici1,ici2,'cu_bm:thesp')
    call getmem2d(ths,jci1,jci2,ici1,ici2,'cu_bm:ths')
    call getmem2d(tthbt,jci1,jci2,ici1,ici2,'cu_bm:tthbt')
    call getmem2d(tthes,jci1,jci2,ici1,ici2,'cu_bm:tthes')
    intall = (jci2-jci1+1)*(ici2-ici1+1)
    call getmem1d(ideep,1,intall,'cu_bm:ideep')
    call getmem1d(ishal,1,intall,'cu_bm:ishal')
    call getmem1d(jdeep,1,intall,'cu_bm:jdeep')
    call getmem1d(jshal,1,intall,'cu_bm:jshal')
    mxarg = -log(epsilon(d_one))
  end subroutine allocate_mod_cu_bm

  subroutine bmpara(m2c)
    implicit none
    type(mod_2_cum), intent(in) :: m2c
    real(rkx), parameter :: h3000 = 3000.0_rkx
    real(rkx), parameter :: stresh = 1.10_rkx
    real(rkx), parameter :: stabs = 1.0_rkx
    real(rkx), parameter :: stabd = 0.90_rkx
    real(rkx), parameter :: rhf = 0.20_rkx
    real(rkx), parameter :: pmn = 6500.0_rkx
    real(rkx), parameter :: epsdn = 1.05_rkx
    real(rkx), parameter :: epsth = 6.0_rkx
    real(rkx), parameter :: pbm = 30000.0_rkx
    real(rkx), parameter :: pqm = 20000.0_rkx
    real(rkx), parameter :: pone = 2500.0_rkx
    real(rkx), parameter :: pfrz = 15000.0_rkx
    real(rkx), parameter :: pshu = 45000.0_rkx
    real(rkx), parameter :: zno = 750.0_rkx
    real(rkx), parameter :: zsh = 3999.0_rkx
    real(rkx), parameter :: fsl = 1.00_rkx
    real(rkx), parameter :: fss = 0.60_rkx
    real(rkx), parameter :: efimn = 0.20_rkx
    real(rkx), parameter :: efmnt = 0.70_rkx
    real(rkx), parameter :: fcc1 = 0.50_rkx
    real(rkx), parameter :: fcp = d_one - fcc1
    !real(rkx), parameter :: dspbfl = -4843.75_rkx
    !real(rkx), parameter :: dsp0fl = -7050.0_rkx
    !real(rkx), parameter :: dsptfl = -2250.0_rkx
    !real(rkx), parameter :: dspbfs = -3875.0_rkx
    !real(rkx), parameter :: dsp0fs = -5875.0_rkx
    !real(rkx), parameter :: dsptfs = -1875.0_rkx
    real(rkx), parameter :: dspbfl = -3875.0_rkx
    real(rkx), parameter :: dsp0fl = -5875.0_rkx
    real(rkx), parameter :: dsptfl = -1875.0_rkx
    real(rkx), parameter :: dspbfs = -3875.0_rkx
    real(rkx), parameter :: dsp0fs = -5875.0_rkx
    real(rkx), parameter :: dsptfs = -1875.0_rkx
    real(rkx), parameter :: dspbsl = dspbfl*fsl
    real(rkx), parameter :: dsp0sl = dsp0fl*fsl
    real(rkx), parameter :: dsptsl = dsptfl*fsl
    real(rkx), parameter :: dspbss = dspbfs*fss
    real(rkx), parameter :: dsp0ss = dsp0fs*fss
    real(rkx), parameter :: dsptss = dsptfs*fss
    ! Control if Deep Convection activated.
    real(rkx), parameter :: epsntp = 1.0e-3_rkx
    real(rkx), parameter :: efifc = 5.0_rkx
    real(rkx), parameter :: avgefi = (efimn+d_one)*d_half
    real(rkx), parameter :: dspc = -3000.0_rkx
    real(rkx), parameter :: epsp = 1.0e-7_rkx
    real(rkx), parameter :: stefi = avgefi
    real(rkx), parameter :: slopbl = (dspbfl-dspbsl)/(d_one-efimn)
    real(rkx), parameter :: slop0l = (dsp0fl-dsp0sl)/(d_one-efimn)
    real(rkx), parameter :: sloptl = (dsptfl-dsptsl)/(d_one-efimn)
    real(rkx), parameter :: slopbs = (dspbfs-dspbss)/(d_one-efimn)
    real(rkx), parameter :: slop0s = (dsp0fs-dsp0ss)/(d_one-efimn)
    real(rkx), parameter :: slopts = (dsptfs-dsptss)/(d_one-efimn)
    real(rkx), parameter :: slope = (d_one-efmnt)/(d_one-efimn)
    real(rkx), parameter :: a23m4l = c3les*(tzero-c4les)*wlhv
    real(rkx), parameter :: cprlg = cpd/(rhoh2o*egrav*wlhv)
    logical, save :: efinit = .false.

    real(rkx) :: avrgt, avrgtl, cthrs, den, dentpy, es, &
               dhdt, difql, diftl, dpkl, dpmix, dqref, drheat, &
               dsp, dsp0k, dspbk, dsptk, dst, dstq, dtdeta,    &
               ee, efi, fefi, fptk, hcorr, otsum, pdiff,       &
               pflag, pk0, pkb, pkt, potsum, pratec,      &
               preck, psum, pthrs, ptpk, qnew,  qotsum, &
               qrfkl, qrftp, qs, qsum, rdpsum, rhh, rhl,       &
               rotsum, rtbar, smix, sumde, sumdp, sumdt, tauk, &
               tcorr, tdpt, thskl, thtpk, thvmkl, tlcl, trfkl, &
               tskl, ztop, dthem, rdp0t, pbot, trel
    integer(ik4) :: i, j, iter, ivi, k, khdeep,    &
               khshal, kk, l, l0, l0m1, lb, lbm1, lbtk, lcor, &
               lqm, lshu, ltp1, ltpk, ltsh, n, &
               ndstn, ndstp
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'bmpara'
    integer(ik4), save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    trel = 3000.0_rkx
    tauk = dt/trel
    cthrs = (0.00635_rkx/secpd)*dt/cprlg
    !
    ! xsm is surface mask: =1 water; =0 land
    !
    do i = ici1, ici2
      do j = jci1, jci2
        if ( m2c%ldmsk(j,i) == 1 ) then
          xsm(j,i) = d_zero
        else
          xsm(j,i) = d_one
        end if
      end do
    end do
    if ( .not. efinit ) then
      do i = ici1, ici2
        do j = jci1, jci2
          cldefi(j,i) = avgefi*xsm(j,i) + stefi*(d_one-xsm(j,i))
        end do
      end do
      efinit = .true.
    end if
    !
    ! lb is currently set to kz-1
    !
    lb = kzm1
    !
    ! find melting level...
    !
    do i = ici1, ici2
      do j = jci1, jci2
        ml(j,i) = kzp1
      end do
    end do
    do k = 1, kz
      do i = ici1, ici2
        do j = jci1, jci2
          t(j,i,k) = m2c%tas(j,i,k)
          q(j,i,k) = m2c%qxas(j,i,k,iqv)
          ape(j,i,k) = d_one/(m2c%pas(j,i,k)/p00)**rovcp
          z0(j,i,k) = m2c%zas(j,i,k)
        end do
      end do
    end do
    do i = ici1, ici2
      do j = jci1, jci2
        do k = 1, kz
          if ( t(j,i,k) > tzero .and. ml(j,i) == kzp1 ) ml(j,i) = k
        end do
      end do
    end do
    do i = ici1, ici2
      do j = jci1, jci2
        cu_kbot(j,i) = kz
        thesp(j,i) = d_zero
        thbt(j,i) = d_zero
        psp(j,i) = 9.5e4_rkx
        tref(j,i,1) = t(j,i,1)
        ! fbuoy = 0 means no positive buoyancy; ifbuoy(j,i) means yes...
        ! p300 is the highest model level in the lowest 300 mb...
        ifbuoy(j,i) = 0
      end do
    end do
    !
    ! search for maximum buoyancy level
    !
    do k = 1, kz
      do i = ici1, ici2
        do j = jci1, jci2
          if ( m2c%pas(j,i,k) >= m2c%psf(j,i)-pbm ) then
            tthbt(j,i) = t(j,i,k)*ape(j,i,k)
            ee = m2c%pas(j,i,k)*q(j,i,k)/(0.622_rkx+q(j,i,k))
            tdpt = d_one/(rtzero-rwat*rwlhv*log(ee/611.0_rkx))
            tdpt = min(tdpt,t(j,i,k))
            tlcl = tdpt - (0.212_rkx+1.571e-3_rkx*(tdpt-tzero) - &
                           4.36e-4_rkx*(t(j,i,k)-tzero))*(t(j,i,k)-tdpt)
            tthes(j,i) = tthbt(j,i)*exp(wlhvocp*q(j,i,k)/tlcl)
            ! check for maximum buoyancy
            if ( tthes(j,i) > thesp(j,i) ) then
              psp(j,i) = p00*(tlcl/tthbt(j,i))**cpovr
              thbt(j,i) = tthbt(j,i)
              thesp(j,i) = tthes(j,i)
            end if
          end if
        end do
      end do
    end do
    !
    ! choose cloud base as model level just below psp
    !
    do k = 1, kzm1
      do i = ici1, ici2
        do j = jci1, jci2
          ! cloud bottom cannot be above pqm
          if ( m2c%pas(j,i,k) < psp(j,i) .and. &
               m2c%pas(j,i,k) >= pqm ) cu_kbot(j,i) = k + 1
        end do
      end do
    end do
    !
    ! warning: cu_kbot must not be gt kz-1 in shallow convection
    ! make sure the cloud base is at least 25 mb above the surface
    !
    do i = ici1, ici2
      do j = jci1, jci2
        pbot = m2c%pas(j,i,cu_kbot(j,i))
        if ( pbot >= m2c%psf(j,i)-pone .or. cu_kbot(j,i) >= kz ) then
          ! cloud bottom is at the surface so recalculate cloud bottom
          do k = 1, kzm1
            if ( m2c%pas(j,i,k) < m2c%psf(j,i)-pone ) cu_kbot(j,i) = k
          end do
        end if
      end do
    end do
    !
    ! cloud top computation
    !
    do i = ici1, ici2
      do j = jci1, jci2
        cu_ktop(j,i) = cu_kbot(j,i)
      end do
    end do
    do ivi = 1, kz
      l = kzp1 - ivi
      ! find environmental saturation equiv pot temp...
      do i = ici1, ici2
        do j = jci1, jci2
          es = aliq*exp((bliq*t(j,i,l)-cliq)/(t(j,i,l)-dliq))
          qs = 0.622_rkx * es/(m2c%pas(j,i,l)-es)
          ths(j,i) = t(j,i,l)*ape(j,i,l)*exp(wlhvocp*qs/t(j,i,l))
        end do
      end do
      ! buoyancy check
      do i = ici1, ici2
        do j = jci1, jci2
          if ( l <= cu_kbot(j,i) ) then
            if ( thesp(j,i) > ths(j,i) ) ifbuoy(j,i) = 1
            if ( thesp(j,i) > ths(j,i)-1.5_rkx .and. &
                 ifbuoy(j,i) == 1 ) cu_ktop(j,i) = l + 1
          end if
        end do
      end do
    end do
    !
    ! define and smooth dsps and cldefi
    !
    do i = ici1, ici2
      do j = jci1, jci2
        efi = cldefi(j,i)
        dspb(j,i) = ((efi-efimn)*slopbs+dspbss)*xsm(j,i) + &
                    ((efi-efimn)*slopbl+dspbsl)*(d_one-xsm(j,i))
        dsp0(j,i) = ((efi-efimn)*slop0s+dsp0ss)*xsm(j,i) + &
                    ((efi-efimn)*slop0l+dsp0sl)*(d_one-xsm(j,i))
        dspt(j,i) = ((efi-efimn)*slopts+dsptss)*xsm(j,i) + &
                    ((efi-efimn)*sloptl+dsptsl)*(d_one-xsm(j,i))
      end do
    end do
    !
    ! initialize changes of t and q due to convection
    !
    do k = 1, kz
      do i = ici1, ici2
        do j = jci1, jci2
          tmod(j,i,k) = d_zero
          qqmod(j,i,k) = d_zero
        end do
      end do
    end do
    !
    ! clean up and gather deep convection points
    !
    khdeep = 0
    do i = ici1, ici2
      do j = jci1, jci2
        if ( cu_ktop(j,i) > cu_kbot(j,i) ) then
          cu_ktop(j,i) = cu_kbot(j,i)
        end if
        cldhgt(j,i) = z0(j,i,cu_ktop(j,i)) - z0(j,i,cu_kbot(j,i))
        ! cloud is less than 90 mb deep or less than 3 sigma layers deep
        if ( cldhgt(j,i) < zno ) then
          cldefi(j,i) = avgefi*xsm(j,i) + stefi*(d_one-xsm(j,i))
        end if
        ! cloud has to be at least 2000m deep
        if ( cldhgt(j,i) >= zsh ) then
          khdeep = khdeep + 1
          ideep(khdeep) = i
          jdeep(khdeep) = j
        end if
      end do
    end do
    !
    ! horizontal loop for deep convection
    !
    do n = 1, khdeep
      i = ideep(n)
      j = jdeep(n)
      pratec = d_zero
      dentpy = d_zero
      avrgt = d_zero
      preck = d_zero
      ltpk = cu_ktop(j,i)
      lbtk = cu_kbot(j,i)
      !
      !dcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd
      !dcdcdcdcdcdc  deep convection   dcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd
      !dcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd
      !
      efi = cldefi(j,i)
      dspbk = dspb(j,i)
      dsp0k = dsp0(j,i)
      dsptk = dspt(j,i)
      !
      ! initialize variables in the convective column
      !
      do k = 1, kz
        dift(k) = d_zero
        difq(k) = d_zero
        tk(k) = t(j,i,k)
        trefk(k) = t(j,i,k)
        qk(k) = q(j,i,k)
        qrefk(k) = q(j,i,k)
        tref(j,i,k) = tpfc(m2c%pas(j,i,k),thesp(j,i),t(j,i,k),ape(j,i,k))
        pk(k) = m2c%pas(j,i,k)
        psk(k) = m2c%pas(j,i,k)
        apek(k) = ape(j,i,k)
        tds(k) = (m2c%pasf(j,i,k+1)-m2c%pasf(j,i,k))/m2c%pasf(j,i,kzp1)
        therk(k) = tref(j,i,k)*ape(j,i,k)
      end do
      !
      ! deep convection reference temperature profile
      !
      ltp1 = min(ltpk + 1,kz)
      lbm1 = max(lb - 1,1)
      pkb = pk(lb)
      pkt = pk(ltpk)
      !
      ! temperature reference profile below freezing level
      !
      l0 = lb
      l0m1 = l0-1
      pk0 = pk(lb)
      rdp0t = d_zero
      dthem = d_zero
      do l = ltpk, lbm1
        ivi = ltpk + lbm1 - l
        if ( trefk(ivi+1) <= 274.16_rkx ) then
          rdp0t = d_one/(pk0-pkt)
          dthem = therk(l0) - trefk(l0)*apek(l0)
          exit
        else
          trefk(ivi) = ((therk(ivi)-therk(ivi+1))*stabd + &
                         trefk(ivi+1)*apek(ivi+1))/apek(ivi)
        end if
        l0 = ivi
        pk0 = pk(l0)
      end do
      if ( l0 /= ltpk ) then
        l0m1 = l0-1
        do l = ltpk, l0m1
          trefk(l) = (therk(l)-(pk(l)-pkt)*dthem*rdp0t)/apek(l)
        end do
      else
        l0m1 = l0-1
      end if
      !
      ! deep convection reference humidity profile
      !
      do l = ltpk, lb
        !
        ! saturation pressure difference
        !
        if ( pkb-pk0 > pfrz ) then
          if ( l < l0 ) then
            dsp = ((pk0-pk(l))*dsptk+(pk(l)-pkt)*dsp0k)/(pk0-pkt)
          else
            dsp = ((pkb-pk(l))*dsp0k+(pk(l)-pk0)*dspbk)/(pkb-pk0)
          end if
        else
          dsp = dspc
        end if
        !
        ! humidity profile
        !
        if ( pk(l) > pqm ) then
          ! pressure must be below 200 mb
          psk(l) = pk(l) + dsp
          apesk(l) = d_one/(psk(l)/p00)**rovcp
          thsk(l) = trefk(l)*apek(l)
          qrefk(l) = pq0/psk(l)*exp(c3les*(thsk(l)-tzero*apesk(l)) / &
                     (thsk(l)-c4les*apesk(l)))
        else
          qrefk(l) = q(j,i,l)
        end if
      end do
      !
      ! enthalpy conservation integral
      !
      do iter = 1, 2
        sumde = d_zero
        sumdp = d_zero
        do l = ltpk, lb
          sumde = ((tk(l)-trefk(l))*cpd + &
                   (qk(l)-qrefk(l))*wlhv) * tds(l) + sumde
          sumdp = sumdp + tds(l)
        end do
        hcorr = sumde/(sumdp-tds(ltpk))
        lcor = ltpk + 1
        !
        ! find lqm
        !
        lqm = 0
        do l = 1, lb
          if ( pk(l) <= pqm ) lqm = l
        end do
        !
        ! above lqm correct temperature only
        !
        if ( lcor <= lqm ) then
          do l = lcor, lqm
            trefk(l) = trefk(l) + hcorr*rcpd
          end do
          lcor = lqm + 1
        end if
        !
        ! below lqm correct both temperature and moisture
        !
        do l = lcor, lb
          tskl = trefk(l)*apek(l)/apesk(l)
          dhdt = qrefk(l)*a23m4l/(tskl-c4les)**2 + cpd
          trefk(l) = hcorr/dhdt + trefk(l)
          thskl = trefk(l)*apek(l)
          qrefk(l) = pq0/psk(l) * exp(c3les*(thskl-tzero*apesk(l)) / &
                     (thskl-c4les*apesk(l)))
        end do
      end do
      do l = 1, kz
        thvref(l) = trefk(l)*apek(l)*(qrefk(l)*ep1+d_one)
      end do
      !
      ! heating, moistening, precipitation
      !
      do l = ltpk, lb
        diftl = (trefk(l)-tk(l))*tauk
        difql = (qrefk(l)-qk(l))*tauk
        avrgtl = (tk(l)+tk(l)+diftl)
        dentpy = (diftl*cpd+difql*wlhv)*tds(l)/avrgtl + dentpy
        avrgt = avrgtl*tds(l) + avrgt
        preck = tds(l)*diftl + preck
        dift(l) = diftl
        difq(l) = difql
      end do
      dentpy = dentpy + dentpy
      avrgt = avrgt/(sumdp+sumdp)
      if ( dentpy < epsntp .or. preck <= d_zero ) then
        cldefi(j,i) = efimn*xsm(j,i) + stefi*(d_one-xsm(j,i))
        ztop = z0(j,i,cu_kbot(j,i)) + zsh - 0.000001_rkx
        do l = 1, lb
          if ( z0(j,i,l) >= ztop ) cu_ktop(j,i) = l + 1
        end do
        !
        ! cloud must be at least 2 layers thick
        !
        if ( cu_kbot(j,i)-cu_ktop(j,i) < 2 ) cu_ktop(j,i) = cu_kbot(j,i) - 2
        cldhgt(j,i) = z0(j,i,cu_ktop(j,i)) - z0(j,i,cu_kbot(j,i))
        cycle
      end if
      !
      ! deep convection otherwise
      !
      total_precip_points = total_precip_points + 1
      ! keep the land value of efi equal to 1 until precip surpasses
      ! a threshold value, currently set to 0.25 inches per 24 hrs
      pthrs = cthrs/m2c%psf(j,i)
      drheat = (preck*xsm(j,i) + &
            max(epsp,preck-pthrs)*(d_one-xsm(j,i)))*cpd/avrgt
      efi = efifc*dentpy/drheat
      !
      efi = (cldefi(j,i)*fcp + efi*fcc1) * xsm(j,i) + (d_one-xsm(j,i))
      if ( efi > d_one ) efi = d_one
      if ( efi < efimn ) efi = efimn
      cldefi(j,i) = efi
      fefi = efmnt + slope*(efi-efimn)
      preck = preck*fefi
      !
      ! update precipitation, temperature & moisture
      !
      pratec = ((m2c%psf(j,i)*preck*cprlg)*d_100)/dtsec
      if ( pratec > dlowval ) then
        ! precipitation rate for surface (mm/s)
        cu_prate(j,i) = cu_prate(j,i) + pratec
      end if
      if ( ichem == 1 ) then
        do k = ltpk, kz
          cu_convpr(j,i,k) = pratec
        end do
      end if
      do l = ltpk, lb
        tmod(j,i,l) = dift(l)*fefi/dt
        qqmod(j,i,l) = difq(l)*fefi/dt
      end do
    end do
    !
    !dcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd
    !dcdcdcdcdcdc  end of deep convection  dcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd
    !dcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd
    !
    !
    ! gather shallow convection points
    !
    khshal = 0
    ndstn = 0
    ndstp = 0
    do i = ici1, ici2
      do j = jci1, jci2
        if ( cldhgt(j,i) >= zno .and. &
             cu_ktop(j,i) <= cu_kbot(j,i)-2 ) then
          if ( cldhgt(j,i) < zsh ) then
            khshal = khshal + 1
            ishal(khshal) = i
            jshal(khshal) = j
          end if
        end if
      end do
    end do
    !
    !scscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscs
    !scscscscscsc  shallow convection  cscscscscscscscscscscscscscscscscscsc
    !scscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscs
    !
    shallow: &
    do n = 1, khshal
      i = ishal(n)
      j = jshal(n)
      do k = 1, kz
        tk(k) = t(j,i,k)
        trefk(k) = t(j,i,k)
        qk(k) = q(j,i,k)
        qrefk(k) = q(j,i,k)
        qsatk(k) = q(j,i,k)
        pk(k) = m2c%pas(j,i,k)
        apek(k) = ape(j,i,k)
        thvmkl = t(j,i,k)*ape(j,i,k)*(q(j,i,k)*ep1+d_one)
        thvref(k) = thvmkl
        pdp(k) = pk(k) - pmn
        tds(k) = (m2c%pasf(j,i,k+1)-m2c%pasf(j,i,k))/m2c%pasf(j,i,kzp1)
      end do
      !
      ! find kdp...kdp(k) is the model level closest to 65 mb (pmn) above k;
      ! this is the depth over which relative humidity drop is measured to
      ! estimate shallow cloud top... see do 545...
      !
      do kk = kz, 1, -1
        pflag = abs(pk(kz)-pdp(kk))
        do k = kz - 1, 1, -1
          pdiff = abs(pk(k)-pdp(kk))
          if ( pdiff < pflag ) then
            pflag = pdiff
            if ( kk == k ) then
              kdp(kk) = k - 1
            else
              kdp(kk) = k
            end if
          end if
        end do
        kdp(kk) = max0(1,kdp(kk))
      end do
      !
      ! search for shallow cloud top
      !
      lbtk = cu_kbot(j,i)
      ltsh = lbtk
      lbm1 = lbtk - 1
      ztop = z0(j,i,cu_kbot(j,i)) + zsh - 0.000001_rkx
      !
      ! cloud top is level just above pbtk-psh
      !
      ltpk = 1
      do l = 1, kz
        if ( z0(j,i,l) >= ztop ) ltpk = l
      end do
      ptpk = pk(ltpk)
      !
      ! highest level allowed is level just below pshu
      !
      if ( ptpk <= pshu ) then
        lshu = 0
        do l = 1, kz
          if ( pk(l) <= pshu ) lshu = l + 1
        end do
        ltpk = lshu
        ptpk = pk(ltpk)
      end if
      ltp1 = ltpk + 1
      do l = ltpk, lbtk
        es = aliq*exp((bliq*tk(l)-cliq)/(tk(l)-dliq))
        qsatk(l) = 0.622_rkx * es/(pk(l)-es)
      end do
      do l = ltp1, lbm1
        rhl = min(max(qk(l)/qsatk(l),rhmin),rhmax)
        rhh = min(max(qk(kdp(l))/qsatk(kdp(l)),rhmin),rhmax)
        if ( rhh+rhf < rhl ) ltsh = l
      end do
      cu_ktop(j,i) = ltsh
      ltp1 = ltsh
      ltpk = ltsh - 1
      cldhgt(j,i) = z0(j,i,cu_ktop(j,i)) - z0(j,i,cu_kbot(j,i))
      ! if cloud is not at least 90 mb or 3 sigma layers deep, then no cloud
      if ( cldhgt(j,i) < zno .or. cu_ktop(j,i) > cu_kbot(j,i)-2 ) then
        cu_ktop(j,i) = cu_kbot(j,i)
        cycle
      end if
      ! scaling potential temperature & table index at top
      thtpk = t(j,i,ltp1)*ape(j,i,ltp1)
      ee = m2c%pas(j,i,ltp1)*q(j,i,ltp1)/(0.622_rkx+q(j,i,ltp1))
      tdpt = d_one/(rtzero-rwat*rwlhv*log(ee/611._rkx))
      tdpt = min(tdpt,t(j,i,ltp1))
      tlcl = tdpt - (0.212_rkx+1.571e-3_rkx*(tdpt-tzero) - &
             4.36e-4_rkx*(t(j,i,ltp1)-tzero))*(t(j,i,ltp1)-tdpt)
      ptpk = p00*(tlcl/thtpk)**cpovr
      dpmix = ptpk - psp(j,i)
      if ( abs(dpmix) < h3000 ) dpmix = -h3000
      !
      ! temperature propfile slope
      !
      smix = (thtpk-thbt(j,i))/dpmix*stabs
      do l = ltp1, lbtk
        ivi = min(ltp1+lbtk-l,kz-1)
        trefk(ivi) = ((pk(ivi)-pk(ivi+1))*smix + &
                     trefk(ivi+1)*apek(ivi+1))/apek(ivi)
      end do
      !
      ! temperature reference profile correction
      !
      sumdt = d_zero
      sumdp = d_zero
      do l = ltp1, lbtk
        sumdt = (tk(l)-trefk(l))*tds(l) + sumdt
        sumdp = sumdp + tds(l)
      end do
      rdpsum = d_one/sumdp
      fpk(lbtk) = trefk(lbtk)
      tcorr = sumdt*rdpsum
      do l = ltp1, lbtk
        trfkl = trefk(l) + tcorr
        trefk(l) = trfkl
        fpk(l) = trfkl
      end do
      !
      ! humidity profile equations
      !
      psum = d_zero
      qsum = d_zero
      potsum = d_zero
      qotsum = d_zero
      otsum = d_zero
      dst = d_zero
      fptk = fpk(ltp1)
      do l = ltp1, lbtk
        dpkl = fpk(l) - fptk
        psum = dpkl*tds(l) + psum
        qsum = qk(l)*tds(l) + qsum
        rtbar = 2.0_rkx/(trefk(l)+tk(l))
        otsum = tds(l)*rtbar + otsum
        potsum = dpkl*rtbar*tds(l) + potsum
        qotsum = qk(l)*rtbar*tds(l) + qotsum
        dst = (trefk(l)-tk(l))*rtbar*tds(l) + dst
      end do

      psum = psum*rdpsum
      qsum = qsum*rdpsum
      rotsum = d_one/otsum
      potsum = potsum*rotsum
      qotsum = qotsum*rotsum
      dst = dst*rotsum*cpowlhv
      !
      ! ensure positive entropy change
      !
      if ( dst > d_zero ) then
        cu_ktop(j,i) = cu_kbot(j,i)
        ndstp = ndstp + 1
        cycle
      else
        dstq = dst*epsdn
      end if
      !
      ! check for isothermal atmosphere
      !
      den = potsum - psum
      if ( -den/psum < 0.00005_rkx ) then
        cu_ktop(j,i) = cu_kbot(j,i)
        cycle
      else
        !
        ! slope of the reference humidity profile
        !
        dqref = (qotsum-dstq-qsum)/den
      end if
      !
      ! humidity doesn`t increase with height
      !
      if ( dqref < d_zero ) then
        cu_ktop(j,i) = cu_kbot(j,i)
        cycle
      end if
      !
      ! humidity at the cloud top
      !
      qrftp = qsum - dqref*psum
      !
      ! humidity profile
      !
      do l = ltp1, lbtk
        qrfkl = (fpk(l)-fptk)*dqref + qrftp
        ! supersaturation not allowed
        qnew = (qrfkl-qk(l))*tauk + qk(l)
        if ( qnew > qsatk(l)*stresh ) then
          cu_ktop(j,i) = cu_kbot(j,i)
          exit shallow
        end if
        thvref(l) = trefk(l)*apek(l)*(qrfkl*ep1+d_one)
        qrefk(l) = qrfkl
      end do
      !
      ! eliminate impossible slopes (betts, dtheta/dq)
      !
      do l = ltp1, lbtk
        dtdeta = (thvref(l-1)-thvref(l))/tds(l)
        if ( dtdeta < epsth ) then
          cu_ktop(j,i) = cu_kbot(j,i)
          exit shallow
        end if
      end do
      if ( dst > d_zero ) then
        ndstp = ndstp + 1
      else
        ndstn = ndstn + 1
      end if
      dentpy = d_zero
      do l = ltp1, lbtk
        dentpy = ((trefk(l)-tk(l))*cpd+(qrefk(l)-qk(l))*wlhv) / &
                  (tk(l)+trefk(l))*tds(l) + dentpy
      end do
      !
      ! relaxation towards reference profiles
      !
      do l = ltp1, lbtk
        tmod(j,i,l) = (trefk(l)-tk(l))/trel
        qqmod(j,i,l) = (qrefk(l)-qk(l))/trel
      end do
    end do shallow
    !
    !scscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscs
    !scscscscscsc  end of shallow convection   scscscscscscscscscscscscscscs
    !scscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscs
    !
    do k = 1, kz
      do i = ici1, ici2
        do j = jci1, jci2
          cu_tten(j,i,k)  = tmod(j,i,k)
          cu_qten(j,i,k,iqv) = qqmod(j,i,k)
        end do
      end do
    end do
#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif
    contains

      !
      ! Calculates tpfc
      !
      pure real(rkx) function tpfc(press,thetae,tgs,pi)
        implicit none
        real(rkx), intent(in) :: pi, press, tgs, thetae
        real(rkx) :: qs, dtx, f1, fo, rp, t1, tguess
        real(rkx) :: es
        real(rkx), parameter :: rl461 = wlhv/rwat
        real(rkx), parameter :: rl1004 = wlhv/cpd
        integer(ik4) :: iloop
        integer(ik4), parameter :: maxiter = 100
        !
        ! iteratively extract temperature from equivalent potential temperature.
        !
        rp = thetae/pi
        es = 611.0_rkx * exp(rl461*(rtzero-d_one/tgs))
        qs = 0.622_rkx*es/(press-es)
        fo = tgs * exp(rl1004*qs/tgs) - rp
        t1 = tgs - d_half*fo
        tguess = tgs
        iloop = 0
        do
          es = 611.0_rkx * exp(rl461*(rtzero-d_one/t1))
          qs = 0.622_rkx*es/(press-es)
          f1 = t1*exp(rl1004*qs/t1) - rp
          if ( abs(f1) < 0.1_rkx ) then
            tpfc = t1
            exit
          else
            dtx = f1 * (t1-tguess)/(f1-fo)
            tguess = t1
            fo = f1
            t1 = t1 - dtx
            iloop = iloop + 1
            if ( iloop > maxiter ) then
              tpfc = t1
              exit
            end if
          end if
        end do
      end function tpfc

  end subroutine bmpara

end module mod_cu_bm
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

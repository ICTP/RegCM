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
 
module mod_cu_bm
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
  use mod_dynparam
  use mod_memutil
  use mod_cu_common

  implicit none

  private

  integer , parameter :: itb = 100
  integer , parameter :: jtb = 150

  real(8) :: pl
  real(8) , pointer , dimension(:,:,:) :: tbase
  real(8) , pointer , dimension(:,:) :: cldefi
  real(8) , pointer , dimension(:,:,:) :: ape , q , qqmod , t , &
         tmod , tref , z0
  real(8) , pointer , dimension(:) :: apek , apesk , difq , dift , ddzq , &
         fpk , pdp , pk , psk , qk , qrefk , qsatk , therk , thsk ,       &
         thvref , tk , trefk
  real(8) , pointer , dimension(:,:) :: cldhgt , dsp0 , dspb , dspt , p , &
         pbot , prtop , psp , xsm , thbt , thesp , ths , tthbt , tthes
  integer , pointer , dimension(:,:) :: ifbuoy , ip300 , lbot , ltop , ml
  integer , pointer , dimension(:) :: kdp , nbotd , nbots , ndpthd ,      &
         ndpths , ntopd , ntops , ideep , ishal , jdeep , jshal

  public :: bmpara , lutbl , allocate_mod_cu_bm
  public :: tbase , cldefi

  contains
!
  subroutine allocate_mod_cu_bm
    implicit none
    call getmem3d(tbase,1,iy,1,kz,1,jxp,'cu_bm:tbase')
    call getmem2d(cldefi,1,iy,1,jxp,'cu_bm:cldefi')
    call getmem1d(apek,1,kz,'cu_bm:apek')
    call getmem1d(apesk,1,kz,'cu_bm:apesk')
    call getmem1d(difq,1,kz,'cu_bm:difq')
    call getmem1d(dift,1,kz,'cu_bm:dift')
    call getmem1d(ddzq,1,kz,'cu_bm:ddzq')
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
    call getmem1d(nbotd,1,kz,'cu_bm:nbotd')
    call getmem1d(nbots,1,kz,'cu_bm:nbots')
    call getmem1d(ndpthd,1,kz,'cu_bm:ndpthd')
    call getmem1d(ndpths,1,kz,'cu_bm:ndpths')
    call getmem1d(ntopd,1,kz,'cu_bm:ntopd')
    call getmem1d(ntops,1,kz,'cu_bm:ntops')
    call getmem3d(ape,1,jxp,1,iy,1,kz,'cu_bm:ape')
    call getmem3d(q,1,jxp,1,iy,1,kz,'cu_bm:q')
    call getmem3d(qqmod,1,jxp,1,iy,1,kz,'cu_bm:qqmod')
    call getmem3d(t,1,jxp,1,iy,1,kz,'cu_bm:t')
    call getmem3d(tmod,1,jxp,1,iy,1,kz,'cu_bm:tmod')
    call getmem3d(tref,1,jxp,1,iy,1,kz,'cu_bm:tref')
    call getmem3d(z0,1,jxp,1,iy,1,kz,'cu_bm:z0')
    call getmem2d(ifbuoy,1,jxp,1,iy,'cu_bm:ifbuoy')
    call getmem2d(ip300,1,jxp,1,iy,'cu_bm:ip300')
    call getmem1d(ideep,1,iy*jxp,'cu_bm:ideep')
    call getmem1d(ishal,1,iy*jxp,'cu_bm:ishal')
    call getmem1d(jdeep,1,iy*jxp,'cu_bm:jdeep')
    call getmem1d(jshal,1,iy*jxp,'cu_bm:jshal')
    call getmem2d(lbot,1,jxp,1,iy,'cu_bm:lbot')
    call getmem2d(ltop,1,jxp,1,iy,'cu_bm:ltop')
    call getmem2d(ml,1,jxp,1,iy,'cu_bm:ml')
    call getmem2d(cldhgt,1,jxp,1,iy,'cu_bm:cldhgt')
    call getmem2d(dsp0,1,jxp,1,iy,'cu_bm:dsp0')
    call getmem2d(dspb,1,jxp,1,iy,'cu_bm:dspb')
    call getmem2d(dspt,1,jxp,1,iy,'cu_bm:dspt')
    call getmem2d(p,1,jxp,1,iy,'cu_bm:p')
    call getmem2d(pbot,1,jxp,1,iy,'cu_bm:pbot')
    call getmem2d(prtop,1,jxp,1,iy,'cu_bm:prtop')
    call getmem2d(psp,1,jxp,1,iy,'cu_bm:psp')
    call getmem2d(xsm,1,jxp,1,iy,'cu_bm:xsm')
    call getmem2d(thbt,1,jxp,1,iy,'cu_bm:thbt')
    call getmem2d(thesp,1,jxp,1,iy,'cu_bm:thesp')
    call getmem2d(ths,1,jxp,1,iy,'cu_bm:ths')
    call getmem2d(tthbt,1,jxp,1,iy,'cu_bm:tthbt')
    call getmem2d(tthes,1,jxp,1,iy,'cu_bm:tthes')
  end subroutine allocate_mod_cu_bm

  subroutine bmpara(jstart,jend,istart,iend,ktau)
    implicit none
    integer , intent(in) :: jstart , jend , istart , iend
    integer(8) , intent(in) :: ktau
!
    real(8) , parameter :: h1 = 1.0D0
    real(8) , parameter :: h3000 = 3000.0D0
    real(8) , parameter :: h10e5 = 100000.0D0
    real(8) , parameter :: d608 = 0.608D0
    real(8) , parameter :: dm2859 = -rgas/cpd
    real(8) , parameter :: epsq = 2.0D-12
    real(8) , parameter :: row = d_1000
    real(8) , parameter :: t1 = tzero+1.0D0
    real(8) , parameter :: stresh = 1.10D0
    real(8) , parameter :: stabs = 1.0D0
    real(8) , parameter :: stabd = 0.90D0
    real(8) , parameter :: rhf = 0.20D0
    real(8) , parameter :: pmn = 6500.0D0
    real(8) , parameter :: epsdn = 1.05D0
    real(8) , parameter :: epsth = 6.0D00
    real(8) , parameter :: pbm = 30000.0D0
    real(8) , parameter :: pqm = 20000.0D0
    real(8) , parameter :: pone = 2500.0D0
    real(8) , parameter :: pfrz = 15000.0D0
    real(8) , parameter :: pshu = 45000.0D0
    real(8) , parameter :: zno = 750.0D0
    real(8) , parameter :: zsh = 3999.0D0
    real(8) , parameter :: fss = 0.60D0
    real(8) , parameter :: efimn = 0.20D0
    real(8) , parameter :: efmnt = 0.70D0
    real(8) , parameter :: fcc1 = 0.50D0
    real(8) , parameter :: fcp = h1 - fcc1
    real(8) , parameter :: dspbfl = -3875.0D0
    real(8) , parameter :: dsp0fl = -5875.0D0
    real(8) , parameter :: dsptfl = -1875.0D0
    real(8) , parameter :: fsl = 1.0D0
    real(8) , parameter :: dspbfs = -3875.0D0
    real(8) , parameter :: dsp0fs = -5875.0D0
    real(8) , parameter :: dsptfs = -1875.0D0
    real(8) , parameter :: dspbsl = dspbfl*fsl
    real(8) , parameter :: dsp0sl = dsp0fl*fsl
    real(8) , parameter :: dsptsl = dsptfl*fsl
    real(8) , parameter :: dspbss = dspbfs*fss
    real(8) , parameter :: dsp0ss = dsp0fs*fss
    real(8) , parameter :: dsptss = dsptfs*fss
    real(8) , parameter :: epsntp = 0.0010D0
    real(8) , parameter :: efifc = 5.0D0
    real(8) , parameter :: avgefi = (efimn+1.0D0)*d_half
    real(8) , parameter :: dspc = -3000.0D0
    real(8) , parameter :: epsp = 1.0D-7
    real(8) , parameter :: stefi = avgefi
    real(8) , parameter :: slopbl = (dspbfl-dspbsl)/(h1-efimn)
    real(8) , parameter :: slop0l = (dsp0fl-dsp0sl)/(h1-efimn)
    real(8) , parameter :: sloptl = (dsptfl-dsptsl)/(h1-efimn)
    real(8) , parameter :: slopbs = (dspbfs-dspbss)/(h1-efimn)
    real(8) , parameter :: slop0s = (dsp0fs-dsp0ss)/(h1-efimn)
    real(8) , parameter :: slopts = (dsptfs-dsptss)/(h1-efimn)
    real(8) , parameter :: slope = (h1-efmnt)/(h1-efimn)
    real(8) , parameter :: a23m4l = c3les*(tzero-c4les)*wlhv
    real(8) , parameter :: cporng = d_one/dm2859
    real(8) , parameter :: elocp = wlhv/cpd
    real(8) , parameter :: cprlg = cpd/(row*egrav*wlhv)
    logical , parameter :: unis = .false.
    logical , parameter :: unil = .true.
    logical , parameter :: oct90 = .true.
!
    real(8) :: ak , akclth , apekl , avrgt , avrgtl , cell , &
               cthrs , den , dentpy , dhdt , difql , diftl , dpkl ,   &
               dpmix , dqref , drheat , dsp , dsp0k , dspbk , dsptk , &
               dst , dstq , dtdeta , dthem , ee , efi , es , fefi ,   &
               fptk , hcorr , otsum , pdiff , pdiffk , pflag , pk0 ,  &
               pkb , pkl , pkt , potsum , pppk , prainx , preck ,     &
               psfck , psum , pthrs , ptpk , qkl , qnew , qotsum ,    &
               qrfkl , qrftp , qs , qsum , qu , rdp0t , rdpsum , rhh ,&
               rhl , rotsum , rtbar , smix , stabdl , sumde , sumdp , &
               sumdt , tauk , tcorr , tdpt , thskl , thtpk , thvmkl , &
               tkl , tlcl , trfkl , tskl , ztop
    integer :: i , j , iconss , iter , ivi , k , kb , kbaseb ,    &
               kclth , khdeep , khshal , kk , l , l0 , l0m1 , lb ,    &
               lbm1 , lbtk , lcor , lqm , lshu , ltp1 , ltpk , ltsh , &
               n , ndeep , ndepth , ndstn , ndstp , nshal , nswap , ll
!
!-----------------------------------------------------------------------
!
    lqm = 0
    lshu = 0
!
    do k = 1 , kz
      do i = istart , iend
        do j = jstart , jend
          rcldlwc(j,i,k) = d_zero
          rcldfra(j,i,k) = d_zero
        end do
      end do
    end do
    if ( lchem ) then
!
!     icumtop = top level of cumulus clouds
!     icumbot = bottom level of cumulus clouds
!     (calculated in cupara and stored for tractend)
      do i = istart , iend
        do j = jstart , jend
          icumtop(j,i) = 0
          icumbot(j,i) = 0
        end do
      end do
    end if
    total_precip_points = 0
    iconss = 0
    tauk = dtcum/trel
    cthrs = (0.006350D0/secpd)*dtcum/cprlg
!
!-----------------------------------------------------------------------
!
!   xsm is surface mask: =1 water; =0 land
!
    do i = istart , iend
      do j = jstart , jend
        if ( lmask(i,j) == 0 ) then
          xsm(j,i) = d_one
        else
          xsm(j,i) = d_zero
        end if
      end do
    end do
    if ( ktau == 0 ) then
      do i = istart , iend
        do j = jstart , jend
          cldefi(i,j) = avgefi*xsm(j,i) + stefi*(h1-xsm(j,i))
        end do
      end do
    end if
!
!   lb is currently set to kz-1
!
    lb = kz - 1
    do k = 1 , kz
      ntopd(k) = 0
      nbotd(k) = 0
      ntops(k) = 0
      nbots(k) = 0
      ndpths(k) = 0
      ndpthd(k) = 0
    end do
!
!   find melting level...
!
    do i = istart , iend
      do j = jstart , jend
        ml(j,i) = kzp1
      end do
    end do

    do i = istart , iend
      do j = jstart , jend
        do k = 1 , kz
          t(j,i,k) = tas(i,k,j)
          if ( t(j,i,k) > tzero .and. ml(j,i) == kzp1 ) ml(j,i) = k
          q(j,i,k) = qvas(i,k,j)
          pppk = (hlev(k)*sfcps(i,j)+ptop)*d_1000
          ape(j,i,k) = (pppk/h10e5)**dm2859
        end do
        lbot(j,i) = kz
        thesp(j,i) = d_zero
        thbt(j,i) = d_zero
        psp(j,i) = 9.5D4
        tref(j,i,1) = t(j,i,1)
        ! fbuoy = 0 means no positive buoyancy; ifbuoy(j,i) means yes...
        ! p300 is the highest model level in the lowest 300 mb...
        ifbuoy(j,i) = 0
        ip300(j,i) = 0
        cell = ptop/sfcps(i,j)
        do k = 1 , kz
          ddzq(k) = rovg*tbase(i,k,j)*dlog((flev(k+1)+cell)/(flev(k)+cell))
        end do
        z0(j,i,kz) = d_half*ddzq(kz)
        do k = kz - 1 , 1 , -1
          z0(j,i,k) = z0(j,i,k+1) + d_half*(ddzq(k)+ddzq(k+1))
        end do
      end do
    end do
!
!   padding specific humidity if too small
!
    do k = 1 , kz
      do i = istart , iend
        do j = jstart , jend
          if ( q(j,i,k) < epsq ) q(j,i,k) = epsq
          pdiff = (d_one-hlev(k))*sfcps(i,j)
          if ( pdiff < 30.0D0 .and. ip300(j,i) == 0 ) ip300(j,i) = k
        end do
      end do
    end do
!
!   search for maximum buoyancy level
!
    do kb = 1 , kz
      do i = istart , iend
        do j = jstart , jend
          pkl = (hlev(kb)*sfcps(i,j)+ptop)*d_1000
          psfck = (hlev(kz)*sfcps(i,j)+ptop)*d_1000
          if ( pkl >= psfck-pbm ) then
            tthbt(j,i) = t(j,i,kb)*ape(j,i,kb)
            ee = pkl*q(j,i,kb)/(ep2+q(j,i,kb))
            tdpt = d_one/(rtzero-rwat/wlhv*dlog(ee/611.D0))
            tdpt = dmin1(tdpt,t(j,i,kb))
            tlcl = tdpt - (0.212D0+1.571D-3*(tdpt-tzero) - &
                           4.36D-4*(t(j,i,kb)-tzero))*(t(j,i,kb)-tdpt)
            tthes(j,i) = tthbt(j,i)*dexp(elocp*q(j,i,kb)/tlcl)
!           check for maximum buoyancy
            if ( tthes(j,i) > thesp(j,i) ) then
              psp(j,i) = h10e5*(tthbt(j,i)/tlcl)**cporng
              thbt(j,i) = tthbt(j,i)
              thesp(j,i) = tthes(j,i)
            end if
          end if
        end do
      end do
    end do
!
!   choose cloud base as model level just below psp
!
    do k = 1 , kzm1
      ak = hlev(k)
      do i = istart , iend
        do j = jstart , jend
          p(j,i) = (ak*sfcps(i,j)+ptop)*d_1000
!         cloud bottom cannot be above 200 mb
          if ( p(j,i) < psp(j,i) .and. p(j,i) >= pqm ) lbot(j,i) = k + 1
        end do
      end do
    end do
!    
!   warning: lbot must not be gt kz-1 in shallow convection
!   make sure the cloud base is at least 25 mb above the surface
!
    do i = istart , iend
      do j = jstart , jend
        pbot(j,i) = (hlev(lbot(j,i))*sfcps(i,j)+ptop)*d_1000
        psfck = (hlev(kz)*sfcps(i,j)+ptop)*d_1000
        if ( pbot(j,i) >= psfck-pone .or. lbot(j,i) >= kz ) then
!         cloud bottom is at the surface so recalculate cloud bottom
          do k = 1 , kzm1
            p(j,i) = (hlev(kz)*sfcps(i,j)+ptop)*d_1000
            if ( p(j,i) < psfck-pone ) lbot(j,i) = k
          end do
          pbot(j,i) = (hlev(lbot(j,i))*sfcps(i,j)+ptop)*d_1000
        end if
      end do
    end do
!    
!   cloud top computation
!    
    do i = istart , iend
      do j = jstart , jend
        prtop(j,i) = pbot(j,i)
        ltop(j,i) = lbot(j,i)
      end do
    end do
    do ivi = 1 , kz
      l = kzp1 - ivi
!     find environmental saturation equiv pot temp...
      do i = istart , iend
        do j = jstart , jend
          p(j,i) = (hlev(l)*sfcps(i,j)+ptop)*d_1000
          es = aliq*dexp((bliq*t(j,i,l)-cliq)/(t(j,i,l)-dliq))
          qs = ep2*es/(p(j,i)-es)
          ths(j,i) = t(j,i,l)*ape(j,i,l)*dexp(elocp*qs/t(j,i,l))
        end do
      end do
!     buoyancy check
      do i = istart , iend
        do j = jstart , jend
          if ( l <= lbot(j,i) ) then
            if ( thesp(j,i) > ths(j,i) ) ifbuoy(j,i) = 1
            if ( thesp(j,i) > ths(j,i)-1.5D0 .and. &
                 ifbuoy(j,i) == 1 ) ltop(j,i) = l + 1
          end if
        end do
      end do
    end do
!    
!   cloud top pressure
!    
    do i = istart , iend
      do j = jstart , jend
        prtop(j,i) = (hlev(ltop(j,i))*sfcps(i,j)+ptop)*d_1000
      end do
    end do
!
!-----------------------------------------------------------------------
!
!   define and smooth dsps and cldefi
!
    if ( unis ) then
      do i = istart , iend
        do j = jstart , jend
          efi = cldefi(i,j)
          dspb(j,i) = (efi-efimn)*slopbs + dspbss
          dsp0(j,i) = (efi-efimn)*slop0s + dsp0ss
          dspt(j,i) = (efi-efimn)*slopts + dsptss
        end do
      end do
    else if ( .not.unil ) then
      do i = istart , iend
        do j = jstart , jend
          efi = cldefi(i,j)
          dspb(j,i) = ((efi-efimn)*slopbs+dspbss)*xsm(j,i) +    &
                    ((efi-efimn)*slopbl+dspbsl)*(h1-xsm(j,i))
          dsp0(j,i) = ((efi-efimn)*slop0s+dsp0ss)*xsm(j,i) +    &
                    ((efi-efimn)*slop0l+dsp0sl)*(h1-xsm(j,i))
          dspt(j,i) = ((efi-efimn)*slopts+dsptss)*xsm(j,i) +    &
                    ((efi-efimn)*sloptl+dsptsl)*(h1-xsm(j,i))
        end do
      end do
    else
      do i = istart , iend
        do j = jstart , jend
          efi = cldefi(i,j)
          dspb(j,i) = ((efi-efimn)*slopbl+dspbsl)
          dsp0(j,i) = ((efi-efimn)*slop0l+dsp0sl)
          dspt(j,i) = ((efi-efimn)*sloptl+dsptsl)
        end do
      end do
    end if
!
!   initialize changes of t and q due to convection
!
    do k = 1 , kz
      do i = istart , iend
        do j = jstart , jend
          tmod(j,i,k) = d_zero
          qqmod(j,i,k) = d_zero
        end do
      end do
    end do
!
!   clean up and gather deep convection points
!
    khdeep = 0
    nswap = 0
    do i = istart , iend
      do j = jstart , jend
        if ( ltop(j,i) > lbot(j,i) ) then
          ltop(j,i) = lbot(j,i)
          prtop(j,i) = pbot(j,i)
        end if
        cldhgt(j,i) = z0(j,i,ltop(j,i)) - z0(j,i,lbot(j,i))
!       cloud is less than 90 mb deep or less than 3 sigma layers deep
        if ( cldhgt(j,i) < zno ) cldefi(i,j) = avgefi*xsm(j,i) &
             + stefi*(h1-xsm(j,i))
!       cloud has to be at least 290 mb deep
        if ( cldhgt(j,i) >= zsh ) then
          khdeep = khdeep + 1
          ideep(khdeep) = i
          jdeep(khdeep) = j
        end if
      end do
    end do
!
!   horizontal loop for deep convection
!
    do n = 1 , khdeep
      i = ideep(n)
      j = jdeep(n)
      dentpy = d_zero
      avrgt = d_zero
      preck = d_zero
      ltpk = ltop(j,i)
      lbtk = lbot(j,i)
!
!dcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd
!dcdcdcdcdcdc  deep convection   dcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd
!dcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd
!
      efi = cldefi(i,j)
      dspbk = dspb(j,i)
      dsp0k = dsp0(j,i)
      dsptk = dspt(j,i)
!
!     initialize variables in the convective column
!
      do k = 1 , kz
        dift(k) = d_zero
        difq(k) = d_zero
        tkl = t(j,i,k)
        tk(k) = tkl
        trefk(k) = tkl
        qkl = q(j,i,k)
        qk(k) = qkl
        qrefk(k) = qkl
        pkl = (hlev(k)*sfcps(i,j)+ptop)*d_1000
        tref(j,i,k) = tpfc(pkl,thesp(j,i),t(j,i,k),wlhv,qu,ape(j,i,k))
        pk(k) = pkl
        psk(k) = pkl
        apekl = ape(j,i,k)
        apek(k) = apekl
        therk(k) = tref(j,i,k)*apekl
      end do
!
!     deep convection reference temperature profile
!
      ltp1 = ltpk + 1
      lbm1 = lb - 1
      pkb = pk(lb)
      pkt = pk(ltpk)
!
!     temperature reference profile below freezing level
!
      l0 = lb
      pk0 = pk(lb)
      do l = ltpk , lbm1
        ivi = ltpk + lbm1 - l
        if ( trefk(ivi+1) <= t1 ) then
!
!         temperature reference profile above freezing level
!
          l0m1 = l0 - 1
          rdp0t = h1/(pk0-pkt)
          dthem = therk(l0) - trefk(l0)*apek(l0)
          do ll = ltpk , l0m1
            trefk(l) = (therk(l)-(pk(l)-pkt)*dthem*rdp0t)/apek(l)
          end do
          go to 50
        else
          stabdl = stabd
          trefk(ivi) = ((therk(ivi)-therk(ivi+1))*stabdl + &
                         trefk(ivi+1)*apek(ivi+1))/apek(ivi)
          l0 = ivi
          pk0 = pk(l0)
        end if
      end do
!
!     freezing level at or above the cloud top
!
      l0m1 = l0 - 1
!
!     deep convection reference humidity profile
!
   50     continue
!
      do l = ltpk , lb
!
!       saturation pressure difference
!
        if ( pkb-pk0 < pfrz ) then
          dsp = dspc
        else if ( l < l0 ) then
          dsp = ((pk0-pk(l))*dsptk+(pk(l)-pkt)*dsp0k)/(pk0-pkt)
        else
          dsp = ((pkb-pk(l))*dsp0k+(pk(l)-pk0)*dspbk)/(pkb-pk0)
        end if
!
!       humidity profile
!
        if ( pk(l) > pqm ) then
!         pressure must be below 200 mb
          psk(l) = pk(l) + dsp
          apesk(l) = (psk(l)/h10e5)**dm2859
          thsk(l) = trefk(l)*apek(l)
          qrefk(l) = pq0/psk(l)*dexp(c3les*(thsk(l)-tzero*apesk(l)) / &
                     (thsk(l)-c4les*apesk(l)))
        else
          qrefk(l) = q(j,i,l)
        end if
      end do
!
!     enthalpy conservation integral
!
      do iter = 1 , 2
        sumde = d_zero
        sumdp = d_zero
        do l = ltpk , lb
          sumde = ((tk(l)-trefk(l))*cpd+(qk(l)-qrefk(l))*wlhv)*dflev(l) + sumde
          sumdp = sumdp + dflev(l)
        end do
        hcorr = sumde/(sumdp-dflev(ltpk))
        lcor = ltpk + 1
!
!       find lqm
!
        do l = 1 , lb
          if ( pk(l) <= pqm ) lqm = l
        end do
!
!       above lqm correct temperature only
!
        if ( lcor <= lqm ) then
          do l = lcor , lqm
            trefk(l) = trefk(l) + hcorr*rcpd
          end do
          lcor = lqm + 1
        end if
!
!       below lqm correct both temperature and moisture
!
        do l = lcor , lb
          tskl = trefk(l)*apek(l)/apesk(l)
          dhdt = qrefk(l)*a23m4l/(tskl-c4les)**d_two + cpd
          trefk(l) = hcorr/dhdt + trefk(l)
          thskl = trefk(l)*apek(l)
          qrefk(l) = pq0/psk(l)*dexp(c3les*(thskl-tzero*apesk(l)) / &
                     (thskl-c4les*apesk(l)))
        end do
      end do
      do l = 1 , kz
        thvref(l) = trefk(l)*apek(l)*(qrefk(l)*d608+h1)
      end do
!
!     heating, moistening, precipitation
!
      do l = ltpk , lb
        tkl = tk(l)
        diftl = (trefk(l)-tkl)*tauk
        difql = (qrefk(l)-qk(l))*tauk
        avrgtl = (tkl+tkl+diftl)
        dentpy = (diftl*cpd+difql*wlhv)*dflev(l)/avrgtl + dentpy
        avrgt = avrgtl*dflev(l) + avrgt
        preck = dflev(l)*diftl + preck
        dift(l) = diftl
        difq(l) = difql
      end do
      dentpy = dentpy + dentpy
      avrgt = avrgt/(sumdp+sumdp)
      if ( dentpy < epsntp .or. preck <= d_zero ) then
        if ( oct90 ) then
          cldefi(i,j) = efimn
        else
          cldefi(i,j) = efimn*xsm(j,i) + stefi*(h1-xsm(j,i))
        end if
        ztop = z0(j,i,lbot(j,i)) + zsh - 0.000001D0
        do l = 1 , lb
          if ( z0(j,i,l) >= ztop ) ltop(j,i) = l + 1
        end do
        prtop(j,i) = pk(ltop(j,i))
!
!       cloud must be at least 2 layers thick
!
        if ( lbot(j,i)-ltop(j,i) < 2 ) ltop(j,i) = lbot(j,i) - 2
        prtop(j,i) = pk(ltop(j,i))
        cldhgt(j,i) = z0(j,i,ltop(j,i)) - z0(j,i,lbot(j,i))
        nswap = nswap + 1
        cycle
      end if
!
!     deep convection otherwise
!
      total_precip_points = total_precip_points + 1
!     keep the land value of efi equal to 1 until precip surpasses
!     a threshold value, currently set to 0.25 inches per 24 hrs
      pthrs = cthrs/sfcps(i,j)
      drheat = (preck*xsm(j,i)+dmax1(epsp,preck-pthrs)*(h1-xsm(j,i)))*cpd/avrgt
      efi = efifc*dentpy/drheat
!
!     unified or separate land/sea conv.
!
      if ( .not.(oct90) ) then
        efi = cldefi(i,j)*fcp + efi*fcc1
      else if ( unis ) then
        efi = cldefi(i,j)*fcp + efi*fcc1
      else if ( .not.unil ) then
        efi = (cldefi(i,j)*fcp+efi*fcc1)*xsm(j,i) + h1 - xsm(j,i)
      else
        efi = h1
      end if
!
      if ( efi > h1 ) efi = h1
      if ( efi < efimn ) efi = efimn
      cldefi(i,j) = efi
!
      fefi = efmnt + slope*(efi-efimn)
!
      preck = preck*fefi
!
!     update precipitation, temperature & moisture
!
      prainx = d_half*((sfcps(i,j)*d_1000*preck*cprlg)*d_100)
      if ( prainx > dlowval ) then
        rainc(i,j) = rainc(i,j) + prainx
!       precipitation rate for bats (mm/s)
        if ( ktau == 0 ) then
          lmpcpc(i,j) = lmpcpc(i,j) + prainx/dtcum
        else
          lmpcpc(i,j) = lmpcpc(i,j) + (prainx*dtcum)/aprdiv
        end if
      end if
      do l = ltpk , lb
        tmod(j,i,l) = dift(l)*fefi/dtcum
        qqmod(j,i,l) = difq(l)*fefi/dtcum
      end do
!
!dcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd
!dcdcdcdcdcdc  end of deep convection  dcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd
!dcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcdcd
!
    end do
    ndeep = 0
    do i = istart , iend
      do j = jstart , jend
        ltpk = ltop(j,i)
        lbtk = lbot(j,i)
        ptpk = prtop(j,i)
        if ( cldhgt(j,i) >= zsh ) then
          ndeep = ndeep + 1
          ndepth = lb - ltpk
          ntopd(ltpk) = ntopd(ltpk) + 1
          nbotd(lb) = nbotd(lb) + 1
          if ( ndepth > 0 ) ndpthd(ndepth) = ndpthd(ndepth) + 1
        end if
      end do
    end do
!
!   gather shallow convection points
!
    khshal = 0
    ndstn = 0
    ndstp = 0
    do i = istart , iend
      do j = jstart , jend
        if ( cldhgt(j,i) >= zno .and. ltop(j,i) <= lbot(j,i)-2 ) then
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
    do n = 1 , khshal
      i = ishal(n)
      j = jshal(n)
      do k = 1 , kz
        tkl = t(j,i,k)
        tk(k) = tkl
        trefk(k) = tkl
        qkl = q(j,i,k)
        qk(k) = qkl
        qrefk(k) = qkl
        qsatk(k) = qkl
        pkl = (hlev(k)*sfcps(i,j)+ptop)*d_1000
        pk(k) = pkl
        apekl = ape(j,i,k)
        apek(k) = apekl
        thvmkl = tkl*apekl*(qkl*d608+h1)
        thvref(k) = thvmkl
        pdp(k) = pk(k) - pmn
      end do
!
!     find kdp...kdp(k) is the model level closest to 65 mb (pmn) above k;
!     this is the depth over which relative humidity drop is measured to
!     estimate shallow cloud top... see do 545...
!
      do kk = kz , 1 , -1
        pflag = dabs(pk(kz)-pdp(kk))
        do k = kz - 1 , 1 , -1
          pdiffk = dabs(pk(k)-pdp(kk))
          if ( pdiffk < pflag ) then
            pflag = pdiffk
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
!     search for shallow cloud top
!
      lbtk = lbot(j,i)
      ltsh = lbtk
      lbm1 = lbtk - 1
      ztop = z0(j,i,lbot(j,i)) + zsh - 0.000001D0
!
!     cloud top is level just above pbtk-psh
!
      do l = 1 , kz
        if ( z0(j,i,l) >= ztop ) ltpk = l
      end do
      ptpk = pk(ltpk)
!
!     highest level allowed is level just below pshu
!
      if ( ptpk <= pshu ) then
        do l = 1 , kz
          if ( pk(l) <= pshu ) lshu = l + 1
        end do
        ltpk = lshu
        ptpk = pk(ltpk)
      end if
      ltp1 = ltpk + 1
      do l = ltpk , lbtk
        if ( l >= ml(j,i) ) then
          es = aliq*dexp((bliq*tk(l)-cliq)/(tk(l)-dliq))
        else
          es = aice*dexp((bice*tk(l)-cice1)/(tk(l)-dice))
        end if
        qsatk(l) = ep2*es/(pk(l)-es)
      end do
      do l = ltp1 , lbm1
        rhl = qk(l)/qsatk(l)
        rhh = qk(kdp(l))/qsatk(kdp(l))
        if ( rhh+rhf < rhl ) ltsh = l
      end do
!
      ltop(j,i) = ltsh
      prtop(j,i) = pk(ltsh)
      ltp1 = ltsh
      ltpk = ltsh - 1
      cldhgt(j,i) = z0(j,i,ltop(j,i)) - z0(j,i,lbot(j,i))
!     if cloud is not at least 90 mb or 3 sigma layers deep, then no cloud
      if ( cldhgt(j,i) < zno .or. ltop(j,i) > lbot(j,i)-2 ) then
        ltop(j,i) = lbot(j,i)
        prtop(j,i) = pbot(j,i)
        cycle
      end if
!     scaling potential temperature & table index at top
      thtpk = t(j,i,ltp1)*ape(j,i,ltp1)
      pkl = (hlev(ltp1)*sfcps(i,j)+ptop)*d_1000
      ee = pkl*q(j,i,ltp1)/(ep2+q(j,i,ltp1))
      tdpt = d_one/(rtzero-rwat/wlhv*dlog(ee/611.D0))
      tdpt = dmin1(tdpt,t(j,i,ltp1))
      tlcl = tdpt - (0.212D0+1.571D-3*(tdpt-tzero) - &
             4.36D-4*(t(j,i,ltp1)-tzero))*(t(j,i,ltp1)-tdpt)
      ptpk = h10e5*(thtpk/tlcl)**cporng
      dpmix = ptpk - psp(j,i)
      if ( dabs(dpmix) < h3000 ) dpmix = -h3000
!
!     temperature propfile slope
!
      smix = (thtpk-thbt(j,i))/dpmix*stabs
      do l = ltp1 , lbtk
        ivi = min(ltp1+lbtk-l,kz-1)
        trefk(ivi) = ((pk(ivi)-pk(ivi+1))*smix + &
                     trefk(ivi+1)*apek(ivi+1))/apek(ivi)
      end do
!
!     temperature reference profile correction
!
      sumdt = d_zero
      sumdp = d_zero
      do l = ltp1 , lbtk
        sumdt = (tk(l)-trefk(l))*dflev(l) + sumdt
        sumdp = sumdp + dflev(l)
      end do
!
      rdpsum = d_one/sumdp
      fpk(lbtk) = trefk(lbtk)
      tcorr = sumdt*rdpsum
      do l = ltp1 , lbtk
        trfkl = trefk(l) + tcorr
        trefk(l) = trfkl
        fpk(l) = trfkl
      end do
!
!     humidity profile equations
!
      psum = d_zero
      qsum = d_zero
      potsum = d_zero
      qotsum = d_zero
      otsum = d_zero
      dst = d_zero
      fptk = fpk(ltp1)
      do l = ltp1 , lbtk
        dpkl = fpk(l) - fptk
        psum = dpkl*dflev(l) + psum
        qsum = qk(l)*dflev(l) + qsum
        rtbar = 2.0D0/(trefk(l)+tk(l))
        otsum = dflev(l)*rtbar + otsum
        potsum = dpkl*rtbar*dflev(l) + potsum
        qotsum = qk(l)*rtbar*dflev(l) + qotsum
        dst = (trefk(l)-tk(l))*rtbar*dflev(l) + dst
      end do
!
      psum = psum*rdpsum
      qsum = qsum*rdpsum
      rotsum = d_one/otsum
      potsum = potsum*rotsum
      qotsum = qotsum*rotsum
      dst = dst*rotsum*cpd/wlhv
!
!     ensure positive entropy change
!
      if ( dst > d_zero ) then
        prtop(j,i) = pbot(j,i)
        ltop(j,i) = lbot(j,i)
        ndstp = ndstp + 1
        cycle
      else
        dstq = dst*epsdn
      end if
!
!     check for isothermal atmosphere
!
      den = potsum - psum
      if ( -den/psum < 0.00005D0 ) then
        ltop(j,i) = lbot(j,i)
        prtop(j,i) = pbot(j,i)
        cycle
      else
!
!       slope of the reference humidity profile
!
        dqref = (qotsum-dstq-qsum)/den
      end if
!
!     humidity doesn`t increase with height
!
      if ( dqref < d_zero ) then
        ltop(j,i) = lbot(j,i)
        prtop(j,i) = pbot(j,i)
        cycle
      end if
!
!     humidity at the cloud top
!
      qrftp = qsum - dqref*psum
!
!     humidity profile
!
      do l = ltp1 , lbtk
        qrfkl = (fpk(l)-fptk)*dqref + qrftp
!       supersaturation not allowed
        qnew = (qrfkl-qk(l))*tauk + qk(l)
        if ( qnew > qsatk(l)*stresh ) then
          ltop(j,i) = lbot(j,i)
          prtop(j,i) = pbot(j,i)
          go to 100
        end if
        thvref(l) = trefk(l)*apek(l)*(qrfkl*d608+h1)
        qrefk(l) = qrfkl
      end do
!
!     eliminate impossible slopes (betts, dtheta/dq)
!
      do l = ltp1 , lbtk
        dtdeta = (thvref(l-1)-thvref(l))/(hlev(l)-hlev(l-1))
        if ( dtdeta < epsth ) then
          ltop(j,i) = lbot(j,i)
          prtop(j,i) = pbot(j,i)
          go to 100
        end if
      end do
      if ( dst > d_zero ) then
        ndstp = ndstp + 1
      else
        ndstn = ndstn + 1
      end if
      dentpy = d_zero
      do l = ltp1 , lbtk
        dentpy = ((trefk(l)-tk(l))*cpd+(qrefk(l)-qk(l))*wlhv) / &
                  (tk(l)+trefk(l))*dflev(l) + dentpy
      end do
!
!     relaxation towards reference profiles
!
      iconss = iconss + 1
      do l = ltp1 , lbtk
        tmod(j,i,l) = (trefk(l)-tk(l))/trel
        qqmod(j,i,l) = (qrefk(l)-qk(l))/trel
      end do
!
!scscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscs
!scscscscscsc  end of shallow convection   scscscscscscscscscscscscscscs
!scscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscscs
!
    end do

   100  continue

    nshal = 0
    do i = istart , iend
      do j = jstart , jend
        ltpk = ltop(j,i)
        lbtk = lbot(j,i)
        ptpk = prtop(j,i)
!       no shallow convection if cloud is not at least 90 mb or 3 sigma
!       layers deep
        if ( cldhgt(j,i) >= zno ) then
          if ( cldhgt(j,i) < zsh ) then
            nshal = nshal + 1
            ntops(ltpk) = ntops(ltpk) + 1
            nbots(lbtk) = nbots(lbtk) + 1
            ndepth = lbtk - ltpk
            if ( ndepth > 0 ) ndpths(ndepth) = ndpths(ndepth) + 1
          end if
!         find cloud fractional cover and liquid water content
          kbaseb = min0(lbtk,kzm2)
          if ( ltpk <= kbaseb ) then
            kclth = kbaseb - ltpk + 1
            akclth = d_one/dble(kclth)
            do k = ltpk , kbaseb
              rcldlwc(j,i,k) = cllwcv
              rcldfra(j,i,k) = d_one - (d_one-clfrcv)**akclth
            end do
          end if
          if ( lchem ) then
            icumtop(j,i) = ltpk
            icumbot(j,i) = kbaseb
          end if
        end if
      end do
    end do
    do k = 1 , kz
      do i = istart , iend
        do j = jstart , jend
          tten(i,k,j)  = tten(i,k,j)  + tmod(j,i,k) *sfcps(i,j)
          qvten(i,k,j) = qvten(i,k,j) + qqmod(j,i,k)*sfcps(i,j)
        end do
      end do
    end do

  end subroutine bmpara
!
! Look up table (calculated version)
!
  subroutine lutbl(ptop)
!
    implicit none
!
    real(8) , parameter :: eps = 2.0D-12 ! little number

!
    real(8) :: ptop
    intent (in) ptop
!
    real(8) :: ape , dp , dqs , dth , dthe , p , pt , qs , qs0k , &
               sqsk , sthek , th , the0k
    real(8) , dimension(jtb) :: pnew , pold , qsnew , qsold ,  &
               thenew , theold , tnew , told , y2p , y2t
    integer :: kp , kpm , kpm1 , kth , kthm , kthm1
    real(8) , parameter :: thl = 210.0D0
    real(8) , parameter :: thh = 385.0D0
    real(8) , parameter :: ph = 105000.0D0
!
!   coarse look-up table for saturation point
!
    pt = ptop*d_1000
!   ptop in pascal
   
    kthm = jtb
    kpm = itb
    kthm1 = kthm - 1
    kpm1 = kpm - 1
!
    pl = pt
!
    dth = (thh-thl)/dble(kthm-1)
    dp = (ph-pl)/dble(kpm-1)
!
    th = thl - dth
!
!-----------------------------------------------------------------------
!
    do kth = 1 , kthm
      th = th + dth
      p = pl - dp
      do kp = 1 , kpm
        p = p + dp
        ape = (100000.0D0/p)**(rovcp)
        qsold(kp) = pq0/p*dexp(c3les*(th-tzero*ape)/(th-c4les*ape))
        pold(kp) = p
      end do
!
      qs0k = qsold(1)
      sqsk = qsold(kpm) - qsold(1)
      qsold(1) = d_zero
      qsold(kpm) = d_one
!
      do kp = 2 , kpm1
        qsold(kp) = (qsold(kp)-qs0k)/sqsk
!       fix due to cyber half prec. limitation
        if ( (qsold(kp)-qsold(kp-1)) < eps ) then
          qsold(kp) = qsold(kp-1) + eps
        end if
      end do
!
      qsnew(1) = d_zero
      qsnew(kpm) = d_one
      dqs = d_one/dble(kpm-1)
!
      do kp = 2 , kpm1
        qsnew(kp) = qsnew(kp-1) + dqs
      end do
!
      y2p(1) = d_zero
      y2p(kpm) = d_zero
!
      call spline(kpm,qsold,pold,y2p,kpm,qsnew,pnew)
!
    end do
!
!   coarse look-up table for t(p) from constant the
!
    p = pl - dp
    do kp = 1 , kpm
      p = p + dp
      th = thl - dth
      do kth = 1 , kthm
        th = th + dth
        ape = (100000.0D0/p)**(rovcp)
        qs = pq0/p*dexp(c3les*(th-tzero*ape)/(th-c4les*ape))
        told(kth) = th/ape
        theold(kth) = th*dexp(eliwv*qs/(cpd*told(kth)))
      end do
!
      the0k = theold(1)
      sthek = theold(kthm) - theold(1)
      theold(1) = d_zero
      theold(kthm) = d_one
!
      do kth = 2 , kthm1
        theold(kth) = (theold(kth)-the0k)/sthek
!       fix due to cyber half prec. limitation
        if ( (theold(kth)-theold(kth-1)) < eps ) then
          theold(kth) = theold(kth-1) + eps
        end if
      end do
!
      thenew(1) = d_zero
      thenew(kthm) = d_one
      dthe = d_one/dble(kthm-1)
!
      do kth = 2 , kthm1
        thenew(kth) = thenew(kth-1) + dthe
      end do
!
      y2t(1) = d_zero
      y2t(kthm) = d_zero
!
      call spline(kthm,theold,told,y2t,kthm,thenew,tnew)
!
    end do
!
  end subroutine lutbl
!
!*****************************************************************
!                                                                *
!  This is a one-dimensional cubic spline fitting routine        *
!  programed for a small scalar machine.                         *
!                                                                *
!  Programer: Z. Janjic, Yugoslav Fed. Hydromet. Inst., Beograd  *
!                                                                *
!  nold - number of given values of the function.  must be ge 3. *
!  xold - locations of the points at which the values of the     *
!         function are given.  must be in ascending order.       *
!  yold - the given values of the function at the points xold.   *
!  y2   - the second derivatives at the points xold.  if natural *
!         spline is fitted y2(1)=0. and y2(nold)=0. must be      *
!         specified.                                             *
!  nnew - number of values of the function to be calculated.     *
!  xnew - locations of the points at which the values of the     *
!         function are calculated.  xnew(k) must be ge xold(1)   *
!         and le xold(nold).                                     *
!  ynew - the values of the function to be calculated.           *
!  p, q - auxiliary vectors of the length nold-2.                *
!                                                                *
!*****************************************************************
!
  subroutine spline(nold,xold,yold,y2,nnew,xnew,ynew)
 
    implicit none
!
    integer :: nnew , nold
    real(8) , dimension(nold) :: xold , yold , y2
    real(8) , dimension(nnew) :: xnew, ynew
    intent (in) nnew , nold , xnew , xold , yold
    intent (out) ynew
    intent (inout) y2
!
    real(8) , dimension(nold-2) :: p , q
    real(8) :: ak , bk , ck , den , dx , dxc , dxl , dxr , dydxl ,    &
               dydxr , rdx , rtdxc , x , xk , xsq , y2k , y2kp1
    integer :: k , k1 , k2 , kold , noldm1
!
!-----------------------------------------------------------------------
!
    ak = d_zero
    bk = d_zero
    ck = d_zero
    noldm1 = nold - 1
!
    dxl = xold(2) - xold(1)
    dxr = xold(3) - xold(2)
    dydxl = (yold(2)-yold(1))/dxl
    dydxr = (yold(3)-yold(2))/dxr
    rtdxc = d_half/(dxl+dxr)
!
    p(1) = rtdxc*(6.0D0*(dydxr-dydxl)-dxl*y2(1))
    q(1) = -rtdxc*dxr
!
    if ( nold == 3 ) then
      k = noldm1
    else
      k = 3
      do
!
        dxl = dxr
        dydxl = dydxr
        dxr = xold(k+1) - xold(k)
        dydxr = (yold(k+1)-yold(k))/dxr
        dxc = dxl + dxr
        den = d_one/(dxl*q(k-2)+dxc+dxc)
!
        p(k-1) = den*(6.0D0*(dydxr-dydxl)-dxl*p(k-2))
        q(k-1) = -den*dxr
!
        k = k + 1
        if ( k >= nold ) then
          k = noldm1
          exit
        end if
      end do
    end if
!
    do
!
      y2(k) = p(k-1) + q(k-1)*y2(k+1)
!
      k = k - 1
      if ( k <= 1 ) then
!-----------------------------------------------------------------------
        k1 = 1
        exit
      end if
    end do
!
 100  continue

    xk = xnew(k1)
!
    do k2 = 2 , nold
      if ( xold(k2) > xk ) then
        kold = k2 - 1
!
        if ( k1 == 1 ) go to 200
        if ( k /= kold ) go to 200
        go to 300
      end if
    end do
    ynew(k1) = yold(nold)
    go to 400
!
   200  continue
    k = kold
!
    y2k = y2(k)
    y2kp1 = y2(k+1)
    dx = xold(k+1) - xold(k)
    rdx = d_one/dx
    ak = (d_five/d_three)*rdx*(y2kp1-y2k)
    bk = d_half*y2k
    ck = rdx*(yold(k+1)-yold(k))-(d_five/d_three)*dx*(y2kp1+y2k+y2k)
!
   300  continue
    x = xk - xold(k)
    xsq = x*x
!
    ynew(k1) = ak*xsq*x + bk*xsq + ck*x + yold(k)
!
   400  continue
!
    k1 = k1 + 1
    if ( k1 <= nnew ) go to 100
!-----------------------------------------------------------------------
  end subroutine spline
!
! Calculates tpfc
!
  function tpfc(press,thetae,tgs,rl,qs,pi)
 
    implicit none
!
    real(8) :: pi , press , qs , rl , tgs , thetae
    real(8) :: tpfc
    intent (in) pi , press , rl , tgs , thetae
    intent (inout) qs
!
    real(8) :: dtx , es , f1 , fo , rlocpd , rlorw , rp , t1 , tguess
!
!   iteratively extract temperature from equivalent potential temperature.
!
    rlorw = rl/rwat
    rlocpd = rl*rcpd
    rp = thetae/pi
    es = 611.0D0*dexp(rlorw*(rtzero-d_one/tgs))
    qs = ep2*es/(press-es)
    fo = tgs*dexp(rlocpd*qs/tgs) - rp
    t1 = tgs - d_half*fo
    tguess = tgs
    do
      es = 611.0D0*dexp(rlorw*(rtzero-d_one/t1))
      qs = ep2*es/(press-es)
      f1 = t1*dexp(rlocpd*qs/t1) - rp
      if ( dabs(f1) < 0.1D0 ) then
        tpfc = t1
        exit
      else
        dtx = f1*(t1-tguess)/(f1-fo)
        tguess = t1
        fo = f1
        t1 = t1 - dtx
      end if
    end do
  end function tpfc
!
end module mod_cu_bm

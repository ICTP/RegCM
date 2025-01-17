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

module mod_pbl_shinhong

  use mod_realkinds , only : rkx
  use mod_regcm_types , only : mod_2_pbl , pbl_2_mod
  use mod_dynparam , only : jci1 , jci2 , ici1 , ici2 , ntr
  use mod_dynparam , only : kz , kzm1 , kzp1 , kzp2 , ntr
  use mod_dynparam , only : idynamic
  use mod_runparams , only : ichem , ipptls , iqv , iqc , iqi , dt , dx
  use mod_constants , only : egrav , regrav , cpd , rcpd , rdry , rwat
  use mod_constants , only : vonkar , ep1 , wlhv , p00 , rovcp , d_one
  use mod_memutil , only : getmem1d , getmem2d

  implicit none

  private

  public :: init_shinhong_pbl , shinhong_pbl

  integer :: numbl , nj , ni , ndiff , ichs

  real(rkx) , parameter :: vconvc = 1.0_rkx
  real(rkx) , parameter :: czo = 0.0185_rkx
  real(rkx) , parameter :: ozo = 1.59e-5_rkx
  real(rkx) , dimension(0:1000) :: psimtb , psihtb

  real(rkx) , dimension(:,:) , pointer :: u2d , v2d , utnp , vtnp
  real(rkx) , dimension(:,:) , pointer :: t2d , ttnp , th2d , p2d , dz2d , pi2d
  real(rkx) , dimension(:,:) , pointer :: p2di
  real(rkx) , dimension(:,:) , pointer :: tke2d
  real(rkx) , dimension(:,:) , pointer :: qtrac , qtnp
  real(rkx) , dimension(:) , pointer :: psfc , hfx , qfx , ust , znt
  real(rkx) , dimension(:) , pointer :: wspd , psim , psih , br
  real(rkx) , dimension(:) , pointer :: hpbl , dusfc , dvsfc , dtsfc , dqsfc
  real(rkx) , dimension(:) , pointer :: wstar , delta , wspd10
  real(rkx) , dimension(:) , pointer :: corf , za
  real(rkx) , dimension(:) , pointer :: govrth , dtg , rah , rpfac
  integer , dimension(:) , pointer :: xland
  integer , dimension(:) , pointer :: kpbl

  contains

  subroutine init_shinhong_pbl( )
    implicit none
    nj = jci2-jci1+1
    ni = ici2-ici1+1
    numbl = nj * ni
    ndiff = 2
    if ( ipptls > 1 ) then
      ndiff = ndiff + 1
    end if
    if ( ichem == 1 ) then
      ichs = ndiff
      ndiff = ndiff + ntr
    else
      ichs = 0
    end if
    call sfclayinit
    call getmem2d(u2d,1,numbl,1,kz,'shinhong_pbl')
    call getmem2d(v2d,1,numbl,1,kz,'shinhong_pbl')
    call getmem2d(t2d,1,numbl,1,kz,'shinhong_pbl')
    call getmem2d(th2d,1,numbl,1,kz,'shinhong_pbl')
    call getmem2d(p2d,1,numbl,1,kz,'shinhong_pbl')
    call getmem2d(dz2d,1,numbl,1,kz,'shinhong_pbl')
    call getmem2d(pi2d,1,numbl,1,kz,'shinhong_pbl')
    call getmem2d(p2di,1,numbl,1,kzp1,'shinhong_pbl')
    call getmem2d(tke2d,1,numbl,1,kz,'shinhong_pbl')
    call getmem2d(qtrac,1,numbl,1,kz*ndiff,'shinhong_pbl')
    call getmem2d(utnp,1,numbl,1,kz,'shinhong_pbl')
    call getmem2d(vtnp,1,numbl,1,kz,'shinhong_pbl')
    call getmem2d(ttnp,1,numbl,1,kz,'shinhong_pbl')
    call getmem2d(qtnp,1,numbl,1,kz*ndiff,'shinhong_pbl')
    call getmem1d(psfc,1,numbl,'shinhong_pbl')
    call getmem1d(hfx,1,numbl,'shinhong_pbl')
    call getmem1d(qfx,1,numbl,'shinhong_pbl')
    call getmem1d(ust,1,numbl,'shinhong_pbl')
    call getmem1d(znt,1,numbl,'shinhong_pbl')
    call getmem1d(wspd,1,numbl,'shinhong_pbl')
    call getmem1d(wspd10,1,numbl,'shinhong_pbl')
    call getmem1d(psim,1,numbl,'shinhong_pbl')
    call getmem1d(psih,1,numbl,'shinhong_pbl')
    call getmem1d(br,1,numbl,'shinhong_pbl')
    call getmem1d(hpbl,1,numbl,'shinhong_pbl')
    call getmem1d(dusfc,1,numbl,'shinhong_pbl')
    call getmem1d(dvsfc,1,numbl,'shinhong_pbl')
    call getmem1d(dtsfc,1,numbl,'shinhong_pbl')
    call getmem1d(dqsfc,1,numbl,'shinhong_pbl')
    call getmem1d(wstar,1,numbl,'shinhong_pbl')
    call getmem1d(delta,1,numbl,'shinhong_pbl')
    call getmem1d(corf,1,numbl,'shinhong_pbl')
    call getmem1d(za,1,numbl,'shinhong_pbl')
    call getmem1d(govrth,1,numbl,'shinhong_pbl')
    call getmem1d(dtg,1,numbl,'shinhong_pbl')
    call getmem1d(rah,1,numbl,'shinhong_pbl')
    call getmem1d(rpfac,1,numbl,'shinhong_pbl')
    call getmem1d(xland,1,numbl,'shinhong_pbl')
    call getmem1d(kpbl,1,numbl,'shinhong_pbl')
  end subroutine init_shinhong_pbl

  subroutine shinhong_pbl(m2p,p2m)
    implicit none
    type(mod_2_pbl) , intent(in) :: m2p
    type(pbl_2_mod) , intent(inout) :: p2m
    !
    ! u2d         3d u-velocity interpolated to theta points (m/s)
    ! v2d         3d v-velocity interpolated to theta points (m/s)
    ! th2d        3d potential temperature (k)
    ! t2d         temperature (k)
    ! qtrac       3d water vapor mixing ratio (kg/kg)
    !             3d cloud mixing ratio (kg/kg)                      +
    !             3d ice mixing ratio (kg/kg) (ipptls > 1)           +
    !             3d tracer mixing ratio (kg/kg) (ichem == 1, ntr)
    ! p2d         3d pressure (pa)
    ! p2di        3d pressure (pa) at interface level
    ! pi2d        3d exner function (dimensionless)
    ! dz2d        dz between full levels (m)
    ! tke2d       3d diagnostic TKE
    ! psfc        pressure at the surface (pa)
    ! hfx	  upward heat flux at the surface (w/m^2)
    ! qfx	  upward moisture flux at the surface (kg/m^2/s)
    ! ust	  u* in similarity theory (m/s)
    ! znt         roughness length (m)
    ! hpbl        pbl height (m)
    ! psim        similarity stability function for momentum
    ! psih        similarity stability function for heat
    ! xland       land mask (1 for land, 0 for water)
    ! wspd        wind speed at lowest model level (m/s)
    ! br          bulk richardson number in surface layer
    ! u10         u-wind speed at 10 m (m/s)
    ! v10         v-wind speed at 10 m (m/s)
    !
    integer :: i , j , k , kk , it , ibin

    if ( idynamic == 3 ) then
      rpfac = 1.0
    else
      do i = ici1 , ici2
        do j = jci1 , jci2
          ibin = (i-ici1)*nj+(j-jci1+1)
          rpfac(ibin) = m2p%psb(j,i)
        end do
      end do
    end if

    do i = ici1 , ici2
      do j = jci1 , jci2
        ibin = (i-ici1)*nj+(j-jci1+1)
        psfc(ibin) = m2p%patmf(j,i,kzp1)
        hfx(ibin) = m2p%hfx(j,i)
        qfx(ibin) = m2p%qfx(j,i)
        ust(ibin) = m2p%ustar(j,i)
        znt(ibin) = max(m2p%zo(j,i),2.0e-4_rkx)
        xland(ibin) = m2p%ldmsk(j,i)
        wspd(ibin) = sqrt(m2p%uxatm(j,i,kz)**2 + m2p%vxatm(j,i,kz)**2)
        wspd10(ibin) = sqrt(m2p%u10m(j,i)**2 + m2p%v10m(j,i)**2)
        za(ibin) = m2p%za(j,i,kz)
        corf(ibin) = m2p%coriol(j,i)
        govrth(ibin) = egrav/m2p%thatm(j,i,kz)
        rah(ibin) = m2p%rah1(j,i)
        dtg(ibin) = m2p%thatm(j,i,kz) - &
          (m2p%tg(j,i)*(1.0_rkx+ep1*m2p%qxatm(j,i,kz,iqv)))
        hpbl(ibin) = za(ibin)
        br(ibin) = govrth(ibin)*za(ibin)*dtg(ibin) / &
                 max(wspd(ibin)*wspd(ibin),0.5_rkx)
      end do
    end do
    do k = 1 , kz
      kk = kzp1 - k
      do i = ici1 , ici2
        do j = jci1 , jci2
          ibin = (i-ici1)*nj+(j-jci1+1)
          p2d(ibin,kk) = m2p%patm(j,i,k)
          u2d(ibin,kk) = m2p%uxatm(j,i,k)
          v2d(ibin,kk) = m2p%vxatm(j,i,k)
          dz2d(ibin,kk) = m2p%dzq(j,i,k)
          t2d(ibin,kk) = m2p%tatm(j,i,k)
          th2d(ibin,kk) = m2p%thatm(j,i,k)
          pi2d(ibin,kk) = (m2p%patm(j,i,k)/p00)**rovcp
          tke2d(ibin,kk) = p2m%tkepbl(j,i,k)
          qtrac(ibin,kk) = m2p%qxatm(j,i,k,iqv)
          qtrac(ibin,kk+kz) = m2p%qxatm(j,i,k,iqc)
        end do
      end do
    end do
    do k = 1 , kzp1
      kk = kzp2 - k
      do i = ici1 , ici2
        do j = jci1 , jci2
          ibin = (i-ici1)*nj+(j-jci1+1)
          p2di(ibin,kk) = m2p%patmf(j,i,k)
        end do
      end do
    end do
    if ( ipptls > 1 ) then
      do k = 1 , kz
        kk = 2*kz + kzp1 - k
        do i = ici1 , ici2
          do j = jci1 , jci2
            ibin = (i-ici1)*nj+(j-jci1+1)
            qtrac(ibin,kk) = m2p%qxatm(j,i,k,iqi)
          end do
        end do
      end do
    end if
    if ( ichem == 1 ) then
      do it = 1 , ntr
        do k = 1 , kz
          kk = (ichs+it-1)*kz + kzp1 - k
          do i = ici1 , ici2
            do j = jci1 , jci2
              ibin = (i-ici1)*nj+(j-jci1+1)
              qtrac(ibin,kk) = m2p%chib(j,i,k,it)
            end do
          end do
        end do
      end do
    end if
    !
    call sfclay(numbl,br,znt,za,ust,govrth,dtg,rah,psim,psih)
    !
    call shinhong2d(numbl,ux=u2d,vx=v2d,tx=t2d,qx=qtrac,p2d=p2d,         &
                    p2di=p2di,pi2d=pi2d,utnp=utnp,vtnp=vtnp,ttnp=ttnp,   &
                    qtnp=qtnp,dz8w2d=dz2d,psfcpa=psfc,znt=znt,ust=ust,   &
                    hpbl=hpbl,psim=psim,psih=psih,xland=xland,hfx=hfx,   &
                    qfx=qfx,wspd=wspd,br=br,dusfc=dusfc,dvsfc=dvsfc,     &
                    dtsfc=dtsfc,dqsfc=dqsfc,dt=dt,rcl=d_one,kpbl1d=kpbl, &
                    wstar=wstar,delta=delta,wspd10=wspd10,tke=tke2d,     &
                    corf=corf,dx=dx,dy=dx)
    do i = ici1 , ici2
      do j = jci1 , jci2
        ibin = (i-ici1)*nj+(j-jci1+1)
        p2m%zpbl(j,i) = hpbl(ibin)
        p2m%kpbl(j,i) = kpbl(ibin)
      end do
    end do
    do k = 1 , kz
      kk = kzp1 - k
      do i = ici1 , ici2
        do j = jci1 , jci2
          ibin = (i-ici1)*nj+(j-jci1+1)
          p2m%tten(j,i,k) = p2m%tten(j,i,k) + &
            (rpfac(ibin)*ttnp(ibin,kk)/pi2d(ibin,kk))
          p2m%uxten(j,i,k) = p2m%uxten(j,i,k)+(rpfac(ibin)*utnp(ibin,kk))
          p2m%vxten(j,i,k) = p2m%vxten(j,i,k)+(rpfac(ibin)*vtnp(ibin,kk))
          p2m%qxten(j,i,k,iqv) = p2m%qxten(j,i,k,iqv) + &
            (rpfac(ibin)*qtnp(ibin,kk))
          p2m%qxten(j,i,k,iqc) = p2m%qxten(j,i,k,iqc) + &
            (rpfac(ibin)*qtnp(ibin,kz+kk))
          p2m%tkepbl(j,i,k) = tke2d(ibin,kk)
        end do
      end do
    end do
    if ( ipptls > 1 ) then
      do k = 1 , kz
        kk = 2*kz + kzp1 - k
        do i = ici1 , ici2
          do j = jci1 , jci2
            ibin = (i-ici1)*nj+(j-jci1+1)
            p2m%qxten(j,i,k,iqi) = p2m%qxten(j,i,k,iqi) + &
              (rpfac(ibin)*qtnp(ibin,kk))
          end do
        end do
      end do
    end if
    if ( ichem == 1 ) then
      do it = 1 , ntr
        do k = 1 , kz
          kk = (ichs+it-1)*kz + kzp1 - k
          do i = ici1 , ici2
            do j = jci1 , jci2
              ibin = (i-ici1)*nj+(j-jci1+1)
              p2m%chiten(j,i,k,it) = p2m%chiten(j,i,k,it) + &
                (rpfac(ibin)*qtnp(ibin,kk))
            end do
          end do
        end do
      end do
    end if
  end subroutine shinhong_pbl

  subroutine shinhong2d(nbl,ux,vx,tx,qx,p2d,p2di,pi2d,utnp,vtnp,ttnp,qtnp, &
         dz8w2d,psfcpa,znt,ust,hpbl,psim,psih,xland,hfx,qfx,wspd,br,dusfc, &
         dvsfc,dtsfc,dqsfc,dt,rcl,kpbl1d,wstar,delta,tke,corf,wspd10,      &
         dx,dy)
    implicit none
    !
    !-----------------------------------------------------------------------
    !
    ! The shinhongpbl (Shin and Hong 2015) is based on the les study of Shin
    ! and Hong (2013). The major ingredients of the shinhongpbl are
    ! 1) The prescribed nonlocal heat transport profile fit to the LES and
    ! 2) Inclusion of explicit scale dependency functions for vertical
    !    transport in convective PBL.
    ! The shinhongpbl works at the gray zone resolution of convective pbl.
    ! Note that Honnert et al. (2011) first suggested explicit scale
    ! dependency function, and Shin and Hong (2013) further classified the
    ! function by stability (u*/w*) in convective pbl and calculated the
    ! function for nonlocal and local transport separately.
    ! Vertical mixing in the stable boundary layer and free atmosphere follows
    ! Hong (2010) and Hong et al. (2006), same as the ysupbl scheme.
    !
    ! shinhongpbl:
    !     coded and implemented by Hyeyum Hailey Shin (NCAR)
    !              Summer 2014
    !
    ! ysupbl:
    !     coded by Song-You Hong (Yonsei university) and implemented by
    !              Song-You Hong (Yonsei university) and Jimy Dudhia (NCAR)
    !              Summer 2002
    !
    ! References:
    !    Shin and Hong (2015) Mon. Wea. Rev.
    !    Shin and Hong (2013) J. Atmos. Sci.
    !    Honnert, Masson, and couvreux (2011) J. Atmos. Sci.
    !    Hong (2010) Quart. J. Roy. Met. Soc
    !    Hong, Noh, and Dudhia (2006), Mon. Wea. Rev.
    !
    !----------------------------------------------------------------------
    !
    integer , intent(in) :: nbl
    real(rkx) , intent(in) :: dt , rcl
    real(rkx) , intent(in) :: dx , dy
    integer , dimension(nbl) , intent(out) :: kpbl1d
    real(rkx) , dimension(nbl,kz) , intent(in) :: dz8w2d , pi2d
    real(rkx) , dimension(nbl,kz) , intent(in) :: ux , vx
    real(rkx) , dimension(nbl,kz) , intent(in) :: tx
    real(rkx) , dimension(nbl,kz*ndiff) , intent(in) :: qx
    real(rkx) , dimension(nbl,kzp1) , intent(in) :: p2di
    real(rkx) , dimension(nbl,kz) , intent(in) :: p2d
    real(rkx) , dimension(nbl,kz) , intent(out) :: utnp , vtnp , ttnp
    real(rkx) , dimension(nbl,kz) , intent(inout) :: tke
    real(rkx) , dimension(nbl,kz*ndiff) , intent(out) :: qtnp
    integer , dimension(nbl) , intent(in) :: xland
    real(rkx) , dimension(nbl) , intent(in) :: hfx , qfx
    real(rkx) , dimension(nbl) , intent(in) :: br , psim , psih , psfcpa
    real(rkx) , dimension(nbl) , intent(in) :: corf
    real(rkx) , dimension(nbl) , intent(in) :: wspd10
    real(rkx) , dimension(nbl) , intent(in) :: ust , znt
    real(rkx) , dimension(nbl) , intent(in) :: wspd
    real(rkx) , dimension(nbl) , intent(out) :: hpbl , wstar , delta
    integer :: i , k , ic , is , nwmass
    integer :: klpbl , kqc , kqi
    integer :: lmh
    real(rkx) :: dt2 , rdt , spdk2 , fm , fh , hol1 , gamfac , vpert
    real(rkx) :: prnum , prnum0 , ss , ri , qmean , tmean , alpha
    real(rkx) :: chi , zk , rl2 , dk , sri , brint , dtodsd , dtodsu
    real(rkx) :: rdz , dsdzt , dsdzq , dsdz2 , rlamdz
    real(rkx) :: utend , vtend , ttend , qtend
    real(rkx) :: dtstep , govrthv
    real(rkx) :: cont , conq , conw , conwrc
    real(rkx) :: delxy , pu1 , pth1 , pq1
    real(rkx) :: zfacdx , dex , hgame_c
    real(rkx) :: amf1 , amf2 , bmf2 , amf3 , bmf3
    real(rkx) :: mlfrac , ezfrac , sfcfracn , sflux0 , snlflux0
    real(rkx) :: uwst , uwstx , csfac
    real(rkx) :: prnumfac , bfx0 , hfx0 , qfx0 , delb , dux , dvx ,    &
      dsdzu , dsdzv , wm3 , dthx , dqx , ross , tem1 , dsig , tvcon ,  &
      conpr , prfac , prfac2 , phim8z , cenlfrac
    integer , dimension(nbl) :: kpbl
    real(rkx) , dimension(nbl) :: rigs , enlfrac2 , cslen , deltaoh
    real(rkx) , dimension(nbl) :: rhox , govrth , zl1 , thermal
    real(rkx) , dimension(nbl) :: wscale , hgamt , hgamq , brdn , brup
    real(rkx) , dimension(nbl) :: phim , phih , dusfc , dvsfc
    real(rkx) , dimension(nbl) :: dtsfc , dqsfc , prpbl , wspd1
    real(rkx) , dimension(nbl) :: ust3 , wstar3 , hgamu , hgamv
    real(rkx) , dimension(nbl) :: wm2 , we , bfxpbl , hfxpbl , qfxpbl
    real(rkx) , dimension(nbl) :: ufxpbl , vfxpbl , dthvx
    real(rkx) , dimension(nbl) :: brcr , sflux , zol1 , brcr_sbro
    real(rkx) , dimension(nbl) :: efxpbl , hpbl_cbl , epshol , ct
    real(rkx) , dimension(nbl,kz) :: xkzm , xkzh , f1 , f2 , r1 , r2
    real(rkx) , dimension(nbl,kz) :: ad , au , cu , al , xkzq , zfac
    real(rkx) , dimension(nbl,kz) :: thx , thvx , del , dza , dzq
    real(rkx) , dimension(nbl,kz) :: xkzom , xkzoh , za
    real(rkx) , dimension(nbl,kz) :: wscalek
    real(rkx) , dimension(nbl,kz) :: xkzml , xkzhl , zfacent , entfac
    real(rkx) , dimension(nbl,kz) :: mf , zfacmf , entfacmf
    real(rkx) , dimension(nbl,kz) :: q2x , hgame2d
    real(rkx) , dimension(nbl,kz) :: tflux_e , qflux_e , tvflux_e
    real(rkx) , dimension(nbl,kzp1) :: zq
    real(rkx) , dimension(nbl,kz,ndiff) :: r3 , f3
    real(rkx) , dimension(kz) :: uxk , vxk , txk , thxk , thvxk
    real(rkx) , dimension(kz) :: q2xk , hgame
    real(rkx) , dimension(kz) :: ps1d , pb1d , eps1d , pt1d
    real(rkx) , dimension(kz) :: xkze1d , eflx_l1d , eflx_nl1d , ptke1
    real(rkx) , dimension(2:kz) :: s2 , gh , rig , el
    real(rkx) , dimension(2:kz) :: akmk , akhk , mfk
    real(rkx) , dimension(2:kz) :: ufxpblk , vfxpblk , qfxpblk
    real(rkx) , dimension(kzp1) :: zqk
    real(rkx) , dimension(kz*ndiff) :: qxk
    logical , dimension(nbl) :: pblflg , sfcflg , stable
    logical , dimension(ndiff) :: ifvmix

    real(rkx) , parameter :: xkzminm = 0.1_rkx
    real(rkx) , parameter :: xkzminh = 0.01_rkx
    real(rkx) , parameter :: xkzmax = 1000.0_rkx
    real(rkx) , parameter :: rimin = -100._rkx
    real(rkx) , parameter :: rlam = 30.0_rkx
    real(rkx) , parameter :: prmin = 0.25_rkx
    real(rkx) , parameter :: prmax = 4.0_rkx
    real(rkx) , parameter :: brcr_ub = 0.0_rkx
    real(rkx) , parameter :: brcr_sb = 0.25_rkx
    real(rkx) , parameter :: afac = 6.8_rkx
    real(rkx) , parameter :: bfac = 6.8_rkx
    real(rkx) , parameter :: pfac = 2.0_rkx
    real(rkx) , parameter :: pfac_q = 2.0_rkx
    real(rkx) , parameter :: phifac = 8.0_rkx
    real(rkx) , parameter :: sfcfrac = 0.1_rkx
    real(rkx) , parameter :: d1 = 0.02_rkx
    real(rkx) , parameter :: d2 = 0.05_rkx
    real(rkx) , parameter :: d3 = 0.001_rkx
    real(rkx) , parameter :: h1 = 0.33333333_rkx
    real(rkx) , parameter :: h2 = 0.6666667_rkx
    real(rkx) , parameter :: zfmin = 1.e-8_rkx
    real(rkx) , parameter :: aphi5 = 5.0_rkx
    real(rkx) , parameter :: aphi16 = 16.0_rkx
    real(rkx) , parameter :: tmin = 1.e-2_rkx
    real(rkx) , parameter :: gamcrt = 3.0_rkx
    real(rkx) , parameter :: gamcrq = 2.e-3_rkx
    integer , parameter :: imvdif = 1
    !
    ! tunable parameters for tke
    !
    real(rkx) , parameter :: epsq2l = 0.01_rkx
    real(rkx) , parameter :: c_1 = 1.0_rkx
    real(rkx) , parameter :: gamcre = 0.224_rkx
    !
    ! tunable parameters for prescribed nonlocal transport profile
    !
    real(rkx) , parameter :: mltop = 1.0_rkx
    real(rkx) , parameter :: sfcfracn1 = 0.075_rkx
    real(rkx) , parameter :: nlfrac = 0.7_rkx
    real(rkx) , parameter :: enlfrac = -0.4_rkx
    real(rkx) , parameter :: a11 = 1.0_rkx
    real(rkx) , parameter :: a12 = -1.15_rkx
    real(rkx) , parameter :: ezfac = 1.5_rkx
    real(rkx) , parameter :: cpent = -0.4_rkx
    real(rkx) , parameter :: rigsmax = 100.0_rkx
    real(rkx) , parameter :: entfmin = 1.0_rkx
    real(rkx) , parameter :: entfmax = 5.0_rkx

    klpbl = kz
    lmh = 1

    cont = cpd*regrav
    conq = wlhv*regrav
    conw = regrav
    conwrc = conw*sqrt(rcl)
    conpr = bfac*vonkar*sfcfrac
    !
    !  k-start index for cloud and rain
    !
    kqc = 1 + kz
    kqi = 1 + kz*2
    nwmass = 3
    ifvmix(:) = .true.

    do k = 1 , kz
      do i = 1 , nbl
        thx(i,k) = tx(i,k)/pi2d(i,k)
      end do
    end do

    do k = 1 , kz
      do i = 1 , nbl
        tvcon = (1.0_rkx + ep1*qx(i,k))
        thvx(i,k) = thx(i,k)*tvcon
      end do
    end do

    do i = 1 , nbl
      tvcon = (1.0_rkx + ep1*qx(i,1))
      rhox(i) = psfcpa(i)/(rdry*tx(i,1)*tvcon)
      govrth(i) = egrav/thx(i,1)
    end do
    !
    !-----compute the height of full- and half-sigma levels above ground
    !     level, and the layer thicknesses.
    !
    do i = 1 , nbl
      zq(i,1) = 0.0_rkx
    end do

    do k = 1 , kz
      do i = 1 , nbl
        zq(i,k+1) = dz8w2d(i,k)+zq(i,k)
      end do
    end do

    do k = 1 , kz
      do i = 1 , nbl
        za(i,k) = 0.5_rkx*(zq(i,k)+zq(i,k+1))
        dzq(i,k) = zq(i,k+1)-zq(i,k)
        del(i,k) = p2di(i,k)-p2di(i,k+1)
      end do
    end do

    do i = 1 , nbl
      dza(i,1) = za(i,1)
    end do

    do k = 2 , kz
      do i = 1 , nbl
        dza(i,k) = za(i,k)-za(i,k-1)
      end do
    end do
    !
    ! initialize vertical tendencies
    !

    do k = 1 , kz
      do i = 1 , nbl
        utnp(i,k) = 0.0_rkx
        vtnp(i,k) = 0.0_rkx
        ttnp(i,k) = 0.0_rkx
      end do
    end do

    do k = 1 , kz*ndiff
      do i = 1 , nbl
        qtnp(i,k) = 0.0_rkx
      end do
    end do

    do i = 1 , nbl
      wspd1(i) = sqrt(ux(i,1)*ux(i,1)+vx(i,1)*vx(i,1))+1.e-9_rkx
    end do
    !
    !---- compute vertical diffusion
    !
    ! - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
    !     compute preliminary variables
    !
    dtstep = dt
    dt2 = 2.0_rkx*dtstep
    rdt = 1.0_rkx/dt2

    do i = 1 , nbl
      bfxpbl(i) = 0.0_rkx
      hfxpbl(i) = 0.0_rkx
      qfxpbl(i) = 0.0_rkx
      ufxpbl(i) = 0.0_rkx
      vfxpbl(i) = 0.0_rkx
      hgamu(i)  = 0.0_rkx
      hgamv(i)  = 0.0_rkx
      delta(i)  = 0.0_rkx
    end do

    do i = 1 , nbl
      efxpbl(i)   = 0.0_rkx
      hpbl_cbl(i) = 0.0_rkx
      epshol(i)   = 0.0_rkx
      ct(i)       = 0.0_rkx
    end do

    do i = 1 , nbl
      deltaoh(i)  = 0.0_rkx
      rigs(i)     = 0.0_rkx
      enlfrac2(i) = 0.0_rkx
      cslen(i)    = 0.0_rkx
    end do

    do k = 1 , klpbl
      do i = 1 , nbl
        wscalek(i,k) = 0.0_rkx
      end do
    end do

    do k = 1 , klpbl
      do i = 1 , nbl
        zfac(i,k) = 0.0_rkx
      end do
    end do

    do k = 1 , kz
      do i = 1 , nbl
        q2x(i,k) = 2.0_rkx*tke(i,k)
      end do
    end do

    do k = 1 , kz
      do i = 1 , nbl
        hgame2d(i,k)  = 0.0_rkx
        tflux_e(i,k)  = 0.0_rkx
        qflux_e(i,k)  = 0.0_rkx
        tvflux_e(i,k) = 0.0_rkx
      end do
    end do

    do k = 1 , kz
      do i = 1 , nbl
        mf(i,k)     = 0.0_rkx
        zfacmf(i,k) = 0.0_rkx
      end do
    end do

    do k = 1 , klpbl-1
      do i = 1 , nbl
        xkzom(i,k) = xkzminm
        xkzoh(i,k) = xkzminh
      end do
    end do

    do i = 1 , nbl
      dusfc(i) = 0.0_rkx
      dvsfc(i) = 0.0_rkx
      dtsfc(i) = 0.0_rkx
      dqsfc(i) = 0.0_rkx
    end do

    do i = 1 , nbl
      hgamt(i)  = 0.0_rkx
      hgamq(i)  = 0.0_rkx
      wscale(i) = 0.0_rkx
      kpbl(i)   = 1
      hpbl(i)   = zq(i,1)
      hpbl_cbl(i) = zq(i,1)
      zl1(i)    = za(i,1)
      thermal(i)= thvx(i,1)
      pblflg(i) = .true.
      sfcflg(i) = .true.
      sflux(i) = hfx(i)/rhox(i)*rcpd + qfx(i)/rhox(i)*ep1*thx(i,1)
      if ( br(i) > 0.0 ) sfcflg(i) = .false.
    end do
    !
    ! compute the first guess of pbl height
    !
    do i = 1 , nbl
      stable(i) = .false.
      brup(i) = br(i)
      brcr(i) = brcr_ub
    end do

    do k = 2 , klpbl
      do i = 1 , nbl
        if ( .not.stable(i) ) then
          brdn(i) = brup(i)
          spdk2   = max(ux(i,k)**2+vx(i,k)**2,1.0_rkx)
          brup(i) = (thvx(i,k)-thermal(i))*(egrav*za(i,k)/thvx(i,1))/spdk2
          kpbl(i) = k
          stable(i) = brup(i) > brcr(i)
        end if
      end do
    end do

    do i = 1 , nbl
      k = kpbl(i)
      if ( brdn(i) >= brcr(i) ) then
        brint = 0.0_rkx
      else if ( brup(i) <= brcr(i) ) then
        brint = 1.0_rkx
      else
        brint = (brcr(i)-brdn(i))/(brup(i)-brdn(i))
      end if
      hpbl(i) = za(i,k-1)+brint*(za(i,k)-za(i,k-1))
      if ( hpbl(i) < zq(i,2) ) kpbl(i) = 1
      if ( kpbl(i) <= 1 ) pblflg(i) = .false.
    end do

    do i = 1 , nbl
      fm = psim(i)
      fh = psih(i)
      zol1(i) = max(br(i)*fm*fm/fh,rimin)
      if ( sfcflg(i) ) then
        zol1(i) = min(zol1(i),-zfmin)
      else
        zol1(i) = max(zol1(i),zfmin)
      end if
      hol1 = zol1(i)*hpbl(i)/zl1(i)*sfcfrac
      epshol(i) = hol1
      if ( sfcflg(i) ) then
        phim(i) = (1.0_rkx-aphi16*hol1)**(-0.25_rkx)
        phih(i) = (1.0_rkx-aphi16*hol1)**(-0.5_rkx)
        bfx0  = max(sflux(i),0.0_rkx)
        hfx0 = max(hfx(i)/rhox(i)*rcpd,0.0_rkx)
        qfx0 = max(ep1*thx(i,1)*qfx(i)/rhox(i),0.0_rkx)
        wstar3(i) = (govrth(i)*bfx0*hpbl(i))
        wstar(i) = (wstar3(i))**h1
      else
        phim(i) = (1.0_rkx+aphi5*hol1)
        phih(i) = phim(i)
        wstar(i)  = 0.0_rkx
        wstar3(i) = 0.0_rkx
      end if
      ust3(i)   = ust(i)**3
      wscale(i) = (ust3(i)+phifac*vonkar*wstar3(i)*0.5_rkx)**h1
      wscale(i) = min(wscale(i),ust(i)*aphi16)
      wscale(i) = max(wscale(i),ust(i)/aphi5)
    end do
    !
    ! compute the surface variables for pbl height estimation
    ! under unstable conditions
    !
    do i = 1 , nbl
      if ( sfcflg(i) .and. sflux(i) > 0.0 ) then
        gamfac   = bfac/rhox(i)/wscale(i)
        hgamt(i) = min(gamfac*hfx(i)*rcpd,gamcrt)
        hgamq(i) = min(gamfac*qfx(i),gamcrq)
        vpert = (hgamt(i)+ep1*thx(i,1)*hgamq(i))/bfac*afac
        thermal(i) = thermal(i)+max(vpert,0.0_rkx) * &
          min(za(i,1)/(sfcfrac*hpbl(i)),1.0_rkx)
        hgamt(i) = max(hgamt(i),0.0_rkx)
        hgamq(i) = max(hgamq(i),0.0_rkx)
        brint    = -15.9_rkx*ust(i)*ust(i)/wspd(i)*wstar3(i)/(wscale(i)**4)
        hgamu(i) = brint*ux(i,1)
        hgamv(i) = brint*vx(i,1)
      else
        pblflg(i) = .false.
      end if
    end do
    !
    ! enhance the pbl height by considering the thermal
    !
    do i = 1 , nbl
      if ( pblflg(i) ) then
        kpbl(i) = 1
        hpbl(i) = zq(i,1)
      end if
    end do

    do i = 1 , nbl
      if ( pblflg(i) ) then
        stable(i) = .false.
        brup(i) = br(i)
        brcr(i) = brcr_ub
      end if
    end do

    do k = 2 , klpbl
      do i = 1 , nbl
        if ( .not. stable(i) .and. pblflg(i) ) then
          brdn(i) = brup(i)
          spdk2   = max(ux(i,k)**2+vx(i,k)**2,1.0_rkx)
          brup(i) = (thvx(i,k)-thermal(i))*(egrav*za(i,k)/thvx(i,1))/spdk2
          kpbl(i) = k
          stable(i) = brup(i) > brcr(i)
        end if
      end do
    end do

    do i = 1 , nbl
      if ( pblflg(i) ) then
        k = kpbl(i)
        if ( brdn(i) >= brcr(i) ) then
          brint = 0.0_rkx
        else if ( brup(i) <= brcr(i) ) then
          brint = 1.0_rkx
        else
          brint = (brcr(i)-brdn(i))/(brup(i)-brdn(i))
        end if
        hpbl(i) = za(i,k-1)+brint*(za(i,k)-za(i,k-1))
        if ( hpbl(i) < zq(i,2) ) kpbl(i) = 1
        if ( kpbl(i) < 1 ) pblflg(i) = .false.
        if ( wstar(i) /= 0 ) then
          uwst  = abs(ust(i)/wstar(i)-0.5_rkx)
          uwstx = -80.0_rkx*uwst+14.0_rkx
          csfac = 0.5_rkx*(tanh(uwstx)+3.0_rkx)
        else
           csfac = 1.0_rkx
        end if
        cslen(i) = csfac*hpbl(i)
      end if
    end do
    !
    ! stable boundary layer
    !
    do i = 1 , nbl
      hpbl_cbl(i) = hpbl(i)
      if ( (.not. sfcflg(i)) .and. hpbl(i) < zq(i,2) ) then
        brup(i) = br(i)
        stable(i) = .false.
      else
        stable(i) = .true.
      end if
    end do

    do i = 1 , nbl
      if ( ( .not. stable(i)) .and. xland(i) > 0 ) then
        ross = wspd10(i) / (max(corf(i),2.546e-5_rkx)*znt(i))
        brcr_sbro(i) = min(0.16_rkx*(1.e-7_rkx*ross)**(-0.18_rkx),0.3_rkx)
      end if
    end do

    do i = 1 , nbl
      if ( .not.stable(i) ) then
        if ( xland(i) > 0 ) then
          brcr(i) = brcr_sbro(i)
        else
          brcr(i) = brcr_sb
        end if
      end if
    end do

    do k = 2 , klpbl
      do i = 1 , nbl
        if ( .not. stable(i) ) then
          brdn(i) = brup(i)
          spdk2   = max(ux(i,k)**2+vx(i,k)**2,1.0_rkx)
          brup(i) = (thvx(i,k)-thermal(i))*(egrav*za(i,k)/thvx(i,1))/spdk2
          kpbl(i) = k
          stable(i) = brup(i) > brcr(i)
        end if
      end do
    end do

    do i = 1 , nbl
      if ( (.not.sfcflg(i)) .and. hpbl(i) < zq(i,2) ) then
        k = kpbl(i)
        if ( brdn(i) >= brcr(i) ) then
          brint = 0.0_rkx
        else if ( brup(i) <= brcr(i) ) then
          brint = 1.0_rkx
        else
          brint = (brcr(i)-brdn(i))/(brup(i)-brdn(i))
        end if
        hpbl(i) = za(i,k-1)+brint*(za(i,k)-za(i,k-1))
        if ( hpbl(i) < zq(i,2) ) kpbl(i) = 1
        if ( kpbl(i) < 1 ) pblflg(i) = .false.
      end if
    end do
    !
    ! scale dependency for nonlocal momentum and moisture transport
    !
    delxy=sqrt(dx*dy)
    do i = 1 , nbl
      pu1 = pu(delxy,cslen(i))
      pq1 = pq(delxy,cslen(i))
      if ( pblflg(i) ) then
        hgamu(i) = hgamu(i)*pu1
        hgamv(i) = hgamv(i)*pu1
        hgamq(i) = hgamq(i)*pq1
      end if
    end do
    !
    ! estimate the entrainment parameters
    !
    delxy = sqrt(dx*dy)
    do i = 1 , nbl
      if ( pblflg(i) ) then
        k = kpbl(i) - 1
        prpbl(i) = 1.0_rkx
        wm3       = wstar3(i) + 5.0_rkx * ust3(i)
        wm2(i)    = wm3**h2
        bfxpbl(i) = -0.15_rkx*thvx(i,1)*regrav*wm3/hpbl(i)
        dthvx(i)  = max(thvx(i,k+1)-thvx(i,k),tmin)
        dthx  = max(thx(i,k+1)-thx(i,k),tmin)
        dqx   = min(qx(i,k+1)-qx(i,k),0.0_rkx)
        we(i) = max(bfxpbl(i)/dthvx(i),-sqrt(wm2(i)))
        hfxpbl(i) = we(i)*dthx
        pq1 = pq(delxy,cslen(i))
        qfxpbl(i) = we(i)*dqx*pq1
        pu1 = pu(delxy,cslen(i))
        dux = ux(i,k+1)-ux(i,k)
        dvx = vx(i,k+1)-vx(i,k)
        if ( dux > tmin ) then
          ufxpbl(i) = max(prpbl(i)*we(i)*dux*pu1,-ust(i)*ust(i))
        else if ( dux < -tmin ) then
          ufxpbl(i) = min(prpbl(i)*we(i)*dux*pu1,ust(i)*ust(i))
        else
          ufxpbl(i) = 0.0_rkx
        end if
        if ( dvx > tmin ) then
          vfxpbl(i) = max(prpbl(i)*we(i)*dvx*pu1,-ust(i)*ust(i))
        else if ( dvx < -tmin ) then
          vfxpbl(i) = min(prpbl(i)*we(i)*dvx*pu1,ust(i)*ust(i))
        else
          vfxpbl(i) = 0.0_rkx
        end if
        delb  = govrth(i)*d3*hpbl(i)
        delta(i) = min(d1*hpbl(i) + d2*wm2(i)/delb,100.0_rkx)
        delb  = govrth(i)*dthvx(i)
        deltaoh(i) = d1*hpbl(i) + d2*wm2(i)/delb
        deltaoh(i) = max(ezfac*deltaoh(i),hpbl(i)-za(i,kpbl(i)-1)-1.0_rkx)
        deltaoh(i) = min(deltaoh(i)      ,hpbl(i))
        if ( (dux /= 0) .or. (dvx /= 0) ) then
          rigs(i) = govrth(i)*dthvx(i)*deltaoh(i)/(dux**2+dvx**2)
        else
          rigs(i) = rigsmax
        end if
        rigs(i)     = max(min(rigs(i), rigsmax),rimin)
        if ( (rigs(i) > 0) .and. (abs(rigs(i)+cpent) <= 1.e-6_rkx) ) then
          cenlfrac = entfmax
        else
          cenlfrac = rigs(i)/(rigs(i)+cpent)
        end if
        cenlfrac = min(cenlfrac,entfmax)
        enlfrac2(i) = max(wm3/wstar3(i)*cenlfrac, entfmin)
        enlfrac2(i) = enlfrac2(i)*enlfrac
      end if
    end do

    do k = 1 , klpbl
      do i = 1 , nbl
        if ( pblflg(i) ) then
          entfacmf(i,k) = sqrt(((zq(i,k+1)-hpbl(i))/deltaoh(i))**2)
        end if
        if ( pblflg(i) .and. k >= kpbl(i) ) then
          entfac(i,k) = ((zq(i,k+1)-hpbl(i))/deltaoh(i))**2
        else
          entfac(i,k) = 1.e30_rkx
        end if
      end do
    end do
    !
    ! compute diffusion coefficients below pbl
    !
    do k = 1 , klpbl
      do i = 1 , nbl
        if ( k < kpbl(i) ) then
          zfac(i,k) = min(max((1.0_rkx-(zq(i,k+1)-zl1(i)) / &
            (hpbl(i)-zl1(i))),zfmin),1.0_rkx)
          zfacent(i,k) = (1.0_rkx-zfac(i,k))**3
          wscalek(i,k) = (ust3(i)+phifac*vonkar*wstar3(i) * &
            (1.0_rkx-zfac(i,k)))**h1
          if ( sfcflg(i) ) then
            prfac = conpr
            prfac2 = 15.9_rkx*wstar3(i)/ust3(i) / &
              (1.0_rkx+4.0_rkx*vonkar*wstar3(i)/ust3(i))
            prnumfac = -3.0_rkx * &
              (max(zq(i,k+1)-sfcfrac*hpbl(i),0.0_rkx))**2/hpbl(i)**2
          else
            prfac = 0.0_rkx
            prfac2 = 0.0_rkx
            prnumfac = 0.0_rkx
            phim8z = 1.0_rkx+aphi5*zol1(i)*zq(i,k+1)/zl1(i)
            wscalek(i,k) = ust(i)/phim8z
            wscalek(i,k) = max(wscalek(i,k),0.001_rkx)
          end if
          prnum0 = (phih(i)/phim(i)+prfac)
          prnum0 = max(min(prnum0,prmax),prmin)
          xkzm(i,k) = wscalek(i,k)*vonkar*zq(i,k+1)*zfac(i,k)**pfac
          prnum =  1.0_rkx + (prnum0-1.0_rkx)*exp(prnumfac)
          xkzq(i,k) = xkzm(i,k)/prnum*zfac(i,k)**(pfac_q-pfac)
          prnum0 = prnum0/(1.0_rkx+prfac2*vonkar*sfcfrac)
          prnum =  1.0_rkx + (prnum0-1.0_rkx)*exp(prnumfac)
          xkzh(i,k) = xkzm(i,k)/prnum
          xkzm(i,k) = xkzm(i,k)+xkzom(i,k)
          xkzh(i,k) = xkzh(i,k)+xkzoh(i,k)
          xkzq(i,k) = xkzq(i,k)+xkzoh(i,k)
          xkzm(i,k) = min(xkzm(i,k),xkzmax)
          xkzh(i,k) = min(xkzh(i,k),xkzmax)
          xkzq(i,k) = min(xkzq(i,k),xkzmax)
        end if
      end do
    end do
    !
    ! compute diffusion coefficients over pbl (free atmosphere)
    !
    do k = 1 , kzm1
      do i = 1 , nbl
        if ( k >= kpbl(i) ) then
          ss = ((ux(i,k+1)-ux(i,k))*(ux(i,k+1)-ux(i,k)) + &
                (vx(i,k+1)-vx(i,k))*(vx(i,k+1)-vx(i,k))) / &
               (dza(i,k+1)*dza(i,k+1))+1.e-9_rkx
          govrthv = egrav/(0.5_rkx*(thvx(i,k+1)+thvx(i,k)))
          ri = govrthv*(thvx(i,k+1)-thvx(i,k))/(ss*dza(i,k+1))
          ! in cloud
          if ( imvdif == 1 .and. nwmass >= 3 ) then
            if ( (qx(i,kqc+k-1)+qx(i,kqi+k-1)) > 0.01e-3_rkx .and. &
                 (qx(i,kqc+k)+qx(i,kqi+k))     > 0.01e-3_rkx ) then
              qmean = 0.5_rkx*(qx(i,k)+qx(i,k+1))
              tmean = 0.5_rkx*(tx(i,k)+tx(i,k+1))
              alpha = wlhv*qmean/rdry/tmean
              chi   = wlhv*wlhv*qmean*rcpd/rwat/tmean/tmean
              ri    = (1.0_rkx+alpha)*(ri-egrav*egrav/ss/ &
                tmean*rcpd*((chi-alpha)/(1.0_rkx+chi)))
            end if
          end if
          zk = vonkar*zq(i,k+1)
          rlamdz = min(max(0.1_rkx*dza(i,k+1),rlam),300.0_rkx)
          rlamdz = min(dza(i,k+1),rlamdz)
          rl2 = (zk*rlamdz/(rlamdz+zk))**2
          dk = rl2*sqrt(ss)
          if ( ri < 0.0_rkx ) then
            ! unstable regime
            ri = max(ri, rimin)
            sri = sqrt(-ri)
            xkzm(i,k) = dk*(1.0_rkx + 8.0_rkx*(-ri)/(1.0_rkx+1.746_rkx*sri))
            xkzh(i,k) = dk*(1.0_rkx + 8.0_rkx*(-ri)/(1.0_rkx+1.286_rkx*sri))
          else
            ! stable regime
            xkzh(i,k) = dk/(1.0_rkx + 5.0_rkx*ri)**2
            prnum = 1.0_rkx + 2.1_rkx*ri
            prnum = min(prnum,prmax)
            xkzm(i,k) = xkzh(i,k)*prnum
          end if

          xkzm(i,k) = xkzm(i,k)+xkzom(i,k)
          xkzh(i,k) = xkzh(i,k)+xkzoh(i,k)
          xkzm(i,k) = min(xkzm(i,k),xkzmax)
          xkzh(i,k) = min(xkzh(i,k),xkzmax)
          xkzml(i,k) = xkzm(i,k)
          xkzhl(i,k) = xkzh(i,k)
        end if
      end do
    end do
    !
    ! prescribe nonlocal heat transport below pbl
    !
    do i = 1 , nbl
      deltaoh(i) = deltaoh(i)/hpbl(i)
    end do

    delxy = sqrt(dx*dy)
    do i = 1 , nbl
      mlfrac      = mltop-deltaoh(i)
      ezfrac      = mltop+deltaoh(i)
      zfacmf(i,1) = min(max((zq(i,2)/hpbl(i)),zfmin),1.0_rkx)
      sfcfracn    = max(sfcfracn1,zfacmf(i,1))
      sflux0      = (a11+a12*sfcfracn)*sflux(i)
      snlflux0    = nlfrac*sflux0
      amf1        = snlflux0/sfcfracn
      if ( pblflg(i) ) then
        amf2        = -snlflux0/(mlfrac-sfcfracn)
        bmf2        = -mlfrac*amf2
      end if
      if ( (deltaoh(i) == 0) .and. (enlfrac2(i) == 0)) then
         amf3       = 0.0_rkx
      else
         amf3       = snlflux0*enlfrac2(i)/deltaoh(i)
      end if
      bmf3        = -amf3*mlfrac
      hfxpbl(i)   = amf3+bmf3
      pth1        = pthnl(delxy,cslen(i))
      hfxpbl(i)   = hfxpbl(i)*pth1

      do k = 1 , klpbl
        zfacmf(i,k) = max((zq(i,k+1)/hpbl(i)),zfmin)
        if ( pblflg(i) .and. k < kpbl(i) ) then
          if ( zfacmf(i,k) <= sfcfracn ) then
            mf(i,k) = amf1*zfacmf(i,k)
          else if ( zfacmf(i,k) <= mlfrac ) then
            mf(i,k) = amf2*zfacmf(i,k)+bmf2
          end if
          mf(i,k) = mf(i,k)+hfxpbl(i)*exp(-entfacmf(i,k))
          mf(i,k) = mf(i,k)*pth1
        end if
      end do
    end do
    !
    ! compute tridiagonal matrix elements for heat
    !
    do k = 1 , kz
      do i = 1 , nbl
        au(i,k) = 0.0_rkx
        al(i,k) = 0.0_rkx
        ad(i,k) = 0.0_rkx
        f1(i,k) = 0.0_rkx
      end do
    end do

    do i = 1 , nbl
      ad(i,1) = 1.0_rkx
      f1(i,1) = thx(i,1)-300.0_rkx+hfx(i)/cont/del(i,1)*dt2
    end do

    delxy = sqrt(dx*dy)
    do k = 1 , kzm1
      do i = 1 , nbl
        dtodsd = dt2/del(i,k)
        dtodsu = dt2/del(i,k+1)
        dsig   = p2d(i,k)-p2d(i,k+1)
        rdz    = 1.0_rkx/dza(i,k+1)
        tem1   = dsig*xkzh(i,k)*rdz
        if ( pblflg(i) .and. k < kpbl(i) ) then
          dsdzt = tem1*(-mf(i,k)/xkzh(i,k))
          f1(i,k)   = f1(i,k)+dtodsd*dsdzt
          f1(i,k+1) = thx(i,k+1)-300.0_rkx-dtodsu*dsdzt
        else if ( pblflg(i) .and. &
                  k >= kpbl(i) .and. entfac(i,k) < 4.6_rkx ) then
          xkzh(i,k) = -we(i)*dza(i,kpbl(i))*exp(-entfac(i,k))
          xkzh(i,k) = sqrt(xkzh(i,k)*xkzhl(i,k))
          xkzh(i,k) = max(xkzh(i,k),xkzoh(i,k))
          xkzh(i,k) = min(xkzh(i,k),xkzmax)
          f1(i,k+1) = thx(i,k+1)-300.0_rkx
        else
          f1(i,k+1) = thx(i,k+1)-300.0_rkx
        end if
        tem1   = dsig*xkzh(i,k)*rdz
        dsdz2     = tem1*rdz
        au(i,k)   = -dtodsd*dsdz2
        al(i,k)   = -dtodsu*dsdz2
        !
        ! scale dependency for local heat transport
        !
        zfacdx = 0.2_rkx*hpbl(i)/zq(i,k+1)
        delxy = sqrt(dx*dy)*max(zfacdx,1.0_rkx)
        pth1 = pthl(delxy,hpbl(i))
        if ( pblflg(i) .and. k < kpbl(i) ) then
          au(i,k) = au(i,k)*pth1
          al(i,k) = al(i,k)*pth1
        end if
        ad(i,k)   = ad(i,k)-au(i,k)
        ad(i,k+1) = 1.0_rkx-al(i,k)
      end do
    end do
    !
    ! copies here to avoid duplicate input args for tridin
    !
    do k = 1 , kz
      do i = 1 , nbl
        cu(i,k) = au(i,k)
        r1(i,k) = f1(i,k)
      end do
    end do

    call tridi1(al,ad,cu,r1,au,f1,nbl)
    !
    ! recover tendencies of heat
    !
    do k = kz , 1 , -1
      do i = 1 , nbl
        ttend = (f1(i,k)-thx(i,k)+300.0_rkx)*rdt*pi2d(i,k)
        ttnp(i,k) = ttnp(i,k)+ttend
        dtsfc(i) = dtsfc(i)+ttend*cont*del(i,k)/pi2d(i,k)
        if ( k == kz ) then
          tflux_e(i,k) = ttend*dz8w2d(i,k)
        else
          tflux_e(i,k) = tflux_e(i,k+1) + ttend*dz8w2d(i,k)
        end if
      end do
    end do
    !
    ! compute tridiagonal matrix elements for moisture, clouds, and gases
    !
    do k = 1 , kz
      do i = 1 , nbl
        au(i,k) = 0.0_rkx
        al(i,k) = 0.0_rkx
        ad(i,k) = 0.0_rkx
      end do
    end do

    do ic = 1 , ndiff
      do i = 1 , nbl
        do k = 1 , kz
          f3(i,k,ic) = 0.0_rkx
        end do
      end do
    end do

    do i = 1 , nbl
      ad(i,1) = 1.0_rkx
      f3(i,1,1) = qx(i,1)+qfx(i)*egrav/del(i,1)*dt2
    end do

    if ( ndiff >= 2 ) then
      do ic = 2 , ndiff
        is = (ic-1) * kz
        do i = 1 , nbl
          f3(i,1,ic) = qx(i,1+is)
        end do
      end do
    end if

    do k = 1  ,kzm1
      do i = 1 , nbl
        if ( k >= kpbl(i) ) then
          xkzq(i,k) = xkzh(i,k)
        end if
      end do
    end do

    do k = 1 , kzm1
      do i = 1 , nbl
        dtodsd = dt2/del(i,k)
        dtodsu = dt2/del(i,k+1)
        dsig   = p2d(i,k)-p2d(i,k+1)
        rdz    = 1.0_rkx/dza(i,k+1)
        tem1   = dsig*xkzq(i,k)*rdz
        if ( pblflg(i) .and. k < kpbl(i) ) then
          dsdzq = tem1*(-qfxpbl(i)*zfacent(i,k)/xkzq(i,k))
          f3(i,k,1) = f3(i,k,1)+dtodsd*dsdzq
          f3(i,k+1,1) = qx(i,k+1)-dtodsu*dsdzq
        else if ( pblflg(i) .and. &
                  k >= kpbl(i) .and. entfac(i,k) < 4.6_rkx ) then
          xkzq(i,k) = -we(i)*dza(i,kpbl(i))*exp(-entfac(i,k))
          xkzq(i,k) = sqrt(xkzq(i,k)*xkzhl(i,k))
          xkzq(i,k) = max(xkzq(i,k),xkzoh(i,k))
          xkzq(i,k) = min(xkzq(i,k),xkzmax)
          f3(i,k+1,1) = qx(i,k+1)
        else
          f3(i,k+1,1) = qx(i,k+1)
        end if
        tem1   = dsig*xkzq(i,k)*rdz
        dsdz2     = tem1*rdz
        au(i,k)   = -dtodsd*dsdz2
        al(i,k)   = -dtodsu*dsdz2
        !
        ! scale dependency for local moisture transport
        !
        zfacdx = 0.2_rkx*hpbl(i)/zq(i,k+1)
        delxy = sqrt(dx*dy)*max(zfacdx,1.0_rkx)
        pq1 = pq(delxy,hpbl(i))
        if ( pblflg(i) .and. k < kpbl(i) ) then
          au(i,k) = au(i,k)*pq1
          al(i,k) = al(i,k)*pq1
        end if
        ad(i,k)   = ad(i,k)-au(i,k)
        ad(i,k+1) = 1.0_rkx-al(i,k)
      end do
    end do

    if ( ndiff >= 2 ) then
      do ic = 2 , ndiff
        is = (ic-1) * kz
        do k = 1 , kzm1
          do i = 1 , nbl
            f3(i,k+1,ic) = qx(i,k+1+is)
          end do
        end do
      end do
    end if
    !
    ! copies here to avoid duplicate input args for tridin
    !
    do k = 1 , kz
      do i = 1 , nbl
        cu(i,k) = au(i,k)
      end do
    end do

    do ic = 1 , ndiff
      do k = 1 , kz
        do i = 1 , nbl
          r3(i,k,ic) = f3(i,k,ic)
        end do
      end do
    end do
    !
    ! solve tridiagonal problem for moisture, clouds, and gases
    !
    call tridin_ysu(al,ad,cu,r3,au,f3,nbl,ndiff)
    !
    ! recover tendencies of heat and moisture
    !
    do k = kz , 1 , -1
      do i = 1 , nbl
        qtend = (f3(i,k,1)-qx(i,k))*rdt
        qtnp(i,k) = qtnp(i,k)+qtend
        dqsfc(i) = dqsfc(i)+qtend*conq*del(i,k)
        if ( k == kz ) then
          qflux_e(i,k) = qtend*dz8w2d(i,k)
        else
          qflux_e(i,k) = qflux_e(i,k+1) + qtend*dz8w2d(i,k)
        end if
        tvflux_e(i,k) = tflux_e(i,k) + qflux_e(i,k)*ep1*thx(i,k)
      end do
    end do

    do k = 1 , kz
      do i = 1 , nbl
        if ( pblflg(i) .and. k < kpbl(i) ) then
          hgame_c = c_1*0.2_rkx*2.5_rkx * &
            (egrav/thvx(i,k))*wstar(i)/(0.25_rkx*(q2x(i,k+1)+q2x(i,k)))
          hgame_c = min(hgame_c,gamcre)
          if ( k == kz ) then
            hgame2d(i,k) = hgame_c*0.5_rkx*tvflux_e(i,k)*hpbl(i)
            hgame2d(i,k) = max(hgame2d(i,k),0.0_rkx)
          else
            hgame2d(i,k) = hgame_c * 0.5_rkx * &
              (tvflux_e(i,k)+tvflux_e(i,k+1))*hpbl(i)
            hgame2d(i,k) = max(hgame2d(i,k),0.0_rkx)
          end if
        end if
      end do
    end do

    if ( ndiff >= 2 ) then
      do ic = 2 , ndiff
        if ( ifvmix(ic) ) then
          is = (ic-1) * kz
          do k = kz , 1 , -1
            do i = 1 , nbl
              qtend = (f3(i,k,ic)-qx(i,k+is))*rdt
              qtnp(i,k+is) = qtnp(i,k+is)+qtend
            end do
          end do
        end if
      end do
    end if
    !
    ! compute tridiagonal matrix elements for momentum
    !
    do i = 1 , nbl
      do k = 1 , kz
        au(i,k) = 0.0_rkx
        al(i,k) = 0.0_rkx
        ad(i,k) = 0.0_rkx
        f1(i,k) = 0.0_rkx
        f2(i,k) = 0.0_rkx
      end do
    end do

    do i = 1 , nbl
      ad(i,1) = 1.0_rkx+ust(i)**2/wspd1(i) * &
                rhox(i)*egrav/del(i,1)*dt2*(wspd1(i)/wspd(i))**2
      f1(i,1) = ux(i,1)
      f2(i,1) = vx(i,1)
    end do

    delxy = sqrt(dx*dy)
    do k = 1 , kzm1
      do i = 1 , nbl
        dtodsd = dt2/del(i,k)
        dtodsu = dt2/del(i,k+1)
        dsig   = p2d(i,k)-p2d(i,k+1)
        rdz    = 1.0_rkx/dza(i,k+1)
        tem1   = dsig*xkzm(i,k)*rdz
        if ( pblflg(i) .and. k < kpbl(i) ) then
          dsdzu     = tem1*(-hgamu(i)/hpbl(i) - &
                      ufxpbl(i)*zfacent(i,k)/xkzm(i,k))
          dsdzv     = tem1*(-hgamv(i)/hpbl(i) - &
                      vfxpbl(i)*zfacent(i,k)/xkzm(i,k))
          f1(i,k)   = f1(i,k)+dtodsd*dsdzu
          f1(i,k+1) = ux(i,k+1)-dtodsu*dsdzu
          f2(i,k)   = f2(i,k)+dtodsd*dsdzv
          f2(i,k+1) = vx(i,k+1)-dtodsu*dsdzv
        else if ( pblflg(i) .and. &
                  k >= kpbl(i) .and. entfac(i,k) < 4.6_rkx ) then
          xkzm(i,k) = prpbl(i)*xkzh(i,k)
          xkzm(i,k) = sqrt(xkzm(i,k)*xkzml(i,k))
          xkzm(i,k) = max(xkzm(i,k),xkzom(i,k))
          xkzm(i,k) = min(xkzm(i,k),xkzmax)
          f1(i,k+1) = ux(i,k+1)
          f2(i,k+1) = vx(i,k+1)
        else
          f1(i,k+1) = ux(i,k+1)
          f2(i,k+1) = vx(i,k+1)
        end if
        tem1   = dsig*xkzm(i,k)*rdz
        dsdz2     = tem1*rdz
        au(i,k)   = -dtodsd*dsdz2
        al(i,k)   = -dtodsu*dsdz2
        !
        ! scale dependency for local momentum transport
        !
        zfacdx = 0.2_rkx*hpbl(i)/zq(i,k+1)
        delxy = sqrt(dx*dy)*max(zfacdx,1.0_rkx)
        pu1 = pu(delxy,hpbl(i))
        if ( pblflg(i) .and. k < kpbl(i) ) then
          au(i,k) = au(i,k)*pu1
          al(i,k) = al(i,k)*pu1
        end if
        ad(i,k)   = ad(i,k)-au(i,k)
        ad(i,k+1) = 1.0_rkx-al(i,k)
      end do
    end do
    !
    ! copies here to avoid duplicate input args for tridin
    !
    do k = 1 , kz
      do i = 1 , nbl
        cu(i,k) = au(i,k)
        r1(i,k) = f1(i,k)
        r2(i,k) = f2(i,k)
      end do
    end do
    !
    ! solve tridiagonal problem for momentum
    !
    call tridi2(al,ad,cu,r1,r2,au,f1,f2,nbl)
    !
    ! recover tendencies of momentum
    !
    do k = kz , 1 , -1
      do i = 1 , nbl
        utend = (f1(i,k)-ux(i,k))*rdt
        vtend = (f2(i,k)-vx(i,k))*rdt
        utnp(i,k) = utnp(i,k)+utend
        vtnp(i,k) = vtnp(i,k)+vtend
        dusfc(i) = dusfc(i) + utend*conwrc*del(i,k)
        dvsfc(i) = dvsfc(i) + vtend*conwrc*del(i,k)
      end do
    end do
    !
    do i = 1 , nbl
      kpbl1d(i) = kpbl(i)
    end do
    !
    !---- calculate sgs tke which is consistent with shinhongpbl algorithm
    !
    tke_calculation: &
    do i = 1 , nbl
      do k = 2 , kz
        s2(k)   = 0.0_rkx
        gh(k)   = 0.0_rkx
        rig(k)  = 0.0_rkx
        el(k)   = 0.0_rkx
        akmk(k) = 0.0_rkx
        akhk(k) = 0.0_rkx
        mfk(k)      = 0.0_rkx
        ufxpblk(k)  = 0.0_rkx
        vfxpblk(k)  = 0.0_rkx
        qfxpblk(k)  = 0.0_rkx
      end do
      do k = 1 , kz
        uxk(k)   = 0.0_rkx
        vxk(k)   = 0.0_rkx
        txk(k)   = 0.0_rkx
        thxk(k)  = 0.0_rkx
        thvxk(k) = 0.0_rkx
        q2xk(k)  = 0.0_rkx
        hgame(k) = 0.0_rkx
        ps1d(k)  = 0.0_rkx
        pb1d(k)  = 0.0_rkx
        eps1d(k) = 0.0_rkx
        pt1d(k)  = 0.0_rkx
        xkze1d(k)    = 0.0_rkx
        eflx_l1d(k)  = 0.0_rkx
        eflx_nl1d(k) = 0.0_rkx
        ptke1(k)     = 1.0_rkx
      end do
      do k = 1 , kzp1
        zqk(k)   = 0.0_rkx
      end do
      do k = 1 , kz*ndiff
        qxk(k) = 0.0_rkx
      end do
      do k = 1 , kz
        uxk(k)   = ux(i,k)
        vxk(k)   = vx(i,k)
        txk(k)   = tx(i,k)
        thxk(k)  = thx(i,k)
        thvxk(k) = thvx(i,k)
        q2xk(k)  = q2x(i,k)
        hgame(k) = hgame2d(i,k)
      end do
      do k = 1 , kzm1
        if ( pblflg(i) .and. k <= kpbl(i) ) then
          zfacdx      = 0.2_rkx*hpbl(i)/za(i,k)
          delxy       = sqrt(dx*dy)*max(zfacdx,1.0_rkx)
          ptke1(k+1)  = ptke(delxy,hpbl(i))
        end if
      end do
      do k = 1 , kzp1
        zqk(k) = zq(i,k)
      end do
      do k = 1 , kz*ndiff
        qxk(k) = qx(i,k)
      end do
      do k = 2 , kz
        akmk(k) = xkzm(i,k-1)
        akhk(k) = xkzh(i,k-1)
        mfk(k)      = mf(i,k-1)/xkzh(i,k-1)
        ufxpblk(k)  = ufxpbl(i)*zfacent(i,k-1)/xkzm(i,k-1)
        vfxpblk(k)  = vfxpbl(i)*zfacent(i,k-1)/xkzm(i,k-1)
        qfxpblk(k)  = qfxpbl(i)*zfacent(i,k-1)/xkzq(i,k-1)
      end do
      if ( pblflg(i) ) then
        k = kpbl(i) - 1
        dex = 0.25_rkx*(q2xk(k+2)-q2xk(k))
        efxpbl(i) = we(i)*dex
      end if
      delxy = sqrt(dx*dy)
      !
      !---- find the mixing length
      !
      call mixlen(lmh,uxk,vxk,txk,thxk,qxk(1),qxk(kzp1),q2xk,zqk,     &
              ust(i),corf(i),epshol(i),s2,gh,rig,el,hpbl(i),kpbl(i),  &
              ct(i),hgamu(i),hgamv(i),hgamq(i),pblflg(i),mfk,ufxpblk, &
              vfxpblk,qfxpblk)
      !
      !---- solve for the production/dissipation of the
      !     turbulent kinetic energy
      !
      call prodq2(dt,ust(i),s2,rig,q2xk,el,zqk,akmk,akhk,     &
              uxk,vxk,thxk,thvxk,hgamu(i),hgamv(i),hgamq(i),delxy, &
              hpbl(i),pblflg(i),kpbl(i),mfk,ufxpblk,vfxpblk,qfxpblk)
      !
      !---- carry out the vertical diffusion of turbulent kinetic energy
      !
      call vdifq(lmh,dt,q2xk,zqk,akhk,ptke1,hgame,hpbl(i), &
              pblflg(i),kpbl(i),efxpbl(i))
      !
      !---- save the new tke and mixing length.
      !
      do k = 1 , kz
        q2x(i,k) = max(q2xk(k),epsq2l)
        tke(i,k) = 0.5_rkx*q2x(i,k)
      end do
    end do tke_calculation
    !
    !---- end of tke calculation
    !
    !---- end of vertical diffusion
    !
  end subroutine shinhong2d

  subroutine tridi1(cl,cm,cu,r1,au,f1,nbl)
    implicit none
    integer , intent(in) :: nbl
    real(rkx) , dimension(nbl,2:kzp1) , intent(in) :: cl
    real(rkx) , dimension(nbl,kz) , intent(in) :: cm
    real(rkx) , dimension(nbl,kz) , intent(in) :: r1
    real(rkx) , dimension(nbl,kz) , intent(inout) :: au , cu
    real(rkx) , dimension(nbl,kz) , intent(inout) :: f1
    real(rkx) :: fk
    integer :: i , k

    do i = 1 , nbl
     fk = 1.0_rkx/cm(i,1)
     au(i,1) = fk*cu(i,1)
     f1(i,1) = fk*r1(i,1)
    end do
    do k = 2 , kzm1
      do i = 1 , nbl
        fk = 1.0_rkx/(cm(i,k)-cl(i,k)*au(i,k-1))
        au(i,k) = fk*cu(i,k)
        f1(i,k) = fk*(r1(i,k)-cl(i,k)*f1(i,k-1))
      end do
    end do
    do i = 1 , nbl
      fk = 1.0_rkx/(cm(i,kz)-cl(i,kz)*au(i,kzm1))
      f1(i,kz) = fk*(r1(i,kz)-cl(i,kz)*f1(i,kzm1))
    end do
    do k = kzm1 , 1 , -1
      do i = 1 , nbl
        f1(i,k) = f1(i,k)-au(i,k)*f1(i,k+1)
      end do
    end do
  end subroutine tridi1

  subroutine tridi2(cl,cm,cu,r1,r2,au,f1,f2,nbl)
    implicit none
    integer , intent(in) :: nbl
    real(rkx) , dimension(nbl,2:kzp1) , intent(in) :: cl
    real(rkx) , dimension(nbl,kz) , intent(in) :: cm , cu , r1 , r2
    real(rkx) , dimension(nbl,kz) , intent(inout) :: au , f1 , f2
    real(rkx) :: fk
    integer :: i , k

    do i = 1 , nbl
      fk = 1.0_rkx/cm(i,1)
      au(i,1) = fk*cu(i,1)
      f1(i,1) = fk*r1(i,1)
      f2(i,1) = fk*r2(i,1)
    end do
    do k = 2 , kzm1
      do i = 1 , nbl
        fk = 1.0_rkx/(cm(i,k)-cl(i,k)*au(i,k-1))
        au(i,k) = fk*cu(i,k)
        f1(i,k) = fk*(r1(i,k)-cl(i,k)*f1(i,k-1))
        f2(i,k) = fk*(r2(i,k)-cl(i,k)*f2(i,k-1))
      end do
    end do
    do i = 1 , nbl
      fk = 1.0_rkx/(cm(i,kz)-cl(i,kz)*au(i,kzm1))
      f1(i,kz) = fk*(r1(i,kz)-cl(i,kz)*f1(i,kzm1))
      f2(i,kz) = fk*(r2(i,kz)-cl(i,kz)*f2(i,kzm1))
    end do
    do k = kzm1 , 1 , -1
      do i = 1 , nbl
        f1(i,k) = f1(i,k)-au(i,k)*f1(i,k+1)
        f2(i,k) = f2(i,k)-au(i,k)*f2(i,k+1)
      end do
    end do
  end subroutine tridi2

  subroutine tridin_ysu(cl,cm,cu,rn,au,fn,nbl,nt)
    implicit none
    integer , intent(in) :: nbl , nt
    real(rkx) , dimension(nbl,2:kzp1) , intent(in) :: cl
    real(rkx) , dimension(nbl,kz) , intent(in) :: cm
    real(rkx) , dimension(nbl,kz,nt) , intent(in) :: rn
    real(rkx) , dimension(nbl,kz) , intent(inout) :: au , cu
    real(rkx) , dimension(nbl,kz,nt) , intent(inout) :: fn
    real(rkx) :: fk
    integer :: i , k , it

    do it = 1 , nt
      do i = 1 , nbl
       fk = 1.0_rkx/cm(i,1)
       au(i,1) = fk*cu(i,1)
       fn(i,1,it) = fk*rn(i,1,it)
      end do
    end do
    do it = 1 , nt
      do k = 2 , kzm1
        do i = 1 , nbl
          fk = 1.0_rkx/(cm(i,k)-cl(i,k)*au(i,k-1))
          au(i,k) = fk*cu(i,k)
          fn(i,k,it) = fk*(rn(i,k,it)-cl(i,k)*fn(i,k-1,it))
        end do
      end do
    end do
    do it = 1 , nt
      do i = 1 , nbl
        fk = 1.0_rkx/(cm(i,kz)-cl(i,kz)*au(i,kzm1))
        fn(i,kz,it) = fk*(rn(i,kz,it)-cl(i,kz)*fn(i,kzm1,it))
      end do
    end do
    do it = 1 , nt
      do k = kzm1 , 1 , -1
        do i = 1 , nbl
          fn(i,k,it) = fn(i,k,it)-au(i,k)*fn(i,k+1,it)
        end do
      end do
    end do
  end subroutine tridin_ysu

  subroutine mixlen(lmh,u,v,t,the,q,cwm,q2,z,ustar,corf,epshol, &
                    s2,gh,ri,el,hpbl,lpbl,ct,hgamu,hgamv,hgamq, &
                    pblflg,mf,ufxpbl,vfxpbl,qfxpbl)
    implicit none
    integer , intent(in) :: lmh , lpbl
    real(rkx) , intent(in) :: hpbl , corf , ustar , hgamu  ,hgamv  ,hgamq
    real(rkx) , intent(inout) :: ct , epshol
    real(rkx) , dimension(kz) , intent(in) :: cwm , q , q2 , t
    real(rkx) , dimension(kz) , intent(in) :: the, u , v
    real(rkx) , dimension(2:kz) , intent(in) :: mf , ufxpbl , vfxpbl
    real(rkx) , dimension(2:kz) , intent(in) :: qfxpbl
    real(rkx) , dimension(kzp1) , intent(in) :: z
    real(rkx) , dimension(2:kz) , intent(out) :: el , ri , gh , s2
    logical , intent(in) :: pblflg
    integer :: k , lpblm
    real(rkx) :: suk , svk , a , aden , b , bden , aubr , bubr
    real(rkx) :: el0 , eloq2x , ghl , s2l , qol2st , qol2un
    real(rkx) :: qdzl , rdz , sq , srel , szq , tem , thm , vkrmz
    real(rkx) :: rlambda , rlb , rln , f , ckp
    real(rkx) , dimension(kz) :: q1 , en2
    real(rkx) , dimension(2:kz) :: dth , elm , rel
    !-----------------------------------------------------------------------
    !  qnse model constants
    !-----------------------------------------------------------------------
    real(rkx) , parameter :: blckdr = 0.0063_rkx
    real(rkx) , parameter :: cn = 0.75_rkx
    real(rkx) , parameter :: eps1 = 1.e-12_rkx
    real(rkx) , parameter :: epsl = 0.32_rkx
    real(rkx) , parameter :: epsru = 1.e-7_rkx
    real(rkx) , parameter :: epsrs = 1.e-7_rkx
    real(rkx) , parameter :: el0max = 1000.0_rkx
    real(rkx) , parameter :: el0min = 1.0_rkx
    real(rkx) , parameter :: elfc = 0.23_rkx * 0.5_rkx
    real(rkx) , parameter :: alph = 0.30_rkx
    real(rkx) , parameter :: beta = 1.0_rkx/273.0_rkx
    real(rkx) , parameter :: btg = beta*egrav
    real(rkx) , parameter :: a1 = 0.659888514560862645_rkx
    real(rkx) , parameter :: a2x = 0.6574209922667784586_rkx
    real(rkx) , parameter :: b1 = 11.87799326209552761_rkx
    real(rkx) , parameter :: b2 = 7.226971804046074028_rkx
    real(rkx) , parameter :: c1 = 0.000830955950095854396_rkx
    real(rkx) , parameter :: adnh = 9.0_rkx*a1*a2x*a2x * &
                                 (12.0_rkx*a1+3.0_rkx*b2)*btg*btg
    real(rkx) , parameter :: adnm = 18.0_rkx*a1*a1*a2x*(b2-3.0_rkx*a2x)*btg
    real(rkx) , parameter :: bdnh = 3.0_rkx*a2x*(7.0_rkx*a1+b2)*btg
    real(rkx) , parameter :: bdnm = 6.0_rkx*a1*a1
    !-----------------------------------------------------------------------
    !  free term in the equilibrium equation for (l/q)**2
    !-----------------------------------------------------------------------
    real(rkx) , parameter :: aeqh = 9.0_rkx*a1*a2x*a2x*b1*btg*btg + &
                                  9.0_rkx*a1*a2x*a2x* &
                                  (12.0_rkx*a1+3.0_rkx*b2)*btg*btg
    real(rkx) , parameter :: aeqm = 3.0_rkx*a1*a2x*b1*(3.0_rkx*a2x + &
                                  3.0_rkx*b2*c1+18.0_rkx*a1*c1-b2)*btg + &
                                  18._rkx*a1*a1*a2x*(b2-3.0_rkx*a2x)*btg
    !-----------------------------------------------------------------------
    !  forbidden turbulence area
    !-----------------------------------------------------------------------
    real(rkx) , parameter :: requ = -aeqh/aeqm
    real(rkx) , parameter :: epsgh = 1.e-9_rkx
    real(rkx) , parameter :: epsgm = requ*epsgh
    !-----------------------------------------------------------------------
    !  near isotropy for shear turbulence, ww/q2 lower limit
    !-----------------------------------------------------------------------
    real(rkx) , parameter :: ubryl = (18.0_rkx*requ*a1*a1*a2x*b2*c1*btg+ &
                                    9.0_rkx*a1*a2x*a2x*b2*btg*btg) / &
                                   (requ*adnm+adnh)
    real(rkx) , parameter :: ubry = (1.0_rkx+epsrs)*ubryl
    real(rkx) , parameter :: ubry3 = 3.0_rkx*ubry
    real(rkx) , parameter :: aubh = 27.0_rkx*a1*a2x*a2x*b2*btg*btg-adnh*ubry3
    real(rkx) , parameter :: aubm = 54.0_rkx*a1*a1*a2x*b2*c1*btg  -adnm*ubry3
    real(rkx) , parameter :: bubh = (9.0_rkx*a1*a2x+3.0_rkx*a2x*b2) * &
                                  btg-bdnh*ubry3
    real(rkx) , parameter :: bubm = 18.0_rkx*a1*a1*c1-bdnm*ubry3
    real(rkx) , parameter :: cubr = 1.0_rkx-ubry3
    real(rkx) , parameter :: rcubr = 1.0_rkx/cubr
    !-----------------------------------------------------------------------
    !  k profile constants
    !-----------------------------------------------------------------------
    real(rkx) , parameter :: elcbl = 0.77_rkx
    real(rkx) , parameter :: elocp = 2.72e6_rkx * cpd

    ct = 0.0_rkx
    do k = 1 , kz
      q1(k) = 0.0_rkx
    end do
    do k = 2 , kz
      dth(k) = the(k)-the(k-1)
    end do
    do k = 3 , kz
      if ( dth(k) > 0.0_rkx .and. dth(k-1) <= 0.0_rkx ) then
        dth(k) = dth(k)+ct
        exit
      end if
    end do
    !
    !  compute local gradient richardson number
    !
    do k = kz , 2 , -1
      rdz = 2.0_rkx/(z(k+1)-z(k-1))
      s2l = ((u(k)-u(k-1))**2+(v(k)-v(k-1))**2)*rdz*rdz ! s**2
      if ( pblflg .and. k <= lpbl ) then
        suk = (u(k)-u(k-1))*rdz
        svk = (v(k)-v(k-1))*rdz
        s2l = (suk-hgamu/hpbl-ufxpbl(k))*suk+(svk-hgamv/hpbl-vfxpbl(k))*svk
      end if
      s2l = max(s2l,epsgm)
      s2(k) = s2l
      tem = (t(k)+t(k-1))*0.5_rkx
      thm = (the(k)+the(k-1))*0.5_rkx
      a = thm*ep1
      b = (elocp/tem-1.0_rkx-ep1)*thm
      ghl = (dth(k)*((q(k)+q(k-1)+cwm(k)+cwm(k-1))*(0.5_rkx*ep1)+1.0_rkx) + &
            (q(k)-q(k-1)+cwm(k)-cwm(k-1))*a + &
            (cwm(k)-cwm(k-1))*b)*rdz    ! dtheta/dz
      if ( pblflg .and. k <= lpbl ) then
        ghl = ghl-mf(k)-(hgamq/hpbl+qfxpbl(k))*a
      end if
      if ( abs(ghl) <= epsgh ) ghl = epsgh
      en2(k) = ghl*egrav/thm ! n**2
      gh(k) = ghl
      ri(k) = en2(k)/s2l
    end do
    !
    !  find maximum mixing lengths and the level of the pbl top
    !
    do k = kz , 2 , -1
      s2l = s2(k)
      ghl = gh(k)
      if ( ghl >= epsgh ) then
        if ( s2l/ghl <= requ ) then
          elm(k) = epsl
        else
          aubr = (aubm*s2l+aubh*ghl)*ghl
          bubr = bubm*s2l+bubh*ghl
          qol2st = (-0.5_rkx*bubr+sqrt(bubr*bubr*0.25_rkx-aubr*cubr))*rcubr
          eloq2x = 1.0_rkx/qol2st
          elm(k) = max(sqrt(eloq2x*q2(k)),epsl)
        end if
      else
        aden = (adnm*s2l+adnh*ghl)*ghl
        bden = bdnm*s2l+bdnh*ghl
        qol2un = -0.5_rkx*bden+sqrt(bden*bden*0.25_rkx-aden)
        eloq2x = 1.0_rkx/(qol2un+epsru)       ! repsr1/qol2un
        elm(k) = max(sqrt(eloq2x*q2(k)),epsl)
      end if
    end do
    do k = lpbl , lmh , -1
      q1(k) = sqrt(q2(k))
    end do
    szq = 0.0_rkx
    sq = 0.0_rkx
    do k = kz , 2 , -1
      qdzl = (q1(k)+q1(k-1))*(z(k)-z(k-1))
      szq = (z(k)+z(k-1)-z(lmh)-z(lmh))*qdzl+szq
      sq = qdzl+sq
    end do
    !
    !  computation of asymptotic l in blackadar formula
    !
    el0 = min(alph*szq*0.5_rkx/sq,el0max)
    el0 = max(el0,el0min)
    !
    !  above the pbl top
    !
    lpblm = min(lpbl+1,kz)
    do k = kz , lpblm , -1
      el(k) = (z(k+1)-z(k-1))*elfc
      rel(k) = el(k)/elm(k)
    end do
    !
    !  inside the pbl
    !
    epshol = min(epshol,0.0_rkx)
    ckp = elcbl*((1.0_rkx-8.0_rkx*epshol)**(1.0_rkx/3.0_rkx))
    if ( lpbl > lmh ) then
      do k = lpbl , lmh+1  ,-1
        vkrmz = (z(k)-z(lmh))*vonkar
        if ( pblflg ) then
          vkrmz = ckp*(z(k)-z(lmh))*vonkar
          el(k)=vkrmz/(vkrmz/el0+1.0_rkx)
        else
          el(k) = vkrmz/(vkrmz/el0+1.0_rkx)
        end if
        rel(k) = el(k)/elm(k)
      end do
    end if
    do k = lpbl-1 , lmh+2 , -1
      srel = min(((rel(k-1)+rel(k+1))*0.5_rkx+rel(k))*0.5_rkx,rel(k))
      el(k) = max(srel*elm(k),epsl)
    end do
    !
    !  mixing length for the qnse model in stable case
    !
    f = max(corf,eps1)
    rlambda = f / (blckdr*ustar)
    do k = kz , 2, -1
      if ( en2(k) >= 0.0_rkx ) then ! stable case
        vkrmz = (z(k)-z(lmh))*vonkar
        rlb = rlambda+1.0_rkx/vkrmz
        rln = sqrt(2.0_rkx*en2(k)/q2(k))/cn
        el(k) = 1.0_rkx/(rlb+rln)
      end if
    end do
  end subroutine mixlen

  subroutine prodq2(dtturbl,ustar,s2,ri,q2,el,z,akm,akh,    &
                    uxk,vxk,thxk,thvxk,hgamu,hgamv,hgamq,delxy, &
                    hpbl,pblflg,kpbl,mf,ufxpbl,vfxpbl,qfxpbl)
    implicit none
    integer , intent(in) :: kpbl
    real(rkx) , intent(in) :: dtturbl , ustar
    real(rkx) , intent(in) :: hgamu , hgamv , hgamq , delxy , hpbl
    logical,  intent(in) :: pblflg
    real(rkx) , dimension(kz) , intent(in) :: uxk , vxk
    real(rkx) , dimension(kz) , intent(in) :: thxk , thvxk
    real(rkx) , dimension(2:kz) , intent(in) :: s2 , ri , akm, akh
    real(rkx) , dimension(2:kz) , intent(in) :: el , mf
    real(rkx) , dimension(2:kz) , intent(in) :: ufxpbl , vfxpbl , qfxpbl
    real(rkx) , dimension(kzp1) , intent(in) :: z
    real(rkx) , dimension(kz) , intent(inout) :: q2
    integer :: k
    real(rkx) :: s2l , q2l , deltaz , akml , akhl , en2 , pr , bpr
    real(rkx) :: dis , suk , svk , gthvk , govrthvk , pru , prv
    real(rkx) :: thm , disel
    real(rkx) , parameter :: epsq2l = 0.01_rkx
    real(rkx) , parameter :: c0 = 0.55_rkx
    real(rkx) , parameter :: ceps = 16.6_rkx
    real(rkx) , parameter :: rc02 = 2.0_rkx/(c0*c0)
    !
    !  start of production/dissipation loop
    !
    main_integration: &
    do k = 2 , kz
      deltaz = 0.5_rkx*(z(k+1)-z(k-1))
      s2l = s2(k)
      q2l = q2(k)
      suk = (uxk(k)-uxk(k-1))/deltaz
      svk = (vxk(k)-vxk(k-1))/deltaz
      gthvk = (thvxk(k)-thvxk(k-1))/deltaz
      govrthvk = egrav/(0.5_rkx*(thvxk(k)+thvxk(k-1)))
      akml = akm(k)
      akhl = akh(k)
      en2 = ri(k)*s2l !n**2
      thm = (thxk(k)+thxk(k-1))*0.5_rkx
      !
      !  turbulence production term
      !
      if ( pblflg .and. k <= kpbl ) then
        pru = (akml*(suk-hgamu/hpbl-ufxpbl(k)))*suk
        prv = (akml*(svk-hgamv/hpbl-vfxpbl(k)))*svk
      else
        pru = akml*suk*suk
        prv = akml*svk*svk
      end if
      pr = pru+prv
      !
      !  buoyancy production
      !
      if ( pblflg .and. k <= kpbl ) then
        bpr = (akhl*(gthvk-mf(k)-(hgamq/hpbl+qfxpbl(k))*ep1*thm))*govrthvk
      else
        bpr = akhl*gthvk*govrthvk
      end if
      !
      !  dissipation
      !
      disel = min(delxy,ceps*el(k))
      dis = (q2l)**1.5_rkx/disel
      q2l = q2l + 2.0_rkx*(pr-bpr-dis)*dtturbl
      q2(k) = max(q2l,epsq2l)
      !
      !  end of production/dissipation loop
      !
    end do main_integration
    !
    !  lower boundary condition for q2
    !
    q2(1) = max(rc02*ustar*ustar,epsq2l)
  end subroutine prodq2

  subroutine vdifq(lmh,dtdif,q2,z,akhk,ptke1,hgame,hpbl, &
                   pblflg,kpbl,efxpbl)
    implicit none
    integer , intent(in) :: lmh , kpbl
    real(rkx) , intent(in) :: dtdif , hpbl , efxpbl
    logical , intent(in) :: pblflg
    real(rkx) , dimension(kz) , intent(in) :: hgame , ptke1
    real(rkx) , dimension(2:kz) , intent(in) :: akhk
    real(rkx) , dimension(kzp1) , intent(in) :: z
    real(rkx) , dimension(kz) , intent(inout) :: q2
    real(rkx) , parameter :: c_k = 1.0_rkx
    real(rkx) , parameter :: esq = 5.0_rkx
    real(rkx) :: akqs , cf , dtozs
    real(rkx) :: esqhf , zak
    real(rkx) , dimension(2:kz) :: zfacentk
    real(rkx) , dimension(3:kz) :: akq , cm , cr , dtoz , rsq2
    integer :: k
    !
    !  vertical turbulent diffusion
    !
    esqhf = 0.5_rkx*esq
    do k = 2 , kz
      zak = 0.5_rkx*(z(k)+z(k-1)) !zak of vdifq = za(k-1) of shinhong2d
      zfacentk(k) = (zak/hpbl)**3
    end do
    do k = kz , 3, -1
      dtoz(k) = (dtdif+dtdif)/(z(k+1)-z(k-1))
      akq(k) = c_k*(akhk(k)/(z(k+1)-z(k-1))+akhk(k-1)/(z(k)-z(k-2)))
      akq(k) = akq(k)*ptke1(k)
      cr(k) = -dtoz(k)*akq(k)
    end do
    akqs = c_k*akhk(2)/(z(3)-z(1))
    akqs = akqs*ptke1(2)
    cm(kz) = dtoz(kz)*akq(kz)+1.0_rkx
    rsq2(kz) = q2(kz)
    do k = kzm1 , 3 , -1
      cf = -dtoz(k)*akq(k+1)/cm(k+1)
      cm(k) = -cr(k+1)*cf+(akq(k+1)+akq(k))*dtoz(k)+1.0_rkx
      rsq2(k) = -rsq2(k+1)*cf+q2(k)
      if ( pblflg .and. k < kpbl ) then
        rsq2(k) = rsq2(k)-dtoz(k)*(2.0_rkx*hgame(k)/hpbl)*akq(k+1) * &
          (z(k+1)-z(k))+dtoz(k)*(2.0_rkx*hgame(k-1)/hpbl)*akq(k)*(z(k)-z(k-1))
        rsq2(k) = rsq2(k)-dtoz(k)*2.0_rkx*efxpbl*zfacentk(k+1) + &
          dtoz(k)*2.0_rkx*efxpbl*zfacentk(k)
      end if
    end do
    dtozs = (dtdif+dtdif)/(z(3)-z(1))
    cf = -dtozs*akq(lmh+2)/cm(lmh+2)
    if ( pblflg .and. ((lmh+1) < kpbl) ) then
      q2(lmh+1) = (dtozs*akqs*q2(lmh)-rsq2(lmh+2)*cf+q2(lmh+1) - &
         dtozs*(2.0_rkx*hgame(lmh+1)/hpbl)*akq(lmh+2)*(z(lmh+2)-z(lmh+1)) + &
         dtozs*(2.0_rkx*hgame(lmh)/hpbl)*akqs*(z(lmh+1)-z(lmh)))
      q2(lmh+1) = q2(lmh+1)-dtozs*2.0_rkx*efxpbl*zfacentk(lmh+2) + &
         dtozs*2.0*efxpbl*zfacentk(lmh+1)
      q2(lmh+1) = q2(lmh+1)/((akq(lmh+2)+akqs)*dtozs-cr(lmh+2)*cf+1.0_rkx)
    else
      q2(lmh+1) = (dtozs*akqs*q2(lmh)-rsq2(lmh+2)*cf+q2(lmh+1)) / &
         ((akq(lmh+2)+akqs)*dtozs-cr(lmh+2)*cf+1.0_rkx)
    end if
    do k = lmh+2 , kz
      q2(k) = (-cr(k)*q2(k-1)+rsq2(k))/cm(k)
    end do
  end subroutine vdifq

  real(rkx) function pu(d,h)
    implicit none
    real(rkx) , intent(in) :: d , h
    real(rkx) , parameter :: pmin = 0.0_rkx
    real(rkx) , parameter :: pmax = 1.0_rkx
    real(rkx) , parameter :: a1 = 1.0_rkx
    real(rkx) , parameter :: a2 = 0.070_rkx
    real(rkx) , parameter :: a3 = 1.0_rkx
    real(rkx) , parameter :: a4 = 0.142_rkx
    real(rkx) , parameter :: a5 = 0.071_rkx
    real(rkx) , parameter :: b1 = 2.0_rkx
    real(rkx) , parameter :: b2 = 0.6666667_rkx
    real(rkx) :: doh , num , den

    if ( abs(h) > tiny(h) ) then
      doh = d/h
      num = a1*(doh)**b1+a2*(doh)**b2
      den = a3*(doh)**b1+a4*(doh)**b2+a5
      pu = num/den
    else
      pu = 1.0_rkx
    end if
    pu = max(pu,pmin)
    pu = min(pu,pmax)
  end function pu

  real(rkx) function pq(d,h)
    implicit none
    real(rkx) , intent(in) :: d , h
    real(rkx) , parameter :: pmin = 0.0_rkx
    real(rkx) , parameter :: pmax = 1.0_rkx
    real(rkx) , parameter :: a1 = 1.0_rkx
    real(rkx) , parameter :: a2 = -0.098_rkx
    real(rkx) , parameter :: a3 = 1.0_rkx
    real(rkx) , parameter :: a4 = 0.106_rkx
    real(rkx) , parameter :: a5 = 0.5_rkx
    real(rkx) , parameter :: b1 = 2.0_rkx
    real(rkx) :: doh , num , den

    if ( abs(h) > tiny(h) ) then
      doh = d/h
      num = a1*(doh)**b1+a2
      den = a3*(doh)**b1+a4
      pq = a5*num/den+(1.0_rkx-a5)
    else
      pq = 1.0_rkx
    end if
    pq = max(pq,pmin)
    pq = min(pq,pmax)
  end function pq

  real(rkx) function pthnl(d,h)
    implicit none
    real(rkx) , intent(in) :: d , h
    real(rkx) , parameter :: pmin = 0.0_rkx
    real(rkx) , parameter :: pmax = 1.0_rkx
    real(rkx) , parameter :: a1 = 1.000_rkx
    real(rkx) , parameter :: a2 = 0.936_rkx
    real(rkx) , parameter :: a3 = -1.110_rkx
    real(rkx) , parameter :: a4 = 1.000_rkx
    real(rkx) , parameter :: a5 = 0.312_rkx
    real(rkx) , parameter :: a6 = 0.329_rkx
    real(rkx) , parameter :: a7 = 0.243_rkx
    real(rkx) , parameter :: b1 = 2.0_rkx
    real(rkx) , parameter :: b2 = 0.875_rkx
    real(rkx) :: doh , num , den

    if ( abs(h) > tiny(h) ) then
      doh = d/h
      num = a1*(doh)**b1+a2*(doh)**b2+a3
      den = a4*(doh)**b1+a5*(doh)**b2+a6
      pthnl = a7*num/den+(1.0_rkx-a7)
    else
      pthnl = 1.0_rkx
    end if
    pthnl = max(pthnl,pmin)
    pthnl = min(pthnl,pmax)
  end function pthnl

  real(rkx) function pthl(d,h)
    implicit none
    real(rkx) , intent(in) :: d , h
    real(rkx) , parameter :: pmin = 0.0_rkx
    real(rkx) , parameter :: pmax = 1.0_rkx
    real(rkx) , parameter :: a1 = 1.000_rkx
    real(rkx) , parameter :: a2 = 0.870_rkx
    real(rkx) , parameter :: a3 = -0.913_rkx
    real(rkx) , parameter :: a4 = 1.000_rkx
    real(rkx) , parameter :: a5 = 0.153_rkx
    real(rkx) , parameter :: a6 = 0.278_rkx
    real(rkx) , parameter :: a7 = 0.280_rkx
    real(rkx) , parameter :: b1 = 2.0_rkx
    real(rkx) , parameter :: b2 = 0.5_rkx
    real(rkx) :: doh , num , den

    if ( abs(h) > tiny(h) ) then
      doh = d/h
      num = a1*(doh)**b1+a2*(doh)**b2+a3
      den = a4*(doh)**b1+a5*(doh)**b2+a6
      pthl = a7*num/den+(1.0_rkx-a7)
    else
      pthl = 1.0_rkx
    end if
    pthl = max(pthl,pmin)
    pthl = min(pthl,pmax)
  end function pthl

  real(rkx) function ptke(d,h)
    implicit none
    real(rkx) , intent(in) :: d , h
    real(rkx) , parameter :: pmin = 0.0_rkx
    real(rkx) , parameter :: pmax = 1.0_rkx
    real(rkx) , parameter :: a1 = 1.000_rkx
    real(rkx) , parameter :: a2 = 0.070_rkx
    real(rkx) , parameter :: a3 = 1.000_rkx
    real(rkx) , parameter :: a4 = 0.142_rkx
    real(rkx) , parameter :: a5 = 0.071_rkx
    real(rkx) , parameter :: b1 = 2.0_rkx
    real(rkx) , parameter :: b2 = 0.6666667_rkx
    real(rkx) :: doh , num , den

    if ( abs(h) > tiny(h) ) then
      doh = d/h
      num = a1*(doh)**b1+a2*(doh)**b2
      den = a3*(doh)**b1+a4*(doh)**b2+a5
      ptke = num/den
    else
      ptke = 1.0_rkx
    end if
    ptke = max(ptke,pmin)
    ptke = min(ptke,pmax)
  end function ptke

  subroutine sfclay(nbl,br,znt,za,ust,govrth,dtg,rah,psim,psih)
    implicit none
    integer , intent(in) :: nbl
    real(rkx) , dimension(nbl) , intent(in) :: br , znt , za , ust
    real(rkx) , dimension(nbl) , intent(in) :: govrth , rah , dtg
    real(rkx) , dimension(nbl) , intent(out) :: psim , psih
    real(rkx) :: gz1oz0 , zol , rzol , mol
    integer :: i , nzol
    !
    ! Diagnose basic parameters for the appropriated stability class:
    !
    ! The stability classes are determined by br (bulk richardson no.)
    ! and hol (height of pbl/monin-obukhov length).
    !
    ! Criteria for the classes are as follows:
    !
    ! 1. br .ge. 0.2;
    !    represents nighttime stable conditions (regime=1),
    !
    ! 2. br .lt. 0.2 .and. br .gt. 0.0;
    !    represents damped mechanical turbulent conditions
    !    (regime=2),
    !
    ! 3. br .eq. 0.0
    !    represents forced convection conditions (regime=3),
    !
    ! 4. br .lt. 0.0
    !    represents free convection conditions (regime=4).
    !
    do i = 1 , nbl
      gz1oz0 = log(za(i)/znt(i))
      if ( br(i) >= 0.2_rkx ) then
        !
        ! Class 1; stable (nighttime) conditions:
        !
        psim(i) = -10.0_rkx*gz1oz0
        psim(i) = max(psim(i),-10.0_rkx)
        psih(i) = psim(i)
      else if ( br(i) > 0.0 ) then
        !
        ! Class 2; damped mechanical turbulence:
        !
        psim(i) = -5.0_rkx*br(i)*gz1oz0/(1.1_rkx-5.0_rkx*br(i))
        psim(i) = max(psim(i),-10.0_rkx)
        psih(i) = psim(i)
      else if ( br(i) == 0.0 ) then
        !
        ! Class 3; forced convection:
        !
        psim(i) = 0.0_rkx
        psih(i) = psim(i)
      else
        !
        ! Class 4; free convection:
        !
        if ( ust(i) < 0.01_rkx ) then
          zol = br(i)*gz1oz0
        else
          mol = dtg(i)/(rah(i)*ust(i))
          zol = vonkar*govrth(i)*za(i)*mol/(ust(i)*ust(i))
        end if
        zol = min(zol,0.0_rkx)
        zol = max(zol,-9.9999_rkx)
        nzol = int(-zol*100.0_rkx)
        rzol = -zol*100.0_rkx-nzol
        psim(i) = psimtb(nzol)+rzol*(psimtb(nzol+1)-psimtb(nzol))
        psih(i) = psihtb(nzol)+rzol*(psihtb(nzol+1)-psihtb(nzol))
        psih(i) = min(psih(i),0.9_rkx*gz1oz0)
        psim(i) = min(psim(i),0.9_rkx*gz1oz0)
      end if
    end do
  end subroutine sfclay

  subroutine sfclayinit
    implicit none
    integer :: n
    real(rkx) :: zoln , x , y
    do n = 0 , 1000
      zoln = -real(n,rkx)*0.01_rkx
      x = (1.0_rkx - 16.0_rkx*zoln)**0.25_rkx
      psimtb(n) = 2.0_rkx*log(0.5_rkx*(1_rkx+x)) + &
        log(0.5_rkx*(1.0_rkx+x*x))-2.0_rkx*atan(x)+2.0_rkx*atan(1.0_rkx)
      y = (1.0_rkx - 16.0_rkx*zoln)**0.5_rkx
      psihtb(n) = 2.0_rkx*log(0.5_rkx*(1.0_rkx+y))
    end do
  end subroutine sfclayinit

end module mod_pbl_shinhong

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

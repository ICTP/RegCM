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
   
  module mod_pbl_holtbl
  !
  ! Holtslag planetary boundary layer scheme
  ! Reference : Holtslag, De Bruijn and Pan - MWR - 8/90
  !
    use mod_dynparam
    use mod_mppparam
    use mod_memutil
    use mod_service
    use mod_mpmessage
    use mod_pbl_common
  !
    private
  !
    public :: allocate_mod_pbl_holtbl , holtbl
  !
    real(dp) , pointer , dimension(:,:,:) :: vv , cgh , kvc , kvh ,   &
                                            kvm , kvq ! , cgq
    real(dp) , pointer, dimension(:,:) :: hfxv , obklen , th10 , &
                                         ustr , xhfx , xqfx , pfcor
  !
    real(dp) , pointer , dimension(:,:,:) :: alphak , betak , &
                          coef1 , coef2 , coef3 , coefe , coeff1 , &
                          coeff2 , tpred1 , tpred2
    real(dp) , pointer , dimension(:,:,:) :: kzm , rc , ttnp , vdep
    real(dp) , pointer , dimension(:,:) :: govrth
    real(dp) , pointer , dimension(:) :: hydf
  !
    real(dp) , pointer , dimension(:,:,:) :: dza , thvx
    real(dp) , pointer , dimension(:,:,:) :: akzz1 , akzz2
    real(dp) , pointer , dimension(:,:,:) :: rhohf
  !
    real(dp) , pointer , dimension(:,:,:) :: ri
    real(dp) , pointer , dimension(:,:) :: therm
  !
  ! minimum eddy diffusivity
    real(dp) , parameter :: kzo = d_one
    real(dp) , parameter :: szkm = 1600.0D0
  ! coef. of proportionality and lower % of bl in sfc layer
    real(dp) , parameter :: fak = 8.5D0
    real(dp) , parameter :: sffrac = 0.1D0
  ! beta coefs. for momentum, stable conditions and heat
    real(dp) , parameter :: betam = 15.0D0
    real(dp) , parameter :: betas = 5.0D0
    real(dp) , parameter :: betah = 15.0D0
    real(dp) , parameter :: mult = 0.61D0
    real(dp) , parameter :: ccon = fak*sffrac*vonkar
    real(dp) , parameter :: gvk = egrav*vonkar
    real(dp) , parameter :: gpcf = egrav/d_1000 ! Grav and pressure conversion
    real(dp) , parameter :: binm = betam*sffrac
    real(dp) , parameter :: binh = betah*sffrac
  ! power in formula for k and critical ri for judging stability
    real(dp) , parameter :: pink = d_two
    real(dp) , parameter :: ricr = d_rfour
  !
    contains
  !
    subroutine allocate_mod_pbl_holtbl(ichem,ichdrdepo)
      implicit none
      integer , intent(in) :: ichem , ichdrdepo

      call getmem3d(alphak,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:alphak')
      call getmem3d(betak,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:betak')
      call getmem3d(coef1,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:coef1')
      call getmem3d(coef2,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:coef2')
      call getmem3d(coef3,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:coef3')
      call getmem3d(coefe,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:coefe')
      call getmem3d(coeff1,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:coeff1')
      call getmem3d(coeff2,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:coeff2')
      call getmem3d(tpred1,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:tpred1')
      call getmem3d(tpred2,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:tpred2')
      call getmem3d(ri,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:ri')
      call getmem3d(kzm,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:kzm')
      call getmem3d(rc,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:rc')
      call getmem3d(ttnp,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:ttnp')
      call getmem2d(govrth,jci1,jci2,ici1,ici2,'mod_holtbl:govrth')
      call getmem1d(hydf,1,kz,'mod_holtbl:hydf')
      call getmem2d(therm,jci1,jci2,ici1,ici2,'mod_holtbl:therm')
      call getmem3d(vv,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:vv')
      call getmem3d(cgh,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:cgh')
      call getmem3d(kvh,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:kvh')
      call getmem3d(kvm,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:kvm')
      call getmem3d(kvq,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:kvq')
      call getmem2d(hfxv,jci1,jci2,ici1,ici2,'mod_holtbl:hfxv')
      call getmem2d(pfcor,jci1,jci2,ici1,ici2,'mod_holtbl:pfcor')
      call getmem2d(obklen,jci1,jci2,ici1,ici2,'mod_holtbl:obklen')
      call getmem2d(th10,jci1,jci2,ici1,ici2,'mod_holtbl:th10')
      call getmem2d(ustr,jci1,jci2,ici1,ici2,'mod_holtbl:ustr')
      call getmem2d(xhfx,jci1,jci2,ici1,ici2,'mod_holtbl:xhfx')
      call getmem2d(xqfx,jci1,jci2,ici1,ici2,'mod_holtbl:xqfx')
      call getmem3d(dza,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:dza')
      call getmem3d(rhohf,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:rhohf')
      call getmem3d(akzz1,jce1-ma%jbl1,jce2, &
                          ice1-ma%ibb1,ice2,1,kz,'mod_holtbl:akzz1')
      call getmem3d(akzz2,jce1-ma%jbl1,jce2, &
                          ice1-ma%ibb1,ice2,1,kz,'mod_holtbl:akzz2')
      call getmem3d(thvx,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:thvx')

      if ( ichem == 1 ) then
        if ( ichdrdepo == 1 ) then
          call getmem3d(vdep,jci1,jci2,ici1,ici2,1,ntr,'mod_holtbl:vdep')
        end if
        call getmem3d(kvc,jci1,jci2,ici1,ici2,1,kz,'mod_holtbl:kvc')
    end if
  end  subroutine allocate_mod_pbl_holtbl
!
  subroutine holtbl(jstart,jend,istart,iend)
!
#ifndef IBM
  use mpi
#endif
  implicit none
#ifdef IBM
  include 'mpif.h'
#endif 
  integer , intent(in) :: jstart , jend , istart , iend
!
  real(dp) :: drgdot , kzmax , oblen , xps , ps2 , ri , &
             sf , sh10 , ss , uflxsf , uflxsfx , vflxsf ,     &
             vflxsfx
  integer :: i , j , k , itr
  character (len=64) :: subroutine_name='holtbl'
  integer :: idindx=0
!
  call time_begin(subroutine_name,idindx)
!
!----------------------------------------------------------------------
!
  ! decouple flux-form variables to give u,v,t,theta,theta-vir.,
  ! t-vir., qv, and qc at cross points and at ktau-1.
  !
  ! *** note ***
  ! the boundary winds may not be adequately affected by friction,
  ! so use only interior values of ubx3d and vbx3d to calculate
  ! tendencies.
  !
  do k = 1 , kz
    hydf(k) = gpcf/dlev(k)
  end do
  do k = 1 , kz
    do i = istart , iend
      do j = jstart , jend
        thvx(j,i,k) = thxatm(j,i,k)*(d_one+ep1*qvatm(j,i,k))
      end do
    end do
  end do
  !
  ! density at surface is stored in rhox2d(j,i), at half levels in rhohf(j,i,k).
  !
  do k = 1 , kzm1
    do i = istart , iend
      do j = jstart , jend
        dza(j,i,k) = za(j,i,k) - za(j,i,k+1)
        xps = (hlev(k)*sfcps(j,i)+ptp)*d_1000
        ps2 = (hlev(k+1)*sfcps(j,i)+ptp)*d_1000
        rhohf(j,i,k) = (ps2-xps)/(egrav*dza(j,i,k))
      end do
    end do
  end do
  do i = istart , iend
    do j = jstart , jend
      govrth(j,i) = egrav/thxatm(j,i,kz)
    end do
  end do
!
! *********************************************************************
!
  !
  !   compute the vertical diffusion term:
  !
  do k = 2 , kz
    do i = istart , iend
      do j = jstart , jend
        rc(j,i,k) = 0.257D0*dzq(j,i,k)**0.175D0
      end do
    end do
  end do
  !
  !   compute the diffusion coefficient:
  !   blackadar scheme above boundary layer top
  !
  do k = 2 , kz
    do i = istart , iend
      do j = jstart , jend
        kzmax = 0.8D0*dza(j,i,k-1)*dzq(j,i,k)*rdtpbl
        vv(j,i,k) = uatm(j,i,k)*uatm(j,i,k) + vatm(j,i,k)*vatm(j,i,k)
        ss = ((uatm(j,i,k-1)-uatm(j,i,k))*   &
              (uatm(j,i,k-1)-uatm(j,i,k))+   &
              (vatm(j,i,k-1)-vatm(j,i,k))*   &
              (vatm(j,i,k-1)-vatm(j,i,k)))/  &
              (dza(j,i,k-1)*dza(j,i,k-1)) + 1.0D-9
        ri = govrth(j,i)*(thvx(j,i,k-1)-thvx(j,i,k))/(ss*dza(j,i,k-1))
        if ( (ri-rc(j,i,k)) >= d_zero ) then
          kzm(j,i,k) = kzo
        else
          kzm(j,i,k) = kzo + dsqrt(ss)*(rc(j,i,k)-ri)*szkm/rc(j,i,k)
        end if
        kzm(j,i,k) = dmin1(kzm(j,i,k),kzmax)
      end do
    end do
  end do
!
! *********************************************************************
!
  !   holtslag pbl
  !
  !   initialize bl diffusion coefficients and counter-gradient terms
  !   with free atmosphere values and make specific humidity
  !
  do k = 2 , kz
    do i = istart , iend
      do j = jstart , jend
        ! eddy diffusivities for momentum, heat and moisture
        kvm(j,i,k) = kzm(j,i,k)
        kvh(j,i,k) = kzm(j,i,k)
        kvq(j,i,k) = kzm(j,i,k)
        if ( lchem ) then
          kvc(j,i,k) = kzm(j,i,k)
        end if
        ! counter gradient terms for heat and moisture
        cgh(j,i,k) = d_zero
      end do
    end do
  end do
 
  do i = istart , iend
    do j = jstart , jend
      ! compute friction velocity
      uflxsfx = uvdrag(j,i)*uatm(j,i,kz)
      vflxsfx = uvdrag(j,i)*vatm(j,i,kz)
      ustr(j,i) = dsqrt(dsqrt(uflxsfx*uflxsfx+vflxsfx*vflxsfx)/rhox2d(j,i))
      ! convert surface fluxes to kinematic units
      xhfx(j,i) = hfx(j,i)/(cpd*rhox2d(j,i))
      xqfx(j,i) = qfx(j,i)/rhox2d(j,i)
      ! compute virtual heat flux at surface
      hfxv(j,i) = xhfx(j,i) + mult*thxatm(j,i,kz)*xqfx(j,i)
      ! limit coriolis parameter to value at 10 deg. latitude
      pfcor(j,i) = dmax1(dabs(coriolis(j,i)),2.546D-5)
    end do
  end do
  !
  !   estimate potential temperature at 10m via log temperature
  !   profile in the surface layer (brutsaert, p. 63).
  !   calculate mixing ratio at 10m by assuming a constant
  !   value from the surface to the lowest model level.
  !
  do i = istart , iend
    do j = jstart , jend
!     th10(j,i) = ((thxatm(j,i,kz)+tg(j,i))*d_half)*(d_one+mult*sh10)
!     th10(j,i) = thvx(j,i,kz) + hfxv(j,i)/(vonkar*ustr(j,i)* &
!                 dlog(za(j,i,kz)*d_r10)
      sh10 = qvatm(j,i,kz)/(qvatm(j,i,kz)+d_one)
      ! "virtual" potential temperature
      if ( hfxv(j,i) >= d_zero ) then
        th10(j,i) = thvx(j,i,kz)
      else
        ! first approximation for obhukov length
        ! th10(j,i) = (0.25*thxatm(j,i,kz)+0.75*tg(j,i))*(d_one+mult*sh10)
        oblen = -d_half*(thxatm(j,i,kz)+tg(j,i)) *  &
                (d_one+mult*sh10)*ustr(j,i)**d_three /  &
                (gvk*(hfxv(j,i)+dsign(1.0D-10,hfxv(j,i))))
        if ( oblen >= za(j,i,kz) ) then
          th10(j,i) = thvx(j,i,kz) + hfxv(j,i)/(vonkar*ustr(j,i))*  &
             (dlog(za(j,i,kz)*d_r10)+d_five/oblen*(za(j,i,kz)-d_10))
        else if ( oblen < za(j,i,kz) .and. oblen > d_10 ) then
          th10(j,i) = thvx(j,i,kz) + hfxv(j,i)/(vonkar*ustr(j,i))*  &
              (dlog(oblen*d_r10)+d_five/oblen*(oblen-d_10)+         &
              6.0D0*dlog(za(j,i,kz)/oblen))
        else if ( oblen <= d_10 ) then
          th10(j,i) = thvx(j,i,kz) + hfxv(j,i)/(vonkar*ustr(j,i)) * &
                      6.0D0*dlog(za(j,i,kz)*d_r10)
        end if
        th10(j,i) = dmax1(th10(j,i),tg(j,i))
!gtb    th10(j,i) = dmin1(th10(j,i),tg(j,i))  ! gtb add to minimize
      end if
      ! obklen compute obukhov length
      obklen(j,i) = -th10(j,i)*ustr(j,i)**d_three / &
              (gvk*(hfxv(j,i)+dsign(1.0D-10,hfxv(j,i))))
    end do
  end do
  !
  ! compute diffusivities and counter gradient terms
  !
  call blhnew(jstart,jend,istart,iend)

  do k = 2 , kz
    do i = istart , iend
      do j = jstart , jend
        akzz1(j,i,k) = rhohf(j,i,k-1)*kvm(j,i,k)/dza(j,i,k-1)
      end do
    end do
  end do
  do k = 1 , kz
    do i = istart , iend
      do j = jstart , jend
        akzz2(j,i,k) = hydf(k)/sfcps(j,i)
      end do
    end do
  end do

  call deco1_exchange_left(akzz1,1,ice1,ice2,1,kz)
  call deco1_exchange_left(akzz2,1,ice1,ice2,1,kz)

  !
  !   calculate coefficients at dot points for u and v wind
  !
  do k = 2 , kz
    do i = istart , iend
      do j = jstart , jend
        betak(j,i,k) = (akzz1(j-1,i,k)+akzz1(j-1,i-1,k)+ &
                        akzz1(j,i,k)  +akzz1(j,i-1,k))*d_rfour
      end do
    end do
  end do
  do k = 1 , kz
    do i = istart , iend
      do j = jstart , jend
        alphak(j,i,k) = (akzz2(j-1,i,k)+akzz2(j-1,i-1,k)+ &
                         akzz2(j,i,k)  +akzz2(j,i-1,k))*d_rfour
      end do
    end do
  end do
!
! **********************************************************************
!
  ! start now procedure for implicit diffusion calculations
  ! performed separately for wind (dot points)
  ! and temperature and water vapor (cross points)
  ! countergradient term is not included in the implicit diffusion
  ! scheme its effect is included as in the old explicit scheme
  ! calculations assume fluxes positive upward, so the sign in front
  ! of uflxsf and vflxsf has been changed in the various terms
  !
  ! wind components
  !
  ! first compute coefficients of the tridiagonal matrix
  !
  ! Atmosphere top
  do i = istart , iend
    do j = jstart , jend
      coef1(j,i,1) = dtpbl*alphak(j,i,1)*betak(j,i,2)
      coef2(j,i,1) = d_one + dtpbl*alphak(j,i,1)*betak(j,i,2)
      coef3(j,i,1) = d_zero
      coefe(j,i,1) = coef1(j,i,1)/coef2(j,i,1)
      coeff1(j,i,1) = udatm(j,i,1)/coef2(j,i,1)
      coeff2(j,i,1) = vdatm(j,i,1)/coef2(j,i,1)
    end do
  end do

  ! top to bottom
  do k = 2 , kz - 1
    do i = istart , iend
      do j = jstart , jend
        coef1(j,i,k) = dtpbl*alphak(j,i,k)*betak(j,i,k+1)
        coef2(j,i,k) = d_one+dtpbl*alphak(j,i,k)*(betak(j,i,k+1)+betak(j,i,k))
        coef3(j,i,k) = dtpbl*alphak(j,i,k)*betak(j,i,k)
        coefe(j,i,k) = coef1(j,i,k)/(coef2(j,i,k)-coef3(j,i,k)*coefe(j,i,k-1))
        coeff1(j,i,k) = (udatm(j,i,k)+coef3(j,i,k)*coeff1(j,i,k-1))/       &
                      (coef2(j,i,k)-coef3(j,i,k)*coefe(j,i,k-1))
        coeff2(j,i,k) = (vdatm(j,i,k)+coef3(j,i,k)*coeff2(j,i,k-1))/       &
                      (coef2(j,i,k)-coef3(j,i,k)*coefe(j,i,k-1))
      end do
    end do
  end do

  ! Nearest to surface
  do i = istart , iend
    do j = jstart , jend
      coef1(j,i,kz) = d_zero
      coef2(j,i,kz) = d_one + dtpbl*alphak(j,i,kz)*betak(j,i,kz)
      coef3(j,i,kz) = dtpbl*alphak(j,i,kz)*betak(j,i,kz)
      drgdot = (uvdrag(j-1,i-1)+uvdrag(j,i-1) + &
                uvdrag(j-1,i)  +uvdrag(j,i))*d_rfour
      uflxsf = drgdot*udatm(j,i,kz)
      vflxsf = drgdot*vdatm(j,i,kz)
      coefe(j,i,kz) = d_zero
      coeff1(j,i,kz) = (udatm(j,i,kz)-dtpbl*alphak(j,i,kz)*uflxsf+      &
                      coef3(j,i,kz)*coeff1(j,i,kz-1))/                  &
                     (coef2(j,i,kz)-coef3(j,i,kz)*coefe(j,i,kz-1))
      coeff2(j,i,kz) = (vdatm(j,i,kz)-dtpbl*alphak(j,i,kz)*vflxsf+      &
                      coef3(j,i,kz)*coeff2(j,i,kz-1))/                  &
                     (coef2(j,i,kz)-coef3(j,i,kz)*coefe(j,i,kz-1))
    end do
  end do
  !
  !   all coefficients have been computed, predict field and put it in
  !   temporary work space tpred
  !
  do i = istart , iend
    do j = jstart , jend
      tpred1(j,i,kz) = coeff1(j,i,kz)
      tpred2(j,i,kz) = coeff2(j,i,kz)
    end do
  end do
 
  do k = kz - 1 , 1 , -1
    do i = istart , iend
      do j = jstart , jend
        tpred1(j,i,k) = coefe(j,i,k)*tpred1(j,i,k+1) + coeff1(j,i,k)
        tpred2(j,i,k) = coefe(j,i,k)*tpred2(j,i,k+1) + coeff2(j,i,k)
      end do
    end do
  end do
  !
  !   calculate tendency due to vertical diffusion using temporary
  !   predicted field
  !
  do k = 1 , kz
    do i = istart , iend
      do j = jstart , jend
        uten(j,i,k) = uten(j,i,k) + &
                      (tpred1(j,i,k)-udatm(j,i,k))*rdtpbl*sfcpd(j,i)
        vten(j,i,k) = vten(j,i,k) + &
                      (tpred2(j,i,k)-vdatm(j,i,k))*rdtpbl*sfcpd(j,i)
      end do
    end do
  end do
  !
  !   Common coefficients.
  !
  do k = 1 , kz
    do i = istart , iend
      do j = jstart , jend
        alphak(j,i,k) = hydf(k)/sfcps(j,i)
      end do
    end do
  end do
  ! 
  !   temperature
  !   calculate coefficients at cross points for temperature
  ! 
  do k = 2 , kz
    do i = istart , iend
      do j = jstart , jend
        betak(j,i,k) = rhohf(j,i,k-1)*kvh(j,i,k)/dza(j,i,k-1)
      end do
    end do
  end do
 
  do i = istart , iend
    do j = jstart , jend
      coef1(j,i,1) = dtpbl*alphak(j,i,1)*betak(j,i,2)
      coef2(j,i,1) = d_one + dtpbl*alphak(j,i,1)*betak(j,i,2)
      coef3(j,i,1) = d_zero
      coefe(j,i,1) = coef1(j,i,1)/coef2(j,i,1)
      coeff1(j,i,1) = thxatm(j,i,1)/coef2(j,i,1)
    end do
  end do

  do k = 2 , kz - 1
    do i = istart , iend
      do j = jstart , jend
        coef1(j,i,k) = dtpbl*alphak(j,i,k)*betak(j,i,k+1)
        coef2(j,i,k) = d_one+dtpbl*alphak(j,i,k)*(betak(j,i,k+1)+betak(j,i,k))
        coef3(j,i,k) = dtpbl*alphak(j,i,k)*betak(j,i,k)
        coefe(j,i,k) = coef1(j,i,k)/(coef2(j,i,k)-coef3(j,i,k)*coefe(j,i,k-1))
        coeff1(j,i,k) = (thxatm(j,i,k)+coef3(j,i,k)*coeff1(j,i,k-1)) / &
                      (coef2(j,i,k)-coef3(j,i,k)*coefe(j,i,k-1))
      end do
    end do
  end do
 
  do i = istart , iend
    do j = jstart , jend
      coef1(j,i,kz) = d_zero
      coef2(j,i,kz) = d_one + dtpbl*alphak(j,i,kz)*betak(j,i,kz)
      coef3(j,i,kz) = dtpbl*alphak(j,i,kz)*betak(j,i,kz)
      coefe(j,i,kz) = d_zero
      coeff1(j,i,kz) = (thxatm(j,i,kz) + dtpbl*alphak(j,i,kz)*hfx(j,i)*rcpd + &
                       coef3(j,i,kz)*coeff1(j,i,kz-1)) / &
                      (coef2(j,i,kz)-coef3(j,i,kz)*coefe(j,i,kz-1))
    end do
  end do
  !
  !   all coefficients have been computed, predict field and put it in
  !   temporary work space tpred
  !
  do i = istart , iend
    do j = jstart , jend
      tpred1(j,i,kz) = coeff1(j,i,kz)
    end do
  end do
 
  do k = kz - 1 , 1 , -1
    do i = istart , iend
      do j = jstart , jend
        tpred1(j,i,k) = coefe(j,i,k)*tpred1(j,i,k+1) + coeff1(j,i,k)
      end do
    end do
  end do
  !
  !   calculate tendency due to vertical diffusion using temporary
  !   predicted field
  !
  do k = 1 , kz
    do i = istart , iend
      do j = jstart , jend
        sf = (tatm(j,i,k)*sfcps(j,i))/thxatm(j,i,k)
        difft(j,i,k) = difft(j,i,k) + &
                       (tpred1(j,i,k)-thxatm(j,i,k))*rdtpbl*sf
      end do
    end do
  end do
  !
  !   water vapor calculate coefficients at cross points for water vapor
  ! 
  do k = 2 , kz
    do i = istart , iend
      do j = jstart , jend
        betak(j,i,k) = rhohf(j,i,k-1)*kvq(j,i,k)/dza(j,i,k-1)
      end do
    end do
  end do
 
  do i = istart , iend
    do j = jstart , jend
      coef1(j,i,1) = dtpbl*alphak(j,i,1)*betak(j,i,2)
      coef2(j,i,1) = d_one + dtpbl*alphak(j,i,1)*betak(j,i,2)
      coef3(j,i,1) = d_zero
      coefe(j,i,1) = coef1(j,i,1)/coef2(j,i,1)
      coeff1(j,i,1) = qvatm(j,i,1)/coef2(j,i,1)
    end do
  end do
 
  do k = 2 , kz - 1
    do i = istart , iend
      do j = jstart , jend
        coef1(j,i,k) = dtpbl*alphak(j,i,k)*betak(j,i,k+1)
        coef2(j,i,k) = d_one+dtpbl*alphak(j,i,k)*(betak(j,i,k+1)+betak(j,i,k))
        coef3(j,i,k) = dtpbl*alphak(j,i,k)*betak(j,i,k)
        coefe(j,i,k) = coef1(j,i,k)/(coef2(j,i,k)-coef3(j,i,k)*coefe(j,i,k-1))
        coeff1(j,i,k) = (qvatm(j,i,k)+coef3(j,i,k)*coeff1(j,i,k-1)) / &
                       (coef2(j,i,k)-coef3(j,i,k)*coefe(j,i,k-1))
      end do
    end do
  end do
 
  do i = istart , iend
    do j = jstart , jend
      coef1(j,i,kz) = d_zero
      coef2(j,i,kz) = d_one + dtpbl*alphak(j,i,kz)*betak(j,i,kz)
      coef3(j,i,kz) = dtpbl*alphak(j,i,kz)*betak(j,i,kz)
      coefe(j,i,kz) = d_zero
      coeff1(j,i,kz) = (qvatm(j,i,kz) + &
               dtpbl*alphak(j,i,kz)*qfx(j,i) + &
               coef3(j,i,kz)*coeff1(j,i,kz-1)) /    &
               (coef2(j,i,kz)-coef3(j,i,kz)*coefe(j,i,kz-1))
    end do
  end do
  !
  !   all coefficients have been computed, predict field and put it in
  !   temporary work space tpred
  ! 
  do i = istart , iend
    do j = jstart , jend
      tpred1(j,i,kz) = coeff1(j,i,kz)
    end do
  end do
 
  do k = kz - 1 , 1 , -1
    do i = istart , iend
      do j = jstart , jend
        tpred1(j,i,k) = coefe(j,i,k)*tpred1(j,i,k+1) + coeff1(j,i,k)
      end do
    end do
  end do
  !
  !   calculate tendency due to vertical diffusion using temporary
  !   predicted field
  !
  if ( ibltyp == 99 ) then
    do k = 1 , kz
      do i = istart , iend
        do j = jstart , jend
          diagqv(j,i,k) = (tpred1(j,i,k)-qvatm(j,i,k))*rdtpbl*sfcps(j,i)
        end do
      end do
    end do
  else
    do k = 1 , kz
      do i = istart , iend
        do j = jstart , jend
          diffq(j,i,k) = diffq(j,i,k) + &
                         (tpred1(j,i,k)-qvatm(j,i,k))*rdtpbl*sfcps(j,i)
        end do
      end do
    end do
  end if
  ! 
  !   calculate coefficients at cross points for cloud vater
  !
  do k = 2 , kz
    do i = istart , iend
      do j = jstart , jend
        betak(j,i,k) = rhohf(j,i,k-1)*kvq(j,i,k)/dza(j,i,k-1)
      end do
    end do
  end do
  do i = istart , iend
    do j = jstart , jend
      coef1(j,i,1) = dtpbl*alphak(j,i,1)*betak(j,i,2)
      coef2(j,i,1) = d_one + dtpbl*alphak(j,i,1)*betak(j,i,2)
      coef3(j,i,1) = d_zero
      coefe(j,i,1) = coef1(j,i,1)/coef2(j,i,1)
      coeff1(j,i,1) = qcatm(j,i,1)/coef2(j,i,1)
    end do
  end do
  do k = 2 , kz - 1
    do i = istart , iend
      do j = jstart , jend
        coef1(j,i,k) = dtpbl*alphak(j,i,k)*betak(j,i,k+1)
        coef2(j,i,k) = d_one+dtpbl*alphak(j,i,k)*(betak(j,i,k+1)+betak(j,i,k))
        coef3(j,i,k) = dtpbl*alphak(j,i,k)*betak(j,i,k)
        coefe(j,i,k) = coef1(j,i,k)/(coef2(j,i,k)-coef3(j,i,k)*coefe(j,i,k-1))
        coeff1(j,i,k) = (qcatm(j,i,k)+coef3(j,i,k)*coeff1(j,i,k-1)) / &
                      (coef2(j,i,k)-coef3(j,i,k)*coefe(j,i,k-1))
      end do
    end do
  end do
  do i = istart , iend
    do j = jstart , jend
      coef1(j,i,kz) = d_zero
      coef2(j,i,kz) = d_one + dtpbl*alphak(j,i,kz)*betak(j,i,kz)
      coef3(j,i,kz) = dtpbl*alphak(j,i,kz)*betak(j,i,kz)
      coefe(j,i,kz) = d_zero
      coeff1(j,i,kz) = (qcatm(j,i,kz)+coef3(j,i,kz)*coeff1(j,i,kz-1)) / &
                     (coef2(j,i,kz)-coef3(j,i,kz)*coefe(j,i,kz-1))
    end do
  end do
  !
  !   all coefficients have been computed, predict field and put it in
  !   temporary work space tpred
  !
  do i = istart , iend
    do j = jstart , jend
      tpred1(j,i,kz) = coeff1(j,i,kz)
    end do
  end do
  do k = kz - 1 , 1 , -1
    do i = istart , iend
      do j = jstart , jend
        tpred1(j,i,k) = coefe(j,i,k)*tpred1(j,i,k+1) + coeff1(j,i,k)
      end do
    end do
  end do
  !
  !   calculate tendency due to vertical diffusion using temporary
  !   predicted field
  !
  if ( ibltyp == 99 ) then
    do k = 1 , kz
      do i = istart , iend
        do j = jstart , jend
          diagqc(j,i,k) = (tpred1(j,i,k)-qcatm(j,i,k))*rdtpbl*sfcps(j,i)
        end do
      end do
    end do
  else
    do k = 1 , kz
      do i = istart , iend
        do j = jstart , jend
          qcten(j,i,k) = qcten(j,i,k) + &
                      (tpred1(j,i,k)-qcatm(j,i,k))*rdtpbl*sfcps(j,i)
        end do
      end do
    end do
  end if
!
! **********************************************************************
!
  !  now add countergradient term to temperature and water vapor equation
  !
  do i = istart , iend
    do j = jstart , jend
      ttnp(j,i,1) = d_zero
    end do
  end do
  do k = 2 , kz
    do i = istart , iend
      do j = jstart , jend
        sf = tatm(j,i,k)/thxatm(j,i,k)
        ttnp(j,i,k) = sf*cpd*rhohf(j,i,k-1)*kvh(j,i,k)*cgh(j,i,k)
      end do
    end do
  end do
  !
  !   compute the tendencies:
  !
  do i = istart , iend
    do j = jstart , jend
      difft(j,i,kz) = difft(j,i,kz) - hydf(kz)*ttnp(j,i,kz)/cpd
    end do
  end do
  do k = 1 , kzm1
    do i = istart , iend
      do j = jstart , jend
        difft(j,i,k) = difft(j,i,k) + hydf(k)*(ttnp(j,i,k+1)-ttnp(j,i,k))/cpd
      end do
    end do
  end do
  if ( lchem .and. lchdrydepo ) then
    !
    !     coef1, coef2, coef3 and coefe are the same as for water vapor
    !     and cloud water so they do not need to be recalculated
    !     recalculation of coef1,2,3  with tracer diffusivity kvc
    ! 
    do k = 2 , kz
      do i = istart , iend
        do j = jstart , jend
          betak(j,i,k) = rhohf(j,i,k-1)*kvc(j,i,k)/dza(j,i,k-1)
        end do
      end do
    end do
    do k = 2 , kz - 1
      do i = istart , iend
        do j = jstart , jend
          coef1(j,i,k) = dtpbl*alphak(j,i,k)*betak(j,i,k+1)
          coef2(j,i,k) = d_one+dtpbl*alphak(j,i,k)*(betak(j,i,k+1)+betak(j,i,k))
          coef3(j,i,k) = dtpbl*alphak(j,i,k)*betak(j,i,k)
        end do
      end do
    end do
    do i = istart , iend
      do j = jstart , jend
        coef1(j,i,1) = dtpbl*alphak(j,i,1)*betak(j,i,2)
        coef2(j,i,1) = d_one + dtpbl*alphak(j,i,1)*betak(j,i,2)
        coef3(j,i,1) = d_zero
        coef1(j,i,kz) = d_zero
        coef2(j,i,kz) = d_one + dtpbl*alphak(j,i,kz)*betak(j,i,kz)
        coef3(j,i,kz) = dtpbl*alphak(j,i,kz)*betak(j,i,kz)
      end do
    end do
    !
    !     set the Vd for in case of prescribed deposition velocities
    !
    do itr = 1 , ntr
      do i = istart , iend
        do j = jstart , jend
          if ( landmsk(j,i) == 0 ) then
            vdep(j,i,itr) = depvel(itr,2)
          else
            vdep(j,i,itr) = depvel(itr,1)
          end if
          !
          !         provisoire test de la routine chdrydep pour les dust
          ! 
!FAB: insert the interactive vdep 
!          if ( chname(itr) == 'DUST' ) vdep(j,i,itr) = d_zero
          vdep(j,i,itr) = d_zero
        end do
      end do
    end do
    do itr = 1 , ntr
      do i = istart , iend
        do j = jstart , jend
          coefe(j,i,1) = coef1(j,i,1)/coef2(j,i,1)
          coeff1(j,i,1) = chmx(j,i,1,itr)/coef2(j,i,1)
        end do
      end do
      do k = 2 , kz - 1
        do i = istart , iend
          do j = jstart , jend
            coefe(j,i,k) = coef1(j,i,k)/(coef2(j,i,k) - &
                           coef3(j,i,k)*coefe(j,i,k-1))
            coeff1(j,i,k) = (chmx(j,i,k,itr)+coef3(j,i,k)*coeff1(j,i,k-1)) / &
                          (coef2(j,i,k)-coef3(j,i,k)*coefe(j,i,k-1))
          end do
        end do
      end do
 
      do i = istart , iend
        do j = jstart , jend
          coefe(j,i,kz) = d_zero
          ! add dry deposition option1
          coeff1(j,i,kz) = (chmx(j,i,kz,itr)-dtpbl*alphak(j,i,kz) * &
                            chmx(j,i,kz,itr)*vdep(j,i,itr)*rhox2d(j,i) + &
                            coef3(j,i,kz)*coeff1(j,i,kz-1)) / &
                           (coef2(j,i,kz)-coef3(j,i,kz)*coefe(j,i,kz-1))
        end do
      end do
      !
      !       all coefficients have been computed, predict field and put
      !       it in temporary work space tpred1
      !
      do i = istart , iend
        do j = jstart , jend
          tpred1(j,i,kz) = coeff1(j,i,kz)
        end do
      end do
      do k = kz - 1 , 1 , -1
        do i = istart , iend
          do j = jstart , jend
            tpred1(j,i,k) = coefe(j,i,k)*tpred1(j,i,k+1) + coeff1(j,i,k)
          end do
        end do
      end do
      !
      !       calculate tendency due to vertical diffusion using temporary
      !       predicted field
      !       Dry deposition option 1 is included
      !
      do k = 1 , kz
        do i = istart , iend
          do j = jstart , jend
            chten(j,i,k,itr) = chten(j,i,k,itr) +  &
                        (tpred1(j,i,k)-chmx(j,i,k,itr))*rdtpbl*sfcps(j,i)
          end do
        end do
      end do
      do i = istart , iend
        do j = jstart , jend
!          if ( chname(itr) /= 'DUST' ) &
! shut down thid diag temporarly
!            drmr(j,i,itr) = drmr(j,i,itr) + chmx(j,i,kz,itr)* &
!                vdep(j,i,itr)*sfcps(j,i)*dtpbl*d_half*rhox2d(j,i)* &
!                hydf(kz)/sfcps(j,i)
 
        end do
      end do
    end do
  end if
  call time_end(subroutine_name,idindx)
 
  end subroutine holtbl
!
! ------------------------------------------------------------
! this routine computes the boundary layer eddy diffusivities
! for momentum, heat and moisture and the counter-gradient
! terms for heat and moisture.
!
! reference : holtslag, de bruijn and pan - mwr - 8/90
!
! input arguments :  j       longitudinal position index
!                    ubx3d   u wind component
!                    vbx3d   v wind component
!                    thx3d   potential temperature
!                    thvx    virtual potential temperature
!                    za      height of half sigma levels
!                    f       coriolis parameter
!                    shum    specific humidity
!                    xhfx    sensible heat flux
!                    xqfx    sfc kinematic moisture flux
!                    th10    virt. pot. temp. at 10m
!                    hfxv    surface virtual heat flux
!                    obklen  monin obukov length
!                    ustr    friction velocity
!
! input/output
! arguments :        therm   thermal temperature excess
!
! output arguments : cgh     counter-gradient term for heat
!                    cgq     counter-gradient term for moisture
!                    kvm     eddy diffusivity for momentum
!                    kvh     eddy diffusivity for heat
!                    kvq     eddy diffusivity for moisture
!                   zpbl     boundary layer height
! ------------------------------------------------------------
!
  subroutine blhnew(jstart,jend,istart,iend)
  implicit none
  integer , intent(in) :: jstart , jend , istart , iend
!
  real(dp) :: fak1 , fak2 , fht , xfmt , pblk , pblk1 , pblk2 , &
             phpblm , pr , therm2 , tkv , tlv , wsc , z , zh , &
             zl , zm , zp , zzh , zzhnew , zzhnew2
  integer :: i , j , k , k2
  character (len=64) :: subroutine_name='blhnew'
  integer :: idindx=0
!
  call time_begin(subroutine_name,idindx)
!
  ! note: kmxpbl, max no. of pbl levels, calculated in param
  ! compute richardson number
  do k = kz , kmxpbl , -1
    do i = istart , iend
      do j = jstart , jend
        ri(j,i,k) = egrav*(thvx(j,i,k)-th10(j,i))*za(j,i,k) / &
                          (th10(j,i)*vv(j,i,k))
      end do
    end do
  end do
  ! first, set bl height to height of lowest model level
  do i = istart , iend
    do j = jstart , jend
      zpbl(j,i) = za(j,i,kz)
      kpbl(j,i) = kz
    end do
  end do
  ! looking for bl top
  do k = kz , kmxpbl + 1 , -1
    k2 = k - 1
    do i = istart , iend
      do j = jstart , jend
        ! bl height lies between this level and the last
        ! use linear interp. of rich. no. to height of ri=ricr
        if ( (ri(j,i,k) < ricr) .and. (ri(j,i,k2) >= ricr) ) then
          zpbl(j,i) = za(j,i,k) + (za(j,i,k2)-za(j,i,k)) &
              *((ricr-ri(j,i,k))/(ri(j,i,k2)-ri(j,i,k)))
          kpbl(j,i) = k
        end if
      end do
    end do
  end do
  do i = istart , iend
    do j = jstart , jend
      ! set bl top to highest allowable model layer
      if ( ri(j,i,kmxpbl) < ricr ) then
        zpbl(j,i) = za(j,i,kmxpbl)
        kpbl(j,i) = kmxpbl
      end if
    end do
  end do
  ! recompute richardson no. at lowest model level
  do i = istart , iend
    do j = jstart , jend
      if ( hfxv(j,i) > d_zero ) then
        ! estimate of convective velocity scale
        xfmt = (d_one-(binm*zpbl(j,i)/obklen(j,i)))**onet
        wsc = ustr(j,i)*xfmt
        ! thermal temperature excess
        therm(j,i) = (xhfx(j,i)+mult*thxatm(j,i,kz)*xqfx(j,i))*fak/wsc
        ri(j,i,kz) = -egrav*therm(j,i)*za(j,i,kz)/(th10(j,i)*vv(j,i,kz))
      else
        therm(j,i) = d_zero
      end if
    end do
  end do
  ! recompute richardson no. at other model levels
  do k = kz - 1 , kmxpbl , -1
    do i = istart , iend
      do j = jstart , jend
        if ( hfxv(j,i) > d_zero ) then
          tlv = th10(j,i) + therm(j,i)
          tkv = thxatm(j,i,k)*(d_one+mult*(qvatm(j,i,k)/(qvatm(j,i,k)+d_one)))
          ri(j,i,k) = egrav*(tkv-tlv)*za(j,i,k)/(th10(j,i)*vv(j,i,k))
        end if
      end do
    end do
  end do
  ! improve estimate of bl height under convective conditions
  ! using convective temperature excess (therm)
  do k = kz , kmxpbl + 1 , -1
    k2 = k - 1
    do i = istart , iend
      do j = jstart , jend
        if ( hfxv(j,i) > d_zero ) then
          ! bl height lies between this level and the last
          ! use linear interp. of rich. no. to height of ri=ricr
          if ( (ri(j,i,k) < ricr) .and. (ri(j,i,k2) >= ricr) ) then
            zpbl(j,i) = za(j,i,k) + (za(j,i,k2)-za(j,i,k)) * &
                        ((ricr-ri(j,i,k))/(ri(j,i,k2)-ri(j,i,k)))
            kpbl(j,i) = k
          end if
        end if
      end do
    end do
  end do
  do i = istart , iend
    do j = jstart , jend
      if ( hfxv(j,i) > d_zero ) then
        ! set bl top to highest allowable model layer
        if ( ri(j,i,kmxpbl) < ricr ) then
          zpbl(j,i) = za(j,i,kmxpbl)
        end if
      end if
    end do
  end do
  ! limit bl height to be at least mech. mixing depth
  do i = istart , iend
    do j = jstart , jend
      ! compute mechanical mixing depth, set to lowest model level if lower
      phpblm = 0.07D0*ustr(j,i)/pfcor(j,i)
      phpblm = dmax1(phpblm,za(j,i,kz))
      zpbl(j,i) = dmax1(zpbl(j,i),phpblm)
      kpbl(j,i) = kz
    end do
  end do
 
  do k = kz , kmxpbl + 1 , -1
    k2 = k - 1
    do i = istart , iend
      do j = jstart , jend
        pblk = d_zero
        zm = za(j,i,k)
        zp = za(j,i,k2)
        if ( zm < zpbl(j,i) ) then
          zp = dmin1(zp,zpbl(j,i))
          z = (zm+zp)*d_half
          zh = z/zpbl(j,i)
          zl = z/obklen(j,i)
          if ( zh <= d_one ) then
            zzh = d_one - zh
            zzh = zzh**pink
!xexp4      zzhnew = zpbl(j,i)*(d_one-zh)*zh**1.5
!xexp5      zzhnew = 0.5*zpbl(j,i)*(d_one-zh)*zh**1.5
!xexp6      zzhnew = d_one - zh
!xexp7      zzhnew =0.5* (d_one - zh)
!Sara
!           zzhnew =0.25* (d_one - zh)
!           zzhnew =0.75* (d_one - zh)
!Sara_
!xexp10     zzhnew =zh * (d_one - zh)**2
            zzhnew = (d_one-zh)*d_rfour
            zzhnew2 = (d_one-zh)**d_two
          else
            zzh = d_zero
            zzhnew = d_zero
            zzhnew2 = d_zero
          end if
          fak1 = ustr(j,i)*zpbl(j,i)*vonkar
          if ( hfxv(j,i) <= d_zero ) then
            ! stable and neutral conditions
            ! igroup = 1
            ! prevent pblk from becoming too small in very stable
            ! conditions
            if ( zl <= d_one ) then
              pblk = fak1*zh*zzh/(d_one+betas*zl)
!xexp5        pblk1 = vonkar * ustr(j,i) / (d_one+betas*zl) * zzhnew
              pblk1 = fak1*zh*zzhnew/(d_one+betas*zl)
              pblk2 = fak1*zh*zzhnew2/(d_one+betas*zl)
            else
              pblk = fak1*zh*zzh/(betas+zl)
!xexp5        pblk1 = vonkar * ustr(j,i) / (betas+zl) * zzhnew
              pblk1 = fak1*zh*zzhnew/(betas+zl)
              pblk2 = fak1*zh*zzhnew2/(betas+zl)
            end if
            ! compute eddy diffusivities
            kvm(j,i,k) = dmax1(pblk,kzo)
            kvh(j,i,k) = kvm(j,i,k)
            kvq(j,i,k) = dmax1(pblk1,kzo)
            ! Erika put k=0 in very stable conditions
            if ( zl <= 0.1D0 ) then
              kvm(j,i,k) = d_zero
              kvh(j,i,k) = kvm(j,i,k)*d_zero
              kvq(j,i,k) = d_zero
            end if
            ! Erika put k=0 in very stable conditions
            if ( lchem ) then
              kvc(j,i,k) = dmax1(pblk2,kzo)
            end if
            ! compute counter-gradient term
            cgh(j,i,k) = d_zero
          else
            ! unstable conditions
            ! compute counter gradient term
            if ( zh >= sffrac ) then
              ! igroup = 2
              xfmt = (d_one-binm*zpbl(j,i)/obklen(j,i))**onet
              fht = dsqrt(d_one-binh*zpbl(j,i)/obklen(j,i))
              wsc = ustr(j,i)*xfmt
              pr = (xfmt/fht) + ccon
              fak2 = wsc*zpbl(j,i)*vonkar
              pblk = fak2*zh*zzh
!xexp5        pblk1 = vonkar * wsc * zzhnew
              pblk1 = fak2*zh*zzhnew
              pblk2 = fak2*zh*zzhnew2
              therm2 = fak/(zpbl(j,i)*wsc)
              cgh(j,i,k) = hfxv(j,i)*therm2
            else
              ! igroup = 3
              pblk = fak1*zh*zzh*(d_one-betam*zl)**onet
!xexp5        pblk1 = vonkar * ustr(j,i) * zzhnew * (d_one-betam*zl)**onet
              pblk1 = fak1*zh*zzhnew*(d_one-betam*zl)**onet
              pblk2 = fak1*zh*zzhnew2*(d_one-betam*zl)**onet
              pr = ((d_one-betam*zl)**onet)/dsqrt(d_one-betah*zl)
              cgh(j,i,k) = d_zero
            end if
            ! compute eddy diffusivities
            kvm(j,i,k) = dmax1(pblk,kzo)
            kvh(j,i,k) = dmax1((pblk/pr),kzo)
            kvq(j,i,k) = dmax1(pblk1,kzo)
            if ( lchem ) then
              kvc(j,i,k) = dmax1(pblk2,kzo)
            end if
          end if
        end if
      end do
    end do
  end do
 
  if ( ibltyp == 99 ) then
    do i = istart , iend
      do j = jstart , jend
        do k = 1 , kz
          kvm(j,i,k) = uwstateb%kzm(j,i,k)
          kvh(j,i,k) = uwstateb%kth(j,i,k)
          kvq(j,i,k) = uwstateb%kth(j,i,k)
          cgh(j,i,k) = d_zero
        end do
      end do
    end do
  end if
  call time_end(subroutine_name,idindx)
!
  end subroutine blhnew
!
end module mod_pbl_holtbl

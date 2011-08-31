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
  use mod_message
  use mod_pbl_common
!
  private
!
  public :: allocate_mod_pbl_holtbl , holtbl
!
  real(8) , pointer , dimension(:,:,:) :: cgh , kvc , kvh ,   &
                                          kvm , kvq ! , cgq
  real(8) , pointer, dimension(:,:) :: hfxv , obklen , th10 , &
                                       ustr , xhfx , xqfx
!
  real(8) , pointer , dimension(:,:) :: alphak , betak , &
                        coef1 , coef2 , coef3 , coefe , coeff1 , &
                        coeff2 , tpred1 , tpred2
  real(8) , pointer , dimension(:,:) :: kzm , rc , ttnp
  real(8) , pointer , dimension(:,:) :: vdep
  real(8) , pointer , dimension(:) :: govrth
!
  real(8) , pointer , dimension(:,:,:) :: dza
  real(8) , pointer , dimension(:,:,:) :: akzz1 , akzz2
  real(8) , pointer , dimension(:) :: wkrecv , wksend
  real(8) , pointer , dimension(:,:,:) :: rhohf
!
  real(8) , pointer , dimension(:,:) :: ri
  real(8) , pointer , dimension(:) :: therm
!
!     minimum eddy diffusivity
  real(8) , parameter :: kzo = d_one
  real(8) , parameter :: szkm = 1600.0D0
!     coef. of proportionality and lower % of bl in sfc layer
  real(8) , parameter :: fak = 8.5D0
  real(8) , parameter :: sffrac = 0.1D0
!     beta coefs. for momentum, stable conditions and heat
  real(8) , parameter :: betam = 15.0D0
  real(8) , parameter :: betas = 5.0D0
  real(8) , parameter :: betah = 15.0D0
  real(8) , parameter :: mult = 0.61D0
  real(8) , parameter :: ccon = fak*sffrac*vonkar
  real(8) , parameter :: binm = betam*sffrac
  real(8) , parameter :: binh = betah*sffrac
!     power in formula for k and critical ri for judging stability
  real(8) , parameter :: pink = d_two
  real(8) , parameter :: ricr = d_rfour
! 
  contains
!
  subroutine allocate_mod_pbl_holtbl(ichem,ichdrdepo)
    implicit none
    integer , intent(in) :: ichem , ichdrdepo

    call getmem2d(alphak,1,iy,1,kz,'mod_holtbl:alphak')
    call getmem2d(betak,1,iy,1,kz,'mod_holtbl:betak')
    call getmem2d(coef1,1,iy,1,kz,'mod_holtbl:coef1')
    call getmem2d(coef2,1,iy,1,kz,'mod_holtbl:coef2')
    call getmem2d(coef3,1,iy,1,kz,'mod_holtbl:coef3')
    call getmem2d(coefe,1,iy,1,kz,'mod_holtbl:coefe')
    call getmem2d(coeff1,1,iy,1,kz,'mod_holtbl:coeff1')
    call getmem2d(coeff2,1,iy,1,kz,'mod_holtbl:coeff2')
    call getmem2d(tpred1,1,iy,1,kz,'mod_holtbl:tpred1')
    call getmem2d(tpred2,1,iy,1,kz,'mod_holtbl:tpred2')
    call getmem2d(ri,1,iy,1,kz,'mod_holtbl:ri')
    call getmem2d(kzm,1,iym1,1,kz,'mod_holtbl:kzm')
    call getmem2d(rc,1,iym1,1,kz,'mod_holtbl:rc')
    call getmem2d(ttnp,1,iym1,1,kz,'mod_holtbl:ttnp')
    call getmem1d(govrth,1,iym1,'mod_holtbl:govrth')
    call getmem1d(therm,1,iy,'mod_holtbl:therm')
    call getmem3d(cgh,1,iy,1,kz,1,jxp,'mod_holtbl:cgh')
!    call getmem3d(cgq,1,iy,1,kz,1,jxp,'mod_holtbl:cgq')
    call getmem3d(kvh,1,iy,1,kz,1,jxp,'mod_holtbl:kvh')
    call getmem3d(kvm,1,iy,1,kz,1,jxp,'mod_holtbl:kvm')
    call getmem3d(kvq,1,iy,1,kz,1,jxp,'mod_holtbl:kvq')
    call getmem2d(hfxv,1,iy,1,jxp,'mod_holtbl:hfxv')
    call getmem2d(obklen,1,iy,1,jxp,'mod_holtbl:obklen')
    call getmem2d(th10,1,iy,1,jxp,'mod_holtbl:th10')
    call getmem2d(ustr,1,iy,1,jxp,'mod_holtbl:ustr')
    call getmem2d(xhfx,1,iy,1,jxp,'mod_holtbl:xhfx')
    call getmem2d(xqfx,1,iy,1,jxp,'mod_holtbl:xqfx')
    call getmem3d(dza,1,iym1,1,kz,1,jxp,'mod_holtbl:dza')
    call getmem3d(akzz1,1,iym1,1,kz,0,jxp+1,'mod_holtbl:akzz1')
    call getmem3d(akzz2,1,iym1,1,kz,0,jxp+1,'mod_holtbl:akzz2')
    call getmem1d(wkrecv,1,2*iym2*kz,'mod_holtbl:wkrecv')
    call getmem1d(wksend,1,2*iym2*kz,'mod_holtbl:wksend')
    call getmem3d(rhohf,1,iy,1,kz,1,jxp,'mod_holtbl:rhohf')
    if ( ichem == 1 ) then
      if ( ichdrdepo == 1 ) then
        call getmem2d(vdep,1,iym1,1,ntr,'mod_holtbl:vdep')
      end if
      call getmem3d(kvc,1,iy,1,kz,1,jxp,'mod_holtbl:kvc')
    end if
  end  subroutine allocate_mod_pbl_holtbl
!
  subroutine holtbl
!
#ifndef IBM
  use mpi
#else 
  include 'mpif.h'
#endif 
!
  implicit none
!
  real(8) :: drgdot , kzmax , oblen , xps , ps2 , ri , &
             sf , sh10 , ss , uflxsf , uflxsfx , vflxsf ,     &
             vflxsfx
  integer :: jdx , jm1
#ifndef BAND
  integer :: jdxm1
#endif
  integer :: i , idx , idxm1 , itr , j , k
  integer :: ierr , ii
!
!----------------------------------------------------------------------
!-----some of the storage spaces for high-resolution pbl
!     are used store the variables in this subroutine.
!     difft(i,k,j)   : temperature tendency (tten)
!     diffq(i,k,j)   : water vapor tendency (qvten)
!
!-----decouple flux-form variables to give u,v,t,theta,theta-vir.,
!     t-vir., qv, and qc at cross points and at ktau-1.
!
!     *** note ***
!     the boundary winds may not be adequately affected by friction,
!     so use only interior values of ubx3d and vbx3d to calculate
!     tendencies.
!
  call mpi_sendrecv(uvdrag(1,jxp),iy,mpi_real8,ieast,1,  &
                    uvdrag(1,0),  iy,mpi_real8,iwest,1,  &
                    mycomm,mpi_status_ignore,ierr)
!
  do j = jbegin , jendx
    do k = 1 , kz
      do i = 2 , iym1
        thvx(i,k,j) = thxatm(i,k,j)*(d_one+ep1*qvatm(i,k,j))
      end do
    end do
!
!.....density at surface is stored in rhox2d(i,j), at half levels in
!       rhohf(i,k,j).
!
    do k = 1 , kzm1
      do i = 2 , iym1
        dza(i,k,j) = za(i,k,j) - za(i,k+1,j)
        xps = (hlev(k)*sfcps(i,j)+ptp)*d_1000
        ps2 = (hlev(k+1)*sfcps(i,j)+ptp)*d_1000
        rhohf(i,k,j) = (ps2-xps)/(egrav*dza(i,k,j))
      end do
    end do
!
    do i = 2 , iym1
      govrth(i) = egrav/thxatm(i,kz,j)
    end do
!
! *********************************************************************
!
!     compute the vertical diffusion term:
!
    do k = 2 , kz
      do i = 2 , iym1
        rc(i,k) = 0.257D0*dzq(i,k,j)**0.175D0
      end do
    end do
!
!-----compute the diffusion coefficient:
!
!       blackadar scheme above boundary layer top
!
    do k = 2 , kz
      do i = 2 , iym1
        kzmax = 0.8D0*dza(i,k-1,j)*dzq(i,k,j)/dtpbl
        ss = ((uatm(i,k-1,j)-uatm(i,k,j))*   &
              (uatm(i,k-1,j)-uatm(i,k,j))+   &
              (vatm(i,k-1,j)-vatm(i,k,j))*   &
              (vatm(i,k-1,j)-vatm(i,k,j)))/  &
              (dza(i,k-1,j)*dza(i,k-1,j)) + 1.0D-9
        ri = govrth(i)*(thvx(i,k-1,j)-thvx(i,k,j))/(ss*dza(i,k-1,j))
        if ( (ri-rc(i,k)) >= d_zero ) then
          kzm(i,k) = kzo
        else
          kzm(i,k) = kzo + dsqrt(ss)*(rc(i,k)-ri)*szkm/rc(i,k)
        end if
        kzm(i,k) = dmin1(kzm(i,k),kzmax)
      end do
    end do
!
! *********************************************************************
!
!       holtslag pbl
!
!       initialize bl diffusion coefficients and counter-gradient terms
!       with free atmosphere values and make specific humidity
!
    do k = 2 , kz
      do i = 2 , iym1
!           eddy diffusivities for momentum, heat and moisture
        kvm(i,k,j) = kzm(i,k)
        kvh(i,k,j) = kzm(i,k)
        kvq(i,k,j) = kzm(i,k)
!chem
        if ( lchem ) then
          kvc(i,k,j) = kzm(i,k)
        end if
!chem_
!           counter gradient terms for heat and moisture
        cgh(i,k,j) = d_zero
      end do
    end do
 
    do i = 2 , iym1
!         compute friction velocity
      idx = i
      jdx = j
#ifndef BAND
      idx = min0(idx,iym1)
      idxm1 = i - 1
      idxm1 = max0(idxm1,2)
      jdxm1 = j - 1
      if ( myid == nproc-1 ) jdx = min0(jdx,jendx)
      if ( myid == 0 ) jdxm1 = max0(jdxm1,2)
#endif
      uflxsfx = uvdrag(idx,jdx)*uatm(i,kz,j)
      vflxsfx = uvdrag(idx,jdx)*vatm(i,kz,j)

      ustr(i,j) = dsqrt(dsqrt(uflxsfx*uflxsfx+vflxsfx*vflxsfx) / &
                         rhox2d(i,j))
 
!         convert surface fluxes to kinematic units
      xhfx(i,j) = hfx(i,j)/(cpd*rhox2d(i,j))
      xqfx(i,j) = qfx(i,j)/rhox2d(i,j)
!         compute virtual heat flux at surface
      hfxv(i,j) = xhfx(i,j) + mult*thxatm(i,kz,j)*xqfx(i,j)
    end do
!
!       estimate potential temperature at 10m via log temperature
!       profile in the surface layer (brutsaert, p. 63).
!       calculate mixing ratio at 10m by assuming a constant
!       value from the surface to the lowest model level.
!
 
    do i = 2 , iym1
      sh10 = qvatm(i,kz,j)/(qvatm(i,kz,j)+d_one)
!     th10(i,j) = ((thxatm(i,kz,j)+tg(i,j))*d_half)*(d_one+mult*sh10)
!     th10(i,j) = thvx(i,kz,j) + hfxv(i,j)/(vonkar*ustr(i,j)* &
!                 dlog(za(i,kz,j)*d_r10)
 
!     "virtual" potential temperature
      if ( hfxv(i,j) >= d_zero ) then
        th10(i,j) = thvx(i,kz,j)
      else
!       th10(i,j) =
!----    (0.25*thxatm(i,kz,j)+0.75*tg(i,j))*(d_one+mult*sh10) first
!       approximation for obhukov length
        oblen = -d_half*(thxatm(i,kz,j)+tg(i,j)) *  &
                (d_one+mult*sh10)*ustr(i,j)**d_three /  &
                (egrav*vonkar*(hfxv(i,j)+dsign(1.0D-10,hfxv(i,j))))
        if ( oblen >= za(i,kz,j) ) then
          th10(i,j) = thvx(i,kz,j) + hfxv(i,j)/(vonkar*ustr(i,j))*  &
             (dlog(za(i,kz,j)*d_r10)+d_five/oblen*(za(i,kz,j)-d_10))
        else if ( oblen < za(i,kz,j) .and. oblen > d_10 ) then
          th10(i,j) = thvx(i,kz,j) + hfxv(i,j)/(vonkar*ustr(i,j))*  &
              (dlog(oblen*d_r10)+d_five/oblen*(oblen-d_10)+         &
              6.0D0*dlog(za(i,kz,j)/oblen))
        else if ( oblen <= d_10 ) then
          th10(i,j) = thvx(i,kz,j) + hfxv(i,j)/(vonkar*ustr(i,j))   &
                      *6.0D0*dlog(za(i,kz,j)*d_r10)
        end if
        th10(i,j) = dmax1(th10(i,j),tg(i,j))
      end if
!gtb      th10(i,j) = dmin1(th10(i,j),tg(i,j))  ! gtb add to minimize
 
!         obklen compute obukhov length
      obklen(i,j) = -th10(i,j)*ustr(i,j)**d_three / &
              (egrav*vonkar*(hfxv(i,j)+dsign(1.0D-10,hfxv(i,j))))
    end do
!
!       compute diffusivities and counter gradient terms
!
  end do

  call blhnew

  do j = jbegin , jendx
#ifndef BAND
  if ( (myid /= nproc-1) .or. (myid == nproc-1 .and. j < jendx)) then
#endif
    do k = 1 , kz
      do i = 2 , iym1
        if ( k > 1 ) then
          akzz1(i,k,j) = rhohf(i,k-1,j)*kvm(i,k,j)/dza(i,k-1,j)
        end if
        akzz2(i,k,j) = egrav/(sfcps(i,j)*d_1000)/dlev(k)
      end do
    end do
#ifndef BAND
  end if
#endif
  end do
  ii = 0
  do k = 1 , kz
    do i = 2 , iym1
      ii = ii + 1
      wksend(ii) = akzz1(i,k,jxp)
    end do
  end do
  do k = 1 , kz
    do i = 2 , iym1
      ii = ii + 1
      wksend(ii) = akzz2(i,k,jxp)
    end do
  end do
  call mpi_sendrecv(wksend,iym2*kz*2,mpi_real8,ieast,1, &
                    wkrecv,iym2*kz*2,mpi_real8,iwest,1, &
                    mycomm,mpi_status_ignore,ierr)
  ii = 0
  do k = 1 , kz
    do i = 2 , iym1
      ii = ii + 1
      akzz1(i,k,0) = wkrecv(ii)
    end do
  end do
  do k = 1 , kz
    do i = 2 , iym1
      ii = ii + 1
      akzz2(i,k,0) = wkrecv(ii)
    end do
  end do
!
  do j = jbegin , jendx
     jm1 = j-1
 
!       calculate coefficients at dot points for u and v wind
 
#ifndef BAND
    if ( myid == 0 .and. j == 2 ) then
      do k = 1 , kz
        do i = 2 , iym1
          idx = i
          idx = min0(idx,iym2)
          idxm1 = i - 1
          idxm1 = max0(idxm1,2)
          if ( k > 1 ) then
            betak(i,k) = d_half*(akzz1(idx,k,j)+akzz1(idxm1,k,j))
          end if
          alphak(i,k) = d_half*(akzz2(idx,k,j)+akzz2(idxm1,k,j))
        end do
      end do
    else if ( myid == nproc-1 .and. j == jendx ) then
      do k = 1 , kz
        do i = 2 , iym1
          idx = i
          idx = min0(idx,iym2)
          idxm1 = i - 1
          idxm1 = max0(idxm1,2)
          if ( k > 1 ) then
            betak(i,k) = d_half*(akzz1(idx,k,jm1)+akzz1(idxm1,k,jm1))
          end if
          alphak(i,k) = d_half*(akzz2(idx,k,jm1)+akzz2(idxm1,k,jm1))
        end do
      end do
    else
#endif
     do k = 1 , kz
       do i = 2 , iym1
         idx = i
         idx = min0(idx,iym2)
         idxm1 = i - 1
         idxm1 = max0(idxm1,2)
         if ( k > 1 ) then
           betak(i,k) = (akzz1(idx,k,jm1)+akzz1(idxm1,k,jm1)+ &
                         akzz1(idx,k,j)+akzz1(idxm1,k,j))*d_rfour
         end if
         alphak(i,k) = (akzz2(idx,k,jm1)+akzz2(idxm1,k,jm1)+ &
                        akzz2(idx,k,j)+akzz2(idxm1,k,j))*d_rfour
       end do
     end do
#ifndef BAND
    end if
#endif
!
! **********************************************************************
!
!       start now procedure for implicit diffusion calculations
!       performed separately for wind (dot points)
!       and temperature and water vapor (cross points)
!       countergradient term is not included in the implicit diffusion
!       scheme its effect is included as in the old explicit scheme
!       calculations assume fluxes positive upward, so the sign in front
!       of uflxsf and vflxsf has been changed in the various terms
!
!       wind components
!
!       first compute coefficients of the tridiagonal matrix
!
! **********************************************************************
!
    do k = 2 , kz - 1
      do i = 2 , iym1
        coef1(i,k) = dtpbl*alphak(i,k)*betak(i,k+1)
        coef2(i,k) = d_one+dtpbl*alphak(i,k)*(betak(i,k+1)+betak(i,k))
        coef3(i,k) = dtpbl*alphak(i,k)*betak(i,k)
      end do
    end do
 
    do i = 2 , iym1
      coef1(i,1) = dtpbl*alphak(i,1)*betak(i,2)
      coef2(i,1) = d_one + dtpbl*alphak(i,1)*betak(i,2)
      coef3(i,1) = d_zero
      coef1(i,kz) = d_zero
      coef2(i,kz) = d_one + dtpbl*alphak(i,kz)*betak(i,kz)
      coef3(i,kz) = dtpbl*alphak(i,kz)*betak(i,kz)
    end do
 
    do i = 2 , iym1
      coefe(i,1) = coef1(i,1)/coef2(i,1)
      coeff1(i,1) = udatm(i,1,j)/coef2(i,1)
      coeff2(i,1) = vdatm(i,1,j)/coef2(i,1)
    end do
 
    do k = 2 , kz - 1
      do i = 2 , iym1
        coefe(i,k) = coef1(i,k)/(coef2(i,k)-coef3(i,k)*coefe(i,k-1))
        coeff1(i,k) = (udatm(i,k,j)+coef3(i,k)*coeff1(i,k-1))/       &
                      (coef2(i,k)-coef3(i,k)*coefe(i,k-1))
        coeff2(i,k) = (vdatm(i,k,j)+coef3(i,k)*coeff2(i,k-1))/       &
                      (coef2(i,k)-coef3(i,k)*coefe(i,k-1))
      end do
    end do
 
    do i = 2 , iym1
      idx = i
      idx = min0(idx,iym1)
      idxm1 = i - 1
      idxm1 = max0(idxm1,2)

#ifdef BAND
      drgdot = (uvdrag(idxm1,jm1)+uvdrag(idxm1,j) + &
                uvdrag(idx,jm1)  +uvdrag(idx,j))*d_rfour
#else
      jdx = j
      jdxm1 = j - 1
      if ( myid == nproc-1 ) jdx = min0(jdx,jendx)
      if ( myid == 0 ) jdxm1 = max0(jdxm1,2)
      drgdot = (uvdrag(idxm1,jdxm1)+uvdrag(idxm1,jdx)+  &
              uvdrag(idx,jdxm1)+uvdrag(idx,jdx))*d_rfour
#endif
      uflxsf = drgdot*udatm(i,kz,j)
      vflxsf = drgdot*vdatm(i,kz,j)
 
      coefe(i,kz) = d_zero
      coeff1(i,kz) = (udatm(i,kz,j)-dtpbl*alphak(i,kz)*uflxsf+          &
                      coef3(i,kz)*coeff1(i,kz-1))/                  &
                     (coef2(i,kz)-coef3(i,kz)*coefe(i,kz-1))
      coeff2(i,kz) = (vdatm(i,kz,j)-dtpbl*alphak(i,kz)*vflxsf+          &
                      coef3(i,kz)*coeff2(i,kz-1))/                  &
                     (coef2(i,kz)-coef3(i,kz)*coefe(i,kz-1))
    end do
!
!       all coefficients have been computed, predict field and put it in
!       temporary work space tpred
!
    do i = 2 , iym1
      tpred1(i,kz) = coeff1(i,kz)
      tpred2(i,kz) = coeff2(i,kz)
    end do
 
    do k = kz - 1 , 1 , -1
      do i = 2 , iym1
        tpred1(i,k) = coefe(i,k)*tpred1(i,k+1) + coeff1(i,k)
        tpred2(i,k) = coefe(i,k)*tpred2(i,k+1) + coeff2(i,k)
      end do
    end do
 
!
!       calculate tendency due to vertical diffusion using temporary
!       predicted field
!
    do k = 1 , kz
      do i = 2 , iym1
        uten(i,k,j) = uten(i,k,j) + &
                        (tpred1(i,k)-udatm(i,k,j))/dtpbl*sfcpd(i,j)
        vten(i,k,j) = vten(i,k,j) + &
                        (tpred2(i,k)-vdatm(i,k,j))/dtpbl*sfcpd(i,j)
      end do
    end do
 
!       temperature
!
!       calculate coefficients at cross points for temperature
 
    do k = 1 , kz
      do i = 2 , iym1
        if ( k > 1 ) then
          betak(i,k) = rhohf(i,k-1,j)*kvh(i,k,j)/dza(i,k-1,j)
        end if
        alphak(i,k) = egrav/(sfcps(i,j)*d_1000)/dlev(k)
      end do
    end do
 
    do k = 2 , kz - 1
      do i = 2 , iym1
        coef1(i,k) = dtpbl*alphak(i,k)*betak(i,k+1)
        coef2(i,k) = d_one+dtpbl*alphak(i,k)*(betak(i,k+1)+betak(i,k))
        coef3(i,k) = dtpbl*alphak(i,k)*betak(i,k)
      end do
    end do
 
    do i = 2 , iym1
      coef1(i,1) = dtpbl*alphak(i,1)*betak(i,2)
      coef2(i,1) = d_one + dtpbl*alphak(i,1)*betak(i,2)
      coef3(i,1) = d_zero
      coef1(i,kz) = d_zero
      coef2(i,kz) = d_one + dtpbl*alphak(i,kz)*betak(i,kz)
      coef3(i,kz) = dtpbl*alphak(i,kz)*betak(i,kz)
    end do
 
    do i = 2 , iym1
      coefe(i,1) = coef1(i,1)/coef2(i,1)
      coeff1(i,1) = thxatm(i,1,j)/coef2(i,1)
    end do
 
    do k = 2 , kz - 1
      do i = 2 , iym1
        coefe(i,k) = coef1(i,k)/(coef2(i,k)-coef3(i,k)*coefe(i,k-1))
        coeff1(i,k) = (thxatm(i,k,j)+coef3(i,k)*coeff1(i,k-1)) / &
                      (coef2(i,k)-coef3(i,k)*coefe(i,k-1))
      end do
    end do
 
    do i = 2 , iym1
      coefe(i,kz) = d_zero
      coeff1(i,kz) = (thxatm(i,kz,j) + &
             dtpbl*alphak(i,kz)*hfx(i,j)*rcpd + &
               coef3(i,kz)*coeff1(i,kz-1)) /       &
               (coef2(i,kz)-coef3(i,kz)*coefe(i,kz-1))
    end do
 
!
!       all coefficients have been computed, predict field and put it in
!       temporary work space tpred
!
 
    do i = 2 , iym1
      tpred1(i,kz) = coeff1(i,kz)
    end do
 
    do k = kz - 1 , 1 , -1
      do i = 2 , iym1
        tpred1(i,k) = coefe(i,k)*tpred1(i,k+1) + coeff1(i,k)
      end do
    end do
 
!
!       calculate tendency due to vertical diffusion using temporary
!       predicted field
!
    do k = 1 , kz
      do i = 2 , iym1
        sf = tatm(i,k,j)/thxatm(i,k,j)
        difft(i,k,j) = difft(i,k,j) + &
                       (tpred1(i,k)-thxatm(i,k,j))/dtpbl*sf
      end do
    end do
!
!       water vapor
!
 
!       calculate coefficients at cross points for water vapor
 
    do k = 1 , kz
      do i = 2 , iym1
        if ( k > 1 ) then
          betak(i,k) = rhohf(i,k-1,j)*kvq(i,k,j)/dza(i,k-1,j)
        end if
        alphak(i,k) = egrav/(sfcps(i,j)*d_1000)/dlev(k)
      end do
    end do
 
    do k = 2 , kz - 1
      do i = 2 , iym1
        coef1(i,k) = dtpbl*alphak(i,k)*betak(i,k+1)
        coef2(i,k) = d_one+dtpbl*alphak(i,k)*(betak(i,k+1)+betak(i,k))
        coef3(i,k) = dtpbl*alphak(i,k)*betak(i,k)
      end do
    end do
 
    do i = 2 , iym1
      coef1(i,1) = dtpbl*alphak(i,1)*betak(i,2)
      coef2(i,1) = d_one + dtpbl*alphak(i,1)*betak(i,2)
      coef3(i,1) = d_zero
      coef1(i,kz) = d_zero
      coef2(i,kz) = d_one + dtpbl*alphak(i,kz)*betak(i,kz)
      coef3(i,kz) = dtpbl*alphak(i,kz)*betak(i,kz)
    end do
 
    do i = 2 , iym1
      coefe(i,1) = coef1(i,1)/coef2(i,1)
      coeff1(i,1) = qvatm(i,1,j)/coef2(i,1)
    end do
 
    do k = 2 , kz - 1
      do i = 2 , iym1
        coefe(i,k) = coef1(i,k)/(coef2(i,k)-coef3(i,k)*coefe(i,k-1))
        coeff1(i,k) = (qvatm(i,k,j)+coef3(i,k)*coeff1(i,k-1)) / &
                       (coef2(i,k)-coef3(i,k)*coefe(i,k-1))
      end do
    end do
 
    do i = 2 , iym1
      coefe(i,kz) = d_zero
      coeff1(i,kz) = (qvatm(i,kz,j) + &
               dtpbl*alphak(i,kz)*qfx(i,j) + &
               coef3(i,kz)*coeff1(i,kz-1)) /    &
               (coef2(i,kz)-coef3(i,kz)*coefe(i,kz-1))
    end do
 
!
!       all coefficients have been computed, predict field and put it in
!       temporary work space tpred
 
    do i = 2 , iym1
      tpred1(i,kz) = coeff1(i,kz)
    end do
 
    do k = kz - 1 , 1 , -1
      do i = 2 , iym1
        tpred1(i,k) = coefe(i,k)*tpred1(i,k+1) + coeff1(i,k)
      end do
    end do
 
!
!       calculate tendency due to vertical diffusion using temporary
!       predicted field
!
    if ( ibltyp == 99 ) then
      do k = 1 , kz
        do i = 2 , iym1
          diagqv(i,k,j) = (tpred1(i,k)-qvatm(i,k,j))/dtpbl*sfcps(i,j)
        end do
      end do
    else
      do k = 1 , kz
        do i = 2 , iym1
          diffq(i,k,j) = diffq(i,k,j) + &
                         (tpred1(i,k)-qvatm(i,k,j))/dtpbl*sfcps(i,j)
        end do
      end do
    end if
 
!       calculate coefficients at cross points for cloud vater
 
    do k = 1 , kz
      do i = 2 , iym1
        if ( k > 1 ) then
          betak(i,k) = rhohf(i,k-1,j)*kvq(i,k,j)/dza(i,k-1,j)
        end if
        alphak(i,k) = egrav/(sfcps(i,j)*d_1000)/dlev(k)
      end do
    end do
 
    do k = 2 , kz - 1
      do i = 2 , iym1
        coef1(i,k) = dtpbl*alphak(i,k)*betak(i,k+1)
        coef2(i,k) = d_one+dtpbl*alphak(i,k)*(betak(i,k+1)+betak(i,k))
        coef3(i,k) = dtpbl*alphak(i,k)*betak(i,k)
      end do
    end do
 
    do i = 2 , iym1
      coef1(i,1) = dtpbl*alphak(i,1)*betak(i,2)
      coef2(i,1) = d_one + dtpbl*alphak(i,1)*betak(i,2)
      coef3(i,1) = d_zero
      coef1(i,kz) = d_zero
      coef2(i,kz) = d_one + dtpbl*alphak(i,kz)*betak(i,kz)
      coef3(i,kz) = dtpbl*alphak(i,kz)*betak(i,kz)
    end do
 
    do i = 2 , iym1
      coefe(i,1) = coef1(i,1)/coef2(i,1)
      coeff1(i,1) = qcatm(i,1,j)/coef2(i,1)
    end do
 
    do k = 2 , kz - 1
      do i = 2 , iym1
        coefe(i,k) = coef1(i,k)/(coef2(i,k)-coef3(i,k)*coefe(i,k-1))
        coeff1(i,k) = (qcatm(i,k,j)+coef3(i,k)*coeff1(i,k-1)) / &
                      (coef2(i,k)-coef3(i,k)*coefe(i,k-1))
      end do
    end do
 
    do i = 2 , iym1
      coefe(i,kz) = d_zero
      coeff1(i,kz) = (qcatm(i,kz,j)+coef3(i,kz)*coeff1(i,kz-1)) / &
                     (coef2(i,kz)-coef3(i,kz)*coefe(i,kz-1))
    end do
 
!
!       all coefficients have been computed, predict field and put it in
!       temporary work space tpred
!
 
    do i = 2 , iym1
      tpred1(i,kz) = coeff1(i,kz)
    end do
 
    do k = kz - 1 , 1 , -1
      do i = 2 , iym1
        tpred1(i,k) = coefe(i,k)*tpred1(i,k+1) + coeff1(i,k)
      end do
    end do
 
!
!       calculate tendency due to vertical diffusion using temporary
!       predicted field
!
    if ( ibltyp == 99 ) then
      do k = 1 , kz
        do i = 2 , iym1
          diagqc(i,k,j) = (tpred1(i,k)-qcatm(i,k,j))/dtpbl*sfcps(i,j)
        end do
      end do
    else
      do k = 1 , kz
        do i = 2 , iym1
          qcten(i,k,j) = qcten(i,k,j) +              &
                           (tpred1(i,k)-qcatm(i,k,j))/dtpbl*sfcps(i,j)
        end do
      end do
    end if
 
!
! **********************************************************************
!
!       now add countergradient term to temperature and water vapor
!       equation
!trapuv
    do i = 2 , iym1
      ttnp(i,1) = d_zero
    end do
!trapuv_
    do k = 2 , kz
      do i = 2 , iym1
        sf = tatm(i,k,j)/thxatm(i,k,j)
        ttnp(i,k) = sf*cpd*rhohf(i,k-1,j)*kvh(i,k,j)*cgh(i,k,j)
      end do
    end do
!
!-----compute the tendencies:
!
    do i = 2 , iym1
      difft(i,kz,j) = difft(i,kz,j) - &
              egrav*ttnp(i,kz)/(d_1000*cpd*dlev(kz))
    end do
!
    do k = 1 , kzm1
      do i = 2 , iym1
        difft(i,k,j) = difft(i,k,j) + &
              egrav*(ttnp(i,k+1)-ttnp(i,k))/(d_1000*cpd*dlev(k))
      end do
    end do
 
!
!chem2
    if ( lchem .and. lchdrydepo ) then
!
!         coef1, coef2, coef3 and coefe are the same as for water vapor
!         and cloud water so they do not need to be recalculated
 
!
!         recalculation of coef1,2,3  with tracer diffusivity kvc
 
      do k = 1 , kz
        do i = 2 , iym1
          if ( k > 1 ) then
            betak(i,k) = rhohf(i,k-1,j)*kvc(i,k,j)/dza(i,k-1,j)
          end if
          alphak(i,k) = egrav/(sfcps(i,j)*d_1000)/dlev(k)
        end do
      end do
 
      do k = 2 , kz - 1
        do i = 2 , iym1
          coef1(i,k) = dtpbl*alphak(i,k)*betak(i,k+1)
          coef2(i,k) = d_one+dtpbl*alphak(i,k)*(betak(i,k+1)+betak(i,k))
          coef3(i,k) = dtpbl*alphak(i,k)*betak(i,k)
        end do
      end do
 
      do i = 2 , iym1
        coef1(i,1) = dtpbl*alphak(i,1)*betak(i,2)
        coef2(i,1) = d_one + dtpbl*alphak(i,1)*betak(i,2)
        coef3(i,1) = d_zero
        coef1(i,kz) = d_zero
        coef2(i,kz) = d_one + dtpbl*alphak(i,kz)*betak(i,kz)
        coef3(i,kz) = dtpbl*alphak(i,kz)*betak(i,kz)
      end do
!
!         set the Vd for in case of prescribed deposition velocities
!
      do itr = 1 , ntr
        do i = 2 , iym1
          if ( landmsk(i,j) == 0 ) then
            vdep(i,itr) = depvel(itr,2)
          else
            vdep(i,itr) = depvel(itr,1)
          end if
!             provisoire test de la routine chdrydep pour les dust
 
          if ( chname(itr) == 'DUST' ) vdep(i,itr) = d_zero
        end do
      end do
!
      do itr = 1 , ntr
!
        do i = 2 , iym1
          coefe(i,1) = coef1(i,1)/coef2(i,1)
          coeff1(i,1) = chmx(i,1,j,itr)/coef2(i,1)
        end do
!
        do k = 2 , kz - 1
          do i = 2 , iym1
            coefe(i,k) = coef1(i,k)/(coef2(i,k)-coef3(i,k)*coefe(i,k-1))
            coeff1(i,k) = (chmx(i,k,j,itr)+coef3(i,k)*coeff1(i,k-1)) / &
                          (coef2(i,k)-coef3(i,k)*coefe(i,k-1))
          end do
        end do
 
        do i = 2 , iym1
          coefe(i,kz) = d_zero
 
!             add dry deposition option1
          coeff1(i,kz) = (chmx(i,kz,j,itr)-dtpbl*alphak(i,kz)*chmx(i,kz,j,itr) &
                         *vdep(i,itr)*rhox2d(i,j)+coef3(i,kz)       &
                         *coeff1(i,kz-1))                           &
                         /(coef2(i,kz)-coef3(i,kz)*coefe(i,kz-1))
        end do
!
!           all coefficients have been computed, predict field and put
!           it in temporary work space tpred1
!
        do i = 2 , iym1
          tpred1(i,kz) = coeff1(i,kz)
        end do
!
        do k = kz - 1 , 1 , -1
          do i = 2 , iym1
            tpred1(i,k) = coefe(i,k)*tpred1(i,k+1) + coeff1(i,k)
          end do
        end do
!
!
!           calculate tendency due to vertical diffusion using temporary
!           predicted field
!           Dry deposition option 1 is included
 
!
        do k = 1 , kz
          do i = 2 , iym1
!qian       chten(i,k,j,itr)=chten(i,k,j,itr)
!CGAFFE     TEST diffusion/10
            chten(i,k,j,itr) = chten(i,k,j,itr) +  &
                        (tpred1(i,k)-chmx(i,k,j,itr))/dtpbl*sfcps(i,j)
          end do
        end do
        do i = 2 , iym1
 
          if ( chname(itr) /= 'DUST' ) &
            drmr(i,j,itr) = drmr(i,j,itr) + chmx(i,kz,j,itr)* &
                vdep(i,itr)*sfcps(i,j)*dtpbl*d_half*rhox2d(i,j)* &
                egrav/(sfcps(i,j)*d_1000*dlev(kz))
 
        end do
      end do
    end if
!chem2_
 
  end do
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
!                    zpbl     boundary layer height
! ------------------------------------------------------------
!
  subroutine blhnew
!
  implicit none
!
  real(8) :: fak1 , fak2 , fht , xfmt , pblk , pblk1 , pblk2 ,  &
             pfcor , phpblm , pr , therm2 , tkv , tlv , vv ,    &
             vvl , wsc , z , zh , zl , zm , zp , zzh , zzhnew , &
             zzhnew2
  integer :: i , j , k , k2
!
  pblk2 = d_zero
  zzhnew2 = d_zero
!
  do j = jbegin , jendx
!       ****note: kmxpbl, max no. of pbl levels, calculated in param
!       ******   compute richardson number
!
    therm(:) = d_zero
 
    do k = kz , kmxpbl , -1
      do i = 2 , iym1
        vv = uatm(i,k,j)*uatm(i,k,j) + &
             vatm(i,k,j)*vatm(i,k,j)
        ri(i,k) = egrav*(thvx(i,k,j)-th10(i,j))*za(i,k,j)/ &
                        (th10(i,j)*vv)
      end do
    end do
!       ******   first, set bl height to height of lowest model level
    do i = 2 , iym1
      zpbl(i,j) = za(i,kz,j)
      kpbl(i,j) = kz
    end do
!       ******   looking for bl top
    do k = kz , kmxpbl + 1 , -1
      k2 = k - 1
      do i = 2 , iym1
!     ******   bl height lies between this level and the last
!     ******   use linear interp. of rich. no. to height of ri=ricr
        if ( (ri(i,k) < ricr) .and. (ri(i,k2) >= ricr) ) then
          zpbl(i,j) = za(i,k,j) + (za(i,k2,j)-za(i,k,j)) &
              *((ricr-ri(i,k))/(ri(i,k2)-ri(i,k)))
          kpbl(i,j) = k
        end if
      end do
    end do
 
    do i = 2 , iym1
!     ******   set bl top to highest allowable model layer
      if ( ri(i,kmxpbl) < ricr ) then
        zpbl(i,j) = za(i,kmxpbl,j)
        kpbl(i,j) = kmxpbl
      end if
    end do
 
!       ******   recompute richardson no. at lowest model level
    do i = 2 , iym1
      if ( hfxv(i,j) > d_zero ) then
!           ******   estimate of convective velocity scale
        xfmt = (d_one-(binm*zpbl(i,j)/obklen(i,j)))**onet
        wsc = ustr(i,j)*xfmt
!           ******   thermal temperature excess
        therm(i) = (xhfx(i,j)+mult*thxatm(i,kz,j)*xqfx(i,j))*fak/wsc
        vvl = uatm(i,kz,j)*uatm(i,kz,j) + vatm(i,kz,j)*vatm(i,kz,j)
        ri(i,kz) = -egrav*therm(i)*za(i,kz,j)/(th10(i,j)*vvl)
      end if
    end do
 
!       ******   recompute richardson no. at other model levels
    do k = kz - 1 , kmxpbl , -1
      do i = 2 , iym1
        if ( hfxv(i,j) > d_zero ) then
          tlv = th10(i,j) + therm(i)
          tkv = thxatm(i,k,j) * &
                (d_one+mult*(qvatm(i,k,j)/(qvatm(i,k,j)+d_one)))
          vv = uatm(i,k,j)*uatm(i,k,j)+ &
               vatm(i,k,j)*vatm(i,k,j)
          ri(i,k) = egrav*(tkv-tlv)*za(i,k,j)/(th10(i,j)*vv)
        end if
      end do
    end do
 
!       ******   improve estimate of bl height under convective
!       conditions ******   using convective temperature excess (therm)
    do k = kz , kmxpbl + 1 , -1
      k2 = k - 1
      do i = 2 , iym1
        if ( hfxv(i,j) > d_zero ) then
!     ******   bl height lies between this level and the last
!     ******   use linear interp. of rich. no. to height of ri=ricr
          if ( (ri(i,k) < ricr) .and. (ri(i,k2) >= ricr) ) then
            zpbl(i,j) = za(i,k,j) + (za(i,k2,j)-za(i,k,j)) &
                           *((ricr-ri(i,k))/(ri(i,k2)-ri(i,k)))
            kpbl(i,j) = k
          end if
        end if
      end do
    end do
 
    do i = 2 , iym1
      if ( hfxv(i,j) > d_zero ) then
!     ******   set bl top to highest allowable model layer
        if ( ri(i,kmxpbl) < ricr ) then
          zpbl(i,j) = za(i,kmxpbl,j)
        end if
      end if
    end do
 
!       ******   limit bl height to be at least mech. mixing depth
    do i = 2 , iym1
!         ******   limit coriolis parameter to value at 10 deg. latitude
      pfcor = dmax1(dabs(coriolis(i,j)),2.546D-5)
!         ******   compute mechanical mixing depth,
!         ******   set to lowest model level if lower
      phpblm = 0.07D0*ustr(i,j)/pfcor
      phpblm = dmax1(phpblm,za(i,kz,j))
      zpbl(i,j) = dmax1(zpbl(i,j),phpblm)
      kpbl(i,j) = kz
    end do
 
    do k = kz , kmxpbl + 1 , -1
      k2 = k - 1
      do i = 2 , iym1
        pblk = d_zero
        zm = za(i,k,j)
        zp = za(i,k2,j)
        if ( zm < zpbl(i,j) ) then
          zp = dmin1(zp,zpbl(i,j))
          z = (zm+zp)*d_half
          zh = z/zpbl(i,j)
          zl = z/obklen(i,j)
          if ( zh <= d_one ) then
            zzh = d_one - zh
            zzh = zzh**pink
!xexp4          zzhnew = zpbl(i,j)*(d_one-zh)*zh**1.5
!xexp5          zzhnew = 0.5*zpbl(i,j)*(d_one-zh)*zh**1.5
!xexp6          zzhnew = d_one - zh
!xexp7          zzhnew =0.5* (d_one - zh)
!Sara
!               zzhnew =0.25* (d_one - zh)
!               zzhnew =0.75* (d_one - zh)
!Sara_
            zzhnew = (d_one-zh)*d_rfour
!xexp10         zzhnew =zh * (d_one - zh)**2
!chem
            if ( lchem ) then
              zzhnew2 = (d_one-zh)**d_two
            end if
!chem_
          else
            zzh = d_zero
            zzhnew = d_zero
!chem
            zzhnew2 = d_zero
!chem_
          end if
          fak1 = ustr(i,j)*zpbl(i,j)*vonkar
          if ( hfxv(i,j) <= d_zero ) then
!**             stable and neutral conditions
!**             igroup = 1
 
!**             prevent pblk from becoming too small in very stable
!               conditions
            if ( zl <= d_one ) then
              pblk = fak1*zh*zzh/(d_one+betas*zl)
!xexp5            pblk1 = vonkar * ustr(i,j) / (d_one+betas*zl) * zzhnew
              pblk1 = fak1*zh*zzhnew/(d_one+betas*zl)
!chem
              if ( lchem ) then
                pblk2 = fak1*zh*zzhnew2/(d_one+betas*zl)
              end if
!chem_
            else
              pblk = fak1*zh*zzh/(betas+zl)
!xexp5            pblk1 = vonkar * ustr(i,j) / (betas+zl) * zzhnew
              pblk1 = fak1*zh*zzhnew/(betas+zl)
!chem
              if ( lchem ) then
                pblk2 = fak1*zh*zzhnew2/(betas+zl)
              end if
!chem_
            end if
!**             compute eddy diffusivities
            kvm(i,k,j) = dmax1(pblk,kzo)
            kvh(i,k,j) = kvm(i,k,j)
            kvq(i,k,j) = dmax1(pblk1,kzo)
! Erika put k=0 in very stable conditions
            if ( zl <= 0.1D0 ) then
              kvm(i,k,j) = d_zero
              kvh(i,k,j) = kvm(i,k,j)*d_zero
              kvq(i,k,j) = d_zero
            end if
! Erika put k=0 in very stable conditions

!chem
            if ( lchem ) then
              kvc(i,k,j) = dmax1(pblk2,kzo)
            end if
!chem_
!**             compute counter-gradient term
            cgh(i,k,j) = d_zero
!               cgq(i,k) = d_zero
          else
!**             unstable conditions
 
!**             compute counter gradient term
            if ( zh >= sffrac ) then
!**               igroup = 2
              xfmt = (d_one-binm*zpbl(i,j)/obklen(i,j))**onet
              fht = dsqrt(d_one-binh*zpbl(i,j)/obklen(i,j))
              wsc = ustr(i,j)*xfmt
              pr = (xfmt/fht) + ccon
              fak2 = wsc*zpbl(i,j)*vonkar
              pblk = fak2*zh*zzh
!xexp5            pblk1 = vonkar * wsc * zzhnew
              pblk1 = fak2*zh*zzhnew
!chem
              if ( lchem ) then
                pblk2 = fak2*zh*zzhnew2
              end if
!chem_
              therm2 = fak/(zpbl(i,j)*wsc)
              cgh(i,k,j) = hfxv(i,j)*therm2
!                 cgq(i,k) = xqfx(i,j)*therm2
            else
!**               igroup = 3
              pblk = fak1*zh*zzh*(d_one-betam*zl)**onet
!xexp5            pblk1 = vonkar * ustr(i,j) * zzhnew *
!                 (d_one-betam*zl)**onet
              pblk1 = fak1*zh*zzhnew*(d_one-betam*zl)**onet
!chem
              if ( lchem ) then
                pblk2 = fak1*zh*zzhnew2*(d_one-betam*zl)**onet
              end if
!chem_
              pr = ((d_one-betam*zl)**onet)/dsqrt(d_one-betah*zl)
              cgh(i,k,j) = d_zero
!                 cgq(i,k) = d_zero
            end if
 
!**             compute eddy diffusivities
            kvm(i,k,j) = dmax1(pblk,kzo)
            kvh(i,k,j) = dmax1((pblk/pr),kzo)
!            kvq(i,k,j) = kvh(i,k,j)
            kvq(i,k,j) = dmax1(pblk1,kzo)
!chem
            if ( lchem ) then
              kvc(i,k,j) = dmax1(pblk2,kzo)
            end if
!chem_
          end if
        end if
      end do
    end do
  end do
 
  if ( ibltyp == 99 ) then
    do j = jbegin , jendx
      do i = 2 , iym1
        do k = 1 , kz
          kvm(i,k,j) = uwstateb%kzm(i,k,j)
          kvh(i,k,j) = uwstateb%kth(i,k,j)
          kvq(i,k,j) = uwstateb%kth(i,k,j)
          cgh(i,k,j) = d_zero
        end do
      end do
    end do
  end if
!
  end subroutine blhnew
!
end module mod_pbl_holtbl

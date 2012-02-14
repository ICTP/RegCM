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
 
module mod_cu_kuo
!
! This module implements Kuo cumulus parameterization scheme.
! The basic method follows Anthes and Keyser (1979) and Kuo (1983).
!
  use mod_dynparam
  use mod_memutil
  use mod_cu_common

  private

  public :: allocate_mod_cu_kuo , cupara , htdiff
!
! qdcrit : the precipitation threshold for moisture convergence.
! pert   : perturbation temperature
! perq   : perturbation mixing ratio
! dlt    : temperature difference used to allow over shooting.
! cdscld : critical cloud depth in delta sigma.
!
  real(dp) , parameter :: qdcrit = 3.0D-7
  real(dp) , parameter :: pert   = 1.0D0
  real(dp) , parameter :: perq   = 1.0D-3
  real(dp) , parameter :: dlt    = 3.0D0
  real(dp) , parameter :: cdscld = 0.3D0
!
  real(dp) , public , pointer , dimension(:,:,:) :: rsheat , rswat
  real(dp) , public , pointer , dimension(:) :: qwght
  real(dp) , public , pointer , dimension(:,:,:) :: twght , vqflx

  integer , public :: k700

  contains
!
  subroutine allocate_mod_cu_kuo
    implicit none
    call getmem3d(rsheat,1,jxp,1,iy,1,kz,'cu_kuo:rsheat')
    call getmem3d(rswat,1,jxp,1,iy,1,kz,'cu_kuo:rswat')
    call getmem1d(qwght,1,kz,'cu_kuo:qwght')
    call getmem3d(twght,1,kz,5,kz,1,kz-3,'cu_kuo:twght')
    call getmem3d(vqflx,1,kz,5,kz,1,kz-3,'cu_kuo:vqflx')
  end subroutine allocate_mod_cu_kuo
!
  subroutine cupara(jstart,jend,istart,iend,ktau)
!
!   All the other arguments are passed from subroutine "tend" and
!   explained in "tend".
!
    implicit none
!
    integer , intent(in) :: jstart , jend , istart , iend
    integer(8) , intent(in) :: ktau
!
    real(dp) :: akclth , apcnt , arh , c301 , dalr ,    &
               deqt , dlnp , dplr , dsc , e1 , eddyf , emax ,   &
               eqt , eqtm , es , plcl , pmax , prainx , psg ,   &
               psx , pux , q , qmax , qs , rh , rsht , rswt ,   &
               sca , siglcl , suma , sumb , t1 , tdmax , tlcl , &
               tmax , tmean , ttconv , ttp , ttsum , xsav , zlcl
    integer :: i , j , k , kbase , kbaseb , kclth , kk , ktop
    real(dp) , dimension(kz) :: seqt
    real(dp) , dimension(kz) :: tmp3
!
!----------------------------------------------------------------------
!
!
    pmax = d_zero
    qmax = d_zero
    tmax = d_zero
!
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
!
!   compute the moisture convergence in a column:
!   at this stage, qvten(j,i,k) only includes horizontal advection.
!   sca: is the amount of total moisture convergence
!
    total_precip_points = 0
    do i = istart , iend
      do j = jstart , jend
!
        sca = d_zero
        do k = 1 , kz
          sca = sca + qvten(j,i,k)*dflev(k)
        end do
!
!       determine if moist convection exists:
!
        if ( sca >= qdcrit ) then
!
!         check for stability
!
!         1) compute eqt (equivalent potential temperature)
!         between surface and 700 mb, with perturbation temperature
!         and moisture added. the maximum eqt will be regarded
!         as the origin of air parcel that produce cloud.
!
          eqtm = d_zero
          do k = k700 , kz
            ttp = ptatm(j,i,k)/psfcps(j,i) + pert
            q = pvqvtm(j,i,k)/psfcps(j,i) + perq
            psg = psfcps(j,i)*hlev(k) + ptop
            t1 = ttp*(d_100/psg)**rovcp
            eqt = t1*dexp(wlhvocp*q/ttp)
            if ( eqt > eqtm ) then
              eqtm = eqt
              tmax = ttp
              qmax = q
              pmax = psg
            end if
          end do
!
!         2) compute lcl, get the sigma and p of lcl
!
          emax = qmax*pmax/(ep2+qmax)
          tdmax = 5418.12D0/(19.84659D0-dlog(emax/0.611D0))
          dalr = egrav*rcpd
          dplr = (egrav*tdmax*tdmax)/(ep2*wlhv*tmax)
          zlcl = (tmax-tdmax)/(dalr-dplr)
          tlcl = tmax - dalr*zlcl
          tmean = (tmax+tlcl)*d_half
          dlnp = (egrav*zlcl)/(rgas*tmean)
          plcl = pmax*dexp(-dlnp)
          siglcl = (plcl-ptop)/psfcps(j,i)
!
!         3) compute seqt (saturation equivalent potential temperature)
!         of all the levels that are above the lcl
!
          do k = 1 , kz
            if ( hlev(k) >= siglcl ) exit
          end do
          kbase = k
          if ( kbase > kz ) kbase = kz
!
!         kbase is the layer where lcl is located.
!
          do k = 1 , kbase
            ttp = ptatm(j,i,k)/psfcps(j,i)
            psg = psfcps(j,i)*hlev(k) + ptop
            es = 0.611D0*dexp(19.84659D0-5418.12D0/ttp)
            qs = ep2*es/(psg-es)
            t1 = ttp*(d_100/psg)**rovcp
            seqt(k) = t1*dexp(wlhvocp*qs/ttp)
          end do
!
!         4) when seqt = eqt + dt, cloud top is reached.
!         eqt is the eqt of cloud (same as lcl eqt).
!
          do kk = 1 , kbase
            k = kbase + 1 - kk
            deqt = seqt(k) - eqtm
            if ( deqt > dlt ) exit
          end do
!
!         cloud top has been reached
!
          ktop = k
!
!         5) check cloud depth
!         if cloud depth is less than critical depth (cdscld = 0.3),
!         the convection is killed
!
          dsc = (siglcl-hlev(ktop))
          if ( dsc >= cdscld ) then
!
!           6) check negative area
!           if negative area is larger than the positive area
!           convection is killed.
!
            ttsum = d_zero
            do k = ktop , kbase
              ttsum = (eqtm-seqt(k))*dflev(k) + ttsum
            end do
            if ( ttsum >= d_zero ) then
!
!             you are here if stability was found.
!
!             if values dont already exist in array twght,vqflx for this
!             kbase/ktop, then flag it, and set kbase/ktop to standard
!
              if ( (kbase < 5) .or. (ktop > kbase-3) ) then
                print 99001 , ktau , i , j , kbase , ktop
                if ( kbase < 5 ) kbase = 5
                if ( ktop > kbase-3 ) ktop = kbase - 3
              end if
!
!             convection exist, compute convective flux of water vapor and
!             latent heating
!             icon   : is a counter which keep track the total points
!             where deep convection occurs.
!             c301   : is the 'b' factor in kuo's scheme.
!
              total_precip_points = total_precip_points + 1
              suma = d_zero
              sumb = d_zero
              arh = d_zero
              psx = psfcps(j,i)
              do k = 1 , kz
                qwght(k) = d_zero
              end do
              do k = ktop , kz
                pux = psx*hlev(k) + ptop
                e1 = 0.611D0*dexp(19.84659D0-5418.12D0/(ptatm(j,i,k)/psx))
                qs = ep2*e1/(pux-e1)
                rh = pvqvtm(j,i,k)/(qs*psx)
                rh = dmin1(rh,d_one)
                xsav = (d_one-rh)*qs
                qwght(k) = xsav
                sumb = sumb + qs*dflev(k)
                arh = arh + rh*qs*dflev(k)
                suma = suma + xsav*dflev(k)
              end do
              arh = arh/sumb
              c301 = d_two*(d_one-arh)
              if ( c301 < d_zero ) c301 = d_zero
              if ( c301 > d_one ) c301 = d_one
              if ( suma <= d_zero ) then
                c301 = d_zero
                suma = d_one
              end if
              do k = ktop , kz
                qwght(k) = qwght(k)/suma
              end do
              do k = 1 , kz
                ttconv = wlhvocp*(d_one-c301)*twght(k,kbase,ktop)*sca
                rsheat(j,i,k) = rsheat(j,i,k) + ttconv*dtcum*d_half
                apcnt = (d_one-c301)*sca/4.3D-3
                eddyf = apcnt*vqflx(k,kbase,ktop)
                qvten(j,i,k) = eddyf
                rswat(j,i,k) = rswat(j,i,k) + c301*qwght(k)*sca*dtcum*d_half
              end do
!
!             find cloud fractional cover and liquid water content
!
              kbaseb = min0(kbase,kzm2)
              if ( ktop <= kbaseb ) then
                kclth = kbaseb - ktop + 1
                akclth = d_one/dble(kclth)
                do k = ktop , kbaseb
                  rcldlwc(j,i,k) = cllwcv
                  rcldfra(j,i,k) = d_one - (d_one-clfrcv)**akclth
                end do
              end if
!             the unit for rainfall is mm.
              prainx = (d_one-c301)*sca*dtmdl*d_1000*regrav
              if ( prainx > dlowval ) then
                rainc(j,i) = rainc(j,i) + prainx
!               instantaneous precipitation rate for use in bats (mm/s)
                if ( ktau == 0 .and. debug_level > 2 ) then
                  lmpcpc(j,i) = lmpcpc(j,i) + prainx/dtmdl
                else
                  lmpcpc(j,i) = lmpcpc(j,i) + prainx/dtmdl*aprdiv
                end if
              end if
!
              if ( lchem ) then
                icumtop(j,i) = ktop
                icumbot(j,i) = kbaseb
              end if
              cycle
            end if
          end if
        end if
!
!       convection not exist, compute the vertical advection term:
!
        tmp3(1) = d_zero
        do k = 2 , kz
          if ( pvqvtm(j,i,k) < 1.0D-15 ) then
            tmp3(k) = d_zero
          else
            tmp3(k) = pvqvtm(j,i,k)*(pvqvtm(j,i,k-1)/pvqvtm(j,i,k))**wlev(k)
          end if
        end do
        qvten(j,i,1) = qvten(j,i,1)-svv(j,i,2)*tmp3(2)/dflev(1)
        do k = 2 , kzm1
          qvten(j,i,k) = qvten(j,i,k)-(svv(j,i,k+1)*tmp3(k+1) - &
                                       svv(j,i,k)*tmp3(k))/dflev(k)
        end do
        qvten(j,i,kz) = qvten(j,i,kz) + svv(j,i,kz)*tmp3(kz)/dflev(kz)
!
      end do
    end do
!
    do k = 1 , kz
      do i = istart , iend
        do j = jstart , jend
          rsheat(j,i,k) = dmax1(rsheat(j,i,k),d_zero)
          rswat(j,i,k) = dmax1(rswat(j,i,k),d_zero)
          rsht = rsheat(j,i,k)/tauht
          rswt = rswat(j,i,k)/tauht
          tten(j,i,k) = tten(j,i,k) + rsht
          qvten(j,i,k) = qvten(j,i,k) + rswt
          rsheat(j,i,k) = rsheat(j,i,k)*(d_one-dtcum/(d_two*tauht))
          rswat(j,i,k) = rswat(j,i,k)*(d_one-dtcum/(d_two*tauht))
        end do
      end do
    end do

99001 format (/,' >>in **cupara**: at ktau=',i8,  &
          ' (i,j)=(',i2,',',i2,'),   ',' kbase/ktop are non-standard:',2I3, &
          ' will be set to closest standard.')
!
  end subroutine cupara
!
  subroutine htdiff(wr1,wr2,dxsq,akht1,jstart,jend,istart,iend)
    implicit none
!
    integer , intent(in) :: jstart , jend , istart , iend
    real(dp) , intent(in) :: akht1 , dxsq
    real(dp) , pointer , intent(in) , dimension(:,:,:) :: wr1
    real(dp) , pointer , intent(in) , dimension(:,:,:) :: wr2
!
    integer :: i , j , k , im1 , ip1 , jm1 , jp1
!
    do k = 1 , kz
      do j = jstart , jend
#ifdef BAND
        jm1 = j - 1
        jp1 = j + 1
#else
        if ( myid == 0 ) then
          jm1 = max0(j-1,jstart)
        else
          jm1 = j - 1
        end if
        if ( myid == nproc-1 ) then
          jp1 = min0(j+1,jend)
        else
          jp1 = j + 1
        end if
#endif
        do i = istart , iend
          im1 = max0(i-1,istart)
          ip1 = min0(i+1,iend)
          rsheat(j,i,k) = rsheat(j,i,k)+akht1*dtmdl/dxsq * &
                   (wr1(j,im1,k)+wr1(j,ip1,k) + &
                    wr1(jm1,i,k)+wr1(jp1,i,k)-d_four*wr1(j,i,k))
        end do
      end do
    end do
!
    do k = 1 , kz
      do j = jstart , jend
#ifdef BAND
        jm1 = j - 1
        jp1 = j + 1
#else
        if ( myid == 0 ) then
          jm1 = max0(j-1,jstart)
        else
          jm1 = j - 1
        end if
        if ( myid == nproc-1 ) then
          jp1 = min0(j+1,jend)
        else
          jp1 = j + 1
        end if
#endif
        do i = istart , iend
          im1 = max0(i-1,istart)
          ip1 = min0(i+1,iend)
          rswat(j,i,k) = rswat(j,i,k)+akht1*dtmdl/dxsq * &
                (wr2(j,im1,k)+wr2(j,ip1,k) + &
                 wr2(jm1,i,k)+wr2(jp1,i,k)-d_four*wr2(j,i,k))
        end do
      end do
    end do
  end subroutine htdiff
!
end module mod_cu_kuo

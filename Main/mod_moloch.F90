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
!    along with ICTP RegCM.  If not, see <http://www.gnu.org/licenses/>
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_moloch

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_constants
  use mod_runparams
  use mod_mppparam
  use mod_mpmessage
  use mod_stdio
  use mod_service
  use mod_memutil
  use mod_atm_interface
  use mod_che_interface
  use mod_cu_interface
  use mod_lm_interface
  use mod_rad_interface
  use mod_pbl_interface
  use mod_micro_interface
  use mod_bdycod
  use mod_slice
  use mod_sun
  use mod_diffusion
  use mod_domain
  use mod_slabocean
  use mod_sound
  use mod_massck
  use mod_zita

  implicit none

  private

  ! generalized vertical velocity
  real(rkx) , pointer , dimension(:,:,:) :: s
  ! nonhydrostatic term in pressure gradient force
  real(rkx) , pointer , dimension(:,:,:) :: deltaw
  ! tridiagonal inversion
  real(rkx) , pointer , dimension(:,:,:) :: wwkw
  real(rkx) , pointer , dimension(:,:,:) :: wx
  real(rkx) , pointer , dimension(:,:,:) :: tkex
  real(rkx) , pointer , dimension(:,:,:) :: wz
  real(rkx) , pointer , dimension(:,:,:) :: p0
  real(rkx) , pointer , dimension(:,:,:) :: zdiv2
  real(rkx) , pointer , dimension(:,:) :: zpby
  real(rkx) , pointer , dimension(:,:) :: zpbw

  real(rkx) , pointer , dimension(:,:,:) :: pf , tf
  real(rkx) , pointer , dimension(:,:,:) :: ten0
  real(rkx) , pointer , dimension(:,:,:) :: qen0
  real(rkx) , pointer , dimension(:,:,:,:) :: chiten0

  public :: allocate_moloch , moloch
  public :: uvstagtox , wstagtox

  integer(ik4) :: nadv = 1
  integer(ik4) :: nsound = 6

  contains

#include <cpmf.inc>

  subroutine allocate_moloch
    implicit none
    call getmem3d(pf,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'moloch:pf')
    call getmem3d(tf,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'moloch:tf')
    call getmem3d(deltaw,jce1ga,jce2ga,ice1ga,ice2ga,1,kzp1,'moloch:deltaw')
    call getmem3d(s,jci1,jci2,ici1,ici2,1,kzp1,'moloch:s')
    call getmem3d(wx,jci1,jci2,ici1,ici2,1,kz,'moloch:wx')
    call getmem3d(zdiv2,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'moloch:zdiv2')
    call getmem3d(wwkw,jci1,jci2,ici1,ici2,1,kzp1,'moloch:wwkw')
    call getmem3d(wz,jci1,jci2,ice1gb,ice2gb,1,kz,'moloch:wz')
    call getmem3d(p0,jce1gb,jce2gb,ici1,ici2,1,kz,'moloch:p0')
    call getmem2d(zpby,jci1,jci2,ici1ga,ice2ga,'moloch:zpby')
    call getmem2d(zpbw,jci1ga,jce2ga,ici1,ici2,'moloch:zpbw')
    if ( ibltyp == 2 ) then
      call getmem3d(tkex,jci1,jci2,ici1,ici2,1,kz,'moloch:tkex')
    end if
    if ( idiag > 0 ) then
      call getmem3d(ten0,jci1,jci2,ici1,ici2,1,kz,'moloch:ten0')
      call getmem3d(qen0,jci1,jci2,ici1,ici2,1,kz,'moloch:qen0')
    end if
    if ( ichem == 1 ) then
      if ( ichdiag > 0 ) then
        call getmem4d(chiten0,jci1,jci2,ici1,ici2,1,kz,1,ntr,'moloch:chiten0')
      end if
    end if
  end subroutine allocate_moloch
  !
  ! Moloch dynamical integration engine
  !
  subroutine moloch
    implicit none
    integer(ik4) :: jadv , jsound
    real(rkx) :: dtsound , dtstepa
    integer(ik4) :: i , j , k
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'moloch'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    dtstepa = dtsec / real(nadv,rkx)
    dtsound = dtsec / real(nsound,rkx)

    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          mo_atm%tetav(j,i,k) =  mo_atm%tvirt(j,i,k)/mo_atm%pai(j,i,k)
          tf(j,i,k) = mo_atm%tetav(j,i,k)
        end do
      end do
    end do

    do jadv = 1 , nadv

      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            pf(j,i,k) = mo_atm%pai(j,i,k)
          end do
        end do
      end do

      call exchange(mo_atm%tetav,2,jce1,jce2,ice1,ice2,1,kz)

      call sound(dtsound)

      call advection(dtstepa)

      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            mo_atm%pai(j,i,k) = mo_atm%pai(j,i,k) - pf(j,i,k)
          end do
        end do
      end do
      call filt3d(mo_atm%pai,0.05_rkx,jci1,jci2,ici1,ici2)
      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            mo_atm%pai(j,i,k) = mo_atm%pai(j,i,k) + pf(j,i,k)
          end do
        end do
      end do

    end do ! Advection loop

    if ( mod(rcmtimer%lcount,25_ik8*int(nadv,ik8)) == 0 ) then
      call filt3d(mo_atm%u,0.06_rkx,jdi1,jdi2,ici1,ici2)
      call filt3d(mo_atm%v,0.06_rkx,jci1,jci2,idi1,idi2)
      call filt3d(mo_atm%w,0.06_rkx,jci1,jci2,ici1,ici2)
      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            mo_atm%tetav(j,i,k) = mo_atm%tetav(j,i,k) - tf(j,i,k)
          end do
        end do
      end do
      call filt3d(mo_atm%tetav,0.06_rkx,jci1,jci2,ici1,ici2)
      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            mo_atm%tetav(j,i,k) = mo_atm%tetav(j,i,k) + tf(j,i,k)
          end do
        end do
      end do
    end if

    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          mo_atm%p(j,i,k) = (mo_atm%pai(j,i,k)**cpovr) * p00
          mo_atm%tvirt(j,i,k) = mo_atm%tetav(j,i,k)*mo_atm%pai(j,i,k)
          mo_atm%t(j,i,k) = mo_atm%tvirt(j,i,k) / &
                   (d_one + ep1*mo_atm%qx(j,i,k,iqv))
        end do
      end do
    end do

    do i = ice1 , ice2
      do j = jce1 , jce2
        sfs%psa(j,i) = mo_atm%p(j,i,kz) !*exp(
      end do
    end do

    !
    ! PHYSICS
    !

    call physical_parametrizations

    !
    ! lateral boundary condition
    !

    call boundary

    !
    ! Next timestep ready : increment elapsed forecast time
    !
    call rcmtimer%advance( )
    if ( islab_ocean == 1 ) xslabtime = xslabtime + dtsec
    !
    ! calculate new solar zenith angle
    !
    call zenitm(mddom%xlat,mddom%xlon,coszrs)

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif

    contains

      subroutine boundary
        implicit none
        call exchange(mo_atm%u,1,jde1,jde2,ice1,ice2,1,kz)
        call exchange(mo_atm%v,1,jce1,jce2,ide1,ide2,1,kz)
        call exchange(mo_atm%t,1,jce1,jce2,ice1,ice2,1,kz)
        call exchange(mo_atm%qx,1,jce1,jce2,ice1,ice2,1,kz,iqv,iqv)
        if ( ichem == 1 ) then
          call exchange(mo_atm%trac,1,jce1,jce2,ice1,ice2,1,kz,1,ntr)
        end if

        if ( iboudy == 1 .or. iboudy == 5 ) then
          if ( idiag > 0 ) then
            ten0 = mo_atm%t(jci1:jci2,ici1:ici2,:)
            qen0 = mo_atm%qx(jci1:jci2,ici1:ici2,:,iqv)
          end if
          call nudge(iboudy,mo_atm%t,xtb)
          call nudge(iboudy,mo_atm%qx,xqb,iqv)
          call nudge(iboudy,mo_atm%u,mo_atm%v,xub,xvb)
          if ( idiag > 0 ) then
            tdiag%bdy = mo_atm%t(jci1:jci2,ici1:ici2,:) - ten0
            qdiag%bdy = mo_atm%qx(jci1:jci2,ici1:ici2,:,iqv) - qen0
          end if
        else if ( iboudy == 4 ) then
          call sponge(mo_atm%t,xtb)
          call sponge(mo_atm%qx,xqb,iqv)
          call sponge(mo_atm%u,mo_atm%v,xub,xvb)
          if ( idiag > 0 ) then
            tdiag%bdy = mo_atm%t(jci1:jci2,ici1:ici2,:) - ten0
            qdiag%bdy = mo_atm%qx(jci1:jci2,ici1:ici2,:,iqv) - qen0
          end if
        end if
        if ( ichem == 1 ) then
          if ( ichdiag > 0 ) then
            chiten0 = mo_atm%trac(jci1:jci2,ici1:ici2,:,:)
          end if
          if ( iboudy == 1 .or. iboudy == 5 ) then
            ! call nudge_chi(kz,mo_atm%trac)
          end if
          if ( ichdiag > 0 ) then
            cbdydiag = mo_atm%trac(jci1:jci2,ici1:ici2,:,:) - chiten0
          end if
        end if
      end subroutine boundary

      subroutine filt2d(p,anu2,j1,j2,i1,i2)
        implicit none
        real(rkx) , dimension(:,:) , pointer , intent(inout) :: p
        real(rkx) , intent(in) :: anu2
        integer(ik4) , intent(in) :: j1 , j2 , i1 , i2
        integer(ik4) :: j , i
        real(rkx) , dimension(j1:j2,i1:i2) :: p2

        call exchange(p,1,j1,j2,i1,i2)

        do i = i1 , i2
          do j = j1 , j2
            p2(j,i) = 0.125_rkx * (p(j,i-1)+p(j-1,i)+p(j+1,i)+p(j,i+1)) - &
                      0.5_rkx*p(j,i)
          end do
        end do
        do i = i1 , i2
          do j = j1 , j2
            p(j,i) = p(j,i) + anu2 * p2(j,i)
          end do
        end do
      end subroutine filt2d

      subroutine filt3d(p,anu2,j1,j2,i1,i2)
        implicit none
        real(rkx) , dimension(:,:,:) , pointer , intent(inout) :: p
        real(rkx) , intent(in) :: anu2
        integer(ik4) , intent(in) :: i1 , i2 , j1 , j2
        integer(ik4) :: j , i , k
        real(rkx) , dimension(j1:j2,i1:i2) :: p2

        call exchange(p,1,j1,j2,i1,i2,1,kz)

        do k = 1  , kz
          do i = i1 , i2
            do j = j1 , j2
              p2(j,i) = 0.125_rkx * (p(j,i-1,k) + p(j-1,i,k) + &
                                     p(j+1,i,k) + p(j,i+1,k)) - &
                        0.5_rkx   * p(j,i,k)
            end do
          end do
          do i = i1 , i2
            do j = j1 , j2
              p(j,i,k) = p(j,i,k) + anu2 * p2(j,i)
            end do
          end do
        end do
      end subroutine filt3d

      subroutine sound(dtsound)
        implicit none
        real(rkx) , intent(in) :: dtsound
        integer(ik4) :: i , j , k , km1 , kp1
        real(rkx) :: zuh , zvh , zcx , zcy , zcyp
        real(rkx) :: zrfmzu , zrfmzum , zrfmzv , zrfmzvp
        real(rkx) :: zup , zum , zvp , zvm , zdiv
        real(rkx) :: zrom1w , zwexpl , zp , zm , zrapp
        real(rkx) :: zzww0 , zzww , zfz
        real(rkx) :: zrom1u , zrom1v
        real(rkx) :: zdtrdx , zdtrdy , zdtrdz , zcs2

        zdtrdx = dtsound/dx
        zdtrdy = dtsound/dx
        zdtrdz = dtsound/mo_dz
        zcs2   = zdtrdz**2*rdrcv

        !  sound waves

        do jsound = 1 , nsound

          ! partial definition of the generalized vertical velocity
          call exchange_lr(mo_atm%u,2,jde1,jde2,ice1,ice2,1,kz)
          call exchange_bt(mo_atm%v,2,jce1,jce2,ide1,ide2,1,kz)

          do i = ici1 , ici2
            do j = jci1 , jci2
              zuh = mo_atm%u(j+1,i,kz)*mddom%hx(j+1,i) + &
                    mo_atm%u(j,i,kz)*mddom%hx(j,i)
              zvh = mo_atm%v(j,i+1,kz)*mddom%hy(j,i+1) + &
                    mo_atm%v(j,i,kz)*mddom%hy(j,i)
              mo_atm%w(j,i,kzp1) = d_half * (zuh+zvh)
              s(j,i,kzp1) = -mo_atm%w(j,i,kzp1)
            end do
          end do

          do k = kzm1, 1 , -1
            do i = ici1 , ici2
              do j = jci1 , jci2
                zuh = (mo_atm%u(j,i,k+1)+mo_atm%u(j,i,k))*mddom%hx(j,i) + &
                      (mo_atm%u(j+1,i,k+1)+mo_atm%u(j+1,i,k))*mddom%hx(j+1,i)
                zvh = (mo_atm%v(j,i,k+1)+mo_atm%v(j,i,k))*mddom%hy(j,i)+ &
                      (mo_atm%v(j,i+1,k+1)+mo_atm%v(j,i+1,k))*mddom%hy(j,i+1)
                s(j,i,k+1) = -0.25_rkx*(zuh+zvh)*gzita(zitah(k))
              end do
            end do
          end do

          ! part of divergence (except w contribution) put in zdiv2

          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
                zcx = zdtrdx*mddom%fmyu(j,i)
                zcy = zdtrdy*mddom%fmyu(j,i)*mddom%clv(j,i)
                zcyp = zdtrdy*mddom%fmyu(j,i)*mddom%clv(j,i+1)
                zrfmzu = d_two/(mo_atm%fmz(j,i,k)+mo_atm%fmz(j+1,i,k))
                zrfmzum = d_two/(mo_atm%fmz(j,i,k)+mo_atm%fmz(j-1,i,k))
                zrfmzv = d_two/(mo_atm%fmz(j,i,k)+mo_atm%fmz(j,i-1,k))
                zrfmzvp = d_two/(mo_atm%fmz(j,i,k)+mo_atm%fmz(j,i+1,k))
                zup = mo_atm%u(j+1,i,k) * zrfmzu
                zum = mo_atm%u(j,i,k) * zrfmzum
                zvp = mo_atm%v(j,i+1,k) * zrfmzvp
                zvm = mo_atm%v(j,i,k) * zrfmzv
                zdiv = zup*zcx - zum*zcx + zvp*zcyp - zvm*zcy + &
                    zdtrdz*(s(j,i,k)-s(j,i,k+1))
                zdiv2(j,i,k) = zdiv*mo_atm%fmz(j,i,k)
              end do
            end do
          end do

          call filt3d(zdiv2,0.8_rkx,jci1,jci2,ici1,ici2)

          ! new w (implicit scheme)

          do k = kz , 2 , -1
            do i = ici1 , ici2
              do j = jci1 , jci2
                deltaw(j,i,k) = -mo_atm%w(j,i,k)
                ! explicit w:
                !    it must be consistent with the initialization of pai
                zrom1w = d_half * cpd * mo_atm%fmzf(j,i,k) * &
                  (mo_atm%tetav(j,i,k-1) + mo_atm%tetav(j,i,k))
                zrom1w = zrom1w - cpd * mo_atm%w(j,i,k) * &
                  mo_atm%fmzf(j,i,k)**2 * jsound * zdtrdz * &
                  (mo_atm%tetav(j,i,k-1) - mo_atm%tetav(j,i,k)) !! GW
                zwexpl = mo_atm%w(j,i,k) - zrom1w * zdtrdz * &
                  (mo_atm%pai(j,i,k-1) - mo_atm%pai(j,i,k)) - egrav*dtsound
                zwexpl = zwexpl + rdrcv * zrom1w * zdtrdz * &
                   (mo_atm%pai(j,i,k-1) * zdiv2(j,i,k-1) - &
                    mo_atm%pai(j,i,k) * zdiv2(j,i,k))
                ! computation of the tridiagonal matrix coefficients
                ! - zp*w(k+1) + (1+zp+zm)*w(k) - zm*w(k-1) = zwexpl
                zp = zcs2 * mo_atm%fmz(j,i,k-1) * zrom1w * &
                  mo_atm%pai(j,i,k-1) + ffilt(k)
                zm = zcs2 * mo_atm%fmz(j,i,k) * zrom1w * &
                  mo_atm%pai(j,i,k) + ffilt(k)
                ! 1st loop for the tridiagonal inversion
                zrapp = d_one/(d_one + zm + zp - zm*wwkw(j,i,k+1))
                mo_atm%w(j,i,k) = zrapp * (zwexpl + zm * mo_atm%w(j,i,k+1))
                wwkw(j,i,k) = zrapp*zp
              end do
            end do
          end do

          ! 2nd loop for the tridiagonal inversion
          do k = 2, kzp1
            do i = ici1 , ici2
              do j = jci1 , jci2
                mo_atm%w(j,i,k) = mo_atm%w(j,i,k) + &
                     wwkw(j,i,k)*mo_atm%w(j,i,k-1)
                deltaw(j,i,k) = deltaw(j,i,k) + mo_atm%w(j,i,k)
              end do
            end do
          end do

          ! new Exner function

          do k = 1, kz
            do i = ici1 , ici2
              do j = jci1 , jci2
                zdiv = zdiv2(j,i,k) + zdtrdz * mo_atm%fmz(j,i,k) * &
                  (mo_atm%w(j,i,k) - mo_atm%w(j,i,k+1))
                mo_atm%pai(j,i,k) = mo_atm%pai(j,i,k) * (d_one - rdrcv*zdiv)
              end do
            end do
          end do

          ! horizontal momentum equations

          call exchange(mo_atm%pai,2,jce1,jce2,ice1,ice2,1,kz)
          call exchange(deltaw,1,jce1,jce2,ice1,ice2,1,kz)
          call exchange(mo_atm%w,1,jce1,jce2,ice1,ice2,1,kz)

          do k = 1, kz
            km1 = max(k-1,1)
            kp1 = min(k+1,kz)
            do i = ici1 , ici2
              do j = jci1 , jci2
                zcx = zdtrdx*mddom%fmyu(j,i)
                zzww0 = mo_atm%fmz(j,i,k) * jsound * zdtrdz * &
                  (mo_atm%w(j,i,k+1) * &
                    (mo_atm%tetav(j,i,km1) - mo_atm%tetav(j,i,k)) + &
                   mo_atm%w(j,i,k) * &
                    (mo_atm%tetav(j,i,k) - mo_atm%tetav(j,i,kp1)))
                zzww = zzww0 + mo_atm%fmz(j-1,i,k) * jsound * zdtrdz * &
                  (mo_atm%w(j-1,i,k+1) * &
                   (mo_atm%tetav(j-1,i,km1) - mo_atm%tetav(j-1,i,k)) + &
                   mo_atm%w(j-1,i,k) * &
                   (mo_atm%tetav(j-1,i,k) - mo_atm%tetav(j-1,i,kp1)))
                zfz = 0.25_rkx * &
                  (deltaw(j-1,i,k) + deltaw(j-1,i,k+1) + &
                   deltaw(j,i,k) + deltaw(j,i,k+1)) + egrav*dtsound
                zrom1u = d_half * cpd * &
                  (mo_atm%tetav(j-1,i,k) + mo_atm%tetav(j,i,k) - d_half*zzww)
                mo_atm%u(j,i,k) = mo_atm%u(j,i,k) - zrom1u * zcx * &
                    (mo_atm%pai(j,i,k) - mo_atm%pai(j-1,i,k)) - &
                    zfz * mddom%hx(j,i) * gzita(zita(k)) + &
                    mddom%coriol(j,i) * mo_atm%v(j,i,k) * dtsound
                zzww = zzww0 + mo_atm%fmz(j,i-1,k) * jsound * zdtrdz * &
                  (mo_atm%w(j,i-1,k+1) * &
                    (mo_atm%tetav(j,i-1,km1) - mo_atm%tetav(j,i-1,k)) + &
                   mo_atm%w(j,i-1,k) * &
                    (mo_atm%tetav(j,i-1,k) - mo_atm%tetav(j,i-1,kp1)))
                zfz = 0.25_rkx * &
                  (deltaw(j,i-1,k) + deltaw(j,i-1,k+1) + &
                   deltaw(j,i,k) + deltaw(j,i,k+1)) + egrav*dtsound
                zrom1v = d_half * cpd * &
                  (mo_atm%tetav(j,i-1,k) + mo_atm%tetav(j,i,k) - d_half*zzww)
                mo_atm%v(j,i,k) = mo_atm%v(j,i,k) - zrom1v * zdtrdy * &
                  (mo_atm%pai(j,i,k) - mo_atm%pai(j,i-1,k)) - &
                   zfz * mddom%hy(j,i) * gzita(zita(k)) - &
                   mddom%coriol(j,i) * mo_atm%u(j,i,k) * dtsound
              end do
            end do
          end do

          if ( ma%has_bdytop ) then
            do k = 1, kz
              km1 = max(k-1,1)
              kp1 = min(k+1,kz)
              do j = jci1 , jci2
                zzww0 = mo_atm%fmz(j,ice2,k) * jsound * zdtrdz * &
                  (mo_atm%w(j,ice2,k+1) * &
                    (mo_atm%tetav(j,ice2,km1) - mo_atm%tetav(j,ice2,k)) + &
                   mo_atm%w(j,ice2,k) * &
                    (mo_atm%tetav(j,ice2,k) - mo_atm%tetav(j,ice2,kp1)))
                zzww = zzww0 + mo_atm%fmz(j,ici2,k) * jsound * zdtrdz * &
                    (mo_atm%w(j,ici2,k+1) * &
                      (mo_atm%tetav(j,ici2,km1) - mo_atm%tetav(j,ici2,k)) + &
                     mo_atm%w(j,ici2,k) * &
                      (mo_atm%tetav(j,ici2,k) - mo_atm%tetav(j,ici2,kp1)))
                zfz = 0.25_rkx * &
                  (deltaw(j,ici2,k) + deltaw(j,ici2,k+1) + &
                   deltaw(j,ice2,k) + deltaw(j,ice2,k+1)) + egrav*dtsound
                zrom1v = d_half * cpd * &
                 (mo_atm%tetav(j,ici2,k) + mo_atm%tetav(j,ice2,k) - d_half*zzww)
                mo_atm%v(j,idi2,k) = mo_atm%v(j,idi2,k) - zrom1v * zdtrdy * &
                    (mo_atm%pai(j,ice2,k) - mo_atm%pai(j,ici2,k)) - &
                     zfz * mddom%hy(j,idi2) * gzita(zita(k)) - &
                     mddom%coriol(j,ice2) * mo_atm%u(j,ice2,k) * dtsound
              end do
            end do
          end if

          if ( ma%has_bdyright ) then
            do k = 1, kz
              km1 = max(k-1,1)
              kp1 = min(k+1,kz)
              do i = ici1 , ici2
                zcx = zdtrdx*mddom%fmyu(jce2,i)
                zzww0 = mo_atm%fmz(jce2,i,k) * jsound * zdtrdz * &
                  (mo_atm%w(jce2,i,k+1) * &
                    (mo_atm%tetav(jce2,i,km1) - mo_atm%tetav(jce2,i,k)) + &
                   mo_atm%w(jce2,i,k) * &
                    (mo_atm%tetav(jce2,i,k) - mo_atm%tetav(jce2,i,kp1)))
                zzww = zzww0 + mo_atm%fmz(jci2,i,k) * jsound * zdtrdz * &
                  (mo_atm%w(jci2,i,k+1) * &
                   (mo_atm%tetav(jci2,i,km1) - mo_atm%tetav(jci2,i,k)) + &
                   mo_atm%w(jci2,i,k) * &
                   (mo_atm%tetav(jci2,i,k) - mo_atm%tetav(jci2,i,kp1)))
                zfz = 0.25_rkx * &
                  (deltaw(jci2,i,k) + deltaw(jci2,i,k+1) + &
                   deltaw(jce2,i,k) + deltaw(jce2,i,k+1)) + egrav*dtsound
                zrom1u = d_half * cpd * &
                 (mo_atm%tetav(jci2,i,k) + mo_atm%tetav(jce2,i,k) - d_half*zzww)
                mo_atm%u(jdi2,i,k) = mo_atm%u(jdi2,i,k) - zrom1u * zcx * &
                    (mo_atm%pai(jce2,i,k) - mo_atm%pai(jci2,i,k)) - &
                    zfz * mddom%hx(jdi2,i) * gzita(zita(k)) + &
                    mddom%coriol(jce2,i) * mo_atm%v(jce2,i,k) * dtsound
              end do
            end do
          end if

        end do ! sound loop

        ! complete computation of generalized vertical velocity

        do k = 2 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              s(j,i,k) = (mo_atm%w(j,i,k)+s(j,i,k)) * mo_atm%fmzf(j,i,k)
            end do
          end do
        end do
        do i = ici1 , ici2
          do j = jci1 , jci2
            s(j,i,1) = d_zero
          end do
        end do
        do i = ici1 , ici2
          do j = jci1 , jci2
            s(j,i,kzp1) = d_zero
          end do
        end do
      end subroutine sound

      subroutine advection(dtstepa)
        implicit none
        integer(ik4) :: i , j , k , n
        real(rkx) :: dtstepa
        real(rkx) , pointer , dimension(:,:,:) :: ptr

        call uvstagtox(mo_atm%u,mo_atm%v,mo_atm%ux,mo_atm%vx)

        ! Compute W (and TKE if required) on zita levels

        call wstagtox(mo_atm%w,wx)

        if ( ibltyp == 2 ) then
          call wstagtox(mo_atm%tke,tkex)
        end if

        call wafone(mo_atm%tetav,mo_atm%u,mo_atm%v,mddom%clv,mddom%fmyu,dtstepa)
        call wafone(mo_atm%pai,mo_atm%u,mo_atm%v,mddom%clv,mddom%fmyu,dtstepa)
        call wafone(mo_atm%ux,mo_atm%u,mo_atm%v,mddom%clv,mddom%fmyu,dtstepa)
        call wafone(mo_atm%vx,mo_atm%u,mo_atm%v,mddom%clv,mddom%fmyu,dtstepa)
        call wafone(wx,mo_atm%u,mo_atm%v,mddom%clv,mddom%fmyu,dtstepa)
        do n = 1 , nqx
          call assignpnt(mo_atm%qx,ptr,n)
          call wafone(ptr,mo_atm%u,mo_atm%v,mddom%clv,mddom%fmyu,dtstepa)
        end do
        if ( ibltyp == 2 ) then
          call wafone(tkex,mo_atm%u,mo_atm%v,mddom%clv,mddom%fmyu,dtstepa)
        end if
        if ( ichem == 1 ) then
          do n = 1 , ntr
            call assignpnt(mo_atm%trac,ptr,n)
            call wafone(ptr,mo_atm%u,mo_atm%v,mddom%clv,mddom%fmyu,dtstepa)
          end do
        end if

        call exchange_lr(mo_atm%ux,2,jce1,jce2,ice1,ice2,1,kz)
        call exchange_bt(mo_atm%vx,2,jce1,jce2,ice1,ice2,1,kz)

        ! Back to wind points: U (fourth order)

        do k = 1 , kz
          do i = ici1 , ici2
            do j = jdii1 , jdii2
              mo_atm%u(j,i,k) = &
                0.5625_rkx * (mo_atm%ux(j,i,k)+mo_atm%ux(j-1,i,k)) - &
                0.0625_rkx * (mo_atm%ux(j+1,i,k)+mo_atm%ux(j-2,i,k))
            end do
          end do
        end do
        if ( ma%has_bdyright ) then
          do k = 1 , kz
            do i = ici1 , ici2
              mo_atm%u(jdi2,i,k) = &
                0.5_rkx * (mo_atm%ux(jci2,i,k)+mo_atm%ux(jce2,i,k))
            end do
          end do
        end if
        if ( ma%has_bdyleft ) then
          do k = 1 , kz
            do i = ici1 , ici2
              mo_atm%u(jdi1,i,k) = &
                0.5_rkx * (mo_atm%ux(jci1,i,k)+mo_atm%ux(jce1,i,k))
            end do
          end do
        end if

        ! Back to wind points: V (fourth order)
        do k = 1 , kz
          do i = idii1 , idii2
            do j = jci1 , jci2
              mo_atm%v(j,i,k) = &
                0.5625_rkx * (mo_atm%vx(j,i,k)+mo_atm%vx(j,i-1,k)) - &
                0.0625_rkx * (mo_atm%vx(j,i+1,k)+mo_atm%vx(j,i-2,k))
            end do
          end do
        end do
        if ( ma%has_bdytop ) then
          do k = 1 , kz
            do j = jci1 , jci2
              mo_atm%v(j,idi2,k) = &
                0.5_rkx * (mo_atm%vx(j,ici2,k)+mo_atm%vx(j,ice2,k))
            end do
          end do
        end if
        if ( ma%has_bdybottom ) then
          do k = 1 , kz
            do j = jci1 , jci2
              mo_atm%v(j,idi1,k) = &
                0.5_rkx * (mo_atm%vx(j,ici1,k)+mo_atm%vx(j,ice1,k))
            end do
          end do
        end if

        ! Back to half-levels

        do k = 3 , kz - 1
          do i = ici1 , ici2
            do j = jci1 , jci2
              mo_atm%w(j,i,k) = 0.5625_rkx * (wx(j,i,k)+wx(j,i,k-1)) - &
                                0.0625_rkx * (wx(j,i,k+1)+wx(j,i,k-2))
            end do
          end do
        end do
        do i = ici1 , ici2
          do j = jci1 , jci2
            mo_atm%w(j,i,2) = 0.5_rkx * (wx(j,i,2)+wx(j,i,1))
            mo_atm%w(j,i,kz) = 0.5_rkx * (wx(j,i,kz)+wx(j,i,kzm1))
          end do
        end do

        if ( ibltyp == 2 ) then
          do k = 3 , kz - 1
            do i = ici1 , ici2
              do j = jci1 , jci2
                mo_atm%tke(j,i,k) = 0.5625_rkx*(tkex(j,i,k)+tkex(j,i,k-1)) - &
                                    0.0625_rkx*(tkex(j,i,k+1)+tkex(j,i,k-2))
              end do
            end do
          end do
          do i = ici1 , ici2
            do j = jci1 , jci2
              mo_atm%tke(j,i,2) = 0.5_rkx * (tkex(j,i,2)+tkex(j,i,1))
              mo_atm%tke(j,i,kz) = 0.5_rkx * (tkex(j,i,kz)+tkex(j,i,kzm1))
            end do
          end do
        end if
      end subroutine advection

      subroutine wafone(p,u,v,clv,fmyu,dtstepa)
        implicit none
        real(rkx) , dimension(:,:,:) , pointer , intent(inout) :: p
        real(rkx) , dimension(:,:,:) , pointer , intent(in) :: u , v
        real(rkx) , dimension(:,:) , pointer , intent(in) :: clv , fmyu
        real(rkx) , intent(in) :: dtstepa
        integer(ik4) :: j , i , k
        integer(ik4) :: k1 , k1m1 , ih , ihm1 , im1 , jh , jhm1 , jm1
        real(rkx) , dimension(jci1:jci2,1:kzp1) :: wfw
        real(rkx) :: zamu , r , b , zphi , is , zdv
        real(rkx) :: zhxvt , zhxvtn
        real(rkx) :: zdtrdx , zdtrdy , zdtrdz

        zdtrdx = dtstepa/dx
        zdtrdy = dtstepa/dx
        zdtrdz = dtstepa/mo_dz

        ! Vertical advection

        do j = jci1 , jci2
          wfw(j,1) = d_zero
          wfw(j,kzp1) = d_zero
        end do

        do i = ici1 , ici2
          do k = 2 , kz
            do j = jci1 , jci2
              zamu = s(j,i,k)*zdtrdz
              if ( zamu >= d_zero ) then
                is = d_one
                k1 = k - 1
                k1m1 = k1 - 1
                k1m1 = max(k1m1,1)
              else
                is = -d_one
                k1 = k + 1
                k1m1 = k1 - 1
                k1 = min(k1,kz)
              end if
              r = rdeno(p(j,i,k1),p(j,i,k1m1),p(j,i,k),p(j,i,k-1))
              b = max(d_zero, min(d_two, max(r, min(d_two*r,d_one))))
              zphi = is + zamu * b - is * b
              wfw(j,k) = d_half * zamu*((d_one+zphi)*p(j,i,k-1) + &
                                        (d_one-zphi)*p(j,i,k))
            end do
          end do
          do k = 1 , kz
            do j = jci1 , jci2
              zdv = (s(j,i,k+1)-s(j,i,k)) * zdtrdz
              wz(j,i,k) = p(j,i,k) + wfw(j,k) - wfw(j,k+1) + p(j,i,k)*zdv
            end do
          end do
        end do

        call exchange_bt(wz,2,jci1,jci2,ici1,ici2,1,kz)

        ! Meridional advection

        do k = 1 , kz
          do i = ici1 , ice2
            do j = jci1 , jci2
              zamu = v(j,i,k)*zdtrdy
              if ( zamu > d_zero ) then
                is = d_one
                ih = max(i-1,icross1+1)
              else
                is = -d_one
                ih = min(i+1,icross2-1)
              end if
              ihm1 = max(ih-1,icross1+1)
              im1 = max(i-1,icross1+1)
              r = rdeno(wz(j,ih,k), wz(j,ihm1,k), wz(j,i,k), wz(j,im1,k))
              b = max(d_zero, min(d_two, max(r, min(d_two*r,d_one))))
              zphi = is+zamu*b - is*b
              zpby(j,i) = d_half*zamu * &
                ((d_one+zphi)*wz(j,im1,k)+(d_one-zphi)*wz(j,i,k))
            end do
          end do

          call exchange_bt(zpby,1,jci1,jci2,ici1,ice2)

          do i = ici1 , ici2
            do j = jci1 , jci2
              zhxvtn = clv(j,i+1)*fmyu(j,i)
              zhxvt  = clv(j,i)*fmyu(j,i)
              zdv = (v(j,i+1,k)*zhxvtn - v(j,i,k)*zhxvt)*zdtrdy
              p0(j,i,k) = wz(j,i,k) + &
                      zpby(j,i)*zhxvt - zpby(j,i+1)*zhxvtn + p(j,i,k)*zdv
            end do
          end do
        end do

        call exchange_lr(p0,2,jci1,jci2,ici1,ici2,1,kz)

        ! Zonal advection

        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jce2
              zamu = u(j,i,k)*zdtrdx*fmyu(j,i)
              if ( zamu > d_zero ) then
                is = d_one
                jh = max(j-1,jcross1+1)
              else
                is = -d_one
                jh = min(j+1,jcross2-1)
              end if
              jhm1 = max(jh-1,jcross1+1)
              jm1 = max(j-1,jcross1+1)
              r = rdeno(p0(jh,i,k), p0(jhm1,i,k), p0(j,i,k), p0(jm1,i,k))
              b = max(d_zero, min(d_two, max(r, min(d_two*r,d_one))))
              zphi = is+zamu*b - is*b
              zpbw(j,i) = d_half*zamu * &
                   ((d_one+zphi)*p0(jm1,i,k)+(d_one-zphi)*p0(j,i,k))
            end do
          end do

          call exchange_lr(zpbw,1,jci1,jce2,ici1,ici2)

          do i = ici1 , ici2
            do j = jci1 , jci2
              zdv = (u(j+1,i,k)-u(j,i,k))*zdtrdx*fmyu(j,i)
              p(j,i,k) = p0(j,i,k) + zpbw(j,i) - zpbw(j+1,i) + p(j,i,k)*zdv
            end do
          end do
        end do
      end subroutine wafone

  end subroutine moloch

  subroutine wstagtox(w,wx)
    implicit none
    real(rkx) , intent(in) , dimension(:,:,:) , pointer :: w
    real(rkx) , intent(inout) , dimension(:,:,:) , pointer :: wx
    integer(ik4) :: i , j , k
    do k = 2 , kz-1
      do i = ici1 , ici2
        do j = jci1 , jci2
          wx(j,i,k) = 0.5625_rkx * (w(j,i,k+1)+w(j,i,k)) - &
                      0.0625_rkx * (w(j,i,k+2)+w(j,i,k-1))
        end do
      end do
    end do
    do i = ici1 , ici2
      do j = jci1 , jci2
        wx(j,i,1)  = 0.5_rkx * (w(j,i,2)+w(j,i,1))
        wx(j,i,kz) = 0.5_rkx * (w(j,i,kzp1)+w(j,i,kz))
      end do
    end do
  end subroutine wstagtox

  subroutine uvstagtox(u,v,ux,vx)
    implicit none
    real(rkx) , intent(inout) , dimension(:,:,:) , pointer :: u , v
    real(rkx) , intent(inout) , dimension(:,:,:) , pointer :: ux , vx
    integer(ik4) :: i , j , k

    call exchange_lr(u,2,jde1,jde2,ice1,ice2,1,kz)
    call exchange_bt(v,2,jce1,jce2,ide1,ide2,1,kz)

    ! Compute U-wind on T points
    do k = 1 , kz
      do i = ice1 , ice2
        do j = jci1 , jci2
          ux(j,i,k) = 0.5625_rkx * (u(j+1,i,k)+u(j,i,k)) - &
                      0.0625_rkx * (u(j+2,i,k)+u(j-1,i,k))
        end do
      end do
    end do
    if ( ma%has_bdyleft ) then
      do k = 1 , kz
        do i = ice1 , ice2
          ux(jci1,i,k) = 0.5_rkx * (u(jde1,i,k)+u(jdi1,i,k))
        end do
      end do
    end if
    if ( ma%has_bdyright ) then
      do k = 1 , kz
        do i = ice1 , ice2
          ux(jce2,i,k) = 0.5_rkx*(u(jde2,i,k) + u(jdi2,i,k))
        end do
      end do
    end if
    ! Compute V-wind on T points
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jce1 , jce2
          vx(j,i,k) = 0.5625_rkx * (v(j,i+1,k)+v(j,i,k)) - &
                      0.0625_rkx * (v(j,i+2,k)+v(j,i-1,k))
        end do
      end do
    end do
    if ( ma%has_bdybottom ) then
      do k = 1 , kz
        do j = jce1 , jce2
          vx(j,ici1,k) = 0.5_rkx * (v(j,ide1,k)+v(j,idi1,k))
        end do
      end do
    end if
    if ( ma%has_bdytop ) then
      do k = 1 , kz
        do j = jce1 , jce2
          vx(j,ice2,k) = 0.5_rkx*(v(j,ide2,k) + v(j,idi2,k))
        end do
      end do
    end if
  end subroutine uvstagtox

  subroutine physical_parametrizations
    implicit none
  end subroutine physical_parametrizations

end module mod_moloch

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

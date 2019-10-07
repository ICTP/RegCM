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
  real(rkx) , pointer , dimension(:,:) :: wfw
  real(rkx) , pointer , dimension(:,:,:) :: p0
  real(rkx) , pointer , dimension(:,:,:) :: zdiv2
  real(rkx) , pointer , dimension(:,:) :: zpby
  real(rkx) , pointer , dimension(:,:) :: zpbw

  real(rkx) , pointer , dimension(:,:,:) :: pf , tf
  real(rkx) , pointer , dimension(:,:,:) :: ten0
  real(rkx) , pointer , dimension(:,:,:) :: qen0
  real(rkx) , pointer , dimension(:,:,:,:) :: chiten0

  real(rkx) , dimension(:,:) , pointer :: p2d
  real(rkx) , dimension(:,:) , pointer :: xlat , xlon , coriol
  real(rkx) , dimension(:,:) , pointer :: clu , fmyu , hx
  real(rkx) , dimension(:,:) , pointer :: clv , fmyv , hy
  real(rkx) , dimension(:,:,:) , pointer :: fmz
  real(rkx) , dimension(:,:,:) , pointer :: fmzf
  real(rkx) , dimension(:,:,:) , pointer :: pai
  real(rkx) , dimension(:,:,:) , pointer :: tetav
  real(rkx) , dimension(:,:,:) , pointer :: tvirt
  real(rkx) , dimension(:,:,:) , pointer :: u , v , w
  real(rkx) , dimension(:,:,:) , pointer :: ux , vx
  real(rkx) , dimension(:,:,:) , pointer :: p , t , qv , qc , qi
  real(rkx) , dimension(:,:,:) , pointer :: tke
  real(rkx) , dimension(:,:,:,:) , pointer :: qx , trac

  public :: allocate_moloch , init_moloch , moloch
  public :: uvstagtox , wstagtox

  integer(ik4) :: nadv = 1
  integer(ik4) :: nsound = 6
  real(rkx) , parameter :: minden = 1.0e-15_rkx

  contains

#include <cpmf.inc>

  subroutine allocate_moloch
    implicit none
    call getmem2d(p2d,jde1,jde2,ide1,ide2,'moloch:p2d')
    call getmem3d(pf,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'moloch:pf')
    call getmem3d(tf,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'moloch:tf')
    call getmem3d(deltaw,jce1ga,jce2ga,ice1ga,ice2ga,1,kzp1,'moloch:deltaw')
    call getmem3d(s,jci1,jci2,ici1,ici2,1,kzp1,'moloch:s')
    call getmem3d(wx,jce1,jce2,ice1,ice2,1,kz,'moloch:wx')
    call getmem3d(zdiv2,jci1ga,jci2ga,ici1ga,ici2ga,1,kz,'moloch:zdiv2')
    call getmem3d(wwkw,jci1,jci2,ici1,ici2,1,kzp1,'moloch:wwkw')
    call getmem3d(wz,jci1,jci2,ice1gb,ice2gb,1,kz,'moloch:wz')
    call getmem2d(wfw,jci1,jci2,1,kzp1,'moloch:wfw')
    call getmem3d(p0,jce1gb,jce2gb,ici1,ici2,1,kz,'moloch:p0')
    call getmem2d(zpby,jci1,jci2,ici1ga,ice2ga,'moloch:zpby')
    call getmem2d(zpbw,jci1ga,jce2ga,ici1,ici2,'moloch:zpbw')
    if ( ibltyp == 2 ) then
      call getmem3d(tkex,jce1,jce2,ice1,ice2,1,kz,'moloch:tkex')
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

  subroutine init_moloch
    implicit none
    call assignpnt(mddom%clu,clu)
    call assignpnt(mddom%clv,clv)
    call assignpnt(mddom%fmyu,fmyu)
    call assignpnt(mddom%fmyv,fmyv)
    call assignpnt(mddom%hx,hx)
    call assignpnt(mddom%hy,hy)
    call assignpnt(mddom%xlat,xlat)
    call assignpnt(mddom%xlon,xlon)
    call assignpnt(mddom%coriol,coriol)
    call assignpnt(mo_atm%fmz,fmz)
    call assignpnt(mo_atm%fmzf,fmzf)
    call assignpnt(mo_atm%pai,pai)
    call assignpnt(mo_atm%tetav,tetav)
    call assignpnt(mo_atm%u,u)
    call assignpnt(mo_atm%ux,ux)
    call assignpnt(mo_atm%v,v)
    call assignpnt(mo_atm%vx,vx)
    call assignpnt(mo_atm%w,w)
    call assignpnt(mo_atm%tvirt,tvirt)
    call assignpnt(mo_atm%p,p)
    call assignpnt(mo_atm%t,t)
    call assignpnt(mo_atm%qx,qv,iqv)
    if ( ipptls > 0 ) then
      call assignpnt(mo_atm%qx,qc,iqc)
      if ( ipptls > 1 ) call assignpnt(mo_atm%qx,qi,iqi)
    end if
    if ( ipptls > 0 ) call assignpnt(mo_atm%qx,qx)
    if ( ibltyp == 2 ) call assignpnt(mo_atm%tke,tke)
    if ( ichem == 1 ) call assignpnt(mo_atm%trac,trac)
  end subroutine init_moloch
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

    if ( ipptls > 0 ) then
      if ( ipptls > 1 ) then
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              tvirt(j,i,k) = t(j,i,k) * (d_one + ep1*qv(j,i,k) - &
                                         qc(j,i,k) - qi(j,i,k))
            end do
          end do
        end do
      else
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              tvirt(j,i,k) = t(j,i,k) * (d_one + ep1*qv(j,i,k) - qc(j,i,k))
            end do
          end do
        end do
      end if
    else
      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            tvirt(j,i,k) = t(j,i,k) * (d_one + ep1*qv(j,i,k))
          end do
        end do
      end do
    end if

    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          tetav(j,i,k) = tvirt(j,i,k)/pai(j,i,k)
          tf(j,i,k) = tetav(j,i,k)
        end do
      end do
    end do

    call exchange_lrbt(tetav,1,jce1,jce2,ice1,ice2,1,kz)

    do jadv = 1 , nadv

      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            pf(j,i,k) = pai(j,i,k)
          end do
        end do
      end do

      call sound(dtsound)

      call advection(dtstepa)

      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            pai(j,i,k) = pai(j,i,k) - pf(j,i,k)
          end do
        end do
      end do
      call filt3d(pai,0.05_rkx,jci1,jci2,ici1,ici2)
      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            pai(j,i,k) = pai(j,i,k) + pf(j,i,k)
          end do
        end do
      end do

    end do ! Advection loop

    if ( mod(rcmtimer%lcount,25_ik8*int(nadv,ik8)) == 0 ) then
      call filt3d(u,0.06_rkx,jdi1,jdi2,ici1,ici2)
      call filt3d(v,0.06_rkx,jci1,jci2,idi1,idi2)
      call filt3d(w,0.06_rkx,jci1,jci2,ici1,ici2)
      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            tetav(j,i,k) = tetav(j,i,k) - tf(j,i,k)
          end do
        end do
      end do
      call filt3d(tetav,0.06_rkx,jci1,jci2,ici1,ici2)
      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            tetav(j,i,k) = tetav(j,i,k) + tf(j,i,k)
          end do
        end do
      end do
    end if

    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          tvirt(j,i,k) = tetav(j,i,k)*pai(j,i,k)
        end do
      end do
    end do

    if ( ipptls > 0 ) then
      if ( ipptls > 1 ) then
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              t(j,i,k) = tvirt(j,i,k) / (d_one + ep1*qv(j,i,k) - &
                             qc(j,i,k) - qi(j,i,k))
            end do
          end do
        end do
      else
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              t(j,i,k) = tvirt(j,i,k) / (d_one + ep1*qv(j,i,k) - qc(j,i,k))
            end do
          end do
        end do
      end if
    else
      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            t(j,i,k) = tvirt(j,i,k) / (d_one + ep1*qv(j,i,k))
          end do
        end do
      end do
    end if

    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          p(j,i,k) = (pai(j,i,k)**cpovr) * p00
        end do
      end do
    end do

    do i = ice1 , ice2
      do j = jce1 , jce2
        sfs%psa(j,i) = p(j,i,kz) !*exp(
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
    call zenitm(xlat,xlon,coszrs)

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif

    contains

      subroutine boundary
        implicit none
        call exchange_lrbt(u,1,jde1,jde2,ice1,ice2,1,kz)
        call exchange_lrbt(v,1,jce1,jce2,ide1,ide2,1,kz)
        call exchange_lrbt(t,1,jce1,jce2,ice1,ice2,1,kz)
        call exchange_lrbt(qv,1,jce1,jce2,ice1,ice2,1,kz)
        if ( ichem == 1 ) then
          call exchange_lrbt(trac,1,jce1,jce2,ice1,ice2,1,kz,1,ntr)
        end if

        if ( iboudy == 1 .or. iboudy == 5 ) then
          if ( idiag > 0 ) then
            ten0 = t(jci1:jci2,ici1:ici2,:)
            qen0 = qv(jci1:jci2,ici1:ici2,:)
          end if
          call nudge(iboudy,t,xtb)
          call nudge(iboudy,qv,xqb)
          call nudge(iboudy,u,v,xub,xvb)
          if ( idiag > 0 ) then
            tdiag%bdy = t(jci1:jci2,ici1:ici2,:) - ten0
            qdiag%bdy = qv(jci1:jci2,ici1:ici2,:) - qen0
          end if
        else if ( iboudy == 4 ) then
          call sponge(t,xtb)
          call sponge(qv,xqb)
          call sponge(u,v,xub,xvb)
          if ( idiag > 0 ) then
            tdiag%bdy = t(jci1:jci2,ici1:ici2,:) - ten0
            qdiag%bdy = qv(jci1:jci2,ici1:ici2,:) - qen0
          end if
        end if
        if ( ichem == 1 ) then
          if ( ichdiag > 0 ) then
            chiten0 = trac(jci1:jci2,ici1:ici2,:,:)
          end if
          if ( iboudy == 1 .or. iboudy == 5 ) then
            ! call nudge_chi(kz,trac)
          end if
          if ( ichdiag > 0 ) then
            cbdydiag = trac(jci1:jci2,ici1:ici2,:,:) - chiten0
          end if
        end if
      end subroutine boundary

      subroutine filt2d(pp,anu2,j1,j2,i1,i2)
        implicit none
        real(rkx) , dimension(:,:) , pointer , intent(inout) :: pp
        real(rkx) , intent(in) :: anu2
        integer(ik4) , intent(in) :: j1 , j2 , i1 , i2
        integer(ik4) :: j , i

        call exchange_lrbt(pp,1,j1,j2,i1,i2)

        do i = i1 , i2
          do j = j1 , j2
            p2d(j,i) = 0.125_rkx * (pp(j,i-1)+pp(j-1,i)+pp(j+1,i)+pp(j,i+1)) - &
                       0.5_rkx * pp(j,i)
          end do
        end do
        do i = i1 , i2
          do j = j1 , j2
            pp(j,i) = pp(j,i) + anu2 * p2d(j,i)
          end do
        end do
      end subroutine filt2d

      subroutine filt3d(pp,anu2,j1,j2,i1,i2)
        implicit none
        real(rkx) , dimension(:,:,:) , pointer , intent(inout) :: pp
        real(rkx) , intent(in) :: anu2
        integer(ik4) , intent(in) :: i1 , i2 , j1 , j2
        integer(ik4) :: j , i , k , k1 , k2

        k1 = lbound(pp,3)
        k2 = ubound(pp,3)
        call exchange_lrbt(pp,1,j1,j2,i1,i2,k1,k2)

        do k = k1 , k2
          do i = i1 , i2
            do j = j1 , j2
              p2d(j,i) = 0.125_rkx * (pp(j-1,i,k) + pp(j+1,i,k) + &
                                      pp(j,i-1,k) + pp(j,i+1,k)) - &
                         0.5_rkx   * pp(j,i,k)
            end do
          end do
          do i = i1 , i2
            do j = j1 , j2
              pp(j,i,k) = pp(j,i,k) + anu2 * p2d(j,i)
            end do
          end do
        end do
      end subroutine filt3d

      subroutine sound(dtsound)
        implicit none
        real(rkx) , intent(in) :: dtsound
        integer(ik4) :: i , j , k , km1 , kp1
        real(rkx) :: zuh , zvh , zcx , zcxp , zcy , zcyp
        real(rkx) :: zrfmzu , zrfmzup , zrfmzv , zrfmzvp
        real(rkx) :: zup , zum , zvp , zvm , zdiv
        real(rkx) :: zrom1w , zwexpl , zp , zm , zrapp
        real(rkx) :: zzww0 , zzww , zfz , gzitak
        real(rkx) :: zrom1u , zrom1v
        real(rkx) :: zdtrdx , zdtrdy , zdtrdz , zcs2

        zdtrdx = dtsound/dx
        zdtrdy = dtsound/dx
        zdtrdz = dtsound/mo_dz
        zcs2   = zdtrdz**2*rdrcv

        !  sound waves

        do jsound = 1 , nsound

          ! partial definition of the generalized vertical velocity
          call exchange_lr(u,1,jde1,jde2,ice1,ice2,1,kz)
          call exchange_bt(v,1,jce1,jce2,ide1,ide2,1,kz)

          do i = ici1 , ici2
            do j = jci1 , jci2
              zuh = u(j,i,kz) * hx(j,i) + u(j+1,i,kz) * hx(j+1,i)
              zvh = v(j,i,kz) * hy(j,i) + v(j,i+1,kz) * hy(j,i+1)
              w(j,i,kzp1) = d_half * (zuh+zvh)
              s(j,i,kzp1) = -w(j,i,kzp1)
            end do
          end do

          ! Equation 10, generalized vertical velocity

          do k = kz , 2 , -1
            gzitak = gzita(zita(k))
            do i = ici1 , ici2
              do j = jci1 , jci2
                zuh = (u(j,i,k)   + u(j,i,k-1))   * hx(j,i) + &
                      (u(j+1,i,k) + u(j+1,i,k-1)) * hx(j+1,i)
                zvh = (v(j,i,k)   + v(j,i,k-1))   * hy(j,i) + &
                      (v(j,i+1,k) + v(j,i+1,k-1)) * hy(j,i+1)
                s(j,i,k) = -0.25_rkx * (zuh+zvh) * gzitak
              end do
            end do
          end do

          ! Part of divergence (except w contribution) put in zdiv2
          ! Equation 16

          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
                zcx  = zdtrdx * fmyu(j,i)
                zcxp = zdtrdx * fmyu(j,i)
                zcy  = zdtrdy * fmyu(j,i) * clv(j,i)
                zcyp = zdtrdy * fmyu(j,i) * clv(j,i+1)
                zrfmzu  = d_two / (fmz(j,i,k) + fmz(j-1,i,k))
                zrfmzv  = d_two / (fmz(j,i,k) + fmz(j,i-1,k))
                zrfmzup = d_two / (fmz(j,i,k) + fmz(j+1,i,k))
                zrfmzvp = d_two / (fmz(j,i,k) + fmz(j,i+1,k))
                zum = u(j,i,k) * zrfmzu
                zup = u(j+1,i,k) * zrfmzup
                zvm = v(j,i,k) * zrfmzv
                zvp = v(j,i+1,k) * zrfmzvp
                zdiv = zup*zcxp - zum*zcx + zvp*zcyp - zvm*zcy + &
                       zdtrdz * (s(j,i,k) - s(j,i,k+1))
                zdiv2(j,i,k) = zdiv * fmz(j,i,k)
              end do
            end do
          end do

          call filt3d(zdiv2,0.8_rkx,jcii1,jcii2,icii1,icii2)

          ! new w (implicit scheme) from Equation 19

          do k = kz , 2 , -1
            do i = ici1 , ici2
              do j = jci1 , jci2
                deltaw(j,i,k) = -w(j,i,k)
                ! explicit w:
                !    it must be consistent with the initialization of pai
                zrom1w = d_half * cpd * fmzf(j,i,k) * &
                        (tetav(j,i,k-1) + tetav(j,i,k))
                zrom1w = zrom1w - cpd * w(j,i,k) * fmzf(j,i,k)**2 * &
                         jsound * zdtrdz * (tetav(j,i,k-1) - tetav(j,i,k)) !! GW
                zwexpl = w(j,i,k) - zrom1w * zdtrdz * &
                         (pai(j,i,k-1) - pai(j,i,k)) - egrav*dtsound
                zwexpl = zwexpl + rdrcv * zrom1w * zdtrdz * &
                         (pai(j,i,k-1) * zdiv2(j,i,k-1) - &
                          pai(j,i,k)   * zdiv2(j,i,k))
                ! computation of the tridiagonal matrix coefficients
                ! - zp*w(k+1) + (1+zp+zm)*w(k) - zm*w(k-1) = zwexpl
                zp = zcs2 * fmz(j,i,k-1) * zrom1w * pai(j,i,k-1) + ffilt(k)
                zm = zcs2 * fmz(j,i,k)   * zrom1w * pai(j,i,k)   + ffilt(k)
                ! 1st loop for the tridiagonal inversion
                zrapp = d_one / (d_one + zm + zp - zm*wwkw(j,i,k+1))
                w(j,i,k) = zrapp * (zwexpl + zm * w(j,i,k+1))
                wwkw(j,i,k) = zrapp * zp
              end do
            end do
          end do

          ! 2nd loop for the tridiagonal inversion
          do k = 2 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
                w(j,i,k) = w(j,i,k) + wwkw(j,i,k)*w(j,i,k-1)
                deltaw(j,i,k) = deltaw(j,i,k) + w(j,i,k)
              end do
            end do
          end do

          ! new Exner function

          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
                zdiv = zdiv2(j,i,k) + zdtrdz * fmz(j,i,k) * &
                      (w(j,i,k) - w(j,i,k+1))
                pai(j,i,k) = pai(j,i,k) * (d_one - rdrcv*zdiv)
              end do
            end do
          end do

          ! horizontal momentum equations

          call exchange_lrbt(pai,1,jce1,jce2,ice1,ice2,1,kz)
          call exchange_lrbt(deltaw,1,jce1,jce2,ice1,ice2,1,kz)
          call exchange_lrbt(w,1,jce1,jce2,ice1,ice2,1,kz)

          do k = 1 , kz
            km1 = max(k-1,1)
            kp1 = min(k+1,kz)
            gzitak = gzita(zitah(k))
            do i = icii1 , icii2
              do j = jcii1 , jcii2
                zcx = zdtrdx*fmyu(j,i)
                zzww0 = fmz(j,i,k) * jsound * zdtrdz * &
                  (w(j,i,k+1) * (tetav(j,i,k)   - tetav(j,i,kp1)) + &
                   w(j,i,k)   * (tetav(j,i,km1) - tetav(j,i,k)))
                zzww = zzww0 + fmz(j-1,i,k) * jsound * zdtrdz * &
                  (w(j-1,i,k+1) * (tetav(j-1,i,k)   - tetav(j-1,i,kp1)) + &
                   w(j-1,i,k)   * (tetav(j-1,i,km1) - tetav(j-1,i,k)))
                zfz = 0.25_rkx * &
                  (deltaw(j-1,i,k) + deltaw(j-1,i,k+1) + &
                   deltaw(j,i,k)   + deltaw(j,i,k+1)) + egrav*dtsound
                zrom1u = d_half * cpd * &
                  (tetav(j-1,i,k) + tetav(j,i,k) - d_half*zzww)
                ! Equation 17
                u(j,i,k) = u(j,i,k) - &
                    zrom1u * zcx * (pai(j,i,k) - pai(j-1,i,k)) - &
                    zfz * hx(j,i) * gzitak + coriol(j,i) * v(j,i,k) * dtsound
                zzww = zzww0 + fmz(j,i-1,k) * jsound * zdtrdz * &
                  (w(j,i-1,k+1) * (tetav(j,i-1,k)   - tetav(j,i-1,kp1)) + &
                   w(j,i-1,k)   * (tetav(j,i-1,km1) - tetav(j,i-1,k)))
                zfz = 0.25_rkx * &
                  (deltaw(j,i-1,k) + deltaw(j,i-1,k+1) + &
                   deltaw(j,i,k)   + deltaw(j,i,k+1)) + egrav*dtsound
                zrom1v = d_half * cpd * &
                  (tetav(j,i-1,k) + tetav(j,i,k) - d_half*zzww)
                ! Equation 18
                v(j,i,k) = v(j,i,k) - &
                    zrom1v * zdtrdy * (pai(j,i,k) - pai(j,i-1,k)) - &
                    zfz * hy(j,i) * gzitak - coriol(j,i) * u(j,i,k) * dtsound
              end do
            end do
          end do

          if ( ma%has_bdytop ) then
            do k = 1 , kz
              km1 = max(k-1,1)
              kp1 = min(k+1,kz)
              gzitak = gzita(zitah(k))
              do j = jcii1 , jcii2
                zzww0 = fmz(j,icii2,k) * jsound * zdtrdz * &
                  (w(j,icii2,k+1) * (tetav(j,icii2,k)  - tetav(j,icii2,kp1)) + &
                   w(j,icii2,k)   * (tetav(j,icii2,km1) - tetav(j,icii2,k)))
                zzww = zzww0 + fmz(j,ici2,k) * jsound * zdtrdz * &
                  (w(j,ici2,k+1) * (tetav(j,ici2,k)   - tetav(j,ici2,kp1)) + &
                   w(j,ici2,k)   * (tetav(j,ici2,km1) - tetav(j,ici2,k)))
                zfz = 0.25_rkx * &
                  (deltaw(j,icii2,k) + deltaw(j,icii2,k+1) + &
                   deltaw(j,ici2,k) + deltaw(j,ici2,k+1)) + egrav*dtsound
                zrom1v = d_half * cpd * &
                 (tetav(j,icii2,k) + tetav(j,ici2,k) - d_half*zzww)
                v(j,idii2,k) = v(j,idii2,k) - &
                     zrom1v * zdtrdy * (pai(j,ici2,k) - pai(j,icii2,k)) - &
                     zfz * hy(j,idii2) * gzitak - &
                     coriol(j,ici2) * u(j,ici2,k) * dtsound
              end do
            end do
          end if

          if ( ma%has_bdyright ) then
            do k = 1 , kz
              km1 = max(k-1,1)
              kp1 = min(k+1,kz)
              gzitak = gzita(zitah(k))
              do i = icii1 , icii2
                zcx = zdtrdx*fmyu(jci2,i)
                zzww0 = fmz(jcii2,i,k) * jsound * zdtrdz * &
                  (w(jcii2,i,k+1) * (tetav(jcii2,i,k)  - tetav(jcii2,i,kp1)) + &
                   w(jcii2,i,k)   * (tetav(jcii2,i,km1) - tetav(jcii2,i,k)))
                zzww = zzww0 + fmz(jci2,i,k) * jsound * zdtrdz * &
                  (w(jci2,i,k+1) * (tetav(jci2,i,k)   - tetav(jci2,i,kp1)) + &
                   w(jci2,i,k)   * (tetav(jci2,i,km1) - tetav(jci2,i,k)))
                zfz = 0.25_rkx * &
                  (deltaw(jcii2,i,k) + deltaw(jcii2,i,k+1) + &
                   deltaw(jci2,i,k) + deltaw(jci2,i,k+1)) + egrav*dtsound
                zrom1u = d_half * cpd * &
                 (tetav(jcii2,i,k) + tetav(jci2,i,k) - d_half*zzww)
                u(jdii2,i,k) = u(jdii2,i,k) - &
                    zrom1u * zcx * (pai(jci2,i,k) - pai(jcii2,i,k)) - &
                    zfz * hx(jdii2,i) * gzitak + &
                    coriol(jci2,i) * v(jci2,i,k) * dtsound
              end do
            end do
          end if

        end do ! sound loop

        ! complete computation of generalized vertical velocity

        do k = 2 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              s(j,i,k) = (w(j,i,k) + s(j,i,k)) * fmzf(j,i,k)
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

        call uvstagtox(u,v,ux,vx)

        ! Compute W (and TKE if required) on zita levels

        call wstagtox(w,wx)

        if ( ibltyp == 2 ) then
          call wstagtox(tke,tkex)
        end if

        call wafone(tetav,dtstepa)
        call wafone(pai,dtstepa)
        call wafone(ux,dtstepa)
        call wafone(vx,dtstepa)
        call wafone(wx,dtstepa)
        call wafone(qv,dtstepa)
        if ( ipptls > 0 ) then
          do n = iqfrst , iqlst
            call assignpnt(qx,ptr,n)
            call wafone(ptr,dtstepa)
          end do
        end if
        if ( ibltyp == 2 ) then
          call wafone(tkex,dtstepa)
        end if
        if ( ichem == 1 ) then
          do n = 1 , ntr
            call assignpnt(trac,ptr,n)
            call wafone(ptr,dtstepa)
          end do
        end if

        call exchange_lr(ux,2,jce1,jce2,ice1,ice2,1,kz)
        call exchange_bt(vx,2,jce1,jce2,ice1,ice2,1,kz)

        ! Back to wind points: U (fourth order)

        do k = 1 , kz
          do i = ici1 , ici2
            do j = jdii1 , jdii2
              u(j,i,k) = 0.5625_rkx * (ux(j,i,k)  +ux(j-1,i,k)) - &
                         0.0625_rkx * (ux(j+1,i,k)+ux(j-2,i,k))
            end do
          end do
        end do
        if ( ma%has_bdyright ) then
          do k = 1 , kz
            do i = ici1 , ici2
              u(jdi2,i,k) = 0.5_rkx * (ux(jci2,i,k)+ux(jce2,i,k))
            end do
          end do
        end if
        if ( ma%has_bdyleft ) then
          do k = 1 , kz
            do i = ici1 , ici2
              u(jdi1,i,k) = 0.5_rkx * (ux(jci1,i,k)+ux(jce1,i,k))
            end do
          end do
        end if

        ! Back to wind points: V (fourth order)
        do k = 1 , kz
          do i = idii1 , idii2
            do j = jci1 , jci2
              v(j,i,k) = 0.5625_rkx * (vx(j,i,k)  +vx(j,i-1,k)) - &
                         0.0625_rkx * (vx(j,i+1,k)+vx(j,i-2,k))
            end do
          end do
        end do
        if ( ma%has_bdytop ) then
          do k = 1 , kz
            do j = jci1 , jci2
              v(j,idi2,k) = 0.5_rkx * (vx(j,ici2,k)+vx(j,ice2,k))
            end do
          end do
        end if
        if ( ma%has_bdybottom ) then
          do k = 1 , kz
            do j = jci1 , jci2
              v(j,idi1,k) = 0.5_rkx * (vx(j,ici1,k)+vx(j,ice1,k))
            end do
          end do
        end if

        ! Back to half-levels

        do k = 3 , kzm1
          do i = ici1 , ici2
            do j = jci1 , jci2
              w(j,i,k) = 0.5625_rkx * (wx(j,i,k)  +wx(j,i,k-1)) - &
                         0.0625_rkx * (wx(j,i,k+1)+wx(j,i,k-2))
            end do
          end do
        end do
        do i = ici1 , ici2
          do j = jci1 , jci2
            w(j,i,2) = 0.5_rkx * (wx(j,i,2)  +wx(j,i,1))
            w(j,i,kz) = 0.5_rkx * (wx(j,i,kz)+wx(j,i,kzm1))
          end do
        end do

        if ( ibltyp == 2 ) then
          do k = 3 , kzm1
            do i = ici1 , ici2
              do j = jci1 , jci2
                tke(j,i,k) = 0.5625_rkx*(tkex(j,i,k)  +tkex(j,i,k-1)) - &
                             0.0625_rkx*(tkex(j,i,k+1)+tkex(j,i,k-2))
              end do
            end do
          end do
          do i = ici1 , ici2
            do j = jci1 , jci2
              tke(j,i,2) = 0.5_rkx * (tkex(j,i,2)  +tkex(j,i,1))
              tke(j,i,kz) = 0.5_rkx * (tkex(j,i,kz)+tkex(j,i,kzm1))
            end do
          end do
        end if
      end subroutine advection

      pure real(rkx) function rdeno(t1,t2,t3,t4)
        implicit none
        real(rkx) , intent(in) :: t1 , t2 , t3 , t4
        real(rkx) :: zzden
        zzden = (t3-t4)
        rdeno = (t1-t2)/sign(max(abs(zzden),minden),zzden)
      end function rdeno

      subroutine wafone(pp,dtstepa)
        implicit none
        real(rkx) , dimension(:,:,:) , pointer , intent(inout) :: pp
        real(rkx) , intent(in) :: dtstepa
        integer(ik4) :: j , i , k
        integer(ik4) :: k1 , k1m1 , ih , ihm1 , im1 , jh , jhm1 , jm1
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
          do k = kzm1 , 1 , -1
            do j = jci1 , jci2
              zamu = s(j,i,k+1)*zdtrdz
              if ( zamu >= d_zero ) then
                is = d_one
                k1 = k + 1
                k1m1 = k1 + 1
                if ( k1m1 > kz ) k1m1 = kz
              else
                is = -d_one
                k1 = k - 1
                if ( k1 < 1 ) k1 = 1
                k1m1 = k1 + 1
              end if
              r = rdeno(pp(j,i,k1),pp(j,i,k1m1),pp(j,i,k),pp(j,i,k+1))
              b = max(d_zero, min(d_two, max(r, min(d_two*r,d_one))))
              zphi = is + zamu * b - is * b
              wfw(j,k+1) = d_half * zamu*((d_one+zphi)*pp(j,i,k+1) + &
                                          (d_one-zphi)*pp(j,i,k))
            end do
          end do
          do k = 1 , kz
            do j = jci1 , jci2
              zdv = (s(j,i,k) - s(j,i,k+1)) * zdtrdz
              wz(j,i,k) = pp(j,i,k) + wfw(j,k+1) - wfw(j,k) + pp(j,i,k)*zdv
            end do
          end do
        end do

        if ( ma%has_bdybottom ) then
          do k = 1 , kz
            do j = jci1 , jci2
              wz(j,ice1,k) = pp(j,ice1,k)
            end do
          end do
        end if
        if ( ma%has_bdytop ) then
          do k = 1 , kz
            do j = jci1 , jci2
              wz(j,ice2,k) = pp(j,ice2,k)
            end do
          end do
        end if

        call exchange_bt(wz,2,jci1,jci2,ice1,ice2,1,kz)

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
              ihm1 = ih-1
              im1 = i-1
              r = rdeno(wz(j,ih,k), wz(j,ihm1,k), wz(j,i,k), wz(j,im1,k))
              b = max(d_zero, min(d_two, max(r, min(d_two*r,d_one))))
              zphi = is + zamu*b - is*b
              zpby(j,i) = d_half*zamu * &
                ((d_one+zphi)*wz(j,im1,k) + (d_one-zphi)*wz(j,i,k))
            end do
          end do

          call exchange_bt(zpby,1,jci1,jci2,ici1,ice2)

          do i = ici1 , ici2
            do j = jci1 , jci2
              zhxvtn = clv(j,i+1)*fmyu(j,i)
              zhxvt  = clv(j,i)*fmyu(j,i)
              zdv = (v(j,i+1,k)*zhxvtn - v(j,i,k)*zhxvt)*zdtrdy
              p0(j,i,k) = wz(j,i,k) + &
                      zpby(j,i)*zhxvt - zpby(j,i+1)*zhxvtn + pp(j,i,k)*zdv
            end do
          end do
        end do

        if ( ma%has_bdyleft ) then
          do k = 1 , kz
            do i = ici1 , ici2
              p0(jce1,i,k) = pp(jce1,i,k)
            end do
          end do
        end if

        if ( ma%has_bdyright ) then
          do k = 1 , kz
            do i = ici1 , ici2
              p0(jce2,i,k) = pp(jce2,i,k)
            end do
          end do
        end if

        call exchange_lr(p0,2,jce1,jce2,ici1,ici2,1,kz)

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
              jhm1 = jh-1
              jm1 = j-1
              r = rdeno(p0(jh,i,k), p0(jhm1,i,k), p0(j,i,k), p0(jm1,i,k))
              b = max(d_zero, min(d_two, max(r, min(d_two*r,d_one))))
              zphi = is + zamu*b - is*b
              zpbw(j,i) = d_half*zamu * &
                   ((d_one+zphi)*p0(jm1,i,k) + (d_one-zphi)*p0(j,i,k))
            end do
          end do

          call exchange_lr(zpbw,1,jci1,jce2,ici1,ici2)

          do i = ici1 , ici2
            do j = jci1 , jci2
              zdv = (u(j+1,i,k) - u(j,i,k))*zdtrdx*fmyu(j,i)
              pp(j,i,k) = p0(j,i,k) + zpbw(j,i) - zpbw(j+1,i) + pp(j,i,k)*zdv
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
    do k = 2 , kzm1
      do i = lbound(wx,2) , ubound(wx,2)
        do j = lbound(wx,1) , ubound(wx,1)
          wx(j,i,k) = 0.5625_rkx * (w(j,i,k+1)+w(j,i,k)) - &
                      0.0625_rkx * (w(j,i,k+2)+w(j,i,k-1))
        end do
      end do
    end do
    do i = lbound(wx,2) , ubound(wx,2)
      do j = lbound(wx,1) , ubound(wx,1)
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
          ux(jce1,i,k) = 0.5_rkx * (u(jde1,i,k)+u(jdi1,i,k))
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
          vx(j,ice1,k) = 0.5_rkx * (v(j,ide1,k)+v(j,idi1,k))
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

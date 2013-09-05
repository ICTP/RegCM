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
 
module mod_tendency

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
  use mod_bdycod
  use mod_precip
  use mod_slice
  use mod_sun
  use mod_advection
  use mod_diffusion
  use mod_mppio
  use mod_domain
  use mod_cloud_s1
  use mod_sladvection
  use mod_slabocean
#ifdef CLM
  use mod_clm
  use mod_mtrxclm
  use clm_varsur
#endif

  private

  public :: allocate_mod_tend , tend

  real(rk8) , pointer , dimension(:,:,:) :: divl
  real(rk8) , pointer , dimension(:,:,:) :: ttld , xkc , xkcf , td , phi , ten0
  real(rk8) , pointer , dimension(:,:,:) :: ps4
  real(rk8) , pointer , dimension(:,:,:) :: ps_4 
  real(rk8) , pointer , dimension(:,:) :: psc , pten

#ifdef DEBUG
  real(rk8) , pointer , dimension(:,:,:) :: wten
#endif

  real(rk8) :: rptn ! Total number of internal points

  contains

  subroutine allocate_mod_tend
    implicit none

    call getmem3d(ps_4,jcross1,jcross2,icross1,icross2,1,4,'tendency:ps_4')
    call getmem3d(divl,jce1,jce2,ice1,ice2,1,kz,'tendency:divl')
    call getmem3d(ttld,jce1,jce2,ice1,ice2,1,kz,'tend:ttld')
    call getmem3d(xkc,jdi1,jdi2,idi1,idi2,1,kz,'tendency:xkc')
    call getmem3d(xkcf,jdi1,jdi2,idi1,idi2,1,kzp1,'tendency:xkcf')
    call getmem3d(ps4,jci1,jci2,ici1,ici2,1,4,'tendency:ps4')
    call getmem3d(td,jce1,jce2,ice1,ice2,1,kz,'tendency:td')
    call getmem2d(psc,jce1,jce2,ice1,ice2,'tendency:psc')
    call getmem2d(pten,jce1,jce2,ice1,ice2,'tendency:pten')
    call getmem3d(phi,jce1-ma%jbl1,jce2,ice1-ma%ibb1,ice2,1,kz,'tendency:phi')
    if ( idiag > 0 ) then 
      call getmem3d(ten0,jce1,jce2,ice1,ice2,1,kz,'tendency:ten0')
    end if
#ifdef DEBUG
    call getmem3d(wten,jde1,jde2,ide1,ide2,1,kz,'tendency:wten')
#endif
    rptn = d_one/dble((jout2-jout1+1)*(iout2-iout1+1))
  end subroutine allocate_mod_tend
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!                                                                     c
! This subroutine computes the tendencies of the prognostic           c
! variables p*, u, v, and t.                                          c
!                                                                     c
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine tend
    implicit none
!
    real(rk8) :: cell , chias , chibs , dudx , dudy , dvdx , dvdy ,  &
               psasum , pt2bar , pt2tot , ptnbar , maxv , lowq ,     &
               ptntot , qxas , qxbs , rovcpm , rtbar , sigpsa , tv , &
               tv1 , tv2 , tv3 , tv4 , tva , tvavg , tvb , tvc ,     &
               xmsf , theta , eccf
    integer(ik4) :: i , itr , j , k , lev , n , ii , jj , kk , iconvec
    logical :: loutrad , labsem
    character (len=32) :: appdat
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'tend'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
    !
    ! Calculate eccentricity factor for radiation calculations
    !
    theta = twopi*calday/dayspy
#ifdef CLM
    eccf  = r2ceccf
#else
    calday = yeardayfrac(idatex)
    eccf = 1.000110D0 + 0.034221D0*dcos(theta) +  &
           0.001280D0 * dsin(theta) + &
           0.000719D0 * dcos(d_two*theta) + &
           0.000077D0 * dsin(d_two*theta)
#endif
    !
    ! multiply ua and va by inverse of mapscale factor at dot point:
    !
    do k = 1 , kz
      do i = ide1 , ide2
        do j = jde1 , jde2
          atm1%u(j,i,k) = atm1%u(j,i,k)*mddom%msfd(j,i)
          atm1%v(j,i,k) = atm1%v(j,i,k)*mddom%msfd(j,i)
        end do
      end do
    end do

    call exchange(sfs%psa,1,jce1,jce2,ice1,ice2)
    call psc2psd(sfs%psa,psdot)
    !
    ! Internal U,V points
    !
    do k = 1 , kz
      do i = idii1 , idii2
        do j = jdii1 , jdii2
          xmsf = mddom%msfd(j,i)
          atmx%u(j,i,k) = atm1%u(j,i,k)/(psdot(j,i)*xmsf)
          atmx%v(j,i,k) = atm1%v(j,i,k)/(psdot(j,i)*xmsf)
        end do
      end do
    end do
    !
    ! Boundary points
    !
    if ( ma%has_bdyleft ) then
      do k = 1 , kz
        do i = idi1 , idi2
          atmx%u(jdi1,i,k) = wui(i,k)
          atmx%v(jdi1,i,k) = wvi(i,k)
        end do
      end do
      do k = 1 , kz
        do i = idi1 , idi2
          atmx%u(jde1,i,k) = wue(i,k)
          atmx%v(jde1,i,k) = wve(i,k)
        end do
      end do
      ! inflow/outflow dependence
      if ( iboudy == 3 .or. iboudy == 4 ) then
        do k = 1 , kz
          do i = idi1 , idi2
            if ( atmx%u(jde1,i,k) <= d_zero ) then
              atmx%u(jde1,i,k) = atmx%u(jdi1,i,k)
              atmx%v(jde1,i,k) = atmx%v(jdi1,i,k)
            end if
          end do
        end do
      end if
    end if
    if ( ma%has_bdyright ) then
      do k = 1 , kz
        do i = idi1 , idi2
          atmx%u(jdi2,i,k) = eui(i,k)
          atmx%v(jdi2,i,k) = evi(i,k)
        end do
      end do
      do k = 1 , kz
        do i = idi1 , idi2
          atmx%u(jde2,i,k) = eue(i,k)
          atmx%v(jde2,i,k) = eve(i,k)
        end do
      end do
      ! inflow/outflow dependence
      if ( iboudy == 3 .or. iboudy == 4 ) then
        do k = 1 , kz
          do i = idi1 , idi2
            if ( atmx%u(jde2,i,k) >= d_zero ) then
              atmx%u(jde2,i,k) = atmx%u(jdi2,i,k)
              atmx%v(jde2,i,k) = atmx%v(jdi2,i,k)
            end if
          end do
        end do
      end if
    end if
    if ( ma%has_bdybottom ) then
      do k = 1 , kz
        do j = jdi1 , jdi2
          atmx%u(j,idi1,k) = sui(j,k)
          atmx%v(j,idi1,k) = svi(j,k)
        end do
      end do
      do k = 1 , kz
        do j = jde1 , jde2
          atmx%u(j,ide1,k) = sue(j,k)
          atmx%v(j,ide1,k) = sve(j,k)
        end do
      end do
      if ( iboudy == 3 .or. iboudy == 4 ) then
        ! inflow/outflow dependence
        do k = 1 , kz
          do j = jde1 , jde2
            if ( atmx%v(j,ide1,k) >= d_zero ) then
              atmx%v(j,ide1,k) = atmx%v(j,idi1,k)
              atmx%u(j,ide1,k) = atmx%u(j,idi1,k)
            end if
          end do
        end do
      end if
    end if
    if ( ma%has_bdytop ) then
      do k = 1 , kz
        do j = jdi1 , jdi2
          atmx%u(j,idi2,k) = nui(j,k)
          atmx%v(j,idi2,k) = nvi(j,k)
        end do
      end do
      do k = 1 , kz
        do j = jde1 , jde2
          atmx%u(j,ide2,k) = nue(j,k)
          atmx%v(j,ide2,k) = nve(j,k)
        end do
      end do
      if ( iboudy == 3 .or. iboudy == 4 ) then
        ! inflow/outflow dependence
        do k = 1 , kz
          do j = jde1 , jde2
            if ( atmx%v(j,ide2,k) <= d_zero ) then
              atmx%v(j,ide2,k) = atmx%v(j,idi2,k)
              atmx%u(j,ide2,k) = atmx%u(j,idi2,k)
            end if
          end do
        end do
      end if
    end if
    !
    ! T , QV , QC decouple
    !
    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          atmx%t(j,i,k) = atm1%t(j,i,k)/sfs%psa(j,i)
        end do
      end do
    end do
    do n = 1 , nqx
      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            atmx%qx(j,i,k,n) = atm1%qx(j,i,k,n)/sfs%psa(j,i)
          end do
        end do
      end do
    end do
    !
    ! call tracer decoupling routine for multiple (ntr) species
    !
    if ( ichem == 1 ) then
      do n = 1 , ntr
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              chi(j,i,k,n) = chia(j,i,k,n)/sfs%psa(j,i)
            end do
          end do
        end do
      end do
    end if
!
!=======================================================================
!
    call exchange(sfs%psb,1,jce1,jce2,ice1,ice2)
    if ( isladvec == 1 ) then
      call exchange(atm1%u,4,jde1,jde2,ide1,ide2,1,kz)
      call exchange(atm1%v,4,jde1,jde2,ide1,ide2,1,kz)
      call exchange(atm1%t,4,jce1,jce2,ice1,ice2,1,kz)
      call exchange(atm1%qx,4,jce1,jce2,ice1,ice2,1,kz,1,nqx)
    else
      call exchange(atm1%u,1,jde1,jde2,ide1,ide2,1,kz)
      call exchange(atm1%v,1,jde1,jde2,ide1,ide2,1,kz)
      call exchange(atm1%t,1,jce1,jce2,ice1,ice2,1,kz)
      call exchange(atm1%qx,1,jce1,jce2,ice1,ice2,1,kz,1,nqx)
    end if
    if ( ibltyp == 2 .or. ibltyp == 99 ) then
      call exchange(atm1%tke,1,jce1,jce2,ice1,ice2,1,kzp1)
    end if
!
    if ( isladvec == 1 ) then
      call exchange(atm2%u,4,jde1,jde2,ide1,ide2,1,kz)
      call exchange(atm2%v,4,jde1,jde2,ide1,ide2,1,kz)
      call exchange(atm2%t,4,jce1,jce2,ice1,ice2,1,kz)
      call exchange(atm2%qx,4,jce1,jce2,ice1,ice2,1,kz,1,nqx)
    else
      call exchange(atm2%u,1,jde1,jde2,ide1,ide2,1,kz)
      call exchange(atm2%v,1,jde1,jde2,ide1,ide2,1,kz)
      call exchange(atm2%t,1,jce1,jce2,ice1,ice2,1,kz)
      call exchange(atm2%qx,1,jce1,jce2,ice1,ice2,1,kz,1,nqx)
    end if
    if ( ibltyp == 2 .or. ibltyp == 99 ) then
      call exchange(atm2%tke,1,jce1,jce2,ice1,ice2,1,kzp1)
    end if
!
    if ( isladvec == 1 ) then
      call exchange(atmx%u,4,jde1,jde2,ide1,ide2,1,kz)
      call exchange(atmx%v,4,jde1,jde2,ide1,ide2,1,kz)
      call exchange(atmx%t,4,jce1,jce2,ice1,ice2,1,kz)
      call exchange(atmx%qx,4,jce1,jce2,ice1,ice2,1,kz,1,nqx)
    else
      call exchange(atmx%u,1,jde1,jde2,ide1,ide2,1,kz)
      call exchange(atmx%v,1,jde1,jde2,ide1,ide2,1,kz)
      call exchange(atmx%t,1,jce1,jce2,ice1,ice2,1,kz)
      call exchange(atmx%qx,1,jce1,jce2,ice1,ice2,1,kz,1,nqx)
    end if
!
    if ( ichem == 1 ) then
      if ( isladvec == 1 ) then
        call exchange(chi,4,jce1,jce2,ice1,ice2,1,kz,1,ntr)
        call exchange(chib,4,jce1,jce2,ice1,ice2,1,kz,1,ntr)
      else
        call exchange(chi,1,jce1,jce2,ice1,ice2,1,kz,1,ntr)
        call exchange(chib,1,jce1,jce2,ice1,ice2,1,kz,1,ntr)
      end if
    end if
    !
    ! Calculate pdot , calculate decoupled variables at B
    !
    call exchange(sfs%psb,1,jce1,jce2,ice1,ice2)
    call psc2psd(sfs%psb,psdot)
    call mkslice
    !
    call exchange(atms%ubd3d,2,jde1,jde2,ide1,ide2,1,kz)
    call exchange(atms%vbd3d,2,jde1,jde2,ide1,ide2,1,kz)
    call exchange(atms%tb3d,2,jce1,jce2,ice1,ice2,1,kz)
    call exchange(atms%ubx3d,2,jce1,jce2,ice1,ice2,1,kz)
    call exchange(atms%vbx3d,2,jce1,jce2,ice1,ice2,1,kz)
    call exchange(atms%qxb3d,2,jce1,jce2,ice1,ice2,1,kz,1,nqx)
    if ( ibltyp == 2 .or. ibltyp == 99 ) then
      call exchange(atms%tkeb3d,2,jce1,jce2,ice1,ice2,1,kzp1)
    end if

    if ( ichem == 1 ) then
      call exchange(atms%chib3d,2,jce1,jce2,ice1,ice2,1,kz,1,ntr)
    end if
    !
    ! compute the pressure tendency
    !
    pten(:,:) = d_zero
    do k = 1 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          divl(j,i,k) = (atm1%u(j+1,i+1,k)+atm1%u(j+1,i,k)- &
                         atm1%u(j,i+1,k)  -atm1%u(j,i,k)) + &
                        (atm1%v(j+1,i+1,k)+atm1%v(j,i+1,k)- &
                         atm1%v(j+1,i,k)  -atm1%v(j,i,k))
          pten(j,i) = pten(j,i) - divl(j,i,k)*dsigma(k) / &
                      (dx2*mddom%msfx(j,i)*mddom%msfx(j,i))
        end do
      end do
    end do
    !
    ! compute vertical sigma-velocity (qdot):
    !
    qdot(:,:,:)  = d_zero
    do k = 2 , kz
      do i = ice1 , ice2
        do j = jce1 , jce2
          qdot(j,i,k) = qdot(j,i,k-1) - (pten(j,i)+divl(j,i,k-1) / &
                       (dx2*mddom%msfx(j,i)*mddom%msfx(j,i)))* &
                       dsigma(k-1)/sfs%psa(j,i)
         end do
      end do
    end do
    call exchange(qdot,1,jce1,jce2,ice1,ice2,1,kzp1)
    !
    ! compute omega
    !
    omega(:,:,:) = d_zero
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          omega(j,i,k) = d_half*sfs%psa(j,i)*(qdot(j,i,k+1)+qdot(j,i,k)) + &
                      hsigma(k)*(pten(j,i)+                         &
                     ((atmx%u(j,i,k)    +atmx%u(j,i+1,k)+           &
                       atmx%u(j+1,i+1,k)+atmx%u(j+1,i,k))*          &
                     (sfs%psa(j+1,i)-sfs%psa(j-1,i))+               &
                     (atmx%v(j,i,k)    +atmx%v(j,i+1,k)+            &
                      atmx%v(j+1,i+1,k)+atmx%v(j+1,i,k))*           &
                     (sfs%psa(j,i+1)-sfs%psa(j,i-1)))/              &
                     (dx8*mddom%msfx(j,i)))
        end do
      end do
    end do
    if ( iboudy == 4 ) then
      call sponge(ba_cr,xpsb,pten)
    else if ( iboudy == 1 .or. iboudy == 5 ) then
      call nudge(ba_cr,xbctime,sfs%psb,iboudy,xpsb,pten)
    end if
    !
    ! psc : forecast pressure
    !
    do i = ici1 , ici2
      do j = jci1 , jci2
        psc(j,i) = sfs%psb(j,i) + pten(j,i)*dt
      end do
    end do
    if ( ma%has_bdyleft ) then
      do i = ici1 , ici2
        psc(jce1,i) = sfs%psb(jce1,i) + xpsb%bt(jce1,i)*dt
      end do
    end if
    if ( ma%has_bdyright ) then
      do i = ici1 , ici2
        psc(jce2,i) = sfs%psb(jce2,i) + xpsb%bt(jce2,i)*dt
      end do
    end if
    if ( ma%has_bdybottom ) then
      do j = jce1 , jce2
        psc(j,ice1) = sfs%psb(j,ice1) + xpsb%bt(j,ice1)*dt
      end do
    end if
    if ( ma%has_bdytop ) then
      do j = jce1 , jce2
        psc(j,ice2) = sfs%psb(j,ice2) + xpsb%bt(j,ice2)*dt
      end do
    end if
    !
    ! compute bleck (1977) noise parameters:
    !
    if ( ktau /= 0 ) then
      ptntot = d_zero
      pt2tot = d_zero
      do i = ici1 , ici2
        do j = jci1 , jci2
          ptntot = ptntot + dabs(pten(j,i))
          pt2tot = pt2tot + dabs((psc(j,i)+sfs%psb(j,i)- &
                   d_two*sfs%psa(j,i))/(dt*dt*d_rfour))
        end do
      end do
    end if
    !
    ! calculate solar zenith angle
    !
    if ( ktau == 0 .or. ichem == 1 .or. &
         mod(ktau+1,ntsrf) == 0 .or. mod(ktau+1,ntrad) == 0 ) then
#ifdef CLM
      call zenit_clm(coszrs)
#else
      call zenitm(coszrs)
#endif
    end if
    !
    ! No diffusion of TKE on lower boundary (kzp1)
    !
    xkc(:,:,:) = d_zero 
    xkcf(:,:,:) = d_zero 
    !
    ! compute the horizontal diffusion coefficient and stored in xkc:
    ! the values are calculated at cross points, but they also used
    ! for dot-point variables.
    !
    do k = 2 , kz
      do i = idi1 , idi2
        do j = jdi1 , jdi2
          dudx = atm2%u(j+1,i,k) + atm2%u(j+1,i+1,k) - &
                 atm2%u(j,i,k)   - atm2%u(j,i+1,k)
          dvdx = atm2%v(j+1,i,k) + atm2%v(j+1,i+1,k) - &
                 atm2%v(j,i,k)   - atm2%v(j,i+1,k)
          dudy = atm2%u(j,i+1,k) + atm2%u(j+1,i+1,k) - &
                 atm2%u(j,i,k)   - atm2%u(j+1,i,k)
          dvdy = atm2%v(j,i+1,k) + atm2%v(j+1,i+1,k) - &
                 atm2%v(j,i,k)   - atm2%v(j+1,i,k)
          cell = (xkhz*hgfact(j,i)+c200 * &
                  dsqrt((dudx-dvdy)*(dudx-dvdy)+(dvdx+dudy)*(dvdx+dudy)))
          xkc(j,i,k) = dmin1(cell,xkhmax)
          if ( k > 1 ) then
            ! TAO: Interpolate the diffusion coefficients to full levels
            ! for use in diffusion of TKE
            cell = twt(k,1)*xkc(j,i,k) + twt(k,2)*xkc(j,i,k-1)
            xkcf(j,i,k) = dmin1(nuk*cell,xkhmax)
          else
            ! TAO: Multiply the horizontal diffusion coefficient by
            ! nuk for TKE.  Without this multiplication, it appears
            ! that TKE does not diffuse fast enough, and stabilities
            ! appear in the TKE field.  While this is different from
            ! Bretherton's treatment, it is consistent with the
            ! scaling of the vertical TKE diffusivity.
            xkcf(j,i,1) = nuk*xkc(j,i,1)
          end if
        end do
      end do
    end do
    !
    ! Initialize the tendencies
    !
    aten%u(:,:,:) = d_zero
    aten%v(:,:,:) = d_zero
    aten%t(:,:,:) = d_zero
    aten%qx(:,:,:,:) = d_zero
#ifdef DEBUG
    wten(:,:,:) = d_zero
#endif
    if ( ibltyp == 2 .or. ibltyp == 99 ) then
      aten%tke(:,:,:) = d_zero
    end if
    if ( ichem == 1 ) then 
      chiten(:,:,:,:)  = d_zero
      if ( ichdiag == 1 ) chiten0(:,:,:,:) = d_zero
    end if 
    if ( idiag > 0 ) ten0(:,:,:) = d_zero
    !
    ! Initialize diffusion terms
    !
    adf%difft(:,:,:) = d_zero
    adf%diffqx(:,:,:,:) = d_zero
    !
    ! compute the horizontal advection term:
    !  
    call hadv(cross,aten%t,atmx%t,kz)
    if ( idiag > 0 ) then
      tdiag%adh = tdiag%adh + (aten%t - ten0) * afdout
      ten0 = aten%t
    end if
#ifdef DEBUG
    call check_temperature_tendency('HADV')
#endif
    !
    ! compute the vertical advection term:
    !
    if ( ibltyp /= 2 .and. ibltyp /= 99 ) then
      call vadv(cross,aten%t,atm1%t,kz,1)
    else
      if ( iuwvadv == 1 ) then
        call vadv(cross,aten%t,atm1%t,kz,3)
      else
        call vadv(cross,aten%t,atm1%t,kz,1)
      end if
    end if
    if ( idiag > 0 ) then
      tdiag%adv = tdiag%adv + (aten%t - ten0) * afdout
      ten0 = aten%t
    end if
#ifdef DEBUG
    call check_temperature_tendency('VADV')
#endif
    !
    ! compute the adiabatic term:
    !
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          rovcpm = rgas/(cpd*(d_one+0.8D0*(atmx%qx(j,i,k,iqv))))
          tv = atmx%t(j,i,k)*(d_one+ep1*(atmx%qx(j,i,k,iqv)))
          aten%t(j,i,k) = aten%t(j,i,k) + (omega(j,i,k)*rovcpm*tv) / &
                          (ptop/sfs%psa(j,i)+hsigma(k))
        end do
      end do
    end do
    if ( idiag > 0 ) then
      tdiag%adi = tdiag%adi + (aten%t - ten0) * afdout
      ten0 = aten%t
    end if
#ifdef DEBUG
    call check_temperature_tendency('DIAB')
#endif

    if ( isladvec == 1 ) then
      call trajcalc_x
    end if
    !
    ! compute the diffusion term for t and store in difft:
    !
    if ( idiag > 0 ) ten0(jci1:jci2,ici1:ici2,:) = adf%difft
    call diffu_x(adf%difft,atms%tb3d,sfs%psb,xkc,kz)
    if ( idiag > 0 ) then
      ! save the h diff diag here
      tdiag%dif(jci1:jci2,ici1:ici2,:) = tdiag%dif(jci1:jci2,ici1:ici2,:) + &
        (adf%difft - ten0(jci1:jci2,ici1:ici2,:)) * afdout
      ! reset ten0 to aten%t 
      ten0 = aten%t
    end if
    !
    ! compute the moisture tendencies for convection
    !
    if ( isladvec == 1 ) then
      call slhadv_x(aten%qx,atm2%qx,iqv)
      call hdvg_x(aten%qx,atm1%qx,iqv)
    else
      call hadv(aten%qx,atmx%qx,kz,iqv)
    end if
    if ( icup /= 1 ) then
      if ( ibltyp /= 2 .and. ibltyp /= 99 ) then
        call vadv(aten%qx,atm1%qx,kz,1,iqv)
      else
        if ( iuwvadv == 1 ) then
          call vadv(aten%qx,atm1%qx,kz,3,iqv)
        else
          call vadv(aten%qx,atm1%qx,kz,1,iqv)
        end if
      end if
    end if
    !
    ! Zero out radiative clouds
    !
    cldfra(:,:,:) = d_zero
    cldlwc(:,:,:) = d_zero
    !
    ! conv tracer diagnostic
    !
    if ( ichem == 1 .and. ichdiag == 1 ) chiten0 = chiten
    !
    ! Call cumulus parametrization
    !
    if ( icup == 1 ) then
      call cupara
    end if
    if ( icup == 2 .or. icup == 99 .or. icup == 98 .or. icup == 96 ) then
      call cuparan
    end if
    if ( icup == 3 ) then
      call bmpara
    end if
    if ( icup == 4 .or. icup == 99 .or. icup == 98 .or. icup == 97 ) then
      call cupemandrv
    end if
    if ( icup == 5 .or. icup == 96 .or. icup == 97 ) then
      ! Tiedtke needs surface heat fluxes to be computed.
      if ( ktau > 0 ) call tiedtkedrv
    end if
 
    if ( idiag > 0 ) then
      tdiag%con = tdiag%con + (aten%t - ten0) * afdout
      ten0 = aten%t
    end if
#ifdef DEBUG
    call check_temperature_tendency('CONV')
    call check_wind_tendency('CONV')
#endif
    !
    ! save cumulus cloud fraction for chemistry before it is
    ! overwritten in cldfrac 
    !
    if ( ichem == 1 .and. ichdiag == 1 ) then
      cconvdiag = cconvdiag + (chiten - chiten0) * cfdout
    end if
    if ( ichem == 1 ) convcldfra(:,:,:) = cldfra(:,:,:)
    !
    ! Large scale precipitation
    !
    if ( ipptls > 0 ) then
      if ( isladvec == 1 ) then
        call slhadv_x(aten%qx,atm2%qx,iqfrst,iqlst)
        call hdvg_x(aten%qx,atm1%qx,iqfrst,iqlst)
      else
        call hadv(aten%qx,atmx%qx,kz,iqfrst,iqlst)
      end if
      if ( ibltyp /= 2 .and. ibltyp /= 99 ) then
        call vadv(aten%qx,atm1%qx,kz,2,iqfrst,iqlst)
      else
        if ( iuwvadv == 1 ) then
          call vadv(aten%qx,atm1%qx,kz,3,iqfrst,iqlst)
        else
          call vadv(aten%qx,atm1%qx,kz,2,iqfrst,iqlst)
        end if
      end if
      if ( ipptls == 2 ) then
        call cldfrac
        call microphys(omega)
      else  
        call pcp
        call cldfrac
      end if
      ! 
      ! compute the diffusion terms:
      ! the diffusion term for qx is stored in diffqx.
      !
      call diffu_x(adf%diffqx,atms%qxb3d,sfs%psb,xkc,kz)
    end if

    if ( idiag > 0 ) then
      ! save tten from pcp (evaporation)
      tdiag%lsc = tdiag%lsc + (aten%t - ten0) * afdout
      ten0 = aten%t
    end if
#ifdef DEBUG
    call check_temperature_tendency('PREC')
#endif
    !
    ! Tracers tendencies
    !
    if ( ichem == 1 ) then
      !
      ! horizontal and vertical advection + diag
      !
      if ( ichdiag == 1 ) chiten0 = chiten
      if ( isladvec == 1 ) then
        call slhadv_x(chiten,chib)
        call hdvg_x(chiten,chia)
      else
        call hadv(chiten,chi,kz)
      end if
      if ( ichdiag == 1 ) then
        cadvhdiag = cadvhdiag + (chiten - chiten0) * cfdout
        chiten0 = chiten
      end if
      if ( icup /= 1 ) then
        if ( ibltyp /= 2 .and. ibltyp /= 99 ) then
          call vadv(chiten,chia,kz,2)
        else
          if ( iuwvadv == 1 ) then
            call vadv(chiten,chia,kz,3)
          else
            call vadv(chiten,chia,kz,2)
          end if
        end if
      end if
      if ( ichdiag == 1 ) then
        cadvvdiag = cadvvdiag + (chiten - chiten0) * cfdout
        chiten0 = chiten
      end if
      !
      ! horizontal diffusion: initialize scratch vars to 0.
      ! need to compute tracer tendencies due to diffusion
      !
      call diffu_x(chiten,chib3d,sfs%psb,xkc,kz)
      if ( ichdiag == 1 ) then
        cdifhdiag = cdifhdiag + (chiten - chiten0) * cfdout 
      end if
      !
      ! Compute chemistry tendencies (other than transport)
      !
      call tractend2(ktau,xyear,xmonth,xday,calday,declin)
    end if ! ichem
    !
    ! call radiative transfer package
    !
    if ( ktau == 0 .or. mod(ktau+1,ntrad) == 0 ) then
      !
      ! calculate albedo
      !
#ifndef CLM
      call albedobats
#else
      call albedoclm
#endif
      loutrad = (ktau == 0 .or. mod(ktau+1,krad) == 0)
      if ( iclimao3 == 1 ) then
        call read_o3data(idatex,scenario,mddom%xlat,mddom%xlon, &
          sfs%psa,ptop,sigma)
      end if
      if ( irrtm == 1 ) then
        call rrtmg_driver(xyear,eccf,loutrad)
      else
        labsem = ( ktau == 0 .or. mod(ktau+1,ntabem) == 0 )
        call colmod3(xyear,eccf,loutrad,labsem)
      end if
    end if

#ifndef CLM
    !
    ! call mtrxbats for surface physics calculations
    !
    if ( ktau == 0 .or. mod(ktau+1,ntsrf) == 0 ) then
      call mtrxbats(ktau)
      if ( islab_ocean == 1 ) call update_slabocean(xslabtime)
    end if
#else
    !
    ! call mtrxclm for surface physics calculations
    !
    if ( ktau == 0 .or. mod(ktau+1,ntrad) == 0 ) then
      r2cdoalb = .true.
    else
      r2cdoalb = .false.
    end if
    if ( ktau == 0 .or. mod(ktau+1,ntsrf) == 0 ) then
      ! Timestep used is the same as for bats
      if ( ktau == 0 ) then
        r2cnstep = 0
      else
        r2cnstep = (ktau+1)/ntsrf
      end if
      call mtrxclm(ktau)
      if ( islab_ocean == 1 ) call update_slabocean(xslabtime)
    end if
#endif
    if ( icup == 1 ) then
      call htdiff(dxsq,akht1)
    end if
    !
    ! Call medium resolution PBL
    !
    if ( ichem == 1 .and. ichdiag == 1 ) chiten0 = chiten
    ! care : pbl update the difft table at this level 
    if ( idiag > 0 ) ten0(jci1:jci2,ici1:ici2,:) = adf%difft
    if ( ibltyp == 2 .or. ibltyp == 99 ) then
      ! Call the Grenier and Bretherton (2001) / Bretherton (2004) TCM
      call uwtcm
      call uvcross2dot(uwten%u,uwten%v,aten%u,aten%v)
      call get_data_from_tcm(uwstateb,uwten,aten,atm1,atm2,.true.)
    end if
    if ( ibltyp == 1 .or. ibltyp == 99 ) then
      ! Call the Holtslag PBL
      call exchange(sfs%psb,1,jce1,jce2,ice1,ice2)
      call psc2psd(sfs%psb,psdot)
      call holtbl
    end if
    if ( ibltyp == 99 ) then
      call check_conserve_qt(holtten%qx,uwten,uwstateb,kz)
      adf%diffqx(:,:,:,:) = adf%diffqx(:,:,:,:) + holtten%qx(:,:,:,:)
    end if
    if ( ichem == 1 .and. ichdiag == 1 ) then
      ctbldiag = ctbldiag + (chiten - chiten0) * cfdout 
    end if
    if ( idiag > 0 ) then
      tdiag%tbl(jci1:jci2,ici1:ici2,:) = tdiag%tbl(jci1:jci2,ici1:ici2,:) + &
                     (adf%difft - ten0(jci1:jci2,ici1:ici2,:)) * afdout    
      ten0 = aten%t 
    end if
#ifdef DEBUG
    call check_temperature_tendency('PBLL')
    call check_wind_tendency('PBLL')
#endif
    !
    ! add ccm radiative transfer package-calculated heating rates to
    ! temperature tendency
    !
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          ! heating rate in deg/sec
          aten%t(j,i,k) = aten%t(j,i,k) + sfs%psb(j,i)*heatrt(j,i,k)
        end do
      end do
    end do
    if ( idiag > 0 ) then
      tdiag%rad = tdiag%rad + (aten%t - ten0) * afdout
      ten0 = aten%t
    end if
#ifdef DEBUG
    call check_temperature_tendency('HEAT')
#endif
    !
    ! apply the sponge boundary conditions on t and qv.
    ! We must do this before adding diffusion terms.
    !
    if ( iboudy == 4 ) then
      call sponge(kz,ba_cr,xtb,aten%t)
      call sponge(kz,iqv,ba_cr,xqb,aten%qx)
      if ( idiag > 0 ) then    
        ! rq : temp condensation tend is added the evap temp tend
        !      calculated in pcp
        tdiag%bdy = tdiag%bdy + (aten%t - ten0) * afdout
        ten0 = aten%t
      end if
    end if

    !
    ! add horizontal diffusion and pbl tendencies for t and qv

    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          aten%t(j,i,k) = aten%t(j,i,k) + adf%difft(j,i,k)
        end do
      end do
    end do
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          aten%qx(j,i,k,iqv) = aten%qx(j,i,k,iqv) + adf%diffqx(j,i,k,iqv)
        end do
      end do
    end do
    ! Rq: the temp tendency diagnostics have been already
    !     saved for tbl and hor. diff but :  
    if ( idiag > 0 ) ten0 = aten%t ! important since aten%t have been updated 
    !
    ! compute the condensation and precipitation terms for explicit
    ! moisture scheme:
    !
    if ( ipptls > 0 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            aten%qx(j,i,k,iqc) = aten%qx(j,i,k,iqc) + adf%diffqx(j,i,k,iqc)
          end do
        end do
      end do
      if ( ipptls == 2 ) then
        !dummystatement
      else
        call condtq(psc)
      end if
      if ( idiag > 0 ) then    
        ! rq : temp condensation tend is added the evap temp tend
        ! calculated in pcp
        tdiag%lsc = tdiag%lsc + (aten%t - ten0) * afdout
        ten0 = aten%t
      end if
    end if
    !
    ! apply the nudging boundary conditions:
    !
    if ( iboudy == 1 .or. iboudy == 5 ) then
      call nudge(kz,ba_cr,xbctime,atm2%t,iboudy,xtb,aten%t)
      call nudge(kz,iqv,ba_cr,xbctime,atm2%qx,iboudy,xqb,aten%qx)
    end if
    if ( ichem == 1 ) then
      if ( ichdiag == 1 ) chiten0 = chiten
      ! keep nudge_chi for now 
      if ( iboudy == 1 .or. iboudy == 5 ) then
        call nudge_chi(kz,xbctime,chib,chiten)
      end if
      if ( ichdiag == 1 ) cbdydiag = cbdydiag + (chiten - chiten0) * cfdout
    end if  
    if ( idiag > 0 ) then    
      tdiag%bdy = tdiag%bdy + (aten%t - ten0) * afdout
      ten0 = aten%t
    end if
#ifdef DEBUG
    call check_temperature_tendency('BDYC')
#endif
    !
    ! forecast t, qv, and qc at tau+1:
    !
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          atmc%t(j,i,k) = atm2%t(j,i,k) + dt*aten%t(j,i,k)
        end do
      end do
    end do
    do n = 1 , nqx
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            atmc%qx(j,i,k,n) = atm2%qx(j,i,k,n) + dt*aten%qx(j,i,k,n)
          end do
        end do
      end do
    end do
    !
    ! forecast tracer chi at at tau+1:
    !
    if ( ichem == 1 ) then
      do itr = 1 , ntr
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              chic(j,i,k,itr) = chib(j,i,k,itr) + dt*chiten(j,i,k,itr)
            end do
          end do
        end do
      end do
    end if
    !
    ! compute weighted p*t (td) for use in ssi:
    !
    if ( ipgf == 1 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            tvc = atmc%t(j,i,k)*(d_one+ep1*(atmc%qx(j,i,k,iqv))/psc(j,i))
            tva = atm1%t(j,i,k)*(d_one+ep1*(atmx%qx(j,i,k,iqv)))
            tvb = atm2%t(j,i,k)*(d_one+ep1* &
                                 (atm2%qx(j,i,k,iqv))/sfs%psb(j,i))
            td(j,i,k) = alpha*(tvc+tvb) + beta*tva
            ttld(j,i,k) = td(j,i,k) - sfs%psa(j,i) * &
                      t00pg*((hsigma(k)*sfs%psa(j,i)+ptop)/p00pg)**pgfaa1
          end do
        end do
      end do
      if ( ma%has_bdyleft ) then
        do k = 1 , kz
          do i = ici1 , ici2
            td(jce1,i,k) = atm1%t(jce1,i,k)*(d_one+ep1*(atmx%qx(jce1,i,k,iqv)))
            ttld(jce1,i,k) = td(jce1,i,k) - sfs%psa(jce1,i) * &
                        t00pg*((hsigma(k)*sfs%psa(jce1,i)+ptop)/p00pg)**pgfaa1
          end do
        end do
      end if
      if ( ma%has_bdyright ) then
        do k = 1 , kz
          do i = ici1 , ici2
            td(jce2,i,k) = atm1%t(jce2,i,k)*(d_one+ep1*(atmx%qx(jce2,i,k,iqv)))
            ttld(jce2,i,k) = td(jce2,i,k) - sfs%psa(jce2,i) * &
                     t00pg*((hsigma(k)*sfs%psa(jce2,i)+ptop)/p00pg)**pgfaa1
          end do
        end do
      end if
      if ( ma%has_bdybottom ) then
        do k = 1 , kz
          do j = jce1 , jce2
            td(j,ice1,k) = atm1%t(j,ice1,k)*(d_one+ep1*(atmx%qx(j,ice1,k,iqv)))
            ttld(j,ice1,k) = td(j,ice1,k) - sfs%psa(j,ice1) * &
                     t00pg*((hsigma(k)*sfs%psa(j,ice1)+ptop)/p00pg)**pgfaa1
          end do
        end do
      end if
      if ( ma%has_bdytop ) then
        do k = 1 , kz
          do j = jce1 , jce2
            td(j,ice2,k) = atm1%t(j,ice2,k)*(d_one+ep1*(atmx%qx(j,ice2,k,iqv)))
            ttld(j,ice2,k) = td(j,ice2,k) - sfs%psa(j,ice2) * &
                     t00pg*((hsigma(k)*sfs%psa(j,ice2)+ptop)/p00pg)**pgfaa1
          end do
        end do
      end if
    else if ( ipgf == 0 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            tvc = atmc%t(j,i,k)*(d_one+ep1*(atmc%qx(j,i,k,iqv))/psc(j,i))
            tva = atm1%t(j,i,k)*(d_one+ep1*(atmx%qx(j,i,k,iqv)))
            tvb = atm2%t(j,i,k)*(d_one+ep1* &
                 (atm2%qx(j,i,k,iqv))/sfs%psb(j,i))
            td(j,i,k) = alpha*(tvc+tvb) + beta*tva
          end do
        end do
      end do
      if ( ma%has_bdyleft ) then
        do k = 1 , kz
          do i = ici1 , ici2
            td(jce1,i,k) = atm1%t(jce1,i,k)*(d_one+ep1*(atmx%qx(jce1,i,k,iqv)))
          end do
        end do
      end if
      if ( ma%has_bdyright ) then
        do k = 1 , kz
          do i = ici1 , ici2
            td(jce2,i,k) = atm1%t(jce2,i,k)*(d_one+ep1*(atmx%qx(jce2,i,k,iqv)))
          end do
        end do
      end if
      if ( ma%has_bdybottom ) then
        do k = 1 , kz
          do j = jce1 , jce2
            td(j,ice1,k) = atm1%t(j,ice1,k)*(d_one+ep1*(atmx%qx(j,ice1,k,iqv)))
          end do
        end do
      end if
      if ( ma%has_bdytop ) then
        do k = 1 , kz
          do j = jce1 , jce2
            td(j,ice2,k) = atm1%t(j,ice2,k)*(d_one+ep1*(atmx%qx(j,ice2,k,iqv)))
          end do
        end do
      end if
    end if
    !
    ! compute the u and v tendencies:
    !
    !   compute the diffusion terms:
    !   put diffusion and pbl tendencies of u and v in difuu and difuv.
    !
    do k = 1 , kz
      do i = idi1 , idi2
        do j = jdi1 , jdi2
          adf%difuu(j,i,k) = aten%u(j,i,k)
          adf%difuv(j,i,k) = aten%v(j,i,k)
        end do
      end do
    end do
    !
    call diffu_d(adf%difuu,atms%ubd3d,psdot,mddom%msfd,xkc,1)
    call diffu_d(adf%difuv,atms%vbd3d,psdot,mddom%msfd,xkc,1)
    !
    ! Reset tendencies for U,V
    !
    aten%u(:,:,:) = d_zero
    aten%v(:,:,:) = d_zero
    !
    ! compute the horizontal advection terms for u and v:
    !
    call hadv(dot,aten%u,atmx%u,kz)
    call hadv(dot,aten%v,atmx%v,kz)
#ifdef DEBUG
    call check_wind_tendency('HADV')
#endif
    !
    ! compute coriolis terms:
    !
    do k = 1 , kz
      do i = idi1 , idi2
        do j = jdi1 , jdi2
          aten%u(j,i,k) = aten%u(j,i,k) + &
                       mddom%coriol(j,i)*atm1%v(j,i,k)/mddom%msfd(j,i)
          aten%v(j,i,k) = aten%v(j,i,k) - &
                       mddom%coriol(j,i)*atm1%u(j,i,k)/mddom%msfd(j,i)
        end do
      end do
    end do
#ifdef DEBUG
    call check_wind_tendency('CORI')
#endif
    !
    ! compute pressure gradient terms:
    !
    if ( ipgf == 1 ) then
      do k = 1 , kz
        do i = idi1 , idi2
          do j = jdi1 , jdi2
            psasum = sfs%psa(j,i) + sfs%psa(j,i-1) + &
                     sfs%psa(j-1,i) + sfs%psa(j-1,i-1)
            sigpsa = psasum
            tv1 = atmx%t(j-1,i-1,k)*(d_one+ep1*(atmx%qx(j-1,i-1,k,iqv)))
            tv2 = atmx%t(j-1,i,k)*(d_one+ep1*(atmx%qx(j-1,i,k,iqv)))
            tv3 = atmx%t(j,i-1,k)*(d_one+ep1*(atmx%qx(j,i-1,k,iqv)))
            tv4 = atmx%t(j,i,k)*(d_one+ep1*(atmx%qx(j,i,k,iqv)))
            rtbar = tv1 + tv2 + tv3 + tv4 - d_four*t00pg*             &
                    ((hsigma(k)*psasum*d_rfour+ptop)/p00pg)**pgfaa1
            rtbar = rgas*rtbar*sigpsa/16.0D0
            aten%u(j,i,k) = aten%u(j,i,k) - rtbar *               &
                  (dlog(d_half*(sfs%psa(j,i)+sfs%psa(j,i-1))*     &
                        hsigma(k)+ptop) -                         &
                   dlog(d_half*(sfs%psa(j-1,i)+sfs%psa(j-1,i-1))* &
                        hsigma(k)+ptop))/(dx*mddom%msfd(j,i))
            aten%v(j,i,k) = aten%v(j,i,k) - rtbar *               &
                  (dlog(d_half*(sfs%psa(j,i)+sfs%psa(j-1,i))*     &
                        hsigma(k)+ptop) -                         &
                   dlog(d_half*(sfs%psa(j-1,i-1)+sfs%psa(j,i-1))* &
                        hsigma(k)+ptop))/(dx*mddom%msfd(j,i))
          end do
        end do
      end do
    else if ( ipgf == 0 ) then
      do k = 1 , kz
        do i = idi1 , idi2
          do j = jdi1 , jdi2
            psasum = sfs%psa(j,i) + sfs%psa(j,i-1) + &
                     sfs%psa(j-1,i) + sfs%psa(j-1,i-1)
            sigpsa = psasum
            tv1 = atmx%t(j-1,i-1,k)*(d_one+ep1*(atmx%qx(j-1,i-1,k,iqv)))
            tv2 = atmx%t(j-1,i,k)*(d_one+ep1*(atmx%qx(j-1,i,k,iqv)))
            tv3 = atmx%t(j,i-1,k)*(d_one+ep1*(atmx%qx(j,i-1,k,iqv)))
            tv4 = atmx%t(j,i,k)*(d_one+ep1*(atmx%qx(j,i,k,iqv)))
            rtbar = rgas*(tv1+tv2+tv3+tv4)*sigpsa/16.0D0
            aten%u(j,i,k) = aten%u(j,i,k) - rtbar *                &
                   (dlog(d_half*(sfs%psa(j,i)+sfs%psa(j,i-1))*     &
                         hsigma(k)+ptop) -                         &
                    dlog(d_half*(sfs%psa(j-1,i)+sfs%psa(j-1,i-1))* &
                         hsigma(k)+ptop))/(dx*mddom%msfd(j,i))
            aten%v(j,i,k) = aten%v(j,i,k) - rtbar *                &
                   (dlog(d_half*(sfs%psa(j,i)+sfs%psa(j-1,i))*     &
                         hsigma(k)+ptop) -                         &
                    dlog(d_half*(sfs%psa(j-1,i-1)+sfs%psa(j,i-1))* &
                         hsigma(k)+ptop))/(dx*mddom%msfd(j,i))
          end do
        end do
      end do
    end if
#ifdef DEBUG
    call check_wind_tendency('PRGR')
#endif
    !
    ! compute geopotential height at half-k levels, cross points:
    !
    if ( ipgf == 1 ) then
      do i = ice1 , ice2
        do j = jce1 , jce2
          tv = (ttld(j,i,kz)/sfs%psa(j,i))/(d_one+atmx%qx(j,i,kz,iqc)/ &
                                       (d_one+atmx%qx(j,i,kz,iqv)))
          phi(j,i,kz) = mddom%ht(j,i) + &
                   rgas*t00pg/pgfaa1*((sfs%psa(j,i)+ptop)/p00pg)**pgfaa1
          phi(j,i,kz) = phi(j,i,kz) - rgas * tv * &
                  dlog((hsigma(kz)+ptop/sfs%psa(j,i))/(d_one+ptop/sfs%psa(j,i)))
        end do
      end do
      do k = 1 , kzm1
        lev = kz - k
        do i = ice1 , ice2
          do j = jce1 , jce2
            tvavg = ((ttld(j,i,lev)*dsigma(lev)+ttld(j,i,lev+1)*   &
                    dsigma(lev+1))/(sfs%psa(j,i)*(dsigma(lev)+     &
                    dsigma(lev+1))))/(d_one+atmx%qx(j,i,lev,iqc)/  &
                                     (d_one+atmx%qx(j,i,lev,iqv)))
            phi(j,i,lev) = phi(j,i,lev+1) - rgas *             &
                   tvavg*dlog((hsigma(lev)+ptop/sfs%psa(j,i))/ &
                             (hsigma(lev+1)+ptop/sfs%psa(j,i)))
          end do
        end do
      end do
    else if ( ipgf == 0 ) then
      do i = ice1 , ice2
        do j = jce1 , jce2
          tv = (td(j,i,kz)/sfs%psa(j,i))/(d_one+atmx%qx(j,i,kz,iqc)/  &
                                     (d_one+atmx%qx(j,i,kz,iqv)))
          phi(j,i,kz) = mddom%ht(j,i) - rgas * tv * &
               dlog((hsigma(kz)+ptop/sfs%psa(j,i))/(d_one+ptop/sfs%psa(j,i)))
        end do
      end do
      do k = 1 , kzm1
        lev = kz - k
        do i = ice1 , ice2
          do j = jce1 , jce2
            tvavg = ((td(j,i,lev)*dsigma(lev)+td(j,i,lev+1)*       &
                    dsigma(lev+1))/(sfs%psa(j,i)*(dsigma(lev)+     &
                    dsigma(lev+1))))/(d_one+atmx%qx(j,i,lev,iqc) / &
                                     (d_one+atmx%qx(j,i,lev,iqv)))
            phi(j,i,lev) = phi(j,i,lev+1) - rgas *              &
                   tvavg*dlog((hsigma(lev)+ptop/sfs%psa(j,i)) / &
                              (hsigma(lev+1)+ptop/sfs%psa(j,i)))
          end do
        end do
      end do
    end if
    call exchange_lb(phi,1,jce1,jce2,ice1,ice2,1,kz)
    !
    ! compute the geopotential gradient terms:
    !
    do k = 1 , kz
      do i = idi1 , idi2
        do j = jdi1 , jdi2
          aten%u(j,i,k) = aten%u(j,i,k) -         &
               (sfs%psa(j-1,i-1)+sfs%psa(j-1,i)+  &
                sfs%psa(j,i-1)+sfs%psa(j,i)) *    &
               (phi(j,i,k)+phi(j,i-1,k)-          &
                phi(j-1,i,k)-phi(j-1,i-1,k)) / (dx8*mddom%msfd(j,i))
          aten%v(j,i,k) = aten%v(j,i,k) -         &
               (sfs%psa(j-1,i-1)+sfs%psa(j-1,i)+  &
                sfs%psa(j,i-1)+sfs%psa(j,i)) *    &
               (phi(j,i,k)+phi(j-1,i,k)-          &
                phi(j,i-1,k)-phi(j-1,i-1,k)) / (dx8*mddom%msfd(j,i))
        end do
      end do
    end do
#ifdef DEBUG
    call check_wind_tendency('GEOP')
#endif
    !
    ! compute the vertical advection terms:
    !
    call vadv(dot,aten%u,atm1%u,kz,2)
    call vadv(dot,aten%v,atm1%v,kz,2)
#ifdef DEBUG
    call check_wind_tendency('VADV')
#endif
    !
    ! apply the sponge boundary condition on u and v:
    !
    if ( iboudy == 4 ) then
      call sponge(kz,ba_dt,xub,aten%u)
      call sponge(kz,ba_dt,xvb,aten%v)
    end if
    !
    ! apply the nudging boundary conditions:
    !
    if ( iboudy == 1 .or. iboudy == 5 ) then
      call nudge(kz,ba_dt,xbctime,atm2%u,iboudy,xub,aten%u)
      call nudge(kz,ba_dt,xbctime,atm2%v,iboudy,xvb,aten%v)
    end if
#ifdef DEBUG
    call check_wind_tendency('BDYC')
#endif
    !
    ! add the diffusion and pbl tendencies to aten%u and aten%v:
    !
    do k = 1 , kz
      do i = idi1 , idi2
        do j = jdi1 , jdi2
          aten%u(j,i,k) = aten%u(j,i,k) + adf%difuu(j,i,k)
          aten%v(j,i,k) = aten%v(j,i,k) + adf%difuv(j,i,k)
        end do
      end do
    end do
#ifdef DEBUG
    call check_wind_tendency('DIFF')
#endif
    !
    ! forecast p*u and p*v at tau+1:
    !
    do k = 1 , kz
      do i = idi1 , idi2
        do j = jdi1 , jdi2
          atmc%u(j,i,k) = atm2%u(j,i,k) + dt*aten%u(j,i,k)
          atmc%v(j,i,k) = atm2%v(j,i,k) + dt*aten%v(j,i,k)
        end do
      end do
    end do
    !
    !  Couple TKE to ps for use in vertical advection
    !
    if ( ibltyp == 2 .or. ibltyp == 99 ) then
      do k = 1 , kzp1
        do i = ice1 , ice2
          do j = jce1 , jce2
            uwstatea%tkeps(j,i,k) = atm1%tke(j,i,k)*sfs%psa(j,i)
            uwstatea%advtke(j,i,k) = d_zero
          end do
        end do
      end do
      do i = ice1 , ice2
        do j = jce1 , jce2
          uwstatea%tkeps(j,i,kz+1) = atm1%tke(j,i,kz+1)*sfs%psa(j,i)
        end do
      end do
    end if
    if ( ibltyp == 2 .or. ibltyp == 99 ) then
      ! Calculate the horizontal advective tendency for TKE
      call hadvtke(uwstatea,atm1,twt,dx4)
      ! Calculate the vertical advective tendency for TKE
      call vadvtke(uwstatea,qdot,2)
      ! Calculate the horizontal, diffusive tendency for TKE
      call diffu_x(uwstatea%advtke,atms%tkeb3d,sfs%psb,xkcf,kzp1)
    end if
    !
    !   store the xxa variables in xxb and xxc in xxa:
    !   perform time smoothing operations.
    !
    do k = 1 , kz
      do i = idi1 , idi2
        do j = jdi1 , jdi2
          atm2%u(j,i,k) = omuhf*atm1%u(j,i,k)/mddom%msfd(j,i) + &
                          gnuhf*(atm2%u(j,i,k)+atmc%u(j,i,k))
          atm2%v(j,i,k) = omuhf*atm1%v(j,i,k)/mddom%msfd(j,i) + &
                          gnuhf*(atm2%v(j,i,k)+atmc%v(j,i,k))
          atm1%u(j,i,k) = atmc%u(j,i,k)
          atm1%v(j,i,k) = atmc%v(j,i,k)
          ! TAO: Once the full loop above is completed, update the TKE
          ! tendency if the UW PBL is running.  NOTE!!! Do not try to
          ! combine these loops with the above loop Advection MUST be
          ! done in a loop separate from the updates.  (I lost 3 days
          ! of working to disocover that this is a problem because I
          ! thought it would be clever to combine loops--TAO)
          if ( ibltyp == 2 .or. ibltyp == 99 ) then
            ! Add the advective tendency to the TKE tendency calculated
            ! by the UW TKE routine
             aten%tke(j,i,k) = aten%tke(j,i,k) + &
                               uwstatea%advtke(j,i,k)/sfs%psa(j,i)
             ! Do a filtered time integration
             atmc%tke(j,i,k) = max(tkemin,atm2%tke(j,i,k) + &
                               dt*aten%tke(j,i,k))
             atm2%tke(j,i,k) = max(tkemin,omuhf*atm1%tke(j,i,k) + &
                               gnuhf*(atm2%tke(j,i,k) + atmc%tke(j,i,k)))
             atm1%tke(j,i,k) = atmc%tke(j,i,k)
          end if ! TKE tendency update
        end do
      end do
    end do
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          atm2%t(j,i,k) = omuhf*atm1%t(j,i,k) + &
                          gnuhf*(atm2%t(j,i,k)+atmc%t(j,i,k))
          atm1%t(j,i,k) = atmc%t(j,i,k)
        end do
      end do
    end do
    do n = 1 , nqx
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            qxas = atmc%qx(j,i,k,n)
            qxbs = omuhf*atm1%qx(j,i,k,n) + &
                   gnuhf*(atm2%qx(j,i,k,n)+atmc%qx(j,i,k,n))
            lowq = minqx*sfs%psa(j,i)
            if ( qxas < lowq ) then
              if ( n == iqv ) then
                qxas = lowq
              else
                qxas = d_zero
              end if
            end if
            lowq = minqx*sfs%psb(j,i)
            if ( qxbs < lowq ) then
              if ( n == iqv ) then
                qxbs = lowq
              else
                qxbs = d_zero
              end if
            end if
            atm1%qx(j,i,k,n) = qxas
            atm2%qx(j,i,k,n) = qxbs
          end do
        end do
      end do
    end do
    if ( ichem == 1 ) then
      do itr = 1 , ntr
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              chias = chic(j,i,k,itr)
              if ( chias < dlowval ) chias = dlowval
              chibs = omu*chia(j,i,k,itr)                             &
                      + gnu*(chib(j,i,k,itr)+chic(j,i,k,itr))
              if ( chibs < dlowval ) chibs = dlowval
              chib(j,i,k,itr) = chibs
              chia(j,i,k,itr) = chias
            end do
          end do
        end do
      end do
    end if
    do i = ici1 , ici2
      do j = jci1 , jci2
        sfs%psb(j,i) = omuhf*sfs%psa(j,i) + gnuhf*(sfs%psb(j,i)+psc(j,i))
        sfs%psa(j,i) = psc(j,i)
      end do
    end do
    !
    ! increment elapsed forecast time:
    !
    ktau = ktau + 1
    xbctime = xbctime + dtsec
    if ( islab_ocean == 1 ) xslabtime = xslabtime + dtsec
    nbdytime = nbdytime + 1
    idatex = idatex + intmdl
    if ( mod(ktau,khour) == 0 ) then
      call split_idate(idatex,xyear,xmonth,xday,xhour)
    end if
    if ( ktau == 2 ) then
      dtbat = dt*dble(ntsrf)
      dt = dt2
      rdt = d_one/dt
      dtsq = dt*dt
      dtcb = dt*dt*dt
    end if
    !
    ! do cumulus transport/mixing  of tracers for the schemes allowing it
    !
    if ( ichem == 1 .and. ichcumtra == 1 .and. &
         ( icup == 4 .or. icup == 5 .or. icup == 97 ) ) then
      call cumtran
    end if
    !
    ! Print out noise parameter
    !
    if ( ktau > 1 ) then
      ! Added a check for nan... The following inequality is wanted.
      if ((ptntot /= ptntot) .or. &
         ((ptntot > d_zero) .eqv. (ptntot <= d_zero))) then
        maxv = dabs(maxval(aten%t))
        if ( (maxv/dtsec) > 0.01D0 ) then ! 50 K per hour
          write(stderr,*) 'MAXVAL ATEN T :', maxv
          maxv = maxv - 0.001D0
          do kk = 1 , kz
            do ii = ici1 , ici2
              do jj = jci1 , jci2
                if ( aten%t(jj,ii,kk) > maxv ) then
                  write(stderr,*) 'II :', global_cross_istart+ii-1 , &
                                ', JJ :', global_cross_jstart+jj-1 , &
                                ', KK :', kk
                end if
              end do
            end do
          end do
        end if
        maxv = dabs(maxval(aten%u))
        if ( (maxv/dtsec) > 0.005D0 ) then  ! 25 m/s per hour
          write(stderr,*) 'MAXVAL ATEN U :', maxv
          maxv = maxv - 0.001D0
          do kk = 1 , kz
            do ii = ici1 , ici2
              do jj = jci1 , jci2
                if ( aten%u(jj,ii,kk) > maxv ) then
                  write(stderr,*) 'II :', global_dot_istart+ii-1 , &
                                ', JJ :', global_dot_jstart+jj-1 , ', KK :', kk
                end if
              end do
            end do
          end do
        end if
        maxv = dabs(maxval(aten%v))
        if ( (maxv/dtsec) > 0.005D0 ) then  ! 25 m/s per hour
          write(stderr,*) 'MAXVAL ATEN V :', maxv
          maxv = maxv - 0.001D0
          do kk = 1 , kz
            do ii = ici1 , ici2
              do jj = jci1 , jci2
                if ( aten%v(jj,ii,kk) > maxv ) then
                  write(stderr,*) 'II :', global_dot_istart+ii-1 , &
                                ', JJ :', global_dot_jstart+jj-1 , ', KK :', kk
                end if
              end do
            end do
          end do
        end if
        maxv = dabs(maxval(aten%qx(:,:,:,iqv)))
        if ( (maxv/dtsec) > 0.001D0 ) then ! 
          write(stderr,*) 'MAXVAL ATEN QV :', maxv
          maxv = maxv - 0.001D0
          do kk = 1 , kz
            do ii = ici1 , ici2
              do jj = jci1 , jci2
                if ( aten%qx(jj,ii,kk,iqv) > maxv ) then
                  write(stderr,*) 'II :', global_cross_istart+ii-1 , &
                                ', JJ :', global_cross_jstart+jj-1 , &
                                ', KK :', kk
                end if
              end do
            end do
          end do
        end if
        maxv = dabs(maxval(aten%qx(:,:,:,iqc)))
        if ( (maxv/dtsec) > 0.001D0 ) then ! 
          write(stderr,*) 'MAXVAL ATEN QC :', maxv
          maxv = maxv - 0.001D0
          do kk = 1 , kz
            do ii = ici1 , ici2
              do jj = jci1 , jci2
                if ( aten%qx(jj,ii,kk,iqc) > maxv ) then
                  write(stderr,*) 'II :', global_cross_istart+ii-1 , &
                                ', JJ :', global_cross_jstart+jj-1 , &
                                ', KK :', kk
                end if
              end do
            end do
          end do
        end if
        write (*,*) 'WHUUUUBBBASAAAGASDDWD!!!!!!!!!!!!!!!!'
        write (*,*) 'No more atmosphere here....'
        write (*,*) 'CFL violation detected, so model STOP'
        write (*,*) '#####################################'
        write (*,*) '#            DECREASE DT !!!!       #'
        write (*,*) '#####################################'
        call fatal(__FILE__,__LINE__,'CFL VIOLATION')
      end if

      if ( mod(ktau,krep) == 0 ) then
        pt2bar = pt2tot
        ptnbar = ptntot
        iconvec = 0
        pt2tot = d_zero
        ptntot = d_zero
        call sumall(total_precip_points,iconvec)
        call sumall(pt2bar,pt2tot)
        call sumall(ptnbar,ptntot)
        if ( myid == italk ) then
          appdat = tochar(idatex)
          ptntot = ptntot*rptn
          pt2tot = pt2tot*rptn
          write(stdout,'(a,a23,a,i16)') ' $$$ ', appdat , ', ktau   = ', ktau
          write(stdout,'(a,2E12.5)') ' $$$ 1st, 2nd time deriv of ps   = ', &
                ptntot , pt2tot
          write(stdout,'(a,i7)') ' $$$  no. of points w/convection = ', iconvec
        end if
      end if
    end if

#ifdef DEBUG
    call time_end(subroutine_name,idindx)

    contains

    ! Check temperature tendency less than 10 K

    subroutine check_temperature_tendency(loc)
      implicit none
      character(len=*) , intent(in) :: loc
      integer(ik4) :: i , j , k , kk , ierr
      real(rk8) :: check_tt , mean_tt
      ierr = 0
      mean_tt = (maxval(aten%t)+minval(aten%t))/d_two
      do k = 1 , kz
        do i = ice1, ice2
          do j = jce1 , jce2
            check_tt = (aten%t(j,i,k)-mean_tt)/psc(j,i)
            if ( dabs(check_tt) > temp_tend_maxval ) then
              write(stderr,*) 'After ', loc, ' at ktau = ', ktau
              write(stderr,*) 'TEMP tendency out of order : ', check_tt
              write(stderr,*) 'At J = ',global_dot_jstart+j
              write(stderr,*) 'At I = ',global_dot_istart+i
              write(stderr,*) 'At K = ',k
              write(stderr,*) 'Surface Temperature : ', sfs%tga(j,i)
              write(stderr,*) 'Vertical PTU profile: '
              do kk = 1 , kz
                write(stderr,'(i2,3f12.7)') kk, &
                        atms%pb3d(j,i,kk), atms%tb3d(j,i,kk), &
                        sqrt(atms%ubx3d(j,i,kk)**2+ &
                             atms%vbx3d(j,i,kk)**2)
              end do
              ierr = ierr + 1
            end if
          end do
        end do
      end do
      if ( ierr /= 0 ) then
        call fatal(__FILE__,__LINE__,'TEMP TENDENCY ERROR')
      end if
    end subroutine check_temperature_tendency

    ! Check wind speed tendency less than 10 m/s

    subroutine check_wind_tendency(loc)
      implicit none
      character(len=*) , intent(in) :: loc
      integer(ik4) :: i , j , k , kk , ierr
      real(rk8) :: check_ww , mean_ww
      ierr = 0
      wten = sqrt(aten%u**2+aten%v**2)
      mean_ww = (maxval(wten)+minval(wten))/d_two
      do k = 1 , kz
        do i = ide1, ide2
          do j = jde1 , jde2
            check_ww = (wten(j,i,k)-mean_ww)/psdot(j,i)
            if ( dabs(check_ww) > wind_tend_maxval ) then
              write(stderr,*) 'After ', loc, ' at ktau = ', ktau
              write(stderr,*) 'WIND tendency out of order : ', check_ww
              write(stderr,*) 'At J = ',global_dot_jstart+j
              write(stderr,*) 'At I = ',global_dot_istart+i
              write(stderr,*) 'At K = ',k
              write(stderr,*) 'Surface Temperature : ', sfs%tga(j,i)
              write(stderr,*) 'Vertical PTU profile: '
              do kk = 1 , kz
                write(stderr,'(i2,3f12.7)') kk, &
                        atms%pb3d(j,i,kk), atms%tb3d(j,i,kk), &
                        sqrt(atms%ubx3d(j,i,kk)**2+ &
                             atms%vbx3d(j,i,kk)**2)
              end do
              ierr = ierr + 1
            end if
          end do
        end do
      end do
      if ( ierr /= 0 ) then
        call fatal(__FILE__,__LINE__,'WIND TENDENCY ERROR')
      end if
    end subroutine check_wind_tendency
#endif

  end subroutine tend
!
end module mod_tendency

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
  use mod_bdycod
  use mod_precip
  use mod_slice
  use mod_sun
  use mod_advection
  use mod_diffusion
  use mod_domain
  use mod_cloud_s1
  use mod_sladvection
  use mod_slabocean
  use mod_sound
  use mod_split
  use mod_timefilter

  implicit none

  private

  public :: allocate_mod_tend , tend
  real(rkx) , pointer , dimension(:,:,:) :: ttld , td , &
         phi , ten0 , qen0 , qcd , qvd , tvfac , ucc , vcc , th , tha
  real(rkx) , pointer , dimension(:,:,:) :: ps4
  real(rkx) , pointer , dimension(:,:,:) :: ps_4
  real(rkx) , pointer , dimension(:,:) :: pten
  real(rkx) , pointer , dimension(:,:,:) :: thten
  real(rkx) , pointer , dimension(:,:) :: dummy , rpsa , rpsb , rpsc
  real(rkx) , pointer , dimension(:,:) :: rpsda

  integer :: ithadv = 1
  integer(ik4) :: iqxvadv , itrvadv
#ifdef DEBUG
  real(rkx) , pointer , dimension(:,:,:) :: wten
#endif

  real(rkx) :: rptn ! Total number of internal points

  contains

#include <cpmf.inc>

  subroutine allocate_mod_tend
    implicit none
    call getmem3d(ps_4,jcross1,jcross2,icross1,icross2,1,4,'tendency:ps_4')
    call getmem3d(ps4,jci1,jci2,ici1,ici2,1,4,'tendency:ps4')
    if ( ipptls == 2 ) then
      call getmem3d(qcd,jce1,jce2,ice1,ice2,1,kz,'tendency:qcd')
    else
      call assignpnt(atmx%qx,qcd,iqc)
    end if
    call getmem3d(tvfac,jce1,jce2,ice1,ice2,1,kz,'tendency:tvfac')
    call assignpnt(atmx%qx,qvd,iqv)
    call getmem2d(pten,jce1,jce2,ice1,ice2,'tendency:pten')
    call getmem2d(dummy,jde1,jde2,ide1,ide2,'tendency:dummy')
    call getmem2d(rpsa,jce1ga,jce2ga,ice1ga,ice2ga,'tendency:rpsa')
    call getmem2d(rpsb,jce1,jce2,ice1,ice2,'tendency:rpsb')
    call getmem2d(rpsda,jde1ga,jde2ga,ide1ga,ide2ga,'tendency:rpsda')
    if ( idynamic == 1 ) then
      ithadv = 0
      call getmem3d(ttld,jce1,jce2,ice1,ice2,1,kz,'tend:ttld')
      call getmem3d(td,jce1,jce2,ice1,ice2,1,kz,'tendency:td')
      call getmem3d(phi,jce1ga,jce2,ice1ga,ice2,1,kz,'tendency:phi')
      call getmem2d(rpsc,jce1,jce2,ice1,ice2,'tendency:rpsc')
      rptn = d_one/real((jout2-jout1+1)*(iout2-iout1+1),rkx)
    else if ( idynamic == 2 ) then
      call getmem3d(ucc,jce1,jce2,ice1,ice2,1,kz,'tendency:ucc')
      call getmem3d(vcc,jce1,jce2,ice1,ice2,1,kz,'tendency:vcc')
      if ( ithadv == 1 ) then
        call getmem3d(thten,jce1,jce2,ice1,ice2,1,kz,'tendency:thten')
        call getmem3d(th,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'tendency:th')
        call getmem3d(tha,jce1,jce2,ice1,ice2,1,kz,'tendency:tha')
      end if
    end if
    if ( idiag > 0 ) then
      call getmem3d(ten0,jce1,jce2,ice1,ice2,1,kz,'tendency:ten0')
      call getmem3d(qen0,jce1,jce2,ice1,ice2,1,kz,'tendency:qen0')
    end if
#ifdef DEBUG
    call getmem3d(wten,jde1,jde2,ide1,ide2,1,kz,'tendency:wten')
#endif
    !
    ! Set number of ghost points for advection for the two schemes
    ! Select advection scheme
    !
    iqxvadv = 1
    itrvadv = 2
    if ( ibltyp == 2 ) then
      if ( iuwvadv == 1 ) then
        itrvadv = 3
        iqxvadv = 3
      end if
    end if
  end subroutine allocate_mod_tend
  !
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !                                                                     c
  ! This subroutine computes the tendencies of the prognostic           c
  ! variables                                                           c
  !                                                                     c
  !cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
  !
  subroutine tend
    implicit none
    real(rkx) :: duv , pt2bar , pt2tot , ptnbar , maxv , ptntot ,  &
                 rovcpm , rtbar , tv , tva , tvavg , tvb , tvc ,   &
                 rho0s , cpm , rofac , uaq , vaq , wabar , amfac , &
                 wadot , wadotp1
    integer(ik4) :: i , itr , j , k , lev , n , ii , jj , kk , iconvec
    logical :: loutrad , labsem
    character (len=32) :: appdat
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'tend'
    integer(ik4) , save :: idindx = 0
    integer(ik4) :: mijx , miiy
    call time_begin(subroutine_name,idindx)
    mijx = (jci2-jci1)/2+1
    miiy = (ici2-ici1)/2+1
#endif
    !
    ! Prepare pressures on two timesteps, on both cross and dot
    ! points, and compute 1/ps.
    !
    call surface_pressures
    !
    ! Decoupling on atmx and apply bdy conditions to U,V
    !
    call decouple
    !
    if ( ipptls == 2 ) then
      qcd(:,:,:) = d_zero
      do n = iqfrst , iqlst
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              qcd(j,i,k) = qcd(j,i,k) + atmx%qx(j,i,k,n)
            end do
          end do
        end do
      end do
    end if
    !
    ! compute vertical sigma-velocity (qdot):
    !
    qdot(:,:,:)  = d_zero
    do i = ice1 , ice2
      do j = jce1 , jce2
        dummy(j,i) = d_one/(dx2*mddom%msfx(j,i)*mddom%msfx(j,i))
      end do
    end do
    if ( idynamic == 1 ) then
      !
      ! compute the pressure tendency
      !
      pten(:,:) = d_zero
      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            !
            ! The surface pressure tendency in the   hydrostatic model:
            ! Eq. 2.1.5 & Eq. 2.4.2 in the MM5 manual
            !
            mdv%cr(j,i,k) = ((atmx%umc(j+1,i+1,k)+atmx%umc(j+1,i,k)- &
                              atmx%umc(j,i+1,k)  -atmx%umc(j,i,k)) + &
                             (atmx%vmc(j+1,i+1,k)+atmx%vmc(j,i+1,k)- &
                              atmx%vmc(j+1,i,k)  -atmx%vmc(j,i,k)))*dummy(j,i)
            pten(j,i) = pten(j,i) - mdv%cr(j,i,k) * dsigma(k)
          end do
        end do
      end do
      do k = 2 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            !
            ! The coordinate vertical velocity in the   hydrostatic model:
            ! Eq. 2.1.6 & Eq. 2.4.3 in the MM5 manual
            !
            qdot(j,i,k) = qdot(j,i,k-1) - (pten(j,i) + &
                       mdv%cr(j,i,k-1)) * dsigma(k-1) * rpsa(j,i)
           end do
        end do
      end do
    else if ( idynamic == 2 ) then
      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            !
            ! Calculate wind components at cross points
            !
            ucc(j,i,k) = (atmx%umd(j,i,k)  + atmx%umd(j,i+1,k) + &
                          atmx%umd(j+1,i,k)+ atmx%umd(j+1,i+1,k))
            vcc(j,i,k) = (atmx%vmd(j,i,k)  + atmx%vmd(j,i+1,k) + &
                          atmx%vmd(j+1,i,k)+ atmx%vmd(j+1,i+1,k))
          end do
        end do
      end do
      do k = 2 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            !
            ! The coordinate vertical velocity in the nonhydrostatic model:
            ! Eq. 2.2.7 & Eq. 2.3.6 in the MM5 manual
            !
            rho0s = twt(k,1)*atm0%rho(j,i,k)+twt(k,2)*atm0%rho(j,i,k-1)
            qdot(j,i,k) = -rho0s*egrav*atmx%w(j,i,k)/atm0%ps(j,i) -  &
              sigma(k) * (dpsdxm(j,i) * (twt(k,1)*ucc(j,i,k) +       &
                                         twt(k,2)*ucc(j,i,k-1)) +    &
                          dpsdym(j,i) * (twt(k,1)*vcc(j,i,k) +       &
                                         twt(k,2)*vcc(j,i,k-1)))
           end do
        end do
      end do
      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            !
            ! The mass divergence term (in cross points) in the
            ! nonhydrostatic model:
            ! Eq. 2.2.6 & Eq. 2.3.5 in the MM5 manual
            !
            mdv%cr(j,i,k) = ((atmx%umc(j+1,i+1,k)+atmx%umc(j+1,i,k) -  &
                              atmx%umc(j,i+1,k)  -atmx%umc(j,i,k))  +  &
                             (atmx%vmc(j+1,i+1,k)+atmx%vmc(j,i+1,k) -  &
                              atmx%vmc(j+1,i,k)  -atmx%vmc(j,i,k))) *  &
                     dummy(j,i) + (qdot(j,i,k+1) - qdot(j,i,k)) *      &
                        sfs%psa(j,i)/dsigma(k)
          end do
        end do
      end do
    end if
    call exchange(mdv%cr,1,jce1,jce2,ice1,ice2,1,kz)
    call exchange(qdot,1,jce1,jce2,ice1,ice2,1,kzp1)
    !
    ! compute omega
    !
    omega(:,:,:) = d_zero
    if ( idynamic == 1 ) then
      do i = ici1 , ici2
        do j = jci1 , jci2
          dummy(j,i) = d_one/(dx8*mddom%msfx(j,i))
        end do
      end do
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            !
            ! omega in the hydrostatic model: Eqs. 2.1.7, 2.1.8 & 2.4.4
            !
            omega(j,i,k) = d_half*(qdot(j,i,k+1)+qdot(j,i,k)) *   &
                        sfs%psa(j,i) + hsigma(k) * (pten(j,i) +   &
                       ((atmx%ud(j,i,k) + atmx%ud(j,i+1,k) +      &
                         atmx%ud(j+1,i+1,k) + atmx%ud(j+1,i,k))*  &
                         (sfs%psa(j+1,i)-sfs%psa(j-1,i)) +        &
                        (atmx%vd(j,i,k) + atmx%vd(j,i+1,k) +      &
                         atmx%vd(j+1,i+1,k) + atmx%vd(j+1,i,k)) * &
                         (sfs%psa(j,i+1)-sfs%psa(j,i-1)))*dummy(j,i))
          end do
        end do
      end do
    else if ( idynamic == 2 ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            !
            ! omega in the non-hydrostatic model: compute from w
            !
            omega(j,i,k) = d_half*egrav*atm0%rho(j,i,k)*rpsb(j,i) * &
                         (atm2%w(j,i,k)+atm2%w(j,i,k+1))
          end do
        end do
      end do
    end if
    !
    ! Prepare fields to be used in physical parametrizations.
    !
    call mkslice
    !
    ! Surface pressure boundary conditions for Hydrostatic
    !
    if ( idynamic == 1 ) then
      if ( iboudy == 4 ) then
        call sponge(xpsb,pten)
      else if ( iboudy == 1 .or. iboudy == 5 ) then
        call nudge(iboudy,sfs%psb,xpsb,pten)
      end if
      !
      ! psc : forecast pressure
      !
      do i = ice1 , ice2
        do j = jce1 , jce2
          sfs%psc(j,i) = sfs%psb(j,i) + pten(j,i)*dt
          rpsc(j,i) = d_one/sfs%psc(j,i)
        end do
      end do
      !
      ! compute bleck (1977) noise parameters:
      !
      if ( ktau /= 0 ) then
        ptntot = d_zero
        pt2tot = d_zero
        do i = ici1 , ici2
          do j = jci1 , jci2
            ptntot = ptntot + abs(pten(j,i))
            pt2tot = pt2tot + abs((sfs%psc(j,i)+sfs%psb(j,i)- &
                     d_two*sfs%psa(j,i))/(dt*dt*d_rfour))
          end do
        end do
      end if
    end if
    !
    ! calculate solar zenith angle
    !
    call zenitm(coszrs)
    !
    ! Compute new diffusion coefficients
    !
    call calc_coeff( )
    !
    ! Initialize the tendencies
    !
    aten%u(:,:,:) = d_zero
    aten%v(:,:,:) = d_zero
    aten%t(:,:,:) = d_zero
    aten%qx(:,:,:,:) = d_zero
    if ( idynamic == 2 ) then
      !
      ! Pressure perturbations and vertical velocity tendencies in the
      ! nonhydrostatic model
      !
      aten%pp(:,:,:) = d_zero
      aten%w(:,:,:) = d_zero
    end if
#ifdef DEBUG
    wten(:,:,:) = d_zero
#endif
    cldfra(:,:,:) = d_zero
    cldlwc(:,:,:) = d_zero
    !
    ! Initialize diffusion terms (temperature, vertical velocity, mixing ratios)
    !
    adf%t(:,:,:)    = d_zero
    adf%qx(:,:,:,:) = d_zero
    adf%u(:,:,:) = d_zero
    adf%v(:,:,:) = d_zero
    if ( idynamic == 2 ) then
      adf%w(:,:,:)    = d_zero
      adf%pp(:,:,:)   = d_zero
    end if
    if ( ibltyp == 2 ) then
      aten%tke(:,:,:) = d_zero
      uwstatea%advtke(:,:,:) = d_zero
    end if
    if ( ichem == 1 ) then
      chiten(:,:,:,:)  = d_zero
      if ( ichdiag == 1 ) chiten0(:,:,:,:) = d_zero
    end if
    if ( idiag > 0 ) then
      ten0(:,:,:) = d_zero
      qen0(:,:,:) = d_zero
    end if
    !
    ! Compute Horizontal advection terms
    !
    call start_advect

    if ( idynamic == 2 ) then
      !
      ! Horizontal and vertical advection of pressure perturbation and vertical
      ! velocity in the nonhydrostatic model: 1st and 2nd term on the RHS of the
      ! Eq. 2.2.3, Eq. 2.2.4, Eq. 2.3.7 and Eq. 2.3.8 in the MM5 manual.
      !
      ! Also, cf. Eq. 2.2.11 of vertical velocity tendency in the MM5 manual.
      !
      call hadv(aten%pp,atmx%pp,0)
      call vadv(aten%pp,atm1%pp,kz,0)
      call hadv(aten%w,atmx%w,1)
      call vadv(aten%w,atm1%w,kzp1,0)
      if ( iboudy == 1 .or. iboudy == 5 ) then
        call nudge(iboudy,atm2%t,xtb,aten%t)
        call nudge(iboudy,atm2%qx,xqb,aten%qx,iqv)
        call nudge(iboudy,atm2%pp,xppb,aten%pp)
        call nudge(iboudy,atm2%w,xwwb,aten%w)
      else if ( iboudy == 4 ) then
        call sponge(xtb,aten%t)
        call sponge(xqb,aten%qx,iqv)
        call sponge(xppb,aten%pp)
        call sponge(xwwb,aten%w)
      end if
      if ( idiag > 0 ) then
        ! rq : temp condensation tend is added the evap temp tend
        !      calculated in pcp
        tdiag%bdy = tdiag%bdy + (aten%t - ten0) * afdout
        ten0 = aten%t
      end if
    end if
    !
    ! compute the horizontal advection term in temperature tendency:
    ! same for hydrostatic and nonhydrostatic models:
    ! in Eqs. 2.1.3, 2.2.5, 2.3.9 (1st RHS term)
    !
    if ( ithadv == 0 ) then
      call hadv(aten%t,atmx%t)
      if ( idiag > 0 ) then
        tdiag%adh = tdiag%adh + (aten%t - ten0) * afdout
        ten0 = aten%t
      end if
#ifdef DEBUG
      call check_temperature_tendency('HADV')
#endif
      !
      ! compute the vertical advection term:
      ! same for hydrostatic and nonhydrostatic models:
      ! in Eqs. 2.1.3, 2.2.5, 2.3.9 (2nd RHS term)
      !
      call vadv(aten%t,atm1%t,kz,1)
      if ( idiag > 0 ) then
        tdiag%adv = tdiag%adv + (aten%t - ten0) * afdout
        ten0 = aten%t
      end if
#ifdef DEBUG
      call check_temperature_tendency('VADV')
#endif
    else
      thten(:,:,:) = d_zero
      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            th(j,i,k) = atmx%t(j,i,k) * (1.0e5_rkx/atm1%pr(j,i,k))**rovcp
            tha(j,i,k) = th(j,i,k) * sfs%psa(j,i)
          end do
        end do
      end do
      call exchange(th,1,jce1,jce2,ice1,ice2,1,kz)
      call hadv(thten,th)
      if ( idiag > 0 ) then
        tdiag%adh = tdiag%adh + (thten - ten0) * afdout
        ten0 = thten
      end if
#ifdef DEBUG
      call check_temperature_tendency('HADV')
#endif
      call vadv(thten,tha,kz,0)
      if ( idiag > 0 ) then
        tdiag%adv = tdiag%adv + (thten - ten0) * afdout
        ten0 = thten
      end if
#ifdef DEBUG
      call check_temperature_tendency('VADV')
#endif
    end if
    !
    ! compute the adiabatic term:
    !
    if ( idynamic == 1 ) then
      !
      ! Adiabatic term in the temperature tendency equation in the
      ! hydrostatic model:    3rd RHS term in Eq. 2.1.3
      !
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            rovcpm = rgas/cpmf(qvd(j,i,k))
            aten%t(j,i,k) = aten%t(j,i,k) + &
                            (omega(j,i,k)*rovcpm*atmx%tv(j,i,k)) / &
                            (ptop*rpsa(j,i)+hsigma(k))
          end do
        end do
      end do
    else if ( idynamic == 2 ) then
      !
      ! Adiabatic term in the temperature tendency equation in the
      ! nonhydrostatic model: 3rd and 4th RHS term in Eq. 2.2.5 and Eq.2.3.9.
      !
      if ( ithadv == 0 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              cpm = cpmf(qvd(j,i,k))
              aten%t(j,i,k) = aten%t(j,i,k) + atmx%t(j,i,k)*mdv%cr(j,i,k) - &
                        (omega(j,i,k)*sfs%psa(j,i) + aten%pp(j,i,k) + &
                         atmx%pp(j,i,k)*mdv%cr(j,i,k)) / (atm1%rho(j,i,k)*cpm)
            end do
          end do
        end do
      else
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              thten(j,i,k) = thten(j,i,k) + th(j,i,k) * mdv%cr(j,i,k)
              aten%t(j,i,k) = aten%t(j,i,k) + &
                   atm1%t(j,i,k)*thten(j,i,k)/tha(j,i,k)
            end do
          end do
        end do
      end if
      !
      ! Divergence term in the pressure perturbation tendency equation in the
      ! nonhydrostatic model: 4th RHS term in Eq. 2.2.4 and Eq. 2.3.8
      !
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            aten%pp(j,i,k) = aten%pp(j,i,k) + atmx%pp(j,i,k)*mdv%cr(j,i,k)
          end do
        end do
      end do
      do n = 1 , nqx
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              aten%qx(j,i,k,n) = aten%qx(j,i,k,n) + &
                                 atmx%qx(j,i,k,n)*mdv%cr(j,i,k)
            end do
          end do
        end do
      end do
      !
      ! Vertical velocity tendency. Following terms are included here:
      ! (1) bouyancy terms: 2nd subterm and part of the 3rd subterm of the
      !     4th RHS term in Eq.2.2.3 and 2.2.11. This is joined into the 5th
      !     RHS term in Eq. 2.3.7.
      ! (2) part of the vertical component of the Coriolis force due to the
      !     horizontal movement (cf. 6th RHS term in Eq. 2.2.11)
      ! (3) vertical curvature term (not explicitly mentioned in the MM5 1994
      !     manual)
      ! (4) mass divergence term (3rd RHS term in Eq. 2.2.3, 2.2.11 and
      !     Eq. 2.3.7)
      !
      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            !
            ! Calculate wind components at cross points
            !
            ucc(j,i,k) = (atmx%uc(j,i,k)  + atmx%uc(j,i+1,k) + &
                          atmx%uc(j+1,i,k)+ atmx%uc(j+1,i+1,k))
            vcc(j,i,k) = (atmx%vc(j,i,k)  + atmx%vc(j,i+1,k) + &
                          atmx%vc(j+1,i,k)+ atmx%vc(j+1,i+1,k))

          end do
        end do
      end do
      do k = 2 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            rofac = ( dsigma(k-1) * atm0%rho(j,i,k) +      &
                      dsigma(k)   * atm0%rho(j,i,k-1) ) /  &
                    ( dsigma(k-1) * atm1%rho(j,i,k) +      &
                      dsigma(k)   * atm1%rho(j,i,k-1) )
            uaq = d_rfour * (twt(k,1) * ucc(j,i,k) + &
                             twt(k,2) * ucc(j,i,k-1))
            vaq = d_rfour * (twt(k,1) * vcc(j,i,k) + &
                             twt(k,2) * vcc(j,i,k-1))
            aten%w(j,i,k) = aten%w(j,i,k) +                &
                  (twt(k,2)*atmx%pr(j,i,k-1) +             &
                   twt(k,1)*atmx%pr(j,i,k)) *              &
                   rofac * egrav * sfs%psa(j,i) +          &
                   mddom%ex(j,i)*(uaq*mddom%crx(j,i) -     &
                                  vaq*mddom%cry(j,i)) +    &
                   (uaq*uaq+vaq*vaq)*rearthrad*rpsa(j,i) + &
                   atmx%w(j,i,k)*(twt(k,1)*mdv%cr(j,i,k) + &
                                  twt(k,2)*mdv%cr(j,i,k-1))
          end do
        end do
      end do
      if ( ipptls > 0 ) then
        do k = 2 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              !
              ! Vertical velocity tendency: water loading term
              ! 5th RHS term in Eq. 2.2.3 & 6th RHS term in Eq. 2.3.7
              !
              aten%w(j,i,k) = aten%w(j,i,k) - egrav * sfs%psa(j,i) * &
                (twt(k,2)*qcd(j,i,k-1) + twt(k,1)*qcd(j,i,k))
            end do
          end do
        end do
      end if
    end if
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
    if ( idiag > 0 ) ten0(jci1:jci2,ici1:ici2,:) = adf%t
    call diffu_x(adf%t,atms%tb3d)
    if ( idiag > 0 ) then
      ! save the h diff diag here
      tdiag%dif(jci1:jci2,ici1:ici2,:) = tdiag%dif(jci1:jci2,ici1:ici2,:) + &
        (adf%t - ten0(jci1:jci2,ici1:ici2,:)) * afdout
      ! reset ten0 to aten%t
      ten0 = aten%t
    end if
    !
    ! compute the diffusion term for vertical velocity w and store in diffw:
    ! compute the diffusion term for perturb pressure pp and store in diffpp:
    !
    if ( idynamic == 2 ) then
      call diffu_x(adf%pp,atms%ppb3d)
      call diffu_x(adf%w,atms%wb3d,d_one)
    end if
    !
    ! compute the moisture tendencies for convection
    !
    if ( isladvec == 1 ) then
      call slhadv_x(aten%qx,atm2%qx,iqv)
      call hdvg_x(aten%qx,atm1%qx,iqv)
    else
      call hadv(aten%qx,atmx%qx,iqv)
    end if
    if ( idiag > 0 ) then
      qdiag%adh = qdiag%adh + (aten%qx(:,:,:,iqv) - qen0) * afdout
      qen0 = aten%qx(:,:,:,iqv)
    end if
    if ( all(icup /= 1) ) then
      call vadv(aten%qx,atm1%qx,iqv)
    end if
    if ( idiag > 0 ) then
      qdiag%adv = qdiag%adv + (aten%qx(:,:,:,iqv) - qen0) * afdout
      qen0 = aten%qx(:,:,:,iqv)
    end if
    if ( ipptls > 0 ) then
      if ( isladvec == 1 ) then
        call slhadv_x(aten%qx,atm2%qx,iqfrst,iqlst)
        call hdvg_x(aten%qx,atm1%qx,iqfrst,iqlst)
      else
        call hadv(aten%qx,atmx%qx,iqfrst,iqlst)
      end if
      call vadv(aten%qx,atm1%qx,iqfrst,iqlst,iqxvadv)
    end if
    if ( ichem == 1 ) then
      !
      ! horizontal and vertical advection + diag
      !
      if ( ichdiag == 1 ) chiten0 = chiten
      if ( isladvec == 1 ) then
        call slhadv_x(chiten,chib)
        call hdvg_x(chiten,chia)
      else
        call hadv(chiten,chi)
      end if
      if ( ichdiag == 1 ) then
        cadvhdiag = cadvhdiag + (chiten - chiten0) * cfdout
        chiten0 = chiten
      end if
      if ( all(icup /= 1) ) then
        call vadv(chiten,chia,1,ntr,itrvadv)
      end if
      if ( ichdiag == 1 ) then
        cadvvdiag = cadvvdiag + (chiten - chiten0) * cfdout
        chiten0 = chiten
      end if
    end if
    !
    ! conv tracer diagnostic
    !
    if ( ichem == 1 .and. ichdiag == 1 ) chiten0 = chiten
    !
    ! Call cumulus parametrization
    !
    call cumulus
    if ( idiag > 0 ) then
      tdiag%con = tdiag%con + (aten%t - ten0) * afdout
      ten0 = aten%t
      qdiag%con = qdiag%con + (aten%qx(:,:,:,iqv) - qen0) * afdout
      qen0 = aten%qx(:,:,:,iqv)
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
      !
      ! Clouds and large scale precipitation
      !
      call cldfrac
      if ( ipptls == 2 ) then
        call microphys
      else
        call pcp
      end if
      !
      ! compute the diffusion terms:
      ! the diffusion term for qx is stored in diffqx.
      !
      if ( idiag > 0 ) then
        qen0(jci1:jci2,ici1:ici2,:) = adf%qx(jci1:jci2,ici1:ici2,:,iqv)
      end if
      if ( ipptls == 1 ) then
        call diffu_x(adf%qx,atms%qxb3d,iqfrst,iqlst,d_one)
      else
        call diffu_x(adf%qx,atms%qxb3d,iqv,d_one)
      end if
      if ( idiag > 0 ) then
        ! save the h diff diag here
        qdiag%dif(jci1:jci2,ici1:ici2,:) = &
                      qdiag%dif(jci1:jci2,ici1:ici2,:) +   &
                     (adf%qx(jci1:jci2,ici1:ici2,:,iqv) -  &
                      qen0(jci1:jci2,ici1:ici2,:)) * afdout
        ! reset qen0 to aten%t
        qen0 = aten%qx(:,:,:,iqv)
      end if
    end if

    if ( idiag > 0 ) then
      ! save tten from pcp (evaporation)
      tdiag%lsc = tdiag%lsc + (aten%t - ten0) * afdout
      ten0 = aten%t
      qdiag%lsc = qdiag%lsc + (aten%qx(:,:,:,iqv) - qen0) * afdout
      qen0 = aten%qx(:,:,:,iqv)
    end if
#ifdef DEBUG
    call check_temperature_tendency('PREC')
#endif
    !
    ! Tracers tendencies
    !
    if ( ichem == 1 ) then
      if ( ichdiag == 1 ) chiten0 = chiten
      !
      ! horizontal diffusion: initialize scratch vars to 0.
      ! need to compute tracer tendencies due to diffusion
      ! Only active if upstream scheme not used.
      !
      call diffu_x(chiten,atms%chib3d,1,ntr,d_one)
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
      call surface_albedo
      ! Update / init Ozone profiles
      if ( iclimao3 == 1 ) then
        call updateo3(idatex,scenario)
      else
        if ( ktau == 0 ) call inito3
      end if
      loutrad = (ktau == 0 .or. mod(ktau+1,krad) == 0)
      labsem = ( ktau == 0 .or. mod(ktau+1,ntabem) == 0 )
      call radiation(xyear,loutrad,labsem)
    end if

    if ( ktau == 0 .or. mod(ktau+1,ntsrf) == 0 ) then
      call surface_model
      if ( islab_ocean == 1 ) call update_slabocean(xslabtime)
    end if
    !
    ! Call medium resolution PBL
    !
    if ( ichem == 1 .and. ichdiag == 1 ) chiten0 = chiten
    ! care : pbl update the difft table at this level
    if ( idiag > 0 ) ten0(jci1:jci2,ici1:ici2,:) = adf%t
    call pblscheme
    if ( ichem == 1 .and. ichdiag == 1 ) then
      ctbldiag = ctbldiag + (chiten - chiten0) * cfdout
    end if
    if ( idiag > 0 ) then
      tdiag%tbl(jci1:jci2,ici1:ici2,:) = tdiag%tbl(jci1:jci2,ici1:ici2,:) + &
                     (adf%t - ten0(jci1:jci2,ici1:ici2,:)) * afdout
      ten0 = aten%t
      qdiag%tbl(jci1:jci2,ici1:ici2,:) = qdiag%tbl(jci1:jci2,ici1:ici2,:) + &
                (adf%qx(jci1:jci2,ici1:ici2,:,iqv) - &
            qen0(jci1:jci2,ici1:ici2,:)) * afdout
      qen0 = aten%qx(:,:,:,iqv)
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
    if ( idynamic == 1 ) then
      if ( iboudy == 4 ) then
        call sponge(xtb,aten%t)
        call sponge(xqb,aten%qx,iqv)
        if ( idiag > 0 ) then
          ! rq : temp condensation tend is added the evap temp tend
          !      calculated in pcp
          tdiag%bdy = tdiag%bdy + (aten%t - ten0) * afdout
          ten0 = aten%t
        end if
      end if
    end if
    !
    ! add horizontal diffusion and pbl tendencies for t and qv
    ! This is the last RHS term in Eqs. 2.1.3 and 2.2.5, 2.3.9
    !
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          aten%t(j,i,k) = aten%t(j,i,k) + adf%t(j,i,k)
        end do
      end do
    end do
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          aten%qx(j,i,k,iqv) = aten%qx(j,i,k,iqv) + adf%qx(j,i,k,iqv)
        end do
      end do
    end do
    !
    ! add horizontal diffusion for vertical velocity w
    ! This is the last RHS term in Eqs. 2.2.3, 2.2.11, 2.3.7
    !
    if ( idynamic == 2 ) then
      do k = 1 , kzp1
        do i = ici1 , ici2
          do j = jci1 , jci2
            aten%w(j,i,k) = aten%w(j,i,k) + adf%w(j,i,k)
          end do
        end do
      end do
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            aten%pp(j,i,k) = aten%pp(j,i,k) + adf%pp(j,i,k)
          end do
        end do
      end do
    end if
    ! Rq: the temp tendency diagnostics have been already
    !     saved for tbl and hor. diff but :
    if ( idiag > 0 ) then
      ten0 = aten%t ! important since aten%t have been updated
      qen0 = aten%qx(:,:,:,iqv) ! important since aten%qx have been updated
    end if
    !
    ! apply the nudging boundary conditions:
    !
    if ( idynamic == 1 ) then
      if ( iboudy == 1 .or. iboudy == 5 ) then
        call nudge(iboudy,atm2%t,xtb,aten%t)
        call nudge(iboudy,atm2%qx,xqb,aten%qx,iqv)
      end if
    end if
    if ( ichem == 1 ) then
      if ( ichdiag == 1 ) chiten0 = chiten
      ! keep nudge_chi for now
      if ( iboudy == 1 .or. iboudy == 5 ) then
        call nudge_chi(kz,chib,chiten)
      end if
      if ( ichdiag == 1 ) cbdydiag = cbdydiag + (chiten - chiten0) * cfdout
    end if
    if ( idiag > 0 ) then
      tdiag%bdy = tdiag%bdy + (aten%t - ten0) * afdout
      ten0 = aten%t
      qdiag%bdy = qdiag%bdy + (aten%qx(:,:,:,iqv) - qen0) * afdout
      qen0 = aten%qx(:,:,:,iqv)
    end if
#ifdef DEBUG
    call check_temperature_tendency('BDYC')
#endif
    ! Rq: the temp tendency diagnostics have been already
    !     saved for tbl and hor. diff but :
    if ( idiag > 0 ) then
      ten0 = aten%t ! important since aten%t have been updated
      qen0 = aten%qx(:,:,:,iqv) ! important since aten%qx have been updated
    end if
    !
    ! compute the condensation and precipitation terms for explicit
    ! moisture schemes
    !
    if ( ipptls > 0 ) then
      do n = iqfrst , iqlst
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              aten%qx(j,i,k,n) = aten%qx(j,i,k,n) + adf%qx(j,i,k,n)
            end do
          end do
        end do
      end do
      if ( ipptls == 1 ) then
        call condtq
      end if
      if ( idiag > 0 ) then
        ! rq : temp condensation tend is added the evap temp tend
        ! calculated in pcp
        tdiag%lsc = tdiag%lsc + (aten%t - ten0) * afdout
        ten0 = aten%t
        qdiag%lsc = qdiag%lsc + (aten%qx(:,:,:,iqv) - qen0) * afdout
        qen0 = aten%qx(:,:,:,iqv)
      end if
    end if
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
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          atmc%qx(j,i,k,iqv) = &
                  max(atm2%qx(j,i,k,iqv) + dt*aten%qx(j,i,k,iqv),minqv)
        end do
      end do
    end do
    do n = iqfrst , iqlst
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            atmc%qx(j,i,k,n) = &
                    max(atm2%qx(j,i,k,n) + dt*aten%qx(j,i,k,n),minqx)
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
              chic(j,i,k,itr) = max(chic(j,i,k,itr),mintr)
            end do
          end do
        end do
      end do
    end if
    if ( idynamic == 1 ) then
      !
      ! compute weighted p*t (td) for use in ssi:
      !
      if ( ipgf == 1 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              tva = atm1%t(j,i,k)*(d_one+ep1*qvd(j,i,k))
              tvb = atm2%t(j,i,k)*(d_one+ep1*atm2%qx(j,i,k,iqv)*rpsb(j,i))
              tvc = atmc%t(j,i,k)*(d_one+ep1*atmc%qx(j,i,k,iqv)*rpsc(j,i))
              td(j,i,k) = alpha_hyd*(tvc+tvb) + beta_hyd*tva
              ttld(j,i,k) = td(j,i,k) - sfs%psa(j,i) * &
                        t00pg*((hsigma(k)*sfs%psa(j,i)+ptop)/p00pg)**pgfaa1
            end do
          end do
        end do
        if ( ma%has_bdyleft ) then
          do k = 1 , kz
            do i = ici1 , ici2
              td(jce1,i,k) = atm1%t(jce1,i,k)*(d_one + ep1*qvd(jce1,i,k))
              ttld(jce1,i,k) = td(jce1,i,k) - sfs%psa(jce1,i) * &
                          t00pg*((hsigma(k)*sfs%psa(jce1,i)+ptop)/p00pg)**pgfaa1
            end do
        end do
        end if
        if ( ma%has_bdyright ) then
          do k = 1 , kz
            do i = ici1 , ici2
              td(jce2,i,k) = atm1%t(jce2,i,k)*(d_one + ep1*qvd(jce2,i,k))
              ttld(jce2,i,k) = td(jce2,i,k) - sfs%psa(jce2,i) * &
                       t00pg*((hsigma(k)*sfs%psa(jce2,i)+ptop)/p00pg)**pgfaa1
            end do
          end do
        end if
        if ( ma%has_bdybottom ) then
          do k = 1 , kz
            do j = jce1 , jce2
              td(j,ice1,k) = atm1%t(j,ice1,k)*(d_one + ep1*qvd(j,ice1,k))
              ttld(j,ice1,k) = td(j,ice1,k) - sfs%psa(j,ice1) * &
                       t00pg*((hsigma(k)*sfs%psa(j,ice1)+ptop)/p00pg)**pgfaa1
            end do
          end do
        end if
        if ( ma%has_bdytop ) then
          do k = 1 , kz
            do j = jce1 , jce2
              td(j,ice2,k) = atm1%t(j,ice2,k)*(d_one + ep1*qvd(j,ice2,k))
              ttld(j,ice2,k) = td(j,ice2,k) - sfs%psa(j,ice2) * &
                       t00pg*((hsigma(k)*sfs%psa(j,ice2)+ptop)/p00pg)**pgfaa1
            end do
          end do
        end if
      else if ( ipgf == 0 ) then
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              tva = atm1%t(j,i,k)*(d_one+ep1*qvd(j,i,k))
              tvb = atm2%t(j,i,k)*(d_one+ep1*atm2%qx(j,i,k,iqv)*rpsb(j,i))
              tvc = atmc%t(j,i,k)*(d_one+ep1*atmc%qx(j,i,k,iqv)*rpsc(j,i))
              td(j,i,k) = alpha_hyd*(tvc+tvb) + beta_hyd*tva
            end do
          end do
        end do
        if ( ma%has_bdyleft ) then
          do k = 1 , kz
            do i = ici1 , ici2
              td(jce1,i,k) = atm1%t(jce1,i,k)*(d_one + ep1*qvd(jce1,i,k))
            end do
          end do
        end if
        if ( ma%has_bdyright ) then
          do k = 1 , kz
            do i = ici1 , ici2
              td(jce2,i,k) = atm1%t(jce2,i,k)*(d_one+ep1*qvd(jce2,i,k))
            end do
          end do
        end if
        if ( ma%has_bdybottom ) then
          do k = 1 , kz
            do j = jce1 , jce2
              td(j,ice1,k) = atm1%t(j,ice1,k)*(d_one+ep1*qvd(j,ice1,k))
            end do
          end do
        end if
        if ( ma%has_bdytop ) then
          do k = 1 , kz
            do j = jce1 , jce2
              td(j,ice2,k) = atm1%t(j,ice2,k)*(d_one+ep1*qvd(j,ice2,k))
            end do
          end do
        end if
      end if
    end if
    !
    ! compute the u and v tendencies:
    !
    ! compute the diffusion terms:
    ! put diffusion and pbl tendencies of u and v in difuu and difuv.
    !
    do k = 1 , kz
      do i = idi1 , idi2
        do j = jdi1 , jdi2
          adf%u(j,i,k) = aten%u(j,i,k)
          adf%v(j,i,k) = aten%v(j,i,k)
        end do
      end do
    end do
    !
    ! Zero again tendencies
    !
    aten%u(:,:,:) = d_zero
    aten%v(:,:,:) = d_zero
    !
    call diffu_d(adf%u,adf%v,atms%ubd3d,atms%vbd3d)
    !
    ! compute the horizontal advection terms for u and v:
    !
    ! compute the horizontal advection term in x and y momentum tendency:
    ! same for hydrostatic and nonhydrostatic models: 1st RHS term in
    ! Eqs. 2.1.1, 2.1.2, 2.2.1, 2.2.2, 2.2.9, 2.2.10, 2.3.3, 2.3.4
    !
    call hadv(aten%u,aten%v,atmx%ud,atmx%vd)
#ifdef DEBUG
    call check_wind_tendency('HADV')
#endif
    !
    ! compute Coriolis and curvature terms:
    !
    if ( idynamic == 1 ) then
      do k = 1 , kz
        do i = idi1 , idi2
          do j = jdi1 , jdi2
            !
            ! Hydrostatic model:
            ! (1) part of the horizontal component of the Coriolis force due
            ! to the horizontal movement (4th RHS term in Eq.2.1.1, Eq.2.1.2)
            !
            aten%u(j,i,k) = aten%u(j,i,k) + mddom%coriol(j,i)*atmx%vc(j,i,k)
            aten%v(j,i,k) = aten%v(j,i,k) - mddom%coriol(j,i)*atmx%uc(j,i,k)
          end do
        end do
      end do
    else if ( idynamic == 2 ) then
      do k = 1 , kz
        do i = idi1 , idi2
          do j = jdi1 , jdi2
            !
            ! Nonhydrostatic model:
            ! (1) part of the horizontal component of the Coriolis force
            !     due to the horizontal movement (5th RHS term in Eq.2.2.1,
            !     Eq.2.2.2, Eq.2.2.9, Eq.2.2.10, Eq.2.3.3, Eq.2.3.4)
            ! (2) part of the horizontal component of the Coriolis force
            !     due to the vertical movement (6th RHS term in Eq.2.2.9,
            !     Eq.2.2.10)
            ! (3) horizontal curvature term
            !     (not explicitly mentioned in the MM5 1994 manual)
            ! (4) vertical curvature term
            !     (not explicitly mentioned in the MM5 1994 manual)
            ! (5) divergence term
            !     (3rd RHS term in Eq.2.2.1, Eq.2.2.2, Eq.2.2.9,
            !      Eq.2.2.10, Eq.2.3.3, Eq.2.3.4)
            !
            wadot   = 0.125_rkx * (atm1%w(j-1,i-1,k) + atm1%w(j-1,i,k)     + &
                                 atm1%w(j,i-1,k)   + atm1%w(j,i,k))
            wadotp1 = 0.125_rkx * (atm1%w(j-1,i-1,k+1) + atm1%w(j-1,i,k+1) + &
                                 atm1%w(j,i-1,k+1)   + atm1%w(j,i,k+1))
            wabar = wadot + wadotp1
            amfac = wabar * rpsda(j,i) * rearthrad
            duv = atmx%uc(j,i,k)*mddom%dmdy(j,i)-atmx%vc(j,i,k)*mddom%dmdx(j,i)
            aten%u(j,i,k) = aten%u(j,i,k) +                      &
                            mddom%coriol(j,i)*atmx%vc(j,i,k) -   & ! H Coriolis
                            mddom%ef(j,i)*mddom%ddx(j,i)*wabar + & ! V Coriolis
                            atmx%vmd(j,i,k)*duv -                & ! H curv
                            atmx%uc(j,i,k)*amfac                   ! V curv
            aten%v(j,i,k) = aten%v(j,i,k) -                      &
                            mddom%coriol(j,i)*atmx%uc(j,i,k) +   & ! H Coriolis
                            mddom%ef(j,i)*mddom%ddy(j,i)*wabar - & ! V Coriolis
                            atmx%umd(j,i,k)*duv -                & ! H curv
                            atmx%vc(j,i,k)*amfac                   ! V curv
          end do
        end do
      end do
    end if
#ifdef DEBUG
    call check_wind_tendency('CORI')
#endif
    if ( idynamic == 1 ) then
      !
      ! compute pressure gradient terms:
      !
      if ( ipgf == 1 ) then
        do k = 1 , kz
          do i = idi1 , idi2
            do j = jdi1 , jdi2
              rtbar = d_rfour * (atmx%tv(j-1,i-1,k) + atmx%tv(j-1,i,k) + &
                                 atmx%tv(j,i-1,k)   + atmx%tv(j,i,k))  - &
                      t00pg*((hsigma(k)*sfs%psdota(j,i)+ptop)/p00pg)**pgfaa1
              rtbar = rgas*rtbar*sfs%psdota(j,i)
              !
              ! Hydrostatic model. The first part of the pressure gradient term:
              ! (1) 3rd RHS term in Eq.2.1.1, Eq.2.1.2., or
              ! (2) 2nd     term in Eq.2.4.1.
              ! (2a) Warning: there is missing sigma in the denominator in the
              !      MM5 manual (cf. Eq.2.4.1 in MM5 manual and Eq.4.2.6
              !      in MM4 manual)
              ! (2b) Hint:
              !      1/[1+p_top/(p* sigma)] dp*/dx = d log(sigma p* + p_top)/dx.
              !      This second form is discretized here.
              !
              aten%u(j,i,k) = aten%u(j,i,k) - rtbar *               &
                    (log(d_half*(sfs%psa(j,i)+sfs%psa(j,i-1))*      &
                          hsigma(k)+ptop) -                         &
                     log(d_half*(sfs%psa(j-1,i)+sfs%psa(j-1,i-1))*  &
                          hsigma(k)+ptop))/(dx*mddom%msfd(j,i))
              aten%v(j,i,k) = aten%v(j,i,k) - rtbar *               &
                    (log(d_half*(sfs%psa(j,i)+sfs%psa(j-1,i))*      &
                          hsigma(k)+ptop) -                         &
                     log(d_half*(sfs%psa(j-1,i-1)+sfs%psa(j,i-1))*  &
                          hsigma(k)+ptop))/(dx*mddom%msfd(j,i))
            end do
          end do
        end do
      else if ( ipgf == 0 ) then
        do k = 1 , kz
          do i = idi1 , idi2
            do j = jdi1 , jdi2
              rtbar = d_rfour * (atmx%tv(j-1,i-1,k) + atmx%tv(j-1,i,k) +   &
                                 atmx%tv(j,i-1,k)   + atmx%tv(j,i,k))
              rtbar = rgas*rtbar*sfs%psdota(j,i)
              !
              ! Hydrostatic model. The first part of the pressure gradient term:
              ! (1) in the 3rd RHS term in Eq.2.1.1, Eq.2.1.2., or
              ! (2)    the 2nd     term in Eq.2.4.1.
              ! (2a) Warning: there is missing sigma in the denominator in
              !      the MM5 manual (cf. Eq.2.4.1 in MM5 manual and Eq.4.2.6
              !      in MM4 manual)
              ! (2b) Hint:
              !      1/[1+p_top/(p* sigma)] dp*/dx = d log(sigma p* + p_top)/dx.
              !      This second form is discretized here.
              !
              aten%u(j,i,k) = aten%u(j,i,k) - rtbar *                &
                     (log(d_half*(sfs%psa(j,i)+sfs%psa(j,i-1))*      &
                           hsigma(k)+ptop) -                         &
                      log(d_half*(sfs%psa(j-1,i)+sfs%psa(j-1,i-1))*  &
                           hsigma(k)+ptop))/(dx*mddom%msfd(j,i))
              aten%v(j,i,k) = aten%v(j,i,k) - rtbar *                &
                     (log(d_half*(sfs%psa(j,i)+sfs%psa(j-1,i))*      &
                           hsigma(k)+ptop) -                         &
                      log(d_half*(sfs%psa(j-1,i-1)+sfs%psa(j,i-1))*  &
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
      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            !
            ! Hydrostatic model: the 2nd part of the Eq.2.4.5
            !
            tvfac(j,i,k) = d_one / (d_one+qcd(j,i,k)/(d_one+qvd(j,i,k)))
          end do
        end do
      end do
      if ( ipgf == 1 ) then
        do i = ice1 , ice2
          do j = jce1 , jce2
            !
            ! Hydrostatic model: the 1st part of the Eq.2.4.5
            !
            tv = ttld(j,i,kz)*rpsa(j,i)*tvfac(j,i,kz)
            phi(j,i,kz) = mddom%ht(j,i) + &
                     rgas*t00pg/pgfaa1*((sfs%psa(j,i)+ptop)/p00pg)**pgfaa1
            phi(j,i,kz) = phi(j,i,kz) - rgas * tv * &
                    log((hsigma(kz)+ptop*rpsa(j,i))/(d_one+ptop*rpsa(j,i)))
          end do
        end do
        do k = 1 , kzm1
          lev = kz - k
          do i = ice1 , ice2
            do j = jce1 , jce2
              !
              ! Hydrostatic model: the 1st part of the Eq.2.4.5
              !    (also, cf. Eq.2.1.9)
              !
              tvavg = ((ttld(j,i,lev)*dsigma(lev)+ttld(j,i,lev+1)*   &
                      dsigma(lev+1))/(sfs%psa(j,i)*(dsigma(lev)+     &
                      dsigma(lev+1))))*tvfac(j,i,lev)
              phi(j,i,lev) = phi(j,i,lev+1) - rgas *                 &
                     tvavg*log((hsigma(lev) + ptop*rpsa(j,i)) /      &
                               (hsigma(lev+1) + ptop*rpsa(j,i)))
            end do
          end do
        end do
      else if ( ipgf == 0 ) then
        do i = ice1 , ice2
          do j = jce1 , jce2
            !
            ! Hydrostatic model: the 1st part of the Eq.2.4.5
            !         (also, cf. Eq.2.1.9)
            !
            tv = td(j,i,kz)*rpsa(j,i)*tvfac(j,i,kz)
            phi(j,i,kz) = mddom%ht(j,i) - rgas * tv * &
                 log((hsigma(kz)+ptop*rpsa(j,i))/(d_one+ptop*rpsa(j,i)))
          end do
        end do
        do k = 1 , kzm1
          lev = kz - k
          do i = ice1 , ice2
            do j = jce1 , jce2
              !
              ! Hydrostatic model: the 1st part of the Eq.2.4.5
              !        (also, cf. Eq.2.1.9)
              !
              tvavg = ((td(j,i,lev)*dsigma(lev)+td(j,i,lev+1)*       &
                      dsigma(lev+1))/(sfs%psa(j,i)*(dsigma(lev)+     &
                      dsigma(lev+1))))*tvfac(j,i,lev)
              phi(j,i,lev) = phi(j,i,lev+1) - rgas *                 &
                     tvavg*log((hsigma(lev)+ptop*rpsa(j,i)) /        &
                                (hsigma(lev+1)+ptop*rpsa(j,i)))
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
            !
            ! Hydrostatic model. The second part of the pressure gradient term:
            ! (1) in the 3rd RHS term in Eq.2.1.1, Eq.2.1.2., or
            ! (2)    the 1st     term in Eq.2.4.1
            !
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
    end if ! idynamic == 1
    !
    ! compute the vertical advection terms:
    !
    ! compute the vertical advection term in x and y momentum tendency:
    ! same for hydrostatic and nonhydrostatic models: 2nd RHS term in
    ! Eqs. 2.1.1, 2.1.2, 2.2.1, 2.2.2, 2.2.9, 2.2.10, 2.3.3, 2.3.4
    !
    call vadv(aten%u,aten%v,atmx%uc,atmx%vc)
#ifdef DEBUG
    call check_wind_tendency('VADV')
#endif
    !
    ! apply the sponge boundary condition on u and v:
    !
    if ( iboudy == 4 ) then
      call sponge(xub,xvb,aten%u,aten%v)
    end if
    !
    ! apply the nudging boundary conditions:
    !
    if ( iboudy == 1 .or. iboudy == 5 ) then
      call nudge(iboudy,atm2%u,atm2%v,xub,xvb,aten%u,aten%v)
    end if
#ifdef DEBUG
    call check_wind_tendency('BDYC')
#endif
    !
    ! add the diffusion and pbl tendencies to aten%u and aten%v:
    ! Last RHS terms in Eq. 2.1.1, 2.1.2, 2.2.1, 2.2.2, 2.2.9, 2.2.10, 2.3.3,
    ! 2.3.4
    !
    do k = 1 , kz
      do i = idi1 , idi2
        do j = jdi1 , jdi2
          aten%u(j,i,k) = aten%u(j,i,k) + adf%u(j,i,k)
          aten%v(j,i,k) = aten%v(j,i,k) + adf%v(j,i,k)
        end do
      end do
    end do
#ifdef DEBUG
    call check_wind_tendency('DIFF')
#endif
    if ( ibltyp == 2 ) then
      !
      !  Couple TKE to ps for use in vertical advection
      !
      do k = 1 , kzp1
        do i = ice1 , ice2
          do j = jce1 , jce2
            uwstatea%tkeps(j,i,k) = atm1%tke(j,i,k)*sfs%psa(j,i)
          end do
        end do
      end do
      ! Calculate the horizontal advective tendency for TKE
      call hadv(uwstatea%advtke,atm1%tke,1)
      ! Calculate the vertical advective tendency for TKE
      call vadv(uwstatea%advtke,uwstatea%tkeps,kzp1,0)
      ! Calculate the horizontal, diffusive tendency for TKE
      ! The multiplication factor was causing instabilities
      ! in the non-hydrostatic model and has been removed.
      call diffu_x(uwstatea%advtke,atm2%tke,d_one)
    end if
    !
    ! Compute future values of t and moisture variables at tau+1:
    !
    if ( idynamic == 1 ) then
      call timefilter_apply(sfs%psa,sfs%psb,sfs%psc,gnu)
    end if
    call timefilter_apply(atm1%t,atm2%t,atmc%t,gnu)
    call timefilter_apply(atm1%qx,atm2%qx,atmc%qx,gnu,sfs%psb)
    if ( idynamic == 1 ) then
      call timefilter_apply(atm1%qx,atm2%qx,atmc%qx, &
                            d_two*gnu,iqfrst,iqlst,minqx)
    else
      call timefilter_apply(atm1%qx,atm2%qx,atmc%qx, &
                            gnu,iqfrst,iqlst,minqx)
    end if
    !
    if ( idynamic == 1 ) then
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
      ! perform time smoothing operations.
      !
      call timefilter_apply(atm1%u,atm2%u,atmc%u, &
                            atm1%v,atm2%v,atmc%v,gnu)
      do i = ice1 , ice2
        do j = jce1 , jce2
          rpsb(j,i) = d_one/sfs%psb(j,i)
        end do
      end do
      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            atm2%pr(j,i,k) = (hsigma(k)*sfs%psb(j,i) + ptop)*d_1000
          end do
        end do
      end do
      !
      ! Hydrostatic split explicit scheme
      !
      call splitf
    else if ( idynamic == 2 ) then
      !
      ! Decouple before calling sound
      !
      do k = 1 , kz
        do i = ide1 , ide2
          do j = jde1 , jde2
            aten%u(j,i,k) = aten%u(j,i,k) * rpsda(j,i)
            aten%v(j,i,k) = aten%v(j,i,k) * rpsda(j,i)
          end do
        end do
      end do
      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            aten%pp(j,i,k) = aten%pp(j,i,k) * rpsa(j,i)
          end do
        end do
      end do
      do k = 1 , kzp1
        do i = ice1 , ice2
          do j = jce1 , jce2
            aten%w(j,i,k) = aten%w(j,i,k) * rpsa(j,i)
          end do
        end do
      end do
      !
      ! Compute u,v,w,pp at ktau+1
      !
      call sound
      !
      ! Recompute new pressure
      !
      do k = 1 , kz
        do i = ice1 , ice2
          do j = jce1 , jce2
            atm2%pr(j,i,k) = atm0%pr(j,i,k) + atm2%pp(j,i,k)*rpsb(j,i)
          end do
        end do
      end do
    end if

    if ( ibltyp == 2 ) then
      ! TAO: Once the full loop above is completed, update the TKE
      ! tendency if the UW PBL is running.  NOTE!!! Do not try to
      ! combine these loops with the above loop Advection MUST be
      ! done in a loop separate from the updates.  (I lost 3 days
      ! of working to disocover that this is a problem because I
      ! thought it would be clever to combine loops--TAO)
      ! Add the advective tendency to the TKE tendency calculated
      ! by the UW TKE routine
      do k = 1 , kz
        do i = idi1 , idi2
          do j = jdi1 , jdi2
             aten%tke(j,i,k) = aten%tke(j,i,k) + &
                               uwstatea%advtke(j,i,k)*rpsa(j,i)
             atmc%tke(j,i,k) = max(tkemin,atm2%tke(j,i,k) + &
                               dt*aten%tke(j,i,k))
          end do
        end do
      end do
      call timefilter_apply(atm1%tke,atm2%tke,atmc%tke,gnu)
    end if ! TKE tendency update
    if ( ichem == 1 ) then
      call timefilter_apply(chia,chib,chic,d_two*gnu,1,ntr,mintr)
      !
      ! do cumulus simple transport/mixing of tracers for the schemes
      ! without explicit convective transport (Grell and KF up to now).
      ! works also in case 2 conv schemes over land and ocean are used,
      ! and if one of them is grell or kf. in this case trac tendency
      ! are not updated in the other convec scheme,
      ! cf mod_cu_em and mod_cu_tiedke.
      !
      if ( ktau > 0 .and. mod(ktau,ntcum) == 0 .and. &
           ichcumtra == 1 .and. any(icup == 2 .or. icup == 6) ) then
        call cumtran
      end if
    end if
    !
    ! Next timestep ready : increment elapsed forecast time
    !
    ktau = ktau + 1
    if ( islab_ocean == 1 ) xslabtime = xslabtime + dtsec
    idatex = idatex + intmdl
    if ( mod(ktau,khour) == 0 ) then
      call split_idate(idatex,xyear,xmonth,xday,xhour)
    end if
    if ( ktau == 2 ) then
      dtbat = dt*real(ntsrf,rkx)
      dt = dt2
      rdt = d_one/dt
      dtsq = dt*dt
      dtcb = dt*dt*dt
    end if
    !
    ! Print out noise parameter
    !
    if ( idynamic == 1 .and. ktau > 1 ) then
      if ( is_nan(ptntot) ) then
        maxv = abs(maxval(aten%t))
        if ( (maxv/dtsec) > 0.01_rkx ) then ! 50 K per hour
          write(stderr,*) 'MAXVAL ATEN T :', maxv
          maxv = maxv - 0.001_rkx
          do kk = 1 , kz
            do ii = ici1 , ici2
              do jj = jci1 , jci2
                if ( aten%t(jj,ii,kk) > maxv ) then
                  write(stderr,*) 'II :', ii, ', JJ :', jj, ', KK :', kk
                end if
              end do
            end do
          end do
        end if
        maxv = abs(maxval(aten%u))
        if ( (maxv/dtsec) > 0.005_rkx ) then  ! 25 m/s per hour
          write(stderr,*) 'MAXVAL ATEN U :', maxv
          maxv = maxv - 0.001_rkx
          do kk = 1 , kz
            do ii = ici1 , ici2
              do jj = jci1 , jci2
                if ( aten%u(jj,ii,kk) > maxv ) then
                  write(stderr,*) 'II :', ii, ', JJ :', jj, ', KK :', kk
                end if
              end do
            end do
          end do
        end if
        maxv = abs(maxval(aten%v))
        if ( (maxv/dtsec) > 0.005_rkx ) then  ! 25 m/s per hour
          write(stderr,*) 'MAXVAL ATEN V :', maxv
          maxv = maxv - 0.001_rkx
          do kk = 1 , kz
            do ii = ici1 , ici2
              do jj = jci1 , jci2
                if ( aten%v(jj,ii,kk) > maxv ) then
                  write(stderr,*) 'II :', ii, ', JJ :', jj, ', KK :', kk
                end if
              end do
            end do
          end do
        end if
        maxv = abs(maxval(aten%qx(:,:,:,iqv)))
        if ( (maxv/dtsec) > 0.001_rkx ) then !
          write(stderr,*) 'MAXVAL ATEN QV :', maxv
          maxv = maxv - 0.001_rkx
          do kk = 1 , kz
            do ii = ici1 , ici2
              do jj = jci1 , jci2
                if ( aten%qx(jj,ii,kk,iqv) > maxv ) then
                  write(stderr,*) 'II :', ii, ', JJ :', jj, ', KK :', kk
                end if
              end do
            end do
          end do
        end if
        maxv = abs(maxval(aten%qx(:,:,:,iqc)))
        if ( (maxv/dtsec) > 0.001_rkx ) then !
          write(stderr,*) 'MAXVAL ATEN QC :', maxv
          maxv = maxv - 0.001_rkx
          do kk = 1 , kz
            do ii = ici1 , ici2
              do jj = jci1 , jci2
                if ( aten%qx(jj,ii,kk,iqc) > maxv ) then
                  write(stderr,*) 'II :', ii, ', JJ :', jj, ', KK :', kk
                end if
              end do
            end do
          end do
        end if
        write (stderr,*) 'WHUUUUBBBASAAAGASDDWD!!!!!!!!!!!!!!!!'
        write (stderr,*) 'No more atmosphere here....'
        write (stderr,*) 'CFL violation detected, so model STOP'
        write (stderr,*) '#####################################'
        write (stderr,*) '#            DECREASE DT !!!!       #'
        write (stderr,*) '#####################################'
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
          if ( any(icup > 0) ) then
            write(stdout,'(a,i7)') &
              ' $$$ no. of points with active convection = ', iconvec
          end if
        end if
      end if
    end if

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif

    contains

#ifdef DEBUG
    ! Check temperature tendency less than 10 K

    subroutine check_temperature_tendency(loc)
      implicit none
      character(len=*) , intent(in) :: loc
      integer(ik4) :: i , j , k , kk , ierr
      real(rkx) :: check_tt , mean_tt
      ierr = 0
      mean_tt = (maxval(aten%t)+minval(aten%t))/d_two
      do k = 1 , kz
        do i = ici1, ici2
          do j = jci1 , jci2
            check_tt = (aten%t(j,i,k)-mean_tt)*rpsb(j,i)
            if ( abs(check_tt) > temp_tend_maxval ) then
              write(stderr,*) 'After ', loc, ' at ktau = ', ktau
              write(stderr,*) 'TEMP tendency out of order : ', check_tt
              write(stderr,*) 'At J = ',j
              write(stderr,*) 'At I = ',i
              write(stderr,*) 'At K = ',k
              write(stderr,*) 'Surface Temperature : ', sfs%tga(j,i)
              write(stderr,*) 'Vertical PTU profile: '
              do kk = 1 , kz
                write(stderr,'(i2,3f12.7)') kk, &
                        atms%pb3d(j,i,kk)*d_r100, atms%tb3d(j,i,kk), &
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
      real(rkx) :: check_ww , mean_ww
      ierr = 0
      wten = sqrt(max(aten%u(jde1:jde2,ide1:ide2,:),epsilon(d_one))**2 + &
                  max(aten%v(jde1:jde2,ide1:ide2,:),epsilon(d_one))**2)
      mean_ww = (maxval(wten)+minval(wten))/d_two
      do k = 1 , kz
        do i = ici1, ici2
          do j = jci1 , jci2
            check_ww = (wten(j,i,k)-mean_ww)/sfs%psdotb(j,i)
            if ( abs(check_ww) > wind_tend_maxval ) then
              write(stderr,*) 'After ', loc, ' at ktau = ', ktau
              write(stderr,*) 'WIND tendency out of order : ', check_ww
              write(stderr,*) 'At J = ',j
              write(stderr,*) 'At I = ',i
              write(stderr,*) 'At K = ',k
              write(stderr,*) 'Surface Temperature : ', sfs%tga(j,i)
              write(stderr,*) 'Vertical PTU profile: '
              do kk = 1 , kz
                write(stderr,'(i2,3f12.7)') kk, &
                        atms%pb3d(j,i,kk)*d_r100, atms%tb3d(j,i,kk), &
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

    subroutine surface_pressures( )
      implicit none
      integer(ik4) :: i , j
      logical , save :: linit = .false.
      !
      ! Compute surface pressure on dot points
      !
      if ( idynamic == 1 ) then
        call exchange(sfs%psa,1,jce1,jce2,ice1,ice2)
        call exchange(sfs%psb,2,jce1,jce2,ice1,ice2)
        do i = ice1ga , ice2ga
          do j = jce1ga , jce2ga
            rpsa(j,i) = d_one/sfs%psa(j,i)
          end do
        end do
        do i = ice1 , ice2
          do j = jce1 , jce2
            rpsb(j,i) = d_one/sfs%psb(j,i)
          end do
        end do
        call psc2psd(sfs%psa,sfs%psdota)
        call psc2psd(sfs%psb,sfs%psdotb)
        call exchange(sfs%psdota,1,jde1,jde2,ide1,ide2)
        call exchange(sfs%psdotb,2,jde1,jde2,ide1,ide2)
      else
        ! Non-hydrostatic pstar pressure is constant == ps0
        if ( .not. linit ) then
          do i = ice1ga , ice2ga
            do j = jce1ga , jce2ga
              rpsa(j,i) = d_one/sfs%psa(j,i)
            end do
          end do
          do i = ice1 , ice2
            do j = jce1 , jce2
              rpsb(j,i) = d_one/sfs%psb(j,i)
            end do
          end do
          do i = ide1ga , ide2ga
            do j = jde1ga , jde2ga
              rpsda(j,i) = d_one/sfs%psdota(j,i)
            end do
          end do
          linit = .true.
        end if
      end if
    end subroutine surface_pressures

    subroutine decouple( )
      implicit none
      integer(ik4) :: i , j , k
      !
      ! Helper
      !
      if ( idynamic == 1 ) then
        do i = ide1ga , ide2ga
          do j = jde1ga , jde2ga
            rpsda(j,i) = d_one/sfs%psdota(j,i)
          end do
        end do
      end if
      !
      ! Exchange ghost points
      !
      call exchange(atm1%u,1,jde1,jde2,ide1,ide2,1,kz)
      call exchange(atm1%v,1,jde1,jde2,ide1,ide2,1,kz)
      call exchange(atm1%t,1,jce1,jce2,ice1,ice2,1,kz)
      call exchange(atm1%qx,1,jce1,jce2,ice1,ice2,1,kz,1,nqx)
      if ( ibltyp == 2 ) then
        call exchange(atm1%tke,1,jce1,jce2,ice1,ice2,1,kzp1)
      end if
      !
      ! Coupled helper
      !
      do k = 1 , kz
        do i = ide1ga , ide2ga
          do j = jde1ga , jde2ga
            atmx%uc(j,i,k) = atm1%u(j,i,k)
            atmx%vc(j,i,k) = atm1%v(j,i,k)
            atmx%umc(j,i,k) = atm1%u(j,i,k)*mddom%msfd(j,i)
            atmx%vmc(j,i,k) = atm1%v(j,i,k)*mddom%msfd(j,i)
          end do
        end do
      end do
      !
      ! Decoupled part with boundary conditions
      !
      do k = 1 , kz
        do i = idi1 , idi2
          do j = jdi1 , jdi2
            atmx%ud(j,i,k) = atm1%u(j,i,k)*rpsda(j,i)
            atmx%vd(j,i,k) = atm1%v(j,i,k)*rpsda(j,i)
          end do
        end do
      end do
      !
      ! Boundary U,V points
      !
      if ( ma%has_bdyleft ) then
        do k = 1 , kz
          do i = idi1 , idi2
            atmx%ud(jdi1,i,k) = wui(i,k)*rpsda(jdi1,i)
            atmx%vd(jdi1,i,k) = wvi(i,k)*rpsda(jdi1,i)
          end do
        end do
        do k = 1 , kz
          do i = idi1 , idi2
            atmx%ud(jde1,i,k) = wue(i,k)*rpsda(jde1,i)
            atmx%vd(jde1,i,k) = wve(i,k)*rpsda(jde1,i)
          end do
        end do
        ! inflow/outflow dependence
        if ( iboudy == 3 .or. iboudy == 4 ) then
          do k = 1 , kz
            do i = idi1 , idi2
              if ( atm1%u(jde1,i,k) <= d_zero ) then
                atmx%ud(jde1,i,k) = atmx%ud(jdi1,i,k)
                atmx%vd(jde1,i,k) = atmx%vd(jdi1,i,k)
              end if
            end do
          end do
        end if
      end if
      if ( ma%has_bdyright ) then
        do k = 1 , kz
          do i = idi1 , idi2
            atmx%ud(jdi2,i,k) = eui(i,k)*rpsda(jdi2,i)
            atmx%vd(jdi2,i,k) = evi(i,k)*rpsda(jdi2,i)
          end do
        end do
        do k = 1 , kz
          do i = idi1 , idi2
            atmx%ud(jde2,i,k) = eue(i,k)*rpsda(jde2,i)
            atmx%vd(jde2,i,k) = eve(i,k)*rpsda(jde2,i)
          end do
        end do
        ! inflow/outflow dependence
        if ( iboudy == 3 .or. iboudy == 4 ) then
          do k = 1 , kz
            do i = idi1 , idi2
              if ( atm1%u(jde2,i,k) >= d_zero ) then
                atmx%ud(jde2,i,k) = atmx%ud(jdi2,i,k)
                atmx%vd(jde2,i,k) = atmx%vd(jdi2,i,k)
              end if
            end do
          end do
        end if
      end if
      if ( ma%has_bdybottom ) then
        do k = 1 , kz
          do j = jdi1 , jdi2
            atmx%ud(j,idi1,k) = sui(j,k)*rpsda(j,idi1)
            atmx%vd(j,idi1,k) = svi(j,k)*rpsda(j,idi1)
          end do
        end do
        do k = 1 , kz
          do j = jde1 , jde2
            atmx%ud(j,ide1,k) = sue(j,k)*rpsda(j,ide1)
            atmx%vd(j,ide1,k) = sve(j,k)*rpsda(j,ide1)
          end do
        end do
        if ( iboudy == 3 .or. iboudy == 4 ) then
          ! inflow/outflow dependence
          do k = 1 , kz
            do j = jde1 , jde2
              if ( atm1%v(j,ide1,k) >= d_zero ) then
                atmx%ud(j,ide1,k) = atmx%ud(j,idi1,k)
                atmx%vd(j,ide1,k) = atmx%vd(j,idi1,k)
              end if
            end do
          end do
        end if
      end if
      if ( ma%has_bdytop ) then
        do k = 1 , kz
          do j = jdi1 , jdi2
            atmx%ud(j,idi2,k) = nui(j,k)*rpsda(j,idi2)
            atmx%vd(j,idi2,k) = nvi(j,k)*rpsda(j,idi2)
          end do
        end do
        do k = 1 , kz
          do j = jde1 , jde2
            atmx%ud(j,ide2,k) = nue(j,k)*rpsda(j,ide2)
            atmx%vd(j,ide2,k) = nve(j,k)*rpsda(j,ide2)
          end do
        end do
        if ( iboudy == 3 .or. iboudy == 4 ) then
          ! inflow/outflow dependence
          do k = 1 , kz
            do j = jde1 , jde2
              if ( atm1%v(j,ide2,k) <= d_zero ) then
                atmx%ud(j,ide2,k) = atmx%ud(j,idi2,k)
                atmx%vd(j,ide2,k) = atmx%vd(j,idi2,k)
              end if
            end do
          end do
        end if
      end if
      if ( isladvec == 1 ) then
        call exchange(atmx%ud,2,jde1,jde2,ide1,ide2,1,kz)
        call exchange(atmx%vd,2,jde1,jde2,ide1,ide2,1,kz)
        do k = 1 , kz
          do i = ide1gb , ide2gb
            do j = jde1gb , jde2gb
              atmx%umd(j,i,k) = atmx%ud(j,i,k)*mddom%msfd(j,i)
              atmx%vmd(j,i,k) = atmx%vd(j,i,k)*mddom%msfd(j,i)
            end do
          end do
        end do
      else
        call exchange(atmx%ud,1,jde1,jde2,ide1,ide2,1,kz)
        call exchange(atmx%vd,1,jde1,jde2,ide1,ide2,1,kz)
        do k = 1 , kz
          do i = ide1ga , ide2ga
            do j = jde1ga , jde2ga
              atmx%umd(j,i,k) = atmx%ud(j,i,k)*mddom%msfd(j,i)
              atmx%vmd(j,i,k) = atmx%vd(j,i,k)*mddom%msfd(j,i)
            end do
          end do
        end do
      end if
      !
      ! T , QV , QC
      !
      do k = 1 , kz
        do i = ice1ga , ice2ga
          do j = jce1ga , jce2ga
            atmx%t(j,i,k) = atm1%t(j,i,k)*rpsa(j,i)
          end do
        end do
      end do
      do k = 1 , kz
        do i = ice1ga , ice2ga
          do j = jce1ga , jce2ga
            atmx%qx(j,i,k,iqv) = max(atm1%qx(j,i,k,iqv)*rpsa(j,i),minqq)
          end do
        end do
      end do
      do n = iqfrst , iqlst
        do k = 1 , kz
          do i = ice1ga , ice2ga
            do j = jce1ga , jce2ga
              atmx%qx(j,i,k,n) = max(atm1%qx(j,i,k,n)*rpsa(j,i),minqx)
            end do
          end do
        end do
      end do
      do k = 1 , kz
        do i = ice1ga , ice2ga
          do j = jce1ga , jce2ga
            atmx%tv(j,i,k) = atmx%t(j,i,k) * (d_one + ep1*atmx%qx(j,i,k,iqv))
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
                chi(j,i,k,n) = chia(j,i,k,n)*rpsa(j,i)
              end do
            end do
          end do
        end do
      end if

      if ( idynamic == 1 ) then
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              atm1%pr(j,i,k) = (hsigma(k)*sfs%psa(j,i) + ptop)*d_1000
              atm1%rho(j,i,k) = atm1%pr(j,i,k) / (rgas*atmx%tv(j,i,k))
            end do
          end do
        end do
      else
        !
        ! Constant reference state and perturbations are defined
        ! for the nonhydrostatic model.
        !
        call exchange(atm1%pp,1,jce1,jce2,ice1,ice2,1,kz)
        call exchange(atm1%w,1,jce1,jce2,ice1,ice2,1,kzp1)
        do k = 1 , kz
          do i = ice1ga , ice2ga
            do j = jce1ga , jce2ga
              atmx%pp(j,i,k) = atm1%pp(j,i,k)*rpsa(j,i)
            end do
          end do
        end do
        do k = 1 , kzp1
          do i = ice1ga , ice2ga
            do j = jce1ga , jce2ga
              atmx%w(j,i,k) = atm1%w(j,i,k)*rpsa(j,i)
            end do
          end do
        end do
        do k = 1 , kz
          do i = ice1ga , ice2ga
            do j = jce1ga , jce2ga
              atm1%pr(j,i,k) = atm0%pr(j,i,k) + atmx%pp(j,i,k)
              atm1%rho(j,i,k) = atm1%pr(j,i,k) / (rgas*atmx%tv(j,i,k))
            end do
          end do
        end do
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              atmx%pr(j,i,k) = (atmx%tv(j,i,k) - atm0%t(j,i,k) - &
                atmx%pp(j,i,k)/(cpd*atm0%rho(j,i,k))) / atmx%t(j,i,k)
            end do
          end do
        end do
      end if
      !
      ! Second timelevel exchange
      !
      call exchange(atm2%u,2,jde1,jde2,ide1,ide2,1,kz)
      call exchange(atm2%v,2,jde1,jde2,ide1,ide2,1,kz)
      call exchange(atm2%t,2,jce1,jce2,ice1,ice2,1,kz)
      if ( isladvec == 1 ) then
        call exchange(atm2%qx,4,jce1,jce2,ice1,ice2,1,kz,1,nqx)
      else
        call exchange(atm2%qx,2,jce1,jce2,ice1,ice2,1,kz,1,nqx)
      end if
      if ( ibltyp == 2 ) then
        call exchange(atm2%tke,2,jce1,jce2,ice1,ice2,1,kzp1)
      end if

      if ( idynamic == 1 ) then
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              atm2%pr(j,i,k) = (hsigma(k)*sfs%psb(j,i) + ptop)*d_1000
            end do
          end do
        end do
      else
        call exchange(atm2%pp,2,jce1,jce2,ice1,ice2,1,kz)
        call exchange(atm2%w,2,jce1,jce2,ice1,ice2,1,kzp1)
        !
        ! Constant reference state and perturbations are defined
        ! for the nonhydrostatic model.
        !
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              atm2%pr(j,i,k) = atm0%pr(j,i,k) + atm2%pp(j,i,k)*rpsb(j,i)
            end do
          end do
        end do
      end if

      if ( ichem == 1 ) then
        call exchange(chi,2,jce1,jce2,ice1,ice2,1,kz,1,ntr)
        if ( isladvec == 1 ) then
          call exchange(chib,4,jce1,jce2,ice1,ice2,1,kz,1,ntr)
        else
          call exchange(chib,2,jce1,jce2,ice1,ice2,1,kz,1,ntr)
        end if
      end if

    end subroutine decouple

  end subroutine tend

end module mod_tendency

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

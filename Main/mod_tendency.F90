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
  use mod_micro_interface
  use mod_bdycod
  use mod_slice
  use mod_sun
  use mod_advection
  use mod_diffusion
  use mod_domain
  use mod_sladvection
  use mod_slabocean
  use mod_sound
  use mod_split
  use mod_timefilter
  use mod_massck

  implicit none

  private

  public :: allocate_mod_tend , tend
  real(rkx) , pointer , dimension(:,:,:) :: ttld , td , &
         phi , qcd , qvd , tvfac , ucc , vcc , th , tha
  real(rkx) , pointer , dimension(:,:,:) :: ps4
  real(rkx) , pointer , dimension(:,:,:) :: ps_4
  real(rkx) , pointer , dimension(:,:) :: pten
  real(rkx) , pointer , dimension(:,:,:) :: thten
  real(rkx) , pointer , dimension(:,:) :: rpsa , rpsb , rpsc
  real(rkx) , pointer , dimension(:,:) :: rpsda

  real(rkx) , pointer , dimension(:,:,:) :: tten
  real(rkx) , pointer , dimension(:,:,:) :: uten
  real(rkx) , pointer , dimension(:,:,:) :: vten
  real(rkx) , pointer , dimension(:,:,:) :: wten
  real(rkx) , pointer , dimension(:,:,:) :: ppten
  real(rkx) , pointer , dimension(:,:,:) :: tketen
  real(rkx) , pointer , dimension(:,:,:,:) :: qxten
  real(rkx) , pointer , dimension(:,:,:,:) :: chiten

  real(rkx) , pointer , dimension(:,:,:) :: tdyn
  real(rkx) , pointer , dimension(:,:,:) :: udyn
  real(rkx) , pointer , dimension(:,:,:) :: vdyn
  real(rkx) , pointer , dimension(:,:,:) :: wdyn
  real(rkx) , pointer , dimension(:,:,:) :: ppdyn
  real(rkx) , pointer , dimension(:,:,:) :: tkedyn
  real(rkx) , pointer , dimension(:,:,:,:) :: qxdyn
  real(rkx) , pointer , dimension(:,:,:,:) :: chidyn

  real(rkx) , pointer , dimension(:,:,:) :: tphy
  real(rkx) , pointer , dimension(:,:,:) :: uphy
  real(rkx) , pointer , dimension(:,:,:) :: vphy
  real(rkx) , pointer , dimension(:,:,:) :: wphy
  real(rkx) , pointer , dimension(:,:,:) :: ppphy
  real(rkx) , pointer , dimension(:,:,:) :: tkephy
  real(rkx) , pointer , dimension(:,:,:,:) :: qxphy
  real(rkx) , pointer , dimension(:,:,:,:) :: chiphy

  real(rkx) , pointer , dimension(:,:,:) :: ten0 , qen0
  real(rkx) , pointer , dimension(:,:,:,:) :: chiten0
  real(rkx) , pointer , dimension(:,:,:) :: tkeps

  integer :: idgq
  integer :: ithadv = 1
  integer(ik4) :: iqxvadv , itrvadv

  real(rkx) :: rptn ! Total number of internal points

  interface ten2diag
    module procedure extracttent
    module procedure extracttenqv
    module procedure extracttenchi
  end interface ten2diag

  contains

#include <cpmf.inc>

  subroutine allocate_mod_tend
    implicit none
    call getmem3d(ps_4,jcross1,jcross2,icross1,icross2,1,4,'tendency:ps_4')
    call getmem3d(ps4,jci1,jci2,ici1,ici2,1,4,'tendency:ps4')
    if ( ipptls > 1 ) then
      call getmem3d(qcd,jce1,jce2,ice1,ice2,1,kz,'tendency:qcd')
    else
      call assignpnt(atmx%qx,qcd,iqc)
    end if
    call getmem3d(tvfac,jce1,jce2,ice1,ice2,1,kz,'tendency:tvfac')
    call assignpnt(atmx%qx,qvd,iqv)
    call getmem2d(pten,jce1,jce2,ice1,ice2,'tendency:pten')
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
        call getmem3d(thten,jci1,jci2,ici1,ici2,1,kz,'tendency:thten')
        call getmem3d(th,jce1ga,jce2ga,ice1ga,ice2ga,1,kz,'tendency:th')
        call getmem3d(tha,jce1,jce2,ice1,ice2,1,kz,'tendency:tha')
      end if
    end if
    !
    ! Set number of ghost points for advection for the two schemes
    ! Select advection scheme
    !
    iqxvadv = 1
    itrvadv = 2
    if ( ibltyp == 2 ) then
      if ( iuwvadv == 1 ) then
        iqxvadv = 3
        itrvadv = 3
      end if
      call getmem3d(tkeps,jce1,jce2,ice1,ice2,1,kzp1,'tendency:tkeps')
    end if

    if ( idiag == 1 ) then
      idgq = iqv
      call getmem3d(ten0,jci1,jci2,ici1,ici2,1,kz,'tendency:ten0')
      call getmem3d(qen0,jci1,jci2,ici1,ici2,1,kz,'tendency:qen0')
    end if

    if ( ichem == 1 ) then
      if ( ichdiag == 1 ) then
        call getmem4d(chiten0,jci1,jci2,ici1,ici2,1,kz,1,ntr,'tendency:chiten0')
      end if
    end if
    call assignpnt(aten%t,tten,pc_total)
    call assignpnt(aten%t,tdyn,pc_dynamic)
    call assignpnt(aten%t,tphy,pc_physic)

    call assignpnt(aten%u,uten,pc_total)
    call assignpnt(aten%u,udyn,pc_dynamic)
    call assignpnt(aten%u,uphy,pc_physic)

    call assignpnt(aten%v,vten,pc_total)
    call assignpnt(aten%v,vdyn,pc_dynamic)
    call assignpnt(aten%v,vphy,pc_physic)

    call assignpnt(aten%qx,qxten,pc_total)
    call assignpnt(aten%qx,qxdyn,pc_dynamic)
    call assignpnt(aten%qx,qxphy,pc_physic)

    if ( idynamic == 2 ) then
      call assignpnt(aten%w,wten,pc_total)
      call assignpnt(aten%w,wdyn,pc_dynamic)
      call assignpnt(aten%w,wphy,pc_physic)
      call assignpnt(aten%pp,ppten,pc_total)
      call assignpnt(aten%pp,ppdyn,pc_dynamic)
      call assignpnt(aten%pp,ppphy,pc_physic)
    end if
    if ( ibltyp == 2 ) then
      call assignpnt(aten%tke,tketen,pc_total)
      call assignpnt(aten%tke,tkedyn,pc_dynamic)
      call assignpnt(aten%tke,tkephy,pc_physic)
    end if
    if ( ichem == 1 ) then
      call assignpnt(aten%chi,chiten,pc_total)
      call assignpnt(aten%chi,chidyn,pc_dynamic)
      call assignpnt(aten%chi,chiphy,pc_physic)
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
    real(rkx) :: pt2bar , pt2tot , ptnbar , maxv , ptntot , &
                 rovcpm , cpm , rofac , uaq , vaq , scr1
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
    ! Compute pressure tendency, mass-divergence, sigma-velocity and
    ! pressure vertical velocity
    !
    call compute_omega
    !
    ! Prepare fields to be used in physical parametrizations.
    !
    call mkslice
    !
    ! Compute forecast surface pressure for hydrostatic
    !
    if ( idynamic == 1 ) then
      call new_pressure
    end if
    !
    ! calculate solar zenith angle
    !
    call zenitm(coszrs)
    !
    ! Compute new diffusion coefficients
    !
    call calc_coeff
    !
    ! Initialize the tendencies
    !
    call init_tendencies
    !
    ! Boundary conditions term
    !
    call boundary
    !
    ! Compute Horizontal advection terms, curvature terms and adiabatic term
    !
    call advection
    call curvature
    call adiabatic
    !
    ! Compute horizontal diffusion term
    !
    call diffusion
    !
    ! Physical parametrizations
    !
    call physical_parametrizations
    !
    ! Compute chemistry tendencies (other than transport)
    !
    if ( ichem == 1 ) then
      call tractend2(ktau,xyear,xmonth,xday,calday,declin)
    end if
    !
    ! Sum up all tendencies for temperature
    !
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          tten(j,i,k) = tten(j,i,k) + tdyn(j,i,k) + tphy(j,i,k)
        end do
      end do
    end do
    !
    ! Sum up all contribution to water vapor tendencies
    ! This is the last RHS term in Eqs. 2.1.3 and 2.2.5, 2.3.9
    !
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          qxten(j,i,k,iqv) = qxten(j,i,k,iqv) + qxdyn(j,i,k,iqv) + &
                             qxphy(j,i,k,iqv)
        end do
      end do
    end do
    if ( idynamic == 2 ) then
      !
      ! Sum up all tendencies for vertical wind and pressure perturbation
      ! This is the last RHS term in Eqs. 2.2.3, 2.2.11, 2.3.7
      !
      do k = 1 , kzp1
        do i = ici1 , ici2
          do j = jci1 , jci2
            wten(j,i,k) = wten(j,i,k) + wdyn(j,i,k) + wphy(j,i,k)

          end do
        end do
      end do
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            ppten(j,i,k) = ppten(j,i,k) + ppdyn(j,i,k) + ppphy(j,i,k)
          end do
        end do
      end do
    end if
    if ( ipptls > 0 ) then
      do n = iqfrst , iqlst
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              qxten(j,i,k,n) = qxten(j,i,k,n) + qxdyn(j,i,k,n) + &
                               qxphy(j,i,k,n)
            end do
          end do
        end do
      end do
      if ( ipptls == 1 ) then
        !
        ! compute the condensation and precipitation terms for SUBEX
        ! moisture schemes
        !
        qxphy(:,:,:,iqv) = d_zero
        qxphy(:,:,:,iqc) = d_zero
        tphy(:,:,:) = d_zero
        call condtq
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              qxten(j,i,k,iqv) = qxten(j,i,k,iqv) + qxphy(j,i,k,iqv)
              qxten(j,i,k,iqc) = qxten(j,i,k,iqc) + qxphy(j,i,k,iqc)
              tten(j,i,k) = tten(j,i,k) + tphy(j,i,k)
            end do
          end do
        end do
      end if
      if ( idiag > 0 ) then
        call ten2diag(aten%t,tdiag%lsc,pc_physic)
        call ten2diag(aten%qx,qdiag%lsc,pc_physic)
      end if
    end if
    if ( idynamic == 2 ) then
      if ( itopnudge == 1 ) then
        call topnudge(atm2%t,tten,xtb)
        call topnudge(atm2%qx,qxten,xqb,iqv)
      end if
      if ( ifrayd == 1 ) then
        call raydamp(atms%za,atm2%t,tten)
        call raydamp(atms%za,atm2%qx,qxten)
      end if
    end if
    !
    ! forecast t, qv, and qc at tau+1:
    !
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          atmc%t(j,i,k) = atm2%t(j,i,k) + dt * tten(j,i,k)
        end do
      end do
    end do
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          atmc%qx(j,i,k,iqv) = max(atm2%qx(j,i,k,iqv) + &
                       dt * qxten(j,i,k,iqv),minqv)
        end do
      end do
    end do
    do n = iqfrst , iqlst
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            atmc%qx(j,i,k,n) = max(atm2%qx(j,i,k,n) + &
                       dt * qxten(j,i,k,n),minqx)
          end do
        end do
      end do
    end do
    !
    ! Pressure Gradient Force for Hydrostatic core
    !
    if ( idynamic == 1 ) then
      call pressure_gradient_force
    end if
    !
    ! Last RHS terms in Eq. 2.1.1, 2.1.2, 2.2.1, 2.2.2, 2.2.9, 2.2.10, 2.3.3,
    ! 2.3.4
    !
    do k = 1 , kz
      do i = idi1 , idi2
        do j = jdi1 , jdi2
          uten(j,i,k) = uten(j,i,k) + udyn(j,i,k) + uphy(j,i,k)
          vten(j,i,k) = vten(j,i,k) + vdyn(j,i,k) + vphy(j,i,k)
        end do
      end do
    end do
    !
    ! Check mass
    !
    if ( debug_level > 0 ) call massck
    !
    ! Compute future values of t and moisture variables at tau+1:
    !
    if ( idynamic == 1 ) then
      call timefilter_apply(sfs%psa,sfs%psb,sfs%psc,gnu1)
    end if
    call timefilter_apply(atm1%t,atm2%t,atmc%t,gnu1)
    call timefilter_apply(atm1%qx,atm2%qx,atmc%qx,gnu1,sfs%psb)
    call timefilter_apply(atm1%qx,atm2%qx,atmc%qx,gnu2,iqfrst,iqlst,minqx)
    !
    if ( idynamic == 1 ) then
      !
      ! forecast p*u and p*v at tau+1:
      !
      do k = 1 , kz
        do i = idi1 , idi2
          do j = jdi1 , jdi2
            atmc%u(j,i,k) = atm2%u(j,i,k) + dt * uten(j,i,k)
            atmc%v(j,i,k) = atm2%v(j,i,k) + dt * vten(j,i,k)
          end do
        end do
      end do
      !
      ! perform time smoothing operations.
      !
      call timefilter_apply(atm1%u,atm2%u,atmc%u, &
                            atm1%v,atm2%v,atmc%v,gnu1)
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
      if ( itopnudge == 1 ) then
        call topnudge(atm2%u,atm2%v,uten,vten)
        call topnudge(atm2%pp,ppten,xppb)
        call topnudge(atm2%w,wten,xwwb)
      end if
      if ( ifrayd == 1 ) then
        call raydamp(atms%za,atm2%u,atm2%v,uten,vten)
        call raydamp(atms%za,atm2%pp,ppten)
        call raydamp(atms%zq,atm2%w,wten,xwwb)
      end if
      do k = 1 , kz
        do i = idi1 , idi2
          do j = jdi1 , jdi2
            uten(j,i,k) = uten(j,i,k) * rpsda(j,i)
            vten(j,i,k) = vten(j,i,k) * rpsda(j,i)
          end do
        end do
      end do
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            ppten(j,i,k) = ppten(j,i,k) * rpsa(j,i)
          end do
        end do
      end do
      do k = 1 , kzp1
        do i = ici1 , ici2
          do j = jci1 , jci2
            wten(j,i,k) = wten(j,i,k) * rpsa(j,i)
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
    !
    ! forecast TKE at at tau+1:
    !
    if ( ibltyp == 2 ) then
      ! TAO: Once the full loop above is completed, update the TKE
      ! tendency if the UW PBL is running.  NOTE!!! Do not try to
      ! combine these loops with the above loop Advection MUST be
      ! done in a loop separate from the updates.  (I lost 3 days
      ! of working to disocover that this is a problem because I
      ! thought it would be clever to combine loops--TAO)
      ! Add the advective tendency to the TKE tendency calculated
      ! by the UW TKE routine
      do k = 1 , kzp1
        do i = ici1 , ici2
          do j = jci1 , jci2
            tketen(j,i,k) = tketen(j,i,k) +  &
                         tkedyn(j,i,k) * rpsa(j,i) + tkephy(j,i,k)
          end do
        end do
      end do
      if ( idynamic == 2 ) then
        if ( ifrayd == 1 ) then
          call raydamp(atms%zq,atm2%tke,tketen,tkemin)
        end if
      end if
      do k = 1 , kzp1
        do i = ici1 , ici2
          do j = jci1 , jci2
             atmc%tke(j,i,k) = max(tkemin,atm2%tke(j,i,k) + &
                               dt * tketen(j,i,k))
          end do
        end do
      end do
      call timefilter_apply(atm1%tke,atm2%tke,atmc%tke,gnu2)
    end if ! TKE tendency update
    !
    ! forecast tracer chi at at tau+1:
    !
    if ( ichem == 1 ) then
      do itr = 1 , ntr
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              chiten(j,i,k,itr) = chiten(j,i,k,itr) + &
                                  chidyn(j,i,k,itr) + chiphy(j,i,k,itr)
            end do
          end do
        end do
      end do
      if ( idynamic == 2 ) then
        if ( ifrayd == 1 ) then
          call raydamp(atms%za,atm2%chi,chiten)
        end if
      end if
      do itr = 1 , ntr
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              atmc%chi(j,i,k,itr) = max(atm2%chi(j,i,k,itr) + &
                       dt * chiten(j,i,k,itr),mintr)
            end do
          end do
        end do
      end do
      call timefilter_apply(atm1%chi,atm2%chi,atmc%chi,gnu2,1,ntr,mintr)
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
        maxv = abs(maxval(tten(:,:,:)))
        if ( (maxv/dtsec) > 0.01_rkx ) then ! 50 K per hour
          write(stderr,*) 'MAXVAL ATEN T :', maxv
          maxv = maxv - 0.001_rkx
          do kk = 1 , kz
            do ii = ici1 , ici2
              do jj = jci1 , jci2
                if ( tten(jj,ii,kk) > maxv ) then
                  write(stderr,*) 'II :', ii, ', JJ :', jj, ', KK :', kk
                end if
              end do
            end do
          end do
        end if
        maxv = abs(maxval(uten(:,:,:)))
        if ( (maxv/dtsec) > 0.005_rkx ) then  ! 25 m/s per hour
          write(stderr,*) 'MAXVAL ATEN U :', maxv
          maxv = maxv - 0.001_rkx
          do kk = 1 , kz
            do ii = ici1 , ici2
              do jj = jci1 , jci2
                if ( uten(jj,ii,kk) > maxv ) then
                  write(stderr,*) 'II :', ii, ', JJ :', jj, ', KK :', kk
                end if
              end do
            end do
          end do
        end if
        maxv = abs(maxval(vten(:,:,:)))
        if ( (maxv/dtsec) > 0.005_rkx ) then  ! 25 m/s per hour
          write(stderr,*) 'MAXVAL ATEN V :', maxv
          maxv = maxv - 0.001_rkx
          do kk = 1 , kz
            do ii = ici1 , ici2
              do jj = jci1 , jci2
                if ( vten(jj,ii,kk) > maxv ) then
                  write(stderr,*) 'II :', ii, ', JJ :', jj, ', KK :', kk
                end if
              end do
            end do
          end do
        end if
        maxv = abs(maxval(qxten(:,:,:,iqv)))
        if ( (maxv/dtsec) > 0.001_rkx ) then !
          write(stderr,*) 'MAXVAL ATEN QV :', maxv
          maxv = maxv - 0.001_rkx
          do kk = 1 , kz
            do ii = ici1 , ici2
              do jj = jci1 , jci2
                if ( qxten(jj,ii,kk,iqv) > maxv ) then
                  write(stderr,*) 'II :', ii, ', JJ :', jj, ', KK :', kk
                end if
              end do
            end do
          end do
        end if
        maxv = abs(maxval(qxten(:,:,:,iqc)))
        if ( (maxv/dtsec) > 0.001_rkx ) then !
          write(stderr,*) 'MAXVAL ATEN QC :', maxv
          maxv = maxv - 0.001_rkx
          do kk = 1 , kz
            do ii = ici1 , ici2
              do jj = jci1 , jci2
                if ( qxten(jj,ii,kk,iqc) > maxv ) then
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

    subroutine check_temperature_tendency(loc,ipc)
      implicit none
      character(len=*) , intent(in) :: loc
      integer(ik4) , intent(in) :: ipc
      integer(ik4) :: i , j , k , kk , ierr
      real(rkx) :: check_tt , mean_tt
      ierr = 0
      mean_tt = (maxval(aten%t(:,:,:,ipc))+minval(aten%t(:,:,:,ipc)))/d_two
      do k = 1 , kz
        do i = ici1, ici2
          do j = jci1 , jci2
            check_tt = (aten%t(j,i,k,ipc)-mean_tt)*rpsb(j,i)
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

    subroutine check_wind_tendency(loc,ipc)
      implicit none
      character(len=*) , intent(in) :: loc
      integer(ik4) , intent(in) :: ipc
      integer(ik4) :: i , j , k , kk , ierr
      real(rkx) :: check_ww , mean_ww
      real(rkx) , dimension(jdi1:jdi2,idi1:idi2,1:kz) :: ww
      ierr = 0
      ww = sqrt(max(aten%u(jdi1:jdi2,idi1:idi2,:,ipc),epsilon(d_one))**2 + &
                max(aten%v(jdi1:jdi2,idi1:idi2,:,ipc),epsilon(d_one))**2)
      mean_ww = (maxval(ww)+minval(ww))/d_two
      do k = 1 , kz
        do i = ici1, ici2
          do j = jci1 , jci2
            check_ww = (ww(j,i,k)-mean_ww)/sfs%psdotb(j,i)
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
        call exchange(sfs%psb,1,jce1,jce2,ice1,ice2)
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
        call exchange(sfs%psdotb,1,jde1,jde2,ide1,ide2)
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
      if ( ichem == 1 ) then
        call exchange(atm1%chi,1,jce1,jce2,ice1,ice2,1,kz,1,ntr)
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
            do i = ice1ga , ice2ga
              do j = jce1ga , jce2ga
                atmx%chi(j,i,k,n) = atm1%chi(j,i,k,n)*rpsa(j,i)
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
      call exchange(atm2%u,1,jde1,jde2,ide1,ide2,1,kz)
      call exchange(atm2%v,1,jde1,jde2,ide1,ide2,1,kz)
      call exchange(atm2%t,1,jce1,jce2,ice1,ice2,1,kz)
      if ( isladvec == 1 ) then
        call exchange(atm2%qx,4,jce1,jce2,ice1,ice2,1,kz,1,nqx)
      else
        call exchange(atm2%qx,1,jce1,jce2,ice1,ice2,1,kz,1,nqx)
      end if
      if ( ibltyp == 2 ) then
        call exchange(atm2%tke,1,jce1,jce2,ice1,ice2,1,kzp1)
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
        call exchange(atm2%pp,1,jce1,jce2,ice1,ice2,1,kz)
        call exchange(atm2%w,1,jce1,jce2,ice1,ice2,1,kzp1)
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
        if ( isladvec == 1 ) then
          call exchange(atm2%chi,4,jce1,jce2,ice1,ice2,1,kz,1,ntr)
        else
          call exchange(atm2%chi,1,jce1,jce2,ice1,ice2,1,kz,1,ntr)
        end if
      end if
      !
      ! Total water load
      !
      if ( ipptls > 1 ) then
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
    end subroutine decouple

    subroutine compute_omega
      implicit none
      integer(ik4) :: i , j , k
      real(rkx) , dimension(jde1:jde2,ide1:ide2) :: dummy

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
              ucc(j,i,k) = (atmx%ud(j,i,k)  + atmx%ud(j,i+1,k) + &
                            atmx%ud(j+1,i,k)+ atmx%ud(j+1,i+1,k))
              vcc(j,i,k) = (atmx%vd(j,i,k)  + atmx%vd(j,i+1,k) + &
                            atmx%vd(j+1,i,k)+ atmx%vd(j+1,i+1,k))
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
              qdot(j,i,k) = -atm0%rhof(j,i,k)*egrav* &
                      atmx%w(j,i,k)/atm0%ps(j,i) - &
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
    end subroutine compute_omega

    subroutine init_tendencies
      implicit none
      aten%u(:,:,:,:) = d_zero
      aten%v(:,:,:,:) = d_zero
      aten%t(:,:,:,:) = d_zero
      aten%qx(:,:,:,:,:) = d_zero
      !
      ! Pressure perturbations and vertical velocity tendencies in the
      ! nonhydrostatic model
      !
      if ( idynamic == 2 ) then
        aten%pp(:,:,:,:) = d_zero
        aten%w(:,:,:,:) = d_zero
      end if
      !
      ! TKE for UW pbl
      !
      if ( ibltyp == 2 ) then
        aten%tke(:,:,:,:) = d_zero
      end if
      !
      ! Chemistry
      !
      if ( ichem == 1 ) then
        aten%chi(:,:,:,:,:)  = d_zero
      end if
      !
      ! Cloud fractions
      !
      cldfra(:,:,:) = d_zero
      cldlwc(:,:,:) = d_zero
    end subroutine init_tendencies

    subroutine advection
      implicit none
      integer(ik4) :: i , j , k

      call start_advect

      if ( idiag > 0 ) then
        ten0 = tdyn
        qen0 = qxdyn(:,:,:,idgq)
      end if
      !
      ! Compute departure points for semi-lagrangian
      !
      if ( isladvec == 1 ) then
        call trajcalc_x
      end if
      !
      ! compute the horizontal advection terms for u and v:
      !
      ! compute the horizontal advection term in x and y momentum tendency:
      ! same for hydrostatic and nonhydrostatic models: 1st RHS term in
      ! Eqs. 2.1.1, 2.1.2, 2.2.1, 2.2.2, 2.2.9, 2.2.10, 2.3.3, 2.3.4
      !
      call hadv(udyn,vdyn,atmx%ud,atmx%vd)
#ifdef DEBUG
      call check_wind_tendency('HADV',pc_dynamic)
#endif
      !
      ! compute the vertical advection terms:
      !
      ! compute the vertical advection term in x and y momentum tendency:
      ! same for hydrostatic and nonhydrostatic models: 2nd RHS term in
      ! Eqs. 2.1.1, 2.1.2, 2.2.1, 2.2.2, 2.2.9, 2.2.10, 2.3.3, 2.3.4
      !
      call vadv(udyn,vdyn,atmx%uc,atmx%vc)
#ifdef DEBUG
      call check_wind_tendency('VADV',pc_dynamic)
#endif
      if ( idynamic == 2 ) then
        !
        ! Horizontal and vertical advection of pressure perturbation and
        ! vertical velocity in the nonhydrostatic model: 1st and 2nd term
        ! on the RHS of the Eq. 2.2.3, Eq. 2.2.4, Eq. 2.3.7 and Eq. 2.3.8
        ! in the MM5 manual.
        !
        ! Also, cf. Eq. 2.2.11 of vertical velocity tendency in the MM5 manual.
        !
        call hadv(ppdyn,atmx%pp,0)
        call vadv(ppdyn,atm1%pp,kz,0)
        call hadv(wdyn,atmx%w,1)
        call vadv(wdyn,atm1%w,kzp1,0)
      end if
      !
      ! compute the horizontal advection term in temperature tendency:
      ! same for hydrostatic and nonhydrostatic models:
      ! in Eqs. 2.1.3, 2.2.5, 2.3.9 (1st RHS term)
      !
      if ( ithadv == 0 ) then
        call hadv(tdyn,atmx%t)
#ifdef DEBUG
        call check_temperature_tendency('HADV',pc_dynamic)
#endif
        if ( idiag > 0 ) then
          call ten2diag(aten%t,tdiag%adh,pc_dynamic,ten0)
          ten0 = tdyn
        end if
        !
        ! compute the vertical advection term:
        ! same for hydrostatic and nonhydrostatic models:
        ! in Eqs. 2.1.3, 2.2.5, 2.3.9 (2nd RHS term)
        !
        call vadv(tdyn,atm1%t,kz,1)
#ifdef DEBUG
        call check_temperature_tendency('VADV',pc_dynamic)
#endif
        if ( idiag > 0 ) then
          call ten2diag(aten%t,tdiag%adv,pc_dynamic,ten0)
        end if
      else
        thten(:,:,:) = d_zero
        do k = 1 , kz
          do i = ice1 , ice2
            do j = jce1 , jce2
              th(j,i,k) = atmx%t(j,i,k) * (p00/atm1%pr(j,i,k))**rovcp
              tha(j,i,k) = th(j,i,k) * sfs%psa(j,i)
            end do
          end do
        end do
        call exchange(th,1,jce1,jce2,ice1,ice2,1,kz)
        call hadv(thten,th)
        call vadv(thten,tha,kz,0)
      end if
      !
      ! compute the moisture tendencies
      !
      if ( isladvec == 1 ) then
        call slhadv_x(qxdyn,atm2%qx,iqv)
        call hdvg_x(qxdyn,atm1%qx,iqv)
      else
        call hadv(qxdyn,atmx%qx,iqv)
      end if
      if ( idiag > 0 .and. idgq == iqv ) then
        call ten2diag(aten%qx,qdiag%adh,pc_dynamic,qen0)
        qen0 = qxdyn(:,:,:,iqv)
      end if
      if ( all(icup /= 1) ) then
        call vadv(qxdyn,atm1%qx,iqv)
      end if
      if ( idiag > 0 .and. idgq == iqv ) then
        call ten2diag(aten%qx,qdiag%adv,pc_dynamic,qen0)
      end if
      if ( ipptls > 0 ) then
        if ( isladvec == 1 ) then
          call slhadv_x(qxdyn,atm2%qx,iqfrst,iqlst)
          call hdvg_x(qxdyn,atm1%qx,iqfrst,iqlst)
        else
          call hadv(qxdyn,atmx%qx,iqfrst,iqlst)
        end if
        if ( idiag > 0 .and. idgq /= iqv ) then
          call ten2diag(aten%qx,qdiag%adh,pc_dynamic)
          qen0 = qxdyn(:,:,:,idgq)
        end if
        call vadv(qxdyn,atm1%qx,iqfrst,iqlst,iqxvadv)
        if ( idiag > 0 .and. idgq /= iqv ) then
          call ten2diag(aten%qx,qdiag%adv,pc_dynamic,qen0)
        end if
      end if
      if ( ichem == 1 ) then
        !
        ! horizontal and vertical advection + diag
        !
        if ( isladvec == 1 ) then
          call slhadv_x(chidyn,atm2%chi)
          call hdvg_x(chidyn,atm1%chi)
        else
          call hadv(chidyn,atmx%chi)
        end if
        if ( ichdiag == 1 ) then
          call ten2diag(aten%chi,cadvhdiag,pc_dynamic)
          chiten0 = chidyn
        end if
        if ( all(icup /= 1) ) then
          call vadv(chidyn,atm1%chi,1,ntr,itrvadv)
        end if
        if ( ichdiag == 1 ) then
          call ten2diag(aten%chi,cadvvdiag,pc_dynamic,chiten0)
        end if
      end if
      if ( ibltyp == 2 ) then
        ! Calculate the horizontal advective tendency for TKE
        call hadv(tkedyn,atm1%tke,1)
        !
        !  Couple TKE to ps for use in vertical advection
        !
        do k = 1 , kzp1
          do i = ice1 , ice2
            do j = jce1 , jce2
              tkeps(j,i,k) = atm1%tke(j,i,k)*sfs%psa(j,i)
            end do
          end do
        end do
        ! Calculate the vertical advective tendency for TKE
        call vadv(tkedyn,tkeps,kzp1,0)
      end if
    end subroutine advection

    subroutine new_pressure
      implicit none
      integer(ik4) :: i , j
      !
      ! Surface pressure boundary conditions
      !
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
      ! Compute bleck (1977) noise parameters:
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
    end subroutine new_pressure

    subroutine boundary
      implicit none
      if ( iboudy == 1 .or. iboudy == 5 ) then
        call nudge(iboudy,atm2%t,xtb,tten)
        call nudge(iboudy,atm2%qx,xqb,qxten,iqv)
        call nudge(iboudy,atm2%u,atm2%v,xub,xvb,uten,vten)
      else if ( iboudy == 4 ) then
        call sponge(xtb,tten)
        call sponge(xqb,qxten,iqv)
        call sponge(xub,xvb,uten,vten)
      end if
#ifdef DEBUG
      call check_temperature_tendency('BDYC',pc_total)
      call check_wind_tendency('BDYC',pc_total)
#endif
      if ( idynamic == 2 ) then
        if ( iboudy == 1 .or. iboudy == 5 ) then
          call nudge(iboudy,atm2%pp,xppb,ppten)
          call nudge(iboudy,atm2%w,xwwb,wten)
        else if ( iboudy == 4 ) then
          call sponge(xppb,ppten)
          call sponge(xwwb,wten)
        end if
      end if
      if ( ichem == 1 ) then
        if ( iboudy == 1 .or. iboudy == 5 ) then
          call nudge_chi(kz,atm2%chi,chiten)
        end if
      end if
      if ( idiag > 0 ) then
        call ten2diag(aten%t,tdiag%bdy,pc_total)
        call ten2diag(aten%qx,qdiag%bdy,pc_total)
        if ( ichem == 1 ) then
          call ten2diag(aten%chi,cbdydiag,pc_total)
        end if
      end if
    end subroutine boundary

    subroutine diffusion
      implicit none
      !
      ! compute the diffusion term for t and qx
      !
      call diffu_d(uphy,vphy,atms%ubd3d,atms%vbd3d)
      if ( idiag > 0 ) then
        ten0 = tphy
        qen0 = qxphy(:,:,:,iqv)
      end if
      call diffu_x(tphy,atms%tb3d)
      call diffu_x(qxphy,atms%qxb3d,1,nqx,d_one)
      if ( idiag > 0 ) then
        call ten2diag(aten%t,tdiag%dif,pc_physic,ten0)
        call ten2diag(aten%qx,qdiag%dif,pc_physic,qen0)
      end if
      if ( idynamic == 2 ) then
        !
        ! compute the diffusion term for vertical velocity w
        ! compute the diffusion term for perturb pressure pp
        !
        call diffu_x(ppphy,atms%ppb3d)
        call diffu_x(wphy,atms%wb3d,d_one)
      end if
      if ( ichem == 1 ) then
        if ( ichdiag > 0 ) chiten0 = chiphy
        call diffu_x(chiphy,atms%chib3d,1,ntr,d_one)
        if ( ichdiag > 0 ) call ten2diag(aten%chi,cdifhdiag,pc_physic,chiten0)
      end if
      if ( ibltyp == 2 ) then
        ! Calculate the horizontal, diffusive tendency for TKE
        ! Here TKE is decoupled , we can pass atm2. Use dynamic component
        ! of tendency, to allow later decoupling to sum up all.
        call diffu_x(tkedyn,atm2%tke,nuk)
      end if
#ifdef DEBUG
      call check_wind_tendency('DIFF',pc_physic)
      call check_temperature_tendency('DIFF',pc_physic)
#endif
    end subroutine diffusion

    subroutine adiabatic
      implicit none
      if ( idiag > 0 ) then
        ten0 = tdyn
        qen0 = qxdyn(:,:,:,iqv)
      end if
      if ( idynamic == 1 ) then
        !
        ! Adiabatic term in the temperature tendency equation in the
        ! hydrostatic model:    3rd RHS term in Eq. 2.1.3
        !
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              rovcpm = rgas/cpmf(qvd(j,i,k))
              tdyn(j,i,k) = tdyn(j,i,k) +  &
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
                scr1 = d_half*egrav*atm0%rho(j,i,k) * &
                             (atm1%w(j,i,k)+atm1%w(j,i,k+1))
                tdyn(j,i,k) = tdyn(j,i,k) + atmx%t(j,i,k)*mdv%cr(j,i,k) - &
                           (scr1 + ppdyn(j,i,k) + ppten(j,i,k) +   &
                           atmx%pp(j,i,k)*mdv%cr(j,i,k))/(atm1%rho(j,i,k)*cpm)
              end do
            end do
          end do
        else
          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
                thten(j,i,k) = thten(j,i,k) + th(j,i,k) * mdv%cr(j,i,k)
                tdyn(j,i,k) = tdyn(j,i,k) + &
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
              ppdyn(j,i,k) = ppdyn(j,i,k) + atmx%pp(j,i,k)*mdv%cr(j,i,k)
            end do
          end do
        end do
        do n = 1 , nqx
          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
                qxdyn(j,i,k,n) = qxdyn(j,i,k,n) + atmx%qx(j,i,k,n)*mdv%cr(j,i,k)
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
              wdyn(j,i,k) = wdyn(j,i,k) + &
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
                wdyn(j,i,k) = wdyn(j,i,k) - egrav * sfs%psa(j,i) * &
                            (twt(k,2)*qcd(j,i,k-1) + twt(k,1)*qcd(j,i,k))
              end do
            end do
          end do
        end if
      end if
      if ( idiag > 0 ) then
        call ten2diag(aten%t,tdiag%adi,pc_dynamic,ten0)
        call ten2diag(aten%qx,qdiag%adi,pc_dynamic,qen0)
      end if
#ifdef DEBUG
      call check_temperature_tendency('DIAB',pc_dynamic)
#endif
    end subroutine adiabatic

    subroutine physical_parametrizations
      implicit none
      !
      ! First guess heating from previous radiation call
      !
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            tphy(j,i,k) = tphy(j,i,k) + sfs%psb(j,i)*heatrt(j,i,k)
          end do
        end do
      end do
      !
      !------------------------------------------------
      !        Call cumulus parametrization
      !------------------------------------------------
      !
      if ( idiag > 0 ) then
        ten0 = tphy
        qen0 = qxphy(:,:,:,idgq)
      end if
      if ( ichem == 1 .and. ichdiag == 1 ) then
        chiten0 = chiphy
      end if
      call cumulus
#ifdef DEBUG
      call check_temperature_tendency('CONV',pc_physic)
      call check_wind_tendency('CONV',pc_physic)
#endif
      if ( idiag > 0 ) then
        call ten2diag(aten%t,tdiag%con,pc_physic,ten0)
        call ten2diag(aten%qx,qdiag%con,pc_physic,qen0)
      end if
      if ( ichem == 1 .and. ichdiag == 1 ) then
        call ten2diag(aten%chi,cconvdiag,pc_physic,chiten0)
      end if
      ! save cumulus cloud fraction for chemistry before it is
      ! overwritten in cldfrac
      if ( ichem == 1 ) convcldfra(:,:,:) = cldfra(:,:,:)
      !
      !------------------------------------------------
      ! Large scale precipitation microphysical schemes
      !------------------------------------------------
      !
      if ( ipptls > 0 ) then
        if ( idiag > 0 ) then
          ten0 = tphy
          qen0 = qxphy(:,:,:,idgq)
        end if
        ! Clouds and large scale precipitation
        call cldfrac
        call microscheme
        if ( idiag > 0 ) then
          call ten2diag(aten%t,tdiag%lsc,pc_physic,ten0)
          call ten2diag(aten%qx,qdiag%lsc,pc_physic,qen0)
        end if
#ifdef DEBUG
        call check_temperature_tendency('MICR',pc_physic)
#endif
      end if
      !
      !------------------------------------------------
      !            Call Surface model
      !------------------------------------------------
      !
      if ( ktau == 0 .or. mod(ktau+1,ntsrf) == 0 ) then
        call surface_model
        if ( islab_ocean == 1 ) call update_slabocean(xslabtime)
      end if
      !
      !------------------------------------------------
      !             Call PBL scheme
      !------------------------------------------------
      !
      if ( idiag > 0 ) then
        ten0 = tphy
        qen0 = qxphy(:,:,:,idgq)
      end if
      if ( ichem == 1 .and. ichdiag == 1 ) then
        chiten0 = chiphy
      end if
      call pblscheme
      if ( idiag > 0 ) then
        call ten2diag(aten%t,tdiag%tbl,pc_physic,ten0)
        call ten2diag(aten%qx,qdiag%tbl,pc_physic,qen0)
      end if
      if ( ichem == 1 .and. ichdiag == 1 ) then
        call ten2diag(aten%chi,ctbldiag,pc_physic,chiten0)
      end if
#ifdef DEBUG
      call check_temperature_tendency('PBLL',pc_physic)
      call check_wind_tendency('PBLL',pc_physic)
#endif
      !
      !------------------------------------------------
      !       Call radiative transfer package
      !------------------------------------------------
      !
      if ( ktau == 0 .or. mod(ktau+1,ntrad) == 0 ) then
        !
        ! Remove first guess heating from previous radiation
        !
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              tphy(j,i,k) = tphy(j,i,k) - sfs%psb(j,i)*heatrt(j,i,k)
            end do
          end do
        end do
        ! calculate albedo
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
        !
        ! Add radiative transfer package-calculated heating rates to
        ! temperature tendency (deg/sec)
        !
        if ( idiag > 0 ) ten0 = tphy
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              tphy(j,i,k) = tphy(j,i,k) + sfs%psb(j,i)*heatrt(j,i,k)
            end do
          end do
        end do
        if ( idiag > 0 ) call ten2diag(aten%t,tdiag%rad,pc_physic,ten0)
#ifdef DEBUG
        call check_temperature_tendency('HEAT',pc_physic)
#endif
      end if
    end subroutine physical_parametrizations

    subroutine curvature
      implicit none
      integer(ik4) :: i , j , k
      real(rkx) :: wadot , wadotp1 , wabar , amfac , duv
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
              udyn(j,i,k) = udyn(j,i,k) + mddom%coriol(j,i)*atmx%vc(j,i,k)
              vdyn(j,i,k) = vdyn(j,i,k) - mddom%coriol(j,i)*atmx%uc(j,i,k)
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
              duv = atmx%uc(j,i,k)*mddom%dmdy(j,i) - &
                    atmx%vc(j,i,k)*mddom%dmdx(j,i)
              udyn(j,i,k) = udyn(j,i,k) + &
                            mddom%coriol(j,i)*atmx%vc(j,i,k) -   & ! H Coriolis
                            mddom%ef(j,i)*mddom%ddx(j,i)*wabar + & ! V Coriolis
                            atmx%vmd(j,i,k)*duv -                & ! H curv
                            atmx%uc(j,i,k)*amfac                   ! V curv
              vdyn(j,i,k) = vdyn(j,i,k) - &
                            mddom%coriol(j,i)*atmx%uc(j,i,k) +   & ! H Coriolis
                            mddom%ef(j,i)*mddom%ddy(j,i)*wabar - & ! V Coriolis
                            atmx%umd(j,i,k)*duv -                & ! H curv
                            atmx%vc(j,i,k)*amfac                   ! V curv
            end do
          end do
        end do
      end if
#ifdef DEBUG
      call check_wind_tendency('CORI',pc_dynamic)
#endif
    end subroutine curvature

    subroutine pressure_gradient_force
      implicit none
      integer(ik4) :: i , j , k
      real(rkx) :: tva , tvb , tvc , rtbar , tv , tvavg
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
              udyn(j,i,k) = udyn(j,i,k) - rtbar * &
                    (log(d_half*(sfs%psa(j,i)+sfs%psa(j,i-1))*      &
                          hsigma(k)+ptop) -                         &
                     log(d_half*(sfs%psa(j-1,i)+sfs%psa(j-1,i-1))*  &
                          hsigma(k)+ptop))/(dx*mddom%msfd(j,i))
              vdyn(j,i,k) = vdyn(j,i,k) - rtbar * &
                    (log(d_half*(sfs%psa(j,i)+sfs%psa(j-1,i))*      &
                          hsigma(k)+ptop) -                         &
                     log(d_half*(sfs%psa(j-1,i-1)+sfs%psa(j,i-1))*  &
                          hsigma(k)+ptop))/(dx*mddom%msfd(j,i))
            end do
          end do
        end do
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
              udyn(j,i,k) = udyn(j,i,k) - rtbar * &
                     (log(d_half*(sfs%psa(j,i)+sfs%psa(j,i-1))*      &
                           hsigma(k)+ptop) -                         &
                      log(d_half*(sfs%psa(j-1,i)+sfs%psa(j-1,i-1))*  &
                           hsigma(k)+ptop))/(dx*mddom%msfd(j,i))
              vdyn(j,i,k) = vdyn(j,i,k) - rtbar * &
                     (log(d_half*(sfs%psa(j,i)+sfs%psa(j-1,i))*      &
                           hsigma(k)+ptop) -                         &
                      log(d_half*(sfs%psa(j-1,i-1)+sfs%psa(j,i-1))*  &
                           hsigma(k)+ptop))/(dx*mddom%msfd(j,i))
            end do
          end do
        end do
      end if
#ifdef DEBUG
      call check_wind_tendency('PRGR',pc_dynamic)
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
            udyn(j,i,k) = udyn(j,i,k) - &
                 (sfs%psa(j-1,i-1)+sfs%psa(j-1,i)+  &
                  sfs%psa(j,i-1)+sfs%psa(j,i)) *    &
                 (phi(j,i,k)+phi(j,i-1,k)-          &
                  phi(j-1,i,k)-phi(j-1,i-1,k)) / (dx8*mddom%msfd(j,i))
            vdyn(j,i,k) = vdyn(j,i,k) - &
                 (sfs%psa(j-1,i-1)+sfs%psa(j-1,i)+  &
                  sfs%psa(j,i-1)+sfs%psa(j,i)) *    &
                 (phi(j,i,k)+phi(j-1,i,k)-          &
                  phi(j,i-1,k)-phi(j-1,i-1,k)) / (dx8*mddom%msfd(j,i))
          end do
        end do
      end do
#ifdef DEBUG
      call check_wind_tendency('GEOP',pc_dynamic)
#endif
    end subroutine pressure_gradient_force

  end subroutine tend

  subroutine extracttent(ten,store,indx,ten0)
    implicit none
    real(rkx) , pointer , dimension(:,:,:,:) :: ten
    real(rkx) , pointer , dimension(:,:,:) :: store
    integer(ik4) , intent(in) :: indx
    real(rkx) , optional , pointer , dimension(:,:,:) :: ten0
    integer :: i , j , k
    if ( present(ten0) ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            store(j,i,k) = store(j,i,k) + (ten(j,i,k,indx)-ten0(j,i,k))*afdout
          end do
        end do
      end do
      return
    end if
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          store(j,i,k) = store(j,i,k) + ten(j,i,k,indx)*afdout
        end do
      end do
    end do
  end subroutine extracttent

  subroutine extracttenqv(qen,store,indx,qen0)
    implicit none
    real(rkx) , pointer , dimension(:,:,:,:,:) :: qen
    real(rkx) , pointer , dimension(:,:,:) :: store
    integer(ik4) , intent(in) :: indx
    real(rkx) , optional , pointer , dimension(:,:,:) :: qen0
    integer :: i , j , k
    if ( present(qen0) ) then
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            store(j,i,k) = store(j,i,k) + &
                    (qen(j,i,k,idgq,indx)-qen0(j,i,k))*afdout
          end do
        end do
      end do
      return
    end if
    do k = 1 , kz
      do i = ici1 , ici2
        do j = jci1 , jci2
          store(j,i,k) = store(j,i,k) + qen(j,i,k,idgq,indx)*afdout
        end do
      end do
    end do
  end subroutine extracttenqv

  subroutine extracttenchi(chiten,store,indx,chiten0)
    implicit none
    real(rkx) , pointer , dimension(:,:,:,:,:) :: chiten
    real(rkx) , pointer , dimension(:,:,:,:) :: store
    integer(ik4) , intent(in) :: indx
    real(rkx) , optional , pointer , dimension(:,:,:,:) :: chiten0
    integer :: i , j , k , n
    if ( present(chiten0) ) then
      do n = 1 , ntr
        do k = 1 , kz
          do i = ici1 , ici2
            do j = jci1 , jci2
              store(j,i,k,n) = store(j,i,k,n) + &
                  (chiten(j,i,k,n,indx)-chiten0(j,i,k,n))*afdout
            end do
          end do
        end do
      end do
      return
    end if
    do n = 1 , ntr
      do k = 1 , kz
        do i = ici1 , ici2
          do j = jci1 , jci2
            store(j,i,k,n) = store(j,i,k,n) + chiten(j,i,k,n,indx)*afdout
          end do
        end do
      end do
    end do
  end subroutine extracttenchi

end module mod_tendency

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

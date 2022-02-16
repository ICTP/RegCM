!:::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
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

module mod_output

  use mod_intkinds
  use mod_realkinds
  use mod_dynparam
  use mod_runparams
  use mod_header
  use mod_mpmessage
  use mod_mppparam
  use mod_service
  use mod_projections
  use mod_atm_interface
  use mod_che_interface
  use mod_che_output
  use mod_lm_interface
  use mod_rad_interface
  use mod_cu_interface
  use mod_pbl_interface
  use mod_micro_interface
  use mod_ncout
  use mod_bdycod
  use mod_split
  use mod_savefile
  use mod_slabocean
  use mod_moloch
  use mod_capecin

  implicit none

  private

  public :: output

  type(rcm_time_and_date) , save , public :: lastout

  type(regcm_projection) , save :: pj

  logical :: rotinit = .false.

  interface uvrot
    module procedure uvrot2d
    module procedure uvrot3d
  end interface uvrot

  contains

  subroutine output
    implicit none
    logical :: ldoatm , ldosrf , ldorad , ldoche, ldoopt
    logical :: ldosav , ldolak , ldosub , ldosts , ldoshf , lnewf
    logical :: ldoslab
    logical :: lstartup
    integer(ik4) :: i , j , k , kk , itr
    real(rkx) , dimension(kz) :: p1d , t1d , rh1d
    real(rkx) :: cell , zz , zz1 , ww , tv
    real(rkx) :: srffac , srafac , radfac , lakfac , subfac , optfac , stsfac
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'output'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif

    lstartup = .false.
    if ( rcmtimer%start( ) .or. doing_restart ) then
      !
      ! Set up static variables (first time in)
      !
      if ( associated(xlon_out) ) then
        xlon_out = mddom%xlon(jci1:jci2,ici1:ici2)
        xlat_out = mddom%xlat(jci1:jci2,ici1:ici2)
        mask_out = mddom%mask(jci1:jci2,ici1:ici2)
        topo_out = mddom%ht(jci1:jci2,ici1:ici2)
        topo_out = topo_out*regrav
      end if
      if ( associated(sub_xlon_out) ) then
        call reorder_subgrid(mdsub%xlon,sub_xlon_out)
        call reorder_subgrid(mdsub%xlat,sub_xlat_out)
        call reorder_subgrid(mdsub%mask,sub_mask_out)
        call reorder_subgrid(mdsub%ht,sub_topo_out)
        sub_topo_out = sub_topo_out*regrav
      end if
      if ( idynamic == 2 ) then
        if ( associated(p0_out) ) then
          p0_out = atm0%ps(jci1:jci2,ici1:ici2) + ptop*d_1000
        end if
      end if
      !
      ! Reset the accumulation arrays
      !
      if ( associated(sts_tgmax_out) )  sts_tgmax_out  = -1.e30_rkx
      if ( associated(sts_tgmin_out) )  sts_tgmin_out  =  1.e30_rkx
      if ( associated(sts_t2max_out) )  sts_t2max_out  = -1.e30_rkx
      if ( associated(sts_t2min_out) )  sts_t2min_out  =  1.e30_rkx
      if ( associated(sts_w10max_out) ) sts_w10max_out = -1.e30_rkx
      if ( associated(sts_psmin_out) )  sts_psmin_out  =  1.e30_rkx
      if ( associated(sts_pcpmax_out) ) sts_pcpmax_out = -1.e30_rkx
      call newoutfiles(rcmtimer%idate)
      lastout = rcmtimer%idate
      lstartup = .true.
      if ( doing_restart ) then
        doing_restart = .false.
#ifdef DEBUG
        call time_end(subroutine_name,idindx)
#endif
        return
      end if
    end if

    lnewf = .false.
    ldoatm = .false.
    ldosrf = .false.
    ldolak = .false.
    ldosub = .false.
    ldorad = .false.
    ldoopt = .false.
    ldoche = .false.
    ldosav = .false.
    ldoslab = .false.
    ldosts = .false.
    ldoshf = .false.

    if ( rcmtimer%integrating( ) ) then
      if ( associated(alarm_out_nwf) ) then
        if ( alarm_out_nwf%act( ) .and. .not. rcmtimer%reached_endtime ) then
          lnewf = .true.
        end if
      end if
      if ( associated(alarm_out_sav) ) then
        if ( savfrq > d_zero ) then
          if ( rcmtimer%reached_endtime .or. alarm_out_sav%act( ) ) then
            ldosav = .true.
          end if
        else
          if ( rcmtimer%reached_endtime .or. &
               alarm_out_sav%act( ) .or. &
               (lfdomonth(rcmtimer%idate) .and. &
                lmidnight(rcmtimer%idate)) ) then
            ldosav = .true.
          end if
        end if
      else
        if ( ( rcmtimer%reached_endtime ) .or. &
             (lfdomonth(rcmtimer%idate) .and. &
              lmidnight(rcmtimer%idate)) ) then
          ldosav = .true.
        end if
      end if
      if ( alarm_out_atm%act( ) ) then
        ldoatm = .true.
      end if
      if ( alarm_out_srf%act( ) ) then
        ldosrf = .true.
      end if
      if ( alarm_out_sts%act( ) ) then
        ldosts = .true.
      end if
      if ( alarm_out_shf%act( ) ) then
        ldoshf = .true.
      end if
      if ( lakemod == 1 ) then
        if ( alarm_out_lak%act( ) ) then
          ldolak= .true.
        end if
      end if
      if ( nsg > 1 ) then
        if ( alarm_out_sub%act( ) ) then
          ldosub= .true.
        end if
      end if
      if ( alarm_out_rad%act( ) ) then
        ldorad = .true.
      end if
      if ( ichem == 1 ) then
        if ( alarm_out_che%act( ) ) then
          ldoche = .true.
        end if
      end if
      if ( alarm_out_opt%act( ) ) then
        ldoopt = .true.
      end if
      if ( rcmtimer%reached_endtime ) then
        ldoslab = .true.
      end if
    end if

    if ( rcmtimer%start( ) ) then
      ldoatm = .true.
      if ( ichem == 1 ) then
        ldoche = .true.
      end if
    end if

    if ( atm_stream > 0 ) then
      if ( ldoatm ) then
        ps_out = sfs%psa(jci1:jci2,ici1:ici2)
        if ( associated(atm_t_out) ) then
          if ( idynamic == 3 ) then
            do k = 1 , kz
              atm_t_out(:,:,k) = mo_atm%t(jci1:jci2,ici1:ici2,k)
            end do
          else
            do k = 1 , kz
              atm_t_out(:,:,k) = atm1%t(jci1:jci2,ici1:ici2,k)/ps_out
            end do
          end if
        end if
        if ( associated(atm_u_out) .and. associated(atm_v_out) ) then
          if ( idynamic == 3 ) then
            call uvstagtox(mo_atm%u,mo_atm%v,mo_atm%ux,mo_atm%vx)
            do k = 1 , kz
              do i = ici1 , ici2
                do j = jci1 , jci2
                  atm_u_out(j,i,k) = mo_atm%ux(j,i,k)
                  atm_v_out(j,i,k) = mo_atm%vx(j,i,k)
                end do
              end do
            end do
          else
            call exchange(atm1%u,1,jde1,jde2,ide1,ide2,1,kz)
            call exchange(atm1%v,1,jde1,jde2,ide1,ide2,1,kz)
            do k = 1 , kz
              do i = ici1 , ici2
                do j = jci1 , jci2
                  atm_u_out(j,i,k) = d_rfour*(atm1%u(j,i,k)+atm1%u(j+1,i,k) + &
                                atm1%u(j,i+1,k)+atm1%u(j+1,i+1,k))/ps_out(j,i)
                  atm_v_out(j,i,k) = d_rfour*(atm1%v(j,i,k)+atm1%v(j+1,i,k) + &
                                atm1%v(j,i+1,k)+atm1%v(j+1,i+1,k))/ps_out(j,i)
                end do
              end do
            end do
          end if
          if ( uvrotate ) then
            call uvrot(atm_u_out,atm_v_out)
          end if
        end if
        if ( associated(atm_omega_out) ) then
          atm_omega_out = omega(jci1:jci2,ici1:ici2,:)*d_10
        end if
        if ( associated(atm_w_out) ) then
          if ( idynamic == 3 ) then
            call wstagtox(mo_atm%w,atm_w_out)
          else
            do k = 1 , kz
              do i = ici1 , ici2
                do j = jci1 , jci2
                  atm_w_out(j,i,k) = d_half * (atm1%w(j,i,k+1) + &
                                               atm1%w(j,i,k)) / ps_out(j,i)
                end do
              end do
            end do
          end if
        end if
        if ( associated(atm_pai_out) ) then
          do k = 1 , kz
            atm_pai_out(:,:,k) = mo_atm%pai(jci1:jci2,ici1:ici2,k)
          end do
        end if
        if ( associated(atm_pp_out) ) then
          do k = 1 , kz
            atm_pp_out(:,:,k) = atm1%pp(jci1:jci2,ici1:ici2,k)/ps_out
          end do
        end if
        if ( associated(atm_qv_out) ) then
          if ( idynamic == 3 ) then
            do k = 1 , kz
              atm_qv_out(:,:,k) = mo_atm%qx(jci1:jci2,ici1:ici2,k,iqv)
            end do
          else
            do k = 1 , kz
              atm_qv_out(:,:,k) = atm1%qx(jci1:jci2,ici1:ici2,k,iqv)/ps_out
            end do
          end if
          ! Specific humidity in the output, not mixing ratio
          atm_qv_out = atm_qv_out/(d_one+atm_qv_out)
        end if
        if ( associated(atm_qc_out) ) then
          if ( idynamic == 3 ) then
            do k = 1 , kz
              atm_qc_out(:,:,k) = mo_atm%qx(jci1:jci2,ici1:ici2,k,iqc)
            end do
          else
            do k = 1 , kz
              atm_qc_out(:,:,k) = atm1%qx(jci1:jci2,ici1:ici2,k,iqc)/ps_out
            end do
          end if
        end if
        if ( associated(atm_qr_out) ) then
          if ( idynamic == 3 ) then
            do k = 1 , kz
              atm_qr_out(:,:,k) = mo_atm%qx(jci1:jci2,ici1:ici2,k,iqr)
            end do
          else
            do k = 1 , kz
              atm_qr_out(:,:,k) = atm1%qx(jci1:jci2,ici1:ici2,k,iqr)/ps_out
            end do
          end if
        end if
        if ( associated(atm_qi_out) ) then
          if ( idynamic == 3 ) then
            do k = 1 , kz
              atm_qi_out(:,:,k) = mo_atm%qx(jci1:jci2,ici1:ici2,k,iqi)
            end do
          else
            do k = 1 , kz
              atm_qi_out(:,:,k) = atm1%qx(jci1:jci2,ici1:ici2,k,iqi)/ps_out
            end do
          end if
        end if
        if ( associated(atm_qs_out) ) then
          if ( idynamic == 3 ) then
            do k = 1 , kz
              atm_qs_out(:,:,k) = mo_atm%qx(jci1:jci2,ici1:ici2,k,iqs)
            end do
          else
            do k = 1 , kz
              atm_qs_out(:,:,k) = atm1%qx(jci1:jci2,ici1:ici2,k,iqs)/ps_out
            end do
          end if
        end if
        if ( associated(atm_rh_out) ) then
          if ( idynamic == 3 ) then
            do k = 1 , kz
              do i = ici1 , ici2
                do j = jci1 , jci2
                  atm_rh_out(j,i,k) = d_100 * &
                    min(rhmax,max(rhmin,(mo_atm%qx(j,i,k,iqv) / &
                           pfwsat(mo_atm%t(j,i,k),mo_atm%p(j,i,k)))))
                end do
              end do
            end do
          else
            do k = 1 , kz
              do i = ici1 , ici2
                do j = jci1 , jci2
                  atm_rh_out(j,i,k) = d_100 * min(rhmax,max(rhmin, &
                     (atm1%qx(j,i,k,iqv)/ps_out(j,i)) / &
                     pfwsat(atm1%t(j,i,k)/ps_out(j,i),atm1%pr(j,i,k))))
                end do
              end do
            end do
          end if
        end if
        if ( associated(atm_pf_out) ) then
          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
                atm_pf_out(j,i,k) = atms%pf3d(j,i,k)
              end do
            end do
          end do
        end if
        if ( idynamic == 1 ) then
          if ( associated(atm_ph_out) ) then
            do k = 1 , kz
              do i = ici1 , ici2
                do j = jci1 , jci2
                  atm_ph_out(j,i,k) = (sigma(k)*sfs%psa(j,i)+ptop)*d_1000
                end do
              end do
            end do
          end if
          if ( associated(atm_zh_out) ) then
            do i = ici1 , ici2
              do j = jci1 , jci2
                cell = ptop / sfs%psa(j,i)
                atm_zh_out(j,i,kz) = rovg * atm1%t(j,i,kz)/sfs%psa(j,i) * &
                     log((sigma(kzp1)+cell)/(sigma(kz)+cell))
                do k = kz-1 , 1 , -1
                  atm_zh_out(j,i,k) = atms%zq(j,i,k+1) +&
                      rovg * atm1%t(j,i,k)/sfs%psa(j,i) *  &
                        log((sigma(k+1)+cell)/(sigma(k)+cell))
                end do
              end do
            end do
          end if
        else if ( idynamic == 2 ) then
          if ( associated(atm_ph_out) ) then
            atm_ph_out(:,:,1) = ptop*d_1000
            do k = 2 , kz
              do i = ici1 , ici2
                do j = jci1 , jci2
                  atm_ph_out(j,i,k) = atm0%pf(j,i,k) + &
                       d_half*(atm1%pp(j,i,k-1)+atm1%pp(j,i,k))/sfs%psa(j,i)
                end do
              end do
            end do
          end if
          if ( associated(atm_zh_out) ) then
            do k = 1 , kz
              do i = ici1 , ici2
                do j = jci1 , jci2
                  atm_zh_out(j,i,k) = atms%za(j,i,k)
                end do
              end do
            end do
          end if
        else if ( idynamic == 3 ) then
          if ( associated(atm_ph_out) ) then
            do k = 1 , kz
              do i = ici1 , ici2
                do j = jci1 , jci2
                  atm_ph_out(j,i,k) = (mo_atm%pai(j,i,k)**cpovr) * p00
                end do
              end do
            end do
          end if
          if ( associated(atm_zh_out) ) then
            do k = 1 , kz
              do i = ici1 , ici2
                do j = jci1 , jci2
                  atm_zh_out(j,i,k) = mo_atm%zeta(j,i,k)
                end do
              end do
            end do
          end if
        end if
        if ( associated(atm_zf_out) ) then
          do k = 1 , kz
            do i = ici1 , ici2
              do j = jci1 , jci2
                atm_zf_out(j,i,k) = atms%zq(j,i,k)
              end do
            end do
          end do
        end if
        if ( ipptls == 2 ) then
          if ( associated(atm_rainls_out) ) then
            do k = 1 , kz
              atm_rainls_out(:,:,k) = rain_ls(jci1:jci2,ici1:ici2,k)
            end do
          end if
        end if
        if ( associated(atm_raincc_out) ) then
          do k = 1 , kz
            atm_raincc_out(:,:,k) = rain_cc(jci1:jci2,ici1:ici2,k)
          end do
        end if
        if ( associated(atm_q_detr_out) ) then
          do k = 1 , kz
            atm_q_detr_out(:,:,k) = q_detr(jci1:jci2,ici1:ici2,k)
          end do
        end if

#ifdef DEBUG
        if ( ipptls == 2 .and. stats ) then
          if ( associated(atm_stats_supw_out) ) then
            do k = 1 , kz
              atm_stats_supw_out(:,:,k) = ngs%statssupw(jci1:jci2,ici1:ici2,k)
            end do
          end if
          if ( associated(atm_stats_supc_out) ) then
            do k = 1 , kz
              atm_stats_supc_out(:,:,k) = ngs%statssupc(jci1:jci2,ici1:ici2,k)
            end do
          end if
          if ( associated(atm_stats_detw_out) ) then
            do k = 1 , kz
              atm_stats_detw_out(:,:,k) = ngs%statsdetrw(jci1:jci2,ici1:ici2,k)
            end do
          end if
          if ( associated(atm_stats_detc_out) ) then
            do k = 1 , kz
              atm_stats_detc_out(:,:,k) = ngs%statsdetrc(jci1:jci2,ici1:ici2,k)
            end do
          end if
          if ( associated(atm_stats_erow_out) ) then
            do k = 1 , kz
              atm_stats_erow_out(:,:,k) = ngs%statserosw(jci1:jci2,ici1:ici2,k)
            end do
          end if
          if ( associated(atm_stats_eroc_out) ) then
            do k = 1 , kz
              atm_stats_eroc_out(:,:,k) = ngs%statserosc(jci1:jci2,ici1:ici2,k)
            end do
          end if
          if ( associated(atm_stats_evw_out) ) then
            do k = 1 , kz
              atm_stats_evw_out(:,:,k) = ngs%statsevapw(jci1:jci2,ici1:ici2,k)
            end do
          end if
          if ( associated(atm_stats_evc_out) ) then
            do k = 1 , kz
              atm_stats_evc_out(:,:,k) = ngs%statsevapc(jci1:jci2,ici1:ici2,k)
            end do
          end if
          if ( associated(atm_stats_con1w_out) ) then
            do k = 1 , kz
              atm_stats_con1w_out(:,:,k) = &
                             ngs%statscond1w(jci1:jci2,ici1:ici2,k)
            end do
          end if
          if ( associated(atm_stats_con1c_out) ) then
            do k = 1 , kz
              atm_stats_con1c_out(:,:,k) = &
                             ngs%statscond1c(jci1:jci2,ici1:ici2,k)
            end do
          end if
          if ( associated(atm_stats_dep_out) ) then
            do k = 1 , kz
              atm_stats_dep_out(:,:,k) = ngs%statsdepos(jci1:jci2,ici1:ici2,k)
            end do
          end if
          if ( associated(atm_stats_melt_out) ) then
            do k = 1 , kz
              atm_stats_melt_out(:,:,k) = ngs%statsmelt(jci1:jci2,ici1:ici2,k)
            end do
          end if
          if ( associated(atm_stats_frz_out) ) then
            do k = 1 , kz
              atm_stats_frz_out(:,:,k) = ngs%statsfrz(jci1:jci2,ici1:ici2,k)
            end do
          end if
          if ( associated(atm_stats_rainev_out) ) then
            do k = 1 , kz
              atm_stats_rainev_out(:,:,k) = &
                             ngs%statsrainev(jci1:jci2,ici1:ici2,k)
            end do
          end if
          if ( associated(atm_stats_snowev_out) ) then
            do k = 1 , kz
              atm_stats_snowev_out(:,:,k) = &
                             ngs%statssnowev(jci1:jci2,ici1:ici2,k)
            end do
          end if
          if ( associated(atm_stats_autocw_out) ) then
            do k = 1 , kz
              atm_stats_autocw_out(:,:,k) = &
                             ngs%statsautocvw(jci1:jci2,ici1:ici2,k)
            end do
          end if
          if ( associated(atm_stats_autocc_out) ) then
            do k = 1 , kz
              atm_stats_autocc_out(:,:,k) = &
                             ngs%statsautocvc(jci1:jci2,ici1:ici2,k)
            end do
          end if
        end if
#endif

        if ( ibltyp == 2 ) then
          if ( associated(atm_tke_out) ) then
            if ( idynamic == 3 ) then
              atm_tke_out = mo_atm%tke(jci1:jci2,ici1:ici2,1:kz)
            else
              atm_tke_out = atm1%tke(jci1:jci2,ici1:ici2,1:kz)
            end if
          end if
          if ( associated(atm_kth_out) ) &
            atm_kth_out = uwstate%kth(jci1:jci2,ici1:ici2,1:kz)
          if ( associated(atm_kzm_out) ) &
            atm_kzm_out = uwstate%kzm(jci1:jci2,ici1:ici2,1:kz)
        else if ( ibltyp == 4 ) then
          if ( associated(atm_tke_out) ) &
            atm_tke_out = atms%tkepbl(jci1:jci2,ici1:ici2,1:kz)
        end if

        if ( ichem == 1 .and. iaerosol == 1 .and. iindirect == 2 ) then
          if ( associated(atm_ccnnum_out) ) then
            do k = 1 , kz
              ! convert to 1/cm3
              atm_ccnnum_out(:,:,k) = ccn(jci1:jci2,ici1:ici2,k) * 1.e-6_rkx
            end do
          end if
          if ( idiag > 0 ) then
            if ( associated(atm_qcrit_out) ) then
              do k = 1 , kz
                atm_qcrit_out(:,:,k) = qdiag%qcr(jci1:jci2,ici1:ici2,k)
              end do
            end if
            if ( associated(atm_qincl_out) ) then
              do k = 1 , kz
                atm_qincl_out(:,:,k) = qdiag%qcl(jci1:jci2,ici1:ici2,k)
              end do
            end if
            if ( associated(atm_autoconvr_out) ) then
              do k = 1 , kz
                atm_autoconvr_out(:,:,k) = qdiag%acr(jci1:jci2,ici1:ici2,k)
              end do
            end if
          end if
        end if

        if ( associated(atm_tpr_out) ) then
          atm_tpr_out = (sfs%rainc+sfs%rainnc)/(atmfrq*secph)
        end if
        if ( associated(atm_tsn_out) ) then
          atm_tsn_out = sfs%snownc/(atmfrq*secph)
        end if
        if ( associated(atm_tgb_out) .and. rcmtimer%lcount == 0 ) then
          atm_tgb_out = sfs%tgbb(jci1:jci2,ici1:ici2)
        end if

        if ( associated(atm_tsw_out) ) then
          if ( rcmtimer%integrating( ) ) then
            where ( mddom%ldmsk == 1 )
              atm_tsw_out = atm_tsw_out / rnsrf_for_atmfrq
            elsewhere
              atm_tsw_out = dmissval
            end where
          else
            atm_tsw_out = dmissval
          end if
        end if

        if ( associated(atm_cape_out) .and. associated(atm_cin_out) ) then
          if ( idynamic == 3 ) then
            do i = ici1 , ici2
              do j = jci1 , jci2
                do k = 1 , kz
                  kk = kzp1 - k
                  p1d(kk) = mo_atm%p(j,i,k)
                  t1d(kk) = mo_atm%t(j,i,k)
                  rh1d(kk) = min(d_one,max(d_zero,(mo_atm%qx(j,i,k,iqv) / &
                      pfwsat(mo_atm%t(j,i,k),mo_atm%p(j,i,k)))))
                end do
                call getcape(kz,p1d,t1d,rh1d,atm_cape_out(j,i),atm_cin_out(j,i))
              end do
            end do
          else
            do i = ici1 , ici2
              do j = jci1 , jci2
                do k = 1 , kz
                  kk = kzp1 - k
                  p1d(kk) = atm1%pr(j,i,k)
                  t1d(kk) = atm1%t(j,i,k)/sfs%psa(j,i)
                  rh1d(kk) = min(d_one,max(d_zero, &
                     (atm1%qx(j,i,k,iqv)/ps_out(j,i)) / &
                     pfwsat(atm1%t(j,i,k)/ps_out(j,i),atm1%pr(j,i,k))))
                end do
                call getcape(kz,p1d,t1d,rh1d,atm_cape_out(j,i),atm_cin_out(j,i))
              end do
            end do
          end if
        end if

        ! FAB add tendency diagnostic here
        if ( idiag > 0 ) then
          if ( associated(atm_tten_adh_out) ) then
            if ( idynamic == 3 ) then
              do k = 1 , kz
                atm_tten_adh_out(:,:,k) = tdiag%adh(jci1:jci2,ici1:ici2,k)
              end do
            else
              do k = 1 , kz
                atm_tten_adh_out(:,:,k) = &
                   tdiag%adh(jci1:jci2,ici1:ici2,k)/ps_out
              end do
            end if
            tdiag%adh = d_zero
          end if
          if ( associated(atm_tten_adv_out) ) then
            if ( idynamic == 3 ) then
              do k = 1 , kz
                atm_tten_adv_out(:,:,k) = tdiag%adv(jci1:jci2,ici1:ici2,k)
              end do
            else
              do k = 1 , kz
                atm_tten_adv_out(:,:,k) = &
                   tdiag%adv(jci1:jci2,ici1:ici2,k)/ps_out
              end do
            end if
            tdiag%adv = d_zero
          end if
          if ( associated(atm_tten_tbl_out) ) then
            if ( idynamic == 3 ) then
              do k = 1 , kz
                atm_tten_tbl_out(:,:,k) = tdiag%tbl(jci1:jci2,ici1:ici2,k)
              end do
            else
              do k = 1 , kz
                atm_tten_tbl_out(:,:,k) = &
                   tdiag%tbl(jci1:jci2,ici1:ici2,k)/ps_out
              end do
            end if
            tdiag%tbl = d_zero
          end if
          if ( associated(atm_tten_dif_out) ) then
            if ( idynamic == 3 ) then
              do k = 1 , kz
                atm_tten_dif_out(:,:,k) = tdiag%dif(jci1:jci2,ici1:ici2,k)
              end do
            else
              do k = 1 , kz
                atm_tten_dif_out(:,:,k) = &
                   tdiag%dif(jci1:jci2,ici1:ici2,k)/ps_out
              end do
            end if
            tdiag%dif = d_zero
          end if
          if ( associated(atm_tten_bdy_out) ) then
            if ( idynamic == 3 ) then
              do k = 1 , kz
                atm_tten_bdy_out(:,:,k) = tdiag%bdy(jci1:jci2,ici1:ici2,k)
              end do
            else
              do k = 1 , kz
                atm_tten_bdy_out(:,:,k) = &
                   tdiag%bdy(jci1:jci2,ici1:ici2,k)/ps_out
              end do
            end if
            tdiag%bdy = d_zero
          end if
          if ( associated(atm_tten_con_out) ) then
            if ( idynamic == 3 ) then
              do k = 1 , kz
                atm_tten_con_out(:,:,k) = tdiag%con(jci1:jci2,ici1:ici2,k)
              end do
            else
              do k = 1 , kz
                atm_tten_con_out(:,:,k) = &
                   tdiag%con(jci1:jci2,ici1:ici2,k)/ps_out
              end do
            end if
            tdiag%con = d_zero
          end if
          if ( associated(atm_tten_adi_out) ) then
            if ( idynamic == 3 ) then
              do k = 1 , kz
                atm_tten_adi_out(:,:,k) = tdiag%adi(jci1:jci2,ici1:ici2,k)
              end do
            else
              do k = 1 , kz
                atm_tten_adi_out(:,:,k) = &
                   tdiag%adi(jci1:jci2,ici1:ici2,k)/ps_out
              end do
            end if
            tdiag%adi = d_zero
          end if
          if ( associated(atm_tten_rad_out) ) then
            if ( idynamic == 3 ) then
              do k = 1 , kz
                atm_tten_rad_out(:,:,k) = tdiag%rad(jci1:jci2,ici1:ici2,k)
              end do
            else
              do k = 1 , kz
                atm_tten_rad_out(:,:,k) = &
                   tdiag%rad(jci1:jci2,ici1:ici2,k)/ps_out
              end do
            end if
            tdiag%rad = d_zero
          end if
          if ( associated(atm_tten_lsc_out) ) then
            if ( idynamic == 3 ) then
              do k = 1 , kz
                atm_tten_lsc_out(:,:,k) = tdiag%lsc(jci1:jci2,ici1:ici2,k)
              end do
            else
              do k = 1 , kz
                atm_tten_lsc_out(:,:,k) = &
                   tdiag%lsc(jci1:jci2,ici1:ici2,k)/ps_out
              end do
            end if
            tdiag%lsc = d_zero
          end if
          if ( associated(atm_qten_adh_out) ) then
            if ( idynamic == 3 ) then
              do k = 1 , kz
                atm_qten_adh_out(:,:,k) = qdiag%adh(jci1:jci2,ici1:ici2,k)
              end do
            else
              do k = 1 , kz
                atm_qten_adh_out(:,:,k) = &
                   qdiag%adh(jci1:jci2,ici1:ici2,k)/ps_out
              end do
            end if
            qdiag%adh = d_zero
          end if
          if ( associated(atm_qten_adv_out) ) then
            if ( idynamic == 3 ) then
              do k = 1 , kz
                atm_qten_adv_out(:,:,k) = qdiag%adv(jci1:jci2,ici1:ici2,k)
              end do
            else
              do k = 1 , kz
                atm_qten_adv_out(:,:,k) = &
                   qdiag%adv(jci1:jci2,ici1:ici2,k)/ps_out
              end do
            end if
            qdiag%adv = d_zero
          end if
          if ( associated(atm_qten_tbl_out) ) then
            if ( idynamic == 3 ) then
              do k = 1 , kz
                atm_qten_tbl_out(:,:,k) = qdiag%tbl(jci1:jci2,ici1:ici2,k)
              end do
            else
              do k = 1 , kz
                atm_qten_tbl_out(:,:,k) = &
                   qdiag%tbl(jci1:jci2,ici1:ici2,k)/ps_out
              end do
            end if
            qdiag%tbl = d_zero
          end if
          if ( associated(atm_qten_dif_out) ) then
            if ( idynamic == 3 ) then
              do k = 1 , kz
                atm_qten_dif_out(:,:,k) = qdiag%dif(jci1:jci2,ici1:ici2,k)
              end do
            else
              do k = 1 , kz
                atm_qten_dif_out(:,:,k) = &
                  qdiag%dif(jci1:jci2,ici1:ici2,k)/ps_out
              end do
            end if
            qdiag%dif = d_zero
          end if
          if ( associated(atm_qten_bdy_out) ) then
            if ( idynamic == 3 ) then
              do k = 1 , kz
                atm_qten_bdy_out(:,:,k) = qdiag%bdy(jci1:jci2,ici1:ici2,k)
              end do
            else
              do k = 1 , kz
                atm_qten_bdy_out(:,:,k) = &
                   qdiag%bdy(jci1:jci2,ici1:ici2,k)/ps_out
              end do
            end if
            qdiag%bdy = d_zero
          end if
          if ( associated(atm_qten_con_out) ) then
            if ( idynamic == 3 ) then
              do k = 1 , kz
                atm_qten_con_out(:,:,k) = qdiag%con(jci1:jci2,ici1:ici2,k)
              end do
            else
              do k = 1 , kz
                atm_qten_con_out(:,:,k) = &
                   qdiag%con(jci1:jci2,ici1:ici2,k)/ps_out
              end do
            end if
            qdiag%con = d_zero
          end if
          if ( associated(atm_qten_adi_out) ) then
            if ( idynamic == 3 ) then
              do k = 1 , kz
                atm_qten_adi_out(:,:,k) = qdiag%adi(jci1:jci2,ici1:ici2,k)
              end do
            else
              do k = 1 , kz
                atm_qten_adi_out(:,:,k) = &
                   qdiag%adi(jci1:jci2,ici1:ici2,k)/ps_out
              end do
            end if
            qdiag%adi = d_zero
          end if
          if ( associated(atm_qten_rad_out) ) then
            if ( idynamic == 3 ) then
              do k = 1 , kz
                atm_qten_rad_out(:,:,k) = qdiag%rad(jci1:jci2,ici1:ici2,k)
              end do
            else
              do k = 1 , kz
                atm_qten_rad_out(:,:,k) = &
                  qdiag%rad(jci1:jci2,ici1:ici2,k)/ps_out
              end do
            end if
            qdiag%rad = d_zero
          end if
          if ( associated(atm_qten_lsc_out) ) then
            if ( idynamic == 3 ) then
              do k = 1 , kz
                atm_qten_lsc_out(:,:,k) = qdiag%lsc(jci1:jci2,ici1:ici2,k)
              end do
            else
              do k = 1 , kz
                atm_qten_lsc_out(:,:,k) = &
                   qdiag%lsc(jci1:jci2,ici1:ici2,k)/ps_out
              end do
            end if
            qdiag%lsc = d_zero
          end if
        end if

        if ( idynamic == 1 ) then
          ps_out = d_1000*(ps_out+ptop)
        else if ( idynamic == 2 ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              ps_out(j,i) = atm0%ps(j,i) + ptop*d_1000 + &
                         atm1%pp(j,i,kz)/sfs%psa(j,i)
            end do
          end do
        end if

        call write_record_output_stream(atm_stream,alarm_out_atm%idate)
        if ( myid == italk ) &
          write(stdout,*) 'ATM variables written at ' , rcmtimer%str( )

        if ( associated(atm_tsw_out) ) atm_tsw_out = d_zero
        sfs%rainc  = d_zero
        sfs%rainnc = d_zero
        if ( ipptls > 1 ) sfs%snownc  = d_zero
        rnsrf_for_atmfrq = d_zero
      end if
    end if

    if ( srf_stream > 0 ) then
      if ( ldosrf ) then
        srffac = d_one / rnsrf_for_srffrq
        srafac = d_one / rnrad_for_srffrq
        if ( idynamic == 1 ) then
          ps_out = d_1000*(sfs%psa(jci1:jci2,ici1:ici2)+ptop)
        else if ( idynamic == 2 ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              ps_out(j,i) = atm0%ps(j,i) + ptop*d_1000 + &
                             atm1%pp(j,i,kz)/sfs%psa(j,i)
            end do
          end do
        else
          ps_out = sfs%psa(jci1:jci2,ici1:ici2)
        end if
        ! Averaged values
        if ( associated(srf_tpr_out) ) &
          srf_tpr_out = srf_tpr_out*srffac
        if ( associated(srf_prcv_out) ) &
          srf_prcv_out = srf_prcv_out*srffac
        if ( associated(srf_zpbl_out) ) &
          srf_zpbl_out = srf_zpbl_out*srffac
        if ( associated(srf_dew_out) .and. associated(srf_evp_out) ) then
          srf_dew_out = -(srf_evp_out*srffac)
          srf_dew_out = max(srf_dew_out, d_zero)
        end if
        if ( associated(srf_evp_out) ) then
          srf_evp_out = srf_evp_out*srffac
          srf_evp_out = max(srf_evp_out, d_zero)
        end if
        if ( associated(srf_scv_out) ) then
          where ( mddom%ldmsk > 0 )
            srf_scv_out = srf_scv_out*srffac
          elsewhere
            srf_scv_out = dmissval
          end where
        end if
        if ( associated(srf_srunoff_out) ) then
          where ( srf_srunoff_out < dlowval )
            srf_srunoff_out = d_zero
          end where
          where ( mddom%ldmsk == 1 )
            srf_srunoff_out = srf_srunoff_out*srffac
          elsewhere
            srf_srunoff_out = dmissval
          end where
        end if
        if ( associated(srf_trunoff_out) ) then
          where ( mddom%ldmsk == 1 )
            where ( srf_trunoff_out < dlowval )
              srf_trunoff_out = d_zero
            end where
            srf_trunoff_out = srf_trunoff_out*srffac
          elsewhere
            srf_trunoff_out = dmissval
          end where
        end if
        if ( associated(srf_sena_out) ) &
          srf_sena_out = srf_sena_out*srffac
        if ( associated(srf_lena_out) ) &
          srf_lena_out = srf_lena_out*srffac
        if ( associated(srf_flw_out) ) &
          srf_flw_out = srf_flw_out*srffac
        if ( associated(srf_fsw_out) ) &
          srf_fsw_out = srf_fsw_out*srffac
        if ( associated(srf_fld_out) ) &
          srf_fld_out = srf_fld_out*srffac
        if ( associated(srf_sina_out) ) &
          srf_sina_out = srf_sina_out*srffac
        if ( associated(srf_uflw_out) ) &
          srf_uflw_out = srf_fld_out - srf_flw_out
        if ( associated(srf_ufsw_out) ) &
          srf_ufsw_out = srf_sina_out - srf_fsw_out
        if ( associated(srf_taux_out) .and. associated(srf_tauy_out) ) then
          srf_taux_out = srf_taux_out*srffac
          srf_tauy_out = srf_tauy_out*srffac
          if ( uvrotate ) then
            call uvrot(srf_taux_out,srf_tauy_out)
          end if
        end if
        if ( associated(srf_snowmelt_out) ) &
          srf_snowmelt_out = srf_snowmelt_out*srffac / alarm_out_srf%dt
        if ( associated(srf_u10m_out) .and. associated(srf_v10m_out) ) then
          if ( associated(srf_wspd_out) ) then
            srf_wspd_out = sqrt(srf_u10m_out(:,:,1)*srf_u10m_out(:,:,1) + &
                                srf_v10m_out(:,:,1)*srf_v10m_out(:,:,1))
          end if
          if ( uvrotate ) then
            call uvrot(srf_u10m_out,srf_v10m_out)
          end if
        end if

        if ( associated(srf_totcf_out) ) then
          srf_totcf_out = srf_totcf_out * srffac * d_100
        end if
        if ( associated(srf_evpot_out) ) then
          srf_evpot_out = srf_evpot_out * srffac
        end if

        if ( associated(srf_ua100_out) .and. &
             associated(srf_va100_out) ) then
          if ( idynamic == 3 ) then
            call uvstagtox(mo_atm%u,mo_atm%v,mo_atm%ux,mo_atm%vx)
            do i = ici1 , ici2
              do j = jci1 , jci2
                zz = mo_atm%zeta(j,i,kz)
                if ( zz > 100.0_rkx ) then
                  srf_ua100_out(j,i,1) = mo_atm%ux(j,i,kz)
                  srf_va100_out(j,i,1) = mo_atm%vx(j,i,kz)
                else
                  vloop1: &
                  do k = kz-1 , 1 , -1
                    zz1 = mo_atm%zeta(j,i,k)
                    if ( zz1 > 100.0_rkx ) then
                      ww = (100.0_rkx-zz)/(zz1-zz)
                      srf_ua100_out(j,i,1) = &
                        ww*mo_atm%ux(j,i,k)+(d_one-ww)*mo_atm%ux(j,i,k+1)
                      srf_va100_out(j,i,1) = &
                        ww*mo_atm%vx(j,i,k)+(d_one-ww)*mo_atm%vx(j,i,k+1)
                      exit vloop1
                    end if
                    zz = zz1
                  end do vloop1
                end if
              end do
            end do
          else if ( idynamic == 2 ) then
            do i = ici1 , ici2
              do j = jci1 , jci2
                zz = atm0%z(j,i,kz)
                if ( zz > 100.0_rkx ) then
                  srf_ua100_out(j,i,1) = &
                    (d_rfour*(atm1%u(j,i,kz)+atm1%u(j+1,i,kz) + &
                              atm1%u(j,i+1,kz)+atm1%u(j+1,i+1,kz)) / &
                              sfs%psa(j,i))
                  srf_va100_out(j,i,1) = &
                    (d_rfour*(atm1%v(j,i,kz)+atm1%v(j+1,i,kz) + &
                              atm1%v(j,i+1,kz)+atm1%v(j+1,i+1,kz)) / &
                              sfs%psa(j,i))
                else
                  vloop2: &
                  do k = kz-1 , 1 , -1
                    zz1 = atm0%z(j,i,k)
                    if ( zz1 > 100.0_rkx ) then
                      ww = (100.0_rkx-zz)/(zz1-zz)
                      srf_ua100_out(j,i,1) = &
                        ww * &
                           (d_rfour*(atm1%u(j,i,k)+atm1%u(j+1,i,k) + &
                                     atm1%u(j,i+1,k)+atm1%u(j+1,i+1,k)) / &
                                     sfs%psa(j,i)) + &
                        (d_one - ww) * &
                           (d_rfour*(atm1%u(j,i,k+1)+atm1%u(j+1,i,k+1) + &
                                     atm1%u(j,i+1,k+1)+atm1%u(j+1,i+1,k+1)) / &
                                     sfs%psa(j,i))
                      srf_va100_out(j,i,1) = &
                        ww * &
                           (d_rfour*(atm1%v(j,i,k)+atm1%v(j+1,i,k) + &
                                     atm1%v(j,i+1,k)+atm1%v(j+1,i+1,k)) / &
                                     sfs%psa(j,i)) + &
                        (d_one - ww) * &
                           (d_rfour*(atm1%v(j,i,k+1)+atm1%v(j+1,i,k+1) + &
                                     atm1%v(j,i+1,k+1)+atm1%v(j+1,i+1,k+1)) / &
                                     sfs%psa(j,i))
                      exit vloop2
                    end if
                    zz = zz1
                  end do vloop2
                end if
              end do
            end do
          else
            do i = ici1 , ici2
              do j = jci1 , jci2
                cell = ptop / sfs%psa(j,i)
                tv = atm1%t(j,i,kz)/sfs%psa(j,i) * &
                            (d_one + ep1*atm1%qx(j,i,kz,iqv)/sfs%psa(j,i))
                zz = rovg * tv * log((sigma(kzp1)+cell)/(sigma(kz)+cell))
                if ( zz > 100.0_rkx ) then
                  srf_ua100_out(j,i,1) = &
                    (d_rfour*(atm1%u(j,i,kz)+atm1%u(j+1,i,kz) + &
                              atm1%u(j,i+1,kz)+atm1%u(j+1,i+1,kz)) / &
                              sfs%psa(j,i))
                  srf_va100_out(j,i,1) = &
                    (d_rfour*(atm1%v(j,i,kz)+atm1%v(j+1,i,kz) + &
                              atm1%v(j,i+1,kz)+atm1%v(j+1,i+1,kz)) / &
                              sfs%psa(j,i))
                else
                  vloop3: &
                  do k = kz-1 , 1 , -1
                    tv = atm1%t(j,i,k)/sfs%psa(j,i) * &
                            (d_one + ep1*atm1%qx(j,i,k,iqv)/sfs%psa(j,i))
                    zz1 = zz + rovg*tv*log((sigma(k+1)+cell)/(sigma(k)+cell))
                    if ( zz1 > 100.0_rkx ) then
                      ww = (100.0_rkx-zz)/(zz1-zz)
                      srf_ua100_out(j,i,1) = &
                        ww * &
                           (d_rfour*(atm1%u(j,i,k)+atm1%u(j+1,i,k) + &
                                     atm1%u(j,i+1,k)+atm1%u(j+1,i+1,k)) / &
                                     sfs%psa(j,i)) + &
                        (d_one - ww) * &
                           (d_rfour*(atm1%u(j,i,k+1)+atm1%u(j+1,i,k+1) + &
                                     atm1%u(j,i+1,k+1)+atm1%u(j+1,i+1,k+1)) / &
                                     sfs%psa(j,i))
                      srf_va100_out(j,i,1) = &
                        ww * &
                           (d_rfour*(atm1%v(j,i,k)+atm1%v(j+1,i,k) + &
                                     atm1%v(j,i+1,k)+atm1%v(j+1,i+1,k)) / &
                                     sfs%psa(j,i)) + &
                        (d_one - ww) * &
                           (d_rfour*(atm1%v(j,i,k+1)+atm1%v(j+1,i,k+1) + &
                                     atm1%v(j,i+1,k+1)+atm1%v(j+1,i+1,k+1)) / &
                                     sfs%psa(j,i))
                      exit vloop3
                    end if
                    zz = zz1
                  end do vloop3
                end if
              end do
            end do
          end if
          if ( uvrotate ) then
            call uvrot(srf_ua100_out,srf_va100_out)
          end if
        end if

        call write_record_output_stream(srf_stream,alarm_out_srf%idate)
        if ( myid == italk ) &
          write(stdout,*) 'SRF variables written at ' , rcmtimer%str( )

        if ( associated(srf_tpr_out) ) srf_tpr_out = d_zero
        if ( associated(srf_prcv_out) ) srf_prcv_out = d_zero
        if ( associated(srf_zpbl_out) ) srf_zpbl_out = d_zero
        if ( associated(srf_evp_out) ) srf_evp_out = d_zero
        if ( associated(srf_scv_out) ) srf_scv_out = d_zero
        if ( associated(srf_srunoff_out) ) srf_srunoff_out = d_zero
        if ( associated(srf_trunoff_out) ) srf_trunoff_out = d_zero
        if ( associated(srf_sena_out) ) srf_sena_out = d_zero
        if ( associated(srf_lena_out) ) srf_lena_out = d_zero
        if ( associated(srf_flw_out) ) srf_flw_out = d_zero
        if ( associated(srf_fsw_out) ) srf_fsw_out = d_zero
        if ( associated(srf_fld_out) ) srf_fld_out = d_zero
        if ( associated(srf_sina_out) ) srf_sina_out = d_zero
        if ( associated(srf_sund_out) ) srf_sund_out = d_zero
        if ( associated(srf_taux_out) ) srf_taux_out = d_zero
        if ( associated(srf_tauy_out) ) srf_tauy_out = d_zero
        if ( associated(srf_snowmelt_out) ) srf_snowmelt_out = d_zero
        if ( associated(srf_totcf_out) ) srf_totcf_out = d_zero
        if ( associated(srf_evpot_out) ) srf_evpot_out = d_zero

        rnsrf_for_srffrq = d_zero
        rnrad_for_srffrq = d_zero

      end if
    end if

    if ( sub_stream > 0 ) then
      if ( ldosub ) then
        subfac = d_one / rnsrf_for_subfrq
        sub_ps_out = sub_ps_out*subfac

        if ( associated(sub_evp_out) ) then
          sub_evp_out = sub_evp_out*subfac
          sub_evp_out = max(sub_evp_out, d_zero)
        end if
        if ( associated(sub_scv_out) ) then
          where ( sub_scv_out < dmissval )
            sub_scv_out = sub_scv_out*subfac
          end where
        end if
        if ( associated(sub_sena_out) ) &
          sub_sena_out = sub_sena_out*subfac
        if ( associated(sub_srunoff_out) ) then
          where ( sub_srunoff_out < dmissval )
            sub_srunoff_out = sub_srunoff_out*subfac
          end where
        end if
        if ( associated(sub_trunoff_out) ) then
          where ( sub_trunoff_out < dmissval )
            sub_trunoff_out = sub_trunoff_out*subfac
          end where
        end if

        if ( associated(sub_u10m_out) .and. associated(sub_v10m_out) ) then
          if ( uvrotate ) then
            call uvrot(sub_u10m_out,sub_v10m_out)
          end if
        end if
        call write_record_output_stream(sub_stream,alarm_out_sub%idate)
        if ( myid == italk ) &
          write(stdout,*) 'SUB variables written at ' , rcmtimer%str( )

        if ( associated(sub_evp_out) ) sub_evp_out = d_zero
        if ( associated(sub_scv_out) ) sub_scv_out = d_zero
        if ( associated(sub_sena_out) ) sub_sena_out = d_zero
        if ( associated(sub_srunoff_out) ) sub_srunoff_out = d_zero
        if ( associated(sub_trunoff_out) ) sub_trunoff_out = d_zero
        rnsrf_for_subfrq = d_zero
      end if
    end if

    if ( lak_stream > 0 ) then
      if ( ldolak ) then
        lakfac = d_one / rnsrf_for_lakfrq

        if ( idynamic == 1 ) then
          ps_out = d_1000*(sfs%psa(jci1:jci2,ici1:ici2)+ptop)
        else if ( idynamic == 2 ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              ps_out(j,i) = atm0%ps(j,i) + ptop*d_1000 + &
                     atm1%pp(j,i,kz) / sfs%psa(j,i)
            end do
          end do
        else
          ps_out = sfs%psa(jci1:jci2,ici1:ici2)
        end if
        if ( associated(lak_tpr_out) ) &
          lak_tpr_out = lak_tpr_out*lakfac
        if ( associated(lak_scv_out) ) then
          where ( mddom%ldmsk > 0 )
            lak_scv_out = lak_scv_out*lakfac
          elsewhere
            lak_scv_out = dmissval
          end where
        end if
        if ( associated(lak_sena_out) ) &
          lak_sena_out = lak_sena_out*lakfac
        if ( associated(lak_flw_out) ) &
          lak_flw_out = lak_flw_out*lakfac
        if ( associated(lak_fsw_out) ) &
          lak_fsw_out = lak_fsw_out*lakfac
        if ( associated(lak_fld_out) ) &
          lak_fld_out = lak_fld_out*lakfac
        if ( associated(lak_sina_out) ) &
          lak_sina_out = lak_sina_out*lakfac
        if ( associated(lak_evp_out) ) then
          lak_evp_out = lak_evp_out*lakfac
          lak_evp_out = max(lak_evp_out, d_zero)
        end if

        call write_record_output_stream(lak_stream,alarm_out_lak%idate)
        if ( myid == italk ) &
          write(stdout,*) 'LAK variables written at ' , rcmtimer%str( )

        if ( associated(lak_tpr_out) )    lak_tpr_out = d_zero
        if ( associated(lak_scv_out) )    lak_scv_out = d_zero
        if ( associated(lak_sena_out) )   lak_sena_out = d_zero
        if ( associated(lak_flw_out) )    lak_flw_out = d_zero
        if ( associated(lak_fsw_out) )    lak_fsw_out = d_zero
        if ( associated(lak_fld_out) )    lak_fld_out = d_zero
        if ( associated(lak_sina_out) )   lak_sina_out = d_zero
        if ( associated(lak_evp_out) )    lak_evp_out = d_zero
        rnsrf_for_lakfrq = d_zero
      end if
    end if

    if ( opt_stream > 0 .and. rcmtimer%integrating( ) ) then
      if ( ldoopt ) then
        optfac = d_one / rnrad_for_optfrq
        if ( idynamic == 1 ) then
          ps_out = d_1000*(sfs%psa(jci1:jci2,ici1:ici2)+ptop)
        else if ( idynamic == 2 ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              ps_out(j,i) = atm0%ps(j,i) + ptop*d_1000 + &
                atm1%pp(j,i,kz)/sfs%psa(j,i)
            end do
          end do
        else
          ps_out = sfs%psa(jci1:jci2,ici1:ici2)
        end if
        if ( associated(opt_pp_out) ) then
          do k = 1 , kz
            opt_pp_out(:,:,k) = atm1%pp(jci1:jci2,ici1:ici2,k) / &
                                  sfs%psa(jci1:jci2,ici1:ici2)
          end do
        end if
        if ( associated(opt_pai_out) ) then
          do k = 1 , kz
            opt_pai_out(:,:,k) = mo_atm%pai(jci1:jci2,ici1:ici2,k)
          end do
        end if
        if ( associated(opt_acstoarf_out) ) &
          opt_acstoarf_out = opt_acstoarf_out * optfac
        if ( associated(opt_acstsrrf_out) ) &
          opt_acstsrrf_out = opt_acstsrrf_out * optfac
        if ( associated(opt_acstalrf_out) ) &
          opt_acstalrf_out = opt_acstalrf_out * optfac
        if ( associated(opt_acssrlrf_out) ) &
          opt_acssrlrf_out = opt_acssrlrf_out * optfac
        if ( associated(opt_aastoarf_out) ) &
          opt_aastoarf_out = opt_aastoarf_out * optfac
        if ( associated(opt_aastsrrf_out) ) &
          opt_aastsrrf_out = opt_aastsrrf_out * optfac
        if ( associated(opt_aastalrf_out) ) &
          opt_aastalrf_out = opt_aastalrf_out * optfac
        if ( associated(opt_aassrlrf_out) ) &
          opt_aassrlrf_out = opt_aassrlrf_out * optfac

        call write_record_output_stream(opt_stream,alarm_out_opt%idate)
        if ( myid == italk ) &
          write(stdout,*) 'OPT variables written at ' , rcmtimer%str( )
        if ( associated(opt_acstoarf_out) ) opt_acstoarf_out = d_zero
        if ( associated(opt_acstsrrf_out) ) opt_acstsrrf_out = d_zero
        if ( associated(opt_acstalrf_out) ) opt_acstalrf_out = d_zero
        if ( associated(opt_acssrlrf_out) ) opt_acssrlrf_out = d_zero
        if ( associated(opt_aastoarf_out) ) opt_aastoarf_out = d_zero
        if ( associated(opt_aastsrrf_out) ) opt_aastsrrf_out = d_zero
        if ( associated(opt_aastalrf_out) ) opt_aastalrf_out = d_zero
        if ( associated(opt_aassrlrf_out) ) opt_aassrlrf_out = d_zero
        rnrad_for_optfrq = d_zero
      end if
    end if

    if ( che_stream > 0 ) then
      if ( ldoche ) then
        if ( idynamic == 1 ) then
          ps_out = d_1000*(sfs%psa(jci1:jci2,ici1:ici2)+ptop)
        else if ( idynamic == 2 ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              ps_out(j,i) = atm0%ps(j,i) + ptop*d_1000 + &
                atm1%pp(j,i,kz)/sfs%psa(j,i)
            end do
          end do
        else
          ps_out = sfs%psa(jci1:jci2,ici1:ici2)
        end if
        if ( associated(che_pp_out) ) then
          do k = 1 , kz
            che_pp_out(:,:,k) = atm1%pp(jci1:jci2,ici1:ici2,k) / &
                                 sfs%psa(jci1:jci2,ici1:ici2)
          end do
        end if
        if ( associated(che_pai_out) ) then
          do k = 1 , kz
            che_pai_out(:,:,k) = mo_atm%pai(jci1:jci2,ici1:ici2,k)
          end do
        end if
        do itr = 1 , ntr
          call fill_chem_outvars(itr)
          call write_record_output_stream(che_stream,alarm_out_che%idate,itr)
        end do
        if ( myid == italk ) then
          write(stdout,*) 'CHE variables written at ' , rcmtimer%str( )
        end if
      end if
    end if

    if ( shf_stream > 0 ) then
      if ( ldoshf ) then
        if ( idynamic == 1 ) then
          ps_out = d_1000*(sfs%psa(jci1:jci2,ici1:ici2)+ptop)
        else if ( idynamic == 2 ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              ps_out(j,i) = atm0%ps(j,i) + ptop*d_1000 + &
                atm1%pp(j,i,kz)/sfs%psa(j,i)
            end do
          end do
        else
          ps_out = sfs%psa(jci1:jci2,ici1:ici2)
        end if

        if ( associated(shf_pcpavg_out) ) &
          shf_pcpavg_out = shf_pcpavg_out * dtsrf/3600.0_rkx
        if ( associated(shf_pcprcv_out) ) &
          shf_pcprcv_out = shf_pcprcv_out * dtsrf/3600.0_rkx
        call write_record_output_stream(shf_stream,alarm_out_shf%idate)
        if ( myid == italk ) &
          write(stdout,*) 'SHF variables written at ' , rcmtimer%str( )

        if ( associated(shf_pcpavg_out) ) shf_pcpavg_out = d_zero
        if ( associated(shf_pcpmax_out) ) shf_pcpmax_out = d_zero
        if ( associated(shf_pcprcv_out) ) shf_pcprcv_out = d_zero
        if ( associated(shf_twetb_out) ) shf_twetb_out = d_zero
      end if
    end if

    if ( sts_stream > 0 ) then
      if ( ldosts ) then
        stsfac = d_one / rnsrf_for_day
        if ( idynamic == 1 ) then
          ps_out = d_1000*(sfs%psa(jci1:jci2,ici1:ici2)+ptop)
        else if ( idynamic == 2 ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              ps_out(j,i) = atm0%ps(j,i) + ptop*d_1000 + &
                atm1%pp(j,i,kz)/sfs%psa(j,i)
            end do
          end do
        else
          ps_out = sfs%psa(jci1:jci2,ici1:ici2)
        end if
        if ( associated(sts_pcpavg_out) ) &
          sts_pcpavg_out = sts_pcpavg_out*stsfac
        if ( associated(sts_t2avg_out) ) &
          sts_t2avg_out = sts_t2avg_out*stsfac
        if ( associated(sts_psavg_out) ) &
          sts_psavg_out = sts_psavg_out*stsfac
        if ( associated(sts_srunoff_out) ) then
          where ( mddom%ldmsk == 1 )
            sts_srunoff_out = sts_srunoff_out*stsfac
          elsewhere
            sts_srunoff_out = dmissval
          end where
        end if
        if ( associated(sts_trunoff_out) ) then
          where ( mddom%ldmsk == 1 )
            sts_trunoff_out = sts_trunoff_out*stsfac
          elsewhere
            sts_trunoff_out = dmissval
          end where
        end if

        call write_record_output_stream(sts_stream,alarm_out_sts%idate)
        if ( myid == italk ) &
          write(stdout,*) 'STS variables written at ' , rcmtimer%str( )

        if ( associated(sts_pcpavg_out) )  sts_pcpavg_out  = d_zero
        if ( associated(sts_t2avg_out) )   sts_t2avg_out   = d_zero
        if ( associated(sts_psavg_out) )   sts_psavg_out   = d_zero
        if ( associated(sts_tgmax_out) )   sts_tgmax_out   = -1.e30_rkx
        if ( associated(sts_tgmin_out) )   sts_tgmin_out   =  1.e30_rkx
        if ( associated(sts_t2max_out) )   sts_t2max_out   = -1.e30_rkx
        if ( associated(sts_t2min_out) )   sts_t2min_out   =  1.e30_rkx
        if ( associated(sts_w10max_out) )  sts_w10max_out  = -1.e30_rkx
        if ( associated(sts_psmin_out) )   sts_psmin_out   =  1.e30_rkx
        if ( associated(sts_pcpmax_out) )  sts_pcpmax_out  = -1.e30_rkx
        if ( associated(sts_sund_out) )    sts_sund_out    = d_zero
        if ( associated(sts_srunoff_out) ) sts_srunoff_out = d_zero
        if ( associated(sts_trunoff_out) ) sts_trunoff_out = d_zero
        rnsrf_for_day = d_zero
      end if
    end if

    if ( rad_stream > 0 ) then
      if ( ldorad ) then
        radfac = d_one / rnrad_for_radfrq
        if ( idynamic == 1 ) then
          ps_out = d_1000*(sfs%psa(jci1:jci2,ici1:ici2)+ptop)
        else if ( idynamic == 2 ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              ps_out(j,i) = atm0%ps(j,i) + ptop*d_1000 + &
                 atm1%pp(j,i,kz)/sfs%psa(j,i)
            end do
          end do
        else
          ps_out = sfs%psa(jci1:jci2,ici1:ici2)
        end if
        if ( associated(rad_pp_out) ) then
          do k = 1 , kz
            rad_pp_out(:,:,k) = atm1%pp(jci1:jci2,ici1:ici2,k)/ &
                             sfs%psa(jci1:jci2,ici1:ici2)
          end do
        end if
        if ( associated(rad_pai_out) ) then
          do k = 1 , kz
            rad_pai_out(:,:,k) = mo_atm%pai(jci1:jci2,ici1:ici2,k)
          end do
        end if
        if ( associated(rad_higcl_out) ) &
          rad_higcl_out = min(rad_higcl_out * radfac,d_one) * d_100
        if ( associated(rad_midcl_out) ) &
          rad_midcl_out = min(rad_midcl_out * radfac,d_one) * d_100
        if ( associated(rad_lowcl_out) ) &
          rad_lowcl_out = min(rad_lowcl_out * radfac,d_one) * d_100
        call write_record_output_stream(rad_stream,alarm_out_rad%idate)
        if ( myid == italk ) &
          write(stdout,*) 'RAD variables written at ' , rcmtimer%str( )
        if ( associated(rad_higcl_out) ) rad_higcl_out = d_zero
        if ( associated(rad_midcl_out) ) rad_midcl_out = d_zero
        if ( associated(rad_lowcl_out) ) rad_lowcl_out = d_zero
        rnrad_for_radfrq = d_zero
      end if
    end if

    if ( slaboc_stream > 0 ) then
      if ( ldoslab ) then
        if ( idynamic == 1 ) then
          ps_out = d_1000*(sfs%psa(jci1:jci2,ici1:ici2)+ptop)
        else if ( idynamic == 2 ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              ps_out(j,i) = atm0%ps(j,i) + ptop*d_1000 + &
                atm1%pp(j,i,kz)/sfs%psa(j,i)
            end do
          end do
        else
          ps_out = sfs%psa(jci1:jci2,ici1:ici2)
        end if
        call fill_slaboc_outvars
        call writevar_output_stream(slaboc_stream,v3dvar_slaboc(slab_qflx))
        if ( myid == italk ) then
          write(stdout,*) 'SOM variables written at ' , rcmtimer%str( )
        end if
      end if
    end if

    if ( ifsave ) then
      if ( ldosav ) then
        if ( do_parallel_save ) then
          if ( idynamic == 3 ) then
            atm_u_io(jde1:jde2,ice1:ice2,:) = mo_atm%u(jde1:jde2,ice1:ice2,:)
            atm_v_io(jce1:jce2,ide1:ide2,:) = mo_atm%v(jce1:jce2,ide1:ide2,:)
            atm_w_io(jce1:jce2,ice1:ice2,:) = mo_atm%w(jce1:jce2,ice1:ice2,:)
            atm_t_io(jce1:jce2,ice1:ice2,:) = mo_atm%t(jce1:jce2,ice1:ice2,:)
            atm_pai_io(jce1:jce2,ice1:ice2,:) = &
                              mo_atm%pai(jce1:jce2,ice1:ice2,:)
            atm_qx_io(jce1:jce2,ice1:ice2,:,:) = &
                               mo_atm%qx(jce1:jce2,ice1:ice2,:,:)
            if ( ibltyp == 2 ) then
              atm_tke_io(jce1:jce2,ice1:ice2,:) = &
                       mo_atm%tke(jce1:jce2,ice1:ice2,:)
            end if
            if ( ichem == 1 ) then
              trac_io(jce1:jce2,ice1:ice2,:,:) = &
                                mo_atm%trac(jce1:jce2,ice1:ice2,:,:)
            end if
            ps_io(jce1:jce2,ice1:ice2) = sfs%psa(jce1:jce2,ice1:ice2)
          else
            atm1_u_io(jde1:jde2,ide1:ide2,:) = atm1%u(jde1:jde2,ide1:ide2,:)
            atm1_v_io(jde1:jde2,ide1:ide2,:) = atm1%v(jde1:jde2,ide1:ide2,:)
            atm1_t_io(jce1:jce2,ice1:ice2,:) = atm1%t(jce1:jce2,ice1:ice2,:)
            atm1_qx_io(jce1:jce2,ice1:ice2,:,:) = &
                               atm1%qx(jce1:jce2,ice1:ice2,:,:)
            atm2_u_io(jde1:jde2,ide1:ide2,:) = atm2%u(jde1:jde2,ide1:ide2,:)
            atm2_v_io(jde1:jde2,ide1:ide2,:) = atm2%v(jde1:jde2,ide1:ide2,:)
            atm2_t_io(jce1:jce2,ice1:ice2,:) = atm2%t(jce1:jce2,ice1:ice2,:)
            atm2_qx_io(jce1:jce2,ice1:ice2,:,:) = &
                               atm2%qx(jce1:jce2,ice1:ice2,:,:)
            if ( ibltyp == 2 ) then
              atm1_tke_io(jce1:jce2,ice1:ice2,:) = &
                       atm1%tke(jce1:jce2,ice1:ice2,:)
              atm2_tke_io(jce1:jce2,ice1:ice2,:) = &
                       atm2%tke(jce1:jce2,ice1:ice2,:)
            end if
            if ( ichem == 1 ) then
              chia_io(jce1:jce2,ice1:ice2,:,:) = &
                                atm1%chi(jce1:jce2,ice1:ice2,:,:)
              chib_io(jce1:jce2,ice1:ice2,:,:) = &
                                atm2%chi(jce1:jce2,ice1:ice2,:,:)
            end if
            psa_io(jce1:jce2,ice1:ice2) = sfs%psa(jce1:jce2,ice1:ice2)
            psb_io(jce1:jce2,ice1:ice2) = sfs%psb(jce1:jce2,ice1:ice2)
          end if
          if ( ibltyp == 2 ) then
            kpbl_io = kpbl
          end if
          if ( ibltyp == 4 ) then
            tke_pbl_io = atms%tkepbl
            kpbl_io = kpbl
            myjsf_uz0_io = sfs%uz0
            myjsf_vz0_io = sfs%vz0
            myjsf_thz0_io = sfs%thz0
            myjsf_qz0_io = sfs%qz0
          end if
          if ( idynamic == 2 ) then
            atm1_pp_io(jce1:jce2,ice1:ice2,:) = atm1%pp(jce1:jce2,ice1:ice2,:)
            atm2_pp_io(jce1:jce2,ice1:ice2,:) = atm2%pp(jce1:jce2,ice1:ice2,:)
            atm1_w_io(jce1:jce2,ice1:ice2,:) = atm1%w(jce1:jce2,ice1:ice2,:)
            atm2_w_io(jce1:jce2,ice1:ice2,:) = atm2%w(jce1:jce2,ice1:ice2,:)
          end if
          hfx_io = sfs%hfx
          qfx_io = sfs%qfx
          tgbb_io = sfs%tgbb
          zo_io = sfs%zo
          uvdrag_io = sfs%uvdrag
          ram_io = sfs%ram1
          rah_io = sfs%rah1
          br_io = sfs%br
          q2m_io = sfs%q2m
          u10m_io = sfs%u10m
          v10m_io = sfs%v10m
          w10m_io = sfs%w10m
          ustar_io = sfs%ustar
          if ( ipptls > 0 ) then
            fcc_io = fcc
          end if
          heatrt_io = heatrt
          o3prof_io = o3prof
          if ( iocnflx == 2 ) then
            zpbl_io = zpbl
          end if
          if ( any(icup == 3) ) then
            cldefi_io = cldefi
          end if
          if ( any(icup == 4) ) then
            cbmf2d_io = cbmf2d
          end if
          if ( any(icup == 6) .or. any(icup == 5) ) then
            cu_avg_ww_io = avg_ww
          end if
          if ( irrtm == 0 ) then
            gasabsnxt_io = gasabsnxt
            gasabstot_io = gasabstot
            gasemstot_io = gasemstot
          end if
          sw_io = lms%sw
#ifdef CLM45
          if ( ichem == 1 ) then
            tsoi_io = tsoi
            swvol_io = sw_vol
          end if
#else
          gwet_io = lms%gwet
          ldew_io = lms%ldew
          taf_io = lms%taf
#endif
          tgrd_io = lms%tgrd
          tgbrd_io = lms%tgbrd
          tlef_io = lms%tlef
          sncv_io = lms%sncv
          sfice_io = lms%sfice
          snag_io = lms%snag
          emisv_io = lms%emisv
          um10_io = lms%um10
          swalb_io = lms%swalb
          lwalb_io = lms%lwalb
          swdiralb_io = lms%swdiralb
          swdifalb_io = lms%swdifalb
          lwdiralb_io = lms%lwdiralb
          lwdifalb_io = lms%lwdifalb
          ldmsk1_io(:,jci1:jci2,ici1:ici2) = mdsub%ldmsk(:,jci1:jci2,ici1:ici2)
          solis_io = solis
          solvs_io = solvs
          solvsd_io = solvsd
          solvl_io = solvl
          solvld_io = solvld
          sabveg_io = sabveg
          flw_io = flw
          flwd_io = flwd
          fsw_io = fsw
          sinc_io = sinc
          ldmsk_io = mddom%ldmsk
#ifndef CLM
          if ( lakemod == 1 ) then
            eta_io = lms%eta
            hi_io = lms%hi
            tlak_io = lms%tlake
          end if
#else
          if ( imask == 2 ) then
            lndcat_io(jce1:jce2,ice1:ice2) = mddom%lndcat(jce1:jce2,ice1:ice2)
          end if
#endif
          if ( idcsst == 1 ) then
            sst_io = lms%sst
            tskin_io = lms%tskin
            deltas_io = lms%deltas
            tdeltas_io = lms%tdeltas
          end if
          if ( idynamic == 1 ) then
            dstor_io(jde1:jde2,ide1:ide2,:) = dstor(jde1:jde2,ide1:ide2,:)
            hstor_io(jde1:jde2,ide1:ide2,:) = hstor(jde1:jde2,ide1:ide2,:)
          end if
          if ( ichem == 1 ) then
            convpr_io = convpr
            rainout_io = rainout
            washout_io = washout
            remdrd_io = remdrd
            if ( igaschem == 1 .and. ichsolver > 0 ) then
              chemall_io = chemall
              taucldsp_io = taucldsp
            end if
            ssw2da_io = ssw2da
#ifdef CLM45
            duflux_io = dustflx_clm
            voflux_io = voc_em_clm
#else
            sdelt_io = sdelt
            sdelq_io = sdelq
            svegfrac2d_io = svegfrac2d
#endif
            sfracv2d_io = sfracv2d
            sfracb2d_io = sfracb2d
            sfracs2d_io = sfracs2d
          end if
          if ( islab_ocean == 1 .and. do_restore_sst ) then
            qflux_restore_sst_io = qflux_restore_sst
          end if
        else
          if ( idynamic == 3 ) then
            call grid_collect(mo_atm%u,atm_u_io,jde1,jde2,ice1,ice2,1,kz)
            call grid_collect(mo_atm%v,atm_v_io,jce1,jce2,ide1,ide2,1,kz)
            call grid_collect(mo_atm%w,atm_w_io,jce1,jce2,ice1,ice2,1,kzp1)
            call grid_collect(mo_atm%t,atm_t_io,jce1,jce2,ice1,ice2,1,kz)
            call grid_collect(mo_atm%pai,atm_pai_io,jce1,jce2,ice1,ice2,1,kz)
            call grid_collect(mo_atm%qx,atm_qx_io, &
                              jce1,jce2,ice1,ice2,1,kz,1,nqx)
            if ( ibltyp == 2 ) then
              call grid_collect(mo_atm%tke,atm_tke_io, &
                                jce1,jce2,ice1,ice2,1,kzp1)
            end if
            if ( ichem == 1 ) then
              call grid_collect(mo_atm%trac,trac_io, &
                                jce1,jce2,ice1,ice2,1,kz,1,ntr)
            end if
            call grid_collect(sfs%psa,ps_io,jce1,jce2,ice1,ice2)
          else
            call grid_collect(atm1%u,atm1_u_io,jde1,jde2,ide1,ide2,1,kz)
            call grid_collect(atm1%v,atm1_v_io,jde1,jde2,ide1,ide2,1,kz)
            call grid_collect(atm1%t,atm1_t_io,jce1,jce2,ice1,ice2,1,kz)
            call grid_collect(atm1%qx,atm1_qx_io,jce1,jce2,ice1,ice2,1,kz,1,nqx)
            call grid_collect(atm2%u,atm2_u_io,jde1,jde2,ide1,ide2,1,kz)
            call grid_collect(atm2%v,atm2_v_io,jde1,jde2,ide1,ide2,1,kz)
            call grid_collect(atm2%t,atm2_t_io,jce1,jce2,ice1,ice2,1,kz)
            call grid_collect(atm2%qx,atm2_qx_io,jce1,jce2,ice1,ice2,1,kz,1,nqx)
            if ( ibltyp == 2 ) then
              call grid_collect(atm1%tke,atm1_tke_io,jce1,jce2,ice1,ice2,1,kzp1)
              call grid_collect(atm2%tke,atm2_tke_io,jce1,jce2,ice1,ice2,1,kzp1)
            end if
            if ( ichem == 1 ) then
              call grid_collect(atm1%chi,chia_io,jce1,jce2,ice1,ice2,1,kz,1,ntr)
              call grid_collect(atm2%chi,chib_io,jce1,jce2,ice1,ice2,1,kz,1,ntr)
            end if
            call grid_collect(sfs%psa,psa_io,jce1,jce2,ice1,ice2)
            call grid_collect(sfs%psb,psb_io,jce1,jce2,ice1,ice2)
          end if
          if ( ibltyp == 2 ) then
            call grid_collect(kpbl,kpbl_io,jci1,jci2,ici1,ici2)
          end if
          if ( ibltyp == 4 ) then
            call grid_collect(atms%tkepbl,tke_pbl_io,jci1,jci2,ici1,ici2,1,kz)
            call grid_collect(kpbl,kpbl_io,jci1,jci2,ici1,ici2)
            call grid_collect(sfs%uz0,myjsf_uz0_io,jci1,jci2,ici1,ici2)
            call grid_collect(sfs%vz0,myjsf_vz0_io,jci1,jci2,ici1,ici2)
            call grid_collect(sfs%thz0,myjsf_thz0_io,jci1,jci2,ici1,ici2)
            call grid_collect(sfs%qz0,myjsf_qz0_io,jci1,jci2,ici1,ici2)
          end if
          if ( idynamic == 2 ) then
            call grid_collect(atm1%pp,atm1_pp_io,jce1,jce2,ice1,ice2,1,kz)
            call grid_collect(atm2%pp,atm2_pp_io,jce1,jce2,ice1,ice2,1,kz)
            call grid_collect(atm1%w,atm1_w_io,jce1,jce2,ice1,ice2,1,kzp1)
            call grid_collect(atm2%w,atm2_w_io,jce1,jce2,ice1,ice2,1,kzp1)
          end if
          call grid_collect(sfs%hfx,hfx_io,jci1,jci2,ici1,ici2)
          call grid_collect(sfs%qfx,qfx_io,jci1,jci2,ici1,ici2)
          call grid_collect(sfs%tgbb,tgbb_io,jci1,jci2,ici1,ici2)
          call grid_collect(sfs%zo,zo_io,jci1,jci2,ici1,ici2)
          call grid_collect(sfs%uvdrag,uvdrag_io,jci1,jci2,ici1,ici2)
          call grid_collect(sfs%ram1,ram_io,jci1,jci2,ici1,ici2)
          call grid_collect(sfs%rah1,rah_io,jci1,jci2,ici1,ici2)
          call grid_collect(sfs%br,br_io,jci1,jci2,ici1,ici2)
          call grid_collect(sfs%q2m,q2m_io,jci1,jci2,ici1,ici2)
          call grid_collect(sfs%u10m,u10m_io,jci1,jci2,ici1,ici2)
          call grid_collect(sfs%v10m,v10m_io,jci1,jci2,ici1,ici2)
          call grid_collect(sfs%w10m,w10m_io,jci1,jci2,ici1,ici2)
          call grid_collect(sfs%ustar,ustar_io,jci1,jci2,ici1,ici2)
          if ( ipptls > 0 ) then
            call grid_collect(fcc,fcc_io,jci1,jci2,ici1,ici2,1,kz)
          end if
          call grid_collect(heatrt,heatrt_io,jci1,jci2,ici1,ici2,1,kz)
          call grid_collect(o3prof,o3prof_io,jci1,jci2,ici1,ici2,1,kzp1)
          if ( iocnflx == 2 ) then
            call grid_collect(zpbl,zpbl_io,jci1,jci2,ici1,ici2)
          end if
          if ( any(icup == 3) ) then
            call grid_collect(cldefi,cldefi_io,jci1,jci2,ici1,ici2)
          end if
          if ( any(icup == 4) ) then
            call grid_collect(cbmf2d,cbmf2d_io,jci1,jci2,ici1,ici2)
          end if
          if ( any(icup == 6) .or. any(icup == 5) ) then
            call grid_collect(avg_ww,cu_avg_ww_io,jci1,jci2,ici1,ici2,1,kz)
          end if
          if ( irrtm == 0 ) then
            call grid_collect(gasabsnxt,gasabsnxt_io, &
                              jci1,jci2,ici1,ici2,1,kz,1,4)
            call grid_collect(gasabstot,gasabstot_io, &
                              jci1,jci2,ici1,ici2,1,kzp1,1,kzp1)
            call grid_collect(gasemstot,gasemstot_io,jci1,jci2,ici1,ici2,1,kzp1)
          end if
          call subgrid_collect(lms%sw,sw_io,jci1,jci2,ici1,ici2, &
                               1,num_soil_layers)
#ifdef CLM45
          if ( ichem == 1 ) then
            call grid_collect(tsoi,tsoi_io,jci1,jci2,ici1,ici2, &
                              1,num_soil_layers)
            call grid_collect(sw_vol,swvol_io,jci1,jci2,ici1,ici2, &
                              1,num_soil_layers)
          end if
#else
          call subgrid_collect(lms%gwet,gwet_io,jci1,jci2,ici1,ici2)
          call subgrid_collect(lms%ldew,ldew_io,jci1,jci2,ici1,ici2)
          call subgrid_collect(lms%taf,taf_io,jci1,jci2,ici1,ici2)
#endif
          call subgrid_collect(lms%tgrd,tgrd_io,jci1,jci2,ici1,ici2)
          call subgrid_collect(lms%tgbrd,tgbrd_io,jci1,jci2,ici1,ici2)
          call subgrid_collect(lms%tlef,tlef_io,jci1,jci2,ici1,ici2)
          call subgrid_collect(lms%sncv,sncv_io,jci1,jci2,ici1,ici2)
          call subgrid_collect(lms%sfice,sfice_io,jci1,jci2,ici1,ici2)
          call subgrid_collect(lms%snag,snag_io,jci1,jci2,ici1,ici2)
          call subgrid_collect(lms%emisv,emisv_io,jci1,jci2,ici1,ici2)
          call subgrid_collect(lms%um10,um10_io,jci1,jci2,ici1,ici2)
          call subgrid_collect(lms%swalb,swalb_io,jci1,jci2,ici1,ici2)
          call subgrid_collect(lms%lwalb,lwalb_io,jci1,jci2,ici1,ici2)
          call subgrid_collect(lms%swdiralb,swdiralb_io,jci1,jci2,ici1,ici2)
          call subgrid_collect(lms%swdifalb,swdifalb_io,jci1,jci2,ici1,ici2)
          call subgrid_collect(lms%lwdiralb,lwdiralb_io,jci1,jci2,ici1,ici2)
          call subgrid_collect(lms%lwdifalb,lwdifalb_io,jci1,jci2,ici1,ici2)
          call subgrid_collect(mdsub%ldmsk,ldmsk1_io,jci1,jci2,ici1,ici2)
          call grid_collect(solis,solis_io,jci1,jci2,ici1,ici2)
          call grid_collect(solvs,solvs_io,jci1,jci2,ici1,ici2)
          call grid_collect(solvsd,solvsd_io,jci1,jci2,ici1,ici2)
          call grid_collect(solvl,solvl_io,jci1,jci2,ici1,ici2)
          call grid_collect(solvld,solvld_io,jci1,jci2,ici1,ici2)
          call grid_collect(sabveg,sabveg_io,jci1,jci2,ici1,ici2)
          call grid_collect(flw,flw_io,jci1,jci2,ici1,ici2)
          call grid_collect(flwd,flwd_io,jci1,jci2,ici1,ici2)
          call grid_collect(fsw,fsw_io,jci1,jci2,ici1,ici2)
          call grid_collect(sinc,sinc_io,jci1,jci2,ici1,ici2)
          call grid_collect(mddom%ldmsk,ldmsk_io,jci1,jci2,ici1,ici2)
#ifndef CLM
          if ( lakemod == 1 ) then
            call subgrid_collect(lms%eta,eta_io,jci1,jci2,ici1,ici2)
            call subgrid_collect(lms%hi,hi_io,jci1,jci2,ici1,ici2)
            call subgrid_collect(lms%tlake,tlak_io,jci1,jci2,ici1,ici2,1,ndpmax)
          end if
#else
          if ( imask == 2 ) then
            call grid_collect(mddom%lndcat,lndcat_io,jci1,jci2,ici1,ici2)
          end if
#endif
          if ( idcsst == 1 ) then
            call subgrid_collect(lms%sst,sst_io,jci1,jci2,ici1,ici2)
            call subgrid_collect(lms%tskin,tskin_io,jci1,jci2,ici1,ici2)
            call subgrid_collect(lms%deltas,deltas_io,jci1,jci2,ici1,ici2)
            call subgrid_collect(lms%tdeltas,tdeltas_io,jci1,jci2,ici1,ici2)
          end if
          if ( idynamic == 1 ) then
            call grid_collect(dstor,dstor_io,jde1,jde2,ide1,ide2,1,nsplit)
            call grid_collect(hstor,hstor_io,jde1,jde2,ide1,ide2,1,nsplit)
          end if
          if ( ichem == 1 ) then
            call grid_collect(convpr,convpr_io,jci1,jci2,ici1,ici2,1,kz)
            call grid_collect(rainout,rainout_io,jci1,jci2,ici1,ici2,1,kz,1,ntr)
            call grid_collect(washout,washout_io,jci1,jci2,ici1,ici2,1,kz,1,ntr)
            call grid_collect(remdrd,remdrd_io,jci1,jci2,ici1,ici2,1,ntr)
            if ( igaschem == 1 .and. ichsolver > 0 ) then
              call grid_collect(chemall,chemall_io,jci1,jci2,ici1,ici2, &
                                1,kz,1,totsp)
              call grid_collect(taucldsp,taucldsp_io,jci1,jci2,ici1,ici2, &
                                0,kz,1,nspi)
            end if
            call grid_collect(ssw2da,ssw2da_io,jci1,jci2,ici1,ici2)
#ifdef CLM45
            call grid_collect(dustflx_clm,duflux_io,jci1,jci2,ici1,ici2,1,4)
            call grid_collect(voc_em_clm,voflux_io,jci1,jci2,ici1,ici2,1,ntr)
#else
            call grid_collect(sdelt,sdelt_io,jci1,jci2,ici1,ici2)
            call grid_collect(sdelq,sdelq_io,jci1,jci2,ici1,ici2)
            call grid_collect(svegfrac2d,svegfrac2d_io,jci1,jci2,ici1,ici2)
#endif
            call grid_collect(sfracv2d,sfracv2d_io,jci1,jci2,ici1,ici2)
            call grid_collect(sfracb2d,sfracb2d_io,jci1,jci2,ici1,ici2)
            call grid_collect(sfracs2d,sfracs2d_io,jci1,jci2,ici1,ici2)
          end if
          if ( islab_ocean == 1 .and. do_restore_sst ) then
            call grid_collect(qflux_restore_sst,qflux_restore_sst_io, &
                              jci1,jci2,ici1,ici2,1,12)
          end if
        end if

        call write_savefile(rcmtimer%idate)

      end if
    end if

    if ( associated(alarm_out_nwf) ) then
      if ( lnewf ) then
        call newoutfiles(rcmtimer%idate)
        call checktime(myid,trim(dirout)//pthsep//trim(prestr)// &
                       trim(domname)//'.'//tochar10(lastout))
        lastout = rcmtimer%idate
      end if
    else
      if ( lfdomonth(rcmtimer%idate) .and. lmidnight(rcmtimer%idate) ) then
        if ( .not. lstartup .and. rcmtimer%idate /= idate2 ) then
          call newoutfiles(rcmtimer%idate)
          call checktime(myid,trim(dirout)//pthsep//trim(prestr)// &
                         trim(domname)//'.'//tochar10(lastout))
          lastout = rcmtimer%idate
        end if
      end if
    end if

#ifdef DEBUG
    call time_end(subroutine_name,idindx)
#endif

    contains

#include <pfesat.inc>
#include <pfwsat.inc>

  end subroutine output

  subroutine vertint(f3,p3,f2,plev)
    implicit none
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: f3
    real(rkx) , pointer , dimension(:,:,:) , intent(in) :: p3
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: f2
    real(rkx) , intent(in) :: plev
    integer(ik4) :: i , j , ik
    real(rkx) , dimension(kz) :: f1 , p1
    real(rkx) :: blw , tlw , dp

    do i = ici1 , ici2
      do j = jci1 , jci2
        f1 = f3(j,i,:)
        p1 = p3(j,i,:)
        ik = findlev()
        if ( ik < 1 ) then
          ! higher than top
          f2(j,i) = f3(j,i,1)
        else if ( ik > kz-1 ) then
          ! lower than bottom
          f2(j,i) = f3(j,i,kz)
        else
          ! in between two levels
          dp = p3(j,i,ik+1) - p3(j,i,ik)
          blw = (plev - p3(j,i,ik)) / dp
          tlw = d_one - blw
          f2(j,i) = (f3(j,i,ik+1)*blw+f3(j,i,ik)*tlw)
        end if
      end do
    end do

    contains

    integer(ik4) function findlev() result(kk)
      implicit none
      integer(ik4) :: k
      kk = 0
      if ( plev >= p1(1) ) then
        do k = 1 , kz
          if ( plev > p1(k) ) then
            kk = k
          end if
        end do
      end if
    end function findlev

  end subroutine vertint
  !
  ! Change U and V from map values (X,Y) to true (N,E)
  !
  subroutine uvrot2d(u,v)
    implicit none
    real(rkx) , pointer , dimension(:,:) , intent(inout) :: u , v
    if ( .not. rotinit ) then
      call alpharot_compute
    end if
    call pj%wind2_antirotate(u,v)
  end subroutine uvrot2d

  !
  ! Change U and V from map values (X,Y) to true (N,E)
  !
  subroutine uvrot3d(u,v)
    implicit none
    real(rkx) , pointer , dimension(:,:,:) , intent(inout) :: u , v
    if ( .not. rotinit ) then
      call alpharot_compute
    end if
    call pj%wind_antirotate(u,v)
  end subroutine uvrot3d

  subroutine alpharot_compute
    implicit none
    type(anyprojparams) :: pjpara
    if ( debug_level > 3 ) then
      if ( myid == italk ) then
        write(stdout,*) 'Computing rotation coefficients'
      end if
    end if
    pjpara%pcode = iproj
    pjpara%ds = dx
    pjpara%clat = clat
    pjpara%clon = clon
    pjpara%plat = plat
    pjpara%plon = plon
    pjpara%trlat1 = truelatl
    pjpara%trlat2 = truelath
    pjpara%nlon = jx
    pjpara%nlat = iy
    pjpara%rotparam = .true.
    pjpara%staggerx = .false.
    pjpara%staggery = .false.
    call pj%initialize(pjpara)
    if ( debug_level > 3 ) then
      if ( myid == italk ) then
        write(stdout,*) 'Done'
      end if
    end if
    rotinit = .true.
  end subroutine alpharot_compute

end module mod_output

! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

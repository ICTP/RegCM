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
 
module mod_output

  use mod_runparams
  use mod_header
  use mod_mpmessage
  use mod_mppparam
  use mod_service
  use mod_atm_interface
  use mod_che_interface
  use mod_che_output
  use mod_lm_interface
  use mod_rad_interface
  use mod_cu_interface
  use mod_pbl_interface
  use mod_ncio
  use mod_ncout
  use mod_bdycod
  use mod_precip
  use mod_split
  use mod_savefile
  use mod_mppio
#ifdef CLM
  use mod_clm
#endif

  private

  public :: output

  contains

  subroutine output
    implicit none
    logical :: ldoatm , ldosrf , ldorad , ldoche
    logical :: ldosav , ldolak , ldosub , ldotmp
    logical :: lstartup
    integer(ik4) :: i , j , k , jp1 , ip1
#ifdef DEBUG
    character(len=dbgslen) :: subroutine_name = 'output'
    integer(ik4) , save :: idindx = 0
    call time_begin(subroutine_name,idindx)
#endif
!
    lstartup = .false.
    if ( ktau == 0 .or. doing_restart ) then
      !
      ! Set up static variables (first time in)
      !
      if ( associated(xlon_out) ) then
        xlon_out = mddom%xlon(jci1:jci2,ici1:ici2)
        xlat_out = mddom%xlat(jci1:jci2,ici1:ici2)
        mask_out = mddom%mask(jci1:jci2,ici1:ici2)
        topo_out = mddom%ht(jci1:jci2,ici1:ici2)
      end if
      if ( associated(sub_xlon_out) ) then
        call reorder_subgrid(xlon1,sub_xlon_out)
        call reorder_subgrid(xlat1,sub_xlat_out)
        call reorder_subgrid(mask1,sub_mask_out)
        call reorder_subgrid(ht1,sub_topo_out)
      end if
      call newoutfiles(idatex)

      ! This must be removed

      if ( ifchem .and. myid == iocpu ) then
        call prepare_chem_out(idatex,ifrest)
      end if

      lstartup = .true.
      if ( doing_restart ) then
        doing_restart = .false.
#ifdef DEBUG
        call time_end(subroutine_name,idindx) 
#endif
        return
      end if
    end if
!
    ldoatm = .false.
    ldosrf = .false.
    ldolak = .false.
    ldosub = .false.
    ldorad = .false.
    ldoche = .false.
    ldosav = .false.
    ldotmp = .false.

    if ( ktau > 0 ) then
      if ( mod(ktau,ksav) == 0 ) then
        ldotmp = .true.
      end if
      if ( ( idatex == idate2 .or. &
           (lfdomonth(idatex) .and. lmidnight(idatex))) ) then
        ldosav = .true.
        ldotmp = .false.
      end if
      if ( mod(ktau,katm) == 0 ) then
        ldoatm = .true.
      end if
      if ( mod(ktau,ksrf) == 0 ) then
        ldosrf = .true.
      end if
      if ( mod(ktau,klak) == 0 ) then
        ldolak= .true.
      end if
      if ( mod(ktau,ksub) == 0 ) then
        ldosub= .true.
      end if
      if ( mod(ktau,krad) == 0 ) then
        ldorad = .true.
      end if
      if ( mod(ktau,kche) == 0 ) then
        ldoche = .true.
      end if
    end if
!
    if ( ktau == 0 ) then
      if ( debug_level > 2 ) then
        ldoatm = .true.
      end if
    end if
!
    if ( atm_stream > 0 ) then
      if ( ldoatm ) then
        ps_out = sfs%psa(jci1:jci2,ici1:ici2)
        if ( associated(atm_t_out) ) then
          do k = 1 , kz
            atm_t_out(:,:,k) = atm1%t(jci1:jci2,ici1:ici2,k)/ps_out
          end do
        end if
        if ( associated(atm_u_out) .and. associated(atm_v_out) ) then
          do k = 1 , kz
            do i = ici1 , ici2
              ip1 = i+1
              if ( ip1 > iy ) ip1 = 1
              do j = jci1 , jci2
                jp1 = j+1
                if ( jp1 > jx ) jp1 = 1
                atm_u_out(j,i,k) = d_rfour*(atm1%u(j,i,k)+atm1%u(jp1,i,k) + &
                                 atm1%u(j,ip1,k)+atm1%u(jp1,ip1,k))/ps_out(j,i)
                atm_v_out(j,i,k) = d_rfour*(atm1%v(j,i,k)+atm1%v(jp1,i,k) + &
                                 atm1%v(j,ip1,k)+atm1%v(jp1,ip1,k))/ps_out(j,i)
              end do
            end do
          end do
        end if
        if ( associated(atm_omega_out) ) &
          atm_omega_out = omega(jci1:jci2,ici1:ici2,:)*d_10
        if ( associated(atm_qv_out) ) then
          do k = 1 , kz
            atm_qv_out(:,:,k) = atm1%qx(jci1:jci2,ici1:ici2,k,iqv)/ps_out
          end do
          atm_qv_out = atm_qv_out/(d_one+atm_qv_out)
        end if
        if ( associated(atm_qc_out) ) then
          do k = 1 , kz
            atm_qc_out(:,:,k) = atm1%qx(jci1:jci2,ici1:ici2,k,iqc)/ps_out
          end do
        end if
        if ( associated(atm_tke_out) ) &
          atm_tke_out = atm1%tke(jci1:jci2,ici1:ici2,1:kz)
        if ( associated(atm_kth_out) ) &
          atm_kth_out = uwstateb%kth(jci1:jci2,ici1:ici2,1:kz)
        if ( associated(atm_kzm_out) ) &
          atm_kzm_out = uwstateb%kzm(jci1:jci2,ici1:ici2,1:kz)

        ps_out = d_10*(ps_out+ptop)

        if ( associated(atm_tpr_out) ) &
          atm_tpr_out = (sfs%rainc+sfs%rainnc)*rsecpd
        if ( associated(atm_tgb_out) ) &
          atm_tgb_out = atm_tgb_out * rsrf_in_atm
        if ( associated(atm_tsw_out) ) &
          atm_tsw_out = atm_tsw_out * rsrf_in_atm

        call write_record_output_stream(atm_stream,idatex)
        if ( myid == italk ) &
          write(stdout,*) 'ATM variables written at ' , tochar(idatex)

        atm_tgb_out = d_zero
        atm_tsw_out = d_zero
        sfs%rainc   = d_zero
        sfs%rainnc  = d_zero

      end if
    end if

    if ( srf_stream > 0 ) then
      if ( ldosrf ) then

        ps_out = d_10*(sfs%psa(jci1:jci2,ici1:ici2)+ptop)
        ! Averaged values
        if ( associated(srf_tpr_out) ) &
          srf_tpr_out = srf_tpr_out*rnsrf_for_srffrq
        if ( associated(srf_prcv_out) ) &
          srf_prcv_out = srf_prcv_out*rnsrf_for_srffrq
        if ( associated(srf_zpbl_out) ) &
          srf_zpbl_out = srf_zpbl_out*rnsrf_for_srffrq
        if ( associated(srf_evp_out) ) &
          srf_evp_out = srf_evp_out*rnsrf_for_srffrq
        if ( associated(srf_scv_out) ) then
          where ( ldmsk > 0 )
            srf_scv_out = srf_scv_out*rnsrf_for_srffrq
          elsewhere
            srf_scv_out = dmissval
          end where
        end if
        if ( associated(srf_runoff_out) ) then
          where ( ldmsk > 0 )
            srf_runoff_out(:,:,1) = srf_runoff_out(:,:,1)*rnsrf_for_srffrq
            srf_runoff_out(:,:,2) = srf_runoff_out(:,:,2)*rnsrf_for_srffrq - &
              srf_runoff_out(:,:,1)
          else where
            srf_runoff_out(:,:,1) = dmissval
            srf_runoff_out(:,:,2) = dmissval
          end where
        end if
        if ( associated(srf_sena_out) ) &
          srf_sena_out = -srf_sena_out*rnsrf_for_srffrq
        if ( associated(srf_flw_out) ) &
          srf_flw_out = srf_flw_out*rnsrf_for_srffrq
        if ( associated(srf_fsw_out) ) &
          srf_fsw_out = srf_fsw_out*rnsrf_for_srffrq
        if ( associated(srf_fld_out) ) &
          srf_fld_out = srf_fld_out*rnsrf_for_srffrq
        if ( associated(srf_sina_out) ) &
          srf_sina_out = srf_sina_out*rnsrf_for_srffrq

        call write_record_output_stream(srf_stream,idatex)
        if ( myid == italk ) &
          write(stdout,*) 'SRF variables written at ' , tochar(idatex)

        if ( associated(srf_tpr_out) ) srf_tpr_out = d_zero
        if ( associated(srf_prcv_out) ) srf_prcv_out = d_zero
        if ( associated(srf_zpbl_out) ) srf_zpbl_out = d_zero
        if ( associated(srf_evp_out) ) srf_evp_out = d_zero
        if ( associated(srf_scv_out) ) srf_scv_out = d_zero
        if ( associated(srf_runoff_out) ) srf_runoff_out = d_zero
        if ( associated(srf_sena_out) ) srf_sena_out = d_zero
        if ( associated(srf_flw_out) ) srf_flw_out = d_zero
        if ( associated(srf_fsw_out) ) srf_fsw_out = d_zero
        if ( associated(srf_fld_out) ) srf_fld_out = d_zero
        if ( associated(srf_sina_out) ) srf_sina_out = d_zero
        if ( associated(srf_sund_out) ) srf_sund_out = d_zero
      end if
    end if

    if ( sub_stream > 0 ) then
      if ( ldosub ) then

        sub_ps_out = sub_ps_out*rnsrf_for_subfrq

        if ( associated(sub_evp_out) ) &
          sub_evp_out = sub_evp_out*rnsrf_for_subfrq
        if ( associated(sub_scv_out) ) then
          where ( sub_scv_out < dmissval )
            sub_scv_out = sub_scv_out*rnsrf_for_subfrq
          end where
        end if
        if ( associated(sub_sena_out) ) &
          sub_sena_out = sub_sena_out*rnsrf_for_subfrq
        if ( associated(sub_runoff_out) ) then
          where ( sub_runoff_out(:,:,1) < dmissval )
            sub_runoff_out(:,:,2) = sub_runoff_out(:,:,2)-sub_runoff_out(:,:,1)
            sub_runoff_out(:,:,1) = sub_runoff_out(:,:,1)*rnsrf_for_subfrq
            sub_runoff_out(:,:,2) = sub_runoff_out(:,:,2)*rnsrf_for_subfrq
          end where
        end if

        call write_record_output_stream(sub_stream,idatex)
        if ( myid == italk ) &
          write(stdout,*) 'SUB variables written at ' , tochar(idatex)

        if ( associated(sub_evp_out) ) sub_evp_out = d_zero
        if ( associated(sub_scv_out) ) sub_scv_out = d_zero
        if ( associated(sub_sena_out) ) sub_sena_out = d_zero
        if ( associated(sub_runoff_out) ) sub_runoff_out = d_zero
      end if
    end if

    if ( lak_stream > 0 ) then
      if ( ldolak ) then

        ps_out = d_10*(sfs%psa(jci1:jci2,ici1:ici2)+ptop)
        if ( associated(lak_tpr_out) ) &
          lak_tpr_out = lak_tpr_out*rnsrf_for_lakfrq
        if ( associated(lak_scv_out) ) then
          where ( ldmsk > 0 )
            lak_scv_out = lak_scv_out*rnsrf_for_lakfrq
          else where
            lak_scv_out = dmissval
          end where
        end if
        if ( associated(lak_sena_out) ) &
          lak_sena_out = lak_sena_out*rnsrf_for_lakfrq
        if ( associated(lak_flw_out) ) &
          lak_flw_out = lak_flw_out*rnsrf_for_lakfrq
        if ( associated(lak_fsw_out) ) &
          lak_fsw_out = lak_fsw_out*rnsrf_for_lakfrq
        if ( associated(lak_fld_out) ) &
          lak_fld_out = lak_fld_out*rnsrf_for_lakfrq
        if ( associated(lak_sina_out) ) &
          lak_sina_out = lak_sina_out*rnsrf_for_lakfrq
        if ( associated(lak_evp_out) ) &
          lak_evp_out = lak_evp_out*rnsrf_for_lakfrq
        if ( associated(lak_aveice_out) ) then
          where ( lak_aveice_out < dmissval )
            lak_aveice_out = lak_aveice_out*rnsrf_for_lakfrq
          end where
        end if

        call write_record_output_stream(lak_stream,idatex)
        if ( myid == italk ) &
          write(stdout,*) 'LAK variables written at ' , tochar(idatex)

        if ( associated(lak_tpr_out) )    lak_tpr_out = d_zero
        if ( associated(lak_scv_out) )    lak_scv_out = d_zero
        if ( associated(lak_sena_out) )   lak_sena_out = d_zero
        if ( associated(lak_flw_out) )    lak_flw_out = d_zero
        if ( associated(lak_fsw_out) )    lak_fsw_out = d_zero
        if ( associated(lak_fld_out) )    lak_fld_out = d_zero
        if ( associated(lak_sina_out) )   lak_sina_out = d_zero
        if ( associated(lak_evp_out) )    lak_evp_out = d_zero
        if ( associated(lak_aveice_out) ) then
          where ( lak_aveice_out < dmissval )
            lak_aveice_out = d_zero
          end where
        end if

      end if
    end if

    if ( opt_stream > 0 ) then
      if ( ldoche ) then
        ps_out = d_10*(sfs%psa(jci1:jci2,ici1:ici2)+ptop)
        if ( iaerosol == 1 ) then
          if ( associated(opt_acstoarf_out) ) &
            opt_acstoarf_out = opt_acstoarf_out * rnrad_for_chem
          if ( associated(opt_acstsrrf_out) ) &
            opt_acstsrrf_out = opt_acstsrrf_out * rnrad_for_chem
          if ( associated(opt_acstalrf_out) ) &
            opt_acstalrf_out = opt_acstalrf_out * rnrad_for_chem
          if ( associated(opt_acssrlrf_out) ) &
            opt_acssrlrf_out = opt_acssrlrf_out * rnrad_for_chem
          if ( associated(opt_aod_out) ) &
            opt_aod_out = opt_aod_out * rnrad_for_chem
          call write_record_output_stream(opt_stream,idatex)
          if ( myid == italk ) &
            write(stdout,*) 'OPT variables written at ' , tochar(idatex)
          if ( associated(opt_acstoarf_out) ) opt_acstoarf_out = d_zero
          if ( associated(opt_acstsrrf_out) ) opt_acstsrrf_out = d_zero
          if ( associated(opt_acstalrf_out) ) opt_acstalrf_out = d_zero
          if ( associated(opt_acssrlrf_out) ) opt_acssrlrf_out = d_zero
          if ( associated(opt_aod_out) ) opt_aod_out = d_zero

          ! For now keep this

          call output_chem(idatex)

        end if
      end if
    end if

    if ( sts_stream > 0 ) then
      if ( mod(ktau+kstsoff,ksts) == 0 .and. ktau > kstsoff+2 ) then

        ps_out = d_10*(sfs%psa(jci1:jci2,ici1:ici2)+ptop)
        if ( associated(sts_pcpavg_out) ) &
          sts_pcpavg_out = sts_pcpavg_out*rnsrf_for_day
        if ( associated(sts_t2avg_out) ) &
          sts_t2avg_out = sts_t2avg_out*rnsrf_for_day

        call write_record_output_stream(sts_stream,idatex)
        if ( myid == italk ) &
          write(stdout,*) 'STS variables written at ' , tochar(idatex)

        if ( associated(sts_pcpavg_out) ) sts_pcpavg_out = d_zero
        if ( associated(sts_t2avg_out) )  sts_t2avg_out  = d_zero
        if ( associated(sts_tgmax_out) )  sts_tgmax_out  = -1.D30
        if ( associated(sts_tgmin_out) )  sts_tgmin_out  =  1.D30
        if ( associated(sts_t2max_out) )  sts_t2max_out  = -1.D30
        if ( associated(sts_t2min_out) )  sts_t2min_out  =  1.D30
        if ( associated(sts_w10max_out) ) sts_w10max_out = -1.D30
        if ( associated(sts_psmin_out) )  sts_psmin_out  =  1.D30
        if ( associated(sts_pcpmax_out) ) sts_pcpmax_out = -1.D30
        if ( associated(sts_sund_out) )   sts_sund_out   = d_zero
      end if
    end if

    if ( rad_stream > 0 ) then
      if ( ldorad ) then
        ps_out = d_10*(sfs%psa(jci1:jci2,ici1:ici2)+ptop)
        if ( associated(rad_frsa_out) ) &
          rad_frsa_out = rad_frsa_out * rnrad_for_radfrq
        if ( associated(rad_frla_out) ) &
          rad_frla_out = rad_frla_out * rnrad_for_radfrq
        if ( associated(rad_clrst_out) ) &
          rad_clrst_out = rad_clrst_out * rnrad_for_radfrq
        if ( associated(rad_clrss_out) ) &
          rad_clrss_out = rad_clrss_out * rnrad_for_radfrq
        if ( associated(rad_clrls_out) ) &
          rad_clrls_out = rad_clrls_out * rnrad_for_radfrq
        if ( associated(rad_clrlt_out) ) &
          rad_clrlt_out = rad_clrlt_out * rnrad_for_radfrq
        if ( associated(rad_sabtp_out) ) &
          rad_sabtp_out = rad_sabtp_out * rnrad_for_radfrq
        if ( associated(rad_firtp_out) ) &
          rad_firtp_out = rad_firtp_out * rnrad_for_radfrq
        call write_record_output_stream(rad_stream,idatex)
        if ( myid == italk ) &
          write(stdout,*) 'RAD variables written at ' , tochar(idatex)
        if ( associated(rad_frsa_out) ) rad_frsa_out = d_zero        
        if ( associated(rad_frla_out) ) rad_frla_out = d_zero        
        if ( associated(rad_clrst_out) ) rad_clrst_out = d_zero        
        if ( associated(rad_clrss_out) ) rad_clrss_out = d_zero        
        if ( associated(rad_clrls_out) ) rad_clrls_out = d_zero        
        if ( associated(rad_clrlt_out) ) rad_clrlt_out = d_zero        
        if ( associated(rad_sabtp_out) ) rad_sabtp_out = d_zero        
        if ( associated(rad_firtp_out) ) rad_firtp_out = d_zero        
      end if
    end if

    if ( ifsave ) then
      if ( ldosav .or. ldotmp ) then
        call grid_collect(atm1%u,atm1_io%u,jde1,jde2,ide1,ide2,1,kz)
        call grid_collect(atm1%v,atm1_io%v,jde1,jde2,ide1,ide2,1,kz)
        call grid_collect(atm1%t,atm1_io%t,jce1,jce2,ice1,ice2,1,kz)
        call grid_collect(atm1%qx,atm1_io%qx,jce1,jce2,ice1,ice2,1,kz,1,nqx)

        call grid_collect(atm2%u,atm2_io%u,jde1,jde2,ide1,ide2,1,kz)
        call grid_collect(atm2%v,atm2_io%v,jde1,jde2,ide1,ide2,1,kz)
        call grid_collect(atm2%t,atm2_io%t,jce1,jce2,ice1,ice2,1,kz)
        call grid_collect(atm2%qx,atm2_io%qx,jce1,jce2,ice1,ice2,1,kz,1,nqx)

        if ( ibltyp == 2 .or. ibltyp == 99 ) then
          call grid_collect(atm1%tke,atm1_io%tke,jce1,jce2,ice1,ice2,1,kzp1)
          call grid_collect(atm2%tke,atm2_io%tke,jce1,jce2,ice1,ice2,1,kzp1)
          call grid_collect(kpbl,kpbl_io,jci1,jci2,ici1,ici2)
        end if

        call grid_collect(sfs%psa,sfs_io%psa,jce1,jce2,ice1,ice2)
        call grid_collect(sfs%psb,sfs_io%psb,jce1,jce2,ice1,ice2)
        call grid_collect(sfs%tga,sfs_io%tga,jce1,jce2,ice1,ice2)
        call grid_collect(sfs%tgb,sfs_io%tgb,jce1,jce2,ice1,ice2)

        call grid_collect(sfs%hfx,sfs_io%hfx,jci1,jci2,ici1,ici2)
        call grid_collect(sfs%qfx,sfs_io%qfx,jci1,jci2,ici1,ici2)
        call grid_collect(sfs%rainc,sfs_io%rainc,jci1,jci2,ici1,ici2)
        call grid_collect(sfs%rainnc,sfs_io%rainnc,jci1,jci2,ici1,ici2)
        call grid_collect(sfs%tgbb,sfs_io%tgbb,jci1,jci2,ici1,ici2)
        call grid_collect(sfs%uvdrag,sfs_io%uvdrag,jci1,jci2,ici1,ici2)

        if ( ipptls == 1 ) then
          call grid_collect(fcc,fcc_io,jci1,jci2,ici1,ici2,1,kz)
        end if
        call grid_collect(heatrt,heatrt_io,jci1,jci2,ici1,ici2,1,kz)
        call grid_collect(o3prof,o3prof_io,jci1,jci2,ici1,ici2,1,kzp1)

        if ( iocnflx == 2 ) then
          call grid_collect(zpbl,zpbl_io,jci1,jci2,ici1,ici2)
        end if
        if ( icup == 1 ) then
          call grid_collect(rsheat,rsheat_io,jci1,jci2,ici1,ici2,1,kz)
          call grid_collect(rswat,rswat_io,jci1,jci2,ici1,ici2,1,kz)
        end if
        if ( icup == 3 ) then
          call grid_collect(tbase,tbase_io,jci1,jci2,ici1,ici2,1,kz)
          call grid_collect(cldefi,cldefi_io,jci1,jci2,ici1,ici2)
        end if
        if ( icup==4 .or. icup==99 .or. icup==98 ) then
          call grid_collect(cbmf2d,cbmf2d_io,jci1,jci2,ici1,ici2)
        end if

        if ( irrtm == 0 ) then
          call grid_collect(gasabsnxt,gasabsnxt_io,jci1,jci2,ici1,ici2,1,kz,1,4)
          call grid_collect(gasabstot,gasabstot_io, &
                            jci1,jci2,ici1,ici2,1,kzp1,1,kzp1)
          call grid_collect(gasemstot,gasemstot_io,jci1,jci2,ici1,ici2,1,kzp1)
        end if

        call subgrid_collect(tlef,tlef_io,jci1,jci2,ici1,ici2)
        call subgrid_collect(ssw,ssw_io,jci1,jci2,ici1,ici2)
        call subgrid_collect(rsw,rsw_io,jci1,jci2,ici1,ici2)
        call subgrid_collect(tgrd,tgrd_io,jci1,jci2,ici1,ici2)
        call subgrid_collect(tgbrd,tgbrd_io,jci1,jci2,ici1,ici2)
        call subgrid_collect(sncv,sncv_io,jci1,jci2,ici1,ici2)
        call subgrid_collect(gwet,gwet_io,jci1,jci2,ici1,ici2)
        call subgrid_collect(snag,snag_io,jci1,jci2,ici1,ici2)
        call subgrid_collect(sfice,sfice_io,jci1,jci2,ici1,ici2)
        call subgrid_collect(ldew,ldew_io,jci1,jci2,ici1,ici2)
        call subgrid_collect(taf,taf_io,jci1,jci2,ici1,ici2)
        call subgrid_collect(tsw,tsw_io,jci1,jci2,ici1,ici2)
        call subgrid_collect(emiss,emiss_io,jci1,jci2,ici1,ici2)
        call subgrid_collect(ldmsk1,ldmsk1_io,jci1,jci2,ici1,ici2)

        call grid_collect(solis,solis_io,jci1,jci2,ici1,ici2)
        call grid_collect(solvd,solvd_io,jci1,jci2,ici1,ici2)
        call grid_collect(solvs,solvs_io,jci1,jci2,ici1,ici2)
        call grid_collect(sabveg,sabveg_io,jci1,jci2,ici1,ici2)
        call grid_collect(flw,flw_io,jci1,jci2,ici1,ici2)
        call grid_collect(flwd,flwd_io,jci1,jci2,ici1,ici2)
        call grid_collect(fsw,fsw_io,jci1,jci2,ici1,ici2)
        call grid_collect(sinc,sinc_io,jci1,jci2,ici1,ici2)
        call grid_collect(ldmsk,ldmsk_io,jci1,jci2,ici1,ici2)

#ifndef CLM
        if ( lakemod == 1 ) then
          call subgrid_collect(eta,eta_io,jci1,jci2,ici1,ici2)
          call subgrid_collect(hi,hi_io,jci1,jci2,ici1,ici2)
          call subgrid_collect(aveice,aveice_io,jci1,jci2,ici1,ici2)
          call subgrid_collect(hsnow,hsnow_io,jci1,jci2,ici1,ici2)
          call subgrid_collect(tlak,tlak_io,jci1,jci2,ici1,ici2,1,ndpmax)
        end if
#else
        call grid_collect(sols2d,sols2d_io,jci1,jci2,ici1,ici2)
        call grid_collect(soll2d,soll2d_io,jci1,jci2,ici1,ici2)
        call grid_collect(solsd2d,solsd2d_io,jci1,jci2,ici1,ici2)
        call grid_collect(solld2d,solld2d_io,jci1,jci2,ici1,ici2)
        call grid_collect(aldirs2d,aldirs2d_io,jci1,jci2,ici1,ici2)
        call grid_collect(aldirl2d,aldirl2d_io,jci1,jci2,ici1,ici2)
        call grid_collect(aldifs2d,aldifs2d_io,jci1,jci2,ici1,ici2)
        call grid_collect(aldifl2d,aldifl2d_io,jci1,jci2,ici1,ici2)
        call grid_collect(lndcat2d,lndcat2d_io,jci1,jci2,ici1,ici2)
#endif
        if ( idcsst == 1 ) then
          call grid_collect(dtskin,dtskin_io,jci1,jci2,ici1,ici2)
          call grid_collect(deltas,deltas_io,jci1,jci2,ici1,ici2)
          call grid_collect(tdeltas,tdeltas_io,jci1,jci2,ici1,ici2)
        end if

        call grid_collect(dstor,dstor_io,jde1,jde2,ide1,ide2,1,nsplit)
        call grid_collect(hstor,hstor_io,jde1,jde2,ide1,ide2,1,nsplit)

        if ( ichem == 1 ) then
          call grid_collect(chia,chia_io,jce1,jce2,ice1,ice2,1,kz,1,ntr)
          call grid_collect(chib,chib_io,jce1,jce2,ice1,ice2,1,kz,1,ntr)
          call grid_collect(remlsc,remlsc_io,jce1,jce2,ice1,ice2,1,kz,1,ntr)
          call grid_collect(remcvc,remcvc_io,jce1,jce2,ice1,ice2,1,kz,1,ntr)
          call grid_collect(remdrd,remdrd_io,jce1,jce2,ice1,ice2,1,ntr)

          call grid_collect(ssw2da,ssw2da_io,jci1,jci2,ici1,ici2)
          call grid_collect(sdeltk2d,sdeltk2d_io,jci1,jci2,ici1,ici2)
          call grid_collect(sdelqk2d,sdelqk2d_io,jci1,jci2,ici1,ici2)
          call grid_collect(sfracv2d,sfracv2d_io,jci1,jci2,ici1,ici2)
          call grid_collect(sfracb2d,sfracb2d_io,jci1,jci2,ici1,ici2)
          call grid_collect(sfracs2d,sfracs2d_io,jci1,jci2,ici1,ici2)
          call grid_collect(svegfrac2d,svegfrac2d_io,jci1,jci2,ici1,ici2)
        end if

        if ( ldosav ) then
          call write_savefile(idatex, .false.)
        else
          call write_savefile(idatex, .true.)
        end if
      end if
    end if

    if ( lfdomonth(idatex) .and. lmidnight(idatex) ) then
      if ( .not. lstartup .and. idatex /= idate2 ) then
        call newoutfiles(idatex)

        ! This must be removed
        if ( ifchem .and. myid == iocpu ) then
          call prepare_chem_out(idatex,ifrest)
        end if

        call checktime(myid)
      end if
    end if
#ifdef DEBUG
    call time_end(subroutine_name,idindx) 
#endif
  end subroutine output
!
end module mod_output

!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::
!
!    This file is part of ICTP RegCM.
!
!    Use of this source code is governed by an MIT-style license that can
!    be found in the LICENSE file or at
!
!         https://opensource.org/licenses/MIT.
!
!    ICTP RegCM is distributed in the hope that it will be useful,
!    but WITHOUT ANY WARRANTY; without even the implied warranty of
!    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.
!
!::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::

module mod_che_interface

  use mod_runparams, only : iaerosol, iindirect, pc_total
  use mod_che_common, only : allocate_mod_che_common, assignpnt, bndp0, &
    bndp1, cadvhdiag, cadvvdiag, cbdydiag, cccn, ccldfra, cconvdiag, cconvpr, &
    cdifhdiag, cdlat, cdlon, cdrydepflx, cdzq, cfcc, cfmz, chem_config, &
    chemall, chemt, chia, chib, chib3d, chifxuw, chiten, cht, clndcat, &
    convcldfra, cpb3d, cpf3d, cps0, cps2d, cpsb, cptrop, cqxb3d, cra, crembc, &
    cremrat, crhb3d, crho2d, crhob3d, csdelqk2d, csdeltk2d, csfracb2d, &
    csfracs2d, csfracv2d, csol2d, cssw2da, ctaucld, ctb3d, ctbldiag, ctg, &
    ctga, cubx3d, custar, cvbx3d, cveg2d, cvegfrac, cw10m, cwetdepflx, &
    cxlai2d, cxlat, cza, czen, czo, czq, drydepv, idust, idynamic, kcumbot, &
    kcumtop, psbb0, psbb1, rainout, remdrd, tvirt0, tvirt1, washout
#if defined(CLM45) || defined(ECLM)
  use mod_che_common, only : cddepv_clm, cdustflx_clm, csw_vol, ctsoi, &
    cvoc_em_clm
#endif
  use mod_che_cumtran, only : cumtran, init_cumtran
  use mod_che_dust, only : allocate_mod_che_dust
  use mod_che_indices, only : totsp
  use mod_che_mppio, only : allocate_mod_che_mppio, ccn_io, chia_io, chib_io, &
    chifxuw_io, cra_io, drydepflx_io, drydepv_io, duflux_io, rainout_io, &
    remdrd_io, sfracb2d_io, sfracs2d_io, ssw2da_io, svegfrac2d_io, &
    sxlai2d_io, trac_io, voflux_io, washout_io, wetdepflx_io
  use mod_che_ncio, only : close_chbc
  use mod_che_bdyco, only : allocate_mod_che_bdyco, cba, chem_bdyin, &
    chem_bdyval, morelax_chiten, nudge_chiten, setup_che_bdycon
  use mod_che_tend, only : tractend2
  use mod_che_start, only : start_chem
  use mod_che_bionit, only : allocate_mod_che_bionit

  implicit none

  private

  public :: start_chem
  public :: init_chem
  public :: cumtran
  public :: tractend2
  public :: nudge_chiten
  public :: morelax_chiten
  public :: chem_config
  public :: setup_che_bdycon
  public :: close_chbc

  public :: allocate_mod_che_common
  public :: allocate_mod_che_mppio
  public :: allocate_mod_che_dust
  public :: allocate_mod_che_bdyco
  public :: allocate_mod_che_bionit

  public :: idust
  public :: totsp
  public :: trac_io
  public :: chia_io
  public :: chib_io
  public :: cra
  public :: convcldfra, cadvhdiag, cadvvdiag, cbdydiag, cconvdiag
  public :: cdifhdiag, ctbldiag
  public :: chem_bdyin, chem_bdyval
  public :: chemall
  public :: washout, washout_io
  public :: remdrd, remdrd_io
  public :: rainout, rainout_io
  public :: sfracb2d_io, sfracs2d_io, ssw2da_io
  public :: duflux_io, voflux_io, svegfrac2d_io
  public :: drydepv, chifxuw
  public :: sxlai2d_io, cra_io, ccn_io, drydepv_io, chifxuw_io
  public :: wetdepflx_io, drydepflx_io

  contains

  subroutine init_chem
    ! this routine define the pointer interface between the chem module and
    ! the rest of the model
    ! It also call startchem which is the chemistry initialisation routine

    use mod_atm_interface, only : aten, atm0, atm1, atm2, atms, ba_cr, ccn, &
      cldfra, convpr, coszrs, drydepflx, fcc, icumbot, icumtop, mddom, &
      mo_atm, nhbh0, nhbh1, ptrop, rembc, remrat, sdelq, sdelt, sfracb2d, &
      sfracs2d, sfracv2d, sfs, solis, ssw2da, svegfrac2d, sxlai2d, wetdepflx, &
      xpsb
#if defined(CLM45) || defined(ECLM)
    use mod_atm_interface, only : ddepv_clm, dustflx_clm, sw_vol, tsoi, &
      voc_em_clm
#endif
    use mod_rad_interface, only : taucldsp
    implicit none

    call assignpnt(icumtop,kcumtop)
    call assignpnt(icumbot,kcumbot)
    call assignpnt(atms%tb3d,ctb3d)
    call assignpnt(atms%qxb3d,cqxb3d)
    call assignpnt(atms%rhob3d,crhob3d)
    call assignpnt(atms%pb3d,cpb3d)
    call assignpnt(atms%pf3d,cpf3d)
    call assignpnt(atms%ps2d,cps2d)
    call assignpnt(atms%rhb3d,crhb3d)
    if ( idynamic == 2 ) then
      call assignpnt(atm0%ps,cps0)
      call assignpnt(nhbh0%ps,bndp0)
      call assignpnt(nhbh1%ps,bndp1)
      call assignpnt(nhbh0%tvirt,tvirt0)
      call assignpnt(nhbh1%tvirt,tvirt1)
    end if
    ! wind at cell center
    call assignpnt(atms%ubx3d,cubx3d)
    call assignpnt(atms%vbx3d,cvbx3d)
    call assignpnt(atms%chib3d,chib3d)
    call assignpnt(atms%rhox2d,crho2d)

    call assignpnt(mddom%lndcat,clndcat)
    call assignpnt(mddom%xlat,cxlat)
    call assignpnt(mddom%dlat,cdlat)
    call assignpnt(mddom%dlon,cdlon)

    if ( idynamic == 3 ) then
      call assignpnt(mo_atm%trac,chemt)
      call assignpnt(mo_atm%chiten,chiten)
      call assignpnt(mo_atm%fmz,cfmz)
      call assignpnt(nhbh0%ps,bndp0)
      call assignpnt(nhbh1%ps,bndp1)
      call assignpnt(nhbh0%tvirt,tvirt0)
      call assignpnt(nhbh1%tvirt,tvirt1)
    else
      call assignpnt(atm1%chi,chia)
      call assignpnt(atm2%chi,chib)
      call assignpnt(aten%chi,chiten,pc_total)
    end if

    call assignpnt(mddom%ht,cht)
    call assignpnt(mddom%iveg,cveg2d)
    call assignpnt(sfs%psb,cpsb)
    call assignpnt(xpsb%b0,psbb0)
    call assignpnt(xpsb%b1,psbb1)
    call assignpnt(sfs%tg,ctga)
    call assignpnt(sfs%tgbb,ctg)

    call assignpnt(sfs%ustar,custar)
    call assignpnt(sfs%w10m,cw10m)
    call assignpnt(sfs%zo,czo)
    call assignpnt(convpr,cconvpr)
    call assignpnt(fcc,cfcc)
    call assignpnt(cldfra,ccldfra)
    call assignpnt(rembc,crembc)
    call assignpnt(remrat,cremrat)
    call assignpnt(solis,csol2d)

    call assignpnt(svegfrac2d,cvegfrac)
    call assignpnt(sfracv2d,csfracv2d)
    call assignpnt(sfracb2d,csfracb2d)
    call assignpnt(sfracs2d,csfracs2d)
    call assignpnt(sxlai2d,cxlai2d)

    call assignpnt(sdelt,csdeltk2d)
    call assignpnt(sdelq,csdelqk2d)

    call assignpnt(atms%za,cza)
    call assignpnt(atms%zq,czq)
    call assignpnt(atms%dzq,cdzq)
    call assignpnt(coszrs,czen)
    call assignpnt(ssw2da,cssw2da)

    call assignpnt(taucldsp, ctaucld)
    call assignpnt(ptrop, cptrop)

    cba%havebound = ba_cr%havebound
    call assignpnt(ba_cr%bsouth, cba%bsouth)
    call assignpnt(ba_cr%bnorth, cba%bnorth)
    call assignpnt(ba_cr%beast, cba%beast)
    call assignpnt(ba_cr%bwest, cba%bwest)
    call assignpnt(ba_cr%ibnd, cba%ibnd)

    cba%ns = ba_cr%ns
    cba%nn = ba_cr%nn
    cba%ne = ba_cr%ne
    cba%nw = ba_cr%nw
    cba%nsp = ba_cr%nspx

    call assignpnt(wetdepflx,cwetdepflx)
    call assignpnt(drydepflx,cdrydepflx)

    if ( iindirect > 0 .and. iaerosol == 1 ) then
      call assignpnt(ccn,cccn)
    end if

#if defined(CLM45) || defined(ECLM)
    call assignpnt(voc_em_clm,cvoc_em_clm)
    call assignpnt(dustflx_clm,cdustflx_clm)
    call assignpnt(ddepv_clm,cddepv_clm)
    call assignpnt(sw_vol,csw_vol)
    call assignpnt(tsoi,ctsoi)
#endif

    call init_cumtran

  end subroutine init_chem

end module mod_che_interface
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

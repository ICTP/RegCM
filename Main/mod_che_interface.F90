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

  use mod_realkinds
  use mod_regcm_types
  use mod_runparams
  use mod_che_common
  use mod_che_cumtran
  use mod_che_dust
  use mod_che_indices
  use mod_che_mppio
  use mod_che_ncio
  use mod_che_param
  use mod_che_drydep
  use mod_che_bdyco
  use mod_che_emission
  use mod_che_carbonaer
  use mod_che_species
  use mod_che_tend
  use mod_che_start
  use mod_che_bionit
  use mod_che_linox

  implicit none

  private

  public :: start_chem
  public :: init_chem
  public :: cumtran
  public :: tractend2
  public :: nudge_chi
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
  public :: convcldfra, cadvhdiag, cadvvdiag, cbdydiag, cconvdiag
  public :: cdifhdiag, ctbldiag
  public :: chem_bdyin, chem_bdyval
  public :: chemall, chemall_io
  public :: washout, washout_io
  public :: remdrd, remdrd_io
  public :: rainout, rainout_io, convpr_io
  public :: sdelq_io, sdelt_io, sfracb2d_io, sfracs2d_io, ssw2da_io
  public :: duflux_io, voflux_io, sfracv2d_io, svegfrac2d_io, taucldsp_io

  contains

  subroutine init_chem
    ! this routine define the pointer interface between the chem module and
    ! the rest of the model
    ! It also call startchem which is the chemistry initialisation routine

    use mod_atm_interface
    use mod_rad_interface
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
    call assignpnt(sfs%ram1,cra)
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
    cba%nsp = ba_cr%nsp

    call assignpnt(wetdepflx,cwetdepflx)
    call assignpnt(drydepflx,cdrydepflx)

    if ( iindirect > 0 .and. iaerosol == 1 ) then
      call assignpnt(ccn,cccn)
    end if

#if (defined CLM45)
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

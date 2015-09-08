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

module mod_che_interface
!
  use mod_realkinds
  use mod_regcm_types
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

  public :: totsp
  public :: chi , chib3d , chic
  public :: chia , chia_io
  public :: chib , chib_io
  public :: chiten , chiten0
  public :: convcldfra , cadvhdiag , cadvvdiag , cbdydiag , cconvdiag
  public :: cdifhdiag , ctbldiag
  public :: chem_bdyin , chem_bdyval
  public :: chemall , chemall_io
  public :: washout , washout_io
  public :: remdrd , remdrd_io
  public :: rainout , rainout_io
  public :: sdelq_io , sdelt_io , sfracb2d_io , sfracs2d_io , ssw2da_io
  public :: sfracv2d_io , svegfrac2d_io , taucldsp_io

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
    call assignpnt(atms%ps2d,cps2d)
    ! wind at cell center
    call assignpnt(atms%ubx3d,cubx3d)
    call assignpnt(atms%vbx3d,cvbx3d)
    call assignpnt(atms%chib3d,chib3d)
    call assignpnt(atms%rhb3d,crhb3d)

    call assignpnt(mddom%lndcat,clndcat)
    call assignpnt(mddom%xlat,cxlat)
    call assignpnt(mddom%ht,cht)
    call assignpnt(mddom%iveg,cveg2d)
    call assignpnt(sfs%psb,cpsb)
    call assignpnt(xpsb%b0,psbb0)
    call assignpnt(xpsb%b1,psbb1)
    call assignpnt(sfs%tgb,ctg)
    call assignpnt(sfs%tga,ctga)
    call assignpnt(sfs%uvdrag,cuvdrag)
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
    call assignpnt(atms%dzq,cdzq)
    call assignpnt(coszrs,czen)
    call assignpnt(ssw2da,cssw2da)

    call assignpnt(taucldsp, ctaucld)

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

#if (defined CLM45)
    call assignpnt(voc_em,cvoc_em)
    call assignpnt(idusts,idust)
#endif
#if defined CLM
#if defined VOC
    call assignpnt(voc_em0,cvoc_em0)
    call assignpnt(voc_em1,cvoc_em1)
    call assignpnt(voc_em2,cvoc_em2)
#endif
    call assignpnt(dep_vels,cdep_vels)
#endif

  end subroutine init_chem
!
end module mod_che_interface
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

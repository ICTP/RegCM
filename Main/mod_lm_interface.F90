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

module mod_lm_interface
!
! Link surface and atmospheric models
!
  use mod_bats_common
  use mod_runparams
  use mod_memutil
  use mod_regcm_types
#ifdef CLM
  use mod_mtrxclm
  use mod_clm
  use mod_bats_mtrxbats
  use clm_varsur , only : landmask
#else
  use mod_bats_param
  use mod_bats_bndry
  use mod_bats_co2
  use mod_bats_drag
  use mod_bats_lake
  use mod_bats_leaftemp
  use mod_bats_mtrxbats
  use mod_bats_zengocn
#endif

  implicit none

  private

  public :: dtbat
  public :: dtlake
  public :: fdaysrf
  public :: cplmsk
  public :: sfice1
  public :: emiss1
  public :: gwet1
  public :: ldew1
  public :: rsw1
  public :: snag1
  public :: sncv1
  public :: ssw1
  public :: taf1
  public :: tgbrd1
  public :: tgrd1
  public :: tlef1
  public :: tsw1
  public :: sfracb2d
  public :: sfracs2d
  public :: sfracv2d
  public :: ssw2da
  public :: svegfrac2d
  public :: sst
  public :: dtskin
  public :: deltas
  public :: tdeltas

#ifndef CLM
  public :: var_aveice
  public :: var_eta
  public :: var_hi
  public :: var_hsnow
  public :: var_tlak
  public :: tlake
  public :: xlake
  public :: llakmsk1
  public :: lakmsk1
  public :: allocate_mod_bats_lake
  public :: lake_fillvar
#endif

  public :: allocate_land_model
  public :: land_albedo
  public :: land_model
  public :: init_land_model
  public :: initbats
  public :: mtrxbats
#ifdef CLM
  public :: voc_em
  public :: voc_em1
  public :: voc_em2
  public :: dep_vels
  public :: initclm
  public :: mtrxclm
  public :: zenit_clm
#endif
  contains

  subroutine allocate_land_model
    implicit none

    rrnnsg = 1.0/real(nnsg)
    rdnnsg = d_one/dble(nnsg)

    call allocate_mod_bats_internal

    if ( iocncpl == 1 ) then
      call getmem2d(cplmsk,jci1,jci2,ici1,ici2,'bats:cplmsk')
      cplmsk(:,:) = 0 
      ! This is for the RTM component
      call getmem3d(dailyrnf,jci1,jci2,ici1,ici2,1,2,'bats:dailyrnf')
    end if
    if ( ichem == 1 ) then
      call getmem2d(ssw2da,jci1,jci2,ici1,ici2,'bats:ssw2da')
      call getmem2d(sfracv2d,jci1,jci2,ici1,ici2,'bats:sfracv2d')
      call getmem2d(sfracb2d,jci1,jci2,ici1,ici2,'bats:sfracb2d')
      call getmem2d(sfracs2d,jci1,jci2,ici1,ici2,'bats:sfracs2d')
      call getmem2d(svegfrac2d,jci1,jci2,ici1,ici2,'bats:svegfrac2d')
    end if

    call getmem3d(gwet1,1,nnsg,jci1,jci2,ici1,ici2,'bats:gwet1')
    call getmem3d(rsw1,1,nnsg,jci1,jci2,ici1,ici2,'bats:rsw1')
    call getmem3d(snag1,1,nnsg,jci1,jci2,ici1,ici2,'bats:snag1')
    call getmem3d(sncv1,1,nnsg,jci1,jci2,ici1,ici2,'bats:sncv1')
    call getmem3d(sfice1,1,nnsg,jci1,jci2,ici1,ici2,'bats:sfice1')
    call getmem3d(ssw1,1,nnsg,jci1,jci2,ici1,ici2,'bats:ssw1')
    call getmem3d(tgrd1,1,nnsg,jci1,jci2,ici1,ici2,'bats:tgrd1')
    call getmem3d(tgbrd1,1,nnsg,jci1,jci2,ici1,ici2,'bats:tgbrd1')
    call getmem3d(tlef1,1,nnsg,jci1,jci2,ici1,ici2,'bats:tlef1')
    call getmem3d(tsw1,1,nnsg,jci1,jci2,ici1,ici2,'bats:tsw1')
    call getmem3d(taf1,1,nnsg,jci1,jci2,ici1,ici2,'bats:taf1')
    call getmem3d(ldew1,1,nnsg,jci1,jci2,ici1,ici2,'bats:ldew1')

    if (idcsst == 1) then
      call getmem2d(deltas,jci1,jci2,ici1,ici2,'bats:deltas')
      call getmem2d(tdeltas,jci1,jci2,ici1,ici2,'bats:tdeltas')
      call getmem2d(dtskin,jci1,jci2,ici1,ici2,'bats:dtskin')
      call getmem2d(sst,jci1,jci2,ici1,ici2,'bats:sst')
    end if

    call getmem3d(sent1,1,nnsg,jci1,jci2,ici1,ici2,'bats:sent1')
    call getmem3d(evpr1,1,nnsg,jci1,jci2,ici1,ici2,'bats:evpr1')
    call getmem3d(drag1,1,nnsg,jci1,jci2,ici1,ici2,'bats:drag1')
    call getmem3d(prcp1,1,nnsg,jci1,jci2,ici1,ici2,'bats:prcp1')
    call getmem3d(q2m,1,nnsg,jci1,jci2,ici1,ici2,'bats:q2m')
    call getmem3d(ps1,1,nnsg,jci1,jci2,ici1,ici2,'bats:ps1')
    call getmem3d(trnof1,1,nnsg,jci1,jci2,ici1,ici2,'bats:trnof1')
    call getmem3d(srnof1,1,nnsg,jci1,jci2,ici1,ici2,'bats:srnof1')
    call getmem3d(t2m,1,nnsg,jci1,jci2,ici1,ici2,'bats:t2m')
    call getmem3d(u10m,1,nnsg,jci1,jci2,ici1,ici2,'bats:u10m')
    call getmem3d(v10m,1,nnsg,jci1,jci2,ici1,ici2,'bats:v10m')
    call getmem3d(taux,1,nnsg,jci1,jci2,ici1,ici2,'bats:taux')
    call getmem3d(tauy,1,nnsg,jci1,jci2,ici1,ici2,'bats:tauy')
    call getmem3d(emiss1,1,nnsg,jci1,jci2,ici1,ici2,'bats:emiss1')

    if ( lakemod == 1 ) then
      call getmem3d(lakmsk1,1,nnsg,jci1,jci2,ici1,ici2,'bats:lakmsk1')
      call getmem3d(llakmsk1,1,nnsg,jci1,jci2,ici1,ici2,'bats:llakmsk1')
      call getmem3d(xlake,1,nnsg,jci1,jci2,ici1,ici2,'bats:xlake')
      call getmem4d(tlake,1,nnsg,jci1,jci2,ici1,ici2,1,ndpmax,'bats:tlake')
    end if

#ifdef CLM
    call getmem2d(r2ctb,1,jxp,1,iyp,'clm:r2ctb')
    call getmem2d(r2cqb,1,jxp,1,iyp,'clm:r2cqb')
    call getmem2d(r2czga,1,jxp,1,iyp,'clm:r2czga')
    call getmem2d(r2cpsb,1,jxp,1,iyp,'clm:r2cpsb')
    call getmem2d(r2cuxb,1,jxp,1,iyp,'clm:r2cuxb')
    call getmem2d(r2cvxb,1,jxp,1,iyp,'clm:r2cvxb')
    call getmem2d(r2crnc,1,jxp,1,iyp,'clm:r2crnc')
    call getmem2d(r2crnnc,1,jxp,1,iyp,'clm:r2crnnc')
    call getmem2d(r2csols,1,jxp,1,iyp,'clm:r2csols')
    call getmem2d(r2csoll,1,jxp,1,iyp,'clm:r2csoll')
    call getmem2d(r2csolsd,1,jxp,1,iyp,'clm:r2csolsd')
    call getmem2d(r2csolld,1,jxp,1,iyp,'clm:r2csolld')
    call getmem2d(r2cflwd,1,jxp,1,iyp,'clm:r2cflwd')
    call getmem2d(r2cxlat,1,jxp,1,iyp,'clm:r2cxlat')
    call getmem2d(r2cxlon,1,jxp,1,iyp,'clm:r2cxlon')
    call getmem2d(r2cxlatd,1,jxp,1,iyp,'clm:r2cxlatd')
    call getmem2d(r2cxlond,1,jxp,1,iyp,'clm:r2cxlond')

    call getmem2d(r2ctb_all,1,jx,1,iy,'clm:r2ctb_all')
    call getmem2d(r2cqb_all,1,jx,1,iy,'clm:r2cqb_all')
    call getmem2d(r2czga_all,1,jx,1,iy,'clm:r2czga_all')
    call getmem2d(r2cpsb_all,1,jx,1,iy,'clm:r2cpsb_all')
    call getmem2d(r2cuxb_all,1,jx,1,iy,'clm:r2cuxb_all')
    call getmem2d(r2cvxb_all,1,jx,1,iy,'clm:r2cvxb_all')
    call getmem2d(r2crnc_all,1,jx,1,iy,'clm:r2crnc_all')
    call getmem2d(r2crnnc_all,1,jx,1,iy,'clm:r2crnnc_all')
    call getmem2d(r2csols_all,1,jx,1,iy,'clm:r2csols_all')
    call getmem2d(r2csoll_all,1,jx,1,iy,'clm:r2csoll_all')
    call getmem2d(r2csolsd_all,1,jx,1,iy,'clm:r2csolsd_all')
    call getmem2d(r2csolld_all,1,jx,1,iy,'clm:r2csolld_all')
    call getmem2d(r2cflwd_all,1,jx,1,iy,'clm:r2cflwd_all')
    call getmem2d(r2ccosz_all,1,jx,1,iy,'clm:r2ccosz_all')
    call getmem2d(r2cxlat_all,1,jx,1,iy,'clm:r2cxlat_all')
    call getmem2d(r2cxlon_all,1,jx,1,iy,'clm:r2cxlon_all')
    call getmem2d(r2cxlatd_all,1,jx,1,iy,'clm:r2cxlatd_all')
    call getmem2d(r2cxlond_all,1,jx,1,iy,'clm:r2cxlond_all')

    call getmem2d(c2rtgb,1,jx,1,iy,'clm:c2rtgb')
    call getmem2d(c2rsenht,1,jx,1,iy,'clm:c2rsenht')
    call getmem2d(c2rlatht,1,jx,1,iy,'clm:c2rlatht')
    call getmem2d(c2ralbdirs,1,jx,1,iy,'clm:c2ralbdirs')
    call getmem2d(c2ralbdirl,1,jx,1,iy,'clm:c2ralbdirl')
    call getmem2d(c2ralbdifs,1,jx,1,iy,'clm:c2ralbdifs')
    call getmem2d(c2ralbdifl,1,jx,1,iy,'clm:c2ralbdifl')
    call getmem2d(c2rtaux,1,jx,1,iy,'clm:c2rtaux')
    call getmem2d(c2rtauy,1,jx,1,iy,'clm:c2rtauy')
    call getmem2d(c2ruvdrag,1,jx,1,iy,'clm:c2ruvdrag')
    call getmem2d(c2rlsmask,1,jx,1,iy,'clm:c2rlsmask')
    call getmem2d(c2rtgbb,1,jx,1,iy,'clm:c2rtgbb')
    call getmem2d(c2rsnowc,1,jx,1,iy,'clm:c2rsnowc')
    call getmem2d(c2rtest,1,jx,1,iy,'clm:c2rtest')
    call getmem2d(c2r2mt,1,jx,1,iy,'clm:c2r2mt')
    call getmem2d(c2r2mq,1,jx,1,iy,'clm:c2r2mq')
    call getmem2d(c2rtlef,1,jx,1,iy,'clm:c2rtlef')
    call getmem2d(c2ru10,1,jx,1,iy,'clm:c2ru10')
    call getmem2d(c2rsm10cm,1,jx,1,iy,'clm:c2rsm10cm')
    call getmem2d(c2rsm1m,1,jx,1,iy,'clm:c2rsm1m')
    call getmem2d(c2rsmtot,1,jx,1,iy,'clm:c2rsmtot')
    call getmem2d(c2rinfl,1,jx,1,iy,'clm:c2rinfl')
    call getmem2d(c2rro_sur,1,jx,1,iy,'clm:c2rro_sur')
    call getmem2d(c2rro_sub,1,jx,1,iy,'clm:c2rro_sub')
    call getmem2d(c2rfracsno,1,jx,1,iy,'clm:c2rfracsno')
    call getmem2d(c2rfvegnosno,1,jx,1,iy,'clm:c2rfvegnosno')
    call getmem2d(c2rprocmap,1,jx,1,iy,'clm:c2rprocmap')
#if (defined VOC)
    call getmem2d(voc_em,1,jx,1,iy,'clm:voc_em')
    call getmem2d(voc_em1,1,jx,1,iy,'clm:voc_em1')
    call getmem2d(voc_em2,1,jx,1,iy,'clm:voc_em2')
#endif
    if ( igaschem == 1 .or. ioxclim == 1 ) then
      call getmem3d(dep_vels,1,jx,1,iy,1,ntr,'clm:dep_vels')
    end if
    call getmem1d(c2rngc,1,nproc,'clm:c2rngc')
    call getmem1d(c2rdisps,1,nproc,'clm:c2rdisps')
    call getmem2d(rs2d,jci1,jci2,ici1,ici2,'clm:rs2d')
    call getmem2d(ra2d,jci1,jci2,ici1,ici2,'clm:ra2d')
#endif
  end subroutine allocate_land_model

  subroutine init_land_model
    use mod_atm_interface
    use mod_che_interface
    implicit none
    ntcpl  = idnint(cpldt/dtsec)
    ntsrf2 = idnint(dtsrf/dtsec)
    if ( idcsst   == 1 ) ldcsst    = .true.
    if ( lakemod  == 1 ) llake     = .true.
    if ( idesseas == 1 ) ldesseas  = .true.
    if ( iseaice  == 1 ) lseaice   = .true.
    call assignpnt(mddom%xlat,xlat)
    call assignpnt(mddom%xlon,xlon)
    call assignpnt(mddom%lndcat,lndcat)
    call assignpnt(mddom%ldmsk,ldmsk)
    call assignpnt(mddom%iveg,iveg)
    call assignpnt(mddom%ht,ht)
    call assignpnt(mddom%snowam,snowam)
    call assignpnt(mdsub%xlat,xlat1)
    call assignpnt(mdsub%xlon,xlon1)
    call assignpnt(mdsub%lndcat,lndcat1)
    call assignpnt(mdsub%ldmsk,ldmsk1)
    call assignpnt(mdsub%ht,ht1)
    call assignpnt(mdsub%dhlake,dhlake1)
    call assignpnt(mdsub%iveg,iveg1)
    call assignpnt(atms%ubx3d,uatm,kz)
    call assignpnt(atms%vbx3d,vatm,kz)
    call assignpnt(atms%tb3d,tatm,kz)
    call assignpnt(atms%qxb3d,qvatm,kz,iqv)
    call assignpnt(atms%thx3d,thatm,kz)
    call assignpnt(atms%rhox2d,rhox)
    call assignpnt(atms%za,hgt,kz)
    call assignpnt(sfs%hfx,hfx)
    call assignpnt(sfs%qfx,qfx)
    call assignpnt(sfs%uvdrag,uvdrag)
    call assignpnt(sfs%tgbb,tgbb)
    call assignpnt(sfs%psb,sfps)
    call assignpnt(sfs%tga,tground1)
    call assignpnt(sfs%tgb,tground2)
    call assignpnt(zpbl,hpbl)
    call assignpnt(pptc,cprate)
    call assignpnt(pptnc,ncprate)
    call assignpnt(coszrs,zencos)
    call assignpnt(fsw,rswf)
    call assignpnt(flw,rlwf)
    call assignpnt(flwd,dwrlwf)
    call assignpnt(sabveg,vegswab)
    call assignpnt(albvl,lwalb)
    call assignpnt(albvs,swalb)
    call assignpnt(aldirs,swdiralb)
    call assignpnt(aldifs,swdifalb)
    call assignpnt(aldirl,lwdiralb)
    call assignpnt(aldifl,lwdifalb)
    call assignpnt(solis,solar)
    call assignpnt(emiss,emissivity)
    call assignpnt(sinc,solinc)
    call assignpnt(solvs,swdir)
    call assignpnt(solvsd,swdif)
    call assignpnt(solvl,lwdir)
    call assignpnt(solvld,lwdif)
    call assignpnt(sdelq,deltaq)
    call assignpnt(sdelt,deltat)
#ifdef CLM
    allocate(landmask(jx,iy))
    call assignpnt(landmask,lmask)
#endif
  end subroutine init_land_model

  subroutine land_model
    implicit none
#ifdef CLM
    if ( ktau == 0 .or. mod(ktau+1,ntrad) == 0 ) then
      r2cdoalb = .true.
    else
      r2cdoalb = .false.
    end if
    ! Timestep used is the same as for bats
    if ( ktau == 0 ) then
      r2cnstep = 0
    else
      r2cnstep = (ktau+1)/ntsrf
    end if
    call mtrxclm
#else
    call mtrxbats
#endif
  end subroutine land_model

  subroutine land_albedo
    implicit none
#ifdef CLM
    call albedoclm
#else
    call albedobats
#endif
  end subroutine land_albedo

end module mod_lm_interface

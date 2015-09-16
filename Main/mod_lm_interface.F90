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
  use mod_runparams
  use mod_memutil
  use mod_regcm_types
  use mod_outvars
  use mod_mppparam
  use mod_mpmessage
  use mod_service
  use mod_bats_common
  use mod_ocn_common
  use mod_che_common
#ifdef CLM
  use mod_clm
  use mod_mtrxclm
  use clm_varsur , only : landmask , numdays
  use clm_varctl , only : filer_rest
  use clm_time_manager , only : get_step_size
  use restFileMod , only : restFile_write, restFile_write_binary
  use restFileMod , only : restFile_filename
  use spmdMod , only : mpicom
  use perf_mod , only : t_prf , t_finalizef
#endif

#ifdef CLM45
  use mod_clm_regcm
#endif

  implicit none

  private

  ! Coupling variables
  real(rk8) :: runoffcount = 0.0D0
  public :: lms

  public :: dtbat

  public :: ocncomm , lndcomm

  public :: import_data_into_surface
  public :: export_data_from_surface

  public :: allocate_surface_model
  public :: surface_albedo
  public :: surface_model
  public :: init_surface_model
  public :: initialize_surface_model
#ifdef CLM
  public :: get_step_size
  public :: filer_rest
  public :: restFile_write
  public :: restFile_write_binary
  public :: restFile_filename
  public :: numdays
  public :: r2ceccf
  public :: mpicom
  public :: t_prf
  public :: t_finalizef
#endif

  type(lm_exchange) :: lm
  type(lm_state) :: lms

  contains

  subroutine allocate_surface_model
    implicit none

    rdnnsg = d_one/dble(nnsg)

    call getmem3d(lms%sent,1,nnsg,jci1,jci2,ici1,ici2,'bats:sent')
    call getmem3d(lms%evpr,1,nnsg,jci1,jci2,ici1,ici2,'bats:evpr')
    call getmem3d(lms%deltat,1,nnsg,jci1,jci2,ici1,ici2,'bats:deltat')
    call getmem3d(lms%deltaq,1,nnsg,jci1,jci2,ici1,ici2,'bats:deltaq')
    call getmem3d(lms%drag,1,nnsg,jci1,jci2,ici1,ici2,'bats:drag')
    call getmem3d(lms%lncl,1,nnsg,jci1,jci2,ici1,ici2,'bats:lncl')
    call getmem3d(lms%prcp,1,nnsg,jci1,jci2,ici1,ici2,'bats:prcp')
    call getmem3d(lms%snwm,1,nnsg,jci1,jci2,ici1,ici2,'bats:snwm')
    call getmem3d(lms%trnof,1,nnsg,jci1,jci2,ici1,ici2,'bats:trnof')
    call getmem3d(lms%srnof,1,nnsg,jci1,jci2,ici1,ici2,'bats:srnof')
    call getmem3d(lms%xlai,1,nnsg,jci1,jci2,ici1,ici2,'bats:xlai')
    call getmem3d(lms%sfcp,1,nnsg,jci1,jci2,ici1,ici2,'bats:sfcp')
    call getmem3d(lms%q2m,1,nnsg,jci1,jci2,ici1,ici2,'bats:q2m')
    call getmem3d(lms%t2m,1,nnsg,jci1,jci2,ici1,ici2,'bats:t2m')
    call getmem3d(lms%u10m,1,nnsg,jci1,jci2,ici1,ici2,'bats:u10m')
    call getmem3d(lms%v10m,1,nnsg,jci1,jci2,ici1,ici2,'bats:v10m')
    call getmem3d(lms%taux,1,nnsg,jci1,jci2,ici1,ici2,'bats:taux')
    call getmem3d(lms%tauy,1,nnsg,jci1,jci2,ici1,ici2,'bats:tauy')
    call getmem3d(lms%wt,1,nnsg,jci1,jci2,ici1,ici2,'bats:wt')
    call getmem3d(lms%swalb,1,nnsg,jci1,jci2,ici1,ici2,'bats:swalb')
    call getmem3d(lms%lwalb,1,nnsg,jci1,jci2,ici1,ici2,'bats:lwalb')
    call getmem3d(lms%swdiralb,1,nnsg,jci1,jci2,ici1,ici2,'bats:swdiralb')
    call getmem3d(lms%lwdiralb,1,nnsg,jci1,jci2,ici1,ici2,'bats:lwdiralb')
    call getmem3d(lms%swdifalb,1,nnsg,jci1,jci2,ici1,ici2,'bats:swdifalb')
    call getmem3d(lms%lwdifalb,1,nnsg,jci1,jci2,ici1,ici2,'bats:lwdifalb')

    call getmem3d(lms%gwet,1,nnsg,jci1,jci2,ici1,ici2,'bats:gwet')
    call getmem3d(lms%ldew,1,nnsg,jci1,jci2,ici1,ici2,'bats:ldew')
    call getmem4d(lms%sw,1,nnsg,jci1,jci2,ici1,ici2,1,num_soil_layers,'bats:sw')
    call assignpnt(lms%sw,lms%ssw,1)
    call assignpnt(lms%sw,lms%rsw,2)
    call getmem3d(lms%tsw,1,nnsg,jci1,jci2,ici1,ici2,'bats:tsw')
    call getmem3d(lms%tgbb,1,nnsg,jci1,jci2,ici1,ici2,'bats:tgbb')
    call getmem3d(lms%tgrd,1,nnsg,jci1,jci2,ici1,ici2,'bats:tgrd')
    call getmem3d(lms%tgbrd,1,nnsg,jci1,jci2,ici1,ici2,'bats:tgbrd')
    call getmem3d(lms%tlef,1,nnsg,jci1,jci2,ici1,ici2,'bats:tlef')
    call getmem3d(lms%taf,1,nnsg,jci1,jci2,ici1,ici2,'bats:taf')
    call getmem3d(lms%sigf,1,nnsg,jci1,jci2,ici1,ici2,'bats:sigf')
    call getmem3d(lms%sfice,1,nnsg,jci1,jci2,ici1,ici2,'bats:sfice')
    call getmem3d(lms%snag,1,nnsg,jci1,jci2,ici1,ici2,'bats:snag')
    call getmem3d(lms%sncv,1,nnsg,jci1,jci2,ici1,ici2,'bats:sncv')
    call getmem3d(lms%scvk,1,nnsg,jci1,jci2,ici1,ici2,'bats:scvk')
    call getmem3d(lms%um10,1,nnsg,jci1,jci2,ici1,ici2,'bats:um10')
    call getmem3d(lms%emisv,1,nnsg,jci1,jci2,ici1,ici2,'bats:emisv')
    call getmem4d(lms%vocemiss,1,nnsg,jci1,jci2,ici1,ici2,1,ntr,'bats:vocemiss')
    call getmem4d(lms%dustemiss,1,nnsg,jci1,jci2,ici1,ici2,1,4,'bats:dustemiss')

#ifdef CLM
    call getmem2d(r2ctb,jde1,jde2,ide1,ide2,'clm:r2ctb')
    call getmem2d(r2cqb,jde1,jde2,ide1,ide2,'clm:r2cqb')
    call getmem2d(r2czga,jde1,jde2,ide1,ide2,'clm:r2czga')
    call getmem2d(r2cpsb,jde1,jde2,ide1,ide2,'clm:r2cpsb')
    call getmem2d(r2cuxb,jde1,jde2,ide1,ide2,'clm:r2cuxb')
    call getmem2d(r2cvxb,jde1,jde2,ide1,ide2,'clm:r2cvxb')
    call getmem2d(r2crnc,jde1,jde2,ide1,ide2,'clm:r2crnc')
    call getmem2d(r2crnnc,jde1,jde2,ide1,ide2,'clm:r2crnnc')
    call getmem2d(r2csols,jde1,jde2,ide1,ide2,'clm:r2csols')
    call getmem2d(r2csoll,jde1,jde2,ide1,ide2,'clm:r2csoll')
    call getmem2d(r2csolsd,jde1,jde2,ide1,ide2,'clm:r2csolsd')
    call getmem2d(r2csolld,jde1,jde2,ide1,ide2,'clm:r2csolld')
    call getmem2d(r2cflwd,jde1,jde2,ide1,ide2,'clm:r2cflwd')
    call getmem2d(r2cxlat,jde1,jde2,ide1,ide2,'clm:r2cxlat')
    call getmem2d(r2cxlon,jde1,jde2,ide1,ide2,'clm:r2cxlon')
    call getmem2d(r2cxlatd,jde1,jde2,ide1,ide2,'clm:r2cxlatd')
    call getmem2d(r2cxlond,jde1,jde2,ide1,ide2,'clm:r2cxlond')

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
    call getmem1d(c2rngc,1,nproc,'clm:c2rngc')
    call getmem1d(c2rdisps,1,nproc,'clm:c2rdisps')
    call getmem2d(rs2d,jci1,jci2,ici1,ici2,'clm:rs2d')
    call getmem2d(ra2d,jci1,jci2,ici1,ici2,'clm:ra2d')
#else
    if ( lakemod == 1 ) then
      call getmem3d(lms%eta,1,nnsg,jci1,jci2,ici1,ici2,'lake:eta')
      call getmem3d(lms%hi,1,nnsg,jci1,jci2,ici1,ici2,'lake:hi')
      call getmem3d(lms%lakmsk,1,nnsg,jci1,jci2,ici1,ici2,'lake:lakmsk')
      call getmem4d(lms%tlake,1,nnsg,jci1,jci2,ici1,ici2,1,ndpmax,'lake:tlake')
    end if
#endif
    if ( idcsst == 1 ) then
      call getmem3d(lms%deltas,1,nnsg,jci1,jci2,ici1,ici2,'sst:deltas')
      call getmem3d(lms%tdeltas,1,nnsg,jci1,jci2,ici1,ici2,'sst:tdeltas')
      call getmem3d(lms%tskin,1,nnsg,jci1,jci2,ici1,ici2,'sst:tskin')
      call getmem3d(lms%sst,1,nnsg,jci1,jci2,ici1,ici2,'sst:sst')
    end if
  end subroutine allocate_surface_model

  subroutine init_surface_model
    use mod_atm_interface
    use mod_che_interface
    implicit none

    call cl_setup(lndcomm,mddom%mask,mdsub%mask)
    call cl_setup(ocncomm,mddom%mask,mdsub%mask,.true.)

#ifdef DEBUG
    write(ndebug+myid,*) 'TOTAL POINTS FOR LAND  IN LNDCOMM : ', &
      lndcomm%linear_npoint_sg(myid+1)
    write(ndebug+myid,*) 'Cartesian p ', lndcomm%cartesian_npoint_g
    write(ndebug+myid,*) 'Cartesian d ', lndcomm%cartesian_displ_g
    write(ndebug+myid,*) 'Linear    p ', lndcomm%linear_npoint_g
    write(ndebug+myid,*) 'Linear    d ', lndcomm%linear_displ_g
    write(ndebug+myid,*) 'Subgrid Cartesian p ', lndcomm%cartesian_npoint_sg
    write(ndebug+myid,*) 'Subgrid Cartesian d ', lndcomm%cartesian_displ_sg
    write(ndebug+myid,*) 'Subgrid Linear    p ', lndcomm%linear_npoint_sg
    write(ndebug+myid,*) 'Subgrid Linear    d ', lndcomm%linear_displ_sg
    write(ndebug+myid,*) 'TOTAL POINTS FOR OCEAN IN OCNCOMM : ', &
      ocncomm%linear_npoint_sg(myid+1)
    write(ndebug+myid,*) 'Cartesian p ', ocncomm%cartesian_npoint_g
    write(ndebug+myid,*) 'Cartesian d ', ocncomm%cartesian_displ_g
    write(ndebug+myid,*) 'Linear    p ', ocncomm%linear_npoint_g
    write(ndebug+myid,*) 'Linear    d ', ocncomm%linear_displ_g
    write(ndebug+myid,*) 'Subgrid Cartesian p ', ocncomm%cartesian_npoint_sg
    write(ndebug+myid,*) 'Subgrid Cartesian d ', ocncomm%cartesian_displ_sg
    write(ndebug+myid,*) 'Subgrid Linear    p ', ocncomm%linear_npoint_sg
    write(ndebug+myid,*) 'Subgrid Linear    d ', ocncomm%linear_displ_sg
#endif

    call allocate_mod_bats_internal(lndcomm)
    call allocate_mod_ocn_internal(ocncomm)

    ntcpl  = idnint(cpldt/dtsec)
    if ( idcsst   == 1 ) ldcsst   = .true.
    if ( lakemod  == 1 ) llake    = .true.
    if ( idesseas == 1 ) ldesseas = .true.
    if ( iseaice  == 1 ) lseaice  = .true.
    if ( iocncpl  == 1 ) lcoup    = .true.
    call assignpnt(mddom%xlat,lm%xlat)
    call assignpnt(mddom%xlon,lm%xlon)
    call assignpnt(mddom%lndcat,lm%lndcat)
    call assignpnt(mddom%ldmsk,lm%ldmsk)
    call assignpnt(mddom%iveg,lm%iveg)
    call assignpnt(mddom%ht,lm%ht)
    call assignpnt(mddom%snowam,lm%snowam)
    call assignpnt(mddom%smoist,lm%smoist)
    call assignpnt(mddom%rmoist,lm%rmoist)
    call assignpnt(mdsub%xlat,lm%xlat1)
    call assignpnt(mdsub%xlon,lm%xlon1)
    call assignpnt(mdsub%lndcat,lm%lndcat1)
    call assignpnt(mdsub%ldmsk,lm%ldmsk1)
    call assignpnt(mdsub%ht,lm%ht1)
    call assignpnt(mdsub%iveg,lm%iveg1)
    call assignpnt(mdsub%dhlake,lm%dhlake1)
    call assignpnt(cplmsk,lm%icplmsk)
    call assignpnt(atms%ubx3d,lm%uatm,kz)
    call assignpnt(atms%vbx3d,lm%vatm,kz)
    call assignpnt(atms%th3d,lm%thatm,kz)
    call assignpnt(atms%tb3d,lm%tatm,kz)
    call assignpnt(atms%pb3d,lm%patm,kz)
    call assignpnt(atms%qxb3d,lm%qvatm,kz,iqv)
    call assignpnt(atms%rhox2d,lm%rhox)
    call assignpnt(atms%ps2d,lm%sfps)
    call assignpnt(atms%ts2d,lm%sfta)
    call assignpnt(atms%za,lm%hgt,kz)
    call assignpnt(sfs%hfx,lm%hfx)
    call assignpnt(sfs%qfx,lm%qfx)
    call assignpnt(sfs%uvdrag,lm%uvdrag)
    call assignpnt(sfs%tgbb,lm%tgbb)
    call assignpnt(sfs%tga,lm%tground1)
    call assignpnt(sfs%tgb,lm%tground2)
    call assignpnt(zpbl,lm%hpbl)
    call assignpnt(pptc,lm%cprate)
    call assignpnt(pptnc,lm%ncprate)
    call assignpnt(coszrs,lm%zencos)
    call assignpnt(fsw,lm%rswf)
    call assignpnt(flw,lm%rlwf)
    call assignpnt(flwd,lm%dwrlwf)
    call assignpnt(sabveg,lm%vegswab)
    call assignpnt(albvl,lm%lwalb)
    call assignpnt(albvs,lm%swalb)
    call assignpnt(aldirs,lm%swdiralb)
    call assignpnt(aldifs,lm%swdifalb)
    call assignpnt(aldirl,lm%lwdiralb)
    call assignpnt(aldifl,lm%lwdifalb)
    call assignpnt(solis,lm%solar)
    call assignpnt(emiss,lm%emissivity)
    call assignpnt(sinc,lm%solinc)
    call assignpnt(solvs,lm%swdir)
    call assignpnt(solvsd,lm%swdif)
    call assignpnt(solvl,lm%lwdir)
    call assignpnt(solvld,lm%lwdif)
    call assignpnt(sdelq,lm%deltaq)
    call assignpnt(sdelt,lm%deltat)
    if ( ichem == 1 ) then
      call assignpnt(sdelt,lm%deltat)
      call assignpnt(sdelq,lm%deltaq)
      call assignpnt(ssw2da,lm%ssw2da)
      call assignpnt(sfracv2d,lm%sfracv2d)
      call assignpnt(sfracb2d,lm%sfracb2d)
      call assignpnt(sfracs2d,lm%sfracs2d)
      call assignpnt(svegfrac2d,lm%svegfrac2d)
      call assignpnt(sxlai2d,lm%sxlai2d)
      call assignpnt(wetdepflx,lm%wetdepflx)
      call assignpnt(drydepflx,lm%drydepflx)
      call assignpnt(idusts,lm%idust)
    end if
    call assignpnt(dailyrnf,lm%dailyrnf)
#ifdef CLM
    allocate(landmask(jx,iy))
    if ( ichem == 1 ) then
      if ( igaschem == 1 .or. ioxclim == 1 ) then
        call assignpnt(dep_vels,lm%dep_vels)
      end if
#ifdef VOC
      call assignpnt(voc_em0,lm%voc_em0)
      call assignpnt(voc_em1,lm%voc_em1)
      call assignpnt(voc_em2,lm%voc_em2)
#endif
    end if
#endif
  end subroutine init_surface_model

  subroutine initialize_surface_model
    implicit none
#ifdef CLM
    integer(ik4) :: i , j , n
#endif
#ifndef CLM45
    call initbats(lm,lms)
#else
    call initclm45(lm,lms)
#endif
    call initocn(lm,lms)
#ifdef CLM
    call initclm(lm,lms)
    if ( ktau == 0 .and. imask == 2 ) then
      ! CLM may have changed the landuse again !
      do i = ici1 , ici2
        do j = jci1 , jci2
          lm%iveg(j,i) = idnint(lm%lndcat(j,i))
        end do
      end do
      ! Correct land/water misalign : set to short grass
      where ( (lm%iveg == 14 .or. lm%iveg == 15) .and. lm%ldmsk == 1 )
        lm%iveg = 2
        lm%lndcat = d_two
      end where
      do i = ici1 , ici2
        do j = jci1 , jci2
          do n = 1 , nnsg
            lm%iveg1(n,j,i) = lm%iveg(j,i)
            lm%lndcat1(n,j,i) = lm%lndcat(j,i)
          end do
        end do
      end do
    end if
#endif
    lm%emissivity = sum(lms%emisv,1) * rdnnsg
  end subroutine initialize_surface_model

  subroutine surface_model
    implicit none
    integer(ik4) :: i , j , n , nn , ierr , i1 , i2
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
      r2cnstep = int((ktau+1)/ntsrf,ik4)
    end if
    call mtrxclm(lm,lms)
#else
#ifdef CLM45
    call runclm45(lm,lms)
    !coupling of biogenic VOC from CLM45 to chemistry
    if ( ichem == 1 ) then
      do n = 1 , ntr
        do i = ici1 , ici2
          do j = jci1 , jci2
            cvoc_em(j,i,n) = sum(lms%vocemiss(:,j,i,n),1) * rdnnsg
          end do
        end do
      end do

      do n = 1 , 4 
        do i = ici1 , ici2
          do j = jci1 , jci2
            cdustflx_clm(j,i,n) = sum(lms%dustemiss(:,j,i,n),1) * rdnnsg
          end do
        end do
    end if
#else
    call vecbats(lm,lms)
#endif
#endif
    call vecocn(lm,lms)
    lm%hfx = sum(lms%sent,1)*rdnnsg
    lm%qfx = sum(lms%evpr,1)*rdnnsg
    lm%uvdrag = sum(lms%drag,1)*rdnnsg
    lm%tgbb = sum(lms%tgbb,1)*rdnnsg
    lm%tground1 = sum(lms%tgrd,1)*rdnnsg
    lm%tground2 = sum(lms%tgrd,1)*rdnnsg
    lm%emissivity = sum(lms%emisv,1) * rdnnsg
    if ( iseaice == 1 .or. lakemod == 1 ) then
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( lm%ldmsk(j,i) /= 1 ) then
            i1 = count(lm%ldmsk1(:,j,i) == 0)
            i2 = count(lm%ldmsk1(:,j,i) == 2)
            if ( i1 > i2 ) then
              lm%ldmsk(j,i) = 0
            else
              lm%ldmsk(j,i) = 2
            end if
          end if
        end do
      end do
    end if
    if ( ichem == 1 ) then
#ifdef CLM
      do i = ici1 , ici2
        do j = jci1 , jci2
          lm%deltat(j,i) = sum(lms%deltat(:,j,i))*rdnnsg
          lm%deltaq(j,i) = sum(lms%deltaq(:,j,i))*rdnnsg
          lm%sfracv2d(j,i) = c2rfvegnosno(j,i)
          lm%sfracb2d(j,i) = d_one - (c2rfvegnosno(j,i)+c2rfracsno(j,i))
          lm%sfracs2d(j,i) = c2rfracsno(j,i)
          lm%ssw2da = sum(lms%ssw,1)*rdnnsg
          lm%sxlai2d = 0.0D0
        end do
      end do
#else
      lm%deltat = sum(lms%deltat,1)*rdnnsg
      lm%deltaq = sum(lms%deltaq,1)*rdnnsg
      lm%sfracv2d = sum(lms%sigf,1)*rdnnsg
      lm%svegfrac2d = sum(lms%lncl,1)*rdnnsg
      lm%sxlai2d = sum(lms%xlai,1)*rdnnsg
      lm%sfracb2d = sum(((d_one-lms%lncl)*(d_one-lms%scvk)),1)*rdnnsg
      lm%sfracs2d = sum((lms%lncl*lms%wt+(d_one-lms%lncl)*lms%scvk),1)*rdnnsg
      lm%ssw2da = sum(lms%ssw,1)*rdnnsg
#endif
    end if
    call collect_output
#ifdef DEBUG
    ! Sanity check of surface temperatures
    ierr = 0
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          if ( lms%tgrd(n,j,i) < 150.0D0 ) then
            write(stderr,*) 'Likely error: Surface temperature too low'
            write(stderr,*) 'MYID = ', myid
            write(stderr,*) 'J    = ',j
            write(stderr,*) 'I    = ',i
            do nn = 1 , nnsg
              write(stderr,*) 'N    = ',nn
              write(stderr,*) 'VAL  = ',lms%tgrd(nn,j,i)
              write(stderr,*) 'MASK = ',lm%ldmsk1(n,j,i)
            end do
            write(stderr,*) 'VAL2 = ',lm%tground1(j,i)
            write(stderr,*) 'VAL2 = ',lm%tground2(j,i)
            ierr = ierr + 1
          end if
        end do
      end do
    end do
    if ( ierr /= 0 ) then
      call fatal(__FILE__,__LINE__,'TEMP CHECK ERROR')
    end if
#endif
  end subroutine surface_model

  subroutine surface_albedo
    implicit none
#ifdef CLM
    logical :: do_call_albedo_bats_for_clm = .false.
    if ( do_call_albedo_bats_for_clm ) then
      call albedobats(lm,lms)
    else
      call albedoclm(lm,lms)
    end if
#else
#ifdef CLM45
    call albedoclm45(lm,lms)
#else
    call albedobats(lm,lms)
#endif
#endif
    call albedoocn(lm,lms)
    lm%swalb = sum(lms%swalb,1)*rdnnsg
    lm%lwalb = sum(lms%lwalb,1)*rdnnsg
    lm%swdiralb = sum(lms%swdiralb,1)*rdnnsg
    lm%lwdiralb = sum(lms%lwdiralb,1)*rdnnsg
    lm%swdifalb = sum(lms%swdifalb,1)*rdnnsg
    lm%lwdifalb = sum(lms%lwdifalb,1)*rdnnsg
  end subroutine surface_albedo

  subroutine export_data_from_surface(expfie)
    implicit none
    type(exp_data) , intent(inout) :: expfie
    integer(ik4) :: j , i
    do i = ici1 , ici2
      do j = jci1 , jci2
        expfie%psfc(j,i) = lm%sfps(j,i)*d_r100
        expfie%tsfc(j,i) = sum(lms%t2m(:,j,i))*rdnnsg
        expfie%qsfc(j,i) = sum(lms%q2m(:,j,i))*rdnnsg
        expfie%swrd(j,i) = lm%rswf(j,i)
        expfie%lwrd(j,i) = lm%rlwf(j,i)
        expfie%dlwr(j,i) = lm%dwrlwf(j,i)
        expfie%dswr(j,i) = lm%swdif(j,i)+lm%swdir(j,i)
        expfie%lhfx(j,i) = sum(lms%evpr(:,j,i))*rdnnsg*wlhv
        expfie%shfx(j,i) = sum(lms%sent(:,j,i))*rdnnsg
        expfie%prec(j,i) = sum(lms%prcp(:,j,i))*rdnnsg
        expfie%wndu(j,i) = sum(lms%u10m(:,j,i))*rdnnsg
        expfie%wndv(j,i) = sum(lms%v10m(:,j,i))*rdnnsg
        expfie%taux(j,i) = sum(lms%taux(:,j,i))*rdnnsg
        expfie%tauy(j,i) = sum(lms%tauy(:,j,i))*rdnnsg
        expfie%sflx(j,i) = (sum(lms%evpr(:,j,i))-sum(lms%prcp(:,j,i)))*rdnnsg
        expfie%snow(j,i) = sum(lms%sncv(:,j,i))*rdnnsg
        expfie%wspd(j,i) = dsqrt(expfie%wndu(j,i)**2+expfie%wndv(j,i)**2)
        expfie%nflx(j,i) = lm%rswf(j,i) - expfie%lhfx(j,i) - &
                           expfie%shfx(j,i) - lm%rlwf(j,i)
      end do
    end do
    if ( mod(ktau+1,kday) == 0 ) then
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( lm%ldmsk(j,i) > 0 ) then
            expfie%rnof(j,i) = lm%dailyrnf(j,i,1)/runoffcount
            expfie%snof(j,i) = lm%dailyrnf(j,i,2)/runoffcount
          else
           expfie%rnof(j,i) = d_zero
           expfie%snof(j,i) = d_zero
          end if
        end do
      end do
      runoffcount = d_zero
      lm%dailyrnf(:,:,:) = d_zero
    end if
  end subroutine export_data_from_surface
!
  subroutine import_data_into_surface(impfie,ldmskb,wetdry,tol)
    use mod_atm_interface
    implicit none
    type(imp_data) , intent(in) :: impfie
    real(rk8) , intent(in) :: tol
    integer(ik4) , pointer , dimension(:,:) , intent(in) :: ldmskb , wetdry
    integer :: i , j , n
    logical :: flag = .false.
    ! real(rk8) :: toth
    ! real(rk8) , parameter :: href = d_two * iceminh
    ! real(rk8) , parameter :: steepf = 1.0D0
    ! integer(ik4) :: ix , jy , imin , imax , jmin , jmax , srad , hveg(22)
    !
    !-----------------------------------------------------------------------
    ! Retrieve information from OCN component
    !-----------------------------------------------------------------------
    !
    do i = ici1, ici2
      do j = jci1, jci2
        if ( lm%iveg(j,i) == 14 .or. lm%iveg(j,i) == 15 ) then
          !
          !--------------------------------------
          ! Update: Sea Surface Temperature (SST)
          !--------------------------------------
          !
          if ( impfie%sst(j,i) < tol ) then
            ! create fixed coupling mask
            if ( mod(ktau+1, ntcpl*2) == 0 ) then
              cplmsk(j,i) = 1
            end if
            lm%tground1(j,i) = impfie%sst(j,i)
            lm%tground2(j,i) = impfie%sst(j,i)
            lm%tgbb(j,i)     = impfie%sst(j,i)
            lms%tgrd(:,j,i)  = impfie%sst(j,i)
            lms%tgbrd(:,j,i) = impfie%sst(j,i)
            lms%tgbb(:,j,i)  = impfie%sst(j,i)
          end if
          !
          !----------------------------------------------------------
          ! Update: Mask and land-use type (based on dynamic wet-dry)
          !----------------------------------------------------------
          !
!         if (importFields%msk(j,i) .lt. tol .and. ldmskb(j,i) == 0) then
!           if (importFields%msk(j,i) .lt. 1.0) then
!             flag = .false.
!             if (lm%ldmsk(j,i) == 0 .or. &
!                 lm%ldmsk(j,i) == 2) flag = .true.
!             ! set land-sea mask
!             lm%ldmsk(j,i) = 1
!             do n = 1, nnsg
!               lm%ldmsk1(n,j,i) = lm%ldmsk(j,i)
!             end do
!             ! count land-use type in a specified search radius (srad)
!             srad = 10
!             jmin = j-srad
!             if (j-srad < jci1) jmin = jci1
!             jmax = j+srad
!             if (j+srad > jci2) jmax = jci2
!             imin = i-srad
!             if (i-srad < ici1) imin = ici1
!             imax = i+srad
!             if (i+srad > ici2) imax = ici2
!             hveg = 0
!             do ix = imin, imax
!               do jy = jmin, jmax
!                 do n = 1, nnsg
!                   hveg(iveg1(n,jy,ix)) = hveg(iveg1(n,jy,ix))+1
!                 end do
!               end do
!             end do
!             hveg(14) = 0
!             hveg(15) = 0
!             ! set array to store change
!             wetdry(j,i) = 1
!             ! write debug info
!             if (flag) then
!               write(*,20) j, i, 'water', 'land ', lm%ldmsk(j,i)
!             end if
!           else
!             if (lm%ldmsk(j,i) == 1 .and. wetdry(j,i) == 1) then
!               flag = .false.
!               if (lm%ldmskb(j,i) /= lm%ldmsk(j,i)) flag = .true.
!               ! set land-sea mask to its original value
!               lm%ldmsk(j,i) = ldmskb(j,i)
!               do n = 1, nnsg
!                 lm%ldmsk1(n,j,i) = lm%ldmsk(j,i)
!               end do
!               ! set array to store change
!               wetdry(j,i) = 0
!               ! write debug info
!               if (flag) then
!                 write(*,20) j, i, 'land ', 'water', lm%ldmsk(j,i)
!               end if
!             end if
!           end if
!         end if
          !
          !------------------------------------------------------------------
          ! Update: Sea-ice, mask and land-use type (based on sea-ice module)
          !------------------------------------------------------------------
          !
          if ( impfie%sit(j,i) < tol .and. lm%ldmsk(j,i) /= 1 ) then
            if ( impfie%sit(j,i) > iceminh ) then
              flag = .false.
              if ( lm%ldmsk(j,i) == 0 ) flag = .true.
              ! set land-sea mask
              lm%ldmsk(j,i) = 2
              do n = 1, nnsg
                lm%ldmsk1(n,j,i) = 2
                ! set sea ice thikness (in meter)
                lms%sfice(n,j,i) = impfie%sit(j,i)
              end do
              ! write debug info
              if ( flag ) then
                write(*,30) j, i, 'water', 'ice  ', &
                   lm%ldmsk(j,i), lms%sfice(1,j,i)
              end if
            else
              if ( ldmskb(j,i) == 0 .and. lm%ldmsk(j,i) == 2 ) then
                do n = 1, nnsg
                  ! reduce to one tenth surface ice: it should melt away
                  lms%sfice(n,j,i) = lms%sfice(n,j,i)*d_r10
                  ! check that sea ice is melted or not
                  if ( lms%sfice(n,j,i) <= iceminh ) then
                    if ( ldmskb(j,i) /= lm%ldmsk(j,i) ) flag = .true.
                    ! set land-sea mask to its original value
                    lm%ldmsk(j,i) = ldmskb(j,i)
                    lm%ldmsk1(n,j,i) = ldmskb(j,i)
                    ! set land-use type to its original value
                    ! set sea ice thikness (in meter)
                    lms%sfice(n,j,i) = d_zero
                  else
                    flag = .false.
                  end if
                end do
                ! write debug info
                if ( flag ) then
                  write(*,40) j, i, 'ice  ', 'water',  &
                    lm%ldmsk(j,i), lms%sfice(1,j,i)
                end if
              end if
            end if
          end if
        end if
      end do
    end do
    !
    !-----------------------------------------------------------------------
    ! Format definition
    !-----------------------------------------------------------------------
    !
! 20 format(' ATM land-sea mask is changed at (',I3,',',I3,') : ',     &
!             A5,' --> ',A5,' [',I2,']')
 30 format(' ATM sea-ice is formed at (',I3,',',I3,') : ',            &
             A5,' --> ',A5,' [',I2,' - ',F12.4,']')
 40 format(' ATM sea-ice is melted at (',I3,',',I3,') : ',       &
             A5,' --> ',A5,' [',I2,' - ',F12.4,']')
  end subroutine import_data_into_surface

  subroutine collect_output
    implicit none
#ifndef CLM
    integer(ik4) :: k
#endif
    integer(ik4) :: i , j , n
    real(rk8) :: qas , tas , ps , qs

    ! Fill accumulators

    if ( ktau > 0 ) then
      if ( ifatm ) then
        if ( associated(atm_tgb_out) ) &
          atm_tgb_out = atm_tgb_out + sum(lms%tgrd,1)*rdnnsg
        if ( associated(atm_tsw_out) ) &
          atm_tsw_out = atm_tsw_out + sum(lms%tsw,1)*rdnnsg
      end if
      if ( ifsrf ) then
        if ( associated(srf_evp_out) ) &
          srf_evp_out = srf_evp_out + lm%qfx
        if ( associated(srf_tpr_out) ) &
          srf_tpr_out = srf_tpr_out + sum(lms%prcp,1)*rdnnsg
        if ( associated(srf_prcv_out) ) &
          srf_prcv_out = srf_prcv_out + lm%cprate*rtsrf
        if ( associated(srf_zpbl_out) ) &
          srf_zpbl_out = srf_zpbl_out + lm%hpbl
        if ( associated(srf_scv_out) ) &
          srf_scv_out = srf_scv_out + sum(lms%sncv,1)*rdnnsg
        if ( associated(srf_sund_out) ) then
          where( lm%rswf > 120.0D0 )
            srf_sund_out = srf_sund_out + dtbat
          end where
        end if
        if ( associated(srf_srunoff_out) ) then
          srf_srunoff_out = srf_srunoff_out+sum(lms%srnof,1)*rdnnsg
        end if
        if ( associated(srf_trunoff_out) ) then
          srf_trunoff_out = srf_trunoff_out+sum(lms%trnof,1)*rdnnsg
        end if
        if ( associated(srf_sena_out) ) then
          srf_sena_out = srf_sena_out + sum(lms%sent,1)*rdnnsg
        end if
        if ( associated(srf_flw_out) ) &
          srf_flw_out = srf_flw_out + lm%rlwf
        if ( associated(srf_fsw_out) ) &
          srf_fsw_out = srf_fsw_out + lm%rswf
        if ( associated(srf_fld_out) ) &
          srf_fld_out = srf_fld_out + lm%dwrlwf
        if ( associated(srf_sina_out) ) &
          srf_sina_out = srf_sina_out + lm%solinc
        if ( associated(srf_snowmelt_out) ) &
          srf_snowmelt_out = srf_snowmelt_out + sum(lms%snwm,1)*rdnnsg
      end if
      if ( ifsub ) then
        call reorder_add_subgrid(lms%sfcp,sub_ps_out)
        if ( associated(sub_evp_out) ) &
          call reorder_add_subgrid(lms%evpr,sub_evp_out)
        if ( associated(sub_scv_out) ) &
          call reorder_add_subgrid(lms%sncv,sub_scv_out,mask=lm%ldmsk1)
        if ( associated(sub_sena_out) ) &
          call reorder_add_subgrid(lms%sent,sub_sena_out)
        if ( associated(sub_srunoff_out) ) &
          call reorder_add_subgrid(lms%srnof,sub_srunoff_out,lm%ldmsk1)
        if ( associated(sub_trunoff_out) ) &
          call reorder_add_subgrid(lms%trnof,sub_trunoff_out,lm%ldmsk1)
      end if
      if ( ifsts ) then
        if ( associated(sts_tgmax_out) ) &
          sts_tgmax_out = max(sts_tgmax_out,sum(lms%tgrd,1)*rdnnsg)
        if ( associated(sts_tgmin_out) ) &
          sts_tgmin_out = min(sts_tgmin_out,sum(lms%tgrd,1)*rdnnsg)
        if ( associated(sts_t2max_out) ) &
          sts_t2max_out(:,:,1) = max(sts_t2max_out(:,:,1),sum(lms%t2m,1)*rdnnsg)
        if ( associated(sts_t2min_out) ) &
          sts_t2min_out(:,:,1) = min(sts_t2min_out(:,:,1),sum(lms%t2m,1)*rdnnsg)
        if ( associated(sts_t2min_out) ) &
          sts_t2avg_out(:,:,1) = sts_t2avg_out(:,:,1) + sum(lms%t2m,1)*rdnnsg
        if ( associated(sts_w10max_out) ) &
          sts_w10max_out(:,:,1) = max(sts_w10max_out(:,:,1), &
            sqrt(sum((lms%u10m**2+lms%v10m**2),1)*rdnnsg))
        if ( associated(sts_pcpmax_out) ) &
          sts_pcpmax_out = max(sts_pcpmax_out,lms%prcp(1,:,:))
        if ( associated(sts_pcpavg_out) ) &
          sts_pcpavg_out = sts_pcpavg_out + lms%prcp(1,:,:)
        if ( associated(sts_psmin_out) ) &
          sts_psmin_out = min(sts_psmin_out,lm%sfps(jci1:jci2,ici1:ici2))
        if ( associated(sts_psavg_out) ) &
          sts_psavg_out = sts_psavg_out + lm%sfps(jci1:jci2,ici1:ici2)
        if ( associated(sts_srunoff_out) ) &
          sts_srunoff_out = sts_srunoff_out+sum(lms%srnof,1)*rdnnsg
        if ( associated(sts_trunoff_out) ) &
          sts_trunoff_out = sts_trunoff_out+sum(lms%trnof,1)*rdnnsg
        if ( associated(sts_sund_out) ) then
          where( lm%rswf > 120.0D0 )
            sts_sund_out = sts_sund_out + dtbat
          end where
        end if
      end if
      if ( iflak ) then
        if ( associated(lak_tpr_out) ) &
          lak_tpr_out = lak_tpr_out + sum(lms%prcp,1)*rdnnsg
        if ( associated(lak_scv_out) ) &
          lak_scv_out = lak_scv_out + sum(lms%sncv,1)*rdnnsg
        if ( associated(lak_sena_out) ) &
          lak_sena_out = lak_sena_out + sum(lms%sent,1)*rdnnsg
        if ( associated(lak_flw_out) ) &
          lak_flw_out = lak_flw_out + lm%rlwf
        if ( associated(lak_fsw_out) ) &
          lak_fsw_out = lak_fsw_out + lm%rswf
        if ( associated(lak_fld_out) ) &
          lak_fld_out = lak_fld_out + lm%dwrlwf
        if ( associated(lak_sina_out) ) &
          lak_sina_out = lak_sina_out + lm%solinc
        if ( associated(lak_evp_out) ) &
          lak_evp_out = lak_evp_out + sum(lms%evpr,1)*rdnnsg
      end if
    end if

    ! Those are for the output, but collected only at POINT in time

    if ( mod(ktau+1,ksrf) == 0 ) then

      if ( ifsrf ) then
        if ( associated(srf_uvdrag_out) ) &
          srf_uvdrag_out = sum(lms%drag,1)*rdnnsg
        if ( associated(srf_tg_out) ) &
          srf_tg_out = sum(lms%tgbb,1)*rdnnsg
        if ( associated(srf_tlef_out) ) then
          where ( lm%ldmsk > 0 )
            srf_tlef_out = sum(lms%tlef,1)*rdnnsg
          elsewhere
            srf_tlef_out = dmissval
          end where
        end if
        if ( associated(srf_aldirs_out) ) &
          srf_aldirs_out = lm%swdiralb
        if ( associated(srf_aldifs_out) ) &
          srf_aldifs_out = lm%swdifalb
        if ( associated(srf_seaice_out) ) &
          srf_seaice_out = sum(lms%sfice,1)*rdnnsg
        if ( associated(srf_t2m_out) ) &
          srf_t2m_out(:,:,1) = sum(lms%t2m,1)*rdnnsg
        if ( associated(srf_q2m_out) ) &
          srf_q2m_out(:,:,1) = sum(lms%q2m,1)*rdnnsg
        if ( associated(srf_rh2m_out) ) then
          srf_rh2m_out = d_zero
          do i = ici1 , ici2
            do j = jci1 , jci2
              do n = 1 , nnsg
                qas = lms%q2m(n,j,i)
                tas = lms%t2m(n,j,i)
                ps = lms%sfcp(n,j,i)
                qs = pfqsat(tas,ps)
                srf_rh2m_out(j,i,1) = srf_rh2m_out(j,i,1)+(qas/qs)*d_100
              end do
            end do
          end do
          srf_rh2m_out = srf_rh2m_out * rdnnsg
        end if
        if ( associated(srf_u10m_out) ) &
          srf_u10m_out(:,:,1) = sum(lms%u10m,1)*rdnnsg
        if ( associated(srf_v10m_out) ) &
          srf_v10m_out(:,:,1) = sum(lms%v10m,1)*rdnnsg
        if ( associated(srf_smw_out) ) then
          do n = 1 , num_soil_layers
            where ( lm%ldmsk > 0 )
              srf_smw_out(:,:,n) = sum(lms%sw(:,:,:,n),1)*rdnnsg
            elsewhere
              srf_smw_out(:,:,n) = dmissval
            end where
          end do
        end if
      end if

      if ( ifsub ) then
        if ( associated(sub_uvdrag_out) ) &
          call reorder_subgrid(lms%drag,sub_uvdrag_out)
        if ( associated(sub_tg_out) ) &
          call reorder_subgrid(lms%tgbb,sub_tg_out)
        if ( associated(sub_tlef_out) ) &
          call reorder_subgrid(lms%tlef,sub_tlef_out,mask=lm%ldmsk1)
        if ( associated(sub_u10m_out) ) &
          call reorder_subgrid(lms%u10m,sub_u10m_out)
        if ( associated(sub_v10m_out) ) &
          call reorder_subgrid(lms%v10m,sub_v10m_out)
        if ( associated(sub_t2m_out) ) &
          call reorder_subgrid(lms%t2m,sub_t2m_out)
        if ( associated(sub_q2m_out) ) &
          call reorder_subgrid(lms%q2m,sub_q2m_out)
        if ( associated(sub_smw_out) ) then
          call reorder_subgrid(lms%sw,sub_smw_out,lm%ldmsk1)
        end if
      end if

#ifndef CLM
      if ( iflak ) then
        if ( associated(lak_tg_out) ) &
          lak_tg_out = sum(lms%tgbb,1)*rdnnsg
        if ( associated(lak_aldirs_out) ) &
          lak_aldirs_out = lm%swdiralb
        if ( associated(lak_aldifs_out) ) &
          lak_aldifs_out = lm%swdifalb
        if ( associated(lak_ice_out) ) &
          lak_ice_out = sum(lms%sfice,1,lms%lakmsk)*rdnnsg
        if ( associated(lak_tlake_out) ) then
          do k = 1 , ndpmax
            lak_tlake_out(:,:,k) = sum(lms%tlake(:,:,:,k),1,lms%lakmsk)*rdnnsg
          end do
        end if
      end if
#endif

    end if ! IF output time

    if ( iocncpl == 1 ) then
      ! Fill for the RTM component
      lm%dailyrnf(:,:,1) = lm%dailyrnf(:,:,1) + sum(lms%srnof,1)*rdnnsg
      lm%dailyrnf(:,:,2) = lm%dailyrnf(:,:,2) + &
        (sum(lms%trnof,1)-sum(lms%srnof,1))*rdnnsg
      runoffcount = runoffcount + d_one
    end if

    ! Reset accumulation from precip and cumulus
    lm%ncprate = d_zero
    lm%cprate  = d_zero

    ! Reset also accumulation for deposition fluxes
    if ( ichem == 1 ) then
      lm%wetdepflx = d_zero
      lm%drydepflx = d_zero
    end if
  end subroutine collect_output

end module mod_lm_interface
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

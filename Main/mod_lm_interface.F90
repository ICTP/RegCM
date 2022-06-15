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
  use mod_dynparam
  use mod_runparams
  use mod_memutil
  use mod_regcm_types
  use mod_outvars
  use mod_constants
  use mod_mppparam
  use mod_mpmessage
  use mod_service
  use mod_bats_common
  use mod_ocn_common
  use mod_stdio
  use mod_slabocean
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
  real(rkx) :: runoffcount = 1.0_rkx
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

  real(rkx) , pointer , dimension(:,:) :: slp , sfp , slp1

  type(lm_exchange) :: lm
  type(lm_state) :: lms

  contains

  subroutine allocate_surface_model
    implicit none

    rdnnsg = d_one/real(nnsg,rkx)

    call getmem3d(lms%sent,1,nnsg,jci1,jci2,ici1,ici2,'lm:sent')
    call getmem3d(lms%evpr,1,nnsg,jci1,jci2,ici1,ici2,'lm:evpr')
    call getmem3d(lms%deltat,1,nnsg,jci1,jci2,ici1,ici2,'lm:deltat')
    call getmem3d(lms%deltaq,1,nnsg,jci1,jci2,ici1,ici2,'lm:deltaq')
    call getmem3d(lms%drag,1,nnsg,jci1,jci2,ici1,ici2,'lm:drag')
    call getmem3d(lms%ustar,1,nnsg,jci1,jci2,ici1,ici2,'lm:ustar')
    call getmem3d(lms%w10m,1,nnsg,jci1,jci2,ici1,ici2,'lm:w10m')
    call getmem3d(lms%zo,1,nnsg,jci1,jci2,ici1,ici2,'lm:zo')
    call getmem3d(lms%rhoa,1,nnsg,jci1,jci2,ici1,ici2,'lm:rho')
    call getmem3d(lms%lncl,1,nnsg,jci1,jci2,ici1,ici2,'lm:lncl')
    call getmem3d(lms%prcp,1,nnsg,jci1,jci2,ici1,ici2,'lm:prcp')
    call getmem3d(lms%snwm,1,nnsg,jci1,jci2,ici1,ici2,'lm:snwm')
    call getmem3d(lms%trnof,1,nnsg,jci1,jci2,ici1,ici2,'lm:trnof')
    call getmem3d(lms%srnof,1,nnsg,jci1,jci2,ici1,ici2,'lm:srnof')
    call getmem3d(lms%xlai,1,nnsg,jci1,jci2,ici1,ici2,'lm:xlai')
    call getmem3d(lms%sfcp,1,nnsg,jci1,jci2,ici1,ici2,'lm:sfcp')
    call getmem3d(lms%q2m,1,nnsg,jci1,jci2,ici1,ici2,'lm:q2m')
    call getmem3d(lms%t2m,1,nnsg,jci1,jci2,ici1,ici2,'lm:t2m')
    call getmem3d(lms%u10m,1,nnsg,jci1,jci2,ici1,ici2,'lm:u10m')
    call getmem3d(lms%v10m,1,nnsg,jci1,jci2,ici1,ici2,'lm:v10m')
    call getmem3d(lms%ram1,1,nnsg,jci1,jci2,ici1,ici2,'lm:ram1')
    call getmem3d(lms%rah1,1,nnsg,jci1,jci2,ici1,ici2,'lm:rah1')
    call getmem3d(lms%br,1,nnsg,jci1,jci2,ici1,ici2,'lm:br')
    call getmem3d(lms%taux,1,nnsg,jci1,jci2,ici1,ici2,'lm:taux')
    call getmem3d(lms%tauy,1,nnsg,jci1,jci2,ici1,ici2,'lm:tauy')
    call getmem3d(lms%wt,1,nnsg,jci1,jci2,ici1,ici2,'lm:wt')
    call getmem3d(lms%swalb,1,nnsg,jci1,jci2,ici1,ici2,'lm:swalb')
    call getmem3d(lms%lwalb,1,nnsg,jci1,jci2,ici1,ici2,'lm:lwalb')
    call getmem3d(lms%swdiralb,1,nnsg,jci1,jci2,ici1,ici2,'lm:swdiralb')
    call getmem3d(lms%lwdiralb,1,nnsg,jci1,jci2,ici1,ici2,'lm:lwdiralb')
    call getmem3d(lms%swdifalb,1,nnsg,jci1,jci2,ici1,ici2,'lm:swdifalb')
    call getmem3d(lms%lwdifalb,1,nnsg,jci1,jci2,ici1,ici2,'lm:lwdifalb')

    call getmem3d(lms%gwet,1,nnsg,jci1,jci2,ici1,ici2,'lm:gwet')
    call getmem4d(lms%sw,1,nnsg,jci1,jci2,ici1,ici2,1,num_soil_layers,'lm:sw')

#ifndef CLM45
    call assignpnt(lms%sw,lms%ssw,1)
    call assignpnt(lms%sw,lms%rsw,2)
    call assignpnt(lms%sw,lms%tsw,3)
#else
    call assignpnt(lms%sw,lms%ssw,1)
    call assignpnt(lms%sw,lms%rsw,2)
    call getmem3d(lms%tsw,1,nnsg,jci1,jci2,ici1,ici2,'lm:tsw')
#endif
    call getmem3d(lms%tgbb,1,nnsg,jci1,jci2,ici1,ici2,'lm:tgbb')
    call getmem3d(lms%tgrd,1,nnsg,jci1,jci2,ici1,ici2,'lm:tgrd')
    call getmem3d(lms%tgbrd,1,nnsg,jci1,jci2,ici1,ici2,'lm:tgbrd')
    call getmem3d(lms%tlef,1,nnsg,jci1,jci2,ici1,ici2,'lm:tlef')
    call getmem3d(lms%taf,1,nnsg,jci1,jci2,ici1,ici2,'lm:taf')
    call getmem3d(lms%sigf,1,nnsg,jci1,jci2,ici1,ici2,'lm:sigf')
    call getmem3d(lms%sfice,1,nnsg,jci1,jci2,ici1,ici2,'lm:sfice')
    call getmem3d(lms%snag,1,nnsg,jci1,jci2,ici1,ici2,'lm:snag')
    call getmem3d(lms%ldew,1,nnsg,jci1,jci2,ici1,ici2,'lm:ldew')
    call getmem3d(lms%sncv,1,nnsg,jci1,jci2,ici1,ici2,'lm:sncv')
    call getmem3d(lms%scvk,1,nnsg,jci1,jci2,ici1,ici2,'lm:scvk')
    call getmem3d(lms%um10,1,nnsg,jci1,jci2,ici1,ici2,'lm:um10')
    call getmem3d(lms%emisv,1,nnsg,jci1,jci2,ici1,ici2,'lm:emisv')
#ifdef CLM45
    call getmem4d(lms%vocemiss,1,nnsg,jci1,jci2,ici1,ici2,1,ntr,'lm:vocemiss')
    call getmem4d(lms%dustemiss,1,nnsg,jci1,jci2,ici1,ici2,1,4,'lm:dustemiss')
    call getmem4d(lms%sw_vol,1,nnsg,jci1,jci2, &
                                    ici1,ici2,1,num_soil_layers,'lm:sw_vol')
    call getmem4d(lms%tsoi,1,nnsg,jci1,jci2, &
                                  ici1,ici2,1,num_soil_layers,'lm:tsoi')
#endif

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
    call getmem2d(sfp,jce1ga,jce2ga,ice1ga,ice2ga,'lm:sfp')
    call getmem2d(slp,jce1ga,jce2ga,ice1ga,ice2ga,'lm:slp')
    call getmem2d(slp1,jce1ga,jce2ga,ice1ga,ice2ga,'lm:slp1')
  end subroutine allocate_surface_model

  subroutine init_surface_model
    use mod_atm_interface
    use mod_che_interface
    implicit none

    call cl_setup(lndcomm,mddom%mask,mdsub%mask)
    call cl_setup(ocncomm,mddom%mask,mdsub%mask,.true.)

#ifdef DEBUG
    write(ndebug,*) 'TOTAL POINTS FOR LAND  IN LNDCOMM : ', &
      lndcomm%linear_npoint_sg(myid+1)
    write(ndebug,*) 'Cartesian p ', lndcomm%cartesian_npoint_g
    write(ndebug,*) 'Cartesian d ', lndcomm%cartesian_displ_g
    write(ndebug,*) 'Linear    p ', lndcomm%linear_npoint_g
    write(ndebug,*) 'Linear    d ', lndcomm%linear_displ_g
    write(ndebug,*) 'Subgrid Cartesian p ', lndcomm%cartesian_npoint_sg
    write(ndebug,*) 'Subgrid Cartesian d ', lndcomm%cartesian_displ_sg
    write(ndebug,*) 'Subgrid Linear    p ', lndcomm%linear_npoint_sg
    write(ndebug,*) 'Subgrid Linear    d ', lndcomm%linear_displ_sg
    write(ndebug,*) 'TOTAL POINTS FOR OCEAN IN OCNCOMM : ', &
      ocncomm%linear_npoint_sg(myid+1)
    write(ndebug,*) 'Cartesian p ', ocncomm%cartesian_npoint_g
    write(ndebug,*) 'Cartesian d ', ocncomm%cartesian_displ_g
    write(ndebug,*) 'Linear    p ', ocncomm%linear_npoint_g
    write(ndebug,*) 'Linear    d ', ocncomm%linear_displ_g
    write(ndebug,*) 'Subgrid Cartesian p ', ocncomm%cartesian_npoint_sg
    write(ndebug,*) 'Subgrid Cartesian d ', ocncomm%cartesian_displ_sg
    write(ndebug,*) 'Subgrid Linear    p ', ocncomm%linear_npoint_sg
    write(ndebug,*) 'Subgrid Linear    d ', ocncomm%linear_displ_sg
#endif

    call allocate_mod_bats_internal(lndcomm)
    call allocate_mod_ocn_internal(ocncomm)

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
    call assignpnt(mddom%itex,lm%itex)
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
    call assignpnt(mdsub%itex,lm%itex1)
    call assignpnt(mdsub%dhlake,lm%dhlake1)
    call assignpnt(cplmsk,lm%icplmsk)
    call assignpnt(atms%ubx3d,lm%uatm,kz)
    call assignpnt(atms%vbx3d,lm%vatm,kz)
    call assignpnt(atms%th3d,lm%thatm,kz)
    call assignpnt(atms%tb3d,lm%tatm,kz)
    call assignpnt(atms%pb3d,lm%patm,kz)
    call assignpnt(atms%qxb3d,lm%qvatm,kz,iqv)
    call assignpnt(atms%za,lm%hgt,kz)
    call assignpnt(atms%rhox2d,lm%rhox)
    call assignpnt(atms%ps2d,lm%sfps)
    call assignpnt(atms%tp2d,lm%sfta)
    call assignpnt(sfs%hfx,lm%hfx)
    call assignpnt(sfs%qfx,lm%qfx)
    call assignpnt(sfs%uvdrag,lm%uvdrag)
    call assignpnt(sfs%ustar,lm%ustar)
    call assignpnt(sfs%w10m,lm%w10m)
    call assignpnt(sfs%zo,lm%zo)
    call assignpnt(sfs%tg,lm%tg)
    call assignpnt(sfs%tgbb,lm%tgbb)
    call assignpnt(sfs%u10m,lm%u10m)
    call assignpnt(sfs%v10m,lm%v10m)
    call assignpnt(sfs%ram1,lm%ram1)
    call assignpnt(sfs%rah1,lm%rah1)
    call assignpnt(sfs%br,lm%br)
    call assignpnt(sfs%q2m,lm%q2m)
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
    call assignpnt(totcf,lm%totcf)
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
#ifdef CLM45
      call assignpnt(sw_vol,lm%sw_vol)
      call assignpnt(tsoi,lm%tsoi)
      call assignpnt(idust,lm%idust)
#endif
    end if
    if ( iocncpl == 1 .or. iwavcpl == 1) then
      call assignpnt(dailyrnf,lm%dailyrnf)
    end if
#ifdef CLM
    allocate(landmask(jx,iy))
#endif
  end subroutine init_surface_model

  subroutine initialize_surface_model
    implicit none
#ifdef CLM
    integer(ik4) :: i , j , n
#endif
    if ( irceideal == 0 ) then
#ifndef CLM45
      call initbats(lm,lms)
#else
      call initclm45(lm,lms)
#endif
    end if
    call initocn(lm,lms)
#ifdef CLM
    if ( irceideal == 0 ) then
      call initclm(lm,lms)
      if ( rcmtimer%start( ) .and. imask == 2 ) then
        ! CLM may have changed the landuse again !
        do i = ici1 , ici2
          do j = jci1 , jci2
            lm%iveg(j,i) = nint(lm%lndcat(j,i))
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
    end if
#endif
    lm%emissivity = sum(lms%emisv,1) * rdnnsg
  end subroutine initialize_surface_model

  subroutine surface_model
#ifdef CLM45
    use mod_atm_interface , only : voc_em_clm , dustflx_clm
#endif
    implicit none
    integer(ik4) :: i , j , n , nn , ierr
#ifdef CLM
    if ( rcmtimer%start( ) .or. syncro_rad%will_act(dtsrf) ) then
      r2cdoalb = .true.
    else
      r2cdoalb = .false.
    end if
    ! Timestep used is the same as for bats
    if ( rcmtimer%start( ) ) then
      r2cnstep = 0
    else
      r2cnstep = syncro_srf%lcount + 1
    end if
    if ( irceideal == 0 ) call mtrxclm(lm,lms)
#else
#ifdef CLM45
    if ( irceideal == 0 ) call runclm45(lm,lms)
    !coupling of biogenic VOC from CLM45 to chemistry
    if ( ichem == 1 ) then
      do n = 1 , ntr
        do i = ici1 , ici2
          do j = jci1 , jci2
            voc_em_clm(j,i,n) = sum(lms%vocemiss(:,j,i,n),1) * rdnnsg
          end do
        end do
      end do
      do n = 1 , 4
        do i = ici1 , ici2
          do j = jci1 , jci2
            dustflx_clm(j,i,n) = sum(lms%dustemiss(:,j,i,n),1) * rdnnsg
          end do
        end do
      end do
    end if
#else
    if ( irceideal == 0 ) call vecbats(lm,lms)
#endif
#endif
!FAB  
    if ( islab_ocean == 1 ) call update_slabocean(xslabtime,lms)

    call vecocn(lm,lms)
    ! Fill land part of this output vars
    do i = ici1, ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          lms%w10m(n,j,i)  = sqrt(lms%u10m(n,j,i)**2 + lms%v10m(n,j,i)**2)
          if ( lm%ldmsk1(n,j,i) == 1 ) then
            lms%rhoa(n,j,i) = lms%sfcp(n,j,i)/(rgas*lms%t2m(n,j,i))
            lms%ustar(n,j,i) = sqrt(sqrt( &
                                  (lms%u10m(n,j,i)*lms%drag(n,j,i))**2 + &
                                  (lms%v10m(n,j,i)*lms%drag(n,j,i))**2) / &
                                  lms%rhoa(n,j,i))
          end if
        end do
      end do
    end do
#ifndef CLM45
    do i = ici1, ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          if ( lm%ldmsk1(n,j,i) == 1 ) then
            lms%taux(n,j,i) = lms%drag(n,j,i) * (lms%u10m(n,j,i)/lm%uatm(j,i))
            lms%tauy(n,j,i) = lms%drag(n,j,i) * (lms%v10m(n,j,i)/lm%vatm(j,i))
          end if
        end do
      end do
    end do
#endif
    lm%hfx = sum(lms%sent,1)*rdnnsg
    lm%qfx = sum(lms%evpr,1)*rdnnsg
    lm%uvdrag = sum(lms%drag,1)*rdnnsg
    lm%ustar = sum(lms%ustar,1)*rdnnsg
    lm%u10m = sum(lms%u10m,1)*rdnnsg
    lm%v10m = sum(lms%v10m,1)*rdnnsg
    lm%ram1 = sum(lms%ram1,1)*rdnnsg
    lm%rah1 = sum(lms%rah1,1)*rdnnsg
    lm%br = sum(lms%br,1)*rdnnsg
    lm%q2m = sum(lms%q2m,1)*rdnnsg
    lm%w10m = sum(lms%w10m,1)*rdnnsg
    lm%zo = sum(lms%zo,1)*rdnnsg
    lm%tgbb = sum(lms%tgbb,1)*rdnnsg
    lm%tg = sum(lms%tgrd,1)*rdnnsg
    lm%emissivity = sum(lms%emisv,1) * rdnnsg
    if ( ichem == 1 ) then
      lm%deltat = sum(lms%deltat,1)*rdnnsg
      lm%deltaq = sum(lms%deltaq,1)*rdnnsg
      lm%sxlai2d = sum(lms%xlai,1)*rdnnsg
#ifdef CLM45
      ! FAB here take humidity of first soil layer, sw should be always defined
      lm%ssw2da = sum(lms%tsw(:,:,:),1)*rdnnsg
      lm%sw_vol = sum(lms%sw_vol(:,:,:,:),1)*rdnnsg
      lm%tsoi = sum(lms%tsoi(:,:,:,:),1)*rdnnsg
      lm%sfracb2d = sum(lms%wt,1)*rdnnsg
#else
      ! FAB here take humidity of first soil layer, sw should be always defined
      lm%ssw2da = sum(lms%ssw(:,:,:),1)*rdnnsg
      lm%sfracv2d = sum(lms%sigf,1)*rdnnsg
      lm%svegfrac2d = sum(lms%lncl,1)*rdnnsg
      lm%sfracs2d = sum((lms%lncl*lms%wt+(d_one-lms%lncl)*lms%scvk),1)*rdnnsg
      lm%sfracb2d = sum(((d_one-lms%lncl)*(d_one-lms%scvk)),1)*rdnnsg
#endif
    end if
    call collect_output
#ifdef DEBUG
    ! Sanity check of surface temperatures
    ierr = 0
    do i = ici1 , ici2
      do j = jci1 , jci2
        do n = 1 , nnsg
          if ( lms%tgrd(n,j,i) < 150.0_rkx ) then
            write(stderr,*) 'Likely error: Surface temperature too low'
            write(stderr,*) 'MYID = ', myid
            write(stderr,*) 'J    = ',j
            write(stderr,*) 'I    = ',i
            do nn = 1 , nnsg
              write(stderr,*) 'N    = ',nn
              write(stderr,*) 'VAL  = ',lms%tgrd(nn,j,i)
              write(stderr,*) 'MASK = ',lm%ldmsk1(n,j,i)
            end do
            write(stderr,*) 'VAL2 = ',lm%tgbb(j,i)
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
      if ( irceideal == 0 ) call albedobats(lm,lms)
    else
      if ( irceideal == 0 ) call albedoclm(lm,lms)
    end if
#else
#ifdef CLM45
    if ( irceideal == 0 ) call albedoclm45(lm,lms)
#else
    if ( irceideal == 0 ) call albedobats(lm,lms)
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

    if ( .not. associated(expfie%psfc) ) then
      call fatal(__FILE__,__LINE__, &
        'RUNNING COUPLED WITHOUT COUPLER INITIALIZATION')
    end if
    do i = ici1 , ici2
      do j = jci1 , jci2
        expfie%psfc(j,i) = lm%sfps(j,i)*d_r100
        expfie%tsfc(j,i) = sum(lms%t2m(:,j,i))*rdnnsg
        expfie%qsfc(j,i) = lm%q2m(j,i)
        expfie%swrd(j,i) = lm%rswf(j,i)
        expfie%lwrd(j,i) = lm%rlwf(j,i)
        expfie%dlwr(j,i) = lm%dwrlwf(j,i)
        expfie%dswr(j,i) = lm%swdif(j,i)+lm%swdir(j,i)
        expfie%lhfx(j,i) = sum(lms%evpr(:,j,i)*wlh(lms%tgrd(:,j,i)))*rdnnsg
        expfie%shfx(j,i) = sum(lms%sent(:,j,i))*rdnnsg
        expfie%prec(j,i) = sum(lms%prcp(:,j,i))*rdnnsg
        expfie%wndu(j,i) = lm%u10m(j,i)
        expfie%wndv(j,i) = lm%v10m(j,i)
        expfie%taux(j,i) = sum(lms%taux(:,j,i))*rdnnsg
        expfie%tauy(j,i) = sum(lms%tauy(:,j,i))*rdnnsg
        expfie%sflx(j,i) = (sum(lms%evpr(:,j,i))-sum(lms%prcp(:,j,i)))*rdnnsg
        expfie%snow(j,i) = sum(lms%sncv(:,j,i))*rdnnsg
        expfie%wspd(j,i) = sqrt(expfie%wndu(j,i)**2+expfie%wndv(j,i)**2)
        expfie%wdir(j,i) = atan2(expfie%wndu(j,i), expfie%wndv(j,i))
        if (expfie%wdir(j,i) < d_zero) then
          expfie%wdir(j,i) = expfie%wdir(j,i)+twopi
        end if
        expfie%ustr(j,i) = sum(lms%ustar(:,j,i))*rdnnsg
        expfie%nflx(j,i) = lm%rswf(j,i) - expfie%lhfx(j,i) - &
                           expfie%shfx(j,i) - lm%rlwf(j,i)
      end do
    end do
    if ( rcmtimer%lcount == 1 .or. alarm_day%will_act(dtsec) ) then
      do i = ici1 , ici2
        do j = jci1 , jci2
          if ( lm%ldmsk(j,i) == 1 ) then
            expfie%rnof(j,i) = lm%dailyrnf(j,i,1)/runoffcount
            expfie%snof(j,i) = lm%dailyrnf(j,i,2)/runoffcount
          else
            expfie%rnof(j,i) = d_zero
            expfie%snof(j,i) = d_zero
          end if
        end do
      end do
      runoffcount = d_one
      lm%dailyrnf(:,:,:) = d_zero
    end if

    contains

#include <wlh.inc>

  end subroutine export_data_from_surface
!
  subroutine import_data_into_surface(impfie,ldmskb,wetdry,tol)
    use mod_atm_interface
    implicit none
    type(imp_data) , intent(in) :: impfie
    real(rkx) , intent(in) :: tol
    integer(ik4) , pointer , dimension(:,:) , intent(in) :: ldmskb , wetdry
    integer :: i , j , n
    logical :: flag = .false.
!    character (len=*) , parameter :: f99001 =                   &
!       "(' ATM land-sea mask is changed at (',I3,',',I3,') : ', &
!          &A5,' --> ',A5,' [',I2,']')"
    character (len=*) , parameter :: f99002 =            &
       "(' ATM sea-ice is formed at (',I3,',',I3,') : ', &
          &A5,' --> ',A5,' [',I2,' - ',F12.4,']')"
    character (len=*) , parameter :: f99003 =            &
       "(' ATM sea-ice is melted at (',I3,',',I3,') : ', &
          &A5,' --> ',A5,' [',I2,' - ',F12.4,']')"
    ! real(rkx) :: toth
    ! real(rkx) , parameter :: href = d_two * iceminh
    ! real(rkx) , parameter :: steepf = 1.0_rkx
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
            if ( syncro_cpl%lcount == 1 ) then
              cplmsk(j,i) = 1
            end if
            lm%tg(j,i)       = impfie%sst(j,i)
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
!               write(*,f99001) j, i, 'water', 'land ', lm%ldmsk(j,i)
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
!                 write(*,f99001) j, i, 'land ', 'water', lm%ldmsk(j,i)
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
                write(*,f99002) j, i, 'water', 'ice  ', &
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
                  write(*,f99003) j, i, 'ice  ', 'water',  &
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
    ! Retrieve information from WAV component
    !-----------------------------------------------------------------------
    !
    do i = ici1, ici2
      do j = jci1, jci2
        if ( lm%iveg(j,i) == 14 .or. lm%iveg(j,i) == 15 ) then
          !
          !--------------------------------------
          ! Update: Surface roughness
          !--------------------------------------
          !
          if ( impfie%zo(j,i) < tol ) then
            lm%zo(j,i) = impfie%zo(j,i)
          else
            lm%zo(j,i) = 1.0e20
          end if
          !
          !--------------------------------------
          ! Update: Friction velocity
          !--------------------------------------
          !
          if ( impfie%ustar(j,i) < tol ) then
            lm%ustar(j,i) = impfie%ustar(j,i)
          else
            lm%ustar(j,i) = 1.0e20
          end if
        end if
      end do
    end do
  end subroutine import_data_into_surface

  subroutine collect_output
    implicit none
#ifndef CLM
    integer(ik4) :: k
#endif
    integer(ik4) :: i , j , n
    real(rkx) :: qas , tas , uas , ps , qs , es , desdt , rh , sws , lws

    ! Fill accumulators

    if ( rcmtimer%integrating( ) ) then
      if ( ifrad ) then
        rnrad_for_srffrq = rnrad_for_srffrq + 1.0_rkx
      end if
      if ( ifatm ) then
        rnsrf_for_atmfrq = rnsrf_for_atmfrq + 1.0_rkx
        if ( associated(atm_tsw_out) ) &
          atm_tsw_out = atm_tsw_out + sum(lms%tsw,1)*rdnnsg
      end if
      if ( ifsrf ) then
        rnsrf_for_srffrq = rnsrf_for_srffrq + 1.0_rkx
        if ( associated(srf_totcf_out) ) &
          srf_totcf_out = srf_totcf_out + lm%totcf
        if ( associated(srf_evp_out) ) &
          srf_evp_out = srf_evp_out + lm%qfx
        if ( associated(srf_tpr_out) ) &
          srf_tpr_out = srf_tpr_out + sum(lms%prcp,1)*rdnnsg
        if ( associated(srf_prcv_out) ) &
          srf_prcv_out = srf_prcv_out + lm%cprate*syncro_srf%rw
        if ( associated(srf_zpbl_out) ) &
          srf_zpbl_out = srf_zpbl_out + lm%hpbl
        if ( associated(srf_scv_out) ) &
          srf_scv_out = srf_scv_out + sum(lms%sncv,1)*rdnnsg
        if ( associated(srf_sund_out) ) then
          where( lm%rswf > 120.0_rkx )
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
        if ( associated(srf_lena_out) ) then
          srf_lena_out = srf_lena_out + sum(lms%evpr*wlh(lms%tgrd),1)*rdnnsg
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
        if ( associated(srf_taux_out) ) &
          srf_taux_out = srf_taux_out + sum(lms%taux,1)*rdnnsg
        if ( associated(srf_tauy_out) ) &
          srf_tauy_out = srf_tauy_out + sum(lms%tauy,1)*rdnnsg
        if ( associated(srf_evpot_out) ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              if ( lm%ldmsk(j,i) == 0 ) then
                srf_evpot_out(j,i) = srf_evpot_out(j,i) + lm%qfx(j,i)
                cycle
              end if
              tas = sum(lms%t2m(:,j,i))*rdnnsg
              ps = sum(lms%sfcp(:,j,i))*rdnnsg
              es = pfesat(tas)
              qs = pfwsat(tas,ps,es)
              qas = lm%q2m(j,i)
              uas = lm%w10m(j,i)
              rh = min(max((qas/qs),d_zero),d_one)
              desdt = pfesdt(tas)
              sws = lm%rswf(j,i)
              lws = lm%rlwf(j,i)
              srf_evpot_out(j,i) = srf_evpot_out(j,i) + &
                   evpt_fao(ps,tas,uas,rh*es,es,desdt,sws,lws)
            end do
          end do
        end if
      end if
      if ( ifsub ) then
        rnsrf_for_subfrq = rnsrf_for_subfrq + 1.0_rkx
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
      if ( ifshf ) then
        if ( associated(shf_pcpmax_out) ) &
          shf_pcpmax_out = max(shf_pcpmax_out,sum(lms%prcp,1)*rdnnsg)
        if ( associated(shf_pcpavg_out) ) &
          shf_pcpavg_out = shf_pcpavg_out + sum(lms%prcp,1)*rdnnsg
        if ( associated(shf_pcprcv_out) ) &
          shf_pcprcv_out = shf_pcprcv_out + lm%cprate*syncro_srf%rw
        if ( associated(shf_twetb_out) ) then
          do i = ici1 , ici2
            do j = jci1 , jci2
              tas = sum(lms%t2m(:,j,i))*rdnnsg
              ps = sum(lms%sfcp(:,j,i))*rdnnsg
              qs = pfwsat(tas,ps)
              qas = lm%q2m(j,i)
              rh = min(max((qas/qs),d_zero),d_one)*100.0_rkx
              shf_twetb_out(j,i) = max(shf_twetb_out(j,i), &
                tas * atan(0.151977_rkx * sqrt(rh + 8.313659_rkx)) + &
                      atan(tas + rh) - atan(rh - 1.676331_rkx) + &
                      0.00391838_rkx * rh**(3.0_rkx/2.0_rkx) * &
                      atan(0.023101_rkx * rh) - 4.686035_rkx)
            end do
          end do
        end if
      end if
      if ( ifsts ) then
        rnsrf_for_day = rnsrf_for_day + 1.0_rkx
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
            sqrt(lm%u10m**2+lm%v10m**2))
        if ( associated(sts_pcpmax_out) ) &
          sts_pcpmax_out = max(sts_pcpmax_out,sum(lms%prcp,1)*rdnnsg)
        if ( associated(sts_pcpavg_out) ) &
          sts_pcpavg_out = sts_pcpavg_out + sum(lms%prcp,1)*rdnnsg
        if ( associated(sts_psmin_out) ) &
          sts_psmin_out = min(sts_psmin_out,lm%sfps(jci1:jci2,ici1:ici2))
        if ( associated(sts_psavg_out) ) &
          sts_psavg_out = sts_psavg_out + lm%sfps(jci1:jci2,ici1:ici2)
        if ( associated(sts_srunoff_out) ) &
          sts_srunoff_out = sts_srunoff_out+sum(lms%srnof,1)*rdnnsg
        if ( associated(sts_trunoff_out) ) &
          sts_trunoff_out = sts_trunoff_out+sum(lms%trnof,1)*rdnnsg
        if ( associated(sts_sund_out) ) then
          where( lm%rswf > 120.0_rkx )
            sts_sund_out = sts_sund_out + dtbat
          end where
        end if
      end if
      if ( iflak ) then
        rnsrf_for_lakfrq = rnsrf_for_lakfrq + 1.0_rkx
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
    if ( alarm_out_atm%will_act(dtsrf) ) then
      if ( ifatm ) then
        if ( associated(atm_tgb_out) ) then
          atm_tgb_out = sum(lms%tgbb,1)*rdnnsg
        end if
      end if
    end if

    if ( alarm_out_srf%will_act(dtsrf) ) then

      if ( ifsrf ) then
        if ( associated(srf_uvdrag_out) ) &
          srf_uvdrag_out = sum(lms%drag,1)*rdnnsg
        if ( associated(srf_ustar_out) ) &
          srf_ustar_out = sum(lms%ustar,1)*rdnnsg
        if ( associated(srf_zo_out) ) &
          srf_zo_out = sum(lms%zo,1)*rdnnsg
        if ( associated(srf_rhoa_out) ) &
          srf_rhoa_out = sum(lms%rhoa,1)*rdnnsg
        if ( associated(srf_tg_out) ) &
          srf_tg_out = sum(lms%tgrd,1)*rdnnsg
        if ( associated(srf_tlef_out) ) then
          where ( lm%ldmsk == 1 )
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
          srf_q2m_out(:,:,1) = lm%q2m
        if ( associated(srf_rh2m_out) ) then
          srf_rh2m_out = d_zero
          do i = ici1 , ici2
            do j = jci1 , jci2
              do n = 1 , nnsg
                qas = lms%q2m(n,j,i)
                tas = lms%t2m(n,j,i)
                ps = lms%sfcp(n,j,i)
                qs = pfwsat(tas,ps)
                srf_rh2m_out(j,i,1) = srf_rh2m_out(j,i,1) + &
                              min(max((qas/qs),rhmin),rhmax)*d_100
              end do
            end do
          end do
          srf_rh2m_out = srf_rh2m_out * rdnnsg
        end if
        if ( associated(srf_u10m_out) ) &
          srf_u10m_out(:,:,1) = lm%u10m
        if ( associated(srf_v10m_out) ) &
          srf_v10m_out(:,:,1) = lm%v10m
        if ( associated(srf_smw_out) ) then
          do n = 1 , num_soil_layers
            where ( lm%ldmsk == 1 )
              srf_smw_out(:,:,n) = sum(lms%sw(:,:,:,n),1)*rdnnsg
            elsewhere
              srf_smw_out(:,:,n) = dmissval
            end where
          end do
        end if
#ifdef CLM45
        if ( associated(srf_tsoil_out) ) then
          do n = 1 , num_soil_layers
            where ( lm%ldmsk == 1 )
              srf_tsoil_out(:,:,n) = sum(lms%tsoi(:,:,:,n),1)*rdnnsg
            elsewhere
              srf_tsoil_out(:,:,n) = dmissval
            end where
          end do
        end if
#endif
        if ( associated(srf_mslp_out) ) then
          call mslp
          srf_mslp_out(jci1:jci2,ici1:ici2) = slp(jci1:jci2,ici1:ici2)
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
            do i = ici1 , ici2
              do j = jci1 , jci2
                if ( lm%iveg(j,i) == 14 ) then
                  lak_tlake_out(j,i,k) = tzero + sum(lms%tlake(:,j,i,k))*rdnnsg
                end if
              end do
            end do
          end do
        end if
      end if
#endif

    end if ! IF output time

    if ( iocncpl == 1 .or. iwavcpl == 1 ) then
      ! Fill for the RTM component
      lm%dailyrnf(:,:,1) = lm%dailyrnf(:,:,1) + sum(lms%srnof,1)*rdnnsg
      lm%dailyrnf(:,:,2) = lm%dailyrnf(:,:,2) + &
        (sum(lms%trnof,1)-sum(lms%srnof,1))*rdnnsg
      runoffcount = runoffcount + d_one
    end if

    ! Reset also accumulation for deposition fluxes
    if ( ichem == 1 ) then
      lm%wetdepflx = d_zero
      lm%drydepflx = d_zero
    end if

    ! Reset accumulation from precip and cumulus
    lm%ncprate(:,:) = d_zero
    lm%cprate(:,:)  = d_zero

    contains

#include <pfesat.inc>
#include <pfwsat.inc>
#include <pfdesatdt.inc>
#include <pqderiv.inc>
#include <wlh.inc>
#include <evpt.inc>

    subroutine mslp
      implicit none
      integer(ik4) :: i , j , n
      real(rkx) :: tstar , hstar , alpha , raval , mval , mall
      integer(ik4) , parameter :: niter = 20
      real(rkx) , dimension(jci1:jci2,ici1:ici2) :: mask

      ! Follow Kallen 1996
      alpha = lrate*rgas/egrav
      do i = ice1 , ice2
        do j = jce1 , jce2
          tstar = lm%tatm(j,i)
          if ( tstar < 255.0_rkx ) then
            tstar = (tstar+255.0_rkx)*0.5_rkx
          else if ( tstar > 290.5_rkx ) then
            tstar = 290.5_rkx + (0.005_rkx*(tstar-290.5_rkx))**2
          end if
          hstar = lm%ht(j,i)/(rgas*tstar)
          raval = d_half*alpha*hstar
          slp(j,i) = lm%sfps(j,i) * &
               exp(hstar*(1.0_rkx - raval + (raval*raval)/3.0_rkx))
        end do
      end do
      ! Gauss Siedel Filtering
      mval = d_half*(maxval(lm%sfps)-minval(lm%sfps))
      call sumall(mval,mall)
      mval = mall/real(nproc,rkx)
      sfp(jce1:jce2,ice1:ice2) = lm%sfps(jce1:jce2,ice1:ice2)
      call exchange(slp,1,jce1,jce2,ice1,ice2)
      call exchange(sfp,1,jce1,jce2,ice1,ice2)
      slp1 = slp
      mask = d_zero
      do i = ici1 , ici2
        do j = jci1 , jci2
          mask(j,i) = (sfp(j,i-1)+sfp(j,i+1) + &
                       sfp(j-1,i)+sfp(j+1,i) - &
                       4.0_rkx*sfp(j,i))/mval
        end do
      end do
      do n = 1 , niter
        do i = ici1 , ici2
          do j = jci1 , jci2
            slp1(j,i) = d_rfour*(slp1(j,i-1)+slp(j,i+1) + &
                                 slp1(j-1,i)+slp(j+1,i)-mask(j,i))
          end do
        end do
        if ( ma%has_bdyleft ) then
          do i = ici1 , ici2
            slp1(jce1,i) = slp1(jci1,i)
          end do
        end if
        if ( ma%has_bdyright ) then
          do i = ici1 , ici2
            slp1(jce2,i) = slp1(jci2,i)
          end do
        end if
        if ( ma%has_bdybottom ) then
          do j = jce1 , jce2
            slp1(j,ice1) = slp1(j,ici1)
          end do
        end if
        if ( ma%has_bdytop ) then
          do j = jce1 , jce2
            slp1(j,ice2) = slp1(j,ici2)
          end do
        end if
        call exchange(slp1,1,jce1,jce2,ice1,ice2)
        slp(:,:) = slp1
      end do
    end subroutine mslp

  end subroutine collect_output

end module mod_lm_interface
! vim: tabstop=8 expandtab shiftwidth=2 softtabstop=2

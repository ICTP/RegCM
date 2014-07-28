module mod_clm_pftvarcon
  use mod_intkinds
  use mod_realkinds
  use mod_mpmessage
  use mod_clm_nchelper
  use mod_dynparam
  use mod_mppparam
  use mod_clm_varpar , only : mxpft , numrad , ivis , inir
  use mod_clm_varctl , only : fpftcon
  use mod_clm_varcon , only : tfrz

  implicit none

  private

  save
  !
  ! Vegetation type constants
  !
  character(len=40) , public , dimension(0:mxpft) :: pftname ! PFT description

  !value for not vegetated
  integer(ik4) , public :: noveg
  !value for Needleleaf evergreen temperate tree
  integer(ik4) , public :: ndllf_evr_tmp_tree
  !value for Needleleaf evergreen boreal tree
  integer(ik4) , public :: ndllf_evr_brl_tree
  !value for Needleleaf deciduous boreal tree
  integer(ik4) , public :: ndllf_dcd_brl_tree
  !value for Broadleaf evergreen tropical tree
  integer(ik4) , public :: nbrdlf_evr_trp_tree
  !value for Broadleaf evergreen temperate tree
  integer(ik4) , public :: nbrdlf_evr_tmp_tree
  !value for Broadleaf deciduous tropical tree
  integer(ik4) , public :: nbrdlf_dcd_trp_tree
  !value for Broadleaf deciduous temperate tree
  integer(ik4) , public :: nbrdlf_dcd_tmp_tree
  !value for Broadleaf deciduous boreal tree
  integer(ik4) , public :: nbrdlf_dcd_brl_tree
  !value for last type of tree
  integer(ik4) , public :: ntree
  !value for Broadleaf evergreen shrub
  integer(ik4) , public :: nbrdlf_evr_shrub
  !value for Broadleaf deciduous temperate shrub
  integer(ik4) , public :: nbrdlf_dcd_tmp_shrub
  !value for Broadleaf deciduous boreal shrub
  integer(ik4) , public :: nbrdlf_dcd_brl_shrub
  !value for C3 arctic grass
  integer(ik4) , public :: nc3_arctic_grass
  !value for C3 non-arctic grass
  integer(ik4) , public :: nc3_nonarctic_grass
  !value for C4 grass
  integer(ik4) , public :: nc4_grass
  !value for first crop
  integer(ik4) , public :: npcropmin
  !value for corn, rain fed (rf)
  integer(ik4) , public :: ncorn
  !value for corn, irrigated (ir)
  integer(ik4) , public :: ncornirrig
  !value for spring temperate cereal (rf)
  integer(ik4) , public :: nscereal
  !value for spring temperate cereal (ir)
  integer(ik4) , public :: nscerealirrig
  !value for winter temperate cereal (rf)
  integer(ik4) , public :: nwcereal
  !value for winter temperate cereal (ir)
  integer(ik4) , public :: nwcerealirrig
  !value for soybean (rf)
  integer(ik4) , public :: nsoybean
  !value for soybean (ir)
  integer(ik4) , public :: nsoybeanirrig
  !value for last prognostic crop in list
  integer(ik4) , public :: npcropmax
  !value for generic crop (rf)
  integer(ik4) , public :: nc3crop
  !value for irrigated generic crop (ir)
  integer(ik4) , public :: nc3irrig

  !characteristic leaf dimension (m)
  real(rk8) , public , dimension(0:mxpft) :: dleaf
  !photosynthetic pathway: 0. = c4, 1. = c3
  real(rk8) , public , dimension(0:mxpft) :: c3psn
  !leaf/stem orientation index
  real(rk8) , public , dimension(0:mxpft) :: xl
  !leaf reflectance: 1=vis, 2=nir
  real(rk8) , public , dimension(0:mxpft,numrad) :: rhol
  !stem reflectance: 1=vis, 2=nir
  real(rk8) , public , dimension(0:mxpft,numrad) :: rhos
  !leaf transmittance: 1=vis, 2=nir
  real(rk8) , public , dimension(0:mxpft,numrad) :: taul
  !stem transmittance: 1=vis, 2=nir
  real(rk8) , public , dimension(0:mxpft,numrad) :: taus
  !ratio of momentum roughness length to canopy top height (-)
  real(rk8) , public , dimension(0:mxpft) :: z0mr
  !ratio of displacement height to canopy top height (-)
  real(rk8) , public , dimension(0:mxpft) :: displar
  !CLM rooting distribution parameter [1/m]
  real(rk8) , public , dimension(0:mxpft) :: roota_par
  !CLM rooting distribution parameter [1/m]
  real(rk8) , public , dimension(0:mxpft) :: rootb_par
  ! crop pft: 0. = not crop, 1. = crop pft
  real(rk8) , public , dimension(0:mxpft) :: crop
  ! irrigated pft: 0. = not, 1. = irrigated
  real(rk8) , public , dimension(0:mxpft) :: irrigated
  !soil water potential at full stomatal opening (mm)
  real(rk8) , public , dimension(0:mxpft) :: smpso
  !soil water potential at full stomatal closure (mm)
  real(rk8) , public , dimension(0:mxpft) :: smpsc
  !foliage nitrogen limitation factor (-)
  real(rk8) , public , dimension(0:mxpft) :: fnitr

  ! begin new pft parameters for CN code
  !SLA at top of canopy [m^2/gC]
  real(rk8) , public , dimension(0:mxpft) :: slatop
  !dSLA/dLAI [m^2/gC]
  real(rk8) , public , dimension(0:mxpft) :: dsladlai
  !leaf C:N [gC/gN]
  real(rk8) , public , dimension(0:mxpft) :: leafcn
  !fraction of leaf N in Rubisco [no units]
  real(rk8) , public , dimension(0:mxpft) :: flnr
  !woody lifeform flag (0 or 1)
  real(rk8) , public , dimension(0:mxpft) :: woody
  !leaf litter C:N (gC/gN)
  real(rk8) , public , dimension(0:mxpft) :: lflitcn
  !fine root C:N (gC/gN)
  real(rk8) , public , dimension(0:mxpft) :: frootcn
  !live wood (phloem and ray parenchyma) C:N (gC/gN)
  real(rk8) , public , dimension(0:mxpft) :: livewdcn
  !dead wood (xylem and heartwood) C:N (gC/gN)
  real(rk8) , public , dimension(0:mxpft) :: deadwdcn
  !growth respiration parameter
  real(rk8) , public , dimension(0:mxpft) :: grperc
  !growth respiration parameter
  real(rk8) , public , dimension(0:mxpft) :: grpnow
  !CLM rooting distribution parameter for C and N inputs [unitless]
  real(rk8) , public , dimension(0:mxpft) :: rootprof_beta

  ! for crop
  !grain C:N (gC/gN)
  real(rk8) , public , dimension(0:mxpft) :: graincn
  !parameter used in accFlds
  real(rk8) , public , dimension(0:mxpft) :: mxtmp
  !parameter used in accFlds
  real(rk8) , public , dimension(0:mxpft) :: baset
  !parameter used in CNAllocation
  real(rk8) , public , dimension(0:mxpft) :: declfact
  !parameter used in CNAllocation
  real(rk8) , public , dimension(0:mxpft) :: bfact
  !parameter used in CNAllocation
  real(rk8) , public , dimension(0:mxpft) :: aleaff
  !parameter used in CNAllocation
  real(rk8) , public , dimension(0:mxpft) :: arootf
  !parameter used in CNAllocation
  real(rk8) , public , dimension(0:mxpft) :: astemf
  !parameter used in CNAllocation
  real(rk8) , public , dimension(0:mxpft) :: arooti
  !parameter used in CNAllocation
  real(rk8) , public , dimension(0:mxpft) :: fleafi
  !parameter used in CNAllocation
  real(rk8) , public , dimension(0:mxpft) :: allconsl
  !parameter used in CNAllocation
  real(rk8) , public , dimension(0:mxpft) :: allconss
  !parameter used in CNVegStructUpdate
  real(rk8) , public , dimension(0:mxpft) :: ztopmx
  !parameter used in CNVegStructUpdate
  real(rk8) , public , dimension(0:mxpft) :: laimx
  !parameter used in CNPhenology
  real(rk8) , public , dimension(0:mxpft) :: gddmin
  !parameter used in CNPhenology
  real(rk8) , public , dimension(0:mxpft) :: hybgdd
  !parameter used in CNPhenology
  real(rk8) , public , dimension(0:mxpft) :: lfemerg
  !parameter used in CNPhenology
  real(rk8) , public , dimension(0:mxpft) :: grnfill
  !parameter used in CNPhenology
  integer(ik4) , public , dimension(0:mxpft) :: mxmat
  !minimum planting date for NorthHemisphere (YYYYMMDD)
  integer(ik4) , public , dimension(0:mxpft) :: mnNHplantdate
  !maximum planting date for NorthHemisphere (YYYYMMDD)
  integer(ik4) , public , dimension(0:mxpft) :: mxNHplantdate
  !minimum planting date for SouthHemisphere (YYYYMMDD)
  integer(ik4) , public , dimension(0:mxpft) :: mnSHplantdate
  !maximum planting date for SouthHemisphere (YYYYMMDD)
  integer(ik4) , public , dimension(0:mxpft) :: mxSHplantdate
  !planting temperature used in CNPhenology (K)
  real(rk8) , public , dimension(0:mxpft) :: planttemp
  !mininum planting temperature used in CNPhenology (K)
  real(rk8) , public , dimension(0:mxpft) :: minplanttemp
  !allocation parameter: new fine root C per new leaf C (gC/gC)
  real(rk8) , public , dimension(0:mxpft) :: froot_leaf
  !allocation parameter: new stem c per new leaf C (gC/gC)
  real(rk8) , public , dimension(0:mxpft) :: stem_leaf
  !allocation parameter: new coarse root C per new stem C (gC/gC)
  real(rk8) , public , dimension(0:mxpft) :: croot_stem
  !allocation parameter: fraction of new wood that is live
  ! (phloem and ray parenchyma) (no units)
  real(rk8) , public , dimension(0:mxpft) :: flivewd
  !allocation parameter: fraction of allocation that goes to
  ! currently displayed growth, remainder to storage
  real(rk8) , public , dimension(0:mxpft) :: fcur
  !alternate fcur for use with cndv
  real(rk8) , public , dimension(0:mxpft) :: fcurdv
  !leaf litter labile fraction
  real(rk8) , public , dimension(0:mxpft) :: lf_flab
  !leaf litter cellulose fraction
  real(rk8) , public , dimension(0:mxpft) :: lf_fcel
  !leaf litter lignin fraction
  real(rk8) , public , dimension(0:mxpft) :: lf_flig
  !fine root litter labile fraction
  real(rk8) , public , dimension(0:mxpft) :: fr_flab
  !fine root litter cellulose fraction
  real(rk8) , public , dimension(0:mxpft) :: fr_fcel
  !fine root litter lignin fraction
  real(rk8) , public , dimension(0:mxpft) :: fr_flig
  !leaf longevity (yrs)
  real(rk8) , public , dimension(0:mxpft) :: leaf_long
  !binary flag for evergreen leaf habit (0 or 1)
  real(rk8) , public , dimension(0:mxpft) :: evergreen
  !binary flag for stress-deciduous leaf habit (0 or 1)
  real(rk8) , public , dimension(0:mxpft) :: stress_decid
  !binary flag for seasonal-deciduous leaf habit (0 or 1)
  real(rk8) , public , dimension(0:mxpft) :: season_decid
  !proportion of deadstem to conversion flux
  real(rk8) , public , dimension(0:mxpft) :: pconv
  !proportion of deadstem to 10-yr product pool
  real(rk8) , public , dimension(0:mxpft) :: pprod10
  !proportion of deadstem to 100-yr product pool
  real(rk8) , public , dimension(0:mxpft) :: pprod100
  !harvest mortality proportion of deadstem to 10-yr pool
  real(rk8) , public , dimension(0:mxpft) :: pprodharv10

  ! pft paraemeters for fire code
  real(rk8) , public , dimension(0:mxpft) :: cc_leaf
  real(rk8) , public , dimension(0:mxpft) :: cc_lstem
  real(rk8) , public , dimension(0:mxpft) :: cc_dstem
  real(rk8) , public , dimension(0:mxpft) :: cc_other
  real(rk8) , public , dimension(0:mxpft) :: fm_leaf
  real(rk8) , public , dimension(0:mxpft) :: fm_lstem
  real(rk8) , public , dimension(0:mxpft) :: fm_dstem
  real(rk8) , public , dimension(0:mxpft) :: fm_other
  real(rk8) , public , dimension(0:mxpft) :: fm_root
  real(rk8) , public , dimension(0:mxpft) :: fm_lroot
  real(rk8) , public , dimension(0:mxpft) :: fm_droot
  real(rk8) , public , dimension(0:mxpft) :: fsr_pft
  real(rk8) , public , dimension(0:mxpft) :: fd_pft

  ! pft parameters for crop code
  !fertilizer
  real(rk8) , public , dimension(0:mxpft) :: fertnitro
  !C:N during grain fill; leaf
  real(rk8) , public , dimension(0:mxpft) :: fleafcn
  !C:N during grain fill; fine root
  real(rk8) , public , dimension(0:mxpft) :: ffrootcn
  !C:N during grain fill; stem
  real(rk8) , public , dimension(0:mxpft) :: fstemcn

  ! pft parameters for CNDV code
  ! from LPJ subroutine pftparameters
  !tree maximum crown area (m2)
  real(rk8) , public , dimension(0:mxpft) :: pftpar20
  !min coldest monthly mean temperature
  real(rk8) , public , dimension(0:mxpft) :: pftpar28
  !max coldest monthly mean temperature
  real(rk8) , public , dimension(0:mxpft) :: pftpar29
  !min growing degree days (>= 5 deg C)
  real(rk8) , public , dimension(0:mxpft) :: pftpar30
  !upper limit of temperature of the warmest month (twmax)
  real(rk8) , public , dimension(0:mxpft) :: pftpar31

  !parameter in allometric equation
  real(rk8) , public , parameter :: reinickerp = 1.6D0
  !cn wood density (gC/m3); lpj:2.0e5
  real(rk8) , public , parameter :: dwood  = 2.5D5
  real(rk8) , public , parameter :: allom1 = 100.0D0   !parameters in
  real(rk8) , public , parameter :: allom2 =  40.0D0   !...allometric
  real(rk8) , public , parameter :: allom3 =   0.5D0   !...equations
  real(rk8) , public , parameter :: allom1s = 250.0D0  !modified for shrubs by
  real(rk8) , public , parameter :: allom2s =   8.0D0  !X.D.Z

  public :: pftconrd ! Read and initialize vegetation (PFT) constants

  contains
  !
  ! Read and initialize vegetation (PFT) constants
  !
  subroutine pftconrd
    implicit none
    integer(ik4) :: i           ! loop indices
    type(clm_filetype) :: ncid  ! file handler
    integer(ik4) :: npft        ! number of pfts on pft-physiology file
    character(len=32) :: subname = 'pftconrd'  ! subroutine name
    !
    ! Expected PFT names: The names expected on the fpftcon file and the
    ! order they are expected to be in.
    ! NOTE: similar types are assumed to be together, first trees (ending
    ! with broadleaf_deciduous_boreal_tree then shrubs, ending with
    ! broadleaf_deciduous_boreal_shrub, then grasses starting with
    ! c3_arctic_grass and finally crops, ending with soybean
    !
    ! DO NOT CHANGE THE ORDER -- WITHOUT MODIFYING OTHER PARTS OF
    ! THE CODE WHERE THE ORDER MATTERS!
    !
    character(len=40) , parameter :: expected_pftnames(0:mxpft) = (/ &
                 'not_vegetated                      '  &
               , 'needleleaf_evergreen_temperate_tree'  &
               , 'needleleaf_evergreen_boreal_tree   '  &
               , 'needleleaf_deciduous_boreal_tree   '  &
               , 'broadleaf_evergreen_tropical_tree  '  &
               , 'broadleaf_evergreen_temperate_tree '  &
               , 'broadleaf_deciduous_tropical_tree  '  &
               , 'broadleaf_deciduous_temperate_tree '  &
               , 'broadleaf_deciduous_boreal_tree    '  &
               , 'broadleaf_evergreen_shrub          '  &
               , 'broadleaf_deciduous_temperate_shrub'  &
               , 'broadleaf_deciduous_boreal_shrub   '  &
               , 'c3_arctic_grass                    '  &
               , 'c3_non-arctic_grass                '  &
               , 'c4_grass                           '  &
               , 'c3_crop                            '  &
               , 'c3_irrigated                       '  &
               , 'corn                               '  &
               , 'irrigated_corn                     '  &
               , 'spring_temperate_cereal            '  &
               , 'irrigated_spring_temperate_cereal  '  &
               , 'winter_temperate_cereal            '  &
               , 'irrigated_winter_temperate_cereal  '  &
               , 'soybean                            '  &
               , 'irrigated_soybean                  '  &
    /)

    ! Set specific vegetation type values

    if (myid == italk) then
      write(stdout,*) 'Attempting to read PFT physiological data .....'
    end if
    call clm_openfile(fpftcon,ncid)
    call clm_inqdim(ncid,'pft',npft)
    call clm_readvar(ncid,'pftname',pftname)
    call clm_readvar(ncid,'z0mr',z0mr)
    call clm_readvar(ncid,'displar',displar)
    call clm_readvar(ncid,'dleaf',dleaf)
    call clm_readvar(ncid,'c3psn',c3psn)
    call clm_readvar(ncid,'rholvis',rhol(:,ivis))
    call clm_readvar(ncid,'rholnir',rhol(:,inir))
    call clm_readvar(ncid,'rhosvis',rhos(:,ivis))
    call clm_readvar(ncid,'rhosnir',rhos(:,inir))
    call clm_readvar(ncid,'taulvis',taul(:,ivis))
    call clm_readvar(ncid,'taulnir',taul(:,inir))
    call clm_readvar(ncid,'tausvis',taus(:,ivis))
    call clm_readvar(ncid,'tausnir',taus(:,inir))
    call clm_readvar(ncid,'xl',xl)
    call clm_readvar(ncid,'roota_par',roota_par)
    call clm_readvar(ncid,'rootb_par',rootb_par)
    call clm_readvar(ncid,'slatop',slatop)
    call clm_readvar(ncid,'dsladlai',dsladlai)
    call clm_readvar(ncid,'leafcn',leafcn)
    call clm_readvar(ncid,'flnr',flnr)
    call clm_readvar(ncid,'smpso',smpso)
    call clm_readvar(ncid,'smpsc',smpsc)
    call clm_readvar(ncid,'fnitr',fnitr)
    call clm_readvar(ncid,'woody',woody)
    call clm_readvar(ncid,'lflitcn',lflitcn)
    call clm_readvar(ncid,'frootcn',frootcn)
    call clm_readvar(ncid,'livewdcn',livewdcn)
    call clm_readvar(ncid,'deadwdcn',deadwdcn)
    call clm_readvar(ncid,'grperc',grperc)
    call clm_readvar(ncid,'grpnow',grpnow)
    call clm_readvar(ncid,'froot_leaf',froot_leaf)
    call clm_readvar(ncid,'stem_leaf',stem_leaf)
    call clm_readvar(ncid,'croot_stem',croot_stem)
    call clm_readvar(ncid,'flivewd',flivewd)
    call clm_readvar(ncid,'fcur',fcur)
    call clm_readvar(ncid,'fcurdv',fcurdv)
    call clm_readvar(ncid,'lf_flab',lf_flab)
    call clm_readvar(ncid,'lf_fcel',lf_fcel)
    call clm_readvar(ncid,'lf_flig',lf_flig)
    call clm_readvar(ncid,'fr_flab',fr_flab)
    call clm_readvar(ncid,'fr_fcel',fr_fcel)
    call clm_readvar(ncid,'fr_flig',fr_flig)
    call clm_readvar(ncid,'leaf_long',leaf_long)
    call clm_readvar(ncid,'evergreen',evergreen)
    call clm_readvar(ncid,'stress_decid',stress_decid)
    call clm_readvar(ncid,'season_decid',season_decid)
    call clm_readvar(ncid,'pftpar20',pftpar20)
    call clm_readvar(ncid,'pftpar28',pftpar28)
    call clm_readvar(ncid,'pftpar29',pftpar29)
    call clm_readvar(ncid,'pftpar30',pftpar30)
    call clm_readvar(ncid,'pftpar31',pftpar31)
    call clm_readvar(ncid,'fertnitro',fertnitro)
    call clm_readvar(ncid,'fleafcn',fleafcn)
    call clm_readvar(ncid,'ffrootcn',ffrootcn)
    call clm_readvar(ncid,'fstemcn',fstemcn)
#ifdef VERTSOILC
    call clm_readvar(ncid,'rootprof_beta',rootprof_beta)
#endif
    call clm_readvar(ncid,'pconv',pconv)
    call clm_readvar(ncid,'pprod10',pprod10)
    call clm_readvar(ncid,'pprodharv10',pprodharv10)
    call clm_readvar(ncid,'pprod100',pprod100)
    call clm_readvar(ncid,'graincn',graincn)
    call clm_readvar(ncid,'mxtmp',mxtmp)
    call clm_readvar(ncid,'baset',baset)
    call clm_readvar(ncid,'declfact',declfact)
    call clm_readvar(ncid,'bfact',bfact)
    call clm_readvar(ncid,'aleaff',aleaff)
    call clm_readvar(ncid,'arootf',arootf)
    call clm_readvar(ncid,'astemf',astemf)
    call clm_readvar(ncid,'arooti',arooti)
    call clm_readvar(ncid,'fleafi',fleafi)
    call clm_readvar(ncid,'allconsl',allconsl)
    call clm_readvar(ncid,'allconss',allconss)
    call clm_readvar(ncid,'crop',crop)
    call clm_readvar(ncid,'irrigated',irrigated)
    call clm_readvar(ncid,'ztopmx',ztopmx)
    call clm_readvar(ncid,'laimx',laimx)
    call clm_readvar(ncid,'gddmin',gddmin)
    call clm_readvar(ncid,'hybgdd',hybgdd)
    call clm_readvar(ncid,'lfemerg',lfemerg)
    call clm_readvar(ncid,'grnfill',grnfill)
    call clm_readvar(ncid,'mxmat',mxmat)
    call clm_readvar(ncid,'cc_leaf',cc_leaf)
    call clm_readvar(ncid,'cc_lstem',cc_lstem)
    call clm_readvar(ncid,'cc_dstem',cc_dstem)
    call clm_readvar(ncid,'cc_other',cc_other)
    call clm_readvar(ncid,'fm_leaf',fm_leaf)
    call clm_readvar(ncid,'fm_lstem',fm_lstem)
    call clm_readvar(ncid,'fm_dstem',fm_dstem)
    call clm_readvar(ncid,'fm_other',fm_other)
    call clm_readvar(ncid,'fm_root',fm_root)
    call clm_readvar(ncid,'fm_lroot',fm_lroot)
    call clm_readvar(ncid,'fm_droot',fm_droot)
    call clm_readvar(ncid,'fsr_pft',fsr_pft)
    call clm_readvar(ncid,'fd_pft',fd_pft)
    call clm_readvar(ncid,'planting_temp',planttemp)
    call clm_readvar(ncid,'min_planting_temp',minplanttemp)
    call clm_readvar(ncid,'min_NH_planting_date',mnNHplantdate)
    call clm_readvar(ncid,'min_SH_planting_date',mnSHplantdate)
    call clm_readvar(ncid,'max_NH_planting_date',mxNHplantdate)
    call clm_readvar(ncid,'max_SH_planting_date',mxSHplantdate)
    call clm_closefile(ncid)

    do i = 0 , mxpft
      if ( trim(adjustl(pftname(i))) /= trim(expected_pftnames(i)) )then
        write(stderr,*)'pftconrd: pftname is NOT what is expected, name = ', &
             trim(pftname(i)), ', expected name = ', trim(expected_pftnames(i))
        call fatal(__FILE__,__LINE__, &
                'pftconrd: bad name for pft on fpftcon dataset' )
      end if
      if ( trim(pftname(i)) == 'not_vegetated'                       ) &
              noveg               = i
      if ( trim(pftname(i)) == 'needleleaf_evergreen_temperate_tree' ) &
              ndllf_evr_tmp_tree  = i
      if ( trim(pftname(i)) == 'needleleaf_evergreen_boreal_tree'    ) &
              ndllf_evr_brl_tree  = i
      if ( trim(pftname(i)) == 'needleleaf_deciduous_boreal_tree'    ) &
              ndllf_dcd_brl_tree  = i
      if ( trim(pftname(i)) == 'broadleaf_evergreen_tropical_tree'   ) &
              nbrdlf_evr_trp_tree  = i
      if ( trim(pftname(i)) == 'broadleaf_evergreen_temperate_tree'  ) &
              nbrdlf_evr_tmp_tree  = i
      if ( trim(pftname(i)) == 'broadleaf_deciduous_tropical_tree'   ) &
              nbrdlf_dcd_trp_tree  = i
      if ( trim(pftname(i)) == 'broadleaf_deciduous_temperate_tree'  ) &
              nbrdlf_dcd_tmp_tree  = i
      if ( trim(pftname(i)) == 'broadleaf_deciduous_boreal_tree'     ) &
              nbrdlf_dcd_brl_tree  = i
      if ( trim(pftname(i)) == 'broadleaf_evergreen_shrub'           ) &
              nbrdlf_evr_shrub     = i
      if ( trim(pftname(i)) == 'broadleaf_deciduous_temperate_shrub' ) &
              nbrdlf_dcd_tmp_shrub = i
      if ( trim(pftname(i)) == 'broadleaf_deciduous_boreal_shrub'    ) &
              nbrdlf_dcd_brl_shrub = i
      if ( trim(pftname(i)) == 'c3_arctic_grass'                     ) &
              nc3_arctic_grass     = i
      if ( trim(pftname(i)) == 'c3_non-arctic_grass'                 ) &
              nc3_nonarctic_grass  = i
      if ( trim(pftname(i)) == 'c4_grass'                            ) &
              nc4_grass            = i
      if ( trim(pftname(i)) == 'c3_crop'                             ) &
              nc3crop              = i
      if ( trim(pftname(i)) == 'c3_irrigated'                        ) &
              nc3irrig               = i
      if ( trim(pftname(i)) == 'corn'                                ) &
              ncorn                = i
      if ( trim(pftname(i)) == 'irrigated_corn'                      ) &
              ncornirrig           = i
      if ( trim(pftname(i)) == 'spring_temperate_cereal'             ) &
              nscereal             = i
      if ( trim(pftname(i)) == 'irrigated_spring_temperate_cereal'   ) &
              nscerealirrig        = i
      if ( trim(pftname(i)) == 'winter_temperate_cereal'             ) &
              nwcereal             = i
      if ( trim(pftname(i)) == 'irrigated_winter_temperate_cereal'   ) &
              nwcerealirrig        = i
      if ( trim(pftname(i)) == 'soybean'                             ) &
              nsoybean             = i
      if ( trim(pftname(i)) == 'irrigated_soybean'                   ) &
              nsoybeanirrig        = i
    end do

    ntree       = nbrdlf_dcd_brl_tree  ! value for last type of tree
    npcropmin   = ncorn                ! first prognostic crop
    npcropmax   = nsoybeanirrig        ! last prognostic crop in list

#if (defined CNDV)
    fcur(:) = fcurdv(:)
#endif
    !
    ! Do some error checking
    !
    if ( npcropmax /= mxpft )then
      call fatal(__FILE__,__LINE__, &
              trim(subname)//' ERROR: npcropmax is NOT the last value' )
    end if
    do i = 0 , mxpft
      if (      irrigated(i) == 1.0D0  .and. (i == nc3irrig .or. &
                                              i == ncornirrig .or. &
                                              i == nscerealirrig .or. &
                                              i == nwcerealirrig .or. &
                                              i == nsoybeanirrig) ) then
        ! correct
      else if ( irrigated(i) == 0.0D0 ) then
        ! correct
      else
        call fatal(__FILE__,__LINE__, &
                trim(subname)//' ERROR: irrigated has wrong values' )
      end if
      if (      crop(i) == 1.0D0 .and. &
              (i >= nc3crop .and. i <= npcropmax) ) then
        ! correct
      else if ( crop(i) == 0.0D0 ) then
        ! correct
      else
        call fatal(__FILE__,__LINE__, &
                trim(subname)//' ERROR: crop has wrong values' )
      end if
      if ( (i /= noveg) .and. (i < npcropmin) .and. &
           abs(pconv(i)+pprod10(i)+pprod100(i) - 1.0D0) > 1.D-7 ) then
        call fatal(__FILE__,__LINE__, &
            trim(subname)//' ERROR: pconv+pprod10+pprod100 do NOT sum to one.' )
      end if
      if ( pprodharv10(i) > 1.0D0 .or. pprodharv10(i) < 0.0D0 ) then
        call fatal(__FILE__,__LINE__, &
                trim(subname)//' ERROR: pprodharv10 outside of range.' )
      end if
    end do

    if (myid == italk) then
      write(stdout,*) 'Successfully read PFT physiological data'
    end if

  end subroutine pftconrd

end module mod_clm_pftvarcon

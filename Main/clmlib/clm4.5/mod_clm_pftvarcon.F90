module mod_clm_pftvarcon

!-----------------------------------------------------------------------
!BOP
!
! !MODULE: pftvarcon
!
! !DESCRIPTION:
! Module containing vegetation constants and method to
! read and initialize vegetation (PFT) constants.
!
! !USES:
  use mod_realkinds
  use mod_mpmessage
  use mod_clm_nchelper
  use mod_dynparam
  use mod_mppparam
  use mod_clm_varpar  , only : mxpft, numrad, ivis, inir
!
! !PUBLIC TYPES:
  implicit none
  save
!
! Vegetation type constants
!
  character(len=40) pftname(0:mxpft) !PFT description

  integer :: noveg                  !value for not vegetated 
  integer :: ndllf_evr_tmp_tree     !value for Needleleaf evergreen temperate tree
  integer :: ndllf_evr_brl_tree     !value for Needleleaf evergreen boreal tree
  integer :: ndllf_dcd_brl_tree     !value for Needleleaf deciduous boreal tree
  integer :: nbrdlf_evr_trp_tree    !value for Broadleaf evergreen tropical tree
  integer :: nbrdlf_evr_tmp_tree    !value for Broadleaf evergreen temperate tree
  integer :: nbrdlf_dcd_trp_tree    !value for Broadleaf deciduous tropical tree
  integer :: nbrdlf_dcd_tmp_tree    !value for Broadleaf deciduous temperate tree
  integer :: nbrdlf_dcd_brl_tree    !value for Broadleaf deciduous boreal tree
  integer :: ntree                  !value for last type of tree
  integer :: nbrdlf_evr_shrub       !value for Broadleaf evergreen shrub
  integer :: nbrdlf_dcd_tmp_shrub   !value for Broadleaf deciduous temperate shrub
  integer :: nbrdlf_dcd_brl_shrub   !value for Broadleaf deciduous boreal shrub
  integer :: nc3_arctic_grass       !value for C3 arctic grass
  integer :: nc3_nonarctic_grass    !value for C3 non-arctic grass
  integer :: nc4_grass              !value for C4 grass
  integer :: npcropmin              !value for first crop
  integer :: ncorn                  !value for corn, rain fed (rf)
  integer :: ncornirrig             !value for corn, irrigated (ir)
  integer :: nscereal               !value for spring temperate cereal (rf)
  integer :: nscerealirrig          !value for spring temperate cereal (ir)
  integer :: nwcereal               !value for winter temperate cereal (rf)
  integer :: nwcerealirrig          !value for winter temperate cereal (ir)
  integer :: nsoybean               !value for soybean (rf)
  integer :: nsoybeanirrig          !value for soybean (ir)
  integer :: npcropmax              !value for last prognostic crop in list
  integer :: nc3crop                !value for generic crop (rf)
  integer :: nc3irrig               !value for irrigated generic crop (ir)

  real(rk8):: dleaf(0:mxpft)       !characteristic leaf dimension (m)
  real(rk8):: c3psn(0:mxpft)       !photosynthetic pathway: 0. = c4, 1. = c3
  real(rk8):: xl(0:mxpft)          !leaf/stem orientation index
  real(rk8):: rhol(0:mxpft,numrad) !leaf reflectance: 1=vis, 2=nir
  real(rk8):: rhos(0:mxpft,numrad) !stem reflectance: 1=vis, 2=nir
  real(rk8):: taul(0:mxpft,numrad) !leaf transmittance: 1=vis, 2=nir
  real(rk8):: taus(0:mxpft,numrad) !stem transmittance: 1=vis, 2=nir
  real(rk8):: z0mr(0:mxpft)        !ratio of momentum roughness length to canopy top height (-)
  real(rk8):: displar(0:mxpft)     !ratio of displacement height to canopy top height (-)
  real(rk8):: roota_par(0:mxpft)   !CLM rooting distribution parameter [1/m]
  real(rk8):: rootb_par(0:mxpft)   !CLM rooting distribution parameter [1/m]
  real(rk8):: crop(0:mxpft)        ! crop pft: 0. = not crop, 1. = crop pft
  real(rk8):: irrigated(0:mxpft)   ! irrigated pft: 0. = not, 1. = irrigated
  real(rk8):: smpso(0:mxpft)       !soil water potential at full stomatal opening (mm)
  real(rk8):: smpsc(0:mxpft)       !soil water potential at full stomatal closure (mm)
  real(rk8):: fnitr(0:mxpft)       !foliage nitrogen limitation factor (-)
  ! begin new pft parameters for CN code
  real(rk8):: slatop(0:mxpft)      !SLA at top of canopy [m^2/gC]
  real(rk8):: dsladlai(0:mxpft)    !dSLA/dLAI [m^2/gC]
  real(rk8):: leafcn(0:mxpft)      !leaf C:N [gC/gN]
  real(rk8):: flnr(0:mxpft)        !fraction of leaf N in Rubisco [no units]
  real(rk8):: woody(0:mxpft)       !woody lifeform flag (0 or 1)
  real(rk8):: lflitcn(0:mxpft)      !leaf litter C:N (gC/gN)
  real(rk8):: frootcn(0:mxpft)      !fine root C:N (gC/gN)
  real(rk8):: livewdcn(0:mxpft)     !live wood (phloem and ray parenchyma) C:N (gC/gN)
  real(rk8):: deadwdcn(0:mxpft)     !dead wood (xylem and heartwood) C:N (gC/gN)
  real(rk8):: grperc(0:mxpft)       !growth respiration parameter
  real(rk8):: grpnow(0:mxpft)       !growth respiration parameter
  real(rk8):: rootprof_beta(0:mxpft)   !CLM rooting distribution parameter for C and N inputs [unitless]

! for crop
  real(rk8):: graincn(0:mxpft)      !grain C:N (gC/gN)
  real(rk8):: mxtmp(0:mxpft)        !parameter used in accFlds
  real(rk8):: baset(0:mxpft)        !parameter used in accFlds
  real(rk8):: declfact(0:mxpft)     !parameter used in CNAllocation
  real(rk8):: bfact(0:mxpft)        !parameter used in CNAllocation
  real(rk8):: aleaff(0:mxpft)       !parameter used in CNAllocation
  real(rk8):: arootf(0:mxpft)       !parameter used in CNAllocation
  real(rk8):: astemf(0:mxpft)       !parameter used in CNAllocation
  real(rk8):: arooti(0:mxpft)       !parameter used in CNAllocation
  real(rk8):: fleafi(0:mxpft)       !parameter used in CNAllocation
  real(rk8):: allconsl(0:mxpft)     !parameter used in CNAllocation
  real(rk8):: allconss(0:mxpft)     !parameter used in CNAllocation
  real(rk8):: ztopmx(0:mxpft)       !parameter used in CNVegStructUpdate
  real(rk8):: laimx(0:mxpft)        !parameter used in CNVegStructUpdate
  real(rk8):: gddmin(0:mxpft)       !parameter used in CNPhenology
  real(rk8):: hybgdd(0:mxpft)       !parameter used in CNPhenology
  real(rk8):: lfemerg(0:mxpft)      !parameter used in CNPhenology
  real(rk8):: grnfill(0:mxpft)      !parameter used in CNPhenology
  integer :: mxmat(0:mxpft)        !parameter used in CNPhenology
  integer :: mnNHplantdate(0:mxpft)!minimum planting date for NorthHemisphere (YYYYMMDD)
  integer :: mxNHplantdate(0:mxpft)!maximum planting date for NorthHemisphere (YYYYMMDD)
  integer :: mnSHplantdate(0:mxpft)!minimum planting date for SouthHemisphere (YYYYMMDD)
  integer :: mxSHplantdate(0:mxpft)!maximum planting date for SouthHemisphere (YYYYMMDD)
  real(rk8):: planttemp(0:mxpft)    !planting temperature used in CNPhenology (K)
  real(rk8):: minplanttemp(0:mxpft) !mininum planting temperature used in CNPhenology (K)
  real(rk8):: froot_leaf(0:mxpft)   !allocation parameter: new fine root C per new leaf C (gC/gC) 
  real(rk8):: stem_leaf(0:mxpft)    !allocation parameter: new stem c per new leaf C (gC/gC)
  real(rk8):: croot_stem(0:mxpft)   !allocation parameter: new coarse root C per new stem C (gC/gC)
  real(rk8):: flivewd(0:mxpft)      !allocation parameter: fraction of new wood that is live (phloem and ray parenchyma) (no units)
  real(rk8):: fcur(0:mxpft)         !allocation parameter: fraction of allocation that goes to currently displayed growth, remainder to storage
  real(rk8):: fcurdv(0:mxpft)       !alternate fcur for use with cndv
  real(rk8):: lf_flab(0:mxpft)      !leaf litter labile fraction
  real(rk8):: lf_fcel(0:mxpft)      !leaf litter cellulose fraction
  real(rk8):: lf_flig(0:mxpft)      !leaf litter lignin fraction
  real(rk8):: fr_flab(0:mxpft)      !fine root litter labile fraction
  real(rk8):: fr_fcel(0:mxpft)      !fine root litter cellulose fraction
  real(rk8):: fr_flig(0:mxpft)      !fine root litter lignin fraction
  real(rk8):: leaf_long(0:mxpft)    !leaf longevity (yrs)
  real(rk8):: evergreen(0:mxpft)    !binary flag for evergreen leaf habit (0 or 1)
  real(rk8):: stress_decid(0:mxpft) !binary flag for stress-deciduous leaf habit (0 or 1)
  real(rk8):: season_decid(0:mxpft) !binary flag for seasonal-deciduous leaf habit (0 or 1)
  real(rk8):: pconv(0:mxpft)        !proportion of deadstem to conversion flux
  real(rk8):: pprod10(0:mxpft)      !proportion of deadstem to 10-yr product pool
  real(rk8):: pprod100(0:mxpft)     !proportion of deadstem to 100-yr product pool
  real(rk8):: pprodharv10(0:mxpft)  !harvest mortality proportion of deadstem to 10-yr pool
  ! pft paraemeters for fire code
  real(rk8):: cc_leaf(0:mxpft)
  real(rk8):: cc_lstem(0:mxpft)
  real(rk8):: cc_dstem(0:mxpft)
  real(rk8):: cc_other(0:mxpft)
  real(rk8):: fm_leaf(0:mxpft)
  real(rk8):: fm_lstem(0:mxpft)
  real(rk8):: fm_dstem(0:mxpft)
  real(rk8):: fm_other(0:mxpft)
  real(rk8):: fm_root(0:mxpft)
  real(rk8):: fm_lroot(0:mxpft)
  real(rk8):: fm_droot(0:mxpft)
  real(rk8):: fsr_pft(0:mxpft)
  real(rk8):: fd_pft(0:mxpft)
  ! pft parameters for crop code
  real(rk8):: fertnitro(0:mxpft)    !fertilizer
  real(rk8):: fleafcn(0:mxpft)      !C:N during grain fill; leaf
  real(rk8):: ffrootcn(0:mxpft)     !C:N during grain fill; fine root
  real(rk8):: fstemcn(0:mxpft)      !C:N during grain fill; stem

  ! pft parameters for CNDV code
  ! from LPJ subroutine pftparameters
  real(rk8) pftpar20(0:mxpft)       !tree maximum crown area (m2)
  real(rk8) pftpar28(0:mxpft)       !min coldest monthly mean temperature
  real(rk8) pftpar29(0:mxpft)       !max coldest monthly mean temperature
  real(rk8) pftpar30(0:mxpft)       !min growing degree days (>= 5 deg C)
  real(rk8) pftpar31(0:mxpft)       !upper limit of temperature of the warmest month (twmax)
  real(rk8), parameter :: reinickerp = 1.6D0 !parameter in allometric equation
  real(rk8), parameter :: dwood  = 2.5D5   !cn wood density (gC/m3); lpj:2.0e5
  real(rk8), parameter :: allom1 = 100.0D0   !parameters in
  real(rk8), parameter :: allom2 =  40.0D0   !...allometric
  real(rk8), parameter :: allom3 =   0.5D0   !...equations
  real(rk8), parameter :: allom1s = 250.0D0  !modified for shrubs by
  real(rk8), parameter :: allom2s =   8.0D0  !X.D.Z
!
! !PUBLIC MEMBER FUNCTIONS:
  public :: pftconrd ! Read and initialize vegetation (PFT) constants
!
! !REVISION HISTORY:
! Created by Sam Levis (put into module form by Mariana Vertenstein)
! 10/21/03, Peter Thornton: Added new variables for CN code
! 06/24/09, Erik Kluzek: Add indices for all pft types, and add expected_pftnames array and comparision
! 09/17/10, David Lawrence: Modified code to read in netCDF pft physiology file
!
!EOP
!-----------------------------------------------------------------------

contains

!-----------------------------------------------------------------------
!BOP
!
! !IROUTINE: pftconrd
!
! !INTERFACE:
  subroutine pftconrd
!
! !DESCRIPTION:
! Read and initialize vegetation (PFT) constants
!
! !USES:
    use mod_clm_varctl, only : fpftcon
    use mod_clm_varcon, only : tfrz
!
! !ARGUMENTS:
    implicit none
!
! !CALLED FROM:
! routine initialize in module initializeMod
!
! !REVISION HISTORY:
! Created by Gordon Bonan
!F. Li and S. Levis (11/06/12)
!
! !LOCAL VARIABLES:
!EOP
    integer :: i,n              ! loop indices
    integer :: ier              ! error code
    type(clm_filetype) :: ncid
    integer :: npft             ! number of pfts on pft-physiology file
    logical :: readv            ! read variable in or not
    character(len=32) :: subname = 'pftconrd'              ! subroutine name
    !
    ! Expected PFT names: The names expected on the fpftcon file and the order they are expected to be in.
    ! NOTE: similar types are assumed to be together, first trees (ending with broadleaf_deciduous_boreal_tree
    !       then shrubs, ending with broadleaf_deciduous_boreal_shrub, then grasses starting with c3_arctic_grass
    !       and finally crops, ending with soybean
    ! DO NOT CHANGE THE ORDER -- WITHOUT MODIFYING OTHER PARTS OF THE CODE WHERE THE ORDER MATTERS!
    !
    character(len=40), parameter :: expected_pftnames(0:mxpft) = (/ &
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
!-----------------------------------------------------------------------

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

    do i = 0, mxpft
       if ( trim(adjustl(pftname(i))) /= trim(expected_pftnames(i)) )then
          write(stderr,*)'pftconrd: pftname is NOT what is expected, name = ', &
                        trim(pftname(i)), ', expected name = ', trim(expected_pftnames(i))
          call fatal(__FILE__,__LINE__,'pftconrd: bad name for pft on fpftcon dataset' )
       end if
       if ( trim(pftname(i)) == 'not_vegetated'                       ) noveg               = i
       if ( trim(pftname(i)) == 'needleleaf_evergreen_temperate_tree' ) ndllf_evr_tmp_tree  = i
       if ( trim(pftname(i)) == 'needleleaf_evergreen_boreal_tree'    ) ndllf_evr_brl_tree  = i
       if ( trim(pftname(i)) == 'needleleaf_deciduous_boreal_tree'    ) ndllf_dcd_brl_tree  = i
       if ( trim(pftname(i)) == 'broadleaf_evergreen_tropical_tree'   ) nbrdlf_evr_trp_tree  = i
       if ( trim(pftname(i)) == 'broadleaf_evergreen_temperate_tree'  ) nbrdlf_evr_tmp_tree  = i
       if ( trim(pftname(i)) == 'broadleaf_deciduous_tropical_tree'   ) nbrdlf_dcd_trp_tree  = i
       if ( trim(pftname(i)) == 'broadleaf_deciduous_temperate_tree'  ) nbrdlf_dcd_tmp_tree  = i
       if ( trim(pftname(i)) == 'broadleaf_deciduous_boreal_tree'     ) nbrdlf_dcd_brl_tree  = i
       if ( trim(pftname(i)) == 'broadleaf_evergreen_shrub'           ) nbrdlf_evr_shrub     = i
       if ( trim(pftname(i)) == 'broadleaf_deciduous_temperate_shrub' ) nbrdlf_dcd_tmp_shrub = i
       if ( trim(pftname(i)) == 'broadleaf_deciduous_boreal_shrub'    ) nbrdlf_dcd_brl_shrub = i
       if ( trim(pftname(i)) == 'c3_arctic_grass'                     ) nc3_arctic_grass     = i
       if ( trim(pftname(i)) == 'c3_non-arctic_grass'                 ) nc3_nonarctic_grass  = i
       if ( trim(pftname(i)) == 'c4_grass'                            ) nc4_grass            = i
       if ( trim(pftname(i)) == 'c3_crop'                             ) nc3crop              = i
       if ( trim(pftname(i)) == 'c3_irrigated'                        ) nc3irrig               = i
       if ( trim(pftname(i)) == 'corn'                                ) ncorn                = i
       if ( trim(pftname(i)) == 'irrigated_corn'                      ) ncornirrig           = i
       if ( trim(pftname(i)) == 'spring_temperate_cereal'             ) nscereal             = i
       if ( trim(pftname(i)) == 'irrigated_spring_temperate_cereal'   ) nscerealirrig        = i
       if ( trim(pftname(i)) == 'winter_temperate_cereal'             ) nwcereal             = i
       if ( trim(pftname(i)) == 'irrigated_winter_temperate_cereal'   ) nwcerealirrig        = i
       if ( trim(pftname(i)) == 'soybean'                             ) nsoybean             = i
       if ( trim(pftname(i)) == 'irrigated_soybean'                   ) nsoybeanirrig        = i
    end do

    ntree                = nbrdlf_dcd_brl_tree  ! value for last type of tree
    npcropmin            = ncorn                ! first prognostic crop
    npcropmax            = nsoybeanirrig        ! last prognostic crop in list

#if (defined CNDV)
    fcur(:) = fcurdv(:)
#endif

    !
    ! Do some error checking
    !
    if ( npcropmax /= mxpft )then
       call fatal(__FILE__,__LINE__,trim(subname)//' ERROR: npcropmax is NOT the last value' )
    end if
    do i = 0, mxpft
       if (      irrigated(i) == 1.0D0  .and. (i == nc3irrig .or. &
                                                i == ncornirrig .or. &
                                                i == nscerealirrig .or. &
                                                i == nwcerealirrig .or. &
                                                i == nsoybeanirrig) )then
          ! correct
       else if ( irrigated(i) == 0.0D0 )then
          ! correct
       else
          call fatal(__FILE__,__LINE__,trim(subname)//' ERROR: irrigated has wrong values' )
       end if
       if (      crop(i) == 1.0D0 .and. (i >= nc3crop .and. i <= npcropmax) )then
          ! correct
       else if ( crop(i) == 0.0D0 )then
          ! correct
       else
          call fatal(__FILE__,__LINE__,trim(subname)//' ERROR: crop has wrong values' )
       end if
       if ( (i /= noveg) .and. (i < npcropmin) .and. &
            abs(pconv(i)+pprod10(i)+pprod100(i) - 1.0D0) > 1.D-7 )then
          call fatal(__FILE__,__LINE__,trim(subname)//' ERROR: pconv+pprod10+pprod100 do NOT sum to one.' )
       end if
       if ( pprodharv10(i) > 1.0D0 .or. pprodharv10(i) < 0.0D0 )then
          call fatal(__FILE__,__LINE__,trim(subname)//' ERROR: pprodharv10 outside of range.' )
       end if
    end do

    if (myid == italk) then
       write(stdout,*) 'Successfully read PFT physiological data'
       write(stdout,*)
    end if

  end subroutine pftconrd

end module mod_clm_pftvarcon
